from dolfin import *
from utilities import *
from energies import *
from data import *
import sys
import numpy as np
import math
import os
import shutil
import pickle
import Borders
import MyProblems
import Solvers

if not has_petsc_tao():
    print("DOLFIN must be compiled with TAO to run this demo.")
    exit(0)

if len(sys.argv) != 3:
	print "\n Incorrect number of input parameters. Usage: \n python inverse_back.py l l0 \n l = TARGET LENGTH \n l0  = INITIAL GUESS. \n"
	exit(0)

# Parsing the command line input parameters
l = float(sys.argv[1])
ltilde = float(sys.argv[2])

os.system("python direct.py %f" %l)

# Boundary instantiation
down = Borders.Down(L, H)
up = Borders.Up(L, H)
right = Borders.Right(L, H)
left = Borders.Left(L, H)

boundaries = FacetFunction("size_t",Th)
boundaries.set_all(0)

down.mark(boundaries, 1)
right.mark(boundaries, 2)
up.mark(boundaries, 3)
left.mark(boundaries, 4)

#Finite Element Spaces
V_u = VectorFunctionSpace(Th, "CG", 1, ndim)
V_alpha = FunctionSpace(Th, "CG", 1)

# Border conditions for alpha
Gamma_alpha_1 = DirichletBC(V_alpha, 0.0, boundaries, 1)
Gamma_alpha_2 = DirichletBC(V_alpha, 0.0, boundaries, 2)
Gamma_alpha_3 = DirichletBC(V_alpha, 0.0, boundaries, 3)
Gamma_alpha_4 = DirichletBC(V_alpha, 0.0, boundaries, 4)

bc_alpha = [Gamma_alpha_1,Gamma_alpha_3,Gamma_alpha_4]

# Border conditions for u
Gamma_u_1 = DirichletBC(V_u.sub(0),  0.0, boundaries, 1)
Gamma_u_3 = DirichletBC(V_u, [0.0, 0.0], boundaries, 3)
Gamma_u_2 = DirichletBC(V_u.sub(0),  0.0, boundaries, 2)
Gamma_u_4 = DirichletBC(V_u.sub(0),  0.0, boundaries, 4)

bc_u = [Gamma_u_1, Gamma_u_3, Gamma_u_4]

# Border conditions for T
Gamma_T_1 = DirichletBC(V_alpha, T0-deltaT, boundaries, 1)
Gamma_T_2 = DirichletBC(V_alpha, T0-deltaT, boundaries, 2)
Gamma_T_3 = DirichletBC(V_alpha, T0-deltaT, boundaries, 3)
Gamma_T_4 = DirichletBC(V_alpha, T0-deltaT, boundaries, 4)

####################################################################################
# TEMPERATURE PROBLEM ##############################################################

Tsol, T, w = Function(V_alpha), TrialFunction(V_alpha), TestFunction(V_alpha)
T_old = interpolate(Expression("0.0"), V_alpha)

####################################################################################
# DIRECT PROBLEM ###################################################################

u, du, v = Function(V_u), TrialFunction(V_u), TestFunction(V_u)
alpha, dalpha, theta = Function(V_alpha), TrialFunction(V_alpha), TestFunction(V_alpha)

# Error variable
alpha_error = Function(V_alpha)
# Recover variable
alphaM = Function(V_alpha)

####################################################################################
# ADJOINT PROBLEM ##################################################################

V, dV, W = Function(V_u), TrialFunction(V_u), TestFunction(V_u)
Beta, dBeta, Z = Function(V_alpha), TrialFunction(V_alpha), TestFunction(V_alpha)

# Error variable
Beta_error = Function(V_alpha)

#####################################################################################

lerror = 1.
P = 1./250
iterinv = 0

One = interpolate(Expression("1.0"), V_alpha)
EPS = interpolate(Expression("0.001"), V_alpha)

# UPDATE ITERATIONS #################################################################

while lerror > 1e-3 and iterinv < 100:

	# Inizialization
	
	alpha_0 = interpolate(Expression("0.0"), V_alpha)
	Beta_0 = interpolate(Expression("1.0"), V_alpha)
	T_old = interpolate(Expression("0.0"), V_alpha)

	lb_Beta = Function(interpolate(Constant(0.),V_alpha))
	ub_Beta = Function(interpolate(Constant(1.),V_alpha))

	lb_alpha = Function(interpolate(Constant(0.),V_alpha))
	ub_alpha = Function(interpolate(Constant(1.),V_alpha))

	t = dt;
	timestep = 0

	# Initialize the temperature solver
	TS = Solvers.TemperatureSolver(T, w, T_old, Gamma_T_2)

	# Initialize the direct solver
	DS = Solvers.DirectSolver(ltilde, u, alpha, v, theta, du, dalpha, Tsol, bc_u, bc_alpha)

	while t<=tf:

		timestep +=1

		TS.solve(Tsol)
		
		DS.solve(u, alpha, alpha_0, lb_alpha, ub_alpha, alpha_error)

		lb_alpha.vector()[:] = alpha.vector()
		T_old.assign(Tsol)

		# Dumping solutions
		f = open('results/alpha%d.pckl'%(timestep), 'w')
    		pickle.dump(alpha.vector().array(), f)
		f.close()

		f = open('results/u%d.pckl'%(timestep), 'w')
    		pickle.dump(u.vector().array(), f)
		f.close()

		f = open('results/T%d.pckl'%(timestep), 'w')
    		pickle.dump(Tsol.vector().array(), f)
		f.close()

		t+=dt
		
		#################################################################
		#################################################################

	timestep = 17	
	t = tf
	der = 0
	IntA = 0.
	kk = 3

	while t> tf - kk*dt:

		timestep -=1

		f = open('results/alphaM%d.pckl'%(timestep))
		saved_data_alphaM = pickle.load(f)
		f.close()

		f = open('results/alpha%d.pckl'%(timestep))
		saved_data_alpha = pickle.load(f)
		f.close()

		f = open('results/u%d.pckl'%(timestep))
		saved_data_u = pickle.load(f)
		f.close()

		f = open('results/T%d.pckl'%(timestep))
		saved_data_T = pickle.load(f)
		f.close()

		alphaM.vector()[:] = saved_data_alphaM
		u.vector()[:] = saved_data_u
		alpha.vector()[:] = saved_data_alpha
		Tsol.vector()[:] = saved_data_T
		
		#################################################################
		# SOLVING ADJOINT PROBLEM #######################################

		IntA = assemble((alpha-alphaM)*dx)/assemble(One*dx)

		IS = Solvers.InverseSolver(ltilde, u, alpha, V, dV, W, Beta, dBeta, Z, Tsol, bc_u, IntA)
		
		IS.solve(V, Beta, Beta_0, lb_Beta, ub_Beta, Beta_error)		

		ub_Beta.vector()[:] = Beta.vector()

		t-=dt

		der += assemble((Beta/(inner(grad(alpha),grad(Beta))+EPS))*dx)/assemble(One*dx)
		
		if abs(IntA) < 2*1e-4:
			der = 0.0
		
		####################################################################

	# SOLVING CONTROL PROBLEM ##################################################
	
	step = (der/kk)
	P = steplength(step, ltilde, iterinv)
	ltildenew = updatel(IntA, ltilde, P, step, iterinv)	
	lerror = abs(ltildenew-ltilde)/ltilde
	ltilde = ltildenew
	iterinv += 1
		

print " "
print "ltilde TARGET: "
print ltilde
print " "
print "Total number of iterations: "
print iterinv
