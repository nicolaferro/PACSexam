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

if len(sys.argv) != 2:
	print "\n Incorrect number of input parameters. Usage: \n python direct.py l \n l = INTERNAL LENGTH \n"
	exit(0)

print "********** EXECUTING DIRECT SIMULATION **********"

# Parsing the command line input parameters
l = float(sys.argv[1])

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

# Finite Element Spaces
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
# TEMPERATURE PROBLEM: VARIABLES ###################################################

Tsol, T, w = Function(V_alpha), TrialFunction(V_alpha), TestFunction(V_alpha)
T_old = interpolate(Expression("0.0"), V_alpha)

####################################################################################
# DIRECT PROBLEM: VARIABLES ########################################################

uM, duM, vM = Function(V_u), TrialFunction(V_u), TestFunction(V_u)
alphaM, dalphaM, thetaM = Function(V_alpha), TrialFunction(V_alpha), TestFunction(V_alpha)

######################################################################################
# SOLVING MODEL PROBLEM: MEASURED FIELD ##############################################

print '*******************************************************************************'
print '*                                                                             *'
print '*                              l = %4f                                   *' %l
print '*                                                                             *'
print '*******************************************************************************'

# Error variable
alphaM_error = Function(V_alpha)

# Inital conditions
alphaM_0 = interpolate(Expression("0.0"), V_alpha)
lbM = Function(interpolate(Constant(0.),V_alpha))
ubM = Function(interpolate(Constant(1.),V_alpha))
t = dt
timestep = 0

# Initialize the temperature solver
TS = Solvers.TemperatureSolver(T, w, T_old, Gamma_T_2)

# Initialize the direct solver
DS = Solvers.DirectSolver(l, uM, alphaM, vM, thetaM, duM, dalphaM, Tsol, bc_u, bc_alpha)

while t<=tf:

	timestep +=1
	
	TS.solve(Tsol)

	DS.solve(uM, alphaM, alphaM_0, lbM, ubM, alphaM_error)

	lbM.vector()[:] = alphaM.vector()
	T_old.assign(Tsol)

	f = open('results/alphaM%d.pckl'%(timestep), 'w')
    	pickle.dump(alphaM.vector().array(), f)
	f.close()

	t+=dt


print "********** END OF THE DIRECT SIMULATION **********"
