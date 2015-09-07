from dolfin import *
from utilities import *
from energies import *
from data import *
import numpy as np
import MyProblems

class TemperatureSolver:

	##
	# @brief Class Constructor
	# It takes in input the temperature variables and the border conditions 
	#
	
	def __init__(self, T, w, T_old, Gamma):

		aT = T*w*dx + dt*inner(kc*0.5*grad(T), grad(w))*dx
		self.LT = (T_old)*w*dx - dt*inner(kc*0.5*grad(T_old), grad(w))*dx

		self.AT = assemble(aT)
		self.Gamma_T = Gamma

	##
	# @brief Solve method
	# It assembles the right hand side of the equation and calls the linear solver
	# @param Tsol the solution function
	#

	def solve(self,Tsol):
		bT = assemble(self.LT)
		self.Gamma_T.apply(self.AT, bT)
		solve(self.AT, Tsol.vector(), bT)


class DirectSolver:

	##
	# @brief Class Constructor
	# It takes in input the primal variables and the border conditions. It instantiates both the solvers for u and alpha
	#

	def __init__(self, l, uM, alphaM, vM, thetaM, duM, dalphaM, Tsol, bc_u, bc_alpha):	
	
		self.total_energyM = total(uM, alphaM, Tsol, mu, Lambda, l, Gc, c_w, beta, T0)
		self.E_U = E_duM(uM, vM, duM, alphaM, Tsol, mu, Lambda, l, Gc, c_w, beta, T0)
		self.E_alphaM = E_dalphaM(uM, alphaM, thetaM, Tsol, mu, Lambda, l, Gc, c_w, beta, T0)
		self.E_alpha_alphaM = E_dalpha_dalphaM(uM, alphaM, dalphaM, thetaM, Tsol, mu, Lambda, l, Gc, c_w, beta, T0)

		a1M = lhs(self.E_U)
		b1M = rhs(self.E_U)
		
		self.problem_uM = LinearVariationalProblem(a1M, b1M, uM, bc_u)
				                       
		self.solver_uM = LinearVariationalSolver(self.problem_uM)

		self.solver_alphaM = PETScTAOSolver()

		# Set some parameters
		self.solver_alphaM.parameters["method"] = "tron"
		self.solver_alphaM.parameters["line_search"] = "gpcg"
		self.solver_alphaM.parameters["monitor_convergence"] = True
		self.solver_alphaM.parameters["report"] = False

		self.BCU = bc_u
		self.BCA = bc_alpha
		

	##
	# @brief solve: Solve method
	# It solves the coupled system with a fixed point scheme
	# @param uM the solution function for the displacements
	# @param alphaM the solution function for the damage
	# @param alphaM_0 the backup function for the previous iteration
	# @param lbM lower bound
	# @param ubM upper bound
	# @param alphaM_error error variable
	#

	def solve(self, uM, alphaM, alphaM_0, lbM, ubM, alphaM_error):
		
		err_alphaM = 1.
		iter = 1

		while err_alphaM>toll and iter<maxiter:
		
			# solve displacement problem
			self.solver_uM.solve()

			# solve damage problem
			self.solver_alphaM.solve(MyProblems.AlphaMProblem(self.total_energyM, self.E_alphaM, self.E_alpha_alphaM, alphaM), alphaM.vector(), lbM.vector(), ubM.vector())

			# test error
			alphaM_error.vector()[:] = alphaM.vector() - alphaM_0.vector()

			# Infinity norm
	       		err_alphaM = np.linalg.norm(alphaM_error.vector().array(), ord = np.Inf)

			# update iteration
			alphaM_0.assign(alphaM)
			iter=iter+1


class InverseSolver:

	##
	# @brief Class Constructor
	# It takes in input the dual variables and the border conditions. It instantiates both the solvers for V and Beta
	#

	def __init__(self, l, u, alpha, V, dV, W, Beta, dBeta, Z, Tsol, bc_u, IntA):
	
		eth = beta*(Tsol-T0)*Identity(2)
		FV = inner(sigma(eps(dV), alpha, mu, Lambda), eps(W))*dx + 2*(alpha-1)*Beta*inner(Stiff(eps(u) - eth, mu, Lambda), eps(W))*dx

		aV = lhs(FV)
		bV = rhs(FV)

		problem_V = LinearVariationalProblem(aV, bV, V, bc_u)
		self.solver_V = LinearVariationalSolver(problem_V)

		self.FBeta = -(IntA)*Z*dx + 2*(alpha-1)*inner(Stiff(eps(u)-eth, mu, Lambda), eps(V))*Z*dx + Beta*inner(Stiff(eps(u)-eth, mu, Lambda), (eps(u)-eth))*Z*dx + Gc/(4.*c_w)*l*inner(grad(Beta), grad(Z))*dx		
		self.JBeta = derivative(self.FBeta, Beta, dBeta)

		# Create the PETScTAOSolver
		self.solver_Beta = PETScSNESSolver()

		# Set some parameters
		self.solver_Beta.parameters["method"] = "vinewtonrsls"
		self.solver_Beta.parameters["report"] = False
		self.solver_Beta.parameters["maximum_iterations"]=800;

	##
	# @brief Solve method
	# It solves the coupled system with a fixed point scheme
	# @param v the solution function for the dual displacements
	# @param Beta the solution function for the dual damage
	# @param Beta_0 the backup function for the previous iteration
	# @param lb_Beta lower bound
	# @param ub_Beta upper bound
	# @param Beta_error error variable
	#	

	def solve(self, V, Beta, Beta_0, lb_Beta, ub_Beta, Beta_error):

		err_Beta = 1.
		iter = 1

		while err_Beta>toll and iter<maxiter:

			# Displacement problem
			self.solver_V.solve()
			
			# solve damage problem
			self.solver_Beta.solve(MyProblems.BetaProblem(self.FBeta, self.JBeta, Beta), Beta.vector(), lb_Beta.vector(), ub_Beta.vector())

			# test error
			Beta_error.vector()[:] = Beta.vector() - Beta_0.vector()

			# Norm
	       		err_Beta = np.linalg.norm(Beta_error.vector().array(), ord = np.Inf)

			# update iteration
			Beta_0.assign(Beta)
			iter=iter+1

