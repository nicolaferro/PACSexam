from dolfin import *

class AlphaMProblem(OptimisationProblem):

	##
	# @brief Class Constructor
       	# @param TE Objective Function
       	# @param EA First derivative of the Objective Function
       	# @param EAA Second derivative of the Objective Function
       	# @param ALP a variable
    	#

    def __init__(self, TE, EA, EAA, ALP):

        OptimisationProblem.__init__(self)
	self.total_energy = TE
	self.E_alpha = EA
	self.E_alpha_alpha = EAA
	self.alphaM = ALP

	##
	# @brief Objective function assembler
       	# @param x the point in which evaluate the function
	#
    
    def f(self, x):

        self.alphaM.vector()[:] = x
        return assemble(self.total_energy)

	##
	# @brief Gradient of the objective function assembler
	# @param b the tensor/vector where to store the result
	# @param x the point in which evaluate the gradient	
	#

    def F(self, b, x):

        self.alphaM.vector()[:] = x
        assemble(self.E_alpha, tensor=b)

	##
	# @brief Hessian of the objective function assembler
	# @param A the tensor/vector where to store the result
       	# @param x the point in which evaluate the Hessian
	#

    def J(self, A, x):

        self.alphaM.vector()[:] = x
        assemble(self.E_alpha_alpha, tensor=A)


class BetaProblem(NonlinearProblem):

	##
	# @brief Class Constructor
       	# @param FB Equation to solve (FB = 0)
       	# @param JB Jacobian of the equation
       	# @param BETA a variable
    	#

    def __init__(self, FB, JB, BETA):

        NonlinearProblem.__init__(self)
	self.FBeta = FB
	self.JBeta = JB
	self.Beta = BETA

	##
	# @brief Equation assembler
       	# @param b the tensor/vector where to store the result
       	# @param x the point in which evaluate the function
    	#

    def F(self, b, x):	

        self.Beta.vector()[:] = x
        assemble(self.FBeta, tensor=b)

	##
	# @brief Jacobian of the equation assembler
       	# @param A the tensor/vector where to store the result
       	# @param x the point in which evaluate the Jacobian
    	#

    def J(self, A, x):	

        self.Beta.vector()[:] = x
        assemble(self.JBeta, tensor=A)

