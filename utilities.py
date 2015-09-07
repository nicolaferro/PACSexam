from dolfin import *
import math

##
# A function that returns the symmetric gradient of a field
# @param u The vector field
#

def eps(u):
    return sym(grad(u))

##
# @brief A function that returns the stiffness tensor
# @param epsilon a tensor field
# @param mu one Lame coefficient
# @param Lambda the other Lame coefficient
#
 

def Stiff(epsilon, mu, Lambda):
   return  2*mu*epsilon + Lambda*tr(epsilon)*Identity(2)

##
# @brief A function that returns the stiffness tensor multiplied by the damage factor
#

def sigma(u, alpha, mu, Lambda):
    return (1-alpha)**2*Stiff(u, mu, Lambda)

##
# @brief A function that returns the steplength
# @param step
# @param ltilde the internal length
# @param iterinv the current iteration number of the inverse procedure
#

def steplength(step, ltilde, iterinv):

	P = 1./250
	if step > 0:
		P = (ltilde**2-1.)/10**(int(math.log10(step))+1)
	if iterinv >=5:
		P /= 4.
	if iterinv >=25:
		P /= 10.
	return P

##
# @brief A function that returns the new internal length
# @param IntA the value of the derivative of the object function
# @param ltilde the internal length
# @param P the steplength
# @param step
# @param iterinv the current iteration number of the inverse procedure
#

def updatel(IntA, ltilde, P, step, iterinv):

	if IntA >= -2*1e-4:

			ltildenew = ltilde**2 - P*(step)

			if ltildenew >= 0:
				ltildenew = sqrt(ltildenew)
			else:
				ltildenew = ltilde**2 - P*(step)/5
				ltildenew = sqrt(ltildenew)

	else:
		
			ltildenew = ltilde - 5*IntA*10**int(-math.log10(-IntA))/(iterinv + 2)**(1./3.)

	return ltildenew
