from dolfin import *
from utilities import *

##
# A function that returns the total elastic energy
#
def total(uM, alphaM, Tsol, mu, Lambda, l, Gc, c_w, beta, T0):

	eth = beta*(Tsol-T0)*Identity(2)

	elastic_energyM = 0.5*inner(sigma(eps(uM),alphaM, mu, Lambda)-sigma(eth,alphaM, mu, Lambda), eps(uM)-eth)*dx
	dissipated_energyM = Gc/(4.*c_w)*(alphaM/l + l*dot(grad(alphaM), grad(alphaM)))*dx

	return elastic_energyM + dissipated_energyM
	
##
# A function that returns the first derivative of the total energy with respect to u
#
def E_duM(uM, vM, duM, alphaM, Tsol, mu, Lambda, l, Gc, c_w, beta, T0):

	eth = beta*(Tsol-T0)*Identity(2)

	elastic_energyM = 0.5*inner(sigma(eps(uM),alphaM, mu, Lambda)-sigma(eth,alphaM, mu, Lambda), eps(uM)-eth)*dx
	dissipated_energyM = Gc/(4.*c_w)*(alphaM/l + l*dot(grad(alphaM), grad(alphaM)))*dx
	total_energyM = elastic_energyM + dissipated_energyM

	E_uM = derivative(total_energyM,uM,vM)
	return replace(E_uM,{uM:duM})

##
# A function that returns the first derivative of the total energy with respect to alpha
#
def E_dalphaM(uM, alphaM, thetaM, Tsol, mu, Lambda, l, Gc, c_w, beta, T0):

	eth=beta*(Tsol-T0)*Identity(2)

	elastic_energyM = 0.5*inner(sigma(eps(uM),alphaM, mu, Lambda)-sigma(eth,alphaM, mu, Lambda), eps(uM)-eth)*dx
	dissipated_energyM = Gc/(4.*c_w)*(alphaM/l + l*dot(grad(alphaM), grad(alphaM)))*dx
	total_energyM = elastic_energyM + dissipated_energyM

	return derivative(total_energyM, alphaM, thetaM)

##
# A function that returns the second derivative of the total energy with respect to alpha
#
def E_dalpha_dalphaM(uM, alphaM, dalphaM, thetaM, Tsol, mu, Lambda, l, Gc, c_w, beta, T0):

	eth=beta*(Tsol-T0)*Identity(2)

	elastic_energyM = 0.5*inner(sigma(eps(uM),alphaM, mu, Lambda)-sigma(eth,alphaM, mu, Lambda), eps(uM)-eth)*dx
	dissipated_energyM = Gc/(4.*c_w)*(alphaM/l + l*dot(grad(alphaM), grad(alphaM)))*dx
	total_energyM = elastic_energyM + dissipated_energyM

	tmp = derivative(total_energyM, alphaM, thetaM)
	return derivative(tmp, alphaM, dalphaM)



