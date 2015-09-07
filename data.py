from dolfin import *

# Algorithm constants: 
# toll: tolerance between iterations in the fixed-point schemes
# maxiter: maximum number of iterations in the fixed-point schemes
toll = 10e-3
maxiter = 100

# Geometry
# The dimensions of the domain
L = 50.;
H = 50.;

# Material properties
# Elastic constants and damage-model-related values
l0 = Constant(0.3);
E = Constant(340.);
nu = 0.22;
Gc = Constant(42.47);
c_w = Constant(2./3.)
mu = E/(2.0*(1.0 + nu))
Lambda = E*nu/(1.0 - nu**2)

# Thermal properties for the heat equation
T0 = 0.;
deltaT = Constant(380.)
beta = sqrt(Gc/(l0*E*deltaT**2))
kc = Constant(1.)

# Time setting
t0 = 0.;
dt = .2;
tf = 3.4;

# Mesh generation
n = 100;
Th = RectangleMesh(0, 0, L, H, n, n);

# Zero and Unit vectors
# Useful definitions
ndim = 2;
zerov = Constant((0.,)*ndim)
e1 = [Constant([1.,0.]),Constant((1.,0.,0.))][ndim-2]
