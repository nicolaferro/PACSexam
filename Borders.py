from dolfin import *

##
# @brief A class for the left border of the domain
#

class Left(SubDomain):

    def __init__(self, LL, HH):
	self.L = LL
	self.H = HH 
	SubDomain.__init__(self)

    ##
    # @param x the point to evaluate
    # \sa near
    # @return a boolean: True if x is on the border, False otherwise

    def inside(self, x, on_boundary):
        return near(x[0], 0.)

##
# @brief A class for the right border of the domain
#
        
class Right(SubDomain):

    def __init__(self, LL, HH):
	self.L = LL
	self.H = HH
	SubDomain.__init__(self)

    ##
    # @param x the point to evaluate
    # \sa near
    # @return a boolean: True if x is on the border, False otherwise


    def inside(self, x, on_boundary):
        return near(x[0], self.L)

##
# @brief A class for the upper border of the domain
#

class Up(SubDomain):

    def __init__(self, LL, HH):
	self.L = LL
	self.H = HH
	SubDomain.__init__(self) 

    ##
    # @param x the point to evaluate
    # \sa near
    # @return a boolean: True if x is on the border, False otherwise

    def inside(self, x, on_boundary):
        return near(x[1], self.H)

##
# @brief A class for the lower border of the domain
#

class Down(SubDomain):

    def __init__(self, LL, HH):
	self.L = LL
	self.H = HH
	SubDomain.__init__(self) 

    ##
    # @param x the point to evaluate
    # \sa near
    # @return a boolean: True if x is on the border, False otherwise

    def inside(self, x, on_boundary):
        return near(x[1], 0.)
