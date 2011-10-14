#!/usr/bin/env python

# integration.py
# Class hierarchy for different numerical integration schemes.

# Quick Reference:
#   Python Scripting for Computational Science
#     Langtangen, 2008 (3 ed.), pp. 384ff

class Integrator:
    """In general, a numerical integration scheme is an
    approximating sum of the form
    
      int_{-1}^1 f(x)dx = sum_{i=1}^n w_i f(x_i).
    

    Usage:
      rule = Integrator()
      integral = rule.eval(lambda x: x**3)
    
    """
    
    def __init__(self):
        self.setup()

    def setup(self):
        # to be overridden in subclasses
        self.weights = None
        self.points  = None

    def eval(self, f):
        sum = 0.0
        for i in range(len(self.points)):
            sum += self.weights[i]*f(self.points[i])
        return sum

class Trapezoidal(Integrator):
    """The trapezoidal rule is
    
      int_{-1}^1 f(x)dx ~ f(-1) + f(1).
    
    """
    
    def setup(self):
        self.weights =  (1, 1)
        self.points  = (-1, 1)

class Simpson(Integrator):
    """The Simpson rule is
    
      int_{-1}^1 f(x)dx ~ [f(-1) + 4 f(0) + f(1)]/3.
    
    """
    
    def setup(self):
        self.weights = (1/3.0, 4/3.0, 1/3.0)
        self.points  = (-1, 0, 1)

class GaussLegendre2(Integrator):
    """The Gauss-Legendre rule is
    
      int_{-1}^1 f(x)dx ~ f(-1/sqrt{3}) + f(1/sqrt{3}).
    
    """
    
    def setup(self):
        p = 1/math.sqrt(3)
        self.weights =  (1, 1)
        self.points  = (-p, p)


class TransFunc:
    """The transition function for the transformation
    
      int_{Omega_j} f(x)dx = int_{-1}^1 g(xi) (h/2) dxi,
    
    is given by g(xi) = f(x(xi)), with
    
      x(xi) = a + (j-1/2)h + (h/2)xi.
    
    
    """
    def __init__(self, f, h, a):
        self.f = f
        self.h = h
        self.a = a
        # Note there is a hidden parameter j labelling
        # which interval.  This appears in coor_mapping().

    def coor_mapping(self, xi):
        """Map local xi in (-1,1) in interval j to global x."""

        # Note that we use self.j without having initialized it
        return self.a + (self.j - 0.5)*self.h + 0.5*self.h*xi

    def __call__(self, xi):
        x = self.coor_mapping(xi)
        return self.f(x)

def integrate(integrator, a, b, f, n):
    """To integrate over [a,b], we may subdivide into n non-overlapping
    intervals Omega_j and transform each \Omega_j to [-1,1]; we then
    integrate and sum:
    
      int_a^b f(x)dx = sum_{j=1}^n int_{Omega_j} f(x) dx.
    
    We have Omega_j = [(j-1)h, jh], with h = (b-a)/n.

    The main trick is the transformation, writing
    
      int_{Omega_j} f(x)dx = int_{-1}^1 g(xi) (h/2) dxi,
    
    with g(xi) = f(x(xi)).  The class TransFunc encapsulates
    this transition function g.
    
    """
    
    # integrator is an instance of a subclass of Integrator
    sum = 0.0
    h = (b-a)/float(n)
    g = TransFunc(f, h, a)
    for j in range(1, n+1):
        g.j = j
        sum += integrator.eval(g)
    return 0.5*h*sum


if __name__ == '__main__':
    def f(x):
        return x**2

    integrator = Simpson()
    value      = integrate(integrator, 0, 1, f, 100)
    # value should be 1**3 / 3 = 0.333
    print("The integral of f(x) from 0 to 1 is %g" % (value))
