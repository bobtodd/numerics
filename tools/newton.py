#!/usr/bin/env python

# newton.py
# A collection of routines for calculating zeros of functions
# using Newton's method.

# Quick reference:
#   Primer on Scientific Programming with Python
#     Langtangen 2009, p. 249

def Newton(f, x, dfdx, epsilon=1.0E-7, N=100, store=False):
    """Find the (local) zero of f given an initial guess x.

    This routine uses Newton's method, which derives straighforwardly
    from the approximate expression for the derivative of f at
    a point x_(n-1):

      f'(x_(n-1)) ~ [f(x_n) - f(x_(n-1))]/[x_n - x_(n-1)]
    
    Rewrite this as

      x_n ~ x_(n-1) + [f(x_n) - f(x_(n-1))]/f'(x_(n-1)),

    then assume f(x_n) = 0 as desired.  Given the starting point
    x_(n-1), the new guess x_n for the zero of f is
    
      x_n ~ x_(n-1) - f(x_(n-1))/f'(x_(n-1)).
    
    """
    
    f_value = f(x)
    n = 0
    if store: info = [(x, f_value)]
    while abs(f_value) > epsilon and n <= N:
        dfdx_value = float(dfdx(x))
        if abs(dfdx_value) < 1E-14:
            raise ValueError("Newton: f'(%g) = %g" % (x, dfdx_value))

        x = x - f_value/dfdx_value

        n += 1
        f_value = f(x)
        if store: info.append((x, f_value))

    if store:
        return x, info
    else:
        return x, n, f_value

def Secant(f, xmin1, xmin2, epsilon=1.0E-7, N=100, store=False):
    """Modification of Newton's method in the case that
    the derivative is unknown.

    The derivation proceeds as with Newton's method, but in place
    of the unknown value f'(x_(n-1)), we use the approximation

      f'(x_(n-1)) ~ [f(x_(n-1)) - f(x_(n-2))]/[x_(n-1) - x_(n-2)].

    This leads to the formula
    
      x_n ~ x_(n-1) - f(x_(n-1))[x_(n-1) - x_(n-2)]/[f(x_(n-1)) - f(x_(n-2))].

    """

    f1 = f(xmin1)
    f2 = f(xmin2)
    n  = 0
    if store: info = [(xmin2, f2), (xmin1, f1)]
    while abs(f1) > epsilon and n <= N:
        dfdx = float((f1 - f2))/float(xmin1 - xmin2)
        if abs(dfdx) < 1E-14:
            raise ValueError("Secant: f'(%g) = %g" % (xmin1, dfdx))

        x = xmin1 - f1/dfdx

        xmin2 = xmin1
        f2    = f(xmin2)
        xmin1 = x
        f1    = f(xmin1)
        
        if store: info.append((xmin1, f1))

    if store:
        return x, info
    else:
        return x, n, f1


if __name__ == '__main__':
    def f(x):
        return x**2 - x

    def df(x):
        return 2*x - 1

    xnew, infonew = Newton(f, 10.5, df, store=True)
    xsec, infosec = Secant(f, 10.5, 10.6, store=True)

    print("Finding zero via Newton's Method:")
    for i in range(len(infonew)):
        print("\tf(%g) = %g" % infonew[i])

    print("\nFinding zero via the Secant Method:")
    for i in range(len(infosec)):
        print("\tf(%g) = %g" % infosec[i])
