#!/usr/bin/env python

import math as m
import numpy as np

class Process:
    def __init__(self, n_steps, total_t, X_0=0.0):
        self.n_steps = int(n_steps)
        self.T       = float(total_t)
        self.X       = np.zeros(self.n_steps)
        self.X[0]    = float(X_0)

    def __getitem__(self, i):
        return self.X[i]

    def __setitem__(self, i, value):
        self.X[i] = float(value)

    def delta_t(self):
        return float(self.T)/float(self.n_steps)

class Stock(Process):
    def __init__(self, n_steps, total_t, X_0, rate_rtn, volatility):
        Process.__init__(self, n_steps, total_t, X_0)
        self.mu    = rate_rtn
        self.sigma = volatility

class Wiener(Process):
    def __init__(self, n_steps, total_t):
        Process.__init__(self, n_steps, total_t)
        self.setup()

    def setup(self):
        dt = self.delta_t()
        for i in range(self.n_steps - 1):
            Z = np.random.normal(0,1)
            self.X[i+1] = self.X[i] + Z * m.sqrt(dt)

class Euler:
    def __init__(self, stock, a, b, wiener):
        self.S = stock
        self.a = a
        self.b = b
        self.W = wiener
    
    def evolve(self):
        dt = self.S.delta_t()
        t  = np.linspace(0, self.S.T, self.S.n_steps)

        dW = np.zeros(self.S.n_steps)
        for i in range(1, len(dW)):
            dW[i] = self.W[i] - self.W[i-1]
        
        for i in range(self.S.n_steps - 1):
            self.S[i+1]  = self.S[i]
            self.S[i+1] += self.a(self.S[i], t[i]) * dt
            self.S[i+1] += self.b(self.S[i], t[i]) * dW[i]

class DriftRate:
    def __init__(self, stock):
        self.mu = stock.mu

    def __call__(self, x, t):
        return self.mu * x

class Bgbm:
    def __init__(self, stock):
        self.sigma = stock.sigma

    def __call__(self, x, t):
        return self.sigma * x

class GBM(Euler):
    def __init__(self, stock, wiener):
        # Note that, b/c of the form of DriftRate and Bgbm,
        # if stock starts at 0, it stays at 0
        a = DriftRate(stock)
        b = Bgbm(stock)
        Euler.__init__(self, stock, a, b, wiener)


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import sys

    visual   = False
    total_t  = 1
    n_steps  = 10000
    rate     = 0.35
    sigma    = 0.3
    S_zero   = 5
    
    while len(sys.argv) > 1:
        option = sys.argv[1]
        del sys.argv[1]

        if option == '-v':
            visual = True
        elif option == '-n':
            n_steps = int(sys.argv[1])
            del sys.argv[1]
        elif option == '-t':
            total_t = float(sys.argv[1])
            del sys.argv[1]
        elif option == '-r':
            rate = float(sys.argv[1])
            del sys.argv[1]
        elif option == '-sig':
            sigma = float(sys.argv[1])
            del sys.argv[1]
        elif option == '-s':
            S_zero = float(sys.argv[1])
            del sys.argv[1]
        else:
            print sys.argv[0], ': Invalid option', option
            sys.exit(1)
    
        

    S = Stock(n_steps, total_t, S_zero, rate, sigma)
    W = Wiener(n_steps, total_t)

    G = GBM(S, W)
    G.evolve()

    
    if visual:
        t = np.linspace(0, total_t, n_steps)

        plt.plot(t, S[:], 'b', label='Stock Price')
        plt.plot(t, W[:], 'r', label='Brownian Motion')
        plt.legend()
        
        plt.show()
