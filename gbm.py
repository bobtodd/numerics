#!/usr/bin/env python

import math as m
import numpy as np

class Process:
    def __init__(self, n_steps, total_t):
        self.n_steps = int(n_steps)
        self.T       = float(total_t)
        self.X       = np.zeros(self.n_steps)

    def delta_t(self):
        return float(self.T)/float(self.n_steps)

class Stock(Process):
    def __init__(self, n_steps, total_t, rate_rtn, volatility):
        Process.__init__(self, n_steps, total_t)
        self.mu    = rate_rtn
        self.sigma = volatility

class Wiener(Process):
    def __init__(self, n_steps, total_t):
        Process.__init__(self, n_steps, total_t)
        self.setup()

    def setup(self):
        dt = self.delta_t()
        self.X[0] = 0.0
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
        t  = np.linspace(0, self.S.total_t, self.S.n_steps)

        dW = np.zeros(self.S.n_steps)
        for i in range(1, len(dW)):
            dW[i] = self.W[i] - self.W[i-1]
        
        for i in range(self.S.n_steps - 1):
            self.S.X[i+1]  = self.S.X[i]
            self.S.X[i+1] += self.a(self.S.X[i], t[i]) * dt
            self.S.X[i+1] += self.b(self.S.X[i], t[i]) * dW[i]

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
        a = DriftRate(stock)
        b = Bgbm(stock)
        Euler.__init__(self, stock, a, b, wiener)
