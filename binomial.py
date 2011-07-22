#!/usr/bin/env python

# A small program to run the binomial pricing algorithm
# for option pricing

import math as m
import numpy as np

class Stock:
    def __init__(self, r, sigma, S_0, T, steps):
        self.r = r
        self.sigma = sigma
        self.S_0 = S_0
        self.T = T
        self.steps = steps
        self.dt = float(T) / float(steps)
        self.S = np.zeros( (steps, steps) )

    def __call__(self, i, j):
        return self.S[i,j]

    def n_steps(self):
        return self.steps

    def rate(self):
        return self.r

    def delta_t(self):
        return self.dt

    def spot(self, i, j):
        return self.S[i,j]

    def beta(self):
        e_minus_rdt = m.exp(-self.r * self.dt)
        e_rdt = 1/e_minus_rdt
        e_sigma2dt = m.exp( (self.sigma**2) * self.dt )
        return (e_minus_rdt + e_rdt*e_sigma2dt) / 2

    def up(self):
        b = self.beta()
        return b + m.sqrt(b**2 - 1)

    def down(self):
        b = self.beta()
        return b - m.sqrt(b**2 - 1)

    def p(self):
        e_rdt = m.exp(self.r * self.dt)
        u = self.up()
        d = self.down()
        return (e_rdt - d)/(u - d)
    
    # the stock needs some way to evolve...
    def create_tree(self):
        u = self.up()
        d = self.down()
        for i in range(self.steps):
            for j in range(i):
                k = i - j
                self.S[j,i] = self.S_0 * u**j * d**k
        


class Option:
    def __init__(self, strike, underlying):
        self.strike = strike
        self.underlying = underlying
        self.V = np.zeros( (self.underlying.steps, self.underlying.steps) )

    def __call__(self, i, j):
        return self.V[i,j]

    def n_steps(self):
        return self.underlying.n_steps()

    def rate(self):
        return self.underlying.rate()

    def delta_t(self):
        return self.underlying.delta_t()

    def fair_price(self):
        return self.V[0,0]



class EuropeanPut(Option):
    def __init__(self, strike, underlying):
        Option.__init__(self, strike, underlying)
    
    def payoff(self, j):
        # calculate the payout
        # depending on the price of the underlying
        return max([self.strike - self.underlying(j, self.underlying.n_steps()-1), 0])
    
    # evolve backwards
    def create_tree(self):
        n_steps = self.n_steps()
        e_minus_rdt = m.exp(-self.rate() * self.delta_t())
        p = self.underlying.p()
        
        for j in range(n_steps):
            self.V[j, n_steps-1] = self.payoff(j)
        
        for j in range(n_steps-2, -1, -1):
            for i in range(n_steps-2, -1, -1):
                self.V[j,i] = e_minus_rdt * ( p*self.V[j+1,i+1] + (1-p)*self.V[j,i+1] )
        


class EuropeanCall(Option):
    def __init__(self, strike, underlying):
        Option.__init__(self, strike, underlying)
    
    def payoff(self, j):
        # calculate the payout
        # depending on the price of the underlying
        return max([self.underlying(j, self.underlying.n_steps()-1) - self.strike, 0])

    # evolve backwards
    def create_tree(self):
        n_steps = self.n_steps()
        e_minus_rdt = m.exp(-self.rate() * self.delta_t())
        p = self.underlying.p()
        
        for j in range(n_steps):
            self.V[j, n_steps-1] = self.payoff(j)
        
        for j in range(n_steps-2, -1, -1):
            for i in range(n_steps-2, -1, -1):
                self.V[j,i] = e_minus_rdt * ( p*self.V[j+1,i+1] + (1-p)*self.V[j,i+1] )




class AmericanPut(Option):
    def __init__(self, strike, underlying):
        Option.__init__(self, strike, underlying)
    
    def payoff(self, j):
        # calculate the payout
        # depending on the price of the underlying
        return max([self.strike - self.underlying(j, self.underlying.n_steps()-1), 0])
    
    # evolve backwards
    def create_tree(self):
        n_steps = self.n_steps()
        e_minus_rdt = m.exp(-self.rate() * self.delta_t())
        p = self.underlying.p()
        
        for j in range(n_steps):
            self.V[j, n_steps-1] = self.payoff(j)
        
        for j in range(n_steps-2, -1, -1):
            for i in range(n_steps-2, -1, -1):
                early_exercise   = max([self.underlying(j,i) - self.strike, 0])
                calculated_value = e_minus_rdt * ( p*self.V[j+1,i+1] + (1-p)*self.V[j,i+1] )
                self.V[j,i]      = max([early_exercise, calculated_value])
        


class AmericanCall(Option):
    def __init__(self, strike, underlying):
        Option.__init__(self, strike, underlying)
    
    def payoff(self, j):
        # calculate the payout
        # depending on the price of the underlying
        return max([self.underlying(j, self.underlying.n_steps()-1) - self.strike, 0])

    # evolve backwards
    def create_tree(self):
        n_steps = self.n_steps()
        e_minus_rdt = m.exp(-self.rate() * self.delta_t())
        p = self.underlying.p()
        
        for j in range(n_steps):
            self.V[j, n_steps-1] = self.payoff(j)
        
        for j in range(n_steps-2, -1, -1):
            for i in range(n_steps-2, -1, -1):
                early_exercise   = max([self.strike - self.underlying(j,i), 0])
                calculated_value = e_minus_rdt * ( p*self.V[j+1,i+1] + (1-p)*self.V[j,i+1] )
                self.V[j,i]      = max([early_exercise, calculated_value])


if __name__ == '__main__':
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    import sys

    visual   = False
    american = False
    put      = False
    total_t  = 1
    n_steps  = 64
    rate     = 0.06
    sigma    = 0.3
    S_zero   = 5
    strike   = 10
    
    while len(sys.argv) > 1:
        option = sys.argv[1]
        del sys.argv[1]

        if option == '-am':
            american = True
        elif option == '-v':
            visual = True
        elif option == '-p':
            put = True
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
        elif option == '-k':
            strike = float(sys.argv[1])
            del sys.argv[1]
        else:
            print sys.argv[0], ': Invalid option', option
            sys.exit(1)
    
    
    times = np.linspace(0, total_t, n_steps)
    

    S = Stock(rate, sigma, S_zero, total_t, n_steps)
    S.create_tree()

    if put:
        if american:
            V = AmericanPut(strike, S)
        else:
            V = EuropeanPut(strike, S)
    else:
        if american:
            V = AmericanCall(strike, S)
        else:
            V = EuropeanCall(strike, S)

    V.create_tree()
    
    outstr = ""
    if american:
        outstr += "American"
    else:
        outstr += "European"

    if put:
        outstr += " put"
    else:
        outstr += " call"

    outstr += " option with"
    print outstr
    print "\tr = %f\n\tsigma = %f\n\tS_0 = %f\n\tK = %f" %(rate,sigma,S_zero,strike)
    print "Fair option price V[0,0] = %f" % V.fair_price()

    if visual:
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        
        x = []
        y = []
        z = []
        for i in range(len(times)):
            for j in range(n_steps):
                if S(j,i) != 0:
                    x.append(times[i])
                    y.append(S(j,i))
                    z.append(V(j,i))
        
        ax.scatter(x, y, 0, zdir='z', c='b', label='S(j,i)')
        ax.scatter(x, y, z, zdir='z', c='r', label='V(j,i)')
        ax.mouse_init()
        ax.set_xlim3d(0,1)
        ax.set_ylim3d(0,max(y))
        ax.set_zlim3d(0,max(z))
        ax.set_xlabel('Time')
        ax.set_ylabel('Stock Price')
        ax.set_zlabel('Option Value')
        
        plt.show()
