#!/usr/bin/env python

# implied_vol.py
# A program for calculating the implied volatility.

import math as m
from scipy.stats import norm
import matplotlib.pyplot as plt
import newton

class Stock:
    def __init__(self, S, sigma, delta):
        self.spot     = S
        self.sigma    = sigma
        self.dividend = delta

class EuropeanOption:
    def __init__(self, stock, K, tau, r):
        self.underlying  = stock
        self.strike      = float(K)
        self.time2mature = float(tau)
        self.interest    = float(r)

    def parameters(self):
        d1  = self.interest - self.underlying.dividend
        d1 += self.underlying.sigma**2 / 2
        d1 *= self.time2mature
        d1 += m.log(self.underlying.spot/self.strike)
        d1 /= self.underlying.sigma * m.sqrt(self.time2mature)

        d2 = d1 - self.underlying.sigma * m.sqrt(self.time2mature)

        return d1, d2

    # to be overridden in subclasses
    def valuation(self):
        value = None
        return value
    
    def implied_volatility(self, data=[], plot=False):
        # take a K, V pair from the data
        # guess a sigma
        # create a stock with that sigma
        # create an option with that stock and K
        # do a valuation of the option
        # compare the valuation with V: f := V - valuation
        # apply SecantX() to minimize f
        
        # so f takes sigma as a variable, but has K as parameter
        # f must create the option as part of its routine
        # it also has tau and r as parameters, in order to create the option
        
        sigma0  = self.underlying.sigma
        sigma1  = self.underlying.sigma + 5
        results = []
        for i in range(len(data)):
            K, V         = data[i]
            self.strike  = K
            f            = Comparison(self, V)
            sigma, n, f1 = newton.SecantX(f, sigma0, sigma1)
            results.append((K, sigma))
        
        if plot:
            Ks = [result[0] for result in results]
            sigmas = [result[1] for result in results]
            plt.plot(Ks, sigmas, 'b')
            plt.show()
        
        return results


class EuropeanCall(EuropeanOption):
    def __init__(self, stock, K, tau, r):
        EuropeanOption.__init__(self, stock, K, tau, r)

    def valuation(self):
        d1, d2 = self.parameters()
        
        x = norm()
        
        s1  = x.cdf(d1)
        s1 *= m.exp(-self.underlying.dividend*self.time2mature)
        s1 *= self.underlying.spot

        s2  = x.cdf(d2)
        s2 *= m.exp(-self.interest*self.time2mature)
        s2 *= self.strike

        return s1 - s2

class EuropeanPut(EuropeanOption):
    def __init__(self, stock, K, tau, r):
        EuropeanOption.__init__(self, stock, K, tau, r)

    def valuation(self):
        d1, d2 = self.parameters()
        
        x = norm()
        
        s1  = x.cdf(-d1)
        s1 *= m.exp(-self.underlying.dividend*self.time2mature)
        s1 *= -self.underlying.spot

        s2  = x.cdf(-d2)
        s2 *= m.exp(-self.interest*self.time2mature)
        s2 *= self.strike

        return s1 + s2

class Comparison:
    def __init__(self, option, V):
        self.option = option
        self.V      = V

    def __call__(self, sigma):
        self.option.underlying.sigma = sigma
        return self.V - self.option.valuation()


if __name__ == '__main__':    
    # default values for some stock and option parameters
    sigma_0 = 5
    delta   = 0.95
    K       = 4600
    # parameters specific to the exercise at hand
    tau = 0.211
    S_0 = 5290.36
    r   = 0.0328

    # (K, V) pairs from DAX on call options for Jan 4, 1999
    data = [(6000, 80.2), \
            (6200, 47.1), \
            (6300, 35.9), \
            (6350, 31.3), \
            (6400, 27.7), \
            (6600, 16.6), \
            (6800, 11.4)]

    stock  = Stock(S_0, sigma_0, delta)
    option = EuropeanCall(stock, K, tau, r)
    results = option.implied_volatility(data, plot=True)
    
    print("Resulting pairs (K, sigma):")
    for item in results:
        print("\t" + str(item))
