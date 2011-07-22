#!/usr/bin/env python

# A short script to simulate an SDE via
# the Euler discretization scheme

import sys
import math as m
import numpy as np
import matplotlib.pyplot as plt

def a(y, t):
    return t*y

def b(y, t):
    return m.sqrt(t) * m.sin(y)

n_steps = 100
delta_t = 0.01
x_0     = 0.5

while len(sys.argv) > 1:
    option = sys.argv[1]
    del sys.argv[1]

    if option == '-n':
        n_steps = int(sys.argv[1])
        del sys.argv[1]
    elif option == '-dt':
        delta_t = float(sys.argv[1])
        del sys.argv[1]
    elif option == '-x':
        x_0 == float(sys.argv[1])
        del sys.argv[1]
    else:
        print sys.argv[0], ': Invalid option', option
        sys.exit(1)

times = np.linspace(0, n_steps * delta_t, n_steps)

y = [x_0]
for t in times[1:]:
    dt = delta_t
    Z = np.random.normal(0,1)
    dW = Z * m.sqrt(dt)
    y.append(y[-1] + a(y[-1],t) * dt + b(y[-1],t) * dW)

plt.plot(times, y)
plt.show()
