#!/usr/bin/env python

# A short script to simulate an SDE via
# the Euler discretization scheme

import sys
import math as m
import numpy as np
import matplotlib.pyplot as plt

def a(y, t):
    # a = t * y
    a = 0.05 * y
    return a

def b(y, t):
    # b = m.sqrt(t) * m.sin(y)
    b = 0.3 * y
    return b

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

t        = np.linspace(0, n_steps * delta_t, n_steps)
y        = np.zeros(n_steps)
y_bar    = np.zeros(n_steps)
y[0]     = x_0
y_bar[0] = x_0
for i in range(n_steps-1):
    dt         = delta_t
    Z          = np.random.normal(0,1)
    dW         = Z * m.sqrt(dt)
    y[i+1]     = y[i] + a(y[i],t[i]) * dt + b(y[i],t[i]) * dW
    y_bar[i+1] = y_bar[i] + a(y_bar[i],t[i]) * dt

plt.plot(t, y, 'b')
plt.plot(t, y_bar, 'r')
plt.show()
