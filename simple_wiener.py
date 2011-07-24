#!/usr/bin/env python

# A short script to simulate a Wiener process

import sys
import math as m
import numpy as np
import matplotlib.pyplot as plt

n_steps = 100
delta_t = 0.01

while len(sys.argv) > 1:
    option = sys.argv[1]
    del sys.argv[1]

    if option == '-n':
        n_steps = int(sys.argv[1])
        del sys.argv[1]
    elif option == '-dt':
        delta_t = float(sys.argv[1])
        del sys.argv[1]
    else:
        print sys.argv[0], ': Invalid option', option
        sys.exit(1)

t    = np.linspace(0, n_steps * delta_t, n_steps)
W    = np.zeros(n_steps)
W[0] = 0.0
for i in range(n_steps-1):
    dt     = delta_t
    Z      = np.random.normal(0,1)
    W[i+1] = W[i] + Z * m.sqrt(dt)

plt.plot(t, list(W))
plt.show()
