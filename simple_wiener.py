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

times = np.linspace(0, n_steps * delta_t, n_steps)

W = [0.0]
for dt in times[1:]:
    z = np.random.normal(0,1)
    W.append(W[-1] + z * m.sqrt(dt))

plt.plot(times, list(W))
plt.show()
