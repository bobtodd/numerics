#!/usr/bin/env python

# A small procedural program to run the binomial pricing algorithm
# for option pricing

import math as m
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

total_t = 1
n_steps = 32
dt = float(total_t) / n_steps

times = np.linspace(0, total_t, n_steps)

rate = 0.06
sigma = 0.3
S_zero = 5
strike = 10

e_minus_rdt = m.exp(-rate*dt)
e_rdt = 1/e_minus_rdt
e_sigma2dt = m.exp( (sigma**2) * dt )

beta = 0.5 * ( e_minus_rdt + e_rdt * e_sigma2dt )
up = beta + m.sqrt( beta**2 - 1 )
down = beta - m.sqrt( beta**2 - 1 )

p = ( e_rdt - down) / (up - down)

S = np.zeros( (n_steps, n_steps) )
for i in range(n_steps):
    for j in range(i):
        k = i - j
        S[j,i] = S_zero * up**j * down**k

V = np.zeros( (n_steps, n_steps) )
for j in range(n_steps):
    V[j, n_steps-1] = max([strike - S[j,n_steps-1], 0])
    
for j in range(n_steps-2, -1, -1):
    for i in range(n_steps-2, -1, -1):
        V[j,i] = e_minus_rdt * ( p*V[j+1,i+1] + (1-p)*V[j,i+1] )
    
print "Put option with"
print "\tr = %f\n\tsigma = %f\n\tS_0 = %f\n\tK = %f" %(rate,sigma,S_zero,strike)
print "Fair option price V[0,0] = %f" % V[0,0]

fig = plt.figure()
ax = fig.gca(projection='3d')

x = []
y = []
z = []
for i in range(len(times)):
    for j in range(n_steps):
        if S[j,i] != 0:
            x.append(times[i])
            y.append(S[j,i])
            z.append(V[j,i])

ax.scatter(x, y, 0, zdir='z', c='b', label='S[j,i]')
ax.scatter(x, y, z, zdir='z', c='r', label='V[j,i]')
ax.mouse_init()
ax.set_xlim3d(0,1)
ax.set_ylim3d(0,max(y))
ax.set_zlim3d(0,max(z))
ax.set_xlabel('Time')
ax.set_ylabel('Stock Price')
ax.set_zlabel('Option Value')

plt.show()
