# *-* coding: utf-8 *-*

"""
Created on  mer 26 mai 2021 14:11:39 CEST

@author : vekemans

"""

import numpy as np
from time import time

from scipy.interpolate import interp2d

import matplotlib
import matplotlib.pyplot as plt
nfig = 1

import cahn_hilliard_sd as chs

tmax = 12e-3
dt = 1e-5
n1 = 128
n2 = 512

x1 = np.arange(0,1,step=1/n1)
x2 = np.arange(0,1,step=1/n2)

c1 = np.ones((n1,n1))
c1[:n1//4,:] = -1
c1[3*n1//4:,:] = -1
c1[:,:n1//4] = -1
c1[:,3*n1//4:] = -1

c2 = np.ones((n2,n2))
c2[:n2//4,:] = -1
c2[3*n2//4:,:] = -1
c2[:,:n2//4] = -1
c2[:,3*n2//4:] = -1

t = np.zeros(int(tmax//dt))
e = np.zeros(int(tmax//dt))
sim1 = iter(chs.simulate(c1,n1,dt))
sim2 = iter(chs.simulate(c2,n2,dt))

for i in range(int(tmax//dt)):
    start = time()
    c1 = next(sim1)
    end = time()

    c2 = next(sim2)
    interp = interp2d(x2,x2,c2, kind='linear')(x1,x1)
    t[i] = end-start
    e[i] = ((c1-interp)**2).mean()**0.5

    print('%d / %d' %(i,int(tmax//dt)))

fig = plt.figure(nfig)
ax1 = fig.add_subplot()
ax2 = ax1.twinx()
ax1.set_xlabel('time')

ax1.plot(np.linspace(dt,tmax,int(tmax//dt)-1),e[1:],c='tab:blue')
ax1.set_ylabel('L2 norm error', color='tab:blue')
ax1.tick_params(axis='y', labelcolor='tab:blue')
ax1.set_yscale('log')

ax2.plot(np.linspace(dt,tmax,int(tmax//dt)-1),t[1:],c='tab:red')
ax2.set_ylabel('Time of integration', color='tab:red')
ax2.tick_params(axis='y', labelcolor='tab:red')
ax2.set_yscale('log')

fig.tight_layout()
nfig += 1

fig.savefig('../figs/efficiency.pdf')
plt.show()
