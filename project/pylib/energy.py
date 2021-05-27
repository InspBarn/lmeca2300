# *-* coding: utf-8 *-*

"""
Created on  mer 26 mai 2021 10:35:15 CEST

@author : vekemans

"""

import numpy as np

import matplotlib
import matplotlib.pyplot as plt
nfig = 1

import cahn_hilliard_sd as chs

tmax = 2e-2
dt = 1e-6
n = 128

c0 = np.random.random(n**2)
c0 = c0 - c0.mean()
c0 = c0.reshape((n,n))

e = np.zeros(int(tmax//dt))

sim = iter(chs.simulate(c0,n,dt))
for i in range(int(tmax//dt)):
    c = next(sim)
    e[i] = (0.25*c**4 - 0.5*c**2 + 1/4).flatten().mean()

fig = plt.figure(nfig)
plt.plot(np.linspace(0,tmax,int(tmax//dt)),e)
plt.ylabel('total free energy')
plt.xlabel('time')
# plt.yscale('log')
# plt.xscale('log')
# plt.xlim([1e-5,3e-2])
plt.grid(which='both')
fig.tight_layout()
nfig += 1

fig2 = plt.figure(nfig)
x = np.linspace(-2,2,100)
plt.plot(x,x**4/4 - x**2/2 + 1/4, c='k')
plt.ylabel(r'total free energy $\mathcal{E}$')
plt.xlabel(r'concentration $c$')
plt.grid()
fig2.tight_layout()

fig.savefig('../figs/energy_time.pdf')
fig2.savefig('../figs/energy_concentration.pdf')

plt.show()
