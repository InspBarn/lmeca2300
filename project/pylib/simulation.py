# *-* coding: utf-8 *-*

"""
Created on  mer 26 mai 2021 11:57:03 CEST

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

fig = plt.figure(nfig,figsize=(6.4,6.4))

sim = iter(chs.simulate(c0,n,dt))
for i in range(int(tmax//dt)):
    c = next(sim)
    if i*dt in [0,4e-3,8e-3,12e-3]:
        im = plt.imshow(c,cmap='binary',vmin=-1,vmax=1,interpolation='bicubic')
        # fig.colorbar(im, format='%.1f')
        # plt.title(r'Time: $t = %.2e$ [s]' %(i*dt))
        plt.xticks([])
        plt.yticks([])
        fig.tight_layout()
        fig.savefig('../figs/sim_%.2e.pdf' %(i*dt))
