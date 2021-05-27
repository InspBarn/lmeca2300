# *-* coding: utf-8 *-*

"""
Created on  dim 23 mai 2021 19:37:03 CEST 

@author : vekemans

"""

import numpy as np
from numpy.fft import \
    rfft,irfft,fftfreq,\
    rfft2,irfft2,rfftfreq
from time import time as now

from scipy import sparse
from scipy.interpolate import interp2d

import matplotlib
import matplotlib.pyplot as plt
nfig = 1

import cahn_hilliard as ch
import cahn_hilliard_sd as chs

range_n = 2**np.arange(2,11)
time_fd = np.zeros(range_n.shape)
time_rk4 = np.zeros(range_n.shape)
time_bdf_ab = np.zeros(range_n.shape)

time = 1e-3

for i,n in enumerate(range_n):
    print(n)
    c0 = np.ones((n,n))
    c0[:n//4,:] = -1
    c0[3*n//4:,:] = -1
    c0[:,:n//4] = -1
    c0[:,3*n//4:] = -1

    if n<256:
        start = now()
        sim = iter(ch.simulate(c0.flatten(),n,time))
        next(sim); next(sim)
        end = now()
    time_fd[i] = end-start if n<256 else None

    if n<256:
        start = now()
        sim = iter(chs.simulate(c0,n,time,temp_scheme='rk4'))
        next(sim); next(sim)
        end = now()
    time_rk4[i] = end-start if n<256 else None

    start = now()
    sim = iter(chs.simulate(c0,n,time))
    next(sim); next(sim)
    end = now()
    time_bdf_ab[i] = end-start


fig = plt.figure(nfig)
plt.plot(range_n,time_fd,ls='-',marker='o',markersize=8, label='Finite Differences')
plt.plot(range_n,time_rk4,ls=':',marker='v',markersize=8, label='RK4')
plt.plot(range_n,time_bdf_ab,ls='--',marker='^',markersize=8, label='BDF/AB')
# plt.plot(range_n,time_bdf_ab,c='k',marker='o',markersize=8)
# plt.xticks([2,4,6,8,10],labels=[r'$2^{%d}$' %i for i in [2,4,6,8,10]])
plt.legend(loc='best'); plt.xscale('log', base=2); plt.yscale('log')
plt.xlabel('n'); plt.ylabel(r'Integration times for 10k iterations')
plt.grid(which='major',axis='both',color='xkcd:grey',linestyle='--',linewidth=0.8)
plt.grid(which='minor',axis='y',color='xkcd:grey',linestyle='--',linewidth=0.8)
fig.tight_layout()

fig.savefig('../figs/time.pdf')

plt.show()
