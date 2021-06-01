# *-* coding: utf-8 *-*

"""
Created on  sam 22 mai 2021 19:11:56 CEST

@author : vekemans

"""

import numpy as np
from numpy.fft import \
    rfft,irfft,fftfreq,\
    rfft2,irfft2,rfftfreq

from scipy import sparse
from scipy.interpolate import interp2d

import matplotlib
import matplotlib.pyplot as plt
nfig = 1

import cahn_hilliard as ch
import cahn_hilliard_sd as chs

range_n = 2**np.arange(2,11)
error_fd = np.zeros(range_n.shape)
error_rk4 = np.zeros(range_n.shape)
error_bdf_ab = np.zeros(range_n.shape)

n_ref = 2*range_n[-1]
time = 1e-4

# c0 = np.random.random(n_ref*n_ref)*2 -1
# c0 = c0.reshape((n_ref,n_ref))
c0 = np.ones((n_ref,n_ref))
c0[:n_ref//4,:] = -1
c0[3*n_ref//4:,:] = -1
c0[:,:n_ref//4] = -1
c0[:,3*n_ref//4:] = -1

sim = iter(chs.simulate(c0,n_ref,time))
next(sim); sol = next(sim)

x = np.arange(0,1,step=1/n_ref)
# c0_interpolation = interp2d(x,x,c0, kind='cubic')
sol_interpolation = interp2d(x,x,sol, kind='linear')

for i,n in enumerate(range_n):
    x = np.arange(0,1,step=1/n)[:n]

    # c0 = c0_interpolation(x,x)
    c0 = np.ones((n,n))
    c0[:n//4,:] = -1
    c0[3*n//4:,:] = -1
    c0[:,:n//4] = -1
    c0[:,3*n//4:] = -1

    sim = iter(ch.simulate(c0.flatten(),n,time))
    next(sim); guess = next(sim)
    error_fd[i] = ((guess-sol_interpolation(x,x))**2).mean()**0.5

    sim = iter(chs.simulate(c0,n,time,temp_scheme='rk4'))
    next(sim); guess = next(sim)
    error_rk4[i] = ((guess-sol_interpolation(x,x))**2).mean()**0.5

    sim = iter(chs.simulate(c0,n,time))
    next(sim); guess = next(sim)
    error_bdf_ab[i] = ((guess-sol_interpolation(x,x))**2).mean()**0.5


fig = plt.figure(nfig)
plt.plot(range_n,error_fd,ls='-',marker='o',markersize=8, label='FD EE')
plt.plot(range_n,error_rk4,ls=':',marker='v',markersize=8, label='Spectral RK4')
plt.plot(range_n,error_bdf_ab,ls='--',marker='^',markersize=8, label='Spectral BDF/AB')
# plt.xticks([2,4,6,8,10],labels=[r'$2^{%d}$' %i for i in [2,4,6,8,10]])
plt.legend(loc='best'); plt.xscale('log', base=2); plt.yscale('log')
plt.xlabel('n'); plt.ylabel(r'L2 norm error at $t=10^{%d}$ [s]' %(round(np.log(time)/np.log(10))))
plt.grid(which='major',axis='both',color='xkcd:grey',linestyle='--',linewidth=0.8)
plt.grid(which='minor',axis='y',color='xkcd:grey',linestyle='--',linewidth=0.8)
fig.tight_layout()

print(np.log(error_bdf_ab[-2]/error_bdf_ab[-1])/np.log(2))

fig.savefig('../figs/convergence.pdf')

plt.show()
