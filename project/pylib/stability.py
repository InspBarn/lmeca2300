# *-* coding: utf-8 *-*

"""
Created on  dim 23 mai 2021 13:46:48 CEST 

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

# import cahn_hilliard as ch
# import cahn_hilliard_sd as chs

a = 0.01
n = 128
c = np.random.random((n,n))*2 - 1
c -= c.mean()
# c /= 10

u = rfftfreq(n, d=1/n).reshape((1,n//2+1)) * 2*np.pi
v = fftfreq(n, d=1/n).reshape((1,n)).T * 2*np.pi
k = np.sqrt(u**2 + v**2)

# vp = -k**2 * rfft(c**2) + k**2 - a**2*k**4
vp = -k**2*rfft(c**2).sum()/2/np.pi + k**2 - a**2*k**4
# vp = -k**2*(rfft(c).sum())**2/2/np.pi + k**2 - a**2*k**4
vp *= 1e-6

fig = plt.figure(nfig)
plt.scatter(vp.real,vp.imag,s=5,c='k')

x = np.linspace(-6,2,100)
y = np.linspace(-4,4,100)
xx,yy = np.meshgrid(x,y)
zz = xx + 1j * yy
rk4 = 1 + zz + zz**2/2 + zz**3/6 + zz**4/24
plt.contour(xx,yy,rk4.imag)

plt.xlabel(r'$\mathcal{R}(\lambda)$')
plt.ylabel(r'$\mathcal{I}(\lambda)$')
plt.grid(which='major',axis='both',color='xkcd:grey',linestyle='--',linewidth=0.8)

fig.tight_layout()
fig.savefig('../figs/stability.pdf')

plt.show()
