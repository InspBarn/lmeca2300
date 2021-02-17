# *-* coding: utf-8 *-*

"""
Created on ven 12 fév 2021 13:18:58 UTC

@author: vekemans

"""

import math as mt
import numpy as np
pi = mt.pi

import matplotlib.pyplot as plt
nfig = 1

def FFT(u,x,k):
	v = np.zeros(k.shape,dtype='complex')
	for (i,ki) in enumerate(k):
		v[i] = np.trapz(u*np.exp(-1j*ki*x),x)
	return v

def backFFT(v,x,k):
	u = np.zeros(x.shape,dtype='complex')
	for (i,xi) in enumerate(x):
		u[i] = np.trapz(v*np.exp(1j*xi*k),k) * 1/(2*pi)
	return u

h = [1e-3,1,2]
x = []
f = []
k = np.arange(-10,10,1e-2)
F = []
for (i,hi) in enumerate(h):
	xi = np.arange(-10,10,hi)
	fi = np.zeros(xi.shape)
	fi[np.where(np.logical_and(xi>=-5,xi<=5))] = 1

	Fi = FFT(fi,xi,k)

	x.append(xi)
	f.append(fi)
	F.append(Fi)

plt.figure(nfig)
for (xi,fi,hi) in zip(x,f,h):
	plt.plot(xi,fi, label='h=%.3f'%hi)

plt.xlabel('$x$')
plt.title('Pulse function')
plt.legend(loc='best')
nfig += 1

plt.figure(nfig)
for (Fi,hi) in zip(F,h):
	plt.plot(k,np.real(Fi), label='h=%.3f'%hi)

plt.xlabel('$k$')
plt.title('Real part of the Fourier Transform of the pulse function')
plt.legend(loc='best')
nfig += 1

plt.show()
