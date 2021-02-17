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

import utils as us

h = [1e-2,1]
x = []
f = []
k = np.arange(-10,10,1e-2)
F = []
for (i,hi) in enumerate(h):
	xi = np.arange(-10,10,hi)
	fi = us.pulse(xi,1,-5,5)
	Fi = us.FFT(fi,xi,k)
	##
	x.append(xi)
	f.append(fi)
	F.append(Fi)

plt.figure(nfig)
plt.plot(x[0],f[0], 'k-')
plt.xlabel('$x$')
plt.title('Pulse function')
nfig += 1

##
ls = ['k-','k--']

plt.figure(nfig)
for (Fi,hi,lsi) in zip(F,h,ls):
	plt.plot(k,np.real(Fi), lsi, label='h=%.2f'%hi)

plt.xlabel('$k$')
plt.title('Real part of the Fourier Transform of the pulse function')
plt.legend(loc='best')
nfig += 1

plt.show()
