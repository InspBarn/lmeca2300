# *-* coding: utf-8 *-*

"""
Created on dim 07 fév 2021 10:47:28 UTC

@author: vekemans

Convergence of periodic spectral method (compare with p1.py)
"""

import math as mt
import numpy as np
import scipy.linalg as lin
import matplotlib.pyplot as plt
nfig = 1

pi = mt.pi
# -----------------------------------------------
# For various N (even), set up a grid as before..
Nvec = np.arange(2,100,2)
plt.figure(nfig)
for (i,N) in enumerate(Nvec):
	h = 2*pi / N
	x = -pi + np.arange(1,N+1)*h
	u = np.exp(np.sin(x))
	uprime = np.cos(x) * u

	# Construct Spectral Differentiation Matrix
	OnetoN = np.arange(1,N)
	c = np.append([0], .5*np.power(-1.,OnetoN)/np.tan(OnetoN*h/2))
	r = c[np.append([0], np.arange(N-1,0,-1))]
	D = lin.toeplitz(c,r)

	error = max(np.abs(D@u - uprime))
	plt.loglog(N, error, 'ko')

	# D = lin.toeplitz(c,-c)
	# error = max(np.abs(D@u - uprime))
	# plt.loglog(N, error, 'o', c='tab:pink', ms=3)

plt.grid(which='both', linestyle='--', linewidth=.5)
plt.xlabel('$N$')
plt.ylabel('Error')
plt.title('Convergence of the spectral')
nfig += 1

# -----------------------------------------------
plt.show()
