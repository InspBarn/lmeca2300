# *-* coding: utf-8 *-*

"""
Created on dim 07 fév 2021 09:27:34 UTC

@author: vekemans

"""

import time
import math as mt
import numpy as np
import scipy.sparse as spr
import matplotlib.pyplot as plt
nfig = 1

pi = mt.pi
# -----------------------------------------------
# Orginal Signal
N = np.power(2,12)
h = 2*pi / N
x = -pi + np.arange(1,N+1)*h

u = np.exp(np.sin(x))
uprime = np.cos(x) * u

plt.figure(nfig)
plt.plot(x,u, 'k--', label=r"$u$")
plt.plot(x,uprime, 'k', label=r"$u'$")
plt.title("Function $e^{sin(x)}$ and its derivative")
plt.legend()
plt.grid()
nfig += 1

# -----------------------------------------------
Nvec = np.power(2, np.arange(3,13))
timevec = np.zeros(Nvec.shape)

plt.figure(nfig)

for (i,N) in enumerate(Nvec):
	t = time.time()
	h = 2*pi / N
	x = -pi + np.arange(1,N+1)*h

	u = np.exp(np.sin(x))
	uprime = np.cos(x)*u

	# Construct the differentiation matrix
	e = np.ones(N)
	row = np.arange(0,N)
	col1 = np.append(np.arange(1,N), [0])
	col2 = np.append(np.arange(2,N), [0,1])
	D = spr.csr_matrix((2.*e/3., (row,col1)), shape=(N,N)) \
		- spr.csr_matrix((e/12., (row,col2)), shape=(N,N))
	D = (D - D.transpose())/h

	# Plot max(abs(D*u - uprime))
	error = max(np.abs(D*u - uprime))
	plt.loglog(N, error, 'ko')
	timevec[i] = time.time() - t

plt.loglog(Nvec, np.float64(Nvec)**(-4), 'k--')
plt.text(100, 1e-10, '$N^{-4}$', fontsize=14)
plt.title("Convergence of fourth-order finite differences")
plt.xlabel('$N$')
plt.ylabel('Error')
plt.grid(which='both', linestyle='--', linewidth=.5)
nfig += 1

# plt.figure(nfig)
# plt.loglog(Nvec, timevec, 'ko')
# plt.title(" ")
# plt.xlabel("$N$")
# plt.ylabel("Time [s]")
# plt.grid()
# nfig += 1

# -----------------------------------------------
plt.show()
