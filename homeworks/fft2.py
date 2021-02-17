# *-* coding: utf-8 *-*

"""
Created on mer 17 fév 2021 09:46:19 UTC

@author: vekemans

"""

import math as mt
import numpy as np
import scipy.linalg as lin
pi = mt.pi

import matplotlib
import matplotlib.pyplot as plt
nfig = 1

import utils as us

## ----------------------------------------------
N = 50
h = 2*pi / N
x = -pi + np.arange(1,N+1)*h

u = np.exp(np.sin(x))
uprime = np.cos(x) * u

OnetoN = np.arange(1,N)
OnetoNSign = np.power(-1,OnetoN)/OnetoN
c = np.append([0], OnetoNSign)
r = -1*c
S = lin.toeplitz(c,r)/h

gprim1 = S@u
error1 = np.abs(S@u-uprime)

c = np.append([0], OnetoNSign - np.flip(OnetoNSign))
r = -1*c
S = lin.toeplitz(c,r)/h

gprim2 = S@u
error2 = np.abs(S@u-uprime)

## ----------------------------------------------
plt.figure(nfig)
plt.plot(x,u, 'k-', label=r"$u$")
plt.plot(x,uprime, 'k--', label=r"$u'$")
plt.plot(x,gprim1, ls=':', c='tab:red') #, label='non-periodic')
plt.plot(x,gprim2, ls=':', c='tab:green') #, label='periodic')
plt.title("Function $e^{sin(x)}$ and its derivative")
plt.legend(loc='best')
plt.grid()
nfig += 1

plt.figure(nfig)
plt.title('Influence of boundary conditions on the error\nof derivation by spectral method')
plt.loglog(x,error1, ':', c='tab:red', label='non-periodic')
plt.loglog(x,error2, ':', c='tab:green', label='periodic')
plt.legend(loc='best')
plt.grid()
nfig += 1

## ----------------------------------------------
Nvec = np.arange(10,100,2)
error = np.empty(Nvec.shape)
for (i,N) in enumerate(Nvec):
	h = 2*pi / N
	x = -pi + np.arange(1,N+1)*h

	u = np.exp(np.sin(x))
	uprime = np.cos(x) * u

	OnetoN = np.arange(1,N)
	OnetoNSign = np.power(-1,OnetoN)/OnetoN

	c = np.append([0], OnetoNSign - np.flip(OnetoNSign))
	r = -1*c
	S = lin.toeplitz(c,r)/h

	error[i] = np.max(np.abs(S@u-uprime))


## ----------------------------------------------
plt.figure(nfig)
plt.plot(Nvec,error, 'k-')
plt.title("Convergence of spectral method")
plt.xlabel(r'$N$')
plt.ylabel('Error')
plt.grid()
nfig += 1

plt.show()
