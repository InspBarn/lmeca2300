# *-* coding: utf-8 *-*

"""
Created on  ven 05 mar 2021 14:40:48 CET 

@author : vekemans

p14.py - Solve non linear BVP u_xx = exp(u), u(-1) = u(1) = 0

"""

from numpy import *
from numpy.linalg import inv,norm

import matplotlib
import matplotlib.pyplot as plt
nfig = 1

N = 16

# D : differentiation matrix
# x : Chebyshev grid
def cheb(N):
	if N==0:
		return 0,1
	x = cos(pi * linspace(0,1, N+1))
	c = concatenate(([2],ones(N-1),[2])) * power(-1, arange(N+1))
	X = x.repeat(N+1).reshape((N+1,N+1))
	dX = X - X.T
	D = (c/c.T) / (dX + eye(N+1))
	D = D - diag(sum(D, axis=1))
	return D,x

D,x = cheb(N)

D2 = D@D
D2 = D2[1:N,1:N]
N2 = len(D2)
u = zeros(N2)
err,n = 1,0
while (err > 1e-15) and (n < 100):
	# fixed-point iteration
	unew = inv(D2)@exp(u)
	err = norm(unew-u, inf)
	u = unew
	n += 1

u = append([0],append(u,[0]))

fig = plt.figure(nfig)
plt.grid(ls='--', lw=.8, c='lightgrey')
plt.scatter(x,u, s=16)
xx = linspace(-1,1, 101)
uu = polyval(polyfit(x,u,N//2),xx)
plt.plot(xx,uu)
plt.title(r'no.steps = %d, $u_0$ = %.5f' %(n,u[N//2]))

plt.show()
