# *-* coding: utf-8 *-*

"""
Created on mer 17 fév 2021 09:54:10 UTC

@author: vekemans

"""

import math as mt
import numpy as np
import scipy.linalg as lin
pi = mt.pi

def pulse(x,u,xbeg,xend):
	f = np.zeros(x.shape)
	f[np.where(np.logical_and(x>=xbeg,x<=xend))] = u
	return f

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

def spectral(x,xc,h,der=0):
	u = np.zeros(x.shape)
	S = lambda x: mt.sin(pi*(x-xc)/h) / (pi*(x-xc)/h)

	if der==0:
		for (i,xi) in enumerate(x):
			if xi==xc:
				u[i] = 1
			else:
				u[i] = S(xi)
	elif der==1:
		for (i,xi) in enumerate(x):
			if xi == xc:
				u[i] = 0
			else:
				u[i] = np.cos(pi*(xi-xc)/h)/(xi-xc) - S(xi)/(xi-xc)
	else:
		raise Exception("We don't know the %dth derivative of the spectral function. Please choose whether 0 or 1." %der)

	return u

def spectral_derivative(N):
	OnetoN = np.arange(1,N)
	c = np.append([0], np.power(-1,OnetoN)/OnetoN)
	r = -1*c

	return lin.toeplitz(c,r)
