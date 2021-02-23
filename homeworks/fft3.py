# *-* coding: utf-8 *-*

"""
Created on mar 23 fév 2021 09:14:21 UTC

@author: vekemans

"""

import math as mt
import numpy as np

from math import pi as π
from numpy.fft import fft,fftshift,ifft

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation
nfig = 1

def dft(func, start,end,N, order=1):
	"""
	Return the 1D-derivative of function 'func', at specified order
	"""

	L = end-start
	h = L/N
	x = np.arange(start,end, step=h)

	# -- Compute the FFT
	u = func(x)
	u_hat = fftshift(fft(u))
	k = np.arange(-N/2,N/2)

	# -- Compute the Derivative
	u_dot_hat = (1j*k)**order * u_hat
	u_dot = ifft(fftshift(u_dot_hat)).real

	return x, u_dot

def pulse(x, lim=1):
	"""
	Pulse function
	"""

	u = np.zeros(x.shape)
	u[(-lim<=x) & (x<=lim)] = 1

	return u

def gaussian(x, μ=0,σ=1):
	"""
	Gaussian function
	"""
	u = 1/mt.sqrt(2*π*σ**2) * np.exp(-(x-μ)**2/(2*σ**2))
	return u

def dgaussian(x, μ=0,σ=1):
	"""
	Derivative of the gaussian function
	"""
	u_dot = -2*(x-μ)/(2*σ**2) * gaussian(x, μ=μ,σ=σ)
	return u_dot

if __name__=='__main__':
	m = np.arange(2,15)
	Nvec = 2**m
	error = np.empty(m.shape)
	
	μ = π
	σ = 0.5
	f = lambda x: gaussian(x, μ=μ,σ=σ)

	for (i,N) in enumerate(Nvec):
		x,u_dot = dft(f, 0,2*π,N, order=1)
		uprime  = -2*(x-μ)/(2*σ**2) * f(x)
		error[i] = np.abs(uprime-u_dot).max()

	# plt.figure(nfig)
	# plt.plot(x,f(x), 'k-')
	# plt.plot(x,uprime, 'k--')
	# nfig += 1

	# plt.figure(nfig)
	# plt.loglog(Nvec,error, 'k-o', markersize=4, lw=0)
	# nfig += 1

	N = 2**9
	x = np.arange(0,2*π, step=2*π/N)
	u = f(x)
	c = 0.2 + np.power(np.sin(x-1),2)

	tmax = 16*π+1
	tplot = .05
	plotgap = int(tplot/0.001)
	dt = tplot/plotgap
	Nplots = int(tmax/tplot)

	data = np.zeros((Nplots,N))
	data[0,:] = u; u_old = u
	time = np.arange(Nplots)/Nplots * tmax
	for i in range(Nplots):
		for n in range(plotgap):
			u_hat = fftshift(fft(u))
			u_dot_hat = (1j*np.arange(-N/2,N/2)) * u_hat
			u_dot = ifft(fftshift(u_dot_hat)).real

			u_new = u_old - dt*c*u_dot
			u_old = u
			u = u_new

		data[i,:] = u

	fig = plt.figure(nfig)
	axe = fig.add_subplot(1,1,1)

	axe.set_xlim([x[0],x[-1]]); axe.set_xlabel(r'$x$')
	axe.set_ylim([0,1]); axe.set_ylabel(r'$u_t$')
	axe.set_title(r'$u_t + c(x) u_x = 0, \qquad u_t^0 = \mathcal{N}(\pi,0.5)$')

	u_plt, = axe.plot([], [], 'k-')
	u_txt  = axe.text(π, 0.9, '', ha='center', fontsize=10)

	def animate(t):
		u_plt.set_data(x, data[t,:])
		u_txt.set_text('Current time : t = %.2f [s]' %(t*tplot))
		return u_plt,u_txt,

	anim = animation.FuncAnimation(fig, animate, Nplots, interval=100*tplot, blit=True)

# -- Show figures
plt.show()
