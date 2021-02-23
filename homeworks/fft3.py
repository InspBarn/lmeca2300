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

	fig = plt.figure(nfig)
	plt.plot(x,f(x), 'k-', label=r'$u$')
	plt.plot(x,uprime, 'k--', label=r"$u'$")
	plt.xlabel(r'$x$')
	plt.title('Gaussian function and its derivative')
	plt.legend()
	plt.tight_layout()
	fig.savefig('../figs/gaussian_function.pdf')
	nfig += 1

	fig = plt.figure(nfig)
	plt.loglog(Nvec,error, 'k-o', markersize=4, lw=0)
	plt.xlabel(r'$N$')
	plt.ylabel('Error')
	plt.title('Convergence of spectral method for differentiation')
	plt.grid(which='major', linestyle='--', linewidth=0.5)
	plt.grid(which='minor', linestyle=':', linewidth=0.25)
	plt.tight_layout()
	fig.savefig('../figs/spectral_convergence.pdf')
	nfig += 1

	N = 2**8
	x = np.arange(0,2*π, step=2*π/N)
	# u = gaussian(x, μ=π,σ=σ)
	u = np.exp(-100*(x-1)**2)
	c = 0.2 + np.power(np.sin(x-1),2)

	tmax = 8.
	tplot = .1
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

			if i==0:
				u_new = u_old - dt*c*u_dot
			else:
				u_new = u_old - 2*dt*c*u_dot
			u_old = u
			u = u_new

		data[i,:] = u

	fig = plt.figure(nfig)
	u_plt, = plt.plot([], [], 'k-')
	u_txt  = plt.text(π, 0.9, '', ha='center', fontsize=10)

	plt.xlim([x[0],x[-1]])
	plt.ylim([0,1])
	plt.xlabel(r'$x$', fontsize=10)
	plt.ylabel(r'$u$', fontsize=10)
	plt.title(r'$u_t + c(x) u_x = 0 \qquad u^0 = \mathcal{N}(\pi,0.5)$', fontsize=12)
	plt.tight_layout()
	nfig += 1

	def animate(t):
		u_plt.set_data(x, data[t,:])
		u_txt.set_text('Current time : t = %.2f [s]' %(t*tplot))
		return u_plt,u_txt,

	anim = animation.FuncAnimation(fig, animate, Nplots, interval=100*tplot, blit=True)

	fig = plt.figure(nfig)
	ax  = plt.axes(projection='3d')
	# for i in range(Nplots):
	# 	t = time[i]*np.ones(x.shape)
	# 	ax.plot3D(x,t, data[i,:], color='tab:blue')
	X,T = np.meshgrid(x,time)
	print(x.shape,time.shape,X.shape,data.shape)
	ax.plot_surface(X,T, data, cmap='w')
	ax.set_xlim([x[0],x[-1]])
	ax.set_ylim([0,tmax])
	ax.set_zlim([0,4])
	ax.set_xlabel('X-axis')
	ax.set_ylabel('Time-axis')
	ax.set_zlabel(r'$u$')
	ax.view_init(40,-70)
	fig.tight_layout()
	fig.savefig('../figs/gaussian_convection.pdf')
	nfig += 1

# -- Show figures
plt.show()
