# *-* coding: utf-8 *-*

"""
Created on  ven 12 mar 2021 11:38:23 CET 

@author : vekemans

"""

import math as mt
import numpy as np

from numpy import pi as π
from numpy.fft import fft,fftshift,ifft
from scipy.interpolate import interp2d

import matplotlib
import matplotlib.pyplot as plt
nfig = 1

# Grid and initial data
N = 24
x = np.cos(π*np.linspace(0,1, N+1))

dt = 6/N**2
xx,yy = np.meshgrid(x,x)
plotgap = int(1/(3*dt))
dt = 1/(plotgap*3)
vv = np.exp(-40*((xx-0.4)**2 + yy**2))
vvold = vv

fig = plt.figure(nfig, figsize=(12.8,9.2))

# Time-stepping by leapfrog formula
ax_pos = 220
for n in range(3*plotgap+1):
    t = n*dt

    if n%plotgap==0:
        # plots at multiples of t=1/3
        ax_pos += 1
        ax = fig.add_subplot(ax_pos, projection='3d')

        x_ax = np.arange(-1,1,step=1/16)
        xxx,yyy = np.meshgrid(x_ax,x_ax)
        vvv = interp2d(xx,yy,vv, kind='cubic')
        ax.plot_wireframe(xxx,yyy,vvv(x_ax,x_ax), color='k')

        ax.set_xlim((-1,1))
        ax.set_ylim((-1,1))
        ax.set_zlim((-0.15,1))
        ax.set_title('t = %.3f' %t)

    uxx = np.zeros((N+1,N+1))
    uyy = np.zeros((N+1,N+1))
    ii = np.arange(1,N)
    for i in range(1,N):
        # 2nd derivs wrt x in each row
        v = vv[i,:]
        V = np.append(v, np.flip(v[ii]))
        U = np.real(fft(V))
        W1 = np.real(ifft(1j*U * np.concatenate((np.arange(N), [0], np.arange(1-N,0)))))   # diff wrt theta
        W2 = np.real(ifft(-U * np.concatenate((np.arange(N+1), np.arange(1-N,0)))**2))     # diff^2 wrt theta
        uxx[i,ii] = W2[ii]/(1-x[ii]**2) - x[ii]*W1[ii]/(1-x[ii]**2)**(3/2)

    for j in range(2,N):
        # 2nd derivs wrt y in each column
        v = vv[:,j]
        V = np.append(v, np.flip(v[ii]))
        U = np.real(fft(V))
        W1 = np.real(ifft(1j*U * np.concatenate((np.arange(N), [0], np.arange(1-N,0)))))    # diff wrt theta
        W2 = np.real(ifft(-U * np.concatenate((np.arange(N+1), np.arange(1-N,0)))**2))      # diff^2 wrt theta
        uyy[ii,j] = W2[ii]/(1-x[ii]**2) - x[ii]*W1[ii]/(1-x[ii]**2)**(3/2)

    vvnew = 2*vv - vvold + (uxx+uyy)*dt**2
    vvold = vv
    vv = vvnew

plt.show()
