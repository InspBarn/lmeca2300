# *-* coding: utf-8 *-*

"""
Created on  ven 12 mar 2021 11:38:23 CET 

@author : vekemans
p20.py - 2nd-order wave eq. in 2D via FFT (compare p19.py)

"""

from numpy import *
from numpy.fft import fft,fftshift,ifft
from scipy.interpolate import interp2d

from matplotlib import pyplot as plt

# Grid and initial data
N = 24
x = cos(pi*linspace(0,1, N+1))

dt = 6/N**2
xx,yy = meshgrid(x,x)
plotgap = int((1./3.)/dt)
dt = (1./3.)/plotgap
vv = exp(-40.*((xx-.4)**2 + yy**2))
vvold = vv

fig = plt.figure(figsize=(6.4*2,4.8*2))

# Time-stepping by leapfrog formula
ax_pos = 220
for n in range(3*plotgap+1):
    t = n*dt

    if remainder(n+.5, plotgap) < 1:
        # plots at multiples of t=1/3
        ax_pos += 1
        ax = fig.add_subplot(ax_pos, projection='3d')

        x4ax = arange(-1,1+1e-5, step=1/16)
        xxx,yyy = meshgrid(x4ax,x4ax)
        vvv = interp2d(x,x, vv, kind='cubic')(x4ax,x4ax)
        # ax.plot_wireframe(xxx,yyy,vvv, color='k', linewidth=0.8)
        # ax.plot_surface(xxx,yyy,vvv, cmap='jet', edgecolor='black')
        # ax.plot_surface(xxx,yyy,vvv, color='white', edgecolor='black', shade=False)

        ax.set_xlim((-1,1)); ax.set_ylim((-1,1)); ax.set_zlim((-0.15,1))
        ax.set_xticks([-1,0,1]); ax.set_yticks([-1,0,1]); ax.set_zticks([0,0.5,1])
        ax.set_title('t = %.3f' %t); ax.view_init(20, -125)

    uxx = zeros((N+1,N+1))
    uyy = zeros((N+1,N+1))
    ii = arange(1,N)

    for i in range(1,N):
        # 2nd derivs wrt x in each row
        v = vv[i,:]; V = concatenate((v, flipud(v[ii]))); U = fft(V).real
        W1 = ifft(1j*concatenate((arange(N), [0], arange(1-N,0))) * U).real     # diff wrt theta
        W2 = ifft(-concatenate((arange(N+1), arange(1-N,0)))**2 * U).real       # diff^2 wrt theta
        uxx[i,ii] = W2[ii]/(1-x[ii]**2) - x[ii]*W1[ii]/(1-x[ii]**2)**(3/2)

    for j in range(1,N):
        # 2nd derivs wrt y in each column
        v = vv[:,j]; V = append(v, flipud(v[ii])); U = fft(V).real
        W1 = ifft(1j*concatenate((arange(N), [0], arange(1-N,0))) * U).real     # diff wrt theta
        W2 = ifft(-concatenate((arange(N+1), arange(1-N,0)))**2  *U).real       # diff^2 wrt theta
        uyy[ii,j] = W2[ii]/(1-x[ii]**2) - x[ii]*W1[ii]/(1-x[ii]**2)**(3/2)

    vvnew = 2*vv - vvold + (uxx+uyy)*dt**2
    vvold = vv
    vv = vvnew

fig.tight_layout()
plt.show()
