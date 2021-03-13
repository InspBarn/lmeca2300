# *-* coding: utf-8 *-*

"""
Created on  ven 12 mar 2021 15:10:09 CET 

@author : vekemans

"""

import math as mt
import numpy as np

from numpy.fft import fft2,ifft2,fftshift,fftfreq

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.ticker import MaxNLocator
nfig = 1

# -- Simulation Parameters
a = 0.01
N = 128

h = 1.0/N
dt = 1e-6/4.
Ntime = int(1e3)

C = np.random.random(N*N)*2-1
C = C.reshape((N,N))

# -- Problem Initialisation
def laplace(f):
    F = fftshift(fft2(f))
    # u = fftfreq(N)
    u = np.repeat(np.arange(-N/2,N/2), N).reshape((N,N))
    v = u.T

    # F_ddot_x = (2j*np.pi*u)**2 * F
    # F_ddot_y = (2j*np.pi*v)**2 * F
    F_ddot_x = (1j*u)**2 * F
    F_ddot_y = (1j*v)**2 * F

    f_ddot_x = ifft2(fftshift(F_ddot_x))
    f_ddot_y = ifft2(fftshift(F_ddot_y))

    return (f_ddot_x+f_ddot_y)/h**2

f = lambda x: laplace(x**3 - x - a**2*laplace(x))

def rk4(C):
    k1 = f(C)
    k2 = f(C+dt*k1/2)
    k3 = f(C+dt*k2/2)
    k4 = f(C+dt*k3)

    return C + (k1 + 2*k2 + 2*k3 + k4)*dt/6

# -- Simulation
# Initialize 'data' vector for final results
data = {}
data[0] = C.reshape((N,N))
for t_idx in range(Ntime):
    print('time : %d / %d' %(t_idx+1,Ntime), end='\r')
    # C = rk4(C)
    C = C + f(C) * dt
    data[t_idx+1] = C.reshape((N,N))

# -- Animation Parameters
# Colorbar
nbins = 200
_bounds = max(-data[t_idx].min(),data[t_idx].max()) # Maximal concentration along the simulation
_levels = MaxNLocator(nbins=nbins).tick_values(-_bounds,_bounds) # Levels of the colorbar
_ticks = np.append(np.flip(np.arange(0,-_bounds,-0.2)),np.arange(0,_bounds,0.2)) # ticks along the colorbar
# Grid coordinates
x = np.arange(0,N)*h + h/2
X,Y = np.meshgrid(x,x) # /!\ symmetrical

# -- Animated Initialization
fig = plt.figure(nfig)
ax = fig.add_subplot()

cahn_hill = [ plt.contourf(X,Y,data[0],_levels, cmap='jet') ]
# fig.colorbar(cahn_hill[0], ticks=_ticks, format='%.1f')

ax.set_title(r'Time : $t = %.5f$ [s]' %0.0)
ax.set_xlim(0,1); ax.set_xticks([])
ax.set_ylim(0,1); ax.set_yticks([])
fig.tight_layout()
nfig += 1

# -- Animation function
Nplots = Ntime//100
def animate(t):
    t_idx = t * (Ntime//Nplots)

    global cahn_hill
    for c in cahn_hill[0].collections:
        c.remove() # clean contourf

    ax.set_title(r'Time : $t = %.5f$ [s]' %(t_idx*dt))
    cahn_hill[0] = plt.contourf(X,Y,data[t_idx],_levels, cmap='jet')

    return cahn_hill[0].collections

if Ntime%Nplots==0:
    anim = animation.FuncAnimation(fig, animate, Nplots+1, interval=10, blit=False) # blit = False → update axis' title
else:
    anim = animation.FuncAnimation(fig, animate, Nplots, interval=10, blit=False)

plt.show()
