# *-* coding: utf-8 *-*

"""
Created on  ven 12 mar 2021 15:10:09 CET 

@author : vekemans

"""

import math as mt
import numpy as np
import copy as cp

import scipy
from scipy.fft import rfft2,irfft2
from scipy.fft import rfftfreq,fftfreq
from scipy import linalg

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.ticker import MaxNLocator
nfig = 1

from time import time

# -- Simulation Parameters
a = 0.01
N = 128

h = 1.0/N
dt_red = 1
dt = 1e-6/dt_red
time_max = 12e-3
Ntime = int(time_max/dt/dt_red)

# -- Problem Initialisation
C = np.random.random(N*N)*2-1
C = C.reshape((N,N))

u = rfftfreq(N, d=1/N).reshape((1,N//2+1))
v = fftfreq(N, d=1/N).reshape((1,N)).T
laplace = -(2*np.pi)**2 * (u**2 + v**2)

C_hat = rfft2(C)

f = lambda c: laplace * (rfft2(irfft2(c)**3) - c - a**2 * laplace*c)
df = lambda c: laplace * (rfft2(irfft2(c)**2)/2/np.pi - 1 - a**2 * laplace)

def rk4(x):
    k1 = f(x)
    k2 = f(x+dt*k1/2)
    k3 = f(x+dt*k2/2)
    k4 = f(x+dt*k3)

    return x + (k1 + 2*k2 + 2*k3 + k4)*dt/6

def implicit_euler(xn):
    x = cp.copy(xn)
    it = 0
    error = 1
    while (error>1e-6) and (it<100):
        dx = (xn+dt*f(x)-x) / (1-dt*df(x))
        x += dx
        error = np.fabs(irfft2(dx)).max(axis=(0,1))
        it += 1
    return x


# -- Simulation
# Initialize 'data' vector for final results
data = {}
data[0] = C
start,t_anl = time(),10
for t_idx in range(Ntime):
    print('time : %d / %d' %(t_idx+1,Ntime), end='\r')

    for _ in range(dt_red):
        # C_hat = rk4(C_hat)
        C_hat = implicit_euler(C_hat)

    data[t_idx+1] = irfft2(C_hat)
print()

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

cahn_hill = plt.imshow(data[0],cmap='jet',vmin=-_bounds,vmax=+_bounds,animated=True,interpolation='bessel')
fig.colorbar(cahn_hill, ticks=_ticks, format='%.1f')

ax.set_title(r'Time : $t = %.5f$ [s]' %0.0)
ax.set_xticks([]); ax.set_yticks([])
nfig += 1

# -- Animation function
Nplots = Ntime//50
def animate(t):
    t_idx = t * (Ntime//Nplots)

    ax.set_title(r'Time : $t = %.5f$ [s]' %(t_idx*dt))
    cahn_hill.set_array(data[t_idx])

    return cahn_hill,

fig.tight_layout()
if Ntime%Nplots==0:
    anim = animation.FuncAnimation(fig, animate, Nplots+1, interval=10, blit=False) # blit = False → update axis' title
else:
    anim = animation.FuncAnimation(fig, animate, Nplots, interval=10, blit=False)

plt.show()
