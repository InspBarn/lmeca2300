# *-* coding: utf-8 *-*

"""
Created on  ven 12 mar 2021 15:10:09 CET 

@author : vekemans

"""

import math as mt
import numpy as np

from numpy.fft import fft2,ifft2,fftfreq
from numpy.fft import rfft2,irfft2,rfftfreq

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
dt = 1e-6/4.
time_max = 6e-3
Ntime = int(time_max/dt)

# -- Problem Initialisation
C = np.random.random(N*N)*2-1
C = C.reshape((N,N))

u = rfftfreq(N, d=1/N).reshape((1,N//2+1))
v = fftfreq(N, d=1/N).reshape((1,N)).T
laplace = -(2*np.pi)**2 * (u**2 + v**2)

C_hat = rfft2(C)

f = lambda x: laplace * (rfft2(irfft2(x)**3) - x - a**2 * laplace*x)

def rk4(x):
    k1 = f(x)
    k2 = f(x+dt*k1/2)
    k3 = f(x+dt*k2/2)
    k4 = f(x+dt*k3)

    return x + (k1 + 2*k2 + 2*k3 + k4)*dt/6

# -- Simulation
# Initialize 'data' vector for final results
data = {}
data[0] = C
start,t_anl = time(),10
for t_idx in range(Ntime):
    print('time : %d / %d' %(t_idx+1,Ntime), end='\r')

    C_hat = rk4(C_hat)
    data[t_idx+1] = irfft2(C_hat)

#     if t_idx%t_anl==0:
#         end = time()
#         print('time for %d it : %.3f' %(t_anl,end-start), end=' ')
#         start = end
#     print('', end='\r')
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

cahn_hill = [ plt.contourf(X,Y,data[0],_levels, cmap='jet') ]
fig.colorbar(cahn_hill[0], ticks=_ticks, format='%.1f')

ax.set_title(r'Time : $t = %.5f$ [s]' %0.0)
ax.set_xlim(0,1); ax.set_xticks([])
ax.set_ylim(0,1); ax.set_yticks([])
fig.tight_layout()
nfig += 1

# -- Animation function
Nplots = Ntime//400
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
