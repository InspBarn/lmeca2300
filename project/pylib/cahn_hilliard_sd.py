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
dt_loc = 1e-6
t_max = 12e-3


def rk4(dt,c,func):
    dt_red = 4
    for _ in range(dt_red):
        k1 = func(c)
        k2 = func(c+(dt/dt_red)*k1/2)
        k3 = func(c+(dt/dt_red)*k2/2)
        k4 = func(c+(dt/dt_red)*k3)
        c = c + (k1 + 2*k2 + 2*k3 + k4)*(dt/dt_red)/6
    return c

def euler_implicit(dt,c,func,dfunc):
    c_next = cp.copy(c)
    it = 0
    error = 1
    while (error>1e-6) and (it<100):
        dx = (c+dt*func(c_next)-c_next) / (1-dt*dfunc(c_next))
        c_next += dx
        error = np.fabs(irfft2(dx)).max(axis=(0,1))
        # error = np.linalg.norm(dx)
        it += 1
    return c_next

def bdf_ab(t,dt,c_hat,c_hat_prev, k):
    c = irfft2(c_hat)
    # c[c > 1] = 1
    # c[c < -1] = -1
    f_hat = rfft2(c**3-c)

    if t==0:
        c_hat_next = (c_hat - dt*k**2*f_hat) / (1+dt*(a*k**2)**2)
    else:
        c_prev = irfft2(c_hat_prev); f_hat_prev = rfft2(c_prev**3-c_prev)
        c_hat_next = (4*c_hat-c_hat_prev - 2*dt*k**2*(2*f_hat - f_hat_prev)) / (3+2*dt*(a*k**2)**2)

    return c_hat_next, c_hat


# -- Simulation
def simulate(c_init,N,dt,temp_scheme='bdf/ab'):
    yield c_init
    c_hat = rfft2(c_init)
    c_hat_previous = cp.copy(c_hat)
    t = 0
    # Spectral frequencies // Spectral laplacian
    u = rfftfreq(N, d=1/N).reshape((1,N//2+1))
    v = fftfreq(N, d=1/N).reshape((1,N)).T
    k = (u**2 + v**2)**0.5 * (2*np.pi)
    # Cahn-Hilliard derivation function && lagrangian
    f = lambda c: -k**2*(-c+rfft2(irfft2(c)**3)) - c*(a*k**2)**2
    # df = lambda c: -k**2*(-1+rfft2(irfft2(c)**2)/(2*np.pi)) - (a*k**2)**2

    while True:
        for _ in range(int(dt//dt_loc)):
            if temp_scheme=='rk4':
                c_hat = rk4(dt_loc, c_hat, f)
            elif temp_scheme=='bdf/ab':
                # c_hat = euler_implicit(dt_loc, c_hat, f,df)
                c_hat, c_hat_previous = bdf_ab(t,dt_loc, c_hat, c_hat_previous, k)
            t += dt_loc
        yield irfft2(c_hat)

if __name__ == "__main__":
    nfig = 1

    N = 128
    h = 1.0/N

    dt = 1e-5
    Nplot = int(t_max/dt)

    # -- Problem Initialisation
    c_init = np.random.random(N*N)*2-1
    c_init = c_init.reshape((N,N))
    # -- Problem Computation
    data = iter(simulate(c_init,N,dt))

    # -- Animation Parameters
    # Colorbar
#     nbins = 200
#     _bounds = max(-data[t_idx].min(),data[t_idx].max()) # Maximal concentration along the simulation
#     _levels = MaxNLocator(nbins=nbins).tick_values(-_bounds,_bounds) # Levels of the colorbar
#     _ticks = np.append(np.flip(np.arange(0,-_bounds,-0.2)),np.arange(0,_bounds,0.2)) # ticks along the colorbar
    # Grid coordinates
    x = np.arange(0,N)*h + h/2
    X,Y = np.meshgrid(x,x) # /!\ symmetrical

    # -- Animated Initialization
    fig = plt.figure(nfig)
    ax = fig.add_subplot()

    cahn_hill = plt.imshow(c_init,cmap='jet',vmin=-1,vmax=+1,animated=True,interpolation='bessel')
    fig.colorbar(cahn_hill, format='%.1f')

    ax.set_title(r'Time : $t = %.5f$ [s]' %0.0)
    ax.set_xticks([]); ax.set_yticks([])
    nfig += 1

    # -- Animation function
    def animate(t):
        ax.set_title(r'Time : $t = %.2f$ [ms]' %(t*dt*1e3))
        # cahn_hill.set_array(data[t_idx])
        cahn_hill.set_array(next(data))

        return cahn_hill,

    fig.tight_layout()
    anim = animation.FuncAnimation(fig, animate, interval=10, blit=False) # blit = False → update axis' title
    plt.show()

    anim = animation.FuncAnimation(fig, animate, frames=Nplot+1, interval=10, blit=False)
    writervideo = animation.FFMpegWriter(fps=60) 
    anim.save('../figs/cahn_hilliard_bdf_ab.mp4', writer=writervideo)
    plt.show()
