# *-* coding: utf-8 *-*

"""
Created on mer 03 mar 2021 22:32:09 UTC

@author: vekemans

"""

import copy as cp
import numpy as np
from scipy import sparse

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.ticker import MaxNLocator

# -- Simulation Parameters
a = 0.01
dt_loc = 1e-6
t_max = 12e-3

# -- Problem Initialization

# -- Simulation
def simulate(c_init,N,dt):
    c = cp.copy(c_init); yield c.reshape((N,N))
    h = 1/N
    # Second Order Centered Difference Methods
    eyeN = np.ones(N)
    D = sparse.block_diag([sparse.spdiags([eyeN, eyeN, -eyeN*4, eyeN, eyeN], [-(N-1),-1,0,1,(N-1)], N,N)]*N) \
        + sparse.spdiags([np.ones(N*N)]*4, [-(N-1)*N,-N,N,(N-1)*N], N*N,N*N)

    while True:
        for _ in range(int(dt//dt_loc)):
            c = c + dt_loc/h**2 * D@(c**3 - c - (a/h)**2 * D@c)
        yield c.reshape((N,N))

if __name__ == "__main__":
    nfig = 1

    N = 128
    h = 1.0/N

    dt = 1e-5
    Nplot = int(t_max/dt)

    # -- Problem Initialisation
    c_init = np.random.random(N*N)*2-1
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

    cahn_hill = plt.imshow(next(data),cmap='binary',vmin=-1,vmax=+1,animated=True,interpolation='bessel')
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
    # anim = animation.FuncAnimation(fig, animate, interval=10, blit=False) # blit = False → update axis' title
    anim = animation.FuncAnimation(fig, animate, frames=Nplot+1, interval=10, blit=False)
    writervideo = animation.FFMpegWriter(fps=60) 
    anim.save('../figs/cahn_hilliard_spat_binary.mp4', writer=writervideo)

    plt.show()
