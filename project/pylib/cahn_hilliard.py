# *-* coding: utf-8 *-*

"""
Created on mer 03 mar 2021 22:32:09 UTC

@author: vekemans

"""

import numpy as np
from scipy import sparse

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.ticker import MaxNLocator
nfig = 1

# -- Simulation Parameters
a = 0.01
N = 128

h = 1.0/N
dt = 1e-6
Ntime = int(12e3)

C = np.random.random(N*N)*2-1

# -- Problem Initialization
# Second Order Centered Difference Methods
eyeN = np.ones(N)
D = sparse.block_diag([sparse.spdiags([eyeN, eyeN, -eyeN*4, eyeN, eyeN], [-(N-1),-1,0,1,(N-1)], N,N)]*N) \
    + sparse.spdiags([np.ones(N*N)]*4, [-(N-1)*N,-N,N,(N-1)*N], N*N,N*N)

# -- Simulation
# Initialize 'data' vector for final results
data = {}
data[0] = C.reshape((N,N))
for t_idx in range(Ntime):
    print('time : %d / %d' %(t_idx+1,Ntime), end='\r')
    F = C**3 - C - (a/h)**2 * D@C
    C = C + dt/h**2 * D@F
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
fig.colorbar(cahn_hill[0], ticks=_ticks, format='%.1f')

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
