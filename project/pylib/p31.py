# *-* coding: utf-8 *-*

"""
Created on  lun 03 mai 2021 14:50:29 CEST 

@author : vekemans
p31.py - gamma function via complex integral, trapezoid rule

"""

from numpy import *
res = finfo('float32').resolution

from matplotlib import pyplot as plt

N = 70
theta = (2*pi/N)*(arange(N)+0.5) - pi
c = -11 # center of circle of integration
r = 16 # radius of circle ofintegration

x = arange(-3.5,4+res,0.1)
y = arange(-2.5,2.5+res,0.1)
xx,yy = meshgrid(x,y)

zz = xx + yy*1j
gaminv = zeros(zz.shape)
for i in range(N):
	t = c + r*exp(1j*theta[i])
	gaminv = gaminv + exp(t)*t**(-zz)*(t-c)
gaminv = gaminv/N
gam = 1./gaminv

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
# ax.plot_wireframe(xx,yy,abs(gam), color='k')
# ax.plot_surface(xx,yy,abs(gam), cmap='jet', edgecolor='black', shade=False)
ax.plot_surface(xx,yy,abs(gam), color='white', edgecolor='black', shade=False)
ax.set_title(r'$|\Gamma(z)|$')

ax.set_xlim((-3.5,4)); ax.set_xlabel(r'$\mathcal{R}(z)$')
ax.set_ylim((-2.5,2.5)); ax.set_ylabel(r'$\mathcal{I}(z)$')
ax.set_zlim((0,6)); ax.view_init(20,-120)

plt.tight_layout()

plt.show()
