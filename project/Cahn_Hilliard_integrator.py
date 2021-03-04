import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import animation
import timeit

myColorMap = matplotlib.cm.jet

def wrapper(func, *args, **kwargs):
    def wrapped():
        return func(*args, **kwargs)
    return wrapped

def Hilliard2(C0,dt,dx,t,a):
  #plus rapide de +/- 1.5 s
  
  m = len(C0)    
  T = np.zeros((len(t),m,m)) # Anim
  F = np.zeros((m+2,m+2))
  C = np.zeros((m+2,m+2))
  C[1:-1,1:-1] = C0
  C[0,1:-1] = C0[-1,:]
  C[1:-1,0] = C0[:,-1]
  C[-1,1:-1] = C0[0,:]
  C[1:-1,-1] = C0[:,0]

   
  for it in range(len(t)):
      
    #T[it,:,:] = C  # Anim
    
    F[1:-1,1:-1] = C[1:-1,1:-1]**3 - C[1:-1,1:-1] - (a**2/dx**2) *(C[2:,1:-1] + C[0:-2,1:-1]   
                            + C[1:-1,2:] + C[1:-1,0:-2] - 4*C[1:-1,1:-1])
    
    C[1:-1,1:-1] += (dt/dx**2) * (F[2:,1:-1] + F[0:-2,1:-1]   
                            + F[1:-1,2:] + F[1:-1,0:-2] - 4*F[1:-1,1:-1])
    
    Cn = C[1:-1,1:-1]
    Fn = F[1:-1,1:-1]
    
    C[0,1:-1] = Cn[-1,:]
    C[1:-1,0] = Cn[:,-1]
    C[-1,1:-1] = Cn[0,:]
    C[1:-1,-1] = Cn[:,0]
    
    F[0,1:-1] = Fn[-1,:]
    F[1:-1,0] = Fn[:,-1]
    F[-1,1:-1] = Fn[0,:]
    F[1:-1,-1] = Fn[:,0]
    
  return Cn#,T

def Hilliard(C0,dt,dx,t,a):
  
  m = len(C0)    
  F = np.zeros((m,m))
  T = np.zeros((len(t),m,m)) # Anim
  C = C0
  
  for it in range(len(t)):
      
    #T[it,:,:] = C  # Anim
    
    # Domaine Periodique
    
    F[0,0] = C[0,0]**3 - C[0,0] - (a**2/dx**2) *(C[-1,0] + C[1,0]   
                            + C[0,-1] + C[0,1] - 4*C[0,0])
        
    C[0,0] += (dt/dx**2) * (F[-1,0] + F[1,0]   
                            + F[0,-1] + F[0,1] - 4*F[0,0])
    ###########################################################################
    F[0,-1] = C[0,-1]**3 - C[0,-1] - (a**2/dx**2) *(C[-1,-1] + C[1,-1]   
                            + C[0,-2] + C[0,0] - 4*C[0,-1])
        
    C[0,-1] += (dt/dx**2) * (F[-1,-1] + F[1,-1]   
                            + F[0,-2] + F[0,0] - 4*F[0,-1])
    ###########################################################################
    F[-1,0] = C[-1,0]**3 - C[-1,0] - (a**2/dx**2) *(C[-2,0] + C[0,0]   
                            + C[-1,-1] + C[-1,1] - 4*C[-1,0])
        
    C[-1,0] += (dt/dx**2) * (F[-2,0] + F[0,0]   
                            + F[-1,-1] + F[-1,1] - 4*F[-1,0])
    ###########################################################################
    F[-1,-1] = C[-1,-1]**3 - C[-1,-1] - (a**2/dx**2) *(C[-2,-1] + C[0,-1]   
                            + C[-1,-2] + C[-1,0] - 4*C[-1,-1])   
    
    C[-1,-1] += (dt/dx**2) * (F[-2,-1] + F[0,-1]   
                            + F[-1,-2] + F[-1,0] - 4*F[-1,-1])
    ###########################################################################
    F[0,1:-1] = C[0,1:-1]**3 - C[0,1:-1] - (a**2/dx**2) *(C[-1,1:-1] + C[1,1:-1]   
                            + C[0,0:-2] + C[0,2:m+1] - 4*C[0,1:-1])
        
    C[0,1:-1] += (dt/dx**2) * (F[-1,1:-1] + F[1,1:-1]   
                            + F[0,0:-2] + F[0,2:m+1] - 4*F[0,1:-1])
    ###########################################################################
    F[1:-1,0] = C[1:-1,0]**3 - C[1:-1,0] - (a**2/dx**2) *(C[0:-2,0] + C[2:m+1,0]   
                            + C[1:-1,-1] + C[1:-1,1] - 4*C[1:-1,0])
    
    
    C[1:-1,0] += (dt/dx**2) * (F[0:-2,0] + F[2:m+1,0]   
                            + F[1:-1,-1] + F[1:-1,1] - 4*F[1:-1,0])
    ###########################################################################
    F[-1,1:-1] = C[-1,1:-1]**3 - C[-1,1:-1] - (a**2/dx**2) *(C[-2,1:-1] + C[0,1:-1]   
                            + C[-1,0:-2] + C[-1,2:m+1] - 4*C[-1,1:-1])
        
    C[-1,1:-1] += (dt/dx**2) * (F[-2,1:-1] + F[0,1:-1]   
                            + F[-1,0:-2] + F[-1,2:m+1] - 4*F[-1,1:-1])
    ###########################################################################
    F[1:-1,-1] = C[1:-1,-1]**3 - C[1:-1,-1] - (a**2/dx**2) *(C[0:-2,-1] + C[2:m,-1]   
                            + C[1:-1,-2] + C[1:-1,0] - 4*C[1:-1,-1])
       
    C[1:-1,-1] += (dt/dx**2) * (F[0:-2,-1] + F[2:m,-1]   
                            + F[1:-1,-2] + F[1:-1,0] - 4*F[1:-1,-1])

    ###########################################################################   
    # Le reste
    F[1:-1,1:-1] = C[1:-1,1:-1]**3 - C[1:-1,1:-1] - (a**2/dx**2) *(C[2:,1:-1] + C[0:-2,1:-1]   
                            + C[1:-1,2:] + C[1:-1,0:-2] - 4*C[1:-1,1:-1])
    
    C[1:-1,1:-1] += (dt/dx**2) * (F[2:,1:-1] + F[0:-2,1:-1]   
                            + F[1:-1,2:] + F[1:-1,0:-2] - 4*F[1:-1,1:-1])
        
  return C#,T


ex = 7
n = 128
dx = 1/n 
x = np.linspace(0,1,n+1) 
a = 0.01

dt = 1e-6
t = np.arange(0,12000*dt+dt,dt)

C0 = np.random.rand(n+1,n+1) - np.random.rand(n+1,n+1)

H2 = Hilliard2(C0,dt,dx,t,a)
#∟H = Hilliard(C0,dt,dx,t,a)

#Plot
plt.figure()
plt.imshow(H2,cmap=myColorMap)
plt.colorbar()
plt.show()
#plt.figure()
#plt.imshow(H,cmap=myColorMap)
#plt.colorbar()
#plt.show()


#Animation

#fig = plt.figure()
#
#im = plt.imshow(C0, animated=True,cmap=myColorMap)
#
#
#def updatefig(i):
#
#    im.set_array(T[i,:,:])
#    return im,
#
#ani = animation.FuncAnimation(fig, updatefig, interval=10, blit=True)
##Writer = animation.writers['ffmpeg']
##writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
##ani.save('im.mp4', writer=writer)
#plt.show()


# Timing
#wrapped = wrapper(Hilliard,C0,dt,dx,t,a)
#time = timeit.timeit(wrapped,number = 20)
#print('Hilliard : ' + str(time))
#
#wrapped2 = wrapper(Hilliard2,C0,dt,dx,t,a)
#time2 = timeit.timeit(wrapped2,number = 20)
#
#print('Hilliard2 : ' + str(time2))