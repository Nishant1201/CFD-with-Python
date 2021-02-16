# Lid driven cavity by vorticity-streamfunction method
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

nx = 17
ny = 17
max_step = 200
visc = 0.1
dt = 0.005
t_current = 0.0
dx = 1.0/(nx-1)

u_wall = 1.0

max_iter = 100
beta = 1.5
max_err = 0.001

psi = np.zeros((nx,ny))
omega = np.zeros((nx,ny))
omega0 = np.zeros((nx,ny))
x = np.zeros((nx,ny))
y = np.zeros((nx,ny))

for i in range(nx):
    for j in range(ny):
        x[i,j] = dx*(i-1)
        y[i,j] = dx*(j-1)

for tstep in range(max_step):
    for iter in range(max_iter):
        omega0 = psi.copy()
        for i in range(1,nx-1):
            for j in range(1,ny-1):
                # Solve for psi by SOR method
                psi[i,j] = 0.25*beta*(psi[i+1,j] + psi[i-1,j] + psi[i,j+1] + psi[i,j-1] + dx*dx*omega[i,j]) + (1.0-beta)*psi[i,j]
        
        err = 0.0
        for i in range(nx):
            for j in range(ny):
                err = err + np.abs(omega0[i,j] - psi[i,j])
        if err <= max_err:
            break

    omega[1:nx-1,0] = -2.0*psi[1:nx-1,1]/(dx*dx)                        # vorticity on bottom wall
    omega[1:nx-1,ny-1] = -2.0*psi[1:nx-1,ny-2]/(dx*dx) - 2.0*u_wall/dx   # vorticity on top moving wall
    omega[0,1:ny-1] = -2.0*psi[1,1:ny-1]/(dx*dx)                        # vorticity on left wall
    omega[nx-1,1:ny-1] = -2.0*psi[nx-2,1:ny-1]/(dx*dx)                  # vorticity on right wall

    omega0 = omega.copy()

    for i in range(1,nx-1):
        for j in range(1,ny-1):
            omega[i,j] = omega0[i,j] + dt*( -0.25*( (psi[i,j+1] - psi[i,j-1])*(omega0[i+1,j]-omega0[i-1,j]) \
                            - (psi[i+1,j] - psi[i-1,j])*(omega0[i,j+1]-omega0[i,j-1]) )/(dx*dx) \
                            + visc*(omega0[i+1,j] + omega0[i-1,j] + omega0[i,j+1] + omega0[i,j-1] -4.0*omega0[i,j])/(dx*dx)  )

    t_current += dt
    
    #plt.subplot(1,2,1)
    #plt.cla()
    #plt.contour(x,y,omega,40)
    #plt.axis('square')
    #plt.pause(0.01)
    
    #plt.subplot(1,2,2)
    #plt.contour(x,y,psi,20)
    #plt.axis('square')
    #plt.pause(0.01)
    


print(err)
plt.subplot(1,2,1)
plt.contour(x,y,omega,40)
plt.title(r'Vorticity, $\Omega$')
plt.axis('square')
plt.subplot(1,2,2)
plt.contour(x,y,psi,20)
plt.axis('square')
plt.title(r'Streamfunction, $\Psi$')
plt.show()

#plt.contour(x,y,omega,40)
#plt.axis('square')
#plt.show()
#plt.contour(x,y,psi,20)
#plt.axis('square')
#plt.show()
 