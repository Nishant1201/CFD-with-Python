# Vorticity by vorticity-streamfunction method
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

nx = 65
ny = 65
max_step = 2000
visc = 0.001
dt = 0.001
t_current = 0.0
dx = 1.0/(nx-1)

max_iter = 100
beta = 1.5
max_err = 0.001

psi = np.zeros((nx+1,ny))
omega = np.zeros((nx+1,ny))
omega0 = np.zeros((nx+1,ny))
x = np.zeros((nx+1,ny))
y = np.zeros((nx+1,ny))

for i in range(nx+1):
    for j in range(ny):
        x[i,j] = dx*(i-1)
        y[i,j] = dx*(j-1)

# Initialize vorticity
omega[1:34,int((ny-1)/2-1)] = 1.0/dx
omega[32:nx-1,int((ny-1)/2+1)] = 1.0/dx

for tstep in range(max_step):
    for iter in range(max_iter):
        psi[0:-1,0] = 0.0
        psi[0:-1,-1] = 0.0
        psi[-1,:] = psi[1,:]
        omega0 = psi.copy()
        for i in range(1,nx):
            for j in range(1,ny-1):
                # Solve for psi by SOR method
                psi[i,j] = 0.25*beta*(psi[i+1,j] + psi[i-1,j] + psi[i,j+1] + psi[i,j-1] + dx*dx*omega[i,j]) + (1.0-beta)*psi[i,j]
        psi[0,:] = psi[nx-1,:]

        # Compute error
        err = 0.0
        for i in range(nx):
            for j in range(ny):
                err = err + np.abs(omega0[i,j] - psi[i,j])
        # Stop if converged
        if err <= max_err:
            break
        
    omega[-1,1:-1] = omega[1,1:-1]
    omega0 = omega.copy()

    for i in range(1,nx):
        for j in range(1,ny-1):
            omega[i,j] = omega0[i,j] + dt*( -0.25*( (psi[i,j+1] - psi[i,j-1])*(omega0[i+1,j]-omega0[i-1,j]) \
                            - (psi[i+1,j] - psi[i-1,j])*(omega0[i,j+1]-omega0[i,j-1]) )/(dx*dx) \
                            + visc*(omega0[i+1,j] + omega0[i-1,j] + omega0[i,j+1] + omega0[i,j-1] -4.0*omega0[i,j])/(dx*dx)  )

    omega[0,1:-1] = omega[nx-1,1:-1]
    t_current += dt
    print(t_current)
    #plt.subplot(1,2,1)
    #plt.cla()
    #plt.contour(x[0:-1,:],y[0:-1,:],omega[0:-1,:],40)
    #plt.axis('square')
    #plt.pause(0.01)
    
    #plt.subplot(1,2,2)
    #plt.contour(x[0:-1,:],y[0:-1,:],psi[0:-1,:],40)
    #plt.axis('square')
    #plt.pause(0.01)
    

print(err)
plt.subplot(1,2,1)
plt.contour(x[0:-1,:],y[0:-1,:],omega[0:-1,:],40)
plt.title(r'Vorticity, $\Omega$')
plt.axis('square')
plt.subplot(1,2,2)
plt.contour(x[0:-1,:],y[0:-1,:],psi[0:-1,:],40)
plt.axis('square')
plt.title(r'Streamfunction, $\Psi$')
plt.show()

#plt.contour(x,y,omega,40)
#plt.axis('square')
#plt.show()
#plt.contour(x,y,psi,20)
#plt.axis('square')
#plt.show()
 