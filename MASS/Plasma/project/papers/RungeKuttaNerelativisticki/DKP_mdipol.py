#Motion of non-relativistic particles in the approximation of a time-independent dipole magnetic field of planet Earth
#The RK4 (Runge-Kutta order 4) algorithm is used for numerical integration
#This program is written following the program from the paper Garcia-Farieta & Hurtado (2019) given in the appendix
#https://rmf.smf.mx/ojs/rmf-e/article/view/524
#Numerical values of physical constants were used from the paper Ozturk (2012)
#https://aapt.scitation.org/doi/10.1119/1.3684537
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
c = 299792458.0 #[m/s]
RZ = 6378137.0 #[m]
m_p = 1.6726219e-27 #[kg]
m_e = 9.10938356e-31 #[kg]
qe = 1.6021766210e-19 #[C]
sinphi = np.sin(11.7*np.pi/180.0) #11.7 degrees is the angle that the magnetic axis makes with the Earth's rotation axis
cosphi = np.cos(11.7*np.pi/180.0)
def BP(R):
    pom = (R[0]**2.0 + R[1]**2.0 + R[2]**2.0)**2.5
    BPx = -7.965626e15*((3.0*R[0]*R[2]*cosphi) + (3.0*R[0]*R[1]*sinphi))/pom
    BPy = -7.965626e15*((3.0*R[1]*R[2]*cosphi) + (2.0*sinphi*(R[1]**2.0)) - (sinphi*(R[0]**2.0)) - (sinphi*(R[2]**2.0)))/pom
    BPz = -7.965626e15*((2.0*cosphi*(R[2]**2.0)) - (cosphi*(R[0]**2.0)) - (cosphi*(R[1]**2.0)) + (3.0*R[2]*R[1]*sinphi))/pom
    BP = np.array([BPx, BPy, BPz])
    return BP #approximation of the dipole magnetic field of planet Earth (see appendix with derivations)
dt = 0.001 #integration step - modify depending on the situation, i.e. initial conditions
ts = 5000.0 #simulation duration
Nkoraka = int(ts/dt) #number of steps
t = np.zeros(Nkoraka) #initializations
rvek = np.zeros((len(t), 3))
vvek = np.zeros((len(t), 3))
m = 1.0*m_p #choice of particle type, by default it is a proton
q = 1.0*qe
#different initial conditions
#dt = 0.001
#ts = 5000.0
rvek[0, :] = np.array([0.0, -7.85, -1.53])*RZ
vvek[0, :] = np.array([0.0, 1.8e6, 1.8e6])
#dt = 0.001
#ts = 50000.0
#rvek[0, :] = np.array([-4.0, -1.0, -6.0])*RZ
#vvek[0, :] = np.array([6e5, 6e5, 6e5])
#dt = 0.01
#ts = 500000.0
#rvek[0, :] = np.array([-30.0, -30.0, -30.0])*RZ
#vvek[0, :] = np.array([-1.2e4, -6e3, 6e4])
#dt = 0.01
#ts = 500000.0
#rvek[0, :] = np.array([-35.0, -30.0, -35.0])*RZ
#vvek[0, :] = np.array([1.2e4, 6e3, 0.0])
for i in range(1, Nkoraka): #RK4 algorithm
    rk1 = rvek[i-1, :]
    vk1 = vvek[i-1, :]
    ak1 = (q/m)*np.cross(vk1, BP(rk1))
    rk2 = rvek[i-1, :] + (0.5*vk1*dt)
    vk2 = vvek[i-1, :] + (0.5*ak1*dt)
    ak2 = (q/m)*np.cross(vk2, BP(rk2))
    rk3 = rvek[i-1, :] + (0.5*vk2*dt)
    vk3 = vvek[i-1, :] + (0.5*ak2*dt)
    ak3 = (q/m)*np.cross(vk3, BP(rk3))
    rk4 = rvek[i-1, :] + (vk3*dt)
    vk4 = vvek[i-1, :] + (ak3*dt)
    ak4 = (q/m)*np.cross(vk4, BP(rk4))
    rvek[i] = rvek[i-1, :] + (dt/6.0)*(vk1 + (2.0*vk2) + (2.0*vk3) + vk4)
    vvek[i] = vvek[i-1, :] + (dt/6.0)*(ak1 + (2.0*ak2) + (2.0*ak3) + ak4)
    t[i] = dt*i
    print(t[i], rvek[i, :]/RZ)
fig = plt.figure() #graphic display
ax = fig.add_subplot(111, projection='3d')
ax.set_aspect('auto')
ax.grid(False)
u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:50j] #drawing a rough sketch of planet Earth as a sphere
x = np.cos(u)*np.sin(v)
y = np.sin(u)*np.sin(v)
z = np.cos(v)
ax.plot_wireframe(x, y, z, color = "blue")
plt.xlabel("$x[R_T]$")
plt.ylabel("$y[R_T]$")
ax.set_zlabel("$z[R_T]$")
plt.axis('on')
ax.set_xlim3d(-75, 75) #modify as needed
ax.set_ylim3d(-75, 75)
ax.set_zlim3d(-75, 75)
plt.plot(rvek[:, 0]/RZ, rvek[:, 1]/RZ, rvek[:, 2]/RZ, color = 'green')
plt.show()
