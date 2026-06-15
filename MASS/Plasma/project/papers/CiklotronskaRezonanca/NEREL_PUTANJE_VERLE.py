# For the given analytical configuration of the relevant fields, it is necessary to display the trajectories of non-relativistic particles.
# Since the validity of the orbital method is assumed, the plasma is considered to be of sufficiently low concentration
# and collision processes are sufficiently rare events so that the motion of the particle collective can be judged
# based on the motion of a single particle. Coulomb interactions between particles are thus neglected.
# This program is inspired by the conceptual code from V. Zekovic & B. Arbutina, Nucl. Part. Phys. Proc. 297-299, 53-57 (2018)
# but with a modified method of numerical integration of the non-relativistic equation of motion
##################################
# Loading necessary libraries
##################################
import numpy as np
import matplotlib.pyplot as plt
import PVREDNOSTI, IZLAZ, DEFPOLJA # In the same directory as the main program
from mpl_toolkits.mplot3d import Axes3D
#############################
# Input data for the simulation
#############################
N = 1 # Number of particles in the considered plasma - e.g. one particle
T_sim = 3000.0 # SIMULATION DURATION [s]
dt = 0.1 # TIME STEP (desirable to be small compared to one gyro-period 2πm/qB)
T_smp = 5.0 # SAMPLING INTERVAL for graphical display
##########################
# Defining variables
##########################
l = 1
q = np.zeros(N)
m = np.zeros(N)
x = 0.0
y = 0.0
z = 0.0
vx = np.zeros(N)
vy = np.zeros(N)
vz = np.zeros(N)
EPx = 0.0
EPy = 0.0
EPz = 0.0
BPx = 0.0
BPy = 0.0
BPz = 0.0
pomx = np.zeros(N)
pomy = np.zeros(N)
pomz = np.zeros(N)
#################
# Initializations
#################
t = 0.0
q, m, x, y, z, vx, vy, vz = PVREDNOSTI.puslovi(N, np)
BPx, BPy, BPz, EPx, EPy, EPz = DEFPOLJA.const_Bpluswave(x, y, z, t, np)
fig, sp1, sp2 = IZLAZ.init_plots(plt)
###############
# Trajectory calculation
###############
# The magnetic Verlet algorithm is used in this example for the non-relativistic equation of motion
# See in detail
# https://aapt.scitation.org/doi/10.1119/10.0001876
# https://arxiv.org/abs/2008.11810
# https://www.compadre.org/PICUP/resources/Numerical-Integration/
#kinetickainit = m*0.5*(vx**2.0 + vy**2.0 + vz**2.0)
xx = np.zeros(30001)
yy = np.zeros(30001)
zz = np.zeros(30001)
ii = 0
xx[ii] = x
yy[ii] = y
zz[ii] = z
while t < (T_sim - dt):
    t = t + dt
    brzina = np.sqrt(vx**2.0 + vy**2.0 + vz**2.0)
    print(brzina)
    pomx = vx + ((q*dt*0.5/m)*(EPx + (vy*BPz) - (BPy*vz)))
    pomy = vy + ((q*dt*0.5/m)*(EPy + (vz*BPx) - (BPz*vx)))
    pomz = vz + ((q*dt*0.5/m)*(EPz + (vx*BPy) - (BPx*vy)))
    x = x + (pomx*dt)
    y = y + (pomy*dt)
    z = z + (pomz*dt)
    BPx, BPy, BPz, EPx, EPy, EPz = DEFPOLJA.const_Bpluswave(x, y, z, t, np)
    pomx = pomx + ((q*dt*0.5/m)*EPx)
    pomy = pomy + ((q*dt*0.5/m)*EPy)
    pomz = pomz + ((q*dt*0.5/m)*EPz)
    vx = ((1.0 + (((q*dt*0.5/m)**2.0)*(BPx**2.0 + BPy**2.0 + BPz**2.0)))**(-1.0))*(pomx + ((q*dt*0.5/m)*(pomy*BPz - pomz*BPy)) + (((q*dt*0.5/m)**2.0)*BPx*(pomx*BPx + pomy*BPy + pomz*BPz)))
    vy = ((1.0 + (((q*dt*0.5/m)**2.0)*(BPx**2.0 + BPy**2.0 + BPz**2.0)))**(-1.0))*(pomy + ((q*dt*0.5/m)*(pomz*BPx - pomx*BPz)) + (((q*dt*0.5/m)**2.0)*BPy*(pomx*BPx + pomy*BPy + pomz*BPz)))
    vz = ((1.0 + (((q*dt*0.5/m)**2.0)*(BPx**2.0 + BPy**2.0 + BPz**2.0)))**(-1.0))*(pomz + ((q*dt*0.5/m)*(pomx*BPy - pomy*BPx)) + (((q*dt*0.5/m)**2.0)*BPz*(pomx*BPx + pomy*BPy + pomz*BPz)))
    #print(t)
    ii = ii + 1
    xx[ii] = x
    yy[ii] = y
    zz[ii] = z
# This could also be written more concisely using np.dot and np.cross for vectors, instead of this scalar way (it is left to the students to improve this code)
###################################
# Graphical representation of the results
###################################
    l += 1
    if l*dt >= T_smp:
        l = 0.0
        IZLAZ.write_plots(sp1, sp2, plt, x, y, z, t, T_sim) # Display in real time with a larger step than calculated
#kinetickafin = m*0.5*(vx**2.0 + vy**2.0 + vz**2.0)
d = 6378137.0
fig = plt.figure() #graphic display
ax = fig.add_subplot(111, projection='3d')
ax.set_aspect('equal')
ax.grid(False)
plt.xlabel("$x[R_T]$")
plt.ylabel("$y[R_T]$")
ax.set_zlabel("$z[R_T]$")
plt.axis('on')
ax.set_xlim3d(-70, 70) #modify as needed
ax.set_ylim3d(-70, 70)
ax.set_zlim3d(-70, 70)
plt.plot(xx/d, yy/d, zz/d, color = 'green')
plt.show()
