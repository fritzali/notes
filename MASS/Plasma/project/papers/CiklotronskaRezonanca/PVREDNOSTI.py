###########################
# Setting initial conditions
###########################
def puslovi(N, np): 
    import random
########################################
# Input data (physical characteristics)
########################################
    RZ = 6378137.0
    m_e = 9.10938356e-31 # Electron mass [kg]
    m_p = 1.6726219e-27 # Proton mass [kg]
    qelem = 1.6021766210e-19 # Elementary charge [C]
    c = 299792458.0 # Speed of light in vacuum [m/s]
    vini = 0.01*c # Initial velocity limit - initial velocities are in the interval (-vini, vini)
    q = np.zeros(N)
    m = np.zeros(N)
    x = 0.0
    y = 0.0
    z = 0.0
    vx = np.zeros(N)
    vy = np.zeros(N)
    vz = np.zeros(N)
    for i in range(0,N):
        q[i] = qelem # Proton motion is observed in this example
        m[i] = m_p
        x = 0.0
        y = 0.0
        z = 0.0
        vx[i] = 0.0
        vy[i] = (random.random() - 0.5)*vini
        vz[i] = (random.random() - 0.5)*vini
    return q, m, x, y, z, vx, vy, vz
