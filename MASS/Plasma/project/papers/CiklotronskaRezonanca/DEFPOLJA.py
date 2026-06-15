#Example of cyclotron resonance
def const_Bpluswave(x, y, z, t, np):
# A homogeneous and stationary B field with Bz component (along the z-axis) and a wave described by the sinusoidal component of the electric field E are given
    BPx = 0.0*1e-10
    BPy = 0.0*1e-10
    BPz = 3.0*1e-10 #Given in micro Gauss [uG] times conversion - 3uG in this example
#For cyclotron resonance it is necessary that omega_c = omega_wave, omega_c = qB/m
    omegac = 1.6021766210e-19 * 3.0*1e-10 / 1.6726219e-27
    EPx = 0.0
    EPy = 0.0 + (5e-5 * np.sin(omegac * t))
    EPz = 0.0
    return BPx, BPy, BPz, EPx, EPy, EPz
