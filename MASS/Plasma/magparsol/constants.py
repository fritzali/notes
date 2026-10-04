"""
Physical constants in SI units.

Fundamental constants are CODATA 2018 values. The Earth's field is the
centred dipole of IGRF-14 at epoch 2025.0.
"""

import numpy as np

# CODATA 2018
C = 299_792_458.0               # speed of light [m/s]
Q_E = 1.602_176_634e-19         # elementary charge [C]
M_E = 9.109_383_7015e-31        # electron mass [kg]
M_P = 1.672_621_923_69e-27      # proton mass [kg]
EPS0 = 8.854_187_8128e-12       # vacuum permittivity [F/m]
MU0 = 1.256_637_062_12e-6       # vacuum permeability [N/A^2]
K_B = 1.380_649e-23             # Boltzmann constant [J/K]

R_EARTH = 6_378_137.0           # equatorial radius, WGS-84 [m]

# Degree-1 Gauss coefficients of IGRF-14 at 2025.0 [nT], reference radius
# a = 6371.2 km. Their magnitude is the equatorial surface field of the
# centred dipole, B0 ≈ 29.7 µT.
_IGRF_G10, _IGRF_G11, _IGRF_H11 = -29_350.0, -1_410.3, 4_545.5
_IGRF_A = 6_371_200.0
_B0 = 1e-9 * np.sqrt(_IGRF_G10**2 + _IGRF_G11**2 + _IGRF_H11**2)

# mu0/(4 pi) times the Earth's dipole moment, about 7.69e15 T m^3
DIPOLE_MOMENT = _B0 * _IGRF_A**3
# Angle between dipole and rotation axis, about 9.21 degrees
DIPOLE_TILT_DEG = float(np.degrees(np.arctan2(np.hypot(_IGRF_G11, _IGRF_H11),
                                              abs(_IGRF_G10))))

# Below this |B|^2 [T^2] Boris C treats the field as zero (avoids 0/0)
B_FLOOR = 1e-30
