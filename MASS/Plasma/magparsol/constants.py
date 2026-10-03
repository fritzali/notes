"""
magparsol/constants.py
-----------------------
Physical constants in SI units and derived convenience values.
Fundamental constants follow CODATA 2018; the geomagnetic dipole follows
IGRF-14 (epoch 2025.0).
"""

import numpy as np

# ── Fundamental constants (CODATA 2018) ───────────────────────────────────────
C          = 299_792_458.0       # Speed of light in vacuum          [m/s]
Q_E        = 1.602_176_634e-19   # Elementary charge                 [C]
M_E        = 9.109_383_7015e-31  # Electron rest mass                [kg]
M_P        = 1.672_621_923_69e-27  # Proton rest mass                [kg]
EPS0       = 8.854_187_8128e-12  # Vacuum permittivity               [F/m]
MU0        = 1.256_637_062_12e-6 # Vacuum permeability               [N/A²]
K_B        = 1.380_649e-23       # Boltzmann constant                [J/K]

# ── Earth parameters ──────────────────────────────────────────────────────────
R_EARTH    = 6_378_137.0         # Earth equatorial radius (WGS-84)  [m]

# Centred-dipole part of IGRF-14 at epoch 2025.0 (Gauss coefficients, nT),
# reference radius a = 6371.2 km:
_IGRF_G10, _IGRF_G11, _IGRF_H11 = -29_350.0, -1_410.3, 4_545.5
_IGRF_A    = 6_371_200.0
_B0        = 1e-9 * np.sqrt(_IGRF_G10**2 + _IGRF_G11**2 + _IGRF_H11**2)  # 29.73 µT

DIPOLE_MOMENT = _B0 * _IGRF_A**3  # = mu0/4pi * m_Earth ≈ 7.69e15    [T·m³]
DIPOLE_TILT_DEG = float(np.degrees(np.arctan2(np.hypot(_IGRF_G11, _IGRF_H11),
                                              abs(_IGRF_G10))))  # ≈ 9.21  [deg]

# ── Numerical safety ──────────────────────────────────────────────────────────
B_FLOOR    = 1e-30               # |B|² threshold for Boris-C (avoid ÷0) [T²]
