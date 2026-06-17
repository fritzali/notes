"""
plasma_box/constants.py
-----------------------
Physical constants in SI units and derived convenience values.
All values follow CODATA 2018 recommendations where applicable.
"""

# ── Fundamental constants ──────────────────────────────────────────────────────
C          = 299_792_458.0       # Speed of light in vacuum          [m/s]
Q_E        = 1.602_176_634e-19   # Elementary (proton) charge        [C]
M_E        = 9.109_383_56e-31    # Electron rest mass                [kg]
M_P        = 1.672_621_9e-27     # Proton rest mass                  [kg]

# ── Earth parameters ──────────────────────────────────────────────────────────
R_EARTH    = 6_378_137.0         # Earth equatorial radius           [m]
DIPOLE_MOMENT = 7.965_626e15     # Earth dipole moment coefficient   [T·m³]
                                  # (= mu0/4pi * m_Earth, signed)
DIPOLE_TILT_DEG = 11.7           # Magnetic axis tilt from rotation axis [deg]

# ── Numerical safety ──────────────────────────────────────────────────────────
B_FLOOR    = 1e-30               # |B|² threshold for Boris-C (avoid ÷0) [T²]
