"""
magparsol/constants.py
-----------------------
Physical constants in SI units and derived convenience values.
All values follow CODATA recommendations where applicable.
"""

### Fundamental Constants ──────────────────────────────────────────────────────────
C               = 299_792_458.0       # speed of light in vacuum          [m/s]
Q_E             = 1.602_176_634e-19   # elementary (proton) charge        [C]
M_E             = 9.109_383_56e-31    # electron rest mass                [kg]
M_P             = 1.672_621_9e-27     # proton rest mass                  [kg]

### Earth Parameters ───────────────────────────────────────────────────────────────
R_EARTH         = 6_378_137.0         # earth equatorial radius           [m]
DIPOLE_MOMENT   = 7.965_626e15        # earth dipole moment coefficient   [T·m³]
DIPOLE_TILT_DEG = 30                  # magnetic axis tilt from rotation  [deg]

### Numerical Safety ───────────────────────────────────────────────────────────────
B_FLOOR         = 1e-30               # guard for division by zero        [T²]
