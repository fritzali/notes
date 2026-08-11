"""
plasma_box/fields.py
--------------------
Electromagnetic field models.

Every FieldModel subclass implements:

    evaluate(r, t) -> (B, E)

where r is an (N, 3) position array and t is a scalar time,
and the return values B, E are (N, 3) arrays in SI units (Tesla, V/m).

A CustomField adapter lets users supply a legacy scalar-signature function::

    def my_field(x, y, z, t):
        return Bx, By, Bz, Ex, Ey, Ez          # scalars or length-N arrays

or a modern vector-signature function::

    def my_field(r, t):                         # r shape (N,3)
        return B, E                              # each shape (N,3)

and wrap it::

    field = CustomField(my_func, vector_api=False)   # legacy
    field = CustomField(my_func, vector_api=True)    # modern
"""

import numpy as np
from abc import ABC, abstractmethod
from plasma_box.constants import (
    Q_E, M_P, C,
    DIPOLE_MOMENT, DIPOLE_TILT_DEG,
)


# ── Abstract base ─────────────────────────────────────────────────────────────

class FieldModel(ABC):
    """Abstract electromagnetic field model.

    Subclasses implement :meth:`evaluate`.  Calling the instance directly
    is equivalent to calling :meth:`evaluate`.

    Class attributes
    ----------------
    is_uniform : bool
        True if the field is spatially uniform (E and B constant in space).
        Used by field-line plotters to skip streamline tracing and draw
        representative arrows instead.
    is_static : bool
        True if the field does not depend on time.
        Used by animation drivers to avoid redundant re-tracing.
    """

    is_uniform: bool = False
    is_static:  bool = True

    @abstractmethod
    def evaluate(self, r: np.ndarray, t: float):
        """Return (B, E) arrays of shape (N, 3) at positions r (N, 3), time t."""

    def __call__(self, r: np.ndarray, t: float):
        return self.evaluate(np.atleast_2d(r), float(t))


# ── Concrete fields ───────────────────────────────────────────────────────────

class UniformB(FieldModel):
    """Homogeneous, static magnetic field; zero electric field.

    Parameters
    ----------
    B : array-like, shape (3,)
        Magnetic field vector [T].  Default = (0, 0, 3e-10) T (3 µG along z).
    """

    is_uniform = True
    is_static  = True

    def __init__(self, B=(0.0, 0.0, 3e-10)):
        self._B = np.asarray(B, dtype=float)  # shape (3,)

    def evaluate(self, r: np.ndarray, t: float):
        N = r.shape[0]
        B = np.broadcast_to(self._B, (N, 3)).copy()
        E = np.zeros((N, 3))
        return B, E


class UniformEB(FieldModel):
    """Homogeneous, static magnetic and electric field.

    Parameters
    ----------
    B : array-like, shape (3,)
        Magnetic field vector [T].
    E : array-like, shape (3,)
        Electric field vector [V/m].
    """

    is_uniform = True
    is_static  = True

    def __init__(self, B=(0.0, 0.0, 3e-10), E=None):
        self._B = np.asarray(B, dtype=float)
        Bz = self._B[2]
        # Default E matches const_EB from source: Ex=Ey=Bz*1e4, Ez=0
        self._E = np.asarray(E if E is not None else [Bz * 1e4, Bz * 1e4, 0.0],
                             dtype=float)

    def evaluate(self, r: np.ndarray, t: float):
        N = r.shape[0]
        B = np.broadcast_to(self._B, (N, 3)).copy()
        E = np.broadcast_to(self._E, (N, 3)).copy()
        return B, E


class CyclotronWaveField(FieldModel):
    """Static uniform B plus a sinusoidal E wave at the cyclotron frequency.

    The resonance condition ω_wave = ω_c = qB/m is enforced automatically.
    B is spatially uniform and static; E is spatially uniform but time-varying.

    Parameters
    ----------
    B : array-like, shape (3,)
        Background magnetic field vector [T].
    q : float
        Particle charge [C].  Used to compute ω_c.
    m : float
        Particle mass [kg].  Used to compute ω_c.
    E_amp : float
        Amplitude of the oscillating electric field [V/m].
    E_axis : int
        Axis index (0=x, 1=y, 2=z) along which the wave oscillates.
    """

    is_uniform = True
    is_static  = False

    def __init__(
        self,
        B=(0.0, 0.0, 3e-10),
        q: float = Q_E,
        m: float = M_P,
        E_amp: float = 5e-5,
        E_axis: int = 1,
    ):
        self._B = np.asarray(B, dtype=float)
        Bmag = np.linalg.norm(self._B)
        self.omega_c = abs(q) * Bmag / m   # cyclotron angular frequency [rad/s]
        self._E_amp = float(E_amp)
        self._E_axis = int(E_axis)

    def evaluate(self, r: np.ndarray, t: float):
        N = r.shape[0]
        B = np.broadcast_to(self._B, (N, 3)).copy()
        E = np.zeros((N, 3))
        E[:, self._E_axis] = self._E_amp * np.sin(self.omega_c * t)
        return B, E


class EarthDipole(FieldModel):
    """Tilted magnetic dipole field of planet Earth; zero electric field.

    The field is derived from the dipole approximation:

        B_x = -M * (3xz cosφ + 3xy sinφ) / r^5
        B_y = -M * (3yz cosφ + (2y²-x²-z²) sinφ) / r^5
        B_z = -M * ((2z²-x²-y²) cosφ + 3zy sinφ) / r^5

    where M = DIPOLE_MOMENT, φ = tilt angle, r = |position|.

    Parameters
    ----------
    tilt_deg : float
        Angle between magnetic axis and Earth's rotation axis [degrees].
    moment : float
        Dipole moment coefficient [T·m³].
    """

    is_uniform = False
    is_static  = True

    def __init__(self, tilt_deg: float = DIPOLE_TILT_DEG, moment: float = DIPOLE_MOMENT):
        phi = np.deg2rad(tilt_deg)
        self._sin_phi = np.sin(phi)
        self._cos_phi = np.cos(phi)
        self._M = float(moment)

    def evaluate(self, r: np.ndarray, t: float):
        x = r[:, 0]
        y = r[:, 1]
        z = r[:, 2]
        r5 = (x**2 + y**2 + z**2) ** 2.5   # |r|^5, shape (N,)
        sp = self._sin_phi
        cp = self._cos_phi
        M  = self._M

        Bx = -M * (3*x*z*cp + 3*x*y*sp) / r5
        By = -M * (3*y*z*cp + (2*y**2 - x**2 - z**2)*sp) / r5
        Bz = -M * ((2*z**2 - x**2 - y**2)*cp + 3*z*y*sp) / r5

        B = np.stack([Bx, By, Bz], axis=1)   # (N, 3)
        E = np.zeros_like(B)
        return B, E


class CustomField(FieldModel):
    """Adapter wrapping a user-supplied field function.

    Two calling conventions are supported:

    **Legacy / scalar API** (``vector_api=False``)::

        def my_field(x, y, z, t):
            return Bx, By, Bz, Ex, Ey, Ez   # scalars or length-N arrays

    **Vector API** (``vector_api=True``)::

        def my_field(r, t):                  # r shape (N, 3)
            return B, E                       # each shape (N, 3)

    Parameters
    ----------
    func : callable
        The user's field function.
    vector_api : bool
        True → func uses the (r, t) → (B, E) convention.
        False → func uses the (x, y, z, t) → 6-tuple scalar convention.
    """

    def __init__(self, func, vector_api: bool = False,
                 is_uniform: bool = False, is_static: bool = True):
        self._func = func
        self._vector_api = vector_api
        self.is_uniform = is_uniform
        self.is_static  = is_static

    def evaluate(self, r: np.ndarray, t: float):
        if self._vector_api:
            B, E = self._func(r, t)
            return np.atleast_2d(B).astype(float), np.atleast_2d(E).astype(float)
        else:
            # Legacy scalar adapter
            x, y, z = r[:, 0], r[:, 1], r[:, 2]
            result = self._func(x, y, z, t)
            Bx, By, Bz, Ex, Ey, Ez = result
            B = np.stack([np.broadcast_to(Bx, x.shape),
                          np.broadcast_to(By, x.shape),
                          np.broadcast_to(Bz, x.shape)], axis=1).astype(float)
            E = np.stack([np.broadcast_to(Ex, x.shape),
                          np.broadcast_to(Ey, x.shape),
                          np.broadcast_to(Ez, x.shape)], axis=1).astype(float)
            return B, E
