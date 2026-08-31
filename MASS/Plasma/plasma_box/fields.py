"""
plasma_box/fields.py
--------------------
Electromagnetic field models.

Every FieldModel subclass implements:

    evaluate(r, t) -> B

where r is an (N, 3) position array and t is a scalar time,
and the return value B is an (N, 3) array in SI units (Tesla).

A CustomField adapter lets users supply a legacy scalar-signature function::

    def my_field(x, y, z, t):
        return Bx, By, Bz                   # scalars or length-N arrays

or a modern vector-signature function::

    def my_field(r, t):                     # r shape (N,3)
        return B                            # shape (N,3)

and wrap it::

    field = CustomField(my_func, vector_api=False)   # legacy
    field = CustomField(my_func, vector_api=True)    # modern
"""

import numpy as np
from abc import ABC, abstractmethod
from plasma_box.constants import (
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
        True if the field is spatially uniform (B constant in space).
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
        """Return B array of shape (N, 3) at positions r (N, 3), time t."""

    def __call__(self, r: np.ndarray, t: float):
        return self.evaluate(np.atleast_2d(r), float(t))


# ── Concrete fields ───────────────────────────────────────────────────────────

class UniformB(FieldModel):
    """Homogeneous, static magnetic field.

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
        return B


class EarthDipole(FieldModel):
    """Tilted magnetic dipole field of planet Earth.

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
        return B


class CustomField(FieldModel):
    """Adapter wrapping a user-supplied field function.

    Two calling conventions are supported:

    **Legacy / scalar API** (``vector_api=False``)::

        def my_field(x, y, z, t):
            return Bx, By, Bz               # scalars or length-N arrays

    **Vector API** (``vector_api=True``)::

        def my_field(r, t):                  # r shape (N, 3)
            return B                         # shape (N, 3)

    Parameters
    ----------
    func : callable
        The user's field function.
    vector_api : bool
        True → func uses the (r, t) → B convention.
        False → func uses the (x, y, z, t) → 3-tuple scalar convention.
    """

    def __init__(self, func, vector_api: bool = False,
                 is_uniform: bool = False, is_static: bool = True):
        self._func = func
        self._vector_api = vector_api
        self.is_uniform = is_uniform
        self.is_static  = is_static

    def evaluate(self, r: np.ndarray, t: float):
        if self._vector_api:
            B = self._func(r, t)
            return np.atleast_2d(B).astype(float)
        else:
            # Legacy scalar adapter
            x, y, z = r[:, 0], r[:, 1], r[:, 2]
            result = self._func(x, y, z, t)
            Bx, By, Bz = result
            B = np.stack([np.broadcast_to(Bx, x.shape),
                          np.broadcast_to(By, x.shape),
                          np.broadcast_to(Bz, x.shape)], axis=1).astype(float)
            return B
