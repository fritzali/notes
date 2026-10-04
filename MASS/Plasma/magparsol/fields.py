"""
Electromagnetic field models.

A field model is called as ``field(r, t)`` with positions ``r`` of shape
(N, 3) in metres and a scalar time ``t`` in seconds, and returns ``(B, E)``,
two (N, 3) arrays in tesla and volt per metre. Subclasses implement
:meth:`FieldModel.evaluate`.

Any other field can be wrapped with :class:`CustomField`, either from a
scalar function ``f(x, y, z, t) -> (Bx, By, Bz, Ex, Ey, Ez)`` or from a
vector function ``f(r, t) -> (B, E)``.
"""

import numpy as np
from abc import ABC, abstractmethod
from magparsol.constants import (
    Q_E, M_P, C,
    DIPOLE_MOMENT, DIPOLE_TILT_DEG,
)


class FieldModel(ABC):
    """Base class of all field models.

    Attributes
    ----------
    is_uniform : bool
        The field does not depend on position. Field line plots then draw
        a grid of arrows instead of tracing lines.
    is_static : bool
        The field does not depend on time. Animations only redraw the field
        of time dependent models.
    display_name : str
        Name used in default plot titles.
    """

    is_uniform: bool = False
    is_static:  bool = True
    display_name: str = "Field"

    @abstractmethod
    def evaluate(self, r: np.ndarray, t: float):
        """Return ``(B, E)``, each of shape (N, 3), at positions ``r`` (N, 3)."""

    def __call__(self, r: np.ndarray, t: float):
        return self.evaluate(np.atleast_2d(r), float(t))


class UniformB(FieldModel):
    """Constant magnetic field, no electric field.

    Parameters
    ----------
    B : array_like, shape (3,)
        Field vector [T]. The default, 3e-10 T along z, is 3 µG.
    """

    is_uniform = True
    is_static  = True
    display_name = r"Uniform $\mathbf{B}$"

    def __init__(self, B=(0.0, 0.0, 3e-10)):
        self._B = np.asarray(B, dtype=float)

    def evaluate(self, r: np.ndarray, t: float):
        N = r.shape[0]
        B = np.broadcast_to(self._B, (N, 3)).copy()
        E = np.zeros((N, 3))
        return B, E


class UniformEB(FieldModel):
    """Constant magnetic and electric field.

    Parameters
    ----------
    B : array_like, shape (3,)
        Magnetic field [T].
    E : array_like, shape (3,) or None
        Electric field [V/m]. By default E = B_z (1, 1, 0) × 10⁴ m/s, which
        gives an E×B drift of √2 × 10⁴ m/s.
    """

    is_uniform = True
    is_static  = True
    display_name = "Crossed Fields"

    def __init__(self, B=(0.0, 0.0, 3e-10), E=None):
        self._B = np.asarray(B, dtype=float)
        Bz = self._B[2]
        self._E = np.asarray(E if E is not None else [Bz * 1e4, Bz * 1e4, 0.0],
                             dtype=float)

    def evaluate(self, r: np.ndarray, t: float):
        N = r.shape[0]
        B = np.broadcast_to(self._B, (N, 3)).copy()
        E = np.broadcast_to(self._E, (N, 3)).copy()
        return B, E


class CyclotronWaveField(FieldModel):
    """Constant B plus a uniform electric field oscillating at the gyrofrequency.

    E points along one axis, E = E_amp sin(ω_c t), with the nonrelativistic
    gyrofrequency ω_c = |q| B / m of the given particle species. The wave is
    therefore resonant with particles of that species as long as v ≪ c.

    Parameters
    ----------
    B : array_like, shape (3,)
        Background magnetic field [T].
    q, m : float
        Charge [C] and mass [kg] that set ω_c.
    E_amp : float
        Amplitude of the electric field [V/m].
    E_axis : int
        Direction of E: 0, 1 or 2 for x, y or z.

    Attributes
    ----------
    omega_c : float
        Wave angular frequency [rad/s].
    """

    is_uniform = True
    is_static  = False
    display_name = "Cyclotron Wave"

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
        self.omega_c = abs(q) * Bmag / m
        self._E_amp = float(E_amp)
        self._E_axis = int(E_axis)

    def evaluate(self, r: np.ndarray, t: float):
        N = r.shape[0]
        B = np.broadcast_to(self._B, (N, 3)).copy()
        E = np.zeros((N, 3))
        E[:, self._E_axis] = self._E_amp * np.sin(self.omega_c * t)
        return B, E


class EarthDipole(FieldModel):
    """Centred, tilted dipole model of the geomagnetic field; no electric field.

    The dipole axis m̂ = (0, sin φ, cos φ) is tilted by φ from the rotation
    axis ẑ towards ŷ, and

        B = -M [3 (m̂·r) r - m̂ r²] / r⁵,

    with M = μ0 m / 4π. The sign makes B point north (along m̂) at the
    magnetic equator, as for the Earth. The defaults are the IGRF-14 values
    from :mod:`magparsol.constants`.

    Parameters
    ----------
    tilt_deg : float
        Tilt φ of the magnetic axis [degrees].
    moment : float
        M = μ0 m / 4π [T m³]; M / r³ is the equatorial field at distance r.

    Attributes
    ----------
    axis : ndarray, shape (3,)
        Unit vector m̂ along the magnetic axis.
    tilt_deg : float
    """

    is_uniform = False
    is_static  = True
    display_name = "Earth Dipole"

    def __init__(self, tilt_deg: float = DIPOLE_TILT_DEG, moment: float = DIPOLE_MOMENT):
        phi = np.deg2rad(tilt_deg)
        self.tilt_deg = float(tilt_deg)
        self._sin_phi = np.sin(phi)
        self._cos_phi = np.cos(phi)
        self._M = float(moment)
        self.axis = np.array([0.0, self._sin_phi, self._cos_phi])

    def B_equator(self, r: float) -> float:
        """Field strength on the magnetic equator at distance ``r`` [T]."""
        return abs(self._M) / r**3

    def evaluate(self, r: np.ndarray, t: float):
        x = r[:, 0]
        y = r[:, 1]
        z = r[:, 2]
        r5 = (x**2 + y**2 + z**2) ** 2.5
        sp = self._sin_phi
        cp = self._cos_phi
        M  = self._M

        # Components of -M [3 (m·r) r - m r²] / r⁵ with m = (0, sp, cp)
        Bx = -M * (3*x*z*cp + 3*x*y*sp) / r5
        By = -M * (3*y*z*cp + (2*y**2 - x**2 - z**2)*sp) / r5
        Bz = -M * ((2*z**2 - x**2 - y**2)*cp + 3*z*y*sp) / r5

        B = np.stack([Bx, By, Bz], axis=1)
        E = np.zeros_like(B)
        return B, E


class CustomField(FieldModel):
    """Field model from a user function.

    With ``vector_api=False`` the function is ``f(x, y, z, t)`` and returns
    the six components ``(Bx, By, Bz, Ex, Ey, Ez)``; each may be a scalar or
    an array of length N. With ``vector_api=True`` it is ``f(r, t)`` with
    ``r`` of shape (N, 3) and returns ``(B, E)`` of shape (N, 3).

    Parameters
    ----------
    func : callable
    vector_api : bool
    is_uniform, is_static : bool
        See :class:`FieldModel`. They only affect plotting and animation.
    name : str or None
        Name used in plot titles.
    """

    display_name = "Custom Field"

    def __init__(self, func, vector_api: bool = False,
                 is_uniform: bool = False, is_static: bool = True,
                 name: str = None):
        if name is not None:
            self.display_name = name
        self._func = func
        self._vector_api = vector_api
        self.is_uniform = is_uniform
        self.is_static  = is_static

    def evaluate(self, r: np.ndarray, t: float):
        if self._vector_api:
            B, E = self._func(r, t)
            return np.atleast_2d(B).astype(float), np.atleast_2d(E).astype(float)
        else:
            x, y, z = r[:, 0], r[:, 1], r[:, 2]
            result = self._func(x, y, z, t)
            Bx, By, Bz, Ex, Ey, Ez = result
            # Broadcast scalar components to length N
            B = np.stack([np.broadcast_to(Bx, x.shape),
                          np.broadcast_to(By, x.shape),
                          np.broadcast_to(Bz, x.shape)], axis=1).astype(float)
            E = np.stack([np.broadcast_to(Ex, x.shape),
                          np.broadcast_to(Ey, x.shape),
                          np.broadcast_to(Ez, x.shape)], axis=1).astype(float)
            return B, E
