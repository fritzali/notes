"""
Trajectory history and gyration diagnostics.

:class:`TrajectoryHistory` collects all states visited during a run; the
functions below compute gyration quantities and check time steps.
"""

import numpy as np
import warnings
from magparsol.constants import C


class TrajectoryHistory:
    """Stored states of a run.

    States are collected in lists while the integrator runs (adaptive
    integrators do not know the number of steps in advance) and converted
    to arrays by :meth:`finalize`. Afterwards the history has

    * ``t``: times, shape (S,)
    * ``r``: positions, shape (S, N, 3)
    * ``v``: velocities, shape (S, N, 3)

    for S stored samples of N particles.

    Parameters
    ----------
    store_dt : float or None
        Sampling interval [s]. States are stored on the grid t0 + k store_dt
        (at the first step on or after each grid time). None stores every
        step.
    """

    def __init__(self, store_dt=None):
        self._t: list = []
        self._r: list = []
        self._v: list = []
        self.store_dt = store_dt
        self._last_stored_t: float = -np.inf
        self._next_store_t: float = -np.inf
        self._finalized: bool = False

        self.t: np.ndarray = None
        self.r: np.ndarray = None
        self.v: np.ndarray = None

    def record(self, state, force: bool = False):
        """Store ``state`` if the next sampling time has been reached.

        ``force=True`` stores it unconditionally and restarts the sampling
        grid at ``state.t``; the integrator does this for the initial state.
        """
        if self._finalized:
            raise RuntimeError("Cannot record into a finalized TrajectoryHistory.")
        if force or self.store_dt is None:
            store = True
            self._next_store_t = state.t + (self.store_dt or 0.0)
        else:
            # The small tolerance absorbs round-off in the accumulated time,
            # which would otherwise skip samples now and then.
            store = state.t >= self._next_store_t - 1e-6 * self.store_dt
            if store:
                while self._next_store_t <= state.t + 1e-6 * self.store_dt:
                    self._next_store_t += self.store_dt
        if store:
            self._t.append(state.t)
            self._r.append(state.r.copy())
            self._v.append(state.v.copy())
            self._last_stored_t = state.t

    def finalize(self):
        """Convert the stored lists to the arrays ``t``, ``r`` and ``v``."""
        if len(self._t) == 0:
            raise RuntimeError("No data recorded; run the integrator first.")
        self.t = np.array(self._t)
        self.r = np.array(self._r)
        self.v = np.array(self._v)
        self._finalized = True

    @property
    def is_finalized(self) -> bool:
        return self._finalized

    def __len__(self):
        return len(self._t) if not self._finalized else len(self.t)

    def kinetic_energy(self, m: np.ndarray, relativistic: bool = True) -> np.ndarray:
        """Kinetic energy [J], shape (S, N), for masses ``m`` of shape (N,).

        (γ - 1) m c² if ``relativistic``, otherwise m v²/2.
        """
        self._check_finalized()
        speed2 = np.sum(self.v**2, axis=2)
        if relativistic:
            beta2 = np.clip(speed2 / C**2, 0.0, 1.0 - 1e-15)
            gamma = 1.0 / np.sqrt(1.0 - beta2)
            return (gamma - 1.0) * m[None, :] * C**2
        else:
            return 0.5 * m[None, :] * speed2

    def speed(self) -> np.ndarray:
        """Speed |v|, shape (S, N)."""
        self._check_finalized()
        return np.linalg.norm(self.v, axis=2)

    def gamma(self) -> np.ndarray:
        """Lorentz factor, shape (S, N)."""
        self._check_finalized()
        beta2 = np.clip(self.speed()**2 / C**2, 0.0, 1.0 - 1e-15)
        return 1.0 / np.sqrt(1.0 - beta2)

    def _check_finalized(self):
        if not self._finalized:
            raise RuntimeError("Call finalize() before accessing trajectory arrays.")


def relative_energy_error(
    history: TrajectoryHistory,
    m: np.ndarray,
    relativistic: bool = True,
) -> np.ndarray:
    """Relative change of the kinetic energy, (K(t) - K(0)) / K(0), shape (S, N).

    In a pure magnetic field K is conserved, so this measures the error
    of the integrator. Particles with K(0) = 0 get the absolute change.
    """
    K = history.kinetic_energy(m, relativistic=relativistic)
    K0 = K[0:1, :]
    return (K - K0) / np.where(np.abs(K0) > 0, np.abs(K0), 1.0)


def gyrofrequency(q: float, m: float, B_mag: float) -> float:
    """Nonrelativistic gyrofrequency ω_c = |q| B / m [rad/s].

    For a relativistic particle use ``m`` = γ m₀.
    """
    return abs(q) * abs(B_mag) / m


def gyroperiod(q: float, m: float, B_mag: float) -> float:
    """Nonrelativistic gyroperiod T_c = 2π m / (|q| B) [s]."""
    return 2.0 * np.pi * m / (abs(q) * abs(B_mag))


def gyroradius(m: float, v_perp: float, q: float, B_mag: float) -> float:
    """Nonrelativistic gyroradius r_L = m v⊥ / (|q| B) [m]."""
    return m * abs(v_perp) / (abs(q) * abs(B_mag))


def check_dt_resolution(
    dt: float,
    q: float,
    m: float,
    B_mag: float,
    warn_threshold: float = 0.1,
) -> float:
    """Return dt / T_c and warn if it exceeds ``warn_threshold``.

    Parameters
    ----------
    dt : float
        Time step [s].
    q, m : float
        Charge [C] and mass [kg].
    B_mag : float
        Typical field strength [T].
    warn_threshold : float
    """
    Tc = gyroperiod(q, m, B_mag)
    ratio = dt / Tc
    if ratio > warn_threshold:
        warnings.warn(
            f"dt/T_c = {ratio:.3f} > {warn_threshold}. "
            "Consider reducing dt for accurate gyration resolution.",
            UserWarning,
            stacklevel=2,
        )
    return ratio


def suggest_dt(
    q: float,
    m: float,
    B_mag: float,
    steps_per_gyration: float = 100.0,
) -> float:
    """Time step [s] that resolves one gyroperiod with ``steps_per_gyration`` steps."""
    Tc = gyroperiod(q, m, B_mag)
    return Tc / steps_per_gyration
