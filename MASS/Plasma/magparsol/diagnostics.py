"""
magparsol/diagnostics.py
-------------------------
Trajectory storage and physical diagnostics.

TrajectoryHistory
~~~~~~~~~~~~~~~~~
Stores particle trajectories from integrator runs using dynamic list
accumulation (required for adaptive step-size integrators where the number
of steps is not known in advance).  Call ``finalize()`` to convert lists to
numpy arrays after the run.

Diagnostics
~~~~~~~~~~~
Standalone functions that operate on a ``TrajectoryHistory`` or on raw
``ParticleState`` data.
"""

import numpy as np
import warnings
from magparsol.constants import C


class TrajectoryHistory:
    """Dynamic-sized storage for particle trajectory data.

    Supports both fixed-step integrators (where every step is stored) and
    adaptive-step integrators (where the time spacing is irregular).  A
    ``store_dt`` threshold can be used to limit storage to a sub-sampled
    sequence, matching the original ``T_smp`` sampling-interval pattern.

    Parameters
    ----------
    store_dt : float or None
        Minimum elapsed simulation time between stored snapshots.
        ``None`` → store every call to :meth:`record`.
    """

    def __init__(self, store_dt=None):
        self._t: list = []
        self._r: list = []
        self._v: list = []
        self.store_dt = store_dt
        self._last_stored_t: float = -np.inf
        self._finalized: bool = False

        # Arrays populated after finalize()
        self.t: np.ndarray = None
        self.r: np.ndarray = None
        self.v: np.ndarray = None

    # ── Recording ─────────────────────────────────────────────────────────────

    def record(self, state, force: bool = False):
        """Append the current particle state if the sampling condition is met.

        Parameters
        ----------
        state : ParticleState
            Current simulation state.
        force : bool
            If True, bypass the ``store_dt`` gate (used for the initial point).
        """
        if self._finalized:
            raise RuntimeError("Cannot record into a finalized TrajectoryHistory.")
        elapsed = state.t - self._last_stored_t
        if force or self.store_dt is None or elapsed >= self.store_dt:
            self._t.append(state.t)
            self._r.append(state.r.copy())   # (N, 3)
            self._v.append(state.v.copy())   # (N, 3)
            self._last_stored_t = state.t

    def finalize(self):
        """Convert internal lists to numpy arrays.

        After calling this, ``self.t``, ``self.r``, ``self.v`` are available:

        - ``t`` : shape (S,)
        - ``r`` : shape (S, N, 3)
        - ``v`` : shape (S, N, 3)

        where S is the number of stored snapshots and N is the particle count.
        """
        if len(self._t) == 0:
            raise RuntimeError("No data recorded; run the integrator first.")
        self.t = np.array(self._t)
        self.r = np.array(self._r)   # (S, N, 3)
        self.v = np.array(self._v)   # (S, N, 3)
        self._finalized = True

    @property
    def is_finalized(self) -> bool:
        return self._finalized

    def __len__(self):
        return len(self._t) if not self._finalized else len(self.t)

    # ── Derived quantities ────────────────────────────────────────────────────

    def kinetic_energy(self, m: np.ndarray, relativistic: bool = True) -> np.ndarray:
        """Kinetic energy along the trajectory.

        Parameters
        ----------
        m : ndarray, shape (N,)
            Particle masses [kg].
        relativistic : bool
            Whether to use the relativistic formula K=(γ-1)mc².

        Returns
        -------
        K : ndarray, shape (S, N) — kinetic energy of each particle at each
            stored time step [J].
        """
        self._check_finalized()
        speed2 = np.sum(self.v**2, axis=2)   # (S, N)
        if relativistic:
            beta2 = np.clip(speed2 / C**2, 0.0, 1.0 - 1e-15)
            gamma = 1.0 / np.sqrt(1.0 - beta2)   # (S, N)
            return (gamma - 1.0) * m[None, :] * C**2
        else:
            return 0.5 * m[None, :] * speed2

    def speed(self) -> np.ndarray:
        """Speed |v| at each snapshot, shape (S, N)."""
        self._check_finalized()
        return np.linalg.norm(self.v, axis=2)

    def gamma(self) -> np.ndarray:
        """Lorentz factor γ at each snapshot, shape (S, N)."""
        self._check_finalized()
        beta2 = np.clip(self.speed()**2 / C**2, 0.0, 1.0 - 1e-15)
        return 1.0 / np.sqrt(1.0 - beta2)

    def _check_finalized(self):
        if not self._finalized:
            raise RuntimeError("Call finalize() before accessing trajectory arrays.")


# ── Physical diagnostic functions ─────────────────────────────────────────────

def relative_energy_error(
    history: TrajectoryHistory,
    m: np.ndarray,
    relativistic: bool = True,
) -> np.ndarray:
    """Relative kinetic-energy deviation Δε/ε₀ over time.

    For pure-magnetic-field simulations the kinetic energy should be constant.
    Deviations indicate integrator error.

    Parameters
    ----------
    history : TrajectoryHistory (finalized)
    m : ndarray, shape (N,)
    relativistic : bool

    Returns
    -------
    err : ndarray, shape (S, N)
        ``(K(t) - K(0)) / K(0)``.
    """
    K = history.kinetic_energy(m, relativistic=relativistic)   # (S, N)
    K0 = K[0:1, :]                                              # (1, N)
    return (K - K0) / np.where(np.abs(K0) > 0, np.abs(K0), 1.0)


def gyrofrequency(q: float, m: float, B_mag: float) -> float:
    """Non-relativistic cyclotron angular frequency [rad/s].

    ω_c = |q| |B| / m
    """
    return abs(q) * abs(B_mag) / m


def gyroperiod(q: float, m: float, B_mag: float) -> float:
    """Non-relativistic gyro-period [s].

    T_c = 2π m / (|q| |B|)
    """
    return 2.0 * np.pi * m / (abs(q) * abs(B_mag))


def gyroradius(m: float, v_perp: float, q: float, B_mag: float) -> float:
    """Non-relativistic Larmor radius [m].

    r_L = m |v_⊥| / (|q| |B|)
    """
    return m * abs(v_perp) / (abs(q) * abs(B_mag))


def check_dt_resolution(
    dt: float,
    q: float,
    m: float,
    B_mag: float,
    warn_threshold: float = 0.1,
) -> float:
    """Check whether dt is small enough relative to the gyro-period.

    Parameters
    ----------
    dt : float
        Integration time step [s].
    q, m : float
        Particle charge and mass.
    B_mag : float
        Representative magnetic field magnitude [T].
    warn_threshold : float
        Emit a warning if dt / T_c > warn_threshold.

    Returns
    -------
    ratio : float
        dt / T_c.
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
    """Suggest a time step giving ``steps_per_gyration`` per gyro-period.

    Parameters
    ----------
    steps_per_gyration : float
        Target number of integration steps per gyro-period.

    Returns
    -------
    dt : float   [s]
    """
    Tc = gyroperiod(q, m, B_mag)
    return Tc / steps_per_gyration
