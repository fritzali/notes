"""
plasma_box/particles.py
-----------------------
Particle state container and initial-condition generators.

Design note on shapes
~~~~~~~~~~~~~~~~~~~~~
All array quantities use shape **(N, 3)** for position/velocity even when N=1.
Scalar particle properties (q, m) use shape **(N,)** and broadcast against
(N, 3) arrays as ``q[:, None]``.  This makes N>1 extension trivial.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional
from plasma_box.constants import C, Q_E, M_P, M_E


@dataclass
class ParticleState:
    """State of N particles at a single instant.

    Attributes
    ----------
    q : ndarray, shape (N,)
        Charges [C].
    m : ndarray, shape (N,)
        Masses [kg].
    r : ndarray, shape (N, 3)
        Positions [m].
    v : ndarray, shape (N, 3)
        Velocities [m/s].
    t : float
        Current simulation time [s].
    """

    q: np.ndarray
    m: np.ndarray
    r: np.ndarray
    v: np.ndarray
    t: float = 0.0

    def __post_init__(self):
        # Ensure correct dtypes and shapes
        self.q = np.atleast_1d(np.asarray(self.q, dtype=float))   # (N,)
        self.m = np.atleast_1d(np.asarray(self.m, dtype=float))   # (N,)
        self.r = np.atleast_2d(np.asarray(self.r, dtype=float))   # (N,3)
        self.v = np.atleast_2d(np.asarray(self.v, dtype=float))   # (N,3)
        N = self.q.shape[0]
        assert self.m.shape == (N,), "m must have shape (N,)"
        assert self.r.shape == (N, 3), "r must have shape (N, 3)"
        assert self.v.shape == (N, 3), "v must have shape (N, 3)"
        self.t = float(self.t)

    @property
    def N(self) -> int:
        """Number of particles."""
        return self.q.shape[0]

    # ── Kinematic helpers ──────────────────────────────────────────────────────

    def speed(self) -> np.ndarray:
        """Speed |v| of each particle, shape (N,)."""
        return np.linalg.norm(self.v, axis=1)

    def gamma(self) -> np.ndarray:
        """Lorentz factor γ of each particle, shape (N,).

        Returns exactly 1.0 where v²/c² is negligibly small (< 1e-15),
        guarding against floating-point issues.
        """
        beta2 = (self.speed() / C) ** 2
        beta2 = np.clip(beta2, 0.0, 1.0 - 1e-15)
        return 1.0 / np.sqrt(1.0 - beta2)

    def kinetic_energy(self, relativistic: bool = True) -> np.ndarray:
        """Kinetic energy of each particle [J], shape (N,).

        Parameters
        ----------
        relativistic : bool
            If True, use K = (γ-1)mc².  If False, use K = ½mv².
        """
        if relativistic:
            return (self.gamma() - 1.0) * self.m * C**2
        else:
            return 0.5 * self.m * np.sum(self.v**2, axis=1)

    def copy(self) -> "ParticleState":
        """Return a deep copy of this state."""
        return ParticleState(
            q=self.q.copy(),
            m=self.m.copy(),
            r=self.r.copy(),
            v=self.v.copy(),
            t=self.t,
        )


# ── Initial-condition factories ────────────────────────────────────────────────

def single_particle(
    q: float = Q_E,
    m: float = M_P,
    r0=(0.0, 0.0, 0.0),
    v0=(0.0, 0.0, 0.0),
    t0: float = 0.0,
) -> ParticleState:
    """Create a single-particle state with explicit initial conditions.

    Parameters
    ----------
    q : float
        Charge [C].
    m : float
        Mass [kg].
    r0 : array-like, shape (3,)
        Initial position [m].
    v0 : array-like, shape (3,)
        Initial velocity [m/s].
    t0 : float
        Initial time [s].

    Returns
    -------
    ParticleState with N=1.
    """
    return ParticleState(
        q=np.array([q]),
        m=np.array([m]),
        r=np.array([r0], dtype=float),
        v=np.array([v0], dtype=float),
        t=t0,
    )


def random_ensemble(
    N: int,
    q: float = Q_E,
    m: float = M_P,
    r0=(0.0, 0.0, 0.0),
    v_max: float = 0.01 * C,
    axes: tuple = (1, 2),
    seed: Optional[int] = None,
    t0: float = 0.0,
) -> ParticleState:
    """Create N particles with randomised initial velocities.

    Velocity components along the specified *axes* are drawn uniformly from
    (-v_max, v_max).  All other velocity components are zero.

    Reproduces the ``PVREDNOSTI.puslovi`` pattern used in the Verlet and
    cyclotron-resonance examples.

    Parameters
    ----------
    N : int
        Number of particles.
    q, m : float
        Charge [C] and mass [kg] (same for all particles).
    r0 : array-like, shape (3,)
        Common initial position [m].
    v_max : float
        Velocity half-range [m/s].
    axes : tuple of int
        Which velocity components (0=x, 1=y, 2=z) are randomised.
    seed : int or None
        RNG seed for reproducibility.
    t0 : float
        Initial time [s].

    Returns
    -------
    ParticleState with shape (N,) / (N, 3).
    """
    rng = np.random.default_rng(seed)
    r = np.tile(np.asarray(r0, dtype=float), (N, 1))
    v = np.zeros((N, 3))
    for ax in axes:
        v[:, ax] = (rng.random(N) - 0.5) * 2.0 * v_max
    return ParticleState(
        q=np.full(N, q),
        m=np.full(N, m),
        r=r,
        v=v,
        t=t0,
    )


def dipole_initial_conditions(
    q: float = Q_E,
    m: float = M_P,
    r0=None,
    v0=None,
    t0: float = 0.0,
) -> ParticleState:
    """Default initial conditions matching the relativistic RK4 Earth-dipole example.

    Position: (2.5 R_E, 0, 0).
    Velocity: relativistic proton at ~0.616 c pitched 60° from y toward z.

    Override r0 and v0 to use different values.
    """
    from plasma_box.constants import R_EARTH
    if r0 is None:
        r0 = [2.5 * R_EARTH, 0.0, 0.0]
    if v0 is None:
        v0 = [0.0, 0.616 * 0.5 * C, 0.616 * 0.866 * C]
    return single_particle(q=q, m=m, r0=r0, v0=v0, t0=t0)
