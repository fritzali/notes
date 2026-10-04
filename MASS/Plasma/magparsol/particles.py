"""
Particle states and initial conditions.

Positions and velocities always have shape (N, 3), also for a single
particle; charges and masses have shape (N,) and broadcast as ``q[:, None]``.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional
from magparsol.constants import C, Q_E, M_P, M_E


@dataclass
class ParticleState:
    """Charges, masses, positions and velocities of N particles at time t.

    Attributes
    ----------
    q, m : ndarray, shape (N,)
        Charges [C] and masses [kg].
    r, v : ndarray, shape (N, 3)
        Positions [m] and velocities [m/s].
    t : float
        Time [s].
    """

    q: np.ndarray
    m: np.ndarray
    r: np.ndarray
    v: np.ndarray
    t: float = 0.0

    def __post_init__(self):
        self.q = np.atleast_1d(np.asarray(self.q, dtype=float))
        self.m = np.atleast_1d(np.asarray(self.m, dtype=float))
        self.r = np.atleast_2d(np.asarray(self.r, dtype=float))
        self.v = np.atleast_2d(np.asarray(self.v, dtype=float))
        N = self.q.shape[0]
        assert self.m.shape == (N,), "m must have shape (N,)"
        assert self.r.shape == (N, 3), "r must have shape (N, 3)"
        assert self.v.shape == (N, 3), "v must have shape (N, 3)"
        self.t = float(self.t)

    @property
    def N(self) -> int:
        """Number of particles."""
        return self.q.shape[0]

    def speed(self) -> np.ndarray:
        """Speeds |v|, shape (N,)."""
        return np.linalg.norm(self.v, axis=1)

    def gamma(self) -> np.ndarray:
        """Lorentz factors, shape (N,).

        β² is capped just below 1 so that γ stays finite for v ≥ c.
        """
        beta2 = (self.speed() / C) ** 2
        beta2 = np.clip(beta2, 0.0, 1.0 - 1e-15)
        return 1.0 / np.sqrt(1.0 - beta2)

    def kinetic_energy(self, relativistic: bool = True) -> np.ndarray:
        """Kinetic energies [J], shape (N,): (γ - 1) m c², or m v²/2."""
        if relativistic:
            return (self.gamma() - 1.0) * self.m * C**2
        else:
            return 0.5 * self.m * np.sum(self.v**2, axis=1)

    def copy(self) -> "ParticleState":
        """Independent copy (integrators modify their state in place)."""
        return ParticleState(
            q=self.q.copy(),
            m=self.m.copy(),
            r=self.r.copy(),
            v=self.v.copy(),
            t=self.t,
        )


def single_particle(
    q: float = Q_E,
    m: float = M_P,
    r0=(0.0, 0.0, 0.0),
    v0=(0.0, 0.0, 0.0),
    t0: float = 0.0,
) -> ParticleState:
    """One particle with charge ``q`` [C], mass ``m`` [kg], position ``r0`` [m],
    velocity ``v0`` [m/s] and start time ``t0`` [s]."""
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
    """N identical particles at ``r0`` with uniformly random velocities.

    The velocity components listed in ``axes`` (0, 1, 2 for x, y, z) are
    drawn uniformly from (-v_max, v_max); the others are zero.

    Parameters
    ----------
    N : int
    q, m : float
        Charge [C] and mass [kg] of every particle.
    r0 : array_like, shape (3,)
        Common start position [m].
    v_max : float
        Largest velocity component [m/s].
    axes : tuple of int
    seed : int or None
        Seed of the random generator.
    t0 : float
        Start time [s].
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


def maxwellian_ensemble(
    N: int,
    T: float,
    m: float = M_P,
    q: float = Q_E,
    r0=(0.0, 0.0, 0.0),
    seed: Optional[int] = None,
    t0: float = 0.0,
) -> ParticleState:
    """N particles at ``r0`` with an isotropic, nonrelativistic Maxwellian.

    Each velocity component is Gaussian with standard deviation
    sqrt(k T / m). ``T`` is the temperature in kelvin; the other parameters
    are as in :func:`random_ensemble`.
    """
    from magparsol.constants import K_B as k_B
    rng  = np.random.default_rng(seed)
    sigma = np.sqrt(k_B * T / m)
    v    = rng.normal(0.0, sigma, size=(N, 3))
    r    = np.tile(np.asarray(r0, dtype=float), (N, 1))
    return ParticleState(q=np.full(N, q), m=np.full(N, m), r=r, v=v, t=t0)


def relativistic_thermal_ensemble(
    N: int,
    theta_e: float,
    m: float = M_P,
    q: float = Q_E,
    r0=(0.0, 0.0, 0.0),
    seed: Optional[int] = None,
    t0: float = 0.0,
) -> ParticleState:
    """N particles at ``r0`` with an isotropic Maxwell–Jüttner distribution.

    ``theta_e`` = k T / (m c²) is the dimensionless temperature. The momentum
    u = γβ follows f(u) ∝ u² exp(-γ/θ); it is drawn by rejection sampling
    from an exponential proposal, and the direction is uniform on the sphere.
    The other parameters are as in :func:`random_ensemble`.
    """
    rng  = np.random.default_rng(seed)
    # Proposal exp(-u/scale); the acceptance weight f(u)/proposal is
    # u² exp(-γ/θ + u/scale), normalised to its largest value in each batch.
    scale   = theta_e + 1.5
    u_samp  = []
    while len(u_samp) < N:
        batch = int((N - len(u_samp)) * 4 + 10)
        u_prop = rng.exponential(scale, size=batch)
        gamma  = np.sqrt(1.0 + u_prop**2)
        log_w  = 2*np.log(u_prop + 1e-30) - gamma/theta_e + u_prop/scale
        log_w -= log_w.max()
        accept = np.log(rng.uniform(size=batch)) < log_w
        u_samp.extend(u_prop[accept].tolist())
    u_mag = np.array(u_samp[:N])

    phi   = rng.uniform(0, 2*np.pi, N)
    costh = rng.uniform(-1, 1, N)
    sinth = np.sqrt(1 - costh**2)
    u_vec = np.column_stack([
        u_mag * sinth * np.cos(phi),
        u_mag * sinth * np.sin(phi),
        u_mag * costh,
    ])
    gamma = np.sqrt(1.0 + u_mag**2)
    v     = u_vec * C / gamma[:, None]   # v = c u / γ

    r = np.tile(np.asarray(r0, dtype=float), (N, 1))
    return ParticleState(q=np.full(N, q), m=np.full(N, m), r=r, v=v, t=t0)


def dipole_initial_conditions(
    q: float = Q_E,
    m: float = M_P,
    r0=None,
    v0=None,
    t0: float = 0.0,
) -> ParticleState:
    """A particle for the Earth dipole, by default a proton at 2.5 R_E.

    Defaults: r0 = (2.5 R_E, 0, 0) and |v0| = 0.616 c in the yz plane,
    60° from y towards z. Pass ``r0`` and ``v0`` to start elsewhere.
    """
    from magparsol.constants import R_EARTH
    if r0 is None:
        r0 = [2.5 * R_EARTH, 0.0, 0.0]
    if v0 is None:
        v0 = [0.0, 0.616 * 0.5 * C, 0.616 * 0.866 * C]
    return single_particle(q=q, m=m, r0=r0, v0=v0, t0=t0)
