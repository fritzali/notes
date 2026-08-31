"""
magparsol/integrators/boris.py
--------------------------------
Relativistic Boris-family integrators.

All three variants follow Zenitani & Umeda (2018),
"On the Boris solver in particle-in-cell simulation",
Physics of Plasmas 25, 112110.

Overview of the three solvers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All three share the same leapfrog structure (operator splitting):

    1. Magnetic rotation: u⁺ = Rotate(u, B, γ)               (variant-specific)
    2. Position:          rⁿ⁺¹ = rⁿ + Δt · u⁺/γ⁺

where u = γv is the relativistic momentum per unit mass.

Boris-A (original Boris 1970, Eq. 7a)
    t = tan(θ/2) b,    θ = qΔt|B|/(mγ)
    Standard cross-product rotation via t-vector (exact angle).

Boris-B (textbook simplification, Eq. 7b)
    t = (θ/2) b  ≈  tan(θ/2) b  for small θ
    Same rotation formula; faster but introduces phase error for large θ.

Boris-C (Zenitani & Umeda 2018, Eqs. 11-12)
    u‖  = (u · b̂) b̂                                          Eq. (11)
    u⁺  = u‖ + (u − u‖) cos θ + (u × b̂) sin θ              Eq. (12)
    Exact analytic rotation — no tan approximation, second-order accurate,
    time-reversible.  Uses B_FLOOR to guard against |B|→0.

Shared rotation kernel
~~~~~~~~~~~~~~~~~~~~~~~
Boris-A and Boris-B both use the same function ``_boris_cross_rotate``
(Eqs. 8-9 of Zenitani & Umeda):

    u′ = u + u × t                                            Eq. (8)
    u⁺ = u + 2/(1+|t|²) u′ × t                               Eq. (9)

Non-relativistic limit
~~~~~~~~~~~~~~~~~~~~~~~
Setting γ = 1 in Boris-A/B recovers the non-relativistic Boris leapfrog.
The ``NonrelBorisVerlet`` class (formerly a separate file) is therefore
subsumed by using Boris-A or Boris-B with a non-relativistic particle
(v ≪ c) — γ ≈ 1 automatically.  No separate class is needed.

Note on shapes
~~~~~~~~~~~~~~
All array operations are (N, 3) so N > 1 particles advance simultaneously.
Scalar particle properties q, m have shape (N,) and broadcast via [:, None].
"""

import numpy as np
from magparsol.integrators.base import Integrator
from magparsol.constants import C, B_FLOOR


# ── Shared rotation kernel (Boris-A and Boris-B) ──────────────────────────────

def _boris_cross_rotate(u: np.ndarray, t_vec: np.ndarray) -> np.ndarray:
    """Apply the Boris cross-product rotation (Eqs. 8-9, Zenitani & Umeda 2018).

    Parameters
    ----------
    u : ndarray, shape (N, 3)
        Pre-rotation relativistic momentum per unit mass [m/s].
    t_vec : ndarray, shape (N, 3)
        Rotation vector  t = tan(θ/2)·b  or  (θ/2)·b.

    Returns
    -------
    u_plus : ndarray, shape (N, 3)
    """
    t2 = np.sum(t_vec**2, axis=1, keepdims=True)    # |t|², shape (N, 1)
    u_prime = u + np.cross(u, t_vec)                # Eq. (8)
    return u + (2.0 / (1.0 + t2)) * np.cross(u_prime, t_vec)   # Eq. (9)


# ── Abstract relativistic Boris base ──────────────────────────────────────────

class _BorisBase(Integrator):
    """Shared position-update logic for all Boris variants.

    Subclasses implement ``_rotate(u, B, gamma)``.
    """

    def __init__(self, state, field, dt, t_max, store_dt=None, relativistic=True):
        super().__init__(state, field, dt, t_max, store_dt=store_dt, relativistic=relativistic)

    # ── Physics helpers ───────────────────────────────────────────────────────

    @staticmethod
    def _gamma_from_u(u: np.ndarray) -> np.ndarray:
        """Lorentz factor from relativistic momentum per unit mass u = γv.

        γ = sqrt(1 + |u|²/c²),  shape (N, 1) for broadcasting.
        """
        u2 = np.sum(u**2, axis=1, keepdims=True)    # (N, 1)
        return np.sqrt(1.0 + u2 / C**2)

    # ── Rotation (to be overridden) ───────────────────────────────────────────

    def _rotate(self, u: np.ndarray, B: np.ndarray, gamma: np.ndarray) -> np.ndarray:
        """Magnetic rotation step.  Must return u_plus, shape (N, 3)."""
        raise NotImplementedError

    # ── Main step ─────────────────────────────────────────────────────────────

    def step(self):
        """One Boris leapfrog step: B-rotate → position."""
        r0 = self.state.r
        v0 = self.state.v
        t0 = self.state.t

        # Evaluate magnetic field at current position
        B = self.field(r0, t0)   # (N, 3)

        # u = γv : relativistic momentum per unit mass
        gamma0 = self.state.gamma()[:, None]   # (N, 1)
        u = gamma0 * v0                        # (N, 3)

        # Magnetic rotation (variant-specific)
        u_plus = self._rotate(u, B, gamma0)

        # New Lorentz factor and velocity
        gamma_new = self._gamma_from_u(u_plus)       # (N, 1)
        v_new = u_plus / gamma_new                   # (N, 3)

        # Position update (mid-point velocity)
        r_new = r0 + self.dt * v_new

        self.state.r = r_new
        self.state.v = v_new
        self.state.t = t0 + self.dt


# ── Boris-A ───────────────────────────────────────────────────────────────────

class BorisA(_BorisBase):
    """Boris-A solver: exact-angle t-vector rotation.

    Rotation vector (Eq. 7a of Zenitani & Umeda 2018):

        θ = qΔt|B| / (mγ)
        t = tan(θ/2) · b̂

    where b̂ = B / |B| is the unit vector along B.

    This is Boris's original (1970) procedure.  The exact ``tan(θ/2)``
    avoids the phase error of Boris-B at large θ (strong fields / large dt),
    but is more expensive due to the ``atan2`` and ``tan`` calls.

    Parameters
    ----------
    state : ParticleState
    field : FieldModel
    dt : float
        Fixed time step [s].  (Boris family is fixed-step only.)
    t_max : float
    store_dt : float or None
    relativistic : bool
        Default True.  Set False to treat γ=1 throughout (non-rel limit).
    """

    def __init__(self, state, field, dt, t_max, store_dt=None, relativistic=True):
        super().__init__(state, field, dt, t_max, store_dt=store_dt, relativistic=relativistic)

    def _rotate(self, u: np.ndarray, B: np.ndarray, gamma: np.ndarray) -> np.ndarray:
        """Exact-angle Boris rotation (Boris-A)."""
        B_mag = np.linalg.norm(B, axis=1, keepdims=True)    # |B|, (N, 1)
        # Safe unit vector: where |B|≈0 rotation is trivial
        safe_mag = np.where(B_mag > 0, B_mag, 1.0)
        b_hat = B / safe_mag                                  # (N, 3)

        # Rotation angle θ = qΔt|B| / (mγ)
        q_dt_over_m = (self.state.q[:, None] * self.dt
                       / self.state.m[:, None])               # (N, 1)
        theta = q_dt_over_m * B_mag / gamma                   # (N, 1)

        # t = tan(θ/2) · b̂
        t_vec = np.tan(theta / 2.0) * b_hat                   # (N, 3)

        # Zero rotation where |B| ≈ 0
        t_vec = np.where(B_mag > 0, t_vec, 0.0)

        return _boris_cross_rotate(u, t_vec)


# ── Boris-B ───────────────────────────────────────────────────────────────────

class BorisB(_BorisBase):
    """Boris-B solver: small-angle (textbook) t-vector rotation.

    Rotation vector (Eq. 7b of Zenitani & Umeda 2018):

        t = (θ/2) · b̂ = qΔt B / (2mγ)

    where the tan approximation tan(θ/2) ≈ θ/2 is used.

    This is the procedure described in Birdsall & Langdon and most PIC
    textbooks.  It introduces a small phase error of order θ² for large θ
    but is cheaper than Boris-A (no trigonometric calls).

    Equivalent to the non-relativistic Verlet algorithm when γ = 1.
    The original scalar formula in ``NEREL_PUTANJE_VERLE.py`` is the
    component-wise expansion of Eqs. (8-9) with t = qΔtB/(2m).

    Parameters
    ----------
    (same as BorisA)
    """

    def __init__(self, state, field, dt, t_max, store_dt=None, relativistic=True):
        super().__init__(state, field, dt, t_max, store_dt=store_dt, relativistic=relativistic)

    def _rotate(self, u: np.ndarray, B: np.ndarray, gamma: np.ndarray) -> np.ndarray:
        """Small-angle Boris rotation (Boris-B, textbook form)."""
        q_dt_over_2m = (self.state.q[:, None] * self.dt
                        / (2.0 * self.state.m[:, None]))      # (N, 1)
        # t = (qΔt / 2mγ) B  — equivalent to (θ/2) b̂ since |B|·b̂ = B
        t_vec = q_dt_over_2m / gamma * B                      # (N, 3)
        return _boris_cross_rotate(u, t_vec)


# ── Boris-C ───────────────────────────────────────────────────────────────────

class BorisC(_BorisBase):
    """Boris-C solver: exact analytic rotation (Zenitani & Umeda 2018).

    Rotation procedure (Eqs. 11-12):

        u‖ = (u · b̂) b̂                                   Eq. (11)
        u⁺ = u‖ + (u − u‖)cosθ + (u × b̂)sinθ           Eq. (12)

    where θ = qΔt|B|/(mγ) and b̂ = B/|B|.

    Unlike Boris-A/B this does NOT use the cross-product rotation kernel
    (Eqs. 8-9).  It directly decomposes u into components parallel and
    perpendicular to B, then applies an exact rotation by θ in the
    perpendicular plane.  This is second-order accurate (Strang splitting),
    time-reversible, and avoids any approximation in the rotation angle.

    A floor on |B|² (``B_FLOOR`` from constants) prevents division by zero
    in field-free regions as recommended by Zenitani & Umeda.

    Parameters
    ----------
    (same as BorisA)
    b_floor : float
        Override the default B² floor (default = ``constants.B_FLOOR``).
    """

    def __init__(self, state, field, dt, t_max, store_dt=None, relativistic=True,
                 b_floor: float = B_FLOOR):
        super().__init__(state, field, dt, t_max, store_dt=store_dt, relativistic=relativistic)
        self._b_floor = float(b_floor)

    def _rotate(self, u: np.ndarray, B: np.ndarray, gamma: np.ndarray) -> np.ndarray:
        """Exact analytic Boris-C rotation (Eqs. 11-12, Zenitani & Umeda 2018)."""
        # |B|² with floor to avoid division by zero (Zenitani & Umeda recommendation)
        B_mag2 = np.maximum(np.sum(B**2, axis=1, keepdims=True), self._b_floor)  # (N, 1)
        B_mag  = np.sqrt(B_mag2)                                                  # (N, 1)
        b_hat  = B / B_mag                                                        # (N, 3)

        # Rotation angle θ = qΔt|B| / (mγ)
        q_dt_over_m = (self.state.q[:, None] * self.dt
                       / self.state.m[:, None])                                   # (N, 1)
        theta = q_dt_over_m * B_mag / gamma                                       # (N, 1)

        # u‖ = (u · b̂) b̂  — parallel component  (Eq. 11)
        u_dot_b = np.sum(u * b_hat, axis=1, keepdims=True)                        # (N, 1)
        u_par   = u_dot_b * b_hat                                                 # (N, 3)

        # u − u‖ = perpendicular component
        u_perp  = u - u_par                                                       # (N, 3)

        # u × b̂  — used for sinθ term (note: cross(u_par, b̂) = 0)
        u_cross_b = np.cross(u, b_hat)                                            # (N, 3)

        # Eq. (12): exact rotation in the perpendicular plane
        u_plus = u_par + u_perp * np.cos(theta) + u_cross_b * np.sin(theta)      # (N, 3)

        return u_plus
