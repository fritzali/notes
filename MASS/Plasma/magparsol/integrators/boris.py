"""
Relativistic Boris integrators A, B and C, after Zenitani & Umeda (2018),
"On the Boris solver in particle-in-cell simulation", Phys. Plasmas 25, 112110.

All three advance u = γv via the same leapfrog sequence:

1. half electric push      u⁻ = uⁿ + (q Δt / 2m) E                  (Eq. 3)
2. magnetic rotation       u⁺ = R(u⁻), by the angle θ = q Δt |B| / (m γ⁻)
3. half electric push      uⁿ⁺¹ = u⁺ + (q Δt / 2m) E                (Eq. 5)
4. position update         rⁿ⁺¹ = rⁿ + Δt uⁿ⁺¹ / γⁿ⁺¹

They differ only in step 2:

* Boris A: t = tan(θ/2) b̂, rotation by the exact angle θ.
* Boris B: t = (θ/2) b̂, the usual textbook form. It rotates by
  2 arctan(θ/2) < θ, i.e. it lags in phase by about θ³/12 per step.
* Boris C: explicit rotation with cos θ and sin θ (Eqs. 11 and 12).

The rotation preserves |u|, so in a pure magnetic field none of them
changes the energy, whatever the step size. The push is always
relativistic; for v ≪ c it reduces to the classical Boris push.
"""

import numpy as np
from magparsol.integrators.base import Integrator
from magparsol.constants import C, B_FLOOR


def _boris_cross_rotate(u_minus: np.ndarray, t_vec: np.ndarray) -> np.ndarray:
    """Boris rotation through the vector t (Eqs. 8 and 9 of Zenitani & Umeda).

    u' = u⁻ + u⁻ × t,  u⁺ = u⁻ + 2/(1 + t²) u' × t. This rotates u⁻ about t
    by the angle 2 arctan |t|. Arrays have shape (N, 3).
    """
    t2 = np.sum(t_vec**2, axis=1, keepdims=True)
    u_prime = u_minus + np.cross(u_minus, t_vec)
    return u_minus + (2.0 / (1.0 + t2)) * np.cross(u_prime, t_vec)


class _BorisBase(Integrator):
    """Electric half pushes and position update shared by all variants.

    Subclasses implement ``_rotate(u_minus, B, gamma_minus)``.
    """

    def __init__(self, state, field, dt, t_max, store_dt=None, relativistic=True):
        super().__init__(state, field, dt, t_max, store_dt=store_dt, relativistic=relativistic)

    @staticmethod
    def _gamma_from_u(u: np.ndarray) -> np.ndarray:
        """γ = sqrt(1 + u²/c²) for u = γv, shape (N, 1) for broadcasting."""
        u2 = np.sum(u**2, axis=1, keepdims=True)
        return np.sqrt(1.0 + u2 / C**2)

    def _half_E_push(self, u: np.ndarray, E: np.ndarray) -> np.ndarray:
        """u + (q Δt / 2m) E."""
        q_over_2m = self.state.q[:, None] / (2.0 * self.state.m[:, None])
        return u + q_over_2m * self.dt * E

    def _rotate(self, u_minus: np.ndarray, B: np.ndarray, gamma_minus: np.ndarray) -> np.ndarray:
        """Magnetic rotation; returns u⁺ of shape (N, 3)."""
        raise NotImplementedError

    def step(self):
        r0 = self.state.r
        v0 = self.state.v
        t0 = self.state.t

        B, E = self.field(r0, t0)

        gamma0 = self.state.gamma()[:, None]
        u = gamma0 * v0

        u_minus = self._half_E_push(u, E)
        gamma_minus = self._gamma_from_u(u_minus)
        u_plus = self._rotate(u_minus, B, gamma_minus)
        u_new = self._half_E_push(u_plus, E)

        gamma_new = self._gamma_from_u(u_new)
        v_new = u_new / gamma_new
        r_new = r0 + self.dt * v_new

        self.state.r = r_new
        self.state.v = v_new
        self.state.t = t0 + self.dt


class BorisA(_BorisBase):
    """Boris A: the original Boris (1970) rotation with t = tan(θ/2) b̂.

    Rotates by exactly θ = q Δt |B| / (m γ⁻), so there is no phase error,
    at the cost of one tangent per step.

    Parameters
    ----------
    state : ParticleState
        Initial state; it is advanced in place.
    field : FieldModel
    dt : float
        Time step [s]. The Boris integrators use a fixed step.
    t_max : float
        End time [s].
    store_dt : float or None
        Sampling interval of the history, see :class:`TrajectoryHistory`.
    relativistic : bool
        Energy formula used by the plotting helpers. The push itself is
        always relativistic.
    """

    def __init__(self, state, field, dt, t_max, store_dt=None, relativistic=True):
        super().__init__(state, field, dt, t_max, store_dt=store_dt, relativistic=relativistic)

    def _rotate(self, u_minus: np.ndarray, B: np.ndarray, gamma_minus: np.ndarray) -> np.ndarray:
        B_mag = np.linalg.norm(B, axis=1, keepdims=True)
        safe_mag = np.where(B_mag > 0, B_mag, 1.0)
        b_hat = B / safe_mag

        q_dt_over_m = (self.state.q[:, None] * self.dt
                       / self.state.m[:, None])
        theta = q_dt_over_m * B_mag / gamma_minus

        t_vec = np.tan(theta / 2.0) * b_hat
        t_vec = np.where(B_mag > 0, t_vec, 0.0)   # no rotation where B = 0

        return _boris_cross_rotate(u_minus, t_vec)


class BorisB(_BorisBase):
    """Boris B: textbook rotation with t = (θ/2) b̂ = q Δt B / (2 m γ⁻).

    This is the form in Birdsall & Langdon and most PIC codes. It needs
    no trigonometric functions but rotates by 2 arctan(θ/2) instead of θ,
    a phase lag of about θ³/12 per step. Parameters as for :class:`BorisA`.
    """

    def __init__(self, state, field, dt, t_max, store_dt=None, relativistic=True):
        super().__init__(state, field, dt, t_max, store_dt=store_dt, relativistic=relativistic)

    def _rotate(self, u_minus: np.ndarray, B: np.ndarray, gamma_minus: np.ndarray) -> np.ndarray:
        q_dt_over_2m = (self.state.q[:, None] * self.dt
                        / (2.0 * self.state.m[:, None]))
        t_vec = q_dt_over_2m / gamma_minus * B
        return _boris_cross_rotate(u_minus, t_vec)


class BorisC(_BorisBase):
    """Boris C: explicit rotation by θ (Eqs. 11 and 12 of Zenitani & Umeda).

    u⁻ is split into its parts along and across b̂ = B/|B|,

        u‖ = (u⁻·b̂) b̂,   u⁺ = u‖ + (u⁻ - u‖) cos θ + (u⁻ × b̂) sin θ,

    which rotates by exactly θ like Boris A. |B|² is floored at ``b_floor``
    so that b̂ stays defined where the field vanishes.

    Parameters
    ----------
    b_floor : float
        Floor on |B|² [T²], by default ``constants.B_FLOOR``.

    The other parameters are as for :class:`BorisA`.
    """

    def __init__(self, state, field, dt, t_max, store_dt=None, relativistic=True,
                 b_floor: float = B_FLOOR):
        super().__init__(state, field, dt, t_max, store_dt=store_dt, relativistic=relativistic)
        self._b_floor = float(b_floor)

    def _rotate(self, u_minus: np.ndarray, B: np.ndarray, gamma_minus: np.ndarray) -> np.ndarray:
        B_mag2 = np.maximum(np.sum(B**2, axis=1, keepdims=True), self._b_floor)
        B_mag  = np.sqrt(B_mag2)
        b_hat  = B / B_mag

        q_dt_over_m = (self.state.q[:, None] * self.dt
                       / self.state.m[:, None])
        theta = q_dt_over_m * B_mag / gamma_minus

        u_dot_b = np.sum(u_minus * b_hat, axis=1, keepdims=True)
        u_par   = u_dot_b * b_hat
        u_perp  = u_minus - u_par
        u_cross_b = np.cross(u_minus, b_hat)   # equals u_perp × b̂

        u_plus = u_par + u_perp * np.cos(theta) + u_cross_b * np.sin(theta)

        return u_plus
