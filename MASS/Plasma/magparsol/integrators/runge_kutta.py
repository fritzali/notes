"""
magparsol/integrators/runge_kutta.py
--------------------------------------
Runge-Kutta integrators for non-relativistic and relativistic particle motion.

Both integrators share the same Dormand-Prince RK45 stepper machinery
(fixed or adaptive step).  The only physics difference is the ``_derivatives``
method:

- ``RKNonrel``      — dv/dt = (q/m)(v × B), non-relativistic
- ``RKRelativistic`` — dv/dt = (q/m)γ⁻¹(v × B)

Dormand-Prince RK45 tableau
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Stages:    7 (FSAL — First Same As Last)
Orders:    5th-order solution; 4th-order embedded for error estimate
Reference: Dormand & Prince (1980), J. Comput. Appl. Math. 6(1), 19-26

Adaptive step-size control
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Uses PI controller (simplified to I-only: Hairer et al.):

    h_new = h * min(f_max, max(f_min, S * err^{-1/5}))

where err is the RMS scaled error norm, S = safety factor (0.9 by default).

Note on multi-particle (N > 1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
All array operations are written in (N, 3) shape so that N particles
advance simultaneously with a single integrator call.  The error norm
for the adaptive controller is taken across all particles.
"""

import numpy as np
from magparsol.integrators.base import Integrator
from magparsol.constants import C

# ── Dormand-Prince RK45 Butcher tableau ──────────────────────────────────────
# c nodes (fraction of dt)
_C = np.array([0.0, 1/5, 3/10, 4/5, 8/9, 1.0, 1.0])

# a coefficients (lower triangular)
_A = [
    [],                                                          # k1 (trivial)
    [1/5],                                                       # k2
    [3/40,       9/40],                                          # k3
    [44/45,     -56/15,      32/9],                              # k4
    [19372/6561, -25360/2187, 64448/6561, -212/729],             # k5
    [9017/3168,  -355/33,    46732/5247,   49/176,  -5103/18656],# k6
    [35/384,      0.0,       500/1113,    125/192, -2187/6784,  11/84],  # k7=k1_next (FSAL)
]

# b5 — 5th-order weights (same as last a row = a[6])
_B5 = np.array([35/384, 0.0, 500/1113, 125/192, -2187/6784, 11/84, 0.0])

# b4 — 4th-order weights (for error estimate e = b5 - b4)
_B4 = np.array([5179/57600, 0.0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40])

# Error coefficients e = b5 - b4
_E = _B5 - _B4


class _RKBase(Integrator):
    """Internal base providing the Dormand-Prince RK45 stepper.

    Subclasses implement ``_derivatives(r, v, t)`` returning (dr/dt, dv/dt)
    arrays of shape (N, 3).
    """

    def __init__(
        self,
        state,
        field,
        dt: float,
        t_max: float,
        adaptive: bool = False,
        dt_min: float = None,
        dt_max: float = None,
        rtol: float = 1e-6,
        atol: float = 1e-9,
        safety: float = 0.9,
        store_dt=None,
        relativistic: bool = False,
    ):
        super().__init__(state, field, dt, t_max, store_dt=store_dt, relativistic=relativistic)
        self.adaptive = adaptive
        self.dt_min = float(dt_min) if dt_min is not None else dt * 1e-6
        self.dt_max = float(dt_max) if dt_max is not None else dt * 1e3
        self.rtol = float(rtol)
        self.atol = float(atol)
        self.safety = float(safety)

    # ── Physics (to be implemented by subclasses) ─────────────────────────────

    def _derivatives(self, r: np.ndarray, v: np.ndarray, t: float):
        """Return (drdt, dvdt), each shape (N, 3)."""
        raise NotImplementedError

    # ── Dormand-Prince RK45 core ──────────────────────────────────────────────

    def _rk45_step(self, r0, v0, t0, h):
        """One RK45 step from (r0, v0) at time t0 with step h.

        Returns
        -------
        r5, v5 : ndarray, shape (N, 3)  — 5th-order solution
        r4, v4 : ndarray, shape (N, 3)  — 4th-order solution (for error)
        """
        # Stage 1
        dr1, dv1 = self._derivatives(r0, v0, t0)

        # Stage 2
        r2 = r0 + h * _A[1][0] * dr1
        v2 = v0 + h * _A[1][0] * dv1
        dr2, dv2 = self._derivatives(r2, v2, t0 + _C[1]*h)

        # Stage 3
        r3 = r0 + h * (_A[2][0]*dr1 + _A[2][1]*dr2)
        v3 = v0 + h * (_A[2][0]*dv1 + _A[2][1]*dv2)
        dr3, dv3 = self._derivatives(r3, v3, t0 + _C[2]*h)

        # Stage 4
        r4 = r0 + h * (_A[3][0]*dr1 + _A[3][1]*dr2 + _A[3][2]*dr3)
        v4 = v0 + h * (_A[3][0]*dv1 + _A[3][1]*dv2 + _A[3][2]*dv3)
        dr4, dv4 = self._derivatives(r4, v4, t0 + _C[3]*h)

        # Stage 5
        r5 = r0 + h * (_A[4][0]*dr1 + _A[4][1]*dr2 + _A[4][2]*dr3 + _A[4][3]*dr4)
        v5 = v0 + h * (_A[4][0]*dv1 + _A[4][1]*dv2 + _A[4][2]*dv3 + _A[4][3]*dv4)
        dr5, dv5 = self._derivatives(r5, v5, t0 + _C[4]*h)

        # Stage 6
        r6 = r0 + h * (_A[5][0]*dr1 + _A[5][1]*dr2 + _A[5][2]*dr3 + _A[5][3]*dr4 + _A[5][4]*dr5)
        v6 = v0 + h * (_A[5][0]*dv1 + _A[5][1]*dv2 + _A[5][2]*dv3 + _A[5][3]*dv4 + _A[5][4]*dv5)
        dr6, dv6 = self._derivatives(r6, v6, t0 + _C[5]*h)

        # 5th-order solution (stage 7 = FSAL, reused as next k1)
        r_out = r0 + h * (_B5[0]*dr1 + _B5[2]*dr3 + _B5[3]*dr4 + _B5[4]*dr5 + _B5[5]*dr6)
        v_out = v0 + h * (_B5[0]*dv1 + _B5[2]*dv3 + _B5[3]*dv4 + _B5[4]*dv5 + _B5[5]*dv6)

        # 4th-order solution for error estimate
        r4_out = r0 + h * (_B4[0]*dr1 + _B4[2]*dr3 + _B4[3]*dr4 + _B4[4]*dr5 + _B4[5]*dr6 + _B4[6]*r_out)  # noqa
        v4_out = v0 + h * (_B4[0]*dv1 + _B4[2]*dv3 + _B4[3]*dv4 + _B4[4]*dv5 + _B4[5]*dv6 + _B4[6]*v_out)  # noqa

        # Correct 4th order: recompute properly without the recursive term
        # (b4[6] multiplies the 7th stage which IS r_out / v_out itself;
        #  simplest and correct: use e = b5 - b4 directly on the stage derivatives)
        dr7, dv7 = self._derivatives(r_out, v_out, t0 + h)

        # Error vector = h * sum_i e_i * k_i
        er = h * (_E[0]*dr1 + _E[2]*dr3 + _E[3]*dr4 + _E[4]*dr5 + _E[5]*dr6 + _E[6]*dr7)
        ev = h * (_E[0]*dv1 + _E[2]*dv3 + _E[3]*dv4 + _E[4]*dv5 + _E[5]*dv6 + _E[6]*dv7)

        return r_out, v_out, er, ev

    # ── Scalar error norm ──────────────────────────────────────────────────────

    def _error_norm(self, r0, v0, r5, v5, er, ev) -> float:
        """RMS scaled error norm (Hairer et al. Eq. 4.11).

        Scale = atol + rtol * max(|y|, |y_new|)
        """
        sc_r = self.atol + self.rtol * np.maximum(np.abs(r0), np.abs(r5))
        sc_v = self.atol + self.rtol * np.maximum(np.abs(v0), np.abs(v5))
        n = er.size + ev.size
        return float(np.sqrt((np.sum((er/sc_r)**2) + np.sum((ev/sc_v)**2)) / n))

    # ── step() ────────────────────────────────────────────────────────────────

    def step(self):
        r0 = self.state.r.copy()
        v0 = self.state.v.copy()
        t0 = self.state.t

        if not self.adaptive:
            # Fixed-step RK45 (use 5th-order solution; ignore error estimate)
            r5, v5, _, _ = self._rk45_step(r0, v0, t0, self.dt)
            self.state.r = r5
            self.state.v = v5
            self.state.t = t0 + self.dt
            return

        # Adaptive step-size control
        # Clamp step so we don't overshoot t_max
        h = min(self.dt, self.t_max - t0)
        while True:
            r5, v5, er, ev = self._rk45_step(r0, v0, t0, h)
            err = self._error_norm(r0, v0, r5, v5, er, ev)

            if err == 0.0:
                # Perfect step — jump to dt_max
                h_new = self.dt_max
                accepted = True
            else:
                factor = self.safety * (1.0 / err) ** 0.2
                factor = max(0.1, min(5.0, factor))
                h_new = np.clip(h * factor, self.dt_min, self.dt_max)
                accepted = err <= 1.0

            if accepted:
                self.state.r = r5
                self.state.v = v5
                self.state.t = t0 + h
                self.dt = h_new
                return
            else:
                # Reject step, try smaller h
                h = h_new
                if h < self.dt_min:
                    # Cannot satisfy tolerance — accept anyway with a warning
                    import warnings
                    warnings.warn(
                        f"Adaptive RK: step size {h:.3e} hit minimum {self.dt_min:.3e}. "
                        "Accepting step with err={err:.3e}.",
                        RuntimeWarning,
                        stacklevel=3,
                    )
                    self.state.r = r5
                    self.state.v = v5
                    self.state.t = t0 + h
                    self.dt = self.dt_min
                    return


# ── Non-relativistic RK integrator ───────────────────────────────────────────

class RKNonrel(_RKBase):
    """Non-relativistic Runge-Kutta integrator (RK45, fixed or adaptive step).

    Equation of motion::

        dr/dt = v
        dv/dt = (q/m)(v × B)

    Consolidates the ``DKP_mdipol.py`` RK4 approach and generalises
    the field model and step-size strategy.

    Parameters
    ----------
    state : ParticleState
    field : FieldModel
    dt : float
        Time step (initial step if adaptive=True) [s].
    t_max : float
    adaptive : bool
        If True, use adaptive Dormand-Prince step control.
    dt_min, dt_max : float or None
        Bounds on adaptive step size [s].
    rtol, atol : float
        Relative and absolute tolerances for adaptive control.
    safety : float
        Safety factor for step-size update (default 0.9).
    store_dt : float or None
        Sub-sampling interval for trajectory storage [s].
    """

    def __init__(self, state, field, dt, t_max, **kwargs):
        kwargs.setdefault("relativistic", False)
        super().__init__(state, field, dt, t_max, **kwargs)

    def _derivatives(self, r, v, t):
        """dr/dt = v,  dv/dt = (q/m)(v × B)."""
        B = self.field(r, t)   # (N, 3)
        q_m = self.state.q[:, None] / self.state.m[:, None]   # (N, 1) broadcasts
        dvdt = q_m * np.cross(v, B)
        return v.copy(), dvdt


# ── Relativistic RK integrator ────────────────────────────────────────────────

class RKRelativistic(_RKBase):
    """Relativistic Runge-Kutta integrator (RK45, fixed or adaptive step).

    Equation of motion (4-force projected onto 3-velocity):

        dr/dt = v
        dv/dt = (q/m) γ⁻¹ [ v × B ]

    This is the standard relativistic EOM in terms of velocity (not
    4-momentum), derived from d(γmv)/dt = q(v × B).

    Parameters
    ----------
    (same as RKNonrel)
    """

    def __init__(self, state, field, dt, t_max, **kwargs):
        kwargs.setdefault("relativistic", True)
        super().__init__(state, field, dt, t_max, **kwargs)

    def _derivatives(self, r, v, t):
        """Relativistic Lorentz force equation."""
        B = self.field(r, t)   # (N, 3)

        speed2 = np.sum(v**2, axis=1, keepdims=True)        # (N, 1)
        gamma_inv = np.sqrt(np.clip(1.0 - speed2 / C**2, 1e-30, 1.0))  # (N, 1)

        q_m = self.state.q[:, None] / self.state.m[:, None]   # (N, 1)

        dvdt = q_m * gamma_inv * np.cross(v, B)
        return v.copy(), dvdt
