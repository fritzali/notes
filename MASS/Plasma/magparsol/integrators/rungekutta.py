"""
Runge–Kutta integrators, nonrelativistic (RKnonrel) and relativistic (RKrel).

Both use the Dormand–Prince 5(4) pair (Dormand & Prince 1980, J. Comput.
Appl. Math. 6, 19) with a fixed or an adaptive step; they differ only in
the acceleration:

* RKnonrel: dv/dt = (q/m) (E + v × B)
* RKrel:    dv/dt = (q / γm) (E + v × B - (v·E) v / c²)

With ``adaptive=True`` a step is accepted if the scaled error norm (root
mean square over all components of all particles) is at most 1, and the
next step is h · clip(S err^(-1/5), 0.1, 5) with the safety factor S,
limited to [dt_min, dt_max].
"""

import numpy as np
from magparsol.integrators.base import Integrator
from magparsol.constants import C

# Dormand–Prince Butcher tableau: nodes c, coefficients a, 5th order
# weights b5 and embedded 4th order weights b4. The last row of a equals
# b5, so the seventh stage is the derivative at the new point.
_C = np.array([0.0, 1/5, 3/10, 4/5, 8/9, 1.0, 1.0])

_A = [
    [],
    [1/5],
    [3/40,       9/40],
    [44/45,     -56/15,      32/9],
    [19372/6561, -25360/2187, 64448/6561, -212/729],
    [9017/3168,  -355/33,    46732/5247,   49/176,  -5103/18656],
    [35/384,      0.0,       500/1113,    125/192, -2187/6784,  11/84],
]

_B5 = np.array([35/384, 0.0, 500/1113, 125/192, -2187/6784, 11/84, 0.0])
_B4 = np.array([5179/57600, 0.0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40])
_E = _B5 - _B4   # error weights


class _RKBase(Integrator):
    """Dormand–Prince stepper; subclasses implement ``_derivatives``.

    Parameters
    ----------
    state, field, dt, t_max, store_dt, relativistic
        See :class:`~magparsol.integrators.base.Integrator`. With
        ``adaptive=True``, ``dt`` is the first trial step.
    adaptive : bool
        Adapt the step to the error tolerances.
    dt_min, dt_max : float or None
        Step limits [s], by default 1e-6 dt and 1e3 dt.
    rtol, atol : float
        Relative and absolute error tolerances, applied to every component
        of r [m] and v [m/s].
    safety : float
        Safety factor S of the step update.
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

    def _derivatives(self, r: np.ndarray, v: np.ndarray, t: float):
        """Return (dr/dt, dv/dt), each of shape (N, 3)."""
        raise NotImplementedError

    def _rk45_step(self, r0, v0, t0, h):
        """One Dormand–Prince step of size h from (r0, v0) at t0.

        Returns
        -------
        r, v : ndarray, shape (N, 3)
            5th order solution at t0 + h.
        er, ev : ndarray, shape (N, 3)
            Error estimate (difference to the embedded 4th order solution).
        """
        dr1, dv1 = self._derivatives(r0, v0, t0)

        r2 = r0 + h * _A[1][0] * dr1
        v2 = v0 + h * _A[1][0] * dv1
        dr2, dv2 = self._derivatives(r2, v2, t0 + _C[1]*h)

        r3 = r0 + h * (_A[2][0]*dr1 + _A[2][1]*dr2)
        v3 = v0 + h * (_A[2][0]*dv1 + _A[2][1]*dv2)
        dr3, dv3 = self._derivatives(r3, v3, t0 + _C[2]*h)

        r4 = r0 + h * (_A[3][0]*dr1 + _A[3][1]*dr2 + _A[3][2]*dr3)
        v4 = v0 + h * (_A[3][0]*dv1 + _A[3][1]*dv2 + _A[3][2]*dv3)
        dr4, dv4 = self._derivatives(r4, v4, t0 + _C[3]*h)

        r5 = r0 + h * (_A[4][0]*dr1 + _A[4][1]*dr2 + _A[4][2]*dr3 + _A[4][3]*dr4)
        v5 = v0 + h * (_A[4][0]*dv1 + _A[4][1]*dv2 + _A[4][2]*dv3 + _A[4][3]*dv4)
        dr5, dv5 = self._derivatives(r5, v5, t0 + _C[4]*h)

        r6 = r0 + h * (_A[5][0]*dr1 + _A[5][1]*dr2 + _A[5][2]*dr3 + _A[5][3]*dr4 + _A[5][4]*dr5)
        v6 = v0 + h * (_A[5][0]*dv1 + _A[5][1]*dv2 + _A[5][2]*dv3 + _A[5][3]*dv4 + _A[5][4]*dv5)
        dr6, dv6 = self._derivatives(r6, v6, t0 + _C[5]*h)

        r_out = r0 + h * (_B5[0]*dr1 + _B5[2]*dr3 + _B5[3]*dr4 + _B5[4]*dr5 + _B5[5]*dr6)
        v_out = v0 + h * (_B5[0]*dv1 + _B5[2]*dv3 + _B5[3]*dv4 + _B5[4]*dv5 + _B5[5]*dv6)

        # Seventh stage at the new point, needed for the error estimate only
        # (it is not reused as the first stage of the next step)
        dr7, dv7 = self._derivatives(r_out, v_out, t0 + h)

        er = h * (_E[0]*dr1 + _E[2]*dr3 + _E[3]*dr4 + _E[4]*dr5 + _E[5]*dr6 + _E[6]*dr7)
        ev = h * (_E[0]*dv1 + _E[2]*dv3 + _E[3]*dv4 + _E[4]*dv5 + _E[5]*dv6 + _E[6]*dv7)

        return r_out, v_out, er, ev

    def _error_norm(self, r0, v0, r5, v5, er, ev) -> float:
        """RMS of the error scaled by atol + rtol max(|y_old|, |y_new|)
        (Hairer, Nørsett & Wanner, Eq. II.4.11)."""
        sc_r = self.atol + self.rtol * np.maximum(np.abs(r0), np.abs(r5))
        sc_v = self.atol + self.rtol * np.maximum(np.abs(v0), np.abs(v5))
        n = er.size + ev.size
        return float(np.sqrt((np.sum((er/sc_r)**2) + np.sum((ev/sc_v)**2)) / n))

    def step(self):
        r0 = self.state.r.copy()
        v0 = self.state.v.copy()
        t0 = self.state.t

        if not self.adaptive:
            r5, v5, _, _ = self._rk45_step(r0, v0, t0, self.dt)
            self.state.r = r5
            self.state.v = v5
            self.state.t = t0 + self.dt
            return

        h = min(self.dt, self.t_max - t0)   # do not step past t_max
        while True:
            r5, v5, er, ev = self._rk45_step(r0, v0, t0, h)
            err = self._error_norm(r0, v0, r5, v5, er, ev)

            if err == 0.0:
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
                h = h_new
                if h <= self.dt_min:
                    # The tolerance cannot be met: take one step of dt_min
                    # anyway and warn.
                    import warnings
                    h = self.dt_min
                    r5, v5, er, ev = self._rk45_step(r0, v0, t0, h)
                    err = self._error_norm(r0, v0, r5, v5, er, ev)
                    warnings.warn(
                        f"Adaptive RK: step size hit minimum {h:.3e} s. "
                        f"Accepting step with err={err:.3e}.",
                        RuntimeWarning,
                        stacklevel=3,
                    )
                    self.state.r = r5
                    self.state.v = v5
                    self.state.t = t0 + h
                    self.dt = h
                    return


class RKnonrel(_RKBase):
    """Nonrelativistic Runge–Kutta integrator, dv/dt = (q/m)(E + v × B).

    Parameters are those of :class:`_RKBase`; ``relativistic`` defaults to
    False here (it only selects the energy formula of the plot helpers).
    """

    def __init__(self, state, field, dt, t_max, **kwargs):
        kwargs.setdefault("relativistic", False)
        super().__init__(state, field, dt, t_max, **kwargs)

    def _derivatives(self, r, v, t):
        B, E = self.field(r, t)
        q_m = self.state.q[:, None] / self.state.m[:, None]
        dvdt = q_m * (E + np.cross(v, B))
        return v.copy(), dvdt


class RKrel(_RKBase):
    """Relativistic Runge–Kutta integrator.

    From d(γmv)/dt = q(E + v × B) and dγ/dt = q (v·E) / (m c²),

        dv/dt = (q / γm) (E + v × B - (v·E) v / c²).

    Parameters are those of :class:`_RKBase`; ``relativistic`` defaults
    to True.
    """

    def __init__(self, state, field, dt, t_max, **kwargs):
        kwargs.setdefault("relativistic", True)
        super().__init__(state, field, dt, t_max, **kwargs)

    def _derivatives(self, r, v, t):
        B, E = self.field(r, t)

        speed2 = np.sum(v**2, axis=1, keepdims=True)
        gamma_inv = np.sqrt(np.clip(1.0 - speed2 / C**2, 1e-30, 1.0))

        vdotE = np.sum(v * E, axis=1, keepdims=True)

        q_m = self.state.q[:, None] / self.state.m[:, None]

        dvdt = q_m * gamma_inv * (
            E + np.cross(v, B) - (vdotE / C**2) * v
        )
        return v.copy(), dvdt
