"""
magparsol/integrators/base.py
------------------------------
Abstract base class for all particle integrators.

The template-method pattern is used throughout:

- ``run()`` owns the time-stepping loop, history recording, and optional
  live plotting — it does not know about physics.
- ``step()`` (abstract) advances ``self.state`` by one ``self.dt``.  For
  adaptive RK integrators the step size used may differ from ``self.dt``
  (which is updated in-place by the adaptive controller).

Subclasses only need to implement ``step()``.
"""

import numpy as np
from abc import ABC, abstractmethod
from magparsol.particles import ParticleState
from magparsol.fields import FieldModel
from magparsol.diagnostics import TrajectoryHistory


class Integrator(ABC):
    """Abstract particle integrator.

    Parameters
    ----------
    state : ParticleState
        Initial particle state.  **Modified in-place** during ``run()``.
        Pass ``state.copy()`` if you need to preserve the original.
    field : FieldModel
        Electromagnetic field callable.
    dt : float
        Time step [s].  For adaptive integrators this is the *initial*
        step size; it is updated after each accepted step.
    t_max : float
        End time of the simulation [s].
    store_dt : float or None
        Minimum elapsed simulation time between stored trajectory snapshots.
        ``None`` → store every step (fine for fixed-step integrators; may
        produce very large arrays for adaptive integrators with small dt).
    relativistic : bool
        Whether to use relativistic kinetic energy in diagnostics/plotting.
    """

    def __init__(
        self,
        state: ParticleState,
        field: FieldModel,
        dt: float,
        t_max: float,
        store_dt=None,
        relativistic: bool = False,
    ):
        self.state = state
        self.field = field
        self.dt = float(dt)
        self.t_max = float(t_max)
        self.store_dt = store_dt
        self.relativistic = relativistic

    # ── Abstract interface ────────────────────────────────────────────────────

    @abstractmethod
    def step(self):
        """Advance ``self.state`` by one time step.

        Must update ``self.state.r``, ``self.state.v``, and ``self.state.t``
        in place.  For adaptive integrators, ``self.dt`` may also be updated.
        """

    # ── Main loop ─────────────────────────────────────────────────────────────

    def run(
        self,
        live_plotter=None,
        live_every: int = 100,
        progress_every: int = 0,
    ) -> TrajectoryHistory:
        """Run the simulation from current state until ``t_max``.

        Parameters
        ----------
        live_plotter : LivePlotter or None
            If provided, ``live_plotter.update(state)`` is called every
            ``live_every`` steps for real-time display.
        live_every : int
            Call ``live_plotter.update`` every this many steps.
        progress_every : int
            Print ``t`` every this many steps.  0 = silent.

        Returns
        -------
        history : TrajectoryHistory (finalized)
        """
        history = TrajectoryHistory(store_dt=self.store_dt)
        history.record(self.state, force=True)   # always store t=0

        step_count = 0
        while self.state.t < self.t_max - 0.5 * self.dt:
            self.step()
            history.record(self.state)
            step_count += 1

            if progress_every and step_count % progress_every == 0:
                print(f"  t = {self.state.t:.6g} s  (step {step_count})")

            if live_plotter is not None and step_count % live_every == 0:
                live_plotter.update(self.state, t_max=self.t_max)

        history.finalize()
        return history

    # ── Convenience diagnostics (delegating to diagnostics / plotting modules) ─

    def gyroperiod_estimate(self, B_mag: float) -> float:
        """Estimate the gyro-period at the current particle state.

        Parameters
        ----------
        B_mag : float
            Representative magnetic field magnitude [T].

        Returns
        -------
        T_c : float [s]
        """
        from magparsol.diagnostics import gyroperiod
        q = float(self.state.q[0])
        m = float(self.state.m[0])
        return gyroperiod(q, m, B_mag)

    def suggest_dt(self, B_mag: float, steps_per_gyration: float = 100.0) -> float:
        """Suggest a dt for the given field strength."""
        from magparsol.diagnostics import suggest_dt
        q = float(self.state.q[0])
        m = float(self.state.m[0])
        return suggest_dt(q, m, B_mag, steps_per_gyration)

    def check_dt(self, B_mag: float, warn_threshold: float = 0.1) -> float:
        """Check current dt resolution against the gyro-period."""
        from magparsol.diagnostics import check_dt_resolution
        q = float(self.state.q[0])
        m = float(self.state.m[0])
        return check_dt_resolution(self.dt, q, m, B_mag, warn_threshold)

    # ── Plotting convenience methods ──────────────────────────────────────────

    def plot_trajectory_3d(self, history: TrajectoryHistory, **kwargs):
        from magparsol.plotting import plot_trajectory_3d
        return plot_trajectory_3d(history, **kwargs)

    def plot_trajectory_2d(self, history: TrajectoryHistory, **kwargs):
        from magparsol.plotting import plot_trajectory_2d
        return plot_trajectory_2d(history, **kwargs)

    def plot_energy(self, history: TrajectoryHistory, **kwargs):
        from magparsol.plotting import plot_energy
        return plot_energy(history, self.state.m, relativistic=self.relativistic, **kwargs)

    def plot_speed(self, history: TrajectoryHistory, **kwargs):
        from magparsol.plotting import plot_speed
        return plot_speed(history, relativistic=self.relativistic, **kwargs)

    # ── Representation ────────────────────────────────────────────────────────

    def __repr__(self):
        cls = type(self).__name__
        return (
            f"{cls}(N={self.state.N}, dt={self.dt:.3g}, "
            f"t_max={self.t_max:.3g}, relativistic={self.relativistic})"
        )
