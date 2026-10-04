"""
Base class shared by all integrators.

:meth:`Integrator.run` owns the time loop, the history and the optional
live plot; subclasses only implement :meth:`Integrator.step`.
"""

import numpy as np
from abc import ABC, abstractmethod
from magparsol.particles import ParticleState
from magparsol.fields import FieldModel
from magparsol.diagnostics import TrajectoryHistory


class Integrator(ABC):
    """Advances a :class:`ParticleState` in a :class:`FieldModel`.

    Parameters
    ----------
    state : ParticleState
        Initial state. It is advanced in place, so pass ``state.copy()`` to
        keep the original.
    field : FieldModel
    dt : float
        Time step [s]; the first trial step for adaptive integrators, which
        update it as they go.
    t_max : float
        End time [s].
    store_dt : float or None
        Sampling interval of the stored history [s]; None stores every step.
    relativistic : bool
        Energy formula used by :meth:`plot_energy` and :meth:`plot_speed`.
        It does not change the equations of motion.
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

    @abstractmethod
    def step(self):
        """Advance ``self.state`` (r, v and t) by one step."""

    def run(
        self,
        live_plotter=None,
        live_every: int = 100,
        progress_every: int = 0,
    ) -> TrajectoryHistory:
        """Integrate from the current state to ``t_max``.

        Parameters
        ----------
        live_plotter : LivePlotter or None
            Updated with the current state every ``live_every`` steps.
        live_every : int
        progress_every : int
            Print the time every this many steps; 0 prints nothing.

        Returns
        -------
        TrajectoryHistory
            Finalized history, starting with the initial state.
        """
        history = TrajectoryHistory(store_dt=self.store_dt)
        history.record(self.state, force=True)

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

    # Shortcuts to the diagnostics, using the first particle's q and m.
    # The gyration quantities are nonrelativistic.

    def gyroperiod_estimate(self, B_mag: float) -> float:
        """Gyroperiod [s] in a field of strength ``B_mag`` [T]."""
        from magparsol.diagnostics import gyroperiod
        q = float(self.state.q[0])
        m = float(self.state.m[0])
        return gyroperiod(q, m, B_mag)

    def suggest_dt(self, B_mag: float, steps_per_gyration: float = 100.0) -> float:
        """Time step [s] giving ``steps_per_gyration`` steps per gyroperiod."""
        from magparsol.diagnostics import suggest_dt
        q = float(self.state.q[0])
        m = float(self.state.m[0])
        return suggest_dt(q, m, B_mag, steps_per_gyration)

    def check_dt(self, B_mag: float, warn_threshold: float = 0.1) -> float:
        """dt / T_c for the current step, with a warning above ``warn_threshold``."""
        from magparsol.diagnostics import check_dt_resolution
        q = float(self.state.q[0])
        m = float(self.state.m[0])
        return check_dt_resolution(self.dt, q, m, B_mag, warn_threshold)

    # Shortcuts to the plotting functions

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

    def __repr__(self):
        cls = type(self).__name__
        return (
            f"{cls}(N={self.state.N}, dt={self.dt:.3g}, "
            f"t_max={self.t_max:.3g}, relativistic={self.relativistic})"
        )
