"""
magparsol/plotting.py
----------------------
Plotting utilities for particle trajectory visualization and diagnostics.

All functions accept a finalized ``TrajectoryHistory`` and produce matplotlib
figures.  Particle index ``pid`` selects which particle to plot for multi-
particle states.

Functions
---------
plot_trajectory_3d    — 3-D trajectory with optional Earth sphere
plot_trajectory_2d    — Side-by-side X-Y and X-Z projection subplots
plot_energy           — Kinetic energy and relative error over time
plot_speed            — Speed (and γ for relativistic) over time
LivePlotter           — Class for real-time animated 2-D display during runs
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D   # noqa: F401 (registers 3d projection)
from magparsol.constants import R_EARTH, Q_E
from magparsol.diagnostics import TrajectoryHistory, relative_energy_error
from magparsol import style


# ── 3-D trajectory ────────────────────────────────────────────────────────────

def plot_trajectory_3d(
    history: TrajectoryHistory,
    pid: int = 0,
    length_unit: float = R_EARTH,
    unit_label: str = r"$R_E$",
    ax_lim: float = 70.0,
    earth_sphere: bool = False,
    color: str = "black",
    title: str = "Trajectory",
    ax=None,
    z_ratio: float = 1.0,
):
    """Plot a 3-D particle trajectory.

    Parameters
    ----------
    history : TrajectoryHistory (finalized)
    pid : int
        Particle index to plot.
    length_unit : float
        Divisor for converting metres to display units.
    unit_label : str
        Axis label unit string.
    ax_lim : float
        Symmetric axis limits in display units.
    earth_sphere : bool
        If True, draw a wireframe unit sphere representing Earth
        (only meaningful when length_unit = R_EARTH).
    color : str
        Trajectory line colour.
    title : str
        Figure title.
    ax : Axes3D or None
        Existing axes to draw into.  If None, a new figure is created.
    z_ratio : float
        Vertical extent relative to the horizontal one; the box aspect is
        adjusted so that the scale stays equal (the Earth stays round).

    Returns
    -------
    fig, ax
    """
    x = history.r[:, pid, 0] / length_unit
    y = history.r[:, pid, 1] / length_unit
    z = history.r[:, pid, 2] / length_unit

    if ax is None:
        fig = plt.figure(figsize=(7.6, 6.8))
        ax = fig.add_subplot(111, projection="3d")
    else:
        fig = ax.get_figure()

    if earth_sphere:
        u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:50j]
        xs = np.cos(u) * np.sin(v)
        ys = np.sin(u) * np.sin(v)
        zs = np.cos(v)
        ax.plot_surface(xs, ys, zs, color=style.EARTH_FACE, alpha=0.4,
                        linewidth=0, shade=True)

    ax.plot(x, y, z, color=color, linewidth=0.6, alpha=0.85,
            rasterized=style.raster(len(x)))
    ax.set_xlabel(f"$x$ [{unit_label}]")
    ax.set_ylabel(f"$y$ [{unit_label}]")
    ax.set_zlabel(f"$z$ [{unit_label}]")

    # Auto ax_lim from data when not specified
    if ax_lim == "auto" or ax_lim is None:
        ax_lim = float(np.max(np.abs([x, y, z]))) * 1.15
        ax_lim = max(ax_lim, 1.1)  # never zero

    ax.set_xlim3d(-ax_lim, ax_lim)
    ax.set_ylim3d(-ax_lim, ax_lim)
    ax.set_zlim3d(-ax_lim * z_ratio, ax_lim * z_ratio)
    # Equal scale on all axes so the Earth renders as a sphere
    ax.set_box_aspect((1, 1, z_ratio))
    ax.set_title(title)
    style.style_3d(ax)
    plt.tight_layout()
    return fig, ax


# ── 2-D projections ───────────────────────────────────────────────────────────

def plot_trajectory_2d(
    history: TrajectoryHistory,
    pid: int = 0,
    length_unit: float = 1.0,
    unit_label: str = "m",
    ax_lim: float = None,
    shared_limits: bool = False,
    title: str = "Trajectory",
):
    """Plot X-Y and X-Z projections side by side.

    Parameters
    ----------
    history : TrajectoryHistory (finalized)
    pid : int
        Particle index to plot.
    length_unit : float
        Divisor for converting metres to display units.
    unit_label : str
        Axis label unit string.
    ax_lim : float or None
        Explicit symmetric limit applied to both subplots.
        None → each subplot is fitted to its own data with equal aspect.
    shared_limits : bool
        If True, both subplots share the same numeric limits (useful when
        comparing absolute orbit extent between planes).  Ignored when ax_lim
        is set explicitly.
    title : str

    Returns
    -------
    fig, (ax1, ax2)
    """
    x = history.r[:, pid, 0] / length_unit
    y = history.r[:, pid, 1] / length_unit
    z = history.r[:, pid, 2] / length_unit

    # Panel widths follow the data aspect ratios so equal-scale panels fill
    # the figure instead of leaving gaps
    asp = [max(np.ptp(x), 1e-12) / max(np.ptp(c), 1e-12) for c in (y, z)]
    asp = [min(max(a, 0.2), 5.0) for a in asp]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(min(3.9 * sum(asp) + 1.3, 11), 4.25),
                                   gridspec_kw={"width_ratios": asp})
    r = style.raster(len(x))

    ax1.plot(x, y, color="steelblue", lw=0.8, rasterized=r)
    ax1.set_xlabel(f"$x$ [{unit_label}]")
    ax1.set_ylabel(f"$y$ [{unit_label}]")
    ax1.set_title("$xy$ Plane")
    ax1.set_aspect("equal")

    ax2.plot(x, z, color="firebrick", lw=0.8, rasterized=r)
    ax2.set_xlabel(f"$x$ [{unit_label}]")
    ax2.set_ylabel(f"$z$ [{unit_label}]")
    ax2.set_title("$xz$ Plane")
    ax2.set_aspect("equal")

    if ax_lim is not None:
        # Explicit user value → apply to both
        for ax in (ax1, ax2):
            ax.set_xlim(-ax_lim, ax_lim)
            ax.set_ylim(-ax_lim, ax_lim)
    elif shared_limits:
        # Same auto-computed limit for both (informative for flat orbits)
        lim = _auto_lim(x, np.concatenate([y, z]))
        for ax in (ax1, ax2):
            ax.set_xlim(-lim, lim)
            ax.set_ylim(-lim, lim)
    else:
        # Fit each panel to its data (centred, equal aspect, 8% margin)
        for ax, a, b in ((ax1, x, y), (ax2, x, z)):
            _fit_box(ax, a, b)

    fig.suptitle(title)
    plt.tight_layout()
    return fig, (ax1, ax2)


def _auto_lim(a, b, margin=1.15):
    """Symmetric axis limit covering both arrays."""
    m = margin * max(float(np.abs(a).max()), float(np.abs(b).max()))
    return max(m, 1e-10)


def _fit_box(ax, a, b, margin=0.08):
    """Data limits with a margin; equal scale, the box shrinks to fit."""
    pad = margin * max(np.ptp(a), np.ptp(b), 1e-12)
    ax.set_xlim(a.min() - pad, a.max() + pad)
    ax.set_ylim(b.min() - pad, b.max() + pad)
    ax.set_aspect("equal", adjustable="box")


# ── Energy diagnostics ────────────────────────────────────────────────────────

def plot_energy(
    history: TrajectoryHistory,
    m: np.ndarray,
    pid: int = 0,
    relativistic: bool = True,
    title: str = "Energy",
):
    """Plot kinetic energy and relative energy error over time.

    Parameters
    ----------
    history : TrajectoryHistory (finalized)
    m : ndarray, shape (N,)
        Particle masses [kg].
    pid : int
        Particle index.
    relativistic : bool
        Whether to compute relativistic kinetic energy.
    title : str

    Returns
    -------
    fig, (ax1, ax2)
    """
    K = history.kinetic_energy(m, relativistic=relativistic)[:, pid]   # (S,)
    err = relative_energy_error(history, m, relativistic=relativistic)[:, pid]  # (S,)
    t = history.t

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.6, 5.1), sharex=True)

    K_eV = K / Q_E
    k_max = float(np.max(np.abs(K_eV))) if len(K_eV) else 0.0
    exp3, prefix = next(((e, p) for e, p in ((9, "G"), (6, "M"), (3, "k"))
                         if k_max >= 10**e), (0, ""))
    ax1.plot(t, K_eV / 10**exp3, color="steelblue")
    ax1.set_ylabel(f"$K$ [{prefix}eV]")
    ax1.set_title(title)

    ax2.plot(t, err, color="firebrick")
    ax2.set_ylabel(r"$\Delta K\,/\,K_0$")
    ax2.set_xlabel("$t$ [s]")
    ax2.axhline(0, color="black", linewidth=0.6, linestyle="--")

    plt.tight_layout()
    return fig, (ax1, ax2)


# ── Speed / gamma ─────────────────────────────────────────────────────────────

def plot_speed(
    history: TrajectoryHistory,
    pid: int = 0,
    relativistic: bool = True,
    title: str = "Speed",
):
    """Plot particle speed (and γ for relativistic) over time.

    Returns
    -------
    fig, axes
    """
    from magparsol.constants import C
    speed = history.speed()[:, pid]
    t = history.t

    if relativistic:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.6, 5.1), sharex=True)
        gamma = history.gamma()[:, pid]
        ax1.plot(t, speed / C, color="steelblue")
        ax1.set_ylabel(r"$|\mathbf{v}|\,/\,c$")
        ax1.set_title(title)

        ax2.plot(t, gamma, color="goldenrod")
        ax2.set_ylabel(r"$\gamma$")
        ax2.set_xlabel("$t$ [s]")

        plt.tight_layout()
        return fig, (ax1, ax2)
    else:
        fig, ax = plt.subplots(figsize=(7.6, 3.0))
        ax.plot(t, speed, color="steelblue")
        ax.set_ylabel(r"$|\mathbf{v}|$ [m/s]")
        ax.set_xlabel("$t$ [s]")
        ax.set_title(title)
        plt.tight_layout()
        return fig, ax


# ── Live (animated) plotter ───────────────────────────────────────────────────

def _in_notebook() -> bool:
    """Return True when running inside a Jupyter kernel."""
    try:
        from IPython import get_ipython
        cfg = get_ipython()
        return cfg is not None and "IPKernelApp" in cfg.config
    except ImportError:
        return False


class LivePlotter:
    """Real-time 2-D trajectory display during an integrator run.

    Mimics the ``IZLAZ.init_plots`` / ``write_plots`` pattern from the
    original code, decoupled from the integrator itself.

    Parameters
    ----------
    ax_lim : float
        Symmetric axis limits for both subplots.
    length_unit : float
        Divisor for converting metres to display units.
    unit_label : str
        Axis label unit string.
    pause : float
        Pause duration per update [s] (passed to ``plt.pause``).
    pid : int
        Particle index to display.
    """

    def __init__(
        self,
        ax_lim: float = 70.0,
        length_unit: float = R_EARTH,
        unit_label: str = r"$R_E$",
        pause: float = 1e-5,
        pid: int = 0,
    ):
        self.length_unit    = length_unit
        self.unit_label     = unit_label
        self.pause          = pause
        self.pid            = pid
        self._ax_lim        = ax_lim
        self._notebook_mode = _in_notebook()
        self._init_figure()

    def _init_figure(self):
        self.fig, (self.ax1, self.ax2) = plt.subplots(1, 2, figsize=(11, 5.1))
        for ax in (self.ax1, self.ax2):
            ax.set_aspect("equal")
            ax.set_xlim(-self._ax_lim, self._ax_lim)
            ax.set_ylim(-self._ax_lim, self._ax_lim)
        self.ax1.set_xlabel(f"$x$ [{self.unit_label}]")
        self.ax1.set_ylabel(f"$y$ [{self.unit_label}]")
        self.ax2.set_xlabel(f"$x$ [{self.unit_label}]")
        self.ax2.set_ylabel(f"$z$ [{self.unit_label}]")
        self.ax1.set_title("$xy$ Plane")
        self.ax2.set_title("$xz$ Plane")
        plt.ion()
        plt.tight_layout()
        plt.show()

    def update(self, state, t_max: float = None):
        """Draw the current particle position.

        Parameters
        ----------
        state : ParticleState
        t_max : float or None
            If provided, show progress in the subplot title.
        """
        d = self.length_unit
        pid = self.pid
        x = float(state.r[pid, 0]) / d
        y = float(state.r[pid, 1]) / d
        z = float(state.r[pid, 2]) / d

        if t_max is not None:
            self.fig.suptitle(f"$t$ = {style.format_time_pair(state.t, t_max)}")

        self.ax1.plot(x, y, ".", color="steelblue", markersize=2)
        self.ax2.plot(x, z, ".", color="firebrick", markersize=2)
        if self._notebook_mode:
            from IPython.display import display, clear_output
            clear_output(wait=True)
            display(self.fig)
        else:
            plt.pause(self.pause)

    def close(self):
        plt.ioff()
        plt.close(self.fig)
