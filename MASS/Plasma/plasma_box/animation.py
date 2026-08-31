"""
plasma_box/animation.py
-----------------------
Overview panel (static multi-quantity summary figure) and GIF animation system.

Overview panel layout (2×3 gridspec)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
┌──────────────┬──────────────┬─────────────────┐
│  X-Y position│  X-Z position│  Field lines     │
├──────────────┼──────────────┼─────────────────┤
│  X-Y velocity│  X-Z velocity│  Spectrum        │
└──────────────┴──────────────┴─────────────────┘

Panel registry
~~~~~~~~~~~~~~
Each panel is a small object with build(ax) → artists and
update(artists, frame_idx) → artists.  This shared structure powers
both standalone per-panel GIFs and the combined overview GIF from one
code path.

GIF generation
~~~~~~~~~~~~~~
make_panel_gif  — single panel standalone GIF
make_overview_gif — all (or selected) panels in the overview layout

Adaptive-RK note: non-uniform time histories are resampled before FFT
inside radiation.spectrum_fft; no special handling needed here.
"""

import warnings
import numpy as np
import matplotlib
matplotlib.use("Agg")   # safe default; caller can switch before import
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.animation import FuncAnimation, PillowWriter

from plasma_box.diagnostics import TrajectoryHistory, relative_energy_error
from plasma_box.constants import C, R_EARTH


# ── Helpers ───────────────────────────────────────────────────────────────────

def _auto_lim(a, b, margin=1.15):
    """Symmetric axis limit covering both arrays with a margin."""
    m = margin * max(float(np.abs(a).max()), float(np.abs(b).max()))
    return m if m > 0 else 1.0


def _trajectory_bbox(history):
    """Combined bounding box across all particles."""
    r = history.r.reshape(-1, 3)
    return r.min(axis=0), r.max(axis=0)


# ── Per-panel classes ─────────────────────────────────────────────────────────

class _PositionPanel:
    """X-Y or X-Z trajectory projection."""

    def __init__(self, plane="xy", length_unit=1.0, unit_label="m",
                 normalize_v=False):
        self.plane       = plane
        self.length_unit = length_unit
        self.unit_label  = unit_label

    def build(self, ax, history, **kwargs):
        lu = self.length_unit
        r  = history.r
        N  = r.shape[1]
        lines, dots = [], []
        colors = plt.cm.tab10(np.linspace(0, 1, N))
        for pid in range(N):
            x = r[:, pid, 0] / lu
            c2 = r[:, pid, 1] / lu if self.plane == "xy" else r[:, pid, 2] / lu
            ln, = ax.plot([], [], lw=0.8, alpha=0.7, color=colors[pid])
            pt, = ax.plot([], [], 'o', ms=4, color=colors[pid])
            lines.append(ln); dots.append(pt)
        # axis limits from full trajectory
        x_all  = r[:, :, 0].ravel() / lu
        c2_all = (r[:, :, 1] if self.plane == "xy" else r[:, :, 2]).ravel() / lu
        lim    = _auto_lim(x_all, c2_all)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_xlabel(f"x [{self.unit_label}]")
        ax.set_ylabel(f"{'y' if self.plane=='xy' else 'z'} [{self.unit_label}]")
        ax.set_title(f"X-{'Y' if self.plane=='xy' else 'Z'} position")
        ax.grid(True, alpha=0.3)
        self._history = history
        return {"lines": lines, "dots": dots}

    def update(self, artists, i):
        lu = self.length_unit
        r  = self._history.r
        for pid, (ln, pt) in enumerate(zip(artists["lines"], artists["dots"])):
            x  = r[:i+1, pid, 0] / lu
            c2 = (r[:i+1, pid, 1] if self.plane == "xy" else r[:i+1, pid, 2]) / lu
            ln.set_data(x, c2)
            pt.set_data([x[-1]], [c2[-1]])
        return list(artists["lines"]) + list(artists["dots"])


class _VelocityPanel:
    """X-Y or X-Z velocity hodograph."""

    def __init__(self, plane="xy", normalize_v=False):
        self.plane       = plane
        self.normalize_v = normalize_v

    def build(self, ax, history, **kwargs):
        N      = history.v.shape[1]
        scale  = C if self.normalize_v else 1.0
        label  = "v/c" if self.normalize_v else "m/s"
        colors = plt.cm.tab10(np.linspace(0, 1, N))
        lines  = []
        for pid in range(N):
            vx = history.v[:, pid, 0] / scale
            c2 = (history.v[:, pid, 1] if self.plane == "xy"
                  else history.v[:, pid, 2]) / scale
            ln, = ax.plot([], [], lw=0.6, alpha=0.5, color=colors[pid])
            lines.append(ln)
        vx_all  = history.v[:, :, 0].ravel() / scale
        c2_all  = (history.v[:, :, 1] if self.plane == "xy"
                   else history.v[:, :, 2]).ravel() / scale
        lim     = _auto_lim(vx_all, c2_all)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_xlabel(f"v_x [{label}]")
        ax.set_ylabel(f"{'v_y' if self.plane=='xy' else 'v_z'} [{label}]")
        ax.set_title(f"V hodograph X-{'Y' if self.plane=='xy' else 'Z'}")
        ax.grid(True, alpha=0.3)
        self._history = history
        self._scale   = scale
        return {"lines": lines}

    def update(self, artists, i):
        s = self._scale
        for pid, ln in enumerate(artists["lines"]):
            vx = self._history.v[:i+1, pid, 0] / s
            c2 = (self._history.v[:i+1, pid, 1] if self.plane == "xy"
                  else self._history.v[:i+1, pid, 2]) / s
            ln.set_data(vx, c2)
        return list(artists["lines"])


class _FieldPanel:
    """Field line / arrow panel (static or periodically refreshed)."""

    def __init__(self, field, history, length_unit=1.0, unit_label="m",
                 components=("B", "E"), density="low",
                 field_update_every=10):
        self.field              = field
        self.history            = history
        self.length_unit        = length_unit
        self.unit_label         = unit_label
        self.components         = components
        self.density            = density
        self.field_update_every = field_update_every
        self._t_series          = None   # time series for animation

    def build(self, ax, history=None, t_series=None, **kwargs):
        from plasma_box.fieldlines import plot_field_lines
        h = history or self.history
        self._t_series = t_series
        # Draw initial field lines at t=0
        plot_field_lines(
            self.field, history=h,
            components=self.components,
            density=self.density,
            t=0.0,
            length_unit=self.length_unit,
            unit_label=self.unit_label,
            projection="xy",
            ax=ax,
        )
        ax.set_title("Field config (B/E)")
        self._ax = ax
        self._h  = h
        return {}   # artists managed internally by re-plot

    def update(self, artists, i, k=None, t_current=0.0):
        """Re-draw field lines if field is time-dependent and on cadence."""
        if self.field.is_static:
            return []
        if k is not None and k % self.field_update_every != 0:
            return []
        from plasma_box.fieldlines import plot_field_lines
        self._ax.cla()
        plot_field_lines(
            self.field, history=self._h,
            components=self.components,
            density=self.density,
            t=t_current,
            length_unit=self.length_unit,
            unit_label=self.unit_label,
            projection="xy",
            ax=self._ax,
        )
        self._ax.set_title(f"Field config (t={t_current:.3g}s)")
        return []


class _SpectrumPanel:
    """Spectrum panel: animated FFT or static retarded integral."""

    def __init__(self, q, m, method="fft", spectrum_update_every=8,
                 show_individual=False, observer=None,
                 store_dt_warn_period=None):
        self.q                     = q
        self.m                     = m
        self.method                = method
        self.spectrum_update_every = spectrum_update_every
        self.show_individual       = show_individual
        self.observer              = observer
        self.store_dt_warn_period  = store_dt_warn_period
        self._prev_spec            = None
        self._converged            = False
        self._ref_freqs            = None

    def build(self, ax, history, **kwargs):
        from plasma_box.radiation import (spectrum_fft, spectrum_retarded,
                                          ensemble_spectrum, _check_spectrum_convergence)
        self._history = history
        self._ax      = ax
        ax.set_xlabel("Frequency [Hz]")
        ax.set_ylabel("Power [arb.]")
        ax.set_title("Spectrum")
        ax.grid(True, alpha=0.3)

        if self.method == "retarded":
            # Compute once from full trajectory
            N = history.r.shape[1]
            if N == 1:
                f, p = spectrum_retarded(history, pid=0, observer=self.observer)
            else:
                f, p, _ = ensemble_spectrum(history, method="retarded",
                                            observer=self.observer)
            self._ref_freqs = f
            ax.semilogy(f, p + 1e-40, color="steelblue", lw=1.2, label="Total")
            ax.legend(fontsize=7)

        # Compute full reference spectrum for convergence checking
        N = history.r.shape[1]
        if N == 1:
            f_full, p_full = spectrum_fft(
                history, pid=0, observer=self.observer,
                store_dt_warn_period=self.store_dt_warn_period)
        else:
            f_full, p_full, _ = ensemble_spectrum(
                history, method="fft", observer=self.observer)
        self._ref_freqs  = f_full
        self._full_spec  = p_full
        self._line_total, = ax.semilogy(f_full, p_full + 1e-40,
                                         color="steelblue", lw=1.2, label="Total",
                                         alpha=0.0 if self.method == "fft" else 1.0)
        return {"line_total": self._line_total}

    def update(self, artists, i, k=None):
        from plasma_box.radiation import (spectrum_fft, ensemble_spectrum,
                                          _check_spectrum_convergence)
        if self.method == "retarded":
            return []   # static — already drawn in build()

        if self._converged:
            return []
        if k is not None and k % self.spectrum_update_every != 0:
            return []

        N = self._history.r.shape[1]
        if N == 1:
            f, p = spectrum_fft(self._history, pid=0, upto=i,
                                 observer=self.observer,
                                 store_dt_warn_period=self.store_dt_warn_period)
        else:
            f, p, _ = ensemble_spectrum(self._history, method="fft",
                                         observer=self.observer,
                                         upto=i)

        # Interpolate onto fixed reference frequency grid
        if self._ref_freqs is not None and len(f) > 1:
            p_interp = np.interp(self._ref_freqs, f, p, left=0.0, right=0.0)
        else:
            p_interp = p

        # Check convergence
        converged = _check_spectrum_convergence(
            self._prev_spec, p_interp, self._ref_freqs)
        if converged:
            self._converged = True
        self._prev_spec = p_interp.copy() if p_interp is not None else None

        line = artists.get("line_total")
        if line is not None and len(p_interp) > 0:
            line.set_alpha(0.9)
            line.set_ydata(p_interp + 1e-40)
            self._ax.relim()
            self._ax.autoscale_view(scaley=True)

        return [line] if line else []


# ── LAYOUT definition ─────────────────────────────────────────────────────────
# Maps panel name → (row, col) in the 2×3 gridspec
LAYOUT = {
    "position_xy": (0, 0),
    "position_xz": (0, 1),
    "field":        (0, 2),
    "velocity_xy":  (1, 0),
    "velocity_xz":  (1, 1),
    "spectrum":     (1, 2),
}


def _build_panel_registry(history, field, q, m,
                           length_unit=1.0, unit_label="m",
                           normalize_v=False,
                           field_components=("B", "E"),
                           field_density="low",
                           spectrum_method="fft",
                           spectrum_update_every=8,
                           show_individual=False,
                           field_update_every=10,
                           store_dt_warn_period=None) -> dict:
    """Instantiate one panel object per panel name."""
    return {
        "position_xy": _PositionPanel("xy", length_unit, unit_label),
        "position_xz": _PositionPanel("xz", length_unit, unit_label),
        "velocity_xy": _VelocityPanel("xy", normalize_v),
        "velocity_xz": _VelocityPanel("xz", normalize_v),
        "field":       _FieldPanel(field, history, length_unit, unit_label,
                                    field_components, field_density,
                                    field_update_every),
        "spectrum":    _SpectrumPanel(q, m, spectrum_method,
                                       spectrum_update_every, show_individual,
                                       store_dt_warn_period=store_dt_warn_period),
    }


# ── Static overview figure ────────────────────────────────────────────────────

def plot_overview(
    history: TrajectoryHistory,
    field,
    q: np.ndarray,
    m: np.ndarray,
    length_unit: float = 1.0,
    unit_label: str = "m",
    normalize_v: bool = False,
    field_components=("B", "E"),
    field_density: str = "low",
    spectrum_method: str = "fft",
    show_individual: bool = False,
    title: str = "Simulation overview",
    store_dt_warn_period: float = None,
):
    """Static 2×3 overview panel.

    Parameters
    ----------
    history : TrajectoryHistory (finalized)
    field : FieldModel
    q, m : ndarray, shape (N,)
    length_unit : float
    unit_label : str
    normalize_v : bool   — show velocity in units of c
    field_components : tuple of "B" and/or "E"
    field_density : "low" | "medium" | "high"
    spectrum_method : "fft" | "retarded"
    show_individual : bool   — overlay per-particle spectra (N>1)
    store_dt_warn_period : float or None  — gyroperiod for Nyquist warning

    Returns
    -------
    fig
    """
    fig = plt.figure(figsize=(16, 9))
    gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.35)
    fig.suptitle(title, fontsize=13)

    registry = _build_panel_registry(
        history, field, q, m, length_unit, unit_label, normalize_v,
        field_components, field_density, spectrum_method,
        store_dt_warn_period=store_dt_warn_period,
    )

    S = len(history.t)
    for name, (row, col) in LAYOUT.items():
        ax  = fig.add_subplot(gs[row, col])
        panel = registry[name]
        artists = panel.build(ax, history)
        i = S - 1
        if name == "field":
            panel.update(artists, i, k=0, t_current=float(history.t[i]))
        elif name == "spectrum":
            panel.update(artists, i, k=0)
        else:
            panel.update(artists, i)

    return fig


# ── GIF generation ────────────────────────────────────────────────────────────

def make_panel_gif(
    history: TrajectoryHistory,
    panel_name: str,
    field=None,
    q: np.ndarray = None,
    m: np.ndarray = None,
    filename: str = None,
    fps: int = 15,
    n_frames: int = 150,
    writer: str = "gif",
    length_unit: float = 1.0,
    unit_label: str = "m",
    normalize_v: bool = False,
    field_components=("B", "E"),
    field_density: str = "low",
    spectrum_method: str = "fft",
    spectrum_update_every: int = 8,
    field_update_every: int = 10,
    store_dt_warn_period: float = None,
    **panel_kwargs,
):
    """Generate a GIF for a single panel.

    Parameters
    ----------
    panel_name : one of "position_xy", "position_xz", "velocity_xy",
                 "velocity_xz", "field", "spectrum"
    filename : output path.  Defaults to f"{panel_name}.gif"
    writer : "gif" (Pillow).  Future: "mp4" (FFMpeg).
    (other params: same as make_overview_gif)
    """
    if panel_name not in LAYOUT:
        raise ValueError(f"Unknown panel '{panel_name}'. Choose from: {list(LAYOUT)}")

    registry = _build_panel_registry(
        history, field, q, m, length_unit, unit_label, normalize_v,
        field_components, field_density, spectrum_method,
        spectrum_update_every, field_update_every=field_update_every,
        store_dt_warn_period=store_dt_warn_period,
    )
    panel = registry[panel_name]

    fig, ax = plt.subplots(figsize=(7, 6))
    artists  = panel.build(ax, history)
    frame_indices = np.round(
        np.linspace(0, len(history.t) - 1, n_frames)).astype(int)

    def update(k):
        i          = frame_indices[k]
        t_current  = float(history.t[i])
        if panel_name == "field":
            return panel.update(artists, i, k=k, t_current=t_current)
        elif panel_name == "spectrum":
            return panel.update(artists, i, k=k)
        else:
            return panel.update(artists, i)

    anim = FuncAnimation(fig, update, frames=len(frame_indices), blit=False)
    outfile = filename or f"{panel_name}.gif"
    if writer == "gif":
        anim.save(outfile, writer=PillowWriter(fps=fps))
    # Future: elif writer == "mp4": anim.save(outfile, writer=FFMpegWriter(fps=fps))
    else:
        raise ValueError(f"Unknown writer '{writer}'. Use 'gif'.")
    plt.close(fig)
    print(f"Saved: {outfile}")
    return outfile


def make_overview_gif(
    history: TrajectoryHistory,
    field,
    q: np.ndarray,
    m: np.ndarray,
    filename: str = "overview.gif",
    fps: int = 15,
    n_frames: int = 150,
    writer: str = "gif",
    panels=None,
    length_unit: float = 1.0,
    unit_label: str = "m",
    normalize_v: bool = False,
    field_components=("B", "E"),
    field_density: str = "low",
    spectrum_method: str = "fft",
    spectrum_update_every: int = 8,
    field_update_every: int = 10,
    store_dt_warn_period: float = None,
):
    """Generate a combined overview GIF with up to 6 panels.

    Parameters
    ----------
    panels : list of str or None
        Subset of panel names to include.  None → all 6.
    filename : output path
    writer : "gif".  Future: "mp4".
    (other params same as make_panel_gif)

    Returns
    -------
    filename : str
    """
    panels = panels or list(LAYOUT.keys())

    registry = _build_panel_registry(
        history, field, q, m, length_unit, unit_label, normalize_v,
        field_components, field_density, spectrum_method,
        spectrum_update_every, field_update_every=field_update_every,
        store_dt_warn_period=store_dt_warn_period,
    )

    fig = plt.figure(figsize=(16, 9))
    gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.35)
    axes_map    = {}
    artists_map = {}

    for name in panels:
        row, col = LAYOUT[name]
        ax = fig.add_subplot(gs[row, col])
        axes_map[name]    = ax
        artists_map[name] = registry[name].build(ax, history)

    frame_indices = np.round(
        np.linspace(0, len(history.t) - 1, n_frames)).astype(int)

    def update(k):
        i         = frame_indices[k]
        t_current = float(history.t[i])
        updated   = []
        for name in panels:
            panel = registry[name]
            if name == "field":
                updated += panel.update(artists_map[name], i,
                                         k=k, t_current=t_current)
            elif name == "spectrum":
                updated += panel.update(artists_map[name], i, k=k)
            else:
                updated += panel.update(artists_map[name], i)
        return updated

    anim = FuncAnimation(fig, update, frames=len(frame_indices), blit=False)
    if writer == "gif":
        anim.save(filename, writer=PillowWriter(fps=fps))
    # Future: elif writer == "mp4": anim.save(filename, writer=FFMpegWriter(fps=fps))
    else:
        raise ValueError(f"Unknown writer '{writer}'. Use 'gif'.")
    plt.close(fig)
    print(f"Saved: {filename}")
    return filename
