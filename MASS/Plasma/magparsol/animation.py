"""
magparsol/animation.py
-----------------------
Overview panel (static multi-quantity summary figure) and GIF animation system.

Overview panel layout (2×3 gridspec)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
┌────────────────┬────────────────┬───────────────┐
│  X-Y Position  │  X-Z Position  │  Field + orbit│
├────────────────┼────────────────┼───────────────┤
│  X-Y Velocity  │  X-Z Velocity  │    Spectrum   │
└────────────────┴────────────────┴───────────────┘

Panel registry
~~~~~~~~~~~~~~
Each panel is a small object with ``build(ax, history)`` → artists and
``update(artists, i, k, t)`` → artists.  This shared structure powers the
static overview, standalone per-panel GIFs and the combined overview GIF
from one code path.

Spectrum animation
~~~~~~~~~~~~~~~~~~
The FFT panel shows the final spectrum as a faint grey reference and, on
top, the spectrum of the trajectory *up to the current frame*.  As the
observation window grows the cyclotron line narrows (Δf ≈ 1/T_window) and
rises out of the window side-lobes.  The retarded-time integral spectrum
is expensive and always shown as the final result only.

Adaptive-RK note: non-uniform time histories are resampled before the FFT
inside ``radiation.spectrum_fft``; no special handling is needed here.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from magparsol.diagnostics import TrajectoryHistory
from magparsol.constants import C, R_EARTH
from magparsol import style


# ── Helpers ───────────────────────────────────────────────────────────────────

def _auto_lim(a, b, margin=1.15):
    """Symmetric axis limit covering both arrays with a margin."""
    m = margin * max(float(np.abs(a).max()), float(np.abs(b).max()))
    return m if m > 0 else 1.0


def _is_dipole(field):
    from magparsol.fields import EarthDipole
    return isinstance(field, EarthDipole)


def _colors(N):
    return style.colors(N)


def mean_gyrofrequency(history, q, m, field):
    """Trajectory-averaged relativistic gyrofrequency f_c = |q|B / (2π γ m) [Hz].

    Used to normalise spectrum axes.  Returns None if no field is given.
    """
    if field is None or q is None or m is None:
        return None
    S, N, _ = history.r.shape
    idx = np.unique(np.linspace(0, S - 1, min(S, 400)).astype(int))
    f = []
    for s in idx:
        B, _ = field(history.r[s], float(history.t[s]))
        Bm = np.linalg.norm(B, axis=1)
        beta2 = np.clip(np.sum(history.v[s]**2, axis=1) / C**2, 0, 1 - 1e-15)
        gamma = 1.0 / np.sqrt(1.0 - beta2)
        f.append(np.abs(q) * Bm / (2*np.pi * gamma * m))
    f_c = float(np.mean(f))
    return f_c if f_c > 0 else None


# ── Per-panel classes ─────────────────────────────────────────────────────────

class _PositionPanel:
    """X-Y or X-Z trajectory projection with a growing trail."""

    def __init__(self, plane="xy", length_unit=1.0, unit_label="m", earth=False):
        self.plane       = plane
        self.length_unit = length_unit
        self.unit_label  = unit_label
        self.earth       = earth

    def _c2(self, r):
        return r[..., 1] if self.plane == "xy" else r[..., 2]

    def build(self, ax, history, **kwargs):
        lu = self.length_unit
        r  = history.r / lu
        N  = r.shape[1]
        lines, dots = [], []
        for pid, col in zip(range(N), _colors(N)):
            ln, = ax.plot([], [], lw=0.7, alpha=0.85, color=col,
                          rasterized=style.raster(len(r)))
            pt, = ax.plot([], [], "o", ms=4, color=col, mec="black", mew=0.5, zorder=5)
            lines.append(ln); dots.append(pt)
        lim = _auto_lim(r[..., 0], self._c2(r))
        if self.earth:
            from magparsol.fieldlines import draw_earth_2d
            draw_earth_2d(ax, lu, zorder=1)
            lim = max(lim, 1.3 * R_EARTH / lu)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_xlabel(f"$x$ [{self.unit_label}]")
        ax.set_ylabel(f"${self.plane[1]}$ [{self.unit_label}]")
        ax.set_title(f"Position (${self.plane}$)")
        self._r = r
        return {"lines": lines, "dots": dots}

    def update(self, artists, i, **kwargs):
        r = self._r
        for pid, (ln, pt) in enumerate(zip(artists["lines"], artists["dots"])):
            x, c2 = r[:i+1, pid, 0], self._c2(r[:i+1, pid])
            ln.set_data(x, c2)
            pt.set_data([x[-1]], [c2[-1]])
        return list(artists["lines"]) + list(artists["dots"])


class _VelocityPanel:
    """X-Y or X-Z velocity hodograph.

    The full hodograph is drawn faintly; the highlighted trail shows only the
    last ``trail`` samples (≈ a few gyrations) so long runs stay readable.
    """

    def __init__(self, plane="xy", normalize_v=False, trail=None):
        self.plane       = plane
        self.normalize_v = normalize_v
        self.trail       = trail

    def build(self, ax, history, **kwargs):
        scale = C if self.normalize_v else 1.0
        unit  = "$c$" if self.normalize_v else "m/s"
        v     = history.v / scale
        N     = v.shape[1]
        c2    = v[..., 1] if self.plane == "xy" else v[..., 2]
        lines, dots = [], []
        for pid, col in zip(range(N), _colors(N)):
            if self.trail is not None:
                ax.plot(v[:, pid, 0], c2[:, pid], lw=0.3, color="lightgrey", zorder=1,
                        rasterized=style.raster(len(v)))
            ln, = ax.plot([], [], lw=0.8, alpha=0.85, color=col)
            pt, = ax.plot([], [], "o", ms=4, color=col, mec="black", mew=0.5, zorder=5)
            lines.append(ln); dots.append(pt)
        lim = _auto_lim(v[..., 0], c2)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_xlabel(f"$v_x$ [{unit}]")
        ax.set_ylabel(f"$v_{self.plane[1]}$ [{unit}]")
        ax.set_title(f"Velocity (${self.plane}$)")
        self._v, self._c2 = v, c2
        return {"lines": lines, "dots": dots}

    def update(self, artists, i, **kwargs):
        i0 = 0 if self.trail is None else max(0, i + 1 - self.trail)
        for pid, (ln, pt) in enumerate(zip(artists["lines"], artists["dots"])):
            vx, c2 = self._v[i0:i+1, pid, 0], self._c2[i0:i+1, pid]
            ln.set_data(vx, c2)
            pt.set_data([vx[-1]], [c2[-1]])
        return list(artists["lines"]) + list(artists["dots"])


class _FieldPanel:
    """Field lines / arrow grid with the particle orbit overlaid.

    Spatially uniform, time-dependent fields (e.g. ``CyclotronWaveField``)
    are animated by updating the arrow vectors.  Non-uniform time-dependent
    fields are re-traced every ``field_update_every`` frames.
    """

    def __init__(self, field, length_unit=1.0, unit_label="m",
                 components=("B", "E"), density="low", field_update_every=10,
                 projection="auto", trail=None):
        self.field              = field
        self.length_unit        = length_unit
        self.unit_label         = unit_label
        self.components         = components
        self.density            = density
        self.field_update_every = field_update_every
        if projection == "auto":
            projection = "xz" if _is_dipole(field) else "xy"
        self.projection         = projection
        self.trail              = trail

    def _draw(self, ax, t):
        from magparsol.fieldlines import plot_field_lines
        comps = tuple(c for c in self.components
                      if not (_is_dipole(self.field) and c == "E"))
        plot_field_lines(
            self.field, components=comps or ("B",), density=self.density,
            t=t, length_unit=self.length_unit, unit_label=self.unit_label,
            ax_lim=self._lim, projection=self.projection, ax=ax,
            ref_times=self._t_all, title="Field and Orbit",
            L_shells=self._L_shells,
        )

    def build(self, ax, history, **kwargs):
        from magparsol.fieldlines import _plane
        e1, e2, _, _ = _plane(self.projection)
        r = history.r / self.length_unit
        self._p1, self._p2 = r @ e1, r @ e2
        self._lim = _auto_lim(self._p1, self._p2)
        self._L_shells = (2, 3, 4, 6, 8)
        if _is_dipole(self.field):
            self._lim = max(self._lim, 1.3 * R_EARTH / self.length_unit)
            # L-shells that fit inside the panel
            lim_RE = self._lim * self.length_unit / R_EARTH
            self._L_shells = tuple(L for L in (1.5, 2, 3, 4, 5, 6, 8, 10, 12)
                                   if L <= 0.95 * lim_RE) or (1.5,)
        self._t_all = history.t
        self._ax    = ax
        self._draw(ax, float(history.t[0]))
        from magparsol.fieldlines import FadingTrail
        N = r.shape[1]
        # Trail ≈ 1/10 of the run (at least 30 samples) so the field stays visible
        n_trail = self.trail or max(30, len(history.t) // 10)
        lines, dots = [], []
        for pid, col in zip(range(N), _colors(N)):
            lines.append(FadingTrail(ax, n_trail, color=style.ORBIT_COLOR,
                                     alpha=0.75, lw=0.9, zorder=6))
            pt, = ax.plot([], [], "o", ms=4.5, color=col, mec="black", mew=0.5, zorder=7)
            dots.append(pt)
        return {"lines": lines, "dots": dots}

    def update(self, artists, i, k=None, t_current=0.0, **kwargs):
        out = []
        for pid, (tr, pt) in enumerate(zip(artists["lines"], artists["dots"])):
            out.append(tr.set_data(self._p1[:i+1, pid], self._p2[:i+1, pid]))
            pt.set_data([self._p1[i, pid]], [self._p2[i, pid]])
            out.append(pt)
        if self.field.is_static:
            return out
        arrows = getattr(self._ax, "_mps_uniform_arrows", [])
        if arrows:
            for arr in arrows:
                out += arr.update(t_current)
        elif k is None or k % self.field_update_every == 0:
            self._ax.cla()
            self._draw(self._ax, t_current)
            for tr, pt in zip(artists["lines"], artists["dots"]):
                self._ax.add_collection(tr.lc); self._ax.add_line(pt)
        return out


class _SpectrumPanel:
    """Spectrum panel: FFT that builds up in time, or static retarded integral."""

    def __init__(self, q, m, field=None, method="fft", spectrum_update_every=5,
                 show_individual=False, observer=None,
                 store_dt_warn_period=None, f_norm=None, f_max_norm=6.0,
                 final_only=False):
        self.q                     = q
        self.m                     = m
        self.field                 = field
        self.method                = method
        self.spectrum_update_every = spectrum_update_every
        self.show_individual       = show_individual
        self.observer              = observer
        self.store_dt_warn_period  = store_dt_warn_period
        self.f_norm                = f_norm
        self.f_max_norm            = f_max_norm
        self.final_only            = final_only   # static figure: no live line

    def _spectrum(self, upto=None, method="fft"):
        from magparsol.radiation import spectrum_fft, spectrum_retarded, ensemble_spectrum
        h = self._history
        if h.r.shape[1] == 1:
            if method == "retarded":
                f, p = spectrum_retarded(h, pid=0, observer=self.observer)
            else:
                f, p = spectrum_fft(h, pid=0, upto=upto, observer=self.observer,
                                    store_dt_warn_period=self.store_dt_warn_period)
            return f, p, []
        kw = {} if method == "retarded" else {"upto": upto}
        return ensemble_spectrum(h, method=method, observer=self.observer, **kw)

    def build(self, ax, history, **kwargs):
        self._history = history
        self._ax      = ax
        f_c = self.f_norm or mean_gyrofrequency(history, self.q, self.m, self.field)
        self._fn = f_c or 1.0
        xlabel = r"$f\,/\,\langle f_c\rangle$" if f_c else r"$f$ [Hz]"

        f, p, ind = self._spectrum(method=self.method)
        x = f / self._fn
        if self.show_individual:
            for pi in ind:
                ax.semilogy(x, np.where(pi > 0, pi, np.nan), color="grey",
                            lw=0.5, alpha=0.3)
        static = self.method == "retarded" or self.final_only
        ax.semilogy(x, np.where(p > 0, p, np.nan),
                    color="steelblue" if static else style.REF_COLOR,
                    lw=1.2 if static else 0.9,
                    label="Retarded" if self.method == "retarded" else "Final")
        live, = ax.semilogy([], [], color="firebrick", lw=1.2, label="Current")
        info  = ax.text(0.97, 0.95, "", transform=ax.transAxes, ha="right",
                        va="top", fontsize=8,
                        bbox=dict(boxstyle="square", fc="white", ec="black",
                                  lw=0.6, alpha=0.92))

        # x range: up to where the final spectrum has fallen 6 decades below
        # its peak (with margin), capped at f_max_norm
        pos   = (x > 0) & (p > 0)
        p_top = float(p[pos].max()) if np.any(pos) else 1.0
        above = np.where(pos & (p > p_top * 1e-6))[0]
        x_max = x[above[-1]] * 1.15 if len(above) else x.max()
        if f_c:
            x_max = min(max(x_max, 1.5), self.f_max_norm)
        x_max = min(x_max, x.max())
        ax.set_xlim(0, x_max)
        # Two decades of headroom keep the legend and time box off the data
        ax.set_ylim(p_top * 1e-6, p_top * 300)
        if f_c:
            ax.axvline(1.0, color="black", lw=0.6, ls="--", alpha=0.6)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Power [arb. units]")
        ax.set_title("Spectrum" + (" (Retarded)" if self.method == "retarded" else ""))
        if not static:
            style.legend(ax, fontsize=8, loc="upper left")
        self._f_c = f_c
        return {"live": live, "info": info}

    def update(self, artists, i, k=None, **kwargs):
        live, info = artists["live"], artists["info"]
        t0, t_end = float(self._history.t[0]), float(self._history.t[-1])
        info.set_text(f"$t$ = {style.format_time(self._history.t[i] - t0, t_end - t0)}")
        if self.method == "retarded" or self.final_only:
            return [info]
        last = i == len(self._history.t) - 1
        if k is not None and k % self.spectrum_update_every != 0 and not last:
            return [info]
        if i < 8:
            live.set_data([], [])
            return [live, info]
        f, p, _ = self._spectrum(upto=i + 1)
        live.set_data(f / self._fn, np.where(p > 0, p, np.nan))
        return [live, info]


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
                          spectrum_update_every=5,
                          show_individual=False,
                          field_update_every=10,
                          store_dt_warn_period=None,
                          field_projection="auto",
                          field_trail=None,
                          final_spectrum_only=False) -> dict:
    """Instantiate one panel object per panel name."""
    earth = _is_dipole(field)
    # Hodograph trail ≈ 3 gyrations when the run covers many more than that
    trail = None
    f_c = mean_gyrofrequency(history, q, m, field)
    if f_c and len(history.t) > 1:
        n_gyr = f_c * float(history.t[-1] - history.t[0])
        if n_gyr > 10:
            trail = max(10, int(round(3 * len(history.t) / n_gyr)))
    return {
        "position_xy": _PositionPanel("xy", length_unit, unit_label, earth),
        "position_xz": _PositionPanel("xz", length_unit, unit_label, earth),
        "velocity_xy": _VelocityPanel("xy", normalize_v, trail),
        "velocity_xz": _VelocityPanel("xz", normalize_v, trail),
        "field":       _FieldPanel(field, length_unit, unit_label,
                                   field_components, field_density,
                                   field_update_every, field_projection,
                                   field_trail),
        "spectrum":    _SpectrumPanel(q, m, field, spectrum_method,
                                      spectrum_update_every, show_individual,
                                      store_dt_warn_period=store_dt_warn_period,
                                      final_only=final_spectrum_only),
    }


def _update_panel(panel, artists, i, k, history):
    return panel.update(artists, i, k=k, t_current=float(history.t[i]))


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
    title: str = "Overview",
    store_dt_warn_period: float = None,
    field_projection: str = "auto",
    field_trail: int = None,
):
    """Static 2×3 overview panel of a finished run.

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
    field_projection : "auto" | "xy" | "xz" | "yz"
    field_trail : int or None
        Length (in stored samples) of the fading orbit trail drawn over the
        field panel.  None → one tenth of the run.

    Returns
    -------
    fig
    """
    fig = plt.figure(figsize=(12.75, 8.1))
    gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.36, wspace=0.34)
    fig.suptitle(title, fontsize=13)

    registry = _build_panel_registry(
        history, field, q, m, length_unit, unit_label, normalize_v,
        field_components, field_density, spectrum_method,
        show_individual=show_individual,
        store_dt_warn_period=store_dt_warn_period,
        field_projection=field_projection, field_trail=field_trail,
        final_spectrum_only=True,
    )
    i = len(history.t) - 1
    for name, (row, col) in LAYOUT.items():
        ax    = fig.add_subplot(gs[row, col])
        panel = registry[name]
        _update_panel(panel, panel.build(ax, history), i, None, history)
    return fig


# ── GIF generation ────────────────────────────────────────────────────────────

def _save(fig, update, n_frames, filename, fps, writer, n_colors=128):
    """Render every frame and write a compact GIF.

    All frames share one palette and are quantised without dithering, so
    pixels that do not change between frames stay identical and the GIF
    encoder only stores the changed region of each frame.  (Per-frame
    adaptive palettes with dithering make every pixel flicker and inflate
    the file by an order of magnitude.)
    """
    from PIL import Image
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    if writer != "gif":
        raise ValueError(f"Unknown writer '{writer}'. Use 'gif'.")
    canvas = FigureCanvasAgg(fig)
    frames = []
    for k in range(n_frames):
        update(k)
        canvas.draw()
        frames.append(Image.fromarray(np.asarray(canvas.buffer_rgba())[..., :3].copy()))
    plt.close(fig)

    # One palette from a mosaic of evenly spaced frames
    pick = np.unique(np.linspace(0, n_frames - 1, min(8, n_frames)).astype(int))
    w, h = frames[0].size
    mosaic = Image.new("RGB", (w, h * len(pick)))
    for j, i in enumerate(pick):
        mosaic.paste(frames[i], (0, h * j))
    pal = mosaic.quantize(colors=n_colors, method=Image.Quantize.MEDIANCUT,
                          dither=Image.Dither.NONE)
    q = [f.quantize(palette=pal, dither=Image.Dither.NONE) for f in frames]
    q[0].save(filename, save_all=True, append_images=q[1:], loop=0,
              duration=int(round(1000 / fps)), disposal=1, optimize=False)
    print(f"Saved: {filename}")
    return filename


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
    spectrum_update_every: int = 3,
    field_update_every: int = 10,
    store_dt_warn_period: float = None,
    field_projection: str = "auto",
    field_trail: int = None,
    dpi: int = 90,
):
    """Generate a GIF for a single panel.

    Parameters
    ----------
    panel_name : one of "position_xy", "position_xz", "velocity_xy",
                 "velocity_xz", "field", "spectrum"
    filename : output path.  Defaults to f"{panel_name}.gif"
    writer : "gif" (Pillow).
    (other params: same as make_overview_gif)

    Returns
    -------
    filename : str
    """
    if panel_name not in LAYOUT:
        raise ValueError(f"Unknown panel '{panel_name}'. Choose from: {list(LAYOUT)}")

    registry = _build_panel_registry(
        history, field, q, m, length_unit, unit_label, normalize_v,
        field_components, field_density, spectrum_method,
        spectrum_update_every, field_update_every=field_update_every,
        store_dt_warn_period=store_dt_warn_period,
        field_projection=field_projection, field_trail=field_trail,
    )
    panel = registry[panel_name]

    fig, ax = plt.subplots(figsize=(5.5, 4.7), dpi=dpi)
    artists = panel.build(ax, history)
    fig.tight_layout()
    frames  = np.round(np.linspace(0, len(history.t) - 1, n_frames)).astype(int)

    def update(k):
        return _update_panel(panel, artists, frames[k], k, history)

    return _save(fig, update, len(frames), filename or f"{panel_name}.gif",
                 fps, writer)


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
    spectrum_update_every: int = 3,
    field_update_every: int = 10,
    store_dt_warn_period: float = None,
    field_projection: str = "auto",
    field_trail: int = None,
    title: str = "Overview",
    dpi: int = 70,
):
    """Generate a combined overview GIF with up to 6 panels.

    Parameters
    ----------
    panels : list of str or None
        Subset of panel names to include.  None → all 6.
    filename : output path
    n_frames : number of GIF frames (evenly spaced over the stored history)
    spectrum_update_every : recompute the FFT every this many frames
    field_update_every : re-trace time-dependent non-uniform fields every
        this many frames (uniform fields update every frame)
    writer : "gif".

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
        field_projection=field_projection, field_trail=field_trail,
    )

    fig = plt.figure(figsize=(12.75, 8.1), dpi=dpi)
    gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.36, wspace=0.34)
    fig.suptitle(title, fontsize=13)
    artists_map = {}
    for name in panels:
        row, col = LAYOUT[name]
        artists_map[name] = registry[name].build(fig.add_subplot(gs[row, col]),
                                                 history)

    frames = np.round(np.linspace(0, len(history.t) - 1, n_frames)).astype(int)

    def update(k):
        i = frames[k]
        out = []
        for name in panels:
            out += _update_panel(registry[name], artists_map[name], i, k, history)
        return out

    return _save(fig, update, len(frames), filename, fps, writer)
