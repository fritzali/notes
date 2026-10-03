"""
magparsol/fieldlines.py
------------------------
Field-line tracing and plotting for any FieldModel.

Rendering strategy (chosen automatically from the field type)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
EarthDipole
    Field lines are traced from their magnetic-equator crossings at the
    requested L-shells, in the magnetic meridian plane(s), down to the
    Earth's surface.  Because the meridian plane contains the (tilted)
    dipole axis, every line is an exact closed loop  r = L cos²λ.
    The Earth disk, the rotation axis and the magnetic axis are drawn.

Spatially uniform fields (``is_uniform=True``)
    A regular grid of equal-length arrows (``quiver``).  Field components
    normal to the viewing plane are shown as ⊙ (out of the plane) and
    ⊗ (into the plane).  With both B and E, the E grid is offset by half
    a cell so the two families never overlap.

Other non-uniform fields
    2-D projections use ``streamplot`` of the in-plane components;
    3-D views (or explicit ``seed_points``) use RK4 line tracing.

Field lines satisfy  dr/ds = F(r)/|F(r)|  with F = B or E.  Tracing is a
fixed-step RK4 on the unit direction field, both ways from every seed,
stopping at ``r_min`` (planet surface), ``r_max`` or ``max_steps``.

Projections
~~~~~~~~~~~
"xy", "xz", "yz" (2-D)  or  "3d".
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgba
from mpl_toolkits.mplot3d import Axes3D   # noqa: F401
from magparsol.constants import R_EARTH, B_FLOOR
from magparsol import style

_DEFAULT_COLORS = {"B": style.B_COLOR, "E": style.E_COLOR}

_PLANES = {
    "xy": (np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]), ("x", "y")),
    "xz": (np.array([1.0, 0.0, 0.0]), np.array([0.0, 0.0, 1.0]), ("x", "z")),
    "yz": (np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0]), ("y", "z")),
}


def _plane(projection):
    """Return (e1, e2, n, axis labels) of a 2-D projection plane."""
    if projection not in _PLANES:
        raise ValueError(f"Unknown projection '{projection}'. "
                         f"Choose from {list(_PLANES) + ['3d']}.")
    e1, e2, labels = _PLANES[projection]
    return e1, e2, np.cross(e1, e2), labels


def _is_dipole(field):
    from magparsol.fields import EarthDipole
    return isinstance(field, EarthDipole)


def _component(field, r, t, comp):
    B, E = field(r, t)
    return B if comp == "B" else E


# ── Seeding helpers ───────────────────────────────────────────────────────────

def _seeds_sphere(n: int, radius: float) -> np.ndarray:
    """n seeds distributed roughly uniformly on a sphere (Fibonacci spiral)."""
    golden = np.pi * (3.0 - np.sqrt(5.0))
    i      = np.arange(n)
    y      = 1.0 - 2.0 * i / (n - 1) if n > 1 else np.array([0.0])
    r_xy   = np.sqrt(np.clip(1.0 - y**2, 0, 1))
    theta  = golden * i
    return radius * np.stack([r_xy*np.cos(theta), y, r_xy*np.sin(theta)], axis=1)


def _seeds_box(n: int, bbox_min: np.ndarray, bbox_max: np.ndarray) -> np.ndarray:
    """About n seeds on a regular grid spanning the bounding box."""
    n_side = max(2, int(round(n ** (1.0/3.0))))
    axes   = [np.linspace(bbox_min[i], bbox_max[i], n_side) for i in range(3)]
    grid   = np.array(np.meshgrid(*axes, indexing="ij")).reshape(3, -1).T
    if len(grid) > n:
        idx  = np.round(np.linspace(0, len(grid)-1, n)).astype(int)
        grid = grid[idx]
    return grid


def _seeds_dipole_lshells(field, L_shells=(2, 3, 4, 5, 6, 8),
                          azimuths=(0.0, np.pi), unit: float = R_EARTH):
    """Seeds on the magnetic equator at radius L·R_E.

    The equator is the plane perpendicular to ``field.axis``.  Azimuth ψ is
    measured from x̂ towards ``axis × x̂`` (the in-equator direction lying in
    the y-z plane), so ψ = 0, π gives the meridian containing x̂ and
    ψ = ±π/2 the meridian lying in the y-z plane.

    Returns
    -------
    seeds : ndarray (n, 3)  [m]
    L_of_seed : ndarray (n,)
    """
    e_a = np.array([1.0, 0.0, 0.0])
    e_b = np.cross(field.axis, e_a)
    e_b /= np.linalg.norm(e_b)
    seeds, Ls = [], []
    for L in L_shells:
        for psi in azimuths:
            d = np.cos(psi) * e_a + np.sin(psi) * e_b
            seeds.append(L * unit * d)
            Ls.append(L)
    return np.array(seeds), np.array(Ls, dtype=float)


def _auto_bbox(history) -> tuple:
    """Bounding box of the trajectory with 20% margin."""
    r = history.r.reshape(-1, 3)
    return r.min(axis=0) * 1.2, r.max(axis=0) * 1.2


# ── Single field-line trace (RK4 on unit direction) ──────────────────────────

def _trace_line(field, r0, t: float, component: str, ds: float,
                max_steps: int, r_max=None, direction: int = 1,
                r_min=None) -> np.ndarray:
    """Trace one field line from r0 along (+1) or against (-1) the field.

    Stops at a null point, when |r| < r_min (the last point is placed on the
    sphere |r| = r_min), when |r| > r_max, or after max_steps.

    Returns
    -------
    points : ndarray, shape (n_points, 3)  [m]
    """
    floor = np.sqrt(B_FLOOR)

    def f_hat(pos):
        fv = _component(field, pos[None, :], t, component)[0]
        fm = np.linalg.norm(fv)
        return direction * fv / fm if fm > floor else None

    r      = np.asarray(r0, dtype=float).copy()
    points = [r.copy()]
    for _ in range(max_steps):
        k1 = f_hat(r)
        if k1 is None:
            break
        k2 = f_hat(r + 0.5*ds*k1)
        k3 = f_hat(r + 0.5*ds*k2) if k2 is not None else None
        k4 = f_hat(r + ds*k3) if k3 is not None else None
        if k4 is None:
            break
        r_new = r + (ds/6.0) * (k1 + 2*k2 + 2*k3 + k4)

        if r_min is not None and np.linalg.norm(r_new) < r_min:
            # Linear interpolation onto the surface |r| = r_min
            a, b = np.linalg.norm(r), np.linalg.norm(r_new)
            s = (a - r_min) / (a - b) if a != b else 1.0
            points.append(r + s * (r_new - r))
            break
        r = r_new
        points.append(r.copy())
        if r_max is not None and np.linalg.norm(r) > r_max:
            break
    return np.array(points)


def _trace_both(field, seed, t, comp, ds, max_steps, r_max, r_min):
    """Trace both ways and join into one line running along the field."""
    fwd = _trace_line(field, seed, t, comp, ds, max_steps, r_max, +1, r_min)
    bwd = _trace_line(field, seed, t, comp, ds, max_steps, r_max, -1, r_min)
    return np.vstack([bwd[::-1], fwd[1:]])


def _add_direction_arrow(ax, pts2d, color, frac=0.5, size=11):
    """Put one arrowhead at fraction ``frac`` of a 2-D polyline's arc length."""
    if len(pts2d) < 3:
        return
    seg = np.linalg.norm(np.diff(pts2d, axis=0), axis=1)
    s   = np.concatenate([[0.0], np.cumsum(seg)])
    if s[-1] <= 0:
        return
    i = int(np.clip(np.searchsorted(s, frac * s[-1]), 1, len(pts2d) - 2))
    ax.annotate("", xy=pts2d[i+1], xytext=pts2d[i-1],
                arrowprops=dict(arrowstyle="-|>", color=color, lw=0,
                                mutation_scale=size, shrinkA=0, shrinkB=0))


# ── Decorations ───────────────────────────────────────────────────────────────

def draw_earth_2d(ax, length_unit, zorder=5):
    """Filled Earth disk of radius R_E (in display units)."""
    ax.add_patch(Circle((0, 0), R_EARTH / length_unit,
                        facecolor=style.earth_fill(), edgecolor="black",
                        lw=0.6, zorder=zorder))


def _draw_earth_2d(ax, field, e1, e2, length_unit, extent):
    """Earth disk, rotation axis (dotted) and magnetic axis (dashed)."""
    draw_earth_2d(ax, length_unit)
    axes = []
    for vec, ls, lbl in ((np.array([0.0, 0.0, 1.0]), ":", "Rotation Axis"),
                         (field.axis, "--", "Magnetic Axis")):
        p = np.array([vec @ e1, vec @ e2])
        if np.linalg.norm(p) < 0.2:
            continue   # axis (nearly) normal to the viewing plane
        axes.append((p / np.linalg.norm(p), ls, lbl))
    if len(axes) == 2 and abs(axes[0][0] @ axes[1][0]) > np.cos(np.radians(0.5)):
        # Both axes project onto the same line: draw and label it once
        axes = [(axes[1][0], "--", "Rotation & Magnetic Axis")]
    for p, ls, lbl in axes:
        p = p * extent
        ax.plot([-p[0], p[0]], [-p[1], p[1]], ls, color=style.AXIS_COLOR,
                lw=0.8, zorder=4, label=lbl)


def _draw_earth_3d(ax, length_unit):
    R = R_EARTH / length_unit
    u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
    ax.plot_surface(R*np.cos(u)*np.sin(v), R*np.sin(u)*np.sin(v), R*np.cos(v),
                    color=style.EARTH_FACE, alpha=0.45, linewidth=0, shade=True)


# ── Fading orbit trail (static and animated) ──────────────────────────────────

class FadingTrail:
    """A line whose opacity ramps from 0 (oldest) to ``alpha`` (newest).

    Only the last ``length`` points are drawn, so an orbit plotted over a
    field configuration does not hide the field.  Call :meth:`set_data` with
    the full path up to the current point.
    """

    def __init__(self, ax, length=200, color=style.ORBIT_COLOR, alpha=0.8,
                 lw=1.0, zorder=6):
        self.length = int(length)
        self.rgba   = np.array(to_rgba(color))
        self.alpha  = alpha
        self.lc     = LineCollection([], linewidths=lw, zorder=zorder,
                                     capstyle="round")
        ax.add_collection(self.lc)

    def set_data(self, x, y):
        x, y = np.asarray(x)[-self.length:], np.asarray(y)[-self.length:]
        if len(x) < 2:
            self.lc.set_segments([])
            return self.lc
        pts  = np.column_stack([x, y])
        segs = np.stack([pts[:-1], pts[1:]], axis=1)
        c    = np.tile(self.rgba, (len(segs), 1))
        c[:, 3] = self.alpha * np.linspace(0.0, 1.0, len(segs))**1.5
        self.lc.set_segments(segs)
        self.lc.set_color(c)
        return self.lc


def field_legend(ax, handles, n_cols=None):
    """Legend in a row above the axes (never covers the field)."""
    if not handles:
        return None
    leg = style.legend(ax, handles=handles, loc="lower center",
                       bbox_to_anchor=(0.5, 1.01),
                       ncol=n_cols or min(len(handles), 4), fontsize=8,
                       handlelength=1.6, columnspacing=1.0, borderpad=0.35)
    ax.set_title(ax.get_title(), pad=26)
    return leg


# ── Uniform-field arrow grids (re-usable by the animation module) ─────────────

class UniformFieldArrows:
    """Quiver grid (+ ⊙/⊗ markers) for one component of a uniform field.

    Arrow length is proportional to |F_in-plane| / ``ref_mag`` so a
    time-dependent field (e.g. ``CyclotronWaveField``) can be animated with
    :meth:`update` while keeping a fixed length scale.
    """

    def __init__(self, ax, field, comp, t, e1, e2, lim, n_grid=7,
                 offset=0.0, color="black", ref_mag=None):
        self.ax, self.field, self.comp = ax, field, comp
        self.e1, self.e2, self.n = e1, e2, np.cross(e1, e2)
        spacing = 2 * lim / n_grid
        g = -lim + spacing * (np.arange(n_grid) + 0.5 + offset)
        g = g[np.abs(g) < lim]
        self.X, self.Y = np.meshgrid(g, g)
        self.arrow_len = 0.75 * spacing
        F0 = _component(field, np.zeros((1, 3)), t, comp)[0]
        self.ref_mag = ref_mag if ref_mag else max(np.linalg.norm(F0), 1e-300)
        zeros = np.zeros_like(self.X)
        self.quiver = ax.quiver(self.X, self.Y, zeros, zeros, color=color,
                                angles="xy", scale_units="xy", scale=1.0,
                                width=0.005, headwidth=4, headlength=5,
                                pivot="middle", zorder=3)
        kw = dict(ls="none", color=color, mfc="none", ms=8, mew=1.1, zorder=3)
        self.ring,  = ax.plot([], [], marker="o", **kw)
        self.dot,   = ax.plot([], [], marker=".", **{**kw, "ms": 4})
        self.cross, = ax.plot([], [], marker="x", **{**kw, "ms": 5})
        self.update(t)

    def update(self, t):
        F = _component(self.field, np.zeros((1, 3)), t, self.comp)[0]
        mag = np.linalg.norm(F)
        f1, f2, fn = F @ self.e1, F @ self.e2, F @ self.n
        s = self.arrow_len / self.ref_mag
        U = np.full_like(self.X, f1 * s)
        V = np.full_like(self.Y, f2 * s)
        self.quiver.set_UVC(U, V)
        # Out-of-plane symbol when the normal component is significant
        show = mag > 0 and abs(fn) > 0.05 * mag and np.hypot(f1, f2) < 0.95 * mag
        x, y = (self.X.ravel(), self.Y.ravel()) if show else ([], [])
        self.ring.set_data(x, y)
        self.dot.set_data(x if fn > 0 else [], y if fn > 0 else [])
        self.cross.set_data(x if fn < 0 else [], y if fn < 0 else [])
        return [self.quiver, self.ring, self.dot, self.cross]


def _ref_magnitude(field, comp, t, ref_times=None):
    """Largest |F| over ``ref_times`` (for a fixed arrow scale), or None."""
    if field.is_static:
        return None
    if ref_times is None:
        w = getattr(field, "omega_c", None)
        ref_times = t + np.linspace(0, 2*np.pi / w, 33) if w else [t]
    ts = np.asarray(ref_times, dtype=float)
    if len(ts) > 200:
        ts = ts[np.linspace(0, len(ts) - 1, 200).astype(int)]
    mags = [np.linalg.norm(_component(field, np.zeros((1, 3)), tt, comp)[0])
            for tt in ts]
    return max(mags) or None


def _uniform_3d(ax, field, comp, t, lim, color, n=4):
    F = _component(field, np.zeros((1, 3)), t, comp)[0]
    mag = np.linalg.norm(F)
    if mag == 0:
        return
    g = np.linspace(-lim, lim, n) * 0.8
    X, Y, Z = np.meshgrid(g, g, g)
    L = 0.6 * (g[1] - g[0])
    f = F / mag * L
    ax.quiver(X, Y, Z, f[0], f[1], f[2], color=color, alpha=0.8,
              arrow_length_ratio=0.3, pivot="middle", linewidth=0.9)


# ── Main public function ──────────────────────────────────────────────────────

def plot_field_lines(
    field,
    history=None,
    components=("B",),
    density: str = "auto",
    seed_points=None,
    seed_strategy: str = "auto",
    n_seeds: int = None,
    t: float = 0.0,
    ds: float = None,
    max_steps: int = 4000,
    r_max: float = None,
    r_min: float = None,
    length_unit: float = 1.0,
    unit_label: str = "m",
    ax_lim=None,
    projection: str = "3d",
    earth_sphere=None,
    color=None,
    ax=None,
    title: str = None,
    L_shells=(2, 3, 4, 6, 8),
    arrows: bool = True,
    legend: bool = True,
    ref_times=None,
    orbit_alpha: float = 0.35,
    orbit_trail: int = None,
):
    """Trace and plot field lines (or arrow grids) for any FieldModel.

    Parameters
    ----------
    field : FieldModel
    history : TrajectoryHistory or None
        Used for automatic axis limits / seeding and overlaid if given.
    components : tuple of "B" and/or "E"
    density : "low" | "medium" | "high" | "auto"
        Number of seeds / arrows / streamline density.
    seed_points : array-like (n, 3) or None
        Explicit seed positions [m] (forces RK4 line tracing).
    seed_strategy : "auto" | "sphere" | "box" | "dipole_Lshells"
        Only used for non-uniform fields in 3-D or with explicit tracing.
    n_seeds : int or None
        Override the density-derived seed count (sphere / box seeding).
    t : float
        Time at which the field is evaluated [s].
    ds : float or None
        Arc-length step for tracing [m].  Auto-estimated if None.
    max_steps : int
    r_max, r_min : float or None
        Outer / inner tracing boundary [m].  For EarthDipole r_min defaults
        to R_EARTH (lines end on the surface) and r_max to 1.5·max(L)·R_E.
    length_unit : float
        Divisor for display units (e.g. R_EARTH).
    unit_label : str
    ax_lim : float or None
        Symmetric axis limit in display units.
    projection : "3d" | "xy" | "xz" | "yz"
    earth_sphere : bool or None
        Draw the Earth.  None → True for EarthDipole.
    color : str, dict or None
        None → B blue, E red.  Dict: {"B": ..., "E": ...}.
    ax : matplotlib Axes or None
    title : str or None
    L_shells : tuple of float
        L values for dipole lines.
    arrows : bool
        Add direction arrowheads to traced lines.
    legend : bool
    ref_times : array-like or None
        Times used to fix the arrow length scale of a time-dependent
        uniform field (arrow length ∝ |F(t)| / max |F(ref_times)|).
        None → one wave period if the field has ``omega_c``, else just t.
    orbit_alpha : float
        Opacity of the overlaid ``history`` orbit (kept low so the field
        stays visible).
    orbit_trail : int or None
        If given (2-D only), draw only the last ``orbit_trail`` samples of the
        orbit as a trail that fades out towards its tail.

    Returns
    -------
    fig, ax
    """
    dens = {"low": 0, "medium": 1, "high": 2, "auto": 1}.get(density, 1)
    if n_seeds is None:
        n_seeds = (6, 12, 24)[dens]
    if color is None:
        color = dict(_DEFAULT_COLORS)
    elif isinstance(color, str):
        color = {c: color for c in components}
    is_dipole = _is_dipole(field)
    if earth_sphere is None:
        earth_sphere = is_dipole
    if is_dipole:
        r_min = R_EARTH if r_min is None else r_min
        if r_max is None:
            r_max = 1.5 * max(L_shells) * R_EARTH

    is3d = projection == "3d"

    # ── Axis limit in display units ───────────────────────────────────────────
    if ax_lim is None:
        if history is not None:
            r_flat = history.r.reshape(-1, 3) / length_unit
            ax_lim = float(np.max(np.abs(r_flat))) * 1.15
        elif is_dipole:
            ax_lim = 1.1 * max(L_shells) * R_EARTH / length_unit
        elif seed_points is not None:
            ax_lim = float(np.max(np.abs(seed_points))) * 1.2 / length_unit
        else:
            ax_lim = 1.0
    ax_lim = float(ax_lim) if ax_lim > 0 else 1.0

    # ── Axes ──────────────────────────────────────────────────────────────────
    own_fig = ax is None
    if own_fig:
        fig = plt.figure(figsize=(6, 5.5))
        ax  = fig.add_subplot(111, projection="3d" if is3d else None)
    else:
        fig = ax.get_figure()

    if not is3d:
        e1, e2, n_hat, labels = _plane(projection)
    handles = []

    # ── Draw each component ───────────────────────────────────────────────────
    ax._mps_uniform_arrows = []
    for k, comp in enumerate(components):
        col = color.get(comp, "grey")
        if field.is_uniform:
            if is3d:
                _uniform_3d(ax, field, comp, t, ax_lim, col, n=(3, 4, 5)[dens])
            else:
                arr = UniformFieldArrows(
                    ax, field, comp, t, e1, e2, ax_lim,
                    n_grid=(5, 7, 9)[dens],
                    offset=0.5 * k if len(components) > 1 else 0.0,
                    color=col,
                    ref_mag=_ref_magnitude(field, comp, t, ref_times))
                ax._mps_uniform_arrows.append(arr)
            handles.append(Line2D([], [], color=col, marker=r"$\rightarrow$",
                                  ms=12, ls="none", label=f"$\\mathbf{{{comp}}}$"))
            continue

        if is_dipole and comp == "E":
            continue   # dipole has no electric field

        if is_dipole and seed_points is None:
            if is3d:
                # Offset so no meridian is seen edge-on from the default view
                az = np.linspace(0, 2*np.pi, (4, 6, 8)[dens], endpoint=False) \
                    + np.pi / 9
            elif projection == "yz":
                az = (np.pi/2, -np.pi/2)
            else:   # xz (true meridian ⊃ x̂) or xy (top view)
                az = (0.0, np.pi) if projection == "xz" else \
                     np.linspace(0, 2*np.pi, 8, endpoint=False)
            seeds, Ls = _seeds_dipole_lshells(field, L_shells, az)
            for seed, L in zip(seeds, Ls):
                step = ds if ds is not None else 0.01 * L * R_EARTH
                pts  = _trace_both(field, seed, t, comp, step, max_steps,
                                   r_max, r_min) / length_unit
                if is3d:
                    ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], color=col,
                            lw=0.9, alpha=0.8)
                else:
                    p2 = np.column_stack([pts @ e1, pts @ e2])
                    ax.plot(p2[:, 0], p2[:, 1], color=col, lw=1.0, zorder=2)
                    if arrows:
                        _add_direction_arrow(ax, p2, col)
            handles.append(Line2D([], [], color=col, label=r"$\mathbf{B}$"))
            continue

        if not is3d and seed_points is None:
            # Generic 2-D view: streamlines of the in-plane components
            n_g = 60
            g = np.linspace(-ax_lim, ax_lim, n_g)
            X, Y = np.meshgrid(g, g)
            P = (X.ravel()[:, None] * e1 + Y.ravel()[:, None] * e2) * length_unit
            F = _component(field, P, t, comp)
            U = (F @ e1).reshape(X.shape)
            V = (F @ e2).reshape(X.shape)
            mag = np.hypot(U, V)
            if r_min is not None:
                inside = np.linalg.norm(P, axis=1).reshape(X.shape) < r_min
                U = np.ma.array(U, mask=inside)
                V = np.ma.array(V, mask=inside)
            if np.nanmax(mag) > 0:
                ax.streamplot(X, Y, U, V, color=col, density=(0.7, 1.1, 1.6)[dens],
                              linewidth=0.9, arrowsize=1.0, zorder=2)
                handles.append(Line2D([], [], color=col,
                                      label=f"$\\mathbf{{{comp}}}$"))
            continue

        # Generic RK4 tracing (3-D, or explicit seeds)
        if seed_points is not None:
            seeds = np.asarray(seed_points, dtype=float)
        else:
            strategy = "sphere" if seed_strategy == "auto" else seed_strategy
            if strategy == "box":
                if history is not None:
                    bmin, bmax = _auto_bbox(history)
                else:
                    s = ax_lim * length_unit
                    bmin, bmax = -s*np.ones(3), s*np.ones(3)
                seeds = _seeds_box(n_seeds, bmin, bmax)
            else:
                seeds = _seeds_sphere(n_seeds, 0.5 * ax_lim * length_unit)
        step = ds if ds is not None else \
            max(float(np.max(np.ptp(seeds, axis=0))), ax_lim*length_unit) / 300.0
        rmax = r_max if r_max is not None else 1.8 * ax_lim * length_unit
        for seed in seeds:
            pts = _trace_both(field, seed, t, comp, step, max_steps, rmax,
                              r_min) / length_unit
            if is3d:
                ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], color=col, lw=0.8,
                        alpha=0.8)
            else:
                p2 = np.column_stack([pts @ e1, pts @ e2])
                ax.plot(p2[:, 0], p2[:, 1], color=col, lw=0.9)
                if arrows:
                    _add_direction_arrow(ax, p2, col)
        handles.append(Line2D([], [], color=col, label=f"$\\mathbf{{{comp}}}$"))

    # ── Earth, trajectory, formatting ─────────────────────────────────────────
    if earth_sphere:
        if is3d:
            _draw_earth_3d(ax, length_unit)
        else:
            _draw_earth_2d(ax, field, e1, e2, length_unit, ax_lim)

    if history is not None:
        r = history.r / length_unit
        for pid in range(r.shape[1]):
            if is3d:
                ax.plot(r[:, pid, 0], r[:, pid, 1], r[:, pid, 2],
                        color=style.ORBIT_COLOR, lw=0.5, alpha=orbit_alpha,
                        rasterized=style.raster(len(r)))
            elif orbit_trail:
                FadingTrail(ax, orbit_trail).set_data(r[:, pid] @ e1, r[:, pid] @ e2)
            else:
                # Below the opaque Earth disk (zorder 5): the orbit is hidden
                # where its projection passes behind the planet
                ax.plot(r[:, pid] @ e1, r[:, pid] @ e2, color=style.ORBIT_COLOR,
                        lw=0.5, alpha=orbit_alpha, zorder=4.5,
                        rasterized=style.raster(len(r)))

    if is3d:
        ax.set_xlim3d(-ax_lim, ax_lim)
        ax.set_ylim3d(-ax_lim, ax_lim)
        ax.set_zlim3d(-ax_lim, ax_lim)
        ax.set_box_aspect((1, 1, 1))
        style.style_3d(ax)
        ax.set_xlabel(f"$x$ [{unit_label}]")
        ax.set_ylabel(f"$y$ [{unit_label}]")
        ax.set_zlabel(f"$z$ [{unit_label}]")
    else:
        ax.set_xlim(-ax_lim, ax_lim)
        ax.set_ylim(-ax_lim, ax_lim)
        ax.set_aspect("equal")
        ax.set_xlabel(f"${labels[0]}$ [{unit_label}]")
        ax.set_ylabel(f"${labels[1]}$ [{unit_label}]")

    if title is None:
        title = getattr(field, "display_name", "Field")
    ax.set_title(title)

    if legend and (len(handles) > 1 or (earth_sphere and not is3d)):
        extra = ax.get_legend_handles_labels()[0]
        if is3d:
            style.legend(ax, handles=handles + extra, loc="upper right", fontsize=8)
        else:
            field_legend(ax, handles + extra)

    if own_fig:
        fig.tight_layout()
    return fig, ax

