"""
magparsol/fieldlines.py
------------------------
Field-line tracing and plotting for any FieldModel.

Tracing
~~~~~~~
Field lines satisfy  dr/ds = F(r)/|F(r)|  where F is B or E.
Integration is a simple fixed-step RK4 on the unit-direction field.
Both directions (+/−) are traced from each seed so closed lines (dipole
loops) are drawn completely.

Seeding strategies
~~~~~~~~~~~~~~~~~~
auto        — bounding-box aspect ratio decides sphere vs box
sphere      — seeds on a sphere of given radius around origin
box         — seeds on a regular grid within a bounding box
dipole_Lshells — (EarthDipole only) seeds at standard L-shell crossing
                  points in the magnetic equatorial plane

Uniform-field handling
~~~~~~~~~~~~~~~~~~~~~~
Spatially uniform fields (is_uniform=True) skip streamline tracing and
instead draw short representative arrow segments at the seed locations.
This is cheaper and avoids meaningless parallel-line clutter.

B and E
~~~~~~~
Pass components=("B",), ("E",), or ("B","E").  Colors default to
steel-blue for B and orange-red for E.  A legend is added when both
components are shown.
"""

import warnings
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D   # noqa: F401
from magparsol.constants import R_EARTH, B_FLOOR


# ── Seeding helpers ───────────────────────────────────────────────────────────

def _seeds_sphere(n: int, radius: float) -> np.ndarray:
    """n seeds distributed roughly uniformly on a sphere of given radius."""
    # Fibonacci / golden-ratio spiral for near-uniform coverage
    golden = np.pi * (3.0 - np.sqrt(5.0))
    i      = np.arange(n)
    y      = 1.0 - 2.0 * i / (n - 1) if n > 1 else np.array([0.0])
    r_xy   = np.sqrt(np.clip(1.0 - y**2, 0, 1))
    theta  = golden * i
    x      = r_xy * np.cos(theta)
    z      = r_xy * np.sin(theta)
    return radius * np.stack([x, y, z], axis=1)   # (n, 3)


def _seeds_box(n: int, bbox_min: np.ndarray, bbox_max: np.ndarray) -> np.ndarray:
    """n seeds on a regular grid spanning the bounding box."""
    n_side = max(2, int(round(n ** (1.0/3.0))))
    axes   = [np.linspace(bbox_min[i], bbox_max[i], n_side) for i in range(3)]
    grid   = np.array(np.meshgrid(*axes, indexing="ij")).reshape(3, -1).T
    # subsample to exactly n if grid is larger
    if len(grid) > n:
        idx  = np.round(np.linspace(0, len(grid)-1, n)).astype(int)
        grid = grid[idx]
    return grid


def _seeds_dipole_lshells(L_shells=(2, 3, 4, 5, 6, 8),
                            n_phi: int = 1,
                            unit: float = R_EARTH) -> np.ndarray:
    """Seeds at L-shell equatorial crossings in the noon-midnight meridian.

    For each L value, place n_phi seeds evenly in azimuth.
    Default n_phi=1 gives the noon meridian only (x-z plane).
    """
    seeds = []
    phis  = np.linspace(0, 2*np.pi, n_phi, endpoint=False)
    for L in L_shells:
        r = L * unit
        for phi in phis:
            seeds.append([r * np.cos(phi), r * np.sin(phi), 0.0])
    return np.array(seeds)


def _auto_strategy(field, history=None) -> str:
    """Choose 'sphere' or 'box' based on bounding-box aspect ratio."""
    from magparsol.fields import EarthDipole
    if isinstance(field, EarthDipole):
        return "sphere"
    if history is not None:
        r = history.r.reshape(-1, 3)
        extents = r.max(axis=0) - r.min(axis=0)
        extents = np.where(extents > 0, extents, 1.0)
        ratio   = extents.max() / extents.min()
        return "box" if ratio > 3.0 else "sphere"
    return "sphere"


def _auto_radius(field, history=None) -> float:
    """Estimate a sensible seed radius from trajectory or a default."""
    if history is not None:
        r_mag = np.linalg.norm(history.r.reshape(-1, 3), axis=1)
        return float(r_mag.mean())
    return R_EARTH * 5.0


def _auto_bbox(history) -> tuple:
    """Bounding box of the trajectory with 20% margin."""
    r    = history.r.reshape(-1, 3)
    bmin = r.min(axis=0) * 1.2
    bmax = r.max(axis=0) * 1.2
    return bmin, bmax


# ── Single field-line trace (RK4 on unit direction) ──────────────────────────

def _trace_line(field, r0: np.ndarray, t: float, component: str,
                ds: float, max_steps: int, r_max: float,
                direction: int = 1) -> np.ndarray:
    """Trace one field line from r0 in the given direction (+1 or -1).

    Parameters
    ----------
    component : "B" | "E"
    direction : +1 (along field) | -1 (against field)

    Returns
    -------
    points : ndarray, shape (n_points, 3)
    """
    r      = np.asarray(r0, dtype=float).copy()
    points = [r.copy()]
    warned = False

    for _ in range(max_steps):
        # Evaluate field at current position
        r2d        = r[None, :]
        B_val = field(r2d, t)
        F          = B_val[0] if component == "B" else E_val[0]
        F_mag      = np.linalg.norm(F)

        if F_mag < np.sqrt(B_FLOOR):
            break   # null point — stop tracing

        f_hat = direction * F / F_mag

        # RK4 on dr/ds = f_hat(r)
        def drdp(pos):
            b_ = field(pos[None, :], t)
            fv     = b_[0] if component == "B" else e_[0]
            fm     = np.linalg.norm(fv)
            return direction * fv / fm if fm > np.sqrt(B_FLOOR) else np.zeros(3)

        k1 = drdp(r)
        k2 = drdp(r + 0.5*ds*k1)
        k3 = drdp(r + 0.5*ds*k2)
        k4 = drdp(r + ds*k3)
        r  = r + (ds/6.0) * (k1 + 2*k2 + 2*k3 + k4)

        points.append(r.copy())

        if r_max is not None and np.linalg.norm(r) > r_max:
            break

    return np.array(points)


# ── Uniform-field arrow drawing ───────────────────────────────────────────────

def _draw_uniform_arrows(ax, field, component: str, seeds: np.ndarray,
                          t: float, color: str, length_unit: float,
                          projection: str, label: str):
    """Draw short arrows for spatially uniform fields."""
    r_probe    = np.zeros((1, 3))
    B_val = field(r_probe, t)
    F          = B_val[0] if component == "B" else E_val[0]
    F_mag      = np.linalg.norm(F)
    if F_mag < 1e-40:
        return   # truly zero field

    # Arrow length: 15% of the seed spread, or a fixed scale
    spread = np.ptp(seeds, axis=0)
    scale  = float(np.max(spread)) * 0.15 if np.max(spread) > 0 else length_unit
    f_hat  = F / F_mag * scale / length_unit

    for seed in seeds:
        s = seed / length_unit
        if projection == "3d":
            ax.quiver(s[0], s[1], s[2], f_hat[0], f_hat[1], f_hat[2],
                      color=color, alpha=0.8, label=label)
        elif projection == "xy":
            ax.annotate("", xy=(s[0]+f_hat[0], s[1]+f_hat[1]), xytext=(s[0], s[1]),
                        arrowprops=dict(arrowstyle="->", color=color))
        else:  # xz
            ax.annotate("", xy=(s[0]+f_hat[0], s[2]+f_hat[2]), xytext=(s[0], s[2]),
                        arrowprops=dict(arrowstyle="->", color=color))
        label = None   # only label the first arrow


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
    max_steps: int = 2000,
    r_max: float = None,
    length_unit: float = 1.0,
    unit_label: str = "m",
    ax_lim=None,
    projection: str = "3d",
    earth_sphere: bool = False,
    color=None,
    ax=None,
    title: str = "Field lines",
):
    """Trace and plot field lines for any FieldModel.

    Parameters
    ----------
    field : FieldModel
    history : TrajectoryHistory or None
        Used for auto-seeding and auto axis limits.
    components : tuple of "B" and/or "E"
    density : "low" | "medium" | "high" | "auto"
        Controls default n_seeds: low=4, medium=10, high=24, auto=medium.
    seed_points : array-like (n, 3) or None
        Explicit seed positions [m].  Overrides seed_strategy/n_seeds.
    seed_strategy : "auto" | "sphere" | "box" | "dipole_Lshells"
    n_seeds : int or None
        Override density-derived seed count.
    t : float
        Time snapshot for field evaluation [s].
    ds : float or None
        Arc-length step for tracing.  Auto-estimated if None.
    max_steps : int
    r_max : float or None
        Stop tracing when |r| > r_max.  None → rely on max_steps only.
    length_unit : float
        Divisor for display units.
    unit_label : str
    ax_lim : float or None
    projection : "3d" | "xy" | "xz"
    earth_sphere : bool
    color : str or dict or None
        None → auto (B: steelblue, E: darkorange).
        Dict: {"B": "...", "E": "..."}.
    ax : matplotlib Axes or None
    title : str

    Returns
    -------
    fig, ax
    """
    # ── Density → n_seeds ────────────────────────────────────────────────────
    _density_map = {"low": 4, "medium": 10, "high": 24, "auto": 10}
    if n_seeds is None:
        n_seeds = _density_map.get(density, 10)

    # ── Color defaults ────────────────────────────────────────────────────────
    if color is None:
        color = {"B": "steelblue", "E": "darkorange"}
    elif isinstance(color, str):
        color = {c: color for c in components}

    # ── Build seeds ───────────────────────────────────────────────────────────
    if seed_points is not None:
        seeds = np.asarray(seed_points, dtype=float)
    else:
        strategy = seed_strategy
        if strategy == "auto":
            strategy = _auto_strategy(field, history)

        if strategy == "dipole_Lshells":
            seeds = _seeds_dipole_lshells()
        elif strategy == "sphere":
            radius = _auto_radius(field, history)
            seeds  = _seeds_sphere(n_seeds, radius)
        else:   # box
            if history is not None:
                bmin, bmax = _auto_bbox(history)
            else:
                s = _auto_radius(field, None)
                bmin, bmax = -s*np.ones(3), s*np.ones(3)
            seeds = _seeds_box(n_seeds, bmin, bmax)

    # ── Auto ds ───────────────────────────────────────────────────────────────
    if ds is None:
        seed_spread = float(np.max(np.ptp(seeds, axis=0)))
        ds = seed_spread / 500.0 if seed_spread > 0 else R_EARTH * 0.01

    # ── Create axes ───────────────────────────────────────────────────────────
    own_fig = ax is None
    if own_fig:
        fig = plt.figure(figsize=(8, 7))
        if projection == "3d":
            ax = fig.add_subplot(111, projection="3d")
        else:
            ax = fig.add_subplot(111)
    else:
        fig = ax.get_figure()

    # ── Earth sphere ──────────────────────────────────────────────────────────
    if earth_sphere and projection == "3d":
        u_e, v_e = np.mgrid[0:2*np.pi:30j, 0:np.pi:30j]
        xs = (R_EARTH/length_unit) * np.cos(u_e) * np.sin(v_e)
        ys = (R_EARTH/length_unit) * np.sin(u_e) * np.sin(v_e)
        zs = (R_EARTH/length_unit) * np.cos(v_e)
        ax.plot_wireframe(xs, ys, zs, color="royalblue", alpha=0.25, linewidth=0.4)

    # ── Warn once if r_max is None ────────────────────────────────────────────
    if r_max is None and not field.is_uniform:
        warnings.warn(
            "r_max is not set. Field-line traces stop only at max_steps. "
            "Set r_max to avoid incomplete or very long traces.",
            UserWarning, stacklevel=2,
        )

    # ── Trace and plot for each component ────────────────────────────────────
    for comp in components:
        col   = color.get(comp, "gray")
        label = f"$\\mathbf{{{'B' if comp=='B' else 'E'}}}$ field"

        if field.is_uniform:
            _draw_uniform_arrows(ax, field, comp, seeds, t, col,
                                 length_unit, projection, label)
        else:
            first_line = True
            for seed in seeds:
                for direction in (+1, -1):
                    pts = _trace_line(field, seed, t, comp,
                                      ds, max_steps, r_max, direction)
                    if len(pts) < 2:
                        continue
                    pts_d = pts / length_unit
                    lbl = label if first_line else None
                    if projection == "3d":
                        ax.plot(pts_d[:,0], pts_d[:,1], pts_d[:,2],
                                color=col, lw=0.8, alpha=0.7, label=lbl)
                    elif projection == "xy":
                        ax.plot(pts_d[:,0], pts_d[:,1],
                                color=col, lw=0.8, alpha=0.7, label=lbl)
                    else:   # xz
                        ax.plot(pts_d[:,0], pts_d[:,2],
                                color=col, lw=0.8, alpha=0.7, label=lbl)
                    first_line = False

    # ── Axis formatting ───────────────────────────────────────────────────────
    if ax_lim is None and history is not None:
        r_flat = history.r.reshape(-1, 3) / length_unit
        ax_lim = float(np.max(np.abs(r_flat))) * 1.2

    if ax_lim is not None:
        if projection == "3d":
            ax.set_xlim3d(-ax_lim, ax_lim)
            ax.set_ylim3d(-ax_lim, ax_lim)
            ax.set_zlim3d(-ax_lim, ax_lim)
            ax.set_box_aspect((1, 1, 1))
            ax.set_xlabel(f"x [{unit_label}]")
            ax.set_ylabel(f"y [{unit_label}]")
            ax.set_zlabel(f"z [{unit_label}]")
        elif projection == "xy":
            ax.set_xlim(-ax_lim, ax_lim); ax.set_ylim(-ax_lim, ax_lim)
            ax.set_xlabel(f"x [{unit_label}]"); ax.set_ylabel(f"y [{unit_label}]")
            ax.set_aspect("equal")
        else:
            ax.set_xlim(-ax_lim, ax_lim); ax.set_ylim(-ax_lim, ax_lim)
            ax.set_xlabel(f"x [{unit_label}]"); ax.set_ylabel(f"z [{unit_label}]")
            ax.set_aspect("equal")

    ax.set_title(title)
    if len(components) > 1:
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(handles, labels, loc="upper right", fontsize=8)

    if own_fig:
        plt.tight_layout()
    return fig, ax
