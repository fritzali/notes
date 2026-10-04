"""
Plot style shared by all plotting functions of this package.

Colours are restricted to named matplotlib colours: black and greys plus
steelblue, firebrick, olivedrab, goldenrod and rebeccapurple. Call
:func:`use` once, e.g. at the top of a notebook, to set the matching
rcParams (colour cycle, legend frame, resolution).
"""

import numpy as np
import matplotlib as mpl
from cycler import cycler

B_COLOR = "steelblue"
E_COLOR = "firebrick"
ORBIT_COLOR = "black"           # orbits drawn over a field configuration
EARTH_FACE = "steelblue"
EARTH_ALPHA = 0.25              # share of EARTH_FACE in the Earth fill
AXIS_COLOR = "dimgrey"          # rotation and magnetic axis
REF_COLOR = "darkgrey"          # reference curves, e.g. the final spectrum

# Drawing order. The grid is drawn at matplotlib's default zorder 1.5,
# so the Earth hides grid and field lines but not the particle.
Z_FIELD = 2
Z_AXIS = 3
Z_EARTH = 4
Z_TRACK = 5
Z_MARKER = 6
Z_LABEL = 10

CYCLE = ["steelblue", "firebrick", "olivedrab", "goldenrod", "rebeccapurple"]

# Lines with more points than this are rasterised inside vector figures
RASTER_THRESHOLD = 4000

LEGEND_KW = dict(frameon=True, fancybox=False, edgecolor="black",
                 framealpha=1.0, borderpad=0.5)


def legend(ax, *args, **kwargs):
    """``ax.legend`` with a thin black square frame."""
    kw = {**LEGEND_KW, **kwargs}
    leg = ax.legend(*args, **kw)
    leg.get_frame().set_linewidth(0.6)
    return leg


def earth_fill():
    """Opaque light blue for the Earth disk: EARTH_FACE mixed into white.

    An opaque fill (rather than a transparent one) hides grid and field
    lines behind the disk.
    """
    from matplotlib.colors import to_rgb
    c = np.array(to_rgb(EARTH_FACE))
    return tuple(EARTH_ALPHA * c + (1 - EARTH_ALPHA) * np.ones(3))


def format_time(t, t_total=None):
    """Time as an integer in seconds, or in milliseconds for runs shorter
    than 10 s (``t_total``), so that short runs do not show "0 s"."""
    ref = t_total if t_total is not None else t
    if ref >= 10:
        return f"{t:.0f} s"
    return f"{1e3 * t:.0f} ms"


def format_time_pair(t, t_total):
    """``"t / total unit"`` with both numbers in the same unit."""
    if t_total >= 10:
        return f"{t:.0f} / {t_total:.0f} s"
    return f"{1e3 * t:.0f} / {1e3 * t_total:.0f} ms"


def raster(n_points):
    """Whether a line of ``n_points`` should be rasterised."""
    return n_points > RASTER_THRESHOLD


def style_3d(ax):
    """Light panes and grid lines for 3D axes."""
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.pane.set_facecolor((1, 1, 1, 0))
        a.pane.set_edgecolor("lightgrey")
        a._axinfo["grid"].update(color=(0.5, 0.5, 0.5, 0.25), linewidth=0.5)


def colors(n):
    """The first ``n`` palette colours, repeating if ``n`` > 5."""
    return [CYCLE[i % len(CYCLE)] for i in range(n)]


def use(dpi=110, savefig_dpi=250):
    """Set matplotlib's rcParams to the package style.

    ``dpi`` is the on-screen resolution, ``savefig_dpi`` that of saved
    figures, including the rasterised parts of saved PDF or SVG files.
    """
    mpl.rcParams.update({
        "axes.prop_cycle":   cycler(color=CYCLE),
        "figure.dpi":        dpi,
        "savefig.dpi":       savefig_dpi,
        "font.size":         10,
        "axes.titlesize":    11,
        "axes.grid":         True,
        "axes.xmargin":      0.0,
        "grid.alpha":        0.3,
        "grid.color":        "grey",
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "legend.fontsize":   9,
        "legend.fancybox":   False,
        "legend.edgecolor":  "black",
        "legend.framealpha": 1.0,
        "patch.linewidth":   0.6,
        "image.cmap":        "Greys",
    })
