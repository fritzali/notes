"""
magparsol/style.py
-------------------
Shared visual style: a small palette of named colours, legend framing and
rasterisation of dense artists, used by every plotting function.

Palette
~~~~~~~
Only named matplotlib colours are used: black and greys, plus

    steelblue · firebrick · olivedrab · goldenrod · rebeccapurple

Call :func:`use` once (e.g. at the top of a notebook) to apply the matching
rcParams: colour cycle, thin black legend frame, high-resolution output.
"""

import matplotlib as mpl
from cycler import cycler

# Semantic colours
B_COLOR     = "steelblue"       # magnetic field
E_COLOR     = "firebrick"       # electric field
ORBIT_COLOR = "black"           # particle trajectories over fields
EARTH_FACE  = "steelblue"
EARTH_ALPHA = 0.25
AXIS_COLOR  = "dimgrey"
REF_COLOR   = "darkgrey"        # reference curves / final spectra

CYCLE = ["steelblue", "firebrick", "olivedrab", "goldenrod", "rebeccapurple"]

# Artists with more points than this are rasterised inside vector output
RASTER_THRESHOLD = 4000

LEGEND_KW = dict(frameon=True, fancybox=False, edgecolor="black",
                 framealpha=0.92, borderpad=0.5)


def legend(ax, *args, **kwargs):
    """``ax.legend`` with the package frame: thin black border, square corners."""
    kw = {**LEGEND_KW, **kwargs}
    leg = ax.legend(*args, **kw)
    leg.get_frame().set_linewidth(0.6)
    return leg


def raster(n_points):
    """True if a line with ``n_points`` should be rasterised in vector output."""
    return n_points > RASTER_THRESHOLD


def colors(n):
    """``n`` palette colours (cycled)."""
    return [CYCLE[i % len(CYCLE)] for i in range(n)]


def use(dpi=110, savefig_dpi=200):
    """Apply the package style to matplotlib's rcParams."""
    mpl.rcParams.update({
        "axes.prop_cycle":   cycler(color=CYCLE),
        "figure.dpi":        dpi,
        "savefig.dpi":       savefig_dpi,
        "font.size":         10,
        "axes.titlesize":    11,
        "axes.grid":         True,
        "grid.alpha":        0.3,
        "grid.color":        "grey",
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "legend.fontsize":   9,
        "legend.fancybox":   False,
        "legend.edgecolor":  "black",
        "legend.framealpha": 0.92,
        "legend.linewidth":  0.6,
        "image.cmap":        "Greys",
    })
