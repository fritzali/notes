#!/usr/bin/env python

import os
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import matplotlib as mpl

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)


ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / 'PLUTO' / 'Simulation_Main'
OUT = ROOT / 'plots'

FORMAT = os.environ.get('PLOT_FORMAT', 'png')
DPI = int(os.environ.get('PLOT_DPI', '600'))
ONLY = [s for s in os.environ.get('PLOT_ONLY', '').split(',') if s]

COL = 82.5 / 25.4
TEXT = 170.0 / 25.4

preamble = r'''
    \usepackage{amsmath}
    \usepackage{mathtools}

    \usepackage{fontspec}
    \usepackage{unicode-math}
    \setmainfont{TeX Gyre Termes}
    \setmathfont{TeX Gyre Termes Math}

    \usepackage[
        locale = US,
        per-mode = power,
        exponent-mode = input,
        print-unity-mantissa = false
    ]{siunitx}

    \def\mathdefault#1{#1}
'''

mpl.rcParams.update({
    "backend": 'pgf',
    "pgf.texsystem": 'lualatex',
    "pgf.rcfonts": False,
    "pgf.preamble": preamble,
    "text.usetex": True,
    "font.family": 'serif',
    "font.size": 9,

    "axes.prop_cycle": "cycler('color', ['steelblue', 'olivedrab', 'goldenrod', 'firebrick', 'rebeccapurple'])",
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "axes.linewidth": 0.6,

    "figure.constrained_layout.use": True,
    "figure.constrained_layout.h_pad": 0.00,
    "figure.constrained_layout.w_pad": 0.00,

    "legend.fontsize": 8,
    "legend.frameon": False,
    "legend.handlelength": 1.5,
    "legend.columnspacing": 1.2,
    "legend.handletextpad": 0.5,
    "legend.borderaxespad": 0.4,

    "lines.linewidth": 1.0,

    "xtick.labelsize": 9,
    "xtick.minor.visible": True,
    "xtick.direction": 'in',
    "xtick.top": True,
    "ytick.labelsize": 9,
    "ytick.minor.visible": True,
    "ytick.direction": 'in',
    "ytick.right": True,

    "savefig.dpi": DPI,
})

import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm, SymLogNorm, TwoSlopeNorm, Normalize
from matplotlib.ticker import FixedLocator, FuncFormatter, MaxNLocator, NullFormatter
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.legend_handler import HandlerBase
from matplotlib.collections import LineCollection, PolyCollection, QuadMesh
from matplotlib.path import Path as MplPath

from scipy import signal
from scipy.ndimage import uniform_filter1d
from scipy.interpolate import CubicSpline

import statsmodels.api as sm

import pyPLUTO as pp


MAIN_PATH = DATA / 'BH_VISC_RAD_HD_TRC_COMP_FIX'
REF_PATH = DATA / 'BH_VISC_RAD_HD_TRC_COMP'

TRACER_RUNS = [('HD', DATA / 'BH_VISC_HD', DATA / 'BH_VISC_HD_TRC'),
               ('MHD', DATA / 'BH_VISC_MHD', DATA / 'BH_VISC_MHD_TRC'),
               ('RMHD', DATA / 'BH_VISC_RES_MHD', DATA / 'BH_VISC_RES_MHD_TRC')]
RAD_TRACER_RUN = DATA / 'BH_VISC_RAD_HD_TRC'

T_VALID = [0.0, 250.0, 500.0, 750.0]
T_EARLY = [100, 200, 300, 400, 500, 600]
FULL_AVG_SPAN = 1.0e4

W_WINDOW = 2500.0
BLOCK = 1000.0
REF_FRAC = 0.2
Z_CRIT = 2.0
TOL_EQ = 0.2

R_WAVE = (4.0, 25.0)
WELCH_SEG = 6
WELCH_TSEG = 2500.0
R_PHASE_REF = 7.0
FIT_FRAC = 0.15
P_FUND_MIN = 0.3

R_PROFILES = [8.0, 12.0, 16.0]
R_ANG = (2.5, 3.0, 4.0, 6.0, 10.0, 20.0)
PROFILE_STRIDE = 4
R_BUDGET = (2.5, 25.0)
R_MIN_HC = 4.0
HEMI = 2.0

CONST_G = 6.6726e-8
CONST_c = 2.99792458e10
CONST_amu = 1.66053886e-24
CONST_kB = 1.3806505e-16
CONST_sigma = 5.67051e-5
CONST_me = 9.1093826e-28
CONST_Msun = 1.989e33


def wanted(name):
    return not ONLY or any(name.startswith(key) for key in ONLY)


LEGEND_LOCS = ['upper right', 'upper left', 'lower left', 'lower right', 'center right', 'center left',
               'lower center', 'upper center', 'center']


def place_legend(ax, handles=None, labels=None, loc=None, **kw):
    if handles is None:
        handles, labels = ax.get_legend_handles_labels()
    elif labels is None:
        labels = [h.get_label() for h in handles]
    ax.__dict__.setdefault('_pending_legends', []).append((handles, labels, loc, kw))


def densify(xy, step=1.5):
    out = []
    good = np.all(np.isfinite(xy), axis=1)
    for a, b, ga, gb in zip(xy[:-1], xy[1:], good[:-1], good[1:]):
        if ga and gb:
            n = int(np.ceil(np.hypot(*(b - a)) / step)) + 1
            out.append(a + np.linspace(0.0, 1.0, n)[:, None] * (b - a))
    if not out:
        return xy[good]
    return np.concatenate(out)


def data_geometry(ax):
    points, polys = [], []
    for a in ax.figure.axes:
        if a.bbox.bounds != ax.bbox.bounds:
            continue
        for line in a.lines:
            if not line.get_visible():
                continue
            xy = line.get_transform().transform(line.get_path().vertices)
            if str(line.get_linestyle()).lower() in ('none', '', ' '):
                points.append(xy[np.all(np.isfinite(xy), axis=1)])
            else:
                points.append(densify(xy))
        for c in a.collections:
            if isinstance(c, QuadMesh) or not c.get_visible():
                continue
            if isinstance(c, PolyCollection):
                polys += [c.get_transform().transform_path(path) for path in c.get_paths()]
            elif isinstance(c, LineCollection):
                points += [densify(c.get_transform().transform(seg)) for seg in c.get_segments() if len(seg)]
            elif len(c.get_offsets()):
                points.append(c.get_offset_transform().transform(c.get_offsets()))
    points = [p for p in points if len(p)]
    return (np.concatenate(points) if points else np.zeros((0, 2))), polys


def legend_rect(ax, loc, w, h, pad):
    A = ax.bbox
    x = {'left': A.x0 + pad, 'right': A.x1 - pad - w, 'center': A.x0 + 0.5 * (A.width - w)}
    y = {'lower': A.y0 + pad, 'upper': A.y1 - pad - h, 'center': A.y0 + 0.5 * (A.height - h)}
    parts = loc.split()
    if len(parts) == 1:
        return x['center'], y['center'], w, h
    return x[parts[1]], y[parts[0]], w, h


def overlap_score(rect, points, polys, others):
    x0, y0, w, h = rect
    x1, y1 = x0 + w, y0 + h
    inside = (points[:, 0] > x0 - 1) & (points[:, 0] < x1 + 1) & (points[:, 1] > y0 - 1) & (points[:, 1] < y1 + 1)
    score = float(inside.sum())
    if polys:
        gx, gy = np.meshgrid(np.linspace(x0, x1, 10), np.linspace(y0, y1, 6))
        grid = np.column_stack([gx.ravel(), gy.ravel()])
        hits = np.zeros(len(grid), dtype=bool)
        for path in polys:
            hits |= path.contains_points(grid)
        score += 2.0 * hits.sum()
    for o in others:
        if not (x1 < o.x0 or x0 > o.x1 or y1 < o.y0 or y0 > o.y1):
            score += 1e6
    return score


def expand_ylim(ax, frac=0.15):
    lo, hi = ax.get_ylim()
    if ax.get_yscale() == 'log':
        ax.set_ylim(lo, hi * (hi / lo)**frac)
    else:
        ax.set_ylim(lo, hi + frac * (hi - lo))


def best_legend(ax, handles, labels, loc, kw, others, name):
    fig = ax.figure
    order = ([loc] if loc else []) + [c for c in LEGEND_LOCS if c != loc]
    for attempt in range(6):
        fig.draw_without_rendering()
        renderer = fig._get_renderer()
        probe = ax.legend(handles, labels, loc=order[0], **kw)
        bb = probe.get_window_extent(renderer)
        pad = probe.borderaxespad * probe._fontsize * fig.dpi / 72.0
        probe.remove()

        points, polys = data_geometry(ax)
        scores = [overlap_score(legend_rect(ax, c, bb.width, bb.height, pad), points, polys, others) for c in order]
        k = int(np.argmin(scores))
        if scores[k] == 0 or attempt == 5:
            if scores[k] > 0:
                print(f'  legend in {name} still covers data, score {scores[k]:.0f}')
            leg = ax.legend(handles, labels, loc=order[k], **kw)
            x0, y0, w, h = legend_rect(ax, order[k], bb.width, bb.height, pad)
            return leg, mpl.transforms.Bbox.from_bounds(x0, y0, w, h)
        expand_ylim(ax)


def finalize_legends(fig, name):
    for ax in fig.axes:
        pending = ax.__dict__.pop('_pending_legends', [])
        placed = []
        for k, (handles, labels, loc, kw) in enumerate(pending):
            leg, box = best_legend(ax, handles, labels, loc, kw, placed, name)
            placed.append(box)
            if k < len(pending) - 1:
                ax.add_artist(leg)


def save(fig, name):
    OUT.mkdir(exist_ok=True)
    finalize_legends(fig, name)
    fig.savefig(OUT / f'{name}.{FORMAT}', dpi=DPI)
    plt.close(fig)
    print(f'wrote {name}.{FORMAT}')


def plain_log(ax, axis='x', mantissas=(1, 2, 3, 4, 6)):
    a = ax.xaxis if axis == 'x' else ax.yaxis
    lo, hi = sorted(ax.get_xlim() if axis == 'x' else ax.get_ylim())
    ticks = [m * 10.0**k for k in range(int(np.floor(np.log10(lo))) - 1, int(np.ceil(np.log10(hi))) + 1) for m in mantissas]
    a.set_major_locator(FixedLocator([t for t in ticks if lo <= t <= hi]))
    a.set_major_formatter(FuncFormatter(lambda v, pos: f'${v:g}$'))
    a.set_minor_formatter(NullFormatter())


def read_definitions(path):
    txt = (path / 'definitions.h').read_text()
    rho0 = re.search(r'#define\s+UNIT_DENSITY\s+([0-9.eE+-]+)', txt)
    mbh = re.search(r'#define\s+BH_MASS\s+\(\s*([0-9.eE+-]+)\s*\*\s*([0-9.eE+-]+)\s*\)', txt)
    rho0 = float(rho0.group(1)) if rho0 else 5.8e-10
    mbh = float(mbh.group(1)) * float(mbh.group(2)) if mbh else 1.0e8 * CONST_Msun
    return rho0, mbh


def read_parameters(path):
    params, inside = {}, False
    for line in (path / 'pluto.ini').read_text().splitlines():
        s = line.strip()
        if s.startswith('['):
            inside = (s == '[Parameters]')
            continue
        if inside and s and not s.startswith('#'):
            key, val = s.split()[:2]
            params[key] = float(val)
    return params


RHO0, M_BH = read_definitions(MAIN_PATH)
PAR = read_parameters(MAIN_PATH)

GAMMA = 5.0 / 3.0
X_H, Z_M = 0.70, 0.02
MU = 1.0 / (2.0 * X_H + 0.75 * (1.0 - X_H - Z_M) + 0.5 * Z_M)

ALPHA = PAR['BETAV']
EPS = PAR['EPS']
DFLOOR = PAR['DFLOOR']

R_G = CONST_G * M_BH / CONST_c**2
T_G = R_G / CONST_c
P0 = RHO0 * CONST_c**2
MDOT0 = RHO0 * CONST_c * R_G**2
L0 = RHO0 * CONST_c**3 * R_G**2

K_T = MU * CONST_amu * CONST_c**2 / CONST_kB
A_RAD = 4.0 * CONST_sigma / CONST_c / P0

KAPPA_ES = 0.2 * (1.0 + X_H)
ETA = 1.0 / 16.0
L_EDD = 4.0 * np.pi * CONST_G * M_BH * CONST_c / KAPPA_ES
MDOT_EDD = L_EDD / (ETA * CONST_c**2)
SCALE = MDOT0 / MDOT_EDD

KAPPA_CODE = RHO0 * R_G
K_FFBF = 4.34e25 * Z_M * (1.0 + X_H) / 10.0 + 3.68e22 * (1.0 - Z_M) * (1.0 + X_H)
T_COMPTON_MAX = 1.0e9


def omega_pw(R):
    return 1.0 / ((R - 2.0) * np.sqrt(R))


def omega_newton(R):
    return R**-1.5


def nu_kin(R):
    return 2.0 / 3.0 * ALPHA * EPS**2 * np.sqrt(R)


def cs0(R):
    return EPS / np.sqrt(R)


def t_orbit(R):
    return 2.0 * np.pi / omega_pw(R)


def t_thermal(R):
    return 1.0 / (ALPHA * omega_newton(R))


def t_viscous(R):
    return R**2 / nu_kin(R)


def t_sound(R):
    return R / cs0(R)


def r_viscous(t):
    return (2.0 * ALPHA * EPS**2 * np.asarray(t, dtype=float) / 3.0)**(2.0 / 3.0)


def read_grid(path):
    lines = [l for l in (path / 'grid.out').read_text().splitlines(True) if not l.startswith('#')]

    edges, k = [], 0
    while k < len(lines) and lines[k].strip():
        n = int(lines[k].split()[0])
        k += 1
        block = np.array([[float(v) for v in lines[k + m].split()[1:3]] for m in range(n)])
        edges.append(np.append(block[:, 0], block[-1, 1]))
        k += n

    return edges


def read_outputs(path):
    entries = {}
    for line in (path / 'dbl.out').read_text().splitlines():
        s = line.split()
        if len(s) > 1:
            entries[int(s[0])] = float(s[1])

    outs = np.array(sorted(entries))
    return outs, np.array([entries[n] for n in outs])


def make_geometry(path):
    R_E, TH_E = read_grid(path)[:2]
    R_C = 0.5 * (R_E[1:] + R_E[:-1])
    TH_C = 0.5 * (TH_E[1:] + TH_E[:-1])
    DOM = 2.0 * np.pi * (np.cos(TH_E[:-1]) - np.cos(TH_E[1:]))
    return dict(R_E=R_E, TH_E=TH_E, R_C=R_C, TH_C=TH_C, NR=R_C.size, NTH=TH_C.size, DTH=np.diff(TH_E), DOM=DOM,
                COT=np.cos(TH_C) / np.sin(TH_C))


def load(n, path, geo, names=('rho', 'vx1', 'vx2', 'vx3', 'prs', 'enr', 'fr1', 'fr2', 'diskfrac', 'nu', 'kappa_abs', 'kappa_scat')):
    D = pp.Load(nout=int(n), path=str(path) + '/')
    shape = (geo['NR'], geo['NTH'])

    out = {}
    for v in names:
        if hasattr(D, v):
            a = np.asarray(getattr(D, v), dtype=float)
            assert a.shape == shape, f"{v} has shape {a.shape}, expected {shape}"
            out[v] = a

    return out


def nearest_out(path, t):
    outs, times = read_outputs(path)
    k = int(np.argmin(np.abs(times - t)))
    return outs[k], times[k]


def gas_temperature(F):
    return K_T * F['prs'] / F['rho']


def rad_temperature(F):
    return (np.clip(F['enr'], 0.0, None) / A_RAD)**0.25


def optical_depth_from_pole(F, g):
    dtau = F['rho'] * (F['kappa_abs'] + F['kappa_scat']) * g['R_C'][:, None] * g['DTH'][None, :]
    return np.cumsum(dtau, axis=1)


def read_cache(path):
    cache = path / 'analysis_cache.npz'
    if not cache.exists():
        sys.exit(f"{cache} is missing, run Full_Analysis.ipynb once to build it")
    with np.load(cache) as Z:
        return {key: Z[key] for key in Z.files}


def frames(span, dt):
    return max(int(round(span / dt)), 1)


def cadences(t):
    dts = np.diff(t)
    dts = dts[dts > 0]
    rounded = np.array([float(f"{x:.3g}") for x in dts])
    values, counts = np.unique(rounded, return_counts=True)
    frequent = values[counts >= 10]
    return frequent if frequent.size else np.array([np.median(dts)])


def long_indices(t, dt_long):
    phase = t / dt_long
    on_grid = np.abs(phase - np.round(phase)) < 0.02
    idx = np.nonzero(on_grid)[0]
    _, first = np.unique(np.round(phase[idx]).astype(int), return_index=True)
    return idx[first]


def dense_indices(t, dt_long):
    fine = cadences(t).min()
    if fine > 0.5 * dt_long:
        return None
    jumps = np.nonzero(np.diff(t) > 1.5 * fine)[0]
    t_dense = t[jumps[-1] + 1] if jumps.size else t[0]
    idx = np.nonzero(t >= t_dense)[0]
    return idx if idx.size > 20 else None


def slice_cache(Cx, idx):
    n = Cx['t'].size
    return {key: (val[idx] if (np.ndim(val) > 0 and val.shape[0] == n and key != 'r') else val) for key, val in Cx.items()}


def mean_se(X, dt):
    X = np.asarray(X, dtype=float)
    n = X.shape[0]
    size = min(frames(BLOCK, dt), max(n // 4, 1))
    nb = n // size

    blocks = X[:nb * size].reshape((nb, size) + X.shape[1:]).mean(axis=1)
    return X.mean(axis=0), blocks.std(axis=0, ddof=1) / np.sqrt(nb)


def window_means(t, X, dt):
    X = np.asarray(X, dtype=float)
    n_w = min(frames(W_WINDOW, dt), X.shape[0])
    cs = np.cumsum(np.concatenate([np.zeros((1,) + X.shape[1:]), X]), axis=0)
    means = (cs[n_w:] - cs[:-n_w]) / n_w
    centres = 0.5 * (t[n_w - 1:] + t[:t.size - n_w + 1])
    return centres, means, n_w


def relaxation(t, X, dt):
    X = np.asarray(X, dtype=float)
    n_ref = max(int(round(REF_FRAC * t.size)), 8)
    tc, m, n_w = window_means(t, X, dt)
    ref, se = mean_se(X[-n_ref:], dt)
    sigma = np.maximum(se * np.sqrt(1.0 + n_ref / n_w), 1e-300)
    D = np.abs(m - ref) / sigma

    t_ref = t[-n_ref]
    judged = tc + 0.5 * W_WINDOW <= t_ref
    if not judged.any():
        return np.full(X.shape[1:], np.nan), tc, D, t_ref

    ok = D[judged] < Z_CRIT
    tail_ok = np.flip(np.logical_and.accumulate(np.flip(ok, axis=0), axis=0), axis=0)
    first = np.argmax(tail_ok, axis=0)
    return np.where(tail_ok.any(axis=0), tc[judged][first], np.nan), tc, D, t_ref


GLOBAL_TEX = {'M_d': r'$M_d$', 'E_rad': r'$E_\mathrm{rad}$', 'mdot_in': r'$\dot M_\mathrm{acc}(r_\mathrm{in})$', 'L_out': r'$L(r_\mathrm{out})$'}


def global_series(Cx):
    return {'M_d': Cx['M_d'], 'E_rad': Cx['E_rad'], 'mdot_in': -Cx['mdot'][:, 0], 'L_out': Cx['L_out']}


def inflow_equilibrium(run):
    t, r, dt = run['t'], run['r'], run['dt']
    tc, wm, _ = window_means(t, -run['C']['mdot'], dt)

    inner = (r > r[0]) & (r < 6.0)
    ref = wm[:, inner].mean(axis=1)
    smooth = uniform_filter1d(wm, 5, axis=1, mode='nearest')
    ok = np.abs(smooth / ref[:, None] - 1.0) < TOL_EQ
    ok[:, r < 6.0] = True
    ok[:, -2:] = True

    flat = np.logical_and.accumulate(ok, axis=1)
    r_eq = np.array([r[np.nonzero(row)[0][-1]] if row.any() else np.nan for row in flat])

    steady = tc - 0.5 * W_WINDOW >= run['t_steady']
    return np.nanpercentile(r_eq[steady], [16, 50, 84]) if steady.any() else np.full(3, np.nan)


def make_run(path, label, color, lw):
    Cx = read_cache(path)
    t_all = Cx['t']
    dt_long = cadences(t_all).max()
    il = long_indices(t_all, dt_long)
    idn = dense_indices(t_all, dt_long)

    run = dict(path=path, label=label, color=color, lw=lw, geo=make_geometry(path), r=Cx['r'], all=Cx,
               C=slice_cache(Cx, il), dt=dt_long)
    run['t'] = run['C']['t']

    times, curves = {}, {}
    for key, X in global_series(run['C']).items():
        tr, tc, D, t_ref = relaxation(run['t'], X, run['dt'])
        times[key] = float(tr)
        curves[key] = (tc, D)
    run['t_ref'] = t_ref
    run['t_relax'], run['D_curves'] = times, curves
    run['t_steady'] = max(v if np.isfinite(v) else t_ref for v in times.values())
    run['sel'] = run['t'] >= run['t_steady']

    if idn is not None:
        D = slice_cache(Cx, idn)
        run['S'] = slice_cache(D, np.nonzero(D['t'] >= run['t_steady'])[0])
        run['dt_s'] = float(np.median(np.diff(t_all[idn])))
    else:
        run['S'] = slice_cache(run['C'], np.nonzero(run['sel'])[0])
        run['dt_s'] = run['dt']

    run['R_eq_stats'] = inflow_equilibrium(run)

    s, dt = run['sel'], run['dt']
    run['mdot_in_mean'], run['mdot_in_se'] = mean_se(-run['C']['mdot'][s, 0], dt)
    run['sig_mean'], run['sig_se'] = mean_se(run['C']['sigma'][s], dt)
    run['hr_mean'], run['hr_se'] = mean_se(run['C']['hr'][s], dt)
    run['cells_H'] = run['hr_mean'] / run['geo']['DTH'][-1]
    run['delta'] = run['S']['sigma'] / run['S']['sigma'].mean(axis=0)[None, :] - 1.0

    return run


class TableHeader:
    def __init__(self, text, color):
        self.text = text
        self.color = color


class TableHeaderHandler(HandlerBase):
    def legend_artist(self, legend, handle, fontsize, handlebox):
        x = 0.5 * handlebox.width - handlebox.xdescent
        y = 0.5 * handlebox.height - handlebox.ydescent
        text = Text(x, y, handle.text, color=handle.color, fontsize=fontsize, ha='center', va='center')
        handlebox.add_artist(text)
        return text


def run_legend(ax, runs, color=False, extra=(), **kw):
    handles = []
    for run in runs:
        if color:
            handles.append(Line2D([], [], color=run['color'], label=run['label']))
        else:
            handles.append(Line2D([], [], color='k', lw=run['lw'], label=run['label']))
    handles += [Line2D([], [], label=label, **{'color': 'k', **style}) for label, style in extra]
    place_legend(ax, handles, **kw)


def run_table(ax, rows, columns, **kw):
    handles, labels = [], []
    for j, (text, color, over) in enumerate(columns):
        last = j == len(columns) - 1
        handles += [TableHeader(text, color)] + [Line2D([], [], **{'color': 'k', **style, **over}) for _, style in rows]
        labels += [''] + [label if last else '' for label, _ in rows]

    opts = dict(ncol=len(columns), columnspacing=0.3, handlelength=1.5, handletextpad=0.4, handleheight=1.8,
                handler_map={TableHeader: TableHeaderHandler()})
    place_legend(ax, handles, labels, **{**opts, **kw})


MAP_FIELDS = {
    'rho':   dict(norm=LogNorm(vmin=DFLOOR, vmax=0.2), cmap='magma', label=r'$\rho\ /\ \rho_0$'),
    'f':     dict(norm=Normalize(vmin=0.0, vmax=1.0), cmap='coolwarm', label=r'$f$'),
    'tr':    dict(norm=Normalize(vmin=0.0, vmax=1.0), cmap='coolwarm', label=r'$f_\mathrm{tr}$'),
    'diff':  dict(norm=Normalize(vmin=-1.0, vmax=1.0), cmap='RdBu_r', label=r'$\Delta f$'),
    'Tg':    dict(norm=LogNorm(vmin=1e4, vmax=1e12), cmap='inferno', label=r'$T_g\ /\ \mathrm{K}$'),
    'ratio': dict(norm=LogNorm(vmin=1e-6, vmax=1e6), cmap='RdBu_r', label=r'$T_g / T_r$'),
    'vr':    dict(norm=SymLogNorm(linthresh=1e-3, vmin=-1.0, vmax=1.0), cmap='RdBu_r', label=r'$v_r\ /\ c$'),
    'vphi':  dict(norm=Normalize(vmin=0.0, vmax=0.8), cmap='viridis', label=r'$v_\phi\ /\ c$'),
}

MAP_LEFT, MAP_BOTTOM, MAP_TOP = 0.40, 0.42, 0.06
MAP_GAP, GROUP_GAP = 0.05, 0.16
CB_PAD, CB_WIDTH, CB_RIGHT = 0.06, 0.08, 0.72


def map_field(F, key):
    if key == 'rho':
        return F['rho']
    if key == 'f':
        return F['diskfrac']
    if key == 'tr':
        return F['tr1']
    if key == 'diff':
        return F['diskfrac'] - F['tr1']
    if key == 'Tg':
        return gas_temperature(F)
    if key == 'ratio':
        return gas_temperature(F) / np.maximum(rad_temperature(F), 1.0)
    if key == 'vr':
        return F['vx1']
    if key == 'vphi':
        return F['vx3']
    raise KeyError(key)


def place_map(fig, left, right, bottom, top):
    spec = fig._map_spec
    groups, ncol, width, aspect = spec['groups'], spec['ncol'], spec['width'], spec['aspect']
    axes, caxes = spec['axes'], spec['caxes']
    w = (width - left - right - (ncol - 1) * MAP_GAP) / ncol
    h = w * aspect
    nrow = sum(groups)
    height = bottom + top + nrow * h + (nrow - len(groups)) * MAP_GAP + (len(groups) - 1) * GROUP_GAP

    fig.set_size_inches(width, height)
    y, i = height - top, 0
    for g, n in enumerate(groups):
        group_top = y
        for k in range(n):
            for j in range(ncol):
                axes[i, j].set_position([(left + j * (w + MAP_GAP)) / width, (y - h) / height, w / width, h / height])
            y -= h
            i += 1
            if k < n - 1:
                y -= MAP_GAP
        cb_left = left + ncol * w + (ncol - 1) * MAP_GAP + CB_PAD
        if caxes[g] is not None:
            caxes[g].set_position([cb_left / width, y / height, CB_WIDTH / width, (group_top - y) / height])
        y -= GROUP_GAP
    spec['margins'] = (left, right, bottom, top)


def tighten_map(fig, pad=0.01, passes=2):
    for _ in range(passes):
        fig.draw_without_rendering()
        bb = fig.get_tightbbox(fig._get_renderer())
        width, height = fig.get_size_inches()
        left, right, bottom, top = fig._map_spec['margins']
        slack = (bb.x0, width - bb.x1, bb.y0, height - bb.y1)
        place_map(fig, *(m - s + pad for m, s in zip((left, right, bottom, top), slack)))


def map_layout(groups, ncol, width, xlim, ylim, cbar=True):
    nrow = sum(groups)
    fig = plt.figure(figsize=[width, width], layout='none')
    axes = np.empty((nrow, ncol), dtype=object)
    for i in range(nrow):
        for j in range(ncol):
            share = axes[0, 0] if i or j else None
            axes[i, j] = fig.add_axes([0, 0, 1, 1], sharex=share, sharey=share)
    caxes = [fig.add_axes([0, 0, 1, 1]) if cbar else None for _ in groups]
    fig._map_spec = dict(groups=groups, ncol=ncol, width=width, aspect=(ylim[1] - ylim[0]) / (xlim[1] - xlim[0]),
                         axes=axes, caxes=caxes)
    place_map(fig, MAP_LEFT, CB_PAD + CB_WIDTH + CB_RIGHT if cbar else 0.1, MAP_BOTTOM, MAP_TOP)

    for i in range(nrow):
        for j in range(ncol):
            ax = axes[i, j]
            ax.tick_params(which='both', direction='out', color='k', top=False, right=False, length=2.5)
            ax.tick_params(which='minor', length=1.5)
            ax.tick_params(labelleft=(j == 0), labelbottom=(i == nrow - 1))
            if j == 0:
                ax.set_ylabel(r'$z\ /\ r_g$')
            if i == nrow - 1:
                ax.set_xlabel(r'$R\ /\ r_g$')
    axes[0, 0].xaxis.set_major_locator(MaxNLocator(5, prune='upper' if ncol > 1 else None))
    axes[0, 0].yaxis.set_major_locator(MaxNLocator(5, prune='upper' if nrow > 1 else None))
    return fig, axes, caxes


def draw_map(ax, field, key, geo, xlim, ylim):
    spec = MAP_FIELDS[key]
    Xe = geo['R_E'][:, None] * np.sin(geo['TH_E'])[None, :]
    Ze = geo['R_E'][:, None] * np.cos(geo['TH_E'])[None, :]
    im = ax.pcolormesh(Xe, Ze, field, norm=spec['norm'], cmap=spec['cmap'], shading='flat', rasterized=True)
    ax.add_patch(plt.Circle((0, 0), geo['R_E'][0], color='k', zorder=3))
    ax.set_facecolor('0.85')
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect('equal')

    Xc = geo['R_C'][:, None] * np.sin(geo['TH_C'])[None, :]
    Zc = geo['R_C'][:, None] * np.cos(geo['TH_C'])[None, :]
    ax._map_cells = ((Xc.ravel() - xlim[0]) / (xlim[1] - xlim[0]), (Zc.ravel() - ylim[0]) / (ylim[1] - ylim[0]),
                     im.to_rgba(np.ma.masked_invalid(field).ravel()))
    return im


def contrast_color(ax, x0, x1, y0, y1):
    fx, fy, rgba = ax._map_cells
    m = (fx >= x0) & (fx <= x1) & (fy >= y0) & (fy <= y1)
    if not m.any():
        lum = 0.85
    else:
        c = rgba[m, :3]
        lum = np.mean(0.2126 * c[:, 0] + 0.7152 * c[:, 1] + 0.0722 * c[:, 2])
    return 'k' if lum > 0.45 else 'w'


def time_label(ax, t):
    color = contrast_color(ax, 0.55, 1.0, 0.75, 1.0)
    ax.text(0.96, 0.93, rf'$t = \num{{{t:.0f}}}\,t_g$', color=color, ha='right', va='top', transform=ax.transAxes, fontsize=8)


def row_tag(ax, text):
    color = contrast_color(ax, 0.0, 0.35, 0.75, 1.0)
    ax.text(0.04, 0.93, text, color=color, ha='left', va='top', transform=ax.transAxes, fontsize=8)


def map_colorbar(fig, im, cax, key):
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(MAP_FIELDS[key]['label'])
    cb.ax.tick_params(which='both', direction='out', right=True, left=False, color='k')
    if key in ('f', 'tr', 'diff'):
        cb.set_ticks([-1.0, -0.5, 0.0, 0.5, 1.0] if key == 'diff' else [0.0, 0.5, 1.0])
        cb.ax.minorticks_off()
    return cb


def map_grid(groups, width, name, xlim=(0, 20), ylim=(0, 10), cbar=True, yaxis=True):
    fig = build_map_grid(groups, width, xlim, ylim, cbar, yaxis)
    tighten_map(fig)
    height = fig.get_size_inches()[1]
    save(fig, name)
    return height


def build_map_grid(groups, width, xlim=(0, 20), ylim=(0, 10), cbar=True, yaxis=True):
    """groups: list of (colorbar key, rows), each row a list of cells (field, key, geo, t, tag[, extra])"""
    ncol = len(groups[0][1][0])
    fig, axes, caxes = map_layout([len(rows) for _, rows in groups], ncol, width, xlim, ylim, cbar=cbar)
    if not yaxis:
        for ax in axes[:, 0]:
            ax.set_ylabel('')
            ax.tick_params(labelleft=False)
    i = 0
    for (cbar_key, rows), cax in zip(groups, caxes):
        im = None
        for row in rows:
            for j, cell in enumerate(row):
                field, key, geo, t, tag = cell[:5]
                im = draw_map(axes[i, j], field, key, geo, xlim, ylim)
                if len(cell) > 5:
                    cell[5](axes[i, j])
                if t is not None:
                    time_label(axes[i, j], t)
                if tag is not None:
                    row_tag(axes[i, j], tag)
            i += 1
        if cax is not None:
            map_colorbar(fig, im, cax, cbar_key)
    return fig


def map_pair(left_fig, right_fig, names, height=None):
    """Tighten two single panel map figures of equal width and align them for placement side by side:
    equal panel size and height, extra horizontal space on the outer sides. With a target height, the
    panels shrink at fixed aspect until both figures have exactly that height."""
    figs = (left_fig, right_fig)
    for fig in figs:
        tighten_map(fig)
    margins = [fig._map_spec['margins'] for fig in figs]
    span = max(l + r for l, r, _, _ in margins)
    bottom = max(m[2] for m in margins)
    top = max(m[3] for m in margins)
    if height is not None:
        spec = left_fig._map_spec
        panel = (height - bottom - top) / spec['aspect']
        span = max(span, spec['width'] - panel)
    (l0, r0, _, _), (l1, r1, _, _) = margins
    place_map(left_fig, l0 + span - (l0 + r0), r0, bottom, top)
    place_map(right_fig, l1, r1 + span - (l1 + r1), bottom, top)
    for fig, name in zip(figs, names):
        save(fig, name)


SNAP_CACHE = {}


def snapshot(path, t_req, names=('rho', 'vx1', 'vx2', 'vx3', 'prs', 'enr', 'diskfrac', 'tr1')):
    n, t = nearest_out(path, t_req)
    key = (str(path), int(n))
    if key not in SNAP_CACHE:
        SNAP_CACHE[key] = (t, make_geometry(path), load(n, path, make_geometry(path), names=names))
    return SNAP_CACHE[key]


def time_row(path, key, times, label_times=True, tag=None):
    row = []
    for j, t_req in enumerate(times):
        t, geo, F = snapshot(path, t_req)
        row.append((map_field(F, key), key, geo, t if label_times else None, tag if j == 0 else None))
    return row


def agreement_series(path):
    geo = make_geometry(path)
    outs, times = read_outputs(path)
    w = geo['DOM'][None, :] * ((geo['R_E'][1:]**3 - geo['R_E'][:-1]**3) / 3.0)[:, None]
    vals = []
    for n in outs:
        F = load(n, path, geo, names=('rho', 'diskfrac', 'tr1'))
        if 'tr1' not in F or 'diskfrac' not in F:
            return None
        diff = np.abs(F['diskfrac'] - F['tr1'])
        mass = F['rho'] * w
        vals.append((1.0 - np.sum(diff * w) / np.sum(w), 1.0 - np.sum(diff * mass) / np.sum(mass)))
    return times, np.array(vals)


def fig_validation():
    diff_rows = []
    for label, tracer_path, df_path in TRACER_RUNS:
        if not (tracer_path / 'dbl.out').exists() or not (df_path / 'dbl.out').exists():
            print(f'skipping validation {label}, outputs missing')
            continue

        tag = label.lower()
        map_grid([('rho', [time_row(tracer_path, 'rho', T_VALID), time_row(df_path, 'rho', T_VALID, False)]),
                  ('f', [time_row(df_path, 'tr', T_VALID, False), time_row(df_path, 'f', T_VALID, False)])],
                 TEXT, f'valid_{tag}_text')
        diff_rows.append(time_row(df_path, 'diff', T_VALID, not diff_rows, rf'$\mathrm{{{label}}}$'))

        for path, key, kind in [(tracer_path, 'rho', 'tracerrun'), (df_path, 'rho', 'dfrun'), (df_path, 'tr', 'tracer'),
                                (df_path, 'f', 'f'), (df_path, 'diff', 'diff')]:
            map_grid([(key, [time_row(path, key, [T_VALID[-1]])])], COL, f'valid_{tag}_{kind}_col')

    if diff_rows:
        map_grid([('diff', diff_rows)], TEXT, 'valid_diff_text')

    if (RAD_TRACER_RUN / 'dbl.out').exists():
        map_grid([('rho', [time_row(RAD_TRACER_RUN, 'rho', T_VALID)]),
                  ('f', [time_row(RAD_TRACER_RUN, 'tr', T_VALID, False), time_row(RAD_TRACER_RUN, 'f', T_VALID, False)])],
                 TEXT, 'valid_rad_text')
        for key, kind in [('tr', 'tracer'), ('f', 'f')]:
            map_grid([(key, [time_row(RAD_TRACER_RUN, key, [T_VALID[-1]])])], COL, f'valid_rad_{kind}_col')

    fig, ax = plt.subplots(figsize=[COL, 0.62 * COL])
    shown = False
    for k, (label, _, df_path) in enumerate(TRACER_RUNS):
        if not (df_path / 'dbl.out').exists():
            continue
        res = agreement_series(df_path)
        if res is None:
            continue
        t, a = res
        ax.plot(t, 100.0 * a[:, 1], color=f'C{k}', label=rf'$\mathrm{{{label}}}$')
        ax.plot(t, 100.0 * a[:, 0], color=f'C{k}', lw=0.6, ls='--')
        shown = True
    if shown:
        ax.set_xlabel(r'$t\ /\ t_g$')
        ax.set_ylabel(r'$1 - \langle|\Delta f|\rangle\ /\ \%$')
        place_legend(ax, ncol=3)
        save(fig, 'valid_agreement_col')
    else:
        plt.close(fig)


def fig_timescales():
    R = np.geomspace(2.5, 30, 300)
    fig, ax = plt.subplots(figsize=[COL, 0.62 * COL])

    ax.axvline(6.0, color='gray', lw=0.6, ls=':')
    ax.plot(R, t_orbit(R), label=r'$t_\mathrm{orb}$')
    ax.plot(R, t_thermal(R), label=r'$t_\mathrm{th}$')
    ax.plot(R, t_sound(R), label=r'$t_{c_s}$')
    ax.plot(R, t_viscous(R), label=r'$t_\nu$')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(2.5, 30)
    plain_log(ax, 'x')

    ax.set_xlabel(r'$R\ /\ r_g$')
    ax.set_ylabel(r'$t\ /\ t_g$')
    place_legend(ax, ncol=2)

    save(fig, 'timescales_col')


def fig_full_domain(run):
    path, geo = run['path'], run['geo']
    lim = (0, geo['R_E'][-1])
    n0 = run['all']['outs'][int(np.argmin(np.abs(run['all']['t'])))]
    F0 = load(n0, path, geo)

    sel = run['t'] >= run['t'][-1] - FULL_AVG_SPAN
    acc = {}
    for n in run['C']['outs'][sel]:
        F = load(n, path, geo)
        for key, val in [('rho', F['rho']), ('Tg', gas_temperature(F)), ('Tr', rad_temperature(F)), ('f', F['diskfrac']),
                         ('vr', F['vx1']), ('vphi', F['vx3']), ('tau', optical_depth_from_pole(F, geo))]:
            acc[key] = acc.get(key, 0.0) + val
    mean = {key: val / sel.sum() for key, val in acc.items()}
    mean['ratio'] = mean['Tg'] / np.maximum(mean['Tr'], 1.0)
    init = {key: map_field(F0, key) for key in ('rho', 'f', 'Tg', 'ratio', 'vr', 'vphi')}
    init['tau'] = optical_depth_from_pole(F0, geo)

    def cells(key):
        return [(init[key], key, geo, None, None), (mean[key], key, geo, None, None)]

    col_height = {}
    for key in ('rho', 'f', 'Tg', 'ratio', 'vr', 'vphi'):
        map_grid([(key, [cells(key)])], TEXT, f'domain_{key}_text', lim, lim)
        for k, label in enumerate(('init', 'mean')):
            col_height[key, label] = map_grid([(key, [[cells(key)[k]]])], COL, f'domain_{key}_{label}_col', lim, lim)

    map_pair(build_map_grid([('rho', [[cells('rho')[0]]])], COL, lim, lim, cbar=False),
             build_map_grid([('rho', [[cells('rho')[1]]])], COL, lim, lim, yaxis=False),
             ('domain_rho_init_alt_col', 'domain_rho_mean_alt_col'), height=col_height['rho', 'init'])

    for name, (a, b) in [('thermal', ('Tg', 'ratio')), ('velocity', ('vr', 'vphi'))]:
        map_grid([(a, [cells(a)]), (b, [cells(b)])], TEXT, f'domain_{name}_text', lim, lim)

    states = [(init['f'], init['tau']), (mean['f'], mean['tau'])]
    outline_figure(states, geo, lim, TEXT, 'domain_outline_text')
    for state, label in zip(states, ('init', 'mean')):
        outline_figure([state], geo, lim, COL, f'domain_outline_{label}_col')

    return mean


def outline_map(ax, f, tau, geo, lim):
    Xc = geo['R_C'][:, None] * np.sin(geo['TH_C'])[None, :]
    Zc = geo['R_C'][:, None] * np.cos(geo['TH_C'])[None, :]
    ax.set_facecolor('w')
    ax.contour(Xc, Zc, f, levels=[0.5], colors='k', linewidths=1.0)
    if np.nanmax(tau) > 1.0:
        ax.contour(Xc, Zc, tau, levels=[1.0], colors='k', linewidths=1.0, linestyles='--')
    th = np.linspace(0.0, 0.5 * np.pi, 200)
    ax.plot(geo['R_E'][-1] * np.sin(th), geo['R_E'][-1] * np.cos(th), 'k-', lw=1.0)
    ax.add_patch(plt.Circle((0, 0), geo['R_E'][0], color='k', zorder=3))
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_aspect('equal')


def outline_figure(states, geo, lim, width, name):
    fig, axes, _ = map_layout([1], len(states), width, lim, lim, cbar=False)
    for ax, (f, tau) in zip(axes[0], states):
        outline_map(ax, f, tau, geo, lim)
    tighten_map(fig)
    save(fig, name)


def fig_early(run, keys=('rho', 'f', 'Tg', 'ratio', 'vr', 'vphi')):
    path = run['path']
    half = len(T_EARLY) // 2

    def rows(key, labels=True):
        return [time_row(path, key, T_EARLY[:half], labels), time_row(path, key, T_EARLY[half:], labels)]

    for key in keys:
        map_grid([(key, rows(key))], TEXT, f'early_{key}_text')

    for name, (a, b) in [('thermal', ('Tg', 'ratio')), ('velocity', ('vr', 'vphi'))]:
        map_grid([(a, rows(a)), (b, rows(b, False))], TEXT, f'early_{name}_text')


def fig_relaxation(runs):
    for width, aspect, suffix in [(TEXT, 0.36, 'text'), (COL, 0.62, 'col')]:
        fig, ax = plt.subplots(figsize=[width, aspect * width])
        tmax = max(run['t'][-1] for run in runs)

        for run in runs:
            Cx, t = run['C'], run['t']
            for k, (key, name) in enumerate([('M', r'$M$'), ('M_d', r'$M_d$'), ('M_c', r'$M_c$')]):
                ax.plot(t, Cx[key], color=f'C{k}', lw=run['lw'], label=name if run is runs[0] else None)

        ax.set_yscale('log')
        ax.set_xscale('symlog', linthresh=11e3, linscale=1.0)
        ax.set_xlim(0, tmax)
        ax.xaxis.set_major_locator(FixedLocator([0, 1e4]))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: '$0$' if v == 0 else rf'$10^{{{int(np.floor(np.log10(v)))}}}$'))
        ax.xaxis.set_minor_locator(FixedLocator([k * 1e3 for k in range(1, 51) if k != 10]))
        ax.xaxis.set_minor_formatter(NullFormatter())

        ax.set_xlabel(r'$t\ /\ t_g$')
        ax.set_ylabel(r'$M\ /\ \rho_0 r_g^3$')

        place_legend(ax, ncol=3, loc='upper right')
        run_legend(ax, runs, ncol=2, loc='lower right')

        save(fig, f'relax_mass_{suffix}')

    fig, ax = plt.subplots(figsize=[COL, 0.62 * COL])
    run = runs[0]
    for k, (key, (tc, D)) in enumerate(run['D_curves'].items()):
        ax.plot(tc, D, color=f'C{k}', label=GLOBAL_TEX[key])
    ax.axhline(Z_CRIT, color='k', lw=0.6, ls='--', zorder=-10)
    ax.axvline(run['t_ref'] - 0.5 * W_WINDOW, color='gray', lw=0.6, ls=':', zorder=-10)
    ax.set_yscale('log')
    ax.set_xlabel(r'$t_c\ /\ t_g$')
    ax.set_ylabel(r'$D(t_c)$')
    place_legend(ax, ncol=2, loc='upper right')
    save(fig, 'relax_distance_col')


def fig_window(runs):
    fig, axs = plt.subplots(3, 1, figsize=[TEXT, 0.8 * TEXT], sharex=True)
    for run in runs:
        s = run['sel']
        Cx, t = run['C'], run['t'][s]
        for k, (key, name) in enumerate([('M', r'$M$'), ('M_d', r'$M_d$'), ('M_c', r'$M_c$')]):
            axs[0].plot(t, Cx[key][s], color=f'C{k}', lw=run['lw'], label=name if run is runs[0] else None)
        for k, (key, name) in enumerate([('E_rad', r'$E_\mathrm{rad}$'), ('E_th', r'$E_\mathrm{th}$'), ('E_kin', r'$E_\mathrm{kin}$')]):
            axs[1].plot(t, Cx[key][s], color=f'C{k}', lw=run['lw'], label=name if run is runs[0] else None)
        axs[2].plot(t, -Cx['mdot'][s, 0] * SCALE, color='C0', lw=run['lw'], label=r'$\dot M_\mathrm{acc}(r_\mathrm{in}) / \dot M_\mathrm{E}$' if run is runs[0] else None)
        axs[2].plot(t, Cx['L_out'][s] * L0 / L_EDD, color='C1', lw=run['lw'], label=r'$L(r_\mathrm{out}) / L_\mathrm{E}$' if run is runs[0] else None)

    axs[0].set_yscale('log')
    axs[1].set_yscale('log')
    axs[0].set_ylabel(r'$M\ /\ \rho_0 r_g^3$')
    axs[1].set_ylabel(r'$E\ /\ \rho_0 c^2 r_g^3$')
    axs[2].set_ylabel(r'$\lambda_\mathrm{E}$')
    axs[2].set_xlabel(r'$t\ /\ t_g$')
    axs[2].set_xlim(min(run['t_steady'] for run in runs), max(run['t'][-1] for run in runs))

    place_legend(axs[0], ncol=3)
    run_legend(axs[0], runs, ncol=2)
    place_legend(axs[1], ncol=3)
    place_legend(axs[2], ncol=2)

    save(fig, 'window_text')


def spacetime(ax, t, data, norm, cmap, r, tlim=None):
    im = ax.pcolormesh(r, t, data, norm=norm, cmap=cmap, shading='nearest', rasterized=True)
    ax.set_xscale('log')
    ax.set_xlabel(r'$r\ /\ r_g$')
    ax.set_ylabel(r'$t\ /\ t_g$')
    if tlim is not None:
        ax.set_ylim(*tlim)
    plain_log(ax, 'x')
    return im


def fig_spacetime(run):
    C, T, r, s = run['C'], run['t'], run['r'], run['sel']
    TS = T[s]

    mdot = -C['mdot']
    lin = np.nanpercentile(np.abs(mdot[s]) * SCALE, 50)
    vmax = np.nanpercentile(np.abs(mdot) * SCALE, 99.5)
    norm = SymLogNorm(linthresh=max(lin, 11e-2), vmin=-vmax, vmax=vmax)

    fig, axs = plt.subplots(1, 2, figsize=[TEXT, 0.4 * TEXT], sharex=True, sharey=True)
    spacetime(axs[0], T, -C['mdot_d'] * SCALE, norm, 'RdBu_r', r, tlim=(TS[0], TS[-1]))
    im = spacetime(axs[1], T, -C['mdot_c'] * SCALE, norm, 'RdBu_r', r, tlim=(TS[0], TS[-1]))
    axs[1].set_ylabel('')
    axs[0].set_title(r'$\dot M_d$')
    axs[1].set_title(r'$\dot M_c$')
    fig.colorbar(im, ax=axs, label=r'$\dot M\ /\ \dot M_\mathrm{E}$', pad=0.02)
    save(fig, 'spacetime_mdot_text')

    delta_long = C['sigma'][s] / C['sigma'][s].mean(axis=0)[None, :] - 1.0
    amp = np.nanpercentile(np.abs(delta_long), 99)
    dnorm = TwoSlopeNorm(vmin=-amp, vcenter=0.0, vmax=amp)
    tlim = (max(run['S']['t'][0], run['S']['t'][-1] - 3000), run['S']['t'][-1])

    fig, axs = plt.subplots(1, 2, figsize=[TEXT, 0.4 * TEXT], sharex=True)
    spacetime(axs[0], TS, delta_long, dnorm, 'RdBu_r', r)
    im = spacetime(axs[1], run['S']['t'], run['delta'], dnorm, 'RdBu_r', r, tlim=tlim)
    for ax in axs:
        ax.axvline(6.0, color='k', lw=0.5, ls=':')
    fig.colorbar(im, ax=axs, label=r'$\delta = \Sigma / \langle\Sigma\rangle - 1$', pad=0.02)
    save(fig, 'spacetime_delta_text')

    for label, t, data, tl in [('long', TS, delta_long, None), ('dense', run['S']['t'], run['delta'], tlim)]:
        fig, ax = plt.subplots(figsize=[COL, 0.7 * COL])
        im = spacetime(ax, t, data, dnorm, 'RdBu_r', r, tlim=tl)
        ax.axvline(6.0, color='k', lw=0.5, ls=':')
        fig.colorbar(im, ax=ax, label=r'$\delta$', pad=0.02)
        save(fig, f'spacetime_delta_{label}_col')


def fig_profiles(runs):
    run_main = runs[0]
    parts = [('mdot', r'$\dot M$'), ('mdot_d', r'$\dot M_d$'), ('mdot_c', r'$\dot M_c$')]
    flows = [('mdot_in', r'$\dot M_\mathrm{in}$'), ('mdot_out', r'$\dot M_\mathrm{out}$')]

    def mdot_panel(ax):
        for run in runs:
            Cx, s, dt, r = run['C'], run['sel'], run['dt'], run['r']
            for k, (key, lab) in enumerate(parts):
                m, e = mean_se(-Cx[key][s] * SCALE, dt)
                ax.plot(r, m, color=f'C{k}', lw=run['lw'], label=lab if run is run_main else None)
                if run is run_main:
                    ax.fill_between(r, m - e, m + e, color=f'C{k}', alpha=0.3, lw=0)
            if run is run_main:
                for k, (key, lab) in enumerate(flows):
                    ax.plot(r, mean_se(-Cx[key][s] * SCALE, dt)[0], color=f'C{k + 3}', ls='--', lw=0.7, label=lab)
        ax.axhline(0, color='k', lw=0.5)
        ax.set_xscale('log')
        ax.set_yscale('symlog', linthresh=1e-3)
        plain_log(ax, 'x')
        ax.set_xlabel(r'$r\ /\ r_g$')
        ax.set_ylabel(r'$\langle\dot M_\mathrm{acc}\rangle\ /\ \dot M_\mathrm{E}$')
        place_legend(ax, ncol=2, loc='lower left')
        run_legend(ax, runs, ncol=2)

    def sigma_panel(ax):
        for run in runs:
            r = run['r']
            ax.plot(r, run['sig_mean'], color=run['color'])
            ax.fill_between(r, run['sig_mean'] - run['sig_se'], run['sig_mean'] + run['sig_se'], color=run['color'], alpha=0.3, lw=0)
            R_ss = r[r > 6.06]
            if run['mdot_in_mean'] > 0:
                ax.plot(R_ss, run['mdot_in_mean'] / (3.0 * np.pi * nu_kin(R_ss)) * (1.0 - np.sqrt(6.0 / R_ss)), color=run['color'],
                        ls='--', lw=0.7, zorder=-10)
        ax.axvline(6.0, color='gray', lw=0.5, ls=':')
        ax.set_xscale('log')
        ax.set_yscale('log')
        plain_log(ax, 'x')
        ax.set_xlabel(r'$r\ /\ r_g$')
        ax.set_ylabel(r'$\Sigma\ /\ \rho_0 r_g$')
        run_table(ax, [(r'$\;\langle\Sigma\rangle$', dict()), (r'$\;\Sigma_\mathrm{ss}$', dict(ls='--', lw=0.7))],
                  [(run['label'], run['color'], dict(color=run['color'])) for run in runs], loc='lower right', handlelength=2.0)

    def alpha_panel(ax):
        for run in runs:
            r = run['r']
            ratio, rel = effective_alpha(run)
            ax.plot(r, ratio, color=run['color'])
            ax.fill_between(r, ratio * (1 - rel), ratio * (1 + rel), color=run['color'], alpha=0.3, lw=0)
            ax.axvline(run['R_eq_stats'][1], color=run['color'], ls=':', lw=0.6, zorder=-10)
        ax.axhline(1.0, color='k', lw=0.5)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(6, run_main['r'][-1])
        plain_log(ax, 'x')
        plain_log(ax, 'y')
        ax.set_xlabel(r'$R\ /\ r_g$')
        ax.set_ylabel(r'$\alpha_\mathrm{eff} / \alpha$')
        run_legend(ax, runs, color=True, loc='lower right')

    def inflow_panel(ax):
        for run in runs:
            r = run['r']
            ax.plot(r, run['mdot_in_mean'] / (2.0 * np.pi * r * run['sig_mean']), color=run['color'])
        R_ss = run_main['r'][run_main['r'] > 6.5]
        ax.plot(R_ss, 1.5 * nu_kin(R_ss) / R_ss / (1.0 - np.sqrt(6.0 / R_ss)), 'k--', lw=0.7, zorder=-10)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(6, run_main['r'][-1])
        plain_log(ax, 'x')
        ax.set_xlabel(r'$R\ /\ r_g$')
        ax.set_ylabel(r'$\bar v_r\ /\ c$')
        run_legend(ax, runs, color=True, extra=[(r'$\bar v_\mathrm{ss}$', dict(ls='--', lw=0.7))], loc='upper right')

    for name, panels in [('mdot_sigma', (mdot_panel, sigma_panel)), ('alpha_inflow', (alpha_panel, inflow_panel))]:
        fig, axs = plt.subplots(1, 2, figsize=[TEXT, 0.34 * TEXT])
        for ax, panel in zip(axs, panels):
            panel(ax)
        save(fig, f'{name}_text')
        for panel, part in zip(panels, name.split('_')):
            fig, ax = plt.subplots(figsize=[COL, 0.62 * COL])
            panel(ax)
            save(fig, f'{part}_col')


def effective_alpha(run):
    r = run['r']
    w = 1.0 - np.sqrt(6.0 / np.maximum(r, 6.0))
    nu_eff = run['mdot_in_mean'] * w / (3.0 * np.pi * run['sig_mean'])
    a_eff = nu_eff / (2.0 / 3.0 * EPS**2 * np.sqrt(r))
    rel = np.sqrt((run['mdot_in_se'] / run['mdot_in_mean'])**2 + (run['sig_se'] / run['sig_mean'])**2)
    return np.where(r > 6.5, a_eff / ALPHA, np.nan), rel


def welch_len(n, dt):
    nper = max(min(int(2 * n / (WELCH_SEG + 1)), frames(WELCH_TSEG, dt)), 8)
    nper = min(nper, n)
    return nper, (n - nper) // (nper // 2) + 1, max(n // nper, 1)


def welch_kw(n, dt):
    nper = welch_len(n, dt)[0]
    return dict(fs=1.0 / dt, window='hann', nperseg=nper, noverlap=nper // 2, detrend='linear', axis=0)


def peak_frequency(f, P):
    k = int(np.argmax(P[1:])) + 1
    if 1 < k < f.size - 1 and np.all(P[k - 1:k + 2] > 0):
        y = np.log(P[k - 1:k + 2])
        den = y[0] - 2.0 * y[1] + y[2]
        shift = 0.5 * (y[0] - y[2]) / den if den != 0 else 0.0
        return f[k] + np.clip(shift, -0.5, 0.5) * (f[1] - f[0]), f[1] - f[0]
    return f[k], f[1] - f[0]


def contiguous(good, i0):
    out = np.zeros_like(good)
    if not good[i0]:
        return out
    i = i0
    while i < good.size and good[i]:
        out[i] = True
        i += 1
    i = i0 - 1
    while i >= 0 and good[i]:
        out[i] = True
        i -= 1
    return out


def line_fit(x, y, s):
    xm = x - x.mean()
    slope = np.sum(xm * y) / np.sum(xm**2)
    resid = y - y.mean() - slope * xm
    se_res = np.sqrt(np.sum(resid**2) / max(x.size - 2, 1) / np.sum(xm**2))
    se_w = np.sqrt(1.0 / np.sum(xm**2 / np.maximum(s, 1e-6)**2))
    return slope, max(se_res, se_w)


def waves(run):
    r, d, dt = run['r'], run['delta'], run['dt_s']
    w = (r >= R_WAVE[0]) & (r <= R_WAVE[1])
    kw = welch_kw(d.shape[0], dt)
    f, P = signal.welch(d[:, w], **kw)
    glob = np.concatenate([[0.0], (P[1:] / P[1:].sum(axis=0, keepdims=True)).mean(axis=1)])
    run['f_peak'], run['df_peak'] = peak_frequency(f, glob)
    run['f_w'], run['P_w'], run['wave'] = f, P, w

    _, nseg, nind = welch_len(d.shape[0], dt)
    i0 = int(np.argmin(np.abs(r - R_PHASE_REF)))
    fc, P0r = signal.csd(d[:, [i0]], d, **kw)
    _, g2 = signal.coherence(d[:, [i0]], d, **kw)
    _, Prr = signal.welch(d, **kw)
    k = int(np.argmin(np.abs(fc - run['f_peak'])))

    lag = -np.angle(P0r[k])
    ph = np.empty_like(lag)
    ph[i0:] = np.unwrap(lag[i0:])
    ph[:i0 + 1] = np.unwrap(lag[:i0 + 1][::-1])[::-1]
    g2k = np.clip(g2[k], 1e-9, 1.0)
    sig = np.sqrt((1.0 - g2k) / (2.0 * nind * g2k))
    share = Prr[k] / Prr[1:].max(axis=0)
    g2_crit = 1.0 - 0.05**(1.0 / (nind - 1)) if nind > 1 else 1.0
    good = (g2k >= g2_crit) & (share >= P_FUND_MIN) & (r >= R_WAVE[0]) & (r <= R_WAVE[1])
    valid = contiguous(good, i0)

    om = 2.0 * np.pi * fc[k]
    v, ev = np.full(r.size, np.nan), np.full(r.size, np.nan)
    for i in np.nonzero(valid)[0]:
        m = valid & (np.abs(np.log(r / r[i])) <= np.log1p(FIT_FRAC))
        if m.sum() >= 4:
            b, eb = line_fit(r[m], ph[m], sig[m])
            v[i], ev[i] = om / b, abs(om / b) * eb / abs(b)

    run['Ph'] = dict(ph=ph - ph[i0], valid=valid, i0=i0, g2_crit=g2_crit)
    run['v_p'], run['v_p_err'] = v, ev
    run['d_rms'] = np.sqrt(np.mean(d**2, axis=0))


def fig_waves(runs):
    for run in runs:
        waves(run)
    main = runs[0]
    S = main['S']
    rho_m, prs_m = S['rho_mid'].mean(axis=0), S['prs_mid'].mean(axis=0)
    pr_m, vr_m = S['enr_mid'].mean(axis=0) / 3.0, S['vr_mid'].mean(axis=0)
    r = main['r']
    speeds = [(np.sqrt((GAMMA * prs_m + 4.0 / 3.0 * pr_m) / rho_m), r'$c_\mathrm{ad}$'),
              (np.sqrt((prs_m + pr_m) / rho_m), r'$c_\mathrm{iso}$'),
              (np.sqrt(GAMMA * prs_m / rho_m), r'$c_\mathrm{gas}$'),
              (np.abs(vr_m), r'$|\langle v_r\rangle|$'),
              (nu_kin(r) / r, r'$\nu/R$')]

    f_w, P_w, w = main['f_w'], main['P_w'], main['wave']
    fig, ax = plt.subplots(figsize=[COL, 0.62 * COL])
    im = ax.pcolormesh(r[w], f_w[1:], P_w[1:] / P_w[1:].max(axis=0, keepdims=True), cmap='viridis', norm=LogNorm(vmin=1e-3, vmax=1.0),
                       shading='nearest', rasterized=True)
    fig.colorbar(im, ax=ax, label=r'$P / P_\mathrm{max}(r)$', pad=0.02)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(f_w[1], f_w[-1])
    plain_log(ax, 'x')
    ax.set_xlabel(r'$r\ /\ r_g$')
    ax.set_ylabel(r'$\nu_f\ /\ t_g^{-1}$')
    save(fig, 'wave_spectrum_col')

    def phase_panel(ax):
        for run in runs:
            rr, Ph = run['r'], run['Ph']
            v = Ph['valid']
            ax.plot(rr[v], Ph['ph'][v], 'o', ms=2, color=run['color'])
            ax.plot(rr[~v], Ph['ph'][~v], 'o', ms=2, mfc='none', mew=0.5, color=run['color'], alpha=0.4)
        ax.axvline(r[main['Ph']['i0']], color='k', lw=0.5, ls=':')
        ax.set_xscale('log')
        ax.set_xlim(*R_WAVE)
        plain_log(ax, 'x')
        ax.set_xlabel(r'$r\ /\ r_g$')
        ax.set_ylabel(r'$\varphi\ /\ \mathrm{rad}$')
        run_table(ax, [(r'$\gamma^2 \geq \gamma^2_\mathrm{crit}$', dict(marker='o', ls='none', ms=2)),
                       (r'$\gamma^2 < \gamma^2_\mathrm{crit}$', dict(marker='o', ls='none', ms=2, mfc='none', mew=0.5, alpha=0.4))],
                  [(run['label'], run['color'], dict(color=run['color'])) for run in runs], loc='upper left', handlelength=1.0)

    def speed_panel(ax):
        for run in runs:
            ax.errorbar(run['r'], run['v_p'], yerr=run['v_p_err'], fmt='o', ms=2, elinewidth=0.5, color=run['color'], ls='none')
        for y, lab in speeds:
            ax.plot(r, y, lw=0.8, label=lab)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(*R_WAVE)
        plain_log(ax, 'x')
        ax.set_xlabel(r'$r\ /\ r_g$')
        ax.set_ylabel(r'$v\ /\ c$')
        place_legend(ax, ncol=3, loc='lower right')
        run_table(ax, [(r'$v_p$', dict(marker='o', ls='none', ms=2))],
                  [(run['label'], run['color'], dict(color=run['color'])) for run in runs], loc='upper left', handlelength=1.0)

    def rms_panel(ax):
        for run in runs:
            ax.plot(run['r'], run['d_rms'], color=run['color'])
        ax.axvline(6.0, color='gray', lw=0.5, ls=':')
        ax.set_xscale('log')
        ax.set_yscale('log')
        plain_log(ax, 'x')
        plain_log(ax, 'y')
        ax.set_xlabel(r'$r\ /\ r_g$')
        ax.set_ylabel(r'$\delta_\mathrm{rms}$')
        run_legend(ax, runs, color=True, loc='lower right')

    def cross_panel(ax):
        drho = S['rho_mid'] / rho_m[None, :] - 1.0
        for k, (b, lab) in enumerate([(S['vr_mid'] - vr_m[None, :], r'$\varphi_{\delta v_r} - \varphi_{\delta\rho}$'),
                                      (S['enr_mid'] / (3.0 * pr_m[None, :]) - 1.0, r'$\varphi_{\delta E_r} - \varphi_{\delta\rho}$')]):
            kw = welch_kw(drho.shape[0], main['dt_s'])
            fc, P = signal.csd(drho, b, **kw)
            _, g2 = signal.coherence(drho, b, **kw)
            i = int(np.argmin(np.abs(fc - main['f_peak'])))
            m = g2[i] >= main['Ph']['g2_crit']
            ax.plot(r[m], np.angle(P[i])[m], 'o', ms=2, color=f'C{k}', label=lab)
        ax.axhline(0, color='k', lw=0.5)
        for y in (-np.pi / 2, np.pi / 2):
            ax.axhline(y, color='gray', lw=0.5, ls=':')
        ax.set_xscale('log')
        ax.set_xlim(*R_WAVE)
        ax.set_ylim(-np.pi, np.pi)
        plain_log(ax, 'x')
        ax.set_xlabel(r'$r\ /\ r_g$')
        ax.set_ylabel(r'$\Delta\varphi(\nu_\mathrm{peak})\ /\ \mathrm{rad}$')
        place_legend(ax, loc='lower left')

    for name, panels in [('phase_speed', (phase_panel, speed_panel)), ('rms_cross', (rms_panel, cross_panel))]:
        fig, axs = plt.subplots(1, 2, figsize=[TEXT, 0.34 * TEXT])
        for ax, panel in zip(axs, panels):
            panel(ax)
        save(fig, f'{name}_text')
        for panel, part in zip(panels, name.split('_')):
            fig, ax = plt.subplots(figsize=[COL, 0.62 * COL])
            panel(ax)
            save(fig, f'{part}_col')


def compton_opacity(Tg, Tr):
    mc2 = CONST_me * CONST_c**2
    Tr4 = Tr**4
    low = KAPPA_ES * 4.0 * CONST_kB * Tr4 / (mc2 * (Tg + Tr) * (Tg**2 + Tr**2))
    high = KAPPA_ES * 4.0 * CONST_kB * Tr4 * (T_COMPTON_MAX - Tr) / (mc2 * (Tg**4 - Tr4))
    k = np.where(Tg <= T_COMPTON_MAX, low, np.where(Tr < T_COMPTON_MAX, high, 0.0))
    return np.minimum(np.nan_to_num(k), KAPPA_ES)


def heating_cooling(run):
    Cx, s, r = run['C'], run['sel'], run['r']
    rho, p = Cx['rho_mid'][s], Cx['prs_mid'][s]
    E = np.clip(Cx['enr_mid'][s], 0.0, None)
    Tg = K_T * p / rho
    Tr = (E / A_RAD)**0.25

    fd = np.clip((Cx['f_mid'][s] - 0.05) / 0.95, 0.0, None)
    beta = np.sqrt(Cx['vr_mid'][s]**2 + Cx['vphi_mid'][s]**2)
    gv = 1.0 / (1.0 + np.exp(np.minimum((beta - 0.9) / 0.03, 50.0)))

    shear = r * np.gradient(Cx['vphi_mid'][s] / r, r, axis=1)
    u_g = p / (GAMMA - 1.0)
    excess = A_RAD * Tg**4 - E

    q_plus = Cx['nu_mid'][s] * shear**2
    L_C = rho * gv * fd * compton_opacity(Tg, Tr) * KAPPA_CODE * excess
    L_ff = rho * gv * fd * K_FFBF * rho * RHO0 * np.maximum(Tg, 1e4)**-3.5 * KAPPA_CODE * excess

    times = dict(t_heat=u_g / q_plus, t_C=u_g / L_C, t_ff=u_g / L_ff, ratio=q_plus / (L_C + L_ff))
    return {key: np.nanpercentile(np.where(np.isfinite(val) & (val > 0), val, np.nan), [16, 50, 84], axis=0)
            for key, val in times.items()}


def fig_heating(runs):
    times = [('t_heat', r'$t_\mathrm{heat}$'), ('t_C', r'$t_C$'), ('t_ff', r'$t_\mathrm{ff}$')]
    for run in runs:
        run['HC'] = heating_cooling(run)

    def times_panel(ax):
        for run in runs:
            r = run['r']
            m = r >= R_MIN_HC
            for k, (key, _) in enumerate(times):
                ax.plot(r[m], run['HC'][key][1][m], color=f'C{k}', lw=1.2 if run is runs[0] else 0.5)
                if run is runs[0]:
                    ax.fill_between(r[m], run['HC'][key][0][m], run['HC'][key][2][m], color=f'C{k}', alpha=0.2, lw=0)
        R = runs[0]['r'][runs[0]['r'] >= R_MIN_HC]
        ax.plot(R, 1.0 / omega_pw(R), 'k--', lw=0.6)
        ax.set_xscale('log')
        ax.set_yscale('log')
        plain_log(ax, 'x')
        ax.set_xlabel(r'$R\ /\ r_g$')
        ax.set_ylabel(r'$t\ /\ t_g$')
        h = [Line2D([], [], color=f'C{k}') for k in range(3)] + [Line2D([], [], color='k', ls='--', lw=0.6)]
        l = [lab for _, lab in times] + [r'$1/\Omega$']
        place_legend(ax, h, l, ncol=2, loc='upper left')
        run_legend(ax, [dict(label=run['label'], lw=1.2 if run is runs[0] else 0.5) for run in runs], loc='lower right')

    def ratio_panel(ax):
        for run in runs:
            r = run['r']
            m = r >= R_MIN_HC
            ax.plot(r[m], run['HC']['ratio'][1][m], color=run['color'])
            ax.fill_between(r[m], run['HC']['ratio'][0][m], run['HC']['ratio'][2][m], color=run['color'], alpha=0.2, lw=0)
        ax.axhline(1.0, color='k', lw=0.5)
        ax.set_xscale('log')
        ax.set_yscale('log')
        plain_log(ax, 'x')
        plain_log(ax, 'y')
        ax.set_xlabel(r'$R\ /\ r_g$')
        ax.set_ylabel(r'$q^+ / (\Lambda_C + \Lambda_\mathrm{ff})$')
        run_legend(ax, runs, color=True, loc='upper right')

    fig, axs = plt.subplots(1, 2, figsize=[TEXT, 0.34 * TEXT])
    times_panel(axs[0])
    ratio_panel(axs[1])
    save(fig, 'heating_text')
    for panel, part in [(times_panel, 'heating_times'), (ratio_panel, 'heating_ratio')]:
        fig, ax = plt.subplots(figsize=[COL, 0.62 * COL])
        panel(ax)
        save(fig, f'{part}_col')


def shell_mass(Cx, a, b):
    return Cx['m_hi'][:, a] + Cx['m_lo'][:, b] + np.sum(Cx['m_lo'][:, a + 1:b] + Cx['m_hi'][:, a + 1:b], axis=1)


def cumulative_trapezoid(t, y):
    return np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(t))])


def mass_budget(Cx, a, b):
    t = Cx['t']
    net = Cx['mdot'][:, a] - Cx['mdot'][:, b]
    flux_int = cumulative_trapezoid(t, net)

    t_e, t_o = t[0::2], t[1::2]
    I_e = np.interp(t_e[-1], t, flux_int) - np.interp(t_e[0], t, flux_int)
    I_o = np.interp(t_o[-1], t, flux_int) - np.interp(t_o[0], t, flux_int)
    err = max(abs(I_e - cumulative_trapezoid(t_e, net[0::2])[-1]), abs(I_o - cumulative_trapezoid(t_o, net[1::2])[-1]))

    M_ab = shell_mass(Cx, a, b)
    dM = M_ab - M_ab[0]
    return dict(t=t, dM=dM, flux=flux_int, floor=dM - flux_int, err=(t - t[0]) / (t[-1] - t[0]) * err)


def fig_budget(run):
    S, r = run['S'], run['r']
    IA = int(np.argmin(np.abs(r - R_BUDGET[0])))
    IB = int(np.argmin(np.abs(r - R_BUDGET[1])))
    budgets = [('restricted', mass_budget(S, IA, IB)), ('complete', mass_budget(S, 0, r.size - 1))]

    def panel(ax, B, legend):
        ax.plot(B['t'], B['dM'], label=r'$\Delta M_{ab}$')
        ax.plot(B['t'], B['flux'], label=r'$\int [\dot M(r_a) - \dot M(r_b)]\,\mathrm{d}t$')
        ax.plot(B['t'], B['floor'], label=r'$\Delta M_\mathrm{floor}$')
        ax.fill_between(B['t'], B['floor'] - B['err'], B['floor'] + B['err'], alpha=0.3, lw=0, color='firebrick', label=r'$\pm\,\sigma_I$')
        ax.axhline(0, color='k', lw=0.5)
        ax.set_xlabel(r'$t\ /\ t_g$')
        ax.set_ylabel(r'$M\ /\ \rho_0 r_g^3$')
        if legend:
            place_legend(ax, loc='upper left')

    fig, axs = plt.subplots(1, 2, figsize=[TEXT, 0.34 * TEXT])
    for k, (ax, (_, B)) in enumerate(zip(axs, budgets)):
        panel(ax, B, k == 1)
    save(fig, 'budget_text')
    for label, B in budgets:
        fig, ax = plt.subplots(figsize=[COL, 0.62 * COL])
        panel(ax, B, True)
        save(fig, f'budget_{label}_col')


def fig_structure(run):
    path, geo = run['path'], run['geo']
    R_C, TH_C, NR, NTH, DTH, COT = geo['R_C'], geo['TH_C'], geo['NR'], geo['NTH'], geo['DTH'], geo['COT']
    outs = run['C']['outs'][run['sel']][::PROFILE_STRIDE]
    st, ct = np.sin(TH_C)[None, :], np.cos(TH_C)[None, :]
    lr = np.log(R_C)

    def vertical(field, R0):
        r_line = R0 / np.sin(TH_C)
        out = np.full(NTH, np.nan)
        for j in np.nonzero((r_line >= R_C[0]) & (r_line <= R_C[-1]))[0]:
            out[j] = np.interp(np.log(r_line[j]), lr, field[:, j])
        return out

    PROF = {R0: {} for R0 in R_PROFILES}
    QP, QM, NQM = np.zeros(NR), np.zeros(NR), np.zeros(NR)
    for n in outs:
        F = load(n, path, geo)
        Fz = F['fr1'] * ct - F['fr2'] * st
        gz = ct / (R_C[:, None] - 2.0)**2
        kap = F['kappa_abs'] + F['kappa_scat']
        fields = dict(rho=F['rho'], p_gas=F['prs'], p_rad=F['enr'] / 3.0, T_g=gas_temperature(F), T_r=rad_temperature(F),
                      Gamma_z=np.where(gz > 0, kap * Fz / np.maximum(gz, 1e-300), np.nan), f=F['diskfrac'],
                      v_z=F['vx1'] * ct - F['vx2'] * st)
        for R0 in R_PROFILES:
            for key, val in fields.items():
                PROF[R0][key] = PROF[R0].get(key, 0.0) + vertical(val, R0)

        Om = F['vx3'] / (R_C[:, None] * st)
        shear = R_C[:, None] * np.gradient(Om, R_C, axis=0)
        QP += np.sum(F['nu'] * shear**2 * R_C[:, None] * DTH[None, :], axis=1)
        tau = optical_depth_from_pole(F, geo)
        for i in range(NR):
            j = np.argmax(tau[i] >= 1.0)
            if tau[i, j] >= 1.0:
                QM[i] += Fz[i, j]
                NQM[i] += 1

    for R0 in R_PROFILES:
        for key in PROF[R0]:
            PROF[R0][key] = PROF[R0][key] / len(outs)
    QP /= len(outs)
    QM = np.where(NQM > 0, QM / np.maximum(NQM, 1), np.nan)

    fig, axs = plt.subplots(2, 3, figsize=[TEXT, 0.5 * TEXT], sharex=True)
    for R0 in R_PROFILES:
        P = PROF[R0]
        line, = axs[0, 0].semilogy(COT, P['rho'], label=rf'$R = {R0:g}$')
        c = line.get_color()
        axs[0, 1].semilogy(COT, P['p_rad'], color=c)
        axs[0, 1].semilogy(COT, P['p_gas'], color=c, ls='--')
        axs[0, 2].semilogy(COT, P['T_r'], color=c)
        axs[0, 2].semilogy(COT, P['T_g'], color=c, ls='--')
        axs[1, 0].plot(COT, P['Gamma_z'], color=c)
        axs[1, 1].plot(COT, P['v_z'], color=c)
        axs[1, 2].plot(COT, P['f'], color=c)
        surface = COT[np.nanargmin(np.abs(P['f'] - 0.5))]
        for ax in axs.flat:
            ax.axvline(surface, color=c, lw=0.5, ls=':')
    axs[1, 0].set_yscale('symlog', linthresh=1e-1)
    axs[0, 0].set_xlim(0, 0.6)
    for ax, lab in zip(axs.flat, [r'$\langle\rho\rangle\ /\ \rho_0$', r'$\langle p\rangle\ /\ \rho_0 c^2$', r'$\langle T\rangle\ /\ \mathrm{K}$',
                                  r'$\langle\Gamma_z\rangle$', r'$\langle v_z\rangle\ /\ c$', r'$\langle f\rangle$']):
        ax.set_ylabel(lab)
    for ax in axs[1]:
        ax.set_xlabel(r'$z/R$')
    place_legend(axs[0, 0])
    place_legend(axs[0, 1], [Line2D([], [], color='k', label=r'$\langle p_r\rangle$'), Line2D([], [], color='k', ls='--', label=r'$\langle p_g\rangle$')],
                     ncol=2, loc='upper right')
    place_legend(axs[0, 2], [Line2D([], [], color='k', label=r'$\langle T_r\rangle$'), Line2D([], [], color='k', ls='--', label=r'$\langle T_g\rangle$')],
                     ncol=2, loc='upper right')
    save(fig, 'vertical_text')

    t_diff = KAPPA_ES * RHO0 * R_G * run['sig_mean'] / 2.0 * run['hr_mean'] * R_C
    fig, axs = plt.subplots(1, 2, figsize=[TEXT, 0.34 * TEXT])
    axs[0].loglog(R_C, QP, label=r'$Q^+$')
    axs[0].loglog(R_C, np.where(QM > 0, QM, np.nan), label=r'$Q^- = F_z|_{\tau = 1}$')
    axs[0].set_xlabel(r'$r\ /\ r_g$')
    axs[0].set_ylabel(r'$Q\ /\ \rho_0 c^3$')
    place_legend(axs[0], loc='lower left')
    axs[1].loglog(R_C, t_diff, label=r'$t_\mathrm{diff}$')
    axs[1].loglog(R_C, t_thermal(R_C), label=r'$t_\mathrm{th}$')
    axs[1].loglog(R_C, t_viscous(R_C), label=r'$t_\nu$')
    axs[1].axhline(1 / run['f_peak'], color='firebrick', ls='--', lw=0.7, label=r'$1/\nu_\mathrm{peak}$')
    axs[1].set_xlabel(r'$r\ /\ r_g$')
    axs[1].set_ylabel(r'$t\ /\ t_g$')
    place_legend(axs[1], ncol=2, loc='lower right')
    for ax in axs:
        plain_log(ax, 'x')
    save(fig, 'balance_text')


def fig_aspect(runs):
    fig, ax = plt.subplots(figsize=[COL, 0.62 * COL])
    axr = ax.twinx()
    for run in runs:
        ax.plot(run['r'], run['hr_mean'], color=run['color'])
        ax.fill_between(run['r'], run['hr_mean'] - run['hr_se'], run['hr_mean'] + run['hr_se'], color=run['color'], alpha=0.3, lw=0)
        axr.plot(run['r'], run['cells_H'], color=run['color'], lw=0.5, alpha=0.6)
    ax.axhline(EPS, color='k', ls=':', lw=0.6)
    ax.set_xscale('log')
    axr.set_ylim(0, None)
    plain_log(ax, 'x')
    ax.set_xlabel(r'$r\ /\ r_g$')
    ax.set_ylabel(r'$\langle H/R\rangle$')
    axr.set_ylabel(r'$\langle H/R\rangle / \Delta\theta$')
    run_legend(ax, runs, color=True, extra=[(r'$\langle H/R\rangle$', dict()), (r'$\langle H/R\rangle / \Delta\theta$', dict(lw=0.5, alpha=0.6)),
                                            (r'$\epsilon$', dict(ls=':', lw=0.6))], loc='upper right')
    save(fig, 'aspect_col')


class NegativeLog(mpl.ticker.LogFormatterSciNotation):
    def __call__(self, x, pos=None):
        s = super().__call__(x, pos)
        return '$-' + s[1:] if s else s


def angular_profiles(run):
    geo, path = run['geo'], run['path']
    outs = run['C']['outs'][run['sel']][::PROFILE_STRIDE]
    idx = np.array([int(np.argmin(np.abs(geo['R_C'] - r0))) for r0 in R_ANG])
    r = geo['R_C'][idx][:, None]
    sth = np.sin(geo['TH_C'])[None, :]

    rec = {key: [] for key in ['mdot', 'rho', 'vr', 'vphi']}
    for n in outs:
        F = {key: val[idx, :] for key, val in load(n, path, geo, names=('rho', 'vx1', 'vx3', 'diskfrac')).items()}
        rec['mdot'].append(-HEMI * 2.0 * np.pi * r**2 * sth * F['rho'] * F['vx1'])
        rec['rho'].append(F['rho'])
        rec['vr'].append(F['vx1'])
        rec['vphi'].append(F['vx3'])

    out = dict(elev=90.0 - np.degrees(geo['TH_C']), dth=geo['DTH'])
    for key, val in rec.items():
        out[key], out[key + '_se'] = mean_se(np.array(val), run['dt'] * PROFILE_STRIDE)
    total = np.sum(out['mdot'] * out['dth'][None, :], axis=1, keepdims=True)
    out['mdot_dist'] = out['mdot'] / total * np.pi / 180.0
    out['mdot_dist_se'] = out['mdot_se'] / np.abs(total) * np.pi / 180.0
    return out


def fig_angular(runs):
    for run in runs:
        run['ANG'] = angular_profiles(run)
    rcol = plt.cm.viridis(np.linspace(0.0, 0.85, len(R_ANG)))
    lws = (1.2, 0.5)
    panels = [('mdot_dist', r'$\mathrm{d}\dot M_\mathrm{acc}/\mathrm{d}\vartheta\ /\ (\dot M_\mathrm{acc}\,\mathrm{deg}^{-1})$', 1.0),
              ('rho', r'$\langle\rho\rangle\ /\ \rho_0$', 1.0),
              ('vr', r'$\langle v_r\rangle\ /\ c$', -1.0),
              ('vphi', r'$\langle v_\phi\rangle\ /\ c$', None)]

    def panel(ax, key, lab, sign, legend):
        for run, lw in zip(runs, lws):
            A = run['ANG']
            for k in range(len(R_ANG)):
                y, e = A[key][k], A[key + '_se'][k]
                lo, hi = y - e, y + e
                if sign is not None:
                    y, lo, hi = sign * y, np.minimum(sign * lo, sign * hi), np.maximum(sign * lo, sign * hi)
                    y, lo, hi = (np.where(v > 0, v, np.nan) for v in (y, lo, hi))
                    lo = np.where(np.isfinite(lo), lo, y)
                ax.plot(A['elev'], y, color=rcol[k], lw=lw)
                ax.fill_between(A['elev'], lo, hi, color=rcol[k], alpha=0.15, lw=0)
        if sign is not None:
            ax.set_yscale('log')
        if sign == -1.0:
            ax.yaxis.set_major_formatter(NegativeLog())
            ax.yaxis.set_minor_formatter(NegativeLog(labelOnlyBase=False))
        if key == 'mdot_dist':
            top = max(np.nanmax(run['ANG']['mdot_dist']) for run in runs)
            ax.set_ylim(top * 1e-4, top * 2.0)
        ax.set_xscale('log')
        ax.set_xlim(0.3, 90.0)
        plain_log(ax, 'x', mantissas=(1, 2, 5))
        ax.set_ylabel(lab)
        ax.set_xlabel(r'$\vartheta\ /\ \mathrm{deg}$')
        if legend:
            place_legend(ax, [Line2D([], [], color=c, label=rf'$r = {r0:g}$') for r0, c in zip(R_ANG, rcol)], ncol=3,
                         handlelength=1.0, columnspacing=0.8, loc='lower left')
            run_legend(ax, [dict(label=run['label'], lw=lw) for run, lw in zip(runs, lws)], ncol=2)

    fig, axs = plt.subplots(2, 2, figsize=[TEXT, 0.66 * TEXT], sharex=True)
    for k, (ax, (key, lab, sign)) in enumerate(zip(axs.flat, panels)):
        panel(ax, key, lab, sign, k == 0)
    for ax in axs[0]:
        ax.set_xlabel('')
    save(fig, 'angular_text')
    for key, lab, sign in panels:
        fig, ax = plt.subplots(figsize=[COL, 0.62 * COL])
        panel(ax, key, lab, sign, key == 'mdot_dist')
        save(fig, f'angular_{key}_col')


def main():
    hr = make_run(MAIN_PATH, r'$\mathrm{HR}$', 'steelblue', 1.2)
    lr = make_run(REF_PATH, r'$\mathrm{LR}$', 'firebrick', 0.5)
    runs = [hr, lr]

    steps = [('timescales', lambda: fig_timescales()),
             ('valid', lambda: fig_validation()),
             ('domain', lambda: fig_full_domain(hr)),
             ('early', lambda: fig_early(hr)),
             ('relax', lambda: fig_relaxation(runs)),
             ('window', lambda: fig_window(runs)),
             ('spacetime', lambda: fig_spacetime(hr)),
             ('profiles', lambda: fig_profiles(runs)),
             ('waves', lambda: fig_waves(runs)),
             ('heating', lambda: fig_heating(runs)),
             ('budget', lambda: fig_budget(hr)),
             ('aspect', lambda: fig_aspect(runs)),
             ('structure', lambda: fig_structure(hr)),
             ('angular', lambda: fig_angular(runs))]

    for name, step in steps:
        if wanted(name):
            if name == 'structure' and 'f_peak' not in hr:
                waves(hr)
            step()


if __name__ == '__main__':
    main()
