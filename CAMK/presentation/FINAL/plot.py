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
DATA = ROOT.parent.parent / 'PLUTO' / 'Simulation_Main'
OUT = ROOT / 'content'

DPI = int(os.environ.get('PLOT_DPI', '600'))
ONLY = [s for s in os.environ.get('SLIDE_ONLY', '').split(',') if s]

SLIDE_W = 5.80
SLIDE_H = 2.36
PLOT_W = 4.80
WIDE_W = 5.20

preamble = r'''
    \usepackage{amsmath}
    \usepackage{mathtools}

    \usepackage{fontspec}
    \usepackage{unicode-math}
    \setmainfont{Fira Sans Light}[BoldFont=Fira Sans Medium]
    \setsansfont{Fira Sans Light}[BoldFont=Fira Sans Medium]
    \setmathfont{Fira Math}

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
    "font.family": 'sans-serif',
    "font.size": 9,

    "axes.prop_cycle": "cycler('color', ['steelblue', 'olivedrab', 'goldenrod', 'firebrick', 'rebeccapurple'])",
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "axes.linewidth": 0.6,

    "figure.constrained_layout.use": True,
    "figure.constrained_layout.h_pad": 0.6 / 72,
    "figure.constrained_layout.w_pad": 0.6 / 72,

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

from scipy import signal
from scipy.ndimage import uniform_filter1d

import pyPLUTO as pp


MAIN_PATH = DATA / 'BH_VISC_RAD_HD_TRC_COMP_FIX'
REF_PATH = DATA / 'BH_VISC_RAD_HD_TRC_COMP'

TRACER_RUNS = [('HD', DATA / 'BH_VISC_HD', DATA / 'BH_VISC_HD_TRC'),
               ('MHD', DATA / 'BH_VISC_MHD', DATA / 'BH_VISC_MHD_TRC'),
               ('RMHD', DATA / 'BH_VISC_RES_MHD', DATA / 'BH_VISC_RES_MHD_TRC')]
RAD_TRACER_RUN = DATA / 'BH_VISC_RAD_HD_TRC'

T_RAD = [0.0, 250.0, 500.0, 750.0]
T_VALID = [0.0, 250.0, 500.0, 750.0]
T_EARLY = [100, 200, 300, 400, 500, 600]
FULL_AVG_SPAN = 1.0e4
FRAME_PAD = mpl.rcParams['axes.linewidth'] / 72
PAIR_GAP = 0.16

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


def write(fig, name):
    OUT.mkdir(exist_ok=True)
    finalize_legends(fig, name)
    fig.savefig(OUT / f'{name}.pdf', dpi=DPI)
    plt.close(fig)
    print(f'wrote {name}.pdf')


LEGEND_LOCS = ['upper right', 'upper left', 'lower left', 'lower right', 'center right', 'center left',
               'lower center', 'upper center', 'center']


PENDING_LEGENDS = {}


def place_legend(ax, handles=None, labels=None, loc=None, **kw):
    if handles is None:
        handles, labels = ax.get_legend_handles_labels()
    elif labels is None:
        labels = [h.get_label() for h in handles]
    if ax not in PENDING_LEGENDS:
        PENDING_LEGENDS[ax] = []
    PENDING_LEGENDS[ax].append((handles, labels, loc, kw))


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
        pending = PENDING_LEGENDS.pop(ax, [])
        placed = []
        for k, (handles, labels, loc, kw) in enumerate(pending):
            leg, box = best_legend(ax, handles, labels, loc, kw, placed, name)
            placed.append(box)
            if k < len(pending) - 1:
                ax.add_artist(leg)


def plain_label(v, pos):
    return f'${v:g}$'


def plain_log(ax, axis='x', mantissas=(1, 2, 3, 4, 6)):
    a = ax.xaxis if axis == 'x' else ax.yaxis
    lo, hi = sorted(ax.get_xlim() if axis == 'x' else ax.get_ylim())
    ticks = [m * 10.0**k for k in range(int(np.floor(np.log10(lo))) - 1, int(np.ceil(np.log10(hi))) + 1) for m in mantissas]
    a.set_major_locator(FixedLocator([t for t in ticks if lo <= t <= hi]))
    a.set_major_formatter(FuncFormatter(plain_label))
    a.set_minor_formatter(NullFormatter())


def power_label(v, pos):
    if v == 0:
        return '$0$'
    return rf'$10^{{{int(np.floor(np.log10(v)))}}}$'


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


GLOBAL_TEX = {'M_d': r'$M_d$', 'E_rad': r'$E_{rad}$', 'mdot_in': r'$\dot M_{acc}(r_{i{\kern-0.25pt}n})$', 'L_out': r'$L(r_{out})$'}


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


def tighten_map(fig, pad=None, passes=2):
    if pad is None:
        pad = 0.5 * FRAME_PAD + 1.0 / DPI
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
    cbars = list(cbar) if isinstance(cbar, (list, tuple)) else [cbar] * len(groups)
    caxes = [fig.add_axes([0, 0, 1, 1]) if c else None for c in cbars]
    fig._map_spec = dict(groups=groups, ncol=ncol, width=width, aspect=(ylim[1] - ylim[0]) / (xlim[1] - xlim[0]),
                         axes=axes, caxes=caxes)
    place_map(fig, MAP_LEFT, CB_PAD + CB_WIDTH + CB_RIGHT if any(cbars) else 0.1, MAP_BOTTOM, MAP_TOP)

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


TIME_COLOR = {'rho': 'w', 'f': 'w', 'tr': 'w', 'ratio': 'w', 'vr': 'w', 'vphi': 'w', 'Tg': 'k'}


def time_label(ax, t, key=None):
    color = TIME_COLOR.get(key) or contrast_color(ax, 0.55, 1.0, 0.75, 1.0)
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


def build_map_grid(groups, width, xlim=(0, 20), ylim=(0, 10), cbar=True, yaxis=True):
    ncol = len(groups[0][1][0])
    cbars = [cbar and key is not None for key, _ in groups]
    fig, axes, caxes = map_layout([len(rows) for _, rows in groups], ncol, width, xlim, ylim, cbar=cbars)
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
                if key == 'outline':
                    outline_map(axes[i, j], *field, geo, xlim)
                    continue
                im = draw_map(axes[i, j], field, key, geo, xlim, ylim)
                if len(cell) > 5:
                    cell[5](axes[i, j])
                if t is not None:
                    time_label(axes[i, j], t, key)
                if tag is not None:
                    row_tag(axes[i, j], tag)
            i += 1
        if cax is not None:
            map_colorbar(fig, im, cax, cbar_key)
    return fig


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


SNAP_CACHE = {}


def snapshot(path, t_req, names=('rho', 'vx1', 'vx2', 'vx3', 'prs', 'enr', 'diskfrac', 'tr1')):
    n, t = nearest_out(path, t_req)
    key = (str(path), int(n))
    if key not in SNAP_CACHE:
        SNAP_CACHE[key] = (t, make_geometry(path), load(n, path, make_geometry(path), names=names))
    return SNAP_CACHE[key]


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


def vertical_profile(field, R0, geo):
    R_C, TH_C = geo['R_C'], geo['TH_C']
    r_line = R0 / np.sin(TH_C)
    out = np.full(geo['NTH'], np.nan)
    for j in np.nonzero((r_line >= R_C[0]) & (r_line <= R_C[-1]))[0]:
        out[j] = np.interp(np.log(r_line[j]), np.log(R_C), field[:, j])
    return out


def spacetime(ax, t, data, norm, cmap, r, tlim=None):
    im = ax.pcolormesh(r, t, data, norm=norm, cmap=cmap, shading='nearest', rasterized=True)
    ax.set_xscale('log')
    ax.set_xlabel(r'$r\ /\ r_g$')
    ax.set_ylabel(r'$t\ /\ t_g$')
    if tlim is not None:
        ax.set_ylim(*tlim)
    plain_log(ax, 'x')
    return im


def mdot_panel(ax, runs):
    run_main = runs[0]
    parts = [('mdot', r'$\dot M$'), ('mdot_d', r'$\dot M_d$'), ('mdot_c', r'$\dot M_c$')]
    flows = [('mdot_in', r'$\dot M_{i{\kern-0.25pt}n}$'), ('mdot_out', r'$\dot M_{out}$')]
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
    ax.set_ylabel(r'$\langle\dot M_{acc}\rangle\ /\ \dot M_{E}$')
    place_legend(ax, ncol=2, loc='lower left')
    run_legend(ax, runs, ncol=2)


def sigma_panel(ax, runs):
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
    run_table(ax, [(r'$\;\langle\Sigma\rangle$', dict()), (r'$\;\Sigma_{ss}$', dict(ls='--', lw=0.7))],
              [(run['label'], run['color'], dict(color=run['color'])) for run in runs], loc='lower right', handlelength=2.0)


def alpha_panel(ax, runs):
    run_main = runs[0]
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
    ax.set_ylabel(r'$\alpha_{ef{\kern-0.5pt}f} / \alpha$')
    run_legend(ax, runs, color=True, loc='lower right')


def inflow_panel(ax, runs):
    run_main = runs[0]
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
    run_legend(ax, runs, color=True, extra=[(r'$\bar v_{ss}$', dict(ls='--', lw=0.7))], loc='upper right')


def wave_speeds(main):
    S, r = main['S'], main['r']
    rho_m, prs_m = S['rho_mid'].mean(axis=0), S['prs_mid'].mean(axis=0)
    pr_m, vr_m = S['enr_mid'].mean(axis=0) / 3.0, S['vr_mid'].mean(axis=0)
    return [(np.sqrt((GAMMA * prs_m + 4.0 / 3.0 * pr_m) / rho_m), r'$c_{ad}$'),
            (np.sqrt((prs_m + pr_m) / rho_m), r'$c_{i{\kern-0.25pt}so}$'),
            (np.sqrt(GAMMA * prs_m / rho_m), r'$c_{gas}$'),
            (np.abs(vr_m), r'$|\langle v_r\rangle|$'),
            (nu_kin(r) / r, r'$\nu/R$')]


def phase_panel(ax, runs):
    main = runs[0]
    for run in runs:
        rr, Ph = run['r'], run['Ph']
        v = Ph['valid']
        ax.plot(rr[v], Ph['ph'][v], 'o', ms=2, color=run['color'])
        ax.plot(rr[~v], Ph['ph'][~v], 'o', ms=2, mfc='none', mew=0.5, color=run['color'], alpha=0.4)
    ax.axvline(main['r'][main['Ph']['i0']], color='k', lw=0.5, ls=':')
    ax.set_xscale('log')
    ax.set_xlim(*R_WAVE)
    plain_log(ax, 'x')
    ax.set_xlabel(r'$r\ /\ r_g$')
    ax.set_ylabel(r'$\varphi\ /\ \mathrm{rad}$')
    run_table(ax, [(r'$\gamma^2 \geq \gamma^2_{crit}$', dict(marker='o', ls='none', ms=2)),
                   (r'$\gamma^2 < \gamma^2_{crit}$', dict(marker='o', ls='none', ms=2, mfc='none', mew=0.5, alpha=0.4))],
              [(run['label'], run['color'], dict(color=run['color'])) for run in runs], loc='upper left', handlelength=1.0)


def speed_panel(ax, runs):
    r = runs[0]['r']
    for run in runs:
        ax.errorbar(run['r'], run['v_p'], yerr=run['v_p_err'], fmt='o', ms=2, elinewidth=0.5, color=run['color'], ls='none')
    for y, lab in wave_speeds(runs[0]):
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


def rms_panel(ax, runs):
    for run in runs:
        ax.plot(run['r'], run['d_rms'], color=run['color'])
    ax.axvline(6.0, color='gray', lw=0.5, ls=':')
    ax.set_xscale('log')
    ax.set_yscale('log')
    plain_log(ax, 'x')
    plain_log(ax, 'y')
    ax.set_xlabel(r'$r\ /\ r_g$')
    ax.set_ylabel(r'$\delta_{rms}$')
    run_legend(ax, runs, color=True, loc='lower right')


def cross_panel(ax, main):
    S, r = main['S'], main['r']
    rho_m = S['rho_mid'].mean(axis=0)
    pr_m, vr_m = S['enr_mid'].mean(axis=0) / 3.0, S['vr_mid'].mean(axis=0)
    drho = S['rho_mid'] / rho_m[None, :] - 1.0
    pairs = [(S['vr_mid'] - vr_m[None, :], r'$\varphi_{\delta v_r} - \varphi_{\delta\rho}$'),
             (S['enr_mid'] / (3.0 * pr_m[None, :]) - 1.0, r'$\varphi_{\delta E_r} - \varphi_{\delta\rho}$')]
    for k, (b, lab) in enumerate(pairs):
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
    ax.set_ylabel(r'$\Delta\varphi(\nu_{peak})\ /\ \mathrm{rad}$')
    place_legend(ax, loc='lower left')


HEATING_TIMES = [('t_heat', r'$t_{h{\kern-0.25pt}eat}$'), ('t_C', r'$t_C$'), ('t_ff', r'$t_{f{\kern-0.5pt}f}$')]


def heating_times_panel(ax, runs):
    for run in runs:
        r = run['r']
        m = r >= R_MIN_HC
        for k, (key, _) in enumerate(HEATING_TIMES):
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
    handles = [Line2D([], [], color=f'C{k}') for k in range(3)] + [Line2D([], [], color='k', ls='--', lw=0.6)]
    labels = [lab for _, lab in HEATING_TIMES] + [r'$1/\Omega$']
    place_legend(ax, handles, labels, ncol=2, loc='upper left')
    run_legend(ax, [dict(label=run['label'], lw=1.2 if run is runs[0] else 0.5) for run in runs], loc='lower right')


def heating_ratio_panel(ax, runs):
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
    ax.set_ylabel(r'$q^+ / (\Lambda_C + \Lambda_{ff})$')
    run_legend(ax, runs, color=True, loc='upper right')


def budget_panel(ax, B, legend):
    ax.plot(B['t'], B['dM'], label=r'$\Delta M_{ab}$')
    ax.plot(B['t'], B['flux'], label=r'$\int [\dot M(r_a) - \dot M(r_b)]\,dt$')
    ax.plot(B['t'], B['floor'], label=r'$\Delta M_{f{\kern-0.25pt}loor}$')
    ax.fill_between(B['t'], B['floor'] - B['err'], B['floor'] + B['err'], alpha=0.3, lw=0, color='firebrick', label=r'$\pm\,\sigma_I$')
    ax.axhline(0, color='k', lw=0.5)
    ax.set_xlabel(r'$t\ /\ t_g$')
    ax.set_ylabel(r'$M\ /\ \rho_0 r_g^3$')
    if legend:
        place_legend(ax, loc='upper left')


ANGULAR_PANELS = [('mdot_dist', r'$d\dot M/ d\vartheta\ /\ (\dot M\,\mathrm{deg}^{-1})$', 1.0),
                  ('rho', r'$\langle\rho\rangle\ /\ \rho_0$', 1.0),
                  ('vr', r'$\langle v_r\rangle\ /\ c$', -1.0),
                  ('vphi', r'$\langle v_\phi\rangle\ /\ c$', None)]


ANGULAR_WIDTHS = (1.2, 0.5)


def angular_panel(ax, runs, key, lab, sign, legend):
    rcol = plt.cm.viridis(np.linspace(0.0, 0.85, len(R_ANG)))
    for run, lw in zip(runs, ANGULAR_WIDTHS):
        A = run['ANG']
        for k in range(len(R_ANG)):
            y, e = A[key][k], A[key + '_se'][k]
            lo, hi = y - e, y + e
            if sign is not None:
                y, lo, hi = sign * y, np.minimum(sign * lo, sign * hi), np.maximum(sign * lo, sign * hi)
                y = np.where(y > 0, y, np.nan)
                lo = np.where(lo > 0, lo, np.nan)
                hi = np.where(hi > 0, hi, np.nan)
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
        place_legend(ax, [Line2D([], [], color=c, label=rf'$r = {r0:g}$') for r0, c in zip(R_ANG, rcol)], ncol=2,
                     fontsize=6.5, handlelength=1.0, columnspacing=0.6, labelspacing=0.3, handletextpad=0.4, loc='lower left')
        run_legend(ax, [dict(label=run['label'], lw=lw) for run, lw in zip(runs, ANGULAR_WIDTHS)], ncol=2,
                   fontsize=6.5, handlelength=1.2, columnspacing=0.6, labelspacing=0.3, handletextpad=0.4)


def is_colorbar(ax):
    return hasattr(ax, '_colorbar')


def main_axes(fig):
    for ax in fig.axes:
        if not is_colorbar(ax):
            return ax


def align_group(items, size):
    boxes = []
    for fig, name in items:
        fig.set_size_inches(*size)
        fig.draw_without_rendering()
        boxes.append(main_axes(fig).get_position())

    x0 = max(box.x0 for box in boxes)
    y0 = max(box.y0 for box in boxes)
    x1 = min(box.x1 for box in boxes)
    y1 = min(box.y1 for box in boxes)

    for (fig, name), old in zip(items, boxes):
        fig.set_layout_engine('none')
        main_axes(fig).set_position([x0, y0, x1 - x0, y1 - y0])
        for cax in fig.axes:
            if is_colorbar(cax):
                p = cax.get_position()
                cax.set_position([p.x0 + x1 - old.x1, y0, p.width, y1 - y0])
        write(fig, name)


def map_group(items, width=SLIDE_W, height=SLIDE_H):
    for fig, name in items:
        tighten_map(fig)

    margins = np.array([fig._map_spec['margins'] for fig, name in items])
    left, right, bottom, top = margins.max(axis=0)

    spec = items[0][0]._map_spec
    groups = spec['groups']
    ncol = spec['ncol']
    nrow = sum(groups)
    h = (height - bottom - top - (nrow - len(groups)) * MAP_GAP - (len(groups) - 1) * GROUP_GAP) / nrow
    extra = width - left - right - ncol * h / spec['aspect'] - (ncol - 1) * MAP_GAP
    extra = max(extra, 0.0)

    for fig, name in items:
        place_map(fig, left + 0.5 * extra, right + 0.5 * extra, bottom, top)
        write(fig, name)


def cell(path, key, t, label=True, tag=None):
    t_out, geo, F = snapshot(path, t)
    if label:
        return (map_field(F, key), key, geo, t_out, tag)
    return (map_field(F, key), key, geo, None, tag)


def slides_tracer():
    path = RAD_TRACER_RUN
    if not (path / 'dbl.out').exists():
        print('skipping radiation tracer, outputs missing')
        return

    items = []
    for k, t in enumerate(T_RAD):
        row = [cell(path, 'tr', t, tag=r'$f_{tr}$'),
               cell(path, 'f', t, False, r'$f_{d{\kern-0.5pt}f}$')]
        fig = build_map_grid([('f', [row])], SLIDE_W)
        items.append((fig, f'rad_tracer_{k}'))
    map_group(items)


def slides_validation():
    items = []
    for label, tracer_path, df_path in TRACER_RUNS:
        if not (tracer_path / 'dbl.out').exists() or not (df_path / 'dbl.out').exists():
            print(f'skipping validation {label}, outputs missing')
            continue
        for k, t in enumerate(T_VALID):
            density_row = [cell(tracer_path, 'rho', t, tag='Tracer Run'),
                           cell(df_path, 'rho', t, False, 'Classifier Run')]
            fraction_row = [cell(df_path, 'tr', t, False, r'$f_{tr}$'),
                            cell(df_path, 'f', t, False, r'$f_{d{\kern-0.5pt}f}$')]
            fig = build_map_grid([('rho', [density_row]), ('f', [fraction_row])], SLIDE_W)
            items.append((fig, f'valid_{label.lower()}_{k}'))
    if items:
        map_group(items)

    runs = []
    for label, tracer_path, df_path in TRACER_RUNS:
        if (df_path / 'dbl.out').exists():
            runs.append((label, df_path))

    if len(runs) == len(TRACER_RUNS):
        diffs = []
        for k, t in enumerate(T_VALID):
            row = []
            for j, (label, path) in enumerate(runs):
                row.append(cell(path, 'diff', t, j == 0, rf'$\mathrm{{{label}}}$'))
            fig = build_map_grid([('diff', [row])], SLIDE_W)
            diffs.append((fig, f'valid_diff_{k}'))
        map_group(diffs)

    fig, ax = plt.subplots(figsize=[PLOT_W, SLIDE_H])
    shown = False
    for k, (label, path) in enumerate(runs):
        res = agreement_series(path)
        if res is None:
            continue
        t, a = res
        ax.plot(t, 100.0 * a[:, 1], color=f'C{k}', label=rf'$\mathrm{{{label}}}$')
        ax.plot(t, 100.0 * a[:, 0], color=f'C{k}', lw=0.6, ls='--')
        shown = True

    if not shown:
        plt.close(fig)
        return None

    ax.set_xlabel(r'$t\ /\ t_g$')
    ax.set_ylabel(r'$1 - \langle|\Delta f|\rangle\ /\ \%$')
    place_legend(ax, ncol=3)
    return fig


def slides_early(run):
    items = []
    for key in ('rho', 'f', 'Tg', 'ratio', 'vr', 'vphi'):
        for k, t in enumerate(T_EARLY):
            fig = build_map_grid([(key, [[cell(run['path'], key, t)]])], SLIDE_W)
            items.append((fig, f'early_{key}_{k}'))
    map_group(items)


def corner_tag(ax, text):
    ax.text(0.96, 0.93, text, color='k', ha='right', va='top', transform=ax.transAxes, fontsize=8)


def pair_map(cells, geo, lim, tags, right=True):
    fig = plt.figure(figsize=[SLIDE_W, SLIDE_H], layout='none')
    ax0 = fig.add_axes([0, 0, 1, 1])
    ax1 = fig.add_axes([0, 0, 1, 1], sharex=ax0, sharey=ax0)

    caxes = []
    for ax, (field, key), text in zip((ax0, ax1), cells, tags):
        if key == 'outline':
            outline_map(ax, field[0], field[1], geo, lim)
            caxes.append(None)
        else:
            im = draw_map(ax, field, key, geo, lim, lim)
            cax = fig.add_axes([0, 0, 1, 1])
            map_colorbar(fig, im, cax, key)
            caxes.append(cax)
        ax.tick_params(which='both', direction='out', color='k', top=False, right=False, length=2.5)
        ax.tick_params(which='minor', length=1.5)
        ax.set_xlabel(r'$R\ /\ r_g$')
        corner_tag(ax, text)

    ax0.set_ylabel(r'$z\ /\ r_g$')
    ax1.tick_params(labelleft=False)
    ax0.xaxis.set_major_locator(MaxNLocator(5))
    ax0.yaxis.set_major_locator(MaxNLocator(5))

    if not right:
        ax1.set_visible(False)
        if caxes[1] is not None:
            caxes[1].set_visible(False)

    fig._pair = dict(axes=(ax0, ax1), caxes=caxes)
    place_pair(fig, 0.4, 0.4, 1.6, 0.6, 0.0)
    return fig


def place_pair(fig, left, bottom, h, lab, offset):
    ax0, ax1 = fig._pair['axes']
    cax0, cax1 = fig._pair['caxes']
    W, H = fig.get_size_inches()

    x0 = offset + left
    x1 = x0 + h + CB_PAD + CB_WIDTH + lab + PAIR_GAP

    ax0.set_position([x0 / W, bottom / H, h / W, h / H])
    ax1.set_position([x1 / W, bottom / H, h / W, h / H])
    if cax0 is not None:
        cax0.set_position([(x0 + h + CB_PAD) / W, bottom / H, CB_WIDTH / W, h / H])
    if cax1 is not None:
        cax1.set_position([(x1 + h + CB_PAD) / W, bottom / H, CB_WIDTH / W, h / H])


def inches(fig, bbox):
    return bbox.transformed(fig.dpi_scale_trans.inverted())


def tight_inches(fig, ax, renderer):
    return inches(fig, ax.get_tightbbox(renderer))


def box_inches(fig, ax):
    return inches(fig, ax.get_position().transformed(fig.transFigure))


def pair_group(items, width=SLIDE_W, height=SLIDE_H):
    pad = 0.5 * FRAME_PAD + 1.0 / DPI

    need = []
    for fig, name in items:
        fig.draw_without_rendering()
        renderer = fig._get_renderer()
        ax0, ax1 = fig._pair['axes']
        cax0, cax1 = fig._pair['caxes']

        visible = []
        for ax in (ax0, ax1, cax0, cax1):
            if ax is not None and ax.get_visible():
                visible.append(ax)

        labels = []
        for cax in (cax0, cax1):
            if cax is not None and cax.get_visible():
                labels.append(tight_inches(fig, cax, renderer).x1 - box_inches(fig, cax).x1)
            else:
                labels.append(0.0)

        left = box_inches(fig, ax0).x0 - tight_inches(fig, ax0, renderer).x0
        bottom = max(box_inches(fig, ax).y0 - tight_inches(fig, ax, renderer).y0 for ax in visible)
        top = max(tight_inches(fig, ax, renderer).y1 - box_inches(fig, ax).y1 for ax in visible)
        need.append((left, bottom, top, labels[0], labels[1]))

    left, bottom, top, lab0, lab1 = np.max(need, axis=0)
    left = left + pad
    bottom = bottom + pad
    top = max(top, 0.0) + pad

    fixed = left + 2.0 * (CB_PAD + CB_WIDTH) + lab0 + PAIR_GAP + lab1 + pad
    h = min(height - bottom - top, 0.5 * (width - fixed))
    offset = 0.5 * (width - fixed - 2.0 * h)

    for fig, name in items:
        place_pair(fig, left, bottom, h, lab0, offset)
        write(fig, name)


def domain_states(run):
    path, geo = run['path'], run['geo']
    n0 = run['all']['outs'][int(np.argmin(np.abs(run['all']['t'])))]
    F0 = load(n0, path, geo)

    sel = run['t'] >= run['t'][-1] - FULL_AVG_SPAN
    acc = {}
    for n in run['C']['outs'][sel]:
        F = load(n, path, geo)
        values = [('rho', F['rho']), ('Tg', gas_temperature(F)), ('Tr', rad_temperature(F)), ('f', F['diskfrac']),
                  ('vr', F['vx1']), ('vphi', F['vx3']), ('tau', optical_depth_from_pole(F, geo))]
        for key, val in values:
            acc[key] = acc.get(key, 0.0) + val

    mean = {}
    for key, val in acc.items():
        mean[key] = val / sel.sum()
    mean['ratio'] = mean['Tg'] / np.maximum(mean['Tr'], 1.0)

    init = {}
    for key in ('rho', 'f', 'Tg', 'ratio', 'vr', 'vphi'):
        init[key] = map_field(F0, key)
    init['tau'] = optical_depth_from_pole(F0, geo)
    return init, mean


def slides_domain(run):
    geo = run['geo']
    lim = (0, geo['R_E'][-1])
    init, mean = domain_states(run)

    density = [(init['rho'], 'rho'), (mean['rho'], 'rho')]
    pairs = [(pair_map(density, geo, lim, ('Initial', 'Average'), right=False), 'domain_rho_init'),
             (pair_map(density, geo, lim, ('Initial', 'Average')), 'domain_rho_both')]

    for label, tag, state in (('init', 'Initial', init), ('mean', 'Average', mean)):
        disk = [(state['f'], 'f'), ((state['f'], state['tau']), 'outline')]
        thermal = [(state['Tg'], 'Tg'), (state['ratio'], 'ratio')]
        velocity = [(state['vr'], 'vr'), (state['vphi'], 'vphi')]
        pairs.append((pair_map(disk, geo, lim, (tag, tag)), f'domain_disk_{label}'))
        pairs.append((pair_map(thermal, geo, lim, (tag, tag)), f'domain_thermal_{label}'))
        pairs.append((pair_map(velocity, geo, lim, (tag, tag)), f'domain_velocity_{label}'))

    pair_group(pairs)


def plot_relax_mass(runs):
    fig, ax = plt.subplots(figsize=[PLOT_W, SLIDE_H])
    tmax = max(run['t'][-1] for run in runs)

    for run in runs:
        Cx, t = run['C'], run['t']
        for k, (key, name) in enumerate([('M', r'$M$'), ('M_d', r'$M_d$'), ('M_c', r'$M_c$')]):
            ax.plot(t, Cx[key], color=f'C{k}', lw=run['lw'], label=name if run is runs[0] else None)

    ax.set_yscale('log')
    ax.set_xscale('symlog', linthresh=11e3, linscale=1.0)
    ax.set_xlim(0, tmax)
    ax.xaxis.set_major_locator(FixedLocator([0, 1e4]))
    ax.xaxis.set_major_formatter(FuncFormatter(power_label))
    ax.xaxis.set_minor_locator(FixedLocator([k * 1e3 for k in range(1, 51) if k != 10]))
    ax.xaxis.set_minor_formatter(NullFormatter())

    ax.set_xlabel(r'$t\ /\ t_g$')
    ax.set_ylabel(r'$M\ /\ \rho_0 r_g^3$')

    place_legend(ax, ncol=3, loc='upper right')
    run_legend(ax, runs, ncol=2, loc='lower right')
    return fig


def plot_relax_distance(run):
    fig, ax = plt.subplots(figsize=[PLOT_W, SLIDE_H])
    for k, (key, (tc, D)) in enumerate(run['D_curves'].items()):
        ax.plot(tc, D, color=f'C{k}', label=GLOBAL_TEX[key])
    ax.axhline(Z_CRIT, color='k', lw=0.6, ls='--', zorder=-10)
    ax.axvline(run['t_ref'] - 0.5 * W_WINDOW, color='gray', lw=0.6, ls=':', zorder=-10)
    ax.set_yscale('log')
    ax.set_xlabel(r'$t_c\ /\ t_g$')
    ax.set_ylabel(r'$D(t_c)$')
    place_legend(ax, ncol=2, loc='upper right')
    return fig


def plot_timescales():
    R = np.geomspace(2.5, 30, 300)
    fig, ax = plt.subplots(figsize=[PLOT_W, SLIDE_H])

    ax.axvline(6.0, color='gray', lw=0.6, ls=':')
    ax.plot(R, t_orbit(R), label=r'$t_{orb}$')
    ax.plot(R, t_thermal(R), label=r'$t_{th}$')
    ax.plot(R, t_sound(R), label=r'$t_{c_s}$')
    ax.plot(R, t_viscous(R), label=r'$t_\nu$')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(2.5, 30)
    plain_log(ax, 'x')

    ax.set_xlabel(r'$R\ /\ r_g$')
    ax.set_ylabel(r'$t\ /\ t_g$')
    place_legend(ax, ncol=2)
    return fig


def plot_panel(panel, *args):
    fig, ax = plt.subplots(figsize=[PLOT_W, SLIDE_H])
    panel(ax, *args)
    return fig


def plot_spectrum(main):
    r = main['r']
    f_w, P_w, w = main['f_w'], main['P_w'], main['wave']
    fig, ax = plt.subplots(figsize=[WIDE_W, SLIDE_H])
    im = ax.pcolormesh(r[w], f_w[1:], P_w[1:] / P_w[1:].max(axis=0, keepdims=True), cmap='viridis', norm=LogNorm(vmin=1e-3, vmax=1.0),
                       shading='nearest', rasterized=True)
    fig.colorbar(im, ax=ax, label=r'$P / P_{max}(r)$', pad=0.02)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(f_w[1], f_w[-1])
    plain_log(ax, 'x')
    ax.set_xlabel(r'$r\ /\ r_g$')
    ax.set_ylabel(r'$\nu_f\ /\ t_g^{-1}$')
    return fig


def plot_rms_cross(runs):
    fig, axs = plt.subplots(1, 2, figsize=[SLIDE_W, SLIDE_H])
    rms_panel(axs[0], runs)
    cross_panel(axs[1], runs[0])
    return fig


def plot_spacetime_delta(run):
    C, r, s = run['C'], run['r'], run['sel']
    delta_long = C['sigma'][s] / C['sigma'][s].mean(axis=0)[None, :] - 1.0
    amp = np.nanpercentile(np.abs(delta_long), 99)
    dnorm = TwoSlopeNorm(vmin=-amp, vcenter=0.0, vmax=amp)
    tlim = (max(run['S']['t'][0], run['S']['t'][-1] - 3000), run['S']['t'][-1])

    fig, ax = plt.subplots(figsize=[WIDE_W, SLIDE_H])
    im = spacetime(ax, run['S']['t'], run['delta'], dnorm, 'RdBu_r', r, tlim=tlim)
    ax.axvline(6.0, color='k', lw=0.5, ls=':')
    fig.colorbar(im, ax=ax, label=r'$\delta$', pad=0.02)
    return fig


def plot_spacetime_mdot(run):
    C, T, r, s = run['C'], run['t'], run['r'], run['sel']
    TS = T[s]

    mdot = -C['mdot']
    lin = np.nanpercentile(np.abs(mdot[s]) * SCALE, 50)
    vmax = np.nanpercentile(np.abs(mdot) * SCALE, 99.5)
    norm = SymLogNorm(linthresh=max(lin, 11e-2), vmin=-vmax, vmax=vmax)

    fig, axs = plt.subplots(1, 2, figsize=[SLIDE_W, SLIDE_H], sharex=True, sharey=True)
    spacetime(axs[0], T, -C['mdot_d'] * SCALE, norm, 'RdBu_r', r, tlim=(TS[0], TS[-1]))
    im = spacetime(axs[1], T, -C['mdot_c'] * SCALE, norm, 'RdBu_r', r, tlim=(TS[0], TS[-1]))
    axs[1].set_ylabel('')
    axs[0].set_title(r'$\dot M_d$')
    axs[1].set_title(r'$\dot M_c$')
    fig.colorbar(im, ax=axs, label=r'$\dot M\ /\ \dot M_E$', pad=0.02)
    return fig


def plot_budget(run):
    S, r = run['S'], run['r']
    IA = int(np.argmin(np.abs(r - R_BUDGET[0])))
    IB = int(np.argmin(np.abs(r - R_BUDGET[1])))

    fig, axs = plt.subplots(1, 2, figsize=[SLIDE_W, SLIDE_H])
    budget_panel(axs[0], mass_budget(S, IA, IB), False)
    budget_panel(axs[1], mass_budget(S, 0, r.size - 1), True)
    return fig


def plot_aspect(runs):
    fig, ax = plt.subplots(figsize=[PLOT_W, SLIDE_H])
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
    extra = [(r'$\langle H/R\rangle$', dict()), (r'$\langle H/R\rangle / \Delta\theta$', dict(lw=0.5, alpha=0.6)),
             (r'$\epsilon$', dict(ls=':', lw=0.6))]
    run_legend(ax, runs, color=True, extra=extra, loc='upper right')
    return fig


def structure_profiles(run):
    path, geo = run['path'], run['geo']
    R_C, TH_C, NR, DTH = geo['R_C'], geo['TH_C'], geo['NR'], geo['DTH']
    outs = run['C']['outs'][run['sel']][::PROFILE_STRIDE]
    st, ct = np.sin(TH_C)[None, :], np.cos(TH_C)[None, :]

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
                PROF[R0][key] = PROF[R0].get(key, 0.0) + vertical_profile(val, R0, geo)

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
    return PROF, QP, QM


def plot_vertical(run, PROF):
    COT = run['geo']['COT']
    fig, axs = plt.subplots(2, 3, figsize=[SLIDE_W, SLIDE_H], sharex=True)
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

    axs[0, 0].set_ylabel(r'$\langle\rho\rangle\ /\ \rho_0$')
    axs[0, 1].set_ylabel(r'$\langle p\rangle\ /\ \rho_0 c^2$')
    axs[0, 2].set_ylabel(r'$\langle T\rangle\ /\ \mathrm{K}$')
    axs[1, 0].set_ylabel(r'$\langle\Gamma_z\rangle$')
    axs[1, 1].set_ylabel(r'$\langle v_z\rangle\ /\ c$')
    axs[1, 2].set_ylabel(r'$\langle f\rangle$')
    for ax in axs[1]:
        ax.set_xlabel(r'$z/R$')

    small = dict(fontsize=6.5, handlelength=1.2, columnspacing=0.6, labelspacing=0.3, handletextpad=0.4)
    place_legend(axs[0, 0], **small)
    place_legend(axs[0, 1], [Line2D([], [], color='k', label=r'$\langle p_r\rangle$'), Line2D([], [], color='k', ls='--', label=r'$\langle p_g\rangle$')],
                 ncol=2, loc='upper right', **small)
    place_legend(axs[0, 2], [Line2D([], [], color='k', label=r'$\langle T_r\rangle$'), Line2D([], [], color='k', ls='--', label=r'$\langle T_g\rangle$')],
                 ncol=2, loc='upper right', **small)
    return fig


def plot_balance(run, QP, QM):
    R_C = run['geo']['R_C']
    t_diff = KAPPA_ES * RHO0 * R_G * run['sig_mean'] / 2.0 * run['hr_mean'] * R_C

    fig, axs = plt.subplots(1, 2, figsize=[SLIDE_W, SLIDE_H])
    axs[0].loglog(R_C, QP, label=r'$Q^+$')
    axs[0].loglog(R_C, np.where(QM > 0, QM, np.nan), label=r'$Q^- = F_z|_{\tau = 1}$')
    axs[0].set_xlabel(r'$r\ /\ r_g$')
    axs[0].set_ylabel(r'$Q\ /\ \rho_0 c^3$')
    place_legend(axs[0], loc='lower left')

    axs[1].loglog(R_C, t_diff, label=r'$t_{d{\kern-0.25pt}i{\kern-0.5pt}f{\kern-0.5pt}f}$')
    axs[1].loglog(R_C, t_thermal(R_C), label=r'$t_{th}$')
    axs[1].loglog(R_C, t_viscous(R_C), label=r'$t_\nu$')
    axs[1].axhline(1 / run['f_peak'], color='firebrick', ls='--', lw=0.7, label=r'$1/\nu_{peak}$')
    axs[1].set_xlabel(r'$r\ /\ r_g$')
    axs[1].set_ylabel(r'$t\ /\ t_g$')
    place_legend(axs[1], ncol=2, loc='lower right')

    plain_log(axs[0], 'x')
    plain_log(axs[1], 'x')
    return fig


def plot_angular(runs):
    fig, axs = plt.subplots(2, 2, figsize=[SLIDE_W, SLIDE_H], sharex=True)
    for k in range(4):
        key, lab, sign = ANGULAR_PANELS[k]
        angular_panel(axs.flat[k], runs, key, lab, sign, k == 0)
    for ax in axs[0]:
        ax.set_xlabel('')
    return fig


def slides_curves(runs, agreement=None):
    hr, lr = runs

    for run in runs:
        run['HC'] = heating_cooling(run)
    for run in runs:
        waves(run)

    group = [(plot_relax_mass(runs), 'relax_mass'),
             (plot_panel(mdot_panel, runs), 'mdot'),
             (plot_panel(sigma_panel, runs), 'sigma'),
             (plot_panel(alpha_panel, runs), 'alpha'),
             (plot_panel(inflow_panel, runs), 'inflow'),
             (plot_panel(heating_times_panel, runs), 'heating_times'),
             (plot_panel(heating_ratio_panel, runs), 'heating_ratio'),
             (plot_panel(phase_panel, runs), 'phase'),
             (plot_panel(speed_panel, runs), 'speed'),
             (plot_relax_distance(hr), 'relax_distance'),
             (plot_timescales(), 'timescales')]
    if agreement is not None:
        group.append((agreement, 'valid_agreement'))
    align_group(group, (PLOT_W, SLIDE_H))

    align_group([(plot_spacetime_delta(hr), 'spacetime_delta'), (plot_spectrum(hr), 'wave_spectrum')], (WIDE_W, SLIDE_H))

    write(plot_spacetime_mdot(hr), 'spacetime_mdot')
    write(plot_rms_cross(runs), 'rms_cross')
    write(plot_budget(hr), 'budget')
    write(plot_aspect(runs), 'aspect')

    PROF, QP, QM = structure_profiles(hr)
    write(plot_vertical(hr, PROF), 'vertical')
    write(plot_balance(hr, QP, QM), 'balance')

    for run in runs:
        run['ANG'] = angular_profiles(run)
    write(plot_angular(runs), 'angular')


def main():
    agreement = None
    if wanted('tracer'):
        slides_tracer()
    if wanted('valid'):
        agreement = slides_validation()

    if wanted('early') or wanted('domain') or wanted('curves'):
        hr = make_run(MAIN_PATH, r'$\mathrm{HR}$', 'steelblue', 1.2)
        lr = make_run(REF_PATH, r'$\mathrm{LR}$', 'firebrick', 0.5)
        runs = [hr, lr]

        if wanted('early'):
            slides_early(hr)
        if wanted('domain'):
            slides_domain(hr)
        if wanted('curves'):
            slides_curves(runs, agreement)
            agreement = None

    if agreement is not None:
        align_group([(agreement, 'valid_agreement')], (PLOT_W, SLIDE_H))


if __name__ == '__main__':
    main()
