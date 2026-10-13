#!/usr/bin/env python

import os
import importlib.util
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
REPORT = ROOT.parent.parent / 'report'
OUT = ROOT / 'content'

ONLY = [s for s in os.environ.get('SLIDE_ONLY', '').split(',') if s]

spec = importlib.util.spec_from_file_location('report_plot', REPORT / 'plot.py')
R = importlib.util.module_from_spec(spec)
spec.loader.exec_module(R)

import matplotlib as mpl
import matplotlib.pyplot as plt


# --- Slide Geometry & Fonts --- #

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
    "pgf.preamble": preamble,
    "font.family": 'sans-serif',
})

T_RAD = [0.0, 250.0, 500.0, 750.0]


def wanted(name):
    return not ONLY or any(name.startswith(key) for key in ONLY)


# --- Saving --- #

CAPTURED = {}


def capture(fig, name):
    CAPTURED[name] = fig


R.save = capture


def take(*names):
    figs = [CAPTURED.pop(name) for name in names]
    for fig in CAPTURED.values():
        plt.close(fig)
    CAPTURED.clear()
    return figs


def write(fig, name):
    OUT.mkdir(exist_ok=True)
    R.finalize_legends(fig, name)
    fig.savefig(OUT / f'{name}.pdf', dpi=R.DPI)
    plt.close(fig)
    print(f'wrote {name}.pdf')


# --- Alignment --- #

def main_axes(fig):
    return [ax for ax in fig.axes if not hasattr(ax, '_colorbar')]


def align_group(items, size):
    """Draw every figure once with constrained layout at the same size, then fix the common inner box
    of the main axes on all of them, so that the axes stay exactly in place between overlays."""
    boxes = []
    for fig, _ in items:
        fig.set_size_inches(*size)
        fig.draw_without_rendering()
        boxes.append(main_axes(fig)[0].get_position())
    x0 = max(b.x0 for b in boxes)
    y0 = max(b.y0 for b in boxes)
    x1 = min(b.x1 for b in boxes)
    y1 = min(b.y1 for b in boxes)
    for (fig, name), old in zip(items, boxes):
        fig.set_layout_engine('none')
        main_axes(fig)[0].set_position([x0, y0, x1 - x0, y1 - y0])
        for cax in fig.axes:
            if hasattr(cax, '_colorbar'):
                p = cax.get_position()
                cax.set_position([p.x0 + x1 - old.x1, y0, p.width, y1 - y0])
        write(fig, name)


def single(fig, name, size):
    fig.set_size_inches(*size)
    write(fig, name)


def map_group(items, width=SLIDE_W, height=SLIDE_H):
    """Tighten all map figures of one layout, take the largest margin on each side and place every figure
    with these margins, so that all panels sit at the same position on canvases of the same size."""
    figs = [fig for fig, _ in items]
    for fig in figs:
        R.tighten_map(fig)
    left, right, bottom, top = np.max([fig._map_spec['margins'] for fig in figs], axis=0)
    spec = figs[0]._map_spec
    groups, ncol = spec['groups'], spec['ncol']
    nrow = sum(groups)
    h = (height - bottom - top - (nrow - len(groups)) * R.MAP_GAP - (len(groups) - 1) * R.GROUP_GAP) / nrow
    extra = max(width - left - right - ncol * h / spec['aspect'] - (ncol - 1) * R.MAP_GAP, 0.0)
    for fig, name in items:
        R.place_map(fig, left + 0.5 * extra, right + 0.5 * extra, bottom, top)
        write(fig, name)


def cell(path, key, t, label=True, tag=None):
    t_out, geo, F = R.snapshot(path, t)
    return (R.map_field(F, key), key, geo, t_out if label else None, tag)


# --- Classifier --- #

def slides_tracer():
    path = R.RAD_TRACER_RUN
    if not (path / 'dbl.out').exists():
        print('skipping radiation tracer, outputs missing')
        return
    items = []
    for k, t in enumerate(T_RAD):
        groups = [('f', [[cell(path, 'tr', t, tag=r'$f_\mathrm{tr}$'), cell(path, 'f', t, False, r'$f_\mathrm{df}$')]])]
        items.append((R.build_map_grid(groups, SLIDE_W), f'rad_tracer_{k}'))
    map_group(items)


def slides_validation():
    items, diffs = [], []
    for label, tracer_path, df_path in R.TRACER_RUNS:
        if not (tracer_path / 'dbl.out').exists() or not (df_path / 'dbl.out').exists():
            print(f'skipping validation {label}, outputs missing')
            continue
        for k, t in enumerate(R.T_VALID):
            groups = [('rho', [[cell(tracer_path, 'rho', t, tag='Tracer Run'), cell(df_path, 'rho', t, False, 'Classifier Run')]]),
                      ('f', [[cell(df_path, 'tr', t, False, r'$f_\mathrm{tr}$'), cell(df_path, 'f', t, False, r'$f_\mathrm{df}$')]])]
            items.append((R.build_map_grid(groups, SLIDE_W), f'valid_{label.lower()}_{k}'))
    if items:
        map_group(items)

    runs = [(label, df_path) for label, _, df_path in R.TRACER_RUNS if (df_path / 'dbl.out').exists()]
    if len(runs) == len(R.TRACER_RUNS):
        for k, t in enumerate(R.T_VALID):
            row = [cell(path, 'diff', t, j == 0, rf'$\mathrm{{{label}}}$') for j, (label, path) in enumerate(runs)]
            diffs.append((R.build_map_grid([('diff', [row])], SLIDE_W), f'valid_diff_{k}'))
        map_group(diffs)

    fig, ax = plt.subplots()
    shown = False
    for k, (label, path) in enumerate(runs):
        res = R.agreement_series(path)
        if res is None:
            continue
        t, a = res
        ax.plot(t, 100.0 * a[:, 1], color=f'C{k}', label=rf'$\mathrm{{{label}}}$')
        ax.plot(t, 100.0 * a[:, 0], color=f'C{k}', lw=0.6, ls='--')
        shown = True
    if shown:
        ax.set_xlabel(r'$t\ /\ t_g$')
        ax.set_ylabel(r'$1 - \langle|\Delta f|\rangle\ /\ \%$')
        R.place_legend(ax, ncol=3)
        return fig
    plt.close(fig)
    return None


# --- Global Maps --- #

def slides_early(run):
    items = []
    for key in ('rho', 'f', 'Tg', 'ratio', 'vr', 'vphi'):
        for k, t in enumerate(R.T_EARLY):
            items.append((R.build_map_grid([(key, [[cell(run['path'], key, t)]])], SLIDE_W), f'early_{key}_{k}'))
    map_group(items)


def slides_domain(run):
    path, geo = run['path'], run['geo']
    lim = (0, geo['R_E'][-1])
    n0 = run['all']['outs'][int(np.argmin(np.abs(run['all']['t'])))]
    F0 = R.load(n0, path, geo)

    sel = run['t'] >= run['t'][-1] - R.FULL_AVG_SPAN
    acc = {}
    for n in run['C']['outs'][sel]:
        F = R.load(n, path, geo)
        for key, val in [('rho', F['rho']), ('Tg', R.gas_temperature(F)), ('Tr', R.rad_temperature(F)), ('f', F['diskfrac']),
                         ('vr', F['vx1']), ('vphi', F['vx3']), ('tau', R.optical_depth_from_pole(F, geo))]:
            acc[key] = acc.get(key, 0.0) + val
    mean = {key: val / sel.sum() for key, val in acc.items()}
    mean['ratio'] = mean['Tg'] / np.maximum(mean['Tr'], 1.0)
    init = {key: R.map_field(F0, key) for key in ('rho', 'f', 'Tg', 'ratio', 'vr', 'vphi')}
    init['tau'] = R.optical_depth_from_pole(F0, geo)

    items = []
    for label, tag, state in (('init', 'Initial', init), ('mean', 'Average', mean)):
        for key in ('rho', 'f', 'Tg', 'ratio', 'vr', 'vphi'):
            fig = R.build_map_grid([(key, [[(state[key], key, geo, None, tag)]])], SLIDE_W, lim, lim)
            items.append((fig, f'domain_{key}_{label}'))
        fig = R.build_map_grid([(None, [[((state['f'], state['tau']), 'outline', geo, None, None)]])], SLIDE_W, lim, lim)
        fig.axes[0].text(0.04, 0.93, tag, color='k', ha='left', va='top', transform=fig.axes[0].transAxes, fontsize=8)
        items.append((fig, f'domain_outline_{label}'))
    map_group(items)


# --- Profiles & Statistics --- #

def slides_curves(runs, agreement=None):
    hr = runs[0]
    R.fig_relaxation(runs)
    R.fig_profiles(runs)
    R.fig_heating(runs)
    R.fig_waves(runs)
    R.fig_timescales()
    curves = take('relax_mass_col', 'mdot_col', 'sigma_col', 'alpha_col', 'inflow_col', 'heating_times_col', 'heating_ratio_col',
                  'phase_col', 'speed_col', 'relax_distance_col', 'timescales_col', 'wave_spectrum_col', 'rms_cross_text')
    *lines, spectrum, rms_cross = curves
    names = ['relax_mass', 'mdot', 'sigma', 'alpha', 'inflow', 'heating_times', 'heating_ratio', 'phase', 'speed',
             'relax_distance', 'timescales']

    group = list(zip(lines, names))
    if agreement is not None:
        group.append((agreement, 'valid_agreement'))
    align_group(group, (PLOT_W, SLIDE_H))

    R.fig_spacetime(hr)
    delta, mdot = take('spacetime_delta_dense_col', 'spacetime_mdot_text')
    align_group([(delta, 'spacetime_delta'), (spectrum, 'wave_spectrum')], (WIDE_W, SLIDE_H))
    single(mdot, 'spacetime_mdot', (SLIDE_W, SLIDE_H))
    single(rms_cross, 'rms_cross', (SLIDE_W, SLIDE_H))

    R.fig_budget(hr)
    single(*take('budget_text'), 'budget', (SLIDE_W, SLIDE_H))
    R.fig_aspect(runs)
    single(*take('aspect_col'), 'aspect', (PLOT_W, SLIDE_H))
    R.fig_structure(hr)
    vertical, balance = take('vertical_text', 'balance_text')
    single(vertical, 'vertical', (SLIDE_W, SLIDE_H))
    single(balance, 'balance', (SLIDE_W, SLIDE_H))
    R.fig_angular(runs)
    single(*take('angular_text'), 'angular', (SLIDE_W, SLIDE_H))


def main():
    agreement = None
    if wanted('tracer'):
        slides_tracer()
    if wanted('valid'):
        agreement = slides_validation()

    if not any(wanted(key) for key in ('early', 'domain', 'curves')):
        if agreement is not None:
            align_group([(agreement, 'valid_agreement')], (PLOT_W, SLIDE_H))
        return
    hr = R.make_run(R.MAIN_PATH, r'$\mathrm{HR}$', 'steelblue', 1.2)
    lr = R.make_run(R.REF_PATH, r'$\mathrm{LR}$', 'firebrick', 0.5)
    runs = [hr, lr]

    if wanted('early'):
        slides_early(hr)
    if wanted('domain'):
        slides_domain(hr)
    if wanted('curves'):
        slides_curves(runs, agreement)
    elif agreement is not None:
        align_group([(agreement, 'valid_agreement')], (PLOT_W, SLIDE_H))


if __name__ == '__main__':
    main()
