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

from matplotlib.ticker import MaxNLocator


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

    rho = [(init['rho'], 'rho'), (mean['rho'], 'rho')]
    pairs = [(pair_map(rho, geo, lim, ('Initial', 'Average'), right=False), 'domain_rho_init'),
             (pair_map(rho, geo, lim, ('Initial', 'Average')), 'domain_rho_both')]
    for label, tag, state in (('init', 'Initial', init), ('mean', 'Average', mean)):
        for name, cells in [('disk', [(state['f'], 'f'), ((state['f'], state['tau']), 'outline')]),
                            ('thermal', [(state['Tg'], 'Tg'), (state['ratio'], 'ratio')]),
                            ('velocity', [(state['vr'], 'vr'), (state['vphi'], 'vphi')])]:
            pairs.append((pair_map(cells, geo, lim, tag), f'domain_{name}_{label}'))
    pair_group(pairs)


def corner_tag(ax, text):
    """The top right corner of the quarter disk lies outside the domain, on the light background."""
    ax.text(0.96, 0.93, text, color='k', ha='right', va='top', transform=ax.transAxes, fontsize=8)


PAIR_GAP = 0.16


def pair_map(cells, geo, lim, tag, right=True):
    """Two square quarter disk panels side by side, each with its own colorbar on its right, where an outline
    panel keeps an empty colorbar slot so that all pairs share the same layout. The tag is one text for both
    panels or one per panel, and without the right panel only the left one is shown at its usual place."""
    tags = (tag, tag) if isinstance(tag, str) else tag
    fig = plt.figure(figsize=[SLIDE_W, SLIDE_H], layout='none')
    ax0 = fig.add_axes([0, 0, 1, 1])
    ax1 = fig.add_axes([0, 0, 1, 1], sharex=ax0, sharey=ax0)
    caxes = []
    for ax, (field, key), text in zip((ax0, ax1), cells, tags):
        if key == 'outline':
            R.outline_map(ax, *field, geo, lim)
            caxes.append(None)
        else:
            im = R.draw_map(ax, field, key, geo, lim, lim)
            cax = fig.add_axes([0, 0, 1, 1])
            R.map_colorbar(fig, im, cax, key)
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
        for a in (ax1, caxes[1]):
            if a is not None:
                a.set_visible(False)
    fig._pair = dict(axes=(ax0, ax1), caxes=caxes)
    place_pair(fig, 0.4, 0.4, 1.6, 0.6, 0.0)
    return fig


def place_pair(fig, left, bottom, h, lab, offset):
    (ax0, ax1), caxes = fig._pair['axes'], fig._pair['caxes']
    W, H = fig.get_size_inches()
    x0 = offset + left
    x1 = x0 + h + R.CB_PAD + R.CB_WIDTH + lab + PAIR_GAP
    for ax, cax, x in ((ax0, caxes[0], x0), (ax1, caxes[1], x1)):
        ax.set_position([x / W, bottom / H, h / W, h / H])
        if cax is not None:
            cax.set_position([(x + h + R.CB_PAD) / W, bottom / H, R.CB_WIDTH / W, h / H])


def pair_group(items, width=SLIDE_W, height=SLIDE_H):
    """Measure label extents of all pair figures, take the largest of each and place every figure with the
    same panel size and positions, centred on canvases of the same size."""
    pad = 0.5 * R.FRAME_PAD + 1.0 / R.DPI
    need = []
    for fig, _ in items:
        fig.draw_without_rendering()
        r = fig._get_renderer()
        (ax0, ax1), caxes = fig._pair['axes'], fig._pair['caxes']
        inch = lambda a: a.get_tightbbox(r).transformed(fig.dpi_scale_trans.inverted())
        box = lambda a: a.get_position().transformed(fig.transFigure).transformed(fig.dpi_scale_trans.inverted())
        parts = [a for a in (ax0, ax1, *caxes) if a is not None and a.get_visible()]
        labs = [inch(c).x1 - box(c).x1 if c is not None and c.get_visible() else 0.0 for c in caxes]
        need.append((box(ax0).x0 - inch(ax0).x0,
                     max(box(a).y0 - inch(a).y0 for a in parts),
                     max(inch(a).y1 - box(a).y1 for a in parts),
                     labs[0], labs[1]))
    left, bottom, top, lab0, lab1 = np.max(need, axis=0)
    left, bottom, top = left + pad, bottom + pad, max(top, 0.0) + pad
    fixed = left + 2.0 * (R.CB_PAD + R.CB_WIDTH) + lab0 + PAIR_GAP + lab1 + pad
    h = min(height - bottom - top, 0.5 * (width - fixed))
    offset = 0.5 * (width - fixed - 2.0 * h)
    for fig, name in items:
        place_pair(fig, left, bottom, h, lab0, offset)
        write(fig, name)


def shrink_legends(fig, size=6.5, ncol=None):
    for ax in fig.axes:
        for handles, _, _, kw in ax.__dict__.get('_pending_legends', []):
            kw.update(fontsize=size, handlelength=min(kw.get('handlelength', 1.5), 1.2), columnspacing=0.6,
                      labelspacing=0.3, handletextpad=0.4)
            if ncol is not None and len(handles) > 4:
                kw['ncol'] = ncol


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
    shrink_legends(vertical)
    single(vertical, 'vertical', (SLIDE_W, SLIDE_H))
    single(balance, 'balance', (SLIDE_W, SLIDE_H))
    R.fig_angular(runs)
    angular, = take('angular_text')
    shrink_legends(angular, ncol=2)
    single(angular, 'angular', (SLIDE_W, SLIDE_H))


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
