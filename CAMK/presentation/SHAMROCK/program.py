import os

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import pyPLUTO as pp

from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
from matplotlib.ticker import ScalarFormatter
from matplotlib.patches import Wedge

plt.rcParams['axes.prop_cycle'] = plt.cycler(
    color=["steelblue", "olivedrab", "goldenrod", "firebrick", "rebeccapurple"]
)
plt.rcParams['figure.figsize'] = [8, 5]
plt.rcParams['figure.constrained_layout.use'] = True
plt.rcParams['legend.frameon'] = False
plt.rcParams["xtick.minor.visible"] = True
plt.rcParams["ytick.minor.visible"] = True


def _sim_times(D_last, path):
    if 'STAR' in path:
        return D_last.timelist / 62.8318
    elif 'BH' in path:
        return D_last.timelist
    else:
        raise ValueError(path)


def resolve_frame(path, nout=None, time=None):
    if (nout is None) == (time is None):
        raise ValueError

    D_last = pp.Load(nout='last', path=path)
    unique_outs = list(D_last.outlist)
    times = _sim_times(D_last, path)

    if nout is not None:
        if nout not in unique_outs:
            raise ValueError(nout)
        frame = nout
    else:
        idx = int(np.argmin(np.abs(np.asarray(times) - time)))
        frame = unique_outs[idx]

    frame_idx = unique_outs.index(frame)
    frame_time = times[frame_idx]
    return frame, frame_time, D_last, times


def _grid_xz(D_last):
    R, Theta = np.meshgrid(D_last.x1, D_last.x2, indexing='ij')
    X = R * np.sin(Theta)
    Z = R * np.cos(Theta)
    return R, Theta, X, Z


def _mask(X, Z):
    return (X >= 0) & (X <= 21) & (Z >= 0) & (Z <= 11)


def _add_center_patch(ax, path):
    if 'STAR' in path:
        ax.add_patch(Wedge((0.0, 0.0), 1.0, 0.0, 90, facecolor='w', edgecolor=None))
    elif 'BH' in path:
        ax.add_patch(Wedge((0.0, 0.0), 2.1, 0.0, 90, facecolor='k', edgecolor=None))


def _time_label(path, t):
    return f'{t:.1f}' if 'STAR' in path else f'{t:.0f}'


def _tag(path):
    return path.rstrip('/').split('/')[-1]


def _frame_stub(path, frame, t):
    return f'{_tag(path)}_n{frame:04d}_t{t:.0f}'


def _save_jpg(fig, outdir, filename, dpi):
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, filename)
    fig.savefig(outpath, format='jpg', dpi=dpi, bbox_inches='tight',
                pil_kwargs={'quality': 95,
                            'optimize': True,
                            'progressive': True
                           })
    plt.close(fig)
    return outpath


def compute_flux_function(D):
    Br = D.Bx1
    theta = D.x2
    r = D.x1
    integrand = Br * (r[:, None] ** 2) * np.sin(theta)[None, :]
    Psi = np.zeros_like(integrand)
    dtheta = np.diff(theta)
    Psi[:, 1:] = np.cumsum(
        0.5 * (integrand[:, 1:] + integrand[:, :-1]) * dtheta[None, :],
        axis=1
    )
    return Psi


def plot_density(path, nout=None, time=None, dpi=300, outdir='content/'):
    frame, frame_time, D_last, _ = resolve_frame(path, nout=nout, time=time)
    Di = pp.Load(nout=frame, path=path)

    _, _, X, Z = _grid_xz(D_last)
    mask = _mask(X, Z)
    rho = np.where(mask, Di.rho, np.nan)

    fig, ax = plt.subplots(figsize=[6, 8])
    ax.set_facecolor('k')

    im = ax.pcolormesh(
        X, Z, rho,
        cmap='magma',
        norm=LogNorm(vmin=np.nanmin(rho), vmax=np.nanmax(rho)),
        shading='auto'
    )

    ax.set_xlabel('R')
    ax.set_ylabel('z')
    ax.set_aspect('equal')
    ax.set_xlim(0, 20)
    ax.set_ylim(0, 10)

    _add_center_patch(ax, path)
    ax.plot([], [], ' ', label=f't = {_time_label(path, frame_time)}')
    ax.legend(loc='upper right', labelcolor='w')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    fig.colorbar(im, cax=cax, label=r'$\rho$')

    return _save_jpg(fig, outdir, f'density_{_frame_stub(path, frame, frame_time)}.jpg', dpi)


def plot_velocity(path, nout=None, time=None, dpi=300, outdir='content/',
                   stride1=1, stride2=1):
    frame, frame_time, D_last, _ = resolve_frame(path, nout=nout, time=time)
    Di = pp.Load(nout=frame, path=path)

    _, Theta, X, Z = _grid_xz(D_last)
    mask = _mask(X, Z)

    rho = np.where(mask, Di.rho, np.nan)

    vr = Di.vx1
    vth = Di.vx2
    vx = vr * np.sin(Theta) + vth * np.cos(Theta)
    vz = vr * np.cos(Theta) - vth * np.sin(Theta)

    vx = np.where(mask, vx, np.nan)
    vz = np.where(mask, vz, np.nan)
    vr_masked = np.where(mask, vr, np.nan)

    fig, ax = plt.subplots(figsize=[6, 8])
    ax.set_facecolor('k')
    im = ax.pcolormesh(
        X, Z, rho,
        cmap='gray',
        norm=LogNorm(vmin=np.nanmin(rho), vmax=np.nanmax(rho)),
        shading='auto'
    )

    Xq = X[::stride1, ::stride2]
    Zq = Z[::stride1, ::stride2]
    Uq = vx[::stride1, ::stride2]
    Wq = vz[::stride1, ::stride2]
    Vrq = vr_masked[::stride1, ::stride2]

    colors = np.where(Vrq >= 0, 'r', 'b').ravel()

    ax.quiver(
        Xq, Zq, Uq, Wq,
        color=colors,
        scale_units='xy',
        angles='xy',
    )
    ax.set_xlabel('R')
    ax.set_ylabel('z')
    ax.set_aspect('equal')
    ax.set_xlim(0, 20)
    ax.set_ylim(0, 10)

    _add_center_patch(ax, path)
    ax.plot([], [], ' ', label=f't = {_time_label(path, frame_time)}')
    ax.legend(loc='upper right', labelcolor='w')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    fig.colorbar(im, cax=cax, label=r'$\rho$')

    return _save_jpg(fig, outdir, f'velocity_{_frame_stub(path, frame, frame_time)}.jpg', dpi)


def plot_field(path, nout=None, time=None, dpi=300, outdir='content/'):
    frame, frame_time, D_last, _ = resolve_frame(path, nout=nout, time=time)
    Di = pp.Load(nout=frame, path=path)

    _, _, X, Z = _grid_xz(D_last)
    mask = _mask(X, Z)

    rho = np.where(mask, Di.rho, np.nan)
    Psi = compute_flux_function(Di)

    absmax = np.nanmax(np.abs(Psi))
    levels_pos = np.linspace(0, absmax, 100)
    levels_neg = -levels_pos[::-1]

    fig, ax = plt.subplots(figsize=[6, 8])
    ax.set_facecolor('k')

    im = ax.pcolormesh(
        X, Z, rho,
        cmap='gray',
        norm=LogNorm(vmin=np.nanmin(rho), vmax=np.nanmax(rho)),
        shading='auto'
    )

    ax.contour(X, Z, Psi, levels=levels_pos, colors='m', linestyle='-', linewidths=0.4)
    ax.contour(X, Z, Psi, levels=levels_neg, colors='m', linestyle='-', linewidths=0.4)

    ax.set_xlabel('R')
    ax.set_ylabel('z')
    ax.set_aspect('equal')
    ax.set_xlim(0, 20)
    ax.set_ylim(0, 10)

    _add_center_patch(ax, path)
    ax.plot([], [], ' ', label=f't = {_time_label(path, frame_time)}')
    ax.legend(loc='upper right', labelcolor='w')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    fig.colorbar(im, cax=cax, label=r'$\rho$')

    return _save_jpg(fig, outdir, f'field_{_frame_stub(path, frame, frame_time)}.jpg', dpi)


def plot_tracer(path, nout=None, time=None, dpi=300, outdir='content/',
                 mode='truth', compare=False, panel=False):
    frame, frame_time, D_last, _ = resolve_frame(path, nout=nout, time=time)
    Di = pp.Load(nout=frame, path=path)

    if panel:
        abs_diff = np.abs(Di.diskfrac - Di.tr1)
        n_theta = abs_diff.shape[1]
        disagreement = np.sum(abs_diff, axis=1) / n_theta

        x1 = D_last.x1
        n1 = D_last.rho.shape[0]
        disagreement_pct = disagreement * 100.0
        total_disagreement_pct = (np.sum(disagreement) / n1) * 100.0
        total_agreement_pct = 100.0 - total_disagreement_pct

        fig, ax = plt.subplots(figsize=[7, 4])
        ax.fill_between(x1, disagreement_pct, -5, color='m', alpha=1 / 4, linewidth=0)
        ax.set_xscale('log')
        ax.set_xticks([2, 4, 6, 8, 10, 20, 30])
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.minorticks_off()
        ax.xaxis.set_minor_formatter(plt.NullFormatter())
        ax.set_yticks([0, 5, 10, 15, 20])
        ax.set_xlabel('R')
        ax.set_ylabel('|A – B| %')
        ax.set_ylim(-5, 20)

        ax.plot([], [], ' ', label=f't = {_time_label(path, frame_time)}')
        ax.legend(loc='upper right')

        ax.text(0.5, 0.05, f'1 – |A – B| = {total_agreement_pct:.1f} %',
                 transform=ax.transAxes, ha='center', va='bottom', fontsize=10)

        return _save_jpg(fig, outdir, f'tracer_comparison_{_frame_stub(path, frame, frame_time)}.jpg', dpi)

    _, _, X, Z = _grid_xz(D_last)
    mask = _mask(X, Z)

    if compare:
        field = Di.diskfrac - Di.tr1
        vmin, vmax = -1, 1
        is_diff = True
        out_filename = f'tracer_diff_{_frame_stub(path, frame, frame_time)}.jpg'
    else:
        field = Di.diskfrac if mode == 'recover' else Di.tr1
        vmin, vmax = 0, 1
        is_diff = False
        out_filename = f'tracer_{mode}_{_frame_stub(path, frame, frame_time)}.jpg'

    field = np.where(mask, field, np.nan)

    fig, ax = plt.subplots(figsize=[6, 8])
    ax.set_facecolor('k')

    im = ax.pcolormesh(
        X, Z, field,
        cmap='coolwarm',
        vmin=vmin, vmax=vmax,
        shading='auto'
    )

    ax.set_xlabel('R')
    ax.set_ylabel('z')
    ax.set_aspect('equal')
    ax.set_xlim(0, 20)
    ax.set_ylim(0, 10)

    _add_center_patch(ax, path)
    ax.plot([], [], ' ', label=f't = {_time_label(path, frame_time)}')
    ax.legend(loc='upper right', labelcolor='w')

    divider = make_axes_locatable(ax)
    if is_diff:
        cax = divider.append_axes("right", size="5%", pad=0.1)
    else:
        cax = divider.append_axes("right", size="5%", pad=0.4)
    cbar = fig.colorbar(im, cax=cax)

    if is_diff:
        cbar.set_label('A – B')
    else:
        cbar.set_ticks([])
        cax.text(0.5, 1.02, 'disk', transform=cax.transAxes,
                  ha='center', va='bottom', color='k', fontsize=10)
        cax.text(0.5, -0.02, 'corona', transform=cax.transAxes,
                  ha='center', va='top', color='k', fontsize=10)

    return _save_jpg(fig, outdir, out_filename, dpi)
