"""
galaxy_zoo_sph.py

Self-contained generator for a small "zoo" of toy galaxies with different
morphologies (spirals, barred spiral, lenticular, ellipticals, an
interacting pair, a merger remnant, and a dwarf irregular), rendered as
SPH-style particle datasets via yt, plus a uniform grid and a 3D pseudo-AMR
representation of the same particles.

Physics summary
----------------
Disk-type galaxies (spirals / barred spiral / lenticular / dwarf) are built
by integrating test-particle orbits in a 2D gravitational potential:
    - an axisymmetric term giving a flat rotation curve
    - zero or more rotating Fourier ("spiral" or "bar") perturbation terms
Spiral arms / bars are NOT painted onto the density -- they emerge from
genuine orbit crowding as particles respond to the (weak, rotating)
non-axisymmetric potential term. See the "2D vs 3D" note at the bottom of
this file for exactly what is and isn't simulated in 3D.

Elliptical galaxies are not integrated dynamically (they are pressure- /
dispersion-supported, not organized by a coherent rotating perturbation);
instead their particles are sampled directly from a triaxial Plummer-sphere
equilibrium density profile.

The interacting pair and merger remnant use a simplified (Toomre & Toomre
-style) restricted encounter: two galaxy centers move under their own
mutual Plummer-softened gravity (a genuine, if tiny, 2-body integration),
and each galaxy's disk particles feel the sum of their own galaxy's
potential (centered on that moving center) plus a Plummer tidal term from
the other galaxy. Tidal bridges/tails emerge from this, they are not
hand-drawn.

The dwarf irregular additionally multiplies the dynamically-computed
density by a fractal (multi-octave smoothed noise) texture, purely as a
post-hoc visual "clumpiness" layer -- this part is NOT dynamical.

Pseudo-AMR (3D): the grid is a true 3D octree-style patch hierarchy with
cubic cells. Only blocks outside a central ellipsoid stay at the coarsest level;
everything inside it is refined at least once, and deeper levels follow an
oblate-ellipsoid region around each nucleus OR a density threshold, so
spiral arms / bars / tidal features keep their resolution. The vertical
extent of the refinement follows the disk scale height rather than a
sphere. One extra central level (finer than the uniform grid by ~20x)
resolves each nucleus via a mass-conserving prolongation informed by the
actual particles.

All three representations (SPH, uniform grid, pseudo-AMR) are built from
the same 3D particle set in the same box and are rendered face-on and
edge-on; the script writes two montages (zoo_montage_faceon.png,
zoo_montage_edgeon.png). The particle simulation is cached in a pickle so
re-rendering never re-runs the integrators.

Requires: numpy, scipy, yt.
"""

import argparse
import pickle
import zlib

import numpy as np
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter, maximum_filter
import yt

yt.set_log_level(50)

rng = np.random.default_rng(7)

# Physical unit calibration, shared by every representation (SPH, uniform
# grid, AMR) so "the same galaxy" means the same physical scale everywhere.
# These are arbitrary but astrophysically reasonable choices: 1 code length
# unit = 10 kpc (so a disk with rd=0.35 code units is a 3.5 kpc scale
# length, and the box radius of 1.0 code unit is 10 kpc -- both sensible
# for a spiral disk), and each simulation particle stands in for a clump of
# stars/gas totalling 1e6 Msun (so a 60,000-particle galaxy has a total
# mass of 6e10 Msun, also a sensible spiral-galaxy mass).
LENGTH_SCALE_KPC = 10.0
MASS_PER_PARTICLE_MSUN = 1.0e6


# =====================================================================
# Fractal noise texture (used for the dwarf-irregular's clumpiness and as
# a mild texture on the grid-based fields)
# =====================================================================

def fractal_noise(shape, rng, octaves=4, base_sigma=6.0, persistence=0.6):
    """Cheap multi-octave 'value noise': white noise Gaussian-filtered at a
    ladder of scales and summed. A lightweight stand-in for Perlin/simplex
    noise that needs nothing beyond numpy + scipy. Works for any ndim."""
    noise = np.zeros(shape)
    amp = 1.0
    amp_sum = 0.0
    sigma = base_sigma
    for _ in range(octaves):
        layer = rng.normal(size=shape)
        layer = gaussian_filter(layer, sigma=sigma)
        layer /= layer.std() + 1e-12
        noise += amp * layer
        amp_sum += amp
        amp *= persistence
        sigma /= 2.0
    return noise / amp_sum


# =====================================================================
# Disk galaxies: axisymmetric potential + rotating Fourier perturbations
# =====================================================================

def make_components(specs):
    """specs: list of dicts with keys m, pitch_deg, amp, r_peak, r_width, omega_p"""
    comps = []
    for s in specs:
        comps.append(dict(
            m=s["m"],
            tan_pitch=np.tan(np.deg2rad(s["pitch_deg"])),
            amp=s["amp"],
            r_peak=s["r_peak"],
            r_width=s["r_width"],
            omega_p=s["omega_p"],
        ))
    return comps


def phi_axisym(r, V0, RC):
    return 0.5 * V0**2 * np.log(r**2 + RC**2)


def phi_components(r, theta, t, comps, V0, r0_ref=0.08):
    phi = 0.0
    r_eff = np.maximum(r, r0_ref)
    for c in comps:
        envelope = c["amp"] * np.exp(-((r - c["r_peak"]) / c["r_width"]) ** 2)
        phase = (c["m"] / c["tan_pitch"]) * np.log(r_eff / r0_ref)
        phi = phi + (-envelope * V0**2 * np.cos(c["m"] * (theta - c["omega_p"] * t) - phase))
    return phi


def accel_disk(x, y, t, comps, V0, RC, c_other=None, GM_other=0.0, eps_other=0.3):
    """Acceleration for a disk centered at the origin (caller offsets by the
    galaxy's own moving center before calling), optionally with a Plummer
    tidal perturber at c_other=(xo,yo) with mass parameter GM_other."""
    h = 1e-4

    def phi_xy(xx, yy):
        r = np.sqrt(xx**2 + yy**2)
        th = np.arctan2(yy, xx)
        val = phi_axisym(r, V0, RC) + phi_components(r, th, t, comps, V0)
        if c_other is not None:
            dx = xx - c_other[0]
            dy = yy - c_other[1]
            val = val - GM_other / np.sqrt(dx**2 + dy**2 + eps_other**2)
        return val

    dphidx = (phi_xy(x + h, y) - phi_xy(x - h, y)) / (2 * h)
    dphidy = (phi_xy(x, y + h) - phi_xy(x, y - h)) / (2 * h)
    return -dphidx, -dphidy


def sample_exponential_disk(n, rd, r0_max, rng):
    r0 = np.empty(0)
    while r0.shape[0] < n:
        cand = rng.uniform(0.02, r0_max, size=n)
        pdf = cand * np.exp(-cand / rd)
        accept = rng.uniform(0, pdf.max(), size=n) < pdf
        r0 = np.concatenate([r0, cand[accept]])
    return r0[:n]


def simulate_disk_galaxy(
    specs, n=60000, t_total=9.0, dt=0.004, fade_time=2.5,
    V0=1.0, RC=0.08, rd=0.35, r0_max=3.0, r_kill=4.0,
    vdisp=0.02, vphi_scatter=0.04,
):
    """Integrate a disk of n test particles in the axisymmetric + spiral/bar
    potential defined by `specs`. Returns real (dynamically evolved) x, y
    positions (2D only -- see module docstring)."""
    comps = make_components(specs)
    r0 = sample_exponential_disk(n, rd, r0_max, rng)
    theta0 = rng.uniform(0, 2 * np.pi, size=n)
    x = r0 * np.cos(theta0)
    y = r0 * np.sin(theta0)

    vphi = V0 * r0 / np.sqrt(r0**2 + RC**2)
    vphi *= 1.0 + rng.normal(0, vphi_scatter, size=n)
    vx = -vphi * np.sin(theta0) + rng.normal(0, vdisp, size=n)
    vy = vphi * np.cos(theta0) + rng.normal(0, vdisp, size=n)

    nsteps = int(t_total / dt)
    t = 0.0
    ax, ay = accel_disk(x, y, t, comps, V0, RC)
    for i in range(nsteps):
        scale = min(1.0, t / fade_time)
        comps_faded = [dict(c, amp=c["amp"] * scale) for c in comps]
        vx += 0.5 * dt * ax
        vy += 0.5 * dt * ay
        x += dt * vx
        y += dt * vy
        t += dt
        ax, ay = accel_disk(x, y, t, comps_faded, V0, RC)
        vx += 0.5 * dt * ax
        vy += 0.5 * dt * ay

    r_final = np.sqrt(x**2 + y**2)
    mask = r_final < r_kill
    return x[mask], y[mask]


# =====================================================================
# Elliptical galaxies: direct equilibrium (Plummer-sphere) sampling
# =====================================================================

def simulate_elliptical(n=40000, a=0.35, qx=1.0, qy=1.0, qz=0.7, r_cap_factor=6.0, rng=rng):
    """Sample n points from a Plummer-sphere density profile, scaled by axis
    ratios (qx, qy, qz) for a triaxial/flattened shape. This is a genuinely
    3D distribution (unlike the disk galaxies -- see module docstring)."""
    u = rng.uniform(0, 0.995, size=n)
    r = a / np.sqrt(u ** (-2.0 / 3.0) - 1.0)
    r = np.minimum(r, a * r_cap_factor)
    costh = rng.uniform(-1, 1, size=n)
    sinth = np.sqrt(1 - costh**2)
    phi = rng.uniform(0, 2 * np.pi, size=n)
    x = r * sinth * np.cos(phi) * qx
    y = r * sinth * np.sin(phi) * qy
    z = r * costh * qz
    return x, y, z


# =====================================================================
# Mergers: two disks on a mutual-gravity (Plummer-softened) encounter
# =====================================================================

def two_body_trajectory(GM_tot, r12_init, v12_init, t_total, dt, eps=0.3):
    """Integrate the relative separation r12 = c2 - c1 of two point masses
    (Plummer-softened) under their mutual gravity. Returns the trajectory
    array (nsteps+1, 2)."""
    nsteps = int(t_total / dt)
    r = np.array(r12_init, dtype=float)
    v = np.array(v12_init, dtype=float)

    def acc(rr):
        d = np.linalg.norm(rr)
        return -GM_tot * rr / (d**2 + eps**2) ** 1.5

    traj = np.zeros((nsteps + 1, 2))
    traj[0] = r
    a = acc(r)
    for i in range(nsteps):
        v = v + 0.5 * dt * a
        r = r + dt * v
        a = acc(r)
        v = v + 0.5 * dt * a
        traj[i + 1] = r
    return traj


def simulate_merger(
    specs1, specs2, n1=20000, n2=20000, t_total=9.0, dt=0.004, fade_time=1.0,
    V0=1.0, RC=0.08, rd1=0.3, rd2=0.3, r0_max1=1.6, r0_max2=1.6, r_kill=8.0,
    GM_tot=3.0, r12_init=(-4.5, 1.0), v12_init=(0.95, -0.05), eps_gal=0.35,
    m1_frac=0.5, spin1=1.0, spin2=1.0,
):
    """Two disk galaxies (specs1, specs2), each with its own potential,
    size (rd1/rd2, r0_max1/r0_max2), particle count (n1/n2), and spin sense
    (spin1/spin2 = +-1), centered on a moving center -- interacting via
    mutual gravity (their centers) plus a Plummer tidal term (felt by each
    other's disk particles). Different sizes/shapes/spins make the two
    galaxies look like genuinely different objects rather than mirror
    images of each other, and prograde vs retrograde spin (relative to the
    orbit) gives the classic asymmetric tidal response real interacting
    pairs show -- one galaxy raises a much stronger tail than the other.
    Returns the combined (x, y) of surviving particles from both galaxies.
    """
    comps1 = make_components(specs1)
    comps2 = make_components(specs2)
    GM1 = GM_tot * m1_frac
    GM2 = GM_tot * (1 - m1_frac)

    traj = two_body_trajectory(GM_tot, r12_init, v12_init, t_total, dt, eps=eps_gal)
    c1_traj = -(1 - m1_frac) * traj  # c1 = COM - (M2/Mtot)*r12 ; COM fixed at origin
    c2_traj = m1_frac * traj         # c2 = COM + (M1/Mtot)*r12

    v12_now = np.array(v12_init, dtype=float)
    v_c1_init = -(1 - m1_frac) * v12_now
    v_c2_init = m1_frac * v12_now

    def make_disk(n, rd_, r0max_, c_init, v_c_init, spin):
        r0 = sample_exponential_disk(n, rd_, r0max_, rng)
        theta0 = rng.uniform(0, 2 * np.pi, size=n)
        x = r0 * np.cos(theta0) + c_init[0]
        y = r0 * np.sin(theta0) + c_init[1]
        vphi = spin * V0 * r0 / np.sqrt(r0**2 + RC**2)
        vx = -vphi * np.sin(theta0) + v_c_init[0] + rng.normal(0, 0.02, size=n)
        vy = vphi * np.cos(theta0) + v_c_init[1] + rng.normal(0, 0.02, size=n)
        return x, y, vx, vy

    x1, y1, vx1, vy1 = make_disk(n1, rd1, r0_max1, c1_traj[0], v_c1_init, spin1)
    x2, y2, vx2, vy2 = make_disk(n2, rd2, r0_max2, c2_traj[0], v_c2_init, spin2)

    x = np.concatenate([x1, x2])
    y = np.concatenate([y1, y2])
    vx = np.concatenate([vx1, vx2])
    vy = np.concatenate([vy1, vy2])
    owner = np.concatenate([np.zeros(n1, dtype=int), np.ones(n2, dtype=int)])

    nsteps = int(t_total / dt)

    def accel_all(x, y, t, step_idx):
        c1 = c1_traj[step_idx]
        c2 = c2_traj[step_idx]
        scale = min(1.0, t / fade_time)
        comps1_f = [dict(c, amp=c["amp"] * scale) for c in comps1]
        comps2_f = [dict(c, amp=c["amp"] * scale) for c in comps2]

        ax = np.empty_like(x)
        ay = np.empty_like(y)

        sel1 = owner == 0
        sel2 = owner == 1
        ax1, ay1 = accel_disk(
            x[sel1] - c1[0], y[sel1] - c1[1], t, comps1_f, V0, RC,
            c_other=(c2[0] - c1[0], c2[1] - c1[1]), GM_other=GM2, eps_other=eps_gal,
        )
        ax2, ay2 = accel_disk(
            x[sel2] - c2[0], y[sel2] - c2[1], t, comps2_f, V0, RC,
            c_other=(c1[0] - c2[0], c1[1] - c2[1]), GM_other=GM1, eps_other=eps_gal,
        )
        ax[sel1], ay[sel1] = ax1, ay1
        ax[sel2], ay[sel2] = ax2, ay2
        return ax, ay

    t = 0.0
    ax, ay = accel_all(x, y, t, 0)
    for i in range(nsteps):
        vx += 0.5 * dt * ax
        vy += 0.5 * dt * ay
        x += dt * vx
        y += dt * vy
        t += dt
        ax, ay = accel_all(x, y, t, i + 1)
        vx += 0.5 * dt * ax
        vy += 0.5 * dt * ay

    r_final = np.sqrt(x**2 + y**2)
    mask = r_final < r_kill
    return x[mask], y[mask]


# =====================================================================
# SPH rendering: add a vertical coordinate, compute kNN smoothing
# length + density, hand off to yt's generic SPH particle loader
# =====================================================================

H_VERT = 0.035  # vertical scale height used for disk-type galaxies


def add_z_disk(x, y, h=H_VERT, zlim=None, rng=rng):
    """Sample a vertical coordinate from a sech^2(z/h) profile (inverse-CDF
    sampling). This is NOT dynamically coupled to the (x, y) orbit
    integration -- see the "2D vs 3D" note below."""
    n = x.shape[0]
    u = rng.uniform(1e-4, 1 - 1e-4, size=n)
    z = h * np.arctanh(2 * u - 1)
    if zlim is not None:
        z = np.clip(z, zlim[0], zlim[1])
    return z


def particles_to_sph(x, y, z, k=32):
    """Real 3D kNN smoothing length + SPH density estimate from the full
    3D particle positions."""
    pos = np.column_stack([x, y, z])
    tree = cKDTree(pos)
    dists, _ = tree.query(pos, k=min(k, pos.shape[0] - 1))
    hsml = np.maximum(dists[:, -1], 1e-4)
    mass = np.ones(x.shape[0])
    density = (k * 1.0) / ((4.0 / 3.0) * np.pi * hsml**3)
    return pos, mass, hsml, density


def z_half_extent(hw, spheroid):
    """Half-thickness of the vertical box, shared by ALL three
    representations: always a cube (same as the in-plane half-width), for
    disks and spheroids alike, so face-on and edge-on views have the same
    layout. (`spheroid` is kept for API compatibility.)"""
    return hw


DENS = ("gas", "density")


def build_sph_ds(x, y, z, xlim, ylim, zlim, weight=None):
    """Build a yt SPH particle dataset.

    x, y, z may extend beyond (xlim, ylim, zlim) -- e.g. the disk galaxies
    are seeded out to a radius well past the visible box (see
    simulate_disk_galaxy's r0_max), so that particles near the edge of the
    plotted view still have real neighbors just outside it. Smoothing
    length / density are computed from the FULL particle set for exactly
    that reason; only the subset that actually falls inside (xlim, ylim,
    zlim) is then handed to yt (which requires every particle to lie
    within the declared bounding box).
    """
    pos, mass, hsml, density = particles_to_sph(x, y, z)
    if weight is not None:
        density = density * weight

    inside = (
        (pos[:, 0] > xlim[0]) & (pos[:, 0] < xlim[1]) &
        (pos[:, 1] > ylim[0]) & (pos[:, 1] < ylim[1]) &
        (pos[:, 2] > zlim[0]) & (pos[:, 2] < zlim[1])
    )
    pos, mass, hsml, density = pos[inside], mass[inside], hsml[inside], density[inside]

    data = {
        "particle_position_x": pos[:, 0],
        "particle_position_y": pos[:, 1],
        "particle_position_z": pos[:, 2],
        "particle_mass": mass,
        "density": density,
        "smoothing_length": hsml,
    }
    bbox = np.array([xlim, ylim, zlim])
    return yt.load_particles(
        data, bbox=bbox,
        length_unit=(LENGTH_SCALE_KPC, "kpc"),
        mass_unit=(MASS_PER_PARTICLE_MSUN, "Msun"),
    )


# =====================================================================
# Shared machinery for the two grid-based representations (uniform grid
# and pseudo-AMR): a single smoothed+denoised "master" density field,
# converted to physical units via exact block-summation so that a coarser
# cell's value is always the true sum of the finer cells it replaces
# (mass-conserving downsampling, not just a resample).
# =====================================================================

def deposit_adaptive(x, y, z, hsml, xlim, ylim, zlim, shape, weights=None,
                     sigma_frac=0.5, min_sigma=0.7, bins_per_octave=2):
    """Adaptive-kernel deposit: every particle is spread as an isotropic
    Gaussian with sigma = sigma_frac * hsml (its own kNN smoothing length,
    the SAME length the SPH representation uses), so grid fields inherit
    SPH's resolution: sharp where particles are dense, smooth where they
    are sparse. Done efficiently by grouping particles into log-spaced
    sigma bins, histogramming each group and filtering it once."""
    cell = (xlim[1] - xlim[0]) / shape[0]
    sig = np.maximum(sigma_frac * hsml / cell, min_sigma)  # in cells
    key = np.round(np.log2(sig) * bins_per_octave).astype(int)
    field = np.zeros(shape)
    for kk in np.unique(key):
        m = key == kk
        Hk, _ = np.histogramdd(np.column_stack([x[m], y[m], z[m]]), bins=shape,
                               range=[xlim, ylim, zlim],
                               weights=None if weights is None else weights[m])
        field += gaussian_filter(Hk, sigma=2.0 ** (kk / bins_per_octave), truncate=3.0)
    return field


def compute_smooth_field_3d(
    x, y, z, xlim, ylim, zlim, shape, weights=None, sigma_fine=2.0,
    sigma_coarse_frac=0.14, coarse_weight=0.6, z_sigma_coarse=None,
    noise_amp=0.15, noise_octaves=4, noise_base_sigma=8.0, rng=rng, hsml=None,
):
    """Real (x, y, z) particle histogram -> blended fine+coarse Gaussian
    smoothing -> correlated noise -> renormalisation. Returns a 3D field
    whose SUM equals the (weighted) number of particles in the box, so
    that sums over any region are a mass-like proxy that later gets
    divided by a cell's physical volume to get a real density.

    The blended smoothing is what keeps sparse regions from looking like
    scattered dots: a narrow pass preserves real structure (arms, bar,
    tidal features) while a much wider pass (a sizeable fraction of the
    box) turns sparse outskirts into a smooth, continuous background.
    Applied isotropically in a thin disk the wide pass would puff the
    disk up vertically, so z_sigma_coarse (in cells) caps its vertical
    width (pass ~ one scale height for disks; None = isotropic, which is
    right for spheroids). The final rescale removes the ~1.6x extra mass
    that summing the two smoothing passes (weight 1 + coarse_weight)
    would otherwise inject."""
    if hsml is not None:
        # adaptive (SPH-matched) smoothing: replaces the fixed fine+coarse blend
        # and the noise texture, which would wash out inter-arm contrast
        field = deposit_adaptive(x, y, z, hsml, xlim, ylim, zlim, shape, weights)
        tot = field.sum()
        field = field + field.max() * 3e-4
        return field * (tot / field.sum())
    H, _ = np.histogramdd(np.column_stack([x, y, z]), bins=shape,
                          range=[xlim, ylim, zlim], weights=weights)
    sc_xy = sigma_coarse_frac * shape[0]
    sc_z = sc_xy if z_sigma_coarse is None else min(sc_xy, z_sigma_coarse)
    H_fine = gaussian_filter(H, sigma=sigma_fine)
    H_coarse = gaussian_filter(H, sigma=(sc_xy, sc_xy, sc_z))
    field = H_fine + coarse_weight * H_coarse
    del H_fine, H_coarse

    noise = fractal_noise(shape, rng, octaves=noise_octaves,
                           base_sigma=noise_base_sigma, persistence=0.55)
    field = field * np.exp(noise_amp * noise)
    del noise
    field = np.clip(field, 0, None)
    field = field + field.max() * 3e-4
    return field * (H.sum() / field.sum())


def block_reduce_sum(field, factor):
    """Downsample an N-d field by an integer factor, summing each
    factor^ndim block. Requires every dimension to be exactly divisible by
    factor (guaranteed by construction everywhere this is used)."""
    shp = []
    for n in field.shape:
        shp += [n // factor, factor]
    return field.reshape(shp).sum(axis=tuple(range(1, 2 * field.ndim, 2)))


def counts_to_density_3d(count_field3d, left_edge, right_edge):
    """Cell 'particle count' -> physical density: cell mass / cell volume,
    in Msun/kpc**3."""
    nx, ny, nz = count_field3d.shape
    vol_kpc3 = 1.0
    for n, lo, hi in zip((nx, ny, nz), left_edge, right_edge):
        vol_kpc3 *= (hi - lo) / n * LENGTH_SCALE_KPC
    return count_field3d * MASS_PER_PARTICLE_MSUN / vol_kpc3


# =====================================================================
# Fixed-resolution uniform grid (NOT adaptive/AMR) built from the same
# particles, for a "what would a different, grid-based code's snapshot
# of this same galaxy look like" comparison. Genuinely 3D with cubic
# cells, so ellipticals are round and disks have a real vertical profile
# (taken from the particles' z, not an imposed analytic sech^2).
# =====================================================================

def build_uniform_ds(
    x, y, z, xlim, ylim, zlim, n_cells=56, supersample=4, h_vert=H_VERT,
    spheroid=False, weights=None, rng=rng, hsml=None, **smooth_kwargs,
):
    """Bin particles onto a single fixed-resolution uniform grid (no
    refinement at all -- every cell is the same size, and deliberately
    modest in count so the cells themselves are visible, like a real
    coarse-ish grid-code run). n_cells is the count across x and y; the
    vertical count follows from cubic cells. The field is computed smooth
    at a finer working resolution first (supersample x) and then block-
    summed down, keeping the smooth, non-dotty properties while still
    giving a visibly blocky final image."""
    dx = (xlim[1] - xlim[0]) / n_cells
    nz = int(round((zlim[1] - zlim[0]) / dx))
    shape = (n_cells * supersample, n_cells * supersample, nz * supersample)
    z_sigma = None if spheroid else h_vert / (dx / supersample)
    fine = compute_smooth_field_3d(x, y, z, xlim, ylim, zlim, shape, weights=weights,
                                   z_sigma_coarse=z_sigma, rng=rng, hsml=hsml, **smooth_kwargs)
    counts = block_reduce_sum(fine, supersample)
    dens = counts_to_density_3d(counts, [xlim[0], ylim[0], zlim[0]], [xlim[1], ylim[1], zlim[1]])
    bbox = np.array([xlim, ylim, list(zlim)])
    return yt.load_uniform_grid(
        {"density": (dens, "Msun/kpc**3")}, dens.shape, bbox=bbox,
        length_unit=(LENGTH_SCALE_KPC, "kpc"),
    )


# =====================================================================
# yt output helpers: labelled PNGs (optional) and raw slice arrays for
# the montages. Face-on = slice normal to z; edge-on = slice normal to y
# (through the midplane).
# =====================================================================

def _style(p, title, unit="Msun/kpc**3", grids=False):
    p.set_unit(DENS, unit)
    p.set_axes_unit("kpc")
    p.set_cmap(DENS, "viridis")
    p.annotate_title(title)
    if grids:
        p.annotate_grids(periodic=False)


def save_yt_plots(ds, base, title, grids=False):
    """Write {base}_faceon.png and {base}_edgeon.png (labelled yt slices)."""
    for axis, tag in (("z", "faceon"), ("y", "edgeon")):
        p = yt.SlicePlot(ds, axis, DENS)
        _style(p, title, grids=grids)
        p.save(f"{base}_{tag}.png")


def slice_arrays(ds, hw, zhalf, res=512):
    """Midplane slices as plain arrays in Msun/kpc**3, oriented for
    imshow(origin='lower'): face-on has x horizontal / y vertical; edge-on
    has x horizontal / z vertical (galaxy lying flat). Face-on is res x res;
    edge-on is res wide and res*zhalf/hw tall (cubic pixels)."""
    q = lambda v: ds.quan(v, "code_length")
    nz_res = max(8, int(round(res * zhalf / hw)))
    face = ds.slice("z", 0.0).to_frb(q(2 * hw), res)[DENS].to("Msun/kpc**3").d
    # yt's image axes for a y-normal slice are (horizontal = z, vertical = x),
    # and frb arrays are indexed [vertical, horizontal]; transpose so the
    # galaxy lies flat (x horizontal, z vertical).
    edge = ds.slice("y", 0.0).to_frb(q(2 * zhalf), (nz_res, res), height=q(2 * hw))[DENS].to("Msun/kpc**3").d
    return np.ascontiguousarray(face), np.ascontiguousarray(edge.T)

# =====================================================================
# 3D pseudo-AMR: cubic cells, refine_by=2 in x, y AND z.
#
# Refinement rules (a block at level L is split into level L+1 if...):
#   L = 0 : its center lies inside an ellipsoid around the box center.
#           Only blocks OUTSIDE that ellipsoid stay at the coarsest level,
#           so the coarse region is an ellipsoidal shell, not a box frame.
#   L >= 1: its center lies inside an OBLATE ELLIPSOID around any
#           identified nucleus (semi-axes R_L in the disk plane and
#           Z_L = min(R_L, z_extent_h * h_vert) vertically), OR its mean
#           density exceeds density_frac[L-1] * (peak fine-cell density).
#   L = max_level : (central level) inside a small oblate ellipsoid around
#           each nucleus, split once more by conservative prolongation.
#
# Why not a sphere? A disk with scale height h << R is ~30x thinner than it
# is wide. A sphere of radius R would spend almost all of its fine cells on
# empty space above/below the disk. The oblate region hugs the disk, and
# the density criterion follows arms / bars / tidal tails wherever they
# are, so detail is kept where there is actually structure. Spheroidal
# galaxies (ellipticals) set spherical=True and get a sphere in a cubic box.
# =====================================================================

def find_centers(
    x, y, xlim, ylim, n_centers=1, resolution=48, smooth_sigma=3.0,
    min_sep_frac=0.2, min_prominence=0.08,
):
    """Locate up to n_centers density peaks (e.g. galaxy nuclei) from the
    particle distribution.

    Two things make this robust against the obvious failure modes:
      - a minimum peak PROMINENCE (relative to the global max) so small
        noise-driven local maxima never count as a peak
      - a minimum SEPARATION between accepted centers, scaled to the box
        size, so a galaxy with a clumpy/double nucleus (e.g. the twin
        star-forming knots at the center of the barred spiral) is never
        mistaken for two different galaxies -- the separation floor is
        comfortably larger than any such sub-galaxy clump spacing, but
        well below the actual galaxy-galaxy separation in the merger
        cases, so real pairs are still found as two centers.
    Peaks are returned strongest-first, deduplicated by the separation
    rule, so the result is always well-separated and never contains two
    (near-)identical centers -- which matters downstream, since a
    pseudo-AMR refinement region built around two coincident centers
    could otherwise produce degenerate, near-zero-area patches that
    break yt's slicing/ray logic.
    """
    H, _, _ = np.histogram2d(x, y, bins=resolution, range=[xlim, ylim])
    Hs = gaussian_filter(H, sigma=smooth_sigma)

    neighborhood = maximum_filter(Hs, size=5, mode="nearest")
    is_local_max = (Hs == neighborhood) & (Hs > Hs.max() * min_prominence)
    coords = np.argwhere(is_local_max)
    vals = Hs[is_local_max]
    order = np.argsort(-vals)
    coords = coords[order]

    dx = (xlim[1] - xlim[0]) / resolution
    dy = (ylim[1] - ylim[0]) / resolution
    min_sep = min_sep_frac * (xlim[1] - xlim[0])

    centers = []
    for (i, j) in coords:
        px = xlim[0] + (i + 0.5) * dx
        py = ylim[0] + (j + 0.5) * dy
        if all(np.hypot(px - cx, py - cy) > min_sep for cx, cy in centers):
            centers.append((px, py))
        if len(centers) >= n_centers:
            break

    if not centers:
        i, j = np.unravel_index(np.argmax(Hs), Hs.shape)
        centers = [(xlim[0] + (i + 0.5) * dx, ylim[0] + (j + 0.5) * dy)]

    return centers


def build_amr_patches_3d(
    fine_field, xlim, ylim, zlim, centers, h_vert, particles=None, base_res=16,
    base_nz=4, block=2, block0=1, refine_by=2, max_level=4, envelope_axes_frac=None,
    radii_frac=(0.55, 0.32, 0.16), density_frac=(0.02, 0.05, 0.12),
    z_extent_h=3.0, spherical=False, central_radius_frac=0.07, central_sigma=1.5,
):
    """Build the list of 3D patches (root + nested refined blocks).

    fine_field must have shape (N, N, Nz) with N = base_res * refine_by**max_level
    and Nz = base_nz * refine_by**max_level, and the cells must be cubic
    (zlim spans exactly base_nz root cells). Patches for levels 0..max_level
    are carved out of fine_field by exact integer block-summation, in
    pixel-index space throughout the recursion, so every patch lands exactly
    on its parent's cell edges (yt's AMR loader requires that) and mass is
    conserved across levels.

    block0 : block size (in root cells) used for the level-0 -> 1 split only.
        Small (1) so the ellipsoid boundary is smooth rather than a coarse staircase.
    envelope_axes_frac : semi-axes (fractions of the half-width) of the
        ELLIPSOID around the box center inside which level-0 blocks are
        refined; blocks outside it stay coarse, so the coarsest level is an
        ellipsoidal shell (elliptical in every coordinate plane), never the
        box faces. Default (0.9, 0.9, 0.9) for spheroids, oblate
        (0.9, 0.9, 0.5) for disks.
    radii_frac / density_frac : per-level (L = 1, 2, 3, ...) nucleus-region
        radius (fraction of half-width) and density threshold (fraction of
        the peak fine-cell value). The last entry is reused for deeper levels.
    z_extent_h : vertical half-thickness of the refinement ellipsoid in
        scale heights (ignored when spherical=True).
    central_radius_frac : if `particles` (an (n, 3) array) is given, ONE
        extra level (max_level + 1) is added inside this (oblate) radius of
        each nucleus. fine_field has no resolution beyond max_level, so the
        extra level is a CONSERVATIVE PROLONGATION: every parent cell's mass
        is split among its children in proportion to a lightly smoothed
        (central_sigma, in child cells) histogram of the actual particles in
        that region. Children therefore sum exactly to their parent (mass
        conserved, no seams between levels) while the sub-cell structure is
        real particle information rather than interpolation. Set to None to
        disable.
    """
    R = refine_by
    master = (base_res * R**max_level, base_res * R**max_level, base_nz * R**max_level)
    assert fine_field.shape == master, f"fine_field must be {master}, got {fine_field.shape}"
    assert base_res % block == 0 and base_nz % block == 0 and base_res % block0 == 0

    dx = (xlim[1] - xlim[0]) / master[0]
    dy = (ylim[1] - ylim[0]) / master[1]
    dz = (zlim[1] - zlim[0]) / master[2]
    assert np.allclose([dx, dy], dz), "cells must be cubic (zlim must span base_nz root cells)"

    hw = 0.5 * (xlim[1] - xlim[0])
    env = hw * np.array(envelope_axes_frac if envelope_axes_frac is not None
                        else ((0.9, 0.9, 0.9) if spherical else (0.9, 0.9, 0.5)))
    box_c = [0.5 * (a + b) for a, b in (xlim, ylim, zlim)]
    ref_value = fine_field.max()
    central = particles is not None and central_radius_frac is not None
    pad = int(np.ceil(4 * central_sigma))  # gaussian_filter truncates at 4 sigma

    patches = [dict(
        left_edge=[xlim[0], ylim[0], zlim[0]], right_edge=[xlim[1], ylim[1], zlim[1]],
        level=0, dims=(base_res, base_res, base_nz),
        density=block_reduce_sum(fine_field, R**max_level),
    )]

    def wants_refinement(level, i0, i1, j0, j1, k0, k1):
        cx = xlim[0] + 0.5 * (i0 + i1) * dx
        cy = ylim[0] + 0.5 * (j0 + j1) * dy
        cz = zlim[0] + 0.5 * (k0 + k1) * dz
        if level == 0:
            return (((cx - box_c[0]) / env[0]) ** 2 + ((cy - box_c[1]) / env[1]) ** 2
                    + ((cz - box_c[2]) / env[2]) ** 2) <= 1.0
        idx = min(level - 1, len(radii_frac) - 1)
        R_L = radii_frac[idx] * hw
        Z_L = R_L if spherical else min(R_L, z_extent_h * h_vert)
        for ccx, ccy in centers:
            if ((cx - ccx) ** 2 + (cy - ccy) ** 2) / R_L**2 + (cz / Z_L) ** 2 <= 1.0:
                return True
        mean_density = fine_field[i0:i1, j0:j1, k0:k1].mean()
        return mean_density > density_frac[min(idx, len(density_frac) - 1)] * ref_value

    R_c = central_radius_frac * hw if central else 0.0
    Z_c = R_c if spherical else min(R_c, z_extent_h * h_vert)

    def in_central(i0, i1, j0, j1, k0, k1):
        cx = xlim[0] + 0.5 * (i0 + i1) * dx
        cy = ylim[0] + 0.5 * (j0 + j1) * dy
        cz = zlim[0] + 0.5 * (k0 + k1) * dz
        return any(((cx - ccx) ** 2 + (cy - ccy) ** 2) / R_c**2 + (cz / Z_c) ** 2 <= 1.0
                   for ccx, ccy in centers)

    def prolong(i0, i1, j0, j1, k0, k1):
        """Split each cell of fine_field[i0:i1, j0:j1, k0:k1] into R^3
        children, weighted by the local particle histogram."""
        parent = fine_field[i0:i1, j0:j1, k0:k1]
        nchild = block * R
        dxf, dyf, dzf = dx / R, dy / R, dz / R
        lo = np.array([xlim[0] + i0 * dx - pad * dxf, ylim[0] + j0 * dy - pad * dyf,
                       zlim[0] + k0 * dz - pad * dzf])
        hi = np.array([xlim[0] + i1 * dx + pad * dxf, ylim[0] + j1 * dy + pad * dyf,
                       zlim[0] + k1 * dz + pad * dzf])
        m = np.all((particles > lo) & (particles < hi), axis=1)
        nb = nchild + 2 * pad
        H, _ = np.histogramdd(particles[m], bins=(nb, nb, nb), range=list(zip(lo, hi)))
        S = gaussian_filter(H, sigma=central_sigma)[pad:-pad, pad:-pad, pad:-pad]
        S = S + 0.1 * S.mean() + 1e-12  # a little floor: empty regions split evenly
        S6 = S.reshape(block, R, block, R, block, R)
        w = S6 / S6.sum(axis=(1, 3, 5), keepdims=True)
        child = w * parent[:, None, :, None, :, None]
        return child.reshape(nchild, nchild, nchild)

    def recurse(level, i0, i1, j0, j1, k0, k1):
        if level > max_level or (level == max_level and not central):
            return
        b = block0 if level == 0 else block
        bpx = b * R ** (max_level - level)  # block size in master pixels
        for pi0 in range(i0, i1, bpx):
            for pj0 in range(j0, j1, bpx):
                for pk0 in range(k0, k1, bpx):
                    pi1, pj1, pk1 = pi0 + bpx, pj0 + bpx, pk0 + bpx
                    if level < max_level:
                        if not wants_refinement(level, pi0, pi1, pj0, pj1, pk0, pk1):
                            continue
                        child_factor = R ** (max_level - level - 1)
                        dchild = block_reduce_sum(fine_field[pi0:pi1, pj0:pj1, pk0:pk1], child_factor)
                    else:  # level == max_level: central prolongation level
                        if not in_central(pi0, pi1, pj0, pj1, pk0, pk1):
                            continue
                        dchild = prolong(pi0, pi1, pj0, pj1, pk0, pk1)
                    patches.append(dict(
                        left_edge=[xlim[0] + pi0 * dx, ylim[0] + pj0 * dy, zlim[0] + pk0 * dz],
                        right_edge=[xlim[0] + pi1 * dx, ylim[0] + pj1 * dy, zlim[0] + pk1 * dz],
                        level=level + 1, dims=dchild.shape, density=dchild,
                    ))
                    if level < max_level:
                        recurse(level + 1, pi0, pi1, pj0, pj1, pk0, pk1)

    recurse(0, 0, master[0], 0, master[1], 0, master[2])
    return patches


def build_amr_ds(
    x, y, z, xlim, ylim, n_centers=1, base_res=16, block=2, refine_by=2,
    max_level=4, h_vert=H_VERT, spherical=False, z_sigma_h=1.0, weights=None,
    rng=rng, hsml=None, **crit_kwargs,
):
    """Deposit particles on a 3D pseudo-AMR hierarchy and return
    (yt dataset, info). Levels 0..max_level come from the smoothed master
    field; one more central level (see build_amr_patches_3d) resolves the
    nuclei beyond it.

    Box: xlim x ylim in the plane; vertically the box is base_nz root cells
    thick (cubic cells): a cube for every galaxy type (see z_half_extent)."""
    R = refine_by
    dx0 = (xlim[1] - xlim[0]) / base_res
    base_nz = base_res
    zhalf = 0.5 * base_nz * dx0
    zlim = (-zhalf, zhalf)
    shape = (base_res * R**max_level, base_res * R**max_level, base_nz * R**max_level)
    dz = (zlim[1] - zlim[0]) / shape[2]
    z_sigma = None if spherical else z_sigma_h * h_vert / dz

    fine_field = compute_smooth_field_3d(
        x, y, z, xlim, ylim, zlim, shape, weights=weights, z_sigma_coarse=z_sigma, rng=rng,
        hsml=hsml)
    centers = find_centers(x, y, xlim, ylim, n_centers=n_centers)

    patches = build_amr_patches_3d(
        fine_field, xlim, ylim, zlim, centers, h_vert,
        particles=np.column_stack([x, y, z]), base_res=base_res,
        base_nz=base_nz, block=block, refine_by=R, max_level=max_level,
        spherical=spherical, **crit_kwargs,
    )
    del fine_field

    grid_data = []
    for p in patches:
        dens = counts_to_density_3d(p["density"], p["left_edge"], p["right_edge"])
        grid_data.append(dict(
            left_edge=p["left_edge"], right_edge=p["right_edge"], level=p["level"],
            dimensions=list(p["dims"]), density=(dens, "Msun/kpc**3"),
        ))

    bbox = np.array([xlim, ylim, list(zlim)])
    ds = yt.load_amr_grids(
        grid_data, [base_res, base_res, base_nz], bbox=bbox, refine_by=R,
        length_unit=(LENGTH_SCALE_KPC, "kpc"),
    )
    counts = np.bincount([p["level"] for p in patches])
    finest_pc = (xlim[1] - xlim[0]) * LENGTH_SCALE_KPC * 1000 / (base_res * R ** (len(counts) - 1))
    return ds, dict(n_grids=len(patches), grids_per_level=counts.tolist(),
                    finest_cell_pc=round(finest_pc, 1))


# =====================================================================
# The 10-galaxy zoo definitions
# =====================================================================

DISK_SPECS = {
    "1_grand_design": dict(
        specs=[dict(m=2, pitch_deg=18, amp=0.55, r_peak=0.5, r_width=1.2, omega_p=0.55)],
        kw=dict(),
    ),
    "2_tightly_wound": dict(
        specs=[dict(m=2, pitch_deg=8, amp=0.55, r_peak=0.5, r_width=1.2, omega_p=0.55)],
        kw=dict(),
    ),
    "3_flocculent": dict(
        specs=[dict(m=6, pitch_deg=25, amp=0.35, r_peak=0.5, r_width=1.0, omega_p=0.55)],
        kw=dict(vdisp=0.05, vphi_scatter=0.08),
    ),
    "4_barred_spiral": dict(
        specs=[
            dict(m=2, pitch_deg=87, amp=0.45, r_peak=0.22, r_width=0.32, omega_p=0.6),
            dict(m=2, pitch_deg=20, amp=0.3, r_peak=0.75, r_width=1.0, omega_p=0.35),
        ],
        kw=dict(fade_time=3.0),
    ),
    "5_lenticular": dict(specs=[], kw=dict()),
}

ELLIPTICALS = {
    "6_elliptical_E0": dict(a=0.3, qx=1.0, qy=1.0, qz=0.8),
    "7_elliptical_E5": dict(a=0.3, qx=1.0, qy=0.45, qz=0.6),
}

MERGER_SPEC_BIG = [dict(m=2, pitch_deg=20, amp=0.4, r_peak=0.4, r_width=0.6, omega_p=0.6)]
MERGER_SPEC_SMALL = [dict(m=4, pitch_deg=15, amp=0.3, r_peak=0.22, r_width=0.3, omega_p=1.0)]
MERGER_SPEC_BARRED_SMALL = [
    dict(m=2, pitch_deg=86, amp=0.4, r_peak=0.12, r_width=0.16, omega_p=1.1),
]

MERGERS = {
    # bigger grand-design primary + smaller, faster-rotating flocculent-ish
    # secondary on a RETROGRADE spin -- prograde/retrograde asymmetry means
    # the two galaxies respond very differently to the same encounter, one
    # raising a much stronger tidal feature than the other (as in real
    # interacting pairs, e.g. the Mice Galaxies)
    "8_interacting_pair": dict(
        specs1=MERGER_SPEC_BIG, specs2=MERGER_SPEC_SMALL,
        rd1=0.32, rd2=0.15, r0_max1=1.5, r0_max2=0.85,
        n1=18000, n2=9000, spin1=1.0, spin2=-1.0, m1_frac=0.65,
        t_total=3.6, GM_tot=4.5, r12_init=(-4.5, 1.2), v12_init=(0.95, -0.75),
        lim=3,
    ),
    # bigger primary + smaller barred secondary, both prograde (drives a
    # deeper, more disruptive encounter toward coalescence)
    "9_merger_remnant": dict(
        specs1=MERGER_SPEC_BIG, specs2=MERGER_SPEC_BARRED_SMALL,
        rd1=0.3, rd2=0.13, r0_max1=1.4, r0_max2=0.7,
        n1=17000, n2=8000, spin1=1.0, spin2=1.0, m1_frac=0.68,
        t_total=19.0, GM_tot=4.0, r12_init=(-4.5, 0.7), v12_init=(0.45, -0.05),
        lim=5,
    ),
}

TITLES = {
    "1_grand_design": "Grand-design spiral",
    "2_tightly_wound": "Tightly-wound spiral",
    "3_flocculent": "Flocculent spiral",
    "4_barred_spiral": "Barred spiral",
    "5_lenticular": "Lenticular (S0)",
    "6_elliptical_E0": "Elliptical E0",
    "7_elliptical_E5": "Elliptical E5",
    "8_interacting_pair": "Interacting pair (bridge)",
    "9_merger_remnant": "Merger remnant (irregular)",
    "10_dwarf_irregular": "Dwarf irregular",
}

ORDER = [
    "1_grand_design", "2_tightly_wound", "3_flocculent", "4_barred_spiral", "5_lenticular",
    "6_elliptical_E0", "7_elliptical_E5", "8_interacting_pair", "9_merger_remnant",
    "10_dwarf_irregular",
]


def generate_all(n_disk=60000, n_ellip=40000, only=None):
    """Run the dynamics/sampling for all 10 galaxies (or just those named in
    `only`). Returns a dict keyed by galaxy name -> dict(x, y, [z], xlim,
    ylim, n_centers, [dwarf])."""
    want = lambda name: only is None or name in only
    results = {}

    for name, cfg in DISK_SPECS.items():
        if not want(name):
            continue
        x, y = simulate_disk_galaxy(cfg["specs"], n=n_disk, t_total=9.0, dt=0.004, **cfg["kw"])
        results[name] = dict(x=x, y=y, xlim=(-1, 1), ylim=(-1, 1), n_centers=1)

    for name, cfg in ELLIPTICALS.items():
        if not want(name):
            continue
        x, y, z = simulate_elliptical(n=n_ellip, **cfg)
        results[name] = dict(x=x, y=y, z=z, xlim=(-1, 1), ylim=(-1, 1), n_centers=1)

    for name, cfg in MERGERS.items():
        if not want(name):
            continue
        cfg = dict(cfg)
        lim = cfg.pop("lim")
        specs1 = cfg.pop("specs1")
        specs2 = cfg.pop("specs2")
        x, y = simulate_merger(specs1, specs2, r_kill=12, **cfg)
        results[name] = dict(x=x, y=y, xlim=(-lim, lim), ylim=(-lim, lim), n_centers=2)

    if want("10_dwarf_irregular"):
        x, y = simulate_disk_galaxy([], n=25000, t_total=9.0, dt=0.004, rd=0.15, r0_max=0.5,
                                     V0=0.5, vdisp=0.05, vphi_scatter=0.15)
        results["10_dwarf_irregular"] = dict(x=x, y=y, xlim=(-1, 1), ylim=(-1, 1),
                                              n_centers=1, dwarf=True)

    return results


METHODS = ("sph", "uniform", "amr")
METHOD_LABELS = {"sph": "SPH", "uniform": "Uniform grid", "amr": "Pseudo-AMR"}


def render_galaxy(name, d, outdir=".", individual=False, grids=False, amr_kwargs=None, res=512):
    """Render one galaxy three ways (SPH, uniform grid, 3D pseudo-AMR), each
    face-on and edge-on. All three share the same particles (x, y from the
    saved simulation, one z per particle) and the same physical box (see
    z_half_extent). Returns {method: {"face": img, "edge": img}, "hw", "zhalf",
    "amr_info"} with images as Msun/kpc**3 midplane-slice arrays.

    Randomness (vertical coordinate, noise textures) is seeded from the
    galaxy name, so re-rendering from the saved particles is reproducible."""
    xlim, ylim = d["xlim"], d["ylim"]
    hw = 0.5 * (xlim[1] - xlim[0])
    spheroid = "z" in d
    zhalf = z_half_extent(hw, spheroid)
    zlim = (-zhalf, zhalf)
    h_vert = H_VERT * hw  # scale vertical thickness with box size
    seed = zlib.crc32(name.encode())
    zrng = np.random.default_rng(seed)

    x, y = d["x"], d["y"]
    z = d["z"] if spheroid else add_z_disk(x, y, h=h_vert, rng=zrng)  # unclipped

    weight = None
    if d.get("dwarf"):  # visual clumpiness layer, applied identically in every method
        nres = 220
        noise = fractal_noise((nres, nres), zrng, octaves=4, base_sigma=6.0, persistence=0.6)
        xi = np.clip(((x - xlim[0]) / (xlim[1] - xlim[0]) * nres).astype(int), 0, nres - 1)
        yi = np.clip(((y - ylim[0]) / (ylim[1] - ylim[0]) * nres).astype(int), 0, nres - 1)
        weight = np.exp(0.9 * noise[xi, yi])

    title = TITLES[name]
    out = dict(hw=hw, zhalf=zhalf)

    ds = build_sph_ds(x, y, z, xlim, ylim, zlim, weight=weight)
    out["sph"] = dict(zip(("face", "edge"), slice_arrays(ds, hw, zhalf, res)))
    if individual:
        save_yt_plots(ds, f"{outdir}/zoo_{name}_sph", title + " (SPH)")
    del ds

    # grid-based representations only need particles inside the visible box
    inside = (x > xlim[0]) & (x < xlim[1]) & (y > ylim[0]) & (y < ylim[1])
    xg, yg, zg = x[inside], y[inside], z[inside]
    wg = None if weight is None else weight[inside]
    hg = particles_to_sph(x, y, z)[2][inside]  # SPH kNN smoothing lengths

    ds = build_uniform_ds(xg, yg, zg, xlim, ylim, zlim, h_vert=h_vert, spheroid=spheroid,
                          weights=wg, hsml=hg, rng=np.random.default_rng(seed + 1))
    out["uniform"] = dict(zip(("face", "edge"), slice_arrays(ds, hw, zhalf, res)))
    if individual:
        save_yt_plots(ds, f"{outdir}/zoo_{name}_uniform", title + " (uniform grid)")
    del ds

    kw = dict(amr_kwargs or {})
    if spheroid:
        kw.setdefault("max_level", 3)  # cubic box: keep the 3D master field small
    ds, info = build_amr_ds(xg, yg, zg, xlim, ylim, n_centers=d.get("n_centers", 1),
                            h_vert=h_vert, spherical=spheroid, weights=wg, hsml=hg,
                            rng=np.random.default_rng(seed + 2), **kw)
    out["amr"] = dict(zip(("face", "edge"), slice_arrays(ds, hw, zhalf, res)))
    out["amr_info"] = info
    if individual:
        save_yt_plots(ds, f"{outdir}/zoo_{name}_amr", title + " (pseudo-AMR)", grids=grids)
    del ds
    return out


def make_montages(frames, outdir=".", dyn_range=3.0):
    """Write zoo_montage_faceon.png and zoo_montage_edgeon.png, both with the
    same layout: rows = SPH / uniform / AMR, two blocks of five galaxies.
    Log color scale with a fixed dynamic range (dyn_range dex), shared by
    the three methods for each galaxy and view."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    names = [n for n in ORDER if n in frames]
    L = LENGTH_SCALE_KPC
    cmap = matplotlib.colormaps["viridis"].copy()
    cmap.set_bad(cmap(0.0))    # empty SPH pixels (zero density) -> floor colour
    cmap.set_under(cmap(0.0))

    def norm(f, key):
        # ONE normalisation per galaxy and view, shared by SPH / uniform / AMR,
        # so the three methods are directly comparable: vmax is the largest
        # 99.9th percentile among them, vmin is dyn_range dex below it.
        vmax = max(np.percentile(f[m][key], 99.9) for m in METHODS)
        return LogNorm(vmin=vmax / 10**dyn_range, vmax=vmax)

    blocks = [names[i:i + 5] for i in range(0, len(names), 5)]
    for key, label in (("face", "faceon"), ("edge", "edgeon")):
        fig, axes = plt.subplots(3 * len(blocks), 5, figsize=(16, 3.25 * 3 * len(blocks)),
                                 squeeze=False)
        for b, block in enumerate(blocks):
            for c in range(5):
                for m, method in enumerate(METHODS):
                    ax = axes[3 * b + m][c]
                    ax.set_xticks([])
                    ax.set_yticks([])
                    if c >= len(block):
                        ax.axis("off")
                        continue
                    f = frames[block[c]]
                    img = f[method][key]
                    e, ez = f["hw"] * L, f["zhalf"] * L
                    ax.imshow(img, origin="lower", extent=[-e, e, -ez, ez], cmap=cmap,
                              norm=norm(f, key), interpolation="nearest")
                    if m == 0:
                        ax.set_title(f"{TITLES[block[c]]}  ({2 * e:.0f} kpc)", fontsize=10)
                    if c == 0:
                        ax.set_ylabel(METHOD_LABELS[method], fontsize=11)
        what = "face-on (z = 0)" if key == "face" else "edge-on (y = 0)"
        fig.suptitle(f"Galaxy zoo, {what} midplane density slices", fontsize=14)
        fig.tight_layout(rect=(0, 0, 1, 0.985))
        fig.savefig(f"{outdir}/zoo_montage_{label}.png", dpi=110)
        plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--cache", default="zoo_cache.pkl",
                    help="pickle of the simulated particles; reused if it exists (integrators are skipped)")
    ap.add_argument("--frb-cache", default="zoo_frames.pkl",
                    help="pickle of rendered slice arrays; lets you re-lay-out montages instantly")
    ap.add_argument("--only", nargs="*", help="galaxy names to (re)generate/render")
    ap.add_argument("--n-disk", type=int, default=60000)
    ap.add_argument("--redo", action="store_true", help="ignore the rendered-slice cache")
    ap.add_argument("--individual", action="store_true", help="also save labelled per-galaxy yt PNGs")
    ap.add_argument("--grids", action="store_true", help="with --individual: AMR grid outlines")
    ap.add_argument("--res", type=int, default=512)
    args = ap.parse_args()

    try:
        with open(args.cache, "rb") as f:
            results = pickle.load(f)
        if args.only and not all(n in results for n in args.only):
            raise FileNotFoundError
        print("loaded cached particles from", args.cache)
    except (FileNotFoundError, EOFError):
        results = generate_all(n_disk=args.n_disk, only=args.only)
        with open(args.cache, "wb") as f:
            pickle.dump(results, f)

    try:
        with open(args.frb_cache, "rb") as f:
            frames = {} if args.redo else pickle.load(f)
    except (FileNotFoundError, EOFError):
        frames = {}

    for name in ORDER:
        if name not in results or (args.only and name not in args.only):
            continue
        if name in frames and not args.individual:
            continue
        frames[name] = render_galaxy(name, results[name], outdir=args.outdir,
                                     individual=args.individual, grids=args.grids, res=args.res)
        print(name, "AMR:", frames[name]["amr_info"], flush=True)
        with open(args.frb_cache, "wb") as f:  # incremental: a crash keeps progress
            pickle.dump(frames, f)

    make_montages(frames, outdir=args.outdir)
    print("wrote montages to", args.outdir)


# =====================================================================
# 2D vs 3D: what is and isn't simulated in 3D
# =====================================================================
# - Disk-type galaxies, mergers, dwarf: orbits are integrated in 2D (x, y)
#   only. The vertical coordinate z is drawn independently from a sech^2
#   profile (add_z_disk) and is NOT dynamically coupled to the in-plane
#   motion -- no vertical oscillations, warps, or disk heating.
# - Ellipticals: genuinely 3D positions (triaxial Plummer sampling), but
#   static -- no orbits are integrated.
# - Dwarf irregular: the fractal-noise weighting is a visual texture only
#   (applied identically in all three representations).
# - All three representations (SPH, uniform grid, pseudo-AMR) are built from
#   the same 3D particle set inside the same box, so edge-on views of disks
#   show the imposed sech^2 thickness, and ellipticals are round.
# =====================================================================
