"""
zoo.py

Toy galaxy zoo rendered with yt: spirals, barred spiral, lenticular,
ellipticals, interacting pair, merger remnant and dwarf irregular. Each
galaxy is shown as three representations built from the same 3D particle set
(SPH, uniform grid, 3D pseudo-AMR), face-on and edge-on, as midplane slices
("slic") and as projections ("proj").

Models
------
Disk types (spirals, bar, lenticular, dwarf):
    Test particles integrated in a 2D potential: flat rotation curve plus
    rotating Fourier perturbations (spiral / bar terms). Arms and bars emerge
    from orbit crowding, they are not painted on.
Ellipticals:
    Sampled directly from a triaxial Plummer sphere, positions and velocities
    (isotropic Plummer distribution function, no dynamics).
Interacting pair / merger remnant:
    Restricted Toomre-type encounter. The two centres follow a softened
    2-body orbit; each disk feels its own potential plus a Plummer tidal term
    from the other galaxy.
Dwarf irregular:
    Dynamical disk whose particle masses carry a fractal noise texture
    (visual only, total mass unchanged).

Physical units
--------------
Every galaxy has its own unit system (PHYS): a length unit L (kpc) and a
velocity unit V (km/s). With G = 1 in code units this fixes the time unit
L / V and the mass unit V^2 L / G. The potentials are in these units, so the
masses follow from the dynamics:
- disks: the disk supplies a fraction f_disk of V0^2 at 2.2 r_d (the peak of
  an exponential disk's rotation curve); the rest is dark matter. Disk mass
  M_d = f_disk V0^2 2.2 r_d / (G 0.645), shared equally by the particles.
- ellipticals: self-gravitating Plummer spheres, V = sqrt(G M_star / L), so
  M_star is exact and consistent with the sampled velocities.
- mergers: each disk's circular speed follows from its own mass in the
  2-body orbit, V0_i = sqrt(G M_i / (2 r0_max_i)); disk masses as above.
Disk scale heights are physical (h_kpc), not tied to the box size.

Kinematics and thermodynamics
-----------------------------
All particles carry 3D velocities. Disk orbits are 2D (x, y); z is drawn
from a sech^2(z / z0) profile. The strong spiral / bar forcing that makes the
arms by orbit crowding also leaves the integrated disk velocities
unphysically hot, so disk velocities are a kinematic model on top of the
simulated positions (disk_kinematics): circular speed of the potential, arm
streaming from the local density contrast, and isotropic random motion
sigma^2 = pi G Sigma z0 of the isothermal sheet matching the sech^2 profile.
Merger tidal debris blends into its simulated (ballistic) velocities.
Ellipticals use their sampled Plummer velocities.

Temperature and pressure are kinetic: each particle's 1D velocity dispersion
sigma^2 is measured over its 32 nearest neighbours after removing a local
linear flow fit (bulk motion, rotation and shear), then
    T = mu m_p sigma^2 / k_B   (mu = 0.6),   P = rho sigma^2,
so P = rho k_B T / (mu m_p) holds exactly. A floor of T_FLOOR applies.

Native yt fields
----------------
All three representations carry ("gas", ...) density, velocity_x/y/z
(hence velocity_magnitude, velocity_cylindrical_theta, ...), temperature and
pressure. The grids deposit the extensive quantities (mass, momentum,
mass * sigma^2) with one shared adaptive kernel and form v = p / m,
sigma^2 = E / m per cell, so mass and momentum are conserved across levels.

Pseudo-AMR
----------
Octree-style patch hierarchy with cubic cells. Blocks outside a central
ellipsoid stay at level 0; everything inside is refined at least once. Deeper
levels follow an oblate ellipsoid around each nucleus (scaled to the disk
scale height) or a density threshold, so arms, bars and tidal features keep
their resolution. One extra central level resolves the nuclei by mass-
conserving prolongation from the actual particles.

Caching
-------
The particle simulation is cached in a pickle (--cache) and the rendered image
arrays in a second one (--frb-cache), so re-rendering and re-laying-out the
montages never re-runs the integrators. Caches written by an older version
(CACHE_VERSION) are regenerated.

Notebook use
------------
    results = load_or_generate("zoo_cache.pkl")
    g = prepare_galaxy("1_grand_design", results["1_grand_design"])
    ds, info = build_ds(g, "amr")            # "sph" | "uniform" | "amr"
    yt_plot(ds, "z", "slic", field="temperature")   # renders inline
    imgs = field_images(ds, g, "slic")       # density, pressure, |v|, T arrays
    prof = radial_profiles(ds, g)            # scaled radial profiles

Requires: numpy, scipy, matplotlib, yt.
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

CACHE_VERSION = 4

# =====================================================================
# Physical constants and unit systems
# =====================================================================

G_KPC = 4.30091e-6           # G in kpc (km/s)^2 / Msun
MYR_PER_KPC_KMS = 977.792    # 1 kpc / (km/s) in Myr
MU = 0.6                     # mean molecular weight (ionised gas)
K_PER_KMS2 = MU * 1.67262192e-27 / 1.380649e-23 * 1e6  # T [K] per sigma^2 [(km/s)^2]
DYN_PER_MSUN_KPC3_KMS2 = 1.98847e33 / 3.0856776e21**3 * 1e10  # rho sigma^2 -> dyn/cm**2
T_FLOOR = 100.0              # K
S2_FLOOR = T_FLOOR / K_PER_KMS2

# Per-galaxy physical calibration. L: code length (kpc); V: code velocity
# (km/s), or M (Msun) for ellipticals, where V = sqrt(G M / L); f_disk: disk
# share of V0^2 at 2.2 r_d; h_kpc: sech^2 scale height z0 of the disk.
PHYS = {
    "1_grand_design": dict(L=10.0, V=220.0, f_disk=0.5, h_kpc=0.35),
    "2_tightly_wound": dict(L=10.0, V=250.0, f_disk=0.55, h_kpc=0.3),
    "3_flocculent": dict(L=8.0, V=170.0, f_disk=0.4, h_kpc=0.35),
    "4_barred_spiral": dict(L=10.0, V=220.0, f_disk=0.55, h_kpc=0.35),
    "5_lenticular": dict(L=8.0, V=230.0, f_disk=0.6, h_kpc=0.5),
    "6_elliptical_E0": dict(L=10.0, M=1.5e11),
    "7_elliptical_E5": dict(L=10.0, M=8.0e10),
    "8_interacting_pair": dict(L=10.0, V=220.0, f_disk=0.5, h_kpc=0.35),
    "9_merger_remnant": dict(L=10.0, V=220.0, f_disk=0.5, h_kpc=0.35),
    "10_dwarf_irregular": dict(L=5.0, V=120.0, f_disk=0.2, h_kpc=0.4),
}

RC_CORE = 0.08            # core radius of the logarithmic potential (code units)
EXP_DISK_FRAC_22 = 0.6454  # mass fraction of an exponential disk inside 2.2 r_d


def make_units(L, V=None, M=None, **_):
    """Unit system from a PHYS entry: L (kpc) and V (km/s), or L and M (Msun)
    for a self-gravitating system of total mass M."""
    if V is None:
        V = np.sqrt(G_KPC * M / L)
    return dict(L_kpc=L, V_kms=V, t_Myr=L / V * MYR_PER_KPC_KMS, M_Msun=V**2 * L / G_KPC)


def yt_units(u):
    """yt load kwargs: the galaxy's unit system, non-periodic (isolated
    galaxy; periodic wrapping also breaks yt's radius fields on stream grids)."""
    return dict(length_unit=(u["L_kpc"], "kpc"), mass_unit=(u["M_Msun"], "Msun"),
                time_unit=(u["t_Myr"], "Myr"), periodicity=(False, False, False))


def disk_mass_code(V0, rd, f_disk):
    """Exponential-disk mass (code units) supplying f_disk of V0^2 at 2.2 rd."""
    return f_disk * V0**2 * 2.2 * rd / EXP_DISK_FRAC_22


def mdyn_log_code(r, V0, RC=RC_CORE):
    """Dynamical mass inside r of the logarithmic potential (code units)."""
    return V0**2 * r**3 / (r**2 + RC**2)


# =====================================================================
# Fractal noise
# =====================================================================

def fractal_noise(shape, rng, octaves=4, base_sigma=6.0, persistence=0.6):
    """Multi-octave smoothed noise: white noise filtered at a ladder of
    scales and summed. Any ndim."""
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
    """Perturbation specs -> internal dicts. Spec keys: m, pitch_deg, amp,
    r_peak, r_width, omega_p."""
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
    """Acceleration for a disk centred on the origin (callers subtract the
    galaxy's own centre first). Optional Plummer tidal perturber at c_other
    with mass parameter GM_other."""
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
    """Radii from an exponential surface-density disk (rejection sampling)."""
    r0 = np.empty(0)
    while r0.shape[0] < n:
        cand = rng.uniform(0.02, r0_max, size=n)
        pdf = cand * np.exp(-cand / rd)
        accept = rng.uniform(0, pdf.max(), size=n) < pdf
        r0 = np.concatenate([r0, cand[accept]])
    return r0[:n]


def simulate_disk_galaxy(
    specs, n=60000, t_total=9.0, dt=0.004, fade_time=2.5,
    V0=1.0, RC=RC_CORE, rd=0.35, r0_max=3.0, r_kill=4.0,
    vdisp=0.02, vphi_scatter=0.04,
):
    """Leapfrog-integrate n test particles in the axisymmetric + `specs`
    potential. Perturbation amplitudes ramp up linearly over fade_time.
    Returns (x, y, vx, vy) of the particles that end inside r_kill."""
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
    return x[mask], y[mask], vx[mask], vy[mask]


# =====================================================================
# Ellipticals: direct Plummer-sphere sampling
# =====================================================================

def simulate_elliptical(n=40000, a=0.35, qx=1.0, qy=1.0, qz=0.7, r_cap_factor=6.0, rng=rng):
    """Sample n particles from a Plummer sphere of unit mass (G = 1, scale a),
    stretched by the axis ratios (qx, qy, qz). Radii are capped at
    r_cap_factor * a. Speeds come from the isotropic Plummer distribution
    function (Aarseth, Henon & Wielen 1974): v = q v_esc(r) with q drawn from
    q^2 (1 - q^2)^(7/2). Velocities are stretched by the same axis ratios,
    which keeps the virial balance only approximately for q != 1.
    Returns (x, y, z, vx, vy, vz)."""
    u = rng.uniform(0, 0.995, size=n)
    r = a / np.sqrt(u ** (-2.0 / 3.0) - 1.0)
    r = np.minimum(r, a * r_cap_factor)

    def iso(size):
        costh = rng.uniform(-1, 1, size=size)
        sinth = np.sqrt(1 - costh**2)
        phi = rng.uniform(0, 2 * np.pi, size=size)
        return sinth * np.cos(phi), sinth * np.sin(phi), costh

    ux, uy, uz = iso(n)
    q = np.empty(0)
    while q.shape[0] < n:
        cand = rng.uniform(0, 1, size=n)
        accept = rng.uniform(0, 0.1, size=n) < cand**2 * (1 - cand**2) ** 3.5
        q = np.concatenate([q, cand[accept]])
    v = q[:n] * np.sqrt(2.0) * (r**2 + a**2) ** -0.25
    wx, wy, wz = iso(n)
    return (r * ux * qx, r * uy * qy, r * uz * qz,
            v * wx * qx, v * wy * qy, v * wz * qz)


# =====================================================================
# Mergers: two disks on a mutual-gravity encounter
# =====================================================================

def two_body_trajectory(GM_tot, r12_init, v12_init, t_total, dt, eps=0.3):
    """Leapfrog-integrate the relative separation r12 = c2 - c1 of two
    Plummer-softened point masses. Returns an (nsteps + 1, 2) array."""
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


def merger_circular_speeds(GM_tot, m1_frac, r0_max1, r0_max2, halo_k=2.0):
    """Circular speed of each disk from its own mass in the 2-body orbit:
    V0_i^2 = G M_i / (halo_k r0_max_i), i.e. the halo mass M_i is enclosed
    at halo_k times the disk's outer radius."""
    GM1 = GM_tot * m1_frac
    GM2 = GM_tot * (1 - m1_frac)
    return np.sqrt(GM1 / (halo_k * r0_max1)), np.sqrt(GM2 / (halo_k * r0_max2))


def simulate_merger(
    specs1, specs2, n1=20000, n2=20000, t_total=9.0, dt=0.004, fade_time=1.0,
    RC=RC_CORE, rd1=0.3, rd2=0.3, r0_max1=1.6, r0_max2=1.6, r_kill=8.0,
    GM_tot=3.0, r12_init=(-4.5, 1.0), v12_init=(0.95, -0.05), eps_gal=0.35,
    m1_frac=0.5, spin1=1.0, spin2=1.0, halo_k=2.0,
):
    """Two disks, each with its own potential, size, particle count and spin
    sense (spin = +-1), centred on moving centres. The centres follow a
    mutual 2-body orbit; each disk also feels a Plummer tidal term from the
    other. Prograde vs retrograde spin gives the asymmetric tidal response of
    real pairs. Each disk's circular speed follows from its mass
    (merger_circular_speeds). Returns (x, y, vx, vy, owner, info) of the
    surviving particles from both disks (owner = 0 / 1); info holds the final
    centres, centre velocities, V0, spins and r0_max of both disks."""
    comps1 = make_components(specs1)
    comps2 = make_components(specs2)
    GM1 = GM_tot * m1_frac
    GM2 = GM_tot * (1 - m1_frac)
    V01, V02 = merger_circular_speeds(GM_tot, m1_frac, r0_max1, r0_max2, halo_k)

    traj = two_body_trajectory(GM_tot, r12_init, v12_init, t_total, dt, eps=eps_gal)
    c1_traj = -(1 - m1_frac) * traj  # centre of mass fixed at the origin
    c2_traj = m1_frac * traj

    v12_now = np.array(v12_init, dtype=float)
    v_c1_init = -(1 - m1_frac) * v12_now
    v_c2_init = m1_frac * v12_now

    def make_disk(n, rd_, r0max_, c_init, v_c_init, spin, V0):
        r0 = sample_exponential_disk(n, rd_, r0max_, rng)
        theta0 = rng.uniform(0, 2 * np.pi, size=n)
        x = r0 * np.cos(theta0) + c_init[0]
        y = r0 * np.sin(theta0) + c_init[1]
        vphi = spin * V0 * r0 / np.sqrt(r0**2 + RC**2)
        vx = -vphi * np.sin(theta0) + v_c_init[0] + rng.normal(0, 0.02, size=n)
        vy = vphi * np.cos(theta0) + v_c_init[1] + rng.normal(0, 0.02, size=n)
        return x, y, vx, vy

    x1, y1, vx1, vy1 = make_disk(n1, rd1, r0_max1, c1_traj[0], v_c1_init, spin1, V01)
    x2, y2, vx2, vy2 = make_disk(n2, rd2, r0_max2, c2_traj[0], v_c2_init, spin2, V02)

    x = np.concatenate([x1, x2])
    y = np.concatenate([y1, y2])
    vx = np.concatenate([vx1, vx2])
    vy = np.concatenate([vy1, vy2])
    owner = np.concatenate([np.zeros(n1, dtype=int), np.ones(n2, dtype=int)])

    nsteps = int(t_total / dt)
    sel1 = owner == 0
    sel2 = owner == 1

    def accel_all(x, y, t, step_idx):
        c1 = c1_traj[step_idx]
        c2 = c2_traj[step_idx]
        scale = min(1.0, t / fade_time)
        comps1_f = [dict(c, amp=c["amp"] * scale) for c in comps1]
        comps2_f = [dict(c, amp=c["amp"] * scale) for c in comps2]

        ax = np.empty_like(x)
        ay = np.empty_like(y)
        ax1, ay1 = accel_disk(
            x[sel1] - c1[0], y[sel1] - c1[1], t, comps1_f, V01, RC,
            c_other=(c2[0] - c1[0], c2[1] - c1[1]), GM_other=GM2, eps_other=eps_gal,
        )
        ax2, ay2 = accel_disk(
            x[sel2] - c2[0], y[sel2] - c2[1], t, comps2_f, V02, RC,
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
    v12_end = (traj[-1] - traj[-2]) / dt
    info = dict(centers=[tuple(c1_traj[-1]), tuple(c2_traj[-1])],
                vcenters=[tuple(-(1 - m1_frac) * v12_end), tuple(m1_frac * v12_end)],
                V0=(V01, V02), spins=(spin1, spin2), r0_max=(r0_max1, r0_max2))
    return x[mask], y[mask], vx[mask], vy[mask], owner[mask], info


# =====================================================================
# Vertical structure and particle kinematics
# =====================================================================

H_VERT = 0.035  # default disk scale height (code units); galaxies use PHYS h_kpc


def add_z_disk(x, y, h=H_VERT, zlim=None, rng=rng):
    """Draw z from a sech^2(z/h) profile (inverse-CDF sampling), independent
    of (x, y)."""
    n = x.shape[0]
    u = rng.uniform(1e-4, 1 - 1e-4, size=n)
    z = h * np.arctanh(2 * u - 1)
    if zlim is not None:
        z = np.clip(z, zlim[0], zlim[1])
    return z


STREAM_AMP = 0.1  # arm streaming: fractional v_phi deficit at saturated density contrast


def disk_kinematics(x, y, mass_code, z0, hosts, rng, k=64, nbins=40):
    """Model 3D velocities (code units) for disk particles.

    The test-particle integrations make the arms and bars by orbit crowding
    under deliberately strong forcing, which leaves their own velocities
    unphysically hot (mean v_phi ~ 0.3 V0, sigma ~ 0.6 V0). The positions
    are kept; the velocities are rebuilt from the potential:
    - ordered flow: circular speed of the host's logarithmic potential,
      v_c = V0 R / sqrt(R^2 + RC^2), about the host centre (moving with it);
    - arm streaming: with delta the local surface density over its
      azimuthal mean at that radius, v_phi -> v_c (1 - A tanh(delta)) and
      v_R = -A/2 |v_c| tanh(delta) (slower, inflowing gas in the arms),
      A = STREAM_AMP;
    - random motion: isotropic, sigma^2 = pi G Sigma z0 with Sigma the local
      surface density (k nearest neighbours), i.e. the isothermal sheet
      whose sech^2(z / z0) profile the z coordinates are drawn from.

    hosts: list of dicts with center, vcenter (code units), V0, spin (+-1),
    sel (bool mask of the host's particles), and optionally r_in, r_out and
    v_sim (n_sel, 2): beyond r_in the ordered in-plane velocity blends
    linearly into v_sim (tidal debris keeps its simulated bulk flow), fully
    at r_out; the random part is kept throughout."""
    n = x.shape[0]
    xy = np.column_stack([x, y])
    d, idx = cKDTree(xy).query(xy, k=min(k, n - 1))
    sig_loc = mass_code[idx].sum(axis=1) / (np.pi * d[:, -1] ** 2)
    vel = rng.normal(size=(n, 3)) * np.sqrt(np.pi * sig_loc * z0)[:, None]
    for h in hosts:
        sel = h["sel"]
        dx, dy = x[sel] - h["center"][0], y[sel] - h["center"][1]
        R = np.maximum(np.hypot(dx, dy), 1e-6)
        cos, sin = dx / R, dy / R
        vc = h["spin"] * h["V0"] * R / np.sqrt(R**2 + RC_CORE**2)
        edges = np.unique(np.quantile(R, np.linspace(0, 1, nbins + 1)))
        ann = np.histogram(R, edges, weights=mass_code[sel])[0] / (np.pi * np.diff(edges**2))
        sig_R = ann[np.clip(np.searchsorted(edges, R, side="right") - 1, 0, len(ann) - 1)]
        t = np.tanh(sig_loc[sel] / sig_R - 1.0)
        vphi = vc * (1.0 - STREAM_AMP * t)
        vR = -0.5 * STREAM_AMP * np.abs(vc) * t
        vxm = vR * cos - vphi * sin + h["vcenter"][0]
        vym = vR * sin + vphi * cos + h["vcenter"][1]
        if "r_in" in h:
            w = np.clip((h["r_out"] - R) / (h["r_out"] - h["r_in"]), 0.0, 1.0)
            vxm = w * vxm + (1 - w) * h["v_sim"][:, 0]
            vym = w * vym + (1 - w) * h["v_sim"][:, 1]
        vel[sel, 0] += vxm
        vel[sel, 1] += vym
    return vel


def particle_sph(g, k=32):
    """kNN smoothing length (code units), SPH density (Msun/kpc**3) and the
    neighbour indices, computed once per prepared galaxy."""
    if g.get("_sph") is None:
        pos = g["pos"]
        d, idx = cKDTree(pos).query(pos, k=min(k, pos.shape[0] - 1))
        hsml = np.maximum(d[:, -1], 1e-4)
        rho = g["mass"][idx].sum(axis=1) / ((4.0 / 3.0) * np.pi * (hsml * g["L"]) ** 3)
        g["_sph"] = (hsml, rho, idx)
    return g["_sph"]


def particle_dispersion(g, chunk=20000):
    """1D velocity dispersion sigma^2 ((km/s)^2) of each particle's kNN
    neighbourhood, after removing a local linear flow v = v0 + B (x - x0)
    (least squares per component, lightly ridge-regularised). The fit takes
    out bulk motion, rotation and shear across the kernel, so only the
    random motion counts as 'thermal'. Floored at T_FLOOR."""
    if g.get("_s2") is None:
        _, _, idx = particle_sph(g)
        pos, vel = g["pos"], g["vel"]
        k = idx.shape[1]
        s2 = np.empty(pos.shape[0])
        eye = np.eye(3)
        for i0 in range(0, pos.shape[0], chunk):
            nb = idx[i0:i0 + chunk]
            dx = pos[nb]
            dv = vel[nb]
            dx = dx - dx.mean(axis=1, keepdims=True)
            dv = dv - dv.mean(axis=1, keepdims=True)
            C = np.einsum("nki,nkj->nij", dx, dx)
            ridge = 1e-3 * np.trace(C, axis1=1, axis2=2)[:, None, None] / 3 + 1e-12
            B = np.linalg.solve(C + ridge * eye, np.einsum("nki,nkj->nij", dx, dv))
            res = dv - np.einsum("nki,nij->nkj", dx, B)
            s2[i0:i0 + chunk] = (res**2).sum(axis=(1, 2)) / (3 * (k - 4))
        g["_s2"] = np.maximum(s2, S2_FLOOR)
    return g["_s2"]


def z_half_extent(hw, spheroid):
    """Vertical half-extent of the box. Always a cube (= in-plane half-width)
    for every galaxy type, so all views share one layout."""
    return hw


DENS = ("gas", "density")


# =====================================================================
# SPH representation
# =====================================================================

def build_sph_ds(g, fields="all"):
    """yt SPH particle dataset of a prepared galaxy.

    Smoothing lengths, densities and dispersions use the full particle set,
    including particles outside the box, so edge particles keep their
    neighbours. Only particles inside the box are passed on, since yt
    requires all particles to lie within the bounding box. fields="all" adds
    velocities, temperature and pressure (exposed as ("gas", ...))."""
    xlim, ylim, zlim = g["xlim"], g["ylim"], g["zlim"]
    pos, vel, mass = g["pos"], g["vel"], g["mass"]
    hsml, rho, _ = particle_sph(g)
    inside = (
        (pos[:, 0] > xlim[0]) & (pos[:, 0] < xlim[1]) &
        (pos[:, 1] > ylim[0]) & (pos[:, 1] < ylim[1]) &
        (pos[:, 2] > zlim[0]) & (pos[:, 2] < zlim[1])
    )
    data = {
        "particle_position_x": pos[inside, 0],
        "particle_position_y": pos[inside, 1],
        "particle_position_z": pos[inside, 2],
        "particle_mass": (mass[inside], "Msun"),
        "density": (rho[inside], "Msun/kpc**3"),
        "smoothing_length": hsml[inside],
    }
    if fields == "all":
        s2 = particle_dispersion(g)[inside]
        for i, ax in enumerate("xyz"):
            data[f"particle_velocity_{ax}"] = (vel[inside, i], "km/s")
        data["temperature"] = (s2 * K_PER_KMS2, "K")
        # the "particle_" prefix makes yt alias it to ("gas", "pressure")
        data["particle_pressure"] = (rho[inside] * s2 * DYN_PER_MSUN_KPC3_KMS2, "dyn/cm**2")
    bbox = np.array([xlim, ylim, zlim])
    return yt.load_particles(data, bbox=bbox, **yt_units(g["units"]))


# =====================================================================
# Grid-based representations: shared extensive fields
# =====================================================================
# Both grids start from smoothed "master" fields of extensive quantities
# (mass, momentum, mass * sigma^2). Coarser cells are exact block sums of
# finer ones, so mass and momentum are conserved across levels; intensive
# fields are formed only at the end (fields_to_yt).

def extensive_weights(g, inside, fields="all"):
    """Per-particle extensive quantities for the grid deposit: mass (Msun),
    momentum px, py, pz (Msun km/s) and e = mass * sigma^2."""
    m = g["mass"][inside]
    ext = {"mass": m}
    if fields == "all":
        v = g["vel"][inside]
        ext.update(px=m * v[:, 0], py=m * v[:, 1], pz=m * v[:, 2],
                   e=m * particle_dispersion(g)[inside])
    return ext


def deposit_fields(pos, hsml, ext, xlim, ylim, zlim, shape, sigma_frac=0.5, min_sigma=0.7,
                   bins_per_octave=2, floor=3e-4):
    """Deposit each particle as an isotropic Gaussian with
    sigma = sigma_frac * hsml (the same kNN length the SPH path uses), so the
    grid inherits SPH's adaptive resolution. Particles are grouped into
    log-spaced sigma bins; each bin is histogrammed and filtered once per
    field, with the same kernel for every field so ratios stay exact.

    A small floor mass (floor * peak cell, capped at 1% of the total mass) at
    zero velocity and T_FLOOR fills empty cells; the total mass is
    renormalised to the deposited mass and the other fields are scaled by
    the same factor."""
    lo = np.array([xlim[0], ylim[0], zlim[0]])
    hi = np.array([xlim[1], ylim[1], zlim[1]])
    cell = (hi - lo) / np.array(shape)
    ijk = np.floor((pos - lo) / cell).astype(np.int64)
    ok = np.all((ijk >= 0) & (ijk < np.array(shape)), axis=1)
    flat = np.ravel_multi_index(ijk[ok].T, shape)
    sig = np.maximum(sigma_frac * hsml[ok] / cell[0], min_sigma)  # in cells
    key = np.round(np.log2(sig) * bins_per_octave).astype(int)
    w_ok = {name: w[ok] for name, w in ext.items()}
    ncell = int(np.prod(shape))

    out = {name: np.zeros(shape) for name in ext}
    for kk in np.unique(key):
        m = key == kk
        s = 2.0 ** (kk / bins_per_octave)
        for name, w in w_ok.items():
            H = np.bincount(flat[m], weights=w[m], minlength=ncell).reshape(shape)
            out[name] += gaussian_filter(H, sigma=s, truncate=3.0)

    mass = out["mass"]
    tot = mass.sum()
    f = min(mass.max() * floor, 0.01 * tot / ncell)
    mass += f
    c = tot / mass.sum()
    mass *= c
    for name in ("px", "py", "pz"):
        if name in out:
            out[name] = (out[name] * c).astype(np.float32)
    if "e" in out:
        out["e"] = ((out["e"] + f * S2_FLOOR) * c).astype(np.float32)
    return out


def block_reduce_sum(field, factor):
    """Downsample an N-d array by summing factor^ndim blocks. All dimensions
    must be divisible by factor."""
    shp = []
    for n in field.shape:
        shp += [n // factor, factor]
    return field.reshape(shp).sum(axis=tuple(range(1, 2 * field.ndim, 2)))


def fields_to_yt(ext, left_edge, right_edge, L):
    """Cell extensive quantities -> yt field dict: density (Msun/kpc**3) and,
    if present, velocity_x/y/z (km/s), temperature (K), pressure (dyn/cm**2)."""
    m = np.asarray(ext["mass"], dtype=float)
    vol = 1.0
    for n, lo, hi in zip(m.shape, left_edge, right_edge):
        vol *= (hi - lo) / n * L
    rho = m / vol
    out = {"density": (rho, "Msun/kpc**3")}
    if "px" in ext:
        for ax in "xyz":
            out[f"velocity_{ax}"] = (ext[f"p{ax}"] / m, "km/s")
        s2 = ext["e"] / m
        out["temperature"] = (s2 * K_PER_KMS2, "K")
        out["pressure"] = (rho * s2 * DYN_PER_MSUN_KPC3_KMS2, "dyn/cm**2")
    return out


# =====================================================================
# Uniform-grid representation
# =====================================================================

def build_uniform_ds(pos, hsml, ext, xlim, ylim, zlim, units, n_cells=56, supersample=4):
    """Fixed-resolution uniform grid with cubic cells. n_cells is the count
    across x and y; the z count follows from the cubic cells. The fields are
    deposited at `supersample` times the resolution and block-summed down."""
    dx = (xlim[1] - xlim[0]) / n_cells
    nz = int(round((zlim[1] - zlim[0]) / dx))
    shape = (n_cells * supersample, n_cells * supersample, nz * supersample)
    fine = deposit_fields(pos, hsml, ext, xlim, ylim, zlim, shape)
    coarse = {k: block_reduce_sum(v, supersample) for k, v in fine.items()}
    del fine
    data = fields_to_yt(coarse, [xlim[0], ylim[0], zlim[0]], [xlim[1], ylim[1], zlim[1]],
                        units["L_kpc"])
    bbox = np.array([xlim, ylim, list(zlim)])
    return yt.load_uniform_grid(data, coarse["mass"].shape, bbox=bbox, **yt_units(units))


# =====================================================================
# yt output helpers
# =====================================================================
# kind: "slic" = midplane slice, "proj" = line-of-sight projection. Density
# is projected without weight (column density); intensive fields are
# density-weighted (mass-weighted line-of-sight averages).
# Face-on = normal to z; edge-on = normal to y.

KINDS = ("slic", "proj")
# slices sit a hair off the midplane: z = 0 and y = 0 are cell faces of every
# grid, where the cell lookup is ambiguous and AMR slices show seams
SLICE_EPS = 1e-6
UNITS = {"slic": "Msun/kpc**3", "proj": "Msun/kpc**2"}
KIND_LABELS = {"slic": "midplane density slices", "proj": "projected density"}

FIELDS = ("density", "pressure", "velocity_magnitude", "temperature")
FIELD_STYLE = {
    "density": dict(label="density", unit={"slic": "Msun/pc**3", "proj": "Msun/pc**2"},
                    tex={"slic": r"$\rho$ [M$_\odot$ pc$^{-3}$]",
                         "proj": r"$\Sigma$ [M$_\odot$ pc$^{-2}$]"},
                    cmap="viridis", log=True, dex=4.0, weight=None),
    "pressure": dict(label="pressure", unit={"slic": "dyn/cm**2", "proj": "dyn/cm**2"},
                     tex={"slic": r"P [dyn cm$^{-2}$]", "proj": r"$\langle P\rangle_\rho$ [dyn cm$^{-2}$]"},
                     cmap="magma", log=True, dex=5.0, weight=DENS),
    "velocity_magnitude": dict(label="|v|", unit={"slic": "km/s", "proj": "km/s"},
                               tex={"slic": "|v| [km/s]", "proj": r"$\langle|v|\rangle_\rho$ [km/s]"},
                               cmap="plasma", log=False, dex=None, weight=DENS),
    "temperature": dict(label="temperature", unit={"slic": "K", "proj": "K"},
                        tex={"slic": "T [K]", "proj": r"$\langle T\rangle_\rho$ [K]"},
                        cmap="inferno", log=True, dex=3.0, weight=DENS),
}


def _style(p, title, kind="slic", grids=False):
    p.set_unit(DENS, UNITS[kind])
    p.set_axes_unit("kpc")
    p.set_cmap(DENS, "viridis")
    p.annotate_title(title)
    if grids:
        p.annotate_grids(periodic=False)


def yt_plot(ds, axis="z", kind="slic", title=None, grids=False, field="density",
            velocity=False):
    """Styled native yt plot (SlicePlot or ProjectionPlot) of one of FIELDS.
    In a notebook, leave it as the last expression of a cell (or call
    .show()) to render it inline; .save(), .zoom() and the other yt plot
    methods work as usual. Edge-on (axis "y") plots are swapped so the disk
    lies flat. velocity=True overlays in-plane velocity arrows."""
    st = FIELD_STYLE[field]
    f = ("gas", field)
    if kind == "slic":
        p = yt.SlicePlot(ds, axis, f, center=ds.arr([SLICE_EPS] * 3, "code_length"))
    else:
        p = yt.ProjectionPlot(ds, axis, f, weight_field=st["weight"])
    p.set_unit(f, st["unit"][kind])
    p.set_log(f, st["log"])
    p.set_cmap(f, st["cmap"])
    p.set_axes_unit("kpc")
    if title:
        p.annotate_title(title)
    if grids:
        p.annotate_grids(periodic=False)
    if velocity:
        p.annotate_velocity(factor=24, normalize=True)
    if axis == "y":
        p.swap_axes()
    return p


def save_yt_plots(ds, base, title, grids=False, kind="slic"):
    """Write {base}_faceon_{kind}.png and {base}_edgeon_{kind}.png (density,
    in the montage units)."""
    for axis, tag in (("z", "faceon"), ("y", "edgeon")):
        if kind == "slic":
            p = yt.SlicePlot(ds, axis, DENS, center=ds.arr([SLICE_EPS] * 3, "code_length"))
        else:
            p = yt.ProjectionPlot(ds, axis, DENS)
        _style(p, title, kind=kind, grids=grids)
        p.save(f"{base}_{tag}_{kind}.png")


def image_arrays(ds, hw, zhalf, res=512, kind="slic", field=DENS, unit=None, weight=None):
    """Face-on and edge-on images as plain arrays for imshow(origin="lower").
    Face-on: res x res, x horizontal / y vertical. Edge-on: res wide and
    res * zhalf / hw tall (cubic pixels), x horizontal / z vertical.
    Defaults give density in the montage units; projections use `weight`."""
    q = lambda v: ds.quan(v, "code_length")
    nz_res = max(8, int(round(res * zhalf / hw)))
    unit = unit or UNITS[kind]
    src = ((lambda ax: ds.slice(ax, SLICE_EPS)) if kind == "slic"
           else (lambda ax: ds.proj(field, ax, weight_field=weight)))
    face = src("z").to_frb(q(2 * hw), res)[field].to(unit).d
    # a y-normal image has axes (horizontal = z, vertical = x) and frb arrays
    # are [vertical, horizontal]; transpose so the galaxy lies flat
    edge = src("y").to_frb(q(2 * zhalf), (nz_res, res), height=q(2 * hw))[field].to(unit).d
    return np.ascontiguousarray(face), np.ascontiguousarray(edge.T)


def field_images(ds, g, kind="slic", fields=FIELDS, res=384):
    """{field: (face, edge)} image arrays of FIELDS for one dataset, in the
    FIELD_STYLE units (slices, or density-weighted projections)."""
    out = {}
    for f in fields:
        st = FIELD_STYLE[f]
        out[f] = image_arrays(ds, g["hw"], g["zhalf"], res, kind, field=("gas", f),
                              unit=st["unit"][kind],
                              weight=st["weight"] if kind == "proj" else None)
    return out


def plot_field_grid(images, g, view="face", kind="slic", fields=FIELDS, methods=None,
                    figsize_per=3.4):
    """Rows = representations, columns = fields. `images` is
    {method: field_images(...)}; view is "face" or "edge". Each column shares
    one colour scale across the representations. Returns the figure."""
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, Normalize

    methods = [m for m in (methods or METHODS) if m in images]
    vi = 0 if view == "face" else 1
    e, ez = g["hw"] * g["L"], g["zhalf"] * g["L"]
    fig, axes = plt.subplots(len(methods), len(fields), squeeze=False, layout="constrained",
                             figsize=(figsize_per * len(fields), figsize_per * len(methods) + 0.9))
    for c, f in enumerate(fields):
        st = FIELD_STYLE[f]
        arrs = [images[m][f][vi] for m in methods]
        if st["log"]:
            vmax = max(np.percentile(a[a > 0], 99.9) for a in arrs if np.any(a > 0))
            norm = LogNorm(vmin=vmax / 10 ** st["dex"], vmax=vmax)
        else:
            norm = Normalize(vmin=0.0, vmax=max(np.percentile(a, 99.5) for a in arrs))
        cmap = matplotlib.colormaps[st["cmap"]].copy()
        cmap.set_bad(cmap(0.0))
        cmap.set_under(cmap(0.0))
        for r, (m, a) in enumerate(zip(methods, arrs)):
            ax = axes[r][c]
            im = ax.imshow(a, origin="lower", extent=[-e, e, -ez, ez], cmap=cmap, norm=norm,
                           interpolation="nearest")
            if c == 0:
                ax.set_ylabel(f"{METHOD_LABELS[m]}\n{'y' if view == 'face' else 'z'} [kpc]")
            else:
                ax.set_yticklabels([])
            if r == len(methods) - 1:
                ax.set_xlabel("x [kpc]")
            else:
                ax.set_xticklabels([])
        fig.colorbar(im, ax=axes[:, c], location="top", shrink=0.9, aspect=25,
                     label=st["tex"][kind])
    what = {"slic": "midplane slices" if view == "face" else "y = 0 slices",
            "proj": "projections (density-weighted)"}[kind]
    fig.suptitle(f"{g['title']}, {'face-on' if view == 'face' else 'edge-on'} {what}",
                 fontsize=13)
    return fig


# =====================================================================
# Radial profiles
# =====================================================================

def profile_centers(g):
    """Profile centres (code units, z = 0): the nuclei found by find_centers
    (both nuclei for the mergers), strongest first."""
    x, y = g["pos"][:, 0], g["pos"][:, 1]
    return find_centers(x, y, g["xlim"], g["ylim"], n_centers=g["d"].get("n_centers", 1),
                        weights=g["mass"])


def profile_targets(g):
    """[(center, r_scale_kpc, label)] to profile: the nucleus with the
    galaxy's scale length, or for the mergers each nucleus with the scale
    length of the disk whose final centre lies nearest to it."""
    centers = profile_centers(g)
    d = g["d"]
    if "pair" not in d:
        return [(centers[0], g["r_scale_kpc"], "")]
    out = []
    for c in centers:
        i = int(np.argmin([np.hypot(c[0] - p[0], c[1] - p[1]) for p in d["pair"]["centers"]]))
        out.append((c, d["r_scale_pair"][i] * g["L"],
                    f"around the {('primary', 'secondary')[i]} nucleus "
                    f"({c[0] * g['L']:.0f}, {c[1] * g['L']:.0f}) kpc"))
    return out


def radial_profiles(ds, g, center=None, nbins=32, rmin_frac=0.05, rmax=None, r_scale_kpc=None):
    """Mass-weighted radial profiles around `center` (code units; default:
    the main nucleus). Disks use cylindrical radius in a slab of +-3 z0
    (widened to the coarse cell size in large boxes; rho is the slab average);
    ellipticals use spherical radius. The frame moves with the mass-weighted
    particle velocity within one r_scale of the centre.

    Returns dict of numpy arrays: r_kpc, x (= r / r_scale), rho (shell mass /
    shell volume, Msun/pc**3), T (K), P (dyn/cm**2), vmag (km/s), vrot
    (v_phi, disks) or vr (v_r, ellipticals) (km/s), menc (Msun, enclosed in
    the slab / sphere), plus r_scale_kpc (default: the galaxy's)."""
    L, rs = g["L"], r_scale_kpc or g["r_scale_kpc"]
    if center is None:
        center = profile_centers(g)[0]
    c3 = np.array([center[0], center[1], 0.0])
    rmax_kpc = rmax if rmax is not None else min(0.95 * g["hw"] * L, 5.0 * rs)
    rmin_kpc = rmin_frac * rs

    # frame velocity: particles near the centre
    d = g["pos"] - c3
    near = (np.hypot(d[:, 0], d[:, 1]) if not g["spheroid"] else np.linalg.norm(d, axis=1)) * L < rs
    if not g["spheroid"]:
        near &= np.abs(d[:, 2]) * L < 3 * g["h_kpc"]
    bulk = np.average(g["vel"][near], axis=0, weights=g["mass"][near])

    # slab half-height: 3 z0, but at least ~one coarse cell (1/32 of the box)
    # so the cell centres of coarse grids next to the midplane are included
    slab = max(3 * g["h_kpc"], 0.6 * 2 * g["hw"] * L / 32)
    c = ds.arr(c3, "code_length")
    if g["spheroid"]:
        src = ds.sphere(c, (rmax_kpc, "kpc"))
        rfield = ("index", "spherical_radius")
        vfield = ("gas", "velocity_spherical_radius")
    else:
        src = ds.disk(c, [0, 0, 1], (rmax_kpc, "kpc"), (slab, "kpc"))
        rfield = ("index", "cylindrical_radius")
        vfield = ("gas", "velocity_cylindrical_theta")
    if getattr(ds, "_sph_ptypes", ()):  # SPH: bin on the particle-based field
        rfield = ("gas", rfield[1])
    src.set_field_parameter("bulk_velocity", ds.arr(bulk, "km/s"))
    src.set_field_parameter("normal", ds.arr([0, 0, 1], ""))

    kw = dict(n_bins=nbins, extrema={rfield: (rmin_kpc, rmax_kpc)}, logs={rfield: True},
              units={rfield: "kpc"})
    pm = yt.create_profile(src, rfield, [("gas", "mass")], weight_field=None, **kw)
    pw = yt.create_profile(src, rfield, [("gas", "temperature"), ("gas", "pressure"),
                                         ("gas", "velocity_magnitude"), vfield],
                           weight_field=("gas", "mass"), **kw)
    edges = pm.x_bins.to("kpc").d
    r = pm.x.to("kpc").d
    mass = pm["gas", "mass"].to("Msun").d
    if g["spheroid"]:
        vol = 4.0 / 3.0 * np.pi * np.diff(edges**3)
    else:
        vol = np.pi * np.diff(edges**2) * 2 * slab
    m_in = float(src.include_inside(rfield, 0, rmin_kpc, units="kpc")["gas", "mass"].sum().to("Msun"))
    out = dict(r_kpc=r, x=r / rs, rho=mass / vol / 1e9,
               T=pw["gas", "temperature"].to("K").d,
               P=pw["gas", "pressure"].to("dyn/cm**2").d,
               vmag=pw["gas", "velocity_magnitude"].to("km/s").d,
               menc=m_in + np.cumsum(mass), edges_kpc=edges, r_scale_kpc=rs)
    out["vr" if g["spheroid"] else "vrot"] = pw[vfield].to("km/s").d
    for k in ("T", "P", "vmag", "vrot", "vr", "rho"):
        if k in out:
            out[k] = np.where(mass > 0, out[k], np.nan)
    return out


def dynamical_mass_profile(g, r_kpc):
    """Total (dynamical) mass enclosed within r in the model potential:
    logarithmic potential for disks, the Plummer sphere itself for
    ellipticals; None for the mergers (two moving potentials)."""
    d = g["d"]
    u = d["units"]
    r = np.asarray(r_kpc) / u["L_kpc"]
    if g["spheroid"]:
        a = d["a"]
        return u["M_Msun"] * r**3 / (r**2 + a**2) ** 1.5
    if "V0" in d:
        return u["M_Msun"] * mdyn_log_code(r, d["V0"])
    return None


def plot_profiles(profiles, g, label=""):
    """3 x 2 panel of radial profiles, one line per representation:
    density, temperature, pressure, |v|, rotation (or radial) velocity and
    enclosed mass (with the model's dynamical mass for comparison). x axis
    is r / r_scale, with kpc on the top axes. `profiles` is {method: dict}
    from radial_profiles. Returns the figure."""
    import matplotlib.pyplot as plt

    rs = next(iter(profiles.values()))["r_scale_kpc"]
    sph = g["spheroid"]
    vkey = "vr" if sph else "vrot"
    panels = [
        ("rho", r"$\rho$ [M$_\odot$ pc$^{-3}$]", True),
        ("T", "T [K] (mass-weighted)", True),
        ("P", r"P [dyn cm$^{-2}$] (mass-weighted)", True),
        ("vmag", "|v| [km/s] (mass-weighted)", False),
        (vkey, ("v$_r$" if sph else r"v$_\phi$") + " [km/s] (mass-weighted)", False),
        ("menc", r"M(<r) [M$_\odot$]", True),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.2))
    styles = {"sph": "-", "uniform": "--", "amr": ":"}
    for ax, (key, ylabel, logy) in zip(axes.flat, panels):
        for m, p in profiles.items():
            ax.plot(p["x"], p[key], styles.get(m, "-"), lw=2, label=METHOD_LABELS[m])
        if key == "menc":
            any_p = next(iter(profiles.values()))
            mdyn = dynamical_mass_profile(g, any_p["r_kpc"])
            if mdyn is not None:
                ax.plot(any_p["x"], mdyn, "k-", lw=1, alpha=0.6, label="dynamical (model)")
        ax.set_xscale("log")
        if logy:
            ax.set_yscale("log")
        ax.set_ylabel(ylabel)
        ax.set_xlabel(f"r / r_scale  (r_scale = {rs:.2f} kpc, "
                      f"{'R_e' if sph else 'r_d'}; {'spherical' if sph else 'cylindrical'})")
        top = ax.secondary_xaxis("top", functions=(lambda x: x * rs, lambda r: r / rs))
        top.set_xlabel("r [kpc]")
        ax.grid(alpha=0.3, which="both")
    axes.flat[0].legend()
    axes.flat[-1].legend()
    fig.suptitle(f"{g['title']} radial profiles {label}".strip(), fontsize=13)
    fig.tight_layout()
    return fig


# =====================================================================
# 3D pseudo-AMR
# =====================================================================
# Cubic cells, refine_by = 2 in x, y and z. A block at level L is split into
# level L + 1 if:
#   L = 0:  its centre lies inside an ellipsoid around the box centre (blocks
#           outside stay coarse, giving an ellipsoidal coarse shell).
#   L >= 1: its centre lies inside an oblate ellipsoid around a nucleus
#           (semi-axes R_L in the plane, Z_L = min(R_L, z_extent_h * h_vert)
#           vertically, measured to the block face nearest the midplane), or
#           the mean density of its densest horizontal layer exceeds
#           density_frac[L-1] * peak fine-cell density.
#   L = max_level (central level): inside a small oblate ellipsoid around a
#           nucleus; split once more by conservative prolongation.
# The oblate region follows the thin disk instead of wasting cells above and
# below it; the density criterion follows arms, bars and tidal tails.
# Spheroids use spherical=True (sphere inside a cubic box).

def find_centers(
    x, y, xlim, ylim, n_centers=1, resolution=48, smooth_sigma=3.0,
    min_sep_frac=0.2, min_prominence=0.08, weights=None,
):
    """Up to n_centers density peaks (galaxy nuclei), strongest first.

    A minimum prominence (relative to the global maximum) rejects noise
    peaks. A minimum separation (fraction of the box width) merges clumpy
    double nuclei into one centre while still resolving real galaxy pairs, and
    guarantees no two centres coincide (which would give degenerate AMR
    patches). Optional particle `weights` (masses)."""
    H, _, _ = np.histogram2d(x, y, bins=resolution, range=[xlim, ylim], weights=weights)
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
    fine_fields, xlim, ylim, zlim, centers, h_vert, particles=None, particle_mass=None,
    base_res=16, base_nz=4, block=2, block0=1, refine_by=2, max_level=4,
    envelope_axes_frac=None, radii_frac=(0.55, 0.32, 0.16), density_frac=(0.02, 0.05, 0.12),
    z_extent_h=3.0, spherical=False, central_radius_frac=0.07, central_sigma=1.5,
):
    """Build the list of 3D patches (root + nested refined blocks).

    fine_fields is a dict of extensive master fields ("mass" drives the
    refinement), each of shape (N, N, Nz) with N = base_res * refine_by**max_level
    and Nz = base_nz * refine_by**max_level (cubic cells). Patches are carved
    out by exact integer block sums in pixel-index space, so each patch sits
    on its parent's cell edges (required by yt) and every extensive quantity
    is conserved. Each patch carries a "fields" dict.

    block0 : block size (root cells) for the level 0 -> 1 split; 1 gives a
        smooth ellipsoid boundary.
    envelope_axes_frac : semi-axes (fractions of the half-width) of the
        level-0 refinement ellipsoid. Default (0.9, 0.9, 0.9) for spheroids,
        (0.9, 0.9, 0.5) for disks.
    radii_frac, density_frac : per-level nucleus-region radius (fraction of
        the half-width) and density threshold (fraction of the peak fine-cell
        value). The last entry is reused for deeper levels.
    z_extent_h : vertical half-thickness of the refinement ellipsoid in scale
        heights (ignored if spherical).
    central_radius_frac : if `particles` (n, 3) is given, one extra level
        (max_level + 1) is added inside this radius of each nucleus. It is a
        conservative prolongation: each parent cell's content (every field
        alike) is split among its children in proportion to a lightly
        smoothed (central_sigma, in child cells) particle mass histogram.
        None disables it.
    """
    R = refine_by
    fine_field = fine_fields["mass"]
    master = (base_res * R**max_level, base_res * R**max_level, base_nz * R**max_level)
    assert fine_field.shape == master, f"fine fields must be {master}, got {fine_field.shape}"
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

    def carve(sl, factor):
        return {k: block_reduce_sum(v[sl], factor) for k, v in fine_fields.items()}

    patches = [dict(
        left_edge=[xlim[0], ylim[0], zlim[0]], right_edge=[xlim[1], ylim[1], zlim[1]],
        level=0, dims=(base_res, base_res, base_nz),
        fields=carve((slice(None),) * 3, R**max_level),
    )]

    def vertical_gap(cz, k0, k1):
        """Distance from the midplane to the nearest face of the block (0 if
        the block straddles it), so blocks taller than a thin disk still
        count as containing it."""
        return 0.0 if spherical else max(abs(cz) - 0.5 * (k1 - k0) * dz, 0.0)

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
        zz = vertical_gap(cz, k0, k1)
        for ccx, ccy in centers:
            if ((cx - ccx) ** 2 + (cy - ccy) ** 2) / R_L**2 + (zz / Z_L) ** 2 <= 1.0:
                return True
        # densest horizontal layer of the block, so a thin disk inside a
        # tall block is not diluted by the empty cells above and below it
        layer_density = fine_field[i0:i1, j0:j1, k0:k1].mean(axis=(0, 1)).max()
        return layer_density > density_frac[min(idx, len(density_frac) - 1)] * ref_value

    R_c = central_radius_frac * hw if central else 0.0
    Z_c = R_c if spherical else min(R_c, z_extent_h * h_vert)

    def in_central(i0, i1, j0, j1, k0, k1):
        cx = xlim[0] + 0.5 * (i0 + i1) * dx
        cy = ylim[0] + 0.5 * (j0 + j1) * dy
        zz = vertical_gap(zlim[0] + 0.5 * (k0 + k1) * dz, k0, k1)
        return any(((cx - ccx) ** 2 + (cy - ccy) ** 2) / R_c**2 + (zz / Z_c) ** 2 <= 1.0
                   for ccx, ccy in centers)

    def prolong(i0, i1, j0, j1, k0, k1):
        """Split each cell of the fine fields [i0:i1, j0:j1, k0:k1] into R^3
        children, weighted by the local particle mass histogram."""
        nchild = block * R
        dxf, dyf, dzf = dx / R, dy / R, dz / R
        lo = np.array([xlim[0] + i0 * dx - pad * dxf, ylim[0] + j0 * dy - pad * dyf,
                       zlim[0] + k0 * dz - pad * dzf])
        hi = np.array([xlim[0] + i1 * dx + pad * dxf, ylim[0] + j1 * dy + pad * dyf,
                       zlim[0] + k1 * dz + pad * dzf])
        m = np.all((particles > lo) & (particles < hi), axis=1)
        nb = nchild + 2 * pad
        H, _ = np.histogramdd(particles[m], bins=(nb, nb, nb), range=list(zip(lo, hi)),
                              weights=None if particle_mass is None else particle_mass[m])
        S = gaussian_filter(H, sigma=central_sigma)[pad:-pad, pad:-pad, pad:-pad]
        S = S + 0.1 * S.mean() + 1e-12  # floor: empty regions split evenly
        S6 = S.reshape(block, R, block, R, block, R)
        w = S6 / S6.sum(axis=(1, 3, 5), keepdims=True)
        out = {}
        for k, v in fine_fields.items():
            parent = v[i0:i1, j0:j1, k0:k1]
            out[k] = (w * parent[:, None, :, None, :, None]).reshape(nchild, nchild, nchild)
        return out

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
                        fchild = carve((slice(pi0, pi1), slice(pj0, pj1), slice(pk0, pk1)),
                                       child_factor)
                    else:  # central prolongation level
                        if not in_central(pi0, pi1, pj0, pj1, pk0, pk1):
                            continue
                        fchild = prolong(pi0, pi1, pj0, pj1, pk0, pk1)
                    patches.append(dict(
                        left_edge=[xlim[0] + pi0 * dx, ylim[0] + pj0 * dy, zlim[0] + pk0 * dz],
                        right_edge=[xlim[0] + pi1 * dx, ylim[0] + pj1 * dy, zlim[0] + pk1 * dz],
                        level=level + 1, dims=fchild["mass"].shape, fields=fchild,
                    ))
                    if level < max_level:
                        recurse(level + 1, pi0, pi1, pj0, pj1, pk0, pk1)

    recurse(0, 0, master[0], 0, master[1], 0, master[2])
    return patches


def build_amr_ds(
    pos, hsml, ext, xlim, ylim, units, n_centers=1, base_res=16, block=2, refine_by=2,
    max_level=4, h_vert=H_VERT, spherical=False, **crit_kwargs,
):
    """Particles -> 3D pseudo-AMR yt dataset. Returns (ds, info).

    Levels 0..max_level come from the deposited master fields; one more
    central level resolves the nuclei (see build_amr_patches_3d). The box is
    a cube: base_nz = base_res root cells in z."""
    R = refine_by
    L = units["L_kpc"]
    dx0 = (xlim[1] - xlim[0]) / base_res
    base_nz = base_res
    zhalf = 0.5 * base_nz * dx0
    zlim = (-zhalf, zhalf)
    shape = (base_res * R**max_level, base_res * R**max_level, base_nz * R**max_level)

    fine_fields = deposit_fields(pos, hsml, ext, xlim, ylim, zlim, shape)
    centers = find_centers(pos[:, 0], pos[:, 1], xlim, ylim, n_centers=n_centers,
                           weights=ext["mass"])

    patches = build_amr_patches_3d(
        fine_fields, xlim, ylim, zlim, centers, h_vert,
        particles=pos, particle_mass=ext["mass"], base_res=base_res,
        base_nz=base_nz, block=block, refine_by=R, max_level=max_level,
        spherical=spherical, **crit_kwargs,
    )
    del fine_fields

    grid_data = []
    for p in patches:
        gd = dict(left_edge=p["left_edge"], right_edge=p["right_edge"], level=p["level"],
                  dimensions=list(p["dims"]))
        gd.update(fields_to_yt(p["fields"], p["left_edge"], p["right_edge"], L))
        grid_data.append(gd)

    bbox = np.array([xlim, ylim, list(zlim)])
    ds = yt.load_amr_grids(grid_data, [base_res, base_res, base_nz], bbox=bbox, refine_by=R,
                           **yt_units(units))
    counts = np.bincount([p["level"] for p in patches])
    finest_pc = (xlim[1] - xlim[0]) * L * 1000 / (base_res * R ** (len(counts) - 1))
    return ds, dict(n_grids=len(patches), grids_per_level=counts.tolist(),
                    finest_cell_pc=round(finest_pc, 1))


# =====================================================================
# Galaxy definitions
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
    "10_dwarf_irregular": dict(
        specs=[], kw=dict(rd=0.15, r0_max=0.5, V0=0.5, vdisp=0.05, vphi_scatter=0.15),
        n=25000, dwarf=True,
    ),
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
    # grand-design primary + smaller, faster secondary on a retrograde spin:
    # the two galaxies respond differently, one raising a much stronger tail
    "8_interacting_pair": dict(
        specs1=MERGER_SPEC_BIG, specs2=MERGER_SPEC_SMALL,
        rd1=0.32, rd2=0.15, r0_max1=1.5, r0_max2=0.85,
        n1=18000, n2=9000, spin1=1.0, spin2=-1.0, m1_frac=0.65,
        t_total=3.6, GM_tot=4.5, r12_init=(-4.5, 1.2), v12_init=(0.95, -0.75),
        lim=3,
    ),
    # primary + smaller barred secondary, both prograde: deeper, more
    # disruptive encounter toward coalescence
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
    """Run the dynamics / sampling for all galaxies, or only those in `only`.

    Returns {name: dict(x, y, [z], vx, vy, [vz], mass, xlim, ylim, n_centers,
    units, h_vert, r_scale, M_dyn_Msun, r_dyn_kpc, t_Myr, ...)}. Positions
    and velocities are in the galaxy's code units, masses in Msun."""
    want = lambda name: only is None or name in only
    results = {}
    dyn_t = 9.0

    for name, cfg in DISK_SPECS.items():
        if not want(name):
            continue
        ph = PHYS[name]
        u = make_units(**ph)
        n = cfg.get("n", n_disk)
        kw = cfg["kw"]
        V0, rd = kw.get("V0", 1.0), kw.get("rd", 0.35)
        x, y, vx, vy = simulate_disk_galaxy(cfg["specs"], n=n, t_total=dyn_t, dt=0.004, **kw)
        m_part = disk_mass_code(V0, rd, ph["f_disk"]) / n * u["M_Msun"]
        results[name] = dict(
            x=x, y=y, vx=vx, vy=vy, mass=np.full(x.shape[0], m_part),
            xlim=(-1, 1), ylim=(-1, 1), n_centers=1, dwarf=cfg.get("dwarf", False),
            units=u, V0=V0, h_vert=ph["h_kpc"] / u["L_kpc"], r_scale=rd,
            M_dyn_Msun=mdyn_log_code(1.0, V0) * u["M_Msun"], r_dyn_kpc=u["L_kpc"],
            t_Myr=dyn_t * u["t_Myr"], version=CACHE_VERSION,
        )

    for name, cfg in ELLIPTICALS.items():
        if not want(name):
            continue
        u = make_units(**PHYS[name])
        x, y, z, vx, vy, vz = simulate_elliptical(n=n_ellip, **cfg)
        a = cfg["a"]
        results[name] = dict(
            x=x, y=y, z=z, vx=vx, vy=vy, vz=vz,
            mass=np.full(x.shape[0], 0.995 * u["M_Msun"] / n_ellip),  # u < 0.995 sampled
            xlim=(-1, 1), ylim=(-1, 1), n_centers=1, units=u, a=a, r_scale=a,
            M_dyn_Msun=u["M_Msun"], r_dyn_kpc=np.inf,  # self-gravitating, no dark matter
            t_Myr=0.0, version=CACHE_VERSION,
        )

    for name, cfg in MERGERS.items():
        if not want(name):
            continue
        ph = PHYS[name]
        u = make_units(**ph)
        cfg = dict(cfg)
        lim = cfg.pop("lim")
        specs1 = cfg.pop("specs1")
        specs2 = cfg.pop("specs2")
        x, y, vx, vy, owner, info = simulate_merger(specs1, specs2, r_kill=12, **cfg)
        V01, V02 = info["V0"]
        m1 = disk_mass_code(V01, cfg["rd1"], ph["f_disk"]) / cfg["n1"]
        m2 = disk_mass_code(V02, cfg["rd2"], ph["f_disk"]) / cfg["n2"]
        results[name] = dict(
            x=x, y=y, vx=vx, vy=vy, owner=owner, pair=info,
            mass=np.where(owner == 0, m1, m2) * u["M_Msun"],
            xlim=(-lim, lim), ylim=(-lim, lim), n_centers=2, units=u,
            V0_pair=(V01, V02), h_vert=ph["h_kpc"] / u["L_kpc"], r_scale=cfg["rd1"],
            r_scale_pair=(cfg["rd1"], cfg["rd2"]),
            M_dyn_Msun=cfg["GM_tot"] * u["M_Msun"], r_dyn_kpc=np.inf,
            t_Myr=cfg["t_total"] * u["t_Myr"], version=CACHE_VERSION,
        )

    return results


def load_or_generate(path, only=None, n_disk=60000):
    """Cached particles from `path`, (re)generating missing galaxies and any
    written by an older CACHE_VERSION."""
    try:
        with open(path, "rb") as f:
            results = pickle.load(f)
    except (FileNotFoundError, EOFError):
        results = {}
    results = {k: v for k, v in results.items() if v.get("version") == CACHE_VERSION}
    wanted = only if only else ORDER
    missing = [n for n in wanted if n not in results]
    if missing:
        results.update(generate_all(n_disk=n_disk, only=missing))
        with open(path, "wb") as f:
            pickle.dump(results, f)
    return results


def mass_budget(results):
    """One row per galaxy: units, particle number and mass, baryonic mass in
    the box, dynamical mass, baryon fraction, scale height and scale length.
    Returns a list of dicts (print with print_mass_budget)."""
    rows = []
    for name in ORDER:
        if name not in results:
            continue
        d = results[name]
        u = d["units"]
        inside = ((d["x"] > d["xlim"][0]) & (d["x"] < d["xlim"][1]) &
                  (d["y"] > d["ylim"][0]) & (d["y"] < d["ylim"][1]))
        M_b = d["mass"][inside].sum()
        rows.append(dict(
            name=name, L_kpc=u["L_kpc"], V_kms=u["V_kms"], t_unit_Myr=u["t_Myr"],
            box_kpc=(d["xlim"][1] - d["xlim"][0]) * u["L_kpc"], N=d["x"].shape[0],
            m_part=float(np.median(d["mass"])), M_baryon=M_b, M_dyn=d["M_dyn_Msun"],
            r_dyn_kpc=d["r_dyn_kpc"], f_b=M_b / d["M_dyn_Msun"],
            h_pc=d.get("h_vert", np.nan) * u["L_kpc"] * 1000,
            r_scale_kpc=d["r_scale"] * u["L_kpc"], age_Myr=d["t_Myr"],
        ))
    return rows


def print_mass_budget(rows):
    hdr = (f"{'galaxy':<20}{'box':>6}{'V':>6}{'t_u':>6}{'N':>7}{'m_part':>9}{'M_bary':>9}"
           f"{'M_dyn':>9}{'(<r)':>6}{'f_b':>6}{'z0':>6}{'r_s':>6}{'age':>6}")
    print(hdr)
    print(f"{'':<20}{'kpc':>6}{'km/s':>6}{'Myr':>6}{'':>7}{'Msun':>9}{'Msun':>9}"
          f"{'Msun':>9}{'kpc':>6}{'':>6}{'pc':>6}{'kpc':>6}{'Myr':>6}")
    for r in rows:
        rd = "all" if not np.isfinite(r["r_dyn_kpc"]) else f"{r['r_dyn_kpc']:.0f}"
        h = "-" if not np.isfinite(r["h_pc"]) else f"{r['h_pc']:.0f}"
        print(f"{r['name']:<20}{r['box_kpc']:>6.0f}{r['V_kms']:>6.0f}{r['t_unit_Myr']:>6.0f}"
              f"{r['N']:>7d}{r['m_part']:>9.2e}{r['M_baryon']:>9.2e}{r['M_dyn']:>9.2e}"
              f"{rd:>6}{r['f_b']:>6.2f}{h:>6}{r['r_scale_kpc']:>6.2f}{r['age_Myr']:>6.0f}")


# =====================================================================
# Rendering
# =====================================================================

METHODS = ("sph", "uniform", "amr")
METHOD_LABELS = {"sph": "SPH", "uniform": "Uniform grid", "amr": "Pseudo-AMR"}


def prepare_galaxy(name, d):
    """Cached simulation result -> full 3D particle set in physical units
    plus box info, for build_ds(). Disks get z from a sech^2 profile and vz
    from the matching isothermal sheet; the dwarf's masses get its clumpy
    texture. Disk velocities come from the kinematic model
    (disk_kinematics); ellipticals keep their sampled velocities. Seeded from
    the galaxy name, so repeated calls give identical particles.

    Positions stay in code units ("pos"), velocities are in km/s ("vel"),
    masses in Msun ("mass")."""
    xlim, ylim = d["xlim"], d["ylim"]
    u = d["units"]
    L, V = u["L_kpc"], u["V_kms"]
    hw = 0.5 * (xlim[1] - xlim[0])
    spheroid = "z" in d
    zhalf = z_half_extent(hw, spheroid)
    zlim = (-zhalf, zhalf)
    h_vert = d.get("h_vert", H_VERT)
    seed = zlib.crc32(name.encode())
    zrng = np.random.default_rng(seed)

    x, y = d["x"], d["y"]
    z = d["z"] if spheroid else add_z_disk(x, y, h=h_vert, rng=zrng)
    mass = d["mass"].copy()

    if d.get("dwarf"):  # clumpiness texture on the masses (total unchanged)
        nres = 220
        noise = fractal_noise((nres, nres), zrng, octaves=4, base_sigma=6.0, persistence=0.6)
        xi = np.clip(((x - xlim[0]) / (xlim[1] - xlim[0]) * nres).astype(int), 0, nres - 1)
        yi = np.clip(((y - ylim[0]) / (ylim[1] - ylim[0]) * nres).astype(int), 0, nres - 1)
        weight = np.exp(0.9 * noise[xi, yi])
        mass *= weight / weight.mean()

    if spheroid:
        vel = np.column_stack([d["vx"], d["vy"], d["vz"]])
    elif "pair" in d:
        p = d["pair"]
        # tidal debris follows the simulated flow, smoothed over 64 neighbours
        # to drop the random motion pumped in by the strong disk forcing
        xy = np.column_stack([x, y])
        _, nb = cKDTree(xy).query(xy, k=64)
        v_flow = np.column_stack([d["vx"][nb].mean(axis=1), d["vy"][nb].mean(axis=1)])
        hosts = [dict(center=p["centers"][i], vcenter=p["vcenters"][i], V0=p["V0"][i],
                      spin=p["spins"][i], sel=d["owner"] == i, r_in=p["r0_max"][i],
                      r_out=2 * p["r0_max"][i], v_sim=v_flow[d["owner"] == i])
                 for i in (0, 1)]
        vel = disk_kinematics(x, y, mass / u["M_Msun"], h_vert, hosts, zrng)
    else:
        hosts = [dict(center=(0.0, 0.0), vcenter=(0.0, 0.0), V0=d["V0"], spin=1.0,
                      sel=np.ones(x.shape[0], bool))]
        vel = disk_kinematics(x, y, mass / u["M_Msun"], h_vert, hosts, zrng)
    pos = np.column_stack([x, y, z])
    vel = vel * V

    return dict(name=name, d=d, pos=pos, vel=vel, mass=mass, units=u, L=L,
                xlim=xlim, ylim=ylim, zlim=zlim, hw=hw, zhalf=zhalf, h_vert=h_vert,
                h_kpc=h_vert * L, spheroid=spheroid, seed=seed, title=TITLES[name],
                r_scale_kpc=d["r_scale"] * L, _sph=None, _s2=None)


def build_ds(g, method, fields="all", amr_kwargs=None):
    """yt dataset for one representation ("sph", "uniform", "amr") of a
    prepared galaxy. fields="all" gives density, velocity, temperature and
    pressure; "density" only the density (faster, used by the montages).
    Returns (ds, info); info is the AMR summary for "amr", else None."""
    if method == "sph":
        return build_sph_ds(g, fields), None

    # grids only need the particles inside the visible box
    pos = g["pos"]
    xlim, ylim, zlim = g["xlim"], g["ylim"], g["zlim"]
    inside = ((pos[:, 0] > xlim[0]) & (pos[:, 0] < xlim[1]) &
              (pos[:, 1] > ylim[0]) & (pos[:, 1] < ylim[1]))
    hsml = particle_sph(g)[0][inside]  # kNN lengths from the full set
    ext = extensive_weights(g, inside, fields)

    if method == "uniform":
        return build_uniform_ds(pos[inside], hsml, ext, xlim, ylim, zlim, g["units"]), None

    kw = dict(amr_kwargs or {})
    if g["spheroid"]:
        kw.setdefault("max_level", 3)  # keeps the 3D master field small
    return build_amr_ds(pos[inside], hsml, ext, xlim, ylim, g["units"],
                        n_centers=g["d"].get("n_centers", 1), h_vert=g["h_vert"],
                        spherical=g["spheroid"], **kw)


def render_galaxy(name, d, outdir=".", individual=False, grids=False, amr_kwargs=None,
                  res=512, kind="slic"):
    """Render one galaxy in all three representations, face-on and edge-on,
    as both density slices and projections (one dataset build serves both).

    Returns {method: {"face_slic", "edge_slic", "face_proj", "edge_proj"},
    "hw", "zhalf", "L_kpc", "amr_info"}; slices in Msun/kpc**3, projections
    in Msun/kpc**2. `kind` selects which one the per-galaxy yt PNGs
    (individual=True) use."""
    g = prepare_galaxy(name, d)
    out = dict(hw=g["hw"], zhalf=g["zhalf"], L_kpc=g["L"], version=CACHE_VERSION)
    for method in METHODS:
        ds, info = build_ds(g, method, fields="density", amr_kwargs=amr_kwargs)
        arrs = {}
        for k in KINDS:
            arrs[f"face_{k}"], arrs[f"edge_{k}"] = image_arrays(ds, g["hw"], g["zhalf"], res, k)
        out[method] = arrs
        if info is not None:
            out["amr_info"] = info
        if individual:
            save_yt_plots(ds, f"{outdir}/zoo_{name}_{method}",
                          f"{g['title']} ({METHOD_LABELS[method]})",
                          grids=grids and method == "amr", kind=kind)
        del ds
    return out


def make_montages(frames, outdir=".", dyn_range=3.0, kind="slic", show=False):
    """Write zoo_montage_faceon_{kind}.png and zoo_montage_edgeon_{kind}.png.

    Rows = SPH / uniform / AMR, in blocks of five galaxies. Log colour scale
    spanning dyn_range dex, shared by the three methods for each galaxy and
    view. Returns {"faceon": fig, "edgeon": fig}; show=True keeps the figures
    open for inline display, otherwise they are closed."""
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    names = [n for n in ORDER if n in frames]
    cmap = matplotlib.colormaps["viridis"].copy()
    cmap.set_bad(cmap(0.0))    # empty SPH pixels -> floor colour
    cmap.set_under(cmap(0.0))

    def norm(f, key):
        # vmax = largest 99.9th percentile among the methods, vmin = dyn_range dex below
        vmax = max(np.percentile(f[m][key], 99.9) for m in METHODS)
        return LogNorm(vmin=vmax / 10**dyn_range, vmax=vmax)

    views = {"face": ("faceon", "face-on (z = 0)", "face-on (along z)"),
             "edge": ("edgeon", "edge-on (y = 0)", "edge-on (along y)")}
    blocks = [names[i:i + 5] for i in range(0, len(names), 5)]
    figs = {}
    for view, (label, what_slic, what_proj) in views.items():
        key = f"{view}_{kind}"
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
                    L = f["L_kpc"]
                    e, ez = f["hw"] * L, f["zhalf"] * L
                    ax.imshow(f[method][key], origin="lower", extent=[-e, e, -ez, ez],
                              cmap=cmap, norm=norm(f, key), interpolation="nearest")
                    if m == 0:
                        ax.set_title(f"{TITLES[block[c]]}  ({2 * e:.0f} kpc)", fontsize=10)
                    if c == 0:
                        ax.set_ylabel(METHOD_LABELS[method], fontsize=11)
        what = what_slic if kind == "slic" else what_proj
        fig.suptitle(f"Galaxy zoo, {what} {KIND_LABELS[kind]}", fontsize=14)
        fig.tight_layout(rect=(0, 0, 1, 0.985))
        fig.savefig(f"{outdir}/zoo_montage_{label}_{kind}.png", dpi=110)
        figs[label] = fig
        if not show:
            plt.close(fig)
    return figs


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")  # script mode only; notebooks keep their own backend

    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--cache", default="zoo_cache.pkl",
                    help="pickle of the simulated particles; reused if present (skips the integrators)")
    ap.add_argument("--frb-cache", default="zoo_frames.pkl",
                    help="pickle of the rendered image arrays (re-layout montages without re-rendering)")
    ap.add_argument("--only", nargs="*", help="galaxy names to (re)generate / render")
    ap.add_argument("--n-disk", type=int, default=60000)
    ap.add_argument("--redo", action="store_true", help="ignore the rendered-image cache")
    ap.add_argument("--individual", action="store_true", help="also save labelled per-galaxy yt PNGs")
    ap.add_argument("--grids", action="store_true", help="with --individual: AMR grid outlines")
    ap.add_argument("--res", type=int, default=512)
    ap.add_argument("--mode", choices=(*KINDS, "both"), default="slic",
                    help="montage content: midplane slices (slic), projections (proj), or both")
    ap.add_argument("--budget", action="store_true", help="print the mass budget table")
    args = ap.parse_args()

    results = load_or_generate(args.cache, only=args.only, n_disk=args.n_disk)
    if args.budget:
        print_mass_budget(mass_budget(results))

    try:
        with open(args.frb_cache, "rb") as f:
            frames = {} if args.redo else pickle.load(f)
    except (FileNotFoundError, EOFError):
        frames = {}

    kinds = KINDS if args.mode == "both" else (args.mode,)
    for name in ORDER:
        if name not in results or (args.only and name not in args.only):
            continue
        if (name in frames and not args.individual
                and frames[name].get("version") == CACHE_VERSION):
            continue
        frames[name] = render_galaxy(name, results[name], outdir=args.outdir,
                                     individual=args.individual, grids=args.grids, res=args.res,
                                     kind=kinds[0])
        print(name, "AMR:", frames[name]["amr_info"], flush=True)
        with open(args.frb_cache, "wb") as f:  # incremental: a crash keeps progress
            pickle.dump(frames, f)

    frames = {k: v for k, v in frames.items() if v.get("version") == CACHE_VERSION}
    for k in kinds:
        make_montages(frames, outdir=args.outdir, kind=k)
    print("wrote montages to", args.outdir)
