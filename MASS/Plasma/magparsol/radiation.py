"""
magparsol/radiation.py
-----------------------
Radiation spectra from single-particle (and ensemble) trajectories.

All computation is post-processing of an already-finalized TrajectoryHistory —
no new simulation is required.

Two spectral methods
--------------------
spectrum_fft :
    FFT of the observer-projected transverse acceleration/velocity, resampled
    onto a uniform time grid when needed (adaptive-RK runs).  Cheap: O(S log S).
    Appropriate for periodic/quasi-periodic gyration (cyclotron, synchrotron).

spectrum_retarded :
    Direct numerical evaluation of Jackson's retarded-time integral at each
    requested frequency.  Accurate for high harmonics / strongly relativistic
    beaming, but O(S × n_ω).  Always shows the final result (never animated
    incrementally).

Ensemble spectra
----------------
ensemble_spectrum :
    Incoherent sum of individual particle spectra (physically correct for
    thermal/uncorrelated ensembles where particle separations ≫ wavelength).

Instantaneous power
-------------------
radiated_power :
    Liénard formula evaluated at each stored trajectory point using stored
    velocity and acceleration (acceleration re-derived from the Lorentz-force
    RHS, not finite-differenced from v).
"""

import warnings
import numpy as np
from magparsol.diagnostics import TrajectoryHistory
from magparsol.constants import C, Q_E

# ── Constants ──────────────────────────────────────────────────────────────────
_EPS0 = 8.854_187_817e-12    # vacuum permittivity [F/m]
_PREFACTOR = Q_E**2 / (6 * np.pi * _EPS0 * C**3)   # Larmor prefactor [W·s²/m²]

# Resampling cap: maximum number of uniform grid points for FFT
_N_RESAMPLE_MAX = 2**18   # ~262 144


# ── Uniform-grid resampling (for adaptive-RK or any non-uniform time series) ──

def _resample_uniform(t: np.ndarray, y: np.ndarray,
                      interp_method: str = "linear") -> tuple:
    """Resample an irregularly-spaced time series onto a uniform grid.

    Parameters
    ----------
    t : ndarray, shape (S,)   — possibly non-uniform times
    y : ndarray, shape (S, …) — signal values
    interp_method : "linear" | "cubic"

    Returns
    -------
    t_uni : ndarray, shape (M,)   — uniform times
    y_uni : ndarray, shape (M, …) — resampled signal
    dt_uni : float                 — uniform step size
    """
    dt_min = float(np.min(np.diff(t)))
    T_total = float(t[-1] - t[0])
    N_naive = int(np.ceil(T_total / dt_min)) + 1

    if N_naive > _N_RESAMPLE_MAX:
        dt_uni = T_total / (_N_RESAMPLE_MAX - 1)
        warnings.warn(
            f"Resampling grid capped at {_N_RESAMPLE_MAX} points "
            f"(dt_resample={dt_uni:.3e} s vs dt_min={dt_min:.3e} s). "
            f"Frequencies above {0.5/dt_uni:.3e} Hz may be underresolved. "
            "Increase N_RESAMPLE_MAX or reduce store_dt for full resolution.",
            UserWarning, stacklevel=3,
        )
    else:
        dt_uni = dt_min

    t_uni = np.linspace(float(t[0]), float(t[-1]),
                        int(round(T_total / dt_uni)) + 1)

    if interp_method == "cubic":
        from scipy.interpolate import CubicSpline
        cs = CubicSpline(t, y, axis=0)
        y_uni = cs(t_uni)
    else:
        # linear — sufficient for smooth trajectories
        shape_out = (len(t_uni),) + y.shape[1:]
        y_uni = np.empty(shape_out)
        for idx in np.ndindex(y.shape[1:]):
            sl = (slice(None),) + idx
            y_uni[sl] = np.interp(t_uni, t, y[sl])

    return t_uni, y_uni, dt_uni


def _check_uniform(t: np.ndarray) -> tuple:
    """Return (is_uniform, dt).  Uniform if max step variation < 1%."""
    diffs = np.diff(t)
    dt = float(diffs.mean())
    is_uni = float(diffs.std()) / dt < 0.01 if dt > 0 else True
    return is_uni, dt


# ── Instantaneous radiated power (Liénard formula) ────────────────────────────

def radiated_power(history: TrajectoryHistory,
                   q: np.ndarray,
                   m: np.ndarray,
                   field,
                   relativistic: bool = True) -> np.ndarray:
    """Instantaneous radiated power at each stored trajectory point.

    Uses the Liénard generalization of Larmor's formula:

        P = (q² γ⁶ / 6πε₀c³) [|a|² − |v×a|²/c²]

    which reduces to P = q²|a|²/(6πε₀c³) for v≪c.

    Acceleration is re-derived analytically from the Lorentz force
    (via the stored field evaluated at stored positions), not finite-
    differenced from velocity — this avoids noise from coarse storage.

    Parameters
    ----------
    history : TrajectoryHistory (finalized)
    q : ndarray, shape (N,)
    m : ndarray, shape (N,)
    field : FieldModel — used to evaluate B at stored positions
    relativistic : bool

    Returns
    -------
    P : ndarray, shape (S, N)   [W]
    """
    history._check_finalized()
    S, N, _ = history.r.shape
    P = np.zeros((S, N))

    for s in range(S):
        r_s = history.r[s]          # (N, 3)
        v_s = history.v[s]          # (N, 3)
        t_s = history.t[s]
        B = field(r_s, t_s)         # (N, 3)

        speed2 = np.sum(v_s**2, axis=1)   # (N,)

        if relativistic:
            beta2  = np.clip(speed2 / C**2, 0.0, 1.0 - 1e-15)
            gamma  = 1.0 / np.sqrt(1.0 - beta2)            # (N,)
            gamma6 = gamma**6
            q_m    = q / m
            a = (q_m * (1.0/gamma))[:, None] * np.cross(v_s, B) # (N, 3)
            vcross_a = np.cross(v_s, a)
            P[s] = (q**2 * gamma6 / (6 * np.pi * _EPS0 * C**3)) * (
                np.sum(a**2, axis=1) - np.sum(vcross_a**2, axis=1) / C**2
            )
        else:
            q_m = q / m
            a   = q_m[:, None] * np.cross(v_s, B)   # (N, 3)
            P[s] = (q**2 / (6 * np.pi * _EPS0 * C**3)) * np.sum(a**2, axis=1)

    return P


def total_radiated_energy(history: TrajectoryHistory,
                           q: np.ndarray, m: np.ndarray,
                           field, relativistic: bool = True) -> np.ndarray:
    """Time-integrated radiated energy for each particle.

    Returns
    -------
    W : ndarray, shape (N,)   [J]
    """
    P = radiated_power(history, q, m, field, relativistic)
    return np.trapz(P, history.t, axis=0)


# ── FFT-based spectrum ────────────────────────────────────────────────────────

def spectrum_fft(history: TrajectoryHistory,
                 pid: int = 0,
                 observer: np.ndarray = None,
                 upto: int = None,
                 interp_method: str = "linear",
                 store_dt_warn_period: float = None) -> tuple:
    """FFT-based single-particle radiated power spectrum.

    Computes the power spectral density of the observer-projected transverse
    velocity (proportional to radiated power spectrum for non-relativistic
    and mildly relativistic cases).

    Parameters
    ----------
    history : TrajectoryHistory (finalized)
    pid : int
        Particle index.
    observer : ndarray, shape (3,) or None
        Observer direction unit vector n̂.  None → (0, 0, 1) (z-axis).
    upto : int or None
        Use only history[:upto] — for animated "building up" panels.
        None → use full trajectory.
    interp_method : "linear" | "cubic"
        Interpolation for non-uniform time grids (adaptive RK).
    store_dt_warn_period : float or None
        Gyroperiod for Nyquist check.  If store_dt > gyroperiod/2, warn.

    Returns
    -------
    freqs : ndarray, shape (M,)   [Hz]
    power : ndarray, shape (M,)   [arb. units, proportional to dI/dω]
    """
    history._check_finalized()

    i_end = upto if upto is not None else len(history.t)
    if i_end < 4:
        return np.array([0.0]), np.array([0.0])

    t = history.t[:i_end]
    v = history.v[:i_end, pid, :]   # (S, 3)

    # Observer projection: transverse velocity components (as vectors, not norm)
    # Taking the norm destroys oscillation info for circular motion.
    # Instead FFT each transverse component separately and sum power spectra.
    if observer is None:
        observer = np.array([0.0, 0.0, 1.0])
    n_hat = np.asarray(observer, dtype=float)
    n_hat = n_hat / np.linalg.norm(n_hat)
    # Transverse velocity vector = v - (v·n̂)n̂  shape (S, 3)
    v_dot_n = np.sum(v * n_hat, axis=1, keepdims=True)
    v_trans = v - v_dot_n * n_hat   # (S, 3) — keep as vector

    # Nyquist check against store_dt
    if store_dt_warn_period is not None and len(t) > 1:
        mean_dt = float(np.mean(np.diff(t)))
        if mean_dt > store_dt_warn_period / 2.0:
            warnings.warn(
                f"store_dt ({mean_dt:.3e} s) > gyroperiod/2 "
                f"({store_dt_warn_period/2:.3e} s). "
                "Spectral content above Nyquist was not stored. "
                "Reduce store_dt for an accurate spectrum.",
                UserWarning, stacklevel=2,
            )

    # Resample to uniform grid if needed (adaptive RK produces uneven spacing)
    is_uni, dt = _check_uniform(t)
    if not is_uni:
        t, v_trans, dt = _resample_uniform(t, v_trans, interp_method)

    # FFT each transverse component; sum power spectra (incoherent combination)
    n     = len(t)
    win   = np.hanning(n)
    freqs = np.fft.rfftfreq(n, d=dt)
    power = np.zeros(len(freqs))
    for ax_idx in range(3):
        sig    = v_trans[:, ax_idx]
        if np.max(np.abs(sig)) < 1e-30:
            continue   # skip zero components (e.g. vz=0 for pure xy gyration)
        fft_c  = np.fft.rfft(sig * win)
        power += (np.abs(fft_c)**2) / n

    # Only positive frequencies
    mask = freqs > 0
    return freqs[mask], power[mask]


def _check_spectrum_convergence(spec_prev: np.ndarray, spec_curr: np.ndarray,
                                 freqs: np.ndarray, rtol: float = 0.05) -> bool:
    """Check if FFT spectrum has converged between two updates.

    Compares peak frequency and its half-power width.
    Returns True if converged (relative change < rtol in both).
    """
    def peak_info(p):
        idx = np.argmax(p)
        f_peak = freqs[idx]
        half  = p[idx] / 2.0
        above = np.where(p >= half)[0]
        fwhm  = freqs[above[-1]] - freqs[above[0]] if len(above) > 1 else 0.0
        return f_peak, fwhm

    if spec_prev is None or len(spec_prev) != len(spec_curr):
        return False
    f0, w0 = peak_info(spec_prev)
    f1, w1 = peak_info(spec_curr)
    if f0 == 0:
        return False
    return (abs(f1 - f0) / abs(f0) < rtol) and (abs(w1 - w0) / (abs(w0) + 1e-30) < rtol)


# ── Retarded-integral spectrum ────────────────────────────────────────────────

def spectrum_retarded(history: TrajectoryHistory,
                       pid: int = 0,
                       omega_array: np.ndarray = None,
                       observer: np.ndarray = None,
                       n_omega: int = 512) -> tuple:
    """Retarded-time Fourier integral spectrum (Jackson Ch. 14).

    dI/dω ∝ |∫ n̂×(n̂×v) exp[iω(t − n̂·r/c)] dt|²

    Accurate for high harmonics and strongly relativistic beaming,
    but O(S × n_ω).  Always computed from the full trajectory.

    Parameters
    ----------
    history : TrajectoryHistory (finalized)
    pid : int
    omega_array : ndarray or None
        Angular frequencies [rad/s] at which to evaluate.
        None → auto-range from 0 to 10× estimated gyrofrequency.
    observer : ndarray, shape (3,) or None
        Observer direction n̂.  None → (0, 0, 1).
    n_omega : int
        Number of frequency points when omega_array is None.

    Returns
    -------
    freqs : ndarray, shape (n_omega,)   [Hz]
    power : ndarray, shape (n_omega,)   [arb. units]
    """
    history._check_finalized()

    t = history.t
    r = history.r[:, pid, :]   # (S, 3)
    v = history.v[:, pid, :]   # (S, 3)

    if observer is None:
        observer = np.array([0.0, 0.0, 1.0])
    n_hat = np.asarray(observer, dtype=float)
    n_hat = n_hat / np.linalg.norm(n_hat)

    # Transverse field: n̂×(n̂×v) = (n̂·v)n̂ − v   (only transverse matters)
    v_dot_n  = np.sum(v * n_hat, axis=1, keepdims=True)   # (S,1)
    v_perp   = v_dot_n * n_hat - v                          # (S, 3)

    # Retardation phase: n̂·r/c
    r_dot_n  = np.sum(r * n_hat, axis=1)   # (S,)

    if omega_array is None:
        # Auto-estimate frequency range from velocity oscillation
        is_uni, dt = _check_uniform(t)
        if not is_uni:
            dt = float(np.mean(np.diff(t)))
        f_max = 0.5 / dt
        omega_array = np.linspace(0, 2 * np.pi * f_max, n_omega + 1)[1:]

    freqs  = omega_array / (2 * np.pi)
    power  = np.zeros(len(omega_array))

    for k, omega in enumerate(omega_array):
        phase      = omega * (t - r_dot_n / C)          # (S,)
        integrand  = v_perp * np.exp(1j * phase)[:, None]   # (S, 3)
        integral   = np.trapz(integrand, t, axis=0)     # (3,)
        power[k]   = float(np.real(np.dot(integral, np.conj(integral))))

    return freqs, power


# ── Ensemble spectrum ─────────────────────────────────────────────────────────

def ensemble_spectrum(history: TrajectoryHistory,
                       method: str = "fft",
                       observer: np.ndarray = None,
                       weights: np.ndarray = None,
                       **kwargs) -> tuple:
    """Incoherent sum of individual particle spectra.

    For a thermal/uncorrelated ensemble the total emitted power spectrum
    is the weighted sum of single-particle spectra.

    Parameters
    ----------
    history : TrajectoryHistory (finalized, N>1)
    method : "fft" | "retarded"
    observer : ndarray, shape (3,) or None
    weights : ndarray, shape (N,) or None
        Per-particle weights.  None → equal weights (1/N each).
    **kwargs : passed to spectrum_fft or spectrum_retarded

    Returns
    -------
    freqs : ndarray   [Hz]
    total_power : ndarray
    individual : list of ndarray — per-particle spectra (same freq grid)
    """
    history._check_finalized()
    N = history.r.shape[1]

    if weights is None:
        weights = np.ones(N) / N

    spec_fn = spectrum_fft if method == "fft" else spectrum_retarded
    individual = []
    ref_freqs  = None

    for pid in range(N):
        f, p = spec_fn(history, pid=pid, observer=observer, **kwargs)
        if ref_freqs is None:
            ref_freqs = f
            total = np.zeros_like(p)
        else:
            # Interpolate onto common frequency grid
            p = np.interp(ref_freqs, f, p, left=0.0, right=0.0)
        total += weights[pid] * p
        individual.append(p)

    return ref_freqs, total, individual
