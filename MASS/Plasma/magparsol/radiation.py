"""
Radiated power and spectra, computed from a finished trajectory history.

* :func:`radiated_power`: Liénard power along the orbit, with the
  acceleration taken from the Lorentz force at each stored point rather
  than from differences of the stored velocities.
* :func:`spectrum_fft`: Fourier spectrum of the transverse motion seen by
  an observer. Fast, but without retardation, so it misses relativistic
  harmonics and beaming.
* :func:`spectrum_retarded`: Jackson's retarded-time integral evaluated
  frequency by frequency. Includes the harmonics, but costs S × n_ω.
* :func:`ensemble_spectrum`: incoherent sum of single-particle spectra,
  appropriate for uncorrelated particles.
"""

import warnings
import numpy as np
from magparsol.diagnostics import TrajectoryHistory
from magparsol.constants import C, EPS0

# NumPy 2 renamed np.trapz to np.trapezoid
_trapezoid = getattr(np, "trapezoid", None) or np.trapz

# Largest uniform grid used when resampling an uneven time series
_N_RESAMPLE_MAX = 2**18


def _resample_uniform(t: np.ndarray, y: np.ndarray,
                      interp_method: str = "linear") -> tuple:
    """Interpolate a time series ``y`` (S, ...) at uneven times ``t`` onto a
    uniform grid.

    The grid step is the smallest step of ``t``, unless that would need
    more than ``_N_RESAMPLE_MAX`` points. ``interp_method`` is "linear" or
    "cubic" (scipy). Returns the uniform times, the resampled series and
    the grid step.
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
        shape_out = (len(t_uni),) + y.shape[1:]
        y_uni = np.empty(shape_out)
        for idx in np.ndindex(y.shape[1:]):
            sl = (slice(None),) + idx
            y_uni[sl] = np.interp(t_uni, t, y[sl])

    return t_uni, y_uni, dt_uni


def _check_uniform(t: np.ndarray) -> tuple:
    """Return ``(is_uniform, mean step)``; uniform means steps vary by < 1 %."""
    diffs = np.diff(t)
    dt = float(diffs.mean())
    is_uni = float(diffs.std()) / dt < 0.01 if dt > 0 else True
    return is_uni, dt


def radiated_power(history: TrajectoryHistory,
                   q: np.ndarray,
                   m: np.ndarray,
                   field,
                   relativistic: bool = True) -> np.ndarray:
    """Radiated power [W] at every stored point, shape (S, N).

    Liénard's formula

        P = q² γ⁶ / (6π ε0 c³) [a² - |v × a|² / c²],

    which for v ≪ c (or ``relativistic=False``) is Larmor's
    P = q² a² / (6π ε0 c³). The acceleration a is computed from the
    Lorentz force with ``field`` at the stored positions and velocities.

    Parameters
    ----------
    history : TrajectoryHistory
    q, m : ndarray, shape (N,)
        Charges [C] and masses [kg].
    field : FieldModel
    relativistic : bool
    """
    history._check_finalized()
    S, N, _ = history.r.shape
    P = np.zeros((S, N))

    for s in range(S):
        r_s = history.r[s]
        v_s = history.v[s]
        t_s = history.t[s]
        B, E = field(r_s, t_s)

        speed2 = np.sum(v_s**2, axis=1)

        if relativistic:
            beta2  = np.clip(speed2 / C**2, 0.0, 1.0 - 1e-15)
            gamma  = 1.0 / np.sqrt(1.0 - beta2)
            gamma6 = gamma**6
            q_m    = q / m
            vdotE  = np.sum(v_s * E, axis=1)
            # a = (q / γm) (E + v × B - (v·E) v / c²)
            a = (q_m * (1.0/gamma))[:, None] * (
                E + np.cross(v_s, B)
                - (vdotE / C**2)[:, None] * v_s
            )
            vcross_a = np.cross(v_s, a)
            P[s] = (q**2 * gamma6 / (6 * np.pi * EPS0 * C**3)) * (
                np.sum(a**2, axis=1) - np.sum(vcross_a**2, axis=1) / C**2
            )
        else:
            q_m = q / m
            a   = q_m[:, None] * (E + np.cross(v_s, B))
            P[s] = (q**2 / (6 * np.pi * EPS0 * C**3)) * np.sum(a**2, axis=1)

    return P


def total_radiated_energy(history: TrajectoryHistory,
                           q: np.ndarray, m: np.ndarray,
                           field, relativistic: bool = True) -> np.ndarray:
    """Radiated energy [J] over the whole run, shape (N,): the time integral
    of :func:`radiated_power` (same parameters)."""
    P = radiated_power(history, q, m, field, relativistic)
    return _trapezoid(P, history.t, axis=0)


def spectrum_fft(history: TrajectoryHistory,
                 pid: int = 0,
                 observer: np.ndarray = None,
                 upto: int = None,
                 interp_method: str = "linear",
                 store_dt_warn_period: float = None,
                 weighting: str = "acceleration") -> tuple:
    """Power spectrum of one particle's transverse motion, by FFT.

    The velocity component perpendicular to the line of sight,
    v⊥ = v - (v·n̂) n̂, is Hann-windowed and Fourier transformed; the power
    of the three components is summed. With ``weighting="acceleration"``
    the result is multiplied by ω², giving the spectrum of the transverse
    acceleration, which is proportional to the dipole radiation spectrum
    dI/dω for slow particles. Retardation is neglected.

    Parameters
    ----------
    history : TrajectoryHistory
    pid : int
        Index of the particle.
    observer : array_like, shape (3,) or None
        Direction of the observer, by default +z.
    upto : int or None
        Use only the first ``upto`` samples (for spectra that build up in
        an animation). None uses the whole run.
    interp_method : "linear" or "cubic"
        Interpolation used to resample uneven time steps (adaptive runs).
    store_dt_warn_period : float or None
        Gyroperiod [s]; warns if the sampling is too coarse to resolve it.
    weighting : "acceleration" or "velocity"

    Returns
    -------
    freqs : ndarray
        Positive frequencies [Hz].
    power : ndarray
        Power in arbitrary units.
    """
    history._check_finalized()

    i_end = upto if upto is not None else len(history.t)
    if i_end < 4:
        return np.array([0.0]), np.array([0.0])

    t = history.t[:i_end]
    v = history.v[:i_end, pid, :]

    if observer is None:
        observer = np.array([0.0, 0.0, 1.0])
    n_hat = np.asarray(observer, dtype=float)
    n_hat = n_hat / np.linalg.norm(n_hat)
    # Keep v⊥ as a vector: its magnitude is constant for circular motion
    # and carries no frequency information.
    v_dot_n = np.sum(v * n_hat, axis=1, keepdims=True)
    v_trans = v - v_dot_n * n_hat

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

    is_uni, dt = _check_uniform(t)
    if not is_uni:
        t, v_trans, dt = _resample_uniform(t, v_trans, interp_method)

    n     = len(t)
    win   = np.hanning(n)
    freqs = np.fft.rfftfreq(n, d=dt)
    power = np.zeros(len(freqs))
    for ax_idx in range(3):
        sig    = v_trans[:, ax_idx]
        if np.max(np.abs(sig)) < 1e-30:
            continue
        fft_c  = np.fft.rfft(sig * win)
        power += (np.abs(fft_c)**2) / n

    if weighting == "acceleration":
        power *= (2 * np.pi * freqs)**2

    mask = freqs > 0
    return freqs[mask], power[mask]


def spectrum_retarded(history: TrajectoryHistory,
                       pid: int = 0,
                       omega_array: np.ndarray = None,
                       observer: np.ndarray = None,
                       n_omega: int = 512,
                       window: str = "hann") -> tuple:
    """Spectrum seen by a distant observer, from the retarded-time integral
    (Jackson, Eq. 14.67):

        d²I / dω dΩ ∝ ω² |∫ n̂ × (n̂ × β) exp[iω (t - n̂·r/c)] dt|²

    The retardation phase n̂·r/c turns relativistic gyration into harmonics
    of the gyrofrequency. The integral is evaluated by the trapezoidal rule
    for every frequency, over the whole run.

    Parameters
    ----------
    history : TrajectoryHistory
    pid : int
        Index of the particle.
    omega_array : ndarray or None
        Angular frequencies [rad/s]. By default ``n_omega`` values up to the
        Nyquist frequency of the stored samples.
    observer : array_like, shape (3,) or None
        Direction of the observer, by default +z.
    n_omega : int
    window : "hann" or None
        Taper over the observation time. The Hann window suppresses the
        leakage caused by the abrupt start and end of the trajectory.

    Returns
    -------
    freqs : ndarray
        Frequencies [Hz].
    power : ndarray
        Power in arbitrary units.
    """
    history._check_finalized()

    t = history.t
    r = history.r[:, pid, :]
    v = history.v[:, pid, :]

    if observer is None:
        observer = np.array([0.0, 0.0, 1.0])
    n_hat = np.asarray(observer, dtype=float)
    n_hat = n_hat / np.linalg.norm(n_hat)

    # n̂ × (n̂ × v) = (n̂·v) n̂ - v
    v_dot_n  = np.sum(v * n_hat, axis=1, keepdims=True)
    v_perp   = v_dot_n * n_hat - v

    r_dot_n  = np.sum(r * n_hat, axis=1)

    if window == "hann":
        tau     = (t - t[0]) / (t[-1] - t[0])
        v_perp  = v_perp * (0.5 - 0.5 * np.cos(2 * np.pi * tau))[:, None]

    if omega_array is None:
        is_uni, dt = _check_uniform(t)
        if not is_uni:
            dt = float(np.mean(np.diff(t)))
        f_max = 0.5 / dt
        omega_array = np.linspace(0, 2 * np.pi * f_max, n_omega + 1)[1:]

    freqs  = omega_array / (2 * np.pi)
    power  = np.zeros(len(omega_array))

    for k, omega in enumerate(omega_array):
        phase      = omega * (t - r_dot_n / C)
        integrand  = v_perp * np.exp(1j * phase)[:, None]
        integral   = _trapezoid(integrand, t, axis=0)
        power[k]   = omega**2 * float(np.real(np.dot(integral, np.conj(integral))))

    return freqs, power


def ensemble_spectrum(history: TrajectoryHistory,
                       method: str = "fft",
                       observer: np.ndarray = None,
                       weights: np.ndarray = None,
                       **kwargs) -> tuple:
    """Weighted incoherent sum of the single-particle spectra of an ensemble.

    Parameters
    ----------
    history : TrajectoryHistory
    method : "fft" or "retarded"
        Uses :func:`spectrum_fft` or :func:`spectrum_retarded`; further
        keyword arguments are passed on.
    observer : array_like, shape (3,) or None
    weights : ndarray, shape (N,) or None
        Particle weights, by default 1/N each.

    Returns
    -------
    freqs : ndarray
        Frequencies [Hz] of the first particle's spectrum; the others are
        interpolated onto them.
    total : ndarray
        Weighted sum.
    individual : list of ndarray
        Unweighted spectrum of each particle.
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
            p = np.interp(ref_freqs, f, p, left=0.0, right=0.0)
        total += weights[pid] * p
        individual.append(p)

    return ref_freqs, total, individual
