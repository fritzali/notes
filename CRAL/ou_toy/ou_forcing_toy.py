"""
Toy Ornstein-Uhlenbeck (OU) turbulence forcing field in Fourier space -- standalone NumPy.

NOT a Shamrock feature. A self-check of the forcing *statistics* before any implementation:
does the generated field have the decorrelation time T we asked for, the expected variance,
and the requested solenoidal/compressive split?

Model (in the spirit of Eswaran & Pope 1988; Federrath et al. 2010, A&A 512, A81):
  f(x, t) = 2 Re sum_k  fhat_k(t) exp(2 pi i k.x / L),   k integer vectors, 1 <= |k| <= 3,
  d fhat_k = -fhat_k dt / T + sigma A(|k|) sqrt(2/T) P_zeta(k) dW_k,
  A(|k|)   = 1 - (|k| - 2)^2   (parabola peaking at |k| = 2, zero at 1 and 3),
  P_zeta   = zeta P_perp + (1 - zeta) P_par,  P_perp = I - kk/k^2 (solenoidal), P_par = kk/k^2.

Two time integrators for the OU coefficients:
  "exact" : fhat <- fhat e^{-dt/T} + sigma A sqrt(1 - e^{-2 dt/T}) P_zeta xi   (exact for ANY dt)
  "euler" : fhat <- fhat (1 - dt/T)  + sigma A sqrt(2 dt/T)          P_zeta xi   (Euler-Maruyama)
Euler-Maruyama has a known, dt-dependent bias:
  effective decorrelation time  T_eff = -dt / ln(1 - dt/T)
  stationary variance ratio     var_EM / var_exact = 1 / (1 - dt/(2T)).
The script measures both from the generated data and compares with these formulas.

Checks:
  1. autocorrelation of the Fourier coefficients vs exp(-tau/T); fitted T and integral time;
  2. the same for the real-space field at a few random points;
  3. stationary variance relative to theory;
  4. divergence of the field for zeta = 1 (should be round-off) and curl for zeta = 0.

usage: python ou_forcing_toy.py [--T 1.0] [--zeta 1.0] [--seed 42] [--out ou_forcing_check.png]
"""

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# --------------------------------------------------------------------------------------
# Forcing modes and projection
# --------------------------------------------------------------------------------------
def make_modes(kmin=1.0, kmax=3.0):
    """Integer wave vectors (units of 2 pi / L) in a half space, 0 < A(|k|)."""
    r = int(np.ceil(kmax))
    ks = []
    for kx in range(-r, r + 1):
        for ky in range(-r, r + 1):
            for kz in range(-r, r + 1):
                # half space: the -k partner is the complex conjugate (real field)
                if (kz > 0) or (kz == 0 and ky > 0) or (kz == 0 and ky == 0 and kx > 0):
                    kn = np.sqrt(kx * kx + ky * ky + kz * kz)
                    if kmin <= kn <= kmax:
                        ks.append((kx, ky, kz))
    ks = np.array(ks, dtype=float)
    kn = np.linalg.norm(ks, axis=1)
    amp = np.clip(1.0 - (kn - 2.0) ** 2, 0.0, None)
    keep = amp > 0
    return ks[keep], amp[keep]


def projection_tensors(ks, zeta):
    """P_zeta(k) = zeta P_perp + (1 - zeta) P_par, shape (nmodes, 3, 3)."""
    khat = ks / np.linalg.norm(ks, axis=1, keepdims=True)
    P_par = np.einsum("mi,mj->mij", khat, khat)
    P_perp = np.eye(3)[None, :, :] - P_par
    return zeta * P_perp + (1.0 - zeta) * P_par


# --------------------------------------------------------------------------------------
# OU integration of the Fourier coefficients
# --------------------------------------------------------------------------------------
def run_ou(ks, amp, P, T, sigma, dt, nsteps, scheme, rng, nburn):
    nm = len(ks)

    def noise():
        xi = (rng.standard_normal((nm, 3)) + 1j * rng.standard_normal((nm, 3))) / np.sqrt(2.0)
        return amp[:, None] * np.einsum("mij,mj->mi", P, xi)

    if scheme == "exact":
        decay = np.exp(-dt / T)
        kick = sigma * np.sqrt(1.0 - np.exp(-2.0 * dt / T))
    elif scheme == "euler":
        decay = 1.0 - dt / T
        kick = sigma * np.sqrt(2.0 * dt / T)
    else:
        raise ValueError(scheme)

    # start from the stationary distribution of the exact process
    f = sigma * noise()
    for _ in range(nburn):
        f = decay * f + kick * noise()
    hist = np.empty((nsteps, nm, 3), dtype=complex)
    for n in range(nsteps):
        f = decay * f + kick * noise()
        hist[n] = f
    return hist


def field_at(points, ks, fhat):
    """Real-space forcing at points (npts,3) in a unit box, fhat (..., nmodes, 3)."""
    phase = np.exp(2j * np.pi * points @ ks.T)  # (npts, nmodes)
    return 2.0 * np.real(np.einsum("pm,...mi->...pi", phase, fhat))


# --------------------------------------------------------------------------------------
# Statistics
# --------------------------------------------------------------------------------------
def acf(series):
    """Normalised autocorrelation along axis 0, averaged over all other axes (real input)."""
    x = series - series.mean(axis=0)
    n = x.shape[0]
    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    X = np.fft.rfft(x, n=nfft, axis=0)
    c = np.fft.irfft(X * np.conj(X), n=nfft, axis=0)[:n]
    c /= (n - np.arange(n)).reshape((-1,) + (1,) * (x.ndim - 1))  # unbiased
    c = c.reshape(n, -1).mean(axis=1)
    return c / c[0]


def fit_T(tau, rho, T_guess):
    """Fit rho = exp(-tau/T) on tau <= 1.5 T_guess (log-linear, through the origin)."""
    m = (tau > 0) & (tau <= 1.5 * T_guess) & (rho > 0.05)
    slope = np.sum(tau[m] * np.log(rho[m])) / np.sum(tau[m] ** 2)
    return -1.0 / slope


def integral_time(tau, rho):
    """Integral of the ACF up to its first zero crossing."""
    cut = np.argmax(rho <= 0) if np.any(rho <= 0) else len(rho)
    return np.trapezoid(rho[:cut], tau[:cut])


# --------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--T", type=float, default=1.0, help="decorrelation (autocorrelation) time")
    ap.add_argument("--sigma", type=float, default=1.0, help="rms amplitude per mode component")
    ap.add_argument("--zeta", type=float, default=1.0, help="1 = solenoidal, 0 = compressive")
    ap.add_argument("--duration", type=float, default=400.0, help="sampled duration in units of T")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default="ou_forcing_check.png")
    args = ap.parse_args()

    T, sigma, zeta = args.T, args.sigma, args.zeta
    rng = np.random.default_rng(args.seed)  # seeded: the run is reproducible

    ks, amp = make_modes()
    P = projection_tensors(ks, zeta)
    print(f"forcing modes (half space, 1<|k|<3): {len(ks)}   T = {T}   zeta = {zeta}   seed = {args.seed}")

    pts = rng.random((8, 3))  # random sample points in the unit box

    cases = [("exact", 0.05), ("exact", 0.5), ("euler", 0.05), ("euler", 0.5)]
    results = {}
    for scheme, dt_over_T in cases:
        dt = dt_over_T * T
        nsteps = int(round(args.duration / dt_over_T))
        hist = run_ou(ks, amp, P, T, sigma, dt, nsteps, scheme, rng, nburn=int(10 / dt_over_T))
        tau = dt * np.arange(nsteps)

        # (1) Fourier coefficients: real and imaginary parts of every non-zero component
        comp = np.concatenate([hist.real, hist.imag], axis=2).reshape(nsteps, -1)
        comp = comp[:, np.abs(comp).mean(axis=0) > 1e-12 * sigma]
        rho_k = acf(comp)

        # (2) real-space field at sample points
        fx = field_at(pts, ks, hist)  # (nsteps, npts, 3)
        rho_x = acf(fx.reshape(nsteps, -1))

        # (3) variance of the coefficients vs the exact stationary value
        var_theory = sigma**2 * np.mean(
            amp[:, None] ** 2 * np.einsum("mij,mij->mi", P, P).real
        ) / 2.0  # per real component: E|P xi|^2 / 2 per Re/Im part
        var_meas = np.mean(np.concatenate([hist.real, hist.imag], axis=2) ** 2)

        T_fit_k = fit_T(tau, rho_k, T)
        T_fit_x = fit_T(tau, rho_x, T)
        T_int_k = integral_time(tau, rho_k)
        if scheme == "exact":
            T_pred, var_ratio_pred = T, 1.0
        else:
            T_pred = -dt / np.log(1.0 - dt_over_T)
            var_ratio_pred = 1.0 / (1.0 - dt_over_T / 2.0)

        results[(scheme, dt_over_T)] = dict(
            tau=tau, rho_k=rho_k, rho_x=rho_x, T_fit_k=T_fit_k, T_fit_x=T_fit_x,
            T_int_k=T_int_k, T_pred=T_pred, var_ratio=var_meas / var_theory,
            var_ratio_pred=var_ratio_pred,
        )
        print(
            f"{scheme:5s} dt/T={dt_over_T:4.2f}: T_fit(coeffs)={T_fit_k:.3f}  T_fit(field)={T_fit_x:.3f}  "
            f"T_integral={T_int_k:.3f}  [predicted {T_pred:.3f}]   "
            f"var/var_exact={var_meas / var_theory:.3f}  [predicted {var_ratio_pred:.3f}]"
        )

    # (4) divergence / curl of the field, spectrally: k . fhat and k x fhat
    hist = run_ou(ks, amp, P, T, sigma, 0.1 * T, 50, "exact", rng, nburn=0)
    fhat = hist[-1]
    knorm = np.linalg.norm(ks, axis=1)
    rel_div = np.abs(np.einsum("mi,mi->m", ks, fhat)) / (knorm * np.linalg.norm(fhat, axis=1))
    rel_curl = np.linalg.norm(np.cross(ks, fhat), axis=1) / (knorm * np.linalg.norm(fhat, axis=1))
    print(f"max |k.fhat|/(|k||fhat|) = {rel_div.max():.2e}   max |k x fhat|/(|k||fhat|) = {rel_curl.max():.2e}"
          f"   (zeta=1 -> first ~0; zeta=0 -> second ~0)")

    # ----------------------------------------------------------------------------------
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.5), dpi=130)
    taus = np.linspace(0, 3 * T, 300)
    for ax, key, title in [(axs[0], "rho_k", "Fourier coefficients"), (axs[1], "rho_x", "real-space field at 8 points")]:
        ax.plot(taus / T, np.exp(-taus / T), "k-", lw=2, label=r"target $e^{-\tau/T}$")
        for (scheme, r), res in results.items():
            m = res["tau"] <= 3 * T
            ls = "o" if scheme == "exact" else "s"
            ax.plot(res["tau"][m] / T, res[key][m], ls, ms=3, alpha=0.7,
                    label=f"{scheme}, dt/T={r}: T_fit={res['T_fit_k' if key == 'rho_k' else 'T_fit_x'] / T:.3f} T")
        ax.set_xlabel(r"lag $\tau / T$")
        ax.set_ylabel("normalised autocorrelation")
        ax.set_title(title)
        ax.set_ylim(-0.1, 1.05)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)

    # slice of the field at z = 0.5 (last snapshot), arrows = in-plane components
    n = 32
    g = (np.arange(n) + 0.5) / n
    X, Y = np.meshgrid(g, g, indexing="xy")
    pts_slice = np.column_stack([X.ravel(), Y.ravel(), 0.5 * np.ones(X.size)])
    fs = field_at(pts_slice, ks, fhat).reshape(n, n, 3)
    im = axs[2].imshow(np.linalg.norm(fs, axis=2), origin="lower", extent=[0, 1, 0, 1], cmap="viridis")
    s = 2
    axs[2].quiver(X[::s, ::s], Y[::s, ::s], fs[::s, ::s, 0], fs[::s, ::s, 1], color="w", scale=None)
    axs[2].set_title(f"|f| at z=0.5, zeta={zeta} (one snapshot)")
    fig.colorbar(im, ax=axs[2])
    fig.suptitle(f"Toy OU forcing: {len(ks)} modes, target T={T}, seed={args.seed}")
    fig.tight_layout()
    fig.savefig(args.out)
    print("saved", args.out)


if __name__ == "__main__":
    main()
