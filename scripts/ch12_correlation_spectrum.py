"""§12.4 (C02, C03): correlations and spectra — the delayed-copy cross-correlation (Fig. 12.4 analogue), r(tau) with its
integral scale, first zero and osculating parabola (Fig. 12.5 analogue), and one spectrum obtained three ways (transform of
the exact correlation, transform of the measured correlation, periodogram), with Parseval's check and S(0) = variance·Lambda/pi.

Run: ``.venv/Scripts/python.exe scripts/ch12_correlation_spectrum.py --no-show [--fast]``
Figures -> outputs/ch12/c02_correlations.png, c03_spectrum_three_ways.png.
Our signal: Gaussian-correlated record, sigma = 1.2 m/s, tau_c = 0.25 s (smooth, so the Taylor microscale exists), seed 0.
"""
from __future__ import annotations

import numpy as np

from ch12_common import COLORS, finish, parse_args, save, setup

from fluidpy.core import turbstats as TS


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    dt, n = 0.01, (2 ** 15 if args.fast else 2 ** 18)
    pair = TS.correlation_spectrum_pair("gaussian", 1.2, 0.25)
    u = TS.smooth_signal(n, dt, pair["S"], seed=0)
    t = np.arange(n) * dt
    lag, R = TS.autocorrelation(u, dt, max_lag=400)
    lag_d, R_d = TS.autocorrelation(u, dt, max_lag=400, method="direct")
    r = R / R[0]
    Lam, tc, lam = TS.integral_scale(lag, r), TS.correlation_time(lag, r), TS.taylor_microscale(lag, r)
    print(f"record: {n} samples, variance {u.var():.4f} (prescribed {pair['variance']:.4f})")
    print(f"FFT vs direct autocorrelation: max difference {np.max(np.abs(R - R_d)):.1e}")
    print(f"integral scale  Lambda_t = {Lam:.4f} s   (exact {pair['Lambda_t']:.4f})")
    print(f"Taylor microscale lambda_t = {lam:.4f} s (exact {pair['lambda_t']:.4f})")
    print(f"first zero of the measured r: t_c = {tc:.3f} s (the exact Gaussian r has none: this is where sampling noise first crosses zero)"
          f" -> about {TS.effective_samples(n * dt, tc):.0f} independent samples in {n * dt:.0f} s by the book's N = T/t_c")

    delay = 0.6
    v = np.roll(u, int(round(delay / dt)))  # v(t) = u(t - delay)
    lag2, Ruv = TS.cross_correlation(u, v, dt, max_lag=200)
    print(f"cross-correlation of u with a copy delayed by {delay} s peaks at lag {lag2[np.argmax(Ruv)]:+.2f} s")
    lag3, Rvu = TS.cross_correlation(v, u, dt, max_lag=200)
    print(f"R_uv(tau) = R_vu(-tau): max difference {np.max(np.abs(Ruv - Rvu[::-1])):.1e}")

    omega = np.linspace(0.0, 40.0, 401)
    S_exact = pair["S"](omega)
    S_corr = TS.spectrum_from_correlation(lag, R, omega)
    om_p, S_p, w_p = TS.periodogram(u, dt, segments=16 if args.fast else 64, return_weights=True)
    om_1, S_1, w_1 = TS.periodogram(u, dt, return_weights=True)
    print(f"Parseval, one boxcar segment: sum(w S) - variance = {TS.spectrum_variance(om_1, S_1, weights=w_1) - u.var():.1e}")
    print(f"area under the exact two-sided spectrum: {TS.spectrum_variance(np.linspace(0, 80, 8001), pair['S'](np.linspace(0, 80, 8001))):.6f}"
          f" (variance {pair['variance']:.6f})")
    print(f"S(0): exact {pair['S0']:.4f}, from the measured correlation {S_corr[0]:.4f}, variance*Lambda/pi = {u.var() * Lam / np.pi:.4f}"
          f"  -> Lambda from S(0): {TS.integral_scale_from_spectrum(pair['S0'], pair['variance']):.4f} s")
    k1, S11 = TS.frequency_to_wavenumber_spectrum(omega, S_exact, U0=5.0)
    print(f"frozen-turbulence conversion at U0 = 5 m/s: variance in k-space {TS.spectrum_variance(k1, S11):.5f}, in omega-space "
          f"{TS.spectrum_variance(omega, S_exact):.5f}")
    for kind in ("exponential", "damped_cosine"):
        p = TS.correlation_spectrum_pair(kind, 1.0, 0.5, omega0=6.0)
        tau = np.linspace(0, 30, 60001)
        print(f"  {kind:14s}: Lambda {TS.integral_scale(tau, p['r'](tau), 'all'):.5f} (exact {p['Lambda_t']:.5f}), "
              f"S(0) {float(TS.spectrum_from_correlation(tau, p['R'](tau), 0.0)):.5f} (exact {p['S0']:.5f})")

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    m = t < 6.0
    ax[0].plot(t[m], u[m], color=COLORS["accent"], lw=1.0, label="u(t)")
    ax[0].plot(t[m], v[m] - 5.0, color=COLORS["teal"], lw=1.0, label=f"v(t) = u(t - {delay} s), offset")
    ax[0].set(xlabel="t [s]", ylabel="signal [m/s]", title="a record and a delayed copy")
    ax[0].legend(fontsize=8)
    ax[1].plot(lag2, Ruv, color=COLORS["teal"], label=r"$\overline{u(t)v(t+\tau)}$")
    ax[1].plot(np.concatenate([-lag[::-1], lag]), np.concatenate([R[::-1], R]), color=COLORS["accent"], label=r"$\overline{u(t)u(t+\tau)}$")
    ax[1].axvline(delay, color=COLORS["orange"], ls=":")
    ax[1].set(xlim=(-2, 2), xlabel=r"lag $\tau$ [s]", ylabel=r"correlation [m$^2$/s$^2$]", title="the peak sits at the aligning lag")
    ax[1].legend(fontsize=8)
    ax[2].plot(lag, r, color=COLORS["accent"], label=r"$r(\tau)$ measured")
    ax[2].plot(lag, pair["r"](lag), "k--", lw=1.0, label="exact")
    ax[2].plot(lag, 1 - (lag / lam) ** 2, color=COLORS["rose"], lw=1.0, label=r"osculating parabola, $\lambda_t$")
    ax[2].plot([0, Lam, Lam], [1, 1, 0], color=COLORS["orange"], ls="--", label=r"equal-area rectangle, $\Lambda_t$")
    ax[2].set(xlim=(0, 1.2), ylim=(-0.2, 1.1), xlabel=r"$\tau$ [s]", ylabel=r"$r(\tau)$", title="memory time and microscale (12.18)-(12.19)")
    ax[2].legend(fontsize=8)
    save(fig, out, "c02_correlations")

    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    ax.loglog(om_p[1:], S_p[1:], color=COLORS["grid"], lw=1.0, label="periodogram (segment average)")
    ax.loglog(omega[1:], S_corr[1:], color=COLORS["teal"], label="transform of the measured correlation (12.20)")
    ax.loglog(omega[1:], S_exact[1:], "k--", lw=1.2, label="exact transform")
    ax.axhline(pair["S0"], color=COLORS["orange"], ls=":", label=r"$S(0)=\overline{u^2}\Lambda_t/\pi$")
    ax.set(xlim=(0.2, 40), ylim=(1e-8, 1), xlabel=r"$\omega$ [rad/s]", ylabel=r"$S_e(\omega)$ [m$^2$/s]", title="one spectrum, three routes")
    ax.legend(fontsize=8)
    save(fig, out, "c03_spectrum_three_ways")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
