"""§12.12 (C16): Taylor's theory of turbulent dispersion — Langevin particles from a point source (Fig. 12.24 analogue), the
mean-square displacement against Taylor's formula (12.119) in closed form, the ballistic (12.121) and diffusive (12.123)
limits, the identity d<X^2>/dt = 2<Xu> (12.116), the eddy diffusivity that grows before it saturates (12.127)-(12.129) with
the printed condition of (12.129) beside the corrected one, the random walk (12.125) (Fig. 12.26 analogue) and the smoke
plume (Fig. 12.27 analogue: linear near the source, square root far away).

Run: ``.venv/Scripts/python.exe scripts/ch12_dispersion.py --no-show [--fast]``
Figures -> outputs/ch12/c16_dispersion.png, c16_random_walk_plume.png.
Our numbers: u_rms = 0.5 m/s, Lambda_t = 2 s, 10000 particles (fast 2000), seed 0; plume in a 5 m/s wind, w_rms = 0.5 m/s, Lambda_t = 10 s.
"""
from __future__ import annotations

import numpy as np

from ch12_common import COLORS, Timer, finish, parse_args, save, setup

from fluidpy import ch12_turbulence as ch12
from fluidpy.core import turbstats as TS


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    urms, Lam = 0.5, 2.0
    u2 = urms ** 2
    npart = 2000 if args.fast else 10000
    t = np.concatenate([[0.0], np.geomspace(0.02, 60.0, 400 if args.fast else 1000)])
    with Timer(f"{npart} Langevin particles x {t.size} steps"):
        Xp, up = ch12.langevin_particles(npart, t, urms, Lam, seed=0)
    X2 = np.mean(Xp ** 2, axis=0)
    exact = ch12.taylor_dispersion_exponential(t, u2, Lam)
    se = np.sqrt(2.0 / npart) * exact
    worst = np.max(np.abs(X2 - exact)[1:] / se[1:])
    print(f"<X^2> of the particles against Taylor's closed form: worst deviation {worst:.2f} standard errors over {t.size - 1} times")
    r = lambda tau: np.exp(-abs(tau) / Lam)  # noqa: E731
    for tt in (0.2, 2.0, 20.0):
        a, b = ch12.taylor_dispersion(tt, r, u2), ch12.taylor_dispersion(tt, r, u2, form="double")
        print(f"  t = {tt:5.1f} s ({ch12.dispersion_regime(tt, Lam):10s}): (12.119) {a:.6f}, (12.118) {b:.6f}, closed form "
              f"{ch12.taylor_dispersion_exponential(tt, u2, Lam):.6f} m2; ballistic u2 t^2 = {u2 * tt ** 2:.4f}, diffusive 2 u2 Lambda t = "
              f"{2 * u2 * Lam * tt:.4f}; local slope {ch12.dispersion_local_slope(tt, Lam):.3f}")
    lhs, rhs = ch12.dispersion_rate_from_particles(t, Xp, up)
    mid = slice(5, -5)
    print(f"d<X^2>/dt vs 2<Xu> (12.115)-(12.116): max difference / max value = {np.max(np.abs(lhs[mid] - rhs[mid]) / np.max(rhs)):.1e}"
          f" (differencing a sampled mean of {npart} rough paths)")
    lagU, RU = TS.autocorrelation(ch12.langevin_particles(1, np.arange(200000) * 0.02, urms, Lam, seed=3)[1][0], 0.02, max_lag=600)
    print(f"Lagrangian integral scale of one long particle record: {TS.integral_scale(lagU, RU):.3f} s (prescribed {Lam})")
    print("eddy diffusivity D_T (12.127), exact / short-time (12.128) / long-time (12.129) / (12.129) under the PRINTED condition:")
    for tt in (0.02, 0.2, 2.0, 20.0, 200.0):
        print(f"  t = {tt:6.2f} s: {ch12.eddy_diffusivity_exponential(tt, u2, Lam):.5f}  {ch12.eddy_diffusivity_asymptote(tt, u2, Lam, 'short'):.5f}  "
              f"{ch12.eddy_diffusivity_asymptote(tt, u2, Lam, 'long'):.5f}  {ch12.eddy_diffusivity_asymptote(tt, u2, Lam, 'long', printed=True):.5f}")
    tg = 1.5
    rg = lambda tau: np.exp(-(tau / tg) ** 2)  # noqa: E731
    print(f"Gaussian correlation, t_c = {tg} s: Lambda_t = sqrt(pi) t_c/2 = {np.sqrt(np.pi) * tg / 2:.4f} s; <X^2>(5 s) = "
          f"{ch12.taylor_dispersion(5.0, rg, u2):.6f} (quad) = {ch12.taylor_dispersion_gaussian(5.0, u2, tg):.6f} (closed form); "
          f"D_T(inf) = {ch12.eddy_diffusivity_taylor(50.0, rg, u2):.5f} m2/s")
    for dim in (1, 2, 3):
        R2 = np.mean(np.sum(ch12.random_walk(400, 4000 if args.fast else 20000, L=1.0, dim=dim, seed=dim) ** 2, axis=2), axis=0)
        print(f"random walk, {dim}-D: <R_n^2>/n at n = 400 -> {R2[-1] / 400:.4f} (expected 1); rms = {np.sqrt(R2[-1]):.2f} = L sqrt(n) = 20")
    R2p = np.mean(np.sum(ch12.random_walk(2000, 4000, dim=2, persistence=0.5, seed=1) ** 2, axis=2), axis=0)
    print(f"  with persistence 0.5: <R_n^2>/n -> {R2p[-1] / 2000:.3f} (expected (1 + p)/(1 - p) = 3)")
    U, wrms, LamP = 5.0, 0.5, 10.0
    x = np.geomspace(1.0, 1e5, 60)
    Z = ch12.smoke_plume_width(x, U, wrms, LamP)
    sl = np.gradient(np.log(Z), np.log(x))
    print(f"smoke plume, U Lambda_t = {U * LamP:.0f} m: Z_rms(10 m) = {ch12.smoke_plume_width(10.0, U, wrms, LamP):.3f} m, "
          f"Z_rms(10 km) = {ch12.smoke_plume_width(1e4, U, wrms, LamP):.1f} m; local slope d ln Z/d ln x = {sl[0]:.3f} near, {sl[-1]:.3f} far")
    zz = np.linspace(-400.0, 400.0, 4001)
    print(f"  plume concentration integrates to Q/U at x = 1 km: {np.trapezoid(ch12.plume_concentration(1e3, zz, 2.0, U, ch12.smoke_plume_width(1e3, U, wrms, LamP)), zz):.5f}"
          f" (Q/U = {2.0 / U:.5f})")
    print(f"  constant-D estimate from the variance of a Gaussian patch sigma^2 = 2 D t, D = 0.3: {ch12.diffusivity_from_variance(t[1:], 0.6 * t[1:])[10]:.4f} m2/s")

    fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.8))
    tl = np.linspace(0.0, 30.0, 601)
    X12, _ = ch12.langevin_particles(12, tl, urms, Lam, seed=5, dim=2)
    for k in range(12):
        ax[0].plot(X12[k, :, 0], X12[k, :, 1], lw=0.9)
    ax[0].plot(0, 0, "ko")
    ax[0].set(aspect="equal", xlabel=r"$X_1$ [m]", ylabel=r"$X_2$ [m]", title="12 particle paths from a point source (30 s)")
    ax[1].loglog(t[1:], np.sqrt(X2[1:]), ".", ms=3, color=COLORS["blue"], label="particles")
    ax[1].loglog(t[1:], np.sqrt(exact[1:]), color=COLORS["accent"], label="Taylor (12.119), exponential r")
    ax[1].loglog(t[1:], urms * t[1:], "k:", label=r"$u_{rms}t$ (12.121)")
    ax[1].loglog(t[1:], urms * np.sqrt(2 * Lam * t[1:]), "k--", lw=1.0, label=r"$u_{rms}\sqrt{2\Lambda_t t}$ (12.123)")
    ax[1].axvline(Lam, color=COLORS["grid"])
    ax[1].set(xlabel="t [s]", ylabel=r"$X_{rms}$ [m]", ylim=(5e-3, 30), title="t first, then the square root of t")
    ax[1].legend(fontsize=8)
    ax[2].semilogx(t[1:], ch12.eddy_diffusivity_exponential(t[1:], u2, Lam), color=COLORS["amber"], label=r"$D_T(t)$ (12.127)")
    ax[2].semilogx(t[1:], 0.5 * rhs[1:], ".", ms=2, color=COLORS["blue"], label=r"$\overline{Xu}$ of the particles")
    ax[2].axhline(u2 * Lam, color=COLORS["muted"], ls="--", label=r"$\overline{u^2}\Lambda_t$ (12.129), $t\gg\Lambda_t$")
    ax[2].semilogx(t[1:40], u2 * t[1:40], "k:", label=r"$\overline{u^2}t$ (12.128)")
    ax[2].set(xlabel="t [s]", ylabel=r"$D_T$ [m$^2$/s]", ylim=(0, 0.7), title="the eddy diffusivity is not a constant")
    ax[2].legend(fontsize=8)
    save(fig, out, "c16_dispersion")

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    w = ch12.random_walk(200, 2000, dim=2, seed=7)
    for k in range(3):
        ax[0].plot(w[k, :, 0], w[k, :, 1], lw=0.8)
    th = np.linspace(0, 2 * np.pi, 200)
    ax[0].plot(np.sqrt(200) * np.cos(th), np.sqrt(200) * np.sin(th), "k--", lw=1.0, label=r"$L\sqrt{n}$, n = 200")
    ax[0].set(aspect="equal", xlabel="x / L", ylabel="y / L", title="three random walks (12.125)")
    ax[0].legend(fontsize=8)
    xs = np.linspace(0.0, 600.0, 400)
    Zs = ch12.smoke_plume_width(xs, U, wrms, LamP)
    ax[1].fill_between(xs, -Zs, Zs, color=COLORS["muted"], alpha=0.35, label=r"$\pm Z_{rms}$")
    ax[1].plot(xs, wrms / U * xs, "k:", label=r"near: $Z_{rms}\propto x$")
    ax[1].plot(xs, wrms * np.sqrt(2 * LamP * xs / U), "k--", lw=1.0, label=r"far: $Z_{rms}\propto\sqrt{x}$")
    ax[1].set(ylim=(-80, 80), xlabel="x [m] downwind", ylabel="z [m]", title="time-averaged smoke plume: a parabola with a pointed vertex")
    ax[1].legend(fontsize=8)
    save(fig, out, "c16_random_walk_plume")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
