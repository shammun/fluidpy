"""§12.7 (C06-C08): the energy cascade and the spectrum — Kolmogorov scales by dimensional analysis (12.50), scale separation
(12.51)-(12.52), the tiers of the cascade, our scale table, the -5/3 law (12.54) with the printed +5/3 beside it, the
relation C1 = (18/55) C, and model spectra in Kolmogorov scaling (Fig. 12.12 analogue).

Run: ``.venv/Scripts/python.exe scripts/ch12_cascade_spectrum.py --no-show [--fast]``
Figures -> outputs/ch12/c07_cascade_scales.png, c08_spectrum_kolmogorov_scaling.png.
"""
from __future__ import annotations

import numpy as np

from ch12_common import COLORS, Timer, finish, parse_args, save, setup

from fluidpy import ch12_turbulence as ch12


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    print("exponents of (12.50) from the Pi theorem (nu^a eps^b):", {k: tuple(str(x) for x in v) for k, v in ch12.kolmogorov_exponents().items()})
    print("our cases (not the book's):")
    for c in ch12.scale_table():
        print(f"  {c['name']:27s} Re_L {c['Re_L']:.2e}  eps {c['eps']:.2e} m2/s3  eta {c['eta']:.2e} m  tau_eta {c['tau_eta']:.2e} s  "
              f"lambda_g {c['lambda_g']:.2e} m  R_lambda {c['R_lambda']:.0f}  decades {c['decades']:.2f}  DNS points {c['grid_points']:.1e}")
    for Re in (1e3, 1e5, 1e7):
        s = ch12.scale_separation(Re)
        print(f"  Re_L = {Re:.0e}: eta/L = {s['eta_over_L']:.2e}, lambda_f/L = {s['lambdaT_over_L']:.2e}, u_K/dU = {s['uK_over_dU']:.3f}, "
              f"decades {ch12.inertial_range_decades(Re):.2f}, ordered {ch12.scale_ordering(Re)['ordered']}")
    tiers = ch12.cascade_tiers(L=1.0, dU=1.0, nu=1.5e-5)
    print(f"cascade from L = 1 m, dU = 1 m/s in air: {tiers['n_tiers']} tiers of ratio 2, eddy Reynolds number from "
          f"{tiers['Re'][0]:.0f} down to {tiers['Re'][-1]:.1f}, eta = {tiers['eta'] * 1e3:.3f} mm")
    kc = ch12.kolmogorov_constants()
    print(f"constants: C = {kc['C']}, one-sided C1 = {kc['C1_one_sided']:.4f}, two-sided {kc['C1_two_sided']:.4f}")
    eps, nu = 0.1, 1.5e-5
    ratio = ch12.one_dimensional_from_3d(50.0, lambda K: ch12.inertial_spectrum_3d(K, eps)) / ch12.inertial_spectrum_1d(50.0, eps, two_sided=False)
    print(f"one-dimensional spectrum of a pure K^-5/3 law / [(18/55) C eps^2/3 k1^-5/3] = {ratio:.10f}")
    k = np.array([10.0, 100.0])
    good, bad = ch12.inertial_spectrum_1d(k, eps), ch12.inertial_spectrum_1d(k, eps, printed=True)
    print(f"slope of (12.54): corrected {np.log(good[1] / good[0]) / np.log(10):+.4f}, as printed {np.log(bad[1] / bad[0]) / np.log(10):+.4f}")
    Kq = np.geomspace(1e-4, 1e7, 200001)
    for L in (None, 1.0):
        Em = ch12.model_spectrum(Kq, eps, nu, L=L)
        print(f"model spectrum (L = {L}): 2 nu int K^2 E dK / eps = {2 * nu * np.trapezoid(Kq ** 2 * Em, Kq) / eps:.6f}"
              + ("" if L is None else f", int E dK = {np.trapezoid(Em, Kq):.4f} m2/s2"))
    print(f"Richardson's law: K(l = 10 m) = {ch12.richardson_diffusivity(10.0, 1e-4):.3f} m2/s, K(100 m) = "
          f"{ch12.richardson_diffusivity(100.0, 1e-4):.3f} m2/s (ratio 10^(4/3) = {10 ** (4 / 3):.2f})")

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    Re = np.geomspace(1e2, 1e9, 50)
    s = ch12.scale_separation(Re)
    ax[0].loglog(Re, s["eta_over_L"], color=COLORS["accent"], label=r"$\eta/L \sim Re_L^{-3/4}$ (12.51)")
    ax[0].loglog(Re, s["lambdaT_over_L"], color=COLORS["teal"], label=r"$\lambda_T/L \propto Re_L^{-1/2}$ (12.52)")
    ax[0].loglog(Re, np.ones_like(Re), color=COLORS["muted"], label="L")
    ax[0].set(xlabel=r"$Re_L=\Delta U L/\nu$", ylabel="scale / L", title="scale separation grows with the Reynolds number")
    ax[0].legend(fontsize=8)
    ax[1].loglog(tiers["size"], tiers["Re"], "o-", color=COLORS["orange"])
    ax[1].axhline(1.0, color=COLORS["muted"], ls=":")
    ax[1].axvline(tiers["eta"], color=COLORS["rose"], ls="--", label=r"$\eta$")
    ax[1].set(xlabel="eddy size l' [m]", ylabel=r"eddy Reynolds number $u'l'/\nu$", title="tiers of the cascade (L = 1 m, air)")
    ax[1].legend(fontsize=8)
    save(fig, out, "c07_cascade_scales")

    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    n_k = 40 if args.fast else 80
    with Timer("one-dimensional spectra from E(K)"):
        for Lm, col in ((0.01, COLORS["muted"]), (0.1, COLORS["teal"]), (1.0, COLORS["orange"]), (10.0, COLORS["accent"])):
            eta = ch12.kolmogorov_scales(nu, eps)[0]
            k1 = np.geomspace(0.3 / Lm, 3.0 / eta, n_k)
            E11 = ch12.one_dimensional_from_3d(k1, lambda K: ch12.model_spectrum(K, eps, nu, L=Lm))
            x, phi = ch12.kolmogorov_normalize_spectrum(k1, E11, nu, eps)
            ax.loglog(x, phi, color=col, label=rf"$L/\eta$ = {Lm / eta:.0f}")
    xx = np.geomspace(1e-5, 0.3, 50)
    ax.loglog(xx, kc["C1_one_sided"] * xx ** (-5 / 3), "k--", lw=1.0, label=r"$C_1 (k_1\eta)^{-5/3}$")
    ax.loglog(xx, 1e3 * kc["C1_one_sided"] * xx ** (5 / 3), color=COLORS["rose"], ls=":", lw=1.0, label="+5/3 as printed (shifted)")
    ax.set(xlabel=r"$k_1\eta$", ylabel=r"$2S_{11}/(u_K^2\eta)$", ylim=(1e-4, 1e9),
           title="model spectra in Kolmogorov scaling (12.53): one curve at high wavenumber")
    ax.legend(fontsize=8)
    save(fig, out, "c08_spectrum_kolmogorov_scaling")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
