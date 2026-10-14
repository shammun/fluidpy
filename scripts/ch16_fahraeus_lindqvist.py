"""Chapter 16, §16.2–§16.3 — the Fåhræus–Lindqvist effect: hematocrit bookkeeping (16.8), the two-layer (core +
cell-free plasma layer) Poiseuille model (16.21)–(16.34) with the corrected core flux (16.30) (slip S3), Charm and
Kurland's core viscosity (16.35)–(16.36), our closure of the model, and the cited in-vitro law of Pries et al. (1992).

Run: ``.venv/Scripts/python.exe scripts/ch16_fahraeus_lindqvist.py --no-show``   Figures -> outputs/ch16/fl_*.png
"""
from __future__ import annotations

import numpy as np
from ch16_common import C_E, C_P, C_REF, C_TAU, C_U, C_WALL, finish, parse_args, save, setup
from scipy.integrate import quad

from fluidpy import ch16_biofluid_mechanics as ch16


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch16.illustrative_inputs()
    mu_p, T = I["mu_plasma"], I["T_body"]
    a, delta, dp, L = 40e-6, I["delta_cfl"], 200.0, 1e-3
    mu_c = ch16.charm_kurland_core_viscosity(0.5, T, mu_p)
    q = ch16.two_layer_flow_rates(a, delta, dp, L, mu_p, mu_c)
    qp = ch16.two_layer_flow_rates(a, delta, dp, L, mu_p, mu_c, printed_16_30=True)
    Qq = quad(lambda r: 2 * np.pi * r * ch16.two_layer_velocity(r, a, delta, dp, L, mu_p, mu_c), 0, a,
              points=[a - delta], epsrel=1e-13)[0]
    print(f"two-layer model (a = {a * 1e6:.0f} µm, δ = {delta * 1e6:.1f} µm, μ_c = {mu_c * 1e3:.3f} mPa s from (16.35)–(16.36) at H_c = 0.5):")
    print(f"  Q_p (16.29) = {q['Qp']:.5e}, Q_c (16.30, corrected) = {q['Qc']:.5e}, sum = {q['Q']:.5e} m³/s; (16.31) = "
          f"{ch16.two_layer_total_flow_rate(a, delta, dp, L, mu_p, mu_c):.5e}; quadrature of the profile = {Qq:.5e}")
    print(f"  printed (16.30) gives Q_c = {qp['Qc']:.5e}: Q_p + Q_c off (16.31) by {qp['Q'] / q['Q'] - 1:+.1%} (slip S3)")
    print(f"  μ_app (16.33) = {ch16.apparent_viscosity_two_layer(a, delta, mu_p, mu_c) * 1e3:.4f} mPa s; thin-layer (16.34) = "
          f"{ch16.apparent_viscosity_thin_layer(a, delta, mu_p, mu_c) * 1e3:.4f} mPa s; μ_c = {mu_c * 1e3:.4f} mPa s")
    cc = ch16.cell_conservation(q["Qc"], q["Qp"], 0.5, a, delta)
    print(f"  cell balance (16.8) with H_c = 0.5: H_F = H_D = {cc['H_F']:.4f}, H_T = {cc['H_T']:.4f} "
          f"(H_T/H_D = {cc['HT_over_HD']:.4f} < 1: the Fåhræus effect)")
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    r = np.linspace(0, a, 400)
    u = ch16.two_layer_velocity(r, a, delta, dp, L, mu_p, mu_c)
    uN = ch16.poiseuille_velocity(r, a, dp, L, mu_c)
    rr = np.concatenate([-r[::-1], r]) * 1e6
    ax[0].plot(rr, np.concatenate([u[::-1], u]) * 1e3, color=C_U, label="two-layer (16.27)–(16.28)")
    ax[0].plot(rr, np.concatenate([uN[::-1], uN]) * 1e3, color=C_REF, ls=":", label="all core fluid (μ_c)")
    for s in (-1, 1):
        ax[0].axvspan(s * (a - delta) * 1e6, s * a * 1e6, color=C_WALL, alpha=0.12)
    ax[0].set(xlabel="r [µm] (full diameter; shaded: cell-free layer)", ylabel="u [mm/s]",
              title="the plasma layer lubricates the core (kink at a − δ)")
    ax[0].legend()
    # μ_app against diameter: our closure (two δ) and the Pries law
    D = np.geomspace(8.0, 500.0, 60)
    ax[1].semilogx(D, ch16.pries_relative_viscosity(D, 0.45), color=C_TAU, label="Pries et al. 1992 (cited), H_D = 0.45")
    for dl, col in ((1.5e-6, C_E), (2.5e-6, C_P)):
        mur = []
        for d in D:
            try:
                mur.append(ch16.fahraeus_lindqvist_model(d / 2 * 1e-6, dl, 0.45, mu_p, T)["mu_rel"])
            except ValueError:
                mur.append(np.nan)
        ax[1].semilogx(D, mur, color=col, ls="--", label=f"two-layer closure (ours), δ = {dl * 1e6:.1f} µm")
    ax[1].set(xlabel="tube diameter D [µm]", ylabel="relative apparent viscosity μ_app/μ_p",
              title="Fåhræus–Lindqvist: thinner tube, lower μ_app (model, ours)")
    ax[1].legend(fontsize=8)
    save(fig, out, "fl_two_layer")
    for d in (300.0, 100.0, 40.0, 20.0):
        m = ch16.fahraeus_lindqvist_model(d / 2 * 1e-6, delta, 0.45, mu_p, T)
        print(f"  D = {d:5.0f} µm: closure (δ = {delta * 1e6:.1f} µm) μ_rel = {m['mu_rel']:.3f}, H_c = {m['H_c']:.3f}, "
              f"H_T = {m['H_T']:.3f}, residual {m['residual']:.1e}; Pries μ_rel = {ch16.pries_relative_viscosity(d, 0.45):.3f}")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
