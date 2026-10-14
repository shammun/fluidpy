"""Chapter 16, §16.3 — flow in a gently curved pipe (Dean): the curvature ratio (16.164), Dean-number conventions
(16.171), (16.178), (16.179) (trap T8), Dean's small-κ series (16.172)–(16.177) solved exactly order by order (ours:
ψ₀, w₁ and the flux coefficient 1541/50400 behind the printed 0.0306), secondary streamlines and the outward skew of
the axial flow, the flux ratio (16.178)/(16.181) and the resistance ratio (16.182).

Run: ``.venv/Scripts/python.exe scripts/ch16_dean.py --no-show``   Figures -> outputs/ch16/dean_*.png
"""
from __future__ import annotations

import numpy as np
from ch16_common import C_E, C_P, C_REF, C_U, C_VISC, Timer, finish, parse_args, save, setup

from fluidpy import ch16_biofluid_mechanics as ch16


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch16.illustrative_inputs()
    with Timer("exact series to order 4"):
        d2 = ch16.dean_series_sympy(2)
        d4 = ch16.dean_series_sympy(4)
    print(f"Dean series (ours, exact rationals): w₀ = {d2['w0']}, ψ₀ = {d2['psi0']}")
    print(f"  w₁ = {d2['w1']}")
    print(f"  flux coefficient of (K/576)²: {d2['flux_coefficient']} = {float(d2['c2']):.7f} (printed 0.0306); "
          f"of (K/576)⁴: {d4['c4']} = {float(d4['c4']):.7f} (printed 0.0120: our exact solution of (16.173)–(16.176) "
          f"is {100 * (float(d4['c4']) / 0.0120 - 1):+.2f} % from the printed 0.0120; the gap is slip S34, the sign of (16.174))")
    print(f"  resistance (16.182): 0.0306² − 0.0120 = {0.0306 ** 2 - 0.0120:+.5f} (printed −0.0110)")
    delta, a, U, nu = I["dean_delta"], 2e-3, 0.1, I["mu"] / I["rho"]
    kap = ch16.dean_number(delta, U, a, nu, convention="kappa")
    print(f"  conventions (ours: δ = {delta}, a = {a * 1e3} mm, W̄ = {U} m/s, ν = {nu:.2e}): κ = De = {kap:.3f}, "
          f"K = 2κ² = {ch16.dean_number(delta, U, a, nu, convention='K'):.3f}, D = 4√2κ = {ch16.dean_number(delta, U, a, nu, convention='D'):.3f}, "
          f"K = D²/16 check {ch16.dean_convert(ch16.dean_number(delta, U, a, nu, convention='D'), from_='D', to='K'):.3f}")
    # cross-section: secondary streamlines and the skewed axial flow
    K = 400.0
    r = np.linspace(0, 1, 161)
    th = np.linspace(0, 2 * np.pi, 241)
    RR, TT = np.meshgrid(r, th)
    Xc, Yc = RR * np.cos(TT), RR * np.sin(TT)  # α_t = 0 (outer bend) on the +x axis
    psi = ch16.dean_secondary_streamfunction(RR, TT)
    w = ch16.dean_axial_velocity(RR, TT, K)
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.4))
    cf = ax[0].contourf(Xc, Yc, w, levels=20, cmap="Blues")
    ax[0].contour(Xc, Yc, psi, levels=14, colors=C_VISC, linewidths=0.9)
    fig.colorbar(cf, ax=ax[0], label="axial speed w₀ + κ²w₁ (scaled by W̄)")
    ax[0].annotate("outer bend (α = 0)", (1.0, 0.0), xytext=(0.55, -1.15), arrowprops=dict(arrowstyle="->"))
    ax[0].set(aspect="equal", xlabel="r cos α", ylabel="r sin α",
              title=f"K = {K:g}: two counter-rotating cells ψ₀ (16.177); axial max pushed outward")
    imax = np.unravel_index(np.argmax(w), w.shape)
    print(f"  K = {K:g}: axial maximum at r = {RR[imax]:.3f}, α = {np.degrees(TT[imax]):.1f}° (0° = outer bend); "
          f"secondary speed at the centre u = {ch16.dean_secondary_velocity(1e-9, 0.0, K)[0]:.4f} W̄ (toward the outer wall)")
    Ks = np.linspace(0, 900, 200)
    ax[1].plot(Ks, ch16.dean_flux_ratio(Ks), color=C_U, label="(16.178) printed coefficients")
    ax[1].plot(Ks, ch16.dean_flux_ratio(Ks, exact=True, exact_c4=True), color=C_P, ls="--", label="our exact coefficients")
    ax[1].plot(Ks, 1 / ch16.curved_resistance_ratio(Ks), color=C_E, ls=":", label="1/λ from (16.182)")
    rel = np.array([ch16.dean_series_reliability(k)["ratio"] for k in Ks])
    ax[1].axvspan(Ks[np.argmax(rel > 0.1)], Ks[-1], color=C_REF, alpha=0.15, label="K⁴ term > 10 % of K² term (ours)")
    ax[1].set(xlabel="Dean number K = 2κ²", ylabel="Q / Q_straight at equal gradient", title="curvature costs flux at O(K²)")
    ax[1].legend(fontsize=8)
    save(fig, out, "dean")
    for k in (100.0, 300.0, 576.0):
        print(f"  K = {k:5.0f}: Q/Q_s (16.178) = {ch16.dean_flux_ratio(k):.5f}, exact coefficients "
              f"{ch16.dean_flux_ratio(k, exact=True, exact_c4=True):.5f}, λ (16.182) = {ch16.curved_resistance_ratio(k):.5f}, "
              f"series reliability ratio {ch16.dean_series_reliability(k)['ratio']:.3f}")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
