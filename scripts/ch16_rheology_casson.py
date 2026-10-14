"""Chapter 16, §16.2 and §16.3 — blood as a Casson fluid: the yield-stress correlation (16.5), the Casson law (16.6),
Whitmore's fit (16.7) (Figs. 16.7–16.8 redrawn with our curves, no data), apparent (secant) against differential
(tangent) viscosity, and Casson flow in a pipe (16.207)–(16.220): the plug, the profile, the flux.

Run: ``.venv/Scripts/python.exe scripts/ch16_rheology_casson.py --no-show``   Figures -> outputs/ch16/casson_*.png
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
    mu_p = I["mu_plasma"]
    print(f"units (slip S27 corrected): 1 cP = {ch16.cP_to_Pa_s(1.0):g} Pa s = 0.01 dyn s/cm²; 1 P = {ch16.poise_to_Pa_s(1.0):g} Pa s")
    for H in (0.40, 0.45, 0.50):
        print(f"  yield stress (16.5): H = {H:.2f}, C_F = {I['C_F']:.2f} g/100 mL -> τ_y = "
              f"{ch16.yield_stress_correlation(H, I['C_F']) * 1e3:.3f} mPa ({ch16.yield_stress_correlation(H, I['C_F']) / 0.1:.4f} dyn/cm²)")
    print(f"  (16.5) over the stated C_F range at H = 0.45: {ch16.yield_stress_correlation(0.45, 0.21) / 0.1:.4f}–"
          f"{ch16.yield_stress_correlation(0.45, 0.46) / 0.1:.4f} dyn/cm²")
    wh = ch16.whitmore_casson(100.0, mu_p)
    print(f"Whitmore (16.7) with μ_p = {mu_p * 1e3:.2f} mPa s: τ_y = 4μ_p = {wh['tau_y'] * 1e3:.2f} mPa, K_c = 1.53√μ_p = "
          f"{wh['Kc']:.4f} (Pa s)^½, K_c² = {wh['Kc'] ** 2 * 1e3:.3f} mPa s = {wh['Kc'] ** 2 / mu_p:.4f} μ_p; "
          f"at γ̇ = 100 s⁻¹: τ = {wh['tau']:.4f} Pa, μ_app = {wh['mu_app'] * 1e3:.3f} mPa s (secant), dτ/dγ̇ = {wh['mu_tangent'] * 1e3:.3f} mPa s (tangent)")
    g = np.geomspace(0.05, 2000.0, 400)
    tau_y, Kc = wh["tau_y"], wh["Kc"]
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    sq = np.linspace(0, 15, 100)
    ax[0].plot(sq, 1.53 * sq + 2.0, color=C_TAU, label="Whitmore (16.7): √(τ/μ_p) = 1.53√γ̇ + 2.0")
    ax[0].plot(sq, 1.53 * sq, color=C_REF, ls=":", label="no yield stress (Newtonian, μ = K_c²)")
    ax[0].set(xlabel="√γ̇ [s^-1/2]", ylabel="√(τ/μ_p) [s^-1/2]", title="Casson form (16.6): a straight line, intercept √(τ_y/μ_p)")
    ax[0].legend()
    ax[1].semilogx(g, ch16.casson_apparent_viscosity(g, tau_y, Kc) * 1e3, color=C_TAU, label="μ_app = τ/γ̇ (secant)")
    ax[1].semilogx(g, ch16.differential_viscosity(g, tau_y, Kc) * 1e3, color=C_E, ls="--", label="dτ/dγ̇ (tangent, trap T27)")
    ax[1].axhline(Kc ** 2 * 1e3, color=C_REF, ls=":", label=f"asymptote K_c² = {Kc ** 2 * 1e3:.2f} mPa s")
    ax[1].set(xlabel="shear rate γ̇ [1/s]", ylabel="viscosity [mPa s]", ylim=(0, 25),
              title="apparent viscosity of the Whitmore fit (ours, Fig. 16.8 style)")
    ax[1].legend()
    save(fig, out, "casson_law")
    # --- Casson pipe flow -----------------------------------------------------------------------------------------------
    a, L = 0.5e-3, 0.05
    ty, Kc2 = I["tau_y"], I["Kc"]
    dpmin = ch16.casson_min_pressure_drop(L, a, ty)
    print(f"Casson pipe (illustrative τ_y = {ty * 1e3:.1f} mPa, K_c = {Kc2} (Pa s)^½, a = {a * 1e3:.1f} mm, L = {L * 100:.0f} cm):")
    print(f"  minimum pressure drop 2Lτ_y/a = {dpmin:.3f} Pa = {ch16.Pa_to_mmHg(dpmin):.5f} mmHg = {dpmin / 0.1:.2f} dyn/cm² "
          "(the book's worked number has a ×1000 slip in its dyn/cm² figure, S26)")
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    r = np.linspace(0, a, 401)
    for dp, col in ((2.0 * dpmin, C_TAU), (5.0 * dpmin, C_E), (20.0 * dpmin, C_U)):
        ac = ch16.casson_plug_radius(dp, L, ty)
        u = ch16.casson_pipe_velocity(r, a, dp, L, ty, Kc2)
        uP = ch16.poiseuille_velocity(r, a, dp, L, Kc2 ** 2)
        Qc = ch16.casson_pipe_flow_rate(a, dp, L, ty, Kc2)
        Qq = np.pi * ac ** 2 * ch16.casson_core_velocity(a, dp, L, ty, Kc2) + quad(
            lambda rr: 2 * np.pi * rr * ch16.casson_pipe_velocity(rr, a, dp, L, ty, Kc2), ac, a, epsrel=1e-12)[0]
        print(f"  Δp = {dp:.3f} Pa: τ_w = {ch16.casson_wall_stress(dp, L, a) * 1e3:.3f} mPa, a_c/a = τ_y/τ_w = {ac / a:.3f}, "
              f"u_c (16.218) = {ch16.casson_core_velocity(a, dp, L, ty, Kc2) * 1e3:.4f} mm/s, Q (16.220) = {Qc:.5e} m³/s "
              f"(quadrature {Qq:.5e}), Q/Q_Pois(μ = K_c²) = {Qc / ch16.poiseuille_flow_rate(a, dp, L, Kc2 ** 2):.4f}")
        ax[0].plot(np.concatenate([-r[::-1], r]) * 1e3, np.concatenate([u[::-1], u]) * 1e3, color=col,
                   label=f"Casson, a_c/a = {ac / a:.2f}")
        ax[0].plot(np.concatenate([-r[::-1], r]) * 1e3, np.concatenate([uP[::-1], uP]) * 1e3, color=col, ls=":", lw=1)
    ax[0].set(xlabel="r [mm] (full diameter)", ylabel="u [mm/s]",
              title="Casson profiles (16.215)–(16.218) and Poiseuille with μ = K_c² (dotted)")
    ax[0].legend()
    xi = np.linspace(0, 1, 201)
    ax[1].plot(xi, ch16.casson_flux_ratio(xi), color=C_WALL, label="1 − 16√ξ/7 + 4ξ/3 − ξ⁴/21 (16.220)")
    ax[1].plot(xi, 1 - 16.0 / 7.0 * np.sqrt(xi), color=C_REF, ls=":", label="first loss term 16√ξ/7")
    ax[1].set(xlabel="ξ = a_c/a = τ_y/τ_w", ylabel="Q / Q_Poiseuille", ylim=(0, 1), title="the plug costs flux")
    ax[1].legend()
    save(fig, out, "casson_pipe")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
