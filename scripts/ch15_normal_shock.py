"""Chapter 15, §15.6 — normal shock waves: the jump relations, the entropy that selects M₁ > 1, the moving shock and
the inside of a shock.

* (15.39)–(15.42) and the three jump laws (15.35)–(15.37) as residuals; M₂ against M₁ (its own inverse);
* entropy jump (15.43), its cubic start (15.44a, b), stagnation-pressure loss (ours), Hugoniot curve against the
  isentrope (ours); the forbidden expansion root; the hydraulic-jump analogy of ch07;
* moving shock: speed of the front and of the gas behind it; strong-shock limits and the crude blast constant;
* shock structure: Becker's profile U(η) for Pr″ = 1 against the ODE, thickness, order-of-magnitude estimate.

Run: ``.venv/Scripts/python.exe scripts/ch15_normal_shock.py --no-show``   Figures -> outputs/ch15/shock_*.png
"""
from __future__ import annotations

import numpy as np
from ch15_common import C_ENT, C_PRES, C_REF, C_SONIC, C_SUB, C_SUP, C_TEMP, COLORS, finish, parse_args, save, setup
from scipy.integrate import solve_ivp

from fluidpy import ch15_compressible_flow as ch15


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch15.illustrative_inputs()
    g, M1 = I["gamma"], I["M1"]
    gas = ch15.PerfectGas(g, I["R"])
    s = ch15.normal_shock(M1, g)
    print(f"normal shock at M₁ = {M1}, γ = {g}: M₂ = {s['M2']:.5f} (15.40), p₂/p₁ = {s['p_ratio']:.5f} (15.39), ρ₂/ρ₁ = {s['rho_ratio']:.5f} (15.41), "
          f"T₂/T₁ = {s['T_ratio']:.5f} (15.42), ΔS/C_v = {s['dS_over_cv']:.5f} (15.43), p₀₂/p₀₁ = {s['p0_ratio']:.5f} (ours)")
    c1 = gas.c(I["T_air"])
    st1 = dict(rho=I["p_air"] / (gas.R * I["T_air"]), u=M1 * c1, p=I["p_air"])
    st2 = dict(rho=st1["rho"] * s["rho_ratio"], u=st1["u"] * s["u_ratio"], p=st1["p"] * s["p_ratio"])
    r = ch15.normal_shock_residuals(st1, st2, g)
    print(f"  jump residuals: mass {r['mass']:.1e}, momentum {r['momentum']:.1e}, energy {r['energy']:.1e}; M₂(M₂(M₁)) = {ch15.normal_shock_mach(s['M2'], g):.6f} (its own inverse)")
    try:
        ch15.normal_shock(s["M2"], g)
    except ValueError as e:
        print(f"  normal_shock({s['M2']:.3f}) refuses: {str(e)[:78]}…")
    fx = ch15.normal_shock(s["M2"], g, allow_expansion=True)
    print(f"  the forbidden root from M₁ = {s['M2']:.4f}: M₂ = {fx['M2']:.4f}, ΔS/C_v = {fx['dS_over_cv']:.5f} < 0, admissible = {fx['admissible']}")
    fl = ch15.normal_shock_fluxes(M1, g)
    print(f"  fluxes downstream/upstream: mass {fl['mass']:.6f}, momentum {fl['momentum']:.6f}, energy {fl['energy']:.6f}; momentum flux is "
          f"{fl['ke1_part']:.1%} ρu² upstream and {fl['p2_part']:.1%} pressure downstream; total enthalpy is {fl['k1_part']:.1%} kinetic upstream, {fl['k2_part']:.1%} downstream")
    print("  weak shocks:  M₁    ΔS/C_v (15.43)   cubic (15.44a)   cubic in Δp/p (15.44b)   1 − p₀₂/p₀₁")
    for m in (1.01, 1.05, 1.1, 1.2, 1.5):
        dpp = ch15.normal_shock_pressure_ratio(m, g) - 1.0
        print(f"                {m:.2f}  {ch15.normal_shock_entropy(m, g):.4e}       {ch15.weak_shock_entropy(m, g):.4e}       {ch15.weak_shock_entropy(dpp, g, form='pressure'):.4e}"
              f"               {1 - ch15.normal_shock_stagnation_pressure_ratio(m, g):.3e}")
    lim = ch15.strong_shock_limits(g)
    print(f"  strong-shock limits: ρ₂/ρ₁ → {lim['rho_ratio']:.3f}, M₂ → {lim['M2']:.5f}; at M₁ = 50: {ch15.normal_shock_density_ratio(50.0, g):.4f}, {ch15.normal_shock_mach(50.0, g):.5f}")
    hj = ch15.hydraulic_jump(1.0, Fr1=2.0)
    print(f"  analogy (ch07): hydraulic jump at Fr₁ = 2: depth ratio {hj['ratio']:.4f}, Fr₂ = {hj['Fr2']:.4f}; shock at M₁ = 2: density ratio {s['rho_ratio']:.4f} "
          f"(the jump loses mechanical energy; the shock conserves energy and loses stagnation pressure)")

    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
    Ma = np.linspace(0.39, 4.0, 600)
    ax[0].plot(Ma[Ma >= 1], ch15.normal_shock_mach(Ma[Ma >= 1], g), color=C_SUP, label="shock (15.40), M₁ > 1")
    ax[0].plot(Ma[Ma <= 1], ch15.normal_shock_mach(Ma[Ma <= 1], g), color=C_ENT, ls="--", label="forbidden: ΔS < 0")
    ax[0].plot([0, 4], [0, 4], color=C_REF, lw=1, label="no change, M₂ = M₁")
    ax[0].axhline(lim["M2"], color=C_SONIC, lw=0.8, ls=":")
    ax[0].set(xlabel="upstream Mach number M₁ [–]", ylabel="M₂ [–]", xlim=(0, 4), ylim=(0, 4), title="M₂ against M₁: a curve that is its own mirror image")
    ax[0].legend(fontsize=8)
    Mb = np.linspace(0.5, 5.0, 500)
    ax[1].plot(Mb, ch15.normal_shock_entropy(Mb, g), color=C_ENT, label=r"$(S_2-S_1)/C_v$ (15.43)")
    Mc = np.linspace(0.8, 1.6, 200)
    ax[1].plot(Mc, ch15.weak_shock_entropy(Mc, g), color=COLORS["ink"], ls="--", lw=1, label="cubic (15.44a)")
    ax[1].plot(Mb[Mb >= 1], 1 - ch15.normal_shock_stagnation_pressure_ratio(Mb[Mb >= 1], g), color=C_PRES, label=r"$1-p_{02}/p_{01}$ (ours)")
    ax[1].axhline(0.0, color=C_REF, lw=0.8)
    ax[1].axvline(1.0, color=C_SONIC, lw=0.8, ls=":")
    ax[1].set(xlabel="M₁ [–]", ylabel="[–]", ylim=(-0.5, 1.5), title="entropy jump: negative below M₁ = 1, third-order contact at 1")
    ax[1].legend(fontsize=8)
    rr = np.linspace(0.4, 5.5, 400)
    ax[2].semilogy(rr, ch15.hugoniot_pressure_ratio(rr, g), color=C_SUP, label="Hugoniot (shock) curve (ours)")
    ax[2].semilogy(rr, rr ** g, color=C_SUB, ls="--", label=r"isentrope $(\rho_2/\rho_1)^\gamma$ (15.2)")
    ax[2].axvline(lim["rho_ratio"], color=C_REF, lw=0.8, ls=":")
    ax[2].set(xlabel=r"$\rho_2/\rho_1$ [–]", ylabel=r"$p_2/p_1$ [–]", ylim=(0.2, 200), title="shock compression against isentropic compression")
    ax[2].legend(fontsize=8)
    save(fig, out, "shock_jump_entropy_hugoniot")

    # --- moving shock ------------------------------------------------------------------------------------------------------
    mv = ch15.moving_shock(I["p_air"], I["T_air"], M_s=I["M_s"], gas=gas)
    print(f"moving shock (ours: M_s = {I['M_s']} into air at {I['T_air']} K): front {mv['u_shock']:.1f} m/s (c₁ = {mv['c1']:.1f}), gas behind follows at {mv['u_gas']:.1f} m/s "
          f"(Mach {mv['M_gas']:.3f} in its own sound speed), p₂ = {mv['p2'] / 1e3:.1f} kPa, T₂ = {mv['T2']:.1f} K; T₀ behind in the laboratory frame {mv['T0_behind_lab']:.1f} K "
          f"(> T₁: stagnation values depend on the frame)")
    wk = ch15.moving_shock(I["p_air"], I["T_air"], p2_over_p1=1.001, gas=gas)
    print(f"  weak limit (Δp/p = 0.001): front {wk['u_shock']:.3f} m/s → c₁; gas speed {wk['u_gas']:.5f} m/s vs acoustic Δp/(ρ₁c₁) = {0.001 * I['p_air'] / (wk['rho1'] * wk['c1']):.5f} m/s (15.16)")
    R_ = ch15.riemann_exact((1.0, 0.0, 1.0), (0.125, 0.0, 0.1), g, 0.0)
    print(f"  exact Riemann solver (ours) on a 10:1 pressure, 8:1 density jump: p* = {R_['p_star']:.5f}, u* = {R_['u_star']:.5f}, shock speed {R_['shock_speed_right']:.5f} "
          f"(= moving_shock: {ch15.moving_shock(0.1, 0.8, p2_over_p1=R_['p_star'] / 0.1, gas=ch15.PerfectGas(g, 1.0))['u_shock']:.5f})")
    print(f"  crude blast constant (Exercise 15.9 model, ours): K = {ch15.blast_constant_crude(g):.4f} (hemisphere), {ch15.blast_constant_crude(g, 'sphere'):.4f} (sphere)")

    # --- shock structure ---------------------------------------------------------------------------------------------------
    Ms = I["M1_structure"]
    U1, U2 = ch15.becker_shock_roots(Ms, g)
    k = ch15.shock_structure_constants(Ms, I["T_air"], I["p_air"], gas)
    eta = np.linspace(-12.0, 12.0, 121 if args.fast else 481)
    U = ch15.becker_shock_profile(eta, Ms, g)
    i0 = eta.size // 2
    fw = solve_ivp(ch15.becker_shock_ode, (0.0, 12.0), [U[i0]], args=(Ms, g), rtol=1e-10, atol=1e-12, t_eval=eta[i0:])
    print(f"shock structure (Becker, Pr″ = 1) at M₁ = {Ms}: plateaus U₁ = {U1:.6f}, U₂ = {U2:.6f} (sum {U1 + U2:.6f} = 2γ/(γ + 1) = {2 * g / (g + 1):.6f}); "
          f"m = {k['m']:.2f} kg/(m² s), V = {k['V']:.2f} m/s, I = {k['I']:.0f} J/kg; implicit profile vs solve_ivp: {np.max(np.abs(fw.y[0] - U[i0:])):.1e}; "
          f"maximum-slope thickness Δη = {ch15.becker_shock_thickness(Ms, g):.4f}")
    du = k["u1"] - k["u2"]
    est = ch15.shock_thickness_estimate(du, 1.5e-5)
    print(f"  estimate (u₁ − u₂)δ/ν ~ 1 with ν = 1.5e-5 m²/s, Δu = {du:.0f} m/s: δ ~ {est:.1e} m = {est / ch15.mean_free_path_air(I['T_air'], I['p_air']):.1f} mean free paths "
          f"(λ = {ch15.mean_free_path_air(I['T_air'], I['p_air']):.2e} m): the continuum model is at its limit")
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    for m, col in zip((1.2, 1.5, 2.0, 3.0), (C_SUB, C_TEMP, C_SUP, C_ENT)):
        a, b = ch15.becker_shock_roots(m, g)
        ax[0].plot(eta, (ch15.becker_shock_profile(eta, m, g) - b) / (a - b), color=col, label=f"M₁ = {m}: Δη = {ch15.becker_shock_thickness(m, g):.1f}")
    ax[0].set(xlabel=r"$\eta = m\int dx/\mu''$ [–]", ylabel=r"$(U-U_2)/(U_1-U_2)$ [–]", title="inside a shock (Becker): stronger shocks are thinner")
    ax[0].legend(fontsize=8)
    ax[1].plot(eta, U, color=C_SUP, label="implicit solution (brentq)")
    ax[1].plot(fw.t, fw.y[0], color=COLORS["ink"], ls="--", lw=1, label="solve_ivp on the structure ODE")
    ax[1].axhline(U1, color=C_REF, lw=0.8, ls=":")
    ax[1].axhline(U2, color=C_REF, lw=0.8, ls=":")
    ax[1].set(xlabel=r"$\eta$ [–]", ylabel=r"$U = u/V$ [–]", title=f"velocity profile at M₁ = {Ms}")
    ax[1].legend(fontsize=8)
    save(fig, out, "shock_structure_becker")
    rs = ch15.normal_shock_sympy()
    ws = ch15.weak_shock_series()
    print(f"sympy: roots of the jump laws ρ₂/ρ₁ ∈ {rs['roots']}; residuals (15.39)–(15.42): {rs['residual_15_39']}, {rs['residual_15_40']}, {rs['residual_15_41']}, {rs['residual_15_42']}; "
          f"series of (15.43) in m = M₁² − 1: coefficients {ws['coefficients'][:3]}, residuals (15.44a) {ws['residual_15_44a']}, (15.44b) {ws['residual_15_44b']}")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
