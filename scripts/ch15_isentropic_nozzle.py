"""Chapter 15, §15.3–§15.5 — the duct laws, the stagnation and sonic states, the area ratio and the area–velocity rule.

* stagnation ratios (15.28)–(15.30) and the compressible correction to a pitot reading; sonic ratios (ours);
* the area ratio (15.31): a minimum of 1 at M = 1, so every area ratio belongs to two Mach numbers; the mass-flux
  function and the choked mass flow (ours);
* the area–velocity rule (15.34) with the corrected (15.32); the printed (15.32) fails the finite-difference check;
* the three solutions in one convergent–divergent duct (subsonic venturi, transonic, supersonic) and the conserved
  quantities along them; wall force (15.24) against ∫p dA;
* a pressure-matched rocket nozzle sized by the chain of Example 15.1 (our inputs).

Run: ``.venv/Scripts/python.exe scripts/ch15_isentropic_nozzle.py --no-show``   Figures -> outputs/ch15/nozzle_*.png
"""
from __future__ import annotations

import numpy as np
from ch15_common import C_PRES, C_REF, C_SONIC, C_SUB, C_SUP, C_TEMP, COLORS, finish, parse_args, save, setup

from fluidpy import ch15_compressible_flow as ch15


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch15.illustrative_inputs()
    g = I["gamma"]
    gas = ch15.PerfectGas(g, I["R"])
    print(f"gas: γ = {g}, R = {gas.R:.3f}, C_p = {gas.cp:.2f}, C_v = {gas.cv:.2f} J/(kg K) (15.1); book-rounded air: R = {ch15.AIR_BOOK.R}, C_p = {ch15.AIR_BOOK.cp:.1f}")
    print("  M     T/T0     p/p0     rho/rho0  A/A*     (p0 − p)/(½ρu²)")
    for M in (0.3, 0.5, 1.0, 2.0, 3.0):
        s = ch15.isentropic_state(M, g)
        print(f"  {M:.1f}   {s['T_over_T0']:.5f}  {s['p_over_p0']:.5f}  {s['rho_over_rho0']:.5f}   {s['A_over_Astar']:.5f}  {ch15.pitot_compressibility_factor(M, g):.5f}")
    sr = ch15.sonic_ratios(g)
    print(f"  sonic state (ours): T*/T0 = {sr['T']:.5f}, p*/p0 = {sr['p']:.5f}, ρ*/ρ0 = {sr['rho']:.5f}; flux function at M = 1: {ch15.mass_flux_function(1.0, g):.5f}")
    print(f"  density changes by 5 % at M = {ch15.mach_for_density_change(0.05, g):.4f} (ch14's ½M² rule: {ch15.compressibility_check(0.3 * gas.c(288.15))['density_change']:.4f} at M = 0.3; exact {ch15.stagnation_density_ratio(0.3, g) - 1:.4f})")
    ar = I["Ae_over_At"]
    print(f"  area ratio {ar}: M_sub = {ch15.mach_from_area_ratio(ar, g, branch='subsonic'):.6f}, M_sup = {ch15.mach_from_area_ratio(ar, g, branch='supersonic'):.6f}")
    mdot = ch15.choked_mass_flow(I["p0"], I["T0"], I["A_throat"], gas)
    print(f"  our nozzle: p0 = {I['p0'] / 1e3:.0f} kPa, T0 = {I['T0']} K, A_t = {1e4 * I['A_throat']:.0f} cm²: choked mass flow {mdot:.4f} kg/s; u_max = √(2C_pT0) = {ch15.max_speed(I['T0'], gas.cp):.1f} m/s")

    # --- influence coefficients and the slip in (15.32) --------------------------------------------------------------------
    print("  area rule du/u per dA/A:  M | (15.34) | finite difference of (15.31) | printed (15.32) used instead")
    for M in (0.5, 0.9, 2.0):
        h = 1e-6
        dlnA = np.log(ch15.area_ratio(M * (1 + h), g) / ch15.area_ratio(M, g))
        c = ch15.influence_coefficients(M, g, dA_over_A=1.0)
        fd = h / dlnA  # dM/M per dA/A
        cp_ = ch15.influence_coefficients(M, g, dA_over_A=1.0, printed_15_32=True)
        print(f"    {M:.1f} | du/u = {c['du_over_u']:+.5f}, dM/M = {c['dM_over_M']:+.5f} | dM/M = {fd:+.5f} | du/u = {cp_['du_over_u']:+.5f}  ({ch15.duct_role(M, -1)} if converging, {ch15.duct_role(M, +1)} if diverging)")

    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
    M = np.linspace(0.05, 3.5, 500)
    s = ch15.isentropic_state(M, g)
    ax[0].plot(M, s["T_over_T0"], color=C_TEMP, label=r"$T/T_0$ (15.28)")
    ax[0].plot(M, s["p_over_p0"], color=C_PRES, label=r"$p/p_0$ (15.29)")
    ax[0].plot(M, s["rho_over_rho0"], color=C_SUB, label=r"$\rho/\rho_0$ (15.30)")
    ax[0].axvline(1.0, color=C_SONIC, lw=1, ls="--")
    ax[0].set(xlabel="Mach number M [–]", ylabel="ratio to the stagnation value [–]", title=f"stagnation ratios, γ = {g}")
    ax[0].legend(fontsize=8)
    ax[1].plot(M[M <= 1], s["A_over_Astar"][M <= 1], color=C_SUB, label="subsonic branch")
    ax[1].plot(M[M >= 1], s["A_over_Astar"][M >= 1], color=C_SUP, label="supersonic branch")
    ax[1].axhline(ar, color=C_REF, lw=0.8, ls=":")
    for br, col in (("subsonic", C_SUB), ("supersonic", C_SUP)):
        ax[1].plot([ch15.mach_from_area_ratio(ar, g, branch=br)], [ar], "o", color=col)
    ax[1].plot(M, ch15.mass_flux_function(M, g) / ch15.mass_flux_function(1.0, g), color=C_SONIC, ls="--", label="mass flux / choked flux (ours)")
    ax[1].set(xlabel="Mach number M [–]", ylabel=r"$A/A^*$ [–]", ylim=(0, 6), title="area ratio (15.31): one area, two Mach numbers")
    ax[1].legend(fontsize=8)

    # --- three solutions in one duct ------------------------------------------------------------------------------------------
    n = 121 if args.fast else 241
    x = np.linspace(0.0, I["L"], n)
    A = lambda xx: ch15.laval_area(xx, I["A_throat"], I["A_inlet"], I["A_exit"], I["x_throat"], I["L"])  # noqa: E731
    sols = {"transonic": ch15.nozzle_isentropic(x, A, g, branch="transonic"),
            "subsonic, sonic throat": ch15.nozzle_isentropic(x, A, g, branch="subsonic"),
            "venturi (A* = 0.8 A_t)": ch15.nozzle_isentropic(x, A, g, branch="subsonic", A_star=0.8 * I["A_throat"]),
            "supersonic throughout": ch15.nozzle_isentropic(x, A, g, branch="supersonic")}
    for (lab, so), col in zip(sols.items(), (COLORS["ink"], C_SUB, C_TEMP, C_SUP)):
        ax[2].plot(x, so["M"], color=col, label=lab)
    ax[2].axvline(I["x_throat"], color=C_SONIC, lw=1, ls="--")
    ax[2].set(xlabel="x [m]", ylabel="M [–]", title="one duct, four isentropic flows: M = 1 only at the throat")
    ax[2].legend(fontsize=8)
    save(fig, out, "nozzle_isentropic_ratios_area")
    so = sols["transonic"]
    rho0 = I["p0"] / (gas.R * I["T0"])
    T = so["T_over_T0"] * I["T0"]
    u = so["M"] * gas.c(T)
    rho, p = so["rho_over_rho0"] * rho0, so["p_over_p0"] * I["p0"]
    md = ch15.mass_flow(rho, u, so["A"])
    h0 = ch15.total_enthalpy(T, u, gas.cp)
    print(f"  transonic solution: ṁ = ρuA constant to {np.ptp(md) / md.mean():.1e} ({md.mean():.4f} kg/s = choked {mdot:.4f}); h₀ constant to {np.ptp(h0) / h0.mean():.1e}; "
          f"exit M = {so['M'][-1]:.4f}, p_e/p0 = {so['p_over_p0'][-1]:.5f}; venturi throat Mach number {sols['venturi (A* = 0.8 A_t)']['M'][so['i_throat']]:.4f} (not sonic)")
    s1 = dict(rho=rho[0], u=u[0], p=p[0])
    s2 = dict(rho=rho[-1], u=u[-1], p=p[-1])
    Fw = ch15.wall_force(s1, s2, so["A"][0], so["A"][-1])
    Fq = float(np.sum(0.5 * (p[1:] + p[:-1]) * np.diff(so["A"])))
    print(f"  wall force on the fluid (15.24): {Fw:.2f} N; ∫p dA along the wall (trapezoid, {n} stations): {Fq:.2f} N; energy residual (15.22): {ch15.energy_residual(s1, s2, 0.0, g):.2e} J/kg")

    # --- rocket nozzle (the chain of Example 15.1, our inputs) --------------------------------------------------------------
    gr = ch15.PerfectGas(I["gamma_rocket"], I["R_rocket"])
    rk = ch15.rocket_nozzle(I["thrust"], I["p0_rocket"], I["T0_rocket"], I["p_exit_rocket"], gr)
    print(f"rocket nozzle (ours: thrust {I['thrust'] / 1e3:.0f} kN, p0 = {I['p0_rocket'] / 1e6:.1f} MPa, T0 = {I['T0_rocket']} K, γ = {gr.gamma}, R = {gr.R}): "
          f"M_e = {rk['M_e']:.4f}, T_e = {rk['T_e']:.1f} K, u_e = {rk['u_e']:.1f} m/s, A_e = {1e4 * rk['A_e']:.2f} cm², A* = {1e4 * rk['A_throat']:.2f} cm² "
          f"(A_e/A* = {rk['Ae_over_At']:.4f}), ṁ = {rk['mdot']:.3f} kg/s = choked {rk['mdot_choked']:.3f}; thrust check ṁu_e = "
          f"{ch15.nozzle_thrust(rk['mdot'], rk['u_e'], I['p_exit_rocket'], I['p_exit_rocket'], rk['A_e']):.1f} N; in vacuum the same nozzle gives "
          f"{ch15.nozzle_thrust(rk['mdot'], rk['u_e'], I['p_exit_rocket'], 0.0, rk['A_e']):.1f} N (pressure term)")
    r = ch15.duct_equations_sympy()
    print(f"sympy: (15.19) {r['eq_15_19']}, (15.25) {r['eq_15_25']}, (15.26) {r['eq_15_26']}, (15.34) {r['eq_15_34']}; printed (15.32) would leave {r['eq_15_33_printed']}")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
