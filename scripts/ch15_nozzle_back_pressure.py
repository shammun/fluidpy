"""Chapter 15, §15.7 — nozzles at different back pressures: choking of a convergent nozzle and the regimes of a
convergent–divergent (Laval) nozzle.  The section has no numbered equation: every formula here is ours, built from
(15.29), (15.31) and the normal-shock relations.

* convergent nozzle: mass flow against back pressure — rises, then stays at the choked value below p*/p₀;
* Laval nozzle: the three dividing back pressures p_b > p_d > p_f, the regime names, the shock position for a given
  back pressure and the pressure distribution along the nozzle;
* a nozzle with a shock at a known Mach number (the chain of Example 15.2, our inputs).

Run: ``.venv/Scripts/python.exe scripts/ch15_nozzle_back_pressure.py --no-show``   Figures -> outputs/ch15/backpressure_*.png
"""
from __future__ import annotations

import numpy as np
from ch15_common import C_ENT, C_PRES, C_REF, C_SONIC, C_SUB, C_SUP, C_TEMP, COLORS, finish, parse_args, save, setup

from fluidpy import ch15_compressible_flow as ch15


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch15.illustrative_inputs()
    g = I["gamma"]
    gas = ch15.PerfectGas(g, I["R"])
    crit = ch15.critical_pressure_ratio(g)
    print(f"convergent nozzle (ours: p0 = {I['p0'] / 1e3:.0f} kPa, T0 = {I['T0']} K, exit {1e4 * I['A_conv_exit']:.0f} cm²): critical ratio p*/p0 = {crit:.5f}")
    print("   p_B/p0   M_exit    p_exit/p0   mdot [kg/s]  choked  matched")
    for r in (1.0, 0.95, 0.8, 0.6, crit, 0.4, 0.1):
        c = ch15.convergent_nozzle(I["p0"], I["T0"], r * I["p0"], I["A_conv_exit"], gas)
        print(f"   {r:.4f}   {c['M_exit']:.5f}   {c['p_exit'] / I['p0']:.5f}     {c['mdot']:.5f}     {c['choked']!s:5}   {c['matched']}")
    e = 1e-6
    sl = [(ch15.mass_flow_vs_back_pressure(crit + s * e, g)["mdot"] - ch15.mass_flow_vs_back_pressure(crit, g)["mdot"]) / (s * e) for s in (+1, -1)]
    print(f"  slope of the dimensionless mass flow at the critical ratio: {sl[0]:.2e} from above, {sl[1]:.2e} from below (the curve arrives with zero slope)")

    ar = I["Ae_over_At"]
    cp = ch15.laval_critical_pressures(ar, g)
    print(f"Laval nozzle with A_e/A_t = {ar}: p_b/p0 = {cp['p_b']:.5f} (subsonic, sonic throat; M_e = {cp['M_sub']:.4f}), p_d/p0 = {cp['p_d']:.5f} (shock at the exit), "
          f"p_f/p0 = {cp['p_f']:.5f} (design; M_e = {cp['M_sup']:.4f})")
    print("   p_B/p0   regime            M_exit   p_exit/p0  M1_shock  A_shock/A_t  p02/p01   mdot/choked")
    for r in (0.99, cp["p_b"], 0.85, I["pB_over_p0_shock"], 0.6, cp["p_d"], 0.3, cp["p_f"], 0.03):
        s = ch15.laval_exit_state(r, ar, g)
        print(f"   {r:.4f}   {s['regime']:16}  {s['M_exit']:.4f}   {s['p_exit_over_p0']:.5f}   {s['M1_shock']:.4f}    {s['A_shock_over_At']:.4f}       {s['p02_over_p01']:.4f}    {s['mdot_over_choked']:.4f}")

    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
    rr = np.linspace(0.0, 1.0, 401)
    mf = ch15.mass_flow_vs_back_pressure(rr, g)
    ax[0].plot(rr, mf["mdot_over_choked"], color=C_SUB, label="mass flow / choked mass flow")
    ax[0].plot(rr, mf["p_exit_over_p0"], color=C_PRES, label=r"exit-plane pressure $p_e/p_0$")
    ax[0].plot(rr, rr, color=C_REF, lw=0.8, ls=":", label=r"$p_e = p_B$")
    ax[0].axvline(crit, color=C_SONIC, ls="--", lw=1, label=rf"$p^*/p_0$ = {crit:.3f}")
    ax[0].invert_xaxis()
    ax[0].set(xlabel=r"back pressure $p_B/p_0$ [–] (lowered →)", ylabel="[–]", title="convergent nozzle: choked below the critical ratio")
    ax[0].legend(fontsize=8)

    n = 161 if args.fast else 321
    x = np.linspace(0.0, I["L"], n)
    A = lambda xx: ch15.laval_area(xx, I["A_throat"], I["A_inlet"], I["A_exit"], I["x_throat"], I["L"])  # noqa: E731
    cases = [(0.985, C_SUB), (cp["p_b"], C_TEMP), (0.85, C_ENT), (I["pB_over_p0_shock"], C_ENT), (0.62, C_ENT), (cp["p_d"], C_PRES), (cp["p_f"], C_SUP)]
    for r, col in cases:
        L_ = ch15.laval_nozzle(x, A, r, g)
        ax[1].plot(x, L_["p_over_p0"], color=col, lw=1.4)
        ax[2].plot(x, L_["M"], color=col, lw=1.4, label=f"{r:.3f}: {L_['regime']}")
        if L_["regime"] == "shock at exit":  # the jump in the exit plane
            ax[1].plot([x[-1], x[-1]], [cp["p_f"], cp["p_d"]], color=col, lw=1.4)
    ax[1].axvline(I["x_throat"], color=C_SONIC, lw=0.8, ls="--")
    ax[1].axhline(crit, color=C_SONIC, lw=0.8, ls=":")
    ax[1].set(xlabel="x [m]", ylabel=r"$p/p_0$ [–]", title="pressure along the Laval nozzle for seven back pressures")
    ax[2].axhline(1.0, color=C_SONIC, lw=0.8, ls=":")
    ax[2].set(xlabel="x [m]", ylabel="M [–]", title="Mach number: the back pressure only moves the shock")
    ax[2].legend(fontsize=7, title=r"$p_B/p_0$: regime")
    save(fig, out, "backpressure_convergent_laval")

    L_ = ch15.laval_nozzle(x, A, I["pB_over_p0_shock"], g)
    md = ch15.mass_flux_function(L_["M"], g) * L_["A"] * np.where(x > L_["x_shock"], L_["p02_over_p01"], 1.0)
    print(f"shock in the nozzle at p_B/p0 = {I['pB_over_p0_shock']}: M₁ = {L_['M1_shock']:.5f}, A_shock/A_t = {L_['A_shock_over_At']:.5f}, x_shock = {L_['x_shock']:.5f} m, "
          f"p02/p01 = {L_['p02_over_p01']:.5f}, exit M = {L_['M_exit']:.5f}, exit p/p0 = {L_['p_over_p0'][-1]:.6f}; mass flow along x constant to {np.ptp(md) / md.mean():.1e}")
    ws = ch15.nozzle_with_shock(L_["M1_shock"], I["A_throat"], I["A_exit"], I["p0"], I["T0"], gas)
    print(f"  the forward chain (Example 15.2's, our inputs): A_shock = {1e4 * ws['A_shock']:.3f} cm², M₂ = {ws['M2']:.5f}, A₂* = {1e4 * ws['A_star2']:.3f} cm² (> A_t = {1e4 * I['A_throat']:.0f}), "
          f"M₃ = {ws['M3']:.5f}, T₃ = {ws['T3']:.2f} K, p₃/p0 = {ws['p3'] / I['p0']:.6f} (= the back pressure); ṁ throat {ws['mdot']:.4f} = exit {ws['mdot_exit']:.4f} kg/s")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
