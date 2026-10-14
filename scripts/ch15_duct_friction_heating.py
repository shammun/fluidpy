"""Chapter 15, §15.8 — friction and heating in a constant-area duct: one quadratic for M₂² (15.46), the differential
forms (15.47), (15.48) with their two printed slips corrected, the temperature paradox and choking.

* (15.45) as residuals on both roots of (15.46); f = q = 0 gives "no change" and the normal shock;
* M₂ against M₁ with f as the parameter (q = 0) and with q as the parameter (f = 0), choking loci (ours);
* (15.47), (15.48): corrected coefficients against finite differences of (15.46); the printed forms fail;
* the usual Fanno and Rayleigh forms (not in the book — ours) and their T–s lines; heating that cools.

Run: ``.venv/Scripts/python.exe scripts/ch15_duct_friction_heating.py --no-show``   Figures -> outputs/ch15/duct_*.png
"""
from __future__ import annotations

import numpy as np
from ch15_common import C_ENT, C_REF, C_SONIC, C_SUB, C_SUP, C_TEMP, COLORS, finish, parse_args, save, setup

from fluidpy import ch15_compressible_flow as ch15


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch15.illustrative_inputs()
    g, f, q = I["gamma"], I["f"], I["q"]
    for M1 in (I["M1_duct_sub"], I["M1_duct_sup"]):
        A = ch15.duct_parameter_A(M1, f, q, g)
        print(f"duct with f = {f}, q = {q} at M₁ = {M1}: A = {A:.6f}, discriminant 1 − 2A(γ + 1) = {1 - 2 * A * (g + 1):.6f}")
        for br in ("continuous", "shocked"):
            M2 = ch15.duct_outlet_mach(M1, f, q, g, branch=br)
            r = ch15.duct_state_ratios(M1, M2, f, q, g)
            s1 = dict(rho=1.0, p=1.0, u=M1 * np.sqrt(g))
            s2 = dict(rho=r["rho_ratio"], p=r["p_ratio"], u=s1["u"] * r["u_ratio"])
            res = ch15.duct_residuals(s1, s2, f, q, g)
            print(f"  {br:10}: M₂ = {M2:.6f}, p₂/p₁ = {r['p_ratio']:.5f}, T₂/T₁ = {r['T_ratio']:.5f}, p₀₂/p₀₁ = {r['p0_ratio']:.5f}, T₀₂/T₀₁ = {r['T0_ratio']:.5f}, "
                  f"ΔS/C_v = {r['dS_over_cv']:+.5f}; residuals of (15.45): {res['mass']:.0e}, {res['momentum']:.0e}, {res['energy']:.0e}")
    print(f"  f = q = 0: continuous root {ch15.duct_outlet_mach([0.5, 2.0], 0, 0, g)}, other root {ch15.duct_outlet_mach([0.5, 2.0], 0, 0, g, branch='shocked')} "
          f"= normal-shock map {ch15.normal_shock_mach([0.5, 2.0], g)}; ΔS/C_v of the supersonic root from M₁ = 0.5: "
          f"{ch15.duct_entropy_change(0.5, ch15.duct_outlet_mach(0.5, 0, 0, g, branch='shocked'), 0, 0, g):+.4f} (excluded)")
    print(f"  printed signs of (15.46) at M₁ = 0.5 and 2 (f = {f}): '+' → {ch15.duct_outlet_mach([0.5, 2.0], f, 0, g, branch='plus')}, '−' → {ch15.duct_outlet_mach([0.5, 2.0], f, 0, g, branch='minus')}")
    print("  choking limits (ours):  M₁    f_max     q_max    M₂ at f_max   M₂ just beyond")
    for M1 in (0.3, 0.5, 0.8, 1.5, 2.0, 3.0):
        lim = ch15.duct_choking_limits(M1, g)
        print(f"                          {M1:.1f}   {lim['f_max']:.5f}   {lim['q_max']:.5f}   {ch15.duct_outlet_mach(M1, lim['f_max'] * (1 - 1e-10), 0, g):.5f}       {ch15.duct_outlet_mach(M1, lim['f_max'] * 1.001, 0, g)}")

    print("  differential forms at three Mach numbers (per unit df or dq; FD = finite difference of (15.46), step 1e-7):")
    e = 1e-7
    for M in (0.5, 0.9, 2.0):
        rf = ch15.duct_state_ratios(M, ch15.duct_outlet_mach(M, e, 0, g), e, 0, g)
        rq = ch15.duct_state_ratios(M, ch15.duct_outlet_mach(M, 0, e, g), 0, e, g)
        c = ch15.influence_coefficients(M, g, df=1.0)
        cq = ch15.influence_coefficients(M, g, dq=1.0)
        print(f"    M = {M}: −dp/p₁/df  FD {-(rf['p_ratio'] - 1) / e:+.5f} | corrected (15.47) {-c['dp_over_p']:+.5f} | printed {-ch15.influence_coefficients(M, g, df=1.0, printed_15_47=True)['dp_over_p']:+.5f}"
              f"   ||   u du/h₁/dq  FD {(g - 1) * M ** 2 * (rq['u_ratio'] - 1) / e:+.5f} | corrected (15.48) {cq['u_du_over_h']:+.5f} | printed "
              f"{ch15.influence_coefficients(M, g, dq=1.0, printed_15_48=True)['u_du_over_h']:+.5f}   ||   dT/T₁/dq  FD {(rq['T_ratio'] - 1) / e:+.5f} | (15.48) {cq['dT_over_T']:+.5f}")
    print(f"  heating cools for 1/√γ = {1 / np.sqrt(g):.4f} < M < 1: sign of dT/dq at M = 0.7, 0.9, 1.5: {ch15.rayleigh_temperature_sign([0.7, 0.9, 1.5], g)}")

    fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
    M1 = np.linspace(0.02, 4.0, 400 if args.fast else 1200)
    sub = M1 < 1.0  # the continuous root does not cross M = 1: draw the two sides as separate curves
    for fv, col in zip((0.0, 0.1, 0.2, 0.3, 0.4), (COLORS["ink"], C_SUB, C_TEMP, C_SONIC, C_SUP)):
        cu = ch15.duct_curve(M1, fv, 0.0, g)
        ax[0].plot(M1, np.where(sub, cu["M2_continuous"], np.nan), color=col, label=f"f = {fv}")
        ax[0].plot(M1, np.where(sub, np.nan, cu["M2_continuous"]), color=col)
        ax[0].plot(M1, cu["M2_shocked"], color=col, ls="--", lw=1)
    ax[0].fill_between([0, 1], 1, 4, color=C_REF, alpha=0.25)
    ax[0].text(0.08, 3.6, "ΔS < 0", fontsize=9)
    ax[0].set(xlabel="inlet Mach number M₁ [–]", ylabel="outlet Mach number M₂ [–]", xlim=(0, 4), ylim=(0, 4), title="friction (q = 0): solid = continuous root, dashed = other root")
    ax[0].legend(fontsize=8)
    for qv, col in zip((-0.4, -0.2, 0.0, 0.2, 0.4), (C_SUB, C_TEMP, COLORS["ink"], C_SONIC, C_SUP)):
        cu = ch15.duct_curve(M1, 0.0, qv, g)
        ax[1].plot(M1, np.where(sub, cu["M2_continuous"], np.nan), color=col, label=f"q = {qv:+.1f}")
        ax[1].plot(M1, np.where(sub, np.nan, cu["M2_continuous"]), color=col)
        ax[1].plot(M1, cu["M2_shocked"], color=col, ls="--", lw=1)
    ax[1].set(xlabel="inlet Mach number M₁ [–]", ylabel="M₂ [–]", xlim=(0, 4), ylim=(0, 4), title="heat (f = 0): heating drives M toward 1, cooling away")
    ax[1].legend(fontsize=8)
    Mm = np.linspace(0.2, 3.0, 400)
    fa, ra = ch15.fanno_line_Ts(Mm, g), ch15.rayleigh_line_Ts(Mm, g)
    ax[2].plot(fa["ds_over_cp"], fa["T_over_Tstar"], color=C_SUB, label="friction (Fanno) line — ours")
    ax[2].plot(ra["ds_over_cp"], ra["T_over_Tstar"], color=C_SUP, label="heating (Rayleigh) line — ours")
    mT = 1.0 / np.sqrt(g)
    ax[2].plot([ch15.rayleigh_line_Ts(mT, g)["ds_over_cp"]], [ch15.rayleigh_ratios(mT, g)["T"]], "o", color=C_ENT, label=rf"max T at $M=1/\sqrt{{\gamma}}$ = {mT:.3f}")
    ax[2].plot([0.0], [1.0], "o", color=C_SONIC, label="M = 1: largest entropy on both")
    ax[2].set(xlabel=r"$(s-s^*)/C_p$ [–]", ylabel=r"$T/T^*$ [–]", xlim=(-1.2, 0.1), title="T–s lines (star = sonic state of the same process)")
    ax[2].legend(fontsize=8)
    save(fig, out, "duct_friction_heat_charts")

    fl = ch15.fanno_length_to_f(0.4, 0.8, g)
    print(f"Fanno form (ours, not in the book): M 0.4 → 0.8 needs 4f_F L/D = {fl['friction_length']:.5f}; the book's lumped f for that duct = {fl['f']:.5f} → duct_outlet_mach returns {ch15.duct_outlet_mach(0.4, fl['f'], 0, g):.6f}")
    a, b = ch15.fanno_length_to_f(0.4, 0.6, g), ch15.fanno_length_to_f(0.6, 0.8, g)
    print(f"  two ducts in series 0.4 → 0.6 → 0.8: friction lengths add ({a['friction_length']:.5f} + {b['friction_length']:.5f} = {a['friction_length'] + b['friction_length']:.5f}); "
          f"their f do not ({a['f']:.5f} + {b['f']:.5f} = {a['f'] + b['f']:.5f} ≠ {fl['f']:.5f}) — each f refers to its own inlet pressure")
    T0r = ch15.rayleigh_ratios(0.6, g)["T0"] / ch15.rayleigh_ratios(0.3, g)["T0"]
    qq = (T0r - 1.0) * ch15.stagnation_temperature_ratio(0.3, g)
    print(f"Rayleigh form (ours): M 0.3 → 0.6 needs T₀₂/T₀₁ = {T0r:.5f}, i.e. q = Q/h₁ = {qq:.5f} → duct_outlet_mach returns {ch15.duct_outlet_mach(0.3, 0, qq, g):.6f}; "
          f"T/T* is largest ({ch15.rayleigh_ratios(mT, g)['T']:.5f}) at M = {mT:.4f}")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
