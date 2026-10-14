"""Chapter 15, §15.10 — thin airfoil theory in supersonic flow: the flat plate at incidence, wave drag, and the lift
slope against Mach number.

* Ackeret's linear theory (15.57)–(15.59) beside the exact shock–expansion result (ours) for a range of α;
* the four waves of Fig. 15.24 computed (fan and shock at each edge) and drawn;
* lift slope: 2π below the speed of sound (ch14), 4/√(M² − 1) above, a gap in between;
* a thin section of any shape (ours, Exercise 15.17): lift depends on incidence only, wave drag adds thickness.

Run: ``.venv/Scripts/python.exe scripts/ch15_supersonic_airfoil.py --no-show``   Figures -> outputs/ch15/airfoil_*.png
"""
from __future__ import annotations

import numpy as np
from ch15_common import C_ENT, C_PRES, C_REF, C_SONIC, C_SUB, C_SUP, COLORS, deg, finish, parse_args, save, setup

from fluidpy import ch15_compressible_flow as ch15


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch15.illustrative_inputs()
    g, Mi, al = I["gamma"], I["M_inf"], I["alpha"]
    lp = ch15.ackeret_flat_plate(al, Mi)
    print(f"flat plate at α = {deg(al):.0f}°, M∞ = {Mi}: Ackeret (15.59) C_L = 4α/√(M² − 1) = {lp['CL']:.5f}, C_D = 4α²/√(M² − 1) = {lp['CD']:.6f}, L/D = 1/α = {lp['L_over_D']:.3f}; "
          f"lift slope {lp['lift_slope']:.4f} per rad against 2π = {2 * np.pi:.4f}; (p − p∞)/p∞ = ±{ch15.ackeret_pressure(al, Mi, g):.5f} (15.57), C_p = ±{ch15.ackeret_cp(al, Mi):.5f}")
    print("   α [deg]   C_L linear  C_L exact   ratio    C_D linear  C_D exact   p_upper/p∞ (lin, exact)   p_lower/p∞ (lin, exact)")
    for a in (0.5, 1.0, 2.0, 5.0, 10.0, 15.0):
        l = ch15.flat_plate_supersonic(np.radians(a), Mi, g, method="linear")
        e = ch15.flat_plate_supersonic(np.radians(a), Mi, g, method="shock_expansion")
        print(f"   {a:5.1f}     {l['CL']:.5f}     {e['CL']:.5f}     {e['CL'] / l['CL']:.4f}   {l['CD']:.6f}    {e['CD']:.6f}    {l['p_upper']:.4f}, {e['p_upper']:.4f}            {l['p_lower']:.4f}, {e['p_lower']:.4f}")
    wp = ch15.flat_plate_wave_pattern(al, Mi, g)
    print(f"wave pattern at α = {deg(al):.0f}°: leading-edge fan from {deg(wp['le_fan_first']):.2f}° to {deg(wp['le_fan_last']):.2f}°, leading-edge shock at {deg(wp['le_shock']):.2f}°, "
          f"trailing-edge shock at {deg(wp['te_shock']):.2f}°, trailing-edge fan from {deg(wp['te_fan_first']):.2f}° to {deg(wp['te_fan_last']):.2f}° (from the free stream, + counter-clockwise); "
          f"M above {wp['M2']:.4f}, below {wp['M3']:.4f}; behind the trailing edge p/p∞ = {wp['p4_upper']:.5f} (upper), {wp['p4_lower']:.5f} (lower) — nearly, not exactly, matched")

    fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
    te = np.array([np.cos(al), -np.sin(al)])
    ax[0].plot([0, te[0]], [0, te[1]], color=COLORS["ink"], lw=3)
    L = 0.9
    for a_ in np.linspace(wp["le_fan_first"], wp["le_fan_last"], 4):
        ax[0].plot([0, L * np.cos(a_)], [0, L * np.sin(a_)], color=C_SUB, lw=0.9)
    ax[0].plot([0, L * np.cos(wp["le_shock"])], [0, L * np.sin(wp["le_shock"])], color=C_SUP, lw=2)
    ax[0].plot([te[0], te[0] + L * np.cos(wp["te_shock"])], [te[1], te[1] + L * np.sin(wp["te_shock"])], color=C_SUP, lw=2)
    for a_ in np.linspace(wp["te_fan_first"], wp["te_fan_last"], 4):
        ax[0].plot([te[0], te[0] + L * np.cos(a_)], [te[1], te[1] + L * np.sin(a_)], color=C_SUB, lw=0.9)
    ax[0].annotate("", xy=(-0.15, 0.3), xytext=(-0.6, 0.3), arrowprops=dict(arrowstyle="->", color=C_REF))
    ax[0].text(-0.6, 0.34, rf"$M_\infty$ = {Mi}", fontsize=9)
    ax[0].text(0.35, 0.12, rf"$p_2/p_\infty$ = {wp['p2']:.3f}", fontsize=8)
    ax[0].text(0.25, -0.3, rf"$p_3/p_\infty$ = {wp['p3']:.3f}", fontsize=8)
    ax[0].set(aspect="equal", xlabel="x / chord", ylabel="y / chord", xlim=(-0.7, 1.9), ylim=(-0.9, 0.8), title=f"flat plate at {deg(al):.0f}°: shocks (thick) and fans (thin)")

    aa = np.radians(np.linspace(0.25, 20.0, 80))
    cl_e = np.array([ch15.flat_plate_supersonic(a, Mi, g, method="shock_expansion")["CL"] for a in aa])
    cd_e = np.array([ch15.flat_plate_supersonic(a, Mi, g, method="shock_expansion")["CD"] for a in aa])
    lin = ch15.ackeret_flat_plate(aa, Mi)
    ax[1].plot(deg(aa), lin["CL"], color=C_SUB, ls="--", label=r"$C_L$ linear (15.59)")
    ax[1].plot(deg(aa), cl_e, color=C_SUB, label=r"$C_L$ shock–expansion (ours)")
    ax[1].plot(deg(aa), lin["CD"], color=C_SUP, ls="--", label=r"$C_D$ linear (15.59)")
    ax[1].plot(deg(aa), cd_e, color=C_SUP, label=r"$C_D$ shock–expansion (ours)")
    ax[1].set(xlabel="angle of attack α [deg]", ylabel="coefficient [–]", title=f"wave drag exists in inviscid flow; M∞ = {Mi}")
    ax[1].legend(fontsize=8)

    M = np.linspace(0.0, 3.5, 701)
    ax[2].plot(M, ch15.lift_slope_vs_mach(M), color=C_SUB, label=r"$2\pi$ (ch14) and $4/\sqrt{M^2-1}$ (15.59)")
    ax[2].plot(M, ch15.lift_slope_vs_mach(M, model="prandtl_glauert"), color=C_PRES, ls="--", lw=1, label=r"$2\pi/\sqrt{1-M^2}$ — not in the book")
    ax[2].axvspan(0.8, 1.2, color=C_REF, alpha=0.25, label="transonic: neither theory holds")
    ax[2].axvline(1.0, color=C_SONIC, lw=0.8, ls=":")
    ax[2].set(xlabel=r"$M_\infty$ [–]", ylabel=r"$dC_L/d\alpha$ [1/rad]", ylim=(0, 12), title="lift slope of a thin flat plate against Mach number")
    ax[2].legend(fontsize=8)
    save(fig, out, "airfoil_flat_plate_waves_slope")
    print(f"lift slope: M = 0.5: {ch15.lift_slope_vs_mach(0.5):.4f} (= 2π; flat-plate Zhukhovsky C_L/α at small α: {ch15.zhukhovsky_lift_coefficient(1e-4, 0.0) / 1e-4:.4f}); "
          f"M = 1.0: {ch15.lift_slope_vs_mach(1.0)}; M = 1.5: {ch15.lift_slope_vs_mach(1.5):.4f}; M = 2: {ch15.lift_slope_vs_mach(2.0):.4f}; M = 3: {ch15.lift_slope_vs_mach(3.0):.4f}; "
          f"slopes equal 2π at M = {np.sqrt(1 + (2 / np.pi) ** 2):.4f}")

    tau = I["thickness_ratio"]
    sections = {"flat plate": (lambda x: 0 * x, lambda x: 0 * x),
                "biconvex (parabolic arcs)": (lambda x: 2 * tau * x * (1 - x), lambda x: -2 * tau * x * (1 - x)),
                "half-sine thickness": (lambda x: 0.5 * tau * np.sin(np.pi * x), lambda x: -0.5 * tau * np.sin(np.pi * x)),
                "double wedge": (lambda x: tau * np.minimum(x, 1 - x), lambda x: -tau * np.minimum(x, 1 - x)),
                "cambered plate": (lambda x: 0.5 * tau * np.sin(np.pi * x), lambda x: 0.5 * tau * np.sin(np.pi * x))}
    b = np.sqrt(Mi ** 2 - 1)
    print(f"thin sections (ours; thickness or camber parameter τ = {tau}) at α = {deg(al):.0f}°, M∞ = {Mi}:   C_L        C_D        = incidence + camber + thickness")
    for name, (yu, yl) in sections.items():
        s = ch15.ackeret_airfoil(yu, yl, al, Mi, 801 if args.fast else 4001)
        print(f"   {name:26} {s['CL']:.6f}   {s['CD']:.6f}   {s['CD_incidence']:.6f} + {s['CD_camber']:.6f} + {s['CD_thickness']:.6f}")
    print(f"   closed forms for the thickness part: biconvex (16/3)τ²/β = {16 / 3 * tau ** 2 / b:.6f}, half-sine (π²/2)τ²/β = {np.pi ** 2 / 2 * tau ** 2 / b:.6f}, double wedge 4τ²/β = {4 * tau ** 2 / b:.6f}")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
