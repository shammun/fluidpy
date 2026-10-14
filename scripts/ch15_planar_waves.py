"""Chapter 15, §15.9 — pressure waves in planar flow: the Mach cone, oblique shocks, the weak-wave rule and the
Prandtl–Meyer expansion.

* wave fronts of a moving source for U < c and U > c; Mach angle (15.49); the delay of the boom (ours);
* the deflection–shock-angle chart (15.53): weak and strong roots, maximum deflection and sonic point (ours),
  detachment; a wedge;
* the weak-wave rule (15.55) against the exact shock and the exact expansion;
* the Prandtl–Meyer function (15.56), the fan at a corner and its staircase of small turns; Mach lines of a curved
  wall.

Run: ``.venv/Scripts/python.exe scripts/ch15_planar_waves.py --no-show``   Figures -> outputs/ch15/waves_*.png
"""
from __future__ import annotations

import numpy as np
from ch15_common import C_ENT, C_PRES, C_REF, C_SONIC, C_SUB, C_SUP, C_TEMP, COLORS, deg, finish, parse_args, save, setup
from scipy.integrate import quad

from fluidpy import ch15_compressible_flow as ch15


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch15.illustrative_inputs()
    g = I["gamma"]
    c = ch15.PerfectGas(g, I["R"]).c(I["T_air"])

    # --- Mach cone ---------------------------------------------------------------------------------------------------------
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
    th = np.linspace(0, 2 * np.pi, 241)
    for a_, M in zip(ax[:2], (0.6, 2.0)):
        w = ch15.source_wavefronts(M * c, c, 6.0, 5, 1.0)
        for xc, r in zip(w["x_centres"], w["radii"]):
            a_.plot((xc + r * np.cos(th)) / 1e3, r * np.sin(th) / 1e3, color=C_SUB, lw=1)
            a_.plot([xc / 1e3], [0], ".", color=C_REF)
        a_.plot([w["x_source"] / 1e3], [0], "o", color=C_SUP)
        if M > 1:
            L = 1.7 * M * c * 1.0 * 5 / 1e3
            for sgn in (+1, -1):
                a_.plot([w["x_source"] / 1e3, w["x_source"] / 1e3 - L * np.cos(w["mu"])], [0, sgn * L * np.sin(w["mu"])], color=C_ENT, lw=1.4)
        a_.set(aspect="equal", xlabel="x [km]", ylabel="y [km]", title=f"source at M = {M}: " + ("fronts run ahead" if M < 1 else f"Mach cone, μ = {deg(w['mu']):.1f}°"))
        print(f"source at M = {M}: the newest front reaches {w['ahead'][0]:+.0f} m ahead of the source ({'ahead' if w['ahead'][0] > 0 else 'none ahead: zone of silence'}); "
              + (f"Mach angle {deg(ch15.mach_angle(M)):.3f}°" if M >= 1 else f"mach_angle({M}) = {ch15.mach_angle(M)} (no cone below M = 1)"))
    Mf, h = I["M_flight"], I["h_flight"]
    print(f"flight at M = {Mf}, height {h / 1e3:.0f} km (c = {c:.1f} m/s): μ = {deg(ch15.mach_angle(Mf)):.2f}°, boom arrives {ch15.mach_cone_delay(h, Mf, c):.2f} s after the aircraft is overhead "
          f"(it is then {Mf * c * ch15.mach_cone_delay(h, Mf, c) / 1e3:.2f} km further on); in zone of action at (x, y) = (2h, h): {ch15.in_zone_of_action(2 * h, h, Mf)}, at (h, h): {ch15.in_zone_of_action(h, h, Mf)}")

    # --- oblique shock chart ---------------------------------------------------------------------------------------------
    for M1, col in zip((1.5, 2.0, 3.0, 5.0), (C_SUB, C_TEMP, C_SUP, C_ENT)):
        cu = ch15.oblique_shock_curve(M1, g, 200 if args.fast else 500)
        ax[2].plot(deg(cu["sigma"]), deg(cu["delta"]), color=col, label=f"M₁ = {M1}: δ_max = {deg(cu['delta_max']):.2f}°")
        ax[2].plot([deg(cu["sigma_max"])], [deg(cu["delta_max"])], "o", color=col, ms=4)
        ax[2].plot([deg(cu["sigma_sonic"])], [deg(cu["delta_sonic"])], "s", color=col, ms=3, mfc="none")
    ax[2].set(xlabel="shock angle σ [deg]", ylabel="deflection δ [deg]", title="(15.53): ● maximum deflection, □ M₂ = 1 (ours)", xlim=(0, 90), ylim=(0, 45))
    ax[2].legend(fontsize=8)
    save(fig, out, "waves_mach_cone_oblique_chart")
    print("  M₁     δ_max [deg]  σ at max [deg]  δ_sonic [deg]  σ_sonic [deg]")
    for M1 in (1.2, 1.5, 2.0, 3.0, 5.0, 1e4):
        md, sd = ch15.max_deflection(M1, g), ch15.sonic_deflection(M1, g)
        print(f"  {M1:<6g} {deg(md['delta_max']):9.4f}    {deg(md['sigma_max']):9.4f}      {deg(sd['delta_sonic']):9.4f}     {deg(sd['sigma_sonic']):9.4f}")
    print(f"  limit δ_max(M₁ → ∞) = asin(1/γ) = {deg(np.arcsin(1 / g)):.4f}°")
    Mw, dw = I["M_wedge"], I["wedge_half_angle"]
    for br in ("weak", "strong"):
        o = ch15.oblique_shock(Mw, delta=dw, branch=br, gamma=g)
        cp_ = ch15.oblique_shock_components(Mw, o["sigma"], g)
        print(f"wedge M₁ = {Mw}, δ = {deg(dw):.0f}°, {br:6} root: σ = {deg(o['sigma']):.4f}°, M_n1 = {o['Mn1']:.4f}, M₂ = {o['M2']:.4f}, p₂/p₁ = {o['p_ratio']:.4f}, ρ₂/ρ₁ = {o['rho_ratio']:.4f} "
              f"= tanσ/tan(σ − δ) = {np.tan(o['sigma']) / np.tan(o['sigma'] - o['delta']):.4f}, p₀₂/p₀₁ = {o['p0_ratio']:.4f}; tangential speed unchanged: v/V₁ = {cp_['v_over_V1']:.4f}")
    wf = ch15.wedge_flow(Mw, dw, g)
    print(f"  wedge_flow: attached = {wf['attached']}; detaches below M₁ = {wf['M_detach']:.4f}; at M₁ = 1.4 the same wedge: attached = {ch15.wedge_flow(1.4, dw, g)['attached']}; "
          f"round trip δ(σ) = {deg(ch15.deflection_angle(Mw, wf['sigma'], g)):.6f}°")
    try:
        ch15.oblique_shock(2.0, delta=np.radians(30.0), gamma=g)
    except ValueError as e:
        print(f"  oblique_shock(2, δ = 30°) raises: {str(e)[:70]}…")

    # --- weak-wave rule -------------------------------------------------------------------------------------------------------
    print("  weak-wave rule (15.55) at M = 2:  δ [deg]   linear     exact      exact/linear")
    for d in (0.5, 1.0, 2.0, 5.0, 10.0, 20.0, -5.0, -10.0):
        w = ch15.weak_wave_error(2.0, np.radians(d), g)
        print(f"                                    {d:+5.1f}    {w['linear']:+.5f}   {w['exact']:+.5f}   {w['ratio']:.4f}  ({w['kind']})")

    # --- Prandtl–Meyer ------------------------------------------------------------------------------------------------------
    for M in (1.5, 2.0, 3.0):
        qv = quad(lambda m: np.sqrt(m * m - 1) / (1 + 0.5 * (g - 1) * m * m) / m, 1.0, M)[0]
        print(f"Prandtl–Meyer ν({M}) = {deg(ch15.prandtl_meyer(M, g)):.4f}° (closed form (15.56)) = {deg(qv):.4f}° (quadrature of the integrand); inverse → {ch15.mach_from_prandtl_meyer(ch15.prandtl_meyer(M, g), g):.6f}")
    print(f"  ν_max = {deg(ch15.prandtl_meyer_max(g)):.3f}°")
    Mfan, turn = I["M_fan"], I["fan_turn"]
    pm = ch15.prandtl_meyer_turn(Mfan, turn, g)
    print(f"fan: M₁ = {Mfan} turned by {deg(turn):.0f}° (expansion; the book's δ = {deg(pm['delta_book']):.0f}°): M₂ = {pm['M2']:.5f}, p₂/p₁ = {pm['p_ratio']:.5f}, T₂/T₁ = {pm['T_ratio']:.5f}; "
          f"first ray at μ₁ = {deg(pm['mu1']):.2f}°, last at μ₂ = {deg(pm['mu2']):.2f}° to the new wall; two half-turns give M₂ = "
          f"{ch15.prandtl_meyer_turn(ch15.prandtl_meyer_turn(Mfan, turn / 2, g)['M2'], turn / 2, g)['M2']:.5f}")
    for N in (5, 50, 500):
        st = ch15.fan_staircase(Mfan, turn, N, g)
        print(f"  staircase of {N:3d} weak turns: M₂ = {st['M'][-1]:.5f} (error {st['M_error']:+.2e}), p₂/p₁ = {st['p_over_p1'][-1]:.5f} by (15.29) (error {st['p_error']:+.2e}), {st['p_over_p1_weak'][-1]:.5f} by repeated (15.55) (error {st['p_weak_error']:+.2e})")

    fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
    fan = ch15.expansion_fan(Mfan, turn, g, n_rays=9, r_start=0.6)
    ax[0].plot([-1.0, 0.0], [0.0, 0.0], color=COLORS["ink"], lw=2)
    ax[0].plot([0.0, 1.6 * np.cos(turn)], [0.0, -1.6 * np.sin(turn)], color=COLORS["ink"], lw=2)
    for ang in fan["ray_angle"]:
        ax[0].plot([0, 1.6 * np.cos(ang)], [0, 1.6 * np.sin(ang)], color=C_SUB, lw=0.8)
    fine = ch15.expansion_fan(Mfan, turn, g, n_rays=60, r_start=0.6)
    y0 = fine["streamline_y"][0]
    xs = np.concatenate([[-1.0], fine["streamline_x"], [fine["streamline_x"][-1] + 0.6 * np.cos(turn)]])
    ys = np.concatenate([[y0], fine["streamline_y"], [fine["streamline_y"][-1] - 0.6 * np.sin(turn)]])
    ax[0].plot(xs, ys, color=C_SUP, lw=1.4, label="a streamline")
    ax[0].set(aspect="equal", xlabel="x [m]", ylabel="y [m]", xlim=(-1, 1.7), title=f"Prandtl–Meyer fan: M {Mfan} → {fan['M2']:.3f} in a {deg(turn):.0f}° turn")
    ax[0].legend(fontsize=8)
    dd = np.radians(np.linspace(-20, 20, 161))
    lin = ch15.weak_wave_pressure(2.0, dd, g)
    ex = np.array([ch15.weak_wave_error(2.0, d, g)["exact"] if d != 0 else 0.0 for d in dd])
    ax[1].plot(deg(dd), ex, color=C_PRES, label="exact: shock (δ > 0), expansion (δ < 0)")
    ax[1].plot(deg(dd), lin, color=COLORS["ink"], ls="--", lw=1, label="weak-wave rule (15.55)")
    ax[1].axvline(0, color=C_REF, lw=0.8)
    ax[1].set(xlabel="deflection δ [deg] (+ into the flow)", ylabel=r"$(p_2-p_1)/p_1$ [–]", title="how small is small? M₁ = 2")
    ax[1].legend(fontsize=8)
    xw = np.linspace(0.0, 1.0, 13)
    for sgn, col, lab in ((+1, C_SUP, "compression: lines converge"), (-1, C_SUB, "expansion: lines diverge")):
        cw = ch15.curved_wall_mach_lines(2.0, lambda s, k=sgn: k * np.radians(12.0) * np.clip(s, 0, 1) ** 2, xw, g)
        off = 0.0 if sgn > 0 else -1.2
        ax[2].plot(cw["x_wall"], cw["y_wall"] + off, color=COLORS["ink"], lw=2)
        for x0, y0_, an in zip(cw["x_wall"], cw["y_wall"], cw["line_angle"]):
            ax[2].plot([x0, x0 + 2.2 * np.cos(an)], [y0_ + off, y0_ + off + 2.2 * np.sin(an)], color=col, lw=0.7)
        if cw["first_crossing"] is not None:
            ax[2].plot([cw["first_crossing"][0]], [cw["first_crossing"][1] + off], "o", color=C_ENT)
            print(f"curved wall (compression, 12° over 1 m, M₁ = 2): neighbouring Mach lines first cross at (x, y) = ({cw['first_crossing'][0]:.3f}, {cw['first_crossing'][1]:.3f}) m — a shock starts there; "
                  f"Mach number at the end of the wall {cw['M'][-1]:.4f}")
        else:
            print(f"curved wall (expansion): no crossing; Mach number at the end of the wall {cw['M'][-1]:.4f}")
        ax[2].text(0.02, off + 0.08, lab, fontsize=8, color=col)
    ax[2].set(aspect="equal", xlabel="x [m]", ylabel="y [m]", xlim=(0, 3), ylim=(-1.5, 1.6), title="Mach lines of a curved wall (simple wave)")
    save(fig, out, "waves_fan_weakrule_curvedwall")
    r1, r2 = ch15.weak_wave_sympy(), ch15.prandtl_meyer_sympy()
    print(f"sympy: (15.54) from (15.50), (15.53): residual {r1['residual_15_54']}; limit slope {r1['slope']} → (15.55) residual {r1['residual_15_55']}; "
          f"d/dM of (15.56) minus the integrand: {r2['derivative_residual']} (numeric {r2['residual_numeric']:.1e}); tangential momentum: v₂ − v₁ = {ch15.tangential_momentum_sympy()['residual']}")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
