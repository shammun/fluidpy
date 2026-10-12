"""Chapter 14, §14.2 and §14.6 — lifting-line theory: downwash (14.13), the sine series (14.18)–(14.22), the solver
of (14.20), and the elliptic loading (14.23)–(14.27).

* planforms: span, area, mean chord, aspect ratio (14.1)–(14.2);
* downwash by direct quadrature of the principal-value integral (14.13) against the series (14.19, last member) —
  and the printed middle member of (14.19), which has the opposite sign (slip #1); Glauert's integral;
* lift and induced drag: (14.21), (14.22) against quadrature of (14.15) and the wake kinetic energy;
* the collocation solver: exact for the elliptic planform, convergence table for a rectangular wing, taper sweep,
  a rectangular wing twisted to an elliptic loading.

Run: ``.venv/Scripts/python.exe scripts/ch14_lifting_line.py --no-show``   Figures -> outputs/ch14/lifting_line_*.png
"""
from __future__ import annotations

import numpy as np
from ch14_common import C_CIRC, C_DRAG, C_FLOW, C_LIFT, C_REF, finish, parse_args, save, setup

from fluidpy import ch14_aerodynamics as ch14


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch14.illustrative_inputs()
    s, A, U, rho, K = I["s"], I["A"], I["U"], I["rho"], I["K"]

    # --- planforms ---------------------------------------------------------------------------------------------------
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
    plans = {"rectangular": ch14.planform("rectangular", s, A / s), "tapered": ch14.planform("tapered", s, 2 * A / (s * (1 + I["taper"])), I["taper"]),
             "elliptic": ch14.planform("elliptic", s, 4 * A / (np.pi * s))}
    print(f"three planforms of the same span s = {s} m and area A = {A} m²:")
    for (name, P), col in zip(plans.items(), (C_FLOW, C_CIRC, C_LIFT)):
        ax[0].plot(P["y"], 0.25 * P["c"], color=col, label=f"{name}, Λ = {P['AR']:.2f}")
        ax[0].plot(P["y"], -0.75 * P["c"], color=col)
        print(f"  {name:12s}: root chord {P['c_root']:.4f} m, A = {P['A']:.4f} m² (trapezoid on 201 points {P['A_trapezoid']:.4f}), "
              f"c̄ = A/s = {P['c_mean']:.4f} m, Λ = s²/A = {P['AR']:.4f} = s/c̄ = {s / P['c_mean']:.4f}")
    ax[0].set(aspect="equal", xlabel="y [m] (span)", ylabel="x [m] (chordwise)", title="planforms (quarter-chord line straight)")
    ax[0].legend(fontsize=8)

    # --- downwash: quadrature of (14.13) against the series ------------------------------------------------------------
    Gn = np.array([11.0, 0.0, 1.4, 0.0, -0.5])
    g1 = np.array([0.4, 1.0, np.pi / 2, 2.3])
    y1 = ch14.span_position(g1, s)
    wq = ch14.downwash(y1, lambda y: ch14.circulation_series(Gn, ch14.span_angle(y, s)), s)
    ws = ch14.downwash_series(Gn, g1, s)
    wp = ch14.downwash_series(Gn, g1, s, printed_middle=True)
    print(f"downwash of the loading Γ_n = {Gn} m²/s at γ₁ = {np.round(g1, 3)}:")
    print(f"  quadrature of (14.13)        : {np.round(wq, 8)} m/s")
    print(f"  series (14.19), last member  : {np.round(ws, 8)} m/s   (max difference {np.max(np.abs(wq - ws)):.1e})")
    print(f"  (14.19) middle member AS PRINTED: {np.round(wp, 8)} m/s   ← opposite sign: the denominator must be cos γ − cos γ₁")
    G1 = 9.0
    we = ch14.downwash(y1, lambda y: ch14.elliptic_circulation(y, s, G1), s)
    print(f"  elliptic loading Γ₁ = {G1}: quadrature {np.round(we, 9)} — uniform, Γ₁/2s = {ch14.elliptic_downwash(G1, s):.9f} (14.24)")
    par = ch14.downwash(0.0, lambda y: G1 * (1.0 - (2.0 * y / s) ** 2), s)
    print(f"  parabolic loading Γ₀(1 − (2y/s)²) with Γ₀ = {G1}: mid-span downwash {par:.6f} m/s = {par / ch14.elliptic_downwash(G1, s):.5f} × the elliptic value "
          f"(4/π = {4 / np.pi:.5f})")
    print("  Glauert's integral PV∫cos nγ dγ/(cos γ − cos γ₁) at γ₁ = 0.9: " + ", ".join(
        f"n={n}: {ch14.glauert_integral(n, 0.9, 'quad'):+.8f} / {ch14.glauert_integral(n, 0.9):+.8f}" for n in range(0, 5)))

    # --- forces: series against quadrature and wake energy ---------------------------------------------------------
    ga = (np.arange(4000) + 0.5) * np.pi / 4000
    Lq, Dq = ch14.lift_and_induced_drag(ch14.span_position(ga, s), ch14.circulation_series(Gn, ga), ch14.downwash_series(Gn, ga, s), rho, U)
    print(f"forces of that loading: L = {ch14.lift_from_series(Gn, s, rho, U):.4f} N (14.21) vs {Lq:.4f} N (trapezoid on (14.15)); "
          f"D_i = {ch14.induced_drag_from_series(Gn, rho):.6f} N (14.22) vs {Dq:.6f} N (14.15) vs wake energy per length "
          f"{ch14.wake_kinetic_energy_per_length(Gn, s, rho):.6f} J/m; span efficiency {ch14.span_efficiency(Gn):.5f}")

    # --- the solver ---------------------------------------------------------------------------------------------------
    alpha, beta = np.deg2rad(4.0), 0.02
    Pe = plans["elliptic"]
    r = ch14.solve_lifting_line(Pe, s, U, alpha, K=K, beta=beta, N=8, rho=rho)
    ew = ch14.elliptic_wing(alpha, Pe["AR"], K=K, beta=beta)
    print(f"elliptic planform, N = 8: Γ_n/Γ₁ for n > 1: max {np.abs(r['Gn'][1:]).max() / r['Gn'][0]:.1e}; C_L = {r['CL']:.8f} "
          f"(closed form K(α + β)/(1 + K/πΛ) = {ew['CL']:.8f}); C_Di = {r['CDi']:.8f} = C_L²/(πΛ) = {ch14.induced_drag_coefficient(r['CL'], Pe['AR']):.8f}; "
          f"slope {r['slope']:.5f} per rad (section {K:.4f}); w/U = {r['w'][0] / U:.6f} uniform")
    print(f"  with K = 2π the slope factor is 1/(1 + 2/Λ): {ch14.finite_wing_lift_slope(2 * np.pi, Pe['AR']) / (2 * np.pi):.6f} = {1 / (1 + 2 / Pe['AR']):.6f}")
    print("rectangular wing Λ = 6, K = 2π, α = 1 rad (convergence of the collocation; symmetric modes):")
    print("     N | span efficiency | slope dC_L/dα | condition number | max residual of (14.17) between stations")
    gm = np.linspace(0.05, np.pi - 0.05, 400)
    for N in (2, 4, 8, 16, 32, 64):
        q = ch14.solve_lifting_line(1.0, 6.0, 1.0, 1.0, K=2 * np.pi, N=N, symmetric=True)
        res = np.max(np.abs(ch14.lifting_line_residual(q["Gn"], gm, 1.0, 6.0, 1.0, 1.0, K=2 * np.pi)))
        print(f"    {N:2d} | {q['e_span']:.7f}       | {q['slope']:.7f}     | {q['cond']:9.1f}        | {res:.2e}")

    Pr = plans["rectangular"]
    rr = ch14.solve_lifting_line(Pr, s, U, alpha, K=K, beta=beta, N=40, symmetric=True, rho=rho)
    gg = np.linspace(0.0, np.pi, 301)
    yy = ch14.span_position(gg, s)
    for (name, P), col in zip(plans.items(), (C_FLOW, C_CIRC, C_LIFT)):
        q = ch14.solve_lifting_line(P, s, U, alpha, K=K, beta=beta, N=40, symmetric=True, rho=rho)
        ax[1].plot(yy, ch14.circulation_series(q["Gn"], gg), color=col, label=f"{name}: e_span = {q['e_span']:.4f}")
        ax[2].plot(yy[1:-1], ch14.downwash_series(q["Gn"], gg[1:-1], s), color=col, label=f"{name}: C_L = {q['CL']:.4f}")
        print(f"our wing as a {name:11s} planform (α = 4°, β = {beta}, K = {K:.4f}): C_L = {q['CL']:.5f}, C_Di = {q['CDi']:.6f}, e_span = {q['e_span']:.5f}, "
              f"L = {q['L']:.1f} N, D_i = {q['Di']:.2f} N, Γ₃/Γ₁ = {q['Gn'][2] / q['Gn'][0]:+.4f}")
    ax[1].set(xlabel="y [m]", ylabel="Γ(y) [m²/s]", title="bound circulation from the solver (same span, area, α)")
    ax[1].legend(fontsize=8)
    ax[2].set(xlabel="y [m]", ylabel="downwash w(y) [m/s]", title="downwash: uniform only for the elliptic loading", ylim=(0, 3.0 * rr["w"][-1]))
    ax[2].legend(fontsize=8)
    save(fig, out, "lifting_line_planforms_loading")

    # --- taper sweep and twist ---------------------------------------------------------------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.4))
    tap = np.linspace(0.0, 1.0, 41)
    for AR, col in zip((5.0, 9.0, 15.0), (C_REF, C_LIFT, C_DRAG)):
        es = [ch14.solve_lifting_line(ch14.planform("tapered", AR, 2.0 / (1.0 + t), t), AR, 1.0, 1.0, K=K, N=40, symmetric=True)["e_span"] for t in tap]
        k = int(np.argmax(es))
        ax[0].plot(tap, es, color=col, label=f"Λ = {AR:.0f}: best taper {tap[k]:.2f}, e_span = {es[k]:.4f}")
        print(f"taper sweep Λ = {AR:4.1f}: e_span = {es[0]:.4f} (pointed tips), max {es[k]:.4f} at taper {tap[k]:.3f}, {es[-1]:.4f} (rectangle)")
    ax[0].set(xlabel="taper ratio (tip chord / root chord) [–]", ylabel=r"span efficiency $e_{span}$ [–]", title="a tapered wing can be nearly elliptic")
    ax[0].legend(fontsize=8)
    c, ar_ = Pr["c_root"], np.deg2rad(3.0)
    T = ch14.elliptic_twist_rectangular(ar_, c, s, K=K)
    q = ch14.solve_lifting_line(c, s, U, T["alpha_t"], K=K, twist=lambda y: ar_ * np.sqrt(np.clip(1 - (2 * y / s) ** 2, 0, None)), N=12, rho=rho)
    print(f"rectangular wing with elliptic twist (root twist 3°): fly at α_t = Kcα_r/(4s) = {np.rad2deg(T['alpha_t']):.4f}°; the solver returns "
          f"Γ_n/Γ₁ (n > 1) ≤ {np.abs(q['Gn'][1:]).max() / q['Gn'][0]:.1e}, e_span = {q['e_span']:.10f}, C_L = {q['CL']:.6f} (closed form {T['CL']:.6f}), "
          f"C_Di = {q['CDi']:.7f} = C_L²/(πΛ) = {T['CDi']:.7f}")
    tw = ch14.twist_for_elliptic_loading(c, yy[1:-1], s, K=K, CL_target=0.5)
    ax[1].plot(yy[1:-1], np.rad2deg(tw["alpha_local"]), color=C_LIFT, label="rectangular planform")
    Pt = plans["tapered"]
    tw2 = ch14.twist_for_elliptic_loading(Pt, yy[1:-1], s, K=K, CL_target=0.5)
    ax[1].plot(yy[1:-1], np.rad2deg(tw2["alpha_local"]), color=C_CIRC, label=f"tapered {I['taper']}")
    ax[1].axhline(np.rad2deg(ch14.elliptic_wing(0.0, Pe["AR"], K=K)["CL"] * 0 + 0.5 / ch14.finite_wing_lift_slope(K, Pe["AR"])), color=C_REF, ls="--", label="elliptic planform (no twist)")
    ax[1].set(xlabel="y [m]", ylabel="local absolute angle of attack [deg]", title="twist that makes the loading elliptic at C_L = 0.5")
    ax[1].legend(fontsize=8)
    save(fig, out, "lifting_line_taper_twist")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
