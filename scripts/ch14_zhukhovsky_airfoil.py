"""Chapter 14, §14.4 — airfoil shapes from the Zhukhovsky transformation z = ζ + b²/ζ (14.3).

* four circles and their images: flat plate, circular arc, symmetric section, cambered airfoil (Figs. 14.11–14.14);
* the circular arc checked against its circle x² + (y + 2b cot 2β)² = (2b csc 2β)², numerically and with sympy;
* the shape rules (14.8) c ≃ 4b, camber ≃ ½βc, t_max/c ≃ (3√3/4)e against the measured contour;
* the cusp at the trailing edge; streamlines round our cambered section at incidence.

Run: ``.venv/Scripts/python.exe scripts/ch14_zhukhovsky_airfoil.py --no-show``   Figures -> outputs/ch14/zhukhovsky_*.png
"""
from __future__ import annotations

import numpy as np
from ch14_common import C_CIRC, C_FLOW, C_LIFT, C_REF, COLORS, finish, parse_args, save, setup, streamlines

from fluidpy import ch14_aerodynamics as ch14


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch14.illustrative_inputs()
    b, e, beta, alpha = I["b"], I["e"], I["beta"], I["alpha"]

    # --- four circles and their images -------------------------------------------------------------------------------
    cases = [("flat plate (e = 0, β = 0)", 0.0, 0.0), (f"circular arc (e = 0, β = {2 * beta:.2f})", 0.0, 2 * beta),
             (f"symmetric (e = {e}, β = 0)", e, 0.0), (f"cambered (e = {e}, β = {2 * beta:.2f})", e, 2 * beta)]
    fig, ax = plt.subplots(2, 4, figsize=(16, 6.4))
    tt = np.linspace(0, 2 * np.pi, 400)
    for k, (name, ee, bb) in enumerate(cases):
        A = ch14.zhukhovsky_airfoil(b, ee, bb, 800)
        mu, a = A["mu"], A["a"]
        ax[0, k].plot(b * np.cos(tt), b * np.sin(tt), color=C_REF, lw=1, ls="--", label="|ζ| = b")
        ax[0, k].plot(A["zeta"].real, A["zeta"].imag, color=C_FLOW, label="circle")
        ax[0, k].plot([mu.real], [mu.imag], "o", color=C_CIRC, ms=5)
        ax[0, k].plot([b, -b], [0, 0], "x", color=COLORS["rose"], ms=7)
        ax[0, k].set(aspect="equal", xlabel="ξ [m]", ylabel="η [m]", title=f"ζ-plane: {name}", xlim=(-1.6 * b, 1.6 * b), ylim=(-1.5 * b, 1.6 * b))
        ax[1, k].plot(np.append(A["x"], A["x"][0]), np.append(A["y"], A["y"][0]), color=C_LIFT)
        ax[1, k].set(aspect="equal", xlabel="x [m]", ylabel="y [m]", title="z-plane", xlim=(-2.6 * b, 2.6 * b), ylim=(-1.2 * b, 1.2 * b))
        sh = ch14.zhukhovsky_shape(b, ee, bb)
        print(f"{name:34s}: centre μ = {mu.real:+.4f}{mu.imag:+.4f}i, a = {a:.4f}; chord {sh['chord']:.5f}, camber "
              f"{sh['camber']:.5f}, t_max {sh['t_max']:.5f} at x = {sh['x_tmax']:+.4f}")
    ax[0, 0].legend(fontsize=7, loc="upper right")
    save(fig, out, "zhukhovsky_four_circles")

    # --- the circular arc --------------------------------------------------------------------------------------------
    bb = 2 * beta
    arc = ch14.circular_arc(b, bb)
    A = ch14.zhukhovsky_airfoil(b, 0.0, bb, 2000)
    res = np.max(np.abs(np.abs(A["z"] - arc["centre"]) - arc["radius"]))
    k = int(np.argmax(A["y"]))
    slope = np.arctan2(A["y"][3] - A["y"][0], -(A["x"][3] - A["x"][0]))
    print(f"circular arc β = {bb}: all 2000 mapped points lie on the circle of centre {arc['centre'].imag:+.5f}i and radius "
          f"{arc['radius']:.5f} to {res:.1e} m; height {A['y'][k]:.6f} = 2b tan β = {arc['height']:.6f}; end tangent "
          f"{slope:.5f} rad ≈ 2β = {arc['end_angle']:.5f}")
    sym = ch14.circular_arc_sym()
    print(f"sympy: residual of (14.5) = {sym['res_14_5']}, of (14.6) = {sym['res_14_6']}, of the arc's circle = "
          f"{sym['res_arc']}  → all zero: {sym['all_zero']}")

    # --- shape rules (14.8) ------------------------------------------------------------------------------------------
    print(f"shape rules (14.8), first-order thickness coefficient 3√3/4 = {ch14.THICKNESS_RULE_COEFFICIENT:.6f} (printed 1.3):")
    es = np.array([0.001, 0.01, 0.03, 0.06, 0.09, 0.12, 0.18, 0.25])
    meas = [ch14.zhukhovsky_shape(b, q, 0.0) for q in es]
    for q, m in zip(es, meas):
        print(f"  e = {q:5.3f}: chord/4b = {m['chord'] / (4 * b):.5f}, t_max/c = {m['thickness_ratio']:.5f}, (t_max/c)/e = "
              f"{m['thickness_ratio'] / q:.4f}, rule/measured − 1 = {ch14.THICKNESS_RULE_COEFFICIENT * q / m['thickness_ratio'] - 1:+.1%}, "
              f"x_tmax = {m['x_tmax']:+.4f} (rule −b)")
    betas = np.array([0.02, 0.06, 0.12, 0.2])
    for q in betas:
        m = ch14.zhukhovsky_shape(b, e, q)
        print(f"  β = {q:4.2f} (e = {e}): camber/c = {m['camber_ratio']:.5f}, rule ½β = {0.5 * q:.5f} ({0.5 * q / m['camber_ratio'] - 1:+.1%}), "
              f"max camber at {m['x_camber_over_c']:.3f} c")
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    ee = np.linspace(0.005, 0.3, 40)
    tr = np.array([ch14.zhukhovsky_shape(b, q, 0.0)["thickness_ratio"] for q in ee])
    ax[0].plot(ee, tr, color=C_LIFT, label="measured on the contour")
    ax[0].plot(ee, ch14.THICKNESS_RULE_COEFFICIENT * ee, color=C_REF, ls="--", label=r"rule (14.8): $(3\sqrt{3}/4)\,e$")
    ax[0].set(xlabel="thickness parameter e [–]", ylabel=r"$t_{\max}/c$ [–]", title="thickness rule against the contour")
    ax[0].legend(fontsize=8)
    bbs = np.linspace(0.0, 0.3, 31)
    cr = np.array([ch14.zhukhovsky_shape(b, e, q)["camber_ratio"] for q in bbs])
    ax[1].plot(bbs, cr, color=C_LIFT, label=f"measured (e = {e})")
    ax[1].plot(bbs, 0.5 * bbs, color=C_REF, ls="--", label=r"rule (14.8): $\frac{1}{2}\beta$")
    ax[1].set(xlabel="camber angle β [rad]", ylabel="camber / chord [–]", title="camber rule against the contour")
    ax[1].legend(fontsize=8)
    for q, col in zip((0.03, e, 0.18), (C_REF, C_LIFT, COLORS["rose"])):
        A = ch14.zhukhovsky_airfoil(b, q, beta, 600)
        ax[2].plot(np.append(A["x"], A["x"][0]), np.append(A["y"], A["y"][0]), color=col, label=f"e = {q}")
    ax[2].axvline(-b, color=C_REF, lw=0.8, ls=":")
    ax[2].set(aspect="equal", xlabel="x [m]", ylabel="y [m]", title=f"sections for β = {np.rad2deg(beta):.0f}°: thickest near x = −b")
    ax[2].legend(fontsize=8)
    save(fig, out, "zhukhovsky_shape_rules")

    # --- the cusp ----------------------------------------------------------------------------------------------------
    print("trailing-edge angle measured from the two neighbouring contour points (a cusp: → 0):")
    for n in (100, 400, 1600, 6400):
        A = ch14.zhukhovsky_airfoil(b, e, beta, n)
        print(f"  n = {n:5d}: {np.rad2deg(ch14.trailing_edge_angle(A['x'], A['y'], 0)):.4f}°")
    for nn in (1.0, 1.1, 1.5):
        w = ch14.trailing_edge_wedge(nn)
        print(f"  corner flow n = {nn}: included angle {np.rad2deg(w['included_angle']):6.2f}°, speed at r = 0.01 d: "
              f"{w['speed'](0.01):.4f} U, C_p at the edge = {w['Cp_edge']}")

    # --- streamlines round our cambered section at incidence -----------------------------------------------------------
    Us = I["U_section"]
    F = ch14.zhukhovsky_flow(Us, b, e, beta, alpha)
    fig, ax = plt.subplots(figsize=(9, 5))
    streamlines(ax, F, (-3.6 * b, 3.6 * b), (-2.2 * b, 2.2 * b))
    st = ch14.stagnation_points(Us, b, e, beta, alpha)
    ax.plot([st["z_front"].real], [st["z_front"].imag], "o", color=COLORS["rose"], ms=6, zorder=6, label="front stagnation point")
    ax.set(xlabel="x [m]", ylabel="y [m]", title=f"Zhukhovsky airfoil e = {e}, β = {np.rad2deg(beta):.0f}°, α = {np.rad2deg(alpha):.0f}°, Kutta circulation")
    ax.legend(fontsize=8, loc="lower left")
    save(fig, out, "zhukhovsky_streamlines")
    print(f"flow: Γ_K = {F.Gamma_kutta:.5f} m²/s (clockwise), C_L exact = {ch14.zhukhovsky_lift_coefficient(alpha, beta, e, True):.4f}, "
          f"2π(α + β) = {ch14.zhukhovsky_lift_coefficient(alpha, beta):.4f}; ψ on the surface varies by "
          f"{np.ptp(F.stream(*[q[1:] for q in (ch14.zhukhovsky_airfoil(b, e, beta, 400)['x'], ch14.zhukhovsky_airfoil(b, e, beta, 400)['y'])])):.1e}")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
