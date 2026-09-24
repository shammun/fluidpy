"""§6.3 (C03, C05): the half-body (6.30)–(6.32) — stream + source, stagnation point a = m/2πU, body ψ = m/2, width
h_max = m/2U, the surface C_p (our reduction) and its zero near 113°; the no-through-flow check (6.16), the far-field
decay (6.17) (1/R for an open body vs 1/R² for a closed one) and the Rankine oval (Exercise 6.19) as the closed twin.

Run: ``.venv/Scripts/python.exe scripts/ch06_half_body.py --no-show``
Figures → outputs/ch06/c05_half_body.png.
"""
from __future__ import annotations

import numpy as np

from ch06_drawings import flow_net, parse_args, save, setup

from fluidpy import ch06_ideal_flow as ch06
from fluidpy.core.style import COLORS


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    U, m = 1.0, 2 * np.pi
    hb = ch06.half_body(U, m)
    num = ch06.half_body_numbers(U, m)
    stag = hb.stagnation_points(box=(-3, 3, -2, 2))
    print(f"half-body U = {U} m/s, m = 2π m²/s: a = {num['a']:.6f} m (Newton: {stag[0].real:.6f}), "
          f"ψ_body = {num['psi_body']:.6f} m²/s, h_max = {num['h_max']:.6f} m, 2Uh_max = {2 * U * num['h_max']:.6f} = m")
    th = np.linspace(0.05, 2 * np.pi - 0.05, 400)
    xb, yb = ch06.half_body_shape(U, m, th)
    print(f"  max |ψ − m/2| on the body: {np.max(np.abs(np.asarray(hb.psi(xb, yb)) - m / 2)):.2e};  "
          f"max |u·n| on the body: {ch06.normal_velocity_on(hb, xb + 1j * yb, closed=False):.2e} m/s")
    cp_flow = np.asarray(hb.cp(xb, yb))
    print(f"  surface C_p: Flow.cp vs reduced formula max diff {np.max(np.abs(cp_flow - ch06.half_body_surface_cp(th))):.2e}")
    t0 = ch06.half_body_cp_zero_angle()
    print(f"  C_p = 0 on the body at θ = {np.degrees(t0):.4f}° (root of sin θ + 2(π − θ)cos θ = 0)")
    R = np.array([10.0, 20.0, 40.0, 80.0])
    ro = ch06.rankine_oval(U, m, 1.0)
    ff_h = [ch06.far_field_check(hb, r) for r in R]
    ff_o = [ch06.far_field_check(ro["flow"], r) for r in R]
    print("  far-field |u − U| slopes (log–log): half-body "
          f"{np.polyfit(np.log(R), np.log(ff_h), 1)[0]:.3f}, Rankine oval {np.polyfit(np.log(R), np.log(ff_o), 1)[0]:.3f}")
    print(f"Rankine oval (source/sink ±1 m): half-length {ro['half_length']:.6f} m, half-width {ro['half_width']:.6f} m, "
          f"stagnation {np.round(ro['flow'].stagnation_points(box=(-4, 4, -3, 3)), 6)}")

    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
    n = 121 if args.fast else 201
    lv = np.linspace(-6, 6, 41)
    flow_net(ax[0], hb, (-3, 6), (-3.5, 3.5), n=n, psi_levels=lv, stag=stag,
             body=xb + 1j * yb, title="half-body: stream + source (6.31), body ψ = m/2")
    ax[0].axhline(num["h_max"], color=COLORS["muted"], ls=":", lw=1)
    thu = np.linspace(0.05, np.pi - 1e-6, 300)
    ax[1].plot(np.degrees(thu), ch06.half_body_surface_cp(thu), color=COLORS["accent"])
    ax[1].axhline(0, color=COLORS["muted"], lw=0.8)
    ax[1].axvline(np.degrees(t0), color=COLORS["orange"], ls="--", label=f"C_p = 0 at {np.degrees(t0):.1f}°")
    ax[1].invert_xaxis()
    ax[1].set_xlabel("θ on the body, from +x at the source [deg] (nose = 180°)")
    ax[1].set_ylabel("C_p = 1 − |u|²/U²  (6.32)")
    ax[1].set_title("half-body surface pressure: excess at the nose, deficit beyond", fontsize=10)
    ax[1].legend(fontsize=8)
    ob = np.linspace(0, 2 * np.pi, 400, endpoint=False)
    # oval outline from the ψ = 0 contour: take ψ on rays and root-find the radius
    from scipy.optimize import brentq
    rr = []
    for t in ob:
        f = lambda r: float(ro["flow"].elements[0].psi(r * np.cos(t), abs(r * np.sin(t)))  # noqa: E731
                            + sum(e.psi(r * np.cos(t), abs(r * np.sin(t))) for e in ro["flow"].elements[1:]))
        rr.append(brentq(f, 0.3, 5.0) if f(0.3) * f(5.0) < 0 else np.nan)
    rr = np.array(rr)
    flow_net(ax[2], ro["flow"], (-4, 4), (-3, 3), n=n, psi_levels=np.linspace(-4, 4, 33),
             stag=ro["flow"].stagnation_points(box=(-4, 4, -3, 3)), body=rr * np.exp(1j * ob),
             title="Rankine oval: Σm = 0 closes the body")
    save(fig, out, "c05_half_body")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
