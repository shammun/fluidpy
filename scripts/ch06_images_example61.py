"""§6.3 (C08): the method of images — vortex images flip sign, source images keep it (Figs. 6.14–6.16), the two-source
flow (6.41)/(6.53) with its implicit streamlines, the circle theorem vs ch05's circle images, and Example 6.1: a vortex
beside a wall drifting at Γ/4πh and the wall-pressure signal (suction, then over-pressure) by the closed form and by
an independent numeric route (point-vortex ODE + ∂φ/∂t by central differences + unsteady Bernoulli).

Run: ``.venv/Scripts/python.exe scripts/ch06_images_example61.py --no-show``
Figures → outputs/ch06/c08_images_example61.png.
"""
from __future__ import annotations

import numpy as np

from ch06_drawings import flow_net, parse_args, save, setup

from fluidpy import ch06_ideal_flow as ch06
from fluidpy.core import biot_savart as BS
from fluidpy.core.style import COLORS


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    x = np.linspace(-3, 3, 13)
    vs = ch06.Flow(ch06.mirror([ch06.Vortex(1.0, (0.4, 0.8))], "y=0"))
    ss = ch06.Flow(ch06.mirror([ch06.Source(1.0, (0.4, 0.8))], "y=0"))
    print(f"wall y = 0: max |v| with vortex images {np.max(np.abs(vs.velocity(x, 0 * x)[1])):.1e}, "
          f"with source images {np.max(np.abs(ss.velocity(x, 0 * x)[1])):.1e} m/s")
    bad = ch06.Flow([ch06.Vortex(1.0, (0.4, 0.8)), ch06.Vortex(1.0, (0.4, -0.8))])
    print(f"  wrong variant (same-sign vortex image): max |v| on the wall {np.max(np.abs(bad.velocity(x, 0 * x)[1])):.3f}")
    m, a = 1.0, 1.0
    ts = ch06.two_sources(m, a)
    xs = np.array([1.5, 2.0, 3.0])
    ys = ch06.two_source_streamline(0.3, m, a, xs)
    print(f"two sources (6.41): ψ on the implicit curve for ψ = 0.3: {np.round(ts.psi(xs, ys), 12)}")
    circ = ch06.circle_theorem(ch06.Flow([ch06.Vortex(1.3, (2.0, 0.5))]), 1.0)
    P, G = BS.circle_image_system(np.array([[2.0], [0.5]]), [1.3], 1.0, inside=False)
    pt = np.array([[1.4, -0.3], [0.2, 2.1]])
    print(f"circle theorem vs ch05 circle_image_system: max |Δu| "
          f"{np.max(np.abs(np.array(circ.velocity(pt[0], pt[1])) - BS.point_vortex_velocity(pt, P, G))):.1e}")

    Gm, h, rho = 1.0, 1.0, 1000.0
    kt = ch06.example_6_1_times(Gm, h, rho)
    t = np.linspace(0, 60, 61 if args.fast else 241)
    cl = ch06.example_6_1(t, Gm, h, rho, route="closed")
    tn = np.array([0.0, 5.0, kt["t_zero"], kt["t_max"], 40.0])
    nu = ch06.example_6_1(tn, Gm, h, rho, route="numeric")
    cn = ch06.example_6_1(tn, Gm, h, rho, route="closed")
    print(f"Example 6.1 (Γ = {Gm} m²/s, h = {h} m, ρ = {rho}): drift Γ/4πh = {Gm / (4 * np.pi * h):.6f} m/s")
    print(f"  p(0,0,0) − p∞ = {kt['p_min']:.6f} Pa; zero at t = 4πh²/Γ = {kt['t_zero']:.4f} s; "
          f"max {kt['p_max']:.6f} Pa at t = 4√3πh²/Γ = {kt['t_max']:.4f} s")
    for ti, pn, pc, yn in zip(tn, nu["p_origin"], cn["p_origin"], nu["xi"][1]):
        print(f"  t = {ti:7.3f} s: p − p∞ closed {pc:+.9f}, numeric {pn:+.9f} Pa; ξ_y numeric {yn:.9f} "
              f"(closed {Gm * ti / (4 * np.pi * h):.9f}) m")

    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
    n = 121 if args.fast else 181
    flow_net(ax[0], ts, (-3, 3), (-3, 3), n=n, psi_levels=np.linspace(-0.49, 0.49, 29),
             stag=ts.stagnation_points(box=(-2, 2, -2, 2)), title="two sources ±a (6.41): x = 0 is a wall")
    ax[0].plot(xs, ys, "s", color=COLORS["amber"], ms=5, label="(6.41) implicit curve, ψ = 0.3")
    ax[0].legend(fontsize=8, loc="lower right")
    ax[1].plot(t, cl["p_origin"], color=COLORS["accent"], label="closed form")
    ax[1].plot(tn, nu["p_origin"], "o", color=COLORS["orange"], label="numeric route")
    ax[1].axhline(0, color=COLORS["muted"], lw=0.8)
    ax[1].set_xlabel("t [s]")
    ax[1].set_ylabel("p(0, 0, t) − p∞ [Pa]")
    ax[1].set_title("Example 6.1: wall pressure at the origin", fontsize=10)
    ax[1].legend(fontsize=8)
    yw = np.linspace(-6, 10, 400)
    for ti, col in ((0.0, COLORS["blue"]), (kt["t_zero"], COLORS["teal"]), (kt["t_max"], COLORS["orange"])):
        wp = ch06.example_6_1_wall_pressure(yw, ti, Gm, h, rho)
        ax[2].plot(wp, yw, color=col, label=f"t = {ti:.1f} s")
    ax[2].axvline(0, color=COLORS["muted"], lw=0.8)
    ax[2].set_xlabel("p − p∞ on the wall x = 0 [Pa]")
    ax[2].set_ylabel("y [m]")
    ax[2].set_title("pressure along the wall (our extension)", fontsize=10)
    ax[2].legend(fontsize=8)
    save(fig, out, "c08_images_example61")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
