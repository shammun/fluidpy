"""§6.3 (C06): the circular cylinder (6.33)–(6.35) and d'Alembert's paradox — streamlines, surface C_p = 1 − 4 sin²θ
against the angle from the front stagnation point, the surface-pressure force (zero drag and lift), parity with ch03's
``cylinder_flow``, and the moving cylinder as a doublet (Fig. 6.11). The "real flow" curve is a labelled qualitative
sketch only (no cited measurement is used).

Run: ``.venv/Scripts/python.exe scripts/ch06_cylinder_cp.py --no-show``
Figures → outputs/ch06/c06_cylinder_cp.png.
"""
from __future__ import annotations

import numpy as np

from ch06_drawings import circle, flow_net, parse_args, pressure_arrows, save, separated_cp_band, setup

from fluidpy import ch03_kinematics as ch03
from fluidpy import ch06_ideal_flow as ch06
from fluidpy.core.style import COLORS


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    U, a, rho = 1.0, 1.0, 1.225
    fl = ch06.cylinder(U, a)
    x = np.array([1.3, -2.0, 0.4, 3.0])
    y = np.array([0.4, 1.1, -1.5, -0.2])
    dv = np.max(np.abs(np.array(fl.velocity(x, y)) - np.array(ch03.cylinder_flow(x, y, U, a))))
    print(f"cylinder U = {U} m/s, a = {a} m: parity with ch03.cylinder_flow max |Δu| = {dv:.2e} m/s")
    th = np.linspace(0, 2 * np.pi, 361)
    cp = ch06.cylinder_surface_cp(th, U, a)
    print(f"  C_p range on the body: max {cp.max():+.6f} (θ = 0, π), min {cp.min():+.6f} (θ = ±90°)")
    print(f"  max |u·n| on r = a: {ch06.normal_velocity_on(fl, circle(a, 256)):.2e} m/s")
    for n in (8, 16, 64):
        F = ch06.surface_pressure_force(fl, R=a, rho=rho, n=n)
        print(f"  pressure force (n = {n:3d} nodes): D = {F.D:+.2e} N/m, L = {F.L:+.2e} N/m (d'Alembert: both 0)")
    print(f"  Blasius on R = 2a: {ch06.blasius_force(fl, R=2 * a, rho=rho)}")

    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
    n = 121 if args.fast else 201
    flow_net(ax[0], fl, (-3, 3), (-2, 2), n=n, psi_levels=np.linspace(-2, 2, 21), body=circle(a),
             stag=fl.stagnation_points(box=(-2, 2, -2, 2)), title="cylinder = stream + doublet d = −2πUa² e_x (6.33)")
    cz = circle(a, 36)
    pressure_arrows(ax[0], cz, 0.5 * rho * U ** 2 * np.asarray(ch06.cylinder_surface_cp(np.angle(cz), U, a)),
                    cz / a, scale=0.25)
    beta = np.linspace(0, np.pi, 181)  # angle from the upstream stagnation point = π − θ
    ax[1].plot(np.degrees(beta), ch06.cylinder_surface_cp(np.pi - beta, U, a), color=COLORS["accent"],
               label="ideal: 1 − 4 sin²θ (6.35)")
    lo, hi = separated_cp_band(np.degrees(beta))
    ax[1].fill_between(np.degrees(beta), lo, hi, color=COLORS["muted"], alpha=0.3,
                       label="real flow — qualitative sketch band, not data (separation ≈ 80°)")
    ax[1].set_xlabel("angle from the front stagnation point, π − θ [deg]")
    ax[1].set_ylabel("C_p")
    ax[1].legend(fontsize=8)
    ax[1].set_title("surface pressure: fore–aft symmetric ⇒ no drag", fontsize=10)
    xg = np.linspace(-3, 3, n)
    yg = np.linspace(-2, 2, n)
    X, Y = np.meshgrid(xg, yg, indexing="xy")
    moving = ch06.Flow([ch06.Doublet((-2 * np.pi * U * a ** 2, 0.0))], inside=fl.inside)
    ax[2].contour(X, Y, np.asarray(moving.psi(X, Y)), levels=np.linspace(-1, 1, 21), colors=COLORS["ink"],
                  linewidths=0.8, negative_linestyles="solid")
    c = circle(a)
    ax[2].fill(c.real, c.imag, color="#d9dce4")
    ax[2].set_aspect("equal")
    ax[2].set_title("cylinder moving through still fluid: a doublet (Fig. 6.11)", fontsize=10)
    ax[2].set_xlabel("x [m]")
    save(fig, out, "c06_cylinder_cp")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
