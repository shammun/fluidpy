"""§8.6 (C15, Fig. 8.20): Oseen's correction (8.53) — fluid-frame streamlines at Re = 0.2, 1, 2 next to Stokes, the
wake asymmetry, the Re → 0 limit to (8.48) in sympy, and the velocity from analytic differentiation checked against
finite differences of ψ.

Run: ``.venv/Scripts/python.exe scripts/ch08_oseen.py --no-show [--fast]``
Figure → outputs/ch08/c15_oseen.png.
"""
from __future__ import annotations

import numpy as np

from ch08_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch08_laminar_flow as ch08
from fluidpy.core.potential import axisym_velocity_spherical


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    a, U = 1.0, 1.0
    n = 100 if args.fast else 200
    xs = np.linspace(-10, 10, n)
    X, Y = np.meshgrid(xs, xs)
    R = np.hypot(X, Y)
    TH = np.arctan2(np.abs(Y), X)
    sgn = np.where(Y >= 0, 1.0, -1.0)
    fig, ax = plt.subplots(1, 4, figsize=(18, 4.6))
    lev = np.array([0.1, 0.3, 0.6, 1.0, 1.5, 2.2, 3.0, 4.0])
    for axi, Re in zip(ax, (0.0, 0.2, 1.0, 2.0)):
        psi = sgn * ch08.oseen_streamfunction(R, TH, U, a, Re, "fluid")
        axi.contour(X, Y, psi, np.concatenate([-lev[::-1], lev]), colors=COLORS["accent"], linewidths=0.8)
        axi.add_patch(plt.Circle((0, 0), a, color=COLORS["muted"]))
        if Re > 0:
            axi.add_patch(plt.Circle((0, 0), 2 * a / Re, fill=False, ls=":", color=COLORS["orange"]))
        axi.set_aspect("equal")
        axi.set_title("Stokes (Re = 0)" if Re == 0 else f"Oseen, Re = 2aU/ν = {Re}")
        axi.set_xlabel("x/a (sphere moving to −x)")
        up = abs(float(ch08.oseen_streamfunction(5.0, 2.6, U, a, Re, "fluid")))
        dn = abs(float(ch08.oseen_streamfunction(5.0, np.pi - 2.6, U, a, Re, "fluid")))
        cd = (f"; C_D Oseen {float(ch08.oseen_drag_coefficient(Re)):.4g} vs Stokes "
              f"{float(ch08.stokes_drag_coefficient(Re)):.4g}") if Re > 0 else ""
        print(f"Re = {Re}: |ψ_fluid| at r = 5a, θ = 149° (upstream) {up:.4f} vs θ = 31° (downstream) {dn:.4f}{cd}")
    ax[0].set_ylabel("y/a")
    save(fig, out, "c15_oseen")
    d = ch08.oseen_limit_sympy()
    print(f"sympy: lim Re→0 ψ_Oseen − ψ_Stokes = {d['residual']}; first-order term {d['first_order']}")
    fr = lambda rr, tt: ch08.oseen_streamfunction(rr, tt, U, a, 1.0)  # noqa: E731
    r_t, th_t = np.array([1.5, 3.0, 8.0]), np.array([0.5, 1.5, 2.5])
    ua = np.array(ch08.oseen_velocity(r_t, th_t, U, a, 1.0))
    un = np.array(axisym_velocity_spherical(fr, r_t, th_t, h=1e-5))
    print(f"oseen_velocity vs finite differences of ψ: max |Δ| = {np.max(np.abs(ua - un)):.2e}")
    ur0, ut0 = ch08.oseen_velocity(a, th_t, U, a, 0.5)
    print(f"no slip only to O(Re): at r = a, Re = 0.5: |u_r| ≤ {np.max(np.abs(ur0)):.3f} U, |u_θ| ≤ "
          f"{np.max(np.abs(ut0)):.3f} U")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
