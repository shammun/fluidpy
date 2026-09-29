"""§8.3 (Example 8.2, Fig. 8.10): Hele-Shaw flow past a disk — gap-averaged streamlines equal 2-D ideal flow round a
cylinder, the parabolic profile across the gap, φ = −z(h − z)p/2μ, and a finite-difference Laplace solve of the
gap-averaged stream function that converges to the ideal one (V3 route).

Run: ``.venv/Scripts/python.exe scripts/ch08_hele_shaw.py --no-show [--fast]``
Figure → outputs/ch08/c06_hele_shaw.png.
"""
from __future__ import annotations

import numpy as np

from ch08_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch08_laminar_flow as ch08
from fluidpy.core import potential as PF  # noqa: F401


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    a, h, U, mu = 0.01, 1e-3, 1e-3, 1e-3  # 1 cm disk, 1 mm gap, 1 mm/s, water
    xs = np.linspace(-4 * a, 4 * a, 241)
    X, Y = np.meshgrid(xs, xs)
    f = ch08.hele_shaw_cylinder(X, Y, h / 2, U, a, h, mu)
    Re_h = U * h / 1e-6
    print(f"gap Reynolds number Uh/ν = {Re_h:.3g}, (h/a)² Re = {(h / a) ** 2 * U * a / 1e-6:.3g} (lubrication valid)")
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.8))
    cs = ax[0].contourf(X / a, Y / a, f["p"], 30, cmap="RdBu_r")
    fig.colorbar(cs, ax=ax[0], label="p − p₀ [Pa]")
    ax[0].contour(X / a, Y / a, f["psi_mean"], np.linspace(-3 * U * a, 3 * U * a, 21), colors="k", linewidths=0.7)
    ax[0].add_patch(plt.Circle((0, 0), 1, color=COLORS["muted"]))
    ax[0].set_aspect("equal")
    ax[0].set_title("gap-averaged streamlines and pressure")
    ax[0].set_xlabel("x/a")
    ax[0].set_ylabel("y/a")
    z = np.linspace(0, h, 50)
    for (xp, yp), c in (((-3 * a, 0.0), COLORS["accent"]), ((0.0, 1.05 * a), COLORS["orange"])):
        fp = ch08.hele_shaw_cylinder(xp, yp, z, U, a, h, mu)
        ax[1].plot(fp["u"] * 1e3, z * 1e3, color=c, lw=2, label=f"(x, y) = ({xp / a:.2f}a, {yp / a:.2f}a)")
        print(f"at ({xp/a:.2f}a, {yp/a:.2f}a): gap mean u = {np.trapezoid(fp['u'], z) / h * 1e3:.5f} mm/s, "
              f"ideal-flow value {float(ch08.hele_shaw_cylinder(xp, yp, 0.0, U, a, h, mu)['u_mean']) * 1e3:.5f} mm/s")
    ax[1].set_xlabel("u [mm/s]")
    ax[1].set_ylabel("z [mm]")
    ax[1].legend(fontsize=8)
    ax[1].set_title("parabolic profile across the gap")
    ns = (33, 65, 129) if args.fast else (33, 65, 129, 257)
    errs, hs = [], []
    for n in ns:
        g = ch08.hele_shaw_streamfunction_grid(n, 1.0, 1.0, 4.0)
        errs.append(g["err_far"])
        hs.append(g["h"])
        print(f"grid n = {n}: max |ψ − ψ_ideal| (r ≥ 2a) = {g['err_far']:.3e}")
    order = np.polyfit(np.log(hs), np.log(errs), 1)[0]
    ax[2].loglog(hs, errs, "o-", color=COLORS["accent"], label=f"observed order {order:.2f}")
    ax[2].set_xlabel("grid spacing h/a")
    ax[2].set_ylabel("max error r ≥ 2a (ψ/Ua)")
    ax[2].legend(fontsize=8)
    ax[2].set_title("finite-difference Laplace solve → ideal flow")
    print(f"observed order of the staircase-disk Laplace solve: {order:.3f}")
    save(fig, out, "c06_hele_shaw")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
