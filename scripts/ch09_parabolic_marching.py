"""§9.1 (N14, N24): the boundary-layer equations are parabolic — x plays the role of time.  Our von Mises marching solver reproduces Blasius,
forgets a wrong inlet profile downstream on a flat plate, and stops at separation for a linearly retarded flow.

Run: ``.venv/Scripts/python.exe scripts/ch09_parabolic_marching.py --no-show``   Figure -> outputs/ch09/n14_parabolic_marching.png.
"""
from __future__ import annotations

import numpy as np

from ch09_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch09_boundary_layers as ch09


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    U, nu = 1.0, 1e-5
    flat = ch09.outer_flow("flat", U=U)
    nx = 150 if args.fast else 300
    x = np.geomspace(1e-3, 1.0, nx)
    r = ch09.march_boundary_layer(flat, x, nu, ny=400)
    exact = ch09.blasius_wall_shear(r["x"], U, 1.0, nu)
    err = r["tau0"][1:] / exact[1:] - 1.0
    print(f"marching vs Blasius: tau0 error at x = 1 m: {err[-1]:+.2e} (BDF2 in x: second order, falls ~4x when nx doubles until the sigma-grid error ~3e-5 is reached)")
    # wrong inlet: exponential profile at x0 = 0.01
    x2 = np.geomspace(1e-2, 1.0, nx)
    d0 = 0.5 * np.sqrt(nu * 1e-2 / U)
    yin = np.linspace(0, 30 * np.sqrt(nu * 1e-2), 800)
    r2 = ch09.march_boundary_layer(flat, x2, nu, u_inlet=(yin, U * (1 - np.exp(-yin / d0))), ny=400)
    e2 = r2["tau0"] / ch09.blasius_wall_shear(r2["x"], U, 1.0, nu) - 1.0
    print(f"memory of a wrong inlet profile: tau0 error {e2[1]:+.2f} right after the inlet, {e2[-1]:+.3f} at x = 1 m")
    of = ch09.outer_flow("linear_retarded", U0=1.0, L=1.0)
    xr = np.linspace(0.01, 0.99, 400)
    rr = ch09.march_boundary_layer(of, xr, 1e-4, ny=400)
    th = ch09.thwaites(np.concatenate([[0.0], xr]), of, 1e-4)
    print(f"retarded flow U_e = U0(1 - x/L): marching separates at x/L = {rr['x_sep']:.4f};  Thwaites (lambda = -0.09): {th['x_sep']:.4f}; "
          f"(l = 0): {th['x_sep_l0']:.4f}")
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].loglog(r["x"][1:], np.abs(err) + 1e-12, color=COLORS["accent"], label="start at the similarity profile")
    ax[0].loglog(r2["x"][1:], np.abs(e2[1:]) + 1e-12, color=COLORS["orange"], label="start at a wrong profile")
    ax[0].set_xlabel("x  [m]")
    ax[0].set_ylabel("|tau0 / tau0,Blasius - 1|")
    ax[0].legend(fontsize=8)
    ax[0].set_title("the inlet profile is forgotten downstream")
    ax[1].plot(rr["x"], rr["tau0"], color=COLORS["teal"], label="marching")
    ax[1].plot(th["x"], th["tau0"], color=COLORS["rose"], ls="--", label="Thwaites (exact-FS closure)")
    ax[1].axhline(0, color=COLORS["muted"], lw=0.8)
    ax[1].set_xlabel("x / L")
    ax[1].set_ylabel("tau0 / rho  [m^2/s^2]")
    ax[1].legend(fontsize=8)
    ax[1].set_title("linearly retarded flow: tau0 falls to zero (separation)")
    save(fig, out, "n14_parabolic_marching")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
