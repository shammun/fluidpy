"""§10.4 (C09): the steady convection–diffusion layer (10.84)–(10.89), centred-scheme wiggles for R_cell > 2 (10.90)–(10.92),
upwind (10.93) and its numerical diffusivity 0.5 R_cell D (10.94); a stretched grid; the overflow of the printed (10.86).

Run: ``.venv/Scripts/python.exe scripts/ch10_steady_cd.py --no-show``   Figure -> outputs/ch10/c09_cell_peclet.png.
"""
from __future__ import annotations

import numpy as np

from ch10_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch10_computational_fluid_dynamics as ch10


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    xf = np.linspace(0, 1, 801)
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
    for ai, (n, R) in enumerate(((20, 20.0), (20, 80.0))):
        ax[ai].plot(xf, ch10.steady_cd_exact(xf, R), color=COLORS["muted"], label="exact (10.86)")
        c = ch10.steady_cd_fd(n, R, "central")
        u = ch10.steady_cd_fd(n, R, "upwind")
        s = ch10.steady_cd_fd(n, R, "central", grid="stretched", beta_s=2.5)
        ax[ai].plot(c["x"], c["T"], "o-", color=COLORS["orange"], ms=4, label="centred (10.91)")
        ax[ai].plot(u["x"], u["T"], "s-", color=COLORS["teal"], ms=4, label="upwind (10.93)")
        ax[ai].plot(s["x"], s["T"], "^-", color=COLORS["blue"], ms=4, label="centred, stretched grid")
        Rc = R / n
        ax[ai].set_title(f"R = {R:g}, n = {n}: R_cell = {Rc:g}")
        ax[ai].set_xlabel("x / L")
        ax[ai].set_ylabel("T")
        ax[ai].legend(fontsize=7)
        w = ch10.wiggle_indicator(c["T"])
        print(f"R = {R:g}, n = {n} (R_cell = {Rc:g}): centred root r = {ch10.discrete_root(Rc):.4g}, wiggles {w['wiggles']} "
              f"(min T = {w['min']:.4f}); upwind r = {ch10.discrete_root(Rc, 'upwind'):.4g}; stretched-grid min T = "
              f"{s['T'].min():.2e}")
        d = ch10.steady_cd_discrete_exact(n, R)
        print(f"   tridiagonal vs closed form T_j = (r^j − 1)/(r^n − 1): max diff {np.abs(d['T'] - c['T']).max():.2e}")
    Rc = np.linspace(0, 6, 601)
    ax[2].plot(Rc, [ch10.discrete_root(r) if abs(r - 2) > 1e-3 else np.nan for r in Rc], color=COLORS["orange"], label="centred r")
    ax[2].plot(Rc, [ch10.discrete_root(r, "upwind") for r in Rc], color=COLORS["teal"], label="upwind r = 1 + R_cell")
    ax[2].axhline(0, color=COLORS["rose"], lw=1)
    ax[2].axvline(2, color=COLORS["muted"], ls=":")
    ax[2].set_ylim(-8, 8)
    ax[2].set_xlabel("R_cell = uΔx/D")
    ax[2].set_ylabel("root r of T_j = r^j")
    ax[2].set_title("r < 0 (alternating signs) ⇔ R_cell > 2")
    ax[2].legend(fontsize=8)
    for R in (10.0, 100.0, 1000.0):
        print(f"layer thickness at e^-1 (R = {R:g}): {ch10.cd_layer_thickness(R):.5f} L (1/R = {1 / R:.5f}); "
              f"at e^-2: {ch10.cd_layer_thickness(R, level=np.exp(-2)):.5f} L")
    with np.errstate(all="ignore"):
        print(f"R = 1000 at x = 0.999: printed form {ch10.steady_cd_exact(0.999, 1000.0, printed=True)}, "
              f"scaled form {ch10.steady_cd_exact(0.999, 1000.0):.6f}")
    print(f"numerical diffusivity of steady upwind, u = 0.1 m/s, Δx = 1e4 m: "
          f"{ch10.numerical_diffusivity(0.1, 1e4, scheme='upwind_steady'):.4g} m²/s")
    for s, p in (("central", 2), ("upwind", 1)):
        errs, hs = [], []
        for n in (20, 40, 80, 160):
            r = ch10.steady_cd_fd(n, 10.0, s)
            errs.append(np.max(np.abs(r["T"] - ch10.steady_cd_exact(r["x"], 10.0))))
            hs.append(1 / n)
        print(f"{s} order against (10.86) at R = 10: {ch10.observed_order(hs, errs):.3f} (expected {p})")
    # upwind vs the exact solution of the modified equation (10.94): the gap shrinks with R_cell
    for n in (5, 10, 20, 40):
        r = ch10.steady_cd_fd(n, 20.0, "upwind")
        Rm = 20.0 / (1 + 0.5 * 20.0 / n)
        print(f"   n = {n:3d}: max |upwind − modified-equation exact| = "
              f"{np.max(np.abs(r['T'] - ch10.steady_cd_exact(r['x'], Rm))):.4f}")
    save(fig, out, "c09_cell_peclet")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
