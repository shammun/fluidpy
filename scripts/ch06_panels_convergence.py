"""Source-panel method (our labelled extension of the chapter's singularity-distribution idea; method of Hess & Smith
1967): constant-strength source panels on a circle (exact to round-off at the control points for every N — a symmetry
of the regular polygon) and on an ellipse (algebraic convergence against the Zhukhovsky-mapped exact surface speed),
with Σλ_jS_j = 0 for the closed body.

Run: ``.venv/Scripts/python.exe scripts/ch06_panels_convergence.py --no-show``
Figures → outputs/ch06/n90_source_panels.png.
"""
from __future__ import annotations

import numpy as np

from ch06_drawings import parse_args, save, setup

from fluidpy import ch06_ideal_flow as ch06
from fluidpy.core.style import COLORS


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    Ns = [8, 16, 32, 64, 128] if not args.fast else [8, 16, 32, 64]
    ec = [ch06.panel_cp_error(N, "circle") for N in Ns]
    ee = [ch06.panel_cp_error(N, "ellipse") for N in Ns]
    print("circle, max |C_p − (1 − 4 sin²θ)| at control points:", [f"{e:.1e}" for e in ec])
    print("ellipse A = 1, B = 0.5, max |ΔC_p|:", [f"{e:.2e}" for e in ee],
          f"→ observed order {-np.polyfit(np.log(Ns[-3:]), np.log(ee[-3:]), 1)[0]:.3f}")
    A, B = 1.0, 0.5
    nu = np.linspace(0, 2 * np.pi, 64, endpoint=False) + np.pi / 64
    r = ch06.source_panels(A * np.cos(nu), B * np.sin(nu), 1.0)
    print(f"  N = 64: Σλ_jS_j = {r['net_source']:.1e} m²/s (closed body), cond = {r['cond']:.2f}")
    u, v = ch06.panel_velocity(np.array([3.0, 0.0]), np.array([0.0, 2.0]), r)
    print(f"  velocity off the body at (3, 0), (0, 2): u = {np.round(u, 6)}, v = {np.round(v, 6)} m/s")

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.2))
    nuc = np.arctan2(r["ym"] / B, r["xm"] / A)
    order = np.argsort(nuc)
    ax[0].plot(np.degrees(nuc[order]), r["cp"][order], "o", color=COLORS["teal"], ms=4, label="64 source panels")
    nn = np.linspace(-np.pi, np.pi, 400)
    ax[0].plot(np.degrees(nn), 1 - np.asarray(ch06.ellipse_surface_speed(nn, 1.0, A, B)) ** 2, color=COLORS["ink"],
               label="exact (mapped flow)")
    ax[0].set_xlabel("ellipse parameter ν [deg]")
    ax[0].set_ylabel("C_p")
    ax[0].legend(fontsize=8)
    ax[0].set_title("ellipse A = 1, B = 0.5 in a stream along x", fontsize=10)
    ax[1].loglog(Ns, ee, "o-", color=COLORS["teal"], label="ellipse")
    ax[1].loglog(Ns, np.maximum(ec, 1e-16), "s--", color=COLORS["muted"], label="circle (round-off)")
    ax[1].set_xlabel("N panels")
    ax[1].set_ylabel("max |ΔC_p|")
    ax[1].legend(fontsize=8)
    ax[1].set_title("panel convergence", fontsize=10)
    save(fig, out, "n90_source_panels")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
