"""§9.7 (C09, D13): OUR separated-flow model of the circular cylinder — ideal C_p = 1 - 4 sin^2(phi) up to the separation angle, a constant base
pressure behind — gives the pressure drag C_D = sin(phi_s)(1 - Cp_base) - (4/3) sin^3(phi_s).  QUALITATIVE: later separation lowers the drag, and
separation at 180 degrees is d'Alembert's paradox (C_D = 0).

Run: ``.venv/Scripts/python.exe scripts/ch09_cylinder_drag_model.py --no-show``   Figure -> outputs/ch09/c09_cylinder_drag_model.png.
"""
from __future__ import annotations

import numpy as np

from ch09_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch09_boundary_layers as ch09


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    print(f"d'Alembert limit: separation at 179.999 deg -> C_D = {ch09.separated_pressure_drag(179.999, None):.2e}")
    for ps, cb in ((82.0, -1.2), (125.0, -0.6)):
        print(f"  phi_sep = {ps:>5.0f} deg, base C_p = {cb:+.1f}: pressure drag {ch09.separated_pressure_drag(ps, cb):.3f} (closed form {ch09.BB._separated_drag_closed(ps, cb):.3f})")
    ps = np.linspace(40, 179, 120)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for cb, col in ((None, COLORS["accent"]), (-1.2, COLORS["orange"]), (-0.6, COLORS["teal"])):
        cd = [ch09.separated_pressure_drag(p, cb) for p in ps]
        ax[0].plot(ps, cd, color=col, label="base C_p = value at separation" if cb is None else f"base C_p = {cb}")
    ax[0].set_xlabel("separation angle from the forward stagnation point  [deg]")
    ax[0].set_ylabel("pressure-drag coefficient (model)")
    ax[0].set_ylim(-0.2, 3)
    ax[0].legend(fontsize=8)
    ax[0].set_title("later separation, less form drag (qualitative)")
    phi = np.linspace(0, 180, 361)
    ax[1].plot(phi, ch09.cp_ideal_cylinder(phi), color=COLORS["muted"], ls="--", label="ideal 1 - 4 sin^2")
    ax[1].plot(phi, ch09.separated_cp(phi, 82.0, -1.2), color=COLORS["orange"], label="separated at 82 deg")
    ax[1].plot(phi, ch09.separated_cp(phi, 125.0, -0.6), color=COLORS["teal"], label="separated at 125 deg")
    ax[1].set_xlabel("phi from the front stagnation point  [deg]")
    ax[1].set_ylabel("C_p")
    ax[1].legend(fontsize=8)
    ax[1].set_title("pressure round the cylinder (model)")
    save(fig, out, "c09_cylinder_drag_model")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
