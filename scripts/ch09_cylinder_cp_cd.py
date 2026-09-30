"""§9.8 (N81, C11): surface pressure round a cylinder before and after the drag crisis (our separated model vs the ideal curve) and a SCHEMATIC C_D(Re)
with the crisis dip.  QUALITATIVE: no dataset was fetched for the book's Fig. 9.21; below Re = 1 the curve is Lamb's asymptote, elsewhere rounded anchor points.

Run: ``.venv/Scripts/python.exe scripts/ch09_cylinder_cp_cd.py --no-show``   Figure -> outputs/ch09/c11_cylinder_cp_cd.png.
"""
from __future__ import annotations

import numpy as np

from ch09_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch09_boundary_layers as ch09


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    dc = ch09.drag_crisis_pair()
    print(f"model pressure drag: subcritical {dc['subcritical']['cd_model']:.3f} -> supercritical {dc['supercritical']['cd_model']:.3f} (ratio {dc['ratio']:.2f})")
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    phi = np.linspace(0, 180, 361)
    ax[0].plot(phi, ch09.cp_ideal_cylinder(phi), color=COLORS["muted"], ls="--", label="ideal flow")
    for st, col in (("subcritical", COLORS["orange"]), ("supercritical", COLORS["teal"])):
        p = dc[st]
        ax[0].plot(phi, ch09.separated_cp(phi, p["phi_sep_deg"], p["cp_base"]), color=col, label=f"{st}: sep {p['phi_sep_deg']:.0f} deg, base {p['cp_base']}")
    ax[0].set_xlabel("phi from the front stagnation point [deg]")
    ax[0].set_ylabel("C_p")
    ax[0].legend(fontsize=8)
    ax[0].set_title("wake pressure plateau (model)")
    Re = np.logspace(-1, 7, 300)
    ax[1].loglog(Re, ch09.cylinder_cd_schematic(Re), color=COLORS["accent"], label="schematic (qualitative)")
    ax[1].axvline(3e5, color=COLORS["rose"], ls=":", label="Re_cr ~ 3e5")
    ax[1].set_xlabel("Re = U d / nu")
    ax[1].set_ylabel("C_D")
    ax[1].legend(fontsize=8)
    ax[1].set_title("drag crisis (schematic)")
    for r in (0.1, 1e3, 3e5, 1e6):
        print(f"  schematic C_D({r:.0e}) = {float(ch09.cylinder_cd_schematic(r)):.3g}")
    save(fig, out, "c11_cylinder_cp_cd")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
