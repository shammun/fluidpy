"""§9.7 (N69, Fig. 9.11 idea): flat-plate drag coefficient vs Re_L — the laminar line 1.328/sqrt(Re_L) (Eq. 9.33), a turbulent line
0.074 Re_L^(-1/5) (NOT in the book: Prandtl's one-seventh-power-law fit, valid 5e5 < Re_L < 1e7) and the transition-patched mixed curve.

Run: ``.venv/Scripts/python.exe scripts/ch09_plate_drag_curve.py --no-show``   Figure -> outputs/ch09/n69_plate_drag_curve.png.
"""
from __future__ import annotations

import numpy as np

from ch09_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch09_boundary_layers as ch09


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    R = np.logspace(5, 9, 300)
    lam = ch09.plate_drag_coefficient(R, "laminar")
    tur = ch09.plate_drag_coefficient(R, "turbulent")
    mix = ch09.plate_drag_coefficient(R, "mixed", Re_tr=5e5)
    for r in (1e5, 5e5, 1e6, 1e7, 1e9):
        print(f"  Re_L = {r:.0e}: laminar {float(ch09.plate_drag_coefficient(r, 'laminar')):.5f}, turbulent {float(ch09.plate_drag_coefficient(r, 'turbulent')):.5f}, "
              f"mixed {float(ch09.plate_drag_coefficient(r, 'mixed')):.5f}")
    print("  laminar layer wins below Re_tr = 5e5; mixed = 0.074 Re^-1/5 - Re_tr(C_t - C_l)/Re above it (one wetted face)")
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.loglog(R, lam, color=COLORS["teal"], label="laminar 1.328/sqrt(Re_L)  (9.33)")
    ax.loglog(R, tur, color=COLORS["rose"], label="turbulent 0.074 Re_L^(-1/5)  (not in book)")
    ax.loglog(R, mix, color=COLORS["accent"], lw=2.4, label="mixed: transition at Re_tr = 5e5")
    ax.set_xlabel("Re_L = U L / nu")
    ax.set_ylabel("C_D per wetted face")
    ax.legend(fontsize=8)
    ax.set_title("laminar, turbulent and transition-patched plate drag")
    save(fig, out, "n69_plate_drag_curve")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
