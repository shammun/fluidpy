"""§9.6 (N63, Fig. 9.8 idea): Thwaites' universal function L(lambda) = 2 l - 2 (2 + H) lambda computed on the exact Falkner-Skan family with our
solver, with the straight line 0.45 - 6.0 lambda of Eq. (9.49) overlaid.  Also the shear correlation l(lambda) and H(lambda) (our closure table).

Run: ``.venv/Scripts/python.exe scripts/ch09_thwaites_L.py --no-show``   Figure -> outputs/ch09/n63_thwaites_L.png.
"""
from __future__ import annotations

import numpy as np

from ch09_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch09_boundary_layers as ch09


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    t = ch09.thwaites_closure_table(30 if args.fast else 60)
    lam = t["lam"]
    L = 2 * t["l"] - 2 * (2 + t["H"]) * lam
    line = 0.45 - 6.0 * lam
    print(f"closure table spans lambda in [{lam[0]:.4f}, {lam[-1]:.4f}] (n = {lam.size} Falkner-Skan members)")
    for lv in (-0.0675, 0.0, 0.05, 0.085):
        Lv, Hv = float(ch09.thwaites_L(lv)), float(ch09.thwaites_H(lv))
        print(f"  lambda = {lv:+.4f}: L = {Lv:.4f} (line {0.45 - 6 * lv:.4f}), l = {float(ch09.thwaites_l(lv)):.4f}, H = {Hv:.3f}")
    fig, ax = plt.subplots(1, 3, figsize=(13, 4))
    ax[0].plot(lam, L, color=COLORS["accent"], label="L(lambda), Falkner-Skan family")
    ax[0].plot(lam, line, color=COLORS["orange"], ls="--", label="0.45 - 6.0 lambda  (9.49)")
    ax[0].set_xlabel("lambda")
    ax[0].set_ylabel("L(lambda)")
    ax[0].legend(fontsize=8)
    ax[0].set_title("Fig. 9.8 idea: L is nearly linear")
    ll = np.linspace(-0.09, 0.15, 200)
    ax[1].plot(lam, t["l"], "o", ms=3, color=COLORS["teal"], label="exact FS members")
    ax[1].plot(ll, ch09.thwaites_l(ll, "falkner_skan"), color=COLORS["teal"], lw=1)
    ax[1].plot(ll, ch09.thwaites_l(ll, "white"), color=COLORS["rose"], ls="--", label="(lambda + 0.09)^0.62")
    ax[1].set_xlabel("lambda")
    ax[1].set_ylabel("l(lambda)")
    ax[1].legend(fontsize=8)
    ax[1].set_title("shear correlation, tau0 = mu (U_e/theta) l")
    ax[2].plot(lam, t["H"], color=COLORS["blue"])
    ax[2].set_xlabel("lambda")
    ax[2].set_ylabel("H = delta*/theta")
    ax[2].set_title("shape factor")
    save(fig, out, "n63_thwaites_L")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
