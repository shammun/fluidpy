"""§9.4 (C05): the Falkner-Skan family f''' + (n+1)/2 f f'' - n f'^2 + n = 0 (Fig. 9.7 with our curves), wall shear vs n, the separation member
n = -0.0904 (f''(0) = 0) and the second (reversed-flow) branch.

Run: ``.venv/Scripts/python.exe scripts/ch09_falkner_skan.py --no-show``   Figures -> outputs/ch09/c05_falkner_skan.png.
"""
from __future__ import annotations

import numpy as np

from ch09_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch09_boundary_layers as ch09


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    sep = ch09.falkner_skan_separation()
    print(f"separation member: n_sep = {sep['m_sep']:.5f}, beta = 2n/(n+1) = {sep['beta_sep']:.5f}")
    ns = [4.0, 1.0, 1 / 3, 1 / 9, 0.0, -0.0654, -0.0904]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    print("       n   f''(0)       H    lambda  inflection eta")
    for n in ns:
        d = ch09.falkner_skan(n)
        s = ch09.falkner_skan_state(n)
        ax[0].plot(np.sqrt((n + 1) / 2) * d["eta"], d["fp"], label=f"n = {n:.4g}")
        infl = "none" if s["inflection_eta"] is None else f"{s['inflection_eta']:.3f}"
        print(f"{n:>8.4f}{s['fpp0']:>10.4f}{s['H']:>8.3f}{s['lam']:>10.4f}{infl:>16}")
    ax[0].set_xlim(0, 4)
    ax[0].set_xlabel("sqrt((n+1)/2) eta")
    ax[0].set_ylabel("u / U_e = f'(eta)")
    ax[0].legend(fontsize=8)
    ax[0].set_title("Fig. 9.7-like family (ours)")
    grid = np.concatenate([np.linspace(-0.09, -0.01, 12), np.linspace(0, 1, 15), np.linspace(1.2, 6, 10)])
    fpp0 = [ch09.falkner_skan_state(float(n))["fpp0"] for n in grid]
    ax[1].plot(grid, fpp0, "o-", color=COLORS["accent"], ms=3, label="attached branch")
    rv = ch09.falkner_skan(-0.08, branch="reversed")
    ax[1].plot([-0.08], [rv["fpp0"]], "s", color=COLORS["rose"], label="reversed branch (n = -0.08)")
    ax[1].axhline(0, color=COLORS["muted"], lw=0.8)
    ax[1].axvline(sep["m_sep"], color=COLORS["orange"], ls="--", label="fold n = -0.0904")
    ax[1].set_xlabel("n")
    ax[1].set_ylabel("f''(0)  (wall shear)")
    ax[1].set_xlim(-0.1, 1.5)
    ax[1].legend(fontsize=8)
    ax[1].set_title("wall shear vs pressure-gradient exponent")
    print(f"reversed branch at n = -0.08: f''(0) = {rv['fpp0']:.4f}, min f' = {rv['fp'].min():.4f} (reverse flow near the wall), success = {rv['success']}")
    save(fig, out, "c05_falkner_skan")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
