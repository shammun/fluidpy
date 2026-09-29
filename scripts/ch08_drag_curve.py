"""§8.6 (C14–C15): sphere drag coefficient at low Re — Stokes 24/Re (8.52), Oseen (24/Re)(1 + 3Re/16), Proudman–Pearson,
and the Morrison correlation (ch04, "experiments") lying between Stokes and Oseen; Re = 2aU/ν throughout.

Run: ``.venv/Scripts/python.exe scripts/ch08_drag_curve.py --no-show``
Figure → outputs/ch08/c14_drag_curve.png.
"""
from __future__ import annotations

import numpy as np

from ch08_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch08_laminar_flow as ch08
from fluidpy.core import similarity as SIM


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    Re = np.logspace(-2, 1, 200)
    st, os_, pp, mo = (ch08.stokes_drag_coefficient(Re), ch08.oseen_drag_coefficient(Re),
                       ch08.proudman_pearson_drag_coefficient(Re), SIM.sphere_drag_coefficient(Re, "morrison"))
    print(f"parity with core.similarity 'stokes': {np.max(np.abs(st - SIM.sphere_drag_coefficient(Re, 'stokes'))):.1e}")
    for R in (0.1, 0.5, 1.0, 2.0, 5.0):
        print(f"Re = {R}: C_D Stokes {float(ch08.stokes_drag_coefficient(R)):.4f}, Oseen "
              f"{float(ch08.oseen_drag_coefficient(R)):.4f}, Proudman–Pearson {float(ch08.proudman_pearson_drag_coefficient(R)):.4f}, "
              f"Morrison {float(SIM.sphere_drag_coefficient(R)):.4f}")
    band = (Re >= 0.1) & (Re <= 5)
    inside = np.all((mo[band] >= st[band]) & (mo[band] <= os_[band]))
    print(f"Morrison correlation between Stokes and Oseen for 0.1 ≤ Re ≤ 5: {inside}")
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
    ax[0].loglog(Re, st, color=COLORS["accent"], lw=2, label="Stokes 24/Re (8.52)")
    ax[0].loglog(Re, os_, color=COLORS["orange"], lw=2, label="Oseen (24/Re)(1 + 3Re/16)")
    ax[0].loglog(Re, pp, color=COLORS["teal"], lw=1.5, ls="--", label="Proudman–Pearson")
    ax[0].loglog(Re, mo, color=COLORS["ink"], lw=1.2, ls=":", label="Morrison correlation (data fit)")
    ax[0].set_xlabel("Re = 2aU/ν")
    ax[0].set_ylabel("C_D")
    ax[0].legend(fontsize=8)
    ax[0].set_title("low-Re sphere drag")
    for y, lab, c in ((os_ / st, "Oseen", COLORS["orange"]), (pp / st, "Proudman–Pearson", COLORS["teal"]),
                      (mo / st, "Morrison", COLORS["ink"])):
        ax[1].semilogx(Re, y, color=c, lw=2, label=lab)
    ax[1].axhline(1, color=COLORS["accent"], lw=1)
    ax[1].set_xlabel("Re")
    ax[1].set_ylabel("C_D / (24/Re)")
    ax[1].legend(fontsize=8)
    ax[1].set_title("ratio to Stokes")
    save(fig, out, "c14_drag_curve")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
