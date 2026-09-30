"""§9.9 (N85): the sign of the Magnus effect as a truth table in (Re_low, Re_high) around the critical Reynolds number — negative Magnus only when
the fast side is past the drag crisis and the slow side is not.  QUALITATIVE (logic table).  The book's sentence prints 'Re < Re_cr' twice; the second
is Re > Re_cr (slip R10).

Run: ``.venv/Scripts/python.exe scripts/ch09_magnus_sign.py --no-show``   Figure -> outputs/ch09/n85_magnus_sign.png.
"""
from __future__ import annotations

import numpy as np

from ch09_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch09_boundary_layers as ch09


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    Re_cr = 3e5
    for lo, hi in ((1e5, 2e5), (1e5, 4e5), (4e5, 6e5), (2e5, 2e5)):
        print(f"  Re_low = {lo:.0e}, Re_high = {hi:.0e}  ->  Magnus sign {ch09.magnus_sign(lo, hi, Re_cr)}")
    g = np.logspace(4.5, 6.2, 60)
    lo, hi = np.meshgrid(g, g)
    val = np.array([[{"+": 1, "−": -1, "none": 0}[ch09.magnus_sign(a, b, Re_cr)] for a, b in zip(ra, rb)] for ra, rb in zip(lo, hi)])
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.pcolormesh(lo, hi, val, cmap="coolwarm_r", shading="auto", vmin=-1, vmax=1)
    ax.axvline(Re_cr, color="k", lw=0.8)
    ax.axhline(Re_cr, color="k", lw=0.8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Re on the slow side")
    ax.set_ylabel("Re on the fast side")
    ax.text(1.2e5, 6e5, "negative", fontsize=10, color=COLORS["ink"])
    ax.text(1e5, 6e4, "none (no spin)", fontsize=8)
    ax.text(4e5, 4.5e5, "positive", fontsize=10)
    ax.set_title("sign of the Magnus force (truth table, qualitative)")
    save(fig, out, "n85_magnus_sign")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
