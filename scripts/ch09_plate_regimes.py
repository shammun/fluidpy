"""§9.7 (N68, Fig. 9.10 idea): regimes along a flat plate — leading edge (Re_x ~ 1), laminar similar region, transition, turbulent — drawn on a
Reynolds-number axis, with the Blasius thickness growth sqrt(x) in the laminar part.  QUALITATIVE: the thresholds are the book's rounded values
(arguments of ``transition_state``); the instability mechanism is Ch. 11.

Run: ``.venv/Scripts/python.exe scripts/ch09_plate_regimes.py --no-show``   Figure -> outputs/ch09/n68_plate_regimes.png.
"""
from __future__ import annotations

import numpy as np

from ch09_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch09_boundary_layers as ch09


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    U, nu = 20.0, 1.5e-5
    x = np.geomspace(1e-4, 10, 300)
    Rex = U * x / nu
    labels = [ch09.transition_state(r) for r in Rex]
    for lab in ("laminar", "transitional", "turbulent"):
        idx = [i for i, l in enumerate(labels) if l == lab]
        print(f"  {lab:<13} x in [{x[idx[0]]:.3g}, {x[idx[-1]]:.3g}] m  (Re_x {Rex[idx[0]]:.2g} - {Rex[idx[-1]]:.2g})")
    fig, ax = plt.subplots(figsize=(9, 3.6))
    d = ch09.blasius_delta99(x, U, nu)
    col = {"laminar": COLORS["teal"], "transitional": COLORS["orange"], "turbulent": COLORS["rose"]}
    for lab in col:
        m = np.array([l == lab for l in labels])
        ax.loglog(x[m], d[m] * (1 if lab == "laminar" else np.nan), color=col[lab], lw=3)
        ax.axvspan(x[m][0], x[m][-1], color=col[lab], alpha=0.12, label=lab)
    ax.loglog(x, d, color=COLORS["ink"], lw=1, ls=":", label="Blasius delta_99 (valid only where laminar)")
    ax.axvline(x[Rex >= 1.0][0], color=COLORS["muted"], lw=0.8)
    ax.text(x[Rex >= 1.0][0] * 1.1, d.max() * 0.3, "Re_x ~ 1: leading edge,\nboundary-layer theory fails", fontsize=8)
    ax.set_xlabel("x  [m]")
    ax.set_ylabel("delta_99  [m]")
    ax.legend(fontsize=8, loc="lower right")
    ax.set_title("regimes of a flat-plate layer (Re_cr thresholds are rounded, qualitative)")
    save(fig, out, "n68_plate_regimes")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
