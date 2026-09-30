"""§9.8 (N77, N80): the flow regimes of a circular cylinder against Re = U d / nu on a log axis, with the Strouhal shedding frequency for a wire.
QUALITATIVE: the thresholds are rounded experimental values (``BB.CYLINDER_THRESHOLDS``), not derived.

Run: ``.venv/Scripts/python.exe scripts/ch09_cylinder_regimes.py --no-show``   Figure -> outputs/ch09/c10_cylinder_regimes.png.
"""
from __future__ import annotations

import numpy as np

from ch09_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch09_boundary_layers as ch09


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    Re = np.logspace(-1, 7, 400)
    labels = [ch09.cylinder_flow_regime(r)["label"] for r in Re]
    edges = [0]
    for i in range(1, len(labels)):
        if labels[i] != labels[i - 1]:
            edges.append(i)
    edges.append(len(labels))
    fig, ax = plt.subplots(figsize=(11, 3.8))
    cols = [COLORS[k] for k in ("accent", "teal", "orange", "blue", "rose", "amber")] * 3
    for j, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
        lo, hi = Re[a], Re[min(b, len(Re) - 1)]
        ax.axvspan(lo, hi, color=cols[j], alpha=0.25)
        ax.text(np.sqrt(lo * hi), 0.5, labels[a].replace(" (", "\n("), rotation=90, ha="center", va="center", fontsize=7)
        print(f"  Re {lo:9.3g} - {hi:9.3g}: {labels[a]}")
    ax.set_xscale("log")
    ax.set_yticks([])
    ax.set_xlabel("Re = U d / nu")
    ax.set_title("cylinder flow regimes (rounded thresholds, qualitative)")
    save(fig, out, "c10_cylinder_regimes")
    d_wire, U = 2e-3, 20.0
    f = float(ch09.shedding_frequency(U, d_wire, 0.2)["f"])
    print(f"wire d = 2 mm in a 20 m/s wind: f = St U/d = {f:.0f} Hz (audible); if Omega is angular: {float(ch09.shedding_angular_frequency(U, d_wire)):.0f} rad/s")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
