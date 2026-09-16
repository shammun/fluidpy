"""Chapter 2, §2.1, Eq. (2.1) and Fig. 2.1 (our drawing): the position vector x = e_i x_i, its three components and
the unit vectors, plus the basis reconstruction of (2.1).

Run: ``.venv/Scripts/python.exe scripts/ch02_fig2_1_position_vector.py --no-show``
Figure → outputs/ch02/fig2_1_position_vector.png.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch02"))
    ap.add_argument("--no-show", action="store_true")
    args = ap.parse_args()
    t0 = time.perf_counter()
    import matplotlib

    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from ch02_drawings import arrow3d, setup_3d
    from fluidpy import ch02_cartesian_tensors as ch02
    from fluidpy.core.style import COLORS, use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    x = np.array([1.0, 1.3, 0.9])  # the running example (any units: metres)
    E = ch02.unit_vectors(3)
    assert np.allclose(ch02.vector_from_components(x, E), x)  # Eq. (2.1): x = e_1 x_1 + e_2 x_2 + e_3 x_3

    fig = plt.figure(figsize=(6.5, 5.5))
    ax = fig.add_subplot(111, projection="3d")
    setup_3d(ax, lim=1.6)
    for e, name in zip(E, ("$\\mathbf{e}_1$", "$\\mathbf{e}_2$", "$\\mathbf{e}_3$")):
        arrow3d(ax, (0, 0, 0), 0.5 * e, color=COLORS["teal"], label=name, lw=1.6)
    arrow3d(ax, (0, 0, 0), x, color=COLORS["accent"], label="$\\mathbf{x}$ = OP", lw=2.5)
    # dashed projections: foot in the (1,2) plane and the three components
    foot = np.array([x[0], x[1], 0.0])
    for a, b in (((0, 0, 0), foot), (foot, x), ((x[0], 0, 0), foot), ((0, x[1], 0), foot)):
        ax.plot(*zip(a, b), color=COLORS["muted"], ls="--", lw=1.0)
    ax.text(x[0] / 2, -0.15, 0, "$x_1$", color=COLORS["ink"])
    ax.text(x[0] + 0.05, x[1] / 2, 0, "$x_2$", color=COLORS["ink"])
    ax.text(x[0] + 0.05, x[1] + 0.05, x[2] / 2, "$x_3$", color=COLORS["ink"])
    ax.set_title("Eq. (2.1): $\\mathbf{x} = \\mathbf{e}_1 x_1 + \\mathbf{e}_2 x_2 + \\mathbf{e}_3 x_3$ — one arrow, three shadows")
    ax.view_init(elev=22, azim=-60)
    fig.savefig(out / "fig2_1_position_vector.png", bbox_inches="tight")

    print(f"x = {x}, |x| = {np.linalg.norm(x):.4f}; e_i rows of unit_vectors():\n{E}")
    print(f"vector_from_components(x, E) = {ch02.vector_from_components(x, E)} (Eq. 2.1 reconstruction)")
    print(f"a·b for (1,2,3)·(4,5,6) = {ch02.dot([1, 2, 3], [4, 5, 6])} (Eq. 2.2); expand: {ch02.expand_indices('a_i b_i')}")
    print(f"saved fig2_1_position_vector.png in {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
