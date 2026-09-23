"""Small matplotlib drawing helpers for the ch03 notebook and scripts, plus re-exports of the four figure builders
(Figs. 3.7, 3.14, 3.17, 3.18). Nothing here is physics; the numbers come from ``fluidpy.ch03_kinematics``.

    from scripts.ch03_drawings import ring_arrows, paddle, flow_lines_figure, shear_elements_frames, …

Run as a script it draws a demo (the relative-velocity ring of a simple shear split into strain and rotation, and a
paddle wheel): ``.venv/Scripts/python.exe scripts/ch03_drawings.py --no-show`` → outputs/ch03/drawings_demo.png.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]  # repo root, so ``import fluidpy`` works when run or imported standalone
for _p in (_ROOT, _ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from fluidpy import ch03_kinematics as ch03  # noqa: E402
from fluidpy.core.style import COLORS  # noqa: E402

_PART_COLORS = {"total": COLORS["accent"], "strain": COLORS["teal"], "rotation": COLORS["orange"]}


def ring_arrows(ax, G, radius: float = 1.0, parts=("total",), n: int = 16, center=(0.0, 0.0), scale: float = 1.0):
    """Relative-velocity arrows du on a ring of neighbours dx = radius·(cos θ, sin θ) around ``center``:
    ``parts`` ⊂ {"total" (G·dx, Eq. 3.10, purple), "strain" (S·dx, teal), "rotation" (½ω × dx, orange)} — the split of
    Eq. (3.19). Arrows are drawn at the neighbour positions, length × ``scale``. Returns the dict of du arrays."""
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    dx = radius * np.stack([np.cos(th), np.sin(th)])
    split = ch03.relative_velocity_split(np.asarray(G, float)[:2, :2], dx)
    arrays = {"total": split.du, "strain": split.du_strain, "rotation": split.du_rot}
    c = np.asarray(center, float)
    ax.plot(c[0] + radius * np.cos(np.linspace(0, 2 * np.pi, 100)), c[1] + radius * np.sin(np.linspace(0, 2 * np.pi, 100)),
            color=COLORS["grid"], lw=1)
    for p in parts:
        du = arrays[p]
        ax.quiver(c[0] + dx[0], c[1] + dx[1], du[0], du[1], color=_PART_COLORS[p], angles="xy", scale_units="xy",
                  scale=1.0 / scale, width=0.005, label=p)
    return {p: arrays[p] for p in parts}


def paddle(ax, xy, angle: float, size: float = 0.15, color=COLORS["ink"], n_blades: int = 4):
    """A paddle wheel at ``xy`` turned by ``angle`` [rad] (counter-clockwise positive): ``n_blades`` blades of length
    ``size`` and a hub. Used for the "does the element spin?" pictures (C10, C13). Returns the list of artists."""
    x0, y0 = xy
    arts = []
    for k in range(n_blades):
        a = angle + 2 * np.pi * k / n_blades
        (ln,) = ax.plot([x0, x0 + size * np.cos(a)], [y0, y0 + size * np.sin(a)], color=color, lw=2)
        arts.append(ln)
    (hub,) = ax.plot([x0], [y0], "o", color=color, ms=4)
    arts.append(hub)
    return arts


try:  # the four figure builders (C.6 of the design); both import styles work
    from scripts.ch03_fig3_7_flow_lines import flow_lines_animation, flow_lines_figure, streamlines_three_instants  # noqa: E402,F401
    from scripts.ch03_fig3_14_shear_elements import shear_elements_figure, shear_elements_frames  # noqa: E402,F401
    from scripts.ch03_fig3_17_leibniz import leibniz_strips_figure  # noqa: E402,F401
    from scripts.ch03_fig3_18_rtt import dt_convergence_figure, rtt_blob_figure  # noqa: E402,F401
except ModuleNotFoundError:  # pragma: no cover - when scripts/ itself is the import root
    from ch03_fig3_7_flow_lines import flow_lines_animation, flow_lines_figure, streamlines_three_instants  # noqa: E402,F401
    from ch03_fig3_14_shear_elements import shear_elements_figure, shear_elements_frames  # noqa: E402,F401
    from ch03_fig3_17_leibniz import leibniz_strips_figure  # noqa: E402,F401
    from ch03_fig3_18_rtt import dt_convergence_figure, rtt_blob_figure  # noqa: E402,F401


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(_ROOT / "outputs" / "ch03"))
    ap.add_argument("--no-show", action="store_true")
    args = ap.parse_args()
    import matplotlib

    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from fluidpy.core.style import use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    G = ch03.velocity_gradient_preset("simple_shear", 1.0)
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.8))
    for ax, parts in zip(axs, (("total",), ("strain",), ("rotation",))):
        d = ring_arrows(ax, G, 1.0, parts, scale=0.5)
        ax.set_aspect("equal")
        ax.set_xlim(-1.8, 1.8)
        ax.set_ylim(-1.8, 1.8)
        ax.set_title({"total": "du = G·dx (3.10)", "strain": "S·dx", "rotation": "½ω × dx"}[parts[0]])
        ax.set_xlabel("$dx_1$ [m]")
        print(f"{parts[0]:8s} du at dx = (0, 1): {d[parts[0]][:, 4].round(6)}")  # θ = 90° is index 4 of 16
    paddle(axs[2], (0, 0), -0.3, 0.4, COLORS["orange"])
    axs[0].set_ylabel("$dx_2$ [m]")
    fig.savefig(out / "drawings_demo.png", bbox_inches="tight")
    print(f"figure → {out / 'drawings_demo.png'}")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
