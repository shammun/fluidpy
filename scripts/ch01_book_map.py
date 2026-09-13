"""Chapter 1, §1.1: a map of the book — which chapters build on which, with the route to Ch. 13 (GFD) highlighted.

Importable helper for the notebook (C01)::

    from scripts.ch01_book_map import draw_book_map
    ax = draw_book_map()

Run: ``.venv/Scripts/python.exe scripts/ch01_book_map.py --no-show`` (figure → outputs/ch01/fig_book_map.png).
The dependency edges are typed here from the book's structure as summarised in the project skill ``fluids-book``
(our summary, not the book's text).
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

#: Prerequisite edges (from, to): chapter "to" leans on chapter "from".
EDGES: list[tuple[str, str]] = [
    ("ch01", "ch02"), ("ch02", "ch03"), ("ch03", "ch04"),
    ("ch04", "ch05"), ("ch04", "ch06"), ("ch04", "ch07"), ("ch04", "ch08"),
    ("ch08", "ch09"), ("ch04", "ch10"), ("ch08", "ch11"), ("ch11", "ch12"), ("ch09", "ch12"),
    ("ch04", "ch13"), ("ch05", "ch13"), ("ch07", "ch13"), ("ch08", "ch13"), ("ch11", "ch13"), ("ch12", "ch13"),
    ("ch06", "ch14"), ("ch09", "ch14"), ("ch04", "ch15"), ("ch08", "ch16"), ("ch09", "ch16"),
]

#: Column (depth) of each chapter in the drawing: foundations left, applications right.
LEVEL: dict[str, int] = {"ch01": 0, "ch02": 1, "ch03": 2, "ch04": 3, "ch05": 4, "ch06": 4, "ch07": 4, "ch08": 4,
                         "ch10": 4, "ch15": 4, "ch09": 5, "ch11": 5, "ch12": 6, "ch14": 6, "ch16": 6, "ch13": 7}


def _ancestors(target: str) -> set[str]:
    found, stack = set(), [target]
    while stack:
        node = stack.pop()
        for a, b in EDGES:
            if b == node and a not in found:
                found.add(a)
                stack.append(a)
    return found


def draw_book_map(ax=None, highlight=("ch13",)):
    """Draw the chapter dependency graph of the book on a matplotlib Axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes, optional
        Axes to draw on (a new 11 × 5.5 in figure if None).
    highlight : sequence of str, optional
        Chapter ids whose prerequisite routes are highlighted (default Ch. 13, geophysical fluid dynamics).

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt

    from fluidpy.core.project import chapters
    from fluidpy.core.style import COLORS

    chaps = {c["id"]: c for c in chapters()}
    if ax is None:
        _, ax = plt.subplots(figsize=(11, 5.5))
    cols: dict[int, list[str]] = {}
    for cid in sorted(chaps):
        cols.setdefault(LEVEL.get(cid, 4), []).append(cid)
    pos = {}
    for lev, ids in cols.items():
        for i, cid in enumerate(ids):
            pos[cid] = (lev, (len(ids) - 1) / 2.0 - i)
    route = set(highlight)
    for h in highlight:
        route |= _ancestors(h)
    for a, b in EDGES:
        if a in pos and b in pos:
            on = a in route and b in route
            ax.annotate("", xy=pos[b], xytext=pos[a],
                        arrowprops=dict(arrowstyle="-|>", lw=2.2 if on else 0.9, shrinkA=18, shrinkB=18,
                                        color=COLORS["orange"] if on else COLORS["grid"]))
    for cid, (x, y) in pos.items():
        on = cid in route
        ax.scatter([x], [y], s=900, color=COLORS["accent"] if cid in highlight else (COLORS["teal"] if on else "white"),
                   edgecolor=COLORS["ink"], zorder=3)
        ax.text(x, y, cid[2:], ha="center", va="center", fontsize=9, zorder=4,
                color="white" if on else COLORS["ink"], weight="bold")
        title = chaps[cid].get("title", "")
        ax.text(x, y - 0.36, title if len(title) < 22 else title[:20] + "…", ha="center", va="top", fontsize=7,
                color=COLORS["muted"])
    ax.set_xlim(-0.6, max(LEVEL.values()) + 0.6)
    ax.set_ylim(-3.2, 3.0)
    ax.axis("off")
    ax.set_title("How the chapters build on each other (orange: the route to " + ", ".join(highlight) + ")")
    return ax


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch01"), help="figure folder")
    ap.add_argument("--no-show", action="store_true", help="do not open a window")
    args = ap.parse_args()
    t0 = time.perf_counter()
    import matplotlib

    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from fluidpy.core.style import use_style

    use_style()
    ax = draw_book_map()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    ax.figure.savefig(out / "fig_book_map.png", bbox_inches="tight")
    print(f"chapters: {len(LEVEL)}  edges: {len(EDGES)}  prerequisites of ch13: {sorted(_ancestors('ch13'))}")
    print(f"saved {out / 'fig_book_map.png'}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
