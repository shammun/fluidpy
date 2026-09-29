"""Shared command line and figure helpers for the ch08 scripts (not a physics module: every number comes from
``fluidpy.ch08_laminar_flow``). The notebook may import it with ``sys.path.insert(0, "scripts")``.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fluidpy.core.style import COLORS  # noqa: E402,F401


def parse_args(doc: str):
    """The common command line of every ch08 script: --out, --no-show, --fast."""
    ap = argparse.ArgumentParser(description=(doc or "").split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch08"))
    ap.add_argument("--no-show", action="store_true")
    ap.add_argument("--fast", action="store_true")
    return ap.parse_args()


def setup(args) -> Path:
    """Non-interactive backend when not showing, house style, output folder."""
    import matplotlib

    if args.no_show:
        matplotlib.use("Agg")
    from fluidpy.core.style import use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    return out


def save(fig, out: Path, name: str) -> Path:
    path = Path(out) / f"{name}.png"
    fig.savefig(path, bbox_inches="tight")
    print(f"  saved {path}")
    return path


def finish(args) -> None:
    import matplotlib.pyplot as plt

    if not args.no_show:
        plt.show()
    plt.close("all")
