"""Shared command line and figure helpers for the ch10 scripts (not a physics module: every number comes from
``fluidpy.ch10_computational_fluid_dynamics``).  The notebook may import it with ``sys.path.insert(0, "scripts")``.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fluidpy.core.style import COLORS  # noqa: E402,F401

REF = ROOT / "reference" / "ch10"


def parse_args(doc: str, extra=None):
    """The common command line of every ch10 script: --out, --no-show, --fast (+ script-specific flags)."""
    ap = argparse.ArgumentParser(description=(doc or "").split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch10"))
    ap.add_argument("--no-show", action="store_true")
    ap.add_argument("--fast", action="store_true")
    if extra:
        extra(ap)
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


class Timer:
    """with Timer("label"): … prints the elapsed wall time."""

    def __init__(self, label: str):
        self.label = label

    def __enter__(self):
        self.t0 = time.time()
        return self

    def __exit__(self, *exc):
        self.dt = time.time() - self.t0
        print(f"  [{self.label}: {self.dt:.1f} s]", flush=True)


def write_csv(path: Path, header_lines, columns: dict, fmt: str = "%.6g") -> Path:
    """Write a small CSV of our own data with '#' header lines (provenance)."""
    import numpy as np

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    names = list(columns)
    data = np.column_stack([np.asarray(columns[k], dtype=float) for k in names])
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for h in header_lines:
            fh.write(f"# {h}\n")
        fh.write(",".join(names) + "\n")
        for row in data:
            fh.write(",".join(fmt % v for v in row) + "\n")
    print(f"  wrote {path}")
    return path
