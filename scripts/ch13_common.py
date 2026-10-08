"""Shared command line and figure helpers for the ch13 scripts (not a physics module: every number comes from
``fluidpy.ch13_geophysical_fluid_dynamics`` and the three core modules it re-exports).  The notebook may import it with
``sys.path.insert(0, "scripts")``.

All inputs are ``ch13.illustrative_inputs()`` — ours, never the book's.
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

REF = ROOT / "reference" / "ch13"
DAY = 86400.0
YEAR = 365.25 * DAY


def parse_args(doc: str, extra=None):
    """The common command line of every ch13 script: --out, --no-show, --fast (+ script-specific flags)."""
    ap = argparse.ArgumentParser(description=(doc or "").split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch13"))
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


def hemispheres(lat_rad: float):
    """The latitude and its southern twin, with labels — every direction-dependent demo runs both."""
    import numpy as np

    return ((abs(lat_rad), f"{np.rad2deg(abs(lat_rad)):.0f}° N"), (-abs(lat_rad), f"{np.rad2deg(abs(lat_rad)):.0f}° S"))


class Timer:
    """``with Timer("label"):`` prints the elapsed wall time."""

    def __init__(self, label: str):
        self.label = label

    def __enter__(self):
        self.t0 = time.time()
        return self

    def __exit__(self, *exc):
        print(f"  [{self.label}: {time.time() - self.t0:.1f} s]")
        return False
