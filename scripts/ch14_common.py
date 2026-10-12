"""Shared command line and figure helpers for the ch14 scripts (not a physics module: every number comes from
``fluidpy.ch14_aerodynamics`` and the two core modules it re-exports).  The notebook may import it with
``sys.path.insert(0, "scripts")``.

All inputs are ``ch14.illustrative_inputs()`` — ours, never the book's.
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

# the chapter's colour code (curation §5): lift teal, drag rose, thrust/drive accent, circulation amber,
# pressure/suction orange, stream and velocities blue, references muted
C_LIFT, C_DRAG, C_DRIVE = COLORS["teal"], COLORS["rose"], COLORS["accent"]
C_CIRC, C_PRES, C_FLOW, C_REF = COLORS["amber"], COLORS["orange"], COLORS["blue"], COLORS["muted"]


def parse_args(doc: str, extra=None):
    """The common command line of every ch14 script: --out, --no-show, --fast (+ script-specific flags)."""
    ap = argparse.ArgumentParser(description=(doc or "").split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch14"))
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


def streamlines(ax, flow, xlim, ylim, n=(260, 200), levels=31, color=None, lw=0.8):
    """Contours of the stream function of a ``ZhukhovskyFlow`` (NaN inside the section) and the filled section."""
    import numpy as np

    from fluidpy.core import airfoils as AF

    x = np.linspace(*xlim, n[0])
    y = np.linspace(*ylim, n[1])
    X, Y = np.meshgrid(x, y, indexing="xy")
    psi = flow.stream(X, Y)
    psi0 = float(flow.stream(4.0 * flow.b, 0.0)) if False else None
    lv = np.linspace(np.nanpercentile(psi, 1), np.nanpercentile(psi, 99), levels)
    ax.contour(X, Y, psi, levels=lv, colors=color or C_FLOW, linewidths=lw, linestyles="solid")
    A = AF.zhukhovsky_airfoil(flow.b, flow.e, flow.beta, 600)
    ax.fill(A["x"], A["y"], color="#d5d9e4", zorder=3)
    ax.plot(np.append(A["x"], A["x"][0]), np.append(A["y"], A["y"][0]), color=COLORS["ink"], lw=1.2, zorder=4)
    ax.set(aspect="equal", xlim=xlim, ylim=ylim)
    return psi0


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
