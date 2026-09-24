"""Shared drawing helpers for the ch07 scripts and notebook (tank outlines, wave surfaces, orbit ghosts, the hydraulic-
jump control volume, beam labels, the two-layer sketch) plus the common command line. Not a physics module: every
number comes from ``fluidpy.ch07_gravity_waves``; the notebook imports this file with ``sys.path.insert(0, "scripts")``.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fluidpy.core.style import COLORS  # noqa: E402

NAVY = "#1e3a8a"  # lower / denser layer (design colour convention)
LIGHT = "#dbeafe"  # upper / lighter layer


def parse_args(doc: str):
    """The common command line of every ch07 script: --out, --no-show, --fast."""
    ap = argparse.ArgumentParser(description=doc.split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch07"))
    ap.add_argument("--no-show", action="store_true")
    ap.add_argument("--fast", action="store_true")
    return ap.parse_args()


def setup(args):
    """Select a non-interactive backend when not showing, apply the house style, create the output folder."""
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


def tank(ax, H: float, L: float, x0: float = 0.0, label: bool = True):
    """Flat bottom at z = −H (hatched), dashed still level z = 0 and the z-axis arrow, for x in [x0, x0 + L]."""
    ax.fill_between([x0, x0 + L], [-H * 1.08] * 2, [-H] * 2, color="#c8b99a", zorder=1)
    ax.plot([x0, x0 + L], [-H, -H], color="#6b5b3e", lw=1.5, zorder=2)
    ax.plot([x0, x0 + L], [0, 0], color=COLORS["muted"], lw=0.8, ls="--", zorder=2)
    if label:
        ax.annotate("", xy=(x0, 0.25 * H), xytext=(x0, -H), arrowprops=dict(arrowstyle="->", color=COLORS["ink"]))
        ax.text(x0 + 0.01 * L, 0.2 * H, "z", color=COLORS["ink"])
        ax.text(x0 + L * 0.98, -H * 0.97, "z = −H", ha="right", va="bottom", fontsize=8, color=COLORS["muted"])
    ax.set_xlim(x0, x0 + L)


def wave_surface(ax, x, eta, fill_to: float | None = None, **kw):
    """The free surface (blue) and, if ``fill_to`` is given, water shading down to that level."""
    kw.setdefault("color", COLORS["blue"])
    kw.setdefault("lw", 2.0)
    ax.plot(x, eta, **kw)
    if fill_to is not None:
        ax.fill_between(x, fill_to, eta, color=COLORS["blue"], alpha=0.08, lw=0)


def orbit_ghosts(ax, x0s, z0s, A, B, n: int = 100, **kw):
    """Ellipses of semi-axes A (horizontal) and B (vertical) round each mean position (teal, orbit ghosts)."""
    kw.setdefault("color", COLORS["teal"])
    kw.setdefault("lw", 1.0)
    s = np.linspace(0, 2 * np.pi, n)
    for x0, z0, a, b in np.broadcast(x0s, z0s, A, B):
        ax.plot(x0 + a * np.cos(s), z0 + b * np.sin(s), **kw)


def cv_box(ax, H1: float, H2: float, x1: float = 0.0, x2: float = 1.0, xj: float = 0.5, width: float = 0.08):
    """A stationary hydraulic jump from depth H1 to H2 around x = xj with the dashed control volume of Fig. 7.20b and
    faces 1 and 2 labelled."""
    xs = np.linspace(x1, x2, 400)
    s = 0.5 * (1 + np.tanh((xs - xj) / width))
    surf = H1 + (H2 - H1) * s
    ax.fill_between(xs, 0, surf, color=COLORS["blue"], alpha=0.15, lw=0)
    ax.plot(xs, surf, color=COLORS["blue"], lw=2)
    ax.plot([x1, x2], [0, 0], color="#6b5b3e", lw=2)
    top = 1.25 * H2
    ax.plot([x1 + 0.05, x1 + 0.05, x2 - 0.05, x2 - 0.05, x1 + 0.05], [0, top, top, 0, 0], color=COLORS["ink"],
            ls="--", lw=1)
    ax.text(x1 + 0.05, top * 1.03, "1", ha="center", fontsize=10)
    ax.text(x2 - 0.05, top * 1.03, "2", ha="center", fontsize=10)


def beam_labels(ax, theta: float, r: float = 1.0):
    """Mark the four St Andrew's-cross beam directions at angle θ from the vertical with dashed rays and a label."""
    for sx in (1, -1):
        for sz in (1, -1):
            ax.plot([0, sx * r * np.sin(theta)], [0, sz * r * np.cos(theta)], color=COLORS["accent"], lw=1, ls="--")
    ax.text(0.05 * r, 0.75 * r, f"θ = {np.degrees(theta):.1f}° from the vertical", color=COLORS["accent"],
            fontsize=8)


def two_layer_sketch(ax, H: float, L: float, eta=None, zeta=None, x=None):
    """Upper layer (light) of thickness H over a deep lower layer (navy), with optional surface η(x) and interface
    −H + ζ(x) curves."""
    if x is None:
        x = np.linspace(0, L, 200)
    top = np.zeros_like(x) if eta is None else eta
    mid = -H + (np.zeros_like(x) if zeta is None else zeta)
    ax.fill_between(x, mid, top, color=LIGHT, lw=0)
    ax.fill_between(x, -3 * H, mid, color=NAVY, alpha=0.35, lw=0)
    ax.plot(x, top, color=COLORS["blue"], lw=2)
    ax.plot(x, mid, color=NAVY, lw=2)
    ax.set_xlim(x[0], x[-1])
