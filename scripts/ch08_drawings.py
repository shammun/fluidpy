"""Drawing helpers for the ch08 scripts and notebook (design Part C.5): channel walls, a pipe section, an annulus, a
tilted bearing pad, a sphere, and profile arrows — plus the common command line re-exported from ``ch08_common``.
Not a physics module: every number comes from ``fluidpy.ch08_laminar_flow``; the notebook imports this file with
``sys.path.insert(0, "scripts")``.
"""
from __future__ import annotations

import numpy as np

from ch08_common import COLORS, ROOT, finish, parse_args, save, setup  # noqa: F401

WALL = "#c8b99a"
EDGE = "#6b5b3e"


def channel(ax, h: float, U: float = 0.0, L: float = 1.0, x0: float = 0.0, label: bool = True):
    """Fixed hatched wall at y = 0 and a wall at y = h (moving at U: arrow), drawn over x0 … x0 + L (data units)."""
    t = 0.06 * h
    ax.fill_between([x0, x0 + L], [-t] * 2, [0, 0], color=WALL, hatch="///", edgecolor=EDGE, lw=0, zorder=1)
    ax.fill_between([x0, x0 + L], [h, h], [h + t] * 2, color=COLORS["muted"] if U else WALL, zorder=1)
    ax.plot([x0, x0 + L], [0, 0], color=EDGE, lw=1.5)
    ax.plot([x0, x0 + L], [h, h], color=COLORS["ink"], lw=1.5)
    if U and label:
        ax.annotate("", xy=(x0 + 0.95 * L, h + 2.5 * t), xytext=(x0 + 0.7 * L, h + 2.5 * t),
                    arrowprops=dict(arrowstyle="->", color=COLORS["orange"], lw=1.8))
        ax.text(x0 + 0.72 * L, h + 3.2 * t, f"U = {U:g}", color=COLORS["orange"], fontsize=8)


def pipe_section(ax, a: float, n: int = 200):
    """Cross-section circle of radius a with the axis marked."""
    th = np.linspace(0, 2 * np.pi, n)
    ax.plot(a * np.cos(th), a * np.sin(th), color=COLORS["ink"], lw=2)
    ax.plot(0, 0, "+", color=COLORS["muted"])
    ax.set_aspect("equal")


def annulus(ax, R1: float, R2: float, n: int = 200):
    """Inner cylinder (filled) and outer wall of a circular-Couette gap; R2 = inf draws only the inner cylinder."""
    th = np.linspace(0, 2 * np.pi, n)
    ax.fill(R1 * np.cos(th), R1 * np.sin(th), color=COLORS["muted"], zorder=2)
    if np.isfinite(R2):
        ax.plot(R2 * np.cos(th), R2 * np.sin(th), color=COLORS["ink"], lw=2)
    ax.set_aspect("equal")


def pad(ax, h0: float, alpha: float, L: float, thickness: float | None = None, U: float | None = None):
    """Tilted bearing pad over a flat floor: gap h = h0(1 + αx/L) for 0 ≤ x ≤ L (Fig. 8.9 idea)."""
    x = np.linspace(0, L, 50)
    hx = h0 * (1 + alpha * x / L)
    tk = thickness if thickness is not None else 0.6 * h0 * (1 + alpha)
    ax.fill_between(x, hx, hx + tk, color=COLORS["muted"], zorder=2)
    ax.fill_between([-0.1 * L, 1.1 * L], [-0.15 * h0] * 2, [0, 0], color=WALL, hatch="///", edgecolor=EDGE, lw=0)
    ax.plot([-0.1 * L, 1.1 * L], [0, 0], color=EDGE, lw=1.5)
    if U:
        ax.annotate("", xy=(0.8 * L, hx[-1] + 1.3 * tk), xytext=(0.5 * L, hx[-1] + 1.3 * tk),
                    arrowprops=dict(arrowstyle="->", color=COLORS["orange"], lw=1.8))


def sphere(ax, a: float = 1.0, center=(0.0, 0.0), color: str | None = None):
    """Filled circle for a sphere of radius a in a meridian plane."""
    import matplotlib.pyplot as plt

    ax.add_patch(plt.Circle(center, a, color=color or COLORS["muted"], zorder=3))
    ax.set_aspect("equal")


def profile_arrows(ax, y, u, x0: float = 0.0, scale: float = 1.0, every: int = 1, color: str | None = None):
    """Velocity profile drawn as horizontal arrows from the line x = x0 (length scale·u) plus the envelope curve."""
    y, u = np.asarray(y, float), np.asarray(u, float)
    c = color or COLORS["accent"]
    for yi, ui in zip(y[::every], u[::every]):
        if abs(ui) > 0:
            ax.annotate("", xy=(x0 + scale * ui, yi), xytext=(x0, yi),
                        arrowprops=dict(arrowstyle="->", color=c, lw=1.0, shrinkA=0, shrinkB=0))
    ax.plot(x0 + scale * u, y, color=c, lw=1.8)
    ax.plot([x0, x0], [y.min(), y.max()], color=COLORS["grid"], lw=0.8)
