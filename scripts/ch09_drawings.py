"""Drawing helpers for the ch09 scripts and notebook (design Part C.5): the boundary layer on a plate, a cylinder with its wake, a free / wall jet,
a teacup in cross-section and profile arrows — plus the common command line re-exported from ``ch09_common``.
Not a physics module: every number comes from ``fluidpy.ch09_boundary_layers``; the notebook imports this file with
``sys.path.insert(0, "scripts")``.  Each figure function takes ``ax=None`` and returns the ``Figure``; all drawings are schematic ("ours").
"""
from __future__ import annotations

import numpy as np

from ch09_common import COLORS, ROOT, finish, parse_args, save, setup  # noqa: F401

NU_AIR, RHO_AIR = 1.5e-5, 1.2  # m²/s, kg/m³ (air, ~20 °C)
NU_W, RHO_W = 1.0e-6, 1000.0  # m²/s, kg/m³ (water)
WALL = "#c8b99a"
EDGE = "#6b5b3e"


def _axes(ax, figsize=(6.0, 2.6)):
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    return fig, ax


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


def plate_layer(ax=None, U: float = 1.0, nu: float = NU_AIR, L: float = 1.0, exaggeration: float = 1.0):
    """Flat plate with the Blasius layer δ₉₉(x) = 4.910 √(νx/U) drawn over 0 < x < L; the thickness is multiplied by ``exaggeration`` (1 = to scale).

    Book: §9.1 (Fig. 9.1 idea, redrawn) with Eq. (9.30).  U [m/s]; nu [m²/s]; L [m].  Schematic: the outer streamlines are the parallel ideal ones."""
    fig, ax = _axes(ax)
    x = np.linspace(1e-4 * L, L, 300)
    d = exaggeration * 4.910 * np.sqrt(nu * x / U)
    dmax = float(d.max())
    ax.fill_between(x, 0, d, color=COLORS["accent"], alpha=0.18, lw=0)
    ax.plot(x, d, color=COLORS["accent"], lw=1.8, label=f"δ₉₉ (×{exaggeration:g})")
    ax.fill_between([0, L], [-0.06 * dmax] * 2, [0, 0], color=WALL, hatch="///", edgecolor=EDGE, lw=0)
    for f in (1.6, 2.3, 3.0):
        ax.plot(x, np.maximum(f * dmax, d + 0.6 * dmax), color=COLORS["muted"], lw=0.8)
    ax.annotate("", xy=(0.12 * L, 2.6 * dmax), xytext=(-0.08 * L, 2.6 * dmax), arrowprops=dict(arrowstyle="->", color=COLORS["orange"], lw=1.6))
    ax.text(-0.08 * L, 2.85 * dmax, f"U = {U:g} m/s", color=COLORS["orange"], fontsize=8)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title("boundary layer on a plate (schematic, ours)", fontsize=9)
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    return fig


def cylinder_sketch(ax=None, phi_sep_deg: float = 82.0, R: float = 1.0):
    """Cylinder in a stream from the left with the separation points at ±φ_s (measured from the forward stagnation point) and a wake region behind.

    Book: §9.7–9.8 (separation angles 82° laminar / 125° turbulent).  Schematic, ours; phi_sep_deg [deg]; R [any length unit]."""
    fig, ax = _axes(ax, (4.2, 3.0))
    th = np.linspace(0, 2 * np.pi, 300)
    ax.fill(R * np.cos(th), R * np.sin(th), color=COLORS["muted"], alpha=0.6, zorder=2)
    ps = np.deg2rad(phi_sep_deg)
    # angle from the forward stagnation point (at −x): point = (−R cos φ, R sin φ)
    for s in (+1, -1):
        px, py = -R * np.cos(ps), s * R * np.sin(ps)
        ax.plot(px, py, "o", color=COLORS["rose"], zorder=3)
        d = np.array([1.0, 0.0])
        ax.plot([px, px + 2.4 * R], [py, py * 0.45], color=COLORS["rose"], lw=1.2, ls="--")
    xw = np.linspace(-R * np.cos(ps), 3.0 * R, 60)
    ax.fill_between(xw, -np.abs(R * np.sin(ps)) * (1 - 0.45 * (xw - xw[0]) / (xw[-1] - xw[0])),
                    np.abs(R * np.sin(ps)) * (1 - 0.45 * (xw - xw[0]) / (xw[-1] - xw[0])), color=COLORS["rose"], alpha=0.08, zorder=1)
    for y0 in (-1.9, -1.4, 1.4, 1.9):
        ax.annotate("", xy=(-1.2 * R, y0 * R), xytext=(-2.6 * R, y0 * R), arrowprops=dict(arrowstyle="->", color=COLORS["blue"], lw=1.0))
    ax.set_aspect("equal")
    ax.set_xlim(-2.7 * R, 3.1 * R)
    ax.set_ylim(-2.2 * R, 2.2 * R)
    ax.set_title(f"separation at ±{phi_sep_deg:g}° (schematic, ours)", fontsize=9)
    ax.axis("off")
    return fig


def jet_sketch(ax=None, kind: str = "free", L: float = 1.0):
    """Sketch of a plane free jet (spreads symmetrically, δ ∝ x^{2/3}) or a wall jet (along a wall, δ ∝ x^{3/4}); L [any length unit].

    Book: §9.10 (Eqs. (9.64), (9.82) for the growth laws).  Schematic, ours."""
    fig, ax = _axes(ax, (5.0, 2.4))
    x = np.linspace(0.02, 1.0, 200) * L
    if kind == "free":
        h = 0.5 * x ** (2.0 / 3.0) * L ** (1.0 / 3.0)
        ax.fill_between(x, -h, h, color=COLORS["teal"], alpha=0.2, lw=0)
        ax.plot(x, h, color=COLORS["teal"], lw=1.5)
        ax.plot(x, -h, color=COLORS["teal"], lw=1.5)
        ax.plot(x, 0 * x, color=COLORS["muted"], lw=0.6, ls=":")
        ax.set_title("free jet: half-width ∝ x^{2/3} (schematic, ours)", fontsize=9)
        ax.set_ylim(-0.7 * L, 0.7 * L)
    elif kind == "wall":
        h = 0.35 * x ** 0.75 * L ** 0.25
        ax.fill_between(x, 0, h, color=COLORS["orange"], alpha=0.2, lw=0)
        ax.plot(x, h, color=COLORS["orange"], lw=1.5)
        ax.fill_between([0, L], [-0.05 * L] * 2, [0, 0], color=WALL, hatch="///", edgecolor=EDGE, lw=0)
        ax.set_title("wall jet: thickness ∝ x^{3/4} (schematic, ours)", fontsize=9)
        ax.set_ylim(-0.08 * L, 0.5 * L)
    else:
        raise ValueError('kind must be "free" or "wall"')
    ax.set_xlim(0, 1.05 * L)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    return fig


def cup_section(ax=None, R: float = 1.0, H: float = 0.8):
    """A teacup in cross-section: swirling water (circle arrows), the meridional secondary loop (in along the floor, up the axis, out along the top, down
    the side wall) and a tea leaf drifting inward on the floor.  Book: §9.11 (teacup).  R, H [any length unit].  Schematic, ours."""
    fig, ax = _axes(ax, (4.0, 3.2))
    ax.plot([-R, -R, R, R], [H, 0, 0, H], color=COLORS["ink"], lw=2)
    ax.fill_between([-R, R], 0, 0.9 * H, color=COLORS["blue"], alpha=0.10, lw=0)
    loop = dict(arrowstyle="->", color=COLORS["accent"], lw=1.4)
    ax.annotate("", xy=(-0.1 * R, 0.08 * H), xytext=(0.75 * R, 0.08 * H), arrowprops=loop)   # in along the floor
    ax.annotate("", xy=(0.0, 0.7 * H), xytext=(0.0, 0.15 * H), arrowprops=loop)              # up the axis
    ax.annotate("", xy=(0.75 * R, 0.8 * H), xytext=(0.1 * R, 0.8 * H), arrowprops=loop)      # out along the top
    ax.annotate("", xy=(0.9 * R, 0.2 * H), xytext=(0.9 * R, 0.75 * H), arrowprops=loop)      # down the side wall
    ax.plot([0.35 * R], [0.03 * H], "s", color=COLORS["amber"], ms=4)
    ax.text(-0.95 * R, 0.95 * H, "swirl about the axis (⊙)", color=COLORS["blue"], fontsize=8)
    ax.set_xlim(-1.1 * R, 1.1 * R)
    ax.set_ylim(-0.1 * H, 1.15 * H)
    ax.set_title("secondary flow in a stirred cup (schematic, ours)", fontsize=9)
    ax.axis("off")
    return fig


if __name__ == "__main__":
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(2, 3, figsize=(11, 6))
    plate_layer(axs[0, 0])
    cylinder_sketch(axs[0, 1])
    jet_sketch(axs[0, 2], "free")
    jet_sketch(axs[1, 0], "wall")
    cup_section(axs[1, 1])
    yy = np.linspace(0, 1, 10)
    profile_arrows(axs[1, 2], yy, np.sin(0.5 * np.pi * yy))
    fig.tight_layout()
    save(fig, out, "drawings")
    finish(args)
