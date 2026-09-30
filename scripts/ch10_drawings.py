"""Schematic drawings for Chapter 10 (our own versions of Figs. 10.1–10.6, 10.9, 10.15, 10.17 and the stencil / characteristic
diagrams) — drawing helpers only, no physics.  Each helper draws on a given matplotlib axis; ``main`` draws them all.

Run: ``.venv/Scripts/python.exe scripts/ch10_drawings.py --no-show``   Figure -> outputs/ch10/drawings.png.
"""
from __future__ import annotations

import numpy as np

from ch10_common import COLORS, finish, parse_args, save, setup


def _clean(ax):
    ax.set_aspect("equal")
    ax.axis("off")


def spacetime_grid(ax, nx: int = 7):
    """Uniform space–time grid x_i = iΔx, t_n = nΔt with three time levels (our Fig. 10.1 analogue; the notation T_i^n)."""
    for k, lab in enumerate(("$t_{n-1}$", "$t_n$", "$t_{n+1}$")):
        ax.plot([0, nx - 1], [k, k], color=COLORS["muted"], lw=1)
        ax.plot(np.arange(nx), np.full(nx, k), "o", color=COLORS["ink"], ms=4)
        ax.text(-0.6, k, lab, ha="right", va="center")
    i = nx // 2
    for d, lab in ((-1, "$x_{i-1}$"), (0, "$x_i$"), (1, "$x_{i+1}$")):
        ax.text(i + d, -0.45, lab, ha="center", va="top")
    ax.text(0, -0.45, "$x_0=0$", ha="center", va="top")
    ax.text(nx - 1, -0.45, "$x_N=L$", ha="center", va="top")
    for d in (-1, 0, 1):  # the FTCS stencil: three nodes at t_n feed one node at t_{n+1}
        ax.plot([i + d, i], [1, 2], color=COLORS["orange"], lw=1.5, zorder=2)
        ax.plot(i + d, 1, "o", color=COLORS["orange"], ms=8, zorder=3)
    ax.plot(i, 2, "s", color=COLORS["orange"], ms=9, zorder=3)
    ax.text(i + 0.25, 2.2, "$T_i^{n+1}$", color=COLORS["orange"])
    ax.annotate("", (i + 1, 0.55), (i, 0.55), arrowprops=dict(arrowstyle="<->", color=COLORS["accent"]))
    ax.text(i + 0.5, 0.62, r"$\Delta x$", ha="center", color=COLORS["accent"])
    ax.set_xlim(-1.5, nx)
    ax.set_ylim(-1, 2.6)
    ax.set_title("space–time grid, $T_i^n\\approx T(x_i,t_n)$")
    _clean(ax)


def stencil_diagram(ax, kind: str = "ftcs"):
    """The points an update formula uses: FTCS (three at t_n → one at t_{n+1}), BTCS (three at t_{n+1}), upwind (two)."""
    pts_old = {"ftcs": (-1, 0, 1), "btcs": (0,), "upwind": (-1, 0), "lax_wendroff": (-1, 0, 1)}[kind]
    pts_new = (-1, 0, 1) if kind == "btcs" else (0,)
    for k in (0, 1):
        ax.plot([-2, 2], [k, k], color=COLORS["grid"], lw=1)
    col = dict(ftcs=COLORS["orange"], btcs=COLORS["blue"], upwind=COLORS["teal"], lax_wendroff=COLORS["accent"])[kind]
    for d in pts_old:
        ax.plot(d, 0, "o", color=col, ms=10)
        for e in pts_new:
            ax.plot([d, e], [0, 1], color=col, lw=1, alpha=0.6)
    for e in pts_new:
        ax.plot(e, 1, "s", color=col, ms=10)
    ax.text(2.2, 0, "$n$", va="center")
    ax.text(2.2, 1, "$n+1$", va="center")
    ax.set_xlim(-2.4, 2.8)
    ax.set_ylim(-0.5, 1.5)
    ax.set_title(f"{kind.upper()} stencil")
    _clean(ax)


def characteristic_diagram(ax, C: float = 0.8):
    """Upwind stencil (i − 1, i at t_n → i at t_{n+1}) and the characteristic x − ut through the new point (CFL: its foot
    must lie inside the stencil, C = uΔt/Δx ≤ 1)."""
    for k in (0, 1):
        ax.plot([-2.5, 1.5], [k, k], color=COLORS["grid"], lw=1)
        ax.plot([-2, -1, 0, 1], [k] * 4, "o", color=COLORS["muted"], ms=5)
    ax.plot([-1, 0], [0, 0], "o", color=COLORS["teal"], ms=10)
    ax.plot(0, 1, "s", color=COLORS["teal"], ms=10)
    ax.fill_between([-1, 0], [0, 0], [0.04, 0.04], color=COLORS["teal"], alpha=0.3)
    ok = C <= 1
    ax.plot([0, -C], [1, 0], color=COLORS["accent"] if ok else COLORS["rose"], lw=2)
    ax.plot(-C, 0, "x", color=COLORS["accent"] if ok else COLORS["rose"], ms=10, mew=2)
    ax.text(-C, -0.25, f"foot at $x_i-C\\Delta x$, C = {C:g}", ha="center", va="top", fontsize=9)
    ax.set_xlim(-2.6, 1.6)
    ax.set_ylim(-0.6, 1.4)
    ax.set_title("CFL: domain of dependence")
    _clean(ax)


def hat_functions(ax, n_el: int = 6):
    """Piecewise-linear hat functions N_0, N_A, N_n on a uniform mesh (our Fig. 10.2 analogue)."""
    from fluidpy.core import fem1d as FEM1

    nodes = np.linspace(0, 1, n_el + 1)
    x = np.linspace(0, 1, 601)
    N = FEM1.hat_basis(nodes, x)
    for A in range(n_el + 1):
        c = COLORS["accent"] if A in (0, n_el) else (COLORS["teal"] if A == n_el // 2 else COLORS["muted"])
        ax.plot(x, N[:, A], color=c, lw=2 if c != COLORS["muted"] else 1)
    ax.plot(nodes, 0 * nodes, "o", color=COLORS["ink"], ms=4)
    ax.text(0, 1.05, "$N_0$", ha="center")
    ax.text(nodes[n_el // 2], 1.05, "$N_A$", ha="center", color=COLORS["teal"])
    ax.text(1, 1.05, "$N_n$", ha="center")
    ax.set_title("hat functions, $N_A(x_B)=\\delta_{AB}$")
    ax.set_xlabel("x / L")
    ax.set_ylim(-0.1, 1.25)


def element_map(ax):
    """Element e = [x_{A−1}, x_A] and the parent element ξ ∈ [−1, 1] with N₁, N₂ (our Fig. 10.3 analogue)."""
    xi = np.linspace(-1, 1, 11)
    ax.plot([0, 2], [0, 0], color=COLORS["ink"])
    ax.plot([0, 2], [1, 0], color=COLORS["teal"])
    ax.plot([0, 2], [0, 1], color=COLORS["accent"])
    ax.text(0, -0.15, "$x_{A-1}$", ha="center", va="top")
    ax.text(2, -0.15, "$x_A$", ha="center", va="top")
    ax.text(0.2, 1.0, "$N_{A-1}$", color=COLORS["teal"])
    ax.text(1.5, 1.0, "$N_A$", color=COLORS["accent"])
    off = 3.9
    ax.plot(off + xi, 0 * xi, color=COLORS["ink"])
    ax.plot(off + xi, 0.5 * (1 - xi), color=COLORS["teal"])
    ax.plot(off + xi, 0.5 * (1 + xi), color=COLORS["accent"])
    ax.text(off - 1, -0.15, "$\\xi=-1$", ha="center", va="top")
    ax.text(off + 1, -0.15, "$\\xi=1$", ha="center", va="top")
    ax.text(off - 0.8, 1.0, "$N_1$", color=COLORS["teal"])
    ax.text(off + 0.6, 1.0, "$N_2$", color=COLORS["accent"])
    ax.annotate("", (off - 1.1, 0.5), (2.1, 0.5), arrowprops=dict(arrowstyle="<->", color=COLORS["muted"]))
    ax.text(1.6 + 0.5 * (off - 1.1 - 2.1) + 0.5, 0.58, "$x(\\xi)$", ha="center", color=COLORS["muted"])
    ax.set_xlim(-0.3, off + 1.4)
    ax.set_ylim(-0.5, 1.3)
    ax.set_title("element and parent element")
    _clean(ax)


def staggered_grid(ax, n: int = 3, highlight=(2, 2)):
    """Staggered (MAC) grid: p at cell centres, u on vertical faces, v on horizontal faces (our Fig. 10.4 analogue)."""
    for k in range(n + 1):
        ax.plot([0, n], [k, k], color=COLORS["grid"])
        ax.plot([k, k], [0, n], color=COLORS["grid"])
    ax.plot([0, n, n, 0, 0], [0, 0, n, n, 0], color=COLORS["ink"], lw=3)  # the boundary Γ: wall faces drawn thick
    hi, hj = highlight
    ax.add_patch(__import__("matplotlib.patches", fromlist=["Rectangle"]).Rectangle((hi - 1, hj - 1), 1, 1,
                                                                                   color=COLORS["amber"], alpha=0.15))
    for j in range(n):
        for i in range(n):
            ax.plot(i + 0.5, j + 0.5, "o", color=COLORS["orange"], ms=6)
    for j in range(n):
        for i in range(n + 1):
            ax.annotate("", (i + 0.18, j + 0.5), (i - 0.18, j + 0.5), arrowprops=dict(arrowstyle="->", color=COLORS["blue"]))
    for j in range(n + 1):
        for i in range(n):
            ax.annotate("", (i + 0.5, j + 0.18), (i + 0.5, j - 0.18), arrowprops=dict(arrowstyle="->", color=COLORS["teal"]))
    ax.text(hi - 0.5, hj - 0.38, "$p_{i,j}$", ha="center", color=COLORS["orange"])
    ax.text(hi + 0.05, hj - 0.62, "$u_{i+1/2,j}$", color=COLORS["blue"], fontsize=9)
    ax.text(hi - 0.45, hj + 0.08, "$v_{i,j+1/2}$", color=COLORS["teal"], fontsize=9)
    ax.set_xlim(-0.3, n + 0.3)
    ax.set_ylim(-0.3, n + 0.3)
    ax.set_title("staggered (MAC / C-) grid")
    _clean(ax)


def p2p1_triangle(ax, iso: bool = False):
    """The P2–P1 (Taylor–Hood) triangle: 6 velocity nodes, 3 pressure vertices (our Fig. 10.5a); ``iso``: iso-P2/P1 (10.5b)."""
    V = np.array([[0, 0], [1, 0], [0.5, np.sqrt(3) / 2]])
    ax.fill(*V.T, color=COLORS["grid"], alpha=0.5)
    ax.plot(*np.vstack([V, V[:1]]).T, color=COLORS["ink"])
    M = 0.5 * (V + np.roll(V, -1, axis=0))
    if iso:
        ax.plot(*np.vstack([M, M[:1]]).T, color=COLORS["ink"], lw=1)
    else:
        ax.plot(*M.T, "o", mfc="white", mec=COLORS["blue"], ms=8)
    ax.plot(*V.T, "o", color=COLORS["orange"], ms=9)
    ax.set_title("iso-P2/P1" if iso else "P2 velocity / P1 pressure")
    _clean(ax)


def _box(ax, x0, y0, w, h, **kw):
    import matplotlib.patches as mp

    ax.add_patch(mp.Rectangle((x0, y0), w, h, **kw))


def cavity_sketch(ax):
    """Lid-driven square cavity (our Fig. 10.6 analogue): side L, lid speed U."""
    _box(ax, 0, 0, 1, 1, fill=False, ec=COLORS["ink"], lw=2)
    ax.plot([0, 1], [1.04, 1.04], color=COLORS["accent"], lw=4)
    ax.annotate("", (0.7, 1.12), (0.3, 1.12), arrowprops=dict(arrowstyle="->", color=COLORS["accent"]))
    ax.text(0.5, 1.16, "U", ha="center", color=COLORS["accent"])
    ax.text(0.5, -0.08, "L", ha="center", va="top")
    ax.plot([0.5, 0.5], [0, 1], "--", color=COLORS["muted"], lw=1)  # the centreline x = L/2 of the benchmark
    ax.text(0.52, 0.03, "x = ½", color=COLORS["muted"], fontsize=8)
    t = np.linspace(0, 2 * np.pi, 100)  # eddy cartoons: the primary eddy and two small bottom-corner eddies
    ax.plot(0.62 + 0.28 * np.cos(t), 0.66 + 0.24 * np.sin(t), color=COLORS["blue"], lw=1.5)
    ax.annotate("", (0.62 + 0.03, 0.90), (0.62 - 0.03, 0.90), arrowprops=dict(arrowstyle="->", color=COLORS["blue"]))
    ax.text(0.62, 0.66, "primary\neddy", ha="center", va="center", fontsize=8, color=COLORS["blue"])
    for cx in (0.06, 0.92):
        ax.plot(cx + 0.045 * np.cos(t), 0.06 + 0.045 * np.sin(t), color=COLORS["teal"], lw=1)
    ax.set_xlim(-0.2, 1.2)
    ax.set_ylim(-0.2, 1.3)
    ax.set_title("lid-driven cavity")
    _clean(ax)


def block_channel_sketch(ax, H: float = 4.0, ahead: float = 8.0, behind: float = 20.0):
    """Square block between two plates sliding at U, uniform inflow U (our Fig. 10.9 analogue, our geometry)."""
    Lx = ahead + 1 + behind
    ax.plot([0, Lx], [0, 0], color=COLORS["ink"], lw=3)
    ax.plot([0, Lx], [H, H], color=COLORS["ink"], lw=3)
    _box(ax, ahead, 0.5 * (H - 1), 1, 1, color=COLORS["muted"])
    for y in np.linspace(0.4, H - 0.4, 5):
        ax.annotate("", (1.5, y), (0, y), arrowprops=dict(arrowstyle="->", color=COLORS["blue"]))
    ax.annotate("", (ahead + 4, -0.35), (ahead, -0.35), arrowprops=dict(arrowstyle="->", color=COLORS["accent"]))
    ax.annotate("", (ahead + 4, H + 0.35), (ahead, H + 0.35), arrowprops=dict(arrowstyle="->", color=COLORS["accent"]))
    ax.text(ahead + 2, H + 0.5, "plates slide at U", ha="center", fontsize=8)
    ax.set_xlim(-0.5, Lx + 0.5)
    ax.set_ylim(-1, H + 1)
    ax.set_title(f"block in a channel (ours: H = {H:g}, {ahead:g} ahead, {behind:g} behind)")
    _clean(ax)


def cylinder_channel_sketch(ax, W: float = 5.0, xmin: float = -8.0, xmax: float = 16.0):
    """Cylinder in a channel with sliding walls; boundaries Γ₁ … Γ₅ (our Fig. 10.15 analogue, our geometry)."""
    import matplotlib.patches as mp

    ax.plot([xmin, xmax], [-W / 2, -W / 2], color=COLORS["ink"], lw=3)
    ax.plot([xmin, xmax], [W / 2, W / 2], color=COLORS["ink"], lw=3)
    ax.plot([xmin, xmin], [-W / 2, W / 2], color=COLORS["blue"], lw=1, ls="--")
    ax.plot([xmax, xmax], [-W / 2, W / 2], color=COLORS["muted"], lw=1, ls="--")
    ax.add_patch(mp.Circle((0, 0), 0.5, color=COLORS["muted"]))
    for lab, (x, y) in {r"$\Gamma_1$": (xmin - 0.8, 0), r"$\Gamma_2$": (xmax + 0.3, 0), r"$\Gamma_3$": (4, W / 2 + 0.4),
                        r"$\Gamma_4$": (4, -W / 2 - 0.8), r"$\Gamma_5$": (0.6, 0.6)}.items():
        ax.text(x, y, lab)
    ax.set_xlim(xmin - 1.5, xmax + 1.5)
    ax.set_ylim(-W / 2 - 1.2, W / 2 + 1.2)
    ax.set_title(f"cylinder in a channel (ours: W = {W:g}d, x from {xmin:g}d to {xmax:g}d)")
    _clean(ax)


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(4, 3, figsize=(14, 13))
    spacetime_grid(axs[0, 0])
    stencil_diagram(axs[0, 1], "ftcs")
    characteristic_diagram(axs[0, 2], 0.8)
    hat_functions(axs[1, 0])
    element_map(axs[1, 1])
    staggered_grid(axs[1, 2])
    p2p1_triangle(axs[2, 0])
    p2p1_triangle(axs[2, 1], iso=True)
    cavity_sketch(axs[2, 2])
    block_channel_sketch(axs[3, 0])
    cylinder_channel_sketch(axs[3, 1])
    stencil_diagram(axs[3, 2], "btcs")
    save(fig, out, "drawings")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
