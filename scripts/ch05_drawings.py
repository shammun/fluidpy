"""Small drawing helpers for the ch05 notebook and scripts (no physics: the numbers come from
``fluidpy.ch05_vorticity_dynamics``). Analogues of the book's sketches drawn by our own code (Figs. 5.1, 5.2, 5.4,
5.6, 5.8, 5.9, 5.10, 5.11–5.14, 5.16), never copies.

    from scripts.ch05_drawings import tube_mesh, tank_section, disc_element, loop_with_element, \
        biot_savart_geometry, helix_with_frame, column_sketch, vortex_marks, sheet_circuit

Run as a script it draws a demo sheet of every helper:
``.venv/Scripts/python.exe scripts/ch05_drawings.py --no-show`` → outputs/ch05/drawings_demo.png (+ the plotly helix
as outputs/ch05/helix_frame.html).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]
for _p in (_ROOT, _ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from fluidpy import ch05_vorticity_dynamics as ch05  # noqa: E402
from fluidpy.core.style import COLORS  # noqa: E402

POS_C, NEG_C = COLORS["teal"], COLORS["rose"]  # counterclockwise (+) and clockwise (−) vortices
P_C, RHO_C = COLORS["orange"], COLORS["blue"]  # isobars, isopycnals
STRETCH_C, TILT_C = COLORS["accent"], COLORS["blue"]


def tube_mesh(ax3d, tube, color=POS_C, every: int = 1, sections=(0, -1)):
    """Draw a traced vortex tube (``ch05.vortex_tube``): its vortex lines and the closed cross-sections ``sections``
    (indices along the lines) — the Fig. 5.1 analogue. Returns the list of artists."""
    L = tube.array
    arts = []
    for j in range(0, L.shape[0], every):
        arts += ax3d.plot(L[j, 0], L[j, 1], L[j, 2], color=color, lw=0.9, alpha=0.85)
    for k in sections:
        P = L[:, :, k]
        P = np.vstack([P, P[:1]])
        arts += ax3d.plot(P[:, 0], P[:, 1], P[:, 2], color=COLORS["ink"], lw=1.6)
    ax3d.set_xlabel("x [m]")
    ax3d.set_ylabel("y [m]")
    ax3d.set_zlabel("z [m]")
    return arts


def tank_section(ax, R: float, H: float, label: str | None = None):
    """Walls and bottom of a cylindrical tank of radius R and height H in the (r, z) plane, mirrored about the axis
    (both sides −R … R), with the axis dashed. Returns the artists."""
    arts = ax.plot([-R, -R, R, R], [H, 0.0, 0.0, H], color=COLORS["ink"], lw=2.0)
    arts += ax.plot([0.0, 0.0], [0.0, H], ls="--", color=COLORS["muted"], lw=0.8)
    if label:
        arts.append(ax.text(0.0, H * 1.02, label, ha="center", va="bottom", color=COLORS["ink"], fontsize=9))
    ax.set_xlabel("r [m]")
    ax.set_ylabel("z [m]")
    return arts


def disc_element(ax, center, radius: float, isobars: float = 0.0, isopycnals: float = 0.5, grad_rho: float = 5.0,
                 rho0: float = 1000.0, g: float = 9.81, n_arrows: int = 16):
    """The Fig. 5.6 analogue: a fluid disc with isobars (orange, horizontal when ``isobars`` = 0 rad) and isopycnals
    (blue, tilted by ``isopycnals`` rad), inward pressure arrows round the rim, the geometric centre, the centre of mass
    G (the true offset R²∇ρ/4ρ₀ is tiny; G is drawn at 0.35R along it) and a curved arrow showing the sense of the
    torque (from ``ch05.baroclinic_element_scenario``). Returns the scenario dict."""
    c = np.asarray(center, float)
    sc = ch05.baroclinic_element_scenario(isopycnals - isobars, grad_rho, radius, rho0, g)
    th = np.linspace(0, 2 * np.pi, 200)
    ax.plot(c[0] + radius * np.cos(th), c[1] + radius * np.sin(th), color=COLORS["ink"], lw=1.8)
    for off in np.linspace(-0.8, 0.8, 5) * radius:
        for ang, col, ls in ((isobars, P_C, "-"), (isopycnals, RHO_C, "--")):
            d = np.array([np.cos(ang), np.sin(ang)])
            nrm = np.array([-d[1], d[0]])
            half = np.sqrt(max(radius ** 2 - off ** 2, 0.0))
            p0, p1 = c + off * nrm - half * d, c + off * nrm + half * d
            ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=col, ls=ls, lw=1.0)
    for a in np.linspace(0, 2 * np.pi, n_arrows, endpoint=False):
        e = np.array([np.cos(a), np.sin(a)])
        ax.annotate("", c + 1.0 * radius * e, c + 1.35 * radius * e,
                    arrowprops=dict(arrowstyle="-|>", color=COLORS["muted"], lw=0.8))
    ax.plot(*c, "o", color=COLORS["ink"], ms=4)
    off = sc["offset"]
    G = c + 0.35 * radius * off / max(np.linalg.norm(off), 1e-300)
    ax.plot(*G, "o", color=COLORS["amber"], ms=6)
    ax.text(G[0] + 0.05 * radius, G[1] + 0.05 * radius, "G", color=COLORS["amber"], fontsize=9)
    if sc["sense"] != "none":
        sgn = 1.0 if sc["sense"] == "counterclockwise" else -1.0
        arc = np.linspace(0.3, 1.3, 30) * sgn
        pts = c[:, None] + 0.45 * radius * np.stack([np.cos(arc + np.pi / 2), np.sin(arc + np.pi / 2)])
        ax.plot(pts[0], pts[1], color=COLORS["rose"], lw=1.6)
        ax.annotate("", pts[:, -1], pts[:, -3], arrowprops=dict(arrowstyle="-|>", color=COLORS["rose"], lw=1.6))
    ax.set_aspect("equal")
    ax.axis("off")
    return sc


def loop_with_element(ax, pts, k: int = 0, u=None, t: float = 0.0, scale: float = 0.3):
    """Draw a closed loop (d, N) and, at point k, the element dx (to the next point, lengthened for visibility) and —
    if a velocity field ``u`` is given — the velocities u and u + du at its two ends (the Fig. 5.4 analogue)."""
    P = np.asarray(pts, float)
    Pc = np.hstack([P, P[:, :1]])
    ax.plot(Pc[0], Pc[1], color=COLORS["accent"], lw=1.6)
    j = (k + max(P.shape[1] // 40, 1)) % P.shape[1]
    ax.annotate("", P[:2, j], P[:2, k], arrowprops=dict(arrowstyle="-|>", color=COLORS["ink"], lw=1.4))
    ax.text(*(0.5 * (P[:2, j] + P[:2, k])), "  dx", color=COLORS["ink"], fontsize=9)
    if u is not None:
        for idx, lab, col in ((k, "u", COLORS["teal"]), (j, "u + du", COLORS["orange"])):
            v = np.asarray(u(P[:, idx], t), float)[:2]
            ax.annotate("", P[:2, idx] + scale * v, P[:2, idx], arrowprops=dict(arrowstyle="-|>", color=col, lw=1.3))
            ax.text(*(P[:2, idx] + scale * v), lab, color=col, fontsize=9)
    ax.set_aspect("equal")
    span = np.ptp(P[:2], axis=1).max()
    ax.set_xlim(P[0].min() - 0.4 * span, P[0].max() + 0.4 * span)
    ax.set_ylim(P[1].min() - 0.4 * span, P[1].max() + 0.4 * span)
    return j


def biot_savart_geometry(ax3d, polyline, x, Gamma: float = 1.0, closed: bool = False, k: int | None = None,
                         scale: float = 1.0):
    """The Fig. 5.8 analogue: a filament (3, M), a field point x, the vector x − x′ from one segment's midpoint
    (index k, default the middle) and that segment's contribution du (5.17) plus the total u at x (scaled).
    Returns (du_k, u) [m/s]."""
    P = np.asarray(polyline, float)
    x = np.asarray(x, float)
    ax3d.plot(P[0], P[1], P[2], color=POS_C, lw=2.0)
    contrib = ch05.filament_contributions(x, P, Gamma, closed)
    k = contrib.shape[1] // 2 if k is None else k
    a, b = P[:, k], P[:, (k + 1) % P.shape[1]]
    mid = 0.5 * (a + b)
    ax3d.plot([mid[0], x[0]], [mid[1], x[1]], [mid[2], x[2]], ls="--", color=COLORS["muted"])
    du = contrib[:, k]
    u = contrib.sum(axis=1)
    for vec, col in ((du / max(np.linalg.norm(du), 1e-300) * 0.3 * scale, COLORS["orange"]),
                     (u / max(np.linalg.norm(u), 1e-300) * 0.5 * scale, COLORS["accent"])):
        ax3d.quiver(*x, *vec, color=col, lw=1.6)
    ax3d.scatter(*x, color=COLORS["ink"], s=12)
    return du, u


def helix_with_frame(fig=None, s0: float = 0.0, a: float = 1.0, c: float = 0.3, turns: float = 2.0,
                     length: float = 0.5):
    """Plotly 3-D figure (the Fig. 5.9 analogue): a helical vortex line (``ch05.helix``) with the natural frame e_s,
    e_n (away from the centre of curvature), e_m (``ch05.helix_frame``) drawn at arc length s₀. Returns the figure."""
    import plotly.graph_objects as go

    fig = go.Figure() if fig is None else fig
    k = np.sqrt(a ** 2 + c ** 2)
    s = np.linspace(0.0, 2 * np.pi * turns * k, 400)
    X = ch05.helix(s, a, c)
    fig.add_trace(go.Scatter3d(x=X[0], y=X[1], z=X[2], mode="lines", line=dict(color=POS_C, width=5),
                               name="vortex line"))
    P = ch05.helix(s0, a, c)
    fr = ch05.helix_frame(s0, a, c)
    for key, col in (("e_s", COLORS["accent"]), ("e_n", COLORS["orange"]), ("e_m", COLORS["blue"])):
        e = fr[key] * length
        fig.add_trace(go.Scatter3d(x=[P[0], P[0] + e[0]], y=[P[1], P[1] + e[1]], z=[P[2], P[2] + e[2]],
                                   mode="lines+text", text=["", key], line=dict(color=col, width=7), name=key))
    fig.update_layout(scene=dict(aspectmode="data", xaxis_title="x [m]", yaxis_title="y [m]", zaxis_title="z [m]"),
                      title=f"natural frame on a helix: κ = {fr['curvature']:.4f} 1/m, τ = {fr['torsion']:.4f} 1/m",
                      margin=dict(l=0, r=0, t=40, b=0))
    return fig


def column_sketch(ax, x, h, zeta, n_cols: int = 5, f: float | None = None):
    """The Fig. 5.10 analogue: the bottom h(x) of a rotating layer (surface at h = max), with ``n_cols`` fluid columns
    whose spin arrows are sized by |ζ| and coloured teal (cyclonic, ζ > 0) or rose (anticyclonic)."""
    x, h, zeta = np.asarray(x, float), np.asarray(h, float), np.asarray(zeta, float)
    top = np.full_like(x, h.max())
    ax.fill_between(x, top - h, top - h.max() * 1.1, color=COLORS["grid"])
    ax.plot(x, top - h, color=COLORS["ink"], lw=1.5)
    ax.plot(x, top, color=COLORS["blue"], lw=1.0)
    zmax = max(np.abs(zeta).max(), 1e-300)
    for i in np.linspace(0, x.size - 1, n_cols).astype(int):
        xb, zb = x[i], top[i] - h[i]
        w = 0.03 * (x[-1] - x[0])
        ax.add_patch(__import__("matplotlib").patches.Rectangle((xb - w / 2, zb), w, h[i], fill=False,
                                                                 ec=COLORS["muted"], lw=1.0))
        col = POS_C if zeta[i] >= 0 else NEG_C
        ax.text(xb, top[i] + 0.04 * h.max(), "↺" if zeta[i] >= 0 else "↻", ha="center", color=col,
                fontsize=8 + 14 * abs(zeta[i]) / zmax)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("height [m]")


def vortex_marks(ax, xv, Gamma, size: float = 0.08, trails=None):
    """Point vortices as circular-arrow marks, teal for counterclockwise Γ > 0, rose for clockwise, with optional
    trails (T, 2, M). Returns the artists."""
    V = np.asarray(xv, float).reshape(2, -1)
    arts = []
    if trails is not None:
        Tr = np.asarray(trails, float)
        for m in range(V.shape[1]):
            arts += ax.plot(Tr[:, 0, m], Tr[:, 1, m], lw=0.8, color=POS_C if Gamma[m] > 0 else NEG_C, alpha=0.6)
    th = np.linspace(0.2, 2 * np.pi - 0.3, 40)
    for m in range(V.shape[1]):
        col = POS_C if Gamma[m] > 0 else NEG_C
        sgn = 1.0 if Gamma[m] > 0 else -1.0
        arc = V[:, m, None] + size * np.stack([np.cos(sgn * th), np.sin(sgn * th)])
        arts += ax.plot(arc[0], arc[1], color=col, lw=1.4)
        arts.append(ax.annotate("", arc[:, -1], arc[:, -3], arrowprops=dict(arrowstyle="-|>", color=col, lw=1.2)))
        arts += ax.plot(*V[:, m], "o", color=col, ms=3)
    ax.set_aspect("equal")
    return arts


def sheet_circuit(ax, ds: float = 1.0, dn: float = 0.4, n_filaments: int = 9, gamma: float = -1.0):
    """The Fig. 5.16 analogue: a row of filaments (marks by the sign of γ) on y = 0, the dn × ds circuit traversed
    counterclockwise, and the tangential velocities u₁ (above) and u₂ (below) of the sheet (∓γ/2 for γ counterclockwise
    positive). Returns (u1, u2, γ = u₂ − u₁)."""
    xs = np.linspace(-0.8 * ds, 0.8 * ds, n_filaments)
    vortex_marks(ax, np.stack([xs, 0 * xs]), np.full(n_filaments, gamma), size=0.07 * ds)
    X = [-0.5 * ds, 0.5 * ds, 0.5 * ds, -0.5 * ds, -0.5 * ds]
    Y = [-0.5 * dn, -0.5 * dn, 0.5 * dn, 0.5 * dn, -0.5 * dn]
    ax.plot(X, Y, color=COLORS["accent"], lw=1.4)
    for (x0, y0, x1, y1) in ((-0.1, -0.5, 0.1, -0.5), (0.5, -0.1, 0.5, 0.1), (0.1, 0.5, -0.1, 0.5), (-0.5, 0.1, -0.5, -0.1)):
        ax.annotate("", (x1 * ds, y1 * dn), (x0 * ds, y0 * dn), arrowprops=dict(arrowstyle="-|>", color=COLORS["accent"]))
    u1, u2 = -0.5 * gamma, 0.5 * gamma
    for y, u, lab in ((0.75 * dn, u1, "u₁"), (-0.75 * dn, u2, "u₂")):
        ax.annotate("", (0.3 * ds * np.sign(u), y), (0.0, y), arrowprops=dict(arrowstyle="-|>", color=COLORS["orange"]))
        ax.text(0.32 * ds * np.sign(u), y, f"{lab} = {u:+.2f}", color=COLORS["orange"], va="center", fontsize=9,
                ha="left" if u > 0 else "right")
    ax.set_aspect("equal")
    ax.set_xlim(-1.0 * ds, 1.0 * ds)
    ax.set_ylim(-1.2 * dn, 1.6 * dn)
    ax.axis("off")
    return u1, u2, float(ch05.vortex_sheet_strength(u1, u2))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(_ROOT / "outputs" / "ch05"))
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
    fig = plt.figure(figsize=(16, 8))
    fig.set_layout_engine("none")  # 3-D axes and constrained layout do not mix
    ax = fig.add_subplot(2, 4, 1, projection="3d")
    tube = ch05.vortex_tube("gaussian_tube", (0, 0, 0.0), (0, 0, 1), 0.08, n_lines=16, s_max=0.6, n=41, both=False,
                            twist=3.0)
    tube_mesh(ax, tube)
    ax.set_title("tube_mesh (narrowing, twisted tube)", fontsize=9)
    ax = fig.add_subplot(2, 4, 2)
    tank_section(ax, 0.1, 0.2, "tank")
    r = np.linspace(-0.1, 0.1, 101)
    ts = ch05.rotating_tank_free_surface(0.1, 0.1, 5.0)
    ax.plot(r, ts["z_vertex"] + 25.0 * r ** 2 / (2 * 9.81), color=COLORS["blue"])
    ax.set_title("tank_section + free surface (Ω = 5 rad/s)", fontsize=9)
    ax = fig.add_subplot(2, 4, 3)
    sc = disc_element(ax, (0, 0), 1.0, 0.0, 0.5)
    ax.set_title(f"disc_element: {sc['sense']}", fontsize=9)
    ax = fig.add_subplot(2, 4, 4)
    s = ch05.kelvin_scenario("cellular")
    loop_with_element(ax, s["pts0"], 10, s["u"])
    ax.set_title("loop_with_element (cellular flow)", fontsize=9)
    ax = fig.add_subplot(2, 4, 5, projection="3d")
    du, u = biot_savart_geometry(ax, ch05.filament_preset("helix"), [0.8, 0.0, 0.2])
    ax.set_title(f"biot_savart_geometry |u| = {np.linalg.norm(u):.3f} m/s", fontsize=9)
    ax = fig.add_subplot(2, 4, 6)
    cs = ch05.column_over_slope(np.linspace(-3e5, 3e5, 121), "ridge")
    column_sketch(ax, cs["x"], cs["h"], cs["zeta"])
    ax.set_title("column_sketch over a ridge", fontsize=9)
    ax = fig.add_subplot(2, 4, 7)
    pr = ch05.point_vortex_preset("unequal_pair")
    tr = ch05.point_vortex_evolve(pr["xv"], pr["Gamma"], np.linspace(0, pr["t_end"] / 2, 80))
    vortex_marks(ax, tr[-1], pr["Gamma"], trails=tr)
    ax.set_title("vortex_marks (Γ₂ = 3Γ₁)", fontsize=9)
    ax = fig.add_subplot(2, 4, 8)
    u1, u2, gm = sheet_circuit(ax)
    ax.set_title(f"sheet_circuit: γ = u₂ − u₁ = {gm:+.2f} m/s", fontsize=9)
    fig.savefig(out / "drawings_demo.png", bbox_inches="tight")
    helix_with_frame(s0=1.0).write_html(str(out / "helix_frame.html"), include_plotlyjs="cdn")
    print(f"drawings demo → {out / 'drawings_demo.png'}; helix → {out / 'helix_frame.html'}")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
