"""Small matplotlib drawing helpers for the ch04 notebook and scripts (no physics: the numbers come from
``fluidpy.ch04_conservation_laws``). Analogues of the book's sketches drawn by our own code (Figs. 4.2–4.3, 4.6–4.7,
4.18–4.19), never copies.

    from scripts.ch04_drawings import cv_box, face_flux_arrows, budget_bars, stream_tube_element, rotating_frames, \
        cone_sweep, pillbox, curved_cap

Run as a script it draws a demo sheet of every helper:
``.venv/Scripts/python.exe scripts/ch04_drawings.py --no-show`` → outputs/ch04/drawings_demo.png.
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

from fluidpy import ch04_conservation_laws as ch04  # noqa: E402
from fluidpy.core.style import COLORS  # noqa: E402

IN_C, OUT_C, STORE_C, FORCE_C = COLORS["blue"], COLORS["orange"], COLORS["accent"], COLORS["rose"]


def cv_box(ax, x0: float, x1: float, y0: float, y1: float, b: float = 0.0, label: str = "CV", ticks: int = 4):
    """A control volume drawn as a dashed rectangle with short outward-normal ticks on each face; ``b`` ≠ 0 adds an
    arrow showing the surface velocity (the CV moving at b along x). Returns the list of artists."""
    arts = [ax.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], ls="--", color=COLORS["ink"], lw=1.2)[0]]
    L = 0.06 * max(x1 - x0, y1 - y0)
    for s in np.linspace(0.15, 0.85, ticks):
        xs, ys = x0 + s * (x1 - x0), y0 + s * (y1 - y0)
        for (px, py, nx, ny) in ((xs, y1, 0, 1), (xs, y0, 0, -1), (x1, ys, 1, 0), (x0, ys, -1, 0)):
            arts.append(ax.annotate("", (px + nx * L, py + ny * L), (px, py),
                                    arrowprops=dict(arrowstyle="->", color=COLORS["muted"], lw=0.8)))
    arts.append(ax.text(x0 + 0.02 * (x1 - x0), y1 - 0.08 * (y1 - y0), label, color=COLORS["ink"], fontsize=9))
    if b:
        xm, ym = 0.5 * (x0 + x1), y1 + 0.12 * (y1 - y0)
        arts.append(ax.annotate("", (xm + np.sign(b) * 0.25 * (x1 - x0), ym), (xm, ym),
                                arrowprops=dict(arrowstyle="-|>", color=STORE_C, lw=1.5)))
        arts.append(ax.text(xm, ym + 0.04 * (y1 - y0), f"b = {b:+.2f} m/s", ha="center", color=STORE_C, fontsize=9))
    return arts


def face_flux_arrows(ax, faces, colors=None, positions=None, scale: float = 1.0, key: str = "mass_flux"):
    """Arrows for per-face fluxes of ``ch04.cv_scenario(...)["faces"]``: inflow (negative) blue, outflow orange, length ∝
    |flux| × ``scale``. ``positions`` maps face name → ((x, y), (nx, ny)) outward normal; default puts inlet left, outlet
    right, top, bottom around the unit square. Returns the artists."""
    colors = colors or {"in": IN_C, "out": OUT_C}
    positions = positions or {"inlet": ((0.0, 0.5), (-1, 0)), "outlet": ((1.0, 0.5), (1, 0)), "top": ((0.5, 1.0), (0, 1)),
                              "bottom": ((0.5, 0.0), (0, -1)), "jet": ((0.0, 0.5), (-1, 0)),
                              "sheet_up": ((0.5, 1.0), (0, 1)), "sheet_down": ((0.5, 0.0), (0, -1)),
                              "nozzle": ((0.5, 0.0), (0, -1))}
    fmax = max(abs(f[key]) for f in faces) or 1.0
    arts = []
    for f in faces:
        if f["name"] not in positions:
            continue
        (px, py), (nx, ny) = positions[f["name"]]
        L = 0.3 * scale * abs(f[key]) / fmax
        out = f[key] > 0
        start = (px, py) if out else (px + nx * L, py + ny * L)
        end = (px + nx * L, py + ny * L) if out else (px, py)
        arts.append(ax.annotate("", end, start, arrowprops=dict(arrowstyle="-|>", lw=2,
                                                                color=colors["out" if out else "in"])))
        arts.append(ax.text(px + nx * (L + 0.05), py + ny * (L + 0.05), f"{f['name']}\n{f[key]:+.3g}", fontsize=8,
                            ha="center", va="center"))
    return arts


def budget_bars(ax, items: dict, total_label: str = "residual"):
    """Waterfall bars: each item (label → value) starts where the previous ended; negatives hatched; a final bar shows
    the running sum (the residual of the budget, ≈ 0). Returns the final sum."""
    run = 0.0
    labels = list(items) + [total_label]
    for i, (k, v) in enumerate(items.items()):
        ax.bar(i, v, bottom=run, color=OUT_C if v >= 0 else IN_C, hatch=None if v >= 0 else "//", edgecolor="white")
        run += v
    ax.bar(len(items), run, color=COLORS["ink"])
    ax.axhline(0, color=COLORS["muted"], lw=0.8)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    return run


def stream_tube_element(ax3d, A0: float = 1.0, A1: float = 1.6, L: float = 3.0, theta: float = 0.35, n: int = 40):
    """A stream-tube element of length L inclined at θ to the horizontal with end areas A0, A1 (Fig. 4.3 analogue), the
    extra pressure force on the conical wall drawn as a rose arrow and gravity as a blue arrow. Returns None."""
    r0, r1 = np.sqrt(A0 / np.pi), np.sqrt(A1 / np.pi)
    s = np.linspace(0, 1, 2)
    ph = np.linspace(0, 2 * np.pi, n)
    S, PH = np.meshgrid(s, ph)
    R = r0 + (r1 - r0) * S
    axis = np.array([np.cos(theta), 0.0, np.sin(theta)])
    e1 = np.array([-np.sin(theta), 0.0, np.cos(theta)])
    e2 = np.array([0.0, 1.0, 0.0])
    P = (S[..., None] * L * axis + R[..., None] * (np.cos(PH)[..., None] * e1 + np.sin(PH)[..., None] * e2))
    ax3d.plot_surface(P[..., 0], P[..., 1], P[..., 2], color=COLORS["teal"], alpha=0.25, linewidth=0)
    for k in (0, 1):
        ax3d.plot(P[:, k, 0], P[:, k, 1], P[:, k, 2], color=COLORS["teal"])
    ax3d.quiver(*(0.5 * L * axis + 0.5 * (r0 + r1) * e1), *(0.6 * axis), color=FORCE_C)
    ax3d.quiver(*(0.5 * L * axis), 0, 0, -0.8, color=IN_C)
    ax3d.set_xlabel("x")
    ax3d.set_zlabel("z")


def rotating_frames(ax, t: float, Omega: float = 0.5, X=(1.5, 0.8), size: float = 0.8):
    """Inertial axes O12 (grey) and a frame O′1′2′ at X turned by Ωt (purple) with the position vector x′ of a point P
    (Fig. 4.6 analogue, plane view). Returns the rotation angle Ωt [rad]."""
    th = Omega * t
    ax.annotate("", (size, 0), (0, 0), arrowprops=dict(arrowstyle="->", color=COLORS["muted"]))
    ax.annotate("", (0, size), (0, 0), arrowprops=dict(arrowstyle="->", color=COLORS["muted"]))
    Xo = np.asarray(X, float)
    for v, lab in ((np.array([np.cos(th), np.sin(th)]), "1′"), (np.array([-np.sin(th), np.cos(th)]), "2′")):
        ax.annotate("", Xo + size * v, Xo, arrowprops=dict(arrowstyle="->", color=STORE_C, lw=1.5))
        ax.text(*(Xo + 1.1 * size * v), lab, color=STORE_C)
    ax.annotate("", Xo, (0, 0), arrowprops=dict(arrowstyle="->", color=COLORS["ink"], ls="--"))
    ax.text(*(0.5 * Xo + np.array([0.02, 0.08])), "X(t)", fontsize=9)
    P = Xo + np.array([0.6, -0.5])
    ax.plot(*P, "o", color=COLORS["orange"])
    ax.annotate("", P, Xo, arrowprops=dict(arrowstyle="->", color=COLORS["orange"]))
    ax.text(*(P + 0.05), "P  (x′)", fontsize=9)
    ax.set_aspect("equal")
    return th


def cone_sweep(ax3d, alpha: float = 0.6, Omega=(0.0, 0.0, 1.0), dt: float = 0.4, n: int = 60):
    """The unit vector e′₁ at angle α from Ω sweeping a cone (Fig. 4.7 analogue): the circle traced, e′₁ at t and t + dt,
    and de′₁ = Ω × e′₁ dt. Returns |de′₁|/dt = sin α |Ω|."""
    Om = np.asarray(Omega, float)
    e0 = np.array([np.sin(alpha), 0.0, np.cos(alpha)])
    ts = np.linspace(0, 2 * np.pi / np.linalg.norm(Om), n)
    ring = np.array([ch04.rotating_basis(Om, t, np.vstack([e0, [0, 1, 0], [0, 0, 1]]))[0] for t in ts])
    ax3d.plot(*ring.T, color=COLORS["grid"])
    e1 = ch04.rotating_basis(Om, dt, np.vstack([e0, [0, 1, 0], [0, 0, 1]]))[0]
    for v, c in ((e0, STORE_C), (e1, COLORS["teal"])):
        ax3d.quiver(0, 0, 0, *v, color=c)
    de = np.cross(Om, e0) * dt
    ax3d.quiver(*e0, *de, color=COLORS["orange"])
    ax3d.quiver(0, 0, 0, *(1.3 * Om / np.linalg.norm(Om)), color=COLORS["ink"])
    return float(np.linalg.norm(np.cross(Om, e0)))


def pillbox(ax, l: float = 0.3, width: float = 1.2, label: bool = True):
    """A pillbox of height l straddling an interface (Fig. 4.18 analogue): the interface line, the box, ±n dA on the end
    faces. Returns the side-wall height l (the term that vanishes as l → 0)."""
    ax.plot([-1, 1], [0, 0], color=COLORS["ink"], lw=1.5)
    ax.fill_between([-width / 2, width / 2], -l / 2, l / 2, color=COLORS["teal"], alpha=0.2)
    ax.plot([-width / 2, width / 2, width / 2, -width / 2, -width / 2], [-l / 2, -l / 2, l / 2, l / 2, -l / 2],
            color=COLORS["teal"])
    ax.annotate("", (0, l / 2 + 0.3), (0, l / 2), arrowprops=dict(arrowstyle="-|>", color=OUT_C))
    ax.annotate("", (0, -l / 2 - 0.3), (0, -l / 2), arrowprops=dict(arrowstyle="-|>", color=OUT_C))
    if label:
        ax.text(0.05, l / 2 + 0.3, "+n dA (medium 1)", fontsize=8)
        ax.text(0.05, -l / 2 - 0.35, "−n dA (medium 2)", fontsize=8)
        ax.text(width / 2 + 0.05, 0.02, f"l = {l:.2f}", fontsize=8)
    ax.set_aspect("equal")
    return l


def curved_cap(ax3d, R1: float = 1.0, R2: float = 2.0, zeta: float = 0.2, n: int = 24, sigma: float = 1.0):
    """The cap z = x²/2R₁ + y²/2R₂ below z = ζ with its rim C and the surface-tension pull t × n round it (Fig. 4.19
    analogue). Returns the net upward pull σ∮(t × n)_z ds (``ch04.cap_surface_tension_force``)."""
    a, b = np.sqrt(2 * R1 * zeta), np.sqrt(2 * R2 * zeta)
    rr, th = np.meshgrid(np.linspace(0, 1, 12), np.linspace(0, 2 * np.pi, 48))
    X, Y = a * rr * np.cos(th), b * rr * np.sin(th)
    ax3d.plot_surface(X, Y, X ** 2 / (2 * R1) + Y ** 2 / (2 * R2), color=COLORS["teal"], alpha=0.3, linewidth=0)
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    x, y = a * np.cos(t), b * np.sin(t)
    tv = np.stack([-a * np.sin(t), b * np.cos(t), 0 * t])
    tv /= np.linalg.norm(tv, axis=0)
    nv = np.stack([-x / R1, -y / R2, np.ones_like(t)])
    nv /= np.linalg.norm(nv, axis=0)
    f = np.cross(tv, nv, axis=0)
    ax3d.plot(np.append(x, x[0]), np.append(y, y[0]), np.full(n + 1, zeta), color=COLORS["ink"])
    ax3d.quiver(x, y, np.full(n, zeta), *(0.3 * f), color=FORCE_C)
    return ch04.cap_surface_tension_force(sigma, R1, R2, zeta)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(_ROOT / "outputs" / "ch04"))
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
    fig = plt.figure(figsize=(14, 8))
    fig.set_layout_engine("none")  # 3-D axes and constrained layout do not mix
    ax = fig.add_subplot(2, 4, 1)
    cv_box(ax, 0, 1, 0, 1, b=-0.5, label="V*")
    sc = ch04.cv_scenario("wake")
    face_flux_arrows(ax, sc["faces"])
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(-0.6, 1.6)
    ax.set_title("cv_box + face_flux_arrows (wake)", fontsize=9)
    ax = fig.add_subplot(2, 4, 2)
    res = budget_bars(ax, {f["name"]: f["mass_flux"] for f in sc["faces"]})
    ax.set_title(f"budget_bars: mass, residual {res:.1e}", fontsize=9)
    ax = fig.add_subplot(2, 4, 3, projection="3d")
    stream_tube_element(ax)
    ax.set_title("stream_tube_element", fontsize=9)
    ax = fig.add_subplot(2, 4, 4)
    rotating_frames(ax, t=1.0)
    ax.set_title("rotating_frames", fontsize=9)
    ax = fig.add_subplot(2, 4, 5, projection="3d")
    r = cone_sweep(ax)
    ax.set_title(f"cone_sweep |de′/dt| = {r:.3f}", fontsize=9)
    ax = fig.add_subplot(2, 4, 6)
    pillbox(ax)
    ax.set_title("pillbox", fontsize=9)
    ax = fig.add_subplot(2, 4, 7, projection="3d")
    F = curved_cap(ax)
    ax.set_title(f"curved_cap σ∮t×n ds = {F:.3f} N", fontsize=9)
    fig.savefig(out / "drawings_demo.png", bbox_inches="tight")
    print(f"drawings demo → {out / 'drawings_demo.png'}")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
