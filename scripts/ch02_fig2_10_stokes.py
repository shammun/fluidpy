"""Chapter 2, §2.13, Eqs. (2.34)–(2.35) and Fig. 2.10 (our drawing): an open cap surface with n, n_c and t = n_c × n
at a boundary point; circulation vs curl flux for solid-body rotation (2|b|πR²), a shear flow and an irrotational
vortex (hypothesis failure), and the shrinking-loop limit Γ/A → n·(∇×u) (Example 2.6).

Run: ``.venv/Scripts/python.exe scripts/ch02_fig2_10_stokes.py --no-show``
Figure → outputs/ch02/fig2_10_stokes.png.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))


def stokes_cap_orientation(shape: str = "hemisphere", R: float = 1.0, phi0: float = -2.5, flip: bool = False) -> dict:
    """The three unit vectors of Fig. 2.10 at the rim point P = R(cos φ₀, sin φ₀, 0) of the cap.

    n = surface normal at P originating from the chosen *outside* (convex side of the hemisphere → horizontal outward;
    upward e₃ for the flat disc); n_c = unit normal to C tangent to A, pointing *into* A (up the cap, or towards the
    disc centre); t = n_c × n (``ch02.boundary_tangent``). ``flip=True`` chooses the other side as outside: n and t
    reverse, n_c stays. Returns ``{"P", "n", "n_c", "t", "det"}`` with det[n_c, n, t] = +1 always (right-handed).
    """
    from fluidpy import ch02_cartesian_tensors as ch02

    radial = np.array([np.cos(phi0), np.sin(phi0), 0.0])
    P = R * radial
    if shape == "hemisphere":
        n = radial.copy()  # outward from the sphere at the rim: horizontal
        n_c = np.array([0.0, 0.0, 1.0])  # tangent to A, ⊥ C, up the cap (into A)
    elif shape == "disc":
        n = np.array([0.0, 0.0, 1.0])
        n_c = -radial  # tangent to A, ⊥ C, towards the centre (into A)
    else:
        raise ValueError('shape must be "hemisphere" or "disc"')
    if flip:
        n = -n
    t = ch02.boundary_tangent(n_c, n)  # t = n_c × n  (Fig. 2.10)
    return {"P": P, "n": n, "n_c": n_c, "t": t, "det": float(np.linalg.det(np.stack([n_c, n, t])))}


def stokes_cap_figure(shape: str = "hemisphere", R: float = 1.0, phi0: float = -2.5, flip: bool = False,
                      n_rim_arrows: int = 6, height: int = 560):
    """Rotatable plotly Fig. 2.10: an open surface A (hemispherical cap by default, or the flat unit disc), its boundary
    C with the tangent t, and at one rim point the right-handed triad (n_c, n, t) with t = n_c × n.

    Book: §2.13, Fig. 2.10 and Eq. (2.34) ∬_A (∇×u)·n dA = ∮_C u·t ds — "the inside and outside of A must be chosen";
    n originates from the outside, n_c is perpendicular to C but tangent to A, and n_c × n = t. Drawing only — the
    vectors come from :func:`stokes_cap_orientation` (``ch02.boundary_tangent``).

    Parameters
    ----------
    shape : "hemisphere" (the book's cap z = √(R² − r²)) or "disc" (A = the unit disc in z = 0, n = e₃)
    R : radius of C [drawing units]
    phi0 : polar angle [rad] of the rim point P where the triad is drawn
    flip : choose the other side as outside — n and t reverse, n_c stays (the storyboard's "change")
    n_rim_arrows : small t arrows along C showing the sense of traversal (counterclockwise seen from the outside)
    height : figure height [px]

    Returns
    -------
    go.Figure: surface A (teal), C (purple) with t arrows, n dA at a patch of A (black), and at P: n (black),
    n_c (orange), t (purple); the title states t = n_c × n and det[n_c, n, t] = +1.
    """
    import plotly.graph_objects as go

    from ch02_drawings import plotly_arrow, plotly_layout_3d
    from fluidpy.core.style import COLORS

    o = stokes_cap_orientation(shape, R, phi0, flip)
    P, n, n_c, t = o["P"], o["n"], o["n_c"], o["t"]
    sgn = -1.0 if flip else 1.0
    fig = go.Figure()
    # the surface A
    r, th = np.meshgrid(np.linspace(0.0, R, 24), np.linspace(0.0, 2.0 * np.pi, 72), indexing="ij")
    X, Y = r * np.cos(th), r * np.sin(th)
    Z = np.sqrt(np.maximum(R ** 2 - r ** 2, 0.0)) if shape == "hemisphere" else np.zeros_like(r)
    fig.add_trace(go.Surface(x=X, y=Y, z=Z, colorscale=[[0, COLORS["teal"]], [1, COLORS["teal"]]], showscale=False,
                             opacity=0.35, name="A", hoverinfo="name", showlegend=True))
    # n dA at a patch (top of the cap / centre of the disc) — the patch is a small square of the surface
    top = np.array([0.0, 0.0, R if shape == "hemisphere" else 0.0])
    n_top = sgn * np.array([0.0, 0.0, 1.0])
    fig.add_traces(plotly_arrow(top, 0.45 * R * n_top, COLORS["ink"], width=4, head=0.3, name="n dA (at a patch of A)",
                                showlegend=True, hovertext="n dA: n from the outside of A"))
    # the boundary C, traversed in the sense of t
    phi = np.linspace(0.0, 2.0 * np.pi, 241)
    fig.add_trace(go.Scatter3d(x=R * np.cos(phi), y=R * np.sin(phi), z=0 * phi, mode="lines",
                               line=dict(color=COLORS["accent"], width=6), name="C", hoverinfo="name"))
    for k in range(int(n_rim_arrows)):
        ph = phi0 + 2.0 * np.pi * (k + 0.5) / n_rim_arrows
        tk = sgn * np.array([-np.sin(ph), np.cos(ph), 0.0])  # counterclockwise about e₃ when the outside is up/out
        fig.add_traces(plotly_arrow(R * np.array([np.cos(ph), np.sin(ph), 0.0]), 0.22 * R * tk, COLORS["accent"],
                                    width=3, head=0.6, name="t along C", showlegend=(k == 0)))
    # the triad at P
    fig.add_traces(plotly_arrow(P, 0.55 * R * n, COLORS["ink"], width=6, head=0.28, name="n (from the outside)",
                                showlegend=True, hovertext=f"n = {np.round(n, 3)}"))
    fig.add_traces(plotly_arrow(P, 0.55 * R * n_c, COLORS["orange"], width=6, head=0.28,
                                name="n<sub>c</sub> (⊥ C, tangent to A, into A)", showlegend=True,
                                hovertext=f"n_c = {np.round(n_c, 3)}"))
    fig.add_traces(plotly_arrow(P, 0.55 * R * t, COLORS["accent"], width=8, head=0.28, name="t = n<sub>c</sub> × n",
                                showlegend=True, hovertext=f"t = n_c × n = {np.round(t, 3)}"))
    for vec, lab, col in ((n, "n", COLORS["ink"]), (n_c, "n<sub>c</sub>", COLORS["orange"]), (t, "t", COLORS["accent"])):
        tip = P + 0.66 * R * vec
        fig.add_trace(go.Scatter3d(x=[tip[0]], y=[tip[1]], z=[tip[2]], mode="text", text=[lab],
                                   textfont=dict(color=col, size=15), hoverinfo="skip", showlegend=False))
    fmt = lambda v: "(" + ", ".join(f"{(0.0 if abs(c) < 1e-9 else c):.2g}" for c in v) + ")"  # noqa: E731
    side = ("inside chosen as outside" if flip else ("convex side is the outside" if shape == "hemisphere"
                                                     else "upper side is the outside"))
    title = (f"Fig. 2.10 idea ({shape}, {side}): t = n<sub>c</sub> × n; (n<sub>c</sub>, n, t) right-handed, "
             f"det = {o['det']:+.0f}<br>P: n = {fmt(n)}, n<sub>c</sub> = {fmt(n_c)}, t = {fmt(t)} · "
             f"t counterclockwise seen from the outside")
    lim = 1.45 * R
    return plotly_layout_3d(fig, ((-lim, lim), (-lim, lim), (-0.6 * R, 1.6 * R)), title=title, height=height,
                            camera_eye=(1.3, -1.7, 1.25))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch02"))
    ap.add_argument("--no-show", action="store_true")
    args = ap.parse_args()
    t0 = time.perf_counter()
    import matplotlib

    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from ch02_drawings import arrow3d, setup_3d
    from fluidpy import ch02_cartesian_tensors as ch02
    from fluidpy.core.style import COLORS, use_style
    from tools.convergence import observed_order

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    # --- Fig. 2.10: a cap z = 1 − r² over the unit disc, boundary circle C in z = 0
    fig = plt.figure(figsize=(14, 4.8))
    ax = fig.add_subplot(131, projection="3d")
    setup_3d(ax, lim=1.3)
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_zlim(-0.2, 1.3)
    r, th = np.meshgrid(np.linspace(0, 1, 15), np.linspace(0, 2 * np.pi, 60))
    ax.plot_surface(r * np.cos(th), r * np.sin(th), 1 - r ** 2, alpha=0.25, color=COLORS["teal"], edgecolor="none")
    phi = np.linspace(0, 2 * np.pi, 200)
    ax.plot(np.cos(phi), np.sin(phi), 0 * phi, color=COLORS["accent"], lw=2, label="C")
    P = np.array([1.0, 0.0, 0.0])
    # Fig. 2.10: n_c is ⊥ C, tangent to A and points *into* A (up the cap; ∂z/∂r = −2 at the rim → direction (−1, 0, 2))
    n_c = np.array([-1.0, 0.0, 2.0]) / np.sqrt(5.0)
    n = np.array([2.0, 0.0, 1.0]) / np.sqrt(5.0)  # surface normal at the rim from the chosen outside (convex side, up)
    t = ch02.boundary_tangent(n_c, n)  # t = n_c × n = e₂ here: counterclockwise seen from above
    arrow3d(ax, P, 0.6 * n, color=COLORS["ink"], label="$\\mathbf{n}$", lw=2)
    arrow3d(ax, P, 0.6 * n_c, color=COLORS["orange"], label="$\\mathbf{n}_c$", lw=2)
    arrow3d(ax, P, 0.6 * t, color=COLORS["accent"], label="$\\mathbf{t} = \\mathbf{n}_c\\times\\mathbf{n}$", lw=2)
    top = np.array([0.0, 0.0, 1.0])
    arrow3d(ax, top, 0.5 * np.array([0, 0, 1.0]), color=COLORS["ink"], label="$\\mathbf{n}\\,dA$", lw=1.5)
    ax.set_title("Fig. 2.10 idea: choose the outside → n;\nt = n_c × n runs counterclockwise seen from outside", fontsize=9.5)
    ax.view_init(elev=20, azim=-70)

    # --- both sides of (2.34) for three plane fields on a disc of radius R
    R = 0.8
    loop = ch02.planar_loop([0.0, 0.0], radius=R)
    disc = ch02.planar_disc([0.0, 0.0], radius=R)
    fields = {"solid-body b × x (b = 1)": ch02.solid_body_rotation_field(1.0, dim=2),
              "simple shear Γx₂ (Γ = 1)": ch02.shear_field(1.0),
              "irrotational vortex K/r": ch02.irrotational_vortex_field(1.0),
              "potential ∇(x² − y²)": ch02.potential_field()}
    results = {name: ch02.stokes_theorem_check(f, loop, disc, curl_fn=f.curl_fn) for name, f in fields.items()}
    ax = fig.add_subplot(132)
    names = list(results)
    w = 0.38
    ax.bar(np.arange(len(names)) - w / 2, [results[k].lhs for k in names], w, color=COLORS["teal"], label="∬ (∇×u)·n dA")
    ax.bar(np.arange(len(names)) + w / 2, [results[k].rhs for k in names], w, color=COLORS["accent"], label="∮ u·t ds")
    ax.set_xticks(np.arange(len(names)), [k.split(" (")[0] for k in names], rotation=12, fontsize=8)
    ax.axhline(0, color=COLORS["grid"])
    ax.set_title(f"Stokes (2.34) on a disc R = {R}:\nboth sides agree unless the core is inside", fontsize=9.5)
    ax.legend(fontsize=8.5)
    ax.text(2, 0.5 * results[names[2]].rhs, "⚠ core inside:\nStokes does not apply", ha="center", fontsize=8, color=COLORS["rose"])

    # --- Example 2.6 / (2.35): shrinking square loop around x0 in the shear flow
    ax = fig.add_subplot(133)
    u3 = ch02.smooth_test_field(3)
    x0 = np.array([0.3, -0.2, 0.5])
    hs = [0.4, 0.2, 0.1, 0.05]
    vals = [ch02.integral_curl_component(u3, x0, [0, 0, 1], h) for h in hs]
    exact = float(u3.curl_fn(*x0)[2])
    err = [abs(v - exact) for v in vals]
    order = observed_order(hs, err)
    ax.loglog(hs, err, "o-", color=COLORS["accent"], label=f"|Γ/A − (∇×u)₃|, slope {order:.2f}")
    ax.loglog(hs, err[0] * (np.array(hs) / hs[0]) ** 2, "--", color=COLORS["muted"], label="∝ h²")
    ax.set_xlabel("loop side h [m]")
    ax.set_ylabel("error [1/s]")
    ax.set_title("Example 2.6 / (2.35):\ncirculation per area → normal curl", fontsize=9.5)
    ax.legend(fontsize=8.5)
    fig.savefig(out / "fig2_10_stokes.png", bbox_inches="tight")

    # the notebook's rotatable Fig. 2.10 (default hemisphere; flipped; flat disc)
    for label, kw in (("hemisphere", {}), ("hemisphere_flip", {"flip": True}), ("disc", {"shape": "disc"})):
        fig_p = stokes_cap_figure(**kw)
        fig_p.write_html(out / f"fig2_10_stokes_cap_{label}.html", include_plotlyjs="cdn")
        o = stokes_cap_orientation(**kw)
        print(f"plotly cap [{label:15s}]: n = {np.round(o['n'], 3)}, n_c = {np.round(o['n_c'], 3)}, t = n_c × n = {np.round(o['t'], 3)}, "
              f"det[n_c, n, t] = {o['det']:+.3f}, {len(fig_p.data)} traces")

    print(f"orientation at P: n_c = {np.round(n_c, 3)}, n = {np.round(n, 3)}, t = {np.round(t, 3)}, det[n_c, n, t] = {np.linalg.det(np.stack([n_c, n, t])):+.3f}")
    for k, res in results.items():
        print(f"{k:28s}: curl flux {res.lhs:+.6f}, circulation {res.rhs:+.6f}, hypothesis ok: {res.hypothesis_ok}")
    print(f"solid body: 2|b|πR² = {2 * np.pi * R ** 2:.6f}; vortex: 2πK = {2 * np.pi:.6f} regardless of R")
    print(f"(2.35) at x0 = {x0}: Γ/A = {[round(v, 6) for v in vals]} → (∇×u)₃ = {exact:.6f}; observed order {order:.3f}")
    print(f"saved fig2_10_stokes.png in {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
