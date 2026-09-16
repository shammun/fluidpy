"""Chapter 2, §2.12, Eqs. (2.30)–(2.32) and Fig. 2.9 (our drawing): a box whose faces are coloured by the outflux
n·Q, the two sides of the divergence theorem converging as the quadrature refines, the sphere benchmark 8π/3, and the
small-box limit (1/V)∮ n·Q dA → ∇·Q (Example 2.5).

Run: ``.venv/Scripts/python.exe scripts/ch02_fig2_9_gauss.py --no-show``
Figure → outputs/ch02/fig2_9_gauss.png.
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


def _components(Q_fn, X, Y, Z) -> np.ndarray:
    """Evaluate a vector field callable on broadcast grids → array (3, ...) (tuples/lists and constants accepted)."""
    out = Q_fn(X, Y, Z)
    if isinstance(out, (list, tuple)):
        return np.stack([np.broadcast_to(np.asarray(c, dtype=float), X.shape) for c in out])
    out = np.asarray(out, dtype=float)
    return np.stack([np.broadcast_to(c, X.shape) for c in out])


def gauss_box_figure(Q_fn, bounds=((0.0, 1.0), (0.0, 1.0), (0.0, 1.0)), n: int = 16, n_face: int = 12, height: int = 560):
    """Rotatable plotly Fig. 2.9: the box whose six faces are coloured by the local outflux n·Q, with the two sides of
    the divergence theorem in the title and each face's total flux written on it.

    Book: §2.12, Eq. (2.30) in its vector form ∭ ∇·Q dV = ∯ n·Q dA ("the volume integral of the divergence of Q is
    equal to the surface integral of the outflux of Q"). Drawing only — the totals come from
    ``ch02.divergence_theorem_box`` and ``ch02.flux_through_box_faces``.

    Parameters
    ----------
    Q_fn : callable ``Q(X, Y, Z) → (3, ...)`` (or a tuple of three arrays) [Q unit]
    bounds : ((x0, x1), (y0, y1), (z0, z1)) [m]
    n : quadrature nodes per direction for the two totals (midpoint rule)
    n_face : vertices per face edge for the colour mesh
    height : figure height [px]

    Returns
    -------
    go.Figure: six ``go.Mesh3d`` faces with intensity n·Q on a shared blue (inflow) – white – orange (outflow) scale,
    six outward-normal arrows, one flux label per face.
    """
    import plotly.graph_objects as go

    from ch02_drawings import plotly_arrow, plotly_layout_3d
    from fluidpy import ch02_cartesian_tensors as ch02
    from fluidpy.core.style import COLORS

    b = np.asarray(bounds, dtype=float)
    lhs, rhs = ch02.divergence_theorem_box(Q_fn, b, n)  # both sides of (2.30)
    face_flux = ch02.flux_through_box_faces(Q_fn, b, n)
    centre = b.mean(axis=1)
    size = b[:, 1] - b[:, 0]

    # sample n·Q on every face first, so all six share one colour scale
    s = np.linspace(0.0, 1.0, n_face)
    A, B = np.meshgrid(s, s, indexing="ij")
    faces = []
    for d in range(3):
        others = [k for k in range(3) if k != d]
        for sign, key in ((+1, "+"), (-1, "-")):
            P = np.empty((3,) + A.shape)
            P[d] = b[d, 1] if sign > 0 else b[d, 0]
            P[others[0]] = b[others[0], 0] + size[others[0]] * A
            P[others[1]] = b[others[1], 0] + size[others[1]] * B
            Q = _components(Q_fn, P[0], P[1], P[2])
            faces.append((d, sign, f"{key}{'xyz'[d]}", P, sign * Q[d]))  # n·Q with n = ±e_d
    vmax = max(1e-300, max(float(np.max(np.abs(fl))) for *_, fl in faces))
    idx = np.arange(n_face * n_face).reshape(n_face, n_face)
    tri_i = np.concatenate([idx[:-1, :-1].ravel(), idx[1:, :-1].ravel()])
    tri_j = np.concatenate([idx[1:, :-1].ravel(), idx[1:, 1:].ravel()])
    tri_k = np.concatenate([idx[:-1, 1:].ravel(), idx[:-1, 1:].ravel()])
    scale = [[0.0, COLORS["blue"]], [0.5, "#ffffff"], [1.0, COLORS["orange"]]]

    fig = go.Figure()
    # the twelve box edges (so the white faces of a zero-flux side still read as a box)
    corners = np.array([[b[0, a], b[1, c], b[2, e]] for a in (0, 1) for c in (0, 1) for e in (0, 1)])
    ex, ey, ez = [], [], []
    for p in range(8):
        for q in range(p + 1, 8):
            if np.sum(np.abs(corners[p] - corners[q]) > 1e-12) == 1:
                ex += [corners[p, 0], corners[q, 0], None]
                ey += [corners[p, 1], corners[q, 1], None]
                ez += [corners[p, 2], corners[q, 2], None]
    fig.add_trace(go.Scatter3d(x=ex, y=ey, z=ez, mode="lines", line=dict(color=COLORS["ink"], width=2),
                               hoverinfo="skip", showlegend=False, name="box"))
    for m, (d, sign, key, P, fl) in enumerate(faces):
        fig.add_trace(go.Mesh3d(x=P[0].ravel(), y=P[1].ravel(), z=P[2].ravel(), i=tri_i, j=tri_j, k=tri_k,
                                intensity=fl.ravel(), colorscale=scale, cmin=-vmax, cmax=vmax, showscale=(m == 0),
                                colorbar=dict(title="n·Q", len=0.6, x=0.98), opacity=0.92, flatshading=False,
                                hovertemplate=f"face {key}: n·Q = %{{intensity:.3f}}<extra></extra>",
                                name=f"{key}: {face_flux[key]:+.3f}", showlegend=True))  # legend lists all six totals
        nvec = np.zeros(3)
        nvec[d] = sign
        c = centre.copy()
        c[d] = b[d, 1] if sign > 0 else b[d, 0]
        fig.add_traces(plotly_arrow(c, 0.28 * float(np.min(size)) * nvec, COLORS["ink"], width=3, head=0.35,
                                    name=f"n = {'+' if sign > 0 else '-'}e{d + 1}"))
        tip = c + 0.5 * float(np.min(size)) * nvec
        fig.add_trace(go.Scatter3d(x=[tip[0]], y=[tip[1]], z=[tip[2]], mode="text",
                                   text=[f"{key}: {face_flux[key]:+.3f}"], textfont=dict(size=12, color=COLORS["ink"]),
                                   hoverinfo="skip", showlegend=False))
    pad = 0.55 * float(np.min(size))
    ranges = tuple((b[k, 0] - pad, b[k, 1] + pad) for k in range(3))
    title = (f"Fig. 2.9 idea: faces coloured by the outflux n·Q (orange out, blue in)"
             f"<br>∭ ∇·Q dV = {lhs:.4f}   ∯ n·Q dA = Σ faces = {rhs:.4f}   (midpoint rule, n = {n})")
    return plotly_layout_3d(fig, ranges, title=title, height=height)


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
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    from ch02_drawings import arrow3d, setup_3d
    from fluidpy import ch02_cartesian_tensors as ch02
    from fluidpy.core.style import COLORS, use_style
    from tools.convergence import observed_order

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    Q = ch02.smooth_test_field(3)  # smooth non-polynomial: the quadrature error is visible
    bounds = ((0.0, 1.0), (0.0, 1.0), (0.0, 1.0))
    ns = [4, 8, 16, 32]
    both = [ch02.divergence_theorem_box(Q, bounds, n, div_fn=Q.div_fn) for n in ns]
    fine = ch02.divergence_theorem_box(Q, bounds, 33, div_fn=Q.div_fn, rule="simpson")
    ref = 0.5 * (fine[0] + fine[1])
    gap = [abs(l - r) for l, r in both]
    faces = ch02.flux_through_box_faces(Q, bounds, 32)
    sphere = ch02.divergence_theorem_sphere(lambda X, Y, Z: (2 * X, Y ** 2, Z ** 2), 1.0, 24)
    # Example 2.5: shrinking box around x0
    x0 = np.array([0.3, -0.2, 0.5])
    hs = [0.4, 0.2, 0.1, 0.05]
    limit = [ch02.integral_divergence(Q, x0, h) for h in hs]
    exact = float(Q.div_fn(*x0))
    order = observed_order(hs, [abs(v - exact) for v in limit])

    fig = plt.figure(figsize=(14, 4.8))
    ax = fig.add_subplot(131, projection="3d")
    setup_3d(ax, lim=1.3)
    ax.set_xlim(-0.1, 1.2)
    ax.set_ylim(-0.1, 1.2)
    ax.set_zlim(-0.1, 1.2)
    vmax = max(abs(v) for v in faces.values())
    cmap = plt.get_cmap("coolwarm")
    for key, val in faces.items():
        d = "xyz".index(key[1])
        v0 = 1.0 if key[0] == "+" else 0.0
        corners = []
        for s_, t_ in ((0, 0), (1, 0), (1, 1), (0, 1)):
            p = [0.0, 0.0, 0.0]
            p[d] = v0
            others = [k for k in range(3) if k != d]
            p[others[0]], p[others[1]] = s_, t_
            corners.append(tuple(p))
        ax.add_collection3d(Poly3DCollection([corners], facecolor=cmap(0.5 + 0.5 * val / vmax), alpha=0.55, edgecolor=COLORS["ink"]))
        n = np.zeros(3)
        n[d] = 1.0 if key[0] == "+" else -1.0
        c = np.full(3, 0.5)
        c[d] = v0
        arrow3d(ax, c, 0.25 * n, color=COLORS["ink"], lw=1.2)
    ax.set_title("Fig. 2.9 idea: faces coloured by the outflux ∬ n·Q dA\n(red out, blue in); their sum is the right side of (2.30)", fontsize=9.5)
    ax.view_init(elev=22, azim=-55)
    ax = fig.add_subplot(132)
    w = 0.38
    ax.bar(np.arange(len(ns)) - w / 2, [b[0] for b in both], w, color=COLORS["teal"], label="∭ ∇·Q dV")
    ax.bar(np.arange(len(ns)) + w / 2, [b[1] for b in both], w, color=COLORS["orange"], label="∯ n·Q dA")
    ax.axhline(ref, color=COLORS["muted"], ls="--", label="Simpson n = 33")
    ax.set_xticks(np.arange(len(ns)), [f"n = {n}" for n in ns])
    ax.set_ylim(min(min(b) for b in both) - 0.02, max(max(b) for b in both) + 0.02)
    ax.set_title("both sides of the divergence theorem\nvs quadrature nodes per direction", fontsize=9.5)
    ax.legend(fontsize=8.5)
    ax = fig.add_subplot(133)
    ax.loglog(hs, [abs(v - exact) for v in limit], "o-", color=COLORS["accent"], label=f"|(1/V)∮n·Q dA − ∇·Q(x₀)|, slope {order:.2f}")
    ax.loglog(hs, [abs(v - exact) for v in limit][0] * (np.array(hs) / hs[0]) ** 2, "--", color=COLORS["muted"], label="∝ h² (D22)")
    ax.set_xlabel("box side h [m]")
    ax.set_ylabel("error [1/s]")
    ax.set_title("Example 2.5: the shrinking box\ngives the divergence at x₀", fontsize=9.5)
    ax.legend(fontsize=8.5)
    fig.savefig(out / "fig2_9_gauss.png", bbox_inches="tight")

    # the notebook's rotatable box (storyboard: F = (2x, y², z²); "change": b × x is solenoidal → net zero)
    F = lambda X, Y, Z: np.stack([2 * X, Y ** 2, Z ** 2])  # noqa: E731
    fig_p = gauss_box_figure(F)
    fig_p.write_html(out / "fig2_9_gauss_box.html", include_plotlyjs="cdn")
    fig_rot = gauss_box_figure(ch02.solid_body_rotation_field((0.0, 0.0, 1.0)), n=8)
    print(f"plotly gauss box F = (2x, y², z²): {len(fig_p.data)} traces; title: {fig_p.layout.title.text.split('<br>')[1]}")
    print(f"plotly gauss box b × x: {fig_rot.layout.title.text.split('<br>')[1]}")
    # 2-D tiling (D25 step 8): the interior faces cancel, the tile sum equals the outer flux
    Q2d = ch02.radial_field(1.0, dim=2)
    tiled = ch02.divergence_theorem_tiled(Q2d, ((0.0, 1.0), (0.0, 1.0)), tiles=4, n=16)
    print(f"tiled 4×4 (Q = (x, y) on the unit square): sum of tile boundaries {tiled.sum_tiles:.15f}, outer {tiled.outer:.15f}, "
          f"interior {tiled.interior:+.1e}; rect2d surface side {ch02.divergence_theorem_rect2d(Q2d, ((0, 1), (0, 1)), 64)[1]:.15f}")

    for n, (l, r) in zip(ns, both):
        print(f"n = {n:3d}: ∭ ∇·Q dV = {l:.8f}, ∯ n·Q dA = {r:.8f}, gap {abs(l - r):.2e}")
    print(f"gap order (midpoint) = {observed_order([1 / n for n in ns], gap):.2f}; Simpson n = 33: {fine[0]:.10f} vs {fine[1]:.10f}")
    print(f"face fluxes (n = 32): {{{', '.join(f'{k}: {v:+.4f}' for k, v in faces.items())}}}")
    print(f"unit sphere, F = (2x, y², z²): volume {sphere[0]:.12f}, surface {sphere[1]:.12f}, 8π/3 = {8 * np.pi / 3:.12f}")
    print(f"Example 2.5 at x0 = {x0}: (1/V)∮ = {[round(v, 6) for v in limit]} → exact {exact:.6f}, observed order {order:.3f}")
    print(f"Q = (x, y, z) on the unit cube: {ch02.divergence_theorem_box(lambda X, Y, Z: (X, Y, Z), bounds, 4)}")
    print(f"saved fig2_9_gauss.png in {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
