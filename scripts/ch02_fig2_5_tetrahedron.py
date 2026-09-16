"""Chapter 2, §2.6, Eq. (2.15) and Fig. 2.5 (our drawing): the tetrahedron whose slanted face has normal n; the
coordinate faces have areas dA_i = n_i dA and the traction on the slanted face is f_i = τ_ji n_j.

Run: ``.venv/Scripts/python.exe scripts/ch02_fig2_5_tetrahedron.py --no-show``
Figure → outputs/ch02/fig2_5_tetrahedron.png.
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


def tetrahedron_geometry(n, size: float = 1.0):
    """Vertices of the slanted face with unit normal ``n`` and of its three coordinate-plane shadows.

    If every |n_i| > 0.1 the slanted face has its vertices on the axes (the book's Fig. 2.5 tetrahedron; a negative n_i
    puts that vertex on the negative axis); otherwise a triangle in the plane ⊥ n is used and the coordinate faces are
    its projections (a zero n_k gives a degenerate shadow of zero area, as dA_k = n_k dA demands).

    Returns ``(slanted (3, 3), shadows [(3, 3)] × 3, dA, dA_i)`` with ``shadows[k]`` the face in the plane x_k = 0.
    """
    from fluidpy import ch02_cartesian_tensors as ch02

    n = np.asarray(n, dtype=float)
    n = n / np.linalg.norm(n)
    if np.min(np.abs(n)) > 0.1:
        d = size * float(np.min(np.abs(n)))  # plane n·x = d; the largest vertex coordinate is then ≤ size
        slanted = np.diag(d / n)  # vertex k = (d / n_k) e_k
    else:
        e1 = np.cross(n, [0.0, 0.0, 1.0]) if abs(n[2]) < 0.9 else np.cross(n, [1.0, 0.0, 0.0])
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(n, e1)
        c = 0.75 * size * n
        phi = np.deg2rad([90.0, 210.0, 330.0])
        slanted = c + 0.55 * size * (np.cos(phi)[:, None] * e1 + np.sin(phi)[:, None] * e2)
    shadows = []
    for k in range(3):
        sh = slanted.copy()
        sh[:, k] = 0.0  # projection onto the plane x_k = 0
        shadows.append(sh)
    dA = 0.5 * float(np.linalg.norm(np.cross(slanted[1] - slanted[0], slanted[2] - slanted[0])))
    dA_i = ch02.tetrahedron_face_areas(n, dA)  # dA_i = n_i dA (D05)
    return slanted, shadows, dA, dA_i


def tetrahedron_figure(tau, n, size: float = 1.0, height: int = 560):
    """Rotatable plotly Fig. 2.5: the slanted face with normal n carries f = n·τ; the three coordinate faces carry the
    reversed rows of τ; the four arrows balance (Cauchy's formula as a force balance).

    Book: §2.6, Eq. (2.15) f_i = τ_ji n_j, derived from the force balance on the tetrahedron whose coordinate faces have
    areas dA_i = n_i dA (D05). Drawing only — the numbers come from ``ch02.traction``, ``ch02.normal_shear_stress`` and
    ``ch02.tetrahedron_face_areas``.

    Parameters
    ----------
    tau : (3, 3) stress tensor [Pa]
    n : outward unit normal of the slanted face (normalised here); n_i ≥ 0 gives the book's picture, zeros allowed
    size : drawing size [drawing units]
    height : figure height [px]

    Returns
    -------
    go.Figure: coordinate faces (teal, orange, blue shadows) each with its traction −τ_k· (grey arrows, lengths ∝ |τ|);
    the slanted face (purple) with f (bold rose), its normal part σ_n n (blue), and n itself (thin black). The title
    reports f, σ_n, τ_s and the balance residual |f dA − Σ_k τ_k· dA_k|.
    """
    import plotly.graph_objects as go

    from ch02_drawings import plotly_arrow, plotly_layout_3d
    from fluidpy import ch02_cartesian_tensors as ch02
    from fluidpy.core.style import COLORS

    tau = np.asarray(tau, dtype=float)
    n = np.asarray(n, dtype=float)
    n = n / np.linalg.norm(n)
    slanted, shadows, dA, dA_i = tetrahedron_geometry(n, size)
    f = ch02.traction(tau, n)  # Eq. (2.15): f_i = τ_ji n_j
    sigma_n, tau_s, _ = ch02.normal_shear_stress(tau, n)
    residual = f * dA - sum(np.sign(n[k]) * tau[k] * abs(dA_i[k]) for k in range(3))  # D05 force balance
    scale = 0.55 * size / max(1e-300, float(np.max(np.abs(tau))))

    fig = go.Figure()
    cols = [COLORS["teal"], COLORS["orange"], COLORS["blue"]]
    for k in range(3):
        v = shadows[k]
        fig.add_trace(go.Mesh3d(x=v[:, 0], y=v[:, 1], z=v[:, 2], i=[0], j=[1], k=[2], color=cols[k], opacity=0.35,
                                name=f"face x{k + 1} = 0: dA{k + 1} = n{k + 1} dA = {dA_i[k]:.3f}", showlegend=True,
                                hoverinfo="name"))
        if abs(n[k]) > 1e-12:
            # outward normal of this face is −sign(n_k) e_k → traction −sign(n_k) (row k of τ)
            vec = -np.sign(n[k]) * tau[k] * scale
            fig.add_traces(plotly_arrow(v.mean(axis=0), vec, COLORS["muted"], width=4,
                                        name=f"−(row {k + 1} of τ) on face x{k + 1} = 0",
                                        hovertext=f"traction on the x{k + 1} = 0 face: {np.round(-np.sign(n[k]) * tau[k], 3)} Pa"))
    fig.add_trace(go.Mesh3d(x=slanted[:, 0], y=slanted[:, 1], z=slanted[:, 2], i=[0], j=[1], k=[2],
                            color=COLORS["accent"], opacity=0.45, name=f"slanted face, dA = {dA:.3f}", showlegend=True,
                            hoverinfo="name"))
    cen = slanted.mean(axis=0)
    fig.add_traces(plotly_arrow(cen, 0.45 * size * n, COLORS["ink"], width=3, head=0.25, name="n",
                                hovertext=f"n = {np.round(n, 3)}", showlegend=True))
    fig.add_traces(plotly_arrow(cen, scale * f, COLORS["rose"], width=8, name=f"f = n·τ = {np.round(f, 3)} Pa",
                                hovertext=f"f_i = τ_ji n_j = {np.round(f, 3)} Pa", showlegend=True))
    fig.add_traces(plotly_arrow(cen, scale * sigma_n * n, COLORS["blue"], width=5,
                                name=f"σ_n n (σ_n = {sigma_n:+.3f} Pa)", hovertext=f"normal part σ_n = {sigma_n:+.4f} Pa",
                                showlegend=True))
    lim = max(1.05 * float(np.max(np.abs(slanted))), 0.5 * size)
    ranges = tuple((min(-0.15 * size, float(slanted[:, k].min()) - 0.1 * size), max(lim, float(slanted[:, k].max()) + 0.1 * size))
                   for k in range(3))
    title = (f"Fig. 2.5 idea: f = n·τ = ({f[0]:+.3g}, {f[1]:+.3g}, {f[2]:+.3g}) Pa balances the three shadows' tractions"
             f"<br>σ<sub>n</sub> = {sigma_n:+.3f} Pa, τ<sub>s</sub> = {tau_s:.3f} Pa, balance residual {np.max(np.abs(residual)):.1e}")
    return plotly_layout_3d(fig, ranges, title=title, height=height, camera_eye=(1.7, 1.5, 1.0))


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

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    # tetrahedron with vertices on the axes: the slanted face has normal ∝ (1/a, 1/b, 1/c)
    a, b, c = 1.2, 1.0, 0.8
    A, B, Cv = np.array([a, 0, 0]), np.array([0, b, 0]), np.array([0, 0, c])
    n = np.cross(B - A, Cv - A)
    dA = 0.5 * np.linalg.norm(n)  # area of the slanted face
    n /= np.linalg.norm(n)
    dA_i = ch02.tetrahedron_face_areas(n, dA)  # dA_i = n_i dA
    coord_face_areas = np.array([0.5 * b * c, 0.5 * a * c, 0.5 * a * b])  # geometry, independently
    assert np.allclose(dA_i, coord_face_areas)  # the vector-area identity behind D05

    tau = np.array([[-1.0, 0.6, 0.2], [0.6, -1.0, 0.0], [0.2, 0.0, -1.0]])  # Pa
    f = ch02.traction(tau, n)  # Eq. (2.15)
    sigma_n, tau_s, s_dir = ch02.normal_shear_stress(tau, n)
    # force balance on the tetrahedron (D05): f dA − Σ_j (row j of τ) dA_j = 0
    residual = f * dA - sum(tau[j] * dA_i[j] for j in range(3))

    fig = plt.figure(figsize=(7.5, 6))
    ax = fig.add_subplot(111, projection="3d")
    setup_3d(ax, lim=1.5)
    faces = [(A, B, Cv), (np.zeros(3), B, Cv), (np.zeros(3), A, Cv), (np.zeros(3), A, B)]
    cols = [COLORS["accent"], COLORS["teal"], COLORS["orange"], COLORS["blue"]]
    for verts, col in zip(faces, cols):
        ax.add_collection3d(Poly3DCollection([list(map(tuple, verts))], alpha=0.25, facecolor=col, edgecolor=COLORS["ink"]))
    cen = (A + B + Cv) / 3
    arrow3d(ax, cen, 0.5 * n, color=COLORS["ink"], label="$\\mathbf{n}$", lw=2)
    arrow3d(ax, cen, 0.6 * f / np.linalg.norm(f), color=COLORS["rose"], label="$\\mathbf{f}(\\mathbf{n}) = \\mathbf{n}\\cdot\\tau$", lw=2.4)
    for k, (lab, col) in enumerate(zip(("$dA_1$", "$dA_2$", "$dA_3$"), cols[1:])):
        centre = np.mean(faces[k + 1], axis=0)
        ax.text(*centre, f"{lab} = n{k + 1} dA = {dA_i[k]:.3f}", color=col, fontsize=9)
    ax.set_title("Fig. 2.5 idea: coordinate faces $dA_i = n_i\\,dA$ carry −(row i of τ); the slanted face carries f\n"
                 f"σ_n = {sigma_n:.3f} Pa, τ_s = {tau_s:.3f} Pa; balance residual {np.max(np.abs(residual)):.1e}")
    ax.view_init(elev=22, azim=35)
    fig.savefig(out / "fig2_5_tetrahedron.png", bbox_inches="tight")

    # the notebook's rotatable version, with the storyboard's numbers (τ = PRESSURE_SHEAR(3, 1), n = (0.6, 0, 0.8))
    tau_nb = np.array([[-3.0, 1.0, 0.0], [1.0, -3.0, 0.0], [0.0, 0.0, -3.0]])
    for label, n_nb in (("0.6_0_0.8", [0.6, 0.0, 0.8]), ("e3", [0.0, 0.0, 1.0]), ("general", [0.5, 0.6, 0.62])):
        fig_p = tetrahedron_figure(tau_nb, n_nb)
        fig_p.write_html(out / f"fig2_5_tetrahedron_{label}.html", include_plotlyjs="cdn")
        nn = np.asarray(n_nb, float) / np.linalg.norm(n_nb)
        _, _, dA_nb, dAi_nb = tetrahedron_geometry(nn)
        print(f"plotly tetrahedron n = {np.round(nn, 3)}: f = {np.round(ch02.traction(tau_nb, nn), 4)} Pa, "
              f"dA_i/dA = {np.round(dAi_nb / dA_nb, 4)} (= n_i), {len(fig_p.data)} traces")

    print(f"n = {np.round(n, 4)}, dA = {dA:.4f}; dA_i = n_i dA = {np.round(dA_i, 4)} (geometry: {np.round(coord_face_areas, 4)})")
    print(f"f = τ_ji n_j = {np.round(f, 4)} Pa; n·τ − τ·n = {np.max(np.abs(f - tau @ n)):.1e} (τ symmetric here)")
    print(f"normal stress {sigma_n:.4f} Pa, shear {tau_s:.4f} Pa along {np.round(s_dir, 3)}; D05 balance residual {np.max(np.abs(residual)):.1e}")
    print(f"saved fig2_5_tetrahedron.png in {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
