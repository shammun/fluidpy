"""Chapter 2, §2.2, Eqs. (2.3)–(2.7) and Fig. 2.2 (our drawing): two frames sharing an origin, one vector x, its two
component sets, the direction-cosine matrix C_ij = e_i·e'_j and the orthogonality check CᵀC = δ.

Run: ``.venv/Scripts/python.exe scripts/ch02_fig2_2_rotated_axes.py --no-show``
Figure → outputs/ch02/fig2_2_rotated_axes.png.
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


def rotated_axes_figure(C, x=(1.0, 2.0, 3.0), height: int = 560):
    """Rotatable plotly Fig. 2.2: two frames sharing an origin, one arrow x, its components in both frames.

    Book: §2.2, Eqs. (2.3)–(2.5) — x = x_i e_i = x'_j e'_j with x'_j = C_ij x_i, C_ij = e_i·e'_j (Fig. 2.2). Drawing only
    — x' comes from ``ch02.transform_vector(x, C)`` (passive convention: the columns of C are the new axes e'_j written
    in the old frame, x' = Cᵀx).

    Parameters
    ----------
    C : (3, 3) direction-cosine matrix (orthogonal; e.g. ``ch02.rotation_matrix_3d(axis, angle)``)
    x : the arrow's components in the old frame
    height : figure height [px]

    Returns
    -------
    go.Figure: teal old frame e_i, orange new frame e'_j (columns of C), the arrow x (purple), and the dashed
    "component staircase" x₁e₁ → x₁e₁ + x₂e₂ → x in each frame (teal / orange); the title lists x_i, x'_j and shows
    that |x| = |x'| and max|CᵀC − I| ≈ 0.
    """
    import plotly.graph_objects as go

    from ch02_drawings import plotly_arrow, plotly_layout_3d
    from fluidpy import ch02_cartesian_tensors as ch02
    from fluidpy.core.style import COLORS

    C = np.asarray(C, dtype=float)
    x = np.asarray(x, dtype=float)
    xp = ch02.transform_vector(x, C)  # Eq. (2.5): x'_j = x_i C_ij
    E_old = np.eye(3)  # rows e_i
    E_new = C.T  # rows e'_j (columns of C)
    L = 1.15 * float(np.max(np.abs(np.concatenate([x, xp]))))
    fig = go.Figure()
    for k in range(3):
        fig.add_traces(plotly_arrow((0, 0, 0), L * E_old[k], COLORS["teal"], width=4, head=0.12,
                                    name="old frame e<sub>i</sub>", showlegend=(k == 0), hovertext=f"e{k + 1}"))
        fig.add_traces(plotly_arrow((0, 0, 0), L * E_new[k], COLORS["orange"], width=4, head=0.12,
                                    name="new frame e′<sub>j</sub> (columns of C)", showlegend=(k == 0),
                                    hovertext=f"e'{k + 1} = {np.round(E_new[k], 3)}"))
    for k, (E, lab, col) in enumerate(((E_old, "", COLORS["teal"]), (E_new, "′", COLORS["orange"]))):
        tip = 1.06 * L * E
        fig.add_trace(go.Scatter3d(x=tip[:, 0], y=tip[:, 1], z=tip[:, 2], mode="text",
                                   text=[f"{i + 1}{lab}" for i in range(3)], textfont=dict(color=col, size=14),
                                   hoverinfo="skip", showlegend=False))
    # component staircases: origin → x₁e₁ → x₁e₁ + x₂e₂ → x, in each frame
    for comps, E, col, lab in ((x, E_old, COLORS["teal"], "x<sub>i</sub> staircase"),
                               (xp, E_new, COLORS["orange"], "x′<sub>j</sub> staircase")):
        pts = np.vstack([np.zeros(3), np.cumsum(comps[:, None] * E, axis=0)])  # partial sums Σ_{k≤m} x_k e_k
        assert np.allclose(pts[-1], x)  # both staircases end at the same arrow tip (Eq. (2.3))
        fig.add_trace(go.Scatter3d(x=pts[:, 0], y=pts[:, 1], z=pts[:, 2], mode="lines",
                                   line=dict(color=col, width=3, dash="dash"), name=lab,
                                   hovertext=[f"{lab}: {np.round(p, 3)}" for p in pts], hoverinfo="text"))
    fig.add_traces(plotly_arrow((0, 0, 0), x, COLORS["accent"], width=8, head=0.15, name="x", showlegend=True,
                                hovertext=f"x = {np.round(x, 3)} = x' = {np.round(xp, 3)} (same arrow)"))
    fmt = lambda v: "(" + ", ".join(f"{c:.3g}" for c in v) + ")"  # noqa: E731
    title = (f"Fig. 2.2 idea: one arrow, two frames — x<sub>i</sub> = {fmt(x)} (teal), "
             f"x′<sub>j</sub> = C<sub>ij</sub> x<sub>i</sub> = {fmt(xp)} (orange)"
             f"<br>|x| = {np.linalg.norm(x):.4f} = |x′| = {np.linalg.norm(xp):.4f};  max|CᵀC − I| = "
             f"{ch02.orthogonality_residual(C):.0e}")
    lim = 1.2 * L
    return plotly_layout_3d(fig, ((-0.35 * lim, lim),) * 3, title=title, height=height, camera_eye=(1.5, -1.8, 0.9))


def frame_rotation_frames(theta_list, x=(1.0, 2.0), degrees: bool = False, figsize=(9.5, 4.2)):
    """Build the C02 frame-rotation animation: a fixed arrow, teal axes fixed, orange axes turning by θ.

    Book: §2.2, Eqs. (2.5) x'_j = x_i C_ij and (2.6) C_ij C_ik = δ_jk (Fig. 2.2 idea in the x₁x₂ plane). Drawing only —
    the numbers come from ``ch02.rotation_matrix_2d`` (passive C, columns = new axes e'_j) and ``ch02.transform_vector``.

    Parameters
    ----------
    theta_list : sequence of frame angles [rad] (``degrees=True`` → degrees); frame ``i`` shows ``theta_list[i]``
    x : the fixed arrow's components in the old frame (default the storyboard's (1, 2))
    degrees : interpret ``theta_list`` in degrees
    figsize : matplotlib figure size [in]

    Returns
    -------
    (fig, update) : ``update(i)`` redraws frame ``i``; feed it to ``fluidpy.core.anim.animate``::

        fig, update = frame_rotation_frames(np.deg2rad(np.linspace(0, 90, 10)))
        show_animation(animate(update, frames=10, fig=fig, interval=500), player="frames")

    Left panel: the arrow (purple, never moves), the fixed axes (teal), the turning axes (orange) with the dashed
    projections x'_j e'_j onto them. Right panel: bars x_i (teal, constant) and x'_j (orange, moving), the live 2 × 2 C
    and max|CᵀC − I|. At θ = 90° the roles swap: x'_1 = x_2, x'_2 = −x_1.
    """
    import matplotlib.pyplot as plt

    from ch02_drawings import arrow2d, draw_frame_2d
    from fluidpy import ch02_cartesian_tensors as ch02
    from fluidpy.core.style import COLORS

    thetas = np.asarray(theta_list, dtype=float)
    if degrees:
        thetas = np.deg2rad(thetas)
    x = np.asarray(x, dtype=float)
    lim = 1.25 * float(np.linalg.norm(x))
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=figsize, gridspec_kw={"width_ratios": [1.15, 1.0]})

    def update(i: int):
        th = float(thetas[int(i)])
        C = ch02.rotation_matrix_2d(th)  # passive C: C_ij = e_i·e'_j, columns = new axes
        xp = ch02.transform_vector(x, C)  # Eq. (2.5): x'_j = x_i C_ij
        ax.clear()
        ax2.clear()
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_xlabel("$x_1$")
        ax.set_ylabel("$x_2$")
        ax.axhline(0, color=COLORS["grid"], lw=0.8)
        ax.axvline(0, color=COLORS["grid"], lw=0.8)
        # teal labels sit below/left of their axis tips so they never overlap the orange ones (they coincide at 90°)
        draw_frame_2d(ax, np.eye(2), COLORS["teal"], length=0.9 * lim, prime="", text_offset=(-0.09 * lim, -0.11 * lim))
        draw_frame_2d(ax, C, COLORS["orange"], length=0.9 * lim, prime="′")
        for j in range(2):  # dashed projections onto the new axes: the foot is x'_j e'_j
            foot = xp[j] * C[:, j]
            ax.plot([x[0], foot[0]], [x[1], foot[1]], ls="--", color=COLORS["orange"], lw=1.0)
        arrow2d(ax, (0, 0), x, color=COLORS["accent"], label="$\\mathbf{x}$", lw=2.6)
        ax.set_title(f"θ = {np.rad2deg(th):.0f}°: the arrow stays, the orange rulers turn", fontsize=10)
        idx = np.arange(2)
        ax2.bar(idx - 0.18, x, width=0.36, color=COLORS["teal"], label="$x_i$ (old frame)")
        ax2.bar(idx + 0.18, xp, width=0.36, color=COLORS["orange"], label="$x'_j = x_i C_{ij}$")
        ax2.set_xticks(idx, ["1", "2"])
        ax2.set_ylim(-lim, lim)
        ax2.axhline(0, color=COLORS["ink"], lw=0.8)
        ax2.set_ylabel("component")
        ax2.legend(loc="upper right", fontsize=8.5)
        txt = ("C = \n" + "\n".join("  ".join(f"{v:+.3f}" for v in row) for row in C)
               + f"\nmax|CᵀC − I| = {ch02.orthogonality_residual(C):.0e}   (2.6)"
               + f"\n|x| = {np.linalg.norm(x):.3f} = |x'| = {np.linalg.norm(xp):.3f}")
        ax2.text(0.02, 0.02, txt, transform=ax2.transAxes, fontsize=8.5, family="monospace", va="bottom")
        ax2.set_title("same arrow, different numbers", fontsize=10)
        return ()

    update(0)
    return fig, update


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch02"))
    ap.add_argument("--no-show", action="store_true")
    ap.add_argument("--angle-deg", type=float, default=40.0, help="rotation angle of the frame about (1,1,1)/sqrt(3)")
    args = ap.parse_args()
    t0 = time.perf_counter()
    import matplotlib

    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from ch02_drawings import arrow3d, draw_frame_3d, setup_3d
    from fluidpy import ch02_cartesian_tensors as ch02
    from fluidpy.core.style import COLORS, use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    axis = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)
    C = ch02.rotation_matrix_3d(axis, np.deg2rad(args.angle_deg))  # passive: columns = new axes e'_j
    E = ch02.unit_vectors(3)
    E_new = C.T  # rows = e'_j in old components
    assert np.allclose(ch02.direction_cosines(E, E_new), C)  # C_ij = e_i·e'_j
    x = np.array([1.0, 1.3, 0.9])
    xp = ch02.transform_vector(x, C)  # Eq. (2.5): x'_j = x_i C_ij
    x_back = ch02.inverse_transform_vector(xp, C)  # Eq. (2.7): x_j = x'_i C_ji
    assert np.allclose(ch02.vector_from_components(xp, E_new), x)  # Eq. (2.3): same x, primed basis

    fig = plt.figure(figsize=(11, 5))
    ax = fig.add_subplot(121, projection="3d")
    setup_3d(ax, lim=1.7)
    draw_frame_3d(ax, E, COLORS["teal"], length=1.4)
    draw_frame_3d(ax, E_new, COLORS["orange"], length=1.4, prime="′")
    arrow3d(ax, (0, 0, 0), x, color=COLORS["accent"], label="$\\mathbf{x}$", lw=2.5)
    ax.set_title(f"Fig. 2.2 idea: O123 (teal) and O1′2′3′ (orange, rotated {args.angle_deg:.0f}° about (1,1,1))")
    ax.view_init(elev=20, azim=-55)

    ax2 = fig.add_subplot(122)
    idx = np.arange(3)
    ax2.bar(idx - 0.18, x, width=0.36, color=COLORS["teal"], label="$x_i$ (old frame)")
    ax2.bar(idx + 0.18, xp, width=0.36, color=COLORS["orange"], label="$x'_j = x_i C_{ij}$ (new frame)")
    ax2.set_xticks(idx, ["1", "2", "3"])
    ax2.set_ylabel("component")
    ax2.set_title(f"same arrow, different numbers   |x| = {np.linalg.norm(x):.4f} = |x'| = {np.linalg.norm(xp):.4f}")
    ax2.legend()
    txt = "C = e_i·e'_j =\n" + "\n".join("  ".join(f"{v:6.3f}" for v in row) for row in C)
    ax2.text(0.02, 0.02, txt + f"\nmax|CᵀC − δ| = {ch02.orthogonality_residual(C):.1e}, det C = {np.linalg.det(C):+.3f}",
             transform=ax2.transAxes, fontsize=8.5, family="monospace", va="bottom")
    fig.savefig(out / "fig2_2_rotated_axes.png", bbox_inches="tight")

    # the notebook's rotatable version (storyboard A.2 #18): rotated_axes_figure(C, x=[1, 2, 3])
    fig_p = rotated_axes_figure(C, x=[1.0, 2.0, 3.0])
    fig_p.write_html(out / "fig2_2_rotated_axes.html", include_plotlyjs="cdn")
    xp3 = ch02.transform_vector([1.0, 2.0, 3.0], C)
    print(f"plotly rotated axes: x = (1, 2, 3) -> x' = {np.round(xp3, 4)}, |x| = {np.sqrt(14):.4f} = |x'| = {np.linalg.norm(xp3):.4f}, "
          f"{len(fig_p.data)} traces")

    # the C02 animation frames (last frame saved as a still): x = (1, 2) at θ = 0°, 30°, 60°, 90°
    fig_a, update = frame_rotation_frames(np.deg2rad([0.0, 30.0, 60.0, 90.0]))
    update(3)
    fig_a.savefig(out / "fig2_2_frame_rotation_90deg.png", bbox_inches="tight")
    C90 = ch02.rotation_matrix_2d(np.pi / 2)
    print(f"frame animation at 90°: x = (1, 2) -> x' = {np.round(ch02.transform_vector([1.0, 2.0], C90), 6)} (roles swap: x'_1 = x_2, x'_2 = -x_1)")

    print(f"C (passive, columns = new axes):\n{np.round(C, 4)}")
    print(f"x = {x} -> x' = Cᵀx = {np.round(xp, 4)} -> back = {np.round(x_back, 4)}")
    print(f"orthogonality residual {ch02.orthogonality_residual(C):.2e}, proper rotation: {ch02.is_proper_rotation(C)}, "
          f"angle recovered {np.rad2deg(ch02.rotation_angle(C)):.2f}°")
    print(f"saved fig2_2_rotated_axes.png in {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
