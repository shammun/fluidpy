"""Chapter 2, §2.4, Fig. 2.4 (our drawing): the stress cube with the nine components on the three +faces (and the
reversed arrows on the hidden −faces), read off with ``cube_face_tractions`` and the sign convention.

Run: ``.venv/Scripts/python.exe scripts/ch02_fig2_4_stress_cube.py --no-show``
Figure → outputs/ch02/fig2_4_stress_cube.png.
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


def stress_cube_figure(tau, show_hidden: bool = False, half: float = 1.0, arrow: float = 0.85, height: int = 540):
    """Rotatable plotly stress cube: on the +x_i face the traction is row i of τ, one arrow per non-zero component.

    Book: §2.4 (Fig. 2.4 idea and the sign convention "first index = face, second = direction of the force"; the
    face letters A–H are the book's, reused by Example 2.5). Drawing only — the (normal, traction) pairs come from
    ``ch02.cube_face_tractions(tau)``.

    Parameters
    ----------
    tau : (3, 3) stress tensor [Pa]
    show_hidden : also draw the three −faces (reversed arrows, grey); 9 cones without, 18 with
    half : half-side of the cube [drawing units]
    arrow : length of the longest arrow [drawing units] (arrows scale with |τ_ij| / max|τ|)
    height : figure height [px]

    Returns
    -------
    go.Figure with ≤ 18 ``go.Cone`` traces: normal stresses (i = j) blue, shears rose, hidden faces grey; a negative
    normal stress points *into* the cube (compression, Ch. 1's pressure).
    """
    import plotly.graph_objects as go

    from ch02_drawings import plotly_arrow, plotly_cube_edges, plotly_layout_3d
    from fluidpy import ch02_cartesian_tensors as ch02
    from fluidpy.core.style import COLORS

    tau = np.asarray(tau, dtype=float)
    faces = ch02.cube_face_tractions(tau)  # {"+1": {"normal", "traction", "letters"}, …}
    scale = arrow / max(1e-300, float(np.max(np.abs(tau))))
    fig = go.Figure()
    fig.add_trace(plotly_cube_edges(half=half, color=COLORS["ink"]))
    for key, f in faces.items():
        plus = key.startswith("+")
        if not plus and not show_hidden:
            continue
        i = int(key[1]) - 1
        n = np.asarray(f["normal"], float)
        centre = half * n
        # a faint face with its book letters on hover
        corners = np.array([centre + half * (a * np.roll(np.eye(3)[i], 1) + b * np.roll(np.eye(3)[i], 2))
                            for a, b in ((-1, -1), (1, -1), (1, 1), (-1, 1))])
        fig.add_trace(go.Mesh3d(x=corners[:, 0], y=corners[:, 1], z=corners[:, 2], i=[0, 0], j=[1, 2], k=[2, 3],
                                color=COLORS["teal"] if plus else COLORS["muted"], opacity=0.12 if plus else 0.05,
                                hovertext=f"face {f['letters']}: outward normal {'+' if plus else '-'}x{i + 1}",
                                hoverinfo="text", name=f"face {key} ({f['letters']})", showlegend=False))
        for j in range(3):
            comp = float(f["traction"][j])  # ±τ_ij: the force per unit area along x_j on this face
            if abs(comp) < 1e-12:
                continue
            vec = np.zeros(3)
            vec[j] = scale * comp
            colour = (COLORS["blue"] if i == j else COLORS["rose"]) if plus else COLORS["muted"]
            label = f"τ{i + 1}{j + 1} = {tau[i, j]:+g} Pa" + ("" if plus else " (reversed on the −face)")
            fig.add_traces(plotly_arrow(centre, vec, colour, name=label, width=6 if plus else 3,
                                        hovertext=f"{label}: face {f['letters']} ({'+' if plus else '-'}x{i + 1}), "
                                                  f"direction {'+' if comp > 0 else '-'}x{j + 1}"))
            if plus:
                # label at the arrow tip when it points away from the cube, else just outside the face (inward arrows
                # end inside the cube, where a label would be hidden)
                outward = float(vec @ n) > 0
                tip = centre + 1.15 * vec if outward else centre + 0.18 * half * n + 0.35 * vec
                fig.add_trace(go.Scatter3d(x=[tip[0]], y=[tip[1]], z=[tip[2]], mode="text",
                                           text=[f"τ<sub>{i + 1}{j + 1}</sub>"], textfont=dict(color=colour, size=13),
                                           hoverinfo="skip", showlegend=False))
    lim = half + arrow + 0.15 * half
    title = ("Fig. 2.4 idea: row i of τ acts on the +x<sub>i</sub> face<br>blue = normal stress (negative points in), "
             "rose = shear" + (", grey = reversed arrows on the −faces" if show_hidden else ""))
    return plotly_layout_3d(fig, ((-lim, lim),) * 3, title=title, height=height, camera_eye=(1.7, 1.4, 1.0))


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

    from ch02_drawings import arrow3d, draw_cube, setup_3d
    from fluidpy import ch02_cartesian_tensors as ch02
    from fluidpy.core.style import COLORS, use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    p, a = 1.0, 0.6  # a compressive pressure p plus a shear a [Pa] — the ch04 form τ = −pδ + viscous part
    tau = np.array([[-p, a, 0.0], [a, -p, 0.0], [0.0, 0.0, -p]])
    faces = ch02.cube_face_tractions(tau)
    assert np.allclose(sum(f["traction"] for f in faces.values()), 0.0)  # opposite faces cancel (Newton III at a point)

    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")
    setup_3d(ax, lim=1.0)
    half = 0.7
    for lim in (ax.set_xlim, ax.set_ylim, ax.set_zlim):
        lim(-1.05, 1.05)
    draw_cube(ax, half=half)
    scale = 0.5 / max(1.0, np.max(np.abs(tau)))
    for key, f in faces.items():
        n = f["normal"]
        centre = half * n
        i = int(key[1]) - 1
        plus = key.startswith("+")
        for j in range(3):
            comp = np.zeros(3)
            comp[j] = f["traction"][j]
            if abs(comp[j]) > 1e-12:
                col = (COLORS["blue"] if i == j else COLORS["orange"]) if plus else COLORS["muted"]
                arrow3d(ax, centre, scale * comp, color=col, label=(f"$\\tau_{{{i + 1}{j + 1}}}$" if plus else None),
                        lw=2.0 if plus else 1.0)
    fig.suptitle("Fig. 2.4 idea: on the +$x_i$ face the traction is row i of τ — first index = face, second = direction\n"
                 "blue = normal (here compressive, −p), orange = shear; grey = the −faces carry the reversed arrows", fontsize=10)
    ax.view_init(elev=18, azim=-58)
    fig.savefig(out / "fig2_4_stress_cube.png", bbox_inches="tight")

    # the notebook's rotatable version (one cone per arrow)
    import plotly.graph_objects as go

    for hidden in (False, True):
        fig_p = stress_cube_figure(tau, show_hidden=hidden)
        n_cones = sum(isinstance(t, go.Cone) for t in fig_p.data)
        assert n_cones <= 18, n_cones
        fig_p.write_html(out / f"fig2_4_stress_cube{'_hidden' if hidden else ''}.html", include_plotlyjs="cdn")
        print(f"plotly stress cube (show_hidden={hidden}): {n_cones} cones, {len(fig_p.data)} traces")

    print(f"τ =\n{tau}")
    for key, f in faces.items():
        print(f"face {key} ({f['letters']}): outward normal {f['normal']}, traction {np.round(f['traction'], 3)} Pa")
    print(ch02.stress_component_meaning(2, 3))
    print(f"saved fig2_4_stress_cube.png in {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
