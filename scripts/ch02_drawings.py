"""Small matplotlib drawing helpers shared by the ch02 scripts (arrows in 2-D and 3-D, frames, cubes).

Nothing here is physics; the numbers come from ``fluidpy.ch02_cartesian_tensors``. Imported by ``scripts/ch02_*.py``
via ``sys.path`` (see the scripts' header).
"""
from __future__ import annotations

import numpy as np

from fluidpy.core.style import COLORS


def arrow2d(ax, start, vec, color=COLORS["accent"], label=None, lw=2.0, ls="-", text_offset=(0.05, 0.05), **kw):
    """Draw an arrow from ``start`` along ``vec`` and (optionally) label its tip."""
    s, v = np.asarray(start, float), np.asarray(vec, float)
    ax.annotate("", xy=s + v, xytext=s, arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, ls=ls, shrinkA=0, shrinkB=0), **kw)
    if label:
        ax.text(*(s + v + np.asarray(text_offset)), label, color=color, fontsize=11)


def arrow3d(ax, start, vec, color=COLORS["accent"], label=None, lw=2.0, ls="-", ratio=0.15):
    """3-D arrow (quiver) from ``start`` along ``vec``."""
    s, v = np.asarray(start, float), np.asarray(vec, float)
    ax.quiver(*s, *v, color=color, lw=lw, arrow_length_ratio=ratio, linestyle=ls)
    if label:
        ax.text(*(s + 1.08 * v), label, color=color, fontsize=11)


def setup_3d(ax, lim=1.6, labels=("$x_1$", "$x_2$", "$x_3$")):
    """Equal-aspect 3-D axes through the origin with light styling."""
    ax.set_xlim(-lim * 0.2, lim)
    ax.set_ylim(-lim * 0.2, lim)
    ax.set_zlim(-lim * 0.2, lim)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    ax.set_zlabel(labels[2])
    ax.grid(False)
    for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
        pane.fill = False
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])


def draw_frame_3d(ax, E, color, names=("1", "2", "3"), length=1.0, prime=""):
    """Draw the three unit vectors (rows of E) as arrows from the origin."""
    for e, nm in zip(E, names):
        arrow3d(ax, (0, 0, 0), length * np.asarray(e), color=color, label=f"{nm}{prime}", lw=1.6)


def draw_frame_2d(ax, C, color, length=1.0, prime="′", names=("1", "2"), text_offset=(0.02, 0.02)):
    """Draw the axes of the frame whose unit vectors are the columns of C."""
    for j in range(2):
        arrow2d(ax, (0, 0), length * C[:, j], color=color, label=f"{names[j]}{prime}", lw=1.6, text_offset=text_offset)


def cube_edges(center=(0, 0, 0), half=0.5):
    """Twelve edges of an axis-aligned cube as (start, end) pairs."""
    c = np.asarray(center, float)
    corners = np.array([[sx, sy, sz] for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)]) * half + c
    edges = []
    for i in range(8):
        for j in range(i + 1, 8):
            if np.sum(np.abs(corners[i] - corners[j]) > 1e-12) == 1:
                edges.append((corners[i], corners[j]))
    return edges


def draw_cube(ax, center=(0, 0, 0), half=0.5, color=COLORS["ink"], lw=1.0):
    for a, b in cube_edges(center, half):
        ax.plot(*zip(a, b), color=color, lw=lw, alpha=0.6)


# ----------------------------------------------------------------------------------------------------------------------
# plotly 3-D helpers (rotatable figures for the notebook; ``go.Cone(sizemode="raw")`` draws a cone whose length equals
# the vector length in axis units, so every arrow is one shaft line + exactly one cone)
# ----------------------------------------------------------------------------------------------------------------------
def plotly_arrow(start, vec, color, name=None, width=5, head=0.3, showlegend=False, hovertext=None):
    """One 3-D arrow as ``[go.Scatter3d shaft, go.Cone head]`` from ``start`` along ``vec`` (head = fraction of |vec|)."""
    import plotly.graph_objects as go

    s, v = np.asarray(start, float), np.asarray(vec, float)
    base = s + (1.0 - head) * v  # the cone's tail sits where the shaft ends
    shaft = go.Scatter3d(x=[s[0], base[0]], y=[s[1], base[1]], z=[s[2], base[2]], mode="lines",
                         line=dict(color=color, width=width), name=name, showlegend=showlegend, hoverinfo="skip",
                         legendgroup=name)
    cone = go.Cone(x=[base[0]], y=[base[1]], z=[base[2]], u=[head * v[0]], v=[head * v[1]], w=[head * v[2]],
                   sizemode="raw", anchor="tail", colorscale=[[0, color], [1, color]], showscale=False, name=name,
                   showlegend=False, hovertext=hovertext, hoverinfo="text" if hovertext else "skip", legendgroup=name)
    return [shaft, cone]


def plotly_cube_edges(center=(0, 0, 0), half=0.5, color=COLORS["ink"], width=2, name="cube"):
    """The twelve cube edges as a single ``go.Scatter3d`` line trace (``None`` breaks between edges)."""
    import plotly.graph_objects as go

    xs, ys, zs = [], [], []
    for a, b in cube_edges(center, half):
        xs += [a[0], b[0], None]
        ys += [a[1], b[1], None]
        zs += [a[2], b[2], None]
    return go.Scatter3d(x=xs, y=ys, z=zs, mode="lines", line=dict(color=color, width=width), name=name,
                        hoverinfo="skip", showlegend=False)


def plotly_layout_3d(fig, ranges, title="", height=520, labels=("x₁", "x₂", "x₃"), camera_eye=(1.6, -1.9, 1.1)):
    """Equal-aspect 3-D layout with the house look: ``ranges`` = ((x0, x1), (y0, y1), (z0, z1))."""
    axis = lambda r, lab: dict(range=list(r), title=lab, showbackground=False, gridcolor=COLORS["grid"],  # noqa: E731
                               zerolinecolor=COLORS["grid"])
    fig.update_layout(title=dict(text=title, font=dict(size=13)), height=height, margin=dict(l=0, r=0, t=60, b=0),
                      scene=dict(xaxis=axis(ranges[0], labels[0]), yaxis=axis(ranges[1], labels[1]),
                                 zaxis=axis(ranges[2], labels[2]), aspectmode="cube",
                                 camera=dict(eye=dict(x=camera_eye[0], y=camera_eye[1], z=camera_eye[2]))),
                      legend=dict(orientation="h", y=-0.02), paper_bgcolor="white")
    return fig
