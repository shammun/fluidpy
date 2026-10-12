"""Schematic drawings and composite figures for chapter 14 — drawn by this code, in our own words; no physics is
computed here beyond calls to the tested functions of ``fluidpy.ch14_aerodynamics``.

Sketches (each takes an optional matplotlib Axes, keyword ``ax``; ``None`` makes a new figure; each returns the Figure):

* ``draw_planform_labelled``     a tapered planform seen from above with span s, chord c(y), root, tip, area A
* ``draw_three_axes``            pitch, roll and yaw axes of an aircraft
* ``draw_section_labelled``      a section with chord line, camber line, camber, thickness and angle of attack
* ``draw_axes_conventions``      the two sets of axes of the chapter: section plane (x chord, y up) and wing (y span, w down)
* ``draw_kelvin_circuit``        the material circuit round the foil and its starting vortex
* ``draw_triangle_OQP``          the triangle O–Q–P of the raised circle (cosine rule behind (14.6))
* ``draw_circle_angles``         the circle with Q, B, α, β and the stream-aligned diameter (the θ of (14.10))
* ``draw_tip_flow_rear_view``    flow round the tips seen from behind: outward below, inward above
* ``draw_horseshoe_staircase``   N nested horseshoe vortices and the staircase Γ(y) they make
* ``draw_force_triangle``        U, w, the tilted stream, dL_e normal to it, split into lift and induced drag
* ``draw_span_semicircle``       the span angle γ: y = −(s/2) cos γ
* ``draw_stalled_section``       a stalled section: separation near the nose, a wide wake
* ``draw_velocity_triangle``     U, V and the relative velocity with lift and its thrust part ("fish", "bird_down", "bird_up", "sail")
* ``draw_sail_polygon``          the sail's lift and drag resolved into drive and lateral force

Composite figures (each returns the Figure): ``fig_four_circles``, ``fig_three_stages``, ``fig_pressure_two_ways``,
``fig_planforms_loadings``.

The notebook imports them with ``sys.path.insert(0, "scripts"); from ch14_drawings import *``.  Colours follow the
chapter code: lift teal, drag rose, thrust / drive purple (accent), circulation amber, pressure / suction orange, stream
and velocities blue, references muted.  Γ is always the clockwise circulation.

Run: ``.venv/Scripts/python.exe scripts/ch14_drawings.py --no-show``   Figures -> outputs/ch14/drawings*.png
"""
from __future__ import annotations

import numpy as np
from ch14_common import (C_CIRC, C_DRAG, C_DRIVE, C_FLOW, C_LIFT, C_PRES, C_REF, COLORS, finish, parse_args, save, setup,
                         streamlines)

from fluidpy import ch14_aerodynamics as ch14

__all__ = ["draw_planform_labelled", "draw_three_axes", "draw_section_labelled", "draw_axes_conventions",
           "draw_kelvin_circuit", "draw_triangle_OQP", "draw_circle_angles", "draw_tip_flow_rear_view",
           "draw_horseshoe_staircase", "draw_force_triangle", "draw_span_semicircle", "draw_stalled_section",
           "draw_velocity_triangle", "draw_sail_polygon", "fig_four_circles", "fig_three_stages",
           "fig_pressure_two_ways", "fig_planforms_loadings"]

INK = COLORS["ink"]


def _ax(ax, figsize=(5.6, 4.2)):
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    return ax.figure, ax


def _arrow(ax, tail, head, color, lw: float = 2.0, style: str = "->"):
    ax.annotate("", head, tail, arrowprops=dict(arrowstyle=style, color=color, lw=lw, shrinkA=0, shrinkB=0))


def _bare(ax, xlim, ylim, title: str):
    ax.set(xlim=xlim, ylim=ylim, aspect="equal", title=title)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(False)


def _arc(ax, c, r, a0, a1, color, lw=1.2):
    t = np.linspace(a0, a1, 40)
    ax.plot(c[0] + r * np.cos(t), c[1] + r * np.sin(t), color=color, lw=lw)


def _section(e=None, beta=None, n=400, b=1.0):
    I = ch14.illustrative_inputs()
    A = ch14.zhukhovsky_airfoil(b, I["e"] if e is None else e, I["beta"] if beta is None else beta, n)
    return A["x"], A["y"]


# ======================================================================================================================
# sketches
# ======================================================================================================================
def draw_planform_labelled(ax=None):
    """A tapered planform from above: span s, chord c(y), root, tips, planform area A, Λ = s²/A (14.1)."""
    fig, ax = _ax(ax, (6.4, 3.6))
    I = ch14.illustrative_inputs()
    s = I["s"]
    P = ch14.planform("tapered", s, 2 * I["A"] / (s * (1 + I["taper"])), I["taper"])
    y, c = P["y"], P["c"]
    ax.fill(np.concatenate([y, y[::-1]]), np.concatenate([0.25 * c, -0.75 * c[::-1]]), color="#d5d9e4", ec=INK, lw=1.2)
    _arrow(ax, (-s / 2, 1.1), (s / 2, 1.1), INK, 1.2, "<->")
    ax.text(0, 1.3, f"wingspan s = {s:g} m", ha="center", fontsize=9)
    yc = 0.3 * s
    cc = float(P["chord"](yc))
    _arrow(ax, (yc, 0.25 * cc), (yc, -0.75 * cc), C_FLOW, 1.4, "<->")
    ax.text(yc + 0.15, -0.3, "chord c(y)", color=C_FLOW, fontsize=9)
    ax.text(0, -0.35, f"planform area A = {P['A']:g} m²", ha="center", fontsize=9)
    ax.text(0, -1.75, "root", ha="center", fontsize=8, color=C_REF)
    ax.text(s / 2, -0.95, "tip", ha="center", fontsize=8, color=C_REF)
    ax.text(-s / 2, -0.95, "tip", ha="center", fontsize=8, color=C_REF)
    _arrow(ax, (-0.42 * s, 2.3), (-0.42 * s, 1.5), C_FLOW, 1.6)
    ax.text(-0.42 * s + 0.2, 2.0, "stream U", color=C_FLOW, fontsize=8)
    _bare(ax, (-0.56 * s, 0.56 * s), (-2.1, 2.6), rf"planform: $\Lambda = s^2/A$ = {P['AR']:.1f}, $\bar c = A/s$ = {P['c_mean']:.2f} m")
    return fig


def draw_three_axes(ax=None):
    """Pitch, roll and yaw axes of an aircraft (a plan-view silhouette with the three axes named)."""
    fig, ax = _ax(ax, (5.6, 4.4))
    ax.fill([-2.2, 1.6, 2.0, 1.6, -2.2], [0.12, 0.18, 0, -0.18, -0.12], color="#d5d9e4", ec=INK)          # fuselage
    ax.fill([0.5, -0.3, -0.6, 0.1], [0.15, 2.2, 2.2, 0.15], color="#d5d9e4", ec=INK)                       # wings
    ax.fill([0.5, -0.3, -0.6, 0.1], [-0.15, -2.2, -2.2, -0.15], color="#d5d9e4", ec=INK)
    ax.fill([-1.8, -2.2, -2.35, -2.0], [0.1, 0.8, 0.8, 0.1], color="#d5d9e4", ec=INK)                      # tail
    ax.fill([-1.8, -2.2, -2.35, -2.0], [-0.1, -0.8, -0.8, -0.1], color="#d5d9e4", ec=INK)
    _arrow(ax, (-2.8, 0), (3.0, 0), C_LIFT, 1.8)
    ax.text(3.05, 0.12, "roll axis\n(nose to tail)", color=C_LIFT, fontsize=8)
    _arrow(ax, (0, -2.6), (0, 2.9), C_DRIVE, 1.8)
    ax.text(0.12, 2.75, "pitch axis (tip to tip)", color=C_DRIVE, fontsize=8)
    ax.plot([0], [0], "o", color=C_CIRC, ms=9)
    ax.plot([0], [0], "o", color="white", ms=3)
    ax.plot([0.1, 0.74], [-0.1, -0.66], color=C_CIRC, lw=0.8)                                             # leader to the label
    ax.text(0.78, -1.05, "yaw axis\n(out of the page)", color=C_CIRC, fontsize=8)
    ax.text(-3.0, -2.75, "ailerons roll · elevator pitches · rudder yaws", fontsize=8, color=C_REF, ha="left", va="top", bbox=dict(fc="white", ec="none", pad=1.5))
    _bare(ax, (-3.1, 4.6), (-3.2, 3.2), "three aircraft-fixed axes")
    return fig


def draw_section_labelled(x=None, y=None, geom=None, ax=None):
    """A section with its chord line, camber line, camber, maximum thickness and the angle of attack.

    ``x``, ``y``: contour (default: our Zhukhovsky section with b = 1); ``geom``: the dict of ``ch14.section_geometry``
    (computed if None)."""
    fig, ax = _ax(ax, (7.0, 3.4))
    if x is None or y is None:
        x, y = _section(e=0.12, beta=0.14)
    if geom is None:
        geom = ch14.section_geometry(x, y)
    le, te, c = geom["le"], geom["te"], geom["chord"]
    d = np.exp(1j * geom["chord_angle"])
    ax.fill(x, y, color="#e9ecf3", ec=INK, lw=1.3)
    ax.plot([le.real, te.real], [le.imag, te.imag], color=C_REF, lw=1.2, ls="--")
    xc, yc = geom["camber_line"]
    zc = le + (xc * c + 1j * yc) * d
    ax.plot(zc.real, zc.imag, color=C_CIRC, lw=1.6)
    k = int(np.argmax(np.abs(yc)))
    foot = le + xc[k] * c * d
    _arrow(ax, (foot.real, foot.imag), (zc[k].real, zc[k].imag), C_CIRC, 1.2, "<->")
    ax.text(zc[k].real + 0.03 * c, 0.5 * (foot.imag + zc[k].imag), "camber", color=C_CIRC, fontsize=8)
    kt = int(np.argmin(np.abs(geom["xc"] - geom["x_tmax"])))
    pu = le + (geom["xc"][kt] * c + 1j * geom["y_upper"][kt]) * d
    pl = le + (geom["xc"][kt] * c + 1j * geom["y_lower"][kt]) * d
    _arrow(ax, (pl.real, pl.imag), (pu.real, pu.imag), C_DRIVE, 1.2, "<->")
    ax.text(pu.real - 0.02 * c, pu.imag + 0.03 * c, r"thickness $t_{\max}$", color=C_DRIVE, fontsize=8, ha="center")
    ax.text(0.5 * (le.real + te.real) + 0.18 * c, le.imag - 0.07 * c, "chord line (length c)", color=C_REF, fontsize=8)
    kc = int(np.argmin(np.abs(xc - 0.22)))
    ax.annotate("camber line", (zc[kc].real, zc[kc].imag), (zc[kc].real - 0.02 * c, le.imag - 0.12 * c), color=C_CIRC, fontsize=8,
                ha="center", arrowprops=dict(arrowstyle="-", color=C_CIRC, lw=0.8, shrinkA=1, shrinkB=0))
    al = np.deg2rad(12.0)
    p0 = le - 0.42 * c * np.exp(1j * al)
    _arrow(ax, (p0.real, p0.imag), (le.real - 0.05 * c, le.imag - 0.05 * c * np.tan(al)), C_FLOW, 1.8)
    ax.plot([le.real - 0.42 * c, le.real], [le.imag, le.imag], color=C_REF, lw=0.8, ls=":")
    _arc(ax, (le.real, le.imag), 0.3 * c, np.pi, np.pi + al, C_FLOW)
    ax.text(le.real - 0.37 * c, le.imag - 0.045 * c, "α", color=C_FLOW, fontsize=10)
    ax.text(p0.real - 0.02 * c, p0.imag - 0.05 * c, "stream U", color=C_FLOW, fontsize=8)
    ax.text(le.real - 0.02 * c, le.imag + 0.06 * c, "leading edge\n(rounded)", fontsize=7, ha="right", color=C_REF)
    ax.text(te.real, te.imag + 0.05 * c, "trailing edge\n(sharp)", fontsize=7, ha="center", color=C_REF)
    _bare(ax, (le.real - 0.5 * c, te.real + 0.12 * c), (le.imag - 0.22 * c, le.imag + 0.3 * c), "airfoil section: the vocabulary")
    return fig


def draw_axes_conventions(ax=None):
    """The two sets of axes of the chapter side by side: the section plane (x along the mapped chord, y up, stream at α)
    and the wing (y along the span, stream along x, downwash w positive downward)."""
    fig, ax = _ax(ax, (8.0, 3.6))
    x, y = _section()
    ax.fill(0.45 * x, 0.45 * y, color="#e9ecf3", ec=INK)
    _arrow(ax, (-1.5, 0), (1.6, 0), INK, 1.2)
    _arrow(ax, (0, -0.9), (0, 1.3), INK, 1.2)
    ax.text(1.62, -0.05, "x (mapped chord)", fontsize=8, va="top", ha="right")
    ax.text(0.06, 1.25, "y (up)", fontsize=8)
    al = np.deg2rad(14)
    _arrow(ax, (-1.6, -0.75), (-1.6 + 0.8 * np.cos(al), -0.75 + 0.8 * np.sin(al)), C_FLOW, 2)
    ax.text(-1.55, -0.95, "stream U at angle α", color=C_FLOW, fontsize=8)
    _arrow(ax, (0.1, 0.15), (0.1 - 0.9 * np.sin(al), 0.15 + 0.9 * np.cos(al)), C_LIFT, 2)
    ax.text(-0.3, 0.72, "lift normal to the\nstream (not along y)", color=C_LIFT, fontsize=8, ha="right")
    ax.text(-0.2, -1.25, "§14.3–§14.5: section plane", fontsize=9, ha="center")
    ox = 4.2
    ax.plot([ox - 1.4, ox + 1.4], [0, 0], color=INK, lw=4)
    _arrow(ax, (ox - 1.7, 0), (ox + 1.8, 0), INK, 1.0)
    ax.text(ox + 1.85, -0.1, "y (span)", fontsize=8, ha="right", va="top")
    _arrow(ax, (ox, 0.0), (ox, 1.2), INK, 1.0)
    ax.text(ox + 0.06, 1.15, "z (up)", fontsize=8)
    for yy in (-1.0, -0.5, 0.0, 0.5, 1.0):
        _arrow(ax, (ox + yy, -0.08), (ox + yy, -0.55), C_FLOW, 1.4)
    ax.text(ox, -0.8, "downwash w (positive down)", color=C_FLOW, fontsize=8, ha="center")
    ax.text(ox - 1.4, 0.15, "−s/2 port", fontsize=7, ha="center", color=C_REF)
    ax.text(ox + 1.45, 0.15, "+s/2 starboard", fontsize=7, ha="right", color=C_REF)
    ax.text(ox, -1.25, "§14.6: wing seen from behind (stream into the page)", fontsize=9, ha="center")
    _bare(ax, (-1.9, 6.2), (-1.45, 1.5), "axes change between sections")
    return fig


def draw_kelvin_circuit(ax=None):
    """The material circuit ABCD round the foil and its starting vortex, split by BD: the two circulations add to zero."""
    fig, ax = _ax(ax, (7.4, 3.4))
    x, y = _section()
    ax.fill(0.35 * x - 2.3, 0.35 * y, color="#e9ecf3", ec=INK)
    ax.plot([-3.6, 3.6, 3.6, -3.6, -3.6], [-1.3, -1.3, 1.3, 1.3, -1.3], color=C_CIRC, lw=2)
    ax.plot([0.4, 0.4], [-1.3, 1.3], color=C_CIRC, lw=1.2, ls="--")
    for name, p in zip("ABCD", ((-3.6, 1.3), (0.4, 1.3), (3.6, -1.3), (0.4, -1.3))):
        ax.text(p[0] + 0.1, p[1] + (0.12 if p[1] > 0 else -0.3), name, color=C_CIRC, fontsize=10)
    ax.text(3.7, 1.35, "(circuit ABCD: A–B–corner–C–D–corner)", fontsize=6, color=C_REF, ha="right", va="bottom")
    _arc(ax, (-2.3, 0), 0.95, np.deg2rad(160), np.deg2rad(20), C_LIFT, 2)
    _arrow(ax, (-1.43, 0.42), (-1.39, 0.3), C_LIFT, 2)
    ax.text(-2.3, 0.75, "Γ clockwise (bound)", color=C_LIFT, fontsize=8, ha="center")
    t = np.linspace(0, 3.2 * np.pi, 200)
    r = 0.07 + 0.045 * t
    ax.plot(2.3 + r * np.cos(t), r * np.sin(t), color=C_DRAG, lw=1.6)
    ax.text(2.3, 0.75, "starting vortex: −Γ\n(counterclockwise)", color=C_DRAG, fontsize=8, ha="center")
    ax.plot(np.linspace(-1.55, 1.7, 30), 0.05 * np.sin(np.linspace(0, 9, 30)), color=C_REF, lw=0.8)
    ax.text(0, -1.75, r"$\Gamma_{ABCD}=0=\Gamma_{\rm foil}+\Gamma_{\rm starting}$ (Kelvin: the circuit started in fluid at rest)", ha="center", fontsize=9)
    _bare(ax, (-3.9, 3.9), (-2.0, 1.7), "material circuit round foil and starting vortex")
    return fig


def draw_triangle_OQP(ax=None, beta: float = 0.35, theta: float = 1.1):
    """The triangle O–Q–P of the raised circle (centre Q on the η-axis, through ζ = ±b): the cosine rule behind (14.6)."""
    fig, ax = _ax(ax, (5.2, 4.6))
    b = 1.0
    mu, a = ch14.zhukhovsky_circle(b, 0.0, beta)
    t = np.linspace(0, 2 * np.pi, 300)
    ax.plot(mu.real + a * np.cos(t), mu.imag + a * np.sin(t), color=C_FLOW)
    # P on the circle at polar angle theta from O: R from R² − 2R b tanβ sinθ − b² = 0
    R = b * np.tan(beta) * np.sin(theta) + np.sqrt((b * np.tan(beta) * np.sin(theta)) ** 2 + b ** 2)
    P = R * np.exp(1j * theta)
    ax.plot([0, mu.real, P.real, 0], [0, mu.imag, P.imag, 0], color=C_CIRC, lw=1.8)
    ax.plot([-b, 0, b], [0, 0, 0], "o", color=INK, ms=4)
    ax.plot([-b, mu.real], [0, mu.imag], color=C_REF, lw=1)
    ax.axhline(0, color=C_REF, lw=0.8)
    ax.axvline(0, color=C_REF, lw=0.8)
    for lab, p, off in (("O", 0j, (-0.13, -0.16)), ("Q", mu, (-0.17, 0.02)), ("P", P, (0.05, 0.05)), ("A (−b)", -b + 0j, (-0.3, -0.17)), ("B (b)", b + 0j, (0.02, -0.17))):
        ax.text(p.real + off[0], p.imag + off[1], lab, fontsize=9)
    ax.text(0.66 * P.real + 0.07, 0.66 * P.imag - 0.06, "OP = R", color=C_CIRC, fontsize=9)
    m = 0.5 * (mu + P)
    ax.text(m.real - 0.07, m.imag + 0.05, "QP = a = b/cos β", color=C_CIRC, fontsize=8, ha="right")
    ax.annotate("OQ = b tan β", (0.0, 0.5 * mu.imag), (-0.62, 0.66), color=C_CIRC, fontsize=8, ha="center",
                arrowprops=dict(arrowstyle="-", color=C_CIRC, lw=0.8, shrinkA=1, shrinkB=1))
    _arc(ax, (0, 0), 0.3, 0, theta, INK)
    ax.text(0.33, 0.12, "θ", fontsize=10)
    _arc(ax, (-b, 0), 0.45, 0, beta, INK)
    ax.text(-b + 0.5, 0.05, "β", fontsize=10)
    _arc(ax, (0, 0), 0.18, theta, np.pi / 2, C_DRAG)
    am = 0.5 * (theta + np.pi / 2)
    ax.annotate("90° − θ", (0.18 * np.cos(am), 0.18 * np.sin(am)), (0.62, 0.36), color=C_DRAG, fontsize=8, va="center",
                arrowprops=dict(arrowstyle="-", color=C_DRAG, lw=0.8, shrinkA=1, shrinkB=0))
    _bare(ax, (-1.5, 1.5), (-0.95, 1.75), "triangle OQP:\nQP² = OP² + OQ² − 2·OQ·OP·cos(90° − θ)")
    return fig


def draw_circle_angles(ax=None, alpha: float = np.deg2rad(18.0), beta: float = np.deg2rad(14.0), e: float = 0.12):
    """The off-centre circle with its centre Q, the trailing-edge point B, the angles α and β and the stream-aligned
    diameter from which the θ of u_θ = −2U sin θ − Γ/(2πa) (14.10) is measured: B sits at θ = −(α + β)."""
    fig, ax = _ax(ax, (5.4, 4.8))
    b = 1.0
    mu, a = ch14.zhukhovsky_circle(b, e, beta)
    t = np.linspace(0, 2 * np.pi, 300)
    ax.plot(mu.real + a * np.cos(t), mu.imag + a * np.sin(t), color=C_FLOW, lw=1.8)
    ax.axhline(0, color=C_REF, lw=0.8)
    ax.axvline(0, color=C_REF, lw=0.8)
    d = a * np.exp(1j * alpha)
    ax.plot([(mu - 1.25 * d).real, (mu + 1.25 * d).real], [(mu - 1.25 * d).imag, (mu + 1.25 * d).imag], color=C_DRIVE, lw=1.4, ls="--")
    ax.text((mu + 1.28 * d).real - 0.25, (mu + 1.28 * d).imag + 0.06, "stream-aligned diameter (θ = 0)", color=C_DRIVE, fontsize=7)
    ax.plot([mu.real, b], [mu.imag, 0], color=C_CIRC, lw=1.6)
    ax.plot([mu.real, mu.real + a], [mu.imag, mu.imag], color=C_REF, lw=0.9, ls=":")
    ax.plot([mu.real, b], [mu.imag, 0], "o", color=INK, ms=5)
    ax.text(mu.real - 0.14, mu.imag + 0.06, "Q", fontsize=10)
    ax.text(b + 0.05, -0.14, "B (ζ = b)", fontsize=9)
    _arc(ax, (mu.real, mu.imag), 0.42, 0, alpha, C_DRIVE)
    ax.text(mu.real + 0.47, mu.imag + 0.06, "α", color=C_DRIVE, fontsize=10)
    _arc(ax, (mu.real, mu.imag), 0.55, -beta, 0, C_CIRC)
    ax.text(mu.real + 0.6, mu.imag - 0.13, "β", color=C_CIRC, fontsize=10)
    p0 = mu - 2.0 * d
    _arrow(ax, (p0.real, p0.imag), ((p0 + 0.55 * d).real, (p0 + 0.55 * d).imag), C_FLOW, 2)
    ax.text(p0.real + 0.05, p0.imag - 0.16, "U", color=C_FLOW, fontsize=10)
    ax.text(mu.real - 0.2, mu.imag - a - 0.3, "at B: θ = −(α + β)  ⇒  Γ = 4πUa sin(α + β)  (14.11), clockwise", fontsize=8, ha="center")
    _bare(ax, (-2.6, 2.0), (-1.75, 1.9), "angles on the circle (ζ-plane)")
    return fig


def draw_tip_flow_rear_view(ax=None):
    """Rear view of a lifting wing: high pressure below, low above; the air moves outward below, inward above and round
    the tips, leaving streamwise vorticity of opposite sign behind the two half-wings."""
    fig, ax = _ax(ax, (6.6, 3.4))
    ax.plot([-2, 2], [0, 0], color=INK, lw=5)
    ax.text(0, 0.42, "low pressure (suction)", color=C_PRES, ha="center", fontsize=8)
    ax.text(0, -0.55, "high pressure", color=C_FLOW, ha="center", fontsize=8)
    for x0, sg in ((0.5, 1), (1.2, 1), (-0.5, -1), (-1.2, -1)):
        _arrow(ax, (x0, -0.2), (x0 + 0.5 * sg, -0.2), C_FLOW, 1.4)
        _arrow(ax, (x0 + 0.5 * sg, 0.2), (x0, 0.2), C_PRES, 1.4)
    for sg in (1, -1):
        t = np.linspace(-0.5 * np.pi, 0.5 * np.pi, 40)
        ax.plot(sg * (2.0 + 0.3 * np.cos(t)), 0.3 * np.sin(t), color=C_CIRC, lw=1.8)
        _arrow(ax, (sg * 2.06, 0.295), (sg * 1.96, 0.3), C_CIRC, 1.8)
        ax.text(sg * 2.75, 0.55, "tip vortex\n" + ("counterclockwise" if sg > 0 else "clockwise"), color=C_CIRC, fontsize=8, ha="center")
    for x0 in (-1.5, -0.75, 0.0, 0.75, 1.5):
        _arrow(ax, (x0, -0.85), (x0, -1.25), C_REF, 1.2)
    ax.text(0, -1.5, "between the tips the air sinks: downwash", color=C_REF, ha="center", fontsize=8)
    ax.text(-2.0, -0.22, "port", fontsize=7, ha="center", color=C_REF)
    ax.text(2.0, -0.22, "starboard", fontsize=7, ha="center", color=C_REF)
    _bare(ax, (-3.5, 3.5), (-1.7, 1.1), "flow round the wing tips (seen from behind)")
    return fig


def draw_horseshoe_staircase(N: int = 7, ax=None):
    """N nested horseshoe vortices (top view) and, above them, the staircase Γ(y) they add up to against the elliptic
    loading; the trailing filament at each step carries the drop of Γ there."""
    fig, ax = _ax(ax, (6.4, 4.6))
    s = 2.0
    ye, G = ch14.staircase_from_loading(lambda q: ch14.elliptic_circulation(q, s, 1.0), s, int(N))
    jump = np.concatenate([[0.0], G]) - np.concatenate([G, [0.0]])
    for yk, jk in zip(ye, jump):
        ax.plot([yk, yk], [0, -1.6], color=C_CIRC, lw=0.6 + 4.0 * abs(jk))
    ax.plot([-1, 1], [0, 0], color=INK, lw=4)
    ax.text(1.05, 0.02, "lifting line", fontsize=8)
    ax.text(0, -1.8, "trailing filaments: strength = the drop of Γ at each step, −(dΓ/dy)dy", color=C_CIRC, fontsize=8, ha="center")
    yy = np.linspace(-1, 1, 200)
    ax.plot(yy, 0.4 + 1.1 * ch14.elliptic_circulation(yy, s, 1.0), color=C_REF, ls="--", lw=1.2)
    ax.step(ye, 0.4 + 1.1 * np.append(G, G[-1]), where="post", color=C_CIRC, lw=1.8)
    ax.plot([ye[0], ye[0]], [0.4, 0.4 + 1.1 * G[0]], color=C_CIRC, lw=1.8)
    ax.plot([ye[-1], ye[-1]], [0.4, 0.4 + 1.1 * G[-1]], color=C_CIRC, lw=1.8)
    ax.plot([-1, 1], [0.4, 0.4], color=C_REF, lw=0.8)
    ax.text(-1.0, 1.55, "Γ(y): staircase and smooth limit", fontsize=8, color=C_CIRC)
    _arrow(ax, (-1.35, -0.4), (-1.35, -1.1), C_FLOW, 1.6)
    ax.text(-1.3, -0.3, "U", color=C_FLOW, fontsize=9)
    _bare(ax, (-1.6, 1.7), (-2.0, 1.75), f"{int(N)} horseshoe vortices side by side")
    return fig


def draw_force_triangle(ax=None, eps_deg: float = 12.0):
    """The stream U and the downwash w give the tilted stream U_e (angle ε); the force dL_e = ρU_eΓ dy is normal to U_e
    and splits into the lift dL (normal to U) and the induced drag dD_i (along U)."""
    fig, ax = _ax(ax, (6.6, 4.0))
    eps = np.deg2rad(eps_deg)
    U = 2.4
    _arrow(ax, (-3.2, 1.2), (-3.2 + U, 1.2), C_FLOW, 2)
    _arrow(ax, (-3.2 + U, 1.2), (-3.2 + U, 1.2 - U * np.tan(eps)), C_FLOW, 1.6)
    _arrow(ax, (-3.2, 1.2), (-3.2 + U, 1.2 - U * np.tan(eps)), C_DRIVE, 2)
    ax.text(-2.1, 1.3, "U", color=C_FLOW, fontsize=10)
    ax.text(-0.72, 0.9, "w", color=C_FLOW, fontsize=10)
    ax.text(-2.3, 0.55, r"$U_e$", color=C_DRIVE, fontsize=10)
    _arc(ax, (-3.2, 1.2), 1.1, -eps, 0, INK)
    ax.text(-2.0, 0.98, "ε", fontsize=10)
    x, y = _section()
    z = (0.3 * (x + 1j * y)) * np.exp(-1j * np.deg2rad(8)) + 1.6 - 0.2j
    ax.fill(z.real, z.imag, color="#e9ecf3", ec=INK)
    o = (1.5, 0.0)
    Le = 2.2
    _arrow(ax, o, (o[0] + Le * np.sin(eps), o[1] + Le * np.cos(eps)), C_DRIVE, 2.2)
    _arrow(ax, o, (o[0], o[1] + Le * np.cos(eps)), C_LIFT, 2.2)
    _arrow(ax, (o[0], o[1] + Le * np.cos(eps)), (o[0] + Le * np.sin(eps), o[1] + Le * np.cos(eps)), C_DRAG, 2.2)
    ax.text(o[0] - 0.42, o[1] + 1.2, "dL", color=C_LIFT, fontsize=10)
    ax.text(o[0] + 0.5, o[1] + 1.2, r"$dL_e$", color=C_DRIVE, fontsize=10)
    ax.text(o[0] + 0.08, o[1] + Le * np.cos(eps) + 0.1, r"$dD_i$", color=C_DRAG, fontsize=10)
    _arc(ax, o, 0.8, np.pi / 2 - eps, np.pi / 2, INK)
    ax.text(o[0] + 0.03, o[1] + 0.88, "ε", fontsize=9)
    ax.text(-0.6, -1.0, r"$\varepsilon=\tan^{-1}(w/U)\simeq w/U$;  $dL\simeq\rho U\Gamma\,dy$,  $dD_i\simeq\rho w\Gamma\,dy$", ha="center", fontsize=9)
    _bare(ax, (-3.5, 3.2), (-1.3, 2.6), "lift tilted back by the downwash")
    return fig


def draw_span_semicircle(ax=None):
    """The span angle γ: a point of the span is the foot of a point on a semicircle of radius s/2, y = −(s/2) cos γ."""
    fig, ax = _ax(ax, (5.6, 3.4))
    t = np.linspace(0, np.pi, 100)
    ax.plot(-np.cos(t), np.sin(t), color=C_REF, lw=1.2)
    ax.plot([-1, 1], [0, 0], color=INK, lw=4)
    for g in np.arange(1, 8) * np.pi / 8:
        ax.plot([-np.cos(g), -np.cos(g)], [0, np.sin(g)], color=C_REF, lw=0.7, ls=":")
        ax.plot([-np.cos(g)], [0], "o", color=C_CIRC, ms=4)
    g = 5 * np.pi / 8
    ax.plot([0, -np.cos(g)], [0, np.sin(g)], color=C_DRIVE, lw=1.6)
    ax.plot([-np.cos(g), -np.cos(g)], [0, np.sin(g)], color=C_DRIVE, lw=1.6)
    _arc(ax, (0, 0), 0.25, np.pi - g, np.pi, C_DRIVE)
    ax.text(-0.38, 0.1, "γ", color=C_DRIVE, fontsize=11)
    ax.text(-np.cos(g) + 0.04, 0.45, "sin γ", color=C_DRIVE, fontsize=8)
    ax.text(-1, -0.2, "γ = 0\nport tip, y = −s/2", ha="center", fontsize=7)
    ax.text(1, -0.2, "γ = π\nstarboard tip, y = +s/2", ha="center", fontsize=7)
    ax.text(0, -0.2, "γ = π/2\nmid-span", ha="center", fontsize=7)
    _bare(ax, (-1.5, 1.5), (-0.45, 1.25), r"span angle: $y=-(s/2)\cos\gamma$ (equal steps of γ crowd the tips)")
    return fig


def draw_stalled_section(ax=None, alpha_deg: float = 20.0):
    """A stalled section (a sketch, not a computation): the suction-side layer leaves near the nose and a wide wake of
    slow, recirculating fluid covers the upper surface."""
    fig, ax = _ax(ax, (6.8, 3.4))
    x, y = _section(e=0.1, beta=0.08)
    z = 0.5 * (x + 1j * y) * np.exp(-1j * np.deg2rad(alpha_deg))
    ax.fill(z.real, z.imag, color="#e9ecf3", ec=INK, zorder=3)
    k = int(np.argmin(z.real))
    nose, tail = z[k], z[0]
    xs = np.linspace(nose.real, 2.9, 60)
    top = nose.imag + 0.25 + 0.42 * np.tanh((xs - nose.real) * 1.4)
    bot = tail.imag + 0.0 * xs - 0.05
    ax.fill_between(xs[xs > nose.real + 0.05], np.maximum(bot, -0.6)[xs > nose.real + 0.05], top[xs > nose.real + 0.05], color=C_DRAG, alpha=0.13, zorder=1)
    ax.plot(xs, top, color=C_DRAG, lw=1.4, ls="--")
    for yy in (-1.0, -0.7, 1.2, 1.5):
        xx = np.linspace(-2.6, 2.9, 80)
        bump = 0.25 * np.exp(-((xx - 0.2) / 1.2) ** 2) * (1 if yy > 0 else -0.6)
        ax.plot(xx, yy + bump, color=C_FLOW, lw=1)
    for c0 in ((0.6, 0.25), (1.5, 0.15)):
        tt = np.linspace(0, 1.8 * np.pi, 40)
        ax.plot(c0[0] + 0.28 * np.cos(tt), c0[1] + 0.18 * np.sin(tt), color=C_DRAG, lw=1)
    ax.plot([nose.real + 0.12], [nose.imag + 0.2], "o", color=C_DRAG, ms=6, zorder=5)
    ax.text(nose.real + 0.0, nose.imag + 0.42, "separation", color=C_DRAG, fontsize=8, ha="center")
    ax.text(1.9, 0.75, "wide wake", color=C_DRAG, fontsize=8)
    _arrow(ax, (-2.6, -0.2), (-1.8, -0.2), C_FLOW, 2)
    ax.text(-2.6, -0.05, "U", color=C_FLOW, fontsize=10)
    _bare(ax, (-2.8, 3.0), (-1.35, 1.85), "stalled section — a sketch (lift falls, drag rises)")
    return fig


def draw_velocity_triangle(kind: str = "fish", ax=None, *, U: float | None = None, V: float | None = None):
    """Velocity triangle of a foil that moves across a stream, with the lift normal to the relative velocity and its
    thrust part.  kind: "fish" (tail moving sideways, top view), "bird_down" / "bird_up" (wing stroke, side view),
    "sail" (true wind, boat velocity, apparent wind).  ``U``, ``V`` (keywords, optional): forward speed and the size of the
    foil's cross-stream speed [m/s] for the three stroke sketches (defaults: the illustrative fish; 10 and 4 m/s for the bird).

    The drawn V arrow is the cross-stream part of the *relative flow* (opposite to the foil's own motion).  Fish and
    downstroke: ``ch14.oscillating_foil`` with V > 0, the lift leans forward (thrust).  Upstroke: the wing moves up, the
    relative wind comes from above-ahead (V < 0 in ``ch14.oscillating_foil``) and a positive lift leans *backward*; a bird
    feathers the wing there (nearly parallel to U_r), so the force is small."""
    fig, ax = _ax(ax, (5.6, 4.2))
    I = ch14.illustrative_inputs()
    if kind == "sail":
        aw = ch14.apparent_wind(I["U_wind"], I["V_boat"], I["theta"])
        sc = 1.0 / I["U_wind"]
        th = I["theta"]
        wind = -np.array([np.cos(th), np.sin(th)]) * I["U_wind"] * sc
        boat = np.array([I["V_boat"], 0.0]) * sc
        app = np.array(aw["vector"]) * sc
        o = np.array([1.2, 1.1])
        _arrow(ax, o, o + wind, C_FLOW, 2.2)
        _arrow(ax, o + wind - boat, o + wind, C_DRIVE, 2.2)
        _arrow(ax, o, o + app, COLORS["amber"], 2.2)
        ax.text(*(o + 0.5 * wind + np.array([0.0, 0.14])), "true wind U", color=C_FLOW, fontsize=9)
        ax.text(*(o + wind - 0.5 * boat + np.array([-0.1, -0.2])), "boat velocity V", color=C_DRIVE, fontsize=9, ha="center")
        ax.text(*(o + 0.5 * app + np.array([-0.95, 0.0])), r"apparent wind $U_r=U-V$", color=COLORS["amber"], fontsize=9)
        ax.text(-0.9, -0.55, f"θ = {np.rad2deg(th):.0f}° off the true wind → β_a = {np.rad2deg(aw['beta_a']):.1f}° off the bow, "
                             f"U_r = {aw['Ur']:.2f} m/s", fontsize=8)
        _bare(ax, (-1.2, 1.6), (-0.8, 1.5), "apparent wind: the wind the boat feels")
        return fig
    if kind == "fish":
        U0, V0, Ldraw, sgn = I["U_fish"], I["V_tail"], 1.0, 1.0
        names = ("U (water past the fish)", "V: water across the tail\n(the tail moves the other way)", "thrust", "side force")
        title = "fish tail (top view):\nthe lift normal to $U_r$ has a forward part"
        tail_note = "thrust power = sideways power (no drag)"
    elif kind == "bird_down":
        U0, V0, Ldraw, sgn = 10.0, 4.0, 1.0, 1.0
        names = ("U (flight speed)", "V: relative wind upward\n(the wing moves down)", "thrust", "upward force")
        title = "bird, downstroke (side view):\nthe lift leans forward — weight support and thrust"
        tail_note = "the forward part of the lift is thrust"
    elif kind == "bird_up":
        U0, V0, Ldraw, sgn = 10.0, 4.0, 0.45, -1.0
        names = ("U (flight speed)", "V: relative wind downward\n(the wing moves up)", "backward part", "upward force")
        title = "bird, upstroke (side view):\na positive lift leans backward — so the wing is feathered"
        tail_note = "small lift (drawn enlarged) leaning backward"
    else:
        raise ValueError('kind must be "fish", "bird_down", "bird_up" or "sail"')
    U = U0 if U is None else float(U)
    V = V0 if V is None else abs(float(V))
    r = ch14.oscillating_foil(U=U, V=sgn * V, L=1.0, D=0.0)     # V < 0: the upstroke (relative wind from above-ahead)
    th = float(r["theta"])                                     # signed: negative on the upstroke
    sc = 1.3 / max(U, V)
    o = np.array([0.0, 0.0 if sgn > 0 else 0.95])
    u = np.array([U, 0.0]) * sc
    v = np.array([0.0, sgn * V]) * sc
    _arrow(ax, o, o + u, C_FLOW, 2.2)
    _arrow(ax, o + u, o + u + v, C_FLOW, 1.6)
    _arrow(ax, o, o + u + v, C_DRIVE, 2.2)
    ax.text(o[0] + 0.5 * u[0], o[1] - sgn * 0.1, names[0], color=C_FLOW, fontsize=8, ha="center", va="top" if sgn > 0 else "bottom")
    ax.text(o[0] + u[0] + 0.06, o[1] + 0.5 * v[1], names[1], color=C_FLOW, fontsize=8, va="center")
    ax.text(o[0] + 0.4 * u[0] - 0.12, o[1] + 0.4 * v[1] + sgn * 0.17, r"$U_r$", color=C_DRIVE, fontsize=11, va="center")
    _arc(ax, o, 0.5, 0, th, INK)
    ax.text(o[0] + 0.56, o[1] + 0.56 * np.tan(0.5 * th), "θ", fontsize=10, va="center")
    f0 = np.array([-1.3 if sgn > 0 else -2.0, 0.2])
    Lv = 1.2 * Ldraw * np.array([-np.sin(th), np.cos(th)])       # normal to U_r: forward (−x) for θ > 0, backward (+x) for θ < 0
    part_col = COLORS["accent"] if sgn > 0 else C_DRAG
    _arrow(ax, f0, f0 + Lv, C_LIFT, 2.4)
    _arrow(ax, f0, f0 + np.array([Lv[0], 0.0]), part_col, 2.0)
    _arrow(ax, f0 + np.array([Lv[0], 0.0]), f0 + Lv, C_REF, 1.4)
    ax.text(*(f0 + Lv + np.array([0.0, 0.08])), "lift L", color=C_LIFT, fontsize=9, ha="center")
    ax.text(f0[0] + 0.5 * Lv[0], f0[1] - 0.08, f"{names[2]}\n= L sin {'θ' if sgn > 0 else '|θ|'}", color=part_col, fontsize=8, ha="center", va="top")
    ax.text(f0[0] + Lv[0] - sgn * 0.08, f0[1] + 0.5 * Lv[1], f"{names[3]}\n= L cos θ", color=C_REF, fontsize=7,
            ha="right" if sgn > 0 else "left", va="center")
    ax.text(0.0, -0.66, f"V/U = tan|θ| = {V / U:.2f} → |θ| = {np.rad2deg(abs(th)):.1f}°; {tail_note}", fontsize=8, ha="center")
    ax.text(-2.8, -0.42, "← forward (the body moves this way)", fontsize=7, color=C_REF)
    _bare(ax, (-2.9, 2.9), (-0.85, 1.75), title)
    return fig


def draw_sail_polygon(sf: dict | None = None, ax=None):
    """Force diagram of a sail on a course along +x: lift normal to the apparent wind, drag along it, their resultant
    split into the drive (along the course) and the lateral force.  ``sf``: the dict of ``ch14.sail_forces`` (default:
    our boat)."""
    fig, ax = _ax(ax, (5.6, 4.6))
    if sf is None:
        I = ch14.illustrative_inputs()
        sf = ch14.sail_forces(I["U_wind"], I["V_boat"], I["theta"], CL=I["CL_sail"], CD=I["CD_sail"], rho=I["rho"], A=I["A_sail"])
    ba = abs(float(sf["beta_a"]))
    sc = 1.0 / float(sf["F"])
    wd = -np.array([np.cos(ba), np.sin(ba)])
    Lv = float(sf["L"]) * sc * np.array([np.sin(ba), -np.cos(ba)])
    Dv = float(sf["D"]) * sc * wd
    o = np.array([0.0, 0.0])
    ax.plot([-0.5, 0.6], [0, 0], color=INK, lw=7, alpha=0.2, solid_capstyle="round")
    _arrow(ax, o - 1.3 * wd, o - 0.75 * wd, C_FLOW, 2.2)
    ax.text(*(o - 1.4 * wd + np.array([0.0, 0.08])), "apparent wind", color=C_FLOW, fontsize=8, ha="center")
    _arrow(ax, o, o + Lv, C_LIFT, 2.2)
    _arrow(ax, o + Lv, o + Lv + Dv, C_DRAG, 2.2)
    _arrow(ax, o, o + Lv + Dv, INK, 1.4)
    _arrow(ax, o, (float(sf["drive"]) * sc, 0.0), C_DRIVE, 2.8)
    _arrow(ax, (float(sf["drive"]) * sc, 0.0), (float(sf["drive"]) * sc, -float(sf["lateral"]) * sc), COLORS["amber"], 2.0)
    ax.text(*(0.86 * Lv + np.array([0.09, 0.0])), "lift (normal to the apparent wind)", color=C_LIFT, fontsize=8, va="center")
    ax.text(*(Lv + 0.5 * Dv + np.array([0.06, -0.13])), "drag (along the apparent wind)", color=C_DRAG, fontsize=8, ha="left")
    ax.text(0.5 * float(sf["drive"]) * sc, 0.07, "drive", color=C_DRIVE, fontsize=9, ha="center")
    ax.text(float(sf["drive"]) * sc + 0.05, -0.27 * float(sf["lateral"]) * sc, "lateral, to leeward\n(the keel resists it)", color=COLORS["amber"], fontsize=8, va="center")
    ax.text(*(0.6 * (Lv + Dv) + np.array([-0.05, 0.0])), "resultant F", color=INK, fontsize=8, ha="right", va="center")
    ax.text(0.68, -0.01, "bow →", color=INK, fontsize=8, va="center")
    ax.text(-1.42, 0.42, "windward side\n(the wind comes from here)", color=C_REF, fontsize=7)
    ax.text(-1.42, -0.5, "leeward side\n(the wind blows toward here)", color=C_REF, fontsize=7)
    _arc(ax, o, 0.5, 0, ba, INK)
    ax.text(0.55, 0.12, r"$\beta_a$", fontsize=10)
    ax.text(-0.1, -1.25, f"drive = L sin β_a − D cos β_a = {float(sf['drive']):.0f} N;  lateral = {float(sf['lateral']):.0f} N;  β_a = {np.rad2deg(ba):.1f}°",
            fontsize=8, ha="center")
    _bare(ax, (-1.5, 1.5), (-1.4, 1.2), "sail forces seen from above (course along +x)")
    return fig


# ======================================================================================================================
# composite figures
# ======================================================================================================================
def fig_four_circles():
    """Four circles and their images under z = ζ + b²/ζ (14.3): flat plate, circular arc, symmetric section, cambered
    airfoil (our e and β, b = 1)."""
    import matplotlib.pyplot as plt

    I = ch14.illustrative_inputs()
    e, beta, b = I["e"], 2.0 * I["beta"], 1.0
    cases = [("centred, a = b → flat plate", 0.0, 0.0), ("raised → circular arc", 0.0, beta), ("shifted back → symmetric section", e, 0.0),
             ("both → cambered airfoil", e, beta)]
    fig, ax = plt.subplots(2, 4, figsize=(15, 6.2))
    tt = np.linspace(0, 2 * np.pi, 300)
    for k, (name, ee, bb) in enumerate(cases):
        A = ch14.zhukhovsky_airfoil(b, ee, bb, 600)
        ax[0, k].plot(np.cos(tt), np.sin(tt), color=C_REF, lw=1, ls="--")
        ax[0, k].plot(np.append(A["zeta"].real, A["zeta"].real[0]), np.append(A["zeta"].imag, A["zeta"].imag[0]), color=C_FLOW)
        ax[0, k].plot([A["mu"].real], [A["mu"].imag], "o", color=C_CIRC, ms=5)
        ax[0, k].plot([1, -1], [0, 0], "x", color=C_DRAG, ms=7)
        ax[0, k].set(aspect="equal", xlim=(-1.6, 1.6), ylim=(-1.5, 1.6), xlabel="ξ / b", ylabel="η / b", title=name)
        ax[1, k].plot(np.append(A["x"], A["x"][0]), np.append(A["y"], A["y"][0]), color=C_LIFT)
        ax[1, k].set(aspect="equal", xlim=(-2.6, 2.6), ylim=(-1.2, 1.2), xlabel="x / b", ylabel="y / b",
                     title=f"e = {ee:g}, β = {np.rad2deg(bb):.0f}°")
    return fig


def fig_three_stages(U: float, b: float, e: float, beta: float, alpha: float):
    """Streamlines round the section for Γ/Γ_K = 0, ½, 1: the rear stagnation point moves to the trailing edge."""
    import matplotlib.pyplot as plt

    a = ch14.zhukhovsky_circle(b, e, beta)[1]
    GK = ch14.kutta_circulation(U, a, alpha, beta)
    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.2))
    for k, f in enumerate((0.0, 0.5, 1.0)):
        F = ch14.zhukhovsky_flow(U, b, e, beta, alpha, Gamma_cw=f * GK)
        streamlines(ax[k], F, (-3.2 * b, 3.6 * b), (-1.9 * b, 1.9 * b), levels=27)
        st = ch14.stagnation_points(U, b, e, beta, alpha, Gamma_cw=f * GK)
        ax[k].plot([st["z_front"].real, st["z_rear"].real], [st["z_front"].imag, st["z_rear"].imag], "o", color=C_DRAG, ms=6, zorder=6)
        ax[k].set(xlabel="x [m]", ylabel="y [m]", title=rf"$\Gamma/\Gamma_K$ = {f:g}: " + ("no lift, edge speed unbounded", "half way", "Kutta: leaves the edge smoothly")[k])
    return fig


def fig_pressure_two_ways(S: dict):
    """The pressure of ``ch14.zhukhovsky_surface`` drawn two ways: arrows normal to the surface (outward = suction) and
    −C_p against x/c with the loop shaded (its area is the normal-force coefficient)."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.4))
    z = S["z"]
    sc = 0.055 * S["chord_x"]
    ax[0].fill(S["x"], S["y"], color="#e9ecf3", ec=INK, lw=1.2)
    zc = np.concatenate([z[-1:], z, z[:1]])
    tang = zc[2:] - zc[:-2]
    nrm = -1j * tang / np.abs(tang)
    step = max(1, z.size // 60)
    for k in range(step // 2 + 1, z.size, step):
        cp = S["cp"][k]
        if not np.isfinite(cp):
            continue
        tip = z[k] + nrm[k] * sc * abs(cp)
        a0, a1 = (z[k], tip) if cp < 0 else (tip, z[k])
        ax[0].annotate("", (a1.real, a1.imag), (a0.real, a0.imag), arrowprops=dict(arrowstyle="->", color=C_PRES if cp < 0 else C_FLOW, lw=1))
    ax[0].set(aspect="equal", xlabel="x [m]", ylabel="y [m]", title="pressure normal to the surface (orange: suction, blue: push)")
    ax[0].margins(0.12, 0.6)
    up, lo = S["upper"], S["lower"]
    ax[1].plot(up["xc"], -up["cp"], color=C_PRES, label="upper surface")
    ax[1].plot(lo["xc"], -lo["cp"], color=C_FLOW, label="lower surface")
    ax[1].fill_between(up["xc"], -up["cp"], np.interp(up["xc"], lo["xc"], -lo["cp"]), color=C_LIFT, alpha=0.18)
    cy = ch14.normal_force_coefficient(up["xc"], up["cp"], lo["xc"], lo["cp"])
    ax[1].set(xlabel="x/c [–]", ylabel=r"$-C_p$ [–]", title=rf"loop area $C_y$ = {cy:.3f}; suction peak $C_p$ = {S['cp_min']:.2f}")
    ax[1].legend(fontsize=8)
    return fig


def fig_planforms_loadings():
    """Three planforms of the same span and area with the loading Γ(y) and the downwash w(y) the lifting-line solver
    gives each (our demo numbers: Λ = 6, K = 2π, α = 5°)."""
    import matplotlib.pyplot as plt

    I = ch14.illustrative_inputs()
    s, c, al, K = I["wd_s"], I["wd_c"], I["wd_alpha"], I["wd_K"]
    A = s * c
    plans = {"rectangular": ch14.planform("rectangular", s, c), "tapered 0.4": ch14.planform("tapered", s, 2 * A / (s * 1.4), 0.4),
             "elliptic": ch14.planform("elliptic", s, 4 * A / (np.pi * s))}
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
    gg = np.linspace(0.0, np.pi, 241)
    yy = ch14.span_position(gg, s)
    for (name, P), col in zip(plans.items(), (C_FLOW, C_CIRC, C_LIFT)):
        r = ch14.solve_lifting_line(P, s, I["wd_U"], al, K=K, N=40, symmetric=True, rho=I["wd_rho"])
        ax[0].plot(P["y"], 0.25 * P["c"], color=col, label=name)
        ax[0].plot(P["y"], -0.75 * P["c"], color=col)
        ax[1].plot(yy, ch14.circulation_series(r["Gn"], gg), color=col, label=rf"{name}: $e_{{span}}$ = {r['e_span']:.3f}")
        ax[2].plot(yy[2:-2], ch14.downwash_series(r["Gn"], gg[2:-2], s), color=col, label=rf"{name}: $C_L$ = {r['CL']:.3f}")
    ax[0].set(aspect="equal", xlabel="y [m]", ylabel="x [m]", title="planforms: same span and area (Λ = 6)")
    ax[1].set(xlabel="y [m]", ylabel="Γ(y) [m²/s]", title="bound circulation")
    ax[2].set(xlabel="y [m]", ylabel="downwash w(y) [m/s]", title="downwash: uniform only for the elliptic loading", ylim=(0, 0.12))
    for a in ax:
        a.legend(fontsize=8)
    return fig


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch14.illustrative_inputs()
    fig, ax = plt.subplots(3, 4, figsize=(22, 13))
    a = ax.ravel()
    draw_planform_labelled(ax=a[0])
    draw_three_axes(ax=a[1])
    draw_section_labelled(ax=a[2])
    draw_axes_conventions(ax=a[3])
    draw_kelvin_circuit(ax=a[4])
    draw_triangle_OQP(ax=a[5])
    draw_circle_angles(ax=a[6])
    draw_tip_flow_rear_view(ax=a[7])
    draw_horseshoe_staircase(7, ax=a[8])
    draw_force_triangle(ax=a[9])
    draw_span_semicircle(ax=a[10])
    draw_stalled_section(ax=a[11])
    save(fig, out, "drawings")
    fig, ax = plt.subplots(1, 5, figsize=(26, 4.8))
    for k, kind in enumerate(("fish", "bird_down", "bird_up", "sail")):
        draw_velocity_triangle(kind, ax=ax[k])
    draw_sail_polygon(ax=ax[4])
    save(fig, out, "drawings_triangles")
    save(fig_four_circles(), out, "drawings_four_circles")
    save(fig_three_stages(I["U_section"], I["b"], I["e"], I["beta"], I["alpha"]), out, "drawings_three_stages")
    S = ch14.zhukhovsky_surface(I["U_section"], I["b"], I["e"], I["beta"], I["alpha"], 1200)
    save(fig_pressure_two_ways(S), out, "drawings_pressure_two_ways")
    save(fig_planforms_loadings(), out, "drawings_planforms_loadings")
    # standalone calls with no Axes (what the notebook does) must return a Figure
    n_ok = sum(f().__class__.__name__ == "Figure" for f in (draw_planform_labelled, draw_three_axes, draw_section_labelled, draw_axes_conventions,
                                                              draw_kelvin_circuit, draw_triangle_OQP, draw_circle_angles, draw_tip_flow_rear_view,
                                                              draw_horseshoe_staircase, draw_force_triangle, draw_span_semicircle, draw_stalled_section,
                                                              draw_velocity_triangle, draw_sail_polygon))
    print(f"14 sketches drawn ({n_ok} return a Figure when called with no arguments) and 4 composite figures; the three stages use "
          f"Γ_K = {ch14.kutta_circulation(I['U_section'], ch14.zhukhovsky_circle(I['b'], I['e'], I['beta'])[1], I['alpha'], I['beta']):.4f} m²/s (clockwise)")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
