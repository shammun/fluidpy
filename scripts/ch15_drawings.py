"""Schematic drawings for chapter 15 — drawn by this code, in our own words; no physics is computed here beyond calls
to the tested functions of ``fluidpy.ch15_compressible_flow`` (the helpers only draw what those functions return).

Sketches (each takes an optional matplotlib Axes, keyword ``ax``; ``None`` makes a new figure; each returns the Figure):

* ``draw_duct_control_volume``    a slowly varying duct with stations 1 and 2, the control volume, wall force and heat
* ``draw_reservoir_nozzle_tank``  a reservoir (p₀, T₀) feeding a nozzle that discharges into a tank at the back pressure
* ``draw_four_ducts``             converging / diverging × subsonic / supersonic: which is a nozzle, which a diffuser
* ``draw_shock_two_frames``       a normal shock at rest and the same shock running into still gas (arrows to scale)
* ``draw_shock_triangles``        the velocity triangles of an oblique shock: (u₁, v) and (u₂, v), angles σ and δ
* ``draw_wedge_streamlines``      a wedge with its attached shock and streamlines turning through it
* ``draw_mach_triangle``          the right triangle behind sin μ = 1/M
* ``draw_plate_waves``            a flat plate at incidence with its two shocks and two fans
* ``draw_bow_shock_sketch``       a detached (bow) shock ahead of a blunt wedge — a sketch, not a computation
* ``draw_jet_sketch``             the jet outside a nozzle exit: over-expanded, design or under-expanded — a sketch
* ``draw_fan``                    a Prandtl–Meyer fan from ``ch15.expansion_fan`` with a streamline
* ``draw_wavefront_panels``       wave fronts of a moving source, one panel per ``AC.source_wavefronts`` result

The notebook imports them with ``sys.path.insert(0, "scripts"); from ch15_drawings import *``.  Colours follow the
chapter code: subsonic blue, supersonic rose, sonic amber, pressure orange, temperature teal, entropy / shocks accent,
references muted.  Station 1 is upstream; σ and δ are the shock angle and the deflection (β and θ elsewhere).

Run: ``.venv/Scripts/python.exe scripts/ch15_drawings.py --no-show``   Figures -> outputs/ch15/drawings*.png
"""
from __future__ import annotations

import numpy as np
from ch15_common import C_ENT, C_PRES, C_REF, C_SONIC, C_SUB, C_SUP, C_TEMP, COLORS, finish, parse_args, save, setup

from fluidpy import ch15_compressible_flow as ch15

__all__ = ["draw_duct_control_volume", "draw_reservoir_nozzle_tank", "draw_four_ducts", "draw_shock_two_frames",
           "draw_shock_triangles", "draw_wedge_streamlines", "draw_mach_triangle", "draw_plate_waves",
           "draw_bow_shock_sketch", "draw_jet_sketch", "draw_fan", "draw_wavefront_panels"]

INK = COLORS["ink"]
WALL = "#d5d9e4"


def _ax(ax, figsize=(5.8, 4.2)):
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    return ax.figure, ax


def _arrow(ax, tail, head, color, lw: float = 2.0, style: str = "->"):
    ax.annotate("", head, tail, arrowprops=dict(arrowstyle=style, color=color, lw=lw, shrinkA=0, shrinkB=0))


def _bare(ax, xlim, ylim, title: str, aspect="equal"):
    ax.set(xlim=xlim, ylim=ylim, aspect=aspect, title=title)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(False)


def _arc(ax, c, r, a0, a1, color, lw=1.2):
    t = np.linspace(a0, a1, 40)
    ax.plot(c[0] + r * np.cos(t), c[1] + r * np.sin(t), color=color, lw=lw)


# ======================================================================================================================
# one-dimensional flow
# ======================================================================================================================
def draw_duct_control_volume(ax=None):
    """A duct whose area changes slowly, the control volume between stations 1 (upstream) and 2, the wall force F on
    the fluid and the heat Q added per unit mass — the picture behind (15.18), (15.22), (15.24)."""
    fig, ax = _ax(ax, (7.0, 3.6))
    x = np.linspace(0.0, 4.0, 200)
    h = 0.9 - 0.35 * np.sin(np.pi * x / 4.0) ** 2 + 0.08 * x  # half-height (our shape)
    ax.fill_between(x, h, h + 0.12, color=WALL)
    ax.fill_between(x, -h, -h - 0.12, color=WALL)
    ax.plot(x, h, color=INK, lw=1.5)
    ax.plot(x, -h, color=INK, lw=1.5)
    x1, x2 = 0.9, 3.1
    for xs, lab in ((x1, "1"), (x2, "2")):
        hs = float(np.interp(xs, x, h))
        ax.plot([xs, xs], [-hs, hs], color=C_ENT, lw=1.6, ls="--")
        ax.text(xs, -hs - 0.36, f"station {lab}", ha="center", fontsize=9, color=C_ENT)
        for yy in np.linspace(-0.6 * hs, 0.6 * hs, 3):
            _arrow(ax, (xs - 0.3, yy), (xs + 0.3, yy), C_SUB, 1.6)
    m = (x >= x1) & (x <= x2)
    ax.plot(x[m], h[m] - 0.03, color=C_ENT, lw=1.2, ls="--")
    ax.plot(x[m], -h[m] + 0.03, color=C_ENT, lw=1.2, ls="--")
    ax.text(x1 - 0.42, 0.0, r"$\rho_1,u_1,p_1,A_1$", ha="right", va="center", fontsize=9)
    ax.text(x2 + 0.42, 0.0, r"$\rho_2,u_2,p_2,A_2$", ha="left", va="center", fontsize=9)
    _arrow(ax, (2.0, 1.45), (2.0, 0.85), C_PRES, 1.8)
    ax.text(2.08, 1.3, "Q  (heat added per unit mass)", fontsize=9, color=C_PRES)
    _arrow(ax, (1.6, -0.62), (2.3, -0.62), C_TEMP, 1.8)
    ax.text(1.95, -0.5, "F  (walls on the fluid, +x)", ha="center", fontsize=9, color=C_TEMP)
    _arrow(ax, (3.6, -1.75), (4.2, -1.75), C_REF, 1.2)
    ax.text(4.25, -1.75, "x", va="center", fontsize=9, color=C_REF)
    _bare(ax, (-1.2, 5.4), (-1.95, 1.7), "control volume in a slowly varying duct")
    return fig


def draw_reservoir_nozzle_tank(ax=None):
    """A large reservoir at the stagnation state (p₀, T₀), a convergent nozzle, and a tank held at the back pressure
    p_B — the setting of §15.7."""
    fig, ax = _ax(ax, (7.0, 3.4))
    ax.fill([0, 2.2, 2.2, 0], [-1.3, -1.3, 1.3, 1.3], color="#eef1f8")
    ax.plot([0, 2.2, 2.2], [1.3, 1.3, 0.55], color=INK, lw=1.6)
    ax.plot([0, 2.2, 2.2], [-1.3, -1.3, -0.55], color=INK, lw=1.6)
    x = np.linspace(2.2, 3.6, 60)
    h = 0.55 - 0.3 * np.sin(0.5 * np.pi * (x - 2.2) / 1.4)
    ax.plot(x, h, color=INK, lw=1.6)
    ax.plot(x, -h, color=INK, lw=1.6)
    ax.plot([3.6, 3.6, 6.2], [h[-1], 1.3, 1.3], color=INK, lw=1.6)
    ax.plot([3.6, 3.6, 6.2], [-h[-1], -1.3, -1.3], color=INK, lw=1.6)
    ax.text(1.0, 0.25, "reservoir", ha="center", fontsize=10)
    ax.text(1.0, -0.3, r"$p_0,\ T_0,\ u\approx 0$", ha="center", fontsize=10, color=C_PRES)
    ax.text(4.9, 0.75, "tank", ha="center", fontsize=10)
    ax.text(4.9, -0.85, r"back pressure $p_B$", ha="center", fontsize=10, color=C_PRES)
    for yy in (-0.12, 0.12):
        _arrow(ax, (2.5, yy), (3.5, yy), C_SUB, 1.6)
    _arrow(ax, (3.7, 0.0), (4.7, 0.0), C_SUP, 2.0)
    ax.plot([3.6, 3.6], [-h[-1], h[-1]], color=C_SONIC, lw=1.4, ls="--")
    ax.text(3.6, -0.6, "exit plane", ha="center", fontsize=8, color=C_SONIC)
    _bare(ax, (-0.2, 6.4), (-1.5, 1.5), "reservoir → nozzle → tank: lower $p_B$ and watch the exit")
    return fig


def draw_four_ducts(ax=None):
    """The four cases of the area rule du/u = −(dA/A)/(1 − M²) (15.34): converging or diverging, subsonic or
    supersonic.  The role of each duct is taken from ``ch15.duct_role``."""
    fig, ax = _ax(ax, (8.0, 4.6))
    cases = [(0.5, -1, 0.0, 1.5), (0.5, +1, 4.4, 1.5), (2.0, -1, 0.0, -1.5), (2.0, +1, 4.4, -1.5)]
    for M, sgn, x0, y0 in cases:
        x = np.linspace(0.0, 3.0, 2)
        h = (0.75 - 0.4 * x / 3.0) if sgn < 0 else (0.35 + 0.4 * x / 3.0)
        col = C_SUB if M < 1 else C_SUP
        ax.fill_between(x0 + x, y0 + h, y0 + h + 0.08, color=WALL)
        ax.fill_between(x0 + x, y0 - h, y0 - h - 0.08, color=WALL)
        ax.plot(x0 + x, y0 + h, color=INK, lw=1.5)
        ax.plot(x0 + x, y0 - h, color=INK, lw=1.5)
        role = ch15.duct_role(M, sgn)
        l_in, l_out = (0.45, 1.0) if role == "nozzle" else (1.0, 0.45)
        _arrow(ax, (x0 + 0.15, y0), (x0 + 0.15 + l_in, y0), col, 2.0)
        _arrow(ax, (x0 + 2.85 - l_out, y0), (x0 + 2.85, y0), col, 2.0)
        ax.text(x0 + 1.5, y0 - 1.1, f"{'M < 1' if M < 1 else 'M > 1'}, {'converging' if sgn < 0 else 'diverging'}: {role}",
                ha="center", fontsize=9, color=col)
    _bare(ax, (-0.4, 7.8), (-2.9, 2.6), "area rule (15.34): the same shape does opposite jobs below and above M = 1")
    return fig


def draw_shock_two_frames(ms: dict | None = None, ax=None):
    """A normal shock seen at rest (gas enters at u₁, leaves at u₂) and the same shock running into still gas (front at
    u₁, gas behind following at u₁ − u₂).  ``ms`` is the dict returned by ``ch15.moving_shock`` (arrows to scale);
    default: our M_s of ``illustrative_inputs``."""
    fig, ax = _ax(ax, (8.4, 3.4))
    if ms is None:
        I = ch15.illustrative_inputs()
        ms = ch15.moving_shock(I["p_air"], I["T_air"], M_s=I["M_s"])
    u1, ug = float(ms["u_shock"]), float(ms["u_gas"])
    u2 = u1 - ug
    s = 1.5 / u1  # arrow length per m/s
    for x0, lab in ((0.0, "frame of the shock"), (5.0, "frame of the still gas")):
        ax.plot([x0 + 2.0, x0 + 2.0], [-1.0, 1.0], color=C_ENT, lw=3)
        ax.text(x0 + 2.0, -1.5, lab, ha="center", fontsize=9)
    _arrow(ax, (0.3, 0.0), (0.3 + s * u1, 0.0), C_SUP, 2.2)
    ax.text(0.3 + 0.5 * s * u1, 0.18, rf"$u_1$ = {u1:.0f} m/s", ha="center", fontsize=8, color=C_SUP)
    _arrow(ax, (2.3, 0.0), (2.3 + s * u2, 0.0), C_SUB, 2.2)
    ax.text(2.3 + 0.5 * s * u2 + 0.35, 0.18, rf"$u_2$ = {u2:.0f} m/s", ha="center", fontsize=8, color=C_SUB)
    ax.text(0.9, 0.75, "ahead (1)", ha="center", fontsize=8, color=C_REF)
    ax.text(3.0, 0.75, "behind (2)", ha="center", fontsize=8, color=C_REF)
    # moving frame: the shock runs to the left into still gas (subtract u₁ everywhere)
    _arrow(ax, (7.0, 1.2), (7.0 - s * u1, 1.2), C_ENT, 2.2)
    ax.text(7.0 - 0.5 * s * u1, 1.38, rf"front: {u1:.0f} m/s", ha="center", fontsize=8, color=C_ENT)
    ax.text(5.9, 0.0, "u = 0", ha="center", fontsize=9, color=C_REF)
    _arrow(ax, (7.3 + s * ug, 0.0), (7.3, 0.0), C_SUB, 2.2)
    ax.text(7.3 + 0.5 * s * ug + 0.2, 0.18, rf"$u_1-u_2$ = {ug:.0f} m/s", ha="center", fontsize=8, color=C_SUB)
    _bare(ax, (-0.2, 9.6), (-1.8, 1.8), "one shock, two frames: subtract $u_1$ everywhere")
    return fig


# ======================================================================================================================
# planar waves
# ======================================================================================================================
def draw_mach_triangle(M: float = 2.0, ax=None):
    """The right triangle behind sin μ = 1/M (15.49): while the source travels U Δt, its first front has grown to
    c Δt; the tangent from the source to that circle makes the Mach angle μ with the path."""
    fig, ax = _ax(ax, (6.4, 3.6))
    mu = float(ch15.mach_angle(M))
    L = 3.0  # U Δt (drawn length)
    r = L / M
    t = np.linspace(0, 2 * np.pi, 200)
    ax.plot(r * np.cos(t), r * np.sin(t), color=C_SUB, lw=1.2)
    ax.plot([0], [0], "o", color=C_REF)
    ax.plot([L], [0], "o", color=C_SUP)
    tp = (r * np.sin(mu), r * np.cos(mu))  # tangent point
    ax.plot([0, L], [0, 0], color=INK, lw=1.4)
    ax.plot([0, tp[0]], [0, tp[1]], color=C_SUB, lw=1.4)
    ax.plot([L, tp[0]], [0, tp[1]], color=C_ENT, lw=1.8)
    ax.plot([L, tp[0]], [0, -tp[1]], color=C_ENT, lw=1.8)
    _arc(ax, (L, 0), 0.7, np.pi - mu, np.pi, C_ENT)
    ax.text(L - 0.95, 0.12, r"$\mu$", fontsize=11, color=C_ENT)
    ax.text(0.5 * L, -0.25, r"$U\,\Delta t$", ha="center", fontsize=10)
    ax.text(0.5 * tp[0] - 0.42, 0.5 * tp[1], r"$c\,\Delta t$", fontsize=10, color=C_SUB)
    ax.text(L + 0.1, 0.25, "source now", fontsize=8, color=C_SUP)
    ax.text(-0.1, -0.3, "source then", fontsize=8, ha="right", color=C_REF)
    _bare(ax, (-2.3, 4.4), (-1.9, 1.9), rf"$\sin\mu = c\,\Delta t/U\,\Delta t = 1/M$:  M = {M:g}, μ = {np.degrees(mu):.1f}°")
    return fig


def draw_wavefront_panels(fronts_list=None, ax=None):
    """Wave fronts of a moving point source, one panel per element of ``fronts_list`` (each the dict returned by
    ``AC.source_wavefronts``); circles labelled with their emission times, Mach lines drawn when M > 1.
    ``ax``: an array of Axes of the same length (``None`` makes a row of panels)."""
    import matplotlib.pyplot as plt

    if fronts_list is None:
        c = 340.0
        fronts_list = [ch15.source_wavefronts(U, c, 4.0, 4, 1.0) for U in (0.0, 0.5 * c, 2.0 * c)]
    n = len(fronts_list)
    if ax is None:
        fig, ax = plt.subplots(1, n, figsize=(5.2 * n, 4.4))
    axes = np.atleast_1d(ax).ravel()
    fig = axes[0].figure
    t = np.linspace(0, 2 * np.pi, 241)
    for a_, w in zip(axes, fronts_list):
        R = float(np.max(w["radii"]))
        for xc, r, te in zip(w["x_centres"], w["radii"], w["t_emit"]):
            a_.plot(xc + r * np.cos(t), r * np.sin(t), color=C_SUB, lw=1.1)
            a_.plot([xc], [0], ".", color=C_REF, ms=5)
            a_.text(xc, r, f"t = {te:g} s", fontsize=7, ha="center", va="bottom", color=C_SUB)
        a_.plot([w["x_source"]], [0], "o", color=C_SUP, ms=7)
        M = float(w["M"])
        if M > 1.0:
            Lc = 1.15 * (w["x_source"] - float(np.min(w["x_centres"])))
            for sg in (+1, -1):
                a_.plot([w["x_source"], w["x_source"] - Lc * np.cos(w["mu"])], [0, sg * Lc * np.sin(w["mu"])], color=C_ENT, lw=1.6)
            ttl = f"M = {M:g}: fronts pile up on the Mach cone, μ = {np.degrees(w['mu']):.1f}°"
        elif M > 0.0:
            ttl = f"M = {M:g}: fronts crowd ahead, but still lead the source"
        else:
            ttl = "source at rest: concentric fronts"
        lo = min(float(np.min(w["x_centres"] - w["radii"])), w["x_source"]) - 0.15 * R
        hi = max(float(np.max(w["x_centres"] + w["radii"])), w["x_source"]) + 0.15 * R
        a_.set(aspect="equal", xlim=(lo, hi), ylim=(-1.25 * R, 1.25 * R), xlabel="x [m]", ylabel="y [m]", title=ttl)
        a_.title.set_fontsize(9)
    return fig


def draw_shock_triangles(M1: float = 2.0, sigma: float | None = None, delta: float | None = None, ax=None):
    """The oblique shock in stream-aligned axes with its velocity triangles: V₁ = (u₁, v) ahead and V₂ = (u₂, v) behind,
    the same shock-parallel v on both sides; σ = shock angle, δ = deflection (β, θ elsewhere).  Give ``sigma`` (the
    deflection is then computed) or ``delta`` (weak root); lengths are to scale in units of V₁."""
    fig, ax = _ax(ax, (6.8, 4.6))
    if sigma is None:
        sigma = ch15.oblique_shock(M1, delta=np.deg2rad(10.0) if delta is None else delta)["sigma"]
    cp = ch15.oblique_shock_components(M1, sigma)
    d = float(cp["delta"]) if delta is None else float(delta)
    es = np.array([np.cos(sigma), np.sin(sigma)])  # along the shock
    en = np.array([np.sin(sigma), -np.cos(sigma)])  # normal to the shock, pointing downstream
    ax.plot([-1.3 * es[0], 1.5 * es[0]], [-1.3 * es[1], 1.5 * es[1]], color=C_ENT, lw=3)
    O = np.zeros(2)
    A = np.array([-1.0, 0.0])
    _arrow(ax, A, O, C_SUP, 2.4)
    ax.text(-0.55, 0.07, r"$V_1$", fontsize=11, color=C_SUP)
    V2 = float(cp["V2_over_V1"]) * np.array([np.cos(d), np.sin(d)])
    _arrow(ax, O, V2, C_SUB, 2.4)
    ax.text(V2[0] * 0.6, V2[1] * 0.6 + 0.1, r"$V_2$", fontsize=11, color=C_SUB)
    # components: V₁ = u₁ e_n + v e_s ; V₂ = u₂ e_n + v e_s
    u1, v, u2 = float(cp["u1_over_V1"]), float(cp["v_over_V1"]), float(cp["u2_over_V1"])
    _arrow(ax, A, A + u1 * en, C_PRES, 1.4)
    _arrow(ax, A + u1 * en, O, C_TEMP, 1.4)
    ax.text(*(A + 0.5 * u1 * en + np.array([-0.14, -0.08])), r"$u_1$", fontsize=10, color=C_PRES)
    ax.text(*(A + u1 * en + 0.5 * v * es + np.array([0.06, -0.1])), r"$v$", fontsize=10, color=C_TEMP)
    _arrow(ax, O, u2 * en, C_PRES, 1.4)
    _arrow(ax, u2 * en, V2, C_TEMP, 1.4)
    ax.text(*(0.5 * u2 * en + np.array([-0.16, -0.1])), r"$u_2$", fontsize=10, color=C_PRES)
    ax.text(*(u2 * en + 0.5 * v * es + np.array([0.06, -0.1])), r"$v$", fontsize=10, color=C_TEMP)
    ax.plot([0, 1.1], [0, 0], color=C_REF, lw=0.8, ls=":")
    _arc(ax, O, 0.42, 0.0, sigma, C_ENT)
    ax.text(0.45 * np.cos(0.6 * sigma), 0.45 * np.sin(0.6 * sigma) + 0.03, r"$\sigma$", fontsize=11, color=C_ENT)
    _arc(ax, O, 0.75, 0.0, d, C_SUB)
    ax.text(0.8, 0.02 + 0.4 * np.sin(d), r"$\delta$", fontsize=11, color=C_SUB)
    _bare(ax, (-1.35, 1.45), (-1.05, 1.25),
          rf"oblique shock, $M_1$ = {M1:g}: σ = {np.degrees(sigma):.1f}°, δ = {np.degrees(d):.1f}°;  $v$ unchanged, $u_2 < u_1$")
    return fig


def draw_wedge_streamlines(M1: float = 2.0, half_angle: float = np.deg2rad(10.0), sigma: float | None = None, ax=None):
    """A symmetric wedge in a supersonic stream with its attached weak shock (angle ``sigma`` from ``ch15.wedge_flow``)
    and streamlines that stay straight and turn by ``half_angle`` at the shock."""
    fig, ax = _ax(ax, (6.6, 4.4))
    wf = ch15.wedge_flow(M1, half_angle)
    if sigma is None:
        sigma = wf["sigma"]
    if not np.isfinite(sigma):
        ax.text(0.5, 0.5, "detached: no attached shock at this angle\n(see draw_bow_shock_sketch)", ha="center", va="center", transform=ax.transAxes)
        _bare(ax, (-1, 2), (-1.2, 1.2), f"wedge of {np.degrees(half_angle):.0f}° at M₁ = {M1:g}: detached")
        return fig
    L = 2.2
    ax.fill([0, L, L], [0, L * np.tan(half_angle), -L * np.tan(half_angle)], color=WALL)
    for sg in (+1, -1):
        ax.plot([0, L], [0, sg * L * np.tan(half_angle)], color=INK, lw=1.6)
        ax.plot([0, L], [0, sg * L * np.tan(sigma)], color=C_ENT, lw=2.4)
        for y0 in np.linspace(0.12, 1.15, 5):
            xs = y0 / np.tan(sigma)  # where this streamline meets the shock
            if xs >= L:
                ax.plot([-1.0, L], [sg * y0, sg * y0], color=C_SUP, lw=1.0)
                continue
            ax.plot([-1.0, xs], [sg * y0, sg * y0], color=C_SUP, lw=1.0)
            ax.plot([xs, L], [sg * y0, sg * (y0 + (L - xs) * np.tan(half_angle))], color=C_SUB if wf["subsonic_behind"] else C_SUP, lw=1.0)
    _arrow(ax, (-0.95, 0.0), (-0.35, 0.0), C_SUP, 2.0)
    ax.text(-0.95, 0.08, rf"$M_1$ = {M1:g}", fontsize=9, color=C_SUP)
    ax.text(1.25, 0.62 * 1.25 * np.tan(sigma) + 0.02, rf"$M_2$ = {wf['M2']:.2f}", fontsize=9)
    _bare(ax, (-1.05, L), (-1.3, 1.3),
          rf"wedge δ = {np.degrees(half_angle):.0f}°: shock at σ = {np.degrees(sigma):.1f}°, $p_2/p_1$ = {wf['p_ratio']:.2f} (δ_max = {np.degrees(wf['delta_max']):.1f}°)")
    return fig


def draw_bow_shock_sketch(ax=None):
    """A detached (bow) shock standing ahead of a wedge too blunt for an attached shock: normal on the axis, strong
    near it (subsonic behind), weak far out — **a sketch**: the shape and the stand-off distance are not computed."""
    fig, ax = _ax(ax, (6.2, 4.4))
    y = np.linspace(-1.5, 1.5, 200)
    xs = -0.55 + 0.42 * y ** 2  # our curve, for illustration only
    ax.plot(xs, y, color=C_ENT, lw=2.6)
    ax.fill([0, 0.9, 2.2, 2.2, 0.9], [0, 0.55, 0.55, -0.55, -0.55], color=WALL)
    ax.plot([2.2, 0.9, 0, 0.9, 2.2], [0.55, 0.55, 0, -0.55, -0.55], color=INK, lw=1.6)
    ys = np.linspace(-0.62, 0.62, 60)
    ax.plot(0.02 + 0.9 * np.abs(ys) ** 0.7, np.where(ys >= 0, 0.62 + 0.25 * (0.62 - np.abs(ys)), -0.62 - 0.25 * (0.62 - np.abs(ys))),
            color=C_SONIC, lw=1.2, ls="--")
    _arrow(ax, (-1.45, 0.0), (-0.75, 0.0), C_SUP, 2.0)
    ax.text(-1.45, 0.1, "M > 1", fontsize=9, color=C_SUP)
    ax.text(-0.42, 0.12, "M < 1", fontsize=9, color=C_SUB)
    ax.text(0.95, 1.05, "M > 1", fontsize=9, color=C_SUP)
    ax.text(-0.5, -0.5, "strong,\nnearly normal", fontsize=8, ha="right", color=C_ENT)
    ax.text(0.55, -1.42, "weak far out", fontsize=8, color=C_ENT)
    ax.text(0.5, 0.82, "M = 1", fontsize=8, color=C_SONIC)
    _bare(ax, (-1.6, 2.3), (-1.65, 1.65), "detached bow shock (δ > δ_max) — sketch, not computed")
    return fig


def draw_fan(fan: dict | None = None, ax=None):
    """A centred Prandtl–Meyer fan from the dict returned by ``ch15.expansion_fan``: the wall and its corner, the rays
    (first at μ₁ to the upstream wall, last at μ₂ to the downstream wall) and one streamline."""
    fig, ax = _ax(ax, (6.8, 4.2))
    if fan is None:
        I = ch15.illustrative_inputs()
        fan = ch15.expansion_fan(I["M_fan"], I["fan_turn"])
    w2 = float(fan["wall_angle2"])
    R = 1.9 * float(np.max(fan["streamline_r"]))
    ax.fill([-R, 0, R * np.cos(w2), R * np.cos(w2), -R], [0, 0, R * np.sin(w2), R * np.sin(w2) - 0.12 * R, -0.12 * R], color=WALL)
    ax.plot([-R, 0, R * np.cos(w2)], [0, 0, R * np.sin(w2)], color=INK, lw=1.8)
    for ang in fan["ray_angle"]:
        ax.plot([0, R * np.cos(ang)], [0, R * np.sin(ang)], color=C_SUB, lw=0.9)
    sx, sy = np.asarray(fan["streamline_x"]), np.asarray(fan["streamline_y"])
    ax.plot(np.concatenate([[-R], sx, [sx[-1] + 0.5 * R * np.cos(w2)]]), np.concatenate([[sy[0]], sy, [sy[-1] + 0.5 * R * np.sin(w2)]]),
            color=C_SUP, lw=1.6)
    M1, M2 = float(fan["M"][0]), float(fan["M2"])
    ax.text(-0.9 * R, sy[0] + 0.06 * R, rf"$M_1$ = {M1:.2f}", fontsize=9, color=C_SUP)
    ax.text(sx[-1] + 0.12 * R, sy[-1] + 0.1 * R, rf"$M_2$ = {M2:.2f}", fontsize=9, color=C_SUP)
    ax.text(0.62 * R * np.cos(fan["ray_angle"][0]), 0.62 * R * np.sin(fan["ray_angle"][0]) + 0.04 * R, r"first ray, $\mu_1$", fontsize=8, color=C_SUB)
    ax.text(0.8 * R * np.cos(fan["ray_angle"][-1]) + 0.02 * R, 0.8 * R * np.sin(fan["ray_angle"][-1]) - 0.09 * R, r"last ray, $\mu_2$", fontsize=8, color=C_SUB)
    _bare(ax, (-R, 1.35 * R), (min(R * np.sin(w2), 0) - 0.15 * R, 0.9 * R),
          rf"Prandtl–Meyer fan: turn {np.degrees(-w2):.0f}°, $p_2/p_1$ = {float(fan['p2_over_p1']):.3f}")
    return fig


def draw_plate_waves(pattern: dict | None = None, ax=None):
    """A flat plate at incidence in a supersonic stream with its four waves, from the dict returned by
    ``ch15.flat_plate_wave_pattern``: thick lines are shocks, thin fans are expansions."""
    fig, ax = _ax(ax, (7.0, 4.6))
    if pattern is None:
        I = ch15.illustrative_inputs()
        pattern = ch15.flat_plate_wave_pattern(I["alpha"], I["M_inf"])
    pa = float(pattern["plate_angle"])
    te = np.array([np.cos(pa), np.sin(pa)])
    ax.plot([0, te[0]], [0, te[1]], color=INK, lw=3.5)
    L = 0.95
    for a_ in np.linspace(pattern["le_fan_first"], pattern["le_fan_last"], 4):
        ax.plot([0, L * np.cos(a_)], [0, L * np.sin(a_)], color=C_SUB, lw=0.9)
    ax.plot([0, L * np.cos(pattern["le_shock"])], [0, L * np.sin(pattern["le_shock"])], color=C_ENT, lw=2.4)
    ax.plot([te[0], te[0] + L * np.cos(pattern["te_shock"])], [te[1], te[1] + L * np.sin(pattern["te_shock"])], color=C_ENT, lw=2.4)
    for a_ in np.linspace(pattern["te_fan_first"], pattern["te_fan_last"], 4):
        ax.plot([te[0], te[0] + L * np.cos(a_)], [te[1], te[1] + L * np.sin(a_)], color=C_SUB, lw=0.9)
    _arrow(ax, (-0.7, 0.32), (-0.2, 0.32), C_SUP, 2.0)
    ax.text(-0.7, 0.38, r"$M_\infty,\ p_\infty$", fontsize=9, color=C_SUP)
    mid = 0.5 * te
    ax.text(mid[0] + 0.02, mid[1] + 0.16, rf"$p_2/p_\infty$ = {float(pattern['p2']):.3f}  (expansion)", fontsize=8, color=C_SUB)
    ax.text(mid[0] - 0.2, mid[1] - 0.24, rf"$p_3/p_\infty$ = {float(pattern['p3']):.3f}  (compression)", fontsize=8, color=C_ENT)
    ax.text(mid[0] - 0.05, mid[1] + 0.03, "b", fontsize=9)
    _bare(ax, (-0.8, 2.0), (-0.95, 0.85), rf"flat plate at α = {np.degrees(-pa):.0f}°: shock below and fan above at the nose, the opposite at the tail")
    return fig


def draw_jet_sketch(kind: str = "underexpanded", ax=None):
    """The jet just outside a nozzle exit for ``kind`` = "overexpanded" (exit pressure below the back pressure: oblique
    shocks turn the jet inward), "design" (matched: a parallel jet) or "underexpanded" (exit pressure above the back
    pressure: expansion fans turn it outward) — **a sketch**: the wave angles are illustrative, not computed."""
    if kind not in ("overexpanded", "design", "underexpanded"):
        raise ValueError('kind must be "overexpanded", "design" or "underexpanded"')
    fig, ax = _ax(ax, (6.2, 3.6))
    x = np.linspace(-2.0, 0.0, 60)
    h = 0.35 + 0.25 * ((x + 2.0) / 2.0) ** 1.3
    ax.fill_between(x, h, h + 0.1, color=WALL)
    ax.fill_between(x, -h, -h - 0.1, color=WALL)
    ax.plot(x, h, color=INK, lw=1.6)
    ax.plot(x, -h, color=INK, lw=1.6)
    he = h[-1]
    turn = {"overexpanded": -0.22, "design": 0.0, "underexpanded": 0.22}[kind]
    for sg in (+1, -1):
        ax.plot([0, 1.6], [sg * he, sg * (he + 1.6 * np.tan(turn))], color=C_REF, lw=1.4, ls="--")  # jet boundary
        if kind == "overexpanded":
            ax.plot([0, 0.9], [sg * he, sg * (he - 0.9 * np.tan(0.75))], color=C_ENT, lw=2.2)
        elif kind == "underexpanded":
            for a_ in np.linspace(0.45, 0.85, 4):
                ax.plot([0, 0.9 * np.cos(a_)], [sg * he, sg * (he - 0.9 * np.sin(a_))], color=C_SUB, lw=0.9)
    _arrow(ax, (-1.6, 0.0), (-0.5, 0.0), C_SUP, 2.0)
    txt = {"overexpanded": r"$p_e < p_B$: oblique shocks, jet contracts", "design": r"$p_e = p_B$: no waves, parallel jet",
           "underexpanded": r"$p_e > p_B$: expansion fans, jet spreads"}[kind]
    ax.text(1.75, 0.0, r"$p_B$", fontsize=10, color=C_PRES, va="center")
    ax.text(-0.3, 0.0, r"$p_e$", fontsize=10, color=C_PRES, va="center")
    _bare(ax, (-2.1, 2.1), (-1.25, 1.25), f"{kind} jet — sketch: {txt}")
    return fig


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch15.illustrative_inputs()
    fig, ax = plt.subplots(3, 4, figsize=(26, 14))
    a = ax.ravel()
    draw_duct_control_volume(ax=a[0])
    draw_reservoir_nozzle_tank(ax=a[1])
    draw_four_ducts(ax=a[2])
    draw_shock_two_frames(ax=a[3])
    draw_mach_triangle(I["M_flight"], ax=a[4])
    wf = ch15.wedge_flow(I["M_wedge"], I["wedge_half_angle"])
    draw_shock_triangles(I["M_wedge"], wf["sigma"], I["wedge_half_angle"], ax=a[5])
    draw_wedge_streamlines(I["M_wedge"], I["wedge_half_angle"], wf["sigma"], ax=a[6])
    draw_bow_shock_sketch(ax=a[7])
    draw_fan(ch15.expansion_fan(I["M_fan"], I["fan_turn"]), ax=a[8])
    draw_plate_waves(ch15.flat_plate_wave_pattern(I["alpha"], I["M_inf"]), ax=a[9])
    draw_jet_sketch("overexpanded", ax=a[10])
    draw_jet_sketch("underexpanded", ax=a[11])
    save(fig, out, "drawings")
    c0 = ch15.PerfectGas(I["gamma"], I["R"]).c(I["T_air"])
    save(draw_wavefront_panels([ch15.source_wavefronts(U, c0, 4.0, 4, 1.0) for U in (0.0, 0.5 * c0, 2.0 * c0)]), out, "drawings_wavefronts")
    # standalone calls with no Axes (what the notebook does) must return a Figure
    fns = (draw_duct_control_volume, draw_reservoir_nozzle_tank, draw_four_ducts, draw_shock_two_frames, draw_shock_triangles,
           draw_wedge_streamlines, draw_mach_triangle, draw_plate_waves, draw_bow_shock_sketch, draw_jet_sketch, draw_fan,
           draw_wavefront_panels)
    n_ok = sum(f().__class__.__name__ == "Figure" for f in fns)
    print(f"{len(fns)} sketches drawn ({n_ok} return a Figure when called with no arguments); wedge M₁ = {I['M_wedge']}, δ = {np.degrees(I['wedge_half_angle']):.0f}°: "
          f"σ = {np.degrees(wf['sigma']):.2f}°; detached case handled: {draw_wedge_streamlines(1.3, np.deg2rad(20.0)).__class__.__name__}")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
