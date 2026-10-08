"""Schematic drawings for chapter 13 — drawn by this code, in our own words; no physics is computed here beyond a few
profile shapes taken from ``fluidpy``.

* ``draw_tangent_plane``        the local east–north–up axes on the rotating sphere and the two components of Ω
* ``draw_parcel_balance``       a parcel between pressure force and Coriolis force, moving along the isobars
* ``draw_thermal_wind_wedge``   tilted density surfaces and the wind that must change with height
* ``draw_taylor_column``        a low obstacle on the floor of a rotating tank and the column of fluid it carries
* ``draw_ekman_layer_sketch``   the turning current under a wind-stressed surface, or the turning wind over the ground
* ``draw_shallow_layer``        one layer of mean depth H with surface displacement η
* ``draw_mode_stack``           a stratified ocean as a stack of shallow layers of different equivalent depth
* ``draw_kelvin_sections``      a wave leaning on a coast: the surface under a crest and under a trough
* ``draw_pv_column``            a column of fluid stretched and squashed while keeping (ζ + f)/h
* ``draw_step_flow_sketch``     a zonal stream meeting a step in depth: eastward (waves) and westward (no waves)
* ``draw_eady_wedge``           sloping density surfaces and the wedge in which exchanging parcels release energy
* ``draw_cascade_arrows``       energy going to large scales and enstrophy to small scales from the injection scale

Every helper takes an optional matplotlib Axes (keyword ``ax``; ``None`` makes a new figure) and returns the Figure; the notebook
imports them with ``sys.path.insert(0, "scripts"); from ch13_drawings import *``.  Colours follow the chapter code:
pressure-gradient force = orange, Coriolis force = teal, friction = rose, velocity / tendency = purple (accent),
density and temperature = blue, vorticity = amber, ghosts = muted.  The hemisphere of every "to the right" is stated.

Run: ``.venv/Scripts/python.exe scripts/ch13_drawings.py --no-show``   Figure -> outputs/ch13/drawings.png.
"""
from __future__ import annotations

import numpy as np
from ch13_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch13_geophysical_fluid_dynamics as ch13

__all__ = ["draw_tangent_plane", "draw_parcel_balance", "draw_thermal_wind_wedge", "draw_taylor_column",
           "draw_ekman_layer_sketch", "draw_shallow_layer", "draw_mode_stack", "draw_kelvin_sections", "draw_pv_column",
           "draw_step_flow_sketch", "draw_eady_wedge", "draw_cascade_arrows"]


def _ax(ax, figsize=(5.2, 4.2)):
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    return ax


def _arrow(ax, tail, head, color, lw: float = 2.0, style: str = "->"):
    ax.annotate("", head, tail, arrowprops=dict(arrowstyle=style, color=color, lw=lw))


def _bare(ax, xlim, ylim, title: str):
    ax.set(xlim=xlim, ylim=ylim, xticks=[], yticks=[], title=title)
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
    return ax.figure


def draw_tangent_plane(ax=None, lat_deg: float = 35.0):
    """The local axes at latitude θ on the rotating Earth, and the rotation vector split into its local-vertical part
    Ω sin θ (which gives f) and its local-horizontal part Ω cos θ."""
    ax = _ax(ax)
    th = np.linspace(0, 2 * np.pi, 200)
    ax.plot(np.cos(th), np.sin(th), color=COLORS["muted"])
    ax.plot([-1, 1], [0, 0], color=COLORS["muted"], lw=0.8, ls=":")
    _arrow(ax, (0, -1.25), (0, 1.45), COLORS["ink"], 1.5)
    ax.text(0.05, 1.42, "Ω", fontsize=12)
    la = np.deg2rad(lat_deg)
    P = np.array([np.cos(la), np.sin(la)])
    up, north = P, np.array([-np.sin(la), np.cos(la)])
    ax.plot([0, P[0]], [0, P[1]], color=COLORS["muted"], lw=0.8)
    ax.text(0.33, 0.07, "θ", fontsize=11)
    _arrow(ax, P, P + 0.5 * up, COLORS["accent"])
    ax.text(*(P + 0.55 * up), "z (up)", color=COLORS["accent"], fontsize=9)
    _arrow(ax, P, P + 0.5 * north, COLORS["accent"])
    ax.text(*(P + 0.55 * north + [-0.25, 0.0]), "y (north)", color=COLORS["accent"], fontsize=9)
    ax.plot(*P, "o", color=COLORS["accent"])
    ax.text(P[0] + 0.05, P[1] - 0.12, "x (east) into the page", color=COLORS["accent"], fontsize=8)
    Q = P + np.array([0.75, -0.35])
    _arrow(ax, Q, Q + [0, 0.6], COLORS["ink"], 1.5)
    _arrow(ax, Q, Q + 0.6 * np.sin(la) * up, COLORS["teal"])
    _arrow(ax, Q, Q + 0.6 * np.cos(la) * north, COLORS["muted"])
    ax.text(*(Q + 0.6 * np.sin(la) * up + [0.03, 0.0]), "Ω sin θ = f/2", color=COLORS["teal"], fontsize=9)
    ax.text(*(Q + 0.6 * np.cos(la) * north + [-0.55, 0.06]), "Ω cos θ", color=COLORS["muted"], fontsize=9)
    ax.set_aspect("equal")
    return _bare(ax, (-1.4, 2.4), (-1.4, 1.7), "local axes on the rotating sphere")


def draw_parcel_balance(hemisphere: int = +1, ax=None):
    """A parcel in geostrophic balance: pressure force toward low pressure, Coriolis force opposite, motion along the
    isobars — low pressure on the left in the northern hemisphere (``hemisphere=+1``), on the right in the southern."""
    ax = _ax(ax)
    for y in np.linspace(-1, 1, 5):
        ax.plot([-2, 2], [y, y], color=COLORS["orange"], lw=0.8)
    ax.text(2.05, 0.0, "isobars", fontsize=8, color=COLORS["orange"], va="center")
    ax.text(-1.95, 1.12, "low pressure", color=COLORS["orange"], fontsize=9)
    ax.text(-1.95, -1.25, "high pressure", color=COLORS["orange"], fontsize=9)
    ax.plot(0, 0, "o", color=COLORS["accent"], ms=10)
    _arrow(ax, (0, 0), (0, 0.9), COLORS["orange"], 2.5)
    ax.text(0.08, 0.8, "pressure force", color=COLORS["orange"], fontsize=9)
    _arrow(ax, (0, 0), (0, -0.9), COLORS["teal"], 2.5)
    ax.text(0.08, -0.85, "Coriolis force", color=COLORS["teal"], fontsize=9)
    s = 1 if hemisphere >= 0 else -1
    _arrow(ax, (0, 0), (1.2 * s, 0), COLORS["accent"], 2.5)
    ax.text(0.6 * s, 0.1, "velocity", color=COLORS["accent"], fontsize=9, ha="center")
    name = "northern, f > 0" if s > 0 else "southern, f < 0"
    ax.text(-1.95, -1.45, "Coriolis force to the " + ("right" if s > 0 else "left") + " of the motion", fontsize=8, color=COLORS["teal"])
    return _bare(ax, (-2.1, 2.9), (-1.6, 1.4), f"geostrophic balance ({name})")


def draw_thermal_wind_wedge(ax=None):
    """Density surfaces sloping down toward the warm side, and the eastward wind that has to increase with height to
    stay in balance (northern hemisphere, cold air to the north)."""
    import matplotlib.patches as mp

    ax = _ax(ax)
    y = np.linspace(-2, 2, 50)
    for z0 in np.linspace(0.2, 1.8, 6):
        ax.plot(y, z0 + 0.25 * y, color=COLORS["blue"], lw=1.0)
    ax.text(-1.95, 2.2, "warm (light)", color=COLORS["blue"], fontsize=9)
    ax.text(1.0, -0.2, "cold (dense)", color=COLORS["blue"], fontsize=9)
    for z0, L in ((0.3, 0.2), (0.9, 0.5), (1.5, 0.8), (2.1, 1.1)):
        ax.add_patch(mp.Circle((0.0, z0), 0.04 + 0.08 * L, color=COLORS["accent"], alpha=0.8))
    ax.text(0.25, 2.05, "wind out of the page (eastward),\nstronger aloft", color=COLORS["accent"], fontsize=9)
    _arrow(ax, (-2.2, -0.2), (-2.2, 2.4), COLORS["ink"], 1.0)
    ax.text(-2.15, 2.4, "z", fontsize=9)
    _arrow(ax, (-2.2, -0.2), (2.3, -0.2), COLORS["ink"], 1.0)
    ax.text(2.2, -0.42, "y (north)", fontsize=9, ha="right")
    return _bare(ax, (-2.4, 2.6), (-0.6, 2.7), "thermal wind: density gradient ⇒ shear")


def draw_taylor_column(ax=None):
    """A short obstacle towed along the floor of a rapidly rotating tank: the fluid above it moves with it as a column,
    and the rest flows round the column at every height."""
    import matplotlib.patches as mp

    ax = _ax(ax)
    ax.plot([-2, 2], [0, 0], color=COLORS["ink"], lw=2)
    ax.plot([-2, 2], [2, 2], color=COLORS["blue"], lw=1.2)
    ax.add_patch(mp.Rectangle((-0.4, 0), 0.8, 0.5, color=COLORS["muted"]))
    ax.add_patch(mp.Rectangle((-0.4, 0.5), 0.8, 1.5, color=COLORS["amber"], alpha=0.25))
    ax.plot([-0.4, -0.4], [0.5, 2], color=COLORS["amber"], ls="--")
    ax.plot([0.4, 0.4], [0.5, 2], color=COLORS["amber"], ls="--")
    ax.text(0, 1.3, "column moves\nwith the obstacle", ha="center", fontsize=8, color=COLORS["amber"])
    for z in (0.3, 0.9, 1.5):
        _arrow(ax, (-1.9, z), (-0.7, z), COLORS["accent"], 1.2)
    ax.text(-1.9, 1.7, "flow goes round, at every height", fontsize=8, color=COLORS["accent"])
    _arrow(ax, (1.5, 0.2), (1.5, 1.8), COLORS["ink"], 1.2)
    ax.text(1.55, 1.75, "Ω", fontsize=12)
    return _bare(ax, (-2.1, 2.1), (-0.2, 2.3), "Taylor column")


def draw_ekman_layer_sketch(kind: str = "surface", hemisphere: int = +1, ax=None):
    """The turning of the flow in an Ekman layer, seen from above: arrows at increasing distance from the boundary.
    ``kind="surface"``: the current under a wind stress (northern hemisphere: 45° to the right at the surface, turning
    further right with depth).  ``kind="bottom"``: the wind over the ground under a geostrophic flow (45° to the left
    at the ground in the northern hemisphere, i.e. toward low pressure)."""
    ax = _ax(ax)
    lat = hemisphere * np.deg2rad(60.0)
    f = ch13.coriolis_parameter(lat)
    if kind == "surface":
        d = ch13.ekman_depth(0.03, f)
        zs = -d * np.array([0.0, 0.4, 0.8, 1.2, 1.8, 2.6])
        u, v = ch13.ekman_surface(zs, 0.07, 0.0, 1027.0, 0.03, f)
        sc = 1.0 / np.hypot(u[0], v[0])
        _arrow(ax, (0, 0), (1.3, 0), COLORS["orange"], 3)
        ax.text(1.3, 0.05, "wind stress", color=COLORS["orange"], fontsize=9)
        lab, ttl = "depth", "Ekman current under the wind"
    elif kind == "bottom":
        d = ch13.ekman_depth(7.0, f)
        zs = d * np.array([0.15, 0.4, 0.8, 1.2, 1.8, 3.0])
        u, v = ch13.ekman_bottom(zs, 12.0, 0.0, 7.0, f)
        sc = 1.0 / 12.0
        _arrow(ax, (0, 0), (1.0, 0), COLORS["muted"], 3)
        ax.text(1.0, -0.1, "geostrophic wind", color=COLORS["muted"], fontsize=9)
        ax.text(-0.9, 0.75 * np.sign(f), "low pressure", color=COLORS["orange"], fontsize=9)
        lab, ttl = "height", "Ekman wind above the ground"
    else:
        raise ValueError('kind must be "surface" or "bottom"')
    for j, (uu, vv) in enumerate(zip(u * sc, v * sc)):
        _arrow(ax, (0, 0), (uu, vv), COLORS["accent"], 2.2 - 0.25 * j)
        ax.text(uu, vv, f" {abs(zs[j]) / d:.1f}δ", fontsize=7, color=COLORS["muted"])
    ax.set_aspect("equal")
    hemi = "N" if hemisphere >= 0 else "S"
    ax.text(-1.05, -1.15, f"labels: {lab} in units of δ", fontsize=7, color=COLORS["muted"])
    return _bare(ax, (-1.1, 1.9), (-1.2, 1.2), f"{ttl} ({hemi} hemisphere)")


def draw_shallow_layer(ax=None):
    """One homogeneous layer: mean depth H, surface displacement η, velocity independent of depth."""
    ax = _ax(ax)
    x = np.linspace(0, 4, 200)
    eta = 0.25 * np.sin(1.6 * x)
    ax.fill_between(x, 0, 2 + eta, color=COLORS["blue"], alpha=0.15)
    ax.plot(x, 2 + eta, color=COLORS["blue"])
    ax.plot([0, 4], [2, 2], color=COLORS["muted"], ls=":")
    ax.plot([0, 4], [0, 0], color=COLORS["ink"], lw=2)
    _arrow(ax, (3.6, 0), (3.6, 2), COLORS["ink"], 1.0, "<->")
    ax.text(3.65, 1.0, "H", fontsize=11)
    _arrow(ax, (1.0, 2.0), (1.0, 2 + 0.25 * np.sin(1.6)), COLORS["blue"], 1.0, "<->")
    ax.text(1.05, 2.1, "η", fontsize=11, color=COLORS["blue"])
    for z in (0.4, 1.0, 1.6):
        _arrow(ax, (1.8, z), (2.5, z), COLORS["accent"], 1.5)
    ax.text(1.8, 1.75, "u, v the same at every depth", fontsize=8, color=COLORS["accent"])
    return _bare(ax, (-0.1, 4.2), (-0.2, 2.6), "shallow-water layer")


def draw_mode_stack(ax=None):
    """A continuously stratified ocean replaced, mode by mode, by shallow layers: the barotropic mode keeps the full
    depth, each baroclinic mode behaves like a layer of much smaller equivalent depth (bars not to scale)."""
    ax = _ax(ax)
    I = ch13.illustrative_inputs()
    m = ch13.modes_uniform_N(I["ocean_N"], I["ocean_H"], n_modes=4, nz=101)
    z = m.z / I["ocean_H"]
    for n in range(4):
        ax.plot(n * 1.6 + 0.5 * m.psi[n], z, color=COLORS["accent"])
        ax.plot([n * 1.6, n * 1.6], [-1, 0], color=COLORS["muted"], lw=0.6)
        ax.text(n * 1.6, 0.08, f"n = {n}", ha="center", fontsize=9)
        ax.text(n * 1.6, -1.18, f"H_e = {m.He[n]:.3g} m\nc = {m.c[n]:.3g} m/s", ha="center", fontsize=7.5, color=COLORS["teal"])
    ax.text(-0.75, -0.5, "ψ_n(z)", rotation=90, fontsize=9, color=COLORS["accent"], va="center")
    return _bare(ax, (-0.9, 5.7), (-1.45, 0.25), "vertical modes: a stack of shallow layers")


def draw_kelvin_sections(hemisphere: int = +1, ax=None):
    """Looking along a coast in the direction a Kelvin wave travels: under a crest the surface slopes up toward the
    coast, under a trough down; the coast is on the right of the direction of travel in the northern hemisphere."""
    ax = _ax(ax)
    y = np.linspace(0, 3, 100)
    s = 1 if hemisphere >= 0 else -1
    yy = 3 - y if s > 0 else y
    ax.plot(yy, 1 + 0.5 * np.exp(-y), color=COLORS["accent"], label="under a crest")
    ax.plot(yy, 1 - 0.5 * np.exp(-y), color=COLORS["rose"], label="under a trough")
    ax.plot([0, 3], [1, 1], color=COLORS["muted"], ls=":")
    xc = 3.0 if s > 0 else 0.0
    ax.plot([xc, xc], [0, 1.7], color=COLORS["ink"], lw=4)
    ax.text(xc, 1.75, "coast", ha="center", fontsize=9)
    _arrow(ax, (xc - s * 1.0, 0.3), (xc, 0.3), COLORS["muted"], 1.0, "<->")
    ax.text(xc - s * 0.5, 0.36, "Λ = c/|f|", ha="center", fontsize=9)
    ax.text(1.5, 1.82, "wave travels into the page", ha="center", fontsize=8)
    ax.legend(fontsize=8, loc="lower center")
    hemi = "northern: coast on the right" if s > 0 else "southern: coast on the left"
    return _bare(ax, (-0.3, 3.3), (-0.05, 1.95), f"Kelvin wave ({hemi})")


def draw_pv_column(ax=None):
    """A column of fluid keeps (ζ + f)/h: stretched it spins up cyclonically, squashed anticyclonically (northern
    hemisphere)."""
    ax = _ax(ax)
    import matplotlib.patches as mp

    for x0, h, w, lab, spin in ((0.0, 1.0, 1.0, "h₀, ζ = 0", 0), (2.2, 1.5, 0.82, "stretched:\nζ > 0", 1),
                                (4.4, 0.6, 1.3, "squashed:\nζ < 0", -1)):
        ax.add_patch(mp.Rectangle((x0 - w / 2, 0), w, h, color=COLORS["blue"], alpha=0.2))
        ax.add_patch(mp.Ellipse((x0, h), w, 0.18, color=COLORS["blue"], alpha=0.5))
        ax.text(x0, -0.35, lab, ha="center", fontsize=8)
        if spin:
            t = np.linspace(0.3, 5.2, 40) * spin
            ax.plot(x0 + 0.3 * w * np.cos(t), h + 0.25 + 0.06 * np.sin(t), color=COLORS["amber"], lw=2)
            ax.annotate("", (x0 + 0.3 * w * np.cos(t[-1]), h + 0.25 + 0.06 * np.sin(t[-1])),
                        (x0 + 0.3 * w * np.cos(t[-3]), h + 0.25 + 0.06 * np.sin(t[-3])),
                        arrowprops=dict(arrowstyle="->", color=COLORS["amber"], lw=2))
    ax.text(2.2, 2.1, "(ζ + f)/h is the same for all three", ha="center", fontsize=9, color=COLORS["amber"])
    return _bare(ax, (-0.9, 5.4), (-0.7, 2.3), "potential vorticity of a column")


def draw_step_flow_sketch(ax=None):
    """A zonal stream crossing a step in depth, seen from above: an eastward stream is left with a standing wave
    downstream, a westward stream bends once and already upstream of the step (curves from the linearised solution of
    ``ch13.flow_over_step`` with our inputs)."""
    ax = _ax(ax, (6.0, 4.2))
    I = ch13.illustrative_inputs()
    f, b = ch13.coriolis_parameter(I["lat"]), ch13.beta_parameter(I["lat"])
    x = np.linspace(-6e6, 1.1e7, 400)
    e = ch13.flow_over_step(x, I["U_mean"], b, f, 4200.0, 4000.0)
    w = ch13.flow_over_step(x, -I["U_mean"], b, f, 4200.0, 4000.0)
    sc = 1.0 / abs(e["Y"]).max()
    for off in (0.0, 0.7):
        ax.plot(x / 1e6, 2.2 + off + sc * e["Y"], color=COLORS["accent"])
        ax.plot(x / 1e6, -1.2 + off + 2 * sc * w["Y"], color=COLORS["rose"])
    ax.axvline(0, color=COLORS["ink"], lw=1.5, ls="--")
    ax.text(0.15, 3.6, "step", fontsize=9)
    _arrow(ax, (-5.5, 3.3), (-3.5, 3.3), COLORS["accent"])
    ax.text(-5.5, 3.45, "eastward: standing Rossby wave", color=COLORS["accent"], fontsize=8)
    _arrow(ax, (10.5, 0.1), (8.5, 0.1), COLORS["rose"])
    ax.text(4.4, 0.25, "westward: no wave, felt upstream", color=COLORS["rose"], fontsize=8)
    return _bare(ax, (-6.2, 11.2), (-2.6, 3.9), "flow over a step on a β-plane (north is up)")


def draw_eady_wedge(ax=None):
    """Sloping density surfaces of a baroclinic current and the wedge between them and the horizontal: parcels
    exchanged along a path inside the wedge move light fluid up and dense fluid down, and release potential energy."""
    ax = _ax(ax)
    y = np.linspace(-2, 2, 50)
    slope = 0.3
    for z0 in np.linspace(0.3, 1.7, 5):
        ax.plot(y, z0 + slope * y, color=COLORS["blue"], lw=1.0)
    ax.fill([0, 1.6, 1.6], [1.0, 1.0, 1.0 + slope * 1.6], color=COLORS["amber"], alpha=0.3)
    ax.text(1.65, 1.2, "wedge of\ninstability", fontsize=8, color=COLORS["amber"])
    _arrow(ax, (0.15, 1.02), (1.4, 1.0 + 0.5 * slope * 1.4), COLORS["accent"], 2)
    _arrow(ax, (-0.15, 0.98), (-1.4, 1.0 - 0.5 * slope * 1.4), COLORS["accent"], 2)
    ax.text(-1.9, 0.35, "dense parcel\nsinks toward the light side", fontsize=7.5, color=COLORS["accent"])
    ax.text(0.2, 0.55, "light parcel rises\ntoward the dense side", fontsize=7.5, color=COLORS["accent"])
    ax.text(-1.95, 2.35, "light", color=COLORS["blue"], fontsize=9)
    ax.text(1.5, -0.1, "dense", color=COLORS["blue"], fontsize=9)
    return _bare(ax, (-2.2, 2.6), (-0.4, 2.7), "baroclinic instability: the wedge")


def draw_cascade_arrows(ax=None):
    """The two inertial ranges of two-dimensional turbulence on logarithmic axes: energy handed to larger scales along
    a −5/3 slope, enstrophy to smaller scales along a −3 slope, from the scale at which the flow is stirred."""
    ax = _ax(ax)
    k1, k2 = np.array([1.0, 10.0]), np.array([10.0, 100.0])
    ax.loglog(k1, k1 ** (-5.0 / 3.0), color=COLORS["accent"], lw=2.5)
    ax.loglog(k2, 10.0 ** (-5.0 / 3.0) * (k2 / 10.0) ** -3.0, color=COLORS["amber"], lw=2.5)
    ax.axvline(10.0, color=COLORS["muted"], ls="--", lw=0.8)
    ax.text(10.5, 0.5, "stirring\nscale K₀", fontsize=8, color=COLORS["muted"])
    ax.annotate("", (1.6, 0.012), (7.0, 0.012), arrowprops=dict(arrowstyle="->", color=COLORS["accent"], lw=2))
    ax.text(1.7, 0.016, "energy → large scales", color=COLORS["accent"], fontsize=8)
    ax.annotate("", (80.0, 0.012), (14.0, 0.012), arrowprops=dict(arrowstyle="->", color=COLORS["amber"], lw=2))
    ax.text(14.5, 0.016, "enstrophy → small scales", color=COLORS["amber"], fontsize=8)
    ax.text(2.0, 0.45, r"$K^{-5/3}$", color=COLORS["accent"], fontsize=11)
    ax.text(32.0, 0.0012, r"$K^{-3}$", color=COLORS["amber"], fontsize=11)
    ax.set(xlabel="wavenumber K (log)", ylabel="energy spectrum (log)", xticks=[], yticks=[], title="two cascades in two dimensions")
    ax.minorticks_off()
    ax.grid(False)
    return ax.figure


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(4, 4, figsize=(22, 18))
    a = ax.ravel()
    draw_tangent_plane(a[0])
    draw_parcel_balance(+1, ax=a[1])
    draw_parcel_balance(-1, ax=a[2])
    draw_thermal_wind_wedge(a[3])
    draw_taylor_column(a[4])
    draw_ekman_layer_sketch("surface", +1, ax=a[5])
    draw_ekman_layer_sketch("surface", -1, ax=a[6])
    draw_ekman_layer_sketch("bottom", +1, ax=a[7])
    draw_shallow_layer(a[8])
    draw_mode_stack(a[9])
    draw_kelvin_sections(+1, ax=a[10])
    draw_kelvin_sections(-1, ax=a[11])
    draw_pv_column(a[12])
    draw_step_flow_sketch(a[13])
    draw_eady_wedge(a[14])
    draw_cascade_arrows(a[15])
    save(fig, out, "drawings")
    solo = [fn() for fn in (draw_tangent_plane, draw_cascade_arrows)]
    print(f"drew {len(__all__)} sketches: {', '.join(__all__)}; called without an Axes they return a {type(solo[0]).__name__}")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
