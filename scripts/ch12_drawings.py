"""Schematic drawings for chapter 12 (our analogues of the book's sketches — drawn by this code, never copied; no physics
is computed here beyond the profile shapes taken from ``fluidpy``):

* ``scales_sketch``          the outer scale L, the velocity difference ΔU and eddies down to η (cf. Fig. 12.11)
* ``parcel_sketch``          a parcel displaced across a mean shear arrives with u < 0, v > 0 (cf. Fig. 12.6)
* ``stress_element_sketch``  viscous and Reynolds shear stress on a fluid element (cf. Fig. 12.7)
* ``fg_geometry_sketch``     longitudinal f(r) and transverse g(r) velocity pairs (cf. Fig. 12.9)
* ``free_shear_sketch``      jet, wake or shear layer with self-similar profiles (cf. Fig. 12.13)
* ``wall_layers_sketch``     sublayer, buffer layer, logarithmic layer and wake region (cf. Figs. 12.16, 12.18)
* ``surface_layer_sketch``   forced convection below |L_M|, free convection above (cf. Fig. 12.21)
* ``plume_sketch``           time-averaged smoke plume: a wedge near the source, a parabola far away (cf. Fig. 12.27)

Every helper takes a matplotlib Axes first and returns it; the notebook imports them with
``sys.path.insert(0, "scripts"); from ch12_drawings import *``.  Colours follow the chapter code: mean = purple (accent),
fluctuation = teal, Reynolds stress = orange, viscous = rose, buoyancy/temperature = blue, scalar = amber, ghosts = muted.
Constants used for the shapes are **labelled illustrative** (ξ½ = 0.10; κ = 0.41, B = 5.0, "a common textbook pair").

Run: ``.venv/Scripts/python.exe scripts/ch12_drawings.py --no-show``   Figure -> outputs/ch12/drawings.png.
"""
from __future__ import annotations

import numpy as np
from ch12_common import COLORS, ILLUSTRATIVE_JET, KAPPA, B_LOG, finish, parse_args, save, setup

from fluidpy import ch12_turbulence as ch12

__all__ = ["scales_sketch", "parcel_sketch", "stress_element_sketch", "fg_geometry_sketch", "free_shear_sketch",
           "wall_layers_sketch", "surface_layer_sketch", "plume_sketch"]


def _arrow(ax, tail, head, color, lw: float = 2.0, style: str = "->"):
    ax.annotate("", head, tail, arrowprops=dict(arrowstyle=style, color=color, lw=lw))


def _bare(ax, xlim, ylim, title: str):
    ax.set(xlim=xlim, ylim=ylim, xticks=[], yticks=[], title=title)
    ax.grid(False)
    return ax


def scales_sketch(ax, n_tiers: int = 5, seed: int = 0):
    """A shear layer of thickness L with velocity difference ΔU and eddies of every size down to η (cf. Fig. 12.11)."""
    import matplotlib.patches as mp

    rng = np.random.default_rng(seed)
    y = np.linspace(-1.0, 1.0, 200)
    ax.plot(0.6 * np.tanh(2.5 * y) + 0.9, y, color=COLORS["accent"], lw=2.5)
    for yy in (-0.8, -0.4, 0.4, 0.8):
        _arrow(ax, (0.9, yy), (0.6 * np.tanh(2.5 * yy) + 0.9, yy), COLORS["accent"], 1.2)
    ax.plot([0.9, 0.9], [-1, 1], color=COLORS["muted"], lw=0.8)
    ax.text(0.25, 1.05, r"$U(y)$", color=COLORS["accent"])
    _arrow(ax, (0.3, -1.12), (1.5, -1.12), COLORS["accent"], 1.2, "<->")
    ax.text(0.9, -1.3, r"$\Delta U$", ha="center", color=COLORS["accent"])
    _arrow(ax, (1.75, -0.7), (1.75, 0.7), COLORS["ink"], 1.2, "<->")
    ax.text(1.82, 0.0, "L", va="center")
    x0, size = 2.9, 0.7
    for k in range(n_tiers):
        m = 2 ** k if k < 3 else 6
        for _ in range(m):
            cx = x0 + rng.uniform(-0.25, 0.25) * (k > 0)
            cy = rng.uniform(-0.85 + size, 0.85 - size) if size < 0.8 else 0.0
            ax.add_patch(mp.Circle((cx, cy), size, fill=False, ec=COLORS["teal"], lw=max(2.2 - 0.4 * k, 0.6)))
        ax.text(x0, -1.3, ("L" if k == 0 else (r"$\eta$" if k == n_tiers - 1 else r"$l'$")), ha="center", color=COLORS["teal"])
        x0 += size + 0.5 * size + 0.35
        size *= 0.5
    _arrow(ax, (2.9, 1.1), (x0 - 0.3, 1.1), COLORS["orange"], 1.5)
    ax.text(0.5 * (2.9 + x0 - 0.3), 1.18, r"energy passed down the scales at the rate $\bar\varepsilon\sim(\Delta U)^3/L$",
            ha="center", fontsize=9, color=COLORS["orange"])
    ax.text(x0 - 0.25, 0.25, "viscosity\nacts here", fontsize=8, color=COLORS["rose"], ha="center")
    ax.set_aspect("equal")
    return _bare(ax, (0.1, x0 + 0.3), (-1.45, 1.45), "outer scale L, eddies down to the Kolmogorov scale (cf. Fig. 12.11)")


def parcel_sketch(ax, shear: float = 1.0):
    """Parcels exchanged across a mean shear: the one moving up arrives slow (u < 0, v > 0), the other fast (cf. Fig. 12.6)."""
    import matplotlib.patches as mp

    y = np.linspace(0.0, 2.0, 50)
    U = 0.4 + 0.8 * shear * (y - 1.0) + 1.0
    ax.plot(U, y, color=COLORS["accent"], lw=2.5)
    ax.plot([0.2, 0.2], [0, 2], color=COLORS["muted"], lw=0.8)
    for yy in (0.25, 0.75, 1.25, 1.75):
        _arrow(ax, (0.2, yy), (0.4 + 0.8 * shear * (yy - 1.0) + 1.0, yy), COLORS["accent"], 1.0)
    ax.text(2.25, 1.9, r"$U(y)$", color=COLORS["accent"])
    ax.axhline(1.0, color=COLORS["muted"], ls="--", lw=0.8)
    for x0, y_from, col, lab in ((3.3, 0.45, COLORS["rose"], r"from below: $v>0$, $u\approx-\ell\,dU/dy<0$"),
                                 (4.6, 1.55, COLORS["blue"], r"from above: $v<0$, $u>0$")):
        ax.add_patch(mp.Circle((x0, y_from), 0.09, color=col, alpha=0.35))
        ax.add_patch(mp.Circle((x0, 1.0), 0.09, color=col))
        _arrow(ax, (x0, y_from), (x0, 1.0 - 0.1 * np.sign(1.0 - y_from)), col, 1.6)
        du = -0.8 * shear * (1.0 - y_from)
        _arrow(ax, (x0, 1.0), (x0 + du, 1.0), COLORS["teal"], 2.0)
        ax.text(x0 + 0.12, y_from, r"$\ell$", color=col, va="center")
        ax.text(x0 - 0.9, 2.12 if y_from > 1 else -0.2, lab, fontsize=8.5, color=col)
    ax.text(5.35, 1.03, r"either way $uv<0$" if shear > 0 else r"either way $uv>0$", fontsize=9, color=COLORS["orange"])
    ax.set_xlabel("x  (mean flow →)")
    ax.set_ylabel("y")
    return _bare(ax, (0.0, 6.6), (-0.3, 2.3), r"why $\overline{uv}<0$ when $dU/dy>0$ (cf. Fig. 12.6)")


def stress_element_sketch(ax):
    """Shear stresses on a fluid element in a mean shear: viscous μ dU/dy and Reynolds −ρ₀ mean(uv) (cf. Fig. 12.7)."""
    import matplotlib.patches as mp

    ax.add_patch(mp.Rectangle((1.0, 1.0), 2.0, 2.0, fill=True, fc=COLORS["grid"], ec=COLORS["ink"], lw=1.5))
    _arrow(ax, (1.2, 3.12), (2.8, 3.12), COLORS["rose"], 2.2)
    _arrow(ax, (2.8, 0.88), (1.2, 0.88), COLORS["rose"], 2.2)
    _arrow(ax, (1.2, 3.34), (2.8, 3.34), COLORS["orange"], 3.0)
    _arrow(ax, (2.8, 0.66), (1.2, 0.66), COLORS["orange"], 3.0)
    ax.text(2.95, 3.08, r"$\mu\,dU/dy$ (viscous)", color=COLORS["rose"], va="center", fontsize=9)
    ax.text(2.95, 3.36, r"$-\rho_0\overline{uv}$ (Reynolds)", color=COLORS["orange"], va="center", fontsize=9)
    ax.text(2.0, 2.0, r"$\bar\tau_{12}=\mu\frac{dU}{dy}-\rho_0\overline{uv}$", ha="center", va="center", fontsize=11)
    _arrow(ax, (0.3, 0.3), (0.9, 0.3), COLORS["ink"], 1.0)
    _arrow(ax, (0.3, 0.3), (0.3, 0.9), COLORS["ink"], 1.0)
    ax.text(0.95, 0.25, "x")
    ax.text(0.22, 0.95, "y")
    ax.set_aspect("equal")
    return _bare(ax, (0.0, 5.4), (0.1, 3.8), "two shear stresses on one element (cf. Fig. 12.7)")


def fg_geometry_sketch(ax):
    """Two points a distance r apart: velocity components along r give f(r), components across r give g(r) (cf. Fig. 12.9)."""
    for y0, lab, col, along in ((1.6, r"longitudinal:  $f(r)=\overline{u_1(\mathbf{x})\,u_1(\mathbf{x}+r\mathbf{e}_1)}/\overline{u_1^2}$",
                                 COLORS["accent"], True),
                                (0.5, r"transverse:  $g(r)=\overline{u_2(\mathbf{x})\,u_2(\mathbf{x}+r\mathbf{e}_1)}/\overline{u_2^2}$",
                                 COLORS["teal"], False)):
        ax.plot([1.0, 3.0], [y0, y0], color=COLORS["muted"], lw=1.0)
        ax.plot([1.0, 3.0], [y0, y0], "o", color=COLORS["ink"])
        for x0 in (1.0, 3.0):
            _arrow(ax, (x0, y0), (x0 + 0.55, y0) if along else (x0, y0 + 0.45), col, 2.4)
        ax.text(2.0, y0 - 0.17, "r", ha="center")
        ax.text(3.75, y0 + 0.05, lab, fontsize=9, color=col, va="center")
    return _bare(ax, (0.6, 8.4), (0.1, 2.2), "the two correlation coefficients of isotropic turbulence (cf. Fig. 12.9)")


def free_shear_sketch(ax, flow: str = "plane_jet", stations=(1.0, 2.0, 3.0, 4.0), constants: dict | None = None):
    """A free turbulent shear flow with its self-similar mean profiles at several stations (cf. Fig. 12.13).

    flow: "plane_jet", "plane_wake" or "shear_layer".  The profile shapes come from ``ch12.free_shear_flow`` with
    **labelled illustrative constants** (not the book's table): the picture shows the growth laws, not measured numbers.
    """
    c = dict(ILLUSTRATIVE_JET) if constants is None else dict(constants)
    c.setdefault("dxi80_coeff", 0.085)
    kw = {"plane_jet": dict(d=0.05, U0=1.0, rho_s=1.0, rho=1.0), "plane_wake": dict(theta=0.02, U_inf=1.0),
          "shear_layer": dict(U1=1.0, U2=0.3)}[flow]
    xs = np.linspace(0.3, 4.6, 60)
    width = np.array([float(ch12.free_shear_flow(flow, x, 0.0, constants=c, **kw)["width"]) for x in xs])
    scale = 2.2 if flow != "plane_wake" else 3.0
    ax.fill_between(xs, -scale * width, scale * width, color=COLORS["grid"], alpha=0.8)
    ax.plot(xs, width, color=COLORS["muted"], ls="--", lw=1.0)
    ax.plot(xs, -width, color=COLORS["muted"], ls="--", lw=1.0)
    half = 1.15 * scale * width[-1]
    y = np.linspace(-half, half, 161)
    ref = {"plane_jet": 0.0, "plane_wake": 1.0, "shear_layer": 0.0}[flow]   # the speed drawn on the station line
    profile = lambda x: np.asarray(ch12.free_shear_flow(flow, x, y, constants=c, **kw)["U"], dtype=float) - ref  # noqa: E731
    amp = 0.55 / float(np.max(np.abs(profile(stations[0]))))   # one drawing scale for all stations: the decay is visible
    for x in stations:
        ax.plot(x + amp * profile(x), y, color=COLORS["accent"], lw=2.0)
        ax.plot([x, x], [-half, half], color=COLORS["muted"], lw=0.6)
    if flow == "plane_jet":
        for sgn in (1, -1):
            for x in (1.5, 3.2):
                w = float(np.interp(x, xs, width))
                _arrow(ax, (x, sgn * 3.6 * w), (x, sgn * 2.5 * w), COLORS["teal"], 1.4)
        ax.text(0.35, 0.92 * half, "ambient fluid is entrained", fontsize=8.5, color=COLORS["teal"])
    ex = ch12.free_shear_exponents(flow)
    name = flow.replace("_", " ")
    ax.text(0.35, -0.95 * half, rf"{name}: width $\propto x^{{{ex['width']}}}$, velocity scale $\propto x^{{{ex['velocity']}}}$",
            fontsize=9)
    ax.set_xlabel("x (downstream)")
    ax.set_ylabel("y")
    return _bare(ax, (0.2, 5.3), (-1.05 * half, 1.05 * half), f"self-similar {name} (illustrative; cf. Fig. 12.13)")


def wall_layers_sketch(ax, Re_tau: float = 2000.0, Pi: float = 0.2):
    """The layers of a wall-bounded turbulent flow on the semi-log profile U⁺(y⁺) (cf. Figs. 12.16 and 12.18).

    Shape: ``ch12.composite_profile_plus`` with κ = 0.41, B = 5.0 (a common textbook pair) and an illustrative wake
    strength; bands at y⁺ = 5, y⁺ = 30 and y/δ = 0.15 (nominal boundaries).
    """
    yp = np.geomspace(0.3, Re_tau, 300)
    U = ch12.composite_profile_plus(yp, Re_tau, kappa=KAPPA, B=B_LOG, Pi=Pi)
    bands = ((0.3, 5.0, COLORS["rose"], "viscous\nsublayer"), (5.0, 30.0, COLORS["amber"], "buffer\nlayer"),
             (30.0, 0.15 * Re_tau, COLORS["teal"], "logarithmic\n(overlap) layer"), (0.15 * Re_tau, Re_tau, COLORS["blue"], "wake\nregion"))
    for lo, hi, col, lab in bands:
        ax.axvspan(lo, hi, color=col, alpha=0.10)
        ax.text(np.sqrt(lo * hi), 1.0, lab, ha="center", va="bottom", fontsize=8, color=col)
    ax.semilogx(yp, U, color=COLORS["accent"], lw=2.5, label="mean profile (composite, illustrative)")
    ys = yp[yp < 14]
    ax.semilogx(ys, ys, color=COLORS["rose"], ls=":", label=r"$U^+=y^+$")
    yl = yp[yp > 8]
    ax.semilogx(yl, ch12.log_law(yl, kappa=KAPPA, B=B_LOG), color=COLORS["ink"], ls="--", lw=1.0,
                label=rf"$U^+=\frac{{1}}{{\kappa}}\ln y^++B$ ($\kappa$ = {KAPPA}, B = {B_LOG}: a common textbook pair)")
    ax.set(xlim=(0.3, Re_tau), ylim=(0, float(U[-1]) + 3), xlabel=r"$y^+=y\,u_*/\nu$", ylabel=r"$U^+=U/u_*$",
           title=rf"layers of a wall flow at $\delta^+$ = {Re_tau:.0f} (cf. Fig. 12.18)")
    ax.legend(fontsize=7.5, loc="upper left")
    return ax


def surface_layer_sketch(ax, L_M: float = -20.0, top: float = 100.0, seed: int = 1):
    """Unstable daytime surface layer: forced convection below |L_M|, free convection with thermal plumes above (cf. Fig. 12.21)."""
    absL = abs(L_M)
    ax.axhspan(0.0, absL, color=COLORS["accent"], alpha=0.08)
    ax.axhspan(absL, top, color=COLORS["blue"], alpha=0.08)
    ax.axhline(absL, color=COLORS["ink"], ls="--", lw=1.0)
    ax.text(0.2, absL * 1.03, r"$z=|L_M|$: buoyancy and shear production equal in size", fontsize=8.5)
    z = np.linspace(0.5, top, 200)
    ax.plot(1.2 + 0.9 * np.log(z / 0.5) / np.log(top / 0.5), z, color=COLORS["accent"], lw=2.5)
    for zz in (4.0, 9.0, 15.0):
        _arrow(ax, (1.2, zz), (1.2 + 0.9 * np.log(zz / 0.5) / np.log(top / 0.5), zz), COLORS["accent"], 1.0)
    ax.text(1.25, 0.55 * absL, "forced convection\n(shear makes the turbulence;\nnear-logarithmic wind)", fontsize=8.5, color=COLORS["accent"])
    rng = np.random.default_rng(seed)
    for x0 in (4.2, 5.6, 7.0, 8.4):
        zz = np.linspace(0.6 * absL, top * rng.uniform(0.8, 0.97), 60)
        w = 0.12 + 0.35 * (zz - zz[0]) / (zz[-1] - zz[0])
        wob = 0.12 * np.sin(zz / 7.0 + rng.uniform(0, 6))
        ax.fill_betweenx(zz, x0 - w + wob, x0 + w + wob, color=COLORS["orange"], alpha=0.35, lw=0)
        _arrow(ax, (x0 + wob[30], zz[30]), (x0 + wob[45], zz[45]), COLORS["orange"], 1.5)
    ax.text(4.0, 0.9 * top, "free convection: thermal plumes", fontsize=8.5, color=COLORS["blue"])
    for x0 in np.linspace(3.6, 9.0, 7):
        _arrow(ax, (x0, 0.5), (x0, 5.0), COLORS["rose"], 1.0)
    ax.text(5.0, 6.0, r"upward heat flux $\overline{wT'}>0$ from the warm ground", fontsize=8, color=COLORS["rose"])
    ax.set(xlim=(0, 9.6), ylim=(0, top), xticks=[], ylabel="z [m] (illustrative)",
           title=rf"unstable surface layer, $L_M$ = {L_M:.0f} m (cf. Fig. 12.21)")
    ax.grid(False)
    return ax


def plume_sketch(ax, U: float = 5.0, w_rms: float = 0.5, Lambda_t: float = 10.0, x_max: float = 400.0):
    """Time-averaged smoke plume behind a chimney: width ∝ x near the source, ∝ √x far away (cf. Fig. 12.27).

    The envelope is ``ch12.smoke_plume_width`` (Taylor's formula with t = x/U) — linear near, square root far: the
    printed caption of the book's figure names the two regimes the other way round (slip #13).
    """
    import matplotlib.patches as mp

    x = np.linspace(0.0, x_max, 300)
    Z = ch12.smoke_plume_width(x, U, w_rms, Lambda_t)
    ax.fill_between(x, -Z, Z, color=COLORS["muted"], alpha=0.35, label=r"time-averaged plume, $\pm Z_{rms}$")
    xn = x[x < 2.5 * U * Lambda_t]
    ax.plot(xn, w_rms / U * xn, color=COLORS["ink"], ls=":", lw=1.2, label=r"near: $Z_{rms}=(w_{rms}/U)\,x$")
    ax.plot(x, w_rms * np.sqrt(2.0 * Lambda_t * x / U), color=COLORS["ink"], ls="--", lw=1.0,
            label=r"far: $Z_{rms}=w_{rms}\sqrt{2\Lambda_t x/U}$")
    ax.axvline(U * Lambda_t, color=COLORS["amber"], lw=1.0)
    ax.text(U * Lambda_t * 1.03, -0.9 * float(Z[-1]) * 1.5, r"$x=U\Lambda_t$", color=COLORS["amber"], fontsize=9)
    zc = 1.5 * float(Z[-1])
    ax.add_patch(mp.Rectangle((-0.02 * x_max, -zc), 0.02 * x_max, zc, color=COLORS["ink"]))
    _arrow(ax, (0.05 * x_max, 0.85 * zc), (0.2 * x_max, 0.85 * zc), COLORS["accent"], 1.5)
    ax.text(0.05 * x_max, 0.9 * zc, "wind U", color=COLORS["accent"], fontsize=9, va="bottom")
    ax.set(xlim=(-0.03 * x_max, x_max), ylim=(-zc, zc), xlabel="x [m] downwind", ylabel="z [m]",
           title="smoke plume: a wedge, then a parabola (cf. Fig. 12.27)")
    ax.legend(fontsize=8, loc="lower right")
    return ax


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(4, 2, figsize=(15, 19))
    scales_sketch(ax[0, 0])
    parcel_sketch(ax[0, 1])
    stress_element_sketch(ax[1, 0])
    fg_geometry_sketch(ax[1, 1])
    free_shear_sketch(ax[2, 0])
    wall_layers_sketch(ax[2, 1])
    surface_layer_sketch(ax[3, 0])
    plume_sketch(ax[3, 1])
    save(fig, out, "drawings")
    fig2, ax2 = plt.subplots(1, 3, figsize=(16, 4.4))
    for a, flow in zip(ax2, ("plane_jet", "plane_wake", "shear_layer")):
        free_shear_sketch(a, flow)
    save(fig2, out, "drawings_free_shear")
    print(f"drew {len(__all__)} sketches: {', '.join(__all__)}")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
