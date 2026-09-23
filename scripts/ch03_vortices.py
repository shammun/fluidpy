"""§3.5 vortices: solid-body rotation vs the irrotational vortex (Figs. 3.15/3.16 by tracking a small element with
``pathline``: it spins in one, keeps its orientation on average in the other), the polar vorticity (3.23) against the
four-leg circulation of a small sector, Γ(r), the δ-core (3.27), and the Rankine/Gaussian profiles (3.28)–(3.29) with
the Gaussian maximum at 1.1209σ (C13, C14, N36–N44).

Run: ``.venv/Scripts/python.exe scripts/ch03_vortices.py --no-show``
Figures → outputs/ch03/fig3_15_16_elements.png, c14_vortex_profiles.png, c13_mean_vorticity_disc.png.
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

from fluidpy import ch03_kinematics as ch03  # noqa: E402
from fluidpy.core.style import COLORS  # noqa: E402
from tools.convergence import observed_order  # noqa: E402


def _pair_turn(P0: np.ndarray, P1: np.ndarray) -> float:
    """Mean signed turn [rad] of the two sides BC and BA of a quadrilateral A B C D between two snapshots."""
    out = 0.0
    for i, j in ((2, 1), (0, 1)):
        s0, s1 = P0[:, i] - P0[:, j], P1[:, i] - P1[:, j]
        out += 0.5 * np.arctan2(s0[0] * s1[1] - s0[1] * s1[0], s0 @ s1)
    return float(out)


def track_element(kind: str, r0: float = 1.0, half: float = 0.08, t_end: float = np.pi / 2, n: int = 5, **params):
    """A small square element A B C D centred at (r0, 0) carried by the vortex ``kind`` (``vortex_profile``
    parameters), at n times up to t_end, and the angle of a rigid paddle wheel riding at its centre: the paddle turns
    at the element's spin ½ω_z (§3.4), so its angle is ∫½ω_z dt along the centre's path line.

    Returns (times, corners (n, 2, 4), paddle_angle (n,) [rad], short-time pair turn rate [rad/s]). The last value is
    the mean turn of two tiny perpendicular material lines through the element centre over dt = 1e-3 s, divided by
    dt — the *instantaneous* spin (½ω₃ of §3.4). Over long times the lines stop being perpendicular and their mean
    turn is **not** the spin (in the irrotational vortex the azimuthal line follows the orbit while the radial line is
    sheared back); lines through different points (e.g. two sides of a finite element) do not measure it either."""
    u = ch03.vortex_velocity_field(kind, **params)
    c0 = np.array([[r0 - half, r0 - half, r0 + half, r0 + half], [half, -half, -half, half]])  # A B C D (as Fig. 3.15)
    ts = np.linspace(0.0, t_end, n)
    corners = np.moveaxis(np.stack([ch03.pathline(u, c0[:, k], 0.0, ts) for k in range(4)], axis=2), 1, 0)  # (n, 2, 4)
    tf = np.linspace(0.0, t_end, 201)
    centre = ch03.pathline(u, [r0, 0.0], 0.0, tf)
    spin = 0.5 * ch03.vorticity(u, centre)[2]  # ½ω_z along the path (h-stencil at each point)
    from scipy.integrate import cumulative_trapezoid

    paddle_angle = np.interp(ts, tf, cumulative_trapezoid(spin, tf, initial=0.0))
    dt, eps = 1e-3, 1e-4  # a tiny cross through the centre: both lines at the same point, as the definition requires
    x0 = np.array([[r0, r0, r0 + eps, r0], [eps, 0.0, 0.0, 0.0]])  # "A" above, "B" = centre, "C" to the right
    x1 = np.stack([ch03.pathline(u, x0[:, k], 0.0, dt) for k in range(4)], axis=1)
    return ts, corners, paddle_angle, _pair_turn(x0, x1) / dt


def elements_figure(t_end: float = np.pi / 2):
    """Figs. 3.15/3.16 side by side: the element at several times in solid-body rotation (it spins as it revolves; S = 0)
    and in the irrotational vortex (it deforms; a paddle wheel riding with it keeps its heading), each with the paddle
    wheel turned by ∫½ω_z dt. Returns (fig, results)."""
    import matplotlib.pyplot as plt

    try:
        from scripts.ch03_drawings import paddle
    except ModuleNotFoundError:  # pragma: no cover
        from ch03_drawings import paddle

    fig, axs = plt.subplots(1, 2, figsize=(10, 4.8))
    out = {}
    for ax, kind, title in ((axs[0], "solid", "solid-body rotation (3.22): the element spins"),
                            (axs[1], "line", "irrotational vortex (3.25): it deforms, no spin")):
        ts, corners, ang, rate = track_element(kind, t_end=t_end, Gamma=2 * np.pi, sigma=1.0)
        out[kind] = (ts, ang, rate)
        for P, a, al in zip(corners, ang, np.linspace(0.35, 1.0, len(ts))):
            ax.fill(*P, color=COLORS["accent"], alpha=0.3 * al)
            ax.plot(*np.hstack([P, P[:, :1]]), color=COLORS["accent"], alpha=al)
            paddle(ax, P.mean(axis=1), a, 0.07, COLORS["orange"])
        th = np.linspace(0, 2 * np.pi, 200)
        for rr in (0.92, 1.08):
            ax.plot(rr * np.cos(th), rr * np.sin(th), color=COLORS["grid"], lw=1)
        ax.plot(0, 0, "+", color=COLORS["ink"])
        ax.set_aspect("equal")
        ax.set_xlim(-0.3, 1.3)
        ax.set_ylim(-0.3, 1.3)
        ax.set_title(title + f"\npaddle turned {np.degrees(ang[-1]):+.1f}° after t = {ts[-1]:.2f} s", fontsize=10)
        ax.set_xlabel("$x$ [m]")
    axs[0].set_ylabel("$y$ [m]")
    return fig, out


def profiles_figure(Gamma: float = 2 * np.pi, sigma: float = 1.0):
    """u_θ(r) and ω_z(r) of the four §3.5 profiles (same Γ, σ) with the Gaussian maximum marked."""
    import matplotlib.pyplot as plt

    r = np.linspace(1e-3, 4.0, 400) * sigma
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4))
    cols = dict(solid=COLORS["muted"], line=COLORS["rose"], rankine=COLORS["teal"], gaussian=COLORS["accent"])
    for kind in ch03.VORTEX_KINDS:
        u, w = ch03.vortex_profile(kind, r, Gamma, sigma)
        a1.plot(r / sigma, u, color=cols[kind], label=kind)
        a2.plot(r / sigma, np.broadcast_to(w, r.shape), color=cols[kind], label=kind)
    rm = ch03.gaussian_vortex_max_radius(sigma)
    um = ch03.gaussian_vortex(rm, Gamma, sigma)[0]
    a1.plot(rm / sigma, um, "o", color=COLORS["accent"])
    a1.annotate(f"max at r = {rm / sigma:.4f}σ", (rm / sigma, um), xytext=(10, 10), textcoords="offset points", fontsize=9)
    a1.set_ylim(0, 1.6 * Gamma / (2 * np.pi * sigma))
    a1.set_xlabel(r"$r/\sigma$")
    a1.set_ylabel(r"$u_\theta$ [m/s]")
    a1.set_title(rf"$\Gamma$ = {Gamma:.3f} m²/s, $\sigma$ = {sigma} m")
    a1.legend(fontsize=8)
    a2.set_ylim(-0.1, 2.5 * Gamma / (np.pi * sigma ** 2) / 2)
    a2.set_xlabel(r"$r/\sigma$")
    a2.set_ylabel(r"$\omega_z$ [1/s]")
    a2.set_title("vorticity: a core, then nothing")
    return fig


def mean_vorticity_figure(B: float = 1.0, omega0: float = 1.0):
    """Mean vorticity in a disc Γ(r)/(πr²) vs r (log–log): 2B/r² for the line vortex (slope −2, the δ-core of (3.27)),
    2ω₀ for solid-body rotation. Returns (fig, slope)."""
    import matplotlib.pyplot as plt

    r = np.logspace(-2, 1, 30)
    mv_line = ch03.mean_vorticity_in_disc(lambda rr: ch03.line_vortex(rr, B), r)
    mv_solid = ch03.mean_vorticity_in_disc(lambda rr: ch03.solid_body_rotation(rr, omega0), r)
    slope = observed_order(r, mv_line)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.loglog(r, mv_line, color=COLORS["rose"], label=f"line vortex (slope {slope:.2f})")
    ax.loglog(r, mv_solid, color=COLORS["muted"], label="solid body (2ω₀)")
    ax.set_xlabel("disc radius $r$ [m]")
    ax.set_ylabel(r"$\Gamma(r)/\pi r^2$ [1/s]")
    ax.set_title("Eq. (3.27): infinite vorticity with a finite integral")
    ax.legend(fontsize=8)
    return fig, slope


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch03"))
    ap.add_argument("--no-show", action="store_true")
    args = ap.parse_args()
    t0 = time.perf_counter()
    import matplotlib

    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from fluidpy.core.style import use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    fig, rots = elements_figure()
    fig.savefig(out / "fig3_15_16_elements.png", bbox_inches="tight")
    for kind, (ts, ang, rate) in rots.items():
        print(f"{kind:5s} element (Γ = 2π, r₀ = 1 m): short-time pair turn rate {rate:+.6f} rad/s (½ω_z = "
              f"{0.5 * ch03.vortex_profile(kind, 1.0)[1]:+.1f}); paddle angle ∫½ω_z dt after {ts[-1]:.3f} s = {ang[-1]:+.6f} rad")
    profiles_figure().savefig(out / "c14_vortex_profiles.png", bbox_inches="tight")
    fig3, slope = mean_vorticity_figure()
    fig3.savefig(out / "c13_mean_vorticity_disc.png", bbox_inches="tight")
    # C13 numbers
    print(f"C13: (3.23) solid body ω₀ = 1: {ch03.polar_vorticity_z(None, lambda r: ch03.solid_body_rotation(r, 1.0), 0.7):.10f}; "
          f"line vortex: {ch03.polar_vorticity_z(None, lambda r: ch03.line_vortex(r, 1.0), 0.7):.1e}")
    for d in (0.1, 0.05, 0.025):
        leg = ch03.annular_sector_circulation("gaussian", 0.6, d, d / 0.6, legs=True)
        rc = 0.6 + d / 2
        print(f"     sector Δr = {d}: Γ/area = {leg['total'] / leg['area']:.8f} vs ω_z(r̄) = {ch03.gaussian_vortex(rc, 2 * np.pi, 1)[1]:.8f}")
    print(f"     line vortex sector Γ_ABCD = {ch03.annular_sector_circulation('line', 1.0, 0.1, 0.2):.1e}; Γ(r) for r = 0.5, 1, 2: "
          f"{[round(ch03.circulation_circle(lambda r: ch03.line_vortex(r, 1.0), rr), 12) for rr in (0.5, 1.0, 2.0)]} (2πB = {2 * np.pi:.12f})")
    print(f"     solid body Γ(r = 1) = {ch03.circulation_circle(lambda r: ch03.solid_body_rotation(r, 1.0), 1.0):.10f}, off-centre "
          f"circle (0.3, 0.2): {ch03.circulation_circle(lambda r: ch03.solid_body_rotation(r, 1.0), 1.0, center=(0.3, 0.2)):.10f}")
    print(f"     line vortex, circle not enclosing the axis: {ch03.circulation_circle(lambda r: ch03.line_vortex(r, 1.0), 0.4, center=(1.0, 0.0)):.1e}")
    print(f"     mean vorticity in a disc: log–log slope {slope:.4f} (line vortex)")
    # C14 numbers
    rm = ch03.gaussian_vortex_max_radius()
    print(f"C14: Γ = 2π, σ = 1: Rankine inside u_θ = r, ω = {ch03.rankine_vortex(0.5, 2 * np.pi, 1)[1]:.1f} s⁻¹, max "
          f"{ch03.rankine_vortex(1.0, 2 * np.pi, 1)[0]:.3f} m/s at r = 1 m; Gaussian max {ch03.gaussian_vortex(rm, 2 * np.pi, 1)[0]:.4f} m/s "
          f"at r = {rm:.10f} m (Lambert W: {ch03.gaussian_vortex_max_radius(method='lambertw'):.10f}); x* = r²/σ² = {rm ** 2:.10f}")
    print(f"     Gaussian Γ(r = 2σ) = {ch03.circulation_circle('gaussian', 2.0):.6f} = 2π(1 − e⁻⁴) = {2 * np.pi * (1 - np.exp(-4)):.6f}")
    print(f"figures → {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
