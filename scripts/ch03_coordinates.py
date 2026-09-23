"""§3.1: one-, two-, three-dimensional flows and coordinate systems — our Fig. 3.1c/d (a developing pipe profile and its
cross-section average, which stays the same: the 1-D description) and our Fig. 3.3 (local unit vectors of plane
polar, cylindrical and spherical coordinates at one point), with the conversion numbers (N03, N05).

Run: ``.venv/Scripts/python.exe scripts/ch03_coordinates.py --no-show``
Figures → outputs/ch03/fig3_1_section_average.png, fig3_3_unit_vectors.png.
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


def section_average_figure(zs=(0.0, 0.5, 1.5, 5.0)):
    """Profiles u(r, z) of ``pipe_profile`` at several stations with their area averages (all equal to U_mean)."""
    import matplotlib.pyplot as plt

    r = np.linspace(-1, 1, 201)
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for z, col in zip(zs, (COLORS["teal"], COLORS["blue"], COLORS["accent"], COLORS["orange"])):
        ubar = ch03.cross_section_average(lambda rr, zz: ch03.pipe_profile(rr, zz), 1.0, z)
        ax.plot(ch03.pipe_profile(r, z), r, color=col, label=f"z = {z} m (ū = {ubar:.6f} m/s)")
    ax.axvline(1.0, color=COLORS["ink"], ls=":", label="1-D description ū = U_mean")
    ax.set_xlabel("$u$ [m/s]")
    ax.set_ylabel("$r/R$")
    ax.set_title("a 2-D flow and its 1-D average (Fig. 3.1c/d idea)")
    ax.legend(fontsize=8)
    return fig


def unit_vectors_figure(P=(1.0, 1.2, 0.9)):
    """3-D axes with the point P and the cylindrical (e_R, e_φ, e_z) and spherical (e_r, e_θ, e_φ) unit vectors there."""
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(6, 5.4))
    ax = fig.add_subplot(111, projection="3d")
    P = np.asarray(P, float)
    R, phi, z = ch03.cylindrical_from_cartesian(*P)
    r, th, ph = ch03.spherical_from_cartesian(*P)
    Ec, Es = ch03.unit_vectors_cylindrical(phi), ch03.unit_vectors_spherical(th, ph)
    for E, names, col, off in ((Ec, ("$e_R$", r"$e_\varphi$", "$e_z$"), COLORS["teal"], 0.0),
                               (Es, ("$e_r$", r"$e_\theta$", r"$e_\varphi$"), COLORS["orange"], 0.02)):
        for e, nm in zip(E, names):
            ax.quiver(*(P + off), *(0.6 * e), color=col, arrow_length_ratio=0.15)
            ax.text(*(P + off + 0.68 * e), nm, color=col)
    ax.plot([0, P[0]], [0, P[1]], [0, P[2]], color=COLORS["muted"], ls=":")
    ax.plot([P[0], P[0]], [P[1], P[1]], [0, P[2]], color=COLORS["muted"], ls=":")
    ax.scatter(*P, color=COLORS["ink"])
    for k, lbl in enumerate("xyz"):
        v = np.zeros(3)
        v[k] = 1.8
        ax.quiver(0, 0, 0, *v, color=COLORS["ink"], arrow_length_ratio=0.05)
        ax.text(*(1.05 * v), f"${lbl}$")
    ax.set_title(f"P: R = {R:.3f}, φ = {np.degrees(phi):.1f}°, r = {r:.3f}, θ = {np.degrees(th):.1f}°")
    return fig


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
    section_average_figure().savefig(out / "fig3_1_section_average.png", bbox_inches="tight")
    unit_vectors_figure().savefig(out / "fig3_3_unit_vectors.png", bbox_inches="tight")
    print(f"N03: Poiseuille u = U(1 − r²/R²) averages to {ch03.cross_section_average(lambda r, z: 1 - r ** 2, 1.0):.12f} U")
    P = np.array([1.0, 1.2, 0.9])
    print(f"N05: P = {P}: cylindrical {np.round(ch03.cylindrical_from_cartesian(*P), 6)}, spherical "
          f"{np.round(ch03.spherical_from_cartesian(*P), 6)}")
    rng = np.random.default_rng(0)
    X = rng.normal(size=(3, 10000))
    back = np.array(ch03.cartesian_from_spherical(*ch03.spherical_from_cartesian(*X)))
    backc = np.array(ch03.cartesian_from_cylindrical(*ch03.cylindrical_from_cartesian(*X)))
    print(f"     round trips on 10⁴ random points: spherical {np.abs(back - X).max():.1e}, cylindrical {np.abs(backc - X).max():.1e}")
    u = ch03.rigid_body_velocity([0, 0, 0], [0, 0, 2.0], P)
    print(f"     solid body Ω = 2e_z at P: (u_R, u_φ, u_z) = {ch03.velocity_components(u, P, 'cylindrical').round(6)} "
          f"(ΩR = {2 * np.hypot(1.0, 1.2):.6f}); spherical {ch03.velocity_components(u, P, 'spherical').round(6)}")
    E = ch03.unit_vectors_spherical(0.7, 1.9)
    print(f"     spherical basis orthonormal: |EEᵀ − I| = {np.abs(E @ E.T - np.eye(3)).max():.1e}, det = {np.linalg.det(E):.12f}")
    print(f"figures → {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
