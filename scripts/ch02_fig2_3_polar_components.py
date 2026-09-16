"""Chapter 2, §2.2, Example 2.1 and Fig. 2.3 (our drawing): a plane vector u resolved in Cartesian (x1, x2) and polar
(r, θ) axes; the polar C = [[cos θ, −sin θ], [sin θ, cos θ]] is Eq. (2.5) with j ∈ {r, θ}.

Run: ``.venv/Scripts/python.exe scripts/ch02_fig2_3_polar_components.py --no-show``
Figure → outputs/ch02/fig2_3_polar_components.png.
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

    from ch02_drawings import arrow2d, draw_frame_2d
    from fluidpy import ch02_cartesian_tensors as ch02
    from fluidpy.core.style import COLORS, use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    u = np.array([1.0, 2.0])
    theta = np.deg2rad(30.0)
    C = ch02.rotation_matrix_2d(theta)
    u_r, u_th = ch02.polar_components(*u, theta)  # Example 2.1
    assert np.allclose([u_r, u_th], ch02.transform_vector(u, C))  # = Eq. (2.5) with the polar C

    fig, axs = plt.subplots(1, 2, figsize=(11, 4.8))
    ax = axs[0]
    ax.set_aspect("equal")
    ax.set_xlim(-0.6, 2.4)
    ax.set_ylim(-0.4, 2.4)
    draw_frame_2d(ax, np.eye(2), COLORS["teal"], length=2.2, prime="", names=("$x_1$", "$x_2$"))
    draw_frame_2d(ax, C, COLORS["orange"], length=2.2, prime="", names=("$r$", "$\\theta$"))
    arrow2d(ax, (0, 0), u, color=COLORS["accent"], label="$\\mathbf{u}$", lw=2.6)
    # Cartesian projections (teal dashed) and polar projections (orange dashed)
    ax.plot([u[0], u[0]], [0, u[1]], color=COLORS["teal"], ls="--", lw=1)
    ax.plot([0, u[0]], [u[1], u[1]], color=COLORS["teal"], ls="--", lw=1)
    pr = u_r * C[:, 0]
    ax.plot([pr[0], u[0]], [pr[1], u[1]], color=COLORS["orange"], ls="--", lw=1)
    ax.plot([0, pr[0]], [0, pr[1]], color=COLORS["orange"], lw=3, alpha=0.4)
    ax.text(pr[0] * 0.5 + 0.1, pr[1] * 0.5 - 0.2, f"$u_r$ = {u_r:.3f}", color=COLORS["orange"])
    ax.text(u[0] + 0.08, u[1] * 0.5, f"$u_2$ = {u[1]:.0f}", color=COLORS["teal"])
    ax.text(u[0] * 0.4, -0.25, f"$u_1$ = {u[0]:.0f}", color=COLORS["teal"])
    ax.set_title(f"Example 2.1: u = (1, 2) at θ = 30°: $u_r$ = {u_r:.3f}, $u_\\theta$ = {u_th:.3f}")
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")

    ax = axs[1]
    th = np.linspace(0, 2 * np.pi, 361)
    ur_all, ut_all = ch02.polar_components(u[0], u[1], th)
    ax.plot(np.rad2deg(th), ur_all, label="$u_r = u_1\\cos\\theta + u_2\\sin\\theta$", color=COLORS["orange"])
    ax.plot(np.rad2deg(th), ut_all, label="$u_\\theta = -u_1\\sin\\theta + u_2\\cos\\theta$", color=COLORS["accent"])
    ax.plot(np.rad2deg(th), np.hypot(ur_all, ut_all), label="$|\\mathbf{u}|$ (unchanged)", color=COLORS["muted"], ls=":")
    ax.axvline(30, color=COLORS["grid"])
    ax.set_xlabel("polar angle θ [deg]")
    ax.set_ylabel("component")
    ax.set_title("components change with the frame, the length does not")
    ax.legend(fontsize=8.5)
    fig.savefig(out / "fig2_3_polar_components.png", bbox_inches="tight")

    print(f"Example 2.1: u = {u}, θ = 30°: u_r = {u_r:.6f}, u_θ = {u_th:.6f}; C =\n{np.round(C, 4)}")
    print(f"round trip: {ch02.cartesian_from_polar(u_r, u_th, theta)}; |u| = {np.linalg.norm(u):.6f} = {np.hypot(u_r, u_th):.6f}")
    print(f"saved fig2_3_polar_components.png in {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
