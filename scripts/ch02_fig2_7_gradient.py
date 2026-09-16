"""Chapter 2, §2.9 and Fig. 2.7 (our drawing): level curves of φ = x² + y²/4, the gradient arrows perpendicular to
them, and ∂φ/∂n = ∇φ·n at a probe point as the direction n turns — maximal along ∇φ, zero along the contour.

Run: ``.venv/Scripts/python.exe scripts/ch02_fig2_7_gradient.py --no-show``
Figure → outputs/ch02/fig2_7_gradient.png.
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
    ap.add_argument("--n", type=int, default=81, help="grid nodes per direction")
    args = ap.parse_args()
    t0 = time.perf_counter()
    import matplotlib

    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from ch02_drawings import arrow2d
    from fluidpy import ch02_cartesian_tensors as ch02
    from fluidpy.core.style import COLORS, use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    g = ch02.grid2d(((-2, 2), (-2, 2)), args.n)
    phi = g.X ** 2 + g.Y ** 2 / 4.0
    grad = ch02.gradient(phi, g.h)  # (∇φ)_i = ∂φ/∂x_i, exact for this quadratic (central differences)
    assert np.allclose(grad[0], 2 * g.X) and np.allclose(grad[1], g.Y / 2)
    probe = np.array([1.0, 1.0])
    gp = np.array([2 * probe[0], probe[1] / 2])  # ∇φ at the probe = (2, 0.5)
    angles = np.linspace(0, 2 * np.pi, 361)
    N = np.stack([np.cos(angles), np.sin(angles)])
    dphidn = ch02.directional_derivative(gp[:, None] * np.ones((1, angles.size)), N)  # ∂φ/∂n = ∇φ·n for every n
    k_max = int(np.argmax(dphidn))

    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    ax = axs[0]
    cs = ax.contour(g.X, g.Y, phi, levels=12, colors=COLORS["muted"], linewidths=1)
    ax.clabel(cs, fmt="%.1f", fontsize=7)
    sk = max(1, args.n // 12)
    ax.quiver(g.X[::sk, ::sk], g.Y[::sk, ::sk], grad[0][::sk, ::sk], grad[1][::sk, ::sk], color=COLORS["teal"], scale=40, width=0.003)
    arrow2d(ax, probe, 0.35 * gp, color=COLORS["accent"], label="$\\nabla\\varphi$", lw=2.5)
    n_probe = np.array([np.cos(np.deg2rad(110)), np.sin(np.deg2rad(110))])
    arrow2d(ax, probe, 0.8 * n_probe, color=COLORS["orange"], label="$\\mathbf{n}$", lw=2)
    ax.plot(*probe, "o", color=COLORS["ink"])
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title("Fig. 2.7 idea: φ = x² + y²/4 — ∇φ (teal) is ⊥ to the level curves")
    ax = axs[1]
    ax.plot(np.rad2deg(angles), dphidn, color=COLORS["accent"])
    ax.axhline(np.linalg.norm(gp), color=COLORS["teal"], ls="--", label="$|\\nabla\\varphi|$ = max")
    ax.axhline(0, color=COLORS["grid"])
    ax.axvline(np.rad2deg(angles[k_max]), color=COLORS["teal"], ls=":", label=f"n ∥ ∇φ at {np.rad2deg(angles[k_max]):.1f}°")
    ax.set_xlabel("direction of n [deg]")
    ax.set_ylabel("∂φ/∂n = ∇φ·n")
    ax.set_title(f"at ({probe[0]:.0f}, {probe[1]:.0f}): ∇φ = ({gp[0]:.1f}, {gp[1]:.1f}), ∂φ/∂n ranges ± {np.linalg.norm(gp):.3f}")
    ax.legend()
    fig.savefig(out / "fig2_7_gradient.png", bbox_inches="tight")

    tang = np.array([-gp[1], gp[0]]) / np.linalg.norm(gp)
    print(f"∇φ at {probe} = {gp}; |∇φ| = {np.linalg.norm(gp):.4f}; ∂φ/∂n along the contour tangent = "
          f"{float(ch02.directional_derivative(gp, tang)):.2e}; along n = (1, 0): {float(ch02.directional_derivative(gp, [1, 0])):.4f}")
    print(f"grid: {g.shape}, h = {g.h}; max |gradient − exact| = {np.max(np.abs(grad[0] - 2 * g.X)):.2e} (quadratic: exact)")
    print(f"saved fig2_7_gradient.png in {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
