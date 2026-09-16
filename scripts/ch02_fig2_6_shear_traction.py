"""Chapter 2, §2.6, Example 2.2 and Fig. 2.6 (our drawing): the channel shear flow τ = [[0, a], [a, 0]], an element
whose normal is 30° from the flow, the traction f, and σ_n(φ) = a sin 2φ, τ_s(φ) = a cos 2φ for every cut angle.

Run: ``.venv/Scripts/python.exe scripts/ch02_fig2_6_shear_traction.py --no-show``
Figure → outputs/ch02/fig2_6_shear_traction.png.
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
    ap.add_argument("--a", type=float, default=1.0, help="shear stress a [Pa] (negative in the other half)")
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

    a = args.a
    ex = ch02.example_2_2(a, np.deg2rad(30.0))
    phi = np.linspace(0, np.pi, 361)
    curves = ch02.stress_vs_angle(ex["tau"], phi)
    assert np.allclose(curves["sigma_n"], a * np.sin(2 * phi)) and np.allclose(curves["tau_s"], a * np.cos(2 * phi))

    fig, axs = plt.subplots(1, 3, figsize=(14, 4.6), gridspec_kw={"width_ratios": [1, 1, 1.3]})
    # panel 1: parabolic profile (Poiseuille-like) with the element
    ax = axs[0]
    y = np.linspace(-1, 1, 41)
    ax.barh(y, 1 - y ** 2, height=0.04, color=COLORS["teal"], alpha=0.6)
    ax.axhline(1, color=COLORS["ink"], lw=3)
    ax.axhline(-1, color=COLORS["ink"], lw=3)
    ax.set_title("channel flow: τ = [[0, a],[a, 0]],\na > 0 above the centreline, a < 0 below", fontsize=9.5)
    ax.set_xlabel("$u_1(x_2)$")
    ax.set_ylabel("$x_2$")
    # panel 2: the element at 30° with n, f and the rotated axes
    ax = axs[1]
    ax.set_aspect("equal")
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.3, 1.3)
    draw_frame_2d(ax, np.eye(2), COLORS["teal"], length=1.1, prime="", names=("1", "2"))
    draw_frame_2d(ax, ex["C"], COLORS["orange"], length=1.1)
    n, f = ex["n"], ex["f"]
    s = np.array([-n[1], n[0]])
    ax.plot([-0.6 * s[0], 0.6 * s[0]], [-0.6 * s[1], 0.6 * s[1]], color=COLORS["ink"], lw=4)  # the cut
    arrow2d(ax, (0, 0), 0.8 * n, color=COLORS["ink"], label="$\\mathbf{n}$ (30°)", lw=1.5)
    arrow2d(ax, (0, 0), 0.9 * f / max(abs(a), 1e-12), color=COLORS["rose"], label="$\\mathbf{f} = \\mathbf{n}\\cdot\\tau$", lw=2.5)
    arrow2d(ax, (0, 0), 0.9 * ex["sigma_n"] / max(abs(a), 1e-12) * n, color=COLORS["blue"], lw=1.5, ls="--")
    arrow2d(ax, 0.9 * ex["sigma_n"] / max(abs(a), 1e-12) * n, 0.9 * ex["tau_s"] / max(abs(a), 1e-12) * s, color=COLORS["accent"], lw=1.5, ls="--")
    ax.set_title(f"Example 2.2 at φ = 30°: f = ({f[0]:.3f}, {f[1]:.3f}) Pa,\n|f| = {ex['magnitude']:.3f}, direction θ = {ex['angle_deg']:.0f}°", fontsize=9.5)
    # panel 3: σ_n and τ_s over all angles
    ax = axs[2]
    ax.plot(np.rad2deg(phi), curves["sigma_n"], color=COLORS["blue"], label="$\\sigma_n = \\tau'_{11} = a\\sin 2\\varphi$")
    ax.plot(np.rad2deg(phi), curves["tau_s"], color=COLORS["accent"], label="$\\tau_s = \\tau'_{12} = a\\cos 2\\varphi$")
    ax.axvline(30, color=COLORS["grid"])
    ax.plot([30], [ex["sigma_n"]], "o", color=COLORS["blue"])
    ax.plot([30], [ex["tau_s"]], "o", color=COLORS["accent"])
    ax.axvline(45, color=COLORS["rose"], ls=":", lw=1)
    ax.text(46, 0.9 * a, "45°: no shear\n(principal axes, Ex. 2.4)", fontsize=8.5, color=COLORS["rose"])
    ax.set_xlabel("angle φ of the normal [deg]")
    ax.set_ylabel("stress [Pa]")
    ax.set_title("the same stress pushes differently\non every plane: σ_n(φ), τ_s(φ)", fontsize=9.5)
    ax.legend(fontsize=8.5, loc="lower left")
    fig.savefig(out / "fig2_6_shear_traction.png", bbox_inches="tight")

    print(f"Example 2.2 (a = {a}): n = {np.round(ex['n'], 4)}, f = {np.round(ex['f'], 4)} Pa, |f| = {ex['magnitude']:.4f}, "
          f"direction {ex['angle_deg']:.1f}°")
    print(f"via (2.12): τ' = Cᵀ τ C =\n{np.round(ex['tau_rot'], 4)} -> normal {ex['sigma_n']:.4f} (√3a/2 = {np.sqrt(3) / 2 * a:.4f}), "
          f"shear {ex['tau_s']:.4f} (a/2 = {a / 2:.4f})")
    neg = ch02.example_2_2(-abs(a), np.deg2rad(30.0))
    print(f"a < 0: f = {np.round(neg['f'], 4)}, direction {neg['angle_deg']:.1f}°")
    print(f"Mohr circle (centre, radius) = {ch02.mohr_circle_2d(ex['tau'])}")
    print(f"saved fig2_6_shear_traction.png in {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
