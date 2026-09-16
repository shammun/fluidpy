"""Chapter 2, §2.11, Example 2.4 and Fig. 2.8 (our drawing): the plane strain rate S = [[0, Γ], [Γ, 0]] of a shear flow
u1(x2), its 45° principal axes b¹, b², a material square deforming under S alone (stretch along b¹, squeeze along b²)
and Mohr's circle for S.

Run: ``.venv/Scripts/python.exe scripts/ch02_fig2_8_principal_axes.py --no-show``
Figure → outputs/ch02/fig2_8_principal_axes.png.
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
    ap.add_argument("--Gamma", type=float, default=0.5, help="S12 = Γ [1/s]")
    ap.add_argument("--t", type=float, default=0.6, help="deformation time [s]")
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

    ex = ch02.example_2_4(args.Gamma)
    S, C, lam = ex["S"], ex["C"], ex["lam"]
    sq0 = ch02.deform_square(S, 0.0, 30, half_width=0.6, boundary_only=True)
    sq_S = ch02.deform_square(S, args.t, 30, half_width=0.6, boundary_only=True)  # pure strain S: no spin
    G = ex["G"]  # the shear flow's velocity gradient [[0, 2Γ], [0, 0]] with S as its symmetric part
    sq_G = ch02.deform_square(G, args.t, 30, half_width=0.6, boundary_only=True)
    centre, radius = ch02.mohr_circle_2d(S)

    fig, axs = plt.subplots(1, 3, figsize=(14, 4.8))
    ax = axs[0]
    ax.set_aspect("equal")
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    draw_frame_2d(ax, np.eye(2), COLORS["teal"], length=1.3, prime="", names=("$x_1$", "$x_2$"))
    draw_frame_2d(ax, C, COLORS["orange"], length=1.3, names=("$x_1'$ (b¹)", "$x_2'$ (b²)"))
    ax.set_title(f"Fig. 2.8 idea: principal axes at {ex['angle_deg']:.0f}°; λ = ({lam[0]:+.2f}, {lam[1]:+.2f}) 1/s")
    ax = axs[1]
    ax.set_aspect("equal")
    ax.plot(*sq0, color=COLORS["muted"], ls="--", label="t = 0")
    ax.plot(*sq_S, color=COLORS["orange"], lw=2, label=f"under S alone, t = {args.t} s")
    ax.plot(*sq_G, color=COLORS["accent"], lw=1.5, ls=":", label="under the full shear G = S + A")
    for k, col in ((0, COLORS["blue"]), (1, COLORS["rose"])):
        arrow2d(ax, (0, 0), 0.9 * np.exp(lam[k] * args.t) * C[:, k], color=col, lw=2)
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.6, 1.6)
    ax.set_title("stretch along b¹ by e^{Γt}, squeeze along b² by e^{−Γt}; G adds a spin")
    ax.legend(fontsize=8)
    ax = axs[2]
    th = np.linspace(0, 2 * np.pi, 361)
    ax.plot(centre + radius * np.cos(th), radius * np.sin(th), color=COLORS["accent"])
    ax.plot([lam.min(), lam.max()], [0, 0], "o", color=COLORS["orange"])
    phi = np.linspace(0, np.pi, 181)
    cur = ch02.stress_vs_angle(S, phi)
    ax.plot(cur["sigma_n"], cur["tau_s"], ".", ms=2, color=COLORS["teal"])
    ax.plot([S[0, 0]], [S[0, 1]], "s", color=COLORS["ink"], label="φ = 0: (S₁₁, S₁₂)")
    ax.set_aspect("equal")
    ax.set_xlabel("normal rate σ_n [1/s]")
    ax.set_ylabel("shear rate τ_s [1/s]")
    ax.set_title(f"Mohr circle: centre {centre:.2f}, radius {radius:.2f} = (λ_max − λ_min)/2")
    ax.legend(fontsize=8)
    fig.savefig(out / "fig2_8_principal_axes.png", bbox_inches="tight")

    print(f"Example 2.4 (Γ = {args.Gamma}): S =\n{S}\nλ = {lam}; b¹ = {np.round(C[:, 0], 4)}, b² = {np.round(C[:, 1], 4)}; "
          f"C = 45° rotation (det {np.linalg.det(C):+.3f})\nS' = Cᵀ S C =\n{np.round(ex['S_prime'], 12)}")
    print(f"principal angle ½ atan2(2S12, S11 − S22) = {np.rad2deg(ch02.principal_angle_2d(S)):.2f}°; "
          f"characteristic polynomial coefficients {ch02.characteristic_polynomial(S)} → roots {np.roots(ch02.characteristic_polynomial(S))}")
    omega3 = ch02.vector_from_antisymmetric(ch02.antisymmetric_part(G))  # angular velocity of the element = ½(∇×u)₃
    vort3 = ch02.vector_from_antisymmetric(ch02.rotation_tensor(G))  # the book's R = G − Gᵀ (3.17): its vector is ∇×u
    print(f"velocity gradient G = [[0, 2Γ],[0, 0]]: symmetric part = S (S12 = ½·2Γ = Γ), antisymmetric part A = ½R with vector "
          f"{omega3:+.3f} = ½(∇×u)3 (element spin); rotation tensor R = G − Gᵀ has vector {vort3:+.3f} = (∇×u)3 (vorticity)")
    print(f"saved fig2_8_principal_axes.png in {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
