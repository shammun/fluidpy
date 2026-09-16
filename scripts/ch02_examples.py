"""Chapter 2, Examples 2.1–2.6 worked numerically and symbolically (the notebook's worked numbers), plus heatmaps of
∇·u and (∇×u)₃ for Example 2.3's two fields.

Run: ``.venv/Scripts/python.exe scripts/ch02_examples.py --no-show``
Figure → outputs/ch02/example2_3_div_curl.png.
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
    import sympy as sp

    from fluidpy import ch02_cartesian_tensors as ch02
    from fluidpy.core.style import use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    print("=== Example 2.1: Cartesian → polar components ===")
    th = sp.Symbol("theta", real=True)
    u1, u2 = sp.symbols("u1 u2", real=True)
    Cs = sp.Matrix([[sp.cos(th), -sp.sin(th)], [sp.sin(th), sp.cos(th)]])
    up = sp.simplify(Cs.T * sp.Matrix([u1, u2]))  # Eq. (2.5): u' = Cᵀ u
    print(f"  u_r = {up[0]},  u_θ = {up[1]}")
    for deg in (0, 30, 90, 137):
        ur, ut = ch02.polar_components(1.0, 2.0, np.deg2rad(deg))
        print(f"  u = (1, 2), θ = {deg:3d}°: u_r = {ur:+.4f}, u_θ = {ut:+.4f}, |u| = {np.hypot(ur, ut):.4f}")

    print("=== Example 2.2: traction in a shear flow ===")
    for a in (1.0, -1.0):
        ex = ch02.example_2_2(a, np.deg2rad(30))
        print(f"  a = {a:+.0f}: f = {np.round(ex['f'], 4)}, |f| = {ex['magnitude']:.4f}, θ = {ex['angle_deg']:.0f}°, "
              f"τ'11 = {ex['sigma_n']:+.4f} (a sin 2φ = {ex['sigma_n_formula']:+.4f}), τ'12 = {ex['tau_s']:+.4f} (a cos 2φ = {ex['tau_s_formula']:+.4f})")
    A = np.array([[1.0, 2.0], [3.0, 4.0]])
    n = np.array([0.6, 0.8])
    print(f"  non-symmetric check: n·τ = {ch02.traction(A, n)} vs τ·n = {A @ n}  (they differ: (2.15) contracts the first index)")

    print("=== Example 2.3: div and curl of a x and b × x ===")
    a_s = sp.Symbol("a", real=True)
    b_s = sp.symbols("b1:4", real=True)
    fr, fb = ch02.radial_field(a_s), ch02.solid_body_rotation_field(b_s)
    print(f"  ∇·(a x) = {fr.div_expr}, ∇×(a x) = {fr.curl_expr};  ∇·(b × x) = {fb.div_expr}, ∇×(b × x) = {fb.curl_expr}")
    ex3 = ch02.example_2_3(2.0, (0.3, -0.2, 1.0), n=16)
    print(f"  grid (a = 2, b = (0.3, −0.2, 1)): ∇·(a x) = {ex3['div_radial_grid']:.12f}, max|∇×(a x)| = {ex3['curl_radial_max']:.1e}, "
          f"max|∇·(b×x)| = {ex3['div_rotation_max']:.1e}, ∇×(b × x) = {np.round(ex3['curl_rotation_grid'], 12)}")
    print(f"  index route: ∂_i(a x_i) = a δ_ii = 3a → {ch02.expand_indices('a delta_ii')};  ε_ijk ε_ljk = 2δ_il → "
          f"{ch02.expand_indices('eps_ijk eps_ljk')[0, 0]}, {ch02.expand_indices('eps_ijk eps_ljk')[0, 1]}")

    print("=== Example 2.4: principal axes of the plane shear strain rate ===")
    ex4 = ch02.example_2_4(1.0)
    print(f"  λ = {ex4['lam']}, b¹ = {np.round(ex4['C'][:, 0], 6)}, b² = {np.round(ex4['C'][:, 1], 6)}, rotation {ex4['angle_deg']:.1f}°")
    print(f"  S' = Cᵀ S C = {np.round(ex4['S_prime'], 12).tolist()}  (Γ taken as S12; for u1(x2), S12 = ½ du1/dx2)")

    print("=== Example 2.5: divergence from the integral definition (2.32) ===")
    conv = ch02.integral_definition_convergence("divergence")
    print(f"  h = {conv['h']}, |(1/V)∮n·Q dA − ∇·Q| = {np.round(conv['err'], 7)}, observed order {conv['order']:.3f}")
    print(f"  linear field a x at x0 = (1, 2, 3), h = 0.5: {ch02.integral_divergence(ch02.radial_field(1.0), [1, 2, 3], 0.5):.12f} (exact 3)")

    print("=== Example 2.6: curl from the integral definition (2.35) ===")
    conv = ch02.integral_definition_convergence("curl_component")
    print(f"  h = {conv['h']}, |Γ/A − (∇×u)₃| = {np.round(conv['err'], 8)}, observed order {conv['order']:.3f}")
    fb1 = ch02.solid_body_rotation_field((0.3, -0.2, 1.0))
    for nvec in ([1, 0, 0], [0, 1, 0], [0, 0, 1]):
        v = ch02.integral_curl_component(fb1, [0.1, 0.2, 0.3], nvec, 0.2)
        print(f"  b × x, n = {nvec}: Γ/A = {v:+.12f} (exact 2 n·b = {2 * np.dot(nvec, [0.3, -0.2, 1.0]):+.1f})")
    print(f"  comma notation (2.36): {ch02.comma_to_partial('eps_ijk u_k,j')} → component 1: {ch02.expand_indices('eps_ijk u_k,j')[0]}")

    # figure: Example 2.3 fields in the plane
    g = ch02.grid2d((-1, 1), 41)
    ur = ch02.evaluate_field(ch02.radial_field(1.0, dim=2), g)
    ub = ch02.evaluate_field(ch02.solid_body_rotation_field(1.0, dim=2), g)
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.8))
    for ax, U, title, val in ((axs[0], ur, "u = a x (a = 1): ∇·u = 2a in the plane, (∇×u)₃ = 0", ch02.divergence(ur, g.h)),
                              (axs[1], ub, "u = b × x (b = e₃): ∇·u = 0, (∇×u)₃ = 2b", ch02.curl(ub, g.h))):
        im = ax.pcolormesh(g.X, g.Y, val, shading="auto", cmap="coolwarm", vmin=-2.2, vmax=2.2)
        sk = 4
        ax.quiver(g.X[::sk, ::sk], g.Y[::sk, ::sk], U[0][::sk, ::sk], U[1][::sk, ::sk], color="k", scale=25, width=0.003)
        ax.set_aspect("equal")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("x₁")
        ax.set_ylabel("x₂")
        fig.colorbar(im, ax=ax, label="∇·u  or  (∇×u)₃  [1/s]")
    fig.savefig(out / "example2_3_div_curl.png", bbox_inches="tight")
    print(f"saved example2_3_div_curl.png in {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
