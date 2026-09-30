"""§9.3 (C03, C04): the Blasius solution with two independent solvers (Toepfer scaling, solve_bvp), Fig. 9.5 (u/U = f'(eta)) and Fig. 9.6
((v/U) sqrt(Re_x) vs eta), the collapse of dimensional profiles at three x, and every Blasius number computed (never typed).

Run: ``.venv/Scripts/python.exe scripts/ch09_blasius.py --no-show``   Figures -> outputs/ch09/c04_blasius_profiles.png, c04_blasius_collapse.png.
"""
from __future__ import annotations

import numpy as np

from ch09_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch09_boundary_layers as ch09


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    c = ch09.blasius_constants()
    bvp = ch09.falkner_skan(0.0, method="bvp")
    print("Blasius numbers (computed):")
    for k in ("fpp0", "eta99", "delta_star", "theta", "H", "v_inf", "cf_coeff", "cd_coeff"):
        print(f"  {k:<11} = {c[k]:.6f}")
    print(f"  bvp vs Toepfer f''(0): {abs(bvp['fpp0'] - c['fpp0']):.1e};  theta - 2 f''(0) = {c['theta'] - 2 * c['fpp0']:.1e}")
    for em in (8.0, 12.0, 16.0):
        print(f"  truncation eta_max = {em:>4}: f''(0) = {ch09.falkner_skan(0.0, eta_max=em, n=11)['fpp0']:.10f}")
    eta = np.linspace(0, 7, 400)
    f, fp, fpp = ch09.blasius_profile(eta)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].plot(fp, eta, color=COLORS["accent"], label="u/U = f'(eta)")
    ax[0].axhline(c["eta99"], color=COLORS["rose"], lw=1, ls="--", label=f"eta_99 = {c['eta99']:.3f} (root)")
    ax[0].axhline(c["delta_star"], color=COLORS["teal"], lw=1, ls="--", label=f"delta* = {c['delta_star']:.3f}")
    ax[0].axhline(c["theta"], color=COLORS["blue"], lw=1, ls="--", label=f"theta = {c['theta']:.3f}")
    ax[0].set_xlabel("u / U")
    ax[0].set_ylabel("eta = y sqrt(U / nu x)")
    ax[0].legend(fontsize=8)
    ax[0].set_title("Fig. 9.5-like: Blasius profile")
    ax[1].plot((eta * fp - f) / 2.0, eta, color=COLORS["orange"])
    ax[1].axvline(c["v_inf"], color=COLORS["muted"], ls=":", label=f"plateau {c['v_inf']:.4f}")
    ax[1].set_xlabel("(v / U) sqrt(Re_x)")
    ax[1].set_ylabel("eta")
    ax[1].legend(fontsize=8)
    ax[1].set_title("Fig. 9.6-like: wall-normal velocity")
    save(fig, out, "c04_blasius_profiles")
    U, nu = 1.0, 1e-5
    fig2, ax2 = plt.subplots(1, 2, figsize=(11, 4))
    y = np.linspace(0, 0.012, 300)
    for x, col in zip((0.1, 0.5, 1.0), (COLORS["accent"], COLORS["teal"], COLORS["orange"])):
        fl = ch09.blasius_fields(x, y, U, nu)
        ax2[0].plot(fl["u"], y * 1e3, color=col, label=f"x = {x} m")
        ax2[1].plot(fl["u"], fl["eta"], color=col, lw=2 - 0.4 * (x > 0.1))
    ax2[0].set_xlabel("u [m/s]")
    ax2[0].set_ylabel("y [mm]")
    ax2[0].legend()
    ax2[0].set_title("dimensional profiles differ with x")
    ax2[1].set_xlabel("u / U")
    ax2[1].set_ylabel("eta")
    ax2[1].set_title("one curve after the similarity scaling")
    spread = ch09.similarity_collapse_error("blasius", 0.0, 0.5)
    print(f"collapse spread at the right exponents: {spread:.1e};  with m = 0.4: {ch09.similarity_collapse_error('blasius', 0.0, 0.4):.2f}")
    save(fig2, out, "c04_blasius_collapse")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
