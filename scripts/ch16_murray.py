"""Chapter 16, §16.3 — Murray's law: the cost of a vessel (16.138), the optimal radius (16.139)–(16.141) with the
metabolic constant in the denominator (slip S15 corrected), the cube law at a fork (16.142)–(16.144), and the optimal
branching angles (16.145)–(16.147) checked against a numerical minimisation of the branch point (ours).

Run: ``.venv/Scripts/python.exe scripts/ch16_murray.py --no-show``   Figures -> outputs/ch16/murray_*.png
"""
from __future__ import annotations

import numpy as np
from ch16_common import C_E, C_P, C_REF, C_U, C_WALL, finish, parse_args, save, setup

from fluidpy import ch16_biofluid_mechanics as ch16


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch16.illustrative_inputs()
    mu, Q, b, L = I["mu"], I["Q_murray"], I["b_metabolic"], 0.05
    a_opt = ch16.murray_radius(Q, mu, b)
    a_pr = ch16.murray_radius(Q, mu, b, printed_16_140=True)
    ag = np.linspace(0.3 * a_opt, 3.0 * a_opt, 200001)
    a_grid = ag[np.argmin(ch16.murray_cost(ag, Q, L, mu, b))]
    tm = ch16.murray_terms(a_opt, Q, L, mu, b)
    print(f"Murray (ours): Q = {Q * 1e6} cm³/s, μ = {mu * 1e3} mPa s, upkeep b = {b} W/m³, L = {L * 100} cm")
    print(f"  optimal radius (16.140, corrected) a = [16μ/(π²b)]^(1/6) Q^(1/3) = {a_opt * 1e3:.5f} mm; grid search "
          f"{a_grid * 1e3:.5f} mm; printed (16.140) would give {a_pr * 1e3:.3f} mm (b in the numerator, slip S15)")
    print(f"  at the optimum: pumping {tm['pumping'] * 1e6:.4f} µW, upkeep {tm['upkeep'] * 1e6:.4f} µW (ratio "
          f"{tm['pumping'] / tm['upkeep']:.6f} = 1/2), total {tm['total'] * 1e6:.4f} µW = (3π/2)bLa² (16.141) "
          f"{ch16.murray_min_cost(Q, L, mu, b) * 1e6:.4f} µW")
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    aa = np.linspace(0.4 * a_opt, 2.5 * a_opt, 400)
    T = ch16.murray_terms(aa, Q, L, mu, b)
    ax[0].plot(aa * 1e3, T["pumping"] * 1e6, color=C_P, label="pumping 8μLQ²/(πa⁴)")
    ax[0].plot(aa * 1e3, T["upkeep"] * 1e6, color=C_WALL, label="upkeep bπa²L")
    ax[0].plot(aa * 1e3, T["total"] * 1e6, color=C_E, lw=2.2, label="total Ê₁ (16.138)")
    ax[0].axvline(a_opt * 1e3, color=C_REF, ls=":", label=f"optimum a = {a_opt * 1e3:.3f} mm")
    ax[0].set(xlabel="radius a [mm]", ylabel="power [µW]", ylim=(0, 3 * tm["total"] * 1e6),
              title="the cheapest vessel: pumping = ½ upkeep")
    ax[0].legend()
    # the fork
    a0 = a_opt
    A, C, D = np.array([0.0, 0.0]), np.array([0.06, 0.035]), np.array([0.055, -0.04])
    for split in (0.5, 0.7, 0.9):
        a1, a2 = ch16.murray_daughters(a0, split)
        ang = ch16.murray_angles(a0, a1, a2)
        bp = ch16.bifurcation_point(A, C, D, a0, a1, a2)
        print(f"  split Q₁/Q₀ = {split}: a₁ = {a1 * 1e3:.4f}, a₂ = {a2 * 1e3:.4f} mm (a₀³ − a₁³ − a₂³ = {a0 ** 3 - a1 ** 3 - a2 ** 3:.1e}); "
              f"(16.146) θ = {np.degrees(ang['theta']):.4f}°, φ = {np.degrees(ang['phi']):.4f}°; numerical optimum of B: "
              f"θ = {np.degrees(bp['theta']):.4f}°, φ = {np.degrees(bp['phi']):.4f}°; (16.147) θ = {np.degrees(ch16.murray_theta(a0, a1)):.4f}°")
        if split == 0.7:
            Bopt, a1p, a2p = bp["B"], a1, a2
    print(f"  symmetric fork: θ = φ = arccos 2^(−1/3) = {np.degrees(ch16.murray_symmetric_angle()):.4f}°")
    for P, a_, col in ((A, a0, C_U), (C, a1p, C_U), (D, a2p, C_U)):
        ax[1].plot([Bopt[0], P[0]], [Bopt[1], P[1]], color=col, lw=2 + 6 * a_ / a0, solid_capstyle="round")
    ax[1].plot(*Bopt, "o", color=C_P, label="optimal B (numerical, ours)")
    for P, lab in ((A, "A"), (C, "C"), (D, "D")):
        ax[1].annotate(lab, P, textcoords="offset points", xytext=(4, 4))
    ax[1].set(xlabel="x [m]", ylabel="y [m]", aspect="equal",
              title="70/30 fork: line width ∝ radius; angles satisfy (16.146)")
    ax[1].legend()
    save(fig, out, "murray")
    Y = ch16.junction_from_geometry(a0, [a0 * 2 ** (-1 / 3)] * 2, 5.0, 5.0, I["rho"])
    print(f"  a Murray-symmetric fork with equal wave speeds reflects R = {Y['R']:.4f} (ΣY/Y₁ = {Y['ratio']:.4f} = 2^(1/3))")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
