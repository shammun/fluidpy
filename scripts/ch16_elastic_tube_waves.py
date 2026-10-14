"""Chapter 16, §16.3 — waves in a viscous liquid inside a thin elastic tube (16.101)–(16.137): the characteristic
quadratic (16.136) in X = k²c₀²/ω², both roots (trap T7), wave speed and attenuation against the Womersley number, the
inviscid limit c → c₀ of (16.62), an independent Bessel-free determinant (ours), and a viscoelastic wall (16.137).
Time convention e^{i(kx − ωt)} (trap T1).

Run: ``.venv/Scripts/python.exe scripts/ch16_elastic_tube_waves.py --no-show``   Figures -> outputs/ch16/viscous_waves_*.png
"""
from __future__ import annotations

import numpy as np
from ch16_common import C_E, C_P, C_REF, C_U, C_WALL, Timer, finish, parse_args, save, setup

from fluidpy import ch16_biofluid_mechanics as ch16


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    A_ = np.geomspace(1.0, 100.0, 60 if args.fast else 150)
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    for nu_, ls in ((0.5, "-"), (0.3, "--")):
        cp = [ch16.elastic_tube_wave(a, 1.0, nu_)["c_over_c0"] for a in A_]
        dp = [ch16.elastic_tube_wave(a, 1.0, nu_)["damping_per_wavelength"] for a in A_]
        ax[0].semilogx(A_, cp, color=C_P, ls=ls, label=f"pressure wave, ν̂ = {nu_}")
        ax[1].semilogx(A_, dp, color=C_P, ls=ls, label=f"pressure wave, ν̂ = {nu_}")
    co = [ch16.elastic_tube_wave(a, 1.0, 0.5, branch="other")["c_over_c0"] for a in A_]
    do = [ch16.elastic_tube_wave(a, 1.0, 0.5, branch="other")["damping_per_wavelength"] for a in A_]
    ax[0].semilogx(A_, co, color=C_REF, ls=":", label="other root (ν̂ = 0.5): fast, wall-dominated")
    ax[1].semilogx(A_, do, color=C_REF, ls=":", label="other root (ν̂ = 0.5)")
    ax[0].axhline(1.0, color=C_WALL, lw=0.8)
    ax[0].set(xlabel="Womersley number α", ylabel="c/c₀", yscale="log",
              title="(16.136): viscosity slows the pulse; c → c₀ (16.62) as α → ∞")
    ax[0].legend(fontsize=8)
    ax[1].set(xlabel="Womersley number α", ylabel="amplitude left after one wavelength, e^{−2πk₂/k₁}",
              title="attenuation per wavelength")
    ax[1].legend(fontsize=8)
    save(fig, out, "viscous_waves")
    for a in (2.0, 5.0, 8.0, 20.0):
        w = ch16.elastic_tube_wave(a, 1.0, 0.5)
        rt = ch16.elastic_tube_roots(a, 0.5)
        print(f"α = {a:5.1f}: pressure root X = {w['X']:.5f}, c/c₀ = {w['c_over_c0']:.4f}, k₂/k₁ = {w['k2_over_k1']:.4f}, "
              f"left after one λ {w['damping_per_wavelength']:.3f}; other root X = {rt['other']:.4f}; "
              f"|residual (16.136)| = {abs(ch16.elastic_tube_characteristic(rt['pressure'], a, 0.5)):.1e}")
    al = np.array([20.0, 40.0, 80.0])
    g1 = [1.0 - ch16.elastic_tube_wave(a, 1.0, 0.5)["c_over_c0"] for a in al]
    g2 = [ch16.elastic_tube_wave(a, 1.0, 0.5)["k2_over_k1"] for a in al]
    print(f"  large α slopes: 1 − c/c₀ ∝ α^{np.polyfit(np.log(al), np.log(g1), 1)[0]:.3f}, k₂/k₁ ∝ α^{np.polyfit(np.log(al), np.log(g2), 1)[0]:.3f} (expect −1)")
    with Timer("Bessel-free determinant at 4 α (×3)"):
        for a in (2.0, 5.0, 10.0, 20.0):
            rt = ch16.elastic_tube_roots(a, 0.5)
            d0 = abs(ch16.elastic_tube_determinant_numeric(rt["pressure"], a, 0.5))
            d1 = abs(ch16.elastic_tube_determinant_numeric(rt["other"], a, 0.5))
            dn = abs(ch16.elastic_tube_determinant_numeric(1.1 * rt["pressure"], a, 0.5))
            print(f"  α = {a:4.1f}: numeric determinant at the pressure root {d0:.1e}, at the other root {d1:.1e}, "
                  f"10 % off the root {dn:.1e}")
    for Er in (1.0, 1.0 - 0.05j, 1.0 - 0.2j):
        w = ch16.elastic_tube_wave(8.0, 1.0, 0.5, E_ratio=Er)
        print(f"  viscoelastic wall (16.137), E*/E = {Er}: c/c₀ = {w['c_over_c0']:.4f}, left after one λ "
              f"{w['damping_per_wavelength']:.3f}")
    c0 = ch16.moens_korteweg_speed(0.45e6, 1.6e-3, 1060.0, 11e-3)
    w = ch16.elastic_tube_wave(8.0, c0, 0.5, omega=2 * np.pi * 66 / 60)
    print(f"  dimensional (ours): c₀ = {c0:.3f} m/s, α = 8 -> c = {w['c']:.3f} m/s, wavelength {2 * np.pi / w['k1']:.2f} m, "
          f"e-folding distance 1/k₂ = {1 / w['k2']:.1f} m")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
