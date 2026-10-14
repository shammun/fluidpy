"""Chapter 16, §16.3 — reflection of pulse waves at an arterial junction (16.148)–(16.163): admittance Y = A/ρc,
reflection and transmission coefficients, the energy balance R² + (ΣY_j/Y₁)T² = 1 (ours), and a pulse run on a grid
through a fork (ours, characteristic transport at Courant number 1) for matched, closed-end-like and open-end-like
junctions.

Run: ``.venv/Scripts/python.exe scripts/ch16_junctions.py --no-show``   Figures -> outputs/ch16/junction_*.png
"""
from __future__ import annotations

import numpy as np
from ch16_common import C_E, C_P, C_REF, C_U, C_WALL, Timer, finish, parse_args, save, setup

from fluidpy import ch16_biofluid_mechanics as ch16


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch16.illustrative_inputs()
    rho = I["rho"]
    ratio = np.geomspace(0.05, 20, 300)
    R = (1 - ratio) / (1 + ratio)
    T = 2 / (1 + ratio)
    rng = np.random.default_rng(0)
    worst = max(abs(ch16.junction_energy(rng.uniform(0.1, 3), rng.uniform(0.1, 3, 2))["total"] - 1) for _ in range(1000))
    print(f"energy check over 1000 random junctions (seed 0): max |R² + (ΣY/Y₁)T² − 1| = {worst:.1e}")
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
    ax[0].semilogx(ratio, R, color=C_P, label="R = (Y₁ − ΣY)/(Y₁ + ΣY) (16.153)")
    ax[0].semilogx(ratio, T, color=C_U, label="T = 2Y₁/(Y₁ + ΣY) = 1 + R")
    ax[0].semilogx(ratio, R ** 2, color=C_E, ls="--", label="energy reflected R²")
    ax[0].axvline(1.0, color=C_REF, ls=":")
    ax[0].set(xlabel="ΣY_j / Y₁", ylabel="coefficient", title="closed-end-like ← matched → open-end-like")
    ax[0].legend(fontsize=8)
    a1, c1 = 5e-3, 5.0
    Y1 = float(ch16.tube_admittance(np.pi * a1 ** 2, rho, c1))
    cases = {"matched": a1 * 2 ** -0.5, "closed-end-like": 2.0e-3, "open-end-like": 4.5e-3}
    with Timer("three fork runs"):
        for (name, ad), col in zip(cases.items(), (C_REF, C_P, C_U)):
            Yd = float(ch16.tube_admittance(np.pi * ad ** 2, rho, c1))
            run = ch16.network_pulse_run(0.5, [0.4, 0.4], c1, [c1, c1], Y1, [Yd, Yd], width=0.02, amplitude=1.0e3,
                                         fast=args.fast)
            en = ch16.junction_energy(Y1, [Yd, Yd])
            print(f"  {name:16s}: daughters a = {ad * 1e3:.2f} mm, ΣY/Y₁ = {2 * Yd / Y1:.3f}: R measured {run['R_measured']:+.6f} "
                  f"(16.153) {run['R_theory']:+.6f}; T measured {run['T_measured']:.6f} (16.153) {run['T_theory']:.6f}; "
                  f"energy reflected {en['reflected']:.4f} + transmitted {en['transmitted']:.4f} = {en['total']:.12f}")
            ax[1].plot(run["probe_t"] * 1e3, run["probe_p"] / 1e3, color=col, label=f"{name}: p")
            ax[2].plot(run["probe_t"] * 1e3, run["probe_Q"] * 1e6, color=col, label=f"{name}: Q")
    ax[1].set(xlabel="t [ms]", ylabel="p at the probe [kPa]", title=f"parent probe at x = {run['x_probe']:.2f} m (16.159)")
    ax[1].legend(fontsize=8)
    ax[2].set(xlabel="t [ms]", ylabel="Q at the probe [cm³/s]", title="flow (16.160): the reflection flips sign in Q")
    ax[2].legend(fontsize=8)
    save(fig, out, "junction")
    w = ch16.junction_waveforms(-0.25, np.array([0.08, 0.16]), lambda s: np.exp(-((s - 0.06) / 0.02) ** 2), c1, 0.4,
                                Y1=Y1, P_I=1.0e3)
    print(f"  closed form (16.159)–(16.160) at x = −0.25 m, t = 80, 160 ms (R = 0.4): p = {np.round(w['p'], 2)} Pa, "
          f"Q = {np.round(w['Q'] * 1e6, 4)} cm³/s")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
