"""§11.5 (C06): double diffusion — the finger criterion (11.46), the regime map in (dT̄/dz, dS̄/dz) with the static-stability
line, and the roots of our free–free cubic (fingers: real σ; diffusive regime: complex σ).

Run: ``.venv/Scripts/python.exe scripts/ch11_double_diffusion.py --no-show``   Figure -> outputs/ch11/c06_double_diffusion.png.
Our seawater: α = 2e-4 1/K, β = 7.6e-4 (g/kg)⁻¹, κ = 1.4e-7, κ_s = 1.5e-9, ν = 1e-6 m²/s, layer d = 5 cm.
"""
from __future__ import annotations

import math

import numpy as np

from ch11_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch11_instability as ch11

AL, BE, KA, KS, NU, D = 2e-4, 7.6e-4, 1.4e-7, 1.5e-9, 1e-6, 0.05


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    sf = ch11.salt_finger_unstable(0.01, 0.002, D)
    print(f"thermocline example dT/dz = 0.01 K/m, dS/dz = 0.002 (g/kg)/m, d = 5 cm: (11.46) lhs = {sf['lhs']:.1f}, "
          f"margin = {sf['margin']:.1f}, density stable = {sf['density_stable']}, R_ρ = {sf['R_rho']:.4f}")
    dmin = __import__("scipy.optimize", fromlist=["brentq"]).brentq(
        lambda d: ch11.salt_finger_unstable(0.01, 0.002, d)["margin"], 1e-3, 0.05)
    print(f"thinnest finger-unstable layer: d = {dmin * 100:.3f} cm")
    print(f"τ = κ_s/κ = {KS / KA:.6f}; RA_FREE_FREE = 27π⁴/4 = {ch11.RA_FREE_FREE:.6f}")
    K2 = math.pi ** 2 / 2
    for Ra, Rs in ((1000.0, 2000.0), (-2e4, -1.9e6)):
        print(f"cubic roots (K² = π²/2, Pr = 7, τ = {KS / KA:.5f}) at (Ra, Rs) = ({Ra:g}, {Rs:g}):",
              np.round(ch11.double_diffusive_sigma(K2, Ra, Rs, 7.0, KS / KA), 5))
    print("τ = 1, Rs = 0, Ra(§11.5) = −2000 vs Bénard free–free:", ch11.double_diffusive_sigma(K2, -2000.0, 0.0, 1.0, 1.0)[0],
          ch11.benard_free_free_sigma(math.sqrt(K2), 2000.0, 1.0)[0])
    Tz = np.linspace(-0.02, 0.02, 81)
    Sz = np.linspace(-0.006, 0.006, 81)
    col = {"stable": 0, "fingers": 1, "diffusive": 2, "overturning": 3}
    M = np.array([[col[ch11.salt_finger_regime(t, s, D)["regime"]] for t in Tz] for s in Sz])
    # the top-heavy but stable strip of (11.46) on the row dS/dz = 0: Rs − Ra = −Ra < 27π⁴/4 ⇔ dT/dz > −27π⁴νκ/(4gαd⁴)
    t_crit = -ch11.RA_FREE_FREE / float(ch11.thermal_rayleigh_signed(1.0, D, AL, NU, KA))
    row0 = [ch11.salt_finger_regime(t, 0.0, D)["regime"] for t in Tz]
    strip = Tz[(Tz < 0) & (np.array(row0) == "stable")]
    print(f"dS/dz = 0 row: top-heavy layers are stable for {t_crit:.5f} < dT/dz < 0 K/m ((11.46): Rs − Ra < 27π⁴/4); "
          f"map points labelled stable there: {len(strip)} (dT/dz = {strip.min():.4f} … {strip.max():.4f}), "
          f"overturning for dT/dz ≤ {Tz[np.array(row0) == 'overturning'].max():.4f}")
    assert all(r == "stable" for t, r in zip(Tz, row0) if t_crit < t < 0)
    assert all(r == "overturning" for t, r in zip(Tz, row0) if t < t_crit)
    fig, ax = plt.subplots(1, 2, figsize=(14, 5))
    from matplotlib.colors import ListedColormap

    cm = ListedColormap([COLORS["teal"], COLORS["amber"], COLORS["accent"], COLORS["rose"]])
    ax[0].pcolormesh(AL * Tz, BE * Sz, M, cmap=cm, vmin=-0.5, vmax=3.5, shading="nearest")
    ax[0].plot(AL * Tz, AL * Tz, color=COLORS["ink"], lw=1.5, label="α dT/dz = β dS/dz (static neutral)")
    lhs_line = (ch11.RA_FREE_FREE * NU / (9.80665 * D ** 4) + AL / KA * Tz) * KS  # β dS/dz on (11.46)
    ax[0].plot(AL * Tz, lhs_line, color=COLORS["blue"], ls="--", label="(11.46) finger threshold")
    ax[0].plot([AL * t_crit, 0.0], [0.0, 0.0], color=COLORS["ink"], lw=3, solid_capstyle="butt",
               label="top-heavy but stable: Rs − Ra < 27π⁴/4")
    ax[0].set(xlim=(AL * Tz[0], AL * Tz[-1]), ylim=(BE * Sz[0], BE * Sz[-1]), xlabel="α dT̄/dz [1/m]", ylabel="β dS̄/dz [1/m]",
              title="regimes: teal stable · amber fingers · purple diffusive · rose overturning")
    ax[0].legend(fontsize=7, loc="lower right")
    Rs = np.linspace(-2e6, 2e6, 200)
    for Ra, c in ((1000.0, COLORS["amber"]), (-2e4, COLORS["accent"])):
        s = np.array([ch11.double_diffusive_sigma(K2, Ra, r, 7.0, KS / KA)[0] for r in Rs])
        ax[1].plot(Rs, s.real, color=c, label=f"Re σ, Ra = {Ra:g}")
    ax[1].axhline(0, color=COLORS["muted"], lw=1)
    ax[1].set(xlabel="Rs", ylabel="max Re σ [κ/d²]", title="leading root of the free–free cubic (ours)", ylim=(-50, 60))
    ax[1].legend(fontsize=8)
    save(fig, out, "c06_double_diffusion")
    # diffusive row: the gradients that give the tested (Ra, Rs) = (−2×10⁴, −1.9×10⁶) (σ = 9.650 ± 70.389i) in this layer —
    # Ra and Rs are linear in the gradients, so divide by their value per unit gradient (≈ −0.228 K/m, −0.0612 (g/kg)/m;
    # the weaker pair (−0.01, −0.0028) used before is damped, σ = −6.54 ± 16.2i: "stable", not diffusive)
    t_dif = -2e4 / float(ch11.thermal_rayleigh_signed(1.0, D, AL, NU, KA))
    s_dif = -1.9e6 / float(ch11.salinity_rayleigh(1.0, D, BE, NU, KS))
    for name, (t, s) in dict(fingers=(0.01, 0.002), diffusive=(t_dif, s_dif), stable=(0.01, -0.001),
                             overturning=(-0.01, 0.0), top_heavy_weak=(-0.005, 0.0)).items():
        r = ch11.salt_finger_regime(t, s, D)
        print(f"  {name:14s} → {r['regime']:12s} R_ρ = {r['R_rho']:.3f}, σ_max = {r['sigma_max']:.4g}  ({r['text']})")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
