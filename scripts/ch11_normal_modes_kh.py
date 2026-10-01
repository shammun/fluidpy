"""§11.1–§11.3 (C01, C02): the four mechanical analogies of Fig. 11.1, normal-mode growth (11.1), Kelvin–Helmholtz (11.18)
with and without surface tension (Ex. 11.1), the vortex sheet (11.20), the static interface (11.19) against Ch. 7 (7.95), the
KH stability boundary ΔU_min(k) and the mixing energy of p. 482.

Run: ``.venv/Scripts/python.exe scripts/ch11_normal_modes_kh.py --no-show``   Figures -> outputs/ch11/c01_potential_wells.png,
c02_kelvin_helmholtz.png.  Our inputs: air over water ρ₁ = 1.2, ρ₂ = 1000 kg/m³, σ_s = 0.074 N/m.
"""
from __future__ import annotations

import numpy as np

from ch11_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch11_instability as ch11
from fluidpy.core import waves as WV

RHO_AIR, RHO_WATER, SIGMA_S = 1.2, 1000.0, 0.074


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    # --- Fig. 11.1 analogue: ball in four potentials, small and large kicks
    fig, ax = plt.subplots(1, 4, figsize=(15, 3.6))
    xs = np.linspace(-1.6, 1.6, 200)
    for a, shape in zip(ax, ("bowl", "cap", "plane", "dimple")):
        V = ch11.potential_well_demo(shape, 0.1)["V"]
        a.plot(xs, V(xs), color=COLORS["muted"])
        for x0, col in ((0.3, COLORS["teal"]), (1.2, COLORS["rose"])):
            r = ch11.potential_well_demo(shape, x0)
            a.plot(r["x"], V(r["x"]) + 0.02, color=col, lw=1.5, label=f"x₀ = {x0}: {'escapes' if r['escaped'] else 'stays'}")
            print(f"{shape:7s} x0 = {x0}: escaped = {r['escaped']}, final x = {r['x'][-1]: .4f}")
        a.set(title=shape, xlabel="x [–]", ylabel="V(x) [–]", ylim=(-1.0, 1.0))
        a.legend(fontsize=7)
    save(fig, out, "c01_potential_wells")

    # --- KH: growth rates, boundary, c-plane
    k = np.geomspace(10.0, 3000.0, 400)
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    for dU, col in ((5.0, COLORS["teal"]), (8.0, COLORS["orange"]), (12.0, COLORS["rose"])):
        ax[0].semilogx(k, ch11.kh_growth_rate(k, dU, 0.0, RHO_AIR, RHO_WATER), color=col, ls="--", label=f"ΔU = {dU} m/s, σ_s = 0")
        ax[0].semilogx(k, ch11.kh_growth_rate(k, dU, 0.0, RHO_AIR, RHO_WATER, surface_tension=SIGMA_S), color=col,
                       label=f"ΔU = {dU} m/s, σ_s = {SIGMA_S}")
        kc = ch11.kh_critical_k(dU, 0.0, RHO_AIR, RHO_WATER)
        band = ch11.kh_unstable_band(dU, RHO_AIR, RHO_WATER, surface_tension=SIGMA_S)
        print(f"ΔU = {dU:4.1f} m/s: k_c (no σ_s) = {kc:8.2f} 1/m (λ = {2 * np.pi / kc * 100:.2f} cm); band with σ_s: "
              f"{band[0]:.1f}–{band[1]:.1f} 1/m")
    ax[0].set(xlabel="k [1/m]", ylabel="growth rate k c_i [1/s]", title="KH air over water (11.18) and Ex. 11.1")
    ax[0].legend(fontsize=7)
    ms = ch11.kh_min_shear(RHO_AIR, RHO_WATER, surface_tension=SIGMA_S)
    print(f"minimum wind shear (Ex. 11.1, h → ∞): ΔU_min = {ms['dU_min']:.3f} m/s at λ = {ms['wavelength'] * 100:.3f} cm")
    print(f"Rayleigh–Taylor cut-off (water over air): λ_c = {ch11.rayleigh_taylor_cutoff(SIGMA_S, RHO_WATER, RHO_AIR) * 100:.3f} cm")
    ax[1].loglog(k, ch11.kh_stability_boundary(k, RHO_AIR, RHO_WATER, surface_tension=SIGMA_S), color=COLORS["accent"],
                 label="ΔU_min(k), σ_s = 0.074")
    ax[1].loglog(k, ch11.kh_stability_boundary(k, RHO_AIR, RHO_WATER), color=COLORS["muted"], ls="--", label="σ_s = 0")
    ax[1].plot(ms["k_star"], ms["dU_min"], "o", color=COLORS["rose"], label=f"min {ms['dU_min']:.2f} m/s")
    ax[1].set(xlabel="k [1/m]", ylabel="ΔU [m/s]", title="stability boundary (unstable above)")
    ax[1].legend(fontsize=8)
    # c-plane: the two roots collide as ΔU passes the boundary at k*
    for dU in np.linspace(0.0, 10.0, 41):
        cp, cm = ch11.kh_phase_speed(ms["k_star"], dU, 0.0, RHO_AIR, RHO_WATER, surface_tension=SIGMA_S)
        col = COLORS["rose"] if np.imag(cp) > 0 else COLORS["teal"]
        ax[2].plot([np.real(cp), np.real(cm)], [np.imag(cp), np.imag(cm)], ".", color=col)
    ax[2].set(xlabel="c_r [m/s]", ylabel="c_i [m/s]", title="c± at k* as ΔU rises 0 → 10 m/s")
    save(fig, out, "c02_kelvin_helmholtz")

    # --- limits and checks
    kk = 50.0
    cp, cm = ch11.kh_phase_speed(kk, 0.0, 0.0, RHO_AIR, RHO_WATER, g=9.81)
    print(f"(11.19) static interface k = {kk}: c = ±{np.real(cp):.6f} m/s; Ch. 7 (7.95) ω/k = "
          f"{WV.interface_omega(kk, RHO_AIR, RHO_WATER, 9.81) / kk:.6f} m/s")
    print(f"(11.20) vortex sheet U1 = 1, U2 = 3: c = {ch11.vortex_sheet_c(1.0, 3.0)}; (11.18) with ρ₁ = ρ₂: "
          f"{ch11.kh_phase_speed(2.0, 1.0, 3.0, 1.0, 1.0)}")
    print("kh_residuals (k = 2, U1 = 3, U2 = −1, ρ = 1, 3):", ch11.kh_residuals(2.0, 3.0, -1.0, 1.0, 3.0))
    mix = ch11.kh_mixing_energy(1.0, 1.0, 1.0)
    print(f"mixing (p. 482): E_f/E_i = {mix['ratio']:.6f} (2/3), momentum {mix['M_i']} → {mix['M_f']}")
    print("normal_mode_growth('kh', 1, U1 = 6, U2 = 0, ρ = 1, 3, g = 10):",
          ch11.normal_mode_growth("kh", 1.0, U1=6, U2=0, rho1=1, rho2=3, g=10))
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
