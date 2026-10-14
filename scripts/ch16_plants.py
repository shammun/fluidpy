"""Chapter 16, §16.4 — water transport in plants: the water potential ψ = p − RTc + ρgz (16.221) (SI, and the book's
bar units with the ρgz/10⁶ conversion, slip S30), the tension in a transpiring xylem column with Poiseuille's viscous
drop (16.17), the xylem Reynolds number (slip S28 handled with our own numbers), capillary rise for comparison (ch01,
ours), and the phloem sieve-tube/sieve-plate pressure drop of Rand et al. (16.222).

Run: ``.venv/Scripts/python.exe scripts/ch16_plants.py --no-show``   Figures -> outputs/ch16/plants_*.png
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
    T, Hh, a, u = I["T_plant"], I["tree_height"], I["a_xylem"], I["u_xylem"]
    mu_w, nu_w = 1.0e-3, 1.0e-6
    print(f"R = {ch16.R_MOLAR:.6f} J/(mol K) = {ch16.R_MOLAR * 10:.4f} cm³ bar/(mol K) (the book's 83.141)")
    t = ch16.water_potential_terms(0.0, 100.0, T, 0.0)
    print(f"  0.1 M solution at the ground, atmospheric: osmotic term −RTc = {t['osmotic'] / 1e6:.4f} MPa; "
          f"book units: {ch16.water_potential_book_units(0.0, 1e-4, T, 0.0):.4f} bar")
    print(f"  pure water lifted 10 m at rest: ρgz = {ch16.water_potential_terms(0.0, 0.0, T, 10.0)['gravity'] / 1e6:.4f} MPa; "
          f"book units with ρgz/10⁶ (S30): {ch16.water_potential_book_units(0.0, 0.0, T, 1000.0, g_cm_s2=980.665):.4f} bar")
    Re = ch16.xylem_reynolds(u, a, nu_w)
    print(f"xylem conduit (ours): a = {a * 1e6:.0f} µm, u = {u * 1e3} mm/s, water ν = {nu_w} m²/s -> Re = u·2a/ν = {Re:.3f} "
          f"(Stokes flow); with a 10× larger ν the same Re would need a 10× faster sap (slip S28)")
    z = np.linspace(0, Hh, 121)
    col = ch16.xylem_column(z, psi_soil=I["psi_soil"], u=u, a=a, mu=mu_w, T=T)
    print(f"  {Hh:.0f} m tree, ψ_soil = {I['psi_soil'] / 1e6} MPa: crown pressure p = {col['p'][-1] / 1e6:.3f} MPa "
          f"(hydrostatic share {col['gravity'][-1] / 1e6:.3f}, viscous {col['viscous'][-1] / 1e6:.3f} MPa, gradient "
          f"{col['grad_viscous'] / 1e6:.4f} MPa/m); ψ at the crown {col['psi'][-1] / 1e6:.3f} MPa; cavitation margin to an "
          f"illustrative {I['p_cav'] / 1e6} MPa: {ch16.cavitation_margin(col['p'][-1], p_cav=I['p_cav']) / 1e6:.3f} MPa")
    hcap = ch16.capillary_rise(0.0728, np.pi / 2, 1000.0, a)
    print(f"  capillary rise in a {a * 1e6:.0f} µm conduit (ch01, full wetting, ours): {hcap:.3f} m — far below the crown")
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.0))
    ax[0].plot(col["p"] / 1e6, z, color=C_P, label="pressure p (tension when < 0)")
    ax[0].plot(col["psi"] / 1e6, z, color=C_U, label="water potential ψ = p − RTc + ρgz")
    ax[0].plot((col["psi"][0] - col["gravity"]) / 1e6, z, color=C_REF, ls=":", label="hydrostatic only (no flow)")
    ax[0].axvline(I["p_cav"] / 1e6, color=C_E, ls="--", label="illustrative cavitation limit (ours)")
    ax[0].axhline(hcap, color=C_WALL, lw=0.8, label=f"capillary rise {hcap:.2f} m")
    ax[0].set(xlabel="[MPa]", ylabel="height z [m]", title="a transpiring xylem column (ours)")
    ax[0].legend(fontsize=8)
    # phloem: sieve plates
    Qs = I["u_phloem"] * np.pi * I["a_sieve"] ** 2
    Ns = np.arange(1, 41)
    shares = []
    for N in Ns:
        r_pore = 0.6 * I["a_sieve"] / np.sqrt(N)  # pores scaled to fit (r ≤ a/√N), ours
        shares.append(ch16.sieve_tube_pressure_drop(Qs, I["mu_phloem"], I["a_sieve"], r_pore, I["L_sieve"],
                                                    I["ell_plate"], N)["plate_share"])
    st = ch16.sieve_tube_pressure_drop(Qs, I["mu_phloem"], I["a_sieve"], I["r_pore"], I["L_sieve"], I["ell_plate"],
                                       I["N_pores"])
    print(f"phloem (ours): a = {I['a_sieve'] * 1e6:.0f} µm, L = {I['L_sieve'] * 1e6:.0f} µm, N = {I['N_pores']:.0f} pores of "
          f"{I['r_pore'] * 1e6} µm: Δp per element (16.222) = {st['total']:.3f} Pa (tube {st['tube']:.3f}, plate viscous "
          f"{st['plate_viscous']:.3f}, plate entry {st['plate_entry']:.3f}); plate share {st['plate_share']:.1%}")
    chk = ch16.sieve_tube_pressure_drop(1e-15, 1e-3, 1e-5, 1e-5, 1e-3, 0.0, 1)
    print(f"  (16.222) with ℓ = 0, N = 1, r = a: {chk['total']:.6e} Pa = Poiseuille {8e-3 * 1e-15 * 1e-3 / (np.pi * 1e-20):.6e} Pa")
    ax[1].plot(Ns, shares, "o-", color=C_WALL, ms=3)
    ax[1].set(xlabel="number of pores N (pore radius 0.6a/√N, ours)", ylabel="share of Δp across the sieve plate",
              title="phloem: sieve plates take much of the pressure drop (16.222)")
    save(fig, out, "plants")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
