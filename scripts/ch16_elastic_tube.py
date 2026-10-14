"""Chapter 16, §16.3 — a distensible (Hookean, thin-walled) tube: hoop stress (16.39), the tube law (16.41) and its
singularity (trap T16), steady Poiseuille flow through it (16.42)–(16.46) with the a₀⁴ of (16.44)–(16.45) restored
(slip S4, checked with pint), the falling effective resistance (Fig. 16.23's mechanism) and the compliance dV/dp.

Run: ``.venv/Scripts/python.exe scripts/ch16_elastic_tube.py --no-show``   Figures -> outputs/ch16/elastic_tube_*.png
"""
from __future__ import annotations

import numpy as np
from ch16_common import C_P, C_REF, C_U, C_WALL, finish, parse_args, save, setup
from scipy.integrate import solve_ivp

from fluidpy import ch16_biofluid_mechanics as ch16


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch16.illustrative_inputs()
    mu = I["mu"]
    a0, h, L = 1.5e-3, 0.15e-3, 0.1            # a soft illustrative vessel (ours)
    E = 0.2e6
    Eh_a0 = E * h / a0
    print(f"soft tube (ours): a₀ = {a0 * 1e3} mm, h = {h * 1e3} mm, E = {E / 1e6} MPa -> Eh/a₀ = {Eh_a0:.0f} Pa "
          f"({ch16.Pa_to_mmHg(Eh_a0):.1f} mmHg): the tube law (16.41) blows up as p − p_e → Eh/a₀")
    p = np.linspace(0, 0.9 * Eh_a0, 200)
    a = ch16.tube_radius_hookean(p, 0.0, a0, E, h)
    print(f"  p − p_e = 0.5Eh/a₀ -> a/a₀ = {ch16.tube_radius_hookean(0.5 * Eh_a0, 0.0, a0, E, h) / a0:.3f}; hoop stress there "
          f"{ch16.hoop_stress(0.5 * Eh_a0, 0.0, ch16.tube_radius_hookean(0.5 * Eh_a0, 0.0, a0, E, h), h) / 1e3:.1f} kPa")
    try:
        ch16.tube_radius_hookean(1.01 * Eh_a0, 0.0, a0, E, h)
    except ValueError as err:
        print(f"  beyond the singularity: ValueError ({str(err)[:60]}…)")
    # pint: the printed (16.45) prefactor has the wrong units (S4)
    Q_ = ch16.Q_
    printed = (8 * Q_(1.0, "Pa*s") / Q_(1.0, "m")).to_base_units()
    correct = (8 * Q_(1.0, "Pa*s") / Q_(1.0, "m") ** 4).to_base_units()
    need = (Q_(1.0, "Pa") / (Q_(1.0, "m^3/s") * Q_(1.0, "m"))).to_base_units()
    print(f"  pint (S4): dp/(Q dx) needs {need.units}; printed 8μ/(πa₀) has {printed.units}; 8μ/(πa₀⁴) has {correct.units}")
    # flow against Δp, rigid line
    dps = np.linspace(0, 0.8 * Eh_a0, 60)
    p2 = 0.0
    Qe = ch16.elastic_tube_flow_rate(p2 + dps, p2, 0.0, a0, E, h, mu, L)
    Qr = ch16.poiseuille_flow_rate(a0, dps, L, mu)
    Q1 = ch16.elastic_tube_flow_rate(0.5 * Eh_a0, 0.0, 0.0, a0, E, h, mu, L)
    sol = solve_ivp(lambda x, y: [-8 * mu * Q1 / (np.pi * ch16.tube_radius_hookean(y[0], 0.0, a0, E, h) ** 4)],
                    (0, L), [0.5 * Eh_a0], rtol=1e-12, atol=1e-9, dense_output=True)
    xq = np.linspace(0, L, 6)
    err = np.max(np.abs(sol.sol(xq)[0] - ch16.elastic_tube_pressure(xq, Q1, 0.5 * Eh_a0, 0.0, a0, E, h, mu)))
    print(f"  Q (16.46) at Δp = 0.5Eh/a₀: {Q1:.5e} m³/s = {Q1 / ch16.poiseuille_flow_rate(a0, 0.5 * Eh_a0, L, mu):.3f}× rigid; "
          f"p(x) closed form vs solve_ivp of (16.43): max |Δp| = {err:.2e} Pa; outlet {sol.y[0, -1]:.3e} Pa")
    print(f"  rigid limit E = 1e12: Q/Q_Pois − 1 = {ch16.elastic_tube_flow_rate(1000.0, 0.0, 0.0, a0, 1e12, h, mu, L) / ch16.poiseuille_flow_rate(a0, 1000.0, L, mu) - 1:.2e}")
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
    ax[0].plot(p / Eh_a0, a / a0, color=C_WALL)
    ax[0].set(xlabel="(p − p_e)/(Eh/a₀)", ylabel="a/a₀", title="tube law (16.41): a → ∞ as a₀(p − p_e)/Eh → 1")
    ax[1].plot(dps, Qe * 1e6, color=C_U, label="elastic tube (16.46)")
    ax[1].plot(dps, Qr * 1e6, color=C_REF, ls=":", label="rigid tube (16.17)")
    ax[1].set(xlabel="Δp = p₁ − p₂ [Pa] (p₂ = p_e = 0)", ylabel="Q [cm³/s]", title="flow is nonlinear in Δp")
    ax[1].legend()
    pm = np.linspace(0, 0.75 * Eh_a0, 60)
    dpf = 0.05 * Eh_a0
    Reff = ch16.elastic_tube_resistance(pm + dpf / 2, pm - dpf / 2, 0.0, a0, E, h, mu, L)
    R0 = ch16.vascular_resistance(mu, L, a0)
    ax[2].plot(pm / ch16.MMHG, Reff / R0, color=C_P)
    ax[2].set(xlabel="mean pressure [mmHg]", ylabel="R_eff / R(a₀)",
              title="a distensible vessel: resistance falls with pressure (ours)")
    save(fig, out, "elastic_tube")
    K0 = ch16.compliance_from_tube(0.0, 0.0, a0, E, h, L)
    print(f"  compliance dV/dp (ours) at p = p_e: {K0:.4e} m³/Pa = {K0 * 1e6 * ch16.MMHG:.4f} mL/mmHg; "
          f"at 0.5Eh/a₀: {ch16.compliance_from_tube(0.5 * Eh_a0, 0.0, a0, E, h, L) / K0:.1f}× larger")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
