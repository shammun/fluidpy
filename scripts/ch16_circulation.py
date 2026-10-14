"""Chapter 16, §16.2–§16.3 (steady) — the circulation in numbers: the ventricle's pressure–volume loop and its work,
SV, EF, CO, CI (16.1)–(16.4), peripheral resistance (16.9)–(16.11), Poiseuille's resistance R = 8μL/(πa⁴) (16.20),
continuity across the vessel tree (Fig. 16.10, ours) and the Reynolds/Womersley scales of each vessel class.

Run: ``.venv/Scripts/python.exe scripts/ch16_circulation.py --no-show``   Figures -> outputs/ch16/circulation_*.png
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
    mm = ch16.MMHG
    # --- the ventricle: PV loop, SV, EF, CO, CI ----------------------------------------------------------------------
    kw = dict(p_min=I["p_min"], p_ED=I["p_ED"], p_open=I["p_open"], p_peak=I["p_peak"], p_ES=I["p_ES"])
    loop = ch16.pv_loop_schematic(I["EDV"], I["ESV"], **kw)
    W = ch16.pv_loop_work(loop["p"], loop["V"])
    SV = ch16.stroke_volume(I["EDV"], I["ESV"])
    EF = ch16.ejection_fraction(SV, I["EDV"])
    CO = ch16.cardiac_output(SV, I["HR"])
    print(f"PV loop (ours): EDV {I['EDV'] * 1e6:.0f} mL, ESV {I['ESV'] * 1e6:.0f} mL -> SV (16.1) = {SV * 1e6:.0f} mL, "
          f"EF (16.2) = {EF:.3f}; HR {I['HR']:.0f}/min -> CO (16.3) = {CO * 6e4:.2f} L/min, "
          f"CI (16.4) = {ch16.cardiac_index(CO, I['BSA']) * 6e4:.2f} L/(min m²)")
    print(f"  work per beat W = -∮p dV = {W['work']:.4f} J ({W['orientation']}); exact for these legs "
          f"{ch16.pv_loop_schematic_work(I['EDV'], I['ESV'], **kw):.4f} J; mean power {ch16.cardiac_power(W['work'], I['HR']):.3f} W")
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    ax[0].fill(loop["V"] * 1e6, loop["p"] / mm, color=C_E, alpha=0.25, label=f"area = work = {W['work']:.2f} J")
    ax[0].plot(loop["V"] * 1e6, loop["p"] / mm, color=C_P)
    for lab, (v, p) in loop["corners"].items():
        ax[0].annotate(lab, (v * 1e6, p / mm), textcoords="offset points", xytext=(6, 4), fontsize=11)
    ax[0].set(xlabel="left-ventricular volume V [mL]", ylabel="pressure p [mmHg]",
              title="PV loop (ours, illustrative): counter-clockwise, W = −∮p dV > 0")
    ax[0].legend(loc="upper center")
    # --- peripheral resistance, the 1/3 rule ----------------------------------------------------------------------------
    kwk = dict(SV=I["SV"], HR=I["HR"], systole_share=I["systole_share"], R=I["R_periph"], K=I["K_compliance"])
    s = ch16.windkessel_summary(**kwk)
    print(f"peripheral resistance R = {I['R_periph']:.3g} Pa s/m³ = {I['R_periph'] / ch16.PRU:.3f} PRU; mean flow "
          f"{s['Q_mean'] * 1e6:.1f} cm³/s -> mean pressure p̄_A = QR (16.9)-(16.10) = {s['p_mean'] / mm:.1f} mmHg")
    print(f"  periodic Windkessel beat: p_S {s['pS'] / mm:.1f}, p_D {s['pD'] / mm:.1f} mmHg; ⅓ rule (16.11) "
          f"{s['map_third'] / mm:.1f} mmHg vs true mean {s['p_mean'] / mm:.1f} mmHg (error {s['map_error'] / mm:+.1f} mmHg)")
    # --- Poiseuille resistance and networks (C04) ------------------------------------------------------------------------
    mu, L = I["mu"], I["L_vessel"]
    a = np.geomspace(5e-6, 5e-3, 200)
    R = ch16.vascular_resistance(mu, L, a)
    print(f"R = 8μL/(πa⁴) (16.20): a = 1 mm, L = 10 cm -> {ch16.vascular_resistance(mu, L, 1e-3):.4g} Pa s/m³; halving a "
          f"multiplies R by {ch16.vascular_resistance(mu, L, 0.5e-3) / ch16.vascular_resistance(mu, L, 1e-3):.1f}")
    Ra = ch16.vascular_resistance(mu, 1e-3, I["a_arteriole"])
    print(f"  one arteriole (a = {I['a_arteriole'] * 1e6:.0f} µm, 1 mm long): {Ra:.3g} Pa s/m³; 10 in parallel: "
          f"{ch16.network_resistance([Ra] * 10, arrangement='parallel'):.3g} (= R/10); 10⁶ in parallel: {Ra / 1e6:.3g}; "
          f"3 in series: {ch16.network_resistance([Ra] * 3, arrangement='series'):.3g} Pa s/m³")
    Q_test = ch16.poiseuille_flow_rate(1e-3, 100.0, L, mu)
    print(f"  check: Δp/Q from (16.17) = {100.0 / Q_test:.6g} = R (16.20) {ch16.vascular_resistance(mu, L, 1e-3):.6g}; "
          f"wall stress (16.19, + sign) {ch16.poiseuille_wall_stress(1e-3, 100.0, L):.4f} Pa")
    ax[1].loglog(a * 1e3, R, color=C_WALL, label="R = 8μL/(πa⁴), L = 10 cm")
    ax[1].loglog(a * 1e3, R[0] * (a / a[0]) ** -4, color=C_REF, ls=":", label="slope −4")
    ax[1].set(xlabel="radius a [mm]", ylabel="resistance R [Pa s/m³]", title="(16.20): halve the radius, ×16 resistance")
    ax[1].legend()
    save(fig, out, "circulation_pv_loop_resistance")
    # --- continuity across the tree (Fig. 16.10, ours) and the scales of each class --------------------------------
    tab = ch16.vessel_scales_table(HR=I["HR"], rho=I["rho"], mu=I["mu"])
    Qco = s["Q_mean"]
    tab["A_total_cm2"] = Qco / tab["u_mean"] * 1e4
    print("vessel classes (ours, illustrative; total area from continuity A = Q/ū):")
    print(tab[["vessel", "a", "u_mean", "A_total_cm2", "Re", "alpha", "regime"]].to_string(
        index=False, float_format=lambda v: f"{v:.4g}"))
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.4))
    xs = np.arange(len(tab))
    ax[0].bar(xs, tab["A_total_cm2"], color=C_WALL, alpha=0.7, label="total area Q/ū [cm²]")
    ax0b = ax[0].twinx()
    ax0b.plot(xs, tab["u_mean"] * 100, "o-", color=C_U, label="mean speed ū [cm/s]")
    ax[0].set_xticks(xs, tab["vessel"], rotation=30)
    ax[0].set(ylabel="total cross-section [cm²]", title="continuity Q = ΣA·ū (ours, illustrative inputs)")
    ax0b.set_ylabel("mean speed [cm/s]")
    ax[1].loglog(tab["Re"], tab["alpha"], "o", color=C_P)
    for _, row in tab.iterrows():
        ax[1].annotate(row["vessel"], (row["Re"], row["alpha"]), textcoords="offset points", xytext=(5, 3), fontsize=9)
    ax[1].axhline(1.0, color=C_REF, ls=":")
    ax[1].set(xlabel="Reynolds number Re = ρū2a/μ", ylabel="Womersley number α = a√(ω/ν) (16.12)",
              title="vessel scales (ours)")
    save(fig, out, "circulation_vessel_scales")
    # --- pressure along resistances in series (Fig. 16.2's reading, ours) ---------------------------------------------
    Rk = np.array([0.03, 0.07, 0.45, 0.25, 0.08, 0.07, 0.05]) * I["R_periph"]
    pk = ch16.series_pressure_drop(Rk, Qco, s["p_mean"])
    print("mean pressure after each class (ours, illustrative resistance shares):",
          ", ".join(f"{v / mm:.0f}" for v in pk), "mmHg")
    print(f"entrance length (16.37): Re = 500, d = 4 mm -> laminar ℓ = {ch16.entrance_length(500.0, 4e-3, regime='laminar') * 100:.1f} cm")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
