"""Chapter 1, §1.10: static atmospheres and oceans — the lapse-rate conventions, the U.S. Standard Atmosphere built from T(z), the isothermal
model and scale height, a synthetic boundary layer (Fig. 1.9 idea) with θ, ρ_θ and N², checks of (1.32) and (1.34),
and the ocean criterion (1.35).

Run: ``.venv/Scripts/python.exe scripts/ch01_atmosphere_profiles.py --no-show``
Figures → outputs/ch01/fig_standard_atmosphere.png, fig1_9_boundary_layer.png, fig_ocean_stability.png,
fig_lapse_rate_conventions.png (C54: the stability criterion in both lapse-rate conventions).
Profiles are synthetic (ours) except USSA-1976.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch01"), help="figure folder")
    ap.add_argument("--no-show", action="store_true", help="do not open a window")
    args = ap.parse_args()
    t_start = time.perf_counter()
    import matplotlib

    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from fluidpy import ch01_introduction as ch01
    from fluidpy.core.style import COLORS, use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    g, R = ch01.G0, ch01.R_AIR

    # --- USSA-1976 vs isothermal model -----------------------------------------------------------------------------------------
    z = np.linspace(0.0, 84000.0, 421)
    T, p, rho = ch01.standard_atmosphere(z)
    p_ode, _, _ = ch01.atmosphere_from_temperature(z[z <= 11000.0], lambda zz: 288.15 - 6.5e-3 * zz, ch01.P_ATM)
    fig, axs = plt.subplots(1, 3, figsize=(12, 4.5), sharey=True)
    axs[0].plot(T, z / 1e3, color=COLORS["accent"])
    axs[0].set_xlabel("T [K]")
    axs[0].set_ylabel("geopotential altitude [km]")
    axs[0].set_title("USSA-1976 temperature")
    axs[1].semilogx(p, z / 1e3, color=COLORS["accent"], label="USSA-1976")
    H = ch01.scale_height(250.0)
    axs[1].semilogx(ch01.isothermal_pressure(z, ch01.P_ATM, 250.0), z / 1e3, color=COLORS["orange"], ls="--",
                    label=f"isothermal 250 K (H = {H / 1e3:.2f} km)")
    for k in (1, 2, 3):
        axs[1].axhline(k * H / 1e3, color=COLORS["grid"], lw=0.8)
    axs[1].set_xlabel("p [Pa]")
    axs[1].legend(fontsize=8)
    axs[1].set_title("Pressure falls by e every H")
    axs[2].plot(ch01.perfect_gas_sound_speed(T), z / 1e3, color=COLORS["teal"])
    axs[2].set_xlabel("c = √(γRT) [m/s]")
    axs[2].set_title("Speed of sound (Eq. 1.27)")
    fig.savefig(out / "fig_standard_atmosphere.png", bbox_inches="tight")

    # --- synthetic boundary layer: T, θ, N² --------------------------------------------------------------------------------------
    zb = np.linspace(0.0, 2000.0, 401)
    col = ch01.synthetic_boundary_layer_column(zb)
    Ga = ch01.adiabatic_lapse_rate()
    fig, axs = plt.subplots(1, 3, figsize=(12, 4.5), sharey=True)
    axs[0].plot(col["T"] - 273.15, zb, color=COLORS["ink"], lw=2, label="environment T(z)")
    for T0 in np.arange(0.0, 40.0, 5.0):  # neutral (dry-adiabatic) reference lines
        axs[0].plot(T0 + Ga * zb, zb, color=COLORS["grid"], lw=0.8)
    axs[0].set_xlim(0, 20)
    axs[0].set_xlabel("T [°C]")
    axs[0].set_ylabel("z [m]")
    axs[0].set_title(f"(a) T(z); grey: Γ_a = {Ga * 1e3:.2f} K/km")
    axs[1].plot(col["theta"] - 273.15, zb, color=COLORS["accent"], lw=2)
    axs[1].set_xlabel("θ [°C]")
    axs[1].set_title("(b) potential temperature (1.31)")
    colors = {"stable": COLORS["teal"], "neutral": COLORS["amber"], "unstable": COLORS["rose"]}
    for lab in ("stable", "neutral", "unstable"):
        mask = col["stability"] == lab
        axs[2].scatter(col["N2"][mask] * 1e4, zb[mask], s=6, color=colors[lab], label=lab)
    axs[2].axvline(0.0, color=COLORS["muted"], lw=0.8)
    axs[2].set_xlabel("N² [10⁻⁴ s⁻²]")
    axs[2].set_title("(c) N² = (g/θ) dθ/dz")
    axs[2].legend(fontsize=8)
    fig.savefig(out / "fig1_9_boundary_layer.png", bbox_inches="tight")

    # --- C54: environment lapse rates vs the dry adiabat, the criterion in both conventions -----------------------------
    zl = np.linspace(0.0, 3000.0, 61)
    T_surf = 288.15
    fig, (ax, axt) = plt.subplots(1, 2, figsize=(14, 5.0), gridspec_kw={"width_ratios": [1.2, 1], "wspace": 0.25})
    ax.plot(ch01.parcel_temperature(T_surf, zl) - 273.15, zl, color=COLORS["ink"], lw=2.5, ls="--",
            label=f"parcel dry adiabat, dT/dz = Γ_a = {Ga * 1e3:.2f} K/km")
    lapse_rows = []
    for dTdz, col_ in ((-6.5e-3, COLORS["teal"]), (-12.0e-3, COLORS["rose"]), (5.0e-3, COLORS["blue"])):
        res_k = ch01.lapse_rate_stability(dTdz, prefix=False)  # Kundu: dT/dz vs Γ_a (negative)
        res_m = ch01.lapse_rate_stability(dTdz, convention="meteorology", prefix=False)  # −dT/dz vs +g/C_p
        assert res_k.code == res_m.code  # the verdict never depends on the convention
        verdict, kundu, met = res_k.verdict, res_k.text, res_m.text
        ax.plot(T_surf + dTdz * zl - 273.15, zl, color=col_, lw=2,
                label=f"{verdict}: dT/dz {kundu}  ⇔  Γ_met {met}")
        lapse_rows.append((dTdz, verdict, kundu, met))
    ax.set_xlabel("T [°C]")
    ax.set_ylabel("z [m]")
    ax.set_title("Environment T(z) vs a parcel rising from the ground")
    ax.legend(fontsize=7.5, loc="lower left")
    axt.axis("off")
    table = [["convention", "Γ ≡", "Γ_a (dry air)", "stable when"],
             ["Kundu (book, fluidpy)", "dT/dz", f"{ch01.lapse_rate_convention(Ga) * 1e3:+.2f} K/km".replace("-", "−"),
              "Γ > Γ_a"],
             ["meteorology", "−dT/dz", f"{ch01.lapse_rate_convention(Ga, 'meteorology') * 1e3:+.2f} K/km",
              "Γ < Γ_a"]]
    tab = axt.table(cellText=table[1:], colLabels=table[0], loc="center", cellLoc="center",
                    colWidths=[0.36, 0.18, 0.26, 0.2])
    tab.auto_set_font_size(False)
    tab.set_fontsize(9)
    tab.scale(1.0, 2.0)
    axt.set_title("Same physics, opposite sign and inequality:\nΓ_met = −Γ_Kundu")
    fig.savefig(out / "fig_lapse_rate_conventions.png", bbox_inches="tight")

    # checks of (1.32) and (1.34) by finite differences of the sampled column (interior points away from the kinks)
    dth_fd = ch01.derivative_2nd_order(col["theta"], zb)
    drt_fd = ch01.derivative_2nd_order(col["rho_theta"], zb)
    inner = (np.abs(zb - 800) > 10) & (np.abs(zb - 1000) > 10) & (zb > 10) & (zb < 1990)
    err132 = np.max(np.abs(dth_fd[inner] - col["dtheta_dz"][inner]))
    err134 = np.max(np.abs(-drt_fd[inner] / col["rho_theta"][inner] - dth_fd[inner] / col["theta"][inner]))
    invariant = np.max(np.abs(col["theta"] * col["rho_theta"] / (ch01.P_REF / R) - 1))

    # --- ocean: in-situ density looks stable, potential density near neutral -------------------------------------------------
    depth = np.linspace(0.0, 4000.0, 401)
    zo = -depth
    T_ocean = ch01.celsius_to_kelvin(2.0 + 16.0 * np.exp(-depth / 500.0))  # thermocline (synthetic)
    S_ocean = 35.0
    c = 1500.0
    rho_surf = ch01.seawater_density_linear(T_ocean, S_ocean)
    # isentropic compression by the hydrostatic pressure: dρ/dh = ρ g/c² at fixed T and S  →  ρ = ρ_s exp(g h/c²)
    rho_insitu = rho_surf * np.exp(g * depth / c ** 2)
    drho_dz = ch01.derivative_2nd_order(rho_insitu, zo)
    crit = ch01.ocean_potential_density_gradient(drho_dz, rho_insitu, c)
    N2_ocean = ch01.brunt_vaisala_sq(rho_insitu, drho_dz, ch01.isentropic_density_gradient(rho_insitu, c))
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)
    axs[0].plot(rho_insitu, depth, color=COLORS["blue"], label="in-situ ρ")
    axs[0].plot(rho_surf, depth, color=COLORS["accent"], ls="--", label="ρ brought to the surface")
    axs[0].invert_yaxis()
    axs[0].set_xlabel("ρ [kg/m³]")
    axs[0].set_ylabel("depth [m]")
    axs[0].legend(fontsize=8)
    axs[0].set_title("Compression adds density with depth")
    axs[1].semilogx(np.abs(drho_dz), depth, color=COLORS["blue"], label="|dρ/dz| (in situ)")
    axs[1].semilogx(np.abs(crit), depth, color=COLORS["teal"], label="|dρ/dz + ρg/c²| (1.35)")
    axs[1].semilogx(np.full_like(depth, 1027.0 * g / c ** 2), depth, color=COLORS["muted"], ls=":", label="ρg/c²")
    axs[1].set_xlabel("[kg/m⁴]")
    axs[1].legend(fontsize=8)
    axs[1].set_title("Deep water is far closer to neutral than ρ(z) suggests")
    fig.savefig(out / "fig_ocean_stability.png", bbox_inches="tight")

    i5, i11 = np.searchsorted(z, 5000.0), np.searchsorted(z, 11000.0)
    print(f"USSA p(5 km) = {p[i5]:.1f} Pa, p(11 km) = {p[i11]:.1f} Pa; ODE from T(z): {p_ode[i5]:.1f}, {p_ode[-1]:.1f} Pa")
    print(f"scale height at 250 K: {H / 1e3:.3f} km; p(H)/p0 = {ch01.isothermal_pressure(H, 1.0, 250.0):.5f} (e^-1 = {np.exp(-1):.5f})")
    print(f"dry adiabatic lapse rate: Kundu Gamma_a = dT_a/dz = {Ga * 1e3:.3f} K/km; meteorology Gamma_a = -dT_a/dz = "
          f"{ch01.lapse_rate_convention(Ga, 'meteorology') * 1e3:.3f} K/km; USSA troposphere dT/dz = -6.5 K/km "
          f"-> N^2 = {ch01.brunt_vaisala_sq_from_lapse(288.15, -6.5e-3):.3e} 1/s^2, "
          f"period {ch01.stability_timescale(ch01.brunt_vaisala_sq_from_lapse(288.15, -6.5e-3))[1] / 60:.1f} min")
    for dTdz, verdict, kundu, met in lapse_rows:
        print(f"  dT/dz = {dTdz * 1e3:+.1f} K/km: {verdict:8s} Kundu: {kundu}   meteorology: {met}   "
              f"(N^2 sign: {ch01.classify_stability(ch01.brunt_vaisala_sq_from_lapse(288.15, dTdz))})")
    print(f"  full texts: {ch01.lapse_rate_stability(-6.5e-3).text!r} | "
          f"{ch01.lapse_rate_stability(-6.5e-3, convention='meteorology').text!r}")
    print(f"theta of air at 500 hPa and 250 K: {ch01.potential_temperature(250.0, 5.0e4):.2f} K")
    print(f"boundary layer: N^2 mixed = {col['N2'][50]:.2e} ({col['stability'][50]}), inversion = {col['N2'][180]:.2e}, "
          f"upper = {col['N2'][300]:.2e} 1/s^2")
    print(f"(1.32) max |FD - formula| = {err132:.2e} K/m; (1.34) max residual = {err134:.2e} 1/m; "
          f"theta*rho_theta/(p_ref/R) - 1 max = {invariant:.1e}")
    print(f"ocean: rho g/c^2 = {1027.0 * g / c**2:.3e} kg/m^4; at 3000 m dρ/dz = {drho_dz[300]:.3e}, (1.35) = {crit[300]:.3e}, "
          f"N^2 = {N2_ocean[300]:.3e} 1/s^2 (sign of (1.35) opposite to N^2: {np.all(np.sign(crit[1:-1]) == -np.sign(N2_ocean[1:-1]))})")
    print(f"saved 4 figures in {out}  ({time.perf_counter() - t_start:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
