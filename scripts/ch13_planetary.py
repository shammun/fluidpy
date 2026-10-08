"""Chapter 13, §13.1–§13.4 — the planet in numbers: f, β, the inertial period, how good the β-plane is, and the size of
every term of the thin-shell equations for synoptic scales.  Also the two-convention lapse-rate table and the
stratification of an idealised ocean and of the standard atmosphere (§13.2).

Run: ``.venv/Scripts/python.exe scripts/ch13_planetary.py --no-show``   Figures -> outputs/ch13/planetary_*.png
"""
from __future__ import annotations

import numpy as np
from ch13_common import COLORS, DAY, finish, hemispheres, parse_args, save, setup

from fluidpy import ch13_geophysical_fluid_dynamics as ch13
from fluidpy.core.stratification import adiabatic_lapse_rate


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch13.illustrative_inputs()
    print(f"Omega (sidereal) = {ch13.OMEGA_EARTH:.6e} rad/s; one turn per solar day = {ch13.OMEGA_SOLAR_DAY:.6e} rad/s "
          f"({100 * (1 - ch13.OMEGA_SOLAR_DAY / ch13.OMEGA_EARTH):.3f} % smaller); sidereal day = {ch13.SIDEREAL_DAY:.1f} s")
    for lat, label in hemispheres(I["lat"]):
        f, b, h = ch13.coriolis_parameter(lat), ch13.beta_parameter(lat), ch13.hemisphere(lat)
        print(f"  {label}: f = {f:+.4e} 1/s, beta = {b:.4e} 1/(m s), inertial period = {ch13.inertial_period(f) / 3600:.2f} h, "
              f"deflection to the {h['turns']}, cyclones turn {h['cyclonic']}")
    y = np.linspace(-3.0e6, 3.0e6, 121)
    err = ch13.beta_plane_error(y, I["lat"])
    for yy in (5.0e5, 1.0e6, 2.0e6):
        print(f"  beta-plane error at y = {yy / 1e3:.0f} km: {ch13.beta_plane_error(yy, I['lat']):+.3e} "
              f"(leading term (y/R)^2/2 = {0.5 * (yy / ch13.EARTH_RADIUS_MEAN) ** 2:.3e})")
    tab = ch13.term_table_thin_layer(I["U_syn"], I["L_syn"], I["atm_H"], I["lat"], nu_H=1.0e5, nu_v=I["nu_v_atm"])
    tl = ch13.thin_layer_terms(I["U_syn"], I["L_syn"], I["atm_H"], I["lat"])
    print(f"thin-shell scales: U = {I['U_syn']} m/s, L = {I['L_syn'] / 1e3:.0f} km, H = {I['atm_H'] / 1e3:.0f} km -> "
          f"W = {tl['W']:.3f} m/s, H/L = {tl['aspect']:.2e}, Ro = {tl['Ro']:.3f}")
    print(tab.to_string(float_format=lambda v: f"{v:.3e}"))
    Ga = adiabatic_lapse_rate()
    print("lapse-rate table for the standard troposphere (both conventions; Gamma_a passed explicitly):")
    print(ch13.lapse_rate_table(-6.5e-3, Gamma_a=Ga).to_string())
    z = np.linspace(0.0, 5.0e4, 501)
    atm = ch13.atmosphere_layers(z)
    zo = np.linspace(-I["ocean_H"], 0.0, 841)
    oc = ch13.ocean_profile_idealized(zo)
    print(f"standard atmosphere: N = {np.sqrt(ch13.atmosphere_layers(5.0e3)['N2']):.4f} rad/s at 5 km (troposphere), "
          f"{np.sqrt(ch13.atmosphere_layers(2.5e4)['N2']):.4f} at 25 km (stratosphere)")
    print(f"idealised ocean: N_max = {oc['N'].max():.4f} rad/s just below the mixed layer -> buoyancy period "
          f"{2 * np.pi / oc['N'].max() / 60:.1f} min; N² from differences of rho_theta agrees with the analytic N² to "
          f"{np.max(np.abs(ch13.buoyancy_frequency_sq(zo, oc['rho_theta'], 1027.0) - oc['N2'])[2:-12]):.1e} 1/s²")

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    lat = np.deg2rad(np.linspace(-90, 90, 181))
    ax[0].plot(np.rad2deg(lat), ch13.coriolis_parameter(lat) * 1e4, color=COLORS["teal"], label="f [10⁻⁴ 1/s]")
    ax[0].plot(np.rad2deg(lat), ch13.beta_parameter(lat) * 1e11, color=COLORS["amber"], label="β [10⁻¹¹ 1/(m s)]")
    ax[0].axhline(0, color=COLORS["muted"], lw=0.8)
    ax[0].set(xlabel="latitude [°]", title="Coriolis parameter and its gradient")
    ax[0].legend()
    ax[1].plot(y / 1e3, 100 * err, color=COLORS["accent"])
    ax[1].set(xlabel="y [km] from 35° N", ylabel="(f₀ + βy − f)/f [%]", title="error of the β-plane")
    names = ["acceleration", "coriolis", "coriolis_w", "pressure", "friction_H", "friction_v"]
    ax[2].barh(names, [tab.loc[n, "x"] for n in names], color=[COLORS[c] for c in ("accent", "teal", "muted", "orange", "rose", "rose")])
    ax[2].set(xscale="log", xlabel="size of the term in the x-equation [m/s²]", title=f"synoptic scales: Ro = {tl['Ro']:.2f}")
    save(fig, out, "planetary_parameters")

    fig2, ax2 = plt.subplots(1, 2, figsize=(10, 4.4))
    ax2[0].plot(atm["T"] - 273.15, z / 1e3, color=COLORS["blue"])
    a2 = ax2[0].twiny()
    a2.plot(np.sqrt(atm["N2"]) * 1e2, z / 1e3, color=COLORS["amber"], ls="--")
    a2.set_xlabel("N [10⁻² rad/s] (dashed)")
    ax2[0].axhline(11, color=COLORS["muted"], lw=0.8)
    ax2[0].text(-55, 12, "tropopause", fontsize=8)
    ax2[0].set(xlabel="T [°C] (standard atmosphere)", ylabel="height [km]")
    ax2[1].plot(oc["T"], zo, color=COLORS["blue"], label="T [°C]")
    a3 = ax2[1].twiny()
    a3.plot(oc["N"] * 1e3, zo, color=COLORS["amber"], ls="--")
    a3.set_xlabel("N [10⁻³ rad/s] (dashed)")
    ax2[1].set(xlabel="T [°C] (idealised ocean, ours)", ylabel="z [m]", ylim=(-2000, 0))
    save(fig2, out, "planetary_stratification")
    print(f"inertial period at 35°: {ch13.inertial_period(ch13.coriolis_parameter(I['lat'])) / DAY:.3f} days")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
