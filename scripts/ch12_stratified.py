"""§12.11 (C14, C15): turbulence in a stratified medium — flux and gradient Richardson numbers (12.107)-(12.109), the gradient
Richardson number from an IN-SITU temperature gradient with the adiabatic lapse rate in both sign conventions, the
Monin-Obukhov length (12.110) and Rf = z/L_M (12.111), the log-linear wind profile for stable / neutral / unstable surface
layers (Figs. 12.21-12.22 analogues), mixing of an unstable column (Fig. 12.10 analogue) and the temperature spectrum with its
-5/3 and -1 ranges (12.113)-(12.114) (Fig. 12.23 analogue).

Lapse-rate convention: computed with Kundu's Gamma = dT/dz (Gamma_a = -g/Cp ~ -9.8 K/km); the meteorological
Gamma = -dT/dz (Gamma_d ~ +9.8 K/km) is printed alongside every verdict.  "Temperature" in the buoyancy terms is potential
temperature.  Heat fluxes are positive upward.

Run: ``.venv/Scripts/python.exe scripts/ch12_stratified.py --no-show``
Figures -> outputs/ch12/c14_richardson.png, c15_surface_layer.png.
Our surface layer: u* = 0.3 m/s, z0 = 3 cm, T = 290 K, rho = 1.2 kg/m3, kappa = 0.40 (the usual atmospheric value), beta = 5.
"""
from __future__ import annotations

import numpy as np

from ch12_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch12_turbulence as ch12
from fluidpy.core import wall_turbulence as WT
from fluidpy.core.stratification import adiabatic_lapse_rate, lapse_rate_convention
from fluidpy.core.thermo import CP_AIR


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    kap, T, rho, cp, us, z0 = 0.40, 290.0, 1.2, CP_AIR, 0.3, 0.03
    Ga = adiabatic_lapse_rate()
    print(f"adiabatic gradient: Kundu Gamma_a = {Ga * 1e3:.2f} K/km; meteorology Gamma_d = {lapse_rate_convention(Ga, 'meteorology') * 1e3:+.2f} K/km")
    print("gradient Richardson number from an in-situ gradient, dU/dz = 0.01 1/s, alpha = 1/288 K:")
    for name, dTdz in (("standard atmosphere", -6.5e-3), ("isothermal layer", 0.0), ("adiabatic layer", Ga), ("superadiabatic", -12e-3),
                       ("night inversion", +20e-3)):
        k = ch12.gradient_richardson_thermal(dTdz, 0.01, 1 / 288.0, Gamma_a=Ga)
        wrong = ch12.gradient_richardson_thermal(dTdz, 0.01, 1 / 288.0)["Ri"]
        print(f"  {name:20s} dT/dz = {dTdz * 1e3:+6.1f} K/km: Ri = {k['Ri']:+8.3f} ({k['verdict']});  {k['verdict_kundu']}  |  {k['verdict_met']}"
              f"   [without Gamma_a one would get {wrong:+.3f}]")
    th = ch12.gradient_richardson_thermal(-6.5e-3 - Ga, 0.01, 1 / 288.0)
    print(f"  same layer from d(theta)/dz = {(-6.5e-3 - Ga) * 1e3:.2f} K/km with Gamma_a = 0: Ri = {th['Ri']:.3f} (identical)")
    print("flux Richardson number (12.107), -uw = 0.09 m2/s2, dU/dz = 0.075 1/s:")
    for wT in (0.10, 0.0, -0.005, -0.02):
        Rf = ch12.flux_richardson(wT, -0.09, 0.075, 1 / T)
        print(f"  mean(wT') = {wT:+.3f} K m/s (H = {ch12.turbulent_heat_flux(1.0, 1.0, wT, rho, cp):+7.1f} W/m2): Rf = {Rf:+.3f} -> {ch12.turbulence_regime(Rf)}"
              f"; Ri at Pr_T = 1, 2: {Rf * 1:+.3f}, {Rf * 2:+.3f}")
    print(f"Rf_cr = 0.25 corresponds to Ri = Pr_T/4: Rf(Ri = 0.5, Pr_T = 2) = {ch12.flux_from_gradient_richardson(0.5, 2.0):.3f}")
    print("surface layer at z = 10 m (book's log-linear profile; unstable side also with Businger-Dyer):")
    for label, H in (("sunny afternoon", 200.0), ("weak heating", 20.0), ("neutral overcast", 0.0), ("clear night", -30.0)):
        s = ch12.surface_layer_state(us, H, T, z0, 10.0, rho, cp, kap)
        bd = ch12.surface_layer_state(us, H, T, z0, 10.0, rho, cp, kap, unstable="businger_dyer")
        print(f"  {label:17s} H = {H:+6.1f} W/m2: wT = {s['wT']:+.4f} K m/s, L_M = {s['L_M']:+9.1f} m, Rf = z/L_M = {s['Rf']:+.3f} ({s['regime']}, {s['layer']}),"
              f" U = {s['U']:.2f} (log-linear), {bd['U']:.2f} (B-D), neutral {s['U_neutral']:.2f} m/s; z(Rf = 1/4) = {s['z_crit']:.1f} m")
        if bd["verdict_kundu"]:
            print(f"  {'':17s} layer at 10 m: {bd['verdict']};  {bd['verdict_kundu']}  |  {bd['verdict_met']}")
    L = ch12.monin_obukhov_from_fluxes(rho * us ** 2, 200.0, rho, cp, T, kappa=kap)
    print(f"L_M from tau = rho u*^2 and H = 200 W/m2: {L:.2f} m; regimes at z = 1, 10, 100 m: {ch12.surface_layer_regime(np.array([1.0, 10.0, 100.0]), L)}")
    z = np.linspace(0.0, 100.0, 201)
    mix = ch12.mixing_potential_energy_change(z, 300.0 - 0.01 * z, 1 / 300.0, 1.2)
    mix2 = ch12.mixing_potential_energy_change(z, 300.0 + 0.01 * z, 1 / 300.0, 1.2)
    print(f"mixing a 100 m column (potential temperature): unstable dtheta/dz = -10 K/km -> dPE = {mix['dPE']:+.1f} J/m2 (released); "
          f"stable +10 K/km -> {mix2['dPE']:+.1f} J/m2 (must be supplied)")
    print(f"temperature-spectrum exponents by the Pi theorem: {ch12.scalar_spectrum_exponents()}")
    nu_w, kap_w, eps, epsT = 1.0e-6, 1.4e-7, 1.0e-6, 1.0e-7
    eta_w = ch12.kolmogorov_scales(nu_w, eps)[0]
    print(f"sea water (our numbers): eta = {eta_w * 1e3:.2f} mm, Batchelor scale = {ch12.batchelor_scale(nu_w, kap_w, eps) * 1e3:.3f} mm "
          f"(= eta/sqrt(Pr), Pr = {nu_w / kap_w:.1f})")

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    zeta = np.linspace(-1.0, 1.0, 401)
    ax[0].plot(zeta, zeta, color=COLORS["accent"], label=r"Rf = $z/L_M$ (12.111)")
    ax[0].plot(zeta, ch12.gradient_richardson_surface_layer(zeta, 1.0, consistent=True), color=COLORS["teal"], ls="--",
               label=r"with the shear of the log-linear profile, $\zeta/(1+5\zeta)$")
    ax[0].axhline(0.25, color=COLORS["rose"], ls=":", label=r"Rf$_{cr}\approx 0.25$ (observed)")
    ax[0].axhspan(-1, 0, color=COLORS["orange"], alpha=0.08)
    ax[0].set(ylim=(-1, 1), xlabel=r"$z/L_M$", ylabel="Rf", title="unstable (shaded) | stable")
    ax[0].legend(fontsize=8)
    K = np.geomspace(1.0, 3e4, 300)
    ax[1].loglog(K, ch12.scalar_spectrum(K, eps, epsT, nu_w, kap_w, cutoff=True), color=COLORS["blue"], label=r"$S_T$: $-5/3$ then $-1$")
    ax[1].loglog(K, ch12.model_spectrum(K, eps, nu_w), color=COLORS["teal"], label="velocity E(K) (model)")
    ax[1].axvline(1 / eta_w, color=COLORS["grid"])
    ax[1].axvline(1 / ch12.batchelor_scale(nu_w, kap_w, eps), color=COLORS["grid"], ls="--")
    ax[1].set(ylim=(1e-16, 1e-4), xlabel="K [rad/m]", ylabel=r"$S_T$ [K$^2$ m], E [m$^3$/s$^2$]", title="temperature and velocity spectra in water")
    ax[1].legend(fontsize=8)
    save(fig, out, "c14_richardson")

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    zz = np.geomspace(z0 * 1.5, 100.0, 200)
    for LM, c, lab in ((30.0, COLORS["blue"], r"stable, $L_M$ = +30 m"), (np.inf, "k", "neutral"), (-30.0, COLORS["orange"], r"unstable, $L_M$ = $-$30 m")):
        ax[0].semilogy(WT.surface_layer_wind(zz, us, z0, LM, kappa=kap, unstable="businger_dyer"), zz, color=c, label=lab)
    ax[0].semilogy(WT.surface_layer_wind(zz[zz < 6], us, z0, -30.0, kappa=kap), zz[zz < 6], color=COLORS["orange"], ls=":",
                   label=r"book's log-linear form, valid for $z < -L_M/5$")
    ax[0].axhline(30.0, color=COLORS["grid"])
    ax[0].set(xlabel="U [m/s]", ylabel="z [m]", title="wind in the surface layer: ln z against U")
    ax[0].legend(fontsize=8)
    Hs = np.linspace(-60.0, 300.0, 181)
    Ls = np.array([ch12.monin_obukhov_from_fluxes(rho * us ** 2, H, rho, cp, T, kappa=kap) for H in Hs])
    ax[1].plot(Hs, 10.0 / Ls, color=COLORS["accent"])
    ax[1].axhline(0.25, color=COLORS["rose"], ls=":", label="Rf = 0.25")
    ax[1].axhline(0, color=COLORS["grid"])
    ax[1].set(xlabel=r"surface heat flux H [W/m$^2$] (upward positive)", ylabel=r"Rf at 10 m = $10/L_M$",
              title=r"$u_*$ = 0.3 m/s: heating makes Rf negative, cooling positive")
    ax[1].legend(fontsize=8)
    save(fig, out, "c15_surface_layer")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
