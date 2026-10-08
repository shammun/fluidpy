"""Chapter 13, §13.8, §13.10–§13.12 — the waves of rotating shallow water.

* the three roots of the cubic dispersion relation (two Poincaré, one Rossby) with the size of each term, for the
  external and the first baroclinic mode;
* Poincaré waves: frequency, group velocity, the current ellipse in both hemispheres; inertial circles;
* Kelvin waves: trapping scale, which way they can travel, residuals of the three equations; a Poincaré and a Kelvin
  wave propagated for one period by the C-grid model.

Run: ``.venv/Scripts/python.exe scripts/ch13_shallow_water_waves.py --no-show``   Figures -> outputs/ch13/waves_*.png
"""
from __future__ import annotations

import numpy as np
from ch13_common import COLORS, DAY, Timer, finish, hemispheres, parse_args, save, setup

from fluidpy import ch13_geophysical_fluid_dynamics as ch13


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch13.illustrative_inputs()
    f, beta, g = ch13.coriolis_parameter(I["lat"]), ch13.beta_parameter(I["lat"]), 9.80665
    c0 = ch13.long_wave_speed(I["ocean_H"])
    c1 = ch13.baroclinic_mode_speed(I["ocean_N"], I["ocean_H"], 1)
    for name, c in (("external", c0), ("first baroclinic", c1)):
        for lam in I["wavelengths"]:
            k = 2.0 * np.pi / lam
            wm, wr, wp = ch13.shallow_water_omega(k, 0.0, c, f, beta)
            t = ch13.dispersion_term_sizes(k, 0.0, c, f, beta, wr)
            big = max(("omega3", "gravity", "rotation", "beta"), key=lambda n: abs(t[n]))
            small = min(("omega3", "gravity", "rotation", "beta"), key=lambda n: abs(t[n]))
            print(f"{name:16s} c = {c:7.2f} m/s, wavelength {lam / 1e3:5.0f} km: omega/f = {wm / f:+.4f}, {wr / f:+.3e}, {wp / f:+.4f}; "
                  f"periods {2 * np.pi / wp / 3600:.2f} h and {2 * np.pi / abs(wr) / DAY:.0f} days; slow root: largest term {big}, "
                  f"smallest {small}; QG value {ch13.rossby_omega(k, 0.0, beta, f, c) / f:+.3e}; regime {ch13.shallow_water_regime(wr, f)['regime']}")
    k = 2.0 * np.pi / I["wavelengths"][1]
    w3 = 3.0 * f
    t3 = ch13.dispersion_term_sizes(np.sqrt(w3 ** 2 - f ** 2) / c0, 0.0, c0, f, beta, w3)
    print(f"slip #7 check: on the fast root with omega = 3f the terms are omega3 {t3['omega3']:.2e}, gravity {t3['gravity']:.2e}, "
          f"rotation {t3['rotation']:.2e}, beta {t3['beta']:.2e} -> the omega^3 term is the largest")

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
    for lat, label in hemispheres(I["lat"]):
        ff = ch13.coriolis_parameter(lat)
        orb = ch13.poincare_orbit(np.linspace(0, 2 * np.pi / ch13.poincare_omega(k, ff, c0), 200), k, 0.5, I["ocean_H"], ff)
        ax[0].plot(orb["u"] * 100, orb["v"] * 100, label=f"{label}: {orb['sense']}")
        cg = ch13.poincare_group_velocity(k, 0.0, ff, c0)
        print(f"Poincaré wave at {label}, wavelength {I['wavelengths'][1] / 1e3:.0f} km: omega/|f| = {orb['omega'] / abs(ff):.3f}, axis ratio "
              f"{orb['axis_ratio']:.3f}, {orb['sense']}; c_phase = {orb['omega'] / k:.1f} m/s, c_group = {cg[0]:.1f} m/s, product/c² = "
              f"{orb['omega'] / k * cg[0] / c0 ** 2:.6f}")
    ax[0].set(aspect="equal", xlabel="u [cm/s]", ylabel="v [cm/s]", title="Poincaré wave: current ellipse (η̂ = 0.5 m)")
    ax[0].legend(fontsize=8)
    q = I["q_inertial"]
    u, v, xi, yi = ch13.inertial_oscillation(np.linspace(0, ch13.inertial_period(f), 200), q, 0.0, f)
    print(f"inertial circle for q = {q} m/s at 35°: radius {ch13.inertial_radius(q, f) / 1e3:.2f} km, period "
          f"{ch13.inertial_period(f) / 3600:.2f} h; path closes to {np.hypot(xi[-1], yi[-1]):.1e} m; speed constant to {np.ptp(np.hypot(u, v)):.1e}")
    ax[1].plot(xi / 1e3, yi / 1e3, color=COLORS["accent"])
    ax[1].set(aspect="equal", xlabel="x [km]", ylabel="y [km]", title="inertial circle (f > 0: clockwise)")

    Lam0, Lam1 = ch13.rossby_radius(c0, f), ch13.rossby_radius(c1, f)
    Lam2 = ch13.rossby_radius_two_layer(I["H1"], I["rho_ocean"] - I["drho"], I["rho_ocean"], f)
    print(f"Rossby radii at 35°: external {Lam0 / 1e3:.0f} km, first baroclinic (uniform N) {Lam1 / 1e3:.1f} km, two-layer "
          f"(H1 = {I['H1']} m, drho = {I['drho']} kg/m³) {Lam2 / 1e3:.1f} km, Eady radius NH/f (atmosphere) "
          f"{ch13.rossby_radius_internal(I['atm_N'], I['atm_H'], f, with_pi=False) / 1e3:.0f} km")
    for lat, label in hemispheres(I["lat"]):
        ff = ch13.coriolis_parameter(lat)
        side = ch13.kelvin_decay_side(ff, +1)
        d = 1 if ff > 0 else -1
        eta, uk = ch13.kelvin_wave(0.0, Lam0, 0.0, 0.5, k, I["ocean_H"], ff)
        res = ch13.kelvin_residuals(1.0e5, 0.7 * Lam0, 3000.0, 0.5, k, I["ocean_H"], ff)
        print(f"Kelvin wave at {label}: travelling toward +x is {'trapped' if side['trapped'] else 'not trapped'} (coast on the "
              f"{side['coast_on']}); default direction {d:+d}; eta(y = Lambda)/eta0 = {eta / 0.5:.4f} (1/e = {np.exp(-1):.4f}), "
              f"u under the crest {uk * 100:+.2f} cm/s; residuals {res['continuity']:.1e}, {res['x_momentum']:.1e}, {res['y_geostrophy']:.1e}")
    y = np.linspace(0, 4 * Lam0, 100)
    ax[2].plot(y / Lam0, ch13.kelvin_wave(0.0, y, 0.0, 0.5, k, I["ocean_H"], f)[0], color=COLORS["accent"])
    ax[2].set(xlabel="distance from the coast y/Λ", ylabel="η [m]", title=f"Kelvin wave: trapped within Λ = {Lam0 / 1e3:.0f} km")
    save(fig, out, "waves_orbits_kelvin")
    save(ch13.fig_poincare_kelvin_dispersion(), out, "waves_dispersion")
    save(ch13.fig_kelvin_sections(), out, "waves_kelvin_sections")

    with Timer("C-grid model: one period of a Poincaré wave and of a Kelvin wave"):
        n = 32 if args.fast else 64
        L = 6.0e6
        md = ch13.ShallowWater(n, n, L, L, I["ocean_H"], f)
        kk, ll = 2 * 2 * np.pi / L, 2 * np.pi / L
        gr = md.grids()
        F = lambda XY, t: ch13.poincare_fields(XY[0], XY[1], t, kk, ll, 0.5, I["ocean_H"], f)  # noqa: E731
        st = ch13.make_state(md, F(gr["center"], 0)[0], F(gr["u"], 0)[1], F(gr["v"], 0)[2])
        T = 2 * np.pi / ch13.poincare_omega(np.hypot(kk, ll), f, c0)
        r = ch13.run(md, st, T, save_every=10 ** 6)
        print(f"  Poincaré wave on {n}² cells after one period ({T / 3600:.2f} h): error {np.max(np.abs(r['eta'][-1] - F(gr['center'], T)[0])) / 0.5:.2e} "
              f"of the amplitude, energy change {r['energy'][-1] / r['energy'][0] - 1:+.1e}, volume change "
              f"{abs(ch13.volume(md, dict(eta=r['eta'][-1])) - ch13.volume(md, st)):.1e} m³ of {0.5 * L * L:.1e}")
        for lat, label in hemispheres(I["lat"]):
            mk = ch13.ShallowWater(n, 2 * n, 8.0e6, 1.6e7, I["ocean_H"], ch13.coriolis_parameter(lat), bc="channel")
            sk = ch13.kelvin_state(mk, 0.5)
            rk = ch13.run(mk, sk, 8.0e6 / c0, save_every=10 ** 6)
            print(f"  Kelvin wave at {label} in a channel after one period: differs from its start by "
                  f"{np.max(np.abs(rk['eta'][-1] - sk['eta'])) / 0.5:.2e} of the amplitude; max |v| {np.max(np.abs(rk['v'][-1])):.1e} m/s")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
