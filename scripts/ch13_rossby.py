"""Chapter 13, §13.15 — Rossby waves.

* frequency, phase speed (always westward) and group velocity (west for long waves, east for short ones) of barotropic
  and first-baroclinic waves; the maximum frequency; the long-wave speed and the time to cross a basin;
* the stationary wave in a westerly flow;
* a wave packet evolved exactly with the linear quasi-geostrophic equation: crests go west, the envelope goes east.

Run: ``.venv/Scripts/python.exe scripts/ch13_rossby.py --no-show``   Figures -> outputs/ch13/rossby_*.png
"""
from __future__ import annotations

import numpy as np
from ch13_common import COLORS, DAY, YEAR, finish, parse_args, save, setup

from fluidpy import ch13_geophysical_fluid_dynamics as ch13


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch13.illustrative_inputs()
    c0 = ch13.long_wave_speed(I["ocean_H"])
    c1 = ch13.baroclinic_mode_speed(I["ocean_N"], I["ocean_H"], 1)
    for lat_name in ("lat_rossby", "lat"):
        lat = I[lat_name]
        f, beta = ch13.coriolis_parameter(lat), ch13.beta_parameter(lat)
        mx0 = ch13.rossby_max_frequency(beta, f, c0)
        print(f"{np.rad2deg(lat):.0f}° N, barotropic (c = {c0:.1f} m/s, Lambda = {ch13.rossby_radius(c0, f) / 1e3:.0f} km): omega_max/f = "
              f"{mx0['omega_max'] / f:.2f} — not small, and Lambda is of the size of the planet: use the non-divergent form (f0/c -> 0) there")
        mx = ch13.rossby_max_frequency(beta, f, c1)
        print(f"{np.rad2deg(lat):.0f}° N, baroclinic (c = {c1:.2f} m/s, Lambda = {ch13.rossby_radius(c1, f) / 1e3:.1f} km): omega_max/f = "
              f"{mx['omega_max'] / f:.2e} -> shortest period {2 * np.pi / mx['omega_max'] / DAY:.0f} days at wavelength "
              f"{2 * np.pi / abs(mx['k']) / 1e3:.0f} km; long-wave speed {ch13.rossby_long_wave_speed(beta, f, c1) * 100:+.2f} cm/s; a 6000 km "
              f"basin is crossed in {ch13.basin_crossing_time(6.0e6, lat, c1) / YEAR:.2f} years")
    lat = I["lat"]
    f, beta = ch13.coriolis_parameter(lat), ch13.beta_parameter(lat)
    for name, c in (("barotropic, non-divergent", np.inf), ("first baroclinic", c1)):
        for lam in I["wavelengths"]:
            k = -2.0 * np.pi / lam
            w = ch13.rossby_omega(k, 0.0, beta, f, c)
            cg = ch13.rossby_group_velocity(k, 0.0, beta, f, c)
            h = 1e-20
            cs = np.imag(ch13.rossby_omega(k + 1j * h, 0.0, beta, f, c)) / h
            print(f"{name} wave at 35°, wavelength {lam / 1e3:.0f} km: omega = {w:+.3e} rad/s (period {2 * np.pi / abs(w) / DAY:.1f} days), "
                  f"c_x = {ch13.rossby_phase_speed(k, 0.0, beta, f, c):+.4f} m/s, c_gx = {cg[0]:+.4f} m/s (complex step {cs:+.4f}): phase west, "
                  f"energy {'east' if cg[0] > 0 else 'west'}")
    U = I["U_mean"]
    lam_s = ch13.stationary_rossby_wavelength(U, beta)
    print(f"stationary wave in a mean flow of {U} m/s at 35°: wavelength {lam_s / 1e3:.0f} km "
          f"({2 * np.pi * ch13.EARTH_RADIUS_MEAN * np.cos(lat) / lam_s:.1f} waves round the latitude circle); check c_x = "
          f"{ch13.rossby_phase_speed(2 * np.pi / lam_s, 0.0, beta, U=U):+.1e} m/s")
    wq = 0.5 * ch13.rossby_max_frequency(beta, f, c1)["omega_max"]
    cir = ch13.rossby_omega_circle(wq, beta, f, c1)
    th = 0.7
    kc, lc = cir["center_k"] + cir["radius"] * np.cos(th), cir["radius"] * np.sin(th)
    print(f"circle of constant omega = omega_max/2 (baroclinic): centre k = {cir['center_k']:.3e}, radius {cir['radius']:.3e} rad/m; a point on "
          f"it returns omega/omega_q = {ch13.rossby_omega(kc, lc, beta, f, c1) / wq:.12f}")
    save(ch13.fig_rossby_dispersion(), out, "rossby_dispersion")

    # a barotropic packet shorter than 2 pi Lambda (energy east) on a periodic strip
    nx, Lx = 512, 6.0e7
    x = np.arange(nx) * Lx / nx
    y = np.arange(8) * 1.0e6
    k0 = -2.0 * np.pi / 3.0e6
    env = np.exp(-((x - 0.3 * Lx) ** 2) / (2 * (4.0e6) ** 2))
    eta0 = (env * np.cos(k0 * x))[None, :] * np.ones((8, 1))
    T = 6.0 * DAY
    eta = ch13.qg_linear_evolve(eta0, x, y, T, beta, f, c0)
    cg = ch13.rossby_group_velocity(k0, 0.0, beta, f, c0)[0]
    cp = ch13.rossby_phase_speed(k0, 0.0, beta, f, c0)
    centre = lambda e: np.sum(x * e ** 2) / np.sum(e ** 2)  # noqa: E731
    moved = centre(eta[0]) - centre(eta0[0])
    print(f"packet of 3000 km waves after {T / DAY:.0f} days: envelope moved {moved / 1e3:+.0f} km (group velocity × time = {cg * T / 1e3:+.0f} km), "
          f"crests moved {cp * T / 1e3:+.0f} km")
    fig, ax = plt.subplots(figsize=(10, 3.8))
    ax.plot(x / 1e6, eta0[0], color=COLORS["muted"], label="t = 0")
    ax.plot(x / 1e6, eta[0], color=COLORS["accent"], label=f"t = {T / DAY:.0f} days")
    ax.set(xlabel="x [1000 km]", ylabel="η (normalised)", xlim=(0, 45), title="Rossby packet: crests west, energy east (short waves)")
    ax.legend()
    save(fig, out, "rossby_packet")
    kk, a = ch13.rossby_packet_spectrum(k0, abs(k0) / 8.0, 33)
    amp = ch13.qg_linear_evolve_1d(a, kk, 0.0, T, beta, f, c0)
    print(f"sum of 33 modes: |sum of amplitudes| at t = 0 is {abs(np.sum(a)):.3f}, eta(x = c_g t, t) = "
          f"{abs(np.sum(amp * np.exp(1j * kk * cg * T))):.3f} (the envelope's peak has moved with the group velocity)")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
