"""Chapter 13, §13.14 — internal waves with rotation (inertia–gravity waves).

* the dispersion relation between f and N, its three limiting forms and their errors; f = 0 reproduces chapter 7;
* group velocity perpendicular to the wavevector (checked by complex-step differentiation), phase up = energy down;
* the WKB solution in a slowly varying N(z) against the numerical solution of the vertical-structure equation;
* the current hodograph in both hemispheres; lee waves over sinusoidal hills.

Run: ``.venv/Scripts/python.exe scripts/ch13_internal_waves.py --no-show``   Figure -> outputs/ch13/internal_waves.png
"""
from __future__ import annotations

import numpy as np
from ch13_common import COLORS, finish, hemispheres, parse_args, save, setup

from fluidpy import ch13_geophysical_fluid_dynamics as ch13
from fluidpy.core.waves import internal_wave_omega


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch13.illustrative_inputs()
    N, f = I["ocean_N"], ch13.coriolis_parameter(I["lat"])
    print(f"band at 35° for N = {N} rad/s: {abs(f):.3e} < omega < {N:.3e} rad/s (N/f = {N / abs(f):.1f}); periods "
          f"{2 * np.pi / N / 60:.1f} min to {2 * np.pi / abs(f) / 3600:.1f} h")
    k = 1.0e-3
    for m in (2.0e-4, 5.0e-3, 5.0e-2):
        w = ch13.inertia_gravity_omega(k, m, N, f)
        reg = ch13.inertia_gravity_regime(w, N, f)
        cg = ch13.inertia_gravity_group_velocity(k, m, N, f)
        h = 1.0e-20
        cs = (np.imag(ch13.inertia_gravity_omega(k + 1j * h, m, N, f)) / h, np.imag(ch13.inertia_gravity_omega(k, m + 1j * h, N, f)) / h)
        print(f"  m/k = {m / k:6.1f}: omega = {w:.4e} rad/s ({w / abs(f):.2f} f, {w / N:.3f} N), regime {reg['regime']:13s} errors in m² "
              f"non-rotating {reg['err_nonrotating']:+.2e}, hydrostatic {reg['err_hydrostatic']:+.2e}; c_g = ({cg[0]:+.3e}, {cg[1]:+.3e}) m/s "
              f"(complex step ({cs[0]:+.3e}, {cs[1]:+.3e})); K.c_g = {k * cg[0] + m * cg[1]:.1e}")
    print(f"  f = 0: {ch13.inertia_gravity_omega(k, 5e-3, N, 0.0):.6e} rad/s = chapter 7's N cos(theta) {internal_wave_omega(k, 5e-3, N):.6e}; "
          f"m² back from omega: {ch13.inertia_gravity_m2(k, 0.0, ch13.inertia_gravity_omega(k, 5e-3, N, f), N, f):.6e} (25.0e-6 put in)")

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.3))
    kk = np.geomspace(1e-5, 1e-1, 300)
    for m, col in ((1e-3, COLORS["accent"]), (1e-2, COLORS["teal"]), (1e-1, COLORS["orange"])):
        ax[0].semilogx(kk, ch13.inertia_gravity_omega(kk, m, N, f) / N, color=col, label=f"m = {m:g} rad/m")
    ax[0].axhline(abs(f) / N, color=COLORS["muted"], ls="--", lw=0.8)
    ax[0].text(1.2e-5, abs(f) / N + 0.02, "ω = |f|", fontsize=8)
    ax[0].set(xlabel="k [rad/m]", ylabel="ω / N", title="inertia–gravity waves: ω between |f| and N")
    ax[0].legend(fontsize=8)

    w = 4.0e-4
    print("WKB against the numerical solution, N = N0 (1 + z/2D) on -D < z < 0, k = 1e-3 rad/m, omega = 4e-4 rad/s:")
    for D in (500.0, 1000.0, 2000.0, 4000.0):
        z = np.linspace(-D, 0.0, 4001)
        e = ch13.wkb_error(z, lambda q, D=D: N * (1.0 + 0.5 * q / D), k, w, f)
        print(f"  D = {D:6.0f} m: Hm = {e['Hm']:5.1f}, largest relative error {e['max_rel_error']:.4f} (error × Hm = {e['max_rel_error'] * e['Hm']:.3f})")
    D = 1000.0
    z = np.linspace(-D, 0.0, 2001)
    Nf = lambda q: N * (1.0 + 0.5 * q / D)  # noqa: E731
    m = np.sqrt(ch13.inertia_gravity_m2(k, 0.0, w, Nf(z), f))
    wk = ch13.wkb_vertical_structure(z, m)
    dm0 = (-3 * m[0] + 4 * m[1] - m[2]) / (z[2] - z[0])
    num = ch13.vertical_structure_solve(z, Nf, k, w, f, w0=wk[0], dw0=(1j * m[0] - 0.5 * dm0 / m[0]) * wk[0])
    ax[1].plot(wk.real, z, color=COLORS["accent"], label="WKB")
    ax[1].plot(num.real, z, color=COLORS["ink"], ls=":", label="numerical")
    ax[1].plot(1 / np.sqrt(m), z, color=COLORS["muted"], lw=0.8, label="envelope m^{-1/2}")
    ax[1].set(xlabel="Re ŵ", ylabel="z [m]", title="vertical structure in a varying N(z)")
    ax[1].legend(fontsize=8)
    u, v, ww = ch13.inertia_gravity_fields(0.0, z, 0.0, k, w, Nf, f)
    print(f"  WKB fields at x = t = 0: max |u| {np.nanmax(np.abs(u)):.2f}, max |v| {np.nanmax(np.abs(v)):.2f}, max |w| {np.nanmax(np.abs(ww)):.3f} "
          f"(per unit A0); |v|/|u| amplitude ratio {abs(f) / w:.3f} = |f|/omega")
    for lat, label in hemispheres(I["lat"]):
        ff = ch13.coriolis_parameter(lat)
        t = np.linspace(0.0, 2 * np.pi / w, 200)
        uh, vh = ch13.inertia_gravity_hodograph(t, w, ff)
        area = 0.5 * np.sum(uh[:-1] * np.diff(vh) - vh[:-1] * np.diff(uh))
        print(f"  hodograph at {label}: axis ratio {np.ptp(vh) / np.ptp(uh):.3f}, {'clockwise' if area < 0 else 'counter-clockwise'}")
        ax[2].plot(uh, vh, label=label)
    ax[2].set(aspect="equal", xlabel="u/û", ylabel="v/û", title="current hodograph (ω = 4.8 |f|)")
    ax[2].legend(fontsize=8)
    save(fig, out, "internal_waves")
    save(ch13.fig_inertia_gravity_orbit(), out, "internal_waves_orbit")

    U, Na, h0 = 12.0, I["atm_N"], 200.0
    print(f"lee waves for U = {U} m/s, N = {Na} rad/s: hills narrower than 2 pi U/N = {2 * np.pi * U / Na / 1e3:.2f} km do not radiate")
    for lam in (3.0e3, 1.2e4, 5.0e4):
        kh = 2 * np.pi / lam
        fld = ch13.lee_wave_field(np.linspace(0, lam, 9), np.array([0.0, 2000.0]), U, Na, kh, h0)
        mm = "evanescent" if fld["m"] is None else f"m = {fld['m']:.3e} rad/m (vertical wavelength {2 * np.pi / fld['m'] / 1e3:.1f} km)"
        print(f"  hills {lam / 1e3:5.1f} km apart: {mm}; tilt {fld['tilt']}; max w at the ground {np.max(np.abs(fld['w'][0])):.2f} m/s, "
              f"at 2 km {np.max(np.abs(fld['w'][1])):.2f} m/s")
    wf = lambda x, y, zz, t: np.cos(2e-3 * x + 5e-3 * zz - ch13.inertia_gravity_omega(2e-3, 5e-3, N, f) * t)  # noqa: E731
    print(f"w-equation residual of a plane wave: {ch13.w_equation_rotating_residual(wf, 10.0, 0.0, -50.0, 100.0, N, f, h=5.0, ht=50.0):.1e} "
          f"against a term size N²k² = {N ** 2 * 4e-6:.1e}")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
