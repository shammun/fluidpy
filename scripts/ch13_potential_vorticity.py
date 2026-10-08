"""Chapter 13, §13.13 — potential vorticity (ζ + f)/h is conserved following the motion.

* a column crossing a step: the relative vorticity it must gain; the linearised streamline of an eastward flow (a
  standing Rossby wave of wavelength 2π sqrt(U/β)) and of a westward flow (no wave, upstream influence) — ours;
* the nonlinear C-grid model: potential vorticity carried by marked particles, with its measured drift.

Run: ``.venv/Scripts/python.exe scripts/ch13_potential_vorticity.py --no-show``   Figure -> outputs/ch13/pv.png
"""
from __future__ import annotations

import numpy as np
from ch13_common import COLORS, Timer, finish, hemispheres, parse_args, save, setup

from fluidpy import ch13_geophysical_fluid_dynamics as ch13


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch13.illustrative_inputs()
    beta = ch13.beta_parameter(I["lat"])
    h0, h1, U = 4200.0, 4000.0, I["U_mean"]
    for lat, label in hemispheres(I["lat"]):
        f = ch13.coriolis_parameter(lat)
        z = ch13.step_vorticity(f, h0, h1)
        print(f"column crossing a step {h0:.0f} -> {h1:.0f} m at {label}: zeta = {z:+.3e} 1/s ({'anticyclonic' if z * f < 0 else 'cyclonic'}); "
              f"q before {ch13.potential_vorticity(0.0, f, h0):+.4e}, after {ch13.potential_vorticity(z, f, h1):+.4e} 1/(m s)")
    f = ch13.coriolis_parameter(I["lat"])
    x = np.linspace(-6.0e6, 1.2e7, 1201)
    east = ch13.flow_over_step(x, U, beta, f, h0, h1)
    west = ch13.flow_over_step(x, -U, beta, f, h0, h1)
    print(f"eastward flow U = {U} m/s: stationary wavelength {east['wavelength'] / 1e3:.0f} km (2 pi sqrt(U/beta) = "
          f"{ch13.stationary_rossby_wavelength(U, beta) / 1e3:.0f} km), largest southward excursion {east['Y'].min() / 1e3:.0f} km, zeta just past "
          f"the step {ch13.flow_over_step(1.0, U, beta, f, h0, h1)['zeta']:+.3e} 1/s")
    print(f"westward flow: no wave; displacement {ch13.flow_over_step(0.0, -U, beta, f, h0, h1)['Y'] / 1e3:+.0f} km already at the step "
          f"(upstream influence over {west['decay_length'] / 1e3:.0f} km), {west['Y'][0] / 1e3:+.0f} km far downstream; monotonic: "
          f"{bool(np.all(np.diff(west['Y']) >= 0) or np.all(np.diff(west['Y']) <= 0))}")
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.3))
    ax[0].plot(x / 1e3, east["Y"] / 1e3, color=COLORS["accent"], label="eastward (flow →)")
    ax[0].plot(x / 1e3, west["Y"] / 1e3, color=COLORS["rose"], label="westward (flow ←)")
    ax[0].axvline(0, color=COLORS["muted"], lw=0.8)
    ax[0].set(xlabel="x [km] (step at 0)", ylabel="streamline displacement Y [km]", title="flow over a step on a β-plane (linearised, ours)")
    ax[0].legend(fontsize=8)

    with Timer("nonlinear shallow-water run with particles"):
        n = 32 if args.fast else 48
        r = ch13.load_reference_run("pv_particles") if not args.fast else None
        src = "cache reference/ch13/pv_particles.npz"
        if r is None:
            r = ch13.reference_run("pv_particles", fast=args.fast)
            src = "live run"
        q = r["q"]
        drift = np.max(np.abs(q / q[0] - 1.0), axis=1)
        print(f"  ({src}; {n if src == 'live run' else r['eta'].shape[-1]}² cells) potential vorticity on {q.shape[1]} particles over "
              f"{r['t'][-1] * abs(float(r['f'])) / (2 * np.pi):.1f} inertial periods: values span {100 * np.ptp(q[0]) / np.mean(q[0]):.0f} % of the mean; "
              f"largest drift along a path {100 * drift.max():.2f} %; particles moved up to "
              f"{np.max(np.hypot(r['xp'] - r['xp'][0], r['yp'] - r['yp'][0])) / 1e3:.0f} km")
    X, Y = np.meshgrid(r["x"], r["y"])
    ax[1].contourf(X / 1e3, Y / 1e3, r["eta"][-1], 12, cmap="Blues")
    ax[1].plot(r["xp"] / 1e3, r["yp"] / 1e3, lw=0.8)
    ax[1].set(aspect="equal", xlabel="x [km]", ylabel="y [km]", title="particle paths round a geostrophic vortex")
    ax[2].plot(r["t"] * abs(float(r["f"])) / (2 * np.pi), q / q[0], lw=0.8)
    ax[2].set(xlabel="time [inertial periods]", ylabel="q(t)/q(0) on each particle", ylim=(0.97, 1.03),
              title="potential vorticity following the motion")
    save(fig, out, "pv")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
