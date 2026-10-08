"""Chapter 13 — geostrophic adjustment to the Rossby radius (ours — not in the book; it ties §13.5, §13.8, §13.11–§13.13
together).  A step in surface height is released: without rotation it flattens completely; with rotation gravity waves
carry away two thirds of the potential energy released and a geostrophic jet of width Λ remains.

The 1-D forward–backward march (``ch13.linear_1d_run``) is compared with the closed-form end state
(``ch13.geostrophic_adjustment_1d``), in both hemispheres.

Run: ``.venv/Scripts/python.exe scripts/ch13_adjustment.py --no-show``   Figure -> outputs/ch13/adjustment.png
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
    g, He, eta0 = 9.80665, 1.3, 0.5          # an equivalent depth of about a metre: a first-baroclinic-like layer
    f = ch13.coriolis_parameter(I["lat"])
    c = ch13.long_wave_speed(He)
    Lam = ch13.rossby_radius(c, f)
    en = ch13.adjustment_energy(eta0, He, f, rho=I["rho_ocean"])
    print(f"layer of equivalent depth {He} m at 35°: c = {c:.3f} m/s, Lambda = {Lam / 1e3:.1f} km, inertial period "
          f"{ch13.inertial_period(f) / 3600:.1f} h")
    print(f"closed form: jet speed at the step {abs(ch13.geostrophic_adjustment_1d(0.0, eta0, He, f)[1]):.3f} m/s = g eta0/c; energy per metre "
          f"of step: released {en['pe_released']:.3e} J/m, kept by the jet {en['ke_jet']:.3e} J/m, radiated {en['radiated']:.3e} J/m; "
          f"ratio kept/released = {en['ratio']:.6f}")
    print(f"  within |x| < 2 Lambda the ratio is {ch13.adjustment_energy(eta0, He, f, L=2 * Lam, rho=I['rho_ocean'])['ratio']:.4f}")
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.3))
    with Timer("1-D marches"):
        for n in ((400,) if args.fast else (400, 800)):
            Lx = 60.0 * Lam
            dx = Lx / n
            x = (np.arange(n) + 0.5 - 0.5 * n) * dx
            dt = 0.5 * dx / c
            Ti = ch13.inertial_period(f)
            nst = int(round(20.0 * Lam / c / dt))
            for lat, label in hemispheres(I["lat"]):
                ff = ch13.coriolis_parameter(lat)
                r = ch13.linear_1d_run(eta0 * np.sign(x), dx, dt, nst, He, ff)
                ex_eta, ex_v = ch13.geostrophic_adjustment_1d(x, eta0, He, ff)
                core = np.abs(x) < 5.0 * Lam
                last = r["t"] >= r["t"][-1] - Ti
                em, vm = r["eta"][last].mean(axis=0), r["v"][last].mean(axis=0)
                print(f"  n = {n}, {label}: after {r['t'][-1] / Ti:.2f} inertial periods the instantaneous eta is within "
                      f"{np.max(np.abs(r['eta'][-1] - ex_eta)[core]) / eta0:.3f} eta0 of the end state, its mean over the last inertial period "
                      f"within {np.max(np.abs(em - ex_eta)[core]) / eta0:.4f} eta0 (jet: {np.max(np.abs(vm - ex_v)[core]) / np.max(np.abs(ex_v)):.4f}); "
                      f"potential vorticity changed by {np.max(np.abs(r['pv'][1] - r['pv'][0])[3:-3]):.1e} 1/s, volume by "
                      f"{abs(np.sum(r['eta'][-1]) * dx):.1e} m²; jet sign {np.sign(vm[n // 2]):+.0f}")
                if n == 400:
                    ax[0].plot(x / Lam, em, label=f"model mean, {label}")
                    ax[1].plot(x / Lam, vm, label=f"model mean, {label}")
            if n == 400:
                ax[0].plot(x / Lam, ex_eta, color=COLORS["ink"], ls=":", label="closed form")
                r0 = ch13.linear_1d_run(eta0 * np.sign(x), dx, dt, nst, He, 0.0)
                ax[0].plot(x / Lam, r0["eta"][-1], color=COLORS["muted"], lw=1, label="f = 0 (all radiated)")
                print(f"  f = 0: largest |eta| within 5 Lambda at the same time {np.max(np.abs(r0['eta'][-1])[core]) / eta0:.3f} eta0 (nothing stays)")
    ax[0].set(xlim=(-6, 6), xlabel="x/Λ", ylabel="η [m]", title="adjusted surface")
    ax[0].legend(fontsize=7)
    ax[1].set(xlim=(-6, 6), xlabel="x/Λ", ylabel="v [m/s]", title="jet along the step (sign follows f)")
    ax[1].legend(fontsize=7)
    ax[2].bar(["released", "jet", "waves"], [en["pe_released"], en["ke_jet"], en["radiated"]],
              color=[COLORS["orange"], COLORS["teal"], COLORS["accent"]])
    ax[2].set(ylabel="energy per metre of step [J/m]", title=f"kept by the jet: {en['ratio']:.3f} of what is released")
    save(fig, out, "adjustment")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
