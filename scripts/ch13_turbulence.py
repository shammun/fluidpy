"""Chapter 13, §13.18 — two-dimensional (geostrophic) turbulence.

* Fjørtoft's argument with our own wavenumber pair: more energy goes to the larger scale, more enstrophy to the smaller;
* the two inertial ranges and their exponents from dimensional analysis (exact rationals);
* a pseudo-spectral run of the barotropic vorticity equation (ours): inviscid invariants, the energy centroid moving to
  large scales, and — with β — the flow organising into zonal bands on the scale of the Rhines length.

The long runs are read from ``reference/ch13`` when present (``scripts/ch13_make_caches.py``); ``--fast`` runs small
live versions instead.

Run: ``.venv/Scripts/python.exe scripts/ch13_turbulence.py --no-show``   Figure -> outputs/ch13/turbulence.png
"""
from __future__ import annotations

from fractions import Fraction

import numpy as np
from ch13_common import COLORS, Timer, finish, parse_args, save, setup

from fluidpy import ch13_geophysical_fluid_dynamics as ch13
from fluidpy.core.dimensional import pi_groups


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch13.illustrative_inputs()
    K0 = 1.0
    ft = ch13.fjortoft_transfer(K0, K0 / 3.0, 1.5 * K0)
    A = np.array([[1.0, 1.0], [(K0 / 3.0) ** 2, (1.5 * K0) ** 2]])
    S = np.linalg.solve(A, [1.0, K0 ** 2])
    print(f"Fjørtoft, K1 = K0/3 and K2 = 3K0/2: energy shares S1 = {ft['S1']:.4f}, S2 = {ft['S2']:.4f} (direct solve {S[0]:.4f}, {S[1]:.4f}); "
          f"energy ratio S1/S2 = {ft['energy_ratio']:.4f}, enstrophy ratio = {ft['enstrophy_ratio']:.4f} -> energy to the large scale, "
          "enstrophy to the small one")
    gE = pi_groups({"S": "m^3/s^2", "eps": "m^2/s^3", "K": "1/m"})[0]     # one group: S eps^(-2/3) K^(5/3)
    gZ = pi_groups({"S": "m^3/s^2", "alpha": "1/s^3", "K": "1/m"})[0]     # one group: S alpha^(-2/3) K^3
    print(f"dimensional analysis (Pi theorem, exact rationals): energy cascade S ~ eps^{-gE['eps']} K^{-gE['K']}; enstrophy "
          f"cascade S ~ alpha^{-gZ['alpha']} K^{-gZ['K']}")
    assert (-gE["eps"], -gE["K"], -gZ["alpha"], -gZ["K"]) == (Fraction(2, 3), Fraction(-5, 3), Fraction(2, 3), Fraction(-3))
    for u, name in zip(I["u_rms"], ("atmosphere", "ocean")):
        print(f"Rhines length for u_rms = {u} m/s at 35° ({name}): {ch13.rhines_length(u, ch13.beta_parameter(I['lat'])) / 1e3:.0f} km")

    with Timer("barotropic runs"):
        runs = {}
        for name in ("turbulence_f", "turbulence_beta"):
            r = None if args.fast else ch13.load_reference_run(name)
            src = "cache"
            if r is None:
                r = ch13.reference_run(name, fast=True)
                src = "live (fast)"
            runs[name] = r
            L = float(r["L"])
            zf = [ch13.zonal_energy_fraction(z, L) for z in (r["zeta"][0], r["zeta"][-1])]
            print(f"  {name} [{src}; {r['zeta'].shape[-1]}², beta = {float(r['beta'])}, nu = {float(r['nu'])}]: to t = {r['t'][-1]:.0f} energy "
                  f"{r['energy'][-1] / r['energy'][0]:.3f}, enstrophy {r['enstrophy'][-1] / r['enstrophy'][0]:.3f} of the initial values; K_E "
                  f"{r['K_E'][0]:.2f} -> {r['K_E'][-1]:.2f}, K_Z {r['K_Z'][0]:.2f} -> {r['K_Z'][-1]:.2f}; zonal energy fraction "
                  f"{zf[0]:.3f} -> {zf[1]:.3f}")
        rb = runs["turbulence_beta"]
        urms = float(np.sqrt(2.0 * rb["energy"][-1]))
        print(f"  beta run: Rhines length sqrt(u_rms/beta) = {ch13.rhines_length(urms, float(rb['beta'])):.3f} box units against 1/K_E = "
              f"{1.0 / rb['K_E'][-1]:.3f} (qualitative: same order)")
        z0 = ch13.random_vorticity(48, 2 * np.pi, seed=1)
        dt = ch13.barotropic_time_step(z0, 2 * np.pi)
        ri = ch13.barotropic_run(48, 2 * np.pi, z0, t_end=1.0, dt=dt, save_every=10 ** 6)
        print(f"  inviscid check (48², t = 1, dt = {dt:.4f}): energy drift {ri['energy'][-1] / ri['energy'][0] - 1:+.1e}, enstrophy drift "
              f"{ri['enstrophy'][-1] / ri['enstrophy'][0] - 1:+.1e}")
    fig, ax = plt.subplots(1, 4, figsize=(19, 4.2))
    for a_, name, ttl in zip(ax[:2], ("turbulence_f", "turbulence_beta"), ("β = 0: isotropic eddies grow", "β ≠ 0: zonal bands")):
        r = runs[name]
        v = np.max(np.abs(r["zeta"][-1]))
        a_.pcolormesh(r["x"], r["x"], r["zeta"][-1], cmap="RdBu_r", vmin=-v, vmax=v, shading="auto")
        a_.set(aspect="equal", xlabel="x", ylabel="y", title=ttl + f" (vorticity, t = {r['t'][-1]:.0f})")
    r = runs["turbulence_f"]
    ax[2].plot(r["t"], r["K_E"], color=COLORS["accent"], label="K_E (energy)")
    ax[2].plot(r["t"], r["K_Z"], color=COLORS["amber"], label="K_Z (enstrophy)")
    ax[2].set(xlabel="t [L/u_rms units]", ylabel="mean wavenumber", title="energy to large scales, enstrophy to small")
    ax[2].legend(fontsize=8)
    for idx, col in ((0, COLORS["muted"]), (-1, COLORS["accent"])):
        K, E = ch13.barotropic_spectrum(r["zeta"][idx], float(r["L"]))
        ax[3].loglog(K, E, color=col, label=f"t = {r['t'][idx]:.0f}")
    Kr = np.array([6.0, 20.0])
    ax[3].loglog(Kr, 2e-2 * (Kr / 6.0) ** -3.0, color=COLORS["rose"], ls="--", label="K⁻³ (guide only)")
    ax[3].set(xlabel="K", ylabel="E(K)  (Σ E ΔK = ½⟨u² + v²⟩)", title="energy spectrum (64²: slope not resolved)")
    ax[3].legend(fontsize=8)
    save(fig, out, "turbulence")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
