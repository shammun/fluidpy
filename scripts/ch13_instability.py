"""Chapter 13, §13.16–§13.17 — barotropic and baroclinic instability.

* Rayleigh–Kuo: β − U″ must change sign; a jet stabilised by β (the criterion, and growth rates from the complex-path
  eigen-solver of chapter 11 with β added);
* the Eady problem: cut-off and fastest-growing wave computed from the dispersion relation, the growth rate in days for
  our atmosphere and ocean inputs, an independent Chebyshev eigenvalue, the westward tilt and the fluxes.

Run: ``.venv/Scripts/python.exe scripts/ch13_instability.py --no-show``   Figures -> outputs/ch13/instability_*.png
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
    U = lambda y: 1.0 / np.cosh(y) ** 2  # noqa: E731
    Up = lambda y: -2.0 * np.tanh(y) / np.cosh(y) ** 2  # noqa: E731
    Upp = lambda y: 4.0 / np.cosh(y) ** 2 - 6.0 / np.cosh(y) ** 4  # noqa: E731
    y = np.linspace(-6.0, 6.0, 1201)
    print("jet U = sech²(y) (non-dimensional): max U'' = 2/3, min U'' = -2")
    for b in (0.0, 0.3, 0.6, 0.7, -2.5):
        rk = ch13.rayleigh_kuo_criterion(y, U(y), b)
        print(f"  beta = {b:+.1f}: beta - U'' in [{rk['min']:+.3f}, {rk['max']:+.3f}], changes sign: {rk['changes_sign']} at y = "
              f"{np.round(rk['y_zero'], 3).tolist()}")
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.3))
    with Timer("Rayleigh–Kuo eigenvalues"):
        betas = (0.0, 0.3, 0.6) if args.fast else (0.0, 0.15, 0.3, 0.45, 0.6, 0.7)
        gr = []
        for b in betas:
            r = ch13.rayleigh_kuo_eigs(1.4, U, Up, Upp, b, bc="decay", y_max=16.0, parity="even")
            gr.append(r["growth_rate"])
            print(f"  k = 1.4, beta = {b:.2f}: c = {r['c']:.5f}, growth rate {r['growth_rate']:.5f}, stable: {r['stable']}")
    ax[0].plot(betas, gr, "o-", color=COLORS["accent"])
    ax[0].axvline(2.0 / 3.0, color=COLORS["rose"], ls="--", lw=0.8)
    ax[0].text(0.5, 0.9 * max(gr), "β = max U″\n(Rayleigh–Kuo)", fontsize=8, color=COLORS["rose"])
    ax[0].set(xlabel="β L²/U₀", ylabel="growth rate k c_i [U₀/L]", title="sech² jet, k L = 1.4: β stabilises")

    ac, fast, wl = ch13.eady_critical(), ch13.eady_fastest(), ch13.eady_wavelengths()
    print(f"Eady problem (computed, not typed): cut-off alpha_c H = {ac:.5f}; fastest wave alpha H = {fast['alphaH']:.5f}, growth "
          f"{fast['sigma_nd']:.5f} f U0/(N H), c = ({fast['cr_over_U0']:.3f} + {fast['ci_over_U0']:.5f} i) U0; wavelengths "
          f"{wl['cutoff_over_Lambda']:.4f} Lambda (cut-off) and {wl['fastest_over_Lambda']:.4f} Lambda (fastest); long-wave limit c_i/U0 = "
          f"{ch13.eady_phase_speed(0.0, 1.0).imag:.5f} = 1/(2 sqrt 3)")
    cases = (("atmosphere", I["atm_N"], I["atm_H"], I["atm_U0"]), ("ocean", I["ocean_N"], 1000.0, 0.11))
    for lat, label in hemispheres(I["lat"]):
        f = ch13.coriolis_parameter(lat)
        for name, N, H, U0 in cases:
            Lam = ch13.rossby_radius_internal(N, H, f, with_pi=False)
            s = ch13.eady_max_growth_rate(f, N, U0 / H)
            kf = fast["alphaH"] / Lam
            ne = ch13.eady_numeric_eigs(kf, 0.0, N, f, H, U0)
            fl = ch13.eady_fluxes(0.5 * H, fast["alphaH"], U0, H, N, f, 1.2 if name == "atmosphere" else I["rho_ocean"])
            print(f"  {name:10s} at {label}: Lambda = NH/|f| = {Lam / 1e3:7.1f} km, fastest wavelength {wl['fastest_over_Lambda'] * Lam / 1e3:7.0f} km, "
                  f"e-folding time {ch13.eady_time_scale(f, N, U0 / H) / DAY:6.2f} days (sigma = {s:.3e} 1/s; closed form at that k "
                  f"{ch13.eady_growth_rate(kf, 0.0, N, f, H, U0):.3e}, Chebyshev {ne['growth_rate']:.3e}); heat flux poleward: "
                  f"{fl['v_T_sign'] > 0}, w'rho' < 0: {fl['w_rho'] < 0}, tilt {np.rad2deg(fl['phase_tilt']):+.1f}° westward with height")
    aH = np.linspace(0.0, 3.2, 321)
    c = np.array([complex(ch13.eady_phase_speed(a, 1.0)) for a in aH])
    ax[1].plot(aH, aH * c.imag, color=COLORS["accent"], label="growth σ N H/(f U₀)")
    ax[1].plot(aH, c.real, color=COLORS["teal"], label="c_r/U₀ (growing or faster wave)")
    ax[1].plot(aH, np.array([complex(ch13.eady_phase_speed(a, 1.0, -1)).real for a in aH]), color=COLORS["teal"], ls=":")
    ax[1].axvline(ac, color=COLORS["rose"], ls="--", lw=0.8)
    ax[1].plot([fast["alphaH"]], [fast["sigma_nd"]], "o", color=COLORS["orange"])
    ax[1].set(xlabel="α H = K N H/|f|", title="Eady: growth rate and phase speed")
    ax[1].legend(fontsize=8)
    H = I["atm_H"]
    z = np.linspace(0.0, H, 81)
    md = ch13.eady_mode(z, fast["alphaH"], I["atm_U0"], H)
    xw = np.linspace(0.0, 2.0, 121)                       # in wavelengths
    P = np.real(md["p_hat"][:, None] * np.exp(2j * np.pi * xw[None, :]))
    ax[2].contourf(xw, z / 1e3, P, 12, cmap="RdBu_r")
    ax[2].plot(-md["phase"] / (2 * np.pi) + 1.0, z / 1e3, color=COLORS["ink"], lw=1.5)
    ax[2].set(xlabel="x / wavelength (east →)", ylabel="z [km]", title="fastest Eady wave: pressure tilts westward with height")
    save(fig, out, "instability_rayleigh_kuo_eady")
    fz = ch13.eady_fluxes(z, fast["alphaH"], I["atm_U0"], H, I["atm_N"], ch13.coriolis_parameter(I["lat"]), 1.2)
    print(f"  fluxes of the fastest wave: w'rho' negative at every interior level: {bool(np.all(fz['w_rho'][1:-1] < 0))}; v'rho' uniform in z to "
          f"{np.ptp(fz['v_rho']) / abs(np.mean(fz['v_rho'])):.1e}; phase difference top − bottom {np.rad2deg(md['phase'][-1] - md['phase'][0]):.1f}°")
    print(f"  beyond the cut-off (alpha H = 3): c = {ch13.eady_phase_speed(3.0, 1.0):.4f} and {ch13.eady_phase_speed(3.0, 1.0, -1):.4f} U0 "
          f"(two neutral edge waves), fluxes {ch13.eady_fluxes(0.5 * H, 3.0, I['atm_U0'], H, I['atm_N'], ch13.coriolis_parameter(I['lat']), 1.2)['w_rho']:.1e}")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
