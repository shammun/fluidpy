"""§12.10 (C13): one- and two-equation models — the k-epsilon eddy viscosity (12.104), the source terms of (12.103) and
(12.105), the decay of homogeneous turbulence (closed form against solve_ivp; C_eps2 fixes the decay exponent), the von
Karman constant the standard constants imply in a log layer, and (optional demo, label qualitative) a channel with wall
functions against the mixing-length channel.

Run: ``.venv/Scripts/python.exe scripts/ch12_k_epsilon.py --no-show``
Figure -> outputs/ch12/c13_k_epsilon.png.   Constants: the standard public set ch12.K_EPSILON_CONSTANTS.
Our decay case: e0 = 0.5 m2/s2, eps0 = 2 m2/s3.
"""
from __future__ import annotations

import numpy as np

from ch12_common import B_LOG, COLORS, KAPPA, Timer, finish, parse_args, save, setup

from fluidpy import ch12_turbulence as ch12


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    K = ch12.K_EPSILON_CONSTANTS
    print(f"standard constants: {K}")
    e0, eps0 = 0.5, 2.0
    print(f"nu_T = C_mu e^2/eps = {ch12.k_epsilon_eddy_viscosity(e0, eps0):.5f} m2/s, length e^1.5/eps = {ch12.k_epsilon_length_scale(e0, eps0):.4f} m")
    print(f"sources with production = dissipation: de/dt = {ch12.k_epsilon_rhs(e0, eps0, eps0)[0]:+.3f}, deps/dt = {ch12.k_epsilon_rhs(e0, eps0, eps0)[1]:+.3f}"
          f" (not zero: (C_eps1 - C_eps2) eps^2/e, balanced by diffusion in the log layer)")
    oc = ch12.one_equation_closure(e0, 0.05, c=0.55, C_eps=0.17, sigma_e=1.0)
    print(f"one-equation closure with illustrative constants c = 0.55, C_eps = 0.17, l_T = 5 cm: {oc}")
    t = np.geomspace(1e-3, 1e2, 60)
    e_c, eps_c, n_dec, t0 = ch12.k_epsilon_decay(e0, eps0, t)
    e_i, eps_i, _, _ = ch12.k_epsilon_decay(e0, eps0, t, method="ivp")
    print(f"decay: n = 1/(C_eps2 - 1) = {n_dec:.4f}, t0 = n e0/eps0 = {t0:.4f} s; closed form vs solve_ivp max relative "
          f"difference {np.max(np.abs(e_c / e_i - 1)):.1e} (e), {np.max(np.abs(eps_c / eps_i - 1)):.1e} (eps)")
    late = slice(-10, None)
    print(f"  late-time slope d ln e/d ln t = {np.polyfit(np.log(t[late]), np.log(e_c[late]), 1)[0]:.4f}")
    for c2 in (1.77, 1.92, 2.0):
        print(f"  C_eps2 = {c2}: decay exponent n = {1 / (c2 - 1):.3f}")
    print(f"log-layer von Karman constant implied by the standard set: {ch12.k_epsilon_loglayer_kappa():.4f} "
          f"(kappa^2 = sqrt(C_mu)(C_eps2 - C_eps1) sigma_eps)")
    with Timer("k-epsilon channel, Re_tau = 2000"):
        ke = ch12.k_epsilon_channel(2000.0, n=200, B=B_LOG)
    ml = ch12.channel_mixing_length(2000.0, KAPPA, 26.0, n=800)
    print(f"channel Re_tau = 2000 (QUALITATIVE): k-epsilon with wall functions converged {ke['converged']} in {ke['iterations']} iterations; "
          f"centreline U+ = {ke['Uplus'][-1]:.2f}, Cf = {ke['Cf']:.5f}; mixing length: {ml['U_cl_plus']:.2f}, {ml['Cf']:.5f}; "
          f"e+ at the wall node {ke['e_plus'][0]:.3f} = 1/sqrt(C_mu), at the centre {ke['e_plus'][-1]:.3f}")

    fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.6))
    ax[0].loglog(t, e_c, color=COLORS["teal"], label=r"$\bar e$ closed form")
    ax[0].loglog(t, e_i, "o", ms=3, color=COLORS["teal"], label=r"$\bar e$ solve_ivp")
    ax[0].loglog(t, eps_c, color=COLORS["rose"], label=r"$\bar\varepsilon$ closed form")
    ax[0].loglog(t, eps_i, "o", ms=3, color=COLORS["rose"])
    ax[0].set(xlabel="t [s]", ylabel=r"$\bar e$ [m$^2$/s$^2$], $\bar\varepsilon$ [m$^2$/s$^3$]", title=f"decay of homogeneous turbulence: n = {n_dec:.2f}")
    ax[0].legend(fontsize=8)
    ax[1].loglog(t, ch12.k_epsilon_eddy_viscosity(e_c, eps_c), color=COLORS["orange"], label=r"$\nu_T=C_\mu \bar e^2/\bar\varepsilon$ [m$^2$/s]")
    ax[1].loglog(t, ch12.k_epsilon_length_scale(e_c, eps_c), color=COLORS["accent"], label=r"$\bar e^{3/2}/\bar\varepsilon$ [m]")
    ax[1].set(xlabel="t [s]", title="the eddies grow while the energy decays")
    ax[1].legend(fontsize=8)
    ax[2].semilogx(ke["yplus"], ke["Uplus"], color=COLORS["accent"], label=r"k-$\varepsilon$ + wall functions")
    ax[2].semilogx(ml["yplus"][1:], ml["Uplus"][1:], color=COLORS["teal"], ls="--", label="mixing length")
    ax[2].set(xlabel=r"$y^+$", ylabel=r"$U^+$", title=r"channel, $Re_\tau$ = 2000 (models, qualitative)")
    ax[2].legend(fontsize=8)
    save(fig, out, "c13_k_epsilon")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
