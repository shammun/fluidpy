"""Chapter 16, §16.3 — Womersley's oscillatory flow in a rigid pipe (16.69)–(16.100): profiles across the full diameter
at several Womersley numbers (our Fig. 16.14), the quasi-steady (16.87) and Stokes-layer (16.88) limits, flow amplitude
and phase lag against α (16.97)–(16.100), the wall shear (16.92, corrected), the small- and large-α forms of F₂ (16.93),
(16.96, corrected — slip S10), a harmonic series (16.83)–(16.86), and an independent Crank–Nicolson solution of (16.72).

Run: ``.venv/Scripts/python.exe scripts/ch16_womersley.py --no-show``   Figures -> outputs/ch16/womersley_*.png
"""
from __future__ import annotations

import numpy as np
from ch16_common import C_E, C_P, C_REF, C_TAU, C_U, C_VISC, Timer, finish, parse_args, save, setup
from scipy.special import jv

from fluidpy import ch16_biofluid_mechanics as ch16


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch16.illustrative_inputs()
    rho, mu = I["rho"], I["mu"]
    nu = mu / rho
    om = 2 * np.pi * I["HR"] / 60.0
    a, G = I["a_artery"], I["G_womersley"]
    al = ch16.womersley_number(a, om, nu)
    sl = ch16.stokes_layer(nu, om)
    print(f"artery (ours): a = {a * 1e3} mm, HR {I['HR']:.0f}/min (ω = {om:.3f} rad/s), ν = {nu:.3g} m²/s -> α (16.12) = {al:.3f} "
          f"({ch16.womersley_regime(al)}); Stokes layer δ_e = √(2ν/ω) = {sl['delta_e'] * 1e3:.3f} mm, α = √2·a/δ_e = "
          f"{np.sqrt(2) * a / sl['delta_e']:.3f}; Re·St/4 check: {ch16.womersley_from_re_st(2 * a * 0.1 / nu, 2 * a * om / 0.1):.3f}")
    r_in = np.array([0.3, 0.6, 0.9]) * a
    lhs = ch16.womersley_ode_residual(r_in, a, om, rho, mu, G)
    print(f"  (16.76) left side at r/a = 0.3, 0.6, 0.9: {np.round(np.real(lhs), 2)} (imag ~{np.max(np.abs(np.imag(lhs))):.1e}); "
          f"−G/μ = {-G / mu:.2f} (the printed +G/μ is slip S6)")
    z = ch16.I32 * 4.0
    print(f"  from-scratch J₀ series vs scipy at z = 4i^(3/2): |Δ| = {abs(ch16.bessel_j0_series(z) - jv(0, z)):.1e}; "
          f"J₁: {abs(ch16.bessel_j1_series(z) - jv(1, z)):.1e}")
    # profiles at four α over half a cycle (full diameter, trap T20)
    fig, axs = plt.subplots(1, 4, figsize=(18, 4.2), sharey=True)
    s = np.linspace(-1, 1, 201)
    phases = np.linspace(0, np.pi, 7)
    cmap = plt.get_cmap("viridis")
    for axk, alpha in zip(axs, (1.0, 3.0, 6.0, 12.0)):
        for j, ph in enumerate(phases):
            axk.plot(ch16.womersley_velocity_nd(np.abs(s), ph, alpha), s, color=cmap(j / len(phases)),
                     label=f"ωt = {np.degrees(ph):.0f}°" if alpha == 1.0 else None)
        axk.axvline(0, color=C_REF, lw=0.6)
        axk.set(title=f"α = {alpha:g} ({ch16.womersley_regime(alpha)})", xlabel="u/(Ga²/4μ)")
    axs[0].set_ylabel("r/a (full diameter)")
    axs[0].legend(fontsize=7)
    fig.suptitle("Womersley profiles (16.82) for −∂p/∂x = G cos ωt (ours, Fig. 16.14 style)")
    save(fig, out, "womersley_profiles")
    # limits
    for alpha in (0.3, 20.0):
        omx = alpha ** 2 * nu / a ** 2
        rr = np.linspace(0, a, 201)
        ex = ch16.womersley_velocity(rr, 0.37, a, omx, rho, mu, G)
        lim = (ch16.womersley_small_alpha if alpha < 1 else ch16.womersley_large_alpha)(rr, 0.37, a, omx, rho, mu, G)
        print(f"  α = {alpha}: max |u − limit| / max|u| = {np.max(np.abs(ex - lim)) / np.max(np.abs(ex)):.2e} "
              f"({'quasi-steady (16.87)' if alpha < 1 else 'Stokes layer (16.88)'})")
    # flow amplitude and phase against α
    A_ = np.geomspace(0.1, 100, 300)
    ratio = ch16.womersley_flow_ratio(A_)
    phi = np.degrees(ch16.womersley_phase_lag(A_))
    X, Y = ch16.womersley_XY(8.0)
    print(f"  α = 8: X = {X:.4f}, Y = {Y:.4f}, |Q|/|Q₀| = 8|Z|/α² = {ch16.womersley_flow_ratio(8.0):.4f}, "
          f"φ = {np.degrees(ch16.womersley_phase_lag(8.0)):.2f}°")
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.4))
    ax[0].loglog(A_, ratio, color=C_U, label="|Q|/|Q₀| = 8|1 − F₂|/α² (16.99)")
    ax[0].loglog(A_, 8 / A_ ** 2, color=C_REF, ls=":", label="8/α² (inertia-limited)")
    ax[0].set(xlabel="α", ylabel="|Q|/|Q₀|", ylim=(1e-3, 2), title="flow amplitude against Poiseuille's")
    ax[0].legend()
    ax[1].semilogx(A_, phi, color=C_P, label="φ = arctan(X/Y) (16.100)")
    ax[1].axhline(90, color=C_REF, ls=":")
    ax[1].set(xlabel="α", ylabel="phase lag of Q behind −∂p/∂x [°]", title="from damper (0°) to mass (90°)")
    ax[1].legend()
    # F2 approximations: error slopes (S10)
    al_s = np.array([0.05, 0.1, 0.2])
    e2 = np.abs(ch16.womersley_F2_series(al_s) - ch16.womersley_F2(al_s))
    al_l = np.array([20.0, 40.0, 80.0])
    ec = np.abs(ch16.womersley_F2_asymptote(al_l) - ch16.womersley_F2(al_l))
    ep = np.abs(ch16.womersley_F2_asymptote(al_l, printed_16_96=True) - ch16.womersley_F2(al_l))
    e0 = np.abs(2 / (ch16.I12 * al_l) - ch16.womersley_F2(al_l))
    sl_ = lambda x, e: np.polyfit(np.log(x), np.log(e), 1)[0]  # noqa: E731
    print(f"  error slopes: (16.93) series {sl_(al_s, e2):.2f} (expect 4); corrected (16.96) {sl_(al_l, ec):.2f} (−3); "
          f"printed (16.96) {sl_(al_l, ep):.2f} (−2); leading term alone {sl_(al_l, e0):.2f} (−2, smaller constant: "
          f"{e0[-1]:.2e} vs printed {ep[-1]:.2e} at α = 80)")
    AA = np.geomspace(5, 200, 50)
    ax[2].loglog(AA, np.abs(ch16.womersley_F2_asymptote(AA) - ch16.womersley_F2(AA)), color=C_U, label="corrected (16.96)")
    ax[2].loglog(AA, np.abs(ch16.womersley_F2_asymptote(AA, printed_16_96=True) - ch16.womersley_F2(AA)), color=C_TAU,
                 label="printed (16.96) (slip S10)")
    ax[2].loglog(AA, np.abs(2 / (ch16.I12 * AA) - ch16.womersley_F2(AA)), color=C_REF, ls=":", label="leading term only")
    ax[2].set(xlabel="α", ylabel="|error in F₂|", title="large-α asymptote of F₂")
    ax[2].legend()
    save(fig, out, "womersley_flow_phase")
    # wall shear, harmonics
    t = np.linspace(0, 2 * np.pi / om, 9)
    h = a * 1e-6
    tw_fd = -mu * (ch16.womersley_velocity(a + h, t, a, om, rho, mu, G) - ch16.womersley_velocity(a - h, t, a, om, rho, mu, G)) / (2 * h)
    print(f"  wall shear (16.92, corrected): |τ_w − (−μ ∂u/∂r)_FD| / max = "
          f"{np.max(np.abs(ch16.womersley_wall_shear(t, a, om, rho, mu, G) - tw_fd)) / np.max(np.abs(tw_fd)):.1e}; amplitude "
          f"(a/2)G|F₂| = {a / 2 * G * abs(ch16.womersley_F2(al)):.3f} Pa (Poiseuille aG/2 = {a * G / 2:.3f} Pa)")
    Gn = [300.0, 600.0, 250.0 - 200.0j, 80.0j]  # ours, a few harmonics
    tq = np.linspace(0, 2 * np.pi / om, 200)
    Qh = ch16.womersley_flow_rate_harmonics(tq, a, om, rho, mu, Gn)
    print(f"  harmonic gradient (ours, G₀…G₃): mean flow {np.mean(Qh[:-1]) * 1e6:.4f} cm³/s = Poiseuille of G₀ "
          f"{np.pi * Gn[0] * a ** 4 / (8 * mu) * 1e6:.4f} cm³/s; peak {Qh.max() * 1e6:.3f}, min {Qh.min() * 1e6:.3f} cm³/s")
    # independent Crank–Nicolson route (V3)
    errs = []
    with Timer("Crank–Nicolson ×3"):
        for nr, ns in ((50, 100), (100, 200), (200, 400)):
            sol = ch16.womersley_pde_solve(a, om, rho, mu, G, nr=nr, steps_per_period=ns)
            R, TT = np.meshgrid(sol["r"], sol["t"])
            ex = ch16.womersley_velocity(R, TT, a, om, rho, mu, G)
            errs.append(np.max(np.abs(sol["u"] - ex)) / np.max(np.abs(ex)))
    print(f"  Crank–Nicolson (ours) vs Bessel solution: errors {', '.join(f'{e:.2e}' for e in errs)}; observed orders "
          f"{', '.join(f'{v:.2f}' for v in np.log2(np.array(errs[:-1]) / np.array(errs[1:])))}")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
