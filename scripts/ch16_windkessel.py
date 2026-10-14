"""Chapter 16, §16.3 — the Windkessel: K dp/dt + p/R = Q(t) (16.49) solved three ways — the integrating-factor
quadrature (16.52), solve_ivp, and our closed-form periodic beat for half-sine ejection — the transient from the first
beat, the diastolic decay with time constant RK, and how compliance sets the pulse pressure (Fig. 16.5's smoothing).

Run: ``.venv/Scripts/python.exe scripts/ch16_windkessel.py --no-show``   Figures -> outputs/ch16/windkessel_*.png
"""
from __future__ import annotations

import numpy as np
from ch16_common import C_E, C_P, C_REF, C_U, C_WALL, Timer, finish, parse_args, save, setup

from fluidpy import ch16_biofluid_mechanics as ch16


def rk4_windkessel(t, Q_fn, R, K, p0):
    """From-scratch RK4 for K dp/dt = Q − p/R (the notebook's transparent version)."""
    p = np.empty_like(t)
    p[0] = p0
    f = lambda tt, pp: (float(Q_fn(tt)) - pp / R) / K  # noqa: E731
    for i in range(len(t) - 1):
        h = t[i + 1] - t[i]
        k1 = f(t[i], p[i])
        k2 = f(t[i] + h / 2, p[i] + h / 2 * k1)
        k3 = f(t[i] + h / 2, p[i] + h / 2 * k2)
        k4 = f(t[i] + h, p[i] + h * k3)
        p[i + 1] = p[i] + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return p


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch16.illustrative_inputs()
    mm = ch16.MMHG
    R, K, HR, SV, sh = I["R_periph"], I["K_compliance"], I["HR"], I["SV"], I["systole_share"]
    T = 60.0 / HR
    Qf = ch16.HalfSineEjection(SV, T, sh * T)
    kw = dict(SV=SV, HR=HR, systole_share=sh, R=R, K=K)
    p0_end = ch16.windkessel_periodic_p0(**kw)
    print(f"Windkessel (ours): R = {R:.3g} Pa s/m³, K = {K:.3g} m³/Pa ({K * 1e6 * mm:.2f} mL/mmHg), RK = {R * K:.3f} s, "
          f"beat {T:.3f} s, systole {sh * T:.3f} s, SV {SV * 1e6:.0f} mL")
    print(f"  periodic end-diastolic pressure (closed form) p₀ = {p0_end / mm:.3f} mmHg; by quadrature "
          f"{ch16.windkessel_periodic_p0_general(Qf, R, K, T) / mm:.3f} mmHg")
    nb = 4 if args.fast else 8
    t = np.linspace(0.0, nb * T, nb * 200 + 1)
    with Timer("three routes"):
        p_quad = ch16.windkessel_pressure(t, Qf, R, K, 0.0)
        p_ode = ch16.windkessel_ode(t, Qf, R, K, 0.0)
        p_rk4 = rk4_windkessel(t, Qf, R, K, 0.0)
    print(f"  from rest: quadrature (16.52) vs solve_ivp {np.max(np.abs(p_quad - p_ode)) / np.max(p_quad):.1e}, vs "
          f"from-scratch RK4 (Δt = T/200) {np.max(np.abs(p_quad - p_rk4)) / np.max(p_quad):.1e} (relative)")
    p_per = ch16.windkessel_periodic_state(t, **kw)
    gap = [np.max(np.abs(p_quad[k * 200:(k + 1) * 200] - p_per[k * 200:(k + 1) * 200])) / mm for k in range(nb)]
    print("  distance to the periodic beat, per beat [mmHg]: " + ", ".join(f"{g:.2f}" for g in gap))
    s = ch16.windkessel_summary(**kw)
    from scipy.integrate import quad

    vol = quad(Qf, 0.0, T, points=[sh * T])[0] - quad(lambda s_: ch16.windkessel_periodic_state(s_, **kw) / R, 0.0, T,
                                                         points=[sh * T], epsrel=1e-13)[0]
    print(f"  periodic beat: p_S = {s['pS'] / mm:.2f}, p_D = {s['pD'] / mm:.2f} mmHg, pulse pressure {s['pulse_pressure'] / mm:.2f} mmHg, "
          f"peak at t = {s['t_peak']:.3f} s; mean = R·Q̄ = {s['p_mean'] / mm:.2f} mmHg (sampled {s['p_mean_sampled'] / mm:.2f}); "
          f"⅓ rule {s['map_third'] / mm:.2f} mmHg; volume in − out over a beat {vol:.2e} m³")
    fig, ax = plt.subplots(1, 2, figsize=(14, 4.6))
    ax[0].plot(t, Qf(t) * 1e6, color=C_U, lw=1, label="ejection Q(t) [cm³/s]")
    ax0b = ax[0].twinx()
    ax0b.plot(t, p_quad / mm, color=C_P, label="aortic p(t) from rest (16.52)")
    ax0b.plot(t, p_per / mm, color=C_REF, ls=":", label="periodic beat (closed form, ours)")
    ax0b.plot(t, R * Qf(t) / mm, color=C_E, lw=0.8, alpha=0.6, label="K → 0: p = QR")
    ax[0].set(xlabel="t [s]", ylabel="Q [cm³/s]", title="Windkessel (16.49): the transient settles in a few RK")
    ax0b.set_ylabel("p [mmHg]")
    ax0b.legend(loc="upper right", fontsize=8)
    Ks = np.geomspace(0.25, 4.0, 25) * K
    pps = [ch16.windkessel_summary(SV=SV, HR=HR, systole_share=sh, R=R, K=k)["pulse_pressure"] / mm for k in Ks]
    ax[1].semilogx(Ks * 1e6 * mm, pps, color=C_WALL)
    ax[1].axvline(K * 1e6 * mm, color=C_REF, ls=":")
    ax[1].set(xlabel="compliance K [mL/mmHg]", ylabel="pulse pressure p_S − p_D [mmHg]",
              title="a stiffer aorta (smaller K): larger pulse pressure, same mean")
    save(fig, out, "windkessel")
    # diastolic decay check: p = p_s e^{−(t − T_s)/RK}
    td = np.linspace(sh * T, T, 50)
    pd_ = ch16.windkessel_periodic_state(td, **kw)
    slope = np.polyfit(td, np.log(pd_), 1)[0]
    print(f"  diastole: d ln p/dt = {slope:.6f} 1/s = −1/RK = {-1 / (R * K):.6f} 1/s")
    tm = ch16.windkessel_terms(0.5 * sh * T, **kw)
    print(f"  mid-systole: inflow {tm['Q'] * 1e6:.1f} cm³/s = storage {tm['storage'] * 1e6:.1f} + run-off {tm['runoff'] * 1e6:.1f} cm³/s")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
