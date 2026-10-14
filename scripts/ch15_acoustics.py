"""Chapter 15, §15.1–§15.2 — sound as the linear limit: p′ = c²ρ′, d'Alembert's two travelling shapes, waves in a
stream, and steepening.

* how small is "small": condensation p′/ρ₀c² (15.7) and the size of the dropped quadratic term; Newton vs Laplace;
* d'Alembert (15.15) with the particle velocity (15.16): a pulse released from rest splits in two; the air moves with
  a compression and against an expansion; (15.16) against the time quadrature (15.14);
* waves in a stream (15.17): the two characteristic families x = x₀ + (U ± c)t for U < c and U > c;
* steepening (ours): the simple wave x = x₀ + (c₀ + ½(γ + 1)u₀)t up to the breaking time, checked by a MacCormack
  run of the Euler equations (ours, optional).

Run: ``.venv/Scripts/python.exe scripts/ch15_acoustics.py --no-show``   Figures -> outputs/ch15/acoustics_*.png
"""
from __future__ import annotations

import numpy as np
from ch15_common import C_ENT, C_PRES, C_REF, C_SUB, C_SUP, COLORS, Timer, finish, parse_args, save, setup

from fluidpy import ch15_compressible_flow as ch15


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch15.illustrative_inputs()
    g, R, T, p0 = I["gamma"], I["R"], I["T_air"], I["p_air"]
    gas = ch15.PerfectGas(g, R)
    c = gas.c(T)
    rho0 = p0 / (R * T)
    print(f"air at T = {T} K, p = {p0} Pa: ρ₀ = {rho0:.4f} kg/m³, c = √(γRT) = {c:.2f} m/s (15.6); Newton's √(RT) = "
          f"{ch15.newton_sound_speed(T, R):.2f} m/s, ratio {ch15.newton_sound_speed(T, R) / c:.4f} = 1/√γ = {1 / np.sqrt(g):.4f}")
    print(f"  ρ₀c² = γp₀ = {rho0 * c ** 2:.1f} Pa; impedance ρ₀c = {ch15.acoustic_impedance(rho0, c):.1f} Pa s/m")
    for pp in (0.02, 1.0, 20.0, 2000.0):
        print(f"  p′ = {pp:7.2f} Pa ({ch15.sound_pressure_level(pp):5.1f} dB): condensation (15.7) = {ch15.condensation(pp, rho0, c):.2e}, "
              f"dropped quadratic term / linear term = {ch15.linearisation_error(pp, rho0, c, g):.2e}, ρ′ = {ch15.acoustic_density(pp, c):.2e} kg/m³, "
              f"T′ = {ch15.acoustic_temperature(pp, p0, T, g):.2e} K, u′ = p′/ρ₀c = {pp / (rho0 * c):.2e} m/s")
    for U in (0.0, 100.0, c, 510.0):
        print(f"  regime at U = {U:6.1f} m/s: M = {ch15.mach_from_speed(U, T, gas):.3f} → {ch15.mach_regime(ch15.mach_from_speed(U, T, gas))}; "
              f"wave speeds (U + c, U − c) = ({ch15.wave_speeds(U, c)[0]:.1f}, {ch15.wave_speeds(U, c)[1]:.1f}) m/s; upstream signal possible: {ch15.upstream_signal_possible(U, c)}")

    # --- d'Alembert: a pulse released from rest -------------------------------------------------------------------------
    F = ch15.pulse("gaussian", I["pulse_width"], I["pulse_amplitude"])
    half = lambda s: 0.5 * F(s)  # noqa: E731  (released from rest: f = g = F/2, Exercise 15.1c with G = 0)
    x = np.linspace(-6.0, 6.0, 1201)
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
    for t_, col in zip((0.0, 0.004, 0.010), (C_REF, C_SUB, C_SUP)):
        ax[0].plot(x, ch15.dalembert(x, t_, c, half, half), color=col, label=f"t = {1e3 * t_:.0f} ms")
        ax[1].plot(x, 1e3 * ch15.plane_wave_velocity(x, t_, c, rho0, half, half), color=col, label=f"t = {1e3 * t_:.0f} ms")
    ax[0].set(xlabel="x [m]", ylabel="p′ [Pa]", title="d'Alembert (15.15): the pulse splits into f(x − ct) + g(x + ct)")
    ax[1].set(xlabel="x [m]", ylabel="u′ [mm/s]", title="particle velocity (15.16): with the wave in a compression")
    ax[0].legend(fontsize=8)
    ax[1].legend(fontsize=8)
    xq, tq = 2.0, 0.006
    uq = ch15.velocity_from_pressure(lambda a, b: ch15.dalembert(a, b, c, half, half), xq, tq, rho0, t0=0.0)
    print(f"d'Alembert pulse (ours: Gaussian, {I['pulse_amplitude']} Pa, width {I['pulse_width']} m): at x = {xq} m, t = {1e3 * tq:.0f} ms "
          f"p′ = {ch15.dalembert(xq, tq, c, half, half):.4f} Pa; u′ by (15.16) = {ch15.plane_wave_velocity(xq, tq, c, rho0, half, half):.6e} m/s, "
          f"by the time integral (15.14) = {uq:.6e} m/s; initial-value form = {ch15.dalembert_initial_value(xq, tq, c, F):.4f} Pa")
    E = [ch15.acoustic_energy(ch15.dalembert(x, t_, c, half, half), ch15.plane_wave_velocity(x, t_, c, rho0, half, half), x, rho0, c) for t_ in (0.004, 0.010)]
    print(f"  acoustic energy per area at 4 ms and 10 ms: {E[0]:.6e}, {E[1]:.6e} J/m² (constant once the halves have separated); "
          f"wave-equation residual (15.13) of the solution: {ch15.wave_equation_residual(lambda a, b: ch15.dalembert(a, b, c, half, half), xq, tq, c, h=1e-3):.2e} Pa/m² "
          f"(scale p′/w² = {I['pulse_amplitude'] / I['pulse_width'] ** 2:.0f})")

    # --- x–t diagram in a stream -------------------------------------------------------------------------------------------
    tt = np.linspace(0.0, 0.01, 2)
    for U, ls, lab in ((I["U_stream_sub"], "-", "U < c"), (I["U_stream_sup"], "--", "U > c")):
        ch_ = ch15.xt_characteristics(U, c, [0.0], tt)
        ax[2].plot(ch_["x_plus"][0], 1e3 * tt, color=C_SUP, ls=ls, label=f"{lab}: dx/dt = U + c = {U + c:.0f} m/s")
        ax[2].plot(ch_["x_minus"][0], 1e3 * tt, color=C_SUB, ls=ls, label=f"{lab}: dx/dt = U − c = {U - c:.0f} m/s")
    ax[2].axvline(0.0, color=C_REF, lw=0.8)
    ax[2].set(xlabel="x [m]", ylabel="t [ms]", title="waves in a stream (15.17): for U > c nothing travels upstream")
    ax[2].legend(fontsize=7)
    save(fig, out, "acoustics_dalembert_stream")
    res = ch15.convected_wave_residual(lambda a, b: ch15.dalembert_in_stream(a, b, c, I["U_stream_sup"], half, half), 3.0, 0.004, c, I["U_stream_sup"], h=1e-3)
    print(f"  convected wave equation (15.12) residual of (15.17) at U = {I['U_stream_sup']} m/s: {res:.2e} Pa/m²")

    # --- steepening ----------------------------------------------------------------------------------------------------------
    lam, ua = I["steep_wavelength"], I["steep_u"]
    n = 201 if args.fast else 401
    xs = np.linspace(0.0, lam, n)[:-1]
    u0 = ua * np.sin(2 * np.pi * xs / lam)
    tb = ch15.breaking_time(u0, xs, g)
    tb_exact = ch15.breaking_time(None, None, g, max_compression_rate=2 * np.pi * ua / lam)
    print(f"steepening (ours): u₀ = {ua} sin(2πx/λ) m/s, λ = {lam} m, M_u = u/c = {ua / c:.3f}, p′ ≈ ρ₀c u = {rho0 * c * ua:.0f} Pa "
          f"(condensation {ua / c:.3f}): breaking time 2/((γ + 1)·2πu/λ) = {1e3 * tb_exact:.4f} ms (from samples {1e3 * tb:.4f} ms) "
          f"= {tb_exact * c / lam:.3f} wavelengths of travel; wavelet speeds from {ch15.simple_wave_speed(-ua, c, g):.1f} to {ch15.simple_wave_speed(ua, c, g):.1f} m/s")
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    for fr, col in zip((0.0, 0.5, 1.0, 1.4), (C_REF, C_SUB, C_PRES, C_ENT)):
        sw = ch15.simple_wave_evolve(u0, xs, fr * tb, c, g)
        xr = sw["x"] - c * fr * tb  # frame moving at c₀
        ax[0].plot(xr / lam, sw["u"], color=col, label=f"t = {fr:.1f} t_b" + ("" if sw["valid"] else "  (multivalued: a shock has formed)"))
    ax[0].set(xlabel="(x − c₀t)/λ", ylabel="u [m/s]", title="simple wave (ours): compressions steepen, expansions flatten")
    ax[0].legend(fontsize=8)
    with Timer("Euler run (MacCormack, ours)"):
        sw = ch15.simple_wave_evolve(u0, xs, 0.6 * tb, c, g)
        cc = c + 0.5 * (g - 1.0) * u0
        e = ch15.euler1d_run(rho0 * (cc / c) ** (2 / (g - 1)), u0, p0 * (cc / c) ** (2 * g / (g - 1)), xs, 0.6 * tb, g, eps_av=0.0)
    xw = np.mod(sw["x"], lam)
    o = np.argsort(xw)
    uref = np.interp(xs, xw[o], sw["u"][o], period=lam)
    ax[1].plot(xs / lam, e["u"], color=C_SUP, label="Euler equations, MacCormack (ours)")
    ax[1].plot(xs / lam, uref, color=COLORS["ink"], ls="--", lw=1, label="simple wave by characteristics")
    ax[1].set(xlabel="x/λ (periodic)", ylabel="u [m/s]", title="at t = 0.6 t_b: the two routes agree")
    ax[1].legend(fontsize=8)
    save(fig, out, "acoustics_steepening")
    print(f"  at 0.6 t_b: max |u_Euler − u_characteristics| = {np.max(np.abs(e['u'] - uref)):.3f} m/s ({np.max(np.abs(e['u'] - uref)) / ua:.2%} of the amplitude), "
          f"{e['n_steps']} steps on {xs.size} nodes; mass drift {e['mass'] - e['mass0']:.1e}, energy drift {(e['energy'] - e['energy0']) / e['energy0']:.1e}")
    r = ch15.linearised_acoustics_sympy()
    print(f"sympy: linear continuity (15.8) = {r['continuity_linear']}; dropped: {r['continuity_dropped']}; cross-differentiation residual {r['cross_residual']}, "
          f"(15.17) in (15.12): {r['dalembert_residual']}")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
