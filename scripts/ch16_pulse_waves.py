"""Chapter 16, §16.3 — pulse waves in an elastic artery: the linearised tube equations (16.56)–(16.58), the wave
equation (16.61), the Moens–Korteweg speed (16.62)–(16.63), d'Alembert's pulse and p̂ = ρcû (16.65)–(16.68), and a
nonlinear one-dimensional run (ours, MacCormack) with the momentum equation (16.54) corrected (slip S5).

Run: ``.venv/Scripts/python.exe scripts/ch16_pulse_waves.py --no-show``   Figures -> outputs/ch16/pulse_*.png
"""
from __future__ import annotations

import numpy as np
from ch16_common import C_E, C_P, C_REF, C_U, C_WALL, Timer, finish, parse_args, save, setup

from fluidpy import ch16_biofluid_mechanics as ch16


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch16.illustrative_inputs()
    rho, a0, h, E = I["rho"], I["a_aorta"], I["h_aorta"], I["E_aorta"]
    A0 = np.pi * a0 ** 2
    c0 = ch16.moens_korteweg_speed(E, h, rho, a0)
    dA = 1e-10
    dpdA = (ch16.tube_law_hookean(A0 + dA, A0, E, h, a0) - ch16.tube_law_hookean(A0 - dA, A0, E, h, a0)) / (2 * dA)
    print(f"aorta-like tube (ours): a₀ {a0 * 1e3} mm, h {h * 1e3} mm, E {E / 1e6} MPa, ρ {rho} kg/m³")
    print(f"  Moens–Korteweg (16.62) c = √(Eh/2ρa₀) = {c0:.4f} m/s; Bramwell–Hill √(A dp/dA/ρ) with (16.55) = "
          f"{ch16.bramwell_hill_speed(A0, dpdA, rho):.4f} m/s; mid-wall (16.63) = {ch16.moens_korteweg_speed(E, h, rho, a0, midwall=True):.4f} m/s")
    print(f"  ∂p/∂A of (16.58) = Eh/(2a₀A₀) = {E * h / (2 * a0 * A0):.6g} Pa/m², numeric from (16.55) {dpdA:.6g}")
    print(f"  inverse: E from c = {ch16.youngs_modulus_from_pwv(c0, h, rho, a0) / 1e6:.4f} MPa; arrival over 0.5 m: {0.5 / c0 * 1e3:.1f} ms")
    lr = ch16.linearisation_ratios(2.5e3, 12e3, 0.5, c0, Eh_over_a0=E * h / a0)
    print(f"  linearisation (ours inputs): p̂/p̄ = {lr['p_ratio']:.3f}, u/c = {lr['u_over_c']:.3f}, strain p̂/(Eh/a₀) = {lr['strain']:.3f}")
    # d'Alembert pulse, particles barely move
    amp, width = 2.0e3, 0.15
    x = np.linspace(-0.5, 2.0, 1001)
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
    for t, col in ((0.0, C_REF), (0.1, C_U), (0.2, C_P)):
        st = ch16.pulse_state(x, t, c=c0, amplitude=amp, width=width, rho=rho, a0=a0, E=E, h=h, ic="right")
        ax[0].plot(x, st["p"] / ch16.MMHG, color=col, label=f"t = {t:.1f} s")
        ax[1].plot(x, (st["a"] - a0) * 1e6, color=col)
    st = ch16.pulse_state(0.0, 0.0, c=c0, amplitude=amp, width=width, rho=rho, a0=a0, E=E, h=h)
    print(f"  pulse peak {amp / ch16.MMHG:.1f} mmHg: u = p/ρc (16.68) = {st['u']:.4f} m/s "
          f"(velocity_from_pulse {ch16.velocity_from_pulse(amp, rho, c0):.4f}) vs c = {c0:.2f} m/s; wall bulge "
          f"η = {ch16.wall_bulge(amp, a0=a0, E=E, h=h) * 1e6:.1f} µm; particle displacement ∫u dt = "
          f"{st['u'] * width * np.sqrt(np.pi) / c0 * 1e3:.2f} mm while the pulse runs on")
    ax[0].set(xlabel="x [m]", ylabel="p [mmHg]", title=f"right-going pulse at c = {c0:.2f} m/s (16.65)")
    ax[0].legend()
    ax[1].set(xlabel="x [m]", ylabel="η = a − a₀ [µm]", title="wall bulge from the tube law (16.41)")
    Es = np.linspace(0.1e6, 2.0e6, 100)
    for hr, col in ((0.08, C_E), (0.15, C_WALL), (0.25, C_P)):
        ax[2].plot(Es / 1e6, ch16.moens_korteweg_speed(Es, hr * a0, rho, a0), color=col, label=f"h/a₀ = {hr}")
    ax[2].set(xlabel="E [MPa]", ylabel="c [m/s]", title="Moens–Korteweg (16.62): stiffer or thicker → faster")
    ax[2].legend()
    save(fig, out, "pulse_linear")
    # nonlinear 1-D model (ours)
    N = 200 if args.fast else 400
    xg = np.linspace(0.0, 4.0, N, endpoint=False)
    g = np.exp(-((xg - 1.0) / 0.2) ** 2)
    eps = 1e-4
    with Timer("small pulse"):
        r = ch16.tube1d_run(xg, A0 * (1 + eps * g), c0 * eps * g, A0, E, h, a0, rho, 2.0 / c0)
    xpk = xg[np.argmax(r["A"])]
    print(f"  nonlinear model, small pulse (ε = {eps}): peak travelled {xpk - 1.0:.4f} m in 2/c₀ (expected 2.0), "
          f"mass Σ AΔx conserved to {r['mass'] / r['mass0'] - 1:.1e}")
    big = 0.3
    A_init = A0 * (1 + big * g)
    # a right-going simple wave: the left-going Riemann invariant u − 4(c₀ − c(A)) keeps its rest value 0, so
    # u = 4(c₀ − c(A)) with c(A) = c₀(A₀/A)^{1/4} for the Hookean law (16.55) (ours)
    cA = c0 * (A0 / A_init) ** 0.25
    u_init = 4.0 * (c0 - cA)
    lam = u_init + cA  # characteristic speed u + c = 4c₀ − 3c(A) (ours)
    t_break = -1.0 / np.min(np.gradient(lam, xg))
    t_end = 0.8 * t_break  # stop before the front breaks: after it, the captured shock overshoots (scheme limit)
    snaps = (0.0, 0.5 * t_end, t_end)
    r1 = ch16.tube1d_run(xg, A_init, u_init, A0, E, h, a0, rho, t_end, nu_art=0.02, snapshots=snaps)
    r2 = ch16.tube1d_run(xg, A_init, u_init, A0, E, h, a0, rho, t_end, nu_art=0.02, snapshots=snaps,
                         printed_16_54=True)
    print(f"  large pulse (30 % area): simple-wave breaking time t_b = −1/min(dλ/dx) = {t_break:.3f} s; at 0.8 t_b the peak is "
          f"{np.max(r1['p']) / ch16.MMHG:.2f} mmHg (initial {np.max(r1['snapshots']['p'][0]) / ch16.MMHG:.2f}), front slope "
          f"×{np.max(np.abs(np.diff(r1['A']))) / np.max(np.abs(np.diff(A_init))):.2f}; corrected vs printed (16.54) "
          f"differ by {np.max(np.abs(r1['A'] - r2['A'])) / A0:.3f} A₀")
    fig, ax = plt.subplots(figsize=(8.5, 4.4))
    for k, (tt, col) in enumerate(zip(snaps, (C_REF, C_U, C_P))):
        ax.plot(xg, r1["snapshots"]["p"][k] / ch16.MMHG, color=col, label=f"corrected (16.54), t = {tt:.3f} s")
        ax.plot(xg, r2["snapshots"]["p"][k] / ch16.MMHG, color=col, ls=":", lw=1)
    ax.set(xlabel="x [m]", ylabel="p − p_e [mmHg]",
           title="nonlinear 1-D tube model (ours): steepening up to breaking; dotted = printed (16.54)")
    ax.legend(fontsize=8)
    save(fig, out, "pulse_nonlinear")
    # an inlet pressure pulse into a resting tube
    xi = np.linspace(0.0, 3.0, 301)
    pin = lambda t: 2.0e3 * np.exp(-((t - 0.08) / 0.03) ** 2)  # noqa: E731
    ri = ch16.tube1d_run(xi, A0, 0.0, A0, E, h, a0, rho, 0.45, inflow=pin, snapshots=(0.25, 0.4))
    sp = (xi[np.argmax(ri["snapshots"]["p"][1])] - xi[np.argmax(ri["snapshots"]["p"][0])]) / 0.15
    st = ch16.tube_radius_hookean(2.0e3, 0.0, a0, E, h)
    u_pk = 4.0 * (c0 - c0 * (a0 / st) ** 0.5)
    print(f"  inlet pulse of 2 kPa: measured peak speed {sp:.3f} m/s vs c₀ = {c0:.3f} and u + c(A) at the peak "
          f"= {u_pk + c0 * (a0 / st) ** 0.5:.3f} m/s (nonlinear)")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
