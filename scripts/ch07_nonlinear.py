"""§7.6 (C11, C12): nonlinear waves — simple-wave steepening of a shallow-water hump to breaking (Fig. 7.19; the Riemann
simple wave is our extension), the Stokes wave (7.82)–(7.83) (Fig. 7.21 remake), Stokes drift from exact path lines
against (7.85)/(7.86) with a dyed vertical line leaning forward (Fig. 7.22 remake), cnoidal and solitary waves (7.88)
(Fig. 7.23 remake), and a KdV (7.87) run of a hump splitting into solitons (IF-RK4 pseudo-spectral, our scheme) with its
invariants.

Run: ``.venv/Scripts/python.exe scripts/ch07_nonlinear.py --no-show [--fast]``
Figures → outputs/ch07/c11_c12_nonlinear.png.
"""
from __future__ import annotations

import numpy as np

from ch07_drawings import parse_args, save, setup

from fluidpy import ch07_gravity_waves as ch07
from fluidpy.core.style import COLORS


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 3, figsize=(16, 8.6))
    # (a) simple-wave steepening
    H, a = 1.0, 0.1
    lam = 20.0
    x = np.linspace(0, lam, 400, endpoint=False)
    eta0 = a * np.sin(2 * np.pi * x / lam)
    sw = ch07.simple_wave_evolve(eta0, x, 0.0, H)
    tb = sw["t_break"]
    tb_formula = 2 * H / (3 * a * (2 * np.pi / lam) * np.sqrt(ch07.G_BOOK * H))
    print(f"simple wave a = {a} m, λ = {lam} m on H = {H} m: breaking time {tb:.3f} s (sine formula 2H/(3akc₀) = "
          f"{tb_formula:.3f} s); Ursell aλ²/H³ = {ch07.ursell_number(a, lam, H):.0f}")
    assert abs(tb / tb_formula - 1) < 0.1 * a / H  # the formula is leading order in a/H
    for t, col in ((0.0, COLORS["muted"]), (0.5 * tb, COLORS["teal"]), (tb, COLORS["blue"]), (1.4 * tb, COLORS["rose"])):
        s = ch07.simple_wave_evolve(eta0, x, t, H)
        ax[0, 0].plot(s["x_points"], s["eta"], color=col, label=f"t = {t / tb:.1f} t_b")
    ax[0, 0].set_xlabel("x [m]")
    ax[0, 0].set_ylabel("η [m]")
    ax[0, 0].set_title("crests outrun troughs: steepening, then multivalued (Fig. 7.19)", fontsize=9)
    ax[0, 0].legend(fontsize=8)
    # (b) Stokes wave
    k = 1.0
    xs = np.linspace(-2 * np.pi, 2 * np.pi, 600)
    for ka, col in ((0.1, COLORS["muted"]), (0.3, COLORS["blue"])):
        ax[0, 1].plot(xs, ch07.stokes_wave_profile(xs, 0.0, ka / k, k) / (ka / k), color=col, label=f"ka = {ka} (7.82)")
    ax[0, 1].plot(xs, np.cos(k * xs), color=COLORS["muted"], ls="--", lw=1, label="linear cos")
    print(f"Stokes wave ka = 0.3: crest {ch07.stokes_wave_profile(0.0, 0.0, 0.3, 1.0):.4f} m, trough "
          f"{ch07.stokes_wave_profile(np.pi, 0.0, 0.3, 1.0):.4f} m; c/√(g/k) = "
          f"{ch07.stokes_wave_speed(1.0, 0.3) / np.sqrt(ch07.G_BOOK):.4f} (7.83); limiting H/λ = {ch07.STOKES_LIMIT_STEEPNESS}")
    se = ch07.stokes_expansion_sympy()
    print(f"  Exercise 7.2 in sympy: third-order α = {se['alpha']}, γ = {se['gamma']} (β = {se['third_order']['beta']},"
          f" δ = {se['third_order']['delta']}); literal set-up (fixed potential amplitude) gives γ = "
          f"{se['exercise_literal']['gamma']}")
    ax[0, 1].set_xlabel("x [m]")
    ax[0, 1].set_ylabel("η/a")
    ax[0, 1].set_title("Stokes wave: peaked crests, flat troughs (Fig. 7.21 remake)", fontsize=9)
    ax[0, 1].legend(fontsize=8)
    # (c) Stokes drift and dyed line (deep water)
    a, k = 0.15, 1.0
    om = ch07.omega_gravity(k)
    T = 2 * np.pi / om
    periods = 6 if args.fast else 10
    for z0 in (-0.2, -0.6):
        num = ch07.stokes_drift_numeric(z0, a, k, periods=periods)
        print(f"drift ka = {a * k}, z₀ = {z0} m: exact path lines {num * 100:.3f} cm/s vs (7.85) "
              f"{ch07.stokes_drift(z0, a, k) * 100:.3f} cm/s")
    p = ch07.particle_path(0.0, -0.2, np.linspace(0, 3 * T, 600), a, k)
    ax[0, 2].plot(p["x"], p["z"], color=COLORS["teal"], lw=1.2, label="one particle, 3 periods (exact (7.33))")
    z0s = np.linspace(-1.5, 0.0, 16)
    for n, col in ((0, COLORS["muted"]), (5, COLORS["blue"]), (10, COLORS["accent"])):
        d = ch07.dyed_line(z0s, n * T, a, k)
        ax[0, 2].plot(d["x"], d["z"], "o-", ms=3, color=col, lw=1, label=f"dyed line after {n} periods")
    ax[0, 2].set_xlabel("x [m]")
    ax[0, 2].set_ylabel("z [m]")
    ax[0, 2].set_title(f"Stokes drift: open orbits, the dye leans forward (ka = {a * k})", fontsize=9)
    ax[0, 2].legend(fontsize=7)
    # (d) cnoidal and solitary waves
    Hs = 1.0
    xx = np.linspace(-40, 40, 1600)
    cw = ch07.cnoidal_wave(xx, 0.0, Hs, 0.2, 0.99)
    ax[1, 0].plot(xx, cw["eta"], color=COLORS["blue"], label=f"cnoidal, m = 0.99, λ = {cw['wavelength']:.1f} m")
    ax[1, 0].plot(xx, ch07.solitary_wave(xx, 0.0, 0.2, Hs), color=COLORS["accent"],
                  label=f"solitary a = 0.2 m, c = {ch07.solitary_wave_speed(0.2, Hs):.3f} m/s (7.88)")
    ax[1, 0].set_xlabel("x [m]")
    ax[1, 0].set_ylabel("η [m]")
    ax[1, 0].set_title("permanent forms of KdV on H = 1 m (Fig. 7.23 remake)", fontsize=9)
    ax[1, 0].legend(fontsize=7)
    print(f"solitary a = 0.2 m on H = 1 m: c = {ch07.solitary_wave_speed(0.2, 1.0):.4f} m/s = "
          f"{ch07.solitary_wave_speed(0.2, 1.0) / np.sqrt(ch07.G_BOOK):.2f}√(gH); sympy residual {ch07.kdv_residual_sympy()['residual']}")
    # (e) KdV: hump splitting into solitons
    N = 256 if args.fast else 512
    L = 200.0
    xk = np.linspace(0, L, N, endpoint=False)
    dxk = xk[1] - xk[0]
    eta_k0 = 0.15 * np.exp(-((xk - 30.0) / 6.0) ** 2)
    tout = np.linspace(0, 120, 4)
    sol = ch07.kdv_solve(eta_k0, xk, tout, 1.0)
    c0 = np.sqrt(ch07.G_BOOK)
    for i, col in enumerate((COLORS["muted"], COLORS["teal"], COLORS["blue"], COLORS["accent"])):
        shift = int(round(c0 * tout[i] / dxk))  # co-moving frame at c₀ (periodic): roll by whole grid cells
        ax[1, 1].plot(xk, np.roll(sol["eta"][i], -shift) + 0.1 * i, color=col, lw=1.2, label=f"t = {tout[i]:.0f} s")
    ax[1, 1].set_xlim(0, 110)
    print(f"  tallest soliton at t = {tout[-1]:.0f} s: {sol['eta'][-1].max():.3f} m (hump 0.15 m)")
    ax[1, 1].set_xlabel("x − c₀t [m] (periodic domain)")
    ax[1, 1].set_ylabel("η (offset per time) [m]")
    ax[1, 1].set_title("KdV (7.87): a hump sorts into solitons, tallest first", fontsize=9)
    ax[1, 1].legend(fontsize=7)
    drift = {n: np.max(np.abs(sol[n] / sol[n][0] - 1)) for n in ("mass", "momentum", "energy")}
    print(f"KdV run N = {N}, {sol['steps']} steps (dt = {sol['dt']:.4f} s): invariant drifts " +
          ", ".join(f"{k_} {v:.1e}" for k_, v in drift.items()))
    assert max(drift.values()) < 1e-6
    # (f) KdV linear phase speed vs (7.29)
    kH = np.logspace(-2, 0, 100)
    ax[1, 2].loglog(kH, np.abs(ch07.kdv_linear_phase_speed(kH, 1.0) / ch07.phase_speed(kH, 1.0) - 1),
                    color=COLORS["accent"], label="|c_KdV/c − 1|")
    ax[1, 2].loglog(kH, kH ** 4 * 19 / 360, color=COLORS["muted"], ls="--", lw=1, label="19(kH)⁴/360 (next Taylor term)")
    ax[1, 2].set_xlabel("kH")
    ax[1, 2].set_ylabel("relative error")
    ax[1, 2].set_title("linearised KdV c₀(1 − k²H²/6) vs (7.29): error ∝ (kH)⁴", fontsize=9)
    ax[1, 2].legend(fontsize=8)
    save(fig, out, "c11_c12_nonlinear")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
