"""§7.5 (C09, C10): group velocity — beats of two nearby waves (7.66) with the printed-slip ghost (½Δω x, whose
envelope would not move) (Fig. 7.13 remake); a Gaussian packet evolved by FFT with its envelope moving at c_g (7.68),
measured with the Hilbert envelope and compared with (7.69) and the closed-form chirped packet; the stone-in-a-pond
train with surface tension (Fig. 7.16 remake: calm centre inside c_g,min·t); and the x–t diagram of crests (slope c)
and of the group (slope c_g) (Fig. 7.17 remake).

Run: ``.venv/Scripts/python.exe scripts/ch07_groups.py --no-show [--fast]``
Figures → outputs/ch07/c09_groups.png.
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

    fig, ax = plt.subplots(2, 2, figsize=(13, 8.5))
    # beats (deep water)
    k1, k2 = 1.0, 1.1
    x = np.linspace(0, 130, 3000)
    b0 = ch07.beat_wave(x, 0.0, k1, k2)
    b1 = ch07.beat_wave(x, 20.0, k1, k2)
    bp = ch07.beat_wave(x, 20.0, k1, k2, printed=True)
    plain = np.cos(k1 * x - ch07.omega_gravity(k1) * 20) + np.cos(k2 * x - ch07.omega_gravity(k2) * 20)
    print(f"beats: (7.66) with ½Δω t vs the plain sum: {np.max(np.abs(b1['eta'] - plain)):.1e}; printed ½Δω x: "
          f"{np.max(np.abs(bp['eta'] - plain)):.2f} (wrong)")
    print(f"  c = {b1['c']:.4f} m/s, Δω/Δk = {b1['cg_finite']:.4f} m/s, dω/dk at the mean k = "
          f"{ch07.group_velocity(b1['k']):.4f} m/s")
    a0 = ax[0, 0]
    a0.plot(x, b1["eta"], color=COLORS["blue"], lw=1, label="η at t = 20 s (7.66)")
    a0.plot(x, np.abs(b0["envelope"]), color=COLORS["muted"], lw=1, ls=":", label="envelope at t = 0")
    a0.plot(x, np.abs(b1["envelope"]), color=COLORS["accent"], lw=1.8, label="envelope at t = 20 s (moved c_g t)")
    a0.plot(x, np.abs(bp["envelope"]), color=COLORS["rose"], lw=1, ls="--", label="printed ½Δω x: envelope frozen")
    a0.set_xlabel("x [m]")
    a0.set_ylabel("η/a")
    a0.set_title("two waves k = 1, 1.1 rad/m (deep): beats", fontsize=10)
    a0.set_ylim(-2.2, 3.4)
    a0.legend(fontsize=7, loc="upper right", ncol=2)
    # packet by FFT, envelope speed
    N = 2048 if args.fast else 8192
    L = 1600.0
    xp = np.linspace(-L / 2, L / 2, N, endpoint=False)
    k0, sx = 0.5, 25.0
    eta0 = np.exp(-xp ** 2 / (2 * sx ** 2)) * np.cos(k0 * xp)
    ts = np.linspace(0, 100, 11)
    ev = ch07.linear_evolve(eta0, xp, ts, lambda k: ch07.omega_gravity(k), direction=1)
    env = ch07.envelope(ev, axis=-1)
    xc = np.array([np.sum(xp * e ** 2) / np.sum(e ** 2) for e in env])
    cg_meas = np.polyfit(ts, xc, 1)[0]
    cg = ch07.group_velocity(k0)
    gp = ch07.gaussian_packet(xp, ts[-1], 1.0, k0, sx)
    print(f"packet k₀ = {k0} rad/m, σ_x = {sx} m (deep): envelope centroid speed {cg_meas:.5f} m/s vs c_g (7.69) = {cg:.5f}"
          f" m/s (rel {cg_meas / cg - 1:.1e}); c = {ch07.phase_speed(k0):.4f} m/s")
    print(f"  FFT vs closed-form chirped packet at t = 100 s: max |Δη| = {np.max(np.abs(ev[-1] - gp['eta'])):.2e}")
    # the energy centroid moves at the spectrum-averaged c_g: c_g(k₀)(1 + O((σ_k/k₀)²)), σ_k = 1/σ_x
    print(f"  (σ_k/k₀)² = {(1 / (sx * k0)) ** 2:.4f} sets the size of that difference")
    assert abs(cg_meas / cg - 1) < (1 / (sx * k0)) ** 2
    a1 = ax[0, 1]
    for i, col in ((0, COLORS["muted"]), (5, COLORS["teal"]), (10, COLORS["blue"])):
        a1.plot(xp, ev[i], color=col, lw=0.8)
        a1.plot(xp, env[i], color=COLORS["accent"], lw=1.5)
    a1.set_xlim(-100, 300)
    a1.set_xlabel("x [m]")
    a1.set_ylabel("η/a")
    a1.set_title(f"Gaussian packet at t = 0, 50, 100 s: envelope at c_g = {cg:.2f} m/s, crests at c = "
                 f"{ch07.phase_speed(k0):.2f} m/s", fontsize=9)
    # stone in a pond
    xs = np.linspace(0, 1.2, 1200)
    tp = [0.5, 1.5, 3.0]
    eta_p = ch07.pond_ripples(xs, tp, width=0.01, n_modes=1024 if not args.fast else 512)
    gm = ch07.min_group_velocity(0.0727, 1000.0)
    a2 = ax[1, 0]
    for i, (t, col) in enumerate(zip(tp, (COLORS["teal"], COLORS["blue"], COLORS["accent"]))):
        a2.plot(xs, eta_p[i] + 0.06 * i, color=col, lw=0.9, label=f"t = {t} s")
        a2.axvline(gm["cg_min"] * t, color=col, ls=":", lw=1)
        calm = np.max(np.abs(eta_p[i][xs < 0.6 * gm["cg_min"] * t])) if t > 0.5 else np.nan
        print(f"stone in a pond t = {t} s: c_g,min·t = {100 * gm['cg_min'] * t:.1f} cm; max |η| inside 0.6 c_g,min t ="
              f" {calm:.2e} (hump amplitude 1)")
    a2.set_xlabel("distance from the impact x [m]")
    a2.set_ylabel("η (offset per time)")
    a2.set_title(f"stone in a pond (1 cm hump, σ = 0.0727 N/m): calm inside c_g,min t, c_g,min = "
                 f"{100 * gm['cg_min']:.1f} cm/s", fontsize=9)
    a2.legend(fontsize=8)
    # x–t diagram: crest paths (slope c) and group path (slope c_g)
    a3 = ax[1, 1]
    tt = np.linspace(0, 100, 101)
    evx = ch07.linear_evolve(eta0, xp, tt, lambda k: ch07.omega_gravity(k), direction=1)
    envx = ch07.envelope(evx, axis=-1)
    sel = (xp > -80) & (xp < 400)
    a3.pcolormesh(xp[sel], tt, evx[:, sel], cmap="RdBu_r", shading="auto", vmin=-1, vmax=1)
    a3.plot(xc, ts, color=COLORS["accent"], lw=2.5, label=f"group: dx/dt = c_g = {cg:.2f} m/s")
    c0 = ch07.phase_speed(k0)
    for n in range(-40, 12):
        x0 = 2 * np.pi * n / k0
        xcr = x0 + c0 * tt
        ev_on = np.array([np.interp(xi, xp, envx[j]) for j, xi in enumerate(xcr)])
        a3.plot(np.where(ev_on > 0.25, xcr, np.nan), tt, color=COLORS["orange"], lw=0.8)  # crest inside the group only
    a3.plot([], [], color=COLORS["orange"], label=f"crests: dx/dt = c = {c0:.2f} m/s")
    a3.set_xlim(-80, 400)
    a3.set_ylim(0, 100)
    a3.set_xlabel("x [m]")
    a3.set_ylabel("t [s]")
    a3.set_title("x–t diagram: crests run through the group (Fig. 7.17 remake)", fontsize=10)
    a3.legend(fontsize=8, loc="upper left")
    save(fig, out, "c09_groups")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
