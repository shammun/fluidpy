"""§7.4 (C08): standing waves and seiches — the sum of a right- and a left-going wave (7.61) pins its nodes; its
streamlines (7.62) (Fig. 7.11 remake); walls at x = 0, L allow λ = 2L/(n + 1) (7.64) with the lake frequencies (7.65)
(Fig. 7.12 remake: u(x) of modes n = 0, 1); a 50 km × 100 m lake's periods with the shallow-water estimate.

Run: ``.venv/Scripts/python.exe scripts/ch07_standing.py --no-show``
Figures → outputs/ch07/c08_standing.png.
"""
from __future__ import annotations

import numpy as np

from ch07_drawings import parse_args, save, setup, wave_surface

from fluidpy import ch07_gravity_waves as ch07
from fluidpy.core.style import COLORS


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    a, k, H = 0.1, 1.0, 1.0
    x = np.linspace(0, 2 * np.pi, 9)
    z = np.linspace(-1, 0, 9)
    s = ch07.standing_wave_fields(x, z, 0.37, a, k, H)
    r1 = ch07.wave_fields(x, z, 0.37, a, k, H)
    r2 = ch07.wave_fields(x, z, 0.37, a, k, H, direction=-1)
    d = max(np.max(np.abs(s[n] - (r1[n] + r2[n]))) for n in ("eta", "psi", "u", "w"))
    print(f"standing wave = sum of the two travelling waves: max |Δ| = {d:.1e}")
    nodes = np.array([np.pi / 2, 3 * np.pi / 2])
    print(f"η at the nodes kx = π/2, 3π/2 over a period: max |η| = "
          f"{max(np.max(np.abs(ch07.standing_wave_fields(nodes, 0, t, a, k, H)['eta'])) for t in np.linspace(0, 3, 13)):.1e}")
    L, Hl = 50e3, 100.0
    for n in range(4):
        st = ch07.seiche_state(L, Hl, n)
        print(f"lake L = 50 km, H = 100 m, mode n = {n}: λ = {st['lam'] / 1e3:.1f} km, T = {st['T_min']:.2f} min"
              f" (shallow 2L/((n+1)√(gH)) = {st['T_shallow'] / 60:.2f} min, error {100 * st['shallow_error']:+.4f} %)")
    xw = np.linspace(0, L, 5)
    for n in range(3):
        kk = float(ch07.seiche_modes(L, Hl, n)["k"])
        uw = ch07.standing_wave_fields(xw[[0, -1]], -50.0, 100.0, 0.1, kk, Hl)["u"]
        assert np.max(np.abs(uw)) < 1e-12 * 0.1 * 1e3

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))
    X, Z = np.meshgrid(np.linspace(0, 2 * np.pi, 300), np.linspace(-H, 0, 120), indexing="xy")
    t = 0.25 * 2 * np.pi / ch07.omega_gravity(k, H)  # quarter period: flow strongest
    f = ch07.standing_wave_fields(X, Z, t, a, k, H)
    m = np.abs(f["psi"]).max()
    ax[0].contour(X, Z, f["psi"], levels=np.linspace(-m, m, 15), colors=COLORS["ink"], linewidths=0.8,
                  negative_linestyles="solid")
    for tt, ls in ((0.0, "-"), (0.5 * 2 * np.pi / ch07.omega_gravity(k, H), "--")):
        wave_surface(ax[0], X[0], ch07.standing_wave_fields(X[0], 0.0, tt, a, k, H)["eta"], ls=ls)
    ax[0].plot(nodes, [0, 0], "o", color=COLORS["orange"], label="nodes (η = 0 always)")
    ax[0].set_title("standing wave: ψ (7.62) at a quarter period, surface at t = 0 and T/2", fontsize=10)
    ax[0].set_xlabel("x [m]")
    ax[0].set_ylabel("z [m]")
    ax[0].legend(fontsize=8)
    xs = np.linspace(0, L, 300)
    for n, col in ((0, COLORS["accent"]), (1, COLORS["teal"]), (2, COLORS["orange"])):
        kk = float(ch07.seiche_modes(L, Hl, n)["k"])
        ax[1].plot(xs / 1e3, np.sin(kk * xs), color=col, label=f"n = {n}: λ = 2L/{n + 1} (7.64)")
    ax[1].axhline(0, color=COLORS["muted"], lw=0.8)
    ax[1].set_xlabel("x [km]")
    ax[1].set_ylabel("u amplitude (∝ sin kx, (7.63))")
    ax[1].set_title("seiche modes: u = 0 at both walls (Fig. 7.12 remake)", fontsize=10)
    ax[1].legend(fontsize=8)
    save(fig, out, "c08_standing")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
