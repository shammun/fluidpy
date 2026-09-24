"""§7.1 (C01): the wave vocabulary — a sinusoid's k, λ, ω, T, c (7.1)–(7.4); a crest tracked by hand moves at ω/k; the
2-D plane wave (7.5) with crest spacing λ = 2π/K along K but 2π/k and 2π/l along the axes (Fig. 7.1 remake) and the
trace velocities (larger than c, not components of c); the Doppler shift (7.9).

Run: ``.venv/Scripts/python.exe scripts/ch07_wave_basics.py --no-show``
Figures → outputs/ch07/c01_wave_basics.png.
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

    wp = ch07.wave_parameters(lam=100.0, T=8.0)
    print("λ = 100 m, T = 8 s → " + ", ".join(f"{k} = {v:.4g}" for k, v in wp.items()))
    # from-scratch crest tracking: argmax of η near the n = 0 crest at several times
    x = np.linspace(-50, 150, 200001)
    ts = np.linspace(0, 4, 5)
    xc = [x[np.argmax(ch07.sinusoid(x, t, 1.0, wp["k"], wp["omega"]) - 1e-9 * np.abs(x - wp["c"] * t))] for t in ts]
    slope = np.polyfit(ts, xc, 1)[0]
    print(f"  crest tracked by hand moves at {slope:.5f} m/s; ω/k = {wp['c']:.5f} m/s (rel diff {slope / wp['c'] - 1:.1e})")
    assert abs(slope / wp["c"] - 1) < 1e-3
    K = np.array([1.0, 1.0])
    om = 1.0
    Kmag = np.linalg.norm(K)
    tr = ch07.trace_velocities(K, om)
    print(f"k = l = 1 rad/m: K = {Kmag:.4f}, λ = 2π/K = {2 * np.pi / Kmag:.3f} m (crests 2π/k = {2 * np.pi:.3f} m apart"
          f" along x); c = ω/K = {om / Kmag:.4f} m/s, trace c_x = c_y = {tr[0]:.4f} m/s")
    print(f"  1/c² = {Kmag ** 2 / om ** 2:.4f} = Σ 1/c_i² = {sum(1 / c ** 2 for c in tr):.4f}")
    w0 = ch07.doppler_frequency(2 * np.pi / 8.0, [-1.0, 0.0], [ch07.wavenumber_from_omega(2 * np.pi / 8.0), 0.0])
    print(f"Doppler: T = 8 s deep-water swell against a 1 m/s current: ω = {2 * np.pi / 8:.4f} → ω₀ = {w0:.4f} rad/s"
          f" (observed period {2 * np.pi / w0:.3f} s)")

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
    X, Y = np.meshgrid(np.linspace(0, 14, 400), np.linspace(0, 14, 400), indexing="xy")
    eta = ch07.plane_wave(np.stack([X, Y], axis=-1), K, om, 1.0)
    ax[0].contourf(X, Y, eta, levels=21, cmap="RdBu_r", alpha=0.6)
    ax[0].contour(X, Y, eta, levels=[0.999], colors=COLORS["orange"], linewidths=1.5)
    ax[0].annotate("", xy=(7 + 1.5 / Kmag, 7 + 1.5 / Kmag), xytext=(7, 7),
                   arrowprops=dict(arrowstyle="->", color=COLORS["ink"], lw=2))
    ax[0].text(7.9, 8.4, "K", fontsize=12)
    ax[0].set_title("plane wave η = a cos(K·x − ωt) (7.5): crests (orange) perpendicular to K", fontsize=10)
    ax[0].set_xlabel("x [m]")
    ax[0].set_ylabel("y [m]")
    ax[0].set_aspect("equal")
    ax[0].text(0.3, 0.4, f"spacing along K: 2π/K = {2 * np.pi / Kmag:.2f} m\nalong x: 2π/k = {2 * np.pi:.2f} m",
               fontsize=8, bbox=dict(fc="white", alpha=0.8))
    xx = np.linspace(0, 300, 1000)
    for t, col in zip((0.0, 2.0, 4.0), (COLORS["blue"], COLORS["teal"], COLORS["accent"])):
        ax[1].plot(xx, ch07.sinusoid(xx, t, 1.0, wp["k"], wp["omega"]), color=col, label=f"t = {t:.0f} s")
        ax[1].plot(ch07.crest_positions(t, wp["k"], wp["omega"], 0), 1.0, "o", color=COLORS["orange"])
    ax[1].set_xlabel("x [m]")
    ax[1].set_ylabel("η/a")
    ax[1].set_title(f"λ = 100 m, T = 8 s: crest n = 0 moves at c = ω/k = {wp['c']:.2f} m/s (7.4)", fontsize=10)
    ax[1].legend(fontsize=8)
    save(fig, out, "c01_wave_basics")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
