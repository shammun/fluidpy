"""§5.6 rotating frame (C09, C11, R18, N39, N54): the fluid column conserving (ζ + f)/h over a ridge, a trough and a
slope; the ring of air carried poleward (5.33); the rotating Kelvin scenario — Γ changes while Γ_a is conserved; the
planetary stretching term 2Ω∂w/∂z and the from-scratch vector area ½∮x × dx.

Run: ``.venv/Scripts/python.exe scripts/ch05_rotating_column.py --no-show``
Figures → outputs/ch05/c11_column.png, c11_absolute_circulation.png.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from fluidpy import ch05_vorticity_dynamics as ch05  # noqa: E402
from fluidpy.core.style import COLORS  # noqa: E402


def column_figure():
    import matplotlib.pyplot as plt
    from ch05_drawings import column_sketch

    x = np.linspace(-3e5, 3e5, 241)
    fig, ax = plt.subplots(2, 3, figsize=(15, 7), sharex=True)
    for j, dep in enumerate(("ridge", "trough", "slope")):
        cs = ch05.column_over_slope(x, dep, 45.0)
        column_sketch(ax[0, j], x, cs["h"], cs["zeta"])
        ax[0, j].set_title(f"{dep}: f = {cs['f']:.3e} 1/s", fontsize=10)
        ax[1, j].plot(x / 1e3, cs["zeta"] / cs["f"], color=COLORS["teal"], label="ζ/f")
        ax[1, j].plot(x / 1e3, (cs["h"] / 1000.0 - 1.0), "--", color=COLORS["muted"], label="h/h₀ − 1")
        ax[1, j].set_xlabel("x [km]")
        ax[1, j].legend(fontsize=8)
    ax[1, 0].set_ylabel("ζ/f")
    return fig


def gamma_a_figure():
    import matplotlib.pyplot as plt

    s = ch05.kelvin_scenario("rotating")
    te = np.linspace(0, s["t_end"], 41)
    G, loops = ch05.material_circulation(s["u"], s["pts0"], te, return_loops=True)
    Ga = [ch05.absolute_circulation(s["u"], L, s["Omega"], tt)[1] for L, tt in zip(loops, te)]
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(te, G, color=COLORS["teal"], label="relative Γ (material loop)")
    ax.plot(te, Ga, color=COLORS["amber"], label="absolute Γ_a = Γ + 2Ω·A (5.33)")
    ax.plot(te, s["Gamma_exact"](te), "--", color=COLORS["muted"], label="2ΩA₀(1 − e^{−αt})")
    ax.set_xlabel("t [s]")
    ax.set_ylabel("[m²/s]")
    ax.legend(fontsize=8)
    ax.set_title("a converging ring in a rotating frame spins up")
    return fig, np.ptp(Ga)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch05"))
    ap.add_argument("--no-show", action="store_true")
    args = ap.parse_args()
    import matplotlib

    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from fluidpy.core.rotating import OMEGA_EARTH, coriolis_parameter
    from fluidpy.core.style import use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    for lat in (0, 30, 45, 60, 90):
        print(f"f({lat:2d}°) = 2Ω sin φ = {coriolis_parameter(np.deg2rad(lat)):.6e} 1/s")
    print(f"column 1000 → 1100 m, f = 1e-4: ζ = {ch05.column_relative_vorticity(1100.0, 1000.0, 0.0, 1e-4):.3e} 1/s "
          f"(cyclonic); 1000 → 900 m: {ch05.column_relative_vorticity(900.0, 1000.0, 0.0, 1e-4):.3e} (anticyclonic)")
    for dep in ("ridge", "trough", "slope"):
        cs = ch05.column_over_slope(np.linspace(-3e5, 3e5, 241), dep, 45.0)
        print(f"   {dep:6s}: ζ/f from {cs['zeta'].min() / cs['f']:+.3f} to {cs['zeta'].max() / cs['f']:+.3f}; "
              f"(ζ + f)/h spread {np.ptp(cs['ratio']):.1e}")
    A = np.pi * 5e5 ** 2
    G1 = ch05.relative_circulation_after_move(0.0, A, 30.0, A, 60.0)
    print(f"ring of air r = 500 km at rest, 30° → 60° N: Γ₁ = {G1:.5e} m²/s, mean ζ = {G1 / A:.5e} 1/s "
          f"(anticyclonic)")
    G = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1e-5]])
    print(f"planetary stretching 2(Ω·∇)u with ∂w/∂z = 1e-5 1/s: "
          f"{ch05.planetary_vorticity_terms(G, [0, 0, OMEGA_EARTH])[2]:.5e} 1/s²")
    pts = ch05.circle_loop_points((0.3, -0.2), 0.7, 128)
    N = pts.shape[1]
    x_next = np.roll(pts, -1, axis=1)
    area_by_hand = 0.5 * np.sum(pts[0] * x_next[1] - pts[1] * x_next[0])  # polygon (shoelace) ½∮x × dx
    print(f"vector area by hand (shoelace) {area_by_hand:.6f} vs spectral {ch05.loop_vector_area(pts)[2]:.6f} vs "
          f"πr² = {np.pi * 0.49:.6f} m²")
    for t in (0.0, 5.0, 10.0):
        print(f"rotating scenario t = {t:4.1f} s: Γ = {ch05.kelvin_scenario_circulation('rotating', t):.6f} m²/s "
              f"(Γ_a = {2 * 0.5 * np.pi:.6f})")
    column_figure().savefig(out / "c11_column.png", bbox_inches="tight")
    fig, spread = gamma_a_figure()
    fig.savefig(out / "c11_absolute_circulation.png", bbox_inches="tight")
    print(f"material loop in the rotating frame: Γ_a constant to {spread:.1e} m²/s")
    print(f"figures → {out}")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
