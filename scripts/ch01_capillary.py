"""Chapter 1, §1.6–1.7: surface tension of water, the Laplace pressure jump (1.5), capillary rise (Example 1.1)
and the pressure along E–F in the tube (Fig. 1.7 idea).

Run: ``.venv/Scripts/python.exe scripts/ch01_capillary.py --no-show``
Figures → outputs/ch01/fig_surface_tension.png, fig_laplace_jump.png, fig_capillary_rise.png.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch01"), help="figure folder")
    ap.add_argument("--no-show", action="store_true", help="do not open a window")
    args = ap.parse_args()
    t_start = time.perf_counter()
    import matplotlib

    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from fluidpy import ch01_introduction as ch01
    from fluidpy.core.style import COLORS, use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    # --- σ(T) -------------------------------------------------------------------------------------------------------------
    T = np.linspace(273.16, 373.15, 101)
    sigma = ch01.surface_tension_water(T)
    T20 = ch01.celsius_to_kelvin(20.0)
    s20 = ch01.surface_tension_water(T20)
    fig, ax = plt.subplots()
    ax.plot(T - 273.15, sigma * 1e3, color=COLORS["blue"])
    ax.plot([20.0], [s20 * 1e3], "o", color=COLORS["orange"])
    ax.annotate(f"{s20 * 1e3:.2f} mN/m at 20 °C", (20.0, s20 * 1e3), xytext=(8, 8), textcoords="offset points")
    ax.set_xlabel("T [°C]")
    ax.set_ylabel("σ water–air [mN/m]")
    ax.set_title("Surface tension of water (IAPWS R1-76)")
    fig.savefig(out / "fig_surface_tension.png", bbox_inches="tight")

    # --- Laplace jump vs radius ----------------------------------------------------------------------------------------------
    R = np.logspace(-7, -2, 100)
    fig, ax = plt.subplots()
    ax.loglog(R, ch01.laplace_pressure_jump(s20, R), color=COLORS["accent"], label="sphere: 2σ/R")
    ax.loglog(R, ch01.laplace_pressure_jump(s20, R, np.inf), color=COLORS["teal"], ls="--", label="cylinder: σ/R")
    ax.axhline(ch01.P_ATM, color=COLORS["muted"], lw=0.8, ls=":")
    ax.text(2e-7, ch01.P_ATM * 1.3, "1 atm", color=COLORS["muted"])
    ax.set_xlabel("radius of curvature R [m]")
    ax.set_ylabel(r"$p_i - p_o$ [Pa]")
    ax.set_title("Small drops and bubbles carry large pressure jumps (Eq. 1.5)")
    ax.legend()
    fig.savefig(out / "fig_laplace_jump.png", bbox_inches="tight")

    # --- capillary rise and pressure along E–F ----------------------------------------------------------------------------
    rho = float(ch01.water_density(T20))
    Rt = np.logspace(-4.5, -2, 100)
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    for a_deg, c in ((90, COLORS["accent"]), (60, COLORS["teal"]), (30, COLORS["orange"])):
        axs[0].loglog(Rt * 1e3, ch01.capillary_rise_deg(s20, a_deg, rho, Rt) * 1e3, color=c,
                      label=f"book's α = {a_deg}° (contact angle {90 - a_deg}°)")
    axs[0].set_xlabel("tube radius R [mm]")
    axs[0].set_ylabel("rise h [mm]")
    axs[0].set_title("h = 2σ sin α/(ρ g R): halve R, double h")
    axs[0].legend(fontsize=8)
    R1, alpha = 0.5e-3, np.radians(90.0)
    h1 = ch01.capillary_rise(s20, alpha, rho, R1)
    z = np.linspace(0.0, h1, 50)  # from F (z = 0, outside level) to E (just under the meniscus)
    p = ch01.hydrostatic_pressure_uniform(z, ch01.P_ATM, rho)  # Eq. (1.9), p(F) = p_atm
    axs[1].plot(ch01.gauge_pressure(p), z * 1e3, color=COLORS["orange"])
    axs[1].axvline(0.0, color=COLORS["muted"], lw=0.8)
    axs[1].set_xlabel("gauge pressure along F→E [Pa]")
    axs[1].set_ylabel("height above outside level [mm]")
    axs[1].set_title(f"R = {R1 * 1e3:g} mm: p(E) is below atmospheric")
    fig.savefig(out / "fig_capillary_rise.png", bbox_inches="tight")

    hmm = ch01.capillary_rise(s20, np.pi / 2, rho, 1e-3)
    pE_laplace = ch01.P_ATM - ch01.laplace_pressure_jump(s20, R1 / np.sin(alpha))
    print(f"sigma(20 °C) = {s20 * 1e3:.3f} mN/m, sigma(25 °C) = {ch01.surface_tension_water(298.15) * 1e3:.3f} mN/m")
    print(f"Laplace jump of a 1 mm radius drop: {ch01.laplace_pressure_jump(s20, 1e-3):.1f} Pa; 1 µm: "
          f"{ch01.laplace_pressure_jump(s20, 1e-6):.3e} Pa")
    print(f"capillary rise of water (20 °C, full wetting) in a 1 mm radius tube: {hmm * 1e3:.2f} mm")
    print(f"R = 0.5 mm: p_E from hydrostatics {ch01.gauge_pressure(p[-1]):.2f} Pa gauge, from Laplace "
          f"{ch01.gauge_pressure(pE_laplace):.2f} Pa gauge (consistency of Ex. 1.1 with Eq. 1.5)")
    print(f"saved 3 figures in {out}  ({time.perf_counter() - t_start:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
