"""Chapter 1, §1.4: pressure from molecular impacts, the continuum averaging volume, and the Knudsen number.

Run: ``.venv/Scripts/python.exe scripts/ch01_continuum.py --no-show [--fast]``
Figures → outputs/ch01/fig_kinetic_pressure.png, fig_continuum_noise.png, fig_knudsen.png.
Random draws are seeded (numpy default_rng(0)).
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
    ap.add_argument("--fast", action="store_true", help="fewer molecules and samples")
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
    rng = np.random.default_rng(0)

    T, p = 300.0, ch01.P_ATM
    m = ch01.molecule_mass(ch01.M_W_AIR)
    n = p / (ch01.K_B * T)  # Eq. (1.21) per unit volume

    # --- kinetic pressure converging with molecule count ---------------------------------------------------------------
    Ns = np.unique(np.logspace(1, 5 if args.fast else 6, 11).astype(int))
    est = np.array([ch01.molecular_pressure(n, m, ch01.maxwellian_velocities(N, T, m, rng)) for N in Ns])
    fig, ax = plt.subplots()
    ax.loglog(Ns, np.abs(est / p - 1), "o-", color=COLORS["accent"], label="|p_kinetic/(n k_B T) − 1|")
    ax.loglog(Ns, np.sqrt(2.0 / (3.0 * Ns)), "--", color=COLORS["muted"], label=r"expected noise $\sqrt{2/(3N)}$")
    ax.set_xlabel("number of sampled molecules N")
    ax.set_ylabel("relative error [-]")
    ax.set_title("Pressure as an average of molecular impacts")
    ax.legend()
    fig.savefig(out / "fig_kinetic_pressure.png", bbox_inches="tight")

    # --- density in a box: noise -> plateau -> macroscopic drift ----------------------------------------------------------
    L = np.logspace(-9, -2, 36 if args.fast else 71)  # box side [m]
    grad, L_flow = 0.5, 1e-3  # 50 % density change over 1 mm (a strong, illustrative gradient; E1's linear model)
    mean, std = ch01.sample_density(L, n, m, rng=rng, n_samples=100 if args.fast else 200, gradient=grad,
                                    L_flow=L_flow)
    rho0 = n * m
    rho_point = rho0  # the linear model's box spans [0, L] from the point, where the density is rho0
    fig, ax = plt.subplots()
    noise = ch01.density_noise_expected(L, n)
    ax.fill_between(L, 1 - noise, 1 + noise, color=COLORS["grid"], label=r"expected $\pm (nL^3)^{-1/2}$")
    ax.semilogx(L, mean / rho_point, ".", color=COLORS["accent"], label="sample mean")
    ax.semilogx(L, ch01.density_expected(L, n, m, grad, L_flow) / rho_point, color=COLORS["orange"],
                label="macroscopic variation seen by the box")
    ax.set_ylim(0.0, 2.0)
    ax.set_xlabel("box side L [m]")
    ax.set_ylabel(r"measured $\delta m/\delta V$ divided by $\rho$ at the point")
    ax.set_title("The continuum value is a plateau between molecular noise and flow variation")
    ax.legend(fontsize=8)
    fig.savefig(out / "fig_continuum_noise.png", bbox_inches="tight")

    # --- Knudsen number vs body size at several altitudes ----------------------------------------------------------------
    body = np.logspace(-9, 1, 200)
    fig, ax = plt.subplots()
    for zkm, c in zip([0, 30, 60, 84], [COLORS["accent"], COLORS["teal"], COLORS["orange"], COLORS["rose"]]):
        Tz, pz, _ = ch01.standard_atmosphere(zkm * 1000.0)
        lz = ch01.mean_free_path_air(Tz, pz)
        ax.loglog(body, ch01.knudsen_number(lz, body), color=c, label=f"z = {zkm} km (l = {lz:.2g} m)")
    for y0, y1, lab, col in ((1e-12, 0.01, "continuum", COLORS["teal"]), (0.01, 0.1, "slip", COLORS["amber"]),
                             (0.1, 10.0, "transition", COLORS["orange"]), (10.0, 1e12, "free molecular", COLORS["rose"])):
        ax.axhspan(y0, y1, color=col, alpha=0.06)
        ax.text(2e-9, np.sqrt(max(y0, 1e-6) * min(y1, 1e6)), lab, fontsize=8, color=col)
    ax.set_ylim(1e-8, 1e6)
    ax.set_xlabel("body size L [m]")
    ax.set_ylabel("Kn = l/L [-]")
    ax.set_title("Where the continuum breaks down (USSA-1976 air)")
    ax.legend(fontsize=8)
    fig.savefig(out / "fig_knudsen.png", bbox_inches="tight")

    l300 = ch01.mean_free_path_air(T, p)
    print(f"air at {T:.0f} K, 1 atm: n = {n:.3e} 1/m^3, m = {m:.3e} kg, mean free path = {l300 * 1e9:.1f} nm, "
          f"collision time = {ch01.collision_time(l300, T, m):.2e} s, mean speed = {ch01.mean_molecular_speed(T, m):.0f} m/s")
    print(f"kinetic pressure with N = {Ns[-1]}: p/(n k_B T) = {est[-1] / p:.5f}; wall-impact form: "
          f"{ch01.wall_impact_pressure(ch01.maxwellian_velocities(Ns[-1], T, m, rng), m, n) / p:.5f}")
    N10um = n * (10e-6) ** 3
    print(f"molecules in a 10 µm cube: {N10um:.2e}, expected relative noise {N10um ** -0.5:.1e}")
    i = np.searchsorted(L, 1e-6)
    print(f"box L = {L[i]:.2e} m: relative std measured {std[i] / mean[i]:.2e}, expected {noise[i]:.2e}")
    j = np.searchsorted(L, 1e-5)
    print(f"box L = {L[j]:.2e} m: mean/rho(point) = {mean[j] / rho_point:.6f} (plateau); box L = {L[-1]:.0e} m: "
          f"{mean[-1] / rho_point:.4f} (drift: density_expected gives "
          f"{ch01.density_expected(L[-1], n, m, grad, L_flow) / rho_point:.4f})")
    print(f"Kn of a 1 µm droplet at sea level: {ch01.knudsen_number(ch01.mean_free_path_air(288.15), 1e-6):.3f}")
    print(f"saved 3 figures in {out}  ({time.perf_counter() - t_start:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
