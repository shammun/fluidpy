"""Chapter 1, §1.7: gauge vs absolute pressure, Pascal's law (1.7), the hydrostatic law (1.8) integrated from scratch
and with ``integrate_hydrostatic``, uniform and layered tanks (1.9), and buoyancy as a net pressure force (D37).

Run: ``.venv/Scripts/python.exe scripts/ch01_hydrostatics.py --no-show``
Figures → outputs/ch01/fig_hydrostatic_tank.png, fig_pascal_slice.png.
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
    g = ch01.G0

    # --- layered tank: 2 m of oil (880 kg/m^3, illustrative) over 3 m of water ------------------------------------------
    z = np.linspace(0.0, -5.0, 501)
    thick, dens = [2.0, 3.0], [880.0, 998.2]
    p_layers = ch01.layered_pressure(z, thick, dens)

    # smooth density rho(z) (a salty bottom layer) — from-scratch Euler vs the library integrator
    def rho_fn(zz, _p=None):
        return 998.2 + 25.0 / (1.0 + np.exp((zz + 3.5) / 0.2))

    p_lib = ch01.integrate_hydrostatic(z, rho_fn, ch01.P_ATM)
    p_euler = np.empty_like(z)
    p_euler[0] = ch01.P_ATM
    for k in range(z.size - 1):  # p[k+1] = p[k] − rho(z_k) g Δz
        p_euler[k + 1] = p_euler[k] - rho_fn(z[k]) * g * (z[k + 1] - z[k])
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)
    axs[0].plot(ch01.gauge_pressure(p_layers) / 1e3, z, color=COLORS["orange"], label="oil over water (1.9 layer by layer)")
    axs[0].plot(ch01.gauge_pressure(ch01.hydrostatic_pressure_uniform(z, ch01.P_ATM, 998.2)) / 1e3, z,
                color=COLORS["blue"], ls="--", label="water only")
    axs[0].axhline(-2.0, color=COLORS["muted"], lw=0.8, ls=":")
    axs[0].set_xlabel("gauge pressure [kPa]")
    axs[0].set_ylabel("z [m] (free surface at 0)")
    axs[0].set_title("Slope dp/dz = −ρg changes at the interface")
    axs[0].legend(fontsize=8)
    axs[1].plot(ch01.gauge_pressure(p_lib) / 1e3, z, color=COLORS["accent"], label="integrate_hydrostatic (1.8)")
    axs[1].plot(ch01.gauge_pressure(p_euler)[::25] / 1e3, z[::25], "o", color=COLORS["teal"], ms=4,
                label="from-scratch Euler march")
    axs[1].set_xlabel("gauge pressure [kPa]")
    axs[1].set_title("Density varying with depth")
    axs[1].legend(fontsize=8)
    fig.savefig(out / "fig_hydrostatic_tank.png", bbox_inches="tight")

    # --- Pascal: a 3-D hydrostatic field has no horizontal gradient -------------------------------------------------------
    x = np.linspace(0.0, 1.0, 21)
    y = np.linspace(0.0, 1.0, 21)
    zz = np.linspace(-1.0, 0.0, 21)
    X, Y, Z = np.meshgrid(x, y, zz, indexing="ij")
    P = ch01.hydrostatic_pressure_uniform(Z, ch01.P_ATM, 998.2)
    dpdx, dpdy, dpdz = np.gradient(P, x, y, zz, edge_order=2)
    fig, ax = plt.subplots()
    cs = ax.contourf(x, zz, ch01.gauge_pressure(P[:, 10, :]).T, levels=15, cmap="viridis")
    fig.colorbar(cs, ax=ax, label="gauge pressure [Pa]")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("z [m]")
    ax.set_title("Isobars are horizontal in a fluid at rest (Eq. 1.7)")
    fig.savefig(out / "fig_pascal_slice.png", bbox_inches="tight")

    # --- buoyancy on a submerged box ----------------------------------------------------------------------------------------
    box = (0.0, 0.2, 0.0, 0.3, -1.5, -1.0)
    F = ch01.net_pressure_force_on_box(lambda xx, yy, zq: ch01.hydrostatic_pressure_uniform(zq, ch01.P_ATM, 998.2), box)
    V = 0.2 * 0.3 * 0.5
    print(f"tyre gauge 2 bar -> absolute {ch01.absolute_pressure(2e5) / 1e3:.1f} kPa; 10 m of water adds "
          f"{998.2 * g * 10 / 1e3:.1f} kPa ({998.2 * g * 10 / ch01.P_ATM:.3f} atm)")
    print(f"layered tank bottom (5 m): {ch01.gauge_pressure(p_layers[-1]) / 1e3:.3f} kPa gauge; "
          f"Euler vs solve_ivp max diff {np.max(np.abs(p_euler - p_lib)):.2f} Pa (dz = 1 cm)")
    print(f"Pascal: max |dp/dx| = {np.abs(dpdx).max():.2e}, max |dp/dy| = {np.abs(dpdy).max():.2e}, "
          f"dp/dz = {dpdz.mean():.2f} Pa/m (−rho g = {-998.2 * g:.2f})")
    print(f"box V = {V:.3f} m^3: net pressure force {F} N; rho g V = {ch01.buoyancy_force(998.2, V):.3f} N")
    print(f"saved 2 figures in {out}  ({time.perf_counter() - t_start:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
