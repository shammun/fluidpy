"""Chapter 1, §1.5: molecular transport — Fick's flux (Fig. 1.2 idea), momentum diffusion between plates (Fig. 1.3
idea: FTCS start-up vs the series solution, wall shear stress), viscosity vs temperature, μ vs ν.

Run: ``.venv/Scripts/python.exe scripts/ch01_transport_relaxation.py --no-show [--fast]``
Figures → outputs/ch01/fig1_2_fick_flux.png, fig1_3_couette_startup.png, fig_viscosity_T.png, fig_mu_vs_nu.png.
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
    ap.add_argument("--fast", action="store_true", help="coarser grid")
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

    # --- Fick: Y(y) increasing upward, flux points down --------------------------------------------------------------------
    y = np.linspace(0.0, 1e-2, 41)  # [m]
    Y = 0.2 + 0.5 * (y / y[-1]) ** 2  # mass fraction (illustrative)
    rho, kappa_m = 1.2, ch01.FLUIDS["air"]["kappa_m"]
    J = ch01.fick_mass_flux(rho, kappa_m, ch01.derivative_2nd_order(Y, y))  # Eq. (1.1), y component
    fig, axs = plt.subplots(1, 2, figsize=(8, 4), sharey=True)
    axs[0].plot(Y, y * 1e3, color=COLORS["accent"])
    axs[0].set_xlabel("mass fraction Y [-]")
    axs[0].set_ylabel("y [mm]")
    axs[1].plot(J, y * 1e3, color=COLORS["orange"])
    for yi, Ji in zip(y[5::10], J[5::10]):
        axs[1].annotate("", xy=(Ji, yi * 1e3 - 1.2), xytext=(Ji, yi * 1e3), arrowprops=dict(arrowstyle="-|>", color=COLORS["orange"]))
    axs[1].set_xlabel(r"flux $J_{m,y}$ [kg m$^{-2}$ s$^{-1}$]")
    axs[0].set_title("Y rises with y …")
    axs[1].set_title("… so the diffusive flux points down (Eq. 1.1)")
    fig.savefig(out / "fig1_2_fick_flux.png", bbox_inches="tight")

    # --- Couette start-up: FTCS (D = ν) vs series ----------------------------------------------------------------------------
    h, U = 1e-2, 1e-2  # gap [m], plate speed [m/s]
    water = ch01.FLUIDS["water"]
    nu, mu = water["nu"], water["mu"]
    N = 51 if args.fast else 101
    yy = np.linspace(0.0, h, N)
    dy = yy[1] - yy[0]
    t_visc = ch01.diffusion_time(h, nu)  # h^2/ν
    dt = ch01.stable_time_step(nu, dy, safety=0.8)  # r = 0.4
    t_end = 1.0 * t_visc
    nsteps = int(np.ceil(t_end / dt))
    dt = t_end / nsteps
    u0 = np.zeros(N)
    u0[-1] = U  # plate starts moving at t = 0
    F = ch01.ftcs_diffusion_1d(u0, nu, dy, dt, nsteps, bc=("dirichlet", "dirichlet"), values=(0.0, U))
    times = np.arange(nsteps + 1) * dt
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    errs = []
    for frac, c in zip([0.02, 0.1, 0.3, 1.0], [COLORS["rose"], COLORS["orange"], COLORS["teal"], COLORS["accent"]]):
        k = int(round(frac * t_visc / dt))
        exact = ch01.couette_startup_profile(yy, times[k], U, h, nu)
        errs.append(float(np.max(np.abs(F[k] - exact)) / U))
        axs[0].plot(exact / U, yy / h, color=c, lw=1, ls="--")
        axs[0].plot(F[k] / U, yy / h, color=c, lw=2, label=f"t = {frac:g} h²/ν")
    axs[0].set_xlabel("u/U [-]")
    axs[0].set_ylabel("y/h [-]")
    axs[0].set_title("Momentum diffuses in from the moving plate (solid FTCS, dashed series)")
    axs[0].legend(fontsize=8)
    tb, tt = ch01.wall_shear_history(F, dy, mu)
    axs[1].plot(times / t_visc, tb / (mu * U / h), color=COLORS["teal"], label="bottom (fixed) wall")
    axs[1].plot(times / t_visc, tt / (mu * U / h), color=COLORS["accent"], label="top (moving) wall")
    axs[1].set_ylim(0, 5)
    axs[1].set_xlabel("t / (h²/ν) [-]")
    axs[1].set_ylabel(r"wall stress $\mu\,\partial u/\partial y$ / ($\mu U/h$)")
    axs[1].set_title("Both wall stresses settle to μU/h (Eq. 1.3)")
    axs[1].legend()
    fig.savefig(out / "fig1_3_couette_startup.png", bbox_inches="tight")

    # --- viscosity vs temperature -----------------------------------------------------------------------------------------
    T = np.linspace(200.0, 400.0, 201)
    mu_s = ch01.sutherland_viscosity(T)
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    axs[0].plot(T, mu_s * 1e6, color=COLORS["accent"], label="Sutherland (USSA-1976)")
    for n_exp, c in ((0.5, COLORS["orange"]), (0.7, COLORS["teal"])):
        axs[0].plot(T, ch01.viscosity_power_law(T, float(ch01.sutherland_viscosity(288.15)), 288.15, n_exp) * 1e6,
                    ls="--", color=c, label=f"μ ∝ T^{n_exp}")
    axs[0].set_xlabel("T [K]")
    axs[0].set_ylabel("μ of air [µPa s]")
    axs[0].set_title("Gas: μ rises with T")
    axs[0].legend()
    Tw = np.linspace(273.15, 373.15, 101)
    axs[1].plot(Tw - 273.15, ch01.water_viscosity(Tw) * 1e3, color=COLORS["blue"])
    axs[1].set_xlabel("T [°C]")
    axs[1].set_ylabel("μ of water [mPa s]")
    axs[1].set_title("Liquid: μ falls with T")
    fig.savefig(out / "fig_viscosity_T.png", bbox_inches="tight")

    # --- μ vs ν bars ------------------------------------------------------------------------------------------------------
    names = list(ch01.FLUIDS)
    fig, ax = plt.subplots(figsize=(8, 4))
    x = np.arange(len(names))
    ax.bar(x - 0.2, [ch01.FLUIDS[k]["mu"] for k in names], 0.4, color=COLORS["rose"], label="μ [Pa s]")
    ax.bar(x + 0.2, [ch01.FLUIDS[k]["nu"] for k in names], 0.4, color=COLORS["teal"], label="ν [m²/s]")
    ax.set_yscale("log")
    ax.set_xticks(x, names)
    ax.set_title("Dynamic vs kinematic viscosity at 20 °C")
    ax.legend()
    fig.savefig(out / "fig_mu_vs_nu.png", bbox_inches="tight")

    air = ch01.FLUIDS["air"]
    print(f"Fick flux at mid-height: {J[20]:.3e} kg/(m^2 s) (negative = downward)")
    print(f"water: nu = {nu:.4e} m^2/s, h^2/nu = {t_visc:.1f} s for h = 1 cm; FTCS N = {N}, steps = {nsteps}, "
          f"r = {nu * dt / dy**2:.3f}")
    print("max |FTCS − series|/U at t/(h^2/nu) = 0.02, 0.1, 0.3, 1: " + ", ".join(f"{e:.2e}" for e in errs))
    print(f"wall stress at t = h^2/nu: bottom {tb[-1] / (mu * U / h):.4f}, top {tt[-1] / (mu * U / h):.4f} (× mu U/h)")
    print(f"mu_air/mu_water = {air['mu'] / water['mu']:.4f}; nu_air/nu_water = {air['nu'] / water['nu']:.2f}")
    print(f"air viscosity at 216.65 K: Sutherland {ch01.sutherland_viscosity(216.65):.4e}, T^0.5 law "
          f"{ch01.viscosity_power_law(216.65, float(ch01.sutherland_viscosity(288.15)), 288.15):.4e} Pa s")
    print(f"saved 4 figures in {out}  ({time.perf_counter() - t_start:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
