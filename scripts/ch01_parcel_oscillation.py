"""Chapter 1, §1.10, Fig. 1.8 idea: a displaced parcel oscillates (N² > 0), stays (N² = 0) or runs away (N² < 0).
Linear solution ζ'' + N²ζ = 0 vs the nonlinear buoyancy equation, in the atmosphere (lapse rates, book sign) and in
an ocean thermocline.

Run: ``.venv/Scripts/python.exe scripts/ch01_parcel_oscillation.py --no-show``
Figures → outputs/ch01/fig1_8_parcel_atmosphere.png, fig1_8_parcel_ocean.png.
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

    T0, zeta0 = 288.15, 100.0
    Ga = ch01.adiabatic_lapse_rate()
    cases = [("standard troposphere, dT/dz = −6.5 K/km", -6.5e-3, COLORS["teal"]),
             (f"dry adiabatic, dT/dz = Γ_a = {Ga * 1e3:.2f} K/km", Ga, COLORS["amber"]),
             ("superadiabatic, dT/dz = −12 K/km", -12e-3, COLORS["rose"])]
    t = np.linspace(0.0, 1500.0, 301)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    report = []
    for label, dTdz, c in cases:
        N2 = ch01.brunt_vaisala_sq_from_lapse(T0, dTdz)
        lin = ch01.parcel_displacement(t, zeta0, N2)
        nl = ch01.parcel_ode_atmosphere(t, zeta0, T0, dTdz, zeta_max=5e3)
        ax.plot(t / 60, lin, color=c, ls="--", lw=1)
        ax.plot(t / 60, nl, color=c, lw=2, label=f"{label} ({ch01.classify_stability(N2, tol=1e-9)})")
        report.append((label, N2, ch01.stability_timescale(N2, tol=1e-9), lin, nl))
    ax.set_ylim(-400, 1500)
    ax.set_xlabel("time after release [min]")
    ax.set_ylabel("displacement ζ [m]")
    ax.set_title("Parcel released 100 m above its rest height (solid nonlinear, dashed linear)")
    ax.legend(fontsize=8)
    fig.savefig(out / "fig1_8_parcel_atmosphere.png", bbox_inches="tight")

    # --- ocean thermocline: linear densities, compressible parcel --------------------------------------------------------
    rho0, c = 1026.0, 1500.0
    drho_a = ch01.isentropic_density_gradient(rho0, c)
    excess = 0.01  # kg/m^4: a strong thermocline (≈ 1 kg/m^3 over 100 m) gives N ≈ 1e-2 1/s, period ≈ 10 min
    drho = drho_a - excess  # environment denser downward than the isentropic gradient by `excess` (synthetic)
    N2o = ch01.brunt_vaisala_sq(rho0, drho, drho_a)
    to = np.linspace(0.0, 3.0 * 2 * np.pi / np.sqrt(N2o), 301)
    fig, ax = plt.subplots(figsize=(8, 3.8))
    ax.plot(to / 60, ch01.parcel_displacement(to, 2.0, N2o), color=COLORS["blue"], ls="--", label="linear")
    ax.plot(to / 60, ch01.parcel_ode_from_gradients(to, 2.0, rho0, drho, drho_a), color=COLORS["accent"],
            label="nonlinear")
    ax.set_xlabel("time [min]")
    ax.set_ylabel("ζ [m]")
    ax.set_title(f"Ocean thermocline: N = {np.sqrt(N2o):.3e} s⁻¹, period {2 * np.pi / np.sqrt(N2o) / 60:.1f} min")
    ax.legend()
    fig.savefig(out / "fig1_8_parcel_ocean.png", bbox_inches="tight")

    # from-scratch Euler–Cromer check for the stable atmosphere case
    N2 = report[0][1]
    dt, zz, ww = 0.5, zeta0, 0.0
    for _ in range(int(600 / dt)):
        ww -= N2 * zz * dt
        zz += ww * dt
    for label, N2, ts, lin, nl in report:
        print(f"{label}: N^2 = {N2:.3e} 1/s^2, {ts[0]} = {ts[1]:.1f} s; zeta(25 min): linear {lin[-1]:.1f} m, "
              f"nonlinear {nl[-1]:.1f} m")
    print(f"Euler–Cromer (dt = 0.5 s) zeta(600 s) = {zz:.3f} m vs analytic {ch01.parcel_displacement(600.0, zeta0, report[0][1]):.3f} m")
    print(f"ocean (strong thermocline, excess gradient {excess} kg/m^4): drho_a/dz = {drho_a:.3e} kg/m^4, "
          f"N^2 = {N2o:.3e} 1/s^2, N = {np.sqrt(N2o):.3e} 1/s, period = {2 * np.pi / np.sqrt(N2o) / 60:.1f} min")
    print(f"saved 2 figures in {out}  ({time.perf_counter() - t_start:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
