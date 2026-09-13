"""Chapter 1, §1.8–1.9: heat and work depend on the path, internal energy and entropy do not; stirring and free
expansion (Clausius–Duhem with the actual heat); isentropes vs isotherms; e(v) of a perfect gas vs van der Waals.

Run: ``.venv/Scripts/python.exe scripts/ch01_process_paths.py --no-show``
Figures → outputs/ch01/fig_pv_paths.png, fig_heat_work_bars.png, fig_isentrope.png, fig_vdw_energy.png.
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
    R, cv = ch01.R_AIR, ch01.CV_AIR

    s1, s2 = (0.8, 300.0), (1.6, 450.0)  # (v [m^3/kg], T [K]) of the two end states (illustrative)
    kinds = ["isochoric-isobaric", "isobaric-isochoric", "isentropic-isochoric", "isothermal-isochoric"]
    colors = [COLORS["accent"], COLORS["teal"], COLORS["orange"], COLORS["blue"]]

    # --- p–v diagram with work areas ----------------------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 4.5))
    vv = np.linspace(0.6, 1.8, 200)
    for Tiso in (250, 300, 350, 400, 450, 500):
        ax.plot(vv, R * Tiso / vv / 1e3, color=COLORS["grid"], lw=1)
    rows = []
    for k, c in zip(kinds, colors):
        P = ch01.process_path(k, s1, s2, n=2001)
        num = ch01.process_heat_work(P["v"], P["T"])
        ex = ch01.path_heat_work_totals(k, *s1, *s2)
        ds = ch01.entropy_change_reversible(num["q"], P["T"])
        rows.append((k, num, ex, ds))
        ax.plot(P["v"], P["p"] / 1e3, color=c, lw=2, label=k)
    ax.plot([s1[0], s2[0]], [R * s1[1] / s1[0] / 1e3, R * s2[1] / s2[0] / 1e3], "ko")
    ax.set_xlabel("specific volume v [m³/kg]")
    ax.set_ylabel("pressure p [kPa]")
    ax.set_title("Four reversible paths between the same two states (grey: isotherms)")
    ax.legend(fontsize=8)
    fig.savefig(out / "fig_pv_paths.png", bbox_inches="tight")

    # --- bars: q, w, Δe, Δs ---------------------------------------------------------------------------------------------------
    stir = ch01.irreversible_process("stirring", 300.0, 0.8, T2=450.0)
    fexp = ch01.irreversible_process("free_expansion", 300.0, 0.8, v2=1.6)
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    x = np.arange(len(kinds))
    wdt = 0.25
    axs[0].bar(x - wdt, [r[2]["q"] / 1e3 for r in rows], wdt, color=COLORS["orange"], label="q (heat in)")
    axs[0].bar(x, [r[2]["w"] / 1e3 for r in rows], wdt, color=COLORS["blue"], label="w (work on gas)")
    axs[0].bar(x + wdt, [r[2]["de"] / 1e3 for r in rows], wdt, color=COLORS["accent"], label="Δe")
    axs[0].set_xticks(x, [k.replace("-", "\n") for k in kinds], fontsize=8)
    axs[0].set_ylabel("[kJ/kg]")
    axs[0].set_title("q and w change with the path, Δe does not (Eq. 1.10)")
    axs[0].legend(fontsize=8)
    labels = ["reversible\n(any path)", "stirring\n(v fixed)", "free expansion\n(T fixed)"]
    ds_vals = [rows[0][2]["ds"], stir["ds"], fexp["ds"]]
    int_vals = [rows[0][2]["int_dq_over_T"], stir["int_dq_over_T"], fexp["int_dq_over_T"]]
    xx = np.arange(3)
    axs[1].bar(xx - 0.18, ds_vals, 0.36, color=COLORS["accent"], label="Δs (from Gibbs 1.18)")
    axs[1].bar(xx + 0.18, int_vals, 0.36, color=COLORS["rose"], label="∫ δq/T (actual heat)")
    axs[1].set_xticks(xx, labels, fontsize=8)
    axs[1].set_ylabel("[J/(kg K)]")
    axs[1].set_title("Clausius–Duhem: Δs ≥ ∫δq/T, equality only if reversible")
    axs[1].legend(fontsize=8)
    fig.savefig(out / "fig_heat_work_bars.png", bbox_inches="tight")

    # --- isentrope vs isotherm (log–log slopes γ and 1), from-scratch ds = 0 integration --------------------------------------
    rho = np.logspace(-0.5, 0.5, 60)
    p_s = ch01.isentropic_pressure(rho, ch01.P_ATM, 1.0)
    p_t = ch01.P_ATM * rho / 1.0
    p_scratch = np.empty_like(rho)
    p_scratch[0] = p_s[0]
    for k in range(rho.size - 1):  # dp/drho = gamma p/rho, midpoint step
        drho = rho[k + 1] - rho[k]
        pm = p_scratch[k] + 0.5 * drho * ch01.GAMMA_AIR * p_scratch[k] / rho[k]
        p_scratch[k + 1] = p_scratch[k] + drho * ch01.GAMMA_AIR * pm / (rho[k] + 0.5 * drho)
    fig, ax = plt.subplots()
    ax.loglog(rho, p_s / 1e3, color=COLORS["orange"], label=r"isentrope $p\propto\rho^{\gamma}$ (1.25)")
    ax.loglog(rho[::5], p_scratch[::5] / 1e3, "o", color=COLORS["teal"], ms=4, label="from-scratch ds = 0 steps")
    ax.loglog(rho, p_t / 1e3, color=COLORS["blue"], ls="--", label=r"isotherm $p\propto\rho$")
    ax.set_xlabel("ρ [kg/m³]")
    ax.set_ylabel("p [kPa]")
    ax.set_title("Compressed without heat exchange, a gas stiffens faster")
    ax.legend(fontsize=8)
    fig.savefig(out / "fig_isentrope.png", bbox_inches="tight")

    # --- e(v) at fixed T: perfect gas flat, van der Waals not -----------------------------------------------------------------
    R_co2 = ch01.gas_constant(ch01.MOLAR_MASS["CO2"])
    cv_co2 = ch01.cv_from_gamma(ch01.GAMMA_BY_ATOMICITY["triatomic"], R_co2)
    v = np.logspace(-3, -1, 100)
    fig, ax = plt.subplots()
    for T, c in ((250.0, COLORS["teal"]), (350.0, COLORS["orange"])):
        ax.semilogx(v, ch01.van_der_waals_internal_energy(T, v, ch01.VDW_CO2["a"], cv_co2) / 1e3, color=c,
                    label=f"van der Waals CO₂, T = {T:g} K")
        ax.semilogx(v, np.full_like(v, ch01.perfect_gas_internal_energy(T, cv_co2) / 1e3), color=c, ls="--",
                    label=f"perfect gas, T = {T:g} K")
    ax.set_xlabel("v [m³/kg]")
    ax.set_ylabel("e [kJ/kg]")
    ax.set_title("Only for a perfect gas is e a function of T alone")
    ax.legend(fontsize=8)
    fig.savefig(out / "fig_vdw_energy.png", bbox_inches="tight")

    for k, num, ex, ds in rows:
        print(f"{k:22s} q = {ex['q'] / 1e3:8.3f} kJ/kg (trapezoid {num['q'][-1] / 1e3:8.3f}), w = {ex['w'] / 1e3:8.3f}, "
              f"de = {ex['de'] / 1e3:7.3f}, ds = {ex['ds']:.4f} (∫dq/T numeric {ds:.4f}) J/(kg K)")
    print(f"stirring 300->450 K: q = 0, w = de = {stir['w'] / 1e3:.2f} kJ/kg, ds = {stir['ds']:.2f} > 0")
    print(f"free expansion v doubles: ds = {fexp['ds']:.3f} J/(kg K) = R ln 2 = {R * np.log(2):.3f}")
    cp = ch01.CP_AIR
    print(f"air: R = {R:.4f}, cp = {cp:.3f}, cv = {cv:.3f}, gamma = {ch01.gamma_from_cp(cp):.4f}; "
          f"c(288.15 K) = {ch01.perfect_gas_sound_speed(288.15):.2f} m/s; "
          f"c from EOS derivative = {ch01.sound_speed_from_eos(lambda r, s: ch01.isentropic_pressure(r, ch01.P_ATM, 1.225), 1.225):.2f}")
    print(f"isentropic compression 1 -> 2 bar: T/T0, rho/rho0 = {ch01.isentropic_ratios(2.0)}; from-scratch vs exact "
          f"max rel diff {np.max(np.abs(p_scratch / p_s - 1)):.1e}")
    print(f"Tait water c at rho0 = {ch01.sound_speed_from_eos(lambda r, s: ch01.tait_pressure(r), 1000.0):.1f} m/s")
    print(f"saved 4 figures in {out}  ({time.perf_counter() - t_start:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
