"""§4.8–§4.9 energy and Bernoulli (C10, C11, C12, E6, E7): Couette heating (steady profile, energy budget, transient
heating), the entropy production, the Bernoulli function of a Rankine vortex along vs across circles, E6's scenarios,
the pitot tube, the draining tank, the U-tube and the accelerating sphere (added mass).

Run: ``.venv/Scripts/python.exe scripts/ch04_energy_bernoulli.py --no-show``
Figures → outputs/ch04/c10_couette_heating.png, c11_rankine_bernoulli.png, c05_tank_drain.png.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from fluidpy import ch04_conservation_laws as ch04  # noqa: E402
from fluidpy.core.style import COLORS  # noqa: E402


def couette_figure(U: float = 5.0, h: float = 1e-3, mu: float = 0.3, k: float = 0.15, rho: float = 880.0,
                   cp: float = 1900.0):
    """An oil bearing (our numbers): u(y), ε(y) and T(y, t) heating to the steady parabola."""
    import matplotlib.pyplot as plt

    y = np.linspace(0, h, 101)
    s = ch04.couette_heating(y, U, h, mu, k)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4))
    a1.plot(s["u"], y * 1e3, color=COLORS["teal"], label="u(y) [m/s]")
    a1b = a1.twiny()
    a1b.plot(s["eps"] / 1e6, y * 1e3, color=COLORS["rose"], label="ρε [MW/m³]")
    a1.set_xlabel("u [m/s]")
    a1b.set_xlabel("ρε = μ(du/dy)² [MW/m³]", color=COLORS["rose"])
    a1.set_ylabel("y [mm]")
    tau = h ** 2 * rho * cp / k
    for f in (0.02, 0.05, 0.1, 0.3, 1.0):
        a2.plot(ch04.couette_heating_transient(y, f * tau, U, h, mu, k, rho, cp) - 293.15, y * 1e3,
                label=f"t = {f:g} h²/κ")
    a2.plot(s["T"] - 293.15, y * 1e3, "k--", label="steady")
    a2.set_xlabel("T − T₀ [K]")
    a2.set_ylabel("y [mm]")
    a2.legend(fontsize=7)
    a2.set_title(f"heat out {s['heat_out']:.0f} W/m² = work in {s['work_in']:.0f} W/m²")
    return fig, s


def rankine_figure(Gamma: float = 1.0, sigma: float = 0.1):
    """B = ½u_θ² + p/ρ across circles of a Rankine vortex: varies in the core (rotational), flat outside."""
    import matplotlib.pyplot as plt

    r = np.linspace(0.0, 4 * sigma, 200)
    d = ch04.rankine_bernoulli(r, Gamma, sigma, 1000.0)
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(r / sigma, d["B"], color=COLORS["accent"], label="B across circles")
    ax.plot(r / sigma, 0.5 * np.asarray(d["u_theta"]) ** 2, color=COLORS["teal"], label="½u_θ²")
    ax.plot(r / sigma, np.asarray(d["p"]) / 1000.0, color=COLORS["orange"], label="p/ρ")
    ax.axvline(1.0, color=COLORS["grid"])
    ax.set_xlabel("r/σ")
    ax.set_ylabel("[m²/s²]")
    ax.set_title("(4.71) holds on each circle; (4.72) only outside the core")
    ax.legend(fontsize=8)
    return fig


def tank_figure():
    """Draining tank with a sharp-edged (C_c = 0.611) and a rounded (1.0) orifice."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.5, 4))
    for Cc, c in ((1.0, COLORS["teal"]), (0.611, COLORS["rose"])):
        r = ch04.tank_drain(1.0, 1.0, 1e-3, Cc)
        ax.plot(r["t"] / 60, r["h"], color=c, label=f"C_c = {Cc}: empty after {r['t_empty'] / 60:.1f} min")
    ax.set_xlabel("t [min]")
    ax.set_ylabel("h [m]")
    ax.set_title("quasi-steady Torricelli: dh/dt = −(C_cA_o/A)√(2gh)")
    ax.legend(fontsize=8)
    return fig


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch04"))
    ap.add_argument("--no-show", action="store_true")
    args = ap.parse_args()
    t0 = time.perf_counter()
    import matplotlib

    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from fluidpy.core.style import use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    fig, s = couette_figure()
    fig.savefig(out / "c10_couette_heating.png", bbox_inches="tight")
    rankine_figure().savefig(out / "c11_rankine_bernoulli.png", bbox_inches="tight")
    tank_figure().savefig(out / "c05_tank_drain.png", bbox_inches="tight")

    c = ch04.couette_heating(np.linspace(0, 1e-3, 5), 1.0, 1e-3, 1e-3, 0.6)
    print(f"C10 water Couette U = 1 m/s, h = 1 mm: ΔT_max {c['dT_max']:.4e} K (μU²/8k = {1e-3 / 4.8:.4e}), heat out {c['heat_out']:.4f}"
          f" = work in {c['work_in']:.4f} W/m²")
    print(f"C10 oil bearing: ΔT_max {s['dT_max']:.2f} K, heat out {s['heat_out']:.1f} W/m²")
    u, _ = ch04.exact_field("couette", U=1.0, h=1e-3, mu=1e-3)
    Tf = lambda X, t: ch04.couette_heating(X[1], 1.0, 1e-3, 1e-3, 0.6)["T"]  # noqa: E731
    res = ch04.internal_energy_residual(1000.0, u, lambda X, t: 4182.0 * Tf(X, t), 0.0, Tf, 1e-3, 0.0, 0.6,
                                        np.array([[0.0, 0.0], [2e-4, 7e-4]]), h=1e-6)
    print(f"C10 (4.60) residual of the steady Couette temperature: {np.abs(res).max():.1e} W/m³ (vs ρε = {1e-3 * 1e6:.0f})")
    print(f"C10 entropy production (∇T = 100 K/m, ε = 1 W/kg): {ch04.entropy_production(100.0, 300.0, 0.6, 1000.0, 1.0):.4e} W/(kg K); "
          f"negative k: {ch04.entropy_production(100.0, 300.0, -0.6, 1000.0, 0.0):.2e} (second law violated)")
    print("C11 energy identities (4.53)…(4.112):", ch04.energy_identities_sym())
    for n in ch04.BERNOULLI_SCENARIOS:
        b = ch04.bernoulli_scenario(n)
        print(f"E6 {n:11s}: valid {b['valid']}, status '{b['status']}', ptp along {np.ptp(b['B_along']):.2e}, "
              f"across {np.ptp(b['B_across']):.2e}")
    print(f"C05 pitot 500 Pa in air: {ch04.pitot_speed(500.0, 0.0, 1.2):.3f} m/s; Torricelli 1 m: {ch04.torricelli_speed(1.0):.4f} m/s")
    ut = ch04.u_tube_column(0.3)
    print(f"C12 U-tube L = 1 m: ω = {ut['omega']:.4f} rad/s, dp = {ut['dp']:.2f} Pa, (4.82) residual {ut['residual_4_82']:.1e}")
    sp_ = ch04.accelerating_sphere_pressure(0.0, 0.1, 2.0)
    print(f"C12 sphere a = 0.1 m, dU/dt = 2: force {sp_['force']:.5f} N vs −m_a dU/dt {-sp_['added_mass'] * 2:.5f} N")
    phi, pf = ch04.accelerating_sphere_fields()
    Xs = np.array([[0.3, 0.5, -0.4], [0.2, -0.3, 0.1], [0.1, 0.2, 0.3]])
    B = ch04.unsteady_bernoulli_B(phi, pf, Xs, 0.4)
    print(f"C12 (4.74) bracket at 3 points: {B} (uniform); gauge: −∫B → {ch04.gauge_absorbed_bracket(2.0, -1)}, "
          f"book's + → {ch04.gauge_absorbed_bracket(2.0, +1)}")
    print(f"N96 stagnation temperature 300 K, 100 m/s: {ch04.stagnation_temperature(300.0, 100.0):.2f} K")
    print(f"figures → {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
