"""§4.2–§4.4: mass and momentum budgets — the expanding flow seen by a fixed box and by a material interval (C01, C02),
the wake drag of Example 4.1 vs box height, the bore of Example 4.3 (with the √(gh) limit and Bélanger's relation),
the rocket of Example 4.4 vs Tsiolkovsky, E1's four control-volume scenarios and the sprinkler of Example 4.6.

Run: ``.venv/Scripts/python.exe scripts/ch04_mass_momentum.py --no-show``
Figures → outputs/ch04/c01_expanding_flow_budget.png, c04_wake_drag_vs_H.png, c04_bore_speed.png, c04_rocket.png.
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
from tools.convergence import observed_order  # noqa: E402


def expanding_budget_figure(a: float = 1.0, rho0: float = 1.0):
    """C01: mass in a fixed interval [1, 2] (falls) vs the material interval that started there (constant), with the
    interval budget terms of (4.5) at t = 0."""
    import matplotlib.pyplot as plt

    ts = np.linspace(0.0, 2.0, 81)
    fixed = [ch04.interval_mass_budget(1.0, 2.0, t)["storage"] for t in ts]
    m_fixed = [rho0 / (1 + a * t) * 1.0 for t in ts]
    m_mat = [ch04.material_mass(1.0, 2.0, t, a, rho0) for t in ts]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4))
    a1.plot(ts, m_fixed, color=COLORS["orange"], label="fixed interval [1, 2] m")
    a1.plot(ts, m_mat, color=COLORS["accent"], label="material interval (moves with the fluid)")
    a1.set_xlabel("t [s]")
    a1.set_ylabel("mass per unit area [kg/m²]")
    a1.set_title("(4.1) vs (4.5): same flow, two boxes")
    a1.legend(fontsize=8)
    r = ch04.interval_mass_budget(1.0, 2.0, 0.0)
    a2.bar(["storage", "−flux in (left)", "flux out (right)", "residual"],
           [r["storage"], -r["flux_left"], r["flux_right"], r["residual"]],
           color=[COLORS["accent"], COLORS["blue"], COLORS["orange"], COLORS["ink"]])
    a2.set_ylabel("[kg/(m² s)]")
    a2.set_title("fixed interval at t = 0: storage + net outflux = 0")
    return fig, fixed


def wake_figure():
    """C04 / Example 4.1: F_D/l from the wake survey vs box height H, converging to the closed form."""
    import matplotlib.pyplot as plt

    Hs = np.linspace(0.05, 1.0, 40)
    F = [ch04.wake_drag_per_span("gaussian", 10.0, 1.2, H) for H in Hs]
    Fc = ch04.wake_drag_gaussian(10.0, 2.0, 0.1, 1.2)
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(Hs, F, color=COLORS["accent"], label=r"$\rho\int U(U_\infty-U)\,dy$ over the box")
    ax.axhline(Fc, color=COLORS["orange"], ls="--", label=f"H → ∞ closed form {Fc:.4f} N/m")
    ax.set_xlabel("box height H [m] (wake half-width b = 0.1 m)")
    ax.set_ylabel(r"$F_D/l$ [N/m]")
    ax.set_title("Example 4.1: a box tall enough captures the whole deficit")
    ax.legend(fontsize=8)
    return fig, Fc


def bore_figure(g: float = 9.81):
    """C04 / Example 4.3: bore speed vs depth ratio with the √(gh) limit."""
    import matplotlib.pyplot as plt

    r = np.linspace(1.0, 2.0, 60)
    U = ch04.bore_speed(1.0, r, g)
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(r, U, color=COLORS["accent"], label=r"$U=\sqrt{g h_{out}(h_{in}+h_{out})/2h_{in}}$")
    ax.axhline(np.sqrt(g), color=COLORS["orange"], ls="--", label=r"$\sqrt{gh}$ (small step)")
    ax.set_xlabel(r"$h_{out}/h_{in}$ ($h_{in}$ = 1 m)")
    ax.set_ylabel("U [m/s]")
    ax.set_title("Example 4.3: a bigger step runs faster")
    ax.legend(fontsize=8)
    return fig


def rocket_figure():
    """Example 4.4: speed during the burn with and without gravity, and Tsiolkovsky's Δv."""
    import matplotlib.pyplot as plt

    M0, mdot, Ve, tb = 1.0, 0.05, 500.0, 10.0
    r0 = ch04.rocket_trajectory(M0, mdot, Ve, tb, g=0.0)
    r1 = ch04.rocket_trajectory(M0, mdot, Ve, tb)
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(r0["t"], r0["b"], color=COLORS["accent"], label="no gravity (ODE)")
    ax.plot(r1["t"], r1["b"], color=COLORS["teal"], label="g = 9.81 m/s² (ODE)")
    ax.plot(r0["t"], ch04.rocket_delta_v(M0, r0["M"], Ve), ":", color=COLORS["orange"], label=r"$V_e\ln(M_0/M)$")
    ax.set_xlabel("t [s]")
    ax.set_ylabel("b [m/s]")
    ax.set_title("Example 4.4: M d²z/dt² = −V_e dM/dt − Mg")
    ax.legend(fontsize=8)
    return fig, r0, r1


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

    fig, _ = expanding_budget_figure()
    fig.savefig(out / "c01_expanding_flow_budget.png", bbox_inches="tight")
    print("C01 interval budget at t = 0:", ch04.interval_mass_budget(1.0, 2.0, 0.0))
    print("C01 material mass at t = 0, 1, 2:", [round(ch04.material_mass(1.0, 2.0, t), 12) for t in (0, 1, 2)])
    rho, u = ch04.expanding_flow_fields(1.0, 1.0, 3)
    box = ch04.MovingBox(lengths=(1, 1, 1), origin=(0.2, 0.1, 0.3))
    mb = ch04.mass_budget(rho, u, box, 0.5)
    print(f"C01 3-D fixed box: storage {mb.storage:+.8f}, outflux {mb.outflux:+.8f}, residual {mb.residual:+.1e} kg/s")
    X = np.array([[0.3, 0.7], [0.1, -0.4], [0.5, 0.2]])
    hs = np.array([0.1, 0.05, 0.025])
    wrong = lambda x, t: 1.0 / (1 + t) ** 2 + 0 * x[0]  # noqa: E731  (dim = 3 needs the cube)
    print("C02 continuity residual (exact, wrong ρ):", ch04.continuity_residual(rho, u, X, 0.5),
          ch04.continuity_residual(wrong, u, X, 0.5))
    errs = [abs(ch04.continuity_residual(lambda x, t: np.sin(x[0]) + 0 * t, lambda x, t: np.stack([np.cos(x[0]) * 0 + x[1],
            0 * x[0], 0 * x[0]]), np.array([0.3, 0.2, 0.1]), 0.0, h=h) - 0.2 * np.cos(0.3)) for h in hs]
    print(f"C02 stencil order (∂(ρu)/∂x of a smooth field): {observed_order(hs, errs):.2f}")

    fig, Fc = wake_figure()
    fig.savefig(out / "c04_wake_drag_vs_H.png", bbox_inches="tight")
    print(f"Ex. 4.1 wake drag (U∞ = 10, Δ = 2, b = 0.1, ρ = 1.2, H = 2): {ch04.wake_drag_per_span('gaussian', 10.0, 1.2, 2.0):.4f}"
          f" N/m, closed form {Fc:.4f}, side outflow {ch04.wake_side_outflow('gaussian', 10.0, 2.0):.5f} m²/s")
    print("Ex. 4.2 sympy result:", ch04.stream_tube_element_balance_sym()["result"])
    bore_figure().savefig(out / "c04_bore_speed.png", bbox_inches="tight")
    hin, hout = 1.0, 1.1
    U = ch04.bore_speed(hin, hout)
    Fr1 = U / np.sqrt(9.81 * hin)
    print(f"Ex. 4.3 bore U = {U:.4f} m/s; Bélanger h2/h1 from Fr1 = {Fr1:.4f}: {(np.sqrt(1 + 8 * Fr1 ** 2) - 1) / 2:.6f}"
          f" (vs {hout / hin}); p_o independence: net = {ch04.bore_pressure_force(hin, hout)['net']:.3f} N, "
          f"{ch04.bore_pressure_force(hin, hout, p_o=0.0)['net']:.3f} N")
    fig, r0, r1 = rocket_figure()
    fig.savefig(out / "c04_rocket.png", bbox_inches="tight")
    print(f"Ex. 4.4 no-g Δb {r0['b'][-1]:.6f} vs Tsiolkovsky {ch04.rocket_delta_v(1.0, 0.5, 500.0):.6f} m/s; with g "
          f"z(t_b) {r1['z'][-1]:.4f} vs closed form {ch04.rocket_closed_form(10.0, 1.0, 0.05, 500.0)[0]:.4f} m")
    for n in ch04.CV_SCENARIOS:
        s = ch04.cv_scenario(n)
        print(f"E1 {n:7s}: result {s['result']:.5g} ({s['result_label']}), residuals mass {s['residual_mass']:+.1e}, "
              f"momentum {s['residual_momentum']:+.1e}")
    print(f"Ex. 4.6 sprinkler torque {ch04.sprinkler_torque(0.2, 1000.0, 1e-4, 5.0, np.pi / 6):.6f} N m, numeric "
          f"{ch04.sprinkler_torque_numeric(0.2, 1000.0, 1e-4, 5.0, np.pi / 6):.6f}; free spin "
          f"{ch04.sprinkler_free_spin_rate(0.2, 5.0, np.pi / 6):.3f} rad/s")
    print(f"figures → {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
