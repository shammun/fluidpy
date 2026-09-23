"""§5.1 pressure and viscous stress in the basic vortices (C02, R02, R03, N05–N09, R07, R08): the rotating tank's
paraboloids (5.6) beside the line vortex's funnel (5.7) and a Rankine tornado; the Bernoulli function across the
streamlines; stress without net force (three routes, D03); the rotating cylinder's torque carried to infinity and its
energy budget; the radial balance (5.5a) of every E2 scenario; (5.5a) integrated outward by hand.

Run: ``.venv/Scripts/python.exe scripts/ch05_vortex_pressure.py --no-show``
Figures → outputs/ch05/c02_isobars.png, c02_profiles.png.
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


def isobar_figure():
    """(r, z) sections: tank paraboloids (ω = 10 s⁻¹) and the Rankine funnel (Γ = 1 m²/s, a = 0.02 m)."""
    import matplotlib.pyplot as plt
    from ch05_drawings import tank_section

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    R, H = 0.1, 0.15
    tank_section(ax[0], R, H, "rotating tank, ω = 10 s⁻¹ (turning at 5 rad/s)")
    r = np.linspace(-R, R, 201)
    ts = ch05.rotating_tank_free_surface(R, 0.1, 5.0)
    for k, dp in enumerate((0.0, -0.03, -0.06)):
        z = ts["z_vertex"] + np.asarray(ch05.isobar_height(r, 10.0, dp))
        ax[0].plot(r, z, color=COLORS["blue"] if k == 0 else COLORS["orange"], ls="-" if k == 0 else "--",
                   label="free surface" if k == 0 else f"isobar, Δp/ρg = {dp:+.2f} m")
    ax[0].legend(fontsize=8)
    r2 = np.concatenate([-np.geomspace(0.3, 1e-3, 200), np.geomspace(1e-3, 0.3, 200)])
    for dp in (0.0, -0.02, -0.04):
        ax[1].plot(r2, np.asarray(ch05.isobar_height(np.abs(r2), 1.0, dp, "rankine", a=0.02)),
                   color=COLORS["blue"] if dp == 0 else COLORS["orange"], ls="-" if dp == 0 else "--")
        ax[1].plot(r2, np.asarray(ch05.isobar_height(np.abs(r2), 1.0, dp, "line")), color=COLORS["muted"], lw=0.8)
    ax[1].set_ylim(-0.4, 0.05)
    ax[1].set_xlabel("r [m]")
    ax[1].set_ylabel("z [m]")
    ax[1].set_title("line vortex funnel (grey) and Rankine core (Γ = 1 m²/s)")
    return fig


def profiles_figure():
    """u_θ, B − B_ref and σ_rθ across the four E2 vortices."""
    import matplotlib.pyplot as plt

    r = np.linspace(0.005, 0.3, 200)
    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    for kind, col in zip(ch05.VORTEX_PRESSURE_KINDS, (COLORS["accent"], COLORS["teal"], COLORS["orange"],
                                                       COLORS["blue"])):
        d = [ch05.vortex_pressure_scenario(kind, rr) for rr in r]
        ax[0].plot(r, [q["u_theta"] for q in d], color=col, label=kind)
        ax[1].plot(r, [q["B"] for q in d], color=col)
        ax[2].plot(r, [q["sigma_rtheta"] for q in d], color=col)
    ax[0].set_ylabel("u_θ [m/s]")
    ax[1].set_ylabel("B − B_ref [J/kg]")
    ax[2].set_ylabel("σ_rθ [Pa]")
    for a in ax:
        a.set_xlabel("r [m]")
    ax[0].legend(fontsize=8)
    ax[1].set_title("B grows outward only where the flow is rotational")
    ax[2].set_title("stress everywhere outside the rigid cores")
    return fig


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch05"))
    ap.add_argument("--no-show", action="store_true")
    args = ap.parse_args()
    import matplotlib

    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.integrate import cumulative_trapezoid

    from fluidpy.core.style import use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    print(f"(5.6) p − p_o at r = 0.1 m, ω = 10 s⁻¹: {ch05.solid_body_pressure(0.1, 0.0, 10.0):.3f} Pa; "
          f"rim − centre height ω²R²/8g = {ch05.isobar_height(0.1, 10.0):.6f} m")
    ts = ch05.rotating_tank_free_surface(0.1, 0.1, 5.0)
    print(f"tank R = 0.1 m, depth 0.1 m, Ω = 5 rad/s: vertex {ts['z_vertex']:.6f} m, rim {ts['z_rim']:.6f} m "
          f"(rise {ts['rim_rise'] * 1e3:.3f} mm = drop {ts['centre_drop'] * 1e3:.3f} mm), volume residual "
          f"{ts['volume_residual']:.1e} m³")
    ts = ch05.rotating_tank_free_surface(0.05, 0.05, 30.0, H=0.08, closed=True)
    print(f"closed tank R = 5 cm, depth 5 cm, H = 8 cm, 30 rad/s: dry bottom radius {ts['r_dry'] * 100:.3f} cm, "
          f"uncovered {ts['uncovered_area'] * 1e4:.3f} cm², meets the lid at r = {ts['r_top'] * 100:.3f} cm")
    print(f"(5.7) Γ = 1 m²/s, r = 0.1 m: p − p_∞ = {ch05.line_vortex_pressure(0.1, 0.0, 1.0):.3f} Pa, funnel "
          f"depth {ch05.isobar_height(0.1, 1.0, kind='line'):.6f} m")
    tor = ch05.vortex_pressure_scenario("rankine", 50.0, rho=1.2, Gamma=1e4, a=50.0)
    print(f"tornado preset (Γ = 1e4 m²/s, a = 50 m, ρ = 1.2): u_max {tor['u_theta']:.3f} m/s, core-edge deficit "
          f"{-tor['p']:.2f} Pa, axis deficit {-ch05.rankine_pressure(0.0, 0.0, 1e4, 50.0, 1.2):.2f} Pa; Γ back from "
          f"the edge pressure {ch05.tornado_circulation(tor['p'], 50.0, 1.2):.1f} m²/s")
    # (5.5a) integrated outward by hand vs (5.6), (5.7)
    r = np.linspace(0.05, 0.3, 2001)
    for kind, param, fn in (("solid", 10.0, lambda q: ch05.solid_body_pressure(q, 0.0, 10.0)),
                            ("line", 1.0, lambda q: ch05.line_vortex_pressure(q, 0.0, 1.0))):
        u = 0.5 * param * r if kind == "solid" else param / (2 * np.pi * r)
        p = fn(r[0]) + cumulative_trapezoid(1000.0 * u ** 2 / r, r, initial=0.0)  # ∂p/∂r = ρu_θ²/r  (5.5a)
        print(f"(5.5a) integrated by hand, {kind:5s}: max |p − formula| = {np.abs(p - fn(r)).max():.2e} Pa")
    print("stress vs net force (N09), r = 5 cm, Γ = 0.01 m²/s:")
    for row in ch05.stress_vs_net_force_table():
        print(f"   {row['name']:12s} σ_rθ = {row['sigma_rtheta']:+.3e} Pa, net force = {row['net_force']:+.3e} N/m³")
    for kind in ("solid", "line", "gaussian"):
        v = ch05.vortex_stress_force(kind, 0.2)
        print(f"   D03 routes {kind:9s}: div {v['force_divergence']:+.2e}, lap {v['force_laplacian']:+.2e}, "
              f"curl {v['force_curl']:+.2e} N/m³")
    print(f"torque per length −2μΓ at r = 0.1, 1, 10 m: {np.round(ch05.torque_per_length([0.1, 1, 10], 1.0), 6)}")
    d = ch05.dissipation_outside_cylinder(0.1, 1.0, 1.0)
    print(f"cylinder energy budget: dissipated {d['dissipated']:.7f} W/m = in {d['power_in']:.7f} − out "
          f"{d['power_out']:.7f} (residual {d['residual']:.1e})")
    for kind in ch05.VORTEX_PRESSURE_KINDS:
        worst = max(abs(ch05.vortex_pressure_scenario(kind, rr)["radial_residual"]) /
                    max(ch05.vortex_pressure_scenario(kind, rr)["centripetal"], 1e-12) for rr in (0.03, 0.1, 0.3))
        print(f"   E2 {kind:8s}: (5.5a) balance relative residual ≤ {worst:.1e}")
    isobar_figure().savefig(out / "c02_isobars.png", bbox_inches="tight")
    profiles_figure().savefig(out / "c02_profiles.png", bbox_inches="tight")
    print(f"figures → {out}")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
