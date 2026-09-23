"""§4.9 Boussinesq and §4.10 boundary conditions (C13, C14, E8, B1): the validity numbers of E8's scenarios, a warm
blob carried and spread by (4.89), the buoyancy/N² link, the kinematic condition of a linear wave (residual O((ka)²)),
the pillbox limit, the curved-cap force balance converging to Laplace's (1.5), and the meniscus of Example 4.7 by the
closed form and by an independent ODE.

Run: ``.venv/Scripts/python.exe scripts/ch04_boussinesq_interfaces.py --no-show``
Figures → outputs/ch04/c13_boussinesq_validity.png, c13_blob.png, c14_kinematic_residual.png, c14_cap_convergence.png,
c14_meniscus.png.
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


def validity_figure():
    """αδT, L/(c²/g) and the viscous-heating ratio for E8's five scenarios on log axes, with the 0.1 threshold."""
    import matplotlib.pyplot as plt

    names = list(ch04.BOUSSINESQ_SCENARIOS)
    keys = [("alpha_dT", "αδT"), ("L_over_Hc", "L/(c²/g)"), ("heating_ratio", "νU/(C_pδT L)")]
    fig, ax = plt.subplots(figsize=(9, 4))
    w = 0.25
    for j, (k, lab) in enumerate(keys):
        vals = []
        for n in names:
            s = ch04.boussinesq_scenario(n)
            v = ch04.boussinesq_validity(s["alpha"], s["dT"], s["L"], s["U"], s["c"], nu=s["nu"], cp=s["cp"])
            vals.append(v[k])
        ax.bar(np.arange(len(names)) + (j - 1) * w, vals, w, label=lab)
    ax.axhline(0.1, color=COLORS["rose"], ls="--", label="≪ 1 threshold (0.1)")
    ax.set_yscale("log")
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=15)
    ax.set_title("Boussinesq validity (§4.9) — the deep atmosphere fails L ≪ c²/g")
    ax.legend(fontsize=8)
    return fig


def blob_figure(U: float = 0.01, kappa: float = 1.4e-7 * 1000, sigma0: float = 0.01):
    """A warm blob carried by a current and spreading, σ² = σ₀² + 2κt (exact solution of (4.89))."""
    import matplotlib.pyplot as plt

    x = np.linspace(-0.05, 0.25, 301)
    fig, ax = plt.subplots(figsize=(7, 4))
    for t in (0.0, 5.0, 10.0, 20.0):
        ax.plot(x, ch04.gaussian_blob_advection_diffusion(x, t, U, kappa, sigma0, 1.0, dim=1), label=f"t = {t:g} s")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("T′ [K]")
    ax.set_title(f"DT/Dt = κ∇²T: moves at U = {U} m/s, spreads with κ = {kappa:.1e} m²/s")
    ax.legend(fontsize=8)
    return fig


def kinematic_figure():
    """B1/C14: the full kinematic condition on z = η has residual O((ka)²) (slope 2); the linearised one is exact."""
    import matplotlib.pyplot as plt

    ka = np.logspace(-3, -1, 8)
    x = np.linspace(0, 2 * np.pi, 13)
    full = [np.abs(ch04.wave_kinematic_residual(x, 0.3, a=a, k=1.0)).max() for a in ka]
    lin = [np.abs(ch04.wave_kinematic_residual(x, 0.3, a=a, k=1.0, full=False)).max() for a in ka]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.loglog(ka, full, "o-", color=COLORS["accent"], label=f"full (4.91) on z = η: slope {observed_order(ka, full):.2f}")
    ax.loglog(ka, np.maximum(lin, 1e-18), "s-", color=COLORS["teal"], label="linearised at z = 0 (exactly 0)")
    ax.set_xlabel("ka")
    ax.set_ylabel("max |Dη/Dt| [m/s]")
    ax.legend(fontsize=8)
    return fig, observed_order(ka, full)


def cap_figure(sigma: float = 0.0728, R1: float = 1e-3, R2: float = 2e-3):
    """C14: the exact rim integral of σ∮t × n ds approaches the book's small-ζ form as ζ → 0 (error O(ζ/R), slope 1)."""
    import matplotlib.pyplot as plt

    zeta = np.logspace(-9, -5, 9)
    err = [abs(ch04.cap_surface_tension_force(sigma, R1, R2, z) / ch04.cap_surface_tension_force(sigma, R1, R2, z, exact=False)
               - 1.0) for z in zeta]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.loglog(zeta, err, "o-", color=COLORS["rose"])
    ax.set_xlabel("cap height ζ [m]")
    ax.set_ylabel("|F_st,exact/F_st,book − 1|")
    ax.set_title(f"(4.98): slope {observed_order(zeta, err):.2f}")
    return fig, observed_order(zeta, err)


def meniscus_figure(sigma: float = 0.0728, rho: float = 998.0):
    """Example 4.7: profiles for several contact angles — closed form (lines) and curvature ODE (dots)."""
    import matplotlib.pyplot as plt

    d = np.sqrt(sigma / (rho * 9.81))
    fig, ax = plt.subplots(figsize=(7, 4))
    worst = 0.0
    for deg in (0, 20, 45, 70):
        th = np.radians(deg)
        x, z = ch04.meniscus_profile_ode(th, 5 * d, sigma, rho)
        ax.plot(x / d, z / d, ".", ms=3)
        zz = np.linspace(z[-1] * 1.0001, z[0], 200) if z[0] > 0 else np.array([])
        if zz.size:
            ax.plot(np.asarray(ch04.meniscus_profile_x(zz, th, sigma, rho)) / d, zz / d, "-", label=f"θ = {deg}°")
            xc = np.asarray(ch04.meniscus_profile_x(z[1:150], th, sigma, rho))
            worst = max(worst, float(np.max(np.abs(xc - x[1:150]) / d)))
    ax.set_xlabel("x/δ")
    ax.set_ylabel("ζ/δ")
    ax.set_title(r"Example 4.7: $h^2 = 2\delta^2(1-\sin\theta)$, $\delta^2 = \sigma/\rho g$")
    ax.legend(fontsize=8)
    return fig, worst


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
    validity_figure().savefig(out / "c13_boussinesq_validity.png", bbox_inches="tight")
    blob_figure().savefig(out / "c13_blob.png", bbox_inches="tight")
    fig, slope_k = kinematic_figure()
    fig.savefig(out / "c14_kinematic_residual.png", bbox_inches="tight")
    fig, slope_c = cap_figure()
    fig.savefig(out / "c14_cap_convergence.png", bbox_inches="tight")
    fig, worst = meniscus_figure()
    fig.savefig(out / "c14_meniscus.png", bbox_inches="tight")

    for n in ch04.BOUSSINESQ_SCENARIOS:
        s = ch04.boussinesq_scenario(n)
        v = ch04.boussinesq_validity(s["alpha"], s["dT"], s["L"], s["U"], s["c"], nu=s["nu"], cp=s["cp"])
        print(f"E8 {n:16s}: αδT {v['alpha_dT']:.1e}, L/(c²/g) {v['L_over_Hc']:.1e}, heating {v['heating_ratio']:.1e} → {v['verdict']}")
    print(f"N108 air at 288 K: c²/g = {ch04.boussinesq_validity(1 / 288, 1, 1, 1, c=np.sqrt(1.4 * 287.05 * 288.0))['H_c'] / 1e3:.2f} km")
    T, u = ch04.gaussian_blob_fields(U=0.01, kappa=1e-4, sigma0=0.01, dim=2)
    X = np.array([[0.01, 0.03], [0.0, -0.01]])
    print(f"C13 (4.89) residual of the moving blob: {np.abs(ch04.temperature_equation_residual(T, u, X, 2.0, kappa=1e-4, h=1e-4, ht=1e-3)).max():.1e} K/s")
    print(f"N137 buoyancy b = −gρ′/ρ₀ for ρ′ = −0.2 kg/m³: {ch04.buoyancy(-0.2, 1000.0):.5f} m/s²; boussinesq_density(303.15 K) = "
          f"{ch04.boussinesq_density(303.15):.4f} kg/m³")
    print(f"C14 kinematic residual slope in ka: {slope_k:.2f}; moving wall: {ch04.kinematic_bc_residual('moving_wall', None, np.array([2.0, 0.3]), 2.0, V=0.5)}")
    print(f"N116 pillbox (side 0.3, volume 2.0 per m): {ch04.pillbox_limit(1.0, 1.0, 0.3, 2.0, np.array([1e-1, 1e-2, 1e-3]))}")
    print(f"C14 cap: exact/book − 1 slope {slope_c:.2f}; Laplace jump (book form) {ch04.laplace_jump_from_balance(0.0728, 1e-3, 1e-3):.4f} Pa, "
          f"exact rim ζ = 1e-9 m {ch04.laplace_jump_from_balance(0.0728, 1e-3, 1e-3, zeta=1e-9, exact=True):.4f} Pa")
    print(f"C14 pressure force: dblquad {ch04.cap_pressure_force(10.0, 1e-3, 2e-3, 1e-5):.6e} N vs closed "
          f"{ch04.cap_pressure_force(10.0, 1e-3, 2e-3, 1e-5, exact=False):.6e} N")
    print(f"Ex. 4.7 meniscus: closed form vs curvature ODE max |Δx|/δ = {worst:.1e}; h(θ = 0) = {ch04.meniscus_height(0.0):.5f} m, "
          f"capillary length {ch04.capillary_length(0.0728, 998.0) * 1e3:.4f} mm")
    print(f"N123 spheroid area at fixed volume: aspect 0.6/1/1.5 → {ch04.spheroid_area(1.0, 0.6):.4f}, {ch04.spheroid_area(1.0, 1.0):.4f}, "
          f"{ch04.spheroid_area(1.0, 1.5):.4f} m²")
    print(f"figures → {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
