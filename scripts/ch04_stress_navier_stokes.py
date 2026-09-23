"""§4.4–§4.6: Cauchy's equation, the Newtonian stress and Navier–Stokes (C06, C07, C08, E3, E4): stress on a rotating
plane for E3's presets (and the μ_v effect on expansion), the spinning-cube argument for τ₁₂ = τ₂₁ (slope −2), the
five terms of (4.39b) for the exact solutions, Stokes' first problem, and the viscous force three ways (4.40).

Run: ``.venv/Scripts/python.exe scripts/ch04_stress_navier_stokes.py --no-show``
Figures → outputs/ch04/c07_stress_on_plane.png, c07_cube_spin.png, c08_ns_term_bars.png, c08_stokes_first.png.
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


def stress_plane_figure(mu: float = 1e-3, p: float = 0.0):
    """σ_n and τ_s on a plane of normal angle θ for the four presets (rate 10 s⁻¹), and the expansion with μ_v = 0, 2μ."""
    import matplotlib.pyplot as plt

    th = np.linspace(0, np.pi, 181)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4))
    cols = dict(shear=COLORS["accent"], extension=COLORS["teal"], rotation=COLORS["muted"], expansion=COLORS["orange"])
    for name in ch04.STRESS_LAB_PRESETS:
        G = ch04.stress_lab_gradient(name, 10.0)
        sn, ts = np.array([ch04.stress_on_plane(G, p, mu, 0.0, t) for t in th]).T
        a1.plot(np.degrees(th), sn, color=cols[name], label=f"{name}: σ_n")
        a1.plot(np.degrees(th), ts, "--", color=cols[name], label=f"{name}: τ_s")
    a1.set_xlabel("plane normal angle θ [°]")
    a1.set_ylabel("stress [Pa]")
    a1.set_title("(4.37) on a rotating plane (μ = 1 mPa s, rate 10/s)")
    a1.legend(fontsize=7, ncol=2)
    G = ch04.stress_lab_gradient("expansion", 10.0)
    for mv, c in ((0.0, COLORS["accent"]), (2e-3, COLORS["rose"])):
        tau = ch04.newtonian_stress(G, 100.0, mu, mu_v=mv)
        a2.bar(f"μ_v = {mv * 1e3:.0f} mPa s", 100.0 - ch04.mean_pressure(tau), color=c)
    a2.set_ylabel("p − p̄ [Pa]")
    a2.set_title("(4.34): p − p̄ = μ_v∇·u (expansion ∇·u = 30/s)")
    return fig


def cube_figure():
    """D08: the angular acceleration of a cube with τ₁₂ ≠ τ₂₁ grows as h⁻² (slope −2) — so τ must be symmetric."""
    import matplotlib.pyplot as plt

    h = np.logspace(-4, -1, 20)
    a = ch04.cube_spin_acceleration(1.0, 0.9, 1000.0, h)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.loglog(h, a, color=COLORS["accent"])
    ax.set_xlabel("cube side h [m]")
    ax.set_ylabel(r"$d\Omega_3/dt$ [rad/s²]")
    ax.set_title(f"τ₁₂ − τ₂₁ = 0.1 Pa: slope {observed_order(h, a):.2f}")
    return fig


def term_bars_figure():
    """E4: the five terms of (4.39b) (x component, per unit mass) at a probe for five exact solutions."""
    import matplotlib.pyplot as plt

    cases = [("poiseuille", 0.0, 2.5e-3, 0.0, {}), ("couette", 0.0, 5e-3, 0.0, dict(G=50.0)),
             ("stokes_first", 0.0, 1e-3, 1.0, {}), ("taylor_green", 0.3, 0.7, 5.0, dict(mu=1e-2, rho=1.0)),
             ("lamb_oseen", 0.02, 0.01, 10.0, dict(mu=1e-3, rho=1.0))]
    keys = ["local", "advective", "pressure", "gravity", "viscous", "residual"]
    cols = [COLORS["blue"], COLORS["teal"], COLORS["orange"], COLORS["muted"], COLORS["rose"], COLORS["ink"]]
    fig, axs = plt.subplots(1, len(cases), figsize=(16, 3.6))
    out = {}
    for ax, (name, x, y, t, kw) in zip(axs, cases):
        d = ch04.ns_terms_preset(name, x, y, t, component=0, **kw)
        out[name] = d
        ax.bar(range(6), [d[k] for k in keys], color=cols)
        ax.set_xticks(range(6))
        ax.set_xticklabels(keys, rotation=45, ha="right", fontsize=7)
        ax.set_title(name, fontsize=9)
    axs[0].set_ylabel("m/s²")
    return fig, out


def stokes_figure(U: float = 1.0):
    """C08/C15: Stokes' first problem at several times and two fluids — one curve u/U vs η = y/√(νt)."""
    import matplotlib.pyplot as plt

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4))
    y = np.linspace(0, 0.02, 200)
    for nu, ls in ((1e-6, "-"), (1.5e-5, "--")):
        for t in (1.0, 10.0, 60.0):
            u = ch04.stokes_first_problem(y, t, U, nu)
            a1.plot(u, y * 1e3, ls, label=f"ν = {nu:.1e}, t = {t:g} s")
            a2.plot(u / U, y / np.sqrt(nu * t), ls, lw=1)
    a1.set_xlabel("u [m/s]")
    a1.set_ylabel("y [mm]")
    a1.legend(fontsize=7)
    a1.set_title("u = U erfc(y/2√(νt))")
    a2.set_ylim(0, 6)
    a2.set_xlabel("u/U")
    a2.set_ylabel(r"$\eta = y/\sqrt{\nu t}$")
    a2.set_title("all curves collapse (similarity)")
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
    stress_plane_figure().savefig(out / "c07_stress_on_plane.png", bbox_inches="tight")
    cube_figure().savefig(out / "c07_cube_spin.png", bbox_inches="tight")
    fig, tb = term_bars_figure()
    fig.savefig(out / "c08_ns_term_bars.png", bbox_inches="tight")
    stokes_figure().savefig(out / "c08_stokes_first.png", bbox_inches="tight")

    G = [[0, 10, 0], [0, 0, 0], [0, 0, 0]]
    print(f"C07 shear 10/s, μ = 1e-3: τ₁₂ = {ch04.newtonian_stress(G, mu=1e-3)[0, 1]:.4f} Pa; plane at 45°: "
          f"{ch04.stress_on_plane(G, 0.0, 1e-3, 0.0, np.pi / 4)}")
    K = ch04.isotropic_fourth_order(1.0, 2.0, 2.0)
    from fluidpy.core.tensors import is_isotropic
    print("D09 K(λ, μ, γ) isotropic under 50 rotations (ok, max residual):", is_isotropic(K))
    print(f"C07 rigid rotation: max |σ| = {np.abs(ch04.viscous_stress(ch04.stress_lab_gradient('rotation', 5.0), 1e-3)).max():.1e} Pa")
    print(f"C10 ε for shear 1000/s in water-like fluid: {ch04.dissipation_rate([[0, 1000, 0], [0, 0, 0], [0, 0, 0]], 1000.0, 1e-3):.4f} W/kg")
    for name, d in tb.items():
        print(f"E4 {name:13s}: " + ", ".join(f"{k} {v:+.3e}" for k, v in d.items()))
    print("E4 contract: poiseuille per volume", ch04.ns_terms_preset("poiseuille", 0.0, 2.5e-4, component=0, per="volume",
                                                                    G=100.0, h=1e-3, mu=1e-3))
    u, p = ch04.exact_field("taylor_green", mu=1e-2, rho=1.0)
    X = np.array([[0.3, 1.1], [0.7, 0.2]])
    lap, d2s, cw = ch04.viscous_force_forms(u, X, 0.5, mu=1e-2)
    print(f"N54 (4.40) three forms at 2 points, max spread {max(np.abs(lap - d2s).max(), np.abs(lap - cw).max()):.1e} N/m³")
    comp = lambda x, t: np.stack([x[0] ** 2, x[0] * x[1]])  # noqa: E731  compressible: ∇·u = 3x
    lap, d2s, cw = ch04.viscous_force_forms(comp, np.array([0.5, 0.3]), 0.0, mu=1.0)
    print(f"N54 compressible field: ∇²u = {lap}, 2∂S_ij/∂x_i = {d2s} (differ by ∇(∇·u) = (3, 0))")
    print(f"figures → {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
