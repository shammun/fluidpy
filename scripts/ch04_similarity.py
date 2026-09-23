"""§4.11 dimensionless forms and dynamic similarity (C15, E9): the coefficients of the scaled Navier–Stokes and energy
equations (sympy), synthetic sphere-drag "experiments" that scatter in dimensional axes and collapse on C_D(Re) (Fig.
4.21 drawn with the Morrison correlation), Example 4.8's ship-model extrapolation, Prandtl numbers and model design.

Run: ``.venv/Scripts/python.exe scripts/ch04_similarity.py --no-show``
Figures → outputs/ch04/c15_sphere_drag_collapse.png, c15_ship_drag.png, c15_prandtl.png.
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


def sphere_figure():
    """Left: F_D vs U for many (d, fluid) — no pattern; right: C_D vs Re — one curve (Morrison + Stokes line)."""
    import matplotlib.pyplot as plt

    d = ch04.synthetic_sphere_drag_data(n=60, seed=0)
    Re = np.logspace(-1, 6, 400)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
    cols = {"air": COLORS["blue"], "water": COLORS["teal"], "glycerine": COLORS["orange"]}
    for f in cols:
        m = d["fluid"] == f
        a1.loglog(d["U"][m], d["F"][m], "o", color=cols[f], ms=4, label=f)
        a2.loglog(d["Re"][m], d["CD"][m], "o", color=cols[f], ms=4)
    a1.set_xlabel("U [m/s]")
    a1.set_ylabel("F_D [N]")
    a1.set_title("dimensional: 5 parameters, no pattern")
    a1.legend(fontsize=8)
    a2.loglog(Re, ch04.sphere_drag_coefficient(Re), color=COLORS["ink"], label="Morrison (2013) correlation")
    a2.loglog(Re[Re < 10], ch04.sphere_drag_coefficient(Re[Re < 10], "stokes"), "--", color=COLORS["rose"], label="24/Re")
    a2.set_xlabel("Re = ρUd/μ")
    a2.set_ylabel(r"$C_D = F_D/(\frac{1}{2}\rho U^2 A)$")
    a2.set_title("dimensionless (4.99): one curve (Fig. 4.21 analogue)")
    a2.legend(fontsize=8)
    return fig


def ship_figure(**kw):
    """Example 4.8: model total → friction + wave → scaled wave + prototype friction; the uncorrected estimate."""
    import matplotlib.pyplot as plt

    r = ch04.ship_drag_extrapolation(**kw)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4))
    a1.bar(["friction", "wave"], [r["D_m_friction"], r["D_m_wave"]], color=[COLORS["teal"], COLORS["blue"]])
    a1.set_ylabel("model drag [N]")
    a1.set_title(f"model at U_m = {r['U_m']:.2f} m/s")
    a2.bar(["wave (scaled)", "friction", "total", "uncorrected"],
           np.array([r["D_p_wave"], r["D_p_friction"], r["D_p_total"], r["D_p_uncorrected"]]) / 1e5,
           color=[COLORS["blue"], COLORS["teal"], COLORS["accent"], COLORS["muted"]])
    a2.set_ylabel("prototype drag [10⁵ N]")
    a2.set_title(f"Re_p/Re_m = {r['Re_ratio']:.0f}: correct the friction")
    return fig, r


def prandtl_figure():
    """Pr(T) of air and water from ch01's properties; Eucken's kinetic-theory values."""
    import matplotlib.pyplot as plt

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.8))
    Ta = np.linspace(250, 350, 21)
    a1.plot(Ta, [ch04.prandtl_of("air", T) for T in Ta], color=COLORS["blue"], label="air")
    a1.axhline(float(ch04.eucken_prandtl(1.4)), color=COLORS["orange"], ls="--", label="Eucken γ = 1.4")
    a1.axhline(float(ch04.eucken_prandtl(5 / 3)), color=COLORS["rose"], ls=":", label="Eucken γ = 5/3 (2/3)")
    a1.set_xlabel("T [K]")
    a1.set_ylabel("Pr")
    a1.legend(fontsize=8)
    Tw = np.linspace(283.15, 303.15, 21)
    a2.plot(Tw, [ch04.prandtl_of("water", T) for T in Tw], color=COLORS["teal"])
    a2.set_xlabel("T [K]")
    a2.set_ylabel("Pr (water)")
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
    sphere_figure().savefig(out / "c15_sphere_drag_collapse.png", bbox_inches="tight")
    fig, r = ship_figure(rho_p=1000.0)
    fig.savefig(out / "c15_ship_drag.png", bbox_inches="tight")
    prandtl_figure().savefig(out / "c15_prandtl.png", bbox_inches="tight")

    print("C15 (4.101) coefficients:", ch04.nondimensional_ns_coefficients())
    print("C15 viscous pressure scaling:", ch04.nondimensional_ns_coefficients("viscous"))
    print("C15 (4.114) coefficients:", ch04.nondimensional_energy_coefficients(), "; (4.110):",
          ch04.nondimensional_continuity_coefficient())
    s = ch04.Scales.from_oscillation(0.1, 2.0, 1000.0, 1e-3)
    print(f"C15 oscillating body l = 0.1 m, Ω = 2 rad/s: groups {s.groups()}")
    print(f"R12 (4.99) Π groups: {ch04.sphere_drag_pi_groups()}")
    print(f"Ex. 4.8 (ρ_p = ρ_m = 1000): U_m {r['U_m']:.2f} m/s, model friction {r['D_m_friction']:.3f} N, wave {r['D_m_wave']:.2f} N, "
          f"prototype wave {r['D_p_wave']:.4e} N + friction {r['D_p_friction']:.4e} N = {r['D_p_total']:.4e} N "
          f"(uncorrected {r['D_p_uncorrected']:.4e} N); Re_p/Re_m = {r['Re_ratio']:.1f}")
    print("E9 model design 1:25 ship in water:", ch04.model_prototype(100.0, 10.0, 1 / 25))
    print(f"N152 Pr: air 300 K {ch04.prandtl_of('air', 300.0):.4f}, water 293.15 K {ch04.prandtl_of('water'):.3f}, "
          f"Eucken monatomic {ch04.eucken_prandtl(5 / 3):.4f}")
    print(f"N147 M(100 m/s, 288.15 K) = {ch04.mach_number(100.0, T=288.15):.4f}; incompressible limit "
          f"{ch04.incompressible_speed_limit():.1f} m/s; U²/c² at M = 0.3: {ch04.compressibility_parameter(0.3, 1.0):.2f}")
    print(f"figures → {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
