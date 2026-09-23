"""§3.2: Lagrangian vs Eulerian descriptions (C01, D01) and the material derivative (C02): a space–time diagram of the
stretching map x = X e^{αt} with its Eulerian twin u = αx; a drifting float and a fixed probe in a moving, warming
temperature front with the three terms of (3.5) along the float's path; the worked numbers; convergence of the
stencils and the independent "d/dt of F along a path line" check.

Run: ``.venv/Scripts/python.exe scripts/ch03_material_derivative.py --no-show``
Figures → outputs/ch03/c01_lagrangian_eulerian.png, c02_thermal_front_terms.png.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from fluidpy import ch03_kinematics as ch03  # noqa: E402
from fluidpy.core.style import COLORS  # noqa: E402
from tools.convergence import observed_order  # noqa: E402


def lagrangian_eulerian_figure(alpha: float = 0.5):
    """Left: particle paths x = X e^{αt} for several labels X (the Lagrangian view); right: the Eulerian velocity
    u(x) = αx at two instants (the same line — the field is steady although every particle accelerates)."""
    import matplotlib.pyplot as plt

    t = np.linspace(0, 2, 100)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.8))
    for X in (0.5, 1.0, 1.5, 2.0):
        a1.plot(t, ch03.lagrangian_map_example(X, t, alpha), label=f"X = {X} m")
    a1.set_xlabel("$t$ [s]")
    a1.set_ylabel("$x = X e^{\\alpha t}$ [m]")
    a1.set_title(rf"Lagrangian: particle paths ($\alpha$ = {alpha} s$^{{-1}}$)")
    a1.legend(fontsize=8)
    x = np.linspace(0, 6, 50)
    for tt, ls in ((0.0, "-"), (1.5, "--")):
        a2.plot(x, alpha * x, ls, color=COLORS["teal"], label=f"u(x, t = {tt} s) = αx")
    a2.set_xlabel("$x$ [m]")
    a2.set_ylabel("$u$ [m/s]")
    a2.set_title("Eulerian: the velocity field (steady)")
    a2.legend(fontsize=8)
    return fig


def thermal_front_figure(G=1e-5, v=10.0, c=5.0, H=1.0 / 3600.0, width=2e5, hours=24):
    """A probe fixed at y = 0 and a float carried by the southerly wind v through a tanh front of width ``width``
    moving north at c and heating at H: T(t) at both, and the three terms of (3.5) at the float vs time."""
    import matplotlib.pyplot as plt

    t = np.linspace(0, hours * 3600.0, 200)
    y_float = v * t - 1.0e5  # the float starts 100 km south of the probe
    T_probe = ch03.thermal_front(0.0, 0.0, t, G, H, c, width)
    T_float = ch03.thermal_front(0.0, y_float, t, G, H, c, width)
    terms = ch03.thermal_front_terms(0.0, y_float, t, 0.0, v, G, H, c, width)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 3.8))
    a1.plot(t / 3600, T_probe, color=COLORS["blue"], label="fixed probe (y = 0): ∂T/∂t")
    a1.plot(t / 3600, T_float, color=COLORS["accent"], label="drifting float: DT/Dt")
    a1.set_xlabel("$t$ [h]")
    a1.set_ylabel("$T$ [K]")
    a1.set_title("two observers of one moving, warming front")
    a1.legend(fontsize=8)
    k = 3600.0
    a2.plot(t / 3600, terms["local"] * k, color=COLORS["blue"], label=r"local $\partial T/\partial t$")
    a2.plot(t / 3600, terms["advective"] * k, color=COLORS["amber"], label=r"advective $\mathbf{u}\cdot\nabla T$")
    a2.plot(t / 3600, terms["total"] * k, color=COLORS["accent"], lw=2.4, label="DT/Dt (sum)")
    dTdt_float = np.gradient(T_float, t) * k
    a2.plot(t[::10] / 3600, dTdt_float[::10], "o", color=COLORS["ink"], ms=3, label="float's measured dT/dt")
    a2.set_xlabel("$t$ [h]")
    a2.set_ylabel("rate [K/h]")
    a2.set_title("Eq. (3.5) at the float")
    a2.legend(fontsize=8)
    return fig


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch03"))
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
    lagrangian_eulerian_figure().savefig(out / "c01_lagrangian_eulerian.png", bbox_inches="tight")
    thermal_front_figure().savefig(out / "c02_thermal_front_terms.png", bbox_inches="tight")

    # C01 worked number and D01
    X, a, t = 2.0, 0.5, 1.0
    x, u, acc = (ch03.lagrangian_map_example(X, t, a, k) for k in range(3))
    uf, af = ch03.lagrangian_velocity_acceleration(lambda s: ch03.lagrangian_map_example(X, s, a), t)
    print(f"C01: x = {x:.4f} m, u = {u:.4f} m/s (= αx {a * x:.4f}), a = {acc:.4f} m/s²; stencils {uf:.6f}, {af:.6f}")
    Xs, ts, al, xs = sp.symbols("X t alpha x", real=True)
    d = ch03.lagrangian_to_eulerian([Xs * sp.exp(al * ts)], [Xs], [xs], ts)
    print(f"D01: label X = {d['label_of_x'][Xs]}, u(x, t) = {d['u'][0]}, Du/Dt − a = {sp.simplify(d['Du_Dt'][0] - d['a_lagrangian'][0])}")

    # C02 worked number: warm advection
    w = ch03.thermal_front_terms(0.0, 0.0, 0.0, 0.0, 10.0, 1e-5, 0.0, 10.0)
    print(f"C02: v ∂T/∂y = {w['advective']:+.1e} K/s, ∂T/∂t = {w['local']:+.1e} K/s ({w['local'] * 3600:+.2f} K/h at a station), "
          f"DT/Dt = {w['total']:+.1e} K/s → {w['regime']}")
    print("N11 index form: DF/Dt − ∂F/∂t =", ch03.expand_indices_str("u_i F_,i"))

    # independent check: DF/Dt from the stencils vs d/dt of F sampled along a path line (C02 from-scratch)
    u_field = lambda X_, T_: np.stack([0.3 * np.sin(X_[1]) + 0.1 * T_, 0.5 * np.cos(X_[0])])  # noqa: E731
    F = lambda X_, T_: np.sin(X_[0]) * np.cos(2 * X_[1]) + 0.2 * T_ ** 2  # noqa: E731
    r0, tp = np.array([0.3, 0.4]), 0.7
    ts_ = tp + np.array([-1e-3, 0.0, 1e-3])
    path = ch03.pathline(u_field, r0, tp, ts_)
    dFdt_path = (F(path[:, 2], ts_[2]) - F(path[:, 0], ts_[0])) / 2e-3
    DF = ch03.material_derivative(F, u_field, r0, tp)
    print(f"DF/Dt: stencil {DF:.8f} vs d/dt F(path line) {dFdt_path:.8f} (diff {abs(DF - dFdt_path):.1e})")
    Fs = sp.sin(xs) * sp.cos(2 * sp.Symbol("y")) + sp.Rational(1, 5) * ts ** 2
    hs = np.array([0.1, 0.05, 0.025, 0.0125])
    exact_expr = ch03.material_derivative_sym(Fs, [0.3 * sp.sin(sp.Symbol("y")) + 0.1 * ts, 0.5 * sp.cos(xs)],
                                              [xs, sp.Symbol("y")], ts)
    exact = float(exact_expr.subs({xs: 0.3, sp.Symbol("y"): 0.4, ts: tp}))
    errs = [abs(ch03.material_derivative(F, u_field, r0, tp, h=h) - exact) for h in hs]
    print(f"stencil vs sympy: {DF - exact:+.1e}; observed order {observed_order(hs, errs):.3f}")
    sw = ch03.streamwise_derivative(F, u_field, r0, tp)
    print(f"(3.6) |u| ∂F/∂s = {sw:.8f} vs advective u·∇F = {ch03.material_derivative_terms(F, u_field, r0, tp).advective:.8f}")
    print(f"figures → {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
