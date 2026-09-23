"""§5.4 and §5.6 the vorticity equation (C06, C09, C10, N19–N22, N30–N38): term bars of (5.13)/(5.30) for the E5/E7
scenes, the sympy derivation checks ((5.12) curls, (B.3.10), (5.30) → (5.13)), the diffusing vortex sheet, Burgers'
stretching–diffusion balance, stretching and tilting of a vortex line in linear flows, a from-scratch (ω·∇)u and
ν∇²ω vs ``vorticity_terms``, ∇·ω = 0 (5.18), the Lamb form (5.25) against (5.20), and the helix frame.

Run: ``.venv/Scripts/python.exe scripts/ch05_vorticity_equation.py --no-show``
Figures → outputs/ch05/c06_budget_bars.png, c06_sheet_diffusion.png, c10_burgers_balance.png, c10_stretch_tilt.png.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from fluidpy import ch05_vorticity_dynamics as ch05  # noqa: E402
from fluidpy.core.style import COLORS  # noqa: E402
from tools.convergence import observed_order  # noqa: E402

TERM_C = dict(local=COLORS["muted"], advective=COLORS["ink"], stretching_tilting=COLORS["accent"],
              planetary=COLORS["amber"], baroclinic=COLORS["orange"], diffusion=COLORS["rose"],
              residual=COLORS["grid"])


def budget_figure():
    import matplotlib.pyplot as plt

    names = ch05.VORTICITY_BUDGET_PRESETS
    fig, axs = plt.subplots(1, len(names), figsize=(18, 3.6))
    for ax, nm in zip(axs, names):
        d = ch05.vorticity_budget_preset(nm, component=1 if nm == "hill" else 2)
        keys = list(d)
        ax.bar(range(len(keys)), [d[k] for k in keys], color=[TERM_C[k] for k in keys])
        ax.set_xticks(range(len(keys)))
        ax.set_xticklabels([k.replace("_", "\n") for k in keys], fontsize=6)
        ax.set_title(nm, fontsize=9)
    axs[0].set_ylabel("term of (5.30) [1/s²]")
    return fig


def sheet_figure():
    import matplotlib.pyplot as plt

    y = np.linspace(-0.02, 0.02, 401)
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
    for t, c in zip((10.0, 100.0, 1000.0), (COLORS["accent"], COLORS["teal"], COLORS["orange"])):
        u, w = ch05.diffusing_vortex_sheet(y, t, 0.1, 1e-6)
        ax[0].plot(y * 1e3, w, color=c, label=f"t = {t:g} s, ∫ω dy = {np.trapezoid(w, y):.4f} m/s")
        ax[1].plot(y * 1e3, u, color=c)
    ax[0].set_xlabel("y [mm]")
    ax[0].set_ylabel("ω_z [1/s]")
    ax[0].legend(fontsize=7)
    ax[1].set_xlabel("y [mm]")
    ax[1].set_ylabel("u [m/s]")
    ax[1].set_title("u = −(γ/2) erf(y/2√(νt)): the jump spreads")
    return fig


def burgers_figure():
    import matplotlib.pyplot as plt

    R = np.linspace(0, 6e-3, 300)
    b = ch05.burgers_balance(R)
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(R * 1e3, b["advective"], color=COLORS["ink"], label="advective u_R∂ω_z/∂R")
    ax.plot(R * 1e3, b["stretching"], color=COLORS["accent"], label="stretching αω_z")
    ax.plot(R * 1e3, b["diffusion"], color=COLORS["rose"], label="diffusion ν∇²ω_z")
    ax.axvline(b["core_radius"] * 1e3, ls="--", color=COLORS["muted"], label="√(4ν/α)")
    ax.set_xlabel("R [mm]")
    ax.set_ylabel("[1/s²]")
    ax.legend(fontsize=8)
    ax.set_title("Burgers' vortex: stretching balanced by diffusion")
    return fig


def stretch_figure():
    import matplotlib.pyplot as plt

    t = np.linspace(0.0, 3.0, 61)
    fig, ax = plt.subplots(figsize=(6.5, 4))
    for pre, ang, c in (("axial_stretch", 0.0, COLORS["accent"]), ("shear_tilt", np.pi / 2, COLORS["blue"]),
                        ("planar", 0.0, COLORS["muted"])):
        w = ch05.uniform_strain_vorticity([np.sin(ang), 0.0, np.cos(ang)], pre, t)
        ax.plot(t, np.linalg.norm(w, axis=1), color=c, label=f"{pre}: |ω|/|ω₀|")
    ax.set_yscale("log")
    ax.set_xlabel("αt or st")
    ax.set_ylabel("|ω(t)|/|ω₀|")
    ax.legend(fontsize=8)
    ax.set_title("stretching grows ω exponentially, tilting linearly, 2-D not at all")
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

    from fluidpy.core.style import use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    for nm in ch05.VORTICITY_BUDGET_PRESETS:
        d = ch05.vorticity_budget_preset(nm, component=1 if nm == "hill" else 2)
        print(f"(5.30) {nm:15s}: " + ", ".join(f"{k} {v:+.4g}" for k, v in d.items()))
    x, y, z, t = sp.symbols("x y z t", real=True)
    nu = sp.Symbol("nu", positive=True)
    tg = [sp.cos(x) * sp.sin(y) * sp.exp(-2 * nu * t), -sp.sin(x) * sp.cos(y) * sp.exp(-2 * nu * t), 0]
    ve = ch05.vorticity_equation_sym(tg, (x, y, z), t, nu, p_expr=sp.Function("p")(x, y, z), Phi_expr=9.81 * z)
    print(f"sympy (5.12): curl(−∇p/ρ) = {list(ve['curl_pressure'])}, curl(−∇Φ) = {list(ve['curl_gravity'])}, "
          f"(B.3.10) residual {list(ve['identity_B310'])}, curl NS − (5.13) = {list(ve['residual_513'])}")
    vb = ch05.vorticity_budget_sym(tg, (x, y, z), t, nu)
    print(f"sympy (5.30): residual {list(vb['residual_530'])}, reduces to (5.13): {vb['reduces_to_513']}")
    # from scratch: (ω·∇)u and ν∇²ω on Burgers' vortex by central differences
    Gam, al, nuB = 1e-3, 1.0, 1e-6
    u = ch05.burgers_vortex_field(Gam, al, nuB)
    p0, h = np.array([1e-3, 0.4e-3, 0.2e-3]), 1e-6
    E = np.eye(3) * h
    Gm = np.stack([(u(p0 + E[j]) - u(p0 - E[j])) / (2 * h) for j in range(3)], axis=1)
    w = lambda q: ch05.vorticity_field(u, h)(q)  # noqa: E731
    lap = sum((w(p0 + E[j]) - 2 * w(p0) + w(p0 - E[j])) / h ** 2 for j in range(3))
    vt = ch05.vorticity_terms(u, p0, 0.0, nuB, h=h)
    print(f"from scratch vs library (Burgers): (ω·∇)u_z {(Gm @ w(p0))[2]:.6f} vs {vt.stretching_tilting[2]:.6f}, "
          f"ν∇²ω_z {nuB * lap[2]:.6f} vs {vt.diffusion[2]:.6f}; residual {vt.residual[2]:+.2e} 1/s²")
    hs = np.array([8e-5, 4e-5, 2e-5, 1e-5])  # h/R from 0.08 to 0.01 (smaller h: round-off of nested stencils)
    err = [abs(ch05.vorticity_terms(u, p0, 0.0, nuB, h=hh).residual[2]) for hh in hs]
    print(f"stencil order of the (5.13) residual on Burgers: {observed_order(hs, err):.2f}")
    print(f"(5.18) ∇·ω of the ABC flow by nested stencils: {ch05.vorticity_divergence(ch05.abc_flow_field(), [0.1, 0.2, 0.3]):.1e}")
    uu = lambda X, tt=0.0: np.stack([np.sin(X[1]) * np.cos(X[2]), np.cos(X[0]) * np.sin(X[2]),  # noqa: E731
                                     0.3 * np.sin(X[0] + X[1])])
    pp = lambda X, tt=0.0: 1000 * np.cos(X[0]) * np.sin(X[1] + X[2])  # noqa: E731
    lf = ch05.rotating_lamb_form_terms(uu, pp, 1000.0, lambda X, tt=0.0: 9.81 * X[2], [0.3, 0.4, 0.5], 0.0, 0.01,
                                       (0, 0, 0.3))
    ns = ch05.rotating_ns_residual(uu, pp, 1000.0, [0.3, 0.4, 0.5], 0.0, 0.01, (0, 0, 0.3), (0, 0, -9.81))
    print(f"(5.25) residual {np.round(lf['residual'], 8)} = (5.20) residual {np.round(ns, 8)} (∇·u = 0 test field)")
    b = ch05.burgers_balance(0.0)
    print(f"Burgers (α = 1, ν = 1e-6, Γ = 1e-3): core {b['core_radius'] * 1e3:.3f} mm, peak ω_z {b['omega_z']:.4f} 1/s")
    for pre, ang in (("axial_stretch", 0.0), ("shear_tilt", np.pi / 2), ("planar", 0.0)):
        sc = ch05.stretching_tilting_scenario(pre, 2.0, angle=ang)
        print(f"{pre:13s} t = 2: |ω|/|ω₀| = {sc['growth']:.6f}, stretch rate {sc['stretch_rate']:+.3f}, tilt rate "
              f"{sc['tilt_rate']:.3f} 1/s — {sc['status']}")
    st = ch05.stretched_tube(2.0, 1.0, 10.0, 1e-4)
    print(f"inviscid tube stretched ×2: ω {st['omega']:.1f} 1/s, A {st['A']:.1e} m², Γ {st['Gamma']:.1e} m²/s")
    hf = ch05.helix_frame(0.0)
    print(f"helix frame (a = 1, c = 0.3): κ = {hf['curvature']:.6f}, τ = {hf['torsion']:.6f}, e_n = {hf['e_n']}")
    u_s, w_s = ch05.diffusing_vortex_sheet(1e-3, 1.0, 1.0, 1e-6)
    print(f"diffusing sheet γ = 1 m/s, ν = 1e-6, t = 1 s: ω(0) = {ch05.diffusing_vortex_sheet(0.0, 1.0, 1.0, 1e-6)[1]:.3f}"
          f" 1/s, u(1 mm) = {u_s:.6f} m/s")
    budget_figure().savefig(out / "c06_budget_bars.png", bbox_inches="tight")
    sheet_figure().savefig(out / "c06_sheet_diffusion.png", bbox_inches="tight")
    burgers_figure().savefig(out / "c10_burgers_balance.png", bbox_inches="tight")
    stretch_figure().savefig(out / "c10_stretch_tilt.png", bbox_inches="tight")
    print(f"figures → {out}")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
