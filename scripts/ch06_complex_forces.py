"""§6.4–6.5 (C09, C10): the complex potential — corner flows w = Azⁿ (6.46) from a 45° wedge to a flat plate with
|dw/dz| ∝ r^{n−1}, Cauchy–Riemann residuals (6.44) for every element and for the non-analytic z* — and forces:
Blasius (6.60) independent of the contour (circles and polygons), the Laurent bars showing only the 1/z term of
(dw/dz)² survives, the sympy Kutta–Zhukhovsky residue (6.61)–(6.62) with the book's printed 1/z² coefficient beside
the correct one, and held singularities (Exercise 6.10).

Run: ``.venv/Scripts/python.exe scripts/ch06_complex_forces.py --no-show``
Figures → outputs/ch06/c09_corners.png, c10_blasius.png.
"""
from __future__ import annotations

import numpy as np

from ch06_drawings import flow_net, parse_args, save, setup

from fluidpy import ch06_ideal_flow as ch06
from fluidpy.core.style import COLORS


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    print("Cauchy–Riemann residuals (6.44) at z = 0.7 + 0.3i (h = 1e-6):")
    for name, p in (("corner", dict(n=2.5)), ("source", {}), ("vortex", {}), ("doublet", {}), ("cylinder", dict(Gamma_cw=1.0)),
                    ("conj", {})):
        r1, r2 = ch06.cauchy_riemann_residual(name, 0.7 + 0.3j, **p)
        pr = ch06.complex_potential_probe(name, 0.7, 0.3, **({"n": 2.5} if name == "corner" else p))
        print(f"  {name:<9s} r1 = {r1: .2e}, r2 = {r2: .2e}  (u, v) = ({pr['u']:.5f}, {pr['v']:.5f})")
    ns = (4.0, 2.0, 2.0 / 3.0, 0.5)
    fig, ax = plt.subplots(1, 5, figsize=(21, 4.2))
    n = 121 if args.fast else 181
    for axi, nn in zip(ax[:4], ns):
        c = ch06.Corner(1.0, nn)
        info = ch06.corner_info(nn)
        flow_net(axi, ch06.Flow([c]), (-2, 2), (-2, 2), n=n, phi=True, title=f"w = zⁿ, n = {nn:.3g}: α = "
                 f"{info['alpha_deg']:.0f}°, {info['regime']}")
        axi.plot([0, 2.2], [0, 0], color=COLORS["rose"], lw=2)
        axi.plot([0, 2.2 * np.cos(info["alpha"])], [0, 2.2 * np.sin(info["alpha"])], color=COLORS["rose"], lw=2)
        r = np.logspace(-4, -1, 7)
        sp_ = np.abs(c.dwdz(r * np.exp(0.5j * info["alpha"])))
        slope = np.polyfit(np.log(r), np.log(sp_), 1)[0]
        zw = r * np.exp(1j * (info["alpha"] - 1e-12))  # just inside the second wall
        wall_psi = max(np.max(np.abs(c.psi(r, 0 * r))), np.max(np.abs(c.psi(zw.real, zw.imag))))
        print(f"  corner n = {nn:.4g}: α = {info['alpha_deg']:.1f}°, log–log slope of |dw/dz| {slope:.6f} "
              f"(n − 1 = {info['exponent']:.6f}); max |ψ| on both walls {wall_psi:.1e}")
        ax[4].loglog(r, sp_, "o-", label=f"n = {nn:.3g}")
    ax[4].set_xlabel("r [m]")
    ax[4].set_ylabel("|dw/dz| [m/s]")
    ax[4].set_title("speed near the corner ∝ r^(n−1)", fontsize=10)
    ax[4].legend(fontsize=8)
    save(fig, out, "c09_corners")

    kz = ch06.kutta_zhukhovsky_sym()
    print(f"sympy (6.61)–(6.62): residue {kz['residue']}, D = {kz['D']}, L = {kz['L']}; 1/z² coefficient correct "
          f"{kz['coeff_z2']} vs printed {kz['coeff_z2_printed']}")
    U, a, G, rho = 1.0, 1.0, 2.0, 1.0
    fl = ch06.cylinder(U, a, Gamma_cw=G)
    Rs = np.linspace(1.05, 5.0, 12)
    Ls = [ch06.blasius_force(fl, R=R, rho=rho).L for R in Rs]
    print(f"cylinder Γ_cw = {G}: Blasius L on 12 circles R ∈ [1.05, 5] — spread {np.ptp(Ls):.1e}, value {Ls[0]:.12f}")
    ns_ = [8, 16, 32, 64]
    errs = [abs(ch06.blasius_force(fl, R=1.2, rho=rho, n=k).L - rho * U * G) for k in ns_]
    print("  circle trapezoid error vs nodes:", ", ".join(f"n = {k}: {e:.1e}" for k, e in zip(ns_, errs)))
    polyerr = []
    for k in (32, 64, 128, 256):
        t = np.linspace(0, 2 * np.pi, k, endpoint=False)
        poly = 2.0 * np.cos(t) + 1j * 1.4 * np.sin(t) + 0.3
        polyerr.append(abs(ch06.blasius_force(fl, contour=poly, rho=rho).L - rho * U * G))
    print("  shifted-ellipse polygon contour error, N = 32…256:", [f"{e:.2e}" for e in polyerr],
          f"→ order {np.polyfit(np.log([32, 64, 128, 256]), np.log(polyerr), 1)[0]:.3f}")
    ef = ch06.elliptic_cylinder_flow(U, 1.2, 1.0, Gamma_cw=G)
    print(f"Zhukhovsky ellipse, same Γ: Blasius {ch06.blasius_force(ef, R=3.0, rho=rho)} (L = ρUΓ = {rho * U * G})")
    lc = ch06.laurent_contributions(fl, 2.0, rho=rho)
    print("Laurent coefficients of dw/dz:", {k: np.round(v, 10) for k, v in lc["coeffs"].items() if abs(v) > 1e-12})
    print("  contributions to ∮(dw/dz)²dz from pairs with j + k = −1:", [(j, k, np.round(v, 10)) for j, k, v in lc["pairs"]])
    print("held source m = 3 in U = 2 (ρ = 1.5):", ch06.force_on_held_singularity("source", 2.0, 3.0, 1.5),
          " held ccw vortex Γ = 3:", ch06.force_on_held_singularity("vortex", 2.0, 3.0, 1.5))

    fig2, ax2 = plt.subplots(1, 2, figsize=(12, 4.2))
    ax2[0].plot(Rs, Ls, "o-", color=COLORS["amber"], label="L by Blasius on |z| = R")
    ax2[0].plot(Rs, [ch06.cv_force_on_body(fl, R, rho=rho).L for R in Rs], "s", color=COLORS["teal"], ms=4,
                label="L by the control volume (6.54)")
    ax2[0].axhline(rho * U * G, color=COLORS["muted"], ls="--", label="ρUΓ (6.62)")
    ax2[0].set_xlabel("contour radius R [m]")
    ax2[0].set_ylabel("L [N/m]")
    ax2[0].set_ylim(0, 3)
    ax2[0].legend(fontsize=8)
    ax2[0].set_title("lift does not depend on the contour", fontsize=10)
    pw = sorted(lc["square"])
    ax2[1].bar([str(p) for p in pw], [abs(lc["square"][p]) for p in pw], color=COLORS["muted"])
    ax2[1].bar(["-1"], [abs(lc["square"][-1])], color=COLORS["amber"], label="only z⁻¹ has a residue")
    ax2[1].set_xlabel("power p of z in (dw/dz)²")
    ax2[1].set_ylabel("|coefficient|")
    ax2[1].legend(fontsize=8)
    ax2[1].set_title("Laurent terms of (dw/dz)² (cylinder, Γ_cw = 2)", fontsize=10)
    save(fig2, out, "c10_blasius")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
