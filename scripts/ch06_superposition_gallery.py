"""§6.2–6.3 (C02–C04): the element kit as flow nets — uniform stream, vortex, source, doublet — with the Laplace
residual of each element, ω = −∇²ψ (6.4) on a Rankine core, the δ-source fluxes (6.6)/(6.13) on circles of every
radius, and the source–sink pair (6.28) shrinking to the doublet (6.29) with its ε² error.

Run: ``.venv/Scripts/python.exe scripts/ch06_superposition_gallery.py --no-show``
Figures → outputs/ch06/c02_element_nets.png, c04_doublet_limit.png.
"""
from __future__ import annotations

import numpy as np

from ch06_drawings import parse_args, save, setup, flow_net

from fluidpy import ch06_ideal_flow as ch06
from fluidpy.core.style import COLORS


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    n = 101 if args.fast else 161
    elems = [("uniform U = 1, V = 0.5 m/s  (6.7), (6.14)", ch06.Flow([ch06.Uniform(1.0, 0.5)])),
             ("vortex Γ = 2π m²/s (ccw)  (6.8)", ch06.Flow([ch06.Vortex(2 * np.pi, cut_angle=0.0)])),
             ("source m = 2π m²/s  (6.15)", ch06.Flow([ch06.Source(2 * np.pi, cut_angle=0.0)])),
             ("doublet d = −1 e_x m³/s  (6.29)", ch06.Flow([ch06.Doublet.from_book_scalar(1.0)]))]
    fig, ax = plt.subplots(1, 4, figsize=(17, 4.4))
    print("Laplace residual ∇²ψ, ∇²φ (4th-order stencil, h = 1e-3) at (0.8, 0.6):")
    for a, (lab, fl) in zip(ax, elems):
        flow_net(a, fl, (-2, 2), (-2, 2), n=n, phi=True, title=lab, cut_ray=fl.elements[0].kind in ("vortex", "source"))
        rp = ch06.laplacian_residual(fl.psi, 0.8, 0.6)
        rf = ch06.laplacian_residual(fl.phi, 0.8, 0.6)
        print(f"  {lab.split('  ')[0]:<34s} ∇²ψ = {rp: .2e}  ∇²φ = {rf: .2e}")
    save(fig, out, "c02_element_nets")

    # ω = −∇²ψ on a Rankine core ψ (u_θ = Γr/2πσ² inside): ω = Γ/πσ² inside, 0 outside
    G, sig = 2 * np.pi, 0.5

    def psi_rankine(x, y):
        return ch06.rankine_vortex_psi(x, y, Gamma=G, a=sig)

    print(f"ω = −∇²ψ, Rankine core Γ = 2π, σ = 0.5: inside {ch06.vorticity_from_psi(psi_rankine, 0.1, 0.2):.6f}"
          f" (exact {G / (np.pi * sig ** 2):.6f}); outside {ch06.vorticity_from_psi(psi_rankine, 1.0, 0.7):.2e}")
    fl_v = ch06.Flow([ch06.Vortex(2.0)])
    fl_s = ch06.Flow([ch06.Source(3.0)])
    radii = [0.01, 1.0, 100.0]
    print("flux of ∇ψ (vortex Γ = 2, expect −Γ):", np.round(ch06.delta_flux_check(fl_v, radii=radii), 12))
    print("flux of ∇φ (source m = 3, expect m):  ", np.round(ch06.delta_flux_check(fl_s, radii=radii, kind="phi"), 12))
    print("flux of ∇φ on a circle not enclosing the source:",
          f"{ch06.delta_flux_check(fl_s, center=(3.0, 0.0), radii=[1.0], kind='phi')[0]:.2e}")

    # the doublet limit (6.28)–(6.29)
    eps = np.array([0.2, 0.1, 0.05, 0.025])
    err = ch06.doublet_limit_error(eps, d=1.0)
    order = np.polyfit(np.log(eps), np.log(err), 1)[0]
    print("source–sink pair vs doublet, relative error:", ", ".join(f"ε = {e:g}: {v:.3e}" for e, v in zip(eps, err)),
          f"→ observed order {order:.3f}")
    fr = ch06.doublet_limit_frames((0.5, 0.2, 0.05), d=1.0, n=n)
    X0, Y0 = fr[0]["X"], fr[0]["Y"]
    psi_d = np.asarray(ch06.Doublet.from_book_scalar(1.0).psi(X0, Y0))
    fig2, ax2 = plt.subplots(1, 4, figsize=(17, 4.2))
    lv = np.linspace(-0.6, 0.6, 25)
    for a, frame in zip(ax2[:3], fr):
        e, P = frame["eps"], frame["psi"]
        a.contour(X0, Y0, P, levels=lv, colors=COLORS["ink"], linewidths=0.8, negative_linestyles="solid")
        a.plot([-e, e], [0, 0], "o", color=COLORS["teal"])
        a.plot([e], [0], "o", color=COLORS["rose"])
        a.set_aspect("equal")
        a.set_title(f"pair, m = d/2ε, ε = {e:g} m", fontsize=10)
        a.set_xlabel("x [m]")
    ax2[3].contour(X0, Y0, psi_d, levels=lv, colors=COLORS["accent"], linewidths=0.8, negative_linestyles="solid")
    ax2[3].set_aspect("equal")
    ax2[3].set_title("doublet (6.29): circles tangent to x at 0", fontsize=10)
    save(fig2, out, "c04_doublet_limit")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
