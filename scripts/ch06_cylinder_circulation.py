"""§6.3 (C07) and §6.5 (C10): the cylinder with clockwise circulation (6.36)–(6.40) — stagnation points moving down
and leaving the body at Γ = 4πaU (Fig. 6.12), lift L = ρUΓ by three independent routes (surface pressure (6.39),
Blasius (6.60) on several contours, control volume (6.54)), the sign convention pinned (Γ_ccw = +… gives L < 0), and
the non-uniqueness family (every Γ satisfies the same boundary conditions).

Run: ``.venv/Scripts/python.exe scripts/ch06_cylinder_circulation.py --no-show``
Figures → outputs/ch06/c07_cylinder_circulation.png.
"""
from __future__ import annotations

import numpy as np

from ch06_drawings import circle, flow_net, parse_args, save, setup

from fluidpy import ch06_ideal_flow as ch06


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    U, a, rho = 10.0, 0.1, 1.225
    G = 2.0
    st = ch06.cylinder_circulation_state(U, a, Gamma_cw=G, rho=rho)
    print(f"U = {U} m/s, a = {a} m, Γ_cw = {G} m²/s (air ρ = {rho}): Γ/4πaU = {st['ratio']:.4f}, {st['regime']}")
    print(f"  stagnation angles {np.round(st['stagnation_angles_deg'], 4)} deg; sin θ = −Γ/4πaU = {-st['ratio']:.6f}")
    fl = ch06.cylinder(U, a, Gamma_cw=G)
    print(f"  L: pressure integral {st['L']:.10f}, ρUΓ {st['rhoUGamma']:.10f} N/m; D = {st['D']:.2e} N/m")
    for R in (0.15, 0.5, 2.0):
        print(f"  Blasius on R = {R} m: {ch06.blasius_force(fl, R=R, rho=rho)};  CV (6.54): "
              f"{ch06.cv_force_on_body(fl, R, rho=rho)}")
    print(f"  sign pin: Gamma_ccw = +2 gives L = {ch06.lift_per_span(rho, U, Gamma_ccw=2.0):+.3f} N/m "
          f"(Blasius {ch06.blasius_force(ch06.cylinder(U, a, Gamma_ccw=2.0), R=0.3, rho=rho).L:+.3f})")
    fam = ch06.circulation_family_check(1.0, 1.0, [0.0, 2 * np.pi, 4 * np.pi, 8 * np.pi])
    print("  non-uniqueness family Γ_cw = 0, 2π, 4π, 8π (U = a = 1): max u·n on r = a",
          [f"{v:.1e}" for v in fam["max_normal"]], "; far field at r = 1000a:", [f"{v:.1e}" for v in fam["far_field"]],
          "; loop circulation (r = 1.5, 3, 10):", [[round(c, 10) for c in cc] for cc in fam["circulation_by_radius"]])
    sp = ch06.cylinder_stagnation_points(1.0, 1.0, Gamma_cw=6 * np.pi, include_inside=True)
    print(f"  Γ = 6πaU: stagnation points {np.round(sp, 6)}, r₊r₋ = {abs(sp[0] * sp[1]):.12f} = a²")

    ratios = [0.0, 0.5, 1.0, 1.5]
    fig, ax = plt.subplots(1, 4, figsize=(18, 4.4))
    n = 121 if args.fast else 181
    for axi, q in zip(ax, ratios):
        Gq = q * 4 * np.pi
        f = ch06.cylinder(1.0, 1.0, Gamma_cw=Gq)
        stq = ch06.cylinder_stagnation_points(1.0, 1.0, Gamma_cw=Gq)
        flow_net(axi, f, (-3, 3), (-3, 2.5), n=n, psi_levels=np.linspace(-3, 3, 25), body=circle(1.0), stag=stq,
                 title=f"Γ_cw/4πaU = {q:g}: L = ρUΓ = {ch06.lift_per_span(1.0, 1.0, Gamma_cw=Gq):.2f}ρ")
    save(fig, out, "c07_cylinder_circulation")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
