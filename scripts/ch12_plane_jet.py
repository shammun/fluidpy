"""§12.8 (C09): the plane turbulent jet — raw profiles at several stations and their collapse on F(y/x) (12.56), (12.66), the
invariant momentum flux (12.62), the growing volume flux (12.68) and constant scalar flux (12.70), the cross-flow and
entrainment velocity (12.58), the stress profile from (12.63), the thin-layer residual (12.61), the similarity coefficients
(12.64), (12.74), and a model energy budget (12.75; qualitative).  The laminar jet of ch09 is drawn as a ghost.

Run: ``.venv/Scripts/python.exe scripts/ch12_plane_jet.py --no-show [--fast]``
Figures -> outputs/ch12/c09_plane_jet_similarity.png, c09_plane_jet_budget.png.
Our jet: air slot d = 10 mm, U0 = 15 m/s; half-widths 0.10 (velocity) and 0.15 (scalar) are ILLUSTRATIVE, C5 from the invariant.
"""
from __future__ import annotations

import numpy as np

from ch12_common import COLORS, ILLUSTRATIVE_JET, finish, parse_args, save, setup

from fluidpy import ch12_turbulence as ch12
from fluidpy.core import jets


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    rho, nu, d, U0 = 1.2, 1.5e-5, 0.010, 15.0
    Js, Ms = ch12.slot_momentum_flux(rho, U0, d), ch12.slot_mass_flux(rho, U0, d)
    xh, xhY = ILLUSTRATIVE_JET["xi_half_U"], ILLUSTRATIVE_JET["xi_half_Y"]
    kw = dict(C5="from_invariant", xi_half=xh)
    pi_ = ch12.profile_integrals(xh, xhY)
    print(f"slot: J_s = {Js:.3f} N/m, M_s = {Ms:.4f} kg/(m s); profile integrals {pi_}; C5 from the invariant = {1 / np.sqrt(pi_['I2']):.4f}")
    y = np.linspace(-3.0, 3.0, 12001 if args.fast else 60001)
    xs = (0.25, 0.5, 1.0, 2.0)
    print("station   U_CL [m/s]  half-width [m]  J [N/m]   V_dot [m2/s]  scalar flux [kg/(m s)]  v_e [m/s]  Re_local")
    for x in xs:
        U = ch12.plane_jet_mean_velocity(x, y, Js, rho, **kw)
        Y = ch12.plane_jet_mass_fraction(x, y, Ms, Js, rho, C6="from_invariant", xi_half_Y=xhY, **kw)
        Ucl = ch12.plane_jet_centerline_velocity(x, Js, rho, **kw)
        print(f"  x = {x:4.2f}  {Ucl:8.3f}   {xh * x:8.4f}      {ch12.jet_momentum_flux_per_span(y, U, rho):.4f}   "
              f"{np.trapezoid(U, y):.4f}        {ch12.scalar_flux_per_span(y, U, Y, rho):.5f}            "
              f"{ch12.plane_jet_entrainment_velocity(x, Js, rho, **kw):.4f}    {Ucl * xh * x / nu:.2e}")
    print(f"volume flux formula (12.68) at x = 2: {ch12.plane_jet_volume_flux(2.0, Js, rho, **kw):.4f} m2/s; "
          f"V(y -> inf) = {ch12.plane_jet_cross_velocity(2.0, 100.0, Js, rho, **kw):+.4f} m/s")
    ex = ch12.free_shear_exponents("plane_jet", return_equations=True)
    print(f"exponents: width x^{ex['width']}, U_CL x^{ex['velocity']}, scalar x^{ex['scalar']}, Reynolds number x^{ex['reynolds']}; from {ex['equations']}")
    yy = np.linspace(-0.6, 0.6, 121)
    tl = ch12.thin_shear_layer_terms(1.0, yy, Js, rho, nu, **kw)
    print(f"thin-layer balance (12.61) at x = 1 m: max residual {np.max(np.abs(tl['residual'])):.1e} against stress gradient "
          f"{np.max(np.abs(tl['stress_gradient'])):.2f} m/s2; viscous/turbulent {tl['viscous_over_stress']:.1e}; max V/U {tl['ratio_V_over_U']:.3f}")
    g1 = ch12.general_similarity_check(lambda x: x, lambda x: x ** -0.5, lambda x: 1.0 / x, 1.0)
    g2 = ch12.general_similarity_check(lambda x: np.exp(0.5 * x), lambda x: np.exp(-0.5 * x), lambda x: np.exp(-0.5 * x), 1.0)
    print(f"(12.64) power family m = 1, n = -1/2: c1 = {g1['c1']:.3f}, c2 = {g1['c2']:.3f}, c3 = {g1['c3']:.3f}, d ln J/d ln x = {g1['momentum_flux_exponent']:.1e}")
    print(f"(12.74) exponential family: c1 = {g2['c1']:.3f}, c2 = {g2['c2']:.1e} (vanishes), c3 = {g2['c3']:.3f}, d ln J/d ln x = {g2['momentum_flux_exponent']:.3f}")
    mom, vol = ch12.wrong_exponent_fluxes(4.0, -0.4, 1.0)
    print(f"wrong decay exponent n = -0.4: momentum flux x {mom:.3f} and volume flux x {vol:.3f} from x = 1 to 4 (momentum exponent 2n + m = {2 * -0.4 + 1:.1f}) - not allowed")
    ev = ch12.plane_jet_eddy_viscosity_profile(yy, xh)
    print(f"uniform eddy viscosity closure: F = sech^2({ev['a']:.3f} xi), nu_T/(U_CL x) = {ev['nu_hat']:.5f}")
    lam1, lam4 = jets.free_jet(1.0, 0.0, Js, rho, nu), jets.free_jet(4.0, 0.0, Js, rho, nu)
    print(f"laminar jet of the same J (ch09) if it could stay laminar: u0(1 m) = {lam1['u0']:.2f} m/s, delta = {lam1['delta'] * 1e3:.2f} mm; "
          f"from x = 1 to 4 m u0 falls by {lam1['u0'] / lam4['u0']:.3f} = 4^(1/3), the turbulent U_CL by "
          f"{ch12.plane_jet_centerline_velocity(1.0, Js, rho, **kw) / ch12.plane_jet_centerline_velocity(4.0, Js, rho, **kw):.3f} = 4^(1/2)")
    s = ch12.sympy_summary("plane_jet_similarity")
    print(f"sympy (12.63): {s['checks']}")

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    cols = (COLORS["muted"], COLORS["teal"], COLORS["orange"], COLORS["accent"])
    for x, c in zip(xs, cols):
        U = ch12.plane_jet_mean_velocity(x, yy, Js, rho, **kw)
        ax[0].plot(yy, U, color=c, label=f"x = {x} m")
        ax[1].plot(yy / x, U / ch12.plane_jet_centerline_velocity(x, Js, rho, **kw), color=c, lw=3 if x == xs[0] else 1.2)
    ax[0].set(xlabel="y [m]", ylabel="U [m/s]", title="raw profiles: slower and wider downstream")
    ax[0].legend(fontsize=8)
    ax[1].plot(yy / 0.25, ch12.plane_jet_eddy_viscosity_profile(yy / 0.25, xh)["F"], "k:", lw=1.0, label=r"sech$^2$ (uniform $\nu_T$; the laminar shape)")
    ax[1].set(xlim=(-0.4, 0.4), xlabel=r"$\xi=y/x$", ylabel=r"$U/U_{CL}$", title="rescaled: one curve F (12.56)")
    ax[1].legend(fontsize=8)
    xx = np.linspace(0.2, 3.0, 40)
    ax[2].plot(xx, [ch12.jet_momentum_flux_per_span(y, ch12.plane_jet_mean_velocity(x, y, Js, rho, **kw), rho) / Js for x in xx],
               color=COLORS["accent"], label=r"momentum flux / $J_s$ (12.62)")
    ax[2].plot(xx, ch12.plane_jet_volume_flux(xx, Js, rho, **kw) / (U0 * d), color=COLORS["teal"], label=r"volume flux / $U_0 d$ (12.68)")
    ax[2].plot(xx, [ch12.scalar_flux_per_span(y, ch12.plane_jet_mean_velocity(x, y, Js, rho, **kw),
                                              ch12.plane_jet_mass_fraction(x, y, Ms, Js, rho, C6="from_invariant", xi_half_Y=xhY, **kw), rho) / Ms for x in xx],
               color=COLORS["amber"], ls="--", label=r"scalar flux / $\dot M_s$ (12.70)")
    ax[2].set(xlabel="x [m]", ylabel="flux / source value", title="what is conserved and what is entrained")
    ax[2].legend(fontsize=8)
    save(fig, out, "c09_plane_jet_similarity")

    xi = np.linspace(0.0, 0.3, 121)
    b = ch12.jet_tke_budget(xi, xh)
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.4))
    ax[0].plot(xi, ch12.gaussian_profile(xi, xh), color=COLORS["accent"], label="F (mean velocity)")
    ax[0].plot(xi, -10 * ch12.plane_jet_stress_profile(xi, xi_half=xh), color=COLORS["orange"], label=r"$10\times\overline{uv}/U_{CL}^2$ from (12.63)")
    ax[0].plot(xi, 10 * b["e"], color=COLORS["teal"], ls="--", label=r"$10\times\bar e/U_{CL}^2$ (assumed shape)")
    ax[0].set(xlabel=r"$\xi=y/x$", ylabel="profile", title="stress is zero on the axis, largest near the steepest shear")
    ax[0].legend(fontsize=8)
    for key, c in (("production", COLORS["orange"]), ("advection", COLORS["accent"]), ("transport", COLORS["teal"]), ("dissipation", COLORS["rose"])):
        ax[1].plot(xi, b[key], color=c, label=key)
    ax[1].axhline(0, color=COLORS["grid"])
    ax[1].set(xlabel=r"$\xi=y/x$", ylabel=r"term / $(U_{CL}^3/x)$", title="energy budget (12.75): a MODEL, qualitative")
    ax[1].legend(fontsize=8)
    save(fig, out, "c09_plane_jet_budget")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
