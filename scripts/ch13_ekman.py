"""Chapter 13, §13.6–§13.7 — Ekman layers.

* the surface spiral in both hemispheres: 45° at the surface, transport at 90° to the wind and independent of the eddy
  viscosity (checked by quadrature for three viscosities and for three K(z) profiles solved numerically);
* the bottom layer: cross-isobar flow toward low pressure, the three-force balance, transport ½Uδ, pumping;
* an Ekman layer in water of finite depth (ours), and the convergence order of the K(z) solver.

Run: ``.venv/Scripts/python.exe scripts/ch13_ekman.py --no-show``   Figures -> outputs/ch13/ekman_*.png
"""
from __future__ import annotations

import numpy as np
from ch13_common import COLORS, finish, hemispheres, parse_args, save, setup

from fluidpy import ch13_geophysical_fluid_dynamics as ch13


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch13.illustrative_inputs()
    tau, rho, nu = I["tau"], I["rho_ocean"], I["nu_v_ocean"]
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    for lat, label in hemispheres(I["lat_ekman"]):
        f = ch13.coriolis_parameter(lat)
        d = ch13.ekman_depth(nu, f)
        z = np.linspace(-14.0 * d, 0.0, 4001)
        u, v = ch13.ekman_surface(z, tau, 0.0, rho, nu, f)
        u0, v0 = ch13.ekman_surface(0.0, tau, 0.0, rho, nu, f)
        M = ch13.ekman_transport(tau, 0.0, rho, f)
        rx, ry = ch13.ekman_residual(z, u, v, nu, f)
        print(f"surface layer at {label}: delta = {d:.2f} m (pi*delta = {ch13.ekman_depth(nu, f, 'pi'):.1f} m), surface current "
              f"{np.hypot(u0, v0) * 100:.2f} cm/s at {np.rad2deg(np.arctan2(v0, u0)):+.1f}° from the wind, transport "
              f"({M[0]:+.3f}, {M[1]:+.3f}) m²/s; quadrature ({np.trapezoid(u, z):+.1e}, {np.trapezoid(v, z):+.4f}); "
              f"max residual {max(np.max(np.abs(rx)), np.max(np.abs(ry))) / (abs(f) * np.hypot(u0, v0)):.1e} of f|V|")
        ax[0].plot(u * 100, v * 100, label=label)
    ax[0].annotate("", (2.0, 0.0), (0.0, 0.0), arrowprops=dict(arrowstyle="->", color=COLORS["orange"], lw=2))
    ax[0].text(2.0, 0.1, "wind stress", color=COLORS["orange"], fontsize=8)
    ax[0].set(aspect="equal", xlabel="u [cm/s]", ylabel="v [cm/s]", title="surface Ekman spiral (hodograph)")
    ax[0].legend(fontsize=8)

    f = ch13.coriolis_parameter(I["lat_ekman"])
    print("transport does not depend on the eddy viscosity (60° N, stress along x):")
    for fac in (0.25, 1.0, 4.0):
        d = ch13.ekman_depth(fac * nu, f)
        z = np.linspace(-16.0 * d, 0.0, 6001)
        u, v = ch13.ekman_surface(z, tau, 0.0, rho, fac * nu, f)
        print(f"  nu_v = {fac * nu:.4f} m²/s: delta = {d:5.2f} m, surface speed {np.hypot(u[-1], v[-1]) * 100:5.2f} cm/s, "
              f"transport by quadrature {np.trapezoid(v, z):+.5f} m²/s (closed form {ch13.ekman_transport(tau, 0.0, rho, f)[1]:+.5f})")
    d = ch13.ekman_depth(nu, f)
    profiles = {"constant K": lambda q: nu + 0.0 * q, "K growing with depth": lambda q: nu * (0.2 - q / d),
                "K decaying with depth": lambda q: nu * np.exp(q / (2.0 * d)) + 1.0e-4}
    zs = -20.0 * d * np.linspace(1.0, 0.0, 1601) ** 2             # clustered near the surface
    for name, K in profiles.items():
        V = ch13.ekman_solve(zs, K, f, tau=tau + 0.0j, rho=rho, bottom="stressfree")
        Mq = np.trapezoid(V, zs)
        print(f"  K(z) solver, {name:22s}: transport ({Mq.real:+.2e}, {Mq.imag:+.6f}) m²/s, surface angle "
              f"{np.rad2deg(np.angle(V[-1])):+.1f}°")
        ax[1].plot(V.real * 100, V.imag * 100, label=name)
    ax[1].set(aspect="equal", xlabel="u [cm/s]", ylabel="v [cm/s]", title="same stress, three K(z): same transport")
    ax[1].legend(fontsize=7)
    errs = []
    for n in (101, 201, 401, 801):
        zz = np.linspace(-14.0 * d, 0.0, n)
        V = ch13.ekman_solve(zz, nu, f, tau=tau + 0.0j, rho=rho)
        errs.append(np.max(np.abs(V - ch13.ekman_surface(zz, tau, 0.0, rho, nu, f, as_complex=True))))
    print(f"  ekman_solve against the closed form: error {errs[-2] / abs(ch13.ekman_surface(0.0, tau, 0.0, rho, nu, f, as_complex=True)):.1e} "
          f"of the surface speed at 401 nodes; observed order {np.round(np.log2(np.array(errs[:-1]) / np.array(errs[1:])), 2)}")

    U, nua = I["U_g"], I["nu_v_atm"]
    for lat, label in hemispheres(I["lat_ekman"]):
        fa = ch13.coriolis_parameter(lat)
        da = ch13.ekman_depth(nua, fa)
        z = np.linspace(0.0, 16.0 * da, 6001)
        ub, vb = ch13.ekman_bottom(z, U, 0.0, nua, fa)
        Mb = ch13.ekman_bottom_transport(U, 0.0, nua, fa)
        fb0, fb = ch13.ekman_force_balance(0.0, U, nua, fa), ch13.ekman_force_balance(0.5 * da, U, nua, fa)
        zeta = 2.0e-5 * np.sign(fa)                               # a cyclone in this hemisphere
        print(f"bottom layer at {label}: delta = {da:.0f} m, cross-isobar transport {Mb[1]:+.0f} m²/s (quadrature "
              f"{np.trapezoid(vb, z):+.0f}; ½Uδ = {0.5 * U * da:.0f}), angle to the isobars {np.rad2deg(fb0['angle_to_isobars']):.0f}° at the "
              f"ground and {np.rad2deg(fb['angle_to_isobars']):.1f}° at z = δ/2, force sum {fb['sum']}, pumping under a cyclone "
              f"(|zeta| = 2e-5 1/s): w = {ch13.ekman_pumping_bottom(zeta, nua, fa) * 100:+.2f} cm/s; u exceeds U by at most "
              f"{(ub.max() / U - 1) * 100:.2f} % (at z = 3 pi delta/4) and by {np.exp(-np.pi) * 100:.2f} % at z = pi delta, where v first vanishes")
        ax[2].plot(ub, vb, label=label)
    ax[2].plot([U], [0], "o", color=COLORS["muted"])
    ax[2].set(aspect="equal", xlabel="u [m/s]", ylabel="v [m/s]", title="Ekman layer on a rigid surface (hodograph)")
    ax[2].legend(fontsize=8)
    save(fig, out, "ekman_spirals")

    H = 4.0 * np.pi * ch13.ekman_depth(nu, f)
    zf = np.linspace(-H, 0.0, 801)
    Vf = ch13.ekman_finite_depth(zf, tau, 0.0, 0.05, 0.0, H, nu, rho, f)
    fig2, ax2 = plt.subplots(1, 2, figsize=(10, 4.2))
    ax2[0].plot(Vf.real * 100, zf, color=COLORS["accent"], label="u")
    ax2[0].plot(Vf.imag * 100, zf, color=COLORS["teal"], label="v")
    ax2[0].set(xlabel="velocity [cm/s]", ylabel="z [m]", title="finite depth: surface layer over bottom layer (ours)")
    ax2[0].legend()
    ax2[1].plot(Vf.real * 100, Vf.imag * 100, color=COLORS["accent"])
    ax2[1].set(aspect="equal", xlabel="u [cm/s]", ylabel="v [cm/s]", title="hodograph")
    save(fig2, out, "ekman_finite_depth")
    print(f"finite depth H = {H:.0f} m with an interior flow of 5 cm/s: V(0) = {Vf[-1] * 100:.2f} cm/s, V(-H) = {abs(Vf[0]):.1e}")
    xg = np.linspace(-1.0e6, 1.0e6, 81)
    Xg, Yg = np.meshgrid(xg, xg)
    tx = -tau * np.cos(np.pi * Yg / 2.0e6)                        # a zonal wind that changes sign across the box (a gyre-like curl)
    wE = ch13.ekman_pumping(tx, 0.0 * tx, xg[1] - xg[0], xg[1] - xg[0], rho, ch13.beta_plane(xg, I["lat"])[:, None])
    curl = -tau * np.pi / 2.0e6 * np.sin(np.pi * 0.25)            # d(-tau_x)/dy at y = L/4
    print(f"Ekman pumping of a sinusoidal wind (ours): w_E between {wE.min() * 86400 * 365.25:+.0f} and {wE.max() * 86400 * 365.25:+.0f} m/yr; "
          f"Sverdrup transport for curl tau = {curl:.2e} N/m³: V = {ch13.sverdrup_transport(curl, rho, ch13.beta_parameter(I['lat'])):+.2f} m²/s")
    print(f"eddy viscosity implied by an Ekman thickness of 600 m at 60°: {ch13.eddy_viscosity_from_depth(600.0, f):.1f} m²/s")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
