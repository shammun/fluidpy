"""Chapter 13, §13.5 — geostrophic balance, how a parcel gets there, the thermal wind and the Taylor–Proudman limit.

* a Gaussian low in both hemispheres: wind along the isobars, the sense of rotation flips with the sign of f;
* a parcel released from rest in a uniform pressure gradient: inertial circles round the geostrophic velocity, and a
  spiral into a cross-isobar flow when drag is added (closed form, ours);
* a front and its jet in thermal-wind balance (westerly in both hemispheres);
* grid convergence of the gridded geostrophic wind (order 2).

Run: ``.venv/Scripts/python.exe scripts/ch13_geostrophic.py --no-show``   Figures -> outputs/ch13/geostrophic_*.png
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
    dp, Rc, rho0 = 1500.0, 5.0e5, 1.2
    x = np.linspace(-2.0e6, 2.0e6, 81)
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    for a, (lat, label) in zip(ax[:2], hemispheres(I["lat"])):
        pc = ch13.pressure_centre(x, x, dp=dp, R_c=Rc, lat_rad=lat, rho0=rho0, kind="low")
        X, Y = np.meshgrid(x, x)
        circ = float(np.sum(X * pc["v"] - Y * pc["u"]))           # > 0 counter-clockwise
        sp = np.hypot(pc["u"], pc["v"])
        ug, vg = ch13.geostrophic_from_field(pc["p"], x[1] - x[0], x[1] - x[0], pc["f"], rho0)
        print(f"low at {label}: f = {pc['f']:+.3e} 1/s, max wind {sp.max():.1f} m/s at r = R_c, rotation "
              f"{'counter-clockwise' if circ > 0 else 'clockwise'}; Ro = {ch13.rossby_number(sp.max(), pc['f'], Rc):.2f}; "
              f"max |u.grad p| = {np.max(np.abs(pc['u'] * (-pc['p'] * X / Rc ** 2) + pc['v'] * (-pc['p'] * Y / Rc ** 2))):.1e}; "
              f"gridded vs analytic wind: {np.max(np.hypot(ug - pc['u'], vg - pc['v'])) / sp.max():.2e}")
        a.contour(X / 1e3, Y / 1e3, pc["p"], 8, colors=COLORS["orange"], linewidths=0.8)
        s = (slice(None, None, 5), slice(None, None, 5))
        a.quiver(X[s] / 1e3, Y[s] / 1e3, pc["u"][s], pc["v"][s], color=COLORS["teal"])
        a.set(aspect="equal", xlabel="x [km]", ylabel="y [km]", title=f"low at {label}: wind along isobars")
    errs = []
    for n in (21, 41, 81, 161):
        xs = np.linspace(-2.0e6, 2.0e6, n)
        pc = ch13.pressure_centre(xs, xs, dp=dp, R_c=Rc, lat_rad=I["lat"], rho0=rho0)
        ug, vg = ch13.geostrophic_from_field(pc["p"], xs[1] - xs[0], xs[1] - xs[0], pc["f"], rho0)
        errs.append(np.max(np.hypot(ug - pc["u"], vg - pc["v"])))
    print("geostrophic_from_field: observed order", np.round(np.log2(np.array(errs[:-1]) / np.array(errs[1:])), 2))

    f = ch13.coriolis_parameter(I["lat"])
    G = 1.0e-3 + 0.0j                                         # pressure-gradient acceleration, pointing east
    t = np.linspace(0.0, 3.0 * ch13.inertial_period(f), 600)
    for r, col, lab in ((0.0, COLORS["accent"], "no drag"), (0.3 * abs(f), COLORS["rose"], "drag r = 0.3 |f|")):
        V = ch13.parcel_adjust(t, G, f, r)
        ax[2].plot(V.real, V.imag, color=col, label=lab)
        Vinf = -G / (r + 1j * f)
        ax[2].plot([Vinf.real], [Vinf.imag], "o", color=col)
        print(f"parcel, {lab}: end state u = {Vinf.real:+.2f}, v = {Vinf.imag:+.2f} m/s; angle to the isobars "
              f"{np.rad2deg(np.arctan2(-Vinf.real, abs(Vinf.imag))):.1f}° toward low pressure")
    ug, vg = ch13.geostrophic_velocity(G.real * rho0, 0.0, f, rho0)
    print(f"  geostrophic velocity for the same gradient: ({ug:+.2f}, {vg:+.2f}) m/s; time-mean of the frictionless parcel "
          f"over 3 inertial periods: ({np.mean(ch13.parcel_adjust(t, G, f).real):+.2f}, {np.mean(ch13.parcel_adjust(t, G, f).imag):+.2f})")
    ax[2].set(aspect="equal", xlabel="u [m/s]", ylabel="v [m/s]", title="a parcel released from rest (pressure falls toward −x)")
    ax[2].legend(fontsize=8)
    save(fig, out, "geostrophic_balance")

    y, z = np.linspace(-2.0e6, 2.0e6, 121), np.linspace(0.0, 1.0e4, 41)
    fig2, ax2 = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
    for a, (lat, label) in zip(ax2, hemispheres(I["lat"])):
        js = ch13.jet_section(y, z, dT=14.0, width=6.0e5, lat_rad=lat, alpha=1.0 / 280.0)
        sh = ch13.thermal_wind_from_temperature(0.0, js["dTdy"][60], ch13.coriolis_parameter(lat), alpha=1.0 / 280.0)
        print(f"jet at {label}: dT/dy at the front = {js['dTdy'][60] * 1e6:+.2f} K per 1000 km, shear du/dz = {sh[0] * 1e3:+.3f} m/s per km, "
              f"wind at 10 km = {js['U'][-1, 60]:+.1f} m/s ({ch13.wind_from_to(js['U'][-1, 60], 0.0)['wind_name']})")
        cs = a.contourf(y / 1e3, z / 1e3, js["U"], 12, cmap="viridis")
        a.contour(y / 1e3, z / 1e3, js["T"], 10, colors=COLORS["blue"], linewidths=0.7)
        a.set(xlabel="y [km] (north →)", title=f"{label}: isotherms (blue) and zonal wind")
        fig2.colorbar(cs, ax=a, label="U [m/s]")
    ax2[0].set_ylabel("z [km]")
    save(fig2, out, "geostrophic_thermal_wind")
    tp = ch13.taylor_proudman_residual(lambda xx, yy, zz: (-yy, xx, 0.0), (1.0, 2.0, 3.0))
    tp2 = ch13.taylor_proudman_residual(lambda xx, yy, zz: (-yy * zz, xx, 0.0), (1.0, 2.0, 3.0))
    print(f"Taylor–Proudman residual: a z-independent vortex {tp}; a sheared one {tp2}")
    print(f"Ekman number of the same low with nu_v = {I['nu_v_atm']} m²/s over a depth of 1 km: "
          f"{ch13.ekman_number(I['nu_v_atm'], f, 1.0e3):.3f}")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
