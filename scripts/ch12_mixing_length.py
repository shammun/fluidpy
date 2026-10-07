"""§12.10 (C12) with §12.7 (C06): the mixing-length closure — the wall profile from the first integral of (12.100) (closed-form
slope, exact U+ without damping), the additive constant B the model implies with and without van Driest damping, the whole
channel from the linear total stress, the two energy budgets (12.46)-(12.47) on that channel, the laminar parabola at the
same pressure gradient, and the general eddy-viscosity solver for (12.99).

Run: ``.venv/Scripts/python.exe scripts/ch12_mixing_length.py --no-show``
Figures -> outputs/ch12/c12_mixing_length.png, c06_channel_energy_budget.png.
kappa = 0.41 (the cited "classical" preset), A+ = 26 (van Driest's damping constant), core cap 0.09 delta.
"""
from __future__ import annotations

import numpy as np

from ch12_common import COLORS, KAPPA, Timer, finish, parse_args, save, setup

from fluidpy import ch12_turbulence as ch12
from fluidpy.core import wall_turbulence as WT


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    print("additive constant B implied by the model (limit of U+ - ln(y+)/kappa):")
    for k in (0.384, 0.41):
        print(f"  kappa = {k}: no damping {ch12.mixing_length_intercept(k):+.3f} = (ln 4k - 1)/k;  A+ = 26: {ch12.mixing_length_intercept(k, 26.0):.3f}")
    for A in (0.0, 15.0, 26.0, 35.0):
        print(f"  kappa = 0.41, A+ = {A:4.1f}: B = {ch12.mixing_length_intercept(KAPPA, A):.3f}")
    yq = np.array([1.0, 5.0, 12.0, 30.0, 100.0, 1000.0])
    und, dmp = ch12.mixing_length_wall_profile(yq, KAPPA), ch12.mixing_length_wall_profile(yq, KAPPA, "van_driest", 26.0)
    for i, y in enumerate(yq):
        s = und["slope"][i]
        print(f"  y+ = {y:6.0f}: undamped U+ = {und['Uplus'][i]:7.3f} (slope {s:.4f}, quadratic residual {s + (KAPPA * y) ** 2 * s * s - 1:+.1e}), "
              f"damped U+ = {dmp['Uplus'][i]:7.3f}, nu_T/nu = {dmp['nuT_over_nu'][i]:8.3f}, log law(0.41, 5.0) = {WT.log_law(y, kappa=KAPPA, B=5.0):6.3f}")
    grid = np.concatenate([[0.0], np.geomspace(0.01, 1000.0, 4000)])
    errs = []
    for n in (250, 500, 1000, 2000):
        g = np.concatenate([[0.0], np.geomspace(0.01, 1000.0, n)])
        errs.append(abs(ch12.mixing_length_wall_profile(g, KAPPA, "van_driest", method="trapezoid")["Uplus"][-1]
                        - ch12.mixing_length_wall_profile(1000.0, KAPPA, "van_driest")["Uplus"]))
    print(f"trapezoid quadrature error at y+ = 1000 for n = 250..2000: {['%.1e' % e for e in errs]} -> order "
          f"{np.log2(errs[0] / errs[1]):.2f}, {np.log2(errs[1] / errs[2]):.2f}, {np.log2(errs[2] / errs[3]):.2f}")
    print("channel from the mixing length + linear total stress:")
    for Ret in (180.0, 550.0, 1000.0, 5200.0):
        b = ch12.channel_energy_budget(Ret, KAPPA, 26.0, n=800)
        print(f"  Re_tau = {Ret:5.0f}: U_cl+ = {b['U_cl_plus']:.2f}, U_bulk+ = {b['U_bulk_plus']:.2f}, Re_bulk = {b['Re_bulk']:.3e}, Cf = {b['Cf']:.5f}, "
              f"production peak {b['production'].max():.4f} at y+ = {b['yplus_peak_production']:.1f}; work {b['integrals']['work']:.3f} = "
              f"direct dissipation {b['integrals']['dissipation']:.3f} + production {b['integrals']['production']:.3f} "
              f"(residual {b['identity_residual']:+.1e}); direct share {100 * b['direct_fraction']:.1f} %")
    b = ch12.channel_energy_budget(1000.0, KAPPA, 26.0, n=800)
    lam_ub = 1000.0 / 3.0   # laminar channel at the same pressure gradient: U+ = y+ (1 - y+/(2 Re_tau)), bulk = Re_tau/3
    print(f"  laminar flow at the same pressure gradient (Re_tau = 1000): U_bulk+ = {lam_ub:.1f} against the turbulent {b['U_bulk_plus']:.1f}"
          f" - turbulence cuts the flow rate by a factor {lam_ub / b['U_bulk_plus']:.1f}")
    # the same channel with the general solver of (12.99), in wall units on the half channel mirrored to the full one
    with Timer("eddy-viscosity BVP (12.99)"):
        Ret = 1000.0
        yh = b["yplus"]
        yfull = np.concatenate([yh, (2 * Ret - yh[::-1])[1:]])
        nuT = lambda y, dU: (np.minimum(KAPPA * np.minimum(y, 2 * Ret - y), 0.09 * Ret)  # noqa: E731
                             * (-np.expm1(-np.minimum(y, 2 * Ret - y) / 26.0))) ** 2 * np.abs(dU)
        sol = ch12.shear_flow_eddy_viscosity_solve(yfull, nuT, dPdx=-1.0 / Ret, rho=1.0, nu=1.0, tol=1e-11, max_iter=2000)
    print(f"  centreline U+ from the BVP: {sol['U'][yh.size - 1]:.3f} ({sol['iterations']} iterations, converged {sol['converged']}) against "
          f"the quadrature {b['U_cl_plus']:.3f}")
    print(f"convective scales (12.102), our layer L = 2 m, dT = 2 K, T = 290 K: w ~ {ch12.convective_velocity_scale(2.0, 2.0, 290.0):.3f} m/s, "
          f"kappa_T ~ {ch12.convective_eddy_diffusivity(2.0, 2.0, 290.0):.3f} m2/s")

    fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.8))
    yy = np.geomspace(0.1, 1000.0, 200)
    for A, c, lab in ((None, COLORS["rose"], "no damping"), (26.0, COLORS["accent"], r"van Driest, $A^+$ = 26")):
        p = ch12.mixing_length_wall_profile(yy, KAPPA, None if A is None else "van_driest", A or 26.0)
        ax[0].semilogx(yy, p["Uplus"], color=c, label=f"{lab}: B = {p['B']:.2f}")
    ax[0].semilogx(yy[yy > 8], WT.log_law(yy[yy > 8], kappa=KAPPA, B=5.0), "k--", lw=1.0, label="log law, B = 5.0")
    ax[0].semilogx(yy[yy < 15], yy[yy < 15], "k:", label=r"$U^+=y^+$")
    ax[0].set(ylim=(0, 25), xlabel=r"$y^+$", ylabel=r"$U^+$", title=r"$\kappa$ sets the slope, the damping sets the intercept")
    ax[0].legend(fontsize=8)
    eta = b["y_over_delta"]
    ax[1].plot(b["Uplus"] / b["U_cl_plus"], eta, color=COLORS["accent"], label="mixing-length model")
    ax[1].plot(eta * (2 - eta), eta, color=COLORS["muted"], ls="--", label="laminar parabola (shape)")
    ax[1].set(xlabel=r"$U/U_{centre}$", ylabel=r"$y/\delta$", title=r"channel, $Re_\tau$ = 1000: blunter than laminar")
    ax[1].legend(fontsize=8)
    ax[2].semilogx(b["yplus"][1:], b["nuT_over_nu"][1:], color=COLORS["orange"], label=r"$\nu_T/\nu$")
    ax[2].semilogx(b["yplus"][1:], b["lT_plus"][1:], color=COLORS["teal"], label=r"$l_T^+$")
    ax[2].set(xlabel=r"$y^+$", ylabel="wall units", yscale="log", ylim=(1e-3, 1e3), title="eddy viscosity and mixing length")
    ax[2].legend(fontsize=8)
    save(fig, out, "c12_mixing_length")

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    yp = b["yplus"]
    ax[0].semilogx(yp[1:], b["mean"]["pressure_work"][1:], color=COLORS["accent"], label="pressure work")
    ax[0].semilogx(yp[1:], b["mean"]["transport"][1:], color=COLORS["muted"], label="transport")
    ax[0].semilogx(yp[1:], b["mean"]["viscous_dissipation"][1:], color=COLORS["rose"], label="direct viscous dissipation")
    ax[0].semilogx(yp[1:], b["mean"]["loss_to_turbulence"][1:], color=COLORS["orange"], label="loss to turbulence")
    ax[0].set(xlabel=r"$y^+$", ylabel=r"term / $(u_*^4/\nu)$", title="mean-flow budget (12.46)")
    ax[0].legend(fontsize=8)
    ax[1].semilogx(yp[1:], b["turbulence"]["production"][1:], color=COLORS["orange"], label="production (the same term, + sign)")
    ax[1].semilogx(yp[1:], b["turbulence"]["dissipation_plus_transport"][1:], color=COLORS["rose"], ls="--", label="dissipation + transport (residual)")
    ax[1].axvline(b["yplus_peak_production"], color=COLORS["grid"])
    ax[1].set(xlabel=r"$y^+$", ylabel=r"term / $(u_*^4/\nu)$", title=f"turbulence budget (12.47): peak 1/4 near y+ = {b['yplus_peak_production']:.0f}")
    ax[1].legend(fontsize=8)
    save(fig, out, "c06_channel_energy_budget")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
