"""§12.9 (C10, C11): wall-bounded turbulence — wall units (12.80)-(12.81), the sublayer (12.82), the logarithmic law (12.88)
and its defect form (12.89), Spalding's inner profile and Coles' wake (Fig. 12.18 analogue), the linear total stress of a
channel and its viscous / Reynolds partition (Fig. 12.17 analogue), the indicator function (12.87), the zero-pressure-gradient
fits with our own plate (Example 12.3 method), skin-friction laws against Blasius, the correlation (12.92), rough walls
(12.93) and the pipe friction law.

Run: ``.venv/Scripts/python.exe scripts/ch12_wall_layers.py --no-show``
Figures -> outputs/ch12/c10_law_of_the_wall.png, c11_friction_laws.png.
Log-law constants: the cited public preset WT.LOG_LAW_CONSTANTS["classical"] (kappa = 0.41, B = 5.0); wake strength 0.3 and
the roughness lengths are illustrative.  Our plate: x = 3 m, U_inf = 40 m/s, air.
"""
from __future__ import annotations

import numpy as np

from ch12_common import B_LOG, COLORS, KAPPA, finish, parse_args, save, setup

from fluidpy.core import wall_turbulence as WT
from fluidpy.core.boundary_layer import blasius_skin_friction
from fluidpy.core.laminar import pipe_friction_factor


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    C = dict(kappa=KAPPA, B=B_LOG)
    tau0, rho, nu = 0.9, 1.2, 1.5e-5
    yp, Up, us, lnu = WT.wall_units(1e-3, 8.0, tau0, rho, nu)
    print(f"wall units for tau0 = {tau0} Pa in air: u* = {us:.4f} m/s, l_nu = {lnu * 1e6:.1f} um; y = 1 mm -> y+ = {yp:.1f}, U = 8 m/s -> U+ = {Up:.2f}")
    print("Pi groups of (12.79):", [{k: str(v) for k, v in g.items()} for g in WT.law_of_the_wall_groups()])
    print("Pi groups of (12.83):", [{k: str(v) for k, v in g.items()} for g in WT.defect_law_groups()])
    print(f"sublayer line meets the log law at y+ = {WT.log_law_crossing(**C):.3f}")
    for y in (1.0, 5.0, 12.0, 30.0, 100.0, 1000.0):
        print(f"  y+ = {y:6.0f}: sublayer {WT.viscous_sublayer(y):8.2f}  log law {WT.log_law(y, **C):6.2f}  Spalding {WT.spalding_uplus(y, **C):6.3f}"
              f"  -> {WT.layer_name(y, y / 5200.0)}")
    ypl = np.geomspace(0.1, 5200.0, 400)
    prof = WT.composite_profile_plus(ypl, 5200.0, Pi=0.3, **C)
    fit0 = WT.fit_log_law(ypl[ypl > 100], WT.log_law(ypl[ypl > 100], **C), window=(100.0, 0.15), Re_tau=5200.0)
    fit = WT.fit_log_law(ypl, prof, window=(100.0, 0.15), Re_tau=5200.0)
    print(f"straight-line fit in 100 < y+ < 0.15 Re_tau: pure log law -> kappa = {fit0['kappa']:.4f}, B = {fit0['B']:.3f}; composite profile "
          f"(Re_tau = 5200, Pi = 0.3) -> kappa = {fit['kappa']:.4f}, B = {fit['B']:.3f} (the wake already bends the profile inside the window)")
    ind = WT.log_law_indicator(ypl, prof)
    print(f"indicator y+ dU+/dy+ at y+ = 300: {np.interp(300.0, ypl, ind):.4f} (1/kappa = {1 / KAPPA:.4f})")
    sp = WT.stress_partition(np.array([0.0, 5.0, 11.0, 30.0, 100.0]), 1000.0, **C)
    print(f"stress partition at Re_tau = 1000, y+ = 0, 5, 11, 30, 100: viscous {np.round(sp['viscous'], 3)}, Reynolds {np.round(sp['reynolds'], 3)}")
    h = 0.05
    print(f"channel h = {h} m, tau0 = {tau0} Pa: dP/dx = {WT.channel_pressure_gradient(tau0, h):.1f} Pa/m (12.90); total stress at y = h/4, h/2: "
          f"{WT.channel_total_stress(h / 4, h, tau0):.3f}, {WT.channel_total_stress(h / 2, h, tau0):.3f} Pa; pipe d = {h} m: "
          f"{WT.pipe_pressure_gradient(tau0, h):.1f} Pa/m (12.91)")
    z = WT.zpg_boundary_layer(3.0, 40.0, nu, kappa=KAPPA)
    print(f"flat plate x = 3 m, U = 40 m/s (Re_x = {z['Re_x']:.2e}): theta = {z['theta'] * 1e3:.2f} mm, delta* = {z['delta_star'] * 1e3:.2f} mm, "
          f"delta99 = {z['delta99'] * 1e3:.1f} mm, H = {z['H']:.3f}, Cf = {z['Cf']:.5f}, u* = {z['u_star']:.3f} m/s, Re_tau = {z['Re_tau']:.0f}")
    print(f"  laminar (Blasius) at the same Re_x would give delta99 = {5.0 * 3.0 / np.sqrt(z['Re_x']) * 1e3:.2f} mm, Cf = {blasius_skin_friction(z['Re_x']):.5f}")
    for Re in (1e6, 1e7, 1e8, 1e9):
        vals = {law: WT.skin_friction_zpg(Re, law, kappa=KAPPA) for law in ("monkewitz", "schultz_grunow", "white")}
        print(f"  Re_x = {Re:.0e}: " + ", ".join(f"{k} {v:.5f}" for k, v in vals.items()) + f", Blasius {blasius_skin_friction(Re):.5f}")
    print(f"correlation (12.92): kappa(B = 5.0) = {WT.nagib_chauhan_kappa(5.0):.4f}; B(kappa = 0.41) = {WT.nagib_chauhan_B(0.41):.3f}; "
          f"B(kappa = 0.384) = {WT.nagib_chauhan_B(0.384):.3f}")
    print(f"rough wall (12.93), u* = 0.3 m/s: U(10 m) over z0 = 0.2 mm (sea) {WT.rough_wall_log_law(10.0, 0.3, 2e-4, kappa=0.40):.2f} m/s, "
          f"over z0 = 3 cm (grass) {WT.rough_wall_log_law(10.0, 0.3, 0.03, kappa=0.40):.2f} m/s; C_D(10 m) = "
          f"{WT.drag_coefficient_neutral(10.0, 2e-4, kappa=0.40):.2e}, {WT.drag_coefficient_neutral(10.0, 0.03, kappa=0.40):.2e}")
    for Re in (1e4, 1e5, 1e6):
        print(f"  pipe Re_d = {Re:.0e}: f = {WT.pipe_friction_factor_turbulent(Re):.5f} (Prandtl's law), "
              f"{WT.pipe_friction_factor_turbulent(Re, **C):.5f} (log law integrated, kappa = 0.41, B = 5.0), laminar 64/Re = {pipe_friction_factor(Re):.5f}")
    slope = np.log(10.0) / (KAPPA * np.sqrt(8.0))
    off = (B_LOG - 1.5 / KAPPA - np.log(2.0 * np.sqrt(8.0)) / KAPPA) / np.sqrt(8.0)
    print(f"  friction-law constants derived from (0.41, 5.0): 1/sqrt(f) = {slope:.4f} log10(Re sqrt f) {off:+.4f}   (Prandtl: 2.0, -0.8)")

    fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.8))
    for Ret, c in ((180.0, COLORS["muted"]), (1000.0, COLORS["teal"]), (5200.0, COLORS["accent"])):
        yq = np.geomspace(0.1, Ret, 300)
        ax[0].semilogx(yq, WT.composite_profile_plus(yq, Ret, Pi=0.3, **C), color=c, label=rf"composite, $Re_\tau$ = {Ret:.0f}")
    yq = np.geomspace(0.1, 5200.0, 300)
    ax[0].semilogx(yq[yq < 20], WT.viscous_sublayer(yq[yq < 20]), "k:", label=r"$U^+=y^+$ (12.82)")
    ax[0].semilogx(yq[yq > 5], WT.log_law(yq[yq > 5], **C), "k--", lw=1.0, label=r"log law (12.88), $\kappa$ = 0.41, B = 5.0")
    for xb in (5.0, 30.0):
        ax[0].axvline(xb, color=COLORS["grid"])
    ax[0].set(ylim=(0, 30), xlabel=r"$y^+$", ylabel=r"$U^+$", title="law of the wall: inner collapse, outer departure")
    ax[0].legend(fontsize=7)
    yq = np.linspace(0.0, 1000.0, 2001)
    sp = WT.stress_partition(yq, 1000.0, **C)
    ax[1].plot(sp["total"], yq / 1000.0, "k", label=r"total $\bar\tau/\tau_0 = 1 - y/\delta$")
    ax[1].plot(sp["viscous"], yq / 1000.0, color=COLORS["rose"], label="viscous")
    ax[1].plot(sp["reynolds"], yq / 1000.0, color=COLORS["orange"], label=r"Reynolds $-\rho\overline{uv}$")
    ax[1].set(xlabel="stress / wall stress", ylabel=r"$y/\delta$", title=r"channel, $Re_\tau$ = 1000: who carries the stress")
    ax[1].legend(fontsize=8)
    ax[2].semilogx(ypl, ind, color=COLORS["accent"])
    ax[2].axhline(1 / KAPPA, color=COLORS["muted"], ls="--", label=r"$1/\kappa$")
    ax[2].set(ylim=(0, 6), xlabel=r"$y^+$", ylabel=r"$y^+\,dU^+/dy^+$", title=r"indicator (12.87), $Re_\tau$ = 5200")
    ax[2].legend(fontsize=8)
    save(fig, out, "c10_law_of_the_wall")

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    Re = np.geomspace(1e6, 1e9, 60)
    for law, c in (("monkewitz", COLORS["accent"]), ("schultz_grunow", COLORS["teal"]), ("white", COLORS["orange"])):
        ax[0].loglog(Re, WT.skin_friction_zpg(Re, law, kappa=KAPPA), color=c, label=law.replace("_", "-"))
    Rel = np.geomspace(1e4, 1e9, 60)
    ax[0].loglog(Rel, blasius_skin_friction(Rel), color=COLORS["muted"], ls="--", label="laminar (Blasius)")
    ax[0].set(xlabel=r"$Re_x$", ylabel=r"$C_f$", title="flat-plate skin friction")
    ax[0].legend(fontsize=8)
    Red = np.geomspace(4e3, 1e7, 50)
    ax[1].loglog(Red, WT.pipe_friction_factor_turbulent(Red), color=COLORS["accent"], label="turbulent (Prandtl's law)")
    ax[1].loglog(Red, WT.pipe_friction_factor_turbulent(Red, **C), color=COLORS["teal"], ls=":", label="log law integrated (0.41, 5.0)")
    Rl = np.geomspace(5e2, 4e3, 20)
    ax[1].loglog(Rl, pipe_friction_factor(Rl), color=COLORS["muted"], ls="--", label="laminar 64/Re")
    ax[1].set(xlabel=r"$Re_d$", ylabel="Darcy friction factor f", title="smooth pipe")
    ax[1].legend(fontsize=8)
    save(fig, out, "c11_friction_laws")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
