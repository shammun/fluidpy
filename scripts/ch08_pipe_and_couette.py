"""§8.2 (C03, C04): Hagen–Poiseuille pipe flow (8.6)–(8.8) with Q ∝ a⁴ and f = 64/Re, and circular Couette flow
(8.10) with its limits (8.11) (rotating cylinder in a bath, with the free-surface dip of Fig. 8.7 from the
irrotational-vortex funnel z_s = −Γ²/(8π²gR²), ours) and (8.12) (solid-body rotation); torque and power = dissipation.

Run: ``.venv/Scripts/python.exe scripts/ch08_pipe_and_couette.py --no-show``
Figures → outputs/ch08/c03_pipe.png, c04_circular_couette.png.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import quad

from ch08_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch03_kinematics as ch03  # noqa: F401
from fluidpy import ch05_vorticity_dynamics as ch05
from fluidpy import ch08_laminar_flow as ch08
from fluidpy.core import navier_stokes as NS
from fluidpy.core import vortices as VX
from fluidpy.core.thermo import G_BOOK


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    a, dpdz, mu, rho = 0.005, -40.0, 1e-3, 1000.0
    R = np.linspace(0, a, 201)
    u = ch08.pipe_poiseuille(R, a, dpdz, mu)
    tau = ch08.pipe_shear_stress(R, dpdz)
    Q, V, umax = ch08.pipe_flow_rate(a, dpdz, mu)
    Qn = quad(lambda r: float(ch08.pipe_poiseuille(r, a, dpdz, mu)) * 2 * np.pi * r, 0, a)[0]
    ref = NS.exact_solution("pipe_poiseuille", np.stack([R, 0 * R, 0 * R]), G=-dpdz, R=a, mu=mu)[0][2]
    Re = V * 2 * a * rho / mu
    tau0 = float(ch08.pipe_wall_stress(a, dpdz))
    print(f"pipe a = {a} m, dp/dz = {dpdz} Pa/m: Q = {Q:.5e} m³/s (quad {Qn:.5e}), V = {V:.4f} m/s, u_max/V = "
          f"{umax / V:.3f}, τ₀ = {tau0:.4f} Pa, Re = {Re:.1f} → {ch08.pipe_flow_regime(V, 2 * a, mu / rho)[1]}, "
          f"f = {float(ch08.pipe_friction_factor(Re)):.4f} (check 8τ₀/ρV² = {8 * abs(tau0) / (rho * V ** 2):.4f}); "
          f"parity with ch04 exact_solution |Δ| = {np.max(np.abs(u - ref)):.1e}")
    print(f"CV force balance: πa²|dp/dz| = {np.pi * a ** 2 * abs(dpdz):.6e} N/m, 2πa|τ₀| = {2 * np.pi * a * abs(tau0):.6e} N/m")
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.3))
    ax[0].plot(u, R * 1e3, color=COLORS["accent"], lw=2.4)
    ax[0].plot(u, -R * 1e3, color=COLORS["accent"], lw=2.4)
    ax[0].axvline(V, color=COLORS["teal"], ls="--", label=f"mean V = {V:.3f} m/s")
    ax[0].set_xlabel("u_z [m/s]")
    ax[0].set_ylabel("R [mm]")
    ax[0].legend(fontsize=8)
    ax[0].set_title("(8.6) parabolic profile")
    ax[1].plot(tau * 1e3, R * 1e3, color=COLORS["orange"], lw=2)
    ax[1].plot(tau * 1e3, -R * 1e3, color=COLORS["orange"], lw=2)
    ax[1].set_xlabel("τ = (R/2) dp/dz [mPa]")
    ax[1].set_title("(8.7) linear stress, wall value (8.8)")
    aa = np.logspace(-4, -2, 30)
    ax[2].loglog(aa * 1e3, [ch08.pipe_flow_rate(r, dpdz, mu)[0] for r in aa], color=COLORS["accent"], lw=2,
                 label="Q ∝ a⁴")
    ax[2].set_xlabel("a [mm]")
    ax[2].set_ylabel("Q [m³/s]")
    ax[2].legend(fontsize=8)
    slope = np.polyfit(np.log(aa), np.log([ch08.pipe_flow_rate(r, dpdz, mu)[0] for r in aa]), 1)[0]
    ax[2].set_title(f"Hagen–Poiseuille: log–log slope {slope:.3f}")
    print(f"Q(a) log–log slope = {slope:.6f}")
    save(fig, out, "c03_pipe")

    # circular Couette
    R1, R2 = 0.05, 0.1
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.3))
    Rr = np.linspace(R1, R2, 200)
    for O2, c in ((0.0, COLORS["accent"]), (1.0, COLORS["teal"]), (-1.0, COLORS["rose"]), (2.0, COLORS["orange"])):
        st = ch08.circular_couette_state(R1, R2, 2.0, O2, mu, rho)
        ax[0].plot(Rr * 100, ch08.circular_couette(Rr, R1, R2, 2.0, O2), color=c, lw=2,
                   label=f"Ω₂/Ω₁ = {O2 / 2:+.1f}, Rayleigh {'stable' if st['rayleigh_stable'] else 'unstable'}")
        print(f"Ω₁ = 2, Ω₂ = {O2}: A = {st['A']:.4f} 1/s, B = {st['B']:.3e} m²/s, torque(in) = {st['torque_inner']:.4e} "
              f"N m/m, power in = {st['power_in']:.4e} W/m = dissipation {st['dissipation']:.4e} W/m")
    ax[0].set_xlabel("R [cm]")
    ax[0].set_ylabel("u_φ [m/s]")
    ax[0].legend(fontsize=7)
    ax[0].set_title(f"(8.10), R₁ = {R1*100:.0f} cm, R₂ = {R2*100:.0f} cm, Ω₁ = 2 rad/s")
    Rb = np.linspace(R1, 6 * R1, 200)
    O1 = 3.0
    ub = ch08.circular_couette(Rb, R1, np.inf, O1, 0.0)
    ax[1].plot(Rb * 100, ub, color=COLORS["accent"], lw=2.4, label="(8.11) R₂ → ∞")
    ax[1].plot(Rb * 100, ch05.rotating_cylinder_flow(Rb, R1, 2 * O1)[0], "o", ms=3, color=COLORS["orange"],
               markevery=12, label="ch05 rotating_cylinder_flow(ω = 2Ω₁)")
    Gam = 2 * np.pi * O1 * R1 ** 2
    zs = -Gam ** 2 / (8 * np.pi ** 2 * G_BOOK * Rb ** 2)
    axb = ax[1].twinx()
    axb.plot(Rb * 100, zs * 1e3, color=COLORS["blue"], ls="--", lw=1.3)
    axb.set_ylabel("free surface z_s [mm] (dashed)", color=COLORS["blue"])
    ax[1].set_xlabel("R [cm]")
    ax[1].set_ylabel("u_φ [m/s]")
    ax[1].legend(fontsize=7)
    ax[1].set_title("cylinder in an infinite bath (Fig. 8.7)")
    pw = ch08.circular_couette_power(R1, np.inf, O1, 0.0, mu)
    print(f"bath: power in = {pw['power_in']:.6e} W/m, dissipation (quad) = {pw['dissipation_quad']:.6e} W/m, "
          f"4πμΩ₁²R₁² = {4 * np.pi * mu * O1 ** 2 * R1 ** 2:.6e}; |Δu| vs ch05 = "
          f"{np.max(np.abs(ub - ch05.rotating_cylinder_flow(Rb, R1, 2 * O1)[0])):.1e}")
    Rt = np.linspace(0, R2, 100)
    ax[2].plot(Rt * 100, ch08.circular_couette(Rt, 0.0, R2, 0.0, 1.5), color=COLORS["accent"], lw=2.4, label="(8.12)")
    ax[2].plot(Rt * 100, VX.solid_body_rotation(Rt, 1.5), ":", color=COLORS["orange"], lw=2,
               label="core.vortices.solid_body_rotation")
    ax[2].set_xlabel("R [cm]")
    ax[2].set_ylabel("u_φ [m/s]")
    ax[2].legend(fontsize=7)
    ax[2].set_title("tank in solid-body rotation, Ω₂ = 1.5 rad/s")
    save(fig, out, "c04_circular_couette")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
