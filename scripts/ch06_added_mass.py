"""§6.9 (C15): the arbitrarily moving sphere — kinematic condition on the surface, (6.104) against ∇φ (and the printed
bracket sign as the wrong variant), the surface pressure (6.105) split into its speed part (zero force, (6.106)) and
its acceleration part (force −M du_s/dt, (6.108)), parity with ch04's accelerating sphere, the added mass by two
independent routes (pressure force, kinetic energy), the cylinder's ρπa², a bubble starting at 2g and a ball with and
without added mass (6.109), and the Rayleigh collapse time (Exercise 6.44).

Run: ``.venv/Scripts/python.exe scripts/ch06_added_mass.py --no-show``
Figures → outputs/ch06/c15_added_mass.png.
"""
from __future__ import annotations

import warnings

import numpy as np

from ch06_drawings import parse_args, save, setup

from fluidpy import ch04_conservation_laws as ch04
from fluidpy import ch06_ideal_flow as ch06
from fluidpy.core.style import COLORS


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    a, rho = 0.1, 1000.0
    xs, us = np.array([0.2, 0.1, -0.3]), np.array([1.0, -0.5, 0.3])
    rng = np.random.default_rng(1)
    E = rng.normal(size=(3, 200))
    E /= np.linalg.norm(E, axis=0)
    u = ch06.moving_sphere_velocity(xs[:, None] + a * E, xs, us, a)
    print(f"kinematic BC on |ξ| = a (200 random points): max |u·n − u_s·n| = "
          f"{np.max(np.abs(np.sum(u * E, 0) - us @ E)):.1e} m/s")
    print(f"(6.104) vs ∇φ: max diff {np.max(np.abs(ch06.moving_sphere_surface_velocity(E, us) - u)):.1e};  printed bracket "
          f"sign: {np.max(np.abs(ch06.moving_sphere_surface_velocity(E, us, printed_bracket=True) - u)):.3f} (wrong)")
    ds = ch06.moving_sphere_dphidt_sym()
    print(f"sympy (6.101): chain − book {ds['chain']}, direct − book {ds['direct']}; (6.104) corrected bracket − ∇φ "
          f"{list(ds['correct_bracket'])}, printed bracket − ∇φ {list(ds['printed_bracket'])}")
    th = np.linspace(0, np.pi, 7)
    for U, dU in ((2.0, 0.0), (0.0, 3.0), (2.0, 3.0)):
        mine = np.array([ch06.moving_sphere_state(t, us=U, dus_dt=dU, a=a, rho=rho)["total"] for t in th])
        with warnings.catch_warnings():  # ch04's force quadrature warns when the force is exactly zero
            warnings.simplefilter("ignore")
            ref = ch04.accelerating_sphere_pressure(th, a=a, dUdt=dU, rho=rho, U=U)["p"]
        print(f"  (6.105) vs ch04.accelerating_sphere_pressure, u_s = {U}, du_s/dt = {dU}: max |Δp| {np.max(np.abs(mine - ref)):.1e} Pa")
    st = ch06.moving_sphere_state(0.4, us=2.0, dus_dt=3.0, accel_angle=0.7, a=a, rho=rho)
    M = ch06.added_mass_sphere(a, rho)
    print(f"oblique case (acceleration 40° to the velocity): F_steady = {np.round(st['F_steady'], 12)} N, "
          f"F_accel = {np.round(st['F_accel'], 10)} N, −M du_s/dt = {np.round(-M * 3 * np.array([np.cos(0.7), np.sin(0.7), 0]), 10)} N")
    Me = ch06.added_mass_by_energy(a, rho)
    print(f"added mass: pressure route M = 2πρa³/3 = {M:.10f} kg, energy route {Me:.10f} kg, M/(ρV) = "
          f"{M / (rho * 4 / 3 * np.pi * a ** 3):.6f}")
    print(f"cylinder added mass per length: closed {ch06.cylinder_added_mass(a, rho):.10f}, energy "
          f"{ch06.cylinder_added_mass(a, rho, 'energy'):.10f} kg/m (= ρπa²)")
    bub = ch06.sphere_motion(0.0, 0.005, rho, mode="bubble", t_eval=np.linspace(0, 0.05, 11))
    print(f"bubble (m → 0, a = 5 mm): initial acceleration {bub['a0']:.6f} m/s² = {bub['a0'] / ch06.G:.6f} g")
    Vb = 4 / 3 * np.pi * a ** 3
    t = np.linspace(0, 1.0, 101 if args.fast else 201)
    wood = ch06.sphere_motion(0.6 * rho * Vb, a, rho, mode="ball", t_eval=t)
    wood0 = ch06.sphere_motion(0.6 * rho * Vb, a, rho, mode="ball", t_eval=t, added_mass=False)
    print(f"wooden ball (ρ_b = 0.6ρ) released under water: a₀ = {wood['a0']:.4f} m/s² with added mass vs "
          f"{wood0['a0']:.4f} without; work − kinetic energy at 1 s: {wood['work'][-1] - wood['kinetic'][-1]:.1e} J")
    tc = ch06.rayleigh_collapse_time(1e-3, rho, 1e5)
    print(f"Rayleigh collapse, R₀ = 1 mm, Δp = 1 bar: t_c = {tc * 1e6:.4f} µs (closed) vs "
          f"{ch06.rayleigh_collapse_time(1e-3, rho, 1e5, 'quad') * 1e6:.4f} µs (quad); t_c/(R₀√(ρ/Δp)) = "
          f"{tc / (1e-3 * np.sqrt(rho / 1e5)):.5f}")

    fig, ax = plt.subplots(1, 3, figsize=(17, 4.4))
    tt = np.linspace(0, np.pi, 181)
    parts = [ch06.moving_sphere_state(q, us=2.0, dus_dt=3.0, a=a, rho=rho) for q in tt]
    ax[0].plot(np.degrees(tt), [p["steady"] for p in parts], color=COLORS["teal"], label="speed part (symmetric)")
    ax[0].plot(np.degrees(tt), [p["acceleration"] for p in parts], color=COLORS["orange"], label="acceleration part")
    ax[0].plot(np.degrees(tt), [p["total"] for p in parts], color=COLORS["ink"], ls="--", label="total (6.105)")
    ax[0].set_xlabel("θ_s from the sphere's velocity [deg]")
    ax[0].set_ylabel("p − p∞ [Pa]")
    ax[0].legend(fontsize=8)
    ax[0].set_title("u_s = 2 m/s, du_s/dt = 3 m/s², a = 0.1 m, water", fontsize=10)
    ax[1].plot(t, wood["u"], color=COLORS["accent"], label="with added mass (6.109)")
    ax[1].plot(t, wood0["u"], color=COLORS["muted"], ls="--", label="without (vacuum inertia)")
    ax[1].set_xlabel("t [s]")
    ax[1].set_ylabel("u [m/s]")
    ax[1].legend(fontsize=8)
    ax[1].set_title("wooden ball rising from rest (no drag)", fontsize=10)
    osc = ch06.sphere_motion(1.0, a, rho, mode="oscillating", amplitude=0.02, omega=2 * np.pi,
                             t_eval=np.linspace(0, 2, 201))
    ax[2].plot(osc["t"], osc["F_s"], color=COLORS["orange"], label="fluid force F_s = −M du/dt")
    ax[2].plot(osc["t"], osc["dudt"] * M, color=COLORS["accent"], ls=":", label="M × acceleration")
    ax[2].set_xlabel("t [s]")
    ax[2].set_ylabel("force [N]")
    ax[2].legend(fontsize=8)
    ax[2].set_title("oscillating sphere: force opposes the acceleration", fontsize=10)
    save(fig, out, "c15_added_mass")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
