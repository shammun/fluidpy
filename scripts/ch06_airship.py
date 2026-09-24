"""§6.8 (C13, C14): axisymmetric ideal flow — the Stokes-stream-function field equation (6.77) (not the plane
Laplacian), the flux 2πdψ (6.78), the sphere (6.89)–(6.91) against the cylinder, the source + line-sink airship
(6.93)–(6.95) in body and fluid frames, and the axial singularity method recovering the exact linear source
distribution of a prolate spheroid as N grows (and its growing condition number).

Run: ``.venv/Scripts/python.exe scripts/ch06_airship.py --no-show``
Figures → outputs/ch06/c13_sphere_airship.png, c14_axial_method.png.
"""
from __future__ import annotations

import numpy as np

from ch06_drawings import parse_args, save, setup

from fluidpy import ch06_ideal_flow as ch06
from fluidpy.core import vortices as VX
from fluidpy.core.style import COLORS


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    U, a = 1.0, 1.0
    sph = ch06.sphere(U, a)
    print(f"sphere: Stokes operator (6.77) of ψ at (R, z) = (1.3, 0.4): {ch06.stokes_operator_residual(sph.psi, 1.3, 0.4):.1e};"
          f" control ψ = Rz (plane-harmonic, not a Stokes solution): {ch06.stokes_operator_residual(lambda R, z: R * z, 1.3, 0.4):.4f}; "
          f"axisym Laplacian of φ: {ch06.axisym_laplacian_residual(sph.phi, 1.3, 0.4):.1e}")
    fq, ff = ch06.axisym_flux_between(sph, (1.5, 0.3), (2.0, -1.0))
    print(f"  flux between stream surfaces (6.78): quadrature {fq:.12f} vs 2πΔψ {ff:.12f} m³/s")
    Rz = (np.array([1.2, 2.0, 0.5]), np.array([0.8, -1.5, 3.0]))
    hill = np.asarray(VX.hill_stream_function(Rz[0], Rz[1], 15.0 * U / (2.0 * a ** 2), a, outside=True))
    print(f"  parity with ch05 Hill exterior (A = 15U/2a², co-moving frame streams along −z, so ψ_Hill = −ψ_sphere): "
          f"max |ψ_sphere + ψ_Hill| = {np.max(np.abs(np.asarray(sph.psi(*Rz)) + hill)):.1e}")
    th = np.linspace(0, np.pi, 7)
    ur, ut = sph.velocity_spherical(a * np.ones_like(th), th)
    print(f"  surface: max |u_r| {np.max(np.abs(ur)):.1e}, max speed {np.max(np.abs(ut)):.6f} U (1.5), "
          f"C_p min {ch06.sphere_surface_cp(np.pi / 2):.4f} (cylinder −3)")
    print(f"  1 % perturbation radius: cylinder r/a = {ch06.perturbation_radius('cylinder'):.3f}, "
          f"sphere {ch06.perturbation_radius('sphere'):.3f}")
    print(f"line sink (6.93) vs (6.94) at (R, z) = (0.7, 1.3): {ch06.line_sink_stream_function(0.7, 1.3, 2.0, 1.0, 'quad'):.12f}"
          f" vs {ch06.line_sink_stream_function(0.7, 1.3, 2.0, 1.0):.12f} m³/s")
    air = ch06.airship(U, 1.0, 1.0)
    print(f"airship U = 1, Q = 1 m³/s, a = 1 m: nose {air['z_nose']:.6f} m, tail {air['z_tail']:.6f} m, length "
          f"{air['length']:.6f} m, max radius {air['R_max']:.6f} m, net strength {air['net_strength']:.1e}")
    zs, Rs = air["contour"]
    print(f"  max |ψ| on the computed body contour: {np.max(np.abs(np.asarray(air['psi'](Rs[1:-1], zs[1:-1])))):.1e}")

    Ns = [4, 8, 16, 32] if args.fast else [4, 8, 16, 32, 64]
    fits = [ch06.axisym_body_fit("ellipsoid", N) for N in Ns]
    for N, f in zip(Ns, fits):
        s = f["solution"]
        print(f"  axial method, prolate spheroid A = 1, B = 0.2, N = {N:3d}: k error {f['k_error']:.3e}, body ψ error "
              f"{f['body_psi_error']:.2e}, cond {s['cond']:.1e}, Σk_nΔξ {s['net_strength']:.1e}")
    ke = [f["k_error"] for f in fits]
    print(f"  observed order of the k error (last three N): {np.polyfit(np.log(Ns[-3:]), np.log(ke[-3:]), 1)[0] * -1:.3f}")
    for kind in ("rankine_oval", "airship", "ellipsoid", "sphere"):
        for N in (10, 20, 40):
            s = ch06.axial_state(kind, N=N)
            print(f"  axial_state {kind:<12s} N = {N:2d} (nose-to-tail segments): body error {s['body_error']:.2e} m, "
                  f"Σk_nΔξ {s['net_strength']:.1e}, cond₁ {s['cond']:.1e}")

    n = 121 if args.fast else 201
    fig, ax = plt.subplots(1, 3, figsize=(18, 4.4))
    zg = np.linspace(-3, 3, n)
    Rg = np.linspace(0, 2.5, n // 2)
    Z, R = np.meshgrid(zg, Rg, indexing="xy")
    ax[0].contour(Z, R, np.asarray(sph.psi(R, Z)), levels=np.linspace(0.05, 3, 16), colors=COLORS["ink"], linewidths=0.8)
    cyl = ch06.cylinder(U, a)
    ax[0].contour(Z, -R, np.asarray(cyl.psi(Z, R)), levels=np.linspace(0.1, 2.5, 13), colors=COLORS["muted"],
                  linewidths=0.8)
    t = np.linspace(0, np.pi, 100)
    ax[0].fill(np.cos(t), np.sin(t), color="#d9dce4")
    ax[0].fill(np.cos(t), -np.sin(t), color="#eceef3")
    ax[0].set_aspect("equal")
    ax[0].set_title("sphere (top, Stokes ψ) vs cylinder (bottom, plane ψ)", fontsize=10)
    ax[0].set_xlabel("z (along the stream) [m]")
    ax[0].set_ylabel("R [m]")
    zg2 = np.linspace(-1.5, 2.5, n)
    Rg2 = np.linspace(1e-3, 1.5, n // 2)
    Z2, R2 = np.meshgrid(zg2, Rg2, indexing="xy")
    ax[1].contour(Z2, R2, np.asarray(air["psi"](R2, Z2)), levels=np.linspace(-0.08, 0.8, 23), colors=COLORS["ink"],
                  linewidths=0.8, negative_linestyles="solid")
    ax[1].contour(Z2, -R2, np.asarray(air["psi_fluid_frame"](R2, Z2)), levels=np.linspace(-0.08, 0.08, 17),
                  colors=COLORS["accent"], linewidths=0.8, negative_linestyles="solid")
    ax[1].fill_between(zs, Rs, -Rs, color="#d9dce4")
    ax[1].plot([0, 1], [0, 0], color=COLORS["rose"], lw=3, label="line sink k = Q/a")
    ax[1].plot([0], [0], "o", color=COLORS["teal"], label="source Q")
    ax[1].set_aspect("equal")
    ax[1].legend(fontsize=8, loc="lower right")
    ax[1].set_title("airship (6.95): body frame (top), fluid frame (bottom)", fontsize=10)
    ax[1].set_xlabel("z [m]")
    f = fits[2]
    xm = f["solution"]["xi_mid"]
    ax[2].step(xm, f["solution"]["k"], where="mid", color=COLORS["teal"], label=f"k_n, N = {Ns[2]}")
    ax[2].plot(xm, f["target"]["exact_k"], "--", color=COLORS["ink"], label="exact k = Kξ")
    ax[2].set_xlabel("ξ along the axis [m]")
    ax[2].set_ylabel("source density k [m²/s]")
    ax[2].legend(fontsize=8)
    ax[2].set_title("axial singularity method: prolate spheroid", fontsize=10)
    save(fig, out, "c13_sphere_airship")

    fig2, axc = plt.subplots(1, 1, figsize=(6, 4))
    axc.loglog(Ns, ke, "o-", color=COLORS["teal"], label="k error")
    axc.loglog(Ns, [f["body_psi_error"] for f in fits], "s-", color=COLORS["accent"], label="body ψ error")
    axc.loglog(Ns, [f["solution"]["cond"] * 1e-12 for f in fits], "^--", color=COLORS["rose"], label="cond × 10⁻¹²")
    axc.set_xlabel("N segments")
    axc.legend(fontsize=8)
    axc.set_title("axial method convergence (spheroid A = 1, B = 0.2)", fontsize=10)
    save(fig2, out, "c14_axial_method")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
