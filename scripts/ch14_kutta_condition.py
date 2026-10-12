"""Chapter 14, §14.3 and §14.5 — the Kutta condition, the Kutta circulation (14.11) and the lift law (14.12).

* three stages of the start (Fig. 14.9): Γ/Γ_K = 0, ½, 1 — the rear stagnation point moves to the trailing edge;
* speed near the trailing edge against Γ: unbounded unless Γ = Γ_K, then U(b/a)cos(α + β) (ours, l'Hôpital on (14.9));
* Γ_K three ways: (14.11), a root finder on "no slip speed at ζ = b", and the contour integral of the velocity;
* the lift curve: 2π(α + β) (14.12), the exact 8π(a/c) sin(α + β), and a generic stalled curve (qualitative model).

Run: ``.venv/Scripts/python.exe scripts/ch14_kutta_condition.py --no-show``   Figures -> outputs/ch14/kutta_*.png
"""
from __future__ import annotations

import numpy as np
from ch14_common import C_CIRC, C_DRAG, C_LIFT, C_REF, COLORS, finish, parse_args, save, setup, streamlines

from fluidpy import ch14_aerodynamics as ch14


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch14.illustrative_inputs()
    b, e, beta, alpha, U = I["b"], I["e"], I["beta"], np.deg2rad(8.0), I["U_section"]
    mu, a = ch14.zhukhovsky_circle(b, e, beta)
    GK = ch14.kutta_circulation(U, a, alpha, beta)
    GKr = ch14.kutta_circulation_root(U, b, e, beta, alpha)
    F = ch14.zhukhovsky_flow(U, b, e, beta, alpha)
    th = np.linspace(0.0, 2 * np.pi, 720, endpoint=False)
    print(f"section b = {b} m, e = {e}, β = {np.rad2deg(beta):.0f}°, α = {np.rad2deg(alpha):.1f}°, U = {U} m/s: circle centre {mu:.4f}, radius a = {a:.5f} m")
    print(f"Kutta circulation (clockwise): (14.11) 4πUa sin(α + β) = {GK:.10f}; root finder on u_θ(ζ = b) = 0: {GKr:.10f} "
          f"(difference {abs(GK - GKr):.1e}); dw/dζ at ζ = b with it: {abs(F.dwdzeta(b)):.1e}")
    for R in (3.0 * b, 6.0 * b, 20.0 * b):
        zc = R * np.exp(1j * th)
        circ = float(np.sum(F.dwdz(zc) * 1j * zc).real * 2 * np.pi / th.size)
        print(f"  ∮u·dx on the circle |z| = {R:4.1f} (counterclockwise): {circ:+.10f}  (−Γ_K = {-GK:+.10f})")
    zf = np.array([400.0 + 300.0j])
    print(f"  far field: dw/dz at |z| = 500 m is {complex(F.dwdz(zf)[0]):.5f}; U e^(−iα) = {U * np.exp(-1j * alpha):.5f}")

    # --- three stages (Fig. 14.9) ------------------------------------------------------------------------------------
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
    for k, frac in enumerate((0.0, 0.5, 1.0)):
        Fk = ch14.zhukhovsky_flow(U, b, e, beta, alpha, Gamma_cw=frac * GK)
        streamlines(ax[k], Fk, (-3.2 * b, 3.6 * b), (-1.9 * b, 1.9 * b), levels=29)
        st = ch14.stagnation_points(U, b, e, beta, alpha, Gamma_cw=frac * GK)
        rs = ch14.rear_stagnation_point(U, a, alpha, beta, Gamma_cw=frac * GK)
        ax[k].plot([st["z_front"].real, st["z_rear"].real], [st["z_front"].imag, st["z_rear"].imag], "o", color=C_DRAG, ms=6, zorder=6)
        ax[k].set(xlabel="x [m]", ylabel="y [m]", title=f"Γ/Γ_K = {frac:.1f}: lift ρUΓ = {ch14.lift_per_span(I['rho'], U, Gamma_cw=frac * GK):.0f} N/m")
        print(f"Γ/Γ_K = {frac:.1f}: rear stagnation point at z = {st['z_rear']:.4f} ({np.rad2deg(rs['offset']):+.2f}° of the circle ahead of the "
              f"edge point B), front at z = {st['z_front']:.4f}; speed 1e-4 b behind the edge: {ch14.edge_speed(U, b, e, beta, alpha, 1e-4 * b, Gamma_cw=frac * GK) / U:.3f} U")
    save(fig, out, "kutta_three_stages")

    # --- edge speed against circulation ----------------------------------------------------------------------------
    fr = np.linspace(0.0, 1.5, 151)
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.4))
    for d, col in zip((1e-2, 1e-3, 1e-4), (C_REF, C_CIRC, C_DRAG)):
        sp = np.array([ch14.edge_speed(U, b, e, beta, alpha, d * b, Gamma_cw=q * GK) for q in fr])
        ax[0].semilogy(fr, sp / U, color=col, label=f"{d:g} b behind the edge")
    te = ch14.trailing_edge_speed(U, a, b, alpha, beta)
    ax[0].axhline(te / U, color=C_LIFT, ls="--", lw=1.2, label="U(b/a)cos(α + β) at the edge, Γ = Γ_K")
    ax[0].set(xlabel=r"$\Gamma/\Gamma_K$ [–]", ylabel="speed / U [–]", title="speed near the trailing edge: finite only at the Kutta value")
    ax[0].legend(fontsize=8)
    print(f"finite edge speed with the Kutta value: U(b/a)cos(α + β) = {te / U:.6f} U = {te:.4f} m/s; approached from behind the edge:")
    for d in (1e-2, 1e-4, 1e-6, 1e-8):
        print(f"  d = {d:.0e} b: {ch14.edge_speed(U, b, e, beta, alpha, d * b) / U:.6f} U   (with 0.9 Γ_K: {ch14.edge_speed(U, b, e, beta, alpha, d * b, Gamma_cw=0.9 * GK) / U:.3f} U)")

    # --- lift curve ---------------------------------------------------------------------------------------------------
    al = np.deg2rad(np.linspace(-8.0, 22.0, 121))
    cl_thin = ch14.zhukhovsky_lift_coefficient(al, beta)
    cl_ex = ch14.zhukhovsky_lift_coefficient(al, beta, e, exact=True)
    cl_gen = ch14.generic_lift_curve(al, K=I["K"], beta=beta, alpha_stall=I["alpha_stall"], CL_post=I["CL_post"], width=I["stall_width"])
    ax[1].plot(np.rad2deg(al), cl_thin, color=C_LIFT, label=r"$2\pi(\alpha+\beta)$ (14.12)")
    ax[1].plot(np.rad2deg(al), cl_ex, color=C_LIFT, ls="--", label=rf"exact $8\pi(a/c)\sin(\alpha+\beta)$, e = {e}")
    ax[1].plot(np.rad2deg(al), cl_gen, color=C_DRAG, label="generic measured curve (qualitative model, our constants)")
    ax[1].axvline(-np.rad2deg(beta), color=C_REF, lw=0.8, ls=":")
    ax[1].axhline(0.0, color=C_REF, lw=0.8)
    ax[1].set(xlabel="angle of attack α [deg]", ylabel=r"$C_L$ [–]", title=f"lift curve: zero lift at α = −β = {-np.rad2deg(beta):.2f}°")
    ax[1].legend(fontsize=8)
    save(fig, out, "kutta_edge_speed_lift_curve")
    sh = ch14.zhukhovsky_shape(b, e, beta)
    print(f"lift at α = {np.rad2deg(alpha):.0f}°: 2π(α + β) = {ch14.zhukhovsky_lift_coefficient(alpha, beta):.5f}; exact = "
          f"{ch14.zhukhovsky_lift_coefficient(alpha, beta, e, True):.5f} = 2Γ_K/(Uc) = {2 * GK / (U * sh['chord']):.5f} (chord {sh['chord']:.5f}); "
          f"exact slope at small angles 8πa/c = {8 * np.pi * a / sh['chord']:.4f} per rad ({8 * np.pi * a / sh['chord'] / (2 * np.pi) - 1:+.1%} above 2π)")
    print(f"zero-lift angle from camber/chord = {sh['camber_ratio']:.5f}: {np.rad2deg(ch14.zero_lift_angle(sh['camber_ratio'])):.3f}° (−β = {-np.rad2deg(beta):.3f}°)")
    fp = ch14.flat_plate_flow(U, b, alpha)
    print(f"flat plate (e = β = 0): Γ_K = {fp.Gamma_kutta:.6f} = 4πUb sin α = {4 * np.pi * U * b * np.sin(alpha):.6f}; C_L = "
          f"{ch14.zhukhovsky_force(U, b, 0.0, 0.0, alpha, 1.0)['CL']:.6f} = 2π sin α = {2 * np.pi * np.sin(alpha):.6f}")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
