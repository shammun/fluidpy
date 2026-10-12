"""Chapter 14, §14.3 and §14.5 — surface pressure, the C_p loop and the force on a Zhukhovsky section.

* C_p round the section (Fig. 14.8, computed): arrows normal to the surface and −C_p against x/c;
* lift four ways: ρUΓ (6.62), the Blasius contour integral (6.60), the surface-pressure integral (with its
  convergence), and the area of the C_p loop;
* circulation as the integral of the surface speed; d'Alembert (zero drag); centre of pressure;
* flat plate: leading-edge suction from a control volume; thin-airfoil theory; a qualitative separation estimate.

Run: ``.venv/Scripts/python.exe scripts/ch14_pressure_lift.py --no-show``   Figures -> outputs/ch14/pressure_*.png
"""
from __future__ import annotations

import numpy as np
from ch14_common import C_DRAG, C_FLOW, C_LIFT, C_PRES, C_REF, COLORS, finish, parse_args, save, setup

from fluidpy import ch14_aerodynamics as ch14


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch14.illustrative_inputs()
    b, e, beta, alpha, U, rho = I["b"], I["e"], I["beta"], I["alpha"], I["U_section"], I["rho"]
    S = ch14.zhukhovsky_surface(U, b, e, beta, alpha, 2000)
    GK = S["Gamma_cw"]

    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
    ax[0].fill(S["x"], S["y"], color="#d5d9e4")
    ax[0].plot(np.append(S["x"], S["x"][0]), np.append(S["y"], S["y"][0]), color=COLORS["ink"], lw=1.2)
    dz = np.gradient(np.append(S["z"], S["z"][:2]))[:-2]
    nrm = -1j * dz / np.abs(dz)                                       # outward normal of a counterclockwise contour
    for k in range(5, S["x"].size, 25):
        cp = S["cp"][k]
        tip = S["z"][k] + nrm[k] * 0.22 * b * abs(cp)
        a0, a1 = (S["z"][k], tip) if cp < 0 else (tip, S["z"][k])   # suction points outward, pressure inward
        ax[0].annotate("", (a1.real, a1.imag), (a0.real, a0.imag), arrowprops=dict(arrowstyle="->", color=C_PRES if cp < 0 else C_FLOW, lw=1))
    ax[0].set(aspect="equal", xlabel="x [m]", ylabel="y [m]", title=f"pressure on the surface (orange = suction), α = {np.rad2deg(alpha):.0f}°", ylim=(-1.2 * b, 1.5 * b))
    up, lo = S["upper"], S["lower"]
    ax[1].plot(up["xc"], -up["cp"], color=C_PRES, label="upper surface")
    ax[1].plot(lo["xc"], -lo["cp"], color=C_FLOW, label="lower surface")
    ax[1].fill_between(up["xc"], -up["cp"], np.interp(up["xc"], lo["xc"], -lo["cp"]), color=C_LIFT, alpha=0.18)
    ax[1].axhline(-1.0, color=C_REF, lw=0.8, ls=":")
    ax[1].set(xlabel="x/c [–]", ylabel=r"$-C_p$ [–]", title="the loop: its area is the normal-force coefficient")
    ax[1].legend(fontsize=8)

    cy = ch14.normal_force_coefficient(up["xc"], up["cp"], lo["xc"], lo["cp"])
    Fb = ch14.zhukhovsky_force(U, b, e, beta, alpha, rho, "blasius")
    q = 0.5 * rho * U ** 2
    print(f"section b = {b} m, e = {e}, β = {np.rad2deg(beta):.0f}°, α = {np.rad2deg(alpha):.1f}°, Γ_K = {GK:.6f} m²/s (clockwise), ρ = {rho} kg/m³, U = {U} m/s")
    print(f"  suction peak C_p = {S['cp_min']:.3f} at x/c = {S['xc_cp_min']:.4f}; trailing-edge speed {S['speed'][0] / U:.5f} U "
          f"(= U(b/a)cos(α + β) = {ch14.trailing_edge_speed(U, S['a'], b, alpha, beta) / U:.5f} U); C_p at the front stagnation point: "
          f"{float(ch14.zhukhovsky_flow(U, b, e, beta, alpha).cp(S['stagnation'][0])):.12f}")
    print(f"  lift per span: ρUΓ = {ch14.lift_per_span(rho, U, Gamma_cw=GK):.8f}; Blasius L = {Fb['L']:.8f}, D = {Fb['D']:+.1e}; "
          f"(Fx, Fy) = ({Fb['Fx']:+.6f}, {Fb['Fy']:+.6f}) — the force is normal to the stream, not to the chord")
    print(f"  C_p loop: C_y = {cy:.6f} on the x-extent {S['chord_x']:.5f} → F_y = {cy * q * S['chord_x']:.6f} (Blasius F_y {Fb['Fy']:.6f}); "
          f"C_L = C_y/cos α · (x-extent/chord) = {ch14.lift_from_cp(up['xc'], up['cp'], lo['xc'], lo['cp'], alpha) * S['chord_x'] / Fb['chord']:.6f}, "
          f"2Γ/(Uc) = {2 * GK / (U * Fb['chord']):.6f}, 2π(α + β) = {ch14.zhukhovsky_lift_coefficient(alpha, beta):.6f}")
    print(f"  circulation from the surface speed: polygon {ch14.surface_circulation(S):.8f}, ζ-plane {ch14.surface_circulation(S, 'circle'):.8f}")
    print(f"  moment about the leading edge {Fb['M_le']:+.6f} N (counterclockwise +); centre of pressure at {Fb['xcp_over_c']:.4f} c")
    print("  surface-pressure route against Blasius (lift error, drag):")
    errs, ns = [], (250, 500, 1000, 2000, 4000)
    for n in ns:
        Fp = ch14.zhukhovsky_force(U, b, e, beta, alpha, rho, "pressure", n=n)
        errs.append(abs(Fp["L"] - Fb["L"]))
        print(f"    n = {n:5d}: |L − L_Blasius| = {errs[-1]:.2e}, D = {Fp['D']:+.2e}")
    print(f"    observed order {np.round(np.log2(np.array(errs[:-1]) / np.array(errs[1:])), 2)}")
    ax[2].loglog(ns, errs, "o-", color=C_LIFT, label="|L − ρUΓ|, pressure integral")
    ax[2].loglog(ns, errs[0] * (ns[0] / np.array(ns)) ** 2, color=C_REF, ls="--", label="slope −2")
    ax[2].set(xlabel="surface points n [–]", ylabel="error [N/m]", title="pressure integral converges to ρUΓ")
    ax[2].legend(fontsize=8)
    save(fig, out, "pressure_cp_loop")

    # --- angle of attack sweep -----------------------------------------------------------------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.4))
    print("angle sweep (loop area against 2Γ/(Uc) and the separation estimate — qualitative, laminar Thwaites on the ideal pressure):")
    xs = []
    alist = np.deg2rad([0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
    for k, al in enumerate(alist):
        Sk = ch14.zhukhovsky_surface(U, b, e, beta, al, 1200)
        ax[0].plot(Sk["upper"]["xc"], -Sk["upper"]["cp"], color=plt.cm.viridis(k / 5), label=f"α = {np.rad2deg(al):.0f}°")
        ax[0].plot(Sk["lower"]["xc"], -Sk["lower"]["cp"], color=plt.cm.viridis(k / 5), lw=1, ls="--")
        sep = ch14.suction_side_separation(al, e, beta, 3.0e5)
        xs.append(sep["x_sep_over_c"])
        cyk = ch14.normal_force_coefficient(Sk["upper"]["xc"], Sk["upper"]["cp"], Sk["lower"]["xc"], Sk["lower"]["cp"])
        Fk = ch14.zhukhovsky_force(U, b, e, beta, al, rho)
        print(f"  α = {np.rad2deg(al):4.1f}°: suction peak {Sk['cp_min']:7.3f} at x/c = {Sk['xc_cp_min']:.4f}; C_y·(x-extent)/(c cos α) = "
              f"{cyk * Sk['chord_x'] / (Fk['chord'] * np.cos(al)):.5f}, C_L = {Fk['CL']:.5f}; laminar separation at x/c = {sep['x_sep_over_c']:.3f}")
    ax[0].set(xlabel="x/c [–]", ylabel=r"$-C_p$ [–]", title="suction peak sharpens with α (dashed: lower surface)", ylim=(-1.2, 6.0))
    ax[0].legend(fontsize=7)
    ax[1].plot(np.rad2deg(alist), xs, "o-", color=C_DRAG)
    ax[1].set(xlabel="angle of attack α [deg]", ylabel="separation position x/c [–]", title="laminar separation estimate — qualitative", ylim=(0, 1))
    save(fig, out, "pressure_alpha_sweep")
    sym_u = ch14.suction_side_separation(0.0, e, 0.0, 3.0e5, side="upper")["x_sep_over_c"]
    sym_l = ch14.suction_side_separation(0.0, e, 0.0, 3.0e5, side="lower")["x_sep_over_c"]
    print(f"  symmetric section at α = 0: separation at x/c = {sym_u:.6f} (upper) and {sym_l:.6f} (lower)")

    # --- flat plate, leading-edge suction, thin-airfoil theory ---------------------------------------------------------
    c, a5 = 4.0 * b, np.deg2rad(5.0)
    fp = ch14.flat_plate_forces(rho, U, c, a5)
    Fpl = ch14.zhukhovsky_force(U, b, 0.0, 0.0, a5, rho)
    cv = ch14.edge_flow_suction(rho, U * np.sin(a5), c, eps=1e-3)
    cv2 = ch14.edge_flow_suction(rho, U * np.sin(a5), c, eps=1e-1)
    print(f"flat plate at 5°: L = {fp['L']:.6f} (Blasius {Fpl['L']:.6f}), normal force N = L cos α = {fp['N']:.6f}, leading-edge suction "
          f"S = L sin α = {fp['S']:.6f}; N cos α + S sin α = {fp['N'] * np.cos(a5) + fp['S'] * np.sin(a5):.6f}, drag N sin α − S cos α = "
          f"{fp['N'] * np.sin(a5) - fp['S'] * np.cos(a5):+.1e}; centre of pressure at {Fpl['xcp_over_c']:.4f} c")
    print(f"  control volume round the edge (U_o²a = U²c sin²α): suction = {-cv['suction'][0]:.6f} with ε = 1e-3, {-cv2['suction'][0]:.6f} "
          f"with ε = 1e-1; closed form πρU_o²a = {cv['closed_form']:.6f}")
    h = ch14.zhukhovsky_shape(b, 0.0, beta)["camber_ratio"]
    ta = ch14.thin_airfoil(lambda x: 4.0 * h * (1.0 - 2.0 * x), alpha)
    print(f"thin-airfoil theory (ours), parabolic camber of height {h:.5f} c: α_L0 = {np.rad2deg(ta['alpha_L0']):.4f}° "
          f"(−2h/c = {np.rad2deg(-2 * h):.4f}°, Zhukhovsky −β = {np.rad2deg(-beta):.4f}°), C_L = {ta['CL']:.5f}, C_M,c/4 = {ta['CM_c4']:.5f}")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
