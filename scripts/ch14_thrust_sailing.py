"""Chapter 14, §14.8–§14.9 — thrust from lift: the fish tail, the bird's wing and sailing against the wind.

* the oscillating foil (Fig. 14.27): thrust L sin θ, side force L cos θ, and the two powers — equal without drag,
  efficiency below 1 with drag (ours);
* the apparent wind U_r = U − V, sail lift and drag resolved into drive and lateral force (Fig. 14.29); the drive
  vanishes at tan β_a = D/L; closest course β_a = ε_A + ε_H; steady boat speed and velocity made good;
* the projectile of Exercise 14.1 (formula against a numerical integration, our numbers).

Run: ``.venv/Scripts/python.exe scripts/ch14_thrust_sailing.py --no-show``   Figures -> outputs/ch14/thrust_*.png
"""
from __future__ import annotations

import numpy as np
from ch14_common import C_DRAG, C_DRIVE, C_FLOW, C_LIFT, C_REF, COLORS, finish, parse_args, save, setup

from fluidpy import ch14_aerodynamics as ch14


def arrow(ax, p0, v, color, label=None, lw=2.0):
    ax.annotate("", (p0[0] + v[0], p0[1] + v[1]), p0, arrowprops=dict(arrowstyle="->", color=color, lw=lw))
    if label:
        ax.text(p0[0] + 1.08 * v[0], p0[1] + 1.08 * v[1], label, color=color, fontsize=8, ha="center", va="center")


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch14.illustrative_inputs()
    # --- fish tail ------------------------------------------------------------------------------------------------
    Uf, Vt, L = I["U_fish"], I["V_tail"], I["L_tail"]
    f0 = ch14.oscillating_foil(Uf, Vt, L)
    f1 = ch14.oscillating_foil(Uf, Vt, L, L / I["LD_tail"])
    print(f"fish tail: forward U = {Uf} m/s, sideways V = {Vt} m/s → θ = tan⁻¹(V/U) = {np.rad2deg(f0['theta']):.2f}°, relative speed {f0['Ur']:.3f} m/s, lift {L} N")
    print(f"  no drag : thrust L sin θ = {f0['thrust']:.4f} N at the rate {f0['power_out']:.4f} W; side force L cos θ = {f0['side_force']:.4f} N at the rate "
          f"{f0['power_in']:.4f} W; efficiency {f0['efficiency']:.12f}")
    print(f"  L/D = {I['LD_tail']}: thrust {f1['thrust']:.4f} N, side force {f1['side_force']:.4f} N, efficiency {f1['efficiency']:.4f}; thrust vanishes at "
          f"θ = tan⁻¹(D/L) = {np.rad2deg(np.arctan(1 / I['LD_tail'])):.2f}°")
    dn = ch14.oscillating_foil(9.0, 3.5, 4.0, 0.3)
    up = ch14.oscillating_foil(9.0, 5.5, 0.4, 0.1)
    print(f"bird wing (vertical strokes, our numbers): downstroke thrust {dn['thrust']:.3f} N, vertical force {dn['side_force']:.3f} N; upstroke (wing nearly "
          f"feathered, small lift) thrust {up['thrust']:+.3f} N, vertical force {up['side_force']:.3f} N")

    fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
    ratio = np.linspace(0.02, 2.0, 200)
    for LD, col in zip((np.inf, 30.0, 12.0, 5.0), (COLORS["ink"], C_LIFT, C_DRIVE, C_DRAG)):
        r = ch14.oscillating_foil(1.0, ratio, 1.0, 0.0 if np.isinf(LD) else 1.0 / LD)
        ax[0].plot(ratio, np.clip(r["efficiency"], 0, None), color=col, label="no drag (book: ideally perfect)" if np.isinf(LD) else f"L/D = {LD:.0f}")
    ax[0].set(xlabel="stroke speed ratio V/U = tan θ [–]", ylabel="thrust power / input power [–]", title="oscillating foil: efficiency with drag (ours)", ylim=(0, 1.05))
    ax[0].legend(fontsize=8)

    # --- sailing ----------------------------------------------------------------------------------------------------
    Uw, VoU, CLs, CDs, As, rho = I["U_wind"], I["V_over_U"], I["CL_sail"], I["CD_sail"], I["A_sail"], I["rho"]
    th = I["theta"]
    aw = ch14.apparent_wind(Uw, VoU * Uw, th)
    sf = ch14.sail_forces(Uw, VoU * Uw, th, CL=CLs, CD=CDs, rho=rho, A=As)
    print(f"sailing: true wind {Uw} m/s, boat speed {VoU * Uw:.2f} m/s on a course {np.rad2deg(th):.0f}° off the wind → apparent wind {aw['Ur']:.3f} m/s at "
          f"{np.rad2deg(aw['beta_a']):.2f}° off the bow (components along/across the course {aw['vector'][0]:+.3f}, {aw['vector'][1]:+.3f} m/s)")
    print(f"  sail lift {sf['L']:.1f} N ⟂ apparent wind, drag {sf['D']:.1f} N ∥ it → drive = L sin β_a − D cos β_a = {sf['drive']:.1f} N, lateral = "
          f"{sf['lateral']:.1f} N; drive²+lateral² − (L²+D²) = {sf['drive'] ** 2 + sf['lateral'] ** 2 - sf['L'] ** 2 - sf['D'] ** 2:.1e}")
    b0 = np.arctan(CDs / CLs)
    thetas = np.deg2rad(np.linspace(0.0, 180.0, 721))
    drive = ch14.sail_drive_curve(thetas, U=Uw, V=VoU * Uw, CL=CLs, CD=CDs, rho=rho, A=As)
    k0 = int(np.argmax(drive > 0))
    ba0 = ch14.apparent_wind(Uw, VoU * Uw, thetas[k0])["beta_a"]
    print(f"  drive changes sign at θ = {np.rad2deg(thetas[k0]):.2f}° off the true wind, where β_a = {np.rad2deg(ba0):.2f}° ≈ tan⁻¹(D/L) = {np.rad2deg(b0):.2f}°")
    lim = ch14.sailboat_course_limit(CLs / CDs, I["LD_hull"])
    print(f"  closest steady course to the apparent wind with a hull of L/D = {I['LD_hull']}: ε_A + ε_H = {np.rad2deg(b0):.2f}° + "
          f"{np.rad2deg(np.arctan(1 / I['LD_hull'])):.2f}° = {np.rad2deg(lim):.2f}°")
    for LD, col in zip((3.0, CLs / CDs, 15.0), (C_DRAG, C_DRIVE, C_LIFT)):
        d = ch14.sail_drive_curve(thetas, U=Uw, V=VoU * Uw, CL=CLs, CD=CLs / LD, rho=rho, A=As)
        ax[1].plot(np.rad2deg(thetas), d, color=col, label=f"sail L/D = {LD:.1f}")
    ax[1].axhline(0.0, color=C_REF, lw=0.8)
    ax[1].axvspan(0.0, np.rad2deg(thetas[k0]), color=C_REF, alpha=0.2)
    ax[1].set(xlabel="course angle off the true wind θ [deg]", ylabel="driving force [N]", title=f"sail drive at V/U = {VoU} (shaded: no-go zone of our sail)")
    ax[1].legend(fontsize=8)

    print("  steady sailing (drive = hull drag; the keel supplies the lateral force):")
    print("    θ [deg] | boat speed [m/s] | V/U   | β_a [deg] | hull L/D needed | keel holds (L/D ≤ 7) | velocity made good [m/s] | polygon closure")
    best = (0.0, 0.0)
    eqs = []
    for t in np.deg2rad(np.arange(20.0, 181.0, 2.0)):
        q = ch14.sailboat_equilibrium(Uw, t, CL=CLs, CD=CDs, rho_air=rho, A_sail=As, rho_water=I["rho_water"], CDA_hull=I["CDA_hull"], LD_hull=I["LD_hull"])
        eqs.append((t, q))
        if q["keel_holds"] and q["vmg"] > best[1]:
            best = (t, q["vmg"])
    for t, q in eqs[::8]:
        print(f"    {np.rad2deg(t):6.1f}  | {q['V']:8.3f}         | {q['V'] / Uw:.3f} | {np.rad2deg(q['beta_a']):7.2f}   | {q['LD_hull_required']:8.2f}        | "
              f"{str(q['keel_holds']):5s}                | {q['vmg']:+8.3f}                 | {q['closure']:+.1e}")
    print(f"  best course to windward that the keel can hold: θ = {np.rad2deg(best[0]):.0f}° with {best[1]:.3f} m/s made good ({best[1] / Uw:.3f} U)")

    # --- force polygon ---------------------------------------------------------------------------------------------
    q = ch14.sailboat_equilibrium(Uw, th, CL=CLs, CD=CDs, rho_air=rho, A_sail=As, rho_water=I["rho_water"], CDA_hull=I["CDA_hull"])
    ba = q["beta_a"]
    sfe = ch14.sail_forces(Uw, q["V"], th, CL=CLs, CD=CDs, rho=rho, A=As)
    sc = 1.0 / sfe["F"]
    wind_dir = -np.array([np.cos(ba), np.sin(ba)])                    # direction the apparent wind blows toward (course along +x)
    Lv = sfe["L"] * sc * np.array([np.sin(ba), -np.cos(ba)])
    Dv = sfe["D"] * sc * wind_dir
    a = ax[2]
    arrow(a, (0, 0), Lv, C_LIFT, "sail lift")
    arrow(a, (0, 0), Dv, C_DRAG, "sail drag")
    arrow(a, (0, 0), (q["drive"] * sc, 0.0), C_DRIVE, "drive", lw=2.6)
    arrow(a, (0, 0), (0.0, -q["lateral"] * sc), COLORS["amber"], "lateral")
    arrow(a, (0, 0), (-q["hull_drag"] * sc, 0.0), C_REF, "hull drag", lw=1.4)
    arrow(a, (0, 0), (0.0, q["keel_lift"] * sc), C_REF, "keel lift", lw=1.4)
    arrow(a, tuple(-1.25 * wind_dir - 0.0), 0.5 * wind_dir, C_FLOW, "apparent wind")
    tw = -np.array([np.cos(th), np.sin(th)])
    arrow(a, tuple(-1.25 * tw), 0.5 * tw, COLORS["blue"], "true wind", lw=1.2)
    a.plot([-0.35, 0.45], [0, 0], color=COLORS["ink"], lw=5, alpha=0.25)
    a.set(aspect="equal", xlim=(-1.5, 1.5), ylim=(-1.3, 1.5), xlabel="along the course (forces / F) [–]", ylabel="across the course [–]",
          title=f"force polygon at θ = {np.rad2deg(th):.0f}°: V = {q['V']:.2f} m/s, β_a = {np.rad2deg(ba):.1f}°")
    save(fig, out, "thrust_fish_sail")

    # --- projectile (Exercise 14.1 formula, our numbers) --------------------------------------------------------------
    print("projectile with constant C_D on a shallow trajectory (our numbers: 40 g, 9 cm², C_D = 0.45, 55 m/s):")
    for deg in (1.0, 3.0, 6.0, 12.0):
        p = ch14.projectile_range_drag(0.040, 1.225, 9.0e-4, 0.45, 55.0, np.deg2rad(deg), numeric=True)
        print(f"  θ₀ = {deg:4.1f}°: formula {p['x_formula']:8.3f} m, numerical {p['x_numeric']:8.3f} m ({p['x_formula'] / p['x_numeric'] - 1:+.2%}), no drag {p['x_vacuum']:8.3f} m")
    lim0 = ch14.projectile_range_drag(0.040, 1.225, 9.0e-4, 1e-9, 55.0, np.deg2rad(3.0), numeric=True)
    print(f"  C_D → 0: formula {lim0['x_formula']:.6f} m → vacuum range {lim0['x_vacuum']:.6f} m")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
