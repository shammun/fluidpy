"""Chapter 14, §14.1–§14.2 and §14.7 — level flight, the drag breakdown, the polar and the best glide.

* is the flow incompressible? Mach number and ½Ma²; level flight: the C_L the wing must fly at; stall speed;
* induced drag coefficient (14.26) against aspect ratio; the polar C_D = C_D0 + C_L²/(πΛ e_span) and its best L/D;
* drag against speed: profile ∝ U², induced ∝ 1/U² (14.25) — the induced share is ½ at the best glide and small at
  high speed (trap T9);
* generic measured lift and drag curves with stall (qualitative models, our constants) and the two stall regimes.

Run: ``.venv/Scripts/python.exe scripts/ch14_drag_polar.py --no-show``   Figures -> outputs/ch14/drag_*.png
"""
from __future__ import annotations

import numpy as np
from ch14_common import C_DRAG, C_FLOW, C_LIFT, C_REF, COLORS, finish, parse_args, save, setup

from fluidpy import ch14_aerodynamics as ch14


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch14.illustrative_inputs()
    W, s, A, AR, U, rho, CD0, es, K = I["W"], I["s"], I["A"], I["AR"], I["U"], I["rho"], I["CD0"], I["e_span"], I["K"]
    cc = ch14.compressibility_check(U)
    print(f"our aircraft: W = {W} N, s = {s} m, A = {A} m², Λ = s²/A = {ch14.aspect_ratio(s, A):.3f}, c̄ = {ch14.mean_chord(A, s):.3f} m, U = {U} m/s")
    print(f"  Mach number {cc['Ma']:.4f} → density changes ≈ ½Ma² = {cc['density_change']:.2%}: incompressible treatment ok = {cc['incompressible_ok']}; "
          f"at 100 m/s: Ma = {ch14.compressibility_check(100.0)['Ma']:.3f}, ½Ma² = {ch14.compressibility_check(100.0)['density_change']:.2%}")
    CL = ch14.level_flight(W, rho, U, A)["CL"]
    CD = ch14.drag_polar(CL, CD0=CD0, AR=AR, e_span=es)
    lf = ch14.level_flight(W, rho, U, A, CD)
    br = ch14.drag_breakdown(CL, CD0=CD0, AR=AR, e_span=es, friction_fraction=0.7)
    print(f"  level flight: C_L = W/(½ρU²A) = {CL:.4f}; C_Di = C_L²/(πΛe_span) = {br['induced']:.5f}, profile C_D0 = {CD0} (friction {br['friction']:.4f} + form "
          f"{br['form']:.4f}, split assumed 70/30), total {br['total']:.5f}; drag {lf['D']:.1f} N, power {lf['power'] / 1000:.2f} kW; induced share {br['induced_share']:.1%}")
    print(f"  Reynolds number on the mean chord {ch14.chord_reynolds(U, A / s, I['nu']):.2e} → '{ch14.stall_regime(ch14.chord_reynolds(U, A / s, I['nu']), Re_laminar=I['Re_laminar'], Re_turbulent=I['Re_turbulent'])}' "
          f"(thresholds ours: {I['Re_laminar']:.0e}, {I['Re_turbulent']:.0e})")
    print(f"  stall speed: {ch14.stall_speed(W, rho, A, I['CL_max']):.2f} m/s with C_L,max = {I['CL_max']}; {ch14.stall_speed(W, rho, A, I['CL_max_flap']):.2f} m/s with flaps (C_L,max = {I['CL_max_flap']})")
    bg = ch14.best_glide(CD0, AR, es)
    Us = ch14.speed_min_drag(W, rho, A, AR, CD0, es)
    cl = np.linspace(0.02, 1.6, 4000)
    ld = cl / ch14.drag_polar(cl, CD0=CD0, AR=AR, e_span=es)
    print(f"  best glide: C_L = √(πΛe C_D0) = {bg['CL']:.5f}, (L/D)_max = ½√(πΛe/C_D0) = {bg['LD_max']:.4f} (brute-force scan {ld.max():.4f} at C_L = {cl[np.argmax(ld)]:.4f}), "
          f"at U = {Us:.2f} m/s; induced share there {ch14.drag_share(bg['CL'], CD0=CD0, AR=AR, e_span=es):.3f}")
    print(f"  profile-drag estimate (approximate: plate friction × form factor 1.25, our choice): laminar {ch14.profile_drag_estimate(1.5e6, 0.12, regime='laminar', form_factor=1.25):.5f}, "
          f"turbulent {ch14.profile_drag_estimate(1.5e6, 0.12, regime='turbulent', form_factor=1.25):.5f} at Re_c = 1.5e6")

    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
    # --- polar for several aspect ratios -----------------------------------------------------------------------------
    for L_, col in zip((5.0, AR, 30.0), (C_REF, C_LIFT, C_FLOW)):
        ax[0].plot(ch14.drag_polar(cl, CD0=CD0, AR=L_, e_span=es), cl, color=col, label=f"Λ = {L_:.1f}: (L/D)max = {ch14.best_glide(CD0, L_, es)['LD_max']:.1f}")
        b_ = ch14.best_glide(CD0, L_, es)
        ax[0].plot([0, b_["CD"]], [0, b_["CL"]], color=col, lw=0.8, ls="--")
        ax[0].plot([b_["CD"]], [b_["CL"]], "o", color=col, ms=5)
    ax[0].axvline(CD0, color=C_REF, lw=0.8, ls=":")
    ax[0].set(xlabel=r"$C_D$ [–]", ylabel=r"$C_L$ [–]", title="drag polar: the tangent from the origin is the best glide", xlim=(0, 0.12))
    ax[0].legend(fontsize=8)
    # --- drag against speed -----------------------------------------------------------------------------------------
    Umin = ch14.stall_speed(W, rho, A, I["CL_max"])
    uu = np.linspace(0.9 * Umin, 75.0, 300)
    dv = ch14.drag_vs_speed(W, rho, A, AR, CD0, uu, es)
    ax[1].plot(uu, dv["induced"], color=C_DRAG, label="induced (14.25) ∝ 1/U²")
    ax[1].plot(uu, dv["profile"], color=C_FLOW, label="profile ∝ U²")
    ax[1].plot(uu, dv["total"], color=COLORS["ink"], label="total")
    ax[1].axvline(Us, color=C_LIFT, ls="--", lw=1.2, label=f"best glide {Us:.1f} m/s")
    ax[1].axvspan(uu[0], Umin, color=C_REF, alpha=0.25, label="below stall speed")
    ax[1].set(xlabel="speed U [m/s]", ylabel="drag [N]", title="level flight: you pay for lift when you fly slowly", ylim=(0, 1.3 * dv["total"].max()))
    ax[1].legend(fontsize=8)
    print("  speed [m/s] | C_L   | induced [N] | profile [N] | induced share | L/D")
    for u in (Umin, Us, U, 50.0, 70.0):
        q = ch14.drag_vs_speed(W, rho, A, AR, CD0, u, es)
        print(f"    {u:6.2f}    | {q['CL']:.3f} | {q['induced']:8.1f}    | {q['profile']:8.1f}    | {q['induced_share']:6.1%}        | {q['LD']:.2f}")
    print("  induced drag coefficient (14.26) at C_L = 0.6: " + ", ".join(f"Λ = {q}: {ch14.induced_drag_coefficient(0.6, q):.5f}" for q in (4, 8, 16, 32, 1e9)))
    print(f"  doubling the span at fixed weight and speed: induced drag × {ch14.elliptic_induced_drag(W, rho, U, 2 * s) / ch14.elliptic_induced_drag(W, rho, U, s):.3f} (14.25)")

    # --- generic measured curves ------------------------------------------------------------------------------------
    al = np.deg2rad(np.linspace(-6.0, 24.0, 241))
    kw = dict(K=K, beta=I["beta"], alpha_stall=I["alpha_stall"], CL_post=I["CL_post"], width=I["stall_width"])
    cl_g = ch14.generic_lift_curve(al, **kw)
    cd_g = ch14.generic_drag_curve(al, CD0=CD0, AR=AR, e_span=es, CD_sep=I["CD_sep"], **kw)
    ax[2].plot(np.rad2deg(al), cl_g, color=C_LIFT, label=r"$C_L$ (model)")
    ax[2].plot(np.rad2deg(al), K * (al + I["beta"]), color=C_LIFT, lw=1, ls="--", label=r"$K(\alpha+\beta)$")
    ax[2].plot(np.rad2deg(al), 10.0 * cd_g, color=C_DRAG, label=r"$10\,C_D$ (model)")
    ax[2].axvline(np.rad2deg(I["alpha_stall"]), color=C_REF, lw=0.8, ls=":")
    ax[2].set(xlabel="angle of attack α [deg]", ylabel="coefficient [–]", title="generic lift and drag curves with stall — qualitative model", ylim=(-0.4, 2.2))
    ax[2].legend(fontsize=8)
    save(fig, out, "drag_polar_speed_stall")
    km = int(np.argmax(cl_g))
    print(f"generic curves (qualitative, our constants): linear slope {K:.4f} per rad = {np.deg2rad(K):.4f} per degree, zero lift at {-np.rad2deg(I['beta']):.2f}°, "
          f"leaves the line at {np.rad2deg(I['alpha_stall']):.0f}°, C_L,max = {cl_g[km]:.3f} at {np.rad2deg(al[km]):.1f}°, plateau {I['CL_post']}; "
          f"C_D minimum {cd_g.min():.4f} at {np.rad2deg(al[np.argmin(cd_g)]):.2f}°")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
