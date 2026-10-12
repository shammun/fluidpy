"""Chapter 14, §14.3 — the starting vortex and Kelvin's theorem (Figs. 14.9–14.10), with a lumped-vortex model (ours).

* an impulsive start: the bound circulation grows toward the Kutta value πUc sin α while the shed circulation grows
  equal and opposite — their sum stays zero at every step (Kelvin);
* the model's circulatory lift ρUΓ_bound / steady lift against R. T. Jones' fit of Wagner's function — a *lift* ratio
  (a cited approximation), drawn on its own lift panel; the exact circulation starts from zero, the exact circulatory
  lift from one half, the model's bound circulation from 1/3;
* first-order convergence in the time step; rigid against free wake; a step up in incidence and a return to zero lift.

Run: ``.venv/Scripts/python.exe scripts/ch14_starting_vortex.py --no-show``   Figures -> outputs/ch14/starting_vortex_*.png
"""
from __future__ import annotations

import numpy as np
from ch14_common import C_CIRC, C_DRAG, C_LIFT, C_REF, COLORS, Timer, finish, parse_args, save, setup

from fluidpy import ch14_aerodynamics as ch14


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    U, c, alpha = 1.0, 1.0, np.deg2rad(6.0)
    Gs = np.pi * U * c * np.sin(alpha)
    print(f"flat plate, chord {c} m, U = {U} m/s, α = 6°: steady (Kutta) circulation πUc sin α = {Gs:.6f} m²/s "
          f"(= (14.11) with a = b = c/4: {ch14.kutta_circulation(U, c / 4, alpha, 0.0):.6f})")
    with Timer("lumped-vortex march, 400 steps, rigid + free wake"):
        r = ch14.impulsive_start_lumped_vortex(U, c, alpha, 20.0, 0.05, history=True)
        rf = ch14.impulsive_start_lumped_vortex(U, c, alpha, 20.0, 0.05, wake="free")
    print(f"  Kelvin: max |Γ_bound + ΣΓ_wake| over the run = {np.abs(r['total']).max():.1e} m²/s (rigid), {np.abs(rf['total']).max():.1e} (free)")
    # Wagner's function is a LIFT ratio (circulatory lift / steady lift), so it is compared with the model's circulatory
    # lift ρUΓ_bound over the steady lift ρUΓ_steady — not with a circulation.  (For a steady stream the two model ratios
    # are the same number; the exact circulation is a different curve, which starts from zero by Kelvin's theorem.)
    rho = 1.0                                                   # the density the march was run with (its default)
    lift_ratio = r["lift_circulatory"] / (rho * U * r["Gamma_steady"])       # model: circulatory lift / steady lift
    lift_ratio_free = rf["lift_circulatory"] / (rho * U * rf["Gamma_steady"])
    phi = ch14.wagner_function(2 * r["s_chords"])               # Jones' fit, argument in semi-chords
    print("  model, bound circulation / steady circulation and circulatory lift / steady lift; Wagner's function is a lift ratio")
    print("  chords travelled | Γ_bound/Γ_steady (rigid) | ρUΓ_bound/steady lift (rigid) | (free wake) | Wagner's lift function, Jones' fit φ(2 chords)")
    for s in (0.05, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0):
        k = int(np.argmin(np.abs(r["s_chords"] - s)))
        print(f"    {r['s_chords'][k]:5.2f}          | {r['ratio'][k]:.4f}                   | {lift_ratio[k]:.4f}                        | {lift_ratio_free[k]:.4f}      | {phi[k]:.4f}")
    k90 = int(np.argmax(r["ratio"] >= 0.9))
    k90w = int(np.argmax(phi >= 0.9))
    late = r["s_chords"] >= 2.0
    print(f"  model: 90 % of the steady value after {r['s_chords'][k90]:.1f} chords; first step of Γ_bound/Γ_steady {r['ratio'][0]:.4f} (model limit 1/3)")
    print(f"  Jones' fit of Wagner's lift function: starts at {float(ch14.wagner_function(0.0)):.3f}, 90 % after {r['s_chords'][k90w]:.1f} chords; "
          f"max |model lift ratio − fit| over 2–20 chords = {np.abs(lift_ratio - phi)[late].max():.4f}")
    print("  exact thin-airfoil start-up (computed in the derivation review, reports/ch14_review.md, not by this script): the circulation "
          "starts at 0 (Kelvin: nothing shed yet) and reaches 90 % after 7.1 chords; the circulatory lift starts at 1/2 and reaches 90 % after 6.4 chords")
    print("  time-step study (bound/steady at 2 chords):")
    vals = []
    for dt in (0.2, 0.1, 0.05, 0.025, 0.0125):
        q = ch14.impulsive_start_lumped_vortex(U, c, alpha, 2.0, dt)
        vals.append(q["ratio"][-1])
        print(f"    dt = {dt:6.4f} s: {vals[-1]:.6f}")
    d = np.diff(vals)
    print(f"    ratios of successive differences {np.round(d[:-1] / d[1:], 2)} (2 = first order)")

    fig, ax = plt.subplots(1, 4, figsize=(23, 4.6))
    # --- panel 0: circulation.  The only reference here is the steady Kutta value (Wagner's curve is NOT a circulation) --
    ax[0].plot(r["s_chords"], r["Gamma_bound_cw"] / Gs, color=C_LIFT, label="bound (foil), clockwise")
    ax[0].plot(r["s_chords"], r["Gamma_wake_cw"] / Gs, color=C_DRAG, label="shed (starting vortex)")
    ax[0].plot(r["s_chords"], r["total"] / Gs, color=C_CIRC, lw=2.5, label="sum (Kelvin: zero)")
    ax[0].axhline(1.0, color=C_REF, lw=1.0, ls=":", label=r"steady Kutta value $\pi U c\sin\alpha$")
    ax[0].set(xlabel="distance travelled [chords]", ylabel=r"circulation / $\pi U c\sin\alpha$ [–]", title="impulsive start: bound + shed = 0")
    ax[0].legend(fontsize=8, loc="center right")
    # --- panel 1: lift.  Wagner's function belongs on this axis: circulatory lift / steady lift -------------------------
    ax[1].plot(r["s_chords"], lift_ratio, color=C_LIFT, label=r"model: $\rho U\Gamma_{bound}$ / steady lift $\rho U\Gamma_{steady}$")
    ax[1].plot(r["s_chords"], phi, color=C_REF, ls="--", label="Wagner's function — circulatory lift / steady lift (Jones' fit)")
    ax[1].axhline(1.0, color=C_REF, lw=0.8, ls=":")
    ax[1].set(xlabel="distance travelled [chords]", ylabel="circulatory lift / steady lift [–]", ylim=(0.0, 1.05),
              title="lift after the start: model against Wagner")
    ax[1].legend(fontsize=8, loc="lower right")
    for kk, col in zip((19, 99, 399), (C_REF, C_CIRC, C_DRAG)):
        x, y = r["wake_x_hist"][kk], r["wake_y_hist"][kk]
        ax[2].plot(x, y, ".", ms=3, color=col, label=f"after {r['s_chords'][kk]:.0f} chords")
    ax[2].plot([0, c], [0, 0], color=COLORS["ink"], lw=3)
    ax[2].set(xlabel="x [m] (plate frame)", ylabel="y [m]", title="shed vortices (rigid wake), plate at 0…c")
    ax[2].legend(fontsize=8)

    # --- events: a step up in incidence, then back to zero incidence (every change of lift sheds a vortex) ----------
    a_of_t = lambda t: alpha if t < 8.0 else (2.0 * alpha if t < 16.0 else 0.0)   # noqa: E731
    ev = ch14.impulsive_start_lumped_vortex(U, c, a_of_t, 30.0, 0.05)
    ax[3].plot(ev["t"], ev["Gamma_bound_cw"] / Gs, color=C_LIFT, label="bound")
    ax[3].plot(ev["t"], ev["Gamma_wake_cw"] / Gs, color=C_DRAG, label="shed (all vortices)")
    ax[3].plot(ev["t"], ev["total"] / Gs, color=C_CIRC, lw=2.5, label="sum")
    ax[3].plot(ev["t"], ev["Gamma_steady"] / Gs, color=C_REF, ls="--", label="steady value for the current U, α")
    ax[3].set(xlabel="time t [s]", ylabel=r"circulation / $\pi U c\sin\alpha_0$ [–]", title="start, α doubled at 8 s, α = 0 from 16 s")
    ax[3].legend(fontsize=8)
    save(fig, out, "starting_vortex_kelvin")
    k8, k16 = int(np.argmin(np.abs(ev["t"] - 8.0))), int(np.argmin(np.abs(ev["t"] - 16.0)))
    print(f"events (bound/steady₀): just before α doubles {ev['Gamma_bound_cw'][k8 - 1] / Gs:.3f}, just before α returns to 0 "
          f"{ev['Gamma_bound_cw'][k16 - 1] / Gs:.3f}, at the end {ev['Gamma_bound_cw'][-1] / Gs:+.4f}; vortex shed at each event "
          f"(clockwise +): start {ev['wake_Gamma_cw'][0] / Gs:+.3f}, α doubled {ev['wake_Gamma_cw'][k8] / Gs:+.3f}, α → 0 "
          f"{ev['wake_Gamma_cw'][k16] / Gs:+.3f} — opposite to the change of bound circulation each time; max |sum| = {np.abs(ev['total']).max():.1e}")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
