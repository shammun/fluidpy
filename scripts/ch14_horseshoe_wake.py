"""Chapter 14, §14.6 — horseshoe vortices, the trailing sheet and its roll-up into tip vortices.

* one horseshoe: centre downwash Γ/(πs), 2× far behind, 0 far ahead, upwash outboard of the tips (formation flight);
* a staircase of horseshoes converging to the downwash integral (14.13);
* downwash ahead of and behind the lifting line (Fig. 14.22); ground effect by an image horseshoe;
* roll-up of the sheet as point vortices (qualitative past the first turn) with its conserved impulse; tip-vortex
  numbers for our aircraft; core size and wake energy.

Run: ``.venv/Scripts/python.exe scripts/ch14_horseshoe_wake.py --no-show [--fast]``   Figures -> outputs/ch14/wake_*.png
"""
from __future__ import annotations

import numpy as np
from ch14_common import C_CIRC, C_DRAG, C_FLOW, C_LIFT, C_REF, COLORS, Timer, finish, parse_args, save, setup

from fluidpy import ch14_aerodynamics as ch14


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch14.illustrative_inputs()
    s, U, rho, W = I["s"], I["U"], I["rho"], I["W"]
    G = 12.0
    w0 = G / (np.pi * s)
    u0 = ch14.horseshoe_velocity(0.0, 0.0, 0.0, G, s)
    print(f"one horseshoe vortex, Γ = {G} m²/s, span {s} m: downwash at the centre {-u0[2]:.8f} m/s = Γ/(πs) = {w0:.8f}; "
          f"far behind (x = 1000 s) {ch14.horseshoe_downwash(1000 * s, 0.0, G, s) / w0:.6f} × that, far ahead {ch14.horseshoe_downwash(-1000 * s, 0.0, G, s) / w0:.2e} ×")
    fin = ch14.horseshoe_velocity(0.3 * s, 0.2 * s, 0.1 * s, G, s, x_end=1.0e6 * s)
    inf = ch14.horseshoe_velocity(0.3 * s, 0.2 * s, 0.1 * s, G, s)
    print(f"  closed half-line legs against three finite segments with legs of 10⁶ spans: max difference {max(abs(a - b) for a, b in zip(fin, inf)):.1e} m/s")
    for yo in (0.25, 0.45, 0.6, 1.0, 2.0):
        print(f"  on the lifting line at y = {yo:4.2f} s: downwash {ch14.horseshoe_downwash(0.0, yo * s, G, s) / w0:+9.4f} × centre value" + ("  (upwash: outboard of the tip)" if yo > 0.5 else ""))
    hw = ch14.horseshoe_wing(G, s, rho, U)
    print(f"  as a 'wing': L = {hw['L']:.1f} N, D_i with the centre downwash everywhere = {hw['Di_centre']:.3f} N = {hw['ratio']:.2f} × the elliptic value "
          f"{hw['Di_elliptic_same_lift']:.3f} N — an artefact ({hw['note']})")

    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
    # --- staircase → integral -----------------------------------------------------------------------------------------
    G1 = 12.0
    load = lambda y: ch14.elliptic_circulation(y, s, G1)                 # noqa: E731
    wex = ch14.elliptic_downwash(G1, s)
    print(f"staircase of N side-by-side horseshoes for the elliptic loading (downwash at mid-span over the exact Γ₁/2s = {wex:.6f} m/s):")
    for N, col in zip((1, 3, 7, 21, 81, 321), (C_REF, C_FLOW, C_CIRC, C_LIFT, C_DRAG, COLORS["ink"])):
        ye, Gs = ch14.staircase_from_loading(load, s, N)
        ym = 0.5 * (ye[:-1] + ye[1:])
        wst = ch14.staircase_downwash(ym, ye, Gs)
        print(f"  N = {N:3d}: w(0)/w_exact = {wst[N // 2] / wex:.6f}")
        if N <= 21:
            ax[0].step(ye, np.append(Gs, Gs[-1]), where="post", color=col, lw=1.2, label=f"N = {N}")
            ax[1].plot(ym, wst / wex, "o-" if N < 21 else "-", ms=3, color=col, label=f"N = {N}")
    yy = np.linspace(-0.5 * s, 0.5 * s, 301)
    ax[0].plot(yy, load(yy), color=COLORS["ink"], ls="--", lw=1, label="Γ(y) elliptic")
    ax[0].set(xlabel="y [m]", ylabel="Γ [m²/s]", title="bound circulation as a staircase of horseshoes")
    ax[0].legend(fontsize=7)
    ax[1].axhline(1.0, color=COLORS["ink"], ls="--", lw=1)
    ax[1].set(xlabel="y [m]", ylabel=r"$w / (\Gamma_1/2s)$ [–]", title="downwash at panel centres → the integral (14.13)", ylim=(0, 3))
    ax[1].legend(fontsize=7)

    # --- Fig. 14.22: downwash along the flight direction -------------------------------------------------------------
    xs = np.linspace(-3.0, 3.0, 121) * s
    wx = ch14.downwash_along_x(xs, [G1], s)
    ax[2].plot(xs / s, wx / wex, color=C_FLOW, label="trailing sheet only")
    xb = xs[np.abs(xs) > 0.12 * s]
    ax[2].plot(xb / s, ch14.downwash_along_x(xb, [G1], s, include_bound=True) / wex, color=C_REF, ls=":", label="sheet + bound vortex (singular at the line)")
    ax[2].axhline(2.0, color=C_REF, lw=0.8)
    ax[2].axhline(0.0, color=C_REF, lw=0.8)
    ax[2].set(xlabel="x / s (downstream +) [–]", ylabel=r"downwash / $w(0)$ [–]", title="mid-span downwash: 0 far ahead, w at the wing, 2w far behind", ylim=(-1.5, 3.5))
    ax[2].legend(fontsize=8)
    save(fig, out, "wake_horseshoe_downwash")
    print("downwash of the trailing sheet on the centre line, over its value at the wing: " + ", ".join(
        f"x = {q:+.1f}s: {ch14.downwash_along_x(q * s, [G1], s) / wex:.4f}" for q in (-20.0, -1.0, -0.2, 0.0, 0.2, 1.0, 20.0)))

    # --- ground effect and the upwash map --------------------------------------------------------------------------
    print("ground effect (image horseshoe of opposite strength): centre downwash over its free-air value")
    for h in (0.05, 0.1, 0.25, 0.5, 1.0, 5.0):
        img = ch14.horseshoe_velocity(0.0, 0.0, 2.0 * h * s, -G, s)[2]
        print(f"  h/s = {h:4.2f}: two horseshoes {(w0 - img) / w0:.6f}, closed form 1 − 1/(1 + (4h/s)²) = {ch14.ground_effect_factor(h * s, s):.6f}")
    X, Y = np.meshgrid(np.linspace(-1.0, 3.0, 161) * s, np.linspace(-1.5, 1.5, 151) * s, indexing="xy")
    Wd = ch14.horseshoe_downwash(X, Y, G, s, z=0.02 * s) / w0
    fig, ax = plt.subplots(1, 2, figsize=(14, 4.8))
    im = ax[0].pcolormesh(X / s, Y / s, np.clip(Wd, -3, 3), cmap="RdBu", shading="auto", vmin=-3, vmax=3)
    fig.colorbar(im, ax=ax[0], label="downwash / centre value [–] (blue = sinking air, red = rising)")
    ax[0].plot([0, 0], [-0.5, 0.5], color=COLORS["ink"], lw=3)
    ax[0].plot([0, 3], [0.5, 0.5], color=COLORS["ink"], lw=1)
    ax[0].plot([0, 3], [-0.5, -0.5], color=COLORS["ink"], lw=1)
    ax[0].set(xlabel="x / s [–]", ylabel="y / s [–]", title="horseshoe vortex: sinking air between the legs, rising outboard")

    # --- roll-up ------------------------------------------------------------------------------------------------------
    with Timer("wake roll-up"):
        R = ch14.wake_rollup([G1], s, fast=args.fast, n_t=41)
    for k, col in zip((0, 10, 20, 40), (C_REF, C_FLOW, C_CIRC, C_DRAG)):
        ax[1].plot(R["y"][k] / s, R["z"][k] / s, ".", ms=3, color=col, label=f"t Γ₀/s² = {R['t'][k] * R['Gamma0'] / s ** 2:.2f}")
    ax[1].set(aspect="equal", xlabel="y / s [–]", ylabel="z / s [–]", title="trailing sheet rolling up (rear view) — qualitative past the first turn")
    ax[1].legend(fontsize=8)
    save(fig, out, "wake_upwash_rollup")
    print(f"roll-up of the elliptic loading ({R['y'].shape[1]} vortices, core {R['eps'] / s:.3f} s): impulse ΣΓy varies by "
          f"{np.ptp(R['impulse']) / abs(R['impulse'][0]):.1e}; starboard centroid y = {R['centroid_y'][0] / s:.5f} s → {R['centroid_y'][-1] / s:.5f} s "
          f"(theory b′/2 = π/8 s = {R['spacing_theory'] / (2 * s):.5f} s); mean descent {-R['centroid_z'][-1] / R['t'][-1]:.5f} m/s "
          f"(pair formula Γ₀/2πb′ = {R['descent_speed_theory']:.5f}); energy drift {abs(R['energy'][-1] - R['energy'][0]) / abs(R['energy'][0]):.1e}")
    tv = ch14.tip_vortex_estimate(W, rho, U, s)
    print(f"our aircraft (W = {W} N, s = {s} m, U = {U} m/s, sea level): tip-vortex circulation {tv['Gamma0']:.2f} m²/s, spacing πs/4 = {tv['spacing']:.2f} m, "
          f"sink rate {tv['descent_speed']:.3f} m/s, downwash at the wing {tv['downwash']:.3f} m/s ({np.rad2deg(tv['downwash'] / U):.2f}°), induced drag {tv['induced_drag']:.1f} N")
    print("  wake kinetic energy per metre for Rankine cores of radius σ (ours): " + ", ".join(
        f"σ = {q:.2f} m: {ch14.vortex_pair_energy_per_length(tv['Gamma0'], tv['spacing'], q, rho):.1f} J/m" for q in (0.05, 0.2, 0.8))
          + f"  — a larger core holds less; (14.25) gives {tv['induced_drag']:.1f} N")
    sig_eq = tv["spacing"] * np.exp(0.25 - np.pi ** 2 / 4.0)
    print(f"  the Rankine core that holds exactly the induced drag: σ = b′exp(¼ − π²/4) = {sig_eq:.3f} m = {sig_eq / s:.4f} s → "
          f"{ch14.vortex_pair_energy_per_length(tv['Gamma0'], tv['spacing'], sig_eq, rho):.2f} J/m")
    rm = ch14.wake_vortex_rolling_moment(rho, U, I["K"], s, I["A"] / s, I["Gamma_v"], I["ell"])
    print(f"wake encounter (our numbers: Γ′ = {I['Gamma_v']} m²/s, ℓ = {I['ell']} m): rolling moment {rm['M']:.0f} N m (quadrature {rm['M_quad']:.0f}), "
          f"{rm['ratio']:.3f} of the ℓ → 0 limit")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
