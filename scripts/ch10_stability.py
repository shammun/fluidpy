"""§10.2 (C04, C05) and §10.4 (C10): von Neumann amplification factors (10.24)–(10.26), Noye's FTCS region (10.27), upwind
and the CFL condition (10.29)–(10.30), MacCormack = Lax–Wendroff on a square pulse (diffusion vs dispersion).

Run: ``.venv/Scripts/python.exe scripts/ch10_stability.py --no-show``   Figures -> outputs/ch10/c04_von_neumann.png,
c05_upwind_cfl.png.
"""
from __future__ import annotations

import numpy as np

from ch10_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch10_computational_fluid_dynamics as ch10


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    circ = np.linspace(0, 2 * np.pi, 200)
    ax[0].plot(np.cos(circ), np.sin(circ), color=COLORS["muted"], lw=1)
    cases = [("ftcs", 0.1, 0.2, COLORS["orange"]), ("ftcs", 0.0, 0.51, COLORS["rose"]), ("btcs", 0.1, 0.2, COLORS["blue"]),
             ("upwind", 0.4, 0.0, COLORS["teal"]), ("lax_wendroff", 0.4, 0.0, COLORS["accent"])]
    for s, a, b, c in cases:
        cv = ch10.amplification_curve(s, a, b)
        v = ch10.stability_verdict(s, a, b)
        ax[0].plot(cv["G"].real, cv["G"].imag, color=c, label=f"{s} α={a}, β={b}")
        ax[1].plot(cv["theta"], cv["modulus"], color=c, label=f"{s}: max|G| = {v['Gmax']:.3f}")
        print(f"{s:12s} α = {a:4.2f} β = {b:4.2f}: max|G| = {v['Gmax']:.4f} at θ = {v['theta_worst']:.3f} rad — {v['reason']}")
    ax[0].set_aspect("equal")
    ax[0].set_xlabel("Re G")
    ax[0].set_ylabel("Im G")
    ax[0].set_title("G(θ) for θ ∈ [0, π] and the unit circle")
    ax[0].legend(fontsize=7)
    ax[1].axhline(1, color=COLORS["muted"], lw=1)
    ax[1].set_xlabel("θ = kπΔx  [rad]")
    ax[1].set_ylabel("|G|")
    ax[1].set_title("amplification per step")
    ax[1].legend(fontsize=7)
    A, B = np.meshgrid(np.linspace(0, 0.6, 121), np.linspace(0, 0.7, 141))
    S = np.vectorize(lambda a, b: ch10.max_amplification("ftcs", a, b, 361))(A, B)
    cs = ax[2].contourf(A, B, np.log10(S), levels=np.linspace(0, 0.5, 11), cmap="Reds", extend="max")
    ax[2].contour(A, B, S, levels=[1 + 1e-9], colors=[COLORS["ink"]])
    a = np.linspace(0, 0.5, 100)
    ax[2].plot(a, 2 * a ** 2, color=COLORS["teal"], lw=2, label=r"$4\alpha^2=2\beta$")
    ax[2].axhline(0.5, color=COLORS["blue"], lw=2, label=r"$2\beta=1$")
    ax[2].set_xlabel("α = uΔt/(2Δx)")
    ax[2].set_ylabel("β = DΔt/Δx²")
    ax[2].set_title("FTCS: log10 max|G| (white = stable, Noye (10.27))")
    ax[2].legend(fontsize=8)
    fig.colorbar(cs, ax=ax[2])
    # brute force vs closed form
    rng = np.random.default_rng(0)
    ab = rng.uniform(0, 0.8, (2000, 2))
    dis = sum(ch10.ftcs_stable(a, b) != ch10.is_von_neumann_stable("ftcs", a, b) for a, b in ab
              if abs(4 * a ** 2 - 2 * b) > 5e-3 and abs(2 * b - 1) > 5e-3)
    print(f"Noye closed form vs brute-force max|G| on 2000 random (α, β): {dis} disagreements")
    save(fig, out, "c04_von_neumann")

    fig2, ax2 = plt.subplots(1, 2, figsize=(12, 4.2))
    for C, ls in ((0.5, "-"), (0.8, "--")):
        for s, c in (("upwind", COLORS["teal"]), ("maccormack", COLORS["accent"])):
            r = ch10.advect_periodic("square", C=C, n_cells=100, n_rev=1.0, scheme=s)
            ax2[0].plot(r["x"], r["T"], color=c, ls=ls, label=f"{s}, C = {r['C']:.2f}")
            print(f"square pulse one revolution, {s:10s} C = {r['C']:.2f}: amplitude ratio {r['amplitude_ratio']:.3f}, "
                  f"rms error {r['rms_error']:.3f}")
    ax2[0].plot(r["x"], r["exact"], color=COLORS["muted"], lw=1, label="exact")
    ax2[0].set_xlabel("x (periodic unit domain)")
    ax2[0].set_ylabel("T")
    ax2[0].set_title("upwind smears (diffusion), MacCormack ripples (dispersion)")
    ax2[0].legend(fontsize=7)
    for C in (0.9, 1.0, 1.1):
        r = ch10.advect_periodic("gauss", C=C, n_cells=100, n_rev=20, scheme="upwind")
        print(f"upwind gauss C = {C}: blew up {r['blew_up']} after {r['steps']} steps; amplitude {r['amplitude_ratio']:.3g}")
    r = ch10.advect_periodic("gauss", C=0.5, n_cells=100, n_rev=20, scheme="ftcs")
    print(f"FTCS pure convection C = 0.5: blew up {r['blew_up']} (step {r['steps']}) — (10.27) is never satisfied when D = 0")
    th = np.linspace(1e-3, np.pi, 300)
    for s, c in (("upwind", COLORS["teal"]), ("lax_wendroff", COLORS["accent"])):
        ax2[1].plot(th, ch10.amplification_modulus(th, 0.4, 0.0, s), color=c, label=f"|G| {s}, C = 0.8")
        ax2[1].plot(th, ch10.phase_error(s, 0.8, th), color=c, ls="--", label=f"relative phase speed {s}")
    ax2[1].axhline(1, color=COLORS["muted"], lw=1)
    ax2[1].set_xlabel("θ  [rad]")
    ax2[1].set_title("amplitude and phase per step (C = 0.8)")
    ax2[1].legend(fontsize=7)
    print(f"numerical diffusivity of upwind, u = 1 m/s, Δx = 0.01 m, C = 0.5: "
          f"{ch10.numerical_diffusivity(1.0, 0.01, scheme='upwind', C=0.5):.4g} m²/s")
    c = np.sqrt(9.81 * 4000.0)
    print(f"climate CFL: external gravity wave √(gH) = {c:.1f} m/s (H = 4000 m): Δt ≤ {ch10.cfl_time_step(c, 100e3):.1f} s "
          f"at Δx = 100 km, {ch10.cfl_time_step(c, 25e3):.1f} s at 25 km")
    save(fig2, out, "c05_upwind_cfl")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
