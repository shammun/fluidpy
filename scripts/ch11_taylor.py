"""§11.6 (C07): Rayleigh's ring interchange and circulation criterion, the Taylor number (11.52), the exact narrow-gap critical
Ta_c(μ) against (11.54), the Galerkin route of Exercise 11.9 (and slip #7), the stability boundary in the (Ω₂, Ω₁) plane for
our radius ratio R₂/R₁ = 1.05 (Fig. 11.17 analogue) and the Taylor vortices.

Run: ``.venv/Scripts/python.exe scripts/ch11_taylor.py --no-show [--fast]``   Figure -> outputs/ch11/c07_taylor.png.
"""
from __future__ import annotations

import numpy as np

from ch11_common import COLORS, Timer, finish, parse_args, save, setup

from fluidpy import ch11_instability as ch11
from fluidpy.core import laminar as LAM


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    r = ch11.ring_interchange_energy(4.0, 2.0, 1.0, 2.0)
    print(f"ring interchange Γ₁ = 4, Γ₂ = 2, r = 1, 2: E_i = {r['E_i']:.5f}, E_f = {r['E_f']:.5f}, ΔE = {r['dE']:.5f} m²/s²")
    R1, R2 = 0.1, 0.105
    for O2 in (0.0, 0.4, 0.46, 0.5):
        Rr = np.linspace(R1, R2, 200)
        rc = ch11.rayleigh_circulation_criterion(Rr, LAM.circular_couette(Rr, R1, R2, 0.5, O2))
        tn = ch11.taylor_number(0.5, O2, R1, R2, 1e-6)
        print(f"Ω₁ = 0.5, Ω₂ = {O2}: Rayleigh-unstable = {rc['unstable']}, Ta = {tn['Ta']:.1f}, rayleigh_stable = "
              f"{tn['rayleigh_stable']} (μ_R = {ch11.couette_rayleigh_line(R1, R2):.4f})")
    print(f"narrow-gap inner-only Ta = {ch11.taylor_number_narrow_inner(0.5, R1, R2 - R1, 1e-6):.1f}")
    mus = np.round(np.arange(-0.5, 1.0001, 0.25 if args.fast else 0.1), 6)
    with Timer("Ta_c(μ)"):
        tab = ch11.taylor_critical_table(mus, N=24 if args.fast else 40, cache=not args.fast)
    for m, T, k, a, e in zip(tab["mu"], tab["Ta_c"], tab["k_c"], tab["approx"], tab["rel_error"]):
        print(f"  μ = {m:5.2f}: Ta_c = {T:9.3f} at k_c = {k:.4f};  (11.54) {a:9.2f} ({e * 100:+.2f} %)")
    print("Galerkin route (Ex. 11.9) at μ = 0, k = 3.1266 vs Chebyshev", f"{ch11.taylor_marginal_Ta(3.1266, 0.0):.4f}:")
    for M in (2, 3, 4, 6, 8):
        print(f"  M = {M}: {ch11.taylor_galerkin_Ta(3.1266, 0.0, M):.4f}   printed (11.93) (slip #7): "
              f"{ch11.taylor_galerkin_Ta(3.1266, 0.0, M, printed=True)}")
    print(f"growth rate at Ta = 3500, μ = 0, k = 3.12: σ = {ch11.taylor_growth_rate(3.12, 3500.0, 0.0):.4f}")
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
    ax[0].plot(tab["mu"], tab["Ta_c"], "o-", color=COLORS["accent"], label="exact narrow gap (11.51)")
    mm = np.linspace(-0.5, 1.0, 100)
    ax[0].plot(mm, [ch11.taylor_critical_approx(m) for m in mm], color=COLORS["muted"], ls="--", label="(11.54) 1708/[½(1 + μ)]")
    ax[0].set(xlabel="μ = Ω₂/Ω₁", ylabel="Ta_c", title="critical Taylor number")
    ax[0].legend(fontsize=8)
    with Timer("stability boundary"):
        sb = ch11.taylor_stability_boundary(1.05, mus=np.linspace(-1.0, 0.86, 10 if args.fast else 25))
    ax[1].plot(sb["x_outer"], sb["y_inner"], color=COLORS["accent"], label="narrow-gap marginal curve")
    ax[1].plot(sb["rayleigh_x"], sb["rayleigh_y"], color=COLORS["muted"], ls="--", label="Rayleigh line Ω₁/Ω₂ = R₂²/R₁²")
    ax[1].set(xlabel="Ω₂R₂²/ν", ylabel="Ω₁R₂²/ν", title="R₂/R₁ = 1.05 (ours): unstable above (cf. Fig. 11.17)",
              ylim=(0, 1.2 * np.max(sb["y_inner"])))
    ax[1].legend(fontsize=8)
    e = ch11.taylor_eigenfunction(tab["k_c"][list(tab["mu"]).index(0.0)] if 0.0 in list(tab["mu"]) else 3.1266, 0.0)
    ax[2].contourf(e["x"], e["z"], e["u_phi"], 21, cmap="RdBu_r")
    ax[2].contour(e["x"], e["z"], e["psi"], 10, colors="k", linewidths=0.7)
    ax[2].set(xlabel="x = (R − R₁)/d", ylabel="z/d", title=f"Taylor vortices μ = 0 (Ta = {e['Ta']:.1f}): ψ over u_φ")
    ax[2].set_aspect("equal")
    save(fig, out, "c07_taylor")
    print(f"boundary: {len(sb['mu'])} points, Ω₁R₂²/ν from {np.min(sb['y_inner']):.0f} to {np.max(sb['y_inner']):.0f}")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
