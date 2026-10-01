"""§11.4 (C03–C05): Rayleigh–Bénard — the Γ sign conventions (slip #10), neutral curves Ra(K) for rigid–rigid, free–free,
rigid–free and the odd mode (Fig. 11.10 analogue), the book's 3 × 3 determinant against Chebyshev (and the printed slip #1),
growth rates σ(K) (exchange of stabilities), the gravest even/odd eigenfunctions (Fig. 11.9) and the planforms
(Figs. 11.11–11.12), spectral convergence in N.

Run: ``.venv/Scripts/python.exe scripts/ch11_benard.py --no-show [--fast]``   Figures -> outputs/ch11/c03_benard_neutral.png,
c04_benard_modes.png.  Our water layer: d = 5 mm, ΔT = 2 K, α = 2.1e-4 1/K, κ = 1.4e-7, ν = 1e-6 m²/s.
"""
from __future__ import annotations

import math

import numpy as np

from ch11_common import COLORS, Timer, finish, parse_args, save, setup

from fluidpy import ch11_instability as ch11


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    gc = ch11.gamma_conventions(2.0, 5e-3)
    sc = ch11.benard_scales(2.1e-4, 2.0, 5e-3, 1.4e-7, 1e-6)
    print(f"water layer d = 5 mm, ΔT = 2 K: Γ (11.21) = {gc['Gamma_11_21']:.0f} K/m, Kundu dT/dz = {gc['dTdz_kundu']:.0f} K/m, "
          f"Γ_met = {gc['Gamma_met']:.0f} K/m; Ra = {sc['Ra']:.1f} (G0), {ch11.rayleigh_number(2.1e-4, 2.0, 5e-3, 1.4e-7, 1e-6, g=9.81):.1f} (g = 9.81); "
          f"Pr = {sc['Pr']:.3f}; w ~ {sc['w_scale']:.2e} m/s, t ~ {sc['t_scale']:.0f} s")
    with Timer("critical points"):
        for name, bc, mode in (("rigid–rigid", ("rigid", "rigid"), "even"), ("free–free", ("free", "free"), "even"),
                               ("rigid–free", ("rigid", "free"), "any"), ("odd, rigid–rigid", ("rigid", "rigid"), "odd")):
            cp = ch11.benard_critical(bc, mode, N=24 if args.fast else 40)
            print(f"  {name:18s} Ra_c = {cp['Ra_c']:.4f} at K_c = {cp['K_c']:.4f} (λ_c = {2 * math.pi / cp['K_c']:.3f} d)")
    ff = ch11.benard_free_free_critical()
    print(f"  free–free exact 27π⁴/4 = {ff['Ra_c']:.6f}, π/√2 = {ff['K_c']:.6f}")
    print("determinant (p. 489) vs Chebyshev, and the printed slip #1:")
    for K in (2.0, 3.1163, 5.0):
        rd, rc, rp = ch11.benard_marginal_Ra_det(K), ch11.benard_marginal_Ra(K), ch11.benard_marginal_Ra_det(K, printed=True)
        print(f"  K = {K}: det {rd:.6f}, Chebyshev {rc:.6f} (diff {abs(rd - rc) / rc:.1e}); printed {rp:.3f} ({(rp / rc - 1) * 100:+.2f} %)")
    print("spectral convergence of Ra(3.1163), rigid–rigid:")
    ref = ch11.benard_marginal_Ra(3.1163, N=60)
    for N in (8, 12, 16, 24, 32, 40):
        print(f"  N = {N:2d}: {ch11.benard_marginal_Ra(3.1163, N=N):.10f} (err {abs(ch11.benard_marginal_Ra(3.1163, N=N) - ref):.1e})")
    Ks = np.linspace(0.8, 9.0, 30 if args.fast else 70)
    with Timer("neutral curves"):
        rr = ch11.benard_neutral_curve(Ks)
        rf = ch11.benard_neutral_curve(Ks, ("rigid", "free"), "any")
        od = ch11.benard_neutral_curve(Ks[Ks > 3.0], mode="odd")
    fig, ax = plt.subplots(1, 2, figsize=(14, 4.8))
    ax[0].plot(rr, Ks, color=COLORS["accent"], label="rigid–rigid")
    ax[0].plot(ch11.benard_free_free_Ra(Ks), Ks, color=COLORS["teal"], label="free–free (11.44)")
    ax[0].plot(rf, Ks, color=COLORS["orange"], label="rigid–free")
    ax[0].plot(od, Ks[Ks > 3.0], color=COLORS["muted"], ls="--", label="odd mode (two rows)")
    ax[0].set(xscale="log", xlabel="Ra", ylabel="K = |K|d", title="marginal curves σ = 0 (cf. Fig. 11.10): unstable to the right")
    ax[0].legend(fontsize=8)
    K = 3.0
    Ras = np.linspace(500, 4000, 60)
    sig = [ch11.benard_growth_rate(K, R, 7.0, N=24 if args.fast else 32) for R in Ras]
    sff = [ch11.benard_free_free_sigma(K, R, 7.0)[0] for R in Ras]
    ax[1].plot(Ras, sig, color=COLORS["accent"], label="rigid–rigid (Chebyshev)")
    ax[1].plot(Ras, sff, color=COLORS["teal"], label="free–free (D12)")
    ax[1].axhline(0, color=COLORS["muted"], lw=1)
    ax[1].set(xlabel="Ra", ylabel="σ d²/κ", title=f"leading growth rate at K = {K}, Pr = 7")
    ax[1].legend(fontsize=8)
    save(fig, out, "c03_benard_neutral")
    allsig = ch11.benard_growth_rate(3.0, 3000.0, 1.0, all=True, return_complex=True)
    print(f"exchange of stabilities: max |σ_i| over {len(allsig)} converged modes = {np.max(np.abs(allsig.imag)):.2e}")
    print(f"σ₁ at the marginal Ra(3) = {ch11.benard_marginal_Ra(3.0):.4f}: {ch11.benard_growth_rate(3.0, ch11.benard_marginal_Ra(3.0), 7.0):.2e}")
    fig, ax = plt.subplots(2, 3, figsize=(15, 7))
    for j, (mode, Kc) in enumerate((("even", 3.1163), ("odd", 5.3647))):
        e = ch11.benard_eigenfunction(Kc, mode=mode)
        ax[j, 0].plot(e["W_profile"], e["z_nodes"], color=COLORS["accent"], label="W(z)")
        ax[j, 0].plot(e["T_profile"] / np.max(np.abs(e["T_profile"])), e["z_nodes"], color=COLORS["rose"], label="T̂(z) (scaled)")
        ax[j, 0].set(title=f"gravest {mode} mode, Ra = {e['Ra']:.1f}", xlabel="amplitude", ylabel="z/d")
        ax[j, 0].legend(fontsize=8)
        cs = ax[j, 1].contourf(e["x"], e["z"], e["T"], 21, cmap="RdBu_r")
        ax[j, 1].contour(e["x"], e["z"], e["psi"], 10, colors="k", linewidths=0.7)
        ax[j, 1].set(title="rolls: ψ (lines) over T′ (colour)", xlabel="x/d", ylabel="z/d")
        fig.colorbar(cs, ax=ax[j, 1])
    X, Y = np.meshgrid(np.linspace(0, 4 * math.pi / 3.117, 120), np.linspace(0, 4 * math.pi / 3.117, 120))
    for a, kind in zip((ax[0, 2], ax[1, 2]), ("rolls", "hexagons")):
        f = ch11.planform(X, Y, 3.117, kind)
        a.contourf(X, Y, f, 21, cmap="RdBu_r")
        a.set(title=f"planform: {kind} (∇_H²f = −K²f)", xlabel="x/d", ylabel="y/d")
        a.set_aspect("equal")
    save(fig, out, "c04_benard_modes")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
