"""§11.8–§11.11 (C11, C13, C14): Orr–Sommerfeld results — Orszag's benchmark and spectral convergence, the plane-Poiseuille
spectrum and neutral curve with its critical point, Squire's theorem checked numerically, the Blasius neutral curve in
frequency form (Fig. 11.26 analogue) with Tollmien's profile, the tanh shear layer (Fig. 11.23), Falkner–Skan loops
(Fig. 11.24), plane Couette, Table 11.1 recomputed and the energy budget (11.88) of a TS wave.

Run: ``.venv/Scripts/python.exe scripts/ch11_neutral_curves.py --no-show [--fast]``   Figures -> outputs/ch11/c11_os_poiseuille.png,
c13_neutral_curves.png, c14_energy_budget.png.  Critical points and neutral curves are cached (outputs/ch11/cache).
"""
from __future__ import annotations

import numpy as np

from ch11_common import COLORS, Timer, finish, parse_args, save, setup

from fluidpy import ch11_instability as ch11


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    U, Upp = (lambda y: 1 - y ** 2), (lambda y: -2.0 + 0 * y)
    print("Orszag (1971) benchmark c = 0.23752649 + 0.00373967i (Re = 1e4, k = 1):")
    for N in (40, 60, 80, 100, 120):
        c = ch11.orr_sommerfeld_eigs(1.0, 1e4, U, Upp, N=N, filter=False)[0]
        print(f"  N = {N:3d}: c = {c.real:.10f} {c.imag:+.10f}i")
    sq = ch11.squire_transform(1.0, 0.5, 1e4)
    c3 = ch11.os_3d_eigs(1.0, 0.5, 1e4, U, lambda y: -2 * y, N=60)[0]
    c2 = ch11.orr_sommerfeld_eigs(sq["kbar"], sq["Rebar"], U, Upp, N=100)[0]
    print(f"Squire: 3-D (k, m, Re) = (1, 0.5, 1e4) → c = {c3:.10f}; 2-D (k̄, Re̅) = ({sq['kbar']:.5f}, {sq['Rebar']:.2f}) → "
          f"c = {c2:.10f}; |diff| = {abs(c3 - c2):.1e}")
    with Timer("Poiseuille critical point (cached)"):
        P = ch11.poiseuille_critical()
    print(f"plane Poiseuille: Re_c = {P['Re_c']:.4f}, k_c = {P['k_c']:.5f}, c_r = {P['c_r']:.5f} (Orszag 5772.22, 1.02056, 0.264)")
    with Timer("Poiseuille neutral curve (cached)"):
        pn = ch11.poiseuille_neutral_curve(fast=args.fast)
    j = int(np.argmin(np.abs(pn["Re"] - 1e4)))
    print(f"  at Re = {pn['Re'][j]:.0f}: unstable k from {pn['k_lower'][j]:.4f} to {pn['k_upper'][j]:.4f}")
    print(f"plane Couette: max c_i over k ∈ (0.5, 1, 2), Re ∈ (1e3, 1e4, 1e5) = "
          f"{ch11.couette_max_growth([0.5, 1, 2], [1e3, 1e4, 1e5]):.4f} (< 0: stable)")
    fig, ax = plt.subplots(1, 2, figsize=(14, 4.8))
    sp = ch11.poiseuille_spectrum(1.0, 1e4)
    ax[0].plot(sp.real, sp.imag, ".", color=COLORS["ink"])
    ax[0].plot(sp[0].real, sp[0].imag, "o", color=COLORS["rose"], label=f"TS mode c = {sp[0]:.5f}")
    ax[0].axhline(0, color=COLORS["muted"], lw=1)
    ax[0].set(xlim=(0, 1), ylim=(-1, 0.1), xlabel="c_r", ylabel="c_i", title="Orr–Sommerfeld spectrum, Poiseuille Re = 10⁴, k = 1")
    ax[0].legend(fontsize=8)
    ax[1].semilogx(pn["Re"], pn["k_lower"], color=COLORS["accent"])
    ax[1].semilogx(pn["Re"], pn["k_upper"], color=COLORS["accent"])
    ax[1].plot(P["Re_c"], P["k_c"], "o", color=COLORS["rose"], label=f"Re_c = {P['Re_c']:.1f}")
    ax[1].set(xlabel="Re (half-width, centreline speed)", ylabel="k", title="plane Poiseuille neutral curve (unstable inside)")
    ax[1].legend(fontsize=8)
    save(fig, out, "c11_os_poiseuille")
    with Timer("Blasius critical point (cached)"):
        B = ch11.blasius_critical()
    print(f"Blasius (parallel): Re_δ*,c = {B['Re_c']:.3f}, αδ* = {B['k_c']:.4f}, c_r = {B['c_r']:.4f}, ω = {B['omega_c']:.4f} "
          f"(Thomas 519.2, 0.303, 0.120; Jordinson 520)")
    for N, ym in ((80, 20.0), (100, 25.0)):
        print(f"  convergence (N, y_max) = ({N}, {ym}): Re_c = {ch11.blasius_critical(N=N, y_max=ym)['Re_c']:.3f}")
    with Timer("Blasius neutral curve (cached)"):
        bn = ch11.blasius_neutral_curve(in_frequency=True, fast=args.fast)
    with Timer("tanh neutral curve (cached)"):
        tn = ch11.tanh_shear_layer_neutral_curve(fast=args.fast)
    print("tanh shear layer upper neutral k_u(Re):", ", ".join(f"{R:.3g}: {k:.4f}" for R, k in zip(tn["Re"], tn["k_upper"])))
    with Timer("Bickley critical (cached)"):
        Bj = ch11.bickley_critical()
    print(f"Bickley jet sinuous: Re_c = {Bj['Re_c']:.4f} at k = {Bj['k_c']:.4f}, c_r = {Bj['c_r']:.4f} (Tatsumi & Kakutani ≈ 4.0 at 0.2)")
    with Timer("Falkner–Skan loops (cached)"):
        fav = ch11.falkner_skan_neutral_curve(0.1, fast=args.fast)
        adv = ch11.falkner_skan_neutral_curve(-0.05, fast=args.fast)
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
    ax[0].plot(bn["Re"], bn["F_lower"] * 1e4, color=COLORS["accent"])
    ax[0].plot(bn["Re"], bn["F_upper"] * 1e4, color=COLORS["accent"])
    ax[0].plot(B["Re_c"], B["omega_c"] / B["Re_c"] * 1e4, "o", color=COLORS["rose"], label=f"Re_δ*,c = {B['Re_c']:.1f}")
    ax[0].set(xlabel="Re_δ* = U∞δ*/ν", ylabel="10⁴ ων/U∞²", title="Blasius neutral curve (cf. Fig. 11.26)")
    ax[0].legend(fontsize=8)
    ax[1].plot(tn["Re"], tn["k_upper"], color=COLORS["accent"])
    ax[1].axhline(1.0, color=COLORS["muted"], ls="--", lw=1)
    ax[1].set(xscale="log", xlabel="Re = U₀L/ν", ylabel="kL", title="tanh shear layer: unstable below (cf. Fig. 11.23)")
    for d, col, lab in ((fav, COLORS["teal"], "favourable m = 0.1"), (adv, COLORS["rose"], "adverse m = −0.05")):
        ax[2].semilogx(d["Re"], d["k_lower"], color=col, label=lab)
        ax[2].semilogx(d["Re"], d["k_upper"], color=col)
    ax[2].set(xlabel="Re_δ*", ylabel="kδ*", title="Falkner–Skan neutral curves (cf. Fig. 11.24)")
    ax[2].legend(fontsize=8)
    save(fig, out, "c13_neutral_curves")
    fin = np.isfinite(adv["k_upper"])
    print(f"Falkner–Skan adverse m = −0.05: upper branch k at the largest Re = {adv['k_upper'][fin][-1]:.4f} (flat, inflection); "
          f"favourable m = 0.1: lowest unstable Re in the sample = {fav['Re'][np.isfinite(fav['k_lower'])][0]:.0f}")
    print("Table 11.1 (ours vs published):")
    for row in ch11.table_11_1():
        print(f"  {row['flow']:18s} U = {row['U']:12s} Re_c ours = {row['Re_c_ours']:<10.4g} benchmark {row['benchmark']:<8g} "
              f"({row['source']}); {row['remark']}")
    m = ch11.ts_mode("poiseuille", 1.0, 1e4)
    b = m["budget"]
    print(f"TS mode energy budget (Re = 1e4, k = 1): E = {b['E']:.5f}, dE/dt = {b['dEdt']:.3e}, P = {b['production']:.3e}, "
          f"Λ = {b['dissipation']:.3e}, P/Λ = {b['ratio']:.4f}, residual {b['residual']:.1e}")
    m2 = ch11.ts_mode("poiseuille", 1.0, 5000.0)
    print(f"  Re = 5000, k = 1: c_i = {m2['c'].imag:.5f}, P/Λ = {m2['budget']['ratio']:.4f} (decays)")
    co = ch11.tollmien_coefficients(0.2)
    print(f"Tollmien profile (our η₁ = 0.2): a = {co['a']:.4f}, b = {co['b']:.4f}; jump at η₁: corrected "
          f"{float(ch11.tollmien_profile(0.2 + 1e-12) - ch11.tollmien_profile(0.2 - 1e-12)):.2e}, printed "
          f"{float(ch11.tollmien_profile(0.2 + 1e-12, printed=True) - ch11.tollmien_profile(0.2 - 1e-12, printed=True)):.3f}")
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    y = m["y"]
    ax[0].plot(np.abs(m["u_hat"]), y, color=COLORS["orange"], label="|û|")
    ax[0].plot(np.abs(m["v_hat"]), y, color=COLORS["blue"], label="|v̂|")
    ax[0].set(xlabel="amplitude", ylabel="y", title=f"TS mode c = {m['c']:.5f}")
    ax[0].legend(fontsize=8)
    ax[1].plot(b["production_density"], y, color=COLORS["orange"], label="production −⟨uv⟩U′")
    ax[1].plot(b["dissipation_density"], y, color=COLORS["blue"], label="dissipation density")
    ax[1].set(xlabel="per unit y", ylabel="y", title=f"P/Λ = {b['ratio']:.3f} (> 1: grows)")
    ax[1].legend(fontsize=8)
    eta = np.linspace(0, 1.4, 300)
    from fluidpy.core.boundary_layer import blasius_constants, blasius_profile

    bc = blasius_constants()
    ax[2].plot(blasius_profile(eta * bc["eta99"])[1], eta, color=COLORS["ink"], label="Blasius f′ (δ = δ₉₉)")
    ax[2].plot(ch11.tollmien_profile(eta), eta, color=COLORS["accent"], label="Tollmien, corrected (η₁ = 0.2)")
    ax[2].plot(ch11.tollmien_profile(eta, printed=True), eta, color=COLORS["muted"], ls=":", label="as printed (slip #6)")
    ax[2].set(xlabel="U/U∞", ylabel="y/δ", title="Tollmien's approximate profile")
    ax[2].legend(fontsize=8)
    save(fig, out, "c14_energy_budget")
    if not args.fast:  # --fast never overwrites the published tables
        with Timer("E8 tables → reference/ch11 (cached)"):
            paths = ch11.neutral_curve_tables()
        for k, v in paths.items():
            print(f"  {k:24s} {v}")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
