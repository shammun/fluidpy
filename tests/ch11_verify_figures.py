"""Verification figures for Chapter 11 (our code only) → outputs/ch11/verify/ (git-ignored).

Run: ``.venv/Scripts/python.exe tests/ch11_verify_figures.py``. Writes one PNG per reproduced figure and prints the numbers
quoted in reports/ch11_verification.md. No book figure is copied: every curve is computed by ``fluidpy`` with our own
parameters (figure numbers are the chapter's; the pictures are ours).
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fluidpy import ch11_instability as ch11  # noqa: E402
from fluidpy.core import stability as ST  # noqa: E402

OUT = ROOT / "outputs" / "ch11" / "verify"
OUT.mkdir(parents=True, exist_ok=True)


def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.png", dpi=90)
    plt.close(fig)
    print("wrote", OUT / f"{name}.png")


def fig_11_10():
    t = ch11.benard_neutral_table(write=False)
    fig, ax = plt.subplots(figsize=(7, 5))
    for key, lab in (("rigid", "rigid–rigid (even)"), ("rigid_free", "rigid–free"), ("free", "free–free (11.44)"),
                     ("odd", "rigid–rigid odd (Ex. 11.7)")):
        ax.semilogy(t["K"], t[key], label=lab)
    for bc, md, mk in ((("rigid", "rigid"), "even", "o"), (("rigid", "free"), "any", "s"), (("free", "free"), "even", "^"),
                       (("rigid", "rigid"), "odd", "d")):
        c = ch11.benard_critical(bc=bc, mode=md)
        ax.plot(c["K_c"], c["Ra_c"], mk, color="k")
        print(f"  Bénard {bc} {md}: Ra_c = {c['Ra_c']:.3f} at K_c = {c['K_c']:.4f}")
    ax.set(xlabel="K (wavenumber × d)", ylabel="marginal Ra", ylim=(400, 1e5), title="Fig. 11.10 analogue: neutral curves")
    ax.legend()
    save(fig, "fig_11_10_benard_neutral_curves")


def kh_figure():
    k = np.geomspace(20, 3000, 400)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    for s, lab in ((0.0, "no surface tension"), (0.074, "σ_s = 0.074 N/m")):
        ax[0].loglog(k, np.maximum(ch11.kh_growth_rate(k, 8.0, 0.0, 1.2, 1000.0, surface_tension=s), 1e-3), label=lab)
    ax[0].axvline(ch11.kh_critical_k(8.0, 0.0, 1.2, 1000.0), ls="--", color="grey")
    ax[0].set(xlabel="k [1/m]", ylabel="k c_i [1/s] (floored at 1e-3)", title="KH growth, air over water, ΔU = 8 m/s")
    ax[0].legend()
    ax[1].loglog(k, ch11.kh_stability_boundary(k, 1.2, 1000.0, surface_tension=0.074), label="ΔU_min(k), deep")
    ax[1].loglog(k, ch11.kh_stability_boundary(k, 1.2, 1000.0, surface_tension=0.074, h=0.005), label="h = 5 mm")
    m = ch11.kh_min_shear(1.2, 1000.0)
    ax[1].plot(m["k_star"], m["dU_min"], "o", color="k", label=f"min {m['dU_min']:.2f} m/s at λ = {100 * m['wavelength']:.2f} cm")
    ax[1].set(xlabel="k [1/m]", ylabel="ΔU [m/s]", title="KH stability boundary (Ex. 11.1)")
    ax[1].legend()
    print(f"  KH min shear {m['dU_min']:.3f} m/s at λ = {m['wavelength'] * 100:.3f} cm")
    save(fig, "kh_growth_and_boundary")


def taylor_figure():
    t = ch11.taylor_critical_table()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(t["mu"], t["Ta_c"], "o-", label="narrow gap, exact (Chebyshev)")
    ax.plot(t["mu"], t["approx"], "--", label="(11.54) 1708/(½(1 + μ))")
    ax.axhline(1707.76, color="grey", lw=0.8)
    ax.set(xlabel="μ = Ω₂/Ω₁", ylabel="Ta_c", title="Taylor–Couette onset vs (11.54)")
    ax.legend()
    print(f"  Taylor: μ = 0 Ta_c = {t['Ta_c'][10]:.2f}, (11.54) error {100 * t['rel_error'][10]:.2f} %")
    save(fig, "taylor_critical_vs_11_54")


def tg_map_figure():
    gm = ch11.tg_growth_map()
    fig, ax = plt.subplots(figsize=(7, 5))
    cs = ax.contourf(gm["k"], gm["J"], gm["kci"], levels=20, cmap="magma")
    fig.colorbar(cs, label="k c_i")
    kk = np.linspace(0, 1, 200)
    ax.plot(kk, kk * (1 - kk), "c-", lw=2, label="exact neutral J = k(1 − k)")
    ax.axhline(0.25, color="w", ls="--", label="Ri = 1/4")
    ax.set(xlabel="k", ylabel="J = Ri_min", title="Taylor–Goldstein growth, U = tanh z, N² = J sech² z (cached table)")
    ax.legend(loc="upper right")
    save(fig, "tg_growth_map")


def semicircle_figure():
    fig, ax = plt.subplots(figsize=(7, 4))
    cr, ci = ST.howard_semicircle(-1.0, 1.0)
    ax.plot(cr, ci, "k-", label="Howard semicircle, U ∈ [−1, 1]")
    for J in (0.0, 0.05, 0.1, 0.15, 0.2):
        pr = ch11.richardson_profiles("tanh", J)
        for k in np.linspace(0.1, 0.9, 9):
            c = ST.taylor_goldstein_eigs(k, pr["U"], pr["Upp"], pr["N2"], domain=(-1, 1), N=80, bc="decay",
                                         y_max=ST.decay_box(k), map_scale=ST.decay_map_scale(k))  # loop 2: as tg_growth
            ax.plot(c.real, c.imag, ".", color="C1")
    bk = ch11.parallel_profile("bickley")
    cr2, ci2 = ST.howard_semicircle(0.0, 1.0)
    ax.plot(cr2, ci2, "b--", label="semicircle, U ∈ [0, 1] (Bickley)")
    for k in np.linspace(0.1, 1.9, 10):
        c = ST.rayleigh_eigs(k, bk["U"], bk["Upp"], N=100, bc="decay", y_max=max(40.0, ST.decay_box(k)), map_scale=1.0,
                              unstable_only=True)
        ax.plot(c.real, c.imag, "x", color="C0")
    ax.set_aspect("equal")
    ax.set(xlabel="c_r", ylabel="c_i", title="Fig. 11.20 analogue: computed unstable c inside the semicircles")
    ax.legend(fontsize=8)
    save(fig, "fig_11_20_howard_semicircle")


def os_figure():
    sp_ = ch11.poiseuille_spectrum(1.0, 1e4)
    pn = ch11.poiseuille_neutral_curve()
    P = ch11.poiseuille_critical()
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    ax[0].plot(sp_.real, sp_.imag, ".")
    ax[0].plot(sp_[0].real, sp_[0].imag, "ro", label=f"TS mode {sp_[0]:.6f}")
    ax[0].set(xlim=(0, 1.05), ylim=(-1, 0.05), xlabel="c_r", ylabel="c_i", title="Orr–Sommerfeld spectrum, Re = 10⁴, k = 1")
    ax[0].legend(fontsize=8)
    ax[1].semilogx(pn["Re"], pn["k_lower"], "b-")
    ax[1].semilogx(pn["Re"], pn["k_upper"], "b-")
    ax[1].plot(P["Re_c"], P["k_c"], "ro", label=f"Re_c = {P['Re_c']:.2f}, k_c = {P['k_c']:.5f}")
    ax[1].set(xlabel="Re", ylabel="k", title="plane Poiseuille neutral curve")
    ax[1].legend(fontsize=8)
    save(fig, "os_poiseuille_spectrum_and_neutral_curve")


def tanh_blasius_figure():
    tn = ch11.tanh_shear_layer_neutral_curve()
    bn = ch11.blasius_neutral_curve(in_frequency=True)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    ax[0].semilogx(tn["Re"], tn["k_upper"], "o-")
    ax[0].axhline(1.0, color="grey", ls="--")
    ax[0].set(xlabel="Re = U₀L/ν", ylabel="kL (upper neutral)", ylim=(0, 1.1), title="Fig. 11.23 analogue: tanh shear layer")
    ax[1].semilogx(bn["Re"], bn["F_lower"], "b-")
    ax[1].semilogx(bn["Re"], bn["F_upper"], "b-")
    B = ch11.blasius_critical()
    ax[1].plot(B["Re_c"], B["omega_c"] / B["Re_c"], "ro", label=f"Re_δ*,c = {B['Re_c']:.1f}")
    ax[1].set(xlabel="Re_δ*", ylabel="F = ων/U∞²", title="Fig. 11.26 analogue: Blasius (parallel) neutral curve")
    ax[1].legend()
    save(fig, "fig_11_23_tanh_and_fig_11_26_blasius")


def budget_figure():
    m = ch11.ts_mode()
    b = m["budget"]
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(b["production_density"], m["y"], label="production −⟨uv⟩U′")
    ax.plot(b["dissipation_density"], m["y"], label="dissipation density")
    ax.set(xlabel="density", ylabel="y", title=f"TS mode Re = 10⁴, k = 1: P/Λ = {b['ratio']:.4f}")
    ax.legend()
    save(fig, "fig_11_25_energy_budget_profiles")


def chaos_figure():
    li = ch11.lorenz_integrate(t_end=40.0)
    sep = ch11.lorenz_separation()
    bd = ch11.bifurcation_diagram(np.linspace(2.8, 4.0, 600), n_transient=600, n_keep=80)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
    ax[0].plot(li["X"], li["Z"], lw=0.4)
    for p in ch11.lorenz_fixed_points(28.0)[1:]:
        ax[0].plot(p[0], p[2], "ro")
    ax[0].set(xlabel="X", ylabel="Z", title="Lorenz attractor (Pr = 10, b = 8/3, r = 28)")
    ax[1].semilogy(sep["t"], sep["sep"])
    ax[1].set(xlabel="t", ylabel="|δ|", title=f"separation, slope {sep['slope']:.3f}")
    ax[2].plot(bd["A"], bd["x"], ",k")
    for a in ch11.period_doubling_points(4)["A_n"]:
        ax[2].axvline(a, color="r", lw=0.5)
    ax[2].set(xlabel="A", ylabel="x", title="logistic map bifurcation tree (Fig. 11.31 analogue)")
    save(fig, "fig_11_29_31_chaos")


def profiles_figure():
    fig, axs = plt.subplots(1, 6, figsize=(16, 3.6), sharey=False)
    for ax, d in zip(axs, ch11.fig_11_21_verdicts()):
        pr = ch11.parallel_profile(d["name"], b=math.pi) if d["name"] == "sin" else ch11.parallel_profile(d["name"])
        lo, hi = pr["domain"] if pr["bc"] == "wall" else (0.0, 8.0)
        y = np.linspace(lo, hi, 300)
        ax.plot(pr["U"](y), y)
        for yi in d["y_I"]:
            ax.axhline(yi, color="r", ls=":")
        ax.set_title(f"({d['panel']}) R:{'Y' if d['rayleigh'] else 'N'} F:{'Y' if d['fjortoft'] else 'N'}", fontsize=9)
    save(fig, "fig_11_21_profiles_and_criteria")


def _csv(name):
    lines = [ln for ln in (ROOT / "reference" / "ch11" / name).read_text(encoding="utf-8").splitlines() if ln and ln[0] != "#"]
    arr = np.array([[float(v) for v in ln.split(",")] for ln in lines[1:]])
    return {h: arr[:, j] for j, h in enumerate(lines[0].split(","))}


def bickley_neutral_figure():
    """Post-review loop 2: the regenerated neutral table of the sinuous Bickley jet (wavelength-scaled box)."""
    nb, nl = _csv("os_neutral_bickley.csv"), _csv("os_neutral_bickley_longwave.csv")
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.loglog(nb["Re"], nb["k_upper"], "o-", ms=3, label="upper branch (main band)")
    ax.loglog(nb["Re"], nb["k_lower"], "s-", ms=3, label="lower edge of the main band (gap edge from Re = 19.3)")
    ax.loglog(nl["Re"], nl["k_long_upper"], "^-", ms=4, label="upper edge of the long-wave band")
    nanl = np.isnan(nb["k_lower"])
    ax.loglog(nb["Re"][nanl], np.full(nanl.sum(), 0.02), "rv", ms=4, label="unstable down to k = 0.02 (edge not resolved)")
    cp = ch11.bickley_critical()
    ax.plot(cp["Re_c"], cp["k_c"], "k*", ms=12, label=f"critical point Re = {cp['Re_c']:.3f}, k = {cp['k_c']:.4f}")
    ax.axhline(2.0, color="gray", ls=":", lw=1)
    ax.set_xlabel("Re = U0 L / nu"), ax.set_ylabel("k L"), ax.legend(fontsize=7), ax.set_title("Bickley jet, sinuous: neutral wavenumbers (ours)")
    save(fig, "bickley_neutral_bands_loop2")


def rayleigh_near_neutral_figure():
    """Post-review loop 2: c_i(k) of the Bickley jet up to the neutral point k = 2 — complex path vs the real-axis solver."""
    B = ch11.parallel_profile("bickley")
    ks = np.linspace(1.2, 2.05, 35)
    con, real = [], []
    for k in ks:
        c = ST.rayleigh_eigs_contour(k, B["U"], B["Up"], B["Upp"], N=120, bc="decay", y_max=40.0, map_scale=1.0, parity="even")
        r = ST.rayleigh_eigs(k, B["U"], B["Upp"], N=120, bc="decay", y_max=40.0, map_scale=1.0, unstable_only=True, tol=1e-4, parity="even")
        con.append(c[0].imag if len(c) else 0.0)
        real.append(r[0].imag if len(r) else 0.0)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(ks, con, "-", label="complex path (rayleigh_eigs_contour)")
    ax.plot(ks, real, "x", label="real axis + N-filter (rayleigh_eigs): 0 = nothing returned")
    ax.plot(ks, 0.1109 * (2.0 - ks), ":", color="gray", label="0.111 (2 - k): the linear approach to the exact neutral point")
    ax.set_ylim(-0.005, 0.12), ax.set_xlabel("k L"), ax.set_ylabel("c_i / U0"), ax.legend(fontsize=7)
    ax.set_title("Bickley jet, sinuous mode: growth up to the neutral wavenumber k = 2")
    save(fig, "rayleigh_jet_near_neutral_loop2")
    print("contour c_i at k = 1.8, 1.9:", [round(v, 5) for v, k in zip(con, ks) if abs(k - 1.8) < 0.013 or abs(k - 1.9) < 0.013])


if __name__ == "__main__":
    if "--loop2" in sys.argv:
        bickley_neutral_figure()
        rayleigh_near_neutral_figure()
        raise SystemExit(0)
    bickley_neutral_figure()
    rayleigh_near_neutral_figure()
    fig_11_10()
    kh_figure()
    taylor_figure()
    tg_map_figure()
    semicircle_figure()
    os_figure()
    tanh_blasius_figure()
    budget_figure()
    chaos_figure()
    profiles_figure()
