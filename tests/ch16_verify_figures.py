"""Reproduce the key plots of chapter 16 with our own code, for the verifier's visual check.

Writes PNGs to ``outputs/ch16/verify/`` (git-ignored).  Every curve is computed by ``fluidpy``; nothing is traced from
the book, and the parameter values are ours.  Run: ``.venv/Scripts/python.exe tests/ch16_verify_figures.py``
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fluidpy import ch16_biofluid_mechanics as ch16  # noqa: E402

OUT = ROOT / "outputs" / "ch16" / "verify"
I = ch16.illustrative_inputs()
PI = math.pi


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.png", dpi=80)
    plt.close(fig)
    print("wrote", OUT / f"{name}.png")


def fig_heart_and_blood():
    """Fig. 16.6 (PV loop), Figs. 16.7–16.8 (Casson), Fig. 16.21 (Casson plug), two-layer model and Pries."""
    fig, ax = plt.subplots(2, 3, figsize=(16, 8.5))
    lp = ch16.pv_loop_schematic(I["EDV"], I["ESV"], p_min=I["p_min"], p_ED=I["p_ED"], p_open=I["p_open"],
                                p_peak=I["p_peak"], p_ES=I["p_ES"])
    w = ch16.pv_loop_work(lp["p"], lp["V"])
    ax[0, 0].plot(lp["V"] * 1e6, ch16.Pa_to_mmHg(lp["p"]), "-")
    for k, (V, p) in lp["corners"].items():
        ax[0, 0].annotate(k, (V * 1e6, ch16.Pa_to_mmHg(p)))
    ax[0, 0].set(xlabel="V [mL]", ylabel="p [mmHg]", title=f"PV loop (ours), W = {w['work']:.3f} J, {w['orientation']}")
    g = np.geomspace(0.05, 2000, 300)
    wc = ch16.whitmore_casson(g, ch16.BLOOD["mu_p"])
    ax[0, 1].plot(np.sqrt(g), np.sqrt(wc["tau"] / ch16.BLOOD["mu_p"]), label="√(τ/μ_p) (16.7)")
    ax[0, 1].set(xlabel="√γ̇ [s^-1/2]", ylabel="√(τ/μ_p) [s^-1/2]", title="Casson: straight line, intercept 2.0")
    ax[0, 1].legend()
    ax[0, 2].loglog(g, wc["mu_app"] * 1e3, label="μ_app (secant)")
    ax[0, 2].loglog(g, ch16.differential_viscosity(g, wc["tau_y"], wc["Kc"]) * 1e3, "--", label="dτ/dγ̇ (tangent)")
    ax[0, 2].axhline(wc["Kc"] ** 2 * 1e3, color="k", lw=0.8, label="K_c² asymptote")
    ax[0, 2].set(xlabel="γ̇ [1/s]", ylabel="viscosity [mPa s]", title="Fig. 16.8 analogue")
    ax[0, 2].legend()
    a, dp, L, ty, Kc = 0.31e-3, 19.0, 0.021, 4.7e-3, 0.057
    r = np.linspace(-a, a, 401)
    ax[1, 0].plot(ch16.casson_pipe_velocity(np.abs(r), a, dp, L, ty, Kc) * 1e3, r * 1e3, label="Casson (16.215)")
    ax[1, 0].plot(ch16.poiseuille_velocity(np.abs(r), a, dp, L, Kc ** 2) * 1e3, r * 1e3, "--", label="Poiseuille, μ = K_c²")
    ac = ch16.casson_plug_radius(dp, L, ty)
    ax[1, 0].axhspan(-ac * 1e3, ac * 1e3, alpha=0.15)
    ax[1, 0].set(xlabel="u [mm/s]", ylabel="r [mm]", title="Fig. 16.21 analogue: plug of radius a_c")
    ax[1, 0].legend()
    rr = np.linspace(0, 20e-6, 300)
    ax[1, 1].plot(ch16.two_layer_velocity(rr, 20e-6, 2.5e-6, 37.0, 1.3e-3, 1.2e-3, 4.0e-3) * 1e3, rr * 1e6)
    ax[1, 1].axhline(17.5, ls=":", color="k")
    ax[1, 1].set(xlabel="u [mm/s]", ylabel="r [µm]", title="two-layer profile (16.27)–(16.28): kink at a − δ")
    D = np.geomspace(20, 500, 60)                         # the closure has no root below about 15 µm at δ = 2.5 µm
    fl = [ch16.fahraeus_lindqvist_model(d * 1e-6 / 2, 2.5e-6, 0.45, 1.2e-3, 310.15)["mu_rel"] for d in D]
    ax[1, 2].semilogx(D, ch16.pries_relative_viscosity(D, 0.45), label="Pries in vitro (Secomb & Pries 2013)")
    ax[1, 2].semilogx(D, fl, "--", label="two-layer + Charm–Kurland, δ = 2.5 µm (ours)")
    ax[1, 2].set(xlabel="D [µm]", ylabel="μ_app/μ_p", title="Fåhræus–Lindqvist effect, H = 0.45")
    ax[1, 2].legend()
    save(fig, "ch16_heart_blood")


def fig_womersley():
    """Fig. 16.14 analogue (profiles over half a cycle), |Q|/|Q0| and phase lag against α, Windkessel beat."""
    fig, ax = plt.subplots(2, 3, figsize=(16, 8.5))
    s = np.linspace(-1, 1, 201)
    for k, al in enumerate((2.5, 5.0, 12.0)):
        axx = ax.flat[k]
        for ph in np.radians(np.arange(0, 181, 30)):
            axx.plot(ch16.womersley_velocity_nd(np.abs(s), ph, al), s, lw=1)
        axx.set(xlabel="u/(Ga²/4μ)", ylabel="r/a", title=f"Womersley profiles, α = {al} (ωt = 0…180°)")
    a = np.geomspace(0.1, 50, 300)
    ax[1, 0].loglog(a, ch16.womersley_flow_ratio(a), label="|Q|/|Q₀|")
    ax[1, 0].loglog(a, 8 / a ** 2, "--", label="8/α²")
    ax2 = ax[1, 0].twinx()
    ax2.semilogx(a, np.degrees(ch16.womersley_phase_lag(a)), "r")
    ax2.set_ylabel("phase lag φ [deg]", color="r")
    ax[1, 0].set(xlabel="α", title="flow amplitude and lag (16.99)–(16.100)")
    ax[1, 0].legend(loc="lower left")
    T = 60 / I["HR"]
    t = np.linspace(0, 3 * T, 1500)
    kw = dict(SV=I["SV"], HR=I["HR"], systole_share=I["systole_share"], R=I["R_periph"], K=I["K_compliance"])
    ax[1, 1].plot(t, ch16.Pa_to_mmHg(ch16.windkessel_periodic_state(t, **kw)), label="p(t) periodic (ours)")
    Qf = ch16.HalfSineEjection(I["SV"], T, I["systole_share"] * T)
    tt = np.linspace(0, 6 * T, 1200)
    ax[1, 1].plot(tt - 3 * T, ch16.Pa_to_mmHg(ch16.windkessel_ode(tt, Qf, I["R_periph"], I["K_compliance"], 0.0)), ":",
                  label="from rest (shifted 3 beats)")
    sm = ch16.windkessel_summary(**kw)
    ax[1, 1].axhline(ch16.Pa_to_mmHg(sm["p_mean"]), ls="--", color="k", lw=0.8, label="true mean = R·CO")
    ax[1, 1].set(xlabel="t [s]", ylabel="p [mmHg]", title="two-element Windkessel (16.49)", xlim=(0, 3 * T))
    ax[1, 1].legend(fontsize=8)
    al = np.geomspace(1, 100, 120)
    cur = ch16.elastic_tube_wave_curve(al)
    ax[1, 2].semilogx(al, cur["c_over_c0"], label="c/c₀ (pressure root)")
    ax[1, 2].semilogx(al, cur["damping_per_wavelength"], label="amplitude kept per wavelength")
    ax[1, 2].set(xlabel="α", title="viscous elastic tube (16.136)", ylim=(0, 1.05))
    ax[1, 2].legend()
    save(fig, "ch16_womersley_windkessel_waves")


def fig_branches_bends_collapse_plants():
    """Fig. 16.15 cost, Figs. 16.17–16.18 Dean, Figs. 16.19–16.20 collapsible tube, §16.4 water potential."""
    fig, ax = plt.subplots(2, 3, figsize=(16, 8.5))
    Q, L, mu, b = 5e-6, 0.05, 3.5e-3, 2000.0
    aa = np.linspace(0.8e-3, 6e-3, 400)
    t = ch16.murray_terms(aa, Q, L, mu, b)
    for k in ("pumping", "upkeep", "total"):
        ax[0, 0].plot(aa * 1e3, t[k], label=k)
    ax[0, 0].axvline(ch16.murray_radius(Q, mu, b) * 1e3, ls="--", color="k")
    ax[0, 0].set(xlabel="a [mm]", ylabel="power [W]", title="Murray cost (16.138), optimum (16.140)", ylim=(0, 3 * t["total"].min()))
    ax[0, 0].legend()
    n = 161
    X, Y = np.meshgrid(np.linspace(-1, 1, n), np.linspace(-1, 1, n))
    R, A = np.hypot(X, Y), np.arctan2(Y, X)
    inside = R <= 1
    psi = np.where(inside, ch16.dean_secondary_streamfunction(np.minimum(R, 1), A), np.nan)
    ax[0, 1].contour(X, Y, psi, 16)
    ax[0, 1].set(aspect="equal", title="Dean ψ₀: two counter-rotating cells (outer bend at x = +1)")
    w = np.where(inside, ch16.dean_axial_velocity(np.minimum(R, 1), A, 576.0), np.nan)
    cs = ax[0, 2].contourf(X, Y, w, 16)
    fig.colorbar(cs, ax=ax[0, 2])
    ax[0, 2].set(aspect="equal", title="axial w₀ + κ²w₁ at K = 576: maximum toward the outer bend")
    al = np.linspace(0.15, 1.6, 300)
    ax[1, 0].plot(al, -ch16.collapsible_tube_law(al), label="shapiro: −P = α⁻ⁿ − 1")
    ax[1, 0].plot(al, -ch16.collapsible_tube_law(al, form="pure"), "--", label="pure: −P = α⁻ⁿ")
    ax[1, 0].set(xlabel="α = A/A₀", ylabel="−P", title="tube law (16.188), n = 3/2", ylim=(-2, 15))
    ax[1, 0].legend()
    prm = dict(L=0.1, A0=PI * 0.004 ** 2, Kp=200.0, n=1.5, rho=1060.0, mu=3.5e-3)
    p2 = np.linspace(-2000, -5, 40)
    fl = ch16.flow_limitation_curve(p2, p1=0.0, pe=0.0, prm=prm)
    ax[1, 1].plot(-p2, fl["Q"] * 1e6, "o-", ms=3)
    ax[1, 1].axvline(-fl["p2_choke"], ls="--", color="k")
    ax[1, 1].set(xlabel="−(p₂ − p_e) [Pa]", ylabel="Q [mL/s]", title="Starling resistor: flow limitation at S = 1")
    z = np.linspace(0, 30, 200)
    col = ch16.xylem_column(z, psi_soil=-0.3e6, u=1e-3, a=20e-6, mu=1e-3)
    ax[1, 2].plot(col["p"] / 1e6, z, label="p")
    ax[1, 2].plot(col["psi"] / 1e6, z, label="ψ")
    ax[1, 2].plot(col["gravity"] / 1e6, z, ":", label="ρgz")
    ax[1, 2].set(xlabel="[MPa]", ylabel="z [m]", title="xylem column (16.221) with Poiseuille drag (ours)")
    ax[1, 2].legend()
    save(fig, "ch16_branches_bends_collapse_plants")


if __name__ == "__main__":
    fig_heart_and_blood()
    fig_womersley()
    fig_branches_bends_collapse_plants()
