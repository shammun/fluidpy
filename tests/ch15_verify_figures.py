"""Reproduce the key plots of chapter 15 with our own code, for the verifier's visual check.

Writes PNGs to ``outputs/ch15/verify/`` (git-ignored).  Every curve is computed by ``fluidpy``; nothing is traced from
the book.  Run: ``.venv/Scripts/python.exe tests/ch15_verify_figures.py``
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

from fluidpy import ch15_compressible_flow as ch15  # noqa: E402
from fluidpy.core import acoustics as AC  # noqa: E402
from fluidpy.core import gasdynamics as GD  # noqa: E402

OUT = ROOT / "outputs" / "ch15" / "verify"
I = ch15.illustrative_inputs()
PI = math.pi


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.png", dpi=80)
    plt.close(fig)
    print("wrote", OUT / f"{name}.png")


def fig_sound_and_nozzles():
    """Fig. 15.1 (pulse), steepening, Table 15.1 as curves, Figs. 15.5–15.6, 15.12, 15.13."""
    fig, ax = plt.subplots(2, 3, figsize=(16, 8.5))
    c = 340.0
    x = np.linspace(-3, 3, 601)
    f = AC.pulse("gaussian", 0.4, 20.0)
    for t, col in ((0.0, "C0"), (2e-3, "C1"), (5e-3, "C2")):
        ax[0, 0].plot(x, AC.dalembert(x, t, c, f=lambda s: 0.5 * f(s), g=lambda s: 0.5 * f(s)), col, label=f"p′, t = {1e3 * t:.0f} ms")
        ax[0, 0].plot(x, 1.2 * c * AC.plane_wave_velocity(x, t, c, 1.2, f=lambda s: 0.5 * f(s), g=lambda s: 0.5 * f(s)), col + "--", lw=1)
    ax[0, 0].set(title="d'Alembert (15.15): a pulse released from rest splits in two; dashed: ρ₀c u′ (15.16)", xlabel="x [m]", ylabel="p′ [Pa]")
    ax[0, 0].legend(fontsize=7)
    xs = np.linspace(0, 1, 801)
    u0 = 30.0 * np.sin(2 * PI * xs)
    tb = AC.breaking_time(None, None, 1.4, max_compression_rate=30.0 * 2 * PI)
    for fr, col in ((0.0, "C0"), (0.5, "C1"), (1.0, "C3"), (1.5, "0.6")):
        s = AC.simple_wave_evolve(u0, xs, fr * tb, c, 1.4)
        ax[0, 1].plot(s["x_points"] - c * fr * tb, s["u"], col, label=f"t = {fr:.1f} t_b")
    ax[0, 1].set(title="steepening of a finite wave (ours): vertical tangent at t_b, multivalued after", xlabel="x − c₀t [m]", ylabel="u [m/s]")
    ax[0, 1].legend(fontsize=7)
    M = np.linspace(0.01, 4, 400)
    s = GD.isentropic_state(M)
    for key, lab in (("p_over_p0", "p/p₀ (15.29)"), ("rho_over_rho0", "ρ/ρ₀ (15.30)"), ("T_over_T0", "T/T₀ (15.28)")):
        ax[0, 2].plot(M, s[key], label=lab)
    ax[0, 2].plot(M, np.minimum(s["A_over_Astar"], 6) / 6, "k--", label="A/A* ÷ 6 (15.31)")
    ax[0, 2].axvline(1.0, color="0.7", lw=0.8)
    ax[0, 2].set(title="isentropic ratios, γ = 1.4; A/A* has its minimum 1 at M = 1", xlabel="M", ylim=(0, 1.05))
    ax[0, 2].legend(fontsize=7)
    xn = np.linspace(0, 1, 501)
    for br, As, col in (("subsonic", 0.8, "C0"), ("subsonic", None, "C2"), ("transonic", None, "C3"), ("supersonic", None, "C1")):
        n = GD.nozzle_isentropic(xn, GD.laval_area, branch=br, A_star=As)
        ax[1, 0].plot(xn, n["M"], col, label=f"{br}" + (" (venturi, A* = 0.8 A_t)" if As else ""))
    ax[1, 0].plot(xn, GD.laval_area(xn), "k:", label="A(x)/A_t")
    ax[1, 0].axhline(1.0, color="0.7", lw=0.8)
    ax[1, 0].set(title="M(x) along one nozzle shape: four isentropic solutions (Figs. 15.5–15.6)", xlabel="x/L", ylabel="M")
    ax[1, 0].legend(fontsize=7)
    cr = GD.laval_critical_pressures(2.0)
    for r in (0.98, cr["p_b"], 0.85, 0.75, 0.6, cr["p_d"], 0.3, cr["p_f"]):
        n = GD.laval_nozzle(np.linspace(0, 1, 801), GD.laval_area, r)
        ax[1, 1].plot(n["x"], n["p_over_p0"], label=f"p_B/p₀ = {r:.3f}: {n['regime']}")
        ax[1, 1].plot([1.0, 1.06], [n["p_over_p0"][-1], r], "k:", lw=0.8)
    ax[1, 1].axhline(cr["p_star"], color="0.7", lw=0.8)
    ax[1, 1].set(title="Laval nozzle, A_e/A_t = 2: p(x)/p₀ for a family of back pressures (Fig. 15.13)", xlabel="x/L", ylabel="p/p₀")
    ax[1, 1].legend(fontsize=6)
    r = np.linspace(0, 1, 401)
    for g in (1.25, 1.4, 5 / 3):
        m = GD.mass_flow_vs_back_pressure(r, g)
        ax[1, 2].plot(r, m["mdot"], label=f"γ = {g:.3g}; p*/p₀ = {GD.critical_pressure_ratio(g):.3f}")
    ax[1, 2].set(title="convergent nozzle (Fig. 15.12): ṁ√(RT₀)/(p₀A_e) against p_B/p₀ — a plateau below p*", xlabel="p_B/p₀")
    ax[1, 2].legend(fontsize=7)
    save(fig, "A_sound_and_nozzles")


def fig_shocks_and_ducts():
    """Figs. 15.8, 15.9, 15.11, 15.15, 15.16 and the friction/heating lines in the T–s plane."""
    fig, ax = plt.subplots(2, 3, figsize=(16, 8.5))
    M1 = np.linspace(0.3780, 4, 600)
    ax[0, 0].plot(M1, GD.normal_shock_mach(M1), "C0", label="(15.40)")
    ax[0, 0].plot([0, 4], [0, 4], "k:", lw=0.8, label="no change")
    ax[0, 0].axhline(GD.strong_shock_limits()["M2"], color="C3", ls="--", lw=0.8, label="√((γ−1)/2γ) = 0.378")
    ax[0, 0].axvspan(0, 1, color="0.92")
    ax[0, 0].set(title="M₂ against M₁ (Fig. 15.8); M₁ < 1 (shaded) is forbidden", xlabel="M₁", ylabel="M₂", xlim=(0, 4), ylim=(0, 4))
    ax[0, 0].legend(fontsize=7)
    M1 = np.linspace(0.5, 5, 500)
    ax[0, 1].plot(M1, GD.normal_shock_entropy(M1), "C0", label="(15.43)")
    Mw = np.linspace(0.7, 1.6, 100)
    ax[0, 1].plot(Mw, GD.weak_shock_entropy(Mw), "C3--", label="(15.44a): cubic in M₁² − 1")
    ax[0, 1].axhline(0, color="0.6", lw=0.8)
    ax[0, 1].axvline(1, color="0.6", lw=0.8)
    ax[0, 1].set(title="entropy jump (S₂ − S₁)/C_v (Fig. 15.9): negative below M₁ = 1, flat contact at 1", xlabel="M₁", ylim=(-0.5, 1.5))
    ax[0, 1].legend(fontsize=7)
    eta = np.linspace(-10, 6, 400)
    for Mv in (1.3, 2.0, 3.5):
        U1, U2 = ch15.becker_shock_roots(Mv)
        ax[0, 2].plot(eta, ch15.becker_shock_profile(eta, Mv), label=f"M₁ = {Mv}: U₁ = {U1:.3f}, U₂ = {U2:.3f}")
    ax[0, 2].set(title="Becker's shock structure, Pr″ = 1 (Fig. 15.11): U = u/V between two plateaus", xlabel="η", ylabel="U", ylim=(0, 0.95))
    ax[0, 2].legend(fontsize=7)
    Mg = np.linspace(0.05, 4, 800)
    for k, (key, vals, a_) in enumerate((("f", (0.0, 0.1, 0.2, 0.3, 0.4), ax[1, 0]), ("q", (-0.4, -0.2, 0.0, 0.1, 0.2, 0.4), ax[1, 1]))):
        for v in vals:
            cur = GD.duct_curve(Mg, **{key: v})
            line, = a_.plot(cur["M1"], cur["M2_continuous"], label=f"{key} = {v}")
            a_.plot(cur["M1"], cur["M2_shocked"], color=line.get_color(), ls=":", lw=0.9)
        a_.plot([0, 4], [1, 1], color="0.7", lw=0.8)
        a_.set(title=f"duct with {'friction' if key == 'f' else 'heat'} (Fig. 15.{15 + k}): M₂ against M₁; dotted = the other root; gaps = choked",
               xlabel="M₁", ylabel="M₂", xlim=(0, 4), ylim=(0, 4))
        a_.legend(fontsize=7)
    Mt = np.linspace(0.25, 3.5, 400)
    fa, ra = GD.fanno_line_Ts(Mt), GD.rayleigh_line_Ts(Mt)
    ax[1, 2].plot(fa["ds_over_cp"], fa["T_over_Tstar"], "C0", label="friction (Fanno) line")
    ax[1, 2].plot(ra["ds_over_cp"], ra["T_over_Tstar"], "C3", label="heating (Rayleigh) line")
    ax[1, 2].plot([0], [1], "ko", ms=4)
    ax[1, 2].set(title="T–s lines (ours): entropy is largest at M = 1; heating line peaks in T at M = 1/√γ", xlabel="(s − s*)/C_p", ylabel="T/T*", xlim=(-1.2, 0.05))
    ax[1, 2].legend(fontsize=7)
    save(fig, "B_shocks_and_ducts")


def fig_planar_waves_and_plate():
    """Figs. 15.17, 15.19, 15.23, 15.24 and the linear-against-exact comparisons."""
    fig, ax = plt.subplots(2, 3, figsize=(16, 8.5))
    for Mv in (1.2, 1.4, 1.6, 2.0, 3.0, 5.0, 1e6):
        c = GD.oblique_shock_curve(Mv, n=600)
        ax[0, 0].plot(np.degrees(c["sigma"]), np.degrees(c["delta"]), label="M₁ → ∞" if Mv > 100 else f"M₁ = {Mv}")
    Ml = np.concatenate([np.linspace(1.001, 6, 200), [10, 30, 1e6]])
    md, sd = GD.max_deflection(Ml), GD.sonic_deflection(Ml)
    ax[0, 0].plot(np.degrees(md["sigma_max"]), np.degrees(md["delta_max"]), "k--", lw=1, label="δ_max locus")
    ax[0, 0].plot(np.degrees(sd["sigma_sonic"]), np.degrees(sd["delta_sonic"]), "k:", lw=1, label="M₂ = 1 locus")
    ax[0, 0].set(title="deflection δ against shock angle σ (15.53), γ = 1.4 (Fig. 15.19)", xlabel="σ [deg]", ylabel="δ [deg]", xlim=(0, 90), ylim=(0, 47))
    ax[0, 0].legend(fontsize=6, ncol=2)
    M = np.linspace(1, 12, 400)
    for g in (1.25, 1.4, 5 / 3):
        line, = ax[0, 1].plot(M, np.degrees(GD.prandtl_meyer(M, g)), label=f"γ = {g:.3g}")
        ax[0, 1].axhline(np.degrees(GD.prandtl_meyer_max(g)), color=line.get_color(), ls=":", lw=0.8)
    ax[0, 1].plot(M, np.degrees(GD.mach_angle(M)), "k--", lw=1, label="Mach angle μ (15.49)")
    ax[0, 1].set(title="Prandtl–Meyer function ν(M) (15.56) and its limits ν_max; μ = sin⁻¹(1/M)", xlabel="M", ylabel="deg")
    ax[0, 1].legend(fontsize=7)
    d = np.radians(np.linspace(0.2, 22, 80))
    for Mv in (1.5, 2.0, 3.0):
        dm = GD.max_deflection(Mv)["delta_max"]
        dd = d[d < 0.98 * dm]
        ax[0, 2].plot(np.degrees(dd), [GD.weak_wave_error(Mv, float(v))["ratio"] for v in dd], label=f"shock, M = {Mv}")
        ax[0, 2].plot(np.degrees(d), [GD.weak_wave_error(Mv, -float(v))["ratio"] for v in d], ls="--", label=f"expansion, M = {Mv}")
    ax[0, 2].axhline(1, color="0.6", lw=0.8)
    ax[0, 2].set(title="exact pressure change ÷ weak-wave rule (15.55): → 1 linearly from both sides", xlabel="|δ| [deg]")
    ax[0, 2].legend(fontsize=6)
    al = np.radians(np.linspace(0.2, 20, 60))
    for Mv, col in ((1.6, "C0"), (2.0, "C1"), (3.0, "C2")):
        ok = al < 0.98 * GD.max_deflection(Mv)["delta_max"]
        ex = [ch15.flat_plate_supersonic(float(a), Mv, method="shock_expansion") for a in al[ok]]
        ax[1, 0].plot(np.degrees(al[ok]), [e["CL"] for e in ex], col, label=f"shock–expansion, M∞ = {Mv}")
        ax[1, 0].plot(np.degrees(al), ch15.ackeret_flat_plate(al, Mv)["CL"], col + "--", lw=1)
        ax[1, 1].plot(np.degrees(al[ok]), [e["CD"] for e in ex], col, label=f"shock–expansion, M∞ = {Mv}")
        ax[1, 1].plot(np.degrees(al), ch15.ackeret_flat_plate(al, Mv)["CD"], col + "--", lw=1)
    ax[1, 0].set(title="flat plate: C_L against α; dashed: 4α/√(M∞² − 1) (15.59)", xlabel="α [deg]", ylabel="C_L")
    ax[1, 1].set(title="flat plate: wave drag C_D; dashed: 4α²/√(M∞² − 1) (15.59)", xlabel="α [deg]", ylabel="C_D")
    ax[1, 0].legend(fontsize=7)
    ax[1, 1].legend(fontsize=7)
    c0 = 340.0
    for k, (Mv, col) in enumerate(((0.5, "C0"), (2.0, "C3"))):
        fr = AC.source_wavefronts(Mv * c0, c0, 6.0, 5, 1.0)
        th = np.linspace(0, 2 * PI, 200)
        off = k * 4200.0
        for xc, r in zip(fr["x_centres"], fr["radii"]):
            ax[1, 2].plot(xc - fr["x_source"] + r * np.cos(th), off + r * np.sin(th), col, lw=0.9)
        ax[1, 2].plot([0], [off], "ko", ms=4)
        if Mv > 1:
            mu = fr["mu"]
            for sgn in (1, -1):
                ax[1, 2].plot([0, -4200], [off, off + sgn * 4200 * math.tan(mu)], "k--", lw=0.9)
    ax[1, 2].set(title="wave fronts of a moving source (Fig. 15.17): M = 0.5 (below), M = 2 with its Mach cone (above)", xlabel="x − x_source [m]", aspect="equal",
                 xlim=(-4500, 2200), ylim=(-2200, 6800))
    save(fig, "C_planar_waves_and_plate")


if __name__ == "__main__":
    fig_sound_and_nozzles()
    fig_shocks_and_ducts()
    fig_planar_waves_and_plate()
