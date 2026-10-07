"""Verification figures for Chapter 12 (our code only) → outputs/ch12/verify/ (git-ignored).

Run: ``.venv/Scripts/python.exe tests/ch12_verify_figures.py``. Writes one PNG per reproduced figure and prints the numbers
quoted in reports/ch12_verification.md. No book figure is copied: every curve is computed by ``fluidpy`` with our own
parameters (figure numbers are the chapter's; the pictures are ours). The DNS points are the public Lee & Moser (2015)
subset in reference/ch12 (see SOURCES.md).
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

from fluidpy import ch12_turbulence as ch12  # noqa: E402
from fluidpy.core import turbstats as TS  # noqa: E402
from fluidpy.core import wall_turbulence as WT  # noqa: E402

OUT = ROOT / "outputs" / "ch12" / "verify"
OUT.mkdir(parents=True, exist_ok=True)
REF = ROOT / "reference" / "ch12"
KAP, B = 0.41, 5.0


def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.png", dpi=90)
    plt.close(fig)
    print("wrote", OUT / f"{name}.png")


def dns(case):
    rows = [ln.split(",") for ln in (REF / "lee_moser_2015_channel_mean.csv").read_text(encoding="utf-8").splitlines()
            if ln and ln[0].isdigit()]
    a = np.array([[float(v) for v in r] for r in rows if int(r[0]) == case])
    return a[0, 1], a[:, 3], a[:, 4], a[:, 5]


def fig_12_18_law_of_the_wall():
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
    for case, mk in ((180, "v"), (1000, "s"), (5200, "o")):
        Re, yp, Up, dU = dns(case)
        ax[0].semilogx(yp[1:], Up[1:], mk, ms=3, label=f"DNS Lee & Moser, Re_τ = {Re:.0f}")
        ax[1].semilogx(yp[1:], yp[1:] * dU[1:], mk, ms=3, label=f"DNS, Re_τ = {Re:.0f}")
    y = np.geomspace(0.2, 5200, 300)
    ax[0].semilogx(y[y < 12], y[y < 12], "k--", label="U⁺ = y⁺ (12.82)")
    ax[0].semilogx(y[y > 8], WT.log_law(y[y > 8], kappa=KAP, B=B), "k-", label="log law (12.88), κ = 0.41, B = 5.0")
    ax[0].semilogx(y, WT.spalding_uplus(y, kappa=KAP, B=B), color="tab:red", lw=1, label="Spalding")
    c = ch12.channel_mixing_length(5185.897, KAP, n=800)
    ax[0].semilogx(c["yplus"][1:], c["Uplus"][1:], color="tab:purple", lw=1, label="mixing-length channel (van Driest)")
    ax[0].set(xlabel="y⁺", ylabel="U⁺", title="Fig. 12.18 analogue: law of the wall")
    ax[0].legend(fontsize=7)
    ax[1].axhline(1 / 0.384, color="k", lw=0.8, label="1/κ, κ = 0.384")
    ax[1].axhline(1 / KAP, color="grey", lw=0.8, ls=":", label="1/κ, κ = 0.41")
    ax[1].set(xlabel="y⁺", ylabel="y⁺ dU⁺/dy⁺", ylim=(0, 6), title="log-law indicator (12.87)")
    ax[1].legend(fontsize=7)
    save(fig, "fig_12_18_law_of_the_wall")
    Re, yp, Up, dU = dns(5200)
    for lo in (30.0, 100.0, 3 * math.sqrt(Re), 350.0):
        f = WT.fit_log_law(yp, Up, window=(lo, 0.15), Re_tau=Re)
        print(f"  DNS 5200 fit from y+ = {lo:6.1f}: kappa = {f['kappa']:.4f}, B = {f['B']:.3f}, n = {f['n_points']}, rms = {f['rms_residual']:.4f}")
    ypl = np.geomspace(1.0, 5000.0, 400)
    comp = WT.composite_profile_plus(ypl, 5000.0, kappa=KAP, B=B, Pi=0.0)
    for lo in (30.0, 100.0, 200.0, 300.0):
        f = WT.fit_log_law(ypl, comp, window=(lo, 0.15), Re_tau=5000.0)
        print(f"  Spalding δ+ = 5000 fit from y+ = {lo:5.0f}: kappa = {f['kappa']:.4f}, B = {f['B']:.3f}")
    ypL = np.geomspace(1.0, 1e5, 2000)
    compL = WT.composite_profile_plus(ypL, 1e5, kappa=KAP, B=B, Pi=0.0)
    for lo in (30.0, 100.0, 300.0, 1000.0):
        f = WT.fit_log_law(ypL, compL, window=(lo, 0.15), Re_tau=1e5)
        print(f"  Spalding δ+ = 1e5 fit from y+ = {lo:5.0f}: kappa = {f['kappa']:.4f}, B = {f['B']:.3f}")
    for case in (180, 550, 1000, 2000, 5200):
        Re, yp, Up, dU = dns(case)
        c = ch12.channel_mixing_length(Re, KAP, n=800)
        dev = np.max(np.abs(np.interp(yp, c["yplus"], c["Uplus"]) - Up))
        yd = yp / Re
        Ub = np.trapezoid(Up, yd) / yd[-1]
        inner = (yp > 0) & (yd < 0.15)
        sdev = np.max(np.abs(WT.spalding_uplus(yp[inner], kappa=KAP, B=B) - Up[inner]))
        print(f"  Re_tau {Re:8.2f}: mixing-length max|dU+| = {dev:.3f}, Cf/Cf_DNS - 1 = {c['Cf'] / (2 / Ub ** 2) - 1:+.4f}, "
              f"Spalding max|dU+| (y/δ < 0.15) = {sdev:.3f}")


def fig_12_12_spectrum():
    fig, ax = plt.subplots(figsize=(7, 5))
    eps = 1.0
    for nu, lab in ((1e-3, "L/η = 180"), (1e-5, "L/η = 5.6e3"), (1e-7, "L/η = 1.8e5")):
        eta = (nu ** 3 / eps) ** 0.25
        k1 = np.geomspace(0.3, 3 / eta, 60)
        S11 = 0.5 * ch12.one_dimensional_from_3d(k1, lambda K: ch12.model_spectrum(K, eps, nu, L=1.0))
        x, Phi = ch12.kolmogorov_normalize_spectrum(k1, S11, nu, eps)
        ax.loglog(x, 2 * Phi, label=lab)
    xx = np.geomspace(1e-5, 1, 50)
    ax.loglog(xx, 18 * 1.5 / 55 * xx ** (-5 / 3), "k--", label="C₁ (k₁η)^(−5/3), C₁ = 18C/55")
    ax.set(xlabel="k₁η", ylabel="2S₁₁/(u_K²η)", ylim=(1e-4, 1e9), title="Fig. 12.12 analogue: model spectra in Kolmogorov scaling")
    ax.legend(fontsize=8)
    save(fig, "fig_12_12_spectrum_collapse")


def fig_12_13_plane_jet():
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    kw = dict(C5="from_invariant", xi_half=0.1)
    for x in (1.0, 2.0, 4.0):
        y = np.linspace(-0.4 * x, 0.4 * x, 801)
        U = ch12.plane_jet_mean_velocity(x, y, 1.0, 1.0, **kw)
        ax[0].plot(y, U, label=f"x = {x:g} m")
        ax[1].plot(y / x, U / U.max(), label=f"x = {x:g} m")
        ax[2].plot(y / x, ch12.plane_jet_reynolds_stress(x, y, 1.0, 1.0, **kw) / U.max() ** 2, label=f"x = {x:g} m")
        print(f"  jet x = {x}: J = {ch12.jet_momentum_flux_per_span(np.linspace(-2 * x, 2 * x, 8001), ch12.plane_jet_mean_velocity(x, np.linspace(-2 * x, 2 * x, 8001), 1.0, 1.0, **kw), 1.0):.10f}"
              f", Vdot = {ch12.plane_jet_volume_flux(x, 1.0, 1.0, **kw):.5f}")
    ax[0].set(xlabel="y [m]", ylabel="U [m/s]", title="plane jet (12.66): wider and slower")
    ax[1].set(xlabel="ξ = y/x", ylabel="U/U_CL", title="the same profiles collapse")
    ax[2].set(xlabel="ξ = y/x", ylabel="−mean(uv)/U_CL²", title="stress profile C₃G = −½F∫F (12.63)")
    for a in ax:
        a.legend(fontsize=8)
    save(fig, "fig_12_13_plane_jet_similarity")


def fig_12_10_energy_budget():
    b = ch12.channel_energy_budget(1000.0, KAP, 26.0, n=800)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    y = b["yplus"]
    for key, lab in (("pressure_work", "pressure work"), ("transport", "transport"), ("viscous_dissipation", "viscous dissipation"),
                     ("loss_to_turbulence", "loss to turbulence")):
        ax[0].semilogx(y[1:], b["mean"][key][1:], label=lab)
    ax[0].set(xlabel="y⁺", ylabel="terms of (12.46) [u_*⁴/ν]", title="mean-flow energy budget, model channel Re_τ = 1000")
    ax[1].semilogx(y[1:], b["turbulence"]["production"][1:], label="production")
    ax[1].semilogx(y[1:], b["turbulence"]["dissipation_plus_transport"][1:], label="dissipation + transport (residual)")
    ax[1].axhline(0.25, color="grey", lw=0.6, ls=":")
    ax[1].set(xlabel="y⁺", ylabel="terms of (12.47)", title="turbulent energy budget (production ≤ ¼)")
    for a in ax:
        a.legend(fontsize=8)
    save(fig, "fig_12_10_channel_energy_budgets")
    I = b["integrals"]
    print(f"  channel 1000: work = {I['work']:.4f}, dissipation = {I['dissipation']:.4f}, production = {I['production']:.4f}, "
          f"identity residual = {b['identity_residual']:.2e}, peak production {np.max(b['production']):.4f} at y+ = {b['yplus_peak_production']:.2f}")


def fig_12_25_dispersion():
    t = np.concatenate([[0.0], np.geomspace(0.02, 300.0, 60)])
    X, u = ch12.langevin_particles(10000, t, 1.0, 10.0, seed=2)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    ax[0].loglog(t[1:], np.mean(X ** 2, axis=0)[1:], "o", ms=3, label="10⁴ Langevin particles")
    ax[0].loglog(t[1:], ch12.taylor_dispersion_exponential(t[1:], 1.0, 10.0), "k-", label="Taylor (12.119), r = e^(−τ/Λ)")
    ax[0].loglog(t[1:], t[1:] ** 2, "--", color="grey", label="u² t² (12.120)")
    ax[0].loglog(t[1:], 20 * t[1:], ":", color="grey", label="2u²Λ t (12.122)")
    ax[0].set(xlabel="t [s]", ylabel="mean(X²) [m²]", ylim=(1e-4, 1e4), title="Fig. 12.25/12.26 analogue")
    ax[0].legend(fontsize=8)
    x = np.linspace(0.0, 2000.0, 400)
    ax[1].plot(x, ch12.smoke_plume_width(x, 5.0, 0.5, 20.0), color="tab:blue")
    ax[1].plot(x, -ch12.smoke_plume_width(x, 5.0, 0.5, 20.0), color="tab:blue")
    ax[1].axvline(5.0 * 20.0, color="grey", ls=":", label="x = UΛ_t")
    ax[1].set(xlabel="x [m]", ylabel="±Z_rms [m]", title="Fig. 12.27 analogue: linear near the source, √x far away")
    ax[1].legend(fontsize=8)
    save(fig, "fig_12_25_taylor_dispersion")


def fig_12_21_surface_layer():
    z = np.geomspace(0.03, 100.0, 200)
    fig, ax = plt.subplots(figsize=(7, 5))
    us, z0, kap = 0.3, 0.03, 0.4
    ax.semilogy(WT.surface_layer_wind(z, us, z0, kappa=kap), z, "k-", label="neutral (12.93)")
    ax.semilogy(WT.surface_layer_wind(z, us, z0, 80.0, kappa=kap), z, color="tab:blue", label="stable, L_M = +80 m (log-linear)")
    ax.semilogy(WT.surface_layer_wind(z, us, z0, -20.0, kappa=kap, unstable="businger_dyer"), z, color="tab:red",
                label="unstable, L_M = −20 m (Businger–Dyer)")
    ax.semilogy(WT.surface_layer_wind(z, us, z0, -20.0, kappa=kap), z, color="tab:red", ls="--", label="unstable, log-linear (outside its range)")
    ax.set(xlabel="U [m/s]", ylabel="z [m]", xlim=(-1, 12), title="Fig. 12.21 analogue: wind in the surface layer, u_* = 0.3 m/s")
    ax.legend(fontsize=8)
    save(fig, "fig_12_21_surface_layer_wind")


def fig_12_5_correlation_spectrum():
    dt, n = 0.05, 2 ** 17
    u = TS.make_ensemble(1, np.arange(n) * dt, 0.0, 1.0, 0.5, seed=3)[0]
    lag, R = TS.autocorrelation(u, dt, max_lag=100)
    om, S = TS.periodogram(u, dt, segments=128)
    p = TS.correlation_spectrum_pair("exponential", 1.0, 0.5)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))
    ax[0].plot(lag, R / R[0], "o", ms=3, label="record, 2¹⁷ samples")
    ax[0].plot(lag, p["r"](lag), "k-", label="e^(−τ/τ_c)")
    ax[0].set(xlabel="τ [s]", ylabel="r₁₁(τ)", title=f"Fig. 12.5 analogue: Λ_t = {TS.integral_scale(lag, R):.3f} s (exact 0.5)")
    ax[1].loglog(om[1:], S[1:], ".", ms=3, label="periodogram, 128 segments")
    ax[1].loglog(om[1:], p["S"](om[1:]), "k-", label="σ²τ_c/[π(1 + ω²τ_c²)]")
    ax[1].set(xlabel="ω [rad/s]", ylabel="S_e(ω)", title="spectrum of the same record (12.20)")
    for a in ax:
        a.legend(fontsize=8)
    save(fig, "fig_12_5_correlation_and_spectrum")


if __name__ == "__main__":
    fig_12_18_law_of_the_wall()
    fig_12_12_spectrum()
    fig_12_13_plane_jet()
    fig_12_10_energy_budget()
    fig_12_25_dispersion()
    fig_12_21_surface_layer()
    fig_12_5_correlation_spectrum()
