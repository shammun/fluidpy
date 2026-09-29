"""Verification figures and convergence metrics for Chapter 8 (our code only) → outputs/ch08/verify/ (git-ignored).

Run: ``.venv/Scripts/python.exe tests/ch08_verify_figures.py``. Prints the observed orders quoted in
reports/ch08_verification.md and writes one PNG per reproduced figure. No book figure is copied: every curve is computed
by ``fluidpy`` with our own parameters.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fluidpy import ch08_laminar_flow as ch08  # noqa: E402
from fluidpy.core import diffusion as DIF  # noqa: E402
from fluidpy.core import similarity as SIM  # noqa: E402
from tools.convergence import observed_order, pairwise_orders  # noqa: E402

OUT = ROOT / "outputs" / "ch08" / "verify"
OUT.mkdir(parents=True, exist_ok=True)


def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.png", dpi=110)
    plt.close(fig)
    print("wrote", OUT / f"{name}.png")


def fig_8_4_couette_poiseuille():  # Fig. 8.4: favourable, adverse (backflow), Couette, Poiseuille; |τ| linear
    h, U, mu = 0.01, 0.1, 1e-3
    y = np.linspace(0, h, 201)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    for g, lab in ((-4.0, "dp/dx = −4 Pa/m (favourable)"), (4.0, "dp/dx = +4 Pa/m (adverse)"), (0.0, "dp/dx = 0 (Couette)")):
        ax[0].plot(ch08.channel_flow(y, h, U, g, mu), y * 1e3, label=lab)
    ax[0].plot(ch08.channel_flow(y, h, 0.0, -4.0, mu), y * 1e3, "k--", label="U = 0 (Poiseuille)")
    ax[0].axvline(0, color="0.6", lw=0.8)
    ax[0].set_xlabel("u [m/s]")
    ax[0].set_ylabel("y [mm]")
    ax[0].legend(fontsize=8)
    ax[0].set_title("(8.5): backflow above 2μU/h² = 2 Pa/m")
    ax[1].plot(np.abs(ch08.channel_shear_stress(y, h, 0.0, -4.0, mu)), y * 1e3)
    ax[1].set_xlabel("|τ| [Pa] (Poiseuille)")
    ax[1].set_ylabel("y [mm]")
    save(fig, "fig8_4_couette_poiseuille")


def fig_8_5_8_6_pipe_and_circular_couette():  # Figs. 8.5–8.6: pipe paraboloid + linear τ; circular Couette profiles
    a, mu = 1e-3, 1e-3
    R = np.linspace(0, a, 101)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    ax[0].plot(ch08.pipe_poiseuille(R, a, -1000.0, mu), R * 1e3, label="u_z (8.6)")
    ax2 = ax[0].twiny()
    ax2.plot(ch08.pipe_shear_stress(R, -1000.0), R * 1e3, "r--", label="τ (8.7)")
    ax[0].set_xlabel("u_z [m/s]")
    ax2.set_xlabel("τ [Pa]", color="r")
    ax[0].set_ylabel("R [mm]")
    Rr = np.linspace(0.01, 0.02, 100)
    for O2 in (0.0, 0.5, 1.0, -0.5):
        ax[1].plot(Rr * 100, ch08.circular_couette(Rr, 0.01, 0.02, 1.0, O2) * 1e3, label=f"Ω₂/Ω₁ = {O2}")
    ax[1].set_xlabel("R [cm]")
    ax[1].set_ylabel("u_φ [mm/s]")
    ax[1].legend(fontsize=8)
    ax[1].set_title("(8.10), R₁ = 1 cm, R₂ = 2 cm, Ω₁ = 1 rad/s")
    save(fig, "fig8_5_8_6_pipe_circular_couette")


def fig_8_9_slider_bearing():  # Fig. 8.9 / Example 8.1: exact vs O(α) vs the printed form; load vs taper
    L, h0, U, mu = 0.05, 50e-6, 5.0, 0.05
    x = np.linspace(0, L, 301)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    for al, c in ((0.2, "C0"), (1.0, "C1")):
        ax[0].plot(x * 100, ch08.slider_bearing(x, h0, al, L, U, mu) / 1e5, c, label=f"exact, α = {al}")
        ax[0].plot(x * 100, ch08.slider_bearing(x, h0, al, L, U, mu, model="linear") / 1e5, c + "--", label=f"O(α), α = {al}")
        ax[0].plot(x * 100, ch08.slider_bearing(x, h0, al, L, U, mu, model="book") / 1e5, c + ":", label=f"as printed, α = {al}")
    ax[0].set_xlabel("x [cm]")
    ax[0].set_ylabel("p − p_e [bar]")
    ax[0].legend(fontsize=7)
    als = np.linspace(0.01, 5, 300)
    W = [ch08.slider_bearing_load(1.0, a_, 1.0, 1.0, 1.0) / 6 for a_ in als]
    ax[1].plot(1 + als, W, label="exact W h₀²/(6μUL²)")
    opt = ch08.slider_optimum_taper()
    ax[1].plot(opt["K_opt"], opt["W_star"], "ko", label=f"optimum K = {opt['K_opt']:.4f}, W* = {opt['W_star']:.4f}")
    ax[1].set_xlabel("K = 1 + α (inlet/exit gap)")
    ax[1].legend(fontsize=8)
    save(fig, "fig8_9_slider_bearing")


def fig_8_11_spreading():  # Fig. 8.11 / Examples 8.3, 8.7: thin-film run vs Huppert; front ∝ t^{1/5}
    rho, g, mu = 1000.0, 9.81, 1.0
    N, X = 400, 0.5
    dx = 2 * X / N
    x = -X + dx * (np.arange(N) + 0.5)
    h0 = np.where(np.abs(x) < 0.02, 0.005, 0.0)
    area = float(np.sum(h0) * dx / 2)
    ts = np.logspace(0, 3, 13)
    res = ch08.thin_film_spread(h0, x, ts, rho, g, mu)
    xN = np.array([float(ch08.viscous_current_similarity(0.0, t, area, rho, g, mu, return_front=True)[1]) for t in ts])
    slope = np.polyfit(np.log(ts[-5:]), np.log(res["x_front"][-5:]), 1)[0]
    print(f"thin film: front slope (last 5) = {slope:.4f}; x_front/x_N(1000 s) = {res['x_front'][-1] / xN[-1]:.4f}; "
          f"max volume drift = {np.max(np.abs(res['volume'] / res['volume'][0] - 1)):.2e}")
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    for k in (0, 6, 12):
        ax[0].plot(x * 100, res["h"][k] * 1e3, label=f"numerical t = {ts[k]:.0f} s")
        ax[0].plot(x * 100, ch08.viscous_current_similarity(x, ts[k], area, rho, g, mu) * 1e3, "k:", lw=1)
    ax[0].set_xlabel("x [cm]")
    ax[0].set_ylabel("h [mm]")
    ax[0].legend(fontsize=8)
    ax[0].set_title("dotted: Huppert similarity")
    ax[1].loglog(ts, res["x_front"], "o", label="numerical front")
    ax[1].loglog(ts, xN, "k-", label="η_N(βA³t)^{1/5}")
    ax[1].set_xlabel("t [s]")
    ax[1].set_ylabel("x_N [m]")
    ax[1].legend(fontsize=8)
    save(fig, "fig8_11_spreading")


def fig_8_13_stokes_first():  # Figs. 8.12–8.13: profiles at three times and their collapse; CN convergence
    nu, U = 1e-6, 1.0
    y = np.linspace(0, 0.05, 300)
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.2))
    for t in (10.0, 100.0, 600.0):
        u = ch08.stokes_first_problem(y, t, U, nu)
        ax[0].plot(u, y * 100, label=f"t = {t:.0f} s")
        ax[1].plot(u, ch08.similarity_variable(y, t, nu, half=True), label=f"t = {t:.0f} s")
    ax[0].set_xlabel("u/U")
    ax[0].set_ylabel("y [cm]")
    ax[0].legend(fontsize=8)
    ax[1].set_ylim(0, 3)
    ax[1].set_xlabel("u/U")
    ax[1].set_ylabel("y/(2√(νt))")
    ax[1].set_title("collapse (8.30)")
    t0, t1 = 100.0, 400.0
    Ly = 12 * np.sqrt(nu * t1)
    yg = np.linspace(0, Ly, 801)
    u0 = ch08.stokes_first_problem(yg, t0, U, nu)
    for theta, lab in ((0.5, "Crank–Nicolson"), (1.0, "backward Euler")):
        uref = DIF.crank_nicolson_1d(u0, yg, (t1 - t0) / 5120, 5120, nu, U, 0.0, startup_be=0, t0=t0, theta=theta)
        dts, errs = [], []
        for n in (10, 20, 40, 80):
            errs.append(np.max(np.abs(DIF.crank_nicolson_1d(u0, yg, (t1 - t0) / n, n, nu, U, 0.0, startup_be=0, t0=t0, theta=theta) - uref)))
            dts.append((t1 - t0) / n)
        print(f"{lab} (time): order {observed_order(dts, errs):.3f}, pairwise {[round(p, 3) for p in pairwise_orders(dts, errs)]}")
        ax[2].loglog(dts, errs, "o-", label=f"{lab}: p = {observed_order(dts, errs):.2f}")
    hs, errs = [], []
    for ny in (41, 81, 161, 321):
        yy = np.linspace(0, Ly, ny)
        u = DIF.crank_nicolson_1d(ch08.stokes_first_problem(yy, t0, U, nu), yy, (t1 - t0) / 2000, 2000, nu, U, 0.0, startup_be=0, t0=t0)
        errs.append(np.max(np.abs(u - ch08.stokes_first_problem(yy, t1, U, nu))))
        hs.append(yy[1] - yy[0])
    print(f"CN (space): order {observed_order(hs, errs):.3f}, pairwise {[round(p, 3) for p in pairwise_orders(hs, errs)]}")
    ax[2].set_xlabel("Δt [s]")
    ax[2].set_ylabel("max error")
    ax[2].legend(fontsize=8)
    save(fig, "fig8_13_stokes_first")


def fig_8_14_8_15_vortex_diffusion():  # Figs. 8.14–8.15: thickening vortex sheet; line-vortex decay
    nu = 1e-6
    y = np.linspace(-0.01, 0.01, 400)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    for t in (5.0, 20.0):
        u, w = ch08.vortex_sheet_diffusion(y, t, 1.0, nu)
        ax[0].plot(w, y * 1e3, label=f"ω_z, t = {t:.0f} s")
    ax[0].set_xlabel("ω_z [1/s]")
    ax[0].set_ylabel("y [mm]")
    ax[0].legend(fontsize=8)
    r = np.linspace(0, 0.02, 300)
    for t in (1.0, 50.0, 200.0, 1000.0):
        ax[1].plot(r * 100, ch08.line_vortex_decay(r, t, 0.01, nu), label=f"νt = {nu * t * 1e6:.0f} mm²")
    ax[1].set_xlabel("r [cm]")
    ax[1].set_ylabel("u_θ [m/s]")
    ax[1].legend(fontsize=8)
    save(fig, "fig8_14_8_15_vortex_diffusion")


def fig_8_16_oscillating_plate():  # Fig. 8.16: profiles at ωt = 0, π/2, π, 3π/2 with the envelope
    nu, om = 1e-6, 2 * np.pi
    s = ch08.stokes_layer(nu, om)
    y = np.linspace(0, 8 * s["delta_e"], 300)
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    for ph in (0, np.pi / 2, np.pi, 3 * np.pi / 2):
        ax.plot(ch08.stokes_second_problem(y, ph / om, 1.0, om, nu), y / s["delta_e"], label=f"ωt = {ph:.2f}")
    env = ch08.stokes_layer_envelope(y, 1.0, om, nu)
    ax.plot(env, y / s["delta_e"], "k:", lw=1)
    ax.plot(-env, y / s["delta_e"], "k:", lw=1)
    ax.axhline(s["delta_book"] / s["delta_e"], color="0.5", lw=0.8)
    ax.set_xlabel("u/U")
    ax.set_ylabel("y/δ_e,  δ_e = √(2ν/ω)")
    ax.legend(fontsize=8)
    save(fig, "fig8_16_oscillating_plate")


def fig_8_17_sphere_surface():  # Fig. 8.17: surface stresses and pressure on Stokes' sphere
    th = np.linspace(0, np.pi, 181)
    srr, srt, tx = ch08.stokes_sphere_surface_stresses(th, 1.0, 1.0, 1.0)
    p = ch08.stokes_sphere_pressure(1.0, th, 1.0, 1.0, 1.0)
    fig, ax = plt.subplots(figsize=(6, 4.2))
    ax.plot(np.degrees(th), p, label="(p − p∞) a/(μU)")
    ax.plot(np.degrees(th), srt, label="σ_rθ a/(μU)")
    ax.plot(np.degrees(th), tx, label="t_x a/(μU) (uniform 1.5)")
    ax.set_xlabel("θ from the downstream axis [deg]")
    ax.legend(fontsize=8)
    save(fig, "fig8_17_sphere_surface")


def fig_8_19_8_20_stokes_oseen():  # Figs. 8.19–8.20: fluid-frame streamlines, Stokes symmetric vs Oseen wake
    xg = np.linspace(-6, 6, 400)
    zg = np.linspace(0.001, 5, 220)
    X, Z = np.meshgrid(xg, zg)
    R = np.hypot(X, Z)
    TH = np.arctan2(Z, X)
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
    for a_, (psi, title) in zip(ax, ((ch08.stokes_sphere_streamfunction(R, TH, 1.0, 1.0, "fluid"), "Stokes (fluid frame)"),
                                     (ch08.oseen_streamfunction(R, TH, 1.0, 1.0, 1.0, "fluid"), "Oseen Re = 1 (fluid frame)"))):
        a_.contour(X, Z, psi, levels=np.linspace(-1.5, 0.05, 24), colors="C0", linewidths=0.8)
        a_.add_patch(plt.Circle((0, 0), 1.0, color="0.7"))
        a_.set_aspect("equal")
        a_.set_title(title + " — sphere moves to −x")
    save(fig, "fig8_19_8_20_stokes_oseen")


def fig_drag_curve():  # C_D(Re): Stokes, Oseen, Proudman–Pearson, Morrison correlation
    Re = np.logspace(-2, 1, 200)
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.loglog(Re, ch08.stokes_drag_coefficient(Re), label="Stokes 24/Re (8.52)")
    ax.loglog(Re, ch08.oseen_drag_coefficient(Re), label="Oseen (24/Re)(1 + 3Re/16)")
    ax.loglog(Re, ch08.proudman_pearson_drag_coefficient(Re), label="Proudman–Pearson")
    ax.loglog(Re, SIM.sphere_drag_coefficient(Re, "morrison"), "k--", label="Morrison (2016) correlation")
    ax.set_xlabel("Re = 2aU/ν")
    ax.set_ylabel("C_D")
    ax.legend(fontsize=8)
    save(fig, "fig_drag_curve")


def metrics():
    hs, es = [], []
    for n in (33, 65, 129, 257):
        g = ch08.hele_shaw_streamfunction_grid(n)
        hs.append(g["h"])
        es.append(g["err_far"])
    print(f"Hele-Shaw grid: errors {[f'{e:.3e}' for e in es]}, order {observed_order(hs, es):.3f}")
    L, h0, al, U, mu = 0.05, 50e-6, 0.5, 5.0, 0.05
    hs, es = [], []
    for N in (21, 41, 81, 161, 321):
        x = np.linspace(0, L, N)
        p, _ = ch08.reynolds_pressure_1d(x, h0 * (1 + al * x / L), U_0=-U, mu=mu)
        pe = np.asarray(ch08.slider_bearing(x, h0, al, L, U, mu))
        es.append(np.max(np.abs(p - pe)) / np.max(pe))
        hs.append(x[1] - x[0])
    print(f"Reynolds 1-D (trapezoid): order {observed_order(hs, es):.3f}")
    nu, om = 1e-6, 2 * np.pi
    de = np.sqrt(2 * nu / om)
    Tp = 2 * np.pi / om
    hs, es = [], []
    for ny in (61, 121, 241, 481):
        y = np.linspace(0, 30 * de, ny)
        u = DIF.crank_nicolson_1d(ch08.stokes_second_problem(y, 0.0, 1.0, om, nu), y, Tp / 4000, 4000, nu,
                                  lambda t: np.cos(om * t), 0.0, startup_be=0)
        es.append(np.max(np.abs(u - ch08.stokes_second_problem(y, Tp, 1.0, om, nu))))
        hs.append(y[1] - y[0])
    print(f"Stokes second CN (space): order {observed_order(hs, es):.3f}")
    r = np.array([50.0, 100.0, 200.0, 500.0])
    print(f"inertia/viscous log–log slope: {np.polyfit(np.log(r), np.log(ch08.inertia_viscous_ratio(r, np.pi / 3, 1, 1, 1)), 1)[0]:.4f}")
    w = [np.max(np.hypot(*ch08.oseen_velocity(1.0, np.linspace(0.1, 3.0, 30), 1.0, 1.0, R_))) for R_ in (1e-3, 1e-2, 1e-1)]
    print(f"Oseen wall slip vs Re: order {observed_order([1e-3, 1e-2, 1e-1], w):.3f}")


if __name__ == "__main__":
    fig_8_4_couette_poiseuille()
    fig_8_5_8_6_pipe_and_circular_couette()
    fig_8_9_slider_bearing()
    fig_8_11_spreading()
    fig_8_13_stokes_first()
    fig_8_14_8_15_vortex_diffusion()
    fig_8_16_oscillating_plate()
    fig_8_17_sphere_surface()
    fig_8_19_8_20_stokes_oseen()
    fig_drag_curve()
    metrics()
