"""Reproduce the key plots of chapter 13 with our own code, for the verifier's visual check.

Writes PNGs to ``outputs/ch13/verify/`` (git-ignored).  Every curve is computed by ``fluidpy``; nothing is traced from
the book.  Run: ``.venv/Scripts/python.exe tests/ch13_verify_figures.py``
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

from fluidpy import ch13_geophysical_fluid_dynamics as ch13  # noqa: E402

OUT = ROOT / "outputs" / "ch13" / "verify"
I = ch13.illustrative_inputs()
G0 = ch13.G0
PI = math.pi


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.png", dpi=110)
    plt.close(fig)
    print("wrote", OUT / f"{name}.png")


def fig_ekman():
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.4))
    for lat, col in ((I["lat_ekman"], "C0"), (-I["lat_ekman"], "C3")):
        f = float(ch13.coriolis_parameter(lat))
        d = float(ch13.ekman_depth(I["nu_v_ocean"], f))
        z = np.linspace(-6 * d, 0, 600)
        u, v = ch13.ekman_surface(z, I["tau"], 0.0, I["rho_ocean"], I["nu_v_ocean"], f)
        ax[0].plot(u, v, col, label=f"f {'>' if f > 0 else '<'} 0")
        for s in (0, PI / 4, PI / 2, PI):
            us, vs = ch13.ekman_surface(-s * d, I["tau"], 0.0, I["rho_ocean"], I["nu_v_ocean"], f)
            ax[0].plot([us], [vs], "o", color=col, ms=4)
    ax[0].axhline(0, color="0.7", lw=0.6)
    ax[0].axvline(0, color="0.7", lw=0.6)
    ax[0].set(xlabel="u [m/s]", ylabel="v [m/s]", title="surface spiral: stress along +x; dots at -z/δ = 0, π/4, π/2, π", aspect="equal")
    ax[0].legend()
    f = float(ch13.coriolis_parameter(I["lat_ekman"]))
    d = float(ch13.ekman_depth(I["nu_v_atm"], f))
    z = np.linspace(0, 7 * d, 700)
    u, v = ch13.ekman_bottom(z, I["U_g"], 0.0, I["nu_v_atm"], f)
    ax[1].plot(u / I["U_g"], v / I["U_g"], "C0")
    ax[1].plot([0, 0.35], [0, 0.35], "k:", lw=0.8, label="45°")
    ax[1].set(xlabel="u/U", ylabel="v/U", title="bottom layer hodograph (f > 0)", aspect="equal")
    ax[1].legend()
    ax[2].plot(u / I["U_g"], z / d, "C0", label="u/U")
    ax[2].plot(v / I["U_g"], z / d, "C1", label="v/U")
    ax[2].axvline(1.0, color="0.6", lw=0.6)
    ax[2].axhline(3 * PI / 4, color="C2", ls="--", lw=0.8, label="3π/4: largest u = 1.067 U")
    ax[2].axhline(PI, color="C3", ls=":", lw=0.8, label="π: v = 0, u = 1.043 U")
    ax[2].set(xlabel="velocity / U", ylabel="z/δ", title="bottom layer profiles")
    ax[2].legend(fontsize=8)
    save(fig, "fig_ekman_layers")


def fig_modes():
    N, H = I["ocean_N"], I["ocean_H"]
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.6))
    X = np.linspace(1e-3, 3.3 * PI, 4000)
    t = np.tan(X)
    t[np.abs(t) > 6] = np.nan
    ax[0].plot(X, t, "C0", label="tan X")
    eps = N * N * H / G0
    ax[0].plot(X, 200 * eps / X, "C1", label="(N²H/g)/X, × 200")
    roots = N * H / ch13.modes_uniform_N(N, H, n_modes=4).c
    ax[0].plot(roots, np.tan(roots), "ko", ms=4, label="roots (true height)")
    ax[0].set(ylim=(-1, 3), xlabel="X = N H / c", title=f"(13.69): X₀ = {roots[0]:.4f}, then just above nπ")
    ax[0].legend(fontsize=8)
    z = np.linspace(-H, 0, 201)
    fd = ch13.vertical_modes(z, N ** 2, n_modes=4)
    ex = ch13.modes_uniform_N(N, H, n_modes=4)
    for n in range(4):
        ax[1].plot(ex.psi[n], ex.z, f"C{n}", label=f"n = {n}, c = {ex.c[n]:.3g} m/s")
        ax[1].plot(fd.psi[n][::10], z[::10], "o", color=f"C{n}", ms=3)
    ax[1].set(xlabel="ψ_n", ylabel="z [m]", title="uniform N, free surface: exact (lines), finite volume (dots)")
    ax[1].legend(fontsize=8)
    zt = np.linspace(-H, 0, 801)
    mt = ch13.vertical_modes(zt, ch13.thermocline_N2(zt), n_modes=4)
    for n in range(4):
        ax[2].plot(mt.psi[n], zt, f"C{n}", label=f"n = {n}, c = {mt.c[n]:.3g} m/s")
    ax[2].set(xlabel="ψ_n", ylabel="z [m]", ylim=(-1500, 0), title="thermocline N(z): modes trapped near the surface")
    ax[2].legend(fontsize=8)
    save(fig, "fig_vertical_modes")


def fig_dispersion():
    f0, beta = float(ch13.coriolis_parameter(I["lat"])), float(ch13.beta_parameter(I["lat"]))
    c = float(ch13.baroclinic_mode_speed(I["ocean_N"], I["ocean_H"], 1))
    Lam = c / f0
    K = np.geomspace(0.02, 200, 400) / Lam
    wm, wr, wp = ch13.shallow_water_omega(-K, 0.0, c, f0, beta)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    ax[0].loglog(K * Lam, np.array(wp) / f0, "C0", label="fast root of the cubic")
    ax[0].loglog(K * Lam, np.array(wr) / f0, "C3", label="slow root (k < 0)")
    ax[0].loglog(K * Lam, ch13.poincare_omega(K, f0, c) / f0, "k--", lw=0.8, label="(13.82)")
    ax[0].loglog(K * Lam, ch13.kelvin_omega(K, c) / f0, "C2:", label="Kelvin ω = cK")
    ax[0].loglog(K * Lam, ch13.rossby_omega(-K, 0.0, beta, f0, c) / f0, "k:", lw=0.8, label="(13.118)")
    ax[0].axhline(1.0, color="0.6", lw=0.6)
    ax[0].set(xlabel="K Λ", ylabel="ω / f", title="first baroclinic mode, 35° N: a gap of four decades", ylim=(1e-5, 300))
    ax[0].legend(fontsize=8)
    k = -np.linspace(0.01, 5, 300) / Lam
    wmax = ch13.rossby_max_frequency(beta, f0, c)["omega_max"]
    ax[1].plot(k * Lam, ch13.rossby_omega(k, 0.0, beta, f0, c) / wmax, "C3")
    ax[1].plot([-1], [1], "ko")
    ax[1].set(xlabel="k Λ", ylabel="ω / ω_max", title="Rossby waves, l = 0: maximum at kΛ = −1")
    save(fig, "fig_dispersion")


def fig_adjustment():
    H, eta0 = 1.3, 0.05
    f = float(ch13.coriolis_parameter(I["lat"]))
    c = math.sqrt(G0 * H)
    Lam = c / f
    n = 600
    dx = 120 * Lam / n
    x = (np.arange(n) + 0.5) * dx - 60 * Lam
    dt = 0.5 * dx / c
    nper = int(round(2 * PI / f / dt))
    r = ch13.linear_1d_run(eta0 * np.sign(x), dx, dt, 6 * nper, H, f)
    ee, ve = ch13.geostrophic_adjustment_1d(x, eta0, H, f)
    em, vm = r["eta"][5 * nper + 1:].mean(0), r["v"][5 * nper + 1:].mean(0)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    ax[0].plot(x / Lam, eta0 * np.sign(x), "0.6", label="start")
    ax[0].plot(x / Lam, ee, "C0", label="end state (D16)")
    ax[0].plot(x[::4] / Lam, em[::4], "k.", ms=3, label="model, mean over the 6th inertial period")
    ax[0].plot(x / Lam, r["eta"][nper // 2], "C1", lw=0.7, label="model at half an inertial period")
    ax[0].set(xlim=(-8, 8), xlabel="x / Λ", ylabel="η [m]", title="geostrophic adjustment of a step")
    ax[0].legend(fontsize=8)
    ax[1].plot(x / Lam, ve, "C0", label="jet (D16)")
    ax[1].plot(x[::4] / Lam, vm[::4], "k.", ms=3, label="model mean")
    _, vs = ch13.geostrophic_adjustment_1d(x, eta0, H, -f)
    ax[1].plot(x / Lam, vs, "C3--", label="f < 0")
    ax[1].set(xlim=(-8, 8), xlabel="x / Λ", ylabel="v [m/s]", title=f"jet along the step, peak gη₀/c = {G0 * eta0 / c:.3f} m/s")
    ax[1].legend(fontsize=8)
    save(fig, "fig_adjustment")


def fig_eady():
    N, H, U0 = I["atm_N"], I["atm_H"], I["atm_U0"]
    f = float(ch13.coriolis_parameter(I["lat"]))
    aH = np.linspace(0.02, 3.6, 500)
    c = ch13.eady_phase_speed(aH, 1.0)
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    ax[0].plot(aH, aH * c.imag, "C0", label="α H c_i / U₀ from (13.141)")
    an = np.linspace(0.2, 3.4, 17)
    num = [ch13.eady_numeric_eigs(a * abs(f) / (N * H), 0.0, N, f, H, U0, n=32)["growth_rate"] * N * H / (abs(f) * U0) for a in an]
    ax[0].plot(an, num, "k.", label="Chebyshev eigenvalues")
    fast, ac = ch13.eady_fastest(), ch13.eady_critical()
    ax[0].axvline(ac, color="C3", ls="--", lw=0.8, label=f"cut-off {ac:.4f}")
    ax[0].plot([fast["alphaH"]], [fast["sigma_nd"]], "C1o", label=f"fastest: {fast['alphaH']:.4f}, {fast['sigma_nd']:.4f}")
    ax[0].set(xlabel="α H", ylabel="σ N H / (f U₀)", title="Eady growth rate")
    ax[0].legend(fontsize=8)
    ax[1].plot(aH, c.real, "C0", label="c_r / U₀ (upper root)")
    ax[1].plot(aH, ch13.eady_phase_speed(aH, 1.0, -1).real, "C0--", label="lower root")
    ax[1].plot(aH, c.imag, "C3", label="c_i / U₀")
    ax[1].axhline(1 / (2 * math.sqrt(3)), color="0.6", lw=0.6)
    ax[1].set(xlabel="α H", title="phase speed: c_r = U₀/2 while growing, two edge waves beyond")
    ax[1].legend(fontsize=8)
    z = np.linspace(0, H, 101)
    md = ch13.eady_mode(z, fast["alphaH"], U0, H)
    ax[2].plot(md["amplitude"], z / H, "C0", label="|p̂|")
    ax[2].plot(md["phase"] / PI, z / H, "C1", label="phase / π")
    ax[2].set(xlabel="", ylabel="z / H", title="fastest mode: phase rises by π/2 (westward tilt)")
    ax[2].legend(fontsize=8)
    save(fig, "fig_eady")


def fig_step_and_kuo():
    beta, f0 = float(ch13.beta_parameter(I["lat"])), float(ch13.coriolis_parameter(I["lat"]))
    U, h0, h1 = I["U_mean"], 4000.0, 3800.0
    Yp = f0 * (h1 - h0) / (beta * h0)
    x = np.linspace(-9e6, 9e6, 2001)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    ax[0].plot(x / 1e6, ch13.flow_over_step(x, U, beta, f0, h0, h1)["Y"] / 1e3, "C0", label="eastward (flows to +x)")
    ax[0].plot(x / 1e6, ch13.flow_over_step(x, -U, beta, f0, h0, h1)["Y"] / 1e3, "C3", label="westward (flows to −x)")
    ax[0].axhline(Yp / 1e3, color="0.5", ls="--", lw=0.8, label="Y_p = f₀(h₁ − h₀)/(β h₀)")
    ax[0].axvline(0, color="k", lw=0.6)
    ax[0].set(xlabel="x [1000 km]", ylabel="northward displacement Y [km]", title="flow over a 5 % step (shallower for x·sign(U) > 0)")
    ax[0].legend(fontsize=8)
    Uf = lambda y: 1 / np.cosh(y) ** 2  # noqa: E731
    Up = lambda y: -2 * np.tanh(y) / np.cosh(y) ** 2  # noqa: E731
    Upp = lambda y: (4 * np.tanh(y) ** 2 - 2 / np.cosh(y) ** 2) / np.cosh(y) ** 2  # noqa: E731
    betas = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.64, 0.7])
    gr = [ch13.rayleigh_kuo_eigs(1.4, Uf, Up, Upp, b, bc="decay", y_max=16.0, parity="even")["growth_rate"] for b in betas]
    ax[1].plot(betas, gr, "o-", color="C0")
    ax[1].axvline(2 / 3, color="C3", ls="--", lw=0.8, label="β = max U″ = 2/3")
    ax[1].set(xlabel="β L²/U₀", ylabel="k c_i", title="sech² jet, kL = 1.4: β stabilises (Rayleigh–Kuo)")
    ax[1].legend(fontsize=8)
    save(fig, "fig_step_and_rayleigh_kuo")


def fig_models():
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    kb = ch13.load_reference_run("kelvin_basin")
    Lam = float(kb["Lambda"])
    for i, col in ((0, "C0"), (2, "C1"), (5, "C2"), (9, "C3")):
        ax[0].contour(kb["x"] / Lam, kb["y"] / Lam, kb["eta"][i], levels=[0.1, 0.25], colors=col)
    ax[0].set(xlabel="x / Λ", ylabel="y / Λ", aspect="equal", title="Kelvin pulse, f > 0: frames 0, 2, 5, 9 (anticlockwise)")
    tf, tb = ch13.load_reference_run("turbulence_f"), ch13.load_reference_run("turbulence_beta")
    ax[1].semilogy(tf["t"], tf["energy"] / tf["energy"][0], "C0", label="energy")
    ax[1].semilogy(tf["t"], tf["enstrophy"] / tf["enstrophy"][0], "C3", label="enstrophy")
    ax[1].semilogy(tf["t"], tf["K_E"] / tf["K_E"][0], "C0--", label="K_E")
    ax[1].semilogy(tb["t"], tb["K_E"] / tb["K_E"][0], "C2--", label="K_E with β")
    ax[1].set(xlabel="t", title="decaying 2-D turbulence: enstrophy falls, energy goes to large scales")
    ax[1].legend(fontsize=8)
    pv = ch13.load_reference_run("pv_particles")
    ax[2].plot(pv["t"] / 86400, (pv["q"] / pv["q"][0] - 1) * 100, lw=0.7)
    ax[2].set(xlabel="t [days]", ylabel="change of q on a particle [%]", title="PV on 24 marked particles (spread of q: 44 %)")
    save(fig, "fig_models")


if __name__ == "__main__":
    for fn in (fig_ekman, fig_modes, fig_dispersion, fig_adjustment, fig_eady, fig_step_and_kuo, fig_models):
        fn()
