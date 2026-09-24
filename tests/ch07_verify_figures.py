"""Verification figures and convergence metrics for Chapter 7 (our code only) → outputs/ch07/verify/ (git-ignored).

Run: ``.venv/Scripts/python.exe tests/ch07_verify_figures.py``. Prints the observed orders quoted in
reports/ch07_verification.md and writes one PNG per reproduced figure. No book figure is copied: every curve is computed
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

from fluidpy import ch07_gravity_waves as ch07  # noqa: E402
from fluidpy.core import waves as W  # noqa: E402
from tools.convergence import observed_order, pairwise_orders  # noqa: E402

OUT = ROOT / "outputs" / "ch07" / "verify"
OUT.mkdir(parents=True, exist_ok=True)
TWO_PI = 2 * np.pi
G = 9.81


def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.png", dpi=110)
    plt.close(fig)
    print("wrote", OUT / f"{name}.png")


def fig_orbits():  # Figs. 7.3–7.4: deep circles, intermediate and shallow ellipses, clockwise
    fig, axs = plt.subplots(1, 3, figsize=(12, 4.2))
    for ax, (kH, title) in zip(axs, ((3.0, "deep, kH = 3"), (1.0, "intermediate, kH = 1"), (0.25, "shallow, kH = 0.25"))):
        k, a = 1.0, 0.12
        H = kH / k
        om = float(W.omega_gravity(k, H, G))
        t = np.linspace(0, TWO_PI / om, 200)
        for z0 in np.linspace(-0.08 * H, -0.92 * H, 5) if kH < 3 else np.linspace(-0.1, -2.4, 5):
            xi, ze = ch07.orbit_linear(0.0, z0, t, a, k, H, G)
            ax.plot(xi, z0 + ze, lw=1.2)
            ax.annotate("", xy=(xi[8], z0 + ze[8]), xytext=(xi[0], z0 + ze[0]),
                        arrowprops=dict(arrowstyle="->", lw=0.8))
        ax.axhline(0, color="k", lw=0.6)
        if np.isfinite(H) and kH < 3:
            ax.axhline(-H, color="saddlebrown", lw=2)
        ax.set_aspect("equal")
        ax.set_title(title)
        ax.set_xlabel("x − x₀ [m]")
        ax.set_ylabel("z [m]")
    save(fig, "fig7_4_orbits")


def fig_dispersion():  # Fig. 7.10: c(λ) with capillary and gravity branches, deep and shallow limits
    lam = np.logspace(-3, 3, 600)
    k = TWO_PI / lam
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for H in (0.1, 1.0, 10.0, np.inf):
        ax.loglog(lam, W.phase_speed(k, H, G, 0.0727, 1000.0), label=f"H = {H} m (σ = 0.0727 N/m)")
    ax.loglog(lam, np.sqrt(G * lam / TWO_PI), "k--", lw=0.8, label="deep √(gλ/2π) (7.45)")
    ax.loglog(lam, np.sqrt(TWO_PI * 0.0727 / (1000 * lam)), "k:", lw=0.8, label="capillary √(2πσ/ρλ) (7.60)")
    m = W.capillary_minimum(0.0727, 1000.0, G)
    ax.plot(m["lam_m"], m["c_min"], "ro", label=f"c_min = {m['c_min']:.4f} m/s at λ_m = {100 * m['lam_m']:.3f} cm")
    ax.set_ylim(0.1, 40)
    ax.set_xlabel("λ [m]")
    ax.set_ylabel("c [m/s]")
    ax.legend(fontsize=7)
    save(fig, "fig7_10_phase_speed")


def fig_residual_scan():  # C02: size of the dropped terms vs ka (abs slope 2, rel slope 1)
    ka = np.logspace(-2.5, -1, 8)
    fig, ax = plt.subplots(figsize=(6, 4.2))
    for kH in (1.0, np.inf):
        s = ch07.free_surface_residual_scan(ka, kH=kH, g=G)
        ax.loglog(ka, s["kinematic_abs"], "o-", label=f"|kinematic (7.16)| kH={kH}: slope {s['slope_abs']:.3f}")
        ax.loglog(ka, s["kinematic_rel"], "s--", label=f"÷ aω: slope {s['slope_rel']:.3f}")
        print("C02 residual slopes kH", kH, s["slope_abs"], s["slope_abs_dynamic"], s["slope_rel"], s["slope_rel_dynamic"])
    ax.set_xlabel("ka")
    ax.set_ylabel("residual")
    ax.legend(fontsize=7)
    save(fig, "c02_residual_scan")


def fig_packet():  # Figs. 7.13–7.15: packet envelope rides at c_g, crests at c
    x = np.linspace(-300, 900, 2 ** 13, endpoint=False)
    k0 = 1.0
    eta0 = np.exp(-x ** 2 / (2 * 20.0 ** 2)) * np.cos(k0 * x)
    ts = np.array([0.0, 100.0, 200.0])
    ev = W.linear_evolve(eta0, x, ts, direction=+1, g=G)
    cg = float(W.group_velocity(k0, np.inf, G))
    c = float(W.phase_speed(k0, np.inf, G))
    fig, ax = plt.subplots(figsize=(9, 4))
    for i, t in enumerate(ts):
        ax.plot(x, ev[i] + 2.5 * i, lw=0.6)
        ax.plot(x, W.envelope(ev[i]) + 2.5 * i, "k", lw=0.8)
        ax.axvline(cg * t, color="purple", ls="--", lw=0.8)
        ax.axvline(c * t, color="orange", ls=":", lw=0.8)
    ax.set_xlim(-100, 500)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("η (offset per time)")
    ax.set_title("packet at t = 0, 100, 200 s: envelope at c_g (purple), c (orange)")
    save(fig, "fig7_15_packet_cg")


def fig_jump():  # Fig. 7.20 + (7.81): conjugate depth ratio and head loss vs Fr₁
    Fr = np.linspace(0.3, 8, 200)
    r = [ch07.hydraulic_jump(1.0, Fr1=f, g=G)["ratio"] for f in Fr]
    dE = [ch07.hydraulic_jump(1.0, Fr1=f, g=G)["dE"] for f in Fr]
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].plot(Fr, r)
    ax[0].axvline(1, color="k", lw=0.6)
    ax[0].set_xlabel("Fr₁")
    ax[0].set_ylabel("H₂/H₁ (7.81)")
    ax[1].plot(Fr, dE)
    ax[1].axhline(0, color="k", lw=0.6)
    ax[1].axvline(1, color="k", lw=0.6)
    ax[1].set_xlabel("Fr₁")
    ax[1].set_ylabel("E₂ − E₁ [J/kg] (H₁ = 1 m)")
    save(fig, "fig7_20_jump")


def fig_stokes():  # Figs. 7.21–7.22: Stokes profile; dyed line leaning forward; drift convergence
    x = np.linspace(0, 2 * TWO_PI, 800)
    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    for order in (1, 2, 3):
        ax[0].plot(x, ch07.stokes_wave_profile(x, 0.0, 0.25, 1.0, G, order=order), label=f"order {order}")
    ax[0].axhline(0, color="k", lw=0.5)
    ax[0].legend()
    ax[0].set_title("Stokes profile ka = 0.25 (7.82)")
    a, k, H = 0.05, 1.0, 3.0
    T = TWO_PI / float(W.omega_gravity(k, H, G))
    z0s = np.linspace(-0.05, -2.5, 25)
    for n in (0, 5, 10, 20):
        d = ch07.dyed_line(z0s, n * T, a, k, H, G)
        ax[1].plot(d["x"], d["z"], "o-", ms=2, label=f"t = {n}T")
    ax[1].set_title("dyed line (exact path lines)")
    ax[1].set_xlabel("x [m]")
    ax[1].legend(fontsize=7)
    ka = np.array([0.01, 0.02, 0.04, 0.08])
    rel = [abs(ch07.stokes_drift_numeric(-0.5, e, 1.0, 3.0, G) / float(ch07.stokes_drift(-0.5, e, 1.0, 3.0, G)) - 1)
           for e in ka]
    ax[2].loglog(ka, rel, "o-")
    ax[2].set_title(f"numeric vs (7.86): order {observed_order(ka, rel):.3f}")
    ax[2].set_xlabel("ka")
    print("Stokes drift order", observed_order(ka, rel), pairwise_orders(ka, rel))
    save(fig, "fig7_22_stokes")


def fig_kdv():  # Fig. 7.23: a hump splitting into solitons; invariants
    H = 1.0
    x = np.linspace(-100, 100, 512, endpoint=False)
    eta0 = 0.3 * np.exp(-(x / 4.0) ** 2)
    ts = np.linspace(0, 20, 5)
    r = ch07.kdv_solve(eta0, x, ts, H, G)
    fig, ax = plt.subplots(figsize=(9, 4))
    for i, t in enumerate(ts):
        ax.plot(x, r["eta"][i] + 0.35 * i, lw=0.9)
    ax.set_xlabel("x [m]")
    ax.set_title("KdV (7.87): hump → solitons (offset per time)")
    print("KdV invariant drift", np.ptp(r["mass"]) / abs(r["mass"][0]), np.ptp(r["momentum"]) / r["momentum"][0],
          np.ptp(r["energy"]) / abs(r["energy"][0]))
    save(fig, "fig7_23_kdv")


def fig_two_layer():  # Figs. 7.27–7.28: barotropic and baroclinic dispersion vs kH; long-wave limit
    k = np.logspace(-4, 0, 300)
    H, r1, r2 = 40.0, 1000.0, 1004.0
    wbt, wbc = W.two_layer_free_surface_omega(k, H, r1, r2, G)
    fig, ax = plt.subplots(figsize=(7, 4.3))
    ax.loglog(k * H, wbt / k, label="barotropic ω/k (7.111)")
    ax.loglog(k * H, wbc / k, label="baroclinic ω/k (7.113)")
    ax.axhline(float(W.two_layer_long_wave_speed(H, r1, r2, G)), color="k", ls="--", lw=0.8, label="√(g′H) (7.116)")
    ax.loglog(k * H, W.interface_omega(k, r1, r2, G) / k, ":", label="deep interface (7.95)")
    ax.set_xlabel("kH")
    ax.set_ylabel("c [m/s]")
    ax.legend(fontsize=8)
    save(fig, "fig7_27_two_layer")


def fig_internal():  # Figs. 7.29, 7.33: St Andrew's cross; c ⟂ c_g
    x = z = np.linspace(-12, 12, 301)
    fig, ax = plt.subplots(1, 2, figsize=(11, 5))
    d = ch07.st_andrews_cross(x, z, 0.0, 0.6, 1.0, 1.0, detail=True)
    ax[0].pcolormesh(x, z, d["field"], cmap="RdBu_r", shading="auto")
    for b in d["beams"]:
        ax[0].plot([0, 11 * b["e_beam"][0]], [0, 11 * b["e_beam"][1]], "k--", lw=0.6)
    ax[0].set_aspect("equal")
    ax[0].set_title(f"beams at {np.degrees(d['theta']):.2f}° from the vertical (ω/N = 0.6)")
    s = ch07.internal_wave_state(0.6, 1.0, 1.0, k_sign=-1.0, m_sign=1.0)
    ax[1].quiver(0, 0, s["k"], s["m"], angles="xy", scale_units="xy", scale=1, color="orange", label="K")
    ax[1].quiver(0, 0, s["cx"], s["cz"], angles="xy", scale_units="xy", scale=1, color="red", label="c")
    ax[1].quiver(0, 0, s["cgx"], s["cgz"], angles="xy", scale_units="xy", scale=1, color="purple", label="c_g")
    ax[1].set_xlim(-1, 1)
    ax[1].set_ylim(-1, 1)
    ax[1].set_aspect("equal")
    ax[1].legend()
    ax[1].set_title("k < 0 (Fig. 7.29 geometry): c up-left, c_g down-left")
    save(fig, "fig7_29_7_33_internal")


def fig_fenton():  # V5: explicit approximations vs our inverse dispersion
    xnd = np.logspace(-2, 2, 2000)
    om = xnd * np.sqrt(G)
    kd = W.wavenumber_from_omega(om, 1.0, G)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.semilogx(xnd, 100 * (kd / W.fenton_mckee_kh(om, 1.0, G) - 1), label="Fenton & McKee (1990), λ error")
    ax.semilogx(xnd, 100 * (kd / W.guo_kh(om, 1.0, G) - 1), label="Guo (β = 5/2), λ error")
    ax.axhline(1.7, color="k", ls=":", lw=0.8)
    ax.axhline(-1.7, color="k", ls=":", lw=0.8)
    ax.set_xlabel("ω√(d/g)")
    ax.set_ylabel("error [%]")
    ax.legend(fontsize=8)
    save(fig, "v5_fenton_guo_errors")


def fig_seiche():  # C08: our Steklov FD vs (7.65); convergence
    sys.path.insert(0, str(ROOT / "tests"))
    import test_ch07 as T

    L, H = 10.0, 2.0
    exact = np.array([float(ch07.seiche_modes(L, H, n, G)["omega"]) for n in range(4)])
    nxs = [20, 40, 80, 160]
    errs = [np.max(np.abs(T._seiche_steklov_fd(L, H, nx) / exact - 1)) for nx in nxs]
    print("seiche FD order", observed_order([L / n for n in nxs[1:]], errs[1:]), pairwise_orders([L / n for n in nxs], errs))
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.loglog([L / n for n in nxs], errs, "o-")
    ax.set_xlabel("Δx [m]")
    ax.set_ylabel("max rel. error of ω₀…ω₃")
    ax.set_title("sloshing eigenproblem (FD, ours) → (7.65)")
    save(fig, "c08_seiche_fd_convergence")


if __name__ == "__main__":
    for f in (fig_orbits, fig_dispersion, fig_residual_scan, fig_packet, fig_jump, fig_stokes, fig_kdv, fig_two_layer,
              fig_internal, fig_fenton, fig_seiche):
        f()
