"""Verification figures and convergence metrics for Chapter 4 (our code only) → outputs/ch04/verify/ (git-ignored).

Run: ``.venv/Scripts/python.exe tests/ch04_verify_figures.py``. Prints the observed orders quoted in
reports/ch04_verification.md and writes one PNG per reproduced figure.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import sympy as sp  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fluidpy import ch04_conservation_laws as ch04  # noqa: E402
from tools.convergence import observed_order, pairwise_orders  # noqa: E402

OUT = ROOT / "outputs" / "ch04" / "verify"
OUT.mkdir(parents=True, exist_ok=True)
G = 9.81


def fig_projectile():
    W, u0 = 0.3, 2.0
    t = np.linspace(0, 6, 61)
    r = ch04.coriolis_projectile(u0, W, t)
    inert, rot = ch04.projectile_paths(u0, W, np.linspace(0, 6, 400))
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    ax[0].plot(inert[0], inert[1], "k-", label="inertial frame: straight line")
    ax[0].set(title="Fig. 4.8 analogue — seen from outer space", xlabel="x [m]", ylabel="y [m]")
    ax[1].plot(rot[0], rot[1], "C1-", label="closed form (u₀t cos Ωt, −u₀t sin Ωt)")
    ax[1].plot(r["rotating"][0], r["rotating"][1], "C0o", ms=3, label="ODE a′ = −2Ω×u′ − Ω×(Ω×x′)")
    ax[1].set(title="rotating frame (Ω = 0.3 rad/s, NH): deflected right", xlabel="x′ [m]", ylabel="y′ [m]")
    ax[0].set_ylim(-1, 1)
    ax[1].set_aspect("equal")
    for a in ax:
        a.legend(fontsize=8)
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "fig4_8_projectile_verify.png", dpi=110)
    plt.close(fig)
    return float(np.max(np.abs(r["rotating"] - ch04.projectile_paths(u0, W, t)[1])))


def fig_effective_gravity():
    lat = np.linspace(0, 90, 361)
    mag, dev = ch04.effective_gravity(np.deg2rad(lat))
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(lat, np.rad2deg(dev) * 60, "C3")
    ax.set(xlabel="latitude φ [°]", ylabel="angle between g_e and the radial [arcmin]",
           title="Fig. 4.9 analogue — deflection of effective gravity")
    ax2 = ax.twinx()
    ax2.plot(lat, mag, "C0--")
    ax2.set_ylabel("|g_e| [m/s²]", color="C0")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "fig4_9_effective_gravity_verify.png", dpi=110)
    plt.close(fig)
    return float(lat[np.argmax(dev)]), float(np.rad2deg(dev.max()) * 60)


def fig_sphere_drag():
    Re = np.logspace(-1, 6, 400)
    d = ch04.synthetic_sphere_drag_data(n=80, seed=0, noise=0.03)
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.loglog(Re, ch04.sphere_drag_coefficient(Re), "k-", label="Morrison (2016) correlation")
    ax.loglog(Re[Re < 50], 24 / Re[Re < 50], "C1--", label="Stokes 24/Re")
    for f, m in (("air", "o"), ("water", "s"), ("glycerine", "^")):
        k = d["fluid"] == f
        ax.loglog(d["Re"][k], d["CD"][k], m, ms=4, label=f"synthetic data: {f}")
    ax.set(xlabel="Re = ρUd/μ", ylabel="C_D = F/(½ρU²πd²/4)", title="Fig. 4.21 analogue — sphere drag collapse")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(OUT / "fig4_21_sphere_drag_verify.png", dpi=110)
    plt.close(fig)
    return float(ch04.sphere_drag_coefficient(1e4)), float(ch04.sphere_drag_coefficient(1e6))


def fig_meniscus():
    sig, rho = 0.0728, 998.0
    dl = np.sqrt(sig / (rho * G))
    fig, ax = plt.subplots(figsize=(6.5, 4))
    worst = 0.0
    for th in (10, 30, 60, 85):
        x, z = ch04.meniscus_profile_ode(np.deg2rad(th), x_max=4 * dl, sigma=sig, rho=rho)
        m = z > 1e-3 * z[0]
        xc = ch04.meniscus_profile_x(z[m], np.deg2rad(th), sig, rho)
        worst = max(worst, float(np.max(np.abs(xc - x[m])) / dl))
        ax.plot(x / dl, z / dl, "-", label=f"θ = {th}° (ODE)")
        ax.plot(xc[::12] / dl, z[m][::12] / dl, "k.", ms=4)
    ax.set(xlabel="x/δ", ylabel="ζ/δ", title="Example 4.7 analogue — meniscus (dots: closed form x(ζ))")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "ex4_7_meniscus_verify.png", dpi=110)
    plt.close(fig)
    return worst


def fig_couette_heating():
    U, h, mu, k, rho, cp = 1.0, 1e-3, 1e-3, 0.6, 1000.0, 4182.0
    y = np.linspace(0, h, 201)
    kap = k / (rho * cp)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    for dpdx, a in ((0.0, ax[0]), (-2e3, ax[1])):
        s = ch04.couette_heating(y, U, h, mu, k, dpdx=dpdx)
        a.plot((s["T"] - 293.15) / s["dT_max"], y / h, "k-", lw=2, label="steady (4.60)")
        for f in (0.0, 0.01, 0.05, 0.2):
            T = ch04.couette_heating_transient(y, f * h ** 2 / kap, U, h, mu, k, rho, cp, dpdx=dpdx)
            a.plot((T - 293.15) / s["dT_max"], y / h, label=f"transient t = {f} h²/κ")
        a.set(xlabel="(T − T₀)/ΔT_max", ylabel="y/h", title=f"Couette heating, dp/dx = {dpdx:g} Pa/m")
        a.legend(fontsize=7)
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "couette_heating_verify.png", dpi=110)
    plt.close(fig)


def fig_rankine_bernoulli():
    r = np.linspace(0.01, 3, 300)
    d = ch04.rankine_bernoulli(r, 2 * np.pi, 1.0)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(r, d["B"], "C3", label="B across circles (varies inside the core)")
    ax.plot(r, d["p"] / 1000.0, "C0--", label="p/ρ")
    ax.plot(r, 0.5 * d["u_theta"] ** 2, "C2:", label="½u_θ²")
    ax.axvline(1.0, color="k", lw=0.8)
    ax.set(xlabel="r/σ", ylabel="[m²/s²]", title="Rankine vortex: B constant on each circle, (4.71) not (4.72)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "rankine_bernoulli_verify.png", dpi=110)
    plt.close(fig)


def orders():
    """The convergence studies of the test suite, re-run for the report (orders and pairwise orders)."""
    X1, X2, X3, TT = sp.symbols("x1 x2 x3 t", real=True)
    out = {}
    # continuity flux-divergence stencil
    rexpr = 2 + X1 ** 2 * sp.sin(TT) + sp.exp(X1) * sp.cos(TT)
    uexpr = sp.sin(X1) * (1 + TT)
    rf = sp.lambdify((X1, TT), rexpr)
    uf = sp.lambdify((X1, TT), uexpr)
    ex = sp.lambdify((X1, TT), sp.diff(rexpr * uexpr, X1))
    Xs = np.linspace(-1, 1, 9)[None, :]
    hs = [0.1, 0.05, 0.025, 0.0125]
    errs = [np.max(np.abs(ch04.continuity_terms(lambda X, t: rf(X[0], t), lambda X, t: np.atleast_2d(uf(X[0], t)),
                                                Xs, 0.3, h=h).flux_divergence - ex(Xs[0], 0.3))) for h in hs]
    out["continuity ∇·(ρu) stencil"] = (hs, errs)
    # stream-function velocity
    x, y = np.array([1.4, -1.8, 0.6]), np.array([0.7, 1.2, -1.9])
    ue, ve = ch04.velocity_preset("cylinder", x, y)
    hs2 = [0.08, 0.04, 0.02, 0.01]
    errs = [max(np.max(np.abs(a - b)) for a, b in zip(ch04.velocity_from_streamfunction_2d("cylinder", x, y, h=h), (ue, ve)))
            for h in hs2]
    out["ψ → (u, v) central differences"] = (hs2, errs)
    # mass budget, moving growing sphere, dt
    rho, u = ch04.expanding_flow_fields(a=1.2, rho0=1.0, dim=3)
    cv = ch04.GrowingSphere(R0=0.5, Rdot=0.3, center=(0.2, -0.1, 0.4), U=(0.7, 0.2, -0.5))
    dts = [0.2, 0.1, 0.05, 0.025]
    out["mass budget (4.5) storage FD, dt"] = (dts, [abs(ch04.mass_budget(rho, u, cv, 0.3, dt=d).residual) for d in dts])
    # NS residual Taylor–Green
    u_fn, p_fn, q = ch04.exact_solution_fields("taylor_green", U0=1.0, k=1.0, rho=1.0, mu=0.1)
    X = np.random.default_rng(5).uniform(-2, 2, (2, 8))
    out["NS (4.39b) residual, Taylor–Green"] = (hs, [np.max(np.abs(ch04.ns_incompressible_terms(
        u_fn, p_fn, X, 0.2, 1.0, 0.1, (0, 0), h=h, ht=1e-4).residual)) for h in hs])
    out["NS (4.38) residual (nested stencil)"] = (hs, [np.max(np.abs(ch04.navier_stokes_residual(
        1.0, u_fn, p_fn, X, 0.2, 0.1, 0.0, (0, 0), h=h, ht=1e-4))) for h in hs])
    # rotating basis
    Om = np.array([0.3, -0.5, 0.8])
    E = ch04.rotating_basis(Om, 1.7)
    exact = ch04.basis_rate_exact(Om, E)
    out["de′/dt = Ω × e′ (central difference)"] = (hs, [np.max(np.abs(ch04.basis_rate(Om, 1.7, h=h) - exact)) for h in hs])
    # kinematic BC amplitude
    k = 1.2
    xx = np.linspace(0, 2 * np.pi / k, 17)
    amps = [0.08, 0.04, 0.02, 0.01]
    out["full (4.91) on a linear wave vs a"] = (amps, [np.max(np.abs(ch04.wave_kinematic_residual(xx, 0.3, a_, k)))
                                                       for a_ in amps])
    # cap force ratio
    zetas = [1e-5, 5e-6, 2.5e-6, 1.25e-6]
    out["cap (4.98) exact / small-ζ − 1 vs ζ"] = (zetas, [abs(ch04.cap_surface_tension_force(0.0728, 1e-3, 3e-3, z_) /
                                                              ch04.cap_surface_tension_force(0.0728, 1e-3, 3e-3, z_, exact=False) - 1)
                                                          for z_ in zetas])
    # bore → √(gh)
    dhs = [0.2, 0.1, 0.05, 0.025]
    out["bore U/√(gh) − 1 vs Δh"] = (dhs, [abs(ch04.bore_speed(2.0, 2.0 + d) / np.sqrt(G * 2.0) - 1) for d in dhs])
    # blob heat-equation residual
    T, uu = ch04.gaussian_blob_fields(U=0.2, kappa=1e-3, sigma0=0.05, dim=2)
    Xb = np.random.default_rng(15).uniform(-0.1, 0.1, (2, 12))
    hb = [0.02, 0.01, 0.005, 0.0025]
    out["(4.89) residual, advected Gaussian"] = (hb, [np.max(np.abs(ch04.temperature_equation_residual(T, uu, Xb, 2.0, 1e-3,
                                                                                                        h=h, ht=1e-5)))
                                                      for h in hb])
    # Coriolis small-time deflection
    ts = [0.4, 0.2, 0.1, 0.05]
    out["deflection ut sin Ωt vs Ωut²"] = (ts, [1 - ch04.coriolis_deflection(2.0, 0.3, t_, exact=True) /
                                                ch04.coriolis_deflection(2.0, 0.3, t_) for t_ in ts])
    fig, ax = plt.subplots(figsize=(7, 5))
    for name, (h_, e_) in out.items():
        p = observed_order(h_, e_)
        print(f"{name}: order {p:.3f}  pairwise {[round(v, 3) for v in pairwise_orders(h_, e_)]}")
        ax.loglog(np.asarray(h_) / h_[0], np.asarray(e_) / e_[0], "o-", ms=3, label=f"{name} ({p:.2f})")
    ax.set(xlabel="step / largest step", ylabel="error / largest error", title="Chapter 4 convergence studies")
    ax.legend(fontsize=6)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(OUT / "convergence_verify.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    print("projectile: max |ODE − closed form| =", fig_projectile())
    print("effective gravity: peak at %.2f°, %.2f arcmin" % fig_effective_gravity())
    print("sphere drag C_D(1e4), C_D(1e6) =", fig_sphere_drag())
    print("meniscus: max |x_closed − x_ODE|/δ =", fig_meniscus())
    fig_couette_heating()
    fig_rankine_bernoulli()
    orders()
    print("figures in", OUT)
