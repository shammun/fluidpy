"""Verification figures and convergence metrics for Chapter 9 (our code only) → outputs/ch09/verify/ (git-ignored).

Run: ``.venv/Scripts/python.exe tests/ch09_verify_figures.py``. Prints the observed orders and key numbers quoted in
reports/ch09_verification.md and writes one PNG per reproduced figure. No book figure is copied: every curve is computed
by ``fluidpy`` with our own parameters (Fig. numbers are the chapter's, the pictures are ours).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fluidpy import ch09_boundary_layers as ch09  # noqa: E402
from fluidpy.core import bluff_body as BB  # noqa: E402
from fluidpy.core import boundary_layer as BL  # noqa: E402
from fluidpy.core import jets as JET  # noqa: E402
from fluidpy.core import similarity as SIM  # noqa: E402
from tools.convergence import observed_order, pairwise_orders  # noqa: E402

OUT = ROOT / "outputs" / "ch09" / "verify"
OUT.mkdir(parents=True, exist_ok=True)
REF = ROOT / "reference" / "ch09" / "benchmarks.json"


def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.png", dpi=110)
    plt.close(fig)
    print("wrote", OUT / f"{name}.png")


def fig_9_5_9_6_blasius():
    c = BL.blasius_constants()
    eta = np.linspace(0, 8, 801)
    f, fp, fpp = BL.blasius_profile(eta)
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.2))
    ax[0].plot(fp, eta, label="f′(η) = u/U")
    ax[0].plot(fpp, eta, "--", label="f″(η)")
    ax[0].axhline(c["eta99"], color="gray", ls=":", label=f"η99 = {c['eta99']:.3f}")
    ax[0].axhline(4.93, color="r", ls=":", label="printed 4.93 (slip R2)")
    ax[0].axhline(c["delta_star"], color="tab:orange", ls="-.", label=f"δ*/δ = {c['delta_star']:.4f}")
    ax[0].axhline(c["theta"], color="tab:green", ls="-.", label=f"θ/δ = {c['theta']:.4f}")
    ax[0].set_xlabel("u/U, f″"), ax[0].set_ylabel("η = y√(U/νx)"), ax[0].set_title("Fig. 9.5 analogue: Blasius profile"), ax[0].legend(fontsize=7)
    U, nu = 1.0, 1.5e-5
    for x, ls in ((0.1, "-"), (0.5, "--"), (2.0, ":")):
        y = np.linspace(0, 0.03, 400)
        ax[1].plot(BL.blasius_fields(x, y, U, nu)["u"], y * 1e3, ls, label=f"x = {x} m")
    ax[1].set_xlabel("u [m/s]"), ax[1].set_ylabel("y [mm]"), ax[1].set_title("dimensional profiles stretch ∝ √x"), ax[1].legend(fontsize=8)
    v = 0.5 * (eta * fp - f)
    ax[2].plot(eta, v)
    ax[2].axhline(c["v_inf"], color="gray", ls=":", label=f"δ*/2 = {c['v_inf']:.4f}")
    ax[2].set_xlabel("η"), ax[2].set_ylabel("v√Re_x/U"), ax[2].set_title("Fig. 9.6 analogue"), ax[2].legend()
    save(fig, "fig_9_5_9_6_blasius")


def fig_9_7_falkner_skan():
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    for n in (4, 1, 1 / 3, 1 / 9, 0, -0.0654, -0.0904):
        d = BL.falkner_skan(n, eta_max=16.0, n=1601)
        ax[0].plot(0.5 * np.sqrt(n + 1) * d["eta"], d["fp"], label=f"n = {n:.4g}")
    ax[0].set_xlim(0, 6), ax[0].set_xlabel("½√(n+1) η"), ax[0].set_ylabel("f′ = u/U_e"), ax[0].set_title("Fig. 9.7 analogue"), ax[0].legend(fontsize=7)
    ms = np.r_[np.linspace(-0.0904, 0, 60), np.linspace(0, 4, 60)[1:]]
    k = [BL.falkner_skan_state(m, eta_max=16.0)["fpp0"] for m in ms]
    ax[1].plot(ms, k, label="attached branch (ours)")
    if REF.exists():
        R = json.loads(REF.read_text(encoding="utf-8"))["belden_falkner_skan"]["points"]
        for p in R:
            m = p["beta"] / (2 - p["beta"])
            ax[1].plot(m, np.sqrt((m + 1) / 2) * p["kappa"], "ks", ms=6, label="Belden et al. (published)" if p is R[0] else None)
    sep = BL.falkner_skan_separation()
    ax[1].plot(sep["m_sep"], 0, "ro", label=f"fold m = {sep['m_sep']:.5f}")
    for mm in (-0.02, -0.04, -0.06):
        d = BL.falkner_skan(mm, branch="reversed")
        ax[1].plot(mm, d["fpp0"], "b^", ms=5, label="reversed branch (ours)" if mm == -0.02 else None)
    ax[1].set_xlabel("n"), ax[1].set_ylabel("f″(0)"), ax[1].set_xlim(-0.1, 1.0), ax[1].legend(fontsize=7), ax[1].axhline(0, color="k", lw=0.5)
    ax[1].set_title("wall shear vs n, with published points")
    lam, L = [], []
    tb = BL.thwaites_closure_table(60)
    ax[2].plot(tb["lam"], tb["L"], "o-", ms=3, label="L(λ) exact FS closure")
    ll = np.linspace(-0.09, 0.12, 50)
    ax[2].plot(ll, 0.45 - 6 * ll, "k--", label="0.45 − 6λ")
    ax[2].set_xlabel("λ"), ax[2].set_ylabel("L"), ax[2].legend(), ax[2].set_title("Fig. 9.8 analogue")
    save(fig, "fig_9_7_9_8_falkner_skan")
    print("fold m_sep =", sep["m_sep"], "beta_sep =", sep["beta_sep"])


def fig_thwaites():
    nu = 1e-6
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    xi = np.linspace(0, 0.3, 601)
    of = BL.outer_flow("diffuser", U1=1.0, L=1.0)
    r = BL.thwaites(xi, of, 1e-5, closure="white", stop_at_separation=False)
    ax[0].plot(xi, r["lam"], label="Thwaites λ(x)")
    ax[0].plot(xi, -(0.45 / 4) * ((1 + xi) ** 4 - 1), "k--", label="closed form −(0.45/4)[(1+x/L)⁴−1]")
    ax[0].axhline(-0.09, color="r", ls=":", label="λ = −0.09"), ax[0].axhline(BL.LAMBDA_SEP_FS, color="m", ls=":", label="λ = −0.0681 (exact FS zero shear)")
    ax[0].set_xlabel("x/L"), ax[0].set_ylabel("λ"), ax[0].set_title("Example 9.2 (diffuser)"), ax[0].legend(fontsize=7)
    ns = np.array([-0.09, -0.05, 0.0, 0.1, 1 / 3, 1.0, 4.0])
    errs = []
    for n in ns:
        st = BL.falkner_skan_state(n)
        x = 0.7
        of2 = BL.outer_flow("wedge", n=n, a=1.0)
        x0 = 1e-6 * x if n < 1.4 else 0.05 * x
        xs = x0 + (x - x0) * np.linspace(0, 1, 4001) ** 2
        th0 = np.sqrt(0.45 * nu * xs[0] ** (1 - n) / (1.0 * (5 * n + 1)))
        rr = BL.thwaites(xs, of2, nu, theta0=th0, stop_at_separation=False)
        errs.append(100 * (rr["theta"][-1] / (st["I_theta"] * np.sqrt(nu * x / x ** n)) - 1))
    ax[1].plot(ns, errs, "o-"), ax[1].axhspan(-3, 3, color="g", alpha=0.15, label="book's ±3 % (favourable)"), ax[1].axhline(8, color="r", ls=":"), ax[1].axhline(-8, color="r", ls=":", label="test bound ±8 %")
    ax[1].set_xlabel("n (U_e = a xⁿ)"), ax[1].set_ylabel("θ error vs exact FS [%]"), ax[1].legend(fontsize=7), ax[1].set_title("Thwaites accuracy on the FS family")
    phi = np.linspace(2, 170, 300)
    ax[2].plot(phi, BL.thwaites_cylinder_closed_form(np.deg2rad(phi)), label="closed form λ(φ)")
    ax[2].plot(phi[::15], BL.thwaites_cylinder(phi[::15]), "k.", label="quadrature")
    ax[2].axhline(-0.09, color="r", ls=":"), ax[2].axvline(BL.thwaites_cylinder_separation(-0.09), color="r", ls=":", label=f"φ_sep = {BL.thwaites_cylinder_separation(-0.09):.2f}°")
    ax[2].set_ylim(-0.3, 0.1), ax[2].set_xlabel("φ from the forward stagnation point [deg]"), ax[2].set_ylabel("λ"), ax[2].legend(fontsize=7), ax[2].set_title("D11: cylinder")
    save(fig, "fig_thwaites")
    print("Thwaites theta errors [%] at n =", dict(zip(np.round(ns, 3), np.round(errs, 2))))


def fig_marching():
    nu = 1e-3
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    U, x1 = 1.0, 1.0
    ex = float(BL.blasius_wall_shear(x1, U, 1.0, nu))
    of = BL.outer_flow("flat", U=U)
    hs, errs = [], []
    for ny in (50, 100, 200, 400, 800):
        r = BL.march_boundary_layer(of, np.linspace(0.1, x1, 161), nu, ny=ny)
        hs.append(1 / ny), errs.append(abs(r["tau0"][-1] / ex - 1))
    print("march plate, Δσ order:", observed_order(hs[:3], errs[:3]), pairwise_orders(hs, errs))
    ax[0].loglog(hs, errs, "o-", label="plate (B = 0)")
    st = BL.falkner_skan_state(0.5)
    exw = st["fpp0"] * nu * np.sqrt(1.0 / nu)
    hs2, errs2 = [], []
    for ny in (50, 100, 200, 400, 800):
        t = BL.march_boundary_layer(BL.outer_flow("wedge", n=0.5, a=1.0), np.linspace(0.1, 1.0, 161), nu, ny=ny)["tau0"][-1]
        hs2.append(1 / ny), errs2.append(abs(t / exw - 1))
    print("march wedge n=0.5, Δσ order:", observed_order(hs2, errs2), pairwise_orders(hs2, errs2))
    ax[0].loglog(hs2, errs2, "s-", label="wedge n = 0.5 (dp/dx ≠ 0)")
    h = np.array(hs)
    ax[0].loglog(h, errs[2] * (h / h[2]) ** 2, "k:", label="slope 2"), ax[0].loglog(h, errs2[2] * (h / h[2]), "k--", label="slope 1")
    ax[0].set_xlabel("Δσ = 1/ny"), ax[0].set_ylabel("|τ₀/τ₀,exact − 1|"), ax[0].legend(fontsize=7), ax[0].set_title("marching: wall-shear error vs Δσ")
    for n, mk in ((0.5, "s-"), (-0.05, "^-")):
        of_n = BL.outer_flow("wedge", n=n, a=1.0)
        ref = BL.march_boundary_layer(of_n, np.linspace(0.1, 1, 1281), nu, ny=200)["tau0"][-1]
        gr = (20, 40, 80, 160) if n < 0 else (10, 20, 40, 80)
        hh, ee = [], []
        for nx in gr:
            hh.append(0.9 / nx), ee.append(abs(BL.march_boundary_layer(of_n, np.linspace(0.1, 1, nx + 1), nu, ny=200)["tau0"][-1] - ref))
        print(f"march BDF2 n={n}, Δx order:", observed_order(hh, ee), pairwise_orders(hh, ee))
        ax[1].loglog(hh, ee, mk, label=f"BDF2, n = {n}")
        ref1 = BL.march_boundary_layer(of_n, np.linspace(0.1, 1, 1281), nu, ny=200, order=1)["tau0"][-1]
        e1 = [abs(BL.march_boundary_layer(of_n, np.linspace(0.1, 1, nx + 1), nu, ny=200, order=1)["tau0"][-1] - ref1) for nx in gr]
        print(f"march BE n={n}, Δx order:", observed_order(hh, e1))
        ax[1].loglog(hh, e1, mk.replace("-", "--"), label=f"backward Euler, n = {n}")
    ax[1].set_xlabel("Δx"), ax[1].set_ylabel("|τ₀ − τ₀,ref|"), ax[1].legend(fontsize=7), ax[1].set_title("marching: Δx convergence")
    x = np.linspace(0.1, 50.0, 1500)
    for n in (-0.089, -0.095, -0.12):
        Ue0 = 0.1 ** n
        uin = (lambda y, Ue0=Ue0: Ue0 * BL.blasius_profile(np.asarray(y) / np.sqrt(nu * 0.1 / Ue0))[1]) if n < -0.0904 else None
        r = BL.march_boundary_layer(BL.outer_flow("wedge", n=n, a=1.0), x, nu, u_inlet=uin)
        ax[2].semilogx(r["x"], r["tau0"] / (nu * np.sqrt(1.0)), label=f"n = {n}: " + (f"separates at x = {r['x_sep']:.2f}" if r["separated"] else "attached to x = 50"))
    ax[2].set_xlabel("x"), ax[2].set_ylabel("τ₀ (scaled)"), ax[2].legend(fontsize=7), ax[2].set_title("wedge flows around the fold n = −0.0904")
    save(fig, "fig_marching")


def fig_bluff():
    fig, ax = plt.subplots(2, 3, figsize=(14, 7.6))
    b = np.linspace(0.05, 1.0, 200)
    ax[0, 0].plot(b, BB.karman_street_growth_closed(b), label="closed form")
    bb = np.array([0.1, 0.2, 0.25, 0.3, 0.35, 0.5, 1.0])
    ax[0, 0].plot(bb, [np.max(BB.karman_street_spectrum_periodic(v, n_pairs=16).real) for v in bb], "ro", label="finite periodic cell (16 pairs)")
    ax[0, 0].plot(bb, [BB.karman_street_growth(v) for v in bb], "k.", label="max over 61 wavenumbers")
    ax[0, 0].plot(bb, [BB.karman_street_growth(v, offset=0.0) for v in bb], "g--", label="facing rows (offset 0): π/4")
    ax[0, 0].axvline(BB.karman_street_ratio(), color="gray", ls=":", label=f"b/a = {BB.karman_street_ratio():.5f}")
    ax[0, 0].set_xlabel("b/a"), ax[0, 0].set_ylabel("growth [Γ/a²]"), ax[0, 0].legend(fontsize=7), ax[0, 0].set_title("Kármán street: growth rate")
    for bv, c in ((0.2, "tab:red"), (BB.karman_street_ratio(), "tab:green"), (0.5, "tab:blue")):
        sp = BB.karman_street_spectrum(bv)
        ax[0, 1].plot(sp.real, sp.imag, "o", color=c, label=f"b/a = {bv:.3f}")
    ax[0, 1].axvline(0, color="k", lw=0.5), ax[0, 1].set_xlabel("Re σ"), ax[0, 1].set_ylabel("Im σ"), ax[0, 1].legend(fontsize=7), ax[0, 1].set_title("spectrum (k = π/a)")
    ph = np.linspace(0, 180, 361)
    ax[0, 2].plot(ph, BB.cp_ideal_cylinder(ph), "k--", label="ideal 1 − 4 sin²φ")
    ax[0, 2].plot(ph, BB.separated_cp(ph, 82.0, cp_base=-1.2), label="model, φs = 82°, Cb = −1.2")
    ax[0, 2].plot(ph, BB.separated_cp(ph, 125.0, cp_base=-0.6), label="model, φs = 125°, Cb = −0.6")
    ax[0, 2].set_ylim(-3.2, 1.2), ax[0, 2].set_xlabel("φ from the forward stagnation point [deg]"), ax[0, 2].set_ylabel("C_p"), ax[0, 2].legend(fontsize=7), ax[0, 2].set_title("Fig. 9.20 analogue (qualitative)")
    Re = np.logspace(-1, 7, 400)
    ax[1, 0].loglog(Re, BB.cylinder_cd_schematic(Re), label="cylinder (qualitative schematic)")
    Rs = np.logspace(-1, 6, 400)
    ax[1, 0].loglog(Rs, SIM.sphere_drag_coefficient(Rs), label="sphere: Morrison correlation (V5)")
    ax[1, 0].loglog(Rs[Rs < 5], ch09.stokes_drag_coefficient(Rs[Rs < 5]), "k:", label="Stokes 24/Re"), ax[1, 0].loglog(Rs[Rs < 5], ch09.oseen_drag_coefficient(Rs[Rs < 5]), "k--", label="Oseen")
    ax[1, 0].set_xlabel("Re (diameter)"), ax[1, 0].set_ylabel("C_D"), ax[1, 0].legend(fontsize=7), ax[1, 0].set_title("Figs. 9.21–9.22 analogues")
    R = np.logspace(4, 7, 200)
    ax[1, 1].semilogx(R, [BB.drag_crisis_state(r)["phi_sep_deg"] for r in R], label="smooth"), ax[1, 1].semilogx(R, [BB.drag_crisis_state(r, rough=True)["phi_sep_deg"] for r in R], "--", label="rough (Re_cr/3)")
    ax[1, 1].set_xlabel("Re"), ax[1, 1].set_ylabel("model φ_sep [deg]"), ax[1, 1].legend(), ax[1, 1].set_title("separation moves aft (illustrative)")
    ax[1, 2].loglog(R, [BB.drag_crisis_state(r)["cd_model"] for r in R], label="model pressure drag C_D,p (illustrative)")
    ax[1, 2].loglog(R, BL.plate_drag_coefficient(np.maximum(R, 1e5), "laminar"), "--", label="plate laminar 1.328/√Re_L")
    ax[1, 2].loglog(R, BL.plate_drag_coefficient(np.maximum(R, 1e5), "turbulent"), ":", label="plate turbulent 0.074 Re^{-1/5}")
    ax[1, 2].loglog(R, BL.plate_drag_coefficient(R, "mixed"), "-.", label="plate mixed (Re_tr = 5e5)")
    ax[1, 2].set_xlabel("Re_L"), ax[1, 2].set_ylabel("C_D"), ax[1, 2].legend(fontsize=7), ax[1, 2].set_title("Fig. 9.11 analogue (plate)")
    save(fig, "fig_bluff_and_street")


def fig_jets():
    fig, ax = plt.subplots(1, 4, figsize=(17, 4.2))
    e = np.linspace(-12, 12, 481)
    p = JET.free_jet_profile(e)
    ax[0].plot(e, p["fp"], label="sech²(η/√6) (closed form)")
    d = JET.free_jet_ode_solve()
    ax[0].plot(d["eta"], d["fp"], "r--", label=f"solve_bvp (max err {d['max_err']:.1e})")
    ax[0].plot(-d["eta"], d["fp"], "r--")
    ax[0].set_xlabel("η"), ax[0].set_ylabel("u/u₀"), ax[0].legend(fontsize=8), ax[0].set_title("free jet profile")
    xs = np.logspace(-2, 1, 30)
    J, rho, nu = 1.0, 1.2, 1.5e-5
    ax[1].loglog(xs, [float(JET.free_jet_centreline(x, J, rho, nu)) for x in xs], label="u₀ ∝ x^{-1/3}")
    ax[1].loglog(xs, [float(JET.free_jet_halfwidth(x, J, rho, nu)) * 1e3 for x in xs], label="h99 [mm] ∝ x^{2/3}")
    ax[1].loglog(xs, [float(JET.free_jet_mass_flux(x, J, rho, nu)) for x in xs], label="ṁ ∝ x^{1/3}")
    ax[1].set_xlabel("x [m]"), ax[1].legend(fontsize=8), ax[1].set_title("exponents from constant J")
    w = JET.wall_jet_ode_solve(1.0 / 72.0, eta_max=120.0, n=3001)
    ax[2].plot(w["eta"], w["fp"] * 12, label="f′ × 12 (IVP, f_∞ = 1)"), ax[2].plot(w["eta"], w["f"], label="f (IVP)")
    pr = JET.wall_jet_profile(w["eta"][::100], 1.0)
    ax[2].plot(w["eta"][::100], pr["fp"] * 12, "k.", label="implicit (9.83) by brentq"), ax[2].set_xlim(0, 40)
    b = JET.wall_jet_ode_solve(1.0 / 72.0, eta_max=120.0, coeff=1.0, n=3001)
    ax[2].plot(b["eta"], b["f"], "r:", label=f"printed ODE: f_∞ = {b['f_inf']:.3f}")
    ax[2].set_xlabel("η"), ax[2].legend(fontsize=7), ax[2].set_title("wall jet: 4f‴+ff″+2f′² = 0")
    xw = np.array([0.25, 0.5, 1.0, 2.0, 4.0])
    inv, mom = [], []
    for x in xw:
        delta = float(JET.wall_jet(x, 0.0, 0.7, 1.3, 1e-3)["delta"])
        y = np.linspace(0, 110 * delta / 1.3, 4001)
        f = JET.wall_jet(x, y, 0.7, 1.3, 1e-3)
        inv.append(JET.wall_jet_invariant(y, f["u"])), mom.append(np.trapezoid(f["u"] ** 2, y))
    ax[3].loglog(xw, np.array(inv) / inv[0], "o-", label="∫u(∫u²) — invariant (9.80)"), ax[3].loglog(xw, np.array(mom) / mom[0], "s-", label="∫u²dy ∝ x^{-1/4}")
    ax[3].set_xlabel("x"), ax[3].legend(fontsize=8), ax[3].set_title("what the wall jet conserves")
    save(fig, "fig_jets")
    print("wall jet default eta_max=40: f_inf =", JET.wall_jet_ode_solve(1 / 72)["f_inf"], "err_vs_9_83 =", JET.wall_jet_ode_solve(1 / 72)["err_vs_9_83"])


if __name__ == "__main__":
    fig_9_5_9_6_blasius()
    fig_9_7_falkner_skan()
    fig_thwaites()
    fig_marching()
    fig_bluff()
    fig_jets()
