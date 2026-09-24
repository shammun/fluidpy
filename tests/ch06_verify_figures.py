"""Verification figures and convergence metrics for Chapter 6 (our code only) → outputs/ch06/verify/ (git-ignored).

Run: ``.venv/Scripts/python.exe tests/ch06_verify_figures.py``. Prints the observed orders quoted in
reports/ch06_verification.md and writes one PNG per reproduced figure.
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

from fluidpy import ch06_ideal_flow as ch06  # noqa: E402
from fluidpy.core import conformal as CM  # noqa: E402
from fluidpy.core import panels as PN  # noqa: E402
from fluidpy.core import potential as PF  # noqa: E402
from tools.convergence import observed_order, pairwise_orders  # noqa: E402

OUT = ROOT / "outputs" / "ch06" / "verify"
OUT.mkdir(parents=True, exist_ok=True)
TWO_PI = 2 * np.pi


def grid(xlim, ylim, n=301):
    x = np.linspace(*xlim, n)
    y = np.linspace(*ylim, n)
    return np.meshgrid(x, y, indexing="xy")


def fig_half_body():
    U, m = 1.0, TWO_PI
    fl = PF.half_body(U, m)
    X, Y = grid((-3, 6), (-4, 4))
    P = fl.psi(X, Y)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].contour(X, Y, P, levels=np.linspace(-6, 10, 33), colors="0.4", linewidths=0.6)
    th = np.linspace(0.03, TWO_PI - 0.03, 600)
    xb, yb = ch06.half_body_shape(U, m, th)
    ax[0].plot(xb, yb, "k", lw=2, label="body ψ = m/2")
    ax[0].plot([-1], [0], "o", color="tab:orange", label="stagnation x = −m/2πU")
    ax[0].axhline(np.pi, ls=":", color="tab:blue", label="h_max = m/2U")
    ax[0].set(xlim=(-3, 6), ylim=(-4, 4), aspect="equal", title="Figs. 6.7 analogue: half-body, U = 1, m = 2π",
              xlabel="x [m]", ylabel="y [m]")
    ax[0].legend(fontsize=7, loc="lower right")
    tt = np.linspace(0.05, np.pi, 400)
    ax[1].plot(np.degrees(tt), ch06.half_body_surface_cp(tt), "k")
    t0 = np.degrees(ch06.half_body_cp_zero_angle())
    ax[1].axvline(t0, color="tab:orange", ls="--", label=f"C_p = 0 at {t0:.1f}°")
    ax[1].axhline(0, color="0.6", lw=0.5)
    ax[1].set(xlabel="θ at the source, from +x [deg]", ylabel="surface C_p", title="Fig. 6.8 analogue: C_p along the body")
    ax[1].legend()
    fig.tight_layout()
    fig.savefig(OUT / "fig6_7_6_8_half_body_verify.png", dpi=110)
    plt.close(fig)
    return t0


def fig_cylinder_cp():
    th = np.linspace(0, np.pi, 181)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(np.degrees(np.pi - th), ch06.cylinder_surface_cp(th), "k", label="ideal 1 − 4 sin²θ (6.35)")
    ax.set(xlabel="angle from the forward stagnation point [deg]", ylabel="C_p",
           title="Fig. 6.10 analogue (ideal curve only; real band qualitative)")
    ax.axhline(0, color="0.6", lw=0.5)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "fig6_10_cylinder_cp_verify.png", dpi=110)
    plt.close(fig)


def fig_circulation():
    U, a = 1.0, 1.0
    crit = 4 * np.pi * a * U
    fig, axs = plt.subplots(1, 4, figsize=(15, 4))
    X, Y = grid((-3, 3), (-3, 3), 301)
    for ax, G in zip(axs, (0.0, crit / 2, crit, 1.5 * crit)):
        fl = PF.cylinder(U, a, Gamma_cw=G)
        P = fl.psi(X, Y)
        ax.contour(X, Y, P, levels=np.linspace(-3, 3, 31), colors="0.4", linewidths=0.6)
        ax.add_patch(plt.Circle((0, 0), a, color="0.8"))
        st = ch06.cylinder_stagnation_points(U, a, Gamma_cw=G)
        ax.plot(st.real, st.imag, "o", color="tab:orange")
        L = ch06.surface_pressure_force(fl, R=a, rho=1.0).L
        ax.set(aspect="equal", xlim=(-3, 3), ylim=(-3, 3), title=f"Γ_cw/4πaU = {G / crit:.1f}, L = {L:.2f} (ρUΓ = {G:.2f})")
    fig.suptitle("Fig. 6.12 analogue: cylinder with clockwise circulation (stagnation points in orange)")
    fig.tight_layout()
    fig.savefig(OUT / "fig6_12_cylinder_circulation_verify.png", dpi=100)
    plt.close(fig)


def fig_example_6_1():
    t = np.linspace(0, 80, 801)
    c = ch06.example_6_1(t)
    n = ch06.example_6_1(t[::40], route="numeric")
    tm = ch06.example_6_1_times()
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(t, c["p_origin"], "k", label="closed form")
    ax.plot(t[::40], n["p_origin"], "o", color="tab:blue", label="moving vortex + image, ∂φ/∂t by differences")
    ax.axvline(tm["t_zero"], color="tab:orange", ls="--", label=f"zero at 4πh²/Γ = {tm['t_zero']:.2f} s")
    ax.axvline(tm["t_max"], color="tab:green", ls=":", label=f"max {tm['p_max']:.2f} Pa at {tm['t_max']:.2f} s")
    ax.set(xlabel="t [s]", ylabel="p(0, 0, t) − p∞ [Pa]", title="Example 6.1: water, Γ = 1 m²/s, h = 1 m")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "example6_1_wall_pressure_verify.png", dpi=110)
    plt.close(fig)


def fig_example_6_2():
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    for axi, refine in zip(ax, (1, 8)):
        r = ch06.example_6_2(Q=1.0, method="direct", refine=refine)
        cs = axi.contour(r["X"], r["Y"], r["psi"], levels=np.linspace(0.1, 0.9, 9), colors="k", linewidths=0.8)
        axi.clabel(cs, fontsize=6)
        axi.fill([5, 9, 9, 5], [0, 0, 2, 2], color="0.8")
        axi.set(aspect="equal", title=f"Fig. 6.25 analogue: ψ/Q, grid Δ = {1 / refine:g} m", xlabel="x [m]",
                ylabel="y [m]")
    fig.tight_layout()
    fig.savefig(OUT / "fig6_25_contraction_verify.png", dpi=110)
    plt.close(fig)


def fig_joukowski():
    U, a, b, G = 1.0, 1.2, 1.0, 2.0
    fl = ch06.elliptic_cylinder_flow(U, a, b, Gamma_cw=G)
    X, Y = grid((-4, 4), (-3, 3), 401)
    P = fl.psi(X, Y)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    ax[0].contour(X, Y, P, levels=np.linspace(-3, 3, 31), colors="0.4", linewidths=0.6)
    th = np.linspace(0, TWO_PI, 400)
    zb = CM.joukowski(a * np.exp(1j * th), b)
    ax[0].fill(zb.real, zb.imag, color="0.8")
    ax[0].plot(fl.stagnation.real, fl.stagnation.imag, "o", color="tab:orange")
    ax[0].set(aspect="equal", title="Fig. 6.21 analogue: ellipse (a = 1.2, b = 1), Γ_cw = 2", xlabel="x", ylabel="y")
    Z = X + 1j * Y
    bad = np.abs(CM.joukowski_inverse(Z, b, "principal")) < a
    ax[1].contourf(X, Y, bad, levels=[0.5, 1.5], colors=["#f4a3a3"])
    ax[1].fill(zb.real, zb.imag, color="0.8")
    ax[1].set(aspect="equal", title="where numpy's principal root lands inside the circle (rose)", xlabel="x")
    fig.tight_layout()
    fig.savefig(OUT / "fig6_21_joukowski_verify.png", dpi=100)
    plt.close(fig)


def fig_sphere_airship():
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    Z, R = grid((-3, 3), (0, 2.5), 301)
    sph = PF.sphere(1.0, 1.0)
    ax[0].contour(Z, R, sph.psi(R, Z), levels=np.linspace(0.02, 3, 16), colors="0.4", linewidths=0.6)
    t = np.linspace(0, np.pi, 100)
    ax[0].fill(np.cos(t), np.sin(t), color="0.8")
    ax[0].set(aspect="equal", title="Fig. 6.27 analogue: sphere (meridian plane)", xlabel="z [m]", ylabel="R [m]")
    s = ch06.airship(1.0, 1.0, 1.0)
    Z2, R2 = grid((-1.5, 2.5), (0, 1.2), 301)
    ax[1].contour(Z2, R2, s["psi"](R2, Z2), levels=np.linspace(0.005, 0.6, 16), colors="0.4", linewidths=0.6)
    zc, Rc = s["contour"]
    ax[1].fill(zc, Rc, color="0.8")
    ax[1].plot([0, 1], [0, 0], color="tab:red", lw=3, label="line sink k = Q/a")
    ax[1].plot([0], [0], "o", color="tab:green", label="point source Q")
    ax[1].set(aspect="equal", title=f"Fig. 6.28 analogue: airship, L = {s['length']:.3f} m", xlabel="z [m]",
              ylabel="R [m]")
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "fig6_27_6_28_sphere_airship_verify.png", dpi=110)
    plt.close(fig)


def fig_added_mass():
    a, rho = 0.1, 1000.0
    th = np.linspace(0, np.pi, 181)
    e = np.stack([np.sin(th), 0 * th, np.cos(th)])
    sp_ = PF.moving_sphere_surface_pressure(e, np.array([0, 0, 1.0]), np.array([0, 0, 2.0]), a, rho, split=True)
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(np.degrees(th), sp_["steady"], color="teal", label="speed part (fore–aft symmetric)")
    ax.plot(np.degrees(th), sp_["acceleration"], color="tab:orange", label="acceleration part")
    ax.plot(np.degrees(th), sp_["total"], "k", label="total (6.105)")
    ax.set(xlabel="θ_s from the velocity [deg]", ylabel="p − p∞ [Pa]",
           title="(6.105): u_s = 1 m/s, du_s/dt = 2 m/s², a = 0.1 m, water")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "eq6_105_sphere_pressure_verify.png", dpi=110)
    plt.close(fig)


def convergence():
    out = {}
    eps = np.array([0.2, 0.1, 0.05, 0.025])
    e = ch06.doublet_limit_error(eps)
    out["doublet_limit"] = (observed_order(eps, e), pairwise_orders(eps, e))
    hs = np.array([0.1, 0.05, 0.025, 0.0125])
    f = lambda X, Y: np.sin(X) * np.sin(Y)  # noqa: E731
    err = [abs(PF.laplacian_residual(f, 0.7, 0.4, h) + 2 * f(0.7, 0.4)) for h in hs]
    out["laplacian_residual_9pt"] = (observed_order(hs, err), pairwise_orders(hs, err))
    errs, hh = [], []
    from fluidpy.core import laplace_solvers as LS
    for n in (16, 32, 64, 128):
        xs = np.linspace(0, 1, n + 1)
        X, Y = np.meshgrid(xs, xs, indexing="xy")
        L = LS.laplacian_5pt(np.sin(np.pi * X) * np.sinh(np.pi * Y), 1.0 / n)
        errs.append(abs(L[n // 2, n // 2]))
        hh.append(1.0 / n)
    out["laplacian_5pt"] = (observed_order(hh, errs), pairwise_orders(hh, errs))
    errs, hh = [], []
    for n in (8, 16, 32, 64):
        xs = np.linspace(0, 1, n + 1)
        X, Y = np.meshgrid(xs, xs, indexing="xy")
        ex = np.sin(np.pi * X) * np.sinh(np.pi * Y) / np.sinh(np.pi)
        mask = np.zeros(X.shape, bool)
        mask[1:-1, 1:-1] = True
        psi, _ = LS.solve_laplace(mask, np.where(mask, 0.0, ex), method="direct", dx=1.0 / n)
        errs.append(abs(psi[n // 2, n // 2] - ex[n // 2, n // 2]))
        hh.append(1.0 / n)
    out["solve_laplace_smooth"] = (observed_order(hh, errs), pairwise_orders(hh, errs))
    rs = [2, 4, 8, 16, 32, 64]
    vals = {"far (1,4)": [], "mid (2,2)": [], "near (5.5,2.5)": []}
    pts = {"far (1,4)": (1.0, 4.0), "mid (2,2)": (2.0, 2.0), "near (5.5,2.5)": (5.5, 2.5)}
    for r in rs:
        res = ch06.example_6_2(method="direct", refine=r)
        for k, (x, y) in pts.items():
            vals[k].append(res["psi"][int(round(y * r)), int(round(x * r))])
    for k in vals:
        d = np.abs(np.diff(vals[k]))
        h = 1.0 / np.array(rs[:-1])
        out[f"example_6_2 {k}"] = (observed_order(h[-3:], d[-3:]), pairwise_orders(h, d))
    Ns = np.array([32, 64, 128, 256])
    e = [ch06.panel_cp_error(int(N)) for N in Ns]
    out["panels_ellipse_cp"] = (observed_order(1 / Ns, e), pairwise_orders(1 / Ns, e))
    errs = []
    for N in Ns:
        th = np.linspace(0, TWO_PI, int(N), endpoint=False) + np.pi / N
        r = PN.source_panels(np.cos(th), np.sin(th))
        P = (np.array([2.0, 0.0, -1.5]), np.array([0.5, 1.8, -1.0]))
        u, v = PN.panel_velocity(r, *P)
        ue, ve = PF.cylinder(1.0, 1.0).velocity(*P)
        errs.append(np.max(np.hypot(u - ue, v - ve)))
    out["panels_off_body_velocity"] = (observed_order(1 / Ns, errs), pairwise_orders(1 / Ns, errs))
    fl = PF.cylinder(10.0, 0.1, Gamma_cw=2.0)
    Nb = np.array([16, 32, 64, 128])
    th = lambda N: 0.3 * np.exp(1j * np.linspace(0, TWO_PI, int(N), endpoint=False))  # noqa: E731
    err = [abs(PF.blasius_force(fl, contour=th(N), rho=1.2).L - 24.0) for N in Nb]
    out["blasius_polygon_circle"] = (observed_order(1 / Nb, err), pairwise_orders(1 / Nb, err))
    ell = ch06.elliptic_cylinder_flow(1.0, 1.2, 1.0, Gamma_cw=2.0)
    errs = []
    for N in (8, 16, 32, 64):
        s = np.linspace(-3.0, 3.0, N, endpoint=False)
        sq = np.concatenate([s - 3j, 3.0 + 1j * s, -s + 3j, -3.0 - 1j * s])
        F = PF.blasius_force(ell, contour=sq, rho=1.2)
        errs.append(np.hypot(F.D, F.L - 2.4))
    out["blasius_polygon_square_ellipse"] = (observed_order(1 / np.array([8, 16, 32, 64]), errs),
                                             pairwise_orders(1 / np.array([8, 16, 32, 64]), errs))
    t = np.array([0.0, 5.0, 12.0, 20.0])
    c = ch06.example_6_1(t)["p_origin"]
    dts = [0.4, 0.2, 0.1, 0.05]
    err = [np.max(np.abs(ch06.example_6_1(t, route="numeric", dt=d)["p_origin"] - c)) for d in dts]
    out["example_6_1_dt"] = (observed_order(dts, err), pairwise_orders(dts, err))
    Ns = np.array([16, 32, 64])
    fits = [ch06.axisym_body_fit("ellipsoid", int(N)) for N in Ns]
    ke = [q["k_error"] for q in fits]
    be = [q["body_psi_error"] for q in fits]
    out["axial_ellipsoid_k"] = (observed_order(1 / Ns, ke), pairwise_orders(1 / Ns, ke))
    out["axial_ellipsoid_body_psi"] = (observed_order(1 / Ns, be), pairwise_orders(1 / Ns, be))
    out["axial_ellipsoid_cond"] = [q["solution"]["cond"] for q in fits]
    xe = np.array([10.0, 100.0, 1000.0])
    D = [abs(ch06.half_body_net_force(x_end=x)["D"]) for x in xe]
    out["half_body_net_force_vs_1/x_end"] = (observed_order(1 / xe, D), pairwise_orders(1 / xe, D))
    d = 1.0
    P = (np.array([0.7, 1.3, 0.2]), np.array([0.9, -0.4, 1.5]))
    ref = PF.Doublet3D(d).phi(*P)
    eps = np.array([0.1, 0.05, 0.025, 0.0125])
    err = []
    for ee in eps:
        Qe = d / (2 * ee)
        pair = PF.PointSource3D(Qe, -ee).phi(*P) + PF.PointSource3D(-Qe, ee).phi(*P)
        err.append(np.max(np.abs(pair - ref)) / np.max(np.abs(ref)))
    out["doublet3d_limit"] = (observed_order(eps, err), pairwise_orders(eps, err))
    for k, v in out.items():
        print(f"{k}: {v}")
    return out


if __name__ == "__main__":
    print("half-body C_p zero angle [deg]:", fig_half_body())
    fig_cylinder_cp()
    fig_circulation()
    fig_example_6_1()
    fig_example_6_2()
    fig_joukowski()
    fig_sphere_airship()
    fig_added_mass()
    convergence()
