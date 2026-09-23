"""Verification figures and convergence metrics for Chapter 5 (our code only) → outputs/ch05/verify/ (git-ignored).

Run: ``.venv/Scripts/python.exe tests/ch05_verify_figures.py``. Prints the observed orders quoted in
reports/ch05_verification.md and writes one PNG per reproduced figure.
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

from fluidpy import ch05_vorticity_dynamics as ch05  # noqa: E402
from fluidpy.core import vortices as VX  # noqa: E402
from tools.convergence import observed_order, pairwise_orders  # noqa: E402

OUT = ROOT / "outputs" / "ch05" / "verify"
OUT.mkdir(parents=True, exist_ok=True)
G = 9.81


def fig_isobars():
    r = np.linspace(0.0, 0.3, 200)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    for dp in (-0.02, 0.0, 0.02):
        ax[0].plot(r, ch05.isobar_height(r, 10.0, dp, "solid"), label=f"Δp/ρg = {dp:+.2f} m")
    t = ch05.rotating_tank_free_surface(0.3, 0.1, 5.0)
    ax[0].plot(r, t["z_vertex"] + 25.0 * r ** 2 / (2 * G), "k--", label="free surface (volume kept)")
    ax[0].set(title="Fig. 5.2 analogue: tank at 5 rad/s (ω = 10 s⁻¹)", xlabel="r [m]", ylabel="z [m]")
    ax[0].legend(fontsize=8)
    rr = np.linspace(0.02, 0.6, 300)
    for dp in (-0.05, -0.2):
        ax[1].plot(rr, ch05.isobar_height(rr, 1.0, dp, "line"), label=f"line vortex, Δp/ρg = {dp}")
        ax[1].plot(rr, ch05.isobar_height(rr, 1.0, dp, "rankine", a=0.1), "--", label=f"Rankine a = 0.1 m, {dp}")
    ax[1].set(title="Fig. 5.3 analogue: funnel (Γ = 1 m²/s)", xlabel="r [m]", ylabel="z [m]", ylim=(-1.0, 0.1))
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "fig5_2_5_3_isobars_verify.png", dpi=110)
    plt.close(fig)
    return t


def fig_tube_budget():
    z = np.linspace(0.0, 2.0, 41)
    s = [ch05.gaussian_tube_section(zz) for zz in z]
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].plot(z, [q["flux"] for q in s], "k-", label="flux Γ(1 − e⁻¹) [m²/s]")
    ax[0].plot(z, [q["area"] * 20 for q in s], "C0-", label="area × 20 [m²]")
    ax[0].plot(z, [q["mean_omega"] / 100 for q in s], "C3-", label="mean ω / 100 [1/s]")
    ax[0].set(title="(5.4): narrowing Gaussian tube (Fig. 5.1 analogue)", xlabel="z [m]")
    ax[0].legend(fontsize=8)
    b = ch05.tube_flux_budget("gaussian_tube", 0.0, 1.0, 0.1)
    bb = ch05.tube_flux_budget("broken", 0.0, 1.0, 0.1)
    lab = ["lower", "side", "upper", "total"]
    x = np.arange(4)
    ax[1].bar(x - 0.2, list(b), 0.4, label="vorticity field (∇·ω = 0)")
    ax[1].bar(x + 0.2, list(bb), 0.4, label="'broken' field (∇·ω ≠ 0)")
    ax[1].set_xticks(x, lab)
    ax[1].set(title="Gauss budget of (5.4) on a tube piece", ylabel="outward flux [m²/s]")
    ax[1].axhline(0, color="k", lw=0.6)
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig5_1_tube_budget_verify.png", dpi=110)
    plt.close(fig)
    return b, bb


def fig_kelvin():
    s = ch05.kelvin_scenario("cellular")
    t = np.linspace(0.0, s["t_end"], 41)
    Gt, loops = ch05.material_circulation(s["u"], s["pts0"], t, return_loops=True)
    lo = ch05.kelvin_scenario("lamb_oseen")
    tl = np.linspace(0.0, 60.0, 13)
    Gl = ch05.material_circulation(lo["u"], lo["pts0"], tl)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    for i, c in ((0, "C0"), (20, "C1"), (40, "C3")):
        L = loops[i]
        ax[0].plot(np.append(L[0], L[0][0]), np.append(L[1], L[1][0]), c, lw=1, label=f"t = {t[i]:.1f} s")
    ax[0].set(title="Fig. 5.4 analogue: material loop in the cellular flow", xlabel="x [m]", ylabel="y [m]",
              aspect="equal")
    ax[0].legend(fontsize=8)
    ax[1].plot(t, (Gt - Gt[0]) / Gt[0], "C0.-", label="cellular (inviscid): (Γ − Γ₀)/Γ₀")
    ax2 = ax[1].twinx()
    ax2.plot(tl, Gl * 1e3, "C3o", label="Lamb–Oseen circle (viscous)")
    ax2.plot(tl, lo["Gamma_exact"](tl) * 1e3, "C3-", lw=0.8)
    ax[1].set(title="Kelvin (5.8) vs (5.11)", xlabel="t [s]", ylabel="relative change of Γ (cellular)")
    ax2.set_ylabel("Γ [10⁻³ m²/s] (Lamb–Oseen)")
    ax[1].legend(loc="upper left", fontsize=8)
    ax2.legend(loc="lower left", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig5_4_kelvin_verify.png", dpi=110)
    plt.close(fig)
    return np.max(np.abs(Gt - Gt[0])), ch05.loop_length(loops[-1]) / ch05.loop_length(loops[0])


def fig_baroclinic():
    Rs = np.array([0.4, 0.2, 0.1, 0.05, 0.025])
    err = [abs(ch05.pressure_torque_on_element(radius=R, grad_rho=(200.0, 50.0), grad_p=(300.0, -9810.0))["ratio"] - 1)
           for R in Rs]
    tilt = np.linspace(-np.pi / 2, np.pi / 2, 31)
    sp_ = [ch05.baroclinic_element_scenario(a, 5.0, 0.05)["spin_up"] for a in tilt]
    lf = ch05.lock_exchange_fields(nx=257, ny=65)
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    ax[0].loglog(Rs, err, "o-")
    ax[0].set(title=f"Fig. 5.6: torque route → (5.28), order {observed_order(Rs, err):.2f}", xlabel="disc radius R [m]",
              ylabel="|spin-up/(∇ρ×∇p/ρ²) − 1|")
    ax[1].plot(np.degrees(tilt), sp_, "C1-")
    ax[1].set(title="baroclinic spin-up vs isopycnal tilt", xlabel="tilt [°]", ylabel="Dω/Dt [1/s²]")
    im = ax[2].pcolormesh(lf["X"], lf["Y"], lf["baroclinic_z"], shading="auto", cmap="Oranges")
    fig.colorbar(im, ax=ax[2], label="(∇ρ×∇p)_z/ρ² [1/s²]")
    ax[2].set(title="Fig. 5.5 analogue: lock exchange t = 0⁺", xlabel="x [m]", ylabel="y [m]")
    fig.tight_layout()
    fig.savefig(OUT / "fig5_5_5_6_baroclinic_verify.png", dpi=110)
    plt.close(fig)
    return observed_order(Rs, err), pairwise_orders(Rs, err)


def fig_biot_savart():
    F = ch05.gaussian_tube_fields()
    b = F["bounds"]
    nodes, w = ch05.cylinder_quadrature(b["radius"], b["z0"], b["z1"], 40, 64, 40)
    r = np.linspace(0.02, 0.8, 30)
    P = np.stack([r, 0 * r, 0 * r])
    ub = ch05.biot_savart_volume(F["omega"], P, nodes, w, eps=0.0)[1]
    up = ch05.velocity_from_curl_omega(F["curl_omega"], P, nodes, w)[1]
    um = ch05.velocity_from_curl_omega(F["curl_omega"], P, nodes, w, sign=-1.0)[1]
    rr = r[r >= 0.3]
    ref = [F["u_theta_reference"](q) for q in rr]
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.plot(r, [F["u_theta_infinite"](q) for q in r], "k-", lw=0.8, label="infinite tube Γ(1−e^{−r²/σ²})/2πr")
    ax.plot(r, ub, "C0o", ms=4, label="Biot–Savart (5.16), L = 4 m")
    ax.plot(r, up, "C1x", ms=5, label="(5.14) with +1/(4π)")
    ax.plot(r, um, "C3s", ms=3, label="(5.14) as printed, −1/(4π)")
    ax.plot(rr, ref, "C2^", ms=4, label="1-D reference (finite tube)")
    ax.axhline(0, color="k", lw=0.5)
    ax.set(title="C07: velocity of a Gaussian tube (Γ = 1 m²/s, σ = 0.1 m)", xlabel="r [m]", ylabel="u_θ [m/s]")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "biot_savart_tube_verify.png", dpi=110)
    plt.close(fig)


def fig_pairs_wall():
    fig, ax = plt.subplots(1, 3, figsize=(13, 4))
    for name, c in (("equal_pair", "C0"), ("unequal_pair", "C1")):
        pr = ch05.point_vortex_preset(name)
        X = ch05.point_vortex_evolve(pr["xv"], pr["Gamma"], np.linspace(0, pr["t_end"] / 2, 300))
        for k in range(2):
            ax[0].plot(X[:, 0, k], X[:, 1, k], c, lw=1 + k)
        g = ch05.centre_of_vorticity(pr["xv"], pr["Gamma"])
        ax[0].plot(*g, c + "*", ms=12, label=f"{name}: G")
    ax[0].set(title="Fig. 5.11 analogue: pairs orbit G", aspect="equal", xlabel="x [m]", ylabel="y [m]")
    ax[0].legend(fontsize=7)
    op = ch05.point_vortex_preset("opposite_pair")
    X = ch05.point_vortex_evolve(op["xv"], op["Gamma"], np.linspace(0, 10, 50))
    ax[1].plot(X[:, 0, 0], X[:, 1, 0], "C0.-", label="+Γ")
    ax[1].plot(X[:, 0, 1], X[:, 1, 1], "C3.-", label="−Γ")
    kb = ch05.point_vortex_preset("knife_bucket")
    Xk = ch05.point_vortex_evolve(kb["xv"], kb["Gamma"], np.linspace(0, kb["t_end"], 300), boundary="circle", **kb["bp"])
    th = np.linspace(0, 2 * np.pi, 200)
    ax[1].plot(0.3 + np.cos(th), np.sin(th), "k-", lw=0.6)
    ax[1].plot(Xk[:, 0, 0], Xk[:, 1, 0], "C3-", lw=0.8)
    ax[1].plot(Xk[:, 0, 1], Xk[:, 1, 1], "C0-", lw=0.8)
    ax[1].set(title="Figs. 5.12–5.13 analogues", aspect="equal", xlabel="x [m]", ylabel="y [m]")
    ax[1].legend(fontsize=7)
    xv, Gm = ch05.wall_image_system(np.array([[0.0], [0.5]]), [1.0])
    xw = np.stack([np.linspace(-2, 2, 21), np.zeros(21)])
    uw = ch05.point_vortex_velocity(xw, xv, Gm)
    ax[2].quiver(xw[0], xw[1], uw[0], uw[1], color="C0", scale=4)
    ax[2].plot(0, 0.5, "C0o")
    ax[2].plot(0, -0.5, "C3o", mfc="none")
    ax[2].axhline(0, color="k")
    ax[2].set(title=f"Fig. 5.14 analogue: max |v| on wall = {np.max(np.abs(uw[1])):.1e}", aspect="equal",
              xlabel="x [m]", ylabel="y [m]", ylim=(-0.8, 0.8))
    fig.tight_layout()
    fig.savefig(OUT / "fig5_11_5_14_point_vortices_verify.png", dpi=110)
    plt.close(fig)


def fig_rings():
    tw = np.linspace(0.0, 20.0, 201)
    w = ch05.ring_dynamics([dict(R=1.0, z=0.0, Gamma=1.0, a=0.1)], tw, wall_z=3.0)
    t = np.linspace(0.0, 30.0, 601)
    d = ch05.ring_dynamics([dict(R=1.0, z=0.0, Gamma=1.0, a=0.1), dict(R=1.0, z=0.5, Gamma=1.0, a=0.1)], t)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].plot(w["R"][:, 0], w["z"][:, 0], "C0-", label="ring")
    ax[0].plot(w["R"][:, 0], 6.0 - w["z"][:, 0], "C3--", label="image ring")
    ax[0].axhline(3.0, color="k")
    ax[0].set(title="Fig. 5.15 analogue: ring toward a wall", xlabel="R [m]", ylabel="z [m]")
    ax[0].legend(fontsize=8)
    ax[1].plot(d["z"][:, 0], d["R"][:, 0], "C0-", lw=0.8, label="ring 1")
    ax[1].plot(d["z"][:, 1], d["R"][:, 1], "C1-", lw=0.8, label="ring 2")
    ax[1].set(title=f"leap-frogging (impulse drift {np.ptp(d['impulse']) / d['impulse'][0]:.1e})", xlabel="z [m]",
              ylabel="R [m]")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig5_15_rings_verify.png", dpi=110)
    plt.close(fig)
    return w["R"][0, 0], w["R"][-1, 0]


def fig_sheet():
    y = np.linspace(-0.1, 0.1, 801)
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    for N, c in ((10, "C0"), (100, "C1"), (1000, "C2")):
        ax.plot(ch05.discrete_sheet_u(0.0, y, 2.0, N), y, c, lw=1, label=f"row of N = {N}")
    ax.plot(ch05.continuous_sheet_velocity(0.0, y, 2.0)[0], y, "k--", lw=1, label="continuous sheet")
    d = ch05.discrete_sheet_convergence(2.0)
    ax.set(title=f"Fig. 5.16 analogue: γ = 2 m/s; L1 order {observed_order(1 / d['N'], d['l1_error']):.3f} in 1/N",
           xlabel="u(0, y) [m/s]", ylabel="y [m]")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig5_16_sheet_verify.png", dpi=110)
    plt.close(fig)
    return d


def fig_column_burgers():
    x = np.linspace(-3e5, 3e5, 301)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    for depth, c in (("ridge", "C3"), ("trough", "C0"), ("slope", "C2")):
        col = ch05.column_over_slope(x, depth=depth)
        ax[0].plot(x / 1e3, col["zeta"] * 1e5, c, label=f"{depth}: ζ [10⁻⁵ s⁻¹]")
    ax[0].axhline(0, color="k", lw=0.5)
    ax[0].set(title="Fig. 5.10 analogue: (ζ + f)/h conserved, 45° N", xlabel="x [km]")
    ax[0].legend(fontsize=8)
    R = np.linspace(0.0, 6e-3, 200)
    bb = ch05.burgers_balance(R)
    ax[1].plot(R * 1e3, bb["advective"], "C0", label="advective u_R ∂ω/∂R")
    ax[1].plot(R * 1e3, bb["stretching"], "C4", label="stretching αω")
    ax[1].plot(R * 1e3, bb["diffusion"], "C3", label="diffusion ν∇²ω")
    ax[1].plot(R * 1e3, bb["residual"], "k--", label="adv − str − diff")
    ax[1].axvline(2.0, color="gray", lw=0.6)
    ax[1].set(title="Burgers balance (D18), core 2 mm", xlabel="R [mm]", ylabel="[1/s²]")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig5_10_column_burgers_verify.png", dpi=110)
    plt.close(fig)


def fig_hill_outside():
    A, a = 1.0, 1.0
    Rg, zg = np.meshgrid(np.linspace(0.05, 2.0, 24), np.linspace(-2.0, 2.0, 40))
    uR, uz, _ = VX.hill_spherical_vortex(Rg, zg, A, a)
    psi = VX.hill_stream_function(Rg, zg, A, a)
    th = np.linspace(0.01, np.pi - 0.01, 200)
    Ro, zo = 1.000001 * np.sin(th), 1.000001 * np.cos(th)
    uo = VX.hill_spherical_vortex(Ro, zo, A, a)
    ui = VX.hill_spherical_vortex(0.999999 * np.sin(th), 0.999999 * np.cos(th), A, a)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.4))
    ax[0].contour(Rg, zg, psi, 21, colors="gray", linewidths=0.6)
    ax[0].quiver(Rg, zg, uR, uz, color="C0", scale=6)
    ax[0].plot(np.sin(th), np.cos(th), "k-")
    ax[0].set(title="Hill: ψ contours (grey) vs coded velocity (blue)", xlabel="R", ylabel="z", aspect="equal")
    ax[1].plot(np.degrees(th), uo[0] * np.sin(th) + uo[1] * np.cos(th), "C3", label="normal u just outside (code)")
    ax[1].plot(np.degrees(th), ui[0] * np.sin(th) + ui[1] * np.cos(th), "C0", label="normal u just inside (code)")
    ax[1].axhline(0, color="k", lw=0.5)
    ax[1].set(title="F2: fluid crosses the sphere outside", xlabel="θ [°]", ylabel="u·e_r on r = a [m/s]")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "hill_outside_F2_verify.png", dpi=110)
    plt.close(fig)
    return float(np.max(np.abs(uo[0] * np.sin(th) + uo[1] * np.cos(th))))


def fig_core_edge():
    r = np.linspace(0.05, 0.3, 501)
    r = np.sort(np.append(r, 0.1))
    d = [ch05.vortex_pressure_scenario("cylinder", q, omega=2.0, a=0.1) for q in r]
    tq = np.array([x["torque"] for x in d])
    fn = np.array([x["net_viscous_force"] for x in d])
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].plot(r, tq * 1e5, "C0.-", ms=2)
    ax[0].axhline(-2e-3 * np.pi * 0.01 * 2.0 * 1e5, color="C3", lw=0.8, label="−2μΓ (torque_per_length)")
    ax[0].set(title="F1: E2 'cylinder' torque 2πr²σ_rθ", xlabel="r [m]", ylabel="torque [10⁻⁵ N m/m]")
    ax[0].legend(fontsize=8)
    ax[1].plot(r, fn, "C3.-", ms=2)
    ax[1].set(title="F1: reported net viscous force (status 'no net viscous force')", xlabel="r [m]",
              ylabel="[N/m³]")
    fig.tight_layout()
    fig.savefig(OUT / "core_edge_F1_verify.png", dpi=110)
    plt.close(fig)


def orders():
    res = {}
    u = ch05.abc_flow(1.0, 0.7, 0.4)
    P = np.random.default_rng(5).uniform(-3, 3, (3, 30))
    hs = [0.1, 0.05, 0.025, 0.0125]
    e = [np.max(np.abs(ch05.vorticity_field(u, h)(P) - u(P))) for h in hs]
    res["curl stencil (ABC)"] = (hs, e)
    ub = ch05.burgers_vortex_field(1e-3, 1.0, 1e-6)
    hs = [4e-5, 2e-5, 1e-5, 5e-6]
    e = [float(np.max(np.abs(ch05.vorticity_terms(ub, np.array([1e-3, 5e-4, 2e-4]), 0.0, 1e-6, h=h).residual)))
         for h in hs]
    res["(5.13) residual (Burgers)"] = (hs, e)
    Ms = [16, 32, 64, 128]
    e = [abs(ch05.filament_velocity_preset("ring", 0, 0, 0, M=M, component=2) - 1.0) for M in Ms]
    res["polygon ring (5.17)"] = ([1 / M for M in Ms], e)
    pts = np.array([[0.2, 0.0, 0.3], [0.8, 0.0, -0.2]]).T
    uR, uz = ch05.ring_ring_velocity(pts[0], pts[2], 0.5, 0.0, 1.0)
    Ms = [32, 64, 128, 256]
    e = []
    for M in Ms:
        uf = ch05.filament_velocity(pts, ch05.filament_preset("ring", M, R=0.5), 1.0)
        e.append(max(np.max(np.abs(uf[0] - uR)), np.max(np.abs(uf[2] - uz))))
    res["ring-ring elliptic vs polygon"] = ([1 / M for M in Ms], e)
    Rs = [0.4, 0.2, 0.1, 0.05]
    e = [abs(ch05.pressure_torque_on_element(radius=R, grad_rho=(200.0, 50.0), grad_p=(300.0, -9810.0))["ratio"] - 1)
         for R in Rs]
    res["pressure torque → (5.28)"] = (Rs, e)
    d = ch05.discrete_sheet_convergence(2.0)
    res["discrete sheet L1"] = (list(1.0 / d["N"]), list(d["l1_error"]))
    u5 = lambda x, t=0.0: np.stack([np.sin(np.asarray(x)[1]) * (1 + 0.3 * t) + 0.2 * np.asarray(x)[0],  # noqa: E731
                                    0.5 * np.sin(np.asarray(x)[0]) - 0.2 * np.asarray(x)[1]])
    pts0 = ch05.circle_loop_points((0.4, 0.1), 0.6, 256)
    lp = ch05.material_loop(u5, pts0, [0.0, 0.7])[-1]
    kr = ch05.kelvin_rate_terms(u5, lp, 0.7, h=1e-5, ht=1e-5).total
    dts, e = [0.2, 0.1, 0.05, 0.025], []
    for dt in dts:
        Gm = ch05.material_circulation(u5, pts0, [0.0, 0.7 - dt, 0.7, 0.7 + dt])
        e.append(abs((Gm[3] - Gm[1]) / (2 * dt) - kr))
    res["(5.9) vs dΓ/dt of a material loop"] = (dts, e)
    lo = lambda nr: abs(ch05.vortex_tube_strength(ch05.lamb_oseen_field(1.0, 1e-3, 2.5), None,  # noqa: E731
                                                  ch05.IT.planar_disc((0, 0), None, 0.1, nr, 64),
                                                  omega=lambda x, t=0.0: np.array([100 / np.pi * np.exp(
                                                      -(np.asarray(x)[0] ** 2 + np.asarray(x)[1] ** 2) / 0.01)])).flux
                        - (1 - np.exp(-1)))
    nrs = [32, 64, 128, 256]
    res["tube strength flux route (disc)"] = ([0.1 / n for n in nrs], [lo(n) for n in nrs])
    fig, ax = plt.subplots(figsize=(7, 5))
    for k, (h, e) in res.items():
        p = observed_order(h, e)
        print(f"{k:40s} observed {p:.3f}  pairwise {[round(v, 3) for v in pairwise_orders(h, e)]}")
        ax.loglog(h, e, "o-", label=f"{k}: {p:.2f}")
    ax.set(title="Ch. 5 convergence studies", xlabel="h (step, 1/M, R, dt, 1/N)", ylabel="error")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "convergence_verify.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    t = fig_isobars()
    print("tank", {k: round(v, 6) for k, v in t.items() if isinstance(v, float)})
    print("tube budget", fig_tube_budget())
    print("kelvin cellular max |ΔΓ|, stretch", fig_kelvin())
    print("baroclinic order", fig_baroclinic())
    fig_biot_savart()
    fig_pairs_wall()
    print("ring toward wall R", fig_rings())
    d = fig_sheet()
    print("sheet L1", d["l1_error"])
    fig_column_burgers()
    print("Hill outside max normal velocity", fig_hill_outside())
    fig_core_edge()
    orders()
    print("figures in", OUT)
