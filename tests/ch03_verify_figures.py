"""Verifier's figure reproductions and convergence numbers for Chapter 3 (Kinematics).

Draws the chapter's key plots with *our* code (never book crops) into ``outputs/ch03/verify/`` (git-ignored) and prints
the observed orders quoted in ``reports/ch03_verification.md``.

Run: ``.venv/Scripts/python.exe tests/ch03_verify_figures.py``
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from fluidpy import ch03_kinematics as ch03  # noqa: E402
from tools.convergence import observed_order, pairwise_orders  # noqa: E402

OUT = ROOT / "outputs" / "ch03" / "verify"
OUT.mkdir(parents=True, exist_ok=True)


def fig_flow_lines(tp: float = math.pi / 6):
    """Fig. 3.7 (Ex. 3.1): closed forms (lines) vs numerical streamline / path line / streak line (markers)."""
    ex = ch03.example_3_1(tp)
    f = ex["field"]
    sl = ch03.streamline(f, [0.0, 0.0], tp, 2.0, n=41)
    pl = ch03.pathline(f, [0.0, 0.0], tp, np.linspace(tp, tp + 2 * np.pi, 41))
    sk = ch03.streakline(f, [0.0, 0.0], tp, np.linspace(tp - 2 * np.pi, tp, 41))
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(*ex["streamline"], color="teal", lw=2, label="streamline $y = x\\tan\\omega t'$ (closed form)")
    ax.plot(*ex["pathline"], color="orange", lw=2, label="path line circle (closed form)")
    ax.plot(*ex["streakline"], color="crimson", lw=2, label="streak line circle (closed form)")
    ax.plot(*sl, "o", ms=3, color="teal", mfc="none", label="streamline (ODE, frozen t')")
    ax.plot(*pl, "s", ms=3, color="orange", mfc="none", label="path line (ODE)")
    ax.plot(*sk, "^", ms=3, color="crimson", mfc="none", label="streak line (ODE ensemble)")
    ax.plot(0, 0, "k*", ms=10, label="dye port / origin")
    ax.set_aspect("equal")
    ax.set_xlabel("x / $\\xi_o$")
    ax.set_ylabel("y / $\\xi_o$")
    ax.set_title(f"Ex. 3.1 flow lines at $\\omega t'$ = {math.degrees(tp):.0f}°")
    ax.legend(fontsize=7, loc="lower left")
    fig.savefig(OUT / "fig3_7_flow_lines_verify.png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    c = ex["path_center"]
    err = max(np.max(np.abs(np.hypot(pl[0] - c[0], pl[1] - c[1]) - 1)),
              np.max(np.abs(np.hypot(sk[0] - ex["streak_center"][0], sk[1] - ex["streak_center"][1]) - 1)))
    print(f"Fig 3.7: max radius error of numerical path/streak circles = {err:.2e}")


def fig_shear_elements(gamma: float = 1.0, t: float = 0.3):
    """Fig. 3.14: aligned element ABCD shears, 45° element PQRS stretches; both turn clockwise."""
    G = ch03.velocity_gradient_preset("simple_shear", gamma)
    M = ch03.linear_flow_map(G, t)
    sq = np.array([[-0.5, 0.5, 0.5, -0.5, -0.5], [-0.5, -0.5, 0.5, 0.5, -0.5]])
    rot = np.array([[math.cos(math.pi / 4), -math.sin(math.pi / 4)], [math.sin(math.pi / 4), math.cos(math.pi / 4)]])
    dia = rot @ sq * 0.7071
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.5))
    for ax, P, name in ((axs[0], sq, "ABCD (aligned)"), (axs[1], dia, "PQRS (45°)")):
        ax.plot(*P, color="gray", lw=1, label="t = 0")
        ax.plot(*(M @ P), color="purple", lw=2, label=f"t = {t} s")
        ax.set_aspect("equal")
        ax.set_title(f"{name}, γ = {gamma} 1/s")
        ax.set_xlabel("x₁ [m]")
        ax.set_ylabel("x₂ [m]")
        ax.legend(fontsize=8)
        sides = np.linalg.norm(np.diff(M @ P, axis=1), axis=0)
        print(f"Fig 3.14 {name}: side lengths at t = {t}: {np.round(sides, 4)}")
    fig.savefig(OUT / "fig3_14_shear_elements_verify.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


def fig_vortex_profiles():
    """(3.28)–(3.29): u_θ and ω_z of the Rankine and Gaussian vortices, peaks at σ and 1.1209σ."""
    Gm, s = 2 * np.pi, 1.0
    r = np.linspace(1e-3, 4, 800)
    ur, wr = ch03.rankine_vortex(r, Gm, s)
    ug, wg = ch03.gaussian_vortex(r, Gm, s)
    rs = ch03.gaussian_vortex_max_radius(s)
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    axs[0].plot(r, ur, color="teal", label="Rankine $u_\\theta$")
    axs[0].plot(r, ug, color="orange", label="Gaussian $u_\\theta$")
    axs[0].plot(r, Gm / (2 * np.pi * r), ":", color="gray", label="line vortex Γ/2πr")
    axs[0].axvline(1.0, color="teal", ls="--", lw=0.8)
    axs[0].axvline(rs, color="orange", ls="--", lw=0.8, label=f"r* = {rs:.4f} σ")
    axs[0].set_ylim(0, 1.2)
    axs[0].set_xlabel("r / σ")
    axs[0].set_ylabel("$u_\\theta$ [m/s] (Γ = 2π m²/s, σ = 1 m)")
    axs[0].legend(fontsize=8)
    axs[1].plot(r, wr, color="teal", label="Rankine $\\omega_z$")
    axs[1].plot(r, wg, color="orange", label="Gaussian $\\omega_z$")
    axs[1].set_xlabel("r / σ")
    axs[1].set_ylabel("$\\omega_z$ [1/s]")
    axs[1].legend(fontsize=8)
    fig.savefig(OUT / "vortex_profiles_verify.png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    print(f"Vortices: Rankine peak {ur.max():.4f} m/s at r = {r[np.argmax(ur)]:.4f}; Gaussian peak {ug.max():.4f} at "
          f"r = {r[np.argmax(ug)]:.4f} (r* = {rs:.10f})")


def fig_elements_in_vortices(t_end: float = math.pi / 2):
    """Figs. 3.15–3.16: a small square carried by solid-body rotation turns; in the line vortex it keeps its orientation."""
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    half = 0.08
    sq = np.array([[-half, half, half, -half, -half], [-half, -half, half, half, -half]]) + np.array([[1.0], [0.0]])
    for ax, kind, fld in ((axs[0], "solid body ω₀ = 1", ch03.vortex_velocity_field(lambda q: ch03.solid_body_rotation(q, 1.0))),
                          (axs[1], "line vortex B = 1", ch03.vortex_velocity_field(lambda q: ch03.line_vortex(q, 1.0)))):
        for tt, col in ((0.0, "gray"), (t_end / 2, "teal"), (t_end, "purple")):
            P = np.stack([ch03.pathline(fld, sq[:, k], 0.0, tt) for k in range(sq.shape[1])], axis=1) if tt > 0 else sq
            ax.plot(*P, color=col, lw=2, label=f"t = {tt:.2f} s")
            edge = P[:, 1] - P[:, 0]  # initially horizontal (angle 0)
            side = P[:, 3] - P[:, 0]  # initially vertical (angle 90°)
            a_b = math.degrees(math.atan2(edge[1], edge[0]))
            a_s = (math.degrees(math.atan2(side[1], side[0])) - 90.0 + 180.0) % 360.0 - 180.0
            print(f"{kind}: t = {tt:.2f}, bottom edge {a_b:7.2f}°, left edge {a_s:7.2f}°, pair average {(a_b + a_s) / 2:7.2f}°")
        for size in (0.08, 0.02):  # the book's statement: instantaneous pair-average turning rate of threads through a point
            tt = 1e-3
            cross = np.array([[1 - size, 1 + size, 1.0, 1.0], [0.0, 0.0, -size, size]])  # two threads through (1, 0)
            P = np.stack([ch03.pathline(fld, cross[:, k], 0.0, tt, rtol=1e-12, atol=1e-14) for k in range(4)], axis=1)
            e1, e2 = P[:, 1] - P[:, 0], P[:, 3] - P[:, 2]
            rate = 0.5 * (math.atan2(e1[1], e1[0]) + math.atan2(e2[1], e2[0]) - math.pi / 2) / tt
            print(f"{kind}: pair-average turning rate (threads through (1,0), half-length {size}, dt = {tt}): {rate:+.6f} rad/s")
        th = np.linspace(0, 2 * np.pi, 200)
        ax.plot(np.cos(th), np.sin(th), ":", color="gray", lw=0.8)
        ax.plot(0, 0, "k+")
        ax.set_aspect("equal")
        ax.set_title(kind)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.legend(fontsize=8)
    fig.savefig(OUT / "fig3_15_16_elements_verify.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


def convergence_studies():
    """Every V3 study of the suite in one log–log figure; prints observed orders."""
    rows = []
    F = lambda x, t: np.sin(np.asarray(x)[0]) * np.exp(0.3 * np.asarray(x)[1]) * np.cos(t) + np.asarray(x)[0] * np.asarray(x)[1] * t  # noqa: E731
    u = lambda x, t: np.stack([np.sin(np.asarray(x)[0]) * np.cos(np.asarray(x)[1]) + 0.3 * t + 0.2 * np.asarray(x)[0],  # noqa: E731
                               -np.cos(np.asarray(x)[0]) * np.sin(np.asarray(x)[1]) + 0.2 * np.asarray(x)[0] * t ** 2 + 0.1 * np.asarray(x)[1] ** 2])
    x0, t0 = np.array([0.3, -0.4]), 0.6
    ref = ch03.material_derivative(F, u, x0, t0, h=1e-3)
    ref = ch03.material_derivative(F, u, x0, t0, h=2.5e-4)
    hs = np.array([0.1, 0.05, 0.025, 0.0125])
    # Richardson-free: compare against a much finer stencil
    errs = [abs(ch03.material_derivative(F, u, x0, t0, h=h) - ref) for h in hs]
    rows.append(("material derivative stencil (h)", hs, errs, 2))
    G_ex = ch03.velocity_gradient_at(u, x0, t0, h=2.5e-4)
    errs = [np.max(np.abs(ch03.velocity_gradient_at(u, x0, t0, h=h) - G_ex)) for h in hs]
    rows.append(("velocity gradient stencil (h)", hs, errs, 2))
    ur = lambda r, th: r ** 2 * np.sin(th)  # noqa: E731
    ut = lambda r, th: r * np.cos(th)  # noqa: E731
    wex = lambda r, th: 2 * np.cos(th) - r * np.cos(th)  # noqa: E731
    ds = np.array([0.2, 0.1, 0.05, 0.025])
    errs = []
    for d in ds:
        lg = ch03.annular_sector_circulation(ut, 1.1 - d / 2, d, d / 1.1, u_r=ur, theta0=0.6 - d / 2.2, legs=True)
        errs.append(abs(lg["total"] / lg["area"] - wex(1.1, 0.6)))
    rows.append(("sector Γ/area → ω_z (size)", ds, errs, 2))
    cv = ch03.GrowingSphere(R0=0.9, Rdot=0.15, center=(0.1, 0, 0), U=(0.2, 0, 0))
    Fs = lambda x, t: (1 + 0.3 * np.asarray(x)[0]) * np.sin(1 + 2 * t)  # noqa: E731
    dFs = lambda x, t: 2 * (1 + 0.3 * np.asarray(x)[0]) * np.cos(1 + 2 * t)  # noqa: E731
    dts = np.array([0.1, 0.05, 0.025, 0.0125])
    rows.append(("RTT: FD of (3.31) vs (3.35) (Δt)", dts, [ch03.rtt_check(Fs, dFs, cv, 0.4, dt=dt).abs_error for dt in dts], 2))
    F3 = lambda x, t: (1 + 0.3 * np.asarray(x)[0] + 0.2 * np.asarray(x)[1] ** 2 + 0.1 * np.asarray(x)[0] * np.asarray(x)[2]) * (1 + 0.5 * t + 0.2 * t ** 2) + 0.05 * np.asarray(x)[2] ** 3 * np.sin(t)  # noqa: E731,E501
    dF3 = lambda x, t: (1 + 0.3 * np.asarray(x)[0] + 0.2 * np.asarray(x)[1] ** 2 + 0.1 * np.asarray(x)[0] * np.asarray(x)[2]) * (0.5 + 0.4 * t) + 0.05 * np.asarray(x)[2] ** 3 * np.cos(t)  # noqa: E731,E501
    dts2 = np.array([0.08, 0.04, 0.02, 0.01])
    sw = [ch03.swept_terms_sphere(0.9, 0.15, F3, dF3, 0.4, dt) for dt in dts2]
    tot = ch03.reynolds_transport(F3, dF3, ch03.GrowingSphere(R0=0.9, Rdot=0.15), 0.4).total
    rows.append(("(3.32) T4 = ∫_ΔV Δt ∂F/∂t (Δt)", dts2, [abs(s["T4"]) for s in sw], 2))
    rows.append(("(3.34) sliver error (Δt)", dts2, [abs(s["sliver_error"]) for s in sw], 2))
    rows.append(("(3.31) one-sided quotient − total (Δt)", dts2, [abs(s["lhs"] - tot) for s in sw], 1))
    Fl = lambda x, t: np.sin(x - t)  # noqa: E731
    ex = ch03.leibniz_example(0.8, "wave")["exact"]
    dtl = np.array([0.2, 0.1, 0.05, 0.025])
    rows.append(("Leibniz measured rate (Δt)", dtl, [abs(ch03.leibniz_check(Fl, lambda t: t / 2, lambda t: 2 + t, 0.8, dt=dt) - ex) for dt in dtl], 2))
    Gr = np.array([[0.3, -0.7], [0.9, -0.2]])
    dtm = np.array([1e-2, 5e-3, 2.5e-3, 1.25e-3])
    rows.append(("tracked-segment stretch → n·S·n (dt)", dtm,
                 [abs(ch03.measured_strain_rates(Gr, [0.8, 0.6], dt)["stretch"] - ch03.linear_strain_rate(Gr, [0.8, 0.6])) for dt in dtm], 1))
    Gs = np.array([[0.7, 0.3], [0.3, -0.4]])
    ts = np.array([0.08, 0.04, 0.02, 0.01])
    rows.append(("ellipse first order vs exact (t)", ts,
                 [np.max(np.abs(ch03.strain_ellipse_axes(Gs, t)[0] - ch03.strain_ellipse_axes(Gs, t, "exact")[0])) for t in ts], 2))
    fig, ax = plt.subplots(figsize=(8, 6))
    for name, h, e, p in rows:
        po = observed_order(h, e)
        print(f"{name:45s} design {p}  observed {po:.3f}  pairwise {np.round(pairwise_orders(h, e), 3)}")
        ax.loglog(h, e, "o-", label=f"{name}: {po:.2f}")
    ax.set_xlabel("step (h, size, Δt or t) [SI]")
    ax.set_ylabel("error")
    ax.set_title("Chapter 3 convergence studies (observed orders in legend)")
    ax.legend(fontsize=7)
    fig.savefig(OUT / "convergence_verify.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    fig_flow_lines()
    fig_shear_elements()
    fig_vortex_profiles()
    fig_elements_in_vortices()
    convergence_studies()
    print("figures in", OUT)
