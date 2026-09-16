"""Chapter 2 verification figures and report metrics (not collected by pytest).

Reproduces the chapter's key plots with our own code into ``outputs/ch02/verify/`` (git-ignored) and prints the numbers
quoted in ``reports/ch02_verification.md`` (convergence orders, benchmark errors, invariant residuals).

Run: ``.venv/Scripts/python.exe tests/ch02_verify_figures.py``
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

from fluidpy import ch02_cartesian_tensors as ch02  # noqa: E402
from fluidpy.core.style import COLORS, use_style  # noqa: E402
from tools.convergence import observed_order, pairwise_orders  # noqa: E402

OUT = ROOT / "outputs" / "ch02" / "verify"
OUT.mkdir(parents=True, exist_ok=True)
use_style()
C1, C2, C3 = COLORS.get("accent", "C0"), COLORS.get("secondary", "C1"), COLORS.get("tertiary", "C2")


def save(fig, name):
    p = OUT / f"{name}.png"
    fig.savefig(p, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {p.relative_to(ROOT)}")


def fig_polar_components():
    """Fig. 2.3 / Example 2.1: u = (1, 2) resolved in (x1, x2) and in (r, θ) at θ = 30° — with our own numbers."""
    u = np.array([1.0, 2.0])
    th = np.deg2rad(30.0)
    ur, ut = ch02.polar_components(u[0], u[1], th)
    er, et = np.array([np.cos(th), np.sin(th)]), np.array([-np.sin(th), np.cos(th)])
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.annotate("", xy=u, xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", lw=2.5, color="k"))
    ax.text(*(u + 0.05), "u", fontsize=13)
    for vec, col, lab in ((ur * er, C1, f"$u_r$ = {ur:.3f}"), (ut * et, C2, f"$u_\\theta$ = {ut:.3f}")):
        ax.annotate("", xy=vec, xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", lw=2, color=col))
        ax.text(*(vec * 1.05), lab, color=col, fontsize=11)
    ax.plot([0, 3 * er[0]], [0, 3 * er[1]], "--", color=C1, lw=1)
    ax.plot([0, 2.2 * et[0]], [0, 2.2 * et[1]], "--", color=C2, lw=1)
    ax.plot([ur * er[0], u[0]], [ur * er[1], u[1]], ":", color="gray")
    ax.plot([ut * et[0], u[0]], [ut * et[1], u[1]], ":", color="gray")
    ax.axhline(0, color="gray", lw=0.8)
    ax.axvline(0, color="gray", lw=0.8)
    ax.set_aspect("equal")
    ax.set_xlim(-1.2, 3)
    ax.set_ylim(-0.5, 2.7)
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")
    ax.set_title("Example 2.1: $(u_1, u_2) = (1, 2)$ at $\\theta = 30^\\circ$ → $(u_r, u_\\theta)$ = "
                 f"({ur:.3f}, {ut:.3f}); |u| = {np.hypot(ur, ut):.3f}")
    save(fig, "fig2_3_polar_components")


def fig_traction_vs_angle_and_mohr():
    """Example 2.2 generalised: σ_n(φ), τ_s(φ) for τ = [[0, a],[a, 0]] and Mohr's circle (a = 1 Pa)."""
    a = 1.0
    tau = ch02.shear_flow_stress(a)
    phi = np.linspace(0, np.pi, 361)
    sv = ch02.stress_vs_angle(tau, phi)
    c, r = ch02.mohr_circle_2d(tau)
    lam, _ = ch02.principal_axes(tau)
    d30 = ch02.example_2_2(a, np.pi / 6)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))
    ax1.plot(np.rad2deg(phi), sv["sigma_n"], color=C1, lw=2, label="$\\sigma_n = a\\sin 2\\varphi$")
    ax1.plot(np.rad2deg(phi), sv["tau_s"], color=C2, lw=2, label="$\\tau_s = a\\cos 2\\varphi$")
    ax1.plot(np.rad2deg(phi), a * np.sin(2 * phi), "k:", lw=1)
    ax1.plot(np.rad2deg(phi), a * np.cos(2 * phi), "k:", lw=1)
    ax1.axvline(30, color="gray", lw=0.8)
    ax1.plot([30, 30], [d30["sigma_n"], d30["tau_s"]], "o", color="k")
    ax1.set_xlabel("angle of the normal φ [deg]")
    ax1.set_ylabel("stress / a")
    ax1.set_title("Example 2.2: normal and shear stress vs plane angle")
    ax1.legend()
    t = np.linspace(0, 2 * np.pi, 200)
    ax2.plot(c + r * np.cos(t), r * np.sin(t), color="k", lw=1.5, label="Mohr circle (centre I₁/2, radius (λ_max−λ_min)/2)")
    ax2.plot(sv["sigma_n"], sv["tau_s"], ".", color=C1, ms=3, label="(σ_n, τ_s) from (2.15) over φ")
    ax2.plot(lam, [0, 0], "s", color=C2, ms=8, label="eigenvalues ±a")
    ax2.set_aspect("equal")
    ax2.set_xlabel("σ_n / a")
    ax2.set_ylabel("τ_s / a")
    ax2.set_title("Mohr's circle: every plane's (σ_n, τ_s) lies on it")
    ax2.legend(fontsize=8, loc="lower left")
    save(fig, "example2_2_traction_mohr")


def fig_gradient_contours():
    """Fig. 2.7: level curves of φ = x² + y²/4 with ∇φ (our grid operator) perpendicular to them."""
    g = ch02.grid2d(((-2, 2), (-2, 2)), 41)
    phi = g.X ** 2 + g.Y ** 2 / 4
    grad = ch02.gradient(phi, g.h)
    fig, ax = plt.subplots(figsize=(6, 5.4))
    cs = ax.contour(g.X, g.Y, phi, levels=10, cmap="viridis")
    ax.clabel(cs, fmt="%.1f", fontsize=7)
    s = slice(2, None, 4)
    ax.quiver(g.X[s, s], g.Y[s, s], grad[0][s, s], grad[1][s, s], color=C1, scale=40, width=0.004)
    n = np.array([np.cos(0.6), np.sin(0.6)])
    p = np.array([1.0, 0.5])
    ax.annotate("", xy=p + 0.8 * n, xytext=p, arrowprops=dict(arrowstyle="-|>", lw=2, color=C2))
    ax.text(*(p + 0.85 * n), "n", color=C2, fontsize=12)
    gp = np.array([2 * p[0], 0.5 * p[1]])
    ax.annotate("", xy=p + 0.4 * gp, xytext=p, arrowprops=dict(arrowstyle="-|>", lw=2, color="k"))
    ax.text(*(p + 0.42 * gp), "∇φ", fontsize=12)
    ax.set_aspect("equal")
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")
    ax.set_title(f"Fig. 2.7: φ = x₁² + x₂²/4, ∇φ ⊥ level curves; ∂φ/∂n = ∇φ·n = {gp @ n:.3f} at (1, 0.5)")
    save(fig, "fig2_7_gradient")


def fig_principal_axes_deforming_square():
    """Fig. 2.8 / Example 2.4: a square in the pure-strain flow S = [[0, Γ],[Γ, 0]] stretches along b¹ = (1,1)/√2."""
    d = ch02.example_2_4(1.0)
    G = d["S"]  # the symmetric tensor itself as a linear flow (A = 0)
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    for t, col, alpha in ((0.0, "k", 1.0), (0.3, C1, 0.9), (0.6, C2, 0.9)):
        ring = ch02.deform_square(G, t, n_side=30, boundary_only=True)
        ax.plot(*np.c_[ring, ring[:, :1]], color=col, lw=2, alpha=alpha, label=f"t = {t}")
    for k, col in ((0, C3), (1, "purple")):
        b = d["B"][:, k]
        ax.annotate("", xy=1.6 * b, xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", lw=2.5, color=col))
        ax.text(*(1.7 * b), f"b{k + 1}: λ = {d['lam'][k]:+.0f}Γ", color=col, fontsize=11)
    ax.set_aspect("equal")
    ax.set_xlim(-2.4, 2.4)
    ax.set_ylim(-2.4, 2.4)
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")
    ax.set_title(f"Example 2.4: principal axes at {d['angle_deg']:.0f}°, S' = diag({d['S_prime'][0, 0]:.0f}, {d['S_prime'][1, 1]:.0f}) Γ")
    ax.legend(loc="lower right")
    save(fig, "fig2_8_principal_axes")


def fig_convergence():
    """V3 evidence: observed order of the grid operators (one-sided and periodic) and the integral definitions."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))
    orders = {}
    for op, col in (("gradient", C1), ("divergence", C2), ("curl", C3), ("laplacian", "purple")):
        r = ch02.operator_convergence(op, ns=(16, 32, 64, 128))
        rp = ch02.operator_convergence(op, ns=(16, 32, 64, 128), bc="periodic")
        orders[op] = (r["order"], rp["order"], pairwise_orders(r["h"], r["err"]))
        ax1.loglog(r["h"], r["err"], "o-", color=col, label=f"{op} one-sided: p = {r['order']:.3f}")
        ax1.loglog(rp["h"], rp["err"], "s--", color=col, alpha=0.6, label=f"{op} periodic: p = {rp['order']:.3f}")
    h = np.array([0.02, 0.15])
    ax1.loglog(h, 0.5 * h ** 2, "k:", label="∝ h²")
    ax1.set_xlabel("grid spacing h")
    ax1.set_ylabel("max |numeric − exact|")
    ax1.set_title("Grid operators (2.22)–(2.25): second order")
    ax1.legend(fontsize=7)
    iorders = {}
    for kind, col in (("divergence", C1), ("curl", C2), ("gradient", C3), ("curl_component", "purple")):
        r = ch02.integral_definition_convergence(kind)
        iorders[kind] = r["order"]
        ax2.loglog(r["h"], r["err"], "o-", color=col, label=f"{kind}: p = {r['order']:.3f}")
    h = np.array([0.05, 0.4])
    ax2.loglog(h, 0.08 * h ** 2, "k:", label="∝ h²")
    ax2.set_xlabel("box / loop size h")
    ax2.set_ylabel("|integral definition − exact|")
    ax2.set_title("Integral definitions (2.31)–(2.33), (2.35): O(h²)")
    ax2.legend(fontsize=8)
    save(fig, "convergence_orders")
    return orders, iorders


def fig_gauss_faces():
    """Fig. 2.9 / Example 2.5 in 2-D: face fluxes of a box for Q = (x², xy) and the tiled cancellation."""
    Q = ch02.smooth_test_field(2)
    tiles, n = 4, 16
    t = ch02.divergence_theorem_tiled(Q, (-1, 1), tiles, n)
    lhs, rhs, faces = ch02.divergence_theorem_rect2d(Q, (-1, 1), tiles * n, div_fn=Q.div_fn)
    fig, ax = plt.subplots(figsize=(6, 5.6))
    g = ch02.grid2d((-1, 1), 21)
    U = Q(*g.coords)
    ax.streamplot(g.X, g.Y, U[0], U[1], color="lightgray", density=0.8)
    xs = np.linspace(-1, 1, tiles + 1)
    for x in xs:
        ax.plot([x, x], [-1, 1], color=C1, lw=0.8)
        ax.plot([-1, 1], [x, x], color=C1, lw=0.8)
    ax.plot([-1, 1, 1, -1, -1], [-1, -1, 1, 1, -1], color="k", lw=2)
    for key, (x, y) in {"+x": (1.08, 0), "-x": (-1.25, 0), "+y": (0, 1.08), "-y": (0, -1.15)}.items():
        ax.text(x, y, f"{key}: {faces[key]:+.3f}", fontsize=9, ha="center")
    ax.set_aspect("equal")
    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-1.3, 1.3)
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")
    ax.set_title(f"Gauss (2.30) on a tiled box: ∬∇·Q = {lhs:.5f}, ∮n·Q = {rhs:.5f}\n"
                 f"Σ tiles = {t.sum_tiles:.5f}; interior faces = {t.interior:+.1e}")
    save(fig, "fig2_9_gauss_tiled")
    return lhs, rhs, t


def fig_stokes_circulation():
    """Stokes (2.34): circulation vs curl flux for u = b × x on discs of radius R, and the singular vortex."""
    b = 1.0
    f = ch02.solid_body_rotation_field(b, dim=2)
    v = ch02.irrotational_vortex_field(1.0)
    Rs = np.linspace(0.2, 2.0, 10)
    circ, flux, circ_v = [], [], []
    for R in Rs:
        L = ch02.planar_loop([0.0, 0.0], radius=R, n=256)
        D = ch02.planar_disc([0.0, 0.0], radius=R)
        circ.append(ch02.circulation(f, L))
        flux.append(ch02.curl_flux(f, D, curl_fn=f.curl_fn))
        circ_v.append(ch02.circulation(v, L))
    fig, ax = plt.subplots(figsize=(6.5, 4.6))
    ax.plot(Rs, circ, "o", color=C1, label="∮ u·t ds, u = b × x")
    ax.plot(Rs, flux, "-", color=C1, label="∬ (∇×u)·n dA = 2bπR²")
    ax.plot(Rs, 2 * b * np.pi * Rs ** 2, "k:", lw=1)
    ax.plot(Rs, circ_v, "s--", color=C2, label="∮ u·t ds, vortex K/r (core inside): 2πK")
    ax.axhline(2 * np.pi, color=C2, lw=0.8, ls=":")
    ax.set_xlabel("loop radius R")
    ax.set_ylabel("circulation / curl flux")
    ax.set_title("Stokes' theorem (2.34): solid-body rotation obeys it; the singular vortex does not")
    ax.legend(fontsize=8)
    save(fig, "fig2_10_stokes_circulation")
    return np.max(np.abs(np.array(circ) - np.array(flux)))


def main():
    print("figures →", OUT.relative_to(ROOT))
    fig_polar_components()
    fig_traction_vs_angle_and_mohr()
    fig_gradient_contours()
    fig_principal_axes_deforming_square()
    orders, iorders = fig_convergence()
    lhs, rhs, t = fig_gauss_faces()
    mism = fig_stokes_circulation()
    print("\n== metrics for the report ==")
    for op, (p1, p2, pw) in orders.items():
        print(f"{op:16s} one-sided order {p1:.3f} (pairwise {[round(x, 3) for x in pw]}), periodic {p2:.3f}")
    for k, p in iorders.items():
        print(f"integral {k:16s} order {p:.3f}")
    print(f"Gauss tiled: lhs {lhs:.12f} rhs {rhs:.12f} sum_tiles {t.sum_tiles:.12f} interior {t.interior:.2e}")
    print(f"Stokes disc max |circ − flux| = {mism:.2e}")
    F = lambda X, Y, Z: (2 * X, Y ** 2, Z ** 2)  # noqa: E731
    l, r = ch02.divergence_theorem_sphere(F, 1.0, 24, div_fn=lambda X, Y, Z: 2 * (1 + Y + Z))
    print(f"sphere 8π/3: lhs rel {abs(l - 8 * np.pi / 3) / (8 * np.pi / 3):.2e}, rhs rel {abs(r - 8 * np.pi / 3) / (8 * np.pi / 3):.2e}")
    print(f"epsilon_delta_residual = {ch02.epsilon_delta_residual()}")
    tau = np.array([[2.0, 1, 0], [1, 3, 1], [0, 1, 4]])
    print(f"invariants {ch02.invariants(tau)}; bounds {ch02.normal_stress_bounds(tau)}; eig {ch02.principal_axes(tau)[0]}")


if __name__ == "__main__":
    main()
