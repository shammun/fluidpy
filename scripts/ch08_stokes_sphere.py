"""§8.6 (C12–C15, Figs. 8.17, 8.19): low-Re scaling coefficients, Stokes' sphere — streamlines in the body and
fluid frames with the ideal-flow sphere alongside, surface pressure (8.50) and stresses, the drag 6πμaU with its ⅓/⅔
split by Gauss–Legendre quadrature, the Stokes residual, and the far-field ratio inertia/viscous ~ Re·r/a.

Run: ``.venv/Scripts/python.exe scripts/ch08_stokes_sphere.py --no-show``
Figures → outputs/ch08/c13_stokes_sphere.png, c14_surface_stresses.png.
"""
from __future__ import annotations

import numpy as np

from ch08_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch08_laminar_flow as ch08


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    sc = ch08.low_re_scaling_sympy()
    print(f"dynamic pressure scale ×Re: {sc['dynamic_times_Re']} → Re → 0: {sc['dynamic_limit']}")
    print(f"viscous pressure scale ×Re: {sc['viscous_times_Re']} → Re → 0: {sc['viscous_limit']}")
    a, U, mu = 1.0, 1.0, 1.0
    xs = np.linspace(-4, 4, 321)
    X, Y = np.meshgrid(xs, xs)
    R = np.hypot(X, Y)
    TH = np.arctan2(np.abs(Y), X)
    sgn = np.where(Y >= 0, 1.0, -1.0)
    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    lev = np.linspace(0.05, 6, 16) ** 1.3 * 0.3
    for axi, frame, title in ((ax[0], "body", "Stokes, body frame (8.48)"),
                              (ax[1], "fluid", "Stokes, fluid frame (Fig. 8.19)")):
        psi = sgn * ch08.stokes_sphere_streamfunction(R, TH, U, a, frame)
        axi.contour(X, Y, psi, np.concatenate([-lev[::-1], lev]), colors=COLORS["accent"], linewidths=0.8)
        axi.add_patch(plt.Circle((0, 0), a, color=COLORS["muted"]))
        axi.set_aspect("equal")
        axi.set_title(title)
        axi.set_xlabel("x/a")
    ideal = sgn * 0.5 * U * R ** 2 * np.sin(TH) ** 2 * (1 - a ** 3 / np.maximum(R, a) ** 3)
    ideal = np.where(R >= a, ideal, np.nan)
    ax[2].contour(X, Y, ideal, np.concatenate([-lev[::-1], lev]), colors=COLORS["teal"], linewidths=0.8)
    ax[2].contour(X, Y, sgn * ch08.stokes_sphere_streamfunction(R, TH, U, a), [0.5, 1.5], colors=COLORS["accent"],
                  linewidths=1.5, linestyles="--")
    ax[2].add_patch(plt.Circle((0, 0), a, color=COLORS["muted"]))
    ax[2].set_aspect("equal")
    ax[2].set_title("ideal-flow sphere (teal) vs Stokes (dashed): ψ = 0.5, 1.5")
    ax[2].set_xlabel("x/a")
    ax[0].set_ylabel("y/a")
    save(fig, out, "c13_stokes_sphere")
    # symmetry, residual
    r_t, th_t = np.array([1.5, 2.0, 3.0]), np.array([0.4, 1.1, 2.3])
    s1 = np.hypot(*ch08.stokes_sphere_velocity(r_t, th_t, U, a, "fluid"))
    s2 = np.hypot(*ch08.stokes_sphere_velocity(r_t, np.pi - th_t, U, a, "fluid"))
    print(f"fluid-frame fore–aft symmetry max ||u|(θ) − |u|(π − θ)| = {np.max(np.abs(s1 - s2)):.1e}")
    rng = np.random.default_rng(0)
    P = rng.uniform(-4, 4, (3, 50))
    P = P[:, np.linalg.norm(P, axis=0) > 1.3]
    uf = lambda Q: np.stack(ch08.stokes_sphere_velocity_xyz(Q[0], Q[1], Q[2], U, a))  # noqa: E731

    def pf(Q):
        rr = np.linalg.norm(Q, axis=0)
        return ch08.stokes_sphere_pressure(rr, np.arccos(Q[0] / rr), U, a, mu)

    res = ch08.stokes_residual(uf, pf, P, mu, 1e-3)
    print(f"Stokes residual ∇p − μ∇²u on {P.shape[1]} random points: max {np.max(np.abs(res)):.2e} "
          f"(scale μU/a² = {mu * U / a ** 2})")
    # surface stresses and drag
    th = np.linspace(0, np.pi, 181)
    srr, srt, tx = ch08.stokes_sphere_surface_stresses(th, U, a, mu)
    p_s = ch08.stokes_sphere_pressure(a, th, U, a, mu)
    fig2, a2 = plt.subplots(1, 2, figsize=(12, 4.3))
    a2[0].plot(np.degrees(th), p_s / (mu * U / a), color=COLORS["accent"], lw=2, label="(p − p∞)/(μU/a), (8.50)")
    a2[0].plot(np.degrees(th), srt / (mu * U / a), color=COLORS["orange"], lw=2, label="σ_rθ/(μU/a)")
    a2[0].plot(np.degrees(th), tx / (mu * U / a), color=COLORS["teal"], lw=2, label="x-traction t_x/(μU/a) (uniform)")
    a2[0].set_xlabel("θ from the downstream axis [deg]")
    a2[0].legend(fontsize=8)
    a2[0].set_title("surface stresses (Fig. 8.17)")
    print(f"pressure at the front θ = π: {float(p_s[-1]):+.4f} μU/a, at the rear θ = 0: {float(p_s[0]):+.4f} μU/a")
    D = 6 * np.pi * mu * a * U
    ns = [2, 4, 8, 16]
    errs = []
    for n in ns:
        Dp = ch08.sphere_drag_quadrature(lambda t: -np.asarray(ch08.stokes_sphere_pressure(a, t, U, a, mu)) * np.cos(t),
                                         a, n)
        Df = ch08.sphere_drag_quadrature(lambda t: -np.asarray(ch08.stokes_sphere_surface_stresses(t, U, a, mu)[1])
                                         * np.sin(t), a, n)
        errs.append(abs(Dp + Df - D) / D)
        print(f"Gauss–Legendre n = {n}: pressure {Dp / D:.12f} D, friction {Df / D:.12f} D, total {(Dp + Df) / D:.12f} D")
    print("stokes_drag(parts=True): " + ", ".join(f"{k} = {v:.6f}" for k, v in ch08.stokes_drag(mu, a, U, True).items())
          + f" (6π = {6 * np.pi:.6f})")
    rr = np.logspace(np.log10(5), np.log10(500), 20)
    for Re_a, c in ((0.01, COLORS["accent"]), (0.1, COLORS["orange"])):
        nu = U * a / Re_a
        ratio = ch08.inertia_viscous_ratio(rr, np.pi / 2, U, a, nu)
        sl = np.polyfit(np.log(rr[-8:]), np.log(ratio[-8:]), 1)[0]
        a2[1].loglog(rr, ratio, color=c, lw=2, label=f"Ua/ν = {Re_a}: slope {sl:.3f}")
        a2[1].axvline(1 / Re_a, color=c, ls=":")
        print(f"Ua/ν = {Re_a}: inertia/viscous at r = 500a: {ratio[-1]:.4f}, log–log slope {sl:.4f}")
    a2[1].axhline(1, color=COLORS["muted"], lw=0.8)
    a2[1].set_xlabel("r/a (θ = 90°)")
    a2[1].set_ylabel("|u·∇u| / |ν∇²u|")
    a2[1].legend(fontsize=8)
    a2[1].set_title("Stokes breaks down at r/a ~ 1/Re")
    save(fig2, out, "c14_surface_stresses")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
