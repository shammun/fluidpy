"""§7.2 refraction and §7.5 kinematic wave theory (C10): rays dx/dt = ∂ω/∂k, dk/dt = −∂ω/∂x keep ω constant (7.79)
while k, c and c_g change — swell turning toward a plane beach (Fig. 7.8 remake; the ray against the closed-form Snell
ray k sin α = const), rays wrapping round a circular island with a sloping beach (Fig. 7.9 remake) and curved x–t rays
over a slope (Fig. 7.18 remake).

Run: ``.venv/Scripts/python.exe scripts/ch07_refraction.py --no-show [--fast]``
Figures → outputs/ch07/c10_refraction.png.
"""
from __future__ import annotations

import numpy as np

from ch07_drawings import parse_args, save, setup

from fluidpy import ch07_gravity_waves as ch07
from fluidpy.core.style import COLORS


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    T, slope, x0 = 10.0, 0.02, 800.0
    alpha0 = np.radians(40.0)
    om = 2 * np.pi / T
    k0 = ch07.wavenumber_from_omega(om, slope * x0)
    Kv = (-k0 * np.cos(alpha0), k0 * np.sin(alpha0))  # toward the shore (−x) and alongshore (+y)
    r = ch07.ray_trace(None, (x0, 0.0), Kv, (0.0, 400.0), H_fn=lambda X: slope * X[0], n_out=400)
    xs = r["x"][0]
    sn = ch07.snell_ray_plane_beach(xs[1:], alpha0, x0, T, slope)
    dy = np.max(np.abs(sn["y"] - r["x"][1][1:]))
    print(f"plane beach 1:{1 / slope:.0f}, T = {T} s, α₀ = 40° at H = {slope * x0:.0f} m: ray reached x = {xs[-1]:.2f} m "
          f"(H = {slope * xs[-1]:.3f} m); ω drift {r['omega_drift']:.1e}; |k| sin α range "
          f"{np.ptp(r['k'][1]):.1e} rad/m; max |y_ray − y_Snell| = {dy:.2e} m")
    print(f"  angle at H = 2 m: α = {ch07.refraction_state(T, alpha0, slope * x0, 2.0)['alpha_deg']:.2f}° (from 40°)")
    assert r["omega_drift"] < 1e-9 and dy < 1e-3

    fig, ax = plt.subplots(1, 3, figsize=(16, 5.2))
    a0 = ax[0]
    for yoff in np.linspace(-600, 600, 7):
        a0.plot(r["x"][0], r["x"][1] + yoff, color=COLORS["accent"], lw=1)
    for j in range(0, len(xs), 25):  # crest segments ⟂ k at points along the rays
        kk = r["k"][:, j] / np.linalg.norm(r["k"][:, j])
        for yoff in np.linspace(-600, 600, 7):
            px, py = xs[j], r["x"][1][j] + yoff
            a0.plot([px - 40 * kk[1], px + 40 * kk[1]], [py + 40 * kk[0], py - 40 * kk[0]], color=COLORS["orange"], lw=1)
    for Hc in (2, 5, 10, 15):
        a0.axvline(Hc / slope, color=COLORS["muted"], ls=":", lw=0.8)
    a0.axvline(0, color="#6b5b3e", lw=3)
    a0.set_xlabel("distance offshore x [m] (dotted: H = 2, 5, 10, 15 m)")
    a0.set_ylabel("alongshore y [m]")
    a0.set_title("rays (purple) and crests (orange) turn toward the beach", fontsize=10)
    a0.set_xlim(-20, x0 + 20)
    a0.set_aspect("equal")
    # circular island
    R0, Hmax, s_isl = 1000.0, 60.0, 0.03

    def H_island(X):
        rr = np.hypot(X[0], X[1])
        return Hmax * np.tanh(s_isl * (rr - R0) / Hmax)

    a1 = ax[1]
    Ti = 12.0
    ki = ch07.wavenumber_from_omega(2 * np.pi / Ti, Hmax)
    ny = 9 if args.fast else 17
    drifts = []
    for y0 in np.linspace(-4000, 4000, ny):
        ri = ch07.ray_trace(None, (-9000.0, y0), (ki, 0.0), (0.0, 1500.0), H_fn=H_island, n_out=300)
        drifts.append(ri["omega_drift"])
        a1.plot(ri["x"][0] / 1e3, ri["x"][1] / 1e3, color=COLORS["accent"], lw=0.9)
    th = np.linspace(0, 2 * np.pi, 200)
    for Hc, col in ((0.0, "#6b5b3e"), (20.0, COLORS["muted"]), (50.0, COLORS["muted"])):
        rc = R0 + Hmax * np.arctanh(Hc / Hmax) / s_isl
        a1.plot(rc * np.cos(th) / 1e3, rc * np.sin(th) / 1e3, color=col, lw=2 if Hc == 0 else 0.8,
                ls="-" if Hc == 0 else ":")
    print(f"island (R = 1 km, beach slope 0.03, H∞ = 60 m), T = {Ti} s: {ny} rays, max ω drift {max(drifts):.1e}")
    a1.set_xlim(-9, 6)
    a1.set_ylim(-5, 5)
    a1.set_aspect("equal")
    a1.set_xlabel("x [km]")
    a1.set_ylabel("y [km]")
    a1.set_title("rays bend round an island (Fig. 7.9 remake, our rays)", fontsize=10)
    # x–t rays over a slope (1-D)
    a2 = ax[2]
    for kx0 in (0.02, 0.03, 0.05, 0.08):
        r1 = ch07.ray_trace(None, 0.0, kx0, (0.0, 600.0), H_fn=lambda x: 5.0 + 0.05 * x, n_out=200)
        a2.plot(r1["x"][0], r1["t"], color=COLORS["accent"], lw=1.2, label=f"k₀ = {kx0} rad/m" if kx0 in (0.02, 0.08) else None)
    a2.set_xlabel("x [m] (depth H = 5 + 0.05x m)")
    a2.set_ylabel("t [s]")
    a2.set_title("x–t rays dx/dt = c_g curve as the depth grows (Fig. 7.18 remake)", fontsize=10)
    a2.legend(fontsize=8)
    save(fig, out, "c10_refraction")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
