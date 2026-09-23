"""§5.5 velocity from vorticity (C07, C08, N24–N29, N52): (5.14) with the corrected +1/(4π) and the printed −1/(4π)
against (5.16) on a finite Gaussian tube; Gauss in curl form (5.15); the periodic Poisson solver; the plane kernel by a
hand-written double loop vs ``biot_savart_2d``; a midpoint sum of (5.17) along a segment vs the closed form; the
segment tending to the infinite line Γ/2πd; the polygon ring converging to ΓR²/2(R² + z²)^{3/2} at order 2; the
elliptic-integral ring against the polygon.

Run: ``.venv/Scripts/python.exe scripts/ch05_biot_savart.py --no-show``
Figures → outputs/ch05/c07_tube_sign.png, c08_segment_ring.png.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from fluidpy import ch05_vorticity_dynamics as ch05  # noqa: E402
from fluidpy.core.style import COLORS  # noqa: E402
from tools.convergence import observed_order  # noqa: E402


def tube_figure(fast: bool):
    """u_θ(r) at the mid-plane of the finite Gaussian tube: (5.16), (5.14) + and − signs, exact references."""
    import matplotlib.pyplot as plt

    tf = ch05.gaussian_tube_fields()
    b = tf["bounds"]
    nodes, w = ch05.cylinder_quadrature(b["radius"], b["z0"], b["z1"], 24 if fast else 32, 32, 48 if fast else 64)
    r = np.linspace(0.3, 1.5, 9 if fast else 13)
    X = np.stack([r, 0 * r, 0 * r])
    u16 = ch05.biot_savart_volume(tf["omega"], X, nodes, w)[1]
    u14 = ch05.velocity_from_curl_omega(tf["curl_omega"], X, nodes, w)[1]
    u14m = ch05.velocity_from_curl_omega(tf["curl_omega"], X, nodes, w, sign=-1.0)[1]
    ref = np.array([tf["u_theta_reference"](rr) for rr in r])
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.plot(r, [tf["u_theta_infinite"](rr) for rr in r], color=COLORS["muted"], ls=":", label="infinite line Γ/2πr")
    ax.plot(r, ref, color=COLORS["ink"], lw=1, label="finite tube (exact)")
    ax.plot(r, u16, "o", color=COLORS["accent"], label="(5.16) Biot–Savart")
    ax.plot(r, u14, "x", color=COLORS["teal"], label="(5.14) with +1/(4π)")
    ax.plot(r, u14m, "s", color=COLORS["rose"], ms=4, label="(5.14) as printed, −1/(4π)")
    ax.axhline(0, color=COLORS["grid"])
    ax.set_xlabel("distance from the tube axis r [m]")
    ax.set_ylabel("u_θ at the mid-plane [m/s]")
    ax.legend(fontsize=8)
    ax.set_title("Γ = 1 m²/s, σ = 0.1 m, length 4 m: the printed sign reverses the swirl")
    return fig, np.abs(u16 - ref).max(), np.abs(u14 - u16).max()


def segment_ring_figure():
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    Ls = np.geomspace(0.5, 200, 30)
    d = 1.0
    sp_ = [ch05.segment_induced_velocity([d, 0, 0], [0, 0, -L], [0, 0, L], 1.0)[1] for L in Ls]
    ax[0].semilogx(Ls / d, sp_, color=COLORS["accent"], label="(Γ/4πd)(cos θ_a − cos θ_b)")
    ax[0].axhline(1 / (2 * np.pi * d), ls="--", color=COLORS["muted"], label="Γ/2πd, (5.2)")
    ax[0].set_xlabel("segment half-length / distance")
    ax[0].set_ylabel("induced speed [m/s]")
    ax[0].legend(fontsize=8)
    Ms = np.array([8, 16, 32, 64, 128, 256])
    err = [abs(ch05.filament_velocity([0, 0, 0.3], ch05.filament_preset("ring", int(M)), 1.0, True)[2]
               - ch05.ring_axis_velocity(0.3, 0.5, 1.0)) for M in Ms]
    ax[1].loglog(Ms, err, "o-", color=COLORS["teal"], label=f"order {observed_order(1.0 / Ms, err):.2f} in 1/M")
    ax[1].set_xlabel("number of segments M")
    ax[1].set_ylabel("|polygon − ΓR²/2(R² + z²)^{3/2}| [m/s]")
    ax[1].legend(fontsize=8)
    return fig, observed_order(1.0 / Ms, err)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch05"))
    ap.add_argument("--no-show", action="store_true")
    ap.add_argument("--fast", action="store_true")
    args = ap.parse_args()
    import matplotlib

    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from fluidpy.core.style import use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    tf = ch05.gaussian_tube_fields()
    print(f"finite Gaussian tube (Γ = 1, σ = 0.1, L = 4) at r = 0.5: exact {tf['u_theta_reference'](0.5):.6f} m/s; "
          f"infinite line {tf['u_theta_infinite'](0.5):.6f} m/s")
    fig, e16, e14 = tube_figure(args.fast)
    fig.savefig(out / "c07_tube_sign.png", bbox_inches="tight")
    print(f"(5.16) vs exact: max error {e16:.1e} m/s; (5.14, +1/4π) vs (5.16): {e14:.1e} m/s; printed sign → −u_θ")
    F = lambda X, Y, Z: np.stack([X * Y ** 2 + np.sin(Z), Y * Z ** 3 + X, X ** 2 * Y * Z])  # noqa: E731
    cc = ch05.curl_theorem_box(F, ((0, 1), (0, 2), (-1, 1)), n=12)
    print(f"(5.15) Gauss in curl form: ∫∇×F dV = {np.round(cc.volume, 10)}, ∮n×F dA = {np.round(cc.surface, 10)}")
    n, Lb = 64, 2 * np.pi
    x = np.arange(n) * Lb / n
    X, Y = np.meshgrid(x, x, indexing="xy")
    uv = ch05.velocity_from_vorticity_fft(-2.0 * np.cos(X) * np.cos(Y), Lb)
    print(f"Poisson (FFT) recovers Taylor–Green: max |u − cos x sin y| = {np.abs(uv[0] - np.cos(X) * np.sin(Y)).max():.1e}")
    # from scratch: the plane kernel by a double loop vs biot_savart_2d on a Gaussian vortex
    g = np.linspace(-1, 1, 41)
    GX, GY = np.meshgrid(g, g, indexing="xy")
    W = 1.0 / (np.pi * 0.2 ** 2) * np.exp(-(GX ** 2 + GY ** 2) / 0.2 ** 2)
    dA = (g[1] - g[0]) ** 2
    px, py = 0.71, 0.23  # off the grid nodes (the kernel is singular on a node)
    u_loop = v_loop = 0.0
    for j in range(g.size):
        for i in range(g.size):
            rx, ry = px - GX[j, i], py - GY[j, i]
            r2 = rx * rx + ry * ry
            u_loop += -W[j, i] * dA * ry / (2 * np.pi * r2)  # e_z × r / 2πr²
            v_loop += W[j, i] * dA * rx / (2 * np.pi * r2)
    ub = ch05.biot_savart_2d(W, GX, GY, [px, py])
    r = np.hypot(px, py)
    print(f"plane kernel by hand ({u_loop:.8f}, {v_loop:.8f}) vs biot_savart_2d ({ub[0]:.8f}, {ub[1]:.8f}); exact "
          f"u_θ = {ch05.gaussian_vortex(r, 1.0, 0.2)[0]:.8f}, |u| = {np.hypot(*ub):.8f}")
    # from scratch: midpoint sum of (5.17) along a straight segment
    Gam, d, a, b = 1.0, 1.0, -1.0, 2.0
    for M in (10, 100, 1000):
        zp = a + (b - a) * (np.arange(M) + 0.5) / M
        dl = (b - a) / M
        du = np.sum(Gam * dl / (4 * np.pi) * d / (d ** 2 + zp ** 2) ** 1.5)  # |e_ω × (x − x′)|/|x − x′|³
        print(f"   (5.17) midpoint sum M = {M:4d}: {du:.10f} m/s")
    print(f"   closed form segment: {ch05.segment_induced_velocity([d, 0, 0], [0, 0, a], [0, 0, b], Gam)[1]:.10f}; "
          f"scalar law {ch05.segment_speed(d, np.arctan2(d, -a) if a < 0 else np.arctan2(d, -a), np.arctan2(d, -b), Gam):.10f}")
    print(f"infinite line Γ/2πd = {1 / (2 * np.pi):.7f}; semi-infinite {ch05.filament_velocity_preset('semi_infinite', 1, 0, 0, length=1e8, component=1):.7f}"
          f"; square loop side 1 at its centre {np.linalg.norm(ch05.filament_velocity([0, 0, 0], ch05.filament_preset('square', 64), 1.0, True)):.6f} "
          f"(2√2/π = {2 * np.sqrt(2) / np.pi:.6f})")
    for M in (8, 16, 64):
        print(f"   polygon ring R = 0.5, M = {M:3d}: centre speed "
              f"{ch05.filament_velocity([0, 0, 0], ch05.filament_preset('ring', M), 1.0, True)[2]:.6f} (exact 1.0)")
    th = 2 * np.pi * np.arange(512) / 512
    poly = np.stack([np.cos(th), np.sin(th), 0.2 + 0 * th])
    up = ch05.filament_velocity([0.7, 0, 0.9], poly, 1.3)
    ue = ch05.ring_ring_velocity(0.7, 0.9, 1.0, 0.2, 1.3)
    print(f"ring by elliptic integrals ({ue[0]:.7f}, {ue[1]:.7f}) vs 512-gon ({up[0]:.7f}, {up[2]:.7f})")
    print(f"Kelvin ring speed Γ/4πR[ln(8R/a) − ¼], R = 1, a = 0.1: {ch05.ring_self_velocity(1.0, 0.1, 1.0):.6f} m/s")
    fig, q = segment_ring_figure()
    fig.savefig(out / "c08_segment_ring.png", bbox_inches="tight")
    print(f"polygon ring axis velocity converges with order {q:.2f} in 1/M")
    print(f"figures → {out}")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
