"""Fig. 3.18 (our 2-D drawing) and the Reynolds transport theorem (3.35) as a budget: a deforming, translating ellipse
over a field F, its boundary velocity b, the swept band coloured by the sign of b·n, and waterfall bars
volume term + surface term = total against the measured rate of ∫F dA; the (3.32) terms for a growing sphere vs Δt
(orders of smallness, slope 2); Example 3.2 three ways (C15, N46–N55).

Run: ``.venv/Scripts/python.exe scripts/ch03_fig3_18_rtt.py --no-show``
Figures → outputs/ch03/fig3_18_rtt.png, rtt_dt_convergence.png.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from fluidpy import ch03_kinematics as ch03  # noqa: E402
from fluidpy.core.style import COLORS  # noqa: E402


def _ellipse_cv(a, b, adot, bdot, c, cdot, t0=0.0):
    c0, cd = np.asarray(c, float), np.asarray(cdot, float)
    return ch03.MovingEllipse2D(a_fn=lambda t: a + adot * (t - t0), b_fn=lambda t: b + bdot * (t - t0),
                                center_fn=lambda t: c0 + cd * (t - t0), rates=lambda t: (adot, bdot, cd[0], cd[1]))


def rtt_blob_figure(a: float = 1.0, b: float = 0.6, adot: float = 0.3, bdot: float = -0.2, cdot=(0.4, 0.0),
                    F_name: str = "warming", dt: float = 0.25, c=(0.0, 0.0), t: float = 0.0):
    """Two panels: (left) F(x, t) heatmap, the ellipse at t (solid) and t + dt (dashed), boundary velocities b and the
    swept band (blue where b·n > 0 — the surface advances and sweeps F in; rose where b·n < 0); (right) waterfall
    bars ∫∂F/∂t dA + ∮F b·n ds = total with the finite-difference d/dt∫F dA as a marker. Returns (fig, terms, fd)."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    F, dFdt = ch03.rtt_field(F_name, dim=2)
    terms = ch03.rtt_ellipse_2d(F, dFdt, a, b, adot, bdot, c, cdot, t)
    cv = _ellipse_cv(a, b, adot, bdot, c, cdot, t)
    fd = ch03.volume_integral_rate_fd(F, cv, t, 1e-4)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11.0, 4.4), gridspec_kw={"width_ratios": [1.5, 1]})
    xs, ys = np.linspace(-2.2, 2.6, 200), np.linspace(-1.6, 1.6, 140)
    X, Y = np.meshgrid(xs, ys, indexing="xy")
    im = ax.pcolormesh(X, Y, F(np.stack([X.ravel(), Y.ravel()]), t).reshape(X.shape), cmap="Greys", alpha=0.35,
                       shading="auto")
    fig.colorbar(im, ax=ax, label=f"F ('{F_name}')")
    phi = np.linspace(0, 2 * np.pi, 97)
    for tt, ls in ((t, "-"), (t + dt, "--")):
        aa, bb, cc = cv.a_fn(tt), cv.b_fn(tt), np.asarray(cv.center_fn(tt))
        ax.plot(cc[0] + aa * np.cos(phi), cc[1] + bb * np.sin(phi), ls, color=COLORS["ink"], lw=1.6)
    Xs, N, _, B = cv.surface_nodes(t, 6)  # 48 boundary points
    bn = np.einsum("ik,ik->k", B, N)
    nb = Xs.shape[1]  # midpoint nodes on a periodic boundary (no duplicated end point)
    for k in range(nb):  # swept band between the boundary at t and at t + dt, closed: last node joins node 0
        k1 = (k + 1) % nb
        seg = np.stack([Xs[:, k], Xs[:, k1], Xs[:, k1] + B[:, k1] * dt, Xs[:, k] + B[:, k] * dt])
        col = COLORS["blue"] if bn[k] > 0 else COLORS["rose"]
        ax.add_patch(Polygon(seg, closed=True, color=col, alpha=0.45, lw=0))
    ax.quiver(Xs[0, ::3], Xs[1, ::3], B[0, ::3], B[1, ::3], color=COLORS["orange"], angles="xy", scale_units="xy",
              scale=1.2, width=0.004, label="b")
    ax.set_aspect("equal")
    ax.set_xlim(xs[0], xs[-1])
    ax.set_ylim(ys[0], ys[-1])
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$y$ [m]")
    ax.set_title("moving control area: blue b·n > 0 (sweeps F in), rose b·n < 0")
    vt, st, tot = terms
    bx.bar([0], [vt], color=COLORS["blue"], label=r"$\int\partial F/\partial t\,dA$")
    bx.bar([1], [st], bottom=[vt], color=COLORS["orange"], label=r"$\oint F\,\mathbf{b}\cdot\mathbf{n}\,ds$")
    bx.bar([2], [tot], color=COLORS["accent"], label="total (3.35)")
    bx.plot([2], [fd], "D", color=COLORS["ink"], ms=8, label=r"measured $d/dt\int F\,dA$")
    bx.set_xticks([0, 1, 2], ["volume", "surface", "total"])
    bx.set_ylabel("rate [F·m²/s]")
    bx.set_title("Reynolds transport budget")
    bx.legend(fontsize=8)
    return fig, terms, fd


def dt_convergence_figure(R: float = 1.0, Rdot: float = 0.1, t: float = 1.0, dts=None):
    """The (3.32) terms for the growing sphere GrowingSphere(R, Rdot) with F = 'warming', vs Δt on log–log axes:
    T2 and T3 fall like Δt, the dropped T4, the sliver error of (3.34) and the Taylor residual like Δt². Returns
    (fig, dict of observed slopes)."""
    import matplotlib.pyplot as plt

    from tools.convergence import observed_order

    dts = np.logspace(-1, -4, 7) if dts is None else np.asarray(dts)
    F, dFdt = ch03.rtt_field("warming")
    rows = [ch03.swept_terms_sphere(R, Rdot, F, dFdt, t, d) for d in dts]
    series = {"T2": [r["T2"] for r in rows], "T3": [r["T3"] for r in rows], "T4 (dropped)": [r["T4"] for r in rows],
              "sliver error (3.34)": [abs(r["sliver_error"]) for r in rows]}
    total = ch03.reynolds_transport(F, dFdt, ch03.GrowingSphere(R, Rdot), t).total
    series["|lhs − RTT total|"] = [abs(r["lhs"] - total) for r in rows]
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    slopes = {}
    for (k, v), col in zip(series.items(), (COLORS["blue"], COLORS["teal"], COLORS["rose"], COLORS["orange"], COLORS["accent"])):
        v = np.abs(np.asarray(v))
        slopes[k] = observed_order(dts, v)
        ax.loglog(dts, v, "o-", color=col, label=f"{k} (slope {slopes[k]:.2f})")
    ax.set_xlabel(r"$\Delta t$ [s]")
    ax.set_ylabel("term [F·m³]  (last: rate error [F·m³/s])")
    ax.set_title("Eq. (3.32): orders of smallness for a growing sphere")
    ax.legend(fontsize=8)
    return fig, slopes


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch03"))
    ap.add_argument("--no-show", action="store_true")
    args = ap.parse_args()
    t0 = time.perf_counter()
    import matplotlib

    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from fluidpy.core.style import use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    fig, terms, fd = rtt_blob_figure()
    fig.savefig(out / "fig3_18_rtt.png", bbox_inches="tight")
    print(f"2-D ellipse RTT: volume {terms.volume_term:.9f} + surface {terms.surface_term:.9f} = {terms.total:.9f}; "
          f"measured {fd:.9f} (diff {abs(fd - terms.total):.1e})")
    fig2, slopes = dt_convergence_figure()
    fig2.savefig(out / "rtt_dt_convergence.png", bbox_inches="tight")
    print("(3.32) slopes vs Δt:", {k: round(v, 3) for k, v in slopes.items()})
    F1, dF1 = (lambda x, t: t + 0 * x[0]), (lambda x, t: 1 + 0 * x[0])
    cvs = ch03.GrowingSphere(1.0 - 0.1, 0.1)  # radius 1 m at t = 1 s, growing at 0.1 m/s
    r = ch03.reynolds_transport(F1, dF1, cvs, 1.0)
    print(f"C15 worked number (F = t, R = 1 m, Ṙ = 0.1 m/s, t = 1 s): volume {r.volume_term:.4f} + surface "
          f"{r.surface_term:.4f} = {r.total:.4f}; measured {ch03.volume_integral_rate_fd(F1, cvs, 1.0):.4f}")
    for cv in (ch03.MovingBox((1, 2, 0.5), (0.1, -0.2, 0.3), (1, 0, 0)), ch03.GrowingCylinder(0.5, 0.1, 2.0, 0.2),
               ch03.GrowingCone(0.5, 0.1, 1.0)):
        chk = ch03.rtt_check(*ch03.rtt_field("warming"), cv, 0.5)
        print(f"  {type(cv).__name__:15s} RTT {chk.rhs_total:.9f} vs FD {chk.lhs_fd:.9f} (rel {chk.rel_error:.1e})")
    ex = ch03.example_3_2(1.0, 0.5, 0.1)
    print("Example 3.2 (h = 1 m, r₀ = 0.5 m, ṙ = 0.1 m/s): " + ", ".join(f"{k} {v:.10f}" for k, v in ex.items()))
    mv = ch03.material_volume_rate(lambda x, t: np.stack([x[0], x[1], x[2]]), ch03.GrowingSphere(0.1, 0.0), 0.0)
    print(f"material volume rate for u = x (∇·u = 3): ∮u·n dA = {mv.surface_flux:.6e}, ∫∇·u dV = {mv.volume_div:.6e}, "
          f"per volume {mv.surface_flux / (4 / 3 * np.pi * 0.1 ** 3):.6f} s⁻¹")
    print(f"figures → {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
