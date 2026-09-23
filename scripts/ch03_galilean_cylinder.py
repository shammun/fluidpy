"""§3.1/§3.3, Fig. 3.2 (our drawing) and Eq. (3.9): the ideal flow past a cylinder seen from the body (steady) and from
the far fluid (unsteady), and the acceleration of the fluid particle at one point for several observers — the local
and advective parts trade places, their sum does not change (C05, N04, N20–N23). Also the C05 wave-frame number and
the linear/quadratic scaling of the two terms (N22).

Run: ``.venv/Scripts/python.exe scripts/ch03_galilean_cylinder.py --no-show``
Figure → outputs/ch03/fig3_2_galilean_cylinder.png.
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


def galilean_cylinder_figure(U: float = 1.0, a: float = 1.0, P=(0.0, 1.5)):
    """Streamlines in the body frame and in the fluid frame (t = 0) with the point P, plus term bars (x and y parts)
    of the acceleration at P for U_frame = 0, U/2, U. Returns (fig, list of term dicts)."""
    import matplotlib.pyplot as plt

    xs, ys = np.linspace(-3, 3, 241), np.linspace(-2.2, 2.2, 177)
    X, Y = np.meshgrid(xs, ys, indexing="xy")
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.2), gridspec_kw={"width_ratios": [1.2, 1.2, 1]})
    for ax, frame, title in ((axs[0], "body", "body frame: steady"), (axs[1], "fluid", "fluid frame: unsteady (t = 0)")):
        u, v = ch03.cylinder_flow(X, Y, U, a, frame)
        ax.streamplot(xs, ys, np.nan_to_num(u), np.nan_to_num(v), density=1.1, color=COLORS["teal"], linewidth=0.8)
        ax.add_patch(plt.Circle((0, 0), a, color=COLORS["muted"], alpha=0.6))
        ax.plot(*P, "o", color=COLORS["orange"])
        ax.set_aspect("equal")
        ax.set_xlim(xs[0], xs[-1])
        ax.set_ylim(ys[0], ys[-1])
        ax.set_title(title)
        ax.set_xlabel("$x$ [m]")
    axs[0].set_ylabel("$y$ [m]")
    rows = [ch03.frame_acceleration_terms(P[0], P[1], U, a, uf) for uf in (0.0, 0.5 * U, U)]
    bx = axs[2]
    w = 0.25
    for k, (r, uf) in enumerate(zip(rows, (0.0, 0.5 * U, U))):
        x0 = k
        for j, comp in enumerate(("x", "y")):
            xx = x0 + (j - 0.5) * w
            bx.bar(xx, r[f"local_{comp}"], w * 0.9, color=COLORS["blue"], label="local" if k == j == 0 else None)
            bx.bar(xx, r[f"advective_{comp}"], w * 0.9, bottom=r[f"local_{comp}"], color=COLORS["amber"],
                   label="advective" if k == j == 0 else None)
            bx.plot(xx, r[f"total_{comp}"], "D", color=COLORS["accent"], label="total" if k == j == 0 else None)
    bx.set_xticks([0, 1, 2], ["$U_f$ = 0\n(fluid)", "$U_f$ = U/2", "$U_f$ = U\n(body)"])
    bx.axhline(0, color=COLORS["grid"])
    bx.set_ylabel("acceleration at P [m/s²] (left bar x, right bar y)")
    bx.set_title("Eq. (3.9): same total, different split")
    bx.legend(fontsize=8)
    return fig, rows


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
    fig, rows = galilean_cylinder_figure()
    fig.savefig(out / "fig3_2_galilean_cylinder.png", bbox_inches="tight")
    for r, name in zip(rows, ("fluid", "halfway", "body")):
        print(f"P = (0, 1.5) m, {name:7s}: local ({r['local_x']:+.6f}, {r['local_y']:+.6f})  advective "
              f"({r['advective_x']:+.6f}, {r['advective_y']:+.6f})  total ({r['total_x']:+.6f}, {r['total_y']:+.6f}) m/s²")
    # C05 wave number: u′ = 0.5 sin x′ (steady in the wave frame), lab frame moving at U = 2 m/s relative to it
    up = lambda x, t: np.array([0.5 * np.sin(x[0])])  # noqa: E731
    lab = ch03.galilean_transform(up, [-2.0])  # the lab frame moves at −2 m/s relative to the wave frame: u = 2 + u′
    x_lab = np.array([np.pi / 4])
    aw = ch03.acceleration(up, x_lab, 0.0)
    al = ch03.acceleration(lab, x_lab, 0.0)
    print(f"C05 wave: wave frame local {aw.local[0]:+.4f}, advective {aw.advective[0]:+.4f}; lab frame local "
          f"{al.local[0]:+.4f}, advective {al.advective[0]:+.4f}; totals {aw.a[0]:.4f} = {al.a[0]:.4f} m/s²")
    # N22: scale u by λ = 2 in the fluid frame: local ×2, advective ×4
    f1 = ch03.cylinder_velocity_field(1.0, 1.0, "fluid")
    f2 = lambda x, t: 2.0 * f1(x, t)  # noqa: E731
    P = np.array([0.4, 1.3])
    a1, a2 = ch03.acceleration(f1, P, 0.0), ch03.acceleration(f2, P, 0.0)
    print(f"N22: λ = 2 → local ratio {a2.local[1] / a1.local[1]:.4f}, advective ratio {a2.advective[1] / a1.advective[1]:.4f}")
    # streamline check: ψ constant along the computed body-frame streamline
    sl = ch03.streamline(ch03.cylinder_velocity_field(1.0, 1.0), [-3.0, 0.4], 0.0, s_max=6.0, both=False)
    psi = ch03.cylinder_streamfunction(sl[0], sl[1], 1.0, 1.0)
    print(f"body-frame streamline from (−3, 0.4): ψ drift {np.nanmax(np.abs(psi - psi[0])):.1e} m²/s over {sl.shape[1]} points")
    print(f"figure → {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
