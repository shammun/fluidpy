"""Example 3.1 / Fig. 3.7 (our drawing): the streamline, path line and streak line through the origin at t = t′ in the
unsteady uniform flow u = ωξ_o cos ωt, v = ωξ_o sin ωt — closed forms overlaid on the numerical
``streamline``/``pathline``/``streakline`` of ``fluidpy.core.kinematics``; plus the streamline pattern at three
instants (C03) and an optional animation (dye + one particle + turning streamline over a period).

Run: ``.venv/Scripts/python.exe scripts/ch03_fig3_7_flow_lines.py --no-show [--anim]``
Figures → outputs/ch03/fig3_7_flow_lines.png, fig3_7_streamlines_three_instants.png (+ fig3_7_flow_lines.gif).
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

C_STREAM, C_PATH, C_STREAK = COLORS["teal"], COLORS["orange"], COLORS["rose"]


def flow_lines_figure(t_prime: float, xi0: float = 1.0, omega: float = 1.0, ax=None, numeric: bool = True):
    """Draw the three flow lines of Example 3.1 at the drawing instant t′ (closed forms, plus the numerical curves as
    dots when ``numeric``), the dye port at the origin and the common tangent. Returns (fig, ax, info) with ``info``
    the max distance between numeric and closed-form curves [m]."""
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(5.6, 5.6))
    else:
        fig = ax.figure
    ex = ch03.example_3_1(t_prime, xi0, omega, n=300)
    ax.plot(*ex["streamline"], color=C_STREAM, lw=2.2, label=r"streamline $y = x\tan\omega t'$")
    ax.plot(*ex["pathline"], color=C_PATH, lw=2.2, label="path line (particle at O at $t'$)")
    ax.plot(*ex["streakline"], color=C_STREAK, lw=2.2, label="streak line (dye from O)")
    for c, col in ((ex["path_center"], C_PATH), (ex["streak_center"], C_STREAK)):
        ax.plot(*c, marker="+", ms=10, color=col)
    info = {}
    if numeric:
        u = ex["field"]
        T = 2 * np.pi / omega
        sl = ch03.streamline(u, [0.0, 0.0], t_prime, s_max=2 * xi0, n=41)
        pl = ch03.pathline(u, [0.0, 0.0], t_prime, np.linspace(t_prime, t_prime + T, 41))
        sk = ch03.streakline(u, [0.0, 0.0], t_prime, np.linspace(t_prime - T, t_prime, 41))
        ax.plot(*sl, ".", color=C_STREAM, ms=4)
        ax.plot(*pl, ".", color=C_PATH, ms=4)
        ax.plot(*sk, ".", color=C_STREAK, ms=4)
        rad = lambda P, c: float(np.max(np.abs(np.hypot(P[0] - c[0], P[1] - c[1]) - xi0)))  # noqa: E731
        info = {"streamline": float(np.max(np.abs(sl[1] * np.cos(omega * t_prime) - sl[0] * np.sin(omega * t_prime)))),
                "pathline": rad(pl, ex["path_center"]), "streakline": rad(sk, ex["streak_center"])}
    ax.plot(0, 0, "s", color=COLORS["ink"], ms=6)
    ax.annotate("dye port O", (0, 0), (0.15 * xi0, -0.35 * xi0), fontsize=9, color=COLORS["ink"])
    ax.set_aspect("equal")
    ax.set_xlim(-2.3 * xi0, 2.3 * xi0)
    ax.set_ylim(-2.3 * xi0, 2.3 * xi0)
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$y$ [m]")
    ax.set_title(rf"Example 3.1 at $\omega t' = {omega * t_prime:.2f}$ rad (dots: numerical)")
    ax.legend(loc="upper left", fontsize=8)
    return fig, ax, info


def flow_lines_animation(frames: int = 60, xi0: float = 1.0, omega: float = 1.0, n_dye: int = 60):
    """Animation over one period: the instantaneous streamline through O turns, one particle released at O at t = 0
    traces its path line, and dye released continuously from O forms the streak line. Returns a ``FuncAnimation``."""
    import matplotlib.pyplot as plt

    from fluidpy.core.anim import animate

    T = 2 * np.pi / omega
    u = ch03.preset_field("ex31", xi0=xi0, omega=omega)
    times = np.linspace(0.0, T, frames)
    path = ch03.pathline(u, [0.0, 0.0], 0.0, times)
    fig, ax = plt.subplots(figsize=(5.0, 5.0))
    ax.set_aspect("equal")
    ax.set_xlim(-2.3 * xi0, 2.3 * xi0)
    ax.set_ylim(-2.3 * xi0, 2.3 * xi0)
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$y$ [m]")
    (ls,) = ax.plot([], [], color=C_STREAM, lw=2, label="streamline now")
    (lp,) = ax.plot([], [], color=C_PATH, lw=2, label="path line (released at t = 0)")
    (pp,) = ax.plot([], [], "o", color=C_PATH, ms=7)
    (lk,) = ax.plot([], [], ".", color=C_STREAK, ms=5, label="dye (streak line)")
    ax.plot(0, 0, "s", color=COLORS["ink"])
    ttl = ax.set_title("")
    ax.legend(loc="upper left", fontsize=8)

    def update(i):
        t = times[i]
        d = np.array([np.cos(omega * t), np.sin(omega * t)])
        s = np.linspace(-2.2 * xi0, 2.2 * xi0, 2)
        ls.set_data(s * d[0], s * d[1])
        lp.set_data(path[0, : i + 1], path[1, : i + 1])
        pp.set_data([path[0, i]], [path[1, i]])
        if t > 0:
            dye = ch03.streakline(u, [0.0, 0.0], t, np.linspace(0.0, t, max(2, int(n_dye * t / T))))
            lk.set_data(dye[0], dye[1])
        ttl.set_text(rf"$\omega t$ = {omega * t:.2f} rad")
        return ls, lp, pp, lk

    return animate(update, frames=frames, fig=fig, interval=60)


def streamlines_three_instants(xi0: float = 1.0, omega: float = 1.0, t_primes=(0.0, np.pi / 4, np.pi / 2)):
    """C03's figure: the Example 3.1 streamline pattern frozen at three instants (it turns), with the path line and
    streak line of the first instant fixed. Returns the figure."""
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(1, len(t_primes), figsize=(4.0 * len(t_primes), 4.0))
    for ax, tp in zip(axs, t_primes):
        ex = ch03.example_3_1(tp, xi0, omega)
        for y0 in np.linspace(-2, 2, 9) * xi0:  # parallel streamlines of the uniform field at this instant
            d = np.array([np.cos(omega * tp), np.sin(omega * tp)])
            s = np.linspace(-4, 4, 2) * xi0
            p0 = np.array([-d[1], d[0]]) * y0
            ax.plot(p0[0] + s * d[0], p0[1] + s * d[1], color=C_STREAM, lw=1.0, alpha=0.6)
        ax.plot(*ex["streamline"], color=C_STREAM, lw=2.4)
        ex0 = ch03.example_3_1(t_primes[0], xi0, omega)
        ax.plot(*ex0["pathline"], color=C_PATH, lw=1.2, ls="--")
        ax.plot(*ex0["streakline"], color=C_STREAK, lw=1.2, ls="--")
        ax.set_aspect("equal")
        ax.set_xlim(-2.2, 2.2)
        ax.set_ylim(-2.2, 2.2)
        ax.set_title(rf"streamlines at $\omega t' = {omega * tp:.2f}$")
        ax.set_xlabel("$x$ [m]")
    axs[0].set_ylabel("$y$ [m]")
    return fig


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch03"))
    ap.add_argument("--no-show", action="store_true")
    ap.add_argument("--anim", action="store_true", help="also save the GIF animation (slower)")
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
    tp = np.pi / 6
    fig, _, info = flow_lines_figure(tp)
    fig.savefig(out / "fig3_7_flow_lines.png", bbox_inches="tight")
    ex = ch03.example_3_1(tp)
    print(f"Example 3.1 at ωt′ = {tp:.4f}: slope tan ωt′ = {ex['slope']:.6f}")
    print(f"  path-line centre {ex['path_center'].round(6)}, streak-line centre {ex['streak_center'].round(6)}, radius {ex['radius']}")
    print(f"  numeric vs closed form (max distance, m): {', '.join(f'{k} {v:.2e}' for k, v in info.items())}")
    fig2 = streamlines_three_instants()
    fig2.savefig(out / "fig3_7_streamlines_three_instants.png", bbox_inches="tight")
    # steady preset: the three lines coincide
    u = ch03.preset_field("steady_vortex", Omega=1.0)
    x0 = np.array([1.0, 0.0])
    pl = ch03.pathline(u, x0, 0.0, np.linspace(0, 2, 30))
    sk = ch03.streakline(u, x0, 2.0, np.linspace(0, 2, 30))
    print(f"steady vortex: path line and streak line radius drift {np.max(np.abs(np.hypot(*pl) - 1)):.1e}, "
          f"{np.max(np.abs(np.hypot(*sk) - 1)):.1e} (both on the streamline r = 1)")
    if args.anim:
        from fluidpy.core.anim import save_gif

        save_gif(flow_lines_animation(frames=40), out / "fig3_7_flow_lines.gif", fps=12)
    print(f"figures → {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
