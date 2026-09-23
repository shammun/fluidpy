"""Fig. 3.17 (our drawing): Leibniz's theorem (3.30) as three strips — the change inside dt∫∂F/∂t dx, the gain at the
moving upper limit db·F(b, t) and the loss at the moving lower limit da·F(a, t) — plus the three terms against the
measured rate of the integral for the ``leibniz_example`` cases (N46, N47, D21).

Run: ``.venv/Scripts/python.exe scripts/ch03_fig3_17_leibniz.py --no-show``
Figure → outputs/ch03/fig3_17_leibniz.png.
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


def leibniz_strips_figure(case: str = "bump", t: float = 0.8, dt: float = 0.25, ax=None):
    """F(x, t) and F(x, t + dt) over the moving interval with the three shaded contributions of (3.30):
    interior (dt·∂F/∂t, blue dots between the curves), gain at b (orange), loss at a (rose). A finite dt is used so
    the strips are visible; their sum differs from the true change by O(dt²). Returns (fig, ax, terms)."""
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(7.0, 4.2))
    else:
        fig = ax.figure
    d0, d1 = ch03.leibniz_example(t, case), ch03.leibniz_example(t + dt, case)
    F = d0["F"]
    a0, b0, a1, b1 = d0["a"], d0["b"], d1["a"], d1["b"]
    lo, hi = min(a0, a1) - 0.3, max(b0, b1) + 0.3
    x = np.linspace(lo, hi, 400)
    ax.plot(x, F(x, t), color=COLORS["ink"], lw=2, label="$F(x, t)$")
    ax.plot(x, F(x, t + dt), color=COLORS["muted"], lw=1.6, ls="--", label=r"$F(x, t + \Delta t)$")
    xi = np.linspace(a0, b0, 200)
    ax.fill_between(xi, F(xi, t), F(xi, t + dt), color=COLORS["blue"], alpha=0.25, label=r"change inside $\approx\Delta t\int\partial F/\partial t\,dx$")
    xb = np.linspace(min(b0, b1), max(b0, b1), 50)
    ax.fill_between(xb, 0, F(xb, t), color=COLORS["orange"], alpha=0.45, label=r"gain at $b$: $\Delta b\,F(b, t)$")
    xa = np.linspace(min(a0, a1), max(a0, a1), 50)
    ax.fill_between(xa, 0, F(xa, t), color=COLORS["rose"], alpha=0.45, hatch="//", label=r"loss at $a$: $\Delta a\,F(a, t)$")
    for xv, lbl in ((a0, "$a(t)$"), (b0, "$b(t)$")):
        ax.axvline(xv, color=COLORS["grid"], lw=1)
        ax.annotate(lbl, (xv, 0), xytext=(2, 4), textcoords="offset points", fontsize=9)
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$F$")
    ax.set_title(rf"Leibniz (3.30), case '{case}', $t$ = {t:g} s, $\Delta t$ = {dt:g} s")
    ax.legend(fontsize=8, loc="upper right")
    return fig, ax, d0


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
    fig, _, _ = leibniz_strips_figure()
    fig.savefig(out / "fig3_17_leibniz.png", bbox_inches="tight")
    for case in ("x2t", "wave", "bump"):
        d = ch03.leibniz_example(1.3, case)
        c = ch03._leibniz_case(case)  # noqa: SLF001 — the case's a(t), b(t) callables for the independent FD check
        fd = ch03.leibniz_check(c["F"], lambda s: float(c["a"](s)), lambda s: float(c["b"](s)), 1.3)
        print(f"{case:5s} t = 1.3: interior {d['interior']:+.6f} + upper {d['upper']:+.6f} − lower {d['lower']:+.6f} "
              f"= {d['total']:+.9f} | sympy {d['exact']:+.9f} | central difference {fd:+.9f}")
    print(f"figure → {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
