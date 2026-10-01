"""Schematic drawings for chapter 11 (our analogues of the book's sketches; no physics computed here beyond the profiles):
KH two-stream set-up (Fig. 11.2), Bénard layer (Fig. 11.8), Taylor–Couette gap (Fig. 11.16), the control volume of (11.88)
(Fig. 11.25), Howard's semicircle (Fig. 11.20) and the six parallel profiles of Fig. 11.21.

Run: ``.venv/Scripts/python.exe scripts/ch11_drawings.py --no-show``   Figure -> outputs/ch11/drawings.png.
The helpers take an Axes and are importable by the notebook (``sys.path.insert(0, "scripts")``).
"""
from __future__ import annotations

import math

import numpy as np

from ch11_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch11_instability as ch11


def kh_sketch(ax, U1: float = 1.0, U2: float = -0.6, amp: float = 0.08):
    """Two streams (U₁ above, U₂ below) separated by a wavy interface ζ = a cos kx (Fig. 11.2 analogue)."""
    x = np.linspace(0, 4 * math.pi, 300)
    ax.fill_between(x, amp * np.cos(x), 1.0, color=COLORS["grid"], alpha=0.6)
    ax.fill_between(x, -1.0, amp * np.cos(x), color=COLORS["blue"], alpha=0.25)
    ax.plot(x, amp * np.cos(x), color=COLORS["ink"])
    for y0, U, lab in ((0.55, U1, r"$U_1,\ \rho_1$"), (-0.55, U2, r"$U_2,\ \rho_2$")):
        for xx in (2.0, 6.0, 10.0):
            ax.annotate("", (xx + 1.2 * U, y0), (xx, y0), arrowprops=dict(arrowstyle="->", color=COLORS["orange"], lw=2))
        ax.text(0.3, y0 + 0.15, lab, fontsize=11)
    ax.text(8.5, 0.15, r"$z=\zeta(x,t)$", fontsize=10)
    ax.set(xlim=(0, 4 * math.pi), ylim=(-1, 1), xticks=[], yticks=[], title="Kelvin–Helmholtz: two streams (cf. Fig. 11.2)")


def benard_sketch(ax):
    """Layer of depth d heated from below, conduction profile T̄ = T₀ − Γ(z + d/2) (Fig. 11.8 analogue)."""
    import matplotlib.patches as mp

    ax.add_patch(mp.Rectangle((0, -0.55), 4, 0.05, color=COLORS["rose"]))
    ax.add_patch(mp.Rectangle((0, 0.5), 4, 0.05, color=COLORS["blue"]))
    z = np.linspace(-0.5, 0.5, 50)
    ax.plot(1.0 - 0.8 * (z + 0.5), z, color=COLORS["ink"])
    ax.text(1.05, -0.45, r"$T_0$")
    ax.text(0.3, 0.38, r"$T_0-\Delta T$")
    ax.text(1.4, 0.0, r"$\bar T=T_0-\Gamma(z+d/2)$, $\Gamma=\Delta T/d$")
    ax.annotate("", (3.6, 0.5), (3.6, -0.5), arrowprops=dict(arrowstyle="<->"))
    ax.text(3.65, 0.0, "d")
    ax.set(xlim=(0, 4), ylim=(-0.7, 0.7), xticks=[], yticks=[], title="Bénard layer (cf. Fig. 11.8)")


def taylor_sketch(ax, kc: float = 3.12):
    """Meridional gap with counter-rotating Taylor vortices of wavelength ≈ 2d (Fig. 11.16 analogue)."""
    e = ch11.taylor_eigenfunction(kc, 0.0)
    ax.contour(e["x"], e["z"], e["psi"], 12, colors=COLORS["accent"], linewidths=1)
    ax.axvline(0, color=COLORS["ink"], lw=3)
    ax.axvline(1, color=COLORS["ink"], lw=3)
    ax.set(xlabel="x = (R − R₁)/d", ylabel="z/d", title="Taylor vortices in the gap (cf. Fig. 11.16)")
    ax.set_aspect("equal")


def cv_sketch(ax):
    """Control volume between walls, an integer number of wavelengths long (Fig. 11.25 analogue)."""
    import matplotlib.patches as mp

    ax.add_patch(mp.Rectangle((0, 0.95), 6, 0.08, color=COLORS["muted"]))
    ax.add_patch(mp.Rectangle((0, -0.03), 6, 0.08, color=COLORS["muted"]))
    ax.add_patch(mp.Rectangle((1, 0.08), 4, 0.84, fill=False, ls="--", ec=COLORS["accent"], lw=2))
    ax.text(2.2, 0.45, "V (integer number\nof wavelengths)")
    ax.text(5.1, 0.5, "u_i periodic")
    ax.set(xlim=(0, 6), ylim=(-0.1, 1.1), xticks=[], yticks=[], title="Control volume of (11.88) (cf. Fig. 11.25)")


def semicircle_sketch(ax, Umin: float = -1.0, Umax: float = 1.0, c=None):
    """Howard's semicircle (p. 507) with optional eigenvalues c."""
    cr, ci = ch11.howard_semicircle(Umin, Umax)
    ax.fill(cr, ci, color=COLORS["accent"], alpha=0.15)
    ax.plot(cr, ci, color=COLORS["accent"])
    if c is not None:
        c = np.atleast_1d(c)
        ax.plot(c.real, c.imag, "o", color=COLORS["rose"])
    ax.set(xlabel="$c_r$", ylabel="$c_i$", title="Howard's semicircle (p. 507)")
    ax.set_aspect("equal")


def six_profiles(axes):
    """Our stand-ins for the six profiles of Fig. 11.21 with their inflection points and verdicts."""
    verdicts = {d["panel"]: d for d in ch11.fig_11_21_verdicts()}
    for ax, (panel, name) in zip(np.ravel(axes), ch11.FIG_11_21.items()):
        pr = ch11.parallel_profile(name, b=math.pi) if name == "sin" else ch11.parallel_profile(name)
        lo, hi = pr["domain"] if pr["bc"] == "wall" else (0.0, 6.0)
        y = np.linspace(lo, hi, 300)
        ax.plot(pr["U"](y), y, color=COLORS["ink"])
        ax.axvline(0, color=COLORS["muted"], lw=0.8)
        v = verdicts[panel]
        for yI in v["y_I"]:
            ax.plot(pr["U"](yI), yI, "o", color=COLORS["rose"])
        ax.set_title(f"({panel}) Rayleigh {'✓' if v['rayleigh'] else '✗'}, Fjørtoft {'✓' if v['fjortoft'] else '✗'}",
                     fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(15, 9))
    gs = fig.add_gridspec(3, 6)
    kh_sketch(fig.add_subplot(gs[0, 0:2]))
    benard_sketch(fig.add_subplot(gs[0, 2:4]))
    taylor_sketch(fig.add_subplot(gs[0, 4:6]))
    cv_sketch(fig.add_subplot(gs[1, 0:2]))
    semicircle_sketch(fig.add_subplot(gs[1, 2:4]), c=[0.0 + 0.42j, 0.3 + 0.2j])
    six_profiles([fig.add_subplot(gs[2, j]) for j in range(6)])
    save(fig, out, "drawings")
    for d in ch11.fig_11_21_verdicts():
        print(f"Fig. 11.21 ({d['panel']}) {d['label']:38s} Rayleigh {d['rayleigh']!s:5s} Fjørtoft {d['fjortoft']}")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
