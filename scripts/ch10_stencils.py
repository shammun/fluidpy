"""§10.2 (C01, C03): finite-difference stencils and their order of accuracy (10.6)–(10.7) — error vs h on log–log axes with
the round-off floor — and the FTCS truncation error (10.17) on the advected Gaussian, plus the convergence rates (10.15).

Run: ``.venv/Scripts/python.exe scripts/ch10_stencils.py --no-show``   Figure -> outputs/ch10/c01_stencil_orders.png.
"""
from __future__ import annotations

import numpy as np

from ch10_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch10_computational_fluid_dynamics as ch10


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    kinds = [("forward", 1, COLORS["orange"]), ("backward", 1, COLORS["amber"]), ("central", 1, COLORS["accent"]),
             ("central2", 2, COLORS["blue"]), ("onesided2", 1, COLORS["teal"])]
    for k, m, _ in kinds:
        st = ch10.stencil_taylor_coefficients(k, m)
        print(f"{k:10s} weights {st['weights']}  order {st['order']}  leading coefficient {st['leading']:.5f}")
    h = np.logspace(-8, -0.5, 60)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
    for k, m, c in kinds:
        e = np.array([abs(ch10.stencil_error("sin", 1.0, hh, k, m)) for hh in h])
        ax[0].loglog(h, e + 1e-18, color=c, label=k)
        hh = np.array([0.1, 0.05, 0.025, 0.0125])
        ee = [abs(ch10.stencil_error("sin", 1.0, x, k, m)) for x in hh]
        print(f"  {k:10s} observed order on h = 0.1 … 0.0125: {ch10.observed_order(hh, ee):.3f}")
    ax[0].loglog(h, h, "k:", lw=1, label="slope 1")
    ax[0].loglog(h, h ** 2, "k--", lw=1, label="slope 2")
    ax[0].set_ylim(1e-14, 1)
    ax[0].set_xlabel("h (spacing, dimensionless)")
    ax[0].set_ylabel("|stencil − exact|  for f = sin x at x = 1")
    ax[0].set_title("order of accuracy, and the round-off floor at small h")
    ax[0].legend(fontsize=8)
    x = np.linspace(0, 1, 401)
    tt = ch10.truncation_terms(0.5, 0.01, 0.01, 0.0025, x, 0.1)
    ax[1].plot(x, tt["time"], color=COLORS["orange"], label=r"$\frac{\Delta t}{2}T_{tt}$")
    ax[1].plot(x, tt["conv"], color=COLORS["accent"], label=r"$u\frac{\Delta x^2}{6}T_{xxx}$")
    ax[1].plot(x, tt["diff"], color=COLORS["blue"], label=r"$-D\frac{\Delta x^2}{12}T_{xxxx}$")
    ax[1].plot(x, tt["measured"], "k--", lw=1, label="one-step residual (measured)")
    ax[1].set_xlabel("x  [m]")
    ax[1].set_ylabel("truncation error  [1/s]")
    ax[1].set_title("FTCS truncation error (10.17), u = 0.5 m/s, D = 0.01 m²/s")
    ax[1].legend(fontsize=8)
    t1 = ch10.truncation_terms(0.5, 0.01, 0.01, 0.0025, 0.35, 0.1)
    print("FTCS truncation terms at x = 0.35 m, t = 0.1 s:", {k: f"{v:.4e}" for k, v in t1.items()})
    for s, r in (("ftcs", "diffusive"), ("btcs", "diffusive"), ("cn", "diffusive"), ("upwind", "advective")):
        st = ch10.convergence_study(s, dt_rule=r)
        print(f"convergence {s:6s} ({r}): order {st['order']:.3f}, pairwise {np.round(st['pairwise'], 3).tolist()}")
    for s in ("forward", "backward", "leapfrog", "trapezoidal"):
        print(f"time difference {s:11s} on y' = -y: order {ch10.ode_scheme_order(s):.3f}")
    save(fig, out, "c01_stencil_orders")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
