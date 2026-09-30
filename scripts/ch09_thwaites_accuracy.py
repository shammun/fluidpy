"""§9.6 (N65, C07): accuracy of Thwaites' method — theta and the wall shear against the exact Falkner-Skan solutions (n > 0 favourable, n < 0 adverse)
and against our marching solver for the linearly retarded flow; the book quotes +-3 % (favourable) and +-10 % (adverse).

Run: ``.venv/Scripts/python.exe scripts/ch09_thwaites_accuracy.py --no-show``   Figure -> outputs/ch09/n65_thwaites_accuracy.png.
"""
from __future__ import annotations

import numpy as np

from ch09_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch09_boundary_layers as ch09


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    ns = np.array([-0.08, -0.06, -0.04, -0.02, 0.0, 0.1, 0.3, 0.6, 1.0, 2.0])
    nu, a = 1e-4, 1.0
    eth, ecf = [], []
    for n in ns:
        of = ch09.outer_flow("wedge", n=float(n), a=a)
        x = np.linspace(1e-3, 1.0, 400)
        # for a pure power law Thwaites' own theta^2 = 0.45 nu x^(1-n)/(a (5n+1)) starts the integration exactly at x0
        th = ch09.thwaites(x, of, nu, theta0=float(np.sqrt(0.45 * nu * 1e-3 ** (1 - n) / (a * (5 * n + 1)))))
        s = ch09.falkner_skan_state(float(n))
        Ue = a * x[-1] ** n
        delta = np.sqrt(nu * x[-1] / Ue)
        eth.append(th["theta"][-1] / (s["I_theta"] * delta) - 1.0)
        ecf.append(th["cf"][-1] / (2 * s["fpp0"] / np.sqrt(Ue * x[-1] / nu)) - 1.0)
    for n, e1, e2 in zip(ns, eth, ecf):
        print(f"  n = {n:+.2f}: theta error {100 * e1:+.2f} %, C_f error {100 * e2:+.2f} %")
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.plot(ns, 100 * np.array(eth), "o-", color=COLORS["accent"], label="theta")
    ax.plot(ns, 100 * np.array(ecf), "s-", color=COLORS["orange"], label="C_f (or wall shear)")
    ax.axhspan(-3, 3, color=COLORS["teal"], alpha=0.12, label="+-3 %")
    ax.axhline(0, color=COLORS["muted"], lw=0.8)
    ax.set_xlabel("Falkner-Skan exponent n  (n < 0 adverse)")
    ax.set_ylabel("Thwaites error vs exact  [%]")
    ax.legend(fontsize=8)
    ax.set_title("Thwaites vs the exact wedge-flow solutions")
    save(fig, out, "n65_thwaites_accuracy")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
