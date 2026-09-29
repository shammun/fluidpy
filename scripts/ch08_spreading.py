"""§8.3–8.4 (C08, Examples 8.3 and 8.7, Fig. 8.11): a viscous bead spreading under gravity — the thin-film equation
∂h/∂t = (ρg/3μ)∂(h³∂h/∂x)/∂x solved by our conservative implicit finite-volume scheme, volume conservation, the front
locking onto x_N ∝ t^{1/5} and the profile onto Huppert's similarity solution; the rescaled collapse error.

Run: ``.venv/Scripts/python.exe scripts/ch08_spreading.py --no-show [--fast]``
Figure → outputs/ch08/c08_spreading.png (and the run cached in outputs/ch08/thin_film_run.npz).
"""
from __future__ import annotations

import numpy as np

from ch08_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch08_laminar_flow as ch08


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    rho, g, mu = 1000.0, 9.80665, 1.0  # glycerol-like, 1 Pa s
    N = 201 if args.fast else 401
    x = np.linspace(-0.5, 0.5, N)
    dx = x[1] - x[0]
    h0 = np.where(np.abs(x) < 0.05, 0.01, 0.0)  # 10 cm wide, 1 cm tall block
    area = float(np.sum(h0) * dx / 2)
    ts = np.logspace(0, 3, 13)
    res = ch08.thin_film_spread(h0, x, ts, rho, g, mu, cache=out / "thin_film_run.npz")
    xN = np.array([float(ch08.viscous_current_similarity(0.0, t, area, rho, g, mu, return_front=True)[1]) for t in ts])
    dv = np.max(np.abs(res["volume"] / res["volume"][0] - 1))
    slope = np.polyfit(np.log(ts[-5:]), np.log(res["x_front"][-5:]), 1)[0]
    print(f"cells {N}, steps {res['n_steps']}, Newton max {res['newton_max']}, rejected {res['rejected']}; "
          f"max relative volume change {dv:.2e}")
    print(f"front x_N(t) fitted slope (last 5 times) = {slope:.4f} (similarity 1/5 = 0.2); x_front/x_N(Huppert) at "
          f"t = {ts[-1]:.0f} s: {res['x_front'][-1] / xN[-1]:.4f}")
    print(f"η_N = {ch08.viscous_current_eta_N():.5f}; state at t = 1000 s: "
          + ", ".join(f"{k} = {v:.4g}" for k, v in ch08.thin_film_state(1000.0, area, rho, g, mu).items()))
    print(f"collapse error of Huppert's solution at (n, m) = (1/5, 1/5): "
          f"{ch08.similarity_collapse_error('spreading', 0.2, 0.2):.1e}; at (1/4, 1/4): "
          f"{ch08.similarity_collapse_error('spreading', 0.25, 0.25):.3f}")
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
    for k in (0, 4, 8, 12):
        ax[0].plot(x * 100, res["h"][k] * 1000, lw=2, label=f"t = {ts[k]:.0f} s")
    ax[0].plot(x * 100, h0 * 1000, color=COLORS["muted"], ls=":", label="t = 0")
    hs = ch08.viscous_current_similarity(x, ts[-1], area, rho, g, mu)
    ax[0].plot(x * 100, hs * 1000, "k--", lw=1, label="Huppert, t = 1000 s")
    ax[0].set_xlabel("x [cm]")
    ax[0].set_ylabel("h [mm]")
    ax[0].legend(fontsize=7)
    ax[0].set_title("spreading bead (numerical)")
    ax[1].loglog(ts, res["x_front"] * 100, "o", color=COLORS["accent"], label="numerical front")
    ax[1].loglog(ts, xN * 100, color=COLORS["teal"], label="similarity x_N ∝ t^{1/5}")
    ax[1].set_xlabel("t [s]")
    ax[1].set_ylabel("x_N [cm]")
    ax[1].legend(fontsize=8)
    ax[1].set_title(f"front: late slope {slope:.3f}")
    for k in (6, 9, 12):
        xi = x / ts[k] ** 0.2
        ax[2].plot(xi, res["h"][k] * ts[k] ** 0.2 * 1000, lw=2, label=f"t = {ts[k]:.0f} s")
    ax[2].set_xlim(0, 0.2)
    ax[2].set_xlabel("x / t^{1/5}")
    ax[2].set_ylabel("h t^{1/5} [mm s^{1/5}]")
    ax[2].legend(fontsize=8)
    ax[2].set_title("rescaled profiles collapse")
    save(fig, out, "c08_spreading")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
