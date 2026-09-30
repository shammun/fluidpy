"""§10.2 (C05, N18, N113): the Lax equivalence theorem made visible — FTCS on a rod heated from both ends (our numbers:
walls at 1, rod at 0, D = 1, L = 1, 20 cells; exact solution (10.199)) at diffusion numbers β = 0.45, 0.50, 0.51.
Consistent at every β, but convergent only while stable (β ≤ ½, (10.28)).

Run: ``.venv/Scripts/python.exe scripts/ch10_lax_demo.py --no-show``   Figure -> outputs/ch10/a1_lax_demo.png.
"""
from __future__ import annotations

import numpy as np

from ch10_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch10_computational_fluid_dynamics as ch10


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    res = ch10.lax_demo(betas=(0.45, 0.50, 0.51), n_cells=20, nsteps=2000)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    for (b, r), c in zip(res.items(), (COLORS["teal"], COLORS["blue"], COLORS["rose"])):
        ax[0].semilogy(np.arange(1, r["max_error_history"].size + 1), r["max_error_history"], color=c, label=f"β = {b}")
        msg = f"blew up at step {r['step_blown']}" if r["blew_up"] else f"final rms error {r['final_error']:.3e}"
        print(f"β = {b:.2f}: {msg}")
    ax[0].set_xlabel("time step n")
    ax[0].set_ylabel("max |T − T_exact|")
    ax[0].set_title("consistent + stable ⇒ convergent (Lax)")
    ax[0].legend()
    # convergence of the stable scheme and of Crank–Nicolson to the series (10.199)
    errs_f, errs_c, hs = [], [], []
    for n in (10, 20, 40, 80):
        dx = 1.0 / n
        x = np.linspace(0, 1, n + 1)
        T0 = np.zeros(n + 1)
        T0[0] = T0[-1] = 1.0
        t_end = 0.05
        dt = 0.25 * dx ** 2
        k = int(round(t_end / dt))
        rf = ch10.solve_transport_1d(T0, x, 0.0, 1.0, t_end / k, k, "ftcs", g=1.0, T_L=1.0)
        kc = int(round(t_end / dx ** 2))
        rc = ch10.solve_transport_1d(T0, x, 0.0, 1.0, t_end / kc, kc, "cn", g=1.0, T_L=1.0)
        ex = ch10.rod_heating_exact(x, t_end, 1.0, 1.0)
        errs_f.append(ch10.error_norm(rf["T"], ex))
        errs_c.append(ch10.error_norm(rc["T"], ex))
        hs.append(dx)
    print(f"FTCS (β = 0.25) order in Δx: {ch10.observed_order(hs, errs_f):.3f};  CN (β = 1) order: "
          f"{ch10.observed_order(hs, errs_c):.3f}")
    ax[1].loglog(hs, errs_f, "o-", color=COLORS["orange"], label="FTCS, β = 0.25")
    ax[1].loglog(hs, errs_c, "s-", color=COLORS["blue"], label="Crank–Nicolson, β = 1")
    ax[1].loglog(hs, np.array(hs) ** 2 * errs_f[0] / hs[0] ** 2, "k--", lw=1, label="slope 2")
    ax[1].set_xlabel("Δx  [m]")
    ax[1].set_ylabel("rms error at t = 0.05 s")
    ax[1].set_title("stable schemes converge to (10.199)")
    ax[1].legend(fontsize=8)
    save(fig, out, "a1_lax_demo")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
