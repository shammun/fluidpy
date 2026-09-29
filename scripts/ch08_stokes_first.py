"""§8.4 (C09): Stokes' first problem — profiles (8.30) at several times (Fig. 8.12), their collapse in η (Fig. 8.13),
the 99 % thickness (8.31), the Crank–Nicolson solution of (8.20) and its observed order, the vorticity content ∫ω dy,
and the stopped plate of Exercise 8.30 (no similarity).

Run: ``.venv/Scripts/python.exe scripts/ch08_stokes_first.py --no-show [--fast]``
Figure → outputs/ch08/c09_stokes_first.png.
"""
from __future__ import annotations

import numpy as np

from ch08_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch08_laminar_flow as ch08
from fluidpy.core import diffusion as DIF


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    U, nu = 1.0, 1e-6
    y = np.linspace(0, 0.02, 401)
    times = (10.0, 40.0, 160.0)
    fig, ax = plt.subplots(1, 4, figsize=(17, 4.4))
    curves = []
    for t in times:
        u = ch08.stokes_first_problem(y, t, U, nu)
        d = float(ch08.diffusion_thickness(t, nu))
        ax[0].plot(u, y * 1000, lw=2, label=f"t = {t:.0f} s, δ₉₉ = {d * 1000:.2f} mm")
        eta = ch08.similarity_variable(y, t, nu, half=True)
        ax[1].plot(u, eta, lw=2)
        curves.append(ch08.stokes_first_problem(np.linspace(0, 4, 50) * 2 * np.sqrt(nu * t), t, U, nu))
    spread = np.max(np.ptp(np.array(curves), axis=0))
    print(f"collapse in η: max spread over t = {times} → {spread:.1e}")
    print(f"99 % thickness coefficient 2 erfc⁻¹(0.01) = {float(ch08.diffusion_thickness(1.0, 1.0)):.5f}; on the "
          f"y/(2√(νt)) axis {float(ch08.diffusion_thickness(1.0, 1.0)) / 2:.5f}")
    for t in times:
        v, e = ch08.vorticity_content(t, U, nu, return_error=True)
        print(f"t = {t:.0f} s: ∫₀^∞ ω dy = {v:.12f} m/s (quad error {e:.1e}); τ_w = "
              f"{ch08.stokes_first_state(t, U, nu)['tau_w']:.4e} Pa")
    ax[0].set_xlabel("u/U")
    ax[0].set_ylabel("y [mm]")
    ax[0].legend(fontsize=7)
    ax[0].set_title("Fig. 8.12 remake")
    ax[1].set_xlabel("u/U")
    ax[1].set_ylabel("y/(2√(νt))")
    ax[1].set_ylim(0, 2.5)
    ax[1].set_title("collapse (Fig. 8.13 axis)")
    # Crank–Nicolson
    t_end = 160.0
    Ns = (100, 200, 400) if args.fast else (100, 200, 400, 800)
    t0 = 10.0
    yg = np.linspace(0, 0.12, 2401)  # ≈ 9.5√(νt_end): the far condition (8.23) is felt < 1e-12
    errs, dts = [], []
    u0 = ch08.stokes_first_problem(yg, t0, U, nu)  # start from the exact profile at t0 (no wall jump)
    n_ref = 16 * Ns[-1]
    # smooth start → pure CN (startup_be=0); final profiles only (return_all=False)
    uref = DIF.crank_nicolson_1d(u0, yg, (t_end - t0) / n_ref, n_ref, nu, U, 0.0, startup_be=0, t0=t0)
    print(f"CN reference ({n_ref} steps) vs (8.30): {np.max(np.abs(uref - ch08.stokes_first_problem(yg, t_end, U, nu))):.2e}"
          " (spatial error floor)")
    for n in Ns:
        ucn = DIF.crank_nicolson_1d(u0, yg, (t_end - t0) / n, n, nu, U, 0.0, startup_be=0, t0=t0)
        errs.append(np.max(np.abs(ucn - uref)))  # temporal error on the same grid
        dts.append((t_end - t0) / n)
    p_t = np.polyfit(np.log(dts), np.log(errs), 1)[0]
    print(f"CN (Δy = {(yg[1] - yg[0]) * 1000:.2f} mm, t = {t0:.0f} → {t_end:.0f} s) errors: " + ", ".join(f"{e:.2e}" for e in errs) + f"; order {p_t:.3f}")
    ax[2].plot(ch08.stokes_first_problem(yg, t_end, U, nu), yg * 1000, color=COLORS["accent"], lw=2.4, label="(8.30)")
    ax[2].plot(ucn[::25], yg[::25] * 1000, "o", color=COLORS["orange"], ms=4, label="Crank–Nicolson")
    ax[2].set_xlabel("u/U")
    ax[2].set_ylabel("y [mm]")
    ax[2].set_ylim(0, 60)
    ax[2].legend(fontsize=8)
    ax[2].set_title(f"t = {t_end:.0f} s: CN order in Δt = {p_t:.2f}")
    # stopped plate (Exercise 8.30)
    T = 40.0
    for t in (20.0, 60.0, 120.0):
        ax[3].plot(ch08.stokes_first_stopped(y, t, T, U, nu), y * 1000, lw=2, label=f"t = {t:.0f} s")
    ax[3].set_xlabel("u/U")
    ax[3].set_ylabel("y [mm]")
    ax[3].legend(fontsize=8)
    ax[3].set_title(f"plate stopped at T = {T:.0f} s (Exercise 8.30)")
    save(fig, out, "c09_stokes_first")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
