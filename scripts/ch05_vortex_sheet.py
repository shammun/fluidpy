"""§5.8 vortex sheet (C14, N44, N60): the tangential jump of a row of filaments converging to a continuous sheet (L1
error ∝ 1/N), the sign convention of the strength (text u₂ − u₁ vs caption u₁ − u₂), the circuit's circulation, a
from-scratch row sum vs ``vortex_sheet_velocity``, and the Krasny-smoothed roll-up of a rippled periodic sheet.

Run: ``.venv/Scripts/python.exe scripts/ch05_vortex_sheet.py --no-show``
Figures → outputs/ch05/c14_sheet_profiles.png, c14_rollup.png.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from fluidpy import ch05_vorticity_dynamics as ch05  # noqa: E402
from fluidpy.core.style import COLORS  # noqa: E402
from tools.convergence import observed_order  # noqa: E402


def profiles_figure(dc):
    import matplotlib.pyplot as plt

    y = np.linspace(-0.1, 0.1, 801)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    for N, c in zip((10, 100, 1000), (COLORS["rose"], COLORS["teal"], COLORS["accent"])):
        ax[0].plot(ch05.discrete_sheet_u(0.0, y, 2.0, N), y, color=c, label=f"N = {N}")
    ax[0].plot(ch05.continuous_sheet_velocity(0.0, y, 2.0)[0], y, "--", color=COLORS["ink"], label="continuous sheet")
    ax[0].set_xlabel("u(0, y) [m/s]")
    ax[0].set_ylabel("y [m]")
    ax[0].legend(fontsize=8)
    ax[0].set_title("γ = 2 m/s: u → ∓γ/2 just above/below")
    ax[1].loglog(dc["N"], dc["l1_error"], "o-", color=COLORS["accent"],
                 label=f"observed order {observed_order(1.0 / dc['N'], dc['l1_error']):.2f}")
    ax[1].set_xlabel("number of filaments N")
    ax[1].set_ylabel("L1 error of u across |y| < 0.1 m [m²/s]")
    ax[1].legend(fontsize=8)
    return fig


def rollup_figure(fast: bool):
    import matplotlib.pyplot as plt

    te = np.linspace(0.0, 4.0, 5)
    sr = ch05.sheet_rollup(N=100 if fast else 200, t_eval=te)
    fig, ax = plt.subplots(1, len(te), figsize=(18, 3.2), sharey=True)
    for i, a in enumerate(ax):
        a.plot(sr["x"][i], sr["y"][i], ".-", color=COLORS["accent"], ms=2, lw=0.6)
        a.set_title(f"t = {te[i]:.1f} s", fontsize=9)
        a.set_xlabel("x [m]")
        a.set_aspect("equal")
    ax[0].set_ylabel("y [m]")
    return fig, sr


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch05"))
    ap.add_argument("--no-show", action="store_true")
    ap.add_argument("--fast", action="store_true")
    args = ap.parse_args()
    import matplotlib

    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from fluidpy.core.style import use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    ua, ub = ch05.discrete_sheet_u(0.0, 1e-3, 2.0, 1000), ch05.discrete_sheet_u(0.0, -1e-3, 2.0, 1000)
    print(f"row of 1000 counterclockwise filaments, γ = 2 m/s: u just above {ua:+.4f}, below {ub:+.4f} m/s → strength "
          f"(text, ccw) {ch05.vortex_sheet_strength(ua, ub):+.4f}, (caption) "
          f"{ch05.vortex_sheet_strength(ua, ub, 'caption'):+.4f} m/s")
    print(f"continuous sheet u(0, 0.05) = {ch05.continuous_sheet_velocity(0.0, 0.05, 2.0)[0]:.6f}; row N = 100: "
          f"{ch05.discrete_sheet_u(0.0, 0.05, 2.0, 100):.6f} m/s")
    # circulation of the dn × ds circuit around the middle of the sheet (counterclockwise)
    S, ds_ = np.stack([-0.5 + (np.arange(1000) + 0.5) / 1000, np.zeros(1000)]), 1e-3
    u_fn = lambda x, t=0.0: ch05.vortex_sheet_velocity(x, S, 2.0, ds_)  # noqa: E731
    circuit = ch05.square_loop_points((0.0, 0.0), 0.1, 2048)
    print(f"circuit 0.1 × 0.1 m around the sheet: Γ = {ch05.loop_circulation(u_fn, circuit):.6f} m²/s = γ·ds = "
          f"{2.0 * 0.1:.6f} (dΓ = (u₂ − u₁)ds)")
    # from scratch: the row sum by hand
    P = np.array([0.013, 0.04])
    by_hand = np.zeros(2)
    for j in range(S.shape[1]):
        rx, ry = P[0] - S[0, j], P[1] - S[1, j]
        by_hand += 2.0 * ds_ / (2 * np.pi) * np.array([-ry, rx]) / (rx * rx + ry * ry)
    print(f"from scratch row sum {by_hand} vs vortex_sheet_velocity {u_fn(P)}")
    dc = ch05.discrete_sheet_convergence(2.0, (10, 30, 100, 300, 1000) if not args.fast else (10, 100, 1000))
    for N, e in zip(dc["N"], dc["l1_error"]):
        print(f"   N = {N:4d}: L1 error {e:.6e} m²/s")
    print(f"observed order in 1/N: {observed_order(1.0 / dc['N'], dc['l1_error']):.3f}")
    profiles_figure(dc).savefig(out / "c14_sheet_profiles.png", bbox_inches="tight")
    fig, sr = rollup_figure(args.fast)
    fig.savefig(out / "c14_rollup.png", bbox_inches="tight")
    print(f"roll-up (N = {sr['N']}, δ = 0.05): mean y conserved to {np.ptp(sr['mean_y']):.1e} m; max |y| at t = 4 s: "
          f"{np.abs(sr['y'][-1]).max():.3f} m")
    print(f"figures → {out}")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
