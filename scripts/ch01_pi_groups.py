"""Chapter 1, §1.11: the dimensional matrix (1.39), its rank by minors, the Π groups (1.37)/(1.40) for the pipe and
Examples 1.2–1.5, invariance under a change of units, and the collapse of synthetic laminar pipe data.

Run: ``.venv/Scripts/python.exe scripts/ch01_pi_groups.py --no-show``
Figures → outputs/ch01/fig_pipe_matrix.png, fig_pi_units_invariance.png, fig_pipe_collapse.png.
Synthetic data: seeded (numpy default_rng(0)); laminar law Δp = 32 μ U Δx/d² (derived in Ch. 8).
"""
from __future__ import annotations

import argparse
import sys
import time
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch01"), help="figure folder")
    ap.add_argument("--no-show", action="store_true", help="do not open a window")
    args = ap.parse_args()
    t_start = time.perf_counter()
    import matplotlib

    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from fluidpy import ch01_introduction as ch01
    from fluidpy.core.style import COLORS, use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)

    # --- the pipe matrix with its minors -------------------------------------------------------------------------------------
    A, names, rows = ch01.dimensional_matrix(ch01.PIPE)
    r, wr, wc = ch01.rank_by_minors(A)
    labels = [f"${ch01.LATEX_SYMBOLS.get(n, n)}$" for n in names]
    fig, ax = plt.subplots(figsize=(7, 3))
    ax.imshow(A, cmap="PuOr", vmin=-3, vmax=3)
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            ax.text(j, i, str(A[i, j]), ha="center", va="center", fontsize=12)
    ax.add_patch(plt.Rectangle((-0.45, -0.45), 2.9, 2.9, fill=False, ec=COLORS["rose"], lw=2))
    ax.add_patch(plt.Rectangle((3.55, -0.45), 2.9, 2.9, fill=False, ec=COLORS["teal"], lw=2))
    ax.set_xticks(range(len(names)), labels)
    ax.set_yticks(range(len(rows)), rows)
    ax.set_title(f"Eq. (1.39): minors {ch01.minor_determinant(A, [0, 1, 2], [0, 1, 2])} (red) and "
                 f"{ch01.minor_determinant(A, [0, 1, 2], [4, 5, 6])} (green) → rank {r}")
    ax.grid(False)
    fig.savefig(out / "fig_pipe_matrix.png", bbox_inches="tight")

    # --- groups for every preset --------------------------------------------------------------------------------------------
    print(f"pipe matrix rows {rows}, rank {r} (witness rows {wr}, cols {wc}), numpy rank {np.linalg.matrix_rank(A)}")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)  # kmol dropped in the scale-height preset (Example 1.2)
        for key, info in ch01.PRESET_INFO.items():
            groups = ch01.pi_groups(info["variables"], info["solution"], info["repeating"])
            Ak, _, _ = ch01.dimensional_matrix(info["variables"])
            n, rk = len(info["variables"]), ch01.rank_by_minors(Ak)[0]
            print(f"{key:12s} n = {n}, r = {rk}, n − r = {n - rk}: " + ",  ".join(ch01.group_latex(gq) for gq in groups)
                  + f"   independent: {ch01.groups_independent(groups)}")
    g_pipe = ch01.pi_groups(ch01.PIPE, "dp", ("U", "d", "rho"))
    print("exponents of Π1 (U, d, rho):", {k: str(v) for k, v in g_pipe[0].items()})
    combined = {"dp": 1, "d": 2, "rho": 1, "mu": -2}  # Π1/Π4^2
    print("{Π1, Π4, Δp d² ρ/μ²} independent?", ch01.groups_independent([g_pipe[0], g_pipe[3], combined]))

    # --- units invariance ------------------------------------------------------------------------------------------------------
    samples = [dict(dp=10 ** rng.uniform(0, 5), dx=rng.uniform(0.1, 10), d=rng.uniform(0.005, 0.5),
                    eps=rng.uniform(1e-6, 1e-3), U=rng.uniform(0.01, 5), rho=rng.uniform(1, 1000),
                    mu=10 ** rng.uniform(-5, -1)) for _ in range(40)]
    fig, ax = plt.subplots()
    worst = 0.0
    for system, c in (("cgs", COLORS["teal"]), ("imperial", COLORS["orange"])):
        si_vals, new_vals = [], []
        for s in samples:
            conv = ch01.rescale_units(s, ch01.PIPE, system)
            for gq in g_pipe:
                si_vals.append(ch01.group_value(gq, s))
                new_vals.append(ch01.group_value(gq, conv))
        si_vals, new_vals = np.array(si_vals), np.array(new_vals)
        worst = max(worst, float(np.max(np.abs(new_vals / si_vals - 1))))
        ax.loglog(si_vals, new_vals, "o", ms=4, color=c, label=f"groups computed in {system} units")
    lim = [1e-8, 1e8]
    ax.loglog(lim, lim, color=COLORS["muted"], lw=0.8)
    ax.set_xlabel("Π computed in SI units")
    ax.set_ylabel("Π computed in another unit system")
    ax.set_title("A dimensionless group has the same value in every unit system")
    ax.legend()
    fig.savefig(out / "fig_pi_units_invariance.png", bbox_inches="tight")

    # --- laminar pipe data collapse --------------------------------------------------------------------------------------------
    N = 200
    d = rng.uniform(0.001, 0.02, N)
    U = rng.uniform(0.001, 0.05, N)
    dx = rng.uniform(0.1, 5.0, N)
    mu = 10 ** rng.uniform(-3, 0, N)
    rho = rng.uniform(800, 1300, N)
    dp = ch01.poiseuille_pressure_drop(mu, U, dx, d)
    Re = rho * U * d / mu
    keep = Re < 2000
    P1 = dp / (rho * U ** 2)
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    axs[0].loglog(dx[keep] / d[keep], P1[keep], "o", ms=3, color=COLORS["rose"])
    axs[0].set_xlabel("Δx/d (μ left out of the variable list)")
    axs[0].set_ylabel("Δp/(ρU²)")
    axs[0].set_title("Missing a variable: no collapse")
    x = (dx / d) * (mu / (rho * U * d))
    axs[1].loglog(x[keep], P1[keep], "o", ms=3, color=COLORS["teal"], label="synthetic laminar data")
    xs = np.logspace(np.log10(x[keep].min()), np.log10(x[keep].max()), 10)
    axs[1].loglog(xs, 32 * xs, color=COLORS["ink"], lw=1, label="Π1 = 32 Π2 Π4")
    axs[1].set_xlabel("(Δx/d)·μ/(ρUd)")
    axs[1].set_title("With the right groups every pipe falls on one line")
    axs[1].legend()
    fig.savefig(out / "fig_pipe_collapse.png", bbox_inches="tight")

    slope = np.polyfit(np.log(x[keep]), np.log(P1[keep]), 1)
    print(f"units invariance: worst relative change of any group = {worst:.1e}")
    print(f"laminar collapse: {keep.sum()} of {N} samples with Re < 2000; log-log fit slope {slope[0]:.6f}, "
          f"prefactor {np.exp(slope[1]):.4f}")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        g_sh = ch01.pi_groups(ch01.SCALE_HEIGHT, "H", ("T0", "Mw", "g", "Ru"))[0]
    vals = dict(H=ch01.scale_height(250.0), T0=250.0, Mw=ch01.M_W_AIR, g=ch01.G0, Ru=ch01.R_U)
    print(f"Example 1.2: Π1 = {ch01.group_latex(g_sh)} = {ch01.group_value(g_sh, vals):.12f} (constant = 1)")
    print(f"saved 3 figures in {out}  ({time.perf_counter() - t_start:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
