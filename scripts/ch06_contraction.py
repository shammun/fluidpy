"""§6.7 (C12): finite-difference Laplace — the five-point rule (6.70)–(6.72) and its second-order accuracy, the
4-point system (6.73) solved exactly and by Gauss–Seidel, Example 6.2's contraction by the book's sweep order (values
after 50 sweeps and converged), Jacobi vs Gauss–Seidel vs SOR residual histories, flux conservation and grid
refinement (reduced order near the 270° corner).

Run: ``.venv/Scripts/python.exe scripts/ch06_contraction.py --no-show``
Figures → outputs/ch06/c12_contraction.png.
"""
from __future__ import annotations

import time

import numpy as np

from ch06_drawings import parse_args, save, setup

from fluidpy import ch06_ideal_flow as ch06
from fluidpy.core.style import COLORS


def book_inputs():
    """The book's flow rate and figure iteration count from the git-ignored private file (tests/book_values_ch06.json);
    without it: our Q = 1 m²/s and no fixed count (ψ scales linearly with Q)."""
    import json
    from pathlib import Path

    f = Path(__file__).resolve().parents[1] / "tests" / "book_values_ch06.json"
    try:
        d = json.loads(f.read_text(encoding="utf-8"))["ex6_2"]
        return float(d["flow_rate_per_depth_m2_s"]), int(d["iterations_quoted_for_figure"])
    except (OSError, KeyError, ValueError):
        return 1.0, None


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    errs, hs = [], []
    for nn in (11, 21, 41, 81):
        x = np.linspace(0, 1, nn)
        X, Y = np.meshgrid(x, x, indexing="xy")
        f = np.sin(np.pi * X) * np.sinh(np.pi * Y)
        lap = ch06.laplacian_5pt(f, x[1] - x[0])
        errs.append(np.nanmax(np.abs(lap)))
        hs.append(x[1] - x[0])
    print(f"five-point Laplacian of the harmonic sin(πx)sinh(πy): max error {['%.2e' % e for e in errs]}, "
          f"observed order {np.polyfit(np.log(hs), np.log(errs), 1)[0]:.3f}")
    bnd = {(1, 2): 1.0, (1, 3): 2.0, (2, 1): 0.0, (3, 1): 0.0, (4, 2): 1.0, (4, 3): 2.0, (2, 4): 3.0, (3, 4): 3.0}
    fp = ch06.four_point_system(bnd, sweeps=12)
    print(f"4-point system (6.73), boundary ψ = y: exact {fp['psi']}, after 12 Gauss–Seidel sweeps "
          f"{np.round(fp['gs_iterates'][-1], 8)}")

    Q, n_book = book_inputs()
    e50 = ch06.example_6_2(Q=Q, n_iter=n_book)
    ec = ch06.example_6_2(Q=Q)
    print(f"Example 6.2 (Δ = 1 m, Q = {Q:g} m²/s{' — the book value, from the private test file' if n_book else ''}):"
          f" after {n_book or 'no fixed count of'} book-order sweeps residual {e50['history']['residual'][-1]:.1e};"
          f" converged to 1e-10 in {ec['iterations']} sweeps; max change → converged "
          f"{np.nanmax(np.abs(e50['psi'] - ec['psi'])):.1e} m²/s")
    print(f"  interior ψ after {e50['iterations']} sweeps (rows y = 1…4 m, columns x = 1…8 m):")
    for j in range(1, 5):
        print("   ", " ".join(f"{v:5.2f}" if np.isfinite(v) else "  ·  " for v in e50["psi"][j, 1:9]))
    print(f"  flux across every interior vertical grid line: min {ec['flux'].min():.12f}, max {ec['flux'].max():.12f}"
          f" (= Q)")
    hist = {}
    for meth in ("jacobi", "gauss_seidel", "sor"):
        t0 = time.time()
        r = ch06.example_6_2(method=meth, refine=2)
        hist[meth] = r["history"]["residual"]
        print(f"  refine ×2, {meth:<12s}: {r['iterations']:5d} sweeps ({time.time() - t0:.2f} s)")
    probe = []
    refs = (1, 2, 4, 8) if not args.fast else (1, 2, 4)
    for rf in refs:
        r = ch06.example_6_2(method="direct", refine=rf)
        probe.append((r["psi"][2 * rf, 2 * rf], r["psi"][3 * rf, 7 * rf]))  # (x, y) = (2, 2) and (7, 3) m
    probe = np.array(probe)
    d = np.abs(np.diff(probe, axis=0))
    print("  refinement (direct): ψ(2 m, 2 m) =", np.round(probe[:, 0], 6), "; ψ(7 m, 3 m) =", np.round(probe[:, 1], 6))
    if len(d) >= 2:
        print(f"  successive-difference ratios (≈ 4 for order 2): {np.round(d[:-1, 0] / d[1:, 0], 3)}, "
              f"{np.round(d[:-1, 1] / d[1:, 1], 3)} (reduced near the re-entrant corner)")

    fig, ax = plt.subplots(1, 2, figsize=(14, 4.4))
    r4 = ch06.example_6_2(method="direct", refine=4)
    X, Y = np.meshgrid(r4["x"], r4["y"], indexing="xy")
    ax[0].contour(X, Y, r4["psi"], levels=np.linspace(0.1, 0.9, 9), colors=COLORS["ink"], linewidths=1.0)
    ax[0].fill([5, 9, 9, 5], [0, 0, 2, 2], color="#d9dce4")
    Xb, Yb = np.meshgrid(e50["x"], e50["y"], indexing="xy")
    ax[0].plot(Xb[e50["mask"]], Yb[e50["mask"]], ".", color=COLORS["orange"], ms=6, label="book grid unknowns")
    ax[0].set_aspect("equal")
    ax[0].set_xlabel("x [m]")
    ax[0].set_ylabel("y [m]")
    ax[0].set_title("Example 6.2 contraction: streamlines (grid ×4)", fontsize=10)
    ax[0].legend(fontsize=8, loc="upper right")
    for meth, col in (("jacobi", COLORS["muted"]), ("gauss_seidel", COLORS["accent"]), ("sor", COLORS["teal"])):
        ax[1].semilogy(hist[meth], color=col, label=meth.replace("_", "–"))
    ax[1].set_xlabel("sweep")
    ax[1].set_ylabel("max defect of the average rule [m²/s]")
    ax[1].set_title("convergence on the ×2 grid", fontsize=10)
    ax[1].legend(fontsize=8)
    save(fig, out, "c12_contraction")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
