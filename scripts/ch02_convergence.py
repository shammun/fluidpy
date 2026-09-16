"""Chapter 2, §2.9 and §2.12–2.13: observed order of the grid operators (gradient, divergence, curl, Laplacian; one-sided
edges and periodic) and of the integral definitions (2.31)–(2.33), (2.35) as the box/loop shrinks — the V3 evidence.

Run: ``.venv/Scripts/python.exe scripts/ch02_convergence.py --no-show [--fast]``
Figure → outputs/ch02/convergence.png.
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


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch02"))
    ap.add_argument("--no-show", action="store_true")
    ap.add_argument("--fast", action="store_true", help="n = 8, 16, 32 only")
    args = ap.parse_args()
    t0 = time.perf_counter()
    import matplotlib

    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from fluidpy import ch02_cartesian_tensors as ch02
    from fluidpy.core.style import CYCLE, use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    ns = (8, 16, 32) if args.fast else (8, 16, 32, 64)

    fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
    ax = axs[0]
    for k, op in enumerate(("gradient", "divergence", "curl", "laplacian", "vector_gradient")):
        c = ch02.operator_convergence(op, ns)
        ax.loglog(c["h"], c["err"], "o-", color=CYCLE[k % len(CYCLE)], label=f"{op} (one-sided edges): slope {c['order']:.2f}")
        print(f"{op:16s} onesided: n = {list(c['n'])}, max err = {np.round(c['err'], 6)}, order {c['order']:.3f}")
    c = ch02.operator_convergence("divergence", ns, bc="periodic")
    ax.loglog(c["h"], c["err"], "s--", color="k", label=f"divergence (periodic): slope {c['order']:.2f}")
    print(f"divergence periodic: max err = {np.round(c['err'], 6)}, order {c['order']:.3f}")
    ax.set_xlabel("grid spacing h")
    ax.set_ylabel("max |numerical − exact|")
    ax.set_title("Eqs. (2.22)–(2.25) on the grid: second order everywhere")
    ax.legend(fontsize=7.5)

    ax = axs[1]
    for k, kind in enumerate(("gradient", "divergence", "curl", "curl_component")):
        c = ch02.integral_definition_convergence(kind)
        eq = {"gradient": "(2.31)", "divergence": "(2.32)", "curl": "(2.33)", "curl_component": "(2.35)"}[kind]
        ax.loglog(c["h"], c["err"], "o-", color=CYCLE[k], label=f"{eq} {kind}: slope {c['order']:.2f}")
        print(f"integral {kind:15s} {eq}: h = {list(c['h'])}, err = {np.round(c['err'], 7)}, order {c['order']:.3f}")
    ax.set_xlabel("box / loop size h")
    ax.set_ylabel("|integral definition − exact derivative|")
    ax.set_title("integral definitions → differential forms as h → 0 (Examples 2.5, 2.6)")
    ax.legend(fontsize=8)
    fig.savefig(out / "convergence.png", bbox_inches="tight")
    print(f"saved convergence.png in {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
