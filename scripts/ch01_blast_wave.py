"""Chapter 1, §1.11, Example 1.4: G. I. Taylor's blast-wave scaling E = K ρ D⁵/t² — recover the energy of a synthetic
explosion from noisy radius–time "photographs" (D ∝ t^{2/5}).

Run: ``.venv/Scripts/python.exe scripts/ch01_blast_wave.py --no-show``
Figure → outputs/ch01/fig1_11_blast_wave.png. The explosion is synthetic (E = 1e14 J, 2 % radius noise, seed 0);
K = 0.856 (Taylor 1950, spherical blast, γ = 1.4, as reported by Díaz, arXiv:2009.05674).
"""
from __future__ import annotations

import argparse
import sys
import time
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

    E_true, rho, K = 1.0e14, 1.2, ch01.TAYLOR_K_GAMMA14
    t = np.logspace(-4, -1.5, 15)  # [s]
    D = ch01.blast_radius(E_true, t, rho, K) * (1.0 + 0.02 * rng.standard_normal(t.size))
    slope, intercept = np.polyfit(np.log(t), np.log(D), 1)
    # with the exponent fixed at 2/5 (dimensional analysis), average E over the photographs
    E_est = np.mean(ch01.blast_energy(D, t, rho, K))
    fig, ax = plt.subplots()
    ax.loglog(t * 1e3, D, "o", color=COLORS["orange"], label="synthetic photographs (2 % noise)")
    tt = np.logspace(-4, -1.5, 100)
    ax.loglog(tt * 1e3, ch01.blast_radius(E_est, tt, rho, K), color=COLORS["accent"],
              label=f"D = (E t²/Kρ)^(1/5), E = {E_est:.3g} J")
    ax.set_xlabel("time t [ms]")
    ax.set_ylabel("blast radius D [m]")
    ax.set_title(f"Fitted slope {slope:.3f} (dimensional analysis: 2/5)")
    ax.legend()
    fig.savefig(out / "fig1_11_blast_wave.png", bbox_inches="tight")

    print(f"log-log slope of D(t): {slope:.4f} (exact 0.4); energy recovered {E_est:.4e} J vs true {E_true:.1e} J "
          f"({(E_est / E_true - 1) * 100:+.2f} %)")
    E_ground = np.mean(ch01.blast_energy(D, t, rho, K, geometry="hemisphere"))
    print(f"same radii read as a ground burst (hemisphere ≡ half of a free sphere of energy 2E): E = {E_ground:.4e} J")
    print(f"round trip blast_energy(blast_radius(E)) / E - 1 ={ch01.blast_energy(ch01.blast_radius(E_true, 0.01, rho), 0.01, rho) / E_true - 1:.1e}")
    print(f"saved fig1_11_blast_wave.png in {out}  ({time.perf_counter() - t_start:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
