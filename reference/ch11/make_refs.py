"""Write chapter 11's reference data into ``reference/ch11/``.

1. ``benchmarks.json`` — published benchmark values with their citations (typed from the sources listed in ``SOURCES.md``,
   each verified on 2026-09-30 as recorded there; analytic values computed).  Nothing here comes from the textbook — the
   book's own rounded numbers live only in the git-ignored ``tests/book_values_ch11.json``.
2. Our own computed tables (labelled "ours"): ``fluidpy.ch11_instability.write_reference_tables`` writes
   benard_neutral_curves.csv, taylor_critical.csv, tg_growth_map.csv, rayleigh_spectra.json, explainer_tables.json,
   critical_points.json and the Orr–Sommerfeld tables os_neutral_*.csv, os_grid_*.csv, os_modes.json.

Run: ``.venv/Scripts/python.exe reference/ch11/make_refs.py [--fast] [--no-tables] [--recompute]`` (idempotent; heavy pieces
are cached in ``outputs/ch11/cache`` and the critical points / TG map fall back to the published files here — after changing a
solver run with ``--recompute``).  Uncached full run ≈ 12 min.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

BENCHMARKS = {
    "_note": "Published benchmarks for chapter 11 (sources and access dates in SOURCES.md). Analytic values computed.",
    "benard_rigid_rigid": {"Ra_c": 1707.762, "K_c": 3.117, "source": "Chandrasekhar (1961), Hydrodynamic and Hydromagnetic "
                           "Stability; quoted in the Nek5000 examples and arXiv:nlin/0302057"},
    "benard_free_free": {"Ra_c": 27 * math.pi ** 4 / 4, "K_c": math.pi / math.sqrt(2), "source": "analytic (Rayleigh 1916)"},
    "benard_rigid_free": {"Ra_c": 1100.65, "K_c": 2.682, "source": "Chandrasekhar (1961) (primary not fetched; our "
                          "Chebyshev value agrees to 6 digits)"},
    "benard_odd_mode": {"Ra_c": 17610.39, "K_c": 5.365, "source": "Chandrasekhar (1961), p. 39"},
    "taylor_narrow_gap_mu1": {"Ta_c": 1707.76, "k_c": 3.117, "source": "narrow-gap limit = rigid Benard; arXiv:2601.14806"},
    "taylor_narrow_gap_corotation": {"Ta_c_mu1": 1707.76, "coefficient": 0.00761,
                                     "formula": "Ta_c = 1707.76 [1 - 0.00761 ((1 - mu)/(1 + mu))^2] as mu -> 1, with "
                                                "Ta = -2 A Omega_1 d^4 (1 + mu)/nu^2 = (1 + mu)/2 x the Ta of Eq. (11.52)",
                                     "source": "Wikipedia 'Taylor-Couette flow', section 'Taylor's criterion' (revision "
                                               "1373953105 of 2026-09-08), read 2026-10-01; asymptotic (leading order in "
                                               "1 - mu), tertiary source"},
    "plane_poiseuille_critical": {"Re_c": 5772.22, "k_c": 1.02056, "c_r": 0.264, "source": "Orszag (1971), J. Fluid Mech. 50, 689"},
    "plane_poiseuille_Re1e4_k1": {"c_r": 0.23752649, "c_i": 0.00373967, "source": "Orszag (1971); eigentools documentation "
                                  "(max c_i = 3.740e-3)"},
    "blasius_critical": {"Re_delta_star": 519.2, "alpha_delta_star": 0.303, "omega": 0.120,
                         "source": "Thomas, via Gallagher, Griffiths & Stephen, Phys. Fluids 28, 074107 (2016); Jordinson "
                                   "(1970): 520"},
    "bickley_jet_sinuous": {"Re_c": 4.0, "k_c": 0.2, "source": "Tatsumi & Kakutani (1958), J. Fluid Mech. 4, 261 "
                            "(approximate, via later papers)"},
    "tanh_shear_layer_inviscid": {"k_max": 0.4446, "kci_max": 0.1897, "neutral_k": 1.0,
                                  "source": "Michalke (1964), J. Fluid Mech. 19, 543 (digits approximate); neutral mode "
                                            "analytic"},
    "bickley_inviscid_neutral": {"sinuous_k": 2.0, "varicose_k": 1.0, "c": 2.0 / 3.0, "source": "analytic (Drazin & Reid 1981)"},
    "piecewise_shear_layer_neutral_kh": {"kh": 1.278465, "source": "root of (kh - 1)^2 = exp(-2kh) (Rayleigh 1880)"},
    "miles_howard_Ri": {"Ri": 0.25, "source": "Miles (1961); Howard (1961)"},
    "tanh_sech2_neutral_curve": {"formula": "J = k(1 - k)", "source": "exact neutral mode |tanh z|^(1-k) sech^k z with c = 0 "
                                 "(checked numerically here; Drazin & Reid 1981)"},
    "lorenz": {"Pr": 10, "b": 8.0 / 3.0, "r": 28, "r_H": 10 * (10 + 8 / 3 + 3) / (10 - 8 / 3 - 1),
               "source": "Lorenz (1963); Wikipedia 'Lorenz system' (r_H = 24.74)"},
    "feigenbaum": {"delta": 4.669201609, "alpha": 2.502907875, "A_2": 1 + math.sqrt(6), "A_3": 3.5440903, "A_4": 3.5644073,
                   "A_inf": 3.5699456, "source": "Feigenbaum (1978); Wikipedia 'Feigenbaum constants'"},
}


def main() -> int:
    ap = argparse.ArgumentParser(description="write reference/ch11")
    ap.add_argument("--fast", action="store_true")
    ap.add_argument("--no-tables", action="store_true")
    ap.add_argument("--recompute", action="store_true",
                    help="ignore every cache and the published tables (needed after changing a solver; ~12 min)")
    args = ap.parse_args()
    (HERE / "benchmarks.json").write_text(json.dumps(BENCHMARKS, indent=1), encoding="utf-8")
    print("wrote", HERE / "benchmarks.json")
    if args.no_tables:
        return 0
    from fluidpy import ch11_instability as ch11

    t0 = time.time()
    paths = ch11.write_reference_tables(HERE, fast=args.fast, cache=not args.recompute)
    for k, v in paths.items():
        print(f"  {k:28s} {v}")
    print(f"tables written in {time.time() - t0:.1f} s")
    crit = json.loads((HERE / "critical_points.json").read_text(encoding="utf-8"))
    print("critical points (ours):", json.dumps({k: v for k, v in crit.items() if k != "note"}, indent=None)[:1500])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
