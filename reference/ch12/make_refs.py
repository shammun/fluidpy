"""Write chapter 12's cited reference data into ``reference/ch12/``.

1. ``lee_moser_2015_channel_mean.csv`` — a small subset (y/δ, y⁺, U⁺, dU⁺/dy⁺) of the public mean-velocity profiles of the
   channel-flow direct numerical simulations of Lee & Moser, J. Fluid Mech. 774, 395–415 (2015), at Re_τ ≈ 180, 550, 1000,
   2000 and 5200.  The full files are downloaded (once) from the authors' server into the git-ignored
   ``data/online/ch12/`` and thinned here; nothing is typed by hand and nothing is interpolated — every row is a row of
   the source file, printed with 8 significant digits.
2. ``benchmarks.json`` — scalar published values with their citations (each verified on the date recorded in
   ``SOURCES.md``).  Nothing here comes from the textbook: the book's own numbers live only in the git-ignored
   ``tests/book_values_ch12.json``.  One entry is not a verified published value and says so: ``businger_dyer_unstable``
   records the coefficients the code uses (literature values, attribution not verified first-hand); the test that reads
   it is an analytic consistency check, not a benchmark.

``explainer_tables.json`` in this folder is *ours* (computed by ``fluidpy.ch12_turbulence.write_reference_tables`` through
``scripts/ch12_tables.py``); this script does not touch it.

Run: ``.venv/Scripts/python.exe reference/ch12/make_refs.py``  (idempotent; needs the network only the first time).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fluidpy.core import refdata  # noqa: E402

BASE = "https://turbulence.oden.utexas.edu/channel2015/data/"
CASES = ("0180", "0550", "1000", "2000", "5200")
N_LOG = 72            # points kept per profile, log-spaced in y+ (plus the wall point)
DENSE = (300.0, 800.0, 3)   # Re_tau = 5200 only: keep every 3rd source point in this y+ band (the logarithmic region)

BENCHMARKS = {
    "_note": "Published values for chapter 12 (sources and access dates in SOURCES.md). No textbook numbers. "
             "Exception: 'businger_dyer_unstable' holds the coefficients the code uses, not a verified published value.",
    "lee_moser_kappa": {"value": 0.384, "uncertainty": 0.004, "Re_tau": 5186,
                        "source": "Lee & Moser, J. Fluid Mech. 774, 395 (2015), abstract (arXiv:1410.7809)"},
    "kolmogorov_C1_one_sided": {"value": 0.53, "std": 0.055, "ci95": 0.03,
                                "source": "Sreenivasan, Phys. Fluids 7, 2778 (1995), as summarised by the ATOMIX wiki "
                                          "'Spectra in the inertial subrange'"},
    "kolmogorov_C1_over_C": {"value": 18.0 / 55.0, "C_implied": 1.5,
                             "source": "isotropic relation C1 = 18C/55 ~ 27/55 (ATOMIX wiki, same page; Pope 2000 §6.5)"},
    "k_epsilon_constants": {"C_mu": 0.09, "C_eps1": 1.44, "C_eps2": 1.92, "sigma_k": 1.0, "sigma_eps": 1.3,
                            "source": "Launder & Sharma (1974) / Launder & Spalding (1974); OpenFOAM LaunderSharmaKE "
                                      "defaults and SimScale k-epsilon documentation"},
    "prandtl_pipe_law": {"slope_log10": 2.0, "intercept": -0.8,
                         "source": "Prandtl (1935), quoted by McKeon, Zagarola & Smits, J. Fluid Mech. 538, 429 (2005)"},
    # NOT a verified published value: the coefficients the code uses (core/wall_turbulence.py).  The test that reads this
    # entry is an analytic consistency check (V1: the wind profile is the integral of the coded phi_m), not a V5 benchmark.
    "businger_dyer_unstable": {"coefficient": 16.0, "exponent": -0.25,
                               "kind": "coefficients as coded - analytic consistency check (V1), not a benchmark (V5)",
                               "verified_first_hand": False,
                               "source": "literature values ('Businger-Dyer' form, after Businger et al. 1971 and Dyer 1974); "
                                         "attribution not verified first-hand (no source was read; the AMS Glossary page "
                                         "returned HTTP 403) - see SOURCES.md"},
}


def read_case(tag: str) -> dict:
    """Download (cached) and parse one mean-profile file; return its header numbers and columns."""
    name = f"LM_Channel_{tag}_mean_prof.dat"
    path = refdata.fetch(BASE + name, ROOT / "data" / "online" / "ch12" / name)
    head = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("%"):
            break
        for key in ("nu", "u_tau", "Re_tau"):
            if f" {key} = " in line and key not in head:
                try:
                    head[key] = float(line.split("=")[-1])
                except ValueError:      # the citation line also contains "Re_tau = 5200," — not a parameter line
                    pass
    a = np.loadtxt(path, comments="%")
    return {"file": name, **head, "y_over_delta": a[:, 0], "yplus": a[:, 1], "Uplus": a[:, 2], "dUdy_plus": a[:, 3]}


def thin(case: dict, dense: bool) -> np.ndarray:
    """Indices of the rows kept: the wall point, ~N_LOG points log-spaced in y+, and (optionally) a denser band."""
    yp = case["yplus"]
    targets = np.geomspace(yp[1], yp[-1], N_LOG)
    idx = {0, len(yp) - 1}
    idx.update(int(np.argmin(np.abs(yp - t))) for t in targets)
    if dense:
        lo, hi, step = DENSE
        band = np.flatnonzero((yp >= lo) & (yp <= hi))
        idx.update(int(i) for i in band[::step])
    return np.array(sorted(idx))


def main() -> None:
    lines = [
        "# Subset of the mean-velocity profiles of Lee & Moser (2015), J. Fluid Mech. 774, 395-415, doi:10.1017/jfm.2015.268",
        "# source files: " + BASE + "LM_Channel_<case>_mean_prof.dat (downloaded by reference/ch12/make_refs.py)",
        "# every row is a row of the source file (thinned, not interpolated); wall units (u_tau, nu); half-channel",
        "# Re_tau = exact value from the file header; dUdy_plus = dU+/dy+ as tabulated by the authors",
        "case,Re_tau,y_over_delta,yplus,Uplus,dUdy_plus",
    ]
    meta = {}
    for tag in CASES:
        c = read_case(tag)
        keep = thin(c, dense=(tag == "5200"))
        meta[tag] = {"Re_tau": c["Re_tau"], "u_tau": c["u_tau"], "nu": c["nu"], "rows_source": int(len(c["yplus"])),
                     "rows_kept": int(len(keep)), "file": c["file"]}
        for i in keep:
            lines.append(f"{int(tag)},{c['Re_tau']:.3f},{c['y_over_delta'][i]:.8e},{c['yplus'][i]:.8e},"
                         f"{c['Uplus'][i]:.8e},{c['dUdy_plus'][i]:.8e}")
        print(f"  Re_tau {c['Re_tau']:9.3f}: kept {len(keep):3d} of {len(c['yplus'])} rows")
    (HERE / "lee_moser_2015_channel_mean.csv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    bm = dict(BENCHMARKS)
    bm["lee_moser_cases"] = meta
    (HERE / "benchmarks.json").write_text(json.dumps(bm, indent=1) + "\n", encoding="utf-8")
    print("wrote", HERE / "lee_moser_2015_channel_mean.csv")
    print("wrote", HERE / "benchmarks.json")


if __name__ == "__main__":
    main()
