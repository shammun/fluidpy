"""Verifier runs for chapter 10 (our computed data, not a cited benchmark): the MacCormack cavity at 128² that the
implementation phase did not cache as a table.

Writes ``reference/ch10/cavity_mck_re100_n128.csv`` in the format of ``scripts/ch10_cavity.py`` (header lines with the run
parameters, the deviation from Ghia et al. (1982) and the mass drift; columns y, u on the vertical centreline x = ½).
Parameters as the script's MacCormack runs: Re = 100, Ma = 0.08, t_end = 30, tol_steady = 1e-6, wall density from
continuity (10.139)–(10.146).  Cost ≈ 6–7 min on a laptop (the result is cached in outputs/ch10 by
``cavity_maccormack``, so a second call is instant).

Run: ``.venv/Scripts/python.exe reference/ch10/verify_runs.py``
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402

from fluidpy import ch10_computational_fluid_dynamics as ch10  # noqa: E402


def main() -> None:
    n = 128
    r = ch10.cavity_maccormack(Re=100.0, Ma=0.08, n=n, t_end=30.0, tol_steady=1e-6)
    e = ch10.cavity_error_vs_ghia(r, 100)
    u = np.asarray(r["u"])
    uc = u[:, n // 2]  # n even: the node column x = ½
    head = [f"# our weakly compressible MacCormack cavity (reference/ch10/verify_runs.py, verifier run 2026-09-30): Re = 100, "
            f"Ma = 0.08, {n}x{n} cells, wall density from continuity (10.139)-(10.146), dt = {r['dt']:.4g}, t = {r['t']:.3f}",
            f"# max deviation from Ghia 1982 = {e['max_dev']:.5f}; relative mass drift {r['mass_drift']:.3e}",
            "y,u"]
    rows = [f"{y:.6g},{v:.6g}" for y, v in zip(r["y"], uc)]
    (HERE / "cavity_mck_re100_n128.csv").write_text("\n".join(head + rows) + "\n", encoding="utf-8")
    print("wrote", HERE / "cavity_mck_re100_n128.csv", e)


if __name__ == "__main__":
    main()
