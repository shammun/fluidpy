"""Write chapter 12's computed tables (ours, computed by fluidpy — no book data, no external benchmark) to
``reference/ch12/explainer_tables.json``:

* ``channel``: the two energy budgets of the mixing-length channel (U+, Reynolds stress, production, pressure work, mean
  dissipation, turbulence sink, integrated budget) at Re_tau = 180, 550, 1000, 5200 — the table of the energy-budget
  explainer (E5) and of the mixing-length explainer's channel view (E8);
* ``van_driest``: the wall profile U+(y+) of the mixing-length model with and without van Driest damping (E8).

The script then reads the file back and proves that it reproduces ``ch12.channel_energy_budget_at`` and
``ch12.mixing_length_wall_profile`` (so the JSON on disk is regenerated, not hand-edited).  It also prints the conventions
table and the list of printed slips.

Run: ``.venv/Scripts/python.exe scripts/ch12_tables.py --no-show [--dest FOLDER]``
"""
from __future__ import annotations

import json

import numpy as np
from ch12_common import REF, Timer, finish, parse_args, setup

from fluidpy import ch12_turbulence as ch12


def main() -> int:
    args = parse_args(__doc__, extra=lambda ap: ap.add_argument("--dest", default=str(REF)))
    setup(args)
    with Timer("write_reference_tables"):
        paths = ch12.write_reference_tables(args.dest)
    for k, v in paths.items():
        print(f"  {k:20s} {v}")
    tab = json.loads(paths["explainer_tables"].read_text(encoding="utf-8"))
    kappa, A_plus = tab["kappa"], tab["A_plus"]
    print(f"  note: {tab['note']}; kappa = {kappa}, A+ = {A_plus}")
    worst = 0.0
    for Re, row in tab["channel"].items():
        print(f"  Re_tau = {Re:>5s}: U_bulk+ = {row['U_bulk_plus']}, U_cl+ = {row['U_cl_plus']}, Cf = {row['Cf']}, production peak "
              f"{row['production_peak']} at y+ = {row['yplus_peak_production']}, work = {row['integrals']['work']} = dissipation "
              f"{row['integrals']['dissipation']} + production {row['integrals']['production']}")
        for i in (10, 40, 70, len(row["yplus"]) - 1):   # parity of the stored rows with the scalar function
            b = ch12.channel_energy_budget_at(row["yplus"][i], float(Re), kappa, A_plus)
            for key, col in (("Uplus", "Uplus"), ("production", "production"), ("mean_dissipation", "mean_dissipation"),
                             ("pressure_work", "pressure_work"), ("uv_plus", "uv_plus"), ("slope", "dUdy_plus")):
                ref = abs(b[key]) if b[key] != 0 else 1.0
                worst = max(worst, abs(row[col][i] - b[key]) / ref)
    vd = tab["van_driest"]
    for i in (0, len(vd["yplus"]) // 2, len(vd["yplus"]) - 1):
        p = ch12.mixing_length_wall_profile(vd["yplus"][i], kappa, damping="van_driest", A_plus=A_plus)
        worst = max(worst, abs(vd["Uplus"][i] - p["Uplus"]) / p["Uplus"])
    print(f"  van Driest table: {len(vd['yplus'])} points, y+ = {vd['yplus'][0]} ... {vd['yplus'][-1]}; intercept B = {vd['B']} "
          f"(undamped {vd['B_undamped']}); U+(y+ = {vd['yplus'][-1]:.0f}) - ln(y+)/kappa = "
          f"{vd['Uplus'][-1] - np.log(vd['yplus'][-1]) / kappa:.4f}")
    print(f"  stored rows against ch12.channel_energy_budget_at / mixing_length_wall_profile: worst relative difference {worst:.1e} "
          f"(6 significant digits are stored)")
    assert worst < 2e-5, "explainer_tables.json does not reproduce the functions it was written from"
    conv = ch12.conventions()
    print(f"conventions table: {conv.shape[0]} rows; first three:")
    print(conv.head(3).to_string(index=False))
    print(f"printed slips carried as corrections: {len(ch12.book_slips())}")
    for s in ch12.book_slips():
        print(f"  #{s['id']:<2d} {s['where']:20s} {s['printed']}  ->  {s['corrected']}   [test: {s['test']}]")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
