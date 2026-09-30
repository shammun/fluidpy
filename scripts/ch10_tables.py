"""Tables of our own computed data for the Chapter 10 explainers and plotly figures (≤ 4 significant figures, labelled "ours"):
reference/ch10/infsup_table.csv (β_h of P1–P1 and P2–P1 on the unit square, F7/E8) and reference/ch10/explainer_tables.json
(cavity centrelines and vortex centres from the cached runs of scripts/ch10_cavity.py for E7; the inf–sup table for E8).

Run: ``.venv/Scripts/python.exe scripts/ch10_tables.py --no-show`` (after scripts/ch10_cavity.py for the cavity tables).
"""
from __future__ import annotations

import json

import numpy as np

from ch10_common import REF, parse_args, setup, write_csv

from fluidpy import ch10_computational_fluid_dynamics as ch10


def _sig(x, n=4):
    return float(f"{float(x):.{n}g}")


def _read(path):
    header = [ln[2:].strip() for ln in open(path, encoding="utf-8") if ln.startswith("#")]
    d = np.loadtxt(path, delimiter=",", comments="#", skiprows=len(header) + 1)
    return header, d


def main() -> int:
    args = parse_args(__doc__)
    setup(args)
    n_list = (2, 3, 4, 6, 8, 12) if args.fast else (2, 3, 4, 6, 8, 12, 16)
    tab = ch10.infsup_table(n_list=n_list)
    rows = [(0 if p == "P1P1" else 1, r["n"], r["h"], r["beta"], r["beta_nonspurious"], r["n_spurious"], r["n_u"], r["n_p"])
            for p, rs in tab.items() for r in rs]
    a = np.array(rows, dtype=float)
    write_csv(REF / "infsup_table.csv",
              ["our discrete inf-sup constants on the unit square (scripts/ch10_tables.py): structured n x n squares cut along one "
               "diagonal, velocity Dirichlet on the whole boundary; pair 0 = P1-P1, 1 = P2-P1; beta = 0 when spurious pressure "
               "modes exist (beta_nonspurious = the smallest nonzero one)"],
              dict(pair=a[:, 0], n=a[:, 1], h=a[:, 2], beta=a[:, 3], beta_nonspurious=a[:, 4], n_spurious=a[:, 5],
                   n_u=a[:, 6], n_p=a[:, 7]))
    for p, rs in tab.items():
        print(p, [(r["n"], round(r["beta"], 4), r["n_spurious"]) for r in rs])
    js = dict(note="ours (fluidpy), 4 significant figures; Ghia et al. (1982) and Hou et al. (1995) are cited public data",
              infsup={p: [dict(n=r["n"], beta=_sig(r["beta"]), beta_nonspurious=_sig(r["beta_nonspurious"]),
                               n_spurious=r["n_spurious"]) for r in rs] for p, rs in tab.items()},
              cavity={}, ghia={str(Re): dict(y=[_sig(v) for v in ch10.ghia_centreline(Re)["y"]],
                                              u=[_sig(v) for v in ch10.ghia_centreline(Re)["u"]]) for Re in (100, 400)},
              ghia_centre={str(Re): ch10.ghia_vortex_centre(Re) for Re in (100, 400)},
              hou_centre={str(k): list(v) for k, v in ch10.hou_centres().items()})
    for f in sorted(REF.glob("cavity_*_re*_n*.csv")):
        header, d = _read(f)
        key = f.stem.replace("cavity_", "")
        js["cavity"][key] = dict(header=header, y=[_sig(v) for v in d[:, 0]], u=[_sig(v) for v in d[:, 1]])
    (REF / "explainer_tables.json").write_text(json.dumps(js, indent=1), encoding="utf-8")
    print(f"  wrote {REF / 'explainer_tables.json'} with cavity runs {sorted(js['cavity'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
