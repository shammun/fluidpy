"""Write chapter 11's computed tables (ours) to ``reference/ch11/``: the Bénard neutral table, the Taylor critical table and
eigenfunctions, the Taylor–Goldstein growth map, Rayleigh spectra, the Lorenz r-sweep, ``explainer_tables.json`` and
``critical_points.json``, plus (unless --no-os) the Orr–Sommerfeld tables of the E8 explainer.

Run: ``.venv/Scripts/python.exe scripts/ch11_tables.py --no-show [--fast] [--no-os]``  (same as reference/ch11/make_refs.py
without the benchmark file).  Heavy pieces are cached in outputs/ch11/cache; the first full run takes ~10–15 min.
"""
from __future__ import annotations

import json

from ch11_common import REF, Timer, finish, parse_args, setup

from fluidpy import ch11_instability as ch11


def main() -> int:
    args = parse_args(__doc__, extra=lambda ap: ap.add_argument("--no-os", action="store_true"))
    setup(args)
    with Timer("write_reference_tables"):
        paths = ch11.write_reference_tables(REF, fast=args.fast, include_os=not args.no_os)
    for k, v in paths.items():
        print(f"  {k:26s} {v}")
    crit = json.loads((REF / "critical_points.json").read_text(encoding="utf-8"))
    for k, v in crit.items():
        if k != "note":
            print(f"  {k}: {v}")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
