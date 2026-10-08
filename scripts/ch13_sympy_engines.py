"""Chapter 13 — the twelve symbolic engines: each re-derives one result of the chapter with sympy and reports a residual
that must be zero.  Two of them can also be run in the form the book prints (``printed=True``) and must then fail:
slip #1 (the sign of the pressure term in the rotating Boussinesq momentum equation) and slip #3 (the sign in the
second tank equation, which breaks the Taylor–Proudman argument).  Also prints the slips and traps tables.

Run: ``.venv/Scripts/python.exe scripts/ch13_sympy_engines.py --no-show``
"""
from __future__ import annotations

import time

from ch13_common import finish, parse_args, setup

from fluidpy import ch13_geophysical_fluid_dynamics as ch13


def main() -> int:
    args = parse_args(__doc__)
    setup(args)
    bad = 0
    for name, fn in ch13.sympy_engines().items():
        t0 = time.time()
        r = fn()
        bad += not r["ok"]
        print(f"  {name:28s} ok = {r['ok']!s:5s} residual = {r['residual']}   [{time.time() - t0:.2f} s]")
    r1 = ch13.boussinesq_rotating_sympy(printed=True)
    r3 = ch13.taylor_proudman_sympy(printed=True)
    print(f"printed forms must fail: slip #1 rest-state residual = {r1['residual']} (ok = {r1['ok']}); slip #3 horizontal divergence = "
          f"{r3['residual']} (ok = {r3['ok']}; what the printed pair gives instead: du/dx - dv/dy = {r3['printed_identity']})")
    d_ok, d_bad = ch13.friction_force_dimensions(), ch13.friction_force_dimensions(printed=True)
    print(f"slip #2: F = (1/rho) d(tau)/dx has dimensions M^{d_ok['M']} L^{d_ok['L']} T^{d_ok['T']} (an acceleration: {d_ok['is_acceleration']}); "
          f"the printed form M^{d_bad['M']} L^{d_bad['L']} T^{d_bad['T']} (an acceleration: {d_bad['is_acceleration']})")
    v = ch13.v_equation_sympy()
    print(f"the only approximation in the v-equation: f² -> f0², relative size {v['neglected_relative']}")
    print(f"operator of the w-equation: {ch13.rotating_internal_wave_set_sympy()['operator']}")
    aud = ch13.reexport_audit()
    print(f"re-export audit: {aud['total']} public names in the three core modules, {aud['same_object']} reachable as ch13.<name>, "
          f"aliased {aud['aliased']}, missing {aud['missing']}")
    slips = ch13.book_slips()
    print(f"book slips carried as corrections: {len(slips)}")
    for key, s in slips.items():
        print(f"  #{key:>2s} {s['where']:42s} test: {s['test']}")
    print(f"traps kept apart: {len(ch13.traps())}; conventions table: {ch13.conventions_table().shape[0]} rows")
    assert bad == 0 and not r1["ok"] and not r3["ok"] and not aud["missing"]
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
