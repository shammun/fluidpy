"""Write chapter 13's cached model runs and computed constants to ``reference/ch13`` — **our own output only** (no book
data, no external benchmark):

* ``kelvin_basin.npz``     a Kelvin pulse going round a closed basin (C-grid shallow-water model, 64²)
* ``turbulence_f.npz``     decaying two-dimensional turbulence, β = 0 (pseudo-spectral, 64²)
* ``turbulence_beta.npz``  the same initial field on a β-plane (zonal bands)
* ``pv_particles.npz``     a nonlinear geostrophic vortex with marked particles and the potential vorticity they carry
* ``explainer_constants.json``  constants the explainers compute themselves, for their parity rows

Frames are stored as float16.  After writing, each file is read back with ``ch13.load_reference_run`` and compared
with the run that produced it (to the float16 quantisation), and the constants are compared with the functions they
came from — so the files on disk are regenerated, never hand-edited.

Run: ``.venv/Scripts/python.exe scripts/ch13_make_caches.py --no-show [--which all|name ...] [--dest FOLDER]``
"""
from __future__ import annotations

import hashlib
import json

import numpy as np
from ch13_common import REF, Timer, finish, parse_args, setup

from fluidpy import ch13_geophysical_fluid_dynamics as ch13


def _extra(ap):
    ap.add_argument("--which", nargs="+", default=["all"], help="'all' or names from ch13.REFERENCE_RUNS")
    ap.add_argument("--dest", default=str(REF))


def main() -> int:
    args = parse_args(__doc__, extra=_extra)
    setup(args)
    which = "all" if args.which == ["all"] else args.which
    with Timer("compute and write"):
        written = ch13.write_reference_runs(which, args.dest)
    total = 0
    for name, info in written.items():
        total += info["bytes"]
        if name == "explainer_constants":
            back = json.loads(info["path"].read_text(encoding="utf-8"))
            assert abs(back["eady_critical"] - ch13.eady_critical()) < 1e-12
            assert abs(back["eady_fastest"]["sigma_nd"] - ch13.eady_fastest()["sigma_nd"]) < 1e-12
            assert abs(back["adjustment_ratio"] - 1.0 / 3.0) < 1e-12
            print(f"  {name:20s} {info['bytes']:8d} bytes: alpha_c H = {back['eady_critical']:.5f}, fastest alpha H = "
                  f"{back['eady_fastest']['alphaH']:.5f}, growth {back['eady_fastest']['sigma_nd']:.5f}, adjustment ratio "
                  f"{back['adjustment_ratio']:.6f} — equal to the functions to 1e-12")
            continue
        fresh, back = info["fresh"], ch13.load_reference_run(name, args.dest)
        assert back is not None and set(back) == set(fresh)
        worst = 0.0
        for key, val in fresh.items():
            a, b = np.asarray(val, dtype=float), np.asarray(back[key], dtype=float)
            scale = max(float(np.max(np.abs(a))), 1e-300)
            worst = max(worst, float(np.max(np.abs(a - b))) / scale)
        digest = hashlib.sha256(b"".join(np.ascontiguousarray(back[k]).tobytes() for k in sorted(back))).hexdigest()[:16]
        frames = [k for k in back if np.ndim(back[k]) == 3][0]
        print(f"  {name:20s} {info['bytes']:8d} bytes: {frames} {back[frames].shape}, t = 0 … {back['t'][-1]:.4g}; read-back differs from the "
              f"fresh run by at most {worst:.1e} of each array's range (float16 frames); sha256 of the arrays {digest}")
        assert worst < 2e-3, f"{name}: the stored file does not reproduce the run"
    print(f"total size of the files written: {total / 1024:.0f} kB (limit 2048 kB)")
    assert total < 2 * 1024 * 1024
    if which == "all":
        kb = ch13.load_reference_run("kelvin_basin", args.dest)
        perim = 4.0 * (kb["x"][-1] - kb["x"][0] + kb["x"][1] - kb["x"][0])
        print(f"kelvin_basin: Lambda = {float(kb['Lambda']) / 1e3:.1f} km, c = {float(kb['c']):.3f} m/s, basin {perim / 4 / float(kb['Lambda']):.0f} Lambda "
              f"wide; run covers {kb['t'][-1] * float(kb['c']) / perim:.2f} circuits; peak height {np.abs(kb['eta'][0]).max():.3f} -> "
              f"{np.abs(kb['eta'][-1]).max():.3f} m")
        j, i = np.unravel_index(np.argmax(kb["eta"][2]), kb["eta"][2].shape)
        j2, i2 = np.unravel_index(np.argmax(kb["eta"][0]), kb["eta"][0].shape)
        print(f"  the pulse starts on the south wall at x = {kb['x'][i2] / 1e3:.0f} km and has moved to x = {kb['x'][i] / 1e3:.0f} km "
              f"(row {j}) by t = {kb['t'][2] / 3600:.1f} h: travel toward +x with the coast on its right (f > 0); c × t = "
              f"{float(kb['c']) * kb['t'][2] / 1e3:.0f} km")
        tf, tb = ch13.load_reference_run("turbulence_f", args.dest), ch13.load_reference_run("turbulence_beta", args.dest)
        print(f"turbulence: K_E {tf['K_E'][0]:.2f} -> {tf['K_E'][-1]:.2f} (beta = 0) and -> {tb['K_E'][-1]:.2f} (beta = {float(tb['beta'])}); zonal "
              f"energy fraction at the end {ch13.zonal_energy_fraction(tf['zeta'][-1], float(tf['L'])):.3f} against "
              f"{ch13.zonal_energy_fraction(tb['zeta'][-1], float(tb['L'])):.3f}")
        pv = ch13.load_reference_run("pv_particles", args.dest)
        print(f"pv_particles: largest drift of q along a path {100 * np.max(np.abs(pv['q'] / pv['q'][0] - 1)):.2f} %")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
