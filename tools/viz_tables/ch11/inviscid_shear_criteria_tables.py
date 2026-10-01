"""Generator of the tables embedded in viz/ch11/inviscid_shear_criteria.html (run from anywhere; writes only next to itself).

1. RAYLEIGH   = reference/ch11/rayleigh_spectra.json (ch11.rayleigh_spectrum_table defaults), read, not recomputed.
2. JET_ODD    = ch11.rayleigh_spectrum_table(["jet_varicose"], cache=False, write=False): phi odd, k = 0.1 ... 2.0.
3. SINB       = U = sin y between walls at +-b, b = 1.6 ... 3.0 step 0.1, k = 0.05 ... 1.0 step 0.05, by
                core.stability.rayleigh_eigs_contour (N = 120, the settings of rayleigh_spectrum_table).
All values rounded to 4 significant figures, c_i = 0 where no mode with c_i > 1e-4 exists.
Also prints off-node check values used to choose honest interpolation tolerances.
"""
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(r"C:\Users\sislam27\Work\Climate Dynamics PHD\Fluid Dynamics")
sys.path.insert(0, str(ROOT))
from fluidpy import ch11_instability as ch11  # noqa: E402

here = Path(__file__).parent
g4 = lambda v: float(f"{v:.4g}")  # noqa: E731
t0 = time.time()
out = {}

ref = json.loads((ROOT / "reference" / "ch11" / "rayleigh_spectra.json").read_text(encoding="utf-8"))["spectra"]
for nm in ("shear_layer", "jet"):
    out[nm] = dict(cr=[0.0 if v is None or abs(v) < 1e-9 else v for v in ref[nm]["c_r"]], ci=ref[nm]["c_i"])

t = ch11.rayleigh_spectrum_table(["jet_varicose", "jet_sinuous"], cache=False, write=False)
out["jet_odd"] = dict(cr=[0.0 if not np.isfinite(v) else g4(v) for v in t["jet_varicose"]["c_r"]],
                      ci=[g4(v) for v in t["jet_varicose"]["c_i"]])
sinu = dict(cr=[0.0 if not np.isfinite(v) else g4(v) for v in t["jet_sinuous"]["c_r"]], ci=[g4(v) for v in t["jet_sinuous"]["c_i"]])
print("jet_sinuous == cache 'jet':", sinu == out["jet"])
print("jet_odd", out["jet_odd"])

ks = np.round(np.arange(0.05, 1.0001, 0.05), 10)
sinb = {}
for b in np.round(np.arange(1.6, 3.0001, 0.1), 10):
    pr = ch11.parallel_profile("sin", b=float(b))
    ci = []
    for k in ks:
        c = ch11.rayleigh_eigs_contour(float(k), pr["U"], pr["Up"], pr["Upp"], domain=pr["domain"], N=120, bc=pr["bc"])
        ci.append(g4(c[0].imag) if len(c) else 0.0)
        if len(c) and abs(c[0].real) > 1e-6:
            print("  note: c_r != 0", b, k, c[0])
    while ci and ci[-1] == 0.0:
        ci.pop()
    sinb[f"{b:.1f}"] = ci
    kn = math.sqrt(1 - (math.pi / (2 * b)) ** 2)
    print(f"sin b={b:.1f} kn={kn:.4f} last node k={ks[len(ci) - 1]:.2f}", ci, flush=True)
out["sinb"] = sinb
(here / "tables_v2.json").write_text(json.dumps(out), encoding="utf-8")

# the JS constants, ready to paste
L = []
L.append("const RAYLEIGH = {")
L.append(f"  shear_layer: {{ cr: Z20, ci: {json.dumps(out['shear_layer']['ci'])} }},")
L.append(f"  jet: {{ cr: {json.dumps(out['jet']['cr'])},\n    ci: {json.dumps(out['jet']['ci'])} }}\n}};")
L.append(f"const JET_ODD = {{ cr: {json.dumps(out['jet_odd']['cr'])},\n  ci: {json.dumps(out['jet_odd']['ci'])} }};")
L.append("const SINB = {")
L.append(",\n".join(f"  '{b}': {json.dumps(v)}" for b, v in sinb.items()))
L.append("};")
(here / "tables_v2.js").write_text("\n".join(L) + "\n", encoding="utf-8")

# off-node check values (for the interpolation tolerances of the selftest rows)
jet = ch11.parallel_profile("jet")
tanh = ch11.parallel_profile("tanh")


def lead(pr, k, **kw):
    ym = max(float(pr["y_max"]), ch11.decay_box(k)) if pr["bc"] == "decay" else None
    c = ch11.rayleigh_eigs_contour(k, pr["U"], pr["Up"], pr["Upp"], domain=pr["domain"], N=120, bc=pr["bc"], y_max=ym,
                                   map_scale=pr.get("map_scale"), **kw)
    return complex(c[0]) if len(c) else None


chk = {}
for k in (1.75, 1.85, 1.95, 1.98):
    chk[f"jet_even_{k}"] = lead(jet, k, parity="even")
for k in (0.85, 0.95, 0.98):
    chk[f"jet_odd_{k}"] = lead(jet, k, parity="odd")
for k in (0.445, 0.95, 0.98):
    chk[f"tanh_{k}"] = lead(tanh, k)
for b, k in ((1.6, 0.12), (1.6, 0.17), (1.7, 0.224), (1.7, 0.33), (1.7, 0.37), (2.0, 0.58), (2.0, 0.61), (3.0, 0.83)):
    chk[f"sin_{b}_{k}"] = lead(ch11.parallel_profile("sin", b=b), k)
for k_, v in chk.items():
    print("check", k_, v)
for b in (1.6, 1.7, 2.0, 3.0):
    print("sin_profile_max_growth", b, ch11.sin_profile_max_growth(b))
print("tanh_max_growth", ch11.tanh_max_growth())
(here / "checks_v2.json").write_text(json.dumps({k_: (None if v is None else [v.real, v.imag]) for k_, v in chk.items()}), encoding="utf-8")
print("total %.0f s" % (time.time() - t0))
