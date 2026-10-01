import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(r"C:\Users\sislam27\Work\Climate Dynamics PHD\Fluid Dynamics")
sys.path.insert(0, str(ROOT))
from fluidpy import ch11_instability as ch11

T = json.load(open(ROOT / "reference/ch11/explainer_tables.json"))["taylor"]
mus = [-0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0]
rs = [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 3, 5, 8, 12, 20, 35, 60]
out = []
for m in mus:
    c = ch11.taylor_critical(m)
    row = []
    for r in rs:
        s = ch11.taylor_growth_rate(c["k_c"], r * c["Ta_c"], m)
        row.append(s)
    print(m, c, [f"{s.real:.4g}{'' if abs(s.imag) < 1e-6 else f'+{s.imag:.3g}j'}" for s in row])
    out.append([float(f"{s.real:.4g}") for s in row])
print(json.dumps(out))
# uR profiles from the table's psi (row at z = lambda/4, sin kz = 1)
for key in ["0", "0.5", "1"]:
    e = T["eigenfunctions"][key]
    psi = np.array(e["psi"]); k = e["k_c"]
    uR = psi[10] * k
    print(key, k, e["Ta"], np.max(uR), json.dumps([float(f"{v:.4g}") for v in uR]))
    ef = ch11.taylor_eigenfunction(k, float(key), x=[0.25, 0.5], z=[0.0])
    print("  uR(0.25,0.5)", ef["u_R"][0], "uphi", ef["u_phi"][0], "max uphi", np.max(np.abs(ef["uphi_profile"])))
print(ch11.taylor_number(0.5, 0.0, 0.1, 0.105, 1e-6), ch11.taylor_number_narrow_inner(0.5, 0.1, 0.005, 1e-6))
print(ch11.taylor_critical(0.999), ch11.taylor_critical(1.0))
for m in (-0.5, 0, 0.5):
    print(m, ch11.taylor_critical_approx(m) / ch11.taylor_critical(m)["Ta_c"] - 1)
