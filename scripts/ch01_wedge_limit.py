"""Chapter 1, §1.7, Eq. (1.6): shrink the triangular element of Fig. 1.5 and watch p2 − p1 → 0 (pressure is a scalar).

Run: ``.venv/Scripts/python.exe scripts/ch01_wedge_limit.py --no-show``
Figure → outputs/ch01/fig1_5_wedge_limit.png.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch01"), help="figure folder")
    ap.add_argument("--no-show", action="store_true", help="do not open a window")
    args = ap.parse_args()
    t_start = time.perf_counter()
    import matplotlib

    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from ch01_drawings import draw_wedge
    from fluidpy import ch01_introduction as ch01
    from fluidpy.core.style import COLORS, use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    rho, p1 = 1000.0, ch01.P_ATM
    dz = np.logspace(-6, 0, 61)  # wedge height [m]
    fig, axs = plt.subplots(1, 2, figsize=(10, 4), gridspec_kw={"width_ratios": [1, 1.4]})
    draw_wedge(axs[0], 1.0, np.radians(35.0))
    axs[0].set_title("forces on the wedge (weight acts at its centre)")
    for th_deg, c, ls in ((20, COLORS["accent"], "-"), (45, COLORS["teal"], "--"), (70, COLORS["orange"], ":")):
        d = ch01.wedge_pressure_difference(rho, dz, np.radians(th_deg))
        axs[1].loglog(dz, d["p2_minus_p1"], color=c, ls=ls, lw=2, label=f"p2 − p1, θ = {th_deg}°")
    ff = ch01.wedge_face_forces(p1, rho, dz, np.radians(35.0))
    axs[1].loglog(dz, ff["weight_to_face_force"], color=COLORS["rose"], label="weight / bottom-face force")
    axs[1].set_xlabel("wedge size dz [m]")
    axs[1].set_ylabel("[Pa]  or  [-]")
    axs[1].set_title("Both vanish ∝ dz: in the limit p1 = p2 = p3 (Eq. 1.6)")
    axs[1].legend(fontsize=8)
    fig.savefig(out / "fig1_5_wedge_limit.png", bbox_inches="tight")

    slope = np.polyfit(np.log(dz), np.log(ch01.wedge_pressure_difference(rho, dz, 0.6)["p2_minus_p1"]), 1)[0]
    print(f"water wedge: p2 − p1 = {ch01.wedge_pressure_difference(rho, 1e-3, 0.6)['p2_minus_p1']:.4f} Pa at dz = 1 mm; "
          f"log-log slope vs dz = {slope:.4f}")
    print(f"balance residuals at dz = 1 mm: res_x = {ff['res_x'][30]:.2e}, res_z = {ff['res_z'][30]:.2e} N/m; "
          f"weight/face force = {ff['weight_to_face_force'][30]:.2e}")
    print(f"saved fig1_5_wedge_limit.png in {out}  ({time.perf_counter() - t_start:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
