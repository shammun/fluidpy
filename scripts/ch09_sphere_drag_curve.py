"""§9.9 (N83, Fig. 9.22 idea): sphere drag coefficient vs Re = U d / nu from Morrison's correlation (ch04, benchmark) with the Stokes and Oseen
asymptotes (ch08) — the curve dips at the drag crisis Re ~ 3-5 e5.

Run: ``.venv/Scripts/python.exe scripts/ch09_sphere_drag_curve.py --no-show``   Figure -> outputs/ch09/n83_sphere_drag_curve.png.
"""
from __future__ import annotations

import numpy as np

from ch09_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch09_boundary_layers as ch09
from fluidpy.core import creeping as CRP
from fluidpy.core import similarity as SIM


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    Re = np.logspace(-2, 6.5, 400)
    cd = SIM.sphere_drag_coefficient(Re, "morrison")
    i = int(np.argmin(np.where(Re > 1e5, cd, np.inf)))
    print(f"Morrison curve: minimum after 1e5 at Re = {Re[i]:.3g}, C_D = {cd[i]:.3f}; C_D(1e4) = {float(SIM.sphere_drag_coefficient(1e4)):.3f}")
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    ax.loglog(Re, cd, color=COLORS["accent"], lw=2, label="Morrison correlation (data fit)")
    lo = Re < 5
    ax.loglog(Re[lo], CRP.stokes_drag_coefficient(Re[lo]), color=COLORS["teal"], ls="--", label="Stokes 24/Re")
    ax.loglog(Re[Re < 10], CRP.oseen_drag_coefficient(Re[Re < 10]), color=COLORS["orange"], ls="--", label="Oseen")
    ax.axvspan(3e5, 5e5, color=COLORS["rose"], alpha=0.15, label="drag crisis")
    ax.set_xlabel("Re = U d / nu")
    ax.set_ylabel("C_D")
    ax.legend(fontsize=8)
    ax.set_title("sphere drag coefficient")
    print("regime at Re = 5e5:", ch09.sphere_flow_regime(6e5)["label"])
    save(fig, out, "n83_sphere_drag_curve")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
