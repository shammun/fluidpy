"""§9.3 (N41): Blasius (spatial) vs the temporally developing boundary layer of ch08 Example 8.5:
C_f sqrt(Re_x) = 0.664 (local Blasius) vs 2/sqrt(pi) = 1.128 (temporal layer with U t -> x) vs 1.328 (plate-averaged C_D sqrt(Re_L)).

Run: ``.venv/Scripts/python.exe scripts/ch09_blasius_vs_temporal.py --no-show``   Figure -> outputs/ch09/n41_blasius_vs_temporal.png.
"""
from __future__ import annotations

import numpy as np

from ch09_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch09_boundary_layers as ch09
from fluidpy.core import laminar as LAM


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    U, nu, rho = 1.0, 1e-5, 1.0
    x = np.geomspace(1e-2, 1.0, 100)
    Rex = U * x / nu
    cf_b = ch09.blasius_skin_friction(Rex)
    cf_t = np.array([LAM.temporal_bl_wall_stress(xx / U, U, nu, rho)["Cf"] for xx in x])
    c = ch09.blasius_constants()
    print(f"Cf sqrt(Re_x): Blasius {c['cf_coeff']:.4f};  temporal {2 / np.sqrt(np.pi):.4f};  ratio {2 / np.sqrt(np.pi) / c['cf_coeff']:.3f}")
    print(f"plate-averaged C_D sqrt(Re_L) = {c['cd_coeff']:.4f} = 2 x local (mean of x^-1/2)")
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.loglog(Rex, cf_b, color=COLORS["accent"], label="Blasius 0.664/sqrt(Re_x)")
    ax.loglog(Rex, cf_t, color=COLORS["orange"], label="temporal layer (t = x/U): 1.128/sqrt(Re_x)")
    ax.loglog(Rex, ch09.blasius_drag_coefficient(Rex), color=COLORS["teal"], ls="--", label="plate mean C_D = 1.328/sqrt(Re_L)")
    ax.set_xlabel("Re_x = U x / nu")
    ax.set_ylabel("skin friction coefficient")
    ax.legend(fontsize=8)
    ax.set_title("why the Galilean map t = x/U overstates the wall stress")
    save(fig, out, "n41_blasius_vs_temporal")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
