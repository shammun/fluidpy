"""§9.2 (N19): streamlines of the Blasius layer are pushed away from the plate; the outer streamlines are displaced by delta* and the
wall-normal velocity approaches v_inf = 0.8604 U / sqrt(Re_x) (Fig. 9.4/9.6 with our solver).

Run: ``.venv/Scripts/python.exe scripts/ch09_displacement_streamlines.py --no-show``   Figure -> outputs/ch09/n19_displacement_streamlines.png.
"""
from __future__ import annotations

import numpy as np

from ch09_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch09_boundary_layers as ch09


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    U, nu = 1.0, 1e-5
    x = np.linspace(0.02, 1.0, 200)
    y = np.linspace(0, 0.02, 300)
    X, Y = np.meshgrid(x, y)
    fl = ch09.blasius_fields(X, Y, U, nu)
    c = ch09.blasius_constants()
    ds = c["delta_star"] * np.sqrt(nu * x / U)
    print(f"delta*/sqrt(nu x/U) = {c['delta_star']:.4f};  v_inf sqrt(Re_x)/U = {c['v_inf']:.4f}")
    v_edge = fl["v"][-1, :]
    print(f"v at y = 2 cm, x = 1 m: {v_edge[-1]:.3e} m/s  vs 0.8604 sqrt(nu U/x)/...: {c['v_inf'] * np.sqrt(nu * U / x[-1]):.3e}")
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.contour(X, Y * 1e3, fl["psi"], levels=np.linspace(0, fl["psi"].max(), 14), colors=COLORS["accent"], linewidths=1)
    ax.plot(x, ds * 1e3, color=COLORS["orange"], lw=2, label="delta* (streamlines above are displaced by this much)")
    ax.plot(x, ch09.blasius_delta99(x, U, nu) * 1e3, color=COLORS["rose"], lw=1.5, ls="--", label="delta_99")
    ax.set_xlabel("x  [m]")
    ax.set_ylabel("y  [mm]")
    ax.legend(loc="upper left", fontsize=8)
    ax.set_title("Blasius streamlines: pushed outward by the slowed fluid")
    save(fig, out, "n19_displacement_streamlines")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
