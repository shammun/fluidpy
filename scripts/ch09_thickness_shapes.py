"""§9.2 (C02): the three thicknesses (delta_99, displacement delta*, momentum theta) and the shape factor H = delta*/theta on five model
profiles (linear, sine, cubic, exponential, Blasius), each integrated numerically and compared with the exact closed forms.

Run: ``.venv/Scripts/python.exe scripts/ch09_thickness_shapes.py --no-show``   Figure -> outputs/ch09/c02_thickness_shapes.png.
"""
from __future__ import annotations

import numpy as np

from ch09_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch09_boundary_layers as ch09


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    y = np.linspace(0.0, 12.0, 4801)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    print(f"{'profile':<12}{'delta*/delta':>14}{'exact':>10}{'theta/delta':>14}{'exact':>10}{'H':>8}")
    for name in ("linear", "sine", "cubic", "exponential", "blasius"):
        p = ch09.profile_shape(name, y, 1.0)
        u = p["u_over_Ue"]
        ds, th = ch09.displacement_thickness(y, u, 1.0), ch09.momentum_thickness(y, u, 1.0)
        print(f"{name:<12}{ds:>14.5f}{p['delta_star']:>10.5f}{th:>14.5f}{p['theta']:>10.5f}{ds / th:>8.3f}")
        ax[0].plot(u, y, label=f"{name}  H = {ds / th:.2f}")
    ax[0].set_ylim(0, 4)
    ax[0].set_xlabel("u / U_e")
    ax[0].set_ylabel("y / delta")
    ax[0].legend(fontsize=8)
    ax[0].set_title("model profiles")
    u = ch09.profile_shape("blasius", y, 1.0)["u_over_Ue"]
    th = ch09.thicknesses(y, u, 1.0)
    ax[1].plot(u, y, color=COLORS["accent"])
    ax[1].fill_betweenx(y, u, 1.0, where=y < th["delta_star"] + 3, color=COLORS["orange"], alpha=0.2, label="deficit")
    for k, c in (("delta99", COLORS["rose"]), ("delta_star", COLORS["teal"]), ("theta", COLORS["blue"])):
        ax[1].axhline(th[k], color=c, lw=1.2, label=f"{k} = {th[k]:.3f}")
    ax[1].set_ylim(0, 6)
    ax[1].set_xlabel("u / U")
    ax[1].set_ylabel("eta = y sqrt(U / nu x)")
    ax[1].legend(fontsize=8)
    ax[1].set_title("Blasius profile with its three thicknesses")
    save(fig, out, "c02_thickness_shapes")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
