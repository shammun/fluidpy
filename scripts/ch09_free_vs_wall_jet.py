"""§9.10 (A6): free jet beside wall jet — exponents of the centre-line speed, thickness and entrained mass flux against x, measured by a log-log fit of
the closed forms, and the entrainment velocity toward the jet.

Run: ``.venv/Scripts/python.exe scripts/ch09_free_vs_wall_jet.py --no-show``   Figure -> outputs/ch09/a6_free_vs_wall_jet.png.
"""
from __future__ import annotations

import numpy as np

from ch09_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch09_boundary_layers as ch09


def slope(x, y):
    return float(np.polyfit(np.log(x), np.log(y), 1)[0])


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    Jm, rho, nu, C, finf = 1.0, 1.0, 1e-4, 0.7, 1.3
    x = np.geomspace(0.2, 20, 40)
    rows = {
        "free jet  u0": (ch09.free_jet_centreline(x, Jm, rho, nu), -1 / 3),
        "free jet  delta": (ch09.free_jet_thickness(x, Jm, rho, nu), 2 / 3),
        "free jet  mdot": (ch09.free_jet_mass_flux(x, Jm, rho, nu), 1 / 3),
        "wall jet  u0": (C * x ** -0.5, -1 / 2),
        "wall jet  delta": (ch09.wall_jet(x, 0.0, C, finf, nu)["delta"], 3 / 4),
        "wall jet  mdot": (ch09.wall_jet_mass_flux(x, C, finf, rho, nu), 1 / 4),
    }
    print(f"{'quantity':<18}{'fitted exponent':>16}{'theory':>9}")
    for k, (v, th) in rows.items():
        print(f"{k:<18}{slope(x, v):>16.4f}{th:>9.4f}")
    Re = np.array([1e2, 1e3, 1e4])
    print("entrainment v/u0 at eta -> +inf:", ", ".join(f"{float(ch09.free_jet_entrainment_velocity(r)):.4f}" for r in Re), "(Re_x = 1e2, 1e3, 1e4)")
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for k, (v, th) in rows.items():
        if "u0" in k or "delta" in k:
            ax[0 if "u0" in k else 1].loglog(x, v, label=f"{k}  ~ x^{th:.3g}")
    ax[0].set_xlabel("x [m]")
    ax[0].set_ylabel("u0 [m/s]")
    ax[0].legend(fontsize=8)
    ax[1].set_xlabel("x [m]")
    ax[1].set_ylabel("delta [m]")
    ax[1].legend(fontsize=8)
    ax[0].set_title("centre-line speed")
    ax[1].set_title("thickness")
    save(fig, out, "a6_free_vs_wall_jet")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
