"""§9.1 (C01): the boundary-layer picture — a thin layer on a flat plate, thickness growing like sqrt(x) (Blasius), profiles at several
stations, and the scaling delta/L ~ Re^(-1/2).  Our figure; every number from ``fluidpy.core.boundary_layer``.

Run: ``.venv/Scripts/python.exe scripts/ch09_boundary_layer_picture.py --no-show``   Figure -> outputs/ch09/c01_boundary_layer_picture.png.
"""
from __future__ import annotations

import numpy as np

from ch09_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch09_boundary_layers as ch09


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    U, nu, L = 1.0, 1.5e-5, 1.0
    sc = ch09.boundary_layer_scales(U, L, nu)
    print(f"air, U = {U} m/s, L = {L} m: Re = {sc['Re']:.3g}, delta/L ~ Re^-1/2 = {sc['delta_over_L']:.4f}, v scale = {sc['v_scale']:.4f} m/s")
    x = np.linspace(1e-3, L, 400)
    d99 = ch09.blasius_delta99(x, U, nu)
    print(f"Blasius delta_99 at x = {L} m: {float(ch09.blasius_delta99(L, U, nu)) * 100:.3f} cm  (eta_99 = {ch09.blasius_constants()['eta99']:.4f})")
    fig, ax = plt.subplots(figsize=(9, 3.6))
    ax.fill_between(x, 0, d99 * 1e3, color=COLORS["accent"], alpha=0.15, label="boundary layer (delta_99)")
    ax.plot(x, d99 * 1e3, color=COLORS["accent"], lw=2)
    for xs in (0.1, 0.3, 0.6, 1.0):
        y = np.linspace(0, 1.4 * float(ch09.blasius_delta99(xs, U, nu)), 60)
        u = ch09.blasius_fields(xs, y, U, nu)["u"]
        ax.plot(xs + 0.12 * u / U, y * 1e3, color=COLORS["teal"], lw=1.6)
        ax.plot([xs, xs], [0, y[-1] * 1e3], color=COLORS["muted"], lw=0.6)
    ax.set_xlabel("x  [m]")
    ax.set_ylabel("y  [mm]")
    ax.set_title("laminar layer on a flat plate: delta grows like sqrt(x); profiles drawn at 4 stations (u/U scaled to 0.12 m)")
    ax.legend(loc="upper left")
    save(fig, out, "c01_boundary_layer_picture")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
