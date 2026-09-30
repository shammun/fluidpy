"""§9.8 (C10, A5): the Karman vortex street as two rows of point vortices.  The staggered street with lateral/longitudinal spacing b/a = arccosh(sqrt2)/pi = 0.2805
is the neutrally stable one: the growth rate of the linearised double row vanishes there (to the eigenvalue noise ~ 1e-8) and grows linearly on either side;
the symmetric (non-staggered) rows are unstable at every spacing.  Also the self-induced speed (Gamma/2a) tanh(pi b/a).

Run: ``.venv/Scripts/python.exe scripts/ch09_karman_street.py --no-show``   Figure -> outputs/ch09/c10_karman_street.png.
"""
from __future__ import annotations

import numpy as np

from ch09_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch09_boundary_layers as ch09


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    r0 = ch09.karman_street_ratio()
    print(f"b/a = arccosh(sqrt2)/pi = {r0:.5f};  cosh^2(pi b/a) = {np.cosh(np.pi * r0) ** 2:.12f}")
    ba = np.linspace(0.05, 0.6, 111)
    g = np.array([ch09.karman_street_growth(b, k=np.pi) for b in ba])  # alternate-vortex mode k = pi/a (exact lattice sums)
    g_sym = np.array([ch09.karman_street_growth(b, offset=0.0, k=np.pi) for b in ba])
    gc = np.array([float(ch09.karman_street_growth_closed(b)) for b in ba])
    print(f"staggered street: growth at b/a = 0.2805: {ch09.karman_street_growth(r0):.1e}; at 0.2: {ch09.karman_street_growth(0.2):.4f}; "
          f"at 0.4: {ch09.karman_street_growth(0.4):.4f}  [units Gamma/a^2, max over 61 wavenumbers]")
    print(f"closed form (pi/2)|1/2 - sech^2(pi b/a)| agrees with the matrix: max diff {np.max(np.abs(g - gc)):.1e}")
    print(f"non-staggered rows: {ch09.karman_street_growth(0.3, offset=0.0):.4f} (= pi/4 = {np.pi / 4:.4f}) for every spacing")
    print(f"street speed for Gamma = 1 m^2/s, a = 1 m: {float(ch09.karman_street_velocity(1.0, r0, 1.0)):.4f} m/s")
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].plot(ba, g, color=COLORS["accent"], label="staggered")
    ax[0].plot(ba, g_sym, color=COLORS["orange"], ls="--", label="non-staggered (symmetric)")
    ax[0].axvline(r0, color=COLORS["rose"], ls=":", label="b/a = 0.2805")
    ax[0].set_xlabel("b/a")
    ax[0].set_ylabel("max growth rate  [Gamma/a^2]")
    ax[0].legend(fontsize=8)
    ax[0].set_title("linear stability of a double row of point vortices")
    for b, col in ((r0, COLORS["teal"]), (0.4, COLORS["rose"])):
        p = ch09.karman_street_positions(4.0, b, eps=0.02, mode="unstable", n_cells=8)
        ax[1].plot(p["zA"].real, p["zA"].imag - (0 if b == r0 else 1.2), "o", color=col, ms=5)
        ax[1].plot(p["zB"].real, p["zB"].imag - (0 if b == r0 else 1.2), "s", color=col, ms=5, label=f"b/a = {b:.3f}, after t = 4 a^2/Gamma")
    ax[1].set_aspect("equal")
    ax[1].set_xlabel("x / a")
    ax[1].set_ylabel("y / a")
    ax[1].legend(fontsize=7, loc="upper right")
    ax[1].set_title("unstable mode after t = 4 a^2/Gamma (0.2805 street; b/a = 0.4 street shifted down)")
    save(fig, out, "c10_karman_street")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
