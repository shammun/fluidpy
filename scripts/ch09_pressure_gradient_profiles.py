"""§9.7 (C08, N71-N73, Fig. 9.12 idea): boundary-layer profiles under a favourable (n = 1), zero (n = 0) and adverse (n = -0.05) pressure gradient;
the sign of d2u/dy2 at the wall follows dp/dx (Eq. 9.9 at y = 0) and an adverse gradient creates an inflection point.

Run: ``.venv/Scripts/python.exe scripts/ch09_pressure_gradient_profiles.py --no-show``   Figure -> outputs/ch09/c08_pressure_gradient_profiles.png.
"""
from __future__ import annotations

import numpy as np

from ch09_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch09_boundary_layers as ch09


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for n, col, lab in ((1.0, COLORS["teal"], "n = 1 (favourable)"), (0.0, COLORS["accent"], "n = 0 (no gradient)"),
                        (-0.05, COLORS["rose"], "n = -0.05 (adverse)")):
        d = ch09.falkner_skan(n)
        s = ch09.falkner_skan_state(n)
        y = d["eta"] * np.sqrt(2 / (n + 1))
        ax[0].plot(d["fp"], y, color=col, label=lab)
        ax[1].plot(d["fppp"], y, color=col, label=lab)
        infl = ch09.profile_inflection(d["eta"], d["fp"])
        print(f"  {lab:<22} f''(0) = {s['fpp0']:.4f}; f'''(0) = {s['fppp0']:+.4f} (= -n); wall curvature sign = "
              f"{np.sign(float(ch09.wall_curvature(-n, 1.0))):+.0f} x dp/dx; inflection at eta = {infl}")
    ax[0].set_ylim(0, 6)
    ax[0].set_xlabel("u / U_e")
    ax[0].set_ylabel("y / delta_(2/(n+1))")
    ax[0].legend(fontsize=8)
    ax[0].set_title("profiles: adverse gradient -> S-shape")
    ax[1].set_ylim(0, 6)
    ax[1].axvline(0, color=COLORS["muted"], lw=0.8)
    ax[1].set_xlabel("d2u/dy2 (scaled) = f'''")
    ax[1].legend(fontsize=8)
    ax[1].set_title("curvature: sign at the wall = sign of dp/dx")
    save(fig, out, "c08_pressure_gradient_profiles")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
