"""§9.11 (C14): secondary flow in a stirred teacup.  Away from the wall the swirl u_e(R) is in radial balance with the pressure gradient
dp/dR = rho u_e^2 / R; in the thin bottom layer dp/dz ~ 0, so the layer feels the same pressure gradient while its slowed swirl u < u_e needs less: the net
inward force rho (u_e^2 - u^2)/R drives fluid toward the axis along the bottom (tea leaves pile up at the centre).  ILLUSTRATIVE profile (not a solution);
the Ekman-layer version of this is Ch. 13.

Run: ``.venv/Scripts/python.exe scripts/ch09_teacup.py --no-show``   Figure -> outputs/ch09/c14_teacup.png.
"""
from __future__ import annotations

import numpy as np

from ch09_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch09_boundary_layers as ch09


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    rho, R, ue, delta = 1000.0, 0.04, 0.15, 2e-3
    z = np.linspace(0, 6e-3, 200)
    u = ch09.secondary_flow_layer_profile(z, delta, ue, "sine")
    F = ch09.secondary_flow_radial_force(ue, u, R, rho)
    print(f"cup R = {R} m, swirl u_e = {ue} m/s: inward force at the floor rho u_e^2/R = {float(ch09.secondary_flow_radial_force(ue, 0.0, R, rho)):.1f} N/m^3 "
          f"(the full centripetal demand is missing there); zero above the layer ({float(F[-1]):.1e})")
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    ax[0].plot(u, z * 1e3, color=COLORS["accent"])
    ax[0].set_xlabel("swirl u_phi [m/s]")
    ax[0].set_ylabel("height above the floor z [mm]")
    ax[0].set_title("bottom layer slows the swirl (illustrative)")
    ax[1].plot(F, z * 1e3, color=COLORS["rose"])
    ax[1].set_xlabel("net inward force rho (u_e^2 - u^2)/R  [N/m^3]")
    ax[1].set_ylabel("z [mm]")
    ax[1].set_title("unbalanced pressure gradient -> inward flow along the floor")
    save(fig, out, "c14_teacup")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
