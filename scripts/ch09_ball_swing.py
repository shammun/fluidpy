"""§9.9 (N84): cricket-ball swing — the seam trips the boundary layer on one side only, the laminar side separates early (~85 deg) and the turbulent side later
(~120 deg), so the wake is pushed sideways and a constant lateral force bends the path into a parabola y = a t^2 / 2.  QUALITATIVE picture from our
separated-pressure model; the deflection is our one-line kinematics with caller-chosen parameters (an assumed side force of 10 % of the weight).

Run: ``.venv/Scripts/python.exe scripts/ch09_ball_swing.py --no-show``   Figure -> outputs/ch09/n84_ball_swing.png.
"""
from __future__ import annotations

import numpy as np

from ch09_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch09_boundary_layers as ch09


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    lam, turb = 85.0, 120.0
    cd_l, cd_t = ch09.separated_pressure_drag(lam, -1.0), ch09.separated_pressure_drag(turb, -0.6)
    print(f"model pressure drag: laminar side (sep {lam:.0f} deg) {cd_l:.3f}, turbulent side (sep {turb:.0f} deg) {cd_t:.3f}")
    d, U, FW = 20.0, 35.0, 0.10
    y = ch09.ball_swing_deflection(FW, d, U)
    print(f"assumed F/W = {FW}: sideways deflection over {d} m at {U} m/s = {y:.3f} m (quadratic in distance: x2 distance -> x{y and ch09.ball_swing_deflection(FW, 2 * d, U) / y:.0f})")
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    phi = np.linspace(0, 180, 361)
    ax[0].plot(phi, ch09.separated_cp(phi, lam, -1.0), color=COLORS["orange"], label="laminar side (seam-free), sep 85 deg")
    ax[0].plot(phi, ch09.separated_cp(phi, turb, -0.6), color=COLORS["teal"], label="turbulent side (seam trips), sep 120 deg")
    ax[0].set_xlabel("phi from the front [deg]")
    ax[0].set_ylabel("C_p (model)")
    ax[0].legend(fontsize=8)
    ax[0].set_title("asymmetric separation")
    dist = np.linspace(0, 20, 100)
    ax[1].plot(dist, ch09.ball_swing_deflection(FW, dist, U), color=COLORS["accent"])
    ax[1].set_xlabel("distance travelled [m]")
    ax[1].set_ylabel("sideways deflection y = a t^2/2 [m]")
    ax[1].set_title("parabolic path (F/W = 0.10 assumed)")
    save(fig, out, "n84_ball_swing")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
