"""§12.8 (C09, N123-N125): the seven self-similar free shear flows — exact exponents from each flow's invariant and growth
law (our derivation of the exponents of Table 12.1), the centreline laws with ILLUSTRATIVE constants (the book's constants
are private and no public table has been confirmed), invariants checked numerically, and the method of Example 12.2 with our
own gas and nozzle.

Run: ``.venv/Scripts/python.exe scripts/ch12_free_shear_table.py --no-show``
Figure -> outputs/ch12/c09_free_shear_flows.png.
Our dilution case: propane (M = 44.10 kg/kmol) from a 5 mm nozzle at 30 m/s into air (M = 28.96), 5 mol O2 per mol fuel.
"""
from __future__ import annotations

import numpy as np

from ch12_common import COLORS, ILLUSTRATIVE_JET, finish, parse_args, save, setup

from fluidpy import ch12_turbulence as ch12


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    print(f"FREE_SHEAR_CONSTANTS has {len(ch12.FREE_SHEAR_CONSTANTS)} entries (no default set); demos use ILLUSTRATIVE {ILLUSTRATIVE_JET}")
    print(f"{'flow':13s} {'width':>6s} {'velocity':>9s} {'scalar':>7s} {'Re':>5s}   invariant")
    for flow in ch12.FREE_SHEAR_FLOWS:
        e = ch12.free_shear_exponents(flow, return_equations=True)
        print(f"{flow:13s} {str(e['width']):>6s} {str(e['velocity']):>9s} {str(e['scalar']):>7s} {str(e['reynolds']):>5s}   {e['invariant']}")
        for q in e["equations"]:
            print(f"{'':15s}{q}")
    C = ILLUSTRATIVE_JET
    src = dict(d=0.01, U0=10.0, rho_s=1.0, rho=1.2)
    r = np.linspace(0.0, 4.0, 40001)
    print("invariants with the illustrative constants (should not depend on x):")
    for x in (0.5, 1.0, 2.0):
        pj = ch12.free_shear_flow("plane_jet", x, r, constants=C, **src)
        rj = ch12.free_shear_flow("round_jet", x, r, constants=C, **src)
        pw = ch12.free_shear_flow("plane_wake", x, r, constants=C, theta=0.01, U_inf=5.0)
        rw = ch12.free_shear_flow("round_wake", x, r, constants=C, theta=0.01, U_inf=5.0)
        print(f"  x = {x}: plane jet 2 int U^2 dy = {2 * np.trapezoid(pj['U'] ** 2, r):.5f}; round jet int U^2 2 pi r dr = "
              f"{np.trapezoid(rj['U'] ** 2 * 2 * np.pi * r, r):.6f}; plane wake 2 int dU dy = {2 * np.trapezoid(5.0 - pw['U'], r):.5f}; "
              f"round wake int dU 2 pi r dr = {np.trapezoid((5.0 - rw['U']) * 2 * np.pi * r, r):.6f}")
    sl = ch12.free_shear_flow("shear_layer", 1.0, np.array([-0.2, 0.0, 0.2]), constants=dict(dxi80_coeff=0.1), U1=10.0, U2=5.0)
    print(f"shear layer U1 = 10, U2 = 5 m/s at x = 1 m: U(y = -0.2, 0, 0.2) = {np.round(sl['U'], 3)}, 80 % span = {sl['dxi80']:.4f} x")
    for f in ch12.FREE_SHEAR_FLOWS:
        print(f"  local Reynolds number of the {f:12s} ~ x^{ch12.local_reynolds_number_exponent(f)}")

    # Example 12.2 method, our inputs, illustrative constant
    st = ch12.stoichiometric_mass_fraction(44.10, 28.96, 5.0, 0.21)
    rho_ratio = 44.10 / 28.96
    x_st = ch12.round_jet_distance_for_mass_fraction(st["Y_fuel"], 0.005, rho_ratio, 1.0, 1.0, C_Y=5.0)
    cl = ch12.free_shear_centerline("round_jet", x_st, constants=dict(C_U=6.0), d=0.005, U0=30.0, rho_s=rho_ratio, rho=1.0)
    print(f"dilution of a propane jet (illustrative C_Y = 5, C_U = 6): stoichiometric Y = {st['Y_fuel']:.4f} "
          f"(fuel volume fraction {st['v_fuel']:.4f}); reached on the axis at x = {x_st:.3f} m = {x_st / 0.005:.0f} d, "
          f"where U_CL = {cl['U_CL']:.2f} m/s")

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    x = np.geomspace(0.1, 10.0, 50)
    for flow, c in zip(("plane_jet", "round_jet", "plane_wake", "round_wake", "round_plume"),
                       (COLORS["accent"], COLORS["teal"], COLORS["orange"], COLORS["rose"], COLORS["blue"])):
        e = ch12.free_shear_exponents(flow)
        ax[0].loglog(x, x ** float(e["velocity"]), color=c, label=f"{flow.replace('_', ' ')}: x^{e['velocity']}")
        ax[1].loglog(x, x ** float(e["width"]), color=c, label=f"{flow.replace('_', ' ')}: x^{e['width']}")
    ax[0].set(xlabel="x / x_ref", ylabel="centreline velocity (or deficit) / value at x_ref", title="decay laws of Table 12.1 (exact exponents)")
    ax[1].set(xlabel="x / x_ref", ylabel="width / value at x_ref", title="growth laws")
    ax[0].legend(fontsize=8)
    ax[1].legend(fontsize=8)
    save(fig, out, "c09_free_shear_flows")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
