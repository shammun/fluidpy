"""§9.10 (C12, C13): the plane laminar free jet (Bickley: sech^2 profile, u0 ~ x^(-1/3), delta ~ x^(2/3), mass flux ~ x^(1/3)) and the wall jet
(Glauert: implicit solution (9.83), u0 ~ x^(-1/2), delta ~ x^(3/4), mass flux ~ x^(1/4)); solve_bvp / solve_ivp cross-checks and the conserved quantities.

Run: ``.venv/Scripts/python.exe scripts/ch09_jets.py --no-show``   Figure -> outputs/ch09/c12_jets.png.
"""
from __future__ import annotations

import numpy as np

from ch09_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch09_boundary_layers as ch09


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    k = ch09.free_jet_constants()
    print(f"free jet: C = {k['C']:.5f} (= 4 sqrt6/3), mass-flux coefficient 36^(1/3) = {k['mdot_coeff']:.4f}, Bickley 0.4543/0.2752 -> {k['u0_coeff']:.4f}/{k['xi_coeff']:.4f}")
    r = ch09.free_jet_ode_solve()
    print(f"  solve_bvp vs sqrt6 tanh(eta/sqrt6): max error {r['max_err']:.1e}")
    print(f"  half-width h99 coefficient: correct {k['h99_coeff']:.4f} (sech^2 = 0.01), printed {k['h99_printed']}")
    Jm, rho, nu = 1.0, 1.0, 1e-4
    y = np.linspace(-1, 1, 4001)
    for x in (0.5, 1.0, 4.0):
        f = ch09.free_jet(x, y, Jm, rho, nu)
        print(f"  x = {x}: momentum flux {ch09.jet_momentum_flux(y, f['u'], rho):.6f} (J = {Jm}), mass flux {rho * np.trapezoid(f['u'], y):.5f} "
              f"(formula {float(ch09.free_jet_mass_flux(x, Jm, rho, nu)):.5f})")
    w = ch09.wall_jet_ode_solve(1.0)
    print(f"wall jet: f_inf = {w['f_inf']:.5f}, f''(0)/f_inf^3 = {w['fpp0_over_f_inf_cubed']:.6f} (1/72 = {1 / 72:.6f}); IVP vs (9.83): {w['err_vs_9_83']:.1e}")
    wp = ch09.wall_jet_ode_solve(1.0, printed=True)
    print(f"  printed ODE (coefficient 1): error vs (9.83) = {wp['err_vs_9_83']:.2f}  -> the printed equation is inconsistent with the book's own integral")
    C, finf = 0.7, 1.3
    yy = np.linspace(0, 3, 6001)
    inv = [ch09.wall_jet_invariant(yy, ch09.wall_jet(x, yy, C, finf, nu)["u"]) for x in (0.5, 1.0, 4.0)]
    print(f"  invariant (9.80) at x = 0.5, 1, 4: {inv[0]:.6e}, {inv[1]:.6e}, {inv[2]:.6e}")
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    eta = np.linspace(-9, 9, 400)
    _pr = ch09.free_jet_profile(eta)
    f, fp = _pr['f'], _pr['fp']
    ax[0].plot(eta, fp, color=COLORS["accent"], label="f' = sech^2(eta/sqrt6)")
    ax[0].plot(eta, f / np.sqrt(6), color=COLORS["teal"], label="f / sqrt6 = tanh")
    ax[0].set_xlabel("eta")
    ax[0].legend(fontsize=8)
    ax[0].set_title("free jet similarity profile")
    e = np.linspace(0, 8, 300)
    _pw = ch09.wall_jet_profile(e, 1.0)
    fw, fpw = _pw['f'], _pw['fp']
    ax[1].plot(e, fpw, color=COLORS["accent"], label="f' (f_inf = 1)")
    ax[1].plot(e, fw, color=COLORS["teal"], label="f")
    ax[1].set_xlabel("eta")
    ax[1].legend(fontsize=8)
    ax[1].set_title("wall jet, from the implicit solution (9.83)")
    x = np.geomspace(0.1, 10, 50)
    ax[2].loglog(x, ch09.free_jet_thickness(x, Jm, rho, nu), color=COLORS["accent"], label="free jet delta ~ x^(2/3)")
    ax[2].loglog(x, ch09.wall_jet(x, 0.5, C, finf, nu)["delta"], color=COLORS["orange"], label="wall jet delta ~ x^(3/4)")
    ax[2].set_xlabel("x [m]")
    ax[2].set_ylabel("delta [m]")
    ax[2].legend(fontsize=8)
    ax[2].set_title("the wall makes the jet spread faster")
    save(fig, out, "c12_jets")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
