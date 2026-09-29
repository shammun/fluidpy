"""§8.4 (C10, Examples 8.5–8.6, Figs. 8.14–8.15): the thickening vortex sheet (ω_z Gaussian, u = U erf) with its
circulation −2U conserved, the temporally developing boundary layer C_f = (2/√π)Re_x^{−1/2}, the decaying line vortex
(Lamb–Oseen) and its spin-up twin (Exercise 8.26), and parity with ch05 / core.vortices; exponents from the
similarity engine and the collapse error.

Run: ``.venv/Scripts/python.exe scripts/ch08_vortex_diffusion.py --no-show``
Figure → outputs/ch08/c10_vortex_diffusion.png.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import quad

from ch08_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch05_vorticity_dynamics as ch05
from fluidpy import ch08_laminar_flow as ch08
from fluidpy.core import vortices as VX


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    U, nu = 0.1, 1e-6
    y = np.linspace(-0.01, 0.01, 401)
    fig, ax = plt.subplots(1, 4, figsize=(17, 4.4))
    for t in (5.0, 20.0, 80.0):
        u, w = ch08.vortex_sheet_diffusion(y, t, U, nu)
        u5, w5 = ch05.diffusing_vortex_sheet(y, t, -2 * U, nu)
        circ = quad(lambda s: float(ch08.vortex_sheet_diffusion(s, t, U, nu)[1]), -np.inf, np.inf)[0]
        print(f"t = {t:.0f} s: ∫ω dy = {circ:.10f} m/s (−2U = {-2 * U}); parity ch05(γ = −2U) |Δu| = "
              f"{np.max(np.abs(u - u5)):.1e}, |Δω| = {np.max(np.abs(w - w5)):.1e}; width(±0.95U) = "
              f"{float(ch08.transition_width(t, nu)) * 1000:.3f} mm")
        ax[0].plot(w, y * 1000, lw=2, label=f"t = {t:.0f} s")
        ax[1].plot(u / U, y / (2 * np.sqrt(nu * t)), lw=2)
    print(f"±0.95U at η = ±{float(ch08.transition_width(1.0, 1.0)) / 2:.4f}, width coefficient "
          f"{float(ch08.transition_width(1.0, 1.0)):.4f}")
    wb = ch08.temporal_bl_wall_stress(10.0, U, nu)
    print("temporal boundary layer: " + ", ".join(f"{k} = {v:.5g}" for k, v in wb.items()))
    ax[0].set_xlabel("ω_z [1/s]")
    ax[0].set_ylabel("y [mm]")
    ax[0].legend(fontsize=8)
    ax[0].set_title("thickening vortex sheet (Fig. 8.14 left)")
    ax[1].set_xlabel("u/U")
    ax[1].set_ylabel("y/(2√(νt))")
    ax[1].set_ylim(-2.5, 2.5)
    ax[1].set_title("velocity in similarity coordinates")
    Gam = 1e-3
    r = np.linspace(0, 0.06, 600)
    for nut in (0.5e-4, 1e-4, 3e-4):
        t = nut / nu
        ud = ch08.line_vortex_decay(r, t, Gam, nu)
        dg = np.max(np.abs(ud - VX.gaussian_vortex(r, Gam, 2 * np.sqrt(nu * t))[0]))
        ax[2].plot(r * 100, ud * 1000, lw=2, label=f"νt = {nut:.1e} m²")
        ax[3].plot(r * 100, ch08.line_vortex_spinup(np.maximum(r, 1e-4), t, Gam, nu) * 1000, lw=2)
        print(f"νt = {nut:.1e} m²: max u_θ = {np.max(ud) * 1000:.4f} mm/s at r = {r[np.argmax(ud)] * 1000:.3f} mm; "
              f"parity with gaussian_vortex(σ = 2√(νt)) {dg:.1e}")
    rr = np.linspace(1e-3, 0.06, 300)
    ax[2].plot(rr * 100, Gam / (2 * np.pi * rr) * 1000, color=COLORS["muted"], ls="--", label="νt = 0 (ideal)")
    ax[2].set_ylim(0, 60)
    ax[2].set_xlabel("r [cm]")
    ax[2].set_ylabel("u_θ [mm/s]")
    ax[2].legend(fontsize=7)
    ax[2].set_title("decay of a line vortex (Fig. 8.15)")
    ax[3].set_ylim(0, 60)
    ax[3].set_xlabel("r [cm]")
    ax[3].set_ylabel("u_θ [mm/s]")
    ax[3].set_title("spin-up (Exercise 8.26)")
    for case, nm in (("vortex_sheet", (0.5, 0.5)), ("line_vortex", (1.0, 0.5))):
        d = ch08.similarity_reduce_sympy(case)
        print(f"{case}: ODE {d['ode']}; exponents n = {d.get('n', 'Γ/2πr prefactor')}; collapse error at the right "
              f"exponents {ch08.similarity_collapse_error(case, *nm):.1e}, at n + 0.1: "
              f"{ch08.similarity_collapse_error(case, nm[0] + 0.1, nm[1]):.3f}")
    save(fig, out, "c10_vortex_diffusion")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
