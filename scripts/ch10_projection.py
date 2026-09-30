"""§10.4 (C10–C12): operator splitting and the Θ-scheme orders, artificial compressibility, the projection on the staggered
grid (divergence before/after, curl-free correction), the checkerboard null space, MAC code verification (Taylor–Green,
Poiseuille), MAC stability limits vs a 2-D von Neumann scan, and weakly compressible MacCormack on Taylor–Green.

Run: ``.venv/Scripts/python.exe scripts/ch10_projection.py --no-show``   Figure -> outputs/ch10/c12_projection.png.
"""
from __future__ import annotations

import numpy as np

from ch10_common import COLORS, Timer, finish, parse_args, save, setup

from fluidpy import ch10_computational_fluid_dynamics as ch10


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    print("splitting orders: Marchuk–Yanenko", round(ch10.splitting_order("marchuk_yanenko"), 3), "| Θ (θ = 1 − 1/√2)",
          round(ch10.splitting_order("theta"), 3), "| Θ (θ = 1/4)", round(ch10.splitting_order("theta", theta_split=0.25), 3))
    th = ch10.theta_scheme_amplification_sympy()
    print("Θ-scheme Δt² coefficient:", th["z2_coeff"], "→ root θ =", th["theta_root"], "=", float(th["theta_root"]))
    with Timer("artificial compressibility"):
        ac = ch10.artificial_compressibility_channel()
    print(f"artificial compressibility (10.95): max|u − Poiseuille| = {ac['max_err']:.2e}, max|∇·u| = {ac['div_final']:.2e} "
          f"after {ac['iterations']} pseudo-time steps")
    for kind in ("divergent", "shear"):
        print(f"projection of the '{kind}' field on 8²:", {k: f"{v:.3e}" for k, v in ch10.mac_projection_summary(8, kind).items()})
    print("gradient null space on 8² periodic: collocated", ch10.gradient_null_space(8, 8, "collocated"),
          "staggered", ch10.gradient_null_space(8, 8, "staggered"))
    tg = ch10.MAC.taylor_green(32, 100.0, 0.5)
    print(f"MAC Taylor–Green 32², Re = 100, t = 0.5: max err u {tg['err_u']:.2e}, p {tg['err_p']:.2e}, order {tg['order']:.2f}")
    po = ch10.MAC.channel_poiseuille(16, 16, 10.0, -1.0)
    print(f"MAC channel Poiseuille (quadratic ghosts): max err {po['max_err']:.1e}")
    for Re, dx in ((100.0, 1 / 32), (100.0, 1 / 64), (10.0, 1 / 32)):
        lim = ch10.MAC.dt_limit(1.0, 0.0, Re, dx, safety=1.0)
        print(f"MAC limit (10.127)–(10.128), Re = {Re:g}, Δx = {dx:.4f}: Δt = {lim:.5g}; 2-D von Neumann max|G| at 0.95×: "
              f"{ch10.ftcs2d_max_amplification(1.0, 0.0, Re, dx, 0.95 * lim):.6f}, at 1.05×: "
              f"{ch10.ftcs2d_max_amplification(1.0, 0.0, Re, dx, 1.05 * lim):.6f}")
    for n in (32, 64):
        r = ch10.MCK.wc_taylor_green(n, 0.05, 100.0, 1.0)
        print(f"weakly compressible MacCormack Taylor–Green {n}², Ma = 0.05: max err u {r['err_u']:.2e}, "
              f"decay {r['decay_measured']:.5f} vs exact {r['decay_exact']:.5f}, mass drift {r['mass_drift']:.1e}")
    # figure: projection stages of a random-free divergent field
    g = ch10.MacGrid(16, 16)
    c = ch10.face_coordinates(g)
    us = np.sin(np.pi * c["xu"]) * np.cos(np.pi * c["yu"]) + 0.5 * np.sin(2 * np.pi * c["xu"])
    vs = np.cos(np.pi * c["xv"]) * np.sin(np.pi * c["yv"])
    us[:, 0] = us[:, -1] = 0.0
    vs[0, :] = vs[-1, :] = 0.0
    st = ch10.projection_stages(us, vs, g, 1.0)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    im = ax[0].imshow(st["div_before"], origin="lower", extent=(0, 1, 0, 1), cmap="RdBu_r")
    fig.colorbar(im, ax=ax[0])
    ax[0].set_title(f"∇·u* before (max {st['max_div_before']:.2f})")
    im = ax[1].imshow(st["p"], origin="lower", extent=(0, 1, 0, 1), cmap="viridis")
    fig.colorbar(im, ax=ax[1])
    ax[1].set_title("pressure from (10.124) (zero mean)")
    im = ax[2].imshow(st["div_after"], origin="lower", extent=(0, 1, 0, 1), cmap="RdBu_r")
    fig.colorbar(im, ax=ax[2])
    ax[2].set_title(f"∇·u after (max {st['max_div_after']:.1e})")
    for a in ax:
        a.set_xlabel("x / L")
        a.set_ylabel("y / L")
    print(f"projection 16²: Σ rhs ΔxΔy = {st['rhs_sum']:.1e}, max|curl of the correction| = {np.abs(st['curl_correction']).max():.1e}")
    save(fig, out, "c12_projection")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
