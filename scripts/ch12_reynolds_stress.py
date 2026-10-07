"""§12.5 (C04): the Reynolds stress — displaced parcels in a mean shear (Fig. 12.6 argument as a Monte-Carlo count), the
(u, v) scatter with its principal axes (Fig. 12.8 analogue), the 2 x 2 stress tensor, the closure count, and the symbolic
Reynolds averaging of the equations of motion (12.27)-(12.35).

Run: ``.venv/Scripts/python.exe scripts/ch12_reynolds_stress.py --no-show [--fast]``
Figure -> outputs/ch12/c04_reynolds_stress.png.
Our numbers: dU/dy = 20 1/s, l_rms = 5 mm, v_rms = 0.3 m/s, air (rho = 1.2 kg/m3, mu = 1.8e-5 Pa s), seed 0.
"""
from __future__ import annotations

import numpy as np

from ch12_common import COLORS, Timer, finish, parse_args, save, setup

from fluidpy import ch12_turbulence as ch12
from fluidpy.core import turbstats as TS


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    dUdy, l_rms, v_rms, rho, mu = 20.0, 5e-3, 0.3, 1.2, 1.8e-5
    n = 20000 if args.fast else 200000
    print("displaced parcels: mean(uv) against -c v_rms l_rms dU/dy")
    for c in (1.0, 0.6, 0.0):
        d = ch12.displaced_parcel_uv(dUdy, l_rms, v_rms, n=n, seed=0, correlation=c, u_extra_rms=0.05)
        print(f"  c = {c:.1f}: mean(uv) = {d['uv']:+.5f} +- {d['stderr']:.5f} m2/s2, expected {d['expected']:+.5f}, r_uv = {d['r_uv']:+.3f}")
    d = ch12.displaced_parcel_uv(dUdy, l_rms, v_rms, n=n, seed=0, correlation=0.6, u_extra_rms=0.05)
    drev = ch12.displaced_parcel_uv(-dUdy, l_rms, v_rms, n=n, seed=0, correlation=0.6, u_extra_rms=0.05)
    print(f"  reversed shear: mean(uv) = {drev['uv']:+.5f} (sign flips)")
    samples = np.stack([d["u"], d["v"]], axis=1)
    tau_R = TS.reynolds_stress(samples, rho0=rho)
    cov = TS.velocity_covariance(samples)
    print(f"Reynolds stress tensor -rho <u_i u_j> [Pa]:\n{tau_R}")
    print(f"  shear component {tau_R[0, 1]:.5f} Pa against the viscous stress mu dU/dy = {mu * dUdy:.2e} Pa "
          f"(ratio {tau_R[0, 1] / (mu * dUdy):.0f}); implied eddy viscosity {d['nu_T']:.2e} m2/s = {d['nu_T'] / (mu / rho):.0f} nu")
    gradU = np.array([[0.0, dUdy], [0.0, 0.0]])
    print(f"mean stress tensor (12.30), P = 0 [Pa]:\n{ch12.mean_stress_tensor(0.0, gradU, mu, rho, cov)}")
    print(f"shear production -<u_i u_j> dU_i/dx_j = {ch12.shear_production(cov, gradU):.4f} m2/s3")
    for row in ch12.closure_count(3):
        print(f"  closure level {row['level']}: {row['equations']} equations, {row['unknowns']} velocity-moment unknowns"
              f" ({row['unknowns_with_extra']} with pressure and dissipation correlations)")
    names = ("rans", "reynolds_stress_budget") if not args.fast else ("rans",)
    for name in names:
        with Timer(f"sympy {name}"):
            s = ch12.sympy_summary(name)
        print(f"  {name}: {s['checks']}")

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
    y = np.linspace(-0.02, 0.02, 9)
    ax[0].plot(dUdy * y, y * 1e3, color=COLORS["accent"])
    sub = slice(0, 60)
    col = np.where(d["u"][sub] * d["v"][sub] < 0, COLORS["orange"], COLORS["muted"])
    ax[0].quiver(np.zeros(60), np.zeros(60), d["u"][sub], d["v"][sub], color=col, angles="xy", scale_units="xy", scale=1.0, width=0.004)
    ax[0].set(xlabel="U(y) and fluctuation u [m/s]", ylabel="y [mm]", title="arrivals at y = 0: orange = uv < 0")
    iso = TS.correlated_pair(4000, 0.0, 0.2, 0.2, seed=2)
    ax[1].scatter(iso[0], iso[1], s=2, color=COLORS["muted"])
    ax[1].set(aspect="equal", xlabel="u [m/s]", ylabel="v [m/s]", title=r"isotropic: $\overline{uv}=0$")
    ax[2].scatter(d["u"][:4000], d["v"][:4000], s=2, color=COLORS["orange"])
    w, vec = np.linalg.eigh(cov)
    for k in range(2):
        ax[2].plot([0, 2 * np.sqrt(w[k]) * vec[0, k]], [0, 2 * np.sqrt(w[k]) * vec[1, k]], color=COLORS["ink"], lw=2)
    ax[2].set(aspect="equal", xlabel="u [m/s]", ylabel="v [m/s]", title=rf"in a shear: $\overline{{uv}}$ = {d['uv']:.4f} m$^2$/s$^2$")
    save(fig, out, "c04_reynolds_stress")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
