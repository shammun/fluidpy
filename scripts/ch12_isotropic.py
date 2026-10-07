"""§12.6 (C05): homogeneous isotropic turbulence on a KINEMATIC three-dimensional synthetic field (random phases, prescribed
spectrum, divergence-free — no cascade): the isotropy report (12.36)-(12.37), f(r) and g(r) (12.38) against
g = f + (r/2) f' (12.41), the scale ratios, and the dissipation rate three ways (12.42)-(12.43).  A 2-D field is shown only
to make the point that it obeys a different relation (g = d(rf)/dr).

Run: ``.venv/Scripts/python.exe scripts/ch12_isotropic.py --no-show [--fast]``
Figure -> outputs/ch12/c05_isotropic_f_g.png.   Field: 64^3 (fast 48^3), box 2*pi m, E(K) ~ K^4 exp(-2(K/K0)^2), K0 = 8, seed 0.
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

    n, Lbox, K0, nu = (48 if args.fast else 64), 2.0 * np.pi, 8.0, 1.0e-3
    E = lambda K: 1e-3 * K ** 4 * np.exp(-2.0 * (K / K0) ** 2)  # noqa: E731
    dx = Lbox / n
    with Timer("3-D synthetic field"):
        u, v, w = TS.synthetic_solenoidal_field(n, Lbox, E, seed=0, dim=3)
    print(f"rms divergence: solenoidal field {ch12.divergence_rms(u, v, dx, w=w):.2e} 1/s, white noise "
          f"{ch12.divergence_rms(*TS.white_noise_field(n, seed=0, dim=3)[:2], dx, w=TS.white_noise_field(n, seed=0, dim=3)[2]):.2e} 1/s")
    rep = ch12.isotropy_report(u, v, w=w, dx=dx)
    print(f"normal-stress ratios {np.round(rep['normal_stress_ratios'], 3)}, shear coefficients "
          f"{ {k: round(val, 3) for k, val in rep['shear_coefficients'].items()} }")
    print(f"<(du_i/dx_j)^2>/<(du_i/dx_i)^2> = {rep['transverse_over_longitudinal']:.3f}  (isotropic incompressible 3-D: 2)")
    u2 = float(np.mean(rep["variances"]))
    r, f, g = TS.longitudinal_transverse_correlation(u, v, dx, w=w)
    g_pred = ch12.transverse_from_longitudinal(r, f, dim=3)
    m = r < 1.2
    print(f"max |g - (f + r f'/2)| for r < 1.2 m: {np.max(np.abs(g - g_pred)[m]):.3f}   (the 2-D form f + r f' misses by "
          f"{np.max(np.abs(g - ch12.transverse_from_longitudinal(r, f, dim=2))[m]):.3f})")
    lam_f, lam_g = TS.taylor_microscale(r, f, fit_points=4), TS.taylor_microscale(r, g, fit_points=4)
    Lf, Lg = TS.integral_scale(r, f), TS.integral_scale(r, g)
    print(f"lambda_f = {lam_f:.4f} m, lambda_g = {lam_g:.4f} m, ratio {lam_g / lam_f:.3f} (1/sqrt 2 = 0.707); "
          f"Lambda_f = {Lf:.4f} m, Lambda_g (to first zero) = {Lg:.4f} m")
    Kq = np.linspace(1e-3, 6 * K0, 20001)
    eps_spec = 2.0 * nu * np.trapezoid(Kq ** 2 * E(Kq), Kq)
    eps_grad = float(ch12.dissipation_isotropic(nu, dudx_sq=float(np.mean(rep["longitudinal"]))))
    eps_full = ch12.dissipation_rate(ch12.velocity_gradient_samples(u, v, w, dx), nu)
    print("dissipation rate [m2/s3]:")
    print(f"  full definition (12.42) on the field      {eps_full:.5e}")
    print(f"  15 nu <(du1/dx1)^2>                        {eps_grad:.5e}")
    print(f"  30 nu u2/lambda_f^2, 15 nu u2/lambda_g^2   {float(ch12.dissipation_isotropic(nu, u2=u2, lambda_f=lam_f)):.5e}, "
          f"{float(ch12.dissipation_isotropic(nu, u2=u2, lambda_g=lam_g)):.5e}")
    print(f"  2 nu int K^2 E dK (prescribed spectrum)    {eps_spec:.5e}")
    print(f"R_lambda (lambda_g) = {ch12.taylor_reynolds_number(u2, lam_g, nu):.1f}")
    sc = ch12.isotropic_scales(np.linspace(0, 6, 6001), lambda rr: np.exp(-rr ** 2 / 0.5 ** 2))
    print(f"Gaussian f: Lambda_g/Lambda_f = {sc['Lambda_ratio']:.6f}, lambda_g/lambda_f = {sc['lambda_ratio']:.6f}")
    gm = ch12.gradient_moments_isotropic_sympy()
    print(f"gradient moments in units of u2/lambda_f^2: {gm['m11']}, {gm['m12']}, {gm['cross']} -> factor {gm['eps_factor']}")
    u2d, v2d = TS.synthetic_solenoidal_field(128, Lbox, E, seed=0, dim=2)
    print(f"2-D solenoidal field for contrast: gradient ratio {ch12.isotropy_report(u2d, v2d, dx=Lbox / 128)['transverse_over_longitudinal']:.2f} (3, not 2)")

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    ax[0].plot(r, f, color=COLORS["accent"], label="f(r) longitudinal")
    ax[0].plot(r, g, color=COLORS["teal"], label="g(r) transverse")
    ax[0].plot(r, g_pred, "k--", lw=1.0, label="f + (r/2) f'  (12.41)")
    ax[0].plot(r, ch12.transverse_from_longitudinal(r, f, dim=2), color=COLORS["rose"], ls=":", label="f + r f'  (2-D relation: wrong here)")
    ax[0].axhline(0, color=COLORS["grid"])
    ax[0].set(xlim=(0, 1.5), xlabel="r [m]", ylabel="correlation coefficient", title="3-D kinematic field (no cascade)")
    ax[0].legend(fontsize=8)
    Ksh, Esh = TS.shell_spectrum((u, v, w), Lbox)
    ax[1].loglog(Ksh, Esh, "o", ms=3, color=COLORS["accent"], label="shell spectrum of the field")
    ax[1].loglog(Ksh, E(Ksh), "k--", lw=1.0, label="prescribed E(K)")
    ax[1].set(ylim=(1e-6, None), xlabel="K [rad/m]", ylabel=r"E(K) [m$^3$/s$^2$]", title="the spectrum the field was built with")
    ax[1].legend(fontsize=8)
    save(fig, out, "c05_isotropic_f_g")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
