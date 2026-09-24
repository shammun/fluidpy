"""§7.2–7.3 (C03, C04, C07): the dispersion relation ω = √(gk tanh kH) (7.28) and what follows from it — cosh, sinh,
tanh (Fig. 7.7 remake), c(λ) with the deep (7.45) and shallow (7.49) limits and the depth-regime errors, the pressure
response cosh k(z + H)/cosh kH (7.31) with its deep (7.48) and hydrostatic (7.52) limits, chord and tangent of ω(k)
(Fig. 7.14), and the capillary–gravity c(λ) with its minimum (7.58) (Fig. 7.10 remake); ocean numbers; the inverse
dispersion checked against the explicit approximations of Fenton & McKee (1990) and Guo (2002).

Run: ``.venv/Scripts/python.exe scripts/ch07_dispersion.py --no-show``
Figures → outputs/ch07/c03_dispersion.png, c07_capillary.png.
"""
from __future__ import annotations

import numpy as np

from ch07_drawings import parse_args, save, setup

from fluidpy import ch01_introduction as ch01
from fluidpy import ch07_gravity_waves as ch07
from fluidpy.core.style import COLORS


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    g = ch07.G_BOOK
    for H in (100.0, 4000.0, np.inf):
        print(f"T = 10 s on H = {H:g} m: λ = {ch07.wavelength_from_period(10.0, H):.2f} m")
    lam_ts = ch07.wavelength_from_period(20 * 60.0, 4000.0)
    print(f"tsunami T = 20 min on H = 4 km: λ = {lam_ts / 1e3:.1f} km, c = {ch07.phase_speed(2 * np.pi / lam_ts, 4000.0):.1f}"
          f" m/s (√(gH) = {np.sqrt(g * 4000):.1f} m/s)")
    for kH in (2.0, 2 * np.pi * 0.07):
        r = ch07.depth_regime(kH, 1.0)
        print(f"kH = {kH:.4f} (H/λ = {r['H_over_lambda']:.4f}): deep-limit error {100 * r['deep_error']:.2f} %,"
              f" shallow-limit error {100 * r['shallow_error']:.2f} %")
    om = np.logspace(-4, 3, 200)
    for H in (1e-3, 1.0, 1e4):
        k = ch07.wavenumber_from_omega(om, H)
        res = np.max(np.abs(ch07.omega_gravity(k, H) / om - 1))
        fm = np.max(np.abs(ch07.fenton_mckee_kh(om, H) / (k * H) - 1))
        gu = np.max(np.abs(ch07.guo_kh(om, H) / (k * H) - 1))
        print(f"H = {H:g} m: inverse residual {res:.1e}; max error Fenton–McKee {100 * fm:.2f} %, Guo {100 * gu:.2f} %")
    kk = 2 * np.pi / 50.0
    print(f"λ = 50 m on H = 10 m: cosh kH = {np.cosh(kk * 10):.3f}; p′ amplitude at the bottom / surface = "
          f"{ch07.pressure_response(kk, -10.0, 10.0):.3f} (ρga = {1000 * g:.0f} Pa → {1000 * g * ch07.pressure_response(kk, -10.0, 10.0):.0f} Pa)")
    print(f"deep water at z = −λ/2: e^(−π) = {ch07.pressure_response(1.0, -np.pi):.4f}")

    fig, ax = plt.subplots(2, 2, figsize=(12, 8.5))
    xs = np.linspace(0, 2.3, 200)
    a0 = ax[0, 0]
    a0.plot(xs, np.cosh(xs), label="cosh x", color=COLORS["accent"])
    a0.plot(xs, np.sinh(xs), label="sinh x", color=COLORS["teal"])
    a0.plot(xs, np.tanh(xs), label="tanh x", color=COLORS["orange"])
    a0.axhline(1, color=COLORS["muted"], ls="--", lw=0.8)
    a0.plot(2.0, np.tanh(2.0), "o", color=COLORS["orange"])
    a0.annotate(f"tanh 2 = {np.tanh(2):.5f}", (2.0, np.tanh(2.0)), (1.2, 0.55), fontsize=8,
                arrowprops=dict(arrowstyle="->", lw=0.8))
    a0.set_ylim(0, 3)
    a0.set_xlabel("x")
    a0.set_title("hyperbolic functions (Fig. 7.7 remake)", fontsize=10)
    a0.legend(fontsize=8)
    a1 = ax[0, 1]
    lam = np.logspace(0, 6, 400)
    for H, col in ((10.0, COLORS["blue"]), (100.0, COLORS["teal"]), (4000.0, COLORS["accent"])):
        a1.loglog(lam, ch07.phase_speed(2 * np.pi / lam, H), color=col, label=f"H = {H:g} m (7.29)")
        a1.axhline(np.sqrt(g * H), color=col, ls=":", lw=0.8)
    a1.loglog(lam, np.sqrt(g * lam / (2 * np.pi)), color=COLORS["muted"], ls="--", label="deep √(gλ/2π) (7.45)")
    a1.set_xlabel("wavelength λ [m]")
    a1.set_ylabel("phase speed c [m/s]")
    a1.set_title("c(λ): dispersive when deep, √(gH) (dotted, (7.49)) when shallow", fontsize=10)
    a1.legend(fontsize=8)
    a2 = ax[1, 0]
    H = 10.0
    z = np.linspace(-H, 0, 200)
    for lam_i, col in ((5.0, COLORS["accent"]), (30.0, COLORS["teal"]), (300.0, COLORS["blue"])):
        k_i = 2 * np.pi / lam_i
        a2.plot(ch07.pressure_response(k_i, z, H), z, color=col, label=f"λ = {lam_i:g} m (kH = {k_i * H:.2f})")
    a2.plot(np.exp(2 * np.pi / 5.0 * z), z, color=COLORS["muted"], ls="--", lw=1, label="deep $e^{kz}$ (7.48), λ = 5 m")
    a2.axvline(1.0, color=COLORS["muted"], ls=":", lw=1, label="hydrostatic (7.52)")
    a2.set_xlabel("p′ amplitude / ρga  (cosh k(z+H)/cosh kH, (7.31))")
    a2.set_ylabel("z [m]")
    a2.set_title("pressure under a wave on H = 10 m", fontsize=10)
    a2.legend(fontsize=8)
    a3 = ax[1, 1]
    k = np.linspace(1e-4, 0.12, 300)
    H = 30.0
    w = ch07.omega_gravity(k, H)
    k0 = 0.07
    w0 = ch07.omega_gravity(k0, H)
    cg0 = ch07.group_velocity(k0, H)
    a3.plot(k, w, color=COLORS["ink"], label="ω(k), H = 30 m (7.28)")
    a3.plot([0, k0], [0, w0], color=COLORS["orange"], label=f"chord: slope c = {w0 / k0:.2f} m/s")
    a3.plot(k, w0 + cg0 * (k - k0), color=COLORS["accent"], ls="--", label=f"tangent: slope c_g = {cg0:.2f} m/s (7.69)")
    a3.plot(k0, w0, "o", color=COLORS["orange"])
    a3.set_ylim(0, w.max() * 1.1)
    a3.set_xlabel("k [rad/m]")
    a3.set_ylabel("ω [rad/s]")
    a3.set_title("phase speed = chord, group speed = tangent (Fig. 7.14 remake)", fontsize=10)
    a3.legend(fontsize=8)
    assert cg0 < w0 / k0
    save(fig, out, "c03_dispersion")

    sig_i = float(ch01.surface_tension_water(293.15))
    for sig, rho in ((0.073, 1000.0), (sig_i, 998.2)):
        m = ch07.capillary_minimum(sig, rho, g)
        gm = ch07.min_group_velocity(sig, rho, g)
        print(f"σ = {sig * 1e3:.2f} mN/m, ρ = {rho:g}: c_min = {100 * m['c_min']:.2f} cm/s at λ_m = {100 * m['lam_m']:.3f} cm;"
              f" c_g,min = {100 * gm['cg_min']:.2f} cm/s at λ = {100 * gm['lam']:.2f} cm")
    from scipy.optimize import minimize_scalar
    r = minimize_scalar(lambda L: ch07.phase_speed(2 * np.pi / L, np.inf, g, 0.073, 1000.0), bounds=(1e-3, 0.1),
                        method="bounded", options={"xatol": 1e-12})
    print(f"  numerical minimum of (7.57): λ = {100 * r.x:.4f} cm, c = {100 * r.fun:.4f} cm/s")
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    lam = np.logspace(-3.3, 2, 500)
    H = 1.0
    ax.loglog(lam, ch07.phase_speed(2 * np.pi / lam, H, g, 0.0727, 1000.0), color=COLORS["accent"],
              label="capillary–gravity, H = 1 m (7.57)")
    ax.loglog(lam, ch07.phase_speed(2 * np.pi / lam, H, g), color=COLORS["blue"], ls="--", lw=1, label="σ = 0 (7.29)")
    ax.loglog(lam, np.sqrt(2 * np.pi * 0.0727 / (1000.0 * lam)), color=COLORS["rose"], ls="--", lw=1,
              label="pure capillary √(2πσ/ρλ) (7.60)")
    ax.axhline(np.sqrt(g * H), color=COLORS["muted"], ls=":", lw=1, label="√(gH)")
    m = ch07.capillary_minimum(0.0727, 1000.0, g)
    ax.plot(m["lam_m"], m["c_min"], "o", color=COLORS["orange"],
            label=f"c_min = {100 * m['c_min']:.1f} cm/s at λ_m = {100 * m['lam_m']:.2f} cm (7.58)")
    ax.set_ylim(0.1, 5)
    ax.set_xlabel("wavelength λ [m]")
    ax.set_ylabel("phase speed c [m/s]")
    ax.set_title("two restoring forces leave a slowest wave (Fig. 7.10 remake, our code)", fontsize=10)
    ax.legend(fontsize=8)
    save(fig, out, "c07_capillary")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
