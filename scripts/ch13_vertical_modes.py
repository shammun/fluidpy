"""Chapter 13, §13.9 — vertical normal modes of a stratified layer.

* uniform N: the exact roots of tan(NH/c) = cN/g against the numerical eigenvalues (finite volumes, order 2; Chebyshev;
  shooting), the barotropic and baroclinic speeds, equivalent depths and internal Rossby radii, the rigid-lid error;
* a thermocline profile: mode shapes, speeds, the WKB estimate, orthogonality with weight 1 (free surface included).

Run: ``.venv/Scripts/python.exe scripts/ch13_vertical_modes.py --no-show``   Figures -> outputs/ch13/modes_*.png
"""
from __future__ import annotations

import numpy as np
from ch13_common import finish, parse_args, save, setup

from fluidpy import ch13_geophysical_fluid_dynamics as ch13


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    I = ch13.illustrative_inputs()
    N, H = I["ocean_N"], I["ocean_H"]
    f = ch13.coriolis_parameter(I["lat"])
    ex = ch13.modes_uniform_N(N, H, n_modes=4, nz=401)
    X = N * H / ex.c
    print(f"uniform N = {N} rad/s, H = {H:.0f} m: roots X_n = NH/c_n = {np.round(X, 5)} (the first is far below 1: slip #6)")
    print(f"  c_0 = {ex.c[0]:.3f} m/s (sqrt(gH) = {ch13.long_wave_speed(H):.3f}, difference {ex.c[0] / ch13.long_wave_speed(H) - 1:+.1e}); "
          f"c_1 = {ex.c[1]:.4f} m/s (NH/pi = {ch13.baroclinic_mode_speed(N, H, 1):.4f}, rigid-lid error {ch13.rigid_lid_error(N, H):+.2e})")
    for n in range(4):
        print(f"  mode {n}: c = {ex.c[n]:9.4f} m/s, equivalent depth {ex.He[n]:10.4f} m, Rossby radius at 35° "
              f"{ch13.rossby_radius(ex.c[n], f) / 1e3:8.1f} km")
    errs = []
    for nz in (101, 201, 401, 801):
        z = np.linspace(-H, 0.0, nz)
        errs.append(abs(ch13.vertical_modes(z, N * N, n_modes=3).c[1] / ex.c[1] - 1))
    print(f"  finite volumes: relative error of c_1 {['%.1e' % e for e in errs]} at 101…801 nodes, observed order "
          f"{np.round(np.log2(np.array(errs[:-1]) / np.array(errs[1:])), 2)}")
    ch_ = ch13.vertical_modes(ex.z, lambda q: N * N + 0.0 * q, n_modes=3, method="cheb")
    sh = [ch13.vertical_modes_shooting(ex.z, lambda q: N * N, n=n)[0] for n in range(3)]
    print(f"  Chebyshev route: relative errors {['%.1e' % abs(a / b - 1) for a, b in zip(ch_.c, ex.c)]}; shooting: "
          f"{['%.1e' % abs(a / b - 1) for a, b in zip(sh, ex.c)]}")
    fd = ch13.vertical_modes(ex.z, N * N, n_modes=4)
    print(f"  orthogonality (weight 1, free surface): numerical modes {np.max(np.abs(ch13.orthogonality_matrix(fd) - np.eye(4))):.1e}, "
          f"exact modes on the same nodes {np.max(np.abs(ch13.orthogonality_matrix(ex) - np.eye(4))):.1e} (trapezoid error only)")

    z = np.linspace(-H, 0.0, 801)
    N2 = ch13.thermocline_N2(z)
    m = ch13.vertical_modes(z, N2, n_modes=4)
    mr = ch13.vertical_modes(z, N2, n_modes=3, lid="rigid")
    print(f"thermocline profile (ours): c = {np.round(m.c, 4)} m/s; rigid lid: {np.round(mr.c, 4)}")
    for n in (1, 2, 3):
        w = ch13.wkb_mode_speed(z, np.sqrt(N2), n)
        print(f"  mode {n}: c = {m.c[n]:.4f} m/s, H_e = {m.He[n]:.4f} m, radius at 35° {ch13.rossby_radius(m.c[n], f) / 1e3:.1f} km; "
              f"WKB estimate {w:.4f} m/s ({100 * (w / m.c[n] - 1):+.0f} %); zero crossings {int(np.sum(m.psi[n][:-1] * m.psi[n][1:] < 0))}")
    a = ch13.project(m, 0.3 * m.psi[1] - 0.1 * m.psi[2])
    print(f"  project/reconstruct: coefficients {np.round(a, 12)}")
    save(ch13.fig_mode_roots(), out, "modes_roots")
    save(ch13.fig_vertical_modes(z, N2), out, "modes_thermocline")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
