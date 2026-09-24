"""§7.2 (C05): particle orbits and streamlines of the linear wave — ellipses (7.36) that are circles in deep water
(7.46) and flat ellipses in shallow water (7.50) (Fig. 7.4 remake for kH = 3, 1, 0.3), one particle's clockwise orbit
from the path-line equations (7.33) against the linearised (7.35) (Fig. 7.3), and the instantaneous streamlines of
(7.37) (Fig. 7.5 remake); parity of the fields with ch04's test field.

Run: ``.venv/Scripts/python.exe scripts/ch07_orbits.py --no-show``
Figures → outputs/ch07/c05_orbits.png, c05_streamlines.png.
"""
from __future__ import annotations

import numpy as np

from ch07_drawings import orbit_ghosts, parse_args, save, setup, tank, wave_surface

from fluidpy import ch04_conservation_laws as ch04
from fluidpy import ch07_gravity_waves as ch07
from fluidpy.core.style import COLORS


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    k, a = 1.0, 0.12
    x = np.linspace(0, 10, 7)
    z = np.linspace(-3, 0, 7)
    d = max(np.max(np.abs(np.asarray(ch07.wave_fields(x, z, 0.7, 0.05, 1.0, 4.0)[n]) -
                          np.asarray(ch04.linear_wave_surface(x, z, 0.7, 0.05, 1.0, 9.81, 4.0)[n])))
            for n in ("eta", "u", "w"))
    print(f"parity with ch04.linear_wave_surface (H = 4 m): max |Δ| = {d:.1e}")
    fig, ax = plt.subplots(1, 3, figsize=(15, 5.2), gridspec_kw={"width_ratios": [1, 1, 1]})
    for axi, kH in zip(ax, (3.0, 1.0, 0.3)):
        H = kH / k
        a = 0.06 * H  # same a/H in every panel (true aspect ratio below)
        lam = 2 * np.pi / k
        xs = np.linspace(0, lam / 2, 300)
        tank(axi, H, lam / 2, label=False)
        wave_surface(axi, xs, ch07.wave_fields(xs, 0.0, 0.0, a, k, H)["eta"], fill_to=-H)
        z0s = -H * np.array([0.1, 0.35, 0.6, 0.85])
        x0s = np.full_like(z0s, lam / 4)
        sa = ch07.orbit_semi_axes(z0s, a, k, H)
        orbit_ghosts(axi, x0s, z0s, sa["A"], sa["B"])
        print(f"kH = {kH}: A = {np.round(sa['A'], 4)}, B = {np.round(sa['B'], 4)} m; focal half-distance "
              f"{sa['focal_half']:.4f} m = √(A² − B²) = {np.sqrt(sa['A'][0] ** 2 - sa['B'][0] ** 2):.4f} m")
        T = 2 * np.pi / ch07.omega_gravity(k, H)
        p = ch07.particle_path(lam / 4, z0s[0], np.linspace(0, T, 200), a, k, H, model="linear")
        xi, ze = ch07.orbit_linear(lam / 4, z0s[0], p["t"], a, k, H)
        err = max(np.max(np.abs(p["x"] - (lam / 4 + xi))), np.max(np.abs(p["z"] - (z0s[0] + ze))))
        area = 0.5 * np.sum(xi[:-1] * ze[1:] - xi[1:] * ze[:-1])
        print(f"  particle_path(linear) vs orbit_linear: {err:.1e} m; signed area {area:+.4f} m² (negative = clockwise)")
        assert area < 0
        axi.set_title(f"kH = {kH}: " + ("deep: circles (7.46)" if kH >= 2 else "shallow: flat ellipses (7.50)"
                                        if kH < 0.5 else "intermediate: ellipses (7.36)")
                      + f"\na = {a:.3g} m, true aspect ratio", fontsize=10)
        axi.set_xlabel("x [m]")
        axi.set_ylabel("z [m]")
        axi.set_ylim(-H * 1.1, 3 * a)
        axi.set_aspect("equal")
    save(fig, out, "c05_orbits")
    fig, axs = plt.subplots(figsize=(9, 4.2))
    a = 0.12
    H = 1.0
    lam = 2 * np.pi / k
    X, Z = np.meshgrid(np.linspace(0, lam, 300), np.linspace(-H, 0, 120), indexing="xy")
    f = ch07.wave_fields(X, Z, 0.0, a, k, H)
    axs.contour(X, Z, f["psi"], levels=np.linspace(-np.abs(f["psi"]).max(), np.abs(f["psi"]).max(), 17),
                colors=COLORS["ink"], linewidths=0.8, negative_linestyles="solid")
    sk = (slice(10, None, 20), slice(15, None, 30))
    axs.quiver(X[sk], Z[sk], f["u"][sk], f["w"][sk], color=COLORS["teal"], scale=6)
    wave_surface(axs, X[0], ch07.wave_fields(X[0], 0.0, 0.0, a, k, H)["eta"])
    axs.set_title("streamlines ψ (7.37) at t = 0, H = 1 m: forward under crests (Fig. 7.5 remake)", fontsize=9)
    axs.set_xlabel("x [m]")
    axs.set_ylabel("z [m]")
    save(fig, out, "c05_streamlines")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
