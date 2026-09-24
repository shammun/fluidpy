"""§7.7 (C13, C14): waves on a density interface — two deep fluids (7.95) with the opposing velocities of a vortex sheet
(Fig. 7.24 remake) and the energy (7.96); a layer over a deep fluid: the barotropic (7.111)–(7.112) and baroclinic
(7.113)–(7.114) modes (Fig. 7.27 remake), their dispersion curves with the (7.95) limit, the long-wave speed √(g′H)
with g′ = g(ρ₂ − ρ₁)/ρ₂ (7.116)–(7.117) (and ch04's ρ₁ convention), and the shallow two-layer mode profiles (Fig. 7.28).

Run: ``.venv/Scripts/python.exe scripts/ch07_interface.py --no-show``
Figures → outputs/ch07/c13_c14_interface.png.
"""
from __future__ import annotations

import numpy as np

from ch07_drawings import NAVY, parse_args, save, setup, two_layer_sketch

from fluidpy import ch05_vorticity_dynamics as ch05
from fluidpy import ch07_gravity_waves as ch07
from fluidpy.core.similarity import reduced_gravity
from fluidpy.core.style import COLORS


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    g = ch07.G_BOOK
    r1, r2 = 1000.0, 1020.0
    k, a = 1.0, 0.1
    f = ch07.interface_fields(0.0, 0.0, 0.0, a, k, r1, r2)
    print(f"two deep fluids ρ₁ = {r1:g}, ρ₂ = {r2:g}: ε² = {ch07.eps2_density(r1, r2):.5f}, ω = {f['omega']:.4f} rad/s vs"
          f" √(gk) = {np.sqrt(g * k):.4f}; A = {f['A']:.4f}, B = {f['B']:.4f} (A = −B = iωa/k)")
    print(f"  sheet strength at a crest γ = u₂ − u₁ = {f['gamma_sheet']:.4f} m/s; ch05.vortex_sheet_strength = "
          f"{ch05.vortex_sheet_strength(f['u1'], f['u2']):.4f} m/s")
    q = ch07.interface_energy(a, k, r1, r2, method="quad")
    print(f"  energy per unit area: quad E_k = {q['Ek']:.6f}, E_p = {q['Ep']:.6f}, E = {q['E']:.6f} J/m² vs"
          f" ½(ρ₂ − ρ₁)ga² = {ch07.wave_energy_density(a, g=g, drho=r2 - r1):.6f}")
    ro1, ro2 = 1025.0, 1027.0
    print(f"reduced gravity ocean ρ₁ = 1025, ρ₂ = 1027: (7.117) {ch07.reduced_gravity_book(ro1, ro2):.5f} m/s², ch04 (ρ₁)"
          f" {reduced_gravity(ro1, ro2, g):.5f} m/s²; oil 800 over water 1000: {ch07.reduced_gravity_book(800, 1000):.3f}"
          f" vs {reduced_gravity(800, 1000, g):.3f} m/s²")
    H = 50.0
    for mode in ("barotropic", "baroclinic"):
        m = ch07.two_layer_modes(2 * np.pi / 5000.0, H, 1025.0, 1027.0, mode=mode)
        print(f"thermocline H = 50 m, λ = 5 km, {mode}: c = {m['c']:.3f} m/s, η/ζ = {m['eta_over_zeta']:.5f}")
    lw = ch07.two_layer_modes(2 * np.pi / 50000.0, H, 1025.0, 1027.0, long_wave=True)
    print(f"  long wave λ = 50 km: c = √(g′H) = {lw['c']:.4f} m/s (exact {lw['omega_exact'] / (2 * np.pi / 50000):.4f});"
          f" η/ζ (7.118) = {lw['eta_over_zeta']:.5f} (exact {lw['eta_over_zeta_exact']:.5f}); p′(−H)/p′(0) = "
          f"{lw['p_prime_check']['interface_over_top']:.5f}")
    res = ch07.two_layer_residuals(np.linspace(0, 5000, 7), 30.0, 2 * np.pi / 5000, H, 1025.0, 1027.0)
    bad = ch07.two_layer_residuals(np.linspace(0, 5000, 7), 30.0, 2 * np.pi / 5000, H, 1025.0, 1027.0,
                                   printed_7_105=True)
    print(f"  residuals of (7.97)–(7.101): max {max(np.max(np.abs(v)) for v in res.values()):.1e}; with (7.105) as printed:"
          f" Laplace {np.max(np.abs(bad['laplace2'])):.1e}, interface kinematic {np.max(np.abs(bad['kin_interface_2'])):.1e}")

    fig, ax = plt.subplots(2, 2, figsize=(13, 8.5))
    lam = 2 * np.pi / k
    X, Z = np.meshgrid(np.linspace(0, lam, 13), np.concatenate([np.linspace(-1.2, -0.15, 5),
                                                                  np.linspace(0.15, 1.2, 5)]), indexing="xy")
    ff = ch07.interface_fields(X, Z, 0.0, a, k, r1, r2)
    xs = np.linspace(0, lam, 200)
    a0 = ax[0, 0]
    a0.fill_between(xs, -1.4, ch07.interface_fields(xs, 0, 0, a, k, r1, r2)["zeta_interface"], color=NAVY, alpha=0.25)
    a0.plot(xs, ch07.interface_fields(xs, 0, 0, a, k, r1, r2)["zeta_interface"], color=NAVY, lw=2)
    a0.quiver(X, Z, ff["u"], ff["w"], color=COLORS["teal"], angles="xy", scale_units="xy", scale=0.1)
    a0.set_ylim(-1.4, 1.4)
    a0.set_title("two deep fluids: u reverses across the interface — a vortex sheet (Fig. 7.24)", fontsize=9)
    a0.set_xlabel("x [m]")
    a0.set_ylabel("z [m]")
    a1 = ax[0, 1]
    kk = np.logspace(-4, -1, 200)
    wbt, wbc = ch07.two_layer_free_surface_omega(kk, H, 1025.0, 1027.0)
    a1.loglog(kk, wbt, color=COLORS["blue"], label="barotropic ω² = gk (7.111)")
    a1.loglog(kk, wbc, color=NAVY, label="baroclinic (7.113), H = 50 m")
    a1.loglog(kk, ch07.interface_omega(kk, 1025.0, 1027.0), color=COLORS["muted"], ls="--", lw=1,
              label="two deep fluids ε√(gk) (7.95)")
    a1.loglog(kk, kk * ch07.two_layer_long_wave_speed(H, 1025.0, 1027.0), color=COLORS["muted"], ls=":", lw=1,
              label="long waves k√(g′H) (7.116)")
    a1.set_xlabel("k [rad/m]")
    a1.set_ylabel("ω [rad/s]")
    a1.set_title("layer over a deep fluid: two modes (ρ₁ = 1025, ρ₂ = 1027)", fontsize=10)
    a1.legend(fontsize=8)
    k2 = 2 * np.pi / 400.0
    xx = np.linspace(0, 400.0, 300)
    for axi, mode in ((ax[1, 0], "barotropic"), (ax[1, 1], "baroclinic")):
        m = ch07.two_layer_modes(k2, H, 1025.0, 1027.0, a=1.0, mode=mode)
        r = m["eta_over_zeta"]
        ph = np.cos(k2 * xx)
        if mode == "barotropic":  # surface 5 m, interface η/r (in phase), both true size
            eta_d, zeta_d, note = 5.0 * ph, 5.0 / r * ph, "true size"
        else:  # interface 10 m; the surface signal is tiny, drawn ×500
            eta_d, zeta_d, note = 500.0 * r * 10.0 * ph, 10.0 * ph, "surface drawn ×500"
        two_layer_sketch(axi, H, 400.0, eta=eta_d, zeta=zeta_d, x=xx)
        axi.set_ylim(-2.2 * H, 15)
        axi.set_xlabel("x [m]")
        axi.set_ylabel("z [m]")
        axi.set_title(f"{mode}: η/ζ = {r:+.5f} " + ("(in phase, (7.112))" if mode == "barotropic"
                                                   else "(antiphase, (7.114))") + f", λ = 400 m; {note}", fontsize=9)
    save(fig, out, "c13_c14_interface")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
