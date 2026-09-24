"""§6.6 (C11): conformal mapping — small elements keep their angles (6.63)–(6.64) except at critical points; the
(φ, ψ) grid of w = z² in the physical plane (Fig. 6.20); w = ln(sin z) by the chain rule; the Zhukhovsky map (6.65)
taking circles to a slit and to ellipses (6.66)–(6.67); the flow round the elliptic cylinder with circulation (6.68)
through the outside root of (6.69), and what numpy's principal root does wrong.

Run: ``.venv/Scripts/python.exe scripts/ch06_conformal_grid.py --no-show``
Figures → outputs/ch06/c11_conformal.png.
"""
from __future__ import annotations

import numpy as np

from ch06_drawings import parse_args, save, setup

from fluidpy import ch06_ideal_flow as ch06
from fluidpy.core.style import COLORS


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    for f in ("square", "exp", "sin", "joukowski"):
        al, be = ch06.angle_preservation(f, z0=0.7 + 0.4j, dz1=1.0, dz2=np.exp(1.1j))
        print(f"  {f:<9s} α = {np.degrees(al):.6f}°, β = {np.degrees(be):.6f}°")
    al, be = ch06.angle_preservation("square", z0=0.0, dz1=1.0, dz2=np.exp(0.5j))
    print(f"  z² at the critical point 0: α = {np.degrees(al):.3f}°, β = {np.degrees(be):.3f}° (doubled)")
    cf = ch06.cot_flow(0.3 + 0.2j)
    print(f"w = ln(sin z): chain rule {cf['dwdz_chain']:.10f}, direct {cf['dwdz_direct']:.10f}")
    b, a = 1.0, 1.2
    E = ch06.joukowski_ellipse(a, b)
    th = np.linspace(0, 2 * np.pi, 400, endpoint=False)
    zc = ch06.joukowski(a * np.exp(1j * th), b)
    print(f"Zhukhovsky b = {b}, a = {a}: semi-axes {E['A']:.6f}, {E['B']:.6f}, foci ±{E['foci']:.6f}; "
          f"max residual of (6.67) on the image {np.max(np.abs(zc.real ** 2 / E['A'] ** 2 + zc.imag ** 2 / E['B'] ** 2 - 1)):.1e}")
    zs = ch06.joukowski(b * np.exp(1j * th), b)
    print(f"  circle |ζ| = b → slit: max |Im z| {np.max(np.abs(zs.imag)):.1e}, Re z ∈ [{zs.real.min():.6f}, {zs.real.max():.6f}]")
    z = -3.0 + 0.5j
    print(f"  inverse at z = {z}: outside root {ch06.joukowski_inverse(z, b):.6f} (|ζ| = {abs(ch06.joukowski_inverse(z, b)):.4f}),"
          f" numpy principal {ch06.joukowski_inverse(z, b, 'principal'):.6f} (|ζ| = "
          f"{abs(ch06.joukowski_inverse(z, b, 'principal')):.4f} < b: inside ⚠)")
    rng = np.random.default_rng(0)
    zeta = (b * 1.001 + 3 * rng.random(10000)) * np.exp(2j * np.pi * rng.random(10000))
    rt = ch06.joukowski_inverse(ch06.joukowski(zeta, b), b)
    rtp = ch06.joukowski_inverse(ch06.joukowski(zeta, b), b, "principal")
    print(f"  round trip on 10⁴ random |ζ| > b (seed 0): outside max error {np.max(np.abs(rt - zeta)):.1e}; "
          f"principal fails at {np.mean(np.abs(rtp - zeta) > 1e-8) * 100:.1f} % of the points")
    G = 1.5
    ef = ch06.elliptic_cylinder_flow(1.0, a, b, Gamma_cw=G)
    ell = E["A"] * np.cos(th) + 1j * E["B"] * np.sin(th)
    print(f"elliptic cylinder Γ_cw = {G}: max |u·n| on the ellipse {ch06.normal_velocity_on(ef, ell):.1e}, "
          f"far field at R = 200: {ch06.far_field_check(ef, 200.0):.1e}, Blasius {ch06.blasius_force(ef, R=3.0)}")

    fig, ax = plt.subplots(1, 3, figsize=(17, 5))
    n = 161 if args.fast else 241
    gi = ch06.grid_image(lambda q: q ** 2, (0, 2), (0, 2), n)
    ax[0].contour(gi["X"], gi["Y"], gi["psi"], levels=np.linspace(0.25, 7.75, 16), colors=COLORS["ink"], linewidths=0.9)
    ax[0].contour(gi["X"], gi["Y"], gi["phi"], levels=np.linspace(-3.75, 3.75, 16), colors=COLORS["teal"],
                  linewidths=0.7, linestyles="--")
    ax[0].set_aspect("equal")
    ax[0].set_title("w = z²: the (φ, ψ) grid in the z-plane (90° corner)", fontsize=10)
    ax[0].set_xlabel("x [m]")
    ax[0].set_ylabel("y [m]")
    for r_, col in ((b, COLORS["rose"]), (1.2, COLORS["accent"]), (1.6, COLORS["teal"]), (2.4, COLORS["blue"])):
        zz = ch06.joukowski(r_ * np.exp(1j * th), b)
        ax[1].plot(np.append(zz.real, zz.real[0]), np.append(zz.imag, zz.imag[0]), color=col, label=f"|ζ| = {r_:g}")
    ax[1].plot([-2 * b, 2 * b], [0, 0], "x", color=COLORS["ink"])
    ax[1].set_aspect("equal")
    ax[1].legend(fontsize=8)
    ax[1].set_title("z = ζ + b²/ζ: circles → slit and ellipses (foci ±2b)", fontsize=10)
    ax[1].set_xlabel("x [m]")
    xg = np.linspace(-4, 4, n)
    yg = np.linspace(-3, 3, n)
    X, Y = np.meshgrid(xg, yg, indexing="xy")
    P = np.asarray(ef.psi(X, Y))
    ax[2].contour(X, Y, P, levels=np.linspace(-3, 3, 31), colors=COLORS["ink"], linewidths=0.8,
                  negative_linestyles="solid")
    ax[2].fill(ell.real, ell.imag, color="#d9dce4")
    st = ef.stagnation
    ax[2].plot(st.real, st.imag, "o", color=COLORS["orange"])
    ax[2].set_aspect("equal")
    ax[2].set_title(f"elliptic cylinder, Γ_cw = {G} m²/s, via (6.68)–(6.69)", fontsize=10)
    ax[2].set_xlabel("x [m]")
    save(fig, out, "c11_conformal")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
