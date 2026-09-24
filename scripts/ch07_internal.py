"""§7.8 (C15, C16): internal waves in a fluid of constant N — ω = N cos θ (7.139) (frequency set by direction only),
K·u = 0 (7.141), c ⟂ c_g with opposite vertical components (7.144)–(7.146) including the book's k < 0 geometry
(Fig. 7.29/7.31 remake), a packet sliding along its crests by 2-D FFT evolution (Fig. 7.32 remake), the St Andrew's
cross (our superposition illustration of Fig. 7.33), equipartition and F = c_g E (7.154)–(7.159), and the interface
limit of the potential energy (7.150)–(7.152).

Run: ``.venv/Scripts/python.exe scripts/ch07_internal.py --no-show [--fast]``
Figures → outputs/ch07/c15_c16_internal.png.
"""
from __future__ import annotations

import numpy as np

from ch07_drawings import beam_labels, parse_args, save, setup

from fluidpy import ch07_gravity_waves as ch07
from fluidpy.core.style import COLORS


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    N = 1.0
    st = ch07.internal_wave_state(0.71, N, 1.0, k_sign=-1.0)
    print(f"ω/N = 0.71: θ_K = {st['theta_K_deg']:.2f}° above the horizontal = beam angle from the vertical; "
          f"c = ({st['cx']:+.4f}, {st['cz']:+.4f}), c_g = ({st['cgx']:+.4f}, {st['cgz']:+.4f}) m/s, c·c_g = {st['dot']:.1e}")
    printed = ch07.internal_wave_velocities(st["k"], st["m"], N, printed=True)["cg"]
    print(f"  k < 0 (Fig. 7.29 geometry): the printed (7.145) gives c_g = ({printed[0]:+.4f}, {printed[1]:+.4f}) — up-right,"
          f" not down-left; ∇_K ω gives ({st['cgx']:+.4f}, {st['cgz']:+.4f})")
    gv = ch07.group_velocity_vector(lambda K: ch07.internal_wave_omega(K[0], K[1], N), [st["k"], st["m"]])
    print(f"  complex-step ∇_K ω = ({gv[0]:+.6f}, {gv[1]:+.6f})")
    assert abs(st["dot"]) < 1e-14 and np.allclose(gv, [st["cgx"], st["cgz"]], atol=1e-12)
    for th_deg in (0.0, 30.0, 60.0):
        th = np.radians(th_deg)
        for K in (0.5, 2.0):
            w = ch07.internal_wave_omega(K * np.cos(th), K * np.sin(th), N)
            print(f"  θ = {th_deg:>4}°, |K| = {K}: ω = {w:.6f} rad/s (N cos θ = {N * np.cos(th):.6f})")
    f = ch07.internal_wave_fields(0.3, 0.2, 0.1, -1.0, 2.0, N, 0.01, residuals=True)
    print(f"plane wave k = −1, m = 2: K·u = {f['K_dot_u']:.1e}; FD residuals " +
          ", ".join(f"{k_} {v:.1e}" for k_, v in f["residuals"].items()))
    e = ch07.internal_wave_energy(-1.0, 2.0, N, 0.01)
    print(f"  E_k = {e['Ek']:.6e}, E_p = {e['Ep']:.6e} J/m³; F = {e['F']} W/m², c_g E = {e['cgE']}")
    assert np.allclose(e["F"], e["cgE"], rtol=1e-12) and abs(e["Ek"] / e["Ep"] - 1) < 1e-12
    for eps in (4.0, 2.0, 1.0, 0.5):
        d = ch07.internal_pe_interface_limit(1.0, 1000.0, 1020.0, eps=eps, k=0.05, detail=True)
        print(f"  δ-function limit (7.152), kε = {0.05 * eps:.3f}: ∫E_p dz = {d['Ep_column']:.4f} J/m² vs ¼Δρga² = "
              f"{d['target']:.4f} (rel {d['rel_error']:+.4f})")

    fig, ax = plt.subplots(1, 3, figsize=(17, 5.4))
    a0 = ax[0]
    s = 1.0 / max(st["c"], st["cg"])
    a0.annotate("", xy=(st["cx"] * s, st["cz"] * s), xytext=(0, 0), arrowprops=dict(arrowstyle="->", color=COLORS["orange"], lw=2.5))
    a0.annotate("", xy=(st["cgx"] * s, st["cgz"] * s), xytext=(0, 0), arrowprops=dict(arrowstyle="->", color=COLORS["accent"], lw=2.5))
    a0.annotate("", xy=(printed[0] * s, printed[1] * s), xytext=(0, 0),
                arrowprops=dict(arrowstyle="->", color=COLORS["rose"], lw=1.2, ls="--"))
    a0.text(st["cx"] * s - 0.25, st["cz"] * s + 0.05, "K and c", color=COLORS["orange"])
    a0.text(st["cgx"] * s - 0.2, st["cgz"] * s - 0.12, "c_g = ∇_K ω", color=COLORS["accent"])
    a0.text(printed[0] * s + 0.02, printed[1] * s, "printed (7.145)\nfor k < 0", color=COLORS["rose"], fontsize=8)
    tt = np.linspace(-1.2, 1.2, 2)
    kd = np.array([st["k"], st["m"]]) / np.hypot(st["k"], st["m"])
    for off in (-0.6, -0.3, 0.0, 0.3, 0.6):
        a0.plot(off * kd[0] + tt * (-kd[1]), off * kd[1] + tt * kd[0], color=COLORS["muted"], lw=0.6)
    a0.set_xlim(-1.2, 1.2)
    a0.set_ylim(-1.2, 1.2)
    a0.set_aspect("equal")
    a0.set_title(f"ω/N = 0.71, K up-left (Fig. 7.29): c ⊥ c_g, θ = {st['theta_K_deg']:.1f}°", fontsize=10)
    a0.set_xlabel("x")
    a0.set_ylabel("z")
    # packet
    n = 128 if args.fast else 256
    L = 200.0
    x = np.linspace(-L / 2, L / 2, n, endpoint=False)
    z = np.linspace(-L / 2, L / 2, n, endpoint=False)
    X, Z = np.meshgrid(x, z, indexing="xy")
    k0, m0 = st["k"], st["m"]
    field0 = np.exp(-(X ** 2 + Z ** 2) / (2 * 12.0 ** 2)) * np.exp(1j * (k0 * X + m0 * Z))
    tp = 40.0
    ft = ch07.linear_evolve_2d(field0, x, z, [0.0, tp], lambda kk, mm: ch07.internal_wave_omega(kk, mm, N))
    cen = [np.array([np.sum(X * np.abs(fr) ** 2), np.sum(Z * np.abs(fr) ** 2)]) / np.sum(np.abs(fr) ** 2) for fr in ft]
    v = (cen[1] - cen[0]) / tp
    print(f"packet (σ = 12 m, K = (−0.71, 0.70) rad/m): centroid velocity ({v[0]:+.4f}, {v[1]:+.4f}) m/s vs c_g "
          f"({st['cgx']:+.4f}, {st['cgz']:+.4f}) m/s")
    a1 = ax[1]
    for fr, alpha in zip(ft, (0.45, 1.0)):
        a1.contour(X, Z, fr, levels=[-0.5, 0.5], colors=[COLORS["blue"], COLORS["rose"]], linewidths=0.8, alpha=alpha)
    a1.annotate("", xy=cen[1], xytext=cen[0], arrowprops=dict(arrowstyle="->", color=COLORS["accent"], lw=2))
    a1.set_xlim(-60, 40)
    a1.set_ylim(-60, 40)
    a1.set_aspect("equal")
    a1.set_title(f"packet at t = 0 (faint) and {tp:.0f} s: phase up-left, energy down-left (Fig. 7.32)", fontsize=9)
    a1.set_xlabel("x [m]")
    a1.set_ylabel("z [m]")
    # St Andrew's cross
    xs = np.linspace(-10, 10, 301)
    fld = ch07.st_andrews_cross(xs, xs, 0.0, 0.71, N, width=1.0)
    a2 = ax[2]
    a2.imshow(fld, extent=(-10, 10, -10, 10), origin="lower", cmap="gray")
    beam_labels(a2, ch07.beam_angle(0.71, N), 9.0)
    a2.set_title("St Andrew's cross, ω/N = 0.71 (our superposition of four beams, cf. Fig. 7.33)", fontsize=9)
    a2.set_xlabel("x [m]")
    a2.set_ylabel("z [m]")
    save(fig, out, "c15_c16_internal")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
