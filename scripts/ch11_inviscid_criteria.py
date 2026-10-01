"""§11.9 (C10, C12): Rayleigh's equation (11.81) and its criteria — the six profiles of Fig. 11.21 through Rayleigh (11.84) and
Fjørtoft (11.86), the tanh shear layer (most-amplified mode, neutral sech mode), the Bickley jet's sinuous/varicose neutral
modes, U = sin y between walls (stable for 2b < π), the piecewise-linear layer of Ex. 11.11, the identities (11.83)–(11.85),
critical layers and the Kelvin cat's eye (Fig. 11.22 analogue).

Run: ``.venv/Scripts/python.exe scripts/ch11_inviscid_criteria.py --no-show [--fast]``   Figure -> outputs/ch11/c12_inviscid.png.
"""
from __future__ import annotations

import math

import numpy as np

from ch11_common import COLORS, Timer, finish, parse_args, save, setup

from fluidpy import ch11_instability as ch11


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    for d in ch11.fig_11_21_verdicts():
        print(f"Fig. 11.21 ({d['panel']}) {d['label']:40s} Rayleigh {d['rayleigh']!s:5s} Fjørtoft {d['fjortoft']}")
    with Timer("tanh max growth"):
        tm = ch11.tanh_max_growth()
    print(f"tanh: k_max = {tm['k']:.5f}, (kc_i)_max = {tm['kci']:.5f} (Michalke 1964: 0.4446, 0.1897); c = {tm['c']:.3e}")
    pr = ch11.parallel_profile("tanh")
    r = ch11.rayleigh_eigs(0.44, pr["U"], pr["Upp"], bc="decay", N=120, return_vectors=True, unstable_only=True)
    ri = ch11.rayleigh_identity_check(0.44, r["c"][0], r["phi"][:, 0], r["grid"].y, pr["U"], pr["Upp"], grid=r["grid"])
    print(f"(11.83)/(11.84) residuals for the tanh mode at k = 0.44: {ri}")
    for k in (0.9, 0.97, 0.99):
        c = ch11.rayleigh_eigs(k, pr["U"], pr["Upp"], bc="decay", N=120, map_scale=0.5, unstable_only=True)
        print(f"  tanh k = {k}: c_i = {c[0].imag if len(c) else 0.0:.5f}  (→ 0 at the neutral k = 1, φ = sech y)")
    bj = ch11.parallel_profile("bickley")
    for k, par in ((1.5, "even"), (1.9, "even"), (0.5, "odd"), (0.9, "odd")):
        c = ch11.rayleigh_eigs(k, bj["U"], bj["Upp"], bc="decay", N=120, parity=par, unstable_only=True)
        print(f"  Bickley {('sinuous' if par == 'even' else 'varicose'):9s} k = {k}: c = {c[0] if len(c) else None}"
              f"  (neutral at k = {2 if par == 'even' else 1}, c = 2/3)")
    for b in (1.4, 1.6, 1.7, 2.0, 3.0):
        print(f"  U = sin y, |y| ≤ {b}: 2b = {2 * b:.2f} (π = {math.pi:.4f}), max kc_i = {ch11.sin_profile_max_growth(b):.5f}")
    print(f"Ex. 11.11 piecewise layer: neutral kh = {ch11.piecewise_neutral_kh():.6f}; c(kh = 0.797) = "
          f"{ch11.piecewise_shear_layer_c(0.797):.5f}")
    print(f"critical layer of Poiseuille at c = 0.5: y_c = {ch11.critical_layer(np.linspace(-1, 1, 41), lambda y: 1 - y ** 2, 0.5)}")
    print(f"cat's-eye half-width A = 0.1, φ_c = 1, U′ = 1: {ch11.cats_eye_width(0.1, 1.0, 1.0):.4f}")
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
    ks = np.linspace(0.02, 1.0, 25 if args.fast else 50)
    gro = [k * (lambda c: c[0].imag if len(c) else 0.0)(ch11.rayleigh_eigs(k, pr["U"], pr["Upp"], bc="decay", N=100,
                                                                           unstable_only=True)) for k in ks]
    ax[0].plot(ks, gro, color=COLORS["rose"], label="tanh y (Rayleigh)")
    kh = np.linspace(0.02, 1.4, 200)
    ax[0].plot(kh, kh * np.imag(ch11.piecewise_shear_layer_c(kh, 2.0)) / 2, color=COLORS["muted"], ls="--",
               label="piecewise layer (Ex. 11.11), ΔU = 2, h = 2")
    ax[0].set(xlabel="k L", ylabel="k c_i L/U₀", title="inviscid growth of free shear layers")
    ax[0].legend(fontsize=8)
    cr, ci = ch11.howard_semicircle(-1.0, 1.0)
    ax[1].fill(cr, ci, color=COLORS["accent"], alpha=0.15)
    for k in (0.2, 0.44, 0.7):
        c = ch11.rayleigh_eigs(k, pr["U"], pr["Upp"], bc="decay", N=100, unstable_only=True)
        ax[1].plot(np.real(c), np.imag(c), "o", color=COLORS["rose"])
    for k in (0.5, 1.0, 1.5):
        c = ch11.rayleigh_eigs(k, bj["U"], bj["Upp"], bc="decay", N=100, unstable_only=True)
        ax[1].plot(np.real(c), np.imag(c), "s", color=COLORS["orange"])
    ax[1].plot(*ch11.howard_semicircle(0.0, 1.0), color=COLORS["orange"], ls=":")
    ax[1].set(xlabel="c_r", ylabel="c_i", title="tanh (●) and Bickley (■) eigenvalues in the semicircles")
    ax[1].set_aspect("equal")
    X, Y = np.meshgrid(np.linspace(0, 4 * math.pi, 300), np.linspace(-1.2, 1.2, 200))
    P = ch11.cats_eye_streamfunction(X, Y, A=0.1, phi_c=1.0, k=1.0, Uy_c=1.0)
    ax[2].contour(X, Y, P, 24, colors=COLORS["accent"], linewidths=0.8)
    ax[2].contour(X, Y, P, [0.1], colors=COLORS["rose"], linewidths=1.5)
    ax[2].set(xlabel="x", ylabel="y − y_c", title="Kelvin cat's eye (p. 514 expansion), separatrix in rose")
    save(fig, out, "c12_inviscid")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
