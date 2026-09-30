"""§10.3 (C06–C08): weak form, Galerkin with hat functions, element matrices and assembly, FE = centred FD for the steady
layer, the consistent-mass transient, and the printed slips R1, R2.

Run: ``.venv/Scripts/python.exe scripts/ch10_fem1d.py --no-show``   Figure -> outputs/ch10/c07_fem1d.png.
"""
from __future__ import annotations

import numpy as np

from ch10_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch10_computational_fluid_dynamics as ch10


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    print("sympy: weak form residual", ch10.weak_form_sympy()["residual"], "| weak ⇒ strong",
          ch10.weak_to_strong_sympy()["residual"], "| Galerkin interior row vs (10.63)",
          ch10.galerkin_equations_sympy(3)["check_interior"])
    m, k = ch10.element_matrices_linear(0.25, 1.0, 0.25)
    print("element matrices (h = 0.25, u = 1, D = 0.25): m =", m.round(5).tolist(), " k =", k.tolist())
    m2, k2 = ch10.element_integrals(0.0, 0.25, 1.0, 0.25)
    print(f"  closed form vs 2-point Gauss: {max(np.abs(m - m2).max(), np.abs(k - k2).max()):.1e}")
    print("interior row (10.63) × h:", ch10.interior_stencil(0.25, 1.0, 0.25)["stiff"].tolist())
    print("connectivity (10.78):", ch10.connectivity(4).tolist(), "  as printed (R2):", ch10.connectivity(4, printed=True).tolist())
    nodes = np.linspace(0, 1, 5)
    s, sp_ = ch10.shape_slopes(0.5, 0.75, A=3), ch10.shape_slopes(0.5, 0.75, printed=True, A=3)
    xm, e = 0.625, 1e-6
    num = {A: (ch10.hat(xm + e, nodes, A) - ch10.hat(xm - e, nodes, A)) / (2 * e) for A in range(5)}
    print(f"R1 check on [0.5, 0.75]: correct {dict(zip(s['nodes'], s['slopes']))}, printed {dict(zip(sp_['nodes'], sp_['slopes']))}, "
          f"numerical {{2: {num[2]:.3f}, 3: {num[3]:.3f}, 4: {num[4]:.3f}}}")
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
    xf = np.linspace(0, 1, 801)
    for n, c in ((4, COLORS["orange"]), (8, COLORS["teal"]), (16, COLORS["blue"])):
        fe = ch10.solve_steady(np.linspace(0, 1, n + 1), 1.0, 0.1, T_L=1.0)
        fd = ch10.steady_cd_fd(n, 10.0)
        ax[0].plot(fe["x"], fe["T"], "o", color=c, label=f"FE, n = {n}")
        ax[0].plot(fd["x"], fd["T"], "x", color="k", ms=5)
        print(f"steady layer R = 10, n = {n}: max |FE − centred FD| = {np.abs(fe['T'] - fd['T']).max():.1e}")
    ax[0].plot(xf, ch10.steady_cd_exact(xf, 10.0), color=COLORS["muted"], label="exact (10.86)")
    ax[0].set_title("steady layer: linear FE = centred FD (×)")
    ax[0].set_xlabel("x / L")
    ax[0].set_ylabel("T")
    ax[0].legend(fontsize=8)
    # weighted hats
    n = 6
    nodes = np.linspace(0, 1, n + 1)
    fe = ch10.solve_steady(nodes, 1.0, 0.25, T_L=1.0)
    N = ch10.hat_basis(nodes, xf)
    for A in range(1, n):
        ax[1].plot(xf, fe["T"][A] * N[:, A], color=COLORS["teal"], lw=1)
    ax[1].plot(xf, ch10.interpolate(fe["T"][1:-1], 0.0, nodes[:-1], xf * (xf <= nodes[-2]) + nodes[-2] * (xf > nodes[-2])) * 0
               + N @ fe["T"], color=COLORS["blue"], lw=2, label="$T^h=\\sum d_A N_A$")
    ax[1].plot(xf, ch10.steady_cd_exact(xf, 4.0), color=COLORS["muted"], ls="--", label="exact")
    ax[1].set_title("weighted hats and their sum (R = 4)")
    ax[1].set_xlabel("x / L")
    ax[1].legend(fontsize=8)
    # transient: consistent vs lumped mass, method of lines
    nodes = np.linspace(0, 1, 21)
    T0 = np.zeros(21)
    for lumped, c in ((False, COLORS["blue"]), (True, COLORS["orange"])):
        r = ch10.solve_transport(nodes, T0, 1.0, 0.1, g=1.0, q=0.0, dt=0.005, nsteps=60, theta=0.5, lumped=lumped)
        ax[2].plot(nodes, r["T"], "o-", color=c, ms=3, label=f"{'lumped' if lumped else 'consistent'} mass, t = {r['t']:.2f} s")
    r2 = ch10.solve_transport(nodes, T0, 1.0, 0.1, g=1.0, q=0.0, dt=0.005, nsteps=60, method="solve_ivp")
    ax[2].plot(nodes, r2["T"], "k:", label="method of lines (Radau)")
    print(f"transient FE: θ-scheme vs method of lines at t = 0.3 s: max diff {np.abs(r2['T'] - ch10.solve_transport(nodes, T0, 1.0, 0.1, g=1.0, q=0.0, dt=0.005, nsteps=60)['T']).max():.2e}")
    ax[2].set_title("transient T(0) = 1, ∂T/∂x(L) = 0 (u = 1, D = 0.1)")
    ax[2].set_xlabel("x  [m]")
    ax[2].legend(fontsize=8)
    save(fig, out, "c07_fem1d")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
