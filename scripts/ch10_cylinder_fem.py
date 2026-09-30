"""§10.4–10.5 (C13, N93–N110, R12, R13): mixed P2–P1 finite elements for a cylinder in a channel with sliding walls (our
geometry: W = 5d, x from −8d to 16d), Newton at every step.  --steady: Re = 1, 10, 40 (streamlines, vorticity, C_D, the
residual history of Newton).  --unsteady: Re = 100 from the (unstable) steady state with a smooth rotation pulse, second-order
time derivative (10.163) with α = 2, β = 1: symmetric wake → Hopf instability → Kármán street; force history → St.
Also the code checks: Poiseuille exact, Kovasznay orders, inf–sup constants of P1–P1 vs P2–P1.

Writes our data to reference/ch10/cylinder_fe_steady.csv and (--unsteady) reference/ch10/cylinder_fe_re100_forces.csv.
The unsteady run takes ~30 min (script-only).

Run: ``.venv/Scripts/python.exe scripts/ch10_cylinder_fem.py --no-show [--steady] [--unsteady] [--fast]``
Figures -> outputs/ch10/c13_cylinder_steady.png, c13_cylinder_re100.png.
"""
from __future__ import annotations

import numpy as np

from ch10_common import COLORS, REF, Timer, finish, parse_args, save, setup, write_csv

from fluidpy import ch10_computational_fluid_dynamics as ch10
from fluidpy.core import fem2d as FEM2


def _extra(ap):
    ap.add_argument("--steady", action="store_true")
    ap.add_argument("--unsteady", action="store_true")
    ap.add_argument("--t-end", type=float, default=150.0)
    ap.add_argument("--dt", type=float, default=0.125)


def checks():
    po = FEM2.poiseuille_test(4)
    pp = FEM2.poiseuille_test(4, printed=True)
    print(f"P2–P1 Poiseuille (exact field): max err u {po['err_u']:.1e}, p {po['err_p']:.1e}; printed (10.172) row: err u "
          f"{pp['err_u']:.1e}")
    k = FEM2.kovasznay_test(8)
    print(f"Kovasznay Re = 40, 8² → 4²: L² err u {k['err_u']:.3e}, p {k['err_p']:.3e}; orders {k['orders']}; Newton residuals "
          f"{[f'{r:.1e}' for r in k['residuals']]}")
    for pair in ("P1P1", "P2P1"):
        print(f"inf–sup {pair}:", [(n, round(FEM2.infsup_constant(n, pair)["beta"], 4), FEM2.infsup_constant(n, pair)["n_spurious"])
                                   for n in (2, 4, 8)])


def steady(out, fast):
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri

    level = "coarse" if fast else "medium"
    res = {}
    for Re in ((40.0,) if fast else (1.0, 10.0, 40.0)):
        with Timer(f"steady Re = {Re:g} ({level} mesh)"):
            r = FEM2.cylinder_steady(Re, level)
        res[Re] = r
        rr = r["residuals"]
        ratios = [rr[i + 1] / rr[i] ** 2 for i in range(len(rr) - 1) if rr[i] > 1e-12]
        print(f"   Newton residuals {[f'{x:.1e}' for x in rr]}; r_(k+1)/r_k² {[f'{x:.2g}' for x in ratios]}; C_D = {r['CD']:.4f}, "
              f"C_L = {r['CL']:.1e} (symmetric mesh ⇒ 0)")
    m = next(iter(res.values()))["mesh"]
    print("mesh:", dict(V=m.V, E=m.E, T=m.T), FEM2.euler_check(m))
    write_csv(REF / "cylinder_fe_steady.csv",
              [f"our steady mixed P2-P1 FE cylinder in a channel (scripts/ch10_cylinder_fem.py): W = 5d, x in [-8d, 16d], "
               f"sliding walls, {level} mesh (V = {m.V}, T = {m.T}); confined: qualitative only against unbounded data"],
              dict(Re=list(res), CD=[r["CD"] for r in res.values()], newton_iterations=[r["iterations"] for r in res.values()]))
    fig, ax = plt.subplots(len(res), 2, figsize=(14, 2.4 * len(res) + 0.6), squeeze=False)
    tri = mtri.Triangulation(m.nodes[:m.V, 0], m.nodes[:m.V, 1], m.tris)
    xg, yg = np.meshgrid(np.linspace(-3, 10, 260), np.linspace(-2.45, 2.45, 100))
    for i, (Re, r) in enumerate(res.items()):
        iu = mtri.LinearTriInterpolator(tri, r["u"][:m.V])(xg, yg)
        iv = mtri.LinearTriInterpolator(tri, r["v"][:m.V])(xg, yg)
        ax[i, 0].streamplot(xg, yg, np.ma.filled(iu, np.nan), np.ma.filled(iv, np.nan), density=1.4, color=COLORS["accent"],
                            linewidth=0.7, arrowsize=0.6)
        ax[i, 1].tricontour(tri, r["omega"], levels=np.linspace(-6, 6, 25), cmap="RdBu_r", linewidths=0.8)
        for a in ax[i]:
            a.add_patch(plt.Circle((0, 0), 0.5, color=COLORS["muted"]))
            a.set_xlim(-3, 10)
            a.set_ylim(-2.5, 2.5)
            a.set_aspect("equal")
        ax[i, 0].set_title(f"streamlines, Re = {Re:g} (C_D = {r['CD']:.3f}, confined)")
        ax[i, 1].set_title(f"vorticity, Re = {Re:g}")
    save(fig, out, "c13_cylinder_steady")


def unsteady(out, t_end, dt):
    import matplotlib.pyplot as plt

    mesh = FEM2.cylinder_channel_mesh(n_theta=48, h_far=0.4, n_r=12, h_wake=0.25, wake_length=10.0)
    bc = FEM2.channel_bcs(mesh)
    print("unsteady mesh:", dict(V=mesh.V, E=mesh.E, T=mesh.T), FEM2.p2_node_counts(mesh.V, mesh.E, mesh.T))
    U = None
    with Timer("steady start Re = 10 → 40 → 70 → 100"):
        for Re in (10.0, 40.0, 70.0, 100.0):
            r = FEM2.newton_solve(mesh, None, Re, bc=bc, U0=U, max_iter=20)
            U = (r["u"], r["v"], r["p"])
    print(f"   steady (symmetric, unstable) Re = 100: C_D = {FEM2.cylinder_forces(mesh, r, 100.0)['CD']:.4f}")
    nsteps = int(round(t_end / dt))
    with Timer(f"unsteady Re = 100, dt = {dt}, {nsteps} steps"):
        mu = FEM2.march_unsteady(mesh, bc, 100.0, dt, nsteps, U0=U, alpha_t=2.0, beta_t=1.0, max_iter=4, tol=1e-7,
                                 kick=1.0, kick_time=2.0, start_be=0, verbose=True)
    t, CD, CL, CM = mu["t"], mu["CD"], mu["CL"], mu["CM"]
    write_csv(REF / "cylinder_fe_re100_forces.csv",
              [f"our mixed P2-P1 FE cylinder in a channel, Re = 100 (scripts/ch10_cylinder_fem.py --unsteady): W = 5d, "
               f"x in [-8d, 16d], sliding walls, V = {mesh.V}, T = {mesh.T} triangles, dt = {dt} d/U, alpha_t = 2, beta_t = 1, "
               "started from the steady state with a smooth rotation pulse (0 < t < 2)",
               "t in d/U; C_D, C_L per unit span / (rho U^2 d/2); C_M = torque / (rho U^2 d^2/2); confined: qualitative only"],
              dict(t=t, CD=CD, CL=CL, CM=CM))
    late = t > t_end - 50.0
    if np.ptp(CL[late]) > 1e-3:
        f = ch10.dominant_frequency(t[late], CL[late])
        f2 = ch10.dominant_frequency(t[late], CL[late], "fft")
        print(f"Re = 100: mean C_D {CD[late].mean():.4f}, C_L amplitude {0.5 * np.ptp(CL[late]):.4f}; shedding frequency "
              f"{f:.4f} (zero crossings) / {f2:.4f} (FFT) → St = {ch10.strouhal_from_period(1 / f):.4f} (period "
              f"{1 / f:.3f} d/U; confined, sliding walls — not the unbounded value)")
    else:
        print("Re = 100: no periodic shedding developed within t_end")
    fig, ax = plt.subplots(1, 1, figsize=(10, 4))
    ax.plot(t, CD, color=COLORS["accent"], label="C_D")
    ax.plot(t, CL, color=COLORS["teal"], label="C_L")
    ax.plot(t, CM, color=COLORS["orange"], label="C_M (torque)")
    ax.set_xlabel("t U/d")
    ax.set_title("Re = 100: symmetric wake → Hopf instability → Kármán street (confined)")
    ax.legend()
    save(fig, out, "c13_cylinder_re100")


def main() -> int:
    args = parse_args(__doc__, _extra)
    out = setup(args)
    checks()
    if args.steady or not args.unsteady:
        steady(out, args.fast)
    if args.unsteady:
        unsteady(out, args.t_end, args.dt)
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
