"""§10.5 (C14, C15): the lid-driven cavity by the MAC projection scheme (16² … 128²) and by explicit MacCormack on the weakly
compressible equations (our Ma = 0.08), against Ghia, Ghia & Shin (1982) Table I and the vortex centres of Ghia (1982) and
Hou et al. (1995) (public data, reference/ch10/); grid convergence with Richardson extrapolation and the GCI.

Writes our cached results to reference/ch10/cavity_mac_re{Re}_n{n}.csv and cavity_mck_re{Re}_n{n}.csv (centreline u(y);
run metadata in the header) — regenerate by re-running this script.  128² and Re = 400 are script-only (minutes).

Run: ``.venv/Scripts/python.exe scripts/ch10_cavity.py --no-show [--fast]``   Figure -> outputs/ch10/c14_cavity.png.
"""
from __future__ import annotations

import numpy as np

from ch10_common import COLORS, REF, Timer, finish, parse_args, save, setup, write_csv

from fluidpy import ch10_computational_fluid_dynamics as ch10


def _mac(Re, n, t_end):
    with Timer(f"MAC Re = {Re:g}, {n}²"):
        st = ch10.cavity(Re=Re, n=n, t_end=t_end, tol_steady=1e-6)
    cl = ch10.cavity_centreline(st)
    err = ch10.cavity_error_vs_ghia(st, int(Re))
    c = st["centre"]
    print(f"   steps {st['steps']}, t = {st['t']:.1f}, steady {st['converged']} (residual {st['residual']:.1e}); "
          f"max dev from Ghia {err['max_dev']:.4f}; centre ({c['x']:.4f}, {c['y']:.4f}), ψ_min = {c['psi_min']:.5f}")
    write_csv(REF / f"cavity_mac_re{int(Re)}_n{n}.csv",
              [f"our MAC projection cavity (scripts/ch10_cavity.py): Re = {Re:g}, {n}x{n} cells, dt = {st['dt']:.4g}, "
               f"t = {st['t']:.3f}, steady = {st['converged']}",
               f"primary vortex centre x = {c['x']:.5f}, y = {c['y']:.5f}, psi_min = {c['psi_min']:.6f}; "
               f"max deviation from Ghia 1982 = {err['max_dev']:.5f}"],
              dict(y=cl["y"], u=cl["u"]))
    return st, cl, err


def _mck(Re, n, t_end, wall_density="continuity"):
    with Timer(f"MacCormack Re = {Re:g}, {n}², Ma = 0.08, {wall_density}"):
        r = ch10.cavity_maccormack(Re=Re, Ma=0.08, n=n, t_end=t_end, tol_steady=1e-6, wall_density=wall_density)
    err = ch10.cavity_error_vs_ghia(r, int(Re))
    print(f"   steps {r['steps']}, t = {r['t']:.1f}, steady {r['converged']} (residual {r['residual']:.1e}), "
          f"mass drift {r['mass_drift']:.2e}; max dev from Ghia {err['max_dev']:.4f}")
    if wall_density == "continuity":
        write_csv(REF / f"cavity_mck_re{int(Re)}_n{n}.csv",
                  [f"our weakly compressible MacCormack cavity (scripts/ch10_cavity.py): Re = {Re:g}, Ma = 0.08, {n}x{n} cells, "
                   f"wall density from continuity (10.139)-(10.146), dt = {r['dt']:.4g}, t = {r['t']:.3f}",
                   f"max deviation from Ghia 1982 = {err['max_dev']:.5f}; relative mass drift {r['mass_drift']:.3e}"],
                  dict(y=r["y"], u=r["u"][:, n // 2]))
    return r, err


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    fast = args.fast
    ns = (16, 32, 64) if fast else (16, 24, 32, 64, 128)
    runs = {n: _mac(100.0, n, 40.0) for n in ns}
    runs400 = {n: _mac(400.0, n, 60.0) for n in ((32, 64) if fast else (32, 64, 128))}
    print("Ghia (1982) centre Re = 100:", ch10.ghia_vortex_centre(100), " Hou et al. (1995):", ch10.hou_centres())
    # grid convergence of the vortex strength ψ_min and of u(½, ½) (refinement ratio 2)
    trio = (32, 64, 128) if not fast else (16, 32, 64)
    psi = [runs[n][0]["centre"]["psi_min"] for n in trio]
    gci = ch10.grid_convergence_index(psi[2], psi[1], psi[0])
    print(f"ψ_min on {trio}: {[round(p, 6) for p in psi]} → observed order {gci['p']:.2f}, Richardson {gci['f_ext']:.6f}, "
          f"GCI_fine {100 * gci['gci_fine']:.3f} %")
    mck = {}
    for n in ((32,) if fast else (32, 64)):
        mck[n] = _mck(100.0, n, 30.0)
    if not fast:
        _mck(100.0, 64, 30.0, "momentum")
        _mck(400.0, 64, 50.0)
    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    nbig = max(ns)
    st = runs[nbig][0]
    g = st["g"]
    X, Y = np.meshgrid(np.arange(g.nx + 1) * g.dx, np.arange(g.ny + 1) * g.dy)
    lev = np.concatenate([np.linspace(st["psi"].min(), -1e-4, 12), [1e-6, 1e-5, 1e-4]])
    ax[0].contour(X, Y, st["psi"], levels=np.sort(lev), colors=[COLORS["accent"]], linewidths=1)
    c = st["centre"]
    ax[0].plot(c["x"], c["y"], "o", color=COLORS["accent"], label="ours")
    gx = ch10.ghia_vortex_centre(100)
    ax[0].plot(gx["x"], gx["y"], "ko", mfc="white", label="Ghia 1982")
    ax[0].plot(*ch10.hou_centres()[100], "k^", mfc="white", label="Hou et al. 1995")
    ax[0].set_aspect("equal")
    ax[0].set_title(f"streamlines, MAC {nbig}², Re = 100")
    ax[0].legend(fontsize=8)
    gh = ch10.ghia_centreline(100)
    for n, col in zip(ns, (COLORS["muted"], COLORS["amber"], COLORS["teal"], COLORS["blue"], COLORS["accent"])):
        cl = runs[n][1]
        ax[1].plot(cl["u"], cl["y"], color=col, label=f"MAC {n}²")
    for n, (r, _) in mck.items():
        ax[1].plot(r["u"][:, n // 2], r["y"], "--", color=COLORS["rose"], label=f"MacCormack {n}², Ma = 0.08")
    ax[1].plot(gh["u"], gh["y"], "ko", mfc="white", label="Ghia 1982")
    ax[1].set_xlabel("u / U on x = L/2")
    ax[1].set_ylabel("y / L")
    ax[1].set_title("centreline, Re = 100")
    ax[1].legend(fontsize=7)
    gh4 = ch10.ghia_centreline(400)
    for n, (s4, cl, _) in runs400.items():
        ax[2].plot(cl["u"], cl["y"], label=f"MAC {n}²")
    ax[2].plot(gh4["u"], gh4["y"], "ko", mfc="white", label="Ghia 1982")
    ax[2].set_xlabel("u / U on x = L/2")
    ax[2].set_title("centreline, Re = 400")
    ax[2].legend(fontsize=7)
    save(fig, out, "c14_cavity")
    hs = [1 / n for n in ns]
    errs = [runs[n][2]["max_dev"] for n in ns]
    print("max deviation from Ghia vs h:", dict(zip(ns, [round(e, 5) for e in errs])),
          f"→ slope {ch10.observed_order(hs[:-1] if len(hs) > 3 else hs, errs[:-1] if len(errs) > 3 else errs):.2f} "
          "(the finest grids approach Ghia's own discretisation error)")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
