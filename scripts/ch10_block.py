"""§10.5 (C15, R11, N86–N91): flow past a square block between sliding plates by explicit MacCormack (our geometry: H = 4,
8 sides ahead, 20 behind; Ma = 0.06, σ = 0.8): drag and lift histories at Re = 20 (steady after an acoustic transient) and
Re = 100 (periodic shedding), and the grid-convergence study of C_D at Re = 20 (Δx = 1/8, 1/16, 1/32) with Richardson and GCI.

Writes our force histories to reference/ch10/block_forces_re{Re}_dx{1/dx}.csv (cached; regenerate by re-running).
The finest grid and Re = 100 take tens of minutes (script-only).

Run: ``.venv/Scripts/python.exe scripts/ch10_block.py --no-show [--fast] [--convergence]``
Figure -> outputs/ch10/c15_block.png.
"""
from __future__ import annotations

import numpy as np

from ch10_common import COLORS, REF, Timer, finish, parse_args, save, setup, write_csv

from fluidpy import ch10_computational_fluid_dynamics as ch10


def _run(Re, dx, t_end, kick=0.0):
    with Timer(f"block Re = {Re:g}, dx = 1/{round(1 / dx)}, t_end = {t_end:g}"):
        r = ch10.block_channel(Re=Re, Ma=0.06, dx=dx, t_end=t_end, kick=kick, verbose=False)
    t, CD, CL = r["t"], r["CD"], r["CL"]
    late = t > t_end - 10.0
    print(f"   {r['nx']}x{r['ny']} nodes, dt = {r['dt']:.3e}, {r['steps']} steps; last 10 time units: mean C_D = "
          f"{CD[late].mean():.4f} (min {CD[late].min():.4f}, max {CD[late].max():.4f}), mean C_L = {CL[late].mean():+.4f}, "
          f"C_L amplitude {0.5 * (CL[late].max() - CL[late].min()):.4f}")
    stride = max(1, int(round(0.1 / (t[1] - t[0]))))
    write_csv(REF / f"block_forces_re{int(Re)}_dx{round(1 / dx)}.csv",
              [f"our explicit MacCormack block-in-channel run (scripts/ch10_block.py): Re = {Re:g}, Ma = 0.06, sigma = 0.8, "
               f"dx = dy = 1/{round(1 / dx)}, H = 4, 8 ahead, 20 behind, FB/BF-BF/FB cycling, kick = {kick:g}",
               "t in units of D/U; C_D, C_L per unit span normalised by rho U^2 D / 2 (confined configuration: qualitative only)"],
              dict(t=t[::stride], CD=CD[::stride], CL=CL[::stride]))
    return r


def main() -> int:
    args = parse_args(__doc__, lambda ap: ap.add_argument("--convergence", action="store_true"))
    out = setup(args)
    import matplotlib.pyplot as plt

    fast = args.fast
    dxs = (1 / 8, 1 / 16) if fast else (1 / 8, 1 / 16, 1 / 32)
    runs = {dx: _run(20.0, dx, 30.0) for dx in dxs}
    cds = {dx: float(r["CD"][r["t"] > 20.0].mean()) for dx, r in runs.items()}
    print("Re = 20, mean C_D over 20 < t < 30:", {f"1/{round(1 / k)}": round(v, 4) for k, v in cds.items()})
    if len(dxs) == 3:
        f3, f2, f1 = (cds[d] for d in dxs)
        g = ch10.grid_convergence_index(f1, f2, f3)
        print(f"   observed order {g['p']:.2f}, Richardson C_D = {g['f_ext']:.4f}, GCI_fine {100 * g['gci_fine']:.2f} %, "
              f"monotone {g['monotone']}")
        if not (1.5 <= g["p"] <= 2.5):
            print("   -> non-asymptotic grid study (observed order far from the scheme's 2): C_D is NOT grid independent; "
                  "treat the Richardson value and the GCI as indicative only")
    shed = None
    if not fast:
        shed = _run(100.0, 1 / 16, 100.0, kick=2.0)
        late = shed["t"] > 60.0
        f = ch10.dominant_frequency(shed["t"][late], shed["CL"][late])
        f2 = ch10.dominant_frequency(shed["t"][late], shed["CL"][late], "fft")
        print(f"Re = 100 (dx = 1/16): shedding frequency {f:.4f} (zero crossings), {f2:.4f} (FFT) → St = f D/U = {f:.4f} "
              f"(confined block, sliding plates: not comparable with unbounded data)")
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
    for (dx, r), c in zip(runs.items(), (COLORS["orange"], COLORS["teal"], COLORS["blue"])):
        ax[0].plot(r["t"], r["CD"], color=c, lw=1, label=f"C_D, Δx = 1/{round(1 / dx)}")
    ax[0].set_ylim(0, 3 * max(cds.values()))
    ax[0].set_xlabel("t U/D")
    ax[0].set_ylabel("C_D")
    ax[0].set_title("Re = 20: acoustic transient, then steady")
    ax[0].legend(fontsize=8)
    h = np.array(dxs)
    ax[1].plot(h, [cds[d] for d in dxs], "o-", color=COLORS["accent"])
    ax[1].set_xlabel("Δx = Δy (block side 1)")
    ax[1].set_ylabel("mean C_D (20 < t < 30)")
    ax[1].set_title("grid convergence of C_D (Re = 20)")
    if shed is not None:
        ax[2].plot(shed["t"], shed["CD"], color=COLORS["accent"], lw=1, label="C_D")
        ax[2].plot(shed["t"], shed["CL"], color=COLORS["teal"], lw=1, label="C_L")
        ax[2].set_xlabel("t U/D")
        ax[2].set_title("Re = 100: periodic shedding (Δx = 1/16)")
        ax[2].legend(fontsize=8)
    r = runs[dxs[-1]]
    if shed is None:
        ax[2].contourf(r["x"], r["y"], np.where(r["solid"], np.nan, r["u"]), 30, cmap="RdBu_r")
        ax[2].set_aspect("equal")
        ax[2].set_xlim(4, 16)
        ax[2].set_title("u at t = 30 (Re = 20)")
    save(fig, out, "c15_block")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
