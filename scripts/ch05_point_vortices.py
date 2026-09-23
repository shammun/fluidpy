"""§5.7 interaction of vortices (C12, C13, N40–N43, N55–N59): E8's point-vortex presets with their invariants, the
pair rates against (Γ₁ + Γ₂)/2πh² and Γ/2πh, the wall drift Γ/4πh with zero normal velocity on the wall, the channel
image series vs its closed form, the knife-blade pair in a bucket, a from-scratch double loop vs
``point_vortex_velocity``, leap-frogging rings and a ring approaching a wall.

Run: ``.venv/Scripts/python.exe scripts/ch05_point_vortices.py --no-show``
Figures → outputs/ch05/c12_pairs.png, c13_images_rings.png.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from fluidpy import ch05_vorticity_dynamics as ch05  # noqa: E402
from fluidpy.core.style import COLORS  # noqa: E402


def pairs_figure(fast: bool):
    import matplotlib.pyplot as plt
    from ch05_drawings import vortex_marks

    names = ("equal_pair", "unequal_pair", "opposite_pair", "three_vortices")
    fig, ax = plt.subplots(1, 4, figsize=(17, 4.2))
    for a, nm in zip(ax, names):
        pr = ch05.point_vortex_preset(nm)
        tr = ch05.point_vortex_evolve(pr["xv"], pr["Gamma"], np.linspace(0, pr["t_end"], 150 if fast else 300))
        vortex_marks(a, tr[-1], pr["Gamma"], size=0.05, trails=tr)
        c = ch05.centre_of_vorticity(pr["xv"], pr["Gamma"])
        if np.all(np.isfinite(c)):
            a.plot(*c, "x", color=COLORS["amber"], ms=8)
        a.set_title(pr["label"], fontsize=9)
    return fig


def images_rings_figure(fast: bool):
    import matplotlib.pyplot as plt
    from ch05_drawings import vortex_marks

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
    pr = ch05.point_vortex_preset("knife_bucket")
    tr = ch05.point_vortex_evolve(pr["xv"], pr["Gamma"], np.linspace(0, pr["t_end"], 200 if fast else 400),
                                  pr["boundary"], **pr["bp"])
    th = np.linspace(0, 2 * np.pi, 200)
    ax[0].plot(0.3 + np.cos(th), np.sin(th), color=COLORS["ink"])
    vortex_marks(ax[0], tr[-1], pr["Gamma"], size=0.04, trails=tr)
    ax[0].set_title(pr["label"], fontsize=9)
    te = np.linspace(0, 30, 121 if fast else 241)
    rd = ch05.ring_dynamics([dict(R=1.0, z=0.0, Gamma=1.0, a=0.1), dict(R=1.0, z=0.5, Gamma=1.0, a=0.1)], te)
    for i, c in enumerate((COLORS["teal"], COLORS["accent"])):
        ax[1].plot(rd["z"][:, i], rd["R"][:, i], color=c, label=f"ring {i + 1}")
    ax[1].set_xlabel("z [m]")
    ax[1].set_ylabel("R [m]")
    ax[1].legend(fontsize=8)
    ax[1].set_title("leap-frogging coaxial rings (meridional plane)", fontsize=9)
    rw = ch05.ring_dynamics([dict(R=1.0, z=0.0, Gamma=1.0, a=0.1)], np.linspace(0, 12, 121), wall_z=2.0)
    ax[2].plot(rw["z"][:, 0], rw["R"][:, 0], color=COLORS["teal"])
    ax[2].axvline(2.0, color=COLORS["ink"], lw=2)
    ax[2].set_xlabel("z [m]")
    ax[2].set_ylabel("R [m]")
    ax[2].set_title("a ring approaching a wall widens and slows", fontsize=9)
    return fig, rd, rw


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch05"))
    ap.add_argument("--no-show", action="store_true")
    ap.add_argument("--fast", action="store_true")
    args = ap.parse_args()
    import matplotlib

    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from fluidpy.core.style import use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    for G1, G2 in ((1.0, 1.0), (1.0, 3.0), (1.0, -1.0)):
        vp = ch05.vortex_pair(G1, G2, 1.0)
        print(f"pair Γ = ({G1:+g}, {G2:+g}), h = 1 m: V1 {vp['V1']:.6f}, V2 {vp['V2']:.6f} m/s, centre from 1 "
              f"{vp['centre_from_1']}, rate {vp['rotation_rate']:.6f} rad/s, translation {vp['translation_speed']:.6f}")
    pr = ch05.point_vortex_preset("equal_pair")
    tr = ch05.point_vortex_evolve(pr["xv"], pr["Gamma"], [0.0, 19.7392088])
    print(f"equal pair after one predicted period 2π/(2Γ/2πh²) = 19.7392 s: displacement "
          f"{np.abs(tr[-1] - tr[0]).max():.1e} m")
    for nm in ch05.POINT_VORTEX_PRESETS:
        pr = ch05.point_vortex_preset(nm)
        tr = ch05.point_vortex_evolve(pr["xv"], pr["Gamma"], np.linspace(0, pr["t_end"], 101), pr["boundary"], **pr["bp"])
        if pr["boundary"] is None:
            inv = ch05.point_vortex_invariants(tr, pr["Gamma"])
            print(f"{nm:15s}: invariants drift P {np.ptp(inv['P'], axis=0).max():.1e}, I {np.ptp(inv['I']):.1e}, "
                  f"H {np.ptp(inv['H']):.1e}")
        else:
            print(f"{nm:15s}: final positions {np.round(tr[-1], 4).tolist()}")
    # wall: normal velocity on the wall and drift
    P, Gs = ch05.wall_image_system([[0.0], [0.5]], [1.0])
    xw = np.linspace(-3, 3, 201)
    uw = ch05.point_vortex_velocity(np.stack([xw, 0 * xw]), P, Gs)
    print(f"wall image: max |v| on the wall {np.abs(uw[1]).max():.1e}; u(0, 0) = {uw[0][100]:.6f} = Γh/π(x² + h²) "
          f"{1.0 * 0.5 / (np.pi * 0.25):.6f}; drift {ch05.vortex_near_wall_speed(1.0, 0.5):.6f} m/s")
    for h in (0.25, 0.5, 0.1):
        print(f"channel H = 1, h = {h}: series {ch05.channel_image_velocity(h, 1.0, 1.0):.10f}, closed form "
              f"{ch05.channel_image_velocity_exact(h, 1.0, 1.0):.10f} m/s")
    # from scratch: double loop over pairs vs point_vortex_velocity
    pr = ch05.point_vortex_preset("three_vortices")
    V, Gm = pr["xv"], pr["Gamma"]
    by_hand = np.zeros_like(V)
    for k in range(V.shape[1]):
        for j in range(V.shape[1]):
            if j == k:
                continue
            rx, ry = V[0, k] - V[0, j], V[1, k] - V[1, j]
            r2 = rx * rx + ry * ry
            by_hand[:, k] += Gm[j] / (2 * np.pi) * np.array([-ry, rx]) / r2
    print(f"from scratch dx_k/dt vs point_vortex_rhs: max difference "
          f"{np.abs(by_hand - ch05.point_vortex_rhs(V, Gm)).max():.1e}")
    fig = pairs_figure(args.fast)
    fig.savefig(out / "c12_pairs.png", bbox_inches="tight")
    fig, rd, rw = images_rings_figure(args.fast)
    fig.savefig(out / "c13_images_rings.png", bbox_inches="tight")
    passes = int(np.sum(np.diff(np.sign(rd["z"][:, 1] - rd["z"][:, 0])) != 0))
    print(f"leap-frog: {passes} pass-throughs in 30 s, impulse ΣΓπR² drift {np.ptp(rd['impulse']):.1e}")
    speed = np.diff(rw["z"][:, 0]) / np.diff(rw["t"])
    print(f"ring toward a wall: R {rw['R'][0, 0]:.3f} → {rw['R'][-1, 0]:.3f} m (monotone: "
          f"{bool(np.all(np.diff(rw['R'][:, 0]) > 0))}), approach speed {speed[0]:.4f} → {speed[-1]:.4f} m/s "
          f"(decreasing: {bool(np.all(np.diff(speed) < 0))})")
    print(f"figures → {out}")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
