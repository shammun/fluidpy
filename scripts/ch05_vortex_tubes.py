"""§5.1 vortex lines and tubes (C01, N02, N03, N45): the narrowing Gaussian tube keeps its strength (5.4) while its
area shrinks and its mean vorticity grows; the "broken" (non-solenoidal) field fails Gauss' budget; a traced, twisted
tube closes its polyhedral budget to round-off; tube strength by the circulation and the flux routes; vortex lines of
the helical swirl keep zR² constant.

Run: ``.venv/Scripts/python.exe scripts/ch05_vortex_tubes.py --no-show``
Figures → outputs/ch05/c01_tube_sections.png, c01_traced_tube.png.
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
from fluidpy.core.integral_theorems import planar_disc, planar_loop  # noqa: E402
from fluidpy.core.style import COLORS  # noqa: E402


def sections_figure():
    """Area, flux and mean vorticity along the narrowing tube (flux flat, mean ω ∝ 1/area) vs the broken field."""
    import matplotlib.pyplot as plt

    z = np.linspace(0.0, 2.0, 41)
    sec = [ch05.gaussian_tube_section(zz) for zz in z]
    broken = [ch05.tube_flux_budget("broken", 0.0, zz, 0.1).upper if zz > 0 else 0.632121 for zz in z]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].plot(z, [s["flux"] for s in sec], color=COLORS["teal"], label="flux Γ through the tube (5.4)")
    ax[0].plot(z, broken, "--", color=COLORS["rose"], label="'broken' field (∇·ω ≠ 0)")
    ax[0].set_xlabel("height z along the tube [m]")
    ax[0].set_ylabel("flux ∫ω·n dA [m²/s]")
    ax[0].legend(fontsize=8)
    ax[0].set_title("tube strength is the same everywhere")
    a2 = ax[1].twinx()
    ax[1].plot(z, [s["area"] for s in sec], color=COLORS["blue"], label="area")
    a2.plot(z, [s["mean_omega"] for s in sec], color=COLORS["accent"], label="mean ω")
    ax[1].set_xlabel("height z [m]")
    ax[1].set_ylabel("cross-section area [m²]", color=COLORS["blue"])
    a2.set_ylabel("mean vorticity [1/s]", color=COLORS["accent"])
    ax[1].set_title("where the tube thins, ω grows")
    return fig


def traced_figure(n_lines: int):
    """3-D view of a traced twisted tube with its end sections; returns (fig, budget dict)."""
    import matplotlib.pyplot as plt
    from ch05_drawings import tube_mesh

    tube = ch05.vortex_tube("gaussian_tube", (0, 0, 0.0), (0, 0, 1), 0.08, n_lines=n_lines, s_max=0.8, n=41,
                            both=False, twist=3.0)
    bud = ch05.tube_flux_budget_traced("gaussian_tube", tube, twist=3.0)
    fig = plt.figure(figsize=(6, 6))
    fig.set_layout_engine("none")
    ax = fig.add_subplot(projection="3d")
    tube_mesh(ax, tube)
    f = bud["flux"]
    ax.set_title(f"traced tube: lower {f.lower:+.6f}, upper {f.upper:+.6f}, side {f.side:.1e}, total {f.total:.1e}",
                 fontsize=9)
    return fig, bud


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

    for zz in (0.0, 1.0):
        s = ch05.gaussian_tube_section(zz)
        print(f"section z = {zz}: radius {s['radius']:.6f} m, area {s['area']:.6f} m², flux {s['flux']:.6f} m²/s, "
              f"mean ω {s['mean_omega']:.3f} 1/s, peak ω {s['peak_omega']:.3f} 1/s")
    for field in ("gaussian_tube", "broken"):
        b = ch05.tube_flux_budget(field, 0.0, 1.0, 0.1)
        print(f"(5.4) budget {field:13s}: lower {b.lower:+.6f}, side {b.side:+.2e}, upper {b.upper:+.6f}, "
              f"total {b.total:+.2e}")
    u = ch05.lamb_oseen_field(1.0, 1e-3, 2.5)  # σ_c = 0.1 m
    ts = ch05.vortex_tube_strength(u, planar_loop((0.0, 0.0), radius=0.1, n=256),
                                   planar_disc((0.0, 0.0), radius=0.1, nr=256, ntheta=128))
    print(f"Lamb–Oseen tube strength at r = σ_c: circulation {ts.circulation:.6f}, flux {ts.flux:.6f} "
          f"(exact 1 − e⁻¹ = {1 - np.exp(-1):.6f})")
    W = ch05.vorticity_field(ch05.helical_swirl_field(1.0), 1e-5)
    _, X = ch05.vortex_line(W, (0.5, 0.0, 1.0), s_max=0.5, n=60)
    zR2 = X[2] * (X[0] ** 2 + X[1] ** 2)
    print(f"helical swirl vortex line: zR² = {zR2.mean():.6f} ± {np.ptp(zR2):.1e} (constant, (5.3))")
    fig = sections_figure()
    fig.savefig(out / "c01_tube_sections.png", bbox_inches="tight")
    fig, bud = traced_figure(16 if args.fast else 32)
    fig.savefig(out / "c01_traced_tube.png", bbox_inches="tight")
    print(f"traced twisted tube: areas {bud['area_lower']:.5f} → {bud['area_upper']:.5f} m², mean ω "
          f"{bud['mean_omega_lower']:.2f} → {bud['mean_omega_upper']:.2f} 1/s, closed-surface total "
          f"{bud['flux'].total:+.1e} m²/s")
    print(f"figures → {out}")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
