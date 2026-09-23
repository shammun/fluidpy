"""§4.3 stream functions (C03, E2): ψ contours at equal Δψ over the speed for E2's presets, the flux through a gate
equal to ψ₂ − ψ₁ whatever its angle, the axisymmetric uniform stream, and the two 3-D stream-function pairs of Fig. 4.1
with the flux (b − a)(d − c) through the patch they bound.

Run: ``.venv/Scripts/python.exe scripts/ch04_stream_functions.py --no-show``
Figures → outputs/ch04/c03_stream_function_presets.png, c03_gate_flux.png.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from fluidpy import ch04_conservation_laws as ch04  # noqa: E402
from fluidpy.core.style import COLORS  # noqa: E402


def presets_figure(n_levels: int = 21):
    """Six planar presets: equal-Δψ contours (the streamlines) over |u| — crowded contours = fast flow."""
    import matplotlib.pyplot as plt

    names = ["uniform", "stagnation", "source_stream", "cylinder", "vortex", "shear"]
    x = np.linspace(-3, 3, 241)
    X, Y = np.meshgrid(x, x)
    fig, axs = plt.subplots(2, 3, figsize=(13, 8.5))
    for ax, name in zip(axs.ravel(), names):
        kw = dict(alpha=0.4) if name == "uniform" else {}
        psi = np.asarray(ch04.streamfunction_preset(name, X, Y, **kw), float)
        u, v = ch04.velocity_preset(name, X, Y, **kw)
        spd = np.hypot(u, v)
        if name == "cylinder":
            inside = X ** 2 + Y ** 2 < 1.0
            psi, spd = np.where(inside, np.nan, psi), np.where(inside, np.nan, spd)
        im = ax.pcolormesh(X, Y, np.clip(spd, 0, 3), shading="auto", cmap="viridis")
        lv = np.linspace(np.nanpercentile(psi, 2), np.nanpercentile(psi, 98), n_levels)
        ax.contour(X, Y, psi, levels=lv, colors="white", linewidths=0.8)
        ax.set_aspect("equal")
        ax.set_title(name, fontsize=10)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        fig.colorbar(im, ax=ax, shrink=0.8, label="|u| [m/s]")
    fig.suptitle(r"$\rho u=\partial\psi/\partial y,\ \rho v=-\partial\psi/\partial x$: equal Δψ between contours")
    return fig


def gate_figure():
    """The flux through a gate of fixed ends but varying angle/shape equals ψ(p2) − ψ(p1)."""
    import matplotlib.pyplot as plt

    p1 = np.array([1.5, 0.2])
    ang = np.linspace(-1.2, 1.2, 25)
    fluxes, dpsi = [], []
    for a in ang:
        p2 = p1 + 1.2 * np.array([np.cos(np.pi / 2 + a), np.sin(np.pi / 2 + a)])
        fluxes.append(ch04.flux_between_streamlines("cylinder", p1, p2))
        dpsi.append(ch04.streamfunction_preset("cylinder", *p2) - ch04.streamfunction_preset("cylinder", *p1))
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(np.degrees(ang), fluxes, "o", color=COLORS["accent"], label="∫u·n ds along the gate (quadrature)")
    ax.plot(np.degrees(ang), dpsi, "-", color=COLORS["orange"], label="ψ(p₂) − ψ(p₁)")
    ax.set_xlabel("gate tilt [°]")
    ax.set_ylabel("flux per unit depth [m²/s]")
    ax.set_title("cylinder flow: the flux between two streamlines is Δψ")
    ax.legend(fontsize=8)
    return fig, np.max(np.abs(np.array(fluxes) - np.array(dpsi)))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch04"))
    ap.add_argument("--no-show", action="store_true")
    args = ap.parse_args()
    t0 = time.perf_counter()
    import matplotlib

    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from fluidpy.core.style import use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    presets_figure().savefig(out / "c03_stream_function_presets.png", bbox_inches="tight")
    fig, err = gate_figure()
    fig.savefig(out / "c03_gate_flux.png", bbox_inches="tight")
    print(f"C03 gate flux vs Δψ: max difference {err:.1e} m²/s over 25 tilts")
    bent = np.array([[1.5, 2.2, 1.9, 1.2], [0.2, 0.6, 1.1, 1.4]])
    print(f"C03 bent gate flux {ch04.flux_along_path('cylinder', bent):.10f} vs Δψ "
          f"{ch04.streamfunction_preset('cylinder', 1.2, 1.4) - ch04.streamfunction_preset('cylinder', 1.5, 0.2):.10f}")
    print("C03 u, v at (2, 0.5) from ψ vs ch03-type closed form:", ch04.velocity_from_streamfunction_2d("cylinder", 2.0, 0.5),
          ch04.velocity_preset("cylinder", 2.0, 0.5))
    print("N17 axisymmetric uniform stream ψ = ½UR² (U = 2): (u_R, u_z) =",
          ch04.velocity_from_streamfunction_axisym("axisym_uniform", 0.5, 0.2, U=2.0))
    num, cl = ch04.stream_tube_mass_flux(a=0.0, b=2.0, c=0.0, d=0.5)
    print(f"N16 parabolic pair χ = y, ψ = z − x²: patch flux {num:.10f} vs (b − a)(d − c) = {cl}")
    ex = ch04.stream_surface_example()
    print(f"N16 second pair χ = x² + y, ψ = z − y²: {ch04.stream_tube_mass_flux(ex['chi'], ex['psi'], patch=ex['patch'])[0]:.10f}"
          f" vs {ex['expected']}")
    ce, pe, cs = ch04.stream_function_pair()
    print("N15 sympy ∇·(∇χ×∇ψ), u·∇χ, u·∇ψ:", ch04.stream_surface_check(ce, pe, cs))
    print(f"figures → {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
