"""Chapter 1, §1.11, Example 1.5: Rayleigh scattering S/I = V² φ₃(n_s)/(d² λ⁴) — why small particles scatter blue.

Run: ``.venv/Scripts/python.exe scripts/ch01_blue_sky.py --no-show``
Figure → outputs/ch01/fig_rayleigh_spectrum.png (curves normalised to 700 nm; φ₃ taken constant across the visible).
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch01"), help="figure folder")
    ap.add_argument("--no-show", action="store_true", help="do not open a window")
    args = ap.parse_args()
    t_start = time.perf_counter()
    import matplotlib

    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from fluidpy import ch01_introduction as ch01
    from fluidpy.core.style import COLORS, use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    lam_nm = np.linspace(380.0, 750.0, 371)
    V, d = (10e-9) ** 3, 1.0  # a 10 nm particle seen from 1 m
    ratio = ch01.rayleigh_scattering_ratio(V, d, lam_nm * 1e-9)
    rel = ratio / ch01.rayleigh_scattering_ratio(V, d, 700e-9)
    r, g, b = ch01.wavelength_to_rgb(lam_nm)
    fig, ax = plt.subplots()
    ax.scatter(lam_nm, rel, c=np.stack([r, g, b], axis=-1), s=12)
    ax.plot(lam_nm, rel, color=COLORS["ink"], lw=0.8)
    ax.set_xlabel("wavelength λ [nm]")
    ax.set_ylabel("scattered intensity relative to 700 nm [-]")
    ax.set_title(r"Rayleigh scattering $\propto \lambda^{-4}$ (Example 1.5)")
    fig.savefig(out / "fig_rayleigh_spectrum.png", bbox_inches="tight")

    blue_red = ch01.rayleigh_scattering_ratio(V, d, 450e-9) / ch01.rayleigh_scattering_ratio(V, d, 700e-9)
    print(f"S/I for a 10 nm particle at 1 m, 550 nm (phi3 = 1): {ch01.rayleigh_scattering_ratio(V, d, 550e-9):.3e}")
    print(f"blue (450 nm) / red (700 nm) = {blue_red:.3f} = (700/450)^4 = {(700 / 450) ** 4:.3f}")
    print(f"doubling the particle volume: x{ch01.rayleigh_scattering_ratio(2 * V, d, 550e-9) / ch01.rayleigh_scattering_ratio(V, d, 550e-9):.1f}; "
          f"doubling the distance: x{ch01.rayleigh_scattering_ratio(V, 2 * d, 550e-9) / ch01.rayleigh_scattering_ratio(V, d, 550e-9):.2f}")
    print(f"saved fig_rayleigh_spectrum.png in {out}  ({time.perf_counter() - t_start:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
