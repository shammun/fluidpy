"""§4.7 noninertial frames (C09, E5): the pole projectile seen from the inertial and the rotating frame (Fig. 4.8), its
deflection vs the Ωut² estimate, a check of the five-term acceleration (4.43) on a path known in both frames, effective
gravity vs latitude (Fig. 4.9), the Earth's oblateness, and Example 4.5's pump equations by sympy.

Run: ``.venv/Scripts/python.exe scripts/ch04_rotating_frames.py --no-show``
Figures → outputs/ch04/c09_coriolis_projectile.png, c09_effective_gravity.png.
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


def projectile_figure(u0: float = 10.0, hours: float = 6.0):
    """The same projectile: a straight line from space, a curve from the rotating Earth (NH), deflection vs Ωut²."""
    import matplotlib.pyplot as plt

    t = np.linspace(0, hours * 3600, 241)
    r = ch04.coriolis_projectile(u0, ch04.OMEGA_EARTH, t)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.5))
    a1.plot(*r["inertial"] / 1e3, color=COLORS["muted"], label="inertial (straight)")
    a1.plot(*r["rotating"] / 1e3, color=COLORS["amber"], label="rotating frame (ODE with Coriolis + centrifugal)")
    a1.set_aspect("equal")
    a1.set_xlabel("x′ [km]")
    a1.set_ylabel("y′ [km]")
    a1.set_title(f"Fig. 4.8 analogue: u₀ = {u0} m/s from the pole, {hours:g} h")
    a1.legend(fontsize=8)
    a2.plot(t / 3600, r["deflection"] / 1e3, color=COLORS["amber"], label="u t sin Ωt (exact)")
    a2.plot(t / 3600, r["deflection_small"] / 1e3, "--", color=COLORS["accent"], label="Ωut² (book)")
    a2.set_xlabel("t [h]")
    a2.set_ylabel("deflection to the right [km]")
    a2.legend(fontsize=8)
    return fig, r


def gravity_figure():
    """Effective gravity magnitude and its deflection from the radial vs latitude (spherical Earth)."""
    import matplotlib.pyplot as plt

    lat = np.radians(np.linspace(0, 90, 181))
    ge, dev = ch04.effective_gravity(lat)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4))
    a1.plot(np.degrees(lat), ge, color=COLORS["blue"])
    a1.set_xlabel("latitude [°]")
    a1.set_ylabel(r"$|g_e|$ [m/s²]")
    a1.set_title(r"$g_e = g + \Omega^2R\,e_R$ (g_n = 9.8 m/s²)")
    a2.plot(np.degrees(lat), np.degrees(dev) * 60, color=COLORS["accent"])
    a2.set_xlabel("latitude [°]")
    a2.set_ylabel("deflection from the radial [arcmin]")
    a2.set_title(f"largest near 45°: {np.degrees(lat[np.argmax(dev)]):.1f}°")
    return fig


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
    fig, r = projectile_figure()
    fig.savefig(out / "c09_coriolis_projectile.png", bbox_inches="tight")
    gravity_figure().savefig(out / "c09_effective_gravity.png", bbox_inches="tight")
    print(f"C09 ODE path vs rotated straight line: max diff {np.abs(r['rotating'] - ch04.projectile_paths(10.0, ch04.OMEGA_EARTH, r['t'])[1]).max():.1e} m")
    c = ch04.coriolis_projectile(10.0, ch04.OMEGA_EARTH, [3600.0])
    print(f"C09 1 h: forward {c['forward'][0]:.0f} m, deflection {c['deflection'][0]:.1f} m, Ωut² {c['deflection_small'][0]:.1f} m, "
          f"angle {np.degrees(c['angle'][0]):.2f}°")
    c2 = ch04.coriolis_projectile(10.0, 1e-3, np.linspace(0, 600, 50), centrifugal=False)
    print(f"C09 Coriolis only: speed spread {np.ptp(c2['speed']):.1e} m/s (does no work)")
    print("C09 Coriolis acceleration for u = (10, 0, 0) at the pole:", ch04.coriolis_acceleration(ch04.OMEGA_EARTH, [10.0, 0, 0]))
    # five-term check on a path known in the inertial frame
    Oz, t, hh = 0.4, 1.1, 1e-4
    Rz = lambda a: np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])  # noqa: E731
    xin = lambda s: np.array([s ** 2, np.sin(s), 0.3 * s])  # noqa: E731
    xr = lambda s: Rz(-Oz * s) @ xin(s)  # noqa: E731
    up = (xr(t + hh) - xr(t - hh)) / (2 * hh)
    ap_ = (xr(t + hh) - 2 * xr(t) + xr(t - hh)) / hh ** 2
    terms = ch04.frame_acceleration_terms(ap_, up, xr(t), [0, 0, Oz])
    print(f"D15 (4.43) five terms vs inertial acceleration: diff {np.abs(terms['total'] - Rz(-Oz * t) @ np.array([2.0, -np.sin(t), 0])).max():.1e}")
    ge, dev = ch04.effective_gravity(np.radians([0.0, 45.0, 90.0]))
    print(f"N64 effective gravity (0°, 45°, 90°): {np.round(ge, 5)} m/s², deflection {np.round(np.degrees(dev) * 60, 3)} arcmin; "
          f"Ω²a = {ch04.OMEGA_EARTH ** 2 * ch04.EARTH_A:.5f} m/s²; 2(a − b) = {ch04.earth_oblateness_diameter() / 1e3:.2f} km")
    t1 = time.perf_counter()
    eqs = ch04.rotating_pump_equations()
    print(f"Ex. 4.5 ({time.perf_counter() - t1:.1f} s): rotation terms {ch04.rotating_pump_terms()['rotation_terms']}")
    for e in eqs:
        print("   ", e)
    print(f"figures → {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
