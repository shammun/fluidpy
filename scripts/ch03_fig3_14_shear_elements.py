"""Fig. 3.14 (our drawing): two small elements in the parallel shear flow u = (γ x₂, 0) — the aligned square ABCD
shears without stretching, the 45° square PQRS stretches along x̄₁ and shrinks along x̄₂ without shearing, and both
turn at the same angular velocity −γ/2 = ω₃/2 (§3.5, N34, N35).

Run: ``.venv/Scripts/python.exe scripts/ch03_fig3_14_shear_elements.py --no-show``
Figure → outputs/ch03/fig3_14_shear_elements.png.
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

from fluidpy import ch03_kinematics as ch03  # noqa: E402
from fluidpy.core.style import COLORS  # noqa: E402


def _corners(kind: str, half: float = 0.5) -> np.ndarray:
    """Corners (2, 4) of the aligned square ABCD (A top-left, B bottom-left, C bottom-right, D top-right) or of the
    45° square PQRS inscribed in it (P left, Q bottom, R right, S top), centred at the origin."""
    if kind == "ABCD":
        return np.array([[-half, -half, half, half], [half, -half, -half, half]])
    return np.array([[-half, 0.0, half, 0.0], [0.0, -half, 0.0, half]])


def _sides_angles(P: np.ndarray):
    """Side lengths (4,) and corner angles (4,) [deg] of the closed quadrilateral with corners P (2, 4)."""
    sides = np.hypot(*(np.roll(P, -1, axis=1) - P))
    ang = []
    for k in range(4):
        a, b = P[:, k - 1] - P[:, k], P[:, (k + 1) % 4] - P[:, k]
        ang.append(np.degrees(np.arccos(np.clip(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)), -1, 1))))
    return sides, np.array(ang)


def shear_elements_frames(gamma: float, times) -> list[dict]:
    """Positions (relative to the element centre) of ABCD and PQRS carried by u = (γ x₂, 0) at each time, with their
    side lengths, corner angles and the turning angle of each side. Units: γ [1/s], times [s], lengths [m].

    Returns a list of dicts: t, ABCD (2, 4), PQRS (2, 4), ABCD_sides, ABCD_angles_deg, PQRS_sides, PQRS_angles_deg.
    """
    G = ch03.velocity_gradient_preset("simple_shear", gamma)
    out = []
    for t in np.atleast_1d(times):
        M = ch03.linear_flow_map(G, t)  # x(t) = e^{Gt} x₀ (exact for the linear shear)
        A, P = M @ _corners("ABCD"), M @ _corners("PQRS")
        sa, aa = _sides_angles(A)
        sp_, ap = _sides_angles(P)
        out.append({"t": float(t), "ABCD": A, "PQRS": P, "ABCD_sides": sa, "ABCD_angles_deg": aa, "PQRS_sides": sp_,
                    "PQRS_angles_deg": ap})
    return out


def shear_elements_figure(gamma: float = 1.0, t: float = 0.3, ax=None):
    """Both elements at t = 0 (thin) and at t (bold) in one axes; principal axes at 45° drawn. Returns (fig, ax)."""
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(6.2, 4.6))
    else:
        fig = ax.figure
    f0, f1 = shear_elements_frames(gamma, [0.0, t])
    for f, lw, al in ((f0, 1.0, 0.45), (f1, 2.4, 1.0)):
        for key, col in (("ABCD", COLORS["accent"]), ("PQRS", COLORS["teal"])):
            P = np.hstack([f[key], f[key][:, :1]])
            ax.plot(*P, color=col, lw=lw, alpha=al, label=f"{key} at t = {f['t']:.2f} s" if lw > 2 else None)
    for lbl, xy in zip("ABCD", f1["ABCD"].T):
        ax.annotate(lbl, xy, xytext=(4, 4), textcoords="offset points", color=COLORS["accent"], fontsize=10)
    for lbl, xy in zip("PQRS", f1["PQRS"].T):
        ax.annotate(lbl, xy, xytext=(4, -10), textcoords="offset points", color=COLORS["teal"], fontsize=10)
    L = 0.75
    ax.plot([-L, L], [-L, L], ls=":", color=COLORS["blue"], lw=1.2, label=r"$\bar x_1$: stretching ($+\gamma/2$)")
    ax.plot([-L, L], [L, -L], ls=":", color=COLORS["rose"], lw=1.2, label=r"$\bar x_2$: compression ($-\gamma/2$)")
    ax.set_aspect("equal")
    ax.set_xlabel("$x_1$ relative to the element centre [m]")
    ax.set_ylabel("$x_2$ [m]")
    ax.set_title(rf"Parallel shear flow, $\gamma$ = {gamma:g} s$^{{-1}}$: shear (ABCD) vs stretch (PQRS)")
    ax.legend(fontsize=8, loc="lower right")
    return fig, ax


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch03"))
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
    gamma, t = 1.0, 0.3
    fig, _ = shear_elements_figure(gamma, t)
    fig.savefig(out / "fig3_14_shear_elements.png", bbox_inches="tight")
    k = ch03.parallel_shear_kinematics(gamma)
    print(f"γ = {gamma}: ω₃ = {k['omega3']:+.3f} s⁻¹, spin = {k['spin']:+.3f} rad/s, AB turns {k['rate_AB']:+.3f}, "
          f"BC {k['rate_BC']:+.3f}; principal rates {k['lam']} at {k['principal_angle_deg']:.1f}°")
    th = np.linspace(0, np.pi, 7)
    print("pair-average spin over θ:", np.round(ch03.perpendicular_pair_rotation_rate(k["G"], th), 12))
    for dt in (0.1, 0.01):
        f = shear_elements_frames(gamma, [dt])[0]
        print(f"after {dt} s: ABCD sides {f['ABCD_sides'].round(6)} angles {f['ABCD_angles_deg'].round(3)}; "
              f"PQRS sides {f['PQRS_sides'].round(6)} angles {f['PQRS_angles_deg'].round(3)}")
    print(f"figure → {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
