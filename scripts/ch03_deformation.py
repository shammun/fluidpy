"""§3.4: the motion of a small fluid element — relative velocity (3.10) and its Taylor remainder (C06), stretching
along any direction measured on tracked segments vs n·S·n (C07, the rose curve), the closing of a right angle (C08),
volume growth (C09, (3.14)), the spin of perpendicular pairs (C10), the strain + rotation split (3.19) (C11) and a
circle becoming an ellipse on the principal axes (C12) — with the curation's worked numbers.

Run: ``.venv/Scripts/python.exe scripts/ch03_deformation.py --no-show``
Figures → outputs/ch03/c07_stretch_rose.png, c10_pair_spin.png, c12_circle_to_ellipse.png.
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
from tools.convergence import observed_order  # noqa: E402


def stretch_rose_figure(G, dt: float = 0.01):
    """Stretching rate (1/ℓ)dℓ/dt measured on tracked unit segments at every angle (dots) vs the formula n·S·n
    (line), plus the shear-closing rate of each perpendicular pair vs 2n₁·S·n₂."""
    import matplotlib.pyplot as plt

    th = np.linspace(0, np.pi, 73)
    n = np.stack([np.cos(th), np.sin(th)])
    formula = ch03.linear_strain_rate(G, n)
    meas = np.array([ch03.measured_strain_rates(G, n[:, k], dt) for k in range(th.size)])
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(np.degrees(th), formula, color=COLORS["blue"], label=r"$\mathbf{n}\cdot\mathbf{S}\cdot\mathbf{n}$")
    ax.plot(np.degrees(th[::3]), [m["stretch"] for m in meas[::3]], "o", color=COLORS["blue"], ms=4,
            label=f"measured stretch (dt = {dt} s)")
    ax.plot(np.degrees(th), [m["closing_formula"] for m in meas], color=COLORS["teal"],
            label=r"$2\,\mathbf{n}_1\cdot\mathbf{S}\cdot\mathbf{n}_2$")
    ax.plot(np.degrees(th[::3]), [m["closing"] for m in meas[::3]], "s", color=COLORS["teal"], ms=4,
            label="measured closing of the right angle")
    ax.axhline(0, color=COLORS["grid"])
    ax.set_xlabel(r"direction of the segment $\theta$ [deg]")
    ax.set_ylabel("rate [1/s]")
    ax.set_title("every entry of S is a measurable rate")
    ax.legend(fontsize=8)
    return fig


def pair_spin_figure(gammas=(0.5, 1.0, 2.0)):
    """Single-line turning rate θ̇(θ) = −γ sin²θ in the parallel shear flow and the flat perpendicular-pair average
    −γ/2 = ω₃/2 for several γ (C10, D14)."""
    import matplotlib.pyplot as plt

    th = np.linspace(0, np.pi, 181)
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for g, col in zip(gammas, (COLORS["teal"], COLORS["accent"], COLORS["orange"])):
        G = ch03.velocity_gradient_preset("simple_shear", g)
        ax.plot(np.degrees(th), ch03.material_line_rotation_rate(G, th), color=col, label=rf"line, $\gamma$ = {g}")
        ax.plot(np.degrees(th), ch03.perpendicular_pair_rotation_rate(G, th), "--", color=col,
                label=rf"pair average = $\omega_3/2$ = {-g / 2}")
    ax.set_xlabel(r"line angle $\theta$ [deg]")
    ax.set_ylabel(r"$d\theta/dt$ [rad/s]")
    ax.set_title("single lines turn at different rates; perpendicular pairs average to ω₃/2")
    ax.legend(fontsize=7, ncol=2)
    return fig


def circle_to_ellipse_figure(G, times=(0.0, 0.25, 0.5, 1.0)):
    """A circle of material points in u = G·x at several times with the principal axes of S (C12) and the exact
    ellipse axes (singular vectors of e^{Gt}) at the last time."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5.4, 5.0))
    for t, al in zip(times, np.linspace(0.3, 1.0, len(times))):
        P = ch03.deform_circle(G, t)
        ax.plot(*np.hstack([P, P[:, :1]]), color=COLORS["accent"], alpha=al, lw=1.6, label=f"t = {t} s")
    lam, C = ch03.principal_strain_rates(G)
    for k, col in ((1, COLORS["blue"]), (0, COLORS["rose"])):
        ax.plot([-1.6 * C[0, k], 1.6 * C[0, k]], [-1.6 * C[1, k], 1.6 * C[1, k]], ":", color=col,
                label=rf"principal axis, $\lambda$ = {lam[k]:+.2f}")
    s, D = ch03.strain_ellipse_axes(G, times[-1], method="exact")
    for k in range(2):
        ax.plot([0, s[k] * D[0, k]], [0, s[k] * D[1, k]], color=COLORS["ink"], lw=1.2)
    ax.set_aspect("equal")
    ax.set_xlabel("$dx_1$ [m]")
    ax.set_ylabel("$dx_2$ [m]")
    ax.set_title("circle → ellipse (black: exact axes at the last time)")
    ax.legend(fontsize=7)
    return fig


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
    Gdemo = np.array([[0.6, 1.0], [-0.2, -0.3]])
    stretch_rose_figure(Gdemo).savefig(out / "c07_stretch_rose.png", bbox_inches="tight")
    pair_spin_figure().savefig(out / "c10_pair_spin.png", bbox_inches="tight")
    circle_to_ellipse_figure(ch03.velocity_gradient_preset("simple_shear", 1.0)).savefig(out / "c12_circle_to_ellipse.png",
                                                                                         bbox_inches="tight")
    # C06
    G6 = np.array([[1.0, 2.0], [0.0, -1.0]])
    print(f"C06: du = G·dx = {ch03.relative_velocity(G6, [0.01, 0.02])} m/s")
    u_nl = lambda x, t: np.stack([np.sin(x[0]) * np.cos(x[1]), x[0] * x[1] ** 2])  # noqa: E731
    x0 = np.array([0.4, 0.7])
    Gx = ch03.velocity_gradient_at(u_nl, x0)
    ds = np.array([0.1, 0.05, 0.025, 0.0125])
    e = [np.linalg.norm(u_nl(x0 + d * np.array([0.6, 0.8]), 0) - u_nl(x0, 0) - ch03.relative_velocity(Gx, d * np.array([0.6, 0.8])))
         for d in ds]
    print(f"     Taylor remainder |u(x+dx) − u(x) − G·dx| slope {observed_order(ds, e):.3f}")
    # C07
    G7 = ch03.velocity_gradient_preset("pure_strain", 2.0)
    l1 = np.linalg.norm(ch03.linear_flow_map(G7, 0.01) @ np.array([0.01, 0.0]))
    print(f"C07: 1 cm segment along x in u = (2x, −2y) after 0.01 s: {l1 * 100:.4f} cm (rate {ch03.linear_strain_rate(G7, [1, 0]):.1f} s⁻¹)")
    dts = np.array([1e-2, 5e-3, 2.5e-3, 1.25e-3])
    err = [abs(ch03.measured_strain_rates(Gdemo, [0.6, 0.8], d)["stretch"] - ch03.linear_strain_rate(Gdemo, [0.6, 0.8])) for d in dts]
    print(f"     measured → formula at order {observed_order(dts, err):.3f} in dt")
    # C08
    G8 = ch03.velocity_gradient_preset("simple_shear", 1.0)
    M = ch03.linear_flow_map(G8, 0.01)
    da = np.arctan2((M @ [0, 1])[0], (M @ [0, 1])[1])
    db = np.arctan2((M @ [1, 0])[1], (M @ [1, 0])[0])
    print(f"C08: γ = 1, dt = 0.01 s: dα = {da:.6f} rad, dβ = {db:.6f} rad → S₁₂ = ½(dα + dβ)/dt = {0.5 * (da + db) / 0.01:.4f} s⁻¹ "
          f"(formula {ch03.shear_strain_rate(G8, [1, 0], [0, 1]):.4f})")
    print(f"     rigid motion: S of U + Ω × x = {np.abs(ch03.strain_rate_tensor(ch03.velocity_gradient_at(lambda x, t: ch03.rigid_body_velocity([1, 2, 3], [0.3, -0.2, 0.5], x), np.array([0.2, 0.1, -0.4])))).max():.1e}")
    # C09
    G9 = np.eye(3)
    print(f"C09: u = (x, y, z): ∇·u = {ch03.volumetric_strain_rate(G9):.0f} s⁻¹; 1 cm³ box after 0.01 s → "
          f"{ch03.material_volume_ratio(G9, 0.01):.4f} cm³ (e^0.03 = {np.exp(0.03):.4f})")
    # C10
    Gs = ch03.velocity_gradient_preset("solid_body_rotation", 1.0)
    print(f"C10: solid body ω₀ = 1: lines turn at {ch03.material_line_rotation_rate(Gs, np.array([0, 0.7, 2.0]))}, ω₃ = "
          f"{ch03.vorticity_from_gradient(Gs)[2]:.1f}, spin {ch03.element_rotation_rate(Gs)[2]:.1f}; "
          f"shear γ = 1 pair average {ch03.perpendicular_pair_rotation_rate(G8, 0.3):.3f}")
    print(f"     rotating observer Ω = ω₃/2 = 1 in solid body: ω′ = {ch03.vorticity_in_rotating_frame(2.0, 1.0):.1f}")
    # C11
    sp_ = ch03.relative_velocity_split(G8, [0.0, 1.0])
    print(f"C11: shear γ = 1, dx = (0, 1): du = {sp_.du} = S·dx {sp_.du_strain} + ½ω × dx {sp_.du_rot}")
    # C12
    ax_fo, _ = ch03.strain_ellipse_axes(G8, 0.1, method="first_order", radius=1e-3)
    ax_s, _ = ch03.strain_ellipse_axes(G8, 0.1, method="strain_only", radius=1e-3)
    ax_e, _ = ch03.strain_ellipse_axes(G8, 0.1, method="exact", radius=1e-3)
    print(f"C12: 1 mm circle after 0.1 s in shear γ = 1: semi-axes first order {ax_fo * 1e3} mm, strain only "
          f"{(ax_s * 1e3).round(4)} mm, exact {(ax_e * 1e3).round(4)} mm; λ = {ch03.principal_strain_rates(G8)[0]}")
    du_bar, dx_bar = ch03.strain_velocity_principal(Gdemo, [0.3, -0.4])
    print(f"     (3.21) in the eigenframe: dū/dx̄ = {du_bar / dx_bar} = λ {ch03.principal_strain_rates(Gdemo)[0]}")
    print(f"figures → {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
