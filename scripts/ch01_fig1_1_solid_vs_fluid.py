"""Chapter 1, §1.3: solid vs fluid (Fig. 1.1 idea), plastic and viscoelastic materials, normal vs shear stress,
molecular spacing in liquids and gases.

Run: ``.venv/Scripts/python.exe scripts/ch01_fig1_1_solid_vs_fluid.py --no-show``
Figures → outputs/ch01/fig1_1_strain_history.png, fig1_1_elements.png. Material parameters are illustrative (ours).
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

    # --- strain histories: same stress on four materials, load removed at t_off --------------------------------------
    t = np.linspace(-0.5, 6.0, 651)          # time [s]
    tau, G, mu, t_off, tau_y = 50.0, 500.0, 100.0, 3.0, 20.0  # Pa, Pa, Pa s, s, Pa (illustrative)
    curves = {
        "elastic solid (τ/G)": ch01.shear_deformation_history(t, tau, "solid", G=G, t_off=t_off),
        "Newtonian fluid (τt/μ)": ch01.shear_deformation_history(t, tau, "fluid", mu=mu, t_off=t_off),
        f"Bingham plastic (τ_y = {tau_y:g} Pa)": ch01.shear_deformation_history(t, tau, "bingham", mu=mu, t_off=t_off,
                                                                              tau_y=tau_y),
        "Maxwell viscoelastic": ch01.shear_deformation_history(t, tau, "maxwell", G=G, mu=mu, t_off=t_off),
    }
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    for (label, g), c in zip(curves.items(), [COLORS["blue"], COLORS["accent"], COLORS["amber"], COLORS["teal"]]):
        ax.plot(t, g, label=label, color=c)
    ax.axvspan(0, t_off, color=COLORS["orange"], alpha=0.08, label="shear stress applied")
    ax.set_xlabel("time t [s]")
    ax.set_ylabel("shear strain γ [rad]")
    ax.set_title("Same shear stress, four materials: only the fluid keeps deforming")
    ax.legend(fontsize=8)
    fig.savefig(out / "fig1_1_strain_history.png", bbox_inches="tight")

    # --- element shapes at a few times (drawn from γ) -----------------------------------------------------------------
    fig2, axs = plt.subplots(1, 2, figsize=(8, 3.2))
    for ax, kind, kw, title in ((axs[0], "solid", dict(G=G), "solid: fixed shape while loaded"),
                                (axs[1], "fluid", dict(mu=mu), "fluid: keeps deforming")):
        for tk, alpha in zip([0.0, 1.0, 2.0, 2.99], [0.25, 0.45, 0.65, 0.9]):
            g = ch01.shear_deformation_history(tk, tau, kind, t_off=t_off, **kw)
            dx = np.tan(g)  # top edge offset of a unit-height element
            ax.plot([0, 1, 1 + dx, dx, 0], [0, 0, 1, 1, 0], color=COLORS["accent"], alpha=alpha)
        ax.set_aspect("equal")
        ax.set_xlim(-0.2, 2.8)
        ax.set_title(title)
        ax.axis("off")
    fig2.savefig(out / "fig1_1_elements.png", bbox_inches="tight")

    # --- key numbers ----------------------------------------------------------------------------------------------------
    s_w = ch01.mean_molecular_spacing(998.2, ch01.MOLAR_MASS["H2O"])
    s_a = ch01.mean_molecular_spacing(1.225, ch01.M_W_AIR)
    sig_n, tau_vec, tau_mag = ch01.traction_components((3.0, 0.0, -4.0), (0.0, 0.0, 1.0))
    print(f"fluid strain at unloading: {curves['Newtonian fluid (τt/μ)'][np.searchsorted(t, t_off) - 1]:.3f} rad; "
          f"solid strain while loaded: {tau / G:.3f} rad")
    print(f"mean molecular spacing: water {s_w * 1e9:.3f} nm, sea-level air {s_a * 1e9:.2f} nm (ratio {s_a / s_w:.1f})")
    print(f"traction (3, 0, -4) Pa on a face with normal +z: sigma_n = {sig_n:.1f} Pa (compression), |tau| = {tau_mag:.1f} Pa")
    print(f"saved fig1_1_strain_history.png, fig1_1_elements.png in {out}  ({time.perf_counter() - t_start:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
