"""§7.6 (C11): the hydraulic jump — the conjugate depth (7.81) from the CV momentum balance (7.80) (Fig. 7.20b remake
with the momentum budget), the energy loss E₂ − E₁ = −g(H₂ − H₁)³/(4H₁H₂) that forbids Fr₁ < 1, the from-scratch
quadratic root, and parity with Example 4.3's bore speed in the moving frame (Fig. 7.20c).

Run: ``.venv/Scripts/python.exe scripts/ch07_hydraulic_jump.py --no-show``
Figures → outputs/ch07/c11_hydraulic_jump.png.
"""
from __future__ import annotations

import numpy as np

from ch07_drawings import cv_box, parse_args, save, setup

from fluidpy import ch04_conservation_laws as ch04
from fluidpy import ch07_gravity_waves as ch07
from fluidpy.core.style import COLORS


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    H1, Fr1 = 0.1, 3.0
    j = ch07.jump_state(H1, Fr1)
    print(f"H₁ = {H1} m, Fr₁ = {Fr1}: u₁ = {j['u1']:.4f} m/s, H₂ = {j['H2']:.4f} m (ratio {j['ratio']:.4f}),"
          f" u₂ = {j['u2']:.4f} m/s, Fr₂ = {j['Fr2']:.4f}; head loss {j['head_loss']:.4f} m, power lost"
          f" {j['power_loss']:.0f} W per metre of width")
    print(f"  momentum budget [N/m]: in {j['mom_in']:.2f}, out {j['mom_out']:.2f}, p₁ {j['p_in']:.2f}, p₂ {j['p_out']:.2f};"
          f" residual {j['residual']:.1e}; (7.80) residual {ch07.jump_momentum_residual(H1, j['H2'], j['Q']):.1e}")
    r = np.roots([1.0, 1.0, -2 * Fr1 ** 2])
    print(f"  from scratch: roots of r² + r − 2Fr₁² = {np.round(r, 6)}; positive root {r.max():.10f} vs (7.81) {j['ratio']:.10f}")
    print(f"  moving bore into H₁: speed u₁ = {j['bore_speed']:.6f} m/s; ch04.bore_speed(H₁, H₂) = "
          f"{ch04.bore_speed(H1, j['H2']):.6f} m/s; flow behind {j['flow_behind']:.4f} m/s")
    jb = ch07.hydraulic_jump(H1, Fr1=0.7)
    print(f"  Fr₁ = 0.7 would give H₂/H₁ = {jb['ratio']:.4f} and E₂ − E₁ = {jb['dE']:+.4f} J/kg > 0: allowed = {jb['allowed']}")
    assert abs(r.max() - j["ratio"]) < 1e-12 and j["dE"] < 0 < jb["dE"]

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    cv_box(ax[0], H1, j["H2"], 0.0, 1.0, 0.5, 0.04)
    ax[0].annotate("", xy=(0.25, 0.5 * H1), xytext=(0.08, 0.5 * H1),
                   arrowprops=dict(arrowstyle="->", color=COLORS["teal"], lw=2))
    ax[0].annotate("", xy=(0.95, 0.5 * j["H2"]), xytext=(0.8, 0.5 * j["H2"]),
                   arrowprops=dict(arrowstyle="->", color=COLORS["teal"], lw=1))
    ax[0].text(0.1, 0.62 * H1, f"u₁ = {j['u1']:.2f} m/s", fontsize=8)
    ax[0].text(0.72, 0.62 * j["H2"], f"u₂ = {j['u2']:.2f} m/s", fontsize=8)
    ax[0].set_xlabel("x [m]")
    ax[0].set_ylabel("depth [m]")
    ax[0].set_title(f"stationary jump, Fr₁ = {Fr1} (Fig. 7.20b remake)", fontsize=10)
    names = ["ρQu₁ in", "−ρQu₂ out", "½ρgH₁² face 1", "−½ρgH₂² face 2"]
    vals = [j["mom_in"], -j["mom_out"], j["p_in"], -j["p_out"]]
    ax[1].barh(names, vals, color=[COLORS["teal"], COLORS["teal"], COLORS["amber"], COLORS["amber"]])
    ax[1].axvline(0, color=COLORS["ink"], lw=0.8)
    ax[1].set_xlabel("x-momentum budget per metre of width [N/m]")
    ax[1].set_title(f"(7.80): the four terms sum to {j['residual']:.0e}", fontsize=10)
    F = np.linspace(0.3, 6, 300)
    ratio = [ch07.hydraulic_jump(1.0, Fr1=f)["ratio"] for f in F]
    loss = [ch07.hydraulic_jump(1.0, Fr1=f)["head_loss"] for f in F]
    ax[2].plot(F, ratio, color=COLORS["accent"], label="H₂/H₁ (7.81)")
    ax[2].plot(F, loss, color=COLORS["rose"], label="head loss/H₁ = (H₂ − H₁)³/(4H₁²H₂)")
    ax[2].axvspan(0.3, 1.0, color=COLORS["rose"], alpha=0.1, label="Fr₁ < 1: would create energy")
    ax[2].axhline(0, color=COLORS["ink"], lw=0.6)
    ax[2].set_xlabel("upstream Froude number Fr₁ (4.104)")
    ax[2].set_ylabel("[–]")
    ax[2].set_title("momentum sets H₂; energy sets the direction", fontsize=10)
    ax[2].legend(fontsize=8)
    save(fig, out, "c11_hydraulic_jump")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
