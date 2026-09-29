"""§8.2 (C02): plane Couette–Poiseuille flow — the four cases of Fig. 8.4 drawn from (8.5) with their signed shear
stress, the backflow threshold dp/dx = 2μU/h² (ours), flow rates, parity with ch04, and a sketch of the entrance region
of Fig. 8.2 with the wall-layer edge estimated by the Stokes-first thickness δ₉₉(x/U) (our estimate, labelled).

Run: ``.venv/Scripts/python.exe scripts/ch08_parallel_flows.py --no-show``
Figures → outputs/ch08/c02_couette_poiseuille.png, c02_entrance_sketch.png.
"""
from __future__ import annotations

import numpy as np

from ch08_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch04_conservation_laws as ch04
from fluidpy import ch08_laminar_flow as ch08
from fluidpy.core import navier_stokes as NS


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    h, U, mu = 0.01, 0.1, 1e-3  # 1 cm gap, 10 cm/s plate, water
    thr = float(ch08.channel_backflow_threshold(U, h, mu))
    print(f"backflow threshold dp/dx* = 2μU/h² = {thr:.4g} Pa/m")
    y = np.linspace(0.0, h, 201)
    cases = [("(a) U > 0, dp/dx < 0", U, -2 * thr), ("(b) U > 0, dp/dx > 0", U, 2 * thr),
             ("(c) plane Couette", U, 0.0), ("(d) plane Poiseuille", 0.0, -2 * thr)]
    fig, ax = plt.subplots(2, 4, figsize=(15, 7), sharey=True)
    for j, (lab, Uc, g) in enumerate(cases):
        u = ch08.channel_flow(y, h, Uc, g, mu)
        tau = ch08.channel_shear_stress(y, h, Uc, g, mu)
        st = ch08.couette_poiseuille_state(h, Uc, g, mu)
        ref = NS.exact_solution("couette", np.stack([0 * y, y]), U=Uc, h=h, G=-g, mu=mu)[0][0]
        print(f"{lab}: Q = {st['Q']:.4e} m²/s, V = {st['V']:.4e} m/s, τ(0) = {st['tau_bottom']:.4g} Pa, "
              f"τ(h) = {st['tau_top']:.4g} Pa, backflow = {st['backflow']}, parity with ch04 |Δu| = "
              f"{np.max(np.abs(u - ref)):.1e}")
        a = ax[0, j]
        a.plot(Uc * y / h * 100, y * 1e3, color=COLORS["muted"], ls="--", lw=1, label="Couette part")
        a.plot((u - Uc * y / h) * 100, y * 1e3, color=COLORS["teal"], ls=":", lw=1.2, label="Poiseuille part")
        a.plot(u * 100, y * 1e3, color=COLORS["accent"], lw=2.4, label="u(y), (8.5)")
        a.fill_betweenx(y * 1e3, 0, np.minimum(u, 0) * 100, color=COLORS["rose"], alpha=0.3)
        a.axvline(0, color=COLORS["ink"], lw=0.6)
        a.set_title(lab, fontsize=10)
        a.set_xlabel("u [cm/s]")
        ax[1, j].plot(tau * 1e3, y * 1e3, color=COLORS["orange"], lw=2)
        ax[1, j].axvline(0, color=COLORS["ink"], lw=0.6)
        ax[1, j].set_xlabel("τ = μ du/dy [mPa]")
    ax[0, 0].set_ylabel("y [mm]")
    ax[1, 0].set_ylabel("y [mm]")
    ax[0, 0].legend(fontsize=7, loc="upper left")
    fig.suptitle(f"Couette–Poiseuille flow, h = {h*1e3:.0f} mm, U = {U*100:.0f} cm/s, μ = {mu} Pa s "
                 f"(|dp/dx| = 2× the backflow threshold {thr:.3g} Pa/m)")
    save(fig, out, "c02_couette_poiseuille")
    # plane Poiseuille parity with ch04.plane_poiseuille (G = −dp/dx)
    d = np.max(np.abs(ch08.channel_flow(y, h, 0.0, -100.0, mu) - ch04.plane_poiseuille(y, 100.0, h, mu)))
    print(f"parity channel_flow(dpdx=−100) vs ch04.plane_poiseuille(G=100): {d:.1e}")
    # onset of backflow by bisection on the sampled profile
    lo, hi = 0.0, 10 * thr
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if np.min(ch08.channel_flow(y[1:-1], h, U, mid, mu)) >= 0 else (lo, mid)
    print(f"bisection onset of u < 0: {hi:.6g} Pa/m (threshold {thr:.6g})")
    # entrance sketch (Fig. 8.2 remake with our estimate)
    fig2, a2 = plt.subplots(figsize=(10, 3.2))
    Um, nu = 0.02, 1e-6
    x = np.linspace(0, 3.0, 300)
    dlt = np.minimum(ch08.diffusion_thickness(np.maximum(x, 1e-9) / Um, nu), h / 2)
    a2.plot(x, dlt * 1e3, color=COLORS["teal"], lw=2, label="wall-layer edge, δ₉₉(t = x/U) (our estimate)")
    a2.plot(x, (h - dlt) * 1e3, color=COLORS["teal"], lw=2)
    a2.axhline(0, color=COLORS["ink"], lw=2)
    a2.axhline(h * 1e3, color=COLORS["ink"], lw=2)
    xm = x[np.argmax(dlt >= h / 2)]
    a2.axvline(xm, color=COLORS["orange"], ls="--", lw=1.2, label=f"layers merge ≈ {xm:.2f} m: fully developed after")
    for xs in (0.2, 1.0, 2.6):
        yy = np.linspace(0, h, 60)
        dd = float(np.minimum(ch08.diffusion_thickness(xs / Um, nu), h / 2))
        prof = np.clip(np.minimum(yy, h - yy) / dd, 0, 1) if xs < xm else 4 * yy * (h - yy) / h ** 2
        a2.plot(xs + 0.15 * prof, yy * 1e3, color=COLORS["accent"], lw=1.5)
    a2.set_xlabel("x [m]")
    a2.set_ylabel("y [mm]")
    a2.legend(fontsize=8, loc="upper right")
    a2.set_title(f"Developing flow in a channel (Fig. 8.2 idea): U = {Um} m/s, h = {h*1e3:.0f} mm, water")
    print(f"entrance estimate: wall layers (δ₉₉ = 2 erfc⁻¹(0.01)√(νx/U)) meet at x ≈ {xm:.3f} m")
    save(fig2, out, "c02_entrance_sketch")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
