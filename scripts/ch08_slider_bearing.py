"""§8.3 (C05–C07): the thin-gap scaling numbers (8.14)–(8.16), the lubrication profile (8.19) in a tilted-pad slider
bearing (Figs. 8.8–8.9, Example 8.1): exact pressure (corrected square), the O(α) form and the printed first-power
ghost, the steady Reynolds-equation quadrature, the load W(α) exact vs linear and the optimum taper (ours, San Andrés).

Run: ``.venv/Scripts/python.exe scripts/ch08_slider_bearing.py --no-show``
Figures → outputs/ch08/c07_slider_pressure.png, c07_slider_load.png, c06_gap_profiles.png.
"""
from __future__ import annotations

import numpy as np

from ch08_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch08_laminar_flow as ch08


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    rho_oil, nu_oil = 880.0, 4e-4  # 30-weight oil, ρ assumed (the book gives ν only)
    s = ch08.lubrication_scales(0.25, 1e-4, 10.0, rho_oil, rho_oil * nu_oil)
    print("oil gap 0.1 mm, surfaces 0.25 m, 10 m/s: " + ", ".join(f"{k} = {v:.4g}" for k, v in s.items()))
    h0, L, U, mu = 50e-6, 0.05, 5.0, 0.1
    for alpha in (0.2, 1.0, 1.5):
        st = ch08.slider_bearing_state(h0, alpha, L, U, mu)
        print(f"α = {alpha}: C₁ = {st['C1']:.4e} m²/s, p_max − p_e = {st['dp_max']:.4e} Pa ({st['p_max_atm']:.2f} atm) "
              f"at x/L = {st['x_pmax'] / L:.4f}, W exact = {st['W_exact']:.4e} N/m, W linear = {st['W_linear']:.4e} "
              f"(error {st['err_linear']:+.1f} %), inlet backflow = {st['inlet_backflow']}")
    x = np.linspace(0, L, 301)
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
    for alpha, c in ((0.2, COLORS["accent"]), (1.0, COLORS["teal"])):
        pe = ch08.slider_bearing(x, h0, alpha, L, U, mu, model="exact")
        pr, q = ch08.reynolds_pressure_1d(x, lambda xx, a_=alpha: h0 * (1 + a_ * xx / L), -U, 0.0, mu)
        print(f"α = {alpha}: Reynolds quadrature vs exact max rel diff = {np.max(np.abs(pr - pe)) / np.max(pe):.2e}, "
              f"q = {q:.4e} m²/s")
        ax[0].plot(x / L, pe / 1e6, color=c, lw=2.4, label=f"exact, α = {alpha}")
        ax[0].plot(x / L, ch08.slider_bearing(x, h0, alpha, L, U, mu, model="linear") / 1e6, color=c, ls="--", lw=1.2,
                   label=f"O(α), α = {alpha}")
        ax[0].plot(x / L, ch08.slider_bearing(x, h0, alpha, L, U, mu, model="book") / 1e6, color=COLORS["rose"],
                   ls=":", lw=1.2, label="printed (first power)" if alpha == 1.0 else None)
    ax[0].set_xlabel("x/L")
    ax[0].set_ylabel("p − p_e [MPa]")
    ax[0].legend(fontsize=7)
    ax[0].set_title(f"pressure under the pad, h₀ = {h0*1e6:.0f} μm, L = {L*100:.0f} cm, U = {U} m/s, μ = {mu} Pa s")
    al = np.linspace(0.01, 4.0, 200)
    We = [ch08.slider_bearing_load(1, a_, 1, 1, 1, "exact") for a_ in al]
    ax[1].plot(al, We, color=COLORS["accent"], lw=2.4, label="exact W h₀²/(μUL²) (ours)")
    ax[1].plot(al, al / 2, color=COLORS["muted"], ls="--", label="linear α/2 (book)")
    opt = ch08.slider_optimum_taper()
    ax[1].plot(opt["alpha_opt"], opt["W_coefficient"], "o", color=COLORS["orange"],
               label=f"optimum α = {opt['alpha_opt']:.4f} (K = {opt['K_opt']:.4f})")
    ax[1].set_ylim(0, 0.3)
    ax[1].set_xlabel("taper α")
    ax[1].set_ylabel("W h₀²/(μUL²)")
    ax[1].legend(fontsize=8)
    ax[1].set_title("load per unit width")
    print(f"optimum taper: α = {opt['alpha_opt']:.5f}, K = h_in/h_out = {opt['K_opt']:.5f}, "
          f"W h₀²/(6μUL²) = {opt['W_star']:.5f}")
    save(fig, out, "c07_slider_pressure")
    # gap profiles (Fig. 8.8 idea) in the pad frame
    fig2, a2 = plt.subplots(figsize=(12, 3.8))
    alpha = 1.5
    for xs in np.linspace(0.05, 0.95, 7) * L:
        hx = h0 * (1 + alpha * xs / L)
        yy = np.linspace(0, hx, 60)
        up = ch08.slider_gap_velocity(xs, yy, h0, alpha, L, U, mu, frame="pad")
        a2.plot(xs / L + 0.02 * up / U, yy * 1e6, color=COLORS["accent"], lw=1.6)
        a2.axvline(xs / L, color=COLORS["grid"], lw=0.8)
    a2.plot(x / L, h0 * (1 + alpha * x / L) * 1e6, color=COLORS["ink"], lw=2)
    a2.axhline(0, color=COLORS["ink"], lw=2)
    a2.set_xlabel("x/L (profiles u − U in the pad frame, scaled)")
    a2.set_ylabel("y [μm]")
    a2.set_title(f"(8.19) Poiseuille + Couette in the local gap, α = {alpha} (flow reversal near the inlet pad)")
    save(fig2, out, "c06_gap_profiles")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
