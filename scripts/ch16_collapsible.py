"""Chapter 16, §16.3 — collapsible tubes and Shapiro's one-dimensional steady theory (16.183)–(16.206): the tube law
(16.188) (Fig. 16.19 redrawn with the two fitted forms only), the wave speed (16.189), the speed index S = u/c
(16.195), the area and speed-index equations (16.204)–(16.206) with their 1 − S² factor, the corrected differential
(16.201) (slip S23, sympy), a smooth throat transition through S = 1 (Fig. 16.20(a)), flow limitation in a Starling
resistor (Case 1), and a giraffe-jugular-like vertical tube (ours).

Run: ``.venv/Scripts/python.exe scripts/ch16_collapsible.py --no-show``   Figures -> outputs/ch16/collapsible_*.png
"""
from __future__ import annotations

import numpy as np
from ch16_common import C_E, C_P, C_REF, C_U, C_WALL, Timer, finish, parse_args, save, setup

from fluidpy import ch16_biofluid_mechanics as ch16


def euler_march(Q, alpha_in, x0, x1, n, prm):
    """From-scratch explicit Euler for dα/dx of (16.204) (the notebook's transparent version); stops near S = 1."""
    x = np.linspace(x0, x1, n + 1)
    al = np.empty_like(x)
    al[0] = alpha_in
    for i in range(n):
        tr = ch16.collapsible_terms(al[i], Q, prm, x=x[i])
        if tr["S2"] > 1 - 1e-3:
            return x[: i + 1], al[: i + 1]
        al[i + 1] = al[i] + (x[i + 1] - x[i]) * tr["dalpha_dx"]
    return x, al


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    I = ch16.illustrative_inputs()
    rho, Kp, n = I["rho_water"], I["Kp_collapse"], 1.5
    A0 = np.pi * I["a_collapse"] ** 2
    al = np.linspace(0.15, 1.0, 300)  # (16.188) is fitted to the collapse branch only: α ≤ 1 (Fig. 16.19)
    P1, dP, d2P = ch16.collapsible_tube_law_derivatives(al, n)
    print(f"tube law (16.188), n = 3/2: M = 3 + αP″/P′ (16.206) = {ch16.shapiro_M(0.5):.3f} = 2 − n "
          f"(finite differences {ch16.shapiro_M(0.5, tube_law=lambda a: 1 - a ** -1.5):.6f})")
    c_fd = ch16.bramwell_hill_speed(0.6 * A0, Kp * (ch16.collapsible_tube_law(0.6 + 1e-6) - ch16.collapsible_tube_law(0.6 - 1e-6))
                                    / (2e-6 * A0), rho)
    print(f"  wave speed at α = 0.6 (16.189): c = √(nK_pα⁻ⁿ/ρ) = {ch16.collapsible_wave_speed(0.6, Kp, rho):.5f} m/s; "
          f"(16.185) by finite differences {c_fd:.5f} m/s")
    sym = ch16.collapsible_differentials_sympy()
    symp = ch16.collapsible_differentials_sympy(printed_16_201=True)
    print(f"  sympy: corrected (16.201) − d(c²)/c² = {sym['dc_difference']}; printed (16.201) leaves {symp['dc_difference']} "
          f"(slip S23); (16.204) residual {sym['eq204_residual']}, (16.205) residual {sym['eq205_residual']}; M(power law) = {sym['M_power_law']}")
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
    ax[0].plot(al, -ch16.collapsible_tube_law(al, n=n), color=C_WALL, label="−P = α⁻ⁿ − 1 (16.188)")
    ax[0].plot(al, -ch16.collapsible_tube_law(al, n=n, form="pure"), color=C_REF, ls="--", label="−P = α⁻ⁿ")
    ax[0].axhline(0, color=C_REF, lw=0.6)
    ax[0].axvline(1.0, color=C_REF, lw=0.6, ls=":")
    ax[0].text(0.98, 9.0, "fit valid for α ≤ 1\n(collapse branch)", ha="right", fontsize=8)
    ax[0].set(xlabel="area ratio α = A/A₀", ylabel="−P = −(p − p_e)/K_p", xlim=(0.1, 1.05), ylim=(-0.5, 18),
              title="collapsible tube law (ours, Fig. 16.19 style)")
    ax[0].legend()
    # smooth throat transition (frictionless): S = 1 exactly at the A₀ minimum.  The whole subcritical branch stays on
    # the fitted collapse branch (α ≤ 1) iff H = K_pP(α_t) + ½ρc_t² = K_p[1 − (1 − n/2)α_t⁻ⁿ] ≤ 0, i.e.
    # α_t ≤ (1 − n/2)^{1/n} = 0.25^{2/3} = 0.397 for n = 3/2; we take α_t = 0.35.
    al_t = 0.35
    al_t_max = (1.0 - n / 2.0) ** (1.0 / n)
    Qs = ch16.choking_flux(al_t, A0, Kp, rho, n)
    H = Kp * ch16.collapsible_tube_law(al_t, n=n) + 0.5 * rho * (Qs / (al_t * A0)) ** 2
    a_far = ch16.alpha_from_tube_law(H / Kp, n)  # far upstream: u → 0, K_pP(α) = H
    ell = 0.05
    A0f = lambda x: A0 * (1 + (x / ell) ** 2)  # noqa: E731
    prm = dict(A0=A0f, Kp=Kp, n=n, rho=rho, mu=I["mu_water"], friction=False)
    x_in = -0.12
    a_in = ch16.collapsible_frictionless_alpha(Qs, H, A0f(x_in), Kp, rho, n, branch="sub")
    m = ch16.collapsible_march(Qs, a_in, (x_in, 0.12), prm, rtol=1e-11, atol=1e-13)
    xb = np.linspace(x_in, 0.12, 121)
    sub = [ch16.collapsible_frictionless_alpha(Qs, H, A0f(x), Kp, rho, n, branch="sub") for x in xb]
    sup = [ch16.collapsible_frictionless_alpha(Qs, H, A0f(x), Kp, rho, n, branch="super") for x in xb]
    xe, ae = euler_march(Qs, a_in, x_in, 0.12, 4000, prm)
    far = m["x"] < -0.02  # away from the square-root behaviour of α(x) at the critical point
    ber_m = np.array([ch16.collapsible_frictionless_alpha(Qs, H, A0f(x), Kp, rho, n) for x in m["x"]])
    fe = xe < -0.02
    ber_e = np.array([ch16.collapsible_frictionless_alpha(Qs, H, A0f(x), Kp, rho, n) for x in xe[fe]])
    print(f"throat mode (ours): α_t = {al_t} (≤ α_t,max = (1 − n/2)^(1/n) = {al_t_max:.4f}), "
          f"H = K_pP + ½ρu² = {H:.2f} Pa ≤ 0; upstream α: {a_in:.4f} at x = {x_in * 100:.0f} cm, "
          f"{a_far:.4f} far upstream (all ≤ 1); max α of the subcritical branch shown = {np.nanmax(sub):.4f}")
    print(f"throat (ours): critical flux Q* = α_tA₀c(α_t) = {Qs * 1e6:.3f} cm³/s at α_t = {al_t}; the march stops at "
          f"x = {m['x_choke']:.2e} m (throat at 0) with S = {m['S'][-1]:.6f}; march vs Bernoulli: "
          f"{np.max(np.abs(ber_m[far] - m['alpha'][far])):.1e} for x < −2 cm, {np.max(np.abs(ber_m - m['alpha'])):.1e} up to "
          f"the throat; explicit Euler (4000 steps) vs Bernoulli for x < −2 cm: {np.max(np.abs(ae[fe] - ber_e)):.1e}")
    ax[1].plot(xb * 100, sub, color=C_U, label="subcritical branch (S < 1)")
    ax[1].plot(xb * 100, sup, color=C_P, label="supercritical branch (S > 1)")
    ax[1].plot(m["x"] * 100, m["alpha"], color=C_E, ls="--", lw=2, label="march of (16.204) to S = 1")
    ax[1].plot(xb * 100, A0f(xb) / A0, color=C_REF, ls=":", label="A₀(x)/A₀,min")
    ax[1].set(xlabel="x [cm]", ylabel="α = A/A₀(x)", ylim=(0, 1.6), title="smooth transition at the throat (Fig. 16.20a)")
    ax[1].legend(fontsize=8)
    # Starling resistor: flow limitation
    prmS = dict(L=I["L_collapse"], A0=A0, Kp=Kp, n=n, rho=rho, mu=I["mu_water"])
    with Timer("flow-limitation curve"):
        p2 = np.linspace(-700.0, -5.0, 10 if args.fast else 24)
        cur = ch16.flow_limitation_curve(p2, p1=0.0, pe=0.0, prm=prmS)
    print(f"Starling resistor (ours: L = {prmS['L']} m, a₀ = {I['a_collapse'] * 1e3} mm, K_p = {Kp} Pa, p₁ = p_e): "
          f"Q_choke = {cur['Q_choke'] * 1e6:.3f} cm³/s at p₂ ≤ {cur['p2_choke']:.1f} Pa; inlet choking flux "
          f"{ch16.choking_flux(1.0, A0, Kp, rho) * 1e6:.3f} cm³/s")
    ax[2].plot(-p2, cur["Q"] * 1e6, "o-", color=C_U, ms=3)
    ax[2].axhline(cur["Q_choke"] * 1e6, color=C_REF, ls=":", label="flow limited: S = 1 at the outlet")
    ax[2].set(xlabel="p₁ − p₂ [Pa]", ylabel="Q [cm³/s]", title="Starling resistor, Case 1 (ours)")
    ax[2].legend()
    save(fig, out, "collapsible")
    st = ch16.starling_state(0.0, -100.0, 0.0, prmS)
    print(f"  p₂ = −100 Pa: Q = {st['Q'] * 1e6:.4f} cm³/s, α_out = {st['alpha_out']:.4f}, S_out = {st['S_out']:.4f}, limited {st['limited']}")
    # friction sign table and the giraffe-jugular-like vertical tube (ours)
    for S_target in (0.3, 2.0):
        a_s = 0.8
        Qx = S_target * ch16.choking_flux(a_s, A0, Kp, rho)
        tr = ch16.collapsible_terms(a_s, Qx, prmS)
        print(f"  friction only, S = {np.sqrt(tr['S2']):.2f}: dα/dx = {tr['dalpha_dx']:+.4f} 1/m "
              f"({'area falls' if tr['dalpha_dx'] < 0 else 'area rises'} downstream)")
    Qg = 0.3 * ch16.choking_flux(1.0, A0, Kp, rho)
    prmS = dict(prmS, mu=I["mu"], rho=I["rho"])  # blood in a vein
    from scipy.optimize import brentq

    def balance(a_):
        return ch16.collapsible_terms(a_, Qg, prmS, dzdx=-1.0)["driver"]

    a_eq = brentq(balance, 1e-4, 1.0)  # bracket kept on the fitted collapse branch α ≤ 1 of (16.188)
    tr = ch16.collapsible_terms(a_eq, Qg, prmS, dzdx=-1.0)
    print(f"  vertical tube draining down (dz/dx = −1), Q = {Qg * 1e6:.2f} cm³/s: uniform state where friction balances gravity "
          f"at α = {a_eq:.4f} (S = {np.sqrt(tr['S2']):.3f}, p − p_e = {Kp * tr['P']:.1f} Pa: partly collapsed)")
    ic = ch16.collapsible_influence_coefficients(0.25, 0.8)
    print("  influence coefficients at S² = 0.25, α = 0.8 (dα/α, dS²/S² per unit driver): "
          + "; ".join(f"{k}: ({v[0]:+.3f}, {v[1]:+.3f})" for k, v in ic.items()))
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
