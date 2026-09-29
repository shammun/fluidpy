"""Elementary lubrication theory: thin-gap scaling, the lubrication velocity profile and flux, the steady 1-D Reynolds
pressure equation, the tilted-pad slider bearing, Hele-Shaw flow and the gravity spreading of a viscous film.

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e (2012), Ch. 8, §8.3, Eqs. (8.13)–(8.19), Examples 8.1 (slider
bearing), 8.2 (Hele-Shaw 1898), 8.3 (spreading bead) and §8.4 Example 8.7 (its similarity form). Transcribed from the
page images chapters/pages/ch08/p345–p353, p363. Reused by Ch. 9 (anisotropic scaling template), Ch. 13 (viscous
gravity currents, depth-averaged thin layers), Ch. 16 (joint and mucus films).

Conventions (analysis ch08 §9)
* Gap problems: x along the gap, y across it, gap h(x, t); lower wall speed U_0, upper wall speed U_h. **Hele-Shaw
  switches to z across the gap** (plates at z = 0 and h) with (x, y) in the plane; φ there is a velocity potential.
* Pressure gradients are dp/dx [Pa/m] (book sign). Thin-film spreading: y up from the plate, g downward.
* Book slips coded in their corrected form (the printed form kept as a labelled option for wrong-variant tests):
  R7 (8.17a) and Example 8.2 print the lubrication balance without ν (dimensionally wrong); R8 Example 8.1's
  intermediate integrals use (1 − αx/L) and its exact p(x) prints (1 + αx/L) to the first power — sympy gives the
  square (``slider_bearing(model="book")`` is the printed form); R8b (8.19) adds U_0 on top of the Couette part so
  u(h) = U_h + U_0 — ``lubrication_velocity(form="consistent")`` (default) has U_h y/h + U_0(1 − y/h).
* Defaults: μ = 1e-3 Pa s (water) for gap flows; the spreading bead uses ρ = 1000 kg/m³, μ = 1 Pa s, g = G0.
"""
from __future__ import annotations

from typing import Callable

import numpy as np
from scipy.integrate import cumulative_trapezoid, quad
from scipy.linalg import solve_banded
from scipy.optimize import minimize_scalar
from scipy.special import gamma as _gamma

from ._util import as_scalar_if_0d, require_positive
from .thermo import G0, P_ATM

__all__ = ["lubrication_scales", "lubrication_term_magnitudes", "lubrication_velocity", "lubrication_flux",
           "reynolds_pressure_1d", "slider_bearing", "slider_bearing_load", "slider_optimum_taper",
           "slider_bearing_state", "slider_gap_velocity", "hele_shaw_velocity", "hele_shaw_potential",
           "hele_shaw_mean_velocity", "hele_shaw_cylinder", "hele_shaw_streamfunction_grid", "thin_film_velocity",
           "thin_film_flux", "thin_film_spread", "viscous_current_similarity", "viscous_current_eta_N",
           "thin_film_state"]

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d


# ======================================================================================================================
# scaling (8.14)–(8.17)
# ======================================================================================================================
def lubrication_scales(L, h, U, rho, mu, p_a: float = P_ATM) -> dict:
    """Dimensionless groups of the thin-gap scaling.

    Book: §8.3, Eq. (8.14) (x* = x/L, y* = y/h, t* = Ut/L, u* = u/U, v* = v/εU, p* = p/P_a) and (8.16a,b): ε = h/L
    (fineness ratio), Re_L = ρUL/μ, the reduced Reynolds number ε²Re_L and the bearing number Λ = μUL/(P_a h²).
    Analysis §9 R24: in the book's own oil example Λ ≈ 10³, not "near unity"; the natural (viscous) pressure scale is
    μUL/h², the one that makes Λ = 1 — returned as ``p_visc``.

    Parameters
    ----------
    L : passage length [m];  h : gap [m];  U : speed scale [m/s];  rho : ρ [kg/m³];  mu : μ [Pa s];
    p_a : pressure scale of (8.14) [Pa] (atmospheric).

    Returns
    -------
    dict of floats: eps, Re_L, eps2_Re_L, eps4_Re_L, Lambda, p_visc [Pa].

    Validation — tests/test_ch08.py: test_lubrication_balance_V2_derivation_and_units,
      test_lubrication_scales_V1_definitions_and_term_magnitudes,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable,
      test_book_V6_section_8_3_lubrication_forms_and_numbers.
    Checks: V1 definitions; V6 the book's ε²Re_L for its oil example (private JSON). Label: analytic.
    """
    L, h, U, rho, mu = float(L), float(h), float(U), float(rho), float(mu)
    eps = h / L
    Re = rho * U * L / mu
    return dict(eps=eps, Re_L=Re, eps2_Re_L=eps ** 2 * Re, eps4_Re_L=eps ** 4 * Re,
                Lambda=mu * U * L / (float(p_a) * h ** 2), p_visc=mu * U * L / h ** 2)


def lubrication_term_magnitudes(L, h, U, rho, mu, p_scale: str = "viscous", p_a: float = P_ATM) -> dict:
    """Coefficient of every term of the scaled momentum equations (8.16a) and (8.16b).

    Book: §8.3, Eqs. (8.16a,b) — x-momentum: ε²Re_L (inertia) = −(1/Λ)∂p*/∂x* + ε²∂²u*/∂x*² + ∂²u*/∂y*²; y-momentum:
    ε⁴Re_L (inertia) = −(1/Λ)∂p*/∂y* + ε⁴∂²v*/∂x*² + ε²∂²v*/∂y*². The pressure coefficient depends on the pressure scale:
    ``"atm"`` (alias ``"atmospheric"``) P_a (the book's (8.14): 1/Λ), ``"viscous"`` μUL/h² (1) or ``"dynamic"`` ρU²
    (ε²Re_L).

    Parameters
    ----------
    L, h : [m];  U : [m/s];  rho : [kg/m³];  mu : [Pa s];  p_scale : "viscous" | "atm" | "dynamic";  p_a : [Pa].

    Returns
    -------
    dict of floats (design Part C 2.2 keys): x_inertia (ε²Re_L), x_pressure, x_diff_along (ε²), x_diff_across (1),
    y_inertia (ε⁴Re_L), y_pressure, y_diff_along (ε⁴), y_diff_across (ε²); the same four diffusion numbers also as
    x/y_diff_streamwise / _cross; plus eps, Re_L, Lambda.

    Validation — tests/test_ch08.py: test_lubrication_scales_V1_definitions_and_term_magnitudes,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V2 equals ``ch08.lubrication_nondim_sympy`` coefficients. Label: analytic.
    """
    s = lubrication_scales(L, h, U, rho, mu, p_a)
    eps, Re = s["eps"], s["Re_L"]
    P = {"atmospheric": float(p_a), "atm": float(p_a), "viscous": s["p_visc"],
         "dynamic": float(rho) * float(U) ** 2}.get(p_scale)
    if P is None:
        raise ValueError('p_scale must be "viscous", "atm" (alias "atmospheric") or "dynamic"')
    cp = P / s["p_visc"]  # = P h²/(μUL); 1/Λ for P = P_a
    return dict(x_inertia=eps ** 2 * Re, x_pressure=cp, x_diff_along=eps ** 2, x_diff_across=1.0,
                y_inertia=eps ** 4 * Re, y_pressure=cp, y_diff_along=eps ** 4, y_diff_across=eps ** 2,
                x_diff_streamwise=eps ** 2, x_diff_cross=1.0, y_diff_streamwise=eps ** 4, y_diff_cross=eps ** 2,
                eps=eps, Re_L=Re, Lambda=s["Lambda"])


# ======================================================================================================================
# the lubrication profile (8.18)–(8.19) and the flux
# ======================================================================================================================
def lubrication_velocity(y, h, dpdx, U_h: float = 0.0, U_0: float = 0.0, mu: float = 1e-3, form: str = "consistent"):
    """Velocity in a thin gap: the local Poiseuille + Couette profile of zeroth-order lubrication theory.

    Book: §8.3, Eq. (8.19) (from (8.18) with u = U_0 on y = 0 and u = U_h on y = h). ``form="consistent"`` (default)
    uses U_h y/h + U_0(1 − y/h), which meets both wall conditions; ``form="book"`` is the printed U_h y/h + U_0, which
    gives u(h) = U_h + U_0 (analysis §9 R8b; identical when U_0 = 0, as in every book example).

    Parameters
    ----------
    y : [m], 0 ≤ y ≤ h;  h : local gap [m];  dpdx : local ∂p/∂x [Pa/m];  U_h, U_0 : upper/lower wall speeds [m/s];
    mu : μ [Pa s];  form : "consistent" | "book".

    Returns
    -------
    u : [m/s].

    Validation — tests/test_ch08.py: test_lubrication_velocity_V1_walls_flux_and_printed_form,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable,
      test_book_V6_section_8_3_lubrication_forms_and_numbers.
    Checks: V1 u(0) = U_0, u(h) = U_h; U_0 = 0 reduces to (8.5); wrong variant ``form="book"``.
    Label: analytic.
    """
    y_, h_ = _F(y), _F(h)
    s = y_ / h_
    pois = -h_ ** 2 / (2.0 * float(mu)) * _F(dpdx) * s * (1.0 - s)
    if form == "consistent":
        # DEVIATION: U_0(1 − y/h) instead of the printed + U_0 — the printed form gives u(h) = U_h + U_0 (analysis §9 R8b)
        return _S(pois + float(U_h) * s + float(U_0) * (1.0 - s))  # Eq. (8.19), U_0 term made consistent
    if form == "book":
        return _S(pois + float(U_h) * s + float(U_0))  # Eq. (8.19) as printed
    raise ValueError('form must be "consistent" or "book"')


def lubrication_flux(h, dpdx, U_h: float = 0.0, U_0: float = 0.0, mu: float = 1e-3):
    """Volume flux per unit width through the gap, q = ∫₀ʰ u dy = −h³(∂p/∂x)/(12μ) + (U_0 + U_h)h/2.

    Book: §8.3 (integral of (8.19) — used in Example 8.1 as C₁ = ∫(u − U)dy; the general flux and the Reynolds equation
    ∂h/∂t + ∂q/∂x = 0 of Exercises 8.19–8.20 are not written in the text: ours).

    Parameters
    ----------
    h : [m];  dpdx : [Pa/m];  U_h, U_0 : [m/s];  mu : [Pa s].

    Returns
    -------
    q : [m²/s].

    Validation — tests/test_ch08.py: test_lubrication_velocity_V1_walls_flux_and_printed_form,
      test_reynolds_equation_V2_derivation, test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V1 ``quad`` of :func:`lubrication_velocity`. Label: analytic.
    """
    h_ = _F(h)
    return _S(-h_ ** 3 * _F(dpdx) / (12.0 * float(mu)) + (float(U_0) + float(U_h)) * h_ / 2.0)


def reynolds_pressure_1d(x, h, U_0: float = 0.0, U_h: float = 0.0, mu: float = 1e-3, p_left: float = 0.0,
                         p_right: float = 0.0):
    """Pressure under a steady 1-D lubricating film of arbitrary gap h(x) (steady Reynolds equation, ours).

    Book: §8.3 (the recipe after (8.19): profile + mass conservation + pressure boundary conditions; Example 8.1 is the
    special case h = h₀(1 + αx/L)). Steady ⇒ q constant, dp/dx = 6μ(U_0 + U_h)/h² − 12μq/h³; q follows from p(x₀) =
    p_left, p(x_end) = p_right. Must be used in the frame in which the gap shape is steady (for Example 8.1 the pad
    frame: U_0 = −U, U_h = 0, and then q = C₁).

    Parameters
    ----------
    x : array_like, shape (N,) — increasing positions [m].
    h : callable x ↦ h [m] (integrals by ``quad`` between nodes) or array of shape (N,) (``cumulative_trapezoid``,
        second order).
    U_0, U_h : wall speeds [m/s];  mu : μ [Pa s];  p_left, p_right : end pressures [Pa].

    Returns
    -------
    (p, q) : pressure at the nodes [Pa] and the constant flux per unit width [m²/s].

    Validation — tests/test_ch08.py: test_reynolds_pressure_1d_V1_exact_slider_and_V3_order_two,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V1 equals the exact slider p(x) to 1e-10 with a callable h; V3 order 2 with an array h.
    Label: analytic, converged.
    """
    x = _F(x)
    mu = float(mu)
    if callable(h):
        f2 = lambda s: 1.0 / float(h(s)) ** 2  # noqa: E731
        f3 = lambda s: 1.0 / float(h(s)) ** 3  # noqa: E731
        seg2 = [quad(f2, a, b, epsabs=0.0, epsrel=1e-13)[0] for a, b in zip(x[:-1], x[1:])]
        seg3 = [quad(f3, a, b, epsabs=0.0, epsrel=1e-13)[0] for a, b in zip(x[:-1], x[1:])]
        I2 = np.concatenate([[0.0], np.cumsum(seg2)])
        I3 = np.concatenate([[0.0], np.cumsum(seg3)])
    else:
        hh = _F(h)
        I2 = cumulative_trapezoid(1.0 / hh ** 2, x, initial=0.0)
        I3 = cumulative_trapezoid(1.0 / hh ** 3, x, initial=0.0)
    S = float(U_0) + float(U_h)
    q = (float(p_left) - float(p_right) + 6.0 * mu * S * I2[-1]) / (12.0 * mu * I3[-1])
    p = float(p_left) + 6.0 * mu * S * I2 - 12.0 * mu * q * I3  # ∫ dp/dx, dp/dx = 6μ(U_0+U_h)/h² − 12μq/h³
    return p, float(q)


# ======================================================================================================================
# Example 8.1 slider bearing
# ======================================================================================================================
def slider_bearing(x, h0, alpha, L, U, mu: float = 1e-3, p_e: float = 0.0, model: str = "exact"):
    """Pressure under a tilted bearing pad (gap h = h₀(1 + αx/L)) moving at U over a flat surface.

    Book: §8.3, Example 8.1. ``model="exact"``: p − p_e = (6μLU/h₀²) α(x/L)(1 − x/L)/[(2 + α)(1 + αx/L)²] — sympy
    (analysis §9 R8; the book prints the denominator (1 + αx/L) to the first power, ``model="book"``, which does not
    satisfy dp/dx = −12μC₁/h³ − 6μU/h²); ``model="linear"``: p − p_e ≅ (3αμLU/h₀²)(x/L)(1 − x/L) (the book's O(α)
    result, exact for α → 0).

    Parameters
    ----------
    x : [m], 0 ≤ x ≤ L;  h0 : gap at x = 0 [m];  alpha : taper α [–] (h(L) = h₀(1 + α));  L : pad length [m];
    U : pad speed [m/s] (+x);  mu : μ [Pa s];  p_e : exterior pressure [Pa];  model : "exact" | "linear" | "book".

    Returns
    -------
    p : [Pa].

    Validation — tests/test_ch08.py: test_reynolds_pressure_1d_V1_exact_slider_and_V3_order_two,
      test_slider_bearing_V1_quadrature_load_and_numbers, test_slider_bearing_V7_reversal_and_ends,
      test_slider_bearing_state_V1_explainer_numbers,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable,
      test_book_V6_section_8_3_lubrication_forms_and_numbers.
    Checks: V2 sympy (ODE and p(0) = p(L) = p_e); V1 = :func:`reynolds_pressure_1d` with h callable (pad
    frame); wrong variant "book". Label: symbolic, analytic.
    """
    s = _F(x) / float(L)
    a, c = float(alpha), 6.0 * float(mu) * float(L) * float(U) / float(h0) ** 2
    if model == "exact":
        # DEVIATION: denominator (1 + αx/L)² — the page prints the first power, which violates dp/dx = −12μC₁/h³ − 6μU/h² (sympy, analysis §9 R8)
        dp = c * a * s * (1.0 - s) / ((2.0 + a) * (1.0 + a * s) ** 2)  # Example 8.1, exact (corrected square)
    elif model == "linear":
        dp = 0.5 * c * a * s * (1.0 - s)  # Example 8.1, O(α): 3αμLU/h₀² (x/L)(1 − x/L)
    elif model == "book":
        dp = c * a * s * (1.0 - s) / ((2.0 + a) * (1.0 + a * s))  # Example 8.1 as printed (first power)
    else:
        raise ValueError('model must be "exact", "linear" or "book"')
    return _S(float(p_e) + dp)


def _load_series(a: float) -> float:
    # W h0²/(μUL²) = α/2 − 3α²/4 + 33α³/40 − 13α⁴/16 + 171α⁵/224 + … (sympy series of the exact load)
    return a / 2.0 - 3.0 * a ** 2 / 4.0 + 33.0 * a ** 3 / 40.0 - 13.0 * a ** 4 / 16.0 + 171.0 * a ** 5 / 224.0


def slider_bearing_load(h0, alpha, L, U, mu: float = 1e-3, model: str = "exact"):
    """Load per unit width carried by the slider bearing, W = ∫₀ᴸ (p − p_e) dx.

    Book: §8.3, Example 8.1: W = αμL²U/(2h₀²) (``model="linear"``, the book's result). ``model="exact"`` (ours, analysis
    D19): W = (6μUL²/(h₀²α²))[ln(1 + α) − 2α/(2 + α)], evaluated by its Taylor series for |α| < 1e-3 (the closed form
    cancels); equal to San Andrés' W(K) with K = 1 + α (inlet/exit gap ratio). ``model="book"``: the integral of the
    printed (first-power) pressure, by ``quad`` (for the ghost curve).

    Parameters
    ----------
    h0 : [m];  alpha : [–] (> −1);  L : [m];  U : [m/s];  mu : [Pa s];  model : "exact" | "linear" | "book".

    Returns
    -------
    W : [N/m] (negative when αU < 0: the pad is sucked down).

    Validation — tests/test_ch08.py: test_slider_bearing_V1_quadrature_load_and_numbers,
      test_slider_bearing_V7_reversal_and_ends, test_slider_optimum_taper_V5_san_andres,
      test_slider_bearing_state_V1_explainer_numbers,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable,
      test_book_V6_section_8_3_lubrication_forms_and_numbers.
    Checks: V1 ``quad`` of :func:`slider_bearing`; V2 series → linear; V5 San Andrés α_opt, W*.
    Label: analytic, symbolic, benchmark.
    """
    a = float(alpha)
    if a <= -1.0:
        raise ValueError("alpha must exceed −1 (positive gap)")
    c = float(mu) * float(U) * float(L) ** 2 / float(h0) ** 2
    if model == "linear":
        return c * a / 2.0  # W = αμL²U/(2h₀²)
    if model == "exact":
        if abs(a) < 1e-3:
            return c * _load_series(a)
        return c * 6.0 / a ** 2 * (np.log1p(a) - 2.0 * a / (2.0 + a))  # exact load (ours)
    if model == "book":
        return float(quad(lambda s: float(slider_bearing(s, h0, a, L, U, mu, 0.0, "book")), 0.0, float(L),
                          epsabs=0.0, epsrel=1e-12)[0])
    raise ValueError('model must be "exact", "linear" or "book"')


def slider_optimum_taper() -> dict:
    """Taper α that maximises the exact slider load at fixed exit gap h₀, U, μ, L (ours; V5 San Andrés 2012).

    Book: §8.3, Example 8.1 stops at O(α); the optimum follows from :func:`slider_bearing_load` ("exact").
    Returns dict(alpha_opt, K_opt = 1 + α (inlet/exit gap ratio), W_star = W h₀²/(6μUL²) (San Andrés' normalisation),
    W_coefficient = W h₀²/(μUL²)). Method: ``scipy.optimize.minimize_scalar`` (bounded, xatol 1e-12).
    Validation — tests/test_ch08.py: test_slider_optimum_taper_V5_san_andres,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V5 K_opt ≈ 2.1889, W* ≈ 0.0267 (reference/ch08). Label: benchmark.
    """
    res = minimize_scalar(lambda a: -slider_bearing_load(1.0, a, 1.0, 1.0, 1.0, "exact"), bounds=(0.05, 10.0),
                          method="bounded", options={"xatol": 1e-12})
    a = float(res.x)
    Wc = float(slider_bearing_load(1.0, a, 1.0, 1.0, 1.0, "exact"))
    return dict(alpha_opt=a, K_opt=1.0 + a, W_star=Wc / 6.0, W_coefficient=Wc)


def slider_bearing_state(h0, alpha, L, U, mu: float = 1e-3, p_e: float = 0.0) -> dict:
    """Every number the slider-bearing explainer shows (curation §9, C07/E3).

    Book: §8.3, Example 8.1: C₁ = ∫₀ʰ(u − U)dy = −(1 + α)Uh₀/(2 + α) (the flux in the pad frame), exact and linear
    pressure and load. The pressure peak is where dp/dx = 0, h* = −2C₁/U = 2(1 + α)h₀/(2 + α), i.e. x = L/(2 + α) (ours).
    Recirculation next to the pad (ours, not the book's). In the pad frame the floor moves at −U, the pad is at rest
    and the flux is C₁ = −U h_m/2 with h_m = h* the gap at the pressure peak. With η = y/h and (8.19),
    (u − U)/U = −(1 − η) + 3(1 − h_m/h) η(1 − η), so near the pad (η → 1) u − U ≈ −U(3h_m/h − 2)(1 − η), while at
    the floor u − U = −U (no reversal there; the profile is a parabola in η, so one sign change at most). The flow next to the pad opposes the bulk flux iff 3h_m/h − 2 < 0, i.e.
    **h > 1.5 h_m**; the reversed velocity grows with h, so it is strongest where the gap is widest. With
    h = h₀(1 + αx/L): for α > 0 the widest gap is h₀(1 + α) at x = L and the criterion is α > 1; for α < 0 it is h₀ at
    x = 0 and the criterion is 1/(1 + α) > 2, i.e. α < −½. Both are "wide-to-narrow gap ratio > 2". Independent of
    the sign of U (C₁ and the profile flip together). For αU < 0 (suction pad) the wide end is downstream.
    Evaluated on a 399-point y grid at the wide end.

    Parameters
    ----------
    h0 : [m];  alpha : [–] (> −1);  L : [m];  U : [m/s];  mu : [Pa s];  p_e : [Pa].

    Returns
    -------
    dict of floats: C1 [m²/s] (the pad-frame flux ∫₀ʰ(u − U)dy, constant; the ground-frame flux C₁ + Uh(x) varies with
    x), p_max [Pa] (p_e + peak), dp_max = p_max − p_e [Pa], x_pmax [m], h_pmax [m], W_exact, W_linear [N/m],
    err_linear = 100(W_linear/W_exact − 1) [%], p_max_atm = dp_max/P_ATM [–], p_visc = μUL/h₀² [Pa] (the viscous
    pressure scale), inlet_backflow (bool; recirculation at the wide end, wherever it is: α > 1 or α < −½, any
    U ≠ 0), backflow_x [m] (station where recirculation is strongest = the wide end, L or 0; NaN if none),
    backflow_any (bool; recirculation anywhere in the gap — equal to inlet_backflow, since the wide end is the worst).

    Validation — tests/test_ch08.py: test_slider_optimum_taper_V5_san_andres,
      test_slider_bearing_state_V1_explainer_numbers, test_slider_backflow_both_ends (α at each threshold ± 0.01, against the sign
      of :func:`slider_gap_velocity` (``frame="pad"``) on an x–y grid),
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable,
      test_book_V6_section_8_3_lubrication_forms_and_numbers.
    Checks: V1 against :func:`slider_bearing` (max over a fine grid) and :func:`slider_bearing_load`.
    Label: analytic.
    """
    h0, a, L, U, mu = float(h0), float(alpha), float(L), float(U), float(mu)
    C1 = -(1.0 + a) * U * h0 / (2.0 + a)  # Example 8.1
    x_star = L / (2.0 + a)
    dpm = float(slider_bearing(x_star, h0, a, L, U, mu, 0.0, "exact"))
    We = float(slider_bearing_load(h0, a, L, U, mu, "exact"))
    Wl = float(slider_bearing_load(h0, a, L, U, mu, "linear"))
    # Recirculation is strongest where the gap is widest (ours, see docstring): x = L for α > 0, x = 0 for α < 0.
    x_wide = L if a > 0.0 else 0.0
    h_wide = h0 * (1.0 + a * x_wide / L)
    dpdx_w = -12.0 * mu * C1 / h_wide ** 3 - 6.0 * mu * U / h_wide ** 2  # Example 8.1
    yy = np.linspace(0.0, h_wide, 401)[1:-1]
    u_pad = _F(lubrication_velocity(yy, h_wide, dpdx_w, 0.0, -U, mu))  # (8.19) pad frame: floor −U, pad 0
    back = bool(C1 != 0.0 and np.any(np.sign(u_pad) == -np.sign(C1))) if U != 0 else False
    return dict(C1=C1, p_max=float(p_e) + dpm, dp_max=dpm, x_pmax=x_star, h_pmax=h0 * (1.0 + a * x_star / L),
                W_exact=We, W_linear=Wl, err_linear=100.0 * (Wl / We - 1.0) if We != 0 else np.nan,
                p_max_atm=dpm / P_ATM, p_visc=mu * U * L / h0 ** 2, inlet_backflow=back,
                backflow_x=x_wide if back else np.nan, backflow_any=back)


def slider_gap_velocity(x, y, h0, alpha, L, U, mu: float = 1e-3, frame: str = "ground"):
    """Velocity u(x, y) in the slider-bearing gap from (8.19) with the exact Example 8.1 pressure gradient.

    Book: §8.3, Example 8.1 (profile u = −(h²/2μ)(dp/dx)(y/h)(1 − y/h) + Uy/h, dp/dx = −12μC₁/h³ − 6μU/h²).
    ``frame="ground"`` (alias ``"lab"``): floor at rest, pad moving at U (the book's); ``frame="pad"``: u − U (steady
    streamlines, tracers).
    Points outside 0 ≤ y ≤ h(x) return NaN.

    Parameters
    ----------
    x, y : [m] (broadcast);  h0 : [m];  alpha : [–];  L : [m];  U : [m/s];  mu : [Pa s];  frame : "ground" | "pad".

    Returns
    -------
    u : [m/s].

    Validation — tests/test_ch08.py: test_slider_bearing_state_V1_explainer_numbers,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V1 ∫₀ʰ(u − U)dy = C₁ at every x (``quad``). Label: analytic.
    """
    x_, y_ = np.broadcast_arrays(_F(x), _F(y))
    h0, a, L, U, mu = float(h0), float(alpha), float(L), float(U), float(mu)
    hx = h0 * (1.0 + a * x_ / L)
    C1 = -(1.0 + a) * U * h0 / (2.0 + a)
    dpdx = -12.0 * mu * C1 / hx ** 3 - 6.0 * mu * U / hx ** 2  # Example 8.1
    u = _F(lubrication_velocity(y_, hx, dpdx, U, 0.0, mu))
    if frame == "pad":
        u = u - U
    elif frame not in ("ground", "lab"):
        raise ValueError('frame must be "ground" or "pad"')
    return _S(np.where((y_ >= 0.0) & (y_ <= hx) & (x_ >= 0.0) & (x_ <= L), u, np.nan))


# ======================================================================================================================
# Example 8.2 Hele-Shaw flow
# ======================================================================================================================
def hele_shaw_velocity(z, h, grad_p, mu: float = 1e-3):
    """In-plane velocity between parallel plates z = 0, h driven by an in-plane pressure gradient.

    Book: §8.3, Example 8.2: u ≅ −(1/2μ)(∂p/∂x) z(h − z), v ≅ −(1/2μ)(∂p/∂y) z(h − z) (from the lubrication equations
    with ν restored, analysis §9 R7).

    Parameters
    ----------
    z : height in the gap [m];  h : gap [m];  grad_p : (∂p/∂x, ∂p/∂y) [Pa/m] (floats or arrays);  mu : [Pa s].

    Returns
    -------
    (u, v) : [m/s].

    Validation — tests/test_ch08.py: test_hele_shaw_V1_velocity_potential_and_mean,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V1 u = v = 0 at z = 0, h; = ∂φ/∂x, ∂φ/∂y of :func:`hele_shaw_potential`. Label: analytic.
    """
    f = -_F(z) * (float(h) - _F(z)) / (2.0 * float(mu))
    return _S(f * _F(grad_p[0])), _S(f * _F(grad_p[1]))


def hele_shaw_potential(p, z, h, mu: float = 1e-3):
    """Velocity potential of Hele-Shaw flow in each plane z = const, φ = −z(h − z)p/(2μ).

    Book: §8.3, Example 8.2 (∇²p = 0 ⇒ ∇²φ = 0: the streamlines of 2-D ideal flow, (6.10), (6.12)).
    Parameters: p [Pa]; z, h [m]; mu [Pa s]. Returns φ [m²/s].
    Validation — tests/test_ch08.py: test_hele_shaw_V1_velocity_potential_and_mean,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V2 ∇²φ = 0 for harmonic p.
    Label: analytic, symbolic.
    """
    return _S(-_F(z) * (float(h) - _F(z)) / (2.0 * float(mu)) * _F(p))


def hele_shaw_mean_velocity(grad_p, h, mu: float = 1e-3):
    """Gap-averaged Hele-Shaw velocity ū = −(h²/12μ)∇p (Darcy-like law; ours, the z-average of Example 8.2).
    Parameters: grad_p (2,) [Pa/m]; h [m]; mu [Pa s]. Returns (ū, v̄) [m/s]. Book: §8.3, Example 8.2. Label: analytic."""
    c = -float(h) ** 2 / (12.0 * float(mu))
    return _S(c * _F(grad_p[0])), _S(c * _F(grad_p[1]))


def hele_shaw_cylinder(x, y, z, U_mean, a, h, mu: float = 1e-3) -> dict:
    """Hele-Shaw flow past a circular disk of radius a trapped between the plates (Fig. 8.10).

    Book: §8.3, Example 8.2 (the gap-averaged flow is 2-D ideal flow; the ideal cylinder (6.35)–(6.36) supplies it).
    Gap-averaged potential φ̄ = U_mean(r + a²/r)cos θ; pressure p = −(12μ/h²)φ̄ (gauge, zero on the y axis); velocity
    at height z: u = (6z(h − z)/h²) ∂φ̄/∂x (parabolic across the gap, mean U_mean far away). Not valid in the O(h)
    layer at the disk (Exercise 8.34) — the no-slip condition there is not met. Points inside the disk return NaN.

    Parameters
    ----------
    x, y : in-plane coordinates [m];  z : height [m] (x, y, z broadcast together; ``None`` → gap averages);
    U_mean : far-field gap-averaged speed [m/s];  a : disk radius [m];  h : gap [m];  mu : [Pa s].

    Returns
    -------
    dict of arrays: p [Pa], u, v [m/s] at height z (gap averages if z is None), phi = −z(h − z)p/2μ [m²/s] (the
    velocity potential of that plane; φ̄ for z = None), u_mean, v_mean [m/s], psi_mean (gap-averaged stream function
    U(r − a²/r) sin θ) [m²/s].

    Validation — tests/test_ch08.py: test_hele_shaw_cylinder_V1_equals_ideal_cylinder,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V1 gap average equals the ``core.potential`` ideal cylinder; V3 grid Laplace solve
    (:func:`hele_shaw_streamfunction_grid`). Label: analytic, converged.
    """
    avg = z is None
    X, Y, Z = np.broadcast_arrays(_F(x), _F(y), _F(0.0 if avg else z))
    U, a, h, mu = float(U_mean), float(a), float(h), float(mu)
    r2 = X ** 2 + Y ** 2
    inside = r2 < a ** 2 * (1.0 - 1e-12)
    with np.errstate(divide="ignore", invalid="ignore"):
        phi = U * X * (1.0 + a ** 2 / r2)
        um = U * (1.0 - a ** 2 * (X ** 2 - Y ** 2) / r2 ** 2)
        vm = -U * a ** 2 * 2.0 * X * Y / r2 ** 2
        psi = U * Y * (1.0 - a ** 2 / r2)
    shape = np.ones_like(Z) if avg else 6.0 * Z * (h - Z) / h ** 2  # u(z)/ū: parabola with unit mean
    nan = lambda A: np.where(inside, np.nan, A)  # noqa: E731
    return dict(p=_S(nan(-12.0 * mu / h ** 2 * phi)), u=_S(nan(shape * um)), v=_S(nan(shape * vm)),
                phi=_S(nan(shape * phi)),  # φ = −z(h − z)p/2μ at height z (= φ̄ for the gap average)
                u_mean=_S(nan(um)), v_mean=_S(nan(vm)), psi_mean=_S(nan(psi)))


def hele_shaw_streamfunction_grid(n: int = 129, a: float = 1.0, U_mean: float = 1.0, box: float = 4.0,
                                  method: str = "direct") -> dict:
    """Gap-averaged Hele-Shaw stream function around a disk from a finite-difference Laplace solve (V3 route).

    Book: §8.3, Example 8.2 (∇²ψ̄ = 0 like 2-D ideal flow) solved with the §6.7 masked five-point solver
    (``core.laplace_solvers.solve_laplace``): ψ̄ = 0 on the disk (staircase boundary: the nodes inside), the exact
    ideal-flow value U(r − a²/r)sin θ on the square box |x|, |y| ≤ box. Our numerical check, not the book's method.

    Parameters
    ----------
    n : nodes per side (odd keeps the axes on the grid);  a : disk radius [m];  U_mean : [m/s];  box : half-width [m];
    method : solver method of ``solve_laplace``.

    Returns
    -------
    dict: x, y (1-D) [m], psi, psi_exact (2-D, [j, i]) [m²/s], err_far (max |error| where r ≥ 2a) [m²/s], h [m].

    Validation — tests/test_ch08.py: test_hele_shaw_grid_V3_staircase_first_order,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V3 err_far decreases under refinement (staircase: order between 1 and 2). Label: converged.
    """
    from .laplace_solvers import solve_laplace
    xs = np.linspace(-box, box, int(n))
    X, Y = np.meshgrid(xs, xs, indexing="xy")
    r2 = X ** 2 + Y ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        exact = np.where(r2 >= a ** 2, U_mean * Y * (1.0 - a ** 2 / r2), 0.0)
    mask = np.ones_like(X, bool)
    mask[0, :] = mask[-1, :] = mask[:, 0] = mask[:, -1] = False
    solid = r2 < a ** 2
    mask &= ~solid
    bc = np.where(mask, 0.0, exact)
    bc[solid] = 0.0
    dx = xs[1] - xs[0]
    psi, _ = solve_laplace(mask, bc, method=method, dx=dx)
    far = r2 >= (2.0 * a) ** 2
    return dict(x=xs, y=xs, psi=psi, psi_exact=exact, err_far=float(np.max(np.abs(psi - exact)[far])), h=float(dx))


# ======================================================================================================================
# Example 8.3 / 8.7 gravity spreading of a viscous bead
# ======================================================================================================================
def thin_film_velocity(y, h, h_x, rho: float = 1000.0, g: float = G0, mu: float = 1.0):
    """Velocity in a spreading viscous film: u ≅ −(ρg/2μ)(∂h/∂x) y(2h − y) (no slip at y = 0, no stress at y = h).

    Book: §8.3, Example 8.3 (hydrostatic p = ρg(h − y) in (8.18)). Parameters: y, h [m]; h_x = ∂h/∂x [–]; rho [kg/m³];
    g [m/s²]; mu [Pa s]. Returns u [m/s].
    Validation — tests/test_ch08.py: test_thin_film_V1_profile_walls_and_flux,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V1 BCs, du/dy(h) = 0. Label: analytic.
    """
    y_ = _F(y)
    return _S(-float(rho) * float(g) / (2.0 * float(mu)) * _F(h_x) * y_ * (2.0 * _F(h) - y_))  # Example 8.3


def thin_film_flux(h, h_x, rho: float = 1000.0, g: float = G0, mu: float = 1.0):
    """Volume flux per unit width in a gravity-spreading film, q = ∫₀ʰ u dy ≅ −(ρg/3μ) h³ ∂h/∂x.

    Book: §8.3, Example 8.3 (then ∂h/∂t = (ρg/3μ)∂(h³∂h/∂x)/∂x). Parameters: h [m]; h_x [–]; rho; g; mu.
    Returns q [m²/s].
    Validation — tests/test_ch08.py: test_thin_film_V1_profile_walls_and_flux,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable,
      test_book_V6_section_8_3_lubrication_forms_and_numbers.
    Checks: V1 ``quad`` of :func:`thin_film_velocity`; wrong variant h³/(2μ).
    Label: analytic.
    """
    return _S(-float(rho) * float(g) / (3.0 * float(mu)) * _F(h) ** 3 * _F(h_x))  # Example 8.3


def thin_film_spread(h0, x, t_out, rho: float = 1000.0, g: float = G0, mu: float = 1.0, h_min: float = 1e-6,
                     cache=None, dt0: float | None = None, growth: float = 1.05, tol: float = 1e-12, max_iter: int = 30,
                     front_level: float | None = None) -> dict:
    """Numerical solution of the thin-film (viscous gravity current) equation ∂h/∂t = β ∂(h³ ∂h/∂x)/∂x, β = ρg/3μ.

    Book: §8.3, Example 8.3 (the book derives the equation and stops; Example 8.7 gives the similarity form). The scheme
    is **ours**: conservative finite volume on uniform cells (face flux β H³(h_{i+1} − h_i)/Δx with H the arithmetic mean
    of the neighbours), zero-flux ends, backward Euler in time with Newton iterations on the tridiagonal Jacobian
    (``solve_banded``), a precursor film h ≥ h_min (the equation is degenerate where h → 0), and time steps that grow
    by ``growth`` from ``dt0`` and land exactly on each output time. Volume Σh Δx is conserved to the Newton tolerance.

    Parameters
    ----------
    h0 : array_like, shape (N,) — initial thickness at the cell centres [m] (values below h_min are raised to h_min).
    x : array_like, shape (N,) — uniformly spaced cell centres [m].
    t_out : sequence of increasing output times [s] (> 0).
    rho, g, mu : ρ [kg/m³], g [m/s²], μ [Pa s].
    h_min : precursor thickness [m].
    cache : optional path of an ``.npz`` file: if it exists and was written for exactly the same inputs (a SHA-256
        key of every argument is stored with it) it is loaded; otherwise the run is computed and saved there.
        Invalidation rule: any change of an input changes the key and forces a recompute.
    dt0 : first time step [s] (default: 1e-3 × the diffusive time Δx²/(βH³) of the initial maximum H).
    growth : step growth factor per step;  tol : Newton tolerance on the max residual × dt [m];  max_iter : per step.
    front_level : thickness defining the front [m] (default max(10 h_min, 1e-3 max h0)).

    Returns
    -------
    dict: t (M,) [s], h (M, N) [m], volume (M,) [m²] (Σ h Δx, precursor included), x_front (M,) [m] (largest x with
    h ≥ front_level, linear interpolation), n_steps (int), newton_max (int), rejected (int: halved steps).

    Validation — tests/test_ch08.py: test_thin_film_spread_V4_volume_conserved,
      test_thin_film_spread_V5_huppert_shape_and_t_one_fifth,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V4 volume conserved (1e-12 relative); V3/V5 front → Huppert x_N ∝ t^{1/5}, shape
    (:func:`viscous_current_similarity`). No grid-refinement study is run, so no convergence claim.
    Label: conserved, benchmark.
    """
    if cache is not None:
        import hashlib
        from pathlib import Path

        key = hashlib.sha256()
        for arr in (_F(h0), _F(x), _F(t_out)):
            key.update(np.ascontiguousarray(arr).tobytes())
        key.update(repr((float(rho), float(g), float(mu), float(h_min), dt0, float(growth), float(tol), int(max_iter),
                         front_level)).encode())
        digest = key.hexdigest()
        path = Path(cache)
        if path.exists():
            with np.load(path, allow_pickle=False) as z:
                if "key" in z and str(z["key"]) == digest:
                    return dict(t=z["t"], h=z["h"], volume=z["volume"], x_front=z["x_front"],
                                n_steps=int(z["n_steps"]), newton_max=int(z["newton_max"]), rejected=int(z["rejected"]))
        res = thin_film_spread(h0, x, t_out, rho, g, mu, h_min, None, dt0, growth, tol, max_iter, front_level)
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(path, key=np.array(digest), **res)
        return res
    x = _F(x)
    h = np.maximum(_F(h0).copy(), float(h_min))
    N = h.size
    dx = float(x[1] - x[0])
    beta = float(rho) * float(g) / (3.0 * float(mu))
    t_out = [float(t) for t in np.atleast_1d(t_out)]
    fl = float(front_level) if front_level is not None else max(10.0 * float(h_min), 1e-3 * float(np.max(h)))
    dt = float(dt0) if dt0 is not None else 1e-3 * dx ** 2 / (beta * float(np.max(h)) ** 3)

    def flux_and_jac(hh):
        H = 0.5 * (hh[1:] + hh[:-1])
        dh = (hh[1:] - hh[:-1]) / dx
        Fc = beta * H ** 3 * dh  # β h³ ∂h/∂x at the faces (the physical flux is −Fc)
        dF_l = beta * (1.5 * H ** 2 * dh - H ** 3 / dx)  # ∂Fc/∂h_i (left cell)
        dF_r = beta * (1.5 * H ** 2 * dh + H ** 3 / dx)  # ∂Fc/∂h_{i+1}
        return Fc, dF_l, dF_r

    def front(hh):
        idx = np.nonzero(hh >= fl)[0]
        if idx.size == 0:
            return np.nan
        i = idx[-1]
        if i + 1 >= N:
            return float(x[-1])
        return float(x[i] + (hh[i] - fl) / (hh[i] - hh[i + 1]) * dx)

    # DEVIATION: Newton (not the Picard iteration named in analysis I19 / design 2.13) — quadratic convergence, same
    # fixed point; the book gives no scheme at all (Example 8.3 stops at the PDE)
    def newton(hn, step):
        """Backward-Euler step by Newton; returns (h, iterations) or (None, iterations) on failure."""
        hk = hn.copy()
        for it in range(max_iter):
            Fc, dl, dr = flux_and_jac(hk)
            div = np.zeros(N)
            div[:-1] += Fc
            div[1:] -= Fc  # div_i·dx = Fc_{i+1/2} − Fc_{i−1/2}, zero-flux ends
            R = (hk - hn) / step - div / dx
            ab = np.zeros((3, N))
            main = np.full(N, 1.0 / step)
            main[:-1] -= dl / dx  # from +Fc_{i+1/2} (depends on h_i via dl)
            main[1:] += dr / dx  # from −Fc_{i−1/2} (depends on h_i via dr)
            ab[1] = main
            ab[0, 1:] = -dr / dx  # ∂R_i/∂h_{i+1}
            ab[2, :-1] = dl / dx  # ∂R_{i+1}/∂h_i
            try:
                d = solve_banded((1, 1), ab, -R)
            except (ValueError, np.linalg.LinAlgError):
                return None, it + 1
            hk = hk + d
            if not np.all(np.isfinite(hk)) or np.min(hk) <= 0.0:
                return None, it + 1
            if np.max(np.abs(d)) <= tol:
                return hk, it + 1
        return None, max_iter

    t = 0.0
    out_t, out_h, out_v, out_f = [], [], [], []
    nsteps, newton_max, rejected = 0, 0, 0
    for T in t_out:
        while t < T * (1.0 - 1e-14):
            step = min(dt, T - t)
            with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
                hnew, its = newton(h, step)
            if hnew is None:  # reject: halve the step and retry
                dt = 0.5 * step
                rejected += 1
                if dt < 1e-14 * max(T, 1.0):
                    raise RuntimeError("thin_film_spread: Newton failed at a vanishing time step")
                continue
            newton_max = max(newton_max, its)
            h = hnew
            t += step
            nsteps += 1
            if step >= dt * (1.0 - 1e-12) and its <= 6:
                dt = step * float(growth)
            elif its > 10:
                dt = 0.7 * step
        out_t.append(t)
        out_h.append(h.copy())
        out_v.append(float(np.sum(h) * dx))
        out_f.append(front(h))
    return dict(t=np.array(out_t), h=np.array(out_h), volume=np.array(out_v), x_front=np.array(out_f),
                n_steps=nsteps, newton_max=newton_max, rejected=rejected)


def viscous_current_eta_N() -> float:
    """Front constant of the 2-D viscous gravity current, η_N = [(1/5)(3/10)^{1/3} π^{1/2} Γ(1/3)/Γ(5/6)]^{−3/5}
    ≈ 1.4112 (Huppert 1982, as summarised by Ball & Huppert). Book: §8.4, Example 8.7 (leaves A and D open).
    Label: benchmark."""
    return float((0.2 * 0.3 ** (1.0 / 3.0) * np.sqrt(np.pi) * _gamma(1.0 / 3.0) / _gamma(5.0 / 6.0)) ** (-0.6))


def viscous_current_similarity(x, t, area, rho: float = 1000.0, g: float = G0, mu: float = 1.0,
                               return_front: bool = False):
    """Similarity solution of the spreading bead (2-D viscous gravity current on a horizontal plate).

    Book: §8.4, Example 8.7 gives the form h = At^{−1/5}F(x/Dt^{1/5}) and stops; the closed form is Huppert's
    (J. Fluid Mech. 121, 43–58, 1982; as written by Ball & Huppert): x_N = η_N(βA³t)^{1/5},
    h = (3/10)^{1/3} η_N^{2/3} (A²/β)^{1/5} t^{−1/5} (1 − x²/x_N²)^{1/3}, β = ρg/3μ, A = ∫₀^{x_N} h dx (half-area).

    Parameters
    ----------
    x : [m] (symmetric about 0);  t : [s] (> 0);  area : half cross-section area A [m²];  rho, g, mu;
    return_front : also return the front position.

    Returns
    -------
    h : thickness [m] (0 beyond the front);  with ``return_front=True``: (h, x_N [m]).

    Validation — tests/test_ch08.py: test_viscous_current_similarity_V1_pde_volume_front,
      test_thin_film_spread_V5_huppert_shape_and_t_one_fifth,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V1 PDE residual and ∫₀^{x_N} h dx = A; V5 η_N; V3 :func:`thin_film_spread` converges to it.
    Label: analytic, benchmark.
    """
    beta = float(rho) * float(g) / (3.0 * float(mu))
    eN = viscous_current_eta_N()
    A = float(area)
    t_ = _F(t)
    xN = eN * (beta * A ** 3 * t_) ** 0.2
    hc = 0.3 ** (1.0 / 3.0) * eN ** (2.0 / 3.0) * (A ** 2 / beta) ** 0.2 * t_ ** -0.2
    s = 1.0 - _F(x) ** 2 / xN ** 2
    h = _S(np.where(s > 0, hc * np.cbrt(np.maximum(s, 0.0)), 0.0))
    return (h, _S(xN)) if return_front else h


def thin_film_state(t, area, rho: float = 1000.0, g: float = G0, mu: float = 1.0) -> dict:
    """Numbers the viscous-gravity-current explainer shows (curation §9, E4), from the similarity solution.

    Book: §8.3 Example 8.3, §8.4 Example 8.7 (t^{1/5}). Parameters: t [s]; area half-area A [m²]; rho, g, mu.
    Returns dict of floats: x_N [m], h_centre [m], beta = ρg/3μ [1/(m s)], eta_N [–], effective_diffusivity = βh_c³
    [m²/s], front_speed = dx_N/dt = x_N/(5t) [m/s], aspect = h_centre/x_N [–].
    Validation — tests/test_ch08.py: test_viscous_current_similarity_V1_pde_volume_front,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V1 against :func:`viscous_current_similarity`. Label: analytic.
    """
    h, xN = viscous_current_similarity(0.0, t, area, rho, g, mu, return_front=True)
    beta = float(rho) * float(g) / (3.0 * float(mu))
    return dict(x_N=float(xN), h_centre=float(h), beta=beta, eta_N=viscous_current_eta_N(),
                effective_diffusivity=beta * float(h) ** 3, front_speed=float(xN) / (5.0 * float(t)),
                aspect=float(h) / float(xN))
