"""Chapter 9 — Boundary layers and related topics: Prandtl scaling, thicknesses, Blasius and Falkner–Skan similarity solutions, the von Kármán
momentum integral, Thwaites' method, separation, bluff-body wakes (cylinder, sphere, sports balls), the free and wall jets, secondary flow.

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e (2012), Ch. 9, §§9.1–9.11, Eqs. (9.1)–(9.85), Examples 9.1–9.2.  Every equation was
transcribed from the rendered page images (chapters/pages/ch09/p388–p435, printed pp. 361–408).

Where the physics lives (every public name is re-exported, so ``ch09.<name>`` reaches every callable)
* ``core.boundary_layer`` (BL) — (9.2)–(9.6), (9.9)–(9.11), (9.16)–(9.17), (9.19)–(9.33), (9.34)–(9.36), (9.43)–(9.50), (9.51)–(9.52), marching.
* ``core.jets`` (JET) — (9.53)–(9.85).
* ``core.bluff_body`` (BB) — §9.7–9.9 (regimes, Kármán street, Strouhal, separated-pressure model, sports balls).
* ``core.similarity_reduce`` — the similarity engine (moved from ch08), now with blasius / falkner_skan / free_jet / wall_jet.
* this module — the sympy derivation engines (``bl_nondim_sympy``, ``momentum_integral_sympy``, ``thwaites_sympy``, ``jet_momentum_sympy``,
  ``wall_jet_invariant_sympy``, ``wall_jet_sympy``), Examples 9.1–9.2, the secondary-flow force, the table of printed slips.

Printed slips (analysis §9) coded corrected, printed forms kept as named options that a test must fail:
R1 (9.7) missing squares (``bl_nondim_sympy(printed_9_7=True)``); R2 (9.30) 4.93 vs 4.910 (``blasius_delta99(printed=True)``);
R3 wall-jet ODE coefficient 1 vs 4 (``similarity_reduce_sympy("wall_jet", printed=True)``, ``wall_jet_ode_solve(printed=True)``);
R4 wall-jet separation of variables (``wall_jet_sympy``); R5 one side of the plate (``sides=``); R6 (9.76) 5.6152 vs 7.3319
(``free_jet_halfwidth(printed=True)``); R7 (9.56) kinematic stress; R10 the Magnus sentence (``magnus_sign``); R11 reverse flow
(``falkner_skan(branch="reversed")``); R14 "Chapter 13" should be 12.  Book-quoted numbers live only in git-ignored ``tests/book_values_ch09.json``.
"""
from __future__ import annotations

import numpy as np
import sympy as sp

from .core import boundary_layer as BL  # noqa: F401
from .core import bluff_body as BB  # noqa: F401
from .core import jets as JET  # noqa: F401
from .core.bluff_body import *  # noqa: F401,F403
from .core.boundary_layer import *  # noqa: F401,F403
from .core.jets import *  # noqa: F401,F403
from .core.creeping import oseen_drag_coefficient, stokes_drag_coefficient  # noqa: F401  (recalled: Fig. 9.22 asymptotes)
from .core.laminar import diffusion_thickness, similarity_variable, stokes_first_problem, temporal_bl_wall_stress  # noqa: F401  (√(νt) bridge)
from .core.similarity import pressure_coefficient, reynolds_number, sphere_drag_coefficient, strouhal_number  # noqa: F401
from .core.similarity_reduce import similarity_collapse_error, similarity_ode_solve, similarity_reduce_sympy  # noqa: F401


def _z(e) -> bool:
    return sp.simplify(e) == 0


# ======================================================================================================================
# §9.1 scaling — sympy
# ======================================================================================================================
def bl_nondim_sympy(printed_9_7: bool = False) -> dict:
    """Substitute the stretched variables (9.6) into the steady 2-D Navier–Stokes equations (9.1) and its y-companion and read off the Re powers (sympy).

    Book: §9.1, Eqs. (9.6)–(9.10): x = Lx*, y = Ly*/√Re, u = Uu*, v = Uv*/√Re, p − p∞ = ρU²p*, ν = UL/Re.  Dividing the x-equation by U²/L gives
    (9.7) with coefficients {u*u*_x*: 1, v*u*_y*: 1, p*_x*: 1, u*_x*x*: 1/Re, u*_y*y*: 1}; dividing the y-equation by U²√Re/L gives (9.8)
    with {1/Re on the advective terms, p*_y*: 1, v*_x*x*: 1/Re², v*_y*y*: 1/Re}.  Re → ∞ leaves (9.9) and (9.10).
    ``printed_9_7=True`` returns the form as PRINTED (slip R1: the second-derivative terms have first-order denominators ∂x*, ∂y*, i.e. lack the squares), whose
    residual against the correct (9.7) is non-zero and which is dimensionally inconsistent.
    Returns dict(x_momentum, y_momentum, continuity [sympy Eq in the starred variables], coefficients {term: power of Re} (ints: x-momentum d²u/dx*² → −1,
    d²u/dy*² → 0, the rest 0; y-momentum inertia → −1, pressure 0, d²v/dx*² → −2, d²v/dy*² → −1; "continuity" → 0), coefficient_exprs {term: coefficient as
    in (9.7)/(9.8)}, limit {…}, residual_correct (0), residual_printed_vs_correct, dimension_check (for the requested form: dict(ok, message, x_terms, printed_ok, correct_ok)).
    Validation: V2 — coefficients exactly {1, 1/Re, 1/Re²}.  Label: symbolic."""
    xs, ys = sp.symbols("xs ys", positive=True)
    L, U, rho, Re = sp.symbols("L U rho Re", positive=True)
    us, vs, ps = (sp.Function(n)(xs, ys) for n in ("us", "vs", "ps"))
    nu = U * L / Re
    dx = lambda q, k=1: sp.diff(q, xs, k) / L ** k  # noqa: E731
    dy = lambda q, k=1: sp.diff(q, ys, k) * sp.sqrt(Re) ** k / L ** k  # noqa: E731
    u, v = U * us, U * vs / sp.sqrt(Re)
    p = rho * U ** 2 * ps
    xm = u * dx(u) + v * dy(u) + dx(p) / rho - nu * (dx(u, 2) + dy(u, 2))  # Eq. (9.1), dimensional residual
    ym = u * dx(v) + v * dy(v) + dy(p) / rho - nu * (dx(v, 2) + dy(v, 2))
    xm_s = sp.expand(sp.simplify(xm * L / U ** 2))
    ym_s = sp.expand(sp.simplify(ym * L / (U ** 2 * sp.sqrt(Re))))
    d = sp.Derivative
    terms_x = {"u*du*/dx*": us * sp.diff(us, xs), "v*du*/dy*": vs * sp.diff(us, ys), "dp*/dx*": sp.diff(ps, xs),
               "d2u*/dx*2": sp.diff(us, xs, 2), "d2u*/dy*2": sp.diff(us, ys, 2)}
    terms_y = {"u*dv*/dx*": us * sp.diff(vs, xs), "v*dv*/dy*": vs * sp.diff(vs, ys), "dp*/dy*": sp.diff(ps, ys),
               "d2v*/dx*2": sp.diff(vs, xs, 2), "d2v*/dy*2": sp.diff(vs, ys, 2)}
    # residual = LHS − RHS, so the pressure and viscous terms carry a minus sign; report the coefficients as they stand in (9.7)/(9.8)
    sgn = {"dp*/dx*": -1, "d2u*/dx*2": -1, "d2u*/dy*2": -1, "dp*/dy*": -1, "d2v*/dx*2": -1, "d2v*/dy*2": -1}
    coef = {k: sgn.get(k, 1) * sp.simplify(xm_s.coeff(t)) for k, t in terms_x.items()}
    coef.update({k: sgn.get(k, 1) * sp.simplify(ym_s.coeff(t)) for k, t in terms_y.items()})
    X = sp.Eq(terms_x["u*du*/dx*"] + terms_x["v*du*/dy*"], -terms_x["dp*/dx*"] + terms_x["d2u*/dx*2"] / Re + terms_x["d2u*/dy*2"])  # Eq. (9.7)
    Y = sp.Eq((terms_y["u*dv*/dx*"] + terms_y["v*dv*/dy*"]) / Re,
              -terms_y["dp*/dy*"] + terms_y["d2v*/dx*2"] / Re ** 2 + terms_y["d2v*/dy*2"] / Re)  # Eq. (9.8)
    correct = sp.simplify(xm_s - (X.lhs - X.rhs))
    powers = {k: int(sp.simplify(sp.expand_log(sp.log(sp.Abs(sp.simplify(v))), force=True) / sp.log(Re))) for k, v in coef.items()}  # |coefficient| = Re^power
    powers["continuity"] = 0  # ∂u*/∂x* + ∂v*/∂y* = 0 has no Re
    # dimension bookkeeping (L, T exponents): u ~ (1, −1), x ~ (1, 0), ν ~ (2, −1); ν ∂²u/∂x² needs a squared length in the denominator
    dm = lambda *t: tuple(sum(c * e[i] for c, e in t) for i in (0, 1))  # noqa: E731
    dx_, du_, dnu_ = (1, 0), (1, -1), (2, -1)
    adv = dm((2, du_), (-1, dx_))
    visc_ok, visc_printed = dm((1, dnu_), (1, du_), (-2, dx_)), dm((1, dnu_), (1, du_), (-1, dx_))
    printed = sp.Eq(X.lhs, -terms_x["dp*/dx*"] + sp.diff(us, xs) / Re + sp.diff(us, ys))  # first-order 'denominators' as printed
    if printed_9_7:
        X = printed
    limit = dict(x_momentum=sp.Eq(terms_x["u*du*/dx*"] + terms_x["v*du*/dy*"], -terms_x["dp*/dx*"] + terms_x["d2u*/dy*2"]),  # Eq. (9.9) scaled
                 y_momentum=sp.Eq(sp.Integer(0), -terms_y["dp*/dy*"]))  # Eq. (9.10)
    chk = dict(printed_ok=visc_printed == adv, correct_ok=visc_ok == adv, x_terms=dict(advective=adv, viscous_correct=visc_ok, viscous_printed=visc_printed))
    chk["ok"] = chk["printed_ok"] if printed_9_7 else chk["correct_ok"]
    chk["message"] = ("ok: every term of (9.7) has dimension L/T²" if chk["ok"] else
                      "second derivative not dimensionless: ν ∂²u/∂x with a first-order denominator has dimension L²/T², not L/T² (slip R1)")
    return dict(x_momentum=X, y_momentum=Y, continuity=sp.Eq(sp.diff(us, xs) + sp.diff(vs, ys), 0),
                coefficients=powers, coefficient_exprs=coef, limit=limit, residual_correct=correct, printed_9_7=printed_9_7,
                residual_printed_vs_correct=sp.simplify(sp.expand((printed.lhs - printed.rhs) - xm_s)), dimension_check=chk)


# ======================================================================================================================
# §9.5 momentum integral — sympy
# ======================================================================================================================
def momentum_integral_sympy() -> dict:
    """The derivation of the von Kármán momentum integral equation (9.43) from (9.37)–(9.42), checked in sympy.

    Book: §9.5, Eqs. (9.37)–(9.43).  Returns dict(steps [list of (label, sympy Eq)], identity_9_38 (0: u·continuity + momentum = ∂(u²)/∂x + ∂(vu)/∂y),
    identity_9_42 (0: the integrands of (9.41) and (9.42) coincide pointwise), deficit_integrals (∫(u² − U_eu)dy = −U_e²θ, ∫(u − U_e)dy = −U_eδ* for
    u = U_eF(y/δ)), result [Eq (9.43)]).  Label: symbolic."""
    x, y = sp.symbols("x y", positive=True)
    rho = sp.Symbol("rho", positive=True)
    u, v, tau = (sp.Function(n)(x, y) for n in ("u", "v", "tau"))
    Ue = sp.Function("U_e")(x)
    dUe = sp.diff(Ue, x)
    cont = sp.diff(u, x) + sp.diff(v, y)
    mom = u * sp.diff(u, x) + v * sp.diff(u, y) - Ue * dUe - sp.diff(tau, y) / rho  # (9.37) residual
    conservative = sp.diff(u ** 2, x) + sp.diff(v * u, y) - Ue * dUe - sp.diff(tau, y) / rho  # (9.38)
    id_938 = sp.simplify(sp.expand(u * cont + mom - conservative))
    integrand_941 = sp.diff(u ** 2, x) - Ue * dUe - Ue * sp.diff(u, x)
    integrand_942 = sp.diff(u ** 2 - Ue * u, x) + dUe * (u - Ue)
    id_942 = sp.simplify(sp.expand(integrand_941 - integrand_942))
    # deficit integrals for a self-similar profile u = U_e F(y/δ(x)); F generic through symbols
    d = sp.Function("delta")(x)
    I_d, I_t = sp.symbols("I_delta I_theta", positive=True)  # ∫(1−F)dη, ∫F(1−F)dη
    theta, dstar = I_t * d, I_d * d
    int_mom = -Ue ** 2 * theta  # ∫(u² − U_e u)dy = −U_e² δ ∫F(1−F)dη
    int_def = -Ue * dstar  # ∫(u − U_e)dy = −U_e δ ∫(1−F)dη
    tau0 = sp.Symbol("tau_0")
    lhs_942 = sp.diff(int_mom, x) + dUe * int_def
    steps = [("(9.37)", sp.Eq(u * sp.diff(u, x) + v * sp.diff(u, y), Ue * dUe + sp.diff(tau, y) / rho)),
             ("(9.38)", sp.Eq(sp.diff(u ** 2, x) + sp.diff(v * u, y), Ue * dUe + sp.diff(tau, y) / rho)),
             ("(9.39)", sp.Eq(sp.Symbol("int_ux_dy"), -sp.Symbol("v_inf"))),
             ("(9.40)", sp.Eq(sp.Symbol("int_d(u2)dx_minus_UeUe'_dy") + Ue * sp.Symbol("v_inf"), -tau0 / rho)),
             ("(9.41)", sp.Eq(sp.Symbol("d/dx int u2 dy") - sp.Symbol("int_UeUe'_dy") - Ue * sp.Symbol("int_ux_dy"), -tau0 / rho)),
             ("(9.42)", sp.Eq(lhs_942, -tau0 / rho))]
    result = sp.Eq(tau0 / rho, sp.diff(Ue ** 2 * theta, x) + Ue * dstar * dUe)  # Eq. (9.43)
    return dict(steps=steps, identity_9_38=id_938, identity_9_42=id_942,
                deficit_integrals=dict(momentum=int_mom, displacement=int_def), lhs_9_42=sp.simplify(lhs_942),
                result=result, residual=sp.simplify(sp.expand(sp.diff(Ue ** 2 * theta, x) + Ue * dstar * dUe + lhs_942)))


# ======================================================================================================================
# §9.6 Thwaites — sympy
# ======================================================================================================================
def thwaites_sympy() -> dict:
    """From (9.47) to Thwaites' closed form (9.50), checked in sympy.

    Book: §9.6, Eqs. (9.47)–(9.50).  Steps: (9.47) ⇒ l(λ) = (2 + H)λ + (U_e/2)d(θ²/ν)/dx; substituting θ² = νλ/U_e′ ⇒ U_e d/dx(λ/U_e′) = 2l − 2(2+H)λ = L(λ) (9.48);
    L ≈ 0.45 − 6λ ⇒ the linear ODE (9.49) for Z = θ²/ν with integrating factor U_e⁶ ⇒ Z U_e⁶ = 0.45∫U_e⁵dx + const (9.50).
    Returns dict(identity_9_47_to_l (0), identity_9_48 (0), ode_9_49 [Eq], integrating_factor_check (0), solution_9_50 [Eq], dsolve [sympy solution of (9.49) for U_e = x^n]).
    Label: symbolic."""
    x, nu, U0 = sp.symbols("x nu U0", positive=True)
    Ue = sp.Function("U_e")(x)
    th = sp.Function("theta")(x)
    H = sp.Function("H")(x)
    dU = sp.diff(Ue, x)
    rhs47 = 2 * th ** 2 * dU / nu + Ue * th * sp.diff(th, x) / nu + th ** 2 / nu * H * dU  # second form of (9.47)
    l_form = (2 + H) * th ** 2 / nu * dU + Ue / 2 * sp.diff(th ** 2 / nu, x)
    id47 = sp.simplify(sp.expand(rhs47 - l_form))
    lam = sp.Function("lambda")(x)
    l = sp.Function("l")(x)
    # θ² = νλ/U_e′ ⇒ (9.48)
    l_of_lam = (2 + H) * lam + Ue / 2 * sp.diff(lam / dU, x)
    L_expr = 2 * l - 2 * (2 + H) * lam
    id48 = sp.simplify(sp.expand(Ue * sp.diff(lam / dU, x) - L_expr.subs(l, l_of_lam)))
    Z = sp.Function("Z")(x)
    ode49 = sp.Eq(sp.diff(Z, x) + 6 * dU / Ue * Z, sp.Rational(45, 100) / Ue)
    # integrating factor U_e^6:  (U_e^6 Z)' = U_e^6 (Z' + 6U_e'/U_e Z) = 0.45 U_e^5
    ifac = sp.simplify(sp.diff(Ue ** 6 * Z, x) - Ue ** 6 * (ode49.lhs) )
    ifac_check = sp.simplify(sp.expand(ifac))
    # closed form for a power law U_e = x^n via dsolve
    n = sp.Symbol("n", positive=True)
    Zn = sp.Function("Zn")(x)
    ode_n = sp.Eq(sp.diff(Zn, x) + 6 * n / x * Zn, sp.Rational(45, 100) * x ** (-n))
    dsol = sp.dsolve(ode_n, Zn)
    I5 = sp.integrate(x ** (5 * n), (x, 0, x), conds="none")
    th0sq = sp.Symbol("theta0_sq_U0_6_over_nu", positive=True)
    s_ = sp.Symbol("s", positive=True)
    sol950 = sp.Eq(Z * Ue ** 6, sp.Rational(45, 100) * sp.Integral(Ue.subs(x, s_) ** 5, (s_, 0, x)) + th0sq)  # Eq. (9.50), Z = theta^2/nu
    return dict(identity_9_47_to_l=id47, identity_9_48=id48, ode_9_49=ode49, integrating_factor_check=ifac_check, solution_9_50=sol950,
                dsolve=dsol, I5_power_law=sp.simplify(I5))


# ======================================================================================================================
# §9.10 jets — sympy
# ======================================================================================================================
def jet_momentum_sympy() -> dict:
    """(9.56)–(9.58): the free jet conserves the momentum flux ∫u²dy; then on the Bickley solution ∫u²dy = J/ρ for every x (sympy).

    Book: §9.10, Eqs. (9.56)–(9.58), (9.61)–(9.64), (9.72).  (9.56): ∫[2u u_x + u v_y + v u_y]dy = ∫ν u_yy dy, i.e. d/dx∫u²dy + [uv] = [ν u_y]; both boundary terms vanish because u and
    u_y → 0 (the book writes τ = ν u_y without ρ — slip R7 — kinematic stress is meant).  Returns dict(pointwise_identity (0), boundary_terms [strings],
    C_free_jet (sympy exact 4√6/3), u0_delta_C (J/ρ, x-independent), residual_x_dependence (0)).  Label: symbolic."""
    x, y, J, rho, nu = sp.symbols("x y J rho nu", positive=True)
    u, v = sp.Function("u")(x, y), sp.Function("v")(x, y)
    cont = sp.diff(u, x) + sp.diff(v, y)
    pointwise = sp.simplify(sp.expand(u * cont + (u * sp.diff(u, x) + v * sp.diff(u, y)) - (sp.diff(u ** 2, x) + sp.diff(u * v, y))))
    eta = sp.Symbol("eta", real=True)
    T = sp.Symbol("T", real=True)
    # ∫sech⁴(η/√6)dη: η = √6 artanh T, dη = √6 dT/(1−T²), sech² = 1 − T²  ⇒  √6 ∫(1 − T²)dT over (−1, 1)
    C = sp.simplify(sp.sqrt(6) * sp.integrate(1 - T ** 2, (T, -1, 1)))
    u0 = (J ** 2 / (C ** 2 * rho ** 2 * nu * x)) ** sp.Rational(1, 3)
    delta = (C * rho * nu ** 2 * x ** 2 / J) ** sp.Rational(1, 3)
    flux = sp.simplify(u0 ** 2 * delta * C)  # ∫u²dy = u0² δ ∫f′²dη  (9.61)
    return dict(pointwise_identity=pointwise, boundary_terms=["[u v] -> 0", "[nu u_y] -> 0"], C_free_jet=C, C_value=sp.simplify(sp.Rational(4, 3) * sp.sqrt(6) - C),
                u0_delta_C=flux, residual_x_dependence=sp.simplify(sp.diff(flux, x)), equals_J_over_rho=sp.simplify(flux - J / rho))


def wall_jet_invariant_sympy() -> dict:
    """(9.79)–(9.82): the conserved 'flux of exterior momentum flux' of the wall jet and the resulting exponents (sympy).

    Book: §9.10, Eqs. (9.80)–(9.82).  With u = u₀f′(η), η = y/δ, δ = (νx/u₀)^{1/2}: ∫₀^∞u(∫_y^∞u²dy′)dy = u₀³δ²K (K = ∫f′∫f′²dη); constancy ⇒ d/dx(νx u₀²) = 0 ⇒ u₀ = Cx^{−1/2}
    (dsolve) and δ = (νx^{3/2}/C)^{1/2} ∝ x^{3/4}; then u₀³δ² = C²ν (Eq. (9.85)).  Returns dict(u0_solution, delta, exponent_delta, u0_cubed_delta_squared, residual (0)).
    Label: symbolic."""
    x, nu, C = sp.symbols("x nu C", positive=True)
    u0 = sp.Function("u0")(x)
    K0 = sp.Symbol("K0", positive=True)
    inv = u0 ** 3 * sp.sqrt(nu * x / u0) ** 2  # = ν x u0² K   (9.81)
    ode = sp.Eq(sp.diff(sp.simplify(inv), x), 0)
    # d/dx(ν x u0²) = 0  ⇒  ν x u0² = K0  ⇒  u0 = √(K0/ν) x^{-1/2}
    u0v = sp.Symbol("u0v", positive=True)
    sol = sp.solve(sp.Eq(nu * x * u0v ** 2, K0), u0v)
    u0s = C * x ** sp.Rational(-1, 2)
    delta_s = sp.sqrt(nu * x / u0s)
    return dict(ode=ode, dsolve=sol, u0_solution=u0s, delta=sp.simplify(delta_s), exponent_delta=sp.simplify(sp.log(delta_s ** 2 / nu * C).expand(force=True)),
                u0_cubed_delta_squared=sp.simplify(u0s ** 3 * delta_s ** 2), residual=sp.simplify(sp.diff(u0s ** 3 * delta_s ** 2, x)),
                delta_matches=sp.simplify(delta_s ** 2 - nu * x ** sp.Rational(3, 2) / C))


def wall_jet_sympy() -> dict:
    """First integrals and the implicit solution (9.83) of the wall-jet ODE 4f‴ + ff″ + 2f′² = 0, checked in sympy.

    Book: §9.10 (text between Eq. (9.82) and (9.84)), Eq. (9.83).  Returns dict with
    first_integral_derivative (0: d/dη[4ff″ − 2f′² + f²f′] = f(4f‴ + ff″ + 2f′²)), second_integral_derivative (0: d/dη[f^{−1/2}f′ + f^{3/2}/6] equals the first integral over 4f^{3/2}),
    integrand_correct (1/(f_∞^{3/2}f^{1/2} − f²)) and integrand_printed (1/(f_∞^{3/2}f − f²), slip R4 — its difference is non-zero), g_integral_derivative (0: d/dg of the LHS of (9.83) = 3/(1 − g³)),
    fpp0_over_finf3 (1/72), K1 (1/40).  Label: symbolic."""
    eta, g, finf = sp.symbols("eta g f_inf", positive=True)
    f = sp.Function("f")(eta)
    fp, fpp, fppp = f.diff(eta), f.diff(eta, 2), f.diff(eta, 3)
    ode = 4 * fppp + f * fpp + 2 * fp ** 2
    first = 4 * f * fpp - 2 * fp ** 2 + f ** 2 * fp
    d1 = sp.simplify(sp.expand(sp.diff(first, eta) - f * ode))
    second = f ** sp.Rational(-1, 2) * fp + f ** sp.Rational(3, 2) / 6
    d2 = sp.simplify(sp.diff(second, eta) - first / (4 * f ** sp.Rational(3, 2)))
    F = sp.Symbol("F", positive=True)
    # f' = F^{1/2}(finf^{3/2} − F^{3/2})/6  ⇒  dF/(finf^{3/2}F^{1/2} − F²) = dη/6
    fprime_expr = sp.sqrt(F) * (finf ** sp.Rational(3, 2) - F ** sp.Rational(3, 2)) / 6
    integrand_correct = 1 / (finf ** sp.Rational(3, 2) * sp.sqrt(F) - F ** 2)
    integrand_printed = 1 / (finf ** sp.Rational(3, 2) * F - F ** 2)
    ok_correct = sp.simplify(integrand_correct - 1 / (6 * fprime_expr))
    ok_printed = sp.simplify(integrand_printed - 1 / (6 * fprime_expr))
    lhs = -sp.log(1 - g) + sp.sqrt(3) * sp.atan((2 * g + 1) / sp.sqrt(3)) + sp.log(1 + g + g ** 2) / 2
    dg = sp.simplify(sp.diff(lhs, g) - 3 / (1 - g ** 3))
    apart = sp.apart(1 / (1 - g ** 3), g)
    return dict(first_integral_derivative=d1, second_integral_derivative=d2, integrand_correct_residual=ok_correct, integrand_printed_residual=ok_printed,
                g_integral_derivative=dg, partial_fractions=apart, fpp0_over_finf3=sp.Rational(1, 72), K1=sp.Rational(1, 40),
                ode_correct=sp.Eq(ode, 0), ode_printed=sp.Eq(fppp + f * fpp + 2 * fp ** 2, 0))


# ======================================================================================================================
# §9.8 Kármán street — sympy
# ======================================================================================================================
def karman_street_sympy() -> dict:
    """The alternate-vortex (k = π/a) stability problem of the staggered street, solved exactly in sympy (Kármán's b/a = 0.2805).

    Book: §9.8 (the book quotes the result only).  Perturbation ζ_{r,n} = P_r(−1)ⁿ, ζ̄_{r,n} = Q_r(−1)ⁿ of the double row (upper row −Γ at na + ib/2, lower row +Γ at
    (n + ½)a − ib/2); the lattice sums are exact: Σ_{m≠0}m^{−2}a^{−2} = π²/(3a²), Σ_{m≠0}(−1)^m(ma)^{−2} = −π²/(6a²), Σ_m(dz − ma)^{−2} = (π/a)²csc²(πdz/a),
    Σ_m(−1)^m(dz − ma)^{−2} = (π/a)²cos(πdz/a)csc²(πdz/a), evaluated at dz = ∓a/2 ± ib (πdz/a = ∓π/2 ± iy, y = πb/a).  The 4×4 matrix is built as in
    :func:`fluidpy.core.bluff_body.karman_street_spectrum`; its characteristic polynomial is μ⁴ + 2(σ² − γ²)μ² + (γ² + σ²)² =
    ((μ − γ)² + σ²)((μ + γ)² + σ²) with γ = (πΓ/2a²)|½ − sech²y|.
    Returns dict(matrix [sympy 4×4 in Γ, a, y = πb/a; state (P₁, Q₁, P₂, Q₂)], char_poly [in μ], factors [the two quadratic factors], gamma [growth rate],
    sigma [oscillation frequency], growth_check (0: expanded product minus char_poly), b_over_a_marginal [= arccosh(√2)/π], marginal_check (0)).
    Validation: V1/V2 vs :func:`fluidpy.core.bluff_body.karman_street_growth_closed` (1e-12) and vs the numerical spectrum.  Label: symbolic."""
    Gam, a, y = sp.symbols("Gamma a y", positive=True)
    mu = sp.Symbol("mu")
    # lattice sums at k = 0 and k = π/a for dz12 = −a/2 + i b, dz21 = a/2 − i b, with π dz/a = ∓π/2 ± i y  (y = π b/a)
    def T_at(phase_z, kind):
        base = (sp.pi / a) ** 2
        if kind == 0:
            return base / sp.sin(phase_z) ** 2
        return base * sp.cos(phase_z) / sp.sin(phase_z) ** 2

    z12 = -sp.pi / 2 + sp.I * y
    z21 = sp.pi / 2 - sp.I * y
    simp = lambda e: sp.simplify(sp.expand_complex(e))  # noqa: E731
    T0 = {12: simp(T_at(z12, 0)), 21: simp(T_at(z21, 0))}
    Tk = {12: simp(T_at(z12, 1)), 21: simp(T_at(z21, 1))}
    conj = lambda e: sp.simplify(sp.conjugate(e))  # noqa: E731
    E0, Ek = sp.pi ** 2 / (3 * a ** 2), -sp.pi ** 2 / (6 * a ** 2)
    G = {1: -Gam, 2: Gam}
    c = 1 / (2 * sp.pi * sp.I)
    M = sp.zeros(4, 4)
    for r, s in ((1, 2), (2, 1)):
        key = 12 if r == 1 else 21
        same = G[r] * c * (E0 - Ek)
        M[2 * (r - 1), 2 * (r - 1) + 1] = same + G[s] * c * conj(T0[key])
        M[2 * (r - 1), 2 * (s - 1) + 1] = -G[s] * c * conj(Tk[key])
        M[2 * (r - 1) + 1, 2 * (r - 1)] = -(same + G[s] * c * T0[key])
        M[2 * (r - 1) + 1, 2 * (s - 1)] = G[s] * c * Tk[key]
    M = M.applyfunc(sp.simplify)
    cp = sp.expand(sp.simplify(M.charpoly(mu).as_expr()))
    coeffs = sp.Poly(cp, mu).all_coeffs()  # μ⁴ + p μ² + q  (odd powers vanish)
    p_, q_ = sp.simplify(coeffs[2]), sp.simplify(coeffs[4])
    sech2 = 1 / sp.cosh(y) ** 2
    gamma = sp.pi * Gam / (2 * a ** 2) * (sp.Rational(1, 2) - sech2)  # signed; the growth rate is |γ|
    g2 = sp.simplify(gamma ** 2)
    # γ² + σ² = √q, σ² − γ² = p/2  ⇒  σ² = (√q + p/2)/2
    sq = sp.simplify(sp.sqrt(sp.factor(q_)))
    sigma2 = sp.simplify((sq + p_ / 2) / 2)
    factors = ((mu - sp.Abs(gamma)) ** 2 + sigma2, (mu + sp.Abs(gamma)) ** 2 + sigma2)
    prod = sp.expand((mu - gamma) ** 2 + sigma2) * sp.expand((mu + gamma) ** 2 + sigma2)
    check = sp.simplify(sp.expand(prod) - cp)
    marginal = sp.acosh(sp.sqrt(2)) / sp.pi
    return dict(matrix=M, char_poly=cp, factors=factors, gamma=sp.Abs(gamma), gamma_signed=gamma, sigma=sp.sqrt(sigma2), sigma_squared=sigma2,
                growth_check=check, b_over_a_marginal=marginal,
                marginal_check=sp.simplify(gamma.subs(y, sp.pi * marginal)))


# ======================================================================================================================
# §9.6 Thwaites on the cylinder — sympy
# ======================================================================================================================
def cylinder_thwaites_sympy() -> dict:
    """Thwaites' λ on the ideal cylinder flow U_e = 2U sin φ, x = aφ, derived in sympy: the integral ∫sin⁵, λ(φ), and its limits.

    Book: §9.6, Eqs. (9.44), (9.50).  Steps: (1) ∫₀^φ U_e⁵dx = 32U⁵a ∫₀^φ sin⁵ψ dψ; with c = cos ψ, sin⁵ψ dψ = −(1 − c²)²dc, so F(φ) = ∫_{cos φ}^{1}(1 − c²)²dc =
    8/15 − cos φ + (2/3)cos³φ − (1/5)cos⁵φ; (2) (9.50): θ² = 0.45νF·32U⁵a/(2U sin φ)⁶ = 0.45(νa/2U)F/sin⁶φ; (3) λ = θ²U_e′/ν with U_e′ = (2U/a)cos φ gives
    λ(φ) = 0.45 F(φ) cos φ/sin⁶φ  (U, a, ν cancel); (4) F ≈ φ⁶/6 as φ → 0, so λ(0) = 0.45/6 = 3/40 (the stagnation-point value); λ(π/2) = 0.
    Returns dict(F, F_derivative_residual (0: dF/dφ − sin⁵φ), lam, lam0 (= 3/40), lam_at_90deg (0), lam_82deg (numeric 0.02630), steps [list of (label, expression)]).
    Validation: V2 residual and limit; V1 vs :func:`fluidpy.core.boundary_layer.thwaites_cylinder_closed_form`.  Label: symbolic."""
    phi, psi, c = sp.symbols("phi psi c", positive=True)
    F = sp.integrate((1 - c ** 2) ** 2, (c, sp.cos(phi), 1))  # substitution c = cos ψ
    F = sp.simplify(F)
    F_direct = sp.integrate(sp.sin(psi) ** 5, (psi, 0, phi))
    lam = sp.Rational(45, 100) * F * sp.cos(phi) / sp.sin(phi) ** 6
    lam0 = sp.limit(lam, phi, 0)
    steps = [("substitute c = cos ψ", sp.Eq(sp.sin(psi) ** 5 * sp.Symbol("dpsi"), -(1 - c ** 2) ** 2 * sp.Symbol("dc"))),
             ("integrate", sp.Eq(sp.Symbol("F"), F)), ("insert in (9.50) and (9.44)", sp.Eq(sp.Symbol("lambda"), lam))]
    return dict(F=F, F_direct_difference=sp.simplify(F - F_direct), F_derivative_residual=sp.simplify(sp.diff(F, phi) - sp.sin(phi) ** 5), lam=lam, lam0=lam0,
                lam_at_90deg=sp.simplify(lam.subs(phi, sp.pi / 2)), lam_82deg=float(lam.subs(phi, 82 * sp.pi / 180)), steps=steps)


# ======================================================================================================================
# Examples 9.1, 9.2
# ======================================================================================================================
def example_9_1(closure: str = "falkner_skan") -> dict:
    """Example 9.1 — Thwaites' method on the flat plate (U_e = U), compared with the exact Blasius numbers.

    Book: §9.6, Eqs. (9.50), (9.45)–(9.46), Example 9.1.  U_e const ⇒ θ² = 0.45νx/U.  Returns (all non-dimensional, scaled by √(νx/U) or √Re_x; our numbers):
    theta_coef = θ/√(νx/U) = √0.45 = 0.6708, theta_err = θ_Thwaites/θ_Blasius − 1 (+0.010), delta_star_coef = H(0)θ/√(νx/U) (1.738), delta_star_err,
    cf_sqrtRex = C_f√Re_x = 2l(0)/θ_coef (0.6575), cf_err = cf/Blasius − 1 (−0.010); also l0, H0 and the ratios theta_ratio_to_blasius, delta_star_ratio,
    cf_ratio_to_blasius.  Print at most 4 significant figures in text (rule R16).
    Validation: V1 exact θ; V5 vs Blasius.  Label: analytic, benchmark."""
    c = BL.blasius_constants()
    th = float(np.sqrt(0.45))
    l0, H0 = float(BL.thwaites_l(0.0, closure)), float(BL.thwaites_H(0.0, closure))
    cf = 2.0 * l0 / th
    return dict(theta_coef=th, theta_err=th / c["theta"] - 1.0, delta_star_coef=H0 * th, delta_star_err=H0 * th / c["delta_star"] - 1.0,
                cf_sqrtRex=cf, cf_err=cf / c["cf_coeff"] - 1.0,
                theta_over_delta=th, theta_ratio_to_blasius=th / c["theta"], delta_star_over_delta=H0 * th, delta_star_ratio=H0 * th / c["delta_star"],
                l0=l0, H0=H0, cf_ratio_to_blasius=cf / c["cf_coeff"])


def example_9_2(theta0: float = 0.0, nu: float = 1e-5, U1: float = 1.0, L: float = 1.0, closure: str = "falkner_skan") -> dict:
    """Example 9.2 — Thwaites' method in a diffuser A(x) = A₁(1 + x/L), U_e = U₁/(1 + x/L).

    Book: §9.6, Eq. (9.50), Example 9.2: λ(x/L) = −(0.45/4)[(1 + x/L)⁴ − 1] − (θ₀²U₁/(νL))(1 + x/L)⁴; separation is predicted where λ reaches the criterion:
    (1 + x/L)⁴ = 1.8 for θ₀ = 0 and λ_sep = −0.09 (x/L = 1.8^{1/4} − 1 = 0.15829), and 1.6053 for the exact-Falkner–Skan value λ_sep = −0.0681 (x/L = 0.12563).
    Parameters: theta0 [m] initial momentum thickness at x = 0; nu [m²/s]; U1 [m/s]; L [m]; closure for the numerical run.
    Returns dict(lam [callable: ξ = x/L → λ, closed form], x_sep_over_L [λ_sep = −0.09], x_sep_over_L_fs [λ_sep = −0.0681] (closed forms; None if λ is already below
    the criterion at ξ = 0), theta0 [m] used, xi, lam_numeric [array, numerical Thwaites on ξ ∈ [0, 0.5]], lam_closed_form [array], max_err, x_sep_numeric [ξ, first crossing
    of the closure's criterion in the numerical run], table (ξ, λ at 0.05…0.2)).
    Validation: V1 λ to 1e-12; x_sep = 1.8^{1/4} − 1.  Label: analytic."""
    kappa = theta0 ** 2 * U1 / (nu * L)
    lam_fn = lambda xi: -(0.45 / 4.0) * ((1.0 + np.asarray(xi, float)) ** 4 - 1.0) - kappa * (1.0 + np.asarray(xi, float)) ** 4  # noqa: E731

    def xsep(lam_sep):
        q = (0.45 / 4.0 - lam_sep) / (0.45 / 4.0 + kappa)  # (1 + ξ)⁴ at which λ = lam_sep
        return None if q <= 1.0 else float(q ** 0.25 - 1.0)

    of = BL.outer_flow("diffuser", U1=U1, L=L)
    xi = np.linspace(0.0, 0.5, 2001)
    r = BL.thwaites(xi * L, of, nu, theta0=theta0, closure=closure, stop_at_separation=False)
    lam_cf = lam_fn(xi)
    tab = np.array([0.05, 0.10, 0.15, 0.20])
    return dict(lam=lam_fn, x_sep_over_L=xsep(LAMBDA_SEP_BOOK), x_sep_over_L_fs=xsep(LAMBDA_SEP_FS), theta0=float(theta0),
                xi=xi, lam_numeric=r["lam"], lam_closed_form=lam_cf, max_err=float(np.max(np.abs(r["lam"] - lam_cf))),
                x_sep_numeric=None if r["x_sep"] is None else r["x_sep"] / L, table=(tab, lam_fn(tab)))


# ======================================================================================================================
# §9.11 secondary flow
# ======================================================================================================================
def secondary_flow_radial_force(u_inviscid, u_layer, R, rho: float = 1000.0):
    """Net inward radial force per unit volume [N/m³] on fluid in the bottom layer of a stirred teacup.

    Book: §9.11: the two largest R-momentum terms give ∂p/∂R = ρu_φ²/R; the thin bottom layer has ∂p/∂z ≈ 0, so it feels the INVISCID pressure gradient ρu_e²/R
    while circular motion at the slower layer speed u would need only ρu²/R.  Net inward force = ρ(u_e² − u²)/R > 0 when the layer is slowed by friction (u < u_e).
    Parameters: u_inviscid u_e [m/s]; u_layer [m/s]; R [m]; rho [kg/m³].  Label: analytic (the sign argument; the Ekman-layer flow itself is Ch. 13)."""
    return float(rho * (np.asarray(u_inviscid, float) ** 2 - np.asarray(u_layer, float) ** 2) / np.asarray(R, float)) if np.ndim(u_layer) == 0 else \
        rho * (np.asarray(u_inviscid, float) ** 2 - np.asarray(u_layer, float) ** 2) / np.asarray(R, float)


def secondary_flow_layer_profile(z, delta: float, u_e: float, shape: str = "exponential", n: float = 1.0 / 7.0):
    """ILLUSTRATIVE swirl profile u_φ(z) in the bottom layer of a teacup: 0 at the floor, u_e far above.  Illustrative, not a solution of the equations.

    Book: §9.11 (the layer is slowed by the floor; its speed profile is not derived).  ``shape``: "exponential" u_e(1 − e^{−z/δ}); "linear" u_e·min(z/δ, 1);
    "power" u_e·min(z/δ, 1)ⁿ (n = 1/7 default); "sine" u_e sin(πz/2δ) for z ≤ δ.  z, delta [m]; u_e [m/s]; n [–].  Label: qualitative."""
    zz = np.asarray(z, float)
    e = zz / delta
    if shape == "exponential":
        return u_e * (1.0 - np.exp(-e))
    if shape == "linear":
        return u_e * np.clip(e, 0.0, 1.0)
    if shape == "power":
        return u_e * np.clip(e, 0.0, 1.0) ** n
    if shape == "sine":
        return u_e * np.sin(np.pi * np.clip(e, 0.0, 1.0) / 2.0)
    raise ValueError('shape must be "exponential", "linear", "power" or "sine"')


# ======================================================================================================================
# printed slips
# ======================================================================================================================
def book_slips() -> list[dict]:
    """The printed slips of Ch. 9 (analysis §9) as data: id, where, printed, correct, evaluator (a callable returning the printed and the correct value).

    Each entry: dict(id, where, printed, correct, printed_value, correct_value) — values are numbers (or None for a structural slip) computed by our functions.
    Book: §9.1 (9.7); §9.3 (9.30), (9.33); §9.4; §9.9; §9.10 (9.56), (9.76), wall-jet ODE and integral, (9.85); §9.10 'Chapter 13'.  Label: analytic."""
    c = BL.blasius_constants()
    return [
        dict(id="R1", where="Eq. (9.7)", printed="second derivatives with denominators ∂x*, ∂y* (no squares)", correct="∂²u*/∂x*², ∂²u*/∂y*²",
             printed_value=None, correct_value=None),
        dict(id="R2", where="Eq. (9.30)", printed="η₉₉ = 4.93 (read from the figure)", correct="η₉₉ = 4.910 (root of f′ = 0.99)",
             printed_value=4.93, correct_value=c["eta99"]),
        dict(id="R3", where="wall-jet ODE below Eq. (9.82)", printed="f‴ + ff″ + 2f′² = 0", correct="4f‴ + ff″ + 2f′² = 0", printed_value=1, correct_value=4),
        dict(id="R4", where="wall-jet separation of variables", printed="∫df/(f_∞^{3/2}f − f²)", correct="∫df/(f_∞^{3/2}f^{1/2} − f²)", printed_value=None, correct_value=None),
        dict(id="R5", where="Eq. (9.33)", printed="C_D = 1.33/√Re_L (looks like the plate)", correct="one side of the plate only (two faces: double)", printed_value=None, correct_value=None),
        dict(id="R6", where="Eq. (9.76)", printed="h₉₉ = 5.6152[…]^{1/3}", correct="h₉₉ = 7.3319[…]^{1/3} (sech² = 0.01)", printed_value=JET.H99_PRINTED,
             correct_value=float(JET.free_jet_at_level(0.01)["coeff"])),
        dict(id="R7", where="Eq. (9.56)", printed="∂τ/∂y without 1/ρ", correct="kinematic stress ν∂u/∂y", printed_value=None, correct_value=None),
        dict(id="R10", where="§9.9 Magnus sentence", printed="'Re < Re_cr' twice", correct="second is Re > Re_cr", printed_value=None, correct_value=None),
        dict(id="R11", where="§9.4", printed="reverse-flow solutions for n < −0.0904", correct="none below; a second branch for −0.0904 < n < 0", printed_value=None, correct_value=None),
        dict(id="R14", where="§9.10", printed="turbulent jets: see Chapter 13", correct="Chapter 12", printed_value=13, correct_value=12),
    ]


def derive_all() -> dict:
    """Run every sympy engine of the chapter once and return their residual checks (all zero) — a one-call self-consistency test.  Label: symbolic."""
    out = {}
    for case in ("blasius", "falkner_skan", "free_jet", "wall_jet"):
        out[f"similarity_{case}"] = similarity_reduce_sympy(case)["residual"]
    ni = bl_nondim_sympy()
    out["bl_nondim"] = ni["residual_correct"]
    mi = momentum_integral_sympy()
    out["momentum_integral"] = (mi["identity_9_38"], mi["identity_9_42"], mi["residual"])
    th = thwaites_sympy()
    out["thwaites"] = (th["identity_9_47_to_l"], th["identity_9_48"], th["integrating_factor_check"])
    jm = jet_momentum_sympy()
    out["jet_momentum"] = (jm["pointwise_identity"], jm["residual_x_dependence"], jm["equals_J_over_rho"])
    wi = wall_jet_invariant_sympy()
    out["wall_jet_invariant"] = (wi["residual"], wi["delta_matches"])
    ws = wall_jet_sympy()
    out["wall_jet"] = (ws["first_integral_derivative"], ws["second_integral_derivative"], ws["g_integral_derivative"], ws["integrand_correct_residual"])
    ct = cylinder_thwaites_sympy()
    out["cylinder_thwaites"] = (ct["F_direct_difference"], ct["F_derivative_residual"], ct["lam0"] - sp.Rational(3, 40), ct["lam_at_90deg"])
    ks = karman_street_sympy()
    out["karman_street"] = (ks["growth_check"], ks["marginal_check"])
    return out
