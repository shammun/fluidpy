"""Similarity-reduction engines (moved from ``ch08_laminar_flow`` in ch09 because two chapters now use them).

Book: §8.4 (Stokes' first problem, Examples 8.4–8.7: ansatz (8.32)) and §9.3–§9.4, §9.10 (Blasius (9.19)–(9.27),
Falkner–Skan (9.34)–(9.36), free jet (9.59)–(9.65), wall jet (9.82)).  The public names are unchanged and are still
re-exported by ``fluidpy.ch08_laminar_flow``: ``similarity_reduce_sympy``, ``similarity_ode_solve``,
``similarity_collapse_error``.
"""
from __future__ import annotations

import numpy as np
import sympy as sp
from scipy.integrate import solve_bvp
from scipy.special import erfc

from .laminar import line_vortex_decay, stokes_first_problem, vortex_sheet_diffusion
from .lubrication import viscous_current_similarity
from ._util import as_scalar_if_0d

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d


def _z(e) -> bool:
    return sp.simplify(e) == 0


# ======================================================================================================================
# §8.4 similarity engine
# ======================================================================================================================
_SIM_KEYS = ("ode", "brackets", "delta", "n", "m", "bcs", "F", "residual")


def _sim_out(d: dict) -> dict:
    """Ensure every design Part C 4.10 key is present (None where it does not apply)."""
    for k in _SIM_KEYS:
        d.setdefault(k, None)
    return d


def similarity_reduce_sympy(case: str = "stokes1", printed: bool = False) -> dict:
    """Reduce a §8.4 PDE to its similarity ODE with the ansatz (8.32a,b) and solve the exponent bookkeeping (sympy).

    Book: §8.4 — "stokes1": u = UF(η), η = y/√(νt) in (8.20) ⇒ (8.26), BCs (8.27)–(8.28), solution (8.30);
    "stokes1_delta": Example 8.4, u = UF(y/δ(t)) ⇒ brackets [δ′/δ], [ν/δ²] proportional ⇒ δ = √(2C₁νt);
    "vortex_sheet": Example 8.5, ω = At^{−n}F(y/√(νt)), the jump constraint ⇒ n = ½, F = De^{−η²/4};
    "line_vortex": Example 8.6, u_θ = (Γ/2πr)F(r/√(νt)) ⇒ (1/η − η/2)F′ = F″, F = 1 − e^{−η²/4};
    "spreading": Example 8.7, h = At^{−n}F(x/Dt^m) ⇒ 3n + 2m = 1 and m = n (volume) ⇒ n = m = 1/5.

    Ch. 9 adds (Book: §9.3–9.4, 9.10): "blasius" ψ = Uδ(x)f(η) (9.19) => brackets (9.25), δ = √(νx/U) (9.26), f''' + ½ff'' = 0 (9.27);
    "falkner_skan" ψ = √(νxU_e)f(η), U_e = axⁿ (9.34) => (9.36); "free_jet" ψ = (Jνx/Cρ)^{1/3}f(η) (9.64) => 3f''' + ff'' + f'² = 0 (unnumbered, sympy);
    "wall_jet" ψ = (νCx^{1/2})^{1/2}f(η) (9.82) => 4f''' + ff'' + 2f'² = 0 (the book prints coefficient 1 — slip R3; ``printed=True`` returns that
    variant, whose ``residual`` is NOT zero).  These return, besides the standard keys: pde_residual_expr, cancelled_terms (Blasius),
    ode_coefficients (dict of the numbers multiplying f''', ff'', f'², 1), prefactor.

    Parameters
    ----------
    case : "stokes1" | "stokes1_delta" | "vortex_sheet" | "line_vortex" | "spreading" | "blasius" | "falkner_skan" | "free_jet" | "wall_jet".
    printed : (ch09 cases) use the book's printed ODE where it differs from the correct one (wall jet).

    Returns
    -------
    dict of sympy objects, always with the design Part C 4.10 keys ode (Eq in F(η)), brackets, delta, n, m, bcs, F
    (closed form), residual (0) — None where a key does not apply — plus per case reduced, solution, pde_residual
    (the dimensional solution in the PDE, 0), exponent_equations …

    Validation — tests/test_ch08.py: test_stokes_first_V2_sympy_residual_and_wrong_variant,
      test_similarity_reduce_V2_all_cases, test_similarity_example_8_4_V2_derivation,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable,
      test_book_V6_section_8_4_to_8_6_numbers_and_slips.
    Checks: V2 reduced ODEs equal (8.26) and the examples' ODEs; exponents; closed forms. Label: symbolic.
    """
    if case in _CH09_CASES:
        return _reduce_ch09(case, printed)
    eta = sp.Symbol("eta", positive=True)
    y, r, xx = sp.symbols("y r x", positive=True)
    t, nu, U, Gam = sp.symbols("t nu U Gamma", positive=True)
    F = sp.Function("F")
    if case == "stokes1":
        u = U * F(y / sp.sqrt(nu * t))
        pde = sp.diff(u, t) - nu * sp.diff(u, y, 2)
        red = sp.simplify((pde * t / U).subs(y, eta * sp.sqrt(nu * t)).doit())
        ode = sp.Eq(-eta / 2 * F(eta).diff(eta), F(eta).diff(eta, 2))  # (8.26)
        sol = 1 - sp.erf(eta / 2)
        ures = sp.simplify(sp.diff(U * (1 - sp.erf(y / (2 * sp.sqrt(nu * t)))), t)
                           - nu * sp.diff(U * (1 - sp.erf(y / (2 * sp.sqrt(nu * t)))), y, 2))
        return _sim_out(dict(case=case, reduced=red, ode=ode, brackets=None, delta=sp.sqrt(nu * t), m=sp.Rational(1, 2), F=sol,
                    matches_8_26=_z(red + (ode.lhs - ode.rhs)) or _z(red - (ode.lhs - ode.rhs)),
                    dsolve=sp.dsolve(ode), solution=sol, A=-1 / sp.sqrt(sp.pi), B=sp.Integer(1),
                    gaussian_integral=sp.integrate(sp.exp(-sp.Symbol("xi") ** 2 / 4), (sp.Symbol("xi"), 0, sp.oo)),
                    bcs=(sol.subs(eta, 0), sp.limit(sol, eta, sp.oo)),
                    residual=sp.simplify(sp.diff(sol, eta, 2) + eta / 2 * sp.diff(sol, eta)), pde_residual=ures,
                    n=sp.Integer(0)))
    if case == "stokes1_delta":
        d = sp.Function("delta")(t)
        C1 = sp.Symbol("C_1", positive=True)
        u = U * F(y / d)
        pde = sp.diff(u, t) - nu * sp.diff(u, y, 2)
        red = sp.expand(sp.simplify((pde / U).subs(y, eta * d).doit()))
        dsol = sp.dsolve(sp.Eq(d.diff(t) / d, C1 * nu / d ** 2), d)
        dsols = dsol if isinstance(dsol, list) else [dsol]
        pos = [s_.rhs for s_ in dsols]
        Cc = sp.Symbol("C1")
        delta = None
        for e in pos:
            e0 = sp.solve(sp.Eq(e.subs(t, 0), 0), Cc)
            if e0:
                cand = sp.simplify(e.subs(Cc, e0[0]))
                if cand.subs({t: 1, nu: 1, C1: 1}).is_positive:
                    delta = cand
        return _sim_out(dict(case=case, reduced=red, brackets=(-d.diff(t) / d, nu / d ** 2), n=sp.Integer(0),
                    ode=sp.Eq(-sp.Symbol("C_1") * eta * F(eta).diff(eta), F(eta).diff(eta, 2)),
                    bracket_t=-d.diff(t) / d, bracket_nu=nu / d ** 2,
                    delta_ode=sp.Eq(d.diff(t) / d, C1 * nu / d ** 2), delta=delta,
                    delta_matches=_z(delta - sp.sqrt(2 * C1 * nu * t)) if delta is not None else False,
                    C1_recovers_8_25=sp.simplify(sp.sqrt(2 * C1 * nu * t).subs(C1, sp.Rational(1, 2)) - sp.sqrt(nu * t)),
                    residual=sp.simplify(sp.sqrt(2 * C1 * nu * t).subs(C1, sp.Rational(1, 2)) - sp.sqrt(nu * t))))
    if case == "vortex_sheet":
        n, A, D = sp.symbols("n A D", real=True)
        w = A * t ** (-n) * F(y / sp.sqrt(nu * t))
        pde = sp.diff(w, t) - nu * sp.diff(w, y, 2)
        red = sp.simplify((pde * t ** (n + 1) / A).subs(y, eta * sp.sqrt(nu * t)).doit())
        ode = sp.Eq(-n * F(eta) - eta / 2 * F(eta).diff(eta), F(eta).diff(eta, 2))
        # jump constraint: −∫ω dy = −A t^{−n} √(νt) ∫F dη = 2U must be t-independent
        tpow = sp.powsimp(t ** (-n) * sp.sqrt(t), force=True)
        n_sol = sp.solve(sp.Eq(sp.log(tpow).expand(force=True).coeff(sp.log(t)), 0), n)[0]
        Fs = D * sp.exp(-eta ** 2 / 4)
        AD = sp.solve(sp.Eq(-A * D * sp.sqrt(nu) * sp.integrate(sp.exp(-eta ** 2 / 4), (eta, -sp.oo, sp.oo)), 2 * U), A)[0] * D
        wz = sp.simplify(AD * t ** (-n_sol) * sp.exp(-y ** 2 / (4 * nu * t)))
        ub = U * sp.erf(y / (2 * sp.sqrt(nu * t)))
        return _sim_out(dict(case=case, reduced=red, ode=ode, n=n_sol, solution=Fs, F=Fs, delta=sp.sqrt(nu * t),
                    m=sp.Rational(1, 2), brackets=(-n / t, -1 / (2 * t), 1 / t), bcs=("F → 0 as η → ∞", "−∫ω dy = 2U"),
                    residual=sp.simplify((ode.lhs - ode.rhs).subs(n, n_sol).subs(F(eta), Fs).doit()),
                    AD=sp.simplify(AD), omega_z=wz, u=ub, omega_matches=_z(wz + sp.diff(ub, y)),
                    pde_residual=sp.simplify(sp.diff(wz, t) - nu * sp.diff(wz, y, 2))))
    if case == "line_vortex":
        u = Gam / (2 * sp.pi * r) * F(r / sp.sqrt(nu * t))
        pde = sp.diff(u, t) - nu * sp.diff(sp.diff(r * u, r) / r, r)
        red = sp.simplify((pde * 2 * sp.pi * r * t / Gam).subs(r, eta * sp.sqrt(nu * t)).doit())
        ode = sp.Eq((1 / eta - eta / 2) * F(eta).diff(eta), F(eta).diff(eta, 2))
        sol = 1 - sp.exp(-eta ** 2 / 4)
        ut = Gam / (2 * sp.pi * r) * (1 - sp.exp(-r ** 2 / (4 * nu * t)))
        us = Gam / (2 * sp.pi * r) * sp.exp(-r ** 2 / (4 * nu * t))
        opr = lambda q: sp.simplify(sp.diff(q, t) - nu * sp.diff(sp.diff(r * q, r) / r, r))  # noqa: E731
        return _sim_out(dict(case=case, reduced=red, ode=ode, bracket=eta ** 2 / 2, brackets=(eta ** 2 / 2,), solution=sol,
                    F=sol, delta=sp.sqrt(nu * t), n=sp.Integer(1), m=sp.Rational(1, 2),
                    residual=sp.simplify((ode.lhs - ode.rhs).subs(F(eta), sol).doit()),
                    bcs=(sol.subs(eta, 0), sp.limit(sol, eta, sp.oo)), pde_residual=opr(ut), spinup_residual=opr(us),
                    u_theta=ut, u_spinup=us))
    if case == "spreading":
        n, m = sp.symbols("n m", real=True)
        eqs = [sp.Eq(-n - 1, -4 * n - 2 * m), sp.Eq(-n + m, 0)]  # ∂h/∂t vs (ρg/3μ)∂(h³h_x)/∂x; volume
        sol = sp.solve(eqs, [n, m], dict=True)[0]
        A, D, beta = sp.symbols("A D beta", positive=True)
        h = A * t ** (-sol[n]) * F(xx / (D * t ** sol[m]))
        pde = sp.diff(h, t) - beta * sp.diff(h ** 3 * sp.diff(h, xx), xx)
        red = sp.simplify((pde * t ** (sol[n] + 1) / A).subs(xx, eta * D * t ** sol[m]).doit())
        return _sim_out(dict(case=case, exponent_equations=eqs, n=sol[n], m=sol[m], reduced=red,
                    ode=sp.Eq(red, 0), delta=D * t ** sol[m], brackets=("t^(-n-1)", "t^(-4n-2m)", "volume t^(m-n)"),
                    t_free=not red.has(t), residual=sp.simplify(sum(sp.Abs((e.lhs - e.rhs).subs(sol)) for e in eqs))))
    raise ValueError('case must be "stokes1", "stokes1_delta", "vortex_sheet", "line_vortex", "spreading", "blasius", "falkner_skan", "free_jet" or "wall_jet"')


_CH09_CASES = ("blasius", "falkner_skan", "free_jet", "wall_jet")


def _reduce_ch09(case: str, printed: bool = False) -> dict:
    """Sympy similarity reduction of u u_x + v u_y = U_e U_e' + nu u_yy for the four Ch. 9 ansaetze."""
    x, y = sp.symbols("x y", positive=True)
    nu, U, a, n, Cc, K = sp.symbols("nu U a n C K", positive=True)
    eta = sp.Symbol("eta", positive=True)
    f = sp.Function("f")
    F0, F1, F2, F3 = sp.symbols("F0 F1 F2 F3")

    def as_symbols(expr):
        expr = expr.doit()
        rep = {sp.Derivative(f(eta), (eta, 3)): F3, sp.Derivative(f(eta), (eta, 2)): F2, sp.Derivative(f(eta), eta): F1}
        for k_, v_ in rep.items():
            expr = expr.subs(k_, v_)
        return expr.subs(f(eta), F0)

    def residual_of(psi, Ue, dUe):
        u = sp.diff(psi, y)
        v = -sp.diff(psi, x)
        return u, v, u * sp.diff(u, x) + v * sp.diff(u, y) - Ue * dUe - nu * sp.diff(u, y, 2)

    if case == "blasius":
        d = sp.Function("delta")(x)
        psi = U * d * f(y / d)
        u, v, R = residual_of(psi, U, 0)
        Rs = sp.expand(as_symbols(R.subs(y, eta * d)))
        c_adv_u = sp.expand(as_symbols((u * sp.diff(u, x)).subs(y, eta * d)))
        c_adv_v = sp.expand(as_symbols((v * sp.diff(u, y)).subs(y, eta * d)))
        canc = (sp.simplify(c_adv_u.coeff(F1 * F2)), sp.simplify(c_adv_v.coeff(F1 * F2)))
        b_ff2 = sp.simplify(-Rs.coeff(F0 * F2))  # U^2 delta'/delta   (Eq. (9.25))
        b_f3 = sp.simplify(-Rs.coeff(F3))  # nu U/delta^2
        dsol = sp.sqrt(nu * x / U)  # Eq. (9.26), C = 2, D = 0
        Rsub = sp.simplify(Rs.subs(d.diff(x), dsol.diff(x)).subs(d, dsol))
        pref = sp.simplify(-U ** 2 / x)
        ode = sp.Eq(F3 + F0 * F2 / 2, 0)
        res = sp.simplify(Rsub - pref * (F3 + F0 * F2 / 2))
        return _sim_out(dict(case=case, brackets=(sp.simplify(b_ff2), sp.simplify(b_f3)), delta=dsol, n=sp.Integer(0), m=sp.Rational(1, 2), ode=ode,
                             cancelled_terms=canc, prefactor=pref, residual=res, pde_residual_expr=Rsub, bcs=("f(0)=0", "f'(0)=0", "f'(inf)=1"),
                             ode_coefficients={"fppp": 1, "f f''": sp.Rational(1, 2)}, bracket_ratio=sp.simplify(b_ff2 / b_f3),
                             delta_check=sp.simplify(sp.diff(dsol, x) * dsol - nu / (2 * U))))
    if case == "falkner_skan":
        Ue = a * x ** n
        delta = sp.sqrt(nu * x / Ue)
        psi = sp.sqrt(nu * x * Ue) * f(y * sp.sqrt(Ue / (nu * x)))
        u, v, R = residual_of(psi, Ue, sp.diff(Ue, x))
        Rs = as_symbols(R.subs(y, eta * delta))
        Rs = sp.simplify(sp.powsimp(sp.expand(Rs), force=True))
        pref = -a ** 2 * x ** (2 * n - 1)
        ode = F3 + (n + 1) / 2 * F0 * F2 - n * F1 ** 2 + n
        res = sp.simplify(sp.powsimp(sp.expand(Rs - pref * ode), force=True))
        return _sim_out(dict(case=case, brackets=None, delta=delta, n=n, m=sp.Rational(1, 2), ode=sp.Eq(ode, 0), prefactor=pref, residual=res,
                             pde_residual_expr=Rs, bcs=("f(0)=0", "f'(0)=0", "f'(inf)=1"),
                             ode_coefficients={"fppp": 1, "f f''": (n + 1) / 2, "f'^2": -n, "1": n}))
    if case == "free_jet":
        delta = (Cc * nu ** 2 * x ** 2 / K) ** sp.Rational(1, 3)
        psi = (K * nu * x / Cc) ** sp.Rational(1, 3) * f(y / delta)
        u, v, R = residual_of(psi, 0, 0)
        Rs = sp.simplify(sp.powsimp(sp.expand(as_symbols(R.subs(y, eta * delta))), force=True))
        Rs = sp.expand(Rs)
        c3 = sp.simplify(Rs.coeff(F3))
        pref = c3 / 3
        ode = 3 * F3 + F0 * F2 + F1 ** 2
        res = sp.simplify(sp.powsimp(sp.expand(Rs - pref * ode), force=True))
        return _sim_out(dict(case=case, brackets=None, delta=delta, n=sp.Rational(1, 3), m=sp.Rational(2, 3), ode=sp.Eq(ode, 0), prefactor=pref, residual=res,
                             pde_residual_expr=Rs, bcs=("f(0)=0", "f'(0)=1", "f'(inf)=0"), F=sp.sqrt(6) * sp.tanh(eta / sp.sqrt(6)),
                             ode_coefficients={"fppp": 3, "f f''": 1, "f'^2": 1}, first_integral="3f'' + f f' = C1 = 0; 3f' + f^2/2 = C2 = 3"))
    if case == "wall_jet":
        delta = sp.sqrt(nu * x ** sp.Rational(3, 2) / Cc)
        psi = sp.sqrt(nu * Cc * sp.sqrt(x)) * f(y / delta)
        u, v, R = residual_of(psi, 0, 0)
        Rs = sp.simplify(sp.powsimp(sp.expand(as_symbols(R.subs(y, eta * delta))), force=True))
        Rs = sp.expand(Rs)
        c3 = sp.simplify(Rs.coeff(F3))
        k3 = 1 if printed else 4
        pref = c3 / 4  # correct: residual = pref times (4 F3 + F0 F2 + 2 F1^2)
        ode = k3 * F3 + F0 * F2 + 2 * F1 ** 2
        res = sp.simplify(sp.powsimp(sp.expand(Rs - pref * ode), force=True))
        return _sim_out(dict(case=case, brackets=None, delta=delta, n=sp.Rational(-1, 2), m=sp.Rational(3, 4), ode=sp.Eq(ode, 0), prefactor=pref, residual=res,
                             pde_residual_expr=Rs, printed=printed, bcs=("f(0)=0", "f'(0)=0", "f'(inf)=0"),
                             ode_coefficients={"fppp": k3, "f f''": 1, "f'^2": 2}))
    raise ValueError("unknown case")


def similarity_ode_solve(case: str = "stokes1", eta_max: float = 12.0, n: int = 400, tol: float = 1e-10,
                         full: bool = False):
    """Solve a §8.4 similarity ODE numerically with ``scipy.integrate.solve_bvp`` (independent of the closed form).

    Book: §8.4 — "stokes1": (8.26) F″ = −(η/2)F′, F(0) = 1, F(η_max) = 0 ((8.27)–(8.28) with ∞ truncated);
    "line_vortex": Example 8.6 F″ = (1/η − η/2)F′, F(η₀) = 0 at η₀ = 1e-4 (the axis; error O(η₀²)), F(η_max) = 1.

    Parameters
    ----------
    case : "stokes1" | "line_vortex";  eta_max : truncation of ∞;  n : initial mesh points;  tol : solve_bvp tolerance.

    full : return the diagnostic dict instead of the (eta, F) pair.

    Returns
    -------
    (eta, F) : arrays (2001 points on [η₀, η_max]);  with ``full=True`` dict(eta, F, F_exact (1 − erf(η/2) or
    1 − e^{−η²/4}), max_err, eta_max, success).

    Validation — tests/test_ch08.py: test_similarity_ode_solve_V3_bvp_matches_closed_forms,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: max_err ≤ 1e-7 against the closed form 1 − erf(η/2) (or 1 − e^{−η²/4}), and the solution is insensitive
    to the truncation η_max ∈ {10, 14}. No mesh-refinement order is measured. Label: analytic (closed-form agreement,
    truncation-insensitive).
    """
    if case == "stokes1":
        e0, bc = 0.0, (1.0, 0.0)
        rhs = lambda e, Y: np.vstack([Y[1], -e / 2.0 * Y[1]])  # noqa: E731
        exact = lambda e: erfc(e / 2.0)  # noqa: E731
    elif case == "line_vortex":
        e0, bc = 1e-4, (0.0, 1.0)
        rhs = lambda e, Y: np.vstack([Y[1], (1.0 / e - e / 2.0) * Y[1]])  # noqa: E731
        exact = lambda e: -np.expm1(-e ** 2 / 4.0)  # noqa: E731
    else:
        raise ValueError('case must be "stokes1" or "line_vortex"')
    e = np.linspace(e0, float(eta_max), int(n))
    guess = np.vstack([bc[0] + (bc[1] - bc[0]) * (1 - np.exp(-e)), (bc[1] - bc[0]) * np.exp(-e)])
    sol = solve_bvp(rhs, lambda a, b: np.array([a[0] - bc[0], b[0] - bc[1]]), e, guess, tol=tol, max_nodes=200000)
    ee = np.linspace(e0, float(eta_max), 2001)
    Fv = sol.sol(ee)[0]
    Fe = exact(ee)
    if not full:
        return ee, Fv
    return dict(eta=ee, F=Fv, F_exact=Fe, max_err=float(np.max(np.abs(Fv - Fe))), eta_max=float(eta_max),
                success=bool(sol.success))


def similarity_collapse_error(case: str = "stokes1", n: float = 0.0, m: float = 0.5, times=(1.0, 2.0, 4.0, 8.0),
                              npts: int = 201) -> float:
    """Spread of rescaled profiles for trial similarity exponents — zero at the right exponents (curation §9, E6, C10).

    Book: §8.4, the ansatz (8.32a) γ = At^{−n}F(ξ/δ), δ ∝ t^m (and (8.32b) γ = Aξ^{−n}F(ξ/δ) for the line vortex).
    Profiles of the exact solutions (ν = 1, U = 1, β = ρg/3μ = 1, unit area) at the given times are rescaled —
    Ch. 9 cases (t <-> x): ``blasius`` u vs y/x^m (n = 0, m = ½), ``free_jet`` u·xⁿ (n = ⅓, m = ⅔), ``wall_jet`` u·xⁿ (n = ½, m = ¾);
    ``stokes1``: u·tⁿ vs y/tᵐ (right: n = 0, m = ½); ``vortex_sheet``: ω·tⁿ vs y/tᵐ (n = m = ½); ``line_vortex``:
    rⁿu_θ vs r/tᵐ (n = 1, m = ½); ``spreading``: h·tⁿ vs x/tᵐ (n = m = 1/5, Huppert's solution) — sampled at the same
    rescaled coordinates ξ ∈ [0, ξ_max] and compared.

    Parameters
    ----------
    case : "stokes1" | "vortex_sheet" | "line_vortex" | "spreading";  n, m : trial exponents;  times : [–];
    npts : samples in ξ.

    Returns
    -------
    spread : max over ξ of (max − min across times) divided by the largest |rescaled profile| [–].

    Validation — tests/test_ch08.py: test_similarity_collapse_V1_right_exponents_only,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V1 ≈ 0 (1e-12) at the right exponents, > 1e-2 away from them. Label: analytic.
    """
    ts = np.asarray(times, dtype=float)
    t0 = ts[0]
    if case == "stokes1":
        f, nat, lo = (lambda X, T: _F(stokes_first_problem(X, T, 1.0, 1.0))), 6.0 * np.sqrt(t0), 0.0
    elif case == "vortex_sheet":
        f, nat, lo = (lambda X, T: _F(vortex_sheet_diffusion(X, T, 1.0, 1.0)[1])), 6.0 * np.sqrt(t0), 0.0
    elif case == "line_vortex":
        f, nat, lo = (lambda X, T: _F(line_vortex_decay(X, T, 2.0 * np.pi, 1.0))), 6.0 * np.sqrt(t0), 0.01
    elif case == "spreading":
        f = lambda X, T: _F(viscous_current_similarity(X, T, 1.0, 3.0, 1.0, 1.0))  # noqa: E731  β = 1
        nat, lo = 1.2 * float(viscous_current_similarity(0.0, t0, 1.0, 3.0, 1.0, 1.0, True)[1]), 0.0
    elif case == "blasius":  # u(y; x) = U f'(y sqrt(U/nu x)), U = nu = 1;  right exponents n = 0, m = 1/2
        from . import boundary_layer as _BL
        f, nat, lo = (lambda X, T: _F(_BL.blasius_fields(T, X, 1.0, 1.0)["u"])), 6.0 * np.sqrt(t0), 0.0
    elif case == "free_jet":  # u(y; x) = u0(x) sech^2(eta/sqrt6), J = rho = nu = 1;  right exponents n = 1/3, m = 2/3
        from . import jets as _J
        f, nat, lo = (lambda X, T: _F(_J.free_jet(T, X, 1.0, 1.0, 1.0)["u"])), 6.0 * float(_J.free_jet_thickness(t0, 1.0, 1.0, 1.0)), 0.0
    elif case == "wall_jet":  # u(y; x) = C x^(-1/2) f'(eta), C = f_inf = nu = 1;  right exponents n = 1/2, m = 3/4
        from . import jets as _J
        f, nat, lo = (lambda X, T: _F(_J.wall_jet(T, X, 1.0, 1.0, 1.0)["u"])), 6.0 * float(np.sqrt(t0 ** 1.5)), 0.0
    else:
        raise ValueError('case must be "stokes1", "vortex_sheet", "line_vortex", "spreading", "blasius", "free_jet" or "wall_jet"')
    xi = np.linspace(lo * nat, nat, int(npts)) / t0 ** float(m)
    P = []
    for T in ts:
        X = xi * T ** float(m)
        prof = f(X, T)
        P.append(X ** float(n) * prof if case == "line_vortex" else T ** float(n) * prof)
    P = np.array(P)
    return float(np.max(np.max(P, axis=0) - np.min(P, axis=0)) / max(np.max(np.abs(P[0])), 1e-300))


