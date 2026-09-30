"""Two-dimensional laminar jets (Ch. 9, §9.10): the free (Bickley–Schlichting) jet and the wall jet (Glauert).

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e, §9.10, Eqs. (9.53)–(9.85).  Both are similarity solutions of the boundary-layer
equation (9.18) with dp/dx = 0; the free jet conserves the momentum flux J (9.57)–(9.58), the wall jet conserves the "flux of
exterior momentum flux" (9.80) because the wall removes ordinary momentum flux.

Coordinates: x along the jet, y across; free jet symmetric about y = 0 (η ∈ (−∞, ∞)); wall jet y = 0 is the wall.  u = ∂ψ/∂y,
v = −∂ψ/∂x.  Symbols: J [N/m] = momentum flux per unit span; u₀(x) centreline (free) or reference speed (wall) [m/s]; Re_x = xu₀/ν.
In code the constants named C in the book are ``C_fj`` (= ∫f′²dη = 4√6/3, free jet) and ``C`` (dimensional, u₀ = Cx^{−1/2}, wall jet).

Printed slips coded corrected (printed forms are named options a test must fail): R3 the wall-jet ODE is 4f‴ + ff″ + 2f′² = 0
(``wall_jet_ode_solve(printed=True)`` integrates the printed coefficient 1); R4 the wall-jet separation of variables; R6 (9.76)
h₉₉ = 7.3319[…]^{1/3}, not 5.6152 (``free_jet_halfwidth(printed=True)``); R7 (9.56) needs the kinematic stress; R8 Ψ of (9.85)
has dimension m⁵/s³ (u²·ν·x), not a force per length.
Finding beyond the book: the wall-jet family f → λf(λη) makes C and f_∞ redundant (one physical constant); see
:func:`wall_jet_constants`.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import cumulative_trapezoid, quad, solve_bvp, solve_ivp
from scipy.optimize import brentq

from ._util import as_scalar_if_0d

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d

_C_FJ = 4.0 * np.sqrt(6.0) / 3.0  # Eq. (9.72), verified by quad in free_jet_constants
H99_PRINTED = 5.6152
"""The coefficient the book prints in Eq. (9.76) (it is the 4 % point, sech² = 0.04).  Label: book value, wrong."""

# ======================================================================================================================
# free jet
# ======================================================================================================================


def free_jet_constants() -> dict:
    """Numerical constants of the free jet, computed (sympy-exact values checked by quad).

    Book: §9.10, Eqs. (9.72), (9.73), (9.62), (9.71), (9.76): C = ∫f′²dη = ∫sech⁴(η/√6)dη = 4√6/3 (9.72); ṁ = (36Jρ²νx)^{1/3} (9.73);
    u₀ = (J²/C²ρ²νx)^{1/3} (9.62); f_∞ = √6; h₉₉: sech²z = 0.01 ⇒ z = arccosh 10.
    Returns dict(C [–], mdot_coeff = 36^{1/3}, u0_coeff = C^{−2/3} (Bickley 0.4543), xi_coeff = C^{−1/3}/√6 (Bickley 0.2752), f_inf, h99_arg,
    h99_coeff = √6 arccosh 10 (7.3319), h99_printed = √6 arccosh 5 = 5.6153, the 4 % point that the book's printed 5.6152 (:data:`H99_PRINTED`) really is).  Validation: V1/V2 quad vs 4√6/3; V5 Bickley (1937) coefficients
    via Wikipedia "Bickley jet".  Label: analytic, benchmark."""
    C_num = quad(lambda e: np.cosh(np.clip(e / np.sqrt(6.0), -300.0, 300.0)) ** -4, -np.inf, np.inf, epsabs=1e-13, epsrel=1e-13)[0]
    assert abs(C_num - _C_FJ) < 1e-10
    arg = float(np.arccosh(10.0))
    return dict(C=_C_FJ, C_quad=C_num, mdot_coeff=36.0 ** (1 / 3), u0_coeff=_C_FJ ** (-2 / 3), xi_coeff=_C_FJ ** (-1 / 3) / np.sqrt(6.0),
                f_inf=float(np.sqrt(6.0)), h99_arg=arg, h99_coeff=float(np.sqrt(6.0) * arg),
                h99_printed=float(np.sqrt(6.0) * np.arccosh(5.0)))  # √6 arccosh 5 = 5.6153: the 4 % point the book's "5.6152" actually is


def free_jet_profile(eta) -> dict:
    """Similarity profile of the free jet: f = √6 tanh(η/√6), f′ = sech²(η/√6), f″ = −(2/√6)sech²tanh  (Eq. (9.70), (9.71)).
    Parameters: eta [–] (η = y/δ(x), array or float).  Returns dict(f, fp, fpp) [–].
    Label: analytic (satisfies 3f‴ + ff″ + f′² = 0; sympy in ch09)."""
    e = _F(eta) / np.sqrt(6.0)
    t = np.tanh(e)
    s2 = 1.0 - t ** 2  # sech²
    return dict(f=_S(np.sqrt(6.0) * t), fp=_S(s2), fpp=_S(-2.0 / np.sqrt(6.0) * s2 * t))  # Eq. (9.70)


def free_jet_profile_table(n: int = 201, eta_max: float = 12.0) -> dict:
    """dict(eta, f, fp, fpp) on n points of [0, η_max] (for explainers)."""
    eta = np.linspace(0.0, eta_max, int(n))
    d = free_jet_profile(eta)
    return dict(eta=eta, f=d["f"], fp=d["fp"], fpp=d["fpp"])


def free_jet_centreline(x, J: float, rho: float, nu: float):
    """Centre-line speed u₀(x) = (J²/(C²ρ²νx))^{1/3} ∝ x^{−1/3}  [m/s].  Book: §9.10, Eq. (9.62)."""
    return _S((J ** 2 / (_C_FJ ** 2 * rho ** 2 * nu * _F(x))) ** (1 / 3))  # Eq. (9.62)


def free_jet_thickness(x, J: float, rho: float, nu: float):
    """Similarity thickness δ(x) = (Cρν²x²/J)^{1/3} ∝ x^{2/3}  [m].  Book: §9.10, Eq. (9.63)."""
    return _S((_C_FJ * rho * nu ** 2 * _F(x) ** 2 / J) ** (1 / 3))  # Eq. (9.63)


def free_jet(x, y, J: float, rho: float, nu: float) -> dict:
    """Velocity field of the plane laminar free jet (Bickley's solution).

    Book: §9.10, Eqs. (9.62)–(9.64), (9.70)–(9.71), (9.74): u = u₀(x) sech²(η/√6) (9.71); ψ = (Jνx/Cρ)^{1/3} f(η), f = √6 tanh(η/√6) (9.64),
    (9.70); v = −ψ_x = −(1/3)(Jν/Cρx²)^{1/3}[f − 2ηf′] (9.74); η = y/δ(x), δ (9.63).
    Parameters: x [m] > 0, y [m] (broadcast); J [N/m] momentum flux per unit span; rho [kg/m³]; nu [m²/s].
    Returns dict(u, v [m/s], psi [m²/s], eta [–], u0 [m/s], delta [m]).
    Assumptions: laminar, dp/dx = 0, x ≫ the slot width (inlet forgotten, Eq. (9.55)), boundary-layer approximation.
    Validation: V2 sympy PDE residual ≡ 0; V4 ∫u²dy = J/ρ for every x; V5 Bickley coefficients.  Label: symbolic, conserved, benchmark."""
    x, y = np.broadcast_arrays(_F(x), _F(y))
    u0 = free_jet_centreline(x, J, rho, nu)
    delta = free_jet_thickness(x, J, rho, nu)
    eta = y / delta
    pr = free_jet_profile(eta)
    f, fp = pr["f"], pr["fp"]
    u = u0 * fp  # Eq. (9.71)
    v = -(1.0 / 3.0) * (J * nu / (_C_FJ * rho * x ** 2)) ** (1 / 3) * (_F(f) - 2.0 * eta * _F(fp))  # Eq. (9.74)
    psi = (J * nu * x / (_C_FJ * rho)) ** (1 / 3) * _F(f)  # Eq. (9.64)
    return dict(u=_S(u), v=_S(v), psi=_S(psi), eta=_S(eta), u0=_S(u0), delta=_S(delta))


def free_jet_ode_solve(eta_max: float = 12.0, n: int = 400, coeff: float = 3.0, tol: float = 1e-10, fast: bool = False,
                       bc: str = "robin") -> dict:
    """Solve coeff·f‴ + ff″ + f′² = 0 (coeff = 3) with f(0) = 0 (9.68), f′(0) = 1 (9.67), f′(∞) = 0 (9.66) by ``solve_bvp`` — independent of the closed form.

    Book: §9.10, the reduced equation (unnumbered; the book skips the differentiation, sympy-checked in ch09) and (9.66)–(9.68).
    Parameters: eta_max [–] truncation of ∞; n initial mesh points; coeff [–] the coefficient of f‴ (3 for the jet; the closed form applies only to 3);
    tol solver tolerance; fast: coarser mesh; bc "robin" (default) applies the far-field condition f″ + (f/coeff)f′ = 0 at η_max (the
    linearisation of the equation for f′ → 0, whose solution decays like e^{−f_∞η/coeff}: an error of order f′² instead of f′) or "dirichlet"
    (f′(η_max) = 0, the book's condition (9.66); the far field decays only like e^{−2η/√6}, so this leaves ≈ 2e-4 at η_max = 12).
    Returns dict(eta, f, fp, max_err [–] vs √6 tanh(η/√6) and sech² (only meaningful for coeff = 3), eta_max, success).
    Validation: V3 max_err ≲ 1e-8 with the Robin condition at η_max = 12.  Label: converged.
    """
    rhs = lambda e, Y: np.vstack([Y[1], Y[2], -(Y[0] * Y[2] + Y[1] ** 2) / coeff])  # noqa: E731
    if bc == "robin":
        bcf = lambda a, b: np.array([a[0], a[1] - 1.0, b[2] + b[0] * b[1] / coeff])  # noqa: E731
    elif bc == "dirichlet":
        bcf = lambda a, b: np.array([a[0], a[1] - 1.0, b[1]])  # noqa: E731
    else:
        raise ValueError('bc must be "robin" or "dirichlet"')
    e = np.linspace(0.0, eta_max, 150 if fast else int(n))
    g = np.vstack([np.sqrt(6) * np.tanh(e / 3.0), 1 / np.cosh(e / 3.0) ** 2 * 0.9, -0.2 * np.exp(-e / 3)])  # crude start, not the answer
    sol = solve_bvp(rhs, bcf, e, g, tol=tol, max_nodes=200000)
    ee = np.linspace(0.0, eta_max, 2001)
    Y = sol.sol(ee)
    ex = free_jet_profile(ee)
    return dict(eta=ee, f=Y[0], fp=Y[1], max_err=float(max(np.max(np.abs(Y[0] - ex["f"])), np.max(np.abs(Y[1] - ex["fp"])))),
                eta_max=float(eta_max), success=bool(sol.success))


def free_jet_mass_flux(x, J: float, rho: float, nu: float):
    """Mass flux ṁ = ρu₀δ·2√6 = (36Jρ²νx)^{1/3} ∝ x^{1/3}  [kg/(m·s)] — the jet entrains ambient fluid.  Book: §9.10, Eq. (9.73)."""
    return _S((36.0 * J * rho ** 2 * nu * _F(x)) ** (1 / 3))  # Eq. (9.73)


def free_jet_entrainment_velocity(Re_x):
    """Entrainment speed v/u₀ → −√6/(3√Re_x) as η → +∞ (ambient fluid drawn toward the jet; +√6/(3√Re_x) at −∞); Re_x = xu₀/ν [–].
    Book: §9.10, Eq. (9.75)."""
    return _S(-np.sqrt(6.0) / (3.0 * np.sqrt(_F(Re_x))))  # Eq. (9.75)


def free_jet_at_level(level: float = 0.01) -> dict:
    """Where the free-jet speed has fallen to ``level``·u₀: sech²(z) = level with z = η/√6 ⇒ z = arccosh(1/√level), η = √6 z.

    Book: §9.10, Eqs. (9.71), (9.76).  Parameters: level [–] in (0, 1) (0.01 → the "99 %" half-width).
    Returns dict(z [–], coeff [–] = √6 z = η_level, the coefficient of (Cρν²x²/J)^{1/3} in the half-width): level 0.01 → 7.3319; 0.04 → 5.6153
    (the coefficient the book prints for 0.01); 0.5 → 2.1589.  Label: analytic."""
    z = float(np.arccosh(1.0 / np.sqrt(level)))
    return dict(z=z, coeff=float(np.sqrt(6.0) * z))


def free_jet_halfwidth(x, J: float, rho: float, nu: float, level: float = 0.01, printed: bool = False):
    """Half-width h at which u = level·u₀: h = η_level·(Cρν²x²/J)^{1/3} ∝ x^{2/3}  [m].

    Book: §9.10, Eq. (9.76) as printed: sech²z = 0.01 ⇒ z = 2.2924 ⇒ h₉₉ = 5.6152[Cρν²x²/J]^{1/3} — but arccosh 5 = 2.2924 is the 4 % point;
    the 1 % point is arccosh 10 = 2.9932, h₉₉ = 7.3319[…]^{1/3} (slip R6).  ``printed=True`` returns the book's coefficient (a named wrong
    variant that a test must fail).  Label: analytic."""
    coeff = H99_PRINTED if printed else free_jet_at_level(level)["coeff"]
    return _S(coeff * (_C_FJ * rho * nu ** 2 * _F(x) ** 2 / J) ** (1 / 3))  # Eq. (9.76)


def free_jet_reynolds(x, J: float, rho: float, nu: float) -> dict:
    """Jet Reynolds numbers Re_x = xu₀/ν = (3Jx/(4√6ρν²))^{2/3} and Re_h99 = u₀h₉₉/ν [–] with the CORRECT 7.3319 (and the book's 5.6152).

    Book: §9.10 (after Eq. (9.76)).  Re_h99 = coeff·(3Jx/(4√6ρν²))^{1/3}.  Returns dict(Re_x, Re_h99, Re_h99_printed)."""
    base = 3.0 * J * _F(x) / (4.0 * np.sqrt(6.0) * rho * nu ** 2)
    return dict(Re_x=_S(base ** (2 / 3)), Re_h99=_S(free_jet_at_level(0.01)["coeff"] * base ** (1 / 3)), Re_h99_printed=_S(H99_PRINTED * base ** (1 / 3)))


def jet_momentum_flux(y, u, rho: float = 1.0):
    """Momentum flux per unit span J = ρ∫u²dy  [N/m] (trapezoid, ``np.trapezoid``) — constant in x for the free jet.  Book: §9.10, Eqs. (9.57)–(9.58)."""
    return float(rho * np.trapezoid(_F(u) ** 2, _F(y)))  # Eq. (9.58)


# ======================================================================================================================
# wall jet
# ======================================================================================================================


def _wall_jet_g(eta, f_inf: float):
    """Solve Eq. (9.83) for q = 1 − g (accurate far out) at given η: returns (g, q)."""
    eta = np.atleast_1d(_F(eta))
    R0 = np.sqrt(3.0) * np.arctan(1.0 / np.sqrt(3.0))

    def T(s):
        q = np.exp(s)
        g = 1.0 - q
        return -s + np.sqrt(3.0) * np.arctan((2 * g + 1) / np.sqrt(3.0)) + 0.5 * np.log(1 + g + g * g)

    out_q = np.empty_like(eta)
    for i, e in enumerate(eta.ravel()):
        R = f_inf * e / 4.0 + R0
        if e <= 0:
            out_q[i] = 1.0
            continue
        s = brentq(lambda s_: T(s_) - R, -(R + 4.0), 0.0, xtol=1e-15, rtol=1e-15)
        out_q[i] = np.exp(s)
    return 1.0 - out_q, out_q


def wall_jet_profile(eta, f_inf: float = 1.0) -> dict:
    """Wall-jet similarity profile from the implicit solution (9.83), g² = f/f_∞.

    Book: §9.10, Eq. (9.83): −ln(1−g) + √3 tan⁻¹((2g+1)/√3) + ½ln(1+g+g²) = (f_∞/4)η + √3 tan⁻¹(1/√3), inverted by ``brentq`` in ln(1−g)
    (so 1 − g ~ e^{−f_∞η/4} stays representable); f = f_∞g², f′ = f_∞²g(1−g³)/6 (from g′ = f_∞(1−g³)/12, Eq. (9.83)'s integrand), f″ from the derivative.
    Parameters: eta [–] ≥ 0 (array or float); f_inf [–] the free scale (the ODE is invariant under f → λf(λη)).
    Returns dict(g, f, fp, fpp) [–] (floats for a scalar η, arrays otherwise).  Assumptions: corrected separation of variables (slip R4: the book's
    ∫df/(f_∞^{3/2}f − f²) lacks f^{1/2}).
    Validation: V2/V3 satisfies 4f‴ + ff″ + 2f′² = 0; agrees with the IVP to 1e-10; f″(0) = f_∞³/72.  Label: symbolic, converged.
    """
    scalar = np.ndim(eta) == 0
    g, q = _wall_jet_g(eta, f_inf)
    P = 1 + g + g * g
    f = f_inf * g ** 2
    gp = f_inf * q * P / 12.0
    fp = f_inf ** 2 * g * q * P / 6.0
    fpp = f_inf ** 2 / 6.0 * gp * (q * P - g * P + g * q * (1 + 2 * g))
    if scalar:
        return dict(g=float(g[0]), f=float(f[0]), fp=float(fp[0]), fpp=float(fpp[0]))
    return dict(g=g, f=f, fp=fp, fpp=fpp)


def wall_jet(x, y, C: float, f_inf: float, nu: float, rho: float = 1.0) -> dict:
    """Velocity field of the plane laminar wall jet.

    Book: §9.10, Eqs. (9.82), (9.83), (9.84) and v = −∂ψ/∂x = −√(νC)(f − 3ηf′)/(4x^{3/4}): u₀ = Cx^{−1/2}, δ = (νx^{3/2}/C)^{1/2} ∝ x^{3/4},
    ψ = (νCx^{1/2})^{1/2}f(η), u = u₀f′(η).
    Parameters: x [m] > 0, y [m] ≥ 0 (broadcast); C [m^{3/2}/s] (u₀ = Cx^{−1/2}); f_inf [–]; nu [m²/s]; rho [kg/m³] (only for ṁ).
    Returns dict(u, v [m/s], psi [m²/s], eta, u0, delta [m], mdot [kg/(m·s)] = ρ√(νC)f_∞x^{1/4}).
    Only the combination C f_∞^{…} is physical (one-parameter family) — see :func:`wall_jet_constants`.  Label: symbolic, conserved."""
    x, y = np.broadcast_arrays(_F(x), _F(y))
    u0 = C * x ** -0.5
    delta = np.sqrt(nu * x ** 1.5 / C)  # Eq. (9.82)
    eta = y / delta
    pr = wall_jet_profile(eta.ravel(), f_inf)
    f, fp = pr["f"].reshape(eta.shape), pr["fp"].reshape(eta.shape)
    u = u0 * fp
    v = -np.sqrt(nu * C) * (f - 3.0 * eta * fp) / (4.0 * x ** 0.75)
    psi = np.sqrt(nu * C * np.sqrt(x)) * f  # ψ = [νCx^{1/2}]^{1/2} f  Eq. (9.82)
    return dict(u=_S(u), v=_S(v), psi=_S(psi), eta=_S(eta), u0=_S(u0), delta=_S(delta), mdot=_S(rho * np.sqrt(nu * C) * f_inf * x ** 0.25))


def wall_jet_mass_flux(x, C: float, f_inf: float, rho: float, nu: float):
    """ṁ = ρ√(νC) f_∞ x^{1/4} ∝ x^{1/4}  [kg/(m·s)] (entrainment weaker than the free jet's x^{1/3}).  Book: §9.10, Eq. (9.84)."""
    return _S(rho * np.sqrt(nu * C) * f_inf * _F(x) ** 0.25)  # Eq. (9.84)


def wall_jet_first_integral_residual(f, fp, fpp):
    """4ff″ − 2f′² + f²f′ — zero along the wall-jet solution (first integral after multiplying 4f‴ + ff″ + 2f′² = 0 by f; constant 0 from η = 0).
    Book: §9.10 (unnumbered, above Eq. (9.83))."""
    f, fp, fpp = _F(f), _F(fp), _F(fpp)
    return _S(4.0 * f * fpp - 2.0 * fp ** 2 + f ** 2 * fp)


def wall_jet_ode_solve(fpp0: float = 1.0, eta_max: float = 40.0, coeff: float = 4.0, n: int = 2001, printed: bool = False) -> dict:
    """Integrate coeff·f‴ + ff″ + 2f′² = 0 (coeff = 4, correct, slip R3) from f(0) = f′(0) = 0, f″(0) = fpp0 with ``solve_ivp`` (DOP853, rtol 1e-12), and
    compare with the closed form (9.83) evaluated at the resulting f_∞.

    Book: §9.10 below Eq. (9.82) (the book prints f‴ + ff″ + 2f′² = 0; ``coeff=1.0`` (or ``printed=True``) integrates that coefficient — a named wrong
    variant, whose solution does NOT satisfy (9.83)), (9.77)–(9.78).  Derived: f_∞³ = 72 f″(0) (from f′ ≈ f^{1/2}f_∞^{3/2}/6 near the wall) and the
    scaling f_∞ ∝ f″(0)^{1/3}.
    Parameters: fpp0 [–] f″(0) (free scale); eta_max [–] truncation of ∞; coeff [–] the coefficient of f‴; n output points; printed [bool] shorthand for coeff = 1.
    Returns dict(eta, f, fp, fpp, f_inf, err_vs_9_83 [–], fpp0_over_finf_cubed [–] (= 1/72 for the correct ODE), first_integral_residual_max).
    Validation: V2/V3/V7.  Label: converged, symbolic.
    """
    k = 1.0 if printed else float(coeff)
    rhs = lambda e, Y: [Y[1], Y[2], -(Y[0] * Y[2] + 2.0 * Y[1] ** 2) / k]  # noqa: E731
    eta = np.linspace(0.0, eta_max, int(n))
    sol = solve_ivp(rhs, (0.0, eta_max), [0.0, 0.0, float(fpp0)], method="DOP853", rtol=1e-12, atol=1e-14, t_eval=eta)
    f, fp, fpp = sol.y
    f_inf = float(f[-1])
    ex = wall_jet_profile(eta, f_inf)
    return dict(eta=eta, f=f, fp=fp, fpp=fpp, f_inf=f_inf, err_vs_9_83=float(max(np.max(np.abs(f - ex["f"])), np.max(np.abs(fp - ex["fp"])))),
                fpp0_over_finf_cubed=float(fpp0 / f_inf ** 3), fpp0_over_f_inf_cubed=float(fpp0 / f_inf ** 3),
                first_integral_residual_max=float(np.max(np.abs(wall_jet_first_integral_residual(f, fp, fpp)))))


def wall_jet_invariant(y, u):
    """The wall-jet invariant ∫₀^∞ u(∫_y^∞u²dy′)dy  [m⁵/s³] — independent of x for the similarity solution.

    Book: §9.10, Eq. (9.80): d/dx∫₀^∞ u(∫_y^∞ u²dy′)dy = 0 (the "flux of exterior momentum flux"); with u = u₀f′ it equals C²ν∫f′(∫_η^∞f′²)dη (9.85).
    Inner tail integral by reversed ``cumulative_trapezoid``, outer by ``np.trapezoid`` (second order).  y [m] increasing from 0; u [m/s].
    Validation: V4 constant in x to 1e-9 on the fields of :func:`wall_jet`.  Label: conserved."""
    y, u = _F(y), _F(u)
    tail = cumulative_trapezoid((u ** 2)[::-1], -y[::-1], initial=0.0)[::-1]  # ∫_y^∞ u² dy′
    return float(np.trapezoid(u * tail, y))  # Eq. (9.80)


def wall_jet_K1(numeric: bool = False) -> float:
    """K₁ = ∫₀^∞ f′(∫_η^∞f′²dη′)dη for f_∞ = 1 [–]; for general f_∞ the double integral is K₁f_∞⁴ (scaling f → λf(λη)).
    Book: §9.10, Eq. (9.85).  Closed form (ours): integrate by parts, K₁ = ∫f f′²dη, and with g (f = g², dη = 12dg/(1−g³)) K₁ = ∫₀¹(g⁴ − g⁷)/3 dg = 1/40.
    ``numeric=True`` recomputes it from the profile by ``np.trapezoid`` on a fine grid (agrees to 1e-6).  Label: analytic."""
    if not numeric:
        return 1.0 / 40.0
    eta = np.linspace(0.0, 60.0, 6001)
    fp = wall_jet_profile(eta, 1.0)["fp"]
    return wall_jet_invariant(eta, fp)


def wall_jet_integrals(f_inf: float = 1.0) -> dict:
    """The integrals that fix the wall-jet constants, by the scaling f → λf(λη) (λ = f_∞ from the f_∞ = 1 solution).

    Book: §9.10, Eqs. (9.77)–(9.85): ∫f′dη = f_∞, ∫f′²dη = f_∞³/18, ∫₀^∞f′(∫_η^∞f′²dη′)dη = f_∞⁴/40 (the K of (9.85)), f″(0) = f_∞³/72.
    The f_∞ = 1 values are computed here in the variable g (f = g², dη = 12dg/(1−g³), f′ = g(1−g³)/6) by ``scipy.integrate.quad`` and then scaled
    (f_∞ multiplies f, f_∞^{−1} multiplies η): ∫f′ dη scales as f_∞, ∫f′² as f_∞³, the double integral as f_∞⁴ and f″(0) as f_∞³.
    Parameters: f_inf [–] the free scale.  Returns dict(int_fp, int_fp2, invariant, fpp0) [–] — for f_∞ = 1: (1, 1/18 = 0.055556, 1/40 = 0.025, 1/72 = 0.013889).
    Validation: V1 exact rationals via quad; V4 the double integral agrees with :func:`wall_jet_invariant` on the profile.  Label: analytic.
    """
    # f_∞ = 1:  dη = 12 dg/(1 − g³), f′ = g(1−g³)/6  ⇒  f′dη = 2g dg;  f′² dη = f′·(2g dg)
    i1 = quad(lambda g: 2.0 * g, 0.0, 1.0, epsabs=1e-14)[0]
    i2 = quad(lambda g: g * (1.0 - g ** 3) / 6.0 * 2.0 * g, 0.0, 1.0, epsabs=1e-14)[0]
    # ∫f′(∫_η^∞ f′² dη′)dη = ∫ f′² · f dη  (integration by parts, f(0) = 0)  = ∫ (g(1−g³)/6)² g² · 12/(1−g³) dg
    i3 = quad(lambda g: (g * (1.0 - g ** 3) / 6.0) ** 2 * g ** 2 * 12.0 / (1.0 - g ** 3), 0.0, 1.0, epsabs=1e-14)[0]
    fpp0 = wall_jet_profile(0.0, 1.0)["fpp"]  # f″(0) of the (9.83) profile at f_∞ = 1 (= 1/72, derived in the docstring's list)
    F = float(f_inf)
    return dict(int_fp=i1 * F, int_fp2=i2 * F ** 3, invariant=i3 * F ** 4, fpp0=fpp0 * F ** 3)


def wall_jet_constants(rho: float, nu: float, *, Psi: float | None = None, mdot: float | None = None, x: float | None = None,
                       f_inf: float = 1.0) -> dict:
    """The wall-jet constant C from ONE physical datum: the invariant Ψ (9.85) or the mass flux ṁ at a station x (9.84).

    Book: §9.10, Eq. (9.85) Ψ = C²ν∫f′∫f′² = C²ν f_∞⁴/40 and Eq. (9.84) ṁ = ρ√(νC) f_∞ x^{1/4}.  FINDING (ours): the family f → λf(λη) makes (C, f_∞)
    redundant — only C f_∞² is physical — so Ψ and ṁ carry the SAME information, related by Ψ = ṁ⁴/(40ρ⁴νx).  Hence: give Ψ, or (ṁ and x), plus the
    gauge f_∞ (free scale, default 1) and get C; the other datum follows.  Giving all of Ψ, ṁ and x returns their consistency residual.
    Parameters: rho [kg/m³]; nu [m²/s]; Psi [m⁵/s³] (keyword); mdot [kg/(m·s)] (keyword, with x); x [m] station of ṁ; f_inf [–] gauge.
    Returns dict(C [m^{3/2}/s] (u₀ = Cx^{−1/2}), C_f_inf_sq [m^{3/2}/s] = C f_∞² (gauge-free), Psi [m⁵/s³], mdot_at_x [kg/(m·s)] (needs x, else None),
    residual [–] = Ψ_given/Ψ_from_ṁ − 1 (only if all three are given, else None)).
    Validation: V1 identity Ψ = ṁ⁴/(40ρ⁴νx) on the fields of :func:`wall_jet`; V4 the invariant.  Label: analytic.
    """
    K = wall_jet_integrals(1.0)["invariant"]  # 1/40
    if Psi is None and mdot is None:
        raise ValueError("give Psi, or mdot together with x")
    if mdot is not None and x is None:
        raise ValueError("mdot is a mass flux AT a station: give x")
    res = None
    if mdot is not None:
        Cf2 = (mdot / rho) ** 2 / (nu * np.sqrt(x))  # C f_∞² from (9.84): ṁ² = ρ²νC f_∞² x^{1/2}
        Psi_m = K * nu * Cf2 ** 2  # (9.85): Ψ = ν(C f_∞²)²/40 = ṁ⁴/(40ρ⁴νx)
        if Psi is not None:
            res = float(Psi / Psi_m - 1.0)
        else:
            Psi = float(Psi_m)
    else:
        Cf2 = float(np.sqrt(Psi / (K * nu)))  # from (9.85)
    mdot_at_x = None if x is None else float(rho * np.sqrt(nu * Cf2) * x ** 0.25)  # (9.84) with C f_∞² = Cf2
    return dict(C=float(Cf2 / f_inf ** 2), C_f_inf_sq=float(Cf2), Psi=float(Psi), mdot_at_x=mdot_at_x, residual=res)


__all__ = [n for n in dir() if not n.startswith("_") and n not in ("annotations",)]
