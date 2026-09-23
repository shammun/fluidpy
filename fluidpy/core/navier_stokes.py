"""The differential conservation laws as residuals and term splits for any candidate field (numeric stencils + sympy).

Book: Kundu, Cohen & Dowling 5e, Ch. 4 — continuity (4.7)–(4.10); momentum (4.22)–(4.24); Navier–Stokes (4.38)–(4.41);
energy (4.53)–(4.60); entropy (4.61)–(4.63); Lamb's identity (4.68); gravity absorbed in p′ (4.84)–(4.85);
Boussinesq (4.86)–(4.89). Every equation was transcribed from the rendered pages chapters/pages/ch04/p125–p164.

This module is the **verifier for every exact solution of later chapters**: plug a velocity/pressure/density field
into a residual and it must vanish to the stencil's O(h²) error.

Conventions
-----------
* Fields are callables ``F(x, t)`` with ``x`` (d,) or (d, N) (coordinates on axis 0; ``core.kinematics``). Scalars
  return () / (N,), vectors (d,) / (d, N), tensors (d, d) / (d, d, N). Constants may be plain numbers.
* Stencils: explicit 2nd-order central differences, length step ``h`` [m] (≪ the flow's length scale) and time step
  ``ht`` [s] (``None`` → 1e-4 s; scale it with the flow's time scale). Nested stencils keep O(h²).
* τ_ij: first index = face normal (§2.4); the net surface force ∂τ_ij/∂x_i in (4.20)–(4.24) contracts the **first**
  index (the book's prose after (4.24), and again in §4.8, writes ∂τ_ij/∂x_j — harmless only because τ is symmetric,
  which (4.25) proves afterwards; analysis §9 item 8). A non-symmetric test τ discriminates.
* z is up and g = (0, 0, −g), g = 9.81 m/s² (the book's value). For 2-D fields in the (x, y) plane the default g has
  no in-plane component (a horizontal plane); pass ``g=(0, -9.81)`` for a vertical plane.
* ``sigma`` arguments are viscous-stress fields σ(x, t) → (d, d, …) [Pa]; :func:`newtonian_viscous_stress_field`
  builds the Newtonian one (4.59) from a velocity field.
"""
from __future__ import annotations

from typing import Callable, NamedTuple

import numpy as np
import sympy as sp
from scipy.special import erfc, exp1

from ._stencil import as_field, ddt, ddx, div, ev, grad, grad_vector, gvec, laplacian, pad3
from ._util import as_scalar_if_0d
from .constitutive import viscous_stress
from .kinematics import acceleration, material_derivative, vorticity
from .statics import integrate_hydrostatic
from .thermo import G_BOOK as G0

__all__ = [
    # continuity
    "continuity_residual", "continuity_residual_sym", "ContinuityTerms", "continuity_terms",
    "ContinuityMaterialTerms", "continuity_material_terms", "density_material_rate", "divergence_free_check",
    # momentum
    "stress_divergence", "divergence_first_index_demo", "momentum_conservative_residual",
    "conservative_to_advective_sym", "CauchyTerms", "cauchy_terms", "cauchy_residual", "newtonian_stress_field",
    "newtonian_viscous_stress_field", "navier_stokes_residual", "navier_stokes_sym", "NSTerms",
    "ns_incompressible_terms", "viscous_force_forms",
    # exact solutions
    "EXACT_SOLUTIONS", "exact_solution", "exact_field", "exact_solution_fields", "ns_terms_preset",
    # Lamb identity
    "lamb_vector", "lamb_identity_terms", "lamb_identity_sym",
    # energy / entropy
    "stress_work_split", "kinetic_energy_budget", "internal_energy_terms", "internal_energy_residual",
    "total_energy_residual_sym", "energy_forms_sym", "energy_identities_sym", "entropy_terms", "entropy_production",
    # §4.9 gravity, Boussinesq
    "perturbation_fields", "buoyancy", "boussinesq_momentum_terms", "heat_equation_terms",
    "temperature_equation_residual",
]

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d


def _pts(x) -> tuple:
    return _F(x).shape[1:]


def _ufield(u: Callable, d: int) -> Callable:
    return lambda X, T: ev(u, X, T, (d,))


def _zero_tensor(d: int) -> Callable:
    return lambda X, T: np.zeros((d, d) + _F(X).shape[1:])


# ======================================================================================================================
# §4.2 continuity
# ======================================================================================================================
def continuity_residual(rho, u: Callable, x, t: float = 0.0, h: float = 1e-4, ht: float | None = None):
    """Residual of the continuity equation ∂ρ/∂t + ∂(ρu_i)/∂x_i, Eq. (4.7) — zero for a mass-conserving field.

    Book: §4.2, Eq. (4.7) (from (4.2) by Gauss (2.30) and the vanishing-integrand argument of (4.6)).

    Parameters
    ----------
    rho : density field ρ(x, t) [kg/m³] (callable or constant);  u : velocity field u(x, t) [m/s]
    x : point(s) [m], (d,) or (d, N);  t : time [s];  h : space step [m];  ht : time step [s] (None → 1e-4 s)

    Returns
    -------
    residual [kg/(m³ s)] — float for one point, (N,) otherwise.

    Assumptions: continuum, smooth fields.
    Validation: V1 the expanding flow u = ax/(1 + at), ρ = ρ₀(1 + at)^(−d) gives 0 and a wrong ρ does not
    (discrimination); V2 equals :func:`continuity_residual_sym`; V3 stencil order 2. Label: analytic, converged.
    """
    return continuity_terms(rho, u, x, t, h, ht).residual


def continuity_residual_sym(rho_expr, u_exprs, coords, t: sp.Symbol):
    """Symbolic Eq. (4.7): ∂ρ/∂t + ∂(ρu_i)/∂x_i for sympy expressions (returned unsimplified).

    Book: §4.2, Eq. (4.7). Label: symbolic.
    """
    rho_e = sp.sympify(rho_expr)
    return sp.diff(rho_e, t) + sum(sp.diff(rho_e * sp.sympify(ui), xi) for ui, xi in zip(u_exprs, coords))  # (4.7)


class ContinuityTerms(NamedTuple):
    """The two terms of (4.7) [kg/(m³ s)]: ``local`` ∂ρ/∂t + ``flux_divergence`` ∇·(ρu) = ``residual`` (0).
    Book: §4.2, Eq. (4.7). Label: analytic.
    """

    local: float | np.ndarray
    flux_divergence: float | np.ndarray
    residual: float | np.ndarray


def continuity_terms(rho, u: Callable, x, t: float = 0.0, h: float = 1e-4, ht: float | None = None) -> ContinuityTerms:
    """∂ρ/∂t and ∇·(ρu) separately — the two terms of the continuity equation (4.7) (C02's mirror-image curves).

    Book: §4.2, Eq. (4.7); the flux-divergence term is "the net loss at a point due to divergence of a flux".
    Parameters: as :func:`continuity_residual`. Returns :class:`ContinuityTerms`. Label: analytic, converged.
    """
    rho_f = as_field(rho)
    x_ = _F(x)
    d = x_.shape[0]
    local = ddt(rho_f, x_, t, ht)  # ∂ρ/∂t
    flux = div(lambda X, T: ev(rho_f, X, T) * ev(u, X, T, (d,)), x_, t, h)  # ∂(ρu_i)/∂x_i
    return ContinuityTerms(_S(local), _S(flux), _S(local + flux))  # Eq. (4.7)


class ContinuityMaterialTerms(NamedTuple):
    """Eq. (4.8) [1/s]: ``rel_rate`` (1/ρ)Dρ/Dt + ``div_u`` ∇·u = 0. Book: §4.2, Eq. (4.8). Label: analytic."""

    rel_rate: float | np.ndarray
    div_u: float | np.ndarray


def continuity_material_terms(rho, u: Callable, x, t: float = 0.0, h: float = 1e-4,
                              ht: float | None = None) -> ContinuityMaterialTerms:
    """The two terms of the advective form of continuity, (1/ρ)Dρ/Dt + ∇·u = 0, Eq. (4.8).

    Book: §4.2, Eq. (4.8) (from (4.7) with ∂(ρu_i)/∂x_i = u_i∂ρ/∂x_i + ρ∂u_i/∂x_i and (3.5)). ∇·u is the volumetric
    strain rate (3.14) — the book's text says "Section 3.6"; it is §3.4 (typo noted in analysis §9 item 9).

    Parameters: as :func:`continuity_residual`. Returns :class:`ContinuityMaterialTerms` (unpacks as a pair) [1/s].
    Validation: V1 sum = residual(4.7)/ρ pointwise; incompressible stratified flow: both 0 with ∇ρ ≠ 0.
    Label: analytic.
    """
    rho_f = as_field(rho)
    x_ = _F(x)
    d = x_.shape[0]
    rel = material_derivative(lambda X, T: ev(rho_f, X, T), u, x_, t, h, ht) / ev(rho_f, x_, t)  # (1/ρ) Dρ/Dt
    dv = div(_ufield(u, d), x_, t, h)  # ∇·u
    return ContinuityMaterialTerms(_S(rel), _S(dv))  # Eq. (4.8)


def density_material_rate(rho, u: Callable, x, t: float = 0.0, h: float = 1e-4, ht: float | None = None):
    """Dρ/Dt = ∂ρ/∂t + u·∇ρ, Eq. (4.9) — zero in *incompressible* flow even when ρ varies in space.

    Book: §4.2, Eq. (4.9) (the incompressibility condition). Units kg/(m³ s).
    Note (analysis §9): §4.11 reuses the label (4.9) for the identity ∇·u = −(1/ρ)Dρ/Dt = −(1/ρc²)Dp/Dt, which is
    general (isentropic), not incompressible.
    Validation: V1 stratified shear flow ρ = ρ₀ + ρ_z z, u = (U₀ + sz, 0, 0): Dρ/Dt = 0 although ∂ρ/∂z ≠ 0; V7 a
    Galilean shift of U₀ leaves it 0. Label: analytic.
    """
    return _S(material_derivative(lambda X, T: ev(as_field(rho), X, T), u, x, t, h, ht))  # Eq. (4.9)


def divergence_free_check(u: Callable, x, t: float = 0.0, h: float = 1e-4, tol: float = 1e-8):
    """Check ∇·u = 0, Eq. (4.10), at the given points; returns (ok, max|∇·u| [1/s]).

    Book: §4.2, Eq. (4.10) (incompressible flow: (4.8) + (4.9)). Label: analytic.
    Validation: V1 cylinder potential flow and every incompressible exact solution pass; the expanding flow does not.
    """
    x_ = _F(x)
    dv = np.abs(_F(div(_ufield(u, x_.shape[0]), x_, t, h)))
    m = float(np.max(dv))
    return bool(m <= tol), m


# ======================================================================================================================
# §4.4 momentum: flux form, Cauchy
# ======================================================================================================================
def stress_divergence(tau: Callable, x, t: float = 0.0, h: float = 1e-4, index: int = 0) -> np.ndarray:
    """Net surface force per unit volume ∂τ_ij/∂x_i (index=0, the book's (4.20b)–(4.24)) or ∂τ_ij/∂x_j (index=1).

    Book: §4.4, Eqs. (4.20b), (4.22), (4.24): ∫ n_iτ_ij dA = ∫ ∂τ_ij/∂x_i dV — the **first** index (the face normal)
    is contracted. ``core.operators.tensor_divergence`` (grid version) defaults to the second index; this point-wise
    version defaults to the first. The two agree only for symmetric τ (4.25).

    Parameters: tau : tensor field τ(x, t) [Pa] ((d, d) per point);  x [m];  t [s];  h [m];  index : 0 or 1.
    Returns (d,) + points [N/m³]. Validation: V1 non-symmetric linear τ: index 0 and 1 differ as predicted; V3 order 2.
    Label: analytic, converged.
    """
    x_ = _F(x)
    d = x_.shape[0]
    tf = lambda X, T: ev(tau, X, T, (d, d))  # noqa: E731
    if index == 0:
        return sum(ddx(lambda X, T, i=i: tf(X, T)[i], x_, t, i, h, (d,)) for i in range(d))  # ∂τ_ij/∂x_i
    return sum(ddx(lambda X, T, j=j: tf(X, T)[:, j], x_, t, j, h, (d,)) for j in range(d))  # ∂τ_ij/∂x_j


def divergence_first_index_demo(tau_fn: Callable, x, h: float = 1e-4, t: float = 0.0):
    """(∂τ_ij/∂x_i, ∂τ_ij/∂x_j) at a point — the book's first-index net surface force (4.20b)/(4.24) and the
    second-index version; they differ when τ is not symmetric (the first-index primer's numeric demo).

    Book: §4.4, Eqs. (4.20b), (4.24) (and the prose slip ∂τ_ij/∂x_j, analysis §9 item 8). ``tau_fn(x, t)`` → (d, d, …)
    [Pa]. Returns (first, second) [N/m³]. Label: analytic.
    """
    return stress_divergence(tau_fn, x, t, h, 0), stress_divergence(tau_fn, x, t, h, 1)


def momentum_conservative_residual(rho, u: Callable, tau: Callable, g=(0.0, 0.0, -G0), x=None, t: float = 0.0,
                                   h: float = 1e-4, ht: float | None = None) -> np.ndarray:
    """Residual of the flux (conservative) form of momentum: ∂(ρu_j)/∂t + ∂(ρu_iu_j)/∂x_i − ρg_j − ∂τ_ij/∂x_i, (4.22).

    Book: §4.4, Eq. (4.22) (from (4.14) with Gauss (4.20a, b) and the vanishing integrand (4.21)).

    Parameters
    ----------
    rho [kg/m³] (callable or constant); u [m/s]; tau : stress field [Pa]; g : body force per mass [m/s²];
    x : point(s) [m]; t [s]; h [m]; ht [s]

    Returns (d,) + points [N/m³].
    Validation: V2 :func:`conservative_to_advective_sym` (flux form = ρDu/Dt + u_j × continuity); V1 rigid rotation with
    τ = −p(r)δ, p = ρΩ²r²/2 gives 0; V3 order 2. Label: analytic, symbolic, converged.
    """
    rho_f = as_field(rho)
    x_ = _F(x)
    d = x_.shape[0]
    rho_u = lambda X, T: ev(rho_f, X, T) * ev(u, X, T, (d,))  # noqa: E731
    local = ddt(rho_u, x_, t, ht, (d,))  # ∂(ρu_j)/∂t
    flux = sum(ddx(lambda X, T, i=i: rho_u(X, T)[i] * ev(u, X, T, (d,)), x_, t, i, h, (d,))
               for i in range(d))  # ∂(ρu_i u_j)/∂x_i
    body = ev(rho_f, x_, t) * gvec(g, d, _pts(x_))
    return local + flux - body - stress_divergence(tau, x_, t, h, 0)  # Eq. (4.22)


def conservative_to_advective_sym(rho_expr, u_exprs, coords, t: sp.Symbol) -> list:
    """The step (4.23): [∂(ρu_j)/∂t + ∂(ρu_iu_j)/∂x_i] − ρDu_j/Dt, simplified — it equals u_j × (4.7)'s left side.

    Book: §4.4, Eq. (4.23) (the bracket u_j[∂ρ/∂t + ∂(ρu_i)/∂x_i] is zero by (4.7), turning (4.22) into (4.24)).
    Returns a list (one per j) of simplified expressions; ``[e − u_j·C for …]`` simplifies to 0 with C =
    :func:`continuity_residual_sym`. Label: symbolic.
    """
    rho_e = sp.sympify(rho_expr)
    u = [sp.sympify(e) for e in u_exprs]
    flux = [sp.diff(rho_e * uj, t) + sum(sp.diff(rho_e * ui * uj, xi) for ui, xi in zip(u, coords)) for uj in u]
    adv = [rho_e * (sp.diff(uj, t) + sum(ui * sp.diff(uj, xi) for ui, xi in zip(u, coords))) for uj in u]
    return [sp.simplify(sp.expand(f - a)) for f, a in zip(flux, adv)]  # Eq. (4.23): = u_j × continuity


class CauchyTerms(NamedTuple):
    """Cauchy's equation (4.24) term by term [N/m³]: ``inertia`` ρDu_j/Dt = ``body`` ρg_j + ``stress_divergence``
    ∂τ_ij/∂x_i; ``residual`` = inertia − body − stress_divergence. Book: §4.4, Eq. (4.24). Label: analytic.
    """

    inertia: np.ndarray
    body: np.ndarray
    stress_divergence: np.ndarray
    residual: np.ndarray


def cauchy_terms(rho, u: Callable, tau: Callable, g=(0.0, 0.0, -G0), x=None, t: float = 0.0, h: float = 1e-4,
                 ht: float | None = None) -> CauchyTerms:
    """Cauchy's equation of motion ρDu_j/Dt = ρg_j + ∂τ_ij/∂x_i, Eq. (4.24), term by term (per unit volume).

    Book: §4.4, Eq. (4.24) — true for any continuum (no constitutive law yet); with (4.7) it has 1 + 3 = 4 equations
    for 13 unknowns (ρ, u_j, τ_ij) (``ch04.closure_count("cauchy")``).

    Parameters: as :func:`momentum_conservative_residual`. Returns :class:`CauchyTerms` [N/m³].
    Validation: V1 rigid rotation (τ = −ρΩ²r²/2 δ) → residual 0; a non-symmetric τ pins the first index (the
    second-index variant fails); V3 order 2. Label: analytic, converged.
    """
    rho_f = as_field(rho)
    x_ = _F(x)
    d = x_.shape[0]
    r = ev(rho_f, x_, t)
    inertia = r * acceleration(u, x_, t, h, ht).a  # ρ Du_j/Dt
    body = r * gvec(g, d, _pts(x_))  # ρ g_j
    surf = stress_divergence(tau, x_, t, h, 0)  # ∂τ_ij/∂x_i (first index)
    return CauchyTerms(inertia, body, surf, inertia - body - surf)  # Eq. (4.24)


def cauchy_residual(rho, u: Callable, tau: Callable, g=(0.0, 0.0, -G0), x=None, t: float = 0.0, h: float = 1e-4,
                    ht: float | None = None) -> np.ndarray:
    """Residual ρDu_j/Dt − ρg_j − ∂τ_ij/∂x_i of Eq. (4.24) [N/m³]. Book: §4.4, Eq. (4.24). Label: analytic, converged."""
    return cauchy_terms(rho, u, tau, g, x, t, h, ht).residual


# ======================================================================================================================
# §4.6 Navier–Stokes
# ======================================================================================================================
def newtonian_viscous_stress_field(u: Callable, mu=1e-3, mu_v=0.0, h: float = 1e-4) -> Callable:
    """Viscous stress field σ(x, t) of a Newtonian fluid (4.59) built from a velocity field (gradient by stencil).

    Book: §4.8, Eq. (4.59) (= the viscous part of (4.37)). ``mu``, ``mu_v`` constants or callables (x, t).
    Returns a callable (x, t) → (d, d) + points [Pa]. Label: analytic.
    """
    muf, mvf = as_field(mu), as_field(mu_v)

    def sigma(X, T):
        X_ = _F(X)
        G = grad_vector(u, X_, T, h)
        return viscous_stress(G, ev(muf, X_, T), ev(mvf, X_, T))
    return sigma


def newtonian_stress_field(u: Callable, p, mu, mu_v=0.0, h: float = 1e-4) -> Callable:
    """Stress field τ(x, t) = −pδ + σ(∇u) of a Newtonian fluid (4.37) built from a velocity field (gradient by stencil).

    Book: §4.5, Eq. (4.37). ``p``, ``mu``, ``mu_v`` constants or callables (x, t) (variable viscosity for (4.38)).
    Returns a callable (x, t) → (d, d) + points [Pa]. Label: analytic.
    """
    pf = as_field(p)
    sig = newtonian_viscous_stress_field(u, mu, mu_v, h)

    def tau(X, T):
        X_ = _F(X)
        d = X_.shape[0]
        return -ev(pf, X_, T) * np.eye(d).reshape((d, d) + (1,) * (X_.ndim - 1)) + sig(X_, T)
    return tau


def navier_stokes_residual(rho, u: Callable, p, x, t: float = 0.0, mu=1e-3, mu_v=0.0, g=(0.0, 0.0, -G0),
                           h: float = 1e-4, ht: float | None = None) -> np.ndarray:
    """Residual of the compressible Navier–Stokes equation with variable viscosity, Eq. (4.38).

    Book: §4.6, Eq. (4.38): ρ(∂u_j/∂t + u_i∂u_j/∂x_i) = −∂p/∂x_j + ρg_j + ∂/∂x_i[μ(∂u_j/∂x_i + ∂u_i/∂x_j)
    + (μ_v − ⅔μ)(∂u_m/∂x_m)δ_ij] — (4.37) substituted into Cauchy's (4.24).

    Parameters
    ----------
    rho [kg/m³]; u [m/s]; p [Pa]; mu, mu_v [Pa s] (constants or callables (x, t)); g [m/s²]; x [m]; t [s]; h [m];
    ht [s]

    Returns (d,) + points [N/m³] — ρDu/Dt − (−∇p + ρg + ∂σ_ij/∂x_i).
    Validation: V1 plane Couette, plane and pipe Poiseuille, Stokes' first problem, Taylor–Green, Lamb–Oseen (the
    :data:`EXACT_SOLUTIONS` set) give ≈ 0; the wrong sign of the viscous term fails; V2 :func:`navier_stokes_sym`
    (4.38) − (4.39a) ≡ 0 for constant μ, μ_v; V3 order 2. Label: analytic, symbolic, converged.
    """
    tau = newtonian_stress_field(u, p, mu, mu_v, h)
    return cauchy_residual(rho, u, tau, g, x, t, h, ht)  # (4.38) = (4.24) with (4.37)


def navier_stokes_sym(rho, u_exprs, p_expr, coords, t: sp.Symbol, mu, mu_v=0, g=None, form: str = "4.38") -> list:
    """Symbolic Navier–Stokes residual, component by component, in one of the book's forms.

    Book: §4.6 — "4.38": ρDu_j/Dt − [−∂p/∂x_j + ρg_j + ∂/∂x_i(μ(∂u_j/∂x_i + ∂u_i/∂x_j) + (μ_v − ⅔μ)∂u_m/∂x_mδ_ij)];
    "4.39a": ρDu_j/Dt − [−∂p/∂x_j + ρg_j + μ∂²u_j/∂x_i² + (μ_v + ⅓μ)∂(∂u_m/∂x_m)/∂x_j] (constant μ, μ_v);
    "4.39b": ρDu_j/Dt − [−∂p/∂x_j + ρg_j + μ∇²u_j] (incompressible). μ, μ_v may be sympy expressions; g a list
    (default zeros).

    Returns a list of sympy expressions (unsimplified). Label: symbolic.
    """
    rho_e = sp.sympify(rho)
    u = [sp.sympify(e) for e in u_exprs]
    p_e = sp.sympify(p_expr)
    mu_e, mv_e = sp.sympify(mu), sp.sympify(mu_v)
    n = len(coords)
    g_ = [0] * n if g is None else [sp.sympify(v) for v in g]
    divu = sum(sp.diff(u[m], coords[m]) for m in range(n))
    out = []
    for j in range(n):
        acc = rho_e * (sp.diff(u[j], t) + sum(u[i] * sp.diff(u[j], coords[i]) for i in range(n)))
        if form == "4.38":
            visc = sum(sp.diff(mu_e * (sp.diff(u[j], coords[i]) + sp.diff(u[i], coords[j]))
                               + (mv_e - sp.Rational(2, 3) * mu_e) * divu * (1 if i == j else 0), coords[i])
                       for i in range(n))  # Eq. (4.38)
        elif form == "4.39a":
            visc = (mu_e * sum(sp.diff(u[j], coords[i], 2) for i in range(n))
                    + (mv_e + mu_e / 3) * sp.diff(divu, coords[j]))  # Eq. (4.39a)
        elif form == "4.39b":
            visc = mu_e * sum(sp.diff(u[j], coords[i], 2) for i in range(n))  # Eq. (4.39b)
        else:
            raise ValueError("form must be '4.38', '4.39a' or '4.39b'")
        out.append(acc - (-sp.diff(p_e, coords[j]) + rho_e * g_[j] + visc))
    return out


class NSTerms(NamedTuple):
    """The five terms of the incompressible Navier–Stokes equation (4.39b): ``local`` ∂u/∂t + ``advective`` (u·∇)u =
    ``pressure`` −∇p/ρ + ``gravity`` g + ``viscous`` ν∇²u; ``residual`` = local + advective − pressure − gravity −
    viscous. Per unit mass [m/s²] (or ×ρ, per unit volume [N/m³]). Book: §4.6, Eq. (4.39b). Label: analytic.
    """

    local: np.ndarray
    advective: np.ndarray
    pressure: np.ndarray
    gravity: np.ndarray
    viscous: np.ndarray
    residual: np.ndarray


def ns_incompressible_terms(u: Callable, p, x, t: float = 0.0, rho: float = 1000.0, mu: float = 1e-3,
                            g=(0.0, 0.0, -G0), h: float = 1e-4, ht: float | None = None, per: str = "mass") -> NSTerms:
    """Term-by-term balance of the incompressible, constant-μ Navier–Stokes equation, Eq. (4.39b).

    Book: §4.6, Eq. (4.39b): ρDu/Dt = −∇p + ρg + μ∇²u; with μ = 0 it is Euler's equation (4.41); with gravity absorbed
    in p′ it is (4.85). The E4 explainer's term bars.

    Parameters
    ----------
    u : velocity field [m/s];  p : pressure field [Pa] (callable or constant);  x : point(s) [m];  t [s];
    rho [kg/m³];  mu [Pa s];  g [m/s²];  h [m] (Laplacian by the 3-point stencil);  ht [s]
    per : "mass" → accelerations [m/s²]; "volume" → forces per volume [N/m³] (×ρ)

    Returns :class:`NSTerms` — arrays (d,) + points.
    Validation: V1 Poiseuille: local = advective = 0 and pressure + viscous = 0; Taylor–Green: local + advective balance
    pressure + viscous; residual ≈ 0 for all :data:`EXACT_SOLUTIONS`; V3 order 2. Label: analytic, converged.
    """
    x_ = _F(x)
    d = x_.shape[0]
    acc = acceleration(u, x_, t, h, ht)
    pres = -grad(as_field(p), x_, t, h) / float(rho)  # −∇p/ρ
    grav = gvec(g, d, _pts(x_))
    visc = (float(mu) / float(rho)) * laplacian(_ufield(u, d), x_, t, h, (d,))  # ν∇²u
    res = acc.local + acc.advective - pres - grav - visc  # Eq. (4.39b) per unit mass
    k = 1.0 if per == "mass" else float(rho) if per == "volume" else None
    if k is None:
        raise ValueError("per must be 'mass' or 'volume'")
    return NSTerms(k * acc.local, k * acc.advective, k * pres, k * grav, k * visc, k * res)


def viscous_force_forms(u: Callable, x, t: float = 0.0, mu: float = 1.0, h: float = 1e-4):
    """The net viscous force per unit volume of an incompressible flow three ways, Eq. (4.40).

    Book: §4.6, Eq. (4.40): (μ∇²u)_j = μ∂²u_j/∂x_i² = 2μ∂S_ij/∂x_i = μ∂/∂x_i(∂u_j/∂x_i + ∂u_i/∂x_j)
    = −με_jik∂ω_k/∂x_i = −μ(∇×ω)_j. Equal only when ∇·u = 0 (for a compressible field 2∂S_ij/∂x_i − ∇²u_j = ∂(∇·u)/∂x_j).

    Parameters: u (2-D or 3-D) [m/s]; x [m]; t [s]; mu [Pa s]; h [m] (the curl of the stencil vorticity is a nested
    stencil).
    Returns (lap, div2S, minus_curl_omega): three (d,) + points arrays [N/m³].
    Validation: V2 sympy: the three forms agree for divergence-free polynomial fields; V1 a compressible field shows the
    extra ∇(∇·u); solid-body rotation gives 0 in all three (the "paradox" resolved). Label: symbolic, analytic.
    """
    x_ = _F(x)
    d = x_.shape[0]
    uf = _ufield(u, d)
    lap = float(mu) * laplacian(uf, x_, t, h, (d,))  # μ∇²u
    twoS = lambda X, T: (lambda G: G + np.swapaxes(G, 0, 1))(grad_vector(uf, X, T, h))  # noqa: E731  2S_ij
    div2S = float(mu) * stress_divergence(twoS, x_, t, h, 0)  # μ ∂(2S_ij)/∂x_i
    om = lambda X, T: vorticity(uf, X, T, h)  # noqa: E731  ω (3 components)
    dom = [ddx(om, x_, t, i, h, (3,)) for i in range(d)]  # ∂ω_k/∂x_i
    if d == 3:
        curl = np.stack([dom[1][2] - dom[2][1], dom[2][0] - dom[0][2], dom[0][1] - dom[1][0]])
    else:  # plane flow: ω = (0, 0, ω₃) → ∇×ω = (∂ω₃/∂y, −∂ω₃/∂x)
        curl = np.stack([dom[1][2], -dom[0][2]])
    return lap, div2S, -float(mu) * curl  # Eq. (4.40)


# ======================================================================================================================
# Exact solutions used as test fields (E4, the C06/C08 figures)
# ======================================================================================================================
EXACT_SOLUTIONS = ("couette", "poiseuille", "pipe_poiseuille", "stokes_first", "taylor_green", "lamb_oseen",
                   "cylinder", "solid_body")

_DEFAULTS = {
    "couette": dict(U=1.0, h=0.01, G=0.0, rho=1000.0, mu=1e-3, p0=0.0),
    "poiseuille": dict(G=100.0, h=0.01, rho=1000.0, mu=1e-3, p0=0.0, g=0.0),
    "pipe_poiseuille": dict(G=100.0, R=0.01, rho=1000.0, mu=1e-3, p0=0.0, g=G0),
    "stokes_first": dict(U=1.0, rho=1000.0, mu=1e-3, p0=0.0),
    "taylor_green": dict(U0=1.0, k=1.0, rho=1000.0, mu=1e-3, p0=0.0),
    "lamb_oseen": dict(Gamma=1.0, rho=1000.0, mu=1e-3, p0=0.0, t0=0.0),
    "cylinder": dict(U=1.0, a=1.0, rho=1000.0, mu=1e-3, p_inf=0.0),
    "solid_body": dict(Omega=1.0, rho=1000.0, mu=1e-3, p0=0.0, g=G0),
}
_ALIASES = {"cylinder_potential": "cylinder"}


def _params(name: str, p: dict) -> tuple[str, dict]:
    name = _ALIASES.get(name, name)
    if name not in _DEFAULTS:
        raise ValueError(f"unknown exact solution {name!r}; choose from {EXACT_SOLUTIONS}")
    q = dict(_DEFAULTS[name])
    if "dpdx" in p:
        p = dict(p)
        p["G"] = -float(p.pop("dpdx"))
    if "dpdz" in p:
        p = dict(p)
        p["G"] = -float(p.pop("dpdz"))
    q.update(p)
    if "nu" in p:  # ν given: μ = ρν
        q["mu"] = float(p["nu"]) * float(q["rho"])
    q["nu"] = float(q["mu"]) / float(q["rho"])
    return name, q


def exact_solution(name: str, x, t: float = 0.0, **p):
    """Velocity and pressure of a standard exact solution of the incompressible Navier–Stokes equations (4.39b).

    Book: §4.6 (these fields satisfy (4.39b); Ch. 8 derives them — standard forms, used here as test fields and
    cross-checked by the residual in the tests, not quoted from the book).

    * "couette" (U, h, G = −dp/dx = 0): walls y = 0 (fixed) and y = h (speed U): u = Uy/h + (G/2μ)y(h − y), p = p₀ − Gx.
    * "poiseuille" (G = −dp/dx, h, mu, rho, g = 0): u = Gy(h − y)/(2μ); p = p₀ − Gx − ρgy (g ≠ 0: gravity along −y, the
      (4.84)–(4.85) demo).  (2-D)
    * "pipe_poiseuille" (G = −dp/dz, R, g): w = G(R² − r²)/(4μ); p = p₀ − (G + ρg)z.  (3-D)
    * "stokes_first" (U, nu): u = U erfc(y/(2√(νt))), p = p₀ (t > 0).  (2-D)
    * "taylor_green" (U0, k, nu, rho): u = U₀ sin kx cos ky e^{−2νk²t}, v = −U₀ cos kx sin ky e^{−2νk²t},
      p = p₀ + (ρU₀²/4)(cos 2kx + cos 2ky)e^{−4νk²t}.  (2-D)
    * "lamb_oseen" (Gamma, nu, t0): u_θ = (Γ/2πr)(1 − e^{−s}), s = r²/a², a² = 4ν(t + t₀) (the Gaussian vortex with
      σ² = 4ν(t + t₀)); p = p₀ − (ρΓ²/8π²a²)[(1 − e^{−s})²/s + 2(E₁(s) − E₁(2s))] (radial balance dp/dr = ρu_θ²/r
      integrated from ∞ in closed form).  (2-D)
    * "cylinder" (U, a, rho, p_inf): ideal flow past a cylinder (ch03 (3.2) field), p = p∞ + ½ρ(U² − |u|²) — an Euler
      solution; its viscous term is exactly 0 ((4.80)).  (2-D)
    * "solid_body" (Omega, rho, g): u = Ω × x (Ω along z), p = p₀ + ρΩ²R²/2 (− ρgz in 3-D).  (2-D or 3-D)

    Parameters
    ----------
    name : one of :data:`EXACT_SOLUTIONS`;  x : point(s) [m];  t : time [s]
    **p : parameters (SI) overriding the defaults (rho [kg/m³], mu [Pa s] or nu [m²/s], U [m/s], h, R, a [m],
        G [Pa/m] (or dpdx), …)

    Returns (u, p) : velocity (shape of x) [m/s] and pressure (points) [Pa].
    Validation: V1 every field gives |``ns_incompressible_terms`` residual| ≈ 0 and ∇·u ≈ 0; V1 form cross-checks with
    the standard Taylor–Green and Lamb–Oseen solutions. Label: analytic.
    """
    name, q = _params(name, p)
    x_ = _F(x)
    rho, mu, nu = float(q["rho"]), float(q["mu"]), float(q["nu"])
    zeros = lambda a: [np.zeros_like(a)] * (x_.shape[0] - 2)  # noqa: E731
    if name == "couette":
        X, Y = x_[0], x_[1]
        U, H, G = float(q["U"]), float(q["h"]), float(q["G"])
        uu = U * Y / H + G / (2.0 * mu) * Y * (H - Y)
        vel = np.stack([uu, np.zeros_like(uu)] + zeros(uu))
        pres = float(q["p0"]) - G * X
    elif name == "poiseuille":
        X, Y = x_[0], x_[1]
        G, H = float(q["G"]), float(q["h"])
        uu = G * Y * (H - Y) / (2.0 * mu)
        vel = np.stack([uu, np.zeros_like(uu)] + zeros(uu))
        pres = float(q["p0"]) - G * X - rho * float(q["g"]) * Y
    elif name == "pipe_poiseuille":
        X, Y, Z = x_[0], x_[1], x_[2]
        G = float(q["G"])
        w = G / (4.0 * mu) * (float(q["R"]) ** 2 - X ** 2 - Y ** 2)
        vel = np.stack([np.zeros_like(w), np.zeros_like(w), w])
        pres = float(q["p0"]) - (G + rho * float(q["g"])) * Z
    elif name == "stokes_first":
        if t <= 0:
            raise ValueError("Stokes' first problem needs t > 0")
        Y = x_[1]
        uu = float(q["U"]) * erfc(Y / (2.0 * np.sqrt(nu * t)))
        vel = np.stack([uu, np.zeros_like(uu)] + zeros(uu))
        pres = float(q["p0"]) + 0.0 * uu
    elif name == "taylor_green":
        X, Y = x_[0], x_[1]
        k, U0 = float(q["k"]), float(q["U0"])
        F = np.exp(-2.0 * nu * k ** 2 * t)
        vel = np.stack([U0 * np.sin(k * X) * np.cos(k * Y) * F, -U0 * np.cos(k * X) * np.sin(k * Y) * F])
        pres = float(q["p0"]) + rho * U0 ** 2 / 4.0 * (np.cos(2 * k * X) + np.cos(2 * k * Y)) * F ** 2
    elif name == "lamb_oseen":
        tt = t + float(q["t0"])
        if tt <= 0:
            raise ValueError("Lamb–Oseen needs t + t0 > 0")
        X, Y = x_[0], x_[1]
        r2 = X ** 2 + Y ** 2
        a2 = 4.0 * nu * tt
        Gam = float(q["Gamma"])
        s = r2 / a2
        with np.errstate(divide="ignore", invalid="ignore"):
            ut_over_r = np.where(r2 > 0, Gam / (2.0 * np.pi) * (-np.expm1(-s)) / np.where(r2 > 0, r2, 1.0),
                                 Gam / (2.0 * np.pi * a2))
            sp_ = np.where(s > 0, s, 1.0)
            bracket = np.where(s > 0, np.expm1(-sp_) ** 2 / sp_ + 2.0 * (exp1(sp_) - exp1(2.0 * sp_)),
                               2.0 * np.log(2.0))
        vel = np.stack([-ut_over_r * Y, ut_over_r * X])
        pres = float(q["p0"]) - rho * Gam ** 2 / (8.0 * np.pi ** 2 * a2) * bracket
    elif name == "cylinder":
        X, Y = x_[0], x_[1]
        U, a = float(q["U"]), float(q["a"])
        r2 = X ** 2 + Y ** 2
        uu = U * (1.0 - a ** 2 * (X ** 2 - Y ** 2) / r2 ** 2)
        vv = -U * a ** 2 * 2.0 * X * Y / r2 ** 2
        vel = np.stack([uu, vv])
        pres = float(q["p_inf"]) + 0.5 * rho * (U ** 2 - uu ** 2 - vv ** 2)
    else:  # solid_body
        Om = float(q["Omega"])
        X, Y = x_[0], x_[1]
        vel = np.stack([-Om * Y, Om * X] + zeros(X))
        pres = float(q["p0"]) + 0.5 * rho * Om ** 2 * (X ** 2 + Y ** 2)
        if x_.shape[0] == 3:
            pres = pres - rho * float(q["g"]) * x_[2]  # hydrostatic part balancing g = (0, 0, −g)
    return vel, _S(pres)


def exact_field(name: str, **p):
    """Callables (u_fn(x, t), p_fn(x, t)) of :func:`exact_solution` for the residual functions. Book: §4.6.
    Label: analytic."""
    nm, q = _params(name, p)
    return (lambda X, T: exact_solution(nm, X, T, **q)[0], lambda X, T: exact_solution(nm, X, T, **q)[1])


def exact_solution_fields(name: str, **p):
    """(u_fn, p_fn, params) — as :func:`exact_field` plus the resolved parameters (rho, mu, nu, …) so a caller passes
    the same values to the residual functions. Book: §4.6. Label: analytic."""
    nm, q = _params(name, p)
    return (lambda X, T: exact_solution(nm, X, T, **q)[0], lambda X, T: exact_solution(nm, X, T, **q)[1], q)


def _scales(name: str, q: dict, t: float) -> tuple[float, float]:
    """(length, time) scales of an exact solution for the default stencil steps."""
    nu = q["nu"]
    L = {"couette": q.get("h", 1.0), "poiseuille": q.get("h", 1.0), "pipe_poiseuille": q.get("R", 1.0),
         "stokes_first": np.sqrt(nu * max(t, 1e-300)), "taylor_green": 1.0 / q.get("k", 1.0),
         "lamb_oseen": np.sqrt(4 * nu * max(t + q.get("t0", 0.0), 1e-300)), "cylinder": q.get("a", 1.0),
         "solid_body": 1.0}[name]
    T = {"stokes_first": max(t, 1e-300), "taylor_green": 1.0 / (nu * q.get("k", 1.0) ** 2),
         "lamb_oseen": max(t + q.get("t0", 0.0), 1e-300)}.get(name, 1.0)
    return float(L), float(T)


def ns_terms_preset(name: str, x: float, y: float, t: float = 0.0, component: int = 0, per: str = "mass",
                    form: str = "laplacian", step: float | None = None, **p) -> dict:
    """The five terms of (4.39b) and the residual for one component of an exact solution at (x, y) — plain floats
    (the parity-friendly wrapper of :func:`ns_incompressible_terms` on :func:`exact_solution`; E4).

    Book: §4.6, Eq. (4.39b). Parameters: name (2-D solutions of :data:`EXACT_SOLUTIONS`); x, y [m]; t [s]; component
    0 (x) or 1 (y); per "mass" [m/s²] or "volume" [N/m³]; form : how the viscous term is computed, "laplacian" (μ∇²u),
    "div2S" (2μ∂S_ij/∂x_i) or "curl" (−μ∇×ω) — the three forms of (4.40); step : stencil step [m] (default 2e-4 × the solution's length
    scale; the time step is 1e-4 × its time scale); **p : the solution's parameters (note: ``h`` is the channel gap).
    Returns dict(local, advective, pressure, gravity, viscous, residual).
    Validation: V1 ``ns_terms_preset("poiseuille", 0.0, 2.5e-4, component=0, per="volume", G=100.0, h=1e-3, mu=1e-3)``
    → pressure +100, viscous −100, others 0 (N/m³). Label: analytic.
    """
    nm, q = _params(name, p)
    L, T = _scales(nm, q, t)
    hs = 2e-4 * L if step is None else float(step)
    u_fn, p_fn = exact_field(nm, **q)
    g = (0.0, -float(q["g"])) if nm == "poiseuille" else (0.0, 0.0)
    X = np.array([float(x), float(y)])
    terms = ns_incompressible_terms(u_fn, p_fn, X, t, q["rho"], q["mu"], g, hs, 1e-4 * T, per)
    out = {k: float(v[component]) for k, v in terms._asdict().items()}
    if form != "laplacian":
        forms = dict(zip(("laplacian", "div2S", "curl"), viscous_force_forms(u_fn, X, t, q["mu"], hs)))
        if form not in forms:
            raise ValueError("form must be 'laplacian', 'div2S' or 'curl'")
        k = 1.0 / q["rho"] if per == "mass" else 1.0
        v = float(forms[form][component]) * k
        out["residual"] += out["viscous"] - v
        out["viscous"] = v
    return out


# ======================================================================================================================
# §4.9 Lamb's identity (4.68)
# ======================================================================================================================
def lamb_vector(u: Callable, x, t: float = 0.0, h: float = 1e-4) -> np.ndarray:
    """Lamb vector u × ω [m/s²] (3 components; a plane field with x of shape (2, …) is the plane z = 0, ω = (0, 0, ω₃)),
    the right side of (4.69)–(4.70).

    Book: §4.9, Eqs. (4.68)–(4.70) (ω = ∇×u, (3.16)). Validation: V1 solid-body rotation u = Ω × x = Ω(−y, x, 0),
    ω = (0, 0, 2Ω): u × ω = +2Ω²(x, y, 0) (outward — balanced by the pressure gradient in ∇B = u × ω); irrotational
    field → 0. Label: analytic.
    """
    x_ = _F(x)
    d = x_.shape[0]
    U3 = pad3(ev(u, x_, t, (d,)))
    om = vorticity(_ufield(u, d), x_, t, h)
    return np.cross(U3, om, axis=0)  # u × ω


def lamb_identity_terms(u: Callable, x, t: float = 0.0, h: float = 1e-4) -> dict:
    """Both sides of Lamb's identity u_i∂u_j/∂x_i = −(u × ω)_j + ∂(½u_i²)/∂x_j, Eq. (4.68).

    Book: §4.9, Eq. (4.68) (Exercise 4.50; our D24 derives it with ε–δ). Returns dict(advective, minus_u_cross_omega,
    grad_ke, residual = advective − minus_u_cross_omega − grad_ke), each (d,) + points [m/s²].
    Validation: V2 :func:`lamb_identity_sym` ≡ 0 for a generic field; V3 order 2; the wrong-sign variant fails.
    Label: symbolic, converged.
    """
    x_ = _F(x)
    d = x_.shape[0]
    adv = acceleration(u, x_, t, h).advective
    ml = -lamb_vector(u, x_, t, h)[:d]
    ke = grad(lambda X, T: 0.5 * np.sum(ev(u, X, T, (d,)) ** 2, axis=0), x_, t, h)
    return {"advective": adv, "minus_u_cross_omega": ml, "grad_ke": ke, "residual": adv - ml - ke}  # Eq. (4.68)


def lamb_identity_sym(u_exprs, coords) -> list:
    """Symbolic residual of (4.68) in 3-D: (u·∇)u − [−u × (∇×u) + ∇(½|u|²)], simplified (→ [0, 0, 0]).

    Book: §4.9, Eq. (4.68). Label: symbolic.
    """
    u = sp.Matrix([sp.sympify(e) for e in u_exprs])
    x1, x2, x3 = coords
    om = sp.Matrix([sp.diff(u[2], x2) - sp.diff(u[1], x3), sp.diff(u[0], x3) - sp.diff(u[2], x1),
                    sp.diff(u[1], x1) - sp.diff(u[0], x2)])
    ke = sp.Rational(1, 2) * (u.T * u)[0]
    adv = [sum(u[i] * sp.diff(u[j], coords[i]) for i in range(3)) for j in range(3)]
    lam = u.cross(om)
    return [sp.simplify(sp.expand(adv[j] - (-lam[j] + sp.diff(ke, coords[j])))) for j in range(3)]  # (4.68)


# ======================================================================================================================
# §4.8 energy
# ======================================================================================================================
def _sigma_or_zero(sigma, d):
    return _zero_tensor(d) if sigma is None else sigma


def stress_work_split(p, sigma: Callable | None, u: Callable, x, t: float = 0.0, h: float = 1e-4) -> dict:
    """The rate of work of the surface stresses split four ways, Eq. (4.54).

    Book: §4.8, Eq. (4.54): ∂(τ_iju_j)/∂x_i = τ_ij∂u_j/∂x_i + u_j∂τ_ij/∂x_i
    = (−p∂u_j/∂x_j + σ_ij∂u_j/∂x_i) + (−u_j∂p/∂x_j + u_j∂σ_ij/∂x_i): deformation work (changes internal energy) +
    force work (changes kinetic energy), with τ = −pδ + σ (4.27).

    Parameters: p [Pa] (callable or constant); sigma : viscous stress field [Pa] (None → 0; Newtonian:
    :func:`newtonian_viscous_stress_field`); u [m/s]; x [m]; t [s]; h [m].
    Returns dict(total, deformation_pressure, deformation_viscous, force_pressure, force_viscous, residual) [W/m³].
    Validation: V1 the four parts sum to the total on smooth fields (residual O(h²)). Label: analytic, converged.
    """
    x_ = _F(x)
    d = x_.shape[0]
    uf = _ufield(u, d)
    pf = as_field(p)
    sig = _sigma_or_zero(sigma, d)
    eye = lambda X: np.eye(d).reshape((d, d) + (1,) * (_F(X).ndim - 1))  # noqa: E731
    tau = lambda X, T: -ev(pf, X, T) * eye(X) + ev(sig, X, T, (d, d))  # noqa: E731  Eq. (4.27)
    total = div(lambda X, T: np.einsum("ij...,j...->i...", tau(X, T), uf(X, T)), x_, t, h)  # ∂(τ_ij u_j)/∂x_i
    G = grad_vector(uf, x_, t, h)  # G[j, i] = ∂u_j/∂x_i
    U = uf(x_, t)
    pdef = -ev(pf, x_, t) * np.trace(G, axis1=0, axis2=1)
    vdef = np.einsum("ij...,ji...->...", ev(sig, x_, t, (d, d)), G)  # σ_ij ∂u_j/∂x_i
    pfor = -np.sum(U * grad(pf, x_, t, h), axis=0)
    vfor = np.sum(U * stress_divergence(sig, x_, t, h, 0), axis=0)  # u_j ∂σ_ij/∂x_i
    res = total - (pdef + vdef + pfor + vfor)  # Eq. (4.54)
    return {k: _S(v) for k, v in dict(total=total, deformation_pressure=pdef, deformation_viscous=vdef,
                                      force_pressure=pfor, force_viscous=vfor, residual=res).items()}


def kinetic_energy_budget(rho, u: Callable, p, sigma: Callable | None, g=(0.0, 0.0, -G0), x=None, t: float = 0.0,
                          h: float = 1e-4, ht: float | None = None) -> dict:
    """Mechanical-energy equation ρD(½u_j²)/Dt = ρg_ju_j − u_j∂p/∂x_j + u_j∂σ_ij/∂x_i, Eq. (4.56).

    Book: §4.8, Eq. (4.56) (Cauchy's (4.24) dotted with u_j; the book says "multiplying (4.22)", from which one must also
    subtract u_j × continuity — our D21 uses (4.24)). Holds only for fields that satisfy the momentum equation.

    Parameters: rho [kg/m³]; u [m/s]; p [Pa]; sigma : viscous stress field [Pa] (None → 0); g [m/s²]; x [m]; t [s];
    h [m]; ht [s].
    Returns dict(lhs = ρD(½|u|²)/Dt, gravity_work, pressure_work, viscous_force_work, residual) [W/m³].
    Validation: V1 exact NS solutions (Taylor–Green: lhs = viscous work < 0; Couette: all 0) give residual ≈ 0.
    Label: analytic, converged.
    """
    x_ = _F(x)
    d = x_.shape[0]
    uf = _ufield(u, d)
    r = ev(as_field(rho), x_, t)
    ke = lambda X, T: 0.5 * np.sum(uf(X, T) ** 2, axis=0)  # noqa: E731
    lhs = r * material_derivative(ke, uf, x_, t, h, ht)
    U = uf(x_, t)
    gw = r * np.sum(U * gvec(g, d, _pts(x_)), axis=0)
    pw = -np.sum(U * grad(as_field(p), x_, t, h), axis=0)
    vw = np.sum(U * stress_divergence(_sigma_or_zero(sigma, d), x_, t, h, 0), axis=0)
    return {k: _S(v) for k, v in dict(lhs=lhs, gravity_work=gw, pressure_work=pw, viscous_force_work=vw,
                                      residual=lhs - gw - pw - vw).items()}  # Eq. (4.56)


def internal_energy_terms(rho, u: Callable, p, sigma: Callable | None, q: Callable | None, x, t: float = 0.0,
                          h: float = 1e-4, ht: float | None = None, e: Callable | None = None) -> dict:
    """The first law for a fluid particle, De/Dt = −p Dv/Dt + (1/ρ)σ_ijS_ij − (1/ρ)∂q_i/∂x_i, Eq. (4.57), per unit mass.

    Book: §4.8, Eq. (4.57) with v = 1/ρ and ε = σ_ijS_ij/ρ (4.58).

    Parameters: rho [kg/m³]; u [m/s]; p [Pa]; sigma : viscous stress field [Pa] (None → 0); q : heat-flux field [W/m²]
    (None → 0); x [m]; t [s]; h [m]; ht [s]; e : internal-energy field [J/kg] (optional, for the residual).
    Returns dict(De_Dt (from e if given, else the sum of the right side), pressure_work, dissipation, conduction,
    residual (NaN without e)) [W/kg].
    Validation: V1 plane Couette with viscous heating: dissipation = νγ², conduction = −dissipation at steady state.
    Label: analytic.
    """
    x_ = _F(x)
    d = x_.shape[0]
    uf = _ufield(u, d)
    rf = as_field(rho)
    r = ev(rf, x_, t)
    pw = -ev(as_field(p), x_, t) * material_derivative(lambda X, T: 1.0 / ev(rf, X, T), uf, x_, t, h, ht)
    G = grad_vector(uf, x_, t, h)
    S = 0.5 * (G + np.swapaxes(G, 0, 1))
    sig = ev(_sigma_or_zero(sigma, d), x_, t, (d, d))
    eps = np.sum(sig * S, axis=(0, 1)) / r  # ε = σ_ij S_ij / ρ (4.58)
    cond = 0.0 * r if q is None else -div(q, x_, t, h) / r
    rhs = pw + eps + cond
    if e is None:
        DeDt, res = rhs, np.full(np.shape(rhs), np.nan)
    else:
        DeDt = material_derivative(lambda X, T: ev(e, X, T), uf, x_, t, h, ht)
        res = DeDt - rhs  # Eq. (4.57)
    return {k: _S(v) for k, v in dict(De_Dt=DeDt, pressure_work=pw, dissipation=eps, conduction=cond,
                                      residual=res).items()}


def internal_energy_residual(rho, u: Callable, e: Callable, p, T: Callable, mu=1e-3, mu_v=0.0, k=0.6, x=None,
                             t: float = 0.0, h: float = 1e-4, ht: float | None = None):
    """Residual of the internal-energy equation with Fourier conduction, Eq. (4.60):
    ρDe/Dt − [−p∂u_m/∂x_m + 2μ(S_ij − ⅓S_mmδ_ij)² + μ_v(∂u_m/∂x_m)² + ∂/∂x_i(k∂T/∂x_i)].

    Book: §4.8, Eq. (4.60) ((4.57) with (4.58) and q = −k∇T (1.2)).

    Parameters: rho [kg/m³]; u [m/s]; e : internal energy field [J/kg]; p [Pa]; T : temperature field [K]; mu, mu_v
    [Pa s]; k [W/(m K)] (constant or callable); x [m]; t [s]; h [m]; ht [s].
    Returns residual [W/m³].
    Validation: V1 steady plane Couette with viscous heating T = T₀ + (μU²/2k)(y/h)(1 − y/h), e = c_vT: residual 0.
    Label: analytic, converged.
    """
    from .constitutive import dissipation_rate
    x_ = _F(x)
    d = x_.shape[0]
    uf = _ufield(u, d)
    r = ev(as_field(rho), x_, t)
    lhs = r * material_derivative(lambda X, T_: ev(e, X, T_), uf, x_, t, h, ht)
    G = grad_vector(uf, x_, t, h)
    divu = np.trace(G, axis1=0, axis2=1)
    phi = r * dissipation_rate(G, r, mu, mu_v)  # 2μ(S − ⅓S_mmδ)² + μ_v S_mm²  (= ρε)
    kf = as_field(k)
    cond = div(lambda X, T_: ev(kf, X, T_) * grad(T, X, T_, h), x_, t, h)  # ∇·(k∇T)
    return _S(lhs - (-ev(as_field(p), x_, t) * divu + phi + cond))  # Eq. (4.60)


def _generic_fields(coords, t):
    """Undefined sympy functions ρ, e, p, T, u_i, q_i of (coords, t) for the identity checks."""
    args = tuple(coords) + (t,)
    n = len(coords)
    rho, e, p, T = (sp.Function(s)(*args) for s in ("rho", "e", "p", "T"))
    u = [sp.Function(f"u{i + 1}")(*args) for i in range(n)]
    q = [sp.Function(f"q{i + 1}")(*args) for i in range(n)]
    return rho, e, p, T, u, q


def _sym_parts(rho, e, p, u, q, coords, t, mu, mu_v, g):
    n = len(coords)
    D = lambda F: sp.diff(F, t) + sum(u[i] * sp.diff(F, coords[i]) for i in range(n))  # noqa: E731  (3.5)
    divu = sum(sp.diff(u[m], coords[m]) for m in range(n))
    sig = [[mu * (sp.diff(u[i], coords[j]) + sp.diff(u[j], coords[i]))
            + (mu_v - sp.Rational(2, 3) * mu) * divu * (1 if i == j else 0) for j in range(n)] for i in range(n)]
    tau = [[-p * (1 if i == j else 0) + sig[i][j] for j in range(n)] for i in range(n)]
    ke = sp.Rational(1, 2) * sum(uj ** 2 for uj in u)
    gu = sum(g[j] * u[j] for j in range(n))
    divq = sum(sp.diff(q[i], coords[i]) for i in range(n))
    return D, divu, sig, tau, ke, gu, divq


def total_energy_residual_sym(rho=None, u=None, e=None, p=None, q=None, coords=None, t=None, mu=None, mu_v=None,
                              g=None):
    """Symbolic residual of the conservative total-energy equation, Eq. (4.53):
    ∂(ρ[e + ½u_j²])/∂t + ∂(ρ[e + ½u_j²]u_i)/∂x_i − ρg_iu_i − ∂(τ_iju_j)/∂x_i + ∂q_i/∂x_i, with Newtonian τ (4.37).

    Book: §4.8, Eq. (4.53). Missing arguments become generic sympy functions/symbols (coords x, y, z; t).
    Returns a sympy expression. Label: symbolic.
    """
    coords = coords or sp.symbols("x y z", real=True)
    t = t if t is not None else sp.Symbol("t", real=True)
    G_rho, G_e, G_p, _, G_u, G_q = _generic_fields(coords, t)
    rho = G_rho if rho is None else rho
    e = G_e if e is None else e
    p = G_p if p is None else p
    u = G_u if u is None else list(u)
    q = G_q if q is None else list(q)
    mu = sp.Symbol("mu") if mu is None else mu
    mu_v = sp.Symbol("mu_v") if mu_v is None else mu_v
    n = len(coords)
    g = [sp.Symbol(f"g{i + 1}") for i in range(n)] if g is None else list(g)
    D, divu, sig, tau, ke, gu, divq = _sym_parts(rho, e, p, u, q, coords, t, mu, mu_v, g)
    E = rho * (e + ke)
    return (sp.diff(E, t) + sum(sp.diff(E * u[i], coords[i]) for i in range(n)) - rho * gu
            - sum(sp.diff(sum(tau[i][j] * u[j] for j in range(n)), coords[i]) for i in range(n)) + divq)  # (4.53)


def energy_identities_sym(coords=None, t=None) -> dict:
    """The chain of energy identities of §4.8 checked with sympy for generic fields ρ, e, p, T, u_i, q_i.

    Book: §4.8, Eqs. (4.53)–(4.57), (4.60), (4.112). With C = ∂ρ/∂t + ∂(ρu_i)/∂x_i (4.7), Cauchy_j = ρDu_j/Dt − ρg_j −
    ∂τ_ij/∂x_i (4.24) and R_N = (left − right) of equation N, each entry expands to exactly 0:

    * "4.53_to_4.55": R53 − R55 − (e + ½u²)C  (expanding (4.53) with continuity gives (4.55); D20)
    * "4.24_dot_u_to_4.56": u_j·Cauchy_j − R56  ((4.56) is u·(4.24); D21)
    * "4.55_minus_4.56_to_4.57": R55 − R56 − ρR57 − (p/ρ)C  ((4.57) with v = 1/ρ and (4.8); D22)
    * "4.54_split": ∂(τ_iju_j)/∂x_i − [(−p∂u_j/∂x_j + σ_ij∂u_j/∂x_i) + (−u_j∂p/∂x_j + u_j∂σ_ij/∂x_i)]  ((4.54))
    * "4.60_to_4.112": R112 − R60 + (p/ρ)C  (enthalpy form with h = e + p/ρ, μ_v = 0, q = −k∇T; analysis a-D58)

    Returns dict name → expanded expression (all 0). ≈ 1–3 s in 2-D coordinates (default x, y). Label: symbolic.
    """
    coords = coords or sp.symbols("x y", real=True)
    t = t if t is not None else sp.Symbol("t", real=True)
    n = len(coords)
    rho, e, p, T, u, q = _generic_fields(coords, t)
    mu, mu_v, kk = sp.symbols("mu mu_v k")
    g = [sp.Symbol(f"g{i + 1}") for i in range(n)]
    D, divu, sig, tau, ke, gu, divq = _sym_parts(rho, e, p, u, q, coords, t, mu, mu_v, g)
    C = continuity_residual_sym(rho, u, coords, t)
    dudx = lambda j, i: sp.diff(u[j], coords[i])  # noqa: E731
    R53 = total_energy_residual_sym(rho, u, e, p, q, coords, t, mu, mu_v, g)
    pdef = -p * divu
    vdef = sum(sig[i][j] * dudx(j, i) for i in range(n) for j in range(n))
    pfor = -sum(u[j] * sp.diff(p, coords[j]) for j in range(n))
    vfor = sum(u[j] * sp.diff(sig[i][j], coords[i]) for i in range(n) for j in range(n))
    R55 = rho * D(e + ke) - (rho * gu + (pdef + vdef) + (pfor + vfor) - divq)  # Eq. (4.55)
    R56 = rho * D(ke) - (rho * gu + pfor + vfor)  # Eq. (4.56)
    R57 = D(e) - (-p * D(1 / rho) + vdef / rho - divq / rho)  # Eq. (4.57)
    cauchy = [rho * D(u[j]) - rho * g[j] - sum(sp.diff(tau[i][j], coords[i]) for i in range(n)) for j in range(n)]
    split = sum(sp.diff(sum(tau[i][j] * u[j] for j in range(n)), coords[i]) for i in range(n)) - (
        pdef + vdef + pfor + vfor)  # Eq. (4.54)
    sub0 = {mu_v: 0}
    cond = sum(sp.diff(kk * sp.diff(T, coords[i]), coords[i]) for i in range(n))
    vdef0 = vdef.subs(sub0)
    R60 = rho * D(e) - (pdef + vdef0 + cond)  # Eq. (4.60) (the squares of (4.58) written as σ_ijS_ij)
    R112 = rho * D(e + p / rho) - (D(p) + vdef0 + cond)  # Eq. (4.112)
    out = {"4.53_to_4.55": R53 - R55 - (e + ke) * C,
           "4.24_dot_u_to_4.56": sum(u[j] * cauchy[j] for j in range(n)) - R56,
           "4.55_minus_4.56_to_4.57": R55 - R56 - rho * R57 - (p / rho) * C,
           "4.54_split": split,
           "4.60_to_4.112": R112 - R60 + (p / rho) * C}
    res = {}
    for k, v in out.items():
        v = sp.expand(v)
        res[k] = v if v == 0 else sp.simplify(v)
    return res


def energy_forms_sym(coords=None, t=None):
    """(4.60) ⇔ (4.112): R112 − R60 + (p/ρ)C for generic fields — the sympy expression 0 (see
    :func:`energy_identities_sym` for the whole chain). Book: §4.8, §4.11, Eqs. (4.60), (4.112). Label: symbolic."""
    return energy_identities_sym(coords, t)["4.60_to_4.112"]


def entropy_terms(rho, T: Callable, q: Callable, eps, x, h: float = 1e-4, t: float = 0.0) -> dict:
    """Entropy of a fluid particle, Ds/Dt = −(1/ρT)∂q_i/∂x_i + ε/T = −(1/ρ)∂(q_i/T)/∂x_i − (q_i/ρT²)∂T/∂x_i + ε/T,
    Eq. (4.62).

    Book: §4.8 (end), Eq. (4.62) ((4.57) with Gibbs (4.61)).

    Parameters: rho [kg/m³]; T [K] field; q : heat flux field [W/m²]; eps : dissipation [W/kg] (constant or field);
    x [m]; h [m]; t [s].
    Returns dict(flux_divergence = −∇·(q/T)/ρ, conduction_production = −q·∇T/(ρT²), dissipation_production = ε/T,
    total, split_residual = [−∇·q/(ρT)] − flux_divergence − conduction_production) [W/(kg K)].
    Validation: V2 the middle and right forms agree (split_residual ≈ 0); with q = −k∇T the production is
    k|∇T|²/(ρT²) ≥ 0 (:func:`entropy_production`). Label: analytic, converged.
    """
    x_ = _F(x)
    d = x_.shape[0]
    r = ev(as_field(rho), x_, t)
    TT = ev(T, x_, t)
    Q = ev(q, x_, t, (d,))
    ep = ev(as_field(eps), x_, t)
    cond = -div(q, x_, t, h) / (r * TT)
    fd = -div(lambda X, T_: ev(q, X, T_, (d,)) / ev(T, X, T_), x_, t, h) / r
    cp_ = -np.sum(Q * grad(T, x_, t, h), axis=0) / (r * TT ** 2)
    visc = ep / TT
    return {k: _S(v) for k, v in dict(flux_divergence=fd, conduction_production=cp_, dissipation_production=visc,
                                      total=cond + visc, split_residual=cond - (fd + cp_)).items()}  # Eq. (4.62)


def entropy_production(gradT, T, k, rho, eps=0.0):
    """Entropy production per unit mass by conduction and viscosity, k|∇T|²/(ρT²) + ε/T, the last two terms of (4.63).

    Book: §4.8, Eq. (4.63): Ds/Dt = (1/ρ)∂(k∂T/∂x_i/T)/∂x_i + (k/ρT²)(∂T/∂x_i)² + ε/T. The second law requires the
    production to be ≥ 0, satisfied by μ, μ_v, k ≥ 0 (the book writes "μ, κ, k > 0": κ there is the bulk viscosity μ_v,
    a 4th-edition leftover, and ≥ 0 is what is required — analysis §9 item 4).

    Parameters: gradT : temperature gradient [K/m] (vector (d,)/(d, N) or a scalar magnitude); T [K]; k [W/(m K)];
    rho [kg/m³]; eps : dissipation rate [W/kg].
    Returns [W/(kg K)]; negative only for unphysical k < 0 or ε < 0 (E7's "second law violated" preset).
    Validation: V1 ≥ 0 over random inputs with k, ε ≥ 0; 0 for uniform T and ε = 0. Label: analytic.
    """
    g_ = _F(gradT)
    g2 = g_ ** 2 if g_.ndim == 0 else np.sum(g_ ** 2, axis=0)
    return _S(_F(k) * g2 / (_F(rho) * _F(T) ** 2) + _F(eps) / _F(T))  # Eq. (4.63), production terms


# ======================================================================================================================
# §4.9 neglect of gravity; Boussinesq
# ======================================================================================================================
def perturbation_fields(p, rho, z, rho_s, g: float = G0, p_s0: float | None = None):
    """Departures from the hydrostatic base state p′ = p − p_s, ρ′ = ρ − ρ_s along a vertical column (before (4.84)).

    Book: §4.9 "Neglect of gravity in constant density flows": 0 = −∇p_s + ρ_sg (≡ (1.8)); subtracting it from (4.39b)
    gives (4.84) with p′, ρ′ (here p′ is a perturbation — not the dummy variable of (4.67)).

    Parameters: p, rho : pressure [Pa] and density [kg/m³] at heights z (arrays); z [m] (up); rho_s : base density
    (callable of z or array on z) [kg/m³]; g [m/s²]; p_s0 : base pressure at z[0] (default p[0]).
    Returns (p_pert, rho_pert) [Pa], [kg/m³]; p_s from ``core.statics.integrate_hydrostatic``.
    Validation: V1 a hydrostatic column gives p′ = ρ′ = 0 (1e-10); constant ρ: p′ = p − (p₀ − ρgz). Label: analytic.
    """
    z_ = _F(z)
    rs = rho_s if callable(rho_s) else (lambda zz, _r=_F(rho_s): np.interp(zz, z_, _r))
    p0 = float(_F(p).ravel()[0]) if p_s0 is None else float(p_s0)
    ps = integrate_hydrostatic(z_, lambda zz, pp: rs(zz), p0, g=g, z0=float(z_[0]))
    return _F(p) - _F(ps), _F(rho) - _F(rs(z_))


def buoyancy(rho_pert, rho0, g: float = G0):
    """Buoyancy (vertical force per unit mass) b = −gρ′/ρ₀ [m/s²], the vertical component of ρ′g/ρ₀ in (4.86).

    Book: §4.9, Eq. (4.86) with g = −g e_z. Positive (upward) for lighter fluid (ρ′ < 0).
    Validation: V1 N² = ∂b_s/∂z = −(g/ρ₀)dρ_s/dz equals ``core.stratification.brunt_vaisala_sq`` (incompressible).
    Label: analytic.
    """
    return _S(-float(g) * _F(rho_pert) / _F(rho0))


def boussinesq_momentum_terms(u: Callable, p_pert, rho_pert, x, t: float = 0.0, rho0: float = 1000.0,
                              nu: float = 1e-6, g: float = G0, h: float = 1e-4, ht: float | None = None) -> dict:
    """Term split of the Boussinesq momentum equation Du/Dt = −∇p′/ρ₀ + (ρ′/ρ₀)g + ν∇²u, Eq. (4.86).

    Book: §4.9, Eq. (4.86) ((4.84) divided by ρ_s with ρ/ρ_s ≈ 1, μ/ρ_s ≈ ν); gravity along −(last coordinate).
    Parameters: u [m/s]; p_pert [Pa] and rho_pert [kg/m³] (callables or constants); x [m]; t [s]; rho0 [kg/m³];
    nu [m²/s]; g [m/s²] (magnitude); h [m]; ht [s].
    Returns dict(local, advective, pressure, buoyancy, viscous, residual, dropped = (ρ′/ρ₀)Du/Dt — the inertia term the
    approximation neglects) [m/s²].
    Validation: V1 a fluid at rest in any base state with p′ = ρ′ = 0 gives all terms 0. Label: analytic.
    """
    x_ = _F(x)
    d = x_.shape[0]
    acc = acceleration(u, x_, t, h, ht)
    pres = -grad(as_field(p_pert), x_, t, h) / float(rho0)
    rp = ev(as_field(rho_pert), x_, t)
    gv = np.zeros((d,) + _pts(x_))
    gv[-1] = -float(g)
    buoy = (rp / float(rho0)) * gv  # (ρ′/ρ₀) g
    visc = float(nu) * laplacian(_ufield(u, d), x_, t, h, (d,))
    res = acc.local + acc.advective - pres - buoy - visc  # Eq. (4.86)
    return {"local": acc.local, "advective": acc.advective, "pressure": pres, "buoyancy": buoy, "viscous": visc,
            "residual": res, "dropped": (rp / float(rho0)) * acc.a}


def heat_equation_terms(T: Callable, u: Callable, x, t: float = 0.0, rho: float = 1000.0, cp: float = 4182.0,
                        k: float = 0.6, eps=0.0, h: float = 1e-4, ht: float | None = None) -> dict:
    """Terms of the Boussinesq heat equation ρC_pDT/Dt = ρε − ∇·q, q = −k∇T (constant k), Eq. (4.88), divided by ρC_p.

    Book: §4.9, Eq. (4.88) (C_p, not C_v, because −p∇·u = −pαDT/Dt = −ρ(C_p − C_v)DT/Dt for a perfect gas); dropping ε
    and with κ = k/ρC_p it is (4.89).
    Returns dict(DT_Dt, dissipation = ε/C_p, conduction = κ∇²T, residual) [K/s].
    Validation: V1 the moving Gaussian blob (ε = 0) gives residual ≈ 0. Label: analytic, converged.
    """
    x_ = _F(x)
    DT = material_derivative(lambda X, T_: ev(T, X, T_), u, x_, t, h, ht)
    ds = ev(as_field(eps), x_, t) / float(cp)
    cond = float(k) / (float(rho) * float(cp)) * laplacian(T, x_, t, h)
    return {k_: _S(v) for k_, v in dict(DT_Dt=DT, dissipation=ds, conduction=cond,
                                        residual=DT - ds - cond).items()}  # Eq. (4.88)/(ρC_p)


def temperature_equation_residual(T: Callable, u: Callable, x, t: float = 0.0, kappa: float = 1.4e-7,
                                  h: float = 1e-4, ht: float | None = None):
    """Residual DT/Dt − κ∇²T of the Boussinesq temperature equation, Eq. (4.89) (κ = k/ρC_p [m²/s]).

    Book: §4.9, Eq. (4.89). Units K/s.
    Validation: V1 the advected, spreading Gaussian (``ch04.gaussian_blob_advection_diffusion``) gives 0; V3 order 2.
    Label: analytic, converged.
    """
    x_ = _F(x)
    DT = material_derivative(lambda X, T_: ev(T, X, T_), u, x_, t, h, ht)
    return _S(DT - float(kappa) * laplacian(T, x_, t, h))  # Eq. (4.89)
