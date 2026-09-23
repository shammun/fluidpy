"""Chapter 4 — Conservation laws: the chapter's worked examples, the fields its figures and explainers use, and
re-exports of the reusable primitives the chapter introduced in ``fluidpy.core``.

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e (2012), Ch. 4, §§4.1–4.11, Eqs. (4.1)–(4.119), Examples 4.1–4.8.
Every equation was transcribed from the rendered page images (chapters/pages/ch04/p124–p178).

Where the physics lives (all public names are re-exported here, so ``ch04.<name>`` reaches every callable)
--------------------------------------------------------------------------------------------------------
* ``core.conservation`` (CB) — integral budgets (4.1)–(4.5), (4.13)–(4.17), (4.46)–(4.48), (4.64)–(4.65).
* ``core.streamfunction`` (SF) — (4.11)–(4.12), plane and axisymmetric stream functions, E2 presets.
* ``core.navier_stokes`` (NS) — residuals and term splits of (4.7)–(4.10), (4.22)–(4.24), (4.38)–(4.41),
  (4.53)–(4.63), (4.68), (4.84)–(4.89); exact solutions for tests.
* ``core.constitutive`` (NW) — Newtonian stress (4.25)–(4.37), dissipation (4.58)–(4.59).
* ``core.rotating`` (RF) — noninertial frames (4.42)–(4.45), Coriolis, centrifugal, effective gravity.
* ``core.curvilinear`` (CU) — Appendix-B operators (Example 4.5) — reached as ``ch04.CU``.
* ``core.bernoulli`` (BE) — the Bernoulli equations (4.19), (4.66)–(4.83) and applications.
* ``core.interfaces`` (IF) — boundary conditions (4.90)–(4.98).
* ``core.similarity`` (SIM) — dimensionless forms and numbers (4.99)–(4.119).
* ``core.thermo.helmholtz_free_energy`` — (4.94).
* this module — the examples (4.1 wake, 4.2 stream-tube element, 4.3 bore, 4.4 rocket, 4.5 pump, 4.6 sprinkler, 4.7
  meniscus, 4.8 ship model), the §4.1 closure ledger, the expanding flow of C01/C02, Couette heating (E7), the
  Boussinesq scenarios (E8), the E1 control-volume scenarios and helpers for the figures.

Defaults: g = 9.81 m/s² (the book's value in Ch. 4; ch01's G0 = 9.80665 is standard gravity), SI units throughout.
Book typos handled (analysis §9): (4.15)'s trailing "= 0" (an identity, not zero); (4.51) dA → dV; (4.74) gauge sign
φ → φ − ∫B dt′; (4.63) "μ, κ, k > 0" → μ, μ_v, k ≥ 0; §4.10 curve C sign (+); Example 4.7's dropped minus and
"η = h/δ" → γ = h/δ; Example 4.2's (½)ρU² → ½U²; the Cauchy prose ∂τ_ij/∂x_j → ∂τ_ij/∂x_i; "Section 3.6" → §3.4;
"(4.100)" → (4.101) and "(4.106), (4.107)" → (4.109), (4.112); Example 4.8's total 9.14e5 N is 9.15e5 N before the
book's rounding. Book-quoted numbers live only in the git-ignored ``tests/book_values_ch04.json``.
"""
from __future__ import annotations

from typing import Callable

import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp, tplquad
from scipy.special import erfc

from .core import bernoulli as BE  # noqa: F401
from .core import conservation as CB  # noqa: F401
from .core import constitutive as NW  # noqa: F401
from .core import curvilinear as CU  # noqa: F401
from .core import interfaces as IF  # noqa: F401
from .core import navier_stokes as NS  # noqa: F401
from .core import rotating as RF  # noqa: F401
from .core import similarity as SIM  # noqa: F401
from .core import streamfunction as SF  # noqa: F401
from .core._util import as_scalar_if_0d
from .core.bernoulli import *  # noqa: F401,F403
from .core.conservation import *  # noqa: F401,F403
from .core.constitutive import *  # noqa: F401,F403
from .core.interfaces import *  # noqa: F401,F403
from .core.navier_stokes import *  # noqa: F401,F403
from .core.rotating import *  # noqa: F401,F403
from .core.similarity import *  # noqa: F401,F403
from .core.streamfunction import *  # noqa: F401,F403
from .core.thermo import G_BOOK, GAMMA_AIR, R_AIR, helmholtz_free_energy  # noqa: F401
from .core.transport import GrowingSphere, MovingBox  # noqa: F401

G = G_BOOK  #: default g [m/s²] of this chapter's functions (the book's 9.81)
_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d


# ======================================================================================================================
# §4.1 the closure ledger
# ======================================================================================================================
_CLOSURE = {
    "empty": dict(equation_names=[], unknown_names=[]),
    "cauchy": dict(equation_names=["continuity (4.7)", "Cauchy (4.24) ×3", "thermodynamic equation ×2"],
                   unknown_names=["ρ", "u_j ×3", "τ_ij ×9"], equations=6, unknowns=13),
    "navier_stokes": dict(equation_names=["continuity (4.7)", "Navier–Stokes (4.38) ×3"],
                          unknown_names=["ρ", "p", "u_j ×3"], equations=4, unknowns=5),
    "barotropic": dict(equation_names=["continuity (4.7)", "Navier–Stokes (4.38) ×3", "p = p(ρ) (or ρ = const)"],
                       unknown_names=["ρ", "p", "u_j ×3"], equations=5, unknowns=5),
    "full": dict(equation_names=["continuity (4.7)", "Navier–Stokes (4.38) ×3", "energy (4.60)",
                                 "thermodynamic equation ×2"], unknown_names=["ρ", "e", "p", "T", "u_j ×3"],
                 equations=7, unknowns=7),
}


def closure_count(stage: str = "full") -> dict:
    """Equations vs unknowns at each stage of the chapter (the counting of §4.4, §4.6, §4.8).

    Book: §4.4 after (4.24): 1 + 3 + 2 = 6 equations for 1 + 3 + 9 = 13 unknowns; §4.6 after (4.38): 4 equations for 5
    unknowns (closed when ρ is constant or p = p(ρ), barotropic); §4.8 after (4.60): 1 + 3 + 1 + 2 = 7 equations for 7
    unknowns (ρ, e, p, T, u_j).

    Parameters: stage ∈ {"empty", "cauchy", "navier_stokes", "barotropic", "full"} ("empty": zeros, the §4.1 table
    skeleton). Returns dict(equations (int), unknowns (int), equation_names, unknown_names, closed). Label: analytic.
    """
    if stage not in _CLOSURE:
        raise ValueError(f"stage must be one of {tuple(_CLOSURE)}")
    s = _CLOSURE[stage]
    out = {"equations": int(s.get("equations", 0)), "unknowns": int(s.get("unknowns", 0)),
           "equation_names": list(s["equation_names"]), "unknown_names": list(s["unknown_names"])}
    out["closed"] = out["equations"] == out["unknowns"] and out["equations"] > 0
    return out


# ======================================================================================================================
# §4.2 conservation of mass: the expanding flow (C01, C02), incompressible but stratified flow
# ======================================================================================================================
def expanding_flow(x, t: float, a: float = 1.0, rho0: float = 1.0, dim: int = 1):
    """An exact compressible flow: uniform expansion u = ax/(1 + at), ρ = ρ₀/(1 + at)^dim (our C01/C02 field).

    Book: §4.2, Eqs. (4.1), (4.5), (4.7) — particles move on x = X(1 + at), so a material interval [X₀, X₁](1 + at) keeps
    its mass ρ₀(X₁ − X₀) while a fixed box loses mass through its faces.

    Parameters: x : position(s) [m] (float or (N,) for dim = 1; (dim,) or (dim, N) otherwise); t [s] (1 + at > 0);
    a : expansion rate at t = 0 [1/s]; rho0 : density at t = 0 [kg/m³]; dim : 1, 2 or 3.
    Returns (u, rho): velocity (shape of x) [m/s] and density [kg/m³] (uniform in space).
    Validation: V2 sympy: (4.7) holds exactly; V4 :func:`material_mass` constant in t; V1 the interval ends follow
    ``core.kinematics.pathline``. Label: analytic, symbolic, conserved.
    """
    x_ = _F(x)
    s = 1.0 + float(a) * float(t)
    if s <= 0:
        raise ValueError("1 + a t must be > 0")
    u = float(a) * x_ / s  # u = a x/(1 + at)
    rho = float(rho0) * s ** (-int(dim)) * np.ones(x_.shape if dim == 1 else x_.shape[1:])  # ρ₀(1 + at)^(−dim)
    return _S(u), _S(rho)


def expanding_flow_fields(a: float = 1.0, rho0: float = 1.0, dim: int = 1):
    """Callables (rho(x, t), u(x, t)) of :func:`expanding_flow` in the ``core.kinematics`` layout ((dim,) or (dim, N)
    points) for ``core.conservation.mass_budget`` and ``core.navier_stokes.continuity_residual``. Book: §4.2, Eqs. (4.1), (4.7). Label: analytic."""
    def rho(x, t):
        return float(rho0) * (1.0 + a * t) ** (-dim) + 0.0 * _F(x)[0]

    def u(x, t):
        return a * _F(x) / (1.0 + a * t)
    return rho, u


def material_interval(X0: float, X1: float, t: float, a: float = 1.0):
    """Ends at time t of the material interval that was [X0, X1] at t = 0 in the expanding flow: X(1 + at) [m].
    Book: §4.2 (the balloon of (4.1) in 1-D). Label: analytic."""
    s = 1.0 + a * t
    return float(X0 * s), float(X1 * s)


def material_mass(X0: float, X1: float, t: float, a: float = 1.0, rho0: float = 1.0, dim: int = 1) -> float:
    """Mass per unit cross-section (dim = 1) of the material interval that was [X0, X1] at t = 0, integrated with
    ``quad`` over its current extent X(1 + at) [kg/m²]; dim > 1: the material cube of side X1 − X0 [kg·m^(dim−3)].

    Book: §4.2, Eq. (4.1): d/dt∫_{V(t)}ρ dV = 0. Validation: V4 (dim = 1) the ``quad`` integral over the moving interval
    equals ρ₀(X1 − X0) at every t (``material_mass(1.0, 2.0, 1.0)`` = 1.0); V1 (dim > 1) the cube's value is a closed-form
    identity ρ₀(1 + at)^(−dim)·((X1 − X0)(1 + at))^dim = ρ₀(X1 − X0)^dim (no independent computation). Label: conserved
    (dim = 1), analytic (dim > 1).
    """
    x0, x1 = material_interval(X0, X1, t, a)
    if dim == 1:
        return float(quad(lambda x: expanding_flow(x, t, a, rho0, 1)[1], x0, x1, epsabs=0.0, epsrel=1e-13)[0])
    s = 1.0 + a * t
    return float(rho0 * s ** (-dim) * ((X1 - X0) * s) ** dim)


def stratified_shear_flow(z, U0: float = 1.0, shear: float = 0.1, rho0: float = 1000.0, drho_dz: float = -0.5):
    """Incompressible but not constant-density flow: u = (U₀ + s z, 0, 0), ρ = ρ₀ + (dρ/dz) z (our C02 example).

    Book: §4.2, Eqs. (4.9)–(4.10): Dρ/Dt = 0 and ∇·u = 0 although ρ varies ("constant-density flows are a subset of
    incompressible flows"). Parameters: z [m]; U0 [m/s]; shear [1/s]; rho0 [kg/m³]; drho_dz [kg/m⁴].
    Returns (u_x, rho). Label: analytic.
    """
    z_ = _F(z)
    return _S(U0 + shear * z_), _S(rho0 + drho_dz * z_)


def stratified_shear_fields(U0: float = 1.0, shear: float = 0.1, rho0: float = 1000.0, drho_dz: float = -0.5):
    """Callables (rho(x, t), u(x, t)) of :func:`stratified_shear_flow` for 3-D points (x, y, z). Book: §4.2, Eqs. (4.9)–(4.10). Label: analytic."""
    def rho(x, t):
        return rho0 + drho_dz * _F(x)[2]

    def u(x, t):
        z = _F(x)[2]
        return np.stack([U0 + shear * z, 0.0 * z, 0.0 * z])
    return rho, u


def is_incompressible_regime(U, T: float = 288.15, threshold: float = 0.3) -> bool:
    """True if a gas flow at speed U is in the nearly incompressible regime M = U/√(γRT) < threshold (0.3).

    Book: §4.2 after (4.10) ("flow speeds less than ~100 m/s, that is, for Mach numbers < 0.3") and §4.11 after (4.111)
    (M = 0.3 → U²/c² ≈ 0.09). Label: analytic.
    """
    return bool(float(mach_number(U, T=T)) < float(threshold))


def incompressible_speed_limit(T: float = 288.15, M: float = 0.3, gamma: float = GAMMA_AIR, R: float = R_AIR) -> float:
    """Speed [m/s] at which a perfect gas at temperature T reaches Mach M (0.3 → ≈ 102 m/s at 288 K). Book: §4.2.
    Label: analytic."""
    return float(M * np.sqrt(gamma * R * T))


def ball_integral(f: Callable, x0, radius: float, epsrel: float = 1e-10) -> float:
    """∫ f dV over the ball of radius r about x0 by ``tplquad`` in spherical coordinates (the C02 localisation demo: as
    r → 0, ∫_ball f dV/V → f(x0), so an integral that vanishes for every ball forces the integrand to vanish, (4.6)).

    Book: §4.2, the argument after (4.6). Parameters: f(X) with X a 3-vector [any units]; x0 (3,) [m]; radius [m].
    Returns [f·m³]. Label: analytic.
    """
    c = _F(x0)

    def integrand(r, th, ph):
        X = c + r * np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])
        return float(f(X)) * r ** 2 * np.sin(th)
    return float(tplquad(integrand, 0.0, 2 * np.pi, 0.0, np.pi, 0.0, radius, epsrel=epsrel)[0])


# ======================================================================================================================
# §4.3 a second concrete 3-D pair of stream functions (Fig. 4.1)
# ======================================================================================================================
def stream_surface_example(y0: float = 0.5, a: float = 1.0, b: float = 2.0, c: float = 0.0, d: float = 1.0) -> dict:
    """A 3-D flow from two stream functions, χ = x² + y and ψ = z − y², ρu = ∇χ × ∇ψ = (1, −2x, −4xy), and a patch on
    the plane y = y₀ bounded by χ = a, b and ψ = c, d (Fig. 4.1's grey surface) — a second pair besides
    ``core.streamfunction.stream_function_pair("parabolic")``.

    Book: §4.3, Eq. (4.12) and ṁ = (b − a)(d − c) (our example). Returns dict(chi, psi (callables of X (3, N)), patch
    (callable (s, r) → X, oriented along −e_y, the local stream direction), expected = (b − a)(d − c)).
    Validation: V1 ``stream_tube_mass_flux(chi, psi, patch=patch)[0]`` = expected to 1e-9. Label: analytic.
    """
    if a - y0 < 0:
        raise ValueError("need a ≥ y0 so that x = √(χ − y0) is real")

    def chi(X):
        return X[0] ** 2 + X[1]

    def psi(X):
        return X[2] - X[1] ** 2

    def patch(s, r):
        chi_v = a + (b - a) * _F(s)  # s ↦ χ
        psi_v = c + (d - c) * _F(r)  # r ↦ ψ
        x = np.sqrt(chi_v - y0)
        return np.stack([x, y0 + 0.0 * x, psi_v + y0 ** 2])
    return {"chi": chi, "psi": psi, "patch": patch, "expected": (b - a) * (d - c)}


# ======================================================================================================================
# §4.4 momentum: body force potential, Examples 4.1–4.4, E1 scenarios
# ======================================================================================================================
def gravity_potential(z, g: float = G):
    """Force potential of gravity Φ = gz [J/kg] (z up), so g = −∇Φ = −g e_z, Eq. (4.18). Book: §4.4, Eq. (4.18). Label: analytic."""
    return _S(float(g) * _F(z))


def body_force_from_potential(Phi, x, h: float = 1e-4, g: float = G) -> np.ndarray:
    """Conservative body force per unit mass g = −∇Φ, Eq. (4.18) (central differences, step h [m]).

    Book: §4.4, Eq. (4.18). Parameters: Phi : callable Φ(X) of points (d, …) [J/kg] or "gravity" (Φ = g × last
    coordinate); x [m]; h [m]; g [m/s²] (for "gravity"). Returns g [m/s²] (shape of x).
    Validation: V1 Φ = gz → (0, 0, −g); the work around a closed loop is 0. Label: analytic.
    """
    from .core._stencil import grad
    if isinstance(Phi, str):
        if Phi != "gravity":
            raise ValueError("Phi must be a callable or 'gravity'")
        f = lambda X, T: float(g) * _F(X)[-1]  # noqa: E731
    else:
        f = lambda X, T: Phi(X)  # noqa: E731
    return -grad(f, _F(x), 0.0, h)  # Eq. (4.18)


def gaussian_wake(y, U_inf: float = 10.0, deficit: float = 2.0, width: float = 0.1):
    """A model wake profile U(y) = U∞ − Δ exp(−y²/b²) [m/s] (our synthetic 'measurement' for Example 4.1).
    Book: §4.4, Example 4.1, Fig. 4.2 (the book leaves U(y) general). Label: analytic."""
    return _S(U_inf - deficit * np.exp(-_F(y) ** 2 / width ** 2))


def _wake_fn(U_of_y, U_inf, p):
    if isinstance(U_of_y, str):
        if U_of_y != "gaussian":
            raise ValueError("U_of_y must be a callable or 'gaussian'")
        q = dict(deficit=p.get("deficit", 2.0), width=p.get("width", 0.1))
        return lambda y: float(gaussian_wake(y, U_inf, q["deficit"], q["width"]))
    return U_of_y


def wake_drag_per_span(U_of_y, U_inf: float = 10.0, rho: float = 1.0, H: float = 1.0, return_error: bool = False,
                       **p):
    """Drag per unit span from a wake survey, F_D/l = ρ∫_{−H/2}^{H/2} U(U∞ − U) dy [N/m] (Example 4.1).

    Book: §4.4, Example 4.1: mass and x-momentum balances of a fixed box whose top and bottom leak fluid ∫(U∞ − U)dy; the
    pressure p∞ integrates to zero (Gauss). F_D > 0 when U(y) < U∞ (drag on the body; −F_D acts on the fluid).

    Parameters: U_of_y : callable U(y) [m/s] or "gaussian" (keywords deficit, width of :func:`gaussian_wake`);
    U_inf [m/s]; rho [kg/m³]; H : box height [m] (≫ wake width); return_error : also return ``quad``'s estimate.
    Validation: V1 Gaussian deficit → ρ[U∞Δ√π b − Δ²√(π/2) b] (:func:`wake_drag_gaussian`) once H ≫ b
    (``wake_drag_per_span("gaussian", 10.0, 1.2, 2.0)`` = 3.6523 N/m); F → 0 as Δ → 0. Label: analytic.
    """
    U = _wake_fn(U_of_y, U_inf, p)
    val, err = quad(lambda y: U(y) * (U_inf - U(y)), -H / 2, H / 2, epsabs=0.0, epsrel=1e-12, limit=200)
    F = rho * val  # F_D/l = ρ∫U(U∞ − U)dy
    return (F, rho * err) if return_error else F


def wake_side_outflow(U_of_y, U_inf: float = 10.0, H: float = 1.0, **p) -> float:
    """Volume flux per span leaving through the top and bottom of the box, ∫_top v dx − ∫_bottom v dx = ∫(U∞ − U)dy
    [m²/s] (the mass balance of Example 4.1; Gaussian: Δ√π b = 0.35449 m²/s for Δ = 2, b = 0.1). Book: §4.4, Example 4.1. Label: analytic."""
    U = _wake_fn(U_of_y, U_inf, p)
    return float(quad(lambda y: U_inf - U(y), -H / 2, H / 2, epsabs=0.0, epsrel=1e-12, limit=200)[0])


def wake_drag_gaussian(U_inf: float, deficit: float, width: float, rho: float = 1.0) -> float:
    """Closed form of :func:`wake_drag_per_span` for :func:`gaussian_wake` and H → ∞:
    F_D/l = ρ[U∞Δ√π b − Δ²√(π/2) b] [N/m] (our sympy result). Scalar-callable (E1 parity). Book: §4.4, Example 4.1. Label: analytic."""
    return float(rho * (U_inf * deficit * np.sqrt(np.pi) * width - deficit ** 2 * np.sqrt(np.pi / 2) * width))


def stream_tube_element_balance_sym() -> dict:
    """Example 4.2 in sympy: the mass and streamwise momentum balances of a stream-tube element of length ds.

    Book: §4.4, Example 4.2, Fig. 4.3: −ρUA + ρ(U + U_s ds)(A + A_s ds) = 0 and the momentum balance with the extra
    pressure force on the conical wall; subtracting U × mass and dropping (ds)² gives U∂U/∂s ds = −g dz − (1/ρ)∂p/∂s ds,
    i.e. d(U²/2) + g dz + dp/ρ = 0 → (4.19).

    Returns dict(mass = O(ds) coefficient of the mass balance, momentum = O(ds) coefficient of (momentum − U·mass),
    dropped = the O(ds²) remainder, result = momentum/(ρA) = UU_s + g sin θ + p_s/ρ, residual_vs_4_19 → 0).
    Label: symbolic.
    """
    rho, U, A, p, g, th, ds = sp.symbols("rho U A p g theta ds", positive=True)
    Us, As, ps = sp.symbols("U_s A_s p_s")
    mass = -rho * U * A + rho * (U + Us * ds) * (A + As * ds)
    lhs = -rho * U ** 2 * A + rho * (U + Us * ds) ** 2 * (A + As * ds)
    rhs = (-rho * g * sp.sin(th) * (A + As * ds / 2) * ds + p * A + (p + ps * ds / 2) * As * ds
           - (p + ps * ds) * (A + As * ds))
    mom = sp.expand(lhs - rhs - U * mass)
    c1 = mom.coeff(ds, 1)
    result = sp.simplify(c1 / (rho * A))
    return {"mass": sp.expand(mass).coeff(ds, 1), "momentum": c1, "dropped": sp.expand(mom - c1 * ds),
            "result": result, "residual_vs_4_19": sp.simplify(result - (U * Us + g * sp.sin(th) + ps / rho))}


def bore_speed(h_in, h_out, g: float = G):
    """Speed of a bore (small solitary wave step) moving into still water, U = √(g h_out(h_in + h_out)/(2h_in)) [m/s].

    Book: §4.4, Example 4.3 (moving CV b = −U e_x; mass U h_in = (U_out + U)h_out; hydrostatic side pressures; p_o
    cancels); U ≈ √(gh) when the depths are close. Parameters: h_in (still side), h_out (disturbed side) [m]; g.
    Scalar-callable (``bore_speed(1.0, 1.1)`` = 3.3661 m/s). Validation: V1 → √(gh) as h_out → h_in (error O(Δh/h)); V1
    equals Bélanger's jump relation h₂/h₁ = (√(1 + 8Fr₁²) − 1)/2 in the wave frame (Wikipedia "Hydraulic jump", form).
    Label: analytic.
    """
    return _S(np.sqrt(float(g) * _F(h_out) * (_F(h_in) + _F(h_out)) / (2.0 * _F(h_in))))


def bore_outlet_velocity(h_in, h_out, g: float = G):
    """Fluid speed behind the bore, U_out = (h_in − h_out)U/h_out [m/s] (Example 4.3's mass balance; negative when
    h_out > h_in: the fluid follows the wave). Book: §4.4, Example 4.3. Label: analytic."""
    return _S((_F(h_in) - _F(h_out)) * _F(bore_speed(h_in, h_out, g)) / _F(h_out))


def bore_pressure_force(h_in: float, h_out: float, rho: float = 1000.0, g: float = G, p_o: float = 1.0e5,
                        width: float = 1.0) -> dict:
    """Horizontal pressure forces on Example 4.3's control volume (hydrostatic p = p_o + ρg(h − y) on the ends, p_o on
    the sloping top) [N].

    Book: §4.4, Example 4.3: left (inlet) ∫(p_o + ρg(h_in − y)) l dy, right (outlet) −∫(p_o + ρg(h_out − y)) l dy, top
    p_o(h_out − h_in)l; net = ρg(h_in²/2 − h_out²/2)l — p_o cancels. Returns dict(left, right, top, net) (``quad`` for
    the end integrals). Validation: V1 net independent of p_o. Label: analytic.
    """
    left = quad(lambda y: p_o + rho * g * (h_in - y), 0.0, h_in)[0] * width
    right = -quad(lambda y: p_o + rho * g * (h_out - y), 0.0, h_out)[0] * width
    top = p_o * (h_out - h_in) * width
    return {"left": left, "right": right, "top": top, "net": left + right + top}


def rocket_delta_v(M0, M1, Ve):
    """Speed gained with no gravity and no surface force, Δb = V_e ln(M₀/M₁) [m/s] (Example 4.4 with g = F_S = 0;
    Tsiolkovsky's equation, Wikipedia form). Book: §4.4, Example 4.4. Label: analytic."""
    return _S(_F(Ve) * np.log(_F(M0) / _F(M1)))


def rocket_closed_form(t, M0: float, mdot: float, Ve: float, g: float = G):
    """Closed-form trajectory of Example 4.4's equation with constant ṁ and F_S = 0 (vertical launch from rest):
    b(t) = −V_e ln(1 − ṁt/M₀) − gt, z(t) = V_e(M₀/ṁ)[(1 − x)ln(1 − x) + x] − ½gt², x = ṁt/M₀.
    Returns (z [m], b [m/s]). Book: §4.4, Example 4.4 (our integration). Label: analytic."""
    t_ = _F(t)
    x = mdot * t_ / M0
    b = -Ve * np.log1p(-x) - g * t_
    z = Ve * (M0 / mdot) * ((1.0 - x) * np.log1p(-x) + x) - 0.5 * g * t_ ** 2
    return _S(z), _S(b)


def rocket_trajectory(M0: float, mdot: float, Ve: float, t_burn: float, g: float = G, Fs: float = 0.0,
                      t_eval=None, rtol: float = 1e-10, atol: float = 1e-12) -> dict:
    """Vertical rocket flight during the burn from M d²z_R/dt² = −V_e dM/dt − Mg + F_S, dM/dt = −ṁ (Example 4.4).

    Book: §4.4, Example 4.4 (accelerating CV b = b(t)e_z; eliminating ρ_eV_eA_e between the mass and momentum balances).
    Our integration: ``solve_ivp`` (RK45) from rest at z = 0 to burn-out t_burn (an event stops it if the mass would
    reach zero).

    Parameters: M0 [kg]; mdot [kg/s] (> 0); Ve : exhaust speed relative to the rocket [m/s]; t_burn [s]; g [m/s²];
    Fs : constant surface force [N]; t_eval (default 201 points); rtol, atol.
    Returns dict(t, z, b, M, a) [s, m, m/s, kg, m/s²].
    Validation: V1 g = F_S = 0: b(t_burn) = V_e ln(M₀/M₁) to 1e-8; with g: equals :func:`rocket_closed_form`; V3 rtol
    sweep converges. Label: analytic, converged.
    """
    if M0 - mdot * t_burn <= 0:
        raise ValueError("the rocket would run out of mass before t_burn")
    t_eval = np.linspace(0.0, t_burn, 201) if t_eval is None else _F(t_eval)

    def rhs(t, y):
        M = M0 - mdot * t
        return [y[1], (Ve * mdot + Fs) / M - g]  # M d²z/dt² = −V_e dM/dt − Mg + F_S, dM/dt = −ṁ

    def empty(t, y):
        return M0 - mdot * t - 1e-9 * M0
    empty.terminal = True
    sol = solve_ivp(rhs, (0.0, float(t_burn)), [0.0, 0.0], t_eval=t_eval, rtol=rtol, atol=atol, events=empty)
    M = M0 - mdot * sol.t
    return {"t": sol.t, "z": sol.y[0], "b": sol.y[1], "M": M, "a": (Ve * mdot + Fs) / M - g}


def jet_plate_force(rho: float, V: float, A: float, theta: float = np.pi / 2) -> float:
    """Force of a free jet (speed V, area A) on a stationary flat plate, normal to the plate: F = ρV²A sin θ [N]
    (θ = angle between jet and plate; π/2 → ρV²A). Our E1 extension of (4.17) (fixed CV, atmospheric pressure all round,
    the sheets leave along the plate carrying no normal momentum). Book: §4.4, Eq. (4.17). Label: analytic."""
    return float(rho * V ** 2 * A * np.sin(theta))


CV_SCENARIOS = ("wake", "bore", "jet", "rocket")


def _face(name, area, un_rel, mass_flux, momentum_flux_x, momentum_flux_along=0.0):
    return dict(name=name, area=float(area), un_rel=float(un_rel), mass_flux=float(mass_flux),
                momentum_flux_x=float(momentum_flux_x), momentum_flux_along=float(momentum_flux_along))


def cv_scenario(name: str, **p) -> dict:
    """E1's four control-volume scenarios: per-face mass and momentum fluxes, storage, forces and the residuals of (4.5)
    and (4.17) (per unit span for the 2-D ones).

    Book: §4.4, Eqs. (4.5), (4.17), Examples 4.1 (wake), 4.3 (bore), 4.4 (rocket); "jet" (jet on a plate) is ours.

    * "wake" (U_inf = 10, deficit = 2, width = 0.1, rho = 1.2, H = 2): fixed box; faces inlet, outlet, top, bottom (the
      leakage ∫(U∞ − U)dy split evenly); result = drag per span F_D/l [N/m].
    * "bore" (h_in = 1, h_out = 1.1, rho = 1000, g, b = CV speed leftward, default the bore speed U): if b ≠ U the jump
      drifts inside the CV and the storage terms pick up the difference; result = bore speed U [m/s].
    * "jet" (rho = 1000, V = 10, A = 1e-3, theta = π/2): fixed CV; faces jet, sheet_up, sheet_down; momentum = plate-normal
      component; result = plate force [N]. With no plate shear and both sheets leaving at V, along-plate momentum splits
      the flow Q₁,₂ = Q(1 ± cos θ)/2 (sheet_up on the side of the jet's along-plate component); the along-plate budget
      (in ρV²A cos θ, out ρV²A(Q₁ − Q₂)/Q) is returned as ``residual_momentum_along`` (0). θ = π/4, ρ = 1000, V = 10,
      A = 1e-3: sheet mass fluxes 8.536 and 1.464 kg/s.
    * "rocket" (M0 = 1, mdot = 0.05, Ve = 500, t = 2, g): CV riding with the rocket; face nozzle; momentum = z-component;
      result = rocket speed b(t) [m/s].

    Returns dict(name, faces (list of dict(name, area, un_rel = (u − b)·n, mass_flux, momentum_flux_x,
    momentum_flux_along)), mass_out, momentum_out_x, storage_mass, storage_x, body_x, surface_x, residual_mass,
    residual_momentum, residual_momentum_along (jet only; 0 otherwise), result, result_label).
    Validation: V1 both residuals vanish (to quad accuracy) for every scenario; b ≠ U in "bore" still closes.
    Label: analytic.
    """
    if name == "wake":
        Ui, De, w = float(p.get("U_inf", 10.0)), float(p.get("deficit", 2.0)), float(p.get("width", 0.1))
        H, rho = float(p.get("H", 2.0)), float(p.get("rho", 1.2))
        U = lambda y: float(gaussian_wake(y, Ui, De, w))  # noqa: E731
        Q_out = quad(U, -H / 2, H / 2, epsrel=1e-12)[0]
        M_out = quad(lambda y: U(y) ** 2, -H / 2, H / 2, epsrel=1e-12)[0]
        leak = wake_side_outflow(U, Ui, H)
        faces = [_face("inlet", H, -Ui, -rho * Ui * H, -rho * Ui ** 2 * H),
                 _face("outlet", H, Q_out / H, rho * Q_out, rho * M_out),
                 _face("top", 1.0, 0.5 * leak, 0.5 * rho * leak, 0.5 * rho * Ui * leak),
                 _face("bottom", 1.0, 0.5 * leak, 0.5 * rho * leak, 0.5 * rho * Ui * leak)]
        FD = wake_drag_per_span(U, Ui, rho, H)
        sm = sx = body = 0.0
        surface = -FD  # the bar pushes the fluid upstream; p∞ integrates to 0
        result, label = FD, "drag per span F_D/l [N/m]"
    elif name == "bore":
        hi, ho = float(p.get("h_in", 1.0)), float(p.get("h_out", 1.1))
        rho, g = float(p.get("rho", 1000.0)), float(p.get("g", G))
        Ub = float(bore_speed(hi, ho, g))
        b = float(p.get("b", Ub))
        Uo = (hi - ho) * Ub / ho
        faces = [_face("inlet", hi, -b, -rho * b * hi, 0.0),
                 _face("outlet", ho, Uo + b, rho * (Uo + b) * ho, rho * Uo * (Uo + b) * ho)]
        sm = rho * (ho - hi) * (Ub - b)  # the jump drifts at (b − U) relative to a CV moving at −b
        sx = rho * Uo * ho * (Ub - b)
        body = 0.0
        surface = rho * g * (hi ** 2 - ho ** 2) / 2.0  # net hydrostatic force (p_o cancels)
        result, label = Ub, "bore speed U [m/s]"
    elif name in ("jet", "jet_plate"):
        rho, V, A = float(p.get("rho", 1000.0)), float(p.get("V", 10.0)), float(p.get("A", 1e-3))
        th = float(p.get("theta", np.pi / 2))
        c = np.cos(th)
        # split fixed by along-plate momentum (no plate shear, both sheets leave at V): Q₁,₂ = Q(1 ± cos θ)/2;
        # sheet_up is the side toward which the jet's along-plate component V cos θ points
        A1, A2 = A * (1.0 + c) / 2.0, A * (1.0 - c) / 2.0
        faces = [_face("jet", A, -V, -rho * V * A, -rho * V ** 2 * A * np.sin(th), -rho * V ** 2 * A * c),
                 _face("sheet_up", A1, V, rho * V * A1, 0.0, rho * V ** 2 * A1),
                 _face("sheet_down", A2, V, rho * V * A2, 0.0, -rho * V ** 2 * A2)]
        sm = sx = body = 0.0
        Fp = jet_plate_force(rho, V, A, th)
        surface = -Fp  # the plate pushes the fluid back along its normal
        result, label = Fp, "force on the plate [N]"
    elif name == "rocket":
        M0, mdot, Ve = float(p.get("M0", 1.0)), float(p.get("mdot", 0.05)), float(p.get("Ve", 500.0))
        t, g = float(p.get("t", 2.0)), float(p.get("g", G))
        M = M0 - mdot * t
        b = float(rocket_closed_form(t, M0, mdot, Ve, g)[1])
        a = Ve * mdot / M - g
        faces = [_face("nozzle", float("nan"), Ve, mdot, mdot * (-Ve + b))]  # (u − b)·n dA = V_e dA, ρ_eV_eA_e = ṁ
        sm = -mdot  # dM/dt
        sx = M * a - mdot * b  # d(Mb)/dt = M db/dt + b dM/dt
        body = -M * g
        surface = 0.0
        result, label = b, "rocket speed b(t) [m/s]"
    else:
        raise ValueError(f"unknown scenario {name!r}; choose from {CV_SCENARIOS}")
    mo = sum(f["mass_flux"] for f in faces)
    po = sum(f["momentum_flux_x"] for f in faces)
    pa = sum(f["momentum_flux_along"] for f in faces)  # along-plate (jet) momentum: no shear, no storage ⇒ 0
    return {"name": name, "faces": faces, "mass_out": float(mo), "momentum_out_x": float(po),
            "storage_mass": float(sm), "storage_x": float(sx), "body_x": float(body), "surface_x": float(surface),
            "residual_mass": float(sm + mo), "residual_momentum": float(sx + po - body - surface),
            "residual_momentum_along": float(pa), "result": float(result), "result_label": label}


def cube_spin_acceleration(tau12, tau21, rho, h):
    """Angular acceleration of a small cube of side h if the shear stresses were not symmetric:
    dΩ₃/dt = 6(τ₁₂ − τ₂₁)/(ρh²) [rad/s²] (torque (τ₁₂ − τ₂₁)h³ over moment of inertia ρh⁵/6).

    Book: §4.5, Eq. (4.25) and Exercise 4.30 (our D08): as h → 0 it would diverge like h⁻² unless τ₁₂ = τ₂₁ — the stress
    tensor is symmetric (no body couples). Scalar-callable (``cube_spin_acceleration(1.0, 0.0, 1000.0, 0.01)`` = 60).
    Label: analytic.
    """
    return _S(6.0 * (_F(tau12) - _F(tau21)) / (_F(rho) * _F(h) ** 2))


# ======================================================================================================================
# §4.6 exact profiles used in figures
# ======================================================================================================================
def stokes_first_problem(y, t, U: float = 1.0, nu: float = 1e-6):
    """Stokes' first problem: a plate at y = 0 set moving at U at t = 0, u = U erfc(y/(2√(νt))) [m/s].

    Book: §4.6 (an exact solution of (4.39b): local ∂u/∂t balances ν∂²u/∂y²; derived in Ch. 8) and §4.11 (u/U = erfc(η/2),
    η = y/√(νt), is the same for every U, ν: similarity). Scalar-callable.
    Validation: V1 NS residual ≈ 0 (``exact_solution("stokes_first")``); two (U, ν, t) collapse on erfc(η/2).
    Label: analytic.
    """
    return _S(float(U) * erfc(_F(y) / (2.0 * np.sqrt(float(nu) * _F(t)))))


def plane_poiseuille(y, G: float = 100.0, h: float = 0.01, mu: float = 1e-3):
    """Plane Poiseuille profile u = Gy(h − y)/(2μ) [m/s] between walls y = 0 and h, G = −dp/dx [Pa/m] (C08 figure).
    Name clash: the parameter ``G`` here is the pressure gradient −dp/dx, **not** the module constant ``ch04.G`` (= g,
    gravity, 9.81 m/s²). Book: §4.6 (pressure and viscous terms of (4.39b) balance: μ d²u/dy² = −G; Ch. 8 derives it).
    Label: analytic."""
    y_ = _F(y)
    return _S(G * y_ * (h - y_) / (2.0 * mu))


# ======================================================================================================================
# §4.7 rotating frames: the pole projectile, Example 4.5, highs and lows
# ======================================================================================================================
def coriolis_projectile(u0: float = 100.0, Omega: float = OMEGA_EARTH, t_eval=None, centrifugal: bool = True,
                        rtol: float = 1e-11, atol: float = 1e-12) -> dict:
    """The pole projectile of Fig. 4.8: forward distance, deflection and turn angle, and the path integrated in the
    rotating frame a′ = −2Ω × u′ − Ω × (Ω × x′) from x′ = 0 with u′ = (u0, 0) (frictionless, gravity balanced).

    Book: §4.7, text before Fig. 4.8: forward distance ut, deflection Ωut², angular deflection Ωut²/ut = Ωt. Exactly
    (our D17) the lateral offset seen in the rotating frame is ut sin Ωt → Ωut² as Ωt → 0.

    Parameters: u0 [m/s]; Omega [rad/s] (> 0 NH, < 0 SH); t_eval [s] (default 201 points over 3600 s; any times ≥ 0);
    centrifugal : include −Ω × (Ω × x′) (then the ODE path is exactly the inertial straight line seen from the turning
    frame; without it the speed is conserved exactly — the Coriolis force does no work); rtol, atol.
    Returns dict(t, forward = u₀t, deflection = u₀t sin|Ω|t (exact, to the right in the NH), deflection_small = |Ω|u₀t²,
    angle = Ωt [rad], rotating (2, N) from the ODE, inertial (2, N) straight path, speed (N,) in the rotating frame).
    Validation: V1 ``coriolis_projectile(10.0, OMEGA_EARTH, [3600.0])`` → forward 36 000 m, deflection 9342.4 m,
    deflection_small 9450.6 m, angle 0.26252 rad; ODE path = projectile_paths to 1e-9; V4 speed constant without
    centrifugal; NH right, SH left. Label: analytic, conserved.
    """
    t_eval = np.linspace(0.0, 3600.0, 201) if t_eval is None else np.atleast_1d(_F(t_eval))
    W = float(Omega)
    c2 = W ** 2 if centrifugal else 0.0

    def rhs(t, y):
        x, yy, u, v = y
        return [u, v, 2.0 * W * v + c2 * x, -2.0 * W * u + c2 * yy]  # −2Ω×u′ − Ω×(Ω×x′), Ω = (0, 0, W)
    t_end = float(np.max(t_eval))
    if t_end > 0:
        sol = solve_ivp(rhs, (0.0, t_end), [0.0, 0.0, float(u0), 0.0], t_eval=np.sort(t_eval), method="DOP853",
                        rtol=rtol, atol=atol)
        rot, spd = sol.y[:2], np.hypot(sol.y[2], sol.y[3])
    else:
        rot, spd = np.zeros((2, t_eval.size)), np.full(t_eval.size, float(u0))
    tt = np.sort(t_eval)
    inertial, _ = projectile_paths(u0, W, tt)
    return {"t": tt, "forward": u0 * tt, "deflection": u0 * tt * np.sin(abs(W) * tt),
            "deflection_small": abs(W) * u0 * tt ** 2, "angle": W * tt, "rotating": rot, "inertial": inertial,
            "speed": spd}


def rotating_pump_terms() -> dict:
    """Example 4.5 by sympy (dict form): ρ(u·∇)u, the right side −∇p + ρ[−2Ω × u − Ω × (Ω × x)] + μ∇²u and the rotation
    terms [ρ(2Ω_zu_φ + Ω_z²R), −2ρΩ_zu_R, 0] in cylindrical components (``core.curvilinear``). Book: §4.7, Example 4.5.
    Returns dict(symbols, lhs, rhs, rotation_terms, residual). Label: symbolic."""
    R, phi, z = CU.coordinates("cylindrical")
    rho, mu, Om = sp.symbols("rho mu Omega_z", positive=True)
    uR, uphi, uz, pf = (sp.Function(n)(R, z) for n in ("u_R", "u_phi", "u_z", "p"))
    u = [uR, uphi, uz]
    lhs = [sp.expand(rho * a) for a in CU.advective_acceleration(u, "cylindrical", simplify=False)]
    Omv = [0, 0, Om]
    x = [R, 0, z]  # position vector R e_R + z e_z in physical components
    cor = CU.cross(Omv, u)
    cen = CU.cross(Omv, CU.cross(Omv, x))
    rot = [sp.expand(rho * (-2 * cor[i] - cen[i])) for i in range(3)]
    gp = CU.gradient(pf, "cylindrical", simplify=False)
    lap = CU.vector_laplacian(u, "cylindrical", simplify=False)
    rhs = [sp.expand(-gp[i] + rot[i] + mu * lap[i]) for i in range(3)]
    return {"symbols": dict(R=R, z=z, rho=rho, mu=mu, Omega_z=Om, u_R=uR, u_phi=uphi, u_z=uz, p=pf),
            "lhs": lhs, "rhs": rhs, "rotation_terms": rot,
            "residual": [sp.simplify(lhs[i] - rhs[i]) for i in range(3)]}


def rotating_pump_equations() -> list:
    """Example 4.5: the radial, azimuthal and axial momentum equations of steady, constant-ρ, constant-μ, axisymmetric
    flow in a von Kármán impeller pump rotating at Ω_z (frame of the disks, no body force), as three sympy Eq.

    Book: §4.7, Example 4.5: ρ(u′·∇′)u′ = −∇′p + ρ[−2Ω × u′ − Ω × (Ω × x′)] + μ∇′²u′ with Ω × u = Ω_zu_Re_φ − Ω_zu_φe_R,
    Ω × (Ω × x′) = −Ω_z²Re_R, Appendix-B cylindrical operators (``core.curvilinear``).
    Returns [Eq_R, Eq_φ, Eq_z] (lhs = ρ(u·∇)u component). Validation: V2 each equals the book's written equation
    (difference simplifies to 0). Label: symbolic.
    """
    d = rotating_pump_terms()
    return [sp.Eq(d["lhs"][i], d["rhs"][i]) for i in range(3)]


def high_low_flow(x, y, U_R: float = 1.0, sense: str = "high"):
    """Radial flow out of a high (sense="high") or into a low ("low") with speed U_R: (u, v) = ±U_R (x, y)/r [m/s] — the
    N62 quiver, before the Coriolis force −2Ω_zu_Re_φ turns it clockwise (high) or counter-clockwise (low) in the NH.
    Book: §4.7, text after Fig. 4.8. Label: analytic."""
    x_, y_ = _F(x), _F(y)
    r = np.hypot(x_, y_)
    s = 1.0 if sense == "high" else -1.0 if sense == "low" else None
    if s is None:
        raise ValueError("sense must be 'high' or 'low'")
    with np.errstate(invalid="ignore", divide="ignore"):
        u = np.where(r > 0, s * U_R * x_ / np.where(r > 0, r, 1.0), 0.0)
        v = np.where(r > 0, s * U_R * y_ / np.where(r > 0, r, 1.0), 0.0)
    return _S(u), _S(v)


# ======================================================================================================================
# §4.8 energy: Couette heating (E7)
# ======================================================================================================================
def couette_heating(y, U: float = 1.0, h: float = 0.01, mu: float = 1e-3, k: float = 0.6, T0: float = 293.15,
                    dpdx: float = 0.0) -> dict:
    """Steady plane Couette–Poiseuille flow with viscous heating, both walls at T₀: velocity, dissipation, temperature,
    wall heat fluxes and the energy budget.

    Book: §4.8, Eqs. (4.58)–(4.60): steady, fully developed, ∇·u = 0 ⇒ 0 = μ(du/dy)² + k d²T/dy² (4.60). Our worked
    field for E7. u = Uy/h + (G/2μ)y(h − y), G = −dp/dx; T from integrating −(μ/k)(du/dy)² twice with T(0) = T(h) = T₀.
    Pure Couette (G = 0): T = T₀ + (μU²/2k)(y/h)(1 − y/h), ΔT_max = μU²/8k, heat out = μU²/h.

    Parameters: y heights [m] in [0, h]; U top-wall speed [m/s]; h gap [m]; mu [Pa s]; k [W/(m K)]; T0 [K]; dpdx [Pa/m].
    Returns dict(u [m/s], dudy [1/s], eps = ρε = μ(du/dy)² [W/m³], T [K], q_bottom, q_top [W/m²] conducted out of the
    fluid into each wall, work_in [W/m²] (shear work of the moving wall + pressure work G∫u dy), heat_out = q_bottom +
    q_top, dT_max [K]).
    Validation: V1 T satisfies (4.60) (``internal_energy_residual`` ≈ 0); V4 work_in = heat_out = ∫ρε dy;
    ``couette_heating(y, U=1, h=1e-3, mu=1e-3, k=0.6)`` → dT_max = 2.0833e-4 K, heat_out = 1.0 W/m². Label: analytic,
    conserved.
    """
    y_ = _F(y)
    Gp = -float(dpdx)
    B = Gp / (2.0 * mu)
    u = U * y_ / h + B * y_ * (h - y_)
    dudy = U / h + B * (h - 2.0 * y_)
    a0, a1 = U / h + B * h, -2.0 * B  # du/dy = a0 + a1 y
    c0, c1, c2 = a0 ** 2, 2 * a0 * a1, a1 ** 2  # (du/dy)² = c0 + c1 y + c2 y²
    F2 = lambda yy: -(c0 * yy ** 2 / 2 + c1 * yy ** 3 / 6 + c2 * yy ** 4 / 12)  # noqa: E731  F'' = −(du/dy)²
    C1 = -F2(h) / h
    Tfun = lambda yy: T0 + (mu / k) * (F2(yy) + C1 * yy)  # noqa: E731
    dTdy0 = (mu / k) * C1
    dTdyh = (mu / k) * (-(c0 * h + c1 * h ** 2 / 2 + c2 * h ** 3 / 3) + C1)
    q_b, q_t = float(k * dTdy0), float(-k * dTdyh)  # q·n_out at y = 0 (n = −e_y) and y = h (n = +e_y)
    tau_top = mu * (a0 + a1 * h)
    flow = U * h / 2.0 + B * h ** 3 / 6.0  # ∫u dy
    work_in = U * tau_top + Gp * flow  # moving-wall shear work + pressure work
    # maximum of T: at a real root in [0, h] of dT/dy ∝ −(c2/3)y³ − (c1/2)y² − c0 y + C1 = 0, or at a wall
    roots = np.roots([-c2 / 3.0, -c1 / 2.0, -c0, C1]) if (c1 or c2) else np.array([C1 / c0 if c0 else 0.0])
    cand = [0.0, h] + [float(r.real) for r in np.atleast_1d(roots) if abs(r.imag) < 1e-12 and 0.0 <= r.real <= h]
    dTmax = float(max(Tfun(np.array(cand))) - T0)
    return {"u": _S(u), "dudy": _S(dudy), "eps": _S(mu * dudy ** 2), "T": _S(Tfun(y_)), "q_bottom": q_b,
            "q_top": q_t, "work_in": float(work_in), "heat_out": q_b + q_t, "dT_max": dTmax}


def couette_heating_transient(y, t, U: float = 1.0, h: float = 0.01, mu: float = 1e-3, k: float = 0.6,
                              rho: float = 1000.0, cp: float = 4182.0, T0: float = 293.15, nterms: int = 200,
                              dpdx: float = 0.0, tol: float = 1e-12):
    """Temperature T(y, t) of the gap heating up from T₀ once the (steady) shear flow is on: ρC_p∂T/∂t = k∂²T/∂y² + φ(y),
    T(0) = T(h) = T₀, T(y, 0) = T₀ — E7's transient.

    Book: §4.8 (4.60) with e = C_pT for a liquid (the Boussinesq (4.88) with ρε kept). Our series solution:
    T = T_s(y) − Σ_{n=1}^{nterms} b_n sin(nπy/h) e^{−κ(nπ/h)²t}, κ = k/ρC_p, with b_n = (2/h)∫₀ʰ(T_s − T₀) sin(nπy/h) dy
    in **closed form**: T_s − T₀ = Σ_m a_m y^m is the quartic of :func:`couette_heating`, and
    I_m(n) = ∫₀ʰ y^m sin(ay) dy (a = nπ/h) follows from I_0 = (1 − (−1)^n)/a, I_m = −h^m(−1)^n/a + (m/a)J_{m−1},
    J_m = ∫₀ʰ y^m cos(ay) dy = −(m/a)I_{m−1} (exact for every n — no quadrature aliasing). All ``nterms`` terms are
    summed (even-n coefficients vanish by symmetry for pure Couette, so the sum never stops at a zero term); terms
    below ``tol``·(T_s − T₀)_max after decay are skipped. Accuracy: b_n ∝ n⁻³, so the truncation error is largest at
    t = 0, about 1/nterms² of ΔT_max (≲ 1e-4 ΔT_max for the default nterms = 200), and falls off like e^{−κ(nπ/h)²t}.

    Parameters: y [m] (array); t [s] (scalar); U, h, mu, k, T0, dpdx as :func:`couette_heating`; rho [kg/m³];
    cp [J/(kg K)]; nterms : number of sine modes; tol : relative size below which a decayed term is skipped.
    Returns T [K] (shape of y).
    Validation: V1 t → ∞ gives :func:`couette_heating`; t = 0 gives T₀ to the truncation error above. Label: analytic.
    """
    y_ = _F(y)
    kap = k / (rho * cp)
    # T_s − T₀ = (μ/k)(F2(y) + C1 y), F2 = −(c0 y²/2 + c1 y³/6 + c2 y⁴/12) — the polynomial of couette_heating
    B = -float(dpdx) / (2.0 * mu)
    a0, a1 = U / h + B * h, -2.0 * B
    c0, c1, c2 = a0 ** 2, 2 * a0 * a1, a1 ** 2
    C1 = (c0 * h ** 2 / 2 + c1 * h ** 3 / 6 + c2 * h ** 4 / 12) / h
    coef = (mu / k) * np.array([0.0, C1, -c0 / 2, -c1 / 6, -c2 / 12])  # a_m, m = 0 … 4
    n = np.arange(1, int(nterms) + 1, dtype=float)
    a = n * np.pi / h
    sgn = (-1.0) ** n
    I = [(1.0 - sgn) / a]  # I_0
    J = [np.zeros_like(a)]  # J_0 = sin(nπ)/a = 0
    for m in range(1, 5):
        I.append(-h ** m * sgn / a + (m / a) * J[m - 1])
        J.append(-(m / a) * I[m - 1])
    bn = (2.0 / h) * sum(coef[m] * I[m] for m in range(5))  # exact sine coefficients
    theta = _F(couette_heating(y_, U, h, mu, k, T0, dpdx)["T"]) - T0
    scale = max(float(np.max(np.abs(np.polyval(coef[::-1], np.linspace(0, h, 201))))), 1e-300)
    decayed = bn * np.exp(-kap * a ** 2 * float(t))
    keep = np.abs(decayed) >= tol * scale
    modes = np.sin(np.multiply.outer(a[keep], y_))  # (n_kept,) + y.shape
    theta = theta - np.tensordot(decayed[keep], modes, axes=(0, 0))
    return _S(T0 + theta)


# ======================================================================================================================
# §4.9 special forms: sprinkler, orifice, Rankine B, U-tube, accelerating sphere, Boussinesq
# ======================================================================================================================
def sprinkler_torque(a: float, rho: float, A: float, U: float, alpha: float) -> float:
    """Torque needed to hold a two-arm lawn sprinkler still, M = 2aρAU² cos α [N m] (Example 4.6, (4.65)).

    Book: §4.9, Example 4.6: stationary CV, atmospheric pressure all round, steady: ∮(r × ρu)(u·n)dA = M; each jet
    contributes (aρU cos α)UA. Parameters: arm radius a [m]; rho [kg/m³]; nozzle area A [m²]; jet speed U [m/s]; jet angle
    α above the horizontal [rad]. Scalar-callable (``sprinkler_torque(0.2, 1000.0, 1e-4, 5.0, π/6)`` = 0.8660 N m).
    Validation: V1 equals :func:`sprinkler_torque_numeric`; α = π/2 → 0. Label: analytic.
    """
    return float(2.0 * a * rho * A * U ** 2 * np.cos(alpha))  # Example 4.6


def sprinkler_torque_numeric(a: float, rho: float, A: float, U: float, alpha: float, n: int = 16) -> float:
    """z-component of ∮(r × ρu)(u·n)dA over two small discs (the nozzle exits) — the numerical route to
    :func:`sprinkler_torque` via ``core.conservation.angular_momentum_flux_surface``. [N m]. Book: §4.9, Eq. (4.65), Example 4.6. Label: analytic."""
    from .core.integral_theorems import planar_disc
    total = 0.0
    radius = np.sqrt(A / np.pi)
    for sgn in (1.0, -1.0):
        c = np.array([sgn * a, 0.0, 0.0])
        e = np.array([0.0, sgn * np.cos(alpha), np.sin(alpha)])  # jet direction (tangential + upward)
        disc = planar_disc(c, e, radius, n, 2 * n)
        ufield = lambda X, T, e=e: np.broadcast_to((U * e)[:, None], _F(X).shape)  # noqa: E731
        total += float(CB.angular_momentum_flux_surface(rho, ufield, disc.points.T, disc.normals.T, disc.dA)[2])
    return total


def sprinkler_free_spin_rate(a: float, U: float, alpha: float) -> float:
    """Rotation rate [rad/s] of a frictionless sprinkler: the jets leave with zero absolute tangential velocity,
    ω = U cos α/a (our extension of Example 4.6: zero torque ⇒ zero angular-momentum outflow). Book: §4.9, Example 4.6. Label: analytic."""
    return float(U * np.cos(alpha) / a)


def orifice_mass_flow(h, A, rho: float = 1000.0, Cc: float = 1.0, g: float = G):
    """Mass flow from a tank orifice at depth h, ṁ = ρC_cA√(2gh) [kg/s] (text after Fig. 4.17; C_c ≈ 0.61 for a sharp
    edge — Wikipedia "Vena contracta" 0.611 — and 1 for a rounded one). Scalar-callable. Book: §4.9 (orifice paragraph, Figs. 4.16–4.17). Label: analytic."""
    return _S(float(rho) * float(Cc) * _F(A) * np.sqrt(2.0 * float(g) * _F(h)))


def tank_drain_time(h0: float, A_tank: float, A_orifice: float, Cc: float = 1.0, g: float = G) -> float:
    """Time to empty a tank through an orifice in the quasi-steady Torricelli model, t = (A_tank/(C_cA_o))√(2h₀/g) [s]
    (our integration of dh/dt = −(C_cA_o/A_tank)√(2gh)). Book: §4.9 (orifice paragraph). Label: analytic."""
    return float(A_tank / (Cc * A_orifice) * np.sqrt(2.0 * h0 / g))


def tank_drain(h0: float, A_tank: float, A_orifice: float, Cc: float = 1.0, g: float = G, t_eval=None,
               rtol: float = 1e-10, atol: float = 1e-14) -> dict:
    """Water level of a draining tank, dh/dt = −(C_cA_o/A_tank)√(2gh) (quasi-steady (4.19) at each instant; valid when
    A_tank ≫ A_o, text before Fig. 4.16), by ``solve_ivp`` with a "tank empty" event.

    Book: §4.9 orifice paragraph (our ODE). Returns dict(t, h, t_empty (event), t_empty_closed).
    Validation: V1 t_empty = :func:`tank_drain_time`; h(t) = (√h₀ − (C_cA_o/A_t)√(g/2)t)²; V7 drain time grows with
    A_tank. Label: analytic, converged.
    """
    tc = tank_drain_time(h0, A_tank, A_orifice, Cc, g)
    t_eval = np.linspace(0.0, 0.999 * tc, 200) if t_eval is None else _F(t_eval)
    kk = Cc * A_orifice / A_tank * np.sqrt(2.0 * g)

    def rhs(t, y):
        return [-kk * np.sqrt(max(y[0], 0.0))]

    def empty(t, y):
        return y[0] - 1e-12 * h0
    empty.terminal = True
    sol = solve_ivp(rhs, (0.0, 1.5 * tc), [h0], t_eval=t_eval[t_eval <= 1.5 * tc], rtol=rtol, atol=atol,
                    events=empty, method="LSODA")
    te = float(sol.t_events[0][0]) if len(sol.t_events[0]) else float("nan")
    return {"t": sol.t, "h": sol.y[0], "t_empty": te, "t_empty_closed": tc}


def rankine_bernoulli(r, Gamma: float = 2 * np.pi, sigma: float = 1.0, rho: float = 1000.0, p_inf: float = 0.0) -> dict:
    """Steady Rankine vortex (core radius σ): u_θ (3.28), pressure from the radial balance dp/dr = ρu_θ²/r, and the
    Bernoulli function B = ½u_θ² + p/ρ — constant on each circle (streamline), varying across circles inside the core
    (rotational: only (4.71) holds), uniform outside ((4.72)).

    Book: §4.9, Eqs. (4.70)–(4.72) (C11 figure; our field). Parameters: r [m]; Gamma [m²/s]; sigma [m]; rho; p_inf [Pa].
    Returns dict(u_theta [m/s], p [Pa], B [m²/s²]). Validation: V1 ``rankine_bernoulli(0.0, 2π, 1.0)["B"]`` = −1.0,
    ``(2.0, …)["B"]`` = 0.0; B inside = p∞/ρ − Γ²/(4π²σ²) + Γ²r²/(4π²σ⁴). Label: analytic.
    """
    r_ = _F(r)
    ut = np.where(r_ <= sigma, Gamma * r_ / (2 * np.pi * sigma ** 2),
                  Gamma / (2 * np.pi * np.where(r_ > 0, r_, 1.0)))
    p = _F(rankine_vortex_pressure(r_, Gamma, sigma, rho, p_inf))
    return {"u_theta": _S(ut), "p": _S(p), "B": _S(0.5 * ut ** 2 + p / rho)}


def u_tube_column(t, L: float = 1.0, h0: float = 0.05, g: float = G, rho: float = 1000.0) -> dict:
    """Frictionless oscillation of a liquid column of length L in a U-tube: surface displacement h = h₀cos ωt,
    ω = √(2g/L), and the unsteady streamline Bernoulli equation (4.82) between the two free surfaces.

    Book: §4.9, Eq. (4.82): ∫₁²∂u/∂t·ds + (½u² + gz + p/ρ)₂ − (…)₁ = 0 with ∫∂u/∂t ds = L dU/dt, p₁ = p₂ = p_atm,
    z₂ − z₁ = 2h ⇒ L dU/dt + 2gh = 0 (our example; column speed U = dh/dt).
    Returns dict(h [m], U [m/s], dUdt [m/s²], dp = ρL dU/dt [Pa] (the pressure difference the column's acceleration
    needs over its length), omega [rad/s], residual_4_82 [J/kg] (≈ 0)); scalar t → floats.
    Validation: V1 residual ≈ 0 (``core.bernoulli.unsteady_streamline_bernoulli``). Label: analytic.
    """
    t_ = _F(t)
    w = np.sqrt(2.0 * g / L)
    x = h0 * np.cos(w * t_)
    U = -h0 * w * np.sin(w * t_)
    dU = -h0 * w ** 2 * np.cos(w * t_)
    s = np.linspace(0.0, L, 11)
    res = np.vectorize(lambda xx, uu, du: unsteady_streamline_bernoulli(
        np.full(s.shape, du), s, dict(U=uu, z=-xx, p=0.0), dict(U=uu, z=xx, p=0.0), rho, g))(x, U, dU)
    return {"h": _S(x), "U": _S(U), "dUdt": _S(dU), "dp": _S(rho * L * dU), "omega": float(w),
            "residual_4_82": _S(res)}


def accelerating_sphere_pressure(theta, a: float = 0.1, dUdt: float = 1.0, rho: float = 1000.0, U: float = 0.0,
                                 p_inf: float = 0.0) -> dict:
    """Surface pressure on a sphere of radius a moving at speed U(t) with acceleration dU/dt through still ideal fluid,
    p = p∞ + ρa(dU/dt)cos θ/2 + ρU²(9cos²θ − 5)/8 [Pa] (θ from the direction of motion), its net force and added mass.

    Book: §4.9, Eq. (4.75) with φ = −U(t)a³cos θ/(2r²) about the moving centre (our example: ∂φ/∂t at a fixed point
    includes −U∂φ/∂x from the moving centre). The dU/dt term integrates to −½(4πa³ρ/3)dU/dt (added mass: half the
    displaced mass); the U² term to zero (d'Alembert).
    Returns dict(p (at theta) [Pa], force = −∮p cos θ dA [N] (quad), added_mass = ½ρ(4/3)πa³ [kg]).
    Validation: V1 force = −added_mass·dU/dt (Wikipedia "Added mass", form); B of (4.74) uniform with
    :func:`accelerating_sphere_fields`. Label: analytic.
    """
    def p_of(th):
        c = np.cos(_F(th))
        return p_inf + 0.5 * rho * a * dUdt * c + rho * U ** 2 * (9.0 * c ** 2 - 5.0) / 8.0
    f = lambda th: float(p_of(th)) * np.cos(th) * 2 * np.pi * a ** 2 * np.sin(th)  # noqa: E731
    force = float(-quad(f, 0.0, np.pi, epsabs=0.0, epsrel=1e-11)[0])
    return {"p": _S(p_of(theta)), "force": force, "added_mass": 0.5 * rho * 4.0 / 3.0 * np.pi * a ** 3}


def accelerating_sphere_fields(a: float = 0.1, U0: float = 1.0, dUdt: float = 0.5, rho: float = 1000.0,
                               p_inf: float = 0.0):
    """(φ(x, t), p(x, t)) of a sphere moving along x with U(t) = U₀ + (dU/dt)t from x = 0 at t = 0, in still fluid:
    φ = −U a³(x − X)/(2|x − X|³) (dipole), p = p∞ − ρ(∂φ/∂t + ½|∇φ|²) (4.75). Points (3, N) outside the sphere.
    Book: §4.9, Eqs. (4.73)–(4.75) (our example). Label: analytic."""
    def X_of(t):
        return U0 * t + 0.5 * dUdt * t ** 2

    def U_of(t):
        return U0 + dUdt * t

    def rel(x, t):
        x_ = _F(x)
        r = x_ - np.array([X_of(t), 0.0, 0.0]).reshape((3,) + (1,) * (x_.ndim - 1))
        return r, np.sqrt(np.sum(r ** 2, axis=0))

    def phi(x, t):
        r, rr = rel(x, t)
        return -U_of(t) * a ** 3 * r[0] / (2.0 * rr ** 3)

    def p(x, t):
        r, rr = rel(x, t)
        c = r[0] / rr
        U = U_of(t)
        ur = U * a ** 3 * c / rr ** 3  # u_r = ∂φ/∂r
        s = np.sqrt(np.maximum(1.0 - c ** 2, 0.0))
        ut = U * a ** 3 * s / (2.0 * rr ** 3)  # u_θ = (1/r)∂φ/∂θ
        ux = ur * c - ut * s
        phit = -dUdt * a ** 3 * c / (2.0 * rr ** 2) - U * ux  # ∂φ/∂t|_x
        return p_inf - rho * (phit + 0.5 * (ur ** 2 + ut ** 2))  # Eq. (4.75)
    return phi, p


def boussinesq_density(T, rho0: float = 1000.0, alpha: float = 2.1e-4, T0: float = 293.15):
    """Linear equation of state of the Boussinesq set, ρ = ρ₀[1 − α(T − T₀)] [kg/m³] (§4.9 summary).
    Book: §4.9, text after (4.89). Validation: V1 equals ``ch01.seawater_density_linear`` with the salinity term off.
    Label: analytic."""
    return _S(rho0 * (1.0 - alpha * (_F(T) - T0)))


def boussinesq_validity(alpha: float, dT: float, L: float, U: float, c: float = 340.0, g: float = G,
                        nu: float = 1.5e-5, cp: float = 1004.5, small: float = 0.1) -> dict:
    """The small parameters of the Boussinesq approximation and the verdict.

    Book: §4.9: (i) (1/ρ)(Dρ/Dt)/∇·u ~ δρ/ρ = αδT ≪ 1; (ii) L ≪ H_c = c²/g (the height over which hydrostatic pressure
    changes the density, ~10 km in air); (iii) viscous heating ρε/(ρC_pDT/Dt) ~ νU/(C_pδT L) ≪ 1 (typically ~1e-7);
    (iv) low Mach number, judged by the §4.2 rule M < 0.3 (U²/c² < 0.09, as :func:`is_incompressible_regime`).

    Parameters: alpha [1/K]; dT [K]; L [m]; U [m/s]; c [m/s]; g [m/s²]; nu [m²/s]; cp [J/(kg K)]; small : threshold
    for (i)–(iii).
    Returns dict(alpha_dT, H_c [m], L_over_Hc, heating_ratio, g_prime = gαδT [m/s²], mach, valid (bool), verdict (text)).
    Validation: V1 ``boussinesq_validity(2e-4, 10.0, 10.0, 0.1, c=1500.0, nu=1e-6, cp=4186.0)`` → alpha_dT 2.0e-3,
    g_prime 0.01962 m/s². Label: analytic.
    """
    Hc = c ** 2 / g
    out = {"alpha_dT": alpha * dT, "H_c": Hc, "L_over_Hc": L / Hc, "heating_ratio": nu * U / (cp * dT * L),
           "g_prime": g * alpha * dT, "mach": U / c}
    names = {"alpha_dT": "αδT not small", "L_over_Hc": "deep layer: L ≈ c²/g", "heating_ratio": "viscous heating",
             "mach": "Mach number not small"}
    fails = [names[k] for k in ("alpha_dT", "L_over_Hc", "heating_ratio") if out[k] >= small]
    if out["mach"] >= 0.3:  # §4.2: M < 0.3 ⇔ nearly incompressible (same rule as is_incompressible_regime)
        fails.append(names["mach"])
    out["valid"] = not fails
    out["verdict"] = "Boussinesq valid" if not fails else "not Boussinesq: " + ", ".join(fails)
    return out


BOUSSINESQ_SCENARIOS = {
    # typical orders of magnitude (ours, not the book's): α [1/K], δT [K], L [m], U [m/s], ν [m²/s], C_p, c [m/s], T0 [K]
    "lake": dict(fluid="water", alpha=1.5e-4, dT=5.0, L=10.0, U=0.05, nu=1.1e-6, cp=4186.0, c=1480.0, T0=288.15),
    "thermocline": dict(fluid="sea water", alpha=2.0e-4, dT=10.0, L=200.0, U=0.1, nu=1.0e-6, cp=3990.0, c=1500.0,
                        T0=283.15),
    "lab_tank": dict(fluid="water", alpha=2.1e-4, dT=2.0, L=0.3, U=0.01, nu=1.0e-6, cp=4182.0, c=1480.0, T0=293.15),
    "abl": dict(fluid="air", alpha=1.0 / 300.0, dT=3.0, L=1000.0, U=5.0, nu=1.5e-5, cp=1005.0, c=340.0, T0=300.0),
    "deep_atmosphere": dict(fluid="air", alpha=1.0 / 250.0, dT=40.0, L=10000.0, U=20.0, nu=1.5e-5, cp=1005.0, c=320.0,
                            T0=250.0),
}


def boussinesq_scenario(name: str) -> dict:
    """E8's five scenarios (typical values, ours): "lake", "thermocline", "lab_tank", "abl" (atmospheric boundary layer),
    "deep_atmosphere" (fails L ≪ c²/g). Returns dict(fluid, alpha, dT, L, U, nu, cp, c, T0) — feed them to
    :func:`boussinesq_validity`. Book: §4.9 (validity conditions). Label: analytic."""
    if name not in BOUSSINESQ_SCENARIOS:
        raise ValueError(f"choose from {tuple(BOUSSINESQ_SCENARIOS)}")
    return dict(BOUSSINESQ_SCENARIOS[name])


def gaussian_blob_advection_diffusion(x, t, U=0.0, kappa: float = 1.4e-7, sigma0: float = 0.01, amp: float = 1.0,
                                      dim: int = 2, T0: float = 0.0, x0=0.0):
    """Exact solution of DT/Dt = κ∇²T (4.89) in a uniform current U: a Gaussian carried by the flow and spreading,
    T′ = T − T₀ = A(σ₀²/σ²)^{dim/2} exp(−|x − x₀ − Ut|²/(2σ²)), σ² = σ₀² + 2κt [K].

    Book: §4.9, Eq. (4.89) (our test field: the heat kernel shifted by the advection). Parameters: x — float/(N,) for
    dim = 1 or (dim,)/(dim, N); t [s]; U [m/s] (scalar along x₁ or a vector); kappa [m²/s]; sigma0 [m]; amp [K]; T0 [K];
    x0 [m]. Validation: V1 ``temperature_equation_residual`` ≈ 0; V4 ∫(T − T₀)dV = A(2πσ₀²)^{dim/2}; the centre moves at
    U. Label: analytic, conserved.
    """
    x_ = _F(x)
    s2 = sigma0 ** 2 + 2.0 * kappa * _F(t)
    if dim == 1:
        r2 = (x_ - _F(x0) - _F(U) * _F(t)) ** 2
    else:
        Uv = np.zeros(dim) if np.ndim(U) == 0 else _F(U)
        if np.ndim(U) == 0:
            Uv[0] = float(U)
        x0v = np.zeros(dim) if np.ndim(x0) == 0 else _F(x0)
        shape = (dim,) + (1,) * (x_.ndim - 1)
        r2 = np.sum((x_ - (x0v + Uv * float(t)).reshape(shape)) ** 2, axis=0)
    return _S(T0 + amp * (sigma0 ** 2 / s2) ** (dim / 2.0) * np.exp(-r2 / (2.0 * s2)))


def blob_rise(t, g_prime: float, tau_d: float = 10.0) -> dict:
    """Rise of a buoyant blob in **our one-line parcel model for E8** (not the book's): dw/dt = g′ − w/τ_d (buoyancy
    g′ = gαδT minus a linear drag with time scale τ_d), from rest at z = 0:
    w = g′τ_d(1 − e^{−t/τ_d}), z = g′τ_d(t − τ_d(1 − e^{−t/τ_d})).

    Book: §4.9, Eq. (4.86) motivates the buoyancy term (ρ′/ρ₀)g; the drag closure is ours. Parameters: t [s];
    g_prime [m/s²]; tau_d [s]. Returns dict(w [m/s], z [m]); scalar t → floats (E8 parity).
    Validation: V1 dw/dt − (g′ − w/τ_d) = 0 by differentiation; w → g′τ_d as t → ∞. Label: analytic.
    """
    t_ = _F(t)
    e = -np.expm1(-t_ / tau_d)  # 1 − e^{−t/τ_d}
    return {"w": _S(g_prime * tau_d * e), "z": _S(g_prime * tau_d * (t_ - tau_d * e))}


def gaussian_blob_fields(U=0.0, kappa: float = 1.4e-7, sigma0: float = 0.01, amp: float = 1.0, dim: int = 2,
                         T0: float = 0.0):
    """Callables (T(x, t), u(x, t)) of :func:`gaussian_blob_advection_diffusion` for (dim, N) points. Book: §4.9, Eq. (4.89). Label: analytic."""
    Uv = np.zeros(dim) if np.ndim(U) == 0 else _F(U)
    if np.ndim(U) == 0:
        Uv[0] = float(U)

    def T(x, t):
        return gaussian_blob_advection_diffusion(x, t, Uv, kappa, sigma0, amp, dim, T0)

    def u(x, t):
        x_ = _F(x)
        return np.broadcast_to(Uv.reshape((dim,) + (1,) * (x_.ndim - 1)), x_.shape).copy()
    return T, u


# ======================================================================================================================
# §4.10 boundary conditions: two-layer problems, slip, a surface wave, the meniscus (Example 4.7)
# ======================================================================================================================
def two_layer_conduction(y, k1: float, k2: float, L1: float, L2: float, T_hot: float, T_cold: float) -> dict:
    """Steady conduction through two layers in contact: T continuous and k₁∂T₁/∂n = k₂∂T₂/∂n at the interface.

    Book: §4.10 (pillbox, Fig. 4.18: continuity of heat flux; no temperature jump T₁ = T₂). Our example: layer 1 on
    0 ≤ y ≤ L₁ (T(0) = T_hot), layer 2 on L₁ ≤ y ≤ L₁ + L₂ (T_cold at the top); q = (T_hot − T_cold)/(L₁/k₁ + L₂/k₂).
    Returns dict(T (at y) [K], q [W/m²], T_interface [K]). Validation: V1 flux continuous, slopes in ratio k₂/k₁.
    Label: analytic.
    """
    y_ = _F(y)
    q = (T_hot - T_cold) / (L1 / k1 + L2 / k2)
    Ti = T_hot - q * L1 / k1
    T = np.where(y_ <= L1, T_hot - q * y_ / k1, Ti - q * (y_ - L1) / k2)
    return {"T": _S(T), "q": float(q), "T_interface": float(Ti)}


def two_fluid_couette(y, mu1: float, mu2: float, h1: float, h2: float, U: float) -> dict:
    """Two immiscible layers sheared between a fixed wall (y = 0) and a wall moving at U (y = h₁ + h₂): velocity and
    shear stress continuous at the interface (no slip between the fluids; n_iτ_ij continuous).

    Book: §4.10 (interface conditions from the pillbox + no slip). τ = U/(h₁/μ₁ + h₂/μ₂); u piecewise linear.
    Returns dict(u (at y) [m/s], tau [Pa], u_interface [m/s]). Label: analytic.
    """
    y_ = _F(y)
    tau = U / (h1 / mu1 + h2 / mu2)
    ui = tau * h1 / mu1
    u = np.where(y_ <= h1, tau * y_ / mu1, ui + tau * (y_ - h1) / mu2)
    return {"u": _S(u), "tau": float(tau), "u_interface": float(ui)}


def navier_slip_couette(y, U: float = 1.0, h: float = 1e-3, slip_length: float = 0.0):
    """Couette flow with a Navier slip condition u = b du/dy on the fixed wall y = 0 (no slip on the moving wall y = h):
    u = U(y + b)/(h + b) [m/s] (our extension of §4.10's no-slip discussion: slip lengths of textured/rarefied
    walls; b = 0 recovers no slip). Book: §4.10 (no-slip and its violations). Label: analytic."""
    y_ = _F(y)
    return _S(U * (y_ + slip_length) / (h + slip_length))


def linear_wave_surface(x, z, t, a: float = 0.05, k: float = 1.0, g: float = G, H: float = np.inf) -> dict:
    """Linear progressive surface wave (a test field for the kinematic condition; Ch. 7 derives it):
    η = a cos(kx − ωt), ω² = gk tanh kH; u = aω[cosh k(z + H)/sinh kH]cos(kx − ωt), w = aω[sinh k(z + H)/sinh kH]
    sin(kx − ωt) (deep water: e^{kz} for both).

    Book: §4.10, Eqs. (4.90)–(4.91) (the moving surface of B1 and the C14 animation). Parameters: x, z [m] (z up, still
    surface z = 0); t [s]; a amplitude [m]; k [1/m]; g; H depth [m] (inf = deep). Returns dict(eta, u, w, omega).
    Validation: V1 the linearised condition ∂η/∂t = w(z = 0) exactly; the full (4.91) on z = η has residual O((ka)²).
    Label: analytic.
    """
    x_, z_, t_ = _F(x), _F(z), _F(t)
    om = float(np.sqrt(g * k * (np.tanh(k * H) if np.isfinite(H) else 1.0)))
    ph = k * x_ - om * t_
    if np.isfinite(H):
        fu = np.cosh(k * (z_ + H)) / np.sinh(k * H)
        fw = np.sinh(k * (z_ + H)) / np.sinh(k * H)
    else:
        fu = fw = np.exp(k * z_)
    return {"eta": _S(a * np.cos(ph)), "u": _S(a * om * fu * np.cos(ph)), "w": _S(a * om * fw * np.sin(ph)),
            "omega": om}


def linear_wave_fields(a: float = 0.05, k: float = 1.0, g: float = G, H: float = np.inf):
    """Callables for ``core.interfaces``: η_imp(x, t) = z − η(x, t) and u(x, t) = (u, w) for points (2, N) = (x, z).
    Book: §4.10 (our test field; = ``core.interfaces.surface_preset("linear_wave")``). Label: analytic."""
    return surface_preset("linear_wave", a=a, k=k, g=g, H=H)


def wave_kinematic_residual(x, t, a: float = 0.05, k: float = 1.0, g: float = G, H: float = np.inf,
                            full: bool = True, h: float = 1e-6):
    """Dη/Dt of the linear wave evaluated on its own surface: full=True → the exact condition (4.91) at z = η(x, t)
    (∂η_imp/∂t + u·∇η_imp with η_imp = z − η), O((ka)²aω); full=False → the linearised condition ∂η/∂t − w at z = 0
    (exactly 0). [m/s]. Book: §4.10, Eqs. (4.90)–(4.91) (B1 parity). Label: analytic."""
    x_ = _F(x)
    if not full:
        d = linear_wave_surface(x_, 0.0, t, a, k, g, H)
        om = d["omega"]
        deta_dt = a * om * np.sin(k * x_ - om * _F(t))
        return _S(deta_dt - _F(d["w"]))
    eta, u = linear_wave_fields(a, k, g, H)
    zs = _F(linear_wave_surface(x_, 0.0, t, a, k, g, H)["eta"])
    X = np.stack([np.atleast_1d(x_), np.atleast_1d(zs)])
    return _S(_F(kinematic_bc_residual(eta, u, X, t, h, 1e-6 / max(1.0, float(np.sqrt(g * k)))))
              .reshape(np.shape(x_)))


def spheroid_area(volume: float, aspect: float) -> float:
    """Surface area of a spheroid of given volume and aspect ratio c/a (polar/equatorial semi-axis) [m²].

    Prolate (c > a): A = 2πa²(1 + (c/(ae)) arcsin e), e² = 1 − a²/c²; oblate (c < a): A = 2πa²(1 + ((1 − e²)/e) artanh e),
    e² = 1 − c²/a²; sphere 4πa²; V = (4/3)πa²c. Book: §4.10 (Batchelor's remark: at fixed volume the sphere has the least
    area — why drops are round). Our illustration. Validation: V1 minimum at aspect 1 (4πr²). Label: analytic.
    """
    a = (3.0 * volume / (4.0 * np.pi * aspect)) ** (1.0 / 3.0)
    c = aspect * a
    if abs(aspect - 1.0) < 1e-12:
        return float(4.0 * np.pi * a ** 2)
    if c > a:
        e = np.sqrt(1.0 - a ** 2 / c ** 2)
        return float(2.0 * np.pi * a ** 2 * (1.0 + c / (a * e) * np.arcsin(e)))
    e = np.sqrt(1.0 - c ** 2 / a ** 2)
    return float(2.0 * np.pi * a ** 2 * (1.0 + (1.0 - e ** 2) / e * np.arctanh(e)))


def meniscus_height(theta, sigma: float = 0.0728, rho: float = 998.0, g: float = G):
    """Rise of a liquid at a vertical wall with contact angle θ, h = √((2σ/ρg)(1 − sin θ)) [m] (Example 4.7).

    Book: §4.10, Example 4.7: first integral (ρg/2σ)ζ² + (1 + ζ′²)^{−1/2} = 1 with ζ′(0) = −cot θ. θ = 90° → flat;
    θ = 0 → √2 δ, δ = √(σ/ρg). Scalar-callable. Label: analytic.
    """
    return _S(np.sqrt(2.0 * sigma / (rho * g) * (1.0 - np.sin(_F(theta)))))


def meniscus_profile_x(zeta, theta, sigma: float = 0.0728, rho: float = 998.0, g: float = G):
    """Closed-form meniscus shape x(ζ) at a vertical wall (Example 4.7), 0 < ζ ≤ h:
    x/δ = cosh⁻¹(2δ/ζ) − (4 − ζ²/δ²)^{1/2} − cosh⁻¹(2δ/h) + (4 − h²/δ²)^{1/2}, δ² = σ/ρg.

    Book: §4.10, Example 4.7 (the book's separated equation drops a minus sign — the slope is negative — and says
    "η = h/δ" for γ = h/δ; the final formula is right: checked by differentiation, analysis §9 item 6).
    Returns x [m]. Validation: V2 sympy dx/dζ = −1/|ζ′| from the curvature ODE; V1 matches :func:`meniscus_profile_ode`
    to 1e-7. Label: analytic, symbolic.
    """
    d = np.sqrt(sigma / (rho * g))
    h = float(meniscus_height(theta, sigma, rho, g))
    z = _F(zeta)
    return _S(d * (np.arccosh(2 * d / z) - np.sqrt(4 - z ** 2 / d ** 2) - np.arccosh(2 * d / h)
                   + np.sqrt(4 - h ** 2 / d ** 2)))


def meniscus_profile_ode(theta, x_max: float = 0.01, sigma: float = 0.0728, rho: float = 998.0, g: float = G,
                         n: int = 200, rtol: float = 1e-11, atol: float = 1e-14):
    """Meniscus profile ζ(x) from the curvature balance ρgζ = σκ integrated from the wall (independent of the closed form).

    Book: §4.10, Example 4.7: ρgζ/σ − ζ″/(1 + ζ′²)^{3/2} = 0 with ζ′(0) = tan(θ + π/2) = −cot θ, ζ → 0 far away.
    # DEVIATION: integrated in arc length s (dx/ds = cos ψ, dζ/ds = sin ψ, dψ/ds = ζ/δ², ψ(0) = θ − π/2, ζ(0) = h)
    instead of the book's dζ/dx separation — the slope −cot θ is infinite for θ → 0 and the arc-length form has no
    singularity (the stiffness risk of analysis §9).

    Parameters: theta [rad]; x_max [m]; sigma [N/m]; rho [kg/m³]; g; n output points; rtol, atol.
    Returns (x, zeta) arrays [m] on x ∈ [0, x_max].
    Validation: V1 equals :func:`meniscus_profile_x` (inverted) to 1e-7 relative; the first integral
    ζ²/(2δ²) + cos ψ = 1 is constant along the solution. Label: analytic, converged.
    """
    d2 = sigma / (rho * g)
    h = float(meniscus_height(theta, sigma, rho, g))
    psi0 = float(theta) - np.pi / 2

    def rhs(s, y):
        return [np.cos(y[2]), np.sin(y[2]), y[1] / d2]

    def reach(s, y):
        return y[0] - x_max
    reach.terminal = True
    s_end = 10.0 * (x_max + h + np.sqrt(d2))
    sol = solve_ivp(rhs, (0.0, s_end), [0.0, h, psi0], events=reach, dense_output=True, rtol=rtol, atol=atol,
                    method="DOP853")
    ss = np.linspace(0.0, sol.t[-1], 20 * n)
    Y = sol.sol(ss)
    xg = np.linspace(0.0, x_max, n)
    order = np.argsort(Y[0])
    return xg, np.interp(xg, Y[0][order], Y[1][order])


# ======================================================================================================================
# §4.11 similarity: Prandtl numbers, the ship model (Example 4.8), sphere-drag data, Π groups
# ======================================================================================================================
def prandtl_of(fluid: str = "air", T: float = 293.15) -> float:
    """Prandtl number Pr = ν/κ = μC_p/k of air or water at temperature T (other ``ch01.FLUIDS`` entries at 20 °C only).

    Book: §4.11, Eq. (4.116) (the book quotes air and water values — private JSON). μ(T) from ``ch01.fluid_properties``
    (Sutherland for air, Vogel for water); C_p and k interpolated linearly in T from the ch01 Incropera rows (air
    250–300 K, water 290–295 K; mild extrapolation outside — k and C_p vary slowly).
    Validation: V5 air 300 K within 0.70–0.73 (Wikipedia "Prandtl number"); water 293.15 K ≈ 7.0 (engineering tables,
    2 %). Label: benchmark.
    """
    from . import ch01_introduction as ch01
    f = fluid.lower()
    if f in ("air", "water"):
        rows = ch01._AIR_A4 if f == "air" else ch01._WATER_A6
        (T1, r1), (T2, r2) = sorted(rows.items())
        w = (T - T1) / (T2 - T1)
        cp = (1 - w) * r1["cp"] + w * r2["cp"]
        k = (1 - w) * r1["k"] + w * r2["k"]
        props = ch01.fluid_properties(f, T)
        return float(prandtl_number(props["nu"], k / (props["rho"] * cp)))
    if abs(T - 293.15) > 1e-9:
        raise ValueError("only air and water have temperature-dependent properties here")
    p = ch01.FLUIDS[f]
    return float(prandtl_number(p["nu"], p["k"] / (p["rho"] * p["cp"])))


def ship_drag_extrapolation(L_p: float = 100.0, U_p: float = 10.0, S_p: float = 300.0, scale: float = 1.0 / 25.0,
                            D_m_total: float = 60.0, CDf_m: float = 0.003, CDf_p: float = 0.0015,
                            rho_m: float = 1000.0, rho_p: float = 1025.0, g: float = G, nu: float = 1e-6) -> dict:
    """Froude-scaled model test with a skin-friction correction (Example 4.8).

    Book: §4.11, Example 4.8: U_m = U_p√(l_m/l_p); model friction ½C_Df,mρ_mU_m²S_m; wave drag = total − friction, which
    scales as ρU²l² between Froude-similar flows: D_wave,p = D_wave,m(ρ_p/ρ_m)(l_p/l_m)²(U_p/U_m)²; plus the prototype
    friction ½C_Df,pρ_pU_p²S_p. The example's inputs are the defaults except ρ_p (the book uses 1000 kg/m³ for both;
    pass rho_p=1000.0 to reproduce it; 1025 is sea water).

    Returns dict(U_m, D_m_friction, D_m_wave, D_p_wave, D_p_friction, D_p_total, Re_ratio = Re_p/Re_m, Re_m, Re_p,
    D_p_uncorrected (model drag scaled as if all wave drag)) [SI].
    Validation: V1 wave-drag ratio (ρ_p/ρ_m)λ⁻³ (λ = scale) exactly; Re_p/Re_m = λ^(−3/2); V6 the book's numbers with
    rho_p = 1000 (private JSON, 0.5 %; the book's total 9.14e5 N is 9.15e5 before rounding). Label: analytic, book-value.
    """
    L_m = L_p * scale
    U_m = float(froude_scaled_speed(U_p, L_p, L_m, g, g))
    S_m = S_p * scale ** 2
    fm = 0.5 * CDf_m * rho_m * U_m ** 2 * S_m
    wm = D_m_total - fm
    k = (rho_p / rho_m) * (L_p / L_m) ** 2 * (U_p / U_m) ** 2
    wp = wm * k
    fp = 0.5 * CDf_p * rho_p * U_p ** 2 * S_p
    Re_m, Re_p = U_m * L_m / nu, U_p * L_p / nu
    return {"U_m": U_m, "D_m_friction": fm, "D_m_wave": wm, "D_p_wave": wp, "D_p_friction": fp, "D_p_total": wp + fp,
            "Re_ratio": Re_p / Re_m, "Re_m": Re_m, "Re_p": Re_p, "D_p_uncorrected": D_m_total * k}


def synthetic_sphere_drag_data(n: int = 40, seed: int = 0, noise: float = 0.03, Re_range=(0.1, 1e6)) -> dict:
    """Synthetic "experiments" for the Fig. 4.21 collapse: spheres of many diameters in air, water and glycerine
    (``ch01.FLUIDS`` at 20 °C) at speeds chosen to span Re_range; drag from the Morrison correlation ×
    (1 + noise·N(0, 1)). Seeded (``seed``).

    Book: §4.11, (4.99) and Fig. 4.21 (data from five dimensional parameters collapse on C_D(Re)).
    Returns dict(fluid, d [m], U [m/s], rho, mu, F [N], Re, CD) arrays. Label: analytic (synthetic data).
    """
    from .ch01_introduction import FLUIDS
    rng = np.random.default_rng(seed)
    Re = 10 ** rng.uniform(np.log10(Re_range[0]), np.log10(Re_range[1]), n)
    names = rng.choice(["air", "water", "glycerine"], n)
    d = 10 ** rng.uniform(-3, -0.5, n)
    rho = np.array([FLUIDS[f]["rho"] for f in names])
    mu = np.array([FLUIDS[f]["mu"] for f in names])
    U = Re * mu / (rho * d)
    CD = _F(sphere_drag_coefficient(Re)) * (1.0 + noise * rng.standard_normal(n))
    F = CD * 0.5 * rho * U ** 2 * np.pi * d ** 2 / 4.0
    return {"fluid": names, "d": d, "U": U, "rho": rho, "mu": mu, "F": F, "Re": Re, "CD": CD}


def sphere_drag_pi_groups() -> dict:
    """The two dimensionless scaling laws of (4.99) from the Π theorem with two repeating sets (``core.dimensional``):
    set 1 (U, D, ρ) → F/(ρU²D²) and μ/(ρUD); set 2 (μ, D, ρ) → Fρ/μ² and ρUD/μ (the book writes its inverse).
    Book: §4.11, Eq. (4.99) (ch01 D28). Returns dict(set1, set2) of group lists (exact Fractions).
    Validation: V1 Fρ/μ² = (F/ρU²D²)·Re² exactly. Label: analytic."""
    from .core.dimensional import SPHERE_DRAG, pi_groups
    return {"set1": pi_groups(SPHERE_DRAG, solution="F", repeating=("U", "D", "rho")),
            "set2": pi_groups(SPHERE_DRAG, solution="F", repeating=("mu", "D", "rho"))}


_NOT_EXPORTED = {"np", "sp", "quad", "solve_ivp", "tplquad", "erfc", "annotations", "Callable", "as_scalar_if_0d"}
