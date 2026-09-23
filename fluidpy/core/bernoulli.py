"""The Bernoulli equations of §4.9 with their hypotheses, and their classic applications (pitot tube, orifice,
stagnation quantities, U-tube, accelerating body).

Book: Kundu, Cohen & Dowling 5e, Ch. 4 — Example 4.2, Eq. (4.19); §4.9 Eqs. (4.66)–(4.83), Figs. 4.13–4.17
(rendered pages chapters/pages/ch04/p132–p133, p155–p161).

"Bernoulli" is **several statements, not one** (the summary after (4.83)):

| eq. | constant along | hypotheses |
|---|---|---|
| (4.19) ½U² + gz + p/ρ | a streamline | steady, inviscid, constant ρ |
| (4.71) ½u² + ∫dp/ρ + gz | streamlines and vortex lines (Lamb surfaces) | steady, inviscid, barotropic |
| (4.72) the same | everywhere | … and irrotational |
| (4.75) ∂φ/∂t + ½|∇φ|² + ∫dp/ρ + gz | everywhere at one instant | inviscid, irrotational, barotropic (unsteady ok) |
| (4.78) h + ½|u|² + gz | streamlines | steady, no viscous stress, no heat conduction (isentropic) |
| (4.82)/(4.83) | a streamline at one instant / everywhere | constant μ and ρ, irrotational (viscous flow) |

z is up, g = 9.81 m/s² by default; ∫dp/ρ is the pressure function of (4.67) from a reference p_o. Book typo handled:
(4.74) prints the gauge change φ → φ + ∫B dt′; absorbing B(t) needs φ → φ − ∫B dt′ (with + the bracket becomes 2B(t))
— analysis §9 item 3 (:func:`gauge_absorbed_bracket`).
"""
from __future__ import annotations

from typing import Callable

import numpy as np
from scipy.integrate import quad

from ._stencil import ddt, ev, grad, laplacian
from ._util import as_scalar_if_0d
from .kinematics import vorticity
from .thermo import G_BOOK as G0, GAMMA_AIR, R_AIR

__all__ = ["bernoulli_head", "bernoulli_solve", "pressure_function", "bernoulli_function", "lamb_surface_check",
           "lamb_surface_terms", "bernoulli_along_line", "unsteady_bernoulli_B", "unsteady_bernoulli_pressure",
           "gauge_absorbed_bracket", "viscous_irrotational_residual", "unsteady_streamline_bernoulli",
           "stagnation_enthalpy", "stagnation_temperature", "stagnation_pressure", "dynamic_pressure", "pitot_speed",
           "pitot_speed_from_heads", "torricelli_speed", "BERNOULLI_FORMS", "which_bernoulli", "which_bernoulli_text",
           "BERNOULLI_SCENARIOS", "bernoulli_scenario", "rankine_vortex_pressure"]

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d
CP_BOOK = 1004.5  #: C_p of air used by the Ch. 4 defaults [J/(kg K)] (ch01's derived CP_AIR is 1004.7)


# ======================================================================================================================
# (4.19) and the pressure function (4.67)
# ======================================================================================================================
def bernoulli_head(U, z, p, rho, g: float = G0):
    """Bernoulli constant of steady, inviscid, constant-density flow, ½U² + gz + p/ρ [J/kg = m²/s²], Eq. (4.19).

    Book: Example 4.2, Eq. (4.19) ("a constant along a streamline"). The example's statement writes (½)ρU² + gz + p/ρ,
    dimensionally inconsistent (J/m³ + J/kg); the derived (4.19) is correct (analysis §9 item 7).
    Parameters: U speed [m/s]; z height [m] (up); p pressure [Pa]; rho [kg/m³]; g [m/s²]. Scalar-callable.
    Validation: V1 constant along streamline traces of the cylinder flow; round trip with :func:`bernoulli_solve`.
    Label: analytic.
    """
    return _S(0.5 * _F(U) ** 2 + float(g) * _F(z) + _F(p) / _F(rho))  # Eq. (4.19)


def bernoulli_solve(state1: dict, state2: dict, unknown: str, rho: float = 1000.0, g: float = G0) -> float:
    """Solve (4.19) between two points of one streamline for the missing quantity of point 2.

    Book: Example 4.2, Eq. (4.19). ``state1`` has keys U [m/s], z [m], p [Pa]; ``state2`` the other two of point 2;
    ``unknown`` ∈ {"U", "z", "p"} names the one solved for. Raises if the requested speed would be imaginary.
    Label: analytic.
    """
    B = bernoulli_head(state1["U"], state1["z"], state1["p"], rho, g)
    s2 = dict(state2)
    if unknown == "U":
        k = 2.0 * (B - g * s2["z"] - s2["p"] / rho)
        if k < 0:
            raise ValueError("no real speed: point 2 has more head in p and z than point 1 in total")
        return float(np.sqrt(k))
    if unknown == "z":
        return float((B - 0.5 * s2["U"] ** 2 - s2["p"] / rho) / g)
    if unknown == "p":
        return float(rho * (B - 0.5 * s2["U"] ** 2 - g * s2["z"]))
    raise ValueError("unknown must be 'U', 'z' or 'p'")


def pressure_function(p, p_o, kind="constant", rho=None, T=None, R: float = R_AIR, gamma: float = GAMMA_AIR,
                      rho_o=None):
    """The barotropic pressure function P(p) = ∫_{p_o}^{p} dp′/ρ(p′) [J/kg], Eq. (4.67).

    Book: §4.9, Eq. (4.67): (1/ρ)∂p/∂x_j = ∂/∂x_j ∫_{p_o}^p dp′/ρ(p′) when ρ = ρ(p) (barotropic). p′ here is a dummy
    integration variable (not the perturbation of (4.84)).

    Parameters
    ----------
    p, p_o : pressure and reference pressure [Pa]
    kind : "constant" (ρ given) → (p − p_o)/ρ; "isothermal" (T given, perfect gas) → RT ln(p/p_o); "isentropic"
        (ρ_o at p_o given) → γ/(γ − 1)(p/ρ − p_o/ρ_o) with ρ = ρ_o(p/p_o)^{1/γ}; or a callable ρ(p) → ``quad``
    rho [kg/m³]; T [K]; R [J/(kg K)]; gamma [-]; rho_o [kg/m³]

    Returns P [J/kg].
    Validation: V1 the three closed forms equal ``quad`` of the matching callable ρ(p) (1e-10); isentropic P equals
    h − h_o = C_p(T − T_o) (link to (4.78) via ``core.thermo.perfect_gas_enthalpy``). Label: analytic.
    """
    p_, po = _F(p), float(p_o)
    if callable(kind):
        f = np.vectorize(lambda pp: quad(lambda s: 1.0 / float(kind(s)), po, float(pp), epsabs=0.0, epsrel=1e-12)[0])
        return _S(f(p_))
    if kind == "constant":
        return _S((p_ - po) / float(rho))
    if kind == "isothermal":
        return _S(float(R) * float(T) * np.log(p_ / po))
    if kind == "isentropic":
        r = float(rho_o) * (p_ / po) ** (1.0 / gamma)
        return _S(gamma / (gamma - 1.0) * (p_ / r - po / float(rho_o)))
    raise ValueError("kind must be 'constant', 'isothermal', 'isentropic' or a callable rho(p)")


def bernoulli_function(speed, p, z, rho=1000.0, g: float = G0, kind="constant", p_o: float = 0.0, **state):
    """Bernoulli function B = ½|u|² + ∫_{p_o}^p dp′/ρ + gz [J/kg], the bracket of Eq. (4.69).

    Book: §4.9, Eqs. (4.69)–(4.72): ∂u/∂t + ∇B = u × ω. Steady ⇒ ∇B = u × ω (4.70), so B is constant on Lamb
    surfaces (streamlines and vortex lines) (4.71), and everywhere if also irrotational (4.72).

    Parameters: speed |u| [m/s]; p [Pa]; z [m]; rho [kg/m³] (for kind="constant"); g [m/s²]; kind, **state: see
    :func:`pressure_function`; p_o : reference pressure [Pa]. Returns B [J/kg]; scalar-callable.
    Validation: V1 Rankine vortex: B constant on each circle, varies across circles inside the core, uniform outside;
    cylinder potential flow: uniform (4.72). Label: analytic.
    """
    P = pressure_function(p, p_o, kind, rho=rho, **state)
    return _S(0.5 * _F(speed) ** 2 + _F(P) + float(g) * _F(z))  # Eq. (4.69) bracket


def lamb_surface_terms(u: Callable, B: Callable, x, t: float = 0.0, h: float = 1e-4) -> dict:
    """(4.70) ∇B = u × ω and its consequences u·∇B = ω·∇B = 0 (B constant along streamlines and vortex lines — Lamb
    surfaces, Fig. 4.13). Returns dict(u_dot_gradB, omega_dot_gradB, residual = ∇B − u × ω (3, …)).
    Book: §4.9, Eqs. (4.70)–(4.71). Label: analytic."""
    x_ = _F(x)
    d = x_.shape[0]
    gB = grad(B, x_, t, h)
    U = ev(u, x_, t, (d,))
    om = vorticity(lambda X, T: ev(u, X, T, (d,)), x_, t, h)
    pad = lambda a: a if d == 3 else np.concatenate([a, np.zeros((1,) + a.shape[1:])])  # noqa: E731
    U3, gB3 = pad(U), pad(gB)
    lam = np.cross(U3, om, axis=0)
    return {"u_dot_gradB": _S(np.sum(U3 * gB3, axis=0)), "omega_dot_gradB": _S(np.sum(om * gB3, axis=0)),
            "residual": gB3 - lam}  # Eq. (4.70)


def lamb_surface_check(u: Callable, B: Callable, x, h: float = 1e-4, t: float = 0.0):
    """|∇B − u × ω| at the point(s) x — zero for steady, inviscid, barotropic flow, Eq. (4.70).

    Book: §4.9, Eqs. (4.70)–(4.71). Parameters: u(x, t) [m/s]; B(x, t) [J/kg]; x [m]; h [m]; t [s]. Returns [m/s²]
    (float for one point). See :func:`lamb_surface_terms` for the parts.
    Validation: V1 steady Rankine vortex (B from :func:`rankine_vortex_pressure`): ≈ 0 everywhere. Label: analytic.
    """
    r = lamb_surface_terms(u, B, x, t, h)["residual"]
    return _S(np.sqrt(np.sum(r ** 2, axis=0)))


def bernoulli_along_line(u: Callable, p: Callable, points, rho: float = 1000.0, g: float = G0, t: float = 0.0,
                         z_axis: int | None = -1) -> np.ndarray:
    """B = ½|u|² + p/ρ + gz at points (d, N) along a line (e.g. a ``core.kinematics.streamline``); z = the coordinate
    ``z_axis`` (default the last; None → no gravity term, a horizontal plane).

    Book: §4.9, Eqs. (4.19), (4.71) (constant density). Returns B [J/kg] (N,). Label: analytic.
    """
    X = _F(points)
    d = X.shape[0]
    U = ev(u, X, t, (d,))
    z = X[z_axis] if z_axis is not None else 0.0 * X[0]
    return 0.5 * np.sum(U ** 2, axis=0) + ev(p, X, t) / rho + g * z  # Eq. (4.19)


def rankine_vortex_pressure(r, Gamma: float, sigma: float, rho: float = 1000.0, p_inf: float = 0.0):
    """Pressure of a steady Rankine vortex (core radius σ) from the radial balance dp/dr = ρu_θ²/r, p → p∞ far away [Pa]:
    outside p = p∞ − ρΓ²/(8π²r²); inside p = p∞ − (ρΓ²/8π²σ²)(2 − r²/σ²) (continuous at r = σ).

    Book: §4.9 (the Bernoulli function of a rotational flow, E6's Rankine scenario; profile (3.28) of ch03). Our field.
    Validation: V1 dp/dr = ρu_θ²/r by central differences; continuity at σ. Label: analytic.
    """
    r_ = _F(r)
    k = rho * Gamma ** 2 / (8.0 * np.pi ** 2)
    with np.errstate(divide="ignore"):
        out = np.where(r_ >= sigma, p_inf - k / np.where(r_ > 0, r_, 1.0) ** 2,
                       p_inf - k / sigma ** 2 * (2.0 - r_ ** 2 / sigma ** 2))
    return _S(out)


# ======================================================================================================================
# unsteady and viscous irrotational forms (4.73)–(4.83)
# ======================================================================================================================
def unsteady_bernoulli_B(phi: Callable, p: Callable, x, t: float = 0.0, rho: float = 1000.0, g: float = G0,
                         h: float = 1e-4, ht: float | None = None, z_axis: int | None = None):
    """The bracket of (4.74) for a constant-density potential flow, ∂φ/∂t + ½|∇φ|² + p/ρ + gz, at the point(s) x — the
    same value B(t) everywhere.

    Book: §4.9, Eq. (4.74): ∇[∂φ/∂t + ½|∇φ|² + ∫dp/ρ + gz] = 0 ⇒ the bracket is B(t) only; the gauge φ → φ − ∫B dt′
    (book prints "+", a sign slip) makes it a constant, (4.75); (4.83) is the same with ρ constant.

    Parameters: phi(x, t) [m²/s]; p(x, t) [Pa]; x [m]; t [s]; rho [kg/m³]; g [m/s²]; h [m]; ht [s]; z_axis (None → no
    gravity). Returns B [J/kg] (float or array; uniform for a valid pair).
    Validation: V1 sphere accelerating in still fluid (``ch04.accelerating_sphere_fields``): B uniform. Label: analytic.
    """
    x_ = _F(x)
    gp = grad(phi, x_, t, h)
    z = x_[z_axis] if z_axis is not None else 0.0 * x_[0]
    return _S(ddt(phi, x_, t, ht) + 0.5 * np.sum(gp ** 2, axis=0) + ev(p, x_, t) / rho + g * z)  # Eq. (4.74)


def gauge_absorbed_bracket(B, sign: float = -1.0):
    """Value of the (4.74) bracket after the gauge change φ → φ + sign·∫B dt′: B + sign·B — 0 for sign = −1 (the correct
    change) and 2B for the book's printed "+" (sign = +1).

    Book: §4.9, (4.74)–(4.75) (analysis §9 item 3). Scalar-callable. Label: analytic.
    """
    return _S(_F(B) + float(sign) * _F(B))


def unsteady_bernoulli_pressure(dphi_dt, speed, z=0.0, rho: float = 1000.0, g: float = G0, C: float = 0.0):
    """Pressure from the unsteady Bernoulli equation (4.75)/(4.83): p = ρ[C − ∂φ/∂t − ½|∇φ|² − gz] [Pa].

    Book: §4.9, Eqs. (4.75), (4.83). Parameters: dphi_dt = ∂φ/∂t [m²/s²]; speed = |∇φ| [m/s] (or a vector (d, …));
    z [m]; rho [kg/m³]; g; C [J/kg] (the constant). Scalar-callable. Label: analytic.
    """
    gp = _F(speed)
    q2 = gp ** 2 if gp.ndim == 0 else np.sum(gp ** 2, axis=0)
    return _S(rho * (C - _F(dphi_dt) - 0.5 * q2 - g * _F(z)))  # Eq. (4.83)


def viscous_irrotational_residual(u: Callable, x, t: float = 0.0, mu: float = 1.0e-3, h: float = 1e-4) -> np.ndarray:
    """The viscous force −μ∇×ω = μ∇²u of an incompressible flow (4.40) — exactly zero for irrotational flow, which is why
    (4.80)–(4.83) hold even with viscosity. Returns (d,) + points [N/m³].

    Book: §4.9, Eqs. (4.79)–(4.80). Validation: V1 any potential flow (cylinder, source, vortex) gives ≈ 0; a shear flow
    with curvature does not. Label: analytic.
    """
    x_ = _F(x)
    d = x_.shape[0]
    return float(mu) * laplacian(lambda X, T: ev(u, X, T, (d,)), x_, t, h, (d,))


def unsteady_streamline_bernoulli(dudt_along, s, state1: dict, state2: dict, rho: float = 1000.0,
                                  g: float = G0) -> float:
    """Residual of the unsteady streamline Bernoulli equation for constant-ρ, constant-μ irrotational flow, (4.82):
    ∫₁²(∂u/∂t)·ds + (½|u|² + gz + p/ρ)₂ − (½|u|² + gz + p/ρ)₁ [J/kg] (0 when the equation holds).

    Book: §4.9, Eqs. (4.80)–(4.82) (points 1, 2 on one streamline at a single instant).
    Parameters: dudt_along : samples of e_s·∂u/∂t along the streamline [m/s²]; s : arc length at the samples [m];
    state1, state2 : dicts U [m/s], z [m], p [Pa]; rho; g. Returns the residual (trapezoid rule for the integral).
    Validation: V1 U-tube column (``ch04.u_tube_column``) gives 0; a straight pipe p₁ − p₂ = ρL dU/dt. Label: analytic.
    """
    I = float(np.trapezoid(_F(dudt_along), _F(s)))
    B1 = bernoulli_head(state1["U"], state1["z"], state1["p"], rho, g)
    B2 = bernoulli_head(state2["U"], state2["z"], state2["p"], rho, g)
    return I + B2 - B1  # Eq. (4.82)


# ======================================================================================================================
# stagnation quantities, pitot, orifice
# ======================================================================================================================
def stagnation_enthalpy(h, U, z=0.0, g: float = G0):
    """h + ½|u|² + gz [J/kg], constant on streamlines of steady inviscid non-conducting flow, Eq. (4.78).

    Book: §4.9, Eqs. (4.76)–(4.78) (from the energy equation (4.55) with σ = q = 0 and steady continuity) — an energy
    statement, not the momentum Bernoulli (4.71). Label: analytic.
    """
    return _S(_F(h) + 0.5 * _F(U) ** 2 + g * _F(z))  # Eq. (4.78)


def stagnation_temperature(T, U, cp: float = CP_BOOK):
    """T₀ = T + U²/(2C_p) [K] — (4.78) with h = C_pT at constant z (a perfect gas slowed to rest adiabatically).

    Book: §4.9, Eq. (4.78). Validation: V1 T₀/T = 1 + (γ − 1)M²/2 with M = U/√(γRT) (ch01 isentropic relations);
    ``stagnation_temperature(300.0, 100.0)`` = 304.98 K. Label: analytic.
    """
    return _S(_F(T) + _F(U) ** 2 / (2.0 * cp))


def stagnation_pressure(p, U, rho):
    """Stagnation (total) pressure p + ½ρ|u|² [Pa] (text after Fig. 4.16). Book: §4.9. Label: analytic."""
    return _S(_F(p) + 0.5 * _F(rho) * _F(U) ** 2)


def dynamic_pressure(U, rho):
    """Dynamic pressure ½ρu² [Pa] (text after Fig. 4.16). Book: §4.9. Label: analytic."""
    return _S(0.5 * _F(rho) * _F(U) ** 2)


def pitot_speed(p_stag, p_static, rho):
    """Flow speed from a pitot-static pair, |u|₁ = √(2(p₂ − p₁)/ρ) [m/s] (Fig. 4.15; (4.19) between a point in the stream
    and the stagnation point on the same streamline). Book: §4.9. Label: analytic."""
    dp = _F(p_stag) - _F(p_static)
    if np.any(dp < 0):
        raise ValueError("stagnation pressure below static pressure")
    return _S(np.sqrt(2.0 * dp / _F(rho)))


def pitot_speed_from_heads(h1, h2, g: float = G0, rho: float = 1000.0, rho_atm: float = 0.0):
    """Pitot speed from the manometer heads, |u|₁ = √(2g(h₂ − h₁)(1 − ρ_atm/ρ)) [m/s] (Fig. 4.15).

    Book: §4.9, text after Fig. 4.15: |u|₁ = √(2g(h₂ − h₁)) when ρ_atm ≪ ρ; the tops of the columns differ in pressure by
    ρ_atm g(h₂ − h₁), which gives the (1 − ρ_atm/ρ) correction (our step). Label: analytic.
    """
    return _S(np.sqrt(2.0 * g * (_F(h2) - _F(h1)) * (1.0 - float(rho_atm) / float(rho))))


def torricelli_speed(h, g: float = G0):
    """Jet speed from a tank orifice at depth h below the free surface, u = √(2gh) [m/s] (text after Fig. 4.17;
    (4.19) from the free surface to section C where p = p_atm). Book: §4.9. Label: analytic."""
    return _S(np.sqrt(2.0 * g * _F(h)))


# ======================================================================================================================
# which Bernoulli? (the summary after (4.83)); E6 scenarios
# ======================================================================================================================
BERNOULLI_FORMS: dict[str, dict] = {
    "4.19": {"tex": r"\tfrac12U^2 + gz + p/\rho = \text{const}", "constant_along": "a streamline",
             "hypotheses": dict(steady=True, viscous=False, constant_density=True)},
    "4.71": {"tex": r"\tfrac12u_i^2 + \int_{p_o}^{p}\frac{dp'}{\rho(p')} + gz = \text{const}",
             "constant_along": "streamlines and vortex lines",
             "hypotheses": dict(steady=True, viscous=False, barotropic=True)},
    "4.72": {"tex": r"\tfrac12u_i^2 + \int_{p_o}^{p}\frac{dp'}{\rho(p')} + gz = \text{const}",
             "constant_along": "everywhere",
             "hypotheses": dict(steady=True, viscous=False, barotropic=True, irrotational=True)},
    "4.75": {"tex": r"\frac{\partial\phi}{\partial t} + \tfrac12|\nabla\phi|^2 + \int_{p_o}^{p}\frac{dp'}{\rho(p')} + gz"
                    r" = \text{const}",
             "constant_along": "everywhere at one instant",
             "hypotheses": dict(viscous=False, barotropic=True, irrotational=True)},
    "4.78": {"tex": r"h + \tfrac12|\mathbf u|^2 + gz = \text{const}", "constant_along": "streamlines",
             "hypotheses": dict(steady=True, viscous=False, isentropic=True)},
    "4.82": {"tex": r"\int_1^2\frac{\partial\mathbf u}{\partial t}\cdot d\mathbf s"
                    r" + \left(\tfrac12|\mathbf u|^2 + gz + \frac{p}{\rho}\right)_2"
                    r" = \left(\tfrac12|\mathbf u|^2 + gz + \frac{p}{\rho}\right)_1",
             "constant_along": "a streamline at one instant (and (4.83) everywhere)",
             "hypotheses": dict(viscous=True, irrotational=True, constant_density=True)},
}


def which_bernoulli(steady: bool, viscous: bool, irrotational: bool, barotropic: bool = False,
                    isentropic: bool = False, constant_density: bool = False) -> list[str]:
    """The Bernoulli equations valid for a flow's hypotheses, in the order (4.19), (4.71), (4.72), (4.75), (4.78), (4.82).

    Book: §4.9, the summary after Eq. (4.83). Rules: constant density implies barotropic; (4.19), (4.71), (4.72),
    (4.75), (4.78) need inviscid flow ((4.78) also isentropic: no heat conduction); (4.82)/(4.83) are *the* Bernoulli
    equation of viscous flow — constant μ and ρ and irrotational (μ∇²u = −μ∇×ω = 0, (4.80)) — so they are listed for
    viscous flows (an inviscid irrotational constant-ρ flow is covered by (4.75)). Inviscid + rotational + unsteady → [].
    Validation: V1 ``which_bernoulli(True, False, True, True, constant_density=True)`` = ["4.19", "4.71", "4.72",
    "4.75"]; without constant density ["4.71", "4.72", "4.75"]. Label: analytic.
    """
    have = dict(steady=bool(steady), viscous=bool(viscous), irrotational=bool(irrotational),
                barotropic=bool(barotropic) or bool(constant_density), isentropic=bool(isentropic),
                constant_density=bool(constant_density))
    return [key for key, form in BERNOULLI_FORMS.items()
            if all(have[k] == v if k == "viscous" else (have[k] or not v) for k, v in form["hypotheses"].items())]


def which_bernoulli_text(steady: bool, viscous: bool, irrotational: bool, barotropic: bool = False,
                         isentropic: bool = False, constant_density: bool = False) -> str:
    """:func:`which_bernoulli` as one string, the labels joined by commas (e.g. "4.19,4.71,4.72,4.75"; "" when none
    applies) — for E6's exact-text parity rows. Book: §4.9, the summary after Eq. (4.83). Label: analytic."""
    return ",".join(which_bernoulli(steady, viscous, irrotational, barotropic, isentropic, constant_density))


BERNOULLI_SCENARIOS = ("pitot", "orifice", "rankine", "cylinder", "u_tube", "hot_nozzle")


def _status(valid: list[str], B_along, B_across) -> str:
    flat = lambda a: float(np.ptp(a)) <= 1e-9 * max(float(np.max(np.abs(a))), 1.0)  # noqa: E731
    if not valid:
        return "no Bernoulli equation applies"
    if "4.72" in valid and flat(B_across):
        return "(4.72): constant everywhere"
    if "4.71" in valid or "4.19" in valid:
        if flat(B_along) and not flat(B_across):
            return "rotational: constant along each streamline, differs across them"
        return "(4.71): constant along this streamline"
    if "4.82" in valid or "4.75" in valid:
        return "unsteady: the ∂φ/∂t term carries the difference ((4.75)/(4.82))"
    if "4.78" in valid:
        return "(4.78): h + ½u² constant on the streamline"
    return ", ".join(valid)


def bernoulli_scenario(name: str, **p) -> dict:
    """E6's six scenarios: hypotheses, valid Bernoulli forms, B along and across streamlines, terms at the probe.

    Book: §4.9 (Figs. 4.13–4.17 and the Bernoulli equations). Our numbers (not the book's):

    * "pitot" (U = 3 m/s, rho, p_static): free stream → stagnation point on one streamline ((4.19));
    * "orifice" (h = 1 m, Cc, A): free surface → vena contracta (Torricelli);
    * "rankine" (Gamma = 1 m²/s, sigma = 0.1 m, rho, r_along = σ/2): B on a circle (flat) and across radii (varies inside
      the core, flat outside) — rotational: (4.71) only;
    * "cylinder" (U, a): ideal flow, B at random points (flat everywhere, (4.72));
    * "u_tube" (L = 1 m, h0 = 0.05 m, t): the two free surfaces — B differs; L dU/dt closes (4.82);
    * "hot_nozzle" (T0 = 600 K, U = 300 m/s, cp): h + ½U² along the nozzle ((4.78)).

    Returns dict(name, hypotheses, valid (labels), status (text), s_along, B_along, s_across, B_across [m²/s²],
    terms (½u², p/ρ, gz, ∂φ/∂t at the probe), numbers (scenario results)).
    Validation: V1 flatness (ptp/scale < 1e-10) where the valid form says so. Label: analytic.
    """
    g = float(p.get("g", G0))
    rho = float(p.get("rho", 1000.0))
    nums: dict = {}
    if name == "pitot":
        U, ps = float(p.get("U", 3.0)), float(p.get("p_static", 101325.0))
        hyp = dict(steady=True, viscous=False, irrotational=True, constant_density=True)
        p0 = float(stagnation_pressure(ps, U, rho))
        s_al = np.array([0.0, 1.0])
        B_al = np.array([bernoulli_head(U, 0.0, ps, rho, g), bernoulli_head(0.0, 0.0, p0, rho, g)])
        s_ac, B_ac = s_al.copy(), B_al.copy()
        terms = dict(ke=0.5 * U ** 2, p_over_rho=ps / rho, gz=0.0, dphi_dt=0.0)
        nums = dict(p_stagnation=p0, dh=(p0 - ps) / (rho * g), U_measured=float(pitot_speed(p0, ps, rho)))
    elif name == "orifice":
        h, Cc, A = float(p.get("h", 1.0)), float(p.get("Cc", 0.611)), float(p.get("A", 1e-3))
        hyp = dict(steady=True, viscous=False, irrotational=True, constant_density=True)
        uj = float(torricelli_speed(h, g))
        s_al = np.array([0.0, 1.0])
        B_al = np.array([bernoulli_head(0.0, h, 0.0, rho, g), bernoulli_head(uj, 0.0, 0.0, rho, g)])
        s_ac, B_ac = s_al.copy(), B_al.copy()
        terms = dict(ke=0.5 * uj ** 2, p_over_rho=0.0, gz=0.0, dphi_dt=0.0)
        nums = dict(jet_speed=uj, mass_flow=rho * Cc * A * uj, jet_area=Cc * A)
    elif name == "rankine":
        Gm, a = float(p.get("Gamma", 1.0)), float(p.get("sigma", 0.1))
        r_al = float(p.get("r_along", 0.5 * a))
        hyp = dict(steady=True, viscous=False, irrotational=False, constant_density=True)

        def ut(r):
            r_ = _F(r)
            return np.where(r_ <= a, Gm * r_ / (2 * np.pi * a ** 2), Gm / (2 * np.pi * np.maximum(r_, 1e-300)))

        def B_at(r):
            return 0.5 * ut(r) ** 2 + _F(rankine_vortex_pressure(r, Gm, a, rho)) / rho
        s_al = np.linspace(0.0, 2 * np.pi * r_al, 9)[:-1]
        B_al = B_at(np.full(s_al.shape, r_al))
        s_ac = np.linspace(0.1 * a, 3.0 * a, 12)
        B_ac = B_at(s_ac)
        terms = dict(ke=float(0.5 * ut(r_al) ** 2), p_over_rho=float(rankine_vortex_pressure(r_al, Gm, a, rho)) / rho,
                     gz=0.0, dphi_dt=0.0)
        nums = dict(B_inside_span=float(np.ptp(B_ac[s_ac <= a])), B_outside_span=float(np.ptp(B_ac[s_ac > a])))
    elif name == "cylinder":
        U, a = float(p.get("U", 1.0)), float(p.get("a", 1.0))
        hyp = dict(steady=True, viscous=False, irrotational=True, constant_density=True)
        rng = np.random.default_rng(int(p.get("seed", 0)))
        r = a * (1.0 + 2.0 * rng.random(24))
        th = 2 * np.pi * rng.random(24)
        X, Y = r * np.cos(th), r * np.sin(th)
        r2 = X ** 2 + Y ** 2
        uu = U * (1 - a ** 2 * (X ** 2 - Y ** 2) / r2 ** 2)
        vv = -2 * U * a ** 2 * X * Y / r2 ** 2
        pp = 0.5 * rho * (U ** 2 - uu ** 2 - vv ** 2)
        B = 0.5 * (uu ** 2 + vv ** 2) + pp / rho
        s_al, B_al = np.arange(8.0), B[:8]
        s_ac, B_ac = np.arange(24.0), B
        terms = dict(ke=float(0.5 * (uu[0] ** 2 + vv[0] ** 2)), p_over_rho=float(pp[0] / rho), gz=0.0, dphi_dt=0.0)
        th_s = np.linspace(0.0, np.pi, 7)
        nums = dict(theta_surface=th_s, Cp_surface=1.0 - 4.0 * np.sin(th_s) ** 2)  # C_p on r = a from (4.72)
    elif name == "u_tube":
        L, h0, t = float(p.get("L", 1.0)), float(p.get("h0", 0.05)), float(p.get("t", 0.3))
        hyp = dict(steady=False, viscous=True, irrotational=True, constant_density=True)
        w = np.sqrt(2.0 * g / L)
        x = h0 * np.cos(w * t)
        U = -h0 * w * np.sin(w * t)
        dUdt = -h0 * w ** 2 * np.cos(w * t)
        s_al = np.array([0.0, L])
        B_al = np.array([bernoulli_head(U, -x, 0.0, rho, g), bernoulli_head(U, x, 0.0, rho, g)])  # the two surfaces
        s_ac, B_ac = s_al.copy(), B_al.copy()
        terms = dict(ke=0.5 * U ** 2, p_over_rho=0.0, gz=g * x, dphi_dt=L * dUdt)
        nums = dict(omega=float(w), x=float(x), U=float(U), dUdt=float(dUdt), dp_over_L=float(rho * dUdt),
                    residual_4_82=float(L * dUdt + B_al[1] - B_al[0]))
    elif name == "hot_nozzle":
        T0, U, cp = float(p.get("T0", 600.0)), float(p.get("U", 300.0)), float(p.get("cp", CP_BOOK))
        hyp = dict(steady=True, viscous=False, irrotational=True, isentropic=True)
        T = T0 - U ** 2 / (2 * cp)
        s_al = np.array([0.0, 1.0])
        B_al = np.array([cp * T0, cp * T + 0.5 * U ** 2])
        s_ac, B_ac = s_al.copy(), B_al.copy()
        terms = dict(ke=0.5 * U ** 2, p_over_rho=0.0, gz=0.0, dphi_dt=0.0, h=cp * T)
        nums = dict(T=T, T0=T0)
    else:
        raise ValueError(f"unknown scenario {name!r}; choose from {BERNOULLI_SCENARIOS}")
    valid = which_bernoulli(hyp.get("steady", False), hyp.get("viscous", False), hyp.get("irrotational", False),
                            barotropic=hyp.get("constant_density", False), isentropic=hyp.get("isentropic", False),
                            constant_density=hyp.get("constant_density", False))
    return dict(name=name, hypotheses=hyp, valid=valid, status=_status(valid, _F(B_al), _F(B_ac)),
                s_along=_F(s_al), B_along=_F(B_al), s_across=_F(s_ac), B_across=_F(B_ac), terms=terms, numbers=nums)
