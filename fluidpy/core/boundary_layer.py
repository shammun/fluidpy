"""Laminar boundary layers (Ch. 9): scales, thicknesses, Blasius and Falkner–Skan similarity solutions, the von Kármán
momentum integral, Thwaites' method, separation diagnostics, a parabolic marching solver, transition and plate drag.

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e, Ch. 9, §9.1–§9.7, Eqs. (9.1)–(9.52), Examples 9.1–9.2.  Every equation
was transcribed from the rendered page images (analysis/ch09.md).  Conventions (analysis §9): x along the wall, y normal;
u = ∂ψ/∂y, v = −∂ψ/∂x; U_e(x) edge speed; dp/dx > 0 is adverse; ``phi`` on a body is measured from the FORWARD stagnation
point; Re_x = U_e x/ν.  All functions are scalar-callable (plain floats in, plain floats out) unless they return arrays
by nature, and use ``np.trapezoid`` (numpy 2).

Printed slips coded corrected (with the printed form available as a named option): R1 (9.7) squares
(``ch09.bl_nondim_sympy(printed_9_7=True)``); R2 (9.30) η₉₉ = 4.93 vs the root 4.910 (``blasius_delta99(printed=True)``);
R11 there is no attached bounded (0 ≤ f′ ≤ 1) solution below the fold; reverse-flow members live on the second branch
(``falkner_skan(branch="reversed")``).  Trap/note (NOT a printed slip): (9.33) is for ONE side of the plate, as the book says (``sides=``).

Numerics: Blasius is solved by Töpfer scaling (one initial-value problem, exact to solver tolerance) AND by
``solve_bvp``; Falkner–Skan by ``solve_bvp`` with continuation in m; the fold at m = −0.0904 is crossed by
parametrising by f″(0).  No long shooting.
"""
from __future__ import annotations

import functools
from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.integrate import cumulative_trapezoid, quad, simpson, solve_bvp, solve_ivp
from scipy.interpolate import CubicSpline, PchipInterpolator
from scipy.linalg import solve_banded
from scipy.optimize import brentq
from scipy.special import erfc

from ._util import as_scalar_if_0d

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d

LAMBDA_SEP_BOOK = -0.09
"""Thwaites' separation value of λ = (θ²/ν)dU_e/dx (Table 9.1: l(−0.09) = 0).  Label: book value (rounded)."""

# ======================================================================================================================
# §9.1 scales and variables
# ======================================================================================================================


def boundary_layer_scales(U, L, nu, rho: float = 1.2) -> dict:
    """Order-of-magnitude sizes of the boundary-layer terms and the thickness δ̄ ~ L/√Re.

    Book: §9.1, Eqs. (9.2)–(9.4): u∂u/∂x ~ U²/L (9.2); ν∂²u/∂y² ~ νU/δ̄² (9.3); equating gives δ̄ ~ √(νL/U), δ̄/L ~ Re^{-1/2}
    (9.4); continuity gives v ~ δ̄U/L = U Re^{-1/2}; the wall stress τ₀ ~ μU/δ̄ gives C_f ~ 2/√Re (unnumbered, §9.1).

    Parameters
    ----------
    U : free-stream speed [m/s].  L : streamwise length [m].  nu : kinematic viscosity [m²/s].
    rho : density [kg/m³] (only for the wall-stress scale in Pa; default air, 1.2).

    Returns
    -------
    dict: Re [–] = UL/ν; delta_over_L [–]; delta [m]; v_scale [m/s]; tau0_scale [Pa] = ρνU/δ̄ = μU/δ̄; cf_estimate [–] = 2/√Re;
    adv [m/s²] = U²/L; visc [m/s²] = νU/δ̄² (equal to adv by construction); visc_x [m/s²] = νU/L² (the streamwise viscous term
    that (9.7) drops).

    Assumptions: order of magnitude only (the factor differs from flow to flow: Blasius has 0.664 not 2).
    Validation: V1 identities δ̄ = L/√Re, adv = visc; parity with core.laminar.diffusion_thickness scaling √(νt), t = L/U.
    Label: analytic.
    """
    U, L, nu = float(U), float(L), float(nu)
    Re = U * L / nu
    delta = np.sqrt(nu * L / U)  # Eq. (9.4)
    return dict(Re=Re, delta_over_L=delta / L, delta=delta, v_scale=delta * U / L,
                adv=U ** 2 / L,  # Eq. (9.2)
                visc=nu * U / delta ** 2,  # Eq. (9.3)
                visc_x=nu * U / L ** 2,
                tau0_scale=rho * nu * U / delta, cf_estimate=2.0 / np.sqrt(Re))


def _re_from(Re, U, L, nu):
    """Re = UL/ν from either ``Re`` or ``nu`` (one of them is required)."""
    if Re is not None:
        return float(Re)
    if nu is None:
        raise ValueError("give Re = UL/nu or nu (Re is not otherwise known)")
    return float(U) * float(L) / float(nu)


def to_bl_variables(x, y, u, v, p, L, U, rho, Re=None, p_inf: float = 0.0, nu=None):
    """Dimensional (x, y, u, v, p) → the stretched, dimensionless boundary-layer variables of Eq. (9.6).

    Book: §9.1, Eq. (9.6): x* = x/L, y* = (y/L)Re^{1/2}, u* = u/U, v* = (v/U)Re^{1/2}, p* = (p − p∞)/(ρU²).
    Parameters: x, y [m]; u, v [m/s]; p [Pa]; L [m]; U [m/s]; rho [kg/m³]; Re = UL/ν [–] (``None`` → computed from ``nu`` [m²/s],
    which is then required); p_inf [Pa].
    Returns (x*, y*, u*, v*, p*) [–].  Label: analytic (round trip with :func:`from_bl_variables`).
    """
    sre = np.sqrt(_re_from(Re, U, L, nu))
    return (_S(_F(x) / L), _S(_F(y) / L * sre), _S(_F(u) / U), _S(_F(v) / U * sre), _S((_F(p) - p_inf) / (rho * U ** 2)))


def from_bl_variables(xs, ys, us, vs, ps, L, U, rho, Re=None, p_inf: float = 0.0, nu=None):
    """Inverse of :func:`to_bl_variables` (Eq. (9.6)); returns dimensional (x, y, u, v, p) [SI]."""
    sre = np.sqrt(_re_from(Re, U, L, nu))
    return (_S(_F(xs) * L), _S(_F(ys) * L / sre), _S(_F(us) * U), _S(_F(vs) * U / sre), _S(_F(ps) * rho * U ** 2 + p_inf))


# ======================================================================================================================
# outer flow, matching (9.11)
# ======================================================================================================================


@dataclass
class OuterFlow:
    """An edge velocity U_e(x) with its derivative and the pressure gradient it imposes, Eq. (9.11).

    Attributes: ``kind`` (str), ``params`` (dict), ``Ue(x)`` [m/s], ``dUe(x)`` [1/s], ``dpdx(x, rho)`` [Pa/m] = −ρU_e dU_e/dx
    (Eq. (9.11); dp/dx > 0 adverse), ``x_stag`` (location where U_e = 0, else None).
    """

    kind: str
    params: dict
    Ue: Callable
    dUe: Callable
    x_stag: float | None = None

    def dpdx(self, x, rho: float = 1.0):
        """dp/dx = −ρ U_e dU_e/dx  [Pa/m], Eq. (9.11)."""
        return -rho * self.Ue(x) * self.dUe(x)  # Eq. (9.11)


def outer_flow(kind: str = "flat", **params) -> OuterFlow:
    """Edge-velocity distributions used in Ch. 9 with the pressure gradient of Eq. (9.11).

    Book: §9.1 Eq. (9.11) −(1/ρ)dp/dx = U_e dU_e/dx; §9.3 flat plate; §9.4 Eq. (9.35) U_e = a xⁿ; §9.6 Example 9.2 diffuser
    (U_e = U₁/(1 + x/L)); cylinder U_e = 2U sin(x/a) (ideal flow (6.x), ch06, x = aφ from the forward stagnation point).

    kinds and parameters (SI): "flat" (U); "wedge" (n, a) U_e = a xⁿ [a in m^{1−n}/s]; "stagnation" (a) = wedge n = 1;
    "diffuser" (U1, L); "cylinder" (U, a); "linear_retarded"/"retarded" (U0, L, c=1) U_e = U₀(1 − c x/L) (Howarth);
    "custom" (Ue callable, dUe optional callable, else a central difference).
    Returns an :class:`OuterFlow`.  Validation: V1 cylinder parity with core.potential.cylinder surface speed (1e-13);
    pressure-gradient identity.  Label: analytic.
    """
    k = kind
    if k == "flat":
        U = float(params.get("U", 1.0))
        return OuterFlow(k, dict(U=U), lambda x: U + 0.0 * _F(x), lambda x: 0.0 * _F(x))
    if k in ("wedge", "stagnation"):
        n = 1.0 if k == "stagnation" else float(params.get("n", 0.0))
        a = float(params.get("a", 1.0))
        return OuterFlow(k, dict(n=n, a=a), lambda x: a * _F(x) ** n,
                         lambda x: a * n * _F(x) ** (n - 1.0) if n != 0 else 0.0 * _F(x),
                         0.0 if n > 0 else None)
    if k == "diffuser":
        U1, L = float(params.get("U1", 1.0)), float(params.get("L", 1.0))
        return OuterFlow(k, dict(U1=U1, L=L), lambda x: U1 / (1.0 + _F(x) / L), lambda x: -U1 / (L * (1.0 + _F(x) / L) ** 2))
    if k == "cylinder":
        U, a = float(params.get("U", 1.0)), float(params.get("a", 1.0))
        return OuterFlow(k, dict(U=U, a=a), lambda x: 2.0 * U * np.sin(_F(x) / a), lambda x: 2.0 * U / a * np.cos(_F(x) / a), 0.0)
    if k in ("linear_retarded", "retarded"):
        U0, L, c = float(params.get("U0", 1.0)), float(params.get("L", 1.0)), float(params.get("c", 1.0))
        return OuterFlow(k, dict(U0=U0, L=L, c=c), lambda x: U0 * (1.0 - c * _F(x) / L), lambda x: -c * U0 / L + 0.0 * _F(x))
    if k == "custom":
        Ue = params["Ue"]
        dUe = params.get("dUe")
        if dUe is None:
            def dUe(x):  # central difference
                x = _F(x)
                h = 1e-6 * (1.0 + np.abs(x))
                return (Ue(x + h) - Ue(x - h)) / (2 * h)
        return OuterFlow(k, dict(), lambda x: _F(Ue(_F(x))), lambda x: _F(dUe(_F(x))), params.get("x_stag"))
    raise ValueError('kind must be "flat", "wedge", "stagnation", "diffuser", "cylinder", "linear_retarded" or "custom"')


def _as_outer(Ue, x=None) -> OuterFlow:
    """Accept an OuterFlow, a callable, or an array sampled on ``x``."""
    if isinstance(Ue, OuterFlow):
        return Ue
    if callable(Ue):
        return outer_flow("custom", Ue=Ue)
    if x is None:
        raise ValueError("an array U_e needs the x grid")
    cs = CubicSpline(_F(x), _F(Ue))
    return OuterFlow("custom", dict(), lambda s: cs(_F(s)), lambda s: cs(_F(s), 1))


# ======================================================================================================================
# residual checks of (9.9)–(9.10)
# ======================================================================================================================


def bl_x_momentum_residual(x, y, u, v, dpdx, nu, rho: float = 1.2, dudx=None):
    """Residual of the boundary-layer momentum equation u u_x + v u_y + (1/ρ)dp/dx − ν u_yy  [m/s²].

    Book: §9.1, Eq. (9.9) u∂u/∂x + v∂u/∂y = −(1/ρ)dp/dx + ν∂²u/∂y² (pressure gradient from Eq. (9.11)).
    Parameters: x [m] (1-D, length nx, or the meshgrid); y [m] (1-D, length ny, or the meshgrid); u, v [m/s] arrays u[j, i] = u(y_j, x_i)
    of shape (ny, nx) (``x`` along the last axis, ``y`` along the second-to-last); dpdx [Pa/m] scalar, (nx,) or (ny, nx); nu [m²/s];
    rho [kg/m³] (default air 1.2).  For a single station pass 1-D u(y), v(y) and ``dudx`` [1/s] (the streamwise derivative there);
    ``x`` is then ignored.
    Returns the residual [m/s²], same shape as u.  Derivatives: second-order central differences, ``np.gradient(edge_order=2)``.
    Validation: O(h²) small on the Blasius / Falkner–Skan fields.  Label: analytic (residual test).
    """
    u, v = _F(u), _F(v)
    y = _F(y)
    if y.ndim == 2:
        y = y[:, 0]
    dpdx = _F(dpdx)
    if u.ndim == 1:
        if dudx is None:
            raise ValueError("1-D input needs dudx")
        uy = np.gradient(u, y, edge_order=2)
        uyy = np.gradient(uy, y, edge_order=2)
        return u * _F(dudx) + v * uy + dpdx / rho - nu * uyy
    x = _F(x)
    if x.ndim == 2:
        x = x[0, :]
    uy = np.gradient(u, y, axis=0, edge_order=2)
    uyy = np.gradient(uy, y, axis=0, edge_order=2)
    ux = np.gradient(u, x, axis=1, edge_order=2)
    if dpdx.ndim == 1:
        dpdx = dpdx[None, :]
    return u * ux + v * uy + dpdx / rho - nu * uyy  # Eq. (9.9)


def bl_dpdy_scaled(u, v, v_x, v_y, v_xx, v_yy, Re):
    """∂p*/∂y* of the scaled y-momentum equation (9.8), solved for the pressure gradient [non-dimensional, in the stretched variables of Eq. (9.6)].

    Book: §9.1, Eq. (9.8): (1/Re)(u*v*_x* + v*v*_y*) = −∂p*/∂y* + (1/Re²)v*_x*x* + (1/Re)v*_y*y*, hence
    ∂p*/∂y* = −(1/Re)(u*v*_x* + v*v*_y*) + (1/Re²)v*_x*x* + (1/Re)v*_y*y*.
    Parameters (all in the starred variables of Eq. (9.6), non-dimensional): u, v velocities; v_x, v_y first derivatives of v; v_xx, v_yy second derivatives; Re = UL/ν [–].
    Returns ∂p*/∂y* [–].  Validation: V1 sign and magnitude on manufactured fields.  Label: analytic."""
    return -(_F(u) * _F(v_x) + _F(v) * _F(v_y)) / Re + _F(v_xx) / Re ** 2 + _F(v_yy) / Re  # Eq. (9.8)


def _bl_pressure_ratio(Re: float, n_x: int, n_y: int, y_max: float = 12.0):
    """(max|p*_y*|, max|advective term|, ∫p*_y* dy* at the last station) on the Blasius field, all in units of ρU² per L (see bl_pressure_variation)."""
    bl = _blasius()
    ys = np.linspace(0.0, y_max, int(n_y))
    xs_all = np.linspace(0.25, 1.0, int(n_x))
    h = 1e-4
    dmax = amax = 0.0
    delta_p = 0.0
    for xs in xs_all:
        def vfun(x):
            e = ys / np.sqrt(x)
            return (e * bl.fp(e) - bl.f(e)) / (2.0 * np.sqrt(x))

        eta = ys / np.sqrt(xs)
        f, fp, fpp = bl.f(eta), bl.fp(eta), bl.fpp(eta)
        u, v = fp, (eta * fp - f) / (2.0 * np.sqrt(xs))
        fppp = -0.5 * f * fpp
        vy = eta * fpp / (2.0 * xs)
        vyy = (fpp + eta * fppp) / (2.0 * xs ** 1.5)
        v_x = (vfun(xs + h) - vfun(xs - h)) / (2 * h)
        v_xx = (vfun(xs + h) - 2 * vfun(xs) + vfun(xs - h)) / h ** 2
        dpdy = bl_dpdy_scaled(u, v, v_x, vy, v_xx, vyy, Re)  # Eq. (9.8) solved for p*_y*
        adv = u * (-eta / (2 * xs)) * fpp + v * fpp / np.sqrt(xs)  # u*u*_x* + v*u*_y* (the terms p*_x* balances in (9.7))
        dmax = max(dmax, float(np.max(np.abs(dpdy))))
        amax = max(amax, float(np.max(np.abs(adv))))
        delta_p = float(np.trapezoid(dpdy, ys))
    return dmax, amax, delta_p


def bl_pressure_variation(Re, n_x: int = 41, n_y: int = 121) -> dict:
    """How much the pressure really changes across a Blasius layer — the O(1/Re) content of Eq. (9.8) ⇒ (9.10).

    Book: §9.1, Eqs. (9.8), (9.10).  Evaluates the right-hand side of the scaled y-momentum equation (9.8),
    ∂p*/∂y* = −(1/Re)(u*v*_x* + v*v*_y*) + (1/Re²)v*_x*x* + (1/Re)v*_y*y*  [(9.8) reads (1/Re)(u*v*_x* + v*v*_y*) = −p*_y* + (1/Re²)v*_x*x* + (1/Re)v*_y*y*], on the Blasius similarity fields at ``n_x`` stations x* ∈ [0.25, 1]
    (x* = x/L, y* = y Re^{1/2}/L, y* ∈ [0, 12] on ``n_y`` points) and compares it with the size of the streamwise terms that ∂p*/∂x* balances in (9.7).
    Parameters: Re = UL/ν [–] ≫ 1; n_x, n_y grid sizes (numerical choices of ours; the answer is insensitive to them).
    Returns dict(Re, ratio [–] = max|∂p*/∂y*| / max|u*u*_x* + v*u*_y*|, dpdy_max [–], adv_max [–], delta_p [–, pressure change across
    the layer in units of ρU²], order [–] = the exponent p in ratio ∝ Re^{−p}, measured between Re and 10Re; ratio_to_advection = ratio).
    DEVIATION: on an exact Blasius field ∂p*/∂x* is zero at leading order (dp/dx = 0), so the denominator is the size of the terms
    u*u*_x* + v*u*_y* that a non-zero ∂p*/∂x* would have to balance in (9.7) — an O(1) scale — rather than the (vanishing) pressure gradient itself.
    Validation: ratio ∝ Re^{−1}, order ≈ 1.  Label: analytic (order of magnitude).
    """
    Re = float(Re)
    d, a, dp = _bl_pressure_ratio(Re, n_x, n_y)
    d10, a10, _ = _bl_pressure_ratio(10.0 * Re, max(int(n_x) // 4, 5), n_y)
    ratio = d / a
    order = float(np.log((d / a) / (d10 / a10)) / np.log(10.0))
    return dict(Re=Re, ratio=ratio, dpdy_max=d, adv_max=a, delta_p=dp, order=order, ratio_to_advection=ratio)


# ======================================================================================================================
# §9.2 thicknesses
# ======================================================================================================================


def _deficit_tail(y, d):
    """Exponential-tail estimate of ∫_{y[-1]}^∞ d dy for a decaying deficit d (0 if d[-1] ≈ 0)."""
    if d[-1] <= 1e-14 or len(d) < 2 or d[-2] <= d[-1]:
        return 0.0
    L = (y[-1] - y[-2]) / np.log(d[-2] / d[-1])
    return float(d[-1] * L)


def displacement_thickness(y, u, Ue, tail: bool = True) -> float:
    """Displacement thickness δ* = ∫₀^∞ (1 − u/U_e) dy  [m].

    Book: §9.2, Eq. (9.16), ∫(U_e − u)dy = U_e δ*: the height by which the outer streamlines are pushed away.
    Parameters: y [m] increasing from 0; u [m/s] on y; Ue [m/s] edge speed; tail: add an exponential-tail estimate beyond y[-1]
    (the deficit is assumed to decay like exp(−y/l) with l from the last two points).  Caveat: the Blasius deficit decays like a Gaussian,
    exp(−η²/4), faster than an exponential, so the tail estimate slightly OVER-estimates the missing integral there (the effect is small when
    y[-1] is well past δ₉₉); pass ``tail=False`` to integrate only the sampled range.
    Method: composite Simpson (scipy.integrate.simpson) on the sampled deficit.
    Validation: V1 exact on u/U = 1 − e^{−y/a} (δ* = a) and the sine profile.  Label: analytic.
    """
    y, u = _F(y), _F(u)
    d = 1.0 - u / Ue  # integrand of Eq. (9.16)
    val = simpson(d, x=y)
    return float(val + (_deficit_tail(y, d) if tail else 0.0))


def momentum_thickness(y, u, Ue, tail: bool = True) -> float:
    """Momentum thickness θ = ∫₀^∞ (u/U_e)(1 − u/U_e) dy  [m].

    Book: §9.2, Eq. (9.17) (ρU_e²θ is the momentum-flux deficit; an angle everywhere else in the book).
    Parameters/method as :func:`displacement_thickness` (including its ``tail=True`` caveat: the exponential-tail estimate over-estimates the
    tail of a Gaussian-decaying Blasius deficit slightly).  Validation: V1 exact (a/2 on 1 − e^{−y/a}).  Label: analytic.
    """
    y, u = _F(y), _F(u)
    r = u / Ue
    g = r * (1.0 - r)  # Eq. (9.17)
    val = simpson(g, x=y)
    return float(val + (_deficit_tail(y, g) if tail else 0.0))


def delta_level(y, u, Ue, level: float = 0.99) -> float:
    """Height where u = level·U_e (δ₉₉ for level = 0.99)  [m].

    Book: §9.2 (δ₉₉ definition; other percentages are equally arbitrary).  Method: PCHIP of u/U_e(y) + ``brentq``.
    Raises ValueError if the profile never reaches ``level``.  Label: analytic (root of an interpolant).
    """
    y, u = _F(y), _F(u)
    r = u / Ue
    if r.max() < level:
        raise ValueError("profile never reaches the requested level")
    p = PchipInterpolator(y, r - level)
    i = int(np.argmax(r >= level))
    return float(brentq(lambda s: float(p(s)), y[max(i - 1, 0)], y[i], xtol=1e-14))


def thicknesses(y, u, Ue, level: float = 0.99) -> dict:
    """δ₉₉ (``level``), δ*, θ [m] and shape factor H = δ*/θ [–] of a sampled profile (Eqs. (9.16), (9.17))."""
    ds, th = displacement_thickness(y, u, Ue), momentum_thickness(y, u, Ue)
    return dict(delta99=delta_level(y, u, Ue, level), delta_star=ds, theta=th, H=ds / th)


def profile_shape(name: str, y, delta: float, p: float | None = None) -> dict:
    """Closed-form model profiles with their exact δ*, θ, H and δ₉₉ (for the thicknesses explainer and primer checks).

    Book: §9.2 Eqs. (9.16)–(9.17) applied to standard shapes.  ``name``: "linear" u/U = y/δ (δ* = δ/2, θ = δ/6, H = 3, δ₉₉ = 0.99δ);
    "sine" sin(πy/2δ) (δ(1−2/π), δ(2/π−½), δ₉₉ = 0.910δ); "cubic" 3η/2 − η³/2 (3δ/8, 39δ/280, δ₉₉ = 0.918δ); "exponential"
    1 − e^{−y/a} with δ ≡ a (δ* = a, θ = a/2, H = 2, δ₉₉ = a ln 100 = 4.605a); "power" u/U = (y/δ)^{1/p}, p ≥ 1 (``p`` default 7, the turbulent
    one-seventh law; δ* = δ/(p+1), θ = pδ/((p+1)(p+2)), H = (p+2)/p, δ₉₉ = 0.99^p δ); "blasius" f′(y/δ) with δ = √(νx/U) (δ* = 1.7208δ,
    θ = 0.6641δ, δ₉₉ = 4.910δ).  Profiles equal 1 beyond y = δ (except "exponential", "blasius").
    Parameters: y [m] (array or float); delta [m] (the scale of the profile); p [–] exponent of "power" only.
    Returns dict(u_over_Ue [–] at y, delta_star [m], theta [m], H [–], delta99 [m]) — the four numbers are exact closed forms (not quadrature).
    Validation: V1 (integrals of the profile vs the closed forms).  Label: analytic.
    """
    y = _F(y)
    eta = y / delta
    e = np.minimum(eta, 1.0)
    if name == "linear":
        r, ds, th, d99 = e, delta / 2, delta / 6, 0.99 * delta
    elif name == "sine":
        r, ds, th, d99 = np.where(eta < 1, np.sin(np.pi * e / 2), 1.0), delta * (1 - 2 / np.pi), delta * (2 / np.pi - 0.5), delta * 2 / np.pi * np.arcsin(0.99)
    elif name == "cubic":
        r, ds, th = 1.5 * e - 0.5 * e ** 3, 3 * delta / 8, 39 * delta / 280
        d99 = delta * brentq(lambda s: 1.5 * s - 0.5 * s ** 3 - 0.99, 0.0, 1.0, xtol=1e-15)
    elif name == "exponential":
        r, ds, th, d99 = 1.0 - np.exp(-eta), delta, delta / 2, delta * np.log(100.0)
    elif name == "power":
        pp = 7.0 if p is None else float(p)
        if pp < 1.0:
            raise ValueError("power profile needs p >= 1")
        r = e ** (1.0 / pp)
        ds, th, d99 = delta / (pp + 1.0), pp * delta / ((pp + 1.0) * (pp + 2.0)), 0.99 ** pp * delta
    elif name == "blasius":
        c = blasius_constants()
        r, ds, th, d99 = _blasius().fp(eta), c["delta_star"] * delta, c["theta"] * delta, c["eta99"] * delta
    else:
        raise ValueError('name must be "linear", "sine", "cubic", "exponential", "power" or "blasius"')
    return dict(u_over_Ue=r, delta_star=float(ds), theta=float(th), H=float(ds / th), delta99=float(d99))


# ======================================================================================================================
# §9.3 Blasius — Töpfer scaling
# ======================================================================================================================


@dataclass
class _BlasiusSolution:
    a: float
    F: Callable
    f: Callable
    fp: Callable
    fpp: Callable
    fpp0: float
    delta_star: float
    eta_end: float


@functools.lru_cache(maxsize=1)
def _blasius() -> _BlasiusSolution:
    """Blasius (9.27)–(9.29) by Töpfer scaling: solve F‴ + ½FF″ = 0, F(0) = F′(0) = 0, F″(0) = 1 as an initial-value problem
    and rescale f(η) = aF(aη), a = F′(∞)^{−1/2}: the equation is invariant under f → aF(aη), so f″(0) = a³ and f′(∞) = 1."""
    xi_max = 16.0
    sol = solve_ivp(lambda s, Y: [Y[1], Y[2], -0.5 * Y[0] * Y[2]], (0.0, xi_max), [0.0, 0.0, 1.0], method="DOP853",
                    rtol=1e-13, atol=1e-15, dense_output=True)
    c = sol.y[1, -1]  # F'(inf)
    a = c ** -0.5
    eta_end = xi_max / a
    ds = (xi_max - a ** 2 * sol.y[0, -1]) / a  # lim (η − f)

    def _ev(eta, k):
        eta = _F(eta)
        xi = np.clip(a * eta, 0.0, xi_max)
        Y = sol.sol(xi.ravel()).reshape((3,) + eta.shape)  # dense output takes 1-D input only
        far = eta > eta_end
        if k == 0:
            return np.where(far, eta - ds, a * Y[0])
        if k == 1:
            return np.where(far, 1.0, a ** 2 * Y[1])
        return np.where(far, 0.0, a ** 3 * Y[2])

    return _BlasiusSolution(a, sol.sol, lambda e: _S(_ev(e, 0)), lambda e: _S(_ev(e, 1)), lambda e: _S(_ev(e, 2)), a ** 3, ds, eta_end)


@functools.lru_cache(maxsize=1)
def _blasius_const_cache() -> dict:
    b = _blasius()
    eta99 = brentq(lambda e: b.fp(e) - 0.99, 2.0, 8.0, xtol=1e-14)  # root of f′ = 0.99  (NOT the figure reading 4.93)
    theta = quad(lambda e: b.fp(e) * (1.0 - b.fp(e)), 0.0, b.eta_end, epsabs=1e-13, epsrel=1e-13, limit=200)[0]  # Eq. (9.17)
    return dict(fpp0=b.fpp0, eta99=eta99, delta_star=b.delta_star, theta=theta, H=b.delta_star / theta, v_inf=0.5 * b.delta_star,
                tau_coeff=b.fpp0, cf_coeff=2.0 * b.fpp0, cd_coeff=4.0 * b.fpp0)


def blasius_constants() -> dict:
    """The numbers of the Blasius solution, COMPUTED (never typed): f″(0), η₉₉, δ*, θ, H, v∞, and the drag coefficients.

    Book: §9.3, Eqs. (9.27)–(9.33): f″(0) = 0.332 (τ₀ = f″(0)ρU²/√Re_x, Eq. (9.31)); η₉₉ (Eq. (9.30) prints 4.93, the root of
    f′ = 0.99 is 4.910 — slip R2); δ*/√(νx/U) = lim(η − f) (1.72); θ/√(νx/U) = ∫f′(1−f′)dη = 2f″(0) (0.664); H = δ*/θ; v∞√Re_x/U
    = ½ lim(ηf′ − f) = δ*/2 (0.86, Fig. 9.6).
    Returns dict(fpp0, eta99, delta_star, theta, H, v_inf, tau_coeff (= f″(0)), cf_coeff (= 2f″(0), Eq. (9.32)),
    cd_coeff (= 4f″(0), Eq. (9.33))), all non-dimensional, scaled by √(νx/U) or √Re_x as stated.
    Validation: V2 residual of (9.27); V3 Töpfer vs solve_bvp (1e-9); V5 f″(0) = 0.332057336215 (Howarth 1938 / Wikipedia
    "Blasius boundary layer"); V7 θ = 2f″(0).  Label: symbolic, converged (Töpfer vs solve_bvp asserted), benchmark.
    """
    return dict(_blasius_const_cache())


def blasius_profile(eta):
    """Blasius f(η), f′(η), f″(η) at arbitrary η (arrays ok): returns (f, fp, fpp).  Book: §9.3 (9.19), (9.23), (9.27).
    Beyond η ≈ 16/a the far field f′ = 1, f = η − δ* is used."""
    b = _blasius()
    return b.f(eta), b.fp(eta), b.fpp(eta)


def blasius_fields(x, y, U: float = 1.0, nu: float = 1e-6) -> dict:
    """Velocity field of the laminar flat-plate boundary layer (Blasius).

    Book: §9.3, Eqs. (9.19), (9.23), (9.24), (9.26): δ = √(νx/U) (9.26); η = y/δ; ψ = Uδ f(η) (9.19); u = ψ_y = U f′ (9.23);
    v = −ψ_x = Uδ′(ηf′ − f) = ½√(νU/x)(ηf′ − f) (9.24); τ₀ = μU f″(0)/δ (9.31).
    Parameters: x [m] > 0 and y [m] ≥ 0 (broadcast against each other, e.g. a meshgrid); U [m/s]; nu [m²/s].
    Returns dict(u, v, psi, eta, delta, tau0_over_rho) — u, v [m/s], ψ [m²/s], η [–], δ [m], τ₀/ρ [m²/s²] = ν U f″(0)/δ.
    Assumptions: laminar, zero pressure gradient, x ≫ ν/U (leading-edge region not resolved).
    Validation: V4/V1 continuity and (9.18) residual O(h²); ∫τ₀dx = ρU²θ; V3 similarity collapse.  Label: analytic (no convergence order asserted).
    """
    x, y = np.broadcast_arrays(_F(x), _F(y))
    if np.any(x <= 0):
        raise ValueError("x must be > 0 (leading edge x = 0 is the singular point)")
    delta = np.sqrt(nu * x / U)  # Eq. (9.26)
    eta = y / delta
    f, fp, fpp = blasius_profile(eta)
    u = U * _F(fp)  # Eq. (9.23)
    v = 0.5 * np.sqrt(nu * U / x) * (eta * _F(fp) - _F(f))  # Eq. (9.24)
    psi = U * delta * _F(f)  # Eq. (9.19)
    return dict(u=_S(u), v=_S(v), psi=_S(psi), eta=_S(eta), delta=_S(delta),
                tau0_over_rho=_S(nu * U * _blasius().fpp0 / delta))


def blasius_delta99(x, U: float, nu: float, printed: bool = False):
    """δ₉₉ = η₉₉ √(νx/U)  [m], with η₉₉ = 4.910 (the root of f′ = 0.99); ``printed=True`` uses the book's figure reading 4.93.

    Book: §9.3, Eq. (9.30) (printed 4.93; slip R2 — the planted wrong variant a test must fail).  Label: analytic.
    """
    eta99 = 4.93 if printed else _blasius_const_cache()["eta99"]
    return _S(eta99 * np.sqrt(nu * _F(x) / U))  # Eq. (9.30)


def blasius_delta_star(x, U: float, nu: float):
    """Displacement thickness δ* = 1.7208 √(νx/U)  [m] (Eq. (9.16) on the Blasius profile; coefficient lim(η − f), computed).

    Book: §9.3, Eqs. (9.16), (9.26), (9.30)–(9.33) (the book quotes the coefficient 1.72).  x [m] > 0; U [m/s]; nu [m²/s].  Label: analytic."""
    return _S(_blasius().delta_star * np.sqrt(nu * _F(x) / U))  # Eq. (9.16) with (9.26)


def blasius_theta(x, U: float, nu: float):
    """Momentum thickness θ = 0.6641 √(νx/U) = 2f″(0)√(νx/U)  [m] (Eq. (9.17) on the Blasius profile; coefficient ∫f′(1−f′)dη = 2f″(0), computed).

    Book: §9.3, Eqs. (9.17), (9.26), (9.31) (coefficient 0.664).  x [m] > 0; U [m/s]; nu [m²/s].  Label: analytic."""
    return _S(blasius_constants()["theta"] * np.sqrt(nu * _F(x) / U))  # Eq. (9.17) with (9.26)


def blasius_wall_shear(x, U: float, rho: float, nu: float):
    """Wall shear stress τ₀ = f″(0) ρU²/√Re_x = 0.332ρU²/√Re_x  [Pa].  Book: §9.3, Eq. (9.31).  Label: analytic/benchmark."""
    Rex = U * _F(x) / nu
    return _S(_blasius().fpp0 * rho * U ** 2 / np.sqrt(Rex))  # Eq. (9.31)


def blasius_skin_friction(Re_x):
    """Local skin-friction coefficient C_f = τ₀/(½ρU²) = 2f″(0)/√Re_x = 0.664/√Re_x  [–].  Book: §9.3, Eq. (9.32)."""
    return _S(2.0 * _blasius().fpp0 / np.sqrt(_F(Re_x)))  # Eq. (9.32)


def blasius_drag(L, U: float, rho: float, nu: float, sides: int = 1):
    """Drag per unit width F_D = ∫₀ᴸ τ₀ dx = 2f″(0) ρU²L/√Re_L = 0.664ρU²L/√Re_L  [N/m] (∝ U^{3/2}), by ``quad`` of Eq. (9.31).

    Book: §9.3 (unnumbered, before Eq. (9.33)).  ``sides``: number of wetted faces — the book says ONE side of the plate (trap/note R5, not a printed slip); a plate
    wetted on both faces has twice the drag.  Label: analytic (quad vs closed form)."""
    val, err = quad(lambda s: float(blasius_wall_shear(s, U, rho, nu)), 0.0, float(L), epsabs=0, epsrel=1e-12, limit=200)
    return sides * val  # ∫ τ0 dx, integrable x^{-1/2} singularity


def blasius_drag_coefficient(Re_L, sides: int = 1):
    """C_D = F_D/(½ρU²L) = 4f″(0)/√Re_L = 1.328/√Re_L per wetted face  [–] (one side, as the book says; trap/note R5).  Book: §9.3, Eq. (9.33)."""
    return _S(sides * 4.0 * _blasius().fpp0 / np.sqrt(_F(Re_L)))  # Eq. (9.33)


@functools.lru_cache(maxsize=1)
def _far_prefactor() -> float:
    b = _blasius()
    e0 = 5.0
    xi = e0 - b.delta_star
    return float(b.fpp(e0) * np.exp(xi ** 2 / 4.0))


def blasius_far_field(eta):
    """Far field f′ − 1 ≈ −A√π erfc(ξ/2), ξ = η − δ*  (Gaussian approach to the free stream; large-ξ form −(2A/ξ)e^{−ξ²/4}).

    Book: §9.3 (unnumbered): f′ − 1 ~ (1/η)exp(−η²/4).  The book omits the prefactor and the shift; the linearisation about f ≈ η − δ* gives
    f‴ + ½ξ f″ = 0 ⇒ f″ = A e^{−ξ²/4} with A = f″(5)e^{(5−δ*)²/4} (≈ 0.234, taken once from the solved profile, NOT typed), and
    integrating from η to ∞ gives −A√π erfc(ξ/2).  Returns f′ − 1 [–] (negative).
    Validation: agrees with the numerical f′ − 1 for 4 ≤ η ≤ 8 (tested); large-ξ form checked.  Label: analytic (asymptotic).
    """
    xi = _F(eta) - _blasius().delta_star
    A = _far_prefactor()
    return _S(-A * np.sqrt(np.pi) * erfc(xi / 2.0))  # f′ − 1 = −∫_η^∞ f″ dη′


# ======================================================================================================================
# §9.4 Falkner–Skan
# ======================================================================================================================
_FS_CACHE: dict = {}
_FOLD = -0.09043  # approximate location of the saddle-node; used only to guard the search


def _fs_rhs_factory(m: float):
    return lambda e, Y: np.vstack([Y[1], Y[2], -(m + 1.0) / 2.0 * Y[0] * Y[2] + m * Y[1] ** 2 - m])


def _bvp_fs(m: float, eta_max: float, tol: float, guess):
    """One solve_bvp of Eq. (9.36) at fixed m with f(0) = f′(0) = 0 (9.28), f′(η_max) = 1 (9.29)."""
    e = np.linspace(0.0, eta_max, 200) if guess is None else guess[0]
    Y0 = None if guess is None else guess[1]
    if Y0 is None:
        b = _blasius()
        Y0 = np.vstack([b.f(e), b.fp(e), b.fpp(e)])
    sol = solve_bvp(_fs_rhs_factory(m), lambda a, b: np.array([a[0], a[1], b[1] - 1.0]), e, Y0, tol=tol, max_nodes=200000)
    return sol


def _fs_solution(m: float, eta_max: float, tol: float = 1e-10):
    """Cached attached-branch BVP solution with continuation from the nearest solved m (bisecting the step on failure)."""
    key = (round(float(m), 12), float(eta_max), float(tol))
    if key in _FS_CACHE:
        return _FS_CACHE[key]
    if abs(m) < 1e-15:
        sol = _bvp_fs(0.0, eta_max, tol, None)
        _FS_CACHE[key] = sol
        return sol
    # nearest cached m with same eta_max/tol (always try Blasius first)
    cands = [(abs(k[0] - m), k[0]) for k in _FS_CACHE if k[1] == float(eta_max) and k[2] == float(tol)]
    if not cands:
        _fs_solution(0.0, eta_max, tol)
        cands = [(abs(m), 0.0)]
    m_c = min(cands)[1]
    # waypoints keep every continuation hop small (a wide hop makes Newton fail and forces deep bisection).  Hops are wide (0.5·max(1, |m|)) beyond
    # m = 0.5 (whatever the sign of the target) and 0.5 / 0.03 for a positive / non-positive target near Blasius.  The step rule depends only on
    # the position of the hop, so the sub-solves below never generate waypoints of their own (no recursion, whatever the cache holds).
    def _step(cur):
        if cur > 0.5:
            return max(0.5, 0.5 * cur)
        return 0.5 if m > 0 else 0.03

    way = []
    cur = m_c
    while abs(m - cur) > _step(cur):
        cur = cur + np.sign(m - cur) * _step(cur)
        way.append(round(float(cur), 12))
    def go(m_from, m_to, depth=0):
        s0 = _FS_CACHE[(round(m_from, 12), float(eta_max), float(tol))]
        guess = (s0.x, s0.y)
        sol = _bvp_fs(m_to, eta_max, tol, guess)
        if sol.success:
            _FS_CACHE[(round(m_to, 12), float(eta_max), float(tol))] = sol
            return sol
        if depth > 18:
            _FS_CACHE[key] = sol
            return sol
        mid = 0.5 * (m_from + m_to)
        go(m_from, mid, depth + 1)
        return go(mid, m_to, depth + 1)

    prev = m_c  # chain the hops iteratively: each starts from the previous solution (a failed hop stores the failure under ``key`` and stops)
    for w_ in way:
        if (round(w_, 12), float(eta_max), float(tol)) not in _FS_CACHE:
            go(prev, w_)
            if key in _FS_CACHE:
                return _FS_CACHE[key]
        prev = w_
    return go(prev, float(m))


def _fs_evaluator(m: float, eta_max: float, tol: float = 1e-10):
    """Callable eta → (f, fp, fpp) for the attached-branch solution (Blasius by Töpfer scaling, far field beyond η_max)."""
    if abs(m) < 1e-15:
        b = _blasius()
        return lambda e: (_F(b.f(e)), _F(b.fp(e)), _F(b.fpp(e)))
    sol = _fs_solution(m, eta_max, tol)
    if not sol.success:
        raise RuntimeError(f"solve_bvp failed for m = {m} (attached branch exists only for m > −0.0904): {sol.message}")

    def ev(e):
        e = _F(e)
        ec = np.clip(e, 0.0, eta_max)
        Y = sol.sol(ec.ravel()).reshape((3,) + e.shape)
        far = e > eta_max
        f = np.where(far, Y[0] + (e - eta_max), Y[0])
        return f, np.where(far, 1.0, Y[1]), np.where(far, 0.0, Y[2])

    return ev


def _fs_by_s(s: float, eta_max: float, guess, tol: float = 1e-10):
    """Falkner–Skan with f″(0) = s prescribed and m an unknown PARAMETER (crosses the fold at m = −0.0904)."""
    def fun(e, Y, p):
        m = p[0]
        return np.vstack([Y[1], Y[2], -(m + 1.0) / 2.0 * Y[0] * Y[2] + m * Y[1] ** 2 - m])

    def bc(a, b, p):
        return np.array([a[0], a[1], a[2] - s, b[1] - 1.0])

    e, Y0, p0 = guess
    sol = solve_bvp(fun, bc, e, Y0, p=[p0], tol=tol, max_nodes=200000)
    return sol


def _fs_track_s(s_values, eta_max: float = 10.0, tol: float = 1e-9):
    """Follow the Falkner–Skan family by decreasing f″(0) from Blasius; returns list of (s, m, solution)."""
    b = _blasius()
    e = np.linspace(0.0, eta_max, 300)
    guess = (e, np.vstack([b.f(e), b.fp(e), b.fpp(e)]), 0.0)
    out = []
    for s in s_values:
        sol = _fs_by_s(float(s), eta_max, guess, tol)
        if not sol.success:
            break
        guess = (sol.x, sol.y, float(sol.p[0]))
        out.append((float(s), float(sol.p[0]), sol))
    return out


def falkner_skan(m: float, eta_max: float | None = None, method: str = "bvp", tol: float = 1e-10, branch: str = "attached",
                 n: int = 400) -> dict:
    """Falkner–Skan similarity solution f‴ + ((m+1)/2) f f″ − m f′² + m = 0 for U_e = a x^m (m is the book's n).

    Book: §9.4, Eq. (9.36) with f(0) = f′(0) = 0 (9.28), f′(∞) = 1 (9.29); m = 0 is the Blasius equation (9.27); m = 1 is the
    stagnation (Hiemenz) flow; f′(η) = u/U_e, η = y/δ(x), δ = √(νx/U_e) (9.34).
    Parameters: m [–] exponent (pressure gradient −dp/dx = m a² x^{2m−1}, Eq. (9.35)); eta_max [–] truncation of ∞ (default None → 12 for m = 0, 16 for m ≤ −0.05 where the layer is thick (m = −0.09: f″(0) = 0.018909 at 10 vs 0.018872 at 16/30),
    10 otherwise; the answer is insensitive to it, tested for 8, 12, 16); method "bvp" (solve_bvp with continuation in m; default), "toepfer" (m = 0 only,
    scaling of an initial-value problem) or "shoot" (brentq on f″(0), m ≥ −0.05 only); tol [–] solver tolerance; branch
    "attached" (default; exists for m > −0.09043) or "reversed" (Stewartson's second branch, −0.0904 < m < 0, f″(0) < 0; the
    book's remark that reverse-flow solutions exist below −0.0904 is not supported: no attached bounded (0 ≤ f′ ≤ 1) solution exists below the fold, R11); n number of output points.
    Returns dict(eta, f, fp, fpp, fppp, fpp0, m, success, method, eta_max) — f‴ from the ODE itself; all non-dimensional.
    DEVIATION: the book gives no method; the BVP/Töpfer/shooting choices and continuation are ours (math-to-python §2).
    Validation: V2 sympy residual of (9.36); V3 bvp vs Töpfer 1e-9 and truncation; V5 f″(0) = 0.332057336 at m = 0, separation
    m_sep = −0.09043; V7 f‴(0) = −m.  Label: symbolic, converged (truncation insensitivity tested), benchmark.
    """
    m = float(m)
    if eta_max is None:  # default truncation: 12 for Blasius (Gaussian tail: f″(0) to 7e-12), 16 near the fold (m ≤ −0.05), 10 otherwise
        eta_max = 12.0 if m == 0.0 else (16.0 if m <= -0.05 else 10.0)
    eta = np.linspace(0.0, float(eta_max), int(n))
    success = True
    if branch == "reversed":
        if not (_FOLD < m < 0.0):
            raise ValueError("the reversed branch exists only for −0.0904 < m < 0")
        sols = _fs_track_reversed(m, eta_max=max(eta_max, 30.0), tol=1e-9)
        f, fp, fpp = sols["ev"](eta)
        success = bool(sols["success"])
        method = "bvp(s-continuation)"
    elif method == "toepfer":
        if m != 0.0:
            raise ValueError("Töpfer scaling applies to the Blasius equation only (m = 0)")
        b = _blasius()
        f, fp, fpp = _F(b.f(eta)), _F(b.fp(eta)), _F(b.fpp(eta))
    elif method == "shoot":
        if m < -0.05:
            raise ValueError("shooting is exponentially sensitive near separation; use method='bvp' for m < −0.05")
        s = _fs_shoot(m, eta_max)
        sol = solve_ivp(lambda e, Y: [Y[1], Y[2], -(m + 1.0) / 2.0 * Y[0] * Y[2] + m * Y[1] ** 2 - m], (0, eta_max), [0, 0, s],
                        method="DOP853", rtol=1e-13, atol=1e-15, t_eval=eta)
        f, fp, fpp = sol.y
        method = "shoot"
    elif method == "bvp":
        if m <= _FOLD:
            return dict(eta=eta, f=np.full_like(eta, np.nan), fp=np.full_like(eta, np.nan), fpp=np.full_like(eta, np.nan),
                        fppp=np.full_like(eta, np.nan), fpp0=float("nan"), m=m, success=False, method="bvp", eta_max=eta_max)
        sol = _fs_solution(m, float(eta_max), tol)
        success = bool(sol.success)
        Y = sol.sol(eta)
        f, fp, fpp = Y
    else:
        raise ValueError('method must be "bvp", "toepfer" or "shoot"')
    fppp = -(m + 1.0) / 2.0 * f * fpp + m * fp ** 2 - m  # Eq. (9.36) solved for f‴
    return dict(eta=eta, f=f, fp=fp, fpp=fpp, fppp=fppp, fpp0=float(fpp[0]), m=m, success=success, method=method, eta_max=float(eta_max))


def _fs_shoot(m: float, eta_max: float) -> float:
    """f″(0) by brentq on f′(η_max) − 1.

    Wrong guesses either overshoot (f′ climbs through 1.6: event, returns +0.6) or, for m > 0, undershoot (f″ turns negative
    before f′ reaches 1: event, returns f′ − 1 < 0 there, i.e. the peak height minus one).  The two events keep every integration
    finite (the overshoot has a finite-η blow-up for m > 0).
    """
    def g(s):
        def ev(e, Y):
            return Y[1] - 1.6
        ev.terminal = True

        def ev2(e, Y):
            return Y[2] if m > 0 else 1.0
        ev2.terminal = True
        ev2.direction = -1
        sol = solve_ivp(lambda e, Y: [Y[1], Y[2], -(m + 1.0) / 2.0 * Y[0] * Y[2] + m * Y[1] ** 2 - m], (0, eta_max), [0, 0, s],
                        method="DOP853", rtol=1e-12, atol=1e-14, events=[ev, ev2])
        if sol.t_events[0].size:
            return 0.6
        return sol.y[1, -1] - 1.0
    lo, hi = 1e-4, 3.0 + 2.0 * abs(m)
    while g(hi) < 0 and hi < 50:
        hi *= 2
    return brentq(g, lo, hi, xtol=1e-15, rtol=1e-14)


@functools.lru_cache(maxsize=8)
def _fs_fold_path(eta_max: float, tol: float):
    """Attached branch followed by decreasing f″(0) from Blasius to the fold f″(0) = 0 (list of (s, m, sol))."""
    b = _blasius()
    s_fwd = np.concatenate([np.linspace(b.fpp0, 0.05, 12), np.linspace(0.05, 0.0, 6)[1:]])
    return _fs_track_s(s_fwd, eta_max=eta_max, tol=tol)


@functools.lru_cache(maxsize=8)
def _fs_reversed_path(eta_max: float, tol: float):
    """Continue past the fold onto the second (reversed-flow) branch: f″(0) < 0, m increasing back towards 0."""
    path = list(_fs_fold_path(eta_max, tol))
    s0, m0, sol0 = path[-1]
    guess = (sol0.x, sol0.y, m0)
    out = []
    for s in -np.arange(1, 200) * 0.002:
        sol = _fs_by_s(float(s), eta_max, guess, tol)
        if not sol.success or sol.y[1].max() > 1.05 or sol.p[0] > 0.0:  # a physical member has f' <= 1 and m < 0 here
            break
        guess = (sol.x, sol.y, float(sol.p[0]))
        out.append((float(s), float(sol.p[0]), sol))
    # the (f″(0), m) curve turns again near m ≈ −0.058: continue with m as the parameter (single-valued from there to m → 0⁻)
    if out:
        sN, mN, solN = out[-1]
        gx, gy = solN.x, solN.y
        for mt in np.arange(mN + 0.002, -0.0057, 0.002):
            sl = _bvp_fs(float(mt), eta_max, tol, (gx, gy))
            if not sl.success or sl.y[2, 0] > 0.0:
                break
            gx, gy = sl.x, sl.y
            out.append((float(sl.y[2, 0]), float(mt), sl))
    return path + out


def _fs_track_reversed(m_target: float, eta_max: float, tol: float) -> dict:
    path = _fs_reversed_path(float(eta_max), float(tol))
    rev = [(s, mm, sol) for (s, mm, sol) in path if s <= 0.0]
    ms = np.array([r[1] for r in rev])
    ok = len(rev) > 2 and ms.max() >= m_target
    if not ok:
        return dict(success=False, ev=lambda e: (np.full_like(_F(e), np.nan),) * 3)
    i = int(np.argmax(ms >= m_target))
    sol = _bvp_fs(m_target, eta_max, tol, (rev[i][2].x, rev[i][2].y))

    def ev(e):
        e = _F(e)
        Y = sol.sol(np.clip(e, 0, eta_max).ravel()).reshape((3,) + e.shape)
        return Y[0], Y[1], Y[2]

    return dict(success=bool(sol.success), ev=ev)


def falkner_skan_separation() -> dict:
    """Falkner–Skan member with zero wall shear: m_sep where f″(0; m) = 0 (the fold of the attached branch).

    Book: §9.4 (n = −0.0904, β = 2n/(n+1) = −0.19884; Fig. 9.7's last curve).  Method: f″(0) is prescribed and m is the
    unknown parameter of ``solve_bvp`` (the fold makes f″(0) a good coordinate, m a bad one), followed from Blasius down to
    f″(0) = 0 (continuation in f″(0)).
    Returns dict(m_sep, beta_sep, fpp0_at_sep).  Validation: V5 m_sep = −0.09043, β = −0.19884 (J. Eng. Math. 2015);
    V7 f″(0) = 0.  Label: benchmark.
    """
    path = _fs_fold_path(16.0, 1e-9)
    hit = [r for r in path if r[0] == 0.0]
    if not hit:
        raise RuntimeError("continuation to f''(0) = 0 failed")
    ms = hit[0][1]
    return dict(m_sep=ms, beta_sep=2 * ms / (ms + 1.0), fpp0_at_sep=0.0)


def falkner_skan_state(m: float, eta_max: float | None = None, tol: float = 1e-10) -> dict:
    """Wall shear, thickness integrals, Thwaites parameters and inflection of one Falkner–Skan member.

    Book: §9.4 (9.34)–(9.36) with the thickness definitions (9.16)–(9.17) and Thwaites' (9.44)–(9.46): for U_e = a xⁿ, δ(x) = √(νx/U_e):
    δ*/δ = I_δ = ∫(1 − f′)dη, θ/δ = I_θ = ∫f′(1 − f′)dη, H = I_δ/I_θ, λ = (θ²/ν)dU_e/dx = m I_θ², l = f″(0)I_θ (τ₀ = μ(U_e/θ)l),
    C_f √Re_x = 2f″(0), and (∂²u/∂y²)_wall ∝ f‴(0) = −m (Eq. (9.36) at η = 0).
    ``eta_max`` [–] truncation of η = ∞: ``None`` (default) picks 10 for m ≥ −0.05 and 16 for m < −0.05 (near separation the layer is thick and
    η_max = 10 leaves a 4 % error in f″(0) at m = −0.0904).
    Returns dict(m, fpp0, fppp0, I_delta, I_theta, H, lam, l, cf_sqrtRex, inflection_eta (0.0 for m = 0, None for m > 0),
    separated).  All non-dimensional.  Validation: V1 momentum integral (9.43) closes for every m (1e-8); V7 f‴(0) = −m.
    Label: analytic, converged.
    """
    m = float(m)
    if eta_max is None:
        eta_max = 16.0 if m <= -0.05 else 10.0
    d = falkner_skan(m, eta_max, tol=tol, n=4001)
    if not d["success"]:
        raise RuntimeError(f"no attached solution for m = {m}")
    eta, f, fp, fpp, fppp = d["eta"], d["f"], d["fp"], d["fpp"], d["fppp"]
    I_d = float(eta[-1] - f[-1])  # ∫(1 − f′)dη = η_max − f(η_max) exactly (f(0) = 0)
    I_t = float(simpson(fp * (1.0 - fp), x=eta))
    infl = None
    if abs(m) < 1e-14:
        infl = 0.0
    elif m < 0:
        sc = np.where(np.diff(np.sign(fppp)) != 0)[0]
        if sc.size:
            i = sc[0]
            infl = float(eta[i] - fppp[i] * (eta[i + 1] - eta[i]) / (fppp[i + 1] - fppp[i]))
    return dict(m=m, fpp0=d["fpp0"], fppp0=float(fppp[0]), I_delta=I_d, I_theta=I_t, H=I_d / I_t, lam=m * I_t ** 2,
                l=d["fpp0"] * I_t, cf_sqrtRex=2.0 * d["fpp0"], inflection_eta=infl, separated=bool(d["fpp0"] <= 0))


def falkner_skan_thickness(x, m: float, a: float, nu: float):
    """Generic thickness δ(x) = √(νx/U_e) = √(νx^{1−m}/a)  [m]; grows for m < 1, constant for m = 1.  Book: §9.4 (unnumbered)."""
    return _S(np.sqrt(nu * _F(x) ** (1.0 - m) / a))


def falkner_skan_fields(x, y, m: float, a: float, nu: float, eta_max: float | None = None) -> dict:
    """Falkner–Skan velocity field for U_e = a xᵐ.

    Book: §9.4, Eq. (9.34): ψ = √(νxU_e) f(η), η = y√(a/ν) x^{(m−1)/2}; u = ψ_y = U_e f′; v = −ψ_x = −√(νa) x^{(m−1)/2}[((m+1)/2) f + ((m−1)/2) ηf′].
    x [m] > 0, y [m] (broadcast); a [m^{1−m}/s]; nu [m²/s]; eta_max [–] truncation (None → 12 for m = 0, 16 for m ≤ −0.05, 10 otherwise, as :func:`falkner_skan`).
    Returns dict(u, v, psi, eta, delta, Ue).  Label: analytic (formulas), converged (truncation insensitivity tested)."""
    x, y = np.broadcast_arrays(_F(x), _F(y))
    if eta_max is None:
        eta_max = 16.0 if float(m) <= -0.05 else 10.0
    ev = _fs_evaluator(float(m), float(eta_max))
    eta = y * np.sqrt(a / nu) * x ** ((m - 1.0) / 2.0)  # Eq. (9.34)
    f, fp, _ = ev(eta)
    Ue = a * x ** m
    u = Ue * fp
    v = -np.sqrt(nu * a) * x ** ((m - 1.0) / 2.0) * (0.5 * (m + 1.0) * f + 0.5 * (m - 1.0) * eta * fp)
    psi = np.sqrt(nu * a) * x ** ((m + 1.0) / 2.0) * f
    return dict(u=_S(u), v=_S(v), psi=_S(psi), eta=_S(eta), delta=_S(np.sqrt(nu * x / Ue)), Ue=_S(Ue))


def _fs_default_grid() -> np.ndarray:
    """41 values of m from just above the fold to 4, dense near the fold (where everything changes fast)."""
    m_fold = -0.09043
    near = m_fold + np.geomspace(1e-4, 0.09043, 14)[:-1]  # −0.09033 … just below 0, geometric spacing from the fold
    mid = np.linspace(-0.01, 1.0, 20)
    hi = np.linspace(1.25, 4.0, 12)
    return np.unique(np.round(np.concatenate([near, mid, hi]), 6))


def falkner_skan_table(m_grid=None, n_eta: int = 161, eta_max: float = 8.0, fast: bool = False) -> dict:
    """Tabulated Falkner–Skan solutions (9.36) for the explainers (6 s.f. are enough there).

    Book: §9.4, Eq. (9.36) with (9.28), (9.29); thickness integrals (9.16)–(9.17); Thwaites' λ, l of (9.44)–(9.45); wall curvature (9.51)–(9.52).
    Parameters: m_grid [–] exponents m of U_e = a xᵐ (default: 41 nodes from the fold −0.09043 + 1e-4 to 4, dense near the fold); n_eta points
    of the η grid on [0, eta_max] (eta_max = 8 is the TABLE range; the solutions themselves are computed with the larger truncation of
    :func:`falkner_skan_state`, so the integrals are not truncated at 8); fast: n_eta ≤ 81, every second default node, and the reversed branch
    skipped (``fpp0_reversed`` = nan) so the call stays within a notebook cell budget.
    Returns dict(m, eta, fp[m, η], fpp[m, η], f[m, η], fpp0, I_delta, I_theta, H, lam, l, inflection_eta (nan where none),
    fpp0_reversed (second branch f″(0) < 0, nan outside −0.0904 < m ≲ −0.008: the BVP continuation of the second branch stops there, where its reversed region becomes very thick)) as float arrays.
    Attached branch only for the profiles.  Precomputed tables are OURS, not the book's.  Validation: rows agree with :func:`falkner_skan_state`
    (identical numbers); f″(0) = 0.33206 at m = 0; V1 momentum integral closes.  Label: analytic (identities); own numbers, no convergence order asserted.
    """
    ms = _fs_default_grid() if m_grid is None else np.asarray(m_grid, float)
    if fast:
        n_eta = min(int(n_eta), 81)
        if m_grid is None:
            ms = ms[::2]
    eta = np.linspace(0.0, float(eta_max), int(n_eta))
    F, FP, FPP = [], [], []
    keys = ("fpp0", "I_delta", "I_theta", "H", "lam", "l")
    cols = {k: [] for k in keys}
    infl, rev = [], []
    for m in ms:
        st = falkner_skan_state(float(m), tol=1e-8)  # 8 digits are ample for a table and 10× cheaper than 1e-10
        ev = _fs_evaluator(float(m), 16.0 if m < -0.05 else 10.0, 1e-8)
        f, fp, fpp = ev(eta)
        F.append(_F(f)), FP.append(_F(fp)), FPP.append(_F(fpp))
        for k in keys:
            cols[k].append(st[k])
        infl.append(np.nan if st["inflection_eta"] is None else st["inflection_eta"])
        r = np.nan
        if (not fast) and _FOLD < m < 0.0:
            try:
                d = falkner_skan(float(m), branch="reversed", n=5)
                r = d["fpp0"] if d["success"] else np.nan
            except (ValueError, RuntimeError):
                r = np.nan
        rev.append(r)
    out = dict(m=ms, eta=eta, f=np.array(F), fp=np.array(FP), fpp=np.array(FPP), inflection_eta=np.array(infl), fpp0_reversed=np.array(rev, float))
    out.update({k: np.array(v, float) for k, v in cols.items()})
    return out


# ======================================================================================================================
# §9.5 momentum integral, Kármán–Pohlhausen
# ======================================================================================================================


def momentum_integral_residual(x, Ue, theta, delta_star, tau0, rho: float = 1.0):
    """Residual of the von Kármán momentum integral equation τ₀/ρ − d(U_e²θ)/dx − U_e δ* dU_e/dx  [m²/s²].

    Book: §9.5, Eq. (9.43): (1/ρ)τ₀ = d/dx[U_e²θ] + U_e δ* dU_e/dx (laminar or time-averaged turbulent).
    Parameters: x [m] increasing; Ue [m/s], theta [m], delta_star [m], tau0 [Pa] sampled on x; rho [kg/m³].
    The x-derivatives use cubic splines through the samples (error O(h³)); avoid the first points near a singular leading edge.
    Validation: ≈ 0 (1e-8 relative) for the exact Blasius / Falkner–Skan thicknesses.  Label: analytic (identity).
    """
    x = _F(x)
    Ue, theta, ds, tau0 = _F(Ue), _F(theta), _F(delta_star), _F(tau0)
    d1 = CubicSpline(x, Ue ** 2 * theta)(x, 1)
    dU = CubicSpline(x, Ue)(x, 1)
    return tau0 / rho - d1 - Ue * ds * dU  # Eq. (9.43)


def _kp_coeffs():
    """Pohlhausen profile shape integrals as functions of Λ (sympy, lambdified once): c_theta, c_delta, c_tau."""
    import sympy as sp
    e, L = sp.symbols("eta Lambda", real=True)
    prof = {
        "cubic": (sp.Rational(3, 2) * e - e ** 3 / 2, False),
        "sine": (sp.sin(sp.pi * e / 2), False),
        "quartic": (2 * e - 2 * e ** 3 + e ** 4 + L / 6 * e * (1 - e) ** 3, True),
    }
    out = {}
    for k, (F, dep) in prof.items():
        ct = sp.integrate(F * (1 - F), (e, 0, 1))
        cd = sp.integrate(1 - F, (e, 0, 1))
        cs = sp.diff(F, e).subs(e, 0)
        out[k] = tuple(sp.lambdify(L, q, "numpy") for q in (ct, cd, cs))
    return out


_KP = None


def karman_pohlhausen(Ue, x, nu, profile: str = "cubic", theta0: float = 0.0, Z0=None, rho: float = 1.0) -> dict:
    """Assumed-profile solution of the momentum integral equation (Pohlhausen's idea; closure of Eq. (9.43)).

    Book: §9.5 — "three unknowns (θ, δ*, τ₀), one equation": choose a profile shape u/U_e = F(y/δ) and solve (9.43) for δ(x).
    Profiles (ours, not Exercise 9.14's numbers): "cubic" 3η/2 − η³/2, "sine" sin(πη/2), "quartic" the Pohlhausen family
    2η − 2η³ + η⁴ + (Λ/6)η(1−η)³ with Λ = δ²U_e′/ν.  With Z = δ², θ = δc_θ(Λ), δ* = δc_δ(Λ), τ₀ = μU_e c_τ(Λ)/δ, Eq. (9.43) is
    solved for Z′ (partial derivatives by central differences) and integrated with ``solve_ivp`` (DOP853).
    Parameters: Ue outer flow (OuterFlow/callable); x [m] increasing with x[0] > 0; nu [m²/s]; theta0 [m] initial momentum thickness at x[0]
    (0 → start from the local flat-plate/similarity balance at x[0], i.e. a layer that began at the leading edge); Z0 initial δ² [m²]
    (overrides theta0; θ₀ maps to δ₀ = θ₀/c_θ(0)); rho [kg/m³].  Returns dict(x, delta, theta, delta_star, tau0, lam) in SI.
    DEVIATION: shape integrals are computed by sympy, not typed; the initial value is our flat-plate estimate.
    Validation: cubic profile on a flat plate: θ/√(νx/U) within 6 % of Blasius (the ODE integration is converged; the assumed profile is not).
    Label: approximate (6 % vs Blasius).
    """
    global _KP
    if _KP is None:
        _KP = _kp_coeffs()
    ct, cd, cs = _KP[profile]
    of = _as_outer(Ue, x)
    x = _F(x)

    def shape(Z, xx):
        Lam = Z * float(of.dUe(xx)) / nu
        c1, c2, c3 = float(ct(Lam)), float(cd(Lam)), float(cs(Lam)) if not callable(cs) else float(cs(Lam))
        return Lam, c1, c2, c3

    def theta_fn(Z, xx):
        return np.sqrt(Z) * float(ct(Z * float(of.dUe(xx)) / nu))

    def rhs(xx, Z):
        Z = Z[0]
        U, dU = float(of.Ue(xx)), float(of.dUe(xx))
        Lam, c1, c2, c3 = shape(Z, xx)
        d = np.sqrt(Z)
        tau0_rho = nu * U * c3 / d
        hz, hx = 1e-6 * Z, 1e-6 * (1 + abs(xx))
        th_Z = (theta_fn(Z + hz, xx) - theta_fn(Z - hz, xx)) / (2 * hz)
        th_x = (theta_fn(Z, xx + hx) - theta_fn(Z, xx - hx)) / (2 * hx)
        th = c1 * d
        rest = tau0_rho - 2 * U * dU * th - U ** 2 * th_x - U * (c2 * d) * dU
        return [rest / (U ** 2 * th_Z)]

    U0 = float(of.Ue(x[0]))
    c1, c2, c3 = float(ct(0.0)), float(cd(0.0)), float(cs(0.0))
    if Z0 is None:
        Z0 = 2.0 * nu * c3 * x[0] / (c1 * U0) if theta0 == 0.0 else (float(theta0) / c1) ** 2
    Z0 = float(Z0)
    sol = solve_ivp(rhs, (x[0], x[-1]), [Z0], t_eval=x, method="DOP853", rtol=1e-8, atol=1e-14)
    Z = sol.y[0]
    xs = sol.t
    Ues, dUs = _F(of.Ue(xs)), _F(of.dUe(xs))
    Lam = Z * dUs / nu
    d = np.sqrt(Z)
    th = np.array([np.sqrt(z) * float(ct(l_)) for z, l_ in zip(Z, Lam)])
    dstar = np.array([np.sqrt(z) * float(cd(l_)) for z, l_ in zip(Z, Lam)])
    tau = rho * nu * Ues * np.array([float(cs(l_)) for l_ in Lam]) / d
    return dict(x=xs, delta=d, theta=th, delta_star=dstar, tau0=tau, lam=Lam)


# ======================================================================================================================
# §9.6 Thwaites
# ======================================================================================================================


def holstein_bohlen(theta, Ue_x, nu):
    """Holstein–Bohlen parameter λ = (θ²/ν) dU_e/dx  [–].  Book: §9.6, Eq. (9.44) (λ < 0 adverse; Thwaites' m = −λ)."""
    return _S(_F(theta) ** 2 / nu * _F(Ue_x))  # Eq. (9.44)


def _lambda_sep_fs() -> float:
    """λ = m I_θ² at the fold of the attached Falkner–Skan branch (l = 0): −0.068148, COMPUTED from :func:`_fs_fold_state` (never typed)."""
    return float(_fs_fold_state()["lam"])


def __getattr__(name):  # lazy module constant: LAMBDA_SEP_FS is computed on first access (the fold continuation costs ~2 s)
    """``LAMBDA_SEP_FS`` — λ at which the exact-Falkner–Skan shear l(λ) vanishes (fold of the attached branch), computed = −0.068148.
    Default separation criterion of :func:`thwaites` for ``closure="falkner_skan"``; the book's rounded criterion is :data:`LAMBDA_SEP_BOOK`."""
    if name == "LAMBDA_SEP_FS":
        return _lambda_sep_fs()
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


@functools.lru_cache(maxsize=1)
def _fs_fold_state() -> dict:
    """(m, λ, l, H) of the member with zero wall shear (the fold), from the f″(0)-continuation path: λ = mI_θ², l = 0, H = I_δ/I_θ."""
    path = _fs_fold_path(16.0, 1e-9)
    hit = [r for r in path if r[0] == 0.0]
    if not hit:
        raise RuntimeError("continuation to f''(0) = 0 failed")
    _, m, sol = hit[0]
    eta = np.linspace(0.0, 16.0, 4001)
    f, fp = sol.sol(eta)[0], sol.sol(eta)[1]
    I_d = float(eta[-1] - f[-1])
    I_t = float(simpson(fp * (1.0 - fp), x=eta))
    return dict(m=float(m), lam=float(m * I_t ** 2), l=0.0, H=I_d / I_t)


@functools.lru_cache(maxsize=4)
def _closure_table(n_points: int):
    """(λ, l, H, m) from the exact Falkner–Skan family, sorted by λ (l = f″(0)I_θ, H = I_δ/I_θ, λ = mI_θ²), ending at the zero-shear fold (l = 0)."""
    n = int(n_points)
    n1 = max(n // 4, 6)
    ms = np.unique(np.concatenate([
        np.linspace(-0.0900, -0.02, n1), np.linspace(-0.02, 0.0, 5), np.linspace(0.0, 1.0, n // 2 + 1),
        np.geomspace(1.0, 50.0, n - n1 - n // 2 - 5 + 5)]))
    rows = []
    for m in ms:
        try:
            s = falkner_skan_state(float(m), tol=1e-8)
        except RuntimeError:
            continue
        rows.append((s["lam"], s["l"], s["H"], float(m)))
    fold = _fs_fold_state()
    rows.append((fold["lam"], 0.0, fold["H"], fold["m"]))
    rows = np.array(sorted(rows))
    keep = np.concatenate([[True], np.diff(rows[:, 0]) > 1e-12])
    return rows[keep]


def thwaites_closure_table(n_points: int = 60, fast: bool = False) -> dict:
    """The shear correlation l(λ), shape factor H(λ) and L(λ) computed from OUR exact Falkner–Skan family (public replacement for Table 9.1).

    Book: §9.6, Table 9.1 (private) and Eqs. (9.45), (9.46), (9.48); λ = m I_θ², l = f″(0)I_θ, H = I_δ/I_θ, L = 2l − 2(2 + H)λ for m from the fold
    (−0.0904, where l = 0) to 50.
    DEVIATION: Thwaites' l, H are an empirical cross-family fit valid for λ ∈ [−0.09, 0.25]; the exact Falkner–Skan values cover only
    λ ∈ [−0.06815, 0.1065] (computed) and give l = 0 at λ = −0.06815 instead of −0.09 — both separation criteria are reported by
    :func:`thwaites`.  ``n_points``: number of family members (60; 24 when ``fast``).
    Accuracy: the members are solved with tolerance 1e-8, so the node values agree with :func:`falkner_skan_state` to ≈ 2e-8 (not 1e-9).
    Returns dict(m, lam, l, H, L) arrays sorted by λ.  Label: converged (own numbers; not asserted equal to the book)."""
    r = _closure_table(24 if fast else int(n_points))  # fast: 24 members instead of 60
    lam, l, H = r[:, 0], r[:, 1], r[:, 2]
    return dict(m=r[:, 3], lam=lam, l=l, H=H, L=2.0 * l - 2.0 * (2.0 + H) * lam)


def _white_l(lam):
    return np.clip(_F(lam) + 0.09, 0.0, None) ** 0.62


def thwaites_l(lam, closure: str = "falkner_skan"):
    """Shear correlation l(λ) with τ₀ = μ(U_e/θ)l(λ).  Book: §9.6, Eq. (9.45), Table 9.1.

    closure "falkner_skan": PCHIP through :func:`thwaites_closure_table` (exact for those flows); below the fold (λ < −0.06815) l = 0 (no attached
    solution — τ₀ and the separation flag agree) and above the table the shape of the "white" law scaled to match.
    closure "white": l ≈ (λ + 0.09)^0.62 (White's fit as quoted in an ANSYS lesson handout — a secondary source, not a benchmark; l(0) = 0.2247 vs Table 9.1 0.220;
    zero at −0.09, the book's criterion).
    Scalar-callable; returns [–].  Label: secondary-sourced fit ("white"); own table, no convergence order asserted ("falkner_skan")."""
    lam = _F(lam)
    if closure == "white":
        return _S(_white_l(lam))
    if closure != "falkner_skan":
        raise ValueError('closure must be "falkner_skan" or "white"')
    t = thwaites_closure_table(60)
    p = PchipInterpolator(t["lam"], t["l"])
    lo, hi = t["lam"][0], t["lam"][-1]
    inside = p(np.clip(lam, lo, hi))
    above = t["l"][-1] * _white_l(lam) / _white_l(hi)
    return _S(np.where(lam < lo, 0.0, np.where(lam > hi, above, inside)))


def thwaites_H(lam, closure: str = "falkner_skan"):
    """Shape factor H(λ) = δ*/θ.  Book: §9.6, Eq. (9.46).  Falkner–Skan closure (PCHIP; linear extrapolation outside the table,
    flagged); the "white" closure has no published H and reuses the exact-FS table.  [–]  Label: converged."""
    lam = _F(lam)
    t = thwaites_closure_table(60)
    p = PchipInterpolator(t["lam"], t["H"])
    lo, hi = t["lam"][0], t["lam"][-1]
    sl, sh = float(p(lo, 1)), float(p(hi, 1))
    v = np.where(lam < lo, t["H"][0] + sl * (lam - lo), np.where(lam > hi, t["H"][-1] + sh * (lam - hi), p(np.clip(lam, lo, hi))))
    return _S(v)


def thwaites_L(lam, closure: str = "falkner_skan"):
    """L(λ) = 2l(λ) − 2(2 + H(λ))λ  [–], Eq. (9.48); nearly linear, L ≈ 0.45 − 6.0λ (Fig. 9.8, Eq. (9.49)).  Label: converged."""
    lam = _F(lam)
    return _S(2.0 * _F(thwaites_l(lam, closure)) - 2.0 * (2.0 + _F(thwaites_H(lam, closure))) * lam)  # Eq. (9.48)


def _gl_cumulative(fun, x, order: int = 10):
    """Cumulative ∫ₓ₀ˣ fun dx′ at the grid points by Gauss–Legendre on each interval (vectorised)."""
    xg, wg = np.polynomial.legendre.leggauss(order)
    a, b = x[:-1], x[1:]
    mid, half = 0.5 * (a + b), 0.5 * (b - a)
    pts = mid[:, None] + half[:, None] * xg[None, :]
    vals = fun(pts)
    seg = half * (vals * wg[None, :]).sum(axis=1)
    return np.concatenate([[0.0], np.cumsum(seg)])


def thwaites(x, Ue, nu: float, theta0: float = 0.0, closure: str = "falkner_skan", rho: float = 1.0,
             lam_sep: float | None = None, stop_at_separation: bool = True) -> dict:
    """Thwaites' method: θ from the closed form of Eq. (9.50), then λ, H, τ₀, δ*, C_f and the separation point.

    Book: §9.6, Eqs. (9.44)–(9.50): θ²U_e⁶/ν = 0.45∫₀ˣU_e⁵dx′ + θ₀²U₀⁶/ν (9.50); λ = (θ²/ν)dU_e/dx (9.44); τ₀ = μ(U_e/θ)l(λ) (9.45);
    δ*/θ = H(λ) (9.46).  Separation is predicted where λ falls to ``lam_sep``: ``None`` → −0.06815 for the exact Falkner–Skan closure (where
    its l = 0, so τ₀ and the flag agree) and −0.09 for ``closure="white"`` (the book's criterion, Table 9.1: l(−0.09) = 0).
    Parameters: x [m] increasing from the start point x[0] (θ = θ₀ there; x[0] = 0 with U_e(0) = 0 is a stagnation point with the
    finite limit θ² = 0.45ν/(6 U_e′(0)), λ(0) = 0.075); Ue OuterFlow / callable / array [m/s]; nu [m²/s]; theta0 [m]; closure "falkner_skan"
    | "white"; rho [kg/m³]; lam_sep [–]; stop_at_separation: cut the arrays just after the first station with λ ≤ lam_sep (the method is
    abandoned after separation; ``False`` keeps the full arrays for plotting).
    Returns dict(x, theta, delta_star, tau0 [Pa], cf, lam, H, l, Ue, separated, x_sep, x_sep_l0) — x_sep None if not reached
    (first crossing of ``lam_sep``, linearly interpolated); x_sep_l0 = first crossing of l = 0.
    Assumptions: attached laminar layer, U_e from ideal flow; abandon the result once separation is predicted.
    DEVIATION: ∫U_e⁵ by 10-point Gauss–Legendre per interval (not the trapezoid) so graded coarse grids stay accurate.
    Validation: V1 Blasius θ = √(0.45νx/U) exactly; diffuser λ(x) = −(0.45/4)[(1+x/L)⁴ − 1] (1e-12); V5 vs exact Falkner–Skan within
    3 % (favourable) / 10 % (adverse).  Label: analytic, benchmark, book-value.
    """
    x = _F(x)
    if lam_sep is None:
        lam_sep = LAMBDA_SEP_BOOK if closure == "white" else _lambda_sep_fs()
    of = _as_outer(Ue, x)
    U = _F(of.Ue(x))
    dU = _F(of.dUe(x))
    I = _gl_cumulative(lambda s: _F(of.Ue(s)) ** 5, x)
    U0 = U[0]
    with np.errstate(divide="ignore", invalid="ignore"):
        th2 = (0.45 * nu * I + theta0 ** 2 * U0 ** 6) / U ** 6  # Eq. (9.50)
    # the analytic stagnation limit is needed only where (9.50) is 0/0 (U_e(x) underflows or is exactly 0); a small but non-zero U_e
    # (e.g. x0 = 1e-6 x on a power law U_e ∝ xⁿ, n ≥ 1.5, with the correct power-law theta0) is evaluated directly
    small = (U < 1e-9 * max(U.max(), 1e-300)) & (~np.isfinite(th2) | (U <= 0.0))
    if small.any():
        if theta0 != 0.0:
            raise ValueError("a stagnation start needs theta0 = 0")
        idx = np.where(small)[0]
        # analytic limit: U_e ≈ U_e′(0)x  ⇒ θ² = 0.45ν/(6U_e′)
        dU0 = float(of.dUe(x[idx]).mean()) if np.all(np.isfinite(of.dUe(x[idx]))) else float(U[1] / x[1])
        if dU0 <= 0.0:
            raise ValueError("θ diverges at a stagnation start with U_e′(0) = 0 (power law n > 1): start at x[0] > 0 with theta0")
        th2[idx] = 0.45 * nu / (6.0 * dU0)
    theta = np.sqrt(th2)
    lam = theta ** 2 / nu * dU  # Eq. (9.44)
    l = _F(thwaites_l(lam, closure))
    H = _F(thwaites_H(lam, closure))
    with np.errstate(divide="ignore", invalid="ignore"):
        tau0 = rho * nu * U / theta * l  # Eq. (9.45), μ = ρν
        cf = 2.0 * tau0 / (rho * U ** 2)
    delta_star = H * theta  # Eq. (9.46)

    def first_crossing(q, thr):
        below = np.where(q <= thr)[0]
        if below.size == 0:
            return None
        i = below[0]
        if i == 0:
            return float(x[0])
        return float(x[i - 1] + (thr - q[i - 1]) * (x[i] - x[i - 1]) / (q[i] - q[i - 1]))

    x_sep = first_crossing(lam, lam_sep)
    x_sep_l0 = first_crossing(l, 0.0)
    out = dict(x=x, theta=theta, delta_star=delta_star, tau0=tau0, cf=cf, lam=lam, H=H, l=l, Ue=U)
    if stop_at_separation and x_sep is not None:
        k = int(np.where(lam <= lam_sep)[0][0]) + 1  # keep the first station past separation, drop the rest
        out = {key: val[:k] for key, val in out.items()}
    return dict(out, separated=x_sep is not None, x_sep=x_sep, x_sep_l0=x_sep_l0)


def thwaites_named(kind: str, x, nu: float, theta0: float = 0.0, closure: str = "falkner_skan", **params) -> dict:
    """One-station Thwaites numbers for a named outer flow (scalar-callable, for explainer parity rows).

    Book: §9.6, Eqs. (9.44)–(9.46), (9.50).  ``kind`` and ``params`` as :func:`outer_flow`: "flat" (U), "diffuser" (U1, L), "wedge" (n, a),
    "cylinder" (U, a; x = aφ from the forward stagnation point), "retarded" (U0, L, c); an optional ``rho`` [kg/m³] (default 1) sets the units of τ₀.
    θ at the station x [m] comes from the marched integral of U_e⁵ on 4000 graded points x·s², s ∈ [0, 1] (10-point Gauss–Legendre per interval).
    For "wedge" (U_e = axⁿ) the march starts at x·1e-6 with the exact power-law value θ² = 0.45νx₀^{1−n}/(a(5n+1)) (the θ² ∝ x^{1−n} solution of (9.49));
    ``theta0`` [m] is the value of θ at x = 0 for the other kinds.
    Returns dict(theta [m], delta_star [m], tau0 [Pa when rho is in kg/m³], cf [–], lam [–], H [–], l [–], x_sep [m or None]) at x; x_sep = first crossing of
    the closure's separation λ (−0.06815 for "falkner_skan", −0.09 for "white") anywhere on [0, x].
    Validation: flat plate θ = √(0.45νx/U) exactly; diffuser λ closed form; cylinder λ = closed form of :func:`thwaites_cylinder_closed_form`.
    Label: analytic.
    """
    rho = float(params.pop("rho", 1.0))
    x = float(x)
    of = outer_flow(kind, **params)
    s = np.linspace(0.0, 1.0, 4001) ** 2
    xs = x * s
    th0 = float(theta0)
    if kind == "wedge":
        n, a = of.params["n"], of.params["a"]
        xs = x * (1e-6 + (1.0 - 1e-6) * s)
        th0 = float(np.sqrt(0.45 * nu * xs[0] ** (1.0 - n) / (a * (5.0 * n + 1.0))))
    r = thwaites(xs, of, nu, theta0=th0, closure=closure, rho=rho, stop_at_separation=False)
    i = -1
    return dict(theta=float(r["theta"][i]), delta_star=float(r["delta_star"][i]), tau0=float(r["tau0"][i]), cf=float(r["cf"][i]),
                lam=float(r["lam"][i]), H=float(r["H"][i]), l=float(r["l"][i]), x_sep=r["x_sep"])


def _cylinder_I5(phi):
    """∫₀^φ sin⁵ψ dψ = 8/15 − cos φ + (2/3)cos³φ − (1/5)cos⁵φ  (closed form; substitution c = cos ψ)."""
    c = np.cos(_F(phi))
    return 8.0 / 15.0 - c + 2.0 / 3.0 * c ** 3 - 0.2 * c ** 5


def thwaites_cylinder_closed_form(phi):
    """Thwaites' λ on the ideal cylinder flow U_e = 2U sin φ, x = aφ (φ from the FORWARD stagnation point) — closed form.

    Derivation: ∫₀^φ sin⁵ = 8/15 − cos φ + (2/3)cos³φ − (1/5)cos⁵φ ≡ F(φ), so by Eq. (9.50) θ² = 0.45(νa/2U)F/sin⁶φ and
    λ = θ²U_e′/ν = 0.45 F cos φ/sin⁶φ  [–] (independent of U, a, ν; λ(0) = 0.45/6 = 0.075; λ(90°) = 0).  phi [rad] in (0, π).
    Book: §9.6 Eq. (9.50) (Exercise 9.21 named only).  Label: analytic."""
    p = _F(phi)
    return _S(0.45 * _cylinder_I5(p) * np.cos(p) / np.sin(p) ** 6)


def _cylinder_lambda_numeric(phi: float) -> float:
    """λ(φ) by quadrature of Eq. (9.50) on U_e = 2 sin ψ (U = a = 1): λ = 0.45 (∫₀^φ U_e⁵dψ) U_e′/U_e⁶ (the ν, U, a cancel)."""
    if phi < 1e-6:
        return 0.45 / 6.0  # stagnation-point limit
    I = quad(lambda t: (2.0 * np.sin(t)) ** 5, 0.0, phi, epsabs=1e-14, epsrel=1e-13)[0]
    return float(0.45 * I * (2.0 * np.cos(phi)) / (2.0 * np.sin(phi)) ** 6)  # Eq. (9.50), (9.44)


def thwaites_cylinder(phi_deg=None, closure=None, U: float = 1.0, a: float = 1.0, nu: float = 1e-6, n: int = 400,
                      phi_max_deg: float = 170.0):
    """Thwaites' method on the ideal-flow cylinder U_e = 2U sin(x/a); angles ``phi_deg`` from the FORWARD stagnation point.

    Book: §9.6 Eqs. (9.44), (9.50).
    * ``phi_deg`` given (float or array [deg]): the NUMERICAL route (quadrature of ∫U_e⁵ by ``scipy.integrate.quad``, not the closed form) to
      λ(φ) [–], for parity with :func:`thwaites_cylinder_closed_form` (equal to 1e-12).  ``closure`` is accepted for signature symmetry only:
      λ does not depend on the l, H closure.
    * ``phi_deg=None``: the full :func:`thwaites` run on ``n`` points up to ``phi_max_deg`` for the given U [m/s], a [m], nu [m²/s] and ``closure``
      (default "falkner_skan"); returns its dict plus phi_deg and phi_sep_deg (first crossing of the closure's λ criterion), phi_sep_l0_deg.
    Label: converged (vs the closed form)."""
    if phi_deg is not None:
        ph = np.deg2rad(_F(phi_deg))
        vals = np.array([_cylinder_lambda_numeric(float(p_)) for p_ in np.atleast_1d(ph)]).reshape(np.shape(ph))
        return _S(vals)
    closure = "falkner_skan" if closure is None else closure
    of = outer_flow("cylinder", U=U, a=a)
    phi = np.linspace(0.0, np.deg2rad(phi_max_deg), int(n))
    r = thwaites(a * phi, of, nu, closure=closure, stop_at_separation=False)
    r["phi_deg"] = np.rad2deg(phi)
    r["phi_sep_deg"] = None if r["x_sep"] is None else float(np.rad2deg(r["x_sep"] / a))
    r["phi_sep_l0_deg"] = None if r["x_sep_l0"] is None else float(np.rad2deg(r["x_sep_l0"] / a))
    return r


def thwaites_cylinder_separation(lam_sep: float = LAMBDA_SEP_BOOK) -> float:
    """Separation angle φ_sep [deg from the forward stagnation point] where Thwaites' λ on the ideal cylinder flow falls to ``lam_sep``.

    Book: §9.6 Eq. (9.50) with U_e = 2U sin φ (Exercise 9.21 named only).  Root of :func:`thwaites_cylinder_closed_form` (λ decreases monotonically from 0
    at 90° to −∞ at 180°) by ``brentq``; independent of U, a, ν.  ``lam_sep`` [–]: −0.09 (the book's criterion) gives 103.1°, −0.06815
    (exact Falkner–Skan zero shear, :data:`LAMBDA_SEP_FS`) gives 100.9°.  These are Thwaites' predictions with the ideal-flow U_e; the measured
    subcritical separation angle is ≈ 82° (rounded experimental value, see core.bluff_body), so the method is only a rough guide on a bluff body.
    Label: analytic."""
    if lam_sep >= 0.0:
        raise ValueError("lam_sep must be negative (adverse gradient, φ > 90°)")
    g = lambda p: float(thwaites_cylinder_closed_form(np.deg2rad(p))) - lam_sep  # noqa: E731
    return float(brentq(g, 90.0, 179.0, xtol=1e-12))


# ======================================================================================================================
# §9.7 wall curvature, inflection, separation
# ======================================================================================================================


def wall_curvature(dpdx, mu):
    """(∂²u/∂y²)_wall = (dp/dx)/μ  [1/(m·s)]: Eq. (9.9) evaluated at the wall where u = v = 0.  Book: §9.7 (unnumbered, before (9.51)).
    dp/dx > 0 (decelerating stream) gives positive curvature — Eq. (9.52); dp/dx < 0 negative — Eq. (9.51)."""
    return _S(_F(dpdx) / mu)


def profile_inflection(y, u):
    """y [m] of the inflection point (∂²u/∂y² changes sign) of a sampled profile, or None; y[0] if the curvature is ≥ 0 at the wall
    only at a point (Blasius: at the wall).  Book: §9.7, Eqs. (9.51)–(9.52).  Uses second differences of a cubic spline; label analytic."""
    y, u = _F(y), _F(u)
    cs = CubicSpline(y, u)
    yy = np.linspace(y[0], y[-1], 20 * len(y))
    raw = cs(yy, 2)
    scale = np.max(np.abs(raw[: int(0.8 * len(raw))]))  # ignore spline end effects
    c2 = np.where(np.abs(raw) < 1e-3 * scale, 0.0, raw)  # spline noise where the curvature is (nearly) zero
    if c2[0] == 0.0 and np.all(c2[1:] <= 0):
        return float(y[0])  # zero curvature at the wall, negative above (Blasius)
    nz = np.where(c2 != 0.0)[0]
    flips = np.where(np.sign(c2[nz][:-1]) * np.sign(c2[nz][1:]) < 0)[0]
    if flips.size == 0:
        return None
    i, j = nz[flips[0]], nz[flips[0] + 1]
    return float(yy[i] - raw[i] * (yy[j] - yy[i]) / (raw[j] - raw[i]))


def separation_point(x, tau0):
    """First x where the wall shear changes sign (τ₀ = 0, (∂u/∂y)_wall = 0), by linear interpolation; None if attached.
    Book: §9.7 (definition of separation).  x [m]; tau0 [Pa].  Label: analytic."""
    x, t = _F(x), _F(tau0)
    idx = np.where((t[:-1] > 0) & (t[1:] <= 0))[0]
    if idx.size == 0:
        return None
    i = idx[0]
    return float(x[i] + t[i] * (x[i + 1] - x[i]) / (t[i] - t[i + 1]))


# ======================================================================================================================
# parabolic marching (von Mises)
# ======================================================================================================================


def march_boundary_layer(Ue, x_grid, nu: float, u_inlet=None, ny: int = 400, y_max_factor: float = 12.0, picard: int = 4,
                         rho: float = 1.0, fast: bool = False, order: int = 2) -> dict:
    """Solve the boundary-layer equations (9.9), (9.18) by marching downstream in x (x plays the role of time — the parabolic character).

    Book: §9.1 (parabolic character, Eq. (9.15) inlet profile), §9.3; the book gives no scheme (Ch. 10).
    Method (OURS): von Mises variables (x, ψ), w = u²:  ∂w/∂x = 2U_eU_e′ + ν√w ∂²w/∂ψ²  (from u u_x|_ψ = U_eU_e′ + νu ∂_ψ(u∂_ψu)).
    ψ = ψ_max σ² (σ uniform on [0,1]) turns the wall singularity w ~ ψ + ψ^{3/2} into a smooth polynomial in σ;  variable-step BDF2
    in x (second order, L-stable; the first step is backward Euler; ``order=1`` gives backward Euler throughout — first order; Crank–Nicolson
    was tried and oscillates when the step is large against the wall-layer diffusion time), central differences in σ, Picard iteration for
    the coefficient √w, tridiagonal solve; τ₀ = μ w_ψ/2 at the wall from the exact two-term extrapolation.  Marching stops when τ₀ ≤ 0
    (separation; the Goldstein singularity follows).
    Parameters: Ue OuterFlow/callable; x_grid [m] increasing, x_grid[0] > 0 (start); nu [m²/s]; u_inlet None (local Falkner–Skan
    profile with m = x U_e′/U_e at x_grid[0]; ValueError if m is below the fold −0.0904) | callable u(y) | (y, u) arrays [m, m/s]; ny nodes in σ; y_max_factor: ψ_max =
    factor·√(νU_max x_end) [m²/s]; picard iterations per step; rho [kg/m³]; fast: ny ≤ 150 (notebook FAST mode); order 2 (BDF2, default) or 1
    (backward Euler) time-like discretisation in x.
    Returns dict(x, psi [m²/s], sigma, y [m] (n_x×ny), u, v [m/s] (n_x×ny), tau0 [Pa], separated, x_sep).
    Validation: V5/V3 τ₀ vs Blasius (error ∝ Δx² for order 2, ∝ Δx for order 1); memory of the inlet forgotten on a flat plate.
    Label: converged.
    """
    x = _F(x_grid)
    if fast:
        ny = min(int(ny), 150)  # notebook FAST mode: coarser wall-normal grid (τ₀ error rises to ~1e-2)
    of = _as_outer(Ue, x)
    Umax = float(np.max(of.Ue(x)))
    psi_max = y_max_factor * np.sqrt(nu * Umax * x[-1])
    N = int(ny)
    sig = np.linspace(0.0, 1.0, N + 1)
    hs = sig[1] - sig[0]
    psi = psi_max * sig ** 2
    Ue0 = float(of.Ue(x[0]))
    if u_inlet is None:
        m0 = x[0] * float(of.dUe(x[0])) / Ue0
        eta = np.linspace(0, 12, 4001)
        if not (m0 > _FOLD):
            raise ValueError(f"march_boundary_layer: the local exponent m0 = x U_e′/U_e = {m0:.4g} at x_grid[0] is at or below the Falkner–Skan "
                             f"separation value {_FOLD}, so no attached similarity inlet profile exists; pass u_inlet (a callable u(y) or (y, u) arrays, "
                             "e.g. the Blasius profile) or start where the layer is attached.")
        d = falkner_skan(m0, 12.0, n=4001)
        f, fp = d["f"], d["fp"]
        # w = u² is interpolated against f (both ∝ η² at the wall, so the interpolant is smooth; interpolating f′ against f would be a √f cusp)
        w = Ue0 ** 2 * np.interp(psi / np.sqrt(nu * x[0] * Ue0), f, fp ** 2)
    else:
        if callable(u_inlet):
            yy = np.linspace(0, 12 * np.sqrt(nu * x[0] / Ue0), 4001)
            uu = _F(u_inlet(yy))
        else:
            yy, uu = _F(u_inlet[0]), _F(u_inlet[1])
        ps = cumulative_trapezoid(uu, yy, initial=0.0)
        w = np.interp(psi, ps, uu ** 2, right=uu[-1] ** 2)
    w[0] = 0.0
    w[-1] = Ue0 ** 2
    coef0 = nu / (4.0 * psi_max ** 2)
    # DEVIATION: series-consistent wall-region operator.  Near the wall w = Aσ² + Bσ³ + Cσ⁴ … with B ∝ dp/dx ≠ 0; the central stencil of
    # L[w] = w_σσ − w_σ/σ is exact on σ² but returns 3jh − h/j instead of 3jh on σ³ at node j (1/3 too small at j = 1).  Multiplying L by
    # a_j = 1/(1 − 1/(3j²)) makes it exact on σ² AND σ³ at every node (a_j → 1 + O(h²/σ²)), which restores second order in Δσ for τ₀.
    jj = np.arange(1, N, dtype=float)
    a_j = 1.0 / (1.0 - 1.0 / (3.0 * jj ** 2))
    W, TAU, XS = [w.copy()], [], []

    def wall_shear(wv):
        A = (8.0 * wv[1] - wv[2]) / (4.0 * hs ** 2)  # w ≈ Aσ² + Bσ³ ;  τ0/μ = A/(2ψ_max)
        return A, A / (2.0 * psi_max)

    A0, _ = wall_shear(w)
    TAU.append(rho * nu * A0 / (2.0 * psi_max))
    x_sep = None
    xs_done = [x[0]]
    sj = sig[1:-1]
    for i in range(1, len(x)):
        h = x[i] - x[i - 1]
        Uen, dUen = float(of.Ue(x[i])), float(of.dUe(x[i]))
        wn = W[-1].copy()
        wit = wn.copy()
        wit[-1] = Uen ** 2
        Uo, dUo = float(of.Ue(x[i - 1])), float(of.dUe(x[i - 1]))
        th_ = 1.0  # implicit in x
        if order == 2 and len(W) >= 2:  # variable-step BDF2: [(1+2ω)/(1+ω) w_n − (1+ω) w_{n−1} + ω²/(1+ω) w_{n−2}]/h = RHS(w_n)
            om = h / (x[i - 1] - x[i - 2])
            a0 = (1.0 + 2.0 * om) / (1.0 + om)
            w_old = (1.0 + om) * W[-1] - om ** 2 / (1.0 + om) * W[-2]
        else:  # backward Euler
            a0, w_old = 1.0, wn
        kap_o = a_j * coef0 * np.sqrt(np.maximum(wn[1:-1], 0.0)) / sj ** 2
        Lw_o = (wn[2:] - 2 * wn[1:-1] + wn[:-2]) / hs ** 2 - (wn[2:] - wn[:-2]) / (2 * hs * sj)
        explicit = (1.0 - th_) * (2.0 * Uo * dUo + kap_o * Lw_o)
        for _ in range(picard):
            kap = a_j * coef0 * np.sqrt(np.maximum(wit[1:-1], 0.0)) / sj ** 2
            lower = -h * th_ * kap * (1.0 / hs ** 2 + 1.0 / (2.0 * hs * sj))
            diag = a0 + h * th_ * kap * 2.0 / hs ** 2
            upper = -h * th_ * kap * (1.0 / hs ** 2 - 1.0 / (2.0 * hs * sj))
            rhs = w_old[1:-1] + h * explicit + th_ * 2.0 * h * Uen * dUen
            rhs[-1] -= upper[-1] * Uen ** 2
            ab = np.zeros((3, N - 1))
            ab[0, 1:] = upper[:-1]
            ab[1, :] = diag
            ab[2, :-1] = lower[1:]
            wnew = np.empty_like(wit)
            wnew[0], wnew[-1] = 0.0, Uen ** 2
            wnew[1:-1] = solve_banded((1, 1), ab, rhs)
            wit = wnew
        A, _ = wall_shear(wit)
        if A <= 0.0:
            x_sep = float(x[i - 1] + TAU[-1] / (TAU[-1] - rho * nu * A / (2 * psi_max)) * h) if TAU[-1] > 0 else float(x[i - 1])
            break
        W.append(wit)
        TAU.append(rho * nu * A / (2.0 * psi_max))
        xs_done.append(x[i])
    W = np.array(W)
    xs = np.array(xs_done)
    Uw = np.sqrt(np.maximum(W, 0.0))
    # y(σ) = ∫ dψ/u = ∫ 2ψ_max σ dσ / u ;  limit at the wall  2ψ_max/√A
    Y = np.zeros_like(W)
    for k in range(len(xs)):
        A = (8.0 * W[k, 1] - W[k, 2]) / (4.0 * hs ** 2)
        integrand = np.empty(N + 1)
        integrand[1:] = 2.0 * psi_max * sig[1:] / Uw[k, 1:]
        integrand[0] = 2.0 * psi_max / np.sqrt(A)
        Y[k] = cumulative_trapezoid(integrand, sig, initial=0.0)
    V = np.zeros_like(Y)
    if len(xs) > 1:
        dydx = np.gradient(Y, xs, axis=0, edge_order=1)
        V = Uw * dydx  # v = u ∂y/∂x|_ψ  (= −ψ_x|_y)
    return dict(x=xs, psi=psi, sigma=sig, y=Y, u=Uw, v=V, tau0=np.array(TAU), separated=x_sep is not None, x_sep=x_sep)


# ======================================================================================================================
# §9.7 transition and flat-plate drag
# ======================================================================================================================


def transition_state(Re_x, Re_cr: float = 5e5, Re_turb: float | None = None) -> str:
    """Regime of a flat-plate layer at local Reynolds number Re_x = Ux/ν: "laminar" (< Re_cr), "transitional" (Re_cr ≤ Re_x < Re_turb), "turbulent".

    Book: §9.7 (Fig. 9.10/9.11): transition begins near Re_x ≈ 5×10⁵ (Re_cr ≈ 10⁶ "within a factor of five" — the book's roundings, slip R9) and
    the layer is fully turbulent some factor further on.  The instability mechanism is Ch. 11.  Both thresholds are arguments (experimental,
    rounded): Re_cr [–] default 5×10⁵; Re_turb [–] default 10·Re_cr (the transitional band is one decade wide).
    Label: qualitative."""
    Re_x = float(Re_x)
    Re_turb = 10.0 * Re_cr if Re_turb is None else Re_turb
    if Re_x < Re_cr:
        return "laminar"
    return "transitional" if Re_x < Re_turb else "turbulent"


def plate_drag_coefficient(Re_L, regime: str = "laminar", Re_tr: float = 5e5, sides: int = 1, turbulent_coeff: float = 0.074):
    """Flat-plate drag coefficient C_D = F/(½ρU²L) per wetted face vs Re_L = UL/ν.

    Book: §9.3 Eq. (9.33) (laminar, 1.328/√Re_L) and Fig. 9.11 (the laminar line, a turbulent line, a transition curve).
    "turbulent": c Re_L^{−1/5} with c = ``turbulent_coeff`` = 0.074 (Prandtl's one-seventh-power-law fit, valid 5×10⁵ < Re_L < 10⁷, NOT in the book;
    Schlichting's refit is c = 0.0725, a 2 % difference — form confirmed by a web search snippet of Schlichting's law, coefficients not re-read
    in a primary source: status "form benchmark, coefficient unverified").  "mixed": laminar up to Re_tr, then C_turb(Re_L) − (Re_tr/Re_L)[C_turb(Re_tr) − C_lam(Re_tr)]
    (turbulent layer starts at Re_tr with the laminar momentum deficit; = 0.074Re^{−1/5} − 1742/Re for Re_tr = 5×10⁵).
    Parameters: Re_L [–]; regime; Re_tr [–] (book value 5×10⁵, rounded); sides wetted faces.
    Label: analytic (laminar); turbulent: form benchmark, coefficient secondary-sourced."""
    R = _F(Re_L)
    lam = 4.0 * _blasius().fpp0 / np.sqrt(R)
    tur = turbulent_coeff * R ** -0.2
    if regime == "laminar":
        c = lam
    elif regime == "turbulent":
        c = tur
    elif regime == "mixed":
        c_tr = turbulent_coeff * Re_tr ** -0.2 - 4.0 * _blasius().fpp0 / np.sqrt(Re_tr)
        c = np.where(R <= Re_tr, lam, tur - Re_tr / R * c_tr)
    else:
        raise ValueError('regime must be "laminar", "turbulent" or "mixed"')
    return _S(sides * c)


__all__ = [n for n in dir() if not n.startswith("_") and n not in ("annotations", "Callable", "dataclass", "functools")]
