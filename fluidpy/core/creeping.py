"""Creeping (Stokes) flow: the Stokes equations, Stokes' sphere (stream function, velocity, pressure, surface stresses,
drag 6πμaU and its ⅓/⅔ split), terminal velocity and Millikan's balance, the far-field breakdown Re·r/a, Oseen's
correction (stream function, velocity, C_D) and the low-Re drag laws.

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e (2012), Ch. 8, §8.6, Eqs. (8.39)–(8.53) (page images
chapters/pages/ch08/p366–p373). Reused by Ch. 13 (sediment and cloud-droplet settling, aerosols), Ch. 16
(micro-swimmers, cells), Ch. 12 (particle-laden flows), Ch. 4 C15 (the low-Re branch of the drag curve, now derived).

Conventions (analysis ch08 §9 R5, R12, R13, R17–R21)
* Spherical (r, θ, φ) with **θ measured from the downstream +x axis** (the direction of the stream U e_x): rear
  stagnation point θ = 0, front θ = π (as ch06 (6.91), not ch03's θ from +z).
* ``frame="body"``: sphere at rest in a stream U e_x (the book's solution); ``frame="fluid"``: the uniform stream
  subtracted (sphere moving to −x through fluid at rest, Figs. 8.19–8.20).
* **Re = 2aU/ν (diameter)** in (8.52), (8.53) and the Oseen C_D; Proudman–Pearson and Wikipedia use the radius
  (Re_a = Re/2: the book's 3/16 is the literature's 3/8).
* Printed slips (coded corrected): R12 the minimum p − p∞ is **−**3μU/2a (rear); R13 Oseen's equation needs
  **−**∂p/∂x_i; (8.44) is the square of the operator E², not the biharmonic.
* Sphere fields return NaN inside the sphere (r < a). g defaults to G0 = 9.80665 m/s² (settling physics).
"""
from __future__ import annotations

import functools
import warnings
from typing import Callable

import numpy as np
import sympy as sp

from . import _stencil as st
from ._util import as_scalar_if_0d, require_positive
from .thermo import G0

__all__ = ["stokes_residual", "E2", "E4_residual", "stokes_sphere_sympy", "stokes_sphere_streamfunction",
           "stokes_sphere_velocity", "stokes_sphere_velocity_xyz", "stokes_sphere_pressure",
           "stokes_sphere_surface_stresses", "stokes_drag", "stokes_drag_running", "side_line_speed",
           "sphere_drag_quadrature", "stokes_drag_coefficient",
           "oseen_drag_coefficient", "proudman_pearson_drag_coefficient", "terminal_velocity",
           "radius_from_terminal_velocity", "millikan_charge", "settling_state", "inertia_viscous_ratio",
           "oseen_streamfunction", "oseen_velocity"]

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d


def _outside(r, a):
    return _F(r) >= float(a) * (1.0 - 1e-12)


# ======================================================================================================================
# (8.43) the Stokes equations; the E² operator (8.44)
# ======================================================================================================================
def stokes_residual(u_fn: Callable, p_fn: Callable, x, mu: float = 1e-3, h: float = 1e-4):
    """Residual ∇p − μ∇²u of the Stokes (creeping-flow) equations at points x.

    Book: §8.6, Eq. (8.43) ∇p = μ∇²u (from (8.42) as Re → 0, with the viscous pressure scale μU/L).

    Parameters
    ----------
    u_fn : callable x ↦ u, x of shape (3,) or (3, N), returning (3,) / (3, N) [m/s].
    p_fn : callable x ↦ p [Pa].
    x : points [m], shape (3,) or (3, N).
    mu : μ [Pa s];  h : central-difference step [m] (second order).

    Returns
    -------
    res : ndarray (3,) or (3, N) [Pa/m].

    Validation — tests/test_ch08.py: test_stokes_residual_V1_sphere_field_and_ideal_flow_fails,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V1 ≈ 0 (relative 1e-6) for the Stokes sphere; ≠ 0 for the ideal-flow sphere.
    Label: analytic.
    """
    U = lambda X, T: u_fn(X)  # noqa: E731
    P = lambda X, T: p_fn(X)  # noqa: E731
    x_ = _F(x)
    d = x_.shape[0]
    return st.grad(P, x_, 0.0, h) - float(mu) * st.laplacian(U, x_, 0.0, h, (d,))  # Eq. (8.43)


def E2(psi, r: sp.Symbol, theta: sp.Symbol):
    """The Stokes operator E²ψ = ∂²ψ/∂r² + (sin θ/r²) ∂/∂θ[(1/sin θ) ∂ψ/∂θ] (sympy).

    Book: §8.6 (ω_φ = −E²ψ/(r sin θ), the display before (8.44); the operator of ch06 (6.77) in spherical form).
    Parameters: psi sympy expression in r, θ; r, theta sympy symbols. Returns sympy expression. Label: symbolic.
    """
    return sp.diff(psi, r, 2) + sp.sin(theta) / r ** 2 * sp.diff(sp.diff(psi, theta) / sp.sin(theta), theta)


def E4_residual(psi, r: sp.Symbol, theta: sp.Symbol):
    """E²(E²ψ), simplified — zero for Stokes-flow stream functions.

    Book: §8.6, Eq. (8.44) (the square of the operator, footnote 2 — not the biharmonic). Returns sympy expression.
    Validation — tests/test_ch08.py: test_stokes_sphere_sympy_V2_engine,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V2 zero for (8.48); nonzero for the biharmonic misreading. Label: symbolic.
    """
    return sp.simplify(E2(E2(psi, r, theta), r, theta))  # Eq. (8.44)


def stokes_sphere_sympy() -> dict:
    """Stokes' sphere solution derived step by step in sympy (the analysis I33 engine; cached — a fresh dict of the
    same sympy objects is returned on every call).

    See :func:`_stokes_sphere_sympy` for the steps and keys. Book: §8.6, Eqs. (8.44)–(8.51). Label: symbolic.
    """
    return dict(_stokes_sphere_sympy())


@functools.lru_cache(maxsize=1)
def _stokes_sphere_sympy() -> dict:
    """Stokes' sphere solution derived step by step in sympy (the analysis I33 engine).

    Book: §8.6, Eqs. (8.44)–(8.51): separable ψ = f(r) sin²θ in (8.44) gives f⁗ − 4f″/r² + 8f′/r³ − 8f/r⁴ = 0, roots
    r⁴, r², r, 1/r; (8.47) ⇒ A = 0, B = U/2; (8.45)–(8.46) ⇒ C = −3Ua/4, D = Ua³/4; (8.48)–(8.49); p from ∇p = μ∇²u
    (8.50) (the integration the book skips); surface stresses; the drag parts (Exercise 8.35) and (8.51).

    Returns
    -------
    dict of sympy objects (design Part C 3.16 keys and longer aliases): symbols, E2_of_fsin2 (E²(f sin²θ)/sin²θ),
    f_ode, roots ([−1, 1, 2, 4]), f_general, constants (A, B, C, D), psi, u_r, u_theta, omega_phi, E2psi, E4psi
    (= E4_psi, 0), biharmonic_psi (scalar ∇⁴ψ, ≠ 0), curl_omega ((r, θ) components of ∇×(ω_φ e_φ)),
    curlcurl_identity (∇p + μ∇×ω, 0), p (p − p∞), p_check_r, p_check_theta (= momentum_residual_r/theta, 0),
    divergence (0), sigma_rr_a (= −p on r = a), sigma_rr_viscous_a (2μ∂u_r/∂r at a, 0), sigma_rtheta_a, t_x_a,
    D_pressure (= drag_pressure, 2πμaU), D_friction (= drag_friction, 4πμaU), drag (6πμaU).

    Validation — tests/test_ch08.py: (via stokes_sphere_sympy) test_stokes_sphere_sympy_V2_engine,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V2 every residual simplifies to 0; drag parts 2πμaU and 4πμaU. Label: symbolic.
    """
    r, th, a, U, mu = sp.symbols("r theta a U mu", positive=True)
    f = sp.Function("f")
    A, B, C, D = sp.symbols("A B C D")
    E2f = sp.simplify(E2(f(r) * sp.sin(th) ** 2, r, th) / sp.sin(th) ** 2)
    ode = sp.expand(sp.simplify(E2(E2(f(r) * sp.sin(th) ** 2, r, th), r, th) / sp.sin(th) ** 2))
    lam = sp.Symbol("lambda")
    char = sp.simplify(ode.subs(f(r), r ** lam).doit() / r ** (lam - 4))
    roots = sorted(sp.solve(sp.expand(char), lam))
    fg = A * r ** 4 + B * r ** 2 + C * r + D / r
    cons = {A: 0, B: U / 2}
    f1 = fg.subs(cons)
    sol = sp.solve([f1.subs(r, a), sp.diff(f1, r).subs(r, a)], [C, D], dict=True)[0]
    cons.update(sol)
    psi = sp.factor(fg.subs(cons)) * sp.sin(th) ** 2
    ur = sp.simplify(sp.diff(psi, th) / (r ** 2 * sp.sin(th)))  # (6.83)
    ut = sp.simplify(-sp.diff(psi, r) / (r * sp.sin(th)))
    om = sp.simplify(-E2(psi, r, th) / (r * sp.sin(th)))
    # spherical vector Laplacian (axisymmetric, no swirl) and ∇p
    lap_r = (sp.diff(r ** 2 * sp.diff(ur, r), r) / r ** 2 + sp.diff(sp.sin(th) * sp.diff(ur, th), th) / (r ** 2 * sp.sin(th))
             - 2 * ur / r ** 2 - 2 * sp.diff(ut * sp.sin(th), th) / (r ** 2 * sp.sin(th)))
    lap_t = (sp.diff(r ** 2 * sp.diff(ut, r), r) / r ** 2 + sp.diff(sp.sin(th) * sp.diff(ut, th), th) / (r ** 2 * sp.sin(th))
             + 2 * sp.diff(ur, th) / r ** 2 - ut / (r ** 2 * sp.sin(th) ** 2))
    # integrate ∂p/∂r = μ lap_r from ∞ (p → p∞ = 0) and check the θ component
    p = sp.simplify(-sp.integrate(sp.simplify(mu * lap_r), (r, r, sp.oo)))
    res_r = sp.simplify(sp.diff(p, r) - mu * lap_r)
    res_t = sp.simplify(sp.diff(p, th) / r - mu * lap_t)
    div = sp.simplify(sp.diff(r ** 2 * ur, r) / r ** 2 + sp.diff(sp.sin(th) * ut, th) / (r * sp.sin(th)))
    s_rr = sp.simplify((-p + 2 * mu * sp.diff(ur, r)).subs(r, a))
    s_rt = sp.simplify((mu * (r * sp.diff(ut / r, r) + sp.diff(ur, th) / r)).subs(r, a))
    tx = sp.simplify(s_rr * sp.cos(th) - s_rt * sp.sin(th))
    dA = 2 * sp.pi * a ** 2 * sp.sin(th)
    Dp = sp.simplify(sp.integrate((-p.subs(r, a)) * sp.cos(th) * dA, (th, 0, sp.pi)))
    Df = sp.simplify(sp.integrate((-s_rt * sp.sin(th)) * dA, (th, 0, sp.pi)))
    # scalar Laplacian applied twice (the biharmonic misreading of (8.44), footnote 2) — must NOT vanish
    lap_s = lambda q: sp.diff(r ** 2 * sp.diff(q, r), r) / r ** 2 + sp.diff(sp.sin(th) * sp.diff(q, th), th) / (  # noqa: E731
        r ** 2 * sp.sin(th))
    bih = sp.simplify(lap_s(lap_s(psi)))
    # ∇×(ω_φ e_φ) in spherical components (r, θ): (1/(r sin θ))∂(sin θ ω)/∂θ, −(1/r)∂(r ω)/∂r
    cw_r = sp.simplify(sp.diff(sp.sin(th) * om, th) / (r * sp.sin(th)))
    cw_t = sp.simplify(-sp.diff(r * om, r) / r)
    # ∇p = μ∇²u = −μ∇×ω (the curl-curl form, ∇·u = 0): residuals of ∇p + μ∇×ω
    cc = sp.simplify(sp.diff(p, r) + mu * cw_r) + sp.simplify(sp.diff(p, th) / r + mu * cw_t)
    s_rr_visc = sp.simplify((2 * mu * sp.diff(ur, r)).subs(r, a))
    E4 = E4_residual(psi, r, th)
    return dict(symbols=(r, th, a, U, mu), E2_of_fsin2=E2f, f_ode=ode, roots=roots, f_general=fg, constants=cons,
                psi=psi, u_r=ur, u_theta=ut, omega_phi=om, E2psi=sp.simplify(E2(psi, r, th)), E4_psi=E4, E4psi=E4,
                biharmonic_psi=bih, curl_omega=(cw_r, cw_t), curlcurl_identity=cc, p=p,
                momentum_residual_r=res_r, momentum_residual_theta=res_t, p_check_r=res_r, p_check_theta=res_t,
                divergence=div, sigma_rr_a=s_rr, sigma_rr_viscous_a=s_rr_visc, sigma_rtheta_a=s_rt, t_x_a=tx,
                drag_pressure=Dp, drag_friction=Df, D_pressure=Dp, D_friction=Df, drag=sp.simplify(Dp + Df))


# ======================================================================================================================
# Stokes' sphere (8.48)–(8.51)
# ======================================================================================================================
def _frame(frame: str) -> bool:
    if frame not in ("body", "fluid"):
        raise ValueError('frame must be "body" or "fluid"')
    return frame == "fluid"


def stokes_sphere_streamfunction(r, theta, U: float = 1.0, a: float = 1.0, frame: str = "body"):
    """Stokes' stream function for creeping flow past a sphere.

    Book: §8.6, Eq. (8.48) ψ = Ur² sin²θ(½ − 3a/4r + a³/4r³) (body frame); fluid frame (sphere moving to −x) the
    display after (8.52): ψ = Ur² sin²θ(−3a/4r + a³/4r³).

    Parameters
    ----------
    r : [m] (≥ a; NaN inside);  theta : θ [rad] from the downstream +x axis;  U : stream speed [m/s];  a : radius [m];
    frame : "body" | "fluid".

    Returns
    -------
    psi : Stokes stream function [m³/s] (u_r = ψ_θ/(r² sin θ), u_θ = −ψ_r/(r sin θ), (6.83)).

    Validation — tests/test_ch08.py: test_stokes_sphere_velocity_V1_walls_far_field_divergence,
      test_stokes_sphere_V7_frames_and_fore_aft_symmetry, test_oseen_V1_velocity_axis_wake_and_no_slip_order,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable,
      test_book_V6_section_8_4_to_8_6_numbers_and_slips.
    Checks: V2 E⁴ψ = 0, BCs (8.45)–(8.47). Label: symbolic, analytic.
    """
    r_, t_ = np.broadcast_arrays(_F(r), _F(theta))
    a, U = float(a), float(U)
    lead = 0.0 if _frame(frame) else 0.5
    with np.errstate(divide="ignore", invalid="ignore"):
        psi = U * r_ ** 2 * np.sin(t_) ** 2 * (lead - 3.0 * a / (4.0 * r_) + a ** 3 / (4.0 * r_ ** 3))  # Eq. (8.48)
    return _S(np.where(_outside(r_, a), psi, np.nan))


def stokes_sphere_velocity(r, theta, U: float = 1.0, a: float = 1.0, frame: str = "body"):
    """Velocity components of Stokes flow past a sphere.

    Book: §8.6, Eq. (8.49): u_r = U cos θ(1 − 3a/2r + a³/2r³), u_θ = −U sin θ(1 − 3a/4r − a³/4r³) (body frame);
    fluid frame: the stream U e_x (u_r = U cos θ, u_θ = −U sin θ) subtracted.

    Parameters
    ----------
    r : [m];  theta : [rad] from +x;  U : [m/s];  a : [m];  frame : "body" | "fluid".

    Returns
    -------
    (u_r, u_theta) : [m/s] (NaN inside the sphere).

    Validation — tests/test_ch08.py: test_stokes_sphere_velocity_V1_walls_far_field_divergence,
      test_stokes_sphere_V7_frames_and_fore_aft_symmetry, test_stokes_drag_V1_quadrature_parts_and_slip_sphere,
      test_oseen_V1_velocity_axis_wake_and_no_slip_order,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable,
      test_book_V6_section_8_4_to_8_6_numbers_and_slips.
    Checks: V1 no slip at r = a, → U e_x far away, fore–aft symmetry of |u| in the fluid frame; V2 from
    (8.48) via (6.83), ∇·u = 0. Label: analytic, symbolic.
    """
    r_, t_ = np.broadcast_arrays(_F(r), _F(theta))
    a, U = float(a), float(U)
    c, s = np.cos(t_), np.sin(t_)
    with np.errstate(divide="ignore", invalid="ignore"):
        q = a / r_
        ur = U * c * (1.0 - 1.5 * q + 0.5 * q ** 3)  # Eq. (8.49)
        ut = -U * s * (1.0 - 0.75 * q - 0.25 * q ** 3)
    if _frame(frame):
        ur, ut = ur - U * c, ut + U * s
    m = _outside(r_, a)
    return _S(np.where(m, ur, np.nan)), _S(np.where(m, ut, np.nan))


def _spherical_from_xyz(x, y, z):
    X, Y, Z = np.broadcast_arrays(_F(x), _F(y), _F(z))
    rc = np.hypot(Y, Z)
    r = np.sqrt(X ** 2 + rc ** 2)
    th = np.arctan2(rc, X)
    with np.errstate(divide="ignore", invalid="ignore"):
        cy = np.where(rc > 0, Y / np.where(rc > 0, rc, 1.0), 1.0)
        cz = np.where(rc > 0, Z / np.where(rc > 0, rc, 1.0), 0.0)
    return r, th, cy, cz


def stokes_sphere_velocity_xyz(x, y, z, U: float = 1.0, a: float = 1.0, frame: str = "body"):
    """Cartesian velocity (u, v, w) of Stokes flow past a sphere centred at the origin, stream along +x.

    Book: §8.6, Eq. (8.49) with e_r = (cos θ, sin θ ĉ), e_θ = (−sin θ, cos θ ĉ), ĉ the unit vector from the x axis in the
    (y, z) plane. Parameters: x, y, z [m] (broadcast); U [m/s]; a [m]; frame. Returns (u, v, w) [m/s] (NaN inside).
    Validation — tests/test_ch08.py: test_stokes_residual_V1_sphere_field_and_ideal_flow_fails,
      test_stokes_sphere_velocity_V1_walls_far_field_divergence,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V1 ∇·u = 0 by finite differences; Stokes residual (8.43) ≈ 0. Label: analytic.
    """
    r, th, cy, cz = _spherical_from_xyz(x, y, z)
    ur, ut = stokes_sphere_velocity(r, th, U, a, frame)
    ur, ut = _F(ur), _F(ut)
    c, s = np.cos(th), np.sin(th)
    radial = ur * s + ut * c  # component along ĉ
    return _S(ur * c - ut * s), _S(radial * cy), _S(radial * cz)


def stokes_sphere_pressure(r, theta, U: float = 1.0, a: float = 1.0, mu: float = 1e-3, p_inf: float = 0.0):
    """Pressure of Stokes flow past a sphere, p − p∞ = −3μaU cos θ/(2r²).

    Book: §8.6, Eq. (8.50) (from ∇p = μ∇²u; the integration is ours). Maximum +3μU/2a at the front stagnation point
    θ = π, minimum **−**3μU/2a at the rear θ = 0 (the page prints the minimum without the minus, analysis §9 R12).
    Frame-independent (the same field seen from the moving sphere).

    Parameters
    ----------
    r : [m];  theta : [rad] from +x;  U : [m/s];  a : [m];  mu : [Pa s];  p_inf : [Pa].

    Returns
    -------
    p : [Pa] (NaN inside).

    Validation — tests/test_ch08.py: test_stokes_residual_V1_sphere_field_and_ideal_flow_fails,
      test_stokes_sphere_V7_frames_and_fore_aft_symmetry, test_stokes_drag_V1_quadrature_parts_and_slip_sphere,
      test_stokes_pressure_V1_extremes_and_printed_minimum,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable,
      test_book_V6_section_8_4_to_8_6_numbers_and_slips.
    Checks: V2 both components of ∇p = μ∇²u; V1 the extremes (wrong variant: +3μaU cos θ/2r²).
    Label: symbolic, analytic.
    """
    r_, t_ = np.broadcast_arrays(_F(r), _F(theta))
    with np.errstate(divide="ignore", invalid="ignore"):
        # (8.50) as printed; only the prose after it slips (the rear minimum is −3μU/2a, analysis §9 R12)
        p = float(p_inf) - 3.0 * float(mu) * float(a) * float(U) * np.cos(t_) / (2.0 * r_ ** 2)  # Eq. (8.50)
    return _S(np.where(_outside(r_, a), p, np.nan))


def stokes_sphere_surface_stresses(theta, U: float = 1.0, a: float = 1.0, mu: float = 1e-3, p_inf: float = 0.0):
    """Stresses on the surface of Stokes' sphere (Fig. 8.17 upper panel; Exercise 8.35), ours from (8.49)–(8.50).

    Book: §8.6, Fig. 8.17 and Exercise 8.35: on r = a the viscous normal stress 2μ∂u_r/∂r vanishes, so
    σ_rr = −p = −p∞ + (3μU/2a) cos θ; σ_rθ = μ[r ∂(u_θ/r)/∂r + (1/r)∂u_r/∂θ] = −(3μU/2a) sin θ; the x-traction on the
    sphere t_x = σ_rr cos θ − σ_rθ sin θ = −p∞ cos θ + 3μU/2a (uniform over the surface apart from p∞).

    Parameters
    ----------
    theta : [rad] from +x;  U : [m/s];  a : [m];  mu : [Pa s];  p_inf : [Pa].

    Returns
    -------
    (sigma_rr, sigma_rtheta, t_x) : [Pa].

    Validation — tests/test_ch08.py: test_stokes_drag_V1_quadrature_parts_and_slip_sphere,
      test_stokes_pressure_V1_extremes_and_printed_minimum,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V2 sympy from the fields; the integral of t_x over the sphere is 6πμaU. Label: symbolic.
    """
    t_ = _F(theta)
    k = 3.0 * float(mu) * float(U) / (2.0 * float(a))
    srr = -float(p_inf) + k * np.cos(t_)
    srt = -k * np.sin(t_)
    tx = srr * np.cos(t_) - srt * np.sin(t_)
    return _S(srr), _S(srt), _S(tx)


def stokes_drag(mu: float, a: float, U: float, parts: bool = False):
    """Stokes' drag on a sphere, D = 6πμaU: one third pressure drag (2πμaU), two thirds skin friction (4πμaU).

    Book: §8.6, Eq. (8.51) (Stokes' law of resistance; the split is stated, derived in Exercise 8.35).

    Parameters
    ----------
    mu : [Pa s];  a : radius [m];  U : speed [m/s];  parts : return dict(pressure, friction, total).

    Returns
    -------
    D : [N] (or dict of floats).

    Validation — tests/test_ch08.py: test_stokes_drag_V1_quadrature_parts_and_slip_sphere,
      test_stokes_law_V1_form_and_drag_coefficient,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable,
      test_book_V6_section_8_4_to_8_6_numbers_and_slips.
    Checks: V2 sympy surface integral; V3 Gauss–Legendre quadrature of the tractions. Label: symbolic,
    analytic.
    """
    D = 6.0 * np.pi * float(mu) * float(a) * float(U)  # Eq. (8.51)
    if parts:
        return dict(pressure=D / 3.0, friction=2.0 * D / 3.0, total=D)
    return D


def stokes_drag_running(theta, mu: float, a: float, U: float) -> dict:
    """Drag on Stokes' sphere collected from the rear stagnation point θ = 0 up to the polar angle θ (design D31).

    Book: §8.6, (8.50)–(8.51) and Exercise 8.35. Integrating the x-tractions of :func:`stokes_sphere_surface_stresses`
    over the cap 0 ≤ θ′ ≤ θ (area element 2πa² sin θ′ dθ′), ours in closed form: pressure πμaU(1 − cos³θ), friction
    πμaU(2 − 3cos θ + cos³θ); at θ = π they are 2πμaU and 4πμaU (⅓ and ⅔ of 6πμaU).

    Parameters
    ----------
    theta : [rad] (0 … π, from the downstream axis);  mu : [Pa s];  a : [m];  U : [m/s].

    Returns
    -------
    dict: pressure, friction, total [N] (floats or arrays like theta).

    Validation — tests/test_ch08.py: test_stokes_drag_V1_quadrature_parts_and_slip_sphere,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V1 quadrature of the tractions; V2 sympy; endpoints ⅓/⅔. Label: analytic.
    """
    c = np.cos(_F(theta))
    k = np.pi * float(mu) * float(a) * float(U)
    pr = k * (1.0 - c ** 3)
    fr = k * (2.0 - 3.0 * c + c ** 3)
    return dict(pressure=_S(pr), friction=_S(fr), total=_S(pr + fr))


def side_line_speed(r, U: float = 1.0, a: float = 1.0, model: str = "stokes", Re: float = 0.1, frame: str = "body"):
    """Speed |u| on the side line θ = π/2 (perpendicular to the stream) for three sphere flows (design C13 number).

    Book: §8.6 — ``"stokes"`` (8.49): u_θ = −U(1 − 3a/4r − a³/4r³); ``"ideal"`` the potential-flow sphere of §6.8
    (6.90): u_θ = −U(1 + a³/2r³); ``"oseen"`` (8.53) via :func:`oseen_velocity`. ``frame="fluid"`` subtracts the stream
    (at θ = π/2 that adds U to u_θ). Parameters: r [m] (≥ a); U [m/s]; a [m]; model; Re (Oseen only); frame.
    Returns |u| [m/s] (NaN inside).
    Validation — tests/test_ch08.py: test_stokes_sphere_V7_frames_and_fore_aft_symmetry,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V1 against the component functions. Label: analytic.
    """
    r_ = _F(r)
    fl = _frame(frame)
    if model == "stokes":
        ur, ut = stokes_sphere_velocity(r_, np.pi / 2, U, a, frame)
    elif model == "oseen":
        ur, ut = oseen_velocity(r_, np.pi / 2, U, a, Re, frame)
    elif model == "ideal":
        with np.errstate(divide="ignore", invalid="ignore"):
            ut = -float(U) * (1.0 + 0.5 * (float(a) / r_) ** 3) + (float(U) if fl else 0.0)
        ur = 0.0 * r_
        ut = np.where(_outside(r_, a), ut, np.nan)
    else:
        raise ValueError('model must be "stokes", "ideal" or "oseen"')
    return _S(np.hypot(_F(ur), _F(ut)))


def sphere_drag_quadrature(stress_fn: Callable, a: float, n: int = 64) -> float:
    """Drag (x-force) on an axisymmetric sphere from its surface x-traction by Gauss–Legendre quadrature in cos θ.

    Book: §8.6 (Exercise 8.35: integrate the surface pressure and shear stress). D = 2πa² ∫₀^π t_x(θ) sin θ dθ =
    2πa² ∫₋₁¹ t_x(arccos μ) dμ (ours; spectrally convergent for smooth tractions).

    Parameters
    ----------
    stress_fn : callable θ ↦ t_x [Pa] (x-component of the traction exerted by the fluid on the sphere).
    a : radius [m];  n : number of Gauss–Legendre nodes.

    Returns
    -------
    D : [N].

    Validation — tests/test_ch08.py: test_stokes_drag_V1_quadrature_parts_and_slip_sphere,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V1 = 6πμaU for the Stokes tractions. The Stokes x-traction is uniform over the sphere, so Gauss–Legendre
    quadrature is exact from n = 1 (no convergence to observe). Label: analytic.
    """
    xg, wg = np.polynomial.legendre.leggauss(int(n))
    return float(2.0 * np.pi * float(a) ** 2 * np.sum(wg * _F(stress_fn(np.arccos(xg)))))


# ======================================================================================================================
# drag coefficients (8.52) and Oseen's C_D
# ======================================================================================================================
def stokes_drag_coefficient(Re):
    """Stokes drag coefficient C_D = D/(½ρU²πa²) = 24/Re, Re = 2aU/ν (diameter).

    Book: §8.6, Eq. (8.52). Parameters: Re [–]. Returns C_D [–].
    Validation — tests/test_ch08.py: test_stokes_law_V1_form_and_drag_coefficient,
      test_drag_laws_V5_oseen_proudman_pearson_morrison,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable,
      test_book_V6_section_8_4_to_8_6_numbers_and_slips.
    Checks: V1 from (8.51); parity
    ``core.similarity.sphere_drag_coefficient(Re, "stokes")``. Label: analytic.
    """
    return _S(24.0 / _F(Re))  # Eq. (8.52)


def oseen_drag_coefficient(Re):
    """Oseen's drag coefficient C_D = (24/Re)(1 + 3Re/16), Re = 2aU/ν (diameter).

    Book: §8.6 (after (8.53); stated). Parameters: Re [–]. Returns C_D [–].
    Validation — tests/test_ch08.py: test_drag_laws_V5_oseen_proudman_pearson_morrison,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable,
      test_book_V6_section_8_4_to_8_6_numbers_and_slips.
    Checks: V5 Wikipedia
    "Oseen equations" radius form (12/Re_a)(1 + 3Re_a/8) with Re_a = Re/2 — an algebraic cross-check of the formula's
    form (no data compared). Label: analytic (form cross-check).
    """
    R = _F(Re)
    return _S(24.0 / R * (1.0 + 3.0 * R / 16.0))


def proudman_pearson_drag_coefficient(Re):
    """Proudman–Pearson (1957) matched-expansion drag, C_D = (24/Re)[1 + (3/8)Re_a + (9/40)Re_a² ln Re_a], Re_a = Re/2.

    Book: §8.6 (named: Kaplun, Proudman & Pearson; the formula is not printed — Wikipedia "Oseen equations").
    Parameters: Re = 2aU/ν [–] (converted to the radius Reynolds number). Returns C_D [–].
    Validation — tests/test_ch08.py: test_drag_laws_V5_oseen_proudman_pearson_morrison (its form against the cited
    radius-Re expression and its Re → 0 approach to Stokes' 24/Re; no data compared). Label: analytic (form cross-check).
    """
    R = _F(Re)
    Ra = R / 2.0
    return _S(24.0 / R * (1.0 + 0.375 * Ra + 0.225 * Ra ** 2 * np.log(Ra)))


# ======================================================================================================================
# settling, Millikan
# ======================================================================================================================
def terminal_velocity(a, rho_p, rho, mu, g: float = G0, warn: bool = True):
    """Terminal (settling) velocity of a small sphere: (4/3)πa³g(ρ′ − ρ) = 6πμaU ⇒ U = 2(ρ′ − ρ)ga²/(9μ).

    Book: §8.6 (after (8.51), Millikan's use of Stokes' law). Positive = downward (ρ′ > ρ).

    Parameters
    ----------
    a : radius [m];  rho_p : particle density ρ′ [kg/m³];  rho : fluid density [kg/m³];  mu : [Pa s];  g : [m/s²];
    warn : warn if Re = 2a|U|ρ/μ > 0.1 (Stokes' law then overestimates U).

    Returns
    -------
    U : [m/s].

    Validation — tests/test_ch08.py: test_stokes_law_V1_form_and_drag_coefficient,
      test_terminal_velocity_V1_force_balance_round_trip_and_warning, test_millikan_V5_synthetic_experiment_recovers_e,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V1 drag = effective weight. Label: analytic.
    """
    a_ = _F(a)
    U = 2.0 * (float(rho_p) - float(rho)) * float(g) * a_ ** 2 / (9.0 * float(mu))
    if warn and np.any(2.0 * a_ * np.abs(U) * float(rho) / float(mu) > 0.1):
        warnings.warn("terminal_velocity: Re = 2aU/ν > 0.1 — outside the creeping-flow range", RuntimeWarning,
                      stacklevel=2)
    return _S(U)


def radius_from_terminal_velocity(U, rho_p, rho, mu, g: float = G0):
    """Radius of a sphere from its Stokes terminal velocity, a = √(9μU/(2(ρ′ − ρ)g)) (Millikan's first step).

    Book: §8.6. Parameters: U [m/s] (> 0 downward); rho_p, rho [kg/m³]; mu [Pa s]; g [m/s²]. Returns a [m].
    Validation — tests/test_ch08.py: test_terminal_velocity_V1_force_balance_round_trip_and_warning,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V1 round trip with :func:`terminal_velocity`. Label: analytic.
    """
    return _S(np.sqrt(9.0 * float(mu) * _F(U) / (2.0 * (float(rho_p) - float(rho)) * float(g))))


def millikan_charge(U_fall, U_rise, rho_p, rho, mu, E, g: float = G0):
    """Charge q = ne on an oil drop from its fall speed (field off) and rise speed U_u (field on).

    Book: §8.6, Millikan (1911): 6πμU_u a + (4/3)πa³g(ρ′ − ρ) = neE (upward electric force balances drag + effective
    weight), a from the fall speed. Since (4/3)πa³g(ρ′ − ρ) = 6πμaU_fall this is q = 6πμa(U_u + U_fall)/E.
    E is the field magnitude [V/m] (the book's E = −V_b/L with the plates switched so the force is upward).

    Parameters
    ----------
    U_fall, U_rise : terminal speeds, both > 0 [m/s];  rho_p, rho : [kg/m³];  mu : [Pa s];  E : [V/m];  g : [m/s²].

    Returns
    -------
    q : [C].

    Validation — tests/test_ch08.py: test_millikan_V5_synthetic_experiment_recovers_e,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V1 identity; V5 synthetic drops recover e (CODATA). Label: analytic.
    """
    a = _F(radius_from_terminal_velocity(U_fall, rho_p, rho, mu, g))
    q = (6.0 * np.pi * float(mu) * _F(U_rise) * a + 4.0 / 3.0 * np.pi * a ** 3 * float(g) * (float(rho_p) - float(rho))) / float(E)
    return _S(q)


def settling_state(a, rho_p, rho, mu, g: float = G0) -> dict:
    """Numbers the settling-sphere explainer shows (curation §9, E9).

    Book: §8.6, (8.51)–(8.52) and the terminal-velocity balance.
    Parameters: a [m]; rho_p, rho [kg/m³]; mu [Pa s]; g [m/s²].
    Returns dict of floats: U_t [m/s], Re = 2aU_tρ/μ [–], D [N] (= effective weight), C_D = 24/Re [–], valid
    (bool, Re < 0.1), D_pressure, D_friction [N], weight_eff [N].
    The ``valid`` cut-off Re < 0.1 is **ours** (a conservative creeping-flow flag, the same as the warning in
    :func:`terminal_velocity`); the book only says the Stokes/Oseen results are "fairly accurate for Re < 5" relative
    to experiment, and Stokes' C_D = 24/Re is already ≈ 2 % low at Re = 0.1 (Oseen's 1 + 3Re/16).
    Validation — tests/test_ch08.py: test_low_re_scaling_V2_engine_and_D26_derivation,
      test_stokes_pressure_V1_extremes_and_printed_minimum,
      test_terminal_velocity_V1_force_balance_round_trip_and_warning,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V1 against the component functions. Label: analytic.
    """
    a, mu, rho = float(a), float(mu), float(rho)
    U = float(terminal_velocity(a, rho_p, rho, mu, g, warn=False))
    Re = 2.0 * a * abs(U) * rho / mu
    d = stokes_drag(mu, a, U, parts=True)
    return dict(U_t=U, Re=Re, D=d["total"], C_D=float(stokes_drag_coefficient(Re)) if Re > 0 else np.inf,
                valid=bool(Re < 0.1), D_pressure=d["pressure"], D_friction=d["friction"],
                weight_eff=4.0 / 3.0 * np.pi * a ** 3 * float(g) * (float(rho_p) - rho))


def inertia_viscous_ratio(r, theta, U: float = 1.0, a: float = 1.0, nu: float = 1e-6, h_rel: float = 1e-4):
    """Ratio |u·∇u| / |ν∇²u| evaluated on Stokes' solution — the far-field breakdown of creeping flow.

    Book: §8.6 (after Fig. 8.19): viscous force/volume ~ μUa/r³, inertia ~ ρU²a/r² ⇒ ratio ~ (ρUa/μ)(r/a) = Re r/a
    as r → ∞: inertia matters beyond r/a ~ 1/Re. Evaluated with second-order central differences of the Cartesian
    field (8.49) in the plane z = 0 (point (r cos θ, r sin θ, 0)), step h = h_rel·r.

    Which Re: the book's estimate uses the radius Reynolds number Re_a = ρUa/μ = Ua/ν; this module's convention
    elsewhere (drag coefficients, :func:`oseen_streamfunction`) is the diameter Re = 2aU/ν = 2 Re_a.

    The O(1) prefactor (ours, sympy on (8.48)–(8.49), body frame): at leading order in a/r the θ-component of u·∇u has
    no O(U²a/r²) term, and the radial one is (3U²a/16r²)(8cos²θ − 4sin²θ) = (3U²a/4r²)(2 − 3sin²θ); the viscous term
    ν∇²u = ∇p/ρ has magnitude (3νUa/2r³)(4cos²θ + sin²θ)^{1/2}. Hence on the axis (θ = 0, π) and at θ = π/2
    ratio → (1/2) Re_a (r/a) = (1/4) Re (r/a), so the crossover (ratio = 1) is at r/a ≈ 2/Re_a = 4/Re there; near
    sin²θ = 2/3 (θ ≈ 54.7°) the leading inertia term vanishes and the ratio is smaller (≈ 0.16 Re_a r/a at θ = π/4).
    The function returns the full finite-difference ratio, not this asymptote.

    Parameters
    ----------
    r : [m] (> a);  theta : [rad] from +x;  U : [m/s];  a : [m];  nu : ν [m²/s];  h_rel : relative step.

    Returns
    -------
    ratio : [–].

    Validation — tests/test_ch08.py: test_inertia_viscous_ratio_V7_linear_growth_in_r_and_Re,
      test_inertia_viscous_ratio_V1_half_Re_a_r_asymptote,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V7 grows linearly in r with slope ∝ Ua/ν (log–log slope 1 ± 0.05 for r/a ∈ [50, 500]).
    Label: analytic.
    """
    r_, t_ = np.broadcast_arrays(_F(r), _F(theta))
    shape = r_.shape
    X = np.stack([(r_ * np.cos(t_)).ravel(), (r_ * np.sin(t_)).ravel(), np.zeros(r_.size)])
    Uf = lambda P, T: np.stack(stokes_sphere_velocity_xyz(P[0], P[1], P[2], U, a, "body"))  # noqa: E731
    hstep = float(h_rel) * float(np.min(r_))
    u = Uf(X, 0.0)
    Gm = st.grad_vector(Uf, X, 0.0, hstep)  # G[i, j] = ∂u_i/∂x_j
    adv = np.einsum("jn,ijn->in", u, Gm)
    lap = float(nu) * st.laplacian(Uf, X, 0.0, hstep, (3,))
    ratio = np.linalg.norm(adv, axis=0) / np.linalg.norm(lap, axis=0)
    return _S(ratio.reshape(shape))


# ======================================================================================================================
# Oseen (8.53)
# ======================================================================================================================
def _phi1(s):
    """(1 − e^{−s})/s, stable for s → 0 (= 1 at s = 0)."""
    s = _F(s)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(np.abs(s) > 1e-300, -np.expm1(-s) / np.where(np.abs(s) > 1e-300, s, 1.0), 1.0)


def oseen_streamfunction(r, theta, U: float = 1.0, a: float = 1.0, Re: float = 0.1, frame: str = "body"):
    """Oseen's stream function for low-Re flow past a sphere.

    Book: §8.6, Eq. (8.53): ψ/(Ua²) = [r²/2a² + a/4r] sin²θ − (3/Re)(1 + cos θ){1 − exp[−(Re/4)(r/a)(1 − cos θ)]},
    Re = 2aU/ν. Evaluated in the algebraically identical form −(3/4)(r/a) sin²θ · (1 − e^{−s})/s, s = (Re/4)(r/a)(1 −
    cos θ), which is stable as Re → 0 (and at Re = 0 gives Stokes' (8.48) exactly). It satisfies no slip only to O(Re)
    (analysis §9 R18). ``frame="fluid"`` subtracts ½Ur² sin²θ (Fig. 8.20).

    Parameters
    ----------
    r : [m];  theta : [rad] from +x (downstream);  U : [m/s];  a : [m];  Re : 2aU/ν [–] (≥ 0);  frame.

    Returns
    -------
    psi : [m³/s] (NaN inside).

    Validation — tests/test_ch08.py: test_oseen_V1_velocity_axis_wake_and_no_slip_order,
      test_oseen_V2_streamfunction_solves_oseen_equation_and_wake_downstream,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable,
      test_book_V6_section_8_4_to_8_6_numbers_and_slips.
    Checks: V2 Re → 0 limit = (8.48); V1 ψ = 0 on the axis; wake asymmetry. Label: symbolic, analytic.
    """
    r_, t_ = np.broadcast_arrays(_F(r), _F(theta))
    a, U, Re = float(a), float(U), float(Re)
    c, s2 = np.cos(t_), np.sin(t_) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        q = r_ / a
        s = Re / 4.0 * q * (1.0 - c)
        psi = U * a ** 2 * ((q ** 2 / 2.0 + 1.0 / (4.0 * q)) * s2 - 0.75 * q * s2 * _phi1(s))  # Eq. (8.53)
    if _frame(frame):
        psi = psi - 0.5 * U * r_ ** 2 * s2
    return _S(np.where(_outside(r_, a), psi, np.nan))


def oseen_velocity(r, theta, U: float = 1.0, a: float = 1.0, Re: float = 0.1, frame: str = "body"):
    """Velocity components of Oseen's solution (8.53), by analytic differentiation through (6.83) (ours).

    Book: §8.6, Eq. (8.53) with u_r = ψ_θ/(r² sin θ), u_θ = −ψ_r/(r sin θ):
    u_r = U[cos θ(1 + a³/2r³) + (3a/4r)((1 − cos θ)φ₁(s) − (1 + cos θ)e^{−s})],
    u_θ = −U sin θ[1 − a³/4r³ − (3a/4r)e^{−s}], s = (Re/4)(r/a)(1 − cos θ), φ₁(s) = (1 − e^{−s})/s.
    Re → 0 gives (8.49). ``frame="fluid"`` subtracts the stream U e_x.

    Parameters
    ----------
    r : [m];  theta : [rad];  U : [m/s];  a : [m];  Re : [–];  frame : "body" | "fluid".

    Returns
    -------
    (u_r, u_theta) : [m/s] (NaN inside).

    Validation — tests/test_ch08.py: test_oseen_V1_velocity_axis_wake_and_no_slip_order,
      test_part_c_V1_every_contract_function_exists_and_is_scalar_callable.
    Checks: V1 finite differences of :func:`oseen_streamfunction`; Re → 0 equals (8.49). Label: analytic.
    """
    r_, t_ = np.broadcast_arrays(_F(r), _F(theta))
    a, U, Re = float(a), float(U), float(Re)
    c, sn = np.cos(t_), np.sin(t_)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        q = a / r_
        s = Re / 4.0 / q * (1.0 - c)
        es = np.exp(-s)
        ur = U * (c * (1.0 + 0.5 * q ** 3) + 0.75 * q * ((1.0 - c) * _phi1(s) - (1.0 + c) * es))
        ut = -U * sn * (1.0 - 0.25 * q ** 3 - 0.75 * q * es)
    if _frame(frame):
        ur, ut = ur - U * c, ut + U * sn
    m = _outside(r_, a)
    return _S(np.where(m, ur, np.nan)), _S(np.where(m, ut, np.nan))
