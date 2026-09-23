"""Kinematics of a fluid: particle and field descriptions, the material derivative, flow lines, frame changes and the
motion of a small fluid element (strain, rotation, volume change, principal axes).

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e, Ch. 3 §§3.2–3.5, Eqs. (3.1)–(3.21) (read from the rendered pages
chapters/pages/ch03/p096–p109). The linear-flow helpers ``velocity_gradient_preset``, ``linear_flow_map``,
``deform_square`` and ``material_line_angle`` were introduced in ch02 (§2.10–2.11) and moved here when Ch. 3 became
their second user; ``fluidpy.ch02_cartesian_tensors`` re-exports them unchanged.

Conventions (project-wide, ``knowledge/notation.md``)
-----------------------------------------------------
* A **velocity field** is a callable ``u(x, t)`` with ``x`` of shape ``(d,)`` (one point) or ``(d, N)`` (N points,
  coordinates on axis 0), ``t`` a float [s]; it returns an array of the same shape as ``x`` [m/s]. A **scalar field**
  ``F(x, t)`` returns a float or an ``(N,)`` array (any unit). ``d`` is 2 (plane flow) or 3.
  (The ch02 integral-theorem routines take ``fn(X, Y[, Z])``; :func:`as_coord_field` adapts.)
* Velocity gradient ``G[i, j] = ∂u_i/∂x_j`` [1/s]; strain rate S = ½(G + Gᵀ) (3.12); the book's rotation tensor
  R = G − Gᵀ (3.13, **no ½**) whose vector is the vorticity ω = ∇×u (3.15); the element spins at ½ω.
* Angles in radians, counter-clockwise positive (so a clockwise spin is negative: ω₃ = −γ in the shear u = (γx₂, 0)).
* Every numerical choice the book leaves open (finite-difference step ``h``, ODE tolerances, resolutions) is an
  argument with a documented default, so convergence studies can vary it.
"""
from __future__ import annotations

from typing import Callable, NamedTuple, Sequence

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.linalg import expm

from ._util import as_scalar_if_0d
from .tensors import levi_civita, principal_axes, rotation_tensor, strain_rate_tensor

__all__ = [
    # promoted from ch02 (linear flows)
    "VELOCITY_GRADIENT_PRESETS", "velocity_gradient_preset", "linear_flow_map", "deform_square", "material_line_angle",
    # field adapters
    "as_coord_field", "from_coord_field",
    # §3.2 Lagrangian / Eulerian, material derivative
    "lagrangian_velocity_acceleration", "lagrangian_velocity_acceleration_sym", "lagrangian_to_eulerian",
    "MaterialDerivativeTerms", "material_derivative", "material_derivative_terms", "material_derivative_sym",
    "streamwise_derivative",
    # §3.3 flow lines, Galilean frames
    "streamline", "pathline", "streakline", "galilean_transform", "AccelerationTerms", "acceleration",
    # §3.4 element kinematics
    "velocity_gradient_at", "relative_velocity", "vorticity_from_gradient", "vorticity", "linear_strain_rate",
    "shear_strain_rate", "volumetric_strain_rate", "material_volume_ratio", "element_rotation_rate",
    "material_line_rotation_rate", "perpendicular_pair_rotation_rate", "vorticity_in_rotating_frame",
    "velocity_potential_2d", "RelativeVelocitySplit", "relative_velocity_split", "principal_strain_rates",
    "strain_velocity_principal", "deform_circle", "deform_sphere", "strain_ellipse_axes", "measured_strain_rates",
]

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731


# ======================================================================================================================
# Linear flows u = G·x (promoted from ch02 §2.10–2.11; used by ch03 §3.4–3.5)
# ======================================================================================================================
VELOCITY_GRADIENT_PRESETS = ("simple_shear", "solid_body_rotation", "pure_strain", "uniaxial_extension", "irrotational_strain")


def velocity_gradient_preset(name: str, Gamma: float = 1.0, dim: int = 2) -> np.ndarray:
    """Velocity-gradient matrices G[i, j] = ∂u_i/∂x_j of the standard linear flows (E3 presets), rate scale Γ [1/s].

    Book: §2.10 (S + A split of ∂u_i/∂x_j), Examples 2.3 and 2.4; the names and matrices are the design's Part C 6.6
    / E3 preset list. Ch. 3 §3.5 uses the ``simple_shear`` preset with Γ = γ = du₁/dx₂.

    * ``simple_shear``: u = (Γ x₂, 0) → G = [[0, Γ], [0, 0]] (Example 2.4's flow; S₁₂ = Γ/2 and A₁₂ = Γ/2 —
      half strain, half rotation; ∇·u = 0)
    * ``solid_body_rotation``: u = Γ e₃ × x = (−Γ x₂, Γ x₁) → G = [[0, −Γ], [Γ, 0]] (Example 2.3 with b = Γ e₃;
      S = 0, ∇·u = 0, (∇×u)₃ = 2Γ)
    * ``pure_strain``: u = (Γ x₁, −Γ x₂) → G = diag(Γ, −Γ) (A = 0, traceless: stretch along x₁, squeeze along x₂,
      area preserved; Example 2.4's S' in its principal frame)
    * ``uniaxial_extension``: u = (Γ x₁, 0) → G = diag(Γ, 0) (A = 0, **∇·u = Γ ≠ 0**: stretch along x₁ only, the
      area grows as e^{Γt})
    * ``irrotational_strain``: u = (Γ x₂, Γ x₁) → G = [[0, Γ], [Γ, 0]] (A = 0, traceless; Example 2.4's S itself as
      a flow — the simple shear with its rotation removed, principal axes at ±45°)

    In 3-D the matrices are embedded in the (1, 2) block with ∂u₃/∂x₃ = 0.

    Validation: V1 S/A split of each preset (A = 0 for the three strain presets, S = 0 for the rotation);
    trace(G) = ∇·u (Γ for uniaxial extension, 0 otherwise); ``linear_flow_map`` limits. Label: analytic.
    """
    if name not in VELOCITY_GRADIENT_PRESETS:
        raise ValueError(f"unknown preset {name!r}; choose from {VELOCITY_GRADIENT_PRESETS}")
    Gam = float(Gamma)
    G2 = {"simple_shear": [[0.0, Gam], [0.0, 0.0]],  # u₁ = Γ x₂
          "solid_body_rotation": [[0.0, -Gam], [Gam, 0.0]],  # u = Γ e₃ × x
          "pure_strain": [[Gam, 0.0], [0.0, -Gam]],  # u = (Γ x₁, −Γ x₂)  (traceless)
          "uniaxial_extension": [[Gam, 0.0], [0.0, 0.0]],  # u = (Γ x₁, 0)  (∇·u = Γ)
          "irrotational_strain": [[0.0, Gam], [Gam, 0.0]]}[name]  # u = (Γ x₂, Γ x₁)  (simple shear minus its rotation)
    G = np.array(G2)
    if dim == 3:
        G3 = np.zeros((3, 3))
        G3[:2, :2] = G
        return G3
    return G


def linear_flow_map(G, t) -> np.ndarray:
    """Position map of the linear flow u = G·x after time t: x(t) = e^{Gt} x₀ (``scipy.linalg.expm``).

    Book: §2.10 (S + A decomposition seen as motion; E3), Ch. 3 §3.4 (tracked material segments, Figs. 3.10–3.13).
    Units: G [1/s], t [s]. Assumes G constant in space and time (a linear, steady flow — or the local linearisation
    (3.10) of any smooth flow for short times).

    Validation: V1 pure rotation G = A gives a rotation matrix by angle ω₃ t; uniaxial strain gives diag(e^{Γt}, e^{−Γt});
    d/dt of the map at t = 0 equals G. Label: analytic.
    """
    return expm(_F(G) * float(t))  # x(t) = e^{Gt} x₀ solves dx/dt = G x


def deform_square(G, t, n_side: int = 10, half_width: float = 1.0, boundary_only: bool = False) -> np.ndarray:
    """Tracer positions of a square of material points carried by the linear flow u = G·x for time t.

    Returns an array ``(2, N)`` (or ``(3, N)`` for a 3 × 3 G — a square in the (1, 2) plane) of positions at time t;
    ``boundary_only`` samples the perimeter (4 n_side points, closed) instead of the filled n_side × n_side lattice.

    Book: E3 "strain vs rotation split", Fig. 2.8's deforming square; Ch. 3 Fig. 3.14 (elements ABCD and PQRS in a
    parallel shear flow). Units: G [1/s], t [s], half_width [m]. Label: analytic.
    """
    G_ = _F(G)
    d = G_.shape[0]
    s = np.linspace(-half_width, half_width, n_side)
    if boundary_only:
        pts = np.concatenate([np.stack([s, np.full_like(s, -half_width)]), np.stack([np.full_like(s, half_width), s]),
                              np.stack([s[::-1], np.full_like(s, half_width)]), np.stack([np.full_like(s, -half_width), s[::-1]])],
                             axis=1)
    else:
        X, Y = np.meshgrid(s, s, indexing="xy")
        pts = np.stack([X.ravel(), Y.ravel()])
    if d == 3:
        pts = np.vstack([pts, np.zeros((1, pts.shape[1]))])
    return linear_flow_map(G_, t) @ pts


def material_line_angle(G, t, theta0=0.0):
    """Angle of a material line element that starts at angle θ₀ after time t in the linear flow u = G·x (2-D).

    Book: §2.10 (S + A split), Ch. 3 §3.4 (Fig. 3.11: the angles α, β of two material lines). Under solid-body
    rotation the line turns uniformly at the angular velocity of the fluid element, which is the vector of the
    antisymmetric part A: Ω = ½(∇×u)₃ = ½ω₃ (ω the book's vorticity, = the vector of ``rotation_tensor``); for the
    ``solid_body_rotation`` preset with G = [[0, −Γ], [Γ, 0]] that is Ω = Γ, so θ(t) = θ₀ + Γt. Under pure strain the
    line tends to the stretching axis; under simple shear it does both (E3 view 3). Units: rad (θ₀ may be an array).

    Validation: V1 solid-body rotation: θ(t) = θ₀ + Γ t (= θ₀ + ½(∇×u)₃ t); pure strain: θ → 0 as t → ∞.
    Label: analytic.
    """
    M = linear_flow_map(G, t)[:2, :2]
    th0 = _F(theta0)
    v = M @ np.stack([np.cos(th0), np.sin(th0)])
    return as_scalar_if_0d(np.arctan2(v[1], v[0]))


# ======================================================================================================================
# Field adapters
# ======================================================================================================================
def as_coord_field(u: Callable, t: float = 0.0) -> Callable:
    """Adapt a kinematics field ``u(x, t)`` to the ch02 convention ``fn(X, Y[, Z])`` (components on axis 0) at time t,
    so that ``core.integral_theorems.circulation``, ``flux_through_sphere`` … can integrate it.

    Label: analytic (pure re-packaging).
    Book: §3.3–3.5 (adapter only; no equation).
    """
    def fn(*coords):
        B = np.broadcast_arrays(*[_F(c) for c in coords])
        shape = B[0].shape
        X = np.stack([b.ravel() for b in B])
        U = _F(u(X, t))
        U = np.broadcast_to(U, X.shape) if U.shape != X.shape else U
        return U.reshape((U.shape[0],) + shape)
    return fn


def from_coord_field(fn: Callable) -> Callable:
    """Adapt a ch02-style callable ``fn(X, Y[, Z])`` (e.g. a :class:`~fluidpy.core.fields.VectorField`) to the
    kinematics convention ``u(x, t)`` (steady: t is ignored). Label: analytic.
    Book: §3.3–3.5 (adapter only; no equation).
    """
    def u(x, t=0.0):
        x_ = _F(x)
        out = _F(fn(*x_))
        return np.broadcast_to(out, x_.shape).copy() if out.shape != x_.shape else out
    return u


def _u_at(u: Callable, x, t) -> np.ndarray:
    """Evaluate a velocity field and make sure the result has the shape of ``x``."""
    x_ = _F(x)
    U = _F(u(x_, t))
    if U.shape != x_.shape:
        U = np.broadcast_to(U.reshape(U.shape + (1,) * (x_.ndim - U.ndim)), x_.shape).copy()
    return U


def _unit(d: int, i: int, ndim: int) -> np.ndarray:
    e = np.zeros(d)
    e[i] = 1.0
    return e.reshape((d,) + (1,) * (ndim - 1))


# ======================================================================================================================
# §3.2 Lagrangian and Eulerian descriptions; the material derivative
# ======================================================================================================================
def lagrangian_velocity_acceleration(r_of_t: Callable, t: float, h: float = 1e-4):
    """Velocity and acceleration of a fluid particle from its trajectory r(t; r_o, t_o), Eq. (3.1).

    Book: §3.2, Eq. (3.1): u = dr(t; r_o, t_o)/dt and a = d²r(t; r_o, t_o)/dt². The labels r_o, t_o are fixed
    parameters of the callable, not variables.

    Parameters
    ----------
    r_of_t : callable t → position [m] (float for a 1-D map, array (d,) otherwise)
    t : time [s]
    h : finite-difference step [s] (explicit 2nd-order central stencils; error O(h²) truncation + O(ε/h²) rounding
        in a — h = 1e-4 gives ~1e-8 relative accuracy for smooth maps)

    Returns
    -------
    (u, a) : velocity [m/s] and acceleration [m/s²] (floats for a scalar map)

    Validation: V1 x = X e^{αt}: u = αx, a = α²x; V3 order 2 in h. Label: analytic, converged.
    """
    rp, r0, rm = _F(r_of_t(t + h)), _F(r_of_t(t)), _F(r_of_t(t - h))
    u = (rp - rm) / (2.0 * h)  # Eq. (3.1): u = dr/dt
    a = (rp - 2.0 * r0 + rm) / h ** 2  # Eq. (3.1): a = d²r/dt²
    return as_scalar_if_0d(u), as_scalar_if_0d(a)


def lagrangian_velocity_acceleration_sym(r_exprs, t: sp.Symbol):
    """Symbolic Eq. (3.1): u_i = dr_i/dt, a_i = d²r_i/dt² for a Lagrangian map given as sympy expressions in t and
    the labels. Returns (u_exprs, a_exprs) as lists. Book: §3.2, Eq. (3.1). Label: symbolic."""
    r = [sp.sympify(e) for e in (r_exprs if isinstance(r_exprs, (list, tuple)) else [r_exprs])]
    return [sp.diff(e, t) for e in r], [sp.diff(e, t, 2) for e in r]  # Eq. (3.1)


def lagrangian_to_eulerian(r_exprs, labels, coords, t: sp.Symbol) -> dict:
    """Turn a Lagrangian map x = r(t; r_o) into the Eulerian velocity u(x, t) by the compatibility condition (3.2).

    Book: §3.2, Eqs. (3.1)–(3.2): F[r(t; r_o, t_o), t] = F(x, t) when x = r(t; r_o, t_o). D01 (our worked case): with
    x = X e^{αt}, u = dr/dt = αX e^{αt}; solving x = X e^{αt} for the label X = x e^{−αt} and substituting gives the
    Eulerian u(x, t) = αx.

    Parameters
    ----------
    r_exprs : list of sympy expressions r_i(t, labels) [m]
    labels : list of label symbols (r_o components)
    coords : list of Eulerian coordinate symbols x_i (same length)
    t : time symbol

    Returns
    -------
    dict(label_of_x={label: expr in x, t}, u=[u_i(x, t)], a_lagrangian=[d²r_i/dt² written in x, t],
    Du_Dt=[∂u_i/∂t + u_j ∂u_i/∂x_j from the Eulerian u]) — ``simplify(Du_Dt − a_lagrangian) == 0`` is the D01 check;
    ``a`` is an alias of ``a_lagrangian``.

    Assumptions: the map is invertible (particles never merge), smooth.
    Validation: V2 Du/Dt computed from the Eulerian u equals d²r/dt² from the map (three maps: stretching, rotation,
    shear). Label: symbolic.
    """
    r = [sp.sympify(e) for e in r_exprs]
    eqs = [sp.Eq(x, ri) for x, ri in zip(coords, r)]
    sol = sp.solve(eqs, list(labels), dict=True)
    if not sol:
        raise ValueError("could not invert the Lagrangian map for the labels (not invertible?)")
    sub = sol[0]
    u = [sp.simplify(sp.diff(ri, t).subs(sub)) for ri in r]  # (3.1) then (3.2): u(x, t) = [dr/dt] at r_o = r_o(x, t)
    a = [sp.simplify(sp.diff(ri, t, 2).subs(sub)) for ri in r]
    DuDt = [sp.simplify(e) for e in material_derivative_sym(u, u, list(coords), t)]  # Eq. (3.5) on the Eulerian twin
    return {"label_of_x": sub, "u": u, "a_lagrangian": a, "a": a, "Du_Dt": DuDt}


class MaterialDerivativeTerms(NamedTuple):
    """The split of Eq. (3.5): ``local`` ∂F/∂t, ``advective`` u·∇F, ``total`` DF/Dt (units of F per second).
    Book: §3.2, Eq. (3.5). Label: analytic.
    """

    local: float | np.ndarray
    advective: float | np.ndarray
    total: float | np.ndarray


def material_derivative_terms(F: Callable, u: Callable, x, t: float, h: float = 1e-4,
                              ht: float | None = None) -> MaterialDerivativeTerms:
    """Local, advective and total rate of change of a field following the fluid, Eqs. (3.4)–(3.5).

    Book: §3.2, Eq. (3.4): d/dt F[r(t; r_o, t_o), t] = (∇F)·u + ∂F/∂t ≡ DF/Dt, and Eq. (3.5):
    DF/Dt = ∂F/∂t + u·∇F = ∂F/∂t + u_i ∂F/∂x_i. The *unsteady* (local) part vanishes for steady F; the *advective*
    part vanishes when F is uniform, u = 0 or u ⟂ ∇F.

    Parameters
    ----------
    F : callable F(x, t) — scalar field (float or (N,) array) or a vector field (same layout as u) — any unit
    u : callable u(x, t) [m/s]
    x : point(s) [m], shape (d,) or (d, N)
    t : time [s]
    h : spatial step [m] of the explicit 2nd-order central stencils (never ``np.gradient``)
    ht : time step [s] (default h)

    Returns
    -------
    MaterialDerivativeTerms(local, advective, total) — [unit of F / s]; floats for a single point and scalar F.

    Assumptions: F, u differentiable at x; one fluid particle per point (continuum).
    Validation: V1/V3 DF/Dt equals d/dt F(pathline(t), t) (an ODE + chain-rule route independent of the stencils);
    V2 equals ``material_derivative_sym``; V3 stencil order 2; local = 0 for steady F, advective = 0 for u ⟂ ∇F.
    Label: analytic, symbolic, converged.
    """
    x_ = _F(x)
    ht = h if ht is None else ht
    d = x_.shape[0]
    local = (_F(F(x_, t + ht)) - _F(F(x_, t - ht))) / (2.0 * ht)  # ∂F/∂t at fixed x
    U = _u_at(u, x_, t)
    adv = 0.0
    for i in range(d):
        e = _unit(d, i, x_.ndim) * h
        dFdxi = (_F(F(x_ + e, t)) - _F(F(x_ - e, t))) / (2.0 * h)  # ∂F/∂x_i
        adv = adv + U[i] * dFdxi  # Eq. (3.5): u_i ∂F/∂x_i
    total = local + adv  # Eq. (3.5): DF/Dt = ∂F/∂t + u·∇F
    return MaterialDerivativeTerms(as_scalar_if_0d(local), as_scalar_if_0d(adv), as_scalar_if_0d(total))


def material_derivative(F: Callable, u: Callable, x, t: float, h: float = 1e-4, ht: float | None = None):
    """DF/Dt = ∂F/∂t + u·∇F at x, t — the total of :func:`material_derivative_terms` (Eqs. (3.4)–(3.5)).

    Book: §3.2, Eqs. (3.4), (3.5). Units: [unit of F / s]. Validation: as ``material_derivative_terms``.
    Label: analytic, symbolic, converged.
    """
    return material_derivative_terms(F, u, x, t, h, ht).total


def material_derivative_sym(F_expr, u_exprs, coords, t: sp.Symbol):
    """Symbolic DF/Dt = ∂F/∂t + u_i ∂F/∂x_i (Eq. (3.5)); ``F_expr`` may be a list (vector F, component-wise).

    Book: §3.2, Eqs. (3.3)–(3.5). Validation: V2 Du/Dt of the Eulerian twin of a Lagrangian map equals d²r/dt²
    (D01/D02). Label: symbolic.
    """
    def one(Fe):
        Fe = sp.sympify(Fe)
        return sp.diff(Fe, t) + sum(ui * sp.diff(Fe, xi) for ui, xi in zip(u_exprs, coords))  # Eq. (3.5)
    if isinstance(F_expr, (list, tuple)):
        return [one(Fe) for Fe in F_expr]
    return one(F_expr)


def streamwise_derivative(F: Callable, u: Callable, x, t: float, h: float = 1e-4):
    """The advective term written along the trajectory: |u| ∂F/∂s with ∂F/∂s = e_u·∇F, Eq. (3.6).

    Book: §3.2, Eq. (3.6), DF/Dt = ∂F/∂t + |u| ∂F/∂s, where s is arc length along the particle path, dr = e_u ds,
    e_u = u/|u|. ⚠️ The book prints the second term as |u| ∂/∂s (the F is missing); both terms must be rates of
    change of F, so we implement |u| ∂F/∂s (analysis §9 item 1).

    Parameters
    ----------
    F, u, x, t, h : as :func:`material_derivative_terms` (x a single point (d,) or points (d, N))

    Returns
    -------
    |u| ∂F/∂s [unit of F / s] — equals the advective term u·∇F.

    Raises
    ------
    ValueError at a stagnation point (|u| = 0: the direction e_u is undefined).

    Validation: V1 equals the advective term of ``material_derivative_terms`` at random points; raises at u = 0.
    Label: analytic.
    """
    x_ = _F(x)
    U = _u_at(u, x_, t)
    speed = np.sqrt(np.sum(U ** 2, axis=0))
    if np.any(speed == 0.0):
        raise ValueError("|u| = 0: the streamwise direction e_u = u/|u| is undefined at a stagnation point")
    e_u = U / speed
    dFds = (_F(F(x_ + h * e_u, t)) - _F(F(x_ - h * e_u, t))) / (2.0 * h)  # ∂F/∂s = e_u·∇F (directional derivative)
    return as_scalar_if_0d(speed * dFds)  # Eq. (3.6) (with the F restored): |u| ∂F/∂s


# ======================================================================================================================
# §3.3 Flow lines
# ======================================================================================================================
def streamline(u: Callable, x0, t_frozen: float = 0.0, s_max: float = 1.0, both: bool = True, n: int = 400,
               rtol: float = 1e-10, atol: float = 1e-12, eps: float = 1e-12) -> np.ndarray:
    """A streamline through x0 at the frozen instant t_frozen: the curve tangent to u everywhere, Eq. (3.7).

    Book: §3.3, Fig. 3.5 and Eq. (3.7): dx/u = dy/v = dz/w (equivalently u × ds = 0). We integrate the equivalent
    arc-length ODE dx/ds = u(x, t_frozen)/|u| (D03: ds ∥ u with |ds| = ds) with ``solve_ivp`` — the clock is frozen,
    so in unsteady flow the pattern belongs to one instant only.

    Parameters
    ----------
    u : callable u(x, t) [m/s]
    x0 : seed point [m], shape (d,)
    t_frozen : the instant t′ [s]
    s_max : arc length [m] integrated in each direction
    both : integrate upstream and downstream (samples s ∈ [−s_max, s_max]) or downstream only (s ∈ [0, s_max])
    n : number of samples (uniform in s)
    rtol, atol : ``solve_ivp`` (DOP853) tolerances
    eps : speed [m/s] below which the curve is stopped (terminal event: a stagnation point is reached)

    Returns
    -------
    ndarray (d, m) — points ordered by s (uniform spacing 2s_max/(n − 1)); m = n unless a stagnation point ended the
    curve, in which case the samples beyond it are dropped (m < n).

    Raises
    ------
    ValueError if x0 itself is a stagnation point.

    Validation: V1 Ex. 3.1 gives the line y = x tan ωt′; solid body gives a circle of constant radius; V4 ψ constant
    along the computed curve (cylinder, stagnation flow xy = c). Label: analytic, conserved.
    """
    x0 = _F(x0)
    if np.linalg.norm(_u_at(u, x0, t_frozen)) < eps:
        raise ValueError("x0 is a stagnation point (|u| < eps): no unique streamline direction")

    def rhs(sign):
        def f(s, x):
            U = _u_at(u, x, t_frozen)
            sp_ = np.linalg.norm(U)
            return sign * U / sp_ if sp_ > eps else np.zeros_like(U)  # dx/ds = ±u/|u|  (tangency, Eq. 3.7)
        return f

    def stop(s, x):
        return np.linalg.norm(_u_at(u, x, t_frozen)) - eps
    stop.terminal = True

    def run(sign, s_eval):
        out = np.full((x0.size, s_eval.size), np.nan)
        sol = solve_ivp(rhs(sign), (0.0, s_max), x0, method="DOP853", rtol=rtol, atol=atol, dense_output=True,
                        events=stop)
        s_end = sol.t[-1]
        ok = s_eval <= s_end + 1e-15
        if np.any(ok):
            out[:, ok] = sol.sol(s_eval[ok])
        return out

    if not both:
        s = np.linspace(0.0, s_max, n)
        pts = run(+1.0, s)
    else:
        s = np.linspace(-s_max, s_max, n)
        pts = np.empty((x0.size, n))
        fwd, bwd = s >= 0.0, s < 0.0
        pts[:, fwd] = run(+1.0, s[fwd])
        pts[:, bwd] = run(-1.0, -s[bwd])
    return pts[:, np.all(np.isfinite(pts), axis=0)]  # drop samples beyond a stagnation stop


def pathline(u: Callable, r0, t0: float, t_eval, rtol: float = 1e-10, atol: float = 1e-12,
             method: str = "DOP853") -> np.ndarray:
    """Path line of the particle that is at r0 at time t0: dr/dt = u(r, t), r(t0) = r0, Eq. (3.8).

    Book: §3.3, Eq. (3.8): dr/dt = [u(x, t)]_{x=r} = u(r, t) with r(t_o) = r_o — the Lagrangian trajectory
    r(t; r_o, t_o) of §3.2 obtained from the Eulerian field. Integration forward and/or backward in time.

    Parameters
    ----------
    u : callable u(x, t) [m/s]
    r0 : initial position [m], shape (d,)
    t0 : label time [s]
    t_eval : float or array of output times [s] (either side of t0)
    rtol, atol, method : ``solve_ivp`` settings (DOP853, 8th order)

    Returns
    -------
    ndarray (d, nt) (or (d,) for scalar t_eval) — positions r(t_eval) [m].

    Validation: V1 linear flows vs ``linear_flow_map`` (e^{G(t−t0)} r0); Ex. 3.1 circle; V3 error decreases
    monotonically with rtol. Label: analytic, converged.
    """
    r0 = _F(r0)
    scalar = np.ndim(t_eval) == 0
    te = np.atleast_1d(_F(t_eval))
    out = np.empty((r0.size, te.size))

    def rhs(t, r):
        return _u_at(u, r, t)  # Eq. (3.8)

    for mask, sign in ((te > t0, 1.0), (te < t0, -1.0)):
        idx = np.nonzero(mask)[0]
        if idx.size == 0:
            continue
        order = idx[np.argsort(sign * te[idx])]
        sol = solve_ivp(rhs, (t0, te[order[-1]]), r0, method=method, rtol=rtol, atol=atol, t_eval=te[order])
        if not sol.success:
            raise RuntimeError(f"pathline integration failed: {sol.message}")
        out[:, order] = sol.y
    out[:, te == t0] = r0[:, None]
    return out[:, 0] if scalar else out


def _u_many(u: Callable, X: np.ndarray, T: np.ndarray) -> np.ndarray:
    """u evaluated at N points X (d, N), each at its own time T (N,); vectorised if the field allows it."""
    try:
        U = _F(u(X, T))
        if U.shape == X.shape:
            return U
    except Exception:  # noqa: BLE001 — fields written for scalar t fall back to a loop
        pass
    return np.stack([_u_at(u, X[:, k], float(T[k])) for k in range(X.shape[1])], axis=1)


def streakline(u: Callable, x0, t: float, t_release, rtol: float = 1e-9, atol: float = 1e-11) -> np.ndarray:
    """Streak line through the fixed point x0 at time t: where every particle released at x0 at the times
    t_release (≤ t) is at time t.

    Book: §3.3 (streak line): integrate (3.8) for all release times t_o with r(t_o) = x_o; at the fixed time t the
    path-line positions x_i = r_i(t; x_o, t_o) give the streak line parametrically in t_o (dye from a port).

    Method (ours — the book prescribes none): each particle k is integrated on its own clock t_k(s) = t_o,k +
    s (t − t_o,k), s ∈ [0, 1], so dr_k/ds = (t − t_o,k) u(r_k, t_k(s)); all particles form one ODE system of size d·N
    (vectorised, one ``solve_ivp`` call).

    Parameters
    ----------
    u : callable u(x, t) [m/s] (if it accepts an array t of shape (N,) the ensemble is evaluated in one call)
    x0 : dye port [m], shape (d,)
    t : drawing time [s]
    t_release : release times t_o [s], each ≤ t
    rtol, atol : ``solve_ivp`` tolerances

    Returns
    -------
    ndarray (d, N) — positions at time t of the particles released at ``t_release`` (same order).

    Validation: V1 Ex. 3.1 circle (centre (ξ_o sin ωt′, −ξ_o cos ωt′), radius ξ_o); V7 steady flow: streak line =
    streamline = path line through x0. Label: analytic.
    """
    x0 = _F(x0)
    tau = np.atleast_1d(_F(t_release))
    if np.any(tau > t + 1e-12):
        raise ValueError("release times must satisfy t_o <= t (a streak line holds particles already released)")
    d, N = x0.size, tau.size
    dur = t - tau

    def rhs(s, y):
        X = y.reshape(d, N)
        return (_u_many(u, X, tau + s * dur) * dur).ravel()  # dr_k/ds = (t − t_o,k) u(r_k, t_o,k + s(t − t_o,k))

    sol = solve_ivp(rhs, (0.0, 1.0), np.repeat(x0[:, None], N, axis=1).ravel(), method="DOP853", rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(f"streakline integration failed: {sol.message}")
    return sol.y[:, -1].reshape(d, N)


# ======================================================================================================================
# §3.3 Galilean transformation and the fluid acceleration
# ======================================================================================================================
def galilean_transform(u: Callable, U, x0p=None) -> Callable:
    """The velocity field seen from a frame O′x′y′z′ that moves at constant velocity U with parallel axes.

    Book: §3.3, Fig. 3.8: u(x, t) = U + u′(x′, t′) with t = t′ and x = x′ + Ut + x′_o, so
    u′(x′, t′) = u(x′ + Ut′ + x′_o, t′) − U.

    Parameters
    ----------
    u : callable u(x, t) in the stationary frame Oxyz [m/s]
    U : constant frame velocity [m/s], shape (d,)
    x0p : x′_o, the vector from O to O′ at t = 0 [m] (default 0)

    Returns
    -------
    callable u_prime(xp, tp) [m/s] with the kinematics field convention.

    Assumptions: U constant (not a rotating frame: that adds Coriolis and centrifugal terms, Ch. 4 §4.7).
    Validation: V1 u′ + U = u at the same physical point; used by the (3.9) tests. Label: analytic.
    """
    U_ = _F(U)
    x0 = np.zeros_like(U_) if x0p is None else _F(x0p)

    def u_prime(xp, tp):
        xp_ = _F(xp)
        shape = (U_.size,) + (1,) * (xp_.ndim - 1)
        x = xp_ + U_.reshape(shape) * tp + x0.reshape(shape)  # x = x′ + Ut + x′_o
        return _u_at(u, x, tp) - U_.reshape(shape)  # u′ = u − U
    return u_prime


class AccelerationTerms(NamedTuple):
    """Fluid-particle acceleration ``a`` = ``local`` ∂u/∂t + ``advective`` (u·∇)u [m/s²] (Eq. 3.9).
    Book: §3.3, Eq. (3.9). Label: analytic.
    """

    a: np.ndarray
    local: np.ndarray
    advective: np.ndarray


def acceleration(u: Callable, x, t: float, h: float = 1e-4, ht: float | None = None) -> AccelerationTerms:
    """Acceleration of the fluid particle at x, t and its local/advective split: Du/Dt = ∂u/∂t + (u·∇)u, Eq. (3.9).

    Book: §3.3, Eq. (3.9): ∂u/∂t + (u·∇)u = (Du/Dt) in Oxyz = (Du′/Dt′) in O′x′y′z′ = ∂u′/∂t′ + (u′·∇′)u′ — the sum
    is the same in every Galilean frame, the split is not (Fig. 3.2). Eq. (3.5) applied component-wise.

    Parameters
    ----------
    u : callable u(x, t) [m/s];  x : point(s) [m] (d,) or (d, N);  t [s];  h, ht : stencil steps [m], [s]

    Returns
    -------
    AccelerationTerms(a, local, advective) — arrays of the shape of x [m/s²].

    Validation: V1 a is equal in both frames of ``galilean_transform`` at 200 points of the cylinder flow and of a
    random unsteady field while local and advective differ; V2 sympy proof of (3.9) for a general field; steady flow
    has local = 0. Label: analytic, symbolic.
    """
    loc, adv, tot = material_derivative_terms(lambda X, T: _u_at(u, X, T), u, x, t, h, ht)  # (3.5) with F = u
    return AccelerationTerms(_F(tot), _F(loc), _F(adv))


# ======================================================================================================================
# §3.4 Motion of a small fluid element
# ======================================================================================================================
def velocity_gradient_at(u: Callable, x, t: float = 0.0, h: float = 1e-4) -> np.ndarray:
    """Velocity gradient G[i, j] = ∂u_i/∂x_j at x by explicit 2nd-order central differences (Eq. (3.10)'s tensor).

    Book: §3.4, Eq. (3.10): du_i = (∂u_i/∂x_j) dx_j — the "velocity gradient tensor". Units: 1/s.

    Parameters
    ----------
    u : callable u(x, t);  x : point (d,) or points (d, N) [m];  t [s];  h : step [m]

    Returns
    -------
    ndarray (d, d) or (d, d, N).

    Validation: V1 vs the sympy Jacobian (``VectorField.grad_fn``); V3 order 2. Label: analytic, converged.
    """
    x_ = _F(x)
    d = x_.shape[0]
    cols = []
    for j in range(d):
        e = _unit(d, j, x_.ndim) * h
        cols.append((_u_at(u, x_ + e, t) - _u_at(u, x_ - e, t)) / (2.0 * h))  # ∂u_i/∂x_j (column j)
    return np.stack(cols, axis=1)


def relative_velocity(G, dx) -> np.ndarray:
    """Velocity of a neighbouring point relative to x to first order: du_i = (∂u_i/∂x_j) dx_j, Eq. (3.10).

    Book: §3.4, Fig. 3.9, Eq. (3.10). Units: G [1/s], dx [m] (shape (d,) or (d, N)) → du [m/s].
    Assumptions: |dx| small (the O(|dx|²) Taylor remainder is dropped).
    Validation: V1/V3 exact u(x + dx) − u(x) − G·dx = O(|dx|²) (slope 2). Label: analytic, converged.
    """
    return np.einsum("ij,j...->i...", _F(G), _F(dx))  # Eq. (3.10): du_i = G_ij dx_j


def vorticity_from_gradient(G):
    """Vorticity ω = vector of the rotation tensor R = G − Gᵀ (Eqs. (3.13), (3.15)): ω_k = −½ ε_ijk R_ij.

    Book: §3.4, Eq. (3.15) (= (2.26), (2.27)): R_ij = −ε_ijk ω_k; components (3.16). Always three components: a 2 × 2
    G (plane flow) returns (0, 0, ω₃) with ω₃ = R₂₁ = ∂u₂/∂x₁ − ∂u₁/∂x₂. Trailing point axes (d, d, N) are allowed
    (→ (3, N)). Units: 1/s.

    Validation: V1 u = b × x gives 2b (discriminates R from ½R); equals the curl of the field at 50 points.
    Label: analytic.
    """
    G_ = _F(G)
    R = rotation_tensor(G_)  # Eq. (3.13): R = G − Gᵀ
    if G_.shape[0] == 2:
        z = np.zeros_like(R[1, 0])
        return np.stack([z, z, R[1, 0]])  # (0, 0, ω₃) with ω₃ = R₂₁  (Eq. 3.16, third component)
    return -0.5 * np.einsum("ijk,ij...->k...", levi_civita(), R)  # Eq. (3.15) inverted: ω_k = −½ ε_ijk R_ij


def vorticity(u: Callable, x, t: float = 0.0, h: float = 1e-4):
    """Vorticity ω = ∇×u at x, Eq. (3.16), via :func:`velocity_gradient_at` and :func:`vorticity_from_gradient`.

    Book: §3.4, Eqs. (3.15)–(3.16). Units: 1/s; three components ((0, 0, ω₃) for a plane field). Label: analytic,
    converged.
    """
    return vorticity_from_gradient(velocity_gradient_at(u, x, t, h))


def _normalise(n) -> np.ndarray:
    n_ = _F(n)
    return n_ / np.sqrt(np.sum(n_ ** 2, axis=0))


def linear_strain_rate(G, n):
    """Rate of stretching per unit length of a material line along the unit direction n: n·S·n.

    Book: §3.4, Fig. 3.10: (1/δx₁) D(δx₁)/Dt = lim (A′B′ − AB)/(AB dt) = ∂u₁/∂x₁ = S₁₁ (no summation on a Greek
    index, S_ηη = ∂u_η/∂x_η); for a general direction n this is n·S·n (D08). n need not be unit (it is normalised);
    shape (d,) or (d, N). Units: 1/s.

    Validation: V1 n = e₁ gives S₁₁; V4 vs tracked material segments (``linear_flow_map``): (1/ℓ)dℓ/dt → n·S·n as
    dt → 0. Label: analytic.
    """
    S = strain_rate_tensor(_F(G))
    nn = _normalise(n)
    return as_scalar_if_0d(np.einsum("i...,ij,j...->...", nn, S, nn))  # (1/ℓ) Dℓ/Dt = n·S·n


def shear_strain_rate(G, n1, n2, tol: float = 1e-9):
    """Half the rate at which two perpendicular material lines along n1 and n2 close their right angle: n1·S·n2.

    Book: §3.4, Fig. 3.11: ½ D(α + β)/Dt = ½(∂u₁/∂x₂ + ∂u₂/∂x₁) = S₁₂ = S₂₁ for n1 = e₁, n2 = e₂; in general
    n1·S·n2 (D09). Units: 1/s.

    Raises
    ------
    ValueError if n1 and n2 are not perpendicular (|cos| > tol).

    Validation: V1 e₁, e₂ gives S₁₂; V4 −½ d(angle)/dt of two tracked segments → n1·S·n2. Label: analytic.
    """
    a, b = _normalise(n1), _normalise(n2)
    if np.any(np.abs(np.sum(a * b, axis=0)) > tol):
        raise ValueError("shear_strain_rate needs perpendicular directions n1 ⟂ n2")
    S = strain_rate_tensor(_F(G))
    return as_scalar_if_0d(np.einsum("i...,ij,j...->...", a, S, b))  # ½ D(α+β)/Dt = n1·S·n2


def volumetric_strain_rate(G) -> float:
    """Volumetric (bulk) strain rate (1/δV) D(δV)/Dt = ∂u_i/∂x_i = S_ii = ∇·u, Eq. (3.14).

    Book: §3.4, Eq. (3.14) (proof Exercise 3.18, written out as D11). The first invariant of S: independent of the
    orientation of the axes. Units: 1/s.

    Validation: V1 trace; V4 tracked box volume rate; V7 invariant under 50 random rotations. Label: analytic.
    """
    G_ = _F(G)
    return float(np.trace(G_))  # Eq. (3.14): S_ii = ∂u_i/∂x_i


def material_volume_ratio(G, t):
    """Volume (area in 2-D) of a material element after time t in the linear flow u = G·x, relative to its start:
    det e^{Gt} (= e^{t tr G}, Jacobi's formula) — the finite-time version of (3.14).

    Book: §3.4, Eq. (3.14) integrated in time (ours; the book states the rate only). Units: dimensionless; t [s]
    (float or array).

    Validation: V1 det expm(Gt) = e^{t tr G}; its log-derivative is tr G = S_ii. Label: analytic.
    """
    G_ = _F(G)
    ts = np.atleast_1d(_F(t))
    out = np.array([np.linalg.det(expm(G_ * ti)) for ti in ts])  # δV(t)/δV(0) = det e^{Gt}
    return as_scalar_if_0d(out[0]) if np.ndim(t) == 0 else out


def element_rotation_rate(G):
    """Angular velocity of a small fluid element = ½ω (half the vorticity).

    Book: §3.4, Fig. 3.11: ½ D(−α + β)/Dt = ½(−∂u₁/∂x₂ + ∂u₂/∂x₁) = −R₁₂/2 = R₂₁/2 — "ω and R_ij represent twice the
    fluid element rotation rate". Three components ((0, 0, ½ω₃) for a 2 × 2 G). Units: rad/s.

    Validation: V1 the average turning rate of any perpendicular pair of material lines equals ½ω₃ (36 angles,
    5 presets); solid body: ½ω = Ω. Label: analytic.
    """
    return 0.5 * _F(vorticity_from_gradient(G))  # spin = ½ω


def material_line_rotation_rate(G, theta):
    """Angular velocity dθ/dt of a single material line at angle θ in a plane flow: e_θ·G·e(θ).

    Book: §3.4–3.5 (Fig. 3.11 lines AB/BC; §3.5 "the angular velocity of line element AB is −γ, and that of BC is
    zero"). With e = (cos θ, sin θ) and e_θ = (−sin θ, cos θ): dθ/dt = e_θ·G·e (ours, D14). For the parallel shear
    G = [[0, γ], [0, 0]] this is −γ sin²θ. Units: rad/s; θ [rad] (float or array); G 2 × 2 (or the (1, 2) block).

    Validation: V1 shear −γ sin²θ; solid body ω₀ for every θ; the perpendicular-pair average = ½ω₃ for all θ.
    Label: analytic.
    """
    G_ = _F(G)[:2, :2]
    th = _F(theta)
    c, s = np.cos(th), np.sin(th)
    rate = -s * (G_[0, 0] * c + G_[0, 1] * s) + c * (G_[1, 0] * c + G_[1, 1] * s)  # e_θ·G·e(θ)
    return as_scalar_if_0d(rate)


def perpendicular_pair_rotation_rate(G, theta):
    """Average turning rate of the perpendicular pair of material lines at θ and θ + π/2 — equal to ½ω₃ for every θ.

    Book: §3.4 (½ D(−α + β)/Dt) and §3.5 ("the average value does not depend on which two mutually perpendicular
    elements … are chosen" — shown in D14: −γ(sin²θ + cos²θ)/2 = −γ/2). Units: rad/s. Label: analytic.
    """
    th = _F(theta)
    return as_scalar_if_0d(0.5 * (_F(material_line_rotation_rate(G, th)) + _F(material_line_rotation_rate(G, th + np.pi / 2))))


def vorticity_in_rotating_frame(omega, Omega):
    """Vorticity seen from a frame rotating at angular velocity Ω: ω′ = ω − 2Ω.

    Book: §3.4 ("ω and R_ij depend on the frame of reference … zero in a co-rotating frame", Exercise 3.19; written
    out as D13): u′ = u − Ω × x and ∇×(Ω × x) = 2Ω (ch02 Example 2.3). Scalars are read as z-components. In
    geophysical language ω′ is the relative vorticity ζ and 2Ω the planetary vorticity (f at the pole), Ch. 13.
    Units: 1/s.

    Validation: V1 numeric curl of u − Ω × x equals ω − 2Ω; zero when Ω = ω/2. Label: analytic.
    """
    return as_scalar_if_0d(_F(omega) - 2.0 * _F(Omega))  # ω′ = ω − 2Ω


def _cumtrapz_from(f: np.ndarray, x: np.ndarray, i0: int, axis: int) -> np.ndarray:
    """∫_{x[i0]}^{x[i]} f dx along ``axis`` (trapezoid), zero at index i0."""
    from scipy.integrate import cumulative_trapezoid

    c = cumulative_trapezoid(f, x, axis=axis, initial=0.0)
    return c - np.take(c, [i0], axis=axis)


def velocity_potential_2d(u: Callable, grid, x_ref=(0.0, 0.0), t: float = 0.0):
    """Velocity potential φ (u = ∇φ, Eq. (3.17)) of a plane field by line integrals of u·ds on a grid, from the node
    nearest ``x_ref`` along two different staircase routes — and how much the two routes disagree.

    Book: §3.4, Eq. (3.17): an irrotational flow can be written u = ∇φ. The converse needs a simply connected region
    (D18): around the line vortex (3.25) the two routes differ by the circulation 2πB and φ = Bθ is multivalued.
    This function therefore **reports** the path dependence instead of asserting irrotational ⇒ potential.

    Parameters
    ----------
    u : callable u(x, t) (plane field) [m/s]
    grid : a ``core.grids.Grid2D`` (``grid2d``) or a pair (x, y) of increasing 1-D node arrays [m] (avoid singular
        points — e.g. offset the grid so that no node sits on a vortex axis)
    x_ref : reference point (φ = 0 at the nearest node) [m];  t : time [s]

    Returns
    -------
    (phi, path_dependence): φ [m²/s] on the grid ``[j, i]`` = (y, x) by route 1 ("along x on the reference row, then
    along y"), and max |φ_route1 − φ_route2| [m²/s] with route 2 "along y on the reference column, then along x".
    Cumulative trapezoid rule (error O(h²)).

    Around a singular axis put ``x_ref`` at a grid corner, far from the core: routes that pass within a few spacings
    of the axis carry a large trapezoid error (the path dependence then over- or under-shoots 2πB by several %).

    Validation: V1 φ = x² − y² recovered up to a constant (1e-6 at 201 nodes); path dependence ≈ 2πB (up to the
    O(h²) error near the axis) for a grid around a line vortex. Label: analytic.
    """
    if hasattr(grid, "x") and hasattr(grid, "y"):
        x_, y_ = _F(grid.x), _F(grid.y)
    else:
        x_, y_ = (_F(v) for v in grid)
    X, Y = np.meshgrid(x_, y_, indexing="xy")
    U = _u_at(u, np.stack([X.ravel(), Y.ravel()]), t).reshape(2, *X.shape)
    i0 = int(np.argmin(np.abs(x_ - float(x_ref[0]))))
    j0 = int(np.argmin(np.abs(y_ - float(x_ref[1]))))
    # route 1: along the reference row y = y_j0 in x, then along each column in y
    phi1 = _cumtrapz_from(U[0, j0, :], x_, i0, 0)[None, :] + _cumtrapz_from(U[1], y_, j0, 0)  # ∫u dx + ∫v dy
    # route 2: along the reference column x = x_i0 in y, then along each row in x
    phi2 = _cumtrapz_from(U[1, :, i0], y_, j0, 0)[:, None] + _cumtrapz_from(U[0], x_, i0, 1)  # ∫v dy + ∫u dx
    return phi1, float(np.nanmax(np.abs(phi1 - phi2)))


class RelativeVelocitySplit(NamedTuple):
    """Eq. (3.19): ``du`` = G·dx = ``du_strain`` S·dx + ``du_rot`` ½ω × dx [m/s].
    Book: §3.4, Eq. (3.19). Label: analytic.
    """

    du: np.ndarray
    du_strain: np.ndarray
    du_rot: np.ndarray


def relative_velocity_split(G, dx) -> RelativeVelocitySplit:
    """Relative velocity of a neighbour = pure deformation + rigid rotation at ω/2, Eq. (3.19).

    Book: §3.4, Eq. (3.19): du_i = (S_ij − ½ε_ijk ω_k) dx_j = S_ij dx_j + ½(ω × dx)_i — from (3.10), (3.11) and
    (3.15) (the book says "(3.11) through (3.14)"; the step uses (3.15), analysis §9 item 5). In 2-D the rotation
    part is ½ω₃ e₃ × dx = ½ω₃(−dx₂, dx₁).

    Parameters
    ----------
    G : (d, d) [1/s];  dx : (d,) or (d, N) [m]

    Returns
    -------
    RelativeVelocitySplit(du, du_strain, du_rot) [m/s]; ``du`` is G·dx computed directly from (3.10), the other two
    from S and ω separately, so ``du_strain + du_rot == du`` is a real check.

    Validation: V1 sum = G·dx; du_rot = ½ω × dx; du_rot ⟂ dx; V2 index identity via ``expand_indices``.
    Label: analytic, symbolic.
    """
    G_, dx_ = _F(G), _F(dx)
    du = relative_velocity(G_, dx_)  # Eq. (3.10)
    du_strain = np.einsum("ij,j...->i...", strain_rate_tensor(G_), dx_)  # S_ij dx_j
    w = vorticity_from_gradient(G_)
    if G_.shape[0] == 2:
        du_rot = 0.5 * float(w[2]) * np.stack([-dx_[1], dx_[0]])  # ½ ω₃ e₃ × dx
    else:
        wv = _F(w).reshape((3,) + (1,) * (dx_.ndim - 1))
        du_rot = 0.5 * np.cross(np.broadcast_to(wv, dx_.shape), dx_, axis=0)  # Eq. (3.19): ½(ω × dx)_i
    return RelativeVelocitySplit(du, du_strain, du_rot)


def principal_strain_rates(G):
    """Principal strain rates (eigenvalues S̄_αα of S, ascending) and principal axes (unit eigenvectors as columns).

    Book: §3.4, Eq. (3.20): in the principal frame (overbar) dū = S̄·dx̄ with S̄ = diag(S̄₁₁, S̄₂₂, S̄₃₃); ch02 §2.11
    ``principal_axes`` (det of the axes matrix = +1, so it is the direction-cosine matrix C of (2.12)).
    Sign convention (ours, for reproducible figures and parity rows): the axis of the largest rate has a positive first
    non-zero component; det = +1 is then restored by flipping the first column if needed. Shear γ = 1: λ = (−0.5, 0.5),
    the +0.5 axis is (1, 1)/√2.
    Units: 1/s. Validation: V1 shear flow ±γ/2 at ±45° (§3.5); S·C = C·diag(λ). Label: analytic.
    """
    lam, C = principal_axes(strain_rate_tensor(_F(G)))
    C = C.copy()
    v = C[:, -1]
    k = int(np.argmax(np.abs(v) > 1e-12))
    if v[k] < 0:
        C[:, -1] = -v
    if np.linalg.det(C) < 0:
        C[:, 0] = -C[:, 0]
    return lam, C


def strain_velocity_principal(G, dx):
    """The strain part of (3.19) in the principal frame: dū_α = S̄_αα dx̄_α (no sum), Eqs. (3.20)–(3.21).

    Book: §3.4, Eqs. (3.20), (3.21). dx̄ = Cᵀ dx and dū = Cᵀ(S·dx) with C the principal axes (ch02 passive
    convention x′ = Cᵀx). Returns (du_bar, dx_bar) [m/s], [m]; ``du_bar == lam * dx_bar`` component by component.

    Validation: V1 du_bar = diag(S̄)·dx_bar to 1e-14 for random G, dx. Label: analytic.
    """
    lam, C = principal_strain_rates(G)
    dx_ = _F(dx)
    S = strain_rate_tensor(_F(G))
    dx_bar = np.einsum("ji,j...->i...", C, dx_)  # dx̄ = Cᵀ dx
    du_bar = np.einsum("ji,j...->i...", C, np.einsum("ij,j...->i...", S, dx_))  # dū = Cᵀ S dx = S̄ dx̄ (3.21)
    return du_bar, dx_bar


def deform_circle(G, t, n: int = 72, radius: float = 1.0, center=(0.0, 0.0)) -> np.ndarray:
    """Material points on a circle of radius ``radius`` (in the x₁–x₂ plane) carried by u = G·x for time t.

    Book: §3.4, Fig. 3.13 (a small sphere becomes an ellipsoid on the principal axes of S). Positions relative to the
    moving centre (the centre itself follows u(center) — omitted: we show the relative motion). Returns (2, n)
    (or (3, n) for a 3 × 3 G, the circle at x₃ = 0), not closed (append the first column to draw).
    Units: G [1/s], t [s], radius [m]. Label: analytic.
    """
    G_ = _F(G)
    ph = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    pts = radius * np.stack([np.cos(ph), np.sin(ph)])
    if G_.shape[0] == 3:
        pts = np.vstack([pts, np.zeros(n)])
    c = _F(center)
    moved = linear_flow_map(G_, t) @ pts
    moved[:2] += c[:2, None]
    return moved


def deform_sphere(G, t, n_theta: int = 24, n_phi: int = 48, radius: float = 1.0) -> np.ndarray:
    """Material points on a sphere of radius ``radius`` carried by u = G·x (3 × 3) for time t: array (3, N) with
    N = n_theta·n_phi (reshape to (3, n_theta, n_phi) for a surface plot; θ polar, φ azimuth, θ varies slowest).
    Book: §3.4, Fig. 3.13. Label: analytic."""
    th = np.linspace(0.0, np.pi, n_theta)
    ph = np.linspace(0.0, 2.0 * np.pi, n_phi)
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    pts = radius * np.stack([np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)]).reshape(3, -1)
    return linear_flow_map(_F(G), t) @ pts


_ELLIPSE_METHODS = {"first_order": "first_order", "strain_only": "strain_only", "strain": "strain_only",
                    "exact": "exact", "full": "exact"}


def strain_ellipse_axes(G, t, method: str = "first_order", radius: float = 1.0):
    """Semi-axes and directions of the ellipse (ellipsoid) that a small circle (sphere) of radius ``radius`` becomes
    after time t.

    Book: §3.4, Fig. 3.13: "a spherical fluid element around O becomes an ellipsoid whose axes are the principal axes
    of the strain-rate tensor S" — a first-order-in-t statement (D16).

    method
        ``"first_order"`` — (1 + S̄_αα t)·radius along the eigenvectors of S (the book's statement, Eq. (3.21) over
        one small step);
        ``"strain_only"`` — the pure-strain part exactly: e^{S̄_αα t}·radius along the eigenvectors of S (expm(St));
        ``"exact"`` — the actual tracer ellipse: singular values and left singular vectors of e^{Gt} (the rotation
        ½ω turns the ellipse; for finite t in shear its axes drift away from the 45° principal axes).
        (Aliases: "strain" = "strain_only", "full" = "exact".)

    Returns
    -------
    (semi_axes, directions) — semi-axes [m] (unit initial radius by default) in **descending** order;
    ``directions[:, k]`` the unit axis of semi-axis k.

    Validation: V1 pure strain: all three methods agree to O(t²); "strain_only" = e^{λt} exactly; product of the
    semi-axes = det e^{Gt}. Label: analytic.
    """
    G_ = _F(G)
    m = _ELLIPSE_METHODS.get(method)
    if m is None:
        raise ValueError('method must be "first_order", "strain_only" or "exact"')
    if m == "exact":
        Umat, s, _ = np.linalg.svd(linear_flow_map(G_, t))
        return s * radius, Umat  # already descending
    lam, C = principal_strain_rates(G_)
    ax = (np.exp(lam * t) if m == "strain_only" else 1.0 + lam * t) * radius  # e^{S̄t} or (1 + S̄t)
    order = np.argsort(ax)[::-1]
    return ax[order], C[:, order]


def _angle(a, b) -> float:
    return float(np.arccos(np.clip(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)), -1.0, 1.0)))


def measured_strain_rates(G, n, dt: float, n2=None) -> dict:
    """What a ruler and a protractor on a moving element measure over a short time dt, next to the formula values.

    Book: §3.4 (Figs. 3.10, 3.11, Eq. (3.14)): the tracked material segments along n and n2 (n2 ⟂ n; in 2-D the
    default is n rotated by +90°) are mapped by e^{G dt}. Measured:

    * ``stretch`` = (ℓ(dt)/ℓ(0) − 1)/dt → n·S·n (``stretch_formula``)
    * ``closing`` = (π/2 − angle between the two segments)/dt → 2 n·S·n2 (``closing_formula``, = D(α + β)/Dt)
    * ``spin`` = mean signed turning rate of the two segments → ½ω₃ (``spin_formula``; 2-D only, NaN in 3-D)
    * ``area`` = (det e^{G dt} − 1)/dt → tr G (``area_formula``; area in 2-D, volume in 3-D)

    All measured values converge to the formulas at first order in dt (the book's lim dt → 0).

    Parameters
    ----------
    G : (d, d) [1/s];  n : direction (d,) (normalised);  dt : measurement interval [s] (> 0);  n2 : second direction

    Returns
    -------
    dict of floats [1/s].

    Validation: V3 |measured − formula| = O(dt) (order 1). Label: converged.
    """
    G_ = _F(G)
    d = G_.shape[0]
    a = _normalise(n)
    if n2 is None:
        if d != 2:
            raise ValueError("give n2 (perpendicular to n) for a 3-D G")
        b = np.array([-a[1], a[0]])
    else:
        b = _normalise(n2)
    M = linear_flow_map(G_, dt)
    ma, mb = M @ a, M @ b
    S = strain_rate_tensor(G_)
    out = {"stretch": (np.linalg.norm(ma) - 1.0) / dt,  # (A′B′ − AB)/(AB dt), Fig. 3.10
           "closing": (0.5 * np.pi - _angle(ma, mb)) / dt,  # D(α + β)/Dt, Fig. 3.11
           "area": (np.linalg.det(M) - 1.0) / dt,  # (1/δV) D(δV)/Dt, Eq. (3.14)
           "stretch_formula": float(a @ S @ a), "closing_formula": float(2.0 * a @ S @ b),
           "area_formula": float(np.trace(G_))}
    if d == 2:
        turn_a = np.arctan2(a[0] * ma[1] - a[1] * ma[0], a @ ma)
        turn_b = np.arctan2(b[0] * mb[1] - b[1] * mb[0], b @ mb)
        out["spin"] = float(0.5 * (turn_a + turn_b) / dt)  # ½ D(−α + β)/Dt
        out["spin_formula"] = float(element_rotation_rate(G_)[2])
    else:
        out["spin"] = out["spin_formula"] = float("nan")
    return {k: float(v) for k, v in out.items()}
