"""Chapter 2 — Cartesian Tensors: worked examples, test fields and the plane-flow kinematics the notebook and
explainers use, plus re-exports of every reusable primitive the chapter introduced in ``fluidpy.core``.

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e (2012), Ch. 2, §§2.1–2.14, Eqs. (2.1)–(2.36), Examples 2.1–2.6.
Every equation was transcribed from the rendered page images (chapters/pages/ch02/p067–p089).

Where the physics lives
-----------------------
* ``core.tensors`` — bases, direction cosines, (2.5)–(2.8), (2.12)–(2.14), (2.15) traction, (2.16)–(2.21) δ, ε,
  cross products, (2.26)–(2.29) symmetric/antisymmetric parts, §2.11 principal axes.
* ``core.index_notation`` — the summation convention as a sympy expander ((2.2), (2.17), (2.19), (2.21), (2.36)).
* ``core.grids`` / ``core.operators`` — the project grid layout and ∇, ∇·, ∇× stencils ((2.22)–(2.25)).
* ``core.fields`` — sympy-defined vector fields with exact div/curl (Example 2.3 and the theorem checks).
* ``core.integral_theorems`` — Gauss (2.30)–(2.33), Stokes (2.34)–(2.35), Examples 2.5–2.6.

Conventions (recorded in ``knowledge/notation.md``; ⚠️ callouts in the notebook)
------------------------------------------------------------------------------
* Passive direction cosines C_ij = e_i·e'_j (columns = new axes): x' = Cᵀx, τ' = CᵀτC. scipy/Wikipedia R = Cᵀ applied.
* Traction contracts the first index of τ (f_i = τ_ji n_j); tensor divergence contracts the second (∂τ_ij/∂x_j).
* Velocity gradient G[i, j] = ∂u_i/∂x_j; the vector of its antisymmetric part is ½∇×u.
* Double dot: book A:B = A_ij B_ji vs Frobenius A_ij B_ij — always passed explicitly.
* Example 2.4's Γ is taken as the off-diagonal element S₁₂ (the book's matrix S = [[0, Γ], [Γ, 0]]); for a flow
  u₁(x₂) alone S₁₂ = ½ du₁/dx₂, so the book's line "2S₁₂ = du₁/dx₂ = Γ" is inconsistent with its matrix
  (analysis §9, item 4). ``example_2_4`` documents the factor.
* Book-quoted numbers are not in this file (git-ignored ``tests/book_values_ch02.json``).
"""
from __future__ import annotations

from typing import Callable, Sequence

import numpy as np
import sympy as sp
from scipy.linalg import expm

from .core._util import as_scalar_if_0d
from .core.fields import ScalarField, VectorField  # noqa: F401
from .core.grids import Grid2D, Grid3D, axis_of_direction, evaluate_field, grid, grid2d, spacing  # noqa: F401
from .core.index_notation import (  # noqa: F401
    classify_indices, comma_to_partial, coordinates, expand_indices, expand_indices_str, field_functions,
    rename_dummy, symbol_matrix, symbol_tensor, symbol_vector, tensor_order,
)
from .core.integral_theorems import (  # noqa: F401
    Loop, StokesCheck, Surface, TiledFlux, boundary_tangent, circulation, curl_flux, divergence_theorem_box,
    divergence_theorem_rect2d, divergence_theorem_sphere, divergence_theorem_tiled, fd_partials,
    flux_through_box_faces, flux_through_faces,
    flux_through_sphere, gauss_gradient_box, gauss_legendre_nodes, integral_curl, integral_curl_component,
    integral_divergence, integral_divergence_2d, integral_gradient, midpoint_nodes, planar_disc, planar_loop,
    planar_rectangle, rectangle_loop, simpson_nodes, stokes_theorem_check, volume_integral_box,
)
from .core.operators import (  # noqa: F401
    curl, curl_components, directional_derivative, divergence, gradient, is_irrotational, is_solenoidal, laplacian,
    partial, second_partial, tensor_divergence, vector_gradient,
)
from .core.tensors import (  # noqa: F401
    STRESS_CUBE_FACES, angle_between, antisymmetric_from_vector, antisymmetric_part, characteristic_polynomial,
    contract, cross, cross_einsum, cube_face_tractions, diagonalize, direction_cosines, dot, dot_tensor_vector,
    double_dot, epsilon_delta_residual, independent_components, inner, invariants, inverse_transform_vector,
    is_antisymmetric, is_isotropic, is_orthogonal, is_proper_rotation, is_symmetric, kronecker_delta, levi_civita,
    mohr_circle_2d, normal_shear_stress, normal_stress_bounds, orthogonality_residual, outer, permutation_sign,
    principal_axes, random_rotation, rotation_angle, rotation_matrix_2d, rotation_matrix_3d, rotation_tensor,
    strain_rate_tensor, stress_component_meaning, symmetric_double_contraction, symmetric_part, tensor_product,
    tetrahedron_face_areas, trace, traction, traction_components, transform_tensor, transform_vector,
    transforms_as_tensor, transforms_as_vector, triple_product, unit_vectors, vector_from_antisymmetric,
    vector_from_components,
)

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731


# ----------------------------------------------------------------------------------------------------------------------
# §2.2 Example 2.1 — Cartesian → polar components
# ----------------------------------------------------------------------------------------------------------------------
def polar_components(u1, u2, theta):
    """Polar components (u_r, u_θ) of the plane vector (u₁, u₂) at polar angle θ.

    Book: §2.2, Example 2.1: u_r = u₁ cos θ + u₂ sin θ, u_θ = −u₁ sin θ + u₂ cos θ — Eq. (2.5) with j ∈ {r, θ} and
    C = [[cos θ, −sin θ], [sin θ, cos θ]] (the polar frame is the Cartesian frame rotated by θ).

    Parameters
    ----------
    u1, u2 : float or array_like — Cartesian components (any unit)
    theta : float or array_like — polar angle of the point [rad]

    Returns
    -------
    (u_r, u_theta) : same unit as the inputs.

    Validation: V1 equals ``transform_vector(u, rotation_matrix_2d(θ))``; θ = 0 identity; θ = π/2 → (u₂, −u₁);
    |u| preserved; V6 Example 2.1 formulas. Label: analytic.
    """
    u1, u2, th = _F(u1), _F(u2), _F(theta)
    u_r = u1 * np.cos(th) + u2 * np.sin(th)  # Example 2.1: u_r = u1 cos θ + u2 cos(π/2 − θ)
    u_th = -u1 * np.sin(th) + u2 * np.cos(th)  # Example 2.1: u_θ = u1 cos(θ + π/2) + u2 cos θ
    return as_scalar_if_0d(u_r), as_scalar_if_0d(u_th)


def cartesian_from_polar(u_r, u_theta, theta):
    """Inverse of :func:`polar_components`: u₁ = u_r cos θ − u_θ sin θ, u₂ = u_r sin θ + u_θ cos θ (Eq. (2.7) form).

    Validation: V1 round trip with ``polar_components`` (1e-14). Label: analytic.
    """
    ur, ut, th = _F(u_r), _F(u_theta), _F(theta)
    return as_scalar_if_0d(ur * np.cos(th) - ut * np.sin(th)), as_scalar_if_0d(ur * np.sin(th) + ut * np.cos(th))


def polar_components_deg(u1, u2, theta_deg):
    """:func:`polar_components` with the angle in degrees (interface convenience)."""
    return polar_components(u1, u2, np.deg2rad(_F(theta_deg)))


# ----------------------------------------------------------------------------------------------------------------------
# §2.6 Example 2.2 — traction in a shear flow; 2-D helpers for the explainer
# ----------------------------------------------------------------------------------------------------------------------
def shear_flow_stress(a) -> np.ndarray:
    """Viscous stress tensor of a parallel channel flow at a point: τ = [[0, a], [a, 0]] (a > 0 in one half, < 0 in the other).

    Book: §2.6, Example 2.2. ``a`` in Pa.
    """
    a = float(a)
    return np.array([[0.0, a], [a, 0.0]])  # Example 2.2


def traction_2d(tau, phi) -> dict:
    """Traction on the plane whose outward normal makes angle φ with the x₁-axis, for any 2 × 2 stress tensor.

    Book: §2.6, Eq. (2.15) with n = (cos φ, sin φ); the split into normal and shear stress as in Example 2.2 and the
    rotated-frame reading τ'₁₁ = f·e'₁, τ'₁₂ = f·e'₂ (Eq. (2.12) with the 2 × 2 C of the example).

    Parameters
    ----------
    tau : array_like, shape (2, 2)  [Pa]
    phi : float — angle of the normal [rad]

    Returns
    -------
    dict with ``n`` (2,), ``f`` (2,) [Pa], ``sigma_n`` = f·n [Pa], ``tau_s`` = f·s (signed shear along s = (−sin φ,
    cos φ), the rotated 2'-axis) [Pa], ``tau_s_mag`` = |f − σ_n n| [Pa], ``angle_rad`` = direction of f in [0, 2π).

    For *any* 2 × 2 τ (symmetric or not) σ_n = τ'₁₁ and τ_s = τ'₁₂ of τ' = CᵀτC with C = ``rotation_matrix_2d(φ)``:
    f_i = τ_ji n_j = τ_ji C_j1, so f·e'₁ = C_i1 C_j1 τ_ji = τ'₁₁ and f·e'₂ = C_i2 C_j1 τ_ji = τ'₁₂ by (2.12) — no symmetry
    is used. (Symmetry would only be needed to also read τ_s as τ'₂₁.)

    Validation: V1 σ_n n + τ_s s == f; σ_n = τ'₁₁ and τ_s = τ'₁₂ of ``transform_tensor`` for symmetric and
    non-symmetric τ; V6 Example 2.2. Label: analytic.
    """
    t = _F(tau)
    ph = float(phi)
    n = np.array([np.cos(ph), np.sin(ph)])  # Example 2.2: n = (cos φ, sin φ)
    s = np.array([-np.sin(ph), np.cos(ph)])  # the rotated 2'-axis, e'_2
    f = traction(t, n)  # Eq. (2.15): f_i = τ_ji n_j
    sigma_n = float(f @ n)
    tau_s = float(f @ s)
    return {"n": n, "s": s, "f": f, "sigma_n": sigma_n, "tau_s": tau_s, "tau_s_mag": abs(tau_s),
            "magnitude": float(np.linalg.norm(f)), "angle_rad": float(np.arctan2(f[1], f[0]) % (2.0 * np.pi))}


def example_2_2(a, phi=np.pi / 6) -> dict:
    """Example 2.2 worked out: traction on an element whose normal is φ from the flow direction in τ = [[0, a], [a, 0]].

    Book: §2.6, Example 2.2 — by (2.15): f = (a sin φ, a cos φ), |f| = |a|, direction θ with sin θ = f₂/|f|,
    cos θ = f₁/|f| (60° for a > 0, 240° for a < 0 at φ = 30°); by (2.12) with C = [[cos φ, −sin φ], [sin φ, cos φ]]:
    τ'₁₁ = 2 sin φ cos φ a = a sin 2φ (normal stress), τ'₁₂ = (cos²φ − sin²φ) a = a cos 2φ (shear stress).

    Parameters
    ----------
    a : float — shear stress level [Pa] (sign selects the half of the channel)
    phi : float — angle of the outward normal from the x₁-axis [rad] (book: 30°)

    Returns
    -------
    dict: ``tau``, ``n``, ``f`` [Pa], ``magnitude`` [Pa], ``angle_rad``, ``angle_deg`` (direction of f in [0, 360°)),
    ``C`` (2 × 2 direction cosines), ``tau_rot`` = CᵀτC, ``sigma_n`` = τ'₁₁, ``tau_s`` = τ'₁₂,
    ``sigma_n_formula`` = a sin 2φ, ``tau_s_formula`` = a cos 2φ.

    Validation: V1 the (2.15) and (2.12) routes agree (f·e'₁ = τ'₁₁, f·e'₂ = τ'₁₂) for 50 φ; |f| = |a|; general forms
    a sin 2φ, a cos 2φ; V6 book values at 30° (private JSON); V5 Mohr formula. Label: analytic.
    """
    a = float(a)
    tau = shear_flow_stress(a)
    d = traction_2d(tau, phi)
    C = rotation_matrix_2d(phi)  # Example 2.2: C_ij = [[cos φ, −sin φ], [sin φ, cos φ]]
    tau_rot = transform_tensor(tau, C)  # Eq. (2.12): τ' = Cᵀ τ C
    return {"tau": tau, "n": d["n"], "f": d["f"], "magnitude": d["magnitude"], "angle_rad": d["angle_rad"],
            "angle_deg": float(np.rad2deg(d["angle_rad"])), "C": C, "tau_rot": tau_rot,
            "sigma_n": float(tau_rot[0, 0]), "tau_s": float(tau_rot[0, 1]),
            "sigma_n_formula": a * np.sin(2.0 * float(phi)), "tau_s_formula": a * np.cos(2.0 * float(phi))}


def stress_vs_angle(tau, phi) -> dict:
    """σ_n(φ) and τ_s(φ) for a symmetric 2 × 2 stress over an array of normal angles (the C05/E2 curves).

    Book: §2.6 Example 2.2 (a sin 2φ, a cos 2φ for pure shear) generalised: σ_n = n·τ·n, τ_s = s·τ·n with
    s = (−sin φ, cos φ). Returns dict(phi, sigma_n, tau_s) [Pa].

    Validation: V1 pure shear reproduces a sin 2φ / a cos 2φ; extremes of σ_n equal the eigenvalues. Label: analytic.
    """
    t = _F(tau)
    ph = _F(phi)
    n = np.stack([np.cos(ph), np.sin(ph)])
    s = np.stack([-np.sin(ph), np.cos(ph)])
    f = np.einsum("ji,j...->i...", t, n)  # Eq. (2.15) for every angle
    return {"phi": ph, "sigma_n": np.einsum("i...,i...->...", f, n), "tau_s": np.einsum("i...,i...->...", f, s)}


# ----------------------------------------------------------------------------------------------------------------------
# §2.9 Example 2.3 and other test fields (sympy twins carried by core.fields.VectorField)
# ----------------------------------------------------------------------------------------------------------------------
def radial_field(a=1.0, dim: int = 3) -> VectorField:
    """u = a x — the field that diverges from the origin: ∇·u = 3a (2a in the plane), ∇×u = 0.

    Book: §2.9, Example 2.3 (the book prints "a x₂ e₂" twice where the last term is a x₃ e₃). ``a`` may be a number
    or a sympy Symbol; the callable uses ``a`` numerically (1.0 for a symbol), ``.div_expr`` keeps it symbolic.

    Validation: V2 div = 3a, curl = 0 symbolically; V1 on the grid with ``core.operators`` (exact: linear field);
    V6 factor 3 (private JSON). Label: symbolic, analytic.
    """
    X = coordinates(dim)
    a_sym = a if isinstance(a, sp.Basic) else sp.Symbol("a", real=True)
    params = {a_sym: 1.0} if isinstance(a, sp.Basic) else {a_sym: float(a)}
    return VectorField([a_sym * x for x in X], X, params=params, name="a x")


def solid_body_rotation_field(b=(0.0, 0.0, 1.0), dim: int = 3) -> VectorField:
    """u = b × x — solid-body rotation about the axis b (angular velocity |b|): ∇·u = 0, ∇×u = 2b.

    Book: §2.9, Example 2.3: u = (b₂x₃ − b₃x₂) e₁ + (b₃x₁ − b₁x₃) e₂ + (b₁x₂ − b₂x₁) e₃. In the plane (``dim=2``) b
    is a scalar b₃ and u = (−b₃ x₂, b₃ x₁) with (∇×u)₃ = 2b₃. Ch. 3: vorticity = 2 × angular velocity.

    Parameters
    ----------
    b : 3-vector (numbers or sympy symbols) or a scalar b₃ for ``dim=2``

    Validation: V2 div = 0, curl = 2b symbolically; V1 Stokes: circulation around a circle of radius R ⊥ b is
    2|b|πR²; V6 factor 2 (private JSON). Label: symbolic, analytic.
    """
    X = coordinates(dim)
    if dim == 2:
        b3 = b if isinstance(b, sp.Basic) else sp.Symbol("b3", real=True)
        params = {b3: 1.0} if isinstance(b, sp.Basic) else {b3: float(np.ravel(b)[-1])}
        return VectorField([-b3 * X[1], b3 * X[0]], X, params=params, name="b × x (plane)")
    bs = list(b)
    syms = [bi if isinstance(bi, sp.Basic) else sp.Symbol(f"b{i + 1}", real=True) for i, bi in enumerate(bs)]
    params = {s_: (1.0 if isinstance(bi, sp.Basic) else float(bi)) for s_, bi in zip(syms, bs)}
    u = sp.Matrix(syms).cross(sp.Matrix(X))  # Example 2.3: u = b × x
    return VectorField(list(u), X, params=params, name="b × x")


def shear_field(Gamma=1.0) -> VectorField:
    """Plane simple shear u = (Γ x₂, 0): ∇·u = 0, (∇×u)₃ = −Γ — straight streamlines with curl (Example 2.4's flow).

    Book: §2.11 Example 2.4 (u₁(x₂) with du₁/dx₂ = Γ); §2.13 (a shear flow has circulation). Units: Γ [1/s].
    Label: symbolic, analytic.
    """
    X = coordinates(2)
    G = sp.Symbol("Gamma", real=True)
    return VectorField([G * X[1], sp.Integer(0)], X, params={G: float(Gamma)}, name="simple shear")


def irrotational_vortex_field(K=1.0) -> VectorField:
    """Plane irrotational vortex u_θ = K/r, i.e. u = K(−x₂, x₁)/(x₁² + x₂²): curl-free everywhere except the core.

    Book: §2.13 (the hypothesis of Stokes' theorem: u must be smooth on A). ``singular_at = (0, 0)`` lets
    ``stokes_theorem_check`` report that the circulation 2πK around any loop enclosing the core is not a curl flux.
    Units: K [m²/s]. Label: symbolic, analytic.
    """
    X = coordinates(2)
    K_ = sp.Symbol("K", real=True)
    r2 = X[0] ** 2 + X[1] ** 2
    return VectorField([-K_ * X[1] / r2, K_ * X[0] / r2], X, params={K_: float(K)}, name="irrotational vortex",
                       singular_at=(0.0, 0.0))


def point_source_field(m=1.0) -> VectorField:
    """Plane point source u = m x/(2π r²) (volume flux m per unit depth): divergence-free except at the origin.

    Book: §2.12 (a box enclosing the source has net outflux m whatever its size). Units: m [m²/s].
    """
    X = coordinates(2)
    m_ = sp.Symbol("m", real=True)
    r2 = X[0] ** 2 + X[1] ** 2
    return VectorField([m_ * X[0] / (2 * sp.pi * r2), m_ * X[1] / (2 * sp.pi * r2)], X, params={m_: float(m)},
                       name="point source", singular_at=(0.0, 0.0))


def potential_field(phi_expr=None, dim: int = 2) -> VectorField:
    """u = ∇φ for a scalar potential (default φ = x₁² − x₂² in 2-D, x₁x₂x₃ in 3-D): irrotational by construction
    (Exercise 2.20: ∇×∇φ = 0; ∮∇φ·t ds = 0).

    Book: §2.9 (gradient), §2.13 (zero circulation). ``phi_expr`` may be a sympy expression in ``coordinates(dim)``.
    """
    X = coordinates(dim)
    if phi_expr is None:
        phi_expr = X[0] ** 2 - X[1] ** 2 if dim == 2 else X[0] * X[1] * X[2]
    phi_expr = sp.sympify(phi_expr)
    return VectorField([sp.diff(phi_expr, x) for x in X], X, name="grad phi")


def smooth_test_field(dim: int = 3) -> VectorField:
    """A smooth non-polynomial field for convergence studies: u = (sin x₁ cos x₂, cos x₂ sin x₃, x₁ x₃ sin x₂)
    (3-D) or (sin x₁ cos x₂, cos x₁ sin x₂ + x₁²) (2-D) — every derivative non-trivial and known exactly.

    Book: tool for the V3 evidence of (2.22)–(2.25), (2.30)–(2.35). Label: symbolic.
    """
    X = coordinates(dim)
    if dim == 3:
        return VectorField([sp.sin(X[0]) * sp.cos(X[1]), sp.cos(X[1]) * sp.sin(X[2]), X[0] * X[2] * sp.sin(X[1])], X,
                           name="smooth test field")
    return VectorField([sp.sin(X[0]) * sp.cos(X[1]), sp.cos(X[0]) * sp.sin(X[1]) + X[0] ** 2], X, name="smooth test field")


def periodic_test_field() -> VectorField:
    """A 2π-periodic smooth field u = (sin x₁ cos x₂, cos x₂ sin x₃, sin x₁ sin x₂ cos x₃) for the ``bc="periodic"``
    stencils on [−π, π]³ (the non-periodic ``smooth_test_field`` would make the wrapped stencil look first order)."""
    X = coordinates(3)
    return VectorField([sp.sin(X[0]) * sp.cos(X[1]), sp.cos(X[1]) * sp.sin(X[2]), sp.sin(X[0]) * sp.sin(X[1]) * sp.cos(X[2])],
                       X, name="periodic test field")


def periodic_scalar_field() -> ScalarField:
    """φ = sin x₁ cos x₂ cos x₃, 2π-periodic, for periodic gradient/Laplacian convergence studies."""
    X = coordinates(3)
    return ScalarField(sp.sin(X[0]) * sp.cos(X[1]) * sp.cos(X[2]), X, name="periodic phi")


def smooth_scalar_field(dim: int = 3) -> ScalarField:
    """φ = sin x₁ cos x₂ (·e^{x₃/2} in 3-D) for gradient convergence studies (Fig. 2.7 uses φ = x² + y²/4 instead)."""
    X = coordinates(dim)
    expr = sp.sin(X[0]) * sp.cos(X[1]) * (sp.exp(X[2] / 2) if dim == 3 else 1)
    return ScalarField(expr, X, name="smooth phi")


def exact_div_curl(field) -> tuple:
    """Exact (sympy) divergence and curl of a field given as a ``VectorField`` or a list of sympy expressions.

    Book: §2.9, Eqs. (2.23), (2.25) evaluated symbolically (Example 2.3's answers 3a, 0, 0, 2b).

    Returns
    -------
    (div_expr, curl_expr) — ``curl_expr`` is a list of three expressions in 3-D, a single expression in 2-D.

    Validation: V2 Example 2.3; V1 the numeric operators converge to these values. Label: symbolic.
    """
    if not isinstance(field, VectorField):
        field = VectorField(list(field))
    return field.div_expr, field.curl_expr


def example_2_3(a=1.0, b=(0.0, 0.0, 1.0), n: int = 16, L: float = 1.0) -> dict:
    """Example 2.3 numerically and symbolically: ∇·(a x) = 3a, ∇×(a x) = 0; ∇·(b × x) = 0, ∇×(b × x) = 2b.

    Book: §2.9, Example 2.3. The grid operators are exact for these linear fields (central differences are exact
    for polynomials of degree ≤ 2), so the check is to round-off.

    Returns
    -------
    dict with the symbolic results (``div_radial``, ``curl_radial``, ``div_rotation``, ``curl_rotation``) and the
    grid maxima ``div_radial_grid`` (→ 3a), ``curl_radial_max`` (→ 0), ``div_rotation_max`` (→ 0), ``curl_rotation_grid``
    (→ 2b), plus the ``grid`` used.

    Validation: V1/V2 as listed; V6 factors 3 and 2 (private JSON). Label: symbolic, analytic.
    """
    g = grid((-L, L), n)
    ur, ub = radial_field(a), solid_body_rotation_field(b)
    Ur, Ub = evaluate_field(ur, g), evaluate_field(ub, g)
    div_r, curl_r = divergence(Ur, g.h), curl(Ur, g.h)
    div_b, curl_b = divergence(Ub, g.h), curl(Ub, g.h)
    return {"div_radial": ur.div_expr, "curl_radial": ur.curl_expr, "div_rotation": ub.div_expr,
            "curl_rotation": ub.curl_expr, "div_radial_grid": float(div_r.mean()),
            "curl_radial_max": float(np.max(np.abs(curl_r))), "div_rotation_max": float(np.max(np.abs(div_b))),
            "curl_rotation_grid": curl_b.reshape(3, -1).mean(axis=1), "grid": g}


# ----------------------------------------------------------------------------------------------------------------------
# §2.11 Example 2.4 — principal axes of a plane shear strain rate; linear-flow kinematics for E3
# ----------------------------------------------------------------------------------------------------------------------
def principal_angle_2d(S):
    """Angle of the first principal axis of a symmetric 2 × 2 tensor: ½ atan2(2 S₁₂, S₁₁ − S₂₂).

    Book: §2.11, Example 2.4 (45° for pure shear); D19's general angle. Returns radians in (−π/2, π/2].

    Validation: V1 45° for S = [[0, Γ], [Γ, 0]]; 0 for diagonal S; the eigenvector of ``principal_axes`` at that angle.
    Label: analytic.
    """
    s = _F(S)
    return float(0.5 * np.arctan2(2.0 * s[0, 1], s[0, 0] - s[1, 1]))


def example_2_4(Gamma=1.0) -> dict:
    """Example 2.4: diagonalise the plane strain rate S = [[0, Γ], [Γ, 0]] by rotating to its principal axes.

    Book: §2.11, Example 2.4 — det|S_ij − λδ_ij| = λ² − Γ² = 0 → λ¹ = Γ, λ² = −Γ; b¹ = (1, 1)/√2, b² = (−1, 1)/√2;
    C = [b¹ b²] = [[1/√2, −1/√2], [1/√2, 1/√2]] (a 45° rotation); S' = CᵀSC = diag(Γ, −Γ): stretching at rate Γ along
    b¹, compression along b², no shear in the principal frame.

    ⚠️ Γ here is the **off-diagonal element** S₁₂, as in the book's matrix. For the flow u = (u₁(x₂), 0) the
    definition S_ij = ½(∂u_i/∂x_j + ∂u_j/∂x_i) gives S₁₂ = ½ du₁/dx₂, so a velocity gradient du₁/dx₂ = γ̇ corresponds to
    Γ = γ̇/2 (the book's line "2S₁₂ = du₁/dx₂ = Γ" carries a stray factor 2 relative to its matrix; analysis §9.4).

    Parameters
    ----------
    Gamma : float — S₁₂ [1/s]

    Returns
    -------
    dict: ``S``, ``lam`` = (Γ, −Γ) in the book's order, ``B`` = [b¹ b²] (columns), ``C`` (= B, det +1), ``S_prime``,
    ``angle_rad`` (π/4), ``angle_deg`` (45), ``G`` = [[0, 2Γ], [0, 0]] (the velocity gradient with S as its symmetric part).

    Validation: V1 S B = B diag(λ); S' diagonal; det C = +1; ``principal_angle_2d(S)`` = 45°; V6 book values (private
    JSON). Label: analytic.
    """
    Gam = float(Gamma)
    S = np.array([[0.0, Gam], [Gam, 0.0]])  # Example 2.4: S = [[0, Γ], [Γ, 0]]
    lam, B = principal_axes(S)  # ascending: (−|Γ|, |Γ|)
    # book ordering: λ¹ = Γ with b¹ along (1, 1)/√2 — pick the column closest to (1, 1) and make b² = rotate(b¹, +90°)
    k = int(np.argmax(np.abs(B.T @ np.array([1.0, 1.0]))))
    b1 = B[:, k] * np.sign(B[0, k] + B[1, k])
    b2 = np.array([-b1[1], b1[0]])  # so that C = [b¹ b²] has det +1 (a rotation, not a reflection)
    C = np.column_stack([b1, b2])
    S_prime = transform_tensor(S, C)  # Eq. (2.12): S' = Cᵀ S C = diag(Γ, −Γ)
    return {"S": S, "lam": np.array([S_prime[0, 0], S_prime[1, 1]]), "B": C, "C": C, "S_prime": S_prime,
            "angle_rad": float(np.arctan2(C[1, 0], C[0, 0])), "angle_deg": float(np.rad2deg(np.arctan2(C[1, 0], C[0, 0]))),
            "G": np.array([[0.0, 2.0 * Gam], [0.0, 0.0]])}


VELOCITY_GRADIENT_PRESETS = ("simple_shear", "solid_body_rotation", "pure_strain", "uniaxial_extension", "irrotational_strain")


def velocity_gradient_preset(name: str, Gamma: float = 1.0, dim: int = 2) -> np.ndarray:
    """Velocity-gradient matrices G[i, j] = ∂u_i/∂x_j of the standard linear flows (E3 presets), rate scale Γ [1/s].

    Book: §2.10 (S + A split of ∂u_i/∂x_j), Examples 2.3 and 2.4; the names and matrices are the design's Part C 6.6
    / E3 preset list.

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

    Book: §2.10 (S + A decomposition seen as motion; E3), previewing Ch. 3 §3.4. Units: G [1/s], t [s].

    Validation: V1 pure rotation G = A gives a rotation matrix by angle ω₃ t; uniaxial strain gives diag(e^{Γt}, e^{−Γt});
    d/dt of the map at t = 0 equals G. Label: analytic.
    """
    return expm(_F(G) * float(t))  # x(t) = e^{Gt} x₀ solves dx/dt = G x


def deform_square(G, t, n_side: int = 10, half_width: float = 1.0, boundary_only: bool = False) -> np.ndarray:
    """Tracer positions of a square of material points carried by the linear flow u = G·x for time t.

    Returns an array ``(2, N)`` (or ``(3, N)`` for a 3 × 3 G — a square in the (1, 2) plane) of positions at time t;
    ``boundary_only`` samples the perimeter (4 n_side points, closed) instead of the filled n_side × n_side lattice.

    Book: E3 "strain vs rotation split" and Fig. 2.8's deforming square. Label: analytic.
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

    Book: §2.10 (S + A split), previewing Ch. 3 §3.4. Under solid-body rotation the line turns uniformly at the
    angular velocity of the fluid element, which is the vector of the antisymmetric part A: Ω = ½(∇×u)₃ = ½ω₃ (ω the
    book's vorticity, = the vector of ``rotation_tensor``); for the ``solid_body_rotation`` preset with G = [[0, −Γ],
    [Γ, 0]] that is Ω = Γ, so θ(t) = θ₀ + Γt. Under pure strain the line tends to the stretching axis; under simple
    shear it does both (E3 view 3). Units: rad (θ₀ may be an array).

    Validation: V1 solid-body rotation: θ(t) = θ₀ + Γ t (= θ₀ + ½(∇×u)₃ t); pure strain: θ → 0 as t → ∞.
    Label: analytic.
    """
    M = linear_flow_map(G, t)[:2, :2]
    th0 = _F(theta0)
    v = M @ np.stack([np.cos(th0), np.sin(th0)])
    return as_scalar_if_0d(np.arctan2(v[1], v[0]))


# ----------------------------------------------------------------------------------------------------------------------
# §2.12–2.13 convergence of the integral definitions (Examples 2.5, 2.6) — used by the notebook and scripts
# ----------------------------------------------------------------------------------------------------------------------
def integral_definition_convergence(kind: str = "divergence", x0=(0.3, -0.2, 0.5), hs: Sequence[float] = (0.4, 0.2, 0.1, 0.05),
                                    field: VectorField | None = None, n_face: int = 8) -> dict:
    """Error of the integral definitions (2.32) (``"divergence"``), (2.33) (``"curl"``), (2.31) on a scalar (``"gradient"``)
    or (2.35) (``"curl_component"``, n = e₃) against the exact sympy value, for a sequence of box/loop sizes h.

    Book: §2.12 Example 2.5 and §2.13 Example 2.6: the Cartesian formulas emerge in the limit h → 0; the error decays
    as h² (D22's neglected Taylor terms).

    Returns
    -------
    dict(h, err, values, exact, order) with ``order`` the least-squares slope of log err vs log h (≈ 2).

    Validation: V3 observed order 2.0 ± 0.15 for the smooth test field. Label: converged.
    """
    from tools.convergence import observed_order  # local import: tools/ is a test-time helper, kept out of the package graph

    f = smooth_test_field(3) if field is None else field
    x0 = _F(x0)
    if kind == "divergence":
        exact = float(f.div_fn(*x0))
        vals = [integral_divergence(f, x0, h, n_face) for h in hs]
    elif kind == "curl":
        exact = f.curl_fn(*x0)
        vals = [integral_curl(f, x0, h, n_face) for h in hs]
    elif kind == "gradient":
        phi = smooth_scalar_field(3)
        exact = phi.grad_fn(*x0)
        vals = [integral_gradient(phi, x0, h, n_face) for h in hs]
    elif kind == "curl_component":
        exact = float(f.curl_fn(*x0)[2])
        vals = [integral_curl_component(f, x0, [0.0, 0.0, 1.0], h, n_face) for h in hs]
    else:
        raise ValueError("kind must be divergence, curl, gradient or curl_component")
    err = np.array([float(np.max(np.abs(np.asarray(v) - np.asarray(exact)))) for v in vals])
    return {"h": np.asarray(hs, dtype=float), "err": err, "values": vals, "exact": exact,
            "order": observed_order(hs, err)}


def operator_convergence(op: str = "gradient", ns: Sequence[int] = (8, 16, 32, 64), bc: str = "onesided", L: float = 1.0) -> dict:
    """Observed order of the grid operators on the smooth test fields: max error vs spacing h for n = ns.

    ``op`` ∈ {"gradient", "divergence", "curl", "laplacian", "vector_gradient"}. Returns dict(n, h, err, order).

    Validation: V3 order 2.0 ± 0.1 including the one-sided edges (or periodic on [−π, π]). Label: converged.
    """
    from tools.convergence import observed_order

    errs, hs = [], []
    for n in ns:
        if bc == "periodic":
            g = grid((-np.pi, np.pi), n, periodic=True)
        else:
            g = grid((-L, L), n)
        if op in ("gradient", "laplacian"):
            phi = periodic_scalar_field() if bc == "periodic" else smooth_scalar_field(3)
            P = phi(*g.coords)
            if op == "gradient":
                num, ex = gradient(P, g.h, bc), phi.grad_fn(*g.coords)
            else:
                X = coordinates(3)
                lap = sp.lambdify(X, sum(sp.diff(phi.expr, x, 2) for x in X), "numpy")
                num, ex = laplacian(P, g.h, bc), lap(*g.coords)
        else:
            f = periodic_test_field() if bc == "periodic" else smooth_test_field(3)
            U = f(*g.coords)
            if op == "divergence":
                num, ex = divergence(U, g.h, bc), f.div_fn(*g.coords)
            elif op == "curl":
                num, ex = curl(U, g.h, bc), f.curl_fn(*g.coords)
            elif op == "vector_gradient":
                num, ex = vector_gradient(U, g.h, bc), f.grad_fn(*g.coords)
            else:
                raise ValueError("op must be gradient, divergence, curl, laplacian or vector_gradient")
        errs.append(float(np.max(np.abs(num - ex))))
        hs.append(g.h[0])
    return {"n": np.asarray(ns), "h": np.asarray(hs), "err": np.asarray(errs), "order": observed_order(hs, errs)}


_NOT_EXPORTED = {"np", "sp", "expm", "annotations", "Callable", "Sequence", "as_scalar_if_0d"}
__all__ = sorted(name for name in dir() if not name.startswith("_") and name not in _NOT_EXPORTED)
