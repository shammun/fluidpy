"""Vector algebra for the AOFD Chapter-0A primer and the rotating-frame chapters: dot and cross products, length,
angle, projection, parallelogram area, scalar triple product, and the horizontal rotation k̂ × u.

Vallis: no numbered equation (background maths); first used in §2.1, Eq. (2.7) v_I = v_R + Ω × r and Eq. (2.12)
(−2Ω × v_R), and in §2.3, Eq. (2.79) (f₀ × u = f₀ k × u).
Kundu: §2.1–§2.8, Eqs. (K 2.2) (a·b = a_i b_i) and (K 2.20)–(K 2.21) (cross product, components and ε_ijk form);
App. B §B.3 (the ∇ operator built on these products) — re-taught in Chapter 0A.

Conventions
-----------
* **Components on the last axis.** A vector is an array of shape (..., n) with n = 2 or 3; any leading axes broadcast
  (numpy rules), so a (N, 3) array of vectors dotted with one (3,) vector gives N numbers. This differs from
  ``fluidpy.core.tensors.cross`` (components on the first axis) — that function is left untouched.
* 2-component vectors are horizontal vectors (u, v) of the tangent plane (x east, y north); where a 3-D operation is
  needed (cross product, triple product) they are padded with a zero vertical component.
* Pure numbers in, pure numbers out: units are those of the inputs (a·b carries the product unit). Scalar results of
  single vectors come back as Python floats (scalar-callable for explainer ``selftest()`` parity rows).
"""
from __future__ import annotations

import numpy as np

from ..core._util import as_scalar_if_0d
from ..core.tensors import dot as _core_dot

__all__ = ["dot", "cross", "norm", "angle_between", "projection", "scalar_projection", "parallelogram_area",
           "triple_scalar", "k_cross"]

_S = as_scalar_if_0d


def _vec(a) -> np.ndarray:
    """float64 array with the components on the last axis; n must be 2 or 3."""
    v = np.asarray(a, dtype=float)
    if v.ndim == 0 or v.shape[-1] not in (2, 3):
        raise ValueError(f"a vector needs 2 or 3 components on the last axis; got shape {v.shape}")
    return v


def _pad3(v: np.ndarray) -> np.ndarray:
    """Pad a horizontal 2-vector (u, v) to (u, v, 0)."""
    if v.shape[-1] == 3:
        return v
    return np.concatenate([v, np.zeros(v.shape[:-1] + (1,))], axis=-1)


def _pair(a, b) -> tuple[np.ndarray, np.ndarray]:
    a, b = _vec(a), _vec(b)
    if a.shape[-1] != b.shape[-1]:
        a, b = _pad3(a), _pad3(b)
    return a, b


def dot(a, b):
    """Dot (scalar) product a·b = Σ_i a_i b_i = |a||b| cos θ.

    Vallis: background maths (used from §2.1, Eq. (2.7) onward; no numbered equation of its own).
    Kundu: §2.1, Eq. (K 2.2); App. B §B.3. Thin wrapper of ``fluidpy.core.tensors.dot`` (one implementation).

    Parameters
    ----------
    a, b : array_like, shape (..., n), n = 2 or 3 — vector components in any consistent unit (e.g. [m/s], [m]);
        leading axes broadcast. A 2-vector meeting a 3-vector is padded with a zero vertical component.

    Returns
    -------
    float or ndarray, shape (...) — a·b in [unit of a × unit of b] (e.g. [m² s⁻²] for two velocities).

    Assumptions: Cartesian (orthonormal) components.
    Validation: V1 against the hand-written sum and ``np.dot``; a·a = |a|²; e_x·e_y = 0; symmetric. Label: analytic.
    """
    a, b = _pair(a, b)
    return _S(_core_dot(a, b))  # Eq. (K 2.2): a_i b_i


def cross(a, b) -> np.ndarray:
    """Cross (vector) product a × b = (a₂b₃ − a₃b₂, a₃b₁ − a₁b₃, a₁b₂ − a₂b₁), perpendicular to a and b, |a × b| =
    |a||b| sin θ, right-handed.

    Vallis: background maths; first used in §2.1, Eq. (2.7) (Ω × r) and Eq. (2.12) (−2Ω × v_R).
    Kundu: §2.8, Eqs. (K 2.20)–(K 2.21); App. B §B.3.

    Parameters
    ----------
    a, b : array_like, shape (..., n), n = 2 or 3 — components [any unit]; 2-vectors are padded to (u, v, 0);
        leading axes broadcast.

    Returns
    -------
    ndarray, shape (..., 3) — a × b in [unit a × unit b] (e.g. Ω × r in [m/s]). Two horizontal vectors give a purely
    vertical result (0, 0, a₁b₂ − a₂b₁).

    Assumptions: right-handed Cartesian components (x east, y north, z up on the tangent plane).
    Validation: V1 equals ``np.cross`` on random pairs; a × b = −b × a; e_x × e_y = e_z; (a × b)·a = 0;
    Lagrange |a × b|² + (a·b)² = |a|²|b|². Label: analytic.
    """
    a, b = _pad3(_vec(a)), _pad3(_vec(b))
    a, b = np.broadcast_arrays(a, b)
    a1, a2, a3 = a[..., 0], a[..., 1], a[..., 2]
    b1, b2, b3 = b[..., 0], b[..., 1], b[..., 2]
    return np.stack([a2 * b3 - a3 * b2, a3 * b1 - a1 * b3, a1 * b2 - a2 * b1], axis=-1)  # Eq. (K 2.20)


def norm(a):
    """Length (Euclidean norm) |a| = √(a·a).

    Vallis: background maths (speeds |u| throughout, e.g. §2.9 inertial circles of radius U/f).
    Kundu: §2.1 (magnitude of a vector, from Eq. (K 2.2)).

    Parameters
    ----------
    a : array_like, shape (..., n), n = 2 or 3 — components [any unit].

    Returns
    -------
    float or ndarray, shape (...) — |a| [same unit as a].

    Validation: V1 |(3, 4)| = 5; |λa| = |λ||a|; equals ``np.linalg.norm``. Label: analytic.
    """
    a = _vec(a)
    return _S(np.sqrt(_core_dot(a, a)))


def angle_between(a, b, degrees: bool = False):
    """Angle θ between two vectors, from a·b = |a||b| cos θ and |a × b| = |a||b| sin θ.

    Vallis: background maths. Kundu: §2.1, Eq. (K 2.2) (a·b = |a||b| cos θ) with §2.8 (|a × b| = |a||b| sin θ).

    Parameters
    ----------
    a, b : array_like, shape (..., n), n = 2 or 3 — components [any unit; only directions matter].
    degrees : bool — return degrees instead of radians (interface convenience; default False).

    Returns
    -------
    float or ndarray — θ ∈ [0, π] [rad] (or [0, 180] [degrees]); NaN where either vector has zero length (no direction).

    Numerics: θ = atan2(|a × b|, a·b) — accurate for nearly parallel and nearly antiparallel vectors, where
    arccos(a·b/(|a||b|)) loses half the digits. (Our choice; mathematically identical.)
    Validation: V1 θ(e_x, e_y) = π/2, θ(a, a) = 0, θ(a, −a) = π; cos θ equals a·b/(|a||b|). Label: analytic.
    """
    a, b = _pair(a, b)
    axb = cross(a, b)
    s = np.sqrt(_core_dot(axb, axb))
    c = _core_dot(a, b)
    zero = (_core_dot(a, a) == 0) | (_core_dot(b, b) == 0)
    th = np.where(zero, np.nan, np.arctan2(s, c))  # θ = atan2(|a × b|, a·b)
    return _S(np.degrees(th) if degrees else th)


def scalar_projection(a, b):
    """Scalar projection (signed component) of a along b: a·b̂ = a·b/|b| = |a| cos θ.

    Vallis: background maths (e.g. the component of Ω along the local vertical, Ω sin ϑ, §2.3 Eq. (2.78)).
    Kundu: §2.1, Eq. (K 2.2).

    Parameters
    ----------
    a : array_like, shape (..., n) — the projected vector [any unit].
    b : array_like, shape (..., n) — the direction (its length does not matter) [any unit], nonzero.

    Returns
    -------
    float or ndarray — a·b̂ [unit of a]; negative when a points against b; NaN where |b| = 0.

    Validation: V1 projection of (3, 4) on e_x = 3; independent of |b|; equals |a| cos θ. Label: analytic.
    """
    a, b = _pair(a, b)
    nb = np.sqrt(_core_dot(b, b))
    with np.errstate(divide="ignore", invalid="ignore"):
        return _S(np.where(nb > 0, _core_dot(a, b) / np.where(nb > 0, nb, 1.0), np.nan))


def projection(a, b) -> np.ndarray:
    """Vector projection of a on b: (a·b/|b|²) b — the part of a that lies along b.

    Vallis: background maths. Kundu: §2.1, Eq. (K 2.2).

    Parameters
    ----------
    a : array_like, shape (..., n) — vector to project [any unit].
    b : array_like, shape (..., n) — direction [any unit], nonzero.

    Returns
    -------
    ndarray, shape (..., n) — the projection [unit of a]; NaN where |b| = 0. The remainder a − proj is ⟂ b.

    Validation: V1 (a − proj)·b = 0; proj of b on b is b; proj(a, λb) = proj(a, b) for λ ≠ 0. Label: analytic.
    """
    a, b = _pair(a, b)
    bb = _core_dot(b, b)
    with np.errstate(divide="ignore", invalid="ignore"):
        k = np.where(bb > 0, _core_dot(a, b) / np.where(bb > 0, bb, 1.0), np.nan)
    return np.asarray(k)[..., None] * b


def parallelogram_area(a, b):
    """Area of the parallelogram spanned by a and b, |a × b| = |a||b| sin θ.

    Vallis: background maths. Kundu: §2.8, Eq. (K 2.20) (magnitude of the cross product).

    Parameters
    ----------
    a, b : array_like, shape (..., n), n = 2 or 3 — edge vectors, e.g. in [m].

    Returns
    -------
    float or ndarray — area [unit a × unit b], e.g. [m²]; ≥ 0, zero for parallel vectors.

    Validation: V1 unit square = 1; |a||b| sin θ; shear invariance A(a, b + λa) = A(a, b). Label: analytic.
    """
    c = cross(a, b)
    return _S(np.sqrt(_core_dot(c, c)))


def triple_scalar(a, b, c):
    """Scalar triple product a·(b × c) — the signed volume of the parallelepiped on a, b, c (= det[a; b; c]).

    Vallis: background maths. Kundu: §2.8, Eq. (K 2.21) (ε_ijk a_i b_j c_k).

    Parameters
    ----------
    a, b, c : array_like, shape (..., n), n = 2 or 3 — edge vectors, e.g. in [m] (2-vectors padded with z = 0, so
        three horizontal vectors give 0).

    Returns
    -------
    float or ndarray — signed volume [product of the units], e.g. [m³]; > 0 for a right-handed triple.

    Validation: V1 equals ``np.linalg.det``; cyclic a·(b × c) = b·(c × a); e_x·(e_y × e_z) = 1; zero for
    coplanar vectors. Label: analytic.
    """
    a = _pad3(_vec(a))
    return _S(_core_dot(a, cross(b, c)))  # Eq. (K 2.21): ε_ijk a_i b_j c_k


def k_cross(u) -> np.ndarray:
    """k̂ × u for a horizontal vector u = (u, v): (−v, u) — u turned 90° to the left (anticlockwise seen from above).

    Vallis: §2.3, Eq. (2.79): the Coriolis term f₀ × u = f₀ k × u, so the Coriolis force per unit mass is −f k × u
    (Vallis §2.9, Eq. (2.233)), i.e. 90° to the right of u when f > 0. Kundu: §13.4, Eq. (K 13.7) (the (−f v, f u)
    components of 2Ω × u).

    Parameters
    ----------
    u : array_like, shape (..., 2) or (..., 3) — horizontal vector components (east, north[, up]) [any unit, e.g. m/s].
        For a 3-vector only the horizontal part is used and the result has a zero vertical component.

    Returns
    -------
    ndarray, same shape as u — k̂ × u [same unit]: (−v, u) or (−v, u, 0).

    Assumptions: k̂ is the local vertical (tangent-plane / traditional approximation).
    Validation: V1 k × e_x = e_y; k × (k × u) = −u_h; |k × u| = |u_h|; (k × u)·u = 0; equals ``cross((0,0,1), u)``.
    Label: analytic.
    """
    u = _vec(u)
    out = np.empty_like(u)
    out[..., 0] = -u[..., 1]  # Eq. (2.79): k × (u, v) = (−v, u)
    out[..., 1] = u[..., 0]
    if u.shape[-1] == 3:
        out[..., 2] = 0.0
    return out
