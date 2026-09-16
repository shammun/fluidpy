"""Cartesian tensor algebra: bases, direction cosines, transformation rules, contractions, the stress tensor and
Cauchy's traction, δ and ε, symmetric/antisymmetric parts, principal axes.

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e, Ch. 2 §§2.1–2.11, Eqs. (2.1)–(2.29). Every equation was read
from the rendered page images (chapters/pages/ch02/p067–p083).

Conventions fixed here for the whole project (see ``knowledge/notation.md``)
------------------------------------------------------------------------------
* **Direction-cosine matrix** ``C[i, j] = e_i · e'_j`` (row = old axis, column = new axis): the *columns* of C are the
  new unit vectors written in old components. Components transform **passively**: ``x' = Cᵀ x`` (2.5), ``x = C x'``
  (2.7), ``τ' = Cᵀ τ C`` (2.12). Wikipedia's R(θ) and ``scipy.spatial.transform.Rotation.as_matrix()`` are *active*
  matrices (they rotate the vector); the book's C for a frame rotated by +θ *equals* R(θ) but is applied as Rᵀ.
* **Traction contracts the FIRST index**: ``f_i = τ_ji n_j`` (2.15); ``f = n·τ`` equals ``τ·n`` only for symmetric τ
  (proved in Ch. 4). The code never assumes symmetry here.
* **Double dot**: the book's ``A:B = A_ij B_ji`` (transpose pairing); the common Frobenius pairing ``A_ij B_ij`` is
  available with ``convention="frobenius"``. They coincide for symmetric operands.
* Stress sign: tensile positive; on the face with outward normal +e_i the positive components point along +e_j.
* Arrays are float64; every function also accepts plain lists (explainer parity rows call them with lists).

Units: this chapter is pure mathematics. Stress inputs are in Pa when physical; direction cosines, unit vectors and
δ, ε are dimensionless. Angles in radians.
"""
from __future__ import annotations

import itertools
from typing import Callable, Sequence

import numpy as np

from ._util import as_scalar_if_0d

__all__ = [
    "unit_vectors", "vector_from_components", "dot", "inner", "outer", "tensor_product", "contract",
    "dot_tensor_vector", "double_dot", "trace", "invariants", "characteristic_polynomial",
    "direction_cosines", "rotation_matrix_2d", "rotation_matrix_3d", "rotation_angle", "random_rotation",
    "orthogonality_residual", "is_orthogonal", "is_proper_rotation",
    "transform_vector", "inverse_transform_vector", "transform_tensor", "transforms_as_vector", "transforms_as_tensor",
    "STRESS_CUBE_FACES", "cube_face_tractions", "stress_component_meaning", "tetrahedron_face_areas",
    "traction", "traction_components", "normal_shear_stress", "mohr_circle_2d",
    "kronecker_delta", "levi_civita", "permutation_sign", "epsilon_delta_residual", "triple_product", "is_isotropic",
    "cross", "cross_einsum", "angle_between",
    "is_symmetric", "is_antisymmetric", "independent_components", "symmetric_part", "antisymmetric_part",
    "strain_rate_tensor", "rotation_tensor", "antisymmetric_from_vector", "vector_from_antisymmetric",
    "symmetric_double_contraction", "principal_axes", "diagonalize", "normal_stress_bounds",
]

_ARR = lambda a: np.asarray(a, dtype=float)  # noqa: E731  (short alias used everywhere below)


# ----------------------------------------------------------------------------------------------------------------------
# §2.1 basis vectors and the summation convention in numpy form
# ----------------------------------------------------------------------------------------------------------------------
def unit_vectors(n: int = 3) -> np.ndarray:
    """The Cartesian unit vectors e_1 … e_n as the rows of an identity matrix.

    Book: §2.1, Eq. (2.1) — x = e_1 x_1 + e_2 x_2 + e_3 x_3; e_1 = (1, 0, 0)ᵀ, e_2 = (0, 1, 0)ᵀ, e_3 = (0, 0, 1)ᵀ.

    Parameters
    ----------
    n : int
        Number of space dimensions (3 in the book; 2 for the plane examples).

    Returns
    -------
    E : ndarray, shape (n, n)
        ``E[i]`` is the unit vector e_i (dimensionless).

    Assumptions: right-handed orthonormal basis. Validation: V1 (E @ E.T = I; e_1 × e_2 = e_3). Label: analytic.
    """
    return np.eye(int(n))


def vector_from_components(components, E) -> np.ndarray:
    """Rebuild a vector from its components in a basis: x = Σ_i x_i e_i.

    Book: §2.1–2.2, Eq. (2.1) (basis e_i) and Eq. (2.3) (the same x in the primed basis e'_j: x = x'_j e'_j).

    Parameters
    ----------
    components : array_like, shape (n,)
        Components x_i (any unit) in the basis whose vectors are the rows of ``E``.
    E : array_like, shape (n, n)
        Basis vectors as rows, expressed in some reference frame (dimensionless).

    Returns
    -------
    x : ndarray, shape (n,)
        The vector in the reference frame (same unit as ``components``).

    Validation: V1 round trip with ``transform_vector`` (x from primed components equals x from unprimed ones).
    Label: analytic.
    """
    c = _ARR(components)
    return np.einsum("i,ij->j", c, _ARR(E))  # Eq. (2.1)/(2.3): x = Σ_i x_i e_i


def dot(a, b):
    """Dot product a·b = a_i b_i (implied sum over the repeated index).

    Book: §2.1, Eq. (2.2); restated §2.8.

    Parameters
    ----------
    a, b : array_like, shape (n,) or (..., n)
        Vector components (same unit each; the product carries the product unit).

    Returns
    -------
    float or ndarray
        a_i b_i, summed over the last axis.

    Validation: V1 against the hand-written loop and ``np.dot``; ``trace(outer(u, v)) == dot(u, v)`` (§2.8).
    Label: analytic.
    """
    return as_scalar_if_0d(np.einsum("...i,...i->...", _ARR(a), _ARR(b)))  # Eq. (2.2)


def inner(A, B) -> np.ndarray:
    """Matrix inner product P_ij = A_ik B_kj (summed over the adjacent index k), i.e. P = A·B.

    Book: §2.3, Eqs. (2.9)–(2.11).

    Parameters
    ----------
    A : array_like, shape (n, m)
    B : array_like, shape (m, p)

    Returns
    -------
    P : ndarray, shape (n, p)

    Validation: V1 equals ``A @ B`` and the hand-written double loop with an inner k-sum. Label: analytic.
    """
    return np.einsum("ik,kj->ij", _ARR(A), _ARR(B))  # Eq. (2.9)


def outer(u, v) -> np.ndarray:
    """Outer product P_ij = u_i v_j — a second-order tensor (Exercise 2.10; §2.4 text).

    Book: §2.4 ("the nine products u_i v_j … transform according to (2.12)").

    Parameters
    ----------
    u, v : array_like, shape (n,)

    Returns
    -------
    P : ndarray, shape (n, n)

    Validation: V1 ``P == np.outer(u, v)``; V7 transforms per (2.12) under random rotations. Label: analytic.
    """
    return np.einsum("i,j->ij", _ARR(u), _ARR(v))


def tensor_product(A, B) -> np.ndarray:
    """Tensor (outer) product of two arrays of any order: P_ij…kl… = A_ij… B_kl… (raises the order).

    Book: §2.5 ("If A and B are two second-order tensors, then the 81 numbers P_ijkl ≡ A_ij B_kl … form a
    fourth-order tensor").

    Returns
    -------
    ndarray, shape ``A.shape + B.shape``.

    Validation: V7 transforms according to Eq. (2.13) when A and B do. Label: analytic.
    """
    return np.multiply.outer(_ARR(A), _ARR(B))


def contract(A, B, pattern: str) -> np.ndarray:
    """A single (or double) contraction of A_ij B_kl written with the book's index letters.

    Book: §2.5, Eq. (2.14): ``"ij,ki->kj"`` gives A_ij B_ki = (B·A)_kj; ``"ij,ik->jk"`` gives (Aᵀ·B)_jk;
    ``"ij,kj->ik"`` gives (A·Bᵀ)_ik; ``"ij,jk->ik"`` gives (A·B)_ik.

    Parameters
    ----------
    A, B : array_like
        Second-order tensors (any consistent unit).
    pattern : str
        An ``np.einsum`` subscript string in the book's letters, e.g. ``"ij,ki->kj"``.

    Returns
    -------
    ndarray

    Validation: V1 the four patterns of (2.14) equal ``B @ A``, ``A.T @ B``, ``A @ B.T``, ``A @ B``. Label: analytic.
    """
    return np.einsum(pattern, _ARR(A), _ARR(B))  # Eq. (2.14)


def dot_tensor_vector(A, u, index: int = 1) -> np.ndarray:
    """Contracted product of a second-order tensor and a vector: A_ij u_j = (A·u)_i or A_ij u_i = (Aᵀ·u)_j.

    Book: §2.5 (text after Eq. (2.14)): the two possibilities.

    Parameters
    ----------
    A : array_like, shape (n, n)
    u : array_like, shape (n,)
    index : {1, 0}
        Which index of A is summed with u: ``1`` (the second, the book's A_ij u_j = A·u) or ``0`` (the first,
        A_ij u_i = Aᵀ·u). Cauchy's traction (2.15) is the ``index=0`` case: f_i = τ_ji n_j.

    Returns
    -------
    ndarray, shape (n,)

    Validation: V1 ``index=1 == A @ u``, ``index=0 == A.T @ u``; they differ for non-symmetric A. Label: analytic.
    """
    A_, u_ = _ARR(A), _ARR(u)
    if index == 1:
        return np.einsum("ij,j->i", A_, u_)  # A_ij u_j = (A·u)_i
    if index == 0:
        return np.einsum("ij,i->j", A_, u_)  # A_ij u_i = (Aᵀ·u)_j
    raise ValueError("index must be 0 (first index summed) or 1 (second index summed)")


def double_dot(A, B, convention: str = "book"):
    """Doubly contracted product of two second-order tensors — a scalar.

    Book: §2.5 (text after (2.14)): ``A:B = A_ij B_ji`` (**book convention**, transpose pairing) and
    ``A:Bᵀ = A_ij B_ij``. Many texts (and Wikipedia) call ``A_ij B_ij`` the double dot ("Frobenius"); the two agree
    when either operand is symmetric, which is the fluid-mechanics case (τ_ij ∂u_i/∂x_j in Ch. 4).

    Parameters
    ----------
    A, B : array_like, shape (n, n)
    convention : {"book", "frobenius"}
        ``"book"`` → A_ij B_ji (= trace(A·B)); ``"frobenius"`` → A_ij B_ij (= trace(A·Bᵀ)). State it where used.

    Returns
    -------
    float

    Validation: V1 ``book == trace(A @ B)``, ``frobenius == trace(A @ B.T)``; V7 both invariant under rotation of both
    operands. Label: analytic.
    """
    A_, B_ = _ARR(A), _ARR(B)
    if convention == "book":
        return as_scalar_if_0d(np.einsum("ij,ji->", A_, B_))  # A:B = A_ij B_ji  (book, §2.5)
    if convention == "frobenius":
        return as_scalar_if_0d(np.einsum("ij,ij->", A_, B_))  # A:Bᵀ = A_ij B_ij (Frobenius)
    raise ValueError('convention must be "book" (A_ij B_ji) or "frobenius" (A_ij B_ij)')


def trace(A):
    """Contraction A_jj = A_11 + A_22 + A_33 — the first invariant of a second-order tensor.

    Book: §2.5 (contraction; "A_jj is an invariant").

    Validation: V1 hand loop; V7 unchanged under random rotations (Exercise 2.9). Label: analytic.
    """
    return as_scalar_if_0d(np.einsum("...ii->...", _ARR(A)))  # A_jj (§2.5)


def invariants(A):
    """The three principal invariants of a second-order tensor (Exercise 2.9; used in §2.11).

    Book: §2.5 ("there are three independent invariants"); the standard set
    I₁ = A_ii, I₂ = ½(I₁² − A_ij A_ji), I₃ = det A, which are the coefficients of the characteristic polynomial
    det(A − λδ) = −λ³ + I₁λ² − I₂λ + I₃ (D18). For a 2 × 2 tensor I₂ = det A and there is no I₃.

    Parameters
    ----------
    A : array_like, shape (n, n), n = 2 or 3

    Returns
    -------
    (I1, I2, I3) : tuple of float
        For n = 2 the tuple is (I1, I2) only.

    Validation: V7 unchanged under 100 random rotations; V1 I₁ = Σλ, I₂ = Σ_{k<l} λ^k λ^l, I₃ = Πλ (Vieta);
    V5 Wikipedia "Cauchy stress tensor" invariants. Label: analytic.
    """
    A_ = _ARR(A)
    I1 = float(np.trace(A_))
    I2 = 0.5 * (I1 ** 2 - float(np.einsum("ij,ji->", A_, A_)))  # ½(I₁² − A_ij A_ji): a closed index chain (D18)
    if A_.shape[0] == 2:
        return I1, I2
    return I1, I2, float(np.linalg.det(A_))


def characteristic_polynomial(tau) -> np.ndarray:
    """Coefficients (highest power first, for ``np.roots``) of det(τ_ij − λδ_ij) = 0 written with the invariants.

    Book: §2.11 fact (1): det|τ_ij − λδ_ij| = 0; expanded (D18) to λ³ − I₁λ² + I₂λ − I₃ = 0 (3 × 3) or
    λ² − I₁λ + I₂ = 0 (2 × 2).

    Returns
    -------
    coeffs : ndarray
        ``[1, -I1, I2, -I3]`` (3 × 3) or ``[1, -I1, I2]`` (2 × 2).

    Validation: V1 ``np.roots`` of the coefficients equals ``np.linalg.eigvalsh`` (sorted) to 1e-10. Label: analytic.
    """
    inv = invariants(tau)
    if len(inv) == 2:
        return np.array([1.0, -inv[0], inv[1]])
    return np.array([1.0, -inv[0], inv[1], -inv[2]])  # λ³ − I₁λ² + I₂λ − I₃ (D18)


# ----------------------------------------------------------------------------------------------------------------------
# §2.2 rotation of axes: direction cosines and transformation rules
# ----------------------------------------------------------------------------------------------------------------------
def direction_cosines(E_old, E_new) -> np.ndarray:
    """Direction-cosine matrix C_ij = e_i · e'_j between two orthonormal bases.

    Book: §2.2, text after Eq. (2.5): "C_ij = e_i·e'_j is a 3 × 3 matrix of direction cosines". Row i = old axis,
    column j = new axis; the columns of C are the new unit vectors in old components.

    Parameters
    ----------
    E_old, E_new : array_like, shape (n, n)
        Basis vectors as **rows**, both written in the same reference frame (dimensionless unit vectors).

    Returns
    -------
    C : ndarray, shape (n, n)
        Passive rotation matrix: x' = Cᵀ x (2.5), x = C x' (2.7).

    Validation: V1 ``direction_cosines(I, E_new) == E_new.T``; orthogonality C Cᵀ = I (Exercise 2.8). Label: analytic.
    """
    return np.einsum("ik,jk->ij", _ARR(E_old), _ARR(E_new))  # C_ij = e_i · e'_j


def rotation_matrix_2d(theta) -> np.ndarray:
    """Direction-cosine matrix of a plane frame rotated by θ (counterclockwise) about the 3-axis.

    Book: §2.2, Example 2.1: ``C = [[cos θ, −sin θ], [sin θ, cos θ]]`` with e'_1 = e_r, e'_2 = e_θ.

    **Passive**: the columns are the new axes (cos θ, sin θ) and (−sin θ, cos θ); components transform as x' = Cᵀ x.
    This C equals the *active* R(θ) of Wikipedia/scipy (which rotates a vector by +θ) — the book applies it as Rᵀ.

    Parameters
    ----------
    theta : float
        Rotation angle of the axes [rad], positive counterclockwise.

    Returns
    -------
    C : ndarray, shape (2, 2)

    Validation: V1 entries against Wikipedia "Rotation matrix" R(θ); C(α)C(β) = C(α+β); CᵀC = I; V6 Example 2.1
    (u_r, u_θ) = Cᵀ(u_1, u_2). Label: analytic.
    """
    c, s = np.cos(float(theta)), np.sin(float(theta))
    return np.array([[c, -s], [s, c]])  # Example 2.1: C_ij = [[e1·e_r, e1·e_θ], [e2·e_r, e2·e_θ]]


def rotation_matrix_3d(axis, angle) -> np.ndarray:
    """Direction-cosine matrix of a frame rotated by ``angle`` about the unit vector ``axis`` (Rodrigues' formula).

    Book: §2.2 (C_ij = e_i·e'_j for a general rotation; the book gives only the 2-D example). Rodrigues:
    ``C = cos θ I + sin θ [k]_× + (1 − cos θ) k kᵀ`` — its columns are the rotated unit vectors e'_j, so C is the
    **passive** matrix of the book (x' = Cᵀ x), numerically equal to the active R of Wikipedia/scipy
    (``Rotation.from_rotvec(θ k).as_matrix()``).

    Parameters
    ----------
    axis : array_like, shape (3,)
        Rotation axis (normalised internally).
    angle : float
        Angle [rad], right-handed about ``axis``.

    Returns
    -------
    C : ndarray, shape (3, 3), det C = +1.

    Validation: V1 ``rotation_matrix_3d(e_3, θ)`` equals Wikipedia R_z(θ) and embeds ``rotation_matrix_2d``; C k = k
    (eigenvalue 1 along the axis); CᵀC = I; V5 agrees with scipy ``Rotation``. Label: analytic.
    """
    k = _ARR(axis)
    k = k / np.linalg.norm(k)
    th = float(angle)
    K = antisymmetric_from_vector(k)  # [k]_× so that K x = k × x
    return np.cos(th) * np.eye(3) + np.sin(th) * K + (1.0 - np.cos(th)) * np.outer(k, k)  # Rodrigues


def rotation_angle(C):
    """Rotation angle encoded by a direction-cosine matrix.

    2 × 2: ``atan2(C_21, C_11)`` (the angle of the new 1-axis, wrapped to (−π, π]); 3 × 3: ``arccos((tr C − 1)/2)``
    (the unsigned angle about the rotation axis).

    Book: §2.2 (inverse of ``rotation_matrix_2d``/``rotation_matrix_3d``; not a numbered equation).

    Validation: V1 ``rotation_angle(rotation_matrix_2d(θ)) == θ`` for θ ∈ (−π, π]; 3-D round trip. Label: analytic.
    """
    C_ = _ARR(C)
    if C_.shape == (2, 2):
        return float(np.arctan2(C_[1, 0], C_[0, 0]))
    return float(np.arccos(np.clip((np.trace(C_) - 1.0) / 2.0, -1.0, 1.0)))


def random_rotation(rng=None, n: int = 3) -> np.ndarray:
    """A uniformly random proper rotation (det +1) for invariance tests, from the QR factorisation of a Gaussian matrix.

    Book: tool for Exercises 2.8, 2.9, 2.11 (orthogonality, invariants, isotropy). Seed with
    ``np.random.default_rng(0)`` for reproducibility.

    Validation: V1 CᵀC = I and det C = +1 on 200 draws. Label: analytic.
    """
    rng = np.random.default_rng(0) if rng is None else rng
    Q, R = np.linalg.qr(rng.standard_normal((n, n)))
    Q = Q * np.sign(np.diag(R))  # make the factorisation unique (Haar-uniform)
    if np.linalg.det(Q) < 0:
        Q[:, -1] *= -1.0  # reflect one axis back: a proper rotation
    return Q


def orthogonality_residual(C):
    """max |CᵀC − δ| — how far a matrix is from satisfying C_ji C_jk = δ_ik.

    Book: §2.2 (used silently by Eq. (2.7)); Exercise 2.8. Zero for every direction-cosine matrix.

    Validation: V1 < 1e-15 for rotation matrices; > 0.1 for a scaled or random matrix. Label: analytic.
    """
    C_ = _ARR(C)
    return float(np.max(np.abs(C_.T @ C_ - np.eye(C_.shape[0]))))  # C_ji C_jk − δ_ik (D02)


def is_orthogonal(C, tol: float = 1e-12) -> bool:
    """True if C Cᵀ = Cᵀ C = δ to within ``tol`` (Exercise 2.8; D02)."""
    C_ = _ARR(C)
    return bool(orthogonality_residual(C_) < tol and np.max(np.abs(C_ @ C_.T - np.eye(C_.shape[0]))) < tol)


def is_proper_rotation(C, tol: float = 1e-12) -> bool:
    """True if C is orthogonal with det C = +1 (a rotation, not a reflection). Book: §2.2 / D02."""
    return bool(is_orthogonal(C, tol) and abs(np.linalg.det(_ARR(C)) - 1.0) < 1e3 * tol)


def transform_vector(u, C) -> np.ndarray:
    """Components of a vector in the rotated frame: u'_j = u_i C_ij (= Cᵀ u).

    Book: §2.2, Eqs. (2.5), (2.6), (2.8) — the formal definition of a vector.

    Parameters
    ----------
    u : array_like, shape (n,)
        Components in the original frame O123.
    C : array_like, shape (n, n)
        Direction cosines C_ij = e_i·e'_j (passive; columns = new axes).

    Returns
    -------
    u_prime : ndarray, shape (n,)
        Components in O1'2'3'.

    Validation: V1 equals ``C.T @ u`` and the hand loop; |u'| = |u|; round trip with ``inverse_transform_vector``;
    V6 Example 2.1. Label: analytic.
    """
    return np.einsum("i,ij->j", _ARR(u), _ARR(C))  # Eq. (2.5)/(2.8): u'_j = u_i C_ij


def inverse_transform_vector(u_prime, C) -> np.ndarray:
    """Components back in the original frame: u_j = u'_i C_ji (= C u').

    Book: §2.2, Eq. (2.7) (Exercise 2.2; D03). Note the summed index of C moves to the second slot.

    Validation: V1 equals ``C @ u_prime``; ``inverse_transform_vector(transform_vector(u, C), C) == u``. Label: analytic.
    """
    return np.einsum("i,ji->j", _ARR(u_prime), _ARR(C))  # Eq. (2.7): x_j = x'_i C_ji


def transform_tensor(T, C) -> np.ndarray:
    """Components of a Cartesian tensor of any order in the rotated frame: one C per index.

    Book: §2.4, Eq. (2.12) τ'_mn = C_im C_jn τ_ij (= Cᵀ τ C) and Eq. (2.13) A'_mnpq = C_im C_jn C_kp C_lq A_ijkl; order 1
    reduces to Eq. (2.8), order 0 is unchanged.

    Parameters
    ----------
    T : array_like, shape (n,)*order
        Components in the original frame.
    C : array_like, shape (n, n)
        Direction cosines (passive).

    Returns
    -------
    T_prime : ndarray, same shape as ``T``.

    Validation: V1 order 2 equals ``C.T @ T @ C``; order 4 equals a four-fold explicit loop; V2 sympy with symbolic
    C and τ; V7 the invariants of ``T`` are unchanged. Label: analytic.
    """
    T_, C_ = _ARR(T), _ARR(C)
    order = T_.ndim
    if order == 0:
        return T_
    if order > 8:
        raise ValueError("transform_tensor supports tensors up to order 8")
    tin = "abcdefgh"[:order]
    tout = "mnopqrst"[:order]
    spec = ",".join(f"{a}{b}" for a, b in zip(tin, tout)) + f",{tin}->{tout}"  # e.g. "am,bn,ab->mn"
    return np.einsum(spec, *([C_] * order), T_)  # Eq. (2.12)/(2.13): one C_ij per index


def _sample_points(points, rng, n: int = 12, dim: int = 3) -> np.ndarray:
    if points is not None:
        return np.atleast_2d(_ARR(points))
    rng = np.random.default_rng(0) if rng is None else rng
    return rng.uniform(-2.0, 2.0, size=(n, dim))


def transforms_as_vector(component_fn: Callable, C, points=None, vector_params: Sequence = (), rng=None) -> float:
    """Residual of the vector test (2.8): does a triple of functions of position transform like the position vector?

    Book: §2.2 ("a vector can be formally defined as any quantity whose components change similarly to the
    components of the position vector under rotation", Eq. (2.8)).

    For each sample point the rule is evaluated twice: ``u = component_fn(x, *params)`` in the old frame and
    ``u' = component_fn(x', *params')`` with x' = Cᵀ x and every vector parameter also rotated (p' = Cᵀ p). A genuine
    vector satisfies u' = Cᵀ u; the residual is max |u' − Cᵀ u| over the points.

    Parameters
    ----------
    component_fn : callable ``fn(x, *vector_params) -> array (n,)``
        The candidate rule, written once, in whatever frame its arguments are given.
    C : array_like, shape (n, n)
        Direction cosines (passive).
    points : array_like, shape (N, n), optional
        Sample points (old components); 12 random points in [−2, 2]³ (seed 0) if omitted.
    vector_params : sequence of array_like
        Constant vectors the rule depends on (e.g. b in b × x); they transform too. Scalars are closed over.

    Returns
    -------
    residual : float
        ≈ 0 (1e-14) for a vector (x, a x, b × x, ∇φ); O(1) for non-vectors ((x₁², x₂², x₃²), (|x|, 0, 0)).

    Validation: V1/V7 as listed. Label: analytic.
    """
    C_ = _ARR(C)
    P = _sample_points(points, rng, dim=C_.shape[0])
    params = [_ARR(p) for p in vector_params]
    params_p = [C_.T @ p for p in params]
    res = 0.0
    for x in P:
        u = _ARR(component_fn(x, *params))
        u_p = _ARR(component_fn(C_.T @ x, *params_p))
        res = max(res, float(np.max(np.abs(u_p - C_.T @ u))))  # Eq. (2.8): u'_j − u_i C_ij
    return res


def transforms_as_tensor(component_fn: Callable, C, points=None, vector_params: Sequence = (), rng=None) -> float:
    """Residual of the tensor test (2.12)/(2.13): does an array-valued rule transform with one C per index?

    Book: §2.4 ("The elements of a matrix represent the components of a second-order tensor only if they obey
    (2.12)"). Same protocol as :func:`transforms_as_vector`; the residual is max |T' − transform_tensor(T, C)|.

    Validation: V1 ≈ 0 for u_i v_j (Exercise 2.10), ∂u_i/∂x_j of a linear field, δ_ij; ≫ 0 for a fixed array declared
    "the same in every frame". Label: analytic.
    """
    C_ = _ARR(C)
    P = _sample_points(points, rng, dim=C_.shape[0])
    params = [_ARR(p) for p in vector_params]
    params_p = [C_.T @ p for p in params]
    res = 0.0
    for x in P:
        T = _ARR(component_fn(x, *params))
        T_p = _ARR(component_fn(C_.T @ x, *params_p))
        res = max(res, float(np.max(np.abs(T_p - transform_tensor(T, C_)))))  # Eq. (2.12)
    return res


# ----------------------------------------------------------------------------------------------------------------------
# §2.4 the stress tensor and its sign convention; §2.6 Cauchy's traction
# ----------------------------------------------------------------------------------------------------------------------
STRESS_CUBE_FACES: dict[str, dict] = {
    # Fig. 2.4 lettering as used by Example 2.5 (EADH/FBCG ⊥ x1, ABCD/EFGH ⊥ x2, ABFE/DCGH ⊥ x3)
    "+1": {"axis": 0, "sign": +1, "letters": "EADH"},
    "-1": {"axis": 0, "sign": -1, "letters": "FBCG"},
    "+2": {"axis": 1, "sign": +1, "letters": "ABCD"},
    "-2": {"axis": 1, "sign": -1, "letters": "EFGH"},
    "+3": {"axis": 2, "sign": +1, "letters": "ABFE"},
    "-3": {"axis": 2, "sign": -1, "letters": "DCGH"},
}


def cube_face_tractions(tau) -> dict[str, dict]:
    """The traction (force per unit area) on each of the six faces of the stress cube of Fig. 2.4.

    Book: §2.4 sign convention: "on a surface whose outward normal points in the positive direction of a coordinate
    axis, the normal and shear stresses are positive if they point in the positive directions of the other axes" —
    so the face with outward normal +e_i carries the traction (τ_i1, τ_i2, τ_i3) = row i of τ, and the opposite face
    (outward normal −e_i) carries −row i (equal and opposite: stresses *at a point*, Newton's third law).

    Parameters
    ----------
    tau : array_like, shape (3, 3)
        Stress tensor [Pa]; ``tau[i, j]``: first index = face normal, second = force direction.

    Returns
    -------
    faces : dict
        Keys ``"+1", "-1", "+2", "-2", "+3", "-3"``; each value ``{"normal": e, "traction": f, "letters": "ABCD"}``
        with ``f`` in Pa (the traction is also ``traction(tau, normal)`` — the same rule).

    Validation: V1 the six tractions sum to zero; face "+2" traction == τ[1, :] (τ21, τ22, τ23 as in the text);
    each face agrees with ``traction(tau, n)``. Label: analytic.
    """
    t = _ARR(tau)
    out = {}
    for key, meta in STRESS_CUBE_FACES.items():
        n = np.zeros(3)
        n[meta["axis"]] = meta["sign"]
        out[key] = {"normal": n, "traction": meta["sign"] * t[meta["axis"], :].copy(), "letters": meta["letters"]}
    return out


def stress_component_meaning(i: int, j: int) -> str:
    """Plain-words meaning of τ_ij (book indices 1..3): first index = face normal, second = force direction.

    Book: §2.4 ("the first (i) index of τ_ij denotes the direction of the surface normal, and the second (j) index
    denotes the force component direction"). Normal stress (i = j) positive when tensile.
    """
    if not (1 <= i <= 3 and 1 <= j <= 3):
        raise ValueError("book indices run 1..3")
    kind = "normal stress (positive = tensile)" if i == j else "shear stress"
    return (f"tau_{i}{j}: force per unit area along +x{j} on the face whose outward normal is +x{i} "
            f"(normal along x{i}, force along x{j}) — {kind}")


def tetrahedron_face_areas(n, dA):
    """Projected areas of the coordinate faces of the tetrahedron of Fig. 2.5: dA_i = n_i dA.

    Book: §2.6 ("The geometry of the tetrahedron requires: dA_i = n_i dA"). This is the vector-area identity for a
    closed surface, n dA − Σ_i e_i dA_i = 0.

    Parameters
    ----------
    n : array_like, shape (3,)
        Unit outward normal of the slanted face (normalised internally).
    dA : float
        Area of the slanted face [m²].

    Returns
    -------
    dA_i : ndarray, shape (3,)  [m²]

    Validation: V1 Σ_i dA_i e_i == n dA. Label: analytic.
    """
    n_ = _ARR(n)
    n_ = n_ / np.linalg.norm(n_)
    return n_ * float(dA)  # dA_i = n_i dA


def traction(tau, n) -> np.ndarray:
    """Cauchy's traction formula: force per unit area on a surface with outward unit normal n, f_i = τ_ji n_j (f = n·τ).

    Book: §2.6, Eq. (2.15). The **first** index of τ (the face-normal index) is contracted with n. ``f = τ·n``
    (second index) gives the same vector only when τ_ij = τ_ji, which Ch. 4 proves for the stress tensor; this
    function does not assume it.

    Parameters
    ----------
    tau : array_like, shape (n, n)
        Stress tensor [Pa] (2 × 2 or 3 × 3).
    n : array_like, shape (n,)
        Outward unit normal (dimensionless; not normalised here — pass a unit vector).

    Returns
    -------
    f : ndarray, shape (n,)  [Pa]

    Validation: V1 symmetric τ: f == τ @ n == n @ τ; non-symmetric τ: f == n @ τ ≠ τ @ n (discrimination); pure
    pressure τ = −pδ gives f = −p n; V2 sympy f_i − τ_ji n_j = 0; V6 Example 2.2; V5 Wikipedia "Cauchy stress
    tensor" T_j = σ_ij n_i. Label: analytic.
    """
    return np.einsum("ji,j->i", _ARR(tau), _ARR(n))  # Eq. (2.15): f_i = τ_ji n_j


def traction_components(traction, normal):
    """Split a traction (stress vector on a surface) into its normal stress and its shear-stress vector.

    Book: §1.3 (normal vs shear stresses; compression and tension; a fluid at rest carries no shear, §1.7);
    used by §2.6 Example 2.2 (normal stress √3a/2, shear a/2). Promoted from ``ch01_introduction`` in ch02.

    Parameters
    ----------
    traction : array_like, shape (n,) or (..., n)
        Force per unit area exerted on the face by the material outside it [Pa] (tuples and lists accepted).
    normal : array_like, shape (n,) or (..., n)
        Outward normal of the face; normalised internally [-].

    Returns
    -------
    (sigma_n, tau_vec, tau_mag) : tuple
        ``sigma_n`` = t·n̂ [Pa], signed: positive = tension (pulling outward), negative = compression;
        ``tau_vec`` = t − sigma_n n̂ [Pa], the tangential (shear) part, orthogonal to n̂;
        ``tau_mag`` = |tau_vec| [Pa].

    Validation: V1 on 50 random tractions and normals (seed 5): tau_vec · n̂ = 0 and sigma_n n̂ + tau_vec = t (1e-12);
    pure pressure on a face gives sigma_n = −p, tau_mag = 0. Label: analytic.
    """
    t = np.asarray(traction, dtype=float)
    n = np.asarray(normal, dtype=float)
    n = n / np.linalg.norm(n, axis=-1, keepdims=True)
    sigma_n = np.sum(t * n, axis=-1)  # dot product: normal component
    tau_vec = t - np.asarray(sigma_n)[..., None] * n  # remaining tangential part
    tau_mag = np.linalg.norm(tau_vec, axis=-1)
    return as_scalar_if_0d(sigma_n), tau_vec, as_scalar_if_0d(tau_mag)


def normal_shear_stress(tau, n):
    """Normal and shear stress on the plane with unit normal n: σ_n = f·n, τ_s = |f − σ_n n|, with f = n·τ (2.15).

    Book: §2.6 (Example 2.2: "The normal stress is therefore √3a/2, and the shear stress is a/2"); §2.11 fact (4)
    (σ_n lies between λ_min and λ_max).

    Parameters
    ----------
    tau : array_like, shape (n, n)  [Pa]
    n : array_like, shape (n,)
        Outward normal (normalised internally).

    Returns
    -------
    (sigma_n, tau_s, shear_dir) : (float [Pa], float [Pa], ndarray unit vector or zeros when τ_s = 0)

    Validation: V1 recomposition σ_n n + τ_s ŝ == f; V6 Example 2.2 at φ = 30°; V7 σ_n ∈ [λ_min, λ_max] over random n.
    Label: analytic.
    """
    n_ = _ARR(n)
    n_ = n_ / np.linalg.norm(n_)
    f = traction(tau, n_)
    sigma_n, tau_vec, tau_mag = traction_components(f, n_)
    shear_dir = tau_vec / tau_mag if tau_mag > 0 else np.zeros_like(tau_vec)
    return sigma_n, tau_mag, shear_dir


def mohr_circle_2d(tau):
    """Centre and radius of Mohr's circle for a symmetric 2 × 2 tensor: (σ_n, τ_s) over all planes lie on it.

    Book: §2.6 Example 2.2 (σ_n = a sin 2φ, τ_s = a cos 2φ trace a circle of radius |a| about 0) and §2.11 fact (4);
    the circle itself is Ch. 4 material previewed. centre = (τ₁₁ + τ₂₂)/2 = I₁/2, radius = √(((τ₁₁ − τ₂₂)/2)² + τ₁₂²)
    = (λ_max − λ_min)/2.

    Returns
    -------
    (center, radius) : floats [Pa]

    Validation: V1 centre ± radius == the eigenvalues (Wikipedia principal-stress formula). Label: analytic.
    """
    t = _ARR(tau)
    center = 0.5 * (t[0, 0] + t[1, 1])
    radius = float(np.hypot(0.5 * (t[0, 0] - t[1, 1]), t[0, 1]))
    return float(center), radius


# ----------------------------------------------------------------------------------------------------------------------
# §2.7 Kronecker delta and the alternating tensor; §2.8 dot and cross products
# ----------------------------------------------------------------------------------------------------------------------
def kronecker_delta(n: int = 3) -> np.ndarray:
    """δ_ij = 1 if i = j, 0 otherwise — the identity matrix. Book: §2.7, Eq. (2.16). Label: analytic."""
    return np.eye(int(n))  # Eq. (2.16)


def permutation_sign(i: int, j: int, k: int) -> int:
    """Value of ε_ijk for three index values: +1 cyclic (123, 231, 312), −1 anti-cyclic (321, 213, 132), 0 if any repeat.

    Book: §2.7, Eq. (2.18). Works with book indices 1..3 or numpy indices 0..2 (only the *order* matters).

    Validation: V1 against ``levi_civita()`` for all 27 triples. Label: analytic.
    """
    if i == j or j == k or i == k:
        return 0  # Eq. (2.18): any two indices equal
    # parity of the permutation that sorts (i, j, k): cyclic order ↔ even
    return 1 if (j - i) * (k - j) * (k - i) > 0 else -1


def levi_civita() -> np.ndarray:
    """The alternating tensor ε_ijk as a 3 × 3 × 3 array (numpy indices: ε[0,1,2] = +1).

    Book: §2.7, Eq. (2.18).

    Validation: V1 ε[0,1,2] = ε[1,2,0] = ε[2,0,1] = 1, ε[2,1,0] = ε[1,0,2] = ε[0,2,1] = −1, zeros elsewhere; moving an
    index two places keeps the value, one place flips the sign; V5 Wikipedia "Levi-Civita symbol". Label: analytic.
    """
    eps = np.zeros((3, 3, 3))
    for i, j, k in itertools.permutations(range(3)):
        eps[i, j, k] = permutation_sign(i, j, k)  # Eq. (2.18)
    return eps


def epsilon_delta_residual() -> float:
    """max over all 81 (i, j, l, m) of |Σ_k ε_ijk ε_klm − (δ_il δ_jm − δ_im δ_jl)| — the epsilon–delta relation check.

    Book: §2.7, Eq. (2.19) ("The reader can verify the validity of this relationship by choosing some values for
    the indices"). Returns exactly 0.0.

    Validation: V1 == 0; V2 sympy enumeration (D09). Label: analytic.
    """
    eps, d = levi_civita(), kronecker_delta()
    lhs = np.einsum("ijk,klm->ijlm", eps, eps)  # Σ_k ε_ijk ε_klm  (adjacent index summed)
    rhs = np.einsum("il,jm->ijlm", d, d) - np.einsum("im,jl->ijlm", d, d)  # δ_il δ_jm − δ_im δ_jl
    return float(np.max(np.abs(lhs - rhs)))  # Eq. (2.19)


def cross(u, v) -> np.ndarray:
    """Cross product written out in components (Exercise 2.14).

    Book: §2.8, Eq. (2.20): u × v = (u₂v₃ − u₃v₂) e₁ + (u₃v₁ − u₁v₃) e₂ + (u₁v₂ − u₂v₁) e₃.

    Validation: V1 == ``np.cross`` == ``cross_einsum`` on 1000 random pairs; u × v = −v × u; e₁ × e₂ = e₃;
    |u × v|² + (u·v)² = |u|²|v|² (Lagrange). Label: analytic.
    """
    u_, v_ = _ARR(u), _ARR(v)
    return np.array([u_[1] * v_[2] - u_[2] * v_[1],
                     u_[2] * v_[0] - u_[0] * v_[2],
                     u_[0] * v_[1] - u_[1] * v_[0]])  # Eq. (2.20)


def cross_einsum(u, v) -> np.ndarray:
    """Cross product in index form: (u × v)_k = ε_ijk u_i v_j (= ε_kij u_i v_j).

    Book: §2.8, Eq. (2.21).

    Validation: V1 equals Eq. (2.20) (``cross``) and ``np.cross``. Label: analytic.
    """
    return np.einsum("ijk,i,j->k", levi_civita(), _ARR(u), _ARR(v))  # Eq. (2.21)


def triple_product(a, b, c) -> np.ndarray:
    """Vector triple product a × (b × c) (Exercise 2.5), equal to (a·c) b − (a·b) c by the epsilon–delta relation.

    Book: §2.7–2.8 (Eq. (2.19) applied to Eq. (2.21); D09).

    Validation: V1 equals (a·c) b − (a·b) c on random vectors. Label: analytic.
    """
    return cross(a, cross(b, c))


def angle_between(u, v):
    """Angle θ between two vectors from u·v = uv cos θ, computed stably as atan2(|u × v|, u·v).

    Book: §2.8 ("u·v = uv cos θ", Exercises 2.12–2.13; "|u × v| = uv sin θ").

    Returns
    -------
    theta : float [rad] in [0, π].

    Validation: V1 equals arccos(u·v/(|u||v|)) away from 0 and π and is exact for parallel/antiparallel vectors.
    Label: analytic.
    """
    u_, v_ = _ARR(u), _ARR(v)
    if u_.shape[-1] == 2:
        s = abs(u_[0] * v_[1] - u_[1] * v_[0])
    else:
        s = np.linalg.norm(np.cross(u_, v_))
    return float(np.arctan2(s, np.dot(u_, v_)))


def is_isotropic(T, rng=None, n_rotations: int = 50, proper: bool = True, tol: float = 1e-12):
    """Test whether a tensor's components are unchanged by every rotation of the frame (isotropic tensor).

    Book: §2.7 ("δ_ij is an isotropic tensor … δ'_ij = δ_ij"; "there is also only one isotropic tensor of third
    order", ε_ijk — Exercise 2.11). ε is isotropic only under **proper** rotations (det C = +1); under a reflection it
    changes sign (a pseudotensor), which ``proper=False`` exposes.

    Parameters
    ----------
    T : array_like, any order
    rng : numpy Generator, optional (seed 0 if omitted)
    n_rotations : int
    proper : bool
        ``True``: random rotations only; ``False``: half of the trial matrices include a reflection.
    tol : float

    Returns
    -------
    (isotropic, max_residual) : (bool, float)

    Validation: V7 δ → True; ε → True (proper), False (improper, residual 2); a random tensor → False; λδ → True.
    Label: analytic.
    """
    T_ = _ARR(T)
    rng = np.random.default_rng(0) if rng is None else rng
    res = 0.0
    for k in range(int(n_rotations)):
        C = random_rotation(rng, T_.shape[0] if T_.ndim else 3)
        if not proper and k % 2 == 1:
            C = C @ np.diag([1.0] * (C.shape[0] - 1) + [-1.0])  # append a reflection: det C = −1
        res = max(res, float(np.max(np.abs(transform_tensor(T_, C) - T_))))
    return bool(res < tol), res


# ----------------------------------------------------------------------------------------------------------------------
# §2.10 symmetric and antisymmetric tensors
# ----------------------------------------------------------------------------------------------------------------------
def is_symmetric(B, tol: float = 1e-12) -> bool:
    """True if B_ij = B_ji within ``tol`` (six independent components). Book: §2.10."""
    B_ = _ARR(B)
    return bool(np.max(np.abs(B_ - B_.T)) <= tol * max(1.0, float(np.max(np.abs(B_)))))


def is_antisymmetric(B, tol: float = 1e-12) -> bool:
    """True if B_ij = −B_ji within ``tol`` (zero diagonal, three independent components). Book: §2.10."""
    B_ = _ARR(B)
    return bool(np.max(np.abs(B_ + B_.T)) <= tol * max(1.0, float(np.max(np.abs(B_)))))


def independent_components(B, tol: float = 1e-12) -> int:
    """Number of independent components: n(n+1)/2 if symmetric, n(n−1)/2 if antisymmetric, n² otherwise (§2.10)."""
    n = _ARR(B).shape[0]
    if is_symmetric(B, tol):
        return n * (n + 1) // 2
    if is_antisymmetric(B, tol):
        return n * (n - 1) // 2
    return n * n


def symmetric_part(B) -> np.ndarray:
    """S_ij = ½(B_ij + B_ji), the symmetric part of any second-order tensor (also for stacked fields B[i, j, ...]).

    Book: §2.10 (unnumbered): B_ij = ½(B_ij + B_ji) + ½(B_ij − B_ji) = S_ij + A_ij.

    Validation: V1 S + A == B, S symmetric; V7 S transforms per (2.12) when B does (D14). Label: analytic.
    """
    B_ = _ARR(B)
    return 0.5 * (B_ + np.swapaxes(B_, 0, 1))  # S_ij = ½(B_ij + B_ji)


def antisymmetric_part(B) -> np.ndarray:
    """A_ij = ½(B_ij − B_ji), the antisymmetric part (zero diagonal; three independent numbers = a vector, (2.27)).

    Book: §2.10 (unnumbered decomposition).

    For a velocity gradient G[i, j] = ∂u_i/∂x_j the antisymmetric part is **half** the book's rotation tensor:
    A = ½R (Ch. 3 (3.17) R_ij = ∂u_i/∂x_j − ∂u_j/∂x_i; see :func:`rotation_tensor`), and
    ``vector_from_antisymmetric(A)`` = ½∇×u = the angular velocity of the fluid element.

    Validation: V1 A antisymmetric, S + A == B; ``vector_from_antisymmetric(antisymmetric_part(∂u_i/∂x_j))`` = ½∇×u
    for u = b × x (= b). Label: analytic.
    """
    B_ = _ARR(B)
    return 0.5 * (B_ - np.swapaxes(B_, 0, 1))  # A_ij = ½(B_ij − B_ji)


def strain_rate_tensor(G) -> np.ndarray:
    """Strain-rate tensor S_ij = ½(∂u_i/∂x_j + ∂u_j/∂x_i) from the velocity gradient G[i, j] = ∂u_i/∂x_j.

    Book: §2.11, Example 2.4 (Ch. 3 material previewed). Units: 1/s.

    Validation: V1 == ``symmetric_part(G)``; Example 2.4's S for G = [[0, Γ], [0, 0]] gives S₁₂ = Γ/2 (see
    ``ch02.example_2_4`` for the book's Γ convention). Label: analytic.
    """
    return symmetric_part(G)  # S_ij = ½(∂u_i/∂x_j + ∂u_j/∂x_i)


def rotation_tensor(G) -> np.ndarray:
    """Rotation tensor R_ij = ∂u_i/∂x_j − ∂u_j/∂x_i = G − Gᵀ — **twice** the antisymmetric part of the velocity gradient.

    Book: §2.10, Eqs. (2.26)–(2.27) ("In Chapter 3, R is recognized as the rotation tensor corresponding to the
    vorticity vector ω"); Ch. 3 Eq. (3.15) R_ij = −ε_ijk ω_k with ω = ∇×u the vorticity, and Eq. (3.17)
    R_ij = ∂u_i/∂x_j − ∂u_j/∂x_i (no ½). The velocity gradient therefore splits as ∂u_i/∂x_j = S_ij + ½R_ij, i.e.
    R = 2A with A = ``antisymmetric_part(G)``; "ω and R represent twice the fluid element rotation rate".

    Parameters
    ----------
    G : array_like, shape (n, n[, ...]) — velocity gradient with ``G[i, j] = ∂u_i/∂x_j``  [1/s]

    Returns
    -------
    R : ndarray, same shape — antisymmetric  [1/s]

    Notes
    -----
    * ``vector_from_antisymmetric(rotation_tensor(G))`` = ∇×u, the vorticity (for u = b × x this is 2b).
    * ``vector_from_antisymmetric(antisymmetric_part(G))`` = ½∇×u, the angular velocity of the fluid element
      (for u = b × x this is b). Use the antisymmetric part A, not R, when you want the rotation *rate*.

    Validation: V1 R = 2·antisymmetric_part(G), R antisymmetric; ``vector_from_antisymmetric(R)`` for
    G = ∂(b × x)_i/∂x_j equals 2b = ∇×u (discriminates R from A); V2 sympy (3.17) entrywise. Label: analytic.
    """
    G_ = _ARR(G)
    return G_ - np.swapaxes(G_, 0, 1)  # Eq. (3.17): R_ij = ∂u_i/∂x_j − ∂u_j/∂x_i  (= 2 A_ij)


def antisymmetric_from_vector(omega) -> np.ndarray:
    """The antisymmetric tensor associated with a vector: R_ij = −ε_ijk ω_k, so that R·x = ω × x.

    Book: §2.10, Eqs. (2.26)–(2.27): R = [[0, −ω₃, ω₂], [ω₃, 0, −ω₁], [−ω₂, ω₁, 0]].
    A scalar ``omega`` is read as ω₃ (plane case) and returns the 2 × 2 block [[0, −ω₃], [ω₃, 0]].

    Parameters
    ----------
    omega : array_like, shape (3,) or scalar

    Returns
    -------
    R : ndarray, shape (3, 3) or (2, 2)

    Validation: V1 entrywise against (2.26) (V2 sympy); R @ x == ω × x (sign discrimination); round trip with
    ``vector_from_antisymmetric``. Label: analytic.
    """
    w = _ARR(omega)
    if w.ndim == 0:
        w3 = float(w)
        return np.array([[0.0, -w3], [w3, 0.0]])  # 2-D block of (2.26)
    return -np.einsum("ijk,k->ij", levi_civita(), w)  # Eq. (2.27): R_ij = −ε_ijk ω_k


def vector_from_antisymmetric(R):
    """The vector associated with an antisymmetric tensor: ω_k = −½ ε_ijk R_ij (inverse of (2.26)).

    Book: §2.10, Eq. (2.27) (the book's printed lower limit "i−1" is a misprint for i = 1). A 2 × 2 input returns
    the scalar ω₃ = R₂₁ (plane case).

    Validation: V1 round trip on random ω; ``vector_from_antisymmetric(antisymmetric_part(G))`` for G = ∂(b × x)_i/∂x_j
    equals b (= ½ ∇×u), while ``vector_from_antisymmetric(rotation_tensor(G))`` = 2b = ∇×u (the book's ω of (3.15)).
    Label: analytic.
    """
    R_ = _ARR(R)
    if R_.shape == (2, 2):
        return float(R_[1, 0])  # ω₃ = R₂₁ = −ε_213 ω₃ … the (2.26) block
    return -0.5 * np.einsum("ijk,ij->k", levi_civita(), R_)  # Eq. (2.27): ω_k = −½ ε_ijk R_ij


def symmetric_double_contraction(tau, B):
    """P = τ_kl B_kl for a symmetric τ split into its symmetric and antisymmetric parts: P = τ_ij S_ij + τ_ij A_ij.

    Book: §2.10, Eqs. (2.28)–(2.29): the A-term is its own negative, hence zero, so τ_ij B_ij = τ_ij S_ij
    = ½ τ_ij (B_ij + B_ji). The pairing here is the index-matched τ_kl B_kl written in (2.28) (Frobenius pairing;
    identical to the book's A:B for symmetric τ).

    Parameters
    ----------
    tau : array_like, shape (n, n) — symmetric (checked; a warning-free ValueError otherwise is *not* raised so the
        discrimination test can show P_A ≠ 0 for a non-symmetric τ).
    B : array_like, shape (n, n)

    Returns
    -------
    (P, P_S, P_A) : floats — total, symmetric-part and antisymmetric-part contributions.

    Validation: V1 P_A == 0 (1e-12) for symmetric τ and random B, P == P_S == double_dot(τ, B, "frobenius");
    non-symmetric τ gives P_A ≠ 0 (discrimination). Label: analytic.
    """
    t, B_ = _ARR(tau), _ARR(B)
    S, A = symmetric_part(B_), antisymmetric_part(B_)
    P_S = float(np.einsum("ij,ij->", t, S))  # τ_ij S_ij
    P_A = float(np.einsum("ij,ij->", t, A))  # τ_ij A_ij  (= 0 when τ is symmetric, Eq. (2.29))
    return P_S + P_A, P_S, P_A  # Eq. (2.28)


# ----------------------------------------------------------------------------------------------------------------------
# §2.11 eigenvalues and eigenvectors of a symmetric tensor
# ----------------------------------------------------------------------------------------------------------------------
def principal_axes(tau, tol: float = 1e-10):
    """Eigenvalues and orthonormal eigenvectors (principal axes) of a real symmetric tensor.

    Book: §2.11 facts (1)–(3): three real λ^k from det|τ_ij − λδ_ij| = 0; mutually orthogonal b^k from
    (τ_ij − λδ_ij) b_j = 0; in the frame of the eigenvectors τ' = diag(λ¹, λ², λ³). The returned B is a proper
    rotation (det B = +1) so it can be used directly as the direction-cosine matrix C of (2.12) (columns = new axes).

    Parameters
    ----------
    tau : array_like, shape (n, n), n = 2 or 3
        Real symmetric tensor [Pa or 1/s]. Raises ``ValueError`` if not symmetric (the facts fail otherwise).
    tol : float
        Symmetry tolerance (relative to max |τ|).

    Returns
    -------
    (lam, B) : (ndarray (n,), ndarray (n, n))
        Eigenvalues in **ascending** order; ``B[:, k]`` is the unit eigenvector of ``lam[k]``; det B = +1.
        For a repeated eigenvalue the eigenvectors are one valid orthonormal choice (compare with projectors).

    Validation: V1 τ B == B diag(λ); Bᵀ B = I, det B = +1; V2/V1 λ are the roots of ``characteristic_polynomial``;
    V5 2-D Mohr formula (σ_x+σ_y)/2 ± √(((σ_x−σ_y)/2)² + τ_xy²) (Wikipedia); V6 Example 2.4. Label: analytic.
    """
    t = _ARR(tau)
    if t.ndim != 2 or t.shape[0] != t.shape[1]:
        raise ValueError("tau must be a square matrix")
    if not is_symmetric(t, tol):
        raise ValueError("principal_axes needs a real symmetric tensor (§2.11 facts (1)–(4) fail otherwise)")
    lam, B = np.linalg.eigh(0.5 * (t + t.T))  # fact (1)–(2): real λ ascending, orthonormal columns
    if np.linalg.det(B) < 0:
        B[:, -1] *= -1.0  # a reflection → flip one axis so that C = B is a rotation (det +1)
    return lam, B


def diagonalize(tau):
    """Rotate a symmetric tensor to its principal axes: C = [b¹ b² b³] and τ' = Cᵀ τ C = diag(λ).

    Book: §2.11 fact (3) and Example 2.4 ("all the components of S in the rotated system can be found by carrying
    out the matrix product Cᵀ·S·C").

    Returns
    -------
    (C, tau_prime) : (ndarray (n, n), ndarray (n, n))
        ``C`` = direction cosines whose columns are the eigenvectors (ascending λ); ``tau_prime`` diagonal.

    Validation: V1 off-diagonal of τ' < 1e-12·|τ|, diagonal == λ; equals ``transform_tensor(tau, C)``. Label: analytic.
    """
    lam, B = principal_axes(tau)
    tau_prime = transform_tensor(tau, B)  # Eq. (2.12): τ'_mn = C_im C_jn τ_ij with C = B
    return B, tau_prime


def normal_stress_bounds(tau, rng=None, n: int = 2000):
    """Monte-Carlo check of §2.11 fact (4): over random unit normals, σ_n = n·τ·n stays within [λ_min, λ_max] and the
    shear |f − σ_n n| never exceeds (λ_max − λ_min)/2 (Mohr).

    Book: §2.11 fact (4) ("they cannot be larger than the largest λ or smaller than the smallest λ") — the sharp
    statement is about the normal stress on any plane (Rayleigh quotient), not the elements.

    Parameters
    ----------
    tau : array_like, shape (n, n) — symmetric [Pa]
    rng : numpy Generator, optional (seed 0)
    n : int — number of random directions

    Returns
    -------
    (sigma_min, sigma_max, tau_s_max) : floats [Pa] — the observed extremes over the sample.

    Validation: V7 σ_min ≥ λ_min − 1e-12, σ_max ≤ λ_max + 1e-12, τ_s_max ≤ (λ_max − λ_min)/2 + 1e-12. Label: analytic.
    """
    t = _ARR(tau)
    rng = np.random.default_rng(0) if rng is None else rng
    N = rng.standard_normal((int(n), t.shape[0]))
    N /= np.linalg.norm(N, axis=1, keepdims=True)
    F = np.einsum("ji,aj->ai", t, N)  # f_i = τ_ji n_j for every sample normal (2.15)
    sigma = np.einsum("ai,ai->a", F, N)
    shear = np.linalg.norm(F - sigma[:, None] * N, axis=1)
    return float(sigma.min()), float(sigma.max()), float(shear.max())
