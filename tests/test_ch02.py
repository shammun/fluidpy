"""Verification suite for Chapter 2 — Cartesian Tensors (Kundu, Cohen & Dowling 5e, §§2.1–2.14).

Evidence levels (``verify-implementation`` skill): V1 analytic · V2 symbolic (sympy) · V3 convergence · V4 conservation /
invariant · V5 published benchmark (``reference/ch02/``, see SOURCES.md) · V6 book value (private,
``tests/book_values_ch02.json``, skipped when absent) · V7 limits / symmetry / invariance. Every test name is
``test_<concept>_<level>_<what>``; the comment on the ``def`` line repeats the level.

A items (CORE, ≥ 2 independent levels): C01 summation convention · C02 direction cosines · C03 vector transformation ·
C04 stress tensor · C05 Cauchy traction · C06 tensor rule · C07 contraction/invariants · C08 alternating tensor ·
C09 gradient · C10 divergence · C11 curl · C12 symmetric/antisymmetric split · C13 principal axes · C14 Gauss ·
C15 integral definitions · C16 Stokes.
Derivations re-derived with sympy (★★/★★★): D02 D05 D06 D09 D15 D17 D18 D21 D22 D25 D26 (+ D01, D12, D14 cheaply).

Pinned conventions with discrimination tests: passive C (x' = Cᵀx; ``rotation_matrix_3d`` == scipy's *active* matrix),
traction contracts the FIRST index, tensor divergence the SECOND, ``double_dot("book")`` = A_ij B_ji, G[i, j] = ∂u_i/∂x_j,
R·x = ω × x, Stokes loops counterclockwise about n (∮(x − c) × t ds = 2A n), Ex. 2.4 Γ ≡ S₁₂, Ex. 2.6 u_y, (2.27) i = 1.

Run: ``.venv/Scripts/python.exe -m pytest tests/test_ch02.py -q -p no:cacheprovider``.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import pytest
import sympy as sp
from scipy.spatial.transform import Rotation

from fluidpy import ch01_introduction as ch01
from fluidpy import ch02_cartesian_tensors as ch02
from fluidpy.core import tensors as T
from tools.convergence import observed_order, pairwise_orders

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference" / "ch02"
BOOK = Path(__file__).resolve().parent / "book_values_ch02.json"
needs_ref = pytest.mark.skipif(not (REF / "benchmarks.json").exists(), reason="run reference/ch02/make_refs.py")
book_only = pytest.mark.skipif(not BOOK.exists(), reason="book values are private; see CLAUDE.md rule 9")

RNG = np.random.default_rng(2)
ORDER_TOL = 0.15  # design order ± this (verify-implementation default)
EXAMPLE_TENSOR = np.array([[2.0, 1.0, 0.0], [1.0, 3.0, 1.0], [0.0, 1.0, 4.0]])  # design Part C 6.8: invariants 9, 24, 18


def book():
    return json.loads(BOOK.read_text(encoding="utf-8"))


def ref_json(name="benchmarks.json"):
    return json.loads((REF / name).read_text(encoding="utf-8"))


def rand_tensor(n=3, symmetric=False, rng=RNG):
    A = rng.standard_normal((n, n))
    return 0.5 * (A + A.T) if symmetric else A


def sym_rot_z(t):
    """Symbolic passive direction-cosine matrix of a frame rotated by t about e3 (columns = new axes)."""
    return sp.Matrix([[sp.cos(t), -sp.sin(t), 0], [sp.sin(t), sp.cos(t), 0], [0, 0, 1]])


def sym_rot_euler(a, b, c):
    Rz = lambda t: sp.Matrix([[sp.cos(t), -sp.sin(t), 0], [sp.sin(t), sp.cos(t), 0], [0, 0, 1]])  # noqa: E731
    Ry = lambda t: sp.Matrix([[sp.cos(t), 0, sp.sin(t)], [0, 1, 0], [-sp.sin(t), 0, sp.cos(t)]])  # noqa: E731
    return Rz(a) * Ry(b) * Rz(c)


# =====================================================================================================================
# C01 — the summation convention (index_notation)                                            N07 N12 N13 N21 N23 N40 N41
# =====================================================================================================================
def test_summation_convention_V1_dot_product_expansion():  # V1
    a1, a2, a3, b1, b2, b3 = sp.symbols("a_1 a_2 a_3 b_1 b_2 b_3", real=True)
    assert sp.expand(ch02.expand_indices("a_i b_i") - (a1 * b1 + a2 * b2 + a3 * b3)) == 0  # Eq. (2.2)
    assert ch02.expand_indices_str("a_i b_i") == "a_1*b_1 + a_2*b_2 + a_3*b_3"  # exact string (B1 parity)
    assert ch02.expand_indices("a_i b_i", values={"a": [1, 2, 3], "b": [4, 5, 6]}) == 32
    assert ch02.expand_indices("a_i b_i", dim=2, values={"a": [1, 2], "b": [4, 5]}) == 14
    assert float(ch02.dot([1, 2, 3], [4, 5, 6])) == 32.0
    assert isinstance(ch02.dot([1.0, 0.0], [0.0, 1.0]), float)  # scalar in → float out


def test_summation_convention_V1_free_dummy_and_errors():  # V1
    assert ch02.classify_indices("x_i C_ij") == (["j"], ["i"])  # §2.2 free j, dummy i
    assert ch02.classify_indices("A_ij B_kl") == (["i", "j", "k", "l"], [])
    assert ch02.tensor_order("A_ij B_kl") == 4 and ch02.tensor_order("a_i b_i") == 0 and ch02.tensor_order("tau_ji n_j") == 1
    assert ch02.rename_dummy("x_i C_ij", "i", "k") == "x_k C_kj"  # Eq. (2.6)
    assert ch02.expand_indices("x_i C_ij") == ch02.expand_indices("x_k C_kj")  # N13: same expansion
    with pytest.raises(ValueError):
        ch02.classify_indices("a_i b_i c_i")  # three repeats
    with pytest.raises(ValueError):
        ch02.rename_dummy("x_i C_ij", "j", "k")  # j is free
    with pytest.raises(ValueError):
        ch02.expand_indices("A_ij + B_ik")  # different free indices
    assert ch02.comma_to_partial("u_i,i") == "∂u_i/∂x_i"  # Eq. (2.36)
    assert "\\partial" in ch02.comma_to_partial("u_i,j", latex=True)


def test_kronecker_delta_V1_substitution_rule():  # V1  (N40, N41; Eqs. (2.16), (2.17))
    assert np.array_equal(ch02.kronecker_delta(), np.eye(3))
    u = ch02.symbol_vector("u")
    res = ch02.expand_indices("delta_ij u_j")
    assert all(sp.simplify(res[i] - u[i]) == 0 for i in range(3))  # Eq. (2.17)
    A = ch02.symbol_matrix("A")
    res2 = ch02.expand_indices("delta_ij A_jk")
    assert all(sp.simplify(res2[i, k] - A[i, k]) == 0 for i in range(3) for k in range(3))
    assert ch02.expand_indices("delta_ii") == 3  # summed: 3, not 1
    assert ch02.expand_indices("delta_11") == 1


def test_matrix_product_V1_inner_is_index_sum():  # V1  (N21, N22, N23; Eqs. (2.9)–(2.11))
    A, B = rand_tensor(), rand_tensor()
    assert np.allclose(ch02.inner(A, B), A @ B, atol=1e-14)
    P = ch02.expand_indices("A_ik B_kj")  # Eq. (2.9)
    As, Bs = ch02.symbol_matrix("A"), ch02.symbol_matrix("B")
    p12 = As[0, 0] * Bs[0, 1] + As[0, 1] * Bs[1, 1] + As[0, 2] * Bs[2, 1]  # Eq. (2.11): P_12
    assert sp.expand(P[0, 1] - p12) == 0
    hand = np.zeros((3, 3))
    for i, j, k in itertools.product(range(3), repeat=3):
        hand[i, j] += A[i, k] * B[k, j]
    assert np.allclose(hand, ch02.inner(A, B), atol=1e-14)


def test_summation_convention_V2_tensor_rule_cross_and_comma_forms():  # V2  (Eqs. (2.12), (2.21), (2.36))
    C, tau = ch02.symbol_matrix("C"), ch02.symbol_matrix("tau")
    res = ch02.expand_indices("C_im C_jn tau_ij")  # Eq. (2.12)
    mat = C.T * tau * C
    assert all(sp.expand(res[m, n] - mat[m, n]) == 0 for m in range(3) for n in range(3))
    u, v = ch02.symbol_vector("u"), ch02.symbol_vector("v")
    cr = ch02.expand_indices("eps_ijk u_i v_j")  # Eq. (2.21)
    assert all(sp.expand(cr[k] - u.cross(v)[k]) == 0 for k in range(3))
    k1 = ch02.expand_indices("eps_ij1 u_i v_j")  # N50: the k = 1 check — only (2,3), (3,2) survive
    assert sp.expand(k1 - (u[1] * v[2] - u[2] * v[1])) == 0
    X = ch02.coordinates(3)
    U = ch02.field_functions("u", 1)
    div = ch02.expand_indices("u_i,i")  # Eq. (2.36)
    assert sp.simplify(div - sum(sp.diff(U[i], X[i]) for i in range(3))) == 0
    curl = ch02.expand_indices("eps_ijk u_k,j")
    assert sp.simplify(curl[0] - (sp.diff(U[2], X[1]) - sp.diff(U[1], X[2]))) == 0  # Eq. (2.25) component 1
    assert sp.simplify(curl[2] - (sp.diff(U[1], X[0]) - sp.diff(U[0], X[1]))) == 0
    om = ch02.expand_indices("-1/2 eps_ijk R_ij")  # Eq. (2.27)
    R = ch02.symbol_matrix("R")
    assert sp.expand(om[2] - sp.Rational(-1, 2) * (R[0, 1] - R[1, 0])) == 0


# =====================================================================================================================
# C02 — direction cosines; N14 N15 (D02, D03)
# =====================================================================================================================
def test_direction_cosines_V1_columns_are_new_axes():  # V1
    th = 0.6
    C = ch02.rotation_matrix_2d(th)
    assert np.allclose(C, [[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]], atol=1e-15)  # Ex. 2.1's C
    E_new = C.T  # rows = new axes in old components
    assert np.allclose(ch02.direction_cosines(np.eye(2), E_new), E_new.T, atol=1e-15)  # C_ij = e_i·e'_j
    assert np.allclose(ch02.direction_cosines(np.eye(2), E_new), C, atol=1e-15)
    # 3-D: rotation about e3 embeds the 2-D matrix (Wikipedia R_z)
    C3 = ch02.rotation_matrix_3d([0, 0, 1], th)
    assert np.allclose(C3[:2, :2], C, atol=1e-15) and np.allclose(C3[2], [0, 0, 1], atol=1e-15)
    # composition and round trip of the angle
    a, b = 0.4, -1.1
    assert np.allclose(ch02.rotation_matrix_2d(a) @ ch02.rotation_matrix_2d(b), ch02.rotation_matrix_2d(a + b), atol=1e-14)
    for t in (-3.0, -1.0, 0.0, 0.3, 2.9):
        assert abs(ch02.rotation_angle(ch02.rotation_matrix_2d(t)) - t) < 1e-13
    assert abs(ch02.rotation_angle(ch02.rotation_matrix_3d([1, 1, 0], 1.3)) - 1.3) < 1e-13
    # the axis is the eigenvector with eigenvalue 1
    k = np.array([0.3, -0.5, 0.8])
    assert np.allclose(ch02.rotation_matrix_3d(k, 0.9) @ k, k, atol=1e-14)
    # random bases: C built from the bases reproduces the transformation of every vector
    Cr = ch02.random_rotation(np.random.default_rng(3))
    x = np.array([1.0, -2.0, 0.5])
    xp = ch02.transform_vector(x, Cr)
    assert np.allclose(ch02.vector_from_components(xp, Cr.T), x, atol=1e-14)  # x = Σ x'_j e'_j  (2.3) with e'_j = column j


def test_direction_cosines_V5_scipy_active_matrix():  # V5 (library cross-check, see SOURCES.md)
    rng = np.random.default_rng(11)
    for _ in range(20):
        k = rng.standard_normal(3)
        th = rng.uniform(-np.pi, np.pi)
        C = ch02.rotation_matrix_3d(k, th)
        R = Rotation.from_rotvec(th * k / np.linalg.norm(k)).as_matrix()  # active rotation by +θ about k
        assert np.allclose(C, R, atol=1e-13)  # the book's passive C *equals* R and is applied as Rᵀ
    # discrimination: the passive components of a fixed vector are Rᵀx, not Rx
    x = np.array([1.0, 0.0, 0.0])
    C = ch02.rotation_matrix_3d([0, 0, 1], np.pi / 2)  # axes turned +90° about e3
    assert np.allclose(ch02.transform_vector(x, C), [0.0, -1.0, 0.0], atol=1e-14)  # e1 seen from turned axes: −e'_2
    assert not np.allclose(ch02.transform_vector(x, C), C @ x)


def test_direction_cosines_V7_orthogonality_and_properness():  # V7  (N15, Exercise 2.8)
    rng = np.random.default_rng(0)
    for _ in range(200):
        C = ch02.random_rotation(rng)
        assert ch02.is_orthogonal(C) and ch02.is_proper_rotation(C)
        assert ch02.orthogonality_residual(C) < 1e-14
    assert not ch02.is_orthogonal(rand_tensor())
    assert not ch02.is_orthogonal(2.0 * ch02.rotation_matrix_2d(0.3))
    assert ch02.is_orthogonal(np.diag([1.0, 1.0, -1.0])) and not ch02.is_proper_rotation(np.diag([1.0, 1.0, -1.0]))
    assert ch02.orthogonality_residual(2.0 * np.eye(3)) > 0.1


def test_orthogonality_V2_derivation():  # V2  D02 (★★): C Cᵀ = Cᵀ C = I, det = +1, from C_ij = e_i·e'_j
    a, b, c = sp.symbols("alpha beta gamma", real=True)
    Ep = sym_rot_euler(a, b, c)  # columns = new orthonormal axes e'_j (any proper rotation)
    e_old = sp.eye(3)
    C = sp.Matrix(3, 3, lambda i, j: e_old[:, i].dot(Ep[:, j]))  # step 0: C_ij = e_i·e'_j
    # step 2–3: e_i expanded in the new basis, dotted with e_k, gives C_ij C_kj
    for i, k in itertools.product(range(3), repeat=2):
        e_i_expanded = sum((C[i, j] * Ep[:, j] for j in range(3)), sp.zeros(3, 1))
        assert sp.simplify(e_i_expanded.dot(e_old[:, k]) - sum(C[i, j] * C[k, j] for j in range(3))) == 0
    assert sp.simplify(C * C.T - sp.eye(3)) == sp.zeros(3)  # step 4
    assert sp.simplify(C.T * C - sp.eye(3)) == sp.zeros(3)  # step 6 (second completeness argument)
    assert sp.simplify(C.det() ** 2 - 1) == 0  # step 7: (det C)² = 1
    assert sp.simplify(C.det() - 1) == 0  # step 8: +1 for a rotation reached continuously from the identity
    assert C.subs({a: 0, b: 0, c: 0}) == sp.eye(3)
    # D03 (★): multiply (2.5) by C_kj, sum on j → x_k = x'_j C_kj
    x = ch02.symbol_vector("x")
    xp = C.T * x  # Eq. (2.5)
    assert sp.simplify(C * xp - x) == sp.zeros(3, 1)  # Eq. (2.7)


def test_inverse_transform_V1_round_trip():  # V1  (N14; Eq. (2.7))
    rng = np.random.default_rng(4)
    for _ in range(20):
        C = ch02.random_rotation(rng)
        x = rng.standard_normal(3)
        xp = ch02.transform_vector(x, C)
        assert np.allclose(ch02.inverse_transform_vector(xp, C), x, atol=1e-14)
        assert np.allclose(ch02.inverse_transform_vector(xp, C), C @ xp, atol=1e-14)
    x = np.array([1.0, 2.0, 3.0])
    assert np.allclose(ch02.vector_from_components(x, ch02.unit_vectors()), x)  # Eq. (2.1)


# =====================================================================================================================
# C03 — transformation of components and the formal definition of a vector (D01); N16 N18 N19
# =====================================================================================================================
def test_vector_transformation_V1_matrix_form_length_and_passive_sign():  # V1
    rng = np.random.default_rng(5)
    for _ in range(20):
        C = ch02.random_rotation(rng)
        x = rng.standard_normal(3)
        xp = ch02.transform_vector(x, C)
        assert np.allclose(xp, C.T @ x, atol=1e-14)  # N16: x' = Cᵀ x
        hand = np.array([sum(x[i] * C[i, j] for i in range(3)) for j in range(3)])  # Eq. (2.5) by hand
        assert np.allclose(xp, hand, atol=1e-14)
        assert abs(np.linalg.norm(xp) - np.linalg.norm(x)) < 1e-13
        assert np.allclose(ch02.transform_tensor(x, C), xp, atol=1e-14)  # order-1 tensor rule == (2.8)
    # passive: axes turned by +30°, a vector along the old 1-axis has a NEGATIVE 2'-component
    C = ch02.rotation_matrix_2d(np.pi / 6)
    xp = ch02.transform_vector([1.0, 0.0], C)
    assert xp[0] > 0 and xp[1] < 0 and np.allclose(xp, [np.cos(np.pi / 6), -np.sin(np.pi / 6)], atol=1e-15)


def test_vector_transformation_V2_derivation():  # V2  D01 (★): x'_j = x·e'_j = x_i C_ij from (2.1), (2.3)
    a, b, c = sp.symbols("alpha beta gamma", real=True)
    Ep = sym_rot_euler(a, b, c)
    x = ch02.symbol_vector("x")
    C = sp.Matrix(3, 3, lambda i, j: sp.eye(3)[:, i].dot(Ep[:, j]))
    for j in range(3):
        proj = x.dot(Ep[:, j])  # Eq. (2.4): x'_j = x·e'_j
        rule = sum(x[i] * C[i, j] for i in range(3))  # Eq. (2.5)
        assert sp.simplify(proj - rule) == 0
    # completeness check: x = Σ x'_j e'_j reconstructs x (Eq. (2.3))
    xp = [x.dot(Ep[:, j]) for j in range(3)]
    assert sp.simplify(sum((xp[j] * Ep[:, j] for j in range(3)), sp.zeros(3, 1)) - x) == sp.zeros(3, 1)


def test_vector_definition_V7_pass_fail_examples():  # V7  (N18; Eq. (2.8))
    C = ch02.random_rotation(np.random.default_rng(6))
    b = np.array([0.3, -1.2, 0.7])
    assert ch02.transforms_as_vector(lambda x: x, C) < 1e-13
    assert ch02.transforms_as_vector(lambda x: 2.5 * x, C) < 1e-13
    assert ch02.transforms_as_vector(lambda x, bb: np.cross(bb, x), C, vector_params=[b]) < 1e-13
    # ∇φ for a scalar written invariantly, φ = (b·x)² + |x|³ → ∇φ = 2(b·x) b + 3|x| x (the rule must be the same formula in both frames)
    grad_phi = lambda x, bb: 2 * (bb @ x) * bb + 3 * np.linalg.norm(x) * x  # noqa: E731
    assert ch02.transforms_as_vector(grad_phi, C, vector_params=[b]) < 1e-12
    # a frame-specific formula for a gradient is (correctly) rejected: ∇(x₁²x₂) written in old components
    assert ch02.transforms_as_vector(lambda x: np.array([2 * x[0] * x[1], x[0] ** 2, 0.0]), C) > 1e-2
    assert ch02.transforms_as_vector(lambda x: x ** 2, C) > 1e-2  # (x1², x2², x3²) is not a vector
    assert ch02.transforms_as_vector(lambda x: np.array([np.linalg.norm(x), 0.0, 0.0]), C) > 1e-2


def test_polar_components_V1_equals_rotation_of_axes():  # V1  (N19, Ex. 2.1 with our own numbers)
    u = np.array([1.0, 2.0])
    for th in np.linspace(-3.0, 3.0, 13):
        ur, ut = ch02.polar_components(u[0], u[1], th)
        assert np.allclose([ur, ut], ch02.transform_vector(u, ch02.rotation_matrix_2d(th)), atol=1e-14)
        assert abs(np.hypot(ur, ut) - np.linalg.norm(u)) < 1e-14
        assert np.allclose(ch02.cartesian_from_polar(ur, ut, th), u, atol=1e-14)
    assert ch02.polar_components(1.0, 2.0, 0.0) == (1.0, 2.0)
    assert np.allclose(ch02.polar_components(1.0, 2.0, np.pi / 2), (2.0, -1.0), atol=1e-15)
    assert np.allclose(ch02.polar_components_deg(1.0, 2.0, 90.0), (2.0, -1.0), atol=1e-15)
    ur, ut = ch02.polar_components(np.array([1.0, 0.0]), np.array([0.0, 1.0]), np.array([0.0, np.pi / 2]))
    assert np.allclose(ur, [1.0, 1.0]) and np.allclose(ut, [0.0, 0.0], atol=1e-15)  # radial vectors have u_θ = 0


@book_only
def test_polar_components_V6_book_example_2_1():  # V6
    ex = book()["example_2_1"]
    u1, u2, th = sp.symbols("u1 u2 theta", real=True)
    ur = sp.sympify(ex["u_r"]).subs({"u1": u1, "u2": u2, "theta": th})
    ut = sp.sympify(ex["u_theta"]).subs({"u1": u1, "u2": u2, "theta": th})
    f = sp.lambdify((u1, u2, th), [ur, ut], "numpy")
    for deg in ex["sample_theta_deg"]:
        t = np.deg2rad(deg)
        ours = ch02.polar_components(1.3, -0.7, t)
        assert np.allclose(ours, f(1.3, -0.7, t), rtol=1e-12, atol=1e-12)
    Cbook = sp.Matrix(2, 2, lambda i, j: sp.sympify(ex["C_matrix"][i][j]).subs("theta", th))
    Cours = sp.Matrix(ch02.rotation_matrix_2d(0.9))
    assert np.allclose(np.array(Cbook.subs(th, 0.9), dtype=float), np.array(Cours, dtype=float), atol=1e-14)


# =====================================================================================================================
# C04 — the stress tensor and its sign convention; N25 N26 R01
# =====================================================================================================================
def test_stress_cube_V1_face_tractions_follow_the_sign_convention():  # V1
    tau = rand_tensor()
    faces = ch02.cube_face_tractions(tau)
    assert set(faces) == {"+1", "-1", "+2", "-2", "+3", "-3"}
    assert np.allclose(faces["+2"]["traction"], tau[1, :])  # τ21, τ22, τ23 as in the text
    assert np.allclose(faces["-2"]["traction"], -tau[1, :])
    assert np.allclose(faces["+2"]["normal"], [0, 1, 0]) and np.allclose(faces["-3"]["normal"], [0, 0, -1])
    for f in faces.values():
        assert np.allclose(f["traction"], ch02.traction(tau, f["normal"]), atol=1e-14)  # the same rule as (2.15)
    assert {f["letters"] for f in faces.values()} == {"EADH", "FBCG", "ABCD", "EFGH", "ABFE", "DCGH"}  # Ex. 2.5 lettering
    s = ch02.stress_component_meaning(2, 3)
    assert "normal along x2" in s and "force along x3" in s and "shear" in s
    assert "normal stress" in ch02.stress_component_meaning(1, 1)
    with pytest.raises(ValueError):
        ch02.stress_component_meaning(0, 1)
    assert ch02.STRESS_CUBE_FACES["+1"]["letters"] == "EADH"


def test_stress_cube_V1_opposite_faces_cancel_at_a_point():  # V1  (R01: Newton III at a point)
    tau = rand_tensor()
    total = sum(f["traction"] for f in ch02.cube_face_tractions(tau).values())
    assert np.allclose(total, 0.0, atol=1e-14)
    for k in ("1", "2", "3"):
        faces = ch02.cube_face_tractions(tau)
        assert np.allclose(faces["+" + k]["traction"] + faces["-" + k]["traction"], 0.0, atol=1e-14)


def test_stress_cube_V7_face_tractions_are_vectors_under_rotation():  # V7
    rng = np.random.default_rng(7)
    tau = rand_tensor(rng=rng)
    for _ in range(10):
        C = ch02.random_rotation(rng)
        tau_p = ch02.transform_tensor(tau, C)
        for f in ch02.cube_face_tractions(tau).values():
            f_p = ch02.traction(tau_p, ch02.transform_vector(f["normal"], C))
            assert np.allclose(f_p, ch02.transform_vector(f["traction"], C), atol=1e-13)


# =====================================================================================================================
# C05 — Cauchy's traction formula (2.15) (D05); N34 N35 N38
# =====================================================================================================================
def test_cauchy_traction_V1_contracts_first_index_discrimination():  # V1
    rng = np.random.default_rng(8)
    n = rng.standard_normal(3)
    n /= np.linalg.norm(n)
    tau_s = rand_tensor(symmetric=True, rng=rng)
    f = ch02.traction(tau_s, n)
    assert np.allclose(f, tau_s @ n, atol=1e-14) and np.allclose(f, n @ tau_s, atol=1e-14)
    tau_a = rand_tensor(rng=rng)  # non-symmetric
    f = ch02.traction(tau_a, n)
    assert np.allclose(f, n @ tau_a, atol=1e-14) and np.allclose(f, tau_a.T @ n, atol=1e-14)  # f_i = τ_ji n_j
    assert not np.allclose(f, tau_a @ n)  # the wrong index would be τ_ij n_j
    assert np.max(np.abs(f - tau_a @ n)) > 0.1
    hand = np.array([sum(tau_a[j, i] * n[j] for j in range(3)) for i in range(3)])
    assert np.allclose(f, hand, atol=1e-14)
    assert np.allclose(ch02.dot_tensor_vector(tau_a, n, index=0), f, atol=1e-14)  # N32: A_ij u_i = (Aᵀ·u)_j
    # pure pressure: f = −p n for every n (ch01's isotropic pressure)
    p = 3.0
    assert np.allclose(ch02.traction(-p * np.eye(3), n), -p * n, atol=1e-14)
    # 2 × 2 and plain lists accepted (explainer parity rows)
    assert np.allclose(ch02.traction([[0, 1], [1, 0]], [0.8660254037844386, 0.5]), [0.5, 0.8660254037844386], atol=1e-12)
    # face check: n = e3 returns row 3 (D05 limit)
    assert np.allclose(ch02.traction(tau_a, [0, 0, 1]), tau_a[2], atol=1e-15)


def test_cauchy_traction_V2_derivation():  # V2  D05 (★★): the tetrahedron balance in the limit h → 0
    tau = ch02.symbol_matrix("tau")  # general (NOT assumed symmetric)
    n = ch02.symbol_vector("n")
    h, rho, k_A, k_V = sp.symbols("h rho k_A k_V", positive=True)
    g = ch02.symbol_vector("g")
    acc = ch02.symbol_vector("a")
    dA = k_A * h ** 2  # slanted face area ∝ h²
    dV = k_V * h ** 3  # volume ∝ h³
    f = ch02.symbol_vector("f")
    dA_j = [n[j] * dA for j in range(3)]  # step 6: dA_j = n_j dA (vector area of a closed surface)
    balance = [f[i] * dA - sum(tau[j, i] * dA_j[j] for j in range(3)) + rho * (g[i] - acc[i]) * dV for i in range(3)]  # step 5–7
    for i in range(3):
        sol = sp.solve(balance[i], f[i])[0]  # step 8: divide by dA
        assert sp.simplify(sol - (sum(tau[j, i] * n[j] for j in range(3)) - rho * (g[i] - acc[i]) * k_V * h / k_A)) == 0
        assert sp.simplify(sp.limit(sol, h, 0) - sum(tau[j, i] * n[j] for j in range(3))) == 0  # step 9: f_i = τ_ji n_j
    # index placement in sympy: f_i − τ_ji n_j == 0 with the expander, and it is NOT τ_ij n_j for general τ
    f_idx = ch02.expand_indices("tau_ji n_j")
    assert all(sp.expand(f_idx[i] - (tau.T * n)[i]) == 0 for i in range(3))
    assert any(sp.expand(f_idx[i] - (tau * n)[i]) != 0 for i in range(3))
    # N35: tetrahedron_face_areas closes the vector area
    nn = np.array([0.6, 0.0, 0.8])
    dAi = ch02.tetrahedron_face_areas(nn, 2.0)
    assert np.allclose(dAi, nn * 2.0) and np.allclose(sum(dAi[i] * np.eye(3)[i] for i in range(3)), nn * 2.0)


@needs_ref
def test_cauchy_traction_V5_wikipedia_index_placement():  # V5  (Cauchy stress tensor: T_j = σ_ij n_i)
    ref = ref_json()["cauchy_stress_tensor"]
    assert "sigma_ij n_i" in ref["traction"]
    tau = rand_tensor(rng=np.random.default_rng(9))
    n = np.array([0.36, 0.48, 0.8])
    assert np.allclose(ch02.traction(tau, n), np.einsum("ij,i->j", tau, n), atol=1e-14)  # same placement as the page


def test_cauchy_traction_V7_covariance_under_rotation():  # V7  (the D06 argument, numerically, for non-symmetric τ)
    rng = np.random.default_rng(10)
    for _ in range(20):
        tau, C = rand_tensor(rng=rng), ch02.random_rotation(rng)
        n = rng.standard_normal(3)
        n /= np.linalg.norm(n)
        f_old = ch02.traction(tau, n)
        f_new = ch02.traction(ch02.transform_tensor(tau, C), ch02.transform_vector(n, C))
        assert np.allclose(f_new, ch02.transform_vector(f_old, C), atol=1e-13)


def test_normal_shear_stress_V1_recomposition_and_ch01_reexport():  # V1
    rng = np.random.default_rng(12)
    for _ in range(20):
        tau = rand_tensor(symmetric=True, rng=rng)
        n = rng.standard_normal(3)
        sig, ts, sdir = ch02.normal_shear_stress(tau, n)
        nn = n / np.linalg.norm(n)
        assert np.allclose(sig * nn + ts * sdir, ch02.traction(tau, nn), atol=1e-13)
        assert abs(sdir @ nn) < 1e-13 or ts == 0.0
    assert ch01.traction_components is T.traction_components  # promoted in ch02; ch01 keeps working
    s, tv, tm = ch02.traction_components([0.5, 0.8660254037844386], [0.8660254037844386, 0.5])
    assert abs(s - 0.8660254037844386) < 1e-12 and abs(tm - 0.5) < 1e-12
    s, tm, sdir = ch02.normal_shear_stress(-2.0 * np.eye(3), [0.0, 0.0, 1.0])
    assert s == -2.0 and tm == 0.0 and np.allclose(sdir, 0.0)


def test_example_2_2_V1_general_forms_and_two_routes():  # V1  (N38: our own inputs)
    for a in (1.0, -2.5, 0.3):
        phis = np.linspace(0.0, 2 * np.pi, 50)
        for ph in phis:
            d = ch02.example_2_2(a, ph)
            assert abs(d["magnitude"] - abs(a)) < 1e-13  # |f| = |a|
            assert abs(d["sigma_n"] - a * np.sin(2 * ph)) < 1e-13  # τ'₁₁ = a sin 2φ
            assert abs(d["tau_s"] - a * np.cos(2 * ph)) < 1e-13  # τ'₁₂ = a cos 2φ
            e1p, e2p = d["C"][:, 0], d["C"][:, 1]
            assert abs(d["f"] @ e1p - d["tau_rot"][0, 0]) < 1e-13  # (2.15) route == (2.12) route
            assert abs(d["f"] @ e2p - d["tau_rot"][0, 1]) < 1e-13
            assert np.allclose(d["tau_rot"], d["C"].T @ d["tau"] @ d["C"], atol=1e-14)
            assert abs(np.trace(d["tau_rot"])) < 1e-13  # I₁ = 0 before and after
            assert 0.0 <= d["angle_rad"] < 2 * np.pi
        t = ch02.traction_2d(ch02.shear_flow_stress(a), 0.4)
        assert np.allclose(t["sigma_n"] * t["n"] + t["tau_s"] * t["s"], t["f"], atol=1e-14)
        assert abs(t["tau_s_mag"] - abs(t["tau_s"])) < 1e-15
    # a > 0 vs a < 0: the traction reverses (angle + 180°)
    dp, dm = ch02.example_2_2(1.0, np.pi / 6), ch02.example_2_2(-1.0, np.pi / 6)
    assert abs(((dm["angle_deg"] - dp["angle_deg"]) % 360.0) - 180.0) < 1e-10
    # stress_vs_angle: extremes of σ_n are the eigenvalues, τ_s vanishes there
    tau = np.array([[3.0, 1.0], [1.0, -1.0]])
    sv = ch02.stress_vs_angle(tau, np.linspace(0, np.pi, 20001))
    lam = np.linalg.eigvalsh(tau)
    assert abs(sv["sigma_n"].max() - lam[1]) < 1e-6 and abs(sv["sigma_n"].min() - lam[0]) < 1e-6
    assert abs(np.abs(sv["tau_s"]).max() - 0.5 * (lam[1] - lam[0])) < 1e-6
    assert abs(sv["tau_s"][np.argmax(sv["sigma_n"])]) < 1e-3


@book_only
def test_example_2_2_V6_book_values():  # V6
    ex = book()["example_2_2"]
    a = 1.0
    d = ch02.example_2_2(a, np.deg2rad(ex["phi_deg"]))
    assert np.allclose(d["tau"] / a, ex["tau_over_a"])
    assert np.allclose(d["n"], ex["n"], atol=1e-12)
    assert np.allclose(d["f"] / a, ex["f_over_a"], atol=1e-12)
    assert abs(d["magnitude"] / abs(a) - ex["magnitude_over_abs_a"]) < 1e-12
    assert abs(d["angle_deg"] - ex["angle_deg_a_positive"]) < 1e-10
    assert abs(ch02.example_2_2(-a, np.deg2rad(ex["phi_deg"]))["angle_deg"] - ex["angle_deg_a_negative"]) < 1e-10
    assert abs(d["tau_rot"][0, 0] / a - ex["tau11_rot_over_a"]) < 1e-12
    assert abs(d["tau_rot"][0, 1] / a - ex["tau12_rot_over_a"]) < 1e-12


@needs_ref
def test_mohr_circle_V5_principal_stress_formula_2d():  # V5  (Cauchy stress tensor page: σ₁,₂ and τ_max)
    ref = ref_json()["cauchy_stress_tensor"]
    assert "sqrt" in ref["principal_2d"]
    rng = np.random.default_rng(13)
    for _ in range(20):
        sx, sy, txy = rng.standard_normal(3)
        tau = np.array([[sx, txy], [txy, sy]])
        c, r = ch02.mohr_circle_2d(tau)
        s1 = 0.5 * (sx + sy) + np.sqrt((0.5 * (sx - sy)) ** 2 + txy ** 2)  # Wikipedia σ₁
        s2 = 0.5 * (sx + sy) - np.sqrt((0.5 * (sx - sy)) ** 2 + txy ** 2)
        assert abs(c + r - s1) < 1e-13 and abs(c - r - s2) < 1e-13
        lam, B = ch02.principal_axes(tau)
        assert np.allclose(lam, [s2, s1], atol=1e-13)
        assert abs(r - 0.5 * abs(s1 - s2)) < 1e-13  # τ_max = ½|σ₁ − σ₂|
        # the Mohr circle actually contains the (σ_n, τ_s) curve of stress_vs_angle
        sv = ch02.stress_vs_angle(tau, np.linspace(0, np.pi, 91))
        assert np.allclose(np.hypot(sv["sigma_n"] - c, sv["tau_s"]), r, atol=1e-12)


# =====================================================================================================================
# C06 — the tensor transformation rule (2.12), (2.13) (D06); N27 N28 N29
# =====================================================================================================================
def test_tensor_transformation_V1_matrix_form_and_higher_orders():  # V1
    rng = np.random.default_rng(14)
    C = ch02.random_rotation(rng)
    tau = rand_tensor(rng=rng)
    assert np.allclose(ch02.transform_tensor(tau, C), C.T @ tau @ C, atol=1e-14)  # Eq. (2.12)
    assert not np.allclose(ch02.transform_tensor(tau, C), C @ tau @ C.T)  # not the active form
    A4 = rng.standard_normal((3, 3, 3, 3))
    ref4 = np.zeros_like(A4)
    for m, n, p, q in itertools.product(range(3), repeat=4):
        for i, j, k, l in itertools.product(range(3), repeat=4):
            ref4[m, n, p, q] += C[i, m] * C[j, n] * C[k, p] * C[l, q] * A4[i, j, k, l]  # Eq. (2.13)
    assert np.allclose(ch02.transform_tensor(A4, C), ref4, atol=1e-12)
    assert ch02.transform_tensor(4.0, C) == 4.0  # order 0 unchanged
    # (2.13) structure: the tensor product of two tensors transforms as the product of their transforms (N28, N30)
    A, B = rand_tensor(rng=rng), rand_tensor(rng=rng)
    assert np.allclose(ch02.transform_tensor(ch02.tensor_product(A, B), C),
                       ch02.tensor_product(ch02.transform_tensor(A, C), ch02.transform_tensor(B, C)), atol=1e-13)
    assert ch02.tensor_product(A, B).shape == (3, 3, 3, 3)
    with pytest.raises(ValueError):
        ch02.transform_tensor(np.zeros((3,) * 9), C)


def test_tensor_transformation_V2_derivation():  # V2  D06 (★★): (2.12) from (2.15) + (2.8) + (2.7), step by step
    th = sp.symbols("theta", real=True)
    C = sym_rot_z(th)  # a passive C (orthogonal; det +1)
    tau = ch02.symbol_matrix("tau")  # no symmetry assumed
    n_p = ch02.symbol_vector("np")  # the new observer's normal n'
    n_old = C * n_p  # step 3: n_j = n'_m C_jm  (Eq. (2.7))
    f_old = tau.T * n_old  # step 2: f_i = τ_ji n_j (2.15)
    f_new_from_vector_rule = C.T * f_old  # step 1: f'_n = f_i C_in  (2.8)
    tau_p = C.T * tau * C  # step 8's matrix form
    f_new_from_cauchy = tau_p.T * n_p  # step 5: f'_n = τ'_mn n'_m
    diff = sp.simplify(f_new_from_vector_rule - f_new_from_cauchy)  # step 6
    assert diff == sp.zeros(3, 1)
    # step 7: "for every n'" ⇒ the coefficient matrix vanishes — compare coefficients of n'_m
    bracket = sp.simplify((f_new_from_vector_rule - f_new_from_cauchy).jacobian(n_p))
    assert bracket == sp.zeros(3)
    # step 8: the two index orders of the dummies are the same expansion (rename i ↔ j)
    r1 = ch02.expand_indices("C_im C_jn tau_ij", free_order="mn")
    r2 = ch02.expand_indices("C_jm C_in tau_ji", free_order="mn")
    assert all(sp.expand(r1[m, n] - r2[m, n]) == 0 for m in range(3) for n in range(3))
    # step 4 vs step 7 written with the expander: C_jm C_in τ_ji n'_m equals τ'_mn n'_m entrywise
    Cs = ch02.symbol_matrix("C")
    lhs = ch02.expand_indices("C_jm C_in tau_ji np_m")
    rhs = (Cs.T * tau * Cs).T * n_p
    assert all(sp.expand(lhs[k] - rhs[k]) == 0 for k in range(3))
    # limit C = I: τ' = τ
    assert tau_p.subs(th, 0) == tau


def test_tensor_definition_V7_pass_fail_examples():  # V7  (N27, N29; Exercise 2.10)
    C = ch02.random_rotation(np.random.default_rng(15))
    v = np.array([0.4, -1.0, 2.0])
    assert ch02.transforms_as_tensor(lambda x, vv: np.outer(vv, x), C, vector_params=[v]) < 1e-13  # u_i v_j
    G = rand_tensor(rng=np.random.default_rng(16))
    # velocity gradient of the linear field u = G x is G itself (a tensor): the rule must return CᵀGC in the new frame
    assert ch02.transforms_as_tensor(lambda x: G, C) > 1e-2  # a *fixed* array declared frame-independent fails
    assert ch02.transforms_as_tensor(lambda x: np.eye(3), C) < 1e-14  # δ passes (isotropic)
    assert ch02.transforms_as_tensor(lambda x: np.outer(x, x), C) < 1e-13  # x_i x_j passes


# =====================================================================================================================
# C07 — contraction, the invariants and the double dot (D18); N30 N31 N32 N33
# =====================================================================================================================
def test_contraction_V1_patterns_of_2_14_and_double_dot_convention():  # V1
    rng = np.random.default_rng(17)
    A, B, u = rand_tensor(rng=rng), rand_tensor(rng=rng), rng.standard_normal(3)
    assert np.allclose(ch02.contract(A, B, "ij,ki->kj"), B @ A, atol=1e-14)  # A_ij B_ki
    assert np.allclose(ch02.contract(A, B, "ij,ik->jk"), A.T @ B, atol=1e-14)  # A_ij B_ik
    assert np.allclose(ch02.contract(A, B, "ij,kj->ik"), A @ B.T, atol=1e-14)  # A_ij B_kj
    assert np.allclose(ch02.contract(A, B, "ij,jk->ik"), A @ B, atol=1e-14)  # A_ij B_jk
    assert np.allclose(ch02.dot_tensor_vector(A, u, 1), A @ u, atol=1e-14)
    assert np.allclose(ch02.dot_tensor_vector(A, u, 0), A.T @ u, atol=1e-14)
    assert not np.allclose(ch02.dot_tensor_vector(A, u, 1), ch02.dot_tensor_vector(A, u, 0))
    with pytest.raises(ValueError):
        ch02.dot_tensor_vector(A, u, 2)
    book_dd = ch02.double_dot(A, B, "book")
    frob = ch02.double_dot(A, B, "frobenius")
    assert abs(book_dd - np.trace(A @ B)) < 1e-13 and abs(frob - np.trace(A @ B.T)) < 1e-13  # N33
    assert abs(book_dd - frob) > 1e-3  # they differ for non-symmetric operands
    S = 0.5 * (A + A.T)
    assert abs(ch02.double_dot(S, B, "book") - ch02.double_dot(S, B, "frobenius")) < 1e-13  # agree when one is symmetric
    assert ch02.double_dot(A, B) == book_dd  # the default is the book's
    with pytest.raises(ValueError):
        ch02.double_dot(A, B, "other")
    assert abs(ch02.trace(ch02.outer(u, u)) - ch02.dot(u, u)) < 1e-13  # N45: u·v = trace of u_i v_j
    assert abs(ch02.trace(A) - (A[0, 0] + A[1, 1] + A[2, 2])) < 1e-15


def test_invariants_V1_vieta_and_characteristic_polynomial():  # V1
    I1, I2, I3 = ch02.invariants(EXAMPLE_TENSOR)
    assert abs(I1 - 9.0) < 1e-13 and abs(I2 - 24.0) < 1e-12 and abs(I3 - 18.0) < 1e-12
    lam = np.linalg.eigvalsh(EXAMPLE_TENSOR)
    assert abs(I1 - lam.sum()) < 1e-12
    assert abs(I2 - (lam[0] * lam[1] + lam[1] * lam[2] + lam[2] * lam[0])) < 1e-12
    assert abs(I3 - np.prod(lam)) < 1e-12
    coeffs = ch02.characteristic_polynomial(EXAMPLE_TENSOR)
    assert np.allclose(coeffs, [1.0, -9.0, 24.0, -18.0], atol=1e-12)
    assert np.allclose(np.sort(np.roots(coeffs)), lam, atol=1e-10)
    for _ in range(10):
        A = rand_tensor(rng=np.random.default_rng(18))
        assert np.allclose(np.sort(np.roots(ch02.characteristic_polynomial(A)).real), np.sort(np.linalg.eigvals(A).real), atol=1e-8)
    S2 = np.array([[3.0, 1.0], [1.0, -1.0]])
    I1, I2 = ch02.invariants(S2)
    assert abs(I1 - 2.0) < 1e-14 and abs(I2 - np.linalg.det(S2)) < 1e-14
    assert np.allclose(np.sort(np.roots(ch02.characteristic_polynomial(S2))), np.linalg.eigvalsh(S2), atol=1e-12)
    assert ch02.invariants(np.eye(3)) == (3.0, 3.0, 1.0)


def test_invariants_V7_unchanged_under_rotation_and_the_chain_remark():  # V7
    rng = np.random.default_rng(19)
    A = rand_tensor(rng=rng)
    inv = np.array(ch02.invariants(A))
    for _ in range(100):
        C = ch02.random_rotation(rng)
        Ap = ch02.transform_tensor(A, C)
        assert np.allclose(np.array(ch02.invariants(Ap)), inv, atol=1e-12)
        assert abs(ch02.trace(Ap) - ch02.trace(A)) < 1e-12
        B = rand_tensor(rng=rng)
        Bp = ch02.transform_tensor(B, C)
        for conv in ("book", "frobenius"):
            assert abs(ch02.double_dot(Ap, Bp, conv) - ch02.double_dot(A, B, conv)) < 1e-12
    # the coefficient of λ in the characteristic polynomial needs the closed chain A_ij A_ji, not A_ij A_ij
    I1 = np.trace(A)
    wrong_I2 = 0.5 * (I1 ** 2 - np.einsum("ij,ij->", A, A))
    assert abs(wrong_I2 - inv[1]) > 1e-3  # for non-symmetric A the two differ …
    assert not np.allclose(np.sort(np.roots([1.0, -I1, wrong_I2, -inv[2]])), np.sort(np.roots(ch02.characteristic_polynomial(A))))
    # … but note: A_ij A_ij = tr(A Aᵀ) is *also* rotation-invariant (the design's D18 remark says otherwise — see report)
    Ap = ch02.transform_tensor(A, ch02.random_rotation(rng))
    assert abs(np.einsum("ij,ij->", Ap, Ap) - np.einsum("ij,ij->", A, A)) < 1e-12


def test_invariants_V2_derivation():  # V2  D18 (★★): I₁, I₂, I₃ invariant; the cubic; Vieta in the principal frame
    th = sp.symbols("theta", real=True)
    C = sym_rot_z(th)
    A = ch02.symbol_matrix("A")
    Ap = C.T * A * C  # Eq. (2.12)
    assert sp.simplify(Ap.trace() - A.trace()) == 0  # steps 1–2: I₁
    assert sp.simplify((Ap * Ap).trace() - (A * A).trace()) == 0  # steps 3–4: the closed chain A_ij A_ji
    assert sp.simplify(Ap.det() - A.det()) == 0  # step 5: I₃
    # step 2 with the expander: C_im C_jm A_ij collapses to A_ii once C_im C_jm = δ_ij
    Cs = ch02.symbol_matrix("C")
    contracted = ch02.expand_indices("C_im C_jm A_ij")
    # the C-pair collapses to δ_ij for any orthogonal C: substitute a numeric rotation, keep A symbolic — every
    # coefficient of the remaining polynomial in the A_ij must vanish to round-off
    Cn = ch02.random_rotation(np.random.default_rng(20))
    subs = {Cs[i, j]: Cn[i, j] for i in range(3) for j in range(3)}
    residual = sp.expand(contracted.subs(subs) - A.trace())
    assert all(abs(float(coef)) < 1e-12 for coef in residual.as_coefficients_dict().values())
    # steps 6–7: det(A − λδ) = −λ³ + I₁λ² − I₂λ + I₃ with I₂ = ½(I₁² − A_ij A_ji) (symbolic 3 × 3, exact)
    lam = sp.symbols("lambda")
    I1 = A.trace()
    I2 = sp.Rational(1, 2) * (I1 ** 2 - (A * A).trace())
    I3 = A.det()
    poly = sp.expand((A - lam * sp.eye(3)).det())
    assert sp.expand(poly - (-lam ** 3 + I1 * lam ** 2 - I2 * lam + I3)) == 0
    M = (A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]) + (A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]) + (A[0, 0] * A[2, 2] - A[0, 2] * A[2, 0])
    assert sp.expand(M - I2) == 0  # step 7: M is the sum of the principal 2 × 2 minors
    # a wrong intermediate would be I₂ with A_ij A_ij: it is not the λ-coefficient for a general A
    I2_wrong = sp.Rational(1, 2) * (I1 ** 2 - sum(A[i, j] ** 2 for i in range(3) for j in range(3)))
    assert sp.expand(I2 - I2_wrong) != 0
    # step 8: principal frame
    l1, l2, l3 = sp.symbols("lambda1 lambda2 lambda3", real=True)
    D = sp.diag(l1, l2, l3)
    assert sp.expand(D.trace() - (l1 + l2 + l3)) == 0
    assert sp.expand(sp.Rational(1, 2) * (D.trace() ** 2 - (D * D).trace()) - (l1 * l2 + l2 * l3 + l3 * l1)) == 0
    assert sp.expand(D.det() - l1 * l2 * l3) == 0


# =====================================================================================================================
# C08 — the alternating tensor, the epsilon–delta relation, cross products (D09); N42 N43 N44 N46–N50
# =====================================================================================================================
def test_alternating_tensor_V1_values_and_index_moves():  # V1  (Eq. (2.18), N43)
    eps = ch02.levi_civita()
    assert eps[0, 1, 2] == eps[1, 2, 0] == eps[2, 0, 1] == 1
    assert eps[2, 1, 0] == eps[1, 0, 2] == eps[0, 2, 1] == -1
    assert np.count_nonzero(eps) == 6 and eps.shape == (3, 3, 3)
    for i, j, k in itertools.product(range(3), repeat=3):
        assert eps[i, j, k] == ch02.permutation_sign(i, j, k) == ch02.permutation_sign(i + 1, j + 1, k + 1)
        assert eps[i, j, k] == float(sp.LeviCivita(i + 1, j + 1, k + 1))
        if len({i, j, k}) < 3:
            assert eps[i, j, k] == 0
    assert np.array_equal(np.transpose(eps, (1, 2, 0)), eps)  # two places: ε_jki = ε_ijk
    assert np.array_equal(np.transpose(eps, (2, 0, 1)), eps)
    assert np.array_equal(np.transpose(eps, (1, 0, 2)), -eps)  # one place: ε_jik = −ε_ijk
    assert np.array_equal(np.transpose(eps, (0, 2, 1)), -eps)
    assert np.array_equal(np.transpose(eps, (2, 1, 0)), -eps)


def test_epsilon_delta_V1_identity_contractions_and_triple_product():  # V1  (N44; Eq. (2.19); Exercises 2.5, 2.7)
    eps, d = ch02.levi_civita(), ch02.kronecker_delta()
    assert ch02.epsilon_delta_residual() == 0.0
    assert np.array_equal(np.einsum("pqi,pqj->ij", eps, eps), 2 * d)  # 2δ_ij
    assert np.einsum("pqr,pqr->", eps, eps) == 6
    assert np.einsum("ii", d) == 3
    rng = np.random.default_rng(21)
    for _ in range(20):
        a, b, c = rng.standard_normal((3, 3))
        assert np.allclose(ch02.triple_product(a, b, c), (a @ c) * b - (a @ b) * c, atol=1e-13)  # BAC − CAB
    a, b, c = np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), np.array([1.0, 0, 0])
    assert np.allclose(ch02.triple_product(a, b, c), [0, 1.0, 0])


def test_epsilon_delta_V2_derivation():  # V2  D09 (★★): the 81 cases, the two contractions, the triple product
    # steps 1–7: full enumeration with sympy's LeviCivita / KroneckerDelta (exact integers)
    for i, j, l, m in itertools.product(range(1, 4), repeat=4):
        L = sum(sp.LeviCivita(i, j, k) * sp.LeviCivita(k, l, m) for k in range(1, 4))
        R = sp.KroneckerDelta(i, l) * sp.KroneckerDelta(j, m) - sp.KroneckerDelta(i, m) * sp.KroneckerDelta(j, l)
        assert L == R  # Eq. (2.19)
        if i == j or l == m:
            assert L == 0  # steps 2–3
        elif {l, m} != {i, j}:
            assert L == 0  # step 5
        elif (l, m) == (i, j):
            assert L == 1  # step 6
        else:
            assert L == -1  # step 7
    # via the expander (both free-index orders agree)
    lhs = ch02.expand_indices("eps_ijk eps_klm", free_order="ijlm")
    rhs = ch02.expand_indices("delta_il delta_jm - delta_im delta_jl", free_order="ijlm")
    assert lhs == rhs
    # step 8: ε_ijk ε_klj = 2δ_il  (δ_jj = 3 summed)
    c1 = ch02.expand_indices("eps_ijk eps_klj")
    assert c1 == sp.Array(2 * sp.eye(3))
    assert ch02.expand_indices("eps_pqi eps_pqj") == sp.Array(2 * sp.eye(3))  # the book's form (two-place moves)
    # step 9: fully contracted → 6
    assert ch02.expand_indices("eps_pqr eps_pqr") == 6
    # steps 10–11: a × (b × c) = (a·c) b − (a·b) c symbolically
    a, b, c = ch02.symbol_vector("a"), ch02.symbol_vector("b"), ch02.symbol_vector("c")
    lhs_v = ch02.expand_indices("eps_mpq a_p eps_qij b_i c_j")
    rhs_v = a.dot(c) * b - a.dot(b) * c
    assert all(sp.expand(lhs_v[k] - rhs_v[k]) == 0 for k in range(3))
    assert sp.simplify(a.cross(b.cross(c)) - rhs_v) == sp.zeros(3, 1)


@needs_ref
def test_levi_civita_V5_wikipedia_identities():  # V5  (Levi-Civita symbol page, three dimensions)
    ref = ref_json()["levi_civita_identities"]
    eps, d = ch02.levi_civita(), ch02.kronecker_delta()
    # (a) ε_ijk ε_imn = δ_jm δ_kn − δ_jn δ_km   (first index summed)
    lhs = np.einsum("ijk,imn->jkmn", eps, eps)
    rhs = np.einsum("jm,kn->jkmn", d, d) - np.einsum("jn,km->jkmn", d, d)
    assert np.array_equal(lhs, rhs)
    # (b) ε_jmn ε_imn = 2 δ_ij ; (c) ε_ijk ε_ijk = 6
    assert np.array_equal(np.einsum("jmn,imn->ij", eps, eps), ref["contraction_2delta_factor"] * d)
    assert np.einsum("ijk,ijk->", eps, eps) == ref["full_contraction"]
    # (d) (a × b)_i = ε_ijk a_j b_k ; (e) det A = ε_ijk a_1i a_2j a_3k
    rng = np.random.default_rng(22)
    a, b = rng.standard_normal(3), rng.standard_normal(3)
    assert np.allclose(np.einsum("ijk,j,k->i", eps, a, b), np.cross(a, b), atol=1e-14)
    assert np.allclose(np.einsum("ijk,j,k->i", eps, a, b), ch02.cross(a, b), atol=1e-14)
    A = rng.standard_normal((3, 3))
    assert abs(np.einsum("ijk,i,j,k->", eps, A[0], A[1], A[2]) - np.linalg.det(A)) < 1e-13
    # pseudotensor: a reflection flips its sign (the reason `is_isotropic(eps, proper=False)` is False)
    P = np.diag([1.0, 1.0, -1.0])
    assert np.array_equal(ch02.transform_tensor(eps, P), -eps)


def test_cross_product_V1_components_and_properties():  # V1  (N46, N47, N49; Eqs. (2.20)–(2.21))
    rng = np.random.default_rng(23)
    U, V = rng.standard_normal((1000, 3)), rng.standard_normal((1000, 3))
    for u, v in zip(U, V):
        w = ch02.cross(u, v)
        assert np.allclose(w, np.cross(u, v), atol=1e-13) and np.allclose(w, ch02.cross_einsum(u, v), atol=1e-13)
        assert np.allclose(ch02.cross(v, u), -w, atol=1e-15)
        assert abs(w @ u) < 1e-12 and abs(w @ v) < 1e-12
        assert abs(w @ w + (u @ v) ** 2 - (u @ u) * (v @ v)) < 1e-11  # Lagrange identity
        th = ch02.angle_between(u, v)
        assert abs(th - np.arccos(np.clip(u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1))) < 1e-12
        assert abs(np.linalg.norm(w) - np.linalg.norm(u) * np.linalg.norm(v) * np.sin(th)) < 1e-12
    e = ch02.unit_vectors()
    assert np.array_equal(ch02.cross(e[0], e[1]), e[2]) and np.array_equal(ch02.cross(e[1], e[2]), e[0])  # right-handed
    assert ch02.angle_between([1, 0, 0], [2, 0, 0]) == 0.0 and abs(ch02.angle_between([1, 0, 0], [-3, 0, 0]) - np.pi) < 1e-15
    assert abs(ch02.angle_between([1.0, 0.0], [0.0, 2.0]) - np.pi / 2) < 1e-15  # 2-D branch
    assert abs(ch02.angle_between([1, 1, 0], [1, 1, 1e-9]) - 1e-9 / np.sqrt(2)) < 1e-15  # stable near 0


def test_isotropic_tensors_V7_delta_epsilon_and_counterexamples():  # V7  (N42; Exercise 2.11)
    ok, res = ch02.is_isotropic(ch02.kronecker_delta())
    assert ok and res < 1e-14
    assert ch02.is_isotropic(2.5 * np.eye(3))[0]  # λδ
    ok, res = ch02.is_isotropic(ch02.levi_civita(), proper=True)
    assert ok and res < 1e-14
    ok, res = ch02.is_isotropic(ch02.levi_civita(), proper=False)
    assert not ok and abs(res - 2.0) < 1e-12  # sign flip under a reflection: max |−ε − ε| = 2
    assert not ch02.is_isotropic(rand_tensor(symmetric=True))[0]
    assert not ch02.is_isotropic(np.array([1.0, 0.0, 0.0]))[0]
    assert not ch02.is_isotropic(np.diag([1.0, 1.0, 2.0]))[0]


# =====================================================================================================================
# grids and operators — C09 gradient, C10 divergence, C11 curl (D12); N51 N52 N53 N54 N55 N74; R02
# =====================================================================================================================
def test_grid_V1_project_layout():  # V1  (math-to-python §3; core.grids)
    g = ch02.grid(((0, 1), (0, 2), (0, 3)), (5, 7, 9))
    assert g.shape == (9, 7, 5) and g.X.shape == (9, 7, 5)  # (nz, ny, nx)
    assert np.all(np.diff(g.X, axis=2) > 0) and np.all(np.diff(g.X, axis=0) == 0) and np.all(np.diff(g.X, axis=1) == 0)
    assert np.all(np.diff(g.Z, axis=0) > 0) and np.all(np.diff(g.Z, axis=2) == 0)
    assert np.allclose(g.h, (1 / 4, 2 / 6, 3 / 8))
    assert g.X[3, 2, 4] == g.x[4] and g.Y[3, 2, 4] == g.y[2] and g.Z[3, 2, 4] == g.z[3]
    gp = ch02.grid((-np.pi, np.pi), 8, periodic=True)
    assert abs(gp.h[0] - 2 * np.pi / 8) < 1e-15 and gp.x[-1] < np.pi  # duplicate end node dropped
    g2 = ch02.grid2d(((0, 1), (0, 2)), (5, 9))
    assert g2.shape == (9, 5) and g2.X[3, 4] == g2.x[4] and g2.Y[3, 4] == g2.y[3]  # u[j, i] = u(y_j, x_i)
    assert ch02.axis_of_direction(3, 0) == 2 and ch02.axis_of_direction(3, 2) == 0 and ch02.axis_of_direction(2, 0) == 1
    assert ch02.spacing((0.1, 0.2, 0.3), 1) == 0.2 and ch02.spacing(0.5, 2) == 0.5
    U = ch02.evaluate_field(lambda X, Y, Z: (X, 0.0, Z), g)
    assert U.shape == (3, 9, 7, 5) and np.all(U[1] == 0.0)
    with pytest.raises(ValueError):
        ch02.grid((0, 1), 2)


def test_partial_V1_exact_on_quadratics_including_edges():  # V1  (R02; the stencil behind every operator)
    g = ch02.grid2d((-1, 1), 17)
    phi = 2.0 * g.X ** 2 - 3.0 * g.X * g.Y + g.Y ** 2 + 0.5 * g.X - 1.0
    dx = ch02.partial(phi, 0, g.h)
    dy = ch02.partial(phi, 1, g.h)
    assert np.max(np.abs(dx - (4.0 * g.X - 3.0 * g.Y + 0.5))) < 1e-12  # edges included (one-sided is exact for quadratics)
    assert np.max(np.abs(dy - (-3.0 * g.X + 2.0 * g.Y))) < 1e-12
    grad = ch02.gradient(phi, g.h)
    assert grad.shape == (2, 17, 17) and np.allclose(grad[0], dx) and np.allclose(grad[1], dy)
    assert np.max(np.abs(ch02.second_partial(phi, 0, g.h) - 4.0)) < 1e-11
    assert np.max(np.abs(ch02.laplacian(phi, g.h) - 6.0)) < 1e-11
    g3 = ch02.grid((-1, 1), 9)
    r2 = g3.X ** 2 + g3.Y ** 2 + g3.Z ** 2
    assert np.max(np.abs(ch02.laplacian(r2, g3.h) - 6.0)) < 1e-11  # ∇²r² = 2d
    with pytest.raises(ValueError):
        ch02.partial(phi, 0, g.h, bc="neumann")
    # leading component axes are allowed
    U = np.stack([phi, -phi])
    assert np.allclose(ch02.partial(U, 0, g.h)[1], -dx)


@pytest.mark.parametrize("op,expected", [("gradient", 1.994), ("divergence", 1.988), ("curl", 2.052),
                                         ("laplacian", 2.026), ("vector_gradient", 1.994)])
def test_operators_V3_second_order_onesided(op, expected):  # V3  (C09, C10, C11; edges with one-sided stencils)
    r = ch02.operator_convergence(op, ns=(16, 32, 64, 128))
    pw = pairwise_orders(r["h"], r["err"])
    print(f"{op} onesided: order {r['order']:.3f} pairwise {[round(p, 3) for p in pw]} err {r['err']}")
    assert abs(r["order"] - 2.0) < ORDER_TOL
    assert abs(pw[-1] - 2.0) < 0.05  # asymptotic pair
    assert np.all(np.diff(r["err"]) < 0)


@pytest.mark.parametrize("op", ["gradient", "divergence", "curl", "laplacian"])
def test_operators_V3_second_order_periodic(op):  # V3  (bc="periodic" with periodic_test_field on [−π, π])
    r = ch02.operator_convergence(op, ns=(16, 32, 64, 128), bc="periodic")
    print(f"{op} periodic: order {r['order']:.3f} err {r['err']}")
    assert abs(r["order"] - 2.0) < ORDER_TOL


def test_gradient_V1_perpendicular_to_level_sets_and_directional_derivative():  # V1  (C09, N51)
    g = ch02.grid2d((-2, 2), 81)
    phi = g.X ** 2 + g.Y ** 2 / 4.0  # Fig. 2.7's ellipses
    grad = ch02.gradient(phi, g.h)
    exact = np.stack([2.0 * g.X, 0.5 * g.Y])
    assert np.max(np.abs(grad - exact)) < 1e-12  # quadratic: exact
    tangent = np.stack([-0.5 * g.Y, 2.0 * g.X])  # tangent to the level curves
    assert np.max(np.abs(np.einsum("i...,i...->...", grad, tangent))) < 1e-11  # ∇φ ⊥ level set
    mag = np.linalg.norm(grad, axis=0)
    rng = np.random.default_rng(24)
    for _ in range(10):
        n = rng.standard_normal(2)
        dn = ch02.directional_derivative(grad, n)
        assert np.all(np.abs(dn) <= mag + 1e-12)  # |∂φ/∂n| ≤ |∇φ|
    mask = mag > 0.1  # away from the origin, where ∇φ = 0 has no direction
    dn_par = ch02.directional_derivative(grad[:, mask], grad[:, mask])  # n ∥ ∇φ → |∇φ|
    assert np.allclose(dn_par, mag[mask], atol=1e-10)
    assert np.max(np.abs(ch02.directional_derivative(grad[:, mask], tangent[:, mask]))) < 1e-10  # n ⊥ ∇φ → 0
    sf = ch02.smooth_scalar_field(2)
    G = ch02.gradient(sf(*g.coords), g.h)
    assert np.max(np.abs(G - sf.grad_fn(*g.coords))) < 5e-3  # truncation only


def test_gradient_V7_rotational_invariance():  # V7  (∇φ is a vector: the rotated field's gradient is Cᵀ∇φ)
    th = 0.7
    C = ch02.rotation_matrix_2d(th)
    phi = ch02.smooth_scalar_field(2)
    errs, hs = [], []
    for n in (33, 65, 129):
        g = ch02.grid2d((-1, 1), n)
        # φ'(x') = φ(C x'): the same scalar field described in the rotated axes
        Xold = C[0, 0] * g.X + C[0, 1] * g.Y
        Yold = C[1, 0] * g.X + C[1, 1] * g.Y
        grad_p = ch02.gradient(phi(Xold, Yold), g.h)  # numerical gradient in the new frame
        exact_old = phi.grad_fn(Xold, Yold)  # exact ∇φ at the same physical points, old components
        exact_new = np.einsum("ij,i...->j...", C, exact_old)  # Cᵀ∇φ
        errs.append(np.max(np.abs(grad_p - exact_new)))
        hs.append(g.h[0])
    assert errs[-1] < 1e-3 and abs(observed_order(hs, errs) - 2.0) < 0.25


def test_divergence_V1_example_2_3_linear_fields_exact():  # V1  (C10, N55; Example 2.3 with our own a, b)
    a, b = 2.0, np.array([1.0, -2.0, 0.5])
    e = ch02.example_2_3(a, b, n=12, L=1.5)
    g = e["grid"]
    Ur = ch02.evaluate_field(ch02.radial_field(a), g)
    Ub = ch02.evaluate_field(ch02.solid_body_rotation_field(b), g)
    assert np.max(np.abs(ch02.divergence(Ur, g.h) - 3.0 * a)) < 1e-12  # ∇·(a x) = 3a everywhere (edges too)
    assert np.max(np.abs(ch02.curl(Ur, g.h))) < 1e-12
    assert np.max(np.abs(ch02.divergence(Ub, g.h))) < 1e-12
    cb = ch02.curl(Ub, g.h)
    for k in range(3):
        assert np.max(np.abs(cb[k] - 2.0 * b[k])) < 1e-12  # ∇×(b × x) = 2b
    assert abs(e["div_radial_grid"] - 3.0 * a) < 1e-12 and np.allclose(e["curl_rotation_grid"], 2.0 * b, atol=1e-12)
    assert e["curl_radial_max"] < 1e-12 and e["div_rotation_max"] < 1e-12
    assert ch02.is_solenoidal(Ub, g.h)[0] and ch02.is_irrotational(Ur, g.h)[0]  # N54
    assert not ch02.is_solenoidal(Ur, g.h)[0] and not ch02.is_irrotational(Ub, g.h)[0]
    g2 = ch02.grid2d((-1, 1), 9)
    U2 = ch02.evaluate_field(ch02.radial_field(a, dim=2), g2)
    assert np.max(np.abs(ch02.divergence(U2, g2.h) - 2.0 * a)) < 1e-12  # plane: 2a
    U2b = ch02.evaluate_field(ch02.solid_body_rotation_field(0.7, dim=2), g2)
    assert np.max(np.abs(ch02.curl(U2b, g2.h) - 1.4)) < 1e-12  # plane scalar curl 2b₃
    with pytest.raises(ValueError):
        ch02.divergence(Ur[:2], g.h)


def test_divergence_V2_example_2_3_symbolic():  # V2  (C10, C11; N55 — symbolic a, b)
    a = sp.Symbol("a", real=True)
    b1, b2, b3 = sp.symbols("b1 b2 b3", real=True)
    fr = ch02.radial_field(a)
    assert sp.simplify(fr.div_expr - 3 * a) == 0 and all(c == 0 for c in fr.curl_expr)
    fb = ch02.solid_body_rotation_field((b1, b2, b3))
    assert fb.div_expr == 0 and all(sp.simplify(fb.curl_expr[k] - 2 * bk) == 0 for k, bk in enumerate((b1, b2, b3)))
    X = ch02.coordinates(3)
    u = sp.Matrix([b1, b2, b3]).cross(sp.Matrix(X))
    assert list(u) == fb.exprs  # u = (b₂x₃ − b₃x₂, b₃x₁ − b₁x₃, b₁x₂ − b₂x₁) as printed (with a x₃ e₃ in the radial one)
    assert fr.exprs[2] == a * X[2]  # the book's "a x₂ e₂" typo written correctly
    d, c = ch02.exact_div_curl([X[0] ** 2, X[1] * X[2], sp.sin(X[0])])
    assert sp.simplify(d - (2 * X[0] + X[2])) == 0 and sp.simplify(c[1] - (-sp.cos(X[0]))) == 0
    d2, c2 = ch02.exact_div_curl(ch02.shear_field(3.0))
    assert d2 == 0 and c2 == -sp.Symbol("Gamma", real=True)  # (∇×u)₃ = −Γ for u₁ = Γx₂
    # the numeric twins agree with the symbolic ones on a grid
    g = ch02.grid((-1, 1), 6)
    f = ch02.smooth_test_field(3)
    assert np.allclose(f.div_fn(*g.coords), sp.lambdify(X, f.div_expr, "numpy")(*g.coords))
    assert f.grad_fn(*g.coords).shape == (3, 3, 6, 6, 6) and f(*g.coords).shape == (3, 6, 6, 6)


def test_vector_gradient_V1_index_order_and_tensor_divergence():  # V1  (N52: G[i, j] = ∂u_i/∂x_j; (∇·τ)_i = ∂τ_ij/∂x_j)
    g = ch02.grid((-1, 1), 9)
    U = np.stack([g.Y, np.zeros_like(g.X), np.zeros_like(g.X)])  # u = (x₂, 0, 0): ∂u₁/∂x₂ = 1
    G = ch02.vector_gradient(U, g.h)
    assert np.allclose(G[0, 1], 1.0, atol=1e-13) and np.allclose(G[1, 0], 0.0, atol=1e-13)  # pins the index order
    assert np.allclose(ch02.divergence(U, g.h), np.einsum("ii...->...", G), atol=1e-13)
    f = ch02.smooth_test_field(3)
    Gs = ch02.vector_gradient(f(*g.coords), g.h)
    assert np.max(np.abs(Gs - f.grad_fn(*g.coords))) < 0.05  # truncation only (h = 0.25)
    assert np.allclose(ch02.divergence(f(*g.coords), g.h), np.einsum("ii...->...", Gs), atol=1e-12)
    # tensor divergence: linear tensor field T_ij = a_ij·x → exact; index 1 vs 0 differ for non-symmetric T
    rng = np.random.default_rng(25)
    A = rng.standard_normal((3, 3, 3))  # T_ij = A_ijk x_k
    Tf = np.einsum("ijk,k...->ij...", A, np.stack([g.X, g.Y, g.Z]))
    div1 = ch02.tensor_divergence(Tf, g.h, index=1)
    div0 = ch02.tensor_divergence(Tf, g.h, index=0)
    ex1 = np.einsum("ijj->i", A)  # ∂T_ij/∂x_j = A_ijj
    ex0 = np.einsum("jij->i", A)  # ∂T_ji/∂x_j = A_jij
    for i in range(3):
        assert np.allclose(div1[i], ex1[i], atol=1e-12) and np.allclose(div0[i], ex0[i], atol=1e-12)
    assert not np.allclose(div1, div0)
    Ts = 0.5 * (Tf + np.swapaxes(Tf, 0, 1))
    assert np.allclose(ch02.tensor_divergence(Ts, g.h, index=1), ch02.tensor_divergence(Ts, g.h, index=0), atol=1e-12)
    # outer-product field vs sympy
    X = ch02.coordinates(3)
    u = [X[0] * X[1], X[2] ** 2, X[0] * X[2]]
    v = [X[1], X[0] + X[2], sp.Integer(1)]
    Tsym = sp.Matrix(3, 3, lambda i, j: u[i] * v[j])
    ex = sp.lambdify(X, [sum(sp.diff(Tsym[i, j], X[j]) for j in range(3)) for i in range(3)], "numpy")
    Tn = np.array(sp.lambdify(X, Tsym.tolist(), "numpy")(*g.coords), dtype=object)
    Tnum = np.stack([np.stack([np.broadcast_to(np.asarray(Tn[i, j], float), g.shape) for j in range(3)]) for i in range(3)])
    exv = np.stack([np.broadcast_to(np.asarray(e, float), g.shape) for e in ex(*g.coords)])
    assert np.max(np.abs(ch02.tensor_divergence(Tnum, g.h) - exv)) < 1e-11  # quadratic → exact


def test_curl_V1_components_equal_index_form():  # V1  (C11, N53; D12; Eqs. (2.24)–(2.25))
    g = ch02.grid((-1, 1), 11)
    f = ch02.smooth_test_field(3)
    U = f(*g.coords)
    assert np.max(np.abs(ch02.curl(U, g.h) - ch02.curl_components(U, g.h))) < 1e-14
    g2 = ch02.grid2d((-1, 1), 21)
    U3 = np.stack([g2.X * g2.Y, -g2.X ** 2, g2.X + g2.Y])  # z-independent 3-component field
    c = ch02.curl(U3, g2.h)
    assert c.shape == (3, 21, 21)
    assert np.allclose(c[0], 1.0, atol=1e-12) and np.allclose(c[1], -1.0, atol=1e-12)  # ∂u₃/∂x₂ − 0, 0 − ∂u₃/∂x₁
    assert np.allclose(c[2], -2.0 * g2.X - g2.X, atol=1e-12)  # ∂u₂/∂x₁ − ∂u₁/∂x₂ = −2x − x
    assert np.allclose(ch02.curl_components(U3, g2.h), c, atol=1e-14)
    assert np.allclose(ch02.curl(U3[:2], g2.h), c[2], atol=1e-14)  # plane scalar curl
    with pytest.raises(ValueError):
        ch02.curl_components(U3[:2], g2.h)
    # D12 (★): the expander produces (2.25) from (2.24)
    X = ch02.coordinates(3)
    Uf = ch02.field_functions("u", 1)
    res = ch02.expand_indices("eps_ijk u_k,j")
    exp = [sp.diff(Uf[2], X[1]) - sp.diff(Uf[1], X[2]), sp.diff(Uf[0], X[2]) - sp.diff(Uf[2], X[0]), sp.diff(Uf[1], X[0]) - sp.diff(Uf[0], X[1])]
    assert all(sp.simplify(res[i] - exp[i]) == 0 for i in range(3))


def test_curl_V1_of_gradient_and_divergence_of_curl_vanish():  # V1  (Exercises 2.19–2.20)
    # The discrete partials along different axes commute exactly (linear stencils on separate axes), so the
    # identities ∇×∇φ = 0 and ∇·(∇×u) = 0 hold to ROUND-OFF on the grid, not merely to truncation.
    f = ch02.smooth_test_field(3)
    phi = ch02.smooth_scalar_field(3)
    for n in (16, 32, 64):
        g = ch02.grid((-1, 1), n)
        assert np.max(np.abs(ch02.curl(ch02.gradient(phi(*g.coords), g.h), g.h))) < 1e-11
        assert np.max(np.abs(ch02.divergence(ch02.curl(f(*g.coords), g.h), g.h))) < 1e-11
        assert np.max(np.abs(ch02.divergence(ch02.curl(f(*g.coords), g.h, bc="onesided"), g.h, bc="onesided"))) < 1e-11
    gp = ch02.grid((-np.pi, np.pi), 32, periodic=True)
    U = ch02.periodic_test_field()(*gp.coords)
    assert np.max(np.abs(ch02.divergence(ch02.curl(U, gp.h, bc="periodic"), gp.h, bc="periodic"))) < 1e-12
    pf = ch02.potential_field(dim=3)
    g = ch02.grid((-1, 1), 9)
    assert ch02.is_irrotational(pf(*g.coords), g.h, tol=1e-11)[0]  # ∇(x₁x₂x₃): quadratic components → exact
    assert all(c == 0 for c in pf.curl_expr)


def test_curl_V7_rotational_invariance_as_a_vector():  # V7  (proper rotation: curl u' = Cᵀ curl u)
    C = ch02.random_rotation(np.random.default_rng(26))
    f = ch02.smooth_test_field(3)
    errs, hs = [], []
    for n in (17, 33, 65):
        g = ch02.grid((-1, 1), n)
        Xp = np.stack(g.coords)  # new-frame coordinates on the grid
        Xold = np.einsum("ij,j...->i...", C, Xp)  # x = C x'
        Uold = f(*Xold)
        Up = np.einsum("ij,i...->j...", C, Uold)  # u' = Cᵀ u
        curl_p = ch02.curl(Up, g.h)
        exact_p = np.einsum("ij,i...->j...", C, f.curl_fn(*Xold))  # Cᵀ (∇×u)
        errs.append(np.max(np.abs(curl_p - exact_p)))
        hs.append(g.h[0])
    assert errs[-1] < 2e-2 and abs(observed_order(hs, errs) - 2.0) < 0.25
    # a reflection flips the sign of the curl (pseudovector): the 2-D scalar curl of a mirrored field is negated
    g2 = ch02.grid2d((-1, 1), 33)
    f2 = ch02.smooth_test_field(2)
    U = f2(*g2.coords)
    Um = f2(-g2.X, g2.Y)
    Um = np.stack([-Um[0], Um[1]])  # mirror x → −x: u₁ flips
    assert np.max(np.abs(ch02.curl(Um, g2.h) + ch02.curl(U, g2.h)[:, ::-1])) < 1e-12


# =====================================================================================================================
# C12 — symmetric + antisymmetric split (D14, D15); N56–N61; linear-flow kinematics (Part C 6.6)
# =====================================================================================================================
def test_decomposition_V1_unique_parts_and_component_counts():  # V1  (N56; D14)
    rng = np.random.default_rng(27)
    for _ in range(10):
        B = rand_tensor(rng=rng)
        S, A = ch02.symmetric_part(B), ch02.antisymmetric_part(B)
        assert np.allclose(S + A, B, atol=1e-15) and ch02.is_symmetric(S) and ch02.is_antisymmetric(A)
        assert np.allclose(np.diag(A), 0.0)
        assert ch02.independent_components(S) == 6 and ch02.independent_components(A) == 3 and ch02.independent_components(B) == 9
        assert np.allclose(ch02.strain_rate_tensor(B), S) and np.allclose(ch02.rotation_tensor(B), A)
    # uniqueness (D14 steps 4–6): if B = S' + A' with S' symmetric and A' antisymmetric then S' = S
    Bs = ch02.symbol_matrix("B")
    Ss = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"S{min(i, j)}{max(i, j)}"))
    As = Bs - Ss
    eqs = [As[i, j] + As[j, i] for i in range(3) for j in range(i, 3)]  # A' antisymmetric
    sol = sp.solve(eqs, list(Ss.free_symbols), dict=True)[0]
    assert sp.simplify(Ss.subs(sol) - (Bs + Bs.T) / 2) == sp.zeros(3)
    assert ch02.independent_components(np.eye(2)) == 3 and ch02.independent_components([[0, 1], [-1, 0]]) == 1
    # stacked fields: parts of a tensor field
    g = ch02.grid((-1, 1), 5)
    Gf = ch02.vector_gradient(ch02.smooth_test_field(3)(*g.coords), g.h)
    assert np.allclose(ch02.symmetric_part(Gf) + ch02.antisymmetric_part(Gf), Gf)


def test_decomposition_V7_parts_transform_as_tensors():  # V7  (D14: the split is frame-independent)
    rng = np.random.default_rng(28)
    for _ in range(10):
        B, C = rand_tensor(rng=rng), ch02.random_rotation(rng)
        Bp = ch02.transform_tensor(B, C)
        assert np.allclose(ch02.symmetric_part(Bp), ch02.transform_tensor(ch02.symmetric_part(B), C), atol=1e-13)
        assert np.allclose(ch02.antisymmetric_part(Bp), ch02.transform_tensor(ch02.antisymmetric_part(B), C), atol=1e-13)
        # the vector of the antisymmetric part transforms as a vector under a proper rotation
        w = ch02.vector_from_antisymmetric(ch02.antisymmetric_part(B))
        wp = ch02.vector_from_antisymmetric(ch02.antisymmetric_part(Bp))
        assert np.allclose(wp, ch02.transform_vector(w, C), atol=1e-13)


def test_antisymmetric_vector_V1_matrix_2_26_cross_product_and_round_trip():  # V1  (N57, N58; Eqs. (2.26)–(2.27))
    w = np.array([1.0, 2.0, 3.0])
    R = ch02.antisymmetric_from_vector(w)
    assert np.allclose(R, [[0, -3, 2], [3, 0, -1], [-2, 1, 0]])  # Eq. (2.26)
    assert ch02.is_antisymmetric(R)
    assert np.allclose(ch02.vector_from_antisymmetric(R), w, atol=1e-15)  # Eq. (2.27), lower limit i = 1
    rng = np.random.default_rng(29)
    for _ in range(20):
        om, x = rng.standard_normal(3), rng.standard_normal(3)
        Rm = ch02.antisymmetric_from_vector(om)
        assert np.allclose(Rm @ x, np.cross(om, x), atol=1e-14)  # R·x = ω × x
        assert not np.allclose(-Rm @ x, np.cross(om, x))  # the negated variant gives −ω × x (sign discrimination)
        assert np.allclose(ch02.vector_from_antisymmetric(Rm), om, atol=1e-14)
    # the vector of the antisymmetric part of ∂(b × x)_i/∂x_j is b = ½∇×u (pins G[i, j] = ∂u_i/∂x_j together with the sign)
    b = np.array([0.3, -1.1, 0.8])
    Gb = ch02.solid_body_rotation_field(b).grad_fn(0.1, 0.2, 0.3)
    assert np.allclose(ch02.vector_from_antisymmetric(ch02.antisymmetric_part(Gb)), b, atol=1e-14)
    assert np.allclose(Gb, ch02.antisymmetric_from_vector(b), atol=1e-14)  # ∂(b × x)/∂x = [b]_×
    g = ch02.grid((-1, 1), 7)
    Gn = ch02.vector_gradient(ch02.evaluate_field(ch02.solid_body_rotation_field(b), g), g.h)
    assert np.allclose(ch02.vector_from_antisymmetric(ch02.antisymmetric_part(Gn[..., 3, 3, 3])), b, atol=1e-13)
    # plane case: a scalar ω₃
    assert np.allclose(ch02.antisymmetric_from_vector(2.0), [[0, -2], [2, 0]])
    assert ch02.vector_from_antisymmetric([[0, -2], [2, 0]]) == 2.0


def test_antisymmetric_vector_V2_derivation():  # V2  D15 (★★): antisymmetry, entries, inversion, R·x = ω × x
    w = ch02.symbol_vector("omega")
    R = ch02.expand_indices("-eps_ijk omega_k")  # start: R_ij = −ε_ijk ω_k
    Rm = sp.Matrix(3, 3, lambda i, j: R[i, j])
    assert sp.simplify(Rm + Rm.T) == sp.zeros(3)  # step 1
    assert Rm[0, 1] == -w[2] and Rm[0, 2] == w[1] and Rm[1, 2] == -w[0]  # step 2: (2.26)
    assert Rm == sp.Matrix([[0, -w[2], w[1]], [w[2], 0, -w[0]], [-w[1], w[0], 0]])
    Rs = ch02.symbol_matrix("R")
    # steps 3–4: ε_ijl R_ij = −2 ω_l when R_ij = −ε_ijk ω_k
    contr = ch02.expand_indices("eps_ijl R_ij")
    subsR = {Rs[i, j]: Rm[i, j] for i in range(3) for j in range(3)}
    assert all(sp.expand(contr[l].subs(subsR) + 2 * w[l]) == 0 for l in range(3))
    # step 5: ω_k = −½ ε_ijk R_ij recovers ω (sum from i = 1, not "i − 1")
    om = ch02.expand_indices("-1/2 eps_ijk R_ij")
    assert all(sp.expand(om[k].subs(subsR) - w[k]) == 0 for k in range(3))
    # steps 6–7: R_ij x_j = ε_ikj ω_k x_j = (ω × x)_i
    x = ch02.symbol_vector("x")
    step6 = ch02.expand_indices("eps_ikj omega_k x_j")
    assert all(sp.expand((Rm * x)[i] - step6[i]) == 0 for i in range(3))
    assert sp.simplify(Rm * x - w.cross(x)) == sp.zeros(3, 1)
    # the opposite sign convention R_ij = +ε_ijk ω_k would give −ω × x, not ω × x
    assert sp.simplify((-Rm) * x - w.cross(x)) != sp.zeros(3, 1)
    assert sp.simplify((-Rm) * x + w.cross(x)) == sp.zeros(3, 1)


def test_symmetric_double_contraction_V1_2_28_2_29():  # V1  (N59, N60, N61)
    rng = np.random.default_rng(30)
    for _ in range(10):
        tau, B = rand_tensor(symmetric=True, rng=rng), rand_tensor(rng=rng)
        P, PS, PA = ch02.symmetric_double_contraction(tau, B)
        assert abs(PA) < 1e-13  # Eq. (2.29): τ_ij A_ij = 0 for symmetric τ
        assert abs(P - PS) < 1e-13 and abs(P - ch02.double_dot(tau, B, "frobenius")) < 1e-13
        assert abs(P - ch02.double_dot(tau, B, "book")) < 1e-13  # coincide for symmetric τ
        Pn, PSn, PAn = ch02.symmetric_double_contraction(rand_tensor(rng=rng), B)
        assert abs(PAn) > 1e-3  # discrimination: a non-symmetric τ sees the antisymmetric part


def test_linear_flow_V1_kinematics_presets_and_material_lines():  # V1  (Part C 6.6; E3)
    Gr = ch02.velocity_gradient_preset("solid_body_rotation", 0.5)
    assert np.allclose(Gr, [[0, -0.5], [0.5, 0]])
    M = ch02.linear_flow_map(Gr, 1.2)
    assert np.allclose(M, ch02.rotation_matrix_2d(0.6), atol=1e-14)  # rotation by ω₃ t
    assert abs(ch02.material_line_angle(Gr, 1.2, 0.3) - 0.9) < 1e-14
    Gu = ch02.velocity_gradient_preset("uniaxial_extension", 0.7)
    assert np.allclose(ch02.linear_flow_map(Gu, 2.0), np.diag([np.exp(1.4), np.exp(-1.4)]), atol=1e-13)
    Gs = ch02.velocity_gradient_preset("simple_shear", 1.0)
    assert np.allclose(Gs, [[0, 1], [0, 0]])
    dt = 1e-6
    assert np.allclose((ch02.linear_flow_map(Gs, dt) - np.eye(2)) / dt, Gs, atol=1e-5)  # d/dt at t = 0 = G
    assert np.allclose(ch02.velocity_gradient_preset("pure_strain", 2.0), ch02.velocity_gradient_preset("irrotational_strain", 2.0))
    G3 = ch02.velocity_gradient_preset("simple_shear", 1.0, dim=3)
    assert G3.shape == (3, 3) and G3[0, 1] == 1.0 and np.all(G3[2] == 0)
    with pytest.raises(ValueError):
        ch02.velocity_gradient_preset("nope")
    pts = ch02.deform_square(Gs, 1.0, n_side=6)
    assert pts.shape == (2, 36)
    ring = ch02.deform_square(Gu, 0.5, n_side=6, boundary_only=True)
    assert ring.shape == (2, 24)
    # area is preserved for a traceless G (det e^{Gt} = e^{t tr G} = 1): shoelace on the boundary
    x, y = ring
    area = 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
    assert abs(area - 4.0) < 1e-12
    assert ch02.deform_square(G3, 0.3).shape[0] == 3
    # Example 2.4's flow: strain_rate_tensor(G) has S₁₂ = Γ/2 for G = [[0, Γ],[0, 0]]
    assert np.allclose(ch02.strain_rate_tensor(Gs), [[0, 0.5], [0.5, 0]])
    assert set(ch02.VELOCITY_GRADIENT_PRESETS) >= {"simple_shear", "solid_body_rotation", "pure_strain", "uniaxial_extension"}


# =====================================================================================================================
# C13 — principal axes of a symmetric tensor (D17); N62 N63
# =====================================================================================================================
def test_principal_axes_V1_eigen_equation_diagonal_form_and_edge_cases():  # V1
    rng = np.random.default_rng(31)
    for _ in range(20):
        tau = rand_tensor(symmetric=True, rng=rng)
        lam, B = ch02.principal_axes(tau)
        assert np.all(np.diff(lam) >= 0)
        assert np.allclose(tau @ B, B @ np.diag(lam), atol=1e-12)  # τ b = λ b
        assert np.allclose(B.T @ B, np.eye(3), atol=1e-13) and abs(np.linalg.det(B) - 1.0) < 1e-12
        C, tp = ch02.diagonalize(tau)
        assert np.allclose(C, B) and np.allclose(tp, np.diag(lam), atol=1e-12)
        assert np.allclose(tp, ch02.transform_tensor(tau, C), atol=1e-13) and np.allclose(tp, C.T @ tau @ C, atol=1e-13)
        assert np.max(np.abs(tp - np.diag(np.diag(tp)))) < 1e-12
    lam, B = ch02.principal_axes(EXAMPLE_TENSOR)
    assert np.allclose(lam, [3 - np.sqrt(3), 3.0, 3 + np.sqrt(3)], atol=1e-12)
    # hydrostatic: all λ = −p, any orthonormal B works
    lam, B = ch02.principal_axes(-2.0 * np.eye(3))
    assert np.allclose(lam, -2.0) and np.allclose(B.T @ B, np.eye(3), atol=1e-14)
    # repeated eigenvalue: compare projectors, not vectors
    lam, B = ch02.principal_axes(np.diag([1.0, 1.0, 2.0]))
    P = B[:, :2] @ B[:, :2].T
    assert np.allclose(P, np.diag([1.0, 1.0, 0.0]), atol=1e-13) and np.allclose(np.abs(B[:, 2]), [0, 0, 1], atol=1e-13)
    with pytest.raises(ValueError):
        ch02.principal_axes(rand_tensor(rng=rng))  # non-symmetric refused
    with pytest.raises(ValueError):
        ch02.principal_axes(np.zeros((2, 3)))


def test_principal_axes_V7_normal_stress_and_shear_bounds():  # V7  (§2.11 fact (4); Mohr)
    rng = np.random.default_rng(32)
    for tau in (EXAMPLE_TENSOR, rand_tensor(symmetric=True, rng=rng), -3.0 * np.eye(3)):
        lam, B = ch02.principal_axes(tau)
        smin, smax, tsmax = ch02.normal_stress_bounds(tau, rng=np.random.default_rng(0), n=4000)
        assert smin >= lam[0] - 1e-12 and smax <= lam[-1] + 1e-12
        assert tsmax <= 0.5 * (lam[-1] - lam[0]) + 1e-12
        if lam[-1] - lam[0] > 1e-6:
            assert smax > lam[-1] - 0.02 * (lam[-1] - lam[0]) and tsmax > 0.45 * (lam[-1] - lam[0])  # bounds are approached
        for k in range(3):  # along a principal axis: σ_n = λ_k and no shear
            s, ts, _ = ch02.normal_shear_stress(tau, B[:, k])
            assert abs(s - lam[k]) < 1e-12 and ts < 1e-12
        # the normal stress is the weighted mean Σ λ_k c_k² (D17 step 13) for random n
        for _ in range(5):
            n = rng.standard_normal(3)
            n /= np.linalg.norm(n)
            c = B.T @ n
            assert abs(n @ tau @ n - np.sum(lam * c ** 2)) < 1e-12 and abs(np.sum(c ** 2) - 1.0) < 1e-13


def test_principal_axes_V2_derivation():  # V2  D17 (★★★): real λ, orthogonal b, CᵀτC diagonal, σ_n = Σλ_k c_k², shear bound
    # steps 1–5 for a general real symmetric 3 × 3: b̄ᵢ τᵢⱼ bⱼ is real for every complex b (so λ = λ̄)
    tau3 = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"t{min(i, j) + 1}{max(i, j) + 1}", real=True))
    xs, ys = sp.symbols("x1:4", real=True), sp.symbols("y1:4", real=True)
    b = sp.Matrix([x + sp.I * y for x, y in zip(xs, ys)])
    Q = sp.expand((b.conjugate().T * tau3 * b)[0])  # step 1's left side
    assert sp.simplify(sp.im(Q)) == 0  # steps 2–4: Q̄ = Q needs τ real (step 3) and symmetric (step 4)
    assert sp.simplify(sp.re((b.conjugate().T * b)[0]) - sum(x ** 2 + y ** 2 for x, y in zip(xs, ys))) == 0  # b̄ᵢbᵢ > 0
    # symmetry is used: for a NON-symmetric real τ the imaginary part does not vanish
    tau_ns = ch02.symbol_matrix("s")
    Qn = sp.expand((b.conjugate().T * tau_ns * b)[0])
    assert sp.simplify(sp.im(Qn)) != 0
    # step 7: b¹·τ·b² = b²·τ·b¹ for symmetric τ (fails without symmetry)
    b1, b2 = ch02.symbol_vector("p"), ch02.symbol_vector("q")
    assert sp.expand((b1.T * tau3 * b2)[0] - (b2.T * tau3 * b1)[0]) == 0
    assert sp.expand((b1.T * tau_ns * b2)[0] - (b2.T * tau_ns * b1)[0]) != 0
    # 2 × 2 construction (the design's check cell, extended): steps 5, 8, 9, 11, 13, 14, 15
    p, q, r, th = sp.symbols("p q r theta", real=True)
    tau = sp.Matrix([[p, q], [q, r]])
    lam = sp.symbols("lambda")
    poly = sp.expand((tau - lam * sp.eye(2)).det())
    disc = sp.expand(poly.coeff(lam, 1) ** 2 - 4 * poly.coeff(lam, 0))
    assert sp.simplify(disc - ((p - r) ** 2 + 4 * q ** 2)) == 0  # step 5: discriminant ≥ 0 → real roots
    lam1, lam2 = [sp.simplify(s) for s in sp.solve(poly, lam)]
    assert sp.simplify(lam1 + lam2 - (p + r)) == 0 and sp.simplify(lam1 * lam2 - (p * r - q ** 2)) == 0  # Vieta (D18 step 8)
    bv1 = sp.Matrix([q, lam1 - p])
    bv2 = sp.Matrix([q, lam2 - p])
    assert sp.simplify((tau * bv1 - lam1 * bv1)) == sp.zeros(2, 1)  # they are eigenvectors
    assert sp.simplify(bv1.dot(bv2)) == 0  # step 8: orthogonal
    unit = lambda b: b / sp.sqrt(b.dot(b))  # noqa: E731  (NOT .norm(): that introduces Abs() and blocks simplify in sympy 1.14)
    C = sp.Matrix.hstack(unit(bv1), unit(bv2))  # step 9
    assert sp.simplify(C.T * C - sp.eye(2)) == sp.zeros(2)
    D = sp.simplify(C.T * tau * C)  # step 11
    assert sp.simplify(D[0, 1]) == 0 and sp.simplify(D[1, 0]) == 0
    assert sp.simplify(D[0, 0] - lam1) == 0 and sp.simplify(D[1, 1] - lam2) == 0
    n = sp.Matrix([sp.cos(th), sp.sin(th)])  # step 12
    c1, c2 = n.dot(C[:, 0]), n.dot(C[:, 1])
    assert sp.simplify(c1 ** 2 + c2 ** 2 - 1) == 0
    sigma_n = (n.T * tau * n)[0]
    assert sp.simplify(sigma_n - (lam1 * c1 ** 2 + lam2 * c2 ** 2)) == 0  # step 13
    # step 14 (bounds) and step 15 (shear = weighted variance ≤ (Δλ/2)²) in the principal frame, symbolic
    l1, l2, ph = sp.symbols("l1 l2 phi", real=True)
    w1, w2 = sp.cos(ph) ** 2, sp.sin(ph) ** 2
    sn = l1 * w1 + l2 * w2
    ts2 = sp.simplify((l1 ** 2 * w1 + l2 ** 2 * w2) - sn ** 2)
    assert sp.simplify(ts2 - (l1 - l2) ** 2 * sp.sin(2 * ph) ** 2 / 4) == 0  # τ_s² = (Δλ)² sin²2φ / 4 ≤ (Δλ/2)²
    assert sp.simplify(sn - l1 * w1 - l2 * (1 - w1)) == 0  # a convex combination of l1 and l2 → between them
    # numeric confirmation of steps 11–15 on a 3 × 3 with the code
    lamn, Bn = ch02.principal_axes(EXAMPLE_TENSOR)
    Cn, tpn = ch02.diagonalize(EXAMPLE_TENSOR)
    assert np.allclose(tpn, np.diag(lamn), atol=1e-12)


def test_example_2_4_V1_pure_shear_principal_axes():  # V1  (N62; Example 2.4 with our own Γ; Γ ≡ S₁₂)
    for Gam in (1.0, 2.5, -0.4):
        d = ch02.example_2_4(Gam)
        assert np.allclose(d["S"], [[0, Gam], [Gam, 0]])
        assert np.allclose(d["lam"], [Gam, -Gam], atol=1e-14)  # the book's order (Γ, −Γ)
        assert np.allclose(d["S"] @ d["B"], d["B"] @ np.diag(d["lam"]), atol=1e-13)
        assert np.allclose(np.abs(d["B"][:, 0]), [1, 1] / np.sqrt(2), atol=1e-14)  # b¹ ∥ (1, 1)
        assert np.allclose(np.abs(d["B"][:, 1]), [1, 1] / np.sqrt(2), atol=1e-14)  # b² ∥ (−1, 1)
        assert abs(np.linalg.det(d["C"]) - 1.0) < 1e-14 and abs(d["angle_deg"] - 45.0) < 1e-12
        assert np.allclose(d["S_prime"], np.diag([Gam, -Gam]), atol=1e-14)
        assert np.allclose(d["S_prime"], d["C"].T @ d["S"] @ d["C"], atol=1e-14)
        assert np.allclose(ch02.strain_rate_tensor(d["G"]), d["S"])  # G = [[0, 2Γ],[0, 0]] has S₁₂ = Γ
    assert abs(ch02.principal_angle_2d([[0, 1], [1, 0]]) - np.pi / 4) < 1e-14
    assert ch02.principal_angle_2d([[3, 0], [0, 1]]) == 0.0
    S = np.array([[3.0, 1.0], [1.0, -1.0]])
    ang = ch02.principal_angle_2d(S)
    lam, B = ch02.principal_axes(S)
    bmax = B[:, 1]
    assert abs(abs(np.cos(ang) * bmax[0] + np.sin(ang) * bmax[1]) - 1.0) < 1e-12  # the angle points along an eigenvector
    assert np.allclose(ch02.rotation_matrix_2d(np.pi / 4), ch02.example_2_4(1.0)["C"], atol=1e-14)


@book_only
def test_example_2_4_V6_book_values():  # V6
    ex = book()["example_2_4"]
    d = ch02.example_2_4(1.0)
    assert np.allclose(d["S"], ex["S_over_Gamma"])
    assert np.allclose(d["lam"], ex["eigenvalues_over_Gamma"], atol=1e-14)
    assert np.allclose(d["B"][:, 0], ex["b1"], atol=1e-12) and np.allclose(d["B"][:, 1], ex["b2"], atol=1e-12)
    assert np.allclose(d["C"], ex["C_matrix"], atol=1e-12)
    assert abs(d["angle_deg"] - ex["rotation_deg"]) < 1e-12
    assert np.allclose(d["S_prime"], ex["S_prime_over_Gamma"], atol=1e-14)


@book_only
def test_example_2_3_V6_book_factors_and_exercise_2_1a():  # V6
    ex = book()["example_2_3"]
    a, b = 1.0, np.array([0.0, 0.0, 1.0])
    e = ch02.example_2_3(a, b)
    assert abs(e["div_radial_grid"] - ex["div_radial_factor"] * a) < 1e-12
    assert np.allclose(e["curl_rotation_grid"], ex["curl_rotation_factor"] * b, atol=1e-12)
    assert e["curl_radial_max"] < 1e-12 and abs(e["div_rotation_max"] - ex["div_rotation"]) < 1e-12
    x = book()["exercise_2_1a"]
    assert ch02.dot(x["b"], x["c"]) == x["b_dot_c"]


# =====================================================================================================================
# C14 — Gauss' theorem (2.30) (D25); N64 N65 — C15 integral definitions (2.31)–(2.33) (D21, D22); N66 N67 N68
# =====================================================================================================================
def test_gauss_theorem_V1_polynomial_fields_on_boxes():  # V1
    lhs, rhs = ch02.divergence_theorem_box(lambda X, Y, Z: (X, Y, Z), (0, 1), 8)
    assert abs(lhs - 3.0) < 1e-12 and abs(rhs - 3.0) < 1e-12
    lhs, rhs = ch02.divergence_theorem_box(lambda X, Y, Z: (X ** 2, 0.0 * Y, 0.0 * Z), (0, 1), 8)
    assert abs(lhs - 1.0) < 1e-12 and abs(rhs - 1.0) < 1e-12
    faces = ch02.flux_through_box_faces(lambda X, Y, Z: (X, Y, Z), (0, 1), 4)
    assert set(faces) == {"+x", "-x", "+y", "-y", "+z", "-z"}
    assert faces["+x"] == 1.0 and faces["-x"] == 0.0 and abs(sum(faces.values()) - 3.0) < 1e-14
    lhs, rhs = ch02.gauss_gradient_box(lambda X, Y, Z: X * Y * Z, (0, 1), 16)  # scalar Q: Eq. (2.30) with a free index
    assert lhs.shape == (3,) and np.allclose(lhs, 0.25, atol=1e-12) and np.allclose(rhs, 0.25, atol=1e-12)
    assert abs(ch02.volume_integral_box(lambda X, Y, Z: X * Y * Z, (0, 1), 8) - 0.125) < 1e-12
    # 2-D
    lhs, rhs, f2 = ch02.divergence_theorem_rect2d(lambda X, Y: (X, Y), (0, 1), 16)
    assert abs(lhs - 2.0) < 1e-12 and abs(rhs - 2.0) < 1e-12 and set(f2) == {"+x", "-x", "+y", "-y"}
    lhs, rhs, _ = ch02.divergence_theorem_rect2d(lambda X, Y: (X ** 2, 0.0 * Y), (0, 1), 16)
    assert abs(lhs - 1.0) < 1e-12 and abs(rhs - 1.0) < 1e-12
    # quadrature building blocks
    x, w = ch02.midpoint_nodes(0.0, 1.0, 4)
    assert np.allclose(x, [0.125, 0.375, 0.625, 0.875]) and np.allclose(w, 0.25)
    x, w = ch02.simpson_nodes(0.0, 1.0, 5)
    assert abs(np.sum(w * x ** 3) - 0.25) < 1e-15  # exact for cubics
    x, w = ch02.gauss_legendre_nodes(0.0, 1.0, 3)
    assert abs(np.sum(w * x ** 5) - 1 / 6) < 1e-15  # exact to degree 5
    with pytest.raises(ValueError):
        ch02.simpson_nodes(0.0, 1.0, 4)


@pytest.mark.xfail(strict=True, raises=ValueError,
                   reason="OPEN ITEM O1: gauss_gradient_box's docstring promises a (3, 3) result for a vector Q, but "
                          "integral_theorems._eval cannot broadcast a rank-2 field output (raises in volume_integral_box); "
                          "Part C 5.1 only requires the scalar-Q case, which passes above")
def test_gauss_theorem_V1_vector_Q_gradient_form_promised_by_docstring():  # V1 (documents a discrepancy)
    lhs, rhs = ch02.gauss_gradient_box(lambda X, Y, Z: (X * Y, Y * Z, Z * X), (0, 1), 8)  # vector Q: 3 × 3, [i, j] = ∂Q_j/∂x_i
    assert lhs.shape == (3, 3) and np.allclose(lhs, rhs, atol=1e-11)
    assert np.allclose(np.diag(lhs), 0.5) and abs(lhs[1, 0] - 0.5) < 1e-11 and abs(lhs[0, 1]) < 1e-11  # ∂(xy)/∂y = x → ∫ = ½


def test_gauss_theorem_V3_quadrature_orders():  # V3  (midpoint 2, Simpson 4)
    Q = ch02.smooth_test_field(3)
    hs, errs = [], []
    for n in (8, 16, 32, 64):
        lhs, rhs = ch02.divergence_theorem_box(Q, (-1, 1), n, div_fn=Q.div_fn, rule="midpoint")
        hs.append(2.0 / n)
        errs.append(abs(lhs - rhs))
    p_mid = observed_order(hs, errs)
    print(f"box midpoint order {p_mid:.3f} errs {errs}")
    assert abs(p_mid - 2.0) < ORDER_TOL
    hs, errs = [], []
    for n in (9, 17, 33, 65):
        lhs, rhs = ch02.divergence_theorem_box(Q, (-1, 1), n, div_fn=Q.div_fn, rule="simpson")
        hs.append(2.0 / (n - 1))
        errs.append(abs(lhs - rhs))
    p_simp = observed_order(hs, errs)
    print(f"box simpson order {p_simp:.3f} errs {errs}")
    assert abs(p_simp - 4.0) < ORDER_TOL
    # 2-D rectangle: order 2
    Q2 = ch02.smooth_test_field(2)
    hs, errs = [], []
    for n in (8, 16, 32, 64):
        lhs, rhs, _ = ch02.divergence_theorem_rect2d(Q2, (-1, 1), n, div_fn=Q2.div_fn)
        hs.append(2.0 / n)
        errs.append(abs(lhs - rhs))
    assert abs(observed_order(hs, errs) - 2.0) < ORDER_TOL
    # the fourth-order FD fallback for ∇·Q agrees with the exact divergence to ~1e-10
    lhs_fd, _ = ch02.divergence_theorem_box(Q, (-1, 1), 16)
    lhs_ex, _ = ch02.divergence_theorem_box(Q, (-1, 1), 16, div_fn=Q.div_fn)
    assert abs(lhs_fd - lhs_ex) < 1e-9
    d = ch02.fd_partials(lambda X, Y, Z: X ** 3 * Y, (np.array([0.5]), np.array([2.0]), np.array([0.0])))
    assert abs(d[0, 0] - 3 * 0.25 * 2.0) < 1e-10 and abs(d[1, 0] - 0.125) < 1e-10


@needs_ref
def test_gauss_theorem_V5_sphere_benchmark_8pi_over_3():  # V5  (Wikipedia divergence-theorem example)
    ref = ref_json()["divergence_theorem_sphere_example"]
    F = lambda X, Y, Z: (2.0 * X, Y ** 2, Z ** 2)  # noqa: E731  the page's field
    lhs, rhs = ch02.divergence_theorem_sphere(F, ref["radius"], n=24, div_fn=lambda X, Y, Z: 2.0 * (1.0 + Y + Z))
    assert abs(rhs - ref["flux"]) < 1e-12 * ref["flux"] and abs(lhs - ref["flux"]) < 1e-12 * ref["flux"]
    assert abs(ch02.flux_through_sphere(F, 1.0) - 8 * np.pi / 3) < 1e-12
    lhs_fd, rhs_fd = ch02.divergence_theorem_sphere(F, 1.0, n=16)  # FD divergence fallback
    assert abs(lhs_fd - 8 * np.pi / 3) < 1e-9 and abs(rhs_fd - 8 * np.pi / 3) < 1e-12
    # independent closed form of the same integral: ∭ 2(1 + y + z) dV over the unit ball = 2·(4π/3)
    assert abs(ref["flux"] - 2.0 * 4.0 * np.pi / 3.0) < 1e-15
    # a shifted sphere sees the same divergence integral only through the constant part: ∭ 2(1+y+z) = 8π/3 (1 + c_y + c_z)
    c = (0.3, -0.2, 0.5)
    lhs_c, rhs_c = ch02.divergence_theorem_sphere(F, 1.0, n=24, center=c, div_fn=lambda X, Y, Z: 2.0 * (1.0 + Y + Z))
    assert abs(rhs_c - 8 * np.pi / 3 * (1 + c[1] + c[2])) < 1e-12 and abs(lhs_c - rhs_c) < 1e-12


def test_gauss_theorem_V4_solenoidal_field_has_zero_net_outflux():  # V4  (N64: ∇·(b × x) = 0 ⇒ ∯ n·Q dA = 0 on every closed surface)
    f = ch02.solid_body_rotation_field([0.4, -1.0, 2.0])
    for bounds in ((-1, 1), ((0, 2), (-1, 0.5), (0.2, 0.7))):
        lhs, rhs = ch02.divergence_theorem_box(f, bounds, 12, div_fn=f.div_fn)
        assert abs(lhs) < 1e-13 and abs(rhs) < 1e-13
    assert abs(ch02.flux_through_sphere(f, 1.3, center=(0.2, 0.1, -0.4))) < 1e-12
    f2 = ch02.solid_body_rotation_field(1.0, dim=2)
    lhs, rhs, faces = ch02.divergence_theorem_rect2d(f2, ((-1, 2), (0, 1)), 16, div_fn=f2.div_fn)
    assert abs(rhs) < 1e-13 and abs(lhs) < 1e-13
    t = ch02.divergence_theorem_tiled(f2, ((-1, 2), (0, 1)), 3, 8)
    assert abs(t.sum_tiles) < 1e-13 and abs(t.outer) < 1e-13 and abs(t.interior) < 1e-13
    # a point source: the net outflux is m for every box enclosing it, 0 otherwise (the field is solenoidal away from 0)
    src = ch02.point_source_field(2.0)
    assert abs(ch02.divergence_theorem_rect2d(src, (-1, 1), 800)[1] - 2.0) < 1e-6  # midpoint, h = 0.0025: O(h²) ≈ 3e-7
    assert abs(ch02.divergence_theorem_rect2d(src, ((-1, 3), (-2, 0.5)), 800)[1] - 2.0) < 1e-6
    assert abs(ch02.divergence_theorem_rect2d(src, ((1, 2), (1, 2)), 256, div_fn=src.div_fn)[1]) < 1e-6  # not enclosing: 0 (O(h²) quadrature)
    assert src.div_expr == 0 and src.singular_at == (0.0, 0.0)


def test_gauss_theorem_V2_derivation():  # V2  D25 (★★): FTC on a box face pair; interior faces cancel under tiling
    X = ch02.coordinates(3)
    a1, b1, a2, b2, a3, b3, c = sp.symbols("a1 b1 a2 b2 a3 b3 c", real=True)
    cs = sp.symbols("k0:6", real=True)
    Q = cs[0] + cs[1] * X[0] ** 2 * X[1] + cs[2] * X[1] * X[2] ** 3 + cs[3] * X[0] * X[2] + cs[4] * X[0] ** 3 + cs[5] * X[1] ** 2
    box = [(X[0], a1, b1), (X[1], a2, b2), (X[2], a3, b3)]
    for i in range(3):
        lhs = sp.integrate(sp.diff(Q, X[i]), *box)  # step 1: iterated volume integral of ∂Q/∂x_i
        others = [box[k] for k in range(3) if k != i]
        inner = sp.integrate(sp.diff(Q, X[i]), box[i])  # step 2: FTC along x_i
        assert sp.expand(inner - (Q.subs(X[i], box[i][2]) - Q.subs(X[i], box[i][1]))) == 0
        rhs = sp.integrate(Q.subs(X[i], box[i][2]), *others) - sp.integrate(Q.subs(X[i], box[i][1]), *others)  # steps 3–5: n_i = ±1 on two faces, 0 on four
        assert sp.expand(lhs - rhs) == 0
    # steps 7–8: split the box at x1 = c: the two sub-boxes' surface terms add to the big box's (the shared face cancels)
    i = 0
    others = [box[1], box[2]]
    face = lambda val: sp.integrate(Q.subs(X[0], val), *others)  # noqa: E731
    left = face(c) - face(a1)
    right = face(b1) - face(c)
    assert sp.expand(left + right - (face(b1) - face(a1))) == 0
    # numerical tiling: interior faces sum to round-off, outer equals the untiled boundary flux (Part C 5.3)
    Q2 = ch02.smooth_test_field(2)
    t = ch02.divergence_theorem_tiled(Q2, (-1, 1), 4, 16)
    assert abs(t.interior) < 1e-13
    lhs, rhs, _ = ch02.divergence_theorem_rect2d(Q2, (-1, 1), 64, div_fn=Q2.div_fn)
    assert abs(t.outer - rhs) < 1e-13 and abs(t.sum_tiles - t.outer - t.interior) < 1e-15
    t1 = ch02.divergence_theorem_tiled(Q2, (-1, 1), 1, 64)
    assert abs(t1.outer - sum(ch02.flux_through_faces(Q2, (-1, 1), 64).values())) < 1e-14 and t1.interior == 0.0
    assert isinstance(t, tuple) and len(t) == 3 and isinstance(t, ch02.TiledFlux)
    with pytest.raises(ValueError):
        ch02.divergence_theorem_tiled(Q2, (-1, 1), 0, 8)


def test_integral_definitions_V1_exact_for_linear_fields_and_index_order():  # V1  (C15; Eqs. (2.31)–(2.33))
    a, b = 1.7, np.array([0.3, -0.8, 1.1])
    x0 = np.array([0.4, -0.3, 0.9])
    for h in (1.0, 0.3, 0.05):
        assert abs(ch02.integral_divergence(ch02.radial_field(a), x0, h) - 3.0 * a) < 1e-12
        assert np.allclose(ch02.integral_curl(ch02.solid_body_rotation_field(b), x0, h), 2.0 * b, atol=1e-12)
        assert abs(ch02.integral_divergence(ch02.solid_body_rotation_field(b), x0, h)) < 1e-12
        lin = lambda X, Y, Z: 2.0 * X - 3.0 * Y + 0.5 * Z + 1.0  # noqa: E731
        assert np.allclose(ch02.integral_gradient(lin, x0, h), [2.0, -3.0, 0.5], atol=1e-12)
        assert abs(ch02.integral_divergence_2d(ch02.radial_field(a, dim=2), x0[:2], h) - 2.0 * a) < 1e-12
    # Q = (x², 0, 0) at x₀ = 1, h = 0.2: front 1.21·0.04 − back 0.81·0.04 = 0.016; /0.008 = 2 = 2x₀ exactly (D21/D22 check)
    Qx2 = lambda X, Y, Z: (X ** 2, 0.0 * Y, 0.0 * Z)  # noqa: E731
    assert abs(ch02.integral_divergence(Qx2, [1.0, 0.0, 0.0], 0.2) - 2.0) < 1e-13
    # Q = (x³, 0, 0): 3 + h²/4 (D22's remainder)
    Qx3 = lambda X, Y, Z: (X ** 3, 0.0 * Y, 0.0 * Z)  # noqa: E731
    assert abs(ch02.integral_divergence(Qx3, [1.0, 0.0, 0.0], 0.2) - (3.0 + 0.04 / 4)) < 1e-13
    # index order of integral_gradient for a vector Q: [i, j] = ∂Q_j/∂x_i (derivative index first, as n_i Q_j)
    Qv = lambda X, Y, Z: (Y, 0.0 * Y, 0.0 * Z)  # noqa: E731  Q₁ = y → ∂Q₁/∂x₂ = 1
    G = ch02.integral_gradient(Qv, x0, 0.1)
    assert G.shape == (3, 3) and abs(G[1, 0] - 1.0) < 1e-12 and abs(G[0, 1]) < 1e-12
    Qs = ch02.smooth_test_field(3)
    assert np.max(np.abs(ch02.integral_gradient(Qs, x0, 0.01) - Qs.grad_fn(*x0).T)) < 1e-5


@pytest.mark.parametrize("kind", ["divergence", "curl", "gradient", "curl_component"])
def test_integral_definitions_V3_second_order_in_h(kind):  # V3  (D21/D22: error ∝ h²)
    r = ch02.integral_definition_convergence(kind)
    print(f"integral {kind}: order {r['order']:.3f} err {r['err']}")
    assert abs(r["order"] - 2.0) < ORDER_TOL
    assert np.all(np.diff(r["err"]) < 0)


def test_integral_definitions_V1_agree_with_grid_operators():  # V1  (the two routes to ∇·, ∇× meet at a node)
    f = ch02.smooth_test_field(3)
    g = ch02.grid((-1, 1), 41)
    k, j, i = 20, 14, 27
    x0 = np.array([g.x[i], g.y[j], g.z[k]])
    div_grid = ch02.divergence(f(*g.coords), g.h)[k, j, i]
    curl_grid = ch02.curl(f(*g.coords), g.h)[:, k, j, i]
    div_int = ch02.integral_divergence(f, x0, g.h[0])
    curl_int = ch02.integral_curl(f, x0, g.h[0])
    exact_d, exact_c = float(f.div_fn(*x0)), f.curl_fn(*x0)
    assert abs(div_grid - exact_d) < 2e-3 and abs(div_int - exact_d) < 2e-3
    assert np.max(np.abs(curl_grid - exact_c)) < 2e-3 and np.max(np.abs(curl_int - exact_c)) < 2e-3
    assert abs(div_grid - div_int) < 3e-3 and np.max(np.abs(curl_grid - curl_int)) < 3e-3


def test_integral_divergence_V2_derivation():  # V2  D21 + D22 (★★): the box outflux per volume is ∂Q_i/∂x_i + O(Δ²)
    X = ch02.coordinates(3)
    d1, d2, d3 = sp.symbols("d1 d2 d3", positive=True)
    D = [d1, d2, d3]
    cs = sp.symbols("c0:10", real=True)
    Qv = [cs[0] + cs[1] * X[0] + cs[2] * X[1] + cs[3] * X[0] ** 2 + cs[4] * X[0] * X[2] + cs[5] * X[0] ** 3,
          cs[6] * X[1] ** 2 + cs[7] * X[0] * X[1],
          cs[8] * X[2] ** 3 + cs[9] * X[1] * X[2]]
    s, t = sp.symbols("s t", real=True)
    V = d1 * d2 * d3

    def face_pair(d):
        o = [k for k in range(3) if k != d]
        total = 0
        for sign in (1, -1):  # steps 1, 3: n = ±e_d at x_d ± Δ_d/2, face integral exact
            sub = {X[d]: X[d] + sign * D[d] / 2, X[o[0]]: X[o[0]] + s, X[o[1]]: X[o[1]] + t}
            total += sign * sp.integrate(Qv[d].subs(sub, simultaneous=True), (s, -D[o[0]] / 2, D[o[0]] / 2), (t, -D[o[1]] / 2, D[o[1]] / 2))
        return sp.expand(total)

    pairs = [face_pair(d) for d in range(3)]
    # step 5: each opposite-face pair gives (∂Q_d/∂x_d) V + higher order; the zeroth-order terms have cancelled
    for d in range(3):
        lead = sp.expand(pairs[d] / V)
        assert sp.expand(lead.subs({d1: 0, d2: 0, d3: 0}) - sp.diff(Qv[d], X[d])) == 0
        assert sp.expand(pairs[d]).subs({d1: 0, d2: 0, d3: 0}) == 0  # no Q(x)·area term survives
    # steps 7–9: sum of six faces / V → ∂Q_i/∂x_i with an O(Δ²) remainder (D21 step 4: the limit is the divergence)
    ratio = sp.expand(sum(pairs) / V)
    div = sum(sp.diff(q, x) for q, x in zip(Qv, X))
    rem = sp.expand(ratio - div)
    assert sp.expand(rem - (cs[5] * d1 ** 2 / 4 + cs[8] * d3 ** 2 / 4)) == 0  # exactly the cubic terms' h²/4
    assert sp.limit(sp.limit(sp.limit(rem, d1, 0), d2, 0), d3, 0) == 0
    # step 6 of D21: n·Q for the vector field; step 7 of D21: ε_kij n_i Q_j = (n × Q)_k
    n, Qs = ch02.symbol_vector("n"), ch02.symbol_vector("Q")
    assert sp.expand(ch02.expand_indices("n_i Q_i") - n.dot(Qs)) == 0
    nxq = ch02.expand_indices("eps_kij n_i Q_j")
    assert all(sp.expand(nxq[k] - n.cross(Qs)[k]) == 0 for k in range(3))
    # D21 step 2 (mean-value theorem) for the one-dimensional model: ∫ f' dx over [x0 − h/2, x0 + h/2] = f' (x*) h with x* inside
    xx, x0, hh = sp.symbols("x x0 h", real=True)
    fpoly = xx ** 3
    mean = sp.integrate(sp.diff(fpoly, xx), (xx, x0 - hh / 2, x0 + hh / 2)) / hh
    assert sp.expand(mean - (3 * x0 ** 2 + hh ** 2 / 4)) == 0  # = f'(x*) for x* = ±h/(2√3) + x0 ∈ V; → f'(x0) as h → 0


# =====================================================================================================================
# C16 — Stokes' theorem (2.34)–(2.35) (D26); N69 N70 N71 N73
# =====================================================================================================================
def test_stokes_V1_orientation_right_hand_rule_and_geometry():  # V1  (N69: t = n_c × n; ∮(x − c) × t ds = 2A n)
    rng = np.random.default_rng(33)
    # Fig. 2.10: n_c is ⊥ C, tangent to A and points INTO A (up the cap); at the rim point (1, 0, 0) of a flat disc
    # with n = e₃, n_c = −e₁ (towards the centre) and t = n_c × n = +e₂ — counterclockwise seen from the tip of n.
    n_c, n = np.array([-1.0, 0, 0]), np.array([0, 0, 1.0])
    t = ch02.boundary_tangent(n_c, n)
    assert np.allclose(t, [0, 1.0, 0]) and abs(np.linalg.det(np.column_stack([n_c, n, t])) - 1.0) < 1e-14  # right-handed triad
    # (an OUTWARD n_c would reverse t — the docstrings of boundary_tangent/planar_loop say "outward": see the report)
    assert np.allclose(ch02.boundary_tangent(-n_c, n), -t)
    assert np.allclose(ch02.boundary_tangent([-1.0, 0.0, 0.0], [0.0, 0.0, 1.0])[1], 1.0)  # the E5 parity row
    for _ in range(10):
        c, nn = rng.standard_normal(3), rng.standard_normal(3)
        nn /= np.linalg.norm(nn)
        L = ch02.planar_loop(c, nn, radius=0.7, n=400)
        moment = np.sum(np.cross(L.points - L.center, L.tangents) * L.ds[:, None], axis=0)
        assert np.allclose(moment, 2.0 * L.area * nn, atol=1e-12)  # counterclockwise about n
        assert abs(L.length - 2 * np.pi * 0.7) < 1e-12 and abs(L.area - np.pi * 0.49) < 1e-15 and L.dim == 3
        assert np.allclose(np.linalg.norm(L.tangents, axis=1), 1.0) and np.allclose(L.tangents @ nn, 0.0, atol=1e-13)
        # the tangent is n_c × n with n_c pointing into A (towards the centre), as in Fig. 2.10
        n_c = -(L.points - L.center) / 0.7
        assert np.allclose(L.tangents, np.cross(n_c, nn), atol=1e-12)
        assert np.allclose(L.tangents, np.cross(nn, (L.points - L.center) / 0.7), atol=1e-12)  # = n × r̂
        # tangents equal the derivative of the points along the loop (central difference of a periodic sample)
        dp = np.roll(L.points, -1, 0) - np.roll(L.points, 1, 0)
        assert np.allclose(dp / np.linalg.norm(dp, axis=1, keepdims=True), L.tangents, atol=1e-4)
        R = ch02.rectangle_loop(c, nn, a=0.5, b=0.8, n=50)
        moment = np.sum(np.cross(R.points - R.center, R.tangents) * R.ds[:, None], axis=0)
        assert np.allclose(moment, 2.0 * 0.4 * nn, atol=1e-12) and abs(R.length - 2.6) < 1e-12 and abs(R.area - 0.4) < 1e-15
        D = ch02.planar_disc(c, nn, radius=0.7)
        assert abs(D.area - np.pi * 0.49) < 1e-12 and np.allclose(D.normals, nn) and D.contains(c) and not D.contains(c + 0.8 * L.tangents[0])
        P = ch02.planar_rectangle(c, nn, a=0.5, b=0.8, n=16)
        assert abs(P.area - 0.4) < 1e-12 and P.contains(c)
    # 2-D convenience: a 2-vector centre → plane z = 0, n = e₃, 2-D points
    L2 = ch02.planar_loop([0.2, -0.1], radius=1.0, n=64)
    assert L2.dim == 2 and L2.points.shape == (64, 2) and np.allclose(L2.normal, [0, 0, 1])
    m2 = np.sum(np.cross(np.c_[L2.points - L2.center, np.zeros(64)], np.c_[L2.tangents, np.zeros(64)]) * L2.ds[:, None], axis=0)
    assert np.allclose(m2, [0, 0, 2 * np.pi], atol=1e-12)
    with pytest.raises(ValueError):
        ch02.planar_loop([0.0, 0.0, 0.0], None)


def test_stokes_V1_solid_body_rotation_and_flipped_normal():  # V1  (C16; the worked number: Γ = 2|b|πR²)
    b = np.array([0.5, -1.0, 2.0])
    f = ch02.solid_body_rotation_field(b)
    nb = b / np.linalg.norm(b)
    for R in (0.5, 1.0, 2.0):
        L = ch02.planar_loop([0.3, 0.1, -0.2], nb, radius=R, n=256)
        D = ch02.planar_disc([0.3, 0.1, -0.2], nb, radius=R)
        expected = 2.0 * np.linalg.norm(b) * np.pi * R ** 2
        assert abs(ch02.circulation(f, L) - expected) < 1e-12 * expected
        assert abs(ch02.curl_flux(f, D, curl_fn=f.curl_fn) - expected) < 1e-12 * expected
        assert abs(ch02.curl_flux(f, D) - expected) < 1e-8 * expected  # FD curl fallback
        chk = ch02.stokes_theorem_check(f, L, D, curl_fn=f.curl_fn)
        assert chk.hypothesis_ok and chk.mismatch < 1e-12 * expected
        lhs, rhs = chk  # iterates as (lhs, rhs)
        assert lhs == chk.lhs and rhs == chk.rhs
        # flipping n flips both sides
        Lm = ch02.planar_loop([0.3, 0.1, -0.2], -nb, radius=R, n=256)
        Dm = ch02.planar_disc([0.3, 0.1, -0.2], -nb, radius=R)
        chk_m = ch02.stokes_theorem_check(f, Lm, Dm, curl_fn=f.curl_fn)
        assert abs(chk_m.lhs + chk.lhs) < 1e-12 * expected and abs(chk_m.rhs + chk.rhs) < 1e-12 * expected
    # 2-D: the plane rotation field on the unit disc: 2πb₃
    f2 = ch02.solid_body_rotation_field(1.0, dim=2)
    chk2 = ch02.stokes_theorem_check(f2, ch02.planar_loop([0.0, 0.0], radius=1.0), ch02.planar_disc([0.0, 0.0], radius=1.0), curl_fn=f2.curl_fn)
    assert abs(chk2.lhs - 2 * np.pi) < 1e-12 and abs(chk2.rhs - 2 * np.pi) < 1e-12
    # shear u₁ = Γx₂ on a unit square: circulation −Γ (straight streamlines with curl)
    sh = ch02.shear_field(1.5)
    assert abs(ch02.circulation(sh, ch02.rectangle_loop([0.0, 0.0], a=1.0, b=1.0, n=32)) + 1.5) < 1e-13


def test_stokes_V1_gradient_field_has_zero_circulation():  # V1  (D26 step 10; Exercise 2.20)
    pf = ch02.potential_field()  # ∇(x₁² − x₂²)
    for c, R in (([0.3, 0.1], 0.7), ([-1.0, 2.0], 1.3)):
        assert abs(ch02.circulation(pf, ch02.planar_loop(c, radius=R, n=128))) < 1e-13
        assert abs(ch02.circulation(pf, ch02.rectangle_loop(c, a=R, b=0.5 * R, n=32))) < 1e-13
    X = ch02.coordinates(3)
    pf3 = ch02.potential_field(sp.sin(X[0]) * X[1] + X[2] ** 2, dim=3)
    L = ch02.planar_loop([0.2, 0.3, 0.1], [1.0, 1.0, 0.5], radius=0.8, n=512)
    assert abs(ch02.circulation(pf3, L)) < 1e-10  # midpoint on a circle is spectrally accurate for smooth fields
    assert all(sp.simplify(c) == 0 for c in pf3.curl_expr)


def test_stokes_V3_rectangle_loop_second_order():  # V3  (loop resolution: midpoint rule)
    f = ch02.smooth_test_field(2)
    hs, errs = [], []
    for n in (8, 16, 32, 64):
        L = ch02.rectangle_loop([0.3, 0.2], a=0.8, b=0.6, n=n)
        S = ch02.planar_rectangle([0.3, 0.2], a=0.8, b=0.6, n=n)
        chk = ch02.stokes_theorem_check(f, L, S, curl_fn=f.curl_fn)
        assert chk.hypothesis_ok
        hs.append(1.0 / n)
        errs.append(chk.mismatch)
    p = observed_order(hs, errs)
    print(f"stokes rectangle order {p:.3f} errs {errs}")
    assert abs(p - 2.0) < ORDER_TOL
    # 3-D rectangle loop in a tilted plane against the exact curl flux
    f3 = ch02.smooth_test_field(3)
    hs, errs = [], []
    for n in (8, 16, 32, 64):
        L = ch02.rectangle_loop([0.1, 0.2, 0.3], [1.0, 2.0, 2.0], a=0.6, b=0.4, n=n)
        S = ch02.planar_rectangle([0.1, 0.2, 0.3], [1.0, 2.0, 2.0], a=0.6, b=0.4, n=n)  # both sides refined together
        errs.append(abs(ch02.circulation(f3, L) - ch02.curl_flux(f3, S, curl_fn=f3.curl_fn)))
        hs.append(1.0 / n)
    assert abs(observed_order(hs, errs) - 2.0) < 0.25


def test_stokes_V1_curl_component_2_35_uses_u_y_in_example_2_6():  # V1  (N71, N73; Ex. 2.6 with the u_y correction)
    b = np.array([0.3, -0.8, 1.1])
    f = ch02.solid_body_rotation_field(b)
    for n, h in (([0, 0, 1.0], 0.4), ([1.0, 0, 0], 0.1), ([1.0, 1.0, 0], 0.05)):
        nn = np.asarray(n) / np.linalg.norm(n)
        assert abs(ch02.integral_curl_component(f, [0.2, -0.1, 0.4], nn, h) - 2.0 * (nn @ b)) < 1e-12  # linear: exact
    # Ex. 2.6: rectangle in the x = const plane, n = e_x: (∇×u)_x = ∂u_z/∂y − ∂u_y/∂z. With u = (0, z, 0) the book's
    # misprinted u_z in the second bracket would give 0; the correct u_y gives −1.
    u = lambda X, Y, Z: (0.0 * X, Z, 0.0 * X)  # noqa: E731
    assert abs(ch02.integral_curl_component(u, [0.0, 0.0, 0.0], [1.0, 0, 0], 0.2) + 1.0) < 1e-13
    u2 = lambda X, Y, Z: (0.0 * X, 0.0 * Y, Y)  # noqa: E731  ∂u_z/∂y = 1
    assert abs(ch02.integral_curl_component(u2, [0.0, 0.0, 0.0], [1.0, 0, 0], 0.2) - 1.0) < 1e-13
    # 2-D: a 2-vector x0 means n = e₃
    sh = ch02.shear_field(2.0)
    assert abs(ch02.integral_curl_component(sh, [0.3, 0.2], None, 0.1) + 2.0) < 1e-13
    # smooth field: order 2 (also covered by integral_definition_convergence("curl_component"))
    f3 = ch02.smooth_test_field(3)
    x0 = np.array([0.3, -0.2, 0.5])
    # (the e₂ component is degenerate for this field — u₁ depends on x₁ only and u₃ is bilinear in the x₂ = const plane,
    #  so the midpoint loop is exact there; use n = e₁, where (∇×u)₁ = x₁x₃cos x₂ − cos x₂ cos x₃ is genuinely nonlinear)
    errs = [abs(ch02.integral_curl_component(f3, x0, [1.0, 0, 0], h) - f3.curl_fn(*x0)[0]) for h in (0.4, 0.2, 0.1)]
    assert abs(observed_order([0.4, 0.2, 0.1], errs) - 2.0) < ORDER_TOL
    assert abs(ch02.integral_curl_component(f3, x0, [0, 1.0, 0], 0.4) - f3.curl_fn(*x0)[1]) < 1e-13  # the degenerate case is exact


def test_stokes_V7_singular_core_breaks_the_hypothesis():  # V7  (the irrotational vortex; explainer ⚠️ status)
    K = 1.5
    v = ch02.irrotational_vortex_field(K)
    assert v.curl_expr == 0 and v.singular_at == (0.0, 0.0)  # curl-free away from the core
    for R in (0.5, 1.0, 3.0):
        chk = ch02.stokes_theorem_check(v, ch02.planar_loop([0.0, 0.0], radius=R, n=256), ch02.planar_disc([0.0, 0.0], radius=R), curl_fn=v.curl_fn)
        assert not chk.hypothesis_ok and "singular" in chk.note
        assert abs(chk.rhs - 2 * np.pi * K) < 1e-12 and chk.lhs == 0.0  # circulation 2πK for every loop; curl flux 0
    # a loop that does not enclose the core: hypothesis OK and both sides ≈ 0
    chk = ch02.stokes_theorem_check(v, ch02.planar_loop([3.0, 0.0], radius=1.0, n=512), ch02.planar_disc([3.0, 0.0], radius=1.0), curl_fn=v.curl_fn)
    assert chk.hypothesis_ok and abs(chk.rhs) < 1e-10 and chk.lhs == 0.0
    # the angular speed is K/r (u_θ = K/r): the field on the unit circle has magnitude K
    U = v(np.array([1.0, 0.0]), np.array([0.0, 1.0]))
    assert np.allclose(np.linalg.norm(U, axis=0), K)


def test_stokes_V2_derivation():  # V2  D26 (★★): one rectangle exactly, the signs of the four sides, the corollary
    X = ch02.coordinates(2)
    x1, x2 = X
    d1, d2 = sp.symbols("d1 d2", positive=True)
    cs = sp.symbols("c0:8", real=True)
    u1 = cs[0] + cs[1] * x1 + cs[2] * x2 + cs[3] * x2 ** 2 + cs[4] * x1 * x2 ** 2
    u2 = cs[5] * x1 + cs[6] * x1 ** 2 + cs[7] * x1 ** 3 * x2
    s = sp.symbols("s", real=True)
    # step 1: the four sides, counterclockwise about e₃: bottom +e₁, right +e₂, top −e₁, left −e₂
    bottom = sp.integrate(u1.subs({x1: x1 + s, x2: x2 - d2 / 2}, simultaneous=True), (s, -d1 / 2, d1 / 2))
    top = -sp.integrate(u1.subs({x1: x1 + s, x2: x2 + d2 / 2}, simultaneous=True), (s, -d1 / 2, d1 / 2))
    right = sp.integrate(u2.subs({x1: x1 + d1 / 2, x2: x2 + s}, simultaneous=True), (s, -d2 / 2, d2 / 2))
    left = -sp.integrate(u2.subs({x1: x1 - d1 / 2, x2: x2 + s}, simultaneous=True), (s, -d2 / 2, d2 / 2))
    A = d1 * d2
    # steps 3–4: bottom + top = −∂u₁/∂x₂ ΔA + O(Δ⁴)
    bt = sp.expand((bottom + top) / A)
    assert sp.expand(bt.subs({d1: 0, d2: 0}) + sp.diff(u1, x2)) == 0
    # step 5: right + left = +∂u₂/∂x₁ ΔA + O(Δ⁴)
    rl = sp.expand((right + left) / A)
    assert sp.expand(rl.subs({d1: 0, d2: 0}) - sp.diff(u2, x1)) == 0
    # step 6: the whole loop / ΔA → (∇×u)₃ = ∂u₂/∂x₁ − ∂u₁/∂x₂ (Eq. (2.35)); polynomial fields: exact Stokes on the rectangle
    circ = sp.expand(bottom + right + top + left)
    curl3 = sp.diff(u2, x1) - sp.diff(u1, x2)
    t = sp.symbols("t", real=True)
    flux = sp.integrate(curl3.subs({x1: x1 + s, x2: x2 + t}, simultaneous=True), (s, -d1 / 2, d1 / 2), (t, -d2 / 2, d2 / 2))
    assert sp.expand(circ - flux) == 0  # Eq. (2.34) exactly for this rectangle
    assert sp.expand(sp.expand(circ / A).subs({d1: 0, d2: 0}) - curl3) == 0  # Eq. (2.35): the limit Δ → 0
    # step 8: two rectangles sharing an edge — the shared edge is walked both ways and cancels
    def loop(xc, yc, a, b):
        bo = sp.integrate(u1.subs({x1: xc + s, x2: yc - b / 2}, simultaneous=True), (s, -a / 2, a / 2))
        to = -sp.integrate(u1.subs({x1: xc + s, x2: yc + b / 2}, simultaneous=True), (s, -a / 2, a / 2))
        ri = sp.integrate(u2.subs({x1: xc + a / 2, x2: yc + s}, simultaneous=True), (s, -b / 2, b / 2))
        le = -sp.integrate(u2.subs({x1: xc - a / 2, x2: yc + s}, simultaneous=True), (s, -b / 2, b / 2))
        return bo + to + ri + le
    two = loop(x1 - d1 / 4, x2, d1 / 2, d2) + loop(x1 + d1 / 4, x2, d1 / 2, d2)
    assert sp.expand(two - circ) == 0
    # step 10: ∮ ∇φ·t ds = 0 for a polynomial φ, hence ∇×∇φ = 0
    phi = cs[0] * x1 ** 3 * x2 + cs[1] * x2 ** 2 + cs[2] * x1 * x2
    g1, g2 = sp.diff(phi, x1), sp.diff(phi, x2)
    circ_phi = (sp.integrate(g1.subs({x1: x1 + s, x2: x2 - d2 / 2}, simultaneous=True), (s, -d1 / 2, d1 / 2))
                - sp.integrate(g1.subs({x1: x1 + s, x2: x2 + d2 / 2}, simultaneous=True), (s, -d1 / 2, d1 / 2))
                + sp.integrate(g2.subs({x1: x1 + d1 / 2, x2: x2 + s}, simultaneous=True), (s, -d2 / 2, d2 / 2))
                - sp.integrate(g2.subs({x1: x1 - d1 / 2, x2: x2 + s}, simultaneous=True), (s, -d2 / 2, d2 / 2)))
    assert sp.expand(circ_phi) == 0
    assert sp.simplify(sp.diff(g2, x1) - sp.diff(g1, x2)) == 0


# =====================================================================================================================
# Part C contract — every callable the notebook and the explainers use exists and is scalar-callable where promised
# =====================================================================================================================
PART_C_NAMES = [
    # C.1 index_notation
    "expand_indices", "expand_indices_str", "classify_indices", "tensor_order", "rename_dummy", "comma_to_partial",
    "coordinates", "symbol_vector", "symbol_matrix", "symbol_tensor", "field_functions",
    # C.2 tensors
    "unit_vectors", "vector_from_components", "dot", "inner", "outer", "tensor_product", "contract", "dot_tensor_vector",
    "double_dot", "trace", "direction_cosines", "rotation_matrix_2d", "rotation_matrix_3d", "rotation_angle",
    "random_rotation", "is_orthogonal", "is_proper_rotation", "orthogonality_residual", "transform_vector",
    "inverse_transform_vector", "transform_tensor", "transforms_as_vector", "transforms_as_tensor", "cube_face_tractions",
    "stress_component_meaning", "tetrahedron_face_areas", "traction", "traction_components", "normal_shear_stress",
    "mohr_circle_2d", "kronecker_delta", "levi_civita", "permutation_sign", "epsilon_delta_residual", "triple_product",
    "is_isotropic", "cross", "cross_einsum", "angle_between", "is_symmetric", "is_antisymmetric", "independent_components",
    "symmetric_part", "antisymmetric_part", "strain_rate_tensor", "rotation_tensor", "antisymmetric_from_vector",
    "vector_from_antisymmetric", "symmetric_double_contraction", "principal_axes", "characteristic_polynomial",
    "diagonalize", "invariants", "normal_stress_bounds", "STRESS_CUBE_FACES",
    # C.3 grids, C.4 operators
    "grid", "grid2d", "evaluate_field", "axis_of_direction", "spacing", "partial", "second_partial", "gradient",
    "directional_derivative", "divergence", "vector_gradient", "tensor_divergence", "curl", "curl_components",
    "is_solenoidal", "is_irrotational", "laplacian",
    # C.5 integral theorems
    "volume_integral_box", "gauss_gradient_box", "divergence_theorem_box", "flux_through_box_faces", "flux_through_sphere",
    "divergence_theorem_sphere", "divergence_theorem_rect2d", "flux_through_faces", "integral_divergence_2d",
    "divergence_theorem_tiled", "TiledFlux", "integral_gradient", "integral_divergence", "integral_curl", "boundary_tangent",
    "planar_loop", "rectangle_loop", "planar_disc", "planar_rectangle", "Loop", "Surface", "circulation", "curl_flux",
    "stokes_theorem_check", "StokesCheck", "integral_curl_component", "fd_partials", "midpoint_nodes", "simpson_nodes",
    "gauss_legendre_nodes",
    # C.6 chapter module
    "polar_components", "cartesian_from_polar", "polar_components_deg", "shear_flow_stress", "example_2_2", "traction_2d",
    "stress_vs_angle", "principal_angle_2d", "example_2_4", "radial_field", "solid_body_rotation_field", "shear_field",
    "irrotational_vortex_field", "potential_field", "point_source_field", "smooth_test_field", "periodic_test_field",
    "periodic_scalar_field", "smooth_scalar_field", "exact_div_curl", "example_2_3", "linear_flow_map", "deform_square",
    "velocity_gradient_preset", "material_line_angle", "operator_convergence", "integral_definition_convergence",
    "VectorField", "ScalarField", "VELOCITY_GRADIENT_PRESETS",
]


def test_part_c_contract_V1_every_callable_exists_and_parity_rows_return_floats():  # V1 (smoke)
    missing = [n for n in PART_C_NAMES if not hasattr(ch02, n)]
    assert not missing, missing
    assert set(PART_C_NAMES) <= set(ch02.__all__)
    # explainer parity rows: plain lists / floats in, one number out
    assert isinstance(ch02.traction_2d([[0, 1], [1, 0]], 0.5236)["sigma_n"], float)
    assert isinstance(ch02.principal_axes([[0, 1], [1, 0]])[0][1], float)
    assert isinstance(ch02.polar_components(1.0, 2.0, 0.5)[0], float)
    assert isinstance(ch02.double_dot([[1, 2], [3, 4]], [[1, 0], [0, 1]], "book"), float)
    assert isinstance(ch02.circulation(ch02.solid_body_rotation_field(1.0, dim=2), ch02.planar_loop([0.0, 0.0], radius=1.0)), float)
    assert isinstance(ch02.mohr_circle_2d([[3, 1], [1, -1]])[1], float)
    assert isinstance(ch02.rotation_angle([[0, -1], [1, 0]]), float)
    assert isinstance(ch02.epsilon_delta_residual(), float)
    assert isinstance(ch02.integral_divergence(ch02.radial_field(1.0), [0.0, 0.0, 0.0], 0.1), float)
    assert isinstance(ch02.vector_from_antisymmetric([[0, -2], [2, 0]]), float)
    assert isinstance(ch02.material_line_angle([[0, -1], [1, 0]], 0.5, 0.0), float)
    # symbolic helpers
    assert ch02.symbol_tensor("A", 3).shape == (3, 3, 3) and ch02.symbol_tensor("A", 0) == sp.Symbol("A", real=True)
    assert len(ch02.field_functions("u", 2)) == 3 and ch02.field_functions("phi", 0).func.__name__ == "phi"
    with pytest.raises(ValueError):
        ch02.field_functions("u", 3)
    assert len(ch02.coordinates(2)) == 2
    # VectorField / ScalarField basics
    X = ch02.coordinates(2)
    with pytest.raises(ValueError):
        ch02.VectorField([X[0]], X)
    v = ch02.VectorField([X[1], -X[0]], X, name="v")
    with pytest.raises(ValueError):
        v(np.zeros(3), np.zeros(3), np.zeros(3))  # 2-D field needs 2 coordinate arrays
    assert "VectorField" in repr(v) and "ScalarField" in repr(ch02.smooth_scalar_field(2))
    pf = ch02.periodic_test_field()
    ps = ch02.periodic_scalar_field()
    assert np.allclose(pf(np.pi, 0.3, 0.2), pf(-np.pi, 0.3, 0.2)) and abs(ps(np.pi, 0.1, 0.2) - ps(-np.pi, 0.1, 0.2)) < 1e-15
    st = ch02.smooth_test_field(3)
    assert st.grad_fn(0.1, 0.2, 0.3).shape == (3, 3) and st.curl_fn(0.1, 0.2, 0.3).shape == (3,)


@needs_ref
def test_reference_V5_files_present_and_documented():  # V5 (bookkeeping)
    ref = ref_json()
    assert {"divergence_theorem_sphere_example", "levi_civita_identities", "cauchy_stress_tensor", "rotation_matrix", "stokes_theorem"} <= set(ref)
    src = (REF / "SOURCES.md").read_text(encoding="utf-8")
    for key in ("Divergence_theorem", "Levi-Civita_symbol", "Cauchy_stress_tensor", "Rotation_matrix", "Stokes"):
        assert key in src
    assert abs(ref["divergence_theorem_sphere_example"]["flux"] - 8 * np.pi / 3) < 1e-15
