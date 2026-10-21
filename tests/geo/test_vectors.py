"""Evidence tests for ``fluidpy.geo.vectors`` (AOFD PoC A and G: dot product, cross product, angle, projection, area,
triple product, k x u).

Levels (verify-implementation ladder): V1 analytic truth on fields of random vectors, V2 sympy derivations of the
identities the explainer teaches, V7 invariance (rotation, reflection, scaling) and limits. ``np.dot``/``np.cross``
comparisons are V1 (an independent implementation of the same closed form), never counted as the only level.
"""
from __future__ import annotations

import mpmath as mp
import numpy as np
import pytest
import sympy as sp

from fluidpy.geo import vectors as V

RNG = np.random.default_rng(20261021)
A = RNG.normal(size=(400, 3))
B = RNG.normal(size=(400, 3))
C = RNG.normal(size=(400, 3))


def _rot(axis, ang):
    """Proper rotation matrix (Rodrigues) — built here, independent of the library."""
    k = np.asarray(axis, float) / np.linalg.norm(axis)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K


# ---------------------------------------------------------------------------------------------------------------- dot
def test_dot_V1_component_sum_on_field():  # V1
    d = V.dot(A, B)
    exact = A[:, 0] * B[:, 0] + A[:, 1] * B[:, 1] + A[:, 2] * B[:, 2]  # Kundu (K 2.2)
    assert d.shape == (400,)
    assert np.max(np.abs(d - exact) / (np.abs(exact) + 1)) < 1e-14
    assert V.dot([1, 2], [3, 4]) == 11.0 and isinstance(V.dot([1, 2], [3, 4]), float)


def test_dot_V1_cosine_form_on_field_of_angles():  # V1: a.b = |a||b| cos(alpha - beta) for planar vectors
    al = RNG.uniform(-np.pi, np.pi, 500)
    be = RNG.uniform(-np.pi, np.pi, 500)
    ra, rb = RNG.uniform(0.1, 10, 500), RNG.uniform(0.1, 10, 500)
    a = np.stack([ra * np.cos(al), ra * np.sin(al)], -1)
    b = np.stack([rb * np.cos(be), rb * np.sin(be)], -1)
    exact = ra * rb * np.cos(al - be)
    assert np.max(np.abs(V.dot(a, b) - exact) / (ra * rb)) < 1e-14


def test_dot_V2_symbolic_cosine_law():  # V2: from components to |a||b| cos(theta) (design D: a.b = |a||b| cos theta)
    ra, rb, al, be = sp.symbols("r_a r_b alpha beta", positive=True)
    a = sp.Matrix([ra * sp.cos(al), ra * sp.sin(al), 0])
    b = sp.Matrix([rb * sp.cos(be), rb * sp.sin(be), 0])
    assert sp.simplify(a.dot(b) - ra * rb * sp.cos(al - be)) == 0
    # law of cosines |a - b|^2 = |a|^2 + |b|^2 - 2 a.b (the explainer's derivation step)
    assert sp.simplify((a - b).dot(a - b) - (ra**2 + rb**2 - 2 * a.dot(b))) == 0


def test_dot_V7_rotation_invariance_symmetry_and_padding():  # V7
    R = _rot([1.0, 2.0, -0.5], 0.7)
    assert np.allclose(V.dot(A @ R.T, B @ R.T), V.dot(A, B), rtol=0, atol=1e-12)
    assert np.array_equal(V.dot(A, B), V.dot(B, A))
    assert V.dot([1.0, 2.0], [3.0, 4.0, 99.0]) == 11.0  # 2-vector padded with w = 0
    assert np.allclose(V.dot(A, [1.0, 0.0, 0.0]), A[:, 0])  # broadcasting
    assert V.dot([1, 0, 0], [0, 1, 0]) == 0.0


def test_dot_V2_dimensions_pint_advection_unit():  # V2 dimensional: the product unit of a velocity and a gradient
    from fluidpy.core.units import Q_
    u = Q_(np.array([10.0, 5.0, 0.0]), "m/s")
    gT = Q_(np.array([-2e-5, 1e-5, 0.0]), "K/m")
    adv = (u * gT).sum()  # what the PoC G "eq:1.7-adv" scratch cell does
    assert adv.check("[temperature]/[time]")
    assert np.isclose(V.dot(u.magnitude, gT.magnitude) * 86400.0, adv.to("K/day").magnitude, rtol=1e-14)


def test_vector_input_validation():
    with pytest.raises(ValueError):
        V.dot([1, 2, 3, 4], [1, 2, 3, 4])
    with pytest.raises(ValueError):
        V.norm(3.0)


# -------------------------------------------------------------------------------------------------------------- cross
def test_cross_V1_components_on_field():  # V1 (K 2.20) written out here, and np.cross
    c = V.cross(A, B)
    exact = np.stack([A[:, 1] * B[:, 2] - A[:, 2] * B[:, 1], A[:, 2] * B[:, 0] - A[:, 0] * B[:, 2],
                      A[:, 0] * B[:, 1] - A[:, 1] * B[:, 0]], -1)
    assert np.max(np.abs(c - exact)) < 1e-14
    assert np.max(np.abs(c - np.cross(A, B))) < 1e-14
    assert np.array_equal(V.cross([1, 0, 0], [0, 1, 0]), [0, 0, 1])  # right-handed


def test_cross_V2_symbolic_epsilon_and_lagrange():  # V2: eps_ijk form (K 2.21), Lagrange identity, perpendicularity
    a = sp.Matrix(sp.symbols("a1:4"))
    b = sp.Matrix(sp.symbols("b1:4"))
    eps = [[sum(sp.LeviCivita(i, j, k) * a[j] * b[k] for j in range(3) for k in range(3))] for i in range(3)]
    assert sp.simplify(sp.Matrix(eps) - a.cross(b)) == sp.zeros(3, 1)
    axb = a.cross(b)
    assert sp.expand(axb.dot(axb) + a.dot(b) ** 2 - a.dot(a) * b.dot(b)) == 0
    assert sp.expand(axb.dot(a)) == 0 and sp.expand(axb.dot(b)) == 0
    # the library obeys the same identities numerically on the field
    c = V.cross(A, B)
    lag = np.sum(c * c, -1) + V.dot(A, B) ** 2 - V.norm(A) ** 2 * V.norm(B) ** 2
    assert np.max(np.abs(lag) / (V.norm(A) ** 2 * V.norm(B) ** 2)) < 1e-13


def test_cross_V7_drawing_invariant_perpendicular_and_antisymmetric():  # V7 + drawing invariant |r.a| < 1e-9
    c = V.cross(A, B)
    assert np.max(np.abs(V.dot(c, A))) < 1e-12 and np.max(np.abs(V.dot(c, B))) < 1e-12
    assert np.array_equal(V.cross(A, B), -V.cross(B, A))
    assert np.allclose(V.cross(A, A), 0.0, atol=0)


def test_cross_V7_rotation_covariance_and_reflection_pseudovector():  # V7: R(a x b) = (Ra) x (Rb); mirror gives -
    R = _rot([0.3, -1.0, 2.0], 1.1)
    assert np.allclose(V.cross(A @ R.T, B @ R.T), V.cross(A, B) @ R.T, rtol=0, atol=1e-12)
    M = np.diag([1.0, 1.0, -1.0])  # mirror in the horizontal plane (det -1)
    assert np.allclose(V.cross(A @ M, B @ M), -(V.cross(A, B) @ M), rtol=0, atol=1e-12)


def test_cross_V1_horizontal_vectors_give_vertical():  # V1 2-vector padding
    c = V.cross([2.0, 1.0], [-1.0, 3.0])
    assert np.array_equal(c, [0.0, 0.0, 7.0])


# -------------------------------------------------------------------------------------------------------- norm, angle
def test_norm_V1_and_V7_scaling():  # V1 + V7
    assert V.norm([3, 4]) == 5.0 and isinstance(V.norm([3, 4]), float)
    assert np.allclose(V.norm(A), np.sqrt(np.sum(A * A, -1)), rtol=1e-15)
    assert np.allclose(V.norm(-3.5 * A), 3.5 * V.norm(A), rtol=1e-15)


def test_angle_between_V1_known_angles_and_field():  # V1
    assert V.angle_between([1, 0, 0], [0, 1, 0]) == pytest.approx(np.pi / 2, abs=1e-15)
    assert V.angle_between([1, 0], [1, 1], degrees=True) == pytest.approx(45.0, abs=1e-13)
    assert V.angle_between([1, 2, 3], [-1, -2, -3]) == pytest.approx(np.pi, abs=1e-15)
    assert V.angle_between([1, 2, 3], [2, 4, 6]) == 0.0
    th = V.angle_between(A, B)
    ref = np.arccos(np.clip(np.sum(A * B, -1) / np.linalg.norm(A, axis=1) / np.linalg.norm(B, axis=1), -1, 1))
    well = np.abs(np.cos(ref)) < 0.99  # arccos is well-conditioned here
    assert np.max(np.abs(th - ref)[well]) < 1e-13
    assert np.all((th >= 0) & (th <= np.pi))


def test_angle_between_V1_near_parallel_against_mpmath():  # V1 at 50-digit precision (design Q5: atan2 accuracy)
    mp.mp.dps = 50
    for eps in (1e-3, 1e-7, 1e-10, 1e-13):
        for sgn in (1.0, -1.0):
            a, b = [1.0, 0.0, 0.0], [sgn, eps, 0.0]
            exact = float(mp.atan2(mp.mpf(eps), mp.mpf(sgn)))
            got = V.angle_between(a, b)
            assert abs(got - exact) <= 4e-16 * max(abs(exact), 1e-300) + 1e-300, (eps, sgn, got, exact)


def test_angle_between_V7_invariance_and_zero_vector():  # V7 rotation and scale invariance; no direction -> NaN
    R = _rot([1, 1, 1], 2.0)
    assert np.allclose(V.angle_between(A @ R.T, 7.0 * B @ R.T), V.angle_between(A, B), rtol=0, atol=1e-12)
    assert np.isnan(V.angle_between([0, 0], [1, 0]))
    assert np.allclose(np.cos(V.angle_between(A, B)), V.dot(A, B) / (V.norm(A) * V.norm(B)), rtol=0, atol=1e-13)


# --------------------------------------------------------------------------------------------------------- projection
def test_projection_V1_remainder_perpendicular_and_idempotent():  # V1
    p = V.projection(A, B)
    assert np.max(np.abs(V.dot(A - p, B)) / (V.norm(A) * V.norm(B))) < 1e-14
    assert np.allclose(V.projection(p, B), p, rtol=0, atol=1e-13)  # projecting twice changes nothing
    assert np.allclose(V.projection(B, B), B, rtol=1e-15)
    assert np.allclose(V.projection(A, -4.0 * B), p, rtol=0, atol=1e-13)  # only the direction of b matters
    assert np.all(np.isnan(V.projection([1.0, 2.0], [0.0, 0.0])))


def test_scalar_projection_V1_and_V7():  # V1 |a| cos(theta); V7 independent of |b|, sign flips with b
    assert V.scalar_projection([3, 4], [5, 0]) == 3.0
    sproj = V.scalar_projection(A, B)
    assert np.allclose(sproj, V.norm(A) * np.cos(V.angle_between(A, B)), rtol=0, atol=1e-12)
    assert np.allclose(V.scalar_projection(A, 0.01 * B), sproj, rtol=1e-14)
    assert np.allclose(V.scalar_projection(A, -B), -sproj, rtol=1e-15)
    assert np.isnan(V.scalar_projection([1.0, 0.0], [0.0, 0.0]))
    p = V.projection(A, B)
    assert np.allclose(V.norm(p), np.abs(sproj), rtol=1e-13)  # the two projections agree


# --------------------------------------------------------------------------------------------- area and triple product
def test_parallelogram_area_V1_base_times_height():  # V1 (0A.2): |a||b| sin theta
    th = RNG.uniform(0, np.pi, 300)
    L1, L2 = RNG.uniform(0.5, 3, 300), RNG.uniform(0.5, 3, 300)
    a = np.stack([L1, 0 * L1], -1)
    b = np.stack([L2 * np.cos(th), L2 * np.sin(th)], -1)
    assert np.max(np.abs(V.parallelogram_area(a, b) - L1 * L2 * np.sin(th))) < 1e-14
    assert V.parallelogram_area([1, 0], [0, 1]) == 1.0


def test_parallelogram_area_V7_shear_invariance_and_parallel_limit():  # V7 (Cavalieri: area unchanged by shear)
    lam = RNG.normal(size=(400, 1))
    assert np.allclose(V.parallelogram_area(A, B + lam * A), V.parallelogram_area(A, B), rtol=1e-11)
    assert np.allclose(V.parallelogram_area(A, 2.5 * A), 0.0, atol=1e-13)


def test_triple_scalar_V1_determinant_and_V7_cyclic():  # V1 det; V7 cyclic/antisymmetric; coplanar -> 0
    det = np.linalg.det(np.stack([A, B, C], axis=1))
    assert np.allclose(V.triple_scalar(A, B, C), det, rtol=1e-12, atol=1e-13)
    assert np.allclose(V.triple_scalar(A, B, C), V.triple_scalar(B, C, A), rtol=1e-12, atol=1e-13)
    assert np.allclose(V.triple_scalar(A, B, C), -V.triple_scalar(B, A, C), rtol=1e-12, atol=1e-13)
    assert V.triple_scalar([1, 0, 0], [0, 1, 0], [0, 0, 1]) == 1.0
    assert np.allclose(V.triple_scalar(A, B, 2 * A - 3 * B), 0.0, atol=1e-12)
    assert V.triple_scalar([1, 2], [3, 4], [5, 6]) == 0.0


# -------------------------------------------------------------------------------------------------------------- k x u
def test_k_cross_V1_equals_cross_with_vertical():  # V1 Vallis (2.79): k x (u, v) = (-v, u)
    U2 = A[:, :2]
    assert np.array_equal(V.k_cross(U2), V.cross([0.0, 0.0, 1.0], U2)[:, :2])
    assert np.array_equal(V.k_cross(A), V.cross([0.0, 0.0, 1.0], A))
    assert np.array_equal(V.k_cross([1.0, 0.0]), [0.0, 1.0])  # east -> north: 90 deg to the left


def test_k_cross_V7_quarter_turn_properties():  # V7: k x (k x u) = -u, |k x u| = |u|, (k x u).u = 0
    U2 = A[:, :2]
    assert np.array_equal(V.k_cross(V.k_cross(U2)), -U2)
    assert np.allclose(V.norm(V.k_cross(U2)), V.norm(U2), rtol=1e-15)
    assert np.max(np.abs(V.dot(V.k_cross(U2), U2))) < 1e-15
    assert np.all(V.cross(U2, V.k_cross(U2))[:, 2] > 0)  # anticlockwise quarter turn seen from above
