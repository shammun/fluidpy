"""Verification suite for Chapter 3 — Kinematics (Kundu, Cohen & Dowling 5e, §§3.1–3.6).

Evidence levels (``verify-implementation`` skill): V1 analytic · V2 symbolic (sympy / dimensions) · V3 convergence ·
V4 conservation / invariant · V5 published benchmark (``reference/ch03/``, see SOURCES.md) · V6 book value (private,
``tests/book_values_ch03.json``, skipped when absent) · V7 limits / symmetry / invariance. Every test name is
``test_<concept>_<level>_<what>``; the comment on the ``def`` line repeats the level.

A items (CORE, ≥ 2 independent levels): C01 Lagrangian/Eulerian bridge (3.2) · C02 material derivative (3.4)–(3.5) ·
C03 streamlines (3.7) · C04 path and streak lines (3.8), Ex. 3.1 · C05 Galilean invariance (3.9) · C06 relative
velocity (3.10) · C07 linear strain rate · C08 shear strain rate · C09 volumetric strain rate (3.14) · C10 spin = ½ω ·
C11 deformation + rotation (3.19) · C12 principal strain axes (sphere → ellipsoid) · C13 polar vorticity (3.23) ·
C14 Rankine (3.28) and Gaussian (3.29) vortices · C15 Reynolds transport theorem (3.35).
Derivations re-derived with sympy (every ★★ and the ★★★ D22, plus the cheap ★ ones): D01–D24, ``test_*_V2_derivation``.

Pinned conventions with discrimination tests: G[i, j] = ∂u_i/∂x_j · book R = G − Gᵀ (no ½), ω = vector(R), spin = ½ω ·
γ = du₁/dx₂ = 2S₁₂ (ch02's Γ ≡ S₁₂ differs by 2) · ω₃ = −γ is clockwise · ω′ = ω − 2Ω in a rotating frame (not ω − Ω)
· Galilean map x = x′ + Ut + x′_o (sign of U) · (3.6) carries |u| ∂F/∂s (the book prints the F-less operator) ·
Leibniz's lower-limit term is subtracted · the RTT sliver is signed · Ex. 3.2: b·n = 0 on the base although b ≠ 0.

Run: ``.venv/Scripts/python.exe -m pytest tests/test_ch03.py -q -p no:cacheprovider``.
"""
from __future__ import annotations

import inspect
import json
import math
from pathlib import Path

import numpy as np
import pytest
import sympy as sp
from scipy.linalg import expm
from scipy.optimize import minimize_scalar
from scipy.spatial.transform import Rotation

from fluidpy import ch02_cartesian_tensors as ch02
from fluidpy import ch03_kinematics as ch03
from fluidpy.core import coords as Cmod
from fluidpy.core import kinematics as K
from fluidpy.core import transport as Rmod
from fluidpy.core import vortices as Xmod
from fluidpy.core.units import Q_
from tools.convergence import observed_order, pairwise_orders

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference" / "ch03"
BOOK = Path(__file__).resolve().parent / "book_values_ch03.json"
needs_ref = pytest.mark.skipif(not (REF / "benchmarks.json").exists(), reason="run reference/ch03/make_refs.py")
book_only = pytest.mark.skipif(not BOOK.exists(), reason="book values are private; see CLAUDE.md rule 9")

RNG = np.random.default_rng(3)
ORDER_TOL = 0.15  # design order ± this (verify-implementation default)
X_STAR = 1.2564312086261697  # root of 1 + 2x = e^x (x > 0), Lambert-W closed form (checked below)


def book():
    return json.loads(BOOK.read_text(encoding="utf-8"))


def ref_json(name="benchmarks.json"):
    return json.loads((REF / name).read_text(encoding="utf-8"))


def rand_rotation(n=3, rng=RNG):
    if n == 2:
        th = rng.uniform(0, 2 * np.pi)
        return np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    return Rotation.random(random_state=int(rng.integers(1_000_000))).as_matrix()


# ---------------------------------------------------------------------------------------------------------------------
# A symbolic test field: exact velocity, gradient, local/advective acceleration, vorticity — independent of the stencils
# ---------------------------------------------------------------------------------------------------------------------
X1, X2, X3, TT = sp.symbols("x1 x2 x3 t", real=True)


class SymField:
    """u(x, t) from sympy expressions, callable in the kinematics convention u(x (d,) or (d, N), t)."""

    def __init__(self, exprs, coords):
        self.exprs = [sp.sympify(e) for e in exprs]
        self.coords = tuple(coords)
        self.d = len(coords)
        args = (*self.coords, TT)
        self._u = sp.lambdify(args, self.exprs, "numpy")
        self.G_expr = [[sp.diff(e, c) for c in self.coords] for e in self.exprs]
        self._G = sp.lambdify(args, self.G_expr, "numpy")
        loc = [sp.diff(e, TT) for e in self.exprs]
        adv = [sum(u * sp.diff(e, c) for u, c in zip(self.exprs, self.coords)) for e in self.exprs]
        self._loc = sp.lambdify(args, loc, "numpy")
        self._adv = sp.lambdify(args, adv, "numpy")
        if self.d == 3:
            u1, u2, u3 = self.exprs
            curl = [sp.diff(u3, X2) - sp.diff(u2, X3), sp.diff(u1, X3) - sp.diff(u3, X1), sp.diff(u2, X1) - sp.diff(u1, X2)]
        else:
            u1, u2 = self.exprs
            curl = [0, 0, sp.diff(u2, X1) - sp.diff(u1, X2)]
        self._curl = sp.lambdify(args, curl, "numpy")
        self._div = sp.lambdify(args, sum(sp.diff(e, c) for e, c in zip(self.exprs, self.coords)), "numpy")

    @staticmethod
    def _stack(vals, shape):
        return np.stack([np.broadcast_to(np.asarray(v, dtype=float), shape) for v in vals])

    def __call__(self, x, t=0.0):
        x_ = np.asarray(x, dtype=float)
        shape = np.broadcast(x_[0], np.asarray(t)).shape
        return self._stack(self._u(*x_, t), shape)

    def G(self, x, t=0.0):
        x_ = np.asarray(x, dtype=float)
        return np.array([[np.broadcast_to(np.asarray(v, float), x_[0].shape) for v in row] for row in self._G(*x_, t)])

    def local(self, x, t=0.0):
        x_ = np.asarray(x, dtype=float)
        return self._stack(self._loc(*x_, t), x_[0].shape)

    def advective(self, x, t=0.0):
        x_ = np.asarray(x, dtype=float)
        return self._stack(self._adv(*x_, t), x_[0].shape)

    def curl(self, x, t=0.0):
        x_ = np.asarray(x, dtype=float)
        return self._stack(self._curl(*x_, t), x_[0].shape)

    def div(self, x, t=0.0):
        x_ = np.asarray(x, dtype=float)
        return np.broadcast_to(np.asarray(self._div(*x_, t), float), x_[0].shape)


# unsteady, nonlinear, compressible test fields (2-D and 3-D) and a scalar field on them
U2 = SymField([sp.sin(X1) * sp.cos(X2) + sp.Rational(3, 10) * TT + sp.Rational(1, 5) * X1,
               -sp.cos(X1) * sp.sin(X2) + sp.Rational(1, 5) * X1 * TT ** 2 + sp.Rational(1, 10) * X2 ** 2], (X1, X2))
U3 = SymField([sp.sin(X2) + sp.Rational(1, 5) * TT * X3 + sp.Rational(1, 10) * X1 ** 2,
               X1 * sp.cos(X3) - sp.Rational(1, 10) * TT + sp.Rational(3, 10) * X2,
               X1 * X2 + sp.sin(TT) - sp.Rational(1, 5) * X3 * X1], (X1, X2, X3))
F2_expr = sp.sin(X1) * sp.exp(sp.Rational(3, 10) * X2) * sp.cos(TT) + X1 * X2 * TT
F2 = sp.lambdify((X1, X2, TT), F2_expr, "numpy")
DF2_expr = sp.diff(F2_expr, TT) + U2.exprs[0] * sp.diff(F2_expr, X1) + U2.exprs[1] * sp.diff(F2_expr, X2)
DF2 = sp.lambdify((X1, X2, TT), DF2_expr, "numpy")
F2_local = sp.lambdify((X1, X2, TT), sp.diff(F2_expr, TT), "numpy")


def F2_field(x, t):
    x_ = np.asarray(x, dtype=float)
    return np.asarray(F2(x_[0], x_[1], t), dtype=float) + 0.0 * x_[0]


def rand_points(d, n, lo=-1.0, hi=1.0, rng=RNG):
    return rng.uniform(lo, hi, size=(d, n))


# =====================================================================================================================
# C01 — Lagrangian and Eulerian descriptions; the bridge (3.2); N07, N08 (3.1); D01
# =====================================================================================================================
def test_lagrangian_eulerian_V1_stretching_map_numbers():  # V1
    # the curation's worked number: X = 2 m, alpha = 0.5 1/s, t = 1 s
    X, al, t = 2.0, 0.5, 1.0
    assert abs(ch03.lagrangian_map_example(X, t, al) - 2.0 * math.exp(0.5)) < 1e-14  # 3.2974 m
    assert abs(ch03.lagrangian_map_example(X, t, al, derivative=1) - 1.0 * math.exp(0.5)) < 1e-14  # 1.6487 m/s
    assert abs(ch03.lagrangian_map_example(X, t, al, derivative=2) - 0.5 * math.exp(0.5)) < 1e-14  # 0.8244 m/s²
    # (3.1) by stencils on the map equals the Eulerian u = alpha x, a = alpha² x on a field of labels and times
    for Xl in np.linspace(-3, 3, 7):
        for tt in (0.0, 0.7, 2.0):
            u, a = ch03.lagrangian_velocity_acceleration(lambda s: ch03.lagrangian_map_example(Xl, s, al), tt)
            x = ch03.lagrangian_map_example(Xl, tt, al)
            assert abs(u - al * x) < 1e-8 * max(1.0, abs(x))
            assert abs(a - al ** 2 * x) < 1e-6 * max(1.0, abs(x))


def test_lagrangian_eulerian_V2_three_maps_Du_Dt_equals_d2r_dt2():  # V2
    x, y, t, X, Y = sp.symbols("x y t X Y", real=True)
    al, w, g = sp.symbols("alpha omega gamma", positive=True)
    cases = [  # (map, labels, coords, expected Eulerian u)
        ([X * sp.exp(al * t)], [X], [x], [al * x]),
        ([X * sp.cos(w * t) - Y * sp.sin(w * t), X * sp.sin(w * t) + Y * sp.cos(w * t)], [X, Y], [x, y], [-w * y, w * x]),
        ([X + g * Y * t, Y], [X, Y], [x, y], [g * y, 0]),
    ]
    for r, labels, coords, u_exp in cases:
        out = ch03.lagrangian_to_eulerian(r, labels, coords, t)
        assert all(sp.simplify(ui - ue) == 0 for ui, ue in zip(out["u"], u_exp))  # (3.2) bridge
        assert all(sp.simplify(d - a) == 0 for d, a in zip(out["Du_Dt"], out["a_lagrangian"]))  # D01 check
        assert out["a"] == out["a_lagrangian"]
        # the label substitution really inverts the map
        assert all(sp.simplify(ri.subs(out["label_of_x"]) - c) == 0 for ri, c in zip(r, coords))
    us, as_ = ch03.lagrangian_velocity_acceleration_sym([X * sp.exp(al * t)], t)
    assert sp.simplify(us[0] - al * X * sp.exp(al * t)) == 0 and sp.simplify(as_[0] - al ** 2 * X * sp.exp(al * t)) == 0


def test_lagrangian_eulerian_V7_compatibility_along_rotation_map():  # V7
    # (3.2): the Eulerian field evaluated on the particle equals the particle's own (Lagrangian) velocity
    w = 0.8
    field = lambda x, t: np.stack([-w * np.asarray(x)[1], w * np.asarray(x)[0]])  # noqa: E731  Eulerian twin
    for X0, Y0 in RNG.uniform(-2, 2, size=(10, 2)):
        r = lambda s: np.array([X0 * np.cos(w * s) - Y0 * np.sin(w * s), X0 * np.sin(w * s) + Y0 * np.cos(w * s)])  # noqa: E731
        for tt in (0.0, 0.9, 3.1):
            u_lag, a_lag = ch03.lagrangian_velocity_acceleration(r, tt)
            assert np.max(np.abs(u_lag - field(r(tt), tt))) < 1e-8
            assert np.max(np.abs(a_lag + w ** 2 * r(tt))) < 1e-6  # centripetal


def test_lagrangian_velocity_V3_order_two():  # V3
    X, al, tt = 2.0, 0.5, 1.0
    hs = [0.2, 0.1, 0.05, 0.025]
    eu, ea = [], []
    for h in hs:
        u, a = ch03.lagrangian_velocity_acceleration(lambda s: X * np.exp(al * s), tt, h=h)
        eu.append(abs(u - al * X * np.exp(al * tt)))
        ea.append(abs(a - al ** 2 * X * np.exp(al * tt)))
    assert abs(observed_order(hs, eu) - 2.0) < ORDER_TOL
    assert abs(observed_order(hs, ea) - 2.0) < ORDER_TOL


def test_lagrangian_map_V2_derivation():  # V2  D01 steps 1-6
    X, al, t, x = sp.symbols("X alpha t x", positive=True)
    path = X * sp.exp(al * t)
    u_lab = sp.diff(path, t)  # step 1
    a_lab = sp.diff(path, t, 2)  # step 2
    X_of_x = sp.solve(sp.Eq(x, path), X)
    assert len(X_of_x) == 1 and sp.simplify(X_of_x[0] - x * sp.exp(-al * t)) == 0  # step 3 (invertible)
    u_x = sp.simplify(u_lab.subs(X, X_of_x[0]))  # steps 4-5
    a_x = sp.simplify(a_lab.subs(X, X_of_x[0]))  # step 6
    assert sp.simplify(u_x - al * x) == 0 and sp.simplify(a_x - al ** 2 * x) == 0
    assert sp.simplify(sp.diff(u_x, t)) == 0  # steady field
    assert sp.simplify(sp.diff(u_x, t) + u_x * sp.diff(u_x, x) - a_x) == 0  # forward check: Du/Dt = a


# =====================================================================================================================
# C02 — the material derivative (3.3)–(3.5); N09–N12; D02
# =====================================================================================================================
def test_material_derivative_V1_equals_rate_along_pathline():  # V1 (independent route: ODE + difference along path)
    # D/Dt by stencils must equal d/dt F(r(t), t) measured along a path line integrated with solve_ivp
    for x0 in rand_points(2, 6).T:
        t0 = float(RNG.uniform(0, 1))
        dd = 1e-3
        pts = ch03.pathline(U2, x0, t0, [t0 - 2 * dd, t0 - dd, t0 + dd, t0 + 2 * dd], rtol=1e-12, atol=1e-14)
        Fv = [float(F2_field(pts[:, k], t0 + s * dd)) for k, s in enumerate((-2, -1, 1, 2))]
        rate_path = (Fv[0] - 8 * Fv[1] + 8 * Fv[2] - Fv[3]) / (12 * dd)  # 4th-order difference along the path
        rate_stencil = ch03.material_derivative(F2_field, U2, x0, t0)
        assert abs(rate_stencil - rate_path) < 1e-7
        assert abs(rate_stencil - float(DF2(x0[0], x0[1], t0))) < 1e-7  # and the exact sympy value


def test_material_derivative_V2_sym_equals_stencils_on_field():  # V2
    F = ch03.material_derivative_sym(F2_expr, U2.exprs, [X1, X2], TT)
    assert sp.simplify(F - DF2_expr) == 0
    Fn = sp.lambdify((X1, X2, TT), F, "numpy")
    P = rand_points(2, 200)
    loc, adv, tot = ch03.material_derivative_terms(F2_field, U2, P, 0.4)
    assert np.max(np.abs(tot - Fn(P[0], P[1], 0.4))) < 1e-7
    assert np.max(np.abs(loc - F2_local(P[0], P[1], 0.4))) < 1e-7
    assert np.max(np.abs(loc + adv - tot)) < 1e-15
    # vector F (component-wise) = the acceleration of the field
    acc = ch03.material_derivative_sym(U2.exprs, U2.exprs, [X1, X2], TT)
    An = sp.lambdify((X1, X2, TT), acc, "numpy")
    assert np.max(np.abs(np.stack(An(P[0], P[1], 0.4)) - (U2.local(P, 0.4) + U2.advective(P, 0.4)))) < 1e-12


def test_material_derivative_V3_stencil_order_two():  # V3
    x0, t0 = np.array([0.3, -0.4]), 0.6
    hs = [0.1, 0.05, 0.025, 0.0125]
    exact = float(DF2(x0[0], x0[1], t0))
    errs = [abs(ch03.material_derivative(F2_field, U2, x0, t0, h=h) - exact) for h in hs]
    assert abs(observed_order(hs, errs) - 2.0) < ORDER_TOL, pairwise_orders(hs, errs)


def test_material_derivative_V7_local_and_advective_vanish():  # V7
    P = rand_points(2, 50)
    steady_F = lambda x, t: np.sin(np.asarray(x)[0]) * np.cos(np.asarray(x)[1])  # noqa: E731
    loc, adv, tot = ch03.material_derivative_terms(steady_F, U2, P, 0.3)
    assert np.max(np.abs(loc)) < 1e-12 and np.max(np.abs(adv)) > 1e-2  # steady F: all advective
    zero_u = lambda x, t: np.zeros_like(np.asarray(x, float))  # noqa: E731
    loc, adv, tot = ch03.material_derivative_terms(F2_field, zero_u, P, 0.3)
    assert np.max(np.abs(adv)) == 0.0  # u = 0: all local
    # u ⟂ ∇F: F = F(y), u = (1, 0)
    Fy = lambda x, t: np.asarray(x)[1] ** 2 * (1 + t)  # noqa: E731
    ux = lambda x, t: np.stack([np.ones_like(np.asarray(x)[0]), np.zeros_like(np.asarray(x)[0])])  # noqa: E731
    loc, adv, tot = ch03.material_derivative_terms(Fy, ux, P, 0.3)
    assert np.max(np.abs(adv)) < 1e-15


def test_thermal_advection_V1_worked_number_and_stencils():  # V1  N10 (climate hook)
    # cold to the north (1 K per 100 km), southerly wind 10 m/s, pattern carried by the wind: DT/Dt = 0
    out = ch03.thermal_front_terms(0.0, 0.0, 0.0, 0.0, 10.0, 1e-5, 0.0, 10.0)
    assert out["advective"] == pytest.approx(-1e-4, rel=1e-12)
    assert out["local"] == pytest.approx(1e-4, rel=1e-12)
    assert abs(out["total"]) < 1e-18
    assert out["local"] * 3600 == pytest.approx(0.36, rel=1e-12)  # K/h at the station
    assert out["regime"] == "warm advection"
    assert ch03.thermal_front_terms(0, 0, 0, 0, -10.0, 1e-5)["regime"] == "cold advection"
    assert ch03.thermal_front_terms(0, 0, 0, 10.0, 0.0, 1e-5)["regime"] == "no advection"
    # exact split == stencil split on the linear and on the tanh front (heating and a moving pattern)
    for width in (None, 2e5):
        G, H, c, u, v = 1e-5, 1.0 / 3600, 5.0, 3.0, 10.0
        F = lambda x, t: ch03.thermal_front(np.asarray(x)[0], np.asarray(x)[1], t, G, H, c, width)  # noqa: E731
        wind = lambda x, t: np.stack([u + 0 * np.asarray(x)[0], v + 0 * np.asarray(x)[0]])  # noqa: E731
        for yy in (-3e5, 0.0, 1.5e5):
            ex = ch03.thermal_front_terms(0.0, yy, 7200.0, u, v, G, H, c, width)
            loc, adv, tot = ch03.material_derivative_terms(F, wind, np.array([0.0, yy]), 7200.0, h=10.0, ht=10.0)
            assert abs(loc - ex["local"]) < 1e-9 and abs(adv - ex["advective"]) < 1e-9
            assert abs(tot - ex["total"]) < 1e-9


def test_streamwise_derivative_V1_equals_advective_term():  # V1  N12 (3.6) corrected
    P = rand_points(2, 100)
    sw = ch03.streamwise_derivative(F2_field, U2, P, 0.2)
    _, adv, _ = ch03.material_derivative_terms(F2_field, U2, P, 0.2)
    assert np.max(np.abs(sw - adv)) < 1e-7
    with pytest.raises(ValueError):
        ch03.streamwise_derivative(F2_field, lambda x, t: np.zeros(2), np.array([0.1, 0.2]), 0.0)


def test_streamwise_derivative_V2_dimensions_expose_book_typo():  # V2  (3.6) as printed lacks F
    # both terms of DF/Dt must be [F]/s: with F = temperature, |u| dF/ds is K/s; the printed |u| d/ds alone is 1/s
    dFdt = Q_(1e-4, "K/s")
    speed, dFds, per_s = Q_(10.0, "m/s"), Q_(-1e-5, "K/m"), Q_(1.0, "1/m")
    corrected = speed * dFds
    printed = speed * per_s
    assert corrected.dimensionality == dFdt.dimensionality
    assert printed.dimensionality != dFdt.dimensionality


def test_material_derivative_V2_index_form():  # V2  N11: u_i dF/dx_i expands to three terms
    expr = ch03.expand_indices("u_i F,i")
    assert len(sp.Add.make_args(sp.expand(expr))) == 3
    assert all(isinstance(term, sp.Mul) for term in sp.Add.make_args(sp.expand(expr)))


def test_material_derivative_V2_derivation():  # V2  D02 steps 1-9 (chain rule along a path)
    t = sp.symbols("t", real=True)
    x, y, z = sp.symbols("x y z", real=True)
    F = sp.Function("F")
    r = [sp.Function(f"r{i}")(t) for i in (1, 2, 3)]
    f = F(*r, t)  # steps 1-2
    dfdt = sp.diff(f, t)  # sympy's chain rule
    xs = sp.symbols("s1 s2 s3 s4")
    G_ = F(*xs)
    chain = sum(sp.diff(G_, xs[i]).subs(dict(zip(xs, [*r, t]))) * sp.diff(r[i], t) for i in range(3)) \
        + sp.diff(G_, xs[3]).subs(dict(zip(xs, [*r, t])))  # step 3, (3.3)
    assert sp.simplify(dfdt - chain) == 0
    # steps 4-8 on a concrete flow whose paths are known: u = (a x, -b y, c), paths r = (X e^{at}, Y e^{-bt}, Z + c t)
    a, b, c, X, Y, Z = sp.symbols("a b c X Y Z", real=True)
    u = [a * x, -b * y, c]
    path = [X * sp.exp(a * t), Y * sp.exp(-b * t), Z + c * t]
    Fc = x ** 2 * y * sp.sin(t) + z * sp.exp(-t) * x
    along = sp.diff(Fc.subs({x: path[0], y: path[1], z: path[2]}), t)
    DFDt = ch03.material_derivative_sym(Fc, u, [x, y, z], t).subs({x: path[0], y: path[1], z: path[2]})
    assert sp.simplify(along - DFDt) == 0  # (3.5) on the particle
    # step 9: D/Dt of a vector acts component by component (the particle acceleration)
    acc = ch03.material_derivative_sym(u, u, [x, y, z], t)
    acc_path = [sp.diff(p, t, 2) for p in path]
    assert all(sp.simplify(ai.subs({x: path[0], y: path[1], z: path[2]}) - ap) == 0 for ai, ap in zip(acc, acc_path))


# =====================================================================================================================
# C03 — streamlines (3.7); N14, N15; D03
# =====================================================================================================================
def test_streamline_V1_ex31_line_and_solid_body_circle():  # V1
    for tp in (0.0, 0.3, 1.1):
        ex = ch03.example_3_1(tp)
        sl = ch03.streamline(ex["field"], [0.0, 0.0], tp, 2.0)
        assert sl.shape == (2, 400)
        # the line y = x tan(ωt'), written without dividing by cos (cross-product form)
        assert np.max(np.abs(sl[1] * np.cos(tp) - sl[0] * np.sin(tp))) < 1e-12
        s = np.hypot(sl[0], sl[1])
        assert abs(s.max() - 2.0) < 1e-9  # arc length s_max each way
    vortex = ch03.preset_field("steady_vortex", Omega=0.7)
    sl = ch03.streamline(vortex, [1.3, 0.0], 0.0, 3.0, rtol=1e-12, atol=1e-14)
    assert np.max(np.abs(np.hypot(sl[0], sl[1]) - 1.3)) < 1e-10  # circles of constant radius (error ∝ rtol)


def test_streamline_V4_streamfunction_constant_along_curve():  # V4
    # cylinder (body frame): psi = U y (1 - a²/r²)
    U, a = 1.0, 1.0
    fld = ch03.cylinder_velocity_field(U, a)
    for seed in ([-4.0, 0.5], [-4.0, 1.5], [3.0, -0.8]):
        sl = ch03.streamline(fld, seed, 0.0, 6.0)
        psi = ch03.cylinder_streamfunction(sl[0], sl[1], U, a)
        assert np.max(np.abs(psi - psi[0])) < 1e-8
    # stagnation-point flow u = (x, -y): psi = x y
    stag = lambda x, t: np.stack([np.asarray(x)[0], -np.asarray(x)[1]])  # noqa: E731
    sl = ch03.streamline(stag, [0.5, 2.0], 0.0, 3.0)
    assert np.max(np.abs(sl[0] * sl[1] - 1.0)) < 1e-8


def test_streamline_V7_stagnation_point_stops_curve():  # V7
    sink = lambda x, t: -np.asarray(x, dtype=float)  # noqa: E731  all streamlines end at the origin
    sl = ch03.streamline(sink, [1.0, 0.0], 0.0, 2.0, both=False, n=201)
    assert sl.shape[1] < 201 and np.all(np.hypot(sl[0], sl[1]) <= 1.0 + 1e-12)
    with pytest.raises(ValueError):
        ch03.streamline(sink, [0.0, 0.0], 0.0, 1.0)


def test_streamline_slope_V1_worked_number():  # V1  N14
    u = lambda x, t: np.array([1.0, 2.0])  # noqa: E731
    assert ch03.streamline_slope(u, [0.3, 0.1], 0.0) == 2.0
    ex = ch03.example_3_1(math.pi / 6)
    assert abs(ch03.streamline_slope(ex["field"], [0.0, 0.0], math.pi / 6) - math.tan(math.pi / 6)) < 1e-14
    assert math.isinf(ch03.streamline_slope(lambda x, t: np.array([0.0, 1.0]), [0, 0], 0.0))


def test_stream_tube_V4_equal_flux_through_two_sections():  # V4  N15
    # axisymmetric stagnation flow u = (-x/2, -y/2, z): div u = 0, stream surfaces R² z = const
    u = lambda X, t: np.stack([-0.5 * X[0], -0.5 * X[1], X[2]])  # noqa: E731
    R1, z1, z2 = 0.8, 1.0, 2.5
    R2 = R1 * math.sqrt(z1 / z2)
    q1 = ch03.flux_through_disc(u, (0, 0, z1), (0, 0, 1), R1)
    q2 = ch03.flux_through_disc(u, (0, 0, z2), (0, 0, 1), R2)
    assert q1 == pytest.approx(math.pi * R1 ** 2 * z1, rel=1e-12)
    assert abs(q1 - q2) < 1e-12 * q1
    # the tube wall is a stream surface: streamlines from the rim stay on R² z = const
    sl = ch03.streamline(lambda x, t: u(np.asarray(x).reshape(3, -1), t).reshape(np.shape(x)), [R1, 0.0, z1], 0.0, 1.5)
    assert np.max(np.abs(np.hypot(sl[0], sl[1]) ** 2 * sl[2] - R1 ** 2 * z1)) < 1e-9


def test_streamline_V2_derivation():  # V2  D03 steps 1-5
    lam, u, v, w, s = sp.symbols("lambda u v w s", positive=True)
    ds = sp.Matrix([lam * u, lam * v, lam * w])  # steps 1-2
    assert ds.cross(sp.Matrix([u, v, w])) == sp.zeros(3, 1)  # u x ds = 0
    assert sp.simplify(ds[0] / u - ds[1] / v) == 0 and sp.simplify(ds[1] / v - ds[2] / w) == 0  # step 3: (3.7)
    lam_s = s / sp.sqrt(u ** 2 + v ** 2 + w ** 2)  # step 4
    step = ds.subs(lam, lam_s)
    assert sp.simplify(step.norm() - s) == 0  # step 5: |dx| = ds, dx/ds = u/|u|
    # check: solid-body rotation u = -y, v = x: dy/dx = -x/y ⇒ x² + y² constant
    x = sp.symbols("x", real=True)
    y = sp.Function("y")(x)
    assert sp.simplify(sp.diff(x ** 2 + y ** 2, x).subs(sp.diff(y, x), -x / y)) == 0


# =====================================================================================================================
# C04 — path lines (3.8), streak lines, Example 3.1 (N13, N16–N19); D04, D05
# =====================================================================================================================
def test_pathline_V1_linear_flows_match_expm():  # V1
    for name in ch03.VELOCITY_GRADIENT_PRESETS:
        G = ch03.velocity_gradient_preset(name, 0.8)
        fld = lambda x, t, G=G: G @ np.asarray(x)  # noqa: E731
        r0 = np.array([0.4, -0.7])
        te = np.array([-1.0, -0.2, 0.5, 1.5])
        pl = ch03.pathline(fld, r0, 0.3, te)
        ex = np.stack([expm(G * (tt - 0.3)) @ r0 for tt in te], axis=1)
        assert np.max(np.abs(pl - ex)) < 1e-9 * max(1.0, np.max(np.abs(ex)))
    # 3-D random linear field
    G = RNG.normal(size=(3, 3)) * 0.5
    r0 = RNG.normal(size=3)
    pl = ch03.pathline(lambda x, t: G @ np.asarray(x), r0, 0.0, 1.0)
    assert pl.shape == (3,) and np.max(np.abs(pl - expm(G) @ r0)) < 1e-9 * np.max(np.abs(expm(G) @ r0))


def test_pathline_V1_ex31_circle():  # V1  D04 result
    for tp in (0.0, 0.6, 2.0):
        ex = ch03.example_3_1(tp, xi0=1.5, omega=0.8)
        te = np.linspace(tp - 3, tp + 9, 60)
        pl = ch03.pathline(ex["field"], [0.0, 0.0], tp, te)
        c = ex["path_center"]
        assert np.max(np.abs(np.hypot(pl[0] - c[0], pl[1] - c[1]) - 1.5)) < 1e-9
        # parametric form of D04 step 4
        px = 1.5 * (np.sin(0.8 * te) - np.sin(0.8 * tp))
        py = 1.5 * (np.cos(0.8 * tp) - np.cos(0.8 * te))
        assert np.max(np.abs(pl - np.stack([px, py]))) < 1e-9
        # the closed-form curve in example_3_1 lies on the same circle
        assert np.max(np.abs(np.hypot(ex["pathline"][0] - c[0], ex["pathline"][1] - c[1]) - 1.5)) < 1e-13


def test_pathline_V3_error_decreases_with_rtol():  # V3
    ex = ch03.example_3_1(0.4)
    te = np.array([7.0])
    exact = np.array([np.sin(7.0) - np.sin(0.4), np.cos(0.4) - np.cos(7.0)])
    errs = []
    for rtol in (1e-4, 1e-6, 1e-8, 1e-10):
        pl = ch03.pathline(ex["field"], [0.0, 0.0], 0.4, te, rtol=rtol, atol=rtol * 1e-2, method="RK45")
        errs.append(float(np.max(np.abs(pl[:, 0] - exact))))
    assert all(e2 < e1 for e1, e2 in zip(errs, errs[1:])), errs
    assert errs[-1] < 1e-8


def test_streakline_V1_ex31_mirror_circle():  # V1  D05 result
    for tp in (0.0, 0.5, 1.9):
        ex = ch03.example_3_1(tp)
        tr = np.linspace(tp - 2 * np.pi, tp, 181)
        sk = ch03.streakline(ex["field"], [0.0, 0.0], tp, tr)
        c = ex["streak_center"]
        assert np.max(np.abs(np.hypot(sk[0] - c[0], sk[1] - c[1]) - 1.0)) < 1e-7
        assert np.allclose(c, -ex["path_center"], atol=1e-15)  # reflected through the origin
        # parametric form in the release time (D05 step 2), column order = release order
        sx = np.sin(tp) - np.sin(tr)
        sy = np.cos(tr) - np.cos(tp)
        assert np.max(np.abs(sk - np.stack([sx, sy]))) < 1e-7
    with pytest.raises(ValueError):
        ch03.streakline(ex["field"], [0.0, 0.0], 0.0, [0.5])


def test_flow_lines_V7_coincide_in_steady_flow():  # V7  N13
    for fld, x0 in ((ch03.preset_field("steady_vortex", Omega=0.9), np.array([1.0, 0.0])),
                    (ch03.preset_field("steady", omega=1.0, xi0=1.0, beta=0.4), np.array([0.0, 0.0]))):
        t = 2.0
        tr = np.linspace(0.0, t, 41)
        sk = ch03.streakline(fld, x0, t, tr)  # dye released from x0
        pl = ch03.pathline(fld, x0, 0.0, t - tr)  # where the particle at x0 at t = 0 is after the same durations
        assert np.max(np.abs(sk - pl)) < 1e-7
        sl = ch03.streamline(fld, x0, 0.0, 3.0, both=False, n=300)
        # every streak/path point lies on the streamline through x0 (distance to the sampled curve)
        from scipy.spatial import cKDTree
        dist, _ = cKDTree(sl.T).query(sk.T)
        assert np.max(dist) < 3.0 / 299  # within one sample spacing
        if np.allclose(fld(np.array([2.0, 1.0]), 0.0), fld(np.array([0.0, 0.0]), 0.0)):
            # uniform flow: all three are the same straight line through x0
            d = fld(x0, 0.0) / np.linalg.norm(fld(x0, 0.0))
            assert np.max(np.abs(sk[0] * d[1] - sk[1] * d[0])) < 1e-9


def test_example_3_1_V1_closed_forms_and_tangency():  # V1  N18, D05 steps 7-8
    for ang in (0.0, 0.5, 1.2, 2.5):
        ex = ch03.example_3_1(ang, xi0=2.0)
        c_p, c_s, R = ex["path_center"], ex["streak_center"], ex["radius"]
        assert R == 2.0 and abs(np.hypot(*c_p) - R) < 1e-14 and abs(np.hypot(*c_s) - R) < 1e-14  # both pass the origin
        # slope of each circle at the origin: -(x - cx)/(y - cy) at (0, 0) = cx/cy... = tan(ωt′)
        sp_path = -(0 - c_p[0]) / (0 - c_p[1])
        sp_streak = -(0 - c_s[0]) / (0 - c_s[1])
        assert abs(sp_path - math.tan(ang)) < 1e-12 and abs(sp_streak - math.tan(ang)) < 1e-12
        assert abs(ex["slope"] - math.tan(ang)) < 1e-12 and ex["angle"] == ang
        np.testing.assert_allclose(ex["pathline_center"], c_p)
        np.testing.assert_allclose(ex["streakline_center"], c_s)
        assert ex["streamline"].shape == ex["pathline"].shape == ex["streakline"].shape == (2, 200)
    assert math.isinf(ch03.example_3_1(math.pi / 2)["slope"])


def test_unsteady_presets_V1_fields():  # V1  E1 fields
    t = np.linspace(0, 3, 7)
    u, v = ch03.unsteady_flow_preset("ex31", 0.3, -0.2, t, xi0=1.5, omega=0.8)
    assert np.allclose(u, 1.2 * np.cos(0.8 * t)) and np.allclose(v, 1.2 * np.sin(0.8 * t))
    u2, v2 = ch03.unsteady_flow_preset("ex31_current", 0.3, -0.2, t, xi0=1.5, omega=0.8, U0=0.4)
    assert np.allclose(u2 - u, 0.4) and np.allclose(v2, v)
    # rotating strain: traceless, irrotational, stretching axis at angle Ω t
    for tt in (0.0, 0.7, 2.2):
        fld = ch03.preset_field("rotating_strain", s=0.6, Omega=0.5)
        G = ch03.velocity_gradient_at(fld, np.array([0.2, 0.1]), tt)
        assert abs(np.trace(G)) < 1e-9 and abs(G[0, 1] - G[1, 0]) < 1e-9
        lam, C = ch03.principal_strain_rates(G)
        assert np.allclose(lam, [-0.6, 0.6], atol=1e-8)
        assert abs(abs(C[:, 1] @ np.array([np.cos(0.5 * tt), np.sin(0.5 * tt)])) - 1.0) < 1e-8
    # ex31 + current: path line is the circle drifting at U0 (a trochoid)
    fld = ch03.preset_field("ex31_current", xi0=1.0, omega=1.0, U0=0.4)
    te = np.linspace(0, 6, 13)
    pl = ch03.pathline(fld, [0.0, 0.0], 0.0, te)
    assert np.max(np.abs(pl[0] - (0.4 * te + np.sin(te)))) < 1e-9 and np.max(np.abs(pl[1] - (1 - np.cos(te)))) < 1e-9
    with pytest.raises(ValueError):
        ch03.unsteady_flow_preset("nope", 0, 0, 0)


def test_pathline_V2_derivation():  # V2  D04 steps 1-7
    t, tp, w, xi = sp.symbols("t tp omega xi0", positive=True)
    cx, cy = sp.symbols("c_x c_y")
    x = sp.integrate(w * xi * sp.cos(w * t), t) + cx  # step 2
    y = sp.integrate(w * xi * sp.sin(w * t), t) + cy
    sol = sp.solve([x.subs(t, tp), y.subs(t, tp)], [cx, cy], dict=True)[0]  # step 3
    assert sp.simplify(sol[cx] + xi * sp.sin(w * tp)) == 0 and sp.simplify(sol[cy] - xi * sp.cos(w * tp)) == 0
    x, y = x.subs(sol), y.subs(sol)  # step 4
    assert sp.simplify(x - xi * (sp.sin(w * t) - sp.sin(w * tp))) == 0
    assert sp.simplify(y - xi * (sp.cos(w * tp) - sp.cos(w * t))) == 0
    circle = (x + xi * sp.sin(w * tp)) ** 2 + (y - xi * sp.cos(w * tp)) ** 2 - xi ** 2  # steps 5-7
    assert sp.simplify(circle) == 0
    assert sp.simplify(sp.diff(x, t) - w * xi * sp.cos(w * t)) == 0  # satisfies (3.8)
    # counterclockwise: (r - c) x v > 0
    cxv, cyv = -xi * sp.sin(w * tp), xi * sp.cos(w * tp)
    cross = (x - cxv) * sp.diff(y, t) - (y - cyv) * sp.diff(x, t)
    assert sp.simplify(cross - w * xi ** 2) == 0


def test_streakline_V2_derivation():  # V2  D05 steps 1-8
    t, tp, to, w, xi = sp.symbols("t tp t_o omega xi0", positive=True)
    x = xi * (sp.sin(w * t) - sp.sin(w * to))  # step 1: released at the origin at t_o
    y = xi * (sp.cos(w * to) - sp.cos(w * t))
    assert x.subs(t, to) == 0 and y.subs(t, to) == 0
    xs, ys = x.subs(t, tp), y.subs(t, tp)  # step 2: snapshot at t'
    circle = (xs - xi * sp.sin(w * tp)) ** 2 + (ys + xi * sp.cos(w * tp)) ** 2 - xi ** 2  # steps 3-5
    assert sp.simplify(circle) == 0
    # step 6: the dye released at t_o sits at polar angle ω t_o (measured from +y, clockwise sense) about the centre,
    # so releases over one period (ω t_o sweeping 2π) cover the whole circle
    assert sp.simplify(xs - xi * sp.sin(w * tp) + xi * sp.sin(w * to)) == 0
    assert sp.simplify(ys + xi * sp.cos(w * tp) - xi * sp.cos(w * to)) == 0
    # steps 7-8: implicit slopes of both circles at the origin equal the streamline slope tan ωt'
    X, Y = sp.symbols("X Y", real=True)
    for cxv, cyv in ((-xi * sp.sin(w * tp), xi * sp.cos(w * tp)), (xi * sp.sin(w * tp), -xi * sp.cos(w * tp))):
        Fc = (X - cxv) ** 2 + (Y - cyv) ** 2
        slope = (-sp.diff(Fc, X) / sp.diff(Fc, Y)).subs({X: 0, Y: 0})
        assert sp.simplify(slope - sp.tan(w * tp)) == 0


# =====================================================================================================================
# C05 — Galilean invariance of the acceleration (3.9); N02, N04, N20–N24; D06
# =====================================================================================================================
def test_galilean_V1_acceleration_equal_in_both_frames():  # V1
    U = np.array([0.7, -0.4])
    x0p = np.array([0.1, 0.2])
    up = ch03.galilean_transform(U2, U, x0p)
    P = rand_points(2, 200)
    tp = 0.35
    X = P + U[:, None] * tp + x0p[:, None]  # x = x' + U t + x'_o
    assert np.max(np.abs(up(P, tp) + U[:, None] - U2(X, tp))) < 1e-14  # u = U + u'
    a_p = ch03.acceleration(up, P, tp)
    a_l = ch03.acceleration(U2, X, tp)
    assert np.max(np.abs(a_p.a - a_l.a)) < 1e-7  # (3.9): the sum agrees
    assert np.max(np.abs(a_l.a - (U2.local(X, tp) + U2.advective(X, tp)))) < 1e-7  # and equals the exact value
    assert np.max(np.abs(a_p.local - a_l.local)) > 0.1  # the split does not
    # discrimination: the wrong-sign map x = x' - U t (a common slip) breaks the equality
    wrong = lambda xp, t: U2(np.asarray(xp) - U[:, None] * t + x0p[:, None], t) - U[:, None]  # noqa: E731
    assert np.max(np.abs(ch03.acceleration(wrong, P, tp).a - a_l.a)) > 1e-2
    # and a rotating observer (u - Ω x x) is not Galilean: its (u·∇)u + ∂u/∂t differs
    Om = 0.5
    rot = lambda x, t: U2(x, t) - Om * np.stack([-np.asarray(x)[1], np.asarray(x)[0]])  # noqa: E731
    assert np.max(np.abs(ch03.acceleration(rot, X, tp).a - a_l.a)) > 1e-2


def test_galilean_V1_sine_wave_worked_number():  # V1  C05 worked number
    # lab frame: u = 2 + 0.5 sin(x - 2t) (a wave carried at 2 m/s); wave frame: u' = 0.5 sin x'
    lab = lambda x, t: np.stack([2 + 0.5 * np.sin(np.asarray(x)[0] - 2 * t), 0 * np.asarray(x)[0]])  # noqa: E731
    P = np.array([math.pi / 4, 0.0])
    a = ch03.acceleration(lab, P, 0.0)
    assert a.local[0] == pytest.approx(-math.cos(math.pi / 4), abs=1e-8)  # -0.707
    assert a.advective[0] == pytest.approx((2 + 0.5 * math.sin(math.pi / 4)) * 0.5 * math.cos(math.pi / 4), abs=1e-8)
    assert a.a[0] == pytest.approx(0.125, abs=1e-8)
    wave = ch03.galilean_transform(lab, [2.0, 0.0])
    aw = ch03.acceleration(wave, P, 0.0)
    assert abs(aw.local[0]) < 1e-9 and aw.advective[0] == pytest.approx(0.125, abs=1e-8)


def test_galilean_V2_derivation():  # V2  D06 (the design's sympy check, general primed field)
    x, y, t, U1, U2_ = sp.symbols("x y t U_1 U_2")
    xp, yp, tp = sp.symbols("xp yp tp")
    f, g = sp.Function("f"), sp.Function("g")
    upr = [f(xp, yp, tp), g(xp, yp, tp)]
    sub = {xp: x - U1 * t, yp: y - U2_ * t, tp: t}
    u = [U1 + upr[0].subs(sub), U2_ + upr[1].subs(sub)]
    a_lab = [sp.diff(q, t) + u[0] * sp.diff(q, x) + u[1] * sp.diff(q, y) for q in u]
    a_pr = [(sp.diff(q, tp) + upr[0] * sp.diff(q, xp) + upr[1] * sp.diff(q, yp)).subs(sub) for q in upr]
    # NOTE (reported to the designer): the design's check_src stops here with sp.simplify(...doit()) and, under
    # sympy 1.14, prints two different Subs forms of the same ∂f/∂t' instead of [0, 0]. Closing the identity for a
    # generic cubic primed field (20 free coefficients per component) removes the representation ambiguity.
    a_, b_, c_ = sp.symbols("a_ b_ c_")
    mons = [a_ ** i * b_ ** j * c_ ** k for i in range(4) for j in range(4) for k in range(4) if i + j + k <= 3]
    Pf = sum(sp.Symbol(f"p{n}") * m for n, m in enumerate(mons))
    Pg = sum(sp.Symbol(f"q{n}") * m for n, m in enumerate(mons))
    repl = {f: sp.Lambda((a_, b_, c_), Pf), g: sp.Lambda((a_, b_, c_), Pg)}
    assert [sp.expand((a - b).subs(repl).doit()) for a, b in zip(a_lab, a_pr)] == [0, 0]
    # step 3 (the trap): the local terms differ by -U·∇'u'
    loc_lab = sp.diff(u[0], t)
    loc_pr = sp.diff(upr[0], tp).subs(sub)
    extra = -(U1 * sp.diff(upr[0], xp) + U2_ * sp.diff(upr[0], yp)).subs(sub)
    assert sp.expand((loc_lab - loc_pr - extra).subs(repl).doit()) == 0
    assert sp.expand((loc_lab - loc_pr).subs(repl).doit()) != 0  # the classic slip ∂u/∂t = ∂u'/∂t' is wrong
    # step 4: space derivatives unchanged
    assert sp.expand((sp.diff(u[0], x) - sp.diff(upr[0], xp).subs(sub)).subs(repl).doit()) == 0


def test_nonlinearity_V1_local_linear_advective_quadratic():  # V1  N22
    P = rand_points(2, 30)
    base = ch03.acceleration(U2, P, 0.2)
    lam = 2.0
    scaled = ch03.acceleration(lambda x, t: lam * U2(x, t), P, 0.2)
    assert np.max(np.abs(scaled.local - lam * base.local)) < 1e-8
    assert np.max(np.abs(scaled.advective - lam ** 2 * base.advective)) < 1e-7


def test_cylinder_flow_V1_boundary_far_field_and_frames():  # V1  N04, N20
    U, a = 1.3, 0.7
    th = np.linspace(0, 2 * np.pi, 200, endpoint=False)
    r = a * (1 + 1e-12)
    u, v = ch03.cylinder_flow(r * np.cos(th), r * np.sin(th), U, a)
    assert np.max(np.abs(u * np.cos(th) + v * np.sin(th))) < 1e-10  # no flow through the body
    uf, vf = ch03.cylinder_flow(1e4, 3e3, U, a)
    assert abs(uf - U) < 1e-7 and abs(vf) < 1e-7
    P = rand_points(2, 100, -3, 3)
    ub, vb = ch03.cylinder_flow(P[0], P[1], U, a, "body")
    ufl, vfl = ch03.cylinder_flow(P[0], P[1], U, a, "fluid")
    ok = np.isfinite(ub)
    assert np.max(np.abs(ub[ok] - ufl[ok] - U)) < 1e-14 and np.max(np.abs(vb[ok] - vfl[ok])) < 1e-14  # u = U + u'
    assert np.all(np.isnan(ub[np.hypot(P[0], P[1]) < a]))
    with pytest.raises(ValueError):
        ch03.cylinder_flow(0, 2, U, a, frame="lab")


def test_cylinder_flow_V2_from_streamfunction_div_curl_zero():  # V2
    x, y, U, a = sp.symbols("x y U a", positive=True)
    psi = U * y * (1 - a ** 2 / (x ** 2 + y ** 2))
    u, v = sp.diff(psi, y), -sp.diff(psi, x)
    assert sp.simplify(sp.diff(u, x) + sp.diff(v, y)) == 0 and sp.simplify(sp.diff(v, x) - sp.diff(u, y)) == 0
    un, vn = sp.lambdify((x, y, U, a), u, "numpy"), sp.lambdify((x, y, U, a), v, "numpy")
    P = rand_points(2, 200, -3, 3)
    P = P[:, np.hypot(P[0], P[1]) > 0.75]
    cu, cv = ch03.cylinder_flow(P[0], P[1], 1.1, 0.7)
    assert np.max(np.abs(cu - un(P[0], P[1], 1.1, 0.7))) < 1e-13 and np.max(np.abs(cv - vn(P[0], P[1], 1.1, 0.7))) < 1e-13
    ps = ch03.cylinder_streamfunction(P[0], P[1], 1.1, 0.7)
    assert np.max(np.abs(ps - sp.lambdify((x, y, U, a), psi)(P[0], P[1], 1.1, 0.7))) < 1e-13


def test_frame_terms_V1_worked_number_and_invariance():  # V1  N23, C05 worked number, N02
    out0 = ch03.frame_acceleration_terms(0.0, 1.5, 1.0, 1.0, 0.0)
    assert out0["local"][1] == pytest.approx(-0.592593, abs=5e-7)  # -16/27 exactly
    assert out0["local_y"] == pytest.approx(-16 / 27, abs=1e-8)
    body = ch03.frame_acceleration_terms(0.0, 1.5, 1.0, 1.0, 1.0)
    assert body["steady"] and np.max(np.abs(body["local"])) < 1e-9
    # exact body-frame acceleration at (0, 1.5): u ∂v/∂x with u = 1 + 1/y², ∂v/∂x = -2/y³
    y = 1.5
    assert body["total"][1] == pytest.approx(-(1 + 1 / y ** 2) * 2 / y ** 3, abs=1e-8)
    for Uf in (0.0, 0.25, 0.5, 0.9, 1.0):
        o = ch03.frame_acceleration_terms(0.0, 1.5, 1.0, 1.0, Uf)
        assert np.max(np.abs(o["total"] - body["total"])) < 1e-8
        assert o["u_x"] == pytest.approx(1 + 1 / y ** 2 - 1.0 + Uf, abs=1e-12)
    # t ≠ 0: the same physical point sits at X + (U - U_f) t in the body frame; totals agree there
    t = 0.3
    for Uf in (0.0, 0.4):
        X = np.array([0.2, 1.5])
        o = ch03.frame_acceleration_terms(X[0], X[1], 1.0, 1.0, Uf, t=t)
        b = ch03.frame_acceleration_terms(X[0] + (1.0 - Uf) * t, X[1], 1.0, 1.0, 1.0, t=0.0)
        assert np.max(np.abs(o["total"] - b["total"])) < 1e-7


# =====================================================================================================================
# C06 — relative velocity (3.10); R01–R03 (3.11)–(3.13); D07
# =====================================================================================================================
def test_velocity_gradient_V1_matches_exact_jacobian():  # V1
    for fld in (U2, U3):
        P = rand_points(fld.d, 50)
        G = ch03.velocity_gradient_at(fld, P, 0.4)
        assert G.shape == (fld.d, fld.d, 50)
        assert np.max(np.abs(G - fld.G(P, 0.4))) < 1e-7
        g1 = ch03.velocity_gradient_at(fld, P[:, 0], 0.4)
        assert g1.shape == (fld.d, fld.d) and np.max(np.abs(g1 - fld.G(P[:, 0], 0.4))) < 1e-7
    # index order pinned: u = (x2, 0) has G[0, 1] = 1
    G = ch03.velocity_gradient_at(lambda x, t: np.stack([np.asarray(x)[1], 0 * np.asarray(x)[0]]), np.array([0.3, 0.2]))
    assert np.allclose(G, [[0, 1], [0, 0]], atol=1e-10)


def test_velocity_gradient_V3_order_two_and_taylor_remainder():  # V3
    x0, t0 = np.array([0.3, -0.2, 0.5]), 0.25
    ex = U3.G(x0, t0)
    hs = [0.1, 0.05, 0.025, 0.0125]
    errs = [np.max(np.abs(ch03.velocity_gradient_at(U3, x0, t0, h=h) - ex)) for h in hs]
    assert abs(observed_order(hs, errs) - 2.0) < ORDER_TOL
    # D07: u(x + dx) - u(x) - G dx = O(|dx|²)
    n = np.array([0.6, -0.48, 0.64])
    ss = [0.1, 0.05, 0.025, 0.0125]
    rem = [np.max(np.abs(U3(x0 + s * n, t0) - U3(x0, t0) - ch03.relative_velocity(ex, s * n))) for s in ss]
    assert abs(observed_order(ss, rem) - 2.0) < ORDER_TOL


def test_relative_velocity_V1_worked_number_and_linear_exactness():  # V1
    G = np.array([[1.0, 2.0], [0.0, -1.0]])
    assert np.allclose(ch03.relative_velocity(G, [0.01, 0.02]), [0.05, -0.02], atol=1e-15)
    # linear field: (3.10) is exact for any dx, vectorised over (d, N)
    dX = rand_points(2, 20)
    assert np.max(np.abs(ch03.relative_velocity(G, dX) - G @ dX)) < 1e-15


def test_strain_rotation_split_V1_recap_R_is_G_minus_GT():  # V1  R01–R03 (book R = G − Gᵀ, no ½)
    for _ in range(10):
        G = RNG.normal(size=(3, 3))
        S, R = ch03.strain_rate_tensor(G), ch03.rotation_tensor(G)
        assert np.allclose(S + 0.5 * R, G, atol=1e-15)  # (3.11)
        assert np.allclose(R, G - G.T, atol=1e-15) and not np.allclose(R, 0.5 * (G - G.T))  # (3.13) factor pinned
        assert np.allclose(S, S.T) and np.allclose(R, -R.T)
        assert np.allclose(ch03.antisymmetric_from_vector(ch03.vorticity_from_gradient(G)), R, atol=1e-14)  # (3.15)


# =====================================================================================================================
# C07 — linear strain rate; D08
# =====================================================================================================================
def test_linear_strain_rate_V1_axes_and_quadratic_form():  # V1
    for _ in range(10):
        G = RNG.normal(size=(3, 3))
        assert abs(ch03.linear_strain_rate(G, [1, 0, 0]) - G[0, 0]) < 1e-15  # S₁₁ = ∂u₁/∂x₁
        n = RNG.normal(size=3)
        nh = n / np.linalg.norm(n)
        assert abs(ch03.linear_strain_rate(G, n) - nh @ G @ nh) < 1e-14  # the antisymmetric part drops out
    Ns = rand_points(3, 40)
    G = RNG.normal(size=(3, 3))
    vec = ch03.linear_strain_rate(G, Ns)
    assert vec.shape == (40,)
    # rigid rotation stretches nothing
    Gr = ch03.velocity_gradient_preset("solid_body_rotation", 1.7)
    assert np.max(np.abs(ch03.linear_strain_rate(Gr, rand_points(2, 30)))) < 1e-15
    # worked number: u = (2x, -2y); a 1 cm thread along x is 1.0202 cm after 0.01 s; at 45° no stretching
    Gs = np.diag([2.0, -2.0])
    assert abs(np.linalg.norm(ch03.linear_flow_map(Gs, 0.01) @ [0.01, 0.0]) - 0.01 * math.exp(0.02)) < 1e-15
    assert ch03.linear_strain_rate(Gs, [1, 0]) == 2.0 and abs(ch03.linear_strain_rate(Gs, [1, 1])) < 1e-15


def test_linear_strain_rate_V3_tracked_segment_converges():  # V3 (geometry: ruler on a tracked segment)
    G = RNG.normal(size=(2, 2))
    n = np.array([0.8, 0.6])
    dts = [1e-2, 5e-3, 2.5e-3, 1.25e-3]
    errs = [abs(ch03.measured_strain_rates(G, n, dt)["stretch"] - ch03.linear_strain_rate(G, n)) for dt in dts]
    assert abs(observed_order(dts, errs) - 1.0) < ORDER_TOL
    m = ch03.measured_strain_rates(G, n, 1e-3)
    assert m["stretch_formula"] == pytest.approx(ch03.linear_strain_rate(G, n), abs=1e-15)


def test_linear_strain_rate_V7_rotation_of_axes_invariant():  # V7
    for _ in range(20):
        G = RNG.normal(size=(3, 3))
        n = RNG.normal(size=3)
        C = rand_rotation()
        assert abs(ch03.linear_strain_rate(C.T @ G @ C, C.T @ n) - ch03.linear_strain_rate(G, n)) < 1e-13


def test_linear_strain_rate_V2_derivation():  # V2  D08 steps 3-6
    dt = sp.symbols("dt", positive=True)
    g = sp.Matrix(2, 2, sp.symbols("g11 g12 g21 g22", real=True))
    th = sp.symbols("theta", real=True)
    n = sp.Matrix([sp.cos(th), sp.sin(th)])
    moved = (sp.eye(2) + g * dt) * n  # first order in dt
    ell = sp.sqrt(moved.dot(moved))
    rate = sp.limit((ell - 1) / dt, dt, 0)
    S = (g + g.T) / 2
    assert sp.simplify(rate - (n.T * S * n)[0]) == 0  # (1/ℓ) Dℓ/Dt = n·S·n
    A = (g - g.T) / 2
    assert sp.simplify((n.T * A * n)[0]) == 0  # the antisymmetric part drops (step 6)
    assert sp.simplify(rate.subs(th, 0) - g[0, 0]) == 0  # step 4: along x₁ it is ∂u₁/∂x₁


# =====================================================================================================================
# C08 — shear strain rate; N27 rigid motion; D09, D10
# =====================================================================================================================
def test_shear_strain_rate_V1_off_diagonal_and_gamma_convention():  # V1
    for _ in range(10):
        G = RNG.normal(size=(3, 3))
        assert abs(ch03.shear_strain_rate(G, [1, 0, 0], [0, 1, 0]) - 0.5 * (G[0, 1] + G[1, 0])) < 1e-15
    with pytest.raises(ValueError):
        ch03.shear_strain_rate(np.eye(2), [1, 0], [1, 1])
    # parallel shear γ = du₁/dx₂: S₁₂ = γ/2 — NOT γ (the ch02 Γ ≡ S₁₂ convention would give γ)
    for gam in (0.5, 1.0, 3.0):
        G = ch03.velocity_gradient_preset("simple_shear", gam)
        assert ch03.shear_strain_rate(G, [1, 0], [0, 1]) == pytest.approx(gam / 2, abs=1e-15)
        assert abs(ch03.shear_strain_rate(G, [1, 0], [0, 1]) - gam) > 0.2
    # worked number: γ = 1, dt = 0.01: the vertical thread tilts by dα ≈ 0.01 rad, the horizontal one stays: S₁₂ = 0.5
    G = ch03.velocity_gradient_preset("simple_shear", 1.0)
    d_alpha = math.pi / 2 - ch03.material_line_angle(G, 0.01, math.pi / 2)
    d_beta = ch03.material_line_angle(G, 0.01, 0.0)
    assert abs(d_alpha - math.atan(0.01)) < 1e-15 and abs(d_beta) < 1e-15
    assert abs(0.5 * (d_alpha + d_beta) / 0.01 - 0.5) < 1e-4


def test_shear_strain_rate_V3_tracked_corner_closing():  # V3
    G = RNG.normal(size=(2, 2))
    n1 = np.array([0.6, 0.8])
    n2 = np.array([-0.8, 0.6])
    dts = [1e-2, 5e-3, 2.5e-3, 1.25e-3]
    target = 2 * ch03.shear_strain_rate(G, n1, n2)
    errs = [abs(ch03.measured_strain_rates(G, n1, dt, n2=n2)["closing"] - target) for dt in dts]
    assert abs(observed_order(dts, errs) - 1.0) < ORDER_TOL


def test_rigid_motion_V1_no_strain_vorticity_twice_omega():  # V1  N27, D10 result (ω = 2Ω)
    for _ in range(20):
        U, Om = RNG.normal(size=3), RNG.normal(size=3)
        fld = lambda x, t, U=U, Om=Om: ch03.rigid_body_velocity(U, Om, x)  # noqa: E731
        P = RNG.normal(size=3)
        G = ch03.velocity_gradient_at(fld, P)
        assert np.max(np.abs(ch03.strain_rate_tensor(G))) < 1e-9
        assert np.max(np.abs(ch03.vorticity_from_gradient(G) - 2 * Om)) < 1e-9
        assert np.max(np.abs(ch03.element_rotation_rate(G) - Om)) < 1e-9  # the element spins at Ω, not 2Ω
    Xs = rand_points(3, 7)
    out = ch03.rigid_body_velocity([1, 2, 3], [0, 0, 1], Xs)
    assert out.shape == (3, 7) and np.allclose(out[2], 3) and np.allclose(out[0], 1 - Xs[1])


def test_rigid_motion_V7_strain_frame_independent():  # V7  D10 step 5
    G = RNG.normal(size=(3, 3))
    for _ in range(10):
        Om = RNG.normal(size=3)
        Grig2 = ch03.velocity_gradient_at(lambda x, t: ch03.rigid_body_velocity([0, 0, 0], Om, x), np.zeros(3))
        # the velocity gradient of Ω × x is the antisymmetric matrix of Ω: G_im = ε_ijm Ω_j = -ε_imj Ω_j
        assert np.allclose(ch03.antisymmetric_from_vector(Om), Grig2, atol=1e-9)
        assert np.allclose(ch03.strain_rate_tensor(G + Grig2), ch03.strain_rate_tensor(G), atol=1e-9)
        assert not np.allclose(ch03.vorticity_from_gradient(G + Grig2), ch03.vorticity_from_gradient(G))


def test_shear_strain_rate_V2_derivation():  # V2  D09 steps 1-8
    dt = sp.symbols("dt", positive=True)
    g11, g12, g21, g22 = sp.symbols("g11 g12 g21 g22", real=True)
    g = sp.Matrix([[g11, g12], [g21, g22]])
    M = sp.eye(2) + g * dt
    vert = M * sp.Matrix([0, 1])  # thread along x₂ (top C)
    horiz = M * sp.Matrix([1, 0])  # thread along x₁
    alpha = sp.atan(vert[0] / vert[1])  # clockwise tilt of the vertical thread (step 2)
    beta = sp.atan(horiz[1] / horiz[0])  # counterclockwise tilt of the horizontal thread (step 4)
    assert sp.limit(alpha / dt, dt, 0) == g12 and sp.limit(beta / dt, dt, 0) == g21  # steps 3-4
    rate = sp.limit((alpha + beta) / (2 * dt), dt, 0)  # steps 5-6
    assert sp.simplify(rate - (g12 + g21) / 2) == 0  # step 7: S₁₂
    # step 8: any perpendicular pair (n1, n2) = rotated axes; the rotated component is n1·S·n2
    th = sp.symbols("theta", real=True)
    C = sp.Matrix([[sp.cos(th), -sp.sin(th)], [sp.sin(th), sp.cos(th)]])
    gr = C.T * g * C  # components of G in axes along n1, n2
    Mr = sp.eye(2) + gr * dt
    a_r = sp.atan((Mr * sp.Matrix([0, 1]))[0] / (Mr * sp.Matrix([0, 1]))[1])
    b_r = sp.atan((Mr * sp.Matrix([1, 0]))[1] / (Mr * sp.Matrix([1, 0]))[0])
    rate_r = sp.limit((a_r + b_r) / (2 * dt), dt, 0)
    n1, n2 = C[:, 0], C[:, 1]
    S = (g + g.T) / 2
    assert sp.simplify(rate_r - (n1.T * S * n2)[0]) == 0


def test_rigid_motion_V2_derivation():  # V2  D10 steps 1-5
    U = sp.Matrix(sp.symbols("U1:4", real=True))
    Om = sp.Matrix(sp.symbols("Omega1:4", real=True))
    x = sp.Matrix(sp.symbols("x1:4", real=True))
    u = U + Om.cross(x)  # step 1
    G = u.jacobian(x)  # step 2
    assert sp.simplify((G + G.T) / 2) == sp.zeros(3, 3)  # step 3: S = 0
    R = G - G.T  # step 4
    w = sp.Matrix([R[2, 1], R[0, 2], R[1, 0]])
    assert sp.simplify(w - 2 * Om) == sp.zeros(3, 1)  # ω = 2Ω
    assert sp.simplify(R - 2 * G) == sp.zeros(3, 3)


# =====================================================================================================================
# C09 — volumetric strain rate (3.14); D11
# =====================================================================================================================
def test_volumetric_strain_rate_V1_trace_and_jacobi():  # V1
    for _ in range(20):
        G = RNG.normal(size=(3, 3))
        assert ch03.volumetric_strain_rate(G) == pytest.approx(np.trace(G), abs=1e-15)
        ts = np.array([0.0, 0.3, 1.1])
        vr = ch03.material_volume_ratio(G, ts)
        assert np.max(np.abs(vr / np.exp(ts * np.trace(G)) - 1)) < 1e-12  # det e^{Gt} = e^{t tr G}
    # worked number: u = (x, y, z): 1 cm³ -> e^{0.03} cm³ after 0.01 s (first order 1.03)
    assert ch03.material_volume_ratio(np.eye(3), 0.01) == pytest.approx(math.exp(0.03), rel=1e-14)
    assert ch03.volumetric_strain_rate(np.eye(3)) == 3.0
    # simple shear keeps volume at every time
    Gs = ch03.velocity_gradient_preset("simple_shear", 2.0, dim=3)
    assert abs(ch03.material_volume_ratio(Gs, 5.0) - 1.0) < 1e-12


def test_volumetric_strain_rate_V3_tracked_area_converges():  # V3
    G = RNG.normal(size=(2, 2))
    dts = [1e-2, 5e-3, 2.5e-3, 1.25e-3]
    errs = [abs(ch03.measured_strain_rates(G, [1, 0], dt)["area"] - np.trace(G)) for dt in dts]
    assert abs(observed_order(dts, errs) - 1.0) < ORDER_TOL
    G3 = RNG.normal(size=(3, 3))
    m = ch03.measured_strain_rates(G3, [1, 0, 0], 1e-4, n2=[0, 1, 0])
    assert abs(m["area"] - np.trace(G3)) < 1e-3 and math.isnan(m["spin"])
    with pytest.raises(ValueError):
        ch03.measured_strain_rates(G3, [1, 0, 0], 1e-3)


def test_volumetric_strain_rate_V7_rotation_invariant():  # V7
    G = RNG.normal(size=(3, 3))
    for _ in range(50):
        C = rand_rotation()
        assert abs(ch03.volumetric_strain_rate(C.T @ G @ C) - ch03.volumetric_strain_rate(G)) < 1e-13


def test_volumetric_strain_rate_V2_derivation():  # V2  D11 steps 1-7
    t, dt = sp.symbols("t dt", positive=True)
    l = [sp.Function(f"l{i}")(t) for i in (1, 2, 3)]
    r = sp.symbols("r1:4", real=True)
    dV = l[0] * l[1] * l[2]
    ddV = sp.diff(dV, t)  # step 1
    sub = {sp.Derivative(l[i], t): r[i] * l[i] for i in range(3)}  # step 2: each edge stretches at its own rate
    assert sp.simplify(ddV.subs(sub) / dV - sum(r)) == 0  # steps 3-4
    # step 5: shear changes volume only at second order: det(I + G dt) = 1 + tr G dt + O(dt²)
    g = sp.Matrix(3, 3, sp.symbols("g1:10", real=True))
    detM = sp.expand((sp.eye(3) + g * dt).det())
    assert sp.simplify(detM.coeff(dt, 1) - g.trace()) == 0 and detM.coeff(dt, 0) == 1
    gs = sp.Matrix([[0, sp.Symbol("gamma")], [0, 0]])  # pure shear: exactly volume preserving
    assert sp.simplify((sp.eye(2) + gs * dt).det() - 1) == 0
    # step 6: the trace is invariant under a rotation of axes
    th = sp.symbols("theta", real=True)
    C = sp.Matrix([[sp.cos(th), -sp.sin(th), 0], [sp.sin(th), sp.cos(th), 0], [0, 0, 1]])
    assert sp.simplify((C.T * g * C).trace() - g.trace()) == 0
    # step 7: Jacobi's formula on a symbolic diagonalisable 2×2 (upper triangular): det e^{Gt} = e^{t tr G}
    a, b, c = sp.symbols("a b c", real=True)
    Gt = sp.Matrix([[a, b], [0, c]])
    E = (Gt * t).exp()
    assert sp.simplify(E.det() - sp.exp(t * (a + c))) == 0


# =====================================================================================================================
# C10 — vorticity is twice the element's spin; R04–R07 (3.15)–(3.18); N28, N29, N34; D12, D13, D14
# =====================================================================================================================
def test_spin_V1_perpendicular_pairs_average_half_vorticity():  # V1
    thetas = np.linspace(0, np.pi, 36, endpoint=False)
    Gs = [ch03.velocity_gradient_preset(n, 1.3) for n in ch03.VELOCITY_GRADIENT_PRESETS] + \
         [RNG.normal(size=(2, 2)) for _ in range(5)]
    for G in Gs:
        w3 = ch03.vorticity_from_gradient(G)[2]
        pair = ch03.perpendicular_pair_rotation_rate(G, thetas)
        assert np.max(np.abs(pair - 0.5 * w3)) < 1e-14
        assert np.allclose(ch03.element_rotation_rate(G), [0, 0, 0.5 * w3], atol=1e-15)
        single = ch03.material_line_rotation_rate(G, thetas)  # single threads differ (unless solid body / pure rotation)
        assert single.shape == thetas.shape
    # parallel shear: θ̇ = -γ sin²θ (design contract: γ = 1, θ = 30° → -0.25)
    G = ch03.velocity_gradient_preset("simple_shear", 1.0)
    assert ch03.material_line_rotation_rate(G, math.radians(30)) == pytest.approx(-0.25, abs=1e-15)
    assert np.allclose(ch03.material_line_rotation_rate(G, thetas), -np.sin(thetas) ** 2, atol=1e-15)
    # solid body: every line turns at ω₀, the element spins at ω₀, the vorticity is 2ω₀ (the factor-2 convention)
    Gr = ch03.velocity_gradient_preset("solid_body_rotation", 0.8)
    assert np.allclose(ch03.material_line_rotation_rate(Gr, thetas), 0.8, atol=1e-15)
    assert ch03.vorticity_from_gradient(Gr)[2] == pytest.approx(1.6, abs=1e-15)
    assert ch03.element_rotation_rate(Gr)[2] == pytest.approx(0.8, abs=1e-15)
    for t in (0.3, 1.0, 2.5):
        assert abs(ch03.material_line_angle(Gr, t, 0.2) - (0.2 + 0.8 * t)) < 1e-13


def test_spin_V3_tracked_pair_converges():  # V3
    G = RNG.normal(size=(2, 2))
    dts = [1e-2, 5e-3, 2.5e-3, 1.25e-3]
    target = 0.5 * ch03.vorticity_from_gradient(G)[2]
    errs = [abs(ch03.measured_strain_rates(G, [0.28, 0.96], dt)["spin"] - target) for dt in dts]
    assert abs(observed_order(dts, errs) - 1.0) < ORDER_TOL
    assert ch03.measured_strain_rates(G, [1, 0], 1e-3)["spin_formula"] == pytest.approx(target, abs=1e-15)


def test_vorticity_V1_equals_curl_of_field():  # V1  R05, R06 (3.15), (3.16)
    P = rand_points(3, 50)
    w = ch03.vorticity(U3, P, 0.3)
    assert w.shape == (3, 50) and np.max(np.abs(w - U3.curl(P, 0.3))) < 1e-7
    w2 = ch03.vorticity(U2, P[:2], 0.3)
    assert np.max(np.abs(w2[:2])) == 0 and np.max(np.abs(w2[2] - U2.curl(P[:2], 0.3)[2])) < 1e-7
    # a single 3 × 3 G and the ε form R_ij = -ε_ijk ω_k
    G = RNG.normal(size=(3, 3))
    om = ch03.vorticity_from_gradient(G)
    assert np.allclose(om, [G[2, 1] - G[1, 2], G[0, 2] - G[2, 0], G[1, 0] - G[0, 1]], atol=1e-15)
    R = -np.einsum("ijk,k->ij", ch03.levi_civita(), om)
    assert np.allclose(R, ch03.rotation_tensor(G), atol=1e-15)
    # the ch02 recap: vector_from_antisymmetric(rotation_tensor) = vorticity; Ex. 2.3 b × x has curl 2b
    b = np.array([0.3, -1.2, 0.7])
    Gb = ch03.velocity_gradient_at(ch03.from_coord_field(ch02.solid_body_rotation_field(b)), np.array([0.2, 0.1, -0.3]))
    assert np.allclose(ch03.vorticity_from_gradient(Gb), 2 * b, atol=1e-9)
    assert np.allclose(ch03.vector_from_antisymmetric(ch03.rotation_tensor(Gb)), 2 * b, atol=1e-9)


def test_vorticity_rotating_frame_V1_curl_of_relative_velocity():  # V1  N28, D13 result (ω′ = ω − 2Ω)
    for _ in range(5):
        Om = RNG.normal(size=3)
        rel = lambda x, t, Om=Om: U3(x, t) - ch03.rigid_body_velocity([0, 0, 0], Om, x)  # noqa: E731
        P = rand_points(3, 20)
        w_rel = ch03.vorticity(rel, P, 0.2)
        w_abs = ch03.vorticity(U3, P, 0.2)
        pred = ch03.vorticity_in_rotating_frame(w_abs, Om[:, None])
        assert np.max(np.abs(w_rel - pred)) < 1e-7
        assert np.max(np.abs(w_rel - (w_abs - Om[:, None]))) > 0.1  # the "ω − Ω" slip is rejected
    # solid body ω₀ = 1 (ω_z = 2) seen from a turntable at Ω = 1: at rest
    assert ch03.vorticity_in_rotating_frame(2.0, 1.0) == 0.0
    assert ch03.vorticity_in_rotating_frame(2.0, 0.0) == 2.0
    # Earth: planetary vorticity 2Ω ≈ 1.46e-4 1/s dwarfs a synoptic ζ ~ 1e-5
    assert ch03.vorticity_in_rotating_frame(1e-5 + 2 * 7.292e-5, 7.292e-5) == pytest.approx(1e-5, rel=1e-9)


def test_parallel_shear_V1_kinematics():  # V1  N34, N35 (Fig. 3.14)
    for gam in (0.5, 1.0, 2.5):
        k = ch03.parallel_shear_kinematics(gam)
        assert k["omega3"] == -gam and k["spin"] == -gam / 2  # clockwise
        assert k["rate_AB"] == pytest.approx(-gam, abs=1e-15) and k["rate_BC"] == 0.0
        assert np.allclose(k["S"], [[0, gam / 2], [gam / 2, 0]]) and np.allclose(k["R"], [[0, gam], [-gam, 0]])
        assert np.allclose(k["lam"], [-gam / 2, gam / 2]) and np.allclose(k["S_bar"], np.diag([gam / 2, -gam / 2]))
        assert k["principal_angle_deg"] == pytest.approx(45.0, abs=1e-12)
        assert np.allclose(k["axes"][:, 1], [1 / math.sqrt(2), 1 / math.sqrt(2)])
        G = k["G"]
        n45, n135 = np.array([1, 1]) / math.sqrt(2), np.array([-1, 1]) / math.sqrt(2)
        # aligned element ABCD: shears (S₁₂ = γ/2) without stretching; 45° element PQRS: stretches ±γ/2 without shear
        assert ch03.shear_strain_rate(G, [1, 0], [0, 1]) == pytest.approx(gam / 2)
        assert ch03.linear_strain_rate(G, [1, 0]) == 0.0 and ch03.linear_strain_rate(G, [0, 1]) == 0.0
        assert ch03.linear_strain_rate(G, n45) == pytest.approx(gam / 2) and ch03.linear_strain_rate(G, n135) == pytest.approx(-gam / 2)
        assert abs(ch03.shear_strain_rate(G, n45, n135)) < 1e-15
        # both elements spin at -γ/2
        assert ch03.perpendicular_pair_rotation_rate(G, 0.0) == pytest.approx(-gam / 2)
        assert ch03.perpendicular_pair_rotation_rate(G, math.pi / 4) == pytest.approx(-gam / 2)


def test_circulation_V1_solid_body_and_line_vortex():  # V1  R07 (3.18), N37 (3.24), N39 (3.26)
    om0 = 1.0
    for r in (0.5, 1.0, 2.0):
        # solid body ω₀ = Γ/(2πσ²) with Γ = 2π, σ = 1 → ω₀ = 1: Γ(r) = 2π r² ω₀, centred and off-centre ("any circuit")
        assert ch03.circulation_circle("solid", r, Gamma=2 * np.pi, sigma=1.0) == pytest.approx(2 * np.pi * r ** 2 * om0, rel=1e-14)
        assert ch03.circulation_circle("solid", r, center=(0.3, -0.2), Gamma=2 * np.pi, sigma=1.0) == \
            pytest.approx(2 * np.pi * r ** 2 * om0, rel=1e-12)
        assert ch03.circulation_circle(lambda rr: 0.7 / rr, r) == pytest.approx(2 * np.pi * 0.7, rel=1e-14)  # 2πB
    assert ch03.circulation_circle("line", 1.0, center=(0.2, 0.1), Gamma=2 * np.pi) == pytest.approx(2 * np.pi, rel=1e-10)
    assert abs(ch03.circulation_circle("line", 0.5, center=(2.0, 0.0), Gamma=2 * np.pi)) < 1e-12  # excludes the axis
    # Stokes cross-check with the ch02 integral theorems for the Rankine core: loop integral = ω × area
    s = 2.0
    Gam = 2 * np.pi
    loop = ch03.planar_loop((0.3, 0.1), radius=0.5, n=256)
    circ = ch03.circulation(ch03.as_coord_field(ch03.vortex_velocity_field("rankine", Gamma=Gam, sigma=s)), loop)
    assert circ == pytest.approx(Gam / (np.pi * s ** 2) * np.pi * 0.25, rel=1e-12)


def test_irrotational_potential_V2_curl_of_gradient_vanishes():  # V2  N29 (3.17)
    phi2 = sp.Function("phi")(X1, X2)
    u, w = ch03.potential_velocity(phi2, [X1, X2])
    assert w == 0 and u == [sp.diff(phi2, X1), sp.diff(phi2, X2)]
    phi3 = sp.Function("phi")(X1, X2, X3)
    u, w = ch03.potential_velocity(phi3, [X1, X2, X3])
    assert w == [0, 0, 0]


def test_velocity_potential_V1_reconstruction_and_line_vortex_path_dependence():  # V1  N29, D18 (simply connected)
    fld = lambda x, t: np.stack([np.cos(np.asarray(x)[0]) * np.exp(np.asarray(x)[1]),  # noqa: E731  u = ∇(sin x e^y)
                                 np.sin(np.asarray(x)[0]) * np.exp(np.asarray(x)[1])])
    xs = np.linspace(-1, 1, 201)
    phi, pd = ch03.velocity_potential_2d(fld, (xs, xs), x_ref=(0.0, 0.0))
    X, Y = np.meshgrid(xs, xs, indexing="xy")
    exact = np.sin(X) * np.exp(Y)
    assert np.max(np.abs(phi - exact)) < 1e-4 and pd < 1e-4  # trapezoid O(h²), h = 0.01
    # line vortex: irrotational off-axis but NOT a potential flow in a region around the axis: the two routes differ by 2πB
    B = 1.0
    xs2 = np.linspace(-2, 2, 201) + 0.0123
    ys2 = np.linspace(-2, 2, 201) + 0.0071
    _, pd2 = ch03.velocity_potential_2d(ch03.vortex_velocity_field("line", Gamma=2 * np.pi * B), (xs2, ys2), x_ref=(-2, -2))
    assert abs(pd2 - 2 * np.pi * B) < 1e-3 * 2 * np.pi * B
    # the ch02 grid operators agree that it is irrotational off the axis
    g = ch03.velocity_gradient_at(ch03.vortex_velocity_field("line", Gamma=2 * np.pi), rand_points(2, 20, 0.5, 2.0))
    assert np.max(np.abs(g[1, 0] - g[0, 1])) < 1e-7


def test_spin_V2_derivation():  # V2  D12 steps 1-8
    dt = sp.symbols("dt", positive=True)
    g = sp.Matrix(2, 2, sp.symbols("g11 g12 g21 g22", real=True))
    M = sp.eye(2) + g * dt
    alpha = sp.atan((M * sp.Matrix([0, 1]))[0] / (M * sp.Matrix([0, 1]))[1])  # clockwise tilt of the vertical thread
    beta = sp.atan((M * sp.Matrix([1, 0]))[1] / (M * sp.Matrix([1, 0]))[0])  # counterclockwise tilt of the horizontal one
    spin = sp.limit((-alpha + beta) / (2 * dt), dt, 0)  # steps 1-4
    R = g - g.T
    assert sp.simplify(spin - (-R[0, 1] / 2)) == 0 and sp.simplify(spin - R[1, 0] / 2) == 0  # step 5
    eps = sp.LeviCivita
    w3 = g[1, 0] - g[0, 1]  # (3.16)
    assert sp.simplify(R[1, 0] - (-eps(2, 1, 3) * w3)) == 0  # step 6: R₂₁ = -ε₂₁₃ ω₃ = ω₃
    assert sp.simplify(spin - w3 / 2) == 0
    # averaging α + β instead (the strain rate) gives S₁₂ — a different quantity unless G is antisymmetric
    assert sp.simplify(sp.limit((alpha + beta) / (2 * dt), dt, 0) - spin) != 0


def test_rotating_frame_V2_derivation():  # V2  D13 steps 1-6
    x = sp.Matrix(sp.symbols("x1:4", real=True))
    u = sp.Matrix([sp.Function(f"u{i}")(*x) for i in (1, 2, 3)])
    Om = sp.symbols("Omega", real=True)
    Omv = sp.Matrix([0, 0, Om])

    def curl(v):
        return sp.Matrix([sp.diff(v[2], x[1]) - sp.diff(v[1], x[2]), sp.diff(v[0], x[2]) - sp.diff(v[2], x[0]),
                          sp.diff(v[1], x[0]) - sp.diff(v[0], x[1])])
    up = u - Omv.cross(x)  # step 1
    assert sp.simplify(curl(Omv.cross(x)) - 2 * Omv) == sp.zeros(3, 1)  # step 3
    assert sp.simplify(curl(up) - (curl(u) - 2 * Omv)) == sp.zeros(3, 1)  # steps 2, 4, 5
    w3 = curl(u)[2]
    assert sp.solve(sp.Eq(w3 - 2 * Om, 0), Om) == [w3 / 2]  # step 6: co-rotating frame


def test_shear_pair_spin_V2_derivation():  # V2  D14 steps 1-6
    gam, th = sp.symbols("gamma theta", real=True)
    G = sp.Matrix([[0, gam], [0, 0]])
    e = sp.Matrix([sp.cos(th), sp.sin(th)])
    et = sp.Matrix([-sp.sin(th), sp.cos(th)])
    Ge = G * e
    assert sp.simplify(Ge - sp.Matrix([gam * sp.sin(th), 0])) == sp.zeros(2, 1)  # step 2
    rate = sp.simplify(et.dot(Ge))
    assert sp.simplify(rate + gam * sp.sin(th) ** 2) == 0  # step 3
    assert rate.subs(th, sp.pi / 2) == -gam and rate.subs(th, 0) == 0  # step 4
    partner = rate.subs(th, th + sp.pi / 2)
    assert sp.simplify(partner + gam * sp.cos(th) ** 2) == 0  # step 5
    assert sp.simplify((rate + partner) / 2 + gam / 2) == 0  # step 6: -γ/2 = ω₃/2 for every θ
    assert rate.subs({gam: 1, th: sp.pi / 6}) == sp.Rational(-1, 4)


# =====================================================================================================================
# C11 — relative velocity = deformation + rigid rotation (3.19); N30, N33; D15
# =====================================================================================================================
def test_relative_velocity_split_V1_parts_sum_and_rotation_is_rigid():  # V1
    for d in (2, 3):
        for _ in range(10):
            G = RNG.normal(size=(d, d))
            dX = rand_points(d, 25)
            du, du_s, du_r = ch03.relative_velocity_split(G, dX)
            assert np.max(np.abs(du - G @ dX)) < 1e-14
            assert np.max(np.abs(du_s + du_r - du)) < 1e-14
            assert np.max(np.abs(np.sum(du_r * dX, axis=0))) < 1e-14  # rotation part ⟂ dx (changes no distance)
            assert np.max(np.abs(du_s - ch03.strain_rate_tensor(G) @ dX)) < 1e-14
            w = ch03.vorticity_from_gradient(G)
            if d == 3:
                assert np.max(np.abs(du_r - 0.5 * np.cross(w[:, None], dX, axis=0))) < 1e-14
    # worked number: shear γ = 1, dx = (0, 1): (0.5, 0) + (0.5, 0) = (1, 0)
    sp_ = ch03.relative_velocity_split(ch03.velocity_gradient_preset("simple_shear", 1.0), [0.0, 1.0])
    assert np.allclose(sp_.du_strain, [0.5, 0]) and np.allclose(sp_.du_rot, [0.5, 0]) and np.allclose(sp_.du, [1, 0])


def test_relative_velocity_split_V2_derivation():  # V2  D15 steps 1-6 (and the sign trap of step 4)
    g = sp.Matrix(3, 3, sp.symbols("g1:10", real=True))
    dx = sp.Matrix(sp.symbols("dx1:4", real=True))
    S, R = (g + g.T) / 2, g - g.T
    w = sp.Matrix([R[2, 1], R[0, 2], R[1, 0]])
    eps = sp.LeviCivita
    Rw = sp.Matrix(3, 3, lambda i, j: -sum(eps(i, j, k) * w[k] for k in range(3)))
    assert sp.simplify(Rw - R) == sp.zeros(3, 3)  # (3.15)
    step2 = sp.Matrix([sum((S[i, j] - sp.Rational(1, 2) * sum(eps(i, j, k) * w[k] for k in range(3))) * dx[j]
                           for j in range(3)) for i in range(3)])
    assert sp.simplify(step2 - g * dx) == sp.zeros(3, 1)  # steps 1-3
    rot = sp.Matrix([sp.Rational(1, 2) * sum(eps(i, k, j) * w[k] * dx[j] for j in range(3) for k in range(3))
                     for i in range(3)])  # step 4: ε_ikj swap absorbs the minus
    assert sp.simplify(rot - w.cross(dx) / 2) == sp.zeros(3, 1)  # step 5
    assert sp.simplify(S * dx + w.cross(dx) / 2 - g * dx) == sp.zeros(3, 1)  # step 6: (3.19)
    assert sp.simplify(S * dx - w.cross(dx) / 2 - g * dx) != sp.zeros(3, 1)  # the sign slip is caught
    # the index machine agrees: -1/2 eps_ijk w_k dx_j expanded equals ½(ω × dx)
    e = ch03.expand_indices("-1/2 eps_ijk w_k dx_j")
    names = {str(s): s for c in e for s in c.free_symbols}
    wsym = sp.Matrix([names[f"w_{i}"] for i in (1, 2, 3)])
    dsym = sp.Matrix([names[f"dx_{i}"] for i in (1, 2, 3)])
    assert all(sp.expand(e[i] - (wsym.cross(dsym) / 2)[i]) == 0 for i in range(3))


# =====================================================================================================================
# C12 — principal strain axes; a small sphere becomes an ellipsoid (3.20)–(3.21); N31, N32, N35; D16
# =====================================================================================================================
def test_principal_strain_V1_eigenframe():  # V1
    for d in (2, 3):
        for _ in range(10):
            G = RNG.normal(size=(d, d))
            lam, C = ch03.principal_strain_rates(G)
            S = ch03.strain_rate_tensor(G)
            assert np.all(np.diff(lam) >= 0)
            assert np.allclose(S @ C, C * lam, atol=1e-13) and np.allclose(C.T @ C, np.eye(d), atol=1e-13)
            assert np.linalg.det(C) == pytest.approx(1.0, abs=1e-12)
            assert np.allclose(C.T @ S @ C, np.diag(lam), atol=1e-13)  # (3.20): S̄ diagonal
            v = C[:, -1]
            assert v[np.argmax(np.abs(v) > 1e-12)] > 0  # documented sign convention
            dX = rand_points(d, 10)
            du_bar, dx_bar = ch03.strain_velocity_principal(G, dX)
            assert np.max(np.abs(du_bar - lam[:, None] * dx_bar)) < 1e-13  # (3.21)
    lam, C = ch03.principal_strain_rates(ch03.velocity_gradient_preset("simple_shear", 1.0))
    assert np.allclose(lam, [-0.5, 0.5]) and np.allclose(C[:, 1], [2 ** -0.5, 2 ** -0.5])


def test_strain_ellipse_V1_exact_axes_and_first_order():  # V1
    # shear γ = 1 at t = 0.1 s: exact semi-axes 1.0513, 0.9513 at 43.6° (D16 check numbers)
    G = ch03.velocity_gradient_preset("simple_shear", 1.0)
    ax, dirs = ch03.strain_ellipse_axes(G, 0.1, "exact")
    assert np.allclose(ax, [math.sqrt(1 + 0.0025) + 0.05, math.sqrt(1 + 0.0025) - 0.05], atol=1e-14)
    assert math.degrees(math.atan2(abs(dirs[1, 0]), abs(dirs[0, 0]))) == pytest.approx(0.5 * math.degrees(math.atan(20)), abs=1e-9)
    ax1, _ = ch03.strain_ellipse_axes(G, 0.1)  # first order: 1 ± 0.05 at ±45°
    assert np.allclose(ax1, [1.05, 0.95])
    ax2, d2 = ch03.strain_ellipse_axes(G, 0.1, "strain_only")
    assert np.allclose(ax2, [math.exp(0.05), math.exp(-0.05)]) and abs(abs(d2[0, 0]) - 2 ** -0.5) < 1e-14
    # product of the exact semi-axes = det e^{Gt}; deformed circle lies on the exact ellipse
    for _ in range(5):
        G = RNG.normal(size=(2, 2))
        t = 0.37
        ax, dirs = ch03.strain_ellipse_axes(G, t, "exact", radius=0.5)
        assert np.prod(ax) == pytest.approx(0.25 * ch03.material_volume_ratio(G, t), rel=1e-12)
        P = ch03.deform_circle(G, t, n=90, radius=0.5)
        M = expm(G * t)
        assert np.max(np.abs(np.linalg.norm(np.linalg.solve(M, P), axis=0) - 0.5)) < 1e-13
        loc = dirs.T @ P  # coordinates along the ellipse axes
        assert np.max(np.abs((loc[0] / ax[0]) ** 2 + (loc[1] / ax[1]) ** 2 - 1)) < 1e-12
    Ps = ch03.deform_sphere(np.diag([0.5, 0.0, -0.5]), 1.0, radius=2.0)
    assert Ps.shape == (3, 24 * 48)
    e = np.exp([0.5, 0.0, -0.5]) * 2.0
    assert np.max(np.abs(np.sum((Ps / e[:, None]) ** 2, axis=0) - 1)) < 1e-12
    with pytest.raises(ValueError):
        ch03.strain_ellipse_axes(G, 0.1, "nope")


def test_strain_ellipse_V3_first_order_error_is_second_order():  # V3  D16: first-order statement, O(dt²) error
    G = np.array([[0.7, 0.3], [0.3, -0.4]])  # pure strain (symmetric): axes do not drift
    ts = [0.08, 0.04, 0.02, 0.01]
    errs = [np.max(np.abs(ch03.strain_ellipse_axes(G, t)[0] - ch03.strain_ellipse_axes(G, t, "exact")[0])) for t in ts]
    assert abs(observed_order(ts, errs) - 2.0) < ORDER_TOL
    assert np.allclose(ch03.strain_ellipse_axes(G, 0.3, "strain_only")[0], ch03.strain_ellipse_axes(G, 0.3, "exact")[0],
                       atol=1e-13)  # G = S: strain alone is exact


def test_deform_square_V1_fig3_14_elements():  # V1  N35 (our Fig. 3.14)
    G = ch03.velocity_gradient_preset("simple_shear", 1.0)
    sq = ch03.deform_square(G, 0.3, n_side=5, boundary_only=True)
    assert sq.shape == (2, 20)
    # every point moved by x ← x + γ t y
    ref = ch03.deform_square(np.zeros((2, 2)), 0.0, n_side=5, boundary_only=True)
    assert np.allclose(sq[0], ref[0] + 0.3 * ref[1]) and np.allclose(sq[1], ref[1])
    assert ch03.deform_square(G, 0.3, n_side=4).shape == (2, 16)
    assert ch03.deform_square(ch03.velocity_gradient_preset("simple_shear", 1.0, dim=3), 0.3).shape == (3, 100)


def test_principal_axes_V2_derivation():  # V2  D16 steps 2-7
    dt, eps_ = sp.symbols("dt epsilon", positive=True)
    s1, s2, s3 = sp.symbols("S1 S2 S3", real=True)
    Sbar = sp.diag(s1, s2, s3)
    dxb = sp.Matrix(sp.symbols("X1:4", real=True))
    dub = Sbar * dxb  # (3.20) → (3.21)
    assert all(sp.simplify(dub[i] - [s1, s2, s3][i] * dxb[i]) == 0 for i in range(3))
    moved = dxb + dub * dt  # step 4
    semi = [eps_ * (1 + s * dt) for s in (s1, s2, s3)]  # step 7
    # a point on the sphere (step 5) lands on the ellipsoid with those semi-axes (step 6)
    th, ph = sp.symbols("theta phi", real=True)
    on_sphere = {dxb[0]: eps_ * sp.sin(th) * sp.cos(ph), dxb[1]: eps_ * sp.sin(th) * sp.sin(ph), dxb[2]: eps_ * sp.cos(th)}
    ell = sum((moved[i] / semi[i]) ** 2 for i in range(3)).subs(on_sphere)
    assert sp.simplify(ell - 1) == 0
    # volume ratio to first order = 1 + S_ii dt (the (3.14) check)
    vol = sp.expand(semi[0] * semi[1] * semi[2] / eps_ ** 3)
    assert sp.simplify(vol.coeff(dt, 1) - (s1 + s2 + s3)) == 0


# =====================================================================================================================
# C13 — vorticity in polar coordinates (3.23); solid body (3.22) vs line vortex (3.25); N36–N41; D17
# =====================================================================================================================
def _polar_test_field():
    ur = lambda r, th: r ** 2 * np.sin(th)  # noqa: E731  non-axisymmetric test pair
    ut = lambda r, th: r * np.cos(th)  # noqa: E731
    w_exact = lambda r, th: 2 * np.cos(th) - r * np.cos(th)  # noqa: E731  (1/r)∂(r u_θ)/∂r - (1/r)∂u_r/∂θ

    def cart(x, t=0.0):
        x_ = np.asarray(x, dtype=float)
        r, th = np.hypot(x_[0], x_[1]), np.arctan2(x_[1], x_[0])
        return ch03.cartesian_components(np.stack([ur(r, th), ut(r, th)]), x_, "polar")
    return ur, ut, w_exact, cart


def test_polar_vorticity_V1_equals_cartesian_curl():  # V1
    ur, ut, w_exact, cart = _polar_test_field()
    P = rand_points(2, 30, 0.3, 1.5)
    r, th = np.hypot(P[0], P[1]), np.arctan2(P[1], P[0])
    wp = ch03.polar_vorticity_z(ur, ut, r, th)
    wc = ch03.vorticity(cart, P)[2]
    assert np.max(np.abs(wp - wc)) < 1e-7 and np.max(np.abs(wp - w_exact(r, th))) < 1e-8
    # the four profiles through the kind-name interface
    rr = np.array([0.3, 0.9, 1.4, 2.5])
    assert np.allclose(ch03.polar_vorticity_z(None, lambda q: 0.8 * q, rr), 1.6, atol=1e-9)  # 2ω₀ everywhere
    assert np.allclose(ch03.polar_vorticity_z(None, lambda q: 0.7 / q, rr), 0.0, atol=1e-9)  # irrotational off-axis
    s, Gm = 1.2, 2 * np.pi
    _, w_rk = ch03.rankine_vortex(rr, Gm, s)
    inside = rr < s - 1e-4
    assert np.allclose(ch03.polar_vorticity_z(None, "rankine", rr[inside], Gamma=Gm, sigma=s), w_rk[inside], atol=1e-8)
    assert np.allclose(ch03.polar_vorticity_z(None, "rankine", rr[rr > s], Gamma=Gm, sigma=s), 0.0, atol=1e-8)
    assert np.allclose(ch03.polar_vorticity_z(None, "gaussian", rr, Gamma=Gm, sigma=s), ch03.gaussian_vortex(rr, Gm, s)[1],
                       atol=1e-8)


def test_polar_vorticity_V2_symbolic_formula():  # V2
    r, th, w0, B = sp.symbols("r theta omega0 B", positive=True)
    assert ch03.polar_vorticity_z_sym(0, w0 * r, r, th) == 2 * w0  # (3.22)
    assert ch03.polar_vorticity_z_sym(0, B / r, r, th) == 0  # (3.25)
    # non-axisymmetric pair vs the Cartesian curl written in e_r, e_θ
    x, y = sp.symbols("x y", real=True)
    ur, ut = r ** 2 * sp.sin(th), r * sp.cos(th)
    rr, tt = sp.sqrt(x ** 2 + y ** 2), sp.atan2(y, x)
    u = (ur * sp.cos(th) - ut * sp.sin(th)).subs({r: rr, th: tt})
    v = (ur * sp.sin(th) + ut * sp.cos(th)).subs({r: rr, th: tt})
    curl = sp.simplify(sp.diff(v, x) - sp.diff(u, y))
    formula = ch03.polar_vorticity_z_sym(ur, ut, r, th).subs({r: rr, th: tt})
    val = sp.lambdify((x, y), curl - formula)(0.7, -0.4)
    assert abs(val) < 1e-12


def test_polar_vorticity_V3_stencil_order_two():  # V3
    ur, ut, w_exact, _ = _polar_test_field()
    hs = [0.1, 0.05, 0.025, 0.0125]
    ug = lambda r, th: np.sin(3 * r) * np.cos(2 * th) + r ** 3  # noqa: E731  cubic-and-trig: nonzero third derivative
    ex = lambda r, th: (np.sin(3 * r) + 3 * r * np.cos(3 * r)) * np.cos(2 * th) / r + 4 * r ** 2 - r * np.cos(th)  # noqa: E731,E501
    urr = lambda r, th: r ** 2 * np.sin(th)  # noqa: E731
    errs = [abs(ch03.polar_vorticity_z(urr, ug, 0.9, 0.4, h=h) - ex(0.9, 0.4)) for h in hs]
    assert abs(observed_order(hs, errs) - 2.0) < ORDER_TOL, pairwise_orders(hs, errs)


def test_sector_circulation_V1_line_vortex_zero_rankine_core_uniform():  # V1  N41, Fig. 3.16
    B = 0.9
    for r0, dr, dth in ((0.5, 0.2, 0.3), (1.0, 0.05, 1.0), (2.0, 1.0, 2.0)):
        assert abs(ch03.annular_sector_circulation(lambda q: B / q, r0, dr, dth)) < 1e-12
        legs = ch03.annular_sector_circulation(lambda q: B / q, r0, dr, dth, legs=True)
        assert legs["outer"] == pytest.approx(B * dth, rel=1e-13) and legs["inner"] == pytest.approx(-B * dth, rel=1e-13)
        assert legs["area"] == pytest.approx(0.5 * ((r0 + dr) ** 2 - r0 ** 2) * dth, rel=1e-14)
    # Rankine core (uniform ω = Γ/πσ²): Γ_sector = ω × area exactly (design check: Γ = 2π, σ = 2 → 0.5 1/s)
    lg = ch03.annular_sector_circulation("rankine", 0.5, 0.3, 0.7, Gamma=2 * np.pi, sigma=2.0, legs=True)
    assert lg["total"] / lg["area"] == pytest.approx(0.5, rel=1e-12)
    # radial legs carry u_r: a pure radial field u_r = f(θ) makes circulation from ∂u_r/∂θ
    lg = ch03.annular_sector_circulation(lambda q: 0.0 * q, 1.0, 0.1, 0.2, u_r=lambda q, th: np.sin(th), legs=True)
    exact = 0.1 * (np.sin(0.0) - np.sin(0.2))
    assert lg["total"] == pytest.approx(exact, rel=1e-12)


def test_sector_circulation_V3_circulation_per_area_converges_to_vorticity():  # V3  D17 in numbers
    ur, ut, w_exact, _ = _polar_test_field()
    r0, th0 = 1.1, 0.6
    ds = [0.2, 0.1, 0.05, 0.025]
    errs = []
    for d in ds:
        lg = ch03.annular_sector_circulation(ut, r0 - d / 2, d, d / r0, u_r=ur, theta0=th0 - d / (2 * r0), legs=True)
        errs.append(abs(lg["total"] / lg["area"] - w_exact(r0, th0)))
    assert abs(observed_order(ds, errs) - 2.0) < ORDER_TOL, pairwise_orders(ds, errs)


def test_mean_vorticity_disc_V1_delta_core_slope_minus_two():  # V1  N40 (3.27)
    B = 0.6
    rs = np.logspace(-3, 0, 12)
    m = ch03.mean_vorticity_in_disc(lambda q: B / q, rs)
    assert np.allclose(m, 2 * B / rs ** 2, rtol=1e-14)
    slope = np.polyfit(np.log(rs), np.log(m), 1)[0]
    assert abs(slope + 2.0) < 1e-12
    assert np.allclose(ch03.mean_vorticity_in_disc("solid", rs, Gamma=2 * np.pi, sigma=1.0), 2.0, rtol=1e-14)
    # Γ = πr² × mean vorticity = 2πB for every r
    assert np.allclose(np.pi * rs ** 2 * m, 2 * np.pi * B, rtol=1e-14)


def test_vortex_fields_V1_cartesian_twins_and_curl():  # V1  N36, N38
    P = rand_points(2, 40, -2, 2)
    P = P[:, np.hypot(P[0], P[1]) > 0.2]
    sb = ch03.vortex_velocity_field(lambda q: ch03.solid_body_rotation(q, 0.7))
    ch2 = ch02.solid_body_rotation_field(0.7, dim=2)
    assert np.max(np.abs(sb(P) - ch2(P[0], P[1]))) < 1e-14
    assert np.max(np.abs(ch03.vorticity(sb, P, h=1e-5)[2] - 1.4)) < 1e-8
    lv = ch03.vortex_velocity_field("line", Gamma=2 * np.pi * 0.5)
    assert np.max(np.abs(lv(P) - ch02.irrotational_vortex_field(0.5)(P[0], P[1]))) < 1e-14
    assert np.max(np.abs(ch03.vorticity(lv, P, h=1e-5)[2])) < 1e-7  # h² u''' ∝ h²/r⁴ near r = 0.2
    gv = ch03.vortex_velocity_field("gaussian", Gamma=2.0, sigma=0.8)
    r = np.hypot(P[0], P[1])
    assert np.max(np.abs(ch03.vorticity(gv, P, h=1e-5)[2] - ch03.gaussian_vortex(r, 2.0, 0.8)[1])) < 1e-7
    assert np.allclose(sb(np.zeros(2)), 0.0) and np.all(np.isnan(lv(np.zeros(2))))


def test_polar_vorticity_V2_derivation():  # V2  D17 steps 1-10 (circulation round a sector, exactly)
    r, th, d = sp.symbols("r theta delta", positive=True)
    a = sp.symbols("a0:5", real=True)
    b = sp.symbols("b0:5", real=True)
    ur = a[0] + a[1] * r + a[2] * th + a[3] * r * th + a[4] * th ** 2
    ut = b[0] + b[1] * r + b[2] * r ** 2 + b[3] * th + b[4] * r * th
    r0, t0 = sp.symbols("r0 theta0", positive=True)
    q = sp.symbols("q", positive=True)
    radial_start = sp.integrate(ur.subs({r: q, th: t0}), (q, r0, r0 + d))  # step 5
    outer = sp.integrate(ut.subs({r: r0 + d, th: q}) * (r0 + d), (q, t0, t0 + d))  # step 2
    radial_end = -sp.integrate(ur.subs({r: q, th: t0 + d}), (q, r0, r0 + d))  # step 6
    inner = -sp.integrate(ut.subs({r: r0, th: q}) * r0, (q, t0, t0 + d))  # step 3
    area = ((r0 + d) ** 2 - r0 ** 2) * d / 2  # step 1 (exact; r dr dθ to leading order)
    Gam = radial_start + outer + radial_end + inner  # step 8
    omega = sp.limit(sp.simplify(Gam / area), d, 0)  # step 9
    formula = (sp.diff(r * ut, r) / r - sp.diff(ur, th) / r).subs({r: r0, th: t0})
    assert sp.simplify(omega - formula) == 0
    # step 4: the arcs combine to ∂(r u_θ)/∂r dr dθ at leading order (the r stays inside the derivative)
    arcs = sp.series(sp.simplify((outer + inner) / d ** 2), d, 0, 1).removeO()
    assert sp.simplify(arcs - sp.diff(r * ut, r).subs({r: r0, th: t0})) == 0
    # step 7: the radial legs give -∂u_r/∂θ dr dθ at leading order
    rad = sp.series(sp.simplify((radial_start + radial_end) / d ** 2), d, 0, 1).removeO()
    assert sp.simplify(rad + sp.diff(ur, th).subs({r: r0, th: t0})) == 0
    # step 10: solid-body rotation
    w0 = sp.symbols("omega0", positive=True)
    assert sp.simplify(sp.diff(r * w0 * r, r) / r) == 2 * w0
    # the code's four legs have the same orientation as the construction (numeric spot check)
    subsn = {**{ai: v for ai, v in zip(a, (0.2, -0.3, 0.5, 0.1, -0.2))}, **{bi: v for bi, v in zip(b, (0.4, 0.1, -0.3, 0.2, 0.6))}}
    urn = sp.lambdify((r, th), ur.subs(subsn))
    utn = sp.lambdify((r, th), ut.subs(subsn))
    legs = ch03.annular_sector_circulation(utn, 0.8, 0.3, 0.3, u_r=urn, theta0=0.2, legs=True)
    ex = {k: float(v.subs(subsn).subs({r0: 0.8, t0: 0.2, d: 0.3}))
          for k, v in (("radial_start", radial_start), ("inner", inner))}
    assert legs["radial_start"] == pytest.approx(ex["radial_start"], rel=1e-12)
    assert legs["inner"] == pytest.approx(ex["inner"], rel=1e-12)


# =====================================================================================================================
# C14 — Rankine (3.28) and Gaussian (3.29) vortices; N42–N44; D18, D19, D20
# =====================================================================================================================
def test_rankine_V1_continuity_peak_circulation():  # V1
    Gm, s = 2 * np.pi, 1.0
    r = np.linspace(1e-6, 5, 50001)
    u, w = ch03.rankine_vortex(r, Gm, s)
    k = np.argmax(u)
    assert abs(r[k] - s) < 2e-4 and u[k] == pytest.approx(Gm / (2 * np.pi * s), rel=1e-4)  # peak at σ: 1 m/s
    ui, _ = ch03.rankine_vortex(s * (1 - 1e-12), Gm, s)
    uo, _ = ch03.rankine_vortex(s * (1 + 1e-12), Gm, s)
    assert abs(ui - uo) < 1e-11  # continuous speed
    assert np.allclose(w[r < s], Gm / (np.pi * s ** 2)) and np.allclose(w[r > s], 0.0)  # vorticity jumps
    for rr in (1.0, 1.5, 3.0):
        assert ch03.circulation_circle("rankine", rr, Gamma=Gm, sigma=s) == pytest.approx(Gm, rel=1e-14)
    assert ch03.circulation_circle("rankine", 0.5, Gamma=Gm, sigma=s) == pytest.approx(Gm * 0.25, rel=1e-14)
    # worked number: Γ = 2π, σ = 1: inside u = r, ω = 2
    assert ch03.rankine_vortex(0.4, Gm, s) == pytest.approx((0.4, 2.0))


def test_gaussian_V1_circulation_limits_and_no_cancellation():  # V1
    Gm, s = 3.0, 0.8
    r = np.linspace(0.01, 6, 300)
    u, w = ch03.gaussian_vortex(r, Gm, s)
    assert np.allclose(2 * np.pi * r * u, Gm * (1 - np.exp(-r ** 2 / s ** 2)), rtol=1e-13)  # Γ(r) (D19 step 3)
    assert np.allclose(w, Gm / (np.pi * s ** 2) * np.exp(-r ** 2 / s ** 2), rtol=1e-14)
    # r ≪ σ: solid body with ω₀ = Γ/2πσ², computed without cancellation (series to 3 terms)
    for rr in (1e-9 * s, 1e-6 * s, 1e-3 * s):
        x = (rr / s) ** 2
        series = Gm / (2 * np.pi * s ** 2) * rr * (1 - x / 2 + x ** 2 / 6)
        assert abs(ch03.gaussian_vortex(rr, Gm, s)[0] / series - 1) < 1e-14
    assert ch03.gaussian_vortex(0.0, Gm, s)[0] == 0.0
    # r ≫ σ: line vortex Γ/2πr
    assert ch03.gaussian_vortex(20 * s, Gm, s)[0] == pytest.approx(Gm / (2 * np.pi * 20 * s), rel=1e-15)
    # worked numbers (Γ = 2π, σ = 1): u(1) = 1 - 1/e, u(5) ≈ 1/5
    assert ch03.gaussian_vortex(1.0, 2 * np.pi, 1.0)[0] == pytest.approx(1 - math.exp(-1), rel=1e-15)
    assert ch03.gaussian_vortex(5.0, 2 * np.pi, 1.0)[0] == pytest.approx(0.2 * (1 - math.exp(-25)), rel=1e-15)


def test_gaussian_max_V1_three_independent_routes():  # V1  N44, D20
    xb = ch03.gaussian_vortex_max_radius(1.0) ** 2
    xl = ch03.gaussian_vortex_max_radius(1.0, method="lambertw") ** 2
    assert abs(xb - xl) < 1e-13 and abs(xb - X_STAR) < 1e-13
    assert abs(1 + 2 * xb - math.exp(xb)) < 1e-13  # residual of 1 + 2x = e^x
    for s in (0.3, 1.0, 4.0):
        res = minimize_scalar(lambda q: -ch03.gaussian_vortex(q, 2 * np.pi, s)[0], bounds=(0.5 * s, 2 * s),
                              method="bounded", options={"xatol": 1e-12 * s})
        assert abs(res.x - ch03.gaussian_vortex_max_radius(s)) < 1e-6 * s  # numerical maximisation of u_θ
    r_star = ch03.gaussian_vortex_max_radius(1.0)
    assert r_star == pytest.approx(1.1209064228, abs=1e-10)
    assert abs(r_star - X_STAR) > 0.1  # r*/σ = √x*, not x* (the D20 trap)
    umax = ch03.gaussian_vortex(r_star, 2 * np.pi, 1.0)[0]
    assert umax == pytest.approx((1 - math.exp(-X_STAR)) / math.sqrt(X_STAR), rel=1e-14)  # 0.6382 Γ/2πσ
    assert abs(umax - 0.6382) < 5e-5
    with pytest.raises(ValueError):
        ch03.gaussian_vortex_max_radius(1.0, method="newton")


@needs_ref
def test_gaussian_max_V5_swirl_lamb_oseen_alpha():  # V5  Canivete Cuissa & Steiner (2022) Eq. (10)
    ref = ref_json()["swirl_lamb_oseen"]
    alpha = ref["alpha"]
    x_star = (ch03.gaussian_vortex_max_radius(1.0)) ** 2
    assert abs(x_star - alpha) <= ref["alpha_half_ulp"]  # 3-digit published value
    # their profile with r_max = √α σ and the normalisation (1 + 1/2α)(1 − e^{−α}) is our (3.29) peaking at v_max
    r, v_max, r_max, a = sp.symbols("r v_max r_max alpha", positive=True)
    v = sp.sympify(ref["velocity"], locals={"r": r, "v_max": v_max, "r_max": r_max, "alpha": a})
    vn = sp.lambdify((r, v_max, r_max, a), v)
    s = 0.7
    rr = np.linspace(0.05, 4, 60)
    ours = ch03.gaussian_vortex(rr, 2 * np.pi, s)[0]
    ours_max = ch03.gaussian_vortex(math.sqrt(x_star) * s, 2 * np.pi, s)[0]
    theirs = vn(rr, ours_max, math.sqrt(x_star) * s, x_star)
    assert np.max(np.abs(theirs / ours - 1)) < 1e-12  # identical profiles when α = x*
    assert abs((1 + 1 / (2 * x_star)) * (1 - math.exp(-x_star)) - 1) < 1e-14  # ⇔ e^α = 1 + 2α
    # the published α is rounded to 3 digits: the exact root of their normalisation lies inside α's rounding interval
    N = lambda al: (1 + 1 / (2 * al)) * (1 - math.exp(-al)) - 1  # noqa: E731
    lo, hi = alpha - ref["alpha_half_ulp"], alpha + ref["alpha_half_ulp"]
    assert N(lo) < 0 < N(hi) and lo <= x_star <= hi


@needs_ref
def test_rankine_V5_wikipedia_form():  # V5
    ref = ref_json()["rankine_vortex"]
    r, G, a = sp.symbols("r Gamma a", positive=True)
    loc = {"r": r, "Gamma": G, "a": a}
    vin = sp.lambdify((r, G, a), sp.sympify(ref["v_theta_inside"], locals=loc))
    vout = sp.lambdify((r, G, a), sp.sympify(ref["v_theta_outside"], locals=loc))
    Om = sp.lambdify((G, a), sp.sympify(ref["Omega_of_Gamma"], locals=loc))
    Gm, s = 5.0, 1.3
    rr = np.linspace(0.01, 5, 200)
    u, w = ch03.rankine_vortex(rr, Gm, s)
    assert np.allclose(u, np.where(rr <= s, vin(rr, Gm, s), vout(rr, Gm, s)), rtol=1e-14)
    assert np.allclose(w, np.where(rr <= s, 2 * Om(Gm, s), 0.0), rtol=1e-14)


@needs_ref
def test_gaussian_V5_lamb_oseen_form():  # V5  (σ² = 4νt, Ch. 5 pointer)
    ref = ref_json()["lamb_oseen_vortex"]
    r, G, nu, t = sp.symbols("r Gamma nu t", positive=True)
    loc = {"r": r, "Gamma": G, "nu": nu, "t": t}
    v = sp.lambdify((r, G, nu, t), sp.sympify(ref["v_theta"], locals=loc))
    w = sp.lambdify((r, G, nu, t), sp.sympify(ref["omega_z"], locals=loc))
    Gm, nun, tn = 0.02, 1.5e-5, 300.0
    s = math.sqrt(4 * nun * tn)
    rr = np.linspace(1e-3, 0.5, 100)
    u_o, w_o = ch03.gaussian_vortex(rr, Gm, s)
    assert np.allclose(u_o, v(rr, Gm, nun, tn), rtol=1e-12) and np.allclose(w_o, w(rr, Gm, nun, tn), rtol=1e-12)


def test_vortex_profile_V1_dispatcher():  # V1
    Gm, s = 2 * np.pi, 1.5
    r = np.array([0.3, 1.5, 3.0])
    us, ws = ch03.vortex_profile("solid", r, Gm, s)
    ur_, _ = ch03.rankine_vortex(r, Gm, s)
    assert us[1] == pytest.approx(ur_[1]) and np.allclose(ws, Gm / (np.pi * s ** 2))  # same speed at σ, ω = 2ω₀
    ul, wl = ch03.vortex_profile("line", r, Gm, s)
    assert np.allclose(ul, Gm / (2 * np.pi * r)) and np.all(wl == 0)
    assert np.allclose(ch03.vortex_profile("rankine", r, Gm, s)[0], ur_)
    assert np.allclose(ch03.vortex_profile("gaussian", r, Gm, s)[0], ch03.gaussian_vortex(r, Gm, s)[0])
    assert ch03.solid_body_rotation(2.0, 0.5) == 1.0 and ch03.line_vortex(2.0, 0.5) == 0.25
    assert math.isnan(ch03.line_vortex(0.0, 1.0))
    assert set(ch03.VORTEX_KINDS) == {"solid", "line", "rankine", "gaussian"}
    with pytest.raises(ValueError):
        ch03.vortex_profile("burgers", r)


def test_rankine_V2_derivation():  # V2  D18 steps 1-5 (and (3.23) applied to each piece of the code's profile)
    r, G, s = sp.symbols("r Gamma sigma", positive=True)
    ui, uo = G * r / (2 * sp.pi * s ** 2), G / (2 * sp.pi * r)
    wi, wo = ch03.polar_vorticity_z_sym(0, ui, r, sp.Symbol("theta")), ch03.polar_vorticity_z_sym(0, uo, r, sp.Symbol("theta"))
    assert sp.simplify(wi - G / (sp.pi * s ** 2)) == 0 and wo == 0  # steps 1-2
    assert sp.simplify(ui.subs(r, s) - uo.subs(r, s)) == 0  # step 3
    assert sp.simplify(2 * sp.pi * r * uo - G) == 0 and sp.simplify(wi * sp.pi * s ** 2 - G) == 0  # step 4
    assert sp.diff(ui, r).is_positive and sp.diff(uo, r).is_negative  # step 5: a kink maximum at σ
    rr = np.array([0.2, 0.7, 1.3, 2.4])
    code_w = ch03.rankine_vortex(rr, 4.0, 1.0)[1]
    assert np.allclose(code_w, np.where(rr <= 1.0, float(wi.subs({G: 4.0, s: 1.0})), 0.0))


def test_gaussian_V2_derivation():  # V2  D19 steps 1-7
    r, rp, G, s = sp.symbols("r r_p Gamma sigma", positive=True)
    w = G / (sp.pi * s ** 2) * sp.exp(-rp ** 2 / s ** 2)
    Gam_r = sp.integrate(w * 2 * sp.pi * rp, (rp, 0, r))  # steps 1-3
    assert sp.simplify(Gam_r - G * (1 - sp.exp(-r ** 2 / s ** 2))) == 0
    u = Gam_r / (2 * sp.pi * r)  # steps 4-5
    assert sp.simplify(u - G / (2 * sp.pi * r) * (1 - sp.exp(-r ** 2 / s ** 2))) == 0
    assert sp.simplify(sp.series(u, r, 0, 2).removeO() - G * r / (2 * sp.pi * s ** 2)) == 0  # step 6
    assert sp.limit(u * r, r, sp.oo) == G / (2 * sp.pi)  # step 7
    assert sp.simplify(ch03.polar_vorticity_z_sym(0, u, r, sp.Symbol("theta")) - w.subs(rp, r)) == 0  # (3.23) closes it


def test_gaussian_max_V2_derivation():  # V2  D20 steps 1-7
    x = sp.symbols("x", positive=True)
    f = (1 - sp.exp(-x)) / sp.sqrt(x)  # step 1
    fp = sp.diff(f, x)
    assert sp.simplify(fp - (sp.exp(-x) / sp.sqrt(x) - (1 - sp.exp(-x)) / (2 * x ** sp.Rational(3, 2)))) == 0  # step 2
    assert sp.simplify(fp * 2 * x ** sp.Rational(3, 2) - (2 * x * sp.exp(-x) - (1 - sp.exp(-x)))) == 0  # step 3
    assert sp.simplify(sp.exp(x) * (2 * x * sp.exp(-x) - (1 - sp.exp(-x))) - (1 + 2 * x - sp.exp(x))) == 0  # step 4
    g = 1 + 2 * x - sp.exp(x)
    assert g.subs(x, 0) == 0  # step 5: the trivial root
    assert float(g.subs(x, 0.5)) == pytest.approx(0.3513, abs=1e-4) and float(g.subs(x, 3)) == pytest.approx(-13.09, abs=1e-2)
    xs = sp.nsolve(g, x, 1.2)  # step 6
    assert abs(float(xs) - X_STAR) < 1e-14
    W = sp.LambertW(-sp.exp(-sp.Rational(1, 2)) / 2, -1)
    assert abs(float(-W - sp.Rational(1, 2)) - X_STAR) < 1e-14  # Lambert-W closed form
    assert abs(math.sqrt(X_STAR) - 1.12091) < 5e-6 and abs(float(f.subs(x, xs)) - 0.6382) < 5e-5  # step 7


# =====================================================================================================================
# C15 — Reynolds transport theorem (3.30)–(3.35), Leibniz, control volumes; N45–N55; D21–D24
# =====================================================================================================================
def F3(x, t):
    x_ = np.asarray(x, dtype=float)
    return (1 + 0.3 * x_[0] + 0.2 * x_[1] ** 2 + 0.1 * x_[0] * x_[2]) * (1 + 0.5 * t + 0.2 * t ** 2) + 0.05 * x_[2] ** 3 * np.sin(t)


def dF3(x, t):
    x_ = np.asarray(x, dtype=float)
    return (1 + 0.3 * x_[0] + 0.2 * x_[1] ** 2 + 0.1 * x_[0] * x_[2]) * (0.5 + 0.4 * t) + 0.05 * x_[2] ** 3 * np.cos(t)


def F2d(x, t):
    x_ = np.asarray(x, dtype=float)
    return (1 + 0.4 * x_[0] - 0.3 * x_[1] ** 2) * np.cos(0.5 * t) + 0.2 * x_[0] * x_[1] * t


def dF2d(x, t):
    x_ = np.asarray(x, dtype=float)
    return -0.5 * (1 + 0.4 * x_[0] - 0.3 * x_[1] ** 2) * np.sin(0.5 * t) + 0.2 * x_[0] * x_[1]


def _cvs():
    return {
        "box": ch03.MovingBox(lengths=(1.0, 0.8, 1.2), rates=(0.1, -0.05, 0.2), velocity=(0.3, 0.0, -0.1),
                              origin=(-0.5, -0.4, -0.6)),
        "sphere": ch03.GrowingSphere(R0=0.9, Rdot=0.15, center=(0.1, -0.2, 0.3), U=(0.2, 0.1, 0.0)),
        "cylinder": ch03.GrowingCylinder(R0=0.6, Rdot=0.1, L=1.0, Ldot=-0.2, base=(0.0, 0.0, -0.5), U=(0.1, 0.0, 0.2)),
        "cone": ch03.GrowingCone(r0=0.5, rdot=0.1, h=1.0),
    }


def test_rtt_V1_equals_definition_for_moving_shapes():  # V1 (two independent routes: RTT terms vs FD of the integral)
    for name, cv in _cvs().items():
        chk = ch03.rtt_check(F3, dF3, cv, 0.4)
        assert chk.abs_error < 1e-7 * max(1.0, abs(chk.rhs_total)), (name, chk)
        assert chk.rel_error < 1e-7
        terms = ch03.reynolds_transport(F3, dF3, cv, 0.4)
        assert terms.total == pytest.approx(terms.volume_term + terms.surface_term, abs=1e-15)
    # 2-D ellipse translating and deforming
    cv2 = ch03.MovingEllipse2D(a_fn=lambda t: 1.0 + 0.3 * t, b_fn=lambda t: 0.6 - 0.2 * t, center_fn=lambda t: (0.4 * t, 0.1))
    lhs = ch03.volume_integral_rate_fd(F2d, cv2, 0.5, dt=1e-4, n=32)
    rt = ch03.rtt_ellipse_2d(F2d, dF2d, 1.15, 0.5, 0.3, -0.2, c=(0.2, 0.1), cdot=(0.4, 0.0), t=0.5)
    assert abs(lhs - rt.total) < 1e-7
    assert abs(ch03.reynolds_transport(F2d, dF2d, cv2, 0.5, n=32).total - rt.total) < 1e-7


def test_rtt_V2_growing_sphere_closed_form_and_worked_number():  # V2
    t, R0, Rd = sp.symbols("t R_0 Rdot", positive=True)
    R = R0 + Rd * t
    I = t * 4 * sp.pi * R ** 5 / 15  # ∫_{ball R} t x₁² dV
    rate = sp.lambdify((t, R0, Rd), sp.diff(I, t))
    cv = ch03.GrowingSphere(R0=0.9, Rdot=0.1)
    Fx = lambda x, tt: tt * np.asarray(x)[0] ** 2  # noqa: E731
    dFx = lambda x, tt: np.asarray(x)[0] ** 2  # noqa: E731
    tot = ch03.reynolds_transport(Fx, dFx, cv, 1.0).total
    assert tot == pytest.approx(rate(1.0, 0.9, 0.1), rel=1e-13)
    # C15 worked number: R = 1 m, Ṙ = 0.1 m/s, F = t at t = 1 s: 4.189 + 1.257 = 5.445 = d/dt(t V)
    terms = ch03.reynolds_transport(lambda x, tt: tt + 0 * np.asarray(x)[0], lambda x, tt: 1 + 0 * np.asarray(x)[0], cv, 1.0)
    assert terms.volume_term == pytest.approx(4 / 3 * np.pi, rel=1e-13)
    assert terms.surface_term == pytest.approx(0.1 * 4 * np.pi, rel=1e-13)
    assert terms.total == pytest.approx(4 / 3 * np.pi + 0.4 * np.pi, rel=1e-13)


def test_rtt_V3_definition_converges_second_order_in_dt():  # V3
    cv = ch03.GrowingSphere(R0=0.9, Rdot=0.15, center=(0.1, 0.0, 0.0), U=(0.2, 0.0, 0.0))
    Fs = lambda x, t: (1 + 0.3 * np.asarray(x)[0]) * np.sin(1 + 2 * t)  # noqa: E731
    dFs = lambda x, t: 2 * (1 + 0.3 * np.asarray(x)[0]) * np.cos(1 + 2 * t)  # noqa: E731
    dts = [0.1, 0.05, 0.025, 0.0125]
    errs = [ch03.rtt_check(Fs, dFs, cv, 0.4, dt=dt).abs_error for dt in dts]
    assert abs(observed_order(dts, errs) - 2.0) < ORDER_TOL, pairwise_orders(dts, errs)


def test_rtt_V3_swept_terms_orders_of_smallness():  # V3  N50–N52, D22 steps 3, 5, 9, 10
    dts = [0.08, 0.04, 0.02, 0.01]
    rows = [ch03.swept_terms_sphere(0.9, 0.15, F3, dF3, 0.4, dt) for dt in dts]
    total = ch03.reynolds_transport(F3, dF3, ch03.GrowingSphere(R0=0.9, Rdot=0.15), 0.4).total
    for key, order in (("T4", 2.0), ("sliver_error", 2.0), ("taylor_residual", 2.0)):
        errs = [abs(r_[key]) for r_ in rows]
        assert abs(observed_order(dts, errs) - order) < ORDER_TOL, (key, pairwise_orders(dts, errs))
    errs = [abs(r_["lhs"] - total) for r_ in rows]
    assert abs(observed_order(dts, errs) - 1.0) < ORDER_TOL  # the one-sided quotient of (3.31) → total at order 1
    r0 = rows[-1]
    assert r0["T1"] + r0["T2"] + r0["T3"] + r0["T4"] + r0["taylor_residual"] == pytest.approx(r0["exact"], rel=1e-14)


def test_rtt_V1_fixed_volume_uniform_F_and_signed_sliver():  # V1  N48, N53, D22 step 1 (signed ΔV)
    fixed = ch03.MovingBox(lengths=(1, 2, 0.5), origin=(-0.5, 0, 0))
    terms = ch03.reynolds_transport(F3, dF3, fixed, 0.7)
    assert terms.surface_term == 0.0 and terms.total == terms.volume_term  # b = 0: d/dt passes inside
    for name, cv in _cvs().items():
        c = 2.5
        tot = ch03.reynolds_transport(lambda x, t: c + 0 * np.asarray(x)[0], lambda x, t: 0 * np.asarray(x)[0], cv, 0.3).total
        assert tot == pytest.approx(c * cv.volume_rate(0.3), rel=1e-12), name
        X, N, dA, b = cv.surface_nodes(0.3)
        assert np.sum(np.einsum("ik,ik->k", b, N) * dA) == pytest.approx(cv.volume_rate(0.3), rel=1e-12)
    # a shrinking balloon sweeps a negative sliver: the swept integral and the surface term are negative
    shrink = ch03.GrowingSphere(R0=1.0, Rdot=-0.2)
    one = lambda x, t: 1 + 0 * np.asarray(x)[0]  # noqa: E731
    assert ch03.surface_flux_term(one, shrink, 0.0) == pytest.approx(-0.2 * 4 * np.pi, rel=1e-13)
    sw = ch03.swept_volume_integral(one, shrink, 0.0, 0.01)
    assert sw == pytest.approx(4 / 3 * np.pi * (0.998 ** 3 - 1), rel=1e-12) and sw < 0


def test_control_volumes_V4_closed_surfaces_and_volumes():  # V4 (geometry: ∮n dA = 0, ∮x·n dA = d·V)
    shapes = dict(_cvs())
    shapes["ellipse"] = ch03.MovingEllipse2D(a_fn=lambda t: 1.0 + 0.3 * t, b_fn=lambda t: 0.6, center_fn=lambda t: (0.2, -0.1))
    one = lambda x, t: 1 + 0 * np.asarray(x)[0]  # noqa: E731
    for name, cv in shapes.items():
        X, N, dA, b = cv.surface_nodes(0.2, 24)
        d = X.shape[0]
        assert np.max(np.abs(N @ dA)) < 1e-12, name
        assert np.allclose(np.linalg.norm(N, axis=0), 1.0, atol=1e-14)
        assert np.sum(np.einsum("ik,ik->k", X, N) * dA) == pytest.approx(d * cv.volume(0.2), rel=1e-12), name
        assert ch03.volume_integral(one, cv, 0.2) == pytest.approx(cv.volume(0.2), rel=1e-12), name
        fd = (cv.volume(0.2 + 1e-5) - cv.volume(0.2 - 1e-5)) / 2e-5
        assert cv.volume_rate(0.2) == pytest.approx(fd, rel=1e-8), name
    assert ch03.volume_integral(one, ch03.GrowingSphere(R0=2.0), 0.0) == pytest.approx(4 / 3 * np.pi * 8, rel=1e-13)
    assert ch03.volume_integral(one, ch03.GrowingCone(r0=0.5, rdot=0.0, h=2.0), 0.0) == pytest.approx(np.pi * 0.25 * 2 / 3, rel=1e-13)
    assert isinstance(shapes["box"], ch03.ControlVolume)
    assert ch03.volume_rate_term(dF3, shapes["sphere"], 0.2) == pytest.approx(ch03.volume_integral(dF3, shapes["sphere"], 0.2))


def test_rtt_V1_one_dimensional_reduction_is_leibniz():  # V1  D22 step 12
    F1 = lambda x, t: np.sin(np.asarray(x)[0] - t) + 0.2 * np.asarray(x)[0] ** 2 * t  # noqa: E731  depends on x₁ only
    dF1 = lambda x, t: -np.cos(np.asarray(x)[0] - t) + 0.2 * np.asarray(x)[0] ** 2  # noqa: E731
    slab = ch03.MovingBox(lengths=(1.5, 1.0, 1.0), rates=(0.4, 0.0, 0.0), velocity=(-0.3, 0.0, 0.0), origin=(0.2, 0, 0))
    t = 0.6
    rtt = ch03.reynolds_transport(F1, dF1, slab, t).total
    a, b = 0.2 - 0.3 * t, 0.2 - 0.3 * t + 1.5 + 0.4 * t
    lz = ch03.leibniz_terms(lambda x, tt: F1(np.asarray(x)[None, ...], tt), lambda x, tt: dF1(np.asarray(x)[None, ...], tt),
                            a, b, -0.3, 0.1, t)
    assert rtt == pytest.approx(lz.total, rel=1e-12)


def test_leibniz_V1_cases_match_closed_form_rates():  # V1  N46 (3.30), D21 check numbers
    for case in ("x2t", "wave", "bump", "power"):
        for t in (0.7, 2.0):
            out = ch03.leibniz_example(t, case)
            assert out["total"] == pytest.approx(out["exact"], rel=1e-10, abs=1e-12), (case, t)
            assert out["total"] == pytest.approx(out["interior"] + out["upper"] - out["lower"], abs=1e-14)
            assert out["exact_rate"] == out["exact"]
    o = ch03.leibniz_example(2.0)
    assert (o["interior"], o["upper"], o["lower"]) == pytest.approx((56 / 3, 128.0, 8.0), rel=1e-13)
    assert o["total"] == pytest.approx((7 * 64 - 4 * 8) / 3, rel=1e-14)  # d/dt (t⁷ − t⁴)/3 at t = 2
    with pytest.raises(ValueError):
        ch03.leibniz_example(1.0, "nope")


def test_leibniz_V3_measured_rate_order_two():  # V3
    F = lambda x, t: np.sin(x - t)  # noqa: E731
    a_fn, b_fn = (lambda t: t / 2), (lambda t: 2 + t)
    t0 = 0.8
    exact = ch03.leibniz_example(t0, "wave")["exact"]
    dts = [0.2, 0.1, 0.05, 0.025]
    errs = [abs(ch03.leibniz_check(F, a_fn, b_fn, t0, dt=dt) - exact) for dt in dts]
    assert abs(observed_order(dts, errs) - 2.0) < ORDER_TOL, pairwise_orders(dts, errs)


@needs_ref
def test_leibniz_V5_boundary_term_signs():  # V5  Wikipedia "Leibniz integral rule"
    ref = ref_json()["leibniz_integral_rule"]
    F = lambda x, t: 1.0 + x ** 2 * t  # noqa: E731
    dF = lambda x, t: x ** 2 + 0 * x  # noqa: E731
    lz = ch03.leibniz_terms(F, dF, 0.5, 2.0, 0.3, -0.4, 1.2)
    assert lz.upper == pytest.approx(-0.4 * F(2.0, 1.2)) and lz.lower == pytest.approx(0.3 * F(0.5, 1.2))
    assert lz.total == pytest.approx(lz.interior + ref["upper_sign"] * lz.upper + ref["lower_sign"] * lz.lower, abs=1e-15)
    # the "+ lower" slip would disagree with the measured rate
    meas = ch03.leibniz_check(F, lambda t: 0.5 + 0.3 * (t - 1.2), lambda t: 2.0 - 0.4 * (t - 1.2), 1.2)
    assert meas == pytest.approx(lz.total, abs=1e-7) and abs(meas - (lz.interior + lz.upper + lz.lower)) > 0.1


@needs_ref
def test_rtt_V5_wikipedia_boundary_term_and_fixed_region():  # V5
    ref = ref_json()["reynolds_transport_theorem"]
    assert "v_b . n" in ref["statement"] and ref["normal"].startswith("outward")
    # sign of the boundary term: a growing sphere with F = 1 (outward b·n > 0) gains volume
    assert ch03.reynolds_transport(lambda x, t: 1 + 0 * np.asarray(x)[0], lambda x, t: 0 * np.asarray(x)[0],
                                   ch03.GrowingSphere(R0=1.0, Rdot=0.1), 0.0).total > 0
    fixed = ch03.GrowingSphere(R0=1.0, Rdot=0.0)
    lhs = ch03.volume_integral_rate_fd(F3, fixed, 0.3)
    assert lhs == pytest.approx(ch03.volume_integral(dF3, fixed, 0.3), rel=1e-7)  # v_b = 0: d/dt passes inside


def test_material_volume_rate_V4_gauss_and_divergence_limit():  # V4  N53, D23
    G = RNG.normal(size=(3, 3))
    lin = lambda x, t: G @ np.asarray(x)  # noqa: E731
    for cv in (ch03.GrowingSphere(R0=0.7, center=(0.2, 0.1, -0.3)), ch03.MovingBox(lengths=(1, 0.5, 2), origin=(0, 0, 0))):
        mv = ch03.material_volume_rate(lin, cv, 0.0)
        assert mv.surface_flux == pytest.approx(mv.volume_div, rel=1e-10)
        assert mv.surface_flux == pytest.approx(np.trace(G) * cv.volume(0.0), rel=1e-10)
    # nonlinear field: Gauss holds, and (1/δV)∮u·n dA → ∇·u at the centre at second order in the radius
    c = np.array([0.2, -0.1, 0.3])
    UD = SymField([U3.exprs[0] + sp.Rational(3, 10) * X1 ** 2 * X2, U3.exprs[1] + sp.sin(X2) * X3,
                   U3.exprs[2] + sp.Rational(1, 5) * X3 ** 2 * X1], (X1, X2, X3))  # div u varies in space
    mv = ch03.material_volume_rate(UD, ch03.GrowingSphere(R0=0.5, center=tuple(c)), 0.4)
    assert mv.surface_flux == pytest.approx(mv.volume_div, rel=1e-8)
    Rs = [0.2, 0.1, 0.05, 0.025]
    errs = []
    for R in Rs:
        sph = ch03.GrowingSphere(R0=R, center=tuple(c))
        errs.append(abs(ch03.material_volume_rate(UD, sph, 0.4).surface_flux / sph.volume(0.0) - float(UD.div(c, 0.4))))
    assert abs(observed_order(Rs, errs) - 2.0) < ORDER_TOL, pairwise_orders(Rs, errs)
    sol = ch03.material_volume_rate(lambda x, t: ch03.rigid_body_velocity([0.1, 0, 0], [0.3, -0.2, 1.0], x),
                                    ch03.GrowingCylinder(R0=0.5, L=1.0), 0.0)
    assert abs(sol.surface_flux) < 1e-12 and abs(sol.volume_div) < 1e-9  # incompressible: volume kept


def test_example_3_2_V1_three_routes_and_base_b_dot_n():  # V1  N54 (Ex. 3.2) with the corrected physics
    for h, r0, rd in ((1.0, 0.5, 0.1), (2.0, 0.3, -0.05), (0.7, 1.2, 0.4)):
        out = ch03.example_3_2(h, r0, rd)
        exact = 2 / 3 * np.pi * h * r0 * rd
        for k in ("direct", "rtt_dblquad", "rtt_cv", "rtt_closed_form"):
            assert out[k] == pytest.approx(exact, rel=1e-12), k
    assert ch03.example_3_2(1.0, 0.5, 0.1)["direct"] == pytest.approx(0.1047197551, rel=1e-9)  # worked number
    cone = ch03.GrowingCone(r0=0.5, rdot=0.1, h=1.0)
    X, N, dA, b = cone.surface_nodes(0.0, 16)
    base = np.isclose(X[2], 1.0) & np.isclose(N[2], 1.0)
    assert base.sum() > 0
    assert np.max(np.abs(np.einsum("ik,ik->k", b[:, base], N[:, base]))) < 1e-15  # b·n = 0 on the base …
    assert np.max(np.linalg.norm(b[:, base], axis=0)) > 0.05  # … although b ≠ 0 there (the book's "b = 0" is too strong)
    assert cone.volume_rate(0.0) == pytest.approx(ch03.example_3_2(1.0, 0.5, 0.1)["direct"], rel=1e-14)


def test_example_3_2_V2_symbolic_cone_integral():  # V2  N54 (a-D33), the "[?]" is the factor z
    z, ph, h, r0, rd = sp.symbols("z phi h r0 rdot", positive=True)
    th = sp.atan(r0 / h)
    integrand = (z / h) * rd * sp.cos(th) * z * sp.tan(th) / sp.cos(th)  # b·n · dA/(dφ dz)
    val = sp.integrate(integrand, (ph, 0, 2 * sp.pi), (z, 0, h))
    assert sp.simplify(val - sp.Rational(2, 3) * sp.pi * h ** 2 * rd * sp.tan(th)) == 0
    assert sp.simplify(val - sp.Rational(2, 3) * sp.pi * h * r0 * rd) == 0
    assert sp.simplify(sp.diff(sp.pi * h * (r0 + rd * sp.Symbol("t")) ** 2 / 3, sp.Symbol("t")).subs(sp.Symbol("t"), 0) - val) == 0
    # without the factor z the integral would be wrong (dimensionally m² instead of m³/s·…): the "[?]" is z
    wrong = sp.integrate((z / h) * rd * sp.cos(th) * sp.tan(th) / sp.cos(th), (ph, 0, 2 * sp.pi), (z, 0, h))
    assert sp.simplify(wrong - val) != 0


def test_rtt_field_V1_definitions_and_carried_pattern():  # V1  E7 fields
    x = np.array([[0.4, -1.0], [0.2, 0.3], [0.0, 1.0]])
    for name, (F_, dF_) in ((n, ch03.rtt_field(n)) for n in ch03.RTT_FIELDS):
        assert np.asarray(F_(x, 2.0)).shape == (2,)
    F, dF = ch03.rtt_field("warming")
    assert F(x, 2.0) == pytest.approx(1 + 0.5 * x[0] + 0.6) and np.allclose(dF(x, 2.0), 0.3)
    F, dF = ch03.rtt_field("carried", dim=2)
    xx = np.array([[0.4], [0.1]])
    h = 1e-6
    dFdx1 = (F(xx + [[h], [0]], 1.0) - F(xx - [[h], [0]], 1.0)) / (2 * h)
    assert abs(np.asarray(dF(xx, 1.0)).item() + ch03.CARRY_SPEED * np.asarray(dFdx1).item()) < 1e-9  # DF/Dt = 0 in u = (0.3, 0)
    rt = ch03.rtt_ellipse_2d(F, dF, 1.0, 0.6, 0.0, 0.0, cdot=(ch03.CARRY_SPEED, 0.0))
    assert abs(rt.total) < 1e-13 and abs(rt.volume_term) > 0.1  # the two terms cancel for a box riding with the pattern
    F1, dF1 = ch03.rtt_field("ramp", dim=1)
    assert ch03.leibniz_terms(F1, dF1, 0.0, 2.0, 0.0, 1.0, 0.0).total == pytest.approx(F1(2.0, 0.0))
    assert ch03.rtt_field("uniform")[0](x, 5.0) == pytest.approx(1.0)
    with pytest.raises(ValueError):
        ch03.rtt_field("nope")


def test_leibniz_V2_derivation():  # V2  D21 steps 1-7 for a general polynomial F and arbitrary moving limits
    x, xp, t, c = sp.symbols("x x_p t c", real=True)
    co = sp.symbols("k0:9", real=True)
    F = sum(co[3 * i + j] * x ** i * t ** j for i in range(3) for j in range(3))
    a, b = sp.Function("a")(t), sp.Function("b")(t)
    Phi = sp.integrate(F.subs(x, xp), (xp, c, x))  # step 1
    assert sp.simplify(sp.diff(Phi, x) - F) == 0
    I = Phi.subs(x, b) - Phi.subs(x, a)  # step 2
    assert sp.simplify(I - sp.integrate(F, (x, a, b))) == 0
    dPhidt = sp.integrate(sp.diff(F, t).subs(x, xp), (xp, c, x))  # step 5
    assert sp.simplify(sp.diff(Phi, t) - dPhidt) == 0
    assert sp.simplify(dPhidt.subs(x, b) - dPhidt.subs(x, a) - sp.integrate(sp.diff(F, t), (x, a, b))) == 0  # step 6
    rhs = sp.integrate(sp.diff(F, t), (x, a, b)) + sp.diff(b, t) * F.subs(x, b) - sp.diff(a, t) * F.subs(x, a)  # step 7
    assert sp.simplify(sp.diff(I, t) - rhs) == 0
    assert sp.simplify(sp.diff(I, t) - (rhs + 2 * sp.diff(a, t) * F.subs(x, a))) != 0  # the lower-limit sign matters


def test_rtt_V2_derivation():  # V2  D22 ★★★ (the design's sympy cell + the intermediate lines)
    t, dt, R0, Rd = sp.symbols("t Delta_t R_0 Rdot", positive=True)
    r, th, ph = sp.symbols("r theta phi", nonnegative=True)
    R = R0 + Rd * t
    F = t * (r * sp.sin(th) * sp.cos(ph)) ** 2
    dV = r ** 2 * sp.sin(th)

    def vol(expr, ra, rb):
        return sp.integrate(expr * dV, (r, ra, rb), (th, 0, sp.pi), (ph, 0, 2 * sp.pi))
    lhs = sp.diff(vol(F, 0, R), t)
    vterm = vol(sp.diff(F, t), 0, R)  # step 7
    sterm = sp.integrate((F * Rd * dV).subs(r, R), (th, 0, sp.pi), (ph, 0, 2 * sp.pi))  # step 10
    assert sp.simplify(lhs - (vterm + sterm)) == 0  # step 11: (3.35)
    T1 = vol(F, 0, R)
    T2 = vol(dt * sp.diff(F, t), 0, R)
    T3 = vol(F, R, R.subs(t, t + dt))  # step 1: the swept shell
    T4 = vol(dt * sp.diff(F, t), R, R.subs(t, t + dt))
    assert sp.series(T4, dt, 0, 2).removeO() == 0  # step 5: T4 = O(Δt²)
    assert sp.simplify(sp.limit(T3 / dt, dt, 0) - sterm) == 0  # steps 9-10: (3.34)
    exact_new = vol(F.subs(t, t + dt), 0, R.subs(t, t + dt))
    # step 3: (3.32) is exact here up to O(Δt²) (F is linear in t, so the Taylor step of F is exact)
    assert sp.simplify(sp.expand(exact_new - (T1 + T2 + T3 + T4))) == 0
    # step 6: (3.33) — the bracket of (3.31) divided by Δt tends to T2/Δt + T3/Δt
    assert sp.simplify(sp.limit((exact_new - T1) / dt, dt, 0) - (vterm + sterm)) == 0
    # step 8: a patch dA moving at b sweeps (bΔt)·n dA — sliding along the surface sweeps nothing
    bvec, nvec = sp.Matrix(sp.symbols("b1:4", real=True)), sp.Matrix([0, 0, 1])
    tang = sp.Matrix([1, 0, 0])
    assert (bvec + 5 * tang).dot(nvec) == bvec.dot(nvec)
    # step 12 (the 1-D reduction) is tested numerically in test_rtt_V1_one_dimensional_reduction_is_leibniz


def test_material_volume_V2_derivation():  # V2  D23 steps 1-5 (2-D material disc, symbolic polynomial u)
    x, y, eps_, th = sp.symbols("x y epsilon theta", positive=True)
    c = sp.symbols("c0:12", real=True)
    u = c[0] + c[1] * x + c[2] * y + c[3] * x ** 2 + c[4] * x * y + c[5] * y ** 2
    v = c[6] + c[7] * x + c[8] * y + c[9] * x ** 2 + c[10] * x * y + c[11] * y ** 2
    on = {x: eps_ * sp.cos(th), y: eps_ * sp.sin(th)}
    flux = sp.integrate(((u * sp.cos(th) + v * sp.sin(th)).subs(on)) * eps_, (th, 0, 2 * sp.pi))  # ∮u·n ds (step 1)
    rr = sp.symbols("rr", positive=True)
    div = sp.diff(u, x) + sp.diff(v, y)
    vol_div = sp.integrate(sp.integrate(div.subs({x: rr * sp.cos(th), y: rr * sp.sin(th)}) * rr, (th, 0, 2 * sp.pi)),
                           (rr, 0, eps_))
    assert sp.simplify(flux - vol_div) == 0  # step 2: Gauss
    rate = sp.limit(flux / (sp.pi * eps_ ** 2), eps_, 0)  # steps 3-5
    assert sp.simplify(rate - div.subs({x: 0, y: 0})) == 0


def test_material_derivative_from_rtt_V2_derivation():  # V2  D24 steps 1-7: the hidden F∇·u term
    x, y, eps_, th, t = sp.symbols("x y epsilon theta t", real=True)
    eps_ = sp.symbols("epsilon", positive=True)
    c = sp.symbols("c0:6", real=True)
    k = sp.symbols("k0:5", real=True)
    u = c[0] + c[1] * x + c[2] * y
    v = c[3] + c[4] * x + c[5] * y + x * y
    F = k[0] + k[1] * x + k[2] * y + k[3] * x * y + k[4] * t * x ** 2
    on = {x: eps_ * sp.cos(th), y: eps_ * sp.sin(th)}
    rr = sp.symbols("rr", positive=True)
    polar = {x: rr * sp.cos(th), y: rr * sp.sin(th)}
    surf = sp.integrate((F * (u * sp.cos(th) + v * sp.sin(th))).subs(on) * eps_, (th, 0, 2 * sp.pi))  # ∮F u·n ds
    vol = sp.integrate(sp.integrate(sp.diff(F, t).subs(polar) * rr, (th, 0, 2 * sp.pi)), (rr, 0, eps_))
    per_area = sp.limit((surf + vol) / (sp.pi * eps_ ** 2), eps_, 0)  # (1/δV) d/dt ∫F dV for a small material disc
    at0 = {x: 0, y: 0}
    DFDt = (sp.diff(F, t) + u * sp.diff(F, x) + v * sp.diff(F, y)).subs(at0)
    divu = (sp.diff(u, x) + sp.diff(v, y)).subs(at0)
    assert sp.simplify(per_area - (DFDt + F.subs(at0) * divu)) == 0  # step 4: DF/Dt + F ∇·u
    assert sp.simplify(per_area - DFDt) != 0  # the book's "extension of (3.5)" hides F∇·u
    # step 3: product rule for the divergence, with generic functions
    f = sp.Function("f")(x, y)
    p, q = sp.Function("p")(x, y), sp.Function("q")(x, y)
    lhs = sp.diff(f * p, x) + sp.diff(f * q, y)
    assert sp.simplify(lhs - (p * sp.diff(f, x) + q * sp.diff(f, y) + f * (sp.diff(p, x) + sp.diff(q, y)))) == 0


# =====================================================================================================================
# §3.1 — N03 section averages, N05 coordinate systems
# =====================================================================================================================
def test_section_average_V1_poiseuille_and_developing_profile():  # V1  N03
    U, R = 2.0, 0.3
    assert ch03.cross_section_average(lambda r, z: U * (1 - r ** 2 / R ** 2), R) == pytest.approx(U / 2, rel=1e-12)
    assert ch03.cross_section_average(lambda r, z: U + 0 * r, R) == pytest.approx(U, rel=1e-12)
    for z in (0.0, 0.3, 1.5, 10.0):  # the developing profile keeps its mean (mass conservation) at every station
        assert ch03.cross_section_average(lambda r, zz: ch03.pipe_profile(r, zz, 1.3, 0.5), 0.5, z) == pytest.approx(1.3, rel=1e-12)
    assert ch03.pipe_profile(0.0, 50.0, 1.0, 1.0) == pytest.approx(2.0, rel=1e-12)  # Poiseuille far downstream: u_max = 2ū


def test_coordinates_V1_round_trips_and_ranges():  # V1  N05
    P = RNG.uniform(-5, 5, size=(3, 10_000))
    R, phi, z = ch03.cylindrical_from_cartesian(*P)
    assert np.max(np.abs(np.stack(ch03.cartesian_from_cylindrical(R, phi, z)) - P)) < 1e-13 * 10
    r, th, ph = ch03.spherical_from_cartesian(*P)
    assert np.max(np.abs(np.stack(ch03.cartesian_from_spherical(r, th, ph)) - P)) < 1e-13 * 10
    assert th.min() >= 0 and th.max() <= np.pi and ph.min() > -np.pi - 1e-15 and ph.max() <= np.pi
    rp, tp = ch03.polar_from_cartesian(P[0], P[1])
    assert np.max(np.abs(np.stack(ch03.cartesian_from_polar_coords(rp, tp)) - P[:2])) < 1e-13 * 10
    # θ from +z: 0 on +z, π on −z, π/2 in the plane; the axis R = 0 is finite
    assert ch03.spherical_from_cartesian(0, 0, 2.0)[1] == 0.0 and ch03.spherical_from_cartesian(0, 0, -2.0)[1] == pytest.approx(np.pi)
    assert ch03.spherical_from_cartesian(1.0, 1.0, 0.0)[1] == pytest.approx(np.pi / 2)
    assert ch03.cylindrical_from_cartesian(0.0, 0.0, 1.0) == (0.0, 0.0, 1.0)


def test_coordinates_V1_unit_vectors_and_projections():  # V1  N05
    for _ in range(20):
        phi = RNG.uniform(-np.pi, np.pi)
        th = RNG.uniform(0, np.pi)
        Ec = ch03.unit_vectors_cylindrical(phi)
        Es = ch03.unit_vectors_spherical(th, phi)
        for E in (Ec, Es):
            assert np.allclose(E @ E.T, np.eye(3), atol=1e-14) and np.linalg.det(E) == pytest.approx(1.0, abs=1e-14)
        assert np.allclose(np.cross(Ec[0], Ec[1]), Ec[2], atol=1e-14) and np.allclose(np.cross(Es[0], Es[1]), Es[2], atol=1e-14)
        Ep = ch03.unit_vectors_polar(phi)
        assert np.allclose(Ep @ Ep.T, np.eye(2), atol=1e-14)
    assert ch03.unit_vectors_cylindrical(np.linspace(0, 1, 5)).shape == (3, 3, 5)
    X = RNG.uniform(-2, 2, size=(3, 50))
    Ucart = RNG.normal(size=(3, 50))
    for system in ("cylindrical", "spherical"):
        uc = ch03.velocity_components(Ucart, X, system)
        assert np.allclose(np.linalg.norm(uc, axis=0), np.linalg.norm(Ucart, axis=0), rtol=1e-13)
        assert np.allclose(ch03.cartesian_components(uc, X, system), Ucart, atol=1e-13)
    up = ch03.velocity_components(Ucart[:2], X[:2], "polar")
    assert np.allclose(np.linalg.norm(up, axis=0), np.linalg.norm(Ucart[:2], axis=0), rtol=1e-13)
    # solid body u = Ω e_z × x: u_φ = ΩR, u_R = u_z = 0
    Om = 0.7
    usb = ch03.rigid_body_velocity([0, 0, 0], [0, 0, Om], X)
    ucyl = ch03.velocity_components(usb, X, "cylindrical")
    assert np.allclose(ucyl[1], Om * np.hypot(X[0], X[1]), atol=1e-14) and np.allclose(ucyl[[0, 2]], 0, atol=1e-14)
    with pytest.raises(ValueError):
        ch03.velocity_components(Ucart, X, "toroidal")


# =====================================================================================================================
# V6 — the book's own statements (private JSON; skipped when absent). The report quotes relative errors only.
# =====================================================================================================================
def _bexpr(s, **extra):
    names = dict(t=sp.Symbol("t"), tp=sp.Symbol("tp"), t0=sp.Symbol("t0"), omega=sp.Symbol("omega"), xi0=sp.Symbol("xi0"),
                 h=sp.Symbol("h"), r0=sp.Symbol("r0"), rdot=sp.Symbol("rdot"), theta=sp.Symbol("theta"),
                 Gamma=sp.Symbol("Gamma"), sigma=sp.Symbol("sigma"), B=sp.Symbol("B"), omega0=sp.Symbol("omega0"),
                 gamma=sp.Symbol("gamma"), r=sp.Symbol("r"), x=sp.Symbol("x"), y=sp.Symbol("y"), z=sp.Symbol("z"))
    names.update(extra)
    return sp.sympify(s, locals=names)


@book_only
def test_example_3_1_V6_book_closed_forms():  # V6
    bk = book()["example_3_1"]
    S = sp.Symbol
    for ang in bk["sample_omega_tp_rad"]:
        w, xi = 1.3, 0.8
        tp = ang / w
        subs = {S("omega"): w, S("xi0"): xi, S("tp"): tp}
        ex = ch03.example_3_1(tp, xi0=xi, omega=w)
        pc = [float(_bexpr(e).subs(subs)) for e in bk["pathline_circle_centre"]]
        sc = [float(_bexpr(e).subs(subs)) for e in bk["streakline_circle_centre"]]
        assert np.allclose(ex["path_center"], pc, atol=1e-14) and np.allclose(ex["streak_center"], sc, atol=1e-14)
        assert float(_bexpr(bk["pathline_radius"]).subs(subs)) == ex["radius"]
        te = np.linspace(tp, tp + 4, 9)
        pl = ch03.pathline(ex["field"], [0.0, 0.0], tp, te)
        bx = sp.lambdify(S("t"), _bexpr(bk["pathline_x"]).subs(subs))
        by = sp.lambdify(S("t"), _bexpr(bk["pathline_y"]).subs(subs))
        assert np.max(np.abs(pl - np.stack([bx(te), by(te)]))) < 1e-9
        tr = np.linspace(tp - 5, tp, 9)
        sk = ch03.streakline(ex["field"], [0.0, 0.0], tp, tr)
        sx = sp.lambdify(S("t0"), _bexpr(bk["streakline_x_param"]).subs(subs))
        sy = sp.lambdify(S("t0"), _bexpr(bk["streakline_y_param"]).subs(subs))
        assert np.max(np.abs(sk - np.stack([sx(tr), sy(tr)]))) < 1e-7
        if abs(math.cos(ang)) > 1e-12:
            assert abs(ex["slope"] - math.tan(ang)) < 1e-12  # streamline y = x tan(ωt')
    assert bk["all_tangent_at_origin"] is True


@book_only
def test_parallel_shear_V6_book_statements():  # V6
    bk = book()["shear_flow_3_5"]
    for gam in (0.5, 1.0, 2.5):
        k = ch03.parallel_shear_kinematics(gam)
        g = {sp.Symbol("gamma"): gam}
        assert k["omega3"] == pytest.approx(float(_bexpr(bk["omega3"]).subs(g)))
        assert k["spin"] == pytest.approx(float(_bexpr(bk["average_angular_velocity"]).subs(g)))
        assert np.allclose(k["S"], [[float(_bexpr(e).subs(g)) for e in row] for row in bk["S"]])
        assert np.allclose(k["S_bar"], [[float(_bexpr(e).subs(g)) for e in row] for row in bk["S_principal"]])
        assert k["principal_angle_deg"] == pytest.approx(bk["principal_angle_deg"])


@book_only
def test_vortices_V6_book_statements():  # V6  (3.22)–(3.29), Exercise 3.26 statement
    bk = book()
    sb = bk["solid_body_3_22_3_24"]
    w0, r = 0.7, 1.3
    sub = {sp.Symbol("omega0"): w0, sp.Symbol("r"): r}
    assert ch03.polar_vorticity_z(None, lambda q: ch03.solid_body_rotation(q, w0), r) == pytest.approx(float(_bexpr(sb["omega_z"]).subs(sub)), abs=1e-9)
    assert ch03.circulation_circle(lambda q: ch03.solid_body_rotation(q, w0), r) == pytest.approx(float(_bexpr(sb["circulation_radius_r"]).subs(sub)), rel=1e-14)
    lv = bk["line_vortex_3_25_3_27"]
    B = 0.9
    subB = {sp.Symbol("B"): B, sp.Symbol("r"): r}
    assert abs(ch03.polar_vorticity_z(None, lambda q: ch03.line_vortex(q, B), r)) < 1e-9 and _bexpr(lv["omega_z_away_from_origin"]) == 0
    assert ch03.circulation_circle(lambda q: ch03.line_vortex(q, B), r) == pytest.approx(float(_bexpr(lv["circulation_about_origin"]).subs(subB)), rel=1e-14)
    assert ch03.mean_vorticity_in_disc(lambda q: ch03.line_vortex(q, B), r) == pytest.approx(float(_bexpr(lv["mean_vorticity_disc_radius_r"]).subs(subB)), rel=1e-14)
    assert abs(ch03.annular_sector_circulation(lambda q: ch03.line_vortex(q, B), 0.5, 0.2, 0.3)) < 1e-12 and _bexpr(lv["circulation_ABCD"]) == 0
    vm = bk["vortex_models_3_28_3_29"]
    rel = abs(ch03.gaussian_vortex_max_radius(1.0) - vm["gaussian_umax_radius_over_sigma"]) / vm["gaussian_umax_radius_over_sigma"]
    assert rel < 5e-6  # the book quotes 6 significant figures; ours differs by rounding only
    rr = np.linspace(0.01, 3, 30001)
    assert rr[np.argmax(ch03.rankine_vortex(rr, 1.0, 1.0)[0])] == pytest.approx(vm["rankine_umax_radius_over_sigma"], abs=1e-4)


@book_only
def test_example_3_2_V6_book_forms():  # V6
    bk = book()["example_3_2"]
    s = bk["sample"]
    sub = {sp.Symbol("h"): s["h"], sp.Symbol("r0"): s["r0"], sp.Symbol("rdot"): s["rdot"],
           sp.Symbol("theta"): math.atan(s["r0"] / s["h"])}
    out = ch03.example_3_2(s["h"], s["r0"], s["rdot"])
    for key in ("dVdt_direct", "dVdt_rtt_intermediate"):
        book_v = float(_bexpr(bk[key]).subs(sub))
        assert abs(out["rtt_cv"] - book_v) / abs(book_v) < 1e-12
    assert float(_bexpr(bk["tan_theta"]).subs(sub)) == pytest.approx(math.tan(sub[sp.Symbol("theta")]))


# =====================================================================================================================
# Design Part C — the contract: every name exists with its signature and is exercised in this file
# =====================================================================================================================
CONTRACT = {
    # C.1 core.kinematics
    "lagrangian_velocity_acceleration": ["r_of_t", "t", "h"], "lagrangian_velocity_acceleration_sym": ["r_exprs", "t"],
    "lagrangian_to_eulerian": ["r_exprs", "labels", "coords", "t"], "material_derivative": ["F", "u", "x", "t", "h", "ht"],
    "material_derivative_terms": ["F", "u", "x", "t", "h"], "material_derivative_sym": ["F_expr", "u_exprs", "coords", "t"],
    "streamwise_derivative": ["F", "u", "x", "t", "h"],
    "streamline": ["u", "x0", "t_frozen", "s_max", "both", "n", "rtol"], "pathline": ["u", "r0", "t0", "t_eval", "rtol", "atol"],
    "streakline": ["u", "x0", "t", "t_release"], "galilean_transform": ["u", "U", "x0p"], "acceleration": ["u", "x", "t", "h"],
    "velocity_gradient_at": ["u", "x", "t", "h"], "relative_velocity": ["G", "dx"], "vorticity_from_gradient": ["G"],
    "vorticity": ["u", "x", "t", "h"], "linear_strain_rate": ["G", "n"], "shear_strain_rate": ["G", "n1", "n2"],
    "volumetric_strain_rate": ["G"], "material_volume_ratio": ["G", "t"], "element_rotation_rate": ["G"],
    "material_line_rotation_rate": ["G", "theta"], "vorticity_in_rotating_frame": ["omega", "Omega"],
    "velocity_potential_2d": ["u", "grid", "x_ref", "t"], "relative_velocity_split": ["G", "dx"],
    "principal_strain_rates": ["G"], "strain_velocity_principal": ["G", "dx"], "deform_circle": ["G", "t", "n", "radius"],
    "deform_sphere": ["G", "t", "n_theta", "n_phi", "radius"], "strain_ellipse_axes": ["G", "t", "method"],
    "measured_strain_rates": ["G", "n", "dt"], "velocity_gradient_preset": ["name", "Gamma", "dim"],
    "linear_flow_map": ["G", "t"], "deform_square": ["G", "t", "n_side", "half_width", "boundary_only"],
    "material_line_angle": ["G", "t", "theta0"], "perpendicular_pair_rotation_rate": ["G", "theta"],
    # C.2 core.coords
    "cylindrical_from_cartesian": ["x", "y", "z"], "cartesian_from_cylindrical": ["R", "phi", "z"],
    "spherical_from_cartesian": ["x", "y", "z"], "cartesian_from_spherical": ["r", "theta", "phi"],
    "unit_vectors_cylindrical": ["phi"], "unit_vectors_spherical": ["theta", "phi"],
    "velocity_components": ["u_cart", "x", "system"], "unit_vectors_polar": ["theta"], "polar_from_cartesian": ["x", "y"],
    "cartesian_components": ["u_curv", "x", "system"],
    # C.3 core.vortices
    "solid_body_rotation": ["r", "omega0"], "line_vortex": ["r", "B"], "rankine_vortex": ["r", "Gamma", "sigma"],
    "gaussian_vortex": ["r", "Gamma", "sigma"], "gaussian_vortex_max_radius": ["sigma", "method"],
    "vortex_profile": ["kind", "r", "Gamma", "sigma"], "vortex_velocity_field": ["profile"],
    "polar_vorticity_z": ["u_r", "u_theta", "r", "theta", "h"], "polar_vorticity_z_sym": ["u_r_expr", "u_theta_expr", "r", "theta"],
    "circulation_circle": ["u_theta", "r", "center", "n"], "mean_vorticity_in_disc": ["u_theta", "r"],
    # C.4 core.transport
    "leibniz_terms": ["F", "dFdt", "a", "b", "dadt", "dbdt", "t", "n"], "leibniz_check": ["F", "a_fn", "b_fn", "t", "dt"],
    "GrowingSphere": ["R0", "Rdot", "center"], "GrowingCylinder": ["R0", "Rdot", "L"],
    "MovingBox": ["lengths", "rates", "velocity", "origin"], "GrowingCone": ["r0", "rdot", "h"],
    "MovingEllipse2D": ["a_fn", "b_fn", "center_fn"],
    "volume_integral": ["F", "cv", "t", "n"], "volume_integral_rate_fd": ["F", "cv", "t", "dt"],
    "swept_volume_integral": ["F", "cv", "t", "dt"], "surface_flux_term": ["F", "cv", "t"],
    "reynolds_transport": ["F", "dFdt", "cv", "t"], "rtt_check": ["F", "dFdt", "cv", "t", "dt"],
    "material_volume_rate": ["u", "cv", "t"], "swept_terms": ["F", "dFdt", "cv", "t", "dt"],
    "rtt_ellipse_2d": ["F", "dFdt", "a", "b", "adot", "bdot", "c", "cdot", "t", "n"],
    "swept_terms_sphere": ["R", "Rdot", "F", "dFdt", "t", "dt"],
    # C.5 ch03_kinematics
    "cross_section_average": ["u_fn", "R", "z"], "cylinder_flow": ["x", "y", "U", "a", "frame"],
    "frame_acceleration_terms": ["x", "y", "U", "a", "U_frame", "t", "h"], "lagrangian_map_example": ["X", "t", "alpha"],
    "streamline_slope": ["u", "x", "t"], "example_3_1": ["t_prime", "xi0", "omega", "n"],
    "unsteady_flow_preset": ["name", "x", "y", "t"], "preset_field": ["name"],
    "thermal_front": ["x", "y", "t", "grad_K_per_m", "heating_K_per_s", "front_speed", "width_m", "T0"],
    "thermal_front_terms": ["x", "y", "t", "u", "v", "grad_K_per_m", "heating_K_per_s", "front_speed", "width_m"],
    "rigid_body_velocity": ["U", "Omega", "x"], "potential_velocity": ["phi_expr", "coords"],
    "parallel_shear_kinematics": ["gamma"], "annular_sector_circulation": ["u_theta", "r", "dr", "dtheta", "u_r"],
    "example_3_2": ["h", "r0", "rdot"], "leibniz_example": ["t", "case"], "rtt_field": ["name", "dim"],
    "flux_through_disc": ["u", "center", "normal", "radius", "t", "nr", "ntheta"], "pipe_profile": ["r", "z", "U_mean", "R"],
    "cylinder_velocity_field": ["U", "a", "frame"], "cylinder_streamfunction": ["x", "y", "U", "a"],
    # C.0 reused (called by the notebook / explainers)
    "strain_rate_tensor": ["G"], "rotation_tensor": ["G"], "antisymmetric_part": [], "antisymmetric_from_vector": [],
    "vector_from_antisymmetric": [], "principal_axes": [], "trace": [], "levi_civita": [], "cross": [],
    "transform_tensor": [], "rotation_matrix_2d": [], "curl": [], "divergence": [], "vector_gradient": [],
    "is_irrotational": [], "circulation": ["u_fn", "loop"], "planar_loop": ["center", "normal", "radius", "n"],
    "rectangle_loop": [], "planar_disc": [], "divergence_theorem_sphere": [], "expand_indices_str": [],
    "as_coord_field": ["u", "t"], "from_coord_field": ["fn"],
}


def test_contract_V1_every_part_c_name_exists_with_signature():  # V1 (contract)
    missing, bad = [], []
    for name, params in CONTRACT.items():
        obj = getattr(ch03, name, None)
        if obj is None:
            missing.append(name)
            continue
        if params:
            got = list(inspect.signature(obj).parameters)
            if got[:len(params)] != params:
                bad.append((name, got[:len(params)], params))
    assert not missing, missing
    assert not bad, bad
    # the core modules own them (re-exported, not copied)
    assert ch03.streamline is K.streamline and ch03.rankine_vortex is Xmod.rankine_vortex
    assert ch03.reynolds_transport is Rmod.reynolds_transport and ch03.velocity_components is Cmod.velocity_components
    assert ch02.velocity_gradient_preset is K.velocity_gradient_preset  # promoted with a re-export (analysis §4)


def test_contract_V1_every_part_c_name_is_exercised_in_this_file():  # V1 (self-check of coverage)
    src = Path(__file__).read_text(encoding="utf-8")
    body = src.split("CONTRACT = {")[0]  # usage before the contract table
    unused = [n for n in CONTRACT if f"ch03.{n}" not in body]
    extra_smoke = {"antisymmetric_part", "principal_axes", "trace", "cross", "transform_tensor", "rotation_matrix_2d",
                   "curl", "divergence", "vector_gradient", "is_irrotational", "rectangle_loop", "planar_disc",
                   "divergence_theorem_sphere", "expand_indices_str", "swept_terms", "unsteady_flow_field"}
    assert set(unused) <= extra_smoke, sorted(set(unused) - extra_smoke)


def test_contract_V1_reused_names_smoke():  # V1 smoke with a physical assertion for each reused callable
    G = RNG.normal(size=(3, 3))
    assert np.allclose(ch03.antisymmetric_part(G), 0.5 * ch03.rotation_tensor(G))  # A = ½R: the factor-2 convention
    lam, C = ch03.principal_axes(ch03.strain_rate_tensor(G))
    assert np.isclose(ch03.trace(ch03.strain_rate_tensor(G)), lam.sum())
    assert np.allclose(ch03.cross([1, 0, 0], [0, 1, 0]), [0, 0, 1])
    Q = rand_rotation()
    assert np.isclose(np.trace(ch03.transform_tensor(G, Q)), np.trace(G))
    assert np.allclose(ch03.rotation_matrix_2d(np.pi / 2) @ [1, 0], [0, 1], atol=1e-15)
    xs = np.linspace(-1, 1, 41)
    Xg, Yg = np.meshgrid(xs, xs, indexing="xy")
    usb = np.stack([-0.7 * Yg, 0.7 * Xg])  # solid body on a grid
    h = xs[1] - xs[0]
    assert np.allclose(ch03.curl(usb, h), 1.4, atol=1e-12) and np.allclose(ch03.divergence(usb, h), 0, atol=1e-12)
    assert ch03.vector_gradient(usb, h).shape[:2] == (2, 2)
    ok, resid = ch03.is_irrotational(usb, h)
    assert (not ok) and resid == pytest.approx(1.4, abs=1e-12)
    ok, resid = ch03.is_irrotational(np.stack([Xg, -Yg]), h)  # pure strain: irrotational
    assert ok and resid < 1e-12
    loop = ch03.rectangle_loop((0.0, 0.0), a=1.0, b=2.0)
    circ = ch03.circulation(lambda X, Y: np.stack([-0.7 * Y, 0.7 * X]), loop)
    assert circ == pytest.approx(1.4 * 2.0, rel=1e-12)  # Γ = ω × area
    disc = ch03.planar_disc((0.0, 0.0, 0.0), (0, 0, 1), radius=1.0)
    assert disc.area == pytest.approx(np.pi, rel=1e-12)
    out = ch03.divergence_theorem_sphere(lambda X, Y, Z: np.stack([X, Y, Z]), 1.0)
    lhs, rhs = out  # ∭ ∇·x dV = ∯ x·n dA = 3 V = 4π for the unit ball
    assert lhs == pytest.approx(4 * np.pi, rel=1e-6) and rhs == pytest.approx(4 * np.pi, rel=1e-12)
    assert ch03.expand_indices_str("u_i F,i").count("+") == 2
    sw = ch03.swept_terms(F3, dF3, ch03.GrowingSphere(R0=0.9, Rdot=0.15), 0.4, 1e-3)
    assert set(sw) >= {"T1", "T2", "T3", "T4", "exact", "lhs", "sliver", "sliver_error", "taylor_residual"}
    assert ch03.unsteady_flow_field is ch03.preset_field
    fn = ch03.as_coord_field(U2, 0.2)
    back = ch03.from_coord_field(fn)
    P = rand_points(2, 5)
    assert np.allclose(back(P), U2(P, 0.2))


def test_scripts_V1_drawing_helpers_run():  # V1 smoke  (C.6)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    from scripts.ch03_drawings import paddle, ring_arrows  # noqa: E402
    from scripts.ch03_fig3_7_flow_lines import flow_lines_figure  # noqa: E402
    from scripts.ch03_fig3_14_shear_elements import shear_elements_frames  # noqa: E402
    from scripts.ch03_fig3_17_leibniz import leibniz_strips_figure  # noqa: E402
    from scripts.ch03_fig3_18_rtt import rtt_blob_figure  # noqa: E402
    fig, ax = plt.subplots()
    ring_arrows(ax, ch03.velocity_gradient_preset("simple_shear", 1.0), parts=("total", "strain", "rotation"))
    paddle(ax, (0.0, 0.0), 0.3)
    plt.close(fig)
    out = flow_lines_figure(0.5)
    plt.close("all")
    frames = shear_elements_frames(1.0, [0.0, 0.3])
    assert len(frames) == 2 and set(frames[0]) >= {"ABCD", "PQRS"}
    # the aligned element ABCD keeps its height (pure shear), PQRS (45°) changes side lengths
    A0, A1 = frames[0]["ABCD"], frames[1]["ABCD"]
    assert np.ptp(A1[1]) == pytest.approx(np.ptp(A0[1]))
    leibniz_strips_figure()
    plt.close("all")
    res = rtt_blob_figure()
    plt.close("all")
    assert out is not None and res is not None
