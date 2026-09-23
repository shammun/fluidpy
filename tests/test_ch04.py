"""Verification suite for Chapter 4 — Conservation Laws (Kundu, Cohen & Dowling 5e, §§4.1–4.11).

Evidence levels (``verify-implementation`` skill): V1 analytic · V2 symbolic (sympy / dimensions) · V3 convergence ·
V4 conservation / invariant · V5 published benchmark (``reference/ch04/``, see SOURCES.md) · V6 book value (private,
``tests/book_values_ch04.json``, skipped when absent) · V7 limits / symmetry / invariance. Every test name is
``test_<concept>_<level>_<what>``; the comment on the ``def`` line repeats the level.

A items (CORE, ≥ 2 independent levels): C01 mass for a moving CV (4.5) · C02 continuity (4.7) · C03 stream function ·
C04 momentum for a moving CV (4.17) · C05 Bernoulli along a streamline (4.19) · C06 Cauchy (4.24) · C07 Newtonian stress
(4.31) · C08 incompressible Navier–Stokes (4.39b) · C09 noninertial frame (4.45) · C10 internal energy (4.57) · C11
Bernoulli function (4.71) · C12 unsteady Bernoulli (4.75) · C13 Boussinesq (4.86) · C14 kinematic BC (4.91) · C15
dimensionless NS (4.101). Derivations D01–D30 re-derived with sympy (``test_*_V2_derivation``); the ★★★ D09 and D15
re-run the construction of the design's check cells.

Pinned conventions with discrimination tests: outward n and the signed (u − b)·n · forces on the fluid (−F_D) · the
**first** index of τ is contracted (∂τ_ij/∂x_i) · 2-D ψ: ρu = ∂ψ/∂y, ρv = −∂ψ/∂x · Coriolis acceleration term +2Ω × u′,
force −2Ω × u′ (factor 2) · centrifugal −Ω × (Ω × x) = +Ω²R e_R · (4.74) gauge φ − ∫B dt (the printed + doubles B) ·
Stokes' assumption μ_v = 0 (not λ = 0) · plane (2 × 2) G is a plane flow of a 3-D fluid (the ⅓ deviator).

Run: ``.venv/Scripts/python.exe -m pytest tests/test_ch04.py -q -p no:cacheprovider``.
"""
from __future__ import annotations

import inspect
import json
import math
import re
from pathlib import Path

import numpy as np
import pytest
import sympy as sp
from scipy.integrate import quad, solve_ivp
from scipy.special import erfc

from fluidpy import ch01_introduction as ch01
from fluidpy import ch03_kinematics as ch03
from fluidpy import ch04_conservation_laws as ch04
from fluidpy.core import kinematics as K
from fluidpy.core import stratification as STRAT
from fluidpy.core import statics as STAT
from fluidpy.core import tensors as TN
from fluidpy.core.units import Q_, dimensional_check
from tools.convergence import observed_order, pairwise_orders

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference" / "ch04"
REF1 = ROOT / "reference" / "ch01"
BOOK = Path(__file__).resolve().parent / "book_values_ch04.json"
needs_ref = pytest.mark.skipif(not (REF / "benchmarks.json").exists(), reason="run reference/ch04/make_refs.py")
book_only = pytest.mark.skipif(not BOOK.exists(), reason="book values are private; see CLAUDE.md rule 9")

RNG = np.random.default_rng(4)
ORDER_TOL = 0.15  # design order ± this (verify-implementation default)
G = 9.81          # the book's (and the chapter's default) g


def book():
    return json.loads(BOOK.read_text(encoding="utf-8"))


def ref_json(name="benchmarks.json", folder=REF):
    return json.loads((folder / name).read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------------------------------------------------
# helpers: sympy expressions → field callables in the core.kinematics layout (components on axis 0, points last)
# ---------------------------------------------------------------------------------------------------------------------
X1, X2, X3, TT = sp.symbols("x1 x2 x3 t", real=True)


def vec_field(exprs, coords):
    """Vector field callable u(x, t) from sympy expressions (x of shape (d,) or (d, N))."""
    f = sp.lambdify((*coords, TT), list(exprs), "numpy")

    def u(x, t=0.0):
        x_ = np.asarray(x, dtype=float)
        shape = x_.shape[1:]
        return np.stack([np.broadcast_to(np.asarray(v, dtype=float), shape) for v in f(*x_, t)])
    return u


def sca_field(expr, coords):
    """Scalar field callable F(x, t) from a sympy expression."""
    f = sp.lambdify((*coords, TT), expr, "numpy")

    def F(x, t=0.0):
        x_ = np.asarray(x, dtype=float)
        return np.broadcast_to(np.asarray(f(*x_, t), dtype=float), x_.shape[1:]).copy()
    return F


def ten_field(rows, coords):
    """Tensor field callable τ(x, t) → (d, d, N) from a nested list of sympy expressions."""
    f = sp.lambdify((*coords, TT), [list(r) for r in rows], "numpy")

    def T(x, t=0.0):
        x_ = np.asarray(x, dtype=float)
        shape = x_.shape[1:]
        return np.array([[np.broadcast_to(np.asarray(v, dtype=float), shape) for v in row] for row in f(*x_, t)])
    return T


def rand_pts(n=40, d=3, lo=-1.0, hi=1.0, rng=RNG):
    return rng.uniform(lo, hi, (d, n))


def rel(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.max(np.abs(a - b)) / max(float(np.max(np.abs(b))), 1e-300))


# =====================================================================================================================
# §4.1 closure ledger (N02)
# =====================================================================================================================
def test_closure_count_V1_equation_unknown_ledger():  # V1
    c = ch04.closure_count("cauchy")
    assert (c["equations"], c["unknowns"], c["closed"]) == (6, 13, False)
    n = ch04.closure_count("navier_stokes")
    assert (n["equations"], n["unknowns"], n["closed"]) == (4, 5, False)
    assert ch04.closure_count("barotropic")["closed"] and ch04.closure_count("barotropic")["equations"] == 5
    f = ch04.closure_count("full")
    assert (f["equations"], f["unknowns"], f["closed"]) == (7, 7, True)
    assert set(f["unknown_names"]) >= {"ρ", "e", "p", "T"}
    e = ch04.closure_count("empty")
    assert (e["equations"], e["unknowns"], e["closed"]) == (0, 0, False)
    with pytest.raises(ValueError):
        ch04.closure_count("nope")


# =====================================================================================================================
# C01 — mass conservation for an arbitrarily moving control volume (4.1)–(4.5); D01
# =====================================================================================================================
def test_mass_cv_V1_interval_budget_worked_numbers():  # V1
    b = ch04.interval_mass_budget(1.0, 2.0, 0.0)  # fixed interval in the expanding flow (a = 1, ρ₀ = 1)
    assert b["storage"] == pytest.approx(-1.0, abs=1e-14)
    assert (b["flux_left"], b["flux_right"], b["net_outflux"]) == pytest.approx((1.0, 2.0, 1.0), abs=1e-14)
    assert abs(b["residual"]) < 1e-14 and b["local"] == pytest.approx(-1.0, abs=1e-14)
    m = ch04.interval_mass_budget(1.0, 2.0, 0.0, dx0dt=1.0, dx1dt=2.0)  # ends move with the fluid: material interval
    assert max(abs(m[k]) for k in ("storage", "flux_left", "flux_right", "net_outflux", "residual")) < 1e-14
    # an arbitrary moving interval at another time and rate: storage = Leibniz, residual 0
    for (x0, x1, t, v0, v1, a, r0) in [(0.3, 1.7, 0.4, -0.2, 0.9, 0.7, 2.0), (-1.0, 0.5, 2.0, 0.3, 0.1, 1.5, 1.2)]:
        b = ch04.interval_mass_budget(x0, x1, t, v0, v1, a=a, rho0=r0)
        s = 1 + a * t
        rho = r0 / s
        assert b["flux_left"] == pytest.approx(rho * (a * x0 / s - v0), rel=1e-13)
        assert b["storage"] == pytest.approx(-a * r0 / s ** 2 * (x1 - x0) + rho * (v1 - v0), rel=1e-12)
        assert abs(b["residual"]) < 1e-13


def test_mass_cv_V1_interval_budget_callable_route_agrees():  # V1 (quad + central difference vs the exact preset)
    a, r0 = 0.8, 1.3
    rho_fn = lambda x, t: r0 / (1 + a * t) + 0.0 * np.asarray(x, float)  # noqa: E731
    u_fn = lambda x, t: a * np.asarray(x, float) / (1 + a * t)  # noqa: E731
    for args in [(0.5, 1.5, 0.2, 0.0, 0.0), (0.5, 1.5, 0.2, -0.4, 0.7)]:
        ex = ch04.interval_mass_budget(*args, a=a, rho0=r0)
        nu = ch04.interval_mass_budget(*args, flow=(rho_fn, u_fn))
        for k in ("storage", "flux_left", "flux_right", "local"):
            assert nu[k] == pytest.approx(ex[k], rel=1e-8, abs=1e-10), k
        assert abs(nu["residual"]) < 1e-8
    # a non-uniform compressible pair: ρ = 2 + x² sin t, u = −x³ cos t/(3ρ) (∂ρ/∂t + ∂(ρu)/∂x = 0 exactly)
    rho2 = lambda x, t: 2 + np.asarray(x, float) ** 2 * np.sin(t)  # noqa: E731
    u2 = lambda x, t: -np.asarray(x, float) ** 3 * np.cos(t) / (3 * rho2(x, t))  # noqa: E731
    for args in [(0.2, 1.1, 0.7, 0.0, 0.0), (0.2, 1.1, 0.7, 0.3, -0.5)]:
        assert abs(ch04.interval_mass_budget(*args, flow=(rho2, u2))["residual"]) < 1e-7
    # a pair that violates continuity leaves a residual = ∫(∂ρ/∂t + ∂(ρu)/∂x) dx
    bad = ch04.interval_mass_budget(0.2, 1.1, 0.7, flow=(rho2, lambda x, t: 0.0 * np.asarray(x, float) + 1.0))
    expected = quad(lambda x: x ** 2 * np.cos(0.7) + 2 * x * np.sin(0.7), 0.2, 1.1)[0]
    assert bad["residual"] == pytest.approx(expected, rel=1e-6)


def test_mass_cv_V4_fixed_box_storage_balances_outflux():  # V4
    rho, u = ch04.expanding_flow_fields(a=0.7, rho0=1.4, dim=3)
    box = ch04.MovingBox(lengths=(1.0, 0.6, 0.8), origin=(0.2, -0.3, 0.5))
    for t in (0.0, 0.5, 2.0):
        b = ch04.mass_budget(rho, u, box, t, dt=1e-4)
        s = 1 + 0.7 * t
        exact = -3 * 0.7 * 1.4 / s ** 4 * (1.0 * 0.6 * 0.8)  # dM/dt = V dρ/dt for the fixed box
        assert b.storage == pytest.approx(exact, rel=1e-7)
        assert b.outflux == pytest.approx(-exact, rel=1e-12)
        assert abs(b.residual) < 1e-7 * abs(exact)
        assert b.mass == pytest.approx(1.4 / s ** 3 * 0.48, rel=1e-13) and math.isnan(b.local)
    lb = ch04.mass_budget(rho, u, box, 0.5, drho_dt=lambda X, T: -3 * 0.7 * 1.4 / (1 + 0.7 * T) ** 4 + 0 * X[0])
    assert lb.local == pytest.approx(lb.storage, rel=1e-7)  # fixed CV: d/dt passes inside (4.2)/(4.6)


def test_mass_cv_V4_material_volume_keeps_its_mass():  # V4 (4.1) with b = u
    rho, u = ch04.expanding_flow_fields(a=0.5, rho0=2.0, dim=3)
    for cv in (ch04.MovingBox(lengths=(1.0, 1.0, 0.5), origin=(0.1, 0.2, -0.4)), ch04.GrowingSphere(0.7, 0.0, (0.3, 0, 0))):
        b = ch04.mass_budget(rho, u, cv, 1.0, material=True)
        assert abs(b.outflux) < 1e-15
        assert abs(b.storage) < 1e-8 * b.mass  # RTT of the coincident material volume vanishes
    # 1-D material interval: mass constant for t ∈ [0, 5] (quad) and ends on the path lines
    ms = [ch04.material_mass(1.0, 2.0, t, a=0.9, rho0=1.7) for t in np.linspace(0, 5, 11)]
    assert np.ptp(ms) < 1e-12 and ms[0] == pytest.approx(1.7, rel=1e-13)
    _, u1 = ch04.expanding_flow_fields(a=0.9, rho0=1.7, dim=1)
    path = K.pathline(u1, np.array([2.0]), 0.0, np.linspace(0, 5, 6))
    ends = [ch04.material_interval(1.0, 2.0, t, a=0.9)[1] for t in np.linspace(0, 5, 6)]
    assert rel(path[0], ends) < 1e-9
    assert ch04.material_mass(1.0, 2.0, 1.0) == pytest.approx(1.0, rel=1e-13)  # the design's contract number
    assert ch04.material_mass(0.0, 1.0, 3.0, a=0.4, rho0=2.0, dim=3) == pytest.approx(2.0, rel=1e-13)


def test_mass_cv_V3_moving_growing_sphere_residual_order_2():  # V3
    rho, u = ch04.expanding_flow_fields(a=1.2, rho0=1.0, dim=3)
    cv = ch04.GrowingSphere(R0=0.5, Rdot=0.3, center=(0.2, -0.1, 0.4), U=(0.7, 0.2, -0.5))
    dts = [0.2, 0.1, 0.05, 0.025]
    errs = [abs(ch04.mass_budget(rho, u, cv, 0.3, dt=d).residual) for d in dts]
    p = observed_order(dts, errs)
    assert abs(p - 2.0) < ORDER_TOL, (p, pairwise_orders(dts, errs))
    assert errs[-1] < 1e-3 * abs(ch04.mass_budget(rho, u, cv, 0.3).storage)


def test_mass_cv_V2_derivation():  # V2 — D01 in 1-D: (4.1) → (4.2) → (4.3) + (4.4) → (4.5)
    x, t = sp.symbols("x t", real=True)
    c = sp.symbols("c0:6", real=True)
    rho = 2 + c[0] * x + c[1] * x ** 2 * t + c[2] * t ** 2  # generic smooth density
    # velocity from continuity: ρu = −∫∂ρ/∂t dx + f(t)
    rho_u = -sp.integrate(sp.diff(rho, t), x) + c[3] * t
    u = rho_u / rho
    assert sp.simplify(sp.diff(rho, t) + sp.diff(rho * u, x)) == 0  # (4.7) holds
    x0, x1 = 0.3 + c[4] * t, 1.2 + c[5] * t ** 2  # arbitrary moving ends: b = dx0/dt, dx1/dt
    storage = sp.diff(sp.integrate(rho, (x, x0, x1)), t)  # d/dt ∫_{V*} ρ dV
    outflux = (rho * (u - sp.diff(x1, t))).subs(x, x1) - (rho * (u - sp.diff(x0, t))).subs(x, x0)
    assert sp.simplify(storage + outflux) == 0  # (4.5)
    # step 2 (4.2): the material interval (b = u) — the coincident volume's d/dt∫ρ equals ∫∂ρ/∂t + [ρu]
    leibniz = sp.integrate(sp.diff(rho, t), (x, x0, x1)) + (rho * sp.diff(x1, t)).subs(x, x1) - (
        rho * sp.diff(x0, t)).subs(x, x0)
    assert sp.simplify(storage - leibniz) == 0  # RTT (3.35) in 1-D = (4.3)
    # (4.4): ∫∂ρ/∂t over V and V* coincide at the instant → subtracting the two RTT statements leaves (4.5)
    material = sp.integrate(sp.diff(rho, t), (x, x0, x1)) + rho_u.subs(x, x1) - rho_u.subs(x, x0)
    assert sp.simplify(material) == 0  # (4.2): the material volume's budget
    # discrimination: dropping b (using u·n instead of (u − b)·n) breaks (4.5)
    wrong = (rho * u).subs(x, x1) - (rho * u).subs(x, x0)
    assert sp.simplify(storage + wrong) != 0


# =====================================================================================================================
# C02 — continuity (4.7)–(4.10); D02, D03
# =====================================================================================================================
def test_continuity_V2_symbolic_expanding_flow_and_discrimination():  # V2
    t, a, r0 = sp.symbols("t a rho0", positive=True)
    xs = sp.symbols("x y z", real=True)
    for dim in (1, 2, 3):
        rho = r0 / (1 + a * t) ** dim
        u = [a * xi / (1 + a * t) for xi in xs[:dim]]
        assert sp.simplify(ch04.continuity_residual_sym(rho, u, xs[:dim], t)) == 0
        wrong = r0 / (1 + a * t) ** (dim + 1)  # a wrong exponent is not mass-conserving
        assert sp.simplify(ch04.continuity_residual_sym(wrong, u, xs[:dim], t)) != 0


def test_continuity_V1_numeric_residual_on_exact_fields():  # V1
    # 3-D expanding flow at 40 random points and three times
    rho, u = ch04.expanding_flow_fields(a=0.9, rho0=1.1, dim=3)
    X = rand_pts(40)
    for t in (0.0, 0.7, 3.0):
        terms = ch04.continuity_terms(rho, u, X, t, h=1e-3, ht=1e-5)
        s = 1 + 0.9 * t
        assert rel(terms.local, np.full(40, -3 * 0.9 * 1.1 / s ** 4)) < 1e-8
        assert rel(terms.flux_divergence, -terms.local) < 1e-8
        assert np.max(np.abs(ch04.continuity_residual(rho, u, X, t, h=1e-3, ht=1e-5))) < 1e-8
    # a non-uniform 1-D field: ρ = 2 + x² sin t, u = −x³ cos t/(3ρ)
    rho1 = sca_field(2 + X1 ** 2 * sp.sin(TT), (X1,))
    u1 = vec_field([-X1 ** 3 * sp.cos(TT) / (3 * (2 + X1 ** 2 * sp.sin(TT)))], (X1,))
    Xs = np.linspace(-1, 1, 21)[None, :]
    r = ch04.continuity_residual(rho1, u1, Xs, 0.4, h=1e-3, ht=1e-4)
    assert np.max(np.abs(r)) < 1e-6
    # discrimination: a wrong density (ρ + x) gives an O(1) residual
    wrong = ch04.continuity_residual(lambda X, T: rho1(X, T) + X[0], u1, Xs, 0.4, h=1e-3, ht=1e-4)
    assert np.max(np.abs(wrong)) > 1e-2


def test_continuity_V3_flux_divergence_stencil_order_2():  # V3
    rho1 = sca_field(2 + X1 ** 2 * sp.sin(TT) + sp.exp(X1) * sp.cos(TT), (X1,))
    uexpr = sp.sin(X1) * (1 + TT)
    u1 = vec_field([uexpr], (X1,))
    exact = sp.lambdify((X1, TT), sp.diff((2 + X1 ** 2 * sp.sin(TT) + sp.exp(X1) * sp.cos(TT)) * uexpr, X1))
    Xs = np.linspace(-1, 1, 9)[None, :]
    hs = [0.1, 0.05, 0.025, 0.0125]
    errs = [np.max(np.abs(ch04.continuity_terms(rho1, u1, Xs, 0.3, h=h).flux_divergence - exact(Xs[0], 0.3)))
            for h in hs]
    p = observed_order(hs, errs)
    assert abs(p - 2.0) < ORDER_TOL, (p, pairwise_orders(hs, errs))


def test_continuity_V1_material_form_and_incompressible_stratified():  # V1 (4.8)–(4.10), N10
    rho1 = sca_field(2 + X1 ** 2 * sp.sin(TT), (X1,))
    u1 = vec_field([sp.cos(X1) * TT], (X1,))  # not mass-conserving on purpose: (4.8) must still equal (4.7)/ρ
    Xs = np.linspace(-1, 1, 11)[None, :]
    mt = ch04.continuity_material_terms(rho1, u1, Xs, 0.5, h=1e-4, ht=1e-5)
    res = ch04.continuity_residual(rho1, u1, Xs, 0.5, h=1e-4, ht=1e-5)
    assert rel(np.asarray(mt.rel_rate) + mt.div_u, res / rho1(Xs, 0.5)) < 1e-6
    # stratified shear flow: incompressible although ρ varies (constant-density ⊂ incompressible)
    r3, u3 = ch04.stratified_shear_fields(U0=2.0, shear=0.3, rho0=1025.0, drho_dz=-0.4)
    X = rand_pts(30)
    assert np.max(np.abs(ch04.density_material_rate(r3, u3, X, 0.0))) < 1e-9
    ok, val = ch04.divergence_free_check(u3, X)
    assert ok and val < 1e-10
    grad_z = (r3(X + np.array([[0], [0], [1e-3]]), 0) - r3(X - np.array([[0], [0], [1e-3]]), 0)) / 2e-3
    assert np.allclose(grad_z, -0.4)  # ∂ρ/∂z ≠ 0
    ux, rz = ch04.stratified_shear_flow(np.array([0.0, 10.0]), 2.0, 0.3, 1025.0, -0.4)
    assert np.allclose(ux, [2.0, 5.0]) and np.allclose(rz, [1025.0, 1021.0])


def test_continuity_V7_galilean_shift_keeps_density_advection_zero():  # V7
    for U0 in (-5.0, 0.0, 12.0):
        r3, u3 = ch04.stratified_shear_fields(U0=U0, shear=0.3, drho_dz=-0.4)
        assert np.max(np.abs(ch04.density_material_rate(r3, u3, rand_pts(20), 1.0))) < 1e-9
    # an expanding flow is not divergence-free; the cylinder potential flow is
    _, u = ch04.expanding_flow_fields(a=0.5, dim=3)
    assert not ch04.divergence_free_check(u, rand_pts(10))[0]
    uc, _ = ch04.exact_field("cylinder", U=1.0, a=1.0)
    Xc = np.stack([np.array([1.5, -2.0, 0.3, 3.0]), np.array([0.5, 1.0, -1.8, 0.1])])
    assert ch04.divergence_free_check(uc, Xc, h=1e-4, tol=1e-7)[0]


@needs_ref
def test_mach_regime_V5_ussa_sound_speed_and_threshold():  # V5
    rows = [ln.split(",") for ln in (REF1 / "ussa1976_table1.csv").read_text().splitlines()[1:]]
    c_sl = float(rows[0][5])  # 340.29 m/s at 288.15 K (USSA-1976)
    c_ours = 100.0 / ch04.mach_number(100.0, T=288.15)
    assert abs(c_ours / c_sl - 1) < 1e-4
    assert ch04.mach_number(100.0, T=288.15) == pytest.approx(0.2939, abs=1e-4)
    assert ch04.is_incompressible_regime(100.0) and not ch04.is_incompressible_regime(105.0)
    assert ch04.incompressible_speed_limit() == pytest.approx(0.3 * c_ours, rel=1e-12)
    assert ch04.compressibility_parameter(0.3 * c_ours, c_ours) == pytest.approx(0.09, rel=1e-12)


def test_localisation_V1_ball_integral_mean_tends_to_point_value():  # V1 (N07, (4.6))
    f = lambda X: 1.0 + X[0] ** 2 + np.sin(X[1]) * X[2]  # noqa: E731
    x0 = np.array([0.3, -0.2, 0.5])
    assert ch04.ball_integral(lambda X: 1.0, x0, 0.4) == pytest.approx(4 / 3 * np.pi * 0.4 ** 3, rel=1e-10)
    rs = [0.4, 0.2, 0.1, 0.05]
    errs = [abs(ch04.ball_integral(f, x0, r) / (4 / 3 * np.pi * r ** 3) - f(x0)) for r in rs]
    assert abs(observed_order(rs, errs) - 2.0) < ORDER_TOL  # mean − point value = O(r²)


def test_continuity_V2_derivation():  # V2 — D02 (Gauss + localisation) and D03 (product rule → (4.8) → (4.10))
    x, y, z, t = sp.symbols("x y z t", real=True)
    c = sp.symbols("c0:9", real=True)
    # D02 step: ∮ρu·n dA over the unit box = ∫∇·(ρu) dV (Gauss) for a generic polynomial mass flux
    F = [c[0] * x ** 2 * y + c[1] * z, c[2] * y ** 2 * z + c[3] * x * t, c[4] * z ** 3 + c[5] * x * y * z]
    flux = sum(sp.integrate(F[i].subs(v, 1) - F[i].subs(v, 0), *[(w, 0, 1) for w in (x, y, z) if w != v])
               for i, v in enumerate((x, y, z)))
    vol = sp.integrate(sum(sp.diff(F[i], v) for i, v in enumerate((x, y, z))), (x, 0, 1), (y, 0, 1), (z, 0, 1))
    assert sp.simplify(flux - vol) == 0
    # D03: ∇·(ρu) = u·∇ρ + ρ∇·u ⇒ (4.7) ⇔ (1/ρ)Dρ/Dt + ∇·u = 0 (4.8); Dρ/Dt = 0 (4.9) ⇒ ∇·u = 0 (4.10)
    rho = sp.Function("rho", positive=True)(x, y, z, t)
    u = [sp.Function(f"u{i}")(x, y, z, t) for i in range(3)]
    C = ch04.continuity_residual_sym(rho, u, (x, y, z), t)
    Drho = sp.diff(rho, t) + sum(ui * sp.diff(rho, v) for ui, v in zip(u, (x, y, z)))
    divu = sum(sp.diff(ui, v) for ui, v in zip(u, (x, y, z)))
    assert sp.simplify(sp.expand(C / rho - (Drho / rho + divu))) == 0
    assert sp.simplify((C / rho).subs(Drho, 0) - divu) == 0 or sp.simplify(sp.expand(C - rho * divu - Drho)) == 0


# =====================================================================================================================
# C03 — stream functions (4.11)–(4.12), plane and axisymmetric; D04
# =====================================================================================================================
def test_streamfunction_V1_presets_velocity_matches_closed_form():  # V1
    rng = np.random.default_rng(11)
    for name, kw in [("uniform", dict(U=2.0, alpha=0.4)), ("stagnation", dict(k=1.5)), ("source_stream", dict(U=1, m=2)),
                     ("cylinder", dict(U=1.3, a=0.8)), ("vortex", dict(Gamma=3.0)), ("shear", dict(gamma_dot=0.7))]:
        r = rng.uniform(1.2, 3.0, 200)
        th = rng.uniform(0, 2 * np.pi, 200)
        x, y = r * np.cos(th), r * np.sin(th)
        u, v = ch04.velocity_from_streamfunction_2d(name, x, y, h=1e-5, **kw)
        ue, ve = ch04.velocity_preset(name, x, y, **kw)
        assert rel(u, ue) < 1e-7 and rel(v, ve) < 1e-7, name
    # the cylinder preset equals ch03's field and ψ = ch03's stream function
    x, y = np.array([1.5, -2.0, 0.4]), np.array([0.3, 1.1, -1.7])
    uc, vc = ch03.cylinder_flow(x, y, 1.3, 0.8)
    assert rel(ch04.velocity_preset("cylinder", x, y, U=1.3, a=0.8)[0], uc) < 1e-13
    assert rel(ch04.streamfunction_preset("cylinder", x, y, U=1.3, a=0.8), ch03.cylinder_streamfunction(x, y, 1.3, 0.8)) < 1e-13
    # sign convention: ψ = Uy is a stream toward +x (ρu = +∂ψ/∂y); ψ = kxy is (kx, −ky); ρ divides
    assert ch04.velocity_from_streamfunction_2d(lambda X, Y: 2.0 * Y, 0.3, 0.1) == pytest.approx((2.0, 0.0), abs=1e-9)
    assert ch04.velocity_from_streamfunction_2d(lambda X, Y: X * Y, 2.0, 3.0) == pytest.approx((2.0, -3.0), abs=1e-8)
    assert ch04.velocity_from_streamfunction_2d(lambda X, Y: 4.0 * Y, 0.3, 0.1, rho=2.0)[0] == pytest.approx(2.0)


def test_streamfunction_V2_symbolic_velocity_and_divergence_free():  # V2
    x, y, R, z = sp.symbols("x y R z", positive=True)
    f = sp.Function("psi")(x, y)
    u, v = ch04.velocity_from_streamfunction_2d_sym(f, x, y)
    assert sp.simplify(sp.diff(u, x) + sp.diff(v, y)) == 0
    rho = sp.Function("rho")(x, y)
    ru, rv = ch04.velocity_from_streamfunction_2d_sym(f, x, y, rho=rho)
    assert sp.simplify(sp.diff(rho * ru, x) + sp.diff(rho * rv, y)) == 0  # (4.11)
    g = sp.Function("psi")(R, z)
    uR, uz = ch04.velocity_from_streamfunction_axisym_sym(g, R, z)
    assert sp.simplify(sp.diff(R * uR, R) / R + sp.diff(uz, z)) == 0  # axisymmetric ∇·u
    uR0, uz0 = ch04.velocity_from_streamfunction_axisym_sym(sp.Rational(1, 2) * 3 * R ** 2, R, z)
    assert (uR0, sp.simplify(uz0)) == (0, 3)


def test_streamfunction_V1_flux_between_points_is_psi_difference():  # V1 / V4
    rng = np.random.default_rng(12)
    for name, kw in [("cylinder", dict(U=1.0, a=1.0)), ("stagnation", dict(k=2.0)), ("vortex", dict(Gamma=2.0)),
                     ("uniform", dict(U=1.5, alpha=0.3))]:
        for _ in range(4):
            r1, r2 = rng.uniform(1.3, 3, 2)
            t1, t2 = rng.uniform(0.1, 1.4, 2)
            p1, p2 = (r1 * np.cos(t1), r1 * np.sin(t1)), (r2 * np.cos(t2), r2 * np.sin(t2))
            dpsi = ch04.streamfunction_preset(name, *p2, **kw) - ch04.streamfunction_preset(name, *p1, **kw)
            assert ch04.flux_between_streamlines(name, p1, p2, **kw) == pytest.approx(dpsi, abs=1e-8), name
    pts = np.array([[1.5, 2.0, 2.5, 3.0], [0.2, 1.4, 0.8, 2.2]])  # a bent gate
    dpsi = ch04.streamfunction_preset("cylinder", 3.0, 2.2) - ch04.streamfunction_preset("cylinder", 1.5, 0.2)
    assert ch04.flux_along_path("cylinder", pts) == pytest.approx(dpsi, abs=1e-8)
    assert ch04.flux_between_streamlines(lambda X, Y: 2 * X * Y, (1, 1), (2, 3), rho=1.5) == pytest.approx(1.5 * (12 - 2), rel=1e-8)


def test_streamfunction_V3_velocity_stencil_order_2():  # V3
    x, y = np.array([1.4, -1.8, 0.6]), np.array([0.7, 1.2, -1.9])
    ue, ve = ch04.velocity_preset("cylinder", x, y)
    hs = [0.08, 0.04, 0.02, 0.01]
    errs = [max(np.max(np.abs(np.subtract(*p))) for p in zip(ch04.velocity_from_streamfunction_2d("cylinder", x, y, h=h),
                                                               (ue, ve))) for h in hs]
    assert abs(observed_order(hs, errs) - 2.0) < ORDER_TOL


def test_streamfunction_V1_axisymmetric_uniform_and_sphere():  # V1
    uR, uz = ch04.velocity_from_streamfunction_axisym("axisym_uniform", np.array([0.5, 2.0]), np.array([0.0, 1.0]), U=3.0)
    assert np.allclose(uR, 0, atol=1e-9) and np.allclose(uz, 3.0, rtol=1e-9)
    with pytest.raises(ValueError):
        ch04.velocity_from_streamfunction_axisym("axisym_uniform", 0.0, 1.0)
    # potential flow past a sphere: ψ = ½UR²(1 − a³/r³), u_z on the equator r = a is 3U/2
    U, a = 2.0, 1.0
    psi = lambda R, z: 0.5 * U * R ** 2 * (1 - a ** 3 / (R ** 2 + z ** 2) ** 1.5)  # noqa: E731
    uR, uz = ch04.velocity_from_streamfunction_axisym(psi, 1.0, 0.0)
    assert (uR, uz) == pytest.approx((0.0, 1.5 * U), abs=1e-7)


def test_stream_functions_3d_V2_two_surfaces_and_flux_patch():  # V2 + V1 (N14–N16)
    x, y, z = sp.symbols("x y z", real=True)
    chi, psi = sp.Function("chi")(x, y, z), sp.Function("psi")(x, y, z)
    chk = ch04.stream_surface_check(chi, psi, (x, y, z))
    assert all(sp.simplify(v) == 0 for v in chk.values())
    Psi = [chi * sp.diff(psi, v) for v in (x, y, z)]  # Ψ = χ∇ψ (4.12)
    m1 = ch04.mass_flux_from_vector_potential(Psi, (x, y, z))
    m2 = ch04.mass_flux_from_stream_functions(chi, psi, (x, y, z))
    assert all(sp.simplify(a_ - b_) == 0 for a_, b_ in zip(m1, m2))
    ce, pe, cs = ch04.stream_function_pair("parabolic")
    assert [sp.simplify(e) for e in ch04.mass_flux_from_stream_functions(ce, pe, cs)] == [1, 0, 2 * cs[0]]
    num, closed = ch04.stream_tube_mass_flux("parabolic", a=0.2, b=1.1, c=-0.5, d=0.9)
    assert num == pytest.approx(closed, rel=1e-9) and closed == pytest.approx(0.9 * 1.4)
    ex = ch04.stream_surface_example()
    num2, _ = ch04.stream_tube_mass_flux(ex["chi"], ex["psi"], patch=ex["patch"], a=1.0, b=2.0, c=0.0, d=1.0)
    assert abs(num2) == pytest.approx(ex["expected"], rel=1e-8)
    X = rand_pts(5, lo=0.2, hi=1.0)
    nchk = ch04.stream_surface_check("parabolic", x=X)
    assert np.max(np.abs(nchk["u_dot_grad_chi"])) < 1e-8 and np.max(np.abs(nchk["u_dot_grad_psi"])) < 1e-8
    assert np.max(np.abs(nchk["div"])) < 1e-6


def test_streamfunction_V2_derivation():  # V2 — D04: χ = −z ⇒ ρu = ∂ψ/∂y, ρv = −∂ψ/∂x; ψ₂ − ψ₁ = flux
    x, y, z, s = sp.symbols("x y z s", real=True)
    psi = sp.Function("psi")(x, y)
    m = ch04.mass_flux_from_stream_functions(-z, psi, (x, y, z))
    assert sp.simplify(m[0] - sp.diff(psi, y)) == 0 and sp.simplify(m[1] + sp.diff(psi, x)) == 0 and m[2] == 0
    # flux across a parametric curve with the right-hand normal: ∫(u n_x + v n_y) ds = ∫ dψ = ψ(end) − ψ(start)
    P = 3 * x ** 2 * y - y ** 3 + x  # a concrete ψ
    X_, Y_ = s ** 2 + 1, sp.sin(s)  # a curved gate, s ∈ [0, 1]
    u, v = sp.diff(P, y), -sp.diff(P, x)
    integrand = (u * sp.diff(Y_, s) - v * sp.diff(X_, s)).subs({x: X_, y: Y_})  # (u, v)·(y′, −x′)
    val = sp.integrate(sp.expand(integrand), (s, 0, 1))
    assert sp.simplify(val - (P.subs({x: X_, y: Y_}).subs(s, 1) - P.subs({x: X_, y: Y_}).subs(s, 0))) == 0


# =====================================================================================================================
# C04 — momentum for an arbitrarily moving control volume (4.13)–(4.18); D05; Examples 4.1, 4.3, 4.4, 4.6
# =====================================================================================================================
def _solid_body_3d(Om=0.8, rho=1000.0, g=G):
    u = lambda X, T: np.stack([-Om * X[1], Om * X[0], 0.0 * X[0]])  # noqa: E731
    p = lambda X: 0.5 * rho * Om ** 2 * (X[0] ** 2 + X[1] ** 2) - rho * g * X[2]  # noqa: E731
    tau = lambda X, T: -p(X) * np.eye(3).reshape(3, 3, *([1] * (np.ndim(X) - 1)))  # noqa: E731
    return u, p, tau


def test_momentum_cv_V1_steady_rotation_fixed_and_moving_boxes():  # V1 (4.17)
    rho = 1000.0
    u, p, tau = _solid_body_3d(rho=rho)
    box = ch04.MovingBox(lengths=(0.6, 0.5, 0.4), origin=(0.3, -0.2, 0.1))
    b = ch04.momentum_budget(rho, u, box, 0.0, g=(0, 0, -G), tau=tau)
    scale = max(np.max(np.abs(b.outflux)), np.max(np.abs(b.body)))
    assert np.max(np.abs(b.storage)) < 1e-10 * scale
    assert np.max(np.abs(b.outflux)) > 1.0  # a nontrivial momentum flux
    assert np.max(np.abs(b.residual)) < 1e-9 * scale
    # the same steady field seen by a translating, growing box: storage ≠ 0, (4.17) still closes
    mv = ch04.MovingBox(lengths=(0.6, 0.5, 0.4), rates=(0.1, -0.05, 0.2), velocity=(0.4, -0.3, 0.2), origin=(0.3, -0.2, 0.1))
    bm = ch04.momentum_budget(rho, u, mv, 0.5, g=(0, 0, -G), tau=tau)
    assert np.max(np.abs(bm.storage)) > 1.0
    assert np.max(np.abs(bm.residual)) < 1e-7 * max(np.max(np.abs(bm.outflux)), np.max(np.abs(bm.body)))
    # traction callable route = tau route
    tr = lambda X, N, T: -p(X) * N  # noqa: E731
    bt = ch04.momentum_budget(rho, u, box, 0.0, g=(0, 0, -G), traction=tr)
    assert rel(bt.surface, b.surface) < 1e-13


def test_momentum_cv_V1_uniform_flow_and_archimedes():  # V1
    box = ch04.MovingBox(lengths=(1.0, 2.0, 0.5), origin=(0.5, 0.0, -1.0))
    b = ch04.momentum_budget(1.2, lambda X, T: np.broadcast_to(np.array([[3.0], [1.0], [-2.0]]), X.shape), box, 0.0,
                             g=(0, 0, 0))
    assert np.max(np.abs(b.outflux)) < 1e-12 and np.max(np.abs(b.residual)) < 1e-12
    rho = 998.0
    ps = lambda X: 1e5 - rho * G * X[2]  # noqa: E731
    rest = ch04.momentum_budget(rho, lambda X, T: 0.0 * X, box, 0.0, g=(0, 0, -G), traction=lambda X, N, T: -ps(X) * N)
    V = 1.0
    assert rest.surface == pytest.approx(-rest.body, abs=1e-8)
    assert rest.surface[2] == pytest.approx(rho * G * V, rel=1e-12)
    F = STAT.net_pressure_force_on_box(lambda x, y, z: 1e5 - rho * 9.81 * z, (0.5, 1.5, 0.0, 2.0, -1.0, -0.5))
    assert rel(rest.surface, F) < 1e-12  # ch01's Archimedes routine


def test_momentum_cv_V3_storage_difference_order_2():  # V3
    # unsteady field (Taylor–Green in x–y, uniform in z) through a moving, growing box: residual → 0 at order 2 in dt
    nu = 0.05
    uexpr = [sp.sin(X1) * sp.cos(X2) * sp.exp(-2 * nu * TT), -sp.cos(X1) * sp.sin(X2) * sp.exp(-2 * nu * TT), 0 * X1]
    pexpr = (sp.cos(2 * X1) + sp.cos(2 * X2)) / 4 * sp.exp(-4 * nu * TT)
    u = vec_field(uexpr, (X1, X2, X3))
    sig = [[nu * (sp.diff(uexpr[i], v) + sp.diff(uexpr[j], w)) for j, v in enumerate((X1, X2, X3))]
           for i, w in enumerate((X1, X2, X3))]
    tau = ten_field([[-pexpr * (1 if i == j else 0) + sig[i][j] for j in range(3)] for i in range(3)], (X1, X2, X3))
    cv = ch04.MovingBox(lengths=(0.7, 0.5, 0.3), rates=(0.3, 0.2, 0.1), velocity=(0.4, -0.2, 0.1), origin=(0.1, 0.3, 0))
    dts = [0.4, 0.2, 0.1, 0.05]
    errs = [np.max(np.abs(ch04.momentum_budget(1.0, u, cv, 0.5, g=(0, 0, 0), tau=tau, dt=d).residual)) for d in dts]
    assert abs(observed_order(dts, errs) - 2.0) < ORDER_TOL, pairwise_orders(dts, errs)


def test_momentum_cv_V4_material_volume_is_newton_second_law():  # V4 (4.13): b = u
    u, p, tau = _solid_body_3d()
    box = ch04.MovingBox(lengths=(0.6, 0.5, 0.4), origin=(0.3, -0.2, 0.1))
    b = ch04.momentum_budget(1000.0, u, box, 0.0, g=(0, 0, -G), tau=tau, material=True,
                             d_rho_u_dt=lambda X, T: 0.0 * X)
    assert np.max(np.abs(b.outflux)) < 1e-12
    # storage = ∮ρu(u·n) (momentum carried by the moving material surface) = body + surface
    assert np.max(np.abs(b.residual)) < 1e-9 * np.max(np.abs(b.body))


def test_body_force_V1_potential_gravity_and_closed_loop():  # V1 (4.18)
    X = rand_pts(10)
    assert np.allclose(ch04.body_force_from_potential("gravity", X), np.array([[0], [0], [-G]]), atol=1e-9)
    Phi = lambda X: X[0] ** 2 * X[1] + np.sin(X[2])  # noqa: E731
    gv = ch04.body_force_from_potential(Phi, X, h=1e-4)
    ex = -np.stack([2 * X[0] * X[1], X[0] ** 2, np.cos(X[2])])
    assert rel(gv, ex) < 1e-7
    assert ch04.gravity_potential(2.0) == pytest.approx(2 * G)
    from fluidpy.core.integral_theorems import circulation, rectangle_loop
    loop = rectangle_loop((0.2, 0.1, 0.3), (0.3, 0.5, 0.8), 0.7, 0.4, n=64)
    work = circulation(lambda x, y, z: ch04.body_force_from_potential(Phi, np.stack([x, y, z])), loop)
    ref_work = circulation(lambda x, y, z: np.stack([-2 * x * y, -x ** 2 + 0 * z, -np.cos(z)]), loop)
    nonc = circulation(lambda x, y, z: np.stack([-y, x, 0 * z]), loop)  # a rotational field does work
    assert abs(nonc) > 0.1
    # a conservative force does no work round a closed loop (to the loop quadrature's error, same for the exact field)
    assert abs(work - ref_work) < 1e-9 and abs(work) < 1e-8 * abs(nonc)


def test_wake_drag_V1_gaussian_closed_form_and_mass_balance():  # V1 (Ex. 4.1)
    for (Ui, De, w, rho) in [(10.0, 2.0, 0.1, 1.2), (5.0, 0.5, 0.3, 1000.0)]:
        F = ch04.wake_drag_per_span("gaussian", Ui, rho, 40 * w, deficit=De, width=w)
        closed = rho * (Ui * De * np.sqrt(np.pi) * w - De ** 2 * np.sqrt(np.pi / 2) * w)
        assert F == pytest.approx(closed, rel=1e-10)
        assert ch04.wake_drag_gaussian(Ui, De, w, rho) == pytest.approx(closed, rel=1e-14)
        assert ch04.wake_side_outflow("gaussian", Ui, 40 * w, deficit=De, width=w) == pytest.approx(De * np.sqrt(np.pi) * w,
                                                                                                    rel=1e-10)
    assert ch04.wake_drag_per_span("gaussian", 10.0, 1.2, 2.0) == pytest.approx(3.6523, abs=5e-5)
    assert ch04.wake_side_outflow("gaussian", 10.0, 2.0) == pytest.approx(0.35449, abs=5e-6)
    F, err = ch04.wake_drag_per_span("gaussian", 10.0, 1.2, 2.0, return_error=True)
    assert err < 1e-9 * F


def test_wake_drag_V7_limits_sign_and_box_height():  # V7
    Fs = [ch04.wake_drag_per_span("gaussian", 10.0, 1.2, 2.0, deficit=d) for d in (1e-3, 1e-2, 1e-1)]
    assert Fs[0] < Fs[1] < Fs[2] and Fs[0] == pytest.approx(1.2 * 10 * 1e-3 * np.sqrt(np.pi) * 0.1, rel=1e-3)
    Hs = [ch04.wake_drag_per_span("gaussian", 10.0, 1.2, H) for H in (1.0, 2.0, 8.0)]
    assert np.ptp(Hs) < 1e-10 * Hs[0]  # independent of H once H ≫ b
    assert ch04.wake_drag_per_span(lambda y: 12.0 - np.exp(-y ** 2 / 0.01), 10.0, 1.2, 2.0) < 0  # U > U∞: thrust
    assert ch04.gaussian_wake(0.0, 10.0, 2.0, 0.1) == pytest.approx(8.0)


def test_wake_drag_V2_derivation_closed_form():  # V2 (N29, a-D09)
    y, Ui, De, b, rho = sp.symbols("y U_inf Delta b rho", positive=True)
    U = Ui - De * sp.exp(-y ** 2 / b ** 2)
    F = rho * sp.integrate(sp.expand(U * (Ui - U)), (y, -sp.oo, sp.oo))
    assert sp.simplify(F - rho * (Ui * De * sp.sqrt(sp.pi) * b - De ** 2 * sp.sqrt(sp.pi / 2) * b)) == 0


def test_bore_V1_speed_limit_and_pressure_force():  # V1 / V7 (Ex. 4.3)
    assert ch04.bore_speed(1.0, 1.1) == pytest.approx(3.36609, abs=5e-6)  # (design Part C prints 3.3660: truncated)
    h = 2.0
    errs, dhs = [], [0.2, 0.1, 0.05, 0.025]
    for dh in dhs:
        errs.append(abs(ch04.bore_speed(h, h + dh) / np.sqrt(G * h) - 1))
    assert abs(observed_order(dhs, errs) - 1.0) < ORDER_TOL  # U → √(gh) with error O(Δh/h)
    for po in (0.0, 1e5, 3e5):
        f = ch04.bore_pressure_force(1.0, 1.3, rho=1000.0, p_o=po)
        assert f["net"] == pytest.approx(1000.0 * G * (1.0 - 1.3 ** 2) / 2, rel=1e-12)  # p_o cancels
    assert ch04.bore_outlet_velocity(1.0, 1.2) == pytest.approx(-0.2 * ch04.bore_speed(1.0, 1.2) / 1.2, rel=1e-14)


@needs_ref
def test_bore_V1_belanger_form_cross_check():  # V1 (form cross-check, not V5)
    ref = ref_json()["belanger_jump"]
    f = sp.lambdify(sp.Symbol("Fr1"), sp.sympify(ref["h2_over_h1"]))
    for (hi, ho) in [(1.0, 1.1), (0.5, 1.5), (2.0, 2.02)]:
        U = ch04.bore_speed(hi, ho)
        Fr1 = U / np.sqrt(G * hi)  # the still water enters the jump at U in the wave's frame
        assert f(Fr1) == pytest.approx(ho / hi, rel=1e-13)


def test_bore_V2_derivation_mass_and_momentum():  # V2 (a-D11): moving CV b = −U e_x
    U, Uo, hi, ho, rho, g = sp.symbols("U U_out h_in h_out rho g", positive=True)
    mass = sp.Eq(U * hi, (Uo + U) * ho)
    Uo_s = sp.solve(mass, Uo)[0]
    # momentum per width, CV riding with the wave: ρ(U_out + U)h_out·U_out − 0 = ρg(h_in² − h_out²)/2 … in the lab
    # x-momentum: outflow of x-momentum ρU_out(U_out + U)h_out = net hydrostatic force
    mom = sp.Eq(rho * Uo_s * (Uo_s + U) * ho, rho * g * (hi ** 2 - ho ** 2) / 2)
    sol = [s for s in sp.solve(mom, U) if s.is_positive is not False]
    target = sp.sqrt(g * ho * (hi + ho) / (2 * hi))
    assert any(sp.simplify(s ** 2 - target ** 2) == 0 for s in sol)
    assert sp.simplify(Uo_s - (hi - ho) * U / ho) == 0


@needs_ref
def test_rocket_V1_tsiolkovsky_and_closed_form_trajectory():  # V1 (form cross-check) / V3
    ref = ref_json()["tsiolkovsky"]
    assert "log(m0/m_f)" in ref["formula"]
    r = ch04.rocket_trajectory(M0=2.0, mdot=0.1, Ve=800.0, t_burn=12.0, g=0.0)
    assert r["b"][-1] == pytest.approx(800.0 * np.log(2.0 / 0.8), rel=1e-8)
    assert ch04.rocket_delta_v(2.0, 0.8, 800.0) == pytest.approx(800.0 * np.log(2.5), rel=1e-15)
    rg = ch04.rocket_trajectory(M0=2.0, mdot=0.1, Ve=800.0, t_burn=12.0)
    z, b = ch04.rocket_closed_form(rg["t"], 2.0, 0.1, 800.0)
    assert rel(rg["b"], b) < 1e-8 and rel(rg["z"], z) < 1e-8
    tols = [1e-4, 1e-6, 1e-8, 1e-10]
    errs = [abs(ch04.rocket_trajectory(2.0, 0.1, 800.0, 12.0, rtol=r_, atol=1e-12)["z"][-1] - z[-1]) for r_ in tols]
    assert errs[-1] < errs[0] and errs[-1] < 1e-7 * abs(z[-1])
    with pytest.raises(ValueError):
        ch04.rocket_trajectory(1.0, 0.1, 500.0, 11.0)


def test_rocket_V2_derivation():  # V2 (a-D12): mass + momentum of the accelerating CV ⇒ closed forms
    t, M0, md, Ve, g = sp.symbols("t M0 mdot V_e g", positive=True)
    M = M0 - md * t
    z, b = ch04.rocket_closed_form(0.5, 2.0, 0.1, 800.0)  # numeric smoke
    assert np.isfinite(z) and b > 0
    b_s = -Ve * sp.log(1 - md * t / M0) - g * t
    assert sp.simplify(M * sp.diff(b_s, t) - (-Ve * sp.diff(M, t) - M * g)) == 0  # M d²z/dt² = −V_e dM/dt − Mg


def test_cv_scenarios_V1_every_budget_closes():  # V1 (E1)
    for name in ch04.CV_SCENARIOS:
        s = ch04.cv_scenario(name)
        scale = max(1.0, max(abs(f["momentum_flux_x"]) for f in s["faces"]))
        assert abs(s["residual_mass"]) < 1e-8 * scale and abs(s["residual_momentum"]) < 1e-8 * scale, name
    for b in (0.0, 1.0, 5.0):  # a CV that does not ride with the bore still closes (storage picks up the drift)
        s = ch04.cv_scenario("bore", b=b)
        assert abs(s["residual_mass"]) < 1e-9 and abs(s["residual_momentum"]) < 1e-7
    s = ch04.cv_scenario("wake")
    assert s["result"] == pytest.approx(ch04.wake_drag_per_span("gaussian", 10.0, 1.2, 2.0), rel=1e-12)
    assert ch04.cv_scenario("jet", theta=np.pi / 6)["result"] == pytest.approx(ch04.jet_plate_force(1000, 10, 1e-3, np.pi / 6))
    assert ch04.jet_plate_force(1000.0, 10.0, 1e-3) == pytest.approx(100.0)
    with pytest.raises(ValueError):
        ch04.cv_scenario("nope")


def test_angular_momentum_V1_sprinkler_numeric_flux_and_free_spin():  # V1 (Ex. 4.6, (4.65))
    for alpha in (0.0, np.pi / 6, np.pi / 3):
        M = ch04.sprinkler_torque(0.2, 1000.0, 1e-4, 5.0, alpha)
        assert ch04.sprinkler_torque_numeric(0.2, 1000.0, 1e-4, 5.0, alpha) == pytest.approx(M, rel=1e-10)
    assert ch04.sprinkler_torque(0.2, 1000.0, 1e-4, 5.0, np.pi / 6) == pytest.approx(0.8660, abs=5e-5)
    assert abs(ch04.sprinkler_torque(0.2, 1000.0, 1e-4, 5.0, np.pi / 2)) < 1e-15
    om = ch04.sprinkler_free_spin_rate(0.2, 5.0, np.pi / 6)
    assert 5.0 * np.cos(np.pi / 6) - om * 0.2 == pytest.approx(0.0, abs=1e-14)  # zero absolute swirl ⇒ zero torque


def test_angular_momentum_V4_budget_closes_for_rigid_rotation():  # V4 (4.65)
    u, p, tau = _solid_body_3d()
    box = ch04.MovingBox(lengths=(0.6, 0.5, 0.4), origin=(0.3, -0.2, 0.1))
    b = ch04.angular_momentum_budget(1000.0, u, box, 0.0, g=(0, 0, -G), tau=tau, origin=(0.1, 0.2, -0.3))
    scale = max(np.max(np.abs(b.outflux)), np.max(np.abs(b.body_torque)))
    assert scale > 1.0 and np.max(np.abs(b.residual)) < 1e-9 * scale
    f = ch04.angular_momentum_flux(1000.0, u, box, 0.0, origin=(0.1, 0.2, -0.3))
    assert rel(f, b.outflux) < 1e-13


def test_momentum_cv_V2_derivation():  # V2 — D05 in 1-D: (4.13) → (4.14) → (4.15) → (4.16) → (4.17)
    x, t = sp.symbols("x t", real=True)
    c = sp.symbols("c0:6", real=True)
    rho = 3 + c[0] * x * t + c[1] * x ** 2
    rho_u = -sp.integrate(sp.diff(rho, t), x) + c[2] * t
    u = rho_u / rho
    # stress chosen so that Cauchy (4.24) holds exactly in 1-D with g: τ_xx from ∂τ/∂x = ρDu/Dt − ρg
    gsym = sp.Symbol("g")
    acc = sp.diff(u, t) + u * sp.diff(u, x)
    tau = sp.integrate(sp.expand(sp.simplify(rho * acc - rho * gsym)), x)
    x0, x1 = c[3] * t, 1 + c[4] * t ** 2
    storage = sp.diff(sp.integrate(rho * u, (x, x0, x1)), t)
    outflux = (rho * u * (u - sp.diff(x1, t))).subs(x, x1) - (rho * u * (u - sp.diff(x0, t))).subs(x, x0)
    body = sp.integrate(rho * gsym, (x, x0, x1))
    surface = tau.subs(x, x1) - tau.subs(x, x0)  # f = n τ: +τ at the right end, −τ at the left
    assert sp.simplify(storage + outflux - body - surface) == 0  # (4.17)
    # (4.15) is an identity (the book's trailing "= 0" is spurious): ∫∂(ρu)/∂t = d/dt∫ρu − [ρu b]
    ident = sp.integrate(sp.diff(rho * u, t), (x, x0, x1)) - (storage - ((rho * u * sp.diff(x1, t)).subs(x, x1)
                                                                        - (rho * u * sp.diff(x0, t)).subs(x, x0)))
    assert sp.simplify(ident) == 0
    assert sp.simplify(sp.integrate(sp.diff(rho * u, t), (x, x0, x1))) != 0  # … and it is not zero by itself


# =====================================================================================================================
# C05 — Bernoulli along a streamline (4.19); Ex. 4.2; pitot, orifice (N103–N105); D06
# =====================================================================================================================
def test_bernoulli_V1_cylinder_surface_pressure_coefficient():  # V1
    th = np.linspace(0, np.pi, 19)
    U, a, rho, pinf = 3.0, 0.5, 1.2, 1e5
    us, vs = ch04.velocity_preset("cylinder", a * np.cos(th), a * np.sin(th) + 0.0, U=U, a=a)
    speed = np.hypot(us, vs)
    B_inf = ch04.bernoulli_head(U, 0.0, pinf, rho, g=0.0)
    p = np.array([ch04.bernoulli_solve(dict(U=U, z=0, p=pinf), dict(U=s_, z=0), "p", rho=rho, g=0.0) for s_ in speed])
    Cp = ch04.pressure_coefficient(p, pinf, rho, U)
    assert rel(Cp, 1 - 4 * np.sin(th) ** 2) < 1e-12
    assert B_inf == pytest.approx(0.5 * 9 + pinf / rho)


def test_bernoulli_V4_constant_along_traced_streamlines():  # V4
    U, a, rho = 1.0, 1.0, 1000.0
    uf, pf = ch04.exact_field("cylinder", U=U, a=a, rho=rho, p_inf=2e5)
    for y0 in (0.3, 0.8, 1.6, -0.5):
        pts = K.streamline(uf, np.array([-4.0, y0]), 0.0, s_max=8.0, both=False, n=200)
        B = ch04.bernoulli_along_line(uf, pf, pts, rho=rho, z_axis=None)
        assert np.ptp(B) / abs(B[0]) < 1e-8
    # a vertical plane with gravity: hydrostatic column at rest keeps B = p/ρ + gz constant
    pts = np.stack([np.zeros(5), np.linspace(-2, 0, 5)])
    Bh = ch04.bernoulli_along_line(lambda X, T: 0.0 * X, lambda X, T: 1e5 - rho * G * X[1], pts, rho=rho)
    assert np.ptp(Bh) < 1e-9


def test_bernoulli_V1_solve_round_trips_and_raises():  # V1
    s1 = dict(U=2.0, z=1.0, p=1.2e5)
    for unknown, s2 in [("U", dict(z=0.4, p=1.1e5)), ("z", dict(U=4.0, p=1.0e5)), ("p", dict(U=0.5, z=3.0))]:
        val = ch04.bernoulli_solve(s1, s2, unknown, rho=998.0)
        full = dict(s2, **{unknown: val})
        assert ch04.bernoulli_head(full["U"], full["z"], full["p"], 998.0) == pytest.approx(
            ch04.bernoulli_head(2.0, 1.0, 1.2e5, 998.0), rel=1e-14)
    with pytest.raises(ValueError):
        ch04.bernoulli_solve(s1, dict(z=100.0, p=2e5), "U")
    with pytest.raises(ValueError):
        ch04.bernoulli_solve(s1, dict(z=1.0, p=2e5), "q")


def test_bernoulli_V2_dimensions_expose_example_typo():  # V2 (dimensional homogeneity; analysis §9 item 7)
    q = dict(U=Q_(2.0, "m/s"), z=Q_(1.0, "m"), p=Q_(1e5, "Pa"), rho=Q_(1000.0, "kg/m**3"), g=Q_(9.81, "m/s**2"))
    dimensional_check(lambda U, z, p, rho, g: 0.5 * U ** 2 + g * z + p / rho, "[length] ** 2 / [time] ** 2", **q)
    with pytest.raises(Exception):
        (0.5 * q["rho"] * q["U"] ** 2 + q["g"] * q["z"]).to("J/kg")  # the example's "(½)ρU² + gz" does not add up


def test_pitot_orifice_V1_speeds_heads_and_drain():  # V1 (N103–N105)
    assert ch04.pitot_speed(500.0, 0.0, 1.2) == pytest.approx(28.868, abs=5e-4)
    assert ch04.pitot_speed(ch04.stagnation_pressure(1e5, 12.0, 1.2), 1e5, 1.2) == pytest.approx(12.0, rel=1e-12)
    assert ch04.stagnation_pressure(1e5, 12.0, 1.2) - 1e5 == pytest.approx(ch04.dynamic_pressure(12.0, 1.2))
    assert ch04.pitot_speed_from_heads(0.1, 0.3) == pytest.approx(np.sqrt(2 * G * 0.2))
    assert ch04.pitot_speed_from_heads(0.1, 0.3, rho=1000.0, rho_atm=1.2) == pytest.approx(np.sqrt(2 * G * 0.2 * (1 - 1.2e-3)))
    with pytest.raises(ValueError):
        ch04.pitot_speed(1.0, 2.0, 1.0)
    assert ch04.torricelli_speed(1.0) == pytest.approx(4.4294, abs=5e-5)
    assert ch04.bernoulli_solve(dict(U=0, z=2.0, p=0), dict(z=0, p=0), "U") == pytest.approx(ch04.torricelli_speed(2.0))
    assert ch04.orifice_mass_flow(2.0, 1e-3, Cc=0.611) == pytest.approx(1000 * 0.611 * 1e-3 * np.sqrt(2 * G * 2.0))
    d = ch04.tank_drain(1.5, 0.5, 1e-3, Cc=0.611)
    kk = 0.611 * 1e-3 / 0.5 * np.sqrt(G / 2)
    assert rel(d["h"], (np.sqrt(1.5) - kk * d["t"]) ** 2) < 1e-8  # closed-form level
    assert d["t_empty"] == pytest.approx(d["t_empty_closed"], rel=1e-6)
    assert ch04.tank_drain_time(1.5, 1.0, 1e-3) > ch04.tank_drain_time(1.5, 0.5, 1e-3)  # V7 monotone in tank area


@needs_ref
def test_orifice_V5_contraction_coefficient_used():  # V5
    Cc = ref_json()["vena_contracta_sharp_orifice"]["Cc"]
    assert ch04.bernoulli_scenario("orifice")["numbers"]["jet_area"] == pytest.approx(Cc * 1e-3, rel=1e-12)
    ratio = ch04.orifice_mass_flow(1.0, 1e-3, Cc=Cc) / ch04.orifice_mass_flow(1.0, 1e-3, Cc=1.0)
    assert ratio == pytest.approx(0.611, rel=1e-12)


def test_bernoulli_V2_derivation():  # V2 — D06 (Ex. 4.2): the stream-tube element, (ds)² dropped → (4.19)
    d = ch04.stream_tube_element_balance_sym()
    assert sp.simplify(d["residual_vs_4_19"]) == 0
    # independent re-run: first-order Taylor of the element's mass and momentum balances
    rho, U, A, p, g, th, ds = sp.symbols("rho U A p g theta ds", positive=True)
    Us, As, ps = sp.symbols("U_s A_s p_s")
    mass = rho * (U + Us * ds) * (A + As * ds) - rho * U * A
    assert sp.simplify(sp.expand(mass).coeff(ds, 1) - rho * (U * As + A * Us)) == 0
    assert sp.simplify(d["mass"] - rho * (U * As + A * Us)) == 0
    mom_out = rho * (U + Us * ds) ** 2 * (A + As * ds) - rho * U ** 2 * A
    forces = p * A - (p + ps * ds) * (A + As * ds) + (p + ps * ds / 2) * As * ds - rho * g * sp.sin(th) * (A + As * ds / 2) * ds
    first = sp.expand(mom_out - U * mass - forces).coeff(ds, 1)
    assert sp.simplify(first / (rho * A) - (U * Us + g * sp.sin(th) + ps / rho)) == 0
    # without the side-pressure term the conical wall force is missing and the result is wrong
    wrong = sp.expand(mom_out - U * mass - (p * A - (p + ps * ds) * (A + As * ds))).coeff(ds, 1)
    assert sp.simplify(wrong / (rho * A) - (U * Us + ps / rho)) != 0
    # U U_s ds + g dz + dp/ρ = 0 with dz = sin θ ds is d(½U² + gz + p/ρ) = 0 (4.19)
    s = sp.Symbol("s")
    Uf, zf, pf = (sp.Function(n)(s) for n in ("U", "z", "p"))
    dB = sp.diff(Uf ** 2 / 2 + g * zf + pf / rho, s)
    assert sp.simplify(dB - (Uf * sp.diff(Uf, s) + g * sp.diff(zf, s) + sp.diff(pf, s) / rho)) == 0


# =====================================================================================================================
# C06 — Cauchy's equation (4.20)–(4.24); D07
# =====================================================================================================================
def test_cauchy_V1_rigid_rotation_both_forms():  # V1
    u, p, tau = _solid_body_3d(Om=1.3, rho=900.0)
    X = rand_pts(25)
    ct = ch04.cauchy_terms(900.0, u, tau, (0, 0, -G), X, 0.0, h=1e-3)
    scale = np.max(np.abs(ct.inertia))
    assert scale > 100 and np.max(np.abs(ct.residual)) < 1e-9 * scale
    assert rel(ct.inertia, 900.0 * np.stack([-1.3 ** 2 * X[0], -1.3 ** 2 * X[1], 0 * X[0]])) < 1e-9
    r2 = ch04.momentum_conservative_residual(900.0, u, tau, (0, 0, -G), X, 0.0, h=1e-3)
    assert np.max(np.abs(r2)) < 1e-9 * scale
    assert np.max(np.abs(ch04.cauchy_residual(900.0, u, tau, (0, 0, -G), X, 0.0, h=1e-3))) < 1e-9 * scale


def test_cauchy_V1_first_index_is_contracted():  # V1 — the non-symmetric discrimination (analysis §9 item 8)
    tau = lambda X, T: np.array([[0 * X[0], 0 * X[0], 0 * X[0]],  # noqa: E731
                                 [X[0], 0 * X[0], 0 * X[0]],  # τ₂₁ = x₁ (non-symmetric)
                                 [0 * X[0], 0 * X[0], 0 * X[0]]])
    X = rand_pts(6)
    first, second = ch04.divergence_first_index_demo(tau, X)
    assert np.allclose(first, 0.0, atol=1e-10)  # ∂τ_i1/∂x_i = ∂τ₂₁/∂x₂ = 0
    assert np.allclose(second[1], 1.0) and np.allclose(second[[0, 2]], 0.0, atol=1e-10)  # ∂τ₂j/∂x_j = 1
    # a fluid at rest with this stress: Cauchy is satisfied with the first index, violated with the second
    rest = lambda X, T: 0.0 * X  # noqa: E731
    assert np.max(np.abs(ch04.cauchy_residual(1.0, rest, tau, (0, 0, 0), X))) < 1e-10
    assert np.max(np.abs(ch04.stress_divergence(tau, X, index=1))) > 0.5
    # symmetric τ: both indices agree
    sym = lambda X, T: np.array([[X[1] ** 2, X[0] * X[2], X[1]], [X[0] * X[2], X[2], X[0]], [X[1], X[0], X[0] * X[1]]])  # noqa
    a_, b_ = ch04.divergence_first_index_demo(sym, X)
    assert rel(a_, b_) < 1e-9


def test_cauchy_V1_flux_form_minus_advective_is_u_times_continuity():  # V1 (4.23) numerically
    rho = sca_field(1 + 0.3 * sp.sin(X1 + TT) + X2 ** 2, (X1, X2, X3))
    u = vec_field([X2 * TT, sp.cos(X1), X1 * X3], (X1, X2, X3))  # does NOT satisfy continuity with this ρ
    tau = ten_field([[X1, 0, X2], [0, X3, 0], [X2, 0, X1 * X2]], (X1, X2, X3))
    X = rand_pts(12)
    d = ch04.momentum_conservative_residual(rho, u, tau, (0, 0, -G), X, 0.4, h=1e-4, ht=1e-5) - \
        ch04.cauchy_residual(rho, u, tau, (0, 0, -G), X, 0.4, h=1e-4, ht=1e-5)
    C = ch04.continuity_residual(rho, u, X, 0.4, h=1e-4, ht=1e-5)
    assert rel(d, u(X, 0.4) * C) < 1e-6


def test_cauchy_V3_residual_stencil_order_2():  # V3
    rho = 1.0
    uexpr = [sp.sin(X1) * sp.cos(X2), -sp.cos(X1) * sp.sin(X2), 0 * X1]
    u = vec_field(uexpr, (X1, X2, X3))
    p = -(sp.cos(2 * X1) + sp.cos(2 * X2)) / 4
    # inviscid steady Taylor–Green-like field: ρ(u·∇)u = −∇p exactly with p = −(cos 2x + cos 2y)/4 ... sign check below
    acc = [sum(uexpr[i] * sp.diff(uexpr[j], v) for i, v in enumerate((X1, X2, X3))) for j in range(3)]
    grad_p = [sp.diff(p, v) for v in (X1, X2, X3)]
    if any(sp.simplify(acc[j] + grad_p[j]) != 0 for j in range(3)):
        p = -p
    tau = ten_field([[-p if i == j else 0 for j in range(3)] for i in range(3)], (X1, X2, X3))
    X = rand_pts(8)
    hs = [0.1, 0.05, 0.025, 0.0125]
    errs = [np.max(np.abs(ch04.cauchy_residual(rho, u, tau, (0, 0, 0), X, 0.0, h=h, ht=1e-3))) for h in hs]
    assert abs(observed_order(hs, errs) - 2.0) < ORDER_TOL, pairwise_orders(hs, errs)


def test_cauchy_V2_conservative_to_advective_symbolic():  # V2 (4.23)
    x, y, z, t = sp.symbols("x y z t", real=True)
    rho = sp.Function("rho")(x, y, z, t)
    u = [sp.Function(f"u{i}")(x, y, z, t) for i in range(3)]
    diffs = ch04.conservative_to_advective_sym(rho, u, (x, y, z), t)
    C = ch04.continuity_residual_sym(rho, u, (x, y, z), t)
    assert all(sp.simplify(sp.expand(d_ - u[j] * C)) == 0 for j, d_ in enumerate(diffs))


def test_cauchy_V2_derivation():  # V2 — D07: Gauss per component pins the FIRST index; (4.20a)–(4.24)
    x, y, z = sp.symbols("x y z", real=True)
    coords = (x, y, z)
    tau = sp.Matrix(3, 3, lambda i, j: (i + 1) * x ** (j + 1) * y + (j + 2) * z ** 2 * x ** i)  # non-symmetric
    # ∮ n_i τ_ij dA over the unit box
    surf = []
    for j in range(3):
        tot = 0
        for i, v in enumerate(coords):
            others = [(w, 0, 1) for w in coords if w != v]
            tot += sp.integrate(tau[i, j].subs(v, 1) - tau[i, j].subs(v, 0), *others)
        surf.append(sp.simplify(tot))
    first = [sp.integrate(sum(sp.diff(tau[i, j], coords[i]) for i in range(3)), (x, 0, 1), (y, 0, 1), (z, 0, 1))
             for j in range(3)]
    second = [sp.integrate(sum(sp.diff(tau[j, i], coords[i]) for i in range(3)), (x, 0, 1), (y, 0, 1), (z, 0, 1))
              for j in range(3)]
    assert all(sp.simplify(a_ - b_) == 0 for a_, b_ in zip(surf, first))  # (4.20b): ∂τ_ij/∂x_i
    assert any(sp.simplify(a_ - b_) != 0 for a_, b_ in zip(surf, second))  # the prose's ∂τ_ij/∂x_j is different
    # (4.20a): ∮ρu_j u·n = ∫∂(ρu_iu_j)/∂x_i and the flux form (4.22) → (4.24) via (4.23) (checked above)
    t = sp.Symbol("t")
    rho = 1 + x * y + t
    u = [y * t, x ** 2, z]
    for j in range(3):
        lhs = 0
        for i, v in enumerate(coords):
            others = [(w, 0, 1) for w in coords if w != v]
            lhs += sp.integrate((rho * u[i] * u[j]).subs(v, 1) - (rho * u[i] * u[j]).subs(v, 0), *others)
        rhs = sp.integrate(sum(sp.diff(rho * u[i] * u[j], coords[i]) for i in range(3)), (x, 0, 1), (y, 0, 1), (z, 0, 1))
        assert sp.simplify(lhs - rhs) == 0


# =====================================================================================================================
# C07 — Newtonian constitutive law (4.25)–(4.37); D08, D09 ★★★, D10
# =====================================================================================================================
def test_newtonian_V1_parallel_shear_rest_and_rotation():  # V1
    for gam, mu in [(10.0, 1e-3), (250.0, 1.8e-5), (-3.0, 0.9)]:
        G_ = ch04.stress_lab_gradient("shear", gam)
        tau = ch04.newtonian_stress(G_, p=2.0, mu=mu)
        assert tau[0, 1] == pytest.approx(ch01.newton_shear_stress(mu, gam), rel=1e-14)  # link to (1.3)
        assert tau[0, 1] == tau[1, 0] and np.allclose(np.diag(tau), -2.0)
        assert ch04.shear_stress_parallel_flow(mu, gam) == pytest.approx(mu * gam)
    assert ch04.newtonian_stress([[0, 10, 0], [0, 0, 0], [0, 0, 0]], mu=1e-3)[0, 1] == pytest.approx(0.010)
    Rg = ch04.stress_lab_gradient("rotation", 4.0)
    assert np.allclose(ch04.viscous_stress(Rg, 1.0, 0.7), 0.0)  # rigid rotation: no viscous stress
    assert np.allclose(ch04.newtonian_stress(np.zeros((3, 3)), p=5.0), ch04.static_stress(5.0))
    assert np.allclose(ch04.static_stress(5.0), -5.0 * np.eye(3))
    assert np.allclose(ch04.total_stress(5.0, np.zeros((3, 3))), ch04.static_stress(5.0))
    n = np.array([0.6, 0.0, 0.8])
    assert np.allclose(TN.traction(ch04.static_stress(3.0), n), -3.0 * n)  # isotropic: traction = −p n


def test_newtonian_V7_isotropic_tensor_invariant_and_linear_law():  # V7 + V1
    K4 = ch04.isotropic_fourth_order(0.3, 1.1, 0.7)
    ok, res = TN.is_isotropic(K4, np.random.default_rng(1), 50)
    assert ok and res < 1e-12
    # (4.29) component by component: λδ_ijδ_mn + μδ_imδ_jn + γδ_inδ_jm
    assert (K4[0, 0, 1, 1], K4[0, 1, 0, 1], K4[0, 1, 1, 0], K4[2, 2, 2, 2]) == pytest.approx((0.3, 1.1, 0.7, 2.1))
    assert K4[0, 1, 2, 2] == 0.0 and K4[0, 0, 0, 1] == 0.0
    # the book's (4.30) argument: K symmetric in (i, j) ⇔ γ = μ
    assert not np.allclose(K4, np.swapaxes(K4, 0, 1))
    Ks = ch04.isotropic_fourth_order(0.3, 1.1, 1.1)
    assert np.allclose(Ks, np.swapaxes(Ks, 0, 1)) and np.allclose(Ks, np.swapaxes(Ks, 2, 3))
    # acting on a non-symmetric tensor A: σ = λ tr A δ + μA + γAᵀ
    A = np.arange(9.0).reshape(3, 3) - 3.0
    assert np.allclose(ch04.linear_stress(K4, A), 0.3 * np.trace(A) * np.eye(3) + 1.1 * A + 0.7 * A.T)
    rng = np.random.default_rng(2)
    for _ in range(20):
        Gm = rng.normal(size=(3, 3))
        S = 0.5 * (Gm + Gm.T)
        lam, mu = rng.uniform(-1, 1), rng.uniform(0.1, 2)
        sig = ch04.linear_stress(ch04.isotropic_fourth_order(lam, mu, mu), S)
        assert np.allclose(sig, ch04.viscous_stress(Gm, mu, ch04.bulk_viscosity(lam, mu)), atol=1e-12)
        # only μ + γ acts on a symmetric S
        assert np.allclose(ch04.linear_stress(ch04.isotropic_fourth_order(lam, 0.4, 2 * mu - 0.4), S), sig, atol=1e-12)
        # (4.31) with λ ≡ (4.37) with μ_v = λ + ⅔μ
        assert np.allclose(ch04.newtonian_stress(Gm, 1.0, mu, lam=lam), ch04.newtonian_stress(Gm, 1.0, mu,
                           mu_v=ch04.bulk_viscosity(lam, mu)), atol=1e-12)
        # rotating the axes rotates the stress: σ(CGCᵀ) = Cσ(G)Cᵀ
        C = TN.random_rotation(rng)
        assert np.allclose(ch04.viscous_stress(C @ Gm @ C.T, mu, 0.3), C @ ch04.viscous_stress(Gm, mu, 0.3) @ C.T, atol=1e-12)


def test_newtonian_V1_pressure_trace_identities():  # V1 (4.32)–(4.36)
    rng = np.random.default_rng(3)
    for _ in range(10):
        Gm = rng.normal(size=(3, 3))
        p, mu, lam = rng.uniform(1, 5), rng.uniform(0.1, 1), rng.uniform(-0.5, 0.5)
        tau = ch04.newtonian_stress(Gm, p, mu, lam=lam)
        div = np.trace(Gm)
        assert ch04.thermodynamic_pressure_from_stress(tau, div, mu, lam) == pytest.approx(p, rel=1e-12)
        assert p - ch04.mean_pressure(tau) == pytest.approx(ch04.pressure_difference(div, mu=mu, lam=lam), rel=1e-10, abs=1e-12)
        assert ch04.pressure_difference(div, mu_v=ch04.bulk_viscosity(lam, mu)) == pytest.approx(
            ch04.pressure_difference(div, mu=mu, lam=lam), rel=1e-12)
        assert ch04.lam_from_bulk(ch04.bulk_viscosity(lam, mu), mu) == pytest.approx(lam, abs=1e-14)
    assert ch04.stokes_assumption_holds(-2 / 3 * 0.9, 0.9) and not ch04.stokes_assumption_holds(0.0, 0.9)
    assert ch04.mean_pressure(ch04.newtonian_stress(rng.normal(size=(3, 3)), 4.0, 0.5)) == pytest.approx(4.0)  # μ_v = 0
    with pytest.raises(ValueError):
        ch04.newtonian_stress(np.eye(3), mu=1.0, incompressible=True)
    with pytest.raises(ValueError):
        ch04.newtonian_stress(np.eye(3), lam=0.1, mu_v=0.2)
    with pytest.raises(ValueError):
        ch04.pressure_difference(1.0, mu=1.0)
    Gs = ch04.stress_lab_gradient("extension", 2.0)
    assert np.allclose(ch04.newtonian_stress(Gs, 1.0, 0.5, incompressible=True), -np.eye(3) + 2 * 0.5 * Gs)  # (4.35)
    A = rng.normal(size=(3, 3))
    assert abs(np.trace(ch04.deviatoric_part(A))) < 1e-14


def test_stress_on_plane_V1_rotatable_plane():  # V1 (E3)
    G_ = [[0, 10, 0], [0, 0, 0], [0, 0, 0]]
    assert ch04.stress_on_plane(G_, 0, 1e-3, 0, np.pi / 4) == pytest.approx((0.010, 0.0), abs=1e-15)
    assert ch04.stress_on_plane(G_, 3.0, 1e-3, 0, 0.0) == pytest.approx((-3.0, 0.010), abs=1e-15)
    for th in np.linspace(0, np.pi, 7):
        sn, ts = ch04.stress_on_plane("rotation", 2.0, 1.0, 0.0, th)
        assert sn == pytest.approx(-2.0) and abs(ts) < 1e-14
        sn, ts = ch04.stress_on_plane(ch04.stress_lab_gradient("expansion", 0.5), 2.0, 1.0, 0.8, th)
        assert sn == pytest.approx(-2.0 + 3 * 0.8 * 0.5) and abs(ts) < 1e-14
        n = np.array([np.cos(th), np.sin(th), 0.0])
        tau = ch04.newtonian_stress(np.array(G_, float), 1.0, 2e-3)
        sn2 = TN.normal_shear_stress(tau, n)[0]
        assert ch04.stress_on_plane(G_, 1.0, 2e-3, 0.0, th)[0] == pytest.approx(sn2, abs=1e-14)
    with pytest.raises(ValueError):
        ch04.stress_lab_gradient("twist")


def test_stress_symmetry_V1_cube_spin_diverges_unless_symmetric():  # V1 / V3 (D08)
    hs = np.array([0.1, 0.05, 0.025, 0.0125])
    acc = ch04.cube_spin_acceleration(1.0, 0.0, 1000.0, hs)
    assert abs(observed_order(hs, acc) + 2.0) < 1e-12  # slope −2: diverges as h → 0
    assert ch04.cube_spin_acceleration(1.0, 0.0, 1000.0, 0.01) == pytest.approx(60.0)
    assert np.all(ch04.cube_spin_acceleration(0.7, 0.7, 1000.0, hs) == 0.0)


def test_stress_symmetry_V2_derivation():  # V2 — D08: torque (τ₁₂ − τ₂₁)h³ over I = ρh⁵/6
    x, y, z, h, rho, t12, t21 = sp.symbols("x y z h rho tau12 tau21", positive=True)
    I = sp.integrate(rho * (x ** 2 + y ** 2), (x, -h / 2, h / 2), (y, -h / 2, h / 2), (z, -h / 2, h / 2))
    assert sp.simplify(I - rho * h ** 5 / 6) == 0
    # torque about z of the shear forces on the four side faces (faces x = ±h/2 carry τ₁₂, y = ±h/2 carry τ₂₁)
    torque = 2 * (h / 2) * (t12 * h ** 2) - 2 * (h / 2) * (t21 * h ** 2)
    alpha = sp.simplify(torque / I)
    assert sp.simplify(alpha - 6 * (t12 - t21) / (rho * h ** 2)) == 0
    assert float(alpha.subs({t12: 1.0, t21: 0.0, rho: 1000.0, h: 0.01})) == pytest.approx(
        ch04.cube_spin_acceleration(1.0, 0.0, 1000.0, 0.01))


def test_newtonian_law_V2_derivation():  # V2 — D09 ★★★, re-running the construction of the design's check cell
    lam, mu, gam, th = sp.symbols("lambda mu gamma theta", real=True)
    d = sp.KroneckerDelta
    idx = range(3)
    K4 = {(i, j, m, n): lam * d(i, j) * d(m, n) + mu * d(i, m) * d(j, n) + gam * d(i, n) * d(j, m)
          for i in idx for j in idx for m in idx for n in idx}  # step 4, (4.29)
    s = sp.symbols("s11 s12 s13 s22 s23 s33", real=True)
    S = sp.Matrix([[s[0], s[1], s[2]], [s[1], s[3], s[4]], [s[2], s[4], s[5]]])  # symmetric S
    sig = sp.Matrix(3, 3, lambda i, j: sum(K4[i, j, m, n] * S[m, n] for m in idx for n in idx))  # steps 2, 5–7 (4.28)
    trS = S.trace()
    step8 = lam * trS * sp.eye(3) + (mu + gam) * S
    assert sp.simplify(sig - step8) == sp.zeros(3, 3)  # step 8: only λ and μ + γ survive
    assert sp.simplify(sig - sig.T) == sp.zeros(3, 3)  # σ symmetric whatever γ is
    assert sp.simplify(sig.subs(gam, mu) - (2 * mu * S + lam * trS * sp.eye(3))) == sp.zeros(3, 3)  # steps 9–10
    # the book's argument (4.30): K symmetric in i, j forces γ = μ
    asym = [sp.expand(K4[i, j, m, n] - K4[j, i, m, n]) for i in idx for j in idx for m in idx for n in idx]
    assert any(a_ != 0 for a_ in asym) and all(sp.expand(a_.subs(gam, mu)) == 0 for a_ in asym)
    # step 4's isotropy: K′_ijmn = C_ip C_jq C_mr C_ns K_pqrs = K_ijmn under a symbolic rotation about z
    C = sp.Matrix([[sp.cos(th), sp.sin(th), 0], [-sp.sin(th), sp.cos(th), 0], [0, 0, 1]])
    rng = np.random.default_rng(9)
    for _ in range(20):
        i, j, m, n = (int(v) for v in rng.integers(0, 3, 4))
        Kp = sum(C[i, p] * C[j, q] * C[m, r] * C[n, w] * K4[p, q, r, w] for p in idx for q in idx for r in idx for w in idx)
        assert sp.simplify(Kp - K4[i, j, m, n]) == 0
    # step 12: parallel shear recovers Newton's law; step 11: τ = −pδ + σ (4.31)
    gd, p = sp.symbols("gammadot p")
    Sh = sp.Matrix([[0, gd / 2, 0], [gd / 2, 0, 0], [0, 0, 0]])
    tau = -p * sp.eye(3) + 2 * mu * Sh + lam * Sh.trace() * sp.eye(3)
    assert sp.simplify(tau[0, 1] - mu * gd) == 0


def test_bulk_viscosity_V2_derivation():  # V2 — D10: (4.31) → (4.37) by adding and subtracting ⅔μS_mmδ_ij
    mu, lam, p = sp.symbols("mu lambda p", real=True)
    S = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"s{min(i, j)}{max(i, j)}"))
    tr = S.trace()
    t431 = -p * sp.eye(3) + 2 * mu * S + lam * tr * sp.eye(3)
    muv = lam + sp.Rational(2, 3) * mu
    t437 = -p * sp.eye(3) + 2 * mu * (S - tr / 3 * sp.eye(3)) + muv * tr * sp.eye(3)
    assert sp.simplify(t431 - t437) == sp.zeros(3, 3)
    assert sp.simplify((-(t431.trace()) / 3) - (p - muv * tr)) == 0  # p̄ = p − μ_v∇·u (4.33)–(4.34)
    assert sp.simplify(t437.subs(lam, -sp.Rational(2, 3) * mu).trace() + 3 * p) == 0  # Stokes (4.36): p = p̄


# =====================================================================================================================
# C08 — Navier–Stokes (4.38)–(4.41); exact solutions; D11–D13
# =====================================================================================================================
def _exact_points(name):
    rng = np.random.default_rng(21)
    if name in ("couette", "poiseuille"):
        return np.stack([rng.uniform(-0.01, 0.01, 12), rng.uniform(0.001, 0.009, 12)])
    if name == "pipe_poiseuille":
        return np.stack([rng.uniform(-0.006, 0.006, 12), rng.uniform(-0.006, 0.006, 12), rng.uniform(-1, 1, 12)])
    if name == "stokes_first":
        return np.stack([rng.uniform(-1, 1, 12), rng.uniform(1e-3, 0.02, 12)])
    if name == "lamb_oseen":
        r, th = rng.uniform(0.005, 0.06, 12), rng.uniform(0, 2 * np.pi, 12)
        return np.stack([r * np.cos(th), r * np.sin(th)])
    if name == "cylinder":
        r, th = rng.uniform(1.2, 3, 12), rng.uniform(0, 2 * np.pi, 12)
        return np.stack([r * np.cos(th), r * np.sin(th)])
    if name == "solid_body":
        return rng.uniform(-1, 1, (3, 12))
    return rng.uniform(-2, 2, (2, 12))  # taylor_green


_EXACT_CASES = {
    "couette": dict(U=1.0, h=0.01, G=50.0), "poiseuille": dict(G=100.0, h=0.01, g=9.81),
    "pipe_poiseuille": dict(G=100.0, R=0.01), "stokes_first": dict(U=1.0, nu=1e-4),
    "taylor_green": dict(U0=1.0, k=1.3, rho=1.0, mu=0.05), "lamb_oseen": dict(Gamma=1e-3, t0=50.0),
    "cylinder": dict(U=1.0, a=1.0), "solid_body": dict(Omega=0.7),
}


@pytest.mark.parametrize("name", list(_EXACT_CASES))
def test_navier_stokes_V1_exact_solutions_have_zero_residual(name):  # V1
    u_fn, p_fn, q = ch04.exact_solution_fields(name, **_EXACT_CASES[name])
    X = _exact_points(name)
    d = X.shape[0]
    L = {"couette": 0.01, "poiseuille": 0.01, "pipe_poiseuille": 0.01, "stokes_first": 0.01, "lamb_oseen": 0.02}.get(name, 1)
    t = {"stokes_first": 0.5, "lamb_oseen": 10.0, "taylor_green": 0.3}.get(name, 0.0)
    g = {"poiseuille": (0.0, -9.81), "pipe_poiseuille": (0.0, 0.0, -G), "solid_body": (0.0, 0.0, -G)}.get(name, (0.0,) * d)
    terms = ch04.ns_incompressible_terms(u_fn, p_fn, X, t, q["rho"], q["mu"], g, h=2e-4 * L, ht=1e-6, per="mass")
    scale = max(np.max(np.abs(getattr(terms, k))) for k in ("local", "advective", "pressure", "viscous", "gravity"))
    assert scale > 0 and np.max(np.abs(terms.residual)) < 1e-6 * scale, name
    ok, dv = ch04.divergence_free_check(u_fn, X, t, h=2e-4 * L, tol=1e-6 / L)
    assert ok, (name, dv)
    # the compressible, variable-μ form (4.38) with constant μ, μ_v = 0 gives the same zero
    r38 = ch04.navier_stokes_residual(q["rho"], u_fn, p_fn, X, t, q["mu"], 0.0, g, h=2e-4 * L, ht=1e-6)
    assert np.max(np.abs(r38)) < 1e-5 * scale * q["rho"], name


def test_navier_stokes_V2_symbolic_residuals_of_exact_solutions():  # V2 (independent sympy fields)
    x, y, z, t = sp.symbols("x y z t", real=True)
    U, h, G_, mu, rho, nu, k, U0, R, Om, g = sp.symbols("U h G mu rho nu k U0 R Omega g", positive=True)
    cases = [
        ([U * y / h + G_ / (2 * mu) * y * (h - y), 0], -G_ * x, rho, mu, [0, 0], (x, y)),
        ([0, 0, G_ / (4 * mu) * (R ** 2 - x ** 2 - y ** 2)], -(G_ + rho * g) * z, rho, mu, [0, 0, -g], (x, y, z)),
        ([U0 * sp.sin(k * x) * sp.cos(k * y) * sp.exp(-2 * nu * k ** 2 * t),
          -U0 * sp.cos(k * x) * sp.sin(k * y) * sp.exp(-2 * nu * k ** 2 * t)],
         rho * U0 ** 2 / 4 * (sp.cos(2 * k * x) + sp.cos(2 * k * y)) * sp.exp(-4 * nu * k ** 2 * t), rho, rho * nu, [0, 0],
         (x, y)),
        ([-Om * y, Om * x, 0], rho * Om ** 2 * (x ** 2 + y ** 2) / 2 - rho * g * z, rho, mu, [0, 0, -g], (x, y, z)),
        ([U * sp.erfc(y / (2 * sp.sqrt(nu * t))), 0], 0, rho, rho * nu, [0, 0], (x, y)),
    ]
    for u, p, r_, m_, gv, cs in cases:
        res = ch04.navier_stokes_sym(r_, u, p, cs, t, m_, 0, g=gv, form="4.39b")
        assert all(sp.simplify(e) == 0 for e in res)
        res38 = ch04.navier_stokes_sym(r_, u, p, cs, t, m_, 0, g=gv, form="4.38")
        assert all(sp.simplify(e) == 0 for e in res38)
    # discrimination: the wrong sign of the viscous term fails for Taylor–Green
    u, p, r_, m_, gv, cs = cases[2]
    wrong = [e + 2 * m_ * sum(sp.diff(u[j], c_, 2) for c_ in cs) for j, e in
             enumerate(ch04.navier_stokes_sym(r_, u, p, cs, t, m_, 0, g=gv, form="4.39b"))]
    assert any(sp.simplify(e) != 0 for e in wrong)


def test_navier_stokes_V2_forms_4_38_4_39a_4_39b_agree():  # V2 — D11, D12
    x, y, z, t = sp.symbols("x y z t", real=True)
    mu, muv = sp.symbols("mu mu_v", positive=True)
    rho = sp.Function("rho")(x, y, z, t)
    p = sp.Function("p")(x, y, z, t)
    u = [sp.Function(f"u{i}")(x, y, z, t) for i in range(3)]
    a = ch04.navier_stokes_sym(rho, u, p, (x, y, z), t, mu, muv, form="4.38")
    b = ch04.navier_stokes_sym(rho, u, p, (x, y, z), t, mu, muv, form="4.39a")
    assert all(sp.simplify(sp.expand(ai - bi)) == 0 for ai, bi in zip(a, b))  # constant μ, μ_v
    c = ch04.navier_stokes_sym(rho, u, p, (x, y, z), t, mu, muv, form="4.39b")
    divu = sum(sp.diff(u[m], v) for m, v in enumerate((x, y, z)))
    assert all(sp.simplify(sp.expand(bi - ci + (muv + mu / 3) * sp.diff(divu, v))) == 0 for bi, ci, v in zip(b, c, (x, y, z)))
    # variable μ(x): (4.38) keeps μ inside the derivative and differs from (4.39a)
    mux = sp.Function("mu")(x)
    a2 = ch04.navier_stokes_sym(rho, u, p, (x, y, z), t, mux, 0, form="4.38")
    b2 = ch04.navier_stokes_sym(rho, u, p, (x, y, z), t, mux, 0, form="4.39a")
    assert any(sp.simplify(sp.expand(ai - bi)) != 0 for ai, bi in zip(a2, b2))
    # D11: (4.38) is (4.24) with (4.37): ∂(−pδ_ij)/∂x_i = −∂p/∂x_j
    tau = [[-p * (1 if i == j else 0) + mu * (sp.diff(u[i], (x, y, z)[j]) + sp.diff(u[j], (x, y, z)[i]))
            + (muv - sp.Rational(2, 3) * mu) * divu * (1 if i == j else 0) for j in range(3)] for i in range(3)]
    cauchy = [rho * (sp.diff(u[j], t) + sum(u[i] * sp.diff(u[j], v) for i, v in enumerate((x, y, z))))
              - sum(sp.diff(tau[i][j], v) for i, v in enumerate((x, y, z))) for j in range(3)]
    assert all(sp.simplify(sp.expand(ci - ai)) == 0 for ci, ai in zip(cauchy, a))


def test_navier_stokes_V3_residual_order_2():  # V3
    u_fn, p_fn, q = ch04.exact_solution_fields("taylor_green", U0=1.0, k=1.0, rho=1.0, mu=0.1)
    X = np.random.default_rng(5).uniform(-2, 2, (2, 8))
    hs = [0.1, 0.05, 0.025, 0.0125]
    errs = [np.max(np.abs(ch04.ns_incompressible_terms(u_fn, p_fn, X, 0.2, 1.0, 0.1, (0, 0), h=h, ht=1e-4).residual))
            for h in hs]
    assert abs(observed_order(hs, errs) - 2.0) < ORDER_TOL, pairwise_orders(hs, errs)
    errs38 = [np.max(np.abs(ch04.navier_stokes_residual(1.0, u_fn, p_fn, X, 0.2, 0.1, 0.0, (0, 0), h=h, ht=1e-4)))
              for h in hs]
    assert abs(observed_order(hs, errs38) - 2.0) < ORDER_TOL, pairwise_orders(hs, errs38)


@needs_ref
def test_navier_stokes_V1_taylor_green_form_cross_check():  # V1 (form cross-check against Wikipedia, not V5)
    ref = ref_json()["taylor_green"]
    x, y, t, U0, k, nu, rho = sp.symbols("x y t U0 k nu rho")
    fu, fv, fp = (sp.lambdify((x, y, t, U0, k, nu, rho), sp.sympify(ref[c]), "numpy") for c in ("u", "v", "p"))
    X = np.random.default_rng(6).uniform(-2, 2, (2, 30))
    for (U0v, kv, rhov, muv, tv) in [(1.0, 1.0, 1000.0, 1e-3, 5.0), (2.0, 3.0, 1.0, 0.2, 0.7)]:
        vel, p = ch04.exact_solution("taylor_green", X, tv, U0=U0v, k=kv, rho=rhov, mu=muv)
        nuv = muv / rhov
        assert rel(vel[0], fu(X[0], X[1], tv, U0v, kv, nuv, rhov)) < 1e-12
        assert rel(vel[1], fv(X[0], X[1], tv, U0v, kv, nuv, rhov)) < 1e-12
        assert rel(p, fp(X[0], X[1], tv, U0v, kv, nuv, rhov)) < 1e-12


def test_navier_stokes_V1_lamb_oseen_pressure_from_radial_balance():  # V1 (quad, independent of the E₁ closed form)
    Gam, rho, nu, t = 2e-3, 1000.0, 1e-6, 30.0
    for r in (0.004, 0.01, 0.03):
        vel, p = ch04.exact_solution("lamb_oseen", np.array([[r], [0.0]]), t, Gamma=Gam, rho=rho, mu=rho * nu)
        ut = lambda s: Gam / (2 * np.pi * s) * (-np.expm1(-s ** 2 / (4 * nu * t)))  # noqa: E731
        pq = -quad(lambda s: rho * ut(s) ** 2 / s, r, np.inf, epsabs=0, epsrel=1e-12, limit=200)[0]
        assert float(np.ravel(p)[0]) == pytest.approx(pq, rel=1e-9)
        assert float(vel[1, 0]) == pytest.approx(ut(r), rel=1e-12)
    # and the Gaussian-vortex profile of ch03 with σ² = 4νt
    from fluidpy.core.vortices import gaussian_vortex
    vel, _ = ch04.exact_solution("lamb_oseen", np.array([[0.01], [0.0]]), t, Gamma=Gam, rho=rho, mu=rho * nu)
    assert float(vel[1, 0]) == pytest.approx(float(gaussian_vortex(0.01, Gam, np.sqrt(4 * nu * t))[0]), rel=1e-12)


def test_ns_terms_preset_V1_poiseuille_balance_and_forms():  # V1 (E4 contract)
    d = ch04.ns_terms_preset("poiseuille", 0.0, 2.5e-4, component=0, per="volume", G=100.0, h=1e-3, mu=1e-3)
    assert d["pressure"] == pytest.approx(100.0, rel=1e-9) and d["viscous"] == pytest.approx(-100.0, rel=1e-6)
    assert abs(d["local"]) < 1e-9 and abs(d["advective"]) < 1e-9 and abs(d["residual"]) < 1e-4
    for form in ("laplacian", "div2S", "curl"):
        e = ch04.ns_terms_preset("taylor_green", 0.3, 0.2, t=0.1, component=0, form=form, rho=1.0, mu=0.05)
        assert abs(e["residual"]) < 1e-6 * abs(e["advective"]), form
    with pytest.raises(ValueError):
        ch04.ns_terms_preset("taylor_green", 0.3, 0.2, form="weird")
    with pytest.raises(ValueError):
        ch04.exact_solution("nope", np.zeros((2, 1)))
    # N107: gravity along −y in the channel is absorbed by the pressure: same velocity with and without g
    a = ch04.exact_solution("poiseuille", np.array([[0.1], [0.004]]), G=100.0, h=0.01, g=9.81)
    b = ch04.exact_solution("poiseuille", np.array([[0.1], [0.004]]), G=100.0, h=0.01, g=0.0)
    assert np.allclose(a[0], b[0]) and float(np.ravel(a[1] - b[1])[0]) == pytest.approx(-1000.0 * 9.81 * 0.004)
    assert ch04.stokes_first_problem(0.01, 2.0, 1.0, 1e-4) == pytest.approx(float(erfc(0.01 / (2 * np.sqrt(2e-4)))))
    assert ch04.plane_poiseuille(0.005, 100.0, 0.01, 1e-3) == pytest.approx(100.0 * 0.005 ** 2 / (2e-3))


def test_viscous_force_V2_three_forms_and_compressible_difference():  # V2 + V1 (4.40); D13
    x, y, z = sp.symbols("x y z", real=True)
    coords = (x, y, z)
    A = [x ** 2 * y * z, sp.sin(x) * z + y ** 3, x * y * z ** 2]
    u = [sp.diff(A[2], y) - sp.diff(A[1], z), sp.diff(A[0], z) - sp.diff(A[2], x), sp.diff(A[1], x) - sp.diff(A[0], y)]
    assert sp.simplify(sum(sp.diff(u[i], coords[i]) for i in range(3))) == 0  # divergence-free (a curl)
    lap = [sum(sp.diff(u[j], c, 2) for c in coords) for j in range(3)]
    div2S = [sum(sp.diff(sp.diff(u[j], coords[i]) + sp.diff(u[i], coords[j]), coords[i]) for i in range(3)) for j in range(3)]
    om = [sp.diff(u[2], y) - sp.diff(u[1], z), sp.diff(u[0], z) - sp.diff(u[2], x), sp.diff(u[1], x) - sp.diff(u[0], y)]
    curl_om = [sp.diff(om[2], y) - sp.diff(om[1], z), sp.diff(om[0], z) - sp.diff(om[2], x), sp.diff(om[1], x) - sp.diff(om[0], y)]
    assert all(sp.simplify(lap[j] - div2S[j]) == 0 and sp.simplify(lap[j] + curl_om[j]) == 0 for j in range(3))
    # numeric: the three stencil forms agree on this field; a compressible field shows ∓μ∇(∇·u)
    uf = vec_field(u, coords)
    X = rand_pts(6, lo=-0.8, hi=0.8)
    l_, d_, c_ = ch04.viscous_force_forms(uf, X, mu=0.7, h=1e-3)
    assert rel(d_, l_) < 1e-5 and rel(c_, l_) < 1e-5
    comp = [x ** 2 * y, y * z ** 2, sp.sin(x) * z]
    cf = vec_field(comp, coords)
    l2, d2, c2 = ch04.viscous_force_forms(cf, X, mu=0.7, h=1e-3)
    gdiv = sp.lambdify(coords, [sp.diff(sum(sp.diff(comp[i], coords[i]) for i in range(3)), c) for c in coords])
    gd = 0.7 * np.array(gdiv(*X))
    assert rel(d2 - l2, gd) < 1e-5 and rel(c2 - l2, -gd) < 1e-5
    # solid-body rotation: all three are 0 (the "paradox": vorticity without viscous force)
    sb = lambda X_, T: np.stack([-X_[1], X_[0], 0 * X_[0]])  # noqa: E731
    assert all(np.max(np.abs(f)) < 1e-7 for f in ch04.viscous_force_forms(sb, X, mu=1.0, h=1e-3))
    # 2-D plane field branch: the forms agree for Taylor–Green
    tg, _ = ch04.exact_field("taylor_green", U0=1.0, k=1.0)
    l3, d3, c3 = ch04.viscous_force_forms(tg, np.array([[0.3, -0.7], [0.2, 1.1]]), mu=1.0, h=1e-3)
    assert rel(d3, l3) < 1e-5 and rel(c3, l3) < 1e-5


def test_viscous_force_V2_derivation():  # V2 — D13 with the ε–δ identity for a generic field
    x, y, z = sp.symbols("x y z", real=True)
    coords = (x, y, z)
    u = [sp.Function(f"u{i}")(x, y, z) for i in range(3)]
    divu = sum(sp.diff(u[i], coords[i]) for i in range(3))
    eps = sp.LeviCivita
    om = [sum(eps(k, m, n) * sp.diff(u[n], coords[m]) for m in range(3) for n in range(3)) for k in range(3)]
    for j in range(3):
        two_dS = sum(sp.diff(sp.diff(u[j], coords[i]) + sp.diff(u[i], coords[j]), coords[i]) for i in range(3))
        lap = sum(sp.diff(u[j], c, 2) for c in coords)
        assert sp.simplify(two_dS - lap - sp.diff(divu, coords[j])) == 0  # 2∂S_ij/∂x_i = ∇²u_j + ∂_j(∇·u)
        curl = sum(eps(j, i, k) * sp.diff(om[k], coords[i]) for i in range(3) for k in range(3))
        assert sp.simplify(sp.expand(-curl - lap + sp.diff(divu, coords[j]))) == 0  # −∇×ω = ∇²u − ∇(∇·u)


def test_navier_stokes_V1_compressible_4_38_with_bulk_viscosity():  # V1 (4.38)/(4.39a) numerics vs sympy
    coords = (X1, X2, X3)
    uexpr = [X1 ** 2 * X2, sp.sin(X2) * X3, X1 * X3]
    rhoexpr = 1 + X1 ** 2
    pexpr = X1 * X2 + X3
    mu, muv = 0.3, 0.2
    u = vec_field(uexpr, coords)
    rhs = ch04.navier_stokes_sym(rhoexpr, uexpr, pexpr, coords, TT, mu, muv, g=[0, 0, -G], form="4.38")
    ex = sp.lambdify(coords, rhs)
    X = rand_pts(8, lo=-0.7, hi=0.7)
    num = ch04.navier_stokes_residual(sca_field(rhoexpr, coords), u, sca_field(pexpr, coords), X, 0.0, mu, muv,
                                      (0, 0, -G), h=1e-3, ht=1e-5)
    assert rel(num, np.array(ex(*X), float)) < 1e-5
    # variable viscosity μ(x) enters inside the derivative
    muf = lambda X_, T: 0.3 + 0.1 * X_[0]  # noqa: E731
    rhs2 = ch04.navier_stokes_sym(rhoexpr, uexpr, pexpr, coords, TT, 0.3 + 0.1 * X1, 0, g=[0, 0, -G], form="4.38")
    num2 = ch04.navier_stokes_residual(sca_field(rhoexpr, coords), u, sca_field(pexpr, coords), X, 0.0, muf, 0.0,
                                       (0, 0, -G), h=1e-3, ht=1e-5)
    assert rel(num2, np.array(sp.lambdify(coords, rhs2)(*X), float)) < 1e-5


# =====================================================================================================================
# C09 — noninertial frames (4.42)–(4.45); Coriolis, centrifugal, effective gravity; D14–D18
# =====================================================================================================================
def test_rotating_basis_V3_turning_rule_order_2():  # V3 + V1 (D14: de′/dt = Ω × e′)
    Om = np.array([0.3, -0.5, 0.8])
    E = ch04.rotating_basis(Om, 1.7)
    assert np.allclose(E @ E.T, np.eye(3), atol=1e-14) and np.linalg.det(E) == pytest.approx(1.0)
    exact = ch04.basis_rate_exact(Om, E)
    hs = [0.1, 0.05, 0.025, 0.0125]
    errs = [np.max(np.abs(ch04.basis_rate(Om, 1.7, h=h) - exact)) for h in hs]
    assert abs(observed_order(hs, errs) - 2.0) < ORDER_TOL
    assert np.allclose(ch04.rotating_basis(0.0, 3.0), np.eye(3))
    assert np.allclose(ch04.rotating_basis(0.5, 2 * np.pi)[0], [-1, 0, 0], atol=1e-14)  # e′₁ turned by Ωt = π about z
    assert np.allclose(ch04.rotating_basis(0.5, np.pi)[0], [0, 1, 0], atol=1e-14)  # … and by π/2
    xp = np.array([1.0, 2.0, 0.5])
    assert np.allclose(ch04.inertial_velocity(np.zeros(3), np.zeros(3), Om, xp), np.cross(Om, xp))  # (4.42)


def _frame_case(seed):
    """A path given in the inertial frame and a frame moving/rotating about a fixed axis; exact derivatives by sympy."""
    rng = np.random.default_rng(seed)
    t = sp.Symbol("t", real=True)
    ax = rng.normal(size=3)
    ax /= np.linalg.norm(ax)
    W0, al = rng.uniform(-1, 1), rng.uniform(-0.5, 0.5)
    th = W0 * t + al * t ** 2 / 2
    n = sp.Matrix(ax)
    Kx = sp.Matrix([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
    R = sp.eye(3) + sp.sin(th) * Kx + (1 - sp.cos(th)) * Kx * Kx  # Rodrigues (columns = e′_i in inertial comps)
    c = rng.uniform(-1, 1, (3, 3))
    xin = sp.Matrix([c[0, 0] + c[0, 1] * t + c[0, 2] * sp.sin(t), c[1, 0] * t ** 2 + c[1, 1], c[2, 0] * sp.cos(2 * t) + c[2, 2] * t])
    Xo = sp.Matrix([0.3 * t ** 2, -0.2 * t, sp.sin(0.5 * t)])
    xp = R.T * (xin - Xo)  # primed components
    tv = float(rng.uniform(0, 2))
    ev = lambda e: np.array(sp.N(e.subs(t, tv)), dtype=float).ravel()  # noqa: E731
    return dict(a_true=ev(R.T * sp.diff(xin, t, 2)), a_p=ev(sp.diff(xp, t, 2)), u_p=ev(sp.diff(xp, t)), x_p=ev(xp),
                Om=ax * float(W0 + al * tv), dOm=ax * al, dU=ev(R.T * sp.diff(Xo, t, 2)))


def test_frame_acceleration_V1_random_paths_reassemble_inertial_acceleration():  # V1 (4.43)
    for seed in range(12):
        c = _frame_case(seed)
        terms = ch04.frame_acceleration_terms(c["a_p"], c["u_p"], c["x_p"], c["Om"], c["dOm"], c["dU"])
        assert rel(terms["total"], c["a_true"]) < 1e-10, seed
        wrong = terms["total"] - 0.5 * terms["coriolis"]  # "Coriolis without the 2"
        if np.linalg.norm(terms["coriolis"]) > 1e-3:
            assert rel(wrong, c["a_true"]) > 1e-3
        forces = ch04.apparent_body_forces(c["u_p"], c["x_p"], c["Om"], c["dOm"], c["dU"], g=(0, 0, -G))
        # (4.45) bracket = g − (frame + Coriolis + angular + centripetal accelerations)
        expect = np.array([0, 0, -G]) - (terms["total"] - terms["relative"])
        assert rel(forces["total"], expect) < 1e-12
    f0 = ch04.apparent_body_forces(np.array([1.0, 2.0, 3.0]), np.array([0.5, 0, 0]), 0.0)
    assert np.allclose(f0["total"], [0, 0, -G])  # inertial frame: gravity only


def test_rotating_frame_V2_derivation():  # V2 — D15 ★★★ as in the design's check cell (and D14 on the way)
    t = sp.Symbol("t", real=True)
    th = sp.Function("theta")(t)
    R = sp.Matrix([[sp.cos(th), -sp.sin(th), 0], [sp.sin(th), sp.cos(th), 0], [0, 0, 1]])  # e′ columns
    Xo = sp.Matrix([sp.Function(f"X{i}")(t) for i in range(3)])
    xp = sp.Matrix([sp.Function(f"x{i}")(t) for i in range(3)])
    x = Xo + R * xp  # D14 step 1
    Om = sp.Matrix([0, 0, sp.diff(th, t)])
    dOm = sp.diff(Om, t)
    cr = lambda a, b: a.cross(b)  # noqa: E731
    # D14 (4.42): Rᵀẋ = RᵀẊ + ẋ′ + Ω × x′
    vel = sp.simplify(R.T * sp.diff(x, t) - (R.T * sp.diff(Xo, t) + sp.diff(xp, t) + cr(Om, xp)))
    assert vel == sp.zeros(3, 1)
    a = R.T * sp.diff(x, t, 2)
    five = R.T * sp.diff(Xo, t, 2) + sp.diff(xp, t, 2) + 2 * cr(Om, sp.diff(xp, t)) + cr(dOm, xp) + cr(Om, cr(Om, xp))
    assert sp.simplify(a - five) == sp.zeros(3, 1)  # (4.43)
    one = five - cr(Om, sp.diff(xp, t))  # Coriolis with factor 1
    assert sp.simplify(a - one - cr(Om, sp.diff(xp, t))) == sp.zeros(3, 1) and sp.simplify(a - one) != sp.zeros(3, 1)
    # step 10: triple product Ω × (Ω × x′) = Ω(Ω·x′) − Ω²x′; for Ω along z it is −Ω²R e_R
    tp = sp.simplify(cr(Om, cr(Om, xp)) - (Om * (Om.dot(xp)) - Om.dot(Om) * xp))
    assert tp == sp.zeros(3, 1)


def test_coriolis_V1_projectile_paths_and_right_deflection():  # V1 + V4 + V7 (Fig. 4.8, D17)
    d = ch04.coriolis_projectile(10.0, ch04.OMEGA_EARTH, [3600.0])
    assert d["forward"][0] == pytest.approx(36000.0) and d["deflection_small"][0] == pytest.approx(9450.6, abs=0.05)
    assert d["deflection"][0] == pytest.approx(9342.4, abs=0.05) and d["angle"][0] == pytest.approx(0.26252, abs=5e-6)
    W, u0 = 0.3, 2.0
    tt = np.linspace(0, 6, 61)
    r = ch04.coriolis_projectile(u0, W, tt)
    _, rot = ch04.projectile_paths(u0, W, tt)
    assert np.max(np.abs(r["rotating"] - rot)) < 1e-9 * u0 * 6  # ODE path = inertial straight line seen turning
    assert np.all(r["rotating"][1, 1:] < 0)  # NH (Ω > 0): deflected to the right (−y′ for motion along +x′)
    sh = ch04.coriolis_projectile(u0, -W, tt)
    assert np.all(sh["rotating"][1, 1:] > 0)  # SH: to the left
    nc = ch04.coriolis_projectile(u0, W, tt, centrifugal=False)
    assert np.ptp(nc["speed"]) < 1e-9 * u0  # the Coriolis force does no work
    cor = ch04.coriolis_acceleration([0, 0, W], [u0, 0, 0])
    assert np.allclose(cor, [0, -2 * W * u0, 0]) and np.allclose(ch04.coriolis_force([0, 0, W], [u0, 0, 0]), cor)
    assert np.allclose(ch04.frame_acceleration_terms(np.zeros(3), [u0, 0, 0], np.zeros(3), W)["coriolis"], -cor)
    # D17: deflection → Ωut² as Ωt → 0 with relative error (Ωt)²/6
    ts = np.array([0.4, 0.2, 0.1, 0.05])
    ratio = np.array([ch04.coriolis_deflection(u0, W, t_, exact=True) / ch04.coriolis_deflection(u0, W, t_) for t_ in ts])
    assert rel(1 - ratio, (W * ts) ** 2 / 6) < 0.02
    assert abs(observed_order(ts, 1 - ratio) - 2.0) < ORDER_TOL


def test_coriolis_V2_derivation_small_time_deflection():  # V2 — D17
    u, W, t = sp.symbols("u Omega t", positive=True)
    exact = u * t * sp.sin(W * t)
    assert sp.simplify(sp.series(exact, t, 0, 4).removeO() - W * u * t ** 2) == 0
    assert sp.simplify((W * u * t ** 2) / (u * t) - W * t) == 0  # angle = deflection / distance


def test_centrifugal_V2_potential_and_effective_gravity():  # V2 + V1 (D18)
    x, y, z, W = sp.symbols("x y z Omega", real=True)
    Om = sp.Matrix([0, 0, W])
    X = sp.Matrix([x, y, z])
    cen = -Om.cross(Om.cross(X))
    Phi = -sp.Rational(1, 2) * W ** 2 * (x ** 2 + y ** 2)
    assert sp.simplify(cen + sp.Matrix([sp.diff(Phi, v) for v in (x, y, z)])) == sp.zeros(3, 1)
    assert ch04.centrifugal_potential(2.0, 0.5) == pytest.approx(-0.5)
    pts = rand_pts(10)
    num = ch04.centrifugal_acceleration([0, 0, 0.7], pts)
    assert rel(num, 0.49 * np.stack([pts[0], pts[1], 0 * pts[0]])) < 1e-14  # +Ω²R e_R, away from the axis
    # effective gravity: independent vector sum in a meridional plane
    for lat in np.deg2rad([0.0, 20.0, 45.0, 70.0, 90.0]):
        rhat = np.array([np.cos(lat), np.sin(lat)])
        ge = -9.8 * rhat + ch04.OMEGA_EARTH ** 2 * ch04.EARTH_A * np.cos(lat) * np.array([1.0, 0.0])
        mag, dev = ch04.effective_gravity(lat)
        assert mag == pytest.approx(np.linalg.norm(ge), rel=1e-14)
        assert dev == pytest.approx(np.arccos(np.clip(-ge @ rhat / np.linalg.norm(ge), -1, 1)), abs=1e-12)
    lats = np.deg2rad(np.linspace(0, 90, 9001))
    peak = np.rad2deg(lats[np.argmax(ch04.effective_gravity(lats)[1])])
    assert 44.0 < peak < 46.0
    assert ch04.effective_gravity(0.0)[0] == pytest.approx(9.8 - ch04.OMEGA_EARTH ** 2 * ch04.EARTH_A)
    assert ch04.coriolis_parameter(np.pi / 2) == pytest.approx(2 * ch04.OMEGA_EARTH)
    assert ch04.coriolis_parameter(-np.pi / 6) == pytest.approx(-ch04.OMEGA_EARTH)


@needs_ref
def test_earth_rotation_V5_wgs84_numbers():  # V5
    w = ref_json()["wgs84"]
    assert ch04.OMEGA_EARTH == pytest.approx(w["omega_rad_s"], rel=1e-12)
    assert ch04.EARTH_A == w["a_m"]
    b = w["a_m"] * (1 - 1 / w["inv_f"])
    assert ch04.EARTH_B == pytest.approx(b, abs=1e-3) and b == pytest.approx(w["b_m"], abs=1e-5)
    assert ch04.earth_oblateness_diameter() == pytest.approx(2 * (w["a_m"] - w["b_m"]), abs=1e-3)
    assert ch04.earth_oblateness_diameter() / 1e3 == pytest.approx(42.77, abs=0.005)
    eq = ch04.centrifugal_acceleration([0, 0, w["omega_rad_s"]], [w["a_m"], 0, 0])
    assert eq[0] == pytest.approx(w["omega_rad_s"] ** 2 * w["a_m"], rel=1e-14)
    assert eq[0] == pytest.approx(0.033916, abs=5e-7)


def test_rotating_ns_V2_derivation():  # V2 — D16: ∇²(U + Ω × x′) = 0 and the forces of (4.45)
    x, y, z = sp.symbols("x y z", real=True)
    Om = sp.Matrix(sp.symbols("W1 W2 W3"))
    U = sp.Matrix(sp.symbols("U1 U2 U3"))
    frame_u = U + Om.cross(sp.Matrix([x, y, z]))
    assert all(sp.simplify(sum(sp.diff(frame_u[j], v, 2) for v in (x, y, z))) == 0 for j in range(3))
    # moving the kinematic terms to the right flips their signs: the bracket of (4.45)
    up, xp = np.array([0.3, -1.0, 0.2]), np.array([1.0, 0.5, -0.4])
    Omv, dOm, dU = np.array([0.1, 0.2, 0.7]), np.array([0.01, 0.0, 0.02]), np.array([0.5, 0.0, -0.1])
    a = ch04.frame_acceleration_terms(np.zeros(3), up, xp, Omv, dOm, dU)
    f = ch04.apparent_body_forces(up, xp, Omv, dOm, dU)
    for ka, kf in (("frame", "frame"), ("coriolis", "coriolis"), ("angular", "angular"), ("centripetal", "centrifugal")):
        assert np.allclose(f[kf], -a[ka])


def test_curvilinear_V1_operators_match_cartesian():  # V1 (Appendix B via core.curvilinear)
    CU = ch04.CU
    R, ph, z = CU.coordinates("cylindrical")
    r, th, ph2 = CU.coordinates("spherical")
    x, y, zz = CU.coordinates("cartesian")
    # a scalar in cylindrical coordinates and its Cartesian twin
    f_cyl = R ** 2 * sp.cos(ph) * z + R ** 3
    f_car = (x ** 2 + y ** 2) * (x / sp.sqrt(x ** 2 + y ** 2)) * zz + (x ** 2 + y ** 2) ** sp.Rational(3, 2)
    lap_car = sum(sp.diff(f_car, v, 2) for v in (x, y, zz))
    pt = dict(R=1.3, phi=0.7, z=0.4)
    xv, yv, zv = pt["R"] * np.cos(pt["phi"]), pt["R"] * np.sin(pt["phi"]), pt["z"]
    assert float(CU.laplacian(f_cyl, "cylindrical").subs({R: 1.3, ph: 0.7, z: 0.4})) == pytest.approx(
        float(lap_car.subs({x: xv, y: yv, zz: zv})), rel=1e-12)
    # solid-body rotation u = Ω × x: cylindrical (0, ΩR, 0): div 0, curl (0, 0, 2Ω), vector Laplacian 0, (u·∇)u = −Ω²R e_R
    W = sp.Symbol("W", positive=True)
    u = [0, W * R, 0]
    assert sp.simplify(CU.divergence(u, "cylindrical")) == 0
    assert [sp.simplify(c) for c in CU.curl(u, "cylindrical")] == [0, 0, 2 * W]
    assert [sp.simplify(c) for c in CU.vector_laplacian(u, "cylindrical")] == [0, 0, 0]
    assert [sp.simplify(c) for c in CU.advective_acceleration(u, "cylindrical")] == [-W ** 2 * R, 0, 0]
    assert sp.simplify(CU.strain_rate(u, "cylindrical")) == sp.zeros(3, 3)  # rigid rotation: no strain
    # spherical: the potential flow past a sphere has zero divergence and curl; ∇²(1/r) = 0
    U, a = sp.symbols("U a", positive=True)
    phi = U * (r + a ** 3 / (2 * r ** 2)) * sp.cos(th)
    us = CU.gradient(phi, "spherical")
    assert sp.simplify(CU.divergence(us, "spherical")) == 0
    assert all(sp.simplify(c) == 0 for c in CU.curl(us, "spherical"))
    assert sp.simplify(CU.laplacian(1 / r, "spherical")) == 0
    tt = sp.Symbol("t")
    ma = CU.material_acceleration([sp.Function("a")(tt), 0, 0], tt, "cylindrical")
    assert sp.simplify(ma[0] - sp.diff(sp.Function("a")(tt), tt)) == 0


def test_rotating_pump_V2_example_4_5_rotation_terms():  # V2 (N65)
    d = ch04.rotating_pump_terms()
    s = d["symbols"]
    rho, Om, uR, uphi, R = s["rho"], s["Omega_z"], s["u_R"], s["u_phi"], s["R"]
    expect = [rho * (2 * Om * uphi + Om ** 2 * R), -2 * rho * Om * uR, 0]
    assert all(sp.simplify(a_ - b_) == 0 for a_, b_ in zip(d["rotation_terms"], expect))
    eqs = ch04.rotating_pump_equations()
    assert len(eqs) == 3
    # the axisymmetric advective terms carry the curvature terms −u_φ²/R and u_Ru_φ/R
    lhs0 = sp.expand(d["lhs"][0] / rho)
    assert sp.simplify(lhs0 - (uR * sp.diff(uR, R) + s["u_z"] * sp.diff(uR, s["z"]) - uphi ** 2 / R)) == 0
    lhs1 = sp.expand(d["lhs"][1] / rho)
    assert sp.simplify(lhs1 - (uR * sp.diff(uphi, R) + s["u_z"] * sp.diff(uphi, s["z"]) + uR * uphi / R)) == 0


def test_high_low_flow_V1_radial_and_coriolis_turning():  # V1 (N62)
    u, v = ch04.high_low_flow(np.array([3.0, 0.0]), np.array([4.0, -2.0]), 2.0, "high")
    assert np.allclose(u, [1.2, 0.0]) and np.allclose(v, [1.6, -2.0])
    ul, vl = ch04.high_low_flow(3.0, 4.0, 2.0, "low")
    assert (ul, vl) == pytest.approx((-1.2, -1.6))
    assert ch04.high_low_flow(0.0, 0.0) == (0.0, 0.0)
    cor = ch04.coriolis_acceleration([0, 0, 1e-4], [1.2, 1.6, 0])  # out of a high in the NH: turned clockwise
    assert np.cross([1.2, 1.6, 0], cor)[2] < 0
    with pytest.raises(ValueError):
        ch04.high_low_flow(1.0, 1.0, sense="middle")


# =====================================================================================================================
# C10 — energy (4.46)–(4.63); dissipation; Couette heating; D19–D23
# =====================================================================================================================
def test_energy_V2_identity_chain():  # V2 — D19–D22, (4.54), (4.60) ⇔ (4.112)
    res = ch04.energy_identities_sym()
    assert set(res) == {"4.53_to_4.55", "4.24_dot_u_to_4.56", "4.55_minus_4.56_to_4.57", "4.54_split", "4.60_to_4.112"}
    assert all(v == 0 for v in res.values()), res
    assert ch04.energy_forms_sym() == 0


def test_energy_V1_total_energy_residual_on_exact_couette_heating():  # V1 (4.53) with an exact field
    y, U, h, mu, k, cv, rho, T0 = sp.symbols("y U h mu k c_v rho T0", positive=True)
    x, z, t = sp.symbols("x z t", real=True)
    T = T0 + mu * U ** 2 / (2 * k) * (y / h) * (1 - y / h)
    res = ch04.total_energy_residual_sym(rho=rho, u=[U * y / h, 0, 0], e=cv * T, p=sp.Symbol("p0"),
                                         q=[0, -k * sp.diff(T, y), 0], coords=(x, y, z), t=t, mu=mu, mu_v=0,
                                         g=[0, 0, -sp.Symbol("g")])
    assert sp.simplify(res) == 0
    wrong = T0 + mu * U ** 2 / k * (y / h) * (1 - y / h)  # twice the heating: not a solution
    res2 = ch04.total_energy_residual_sym(rho=rho, u=[U * y / h, 0, 0], e=cv * wrong, p=sp.Symbol("p0"),
                                          q=[0, -k * sp.diff(wrong, y), 0], coords=(x, y, z), t=t, mu=mu, mu_v=0,
                                          g=[0, 0, 0])
    assert sp.simplify(res2) != 0


def test_couette_heating_V1_profile_satisfies_4_60():  # V1
    U, h, mu, k, rho, cv = 2.0, 1e-3, 1e-3, 0.6, 1000.0, 4180.0
    Tf = lambda X, t: ch04.couette_heating(X[1], U, h, mu, k)["T"]  # noqa: E731
    uf = lambda X, t: np.stack([U * X[1] / h, 0 * X[1]])  # noqa: E731
    X = np.stack([np.linspace(-1e-3, 1e-3, 7), np.linspace(1e-4, 9e-4, 7)])
    r = ch04.internal_energy_residual(rho, uf, lambda X, t: cv * Tf(X, t), 1e5, Tf, mu, 0.0, k, X, 0.0, h=5e-5, ht=1e-4)
    scale = mu * (U / h) ** 2
    assert np.max(np.abs(r)) < 1e-6 * scale
    c = ch04.couette_heating(np.linspace(0, h, 5), U, h, mu, k)
    assert c["dT_max"] == pytest.approx(mu * U ** 2 / (8 * k), rel=1e-12)
    assert c["T"][2] - 293.15 == pytest.approx(mu * U ** 2 / (8 * k), rel=1e-12)
    assert np.allclose(c["eps"], mu * (U / h) ** 2)
    d = ch04.couette_heating(np.linspace(0, 1e-3, 5), U=1, h=1e-3, mu=1e-3, k=0.6)
    assert d["dT_max"] == pytest.approx(2.0833e-4, rel=1e-4) and d["heat_out"] == pytest.approx(1.0, rel=1e-12)


def test_couette_heating_V4_work_in_equals_heat_out():  # V4 (also with a pressure gradient)
    for dpdx in (0.0, -2e3, 5e3):
        U, h, mu, k = 1.5, 2e-3, 2e-3, 0.5
        c = ch04.couette_heating(np.linspace(0, h, 3), U, h, mu, k, dpdx=dpdx)
        diss = quad(lambda y: float(ch04.couette_heating(y, U, h, mu, k, dpdx=dpdx)["eps"]), 0, h, epsrel=1e-13)[0]
        assert c["work_in"] == pytest.approx(c["heat_out"], rel=1e-10)
        assert c["heat_out"] == pytest.approx(diss, rel=1e-10)
    # energy_budget on a box spanning the gap: shear work in at the moving wall = heat conducted out
    U, h, mu, k, rho, cv = 1.0, 1e-2, 1e-3, 0.6, 1000.0, 4180.0
    Tfun = lambda yy: ch04.couette_heating(yy, U, h, mu, k)["T"]  # noqa: E731
    u = lambda X, t: np.stack([U * X[1] / h, 0 * X[1], 0 * X[1]])  # noqa: E731
    e = lambda X, t: cv * np.asarray(Tfun(X[1]))  # noqa: E731
    q = lambda X, t: np.stack([0 * X[1], -k * (mu / k) * (U / h) ** 2 * (h / 2 - X[1]), 0 * X[1]])  # noqa: E731
    tau = lambda X, t: np.array([[-1e5 + 0 * X[1], mu * U / h + 0 * X[1], 0 * X[1]],  # noqa: E731
                                 [mu * U / h + 0 * X[1], -1e5 + 0 * X[1], 0 * X[1]],
                                 [0 * X[1], 0 * X[1], -1e5 + 0 * X[1]]])
    box = ch04.MovingBox(lengths=(1.0, h, 1.0), origin=(0.0, 0.0, 0.0))
    b = ch04.energy_budget(rho, u, e, box, 0.0, g=(0, 0, -G), q=q, tau=tau)
    assert b.surface_work == pytest.approx(mu * U ** 2 / h, rel=1e-10)
    assert b.heat_out == pytest.approx(mu * U ** 2 / h, rel=1e-10)
    assert abs(b.residual) < 1e-9 * b.heat_out


def test_energy_budget_V1_uniform_adiabatic_flow_and_material_cv():  # V1
    box = ch04.MovingBox(lengths=(1.0, 0.5, 0.5), origin=(0, 0, 0))
    uu = lambda X, t: np.broadcast_to(np.array([[2.0], [0.0], [0.0]]), X.shape)  # noqa: E731
    b = ch04.energy_budget(1.2, uu, 2e5, box, 0.0, g=(0, 0, -G), tau=lambda X, t: -1e5 * np.eye(3)[..., None] + 0 * X[0])
    assert max(abs(b.storage), abs(b.outflux), abs(b.body_work), abs(b.residual)) < 1e-6
    bm = ch04.energy_budget(1.2, uu, 2e5, box, 0.0, g=(0, 0, 0), material=True)
    assert abs(bm.outflux) < 1e-12


def test_couette_heating_V1_transient_satisfies_heat_equation_and_long_time_limit():  # V1
    U, h, mu, k, rho, cp = 1.0, 1e-3, 1e-3, 0.6, 1000.0, 4182.0
    y = np.linspace(0, h, 11)
    for dpdx in (0.0, -2e3):
        Ts = ch04.couette_heating(y, U, h, mu, k, dpdx=dpdx)["T"]
        T_late = ch04.couette_heating_transient(y, 50.0, U, h, mu, k, rho, cp, dpdx=dpdx)
        assert rel(T_late - 293.15, Ts - 293.15) < 1e-10  # t → ∞: the steady profile
        # interior PDE residual ρC_p∂T/∂t − k∂²T/∂y² − φ(y) by central differences, at t = 0.2 h²/κ
        kap = k / (rho * cp)
        tt, dy, dt = 0.2 * h ** 2 / kap, 1e-5, 1e-3 * h ** 2 / kap
        yi = np.linspace(0.2 * h, 0.8 * h, 5)
        f = lambda yy, t_: ch04.couette_heating_transient(yy, t_, U, h, mu, k, rho, cp, dpdx=dpdx)  # noqa: E731
        Tt = (f(yi, tt + dt) - f(yi, tt - dt)) / (2 * dt)
        Tyy = (f(yi + dy, tt) - 2 * f(yi, tt) + f(yi - dy, tt)) / dy ** 2
        phi = ch04.couette_heating(yi, U, h, mu, k, dpdx=dpdx)["eps"]
        res = rho * cp * Tt - k * Tyy - phi
        assert np.max(np.abs(res)) < 1e-4 * np.max(phi), dpdx


def test_couette_heating_V1_transient_equals_independent_quad_series():  # V1 (loop 2: F1 fix checked on our terms)
    from scipy.integrate import quad as _quad
    U, h, mu, k, rho, cp = 1.0, 1e-3, 1e-3, 0.6, 1000.0, 4182.0
    kap = k / (rho * cp)
    y = np.linspace(0, h, 23)
    for dpdx in (0.0, -2e3, 5e3):
        Ts = lambda yy: float(ch04.couette_heating(yy, U, h, mu, k, dpdx=dpdx)["T"]) - 293.15  # noqa: E731
        bn = []
        for n in range(1, 41):  # adaptive quad, independent of the implementer's recursion
            bn.append(2 / h * _quad(lambda yy: Ts(yy) * np.sin(n * np.pi * yy / h), 0, h, limit=200,
                                    epsabs=1e-16, epsrel=1e-10)[0])
        if dpdx == 0.0:
            assert max(abs(b) for b in bn[1::2]) < 1e-12 * max(abs(b) for b in bn)  # even modes vanish by symmetry
        dT = ch04.couette_heating(y, U, h, mu, k, dpdx=dpdx)["dT_max"]
        for f in (0.002, 0.02, 0.2):  # at these times 40 modes are converged far below 1e-10 ΔT_max
            t = f * h ** 2 / kap
            ours = np.array([Ts(yy) for yy in y]) - sum(b * np.sin((n + 1) * np.pi * y / h) *
                                                        np.exp(-kap * ((n + 1) * np.pi / h) ** 2 * t)
                                                        for n, b in enumerate(bn))
            theirs = ch04.couette_heating_transient(y, t, U, h, mu, k, rho, cp, dpdx=dpdx) - 293.15
            assert np.max(np.abs(theirs - ours)) < 1e-9 * dT, (dpdx, f)
        late = ch04.couette_heating_transient(y, 100 * h ** 2 / kap, U, h, mu, k, rho, cp, dpdx=dpdx) - 293.15
        assert np.max(np.abs(late - np.array([Ts(yy) for yy in y]))) < 1e-12 * dT  # t → ∞: the steady profile
    # truncation error at t = 0 falls like nterms⁻² (b_n ∝ n⁻³)
    errs = [np.max(np.abs(ch04.couette_heating_transient(y, 0.0, U, h, mu, k, rho, cp, nterms=N, dpdx=-2e3) - 293.15))
            for N in (25, 50, 100, 200)]
    assert errs[-1] < errs[0] and errs[-1] < 1e-4 * ch04.couette_heating(y, U, h, mu, k, dpdx=-2e3)["dT_max"]


def test_couette_heating_V1_transient_starts_from_T0_at_documented_accuracy():  # V1 (F1 in loop 1)
    # docstring: "t = 0 is reproduced to ~1/nterms² relative" (default nterms = 200 → ~2.5e-5); we allow 10/nterms².
    U, h, mu, k, rho, cp = 1.0, 1e-3, 1e-3, 0.6, 1000.0, 4182.0
    y = np.linspace(0, h, 11)
    errs = {}
    for dpdx in (0.0, -2e3):
        dT = ch04.couette_heating(y, U, h, mu, k, dpdx=dpdx)["dT_max"]
        T0 = ch04.couette_heating_transient(y, 0.0, U, h, mu, k, rho, cp, dpdx=dpdx)
        errs[dpdx] = float(np.max(np.abs(T0 - 293.15)) / dT)
    assert max(errs.values()) < 10 / 200 ** 2, errs


def test_dissipation_V1_two_routes_signs_and_special_flows():  # V1 + V7 (4.58)
    rng = np.random.default_rng(7)
    Gs = rng.normal(size=(3, 3, 1000))
    sq, con = ch04.dissipation_rate(Gs, 1.3, 0.8, 0.4, form="both")
    assert rel(sq, con) < 1e-12 and np.all(sq >= 0)
    neg = ch04.dissipation_rate(Gs, 1.3, -0.8, 0.0)
    assert np.any(neg < 0)  # μ < 0 would violate the second law
    assert ch04.dissipation_rate([[0, 1000, 0], [0, 0, 0], [0, 0, 0]], 1000.0, 1e-3) == pytest.approx(1.0)  # νγ²
    assert abs(ch04.dissipation_rate(ch04.stress_lab_gradient("rotation", 3.0), 1.0, 1.0)) < 1e-15
    Ginc = ch04.stress_lab_gradient("extension", 2.0)
    S = 0.5 * (Ginc + Ginc.T)
    assert ch04.dissipation_rate(Ginc, 1.0, 0.7) == pytest.approx(2 * 0.7 * np.sum(S * S))  # 2νS:S incompressible
    C = TN.random_rotation(rng)
    G0_ = rng.normal(size=(3, 3))
    assert ch04.dissipation_rate(C @ G0_ @ C.T, 1.0, 0.5, 0.2) == pytest.approx(ch04.dissipation_rate(G0_, 1.0, 0.5, 0.2))
    # a plane (2 × 2) G is the plane flow of a 3-D fluid: same ε as the padded 3 × 3
    G2 = rng.normal(size=(2, 2))
    G3 = np.zeros((3, 3))
    G3[:2, :2] = G2
    assert ch04.dissipation_rate(G2, 1.0, 0.5, 0.2) == pytest.approx(ch04.dissipation_rate(G3, 1.0, 0.5, 0.2), rel=1e-13)
    assert ch04.dissipation_rate(G2, 1.0, 0.5, 0.2, form="contraction") == pytest.approx(
        ch04.dissipation_rate(G3, 1.0, 0.5, 0.2), rel=1e-13)
    assert np.allclose(ch04.viscous_stress(G2, 0.5, 0.2), ch04.viscous_stress(G3, 0.5, 0.2)[:2, :2])
    with pytest.raises(ValueError):
        ch04.dissipation_rate(G2, 1.0, 0.5, form="bad")


def test_dissipation_V2_derivation():  # V2 — D23: σ_ijS_ij = 2μ(S − ⅓S_mmδ)² + μ_vS_mm²
    mu, muv = sp.symbols("mu mu_v", real=True)
    s = sp.symbols("s11 s12 s13 s22 s23 s33", real=True)
    S = sp.Matrix([[s[0], s[1], s[2]], [s[1], s[3], s[4]], [s[2], s[4], s[5]]])
    tr = S.trace()
    sig = 2 * mu * (S - tr / 3 * sp.eye(3)) + muv * tr * sp.eye(3)
    D = S - tr / 3 * sp.eye(3)
    lhs = sum(sig[i, j] * S[i, j] for i in range(3) for j in range(3))
    rhs = 2 * mu * sum(D[i, j] ** 2 for i in range(3) for j in range(3)) + muv * tr ** 2
    assert sp.expand(lhs - rhs) == 0
    assert sp.expand(sum(sp.eye(3)[i, j] ** 2 for i in range(3) for j in range(3)) - 3) == 0  # δ_ijδ_ij = 3
    # σ_ij ∂u_j/∂x_i = σ_ij S_ij for symmetric σ (R07): the rotation part drops out
    w = sp.symbols("w1 w2 w3")
    Rm = sp.Matrix([[0, w[0], w[1]], [-w[0], 0, w[2]], [-w[1], -w[2], 0]])
    assert sp.expand(sum(sig[i, j] * (S + Rm)[j, i] for i in range(3) for j in range(3)) - lhs) == 0


def test_energy_terms_V1_stress_work_kinetic_and_internal_budgets():  # V1 (4.54), (4.56), (4.57)
    coords = (X1, X2, X3)
    uexpr = [sp.sin(X1) * X2, X3 ** 2, sp.cos(X2) * X1]
    u = vec_field(uexpr, coords)
    p = sca_field(1 + X1 * X2 * X3, coords)
    sig = ch04.newtonian_viscous_stress_field(u, mu=0.4, mu_v=0.1, h=1e-4)
    X = rand_pts(6, lo=-0.7, hi=0.7)
    sw = ch04.stress_work_split(p, sig, u, X, 0.0, h=1e-3)
    assert np.max(np.abs(sw["residual"])) < 1e-5 * np.max(np.abs(sw["total"]))
    # Taylor–Green: the mechanical-energy budget closes and the kinetic energy decays by viscous work
    tg, tp, q = ch04.exact_solution_fields("taylor_green", U0=1.0, k=1.0, rho=1.0, mu=0.1)
    Xt = np.random.default_rng(8).uniform(-2, 2, (2, 6))
    sigt = ch04.newtonian_viscous_stress_field(tg, mu=0.1, h=1e-4)
    ke = ch04.kinetic_energy_budget(1.0, tg, tp, sigt, (0, 0), Xt, 0.3, h=1e-3, ht=1e-5)
    assert np.max(np.abs(ke["residual"])) < 1e-5 * np.max(np.abs(ke["lhs"]))
    # Couette: everything 0 in the mechanical budget (steady, u ⟂ ∇u)
    uc, pc = ch04.exact_field("couette", U=1.0, h=0.01)
    kc = ch04.kinetic_energy_budget(1000.0, uc, pc, ch04.newtonian_viscous_stress_field(uc, 1e-3), (0, 0),
                                    np.array([[0.1, 0.3], [0.002, 0.007]]), 0.0, h=1e-5)
    assert max(np.max(np.abs(v)) for v in kc.values()) < 1e-9
    # internal energy (4.57) on Couette with heating: dissipation = νγ², conduction = −dissipation
    U, h, mu, k = 1.0, 1e-3, 1e-3, 0.6
    uf = lambda X_, t: np.stack([U * X_[1] / h, 0 * X_[1]])  # noqa: E731
    qf = lambda X_, t: np.stack([0 * X_[1], -k * (mu / k) * (U / h) ** 2 * (h / 2 - X_[1])])  # noqa: E731
    it = ch04.internal_energy_terms(1000.0, uf, 1e5, ch04.newtonian_viscous_stress_field(uf, mu, h=1e-6), qf,
                                    np.array([[0.0, 0.1], [3e-4, 6e-4]]), 0.0, h=1e-6,
                                    e=lambda X_, t: 4180.0 * ch04.couette_heating(X_[1], U, h, mu, k)["T"])
    assert np.allclose(it["dissipation"], mu * (U / h) ** 2 / 1000.0, rtol=1e-8)
    assert np.allclose(it["conduction"], -it["dissipation"], rtol=1e-6)
    assert np.max(np.abs(it["residual"])) < 1e-6 * np.max(it["dissipation"])


def test_entropy_V1_split_and_production_sign():  # V1 + V2 (4.62)–(4.63)
    coords = (X1, X2, X3)
    Texpr = 300 + 10 * sp.sin(X1) * X2 + X3 ** 2
    T = sca_field(Texpr, coords)
    kk = 0.6
    q = vec_field([-kk * sp.diff(Texpr, v) for v in coords], coords)
    X = rand_pts(8)
    et = ch04.entropy_terms(1.2, T, q, 0.3, X, h=1e-4)
    assert np.max(np.abs(et["split_residual"])) < 1e-8 * np.max(np.abs(et["flux_divergence"]))
    gT = np.array(sp.lambdify(coords, [sp.diff(Texpr, v) for v in coords])(*X), float)
    prod = ch04.entropy_production(gT, T(X), kk, 1.2, 0.3)
    assert rel(et["conduction_production"] + et["dissipation_production"], prod) < 1e-7
    rng = np.random.default_rng(10)
    assert np.all(np.asarray(ch04.entropy_production(rng.normal(size=(3, 100)), rng.uniform(200, 400, 100),
                                                     rng.uniform(0, 1, 100), 1.0, rng.uniform(0, 1, 100))) >= 0)
    assert ch04.entropy_production(0.0, 300.0, 0.6, 1.0, 0.0) == 0.0
    assert ch04.entropy_production(2.0, 300.0, -0.6, 1.0, 0.0) < 0  # k < 0 would violate the second law
    # sympy: the middle and right forms of (4.62) are identical
    x, y, z = sp.symbols("x y z")
    Tf = sp.Function("T")(x, y, z)
    qf = [sp.Function(f"q{i}")(x, y, z) for i in range(3)]
    rho = sp.Symbol("rho", positive=True)
    mid = -sum(sp.diff(qi, v) for qi, v in zip(qf, (x, y, z))) / (rho * Tf)
    right = -sum(sp.diff(qi / Tf, v) for qi, v in zip(qf, (x, y, z))) / rho - sum(
        qi * sp.diff(Tf, v) for qi, v in zip(qf, (x, y, z))) / (rho * Tf ** 2)
    assert sp.simplify(mid - right) == 0


# =====================================================================================================================
# C11 — the Bernoulli function (4.66)–(4.72), (4.78); D24, D25
# =====================================================================================================================
def test_lamb_identity_V2_symbolic_and_sign():  # V2 (4.68)
    x, y, z = sp.symbols("x y z", real=True)
    u = [sp.Function(f"u{i}")(x, y, z) for i in range(3)]
    assert ch04.lamb_identity_sym(u, (x, y, z)) == [0, 0, 0]
    # wrong sign of u × ω fails
    U_ = sp.Matrix([y * z, x ** 2, sp.sin(y)])
    om = sp.Matrix([sp.diff(U_[2], y) - sp.diff(U_[1], z), sp.diff(U_[0], z) - sp.diff(U_[2], x),
                    sp.diff(U_[1], x) - sp.diff(U_[0], y)])
    ke = (U_.T * U_)[0] / 2
    adv = [sum(U_[i] * sp.diff(U_[j], v) for i, v in enumerate((x, y, z))) for j in range(3)]
    lam = U_.cross(om)
    assert any(sp.simplify(adv[j] - (lam[j] + sp.diff(ke, v))) != 0 for j, v in enumerate((x, y, z)))


def test_lamb_identity_V3_stencil_order_2_and_solid_body():  # V3 + V1
    u = vec_field([sp.sin(X2) * X3, X1 ** 2 * X3, sp.cos(X1 + X2)], (X1, X2, X3))
    X = rand_pts(6)
    hs = [0.1, 0.05, 0.025, 0.0125]
    errs = [np.max(np.abs(ch04.lamb_identity_terms(u, X, h=h)["residual"])) for h in hs]
    assert abs(observed_order(hs, errs) - 2.0) < ORDER_TOL, pairwise_orders(hs, errs)
    sb = lambda X_, T: np.stack([-0.5 * X_[1], 0.5 * X_[0]])  # noqa: E731
    Xp = rand_pts(5, d=2)
    L = ch04.lamb_vector(sb, Xp, h=1e-3)
    assert rel(L, np.stack([2 * 0.25 * Xp[0], 2 * 0.25 * Xp[1], 0 * Xp[0]])) < 1e-9  # u × ω = +2Ω²(x, y, 0)
    # consistency with (4.68): (u·∇)u = −u × ω + ∇(½|u|²) = −Ω²(x, y) for solid-body rotation
    lt = ch04.lamb_identity_terms(sb, Xp, h=1e-3)
    assert rel(lt["advective"], -0.25 * Xp) < 1e-9 and rel(lt["minus_u_cross_omega"], -2 * 0.25 * Xp) < 1e-9


def test_pressure_function_V1_closed_forms_vs_quad():  # V1 (4.67)
    po, p = 1.0e5, 1.4e5
    assert ch04.pressure_function(p, po, "constant", rho=998.0) == pytest.approx(0.4e5 / 998.0, rel=1e-14)
    assert ch04.pressure_function(p, po, lambda s: 998.0) == pytest.approx(0.4e5 / 998.0, rel=1e-10)
    R_, T_ = 287.058, 250.0
    iso = ch04.pressure_function(p, po, "isothermal", T=T_, R=R_)
    assert iso == pytest.approx(ch04.pressure_function(p, po, lambda s: s / (R_ * T_)), rel=1e-10)
    rho_o, gam = 1.2, 1.4
    ise = ch04.pressure_function(p, po, "isentropic", rho_o=rho_o, gamma=gam)
    assert ise == pytest.approx(ch04.pressure_function(p, po, lambda s: rho_o * (s / po) ** (1 / gam)), rel=1e-10)
    # isentropic P = h − h_o = C_p(T − T_o) for a perfect gas (link to (4.78))
    To = po / (rho_o * R_)
    T1 = To * (p / po) ** ((gam - 1) / gam)
    assert ise == pytest.approx(gam * R_ / (gam - 1) * (T1 - To), rel=1e-12)
    with pytest.raises(ValueError):
        ch04.pressure_function(p, po, "polytropic")


def test_bernoulli_function_V1_rankine_rotational_and_cylinder_uniform():  # V1 + V4 (4.70)–(4.72)
    Gam, sig, rho = 2 * np.pi, 1.0, 1000.0
    assert ch04.rankine_bernoulli(0.0, Gam, sig)["B"] == pytest.approx(-1.0)
    assert ch04.rankine_bernoulli(2.0, Gam, sig)["B"] == pytest.approx(0.0, abs=1e-15)
    r = np.linspace(0.05, 0.95, 10)
    Bin = ch04.rankine_bernoulli(r, Gam, sig)["B"]
    assert np.allclose(Bin, -Gam ** 2 / (4 * np.pi ** 2 * sig ** 2) + Gam ** 2 * r ** 2 / (4 * np.pi ** 2 * sig ** 4))
    assert np.ptp(Bin) > 0.5  # varies across circles inside the core (rotational)
    assert np.ptp(ch04.rankine_bernoulli(np.linspace(1.0, 5, 10), Gam, sig)["B"]) < 1e-14  # uniform outside
    # dp/dr = ρu_θ²/r and continuity at σ
    rr = np.array([0.3, 0.9, 1.5, 3.0])
    dp = (ch04.rankine_vortex_pressure(rr + 1e-6, Gam, sig, rho) - ch04.rankine_vortex_pressure(rr - 1e-6, Gam, sig, rho)) / 2e-6
    ut = ch04.rankine_bernoulli(rr, Gam, sig)["u_theta"]
    assert rel(dp, rho * ut ** 2 / rr) < 1e-7
    assert ch04.rankine_vortex_pressure(1 - 1e-12, Gam, sig) == pytest.approx(ch04.rankine_vortex_pressure(1.0, Gam, sig))
    # ∇B = u × ω everywhere for the steady Rankine vortex (Lamb surfaces)
    def u2(X, T):
        r_ = np.hypot(X[0], X[1])
        ut_ = ch04.rankine_bernoulli(r_, Gam, sig)["u_theta"]
        return np.stack([-ut_ * X[1] / r_, ut_ * X[0] / r_])
    Bf = lambda X, T: ch04.rankine_bernoulli(np.hypot(X[0], X[1]), Gam, sig)["B"]  # noqa: E731
    Xp = np.stack([np.array([0.3, -0.5, 1.8, 2.5]), np.array([0.2, 0.4, -1.0, 0.7])])
    assert np.max(ch04.lamb_surface_check(u2, Bf, Xp, h=1e-5)) < 1e-6
    lt = ch04.lamb_surface_terms(u2, Bf, Xp, h=1e-5)
    assert np.max(np.abs(lt["u_dot_gradB"])) < 1e-6
    # cylinder potential flow: B uniform everywhere (4.72) at 1000 points
    uc, pc = ch04.exact_field("cylinder", U=2.0, a=1.0, rho=rho, p_inf=1e5)
    rngp = np.random.default_rng(13)
    rp, tp = rngp.uniform(1.0, 4.0, 1000), rngp.uniform(0, 2 * np.pi, 1000)
    Xc = np.stack([rp * np.cos(tp), rp * np.sin(tp)])
    Bc = ch04.bernoulli_function(np.hypot(*uc(Xc, 0)), pc(Xc, 0), 0.0, rho=rho)
    assert np.ptp(Bc) / abs(np.mean(Bc)) < 1e-12


def test_stagnation_V1_enthalpy_temperature_and_isentropic_relation():  # V1 (4.78)
    from fluidpy.core.thermo import CP_AIR, GAMMA_AIR, R_AIR
    assert ch04.stagnation_temperature(300.0, 100.0) == pytest.approx(304.98, abs=5e-3)
    for T, U in [(288.15, 50.0), (220.0, 300.0)]:
        T0 = ch04.stagnation_temperature(T, U, cp=CP_AIR)
        M = U / np.sqrt(GAMMA_AIR * R_AIR * T)
        assert T0 / T == pytest.approx(1 + (GAMMA_AIR - 1) / 2 * M ** 2, rel=1e-13)
    assert ch04.stagnation_enthalpy(3e5, 20.0, 5.0) == pytest.approx(3e5 + 200 + 5 * G)
    s = ch04.bernoulli_scenario("hot_nozzle")
    assert np.ptp(s["B_along"]) < 1e-9 and "4.78" in s["valid"]


def test_which_bernoulli_V1_hypotheses_table():  # V1 (N102)
    assert ch04.which_bernoulli(True, False, True, True, constant_density=True) == ["4.19", "4.71", "4.72", "4.75"]
    assert ch04.which_bernoulli(True, False, True, True) == ["4.71", "4.72", "4.75"]
    assert ch04.which_bernoulli(False, False, False, True) == []  # inviscid, rotational, unsteady → none
    assert ch04.which_bernoulli(True, False, False, True, isentropic=True) == ["4.71", "4.78"]
    assert ch04.which_bernoulli(False, True, True, constant_density=True) == ["4.82"]
    assert ch04.which_bernoulli(True, True, False, True) == []
    assert ch04.which_bernoulli_text(True, False, True, True) == "4.71,4.72,4.75"
    assert ch04.which_bernoulli_text(False, False, False) == ""
    assert set(ch04.BERNOULLI_FORMS) == {"4.19", "4.71", "4.72", "4.75", "4.78", "4.82"}


def test_bernoulli_scenarios_V1_flatness_matches_valid_forms():  # V1 (E6)
    for name in ch04.BERNOULLI_SCENARIOS:
        s = ch04.bernoulli_scenario(name)
        assert s["name"] == name and s["valid"] == ch04.which_bernoulli(
            s["hypotheses"].get("steady", False), s["hypotheses"].get("viscous", False),
            s["hypotheses"].get("irrotational", False), barotropic=s["hypotheses"].get("constant_density", False),
            isentropic=s["hypotheses"].get("isentropic", False),
            constant_density=s["hypotheses"].get("constant_density", False))
    for name in ("pitot", "orifice", "cylinder"):
        s = ch04.bernoulli_scenario(name)
        assert np.ptp(s["B_along"]) < 1e-9 * max(1, np.max(np.abs(s["B_along"])))
    rk = ch04.bernoulli_scenario("rankine")
    assert np.ptp(rk["B_along"]) < 1e-12 and rk["numbers"]["B_inside_span"] > 1e-3 and rk["numbers"]["B_outside_span"] < 1e-12
    assert rk["status"].startswith("rotational")
    ut = ch04.bernoulli_scenario("u_tube")
    assert abs(ut["numbers"]["residual_4_82"]) < 1e-12 and np.ptp(ut["B_along"]) > 0
    pt = ch04.bernoulli_scenario("pitot", U=3.0)
    assert pt["numbers"]["U_measured"] == pytest.approx(3.0)
    cy = ch04.bernoulli_scenario("cylinder")
    assert rel(cy["numbers"]["Cp_surface"], 1 - 4 * np.sin(cy["numbers"]["theta_surface"]) ** 2) < 1e-15
    with pytest.raises(ValueError):
        ch04.bernoulli_scenario("nope")


def test_bernoulli_function_V2_derivation():  # V2 — D24 (4.66) + (4.67) + (4.68) → (4.69); D25
    x, y, z, t = sp.symbols("x y z t", real=True)
    coords = (x, y, z)
    u = [sp.Function(f"u{i}")(x, y, z, t) for i in range(3)]
    p = sp.Function("p")(x, y, z, t)
    rho_of = sp.Function("rho")
    s = sp.Symbol("s")
    po = sp.Symbol("p_o")
    P = sp.Integral(1 / rho_of(s), (s, po, p))  # (4.67) pressure function
    for v in coords:
        assert sp.simplify(sp.diff(P, v) - sp.diff(p, v) / rho_of(p)) == 0  # FTC + chain rule
    Phi = sp.Function("Phi")(x, y, z)
    U = sp.Matrix(u)
    om = sp.Matrix([sp.diff(u[2], y) - sp.diff(u[1], z), sp.diff(u[0], z) - sp.diff(u[2], x), sp.diff(u[1], x) - sp.diff(u[0], y)])
    B = sum(ui ** 2 for ui in u) / 2 + P + Phi
    lam = U.cross(om)
    for j, v in enumerate(coords):
        euler = sp.diff(u[j], t) + sum(u[i] * sp.diff(u[j], coords[i]) for i in range(3)) + sp.diff(p, v) / rho_of(p) + sp.diff(Phi, v)
        form69 = sp.diff(u[j], t) + sp.diff(B, v) - lam[j]
        assert sp.simplify(sp.expand(euler - form69)) == 0  # (4.66) ≡ (4.69)
    assert sp.simplify(U.dot(lam)) == 0 and sp.simplify(om.dot(lam)) == 0  # D25: u·∇B = ω·∇B = 0


# =====================================================================================================================
# C12 — unsteady Bernoulli (4.73)–(4.83); D26
# =====================================================================================================================
def test_unsteady_bernoulli_V1_accelerating_sphere_bracket_uniform():  # V1 (4.74)/(4.75)
    a, U0, dUdt, rho = 0.1, 1.0, 0.5, 1000.0
    phi, p = ch04.accelerating_sphere_fields(a, U0, dUdt, rho, p_inf=1e5)
    rng = np.random.default_rng(14)
    r, th, ph = rng.uniform(1.2 * a, 4 * a, 30), rng.uniform(0, np.pi, 30), rng.uniform(0, 2 * np.pi, 30)
    t = 0.4
    Xc = U0 * t + 0.5 * dUdt * t ** 2
    X = np.stack([Xc + r * np.cos(th), r * np.sin(th) * np.cos(ph), r * np.sin(th) * np.sin(ph)])
    B = ch04.unsteady_bernoulli_B(phi, p, X, t, rho=rho, h=1e-6, ht=1e-6)
    assert np.ptp(B) < 1e-6 * np.max(np.abs(B))  # the same B(t) everywhere (here p∞/ρ)
    # and the field satisfies Euler's equation (the independent route): u = ∇φ, p from (4.75)
    u = lambda X_, T: np.stack([(phi(X_ + e[:, None] * 1e-6, T) - phi(X_ - e[:, None] * 1e-6, T)) / 2e-6  # noqa: E731
                                for e in np.eye(3)])
    terms = ch04.ns_incompressible_terms(u, p, X, t, rho, 0.0, (0, 0, 0), h=1e-4, ht=1e-5)
    scale = np.max(np.abs(terms.local)) + np.max(np.abs(terms.advective))
    assert np.max(np.abs(terms.residual)) < 1e-4 * scale


@needs_ref
def test_accelerating_sphere_V1_added_mass_form_cross_check():  # V1 (form cross-check, not V5)
    frac = ref_json()["added_mass_sphere"]["fraction_of_displaced"]
    for (a, dUdt, rho, U) in [(0.1, 1.0, 1000.0, 0.0), (0.05, -3.0, 1.2, 2.0)]:
        d = ch04.accelerating_sphere_pressure(0.3, a, dUdt, rho, U=U)
        assert d["added_mass"] == pytest.approx(frac * rho * 4 / 3 * np.pi * a ** 3, rel=1e-14)
        assert d["force"] == pytest.approx(-d["added_mass"] * dUdt, rel=1e-10)  # U² part integrates to 0 (d'Alembert)
    # surface pressure formula vs the field's p at r = a (fixed point on the sphere at t)
    a, U0, dUdt, rho, t = 0.1, 1.0, 0.5, 1000.0, 0.4
    phi, p = ch04.accelerating_sphere_fields(a, U0, dUdt, rho)
    for th in (0.2, 1.1, 2.5):
        Xc = U0 * t + 0.5 * dUdt * t ** 2
        X = np.array([[Xc + a * np.cos(th)], [a * np.sin(th)], [0.0]])
        pf = float(p(X, t)[0])
        pd = ch04.accelerating_sphere_pressure(th, a, dUdt, rho, U=U0 + dUdt * t)["p"]
        assert pf == pytest.approx(pd, rel=1e-12, abs=1e-9)


def test_unsteady_bernoulli_V1_gauge_utube_and_viscous_irrotational():  # V1 (4.74)–(4.82)
    assert ch04.gauge_absorbed_bracket(3.0) == 0.0 and ch04.gauge_absorbed_bracket(3.0, sign=+1.0) == 6.0
    assert ch04.unsteady_bernoulli_pressure(2.0, 3.0, 1.0, rho=1000.0, C=20.0) == pytest.approx(
        1000 * (20 - 2 - 4.5 - G))
    assert ch04.unsteady_bernoulli_pressure(0.0, np.array([[3.0], [4.0]]), 0.0, rho=1.0) == pytest.approx(-12.5)
    for t in (0.0, 0.3, 1.7):
        u = ch04.u_tube_column(t, L=1.2, h0=0.05)
        assert abs(u["residual_4_82"]) < 1e-12 and u["dp"] == pytest.approx(1000 * 1.2 * u["dUdt"])
        assert u["dUdt"] == pytest.approx(-2 * G * u["h"] / 1.2, rel=1e-12)  # L dU/dt + 2gh = 0
    # straight pipe of length L: p₁ − p₂ = ρL dU/dt from (4.82)
    s = np.linspace(0, 3.0, 31)
    r = ch04.unsteady_streamline_bernoulli(np.full(31, 0.7), s, dict(U=2, z=0, p=1e5), dict(U=2, z=0, p=1e5 - 1000 * 3 * 0.7))
    assert abs(r) < 1e-10
    # viscous force of potential flows vanishes (4.80); a curved shear flow's does not
    for name, kw in [("cylinder", dict(U=1.0, a=1.0)), ("source_stream", dict(U=1.0, m=1.0)), ("vortex", dict(Gamma=2.0))]:
        uf = lambda X, T, n=name, k=kw: np.stack(ch04.velocity_preset(n, X[0], X[1], **k))  # noqa: E731
        X = np.stack([np.array([1.5, -2.0, 0.4]), np.array([0.7, 1.1, -1.9])])
        assert np.max(np.abs(ch04.viscous_irrotational_residual(uf, X, mu=1.0, h=1e-3))) < 1e-5, name
    par = lambda X, T: np.stack([X[1] ** 2, 0 * X[1]])  # noqa: E731
    assert np.max(np.abs(ch04.viscous_irrotational_residual(par, np.array([[0.0], [0.5]]), mu=1.0, h=1e-3))) > 1.0


def test_unsteady_bernoulli_V2_derivation():  # V2 — D26: (4.73) → (4.74) → (4.75), and the gauge sign
    x, y, z, t = sp.symbols("x y z t", real=True)
    phi = sp.Function("phi")(x, y, z, t)
    u = [sp.diff(phi, v) for v in (x, y, z)]
    rho = sp.Symbol("rho", positive=True)
    p, Phi = sp.Function("p")(x, y, z, t), sp.Function("Phi")(x, y, z)
    bracket = sp.diff(phi, t) + sum(ui ** 2 for ui in u) / 2 + p / rho + Phi
    for j, v in enumerate((x, y, z)):
        euler = sp.diff(u[j], t) + sum(u[i] * sp.diff(u[j], w) for i, w in enumerate((x, y, z))) + sp.diff(p, v) / rho + sp.diff(Phi, v)
        assert sp.simplify(sp.expand(euler - sp.diff(bracket, v))) == 0  # Euler ≡ ∇[bracket] (Schwarz)
    Bt = sp.Function("B")(t)
    tp = sp.Symbol("tp")
    minus = sp.diff(-sp.Integral(Bt.subs(t, tp), (tp, 0, t)), t)  # φ_new = φ − ∫B dt′ adds −B to ∂φ/∂t
    plus = sp.diff(sp.Integral(Bt.subs(t, tp), (tp, 0, t)), t)
    assert sp.simplify(Bt + minus) == 0 and sp.simplify(Bt + plus - 2 * Bt) == 0


# =====================================================================================================================
# C13 — Boussinesq (4.84)–(4.89); D27, D28
# =====================================================================================================================
def test_boussinesq_V1_base_state_perturbations_and_buoyancy():  # V1 (4.84), N², ch01 link
    z = np.linspace(0, 100, 201)
    rho_s = lambda zz: 1025.0 - 0.02 * np.asarray(zz)  # noqa: E731
    ps = STAT.integrate_hydrostatic(z, lambda zz, pp: rho_s(zz), 2e5, g=G, z0=0.0)
    pp, rp = ch04.perturbation_fields(ps, rho_s(z), z, rho_s, g=G)
    assert np.max(np.abs(pp)) < 1e-8 and np.max(np.abs(rp)) < 1e-12
    pp2, rp2 = ch04.perturbation_fields(1e5 - 1000.0 * G * z + 50.0, np.full(z.shape, 1000.3), z, lambda zz: 1000.0 + 0 * zz,
                                        g=G, p_s0=1e5)
    assert np.allclose(pp2, 50.0) and np.allclose(rp2, 0.3)
    assert ch04.buoyancy(-0.5, 1000.0) == pytest.approx(G * 0.5e-3)  # lighter fluid: upward
    # N² = ∂b_s/∂z = −(g/ρ₀)dρ_s/dz equals ch01's incompressible N² (the adiabatic gradient is 0 for a liquid)
    N2_b = (ch04.buoyancy(rho_s(10.0 + 1e-3) - 1025.0, 1025.0) - ch04.buoyancy(rho_s(10.0 - 1e-3) - 1025.0, 1025.0)) / 2e-3
    assert N2_b == pytest.approx(STRAT.brunt_vaisala_sq(1025.0, -0.02, 0.0, g=G), rel=1e-9)
    assert ch04.boussinesq_density(303.15, 1000.0, 2e-4, 293.15) == pytest.approx(998.0)
    sw = ch01.seawater_density_linear(283.15 + 5.0, 35.0)
    assert ch04.boussinesq_density(288.15, 1027.0, 1.67e-4, 283.15) == pytest.approx(sw, rel=1e-14)


def test_boussinesq_V1_momentum_terms_rest_and_hydrostatic_perturbation():  # V1 (4.86)
    X = rand_pts(6)
    rest = lambda X_, T: 0.0 * X_  # noqa: E731
    d = ch04.boussinesq_momentum_terms(rest, 0.0, 0.0, X, 0.0)
    assert all(np.max(np.abs(v)) == 0 for v in d.values())
    # a static density anomaly ρ′(z) balanced by p′ = −g∫ρ′dz: no motion, residual 0
    rp = lambda X_, T: 0.5 * np.sin(X_[2])  # noqa: E731
    pp = lambda X_, T: -G * 0.5 * (1 - np.cos(X_[2]))  # noqa: E731
    d2 = ch04.boussinesq_momentum_terms(rest, pp, rp, X, 0.0, rho0=1000.0, g=G, h=1e-4)
    assert np.max(np.abs(d2["residual"])) < 1e-9 and np.max(np.abs(d2["buoyancy"])) > 1e-5
    assert np.max(np.abs(d2["dropped"])) == 0
    # an unbalanced anomaly accelerates the fluid (sign: ρ′ > 0 sinks)
    d3 = ch04.boussinesq_momentum_terms(rest, 0.0, 0.2, X, 0.0, rho0=1000.0)
    assert np.all(d3["residual"][2] > 0)  # local acc (0) − buoyancy (−gρ′/ρ₀) > 0 ⇒ the equation wants Du/Dt < 0


def test_boussinesq_heat_V1_blob_residual_conservation_and_order():  # V1 + V4 + V3 (4.89)
    kap, U, s0 = 1e-3, 0.2, 0.05
    T, u = ch04.gaussian_blob_fields(U=U, kappa=kap, sigma0=s0, dim=2)
    X = np.random.default_rng(15).uniform(-0.1, 0.1, (2, 12))
    t = 2.0
    r = ch04.temperature_equation_residual(T, u, X + np.array([[U * t], [0]]), t, kap, h=2e-5, ht=1e-5)
    scale = np.max(np.abs(ch04.heat_equation_terms(T, u, X + np.array([[U * t], [0]]), t, 1000.0, 4182.0, kap * 4.182e6,
                                                   h=2e-5, ht=1e-5)["conduction"]))
    assert np.max(np.abs(r)) < 1e-6 * scale
    ht_ = ch04.heat_equation_terms(T, u, X, t, 1000.0, 4182.0, kap * 4.182e6, eps=0.0, h=2e-5, ht=1e-5)
    assert np.max(np.abs(ht_["residual"])) < 1e-6 * max(scale, 1e-12) and np.all(ht_["dissipation"] == 0)
    hs = [0.02, 0.01, 0.005, 0.0025]
    errs = [np.max(np.abs(ch04.temperature_equation_residual(T, u, X, t, kap, h=h, ht=1e-5))) for h in hs]
    assert abs(observed_order(hs, errs) - 2.0) < ORDER_TOL, pairwise_orders(hs, errs)
    # ∫(T − T₀)dA conserved = A(2πσ₀²), the centre moves at U
    for tt in (0.0, 3.0, 10.0):
        from scipy.integrate import dblquad
        tot = dblquad(lambda yy, xx: float(ch04.gaussian_blob_advection_diffusion(np.array([xx, yy]), tt, U, kap, s0)),
                      -1 + U * tt, 1 + U * tt, -1, 1, epsabs=1e-12)[0]
        assert tot == pytest.approx(2 * np.pi * s0 ** 2, rel=1e-8)
    xs = np.linspace(-0.5, 2.5, 30001)
    prof = ch04.gaussian_blob_advection_diffusion(xs, 5.0, U, kap, s0, dim=1)
    assert xs[np.argmax(prof)] == pytest.approx(U * 5.0, abs=2e-4)


def test_boussinesq_validity_V1_numbers_and_scenarios():  # V1 (N108, N112)
    b = ch04.boussinesq_validity(2e-4, 10.0, 10.0, 0.1, c=1500.0, nu=1e-6, cp=4186.0)
    assert b["alpha_dT"] == pytest.approx(2.0e-3) and b["g_prime"] == pytest.approx(0.01962) and b["valid"]
    assert b["H_c"] == pytest.approx(1500.0 ** 2 / G) and b["heating_ratio"] == pytest.approx(1e-7 / (4186.0 * 100))
    c_air = 100.0 / ch04.mach_number(100.0, T=288.15)
    assert c_air ** 2 / G / 1e3 == pytest.approx(11.8, abs=0.05)  # c²/g for air at 288 K ≈ 11.8 km
    for name in ("lake", "thermocline", "lab_tank", "abl"):
        assert ch04.boussinesq_validity(**{k: v for k, v in ch04.boussinesq_scenario(name).items()
                                           if k in ("alpha", "dT", "L", "U", "c", "nu", "cp")})["valid"], name
    deep = ch04.boussinesq_validity(**{k: v for k, v in ch04.boussinesq_scenario("deep_atmosphere").items()
                                       if k in ("alpha", "dT", "L", "U", "c", "nu", "cp")})
    assert not deep["valid"] and "deep layer" in deep["verdict"]
    with pytest.raises(ValueError):
        ch04.boussinesq_scenario("mars")
    br = ch04.blob_rise(np.array([0.0, 5.0, 1e4]), 0.02, 4.0)
    w, z = br["w"], br["z"]
    assert w[0] == 0 and z[0] == 0 and w[-1] == pytest.approx(0.02 * 4.0)  # terminal speed g′τ_d
    tt = np.array([1.0, 7.0])
    dw = (ch04.blob_rise(tt + 1e-5, 0.02, 4.0)["w"] - ch04.blob_rise(tt - 1e-5, 0.02, 4.0)["w"]) / 2e-5
    assert rel(dw, 0.02 - ch04.blob_rise(tt, 0.02, 4.0)["w"] / 4.0) < 1e-8  # dw/dt = g′ − w/τ_d
    dz = (ch04.blob_rise(tt + 1e-5, 0.02, 4.0)["z"] - ch04.blob_rise(tt - 1e-5, 0.02, 4.0)["z"]) / 2e-5
    assert rel(dz, ch04.blob_rise(tt, 0.02, 4.0)["w"]) < 1e-8


def test_boussinesq_V2_derivation():  # V2 — D27 (subtract the base state) and D28 (C_v → C_p)
    x, y, z, t, g = sp.symbols("x y z t g", real=True)
    ps = sp.Function("p_s")(z)
    rhos = -sp.diff(ps, z) / g  # hydrostatic base: dp_s/dz = −ρ_s g
    pp, rp = sp.Function("pp")(x, y, z, t), sp.Function("rp")(x, y, z, t)
    p, rho = ps + pp, rhos + rp
    full = [-sp.diff(p, v) for v in (x, y, z)]
    full[2] += -rho * g
    pert = [-sp.diff(pp, v) for v in (x, y, z)]
    pert[2] += -rp * g
    assert all(sp.simplify(a_ - b_) == 0 for a_, b_ in zip(full, pert))  # −∇p + ρg = −∇p′ + ρ′g (4.84)
    # D28: perfect gas at (nearly) constant p: ρ = p/(RT) ⇒ (1/ρ)Dρ/Dt = −(1/T)DT/Dt; p∇·u = ρR DT/Dt
    R, cv, P = sp.symbols("R c_v P", positive=True)
    Tt = sp.Function("T")(t)
    rho_t = P / (R * Tt)
    divu = -sp.diff(rho_t, t) / rho_t  # continuity (4.8)
    lhs = rho_t * cv * sp.diff(Tt, t) + P * divu  # ρ De/Dt + p∇·u with e = C_vT
    assert sp.simplify(lhs - rho_t * (cv + R) * sp.diff(Tt, t)) == 0  # = ρC_p DT/Dt with C_p = C_v + R (4.88)
    alpha = sp.simplify(-sp.diff(P / (R * sp.Symbol("T", positive=True)), sp.Symbol("T", positive=True)) /
                        (P / (R * sp.Symbol("T", positive=True))))
    assert sp.simplify(alpha - 1 / sp.Symbol("T", positive=True)) == 0  # α = 1/T


# =====================================================================================================================
# C14 — kinematic boundary condition (4.90)–(4.93); interface conditions; surface tension; Example 4.7; D29
# =====================================================================================================================
def test_kinematic_bc_V1_moving_wall_and_moving_front():  # V1
    X = np.stack([np.full(5, 0.6 * 1.3), np.linspace(-1, 1, 5)])
    assert np.allclose(ch04.kinematic_bc_residual("moving_wall", None, X, 0.6, V=1.3), 0, atol=1e-10)
    assert np.allclose(ch04.surface_normal_speed("moving_wall", X, 0.6, V=1.3), 1.3)
    assert np.allclose(ch04.relative_normal_velocity("moving_wall", None, X, 0.6, V=1.3), 0, atol=1e-10)
    # a plane front η = x − Vt in a uniform flow U: relative normal velocity U − V, mass flux ρ(U − V)
    eta = lambda X_, t: X_[0] - 0.4 * t  # noqa: E731
    uf = lambda X_, t: np.stack([2.0 + 0 * X_[0], 0.5 + 0 * X_[0]])  # noqa: E731
    assert np.allclose(ch04.relative_normal_velocity(eta, uf, X, 1.0), 1.6)
    assert np.allclose(ch04.interface_mass_flux(eta, uf, 1.2, X, 1.0), 1.2 * 1.6)
    assert np.allclose(ch04.interface_mass_flux(eta, uf, lambda X_, t: 2.0 + 0 * X_[0], X, 1.0), 3.2)
    # an oblique front: only the normal component counts; scaling η does not change (4.92)
    eta2 = lambda X_, t: 3.0 * (X_[0] + X_[1] - 0.2 * t)  # noqa: E731
    n = np.array([1, 1]) / np.sqrt(2)
    assert np.allclose(ch04.relative_normal_velocity(eta2, uf, X, 0.5), 2.5 / np.sqrt(2) - 0.2 / np.sqrt(2))
    assert np.allclose(ch04.surface_normal_speed(eta2, X, 0.5), 0.2 / np.sqrt(2)) and n @ n == pytest.approx(1)
    with pytest.raises(ValueError):
        ch04.surface_preset("drum")


def test_kinematic_bc_V3_linear_wave_residual_is_second_order_in_amplitude():  # V3
    k = 1.2
    x = np.linspace(0, 2 * np.pi / k, 17)
    amps = np.array([0.08, 0.04, 0.02, 0.01])
    errs = [np.max(np.abs(ch04.wave_kinematic_residual(x, 0.3, a, k))) for a in amps]
    assert abs(observed_order(amps, errs) - 2.0) < ORDER_TOL, pairwise_orders(amps, errs)
    om = np.sqrt(G * k)
    assert errs[-1] / (0.01 * om) == pytest.approx(k * 0.01, rel=0.1)  # residual/(aω) ~ ka
    assert np.max(np.abs(ch04.wave_kinematic_residual(x, 0.3, 0.05, k, full=False))) < 1e-15  # linearised: exact
    assert np.max(np.abs(ch04.wave_kinematic_residual(x, 0.3, 0.05, k, H=2.0, full=False))) < 1e-15
    w = ch04.linear_wave_surface(x, 0.0, 0.3, 0.05, k, H=2.0)
    assert w["omega"] == pytest.approx(np.sqrt(G * k * np.tanh(2 * k)))
    eta, u = ch04.linear_wave_fields(0.05, k, H=2.0)
    Xs = np.stack([x, 0.05 * np.cos(k * x - w["omega"] * 0.3)])
    assert np.allclose(eta(Xs, 0.3), 0.0, atol=1e-15)


def test_kinematic_bc_V2_derivation():  # V2 — D29
    t, x, z, a, k, om = sp.symbols("t x z a k omega", positive=True)
    X_, Z_ = sp.Function("X")(t), sp.Function("Z")(t)
    eta = sp.Function("eta")
    ddt_path = sp.diff(eta(X_, Z_, t), t)  # a particle riding on the surface: dη/dt along its path
    u, w = sp.symbols("u w")
    chain = ddt_path.subs({sp.Derivative(X_, t): u, sp.Derivative(Z_, t): w})
    expected = (sp.Subs(sp.Derivative(eta(x, Z_, t), x), x, X_) * u).doit()
    assert sp.simplify(chain.doit() - (expected + (sp.Subs(sp.Derivative(eta(X_, z, t), z), z, Z_) * w).doit()
                                       + sp.Subs(sp.Derivative(eta(X_, Z_, sp.Symbol("s")), sp.Symbol("s")),
                                                 sp.Symbol("s"), t).doit())) == 0
    # linear deep-water wave: the full (4.91) at z = η starts at O(a²); the linearised form is exact
    ph = k * x - om * t
    etaexpr = z - a * sp.cos(ph)
    uu, ww = a * om * sp.exp(k * z) * sp.cos(ph), a * om * sp.exp(k * z) * sp.sin(ph)
    res = (sp.diff(etaexpr, t) + uu * sp.diff(etaexpr, x) + ww * sp.diff(etaexpr, z)).subs(z, a * sp.cos(ph))
    ser = sp.series(res, a, 0, 3).removeO()
    assert sp.simplify(ser.coeff(a, 1)) == 0 and sp.simplify(ser.coeff(a, 2)) != 0
    lin = sp.diff(a * sp.cos(ph), t) - ww.subs(z, 0)
    assert sp.simplify(lin) == 0


def test_interfaces_V1_pillbox_two_layers_and_slip():  # V1 + V3 (N116, N117)
    ls = np.array([0.1, 0.05, 0.025, 0.0125])
    d = ch04.pillbox_limit(3.0, 3.0, 2.0, -0.5, ls)
    assert abs(observed_order(ls, np.abs(d["residual"])) - 1.0) < 1e-12  # side + volume ∝ l
    assert d["jump"] == 0.0
    tl = ch04.two_layer_conduction(np.array([0.0, 0.1, 0.1 - 1e-12, 0.1 + 1e-12, 0.3]), 1.0, 0.25, 0.1, 0.2, 400.0, 300.0)
    assert tl["q"] == pytest.approx(100.0 / (0.1 + 0.8)) and tl["T"][2] == pytest.approx(tl["T"][3], abs=1e-8)
    assert tl["T"][0] == 400.0 and tl["T"][-1] == pytest.approx(300.0) and tl["T_interface"] == pytest.approx(tl["T"][1])
    s1, s2 = -tl["q"] / 1.0, -tl["q"] / 0.25
    assert s2 / s1 == pytest.approx(1.0 / 0.25)  # slopes in ratio k₁/k₂ … flux k∂T/∂n continuous
    tf = ch04.two_fluid_couette(np.array([0.0, 0.3, 0.3 + 1e-12, 1.0]), 1e-3, 5e-3, 0.3, 0.7, 2.0)
    assert tf["u"][0] == 0.0 and tf["u"][-1] == pytest.approx(2.0) and tf["u"][1] == pytest.approx(tf["u"][2], abs=1e-9)
    assert 1e-3 * tf["u_interface"] / 0.3 == pytest.approx(tf["tau"]) and 5e-3 * (2.0 - tf["u_interface"]) / 0.7 == pytest.approx(tf["tau"])
    y = np.linspace(0, 1e-3, 5)
    assert np.allclose(ch04.navier_slip_couette(y, 1.0, 1e-3, 0.0), y / 1e-3)
    us = ch04.navier_slip_couette(y, 1.0, 1e-3, 2e-4)
    assert us[0] == pytest.approx(2e-4 * (us[1] - us[0]) / (y[1] - y[0]))  # u = b du/dy at the wall


def test_surface_tension_V1_cap_forces_and_laplace_jump():  # V1 + V3 (4.97)–(4.98), N127
    for (R1, R2, zeta, dp) in [(1e-3, 2e-3, 1e-5, 50.0), (0.5, 0.2, 1e-3, 3.0)]:
        Fz = ch04.cap_pressure_force(dp, R1, R2, zeta)
        assert Fz == pytest.approx(ch04.cap_pressure_force(dp, R1, R2, zeta, exact=False), rel=1e-10)
        comps = ch04.cap_pressure_force(dp, R1, R2, zeta, components=True)
        assert np.max(np.abs(comps[:2])) < 1e-10 * abs(Fz)
        st = ch04.cap_surface_tension_force(0.07, R1, R2, zeta, components=True)
        assert np.max(np.abs(st[:2])) < 1e-12 * abs(st[2])
    R1, R2, sig = 1e-3, 3e-3, 0.0728
    zetas = np.array([1e-5, 5e-6, 2.5e-6, 1.25e-6])
    ratio = np.array([ch04.cap_surface_tension_force(sig, R1, R2, zz) / ch04.cap_surface_tension_force(sig, R1, R2, zz, exact=False)
                      for zz in zetas])
    assert abs(observed_order(zetas, np.abs(ratio - 1)) - 1.0) < ORDER_TOL  # (ratio − 1) = O(ζ/R)
    assert ch04.laplace_jump_from_balance(0.0728, 1e-3, 1e-3) == pytest.approx(145.6, rel=1e-12)
    assert ch04.laplace_jump_from_balance(sig, R1, R2) == pytest.approx(ch01.laplace_pressure_jump(sig, R1, R2), rel=1e-12)
    assert ch04.laplace_jump_from_balance(sig, R1, R2, zeta=1e-9, exact=True) == pytest.approx(
        ch01.laplace_pressure_jump(sig, R1, R2), rel=1e-5)


def test_surface_tension_V2_quarter_path_integral():  # V2 — the small-ζ evaluation of (4.98) → Laplace (1.5)
    R1, R2, zeta, xi = sp.symbols("R1 R2 zeta xi", positive=True)
    a, b = sp.sqrt(2 * R1 * zeta), sp.sqrt(2 * R2 * zeta)
    # small cap: n ≈ (−x/R₁, −y/R₂, 1), t horizontal ⇒ (t × n)_z ds = x dy/R₁ − y dx/R₂ on the rim ellipse
    x, y = a * sp.cos(xi), b * sp.sin(xi)
    integrand = x * sp.diff(y, xi) / R1 - y * sp.diff(x, xi) / R2
    Fz = 4 * sp.integrate(sp.expand(integrand), (xi, 0, sp.pi / 2))  # four equal quarter paths
    target = sp.pi * a * b * (1 / R1 + 1 / R2)
    assert sp.simplify(Fz - target) == 0
    dp, sig = sp.symbols("dp sigma", positive=True)
    sol = sp.solve(sp.Eq(-dp * sp.pi * a * b + sig * target, 0), dp)[0]  # F_p + F_st = 0
    assert sp.simplify(sol - sig * (1 / R1 + 1 / R2)) == 0  # Laplace (1.5)
    # the numeric closed form of the code is this expression
    assert ch04.cap_surface_tension_force(0.07, 1e-3, 2e-3, 1e-6, exact=False) == pytest.approx(
        float(target.subs({R1: 1e-3, R2: 2e-3, zeta: 1e-6})) * 0.07, rel=1e-14)


@needs_ref
def test_capillary_length_V5_water_20C():  # V5
    ref = ref_json()["capillary_length_water_20C"]
    rows = [ln.split(",") for ln in (REF1 / "iapws_sigma.csv").read_text().splitlines()[1:]]
    sig20 = [float(r[3]) for r in rows if float(r[0]) == 20.0][0] * 1e-3
    lc = ch04.capillary_length(sig20, 998.2, g=9.80665, rho_other=1.2)
    half = 0.5 * 10 ** (np.floor(np.log10(ref["length_m"])) - (ref["digits"] - 1))
    assert abs(lc - ref["length_m"]) <= half + 0.01 * ref["length_m"]  # 1 % plus the source's 3-digit rounding
    assert ch04.capillary_length(0.0728, 998.0) == pytest.approx(2.72688e-3, abs=5e-9)  # design prints 2.7266 mm (slip)
    assert ch04.bond_number(998.0, G, 0.01, 0.0728) == pytest.approx((0.01 / ch04.capillary_length(0.0728, 998.0)) ** 2)


def test_meniscus_V1_closed_form_vs_ode_and_limits():  # V1 (Ex. 4.7)
    sig, rho = 0.0728, 998.0
    d = np.sqrt(sig / (rho * G))
    for th in np.deg2rad([10.0, 30.0, 60.0, 85.0]):
        h = ch04.meniscus_height(th, sig, rho)
        assert h ** 2 == pytest.approx(2 * d ** 2 * (1 - np.sin(th)), rel=1e-14)
        x, zeta = ch04.meniscus_profile_ode(th, x_max=4 * d, sigma=sig, rho=rho)
        assert zeta[0] == pytest.approx(h, rel=1e-12)
        m = zeta > 1e-3 * h
        xc = ch04.meniscus_profile_x(zeta[m], th, sig, rho)
        assert np.max(np.abs(xc - x[m])) < 1e-7 * d * 10  # linear interpolation of the dense ODE output
        assert np.all(np.diff(zeta) < 0)
    assert ch04.meniscus_height(np.pi / 2) == pytest.approx(0.0, abs=1e-12)
    assert ch04.meniscus_height(0.0, sig, rho) == pytest.approx(np.sqrt(2) * d)


def test_meniscus_V2_closed_form_satisfies_first_integral():  # V2 (Ex. 4.7; the book's sign slip)
    zeta, d = sp.symbols("zeta delta", positive=True)
    h = sp.Symbol("h", positive=True)
    x = d * (sp.acosh(2 * d / zeta) - sp.sqrt(4 - zeta ** 2 / d ** 2) - sp.acosh(2 * d / h) + sp.sqrt(4 - h ** 2 / d ** 2))
    dxdz = sp.diff(x, zeta)
    slope2 = 1 / dxdz ** 2  # ζ′² = (dζ/dx)²
    first = zeta ** 2 / (2 * d ** 2) + 1 / sp.sqrt(1 + slope2)  # (ρg/2σ)ζ² + (1 + ζ′²)^{−1/2}
    f = sp.lambdify((zeta, d), sp.simplify(first))
    for zz in (0.1, 0.5, 1.0, 1.3):
        assert f(zz, 1.0) == pytest.approx(1.0, abs=1e-12)
        assert float(dxdz.subs({zeta: zz, d: 1.0})) < 0  # the slope is negative (the book's separated form drops a −)
    assert x.subs(zeta, h) == 0


def test_free_energy_V2_isothermal_work_and_spheroid_minimum():  # V2 + V1 (N121–N123)
    T, v, R, cv = sp.symbols("T v R c_v", positive=True)
    e = cv * T
    s = cv * sp.log(T) + R * sp.log(v)
    f = e - T * s
    assert sp.simplify(sp.diff(f, v) + R * T / v) == 0  # (∂f/∂v)_T = −p for a perfect gas
    assert ch04.helmholtz_free_energy(2.0e5, 300.0, 100.0) == pytest.approx(2.0e5 - 3.0e4)
    V = 1.0
    r = (3 * V / (4 * np.pi)) ** (1 / 3)
    assert ch04.spheroid_area(V, 1.0) == pytest.approx(4 * np.pi * r ** 2, rel=1e-13)
    for asp in (0.5, 0.9, 1.1, 2.0):
        assert ch04.spheroid_area(V, asp) > ch04.spheroid_area(V, 1.0)
        a = (3 * V / (4 * np.pi * asp)) ** (1 / 3)
        c = asp * a
        num = quad(lambda t_: 2 * np.pi * a * np.sin(t_) * np.sqrt(a ** 2 * np.cos(t_) ** 2 + c ** 2 * np.sin(t_) ** 2),
                   0, np.pi, epsrel=1e-13)[0]
        assert ch04.spheroid_area(V, asp) == pytest.approx(num, rel=1e-10)


# =====================================================================================================================
# C15 — dimensionless equations and dynamic similarity (4.99)–(4.119); D30
# =====================================================================================================================
def test_similarity_V2_ns_coefficients():  # V2
    l, U, rho, mu, g, Om = sp.symbols("l U rho mu g Omega", positive=True)
    c = ch04.nondimensional_ns_coefficients()
    expect = dict(unsteady=Om * l / U, advective=1, pressure=1, gravity=g * l / U ** 2, viscous=mu / (rho * U * l))
    assert all(sp.simplify(c[k] - v) == 0 for k, v in expect.items())
    cv = ch04.nondimensional_ns_coefficients("viscous")
    assert sp.simplify(cv["pressure"] - mu / (rho * U * l)) == 0
    ch = ch04.nondimensional_ns_coefficients("hydrostatic", "advective")
    assert sp.simplify(ch["pressure"] - g * l / U ** 2) == 0 and sp.simplify(ch["unsteady"] - 1) == 0
    eq = ch04.nondimensional_ns_sym()
    assert isinstance(eq, sp.Eq) and eq.has(Om * l / U)
    e = ch04.nondimensional_energy_coefficients()
    Uo, lo, ro, muo, ko, cp, To, Tw = sp.symbols("U l rho_o mu_o k_o C_p T_o T_w", positive=True)
    Ec = Uo ** 2 / (cp * (Tw - To))
    Re = ro * Uo * lo / muo
    Pr = muo * cp / ko
    assert sp.simplify(e["pressure_work"] - Ec) == 0 and sp.simplify(e["dissipation"] - Ec / Re) == 0
    assert sp.simplify(e["conduction"] - 1 / (Pr * Re)) == 0
    c_ = sp.Symbol("c", positive=True)
    assert sp.simplify(ch04.nondimensional_continuity_coefficient() - Uo ** 2 / c_ ** 2) == 0


def test_similarity_V2_derivation():  # V2 — D30, independent route: a concrete scaled field substituted into (4.39b)
    l, U, rho, mu, g, Om = sp.symbols("l U rho mu g Omega", positive=True)
    x, t, xs, ts = sp.symbols("x t x_s t_s", real=True)
    f = sp.sin(xs) * sp.exp(-ts) + xs ** 2  # u* = f(x*, t*)
    q = sp.cos(2 * xs) * ts  # p* = q(x*, t*)
    sub = {xs: x / l, ts: Om * t}
    u = U * f.subs(sub)
    p = sp.Symbol("p_inf") + rho * U ** 2 * q.subs(sub)
    dim = rho * sp.diff(u, t) + rho * u * sp.diff(u, x) + sp.diff(p, x) + rho * g - mu * sp.diff(u, x, 2)
    star = (Om * l / U) * sp.diff(f, ts) + f * sp.diff(f, xs) + sp.diff(q, xs) + g * l / U ** 2 - mu / (rho * U * l) * sp.diff(f, xs, 2)
    assert sp.simplify(dim / (rho * U ** 2 / l) - star.subs(sub)) == 0  # ∂/∂t = Ω∂/∂t*, ∇ = ∇*/l, ÷ ρU²/l


def test_similarity_V1_collapse_of_dimensional_solutions():  # V1 (§4.11: same equation, same solution)
    eta = np.linspace(0, 4, 41)
    curves = []
    for (U, nu, t) in [(1.0, 1e-6, 10.0), (3.5, 1.5e-5, 0.7), (0.2, 1e-4, 100.0)]:
        y = eta * np.sqrt(nu * t)
        curves.append(ch04.stokes_first_problem(y, t, U, nu) / U)
    assert max(np.max(np.abs(c - curves[0])) for c in curves) < 1e-12
    Y = np.linspace(0, 1, 21)
    P = [ch04.plane_poiseuille(Y * h, G_, h, mu) / (G_ * h ** 2 / mu) for (G_, h, mu) in [(100, 0.01, 1e-3), (5, 0.3, 1.8e-5)]]
    assert np.max(np.abs(P[0] - P[1])) < 1e-14
    # cylinder C_p is the same for any (U, a): dynamic similarity of potential flow
    th = np.linspace(0, np.pi, 13)
    cps = []
    for (U, a) in [(1, 1), (7, 0.2), (0.3, 5)]:
        uu, vv = ch04.velocity_preset("cylinder", 1.3 * a * np.cos(th), 1.3 * a * np.sin(th), U=U, a=a)
        cps.append(1 - (uu ** 2 + vv ** 2) / U ** 2)
    assert max(np.max(np.abs(c - cps[0])) for c in cps) < 1e-13


def test_named_numbers_V2_dimensionless_and_identities():  # V2 (pint) + V1 identities
    q = dict(U=Q_(2.0, "m/s"), l=Q_(0.5, "m"), nu=Q_(1e-6, "m**2/s"), g=Q_(9.81, "m/s**2"), Om=Q_(3.0, "1/s"),
             rho=Q_(1000.0, "kg/m**3"), mu=Q_(1e-3, "Pa*s"), sigma=Q_(0.07, "N/m"), cp=Q_(4182.0, "J/(kg*K)"),
             dT=Q_(10.0, "K"), c=Q_(1500.0, "m/s"), kap=Q_(1.4e-7, "m**2/s"))
    groups = [(lambda Om, l, U: Om * l / U, dict(Om=q["Om"], l=q["l"], U=q["U"])),
              (lambda U, l, nu: U * l / nu, dict(U=q["U"], l=q["l"], nu=q["nu"])),
              (lambda U, l, g: U / (g * l) ** 0.5, dict(U=q["U"], l=q["l"], g=q["g"])),
              (lambda U, c: U / c, dict(U=q["U"], c=q["c"])),
              (lambda U, cp, dT: U ** 2 / (cp * dT), dict(U=q["U"], cp=q["cp"], dT=q["dT"])),
              (lambda nu, kap: nu / kap, dict(nu=q["nu"], kap=q["kap"])),
              (lambda rho, U, l, sigma: rho * U ** 2 * l / sigma, dict(rho=q["rho"], U=q["U"], l=q["l"], sigma=q["sigma"])),
              (lambda rho, g, l, sigma: rho * g * l ** 2 / sigma, dict(rho=q["rho"], g=q["g"], l=q["l"], sigma=q["sigma"])),
              (lambda mu, U, sigma: mu * U / sigma, dict(mu=q["mu"], U=q["U"], sigma=q["sigma"])),
              (lambda U, Om, l: U / (Om * l), dict(U=q["U"], Om=q["Om"], l=q["l"]))]
    for fn, kw in groups:
        val = fn(**kw)
        assert val.dimensionless, (val, val.dimensionality)  # every group of (4.102)–(4.119) is a pure number
    assert ch04.strouhal_number(3.0, 0.5, 2.0) == pytest.approx(0.75)
    assert ch04.reynolds_number(1.0, 0.01, 1e-6) == pytest.approx(1e4) and ch04.reynolds_number(1.0, 0.01, rho=1000.0, mu=1e-3) == pytest.approx(1e4)
    with pytest.raises(ValueError):
        ch04.reynolds_number(1.0, 1.0)
    Re, We = ch04.reynolds_number(2.0, 0.5, rho=1000.0, mu=1e-3), ch04.weber_number(1000.0, 2.0, 0.5, 0.07)
    assert ch04.capillary_number(1e-3, 2.0, 0.07) == pytest.approx(We / Re, rel=1e-14)
    gp = ch04.reduced_gravity(1000.0, 1002.0)
    assert gp == pytest.approx(G * 2e-3)
    Fr_i = ch04.internal_froude_number(0.1, gp, 10.0)
    assert ch04.richardson_number(gp, 10.0, 0.1) == pytest.approx(1 / Fr_i ** 2, rel=1e-14)
    assert ch04.richardson_number(9.81e-3, 100.0, 0.1) == pytest.approx(98.1)
    N = 0.02
    assert ch04.internal_froude_number(0.1, N=N, l=10.0) == pytest.approx(ch04.internal_froude_number(0.1, N ** 2 * 10.0, 10.0))
    with pytest.raises(ValueError):
        ch04.internal_froude_number(0.1)
    assert ch04.gradient_richardson_number(1e-4, 0.01) == pytest.approx(1.0)
    assert ch04.froude_number(3.0, 1.0, 9.0) == pytest.approx(1.0)
    assert ch04.eckert_number(10.0, 1000.0, 5.0) == pytest.approx(0.02) and ch04.prandtl_number(1e-6, 1.4e-7) == pytest.approx(7.142857, rel=1e-6)
    assert ch04.rossby_number(10.0, 7.292115e-5, 1.0e6) == pytest.approx(0.06857, abs=5e-6)
    assert ch04.rossby_number(10.0, 1e-4, 1e5, factor=1.0) == pytest.approx(1.0)
    assert ch04.mach_number(170.0, c=340.0) == 0.5
    with pytest.raises(ValueError):
        ch04.mach_number(1.0)


def test_coefficients_V1_pressure_drag_lift_and_areas():  # V1 (4.106)–(4.108)
    assert ch04.pressure_coefficient(1e5 + 50, 1e5, 1.0, 10.0) == pytest.approx(1.0)
    F = ch04.drag_coefficient(12.0, 1.2, 5.0, 0.4) * 0.5 * 1.2 * 25 * 0.4
    assert F == pytest.approx(12.0) and ch04.lift_coefficient(3.0, 1.0, 1.0, 6.0) == pytest.approx(1.0)
    assert ch04.reference_area("sphere", d=0.2) == pytest.approx(np.pi * 0.01)
    assert ch04.reference_area("plate", b=2.0, c=0.5) == 1.0 and ch04.reference_area("airfoil", s=3.0, l=0.2) == pytest.approx(0.6)
    assert ch04.reference_area("cylinder", b=2.0, d=0.1) == pytest.approx(0.2)
    with pytest.raises(ValueError):
        ch04.reference_area("cube", d=1.0)


@needs_ref
def test_sphere_drag_V5_morrison_form_and_limits():  # V5 + V1
    ref = ref_json()["morrison_sphere_drag"]
    f = sp.lambdify(sp.Symbol("Re"), sp.sympify(ref["formula"]), "numpy")
    Re = np.logspace(-1, 6, 200)
    assert rel(ch04.sphere_drag_coefficient(Re), f(Re)) < 1e-12
    assert ch04.sphere_drag_coefficient(1e4) == pytest.approx(0.3926, abs=5e-5)
    assert ch04.sphere_drag_coefficient(0.1) / (24 / 0.1) == pytest.approx(1.0, abs=0.01)  # Stokes limit within 1 %
    assert ch04.sphere_drag_coefficient(0.1, model="stokes") == pytest.approx(240.0)
    with pytest.raises(ValueError):
        ch04.sphere_drag_coefficient(1.0, model="oseen")
    data = ch04.synthetic_sphere_drag_data(n=60, seed=1, noise=0.0)
    CD = ch04.drag_coefficient(data["F"], data["rho"], data["U"], np.pi * data["d"] ** 2 / 4)
    Re_ = ch04.reynolds_number(data["U"], data["d"], rho=data["rho"], mu=data["mu"])
    assert rel(Re_, data["Re"]) < 1e-12 and rel(CD, ch04.sphere_drag_coefficient(Re_)) < 1e-12  # collapse on C_D(Re)


def test_pi_groups_V1_sphere_drag_two_repeating_sets():  # V1 (R12, (4.99))
    g = ch04.sphere_drag_pi_groups()
    from fractions import Fraction
    s1 = {k: v for k, v in g["set1"][0].items() if v != 0}
    s2 = {k: v for k, v in g["set2"][0].items() if v != 0}
    assert s1 == {"F": 1, "U": -2, "D": -2, "rho": -1} and s2 == {"F": 1, "mu": -2, "rho": 1}
    # Fρ/μ² = (F/ρU²D²)·Re² exactly: exponent vectors add
    re = {"rho": Fraction(1), "U": Fraction(1), "D": Fraction(1), "mu": Fraction(-1)}
    combo = {k: s1.get(k, 0) + 2 * re.get(k, 0) for k in set(s1) | set(re)}
    assert {k: v for k, v in combo.items() if v != 0} == s2


def test_ship_model_V1_froude_scaling_ratios():  # V1 (Ex. 4.8, N156)
    assert ch04.froude_scaled_speed(10.0, 100.0, 4.0) == pytest.approx(2.0)
    d = ch04.ship_drag_extrapolation(scale=1 / 25, rho_p=1025.0)
    lam = 1 / 25
    assert d["D_p_wave"] / d["D_m_wave"] == pytest.approx(1025.0 / 1000.0 * lam ** -3, rel=1e-13)
    assert d["Re_ratio"] == pytest.approx(lam ** -1.5, rel=1e-13)
    assert d["D_p_total"] == pytest.approx(d["D_p_wave"] + d["D_p_friction"])
    m = ch04.model_prototype(100.0, 10.0, 1 / 25)
    assert m["matched"] == ["Fr"] and m["Re_ratio"] == pytest.approx(125.0, rel=1e-12)
    mr = ch04.model_prototype(100.0, 10.0, 1 / 25, match="Re")
    assert mr["matched"] == ["Re"] and mr["U_m"] == pytest.approx(250.0, rel=1e-12)
    mix = ch04.model_prototype(1.0, 30.0, 1.0, fluid_p="air", fluid_m="water", match="Re")
    assert mix["Re_p"] == pytest.approx(mix["Re_m"], rel=1e-12)
    with pytest.raises(ValueError):
        ch04.model_prototype(1.0, 1.0, 0.1, match="We")


def test_scales_V1_round_trip_and_oscillating_body():  # V1 + V7
    s = ch04.Scales(l=0.2, U=3.0, rho=1000.0, mu=1e-3, Omega=5.0, p_inf=1e5, c=1500.0, T_o=280.0, T_w=300.0, cp=4182.0, k=0.6)
    vals = dict(x=np.array([0.1, 0.4]), t=2.0, u=1.5, p=1.2e5, T=290.0)
    back = s.redimensionalise(**s.nondimensionalise(**vals))
    assert all(np.allclose(back[k], v, rtol=1e-15) for k, v in vals.items())
    assert s.St == pytest.approx(5.0 * 0.2 / 3.0) and s.Re == pytest.approx(3.0 * 0.2 / 1e-6)
    assert s.Fr == pytest.approx(3.0 / np.sqrt(G * 0.2)) and s.M == pytest.approx(2e-3) and s.Ec == pytest.approx(9 / (4182 * 20))
    assert s.Pr == pytest.approx(1e-3 * 4182.0 / 0.6) and set(s.groups()) == {"St", "Re", "Fr", "M", "Ec", "Pr"}
    assert s.t_ref == pytest.approx(0.2) and s.p_ref == pytest.approx(9000.0)
    o = ch04.Scales.from_oscillation(0.1, 4.0, 1000.0, 1e-3)
    assert o.St == pytest.approx(1.0) and o.Re == pytest.approx(4.0 * 0.01 / 1e-6)
    assert o.Fr == pytest.approx(4.0 * np.sqrt(0.1 / G))
    adv = ch04.Scales(l=1.0, U=2.0, rho=1.0, mu=1.0, pressure_scale="viscous")
    assert adv.time_scale == "advective" and adv.t_ref == 0.5 and adv.p_ref == 2.0 and adv.St is None
    with pytest.raises(ValueError):
        ch04.Scales(l=1.0, U=1.0, rho=1.0, mu=1.0, time_scale="omega")
    with pytest.raises(ValueError):
        ch04.Scales(l=1.0, U=1.0, rho=1.0, mu=1.0, pressure_scale="magic")


@needs_ref
def test_prandtl_V5_air_water_and_kinetic_theory():  # V5
    pr = ref_json()["prandtl"]
    lo, hi = pr["air_range_250_1000K"]
    for T in (250.0, 280.0, 300.0):
        assert lo <= ch04.prandtl_of("air", T) <= hi, T
    w = pr["water_300K"]
    assert abs(ch04.prandtl_of("water", 300.0) - w) <= 0.05 + 0.01 * w  # 2-digit value: rounding band + 1 %
    assert ch04.eucken_prandtl(5 / 3) == pytest.approx(float(sp.Rational(pr["monatomic"])), rel=1e-14)
    assert ch04.eucken_prandtl(1.4) == pytest.approx(5.6 / 7.6)
    assert ch04.prandtl_of("water", 293.15) > ch04.prandtl_of("water", 300.0)  # V7: Pr of water falls with T
    with pytest.raises(ValueError):
        ch04.prandtl_of("glycerine", 300.0)


# =====================================================================================================================
# V6 — the book's own numbers (private JSON; skipped when absent)
# =====================================================================================================================
@book_only
def test_book_values_V6_example_4_8_ship_model():  # V6
    b = book()["example_4_8_ship"]
    d = ch04.ship_drag_extrapolation(L_p=b["L_p_m"], U_p=b["U_p_m_s"], S_p=b["wetted_area_p_m2"], scale=b["scale"],
                                     D_m_total=b["model_total_drag_N"], CDf_m=b["CD_friction_model"],
                                     CDf_p=b["CD_friction_prototype"], rho_m=b["rho_kg_m3"], rho_p=b["rho_kg_m3"],
                                     nu=b["nu_m2_s"])
    pairs = [("U_m", "U_m_m_s"), ("Re_m", "Re_model"), ("Re_p", "Re_prototype"), ("D_m_friction", "friction_model_N"),
             ("D_m_wave", "wave_model_N"), ("D_p_wave", "wave_prototype_N"), ("D_p_friction", "friction_prototype_N"),
             ("D_p_uncorrected", "uncorrected_prototype_N")]
    for ours, key in pairs:
        assert d[ours] == pytest.approx(b[key], rel=5e-3), key
    # the printed total adds the book's rounded wave drag (analysis §9 item 12; explained in the report)
    assert d["D_p_total"] == pytest.approx(b["wave_prototype_N"] + b["friction_prototype_N"], rel=5e-3)
    assert abs(d["D_p_total"] / b["total_prototype_N"] - 1) < 1.5e-3


@book_only
def test_book_values_V6_closed_forms_of_examples():  # V6 (the examples' formulas, parsed from the private JSON)
    B = book()
    y, H, Ui, rho = sp.symbols("y H U_inf rho", positive=True)
    Uw = sp.Function("U")(y)
    wake = sp.sympify(B["example_4_1_wake_drag"]["drag_per_span"], locals={"U": Uw, "y": y, "H": H, "U_inf": Ui, "rho": rho})
    Fb = float(wake.subs(Uw, 10.0 - 2.0 * sp.exp(-y ** 2 / sp.Rational(1, 100))).subs({H: 2, Ui: 10, rho: sp.Rational(6, 5)}).evalf())
    assert Fb == pytest.approx(ch04.wake_drag_per_span("gaussian", 10.0, 1.2, 2.0), rel=1e-10)
    g, hi, ho = sp.symbols("g h_in h_out", positive=True)
    spd = sp.sympify(B["example_4_3_bore"]["speed"], locals=dict(g=g, h_in=hi, h_out=ho))
    assert float(spd.subs({g: G, hi: 1.0, ho: 1.1})) == pytest.approx(ch04.bore_speed(1.0, 1.1), rel=1e-14)
    net = sp.sympify(B["example_4_3_bore"]["net_side_pressure_force_per_width"], locals=dict(g=g, h_in=hi, h_out=ho, rho=rho))
    assert float(net.subs({g: G, hi: 1.0, ho: 1.3, rho: 1000})) == pytest.approx(ch04.bore_pressure_force(1.0, 1.3)["net"], rel=1e-12)
    a, A, U, al = sp.symbols("a A U alpha", positive=True)
    tq = sp.sympify(B["example_4_6_sprinkler"]["torque"], locals=dict(a=a, rho=rho, A=A, U=U, alpha=al))
    assert float(tq.subs({a: 0.2, rho: 1000, A: 1e-4, U: 5, al: sp.pi / 6})) == pytest.approx(0.8660254, rel=1e-7)
    th, sg = sp.symbols("theta sigma", positive=True)
    h2 = sp.sympify(B["example_4_7_meniscus"]["height_squared"], locals=dict(sigma=sg, rho=rho, g=g, theta=th))
    assert float(sp.sqrt(h2).subs({sg: 0.0728, rho: 998, g: G, th: 0.5})) == pytest.approx(ch04.meniscus_height(0.5), rel=1e-12)
    R1, R2, ze, dp = sp.symbols("R1 R2 zeta dp", positive=True)
    Fst = sp.sympify(B["surface_tension_cap_4_10"]["Fst_z"], locals=dict(sigma=sg, R1=R1, R2=R2, zeta=ze))
    assert float(Fst.subs({sg: 0.07, R1: 1e-3, R2: 2e-3, ze: 1e-6})) == pytest.approx(
        ch04.cap_surface_tension_force(0.07, 1e-3, 2e-3, 1e-6, exact=False), rel=1e-12)
    W, u, t = sp.symbols("Omega u t", positive=True)
    defl = sp.sympify(B["coriolis_projectile_4_7"]["deflection"], locals=dict(Omega=W, u=u, t=t))
    assert float(defl.subs({W: ch04.OMEGA_EARTH, u: 10, t: 3600})) == pytest.approx(
        ch04.coriolis_projectile(10.0, ch04.OMEGA_EARTH, [3600.0])["deflection_small"][0], rel=1e-12)


@book_only
def test_book_values_V6_quoted_numbers():  # V6
    B = book()
    assert ch04.earth_oblateness_diameter() / 1e3 == pytest.approx(B["earth_oblateness_4_7"]["equatorial_minus_polar_diameter_km"], rel=0.02)
    cap = B["capillary_scale_4_10"]
    assert ch04.capillary_length(cap["sigma_N_per_m"], cap["rho"], cap["g"]) == pytest.approx(cap["length_m"], rel=5e-3)
    th = B["incompressible_air_threshold_4_2"]  # a round speed and a Mach bound
    assert ch04.mach_number(th["speed_m_s"], T=288.15) < th["mach"]
    assert abs(ch04.incompressible_speed_limit(M=th["mach"]) - th["speed_m_s"]) < 0.05 * th["speed_m_s"]  # 102 m/s
    m = B["mach_threshold_4_11"]  # "M < 0.3 ⇒ departures below 10 %" is a bound: M² = 0.09 < 0.10
    assert ch04.compressibility_parameter(m["M_incompressible"], 1.0) < m["departure_fraction"]
    assert ch04.eucken_prandtl(5 / 3) == pytest.approx(B["prandtl_4_11"]["kinetic_theory_hard_sphere_monatomic"], rel=5e-3)
    assert B["orifice_4_9"]["sharp_edge_contraction_fraction"] == pytest.approx(0.611, rel=0.02)  # book value vs Wikipedia's 0.611 (pre-registered 2 %)
    c_air = 100.0 / ch04.mach_number(100.0, T=288.15)
    Hc_km = c_air ** 2 / G / 1e3  # 11.8 km; the book quotes an order of magnitude
    assert round(Hc_km, -1) == B["boussinesq_4_9"]["c2_over_g_air_km"]
    p = ch04.rotating_pump_terms()
    s = p["symbols"]
    rt = sp.sympify(B["example_4_5_pump"]["radial_rhs_rotation_terms"],
                    locals=dict(rho=s["rho"], Omega_z=s["Omega_z"], u_phi=s["u_phi"], R=s["R"]))
    assert sp.simplify(p["rotation_terms"][0] - rt) == 0


# =====================================================================================================================
# Part C contract: every function the notebook and explainers call exists, and is exercised in this file
# =====================================================================================================================
PART_C = """
MassBudget mass_budget MomentumBudget momentum_budget EnergyBudget energy_budget interval_mass_budget angular_momentum_flux
velocity_from_streamfunction_2d velocity_from_streamfunction_2d_sym flux_between_streamlines flux_along_path
velocity_from_streamfunction_axisym velocity_from_streamfunction_axisym_sym mass_flux_from_vector_potential
mass_flux_from_stream_functions stream_function_pair stream_surface_check stream_tube_mass_flux streamfunction_preset
velocity_preset STREAMFUNCTION_PRESETS static_stress total_stress linear_stress isotropic_fourth_order newtonian_stress
viscous_stress mean_pressure thermodynamic_pressure_from_stress pressure_difference bulk_viscosity lam_from_bulk
stress_on_plane dissipation_rate deviatoric_part continuity_residual continuity_residual_sym ContinuityTerms
continuity_terms continuity_material_terms density_material_rate divergence_free_check momentum_conservative_residual
conservative_to_advective_sym CauchyTerms cauchy_terms cauchy_residual navier_stokes_residual navier_stokes_sym NSTerms
ns_incompressible_terms viscous_force_forms exact_solution exact_field ns_terms_preset lamb_vector lamb_identity_terms
lamb_identity_sym stress_work_split kinetic_energy_budget internal_energy_terms internal_energy_residual
total_energy_residual_sym energy_forms_sym entropy_terms entropy_production perturbation_fields buoyancy
boussinesq_momentum_terms heat_equation_terms temperature_equation_residual divergence_first_index_demo OMEGA_EARTH
EARTH_RADIUS rotating_basis basis_rate inertial_velocity frame_acceleration_terms apparent_body_forces
coriolis_acceleration coriolis_parameter centrifugal_acceleration centrifugal_potential effective_gravity projectile_paths
bernoulli_head bernoulli_solve pitot_speed pitot_speed_from_heads stagnation_pressure dynamic_pressure torricelli_speed
pressure_function bernoulli_function bernoulli_along_line lamb_surface_check unsteady_bernoulli_pressure
unsteady_bernoulli_B gauge_absorbed_bracket unsteady_streamline_bernoulli viscous_irrotational_residual
stagnation_enthalpy stagnation_temperature BERNOULLI_FORMS which_bernoulli which_bernoulli_text bernoulli_scenario
surface_normal_speed kinematic_bc_residual relative_normal_velocity interface_mass_flux pillbox_limit cap_pressure_force
cap_surface_tension_force laplace_jump_from_balance capillary_length Scales nondimensional_ns_coefficients
nondimensional_ns_sym nondimensional_energy_coefficients strouhal_number reynolds_number froude_number reduced_gravity
internal_froude_number richardson_number gradient_richardson_number mach_number compressibility_parameter eckert_number
prandtl_number eucken_prandtl weber_number bond_number capillary_number rossby_number pressure_coefficient
drag_coefficient lift_coefficient reference_area sphere_drag_coefficient froude_scaled_speed model_prototype
helmholtz_free_energy closure_count expanding_flow material_mass material_interval stratified_shear_flow
is_incompressible_regime ball_integral body_force_from_potential gravity_potential gaussian_wake wake_side_outflow
wake_drag_per_span bore_speed bore_pressure_force rocket_trajectory rocket_delta_v jet_plate_force cv_scenario
stream_tube_element_balance_sym orifice_mass_flow tank_drain sprinkler_torque sprinkler_free_spin_rate
cube_spin_acceleration stokes_first_problem plane_poiseuille coriolis_projectile rotating_pump_equations high_low_flow
couette_heating couette_heating_transient rankine_bernoulli u_tube_column accelerating_sphere_pressure
boussinesq_validity boussinesq_scenario boussinesq_density gaussian_blob_advection_diffusion blob_rise
two_layer_conduction two_fluid_couette navier_slip_couette linear_wave_surface wave_kinematic_residual spheroid_area
meniscus_height meniscus_profile_x meniscus_profile_ode prandtl_of ship_drag_extrapolation
""".split()
PART_C_CU = ["gradient", "divergence", "curl", "laplacian", "vector_laplacian", "advective_acceleration", "strain_rate",
             "coordinates"]


def test_contract_V1_every_part_c_name_exists_and_is_reexported():  # V1 (design Part C)
    missing = [n for n in PART_C if not hasattr(ch04, n)]
    assert not missing, missing
    assert all(callable(getattr(ch04.CU, n)) for n in PART_C_CU)
    from fluidpy.core import bernoulli, constitutive, navier_stokes, rotating, similarity
    assert ch04.newtonian_stress is constitutive.newtonian_stress and ch04.which_bernoulli is bernoulli.which_bernoulli
    assert ch04.cauchy_terms is navier_stokes.cauchy_terms and ch04.coriolis_acceleration is rotating.coriolis_acceleration
    assert ch04.Scales is similarity.Scales
    # the contracted leading parameters
    sig = {n: list(inspect.signature(getattr(ch04, n)).parameters)[:3] for n in
           ("mass_budget", "stress_on_plane", "ns_terms_preset", "frame_acceleration_terms", "couette_heating")}
    assert sig["mass_budget"] == ["rho", "u", "cv"] and sig["stress_on_plane"] == ["G", "p", "mu"]
    assert sig["ns_terms_preset"] == ["name", "x", "y"] and sig["frame_acceleration_terms"] == ["a_prime", "u_prime", "x_prime"]
    assert sig["couette_heating"] == ["y", "U", "h"]


def test_contract_V1_every_part_c_name_is_exercised_in_this_file():  # V1
    full = Path(__file__).read_text(encoding="utf-8")
    start = full.index('PART_C = """')
    src = full[:start] + full[full.index('""".split()', start):]  # everything except the name list itself
    unused = [n for n in PART_C if not re.search(rf"ch04\.{n}\b", src)]
    assert not unused, unused
    unused_cu = [n for n in PART_C_CU if not re.search(rf"CU\.{n}\b", src)]
    assert not unused_cu, unused_cu


def test_result_types_V1_named_fields_and_constants():  # V1 (the result classes the notebook teaches)
    rho, u = ch04.expanding_flow_fields(a=1.0, rho0=1.0, dim=3)
    box = ch04.MovingBox()
    mb = ch04.mass_budget(rho, u, box, 0.0)
    assert isinstance(mb, ch04.MassBudget) and mb.residual == pytest.approx(mb.storage + mb.outflux)
    mo = ch04.momentum_budget(1.0, lambda X, T: 0.0 * X, box, 0.0)
    assert isinstance(mo, ch04.MomentumBudget) and np.allclose(mo.residual, [0, 0, 9.81])  # at rest, no surface force
    eb = ch04.energy_budget(1.0, lambda X, T: 0.0 * X, 1e3, box, 0.0)
    assert isinstance(eb, ch04.EnergyBudget) and eb.residual == 0.0
    ct = ch04.continuity_terms(rho, u, np.zeros(3), 0.0)
    assert isinstance(ct, ch04.ContinuityTerms) and ct._fields == ("local", "flux_divergence", "residual")
    ca = ch04.cauchy_terms(1.0, lambda X, T: 0.0 * X, lambda X, T: np.zeros((3, 3) + np.shape(X)[1:]), (0, 0, -9.81),
                           np.zeros(3))
    assert isinstance(ca, ch04.CauchyTerms) and np.allclose(ca.residual, [0, 0, 9.81])
    nt = ch04.ns_incompressible_terms(lambda X, T: 0.0 * X, lambda X, T: 1e5 - 1000 * 9.81 * X[2], np.zeros(3))
    assert isinstance(nt, ch04.NSTerms) and np.allclose(nt.residual, 0.0, atol=1e-8)  # hydrostatic rest
    assert set(ch04.STREAMFUNCTION_PRESETS) == {"uniform", "stagnation", "source_stream", "cylinder", "vortex", "shear",
                                                "axisym_uniform"}
    assert ch04.EARTH_RADIUS == ch04.EARTH_A == 6378137.0
    uu, rr = ch04.expanding_flow(np.array([1.0, 2.0]), 1.0, a=1.0, rho0=2.0)
    assert np.allclose(uu, [0.5, 1.0]) and np.allclose(rr, 1.0)
    u3, r3 = ch04.expanding_flow(np.ones((3, 4)), 1.0, dim=3)
    assert u3.shape == (3, 4) and np.allclose(r3, 1 / 8)
    with pytest.raises(ValueError):
        ch04.expanding_flow(1.0, -2.0)


def test_energy_V2_derivation_one_dimensional_chain():  # V2 — D19–D22 re-derived independently in 1-D
    x, t = sp.symbols("x t", real=True)
    rho, u, e, p, q = (sp.Function(n)(x, t) for n in ("rho", "u", "e", "p", "q"))
    mu, muv, g = sp.symbols("mu mu_v g", real=True)
    ux = sp.diff(u, x)
    sig = sp.Rational(4, 3) * mu * ux + muv * ux  # 1-D Newtonian σ₁₁ = 2μ(u_x − ⅓u_x) + μ_v u_x (4.59)
    tau = -p + sig
    D = lambda F: sp.diff(F, t) + u * sp.diff(F, x)  # noqa: E731
    C = sp.diff(rho, t) + sp.diff(rho * u, x)  # continuity (4.7)
    E = e + u ** 2 / 2
    # D19: (4.46) over a fixed interval [a, b] by the RTT + Gauss ⇒ the local form (4.53)
    a_, b_ = sp.symbols("a b", real=True)
    integral = (sp.integrate(sp.diff(rho * E, t), (x, a_, b_)) + (rho * E * u).subs(x, b_) - (rho * E * u).subs(x, a_)
                - sp.integrate(rho * g * u, (x, a_, b_)) - ((tau * u).subs(x, b_) - (tau * u).subs(x, a_))
                + (q.subs(x, b_) - q.subs(x, a_)))
    R53 = sp.diff(rho * E, t) + sp.diff(rho * E * u, x) - rho * g * u - sp.diff(tau * u, x) + sp.diff(q, x)
    # (4.49)–(4.52): all under one integral ⇔ d/db of the interval statement is the integrand at b (Gauss in 1-D)
    assert sp.simplify(sp.expand(sp.diff(integral, b_) - R53.subs(x, b_))) == 0
    assert sp.simplify(sp.expand(sp.diff(integral, a_) + R53.subs(x, a_))) == 0
    # D20: (4.53) − E × continuity = ρDE/Dt − … (4.55)
    R55 = rho * D(E) - rho * g * u - sp.diff(tau * u, x) + sp.diff(q, x)
    assert sp.simplify(sp.expand(R53 - R55 - E * C)) == 0
    # D21: u × Cauchy (4.24) = the mechanical-energy equation (4.56)
    cauchy = rho * D(u) - rho * g - sp.diff(tau, x)
    R56 = rho * D(u ** 2 / 2) - rho * g * u + u * sp.diff(p, x) - u * sp.diff(sig, x)
    assert sp.simplify(sp.expand(u * cauchy - R56)) == 0
    # D22: (4.55) − (4.56) = ρ[De/Dt + p Dv/Dt − σ u_x/ρ + q_x/ρ] + (p/ρ) × continuity (4.57)
    R57 = D(e) + p * D(1 / rho) - sig * ux / rho + sp.diff(q, x) / rho
    assert sp.simplify(sp.expand(R55 - R56 - rho * R57 - p / rho * C)) == 0
    # the force-work terms cancel, the deformation work stays: −p u_x + σ u_x (4.54)
    split = sp.diff(tau * u, x) - ((-p * ux + sig * ux) + (-u * sp.diff(p, x) + u * sp.diff(sig, x)))
    assert sp.simplify(sp.expand(split)) == 0


def test_energy_budget_V2_derivation_box_gauss():  # V2 — D19 in 3-D: the surface integrals of (4.48) by Gauss
    x, y, z = sp.symbols("x y z", real=True)
    coords = (x, y, z)
    u = [x * y, z ** 2, x + y * z]
    rhoE = 1 + x ** 2 * z
    tau = sp.Matrix(3, 3, lambda i, j: (i + 2 * j + 1) * x ** i * y + z * (j + 1))
    qv = [y * z, x ** 2, z * x]
    surf = 0
    for i, v in enumerate(coords):
        others = [(w, 0, 1) for w in coords if w != v]
        flux_i = rhoE * u[i] - sum(tau[i, j] * u[j] for j in range(3)) + qv[i]  # (ρE u − τ·u + q)·e_i
        surf += sp.integrate(flux_i.subs(v, 1) - flux_i.subs(v, 0), *others)
    vol = sp.integrate(sum(sp.diff(rhoE * u[i] - sum(tau[i, j] * u[j] for j in range(3)) + qv[i], coords[i])
                           for i in range(3)), (x, 0, 1), (y, 0, 1), (z, 0, 1))
    assert sp.simplify(surf - vol) == 0  # (4.49)–(4.51) with dV (the book's (4.51) prints dA — typo)


def test_scripts_V1_drawing_helpers_run():  # V1 smoke (C.12)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    from scripts.ch04_drawings import (budget_bars, cone_sweep, curved_cap, cv_box, face_flux_arrows, pillbox,  # noqa
                                       rotating_frames, stream_tube_element)
    fig = plt.figure()
    ax = fig.add_subplot(121)
    ax3 = fig.add_subplot(122, projection="3d")
    cv_box(ax, 0, 1, 0, 1, b=0.3)
    s = ch04.cv_scenario("wake")
    face_flux_arrows(ax, s["faces"])
    budget_bars(ax, {"storage": 0.0, "outflux": s["mass_out"]})
    rotating_frames(ax, 0.5)
    pillbox(ax, 0.2)
    stream_tube_element(ax3)
    cone_sweep(ax3)
    curved_cap(ax3)
    plt.close("all")


def test_scalar_callable_V1_explainer_parity_functions_return_floats():  # V1 (selftest parity rows need floats)
    vals = [ch04.bore_speed(1.0, 1.1), ch04.torricelli_speed(1.0), ch04.pitot_speed(500.0, 0.0, 1.2),
            ch04.cube_spin_acceleration(1.0, 0.0, 1000.0, 0.01), ch04.stagnation_temperature(300.0, 100.0),
            ch04.capillary_length(0.0728, 998.0), ch04.meniscus_height(0.3), ch04.sphere_drag_coefficient(1e4),
            ch04.reynolds_number(1.0, 0.01, 1e-6), ch04.froude_scaled_speed(10.0, 100.0, 4.0),
            ch04.flux_between_streamlines("cylinder", (1.5, 0.2), (2.0, 1.0)), ch04.gravity_potential(1.0),
            ch04.gauge_absorbed_bracket(2.0), ch04.coriolis_parameter(0.5), ch04.centrifugal_potential(1.0, 1.0),
            ch04.wake_drag_gaussian(10.0, 2.0, 0.1, 1.2), ch04.dissipation_rate([[0, 1, 0], [0, 0, 0], [0, 0, 0]], 1.0, 1.0),
            ch04.buoyancy(-0.5, 1000.0), ch04.bernoulli_function(2.0, 1e5, 1.0), ch04.effective_gravity(0.5)[0],
            ch04.stress_on_plane("shear", 0.0, 1e-3, 0.0, 0.3)[1], ch04.laplace_jump_from_balance(0.07, 1e-3, 2e-3)]
    assert all(isinstance(v, float) for v in vals), [type(v) for v in vals]
