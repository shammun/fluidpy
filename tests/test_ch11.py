"""Verification suite for Chapter 11 — Instability (Kundu, Cohen & Dowling 5e, §§11.1–11.14).

Evidence levels (``verify-implementation`` skill): V1 analytic · V2 symbolic (sympy / dimensions) · V3 convergence ·
V4 conservation / invariant / integral identity · V5 published benchmark (``reference/ch11/``, see SOURCES.md) · V6 book
value (private, ``tests/book_values_ch11.json``, skipped when absent) · V7 limits / symmetry / invariance / theorems.
Every test name is ``test_<concept>_<level>_<what>``; the comment on the ``def`` line repeats the level and the curation ID.

A items (CORE, ≥ 2 independent levels): C01 normal modes (11.1) · C02 Kelvin–Helmholtz (11.18) · C03 Bénard amplitude
problem (11.36)–(11.37) · C04 rigid–rigid neutral curve, Ra_c · C05 free–free (11.44) · C06 salt fingers (11.46) · C07 Taylor
(11.51)–(11.54) · C08 Taylor–Goldstein (11.61) · C09 Miles–Howard (11.67) · C10 Howard's semicircle (11.72) · C11
Orr–Sommerfeld (11.79) · C12 Rayleigh's inflection-point theorem (11.83)–(11.84) · C13 plane Poiseuille / Table 11.1 · C14 the
disturbance-energy equation (11.88) · C15 Lorenz (11.91).
Derivations: every ★★ / ★★★ D row (D02 D03 D04 D06 D08 D10 D11 D12 D13 D14 D17 D18 D19 D21 D22 D23 D24 D25) is re-derived
with sympy in a ``test_*_V2_derivation`` test, independently of the chapter's own sympy engines (which are tested
separately), and the intermediate lines of design Part F are checked where cheap.

Printed slips kept as wrong variants that must FAIL: S1 ``benard_determinant(printed=True)`` · S2
``benard_free_free_mode(printed=True)`` · S3 ``benard_free_free_sympy()["printed_root"]`` · S4/S11
``stratified_shear_sympy`` · S6 ``tollmien_profile(printed=True)`` · S7 ``taylor_galerkin_Ta(printed=True)`` · S12
``taylor_perturbation_sympy()["printed_continuity_units_ok"]``.

Caches (ch10 lesson: a cached table written by the same code can hide a bug): every critical point, neutral curve and
table in ``reference/ch11`` is checked against an **independent** analytic/published value and against a **live**
recomputation (``cache=False``) at a few points; the heavy table writers run in the slow tests.

Run: ``.venv/Scripts/python.exe -m pytest tests/test_ch11.py -q -p no:cacheprovider``  (``-m "not slow"`` skips the script
smoke runs and the table writers).
"""
from __future__ import annotations

import cmath
import json
import math
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import scipy.linalg as sla
import sympy as sp
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq, minimize_scalar

from fluidpy import ch11_instability as ch11
from fluidpy.core import laminar as LAM
from fluidpy.core import stability as ST
from fluidpy.core import waves as WV
from fluidpy.core.units import Q_, dimensional_check
from tools.convergence import observed_order

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference" / "ch11"
BOOK = Path(__file__).resolve().parent / "book_values_ch11.json"
book_only = pytest.mark.skipif(not BOOK.exists(), reason="book values are private; see CLAUDE.md rule 9")
slow = pytest.mark.slow
G0 = ch11.G0
PI = math.pi
TAU_SW = 1.5e-9 / 1.4e-7  # κ_s/κ for the design's seawater numbers (0.0107143)


def book():
    return json.loads(BOOK.read_text(encoding="utf-8"))


def rel(a, b):
    return abs(complex(a) - complex(b)) / max(abs(complex(b)), 1e-300)


def z0(expr) -> bool:
    return sp.simplify(sp.expand(expr)) == 0


def read_csv(name: str) -> dict:
    lines = [ln for ln in (REF / name).read_text(encoding="utf-8").splitlines() if ln.strip() and not ln.startswith("#")]
    hdr = lines[0].split(",")
    data = np.array([[float(v) for v in ln.split(",")] for ln in lines[1:]])
    return {h: data[:, j] for j, h in enumerate(hdr)}


POIS = ch11.parallel_profile("poiseuille")
TANH = ch11.parallel_profile("tanh")
BICK = ch11.parallel_profile("bickley")
ORSZAG_C = complex(0.23752649, 0.00373967)  # Orszag (1971) Table 2/3, Re = 1e4, k = 1 (reference/ch11/SOURCES.md)


# ======================================================================================================================
# Chebyshev toolkit (core.stability) — used by every solver
# ======================================================================================================================
def test_cheb_V1_exact_for_polynomials_and_nodes():  # V1 (C03 primer)
    D, x = ST.cheb(4)
    assert np.allclose(x, [1, math.sqrt(0.5), 0, -math.sqrt(0.5), -1], atol=1e-15)
    assert np.max(np.abs(D @ x ** 3 - 3 * x ** 2)) < 1e-13
    for N in (6, 12, 20):
        D, x = ST.cheb(N, (0.0, 2.0))
        for p in range(1, N + 1):
            err = np.max(np.abs(D @ x ** p - p * x ** (p - 1))) / max(1.0, p * 2.0 ** (p - 1))
            assert err < 1e-11, (N, p, err)
        assert np.max(np.abs(D @ np.ones(N + 1))) < 1e-12  # rows annihilate constants
    M = ST.cheb_matrices(10, orders=(1, 2, 4))
    D, _ = ST.cheb(10)
    assert np.allclose(M["D2"], D @ D) and np.allclose(M["D4"], D @ D @ D @ D)
    with pytest.raises(ValueError):
        ST.cheb(0)


def test_cheb_V3_dirichlet_eigenvalues_converge_spectrally():  # V3 (C03 primer: digits gained, not an order)
    errs = []
    for N in (8, 16, 32):
        D, _ = ST.cheb(N)
        lam = np.sort(np.linalg.eigvals(-(D @ D)[1:-1, 1:-1]).real)[:3]
        errs.append(np.max(np.abs(lam - (np.arange(1, 4) * PI / 2) ** 2)))
    assert errs[0] > 1e-6 and errs[2] < 1e-10 and errs[1] < errs[0] * 1e-3, errs


def test_clenshaw_curtis_V1_exact_up_to_degree_N():  # V1
    for N in (8, 9, 16):
        w = ST.clenshaw_curtis_weights(N)
        _, x = ST.cheb(N)
        assert abs(w.sum() - 2.0) < 1e-14
        for p in range(0, N + 1):
            exact = 0.0 if p % 2 else 2.0 / (p + 1)
            assert abs(np.sum(w * x ** p) - exact) < 1e-13, (N, p)
    w = ST.clenshaw_curtis_weights(12, (0.0, 3.0))
    _, x = ST.cheb(12, (0.0, 3.0))
    assert abs(np.sum(w * x ** 2) - 9.0) < 1e-12


def test_cheb_grid_V1_mapped_derivatives_and_quadrature():  # V1 (maps of the unbounded / semi-infinite solvers)
    g = ST.cheb_grid(120, (-1, 1), "tan", y_max=30.0, s=1.0)
    assert np.max(np.abs(g.D1 @ np.tanh(g.y) - 1 / np.cosh(g.y) ** 2)) < 1e-8
    assert np.max(np.abs(g.D2 @ np.tanh(g.y) + 2 * np.tanh(g.y) / np.cosh(g.y) ** 2)) < 1e-7
    assert abs(g.integrate(1 / np.cosh(g.y) ** 2) - 2 * math.tanh(30.0)) < 1e-9
    a = ST.cheb_grid(80, (0.0, 1.0), "algebraic", y_max=40.0, s=2.0)
    assert abs(a.y[0] - 40.0) < 1e-12 and abs(a.y[-1]) < 1e-12
    assert np.max(np.abs(a.D1 @ np.exp(-a.y) + np.exp(-a.y))) < 1e-8
    assert abs(a.integrate(np.exp(-a.y)) - (1 - math.exp(-40.0))) < 1e-8
    with pytest.raises(ValueError):
        ST.cheb_grid(10, (0, 1), "tan")


def _os_rowreplace(k, Re, N):
    """From scratch (the C03/C11 teaching route): row-replacement BCs + generalized_eigs — independent of constrained_eig."""
    D, y = ST.cheb(N)
    I = np.eye(N + 1)
    D2 = D @ D
    D4 = D2 @ D2
    L = D2 - k ** 2 * I
    U, Upp = 1 - y ** 2, -2.0 * np.ones_like(y)
    A = np.diag(U) @ L - np.diag(Upp) - (D4 - 2 * k ** 2 * D2 + k ** 4 * I) / (1j * k * Re)
    rows = [(0, I[0]), (N, I[N]), (1, D[0]), (N - 1, D[N])]
    A, B = ST.apply_bc_rows(A, L, rows)
    assert np.all(B[[0, N, 1, N - 1]] == 0)
    return ST.generalized_eigs(A, B, sort="imag", c_max=5.0)


def test_apply_bc_rows_V1_row_replacement_route_reproduces_orszag():  # V1 independent BC route (C11 from scratch)
    c = _os_rowreplace(1.0, 1e4, 80)
    assert abs(c[0] - ORSZAG_C) < 2e-8
    c_constrained = ST.orr_sommerfeld_eigs(1.0, 1e4, POIS["U"], POIS["Upp"], N=80, filter=False)
    assert abs(c[0] - c_constrained[0]) < 1e-9


def test_converged_eigs_V3_keeps_the_TS_mode_and_drops_spurious_ones():  # V3 (N-convergence filter, design C.1 1.6)
    fn = lambda N: _os_rowreplace(1.0, 1e4, N)  # noqa: E731
    raw = fn(60)
    kept = ST.converged_eigs(fn, 60, factor=1.5, tol=1e-6)
    assert abs(kept[0] - ORSZAG_C) < 1e-7
    assert 0 < len(kept) < len(raw)
    mask = ST.converged_mask(np.array([1.0, 2.0 + 1j]), np.array([1.0 + 1e-9, 5.0]))
    assert list(mask) == [True, False]
    assert not ST.converged_mask(np.array([1.0]), np.array([])).any()


def test_constrained_eig_V7_singular_constraint_raises():  # V7 guard
    with pytest.raises(np.linalg.LinAlgError):
        ST.constrained_eig(np.eye(3), np.eye(3), np.array([[1.0, 0.0, 0.0]]), [1])


# ======================================================================================================================
# C01 — normal modes (11.1) and the stability vocabulary
# ======================================================================================================================
def test_normal_mode_V1_two_forms_of_11_1_agree_on_a_field():  # V1 (C01, D01)
    rng = np.random.default_rng(0)
    x, y, t = np.meshgrid(np.linspace(0, 5, 21), np.linspace(-2, 2, 11), np.linspace(0, 3, 7), indexing="ij")
    for _ in range(5):
        k, m = rng.uniform(0.1, 2.0, 2)
        c = complex(*rng.uniform(-1, 1, 2))
        uh = complex(*rng.normal(size=2))
        K = math.hypot(k, m)
        s = ST.sigma_from_c(K, c)
        assert abs(s - (-1j * K * c)) < 1e-15  # (11.1): σ = −i|K|c
        assert abs(s.real - K * c.imag) < 1e-14 and abs(s.imag + K * c.real) < 1e-14  # D01 steps 4–5
        assert abs(ST.c_from_sigma(K, s) - c) < 1e-14
        f1 = ST.normal_mode(x, y, t, uh, k, m, sigma=s)
        f2 = ST.normal_mode(x, y, t, uh, k, m, c=c)
        assert np.max(np.abs(f1 - f2)) < 1e-12 * max(1.0, np.max(np.abs(f1)))
        assert np.allclose(f1, WV.real_field(uh * np.exp(s.real * t), k * x + m * y + s.imag * t), atol=1e-12)  # ch07 parity
    assert ST.sigma_from_c(2.0, 1 + 0.5j) == pytest.approx(1 - 2j)  # design C.1 1.8 example
    with pytest.raises(ValueError):
        ST.c_from_sigma(0.0, 1.0)
    with pytest.raises(ValueError):
        ST.normal_mode(0, 0, 0, 1, 1.0)


def test_stability_vocabulary_V7_class_verdict_and_marginal_type():  # V7 (C01, N04, N06)
    assert ST.stability_class(sigma=-0.1) == "stable" and ST.stability_class(sigma=0.0) == "neutral"
    assert ST.stability_class(sigma=1e-6) == "unstable" and ST.stability_class(c=0.3 + 1e-3j) == "unstable"
    assert ST.stability_class(c=[0.3 - 0.1j, 0.2 - 1e-3j]) == "stable"
    v = ST.stability_verdict([-0.3, -0.1, 0.02, -0.5])
    assert v["verdict"] == "unstable" and v["k_index_max"] == 2 and v["max_growth"] == pytest.approx(0.02)
    assert ST.stability_verdict([-1.0, -0.2])["verdict"] == "stable"
    assert ST.stability_verdict([-1.0, 0.0])["verdict"] == "neutral"
    assert ST.marginal_type(0.0 + 0.0j) == "stationary" and ST.marginal_type(1e-12 + 3.0j) == "oscillatory"
    with pytest.raises(ValueError):
        ST.marginal_type(0.1 + 1j)
    with pytest.raises(ValueError):
        ST.stability_verdict([])


def test_potential_well_V4_energy_never_increases_and_the_dimple_escapes():  # V4 (N02, Fig. 11.1)
    for shape, x0 in (("bowl", 1.0), ("cap", 0.01), ("plane", 0.3), ("dimple", 0.3), ("dimple", 1.2)):
        r = ch11.potential_well_demo(shape, x0, v0=0.0, damping=0.3, t_end=20.0)
        E = 0.5 * r["v"] ** 2 + r["V"](r["x"])
        assert np.all(np.diff(E) <= 1e-9), shape
    assert not ch11.potential_well_demo("dimple", 0.3)["escaped"]
    assert ch11.potential_well_demo("dimple", 1.2)["escaped"]
    assert ch11.potential_well_demo("cap", 0.01, t_end=40.0)["escaped"]
    b = ch11.potential_well_demo("bowl", 1.0, damping=0.0, t_end=2 * PI, n=201)
    assert abs(b["x"][-1] - 1.0) < 1e-8  # undamped bowl: period 2π


def test_normal_mode_growth_V1_kh_interface_and_benard_branches():  # V1 (C01, E1 mirror)
    s = ch11.normal_mode_growth("kh", 1.0, U1=6.0, U2=0.0, rho1=1.0, rho2=3.0, g=10.0)
    assert abs(s - complex(math.sqrt(1.75), -1.5)) < 1e-14  # σ = −ik(1.5 + i√1.75)
    s2 = ch11.normal_mode_growth("interface", 20.0, rho1=1.2, rho2=1000.0, g=9.81)
    assert abs(abs(s2.imag) - WV.interface_omega(20.0, 1.2, 1000.0, g=9.81)) < 1e-12 and s2.real == 0.0
    s3 = ch11.normal_mode_growth("interface", 20.0, rho1=1000.0, rho2=1.2, g=9.81)  # Rayleigh–Taylor
    assert s3.real > 0 and abs(s3.real - math.sqrt(9.81 * 20 * (1000 - 1.2) / 1001.2)) < 1e-10
    s4 = ch11.normal_mode_growth("benard_free", PI / math.sqrt(2), Ra=2000.0, Pr=1.0)
    assert abs(s4 - ch11.benard_free_free_sigma(PI / math.sqrt(2), 2000.0, 1.0)[0]) < 1e-14
    with pytest.raises(ValueError):
        ch11.normal_mode_growth("bogus", 1.0)


# ======================================================================================================================
# C02 — Kelvin–Helmholtz (11.18)
# ======================================================================================================================
def test_kh_phase_speed_V1_both_roots_solve_the_quadratic_on_a_random_field():  # V1 (C02, N15)
    rng = np.random.default_rng(1)
    k = rng.uniform(0.05, 50.0, 400)
    U1, U2 = rng.uniform(-10, 10, 2)
    for r1, r2 in ((1.2, 1000.0), (1.0, 3.0), (2.0, 2.0), (3.0, 1.0)):
        cp, cm = ch11.kh_phase_speed(k, U1, U2, r1, r2, g=9.81)
        for c in (cp, cm):
            res = r1 * (U1 - c) ** 2 + r2 * (U2 - c) ** 2 - (9.81 / k) * (r2 - r1)  # p. 479
            scale = r1 * (U1 - c.real) ** 2 + r2 * (U2 - c.real) ** 2 + (9.81 / k) * abs(r2 - r1) + 1.0
            assert np.max(np.abs(res) / scale) < 1e-12
        assert np.all(np.imag(cp) >= 0) and np.allclose(cp + cm, 2 * (r1 * U1 + r2 * U2) / (r1 + r2))
    cp, cm = ch11.kh_phase_speed(1.0, 6.0, 0.0, 1.0, 3.0, g=10.0)
    assert abs(cp - (1.5 + 1j * math.sqrt(1.75))) < 1e-14 and abs(cm - (1.5 - 1j * math.sqrt(1.75))) < 1e-14
    cp, cm = ch11.kh_phase_speed(1.0, 4.0, 0.0, 1.0, 3.0, g=10.0)  # neutral pair (design example)
    assert abs(cp - (1 + math.sqrt(2))) < 1e-14 and abs(cm - (1 - math.sqrt(2))) < 1e-14


def test_kh_phase_speed_V7_limits_galilean_shift_and_signs():  # V7 (C02, R09, R10, N19)
    k = np.geomspace(0.5, 500, 40)
    cp, cm = ch11.kh_phase_speed(k, 0.0, 0.0, 1.2, 1000.0, g=9.81)  # (11.19) = ch07 (7.95)
    assert np.allclose(cp.real, WV.interface_omega(k, 1.2, 1000.0, g=9.81) / k, rtol=1e-13)
    assert np.allclose(cp.real, -cm.real) and np.all(cp.imag == 0)
    cp, _ = ch11.kh_phase_speed(k, 1.0, 3.0, 5.0, 5.0, g=9.81)  # ρ₁ = ρ₂ → (11.20)
    assert np.allclose(cp, ch11.vortex_sheet_c(1.0, 3.0)[0], atol=1e-13)
    V = 7.3  # Galilean: adding V to both streams shifts c by V
    a = np.array(ch11.kh_phase_speed(k, 2.0, -1.0, 1.0, 3.0, g=9.81))
    b = np.array(ch11.kh_phase_speed(k, 2.0 + V, -1.0 + V, 1.0, 3.0, g=9.81))
    assert np.allclose(b - a, V, atol=1e-12)
    # sign matters: moving the upper vs the lower stream flips the sign of the drift, never the growth
    up = ch11.kh_phase_speed(300.0, 5.0, 0.0, 1.0, 1.0, g=9.81)[0]
    lo = ch11.kh_phase_speed(300.0, 0.0, 5.0, 1.0, 1.0, g=9.81)[0]
    assert up.real == pytest.approx(2.5) and lo.real == pytest.approx(2.5)
    up2 = ch11.kh_phase_speed(300.0, 5.0, 0.0, 1.0, 3.0, g=9.81)[0]
    lo2 = ch11.kh_phase_speed(300.0, -5.0, 0.0, 1.0, 3.0, g=9.81)[0]
    assert up2.real > 0 > lo2.real and up2.imag == pytest.approx(lo2.imag, rel=1e-14)
    vs = ch11.vortex_sheet_c(3.0, 1.0)
    assert vs[0] == 2 + 1j and vs[1] == 2 - 1j  # growing root first whichever stream is faster
    g = ch11.kh_growth_rate(k, 1.0, 3.0, 1.0, 1.0, g=9.81)
    assert np.allclose(g, k * 1.0)  # vortex sheet: kΔU/2 at every k
    with pytest.raises(ValueError):
        ch11.kh_phase_speed(0.0, 1, 0, 1, 2)


def test_kh_critical_k_V1_neutral_at_kc_unstable_just_above():  # V1 (C02, N16, N18)
    for U1, U2, r1, r2 in ((5.0, 0.0, 1.2, 1000.0), (0.0, -3.0, 1.0, 3.0), (-2.0, 4.0, 900.0, 1000.0)):
        kc = ch11.kh_critical_k(U1, U2, r1, r2, g=9.81)
        assert abs(kc - 9.81 * (r2 ** 2 - r1 ** 2) / (r1 * r2 * (U2 - U1) ** 2)) < 1e-12 * kc
        assert ch11.kh_discriminant_terms(kc, U1, U2, r1, r2, g=9.81)["total"] == pytest.approx(0.0, abs=1e-12 * 9.81 / kc)
        assert ch11.kh_growth_rate(kc * (1 - 1e-6), U1, U2, r1, r2, g=9.81) == 0.0
        assert ch11.kh_growth_rate(kc * (1 + 1e-6), U1, U2, r1, r2, g=9.81) > 0.0
    assert ch11.kh_critical_k(5.0, 0.0, 1.2, 1000.0) == pytest.approx(326.89, abs=0.01)  # design: 327 m⁻¹ (λ_c = 1.92 cm)
    assert 2 * PI / ch11.kh_critical_k(5.0, 0.0, 1.2, 1000.0) == pytest.approx(0.0192, abs=5e-5)
    assert ch11.kh_critical_k(1.0, 1.0, 1.0, 2.0) == math.inf and ch11.kh_critical_k(1.0, 0.0, 2.0, 1.0) == 0.0


def test_kh_discriminant_terms_V1_tiny_example_and_sum():  # V1 (C02, E2 term bars)
    t = ch11.kh_discriminant_terms(1.0, 6.0, 0.0, 1.0, 3.0, g=10.0)
    assert (t["gravity"], t["tension"], t["shear"], t["total"], t["mean"]) == pytest.approx((5.0, 0.0, -6.75, -1.75, 1.5))
    t2 = ch11.kh_discriminant_terms(300.0, 6.0, 0.0, 1.2, 1000.0, g=9.81, surface_tension=0.074, h=0.01)
    assert t2["total"] == pytest.approx(t2["gravity"] + t2["tension"] + t2["shear"], rel=1e-14)
    cp, _ = ch11.kh_phase_speed(300.0, 6.0, 0.0, 1.2, 1000.0, g=9.81, surface_tension=0.074, h=0.01)
    assert abs(cp - (t2["mean"] + cmath.sqrt(t2["total"]))) < 1e-12


def test_kh_depth_and_tension_V1_boundary_min_shear_band_and_RT_cutoff():  # V1 (N120, N121, D04)
    r1, r2, s = 1.2, 1000.0, 0.074
    k = np.geomspace(10, 3000, 300)
    for h in (None, 0.01, 0.1):
        dU = ch11.kh_stability_boundary(k, r1, r2, g=9.81, surface_tension=s, h=h)
        tot = ch11.kh_discriminant_terms(k, dU, 0.0, r1, r2, g=9.81, surface_tension=s, h=h)["total"]
        assert np.max(np.abs(tot) / (9.81 / k)) < 1e-11  # zero discriminant on the boundary
        assert np.all(ch11.kh_growth_rate(k, 1.001 * dU, 0.0, r1, r2, g=9.81, surface_tension=s, h=h) > 0)
        assert np.all(ch11.kh_growth_rate(k, 0.999 * dU, 0.0, r1, r2, g=9.81, surface_tension=s, h=h) == 0)
    ms = ch11.kh_min_shear(r1, r2, g=9.81, surface_tension=s)
    r = minimize_scalar(lambda kk: float(ch11.kh_stability_boundary(kk, r1, r2, g=9.81, surface_tension=s)),
                        bounds=(50, 2000), method="bounded", options=dict(xatol=1e-8))
    assert ms["dU_min"] == pytest.approx(r.fun, rel=1e-10) and ms["k_star"] == pytest.approx(r.x, rel=1e-5)
    m0 = ch11.kh_min_shear(r1, r2)  # default g = G0
    assert m0["dU_min"] == pytest.approx(6.70, abs=0.005) and m0["wavelength"] == pytest.approx(0.0173, abs=5e-5)
    k1, k2 = ch11.kh_unstable_band(8.0, r1, r2, surface_tension=s)
    assert (k1, k2) == pytest.approx((149.2, 887.4), abs=0.1)  # design C.2 2.9
    for kk in (k1, k2):
        assert ch11.kh_discriminant_terms(kk, 8.0, 0.0, r1, r2, surface_tension=s)["total"] == pytest.approx(0, abs=1e-12)
    assert np.isnan(ch11.kh_unstable_band(5.0, r1, r2, surface_tension=s)[0])
    assert ch11.kh_unstable_band(5.0, r1, r2)[0] == pytest.approx(ch11.kh_critical_k(5.0, 0.0, r1, r2))
    lam = ch11.rayleigh_taylor_cutoff(s, r2, r1, g=9.81)  # water over air
    kc = 2 * PI / lam
    assert ch11.kh_growth_rate(kc * 1.001, 0, 0, r2, r1, g=9.81, surface_tension=s) == 0.0  # shorter: stable
    assert ch11.kh_growth_rate(kc * 0.999, 0, 0, r2, r1, g=9.81, surface_tension=s) > 0.0  # longer: Rayleigh–Taylor
    assert ch11.rayleigh_taylor_cutoff(s, r2, r1) == pytest.approx(0.0173, abs=5e-5)
    with pytest.raises(ValueError):
        ch11.kh_min_shear(r1, r2, surface_tension=0.0)
    with pytest.raises(ValueError):
        ch11.rayleigh_taylor_cutoff(s, 1.0, 2.0)


def test_kh_ex11_1_V7_deep_limit_returns_11_18_plus_tension():  # V7 (N120)
    k = np.geomspace(1, 1000, 50)
    a = np.array(ch11.kh_phase_speed(k, 3.0, -1.0, 1.0, 2.0, g=9.81, surface_tension=0.05, h=200.0 / k[0]))
    t = ch11.kh_discriminant_terms(k, 3.0, -1.0, 1.0, 2.0, g=9.81)
    b = np.array([t["mean"] + s * np.lib.scimath.sqrt(t["total"] + 0.05 * k / 3.0) for s in (1, -1)])
    assert np.max(np.abs(a - b)) < 1e-10


def test_kh_fields_V1_residuals_and_mode_structure():  # V1 (N09, N12, N14, R08)
    for args in ((2.0, 3.0, -1.0, 1.0, 3.0), (50.0, 6.0, 0.0, 1.2, 1000.0), (1.0, 0.5, 0.0, 1.0, 3.0)):
        r = ch11.kh_residuals(*args, g=9.81)
        assert max(r.values()) < 1e-12, r
    k, U1, U2, r1, r2 = 2.0, 3.0, -1.0, 1.0, 3.0
    c = ch11.kh_phase_speed(k, U1, U2, r1, r2, g=9.81)[0]
    Am, Ap = ch11.kh_amplitudes(k, c, U1, U2, 1.0)
    assert Am == -1j * (U1 - c) and Ap == 1j * (U2 - c)  # (11.16) solved
    assert abs(-1j * U1 * k - k * Am - (-1j * k * c)) < 1e-12 and abs(-1j * U2 * k + k * Ap - (-1j * k * c)) < 1e-12
    x = np.linspace(0, PI, 201)
    z = 0.3
    h = 1e-5
    f = ch11.kh_fields(x, z, 0.4, k, U1, U2, r1, r2, g=9.81, zeta0=0.01)
    fx = ch11.kh_fields(x + h, z, 0.4, k, U1, U2, r1, r2, g=9.81, zeta0=0.01)
    fz = ch11.kh_fields(x, z + h, 0.4, k, U1, U2, r1, r2, g=9.81, zeta0=0.01)
    assert np.allclose((fx["phi1"] - f["phi1"]) / h, f["u1"], atol=1e-6)
    assert np.allclose((fz["phi1"] - f["phi1"]) / h, f["w1"], atol=1e-6)
    f0 = ch11.kh_fields(0.0, 0.0, 0.0, k, U1, U2, r1, r2, g=9.81, zeta0=0.01)
    f1 = ch11.kh_fields(0.0, 0.0, 1.0, k, U1, U2, r1, r2, g=9.81, zeta0=0.01)
    amp = lambda tt: abs(0.01 * np.exp(-1j * k * c * tt))  # noqa: E731
    assert amp(1.0) / amp(0.0) == pytest.approx(math.exp(k * c.imag))  # the growing root grows
    assert abs(f0["zeta"]) <= 0.01 + 1e-15 and f1["c"] == c
    assert abs(ch11.kh_fields(0.0, 2.0, 0.0, k, U1, U2, r1, r2, g=9.81)["phi1"]) < abs(
        ch11.kh_fields(0.0, 0.1, 0.0, k, U1, U2, r1, r2, g=9.81)["phi1"]) + 1e-15  # decays upward


def test_kh_mixing_energy_V1_two_thirds_and_V4_momentum_and_smoothing():  # V1 + V4 (N22, N23)
    m = ch11.kh_mixing_energy(3.0, 0.5, 1000.0)
    assert m["ratio"] == pytest.approx(2 / 3, rel=1e-14) and m["M_i"] == m["M_f"]
    assert m["E_f_quad"] == pytest.approx(m["E_f"], rel=1e-12)
    rng = np.random.default_rng(2)
    for _ in range(10):  # any momentum-conserving smoothing of the step lowers the energy (Cauchy–Schwarz)
        w = rng.uniform(0.05, 0.5)
        prof = lambda z, w=w: 3.0 * 0.5 * (1 + math.erf(z / w))  # noqa: E731
        r = ch11.kh_mixing_energy(3.0, 50.0, 1000.0, profile=prof)
        assert r["E_f"] < r["E_i"] and r["M_f"] == pytest.approx(r["M_i"], rel=1e-9)


def test_kh_V2_derivation_D02_linearised_interface_conditions():  # V2 (D02 ★★, steps 1–12)
    x, z, t, eps = sp.symbols("x z t epsilon", real=True)
    U1, U2, g, r1, r2, C1, C2 = sp.symbols("U1 U2 g rho1 rho2 C1 C2", real=True)
    Z = sp.Function("Z")(x, t)
    p1, p2 = sp.Function("p1")(x, z, t), sp.Function("p2")(x, z, t)
    phit1, phit2 = U1 * x + eps * p1, U2 * x + eps * p2  # (11.2)
    zeta = eps * Z
    # steps 1–5: unit normal of f = z − ζ, three dot products, the square root cancels
    zx = sp.diff(zeta, x)
    root = sp.sqrt(1 + zx ** 2)
    n = sp.Matrix([-zx, 1]) / root
    Us = sp.Matrix([0, sp.diff(zeta, t)])
    for pt in (phit1, phit2):
        grad = sp.Matrix([sp.diff(pt, x), sp.diff(pt, z)])
        lhs = sp.simplify((n.dot(grad) - n.dot(Us)) * root)
        assert z0(lhs - (-sp.diff(pt, x) * zx + sp.diff(pt, z) - sp.diff(zeta, t)))  # step 5
    # steps 6–7: O(ε) on z = ζ (Taylor transfer to z = 0)
    kin = (-sp.diff(phit1, x) * zx + sp.diff(phit1, z) - sp.diff(zeta, t)).subs(z, zeta)
    kin1 = sp.diff(kin, eps).subs(eps, 0).doit()
    target = (-U1 * sp.diff(Z, x) + sp.diff(p1, z) - sp.diff(Z, t)).subs(z, 0).doit()
    assert z0(kin1 - target)  # (11.9)
    assert sp.diff(kin, eps, 2).subs(eps, 0).doit() != 0  # the dropped terms are second order, not zero
    # steps 8–12: Bernoulli pressures equal on z = ζ, minus the undisturbed balance (11.12)
    P = lambda rho, C, pt: rho * (C - sp.diff(pt, t) - (sp.diff(pt, x) ** 2 + sp.diff(pt, z) ** 2) / 2 - g * z)  # noqa: E731
    jump = (P(r1, C1, phit1) - P(r2, C2, phit2)).subs(z, zeta)
    base = jump.subs(eps, 0)
    assert z0(base - (r1 * (C1 - U1 ** 2 / 2) - r2 * (C2 - U2 ** 2 / 2)))  # (11.12) is the O(1) balance
    dyn1 = sp.diff(jump, eps).subs(eps, 0).doit()
    target = -(r1 * (sp.diff(p1, t) + U1 * sp.diff(p1, x) + g * Z) - r2 * (sp.diff(p2, t) + U2 * sp.diff(p2, x) + g * Z))
    assert z0(dyn1 - target.subs(z, 0).doit())  # (11.13): U∂φ/∂x survives, not ½(∂φ/∂x)²


def test_kh_V2_derivation_D03_dispersion_relation_from_normal_modes():  # V2 (D03 ★★, steps 1–14)
    z, k, g, r1, r2, z_0 = sp.symbols("z k g rho1 rho2 zeta0", positive=True)
    U1, U2, c, Am, Ap = sp.symbols("U1 U2 c A_m A_p")
    A = sp.Function("A")
    sol = sp.dsolve(sp.Eq(A(z).diff(z, 2) - k ** 2 * A(z), 0))  # steps 1–2
    assert {sp.exp(k * z), sp.exp(-k * z)} <= set(sol.rhs.atoms(sp.exp))
    f1, f2 = Am * sp.exp(-k * z), Ap * sp.exp(k * z)  # step 3 (11.15): decaying on each side
    kin = [sp.Eq(-sp.I * U1 * k * z_0 + sp.diff(f1, z).subs(z, 0), -sp.I * k * c * z_0),
           sp.Eq(-sp.I * U2 * k * z_0 + sp.diff(f2, z).subs(z, 0), -sp.I * k * c * z_0)]  # step 5 (11.16)
    amp = sp.solve(kin, [Am, Ap], dict=True)[0]
    assert z0(amp[Am] + sp.I * (U1 - c) * z_0) and z0(amp[Ap] - sp.I * (U2 - c) * z_0)  # step 6
    dyn = (r1 * (-sp.I * k * c * Am + sp.I * k * U1 * Am + g * z_0) - r2 * (-sp.I * k * c * Ap + sp.I * k * U2 * Ap + g * z_0))
    step8 = r1 * (k * (U1 - c) ** 2 + g) * z_0 - r2 * (-k * (U2 - c) ** 2 + g) * z_0
    assert z0(dyn.subs(amp) - step8)  # step 8
    quad_ = sp.expand(step8 / (k * z_0))
    assert z0(quad_ - (r1 * (U1 - c) ** 2 + r2 * (U2 - c) ** 2 - g / k * (r2 - r1)))  # step 9
    assert z0((r1 * U1 + r2 * U2) ** 2 - (r1 + r2) * (r1 * U1 ** 2 + r2 * U2 ** 2) + r1 * r2 * (U1 - U2) ** 2)  # step 12
    disc = (r2 - r1) / (r2 + r1) * g / k - r1 * r2 * (U2 - U1) ** 2 / (r1 + r2) ** 2
    for s in (1, -1):
        root = (r1 * U1 + r2 * U2) / (r1 + r2) + s * sp.sqrt(disc)  # (11.18)
        assert z0(quad_.subs(c, root))
    assert z0(k * (r1 + r2) ** 2 * disc - (g * (r2 ** 2 - r1 ** 2) - k * r1 * r2 * (U2 - U1) ** 2))  # step 13
    # the coded root equals the derived root at random numbers (independent of kh_sympy)
    rng = np.random.default_rng(3)
    for _ in range(20):
        vals = dict(zip((k, g, r1, r2), rng.uniform(0.1, 5, 4)))
        vals.update({U1: rng.uniform(-3, 3), U2: rng.uniform(-3, 3)})
        der = complex(sp.N(((r1 * U1 + r2 * U2) / (r1 + r2)).subs(vals))) + cmath.sqrt(complex(sp.N(disc.subs(vals))))
        code = ch11.kh_phase_speed(float(vals[k]), float(vals[U1]), float(vals[U2]), float(vals[r1]), float(vals[r2]),
                                   g=float(vals[g]))[0]
        assert abs(code - der) < 1e-12 * max(1, abs(der))


def test_kh_V2_derivation_D04_depth_and_surface_tension():  # V2 (D04 ★★, steps 1–9)
    z, k, g, r1, r2, sg, h, z_0 = sp.symbols("z k g rho1 rho2 sigma_s h zeta0", positive=True)
    U1, U2, c, Am, Ap = sp.symbols("U1 U2 c A_m A_p")
    f2 = Ap * sp.cosh(k * (z + h))  # step 1
    assert z0(sp.diff(f2, z, 2) - k ** 2 * f2) and z0(sp.diff(f2, z).subs(z, -h))
    sol = sp.solve([sp.Eq(-sp.I * U1 * k * z_0 - k * Am, -sp.I * k * c * z_0),
                    sp.Eq(-sp.I * U2 * k * z_0 + sp.diff(f2, z).subs(z, 0), -sp.I * k * c * z_0)], [Am, Ap], dict=True)[0]
    assert z0(sol[Ap] - sp.I * (U2 - c) * z_0 / sp.sinh(k * h))  # step 3
    phi2_0 = f2.subs(z, 0).subs(sol)
    p1 = -r1 * (-sp.I * k * c * sol[Am] + sp.I * k * U1 * sol[Am] + g * z_0)
    p2 = -r2 * (-sp.I * k * c * phi2_0 + sp.I * k * U2 * phi2_0 + g * z_0)
    cond = sp.expand((p1 - p2 + sg * k ** 2 * z_0) / (k * z_0))  # step 4–5: p₁ − p₂ = σ_sζ_xx = −σ_sk²ζ
    step6 = r1 * (U1 - c) ** 2 + r2 * sp.cosh(k * h) / sp.sinh(k * h) * (U2 - c) ** 2 - (g / k * (r2 - r1) + sg * k)
    assert z0((cond + step6).rewrite(sp.exp))  # step 6 (up to the overall sign)
    rng = np.random.default_rng(4)
    for _ in range(15):  # step 7: the coded Ex. 11.1 root solves step 6
        kk, hh, ss, a1, a2 = rng.uniform(0.2, 5), rng.uniform(0.05, 3), rng.uniform(0, 2), rng.uniform(0.5, 2), rng.uniform(2, 4)
        u1, u2 = rng.uniform(-3, 3, 2)
        for cc in ch11.kh_phase_speed(kk, u1, u2, a1, a2, g=9.81, surface_tension=ss, h=hh):
            res = a1 * (u1 - cc) ** 2 + a2 / math.tanh(kk * hh) * (u2 - cc) ** 2 - (9.81 / kk * (a2 - a1) + ss * kk)
            assert abs(res) < 1e-11 * (1 + abs(a1 * (u1 - cc) ** 2))
    kk = sp.symbols("kk", positive=True)
    dr = sp.symbols("Delta_rho", positive=True)
    f = g / kk * dr + sg * kk  # step 8
    ks = sp.solve(sp.diff(f, kk), kk)
    assert len(ks) == 1 and z0(ks[0] - sp.sqrt(g * dr / sg)) and z0(f.subs(kk, ks[0]) - 2 * sp.sqrt(g * dr * sg))
    dU2 = sp.simplify(f.subs(kk, ks[0]) * (r1 + r2) / (r1 * r2))  # step 9
    assert z0(dU2 - 2 * sp.sqrt(g * dr * sg) * (r1 + r2) / (r1 * r2))
    assert sp.limit(sp.cosh(k * h) / sp.sinh(k * h), h, sp.oo) == 1


def test_kh_engines_V2_sympy_residuals_are_zero():  # V2 (engine check: kh_sympy, kh_depth_tension_sympy)
    r = ch11.kh_sympy()
    for key in ("kinematic_lin", "dynamic_lin", "identity", "criterion"):
        assert z0(r[key]), key
    assert all(z0(v) for v in r["residual_11_18"] + r["vortex_sheet"] + r["static"])
    d = ch11.kh_depth_tension_sympy()
    assert all(z0(v) for v in d["residual"]) and z0(d["limit_residual"]) and d["coth_limit"] == 1


# ======================================================================================================================
# C03 — the Bénard amplitude problem (11.36)–(11.37), Ra (11.21), exchange of stabilities
# ======================================================================================================================
def test_rayleigh_number_V2_dimensionless_and_V7_unit_invariance_and_sign():  # V2 + V7 (N25, slip S10)
    q = (Q_(9.8, "m/s**2") * Q_(2.1e-4, "1/K") * Q_(2.0, "K") * Q_(5, "mm") ** 3 / (Q_(1.4e-7, "m**2/s") * Q_(1e-6, "m**2/s")))
    assert q.dimensionless and q.to("dimensionless").magnitude == pytest.approx(9.8 * 2.1e-4 * 2 * 125e-9 / 1.4e-13)  # (11.21)
    dimensional_check(lambda alpha, dT, d: alpha * dT / d, "[length] ** -1", alpha=Q_(2.1e-4, "1/K"), dT=Q_(2.0, "K"),
                      d=Q_(5, "mm"))  # αΓ (Γ = ΔT/d) is an inverse length
    base = ch11.rayleigh_number(2.1e-4, 2.0, 5e-3, 1.4e-7, 1e-6, g=9.81)
    assert base == pytest.approx(9.81 * 2.1e-4 * 2.0 * 5e-3 ** 3 / (1.4e-7 * 1e-6), rel=1e-14)
    lam, tau, th = 1000.0, 60.0, 1.8  # lengths in mm, time in minutes, temperature in °F-sized units: Ra must not change
    alt = ch11.rayleigh_number(2.1e-4 / th, 2.0 * th, 5e-3 * lam, 1.4e-7 * lam ** 2 * tau, 1e-6 * lam ** 2 * tau,
                               g=9.81 * lam * tau ** 2)
    assert alt == pytest.approx(base, rel=1e-13)
    assert ch11.rayleigh_number(2.1e-4, 2.0, 5e-3, 1.4e-7, 1e-6) == pytest.approx(3677.49, abs=0.01)  # G0 (design: 3677)
    assert base == pytest.approx(3678.75, abs=0.01)  # g = 9.81 (design: 3679)
    assert ch11.rayleigh_number(2.1e-4, -2.0, 5e-3, 1.4e-7, 1e-6) < 0  # heated from above: Ra < 0 (sign matters)
    gc = ch11.gamma_conventions(2.0, 5e-3)
    assert (gc["Gamma_11_21"], gc["dTdz_kundu"], gc["Gamma_met"]) == pytest.approx((400.0, -400.0, 400.0))
    T = ch11.benard_base_state(np.array([-2.5e-3, 2.5e-3]), 300.0, 2.0, 5e-3)["T"]
    assert (T[1] - T[0]) / 5e-3 == pytest.approx(gc["dTdz_kundu"])  # Kundu's Γ = dT/dz of the actual profile
    with pytest.raises(ValueError):
        ch11.rayleigh_number(2.1e-4, 2.0, 0.0, 1.4e-7, 1e-6)


def test_benard_base_state_V1_hydrostatic_conduction_and_scales():  # V1 (R13, N26, N28, N31)
    z = np.linspace(-2.5e-3, 2.5e-3, 41)
    T0, dT, d, rho0, al = 300.0, 2.0, 5e-3, 1000.0, 2.1e-4
    b = ch11.benard_base_state(z, T0, dT, d, rho0=rho0, alpha=al, g=9.81, P0=1e5)
    assert b["T"][0] == pytest.approx(T0) and b["T"][-1] == pytest.approx(T0 - dT)  # bottom T₀, top T₀ − ΔT
    assert np.allclose(np.diff(b["T"], 2), 0, atol=1e-12)  # ∂²T̄/∂z² = 0
    h = 1e-7
    zz = z[5:-5]
    dP = (ch11.benard_base_state(zz + h, T0, dT, d, rho0, al, 9.81)["P"]
          - ch11.benard_base_state(zz - h, T0, dT, d, rho0, al, 9.81)["P"]) / (2 * h)
    Tz = ch11.benard_base_state(zz, T0, dT, d, rho0, al, 9.81)["T"]
    assert np.allclose(dP, -rho0 * 9.81 * (1 - al * (Tz - T0)), rtol=1e-7)  # (11.23)
    s = ch11.benard_scales(2.1e-4, 2.0, 5e-3, 1.4e-7, 1e-6, g=9.81)
    assert s["Pr"] == pytest.approx(1e-6 / 1.4e-7) and s["w_scale"] == pytest.approx(2.8e-5)
    assert s["t_scale"] == pytest.approx(178.571, rel=1e-5) and s["Gamma"] == pytest.approx(400.0)


def test_benard_growth_rate_V1_free_free_equals_the_D12_quadratic():  # V1 (C03, C05, D12)
    for K, Ra, Pr in ((PI / math.sqrt(2), 2000.0, 1.0), (1.5, 900.0, 7.0), (3.0, 5000.0, 0.1), (2.0, 300.0, 1.0)):
        sp_, sm = ch11.benard_free_free_sigma(K, Ra, Pr)
        assert ch11.benard_growth_rate(K, Ra, Pr, bc=("free", "free"), N=40) == pytest.approx(sp_, abs=1e-8)
        allc = ch11.benard_growth_rate(K, Ra, Pr, bc=("free", "free"), N=40, all=True)
        assert np.min(np.abs(allc - sm)) < 1e-6 * max(1, abs(sm))
    a = ch11.benard_free_free_sigma(PI / math.sqrt(2), 2000.0, 1.0)
    assert a == pytest.approx((11.0155, -40.6243), abs=1e-4)
    assert ch11.benard_free_free_sigma(PI / math.sqrt(2), ch11.RA_FREE_FREE, 1.0)[0] == pytest.approx(0, abs=1e-10)


def test_benard_growth_rate_V7_real_spectrum_and_sign_change_at_the_margin():  # V7 (C03, N40 exchange of stabilities)
    for K, Ra, Pr in ((3.0, 2500.0, 0.7), (2.0, 1000.0, 7.0), (4.0, 10000.0, 0.025)):
        s = ch11.benard_growth_rate(K, Ra, Pr, all=True, return_complex=True)
        assert len(s) > 10 and np.max(np.abs(s.imag)) < 1e-10 * max(1.0, np.max(np.abs(s)))
    for K in (2.0, 3.1163, 5.0):
        Ram = ch11.benard_marginal_Ra(K)
        assert abs(ch11.benard_growth_rate(K, Ram, 1.0)) < 1e-6
        assert ch11.benard_growth_rate(K, 1.01 * Ram, 1.0) > 0 > ch11.benard_growth_rate(K, 0.99 * Ram, 1.0)
    # Pr changes the rate, not the margin
    Ram = ch11.benard_marginal_Ra(3.0)
    assert all(abs(ch11.benard_growth_rate(3.0, Ram, Pr)) < 1e-6 for Pr in (0.01, 1.0, 100.0))
    # heated from above: every mode decays
    assert ch11.benard_growth_rate(3.0, -3000.0, 1.0) < 0


def test_benard_growth_rate_V3_spectral_convergence_in_N():  # V3 (C03)
    ref = ch11.benard_growth_rate(3.0, 2500.0, 0.7, N=30, filter=False)
    errs = [abs(ch11.benard_growth_rate(3.0, 2500.0, 0.7, N=N, filter=False) - ref) for N in (12, 16, 20, 24)]
    assert errs[0] > 1e-7 and max(errs[1:]) < 1e-9, errs  # 1.2e-6 at N = 12, round-off floor ~1e-12 from N = 16
    assert abs(ch11.benard_growth_rate(3.0, 2500.0, 0.7) - ref) < 1e-8  # the default N = 40 (filtered) sits on the floor


def test_benard_V2_derivation_D06_pressure_elimination():  # V2 (D06 ★★, steps 1–7; planted ∇² fails)
    x, y, z, t = sp.symbols("x y z t", real=True)
    rho0, g, al, nu = sp.symbols("rho0 g alpha nu", positive=True)
    A1, A2, A3, p, T = [sp.Function(n)(x, y, z, t) for n in ("A1", "A2", "A3", "p", "T")]
    u = sp.diff(A3, y) - sp.diff(A2, z)  # a divergence-free field from a vector potential (11.25 by construction)
    v = sp.diff(A1, z) - sp.diff(A3, x)
    w = sp.diff(A2, x) - sp.diff(A1, y)
    assert z0(sp.diff(u, x) + sp.diff(v, y) + sp.diff(w, z))
    lap = lambda f: sp.diff(f, x, 2) + sp.diff(f, y, 2) + sp.diff(f, z, 2)  # noqa: E731
    R = [sp.diff(q, t) + sp.diff(p, s) / rho0 - nu * lap(q) for q, s in ((u, x), (v, y), (w, z))]  # (11.26) residuals
    R[2] -= g * al * T
    divR = sp.diff(R[0], x) + sp.diff(R[1], y) + sp.diff(R[2], z)
    assert z0(divR - (lap(p) / rho0 - g * al * sp.diff(T, z)))  # steps 2–3: only the Poisson part survives
    e29 = sp.diff(lap(w), t) - g * al * (sp.diff(T, x, 2) + sp.diff(T, y, 2)) - nu * lap(lap(w))  # (11.29)
    assert z0(lap(R[2]) - sp.diff(divR, z) - e29)  # steps 1, 4–7
    planted = sp.diff(lap(w), t) - g * al * lap(T) - nu * lap(lap(w))
    assert not z0(lap(R[2]) - sp.diff(divR, z) - planted)


def test_benard_V2_derivation_D08_exchange_of_stabilities():  # V2 (D08 ★★, steps 1–12) + numeric check on computed modes
    z = sp.symbols("z", real=True)
    K, Pr, Ra = sp.symbols("K Pr Ra", positive=True)
    rng = np.random.default_rng(5)
    cT = [sp.Rational(int(a), 7) + sp.I * sp.Rational(int(b), 5) for a, b in rng.integers(-5, 6, (3, 2))]
    cW = [sp.Rational(int(a), 3) + sp.I * sp.Rational(int(b), 4) for a, b in rng.integers(-5, 6, (3, 2))]
    T = (z ** 2 - sp.Rational(1, 4)) * sum(cc * z ** j for j, cc in enumerate(cT))  # T̂ = 0 at the walls
    W = (z ** 2 - sp.Rational(1, 4)) ** 2 * sum(cc * z ** j for j, cc in enumerate(cW))  # W = W′ = 0 at the walls
    I = lambda f: sp.integrate(sp.expand(f), (z, -sp.Rational(1, 2), sp.Rational(1, 2)))  # noqa: E731
    cj = sp.conjugate
    assert z0(-I(cj(T) * sp.diff(T, z, 2)) - I(cj(sp.diff(T, z)) * sp.diff(T, z)))  # step 2
    L = lambda f: sp.diff(f, z, 2) - K ** 2 * f  # noqa: E731
    J1 = I(cj(sp.diff(W, z)) * sp.diff(W, z) + K ** 2 * cj(W) * W)
    assert z0(-I(cj(W) * L(W)) - J1)  # step 6
    assert z0(I(cj(W) * L(L(W))) - I(cj(L(W)) * L(W)))  # steps 7–8: J₂ = ∫|(D² − K²)W|²
    s_r, s_i = sp.symbols("sigma_r sigma_i", real=True)
    i1, i2, j1, j2 = sp.symbols("I1 I2 J1 J2", positive=True)
    s = s_r + sp.I * s_i
    step10 = s / Pr * j1 + j2 - Ra * K ** 2 * (sp.conjugate(s) * i1 + i2)
    assert z0(sp.im(sp.expand(step10)) - s_i * (j1 / Pr + Ra * K ** 2 * i1))  # step 11
    assert sp.solve(sp.Eq(sp.im(sp.expand(step10)), 0), s_i) == [0]  # step 12 (Ra > 0)
    # the two energy relations hold for actual computed eigenmodes of (11.36)–(11.38)
    K_, Ra_, Pr_ = 3.0, 2500.0, 0.7
    D, D2, Lm, Im, zz, n = ch11._benard_ops(K_, 40)
    A = np.block([[Lm @ Lm, -Ra_ * K_ ** 2 * Im], [Im, Lm]])
    B = np.block([[Lm / Pr_, np.zeros_like(Im)], [np.zeros_like(Im), Im]])
    C, elim = ch11._benard_constraints(("rigid", "rigid"), D, D2, Im, n, 40)
    w, V = ST.constrained_eig(A, B, C, elim, return_vectors=True, sort="real")
    wq = ST.clenshaw_curtis_weights(40, (-0.5, 0.5))
    for j in range(3):
        sg, Wv, Tv = w[j], V[:n, j], V[n:, j]
        Iq = lambda f: np.sum(wq * f)  # noqa: E731
        I1, I2 = Iq(abs(Tv) ** 2), Iq(abs(D @ Tv) ** 2 + K_ ** 2 * abs(Tv) ** 2)
        J1n, J2n = Iq(abs(D @ Wv) ** 2 + K_ ** 2 * abs(Wv) ** 2), Iq(abs(Lm @ Wv) ** 2)
        lhs, rhs = sg / Pr_ * J1n + J2n, Ra_ * K_ ** 2 * (np.conj(sg) * I1 + I2)
        assert abs(lhs - rhs) < 1e-6 * abs(lhs)  # step 10 on a computed mode
        assert abs(sg.imag) < 1e-10 * max(1, abs(sg))


def test_benard_engines_V2_perturbation_and_exchange_sympy():  # V2 (engines: benard_perturbation_sympy, exchange…)
    r = ch11.benard_perturbation_sympy()
    for key in ("eq_11_27", "eq_11_29_residual", "eq_11_32", "eq_11_36", "eq_11_37", "eq_11_40_residual"):
        assert z0(r[key]), key
    assert not z0(r["planted_residual"])
    e = ch11.exchange_of_stabilities_sympy()
    assert z0(e["relation_T"]) and z0(e["relation_W"]) and z0(e["imag_identity"]) and e["sigma_i_solutions"] == [0]


# ======================================================================================================================
# C04 — rigid–rigid neutral curve and Ra_c
# ======================================================================================================================
def test_benard_critical_V5_chandrasekhar_rigid_rigid_and_rigid_free():  # V5 (C04, N50; Nek5000/arXiv tables)
    rr = ch11.benard_critical()
    assert rel(rr["Ra_c"], 1707.762) < 1e-5 and rel(rr["K_c"], 3.117) < 1e-3  # source rounds 3.1163 to 3.117
    rf = ch11.benard_critical(bc=("rigid", "free"), mode="any")
    assert rel(rf["Ra_c"], 1100.65) < 1e-5 and rel(rf["K_c"], 2.682) < 1e-3
    ff = ch11.benard_critical(bc=("free", "free"))
    assert rel(ff["Ra_c"], 27 * PI ** 4 / 4) < 1e-8 and abs(ff["K_c"] - PI / math.sqrt(2)) < 1e-4  # V1 analytic
    assert 2 * PI / rr["K_c"] == pytest.approx(2.016, abs=1e-3)  # cells about twice as wide as the layer is deep


def test_benard_marginal_V1_chebyshev_equals_the_determinant_route():  # V1 independent route (C04, D10, N123)
    for K in (2.0, 3.1163, 5.0, 7.5):
        a, b = ch11.benard_marginal_Ra(K), ch11.benard_marginal_Ra_det(K)
        assert rel(a, b) < 1e-8, (K, a, b)
    assert ch11.benard_marginal_Ra(2.0) == pytest.approx(2177.41, abs=0.01)
    assert ch11.benard_marginal_Ra(5.0) == pytest.approx(2439.32, abs=0.01)
    for K in (4.5, 5.3647, 7.0):  # odd mode (Ex. 11.7): sin/sinh determinant vs Chebyshev parity filter
        assert rel(ch11.benard_marginal_Ra(K, mode="odd"), ch11.benard_marginal_Ra_det(K, mode="odd")) < 1e-8
    odd = ch11.benard_critical(mode="odd")
    assert rel(odd["Ra_c"], 17610.39) < 1e-6 and abs(odd["K_c"] - 5.365) < 6e-4  # Chandrasekhar (1961) p. 39 (not fetched)
    with pytest.raises(ValueError):
        ch11.benard_marginal_Ra(3.0, bc=("rigid", "free"), mode="even")


def test_benard_guards_V7_unresolved_or_unbracketed_cases_raise_value_error():  # V7 (loop 2: no IndexError, no NaN)
    with pytest.raises(ValueError, match="filter"):  # nothing survives the N-filter at N = 8 (was a bare IndexError)
        ch11.benard_growth_rate(3.0, 2500.0, 0.7, N=8)
    ref = ch11.benard_growth_rate(3.0, 2500.0, 0.7)
    assert ch11.benard_growth_rate(3.0, 2500.0, 0.7, N=8, filter=False) == pytest.approx(ref, rel=1e-3)  # the way out works
    assert ch11.benard_growth_rate(3.0, 2500.0, 0.7, N=12) == pytest.approx(ref, abs=1e-5)  # N = 12 already resolves it
    assert len(ch11.benard_growth_rate(3.0, 2500.0, 0.7, N=8, all=True)) == 0  # all=True hands back the (empty) list
    for K, Pr in ((0.0, 1.0), (-1.0, 1.0), (3.0, 0.0)):
        with pytest.raises(ValueError):
            ch11.benard_growth_rate(K, 2500.0, Pr)
    # determinant route: the root lies above the default scan end Ra_max = 2e4 for the even mode at K ≤ 0.5 and for the odd
    # mode at K < 4.0 (was a silent NaN) — and is found, equal to the Chebyshev value, once Ra_max is raised
    for K, mode in ((0.5, "even"), (1.0, "odd"), (3.0, "odd")):
        with pytest.raises(ValueError, match="Ra_max"):
            ch11.benard_marginal_Ra_det(K, mode=mode)
        a, b = ch11.benard_marginal_Ra_det(K, mode=mode, Ra_max=1e6), ch11.benard_marginal_Ra(K, mode=mode)
        assert a > 2e4 and rel(a, b) < 1e-8, (K, mode, a, b)
    for bad in (0.0, -1.0, float("nan")):
        with pytest.raises(ValueError):
            ch11.benard_marginal_Ra_det(bad)
    with pytest.raises(ValueError):
        ch11.benard_marginal_Ra_det(3.0, mode="both")
    assert all(np.isfinite(ch11.benard_marginal_Ra_det(K)) for K in (0.6, 2.0, 9.0))  # default bracket: even K = 0.6 … 9
    assert all(np.isfinite(ch11.benard_marginal_Ra_det(K, mode="odd")) for K in (4.0, 5.3647, 9.0))  # odd K = 4 … 9


def test_benard_marginal_V3_converges_in_N():  # V3 (C04)
    ref = ch11.benard_marginal_Ra_det(3.1163)
    errs = [abs(ch11.benard_marginal_Ra(3.1163, N=N) - ref) / ref for N in (10, 14, 20, 30)]
    assert errs[0] > 1e-6 and errs[1] < 1e-8 and errs[-1] < 1e-9, errs  # 3e-6, 2e-10, then the 1e-12 round-off floor


def test_benard_V2_derivation_D10_roots_operator_and_imaginary_determinant():  # V2 (D10 ★★★, steps 1–14)
    z, K, Ra = sp.symbols("z K Ra", positive=True)
    s = (Ra / K ** 4) ** sp.Rational(1, 3)
    for q2 in (-K ** 2 * (s - 1), K ** 2 * (1 + s * (1 + sp.sqrt(3) * sp.I) / 2), K ** 2 * (1 + s * (1 - sp.sqrt(3) * sp.I) / 2)):
        assert z0((q2 - K ** 2) ** 3 + Ra * K ** 2)  # steps 1–4: (11.42) solves (q² − K²)³ = −RaK²
    q, q0 = sp.symbols("q q0")
    L = lambda f: sp.diff(f, z, 2) - K ** 2 * f  # noqa: E731
    assert z0(L(L(sp.cosh(q * z))) - (q ** 2 - K ** 2) ** 2 * sp.cosh(q * z))  # step 9
    assert z0(L(L(sp.cos(q0 * z))) - (q0 ** 2 + K ** 2) ** 2 * sp.cos(q0 * z))
    assert z0(L(L(L(sp.cosh(q * z)))) - (q ** 2 - K ** 2) ** 3 * sp.cosh(q * z))  # (11.40) on each exponential
    # step 12: the determinant is purely imaginary (conjugate columns) — and the printed slip #1 breaks it
    for Ra_, K_ in ((1000.0, 2.0), (1707.762, 3.1163), (5000.0, 5.0), (3000.0, 3.0)):
        q0n, qn, qs = ch11.benard_char_roots(Ra_, K_)
        assert abs(qs - qn.conjugate()) < 1e-15 and qn.real > 0
        assert abs((qn ** 2 - K_ ** 2) ** 3 + Ra_ * K_ ** 2) < 1e-9 * Ra_ * K_ ** 2
        assert abs((-q0n ** 2 - K_ ** 2) ** 3 + Ra_ * K_ ** 2) < 1e-9 * Ra_ * K_ ** 2  # q = iq₀
        d = ch11.benard_determinant(Ra_, K_)
        assert abs(d.real) < 1e-9 * max(abs(d), 1.0)
        dp = ch11.benard_determinant(Ra_, K_, printed=True)
        assert abs(dp.real) > 1e-3 * abs(dp)  # slip #1 matrix has no conjugate symmetry
    q0n, qn, _ = ch11.benard_char_roots(1000.0, 2.0)
    assert q0n == pytest.approx(3.4459, abs=1e-4) and abs(qn - (3.8822 + 1.7705j)) < 1e-4  # design example
    # step 13–14: Im det changes sign at the Chebyshev marginal Ra
    Ram = ch11.benard_marginal_Ra(3.1163)
    assert ch11.benard_determinant(0.999 * Ram, 3.1163).imag * ch11.benard_determinant(1.001 * Ram, 3.1163).imag < 0
    # slip #1 moves the marginal Ra: 1.6 % at K = 2, 2.4 % at K = 5 (only 0.44 % near K_c)
    for K_, lo in ((2.0, 0.01), (5.0, 0.01), (3.1163, 0.003)):
        shift = rel(ch11.benard_marginal_Ra_det(K_, printed=True), ch11.benard_marginal_Ra_det(K_))
        assert shift > lo, (K_, shift)


def test_benard_eigenfunction_V1_boundary_conditions_parity_and_roll_kinematics():  # V1 (N44, Fig. 11.9)
    for mode in ("even", "odd"):
        Kx = 3.1163 if mode == "even" else 5.3647
        Ra, m = ch11.benard_marginal_Ra(Kx, mode=mode, return_mode=True)
        D, zn = ST.cheb(40, (-0.5, 0.5))
        W, T = m["W"], m["T"]
        assert max(abs(W[0]), abs(W[-1]), abs((D @ W)[0]), abs((D @ W)[-1]), abs(T[0]), abs(T[-1])) < 1e-10
        sgn = 1 if mode == "even" else -1
        assert np.allclose(W, sgn * W[::-1], atol=1e-7) and np.max(np.abs(W)) == pytest.approx(1.0)  # observed 5e-9
        L = D @ D - Kx ** 2 * np.eye(41)
        assert np.max(np.abs((L @ T + W)[1:-1])) < 1e-7  # (11.39) first equation (interior)
    e = ch11.benard_eigenfunction(3.1163, x=np.linspace(0, 2, 201), z=np.linspace(-0.45, 0.45, 91))
    h = e["x"][0, 1] - e["x"][0, 0]
    hz = e["z"][1, 0] - e["z"][0, 0]
    dpsi_dx = np.gradient(e["psi"], h, axis=1)
    dpsi_dz = np.gradient(e["psi"], hz, axis=0)
    sl = (slice(2, -2), slice(2, -2))
    assert np.max(np.abs((dpsi_dx - e["w"])[sl])) < 1e-3 and np.max(np.abs((-dpsi_dz - e["u"])[sl])) < 1e-3  # §11.14 sign
    dux = np.gradient(e["u"], h, axis=1)
    div = dux + np.gradient(e["w"], hz, axis=0)
    assert np.max(np.abs(div[sl])) < 2e-3 * np.max(np.abs(dux))  # continuity, up to the O(h²) FD error
    assert e["Ra"] == pytest.approx(1707.762, rel=1e-5)


def test_planform_V1_helmholtz_equation():  # V1 (N51)
    x, y = np.meshgrid(np.linspace(0, 3, 61), np.linspace(0, 3, 61))
    h = 1e-4
    for kind in ("rolls", "squares", "hexagons"):
        f = ch11.planform(x, y, 3.1, kind)
        lap = (ch11.planform(x + h, y, 3.1, kind) + ch11.planform(x - h, y, 3.1, kind) + ch11.planform(x, y + h, 3.1, kind)
               + ch11.planform(x, y - h, 3.1, kind) - 4 * f) / h ** 2
        assert np.max(np.abs(lap + 3.1 ** 2 * f)) < 1e-4
    with pytest.raises(ValueError):
        ch11.planform(x, y, 3.1, "triangles")


def test_benard_tables_V1_cached_neutral_table_equals_live_and_analytic():  # V1 anti-cache (E3/F2 table)
    t = read_csv("benard_neutral_curves.csv")
    idx = [20, 50, 95]
    for i in idx:
        K = t["K"][i]
        assert rel(t["Ra_rigid_rigid"][i], ch11.benard_marginal_Ra_det(K)) < 1e-5
        assert rel(t["Ra_free_free"][i], (PI ** 2 + K ** 2) ** 3 / K ** 2) < 1e-5
        assert rel(t["Ra_rigid_free"][i], ch11.benard_marginal_Ra(K, bc=("rigid", "free"), mode="any")) < 1e-5
        assert rel(t["Ra_odd_rigid_rigid"][i], ch11.benard_marginal_Ra_det(K, mode="odd", Ra_max=1e6)) < 1e-5
    live = ch11.benard_neutral_table(Ks=[1.0, 3.0], write=False, cache=False)
    assert rel(live["rigid"][1], ch11.benard_marginal_Ra_det(3.0)) < 1e-8
    assert np.allclose(ch11.benard_neutral_curve([1.0, 3.0]), live["rigid"])
    assert np.all(t["Ra_rigid_rigid"] > t["Ra_rigid_free"]) and np.all(t["Ra_rigid_free"] > t["Ra_free_free"])


def test_critical_points_json_V1_equals_live_recomputation():  # V1 anti-cache (reference/ch11/critical_points.json)
    cp = json.loads((REF / "critical_points.json").read_text(encoding="utf-8"))
    b = cp["benard"]
    assert rel(b["rigid_rigid"]["Ra_c"], ch11.benard_critical()["Ra_c"]) < 1e-10
    assert rel(b["rigid_rigid_odd"]["Ra_c"], ch11.benard_critical(mode="odd")["Ra_c"]) < 1e-10
    assert rel(cp["taylor_mu0"]["Ta_c"], ch11.taylor_critical(0.0)["Ta_c"]) < 1e-10
    assert rel(cp["poiseuille"]["Re_c"], ch11.poiseuille_critical(N=100, cache=False)["Re_c"]) < 1e-7
    assert rel(cp["lorenz_r_H"], 10 * (10 + 8 / 3 + 3) / (10 - 8 / 3 - 1)) < 1e-14
    assert rel(cp["piecewise_layer_neutral_kh"], ch11.piecewise_neutral_kh()) < 1e-14
    assert np.allclose(cp["feigenbaum_ratios"], ch11.feigenbaum_estimate(8), rtol=1e-12)
    bl = ch11.blasius_critical()  # the cached value used by the notebook
    assert rel(bl["Re_c"], cp["blasius"]["Re_c"]) < 1e-12


# ======================================================================================================================
# C05 — free–free boundaries (11.44)
# ======================================================================================================================
def test_free_free_V1_formula_examples_and_chebyshev_route():  # V1 (C05)
    assert [round(float(ch11.benard_free_free_Ra(K)), 2) for K in (1, 2, 3, 4)] == [1284.23, 667.01, 746.53, 1082.06]
    for K in (0.8, 2.2214, 4.0, 7.0):
        assert rel(ch11.benard_marginal_Ra(K, bc=("free", "free")), (PI ** 2 + K ** 2) ** 3 / K ** 2) < 1e-8  # observed ≤ 1.6e-9
    c = ch11.benard_free_free_critical()
    assert c["Ra_c"] == 27 * PI ** 4 / 4 and c["K_c"] == PI / math.sqrt(2) and c["K_c2"] == PI ** 2 / 2
    assert ch11.benard_free_free_Ra(2.0, n=2) > ch11.benard_free_free_Ra(2.0, n=1)  # higher n only raise Ra
    zz = np.linspace(-0.5, 0.5, 11)
    W = ch11.benard_free_free_mode(zz)
    assert abs(W[0]) < 1e-15 and abs(W[-1]) < 1e-15 and np.allclose(W, np.cos(PI * zz), atol=1e-15)
    Wp = ch11.benard_free_free_mode(np.array([-0.5, 0.5]), printed=True)
    assert np.allclose(Wp, [-1.0, 1.0])  # slip #2: sin(πz) ≠ 0 on the walls


def test_free_free_V2_derivation_D11_sine_modes_and_the_minimum():  # V2 (D11 ★★, steps 4–11; slips #2, #3)
    z, K2 = sp.symbols("z K2", positive=True)
    for n_ in (1, 2, 3):
        W = sp.sin(n_ * sp.pi * (z + sp.Rational(1, 2)))
        for m in (0, 2, 4, 6):
            for s in (sp.Rational(1, 2), -sp.Rational(1, 2)):
                assert sp.simplify(sp.diff(W, z, m).subs(z, s)) == 0  # steps 5–7: every even derivative vanishes
        L = lambda f: sp.diff(f, z, 2) - K2 * f  # noqa: E731
        Ra_n = (n_ ** 2 * sp.pi ** 2 + K2) ** 3 / K2
        assert z0(L(L(L(W))) + Ra_n * K2 * W)  # step 8 (11.44)
    assert sp.sin(sp.pi * sp.Rational(1, 2)) != 0  # slip #2 for n = 1
    Ra1 = (sp.pi ** 2 + K2) ** 3 / K2
    dRa = sp.diff(Ra1, K2)
    assert z0(dRa - (3 * (sp.pi ** 2 + K2) ** 2 / K2 - (sp.pi ** 2 + K2) ** 3 / K2 ** 2))  # step 9
    roots = [r for r in sp.solve(dRa, K2) if r.is_positive]
    assert roots == [sp.pi ** 2 / 2] and z0(Ra1.subs(K2, roots[0]) - sp.Rational(27, 4) * sp.pi ** 4)  # steps 10–11
    printed = 3 * (sp.pi ** 2 + K2) ** 2 / K2 - 3 * (sp.pi ** 2 + K2) ** 3 / K2 ** 2
    assert [r for r in sp.solve(printed, K2) if r.is_positive] == []  # slip #3: no root
    e = ch11.benard_free_free_sympy()
    assert z0(e["residual"]) and all(z0(b) for b in e["bc_W4"]) and e["root"] == sp.pi ** 2 / 2
    assert e["printed_root"] == [] and z0(e["Ra_c"] - sp.Rational(27, 4) * sp.pi ** 4)


def test_free_free_V2_derivation_D12_growth_rate_quadratic():  # V2 (D12 ★★, steps 1–8)
    z, s, K, Ra, Pr = sp.symbols("z sigma K Ra Pr", positive=True)
    W0, T0 = sp.symbols("W0 T0")
    S = sp.sin(sp.pi * (z + sp.Rational(1, 2)))
    a2 = sp.pi ** 2 + K ** 2
    e36 = (s + K ** 2) * T0 * S - sp.diff(T0 * S, z, 2) - W0 * S  # (11.36)
    LW = sp.diff(W0 * S, z, 2) - K ** 2 * W0 * S
    e37 = (s / Pr + K ** 2) * LW - sp.diff(LW, z, 2) + Ra * K ** 2 * T0 * S  # (11.37)
    assert z0(e36 / S - ((s + a2) * T0 - W0))  # step 3
    assert z0(e37 / S - (-(s / Pr + a2) * a2 * W0 + Ra * K ** 2 * T0))  # step 4
    quad_ = sp.expand((s + a2) * (s / Pr + a2) * a2 - Ra * K ** 2)  # step 5
    assert z0(sp.expand(quad_ / a2) - (s ** 2 / Pr + a2 * (1 + 1 / Pr) * s + a2 ** 2 - Ra * K ** 2 / a2))  # step 6
    disc = sp.discriminant(sp.expand(quad_ / a2), s)
    assert z0(disc - (a2 ** 2 * (1 - 1 / Pr) ** 2 + 4 * Ra * K ** 2 / (Pr * a2)))  # step 7: > 0 for Ra > 0
    assert z0(quad_.subs(s, 0).subs(Ra, a2 ** 3 / K ** 2))  # step 8: σ = 0 on (11.44)
    rng = np.random.default_rng(6)
    for _ in range(10):
        Kn, Ran, Prn = rng.uniform(0.5, 5), rng.uniform(0, 5000), rng.uniform(0.05, 10)
        roots = np.roots([1 / Prn, (PI ** 2 + Kn ** 2) * (1 + 1 / Prn), (PI ** 2 + Kn ** 2) ** 2 - Ran * Kn ** 2 / (PI ** 2 + Kn ** 2)])
        code = ch11.benard_free_free_sigma(Kn, Ran, Prn)
        assert np.allclose(sorted(np.array(code, dtype=complex), key=lambda v: -v.real),
                           sorted(roots, key=lambda v: -v.real), rtol=1e-10, atol=1e-10)


# ======================================================================================================================
# C06 — double diffusion (11.45)–(11.46)
# ======================================================================================================================
def test_salt_finger_V1_criterion_is_Rs_minus_Ra_and_flips_at_27pi4_over_4():  # V1 (C06, N56)
    sf = ch11.salt_finger_unstable(0.01, 0.002, 0.05)
    assert sf["lhs"] == pytest.approx(sf["Rs"] - sf["Ra"], rel=1e-13)  # (11.46) ⇔ Rs − Ra (§11.5 signs)
    assert sf["margin"] == pytest.approx(ch11.double_diffusive_margin(sf["Ra"], sf["Rs"]), rel=1e-13)
    assert sf["density_stable"] and sf["unstable"] and sf["R_rho"] == pytest.approx(2e-4 * 0.01 / (7.6e-4 * 0.002))
    assert sf["lhs"] == pytest.approx(61233.19, abs=0.01)  # g = G0 (design: 61 254 with g = 9.81)
    assert ch11.salt_finger_unstable(0.01, 0.002, 0.05, g=9.81)["lhs"] == pytest.approx(61254.1, abs=0.1)
    d_thin = brentq(lambda d: ch11.salt_finger_unstable(0.01, 0.002, d)["margin"], 1e-3, 0.05, xtol=1e-12)
    assert d_thin == pytest.approx(0.0161, abs=5e-5)  # design: thinnest finger-unstable layer 1.61 cm
    assert not ch11.salt_finger_unstable(0.01, 0.002, 0.999 * d_thin)["unstable"]
    assert ch11.salt_finger_unstable(0.01, 0.002, 1.001 * d_thin)["unstable"]
    assert ch11.double_diffusive_margin(0.0, ch11.RA_FREE_FREE) == pytest.approx(0, abs=1e-12)
    assert ch11.thermal_rayleigh_signed(0.01, 0.05, 2e-4, 1e-6, 1.4e-7) == pytest.approx(875.59, abs=0.01)
    assert ch11.thermal_rayleigh_signed(-0.01, 0.05, 2e-4, 1e-6, 1.4e-7) < 0  # heated from below: §11.5 Ra < 0
    assert ch11.salinity_rayleigh(0.002, 0.05, 7.6e-4, 1e-6, 1.5e-9) == pytest.approx(
        ch11.salinity_rayleigh_prime(0.002, 0.05, 7.6e-4, 1e-6, 1.4e-7) * 1.4e-7 / 1.5e-9, rel=1e-13)  # Rs = (κ/κ_s)Rs′


def test_salt_finger_V7_single_component_limit():  # V7 (C06: κ_s = κ, no salt → Bénard's 27π⁴/4)
    d, al, nu, ka = 0.01, 2e-4, 1e-6, 1.4e-7
    dTdz = -ch11.RA_FREE_FREE * nu * ka / (G0 * al * d ** 4)  # exactly marginal, heated from below
    sf = ch11.salt_finger_unstable(dTdz, 0.0, d, alpha=al, nu=nu, kappa=ka, kappa_s=ka)
    assert sf["margin"] == pytest.approx(0, abs=1e-9)
    assert ch11.rayleigh_number(al, -dTdz * d, d, ka, nu) == pytest.approx(ch11.RA_FREE_FREE, rel=1e-13)  # (11.21) sign
    assert not sf["density_stable"]


def test_double_diffusive_sigma_V1_reduces_to_benard_and_roots_solve_the_cubic():  # V1 + V7 (N54, D13)
    for K2 in (PI ** 2 / 2, 1.0, 9.0):
        a2 = PI ** 2 + K2
        r = ch11.double_diffusive_sigma(K2, -2000.0, 0.0, 1.0, 1.0)
        bf = ch11.benard_free_free_sigma(math.sqrt(K2), 2000.0, 1.0)
        assert np.min(np.abs(r - bf[0])) < 1e-9 * max(1, abs(bf[0])) and np.min(np.abs(r - bf[1])) < 1e-8 * abs(bf[1])
        assert np.min(np.abs(r + a2)) < 1e-9  # the extra root of τ = 1 is σ = −a²
    for Ra, Rs, tau in ((1000.0, 2000.0, TAU_SW), (-2e4, -1.9e6, TAU_SW), (500.0, -3000.0, 0.3)):
        a2 = PI ** 2 / 2 + PI ** 2
        for s in ch11.double_diffusive_sigma(PI ** 2 / 2, Ra, Rs, 7.0, tau):
            lhs = (s / 7.0 + a2) * a2 * (s + a2) * (s + tau * a2)
            rhs = PI ** 2 / 2 * (-Ra * (s + tau * a2) + tau * Rs * (s + a2))
            assert abs(lhs - rhs) < 1e-9 * max(abs(lhs), 1)
    s = ch11.double_diffusive_sigma(PI ** 2 / 2, 1000.0, 2000.0, 7.0, TAU_SW)
    assert abs(s[0] - 0.0330081) < 1e-6 and abs(s[0].imag) == 0  # fingers: real positive root (NOT the design's 0.0308)
    s = ch11.double_diffusive_sigma(PI ** 2 / 2, -2e4, -1.9e6, 7.0, TAU_SW)
    assert abs(s[0].real - 9.650) < 1e-3 and abs(abs(s[0].imag) - 70.389) < 1e-3  # diffusive: oscillatory growth


def test_double_diffusion_V2_derivation_D13_cubic_from_the_dimensional_equations():  # V2 (D13 ★★, steps 1–12)
    sd, kd, d, g, al, be, nu, ka, ks, Tz, Sz = sp.symbols("sigma_d k_d d g alpha beta nu kappa kappa_s T_z S_z")
    m = sp.pi / d  # free–free sine mode sin(π(z/d + ½)): ∂²/∂z² → −m², ∇_H² → −k_d²
    lap = -(kd ** 2 + m ** 2)
    wh, Th, sh = sp.symbols("w_h T_h s_h")
    eqT = sd * Th + wh * Tz - ka * lap * Th  # step 1 (heat), dimensional
    eqS = sd * sh + wh * Sz - ks * lap * sh  # step 1 (salt)
    eqW = sd * lap * wh - g * (-kd ** 2) * (al * Th - be * sh) - nu * lap ** 2 * wh  # step 3: ∂_t∇²w = g∇_H²(αT − βs) + ν∇⁴w
    M = sp.Matrix([[sp.diff(e, v) for v in (wh, Th, sh)] for e in (eqT, eqS, eqW)])
    det_dim = sp.expand(M.det())
    s, K2, Pr, tau, Ra, Rs = sp.symbols("sigma K2 Pr tau Ra Rs")
    sub = {sd: s * ka / d ** 2, kd: sp.sqrt(K2) / d, nu: Pr * ka, ks: tau * ka,
           Tz: Ra * nu * ka / (g * al * d ** 4), Sz: Rs * nu * ks / (g * be * d ** 4)}  # step 4–5 scales, §11.5 signs
    sub[Tz] = sub[Tz].subs(nu, Pr * ka)
    sub[Sz] = sub[Sz].subs({nu: Pr * ka, ks: tau * ka})
    nd = sp.simplify(det_dim.subs(sub))
    a2 = sp.pi ** 2 + K2
    cubic = (s / Pr + a2) * a2 * (s + a2) * (s + tau * a2) - K2 * (-Ra * (s + tau * a2) + tau * Rs * (s + a2))  # ours
    ratio = sp.simplify(nd / cubic)
    assert ratio.free_symbols <= {ka, d, Pr, g, al, be}, ratio  # same polynomial in σ up to a σ-free factor
    assert z0(cubic.subs(s, 0) - tau * a2 * (a2 ** 3 - K2 * (Rs - Ra)))  # steps 9–12: σ = 0 ⇒ Rs − Ra = a⁶/K²
    assert z0(((sp.pi ** 2 + K2) ** 3 / K2).subs(K2, sp.pi ** 2 / 2) - sp.Rational(27, 4) * sp.pi ** 4)  # (11.46)
    vals = {K2: sp.pi ** 2 / 2, Pr: 7, tau: sp.Rational(15, 1400), Ra: 1000, Rs: 2000}
    roots = sp.Poly(sp.expand(cubic.subs(vals)), s).nroots(n=30)
    real_pos = [complex(r) for r in roots if abs(sp.im(r)) < 1e-20 and sp.re(r) > 0]
    assert len(real_pos) == 1 and abs(real_pos[0] - 0.0330081) < 1e-6  # decides 0.03301 (code) vs 0.0308 (design)
    alt = sp.Poly(sp.expand(cubic.subs({**vals, tau: sp.Rational(1, 100)})), s).nroots(n=30)
    assert any(abs(complex(r) - 0.0307930) < 1e-6 for r in alt)  # the design's 0.0308 is τ = 0.01, not 0.0107


def test_salt_finger_regime_V7_four_regimes_and_eos():  # V7 (N53, N54)
    assert ch11.salt_finger_regime(0.01, 0.002, 0.05)["regime"] == "fingers"
    assert ch11.salt_finger_regime(-0.01, 0.0, 0.05)["regime"] == "overturning"
    assert ch11.salt_finger_regime(0.01, -0.002, 0.05)["regime"] == "stable"
    # cold fresh over hot salty, statically stable: R_ρ just above 1 → oscillatory growth (diffusive regime)
    d, al, be = 0.05, 2e-4, 7.6e-4
    Tz = -2e4 * 1e-6 * 1.4e-7 / (G0 * al * d ** 4)
    Sz = -1.9e6 * 1e-6 * 1.5e-9 / (G0 * be * d ** 4)
    r = ch11.salt_finger_regime(Tz, Sz, d)
    assert r["density_stable"] and r["regime"] == "diffusive" and r["sigma_max"].real > 0 and abs(r["sigma_max"].imag) > 1
    rho = ch11.linear_eos(np.array([10.0, 12.0]), np.array([35.0, 35.5]))
    assert rho == pytest.approx(1027.0 * np.array([1.0, 1 - 2e-4 * 2 + 7.6e-4 * 0.5]), rel=1e-14)


# ======================================================================================================================
# C07 — Taylor–Couette (11.47)–(11.54)
# ======================================================================================================================
def test_taylor_critical_V1_mu_to_1_is_rigid_benard():  # V1 (C07, D15) + V5-lite (arXiv:2601.14806)
    t1 = ch11.taylor_critical(1.0)
    b = ch11.benard_critical()
    assert rel(t1["Ta_c"], b["Ra_c"]) < 1e-8 and abs(t1["k_c"] - b["K_c"]) < 1e-5
    for k in (2.5, 3.1, 4.0):
        assert rel(ch11.taylor_marginal_Ta(k, 1.0), ch11.benard_marginal_Ra_det(k)) < 1e-7  # same eigenproblem (observed ≤ 2.5e-8)


def test_taylor_marginal_V1_galerkin_route_agrees_and_slip7_fails():  # V1 independent route (N124, Ex. 11.9)
    for k, mu in ((3.0, 0.0), (3.1266, 0.0), (3.2, 0.5), (3.1, 1.0), (3.3, -0.25)):
        a, b = ch11.taylor_galerkin_Ta(k, mu, n_modes=4), ch11.taylor_marginal_Ta(k, mu)
        assert rel(a, b) < 1e-3, (k, mu, a, b)
    g6 = ch11.taylor_galerkin_Ta(3.1266, 0.0, n_modes=8)
    assert rel(g6, ch11.taylor_marginal_Ta(3.1266, 0.0)) < 2e-5  # more modes → closer
    assert math.isnan(ch11.taylor_galerkin_Ta(3.1266, 0.0, printed=True))  # slip #7: no positive Ta at all


def test_taylor_marginal_V3_converges_in_N():  # V3 (C07)
    ref = ch11.taylor_marginal_Ta(3.1266, 0.0, N=24)
    errs = [abs(ch11.taylor_marginal_Ta(3.1266, 0.0, N=N) - ref) / ref for N in (10, 14, 20)]
    assert errs[0] > 1e-6 and errs[1] < 1e-8 and errs[2] < 1e-9, errs  # 9.5e-6, 1.5e-9, 3e-12
    assert abs(ch11.taylor_marginal_Ta(3.1266, 0.0) - ref) / ref < 1e-9  # default N = 40 (round-off grows ~N⁴: 6.5e-9 at N = 60)


def test_taylor_critical_V1_11_54_accuracy_and_cached_table_is_live():  # V1 + anti-cache (C07, N68)
    for mu, Ta, kc in ((0.0, 3389.90, 3.1266), (0.5, 2275.09, 3.1175), (-0.5, 6413.72, 3.1985)):
        r = ch11.taylor_critical(mu)
        assert abs(r["Ta_c"] - Ta) < 0.01 and abs(r["k_c"] - kc) < 1e-4
    assert ch11.taylor_critical_approx(0.0) / ch11.taylor_critical(0.0)["Ta_c"] - 1 == pytest.approx(0.0077, abs=1e-4)
    assert ch11.taylor_critical_approx(0.5) / ch11.taylor_critical(0.5)["Ta_c"] - 1 == pytest.approx(0.0010, abs=1e-4)
    assert ch11.taylor_critical_approx(1.0) == 1708.0
    t = read_csv("taylor_critical.csv")
    for i in (0, 10, 30):
        live = ch11.taylor_critical(t["mu"][i])
        assert rel(t["Ta_c"][i], live["Ta_c"]) < 1e-5 and abs(t["k_c"][i] - live["k_c"]) < 1e-4
        assert rel(t["Ta_11_54"][i], ch11.taylor_critical_approx(t["mu"][i])) < 1e-5
    tab = ch11.taylor_critical_table([0.0, 1.0], cache=False)
    assert rel(tab["Ta_c"][1], 1707.762) < 1e-5 and tab["rel_error"][0] == pytest.approx(0.0077, abs=1e-4)
    with pytest.raises(ValueError):
        ch11.taylor_critical_approx(-1.0)


def test_taylor_growth_rate_V7_zero_at_the_margin_and_real_for_corotation():  # V7 (N67)
    for mu in (0.0, 0.5):
        Tam = ch11.taylor_marginal_Ta(3.12, mu)
        assert abs(ch11.taylor_growth_rate(3.12, Tam, mu)) < 1e-6
        assert ch11.taylor_growth_rate(3.12, 1.05 * Tam, mu).real > 0 > ch11.taylor_growth_rate(3.12, 0.95 * Tam, mu).real
        s = ch11.taylor_growth_rate(3.12, 1.2 * Tam, mu, all=True)
        assert np.max(np.abs(s[:10].imag) / np.maximum(1.0, np.abs(s[:10]))) < 1e-8  # observed ≤ 1.9e-9 (round-off)
    assert ch11.taylor_growth_rate(3.12, 3500.0, 0.0).real == pytest.approx(0.4214, abs=1e-4)


def test_taylor_eigenfunction_V1_walls_and_narrow_gap_continuity():  # V1 (N66, Fig. 11.16)
    e = ch11.taylor_eigenfunction(3.1266, 0.0, x=np.linspace(0, 1, 101), z=np.linspace(0, 2, 201))
    assert max(abs(e["uR_profile"][0]), abs(e["uR_profile"][-1]), abs(e["uphi_profile"][0]), abs(e["uphi_profile"][-1])) < 1e-10
    D, _ = ST.cheb(40, (0.0, 1.0))
    assert max(abs((D @ e["uR_profile"])[0]), abs((D @ e["uR_profile"])[-1])) < 1e-9  # (11.53)
    hx, hz = 0.01, 0.01
    dux = np.gradient(e["u_R"], hx, axis=1)
    div = dux + np.gradient(e["u_z"], hz, axis=0)
    sc = np.max(np.abs(dux))
    assert np.max(np.abs(div[2:-2, 2:-2])) < 2e-3 * sc  # narrow-gap continuity (second-order FD error only)
    assert np.max(np.abs(np.gradient(e["psi"], hz, axis=0) - e["u_R"])[2:-2, 2:-2]) < 2e-3  # u_R = ∂ψ/∂z
    assert np.max(np.abs(-np.gradient(e["psi"], hx, axis=1) - e["u_z"])[2:-2, 2:-2]) < 2e-3 * sc  # u_z = −∂ψ/∂x
    assert e["Ta"] == pytest.approx(3389.90, abs=0.01)


def test_taylor_number_V1_examples_couette_parity_and_V7_dimensionless():  # V1 + V7 (N65, D15, R16)
    t = ch11.taylor_number(0.5, 0.0, 0.1, 0.105, 1e-6)
    assert t["Ta"] == pytest.approx(6097.56, abs=0.01) and not t["rayleigh_stable"] and t["mu"] == 0.0
    assert ch11.taylor_number_narrow_inner(0.5, 0.1, 0.005, 1e-6) == pytest.approx(6250.0)
    for dR in (1e-3, 1e-4):  # narrow-gap limit: relative difference O(d/R₁)
        exact = ch11.taylor_number(0.5, 0.0, 1.0, 1.0 + dR, 1e-8)["Ta"]
        assert abs(ch11.taylor_number_narrow_inner(0.5, 1.0, dR, 1e-8) / exact - 1) < 1.5 * dR
    _, A, _ = LAM.circular_couette(0.1, 0.1, 0.105, 0.5, 0.2, return_coeffs=True)
    assert ch11.taylor_number(0.5, 0.2, 0.1, 0.105, 1e-6)["Ta"] == pytest.approx(-4 * A * 0.5 * 0.005 ** 4 / 1e-12, rel=1e-12)
    lam, tau = 100.0, 7.0  # lengths ×100, time ×7: a dimensionless number does not change
    a = ch11.taylor_number(0.5, 0.1, 0.1, 0.105, 1e-6)["Ta"]
    b = ch11.taylor_number(0.5 / tau, 0.1 / tau, 0.1 * lam, 0.105 * lam, 1e-6 * lam ** 2 / tau)["Ta"]
    assert b == pytest.approx(a, rel=1e-12)
    mR = ch11.couette_rayleigh_line(0.1, 0.105)
    assert mR == pytest.approx(0.9070, abs=1e-4)
    assert ch11.taylor_number(1.0, 1.001 * mR, 0.1, 0.105, 1e-6)["Ta"] < 0  # beyond Rayleigh's line: no forcing (sign)
    assert ch11.taylor_number(1.0, 0.999 * mR, 0.1, 0.105, 1e-6)["Ta"] > 0
    with pytest.raises(ValueError):
        ch11.taylor_number(1.0, 0.0, 0.1, 0.1, 1e-6)


def test_ring_interchange_and_rayleigh_criterion_V1_V7():  # V1 + V7 (N59, R14)
    r = ch11.ring_interchange_energy(4.0, 2.0, 1.0, 2.0)
    assert r["dE"] == pytest.approx(r["E_f"] - r["E_i"], rel=1e-14)
    assert (r["E_i"], r["E_f"], r["dE"]) == pytest.approx((0.21531, 0.10132, -0.11399), abs=1e-5)
    assert ch11.ring_interchange_energy(2.0, 4.0, 1.0, 2.0)["dE"] > 0  # Γ² increasing outward: costs energy (stable)
    R = np.linspace(1.0, 1.5, 200)
    assert not ch11.rayleigh_circulation_criterion(R, lambda rr: 0.7 * rr)["unstable"]  # solid body
    out_rest = LAM.circular_couette(R, 1.0, 1.5, 1.0, 0.0)
    crit = ch11.rayleigh_circulation_criterion(R, out_rest)
    assert crit["unstable"] and crit["where"].size == R.size - 1
    mR = ch11.couette_rayleigh_line(1.0, 1.5)
    for f, unstable in ((1.02, False), (0.98, True)):
        u = LAM.circular_couette(R, 1.0, 1.5, 1.0, f * mR)
        assert ch11.rayleigh_circulation_criterion(R, u)["unstable"] is unstable
    with pytest.raises(ValueError):
        ch11.ring_interchange_energy(1, 1, 2.0, 1.0)


def test_taylor_stability_boundary_V1_points_reproduce_ta_c():  # V1 (N68, Fig. 11.17 analogue)
    sb = ch11.taylor_stability_boundary(1.05, mus=[-0.5, 0.0, 0.5])
    R1, R2 = 1.0, 1.05
    for xo, yi, mu, Tac in zip(sb["x_outer"], sb["y_inner"], sb["mu"], sb["Ta_c"]):
        Om1, Om2 = yi / R2 ** 2, xo / R2 ** 2  # ν = 1 units
        assert Om2 / Om1 == pytest.approx(mu)
        assert ch11.taylor_number(Om1, Om2, R1, R2, 1.0)["Ta"] == pytest.approx(Tac, rel=1e-10)
    assert np.allclose(sb["rayleigh_y"][1:] / sb["rayleigh_x"][1:], (R2 / R1) ** 2)


def test_taylor_V2_derivation_D14_linearisation_elimination_and_narrow_gap():  # V2 (D14 ★★★, steps 1–15; slip #12)
    R, z, t = sp.symbols("R z t", positive=True)
    A, B, rho, nu, eps = sp.symbols("A B rho nu epsilon")
    U = A * R + B / R  # (11.49)
    uR, uph, uz, p = [sp.Function(n)(R, z, t) for n in ("u_R", "u_phi", "u_z", "p")]
    P = sp.Function("P")(R)
    lap = lambda f: sp.diff(f, R, 2) + sp.diff(f, R) / R + sp.diff(f, z, 2)  # noqa: E731
    tR, tph, tz = eps * uR, U + eps * uph, eps * uz
    adv = lambda f: sp.diff(f, t) + tR * sp.diff(f, R) + tz * sp.diff(f, z)  # noqa: E731
    eR = adv(tR) - tph ** 2 / R + sp.diff(P + eps * p, R) / rho - nu * (lap(tR) - tR / R ** 2)  # (11.47) radial
    eph = adv(tph) + tR * tph / R - nu * (lap(tph) - tph / R ** 2)
    eR1 = sp.diff(eR, eps).subs(eps, 0)
    assert z0(eR1 - (sp.diff(uR, t) - 2 * U * uph / R + sp.diff(p, R) / rho - nu * (lap(uR) - uR / R ** 2)))  # steps 2–3
    assert z0(sp.diff(eph, eps).subs(eps, 0) - (sp.diff(uph, t) + (sp.diff(U, R) + U / R) * uR - nu * (lap(uph) - uph / R ** 2)))
    assert z0(sp.diff(U, R) + U / R - 2 * A)  # step 9
    assert z0(sp.expand((U + eps) ** 2 / R) - sp.expand(U ** 2 / R) - 2 * U * eps / R - eps ** 2 / R)  # step 2: the factor 2
    # steps 5–10 with normal modes e^{ikz + σt}: û_z and p̂ eliminated
    k, s = sp.symbols("k sigma", positive=True)
    f, h = sp.Function("f")(R), sp.Function("h")(R)  # f = û_R, h = û_φ
    D = lambda q: sp.diff(q, R)  # noqa: E731
    Ds = lambda q: sp.diff(q, R) + q / R  # noqa: E731
    assert z0(D(Ds(f)) - (f.diff(R, 2) + f.diff(R) / R - f / R ** 2)) and z0(Ds(D(f)) - (f.diff(R, 2) + f.diff(R) / R))  # step 6
    assert z0(D(Ds(D(f)) - k ** 2 * f) - (D(Ds(D(f))) - k ** 2 * D(f)))  # step 10: D(D_*D − k²) = (DD_* − k²)D
    Lc = lambda q: D(Ds(q)) - k ** 2 * q  # noqa: E731
    uzh = sp.I / k * Ds(f)  # step 7
    assert z0(Ds(f) + sp.I * k * uzh)
    ph = (nu * (Ds(D(Ds(f))) - k ** 2 * Ds(f)) - s * Ds(f)) / k ** 2  # step 8 (p̂/ρ)
    axial = s * uzh + sp.I * k * ph - nu * (Ds(D(uzh)) - k ** 2 * uzh)
    assert z0(axial)
    radial = s * f - 2 * U / R * h + D(ph) - nu * Lc(f)
    step10 = nu * Lc(Lc(f)) - s * Lc(f) - 2 * k ** 2 * U / R * h
    assert z0(k ** 2 * radial - step10)  # step 10
    # steps 11–15: narrow gap, scale by d, rescale û_φ, read off Ta
    x, d, Om1, al, kt, st = sp.symbols("x d Omega1 alpha k_t sigma_t", positive=True)
    F, Hh = sp.Function("F")(x), sp.Function("H")(x)
    Lx = lambda q: q.diff(x, 2) - kt ** 2 * q  # noqa: E731
    first = (nu / d ** 2) * (Lx(Lx(F)) - st * Lx(F)) / d ** 2 - 2 * (kt / d) ** 2 * Om1 * (1 + al * x) * Hh
    cst = nu / (2 * kt ** 2 * d ** 2 * Om1)
    resc1 = sp.expand(first.subs(Hh, cst * Hh).doit() * d ** 4 / nu)
    assert z0(resc1 - (Lx(Lx(F)) - st * Lx(F) - (1 + al * x) * Hh))  # (11.51) first line
    second = (nu / d ** 2) * (Lx(Hh) - st * Hh) - 2 * A * F
    resc2 = sp.expand(second.subs(Hh, cst * Hh).doit() * d ** 2 / nu / cst)
    Ta = -4 * A * Om1 * d ** 4 / nu ** 2
    assert z0(resc2 - (Lx(Hh) - st * Hh + Ta * kt ** 2 * F))  # (11.51) second line, Ta = −4AΩ₁d⁴/ν²
    O1, O2, R1, R2, nv = sp.symbols("Omega1 Omega2 R1 R2 nu_v", positive=True)
    A_ = (O2 * R2 ** 2 - O1 * R1 ** 2) / (R2 ** 2 - R1 ** 2)
    assert z0(-4 * A_ * O1 * (R2 - R1) ** 4 / nv ** 2 - 4 * (O1 * R1 ** 2 - O2 * R2 ** 2) / (R2 ** 2 - R1 ** 2) * O1 * (R2 - R1) ** 4 / nv ** 2)
    # slip #12: with u_R = V g(R/L), u_z = V q(z/L) the printed ∂(Rũ_R)/∂R carries L⁰ and ∂ũ_z/∂z carries L⁻¹
    Lsc, V, sv = sp.symbols("L V s", positive=True)
    gfun, qfun = sp.Function("g"), sp.Function("q")
    t1_printed = sp.diff(R * V * gfun(R / Lsc), R).subs(R, Lsc * sv)
    t1_correct = (sp.diff(R * V * gfun(R / Lsc), R) / R).subs(R, Lsc * sv)
    t2 = sp.diff(V * qfun(z / Lsc), z).subs(z, Lsc * sv)
    deg = lambda e: sp.simplify(Lsc * sp.diff(e, Lsc) / e)  # noqa: E731  (the power of L in a monomial-in-L term)
    assert deg(t1_printed) == 0 and deg(t2) == -1 and deg(t1_correct) == -1
    tp = ch11.taylor_perturbation_sympy()
    for key in ("lin_R", "lin_phi", "lin_z", "base_R", "base_phi", "term_2A", "operator_identity", "Ta_identity",
                "continuity_correct"):
        assert z0(tp[key]), key
    assert tp["printed_continuity_units_ok"] is False and tp["correct_continuity_units_ok"] is True
    assert not z0(tp["continuity_printed"])


# ======================================================================================================================
# C08 — Taylor–Goldstein (11.61); C09 — Miles–Howard (11.67); C10 — Howard's semicircle (11.72)
# ======================================================================================================================
def _tg_mode(k=0.4, J=0.1, N=100):
    pr = ch11.richardson_profiles("tanh", J)
    res = ST.taylor_goldstein_eigs(k, pr["U"], pr["Upp"], pr["N2"], domain=(-1, 1), N=N, bc="decay", return_vectors=True)
    return pr, res


def test_taylor_goldstein_V1_N2_zero_equals_rayleigh_and_V5_michalke():  # V1 (C08, N91) + V5 (Michalke 1964, approx.)
    pr = ch11.richardson_profiles("tanh", 0.0)
    for k in (0.2, 0.4449, 0.7):
        ym = ST.decay_box(k)  # loop 1: both solvers in the converged box (60 at k = 0.2; the fixed 30 is 2.3e-6 off in c there)
        ctg = ST.taylor_goldstein_eigs(k, pr["U"], pr["Upp"], pr["N2"], domain=(-1, 1), N=100, bc="decay", y_max=ym)
        cr = ST.rayleigh_eigs(k, TANH["U"], TANH["Upp"], N=100, bc="decay", y_max=ym, map_scale=0.5, unstable_only=True)
        assert abs(ctg[0] - cr[0]) < 1e-8 and abs(ctg[0].real) < 1e-10
    m = ch11.tanh_max_growth()
    assert abs(m["k"] - 0.4446) / 0.4446 < 1e-3 and abs(m["kci"] - 0.1897) / 0.1897 < 1e-3
    assert ch11.tg_growth(0.4449, 0.0) == pytest.approx(0.18970, abs=1e-5)
    assert ch11.tg_growth(0.4, 0.1) == pytest.approx(0.1245, abs=1e-4)


def test_taylor_goldstein_V1_exact_neutral_curve_J_equals_k_one_minus_k():  # V1 (N78; Drazin/Hazel neutral mode)
    zz = np.linspace(0.2, 6, 400)  # z > 0 (the mode is |tanh z|^{1−k} sech^k z)
    for k in (0.3, 0.5, 0.7):
        J = float(ch11.tg_tanh_neutral_J(k))
        assert J == pytest.approx(k * (1 - k))
        h = 1e-4
        psi = lambda q: ch11.tg_tanh_neutral_mode(q, k)  # noqa: E731
        d2 = (psi(zz + h) - 2 * psi(zz) + psi(zz - h)) / h ** 2
        U, Upp, N2 = np.tanh(zz), -2 * np.tanh(zz) / np.cosh(zz) ** 2, J / np.cosh(zz) ** 2
        res = U * (d2 - k ** 2 * psi(zz)) - Upp * psi(zz) + N2 * psi(zz) / U  # (11.61) with c = 0
        assert np.max(np.abs(res)) < 1e-6
        assert ch11.tg_growth(k, J - 0.05) > 0 and ch11.tg_growth(k, J + 0.01) == 0.0  # the tongue's edge
    zs, ks = sp.symbols("z k", positive=True)
    psi_s = sp.tanh(zs) ** (1 - ks) / sp.cosh(zs) ** ks
    U = sp.tanh(zs)
    tg = U * (sp.diff(psi_s, zs, 2) - ks ** 2 * psi_s) - sp.diff(U, zs, 2) * psi_s + ks * (1 - ks) / sp.cosh(zs) ** 2 * psi_s / U
    for kv in (sp.Rational(3, 10), sp.Rational(3, 5)):  # exact identity checked with 40-digit arithmetic (V2)
        for zv in (sp.Rational(1, 7), sp.Rational(9, 10), 3):
            assert abs(sp.N((tg / psi_s).subs({ks: kv, zs: zv}), 40)) < 1e-30


def test_taylor_goldstein_V7_conjugate_pairs_and_miles_howard_no_growth():  # V7 (N73, C09)
    pr = ch11.richardson_profiles("tanh", 0.1)
    c = ST.taylor_goldstein_eigs(0.4, pr["U"], pr["Upp"], pr["N2"], domain=(-1, 1), N=80, bc="decay", unstable_only=False,
                                 filter=False)
    big = c[np.abs(c.imag) > 1e-3]
    for v in big:
        assert np.min(np.abs(c - v.conjugate())) < 1e-6 * max(1, abs(v))  # real coefficients ⇒ c* is an eigenvalue
    for J in (0.26, 0.3, 0.5):
        assert all(ch11.tg_growth(k, J, N=80) == 0.0 for k in (0.1, 0.3, 0.5, 0.7, 0.9))
    rng = np.random.default_rng(7)
    for _ in range(6):  # random smooth profiles with Ri_min > 1/4 on a walled channel
        a, b_, J = rng.uniform(0.5, 2.0), rng.uniform(-0.5, 0.5), rng.uniform(0.27, 1.0)
        U = lambda zz, a=a, b_=b_: np.tanh(a * zz) + b_ * zz  # noqa: E731
        Up = lambda zz, a=a, b_=b_: a / np.cosh(a * zz) ** 2 + b_  # noqa: E731
        Upp = lambda zz, a=a: -2 * a ** 2 * np.tanh(a * zz) / np.cosh(a * zz) ** 2  # noqa: E731
        N2 = lambda zz, J=J: J * Up(zz) ** 2 + 0.0 * zz  # noqa: E731
        assert ch11.miles_howard_stable(np.linspace(-1, 1, 201), dUdz=Up, N2=N2)["guaranteed_stable"]
        for k in (0.5, 1.5, 3.0):
            assert len(ST.taylor_goldstein_eigs(k, U, Upp, N2, domain=(-1, 1), N=60)) == 0


def test_taylor_goldstein_V4_integral_identities_of_computed_modes():  # V4 (C09 (11.65), C10 (11.69)–(11.70), N125)
    for k, J in ((0.4, 0.1), (0.3, 0.05), (0.5, 0.15)):
        pr, res = _tg_mode(k, J)
        c, psi = res["c"][0], res["psi"][:, 0]
        ri = ch11.richardson_identity_check(k, c, psi, res["y"], pr["U"], pr["Up"], pr["Upp"], pr["N2"], grid=res["grid"])
        assert ri["residual"] < 1e-8 and ri["imag_lhs"] == pytest.approx(ri["imag_rhs"], rel=1e-7)
        assert ri["imag_lhs"] < 0  # = −c_i × a positive integral
        hw = ch11.howard_identity_check(k, c, psi, res["y"], pr["U"], pr["N2"], grid=res["grid"], Up=pr["Up"])
        assert hw["res_69"] < 1e-8 and hw["res_70"] < 1e-8 and hw["cr_mean"] == pytest.approx(c.real, abs=1e-8)
        assert hw["int_11_72"] >= 0 and hw["in_semicircle"]
        eb = ch11.stratified_energy_budget(k, c, psi, res["y"], pr["U"], pr["N2"], grid=res["grid"], Up=pr["Up"])
        assert abs(eb["residual"]) < 1e-9 * max(abs(eb["production"]), 1e-12)  # Ex. 11.10: dK/dt + dP/dt = production
        assert eb["APE"] > 0 and eb["production"] > 0


def test_richardson_V1_closed_form_and_miles_howard_verdict():  # V1 (R20, C09, N78)
    z = np.linspace(-3, 3, 121)
    pr = ch11.richardson_profiles("tanh", 0.2)
    Ri = ch11.gradient_richardson(z, U=pr["U"], N2=pr["N2"])
    assert np.allclose(Ri, 0.2 * np.cosh(z) ** 2, rtol=1e-6)
    Ri2 = ch11.gradient_richardson(z, N2=pr["N2"], dUdz=pr["Up"])
    assert np.allclose(Ri2, 0.2 * np.cosh(z) ** 2, rtol=1e-12)
    Ri3 = ch11.gradient_richardson(z, U=np.tanh(z), N2=pr["N2"](z))
    assert np.allclose(Ri3[10:-10], 0.2 * np.cosh(z[10:-10]) ** 2, rtol=1e-2)
    assert np.allclose(pr["Ri"](z), 0.2 * np.cosh(z) ** 2)
    v = ch11.miles_howard_stable(z, U=pr["U"], N2=pr["N2"])
    assert not v["guaranteed_stable"] and v["Ri_min"] == pytest.approx(0.2, rel=1e-6) and abs(v["z_min"]) < 1e-12
    assert ch11.miles_howard_stable(z, dUdz=pr["Up"], N2=ch11.richardson_profiles("tanh", 0.26)["N2"])["guaranteed_stable"]
    assert ch11.gradient_richardson(np.array([0.0]), dUdz=lambda q: 0 * q, N2=lambda q: 1 + 0 * q)[0] == np.inf
    with pytest.raises(ValueError):
        ch11.gradient_richardson(z, U=pr["U"])


def test_tg_growth_map_V1_cached_csv_matches_live_and_the_exact_tongue():  # V1 anti-cache (E6/F5 table)
    gm = ch11.tg_growth_map()
    neutral = gm["k"] * (1 - gm["k"])
    KK, JJ = np.meshgrid(gm["k"], gm["J"])
    assert np.all(gm["kci"][JJ >= 0.25] == 0)  # Miles–Howard
    assert np.all(gm["kci"] >= 0)
    ii, jj = np.nonzero(gm["kci"] > 0)
    outside = [(float(gm["k"][j]), float(gm["J"][i]), float(gm["kci"][i, j]), float(neutral[j]))
               for i, j in zip(ii, jj) if gm["J"][i] >= neutral[j]]
    assert len(ii) > 50 and not outside, f"growth reported outside the exact tongue J < k(1 − k): (k, J, kc_i, k(1−k)) = {outside}"
    assert np.all(KK.shape == gm["kci"].shape)
    for (i, j) in ((0, 8), (10, 7), (5, 3), (20, 9), (2, 15)):
        assert gm["kci"][i, j] == pytest.approx(ch11.tg_growth(gm["k"][j], gm["J"][i]), abs=2e-4)  # 4 s.f. table
    small = ch11.tg_growth_map(ks=[0.4], Js=[0.1], cache=False)
    assert small["kci"][0, 0] == pytest.approx(0.1245, abs=1e-4)


def test_tg_growth_V3_converged_in_the_decay_box_at_small_k():  # V3 (C09 table E6/F5: box truncation y_max)
    """kc_i of the tanh / J sech² layer must not depend on where the decaying condition is imposed.  The long modes decay
    like e^{−k|z|}, so a box y_max = 30 is only 1.5 e-folds wide at k = 0.05."""
    for k, J in ((0.05, 0.0), (0.05, 0.04), (0.1, 0.05), (0.3, 0.1)):
        ref = ch11.tg_growth(k, J, y_max=240.0)
        assert ch11.tg_growth(k, J, y_max=120.0) == pytest.approx(ref, rel=1e-3)  # the reference is converged
        assert ch11.tg_growth(k, J) == pytest.approx(ref, rel=1e-2), (k, J, ch11.tg_growth(k, J), ref)  # default box
    assert ch11.tg_growth(0.05, 0.05) == 0.0  # J = 0.05 > k(1 − k) = 0.0475: outside the exact tongue


def test_tg_growth_V3_no_false_zero_under_the_neutral_curve_at_large_k():  # V3 (C08 table E6/F5: resolution in N)
    """The table and ``tg_growth`` document "approximate within ≈ 0.02 of the exact neutral curve J = k(1 − k)" (modes dropped
    by the N-filter are reported as 0).  So the zero strip under the neutral curve may be at most 0.02 + one J step (0.01)
    wide in every column, and a point further inside must carry the converged growth (N = 160 and 240 agree to 1 %)."""
    gm = ch11.tg_growth_map()
    strip = {}
    for j, k in enumerate(gm["k"]):
        pos = np.nonzero(gm["kci"][:, j] > 0)[0]
        if k * (1 - k) > 0.03 + 1e-9:
            strip[round(float(k), 2)] = float(k * (1 - k) - (gm["J"][pos[-1]] if len(pos) else -np.inf))
    too_wide = {k: round(w, 4) for k, w in strip.items() if w > 0.03 + 1e-9}
    assert not too_wide, f"zero strip under J = k(1 − k) wider than the documented 0.02 (+ 0.01 grid step): {too_wide}"
    for k, J in ((0.9, 0.04), (0.8, 0.13), (0.95, 0.01)):  # 0.05, 0.03, 0.0375 inside the tongue
        ref = ch11.tg_growth(k, J, N=240)
        assert ref > 0.02 and ch11.tg_growth(k, J, N=160) == pytest.approx(ref, rel=1e-2)  # the reference is converged
        assert ch11.tg_growth(k, J) == pytest.approx(ref, rel=1e-2), (k, J, ch11.tg_growth(k, J), ref)


def test_decay_box_V1_formula_guards_and_efolds():  # V1 (loop 1: core.stability.decay_box = max(30, 12/k))
    assert ST.decay_box(0.05) == pytest.approx(240.0, rel=1e-14) and ST.decay_box(0.1) == pytest.approx(120.0, rel=1e-14)
    assert ST.decay_box(0.4) == 30.0 and ST.decay_box(1.0) == 30.0 and ST.decay_box(0.2) == pytest.approx(60.0, rel=1e-14)
    assert ST.decay_box(0.5, y_min=10.0, n_efold=8.0) == 16.0 and ST.decay_box(0.1, 5.0, 3.0) == pytest.approx(30.0)
    assert ch11.decay_box is ST.decay_box and np.ndim(ST.decay_box(0.3)) == 0  # one definition, scalar-callable
    ks = np.geomspace(1e-3, 5.0, 60)
    box = np.array([ST.decay_box(k) for k in ks])
    assert np.all(box >= 30.0) and np.all(np.exp(-ks * box) <= math.exp(-12.0) * (1 + 1e-12))  # ≥ 12 e-folds of e^{−k|y|}
    assert np.all(np.diff(box) <= 0)  # longer waves never get a smaller box
    for bad in (0.0, -0.3):
        with pytest.raises(ValueError):
            ST.decay_box(bad)


def test_decay_box_V3_tg_growth_unchanged_by_efold_box_and_N_doubling():  # V3 (loop 1: 12 e-folds are enough)
    """kc_i with the default rule against 8 and 16 e-folds, twice and four times the box, and N = 160 (observed: ≤ 2e-7
    absolute away from the neutral curve; 1.3e-6 absolute = 2e-4 relative at (0.15, 0.125), 0.0025 below J = k(1 − k))."""
    for k, J in ((0.05, 0.04), (0.05, 0.045), (0.1, 0.08), (0.2, 0.15)):
        ref = ch11.tg_growth(k, J)
        assert ref > 0
        for ne in (8.0, 16.0, 24.0, 48.0):
            assert ch11.tg_growth(k, J, y_max=ST.decay_box(k, n_efold=ne)) == pytest.approx(ref, abs=1e-6), (k, J, ne)
        assert ch11.tg_growth(k, J, N=160, y_max=2 * ST.decay_box(k)) == pytest.approx(ref, abs=1e-6), (k, J)
    assert ch11.tg_growth(0.15, 0.125, N=160) == pytest.approx(ch11.tg_growth(0.15, 0.125), rel=1e-3)  # next to neutral
    for J in (0.05, 0.06):  # outside the tongue at k = 0.05 in every box that is long enough
        assert all(ch11.tg_growth(0.05, J, y_max=ym) == 0.0 for ym in (160.0, 240.0, 480.0))


def test_tg_growth_V7_default_is_the_decay_box_and_explicit_y_max_overrides():  # V7 (loop 1: wiring of the default)
    for k, J in ((0.05, 0.04), (0.1, 0.08), (0.4, 0.1)):
        assert ch11.tg_growth(k, J) == ch11.tg_growth(k, J, y_max=ST.decay_box(k))  # same call, bit for bit
    short = ch11.tg_growth(0.05, 0.04, y_max=30.0)  # the loop-0 defect, reproduced on request only
    assert short == pytest.approx(0.0205913, abs=1e-6) and short > 1.2 * ch11.tg_growth(0.05, 0.04)
    assert ch11.tg_growth(0.05, 0.05, y_max=30.0) > 5e-3  # the spurious mode of the short box
    live = ch11.tg_growth_map(ks=[0.05, 0.4], Js=[0.04, 0.05], cache=False)
    assert live["kci"][0, 0] == ch11.tg_growth(0.05, 0.04) and live["kci"][1, 0] == 0.0
    assert live["kci"][0, 1] == ch11.tg_growth(0.4, 0.04)
    a = ch11.tg_growth_map(ks=[0.05], Js=[0.04], N=100)  # cached twice: the cache key must carry the box rule
    b = ch11.tg_growth_map(ks=[0.05], Js=[0.04], N=100, y_max=30.0)
    assert a["kci"][0, 0] == pytest.approx(0.0163047, abs=1e-6) and b["kci"][0, 0] == pytest.approx(0.0205913, abs=1e-6)


def test_tg_tables_V1_csv_equals_explainer_json_and_small_k_columns_are_live():  # V1 anti-cache (loop 1: E6 table)
    gm = ch11.tg_growth_map()
    tj = json.loads((REF / "explainer_tables.json").read_text(encoding="utf-8"))["tg_map"]
    assert np.array_equal(np.array(tj["kci"], float), gm["kci"])  # the explainer's copy is the CSV, entry for entry
    assert np.allclose(tj["k"], gm["k"], atol=1e-12) and np.allclose(tj["J"], gm["J"], atol=1e-12)
    assert np.allclose(tj["neutral_J"], gm["k"] * (1 - gm["k"]), atol=1e-4)
    KK, JJ = np.meshgrid(gm["k"], gm["J"])
    # Loop 2: the count is derived from the grid, not pinned.  Loop 1 pinned 299 here — that number encoded the defect
    # (36 unstable grid points at k ≥ 0.65 reported as 0).  The exact tongue J < k(1 − k) holds 335 of the 620 grid points
    # (the 1e-12 keeps the points that lie on the curve, e.g. k = 0.4, J = 0.24, out: kc_i = 0 exactly there).
    under = JJ < KK * (1 - KK) - 1e-12
    assert int(under.sum()) == 335 and np.array_equal(gm["kci"] > 0, under)  # growth at every point under the curve, only there
    assert not np.any((gm["kci"] > 0) & (JJ >= KK * (1 - KK)))
    for j in (0, 1, 2):  # k = 0.05, 0.1, 0.15: the columns the short box spoiled — every entry recomputed
        assert gm["k"][j] == pytest.approx(0.05 * (j + 1))
        live = np.array([ch11.tg_growth(gm["k"][j], J) for J in gm["J"]])
        assert np.array_equal(live > 0, gm["kci"][:, j] > 0), j
        assert np.allclose(live, gm["kci"][:, j], rtol=6e-4, atol=0), (j, live, gm["kci"][:, j])  # 4 s.f. rounding


def _tg_shoot_kci(k, J, ci0, L=18.0):
    """Growth rate of the tanh / J sech² layer by a route that shares nothing with the Chebyshev solver: no box, no map,
    no N-filter.  Integrate (11.61), ψ̂″ = [k² + U″/(U − c) − N²/(U − c)²]ψ̂, from z = ∓L inward with the exact far-field
    solutions e^{±kz} (sech²18 ~ 1e-15) and find the purely imaginary c = i c_i (c_r = 0 by the symmetry U(−z) = −U(z)) at
    which the two solutions match at z = 0 (Wronskian = 0), by the secant method from the seed ci0.
    Returns (k c_i, |Im c|, |normalised Wronskian|)."""
    def wronskian(c):
        def rhs(z, y):
            U, s2 = np.tanh(z), 1.0 / np.cosh(z) ** 2
            return [y[1], (k * k - 2.0 * U * s2 / (U - c) - J * s2 / (U - c) ** 2) * y[0]]
        a = solve_ivp(rhs, (-L, 0.0), [1.0 + 0j, k + 0j], method="DOP853", rtol=1e-11, atol=1e-14).y[:, -1]
        b = solve_ivp(rhs, (L, 0.0), [1.0 + 0j, -k + 0j], method="DOP853", rtol=1e-11, atol=1e-14).y[:, -1]
        return (a[0] * b[1] - a[1] * b[0]) / (abs(a[0]) * abs(b[0]) + abs(a[1]) * abs(b[1]))

    x0, x1 = complex(ci0), complex(1.05 * ci0 + 1e-4)
    f0, f1 = wronskian(1j * x0), wronskian(1j * x1)
    for _ in range(60):
        if f1 == f0 or abs(x1 - x0) < 1e-12:
            break
        x0, f0, x1 = x1, f1, x1 - f1 * (x1 - x0) / (f1 - f0)
        f1 = wronskian(1j * x1)
    return k * x1.real, abs(x1.imag), abs(f1)


def test_decay_map_scale_V1_formula_guards_and_limits():  # V1 (loop 2: core.stability.decay_map_scale)
    """s(k) = min(0.5, max(0.035, 0.025/k)): 0.5 for k ≤ 0.05, s·k = 0.025 in between, 0.035 for k ≥ 5/7."""
    for k, s in ((0.01, 0.5), (0.05, 0.5), (0.1, 0.25), (0.2, 0.125), (0.5, 0.05), (0.7, 0.025 / 0.7), (5 / 7, 0.035),
                 (0.9, 0.035), (1.0, 0.035), (3.0, 0.035)):
        assert ST.decay_map_scale(k) == pytest.approx(s, rel=1e-14), k
    ks = np.geomspace(1e-3, 5.0, 80)
    s = np.array([ST.decay_map_scale(k) for k in ks])
    assert np.all((s >= 0.035) & (s <= 0.5)) and np.all(np.diff(s) <= 0)  # bounded; shorter waves never get a coarser centre
    mid = (ks > 0.05) & (ks < 5 / 7)
    assert np.allclose(s[mid] * ks[mid], 0.025, rtol=1e-14)  # between the limits the grid is the same in the wave's own ky
    assert np.allclose([ST.decay_box(k) / ST.decay_map_scale(k) for k in ks[(ks > 0.05) & (ks < 0.4)]], 480.0, rtol=1e-12)
    assert ST.decay_map_scale(0.2, s_max=0.3, s_min=0.01, sk=0.1) == 0.3 and ST.decay_map_scale(2.0, 0.3, 0.01, 0.1) == 0.05
    assert ST.decay_map_scale(20.0, 0.3, 0.01, 0.1) == 0.01
    assert ch11.decay_map_scale is ST.decay_map_scale and np.ndim(ST.decay_map_scale(0.3)) == 0  # one definition, scalar
    for bad in (0.0, -0.3):
        with pytest.raises(ValueError):
            ST.decay_map_scale(bad)
    # what the rule does to the grid: half of the nodes lie inside |y| < s (tan map), so the centre is 14× finer at k = 0.9
    g_new = ST.cheb_grid(100, (-1, 1), "tan", y_max=ST.decay_box(0.9), s=ST.decay_map_scale(0.9))
    g_old = ST.cheb_grid(100, (-1, 1), "tan", y_max=ST.decay_box(0.9), s=0.5)
    assert int(np.sum(np.abs(g_new.y) < 0.05)) >= 2.5 * int(np.sum(np.abs(g_old.y) < 0.05))
    assert g_new.y[0] == pytest.approx(30.0) and g_old.y[0] == pytest.approx(30.0)  # the box is not changed by the scale


def test_tg_growth_V1_shooting_on_the_unbounded_layer_agrees_incl_near_neutral():  # V1 independent route (loop 2, C08)
    """`tg_growth` (Chebyshev, tan map with the k-dependent scale, box 12/k, N-filter) against shooting on the unbounded
    layer.  Observed |Δkc_i| ≤ 1.7e-6 (4 s.f. table); asserted 3e-6 absolute.  The first six points are the ones the
    loop-1 table reported as 0."""
    for k, J in ((0.9, 0.04), (0.9, 0.08), (0.8, 0.13), (0.95, 0.04), (0.65, 0.22), (0.5, 0.24), (0.4, 0.1), (0.4449, 0.0),
                 (0.05, 0.04), (0.15, 0.125)):
        ours = ch11.tg_growth(k, J)
        kci, im, res = _tg_shoot_kci(k, J, ours / k)
        assert ours > 0 and res < 1e-9 and im < 1e-9, (k, J, ours, kci, im, res)
        assert ours == pytest.approx(kci, abs=3e-6), (k, J, ours, kci)


def test_decay_map_scale_V3_tg_growth_unchanged_by_N_doubling_and_scale_halving():  # V3 (loop 2: the rule is converged)
    """Default (N = 100, s = decay_map_scale(k)) against N = 200 with the rule and with half the scale (observed ≤ 1.6e-6
    absolute; asserted 3e-6), next to the neutral curve at large k and at the loop-1 small-k point (0.15, 0.125).  The
    loop-1 fixed scale 0.5 still returns 0 at the large-k points (kept reachable through ``map_scale=0.5``)."""
    for k, J in ((0.9, 0.04), (0.95, 0.04), (0.65, 0.22), (0.15, 0.125)):
        ref = ch11.tg_growth(k, J, N=200)
        assert ref > 4e-3 and ch11.tg_growth(k, J) == pytest.approx(ref, abs=3e-6), (k, J)
        assert ch11.tg_growth(k, J, N=200, map_scale=0.5 * ST.decay_map_scale(k)) == pytest.approx(ref, abs=3e-6), (k, J)
        if k >= 0.65:
            assert ch11.tg_growth(k, J, map_scale=0.5) == 0.0  # the loop-1 defect, on request only


def test_tg_growth_V7_no_growth_at_or_above_the_neutral_curve_off_grid():  # V7 (loop 2: no spurious mode from the finer centre)
    """Nothing at J ≥ k(1 − k) (exact neutral curve), nothing at J ≥ ¼ (Miles–Howard (11.67)), nothing for k ≥ 1 (the
    unstratified layer is neutral at k = 1), at wavenumbers that are not on the table's grid."""
    bad = []
    for k in (0.07, 0.23, 0.5, 0.62, 0.78, 0.905, 0.955, 0.99, 1.0, 1.1):
        nj = max(k * (1 - k), 0.0)
        for J in (nj, nj + 1e-4, nj + 1e-3, nj + 0.01, 0.25, 0.2501, 0.3):
            if J >= nj and ch11.tg_growth(k, J) != 0.0:
                bad.append((k, J, ch11.tg_growth(k, J)))
    assert not bad, bad
    assert ch11.tg_growth(1.1, 0.0) == 0.0 and ch11.tg_growth(0.99, 0.0) > 0  # tanh layer: unstable only for k < 1


def test_tg_growth_V3_near_neutral_sliver_is_as_documented():  # V3 (loop 2: the stated limitation, measured)
    """Documented: a 0 is returned although the layer is unstable only within 0.0002 (k = 0.05) … 0.006 (k = 0.95) of the
    neutral curve, where kc_i < 0.004.  Checked at five k: growth is reported at 1.25 × the documented width below the curve,
    the growth lost at 0.5 × that width (shooting) is below 0.004 and above 0 (so the 0 there is a resolution limit, not
    stability), and a reported growth just outside the sliver equals the shooting value to 3e-6."""
    for k, width in ((0.05, 0.0002), (0.3, 0.0011), (0.5, 0.0021), (0.8, 0.0043), (0.95, 0.0059)):
        nj = k * (1 - k)
        g = ch11.tg_growth(k, nj - 1.25 * width)
        assert 0 < g < 0.006, (k, g)
        assert g == pytest.approx(_tg_shoot_kci(k, nj - 1.25 * width, g / k)[0], abs=3e-6), k
        lost, im, res = _tg_shoot_kci(k, nj - 0.5 * width, 0.5 * g / k)
        assert res < 1e-9 and im < 1e-9 and 0 < lost < 0.004, (k, lost)
    assert "0.006" in ch11.tg_growth.__doc__ and "0.006" in ch11.tg_growth_map.__doc__  # the limitation is stated
    assert "0.006" in (REF / "tg_growth_map.csv").read_text(encoding="utf-8").splitlines()[1]


def test_tg_growth_V7_default_is_the_map_scale_rule_and_explicit_scale_overrides():  # V7 (loop 2: wiring of the default)
    for k, J in ((0.9, 0.04), (0.3, 0.1), (0.05, 0.04)):
        assert ch11.tg_growth(k, J) == ch11.tg_growth(k, J, map_scale=ST.decay_map_scale(k))  # same call, bit for bit
    assert ch11.tg_growth(0.9, 0.04) == pytest.approx(0.0332469, abs=1e-6) and ch11.tg_growth(0.9, 0.04, map_scale=0.5) == 0.0
    live = ch11.tg_growth_map(ks=[0.9, 0.3], Js=[0.04, 0.1], cache=False)
    assert live["kci"][0, 0] == ch11.tg_growth(0.9, 0.04) and live["kci"][1, 1] == ch11.tg_growth(0.3, 0.1)
    assert live["kci"][1, 0] == 0.0  # J = 0.1 > k(1 − k) = 0.09 at k = 0.9
    a = ch11.tg_growth_map(ks=[0.9], Js=[0.04], N=100)  # cached twice: the cache key must carry the map-scale rule
    b = ch11.tg_growth_map(ks=[0.9], Js=[0.04], N=100, map_scale=0.5)
    assert a["kci"][0, 0] == pytest.approx(0.0332469, abs=1e-6) and b["kci"][0, 0] == 0.0
    assert "decay_map_scale" in ch11.tg_growth.__doc__ and "decay_map_scale" in ST.taylor_goldstein_eigs.__doc__


def test_taylor_goldstein_V1_N2_zero_equals_rayleigh_on_the_rule_grid():  # V1 (loop 2: C08 parity survives the new scale)
    pr = ch11.richardson_profiles("tanh", 0.0)
    for k in (0.05, 0.2, 0.7, 0.9, 0.95):
        s, ym = ST.decay_map_scale(k), ST.decay_box(k)
        ctg = ST.taylor_goldstein_eigs(k, pr["U"], pr["Upp"], pr["N2"], domain=(-1, 1), N=100, bc="decay", y_max=ym, map_scale=s)
        cr = ST.rayleigh_eigs(k, TANH["U"], TANH["Upp"], N=100, bc="decay", y_max=ym, map_scale=s, unstable_only=True)
        assert len(ctg) == 1 and abs(ctg[0] - cr[0]) < 1e-8 and abs(ctg[0].real) < 1e-9, (k, ctg, cr)
        assert ch11.tg_growth(k, 0.0) == pytest.approx(k * ctg[0].imag, abs=1e-12)  # tg_growth is this call


def test_parallel_profile_V7_sech2_far_field_is_finite_and_the_clip_is_invisible():  # V7 (loop 2: _sech2 clip at ±350)
    y = np.linspace(-20, 20, 401)
    t, s2 = np.tanh(y), 1.0 / np.cosh(y) ** 2
    exact = dict(shear_layer=(t, -2 * t * s2), jet=(s2, 4 * s2 * t ** 2 - 2 * s2 ** 2))
    far = np.array([-1000.0, -480.0, -356.0, 356.0, 480.0, 1000.0])  # boxes decay_box(k) for k < 0.034 reach these
    for nm, (U, Upp) in exact.items():
        pr = ch11.parallel_profile(nm)
        assert np.allclose(pr["U"](y), U, rtol=1e-13, atol=0) and np.allclose(pr["Upp"](y), Upp, rtol=1e-12, atol=1e-300)
        with np.errstate(over="raise", invalid="raise"):
            u, upp = pr["U"](far), pr["Upp"](far)
        assert np.all(np.isfinite(u)) and np.all(np.isfinite(upp)) and np.all(np.abs(upp) < 1e-300)
    src = (ROOT / "scripts" / "ch11_inviscid_criteria.py").read_text(encoding="utf-8")
    assert "y_max=ch11.decay_box(k)" in src  # the tanh sweep from k = 0.02 runs in the decay box (loop-1 Open item 2)
    c = [ST.rayleigh_eigs(0.02, TANH["U"], TANH["Upp"], bc="decay", N=100, y_max=ym, unstable_only=True)[0].imag
         for ym in (ST.decay_box(0.02), 2 * ST.decay_box(0.02), 30.0)]
    assert c[0] == pytest.approx(c[1], rel=1e-6) and abs(c[2] - c[1]) / c[1] > 0.01  # box 600 converged; box 30 is 1.9 % low


def test_richardson_profiles_V7_far_field_is_finite_and_the_clip_is_invisible():  # V7 (loop 1: cosh clip at ±350)
    z = np.linspace(-20, 20, 401)
    for J, R in ((0.1, 1.0), (0.2, 2.0), (0.05, 0.5)):
        pr = ch11.richardson_profiles("tanh", J, R)
        assert np.allclose(pr["N2"](z), J / np.cosh(R * z) ** 2, rtol=1e-14, atol=0)
        assert np.allclose(pr["Upp"](z), -2 * np.tanh(z) / np.cosh(z) ** 2, rtol=1e-14, atol=0)
        far = np.array([-1000.0, -480.0, -356.0, 356.0, 480.0, 1000.0])
        with np.errstate(over="raise", invalid="raise"):
            n2, upp = pr["N2"](far / R), pr["Upp"](far)
        assert np.all(np.isfinite(n2)) and np.all(np.isfinite(upp)) and np.all(n2 >= 0) and np.all(np.abs(upp) < 1e-300)
        assert np.all(n2 < 1e-300)


def test_decay_callers_V7_scripts_never_use_the_fixed_box_for_long_tg_waves():  # V7 (loop 1: structural check)
    """Every decay-BC Taylor–Goldstein call in scripts/ch11_*.py and in the figure script either passes y_max or has a
    literal k for which the default box 30 already is decay_box(k) (k ≥ 0.4)."""
    import ast

    bad = []
    files = sorted((ROOT / "scripts").glob("ch11_*.py")) + [ROOT / "tests" / "ch11_verify_figures.py"]
    n_calls = 0
    for f in files:
        for node in ast.walk(ast.parse(f.read_text(encoding="utf-8"))):
            if not (isinstance(node, ast.Call) and getattr(node.func, "attr", "") == "taylor_goldstein_eigs"):
                continue
            kw = {k.arg: k.value for k in node.keywords}
            if not (isinstance(kw.get("bc"), ast.Constant) and kw["bc"].value == "decay"):
                continue
            n_calls += 1
            k0 = node.args[0]
            literal_ok = isinstance(k0, ast.Constant) and ST.decay_box(float(k0.value)) == 30.0
            if "y_max" not in kw and not literal_ok:
                bad.append((f.name, node.lineno))
    assert n_calls >= 4 and not bad, bad
    assert ST.taylor_goldstein_eigs.__doc__.count("decay_box") >= 1 and "decay_box" in ch11.tg_growth.__doc__


def test_rayleigh_decay_box_V3_fixed_profile_boxes_are_converged_for_k_from_0_2():  # V3 (E7 table, k ≥ 0.2 rows)
    """The profiles' own boxes, 30 (tanh) and 40 (Bickley), against decay_box(k): converged to 5e-6 relative in c at k = 0.2
    and to 1e-7 for k ≥ 0.3 (observed 3.3e-6, 1.9e-7; 1.6e-8, 1e-10).  At k = 0.1 they are 5e-4 off — since loop 2 the
    Rayleigh table (E7/F6) therefore uses max(profile box, decay_box(k)) (`test_rayleigh_spectra_V1_…` checks its k = 0.1
    and 0.2 rows against twice the box); the fixed boxes remain the defaults of `parallel_profile` for k ≥ 0.2."""
    for nm in ("shear_layer", "jet"):
        pr = ch11.parallel_profile(nm)
        for k, tol in ((0.2, 5e-6), (0.3, 1e-7), (0.6, 1e-7)):
            c = [ST.rayleigh_eigs(k, pr["U"], pr["Upp"], N=120, bc="decay", y_max=ym, map_scale=pr["map_scale"],
                                  unstable_only=True, tol=1e-4)[0] for ym in (pr["y_max"], 2 * ST.decay_box(k))]
            assert abs(c[0] - c[1]) < tol * abs(c[1]), (nm, k, c)


@slow
def test_decay_boxes_V3_other_unbounded_solvers_are_converged_in_y_max():  # V3 (DEVIATION: ∞ truncated at y_max)
    tn = [ch11.tanh_shear_layer_neutral_curve(Re_values=[2.0], cache=False, y_max=ym)["k_upper"][0] for ym in (40.0, 80.0)]
    assert abs(tn[0] - tn[1]) < 1e-8  # Fig. 11.23 (default 40)
    for nm in ("shear_layer", "jet"):
        pr = ch11.parallel_profile(nm)
        c = [ST.rayleigh_eigs(0.1, pr["U"], pr["Upp"], N=120, bc="decay", y_max=ym, map_scale=pr["map_scale"], unstable_only=True,
                              tol=1e-4)[0] for ym in (pr["y_max"], 120.0, 480.0)]
        assert abs(c[0] - c[1]) < 1e-3 * abs(c[1]), (nm, c)  # the profile's fixed box at k = 0.1 (observed ≤ 5e-4)
        assert abs(c[1] - c[2]) < 1e-8 * abs(c[2]), (nm, c)  # loop 2: the table's box decay_box(0.1) = 120 is converged
    b = [ch11.blasius_critical(N=100, y_max=ym, cache=False)["Re_c"] for ym in (20.0, 40.0)]
    assert abs(b[0] - b[1]) < 1e-4 * b[1]  # observed 4e-5 (a linear map on [0, 40] needs N ≈ 100: N = 60 gives 522.0)
    bk = [ch11.bickley_critical(N=60, y_max=ym, cache=False)["Re_c"] for ym in (60.0, 100.0)]
    assert abs(bk[0] - bk[1]) < 1e-5 * bk[1]


def test_howard_semicircle_V1_boundary_and_membership():  # V1 (C10)
    cr, ci = ST.howard_semicircle(-1.0, 3.0, n=101)
    assert np.allclose((cr - 1.0) ** 2 + ci ** 2, 4.0) and np.all(ci >= -1e-15)
    assert ST.in_howard_semicircle(1 + 1.99j, -1, 3) and not ST.in_howard_semicircle(1 + 2.01j, -1, 3)
    assert not ST.in_howard_semicircle(3.2 + 0.01j, -1, 3) and ST.in_howard_semicircle(-0.99 + 0.0j, -1, 3)
    with pytest.raises(ValueError):
        ST.howard_semicircle(1.0, 1.0)


def test_howard_semicircle_V7_every_computed_unstable_eigenvalue_obeys_it():  # V7 theorem check (C10, D19)
    cases = [(TANH, (-1, 1)), (BICK, (0, 1))]
    for pr, (umin, umax) in cases:
        for k in (0.2, 0.5, 0.9):
            c = ST.rayleigh_eigs(k, pr["U"], pr["Upp"], N=100, bc="decay", y_max=30.0, map_scale=1.0, unstable_only=True)
            for v in c:
                assert ST.in_howard_semicircle(v, umin, umax) and k * v.imag <= k * (umax - umin) / 2
    for k, J in ((0.3, 0.0), (0.4, 0.1), (0.6, 0.05)):
        _, res = _tg_mode(k, J)
        assert all(ST.in_howard_semicircle(v, -1, 1) for v in res["c"])
    sinp = ch11.parallel_profile("sin", b=PI)
    c = ST.rayleigh_eigs(0.5, sinp["U"], sinp["Upp"], domain=(-PI, PI), N=80, unstable_only=True, tol=1e-4)
    assert len(c) and all(ST.in_howard_semicircle(v, -1, 1) for v in c)


def test_tg_V2_derivation_D17_taylor_goldstein_from_11_58_to_11_60():  # V2 (D17 ★★, steps 1–9; slip #4)
    x, z, t, k, g, rho0 = sp.symbols("x z t k g rho0", positive=True)
    c = sp.symbols("c")
    U, N2 = sp.Function("U")(z), sp.Function("N2")(z)
    ps, pp, rh = sp.Function("psi")(z), sp.Function("p")(z), sp.Function("rho")(z)
    E = sp.exp(sp.I * k * (x - c * t))
    psi, p, rho = ps * E, pp * E, rh * E
    u, w = sp.diff(psi, z), -sp.diff(psi, x)  # §11.7: u = ∂ψ/∂z, w = −∂ψ/∂x
    xm = sp.diff(u, t) + w * sp.diff(U, z) + U * sp.diff(u, x) + sp.diff(p, x) / rho0  # (11.55) x
    zm = sp.diff(w, t) + U * sp.diff(w, x) + sp.diff(p, z) / rho0 + g * rho / rho0  # (11.55) z, corrected
    de = sp.diff(rho, t) + U * sp.diff(rho, x) - rho0 * N2 * w / g  # (11.56)
    e58 = sp.simplify(xm / (sp.I * k * E))
    assert z0(e58 - ((U - c) * ps.diff(z) - U.diff(z) * ps + pp / rho0))  # step 2 (11.58)
    e59 = sp.simplify(zm / E)
    assert z0(e59 - (k ** 2 * (U - c) * ps + g * rh / rho0 + pp.diff(z) / rho0))  # step 3 (11.59)
    e60 = sp.simplify(de / (sp.I * k * E))
    assert z0(e60 - ((U - c) * rh + rho0 * N2 * ps / g))  # step 4 (11.60)
    p_sol = sp.solve(e58, pp)[0]
    r_sol = sp.solve(e60, rh)[0]  # step 7 (needs U ≠ c)
    tg = sp.simplify(e59.subs(rh, r_sol).subs(pp, p_sol).doit())
    target = (U - c) * (ps.diff(z, 2) - k ** 2 * ps) - U.diff(z, 2) * ps + N2 * ps / (U - c)  # (11.61)
    assert z0(tg + target)  # steps 5–9 (the U′ψ̂′ terms cancel)
    zm_pr = sp.diff(w, t) + U * sp.diff(w, x) + sp.diff(p, x) / rho0 + g * rho / rho0  # slip #4: ∂p/∂x
    tg_pr = sp.simplify((zm_pr / E).subs(rh, r_sol).subs(pp, p_sol).doit())
    assert not z0(tg_pr + target)
    s = ch11.stratified_shear_sympy()
    assert z0(s["tg_residual"]) and all(z0(v) for v in s["normal_modes"]) and s["printed_slip4_reaches_tg"] is False


def test_miles_howard_V2_derivation_D18():  # V2 (D18 ★★★, steps 2–14)
    z, k = sp.symbols("z k")
    c = sp.symbols("c")
    U, N2, phi = sp.Function("U")(z), sp.Function("N2")(z), sp.Function("phi")(z)
    s = sp.sqrt(U - c)
    psi = s * phi  # step 2 (11.63)
    step3 = s * phi.diff(z) + phi * U.diff(z) / (2 * s)
    assert z0(psi.diff(z) - step3)
    step4 = s * phi.diff(z, 2) + (phi.diff(z) * U.diff(z) + phi * U.diff(z, 2) / 2) / s - phi * U.diff(z) ** 2 / (4 * (U - c) ** sp.Rational(3, 2))
    assert z0(psi.diff(z, 2) - step4)
    TG = (U - c) * (psi.diff(z, 2) - k ** 2 * psi) - U.diff(z, 2) * psi + N2 * psi / (U - c)
    step6 = ((U - c) ** sp.Rational(3, 2) * phi.diff(z, 2) + s * phi.diff(z) * U.diff(z) - s * U.diff(z, 2) * phi / 2
             - k ** 2 * (U - c) ** sp.Rational(3, 2) * phi - (U.diff(z) ** 2 / 4 - N2) * phi / s)
    assert z0(TG - step6)  # step 6 ("after some rearrangement", written out)
    step8 = sp.diff((U - c) * phi.diff(z), z) - (k ** 2 * (U - c) + U.diff(z, 2) / 2 + (U.diff(z) ** 2 / 4 - N2) / (U - c)) * phi
    assert z0(TG / s - step8)  # (11.64)
    Ur, cr, ci = sp.symbols("U_r c_r c_i", real=True)
    assert z0(sp.im(1 / (Ur - cr - sp.I * ci)) - ci / ((Ur - cr) ** 2 + ci ** 2))  # step 12
    assert z0(sp.im(Ur - (cr + sp.I * ci)) + ci)
    # step 10 for an actual mode: the identity (11.65) holds and the integrand of step 13 is not positive
    pr, res = _tg_mode(0.4, 0.1)
    cc, ps = res["c"][0], res["psi"][:, 0]
    ri = ch11.richardson_identity_check(0.4, cc, ps, res["y"], pr["U"], pr["Up"], pr["Upp"], pr["N2"], grid=res["grid"])
    assert ri["imag_lhs"] / cc.imag < 0  # step 13: the left integral is negative ⇒ Ri < 1/4 somewhere
    st = ch11.stratified_shear_sympy()
    assert z0(st["self_adjoint_residual"])


def test_howard_V2_derivation_D19():  # V2 (D19 ★★★, steps 1–14; slips #8, #11)
    z, k = sp.symbols("z k")
    c = sp.symbols("c")
    U, N2, F = sp.Function("U")(z), sp.Function("N2")(z), sp.Function("F")(z)
    psi = (U - c) * F  # step 1
    assert z0(psi.diff(z, 2) - ((U - c) * F.diff(z, 2) + 2 * U.diff(z) * F.diff(z) + U.diff(z, 2) * F))  # step 2
    TG = (U - c) * (psi.diff(z, 2) - k ** 2 * psi) - U.diff(z, 2) * psi + N2 * psi / (U - c)
    step4 = (U - c) * ((U - c) * F.diff(z, 2) + 2 * U.diff(z) * F.diff(z) - k ** 2 * (U - c) * F) + N2 * F
    assert z0(TG - step4)  # steps 3–4: the U″ terms cancel
    div = sp.diff((U - c) ** 2 * F.diff(z), z) - k ** 2 * (U - c) ** 2 * F + N2 * F  # step 5
    assert z0(TG - div)
    div_pr = sp.diff((U - c) ** 2 * F.diff(z), z) - k ** 2 * (U - c) * F + N2 * F  # slip #11 as printed
    assert not z0(TG - div_pr)
    Ur, cr, ci = sp.symbols("U_r c_r c_i", real=True)
    assert z0(sp.expand((Ur - cr - sp.I * ci) ** 2) - ((Ur - cr) ** 2 - ci ** 2 - 2 * sp.I * ci * (Ur - cr)))  # step 8
    IU2, IU, IQ, IN = sp.symbols("I_U2 I_U I_Q I_N")
    real_part = IU2 - 2 * cr * IU + (cr ** 2 - ci ** 2) * IQ - IN  # (11.69) written with moments of Q
    assert z0(real_part.subs(IU, cr * IQ) - (IU2 - (cr ** 2 + ci ** 2) * IQ - IN))  # steps 9–11 (11.72)
    Umax, Umin = sp.symbols("U_max U_min", real=True)
    step13 = cr ** 2 + ci ** 2 - cr * (Umax + Umin) + Umax * Umin
    assert z0((cr - (Umax + Umin) / 2) ** 2 + ci ** 2 - ((Umax - Umin) / 2) ** 2 - step13)  # step 14
    st = ch11.stratified_shear_sympy()
    assert z0(st["divergence_residual"]) and not z0(st["printed_slip11_residual"])


# ======================================================================================================================
# C11 — Orr–Sommerfeld (11.79), Squire (11.78)
# ======================================================================================================================
def test_orr_sommerfeld_V5_orszag_eigenvalue():  # V5 (C11, C13; Orszag 1971)
    c = ST.orr_sommerfeld_eigs(1.0, 1e4, POIS["U"], POIS["Upp"], N=100)
    assert abs(c[0] - ORSZAG_C) < 1e-8
    assert ch11.poiseuille_spectrum(1.0, 1e4)[0] == pytest.approx(c[0], abs=1e-12)
    assert np.all(ch11.poiseuille_spectrum(1.0, 1e4)[1:].imag < 0)  # exactly one unstable mode


def test_orr_sommerfeld_V3_spectral_convergence():  # V3 (C11)
    ref = ST.orr_sommerfeld_eigs(1.0, 1e4, POIS["U"], POIS["Upp"], N=140, filter=False)[0]
    errs = [abs(ST.orr_sommerfeld_eigs(1.0, 1e4, POIS["U"], POIS["Upp"], N=N, filter=False)[0] - ref) for N in (30, 40, 50, 60, 80)]
    assert errs[0] > 1e-5 and errs[3] < 1e-7 and errs[4] < 1e-9, errs
    rates = np.diff(np.log(errs[:4])) / 10  # exponential: log error falls linearly in N
    assert np.all(rates < -0.2), rates


def test_orr_sommerfeld_V7_large_Re_approaches_rayleigh_and_guards():  # V7 (N91, N92)
    cr = ST.rayleigh_eigs(0.4449, TANH["U"], TANH["Upp"], N=120, bc="decay", y_max=30.0, map_scale=1.0, unstable_only=True)[0]
    errs = []
    for Re in (1e3, 1e4, 1e5):
        c = ST.orr_sommerfeld_eigs(0.4449, Re, TANH["U"], TANH["Upp"], N=120, bc="decay", y_max=30.0, map_scale=1.0,
                                   filter=False)
        c = c[np.abs(c) < 5]
        errs.append(abs(c[0] - cr))
    assert errs[0] > errs[1] > errs[2] and errs[2] < 1e-4 * 1.1
    full = ST.orr_sommerfeld_eigs(1.0, 1e4, POIS["U"], POIS["Upp"], N=80, filter=False)
    for v in full[np.abs(full.imag) > 0.05][:5]:  # viscosity breaks the c ↔ c* symmetry of Rayleigh's equation
        assert np.min(np.abs(full - v.conjugate())) > 1e-3
    with pytest.raises(ValueError):
        ST.orr_sommerfeld_eigs(0.0, 1e4, POIS["U"], POIS["Upp"])
    with pytest.raises(ValueError):
        ST.orr_sommerfeld_eigs(1.0, -1.0, POIS["U"], POIS["Upp"])


def test_os_mode_V1_normalisation_and_stream_function_convention():  # V1 (N90, R22)
    m = ST.os_mode(1.0, 1e4, POIS["U"], POIS["Upp"], N=100)
    j = int(np.argmax(np.abs(m["u_hat"])))
    assert abs(m["u_hat"][j] - 1.0) < 1e-12 and np.max(np.abs(m["u_hat"])) == pytest.approx(1.0)
    assert np.allclose(m["v_hat"], -1j * m["k"] * m["phi"]) and np.allclose(m["u_hat"], m["grid"].D1 @ m["phi"])
    assert max(abs(m["phi"][0]), abs(m["phi"][-1]), abs(m["u_hat"][0]), abs(m["u_hat"][-1])) < 1e-10  # (11.80)
    g = m["grid"]
    lhs = (np.diag(POIS["U"](g.y)) - m["c"] * np.eye(len(g.y))) @ (g.D2 @ m["phi"] - m["phi"]) - POIS["Upp"](g.y) * m["phi"]
    rhs = (g.D4 @ m["phi"] - 2 * g.D2 @ m["phi"] + m["phi"]) / (1j * 1.0 * 1e4)
    assert np.max(np.abs((lhs - rhs)[5:-5])) < 1e-6  # (11.79) residual at interior nodes
    t = ch11.ts_mode()
    assert abs(t["c"] - ORSZAG_C) < 1e-8 and t["budget"]["ratio"] == pytest.approx(1.6163, abs=1e-4)


def test_squire_V1_transform_and_numerical_3d_equals_2d():  # V1 + V7 (N88, N89; a-D24)
    sq = ST.squire_transform(1.0, math.tan(math.radians(30.0)), 1e4)
    assert sq["Rebar"] == pytest.approx(1e4 * math.cos(math.radians(30.0)), rel=1e-14) and sq["angle_deg"] == pytest.approx(30.0)
    assert sq["kbar"] * sq["Rebar"] == pytest.approx(1.0 * 1e4)  # k̄Re̅ = kRe
    assert sq["Rebar"] <= 1e4 and sq["growth_factor"] >= 1
    s0 = ST.squire_transform(0.7, 0.0, 5000.0)
    assert s0["kbar"] == 0.7 and s0["Rebar"] == 5000.0
    sq2 = ST.squire_transform(1.0, 0.5, 1e4, uhat=np.array([1.0, 2.0]), what=np.array([3.0, -1.0]), phat=np.array([2.0]))
    assert np.allclose(sq2["ubar"], (1.0 * np.array([1, 2]) + 0.5 * np.array([3, -1])) / math.hypot(1, 0.5))
    assert np.allclose(sq2["pbar"], [2.0 * math.hypot(1, 0.5)])
    for k, m, Re in ((1.0, 0.5, 1e4), (0.8, 1.2, 2e4)):
        c3 = ch11.os_3d_eigs(k, m, Re, POIS["U"], POIS["Up"], N=80)
        s = ST.squire_transform(k, m, Re)
        c2 = ST.orr_sommerfeld_eigs(s["kbar"], s["Rebar"], POIS["U"], POIS["Upp"], N=100)
        assert np.min(np.abs(c3 - c2[0])) < 1e-8  # the 2-D twin is in the 3-D spectrum (plus damped Squire modes)
    # the contract slip (Part C wrote Upp for os_3d_eigs): passing U″ is wrong physics and gives a different c
    c3a = ch11.os_3d_eigs(1.0, 0.5, 1e4, POIS["U"], POIS["Up"], N=60)
    assert abs(ch11.os_3d_eigs(1.0, 0.5, 1e4, POIS["U"], POIS["Upp"], N=60)[0] - c3a[0]) > 0.1
    assert abs(ch11.os_3d_eigs(1.0, 0.5, 1e4, POIS["U"], None, N=60)[0] -
               ch11.os_3d_eigs(1.0, 0.5, 1e4, POIS["U"], POIS["Up"], N=60)[0]) < 1e-9
    with pytest.raises(ValueError):
        ST.squire_transform(0.0, 1.0, 1.0)


def test_os_V2_derivation_D21_orr_sommerfeld_from_11_77():  # V2 (D21 ★★, steps 1–10; planted v̂ = +ikφ)
    y = sp.symbols("y", real=True)
    k, Re = sp.symbols("k Re", positive=True)
    c = sp.symbols("c")
    U, phi, pfun = sp.Function("U")(y), sp.Function("phi")(y), sp.Function("p")(y)
    OS = (U - c) * (phi.diff(y, 2) - k ** 2 * phi) - U.diff(y, 2) * phi - (phi.diff(y, 4) - 2 * k ** 2 * phi.diff(y, 2)
                                                                        + k ** 4 * phi) / (sp.I * k * Re)  # (11.79)

    def eliminate(sign):
        uh, vh = phi.diff(y), sign * sp.I * k * phi  # step 1
        xeq = sp.I * k * (U - c) * uh + vh * U.diff(y) + sp.I * k * pfun - (uh.diff(y, 2) - k ** 2 * uh) / Re  # (11.77) x
        yeq = sp.I * k * (U - c) * vh + pfun.diff(y) - (vh.diff(y, 2) - k ** 2 * vh) / Re  # (11.77) y
        assert z0(sp.I * k * uh + vh.diff(y)) if sign == -1 else True  # continuity automatically
        step4 = sp.expand(xeq.diff(y))  # differentiate the x-equation
        step5 = sp.expand(sp.I * k * yeq)  # multiply the y-equation by ik
        diff_ = sp.expand(step4 - step5)  # the pressure cancels
        assert pfun.diff(y) not in diff_.atoms(sp.Derivative)
        return diff_
    good = eliminate(-1)
    assert z0(good / (sp.I * k) - OS)  # steps 6–8
    step6 = sp.I * k * (U - c) * (phi.diff(y, 2) - k ** 2 * phi) - sp.I * k * U.diff(y, 2) * phi - (phi.diff(y, 4) - 2 * k ** 2 * phi.diff(y, 2) + k ** 4 * phi) / Re
    assert z0(good - step6)
    assert not z0(eliminate(+1) / (sp.I * k) - OS)
    r = ch11.os_derivation_sympy()
    assert z0(r["residual"]) and not z0(r["planted_residual"]) and r["factor"] == -1


# ======================================================================================================================
# C12 — Rayleigh's inflection-point theorem, Fjørtoft, critical layers, profiles
# ======================================================================================================================
def test_parallel_profiles_V1_derivatives_are_consistent():  # V1 (every base flow handed to a solver)
    h = 1e-4
    for name, kw, ys in (("poiseuille", {}, np.linspace(-0.9, 0.9, 13)), ("couette", {}, np.linspace(-0.9, 0.9, 13)),
                         ("tanh", {}, np.linspace(-4, 4, 17)), ("bickley", {}, np.linspace(-4, 4, 17)),
                         ("sin", dict(b=2.0), np.linspace(-1.9, 1.9, 13)), ("wall_vorticity_max", {}, np.linspace(-0.9, 0.9, 13)),
                         ("shear_layer_walls", {}, np.linspace(-0.9, 0.9, 13)), ("blasius", {}, np.linspace(0.2, 6, 15)),
                         ("falkner_skan", dict(m=0.1), np.linspace(0.2, 6, 15))):
        p = ch11.parallel_profile(name, **kw)
        U, Up, Upp = p["U"], p["Up"], p["Upp"]
        fd1 = (U(ys + h) - U(ys - h)) / (2 * h)
        fd2 = (Up(ys + h) - Up(ys - h)) / (2 * h)
        assert np.max(np.abs(fd1 - Up(ys))) < 1e-6 * max(1, np.max(np.abs(Up(ys)))), name
        assert np.max(np.abs(fd2 - Upp(ys))) < 2e-5 * max(1, np.max(np.abs(Upp(ys)))), name
    bb = ch11.blasius_base()
    yy = np.linspace(0, 20, 4001)
    assert abs(np.trapezoid(1 - bb["U"](yy), yy) - 1.0) < 1e-5  # δ* units: ∫(1 − U)dy = 1
    assert bb["delta_star_eta"] == pytest.approx(1.7208, abs=1e-4) and abs(bb["U"](20.0) - 1) < 1e-12
    assert ch11.inviscid_profile("jet")["name"] == "bickley" and ch11.inviscid_profile("shear_layer")["name"] == "tanh"
    with pytest.raises(ValueError):
        ch11.parallel_profile("nonsense")


def test_rayleigh_V1_analytic_neutral_modes_tanh_and_bickley():  # V1 (C12, R23, N127; Drazin & Reid)
    y = sp.symbols("y", real=True)
    for U, phi, k, c in ((sp.tanh(y), 1 / sp.cosh(y), 1, 0),  # tanh: k = 1, c = 0, φ = sech y
                         (1 / sp.cosh(y) ** 2, 1 / sp.cosh(y) ** 2, 2, sp.Rational(2, 3)),  # Bickley sinuous
                         (1 / sp.cosh(y) ** 2, sp.tanh(y) / sp.cosh(y), 1, sp.Rational(2, 3))):  # Bickley varicose
        res = (U - c) * (sp.diff(phi, y, 2) - k ** 2 * phi) - sp.diff(U, y, 2) * phi  # (11.81)
        for yv in (sp.Rational(-13, 10), sp.Rational(1, 3), 2):
            assert abs(sp.N(res.subs(y, yv), 40)) < 1e-30  # exact identity at 40 digits
    # numerically the unstable branch ends at those neutral points (c_r → 2/3, c_i → 0)
    for par, k, kn in (("even", 1.95, 2.0), ("odd", 0.95, 1.0)):
        c = ST.rayleigh_eigs(k, BICK["U"], BICK["Upp"], N=160, bc="decay", y_max=40.0, map_scale=1.0, parity=par,
                             unstable_only=False, filter=False)
        c = c[(np.abs(c) < 2) & (c.imag > 1e-3)]
        assert abs(c[0].real - 2 / 3) < 0.01 and 0 < c[0].imag < 0.02
    c = ST.rayleigh_eigs(0.9, TANH["U"], TANH["Upp"], N=120, bc="decay", y_max=30.0, map_scale=1.0, unstable_only=True)
    assert abs(c[0].real) < 1e-10 and 0 < c[0].imag < 0.1  # c_i → 0 as k → 1, c_r = 0


def test_rayleigh_V4_identity_conjugate_pairs_and_V7_no_inflection_no_growth():  # V4 + V7 (C12, N92, D22)
    for k in (0.3, 0.6):
        r = ST.rayleigh_eigs(k, TANH["U"], TANH["Upp"], N=120, bc="decay", y_max=30.0, map_scale=1.0, unstable_only=True,
                             return_vectors=True)
        ck = ch11.rayleigh_identity_check(k, r["c"][0], r["phi"][:, 0], r["y"], TANH["U"], TANH["Upp"], grid=r["grid"])
        assert ck["res_real"] < 1e-9 and ck["res_imag"] < 1e-9
    full = ST.rayleigh_eigs(0.5, BICK["U"], BICK["Upp"], N=100, bc="decay", y_max=40.0, map_scale=1.0, filter=False)
    for v in full[full.imag > 1e-2]:
        assert np.min(np.abs(full - v.conjugate())) < 1e-8
    for name in ("poiseuille", "couette", "blasius"):
        p = ch11.parallel_profile(name)
        for k in (0.3, 1.0):
            kw = dict(domain=p["domain"], bc=p["bc"], y_max=p.get("y_max"))
            assert len(ST.rayleigh_eigs(k, p["U"], p["Upp"], N=80, unstable_only=True, **kw)) == 0, name
    with pytest.raises(ValueError):
        ST.rayleigh_eigs(0.0, TANH["U"], TANH["Upp"])


def test_rayleigh_V2_derivation_D22_inflection_point_theorem():  # V2 (D22 ★★, steps 1–9)
    y = sp.symbols("y", real=True)
    rng = np.random.default_rng(8)
    coef = [sp.Rational(int(a), 3) + sp.I * sp.Rational(int(b), 2) for a, b in rng.integers(-4, 5, (4, 2))]
    phi = (1 - y ** 2) * sum(cc * y ** j for j, cc in enumerate(coef))  # φ = 0 at the walls ±1
    I = lambda f: sp.integrate(sp.expand(f), (y, -1, 1))  # noqa: E731
    assert z0(I(sp.conjugate(phi) * sp.diff(phi, y, 2)) + I(sp.conjugate(sp.diff(phi, y)) * sp.diff(phi, y)))  # step 4
    Ur, cr, ci = sp.symbols("U_r c_r c_i", real=True)
    assert z0(1 / (Ur - cr - sp.I * ci) - (Ur - cr + sp.I * ci) / ((Ur - cr) ** 2 + ci ** 2))  # step 6
    # steps 7–9 on a computed mode: c_i∫U″|φ|²/|U − c|² = 0 with U″ of both signs
    r = ST.rayleigh_eigs(0.4, TANH["U"], TANH["Upp"], N=120, bc="decay", y_max=30.0, map_scale=1.0, unstable_only=True,
                         return_vectors=True)
    g, c, ph = r["grid"], r["c"][0], r["phi"][:, 0]
    wgt = np.abs(ph) ** 2 / np.abs(TANH["U"](g.y) - c) ** 2
    pos = g.integrate(np.maximum(TANH["Upp"](g.y), 0) * wgt)
    neg = g.integrate(np.minimum(TANH["Upp"](g.y), 0) * wgt)
    assert pos > 0 > neg and abs(pos + neg) < 1e-9 * pos


def test_inflection_and_fjortoft_V1_fig_11_21_verdicts():  # V1 + V7 (C12, N93, N94)
    v = {d["panel"]: (d["rayleigh"], d["fjortoft"]) for d in ch11.fig_11_21_verdicts()}
    assert v == {"a": (False, False), "b": (False, False), "c": (False, False), "d": (True, False), "e": (True, True),
                 "f": (True, True)}
    y = np.linspace(-3, 3, 301)
    assert np.allclose(ST.inflection_points(y, Upp=TANH["Upp"]), [0.0], atol=1e-12)
    assert np.allclose(ST.inflection_points(y, U=np.tanh), [0.0], atol=1e-6)
    assert np.allclose(ST.inflection_points(y, U=np.tanh(y)), [0.0], atol=1e-3)
    yI = ST.inflection_points(np.linspace(-4, 4, 401), Upp=BICK["Upp"])
    assert np.allclose(np.abs(yI), math.atanh(1 / math.sqrt(3)), atol=1e-10)  # sech² y: tanh² y = 1/3
    assert ST.inflection_points(np.linspace(-1, 1, 51), Upp=POIS["Upp"]).size == 0
    fj = ch11.fjortoft_criterion(np.linspace(-4, 4, 401), BICK["U"], BICK["Upp"])
    assert fj["satisfied"] and fj["min_product"] < 0
    fw = ch11.fjortoft_criterion(np.linspace(-1, 1, 201), ch11.parallel_profile("wall_vorticity_max")["U"])
    assert ch11.rayleigh_criterion(np.linspace(-1, 1, 201), Upp=ch11.parallel_profile("wall_vorticity_max")["Upp"])["has_inflection"]
    assert not fw["satisfied"]
    with pytest.raises(ValueError):
        ST.inflection_points(y)


def test_sin_profile_V7_inflected_yet_stable_for_2b_below_pi():  # V7 (N95)
    assert ch11.sin_profile_max_growth(1.5, N=80) == 0.0
    vals = [ch11.sin_profile_max_growth(b, N=80) for b in (1.7, 2.0, 3.0)]
    assert 0 < vals[0] < vals[1] < vals[2]
    assert vals[1] == pytest.approx(0.0596, abs=1e-3)


def test_critical_layer_and_cats_eye_V1():  # V1 (N96, N97)
    assert np.allclose(ch11.critical_layer(np.linspace(-1, 1, 41), POIS["U"], 0.5 + 0.01j), [-math.sqrt(0.5), math.sqrt(0.5)])
    U = lambda q: np.sin(q) + 0.3  # noqa: E731
    yc = 0.2
    cc = float(U(yc))
    x = np.linspace(0, 2 * PI, 9)
    errs = []
    for dy in (0.04, 0.02, 0.01):
        yy = np.full_like(x, yc + dy)
        ex = ch11.cats_eye_streamfunction(x, yy, y_c=yc, A=0.01, U=U, exact=True)
        ap = ch11.cats_eye_streamfunction(x, yy, y_c=yc, A=0.01, U=U)
        errs.append(np.max(np.abs(ex - ap)))
    assert abs(observed_order([0.04, 0.02, 0.01], errs) - 3.0) < 0.15  # the dropped term is O((y − y_c)³)
    A_, ph, Uc = 0.1, 1.0, 1.0
    w = ch11.cats_eye_width(A_, ph, Uc)
    assert w == pytest.approx(0.6325, abs=1e-4)
    sep = ch11.cats_eye_streamfunction(0.0, w / 2 * 0 + w, A=A_, phi_c=ph, Uy_c=Uc)  # x = 0 (cos = 1)
    sep2 = ch11.cats_eye_streamfunction(PI, w, A=A_, phi_c=ph, Uy_c=Uc)  # x = π/k: the separatrix level Aφ_c
    assert sep2 == pytest.approx(A_ * ph, rel=1e-12) and sep > sep2


def test_piecewise_shear_layer_V1_neutral_kh_and_small_kh_expansion():  # V1 (N126, Ex. 11.11)
    kh = ch11.piecewise_neutral_kh()
    assert kh == pytest.approx(1.278465, abs=1e-6) and abs((kh - 1) ** 2 - math.exp(-2 * kh)) < 1e-14
    assert abs(ch11.piecewise_shear_layer_c(kh * 1.0001).imag) == 0 and ch11.piecewise_shear_layer_c(kh * 0.999).imag > 0
    for a in (1e-3, 1e-2):
        c = ch11.piecewise_shear_layer_c(a, dU=2.0)
        approx = 1j * (2.0 / 2) * math.sqrt(1 - 4 / 3 * a)
        assert abs(c - approx) < 3 * a ** 2  # Ex. 11.11(d)
    r = minimize_scalar(lambda a: -a * ch11.piecewise_shear_layer_c(a).imag, bounds=(0.3, 1.2), method="bounded",
                        options=dict(xatol=1e-10))
    assert r.x == pytest.approx(0.7968, abs=1e-3) and -r.fun == pytest.approx(0.2012, abs=1e-4)
    mg = ST.max_growth(lambda a: ch11.piecewise_shear_layer_c(a), (0.05, 1.25))  # the sweep tool on an analytic c(k)
    assert mg["k"] == pytest.approx(r.x, abs=1e-6) and mg["growth"] == pytest.approx(-r.fun, rel=1e-10)
    assert mg["c"] == pytest.approx(ch11.piecewise_shear_layer_c(mg["k"]))
    assert ST.max_growth(lambda a: ch11.piecewise_shear_layer_c(a), (0.05, 1.25), measure="ci")["k"] < 0.1  # c_i largest as k → 0


def test_rayleigh_spectra_V1_cached_json_matches_live():  # V1 anti-cache (E7/F6 table)
    js = json.loads((REF / "rayleigh_spectra.json").read_text(encoding="utf-8"))["spectra"]
    live = ch11.rayleigh_spectrum_table(names=["shear_layer", "poiseuille"], ks=[0.4, 0.8], cache=False, write=False)
    for nm in ("shear_layer", "poiseuille"):
        ks = js[nm]["k"]
        for j, kk in enumerate((0.4, 0.8)):
            i = ks.index(kk)
            assert js[nm]["c_i"][i] == pytest.approx(live[nm]["c_i"][j], abs=1e-3 * max(1, abs(live[nm]["c_i"][j])))
    assert all(v == 0 for v in js["poiseuille"]["c_i"]) and all(v == 0 for v in js["couette"]["c_i"])
    assert max(js["shear_layer"]["c_i"]) > 0.4
    # Loop 2: the long-wave rows (k = 0.1, 0.2) of the two unbounded profiles are computed in max(profile box, decay_box(k)).
    # The file (and the explainer's copy) is the 4-s.f. rounding of a live recomputation, and the live value is the
    # box-converged one (twice the box: ≤ 1e-8), which the fixed boxes 30 / 40 missed by 5.1e-4 / 4.4e-4 at k = 0.1.
    live = ch11.rayleigh_spectrum_table(names=["shear_layer", "jet"], ks=[0.1, 0.2], cache=False, write=False)
    ex = json.loads((REF / "explainer_tables.json").read_text(encoding="utf-8"))["rayleigh"]
    r4 = lambda v: float(f"{v:.4g}")  # noqa: E731
    for nm in ("shear_layer", "jet"):
        assert ex[nm] == js[nm]  # the explainer's copy is the file, entry for entry
        pr = ch11.parallel_profile(nm)
        for j, k in enumerate((0.1, 0.2)):
            i = js[nm]["k"].index(k)
            assert js[nm]["c_i"][i] == r4(live[nm]["c_i"][j]), (nm, k)
            if nm == "jet":
                assert js[nm]["c_r"][i] == r4(live[nm]["c_r"][j]), (nm, k)
            big, fixed = (ST.rayleigh_eigs(k, pr["U"], pr["Upp"], N=120, bc="decay", y_max=ym, map_scale=pr["map_scale"],
                                           unstable_only=True, tol=1e-4)[0] for ym in (2 * ST.decay_box(k), pr["y_max"]))
            c_live = complex(live[nm]["c_r"][j], live[nm]["c_i"][j])
            assert abs(c_live - big) < 1e-8 * abs(big), (nm, k, c_live, big)
            if k == 0.1:
                assert 3e-4 < abs(fixed - big) / abs(big) < 7e-4  # what the fixed box cost (loop-1 Open item 2)
    assert js["shear_layer"]["c_i"][0] == 0.8364 and (js["jet"]["c_r"][0], js["jet"]["c_i"][0]) == (0.0928, 0.215)


# ======================================================================================================================
# C13 — plane Poiseuille, Table 11.1, Blasius, Bickley, tanh, Couette, Tollmien
# ======================================================================================================================
def test_poiseuille_critical_V5_orszag_live_and_V3_in_N():  # V5 + V3 (C13; live, never the cache)
    out = {N: ch11.poiseuille_critical(N=N, cache=False) for N in (60, 80)}
    for r in out.values():
        assert rel(r["Re_c"], 5772.22) < 1e-5 and abs(r["k_c"] - 1.02056) < 3e-5 and abs(r["c_r"] - 0.26400) < 1e-5
    assert rel(out[60]["Re_c"], out[80]["Re_c"]) < 1e-6
    cached = ch11.poiseuille_critical()
    assert rel(cached["Re_c"], out[80]["Re_c"]) < 1e-6  # the committed value is the live one
    assert len(ST.rayleigh_eigs(1.0, POIS["U"], POIS["Upp"], N=100, unstable_only=True)) == 0  # inviscidly stable


def test_poiseuille_neutral_curve_V1_live_slice_and_thumb_shape():  # V1 anti-cache (C13, F7)
    fn = lambda k, Re: (lambda c: c[np.abs(c) < 5][:1])(ST.orr_sommerfeld_eigs(k, Re, POIS["U"], POIS["Upp"], N=80, filter=False))  # noqa: E731
    live = ST.neutral_curve(fn, [1e4], (0.6, 1.2), n_k=13)
    assert live["k_lower"][0] == pytest.approx(0.797, abs=2e-3) and live["k_upper"][0] == pytest.approx(1.0947, abs=2e-3)
    for kb, sgn in ((live["k_lower"][0], 1), (live["k_upper"][0], -1)):  # independent row-replacement route brackets them
        assert sgn * _os_rowreplace(kb + 2e-3, 1e4, 100)[0].imag > 0 > sgn * _os_rowreplace(kb - 2e-3, 1e4, 100)[0].imag
    pn = ch11.poiseuille_neutral_curve()
    j = int(np.argmin(np.abs(pn["Re"] - 1e4)))
    kl = ST.neutral_curve(fn, [float(pn["Re"][j])], (0.6, 1.2), n_k=13)
    assert abs(kl["k_lower"][0] - pn["k_lower"][j]) < 1e-6 and abs(kl["k_upper"][0] - pn["k_upper"][j]) < 1e-6
    width = pn["k_upper"] - pn["k_lower"]
    assert np.all(width > 0) and width[-1] < np.max(width)  # the thumb closes as Re → ∞
    cp = ST.critical_point(fn, (0.95, 1.1), (5000.0, 7000.0), k_guess=1.02)
    assert rel(cp["Re_c"], 5772.22) < 1e-5
    with pytest.raises(RuntimeError):
        ST.critical_point(fn, (0.95, 1.1), (100.0, 200.0))


def test_os_tables_V1_cached_csvs_match_live_modes():  # V1 anti-cache (os_grid_*.csv, os_neutral_*.csv)
    g = read_csv("os_grid_poiseuille.csv")
    for i in (100, 350, 550):
        Re, k = g["Re"][i], g["k"][i]
        m = ST.os_mode(float(k), float(Re), POIS["U"], POIS["Upp"], N=80)
        b = ST.disturbance_energy_budget(float(k), m["c"], m["phi"], m["y"], POIS["Up"], float(Re), grid=m["grid"])
        assert g["c_i"][i] == pytest.approx(m["c"].imag, abs=2e-3 * max(abs(m["c"].imag), 1e-3) + 1e-6)
        assert g["P"][i] == pytest.approx(b["production"], rel=2e-3, abs=1e-8)
    nb = read_csv("os_neutral_poiseuille.csv")
    assert np.all(nb["k_lower"][np.isfinite(nb["k_lower"])] < nb["k_upper"][np.isfinite(nb["k_lower"])])
    modes = json.loads((REF / "os_modes.json").read_text(encoding="utf-8"))["modes"]
    assert modes[0]["flow"] == "poiseuille" and abs(complex(*modes[0]["c"]) - ORSZAG_C) < 1e-5


def test_couette_V7_linearly_stable_everywhere():  # V7 (N100)
    assert ch11.couette_max_growth(1.0, 1e3) == pytest.approx(-0.1192, abs=1e-4)
    assert ch11.couette_max_growth(1.0, 1e4) == pytest.approx(-0.0521, abs=1e-4)
    assert ch11.couette_max_growth([0.1, 0.5, 1.0, 3.0], [1e2, 1e4, 1e6]) < 0


def test_blasius_critical_V5_thomas_live_and_V3():  # V5 + V3 (N108; Thomas via Gallagher et al. 2016)
    live = ch11.blasius_critical(N=60, cache=False)
    cached = ch11.blasius_critical()
    for r in (live, cached):
        assert rel(r["Re_c"], 519.2) < 1e-3 and rel(r["k_c"], 0.303) < 5e-3 and rel(r["omega_c"], 0.120) < 5e-3
    assert rel(live["Re_c"], cached["Re_c"]) < 5e-4  # N = 60 vs N = 100 (observed 1.4e-4)
    assert cached["omega_c"] == pytest.approx(cached["k_c"] * cached["c_r"], rel=1e-12)
    bn = ch11.blasius_neutral_curve(in_frequency=True)
    ok = np.isfinite(bn["F_lower"])
    assert np.allclose(bn["F_lower"][ok], bn["k_lower"][ok] * np.real(bn["c_lower"][ok]) / bn["Re"][ok])
    assert np.allclose(bn["lower"][ok], bn["F_lower"][ok]) and np.all(bn["Re"] > live["Re_c"])


def test_bickley_and_tanh_V5_table_11_1_rows_and_V7_tanh_limit():  # V5 (R23; Tatsumi & Kakutani, approx.) + V7 (N99)
    live = ch11.bickley_critical(N=60, cache=False)
    cached = ch11.bickley_critical()
    assert rel(live["Re_c"], 4.0) < 0.01 and rel(cached["Re_c"], 4.0) < 0.01
    assert rel(live["Re_c"], cached["Re_c"]) < 1e-6 and abs(live["k_c"] - cached["k_c"]) < 1e-3
    tn = ch11.tanh_shear_layer_neutral_curve(Re_values=[2.0, 10.0, 200.0], cache=False)
    assert np.allclose(tn["k_upper"], [0.2732, 0.6544, 0.9702], atol=1e-4)
    assert np.all(np.diff(tn["k_upper"]) > 0) and tn["k_upper"][-1] < 1.0  # → 1/L as Re → ∞ (Re_c = 0)
    cached_t = ch11.tanh_shear_layer_neutral_curve()
    assert np.all(cached_t["k_upper"][np.isfinite(cached_t["k_upper"])] < 1.0)
    rows = {r["flow"]: r for r in ch11.table_11_1()}
    assert rows["plane Poiseuille"]["Re_c_ours"] == pytest.approx(5772.22, rel=1e-5)
    assert rows["Blasius"]["Re_c_ours"] == pytest.approx(519.2, rel=1e-3) and rows["plane Couette"]["Re_c_ours"] == math.inf


def test_bickley_V1_half_domain_route_confirms_k_c_below_0_2():  # V1 independent route (R23: our k_c vs the 1958 0.2)
    """Solve the sinuous Bickley OS problem on y ∈ [0, y_max] with symmetry rows φ′(0) = φ‴(0) = 0 (an independent
    discretisation of the even modes) and compare the neutral Re at our k_c with that at the 1958 k ≈ 0.2."""
    def ci_half(k, Re, N=140, ym=80.0):
        g = ST.cheb_grid(N, (0.0, ym), "linear")
        n = len(g.y)
        I = np.eye(n)
        L = g.D2 - k ** 2 * I
        U, Upp = BICK["U"](g.y), BICK["Upp"](g.y)
        A = np.diag(U) @ L - np.diag(Upp) - (g.D4 - 2 * k ** 2 * g.D2 + k ** 4 * I) / (1j * k * Re)
        C = np.vstack([I[0], g.D1[0], g.D1[n - 1], g.D3[n - 1]])  # φ = φ′ = 0 far; φ′ = φ‴ = 0 at the centre
        w = ST.constrained_eig(A, L.astype(complex), C, [0, 1, n - 1, n - 2])
        w = w[np.abs(w) < 2]
        return w[0].imag
    Re_n = lambda k: brentq(lambda R: ci_half(k, R), 2.5, 15.0, xtol=1e-9)  # noqa: E731
    kc = ch11.bickley_critical()["k_c"]
    a, b = Re_n(kc), Re_n(0.2)
    assert a == pytest.approx(ch11.bickley_critical()["Re_c"], rel=1e-4) and a < b  # observed 1.3e-5; Re_n(0.2) = 4.037


def test_tollmien_profile_V1_continuity_and_slip6():  # V1 (N105)
    co = ch11.tollmien_coefficients(0.2)
    assert (co["a"], co["b"]) == pytest.approx((5 / 3, 1 / 0.96))
    for e1 in (0.15, 0.2, 0.3):
        h = 1e-7
        lo, hi = ch11.tollmien_profile(e1 - h, e1), ch11.tollmien_profile(e1 + h, e1)
        assert abs(hi - lo) < 1e-6 and abs(ch11.tollmien_profile(e1 + 2 * h, e1) - hi - (lo - ch11.tollmien_profile(e1 - 2 * h, e1))) < 1e-9
        assert ch11.tollmien_profile(1.0, e1) == pytest.approx(1.0) and ch11.tollmien_profile(2.0, e1) == 1.0
        pr = ch11.tollmien_profile(np.array([e1 - h, e1 + h]), e1, printed=True)
        assert abs(pr[1] - pr[0]) > 0.1  # slip #6: the printed middle branch jumps


# ======================================================================================================================
# C14 — the disturbance-energy equation (11.88)
# ======================================================================================================================
def test_energy_budget_V4_closes_for_computed_modes_and_V7_sign_of_growth():  # V4 + V7 (C14, N104)
    cases = [(POIS, 1.0, 1e4, {}), (POIS, 1.0, 5000.0, {}), (POIS, 0.6, 2e4, {}),
             (TANH, 0.45, 50.0, dict(bc="decay", y_max=30.0, map_scale=1.0)),
             (ch11.parallel_profile("blasius"), 0.25, 1000.0, dict(bc="semi_infinite", y_max=20.0))]
    for pr, k, Re, kw in cases:
        m = ST.os_mode(k, Re, pr["U"], pr["Upp"], domain=pr["domain"], N=100, **kw)
        b = ST.disturbance_energy_budget(k, m["c"], m["phi"], m["y"], pr["Up"], Re, grid=m["grid"])
        assert abs(b["residual"]) < 1e-8 * max(b["production"], b["dissipation"]), (k, Re, b["residual"])
        assert np.sign(b["production"] - b["dissipation"]) == np.sign(m["c"].imag)  # P > Λ iff c_i > 0
        assert b["dissipation"] > 0 and b["E"] > 0
    m = ST.os_mode(1.0, 5000.0, POIS["U"], POIS["Upp"], N=100)
    b = ST.disturbance_energy_budget(1.0, m["c"], m["phi"], m["y"], POIS["Up"], 5000.0, grid=m["grid"])
    assert b["ratio"] == pytest.approx(0.7654, abs=1e-4)  # decays (design C.1 1.22)
    g = read_csv("os_grid_poiseuille.csv")
    ok = np.isfinite(g["c_i"]) & (np.abs(g["c_i"]) > 1e-4)
    assert np.all(np.sign(g["P"][ok] - g["Lambda"][ok]) == np.sign(g["c_i"][ok]))  # whole cached grid
    with pytest.raises(ValueError):
        ST.disturbance_energy_budget(1.0, m["c"], m["phi"][:-1], m["y"][:-1], POIS["Up"], 5000.0)


def test_energy_V2_derivation_D23_and_wavelength_averages():  # V2 (D23 ★★, steps 1–12)
    x, y, t, rho, nu = sp.symbols("x y t rho nu", positive=True)
    psi, p, U = sp.Function("psi")(x, y, t), sp.Function("p")(x, y, t), sp.Function("U")(y)
    uu = [psi.diff(y), -psi.diff(x)]  # ∂_ju_j = 0 by construction
    UU = [U, 0]
    X = [x, y]
    d = lambda f, j: sp.diff(f, X[j])  # noqa: E731
    R = [sp.diff(uu[i], t) + sum(UU[j] * d(uu[i], j) + uu[j] * d(UU[i], j) + uu[j] * d(uu[i], j) for j in range(2))
         + d(p, i) / rho - nu * sum(d(d(uu[i], j), j) for j in range(2)) for i in range(2)]  # step 1 (index form)
    q = sum(ui ** 2 for ui in uu) / 2
    flux = [q * UU[j] + q * uu[j] + p * uu[j] / rho - nu * sum(uu[i] * d(uu[i], j) for i in range(2)) for j in range(2)]
    prod = sum(uu[i] * uu[j] * d(UU[i], j) for i in range(2) for j in range(2))
    diss = nu * sum(d(uu[i], j) ** 2 for i in range(2) for j in range(2))
    step7 = sp.diff(q, t) + sum(d(flux[j], j) for j in range(2)) + prod + diss
    assert z0(sum(uu[i] * R[i] for i in range(2)) - step7)  # steps 2–7
    assert z0(prod - uu[0] * uu[1] * sp.diff(U, y))  # step 12: only uvU′ survives
    k, th = sp.symbols("k theta", positive=True)
    a1, a2, b1, b2 = sp.symbols("a1 a2 b1 b2", real=True)
    uh, vh = a1 + sp.I * a2, b1 + sp.I * b2
    avg = sp.integrate(sp.re(uh * sp.exp(sp.I * th)) * sp.re(vh * sp.exp(sp.I * th)), (th, 0, 2 * sp.pi)) / (2 * sp.pi)
    assert z0(avg - sp.re(uh * sp.conjugate(vh)) / 2)  # ⟨uv⟩ = ½Re(ûv̂*)
    assert z0(sp.integrate(sp.diff(sp.cos(k * x) ** 3 * sp.sin(k * x), x), (x, 0, 2 * sp.pi / k)))  # step 10 (P273)
    assert z0(ch11.energy_equation_sympy()["residual"])
    m = ST.os_mode(1.0, 1e4, POIS["U"], POIS["Upp"], N=100)  # the closed form of ts_wave_fields vs a numerical x-average
    xs = np.linspace(0, 2 * PI, 400, endpoint=False)
    f = ch11.ts_wave_fields(xs, m["y"], 0.0, 1.0, m["c"], m["phi"], m["u_hat"], amp=0.05, U=POIS["U"])
    assert np.allclose((f["u"] * f["v"]).mean(axis=1), f["uv_mean"], atol=1e-14)
    assert np.allclose(f["u_total"], POIS["U"](m["y"])[:, None] + f["u"])
    h = xs[1] - xs[0]
    assert np.max(np.abs(-np.gradient(f["psi"], h, axis=1) - f["v"])[:, 3:-3]) < 2e-4 * np.max(np.abs(f["v"]))  # v = −∂ψ/∂x


# ======================================================================================================================
# C15 — chaos: pendulum, Hopf, Lorenz (11.91), logistic map
# ======================================================================================================================
def test_lorenz_V1_fixed_points_eigenvalues_and_V5_hopf_point():  # V1 + V5 (C15, N116; Wikipedia "Lorenz system")
    for r in (0.5, 10.0, 28.0):
        for pnt in ch11.lorenz_fixed_points(r):
            assert np.max(np.abs(ch11.lorenz_rhs(0, pnt, r=r))) < 1e-12
    assert len(ch11.lorenz_fixed_points(0.5)) == 1
    c = ch11.lorenz_fixed_points(28.0)[1]
    assert c == pytest.approx((8.4853, 8.4853, 27.0), abs=1e-4)
    e = ch11.lorenz_eigs(28.0)
    assert e[0].real == pytest.approx(0.0940, abs=1e-4) and abs(e[0].imag) == pytest.approx(10.1945, abs=1e-4)
    assert e[2].real == pytest.approx(-13.8546, abs=1e-4)
    eo = ch11.lorenz_eigs(28.0, which="O")
    assert eo.real == pytest.approx([11.8277, -2.6667, -22.8277], abs=1e-4)
    rH = ch11.lorenz_hopf_r()
    assert rH == pytest.approx(24.7368, abs=1e-4) and rel(rH, 24.74) < 2e-4  # Wikipedia: 24.74
    # independent of the formula: the real part of the complex pair crosses zero at r_H
    re_ = lambda r: ch11.lorenz_eigs(r)[0].real  # noqa: E731
    assert brentq(re_, 20.0, 28.0, xtol=1e-12) == pytest.approx(rH, abs=1e-8)
    assert re_(20.0) < 0 < re_(28.0) and ch11.lorenz_eigs(20.0)[0].real == pytest.approx(-0.155, abs=1e-3)
    with pytest.raises(ValueError):
        ch11.lorenz_hopf_r(Pr=3.0)
    with pytest.raises(ValueError):
        ch11.lorenz_eigs(0.5)


def test_lorenz_V4_volume_contraction_and_V3_rk4_order():  # V4 + V3 (C15, E9 parity)
    rng = np.random.default_rng(9)
    for s in rng.normal(0, 10, (5, 3)):
        assert np.trace(ch11.lorenz_jacobian(s)) == pytest.approx(ch11.lorenz_divergence(), rel=1e-14)
    assert ch11.lorenz_divergence() == pytest.approx(-(10 + 1 + 8 / 3))
    h = 1e-6
    s0 = np.array([1.0, 2.0, 3.0])
    Jfd = np.array([(np.array(ch11.lorenz_rhs(0, s0 + h * e)) - np.array(ch11.lorenz_rhs(0, s0 - h * e))) / (2 * h)
                    for e in np.eye(3)]).T
    assert np.allclose(Jfd, ch11.lorenz_jacobian(s0), atol=1e-8)
    assert np.allclose(ch11.lorenz_rhs(0, (1, 1, 1)), (0, 26, -5 / 3))
    one = ch11.lorenz_rk4((1, 1, 1), 0.01, 1)[1]
    assert one == pytest.approx((1.012567, 1.259918, 0.984891), abs=1e-6)
    ref = solve_ivp(ch11.lorenz_rhs, (0, 0.5), [1, 1, 1], method="DOP853", rtol=1e-13, atol=1e-14).y[:, -1]
    dts = [0.02, 0.01, 0.005]
    errs = [np.linalg.norm(ch11.lorenz_rk4((1, 1, 1), dt, int(round(0.5 / dt)))[-1] - ref) for dt in dts]
    assert abs(observed_order(dts, errs) - 4.0) < 0.15, errs
    li = ch11.lorenz_integrate(t_end=5.0, n=501)
    rk = ch11.lorenz_rk4((1, 1, 1), 0.005, 1000)
    assert np.allclose(rk[-1], [li["X"][-1], li["Y"][-1], li["Z"][-1]], atol=1e-6)
    long = ch11.lorenz_integrate(t_end=40.0, fast=True)
    assert long["t"][-1] == 25.0 and np.max(np.abs(long["X"])) < 25 and np.min(long["Z"]) > -1e-9  # bounded attractor


def test_lorenz_V7_sensitivity_and_regimes():  # V7 (N112, N116; Lyapunov qualitative)
    sep = ch11.lorenz_separation()
    assert 0.8 < sep["slope"] < 1.0  # ≈ 0.906 (literature, qualitative)
    T = ch11.lorenz_predictability_time()
    assert 26.0 < T < 32.0 and np.isnan(ch11.lorenz_predictability_time(threshold=1e3))
    lam = ch11.lorenz_largest_lyapunov(t_end=60.0, cache=False)
    assert 0.7 < lam < 1.2
    sw = ch11.lorenz_r_sweep(r_values=[0.5, 15.0], t_end=30.0, n=600, cache=False)
    assert abs(sw["X"][0][-1]) < 1e-3  # r < 1: conduction (origin attracts)
    assert abs(abs(sw["X"][1][-1]) - math.sqrt(8 / 3 * 14)) < 1e-2  # 1 < r < r_H: steady rolls C±
    cached = ch11.lorenz_r_sweep()
    j = list(cached["r"]).index(28.0)
    assert np.ptp(cached["X"][j][len(cached["t"]) // 2:]) > 20  # r = 28: switching between the two lobes
    assert ch11.lorenz_b(PI / math.sqrt(2)) == pytest.approx(8 / 3, rel=1e-15)
    for k in (1.0, PI / math.sqrt(2), 4.0):
        assert ch11.lorenz_r(ch11.benard_free_free_Ra(k), k) == pytest.approx(1.0, rel=1e-14)  # r = Ra/Ra_c(k)


def test_lorenz_fields_V1_truncation_kinematics():  # V1 (N115; §11.14 sign u = −∂ψ/∂z, w = ∂ψ/∂x)
    x, z = np.meshgrid(np.linspace(0, 4, 401), np.linspace(-0.5, 0.5, 161))
    f = ch11.lorenz_fields(x, z, 2.0, 1.5, -0.7)
    hx, hz = x[0, 1] - x[0, 0], z[1, 0] - z[0, 0]
    sl = (slice(2, -2), slice(2, -2))
    assert np.max(np.abs((-np.gradient(f["psi"], hz, axis=0) - f["u"])[sl])) < 2e-3
    assert np.max(np.abs((np.gradient(f["psi"], hx, axis=1) - f["w"])[sl])) < 2e-3
    div = np.gradient(f["u"], hx, axis=1) + np.gradient(f["w"], hz, axis=0)
    assert np.max(np.abs(div[sl])) < 1e-2
    assert np.allclose(f["w"][[0, -1]], 0, atol=1e-12) and np.allclose(f["T"][[0, -1]], 0, atol=1e-12)  # free–free walls


def test_lorenz_V2_derivation_D24_galerkin_projection():  # V2 (D24 ★★★, steps 1–15)
    x, z, t = sp.symbols("x z t", real=True)
    k, Pr, Ra = sp.symbols("k Pr Ra", positive=True)
    A, B, C = [sp.Function(n)(t) for n in ("A", "B", "C")]
    psi = A * sp.cos(sp.pi * z) * sp.sin(k * x)  # (11.90)
    th = B * sp.cos(sp.pi * z) * sp.cos(k * x) + C * sp.sin(2 * sp.pi * z)
    J = lambda f, g: sp.diff(f, x) * sp.diff(g, z) - sp.diff(f, z) * sp.diff(g, x)  # noqa: E731  step 3
    lap = lambda f: sp.diff(f, x, 2) + sp.diff(f, z, 2)  # noqa: E731
    a2 = sp.pi ** 2 + k ** 2
    assert z0(lap(psi) + a2 * psi) and z0(J(psi, lap(psi)))  # step 4
    step8 = -sp.pi * k / 2 * A * B * sp.sin(2 * sp.pi * z) + 2 * sp.pi * k * A * C * sp.cos(sp.pi * z) * sp.cos(2 * sp.pi * z) * sp.cos(k * x)
    assert z0((J(psi, th) - step8).rewrite(sp.exp))  # step 8
    vort = (sp.diff(lap(psi), t) + J(psi, lap(psi))) / Pr - Ra * sp.diff(th, x) - lap(lap(psi))  # step 5
    heat = sp.diff(th, t) + J(psi, th) - sp.diff(psi, x) - lap(th)
    proj = lambda e, s: sp.integrate(sp.integrate(sp.expand(e * s), (x, 0, 2 * sp.pi / k)), (z, -sp.Rational(1, 2), sp.Rational(1, 2)))  # noqa: E731
    dA = sp.solve(proj(vort, sp.cos(sp.pi * z) * sp.sin(k * x)), sp.diff(A, t))[0]
    dB = sp.solve(proj(heat, sp.cos(sp.pi * z) * sp.cos(k * x)), sp.diff(B, t))[0]
    dC = sp.solve(proj(heat, sp.sin(2 * sp.pi * z)), sp.diff(C, t))[0]
    assert z0(dA - Pr * (Ra * k / a2 * B - a2 * A))  # step 7
    assert z0(dB - (-a2 * B + k * A - sp.pi * k * A * C))  # step 9
    assert z0(dC - (sp.pi * k / 2 * A * B - 4 * sp.pi ** 2 * C))  # step 10
    X, Y, Z, r = sp.symbols("X Y Z r", positive=True)
    aL = sp.pi * k / (sp.sqrt(2) * a2)  # step 13
    rr = Ra * k ** 2 / a2 ** 3  # step 14
    bL = aL * Ra * k / a2 ** 2
    gL = sp.pi * rr
    assert z0(bL * k / (aL * a2) - rr)
    sub = {A: X / aL, B: Y / bL, C: Z / gL}
    dX = sp.simplify(aL * dA.subs(sub) / a2)  # step 11: τ = a²t
    dY = sp.simplify(bL * dB.subs(sub) / a2)
    dZ = sp.simplify(gL * dC.subs(sub) / a2)
    b = 4 * sp.pi ** 2 / a2
    assert z0(dX - Pr * (Y - X)) and z0(dY - (-X * Z + rr * X - Y)) and z0(dZ - (X * Y - b * Z))  # step 15 (11.91)
    assert b.subs(k, sp.pi / sp.sqrt(2)) == sp.Rational(8, 3)
    e = ch11.lorenz_sympy()
    assert all(z0(v) for v in e["scaled_residuals"]) and z0(e["b_at_kc"] - sp.Rational(8, 3))


def test_lorenz_V2_derivation_D25_fixed_points_and_hopf():  # V2 (D25 ★★, steps 1–11)
    X, Y, Z = sp.symbols("X Y Z", real=True)
    r, Pr, b, lam = sp.symbols("r Pr b lambda")
    rhs = [Pr * (Y - X), -X * Z + r * X - Y, X * Y - b * Z]
    sols = sp.solve(rhs, [X, Y, Z], dict=True)
    assert any(z0(s[X] ** 2 - b * (r - 1)) and z0(s[Z] - (r - 1)) and z0(s[Y] - s[X]) for s in sols)  # steps 1–4
    assert any(s[X] == 0 and s[Y] == 0 and s[Z] == 0 for s in sols) and len(sols) == 3
    Jm = sp.Matrix(rhs).jacobian([X, Y, Z])
    assert Jm == sp.Matrix([[-Pr, Pr, 0], [r - Z, -1, -X], [Y, X, -b]])  # step 5
    p0 = sp.factor((lam * sp.eye(3) - Jm.subs({X: 0, Y: 0, Z: 0})).det())
    assert z0(p0 - (lam + b) * (lam ** 2 + (Pr + 1) * lam + Pr * (1 - r)))  # step 6
    Xb = sp.sqrt(b * (r - 1))
    pc = sp.expand((lam * sp.eye(3) - Jm.subs({X: Xb, Y: Xb, Z: r - 1})).det())
    assert z0(pc - (lam ** 3 + (Pr + b + 1) * lam ** 2 + b * (r + Pr) * lam + 2 * b * Pr * (r - 1)))  # step 9
    rH = sp.solve(sp.Eq((Pr + b + 1) * b * (r + Pr), 2 * b * Pr * (r - 1)), r)  # steps 10–11
    assert len(rH) == 1 and z0(rH[0] - Pr * (Pr + b + 3) / (Pr - b - 1))
    w = sp.symbols("omega", positive=True)
    a2, a1, a0 = sp.symbols("a2 a1 a0", positive=True)
    cub = (sp.I * w) ** 3 + a2 * (sp.I * w) ** 2 + a1 * sp.I * w + a0
    assert z0(sp.re(sp.expand(cub)).subs(w, sp.sqrt(a1)).subs(a0, a2 * a1)) and z0(sp.im(sp.expand(cub)).subs(w, sp.sqrt(a1)))


def test_pendulum_hopf_V4_energy_and_V1_limit_cycle():  # V4 + V1 (N113, N114)
    tr = ch11.phase_portrait(ch11.pendulum_rhs, [(1.0, 0.0), (2.5, 0.3), (0.2, 2.2)], t_end=20.0, n=400)
    for d in tr:
        E = ch11.pendulum_energy((d["X"], d["Y"]))
        assert np.ptp(E) < 1e-9
    damped = ch11.phase_portrait(ch11.pendulum_rhs, [(1.0, 0.0)], t_end=20.0, damping=0.2)[0]
    Ed = ch11.pendulum_energy((damped["X"], damped["Y"]))
    assert np.all(np.diff(Ed) <= 1e-10)
    for mu in (0.04, 0.25):
        hp = ch11.phase_portrait(ch11.hopf_normal_form, [(0.05, 0.0), (1.5, 0.0)], t_end=200.0, n=200, mu=mu)
        for d in hp:
            assert math.hypot(d["X"][-1], d["Y"][-1]) == pytest.approx(ch11.limit_cycle_amplitude(mu), rel=1e-6)
    assert ch11.limit_cycle_amplitude(-0.3) == 0.0
    hp = ch11.phase_portrait(ch11.hopf_normal_form, [(0.5, 0.0)], t_end=60.0, n=50, mu=-0.2)[0]
    assert math.hypot(hp["X"][-1], hp["Y"][-1]) < 1e-4


def test_logistic_V1_fixed_points_and_V5_feigenbaum():  # V1 + V5 (N117, N129; Wikipedia "Feigenbaum constants")
    for A in (1.5, 2.8, 3.3):
        fp = ch11.logistic_fixed_point(A)
        xs = fp["x_star"]
        assert A * xs * (1 - xs) == pytest.approx(xs) and fp["multiplier"] == pytest.approx(2 - A)
        assert fp["stable"] is (1 < A < 3)
    x = ch11.logistic_map(2.8, 0.3, 400)
    assert x[0] == 0.3 and abs(x[-1] - (1 - 1 / 2.8)) < 1e-12
    cw = ch11.cobweb(2.8, 0.3, 5)
    assert len(cw["xs"]) == 11 and cw["ys"][0] == 0 and cw["xs"][1] == 0.3 and cw["ys"][1] == pytest.approx(x[1])
    bd = ch11.bifurcation_diagram([2.8, 3.2], n_transient=2000, n_keep=50)
    assert len(np.unique(np.round(bd["x"][bd["A"] == 2.8], 8))) == 1
    assert len(np.unique(np.round(bd["x"][bd["A"] == 3.2], 8))) == 2
    pd_ = ch11.period_doubling_points(5)
    An = pd_["A_n"]
    assert An[0] == 3.0 and An[1] == pytest.approx(1 + math.sqrt(6), abs=1e-12)
    assert An[2] == pytest.approx(3.5440903, abs=2e-7) and An[3] == pytest.approx(3.5644073, abs=2e-7)
    assert An[4] == pytest.approx(3.5687594, abs=2e-7)
    S = ch11.superstable_points(6)
    assert S[1] == pytest.approx(1 + math.sqrt(5)) and np.all(np.diff(S) > 0)
    for n_, Sn in enumerate(S):  # x = ½ lies on the superstable 2ⁿ cycle
        assert abs(ch11._f_iter(Sn, 0.5, 2 ** n_) - 0.5) < 1e-10
    est = ch11.feigenbaum_estimate(8)
    assert abs(est[-1] - 4.669201609) < 1e-3 and abs(est[-1] - 4.669201609) < abs(est[0] - 4.669201609)
    assert est[:5] == pytest.approx([4.709, 4.681, 4.663, 4.668, 4.669], abs=1e-3)


# ======================================================================================================================
# Slips, engines, contract
# ======================================================================================================================
def test_book_slips_V7_every_coded_slip_fails_its_check():  # V7 (analysis §9 S1–S12)
    rows = {r["key"]: r for r in ch11.book_slips()}
    assert sorted(rows, key=lambda s: int(s[1:])) == [f"S{i}" for i in range(1, 13)]
    assert rel(ch11.benard_marginal_Ra_det(2.0, printed=True), ch11.benard_marginal_Ra_det(2.0)) > 0.01  # S1
    assert abs(ch11.benard_free_free_mode(0.5, printed=True)) > 0.5  # S2
    assert ch11.benard_free_free_sympy()["printed_root"] == []  # S3
    st = ch11.stratified_shear_sympy(printed=True)
    assert st["printed_slip4_reaches_tg"] is False and not z0(st["residual"][1])  # S4, S11
    assert abs(ch11.tollmien_profile(0.2 + 1e-9, printed=True) - ch11.tollmien_profile(0.2 - 1e-9, printed=True)) > 0.1  # S6
    assert math.isnan(ch11.taylor_galerkin_Ta(3.12, 0.0, printed=True))  # S7
    assert ch11.taylor_perturbation_sympy()["printed_continuity_units_ok"] is False  # S12
    for r in rows.values():
        assert r["printed"] and r["correct"] and r["where"]


def test_derive_all_V2_every_engine_reports_ok():  # V2 (engine roll-call)
    d = ch11.derive_all()
    assert set(d) == {"kh_sympy", "kh_depth_tension_sympy", "benard_perturbation_sympy", "exchange_of_stabilities_sympy",
                      "benard_free_free_sympy", "taylor_perturbation_sympy", "stratified_shear_sympy", "os_derivation_sympy",
                      "energy_equation_sympy", "lorenz_sympy"}
    assert all(v["ok"] for v in d.values()), {k: v for k, v in d.items() if not v["ok"]}


def test_reference_files_V7_benchmarks_json_matches_make_refs_and_sources():  # V7 (reference/ch11 consistency)
    bj = json.loads((REF / "benchmarks.json").read_text(encoding="utf-8"))
    assert bj["plane_poiseuille_critical"]["Re_c"] == 5772.22 and bj["plane_poiseuille_Re1e4_k1"]["c_i"] == 0.00373967
    assert bj["benard_rigid_rigid"]["Ra_c"] == 1707.762 and bj["benard_free_free"]["Ra_c"] == 27 * PI ** 4 / 4
    assert bj["lorenz"]["r_H"] == pytest.approx(ch11.lorenz_hopf_r(), rel=1e-15)
    assert bj["feigenbaum"]["A_2"] == pytest.approx(1 + math.sqrt(6))
    src = (REF / "SOURCES.md").read_text(encoding="utf-8")
    for key in ("Orszag", "Chandrasekhar", "Gallagher", "Tatsumi", "Michalke", "Feigenbaum", "Lorenz"):
        assert key in src
    ex = json.loads((REF / "explainer_tables.json").read_text(encoding="utf-8"))
    bt = read_csv("benard_neutral_curves.csv")
    assert ex["benard"]["rigid"][20] == pytest.approx(bt["Ra_rigid_rigid"][20], rel=1e-3)
    assert ex["taylor"]["eigenfunctions"]["0"]["Ta"] == pytest.approx(3389.90, abs=0.01)


PART_C_EXTRA = ["parallel_profile", "fig_11_21_verdicts", "ts_mode", "tg_tanh_neutral_J", "tg_tanh_neutral_mode",
                "superstable_points", "logistic_fixed_point", "tollmien_coefficients", "lorenz_divergence"]


def test_part_c_V7_every_contract_function_exists_and_is_exercised():  # V7 (design Part C: 149 names)
    design = (ROOT / "analysis" / "ch11_design.md").read_text(encoding="utf-8").splitlines()
    a = next(i for i, ln in enumerate(design) if ln.startswith("### C.1"))
    b = next(i for i, ln in enumerate(design) if ln.startswith("### C.3"))
    names = []
    for ln in design[a:b]:
        if ln.startswith("| ") and not ln.startswith("| #") and not ln.startswith("|---"):
            for m in re.finditer(r"`([A-Za-z_][A-Za-z0-9_]*)\(", ln.split("|")[2]):
                if m.group(1) not in names:
                    names.append(m.group(1))
    assert len(names) >= 140
    src = Path(__file__).read_text(encoding="utf-8")
    missing = [n for n in names if not hasattr(ch11, n)]
    untested = [n for n in names + PART_C_EXTRA if not re.search(r"\b(?:ch11|ST)\." + re.escape(n) + r"\b", src)]
    assert not missing and not untested, (missing, untested)
    for n in ("sigma_from_c", "kh_phase_speed", "benard_free_free_Ra", "benard_free_free_sigma", "lorenz_b", "lorenz_r",
              "taylor_critical_approx", "rayleigh_number", "gradient_richardson", "piecewise_shear_layer_c", "tollmien_profile",
              "cats_eye_width"):
        args = {"sigma_from_c": (1.0, 0.5), "kh_phase_speed": (1.0, 1.0, 0.0, 1.0, 2.0), "benard_free_free_Ra": (2.0,),
                "benard_free_free_sigma": (2.0, 1000.0, 1.0), "lorenz_b": (1.0,), "lorenz_r": (1000.0, 2.0),
                "taylor_critical_approx": (0.5,), "rayleigh_number": (2e-4, 1.0, 0.01, 1.4e-7, 1e-6),
                "gradient_richardson": (0.3, None, 0.1, 0.5), "piecewise_shear_layer_c": (0.5,), "tollmien_profile": (0.3,),
                "cats_eye_width": (0.1, 1.0, 1.0)}[n]
        v = getattr(ch11, n)(*args)  # scalar-callable for the explainers' selftest parity rows
        vv = v if isinstance(v, tuple) else (v,)
        assert all(np.ndim(e) == 0 for e in vv), n


# ======================================================================================================================
# V6 — numbers printed by the book (private; skipped when the JSON is absent)
# ======================================================================================================================
@book_only
def test_book_V6_benard_taylor_and_double_diffusion_numbers():  # V6 (C04, C05, C06, C07, N50, N123)
    bv = book()
    b4 = bv["sec_11_4_benard"]
    rr = ch11.benard_critical()
    assert rel(rr["Ra_c"], b4["Ra_cr_rigid_rigid"]) < 5e-3 and rel(rr["K_c"], b4["K_cr_rigid_rigid"]) < 5e-3
    assert rel(2 * PI / rr["K_c"], b4["lambda_cr_over_d"]) < 0.01
    assert rel(ch11.RA_FREE_FREE, b4["Ra_cr_free_free"]) < 5e-3
    rf = ch11.benard_critical(bc=("rigid", "free"), mode="any")
    assert rel(rf["Ra_c"], b4["Ra_cr_rigid_free"]) < 5e-3 and rel(rf["K_c"], b4["K_cr_rigid_free"]) < 5e-3
    ex = bv["exercises"]
    od = ch11.benard_critical(mode="odd")
    assert rel(od["Ra_c"], ex["ex_11_7_odd_mode_Ra"]) < 5e-3 and rel(od["K_c"], ex["ex_11_7_odd_mode_K"]) < 5e-3
    b6 = bv["sec_11_6_taylor"]
    assert rel(ch11.taylor_critical(1.0)["Ta_c"], b6["Ta_cr_numerator_11_54"]) < 5e-3
    assert rel(ch11.taylor_critical(1.0)["k_c"], b6["k_cr"]) < 5e-3
    assert rel(ch11.taylor_critical(0.0)["k_c"], ex["ex_11_9_k_cr"]) < 5e-3
    sb = ch11.taylor_stability_boundary(b6["fig_11_17_radius_ratio_R2_over_R1"], mus=[-0.5, 0.0])
    assert np.all(sb["y_inner"] > 0)
    b5 = bv["sec_11_5_double_diffusion"]
    assert rel(ch11.RA_FREE_FREE, b5["Rs_minus_Ra_critical"]) < 5e-3
    b3 = bv["sec_11_3_kelvin_helmholtz"]
    m = ch11.kh_mixing_energy(2.0, 1.0, 1.0)
    assert m["E_i"] / 4.0 == pytest.approx(b3["E_initial_over_rho_U1sq_h"]) and m["E_f"] / 4.0 == pytest.approx(b3["E_final_over_rho_U1sq_h"])


@book_only
def test_book_V6_viscous_table_11_1_tollmien_stratified_and_chaos():  # V6 (C13, N105, N78, C15)
    bv = book()
    t = bv["sec_11_10_viscous"]["table_11_1"]
    rows = {r["flow"]: r for r in ch11.table_11_1()}
    assert rel(rows["jet (Bickley)"]["Re_c_ours"], t["jet_sech2_Re_cr"]) < 5e-3
    assert rel(rows["Blasius"]["Re_c_ours"], t["blasius_Re_cr_delta_star"]) < 5e-3
    assert rel(rows["plane Poiseuille"]["Re_c_ours"], t["plane_poiseuille_Re_cr"]) < 5e-3  # the book rounds 5772.22 up
    assert rows["shear layer"]["Re_c_ours"] == t["shear_layer_tanh_Re_cr"]
    assert t["pipe_Re_cr"] == "infinity" and rows["pipe"]["Re_c_ours"] == math.inf and rows["plane Couette"]["Re_c_ours"] == math.inf
    tp = bv["sec_11_11_tollmien_schlichting"]["tollmien_profile_printed"]
    a_ = float(re.search(r"([\d.]+)\s*\(y/delta\)", tp["inner"]).group(1))
    e1 = float(re.search(r"<=\s*([\d.]+)", tp["inner"]).group(1))
    b_ = float(re.search(r"1 - ([\d.]+)", tp["outer"]).group(1))
    jump = abs(ch11.tollmien_profile(e1 + 1e-12, e1, a=a_, b=b_) - ch11.tollmien_profile(e1 - 1e-12, e1, a=a_, b=b_))
    assert jump / (a_ * e1) < 5e-3  # corrected form continuous to the book's rounding
    jp = abs(ch11.tollmien_profile(e1 + 1e-12, e1, printed=True, a=a_, b=b_) - ch11.tollmien_profile(e1 - 1e-12, e1, printed=True, a=a_, b=b_))
    assert jp > 0.2
    s7 = bv["sec_11_7_stratified"]
    assert ch11.miles_howard_stable(np.linspace(-3, 3, 61), dUdz=TANH["Up"], N2=lambda z: 1.01 * s7["Ri_critical"] / np.cosh(z) ** 2)["guaranteed_stable"]
    lam_over_h = 2 * PI / ch11.tanh_max_growth()["k"] / 2.0  # tanh(z/L): vorticity thickness h = 2L
    assert abs(lam_over_h - s7["most_unstable_wavelength_over_h"]) / s7["most_unstable_wavelength_over_h"] < 0.015
    c = bv["sec_11_14_chaos"]
    assert ch11.lorenz_hopf_r(c["lorenz_Pr"], c["lorenz_b"]) == pytest.approx(c["derived_by_analyst_lorenz_r_hopf"], rel=1e-6)
    assert rel(ch11.feigenbaum_estimate(8)[-1], c["feigenbaum_delta_printed"]) < 5e-4
    ex = bv["exercises"]
    assert ch11.piecewise_neutral_kh() == pytest.approx(ex["derived_by_analyst_ex_11_11_neutral_kh"], abs=1e-6)
    d = bv["sec_11_10_viscous"]["derived_by_analyst"]
    assert rel(ch11.poiseuille_critical()["Re_c"], d["poiseuille_Re_c"]) < 1e-6


# ======================================================================================================================
# slow: the table writers, the scripts, the reference builder
# ======================================================================================================================
@slow
def test_neutral_curve_tables_and_reference_writer_V3_fast_rebuild(tmp_path):  # V3 (table writers run end to end)
    paths = ch11.neutral_curve_tables(tmp_path, fast=True, cache=False)
    assert set(paths) >= {"os_neutral_poiseuille", "os_grid_blasius", "os_modes"}
    nb = [ln for ln in Path(paths["os_neutral_poiseuille"]).read_text(encoding="utf-8").splitlines() if ln[:1].isdigit()]
    assert len(nb) == 8
    fs = ch11.falkner_skan_neutral_curve(-0.05, Re_values=[300.0, 1000.0], cache=False, N=60)
    fp = ch11.falkner_skan_neutral_curve(0.05, Re_values=[300.0, 1000.0], cache=False, N=60)
    assert np.isfinite(fs["k_upper"]).all() and not np.isfinite(fp["k_lower"][0])  # adverse unstable earlier
    gm = ch11.tg_growth_map(fast=True, cache=False)
    assert np.all(gm["kci"][gm["J"] >= 0.25] == 0)


@slow
def test_scripts_and_make_refs_V7_run_clean(tmp_path):  # V7 smoke (scripts/ch11_*.py --no-show --fast)
    py = sys.executable
    for s in sorted((ROOT / "scripts").glob("ch11_*.py")):
        if s.name in ("ch11_common.py", "ch11_tables.py"):  # ch11_tables rewrites reference/ch11 (run by make_refs)
            continue
        r = subprocess.run([py, str(s), "--no-show", "--fast", "--out", str(tmp_path)], cwd=str(ROOT), capture_output=True,
                           text=True, encoding="utf-8", errors="replace", timeout=900)
        assert r.returncode == 0, (s.name, r.stderr[-2000:])
    before = (REF / "benchmarks.json").read_bytes()
    r = subprocess.run([py, str(REF / "make_refs.py"), "--no-tables"], cwd=str(ROOT), capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=300)
    assert r.returncode == 0 and (REF / "benchmarks.json").read_bytes() == before
