"""Verification suite for Chapter 10 — Computational Fluid Dynamics (Kundu, Cohen & Dowling 5e, §§10.1–10.6).

Evidence levels (``verify-implementation`` skill): V1 analytic · V2 symbolic (sympy / dimensions) · V3 convergence ·
V4 conservation / invariant · V5 published benchmark (``reference/ch10/``, see SOURCES.md) · V6 book value (private,
``tests/book_values_ch10.json``, skipped when absent) · V7 limits / symmetry / invariance. Every test name is
``test_<concept>_<level>_<what>``; the comment on the ``def`` line repeats the level and the curation ID.

A items (CORE, ≥ 2 independent levels): C01 stencils (10.4)–(10.8) · C02 FTCS/BTCS (10.9)–(10.13) · C03 consistency and the
truncation error (10.14)–(10.17) · C04 von Neumann (10.18)–(10.28) · C05 upwind, CFL, Lax (10.29)–(10.31) · C06 weak form
(10.32)–(10.38) · C07 Galerkin (10.39)–(10.63) · C08 element matrices and assembly (10.64)–(10.78) · C09 cell Péclet
(10.84)–(10.94) · C10 MacCormack (10.95)–(10.110) · C11 splitting and projection (10.111)–(10.118) · C12 staggered grid
(10.119)–(10.128) · C13 mixed FE and LBB (10.134)–(10.137), cylinder (10.156)–(10.198) · C14 lid-driven cavity (Ghia 1982,
Hou 1995) · C15 grid convergence, Richardson, GCI.
Derivations: every D row D01–D23 (the ★★ ones and the ★★★ D18 step by step) is re-derived with sympy in a
``test_*_V2_derivation`` test, independently of the chapter's own sympy engines (which are tested separately).

Printed slips kept as wrong variants that must FAIL: R1 ``shape_slopes(printed=True)`` · R2 ``connectivity(printed=True)`` ·
R3 ``steady_cd_fd(scheme="forward")`` · R5 ``cavity_maccormack(printed_step5=True)`` · R6 ``poiseuille_test(printed=True)``
· R11 ``scheme="upwind_printed"`` with u < 0 · R12 ``maccormack_dt_asymptotic`` at Ma = 0.5.

Caches: tests that check *code* run live (``cache=False``); tests that use the heavy runs (128² cavities, fine block grids,
the unsteady cylinder) read our committed result tables in ``reference/ch10/`` and say so — they are evidence about the
cached data, re-checked by the slow tests that re-run the scripts.

Run: ``.venv/Scripts/python.exe -m pytest tests/test_ch10.py -q -p no:cacheprovider``  (``-m "not slow"`` skips the script
smoke test and the long live runs).
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import sympy as sp

from fluidpy import ch10_computational_fluid_dynamics as ch10
from fluidpy.core import fd as FD
from fluidpy.core import fem1d as FEM1
from fluidpy.core import fem2d as FEM2
from fluidpy.core import mac as MAC
from fluidpy.core import maccormack as MCK
from fluidpy.core.diffusion import ftcs_diffusion_1d
from fluidpy.core.units import Q_, dimensional_check
from tools.convergence import grid_convergence_index, observed_order, pairwise_orders, richardson

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference" / "ch10"
BOOK = Path(__file__).resolve().parent / "book_values_ch10.json"
needs_ref = pytest.mark.skipif(not (REF / "ghia1982_table1.csv").exists(), reason="run reference/ch10/make_refs.py")
book_only = pytest.mark.skipif(not BOOK.exists(), reason="book values are private; see CLAUDE.md rule 9")
slow = pytest.mark.slow

ORDER_TOL = 0.15  # design order ± this (verify-implementation default)


def book():
    return json.loads(BOOK.read_text(encoding="utf-8"))


def rel(a, b):
    return abs(float(a) - float(b)) / max(abs(float(b)), 1e-300)


def z0(expr) -> bool:
    return sp.simplify(sp.expand(expr)) == 0


def read_csv(name: str) -> dict:
    """Our result tables in reference/ch10: '#' comment lines, one header line, numbers."""
    lines = [ln for ln in (REF / name).read_text(encoding="utf-8").splitlines() if ln.strip() and not ln.startswith("#")]
    hdr = lines[0].split(",")
    data = np.array([[float(v) for v in ln.split(",")] for ln in lines[1:]])
    return {h: data[:, k] for k, h in enumerate(hdr)}


def header_value(name: str, key: str) -> float:
    txt = "\n".join(ln for ln in (REF / name).read_text(encoding="utf-8").splitlines() if ln.startswith("#"))
    m = re.search(rf"{re.escape(key)}\s*=\s*(-?[0-9.eE+-]+)", txt)
    assert m, f"{key} not in the header of {name}"
    return float(m.group(1))


# ======================================================================================================================
# C01 — stencils from Taylor series (10.4)–(10.8)
# ======================================================================================================================
def test_stencils_V1_exact_rational_weights_of_the_book_stencils():  # V1 (C01, N80, N88): (10.6), (10.7), p. 450, (10.149)
    R = sp.Rational
    assert FD.fd_weights([0, 1], 1) == [-1, 1]
    assert FD.fd_weights([-1, 0], 1) == [-1, 1]
    assert FD.fd_weights([-1, 0, 1], 1) == [R(-1, 2), 0, R(1, 2)]
    assert FD.fd_weights([-1, 0, 1], 2) == [1, -2, 1]
    assert FD.fd_weights([0, 1, 2], 1) == [R(-3, 2), 2, R(-1, 2)]
    assert FD.fd_weights([0, -1, -2, -3], 2) == [2, -5, 4, -1]
    assert np.allclose(FD.fd_weights_float([-1, 0, 1], 2), [1, -2, 1], atol=0)


@pytest.mark.parametrize("offsets,m,order,coef", [((0, 1), 1, 1, sp.Rational(1, 2)), ((-1, 0), 1, 1, sp.Rational(-1, 2)),
                                                   ((-1, 0, 1), 1, 2, sp.Rational(1, 6)), ((-1, 0, 1), 2, 2, sp.Rational(1, 12)),
                                                   ((0, 1, 2), 1, 2, sp.Rational(-1, 3)), ((0, -1, -2, -3), 2, 2, sp.Rational(-11, 12))])
def test_stencils_V2_leading_error_independently_from_exponential_series(offsets, m, order, coef):  # V2 (C01, D01)
    # independent of the code's Taylor table: apply the stencil to e^{a x} at x = 0 and expand in h
    a, h = sp.symbols("a h", positive=True)
    w = FD.fd_weights(offsets, m)
    stencil = sum(wk * sp.exp(a * s * h) for wk, s in zip(w, offsets)) / h ** m
    err = sp.series(stencil - a ** m, h, 0, order + 1).removeO()
    lead = sp.simplify(err.coeff(h, order) / a ** (m + order))
    assert lead == coef
    assert all(sp.simplify(err.coeff(h, k)) == 0 for k in range(order))
    p, c = FD.fd_leading_error(offsets, m)
    assert (p, c) == (order, coef)


def test_stencil_taylor_coefficients_V1_design_rows_and_taylor_table():  # V1 (C01, R02): design C.1 rows 1.2, 1.3
    fw = FD.stencil_taylor_coefficients("forward", 1)
    assert np.allclose(fw["coeffs"], [0, 1, 0.5, 1 / 6, 1 / 24, 1 / 120], atol=1e-15) and fw["order"] == 1 and fw["leading"] == 0.5
    ce = FD.stencil_taylor_coefficients("central", 1)
    assert np.allclose(ce["coeffs"], [0, 1, 0, 1 / 6, 0, 1 / 120], atol=1e-15) and ce["order"] == 2
    c2 = FD.stencil_taylor_coefficients("central2")
    assert np.allclose(c2["coeffs"], [0, 0, 1, 0, 1 / 12, 0], atol=1e-15) and c2["order"] == 2 and c2["leading"] == pytest.approx(1 / 12)
    tab = FD.taylor_table_sympy([1, -1], 4)  # the rows of (10.4) and (10.5)
    R = sp.Rational
    assert list(tab.row(0)) == [1, 1, R(1, 2), R(1, 6), R(1, 24)]
    assert list(tab.row(1)) == [1, -1, R(1, 2), R(-1, 6), R(1, 24)]


def test_stencils_V3_observed_orders_on_sin_and_design_numbers():  # V3 (C01): (10.6)–(10.7), p. 450 on sin x at x0 = 1
    # design Part C 1.5 numbers (h = 0.1)
    assert FD.stencil_error("sin", 1.0, 0.1, "forward") == pytest.approx(-4.294e-2, abs=5e-6)
    assert FD.stencil_error("sin", 1.0, 0.1, "backward") == pytest.approx(4.114e-2, abs=5e-6)
    assert FD.stencil_error("sin", 1.0, 0.1, "central") == pytest.approx(-9.001e-4, abs=5e-8)
    assert FD.stencil_error("sin", 1.0, 0.1, "central2") == pytest.approx(7.010e-4, abs=5e-8)
    assert FD.stencil_error("sin", 1.0, 0.1, "onesided2") == pytest.approx(1.585e-3, abs=5e-7)
    hs = 0.1 / 2.0 ** np.arange(5)
    for kind, p in (("forward", 1), ("backward", 1), ("central", 2), ("central2", 2), ("onesided2", 2)):
        e = [abs(FD.stencil_error("sin", 1.0, h, kind)) for h in hs]
        assert abs(observed_order(hs, e) - p) < ORDER_TOL, (kind, observed_order(hs, e))
    # leading term predicts the error (D01 check): forward ≈ (h/2) f''
    assert FD.stencil_error("sin", 1.0, 1e-3, "forward") == pytest.approx(0.5e-3 * -np.sin(1.0), rel=2e-3)
    # the other test functions exist with exact derivatives
    assert FD.test_function("exp", 0.3, 4) == pytest.approx(np.exp(0.3))
    assert FD.test_function("gauss", 0.2, 2) == pytest.approx((4 * 0.04 - 2) * np.exp(-0.04))


def test_fd_derivative_V1_exact_on_polynomials_and_V3_periodic_order():  # V1+V3 (C01): (10.6)–(10.7) as array operators
    x = np.linspace(0, 1, 11)
    dx = x[1] - x[0]
    q = 3 * x ** 2 - x + 2
    assert np.allclose(FD.fd_derivative(q, dx, "central")[1:-1], (6 * x - 1)[1:-1], atol=1e-12)
    c = x ** 3
    assert np.allclose(FD.fd_derivative(c, dx, "central2", 2)[1:-1], 6 * x[1:-1], atol=1e-11)
    assert np.allclose(FD.fd_derivative(q, dx, "onesided2")[:-2], (6 * x - 1)[:-2], atol=1e-11)
    assert np.allclose(FD.fd_derivative(q, dx, "onesided2_backward")[2:], (6 * x - 1)[2:], atol=1e-11)
    d = FD.fd_derivative(q, dx, "central")
    assert np.isnan(d[0]) and np.isnan(d[-1])
    errs, hs = [], []
    for n in (16, 32, 64, 128):
        xx = np.arange(n) * 2 * np.pi / n
        errs.append(np.max(np.abs(FD.fd_derivative(np.sin(xx), xx[1], "central2", 2, periodic=True) + np.sin(xx))))
        hs.append(xx[1])
    assert abs(observed_order(hs, errs) - 2) < ORDER_TOL


def test_mixed_derivative_V1_exact_on_biquadratics_N88():  # V1 (C15/N88): (10.150) one-sided in x, centred in y
    x = np.linspace(0, 1, 9)
    y = np.linspace(0, 2, 11)
    X, Y = np.meshgrid(x, y, indexing="xy")
    f = 2 * X ** 2 * Y ** 2 + 3 * X * Y - X ** 2 + Y
    exact = 8 * X * Y + 3
    for side in ("backward_x", "forward_x", "backward_y", "forward_y"):
        m = FD.mixed_derivative_onesided(f, x[1] - x[0], y[1] - y[0], side)
        ok = np.isfinite(m)
        assert ok.sum() > 30
        assert np.allclose(m[ok], exact[ok], atol=1e-10), side


def test_time_differences_V3_orders_of_10_8_N07():  # V3 (C01/N07): (10.8) on y' = λy
    assert abs(FD.ode_scheme_order("forward") - 1) < 0.1
    assert abs(FD.ode_scheme_order("backward") - 1) < 0.1
    assert abs(FD.ode_scheme_order("leapfrog") - 2) < 0.1
    assert abs(FD.ode_scheme_order("trapezoidal") - 2) < 0.1


def test_stencils_D01_V2_derivation():  # V2 (C01/D01): every line of Part F D01
    h = sp.symbols("Delta_x", positive=True)
    T = sp.symbols("T0:6")  # T, T_x, T_xx, ...
    Tp = sum(h ** k / sp.factorial(k) * T[k] for k in range(6))  # (10.4)
    Tm = sum((-h) ** k / sp.factorial(k) * T[k] for k in range(6))  # (10.5)
    fw = sp.expand((Tp - T[0]) / h - T[1])  # steps 2–3
    assert fw.coeff(h, 1) == T[2] / 2 and fw.coeff(h, 0) == 0
    bw = sp.expand((T[0] - Tm) / h - T[1])  # step 4
    assert bw.coeff(h, 1) == -T[2] / 2
    dif = sp.expand(Tp - Tm)  # step 5
    assert dif.coeff(h, 1) == 2 * T[1] and dif.coeff(h, 2) == 0 and dif.coeff(h, 3) == T[3] / 3
    ce = sp.expand((Tp - Tm) / (2 * h) - T[1])  # step 6
    assert ce.coeff(h, 2) == T[3] / 6 and ce.coeff(h, 1) == 0 and ce.coeff(h, 3) == 0
    sm = sp.expand(Tp + Tm)  # step 7
    assert sm.coeff(h, 0) == 2 * T[0] and sm.coeff(h, 2) == T[2] and sm.coeff(h, 4) == T[4] / 12 and sm.coeff(h, 5) == 0
    c2 = sp.expand((Tp - 2 * T[0] + Tm) / h ** 2 - T[2])  # step 8
    assert c2.coeff(h, 2) == T[4] / 12 and c2.coeff(h, 0) == 0 and c2.coeff(h, 1) == 0
    # check numbers of D01 (sin at 1, Δx = 0.1): predicted −0.0421 and −9.005e-4
    assert 0.05 * -np.sin(1.0) == pytest.approx(-0.0421, abs=5e-5)
    assert (0.01 / 6) * -np.cos(1.0) == pytest.approx(-9.005e-4, abs=5e-8)


# ======================================================================================================================
# C02 — explicit FTCS and implicit BTCS (10.9)–(10.13)
# ======================================================================================================================
def test_ftcs_V1_design_five_point_example_and_coefficients():  # V1 (C02): (10.10)–(10.11), design C.1 rows 1.9, 1.11
    assert FD.ftcs_coefficients(0.1, 0.01, 0.1, 0.2) == pytest.approx((0.1, 0.2))
    out = FD.transport_1d_step([0, 0, 1, 0, 0], 0.1, 0.2, "ftcs")
    assert np.allclose(out, [0, 0.1, 0.6, 0.3, 0], atol=1e-15)
    assert out.sum() == pytest.approx(1.0, abs=1e-15)
    assert FD.cfl_number(-2.0, 0.1, 0.5) == pytest.approx(0.4)
    assert FD.diffusion_number(0.01, 0.2, 0.1) == pytest.approx(0.2)
    assert FD.cell_peclet(1.0, 0.01, 0.005) == pytest.approx(2.0)  # N19 number
    assert FD.cell_peclet(-1.0, 0.01, 0.005, signed=True) == pytest.approx(-2.0)


def test_ftcs_V1_from_scratch_np_roll_update_agrees():  # V1 (C02): independent one-line (10.10)
    rng = np.random.default_rng(0)
    T = rng.normal(size=37)
    a, b = 0.13, 0.31
    mine = T - a * (np.roll(T, -1) - np.roll(T, 1)) + b * (np.roll(T, -1) - 2 * T + np.roll(T, 1))
    assert np.allclose(FD.transport_1d_step(T, a, b, "ftcs"), mine, atol=1e-14)


def test_ftcs_V1_u_zero_reproduces_the_ch01_diffusion_solver():  # V1 (C02/R05): u = 0 → ch01 FTCS with r = β
    x = np.linspace(0, 1, 21)
    T0 = np.sin(np.pi * x) + 0.3 * x
    dx = x[1] - x[0]
    D, dt = 0.02, 0.4 * dx ** 2 / 0.02
    r = FD.solve_transport_1d(T0, x, 0.0, D, dt, 60, "ftcs")
    F = ftcs_diffusion_1d(T0, D, dx, dt, 60)
    assert np.allclose(r["T"], F[-1], atol=1e-14)


def test_ftcs_V4_periodic_sum_is_conserved():  # V4 (C02): the weights (α+β) + (1−2β) + (β−α) = 1
    rng = np.random.default_rng(1)
    T = rng.random(50)
    s0 = T.sum()
    for _ in range(500):
        T = FD.transport_1d_step(T, 0.05, 0.3, "ftcs")
    assert abs(T.sum() - s0) / s0 < 1e-12


def test_ftcs_V3_second_order_in_dx_with_diffusive_time_step():  # V3 (C02/C03/N10): ‖e‖ ≤ K Δx^a Δt^b, Δt ∝ Δx²
    st = FD.convergence_study("ftcs")
    assert abs(st["order"] - 2) < ORDER_TOL, st
    assert abs(st["pairwise"][-1] - 2) < ORDER_TOL
    cn = FD.convergence_study("cn")
    assert abs(cn["order"] - 2) < ORDER_TOL


def test_btcs_V3_second_order_in_the_asymptotic_range_N08():  # V3 (C02/N08)
    # the default grids 20–160 give a least-squares 1.80 (pairwise 1.53, 1.88, 1.95): pre-asymptotic; refined grids give 2
    d = FD.convergence_study("btcs")
    assert d["pairwise"][0] < d["pairwise"][1] < d["pairwise"][2]  # drifting up towards 2
    assert abs(d["pairwise"][-1] - 2) < ORDER_TOL
    f = FD.convergence_study("btcs", n_list=(80, 160, 320, 640))
    assert abs(f["order"] - 2) < ORDER_TOL and abs(f["pairwise"][-1] - 2) < 0.05


def test_btcs_V1_one_step_solves_the_tridiagonal_system_10_13():  # V1 (C02/N08): residual of (10.13) after one solve
    rng = np.random.default_rng(2)
    Told = rng.random(30)
    a, b = 0.7, 3.0
    Tn = FD.transport_1d_step(Told, a, b, "btcs")
    res = Tn + a * (np.roll(Tn, -1) - np.roll(Tn, 1)) - b * (np.roll(Tn, -1) - 2 * Tn + np.roll(Tn, 1)) - Told
    assert np.max(np.abs(res)) < 1e-13


def test_neumann_ghost_V1_linear_profile_is_steady_N06():  # V1 (C02/N06): ghost T_N = T_{N−2} + 2Δx q is exact for linear T
    x = np.linspace(0, 1, 26)
    q, g = 0.5, 1.0
    T0 = g + q * x
    dx = x[1] - x[0]
    for scheme, dt in (("ftcs", 0.4 * dx ** 2), ("btcs", 50 * dx ** 2), ("cn", 5 * dx ** 2)):
        r = FD.solve_transport_1d(T0, x, 0.0, 1.0, dt, 100, scheme, g=g, q=q)
        assert np.max(np.abs(r["T"] - T0)) < 1e-12, scheme


def test_solve_transport_V7_stability_guard_and_growth_cap():  # V7 (C02/C04): raises outside (10.27); cap stops a blow-up
    x = np.linspace(0, 1, 21)
    T0 = np.sin(np.pi * x)
    dx = x[1] - x[0]
    with pytest.raises(ValueError):
        FD.solve_transport_1d(T0, x, 0.0, 1.0, 0.51 * dx ** 2, 10, "ftcs")
    r = FD.solve_transport_1d(T0 + 1e-8 * (-1.0) ** np.arange(21), x, 0.0, 1.0, 0.6 * dx ** 2, 5000, "ftcs",
                              check_stability=False)
    assert r["blew_up"] and r["step_blown"] < 5000 and np.all(np.isfinite(r["T"]))
    with pytest.raises(ValueError):
        FD.solve_transport_1d(T0, x, 2.0, 0.0, 0.6 * dx, 10, "upwind", periodic=True)  # C = 1.2


def test_ftcs_D02_V2_derivation():  # V2 (C02/D02): (10.9) solved for T^{n+1} is (10.10) with (10.11)
    Tn, Tp, Tm, Tnew, u, D, dt, dx = sp.symbols("T_i T_ip T_im T_new u D Delta_t Delta_x")
    eq = sp.Eq((Tnew - Tn) / dt + u * (Tp - Tm) / (2 * dx), D * (Tp - 2 * Tn + Tm) / dx ** 2)  # (10.9) without the O(·)
    sol = sp.solve(eq, Tnew)[0]
    al, be = u * dt / (2 * dx), D * dt / dx ** 2  # (10.11)
    assert z0(sol - (Tn - al * (Tp - Tm) + be * (Tp - 2 * Tn + Tm)))  # (10.10)
    coeffs = sp.Poly(sp.expand(sol), Tm, Tn, Tp)
    assert z0(sum(coeffs.coeffs()) - 1)  # weights sum to one (D02 check)
    # dimensions (V2): α and β dimensionless
    assert (Q_(1.0, "m/s") * Q_(0.1, "s") / (2 * Q_(0.01, "m"))).to_reduced_units().dimensionless
    assert (Q_(1e-3, "m**2/s") * Q_(0.1, "s") / Q_(0.01, "m") ** 2).to_reduced_units().dimensionless
    dimensional_check(lambda u, dx: 0.5 * u * dx, "kinematic_viscosity", u=Q_(1.0, "m/s"), dx=Q_(0.01, "m"))  # D_num (10.94)


# ======================================================================================================================
# C03 — consistency and the truncation error (10.14)–(10.17)
# ======================================================================================================================
def test_truncation_V2_code_returns_10_17_and_the_other_schemes():  # V2 (C03): truncation_error_sympy
    r = FD.truncation_error_sympy("ftcs")
    s = r["symbols"]
    dt, dx, u, D = s["dt"], s["dx"], s["u"], s["D"]
    E = dt / 2 * s["T_tt"] + u * dx ** 2 / 6 * s["T_xxx"] - D * dx ** 2 / 12 * s["T_xxxx"]
    assert z0(r["E"] - E) and (r["order_t"], r["order_x"]) == (1, 2)
    b = FD.truncation_error_sympy("btcs")
    assert z0(b["E"] - (-dt / 2 * s["T_tt"] + u * dx ** 2 / 6 * s["T_xxx"] - D * dx ** 2 / 12 * s["T_xxxx"]))
    up = FD.truncation_error_sympy("upwind")
    assert (up["order_t"], up["order_x"]) == (1, 1)
    us = FD.truncation_error_sympy("upwind_steady")
    Tx, Txx = us["symbols"]["T_x"], us["symbols"]["T_xx"]
    assert z0(us["modified_equation"].rhs - (D + u * dx / 2) * Txx) and us["modified_equation"].lhs == u * Tx  # (10.94)
    cs = FD.truncation_error_sympy("central_steady")
    assert cs["order_x"] == 2 and z0(cs["modified_equation"].rhs - D * Txx)
    # upwind with the PDE substituted: the Δt u²/2 T_xx and −uΔx/2 T_xx combine into −(uΔx/2)(1 − C) T_xx (D08 note)
    sub = FD.truncation_error_sympy("upwind", substitute_pde=True)
    Dsub = sub["E"].subs(D, 0)
    C = sp.symbols("C", positive=True)
    assert z0(Dsub.subs(dt, C * dx / u) - (-(u * dx / 2) * (1 - C) * Txx))


def test_truncation_D03_V2_derivation_on_an_exponential_mode():  # V2 (C03/D03): independent of the code's Taylor table
    a, b, u, D = sp.symbols("a b u D")
    dt, dx, beta = sp.symbols("Delta_t Delta_x beta", positive=True)
    # T = e^{ax + bt} makes every difference a closed form; divide the residual by T
    time_diff = (sp.exp(b * dt) - 1) / dt  # step 3
    conv = (sp.exp(a * dx) - sp.exp(-a * dx)) / (2 * dx)  # step 4
    diff = (sp.exp(a * dx) - 2 + sp.exp(-a * dx)) / dx ** 2  # step 5
    assert z0(sp.series(time_diff, dt, 0, 2).removeO() - (b + dt / 2 * b ** 2))
    assert z0(sp.series(conv, dx, 0, 3).removeO() - (a + dx ** 2 / 6 * a ** 3))
    assert z0(sp.series(diff, dx, 0, 3).removeO() - (a ** 2 + dx ** 2 / 12 * a ** 4))
    E = (sp.series(time_diff, dt, 0, 2).removeO() - b) + u * (sp.series(conv, dx, 0, 3).removeO() - a) \
        - D * (sp.series(diff, dx, 0, 3).removeO() - a ** 2)  # steps 6–7: E/T
    assert z0(E - (dt / 2 * b ** 2 + u * dx ** 2 / 6 * a ** 3 - D * dx ** 2 / 12 * a ** 4))  # (10.17) with T_tt = b²T …
    # step 11: pure diffusion, b = D a², Δt = βΔx²/D ⇒ E = DΔx²(β/2 − 1/12) T_xxxx
    E11 = E.subs(u, 0).subs(b, D * a ** 2).subs(dt, beta * dx ** 2 / D)
    assert z0(E11 - D * dx ** 2 * (beta / 2 - sp.Rational(1, 12)) * a ** 4)
    # the D03 check number: π⁴·0.01·(0.125 − 1/12) = 0.0406
    assert float(np.pi ** 4 * 0.01 * (0.125 - 1 / 12)) == pytest.approx(0.0406, abs=5e-5)


def test_truncation_terms_V3_measured_residual_minus_10_17_is_fourth_order():  # V3 (C03): O(Δt², Δx⁴) remainder
    x = np.linspace(0, 1, 101)
    rem, tot, hs = [], [], []
    for s in (1, 2, 4, 8):
        dx, dt = 0.01 / s, 0.001 / s ** 2
        tt = FD.truncation_terms(0.5, 0.01, dx, dt, x, 0.05)
        rem.append(np.max(np.abs(tt["measured"] - tt["total"])))
        tot.append(np.max(np.abs(tt["total"])))
        hs.append(dx)
    assert abs(observed_order(hs, rem) - 4) < ORDER_TOL  # remainder O(Δx⁴) with Δt ∝ Δx²
    assert abs(observed_order(hs, tot) - 2) < ORDER_TOL  # (10.17) itself O(Δx²) on this path: consistency
    one = FD.truncation_terms(0.5, 0.01, 0.01, 0.001, 0.37, 0.05)
    assert one["total"] == pytest.approx(one["time"] + one["conv"] + one["diff"])
    assert rel(one["measured"], one["total"]) < 0.05  # "within 5 %" (D03 check)


def test_error_norms_V1_kinds_and_convergence_rates_N10():  # V1 (C03/N10): (10.14)–(10.15)
    a = np.array([1.0, -2.0, 3.0])
    assert FD.error_norm(a, 0 * a) == pytest.approx(np.sqrt(14 / 3))
    assert FD.error_norm(a, 0 * a, "max") == 3.0 and FD.error_norm(a, 0 * a, "l1") == 2.0
    r = FD.convergence_rates(lambda h: 3 * h ** 2, [0.1, 0.05, 0.025])
    assert r["order"] == pytest.approx(2.0) and np.allclose(r["pairwise"], 2.0)
    assert FD.observed_order([0.1, 0.05], [1e-3, 2.5e-4]) == pytest.approx(2.0)


# ======================================================================================================================
# C04 — von Neumann stability (10.18)–(10.28)
# ======================================================================================================================
def test_von_neumann_V1_G_modulus_equals_10_26_and_design_values():  # V1 (C04, N14): |G(10.24)|² = (10.26)
    rng = np.random.default_rng(3)
    th = np.linspace(0, np.pi, 181)
    for a, b in rng.uniform([-0.8, 0.0], [0.8, 0.8], size=(20, 2)):
        G = FD.amplification_factor(th, a, b, "ftcs")
        assert np.allclose(np.abs(G) ** 2, FD.ftcs_amplification_modulus2(th, a, b), atol=1e-12)
        assert np.allclose(FD.amplification_modulus(th, a, b), np.abs(G), atol=1e-15)
    assert FD.amplification_factor(np.pi / 2, 0.1, 0.2) == pytest.approx(0.6 - 0.2j)
    assert FD.ftcs_amplification_modulus2(np.pi / 2, 0.1, 0.2) == pytest.approx(0.40)
    assert FD.amplification_factor(np.pi, 0.0, 100.0, "btcs") == pytest.approx(1 / 401)
    curve = FD.amplification_curve("ftcs", 0.1, 0.2, 91)
    assert curve["theta"][0] == 0 and curve["theta"][-1] == pytest.approx(np.pi) and np.allclose(curve["modulus"], np.abs(curve["G"]))


@pytest.mark.parametrize("scheme,a,b", [("ftcs", 0.1, 0.3), ("ftcs", 0.0, 0.51), ("btcs", 0.4, 2.0), ("cn", 0.3, 1.5),
                                        ("upwind", 0.3, 0.1), ("upwind", -0.3, 0.1), ("lax_wendroff", 0.4, 0.0)])
def test_von_neumann_V1_a_single_fourier_mode_is_multiplied_by_G_each_step(scheme, a, b):  # V1 (C04/D04 step 10)
    N, k, n = 64, 5, 40
    th = 2 * np.pi * k / N
    xi0 = np.cos(th * np.arange(N))
    r = FD.propagate_error(xi0, a, b, n, scheme, growth_cap=1e30)
    rms = np.sqrt(np.mean(r["history"] ** 2, axis=1))
    G = abs(FD.amplification_factor(th, a, b, scheme))
    assert np.allclose(rms / rms[0], G ** np.arange(n + 1), rtol=1e-9)


def test_noye_V1_closed_form_equals_brute_force_scan():  # V1 (C04/N15): (10.27) ⇔ max_θ |G| ≤ 1
    mism = 0
    for a in np.linspace(-0.6, 0.6, 31):
        for b in np.linspace(0.0, 0.7, 29):
            edge = min(abs(4 * a * a - 2 * b), abs(2 * b - 1))
            if edge < 1e-6:
                continue
            mism += FD.ftcs_stable(a, b) != FD.is_von_neumann_stable("ftcs", a, b)
    assert mism == 0
    assert FD.ftcs_stable(0.1, 0.2) and not FD.ftcs_stable(0.4, 0.2) and not FD.ftcs_stable(0.0, 0.51)


def test_noye_V7_blowup_just_outside_decay_just_inside():  # V7 (C04, N11–N13): both edges of (10.27)
    rng = np.random.default_rng(4)
    xi0 = 1e-10 * rng.normal(size=64)
    grow = lambda a, b, n: FD.propagate_error(xi0, a, b, n)["max_abs"][-1] / np.max(np.abs(xi0))  # noqa: E731
    assert grow(0.0, 0.51, 2000) > 1e3 and grow(0.0, 0.49, 2000) < 1.0  # zigzag edge 2β = 1
    assert grow(0.4, 0.2, 3000) > 1e3 and grow(0.1, 0.2, 3000) < 1.0  # long-wave edge 4α² = 2β
    assert FD.max_amplification("ftcs", 0.0, 0.51) == pytest.approx(1.04)
    assert FD.worst_theta("ftcs", 0.0, 0.51) == pytest.approx(np.pi)
    assert FD.worst_theta("ftcs", 0.4, 0.2) < np.pi / 2  # long waves fail first


def test_btcs_V1_unconditionally_stable_N17():  # V1 (C04/N17, D07): |G_BTCS| ≤ 1 for every α, β ≥ 0
    for a in np.linspace(-5, 5, 21):
        for b in (0.0, 0.5, 2.0, 100.0):
            assert FD.max_amplification("btcs", a, b) <= 1 + 1e-14
            assert FD.max_amplification("cn", a, b) <= 1 + 1e-14
    x = np.linspace(0, 1, 41)
    r = FD.solve_transport_1d(np.sin(np.pi * x), x, 0.0, 1.0, 100 * (x[1] - x[0]) ** 2, 50, "btcs")
    assert np.max(np.abs(r["T"])) <= 1.0 and not r["blew_up"]


def test_pure_convection_V1_ftcs_never_stable():  # V1 (C04/R05 new half, D06 step 13)
    th = np.linspace(0.01, np.pi - 0.01, 50)
    for a in (0.01, 0.2, 1.0):
        assert np.allclose(FD.ftcs_amplification_modulus2(th, a, 0.0), 1 + 4 * a * a * np.sin(th) ** 2)
        assert not FD.is_von_neumann_stable("ftcs", a, 0.0)
    r = FD.advect_periodic("gauss", 0.2, 100, 20.0, "ftcs")
    assert r["blew_up"] or r["amplitude_ratio"] > 10


def test_stability_verdict_V1_reason_strings():  # V1 (C04/C05): the explainers' status lines
    assert FD.stability_verdict("ftcs", 0.1, 0.2)["reason"] == "stable"
    assert FD.stability_verdict("ftcs", 0.0, 0.51)["reason"] == "unstable: 2b > 1, the zigzag grows"
    assert FD.stability_verdict("ftcs", 0.4, 0.2)["reason"] == "unstable: 4a^2 > 2b, long waves grow"
    assert FD.stability_verdict("ftcs", 0.3, 0.0)["reason"] == "unstable: pure convection, every wave grows"
    assert FD.stability_verdict("upwind", 0.55, 0.0)["reason"] == "unstable: C > 1, the characteristic leaves the stencil"
    assert FD.stability_verdict("btcs", 3.0, 9.0)["reason"] == "stable for every dt (implicit)"
    assert FD.stability_verdict("upwind_printed", -0.25, 0.0)["stable"] is False


def test_fourier_mode_V1_book_and_standard_conventions_R04():  # V1 (C04/R04): e^{iπkx} ↔ θ = kπΔx
    x = np.linspace(0, 1, 11)
    k = 3.0
    assert np.allclose(FD.fourier_mode(x, k, "book"), FD.fourier_mode(x, np.pi * k, "standard"))
    m = FD.fourier_mode(x, k)
    th = k * np.pi * (x[1] - x[0])
    assert np.allclose(m[1:] / m[:-1], np.exp(1j * th))


def test_von_neumann_D04_V2_derivation():  # V2 (C04/D04)
    al, be, th = sp.symbols("alpha beta theta", real=True)
    xm, x0, xp = sp.symbols("xi_m xi_0 xi_p")
    step3 = x0 - al * (xp - xm) + be * (xp - 2 * x0 + xm)
    step4 = (al + be) * xm + (1 - 2 * be) * x0 + (be - al) * xp  # (10.19)
    assert z0(step3 - step4)
    G = step4.subs({xm: sp.exp(-sp.I * th), x0: 1, xp: sp.exp(sp.I * th)})  # steps 7–9
    G1024 = (al + be) * sp.exp(-sp.I * th) + (1 - 2 * be) + (be - al) * sp.exp(sp.I * th)  # (10.24)
    assert z0(G - G1024)
    assert z0(G.subs(th, 0) - 1)  # a constant error neither grows nor decays
    assert complex(G.subs({al: 0.1, be: 0.2, th: sp.pi / 2})) == pytest.approx(0.6 - 0.2j)
    assert z0(sp.conjugate(G) - G.subs(th, -th))  # G(−θ) = conj G(θ): θ ∈ [0, π] suffices (step 11)


def test_von_neumann_D05_V2_derivation():  # V2 (C04/D05)
    al, be, th = sp.symbols("alpha beta theta", real=True)
    G = ((al + be) * (sp.cos(th) - sp.I * sp.sin(th)) + (1 - 2 * be) + (be - al) * (sp.cos(th) + sp.I * sp.sin(th)))  # step 1
    re, im = sp.re(sp.expand(G)), sp.im(sp.expand(G))
    assert z0(re - (1 - 2 * be * (1 - sp.cos(th))))  # step 2
    assert z0(im - (-2 * al * sp.sin(th)))  # step 3
    assert z0(sp.expand_trig((1 - sp.cos(th)) - 2 * sp.sin(th / 2) ** 2))  # step 4 (half angle)
    mod2 = (1 - 4 * be * sp.sin(th / 2) ** 2) ** 2 + (2 * al * sp.sin(th)) ** 2  # (10.26)
    assert z0(sp.expand_trig(re ** 2 + im ** 2 - mod2))  # steps 5–7


def test_noye_D06_V2_derivation():  # V2 (C04/D06): s-substitution, factorisation, endpoint rule
    al, be, s, th = sp.symbols("alpha beta s theta", real=True)
    hh = sp.symbols("h", real=True)  # θ = 2h
    assert z0(sp.expand_trig(sp.sin(2 * hh) ** 2 - 4 * sp.sin(hh) ** 2 * (1 - sp.sin(hh) ** 2)))  # step 2
    g2 = (1 - 4 * be * s) ** 2 + 16 * al ** 2 * s * (1 - s)  # step 3
    assert z0(sp.expand(g2) - (1 - 8 * be * s + 16 * be ** 2 * s ** 2 + 16 * al ** 2 * s - 16 * al ** 2 * s ** 2))  # step 4
    f = 16 * al ** 2 - 8 * be + s * (16 * be ** 2 - 16 * al ** 2)
    assert z0(g2 - 1 - s * f)  # step 5
    assert z0(f.subs(s, 0) - 4 * (4 * al ** 2 - 2 * be))  # step 9: 4α² ≤ 2β
    assert z0(f.subs(s, 1) - 8 * be * (2 * be - 1))  # step 10: 0 ≤ β ≤ ½
    assert z0(sp.expand_trig(g2.subs(be, 0).subs(s, sp.sin(th / 2) ** 2) - 1 - 4 * al ** 2 * sp.sin(th) ** 2))  # step 13
    # D06 numbers: (0.1, 0.2) stable, (0.4, 0.2) and (0, 0.51) unstable; physics form u²Δt ≤ 2D
    assert 4 * 0.1 ** 2 <= 2 * 0.2 <= 1 and not (4 * 0.4 ** 2 <= 2 * 0.2)


def test_btcs_D07_V2_derivation():  # V2 (C04/D07)
    al, be, th = sp.symbols("alpha beta theta", real=True)
    bracket = 1 + al * (sp.exp(sp.I * th) - sp.exp(-sp.I * th)) - be * (sp.exp(sp.I * th) - 2 + sp.exp(-sp.I * th))  # step 2
    target = 1 + 4 * be * sp.sin(th / 2) ** 2 + 2 * sp.I * al * sp.sin(th)  # steps 3–4
    assert z0(sp.expand_trig(sp.expand(bracket.rewrite(sp.cos)) - target))
    # step 6: |denominator|² − (real part)² = 4α² sin²θ ≥ 0, real part ≥ 1 for β ≥ 0
    assert z0(sp.expand(sp.re(target) ** 2 + sp.im(target) ** 2 - (1 + 4 * be * sp.sin(th / 2) ** 2) ** 2) - 4 * al ** 2 * sp.sin(th) ** 2)
    G = 1 / target  # step 5
    assert complex(G.subs({al: 0, be: 100, th: sp.pi})) == pytest.approx(1 / 401)


# ======================================================================================================================
# C05 — upwind differencing, the CFL condition, Lax equivalence (10.29)–(10.31), (10.199)
# ======================================================================================================================
def test_upwind_V1_design_example_and_exact_shift_at_C_equal_one():  # V1 (C05/N16, D08 step 8)
    assert np.allclose(FD.transport_1d_step([0, 0, 1, 0, 0], 0.25, 0.0, "upwind"), [0, 0, 0.5, 0.5, 0])
    for scheme in ("upwind", "maccormack", "maccormack_bf", "lax_wendroff"):
        for prof in ("square", "gauss"):
            r = FD.advect_periodic(prof, 1.0, 50, 1.0, scheme)
            assert np.max(np.abs(r["T"] - r["T0"])) < 1e-14 and r["steps"] == 50, (scheme, prof)  # one revolution: T = T0
        r = FD.advect_periodic("square", 1.0, 50, 0.3, scheme)  # 15 cells: an exact shift by 15 nodes
        assert np.max(np.abs(r["T"] - np.roll(r["T0"], 15))) < 1e-14


@pytest.mark.parametrize("n_rev", [1.0, 0.3, 2.0])
def test_advect_periodic_V1_reported_error_is_zero_for_an_exact_shift(n_rev):  # V1 (C05, design C.10.1: E3 summary card)
    # At C = 1 upwind is an exact shift (previous test), so the reported rms_error must be 0. (F2, fixed in loop 1: failed for
    # the square pulse: the reference profile f(x − k·C/N) is evaluated in floating point, and (0.2 − 1.0) % 1 =
    # 0.19999999999999996 < 0.2 moves the pulse edge by one node — rms_error = 0.1414 (= √(1/50)) after one revolution.
    r = FD.advect_periodic("square", 1.0, 50, n_rev, "upwind")
    assert r["rms_error"] < 1e-14, f"rms_error = {r['rms_error']:.4g} for an exact shift (n_rev = {n_rev})"


def test_upwind_V1_modulus_formula_of_D08():  # V1 (C05/D08): |G|² = 1 − 2C(1 − C)(1 − cos θ)
    th = np.linspace(0, np.pi, 91)
    for C in (0.2, 0.5, 0.9, 1.0, 1.1):
        G = FD.amplification_factor(th, C / 2, 0.0, "upwind")
        assert np.allclose(np.abs(G) ** 2, 1 - 2 * C * (1 - C) * (1 - np.cos(th)), atol=1e-13)
    assert FD.max_amplification("upwind", 0.55, 0.0) == pytest.approx(1.2)  # C = 1.1: |G(π)| = 2C − 1


def test_cfl_V7_upwind_blows_up_just_above_one_and_is_bounded_below():  # V7 (C05): (10.30)
    hi = FD.advect_periodic("gauss", 1.1, 100, 5.0, "upwind")
    lo = FD.advect_periodic("gauss", 0.9, 100, 5.0, "upwind")
    assert hi["blew_up"]
    assert not lo["blew_up"] and lo["amplitude_ratio"] <= 1.0
    assert FD.is_von_neumann_stable("upwind", 0.45, 0.0) and not FD.is_von_neumann_stable("upwind", 0.55, 0.0)


def test_upwind_V3_first_order_in_the_asymptotic_range():  # V3 (C05): advective Δt = CΔx
    # default grids 20–160: least squares 0.76 (pairwise 0.75, 0.71, 0.83) — pre-asymptotic (Gaussian width 0.05)
    fine = FD.convergence_study("upwind", n_list=(160, 320, 640, 1280), dt_rule="advective")
    assert abs(fine["order"] - 1) < ORDER_TOL
    assert abs(fine["pairwise"][-1] - 1) < 0.05
    assert fine["pairwise"][0] < fine["pairwise"][1] < fine["pairwise"][2]  # approaching 1 from below


def test_upwind_V1_solution_follows_its_modified_equation():  # V1 (C05, D08 note): D_num = |u|Δx(1 − C)/2
    N, C = 200, 0.5
    dx = 1.0 / N
    x = np.arange(N) * dx
    T = FD.advected_gaussian(x, 0.0, 1.0, 0.0)
    steps = int(round(N / C))
    for _ in range(steps):
        T = FD.transport_1d_step(T, C / 2, 0.0, "upwind")
    Dnum = FD.numerical_diffusivity(1.0, dx, scheme="upwind", C=C)
    assert Dnum == pytest.approx(0.5 * dx * (1 - C))
    e_mod = FD.error_norm(T, FD.advected_gaussian(x, 1.0, 1.0, Dnum))
    e_pure = FD.error_norm(T, FD.advected_gaussian(x, 1.0, 1.0, 0.0))
    assert e_mod < 0.1 * e_pure
    assert FD.numerical_diffusivity(1.0, 0.01, C=0.5, scheme="upwind") == pytest.approx(0.0025)  # design example
    assert FD.numerical_diffusivity(1.0, 0.01, C=0.5, scheme="ftcs") < 0  # anti-diffusion


def test_upwind_R11_V7_printed_stencil_fails_for_negative_u_and_mirror_symmetry():  # V7 (C05/N16, slip R11)
    assert FD.max_amplification("upwind_printed", -0.25, 0.0) > 1.4  # C = −0.5 with the i − 1 side: unstable
    assert FD.max_amplification("upwind", -0.25, 0.0) <= 1 + 1e-14
    rng = np.random.default_rng(5)
    T0 = rng.random(40)
    a, b = FD.transport_1d_step(T0, -0.3, 0.1, "upwind"), FD.transport_1d_step(T0[::-1], 0.3, 0.1, "upwind")
    assert np.allclose(a, b[::-1], atol=1e-15)  # u → −u mirrors the solution
    Tp = T0.copy()
    for _ in range(300):
        Tp = FD.transport_1d_step(Tp, -0.3, 0.0, "upwind_printed")
    assert np.max(np.abs(Tp)) > 1e3  # the printed stencil explodes
    assert FD.cfl_number(-1.0, 0.5, 1.0) == 0.5  # |u| in (10.30)


def test_cfl_time_step_V1_climate_numbers():  # V1 (C05 climate note): Δt = CΔx/c with c = √(gH)
    c = np.sqrt(9.81 * 4000.0)
    assert c == pytest.approx(198.1, abs=0.05)
    assert ch10.cfl_time_step(c, 100e3) == pytest.approx(504.8, abs=0.05)
    assert ch10.cfl_time_step(c, 25e3) == pytest.approx(126.2, abs=0.05)
    assert ch10.cfl_time_step(320.0, 25e3) == pytest.approx(78.1, abs=0.05)


def test_lax_equivalence_V7_demo_N18():  # V7 (C05/N18): consistent + stable converges; β = 0.51 blows up
    d = FD.lax_demo()
    assert not d[0.45]["blew_up"] and not d[0.5]["blew_up"]
    assert d[0.45]["final_error"] < 1e-8 and d[0.5]["final_error"] < 1e-8
    assert d[0.51]["blew_up"] and 300 < d[0.51]["step_blown"] < 2000


def test_rod_heating_V1_series_solves_the_heat_equation_N113():  # V1 (C05/N113): (10.199) with our numbers
    D, Tw, L = 0.7, 2.0, 1.5
    x = np.linspace(0, L, 31)
    t = 0.05
    T = ch10.rod_heating_exact(x, t, D, Tw, 0.5, L)
    assert T[0] == pytest.approx(Tw, abs=1e-12) and T[-1] == pytest.approx(Tw, abs=1e-12)
    assert np.allclose(T, T[::-1], atol=1e-12)  # symmetric rod
    h, k = 1e-3, 1e-5
    xi = x[5:-5]
    Tt = (ch10.rod_heating_exact(xi, t + k, D, Tw, 0.5, L) - ch10.rod_heating_exact(xi, t - k, D, Tw, 0.5, L)) / (2 * k)
    Txx = (ch10.rod_heating_exact(xi + h, t, D, Tw, 0.5, L) - 2 * ch10.rod_heating_exact(xi, t, D, Tw, 0.5, L)
           + ch10.rod_heating_exact(xi - h, t, D, Tw, 0.5, L)) / h ** 2
    assert np.max(np.abs(Tt - D * Txx)) < 1e-4 * np.max(np.abs(Tt))
    assert np.allclose(ch10.rod_heating_exact(x, 0.0, D, Tw, 0.5, L)[1:-1], 0.5)


def test_rod_heating_V3_ftcs_converges_to_the_series_at_second_order():  # V3 (C05/N113)
    errs, hs = [], []
    for n in (10, 20, 40, 80):
        x = np.linspace(0, 1, n + 1)
        dx = x[1] - x[0]
        steps = int(round(0.05 / (0.25 * dx ** 2)))
        dt = 0.05 / steps
        T0 = np.zeros(n + 1)
        T0[0] = T0[-1] = 1.0
        r = FD.solve_transport_1d(T0, x, 0.0, 1.0, dt, steps, "ftcs", g=1.0, T_L=1.0)
        errs.append(FD.error_norm(r["T"], ch10.rod_heating_exact(x, 0.05, 1.0, 1.0)))
        hs.append(dx)
    assert abs(observed_order(hs, errs) - 2) < ORDER_TOL


def test_upwind_D08_V2_derivation():  # V2 (C05/D08)
    C, th = sp.symbols("C theta", real=True)
    G = 1 - C * (1 - sp.exp(-sp.I * th))  # step 2
    re, im = sp.re(sp.expand(G.rewrite(sp.cos))), sp.im(sp.expand(G.rewrite(sp.cos)))
    assert z0(re - (1 - C + C * sp.cos(th))) and z0(im + C * sp.sin(th))  # step 3
    g2 = sp.expand(re ** 2 + im ** 2)
    assert z0(sp.simplify(g2 - ((1 - C) ** 2 + 2 * C * (1 - C) * sp.cos(th) + C ** 2)))  # step 5
    assert z0(sp.simplify(g2 - (1 - 2 * C * (1 - C) * (1 - sp.cos(th)))))  # step 6
    assert z0(sp.simplify(G.subs(C, 1) - sp.exp(-sp.I * th)))  # step 8: exact shift
    assert z0(sp.simplify(g2.subs(th, sp.pi) - (2 * C - 1) ** 2))
    assert float(sp.sqrt(g2.subs({th: sp.pi, C: 1.1}))) == pytest.approx(1.2)
    # step 9: printed stencil with C < 0 ⇒ C(1 − C) < 0 ⇒ growth at every θ ≠ 0
    assert float(g2.subs({C: -0.5, th: 1.0})) > 1


def test_phase_error_V1_design_values():  # V1 (C05/C10): relative phase speed −arg G/(Cθ)
    assert FD.phase_error("lax_wendroff", 0.5, np.pi / 2) == pytest.approx(0.7487, abs=5e-5)
    assert FD.phase_error("lax_wendroff", 0.8, np.pi / 4) == pytest.approx(0.9679, abs=5e-5)
    th = np.linspace(0.1, 3.0, 20)
    assert np.allclose(FD.phase_error("upwind", 0.5, th), 1.0)  # C = ½: G = e^{−iθ/2} cos(θ/2)
    assert np.allclose(FD.phase_error("lax_wendroff", 1.0, th), 1.0)


# ======================================================================================================================
# C06 — the weak form (10.32)–(10.38)
# ======================================================================================================================
def test_weak_form_V2_chapter_sympy_engines_vanish():  # V2 (C06, N22): D09, D10 engines
    r = ch10.weak_form_sympy()
    assert r["residual"] == 0 and r["by_parts"] == 0 and r["weak_minus_strong"] == 0
    s = ch10.weak_to_strong_sympy()
    assert s["identity"] == 0 and s["residual"] == 0


def test_weak_form_V1_exact_steady_solution_passes_every_test_function():  # V1 (C06): (10.36) with ∂T/∂t = 0
    u, D, g, q, L = 1.0, 0.5, 0.3, 0.7, 1.0
    B = q * D / (u * np.exp(u * L / D))
    T = lambda x: g + B * (np.exp(u * np.asarray(x) / D) - 1)  # noqa: E731  exact: u T_x = D T_xx, T(0) = g, T_x(L) = q
    for w in (lambda x: np.sin(np.pi * np.asarray(x) / 2), lambda x: np.asarray(x) ** 2, lambda x: np.asarray(x) * (1.3 - np.asarray(x))):
        assert abs(FEM1.weak_residual(T, w, u, D, q, L, g)) < 1e-9
    Twrong = lambda x: g + 1.2 * B * (np.exp(u * np.asarray(x) / D) - 1)  # noqa: E731  wrong flux at L
    assert abs(FEM1.weak_residual(Twrong, lambda x: np.asarray(x), u, D, q, L, g)) > 1e-3


def test_weak_form_D09_V2_derivation():  # V2 (C06/D09): integration by parts with w(0) = 0 and T_x(L) = q
    x, L, D, q = sp.symbols("x L D q", positive=True)
    c = sp.symbols("c0:4")
    b1, b2 = sp.symbols("b1 b2")
    T = sum(c[k] * x ** k for k in range(4))
    w = x * (b1 + b2 * x)
    lhs = sp.integrate(sp.diff(T, x, 2) * w, (x, 0, L))  # step 5 left side
    rhs = (sp.diff(T, x) * w).subs(x, L) - (sp.diff(T, x) * w).subs(x, 0) - sp.integrate(sp.diff(T, x) * sp.diff(w, x), (x, 0, L))
    assert z0(lhs - rhs)
    assert z0((sp.diff(T, x) * w).subs(x, 0))  # step 7: w(0) = 0 removes the Dirichlet end
    # tiny example T = x², w = x on [0, 1]: 1 = 2 − 1
    X = sp.symbols("X")
    assert sp.integrate(2 * X, (X, 0, 1)) == 1 and (2 * X * X).subs(X, 1) - sp.integrate(2 * X, (X, 0, 1)) == 1


def test_weak_form_D10_V2_derivation():  # V2 (C06/D10): weak ⇒ strong and the natural condition
    x, L, D, q, u = sp.symbols("x L D q u", positive=True)
    c = sp.symbols("c0:4")
    T = sum(c[k] * x ** k for k in range(4))
    w = x / L  # step 8: w ∈ V with w(L) = 1
    weak = sp.integrate(u * sp.diff(T, x) * w + D * sp.diff(T, x) * sp.diff(w, x), (x, 0, L)) - D * q * w.subs(x, L)
    strong = sp.integrate((u * sp.diff(T, x) - D * sp.diff(T, x, 2)) * w, (x, 0, L)) + D * (sp.diff(T, x).subs(x, L) - q) * w.subs(x, L)
    assert z0(weak - strong)  # (10.37)
    # if the PDE holds (integral term 0) the remainder is D[T_x(L) − q] (10.38)
    assert z0((strong - sp.integrate((u * sp.diff(T, x) - D * sp.diff(T, x, 2)) * w, (x, 0, L))) - D * (sp.diff(T, x).subs(x, L) - q))


# ======================================================================================================================
# C07 — Galerkin with hat functions (10.39)–(10.63)
# ======================================================================================================================
def test_galerkin_V2_sympy_matrices_and_interior_row():  # V2 (C07, N23–N28, N32)
    r = FEM1.galerkin_equations_sympy(3)
    assert r["check_interior"] == 0
    assert r["check_M"] == sp.zeros(3, 3)
    s = r["symbols"]
    assert r["F"][2] == s["D"] * s["q"]  # the Neumann datum in the last row only (10.57)
    assert not z0(r["K"][0, 1] - r["K"][1, 0])  # K not symmetric with convection (D11 step 12)


def test_galerkin_V1_steady_FE_equals_centred_FD_N33():  # V1 (C07/N33): the mass matrix drops out
    st = FEM1.solve_steady(np.linspace(0, 1, 5), 1.0, 0.25, T_L=1.0)
    assert np.allclose(st["T"], [0, 0.025, 0.1, 0.325, 1.0], atol=1e-14)
    for n, R in ((4, 4.0), (16, 10.0), (20, 60.0)):
        fe = FEM1.solve_steady(np.linspace(0, 1, n + 1), R, 1.0, T_L=1.0)["T"]
        fdv = FD.steady_cd_fd(n, R)["T"]
        assert np.max(np.abs(fe - fdv)) < 1e-12


def test_hat_basis_V1_partition_of_unity_kronecker_and_interpolation_N30_N31():  # V1 (C07/N25–N31)
    nodes = np.array([0.0, 0.2, 0.5, 0.6, 1.0])
    x = np.linspace(0, 1, 101)
    N = FEM1.hat_basis(nodes, x)
    assert np.allclose(N.sum(axis=1), 1.0, atol=1e-14)
    assert np.allclose(FEM1.hat_basis(nodes, nodes), np.eye(5), atol=1e-14)
    d = np.array([1.0, -2.0, 0.5, 3.0])
    assert np.allclose(FEM1.interpolate(d, 0.7, nodes, nodes), np.concatenate([[0.7], d]), atol=1e-14)
    assert FEM1.hat(0.35, nodes, 2) == pytest.approx(0.5)  # halfway up the ramp [0.2, 0.5]


def test_interior_stencil_V1_assembled_rows_are_10_63_times_h():  # V1 (C07/C08, N32): assembly vs closed form
    n, u, D = 8, 1.3, 0.4
    h = 1 / n
    M, K, F = FEM1.assemble_1d(np.linspace(0, 1, n + 1), u, D, sparse=False)
    st = FEM1.interior_stencil(h, u, D)
    for k in range(1, n - 2):  # unknown k ↔ node k+1 (interior)
        assert np.allclose(M[k, k - 1:k + 2], st["M"], atol=1e-15)
        assert np.allclose(K[k, k - 1:k + 2], st["K"], atol=1e-14)
    s = FEM1.interior_stencil(0.25, 1.0, 0.25)
    assert np.allclose(s["stiff"], [-1.5, 2.0, -0.5]) and np.allclose(s["M_over_h"], [1 / 6, 2 / 3, 1 / 6])
    assert abs(s["K"].sum()) < 1e-15  # a constant has no convection or diffusion


def test_fe_steady_V3_dirichlet_and_natural_neumann_converge():  # V3 (C06/C07): T(0) = g essential, T_x(L) = q natural
    u, D, g, q = 1.0, 0.5, 0.3, 0.7
    B = q * D / (u * np.exp(u / D))
    errs, slope_err, hs = [], [], []
    for n in (8, 16, 32, 64):
        x = np.linspace(0, 1, n + 1)
        T = FEM1.solve_steady(x, u, D, g=g, q=q)["T"]
        ex = g + B * (np.exp(u * x / D) - 1)
        assert T[0] == g  # essential: exact on every mesh
        errs.append(np.max(np.abs(T - ex)))
        slope_err.append(abs((T[-1] - T[-2]) / (x[1] - x[0]) - q))
        hs.append(1 / n)
    assert abs(observed_order(hs, errs) - 2) < ORDER_TOL
    assert slope_err[-1] < slope_err[0] and observed_order(hs, slope_err) > 0.85  # the natural slope approached (order ≥ 1)


def test_fe_transient_V3_half_rod_against_the_series_solution():  # V3 (C07/N29, N113): θ-scheme on M ḋ + K d = F
    errs, hs = [], []
    for n in (8, 16, 32):
        x = np.linspace(0, 0.5, n + 1)  # half rod, insulated centre q = 0 (symmetry)
        T0 = np.zeros(n + 1)
        T0[0] = 1.0
        r = FEM1.solve_transport(x, T0, 0.0, 1.0, g=1.0, q=0.0, dt=2e-5, nsteps=2500, theta=0.5)
        errs.append(np.max(np.abs(r["T"] - ch10.rod_heating_exact(x, 0.05, 1.0, 1.0))))
        hs.append(x[1])
    assert abs(observed_order(hs, errs) - 2) < 0.25  # ±0.25: < 2 decades


def test_solve_transport_V1_theta_scheme_agrees_with_method_of_lines_N29():  # V1 (C07/N29)
    x = np.linspace(0, 1, 21)
    T0 = np.exp(-((x - 0.4) / 0.1) ** 2)
    T0[0] = 0.0
    a = FEM1.solve_transport(x, T0, 0.5, 0.02, 0.0, 0.0, dt=1e-3, nsteps=200, theta=0.5)
    b = FEM1.solve_transport(x, T0, 0.5, 0.02, 0.0, 0.0, dt=1e-3, nsteps=200, method="solve_ivp")
    assert np.max(np.abs(a["T"] - b["T"])) < 1e-4
    lumped = FEM1.solve_transport(x, T0, 0.5, 0.02, 0.0, 0.0, dt=1e-3, nsteps=200, lumped=True)
    assert 0 < np.max(np.abs(lumped["T"] - a["T"])) < 0.05


def test_galerkin_D11_V2_derivation():  # V2 (C07/D11): my own element-by-element integrals of three hats
    x, h, u, D, g, q = sp.symbols("x h u D g q", positive=True)
    n = 3
    hat = lambda A, e: ((e * h - x) / h if A == e - 1 else (x - (e - 1) * h) / h if A == e else sp.Integer(0))  # noqa: E731
    I = lambda A, B, da, db: sum(sp.integrate((sp.diff(hat(A, e), x) if da else hat(A, e)) *  # noqa: E731
                                              (sp.diff(hat(B, e), x) if db else hat(B, e)), (x, (e - 1) * h, e * h)) for e in range(1, n + 1))
    M = sp.Matrix(n, n, lambda i, j: I(i + 1, j + 1, False, False))  # (10.55)
    K = sp.Matrix(n, n, lambda i, j: u * I(i + 1, j + 1, False, True) + D * I(i + 1, j + 1, True, True))  # (10.56)
    F = sp.Matrix(n, 1, lambda i, j: D * q * (1 if i + 1 == n else 0) - g * (u * I(i + 1, 0, False, True) + D * I(i + 1, 0, True, True)))
    r = FEM1.galerkin_equations_sympy(3)
    assert sp.simplify(M - r["M"]) == sp.zeros(3, 3)
    assert sp.simplify(K - r["K"]) == sp.zeros(3, 3)
    assert sp.simplify(F - r["F"]) == sp.zeros(3, 1)
    assert sp.simplify(K.subs(u, 0) - K.subs(u, 0).T) == sp.zeros(3, 3)  # pure diffusion: symmetric
    # the numbers of the code match (h = 0.25, u = 1, D = 0.25, g = 0.4, q = 0.3)
    Mn, Kn, Fn = FEM1.assemble_1d(np.linspace(0, 0.75, 4), 1.0, 0.25, g=0.4, q=0.3, sparse=False)
    sub = {h: 0.25, u: 1.0, D: 0.25, g: 0.4, q: 0.3}
    assert np.allclose(np.array(M.subs(sub), float), Mn) and np.allclose(np.array(K.subs(sub), float), Kn)
    assert np.allclose(np.array(F.subs(sub), float).ravel(), Fn)


def test_galerkin_D12_V2_derivation():  # V2 (C07/D12): the hat integrals on a uniform mesh
    s, h, u, D = sp.symbols("s h u D", positive=True)
    rise, fall = s / h, 1 - s / h  # on one element [0, h]
    assert sp.integrate(rise ** 2, (s, 0, h)) * 2 == 2 * h / 3  # step 2
    assert sp.integrate(rise * fall, (s, 0, h)) == h / 6  # step 3
    assert sp.integrate(sp.diff(rise, s) * fall, (s, 0, h)) == sp.Rational(1, 2)  # step 5: ∫N_{A+1,x}N_A on [x_A, x_{A+1}]
    assert sp.integrate(sp.diff(fall, s) * rise, (s, 0, h)) == -sp.Rational(1, 2)  # ∫N_{A−1,x}N_A
    assert sp.integrate(sp.diff(rise, s) * rise, (s, 0, h)) + sp.integrate(sp.diff(fall, s) * fall, (s, 0, h)) == 0
    assert 2 * sp.integrate(sp.diff(rise, s) ** 2, (s, 0, h)) == 2 / h  # step 6
    assert sp.integrate(sp.diff(rise, s) * sp.diff(fall, s), (s, 0, h)) == -1 / h
    row_K = (-u / 2 - D / h, 2 * D / h, u / 2 - D / h)  # step 8
    assert z0(sum(row_K))


# ======================================================================================================================
# C08 — element matrices and assembly (10.64)–(10.78)
# ======================================================================================================================
def test_element_matrices_V1_closed_form_equals_quadrature_and_design_numbers():  # V1 (C08, N37): (10.74)–(10.76)
    for xa, xb, u, D in ((0.0, 0.25, 1.0, 0.25), (0.3, 0.9, -2.0, 0.1), (1.0, 1.01, 5.0, 3.0)):
        m1, k1 = FEM1.element_matrices_linear(xb - xa, u, D)
        m2, k2 = FEM1.element_integrals(xa, xb, u, D)
        assert np.allclose(m1, m2, atol=1e-15) and np.allclose(k1, k2, atol=1e-13)
    m, k = FEM1.element_matrices_linear(0.25, 1.0, 0.25)
    assert np.allclose(m, [[0.25 / 3, 0.25 / 6], [0.25 / 6, 0.25 / 3]]) and np.allclose(k, [[0.5, -0.5], [-1.5, 1.5]])
    assert np.allclose(k.sum(axis=1), 0)  # rows of k sum to zero


def test_shape_slopes_V1_chain_rule_and_R1_printed_labels_fail():  # V1 (C08/N35, D13; slip R1)
    nodes = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    A = 3
    xa, xb = nodes[A - 1], nodes[A]
    xm = np.linspace(xa + 0.01, xb - 0.01, 7)
    h = 1e-6
    num = lambda B: (FEM1.hat(xm + h, nodes, B) - FEM1.hat(xm - h, nodes, B)) / (2 * h)  # noqa: E731
    good = FEM1.shape_slopes(xa, xb, A=A)
    assert good["labels"] == ("A-1", "A") and good["nodes"] == (A - 1, A)
    for node, slope in zip(good["nodes"], good["slopes"]):
        assert np.allclose(num(node), slope, atol=1e-6)
    bad = FEM1.shape_slopes(xa, xb, printed=True, A=A)
    assert not all(np.allclose(num(node), slope, atol=1e-6) for node, slope in zip(bad["nodes"], bad["slopes"]))
    assert FEM1.shape_slopes(0.5, 0.75)["slopes"] == pytest.approx((-4.0, 4.0))


def test_parent_map_V1_round_trip_and_parent_shapes_N34():  # V1 (C08/N34): (10.64)–(10.66)
    xi = np.linspace(-1, 1, 9)
    x = FEM1.map_to_element(xi, 0.3, 0.8)
    assert np.allclose(FEM1.map_to_parent(x, 0.3, 0.8), xi) and x[0] == pytest.approx(0.3) and x[-1] == pytest.approx(0.8)
    N1, N2 = FEM1.parent_shapes(xi)
    assert np.allclose(N1 + N2, 1) and N1[0] == 1 and N2[-1] == 1
    assert np.allclose(N1 * 0.3 + N2 * 0.8, x)  # the map is built from the shapes (isoparametric)


def test_connectivity_V1_and_R2_printed_rows_fail_N38():  # V1 (C08/N38; slip R2)
    c = FEM1.connectivity(4)
    assert c.tolist() == [[0, 1], [1, 2], [2, 3], [3, 4]]
    p = FEM1.connectivity(4, printed=True)
    assert 0 not in p and p.max() == 5  # never touches the Dirichlet node, runs past node 4


def test_element_force_and_assembly_trace_V1_rebuild_assemble_1d():  # V1 (C08/N36–N37): (10.69)–(10.78)
    n, u, D, g, q = 5, 1.2, 0.3, 0.4, 0.9
    h = 1 / n
    f1 = FEM1.element_force(1, n, h, u, D, g, q)
    _, k = FEM1.element_matrices_linear(h, u, D)
    assert np.allclose(f1, -g * k[:, 0])
    assert np.allclose(FEM1.element_force(n, n, h, u, D, g, q), [0, D * q])
    assert np.allclose(FEM1.element_force(3, n, h, u, D, g, q), 0)
    assert np.allclose(FEM1.element_force(1, 1, h, u, D, g, q), -g * k[:, 0] + [0, D * q])
    tr = FEM1.assembly_trace(n, u, D, 1.0, g, q)
    assert [t["nodes"] for t in tr] == FEM1.connectivity(n).tolist()
    M, K, F = FEM1.assemble_1d(np.linspace(0, 1, n + 1), u, D, g, q, sparse=False)
    Mf, Kf = np.array(tr[-1]["M_after"]), np.array(tr[-1]["K_after"])
    assert np.allclose(Mf[1:, 1:], M) and np.allclose(Kf[1:, 1:], K)
    Ff = np.zeros(n + 1)
    for t in tr:
        Ff[t["nodes"]] += t["f"]
    assert np.allclose(Ff[1:], F)


def test_bilinear_form_V1_equals_the_stiffness_entries():  # V1 (C07/N23): a(N_A, N_B) = K_AB (10.42), (10.56)
    n, u, D = 4, 0.8, 0.3
    x = np.linspace(0, 1, n + 1)
    Kf = np.array(FEM1.assembly_trace(n, u, D)[-1]["K_after"])
    E = np.eye(n + 1)
    for A in range(n + 1):
        for B in range(n + 1):
            assert FEM1.bilinear_form(E[A], E[B], x, u, D) == pytest.approx(Kf[A, B], abs=1e-14)


def test_element_D13_V2_derivation():  # V2 (C08/D13)
    xi, x, xa, xb = sp.symbols("xi x x_a x_b", real=True)
    N1, N2 = (1 - xi) / 2, (1 + xi) / 2  # (10.64)
    xmap = N1 * xa + N2 * xb
    assert z0(xmap - ((xb - xa) * xi + xb + xa) / 2)  # (10.65)
    inv = sp.solve(sp.Eq(x, xmap), xi)[0]
    assert z0(inv - (2 * x - xb - xa) / (xb - xa))  # (10.66)
    assert z0(sp.diff(inv, x) - 2 / (xb - xa))  # step 4
    assert z0(sp.diff(N1, xi) * sp.diff(inv, x) + 1 / (xb - xa)) and z0(sp.diff(N2, xi) * sp.diff(inv, x) - 1 / (xb - xa))  # step 6


def test_element_D14_V2_derivation():  # V2 (C08/D14): parent-element integrals and the assembled row
    xi, h, u, D = sp.symbols("xi h u D", positive=True)
    N = [(1 - xi) / 2, (1 + xi) / 2]
    dN = [-1 / h, 1 / h]  # D13
    J = h / 2
    m = sp.Matrix(2, 2, lambda a, b: sp.integrate(N[a] * N[b] * J, (xi, -1, 1)))
    k = sp.Matrix(2, 2, lambda a, b: sp.integrate(u * dN[b] * N[a] * J + D * dN[b] * dN[a] * J, (xi, -1, 1)))
    assert sp.simplify(m - h / 6 * sp.Matrix([[2, 1], [1, 2]])) == sp.zeros(2, 2)
    assert sp.simplify(k - (u / 2 * sp.Matrix([[-1, 1], [-1, 1]]) + D / h * sp.Matrix([[1, -1], [-1, 1]]))) == sp.zeros(2, 2)
    row = (k[1, 0], k[1, 1] + k[0, 0], k[0, 1])  # step 10: node A = local 2 of element A, local 1 of element A + 1
    assert all(z0(a - b) for a, b in zip(row, (-u / 2 - D / h, 2 * D / h, u / 2 - D / h)))  # step 11
    assert z0(m[1, 1] + m[0, 0] - 2 * h / 3) and z0(m[1, 0] - h / 6)
    mn, kn = FEM1.element_matrices_linear(0.3, 1.7, 0.2)
    assert np.allclose(np.array(m.subs(h, 0.3), float), mn) and np.allclose(np.array(k.subs({h: 0.3, u: 1.7, D: 0.2}), float), kn)


# ======================================================================================================================
# C09 — the steady layer, cell Péclet wiggles and upwind numerical diffusion (10.84)–(10.94)
# ======================================================================================================================
def test_steady_cd_exact_V2_solves_10_84_and_V1_overflow_safe_form():  # V2+V1 (C09/N43–N45): (10.86)–(10.88)
    x, R = sp.symbols("x R", positive=True)
    T = (sp.exp(R * x) - 1) / (sp.exp(R) - 1)  # (10.86), L = 1
    assert z0(R * sp.diff(T, x) - sp.diff(T, x, 2)) and T.subs(x, 0) == 0 and sp.simplify(T.subs(x, 1)) == 1
    assert sp.limit(T, R, 0) == x  # conduction limit
    xx = np.linspace(0, 1, 41)
    for Rv in (-30.0, -1.0, 0.5, 4.0, 50.0):
        assert np.allclose(FD.steady_cd_exact(xx, Rv), FD.steady_cd_exact(xx, Rv, printed=True), atol=1e-13, rtol=1e-12)
    big = FD.steady_cd_exact(xx, 1e4)
    assert np.all(np.isfinite(big)) and big[-1] == pytest.approx(1.0)
    assert not np.all(np.isfinite(FD.steady_cd_exact(xx, 1e4, printed=True)))  # the printed form overflows
    assert np.allclose(FD.steady_cd_exact([0.25, 0.5, 0.75], 4.0), [0.03206, 0.1192, 0.3561], atol=5e-5)
    assert np.allclose(FD.steady_cd_exact(xx, 1e-12), xx, atol=1e-9)
    lg = FD.steady_cd_exact(xx, 40.0, limit="large_R")
    assert np.max(np.abs(lg - FD.steady_cd_exact(xx, 40.0))) < 2 * np.exp(-40.0) + 1e-15  # (10.88)


def test_layer_thickness_V1_10_89_and_the_e_folds():  # V1 (C09/N45): δ/L = O(1/R); T(1 − 1/R) → e⁻¹
    assert FD.cd_layer_thickness(100.0) == pytest.approx(0.01, rel=1e-12)
    for R in (50.0, 200.0, 1000.0):
        assert FD.steady_cd_exact(1 - 1 / R, R) == pytest.approx(np.exp(-1), rel=1e-6)
        assert FD.steady_cd_exact(1 - 2 / R, R) == pytest.approx(np.exp(-2), rel=1e-6)
        assert FD.cd_layer_thickness(R) * R == pytest.approx(1.0, rel=1e-6)
    assert FD.cd_layer_thickness(100.0, level=np.exp(-2)) == pytest.approx(0.02, rel=1e-9)


def test_steady_cd_fd_V1_tridiagonal_solve_equals_the_discrete_closed_form():  # V1 (C09/N46, N48, D16)
    assert np.allclose(FD.steady_cd_fd(4, 4.0)["T"], [0, 0.025, 0.1, 0.325, 1], atol=1e-14)
    assert np.allclose(FD.steady_cd_fd(4, 16.0)["T"], [0, -0.05, 0.1, -0.35, 1], atol=1e-14)
    assert np.allclose(FD.steady_cd_fd(4, 16.0, "upwind")["T"], [0, 0.00641, 0.03846, 0.1987, 1], atol=5e-5)
    for n in (4, 10, 25):
        for R in (0.5, 5.0, 2.0 * n, 3.0 * n, 40.0 * n, -7.0):
            for sch in ("central", "upwind"):
                a = FD.steady_cd_fd(n, R, sch)["T"]
                b = FD.steady_cd_discrete_exact(n, R, sch)["T"]
                assert np.max(np.abs(a - b)) < 1e-11, (n, R, sch)
    assert FD.discrete_root(1.0) == pytest.approx(3.0) and FD.discrete_root(4.0) == pytest.approx(-3.0)
    assert FD.discrete_root(4.0, "upwind") == pytest.approx(5.0) and FD.discrete_root(2.0) == float("inf")


def test_cell_peclet_V7_wiggles_if_and_only_if_Rcell_above_two():  # V7 (C09/N19, (10.31), (10.92))
    n = 20
    for Rc in (0.5, 1.0, 1.9, 2.1, 3.0, 10.0, 50.0):
        c = FD.wiggle_indicator(FD.steady_cd_fd(n, Rc * n)["T"])
        up = FD.wiggle_indicator(FD.steady_cd_fd(n, Rc * n, "upwind")["T"])
        assert c["wiggles"] == (Rc > 2), Rc
        assert not up["wiggles"] and not up["has_negative"]
        assert (FD.discrete_root(Rc) < 0) == (Rc > 2)


def test_steady_cd_V3_orders_centred_two_upwind_one():  # V3 (C09)
    ns = (10, 20, 40, 80, 160)
    ec, eu = [], []
    for n in ns:
        x = np.arange(n + 1) / n
        ex = FD.steady_cd_exact(x, 10.0)
        ec.append(np.max(np.abs(FD.steady_cd_fd(n, 10.0)["T"] - ex)))
        eu.append(np.max(np.abs(FD.steady_cd_fd(n, 10.0, "upwind")["T"] - ex)))
    hs = [1 / n for n in ns]
    assert abs(observed_order(hs, ec) - 2) < ORDER_TOL
    assert abs(pairwise_orders(hs, eu)[-1] - 1) < ORDER_TOL


def test_upwind_steady_V1_nodes_approach_the_modified_equation_solution():  # V1 (C09/N49, D17): (10.94)
    # the gap to the exact solution of u T_x = (D + D_num) T_xx shrinks (it is not zero) and faster than the error itself
    gaps, errs, hs = [], [], []
    for n in (10, 20, 40, 80, 160):
        x = np.arange(n + 1) / n
        up = FD.steady_cd_fd(n, 10.0, "upwind")["T"]
        Dnum = FD.numerical_diffusivity(10.0, 1 / n, D=1.0)  # u/D = R: D = 1, u = 10
        assert Dnum == pytest.approx(0.5 * (10.0 / n) * 1.0)  # 0.5 R_cell D
        me = FD.steady_cd_exact(x, 10.0 / (1 + Dnum))
        gaps.append(np.max(np.abs(up - me)))
        errs.append(np.max(np.abs(up - FD.steady_cd_exact(x, 10.0))))
        hs.append(1 / n)
    assert gaps[-1] > 0 and all(g < e for g, e in zip(gaps, errs))
    assert observed_order(hs, gaps) > 1.5 and abs(observed_order(hs, errs) - 1) < 0.25
    assert gaps[-1] / errs[-1] < 0.02


def test_stretched_grid_V7_shrinks_the_wiggles_by_four_orders_N47():  # V7 (C09/N47): the same n, nodes clustered at x = L
    # The stretched grid does NOT remove the sign alternation (the coarse upstream cells now have a local R_cell up to 8),
    # but it shrinks its amplitude from 0.35 to < 1e-4 and the error from 0.39 to 0.05 — the curation's "removes the
    # wiggles" holds for the eye only (documentation correction in the report).
    n, R = 20, 80.0  # R_cell = 4 on the uniform grid
    uni = FD.steady_cd_fd(n, R)["T"]
    assert FD.wiggle_indicator(uni)["wiggles"] and uni.min() < -0.1
    xs = FD.stretched_grid(n, 1.0, 2.0)
    assert np.diff(xs)[-1] < np.diff(xs)[0] and xs[0] == 0 and xs[-1] == pytest.approx(1.0)
    assert np.allclose(FD.stretched_grid(n, 1.0, 0.0), np.linspace(0, 1, n + 1))
    st = FD.steady_cd_fd(n, R, grid="stretched")
    assert st["T"].min() > -1e-4
    assert np.max(np.abs(st["T"] - FD.steady_cd_exact(st["x"], R))) < 0.05 < np.max(np.abs(uni - FD.steady_cd_exact(np.linspace(0, 1, n + 1), R)))


def test_forward_scheme_R3_V1_downwind_root_wiggles_above_Rcell_one():  # V1 (C09/N48, slip R3)
    # R3: the book calls T_j − T_{j−1} "forward"; the true forward (downwind) difference has the root r = 1/(1 − R_cell):
    # monotone for R_cell < 1, oscillating for R_cell > 1 (it is anti-diffusive, D(1 − 0.5R_cell)) — not "every R_cell"
    n = 10
    for Rc, wig in ((0.5, False), (0.9, False), (1.5, True), (3.0, True)):
        assert FD.wiggle_indicator(FD.steady_cd_fd(n, Rc * n, "forward")["T"])["wiggles"] == wig, Rc
    P = sp.symbols("P", positive=True)
    r = sp.symbols("r")
    roots = sp.solve(P * (r ** 2 - r) - (r ** 2 - 2 * r + 1), r)
    assert set(sp.simplify(x) for x in roots) == {1, sp.simplify(1 / (1 - P))}


def test_steady_cd_D15_V2_derivation():  # V2 (C09/D15)
    x, u, D, L, R = sp.symbols("x u D L R", positive=True)
    m = sp.symbols("m")
    Tf = sp.Function("T")
    sol = sp.dsolve(sp.Eq(D * Tf(x).diff(x, 2) - u * Tf(x).diff(x), 0), Tf(x), ics={Tf(0): 0, Tf(L): 1}).rhs
    assert z0(sp.simplify(sol - (sp.exp(u * x / D) - 1) / (sp.exp(u * L / D) - 1)))  # (10.86) with R = uL/D
    assert set(sp.solve(D * m ** 2 - u * m, m)) == {0, u / D}  # step 3
    ratio = sp.simplify(((sp.exp(R * sp.Rational(7, 10)) - 1) / (sp.exp(R) - 1)) / sp.exp(-R * sp.Rational(3, 10)))
    assert sp.limit(ratio, R, sp.oo) == 1  # (10.88) as R → ∞ at fixed ξ = x/L = 0.7 (step 7)
    assert float(sp.exp(-1)) == pytest.approx(0.368, abs=5e-4) and float(sp.exp(-2)) == pytest.approx(0.135, abs=5e-4)


def test_steady_cd_D16_V2_derivation():  # V2 (C09/D16)
    P, r = sp.symbols("P r")
    j, n = sp.symbols("j n", integer=True, positive=True)
    quad = (1 - P / 2) * r ** 2 - 2 * r + (1 + P / 2)  # step 6
    assert z0(sp.expand(sp.Rational(1, 2) * P * (r ** 2 - 1) - (r ** 2 - 2 * r + 1)) + quad)  # step 5 ⇔ step 6
    roots = sp.solve(quad, r)
    assert set(sp.simplify(x) for x in roots) == {1, sp.simplify((1 + P / 2) / (1 - P / 2))}  # step 7
    r2 = (1 + P / 2) / (1 - P / 2)
    Z, rn = sp.symbols("Z r_n")  # Z = r^{j−1}, r_n = r^n: T_{j+k} = (Z r^{k+1} − 1)/(r_n − 1)
    Tk = lambda k: (Z * r2 ** (k + 1) - 1) / (rn - 1)  # noqa: E731  step 8 with k = −1, 0, 1 ↔ j − 1, j, j + 1
    rec = sp.Rational(1, 2) * P * (Tk(1) - Tk(-1)) - (Tk(1) - 2 * Tk(0) + Tk(-1))
    assert z0(sp.simplify(rec))  # the closed form satisfies (10.91) at every interior j
    Tj = lambda k: (r2 ** k - 1) / (r2 ** n - 1)  # noqa: E731
    vals = [float(Tj(k).subs({P: 4, n: 4})) for k in (1, 2, 3)]
    assert np.allclose(vals, [-0.05, 0.1, -0.35])  # tiny example n = 4, R = 16
    up = sp.solve(P * (r - 1) - (r - 1) ** 2, r)  # step 12
    assert set(up) == {1, P + 1}


def test_upwind_steady_D17_V2_derivation():  # V2 (C09/D17): modified equation (10.94)
    a, dx, u, D = sp.symbols("a Delta_x u D", positive=True)
    rhs = sp.exp(a * dx) - 2 + sp.exp(-a * dx)  # (T_{j+1} − 2T_j + T_{j−1}) / T_j for T = e^{ax}, step 3
    # keep R_cell = P fixed (steps 5–7): the balance per Δx² is (P/Δx) a − (1 + P/2) a² + O(Δx)
    Pf = sp.symbols("P", positive=True)
    bal2 = sp.series((Pf * (1 - sp.exp(-a * dx)) - rhs) / dx ** 2, dx, 0, 1).removeO()  # R_cell(T_j − T_{j−1}) − (…)
    assert z0(sp.expand(bal2 - (Pf * a / dx - (1 + Pf / 2) * a ** 2)))
    # with P/Δx = u/D: u T_x = D(1 + 0.5 R_cell) T_xx, i.e. D_num = uΔx/2 (step 8)
    assert z0((Pf * a / dx - (1 + Pf / 2) * a ** 2).subs(Pf, u * dx / D) * D - (u * a - D * (1 + u * dx / (2 * D)) * a ** 2))
    assert FD.numerical_diffusivity(2.0, 0.1) == pytest.approx(0.1)  # D_num = uΔx/2 = 0.5 R_cell D


# ======================================================================================================================
# C10 — weakly compressible NS and MacCormack (10.95)–(10.110)
# ======================================================================================================================
def test_maccormack_V1_equals_lax_wendroff_every_step():  # V1 (C10, D18 steps 1–5)
    assert np.allclose(MCK.maccormack_advection_1d([0, 0, 1, 0, 0], 0.5), [0, -0.125, 0.75, 0.375, 0], atol=1e-15)
    rng = np.random.default_rng(6)
    T0 = rng.normal(size=33)
    for C in (0.3, 0.8, 1.0):
        lw = MCK.lax_wendroff_advection_1d(T0, C, 50)
        assert np.allclose(MCK.maccormack_advection_1d(T0, C, 50, "FB"), lw, atol=1e-11)
        assert np.allclose(MCK.maccormack_advection_1d(T0, C, 50, "BF"), lw, atol=1e-11)
        assert np.allclose(FD.transport_1d_step(T0, C / 2, 0.0, "lax_wendroff"), MCK.lax_wendroff_advection_1d(T0, C, 1), atol=1e-14)


def test_maccormack_step_V1_generic_flux_reduces_to_the_scalar_scheme():  # V1 (C10): (10.101)–(10.102) with E = uT
    rng = np.random.default_rng(7)
    T0 = rng.normal(size=24)
    C = 0.6
    U = MCK.maccormack_step(T0[None, :], lambda V: 2.0 * V, None, dt=0.3, dx=1.0, arrangement="F/B")
    assert np.allclose(U[0], MCK.maccormack_advection_1d(T0, C, 1, "FB"), atol=1e-14)
    # 2-D field depending on y only, advected in y: each column follows the 1-D scheme
    F = np.tile(T0[:, None], (1, 5))[None]
    U2 = MCK.maccormack_step(F, lambda V: 0 * V, lambda V: 0.5 * V, dt=1.2, dx=1.0, dy=1.0, arrangement="FF/BB")
    assert np.allclose(U2[0][:, 3], MCK.maccormack_advection_1d(T0, C, 1, "FB"), atol=1e-14)


def test_maccormack_V2_chapter_engine_lax_wendroff_error_and_G():  # V2 (C10/D18): maccormack_linear_sympy
    r = ch10.maccormack_linear_sympy()
    assert r["lw_difference"] == 0 and r["G2_difference"] == 0
    C = sp.symbols("C", real=True)
    x, h, u = sp.symbols("x h u", positive=True)
    f = sp.Function("f")
    target = h ** 2 * u * (C - 1) * (C + 1) * sp.Derivative(f(x), (x, 3)) / 6  # −(uΔx²/6)(1 − C²) T_xxx per unit time
    assert z0(r["local_error"] - target)


def test_maccormack_D18_V2_derivation_step_by_step():  # V2 (C10/D18 ★★★): every line of Part F D18
    C, th, a = sp.symbols("C theta a", real=True)
    Tm, T0, Tp = sp.symbols("T_m T_0 T_p")
    pred = lambda left, right: left - C * (right - left)  # noqa: E731  step 1
    Ts0, Tsm = pred(T0, Tp), pred(Tm, T0)
    assert z0(Ts0 - Tsm - ((T0 - Tm) - C * (Tp - 2 * T0 + Tm)))  # step 3
    corr = sp.Rational(1, 2) * (T0 + Ts0 - C * (Ts0 - Tsm))  # step 2
    step4 = sp.Rational(1, 2) * (2 * T0 - C * (Tp - T0) - C * (T0 - Tm) + C ** 2 * (Tp - 2 * T0 + Tm))
    assert z0(corr - step4)
    lw = T0 - C / 2 * (Tp - Tm) + C ** 2 / 2 * (Tp - 2 * T0 + Tm)  # step 5
    assert z0(corr - lw)
    xx, dx, u, dt = sp.symbols("x dx u dt", positive=True)
    f = sp.Function("f")
    # steps 6–8: exact step in space derivatives (T_t = −uT_x etc. for T = f(x − ut))
    exact = f(xx - u * dt)
    ser_exact = sp.series(exact, dt, 0, 4).removeO().doit()
    tgt8 = f(xx) - u * dt * f(xx).diff(xx) + u ** 2 * dt ** 2 / 2 * f(xx).diff(xx, 2) - u ** 3 * dt ** 3 / 6 * f(xx).diff(xx, 3)
    assert z0(sp.simplify(ser_exact - tgt8))
    # steps 9–12: the scheme's step and the difference, C held fixed (u dt = C dx)
    scheme = lw.subs({Tm: f(xx - dx), T0: f(xx), Tp: f(xx + dx)})
    err = sp.series(scheme - f(xx - C * dx), dx, 0, 4).removeO().doit()
    assert z0(sp.simplify(err - C * (C - 1) * (C + 1) * dx ** 3 * f(xx).diff(xx, 3) / 6))  # = −(uΔtΔx²/6)(1 − C²)T_xxx
    # steps 13–14: G and |G|²
    G = 1 - sp.I * C * sp.sin(th) - C ** 2 * (1 - sp.cos(th))
    Gw = sp.expand(lw.subs({Tm: sp.exp(-sp.I * th), T0: 1, Tp: sp.exp(sp.I * th)}).rewrite(sp.cos))
    assert z0(sp.expand(Gw - G))
    s = sp.sin(th / 2) ** 2
    G2 = sp.expand(sp.re(G) ** 2 + sp.im(G) ** 2)
    assert z0(sp.expand_trig((G2 - ((1 - 2 * C ** 2 * s) ** 2 + 4 * C ** 2 * s * (1 - s))).subs(th, 2 * a)))
    assert z0(sp.expand_trig((G2 - (1 - 4 * C ** 2 * (1 - C ** 2) * sp.sin(th / 2) ** 4)).subs(th, 2 * a)))
    assert z0(sp.simplify(G.subs(C, 1) - sp.exp(-sp.I * th).rewrite(sp.cos)))  # C = 1: the exact shift


def test_maccormack_V3_second_order_on_a_smooth_pulse():  # V3 (C10): (10.100)–(10.102)
    errs, hs = [], []
    for n in (100, 200, 400, 800):  # N = 50 puts only 2.5 cells in the Gaussian's width (pairwise 1.47 there)
        r = FD.advect_periodic("gauss", 0.5, n, 1.0, "maccormack")
        errs.append(r["rms_error"])
        hs.append(1 / n)
    assert abs(observed_order(hs, errs) - 2) < ORDER_TOL
    lw = FD.convergence_study("lax_wendroff", dt_rule="advective")
    assert abs(lw["order"] - 2) < ORDER_TOL


def test_maccormack_V7_courant_limit_one():  # V7 (C10/D18): stable iff C ≤ 1
    assert FD.advect_periodic("gauss", 1.05, 100, 10.0, "maccormack")["blew_up"]
    ok = FD.advect_periodic("gauss", 0.95, 100, 10.0, "maccormack")
    assert not ok["blew_up"] and ok["amplitude_ratio"] <= 1.0 + 1e-12
    assert FD.max_amplification("maccormack", 0.5, 0.0) <= 1 + 1e-14
    assert FD.max_amplification("maccormack", 0.525, 0.0) > 1


def test_ns_fluxes_V1_isothermal_state_10_99():  # V1 (C10/N53, R09): E = (ρu, ρu² + c²ρ, ρuv), F likewise
    rho, u, v, c = 1.2, 3.0, -2.0, 10.0
    E, F = MCK.ns_fluxes(np.array([rho]), np.array([rho * u]), np.array([rho * v]), c)
    assert np.allclose(E[:, 0], [rho * u, rho * u * u + c * c * rho, rho * u * v])
    assert np.allclose(F[:, 0], [rho * v, rho * u * v, rho * v * v + c * c * rho])
    assert MCK.pressure_isothermal(1.2, 10.0) == pytest.approx(120.0)


def test_ns_coefficients_V1_10_109_and_cavity_coefficients_are_the_same_numbers():  # V1 (C10/N56, N85)
    dt, dx, dy, mu = 1e-3, 0.02, 0.03, 0.4
    c = MCK.ns_coefficients(dt, dx, dy, mu)
    assert c == pytest.approx(dict(c1=dt / dx, c2=dt / dy, c3=mu * dt / dx ** 2, c4=mu * dt / dy ** 2, c5=mu * dt / (12 * dx * dy)))
    Ma, Re = 0.08, 100.0
    a = MCK.cavity_coefficients(dt, dx, dy, Ma, Re)
    b = MCK.coefficients_from_ns(MCK.ns_coefficients(dt, dx, dy, 1 / Re, c=1 / Ma))
    for k in ("a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8", "a9", "a10", "a11"):
        assert a[k] == pytest.approx(b[k], rel=1e-14), k
    assert a["a5"] == pytest.approx(4 * dt / (3 * Re * dx ** 2)) and a["a3"] == pytest.approx(dt / (dx * Ma ** 2))
    assert a["a10"] == pytest.approx(2 * (a["a5"] + a["a6"])) and a["a9"] == pytest.approx(dt / (12 * Re * dx * dy))


def _ns_manufactured(N, Ma, Re):
    X, Y = sp.symbols("x y")
    rho = 1 + sp.Rational(1, 10) * sp.sin(X) * sp.cos(2 * Y)
    u = sp.sin(X) * sp.cos(Y) + sp.Rational(1, 5)
    v = -sp.cos(X) * sp.sin(Y) + sp.Rational(1, 10) * sp.cos(X + Y)
    m, n = rho * u, rho * v
    p = rho / Ma ** 2
    div = sp.diff(u, X) + sp.diff(v, Y)
    rhs = dict(Rr=-(sp.diff(m, X) + sp.diff(n, Y)),  # (10.96)
               Rm=-(sp.diff(m * u + p, X) + sp.diff(m * v, Y)) + (sp.diff(u, X, 2) + sp.diff(u, Y, 2) + sp.diff(div, X) / 3) / Re,  # (10.97)
               Rn=-(sp.diff(m * v, X) + sp.diff(n * v + p, Y)) + (sp.diff(v, X, 2) + sp.diff(v, Y, 2) + sp.diff(div, Y) / 3) / Re)  # (10.98)
    fs = {k: sp.lambdify((X, Y), e, "numpy") for k, e in dict(rho=rho, m=m, n=n, **rhs).items()}
    h = 2 * np.pi / N
    xg = np.arange(N) * h
    XX, YY = np.meshgrid(xg, xg, indexing="xy")
    st = (fs["rho"](XX, YY), fs["m"](XX, YY), fs["n"](XX, YY))
    return st, [fs[k](XX, YY) for k in ("Rr", "Rm", "Rn")], h


def test_ns_predictor_V3_consistent_with_the_compressible_equations_10_96_to_10_98():  # V3 (C10/N54–N56)
    # independent of the code's coefficients: (U* − Uⁿ)/Δt → the right side of (10.96)–(10.98) (μ_v = 0, p = ρ/Ma²);
    # one-sided predictors are first order, the average of FF and BB is second order
    Ma, Re, dt = 0.3, 50.0, 1e-6
    single, avg, hs = [], [], []
    for N in (32, 64, 128):
        st, R, h = _ns_manufactured(N, Ma, Re)
        a = MCK.cavity_coefficients(dt, h, h, Ma, Re)
        sF = MCK.ns_predictor(st, a, "FF/BB", periodic=True)
        sB = MCK.ns_predictor(st, a, "BB/FF", periodic=True)
        single.append(max(np.max(np.abs((sF[k] - st[k]) / dt - R[k])) for k in range(3)))
        avg.append(max(np.max(np.abs(((sF[k] + sB[k]) / 2 - st[k]) / dt - R[k])) for k in range(3)))
        hs.append(h)
    assert abs(observed_order(hs, single) - 1) < ORDER_TOL
    assert abs(observed_order(hs, avg) - 2) < ORDER_TOL


def test_ns_corrector_V1_equals_the_six_substep_periodic_step_and_rho_prime_storage():  # V1+V7 (C10/N55, N85, N90)
    st, _, h = _ns_manufactured(32, 0.3, 50.0)
    a = MCK.cavity_coefficients(2e-3, h, h, 0.3, 50.0)
    full = MCK.ns_corrector(st, MCK.ns_predictor(st, a, "FF/BB"), a, "FF/BB")
    w = MCK.weakly_compressible_step(st, a, "periodic", "FF/BB", perturbation=False)
    assert all(np.allclose(p, q, atol=1e-14) for p, q in zip(full, w))
    prime = MCK.weakly_compressible_step((st[0] - 1.0, st[1], st[2]), a, "periodic", "FF/BB", perturbation=True)
    assert np.allclose(prime[0] + 1.0, w[0], atol=1e-13) and np.allclose(prime[1], w[1], atol=1e-13)
    d = dict(rho=st[0], rhou=st[1], rhov=st[2])
    assert isinstance(MCK.ns_predictor(d, a), dict)  # dict states round-trip


def test_arrangements_V3_ffbb_and_bbff_differ_by_second_order_N57():  # V3 (C10/N57)
    diffs, hs = [], []
    for N in (32, 64, 128):
        st, _, h = _ns_manufactured(N, 0.3, 50.0)
        dt = 0.2 * h * 0.3
        a = MCK.cavity_coefficients(dt, h, h, 0.3, 50.0)
        A = MCK.weakly_compressible_step(st, a, "periodic", "FF/BB", perturbation=False)
        B = MCK.weakly_compressible_step(st, a, "periodic", "BB/FF", perturbation=False)
        diffs.append(max(np.max(np.abs(p - q)) for p, q in zip(A, B)) / dt)
        hs.append(h)
    assert observed_order(hs, diffs) > 2 - ORDER_TOL


def test_compressible_ns_V2_R09_against_the_ch04_equations():  # V2 (C10/R09): (10.96)–(10.98) vs (4.38) with μ_v = 0
    r = ch10.compressible_ns_sympy()
    assert r["difference_x"] == 0 and r["difference_y"] == 0 and r["residual"] == 0


def test_weakly_compressible_V4_mass_conserved_and_V3_taylor_green_convergence_N53():  # V4+V3 (C10/N53)
    res = [MCK.wc_taylor_green(n, Ma=0.05) for n in (16, 32, 64)]
    assert all(abs(r["mass_drift"]) < 1e-14 for r in res)
    e = [r["err_u"] for r in res]
    hs = [2 * np.pi / n for n in (16, 32, 64)]
    assert observed_order(hs, e) > 2 - ORDER_TOL  # measured 2.85 (pre-asymptotic, higher than the design 2)
    # low-Mach behaviour on a fixed grid (reported, not a defect of the transcription): the MacCormack error grows like 1/Ma
    # because Δt ∝ Ma Δx and the acoustic truncation error ∝ c Δx² accumulates — the compressibility error O(Ma²) is hidden
    e_ma = [MCK.wc_taylor_green(32, Ma=m)["err_u"] for m in (0.1, 0.05, 0.025)]
    assert e_ma[0] < e_ma[1] < e_ma[2]
    assert MCK.wc_taylor_green(32, Ma=0.05)["decay_measured"] < MCK.wc_taylor_green(32, Ma=0.05)["decay_exact"]


def test_maccormack_dt_V1_design_value_asymptotic_limit_and_R12():  # V1 (C10/N58, C15/N91, slip R12): (10.110), (10.155)
    assert MCK.maccormack_dt(1, 0, 12.5, 1 / 64, 1 / 64, 1, 0.01) == pytest.approx(2.935e-4, abs=5e-8)
    assert MCK.maccormack_dt_asymptotic(0.08, 1 / 64) == pytest.approx(7.071e-4, abs=5e-8)
    assert MCK.maccormack_dt_asymptotic(0.08, 1 / 32) == pytest.approx(1.414e-3, abs=5e-7)
    dx = 1 / 64
    for Ma in (0.001, 0.01, 0.1, 0.5):
        full = MCK.maccormack_dt(1.0, 1.0, 1 / Ma, dx, dx, 1.0, 1e-14)  # Re_Δ → ∞
        ratio = MCK.maccormack_dt_asymptotic(Ma, dx) / full
        assert ratio == pytest.approx(1 + np.sqrt(2) * Ma, rel=1e-9)  # → 1 only as Ma → 0 (R12)
    assert MCK.maccormack_dt_asymptotic(0.5, dx) / MCK.maccormack_dt(1.0, 1.0, 2.0, dx, dx, 1.0, 1e-14) > 1.7


def test_cavity_density_bc_V1_one_sided_continuity_exact_for_quadratics_N79_N84():  # V1 (C14/N79–N84): (10.138)–(10.146)
    n, dt, U = 12, 1e-3, 1.0
    s = np.linspace(0, 1, n + 1)
    X, Y = np.meshgrid(s, s, indexing="xy")
    h = s[1]
    rho = 1 + 0.2 * X ** 2 - 0.1 * Y + 0.3 * X * Y
    m = 0.5 * X ** 2 - 0.3 * X + 0.2 * Y ** 2  # ρu
    nn = -0.4 * Y ** 2 + 0.7 * Y + 0.1 * X  # ρv
    st = (rho, m, nn)
    mx, ny, rx = X - 0.3, -0.8 * Y + 0.7, 0.4 * X + 0.3 * Y
    left = MCK.cavity_density_bc(st, "left", "predictor", dt, h, h, U)
    assert np.allclose(left, rho[:, 0] - dt * mx[:, 0], atol=1e-13)  # ρ_t = −(ρu)_x (v = 0 along the wall)
    right = MCK.cavity_density_bc(st, "right", "predictor", dt, h, h, U)
    assert np.allclose(right, rho[:, -1] - dt * mx[:, -1], atol=1e-13)
    bot = MCK.cavity_density_bc(st, "bottom", "predictor", dt, h, h, U)
    assert np.allclose(bot, rho[0, 1:-1] - dt * ny[0, 1:-1], atol=1e-13)
    top = MCK.cavity_density_bc(st, "top", "predictor", dt, h, h, U)
    assert np.allclose(top, rho[-1, 1:-1] - dt * (U * rx[-1, 1:-1] + ny[-1, 1:-1]), atol=1e-13)  # (ρu)_x = Uρ_x on the lid
    cor = MCK.cavity_density_bc(st, "left", "corrector", dt, h, h, U, star=st)
    assert np.allclose(cor, 0.5 * (rho[:, 0] + rho[:, 0] - dt * mx[:, 0]), atol=1e-13)


def test_cavity_maccormack_V4_mass_and_R5_printed_step5_breaks_it():  # V4 (C14/N85, slip R5)
    good = MCK.cavity_maccormack(100, n=16, t_end=1.0, cache=False, tol_steady=0)
    assert abs(good["mass_drift"]) < 1e-3 and np.all(np.isfinite(good["u"]))
    with np.errstate(all="ignore"):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            bad = MCK.cavity_maccormack(100, n=16, t_end=1.0, cache=False, tol_steady=0, printed_step5=True)
    assert not np.isfinite(bad["mass_drift"]) or abs(bad["mass_drift"]) > 0.1


def test_block_wall_density_V1_exact_on_manufactured_fields_10_147_to_10_154():  # V1 (C15/N87–N89)
    Ma, Re, h = 0.3, 20.0, 0.1
    k = Ma ** 2 / Re
    N = 9
    i0 = j0 = 4
    idx = np.arange(N)
    J, I = np.meshgrid(idx, idx, indexing="ij")
    for side in ("front", "back", "top", "bottom"):
        if side in ("front", "back"):
            dn = (I - i0) * h  # signed distance along the normal axis (x)
            dt_ = (J - j0) * h
        else:
            dn = (J - j0) * h
            dt_ = (I - i0) * h
        # normal velocity: cubic in dn, zero on the face; tangential velocity: quadratic in dn × linear in the tangent
        un = dn * (0.3 + 0.8 * dn - 0.5 * dn ** 2) * (1 + 0.2 * dt_)
        ut = (0.2 + 0.4 * dn + 0.6 * dn ** 2) * (0.5 + 0.7 * dt_)
        unn = 2 * 0.8 * (1 + 0.2 * dt_)  # ∂²u_n/∂n² at dn = 0
        utnt = (0.4) * 0.7  # ∂²u_t/∂n∂t at dn = 0
        s = k * (4 / 3 * unn + utnt / 3)  # (10.147): ∂ρ/∂n = (Ma²/Re)(4/3 ∂²u_n/∂n² + 1/3 ∂²u_t/∂n∂t)
        rho0 = 0.05 * (1 + 0.3 * dt_)
        rho = rho0 + s * dn + 0.7 * dn ** 2
        u, v = (un, ut) if side in ("front", "back") else (ut, un)
        R = rho + 1.0
        st = (rho, u * R, v * R)
        if side in ("front", "back"):
            out = MCK.block_wall_density(st, side, h, h, Ma, Re, i=np.array([i0]), j=np.array([j0]))
        else:
            out = MCK.block_wall_density(st, side, h, h, Ma, Re, i=np.array([i0]), j=np.array([j0]))
        assert out["rho"][0] == pytest.approx(rho[j0, i0], abs=1e-12), side


def test_body_forces_V1_pressure_drag_sign_and_zero_for_uniform_state():  # V1 (C15/R11): C_D (4.107), C_L (4.108)
    (st0, geo) = MCK.block_init(0.125, 4.0, 3.0, 4.0)
    zero = (np.zeros_like(st0[0]) + 0.01, np.zeros_like(st0[0]), np.zeros_like(st0[0]))
    f = MCK.body_forces(zero, geo, 20.0, 0.1)
    assert abs(f["CD"]) < 1e-12 and abs(f["CL"]) < 1e-12  # uniform pressure: no net force on a closed surface
    Ma = 0.1
    X = np.arange(geo.nx) * geo.dx
    Y = np.arange(geo.ny) * geo.dx
    XX, YY = np.meshgrid(X, Y, indexing="xy")
    lin = (-0.002 * XX, np.zeros_like(XX), np.zeros_like(XX))  # p = ρ′/Ma² falls with x: pushes the block downstream
    f = MCK.body_forces(lin, geo, 20.0, Ma)
    assert f["CD"] == pytest.approx(2 * 0.002 * 1.0 / Ma ** 2, rel=1e-12) and abs(f["CL"]) < 1e-12
    liny = (-0.003 * YY, np.zeros_like(XX), np.zeros_like(XX))
    assert MCK.body_forces(liny, geo, 20.0, Ma)["CL"] == pytest.approx(2 * 0.003 / Ma ** 2, rel=1e-12)
    g = MCK.block_geometry(0.125, 4.0, 8.0, 20.0)
    assert g["nodes_across_block"] == 9 and g["ny"] == 33


def test_block_channel_smoke_V7_short_live_run_is_finite_and_symmetric_at_start():  # V7 (C15/N86): a short live run
    r = MCK.block_channel(20.0, 0.1, 0.25, 2.0, H=4.0, ahead=3.0, behind=6.0, cache=False)
    assert np.all(np.isfinite(r["CD"])) and len(r["t"]) > 5
    assert np.max(np.abs(r["CL"][-5:])) < 0.05 * np.max(np.abs(r["CD"][-5:]))  # symmetric geometry: no lift yet
    assert r["dt"] == pytest.approx(MCK.maccormack_dt_asymptotic(0.1, 0.25, 0.8))


def test_maccormack_D18_supports_the_lax_wendroff_sympy_via_derive_all():  # V2 (C10, C06, C11, C13): derive_all
    d = ch10.derive_all()
    assert all(v == 0 for v in d.values()), d


# ======================================================================================================================
# C11 — operator splitting, Θ-scheme and the projection (10.111)–(10.118), (10.129)–(10.133)
# ======================================================================================================================
def test_splitting_V3_marchuk_yanenko_first_order_N60_N61():  # V3 (C11/N60–N61)
    s = ch10.split_linear_system()
    assert np.max(np.abs(s["commutator"])) > 0  # non-commuting parts
    assert abs(ch10.splitting_order("marchuk_yanenko") - 1) < 0.1
    assert abs(ch10.splitting_order("marchuk_yanenko", dt_list=(0.02, 0.01, 0.005, 0.0025)) - 1) < 0.05
    # V1: with commuting parts the split is exact per substep pair up to the implicit Euler error of each part
    ex = s["exact"](0.0)
    assert np.allclose(ex, s["phi0"])


def test_theta_scheme_V3_second_order_only_at_one_minus_one_over_root_two_N74():  # V3 (C11/N74)
    assert abs(ch10.splitting_order("theta") - 2) < ORDER_TOL
    # θ = ¼: the design expects 1.0 on the default Δt list, measured 0.86 there (pre-asymptotic); finer Δt gives 1
    assert abs(ch10.splitting_order("theta", dt_list=(0.004, 0.002, 0.001, 0.0005), theta_split=0.25) - 1) < 0.05


def test_theta_scheme_V2_the_dt2_coefficient_vanishes_only_at_the_root():  # V2 (C11/N74, ★★★ demoted)
    r = ch10.theta_scheme_amplification_sympy()
    th, l1, l2 = sp.symbols("theta lambda1 lambda2", positive=True)  # the engine's symbols
    assert r["z2_at_root"] == 0 and sp.simplify(r["theta_root"] - (1 - 1 / sp.sqrt(2))) == 0
    c2 = r["z2_coeff"]
    assert th in c2.free_symbols
    assert sp.simplify(c2.subs({th: sp.Rational(1, 4), l1: 1, l2: 0})) != 0  # θ = ¼: the Δt² error survives
    assert sp.simplify(c2.subs({th: 1 - 1 / sp.sqrt(2), l1: 3, l2: 7})) == 0  # the root works for every (λ₁, λ₂)
    assert sp.simplify(r["z3_at_root"] - (106 - 75 * sp.sqrt(2)) / 6) == 0


def test_projection_V4_divergence_removed_to_round_off():  # V4 (C11/C12): (10.117)–(10.124) on a non-square walled grid
    g = MAC.MacGrid(20, 14, 1.0, 0.7)
    rng = np.random.default_rng(8)
    us, vs = rng.normal(size=g.u_shape), rng.normal(size=g.v_shape)
    r = MAC.project(us, vs, g, 0.01, return_parts=True)
    assert np.max(np.abs(r["div_after"])) < 1e-12 * np.max(np.abs(r["div_before"]))
    assert abs(r["rhs_sum"]) < 1e-9
    assert np.max(np.abs(r["curl_correction"])) < 1e-10 * np.max(np.abs(r["corr_u"])) / g.dx  # irrotational correction (N72)
    gp = MAC.MacGrid(16, 16, 2 * np.pi, 2 * np.pi, (True, True))
    up, vp = rng.normal(size=gp.u_shape), rng.normal(size=gp.v_shape)
    u2, v2, p2 = MAC.project(up, vp, gp, 0.1, method="fft")
    u3, v3, p3 = MAC.project(up, vp, gp, 0.1, method="splu")
    assert np.max(np.abs(MAC.divergence(u2, v2, gp))) < 1e-11 and np.allclose(u2, u3, atol=1e-10) and np.allclose(p2, p3, atol=1e-10)


def test_projection_V1_independent_divergence_and_summary_numbers():  # V1 (C12/N66): (10.123) by hand, E6 parity
    g = MAC.MacGrid(8, 8)
    s = MAC.mac_projection_summary(8, "divergent")
    assert s["div_after_max"] < 1e-12 and abs(s["rhs_sum"]) < 1e-12 and s["curl_correction_max"] < 1e-12
    assert s["div_before_max"] > 1 and s["p_max"] > 0 > s["p_min"]
    rng = np.random.default_rng(9)
    u, v = rng.normal(size=g.u_shape), rng.normal(size=g.v_shape)
    d = np.zeros((8, 8))
    for j in range(8):
        for i in range(8):
            d[j, i] = (u[j, i + 1] - u[j, i]) / g.dx + (v[j + 1, i] - v[j, i]) / g.dy
    assert np.allclose(MAC.divergence(u, v, g), d, atol=1e-13)
    st = MAC.projection_stages(u, v, g, 0.5)
    assert st["max_div_after"] < 1e-12 and set(st) >= {"div_before", "rhs", "p", "gx", "gy", "u", "v", "div_after"}


def test_projection_V1_no_pressure_boundary_condition_needed_N70():  # V1 (C12/N70, D20): boundary u* does not matter
    g = MAC.MacGrid(10, 10)
    rng = np.random.default_rng(10)
    us, vs = rng.normal(size=g.u_shape), rng.normal(size=g.v_shape)
    a = MAC.project(us, vs, g, 0.1)
    us2, vs2 = us.copy(), vs.copy()
    us2[:, 0] += 5.0
    vs2[-1, :] -= 3.0
    b = MAC.project(us2, vs2, g, 0.1)
    assert np.allclose(a[0], b[0], atol=1e-12) and np.allclose(a[1], b[1], atol=1e-12)


def test_projection_V2_chapter_engine_D19():  # V2 (C11/D19): projection_sympy
    r = ch10.projection_sympy()
    assert r["divergence_after"] == 0 and r["curl"] == 0


def test_projection_D19_V2_derivation():  # V2 (C11/D19): Poisson equation, curl-free correction, the tiny example
    x, y, dt = sp.symbols("x y Delta_t", positive=True)
    p = sp.Function("p")(x, y)
    us, vs = sp.Function("u_s")(x, y), sp.Function("v_s")(x, y)
    u, v = us - dt * p.diff(x), vs - dt * p.diff(y)  # step 1
    div = sp.expand(u.diff(x) + v.diff(y))
    assert z0(div - (us.diff(x) + vs.diff(y) - dt * (p.diff(x, 2) + p.diff(y, 2))))  # steps 2–3
    assert z0((v.diff(x) - u.diff(y)) - (vs.diff(x) - us.diff(y)))  # step 8: vorticity unchanged
    # tiny example: u* = (x + y, 0), Δt = 1 ⇒ p = x²/2, u = (y, 0)
    P = x ** 2 / 2
    un, vn = (x + y) - P.diff(x), 0 - P.diff(y)
    assert sp.simplify(un - y) == 0 and vn == 0
    assert sp.diff(P, x, 2) + sp.diff(P, y, 2) == sp.diff(x + y, x)  # ∇²p = ∇·u*/Δt
    assert sp.diff(0, x) - sp.diff(x + y, y) == -1 and sp.diff(vn, x) - sp.diff(un, y) == -1


# ======================================================================================================================
# C12 — the staggered grid, the discrete Poisson equation, the checkerboard (10.119)–(10.128)
# ======================================================================================================================
def test_checkerboard_V1_collocated_gradient_zero_staggered_two_over_dx_D21():  # V1 (C12/N68, N69, D21)
    p = ch10.checkerboard(8, 6)
    assert p[0, 0] == 1 and p[0, 1] == -1 and p[1, 0] == -1 and p[1, 1] == 1
    gx, gy = ch10.collocated_gradient(p, 0.5, 0.25)
    assert np.max(np.abs(gx)) == 0 and np.max(np.abs(gy)) == 0
    g = MAC.MacGrid(8, 6, 4.0, 1.5)
    sx, sy = MAC.gradient(p, g)
    assert np.allclose(np.abs(sx[:, 1:-1]), 2 / g.dx) and np.allclose(np.abs(sy[1:-1, :]), 2 / g.dy)
    assert ch10.gradient_null_space(8, 8, "collocated") == 4 and ch10.gradient_null_space(8, 8, "staggered") == 1


def test_checkerboard_D21_V2_derivation():  # V2 (C12/D21)
    i, j = sp.symbols("i j", integer=True)
    p = lambda a, b: (-1) ** (a + b)  # noqa: E731
    dx = sp.symbols("Delta_x", positive=True)
    assert sp.simplify(p(i + 1, j) + p(i, j)) == 0  # step 2: a neighbour has the opposite sign
    assert sp.simplify((p(i + 1, j) - p(i - 1, j)) / (2 * dx)) == 0  # step 3: (10.125) vanishes
    assert sp.simplify((p(i + 1, j) - p(i, j)) / dx + 2 * p(i, j) / dx) == 0  # step 5: (10.126) = ∓2/Δx
    uu = lambda a: (-1) ** a  # noqa: E731  step 7: the same for a zigzag velocity
    assert sp.simplify(uu(i + 1) - uu(i - 1)) == 0 and sp.simplify(uu(i + 1) - uu(i) + 2 * uu(i)) == 0


def test_pressure_poisson_V1_rows_sum_to_zero_rank_and_symmetry():  # V1 (C12/N67): (10.124)
    for nx, ny in ((8, 8), (5, 7), (2, 3)):
        g = MAC.MacGrid(nx, ny, 1.0, 1.3)
        A = MAC.pressure_poisson_matrix(g).toarray()
        assert np.max(np.abs(A @ np.ones(nx * ny))) < 1e-10 * np.max(np.abs(A))
        assert np.allclose(A, A.T)
        assert np.linalg.matrix_rank(A) == nx * ny - 1
        Ap = MAC.pressure_poisson_matrix(g, pin="corner").toarray()
        assert np.linalg.matrix_rank(Ap) == nx * ny
    A8 = MAC.pressure_poisson_matrix(MAC.MacGrid(8, 8)).toarray()
    assert np.linalg.matrix_rank(A8) == 63  # D20 check


def test_pressure_poisson_V1_single_cell_direction_rows_sum_to_zero_D20_pipe():  # V1 (C12/D20 tiny example)
    # The 3-cell pipe of D20: walls above and below one row of cells. (10.124) must reduce to the rows (−1, 1, 0),
    # (1, −2, 1), (0, 1, −1) (× 1/Δx²) with a constant null vector. (F1, fixed in loop 1: _lap1d with n = 1 used to set the
    # single wall-bounded cell to −1/Δy² (a Dirichlet ghost) instead of 0 — the y-part of every row is −1/Δy².
    g = MAC.MacGrid(3, 1, 3.0, 1.0)
    A = MAC.pressure_poisson_matrix(g).toarray()
    assert np.allclose(A, [[-1, 1, 0], [1, -2, 1], [0, 1, -1]]), A
    p = MAC.solve_pressure(np.array([[1.0, 0.0, -1.0]]), g, pin="corner")
    assert np.allclose(p, [[0.0, 1.0, 2.0]])


@pytest.mark.parametrize("nx,ny,per", [(4, 1, (False, False)), (1, 4, (False, False)), (4, 1, (True, False)),
                                       (1, 4, (True, False)), (2, 4, (True, False)), (4, 2, (False, True))])
def test_pressure_poisson_V1_constant_null_vector_on_thin_grids(nx, ny, per):  # V1 (C12/N67, D20 step 10)
    # A·1 = 0 must hold on every grid (a constant pressure has no gradient). FAILS for a direction with one cell (wall or
    # periodic: the diagonal gets −1/h² instead of 0) and for a periodic direction with two cells (the two neighbours are
    # the same cell: the off-diagonal must be 2/h², it is 1/h²).
    A = MAC.pressure_poisson_matrix(MAC.MacGrid(nx, ny, 1.0, 1.0, per)).toarray()
    rowsum = np.max(np.abs(A @ np.ones(nx * ny)))
    assert rowsum < 1e-9 * np.max(np.abs(A)), f"max |A·1| = {rowsum:.3g} (|A|max = {np.max(np.abs(A)):.3g})"


def test_pressure_poisson_D20_V2_derivation():  # V2 (C12/D20): the cell budget → (10.124), boundary rows, compatibility
    dt, dx, dy = sp.symbols("Delta_t Delta_x Delta_y", positive=True)
    ue, uw, vn, vs_ = sp.symbols("u_e u_w v_n v_s")  # predicted face velocities
    pc, pe, pw, pn, ps = sp.symbols("p_c p_e p_w p_n p_s")
    ue1, uw1 = ue - dt / dx * (pe - pc), uw - dt / dx * (pc - pw)  # step 2 (10.121)
    vn1, vs1 = vn - dt / dy * (pn - pc), vs_ - dt / dy * (pc - ps)  # (10.122)
    budget = (ue1 - uw1) / dx + (vn1 - vs1) / dy  # (10.123)
    lap = (pe - 2 * pc + pw) / dx ** 2 + (pn - 2 * pc + ps) / dy ** 2
    assert z0(budget - (((ue - uw) / dx + (vn - vs_) / dy) - dt * lap))  # steps 3–6: ∇²p = ∇·u*/Δt
    uwall = sp.symbols("u_wall")
    wall_budget = (ue1 - uwall) / dx + (vn1 - vs1) / dy  # steps 7–9: the wall face is not corrected
    assert not wall_budget.has(pw)
    # step 10–12 on the 3-cell pipe (numbers of the D20 check): rows sum to zero, Σ rhs = 0, p = (0, 1, 2)
    A = np.array([[-1.0, 1, 0], [1, -2, 1], [0, 1, -1]])
    rhs = np.array([1.0, 0, -1])
    assert np.allclose(A.sum(axis=1), 0) and rhs.sum() == 0
    Ap = A.copy()
    Ap[0] = [1, 0, 0]
    p = np.linalg.solve(Ap, np.array([0.0, 0, -1]))
    assert np.allclose(p, [0, 1, 2]) and np.allclose(A @ p, rhs)


def test_solve_pressure_V1_methods_agree_and_compatibility_is_enforced():  # V1 (C12/N67, N41)
    g = MAC.MacGrid(12, 10)
    rng = np.random.default_rng(11)
    rhs = rng.normal(size=(10, 12))
    rhs -= rhs.mean()
    a = MAC.solve_pressure(rhs, g, "splu")
    b = MAC.solve_pressure(rhs, g, "sor")
    assert np.max(np.abs(a - b)) < 1e-8 and abs(a.mean()) < 1e-12
    A = MAC.pressure_poisson_matrix(g)
    assert np.max(np.abs((A @ a.ravel()).reshape(10, 12) - rhs)) < 1e-9
    with pytest.raises(ValueError):
        MAC.solve_pressure(rhs + 1.0, g)
    c = MAC.solve_pressure(rhs, g, pin="corner")
    assert c[0, 0] == 0 and np.allclose(c - c.mean(), a)  # the constant is a gauge (N41)
    assert np.allclose(MAC.pin_pressure(a + 7.0), a)


def test_mac_V1_channel_poiseuille_exact_with_quadratic_ghosts_and_order_two_linear():  # V1+V3 (C12/C11)
    r = MAC.channel_poiseuille()
    assert r["max_err"] < 1e-11
    errs, hs = [], []
    for ny in (8, 16, 32):
        errs.append(MAC.channel_poiseuille(ny=ny, ghost="linear")["max_err"])
        hs.append(1 / ny)
    assert abs(observed_order(hs, errs) - 2) < ORDER_TOL


def test_mac_V3_taylor_green_decays_at_the_exact_rate():  # V1+V3 (C11/C12): periodic MAC vs the exact vortex
    r = MAC.taylor_green(32, 100.0, 1.0)
    assert r["err_u"] < 1e-4 and r["err_p"] < 5e-3
    assert abs(r["order"] - 2) < ORDER_TOL


def test_mac_stability_V1_von_neumann_scan_of_10_127_and_10_128_N73():  # V1 (C12/N73)
    for Re, dx in ((100.0, 1 / 32), (10.0, 1 / 32), (1000.0, 1 / 16)):
        lim = MAC.dt_limit(1.0, 0.0, Re, dx, safety=1.0)
        assert FD.ftcs2d_max_amplification(1.0, 0.0, Re, dx, 0.95 * lim) <= 1 + 1e-12
        assert FD.ftcs2d_max_amplification(1.0, 0.0, Re, dx, 1.05 * lim) > 1
    rng = np.random.default_rng(12)
    for ang in rng.uniform(0, 2 * np.pi, 6):
        u, v = np.cos(ang), np.sin(ang)
        lim = MAC.dt_limit(abs(u), abs(v), 100.0, 1 / 32, safety=1.0)
        assert FD.ftcs2d_max_amplification(u, v, 100.0, 1 / 32, 0.95 * lim) <= 1 + 1e-12  # sufficient in every direction
    assert MAC.dt_limit(1, 0, 100, 1 / 32, safety=1) == pytest.approx(0.02)
    assert MAC.dt_limit(1, 0, 100, 1 / 64, safety=1) == pytest.approx(0.0061035, abs=5e-8)
    g = MAC.MacGrid(8, 8)
    with pytest.raises(ValueError):
        MAC.run(MAC.new_state(g), g, 100.0, 0.05, 2, lid=1.0)  # convective limit 0.02 < 0.05 (lid default is now the grid's)


def test_mac_operators_V3_laplacian_and_convection_forms_agree_on_a_divergence_free_field():  # V3 (C12/R08, N64)
    e_lap, e_conv, hs = [], [], []
    for n in (16, 32, 64):
        g = MAC.MacGrid(n, n, 2 * np.pi, 2 * np.pi, (True, True))
        c = MAC.face_coordinates(g)
        u = np.sin(c["xu"]) * np.cos(c["yu"])
        v = -np.cos(c["xv"]) * np.sin(c["yv"])
        lu, lv = MAC.laplacian_faces(u, v, g)
        e_lap.append(np.max(np.abs(lu + 2 * u)))
        ca, _ = MAC.convective_terms(u, v, g, "advective")
        cc, _ = MAC.convective_terms(u, v, g, "conservative")
        exact = 0.5 * np.sin(2 * c["xu"])  # (u·∇)u for Taylor–Green, x component
        e_conv.append(max(np.max(np.abs(ca - exact)), np.max(np.abs(cc - exact))))
        hs.append(g.dx)
    assert abs(observed_order(hs, e_lap) - 2) < ORDER_TOL and abs(observed_order(hs, e_conv) - 2) < ORDER_TOL
    w = MAC.vorticity(u, v, g)
    assert np.allclose(w, MAC.curl(u, v, g))


# ======================================================================================================================
# C13 — mixed finite elements, LBB, the FE cylinder (10.134)–(10.137), (10.156)–(10.198)
# ======================================================================================================================
def test_p2_shapes_V1_kronecker_partition_and_gradients_N105():  # V1 (C13/N105): (10.185), (10.187)
    nodes = np.array([(0, 0), (1, 0), (0, 1), (0.5, 0), (0.5, 0.5), (0, 0.5)])
    Phi = FEM2.p2_shape(nodes[:, 0], nodes[:, 1])
    assert np.allclose(Phi, np.eye(6), atol=1e-15)
    xi, eta = np.random.default_rng(13).random((2, 20)) * 0.5
    assert np.allclose(FEM2.p2_shape(xi, eta).sum(axis=0), 1.0)
    h = 1e-6
    num_x = (FEM2.p2_shape(xi + h, eta) - FEM2.p2_shape(xi - h, eta)) / (2 * h)
    num_y = (FEM2.p2_shape(xi, eta + h) - FEM2.p2_shape(xi, eta - h)) / (2 * h)
    G = FEM2.p2_shape_grad(xi, eta)
    assert np.allclose(G[:, 0], num_x, atol=1e-8) and np.allclose(G[:, 1], num_y, atol=1e-8)
    P1 = FEM2.p1_shape(nodes[:3, 0], nodes[:3, 1])
    assert np.allclose(P1, np.eye(3))


def test_quadrature_V1_seven_point_rule_exact_to_degree_five_not_six_N108():  # V1 (C13/N108): (10.198)
    P, W = FEM2.tri_quad_7pt()
    assert W.sum() == pytest.approx(1.0, abs=1e-15)
    from math import factorial
    exact = lambda p, q: factorial(p) * factorial(q) / factorial(p + q + 2)  # noqa: E731  ∫ ξ^p η^q over the parent triangle
    for p in range(6):
        for q in range(6 - p):
            assert 0.5 * np.sum(W * P[:, 0] ** p * P[:, 1] ** q) == pytest.approx(exact(p, q), abs=1e-15), (p, q)
    worst = max(abs(0.5 * np.sum(W * P[:, 0] ** p * P[:, 1] ** (6 - p)) - exact(p, 6 - p)) for p in range(7))
    assert worst > 1e-6


def test_iso_map_V1_affine_jacobian_and_curved_integration_N104():  # V1 (C13/N104, N108)
    xe = np.array([1.0, 3.0, 1.0, 2.0, 2.0, 1.0])
    ye = np.array([0.0, 0.0, 2.0, 0.0, 1.0, 1.0])  # straight triangle (1,0), (3,0), (1,2), mid-nodes on the sides
    assert FEM2.jacobian(xe, ye, 0.2, 0.3) == pytest.approx(2 * 2.0)  # J = 2 × area
    x, y = FEM2.iso_map(xe, ye, np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0]))
    assert np.allclose(x, [1, 3, 1]) and np.allclose(y, [0, 0, 2])
    assert FEM2.integrate_element(lambda X, Y: X * Y + Y ** 2, xe, ye) == pytest.approx(10 / 3, rel=1e-13)  # by hand
    # a quarter of the unit disc's sector: vertices (0,0), (1,0), (0,1), mid-node of the arc on the circle → area ≈ π/4
    s = np.sqrt(0.5)
    xc = np.array([0, 1, 0, 0.5, s, 0])
    yc = np.array([0, 0, 1, 0, s, 0.5])
    area = FEM2.integrate_element(lambda X, Y: np.ones_like(X), xc, yc)
    assert abs(area - np.pi / 4) < 0.01 and abs(area - 0.5) > 0.2


def test_meshes_V1_counts_and_euler_N94():  # V1 (C13/N94): N_u = V + E, N_p = V, V − E + T = 1 − holes
    m = FEM2.structured_square_mesh(4)
    c = m["counts"]
    assert (c["V"], c["E"], c["T"], c["N_u"], c["N_p"], c["euler"]) == (25, 56, 32, 81, 25, 1)
    assert FEM2.p2_node_counts(25, 56, 32) == dict(N_u=81, N_p=25, euler=1)
    cyl = FEM2.cylinder_channel_mesh()
    assert FEM2.euler_check(cyl)["euler"] == 0  # one hole
    on = np.abs(np.hypot(*cyl.nodes[cyl.boundary["cylinder"]].T) - 0.5) < 1e-9
    assert on.all() and cyl.curved
    y = cyl.nodes[:, 1]
    assert np.allclose(np.sort(y), np.sort(-y), atol=1e-9)  # up–down symmetric mesh
    p1 = FEM2.structured_square_mesh(3, "P1P1")
    assert p1["counts"]["N_u"] == 16


def test_mixed_fe_V1_poiseuille_exact_R6_printed_fails_N109():  # V1 (C13/N109, slip R6)
    a = FEM2.poiseuille_test()
    assert a["err_u"] < 1e-11 and a["err_v"] < 1e-11 and a["err_p"] < 1e-10
    b = FEM2.poiseuille_test(method="penalty")
    assert b["err_u"] < 1e-8 and b["err_p"] < 1e-7  # both essential-BC methods give the same flow
    with np.errstate(all="ignore"):
        bad = FEM2.poiseuille_test(printed=True)
    assert not np.isfinite(bad["err_u"]) or bad["err_u"] > 1e-2


def test_mixed_fe_V3_kovasznay_orders_three_and_two():  # V3 (C13/N95–N108): Newton + P2–P1 on Kovasznay's exact flow
    r16 = FEM2.kovasznay_test(16)
    r32 = FEM2.kovasznay_test(32)
    assert abs(r32["orders"]["u"] - 3) < 0.25 and r32["orders"]["p"] > 2 - ORDER_TOL
    assert r32["err_u"] < r16["err_u"] < 0.01
    res = r16["residuals"]
    assert res[-1] < 1e-9 and all(res[k + 1] < 5 * res[k] ** 2 for k in range(1, len(res) - 1))  # quadratic Newton (N98)


def test_kovasznay_V2_exact_field_solves_steady_navier_stokes():  # V2 (C13): the exact field used above
    x, y, Re = sp.symbols("x y Re", positive=True)
    lam = Re / 2 - sp.sqrt(Re ** 2 / 4 + 4 * sp.pi ** 2)
    u = 1 - sp.exp(lam * x) * sp.cos(2 * sp.pi * y)
    v = lam / (2 * sp.pi) * sp.exp(lam * x) * sp.sin(2 * sp.pi * y)
    p = (1 - sp.exp(2 * lam * x)) / 2
    lap = lambda f: sp.diff(f, x, 2) + sp.diff(f, y, 2)  # noqa: E731
    rx = u * u.diff(x) + v * u.diff(y) + p.diff(x) - lap(u) / Re
    ry = u * v.diff(x) + v * v.diff(y) + p.diff(y) - lap(v) / Re
    for R in (10, 40):
        for pt in ((0.1, 0.2), (0.7, -0.3)):
            sub = {Re: R, x: pt[0], y: pt[1]}
            assert abs(float(rx.subs(sub))) < 1e-12 and abs(float(ry.subs(sub))) < 1e-12
    assert z0(u.diff(x) + v.diff(y))
    un, vn, pn = FEM2.kovasznay_exact(0.1, 0.2, 40.0)
    assert un == pytest.approx(float(u.subs({Re: 40, x: 0.1, y: 0.2})))


def test_lbb_V1_taylor_hood_bounded_equal_order_spurious_C13():  # V1+V3 (C13): inf–sup constants
    betas = [FEM2.infsup_constant(n, "P2P1")["beta"] for n in (2, 4, 8)]
    assert min(betas) > 0.36 and max(betas) - min(betas) < 0.01  # mesh independent
    p1 = [FEM2.infsup_constant(n, "P1P1") for n in (4, 8)]
    assert all(r["n_spurious"] > 0 and r["beta"] == 0 for r in p1)
    assert p1[1]["beta_nonspurious"] < p1[0]["beta_nonspurious"]  # falling with refinement
    tab = read_csv("infsup_table.csv")
    p2rows = tab["pair"] == 1
    assert np.all(tab["beta"][p2rows] > 0.36) and np.all(tab["n_spurious"][~p2rows] > 0)
    t = FEM2.infsup_table(("P2P1",), (2, 3))
    assert t["P2P1"][0]["beta"] == pytest.approx(0.36657, abs=5e-5)


def test_stokes_cavity_V1_equal_order_pressure_is_rank_deficient():  # V1 (C13): E8's two pictures
    a = FEM2.stokes_cavity(4, "P2P1")
    b = FEM2.stokes_cavity(4, "P1P1")
    assert a["rank_deficiency"] == 0 and a["residual"] < 1e-10
    assert b["rank_deficiency"] > 0
    assert b["p_std_checker"] > 2 * a["p_std_checker"]


def test_weak_ns_V2_chapter_engines_D22_N95():  # V2 (C13/D22, N95)
    assert ch10.weak_ns_identity_sympy()["residual"] == 0
    r = ch10.weak_ns_components_sympy()
    assert r["difference"] == 0 and r["book_bracket"] == 0


def test_weak_ns_D22_V2_derivation():  # V2 (C13/D22): steps 2, 4, 5, 8 independently
    x, y = sp.symbols("x y")
    p, ut, vt, u, v = (sp.Function(n)(x, y) for n in ("p", "ut", "vt", "u", "v"))
    assert z0((p.diff(x) * ut + p.diff(y) * vt) - ((p * ut).diff(x) + (p * vt).diff(y) - p * (ut.diff(x) + vt.diff(y))))  # step 2
    G = sp.Matrix([[u.diff(x), u.diff(y)], [v.diff(x), v.diff(y)]])  # G_ij = ∂u_i/∂x_j
    divT = [sum(sp.diff((G + G.T)[i, jj], (x, y)[jj]) for jj in range(2)) for i in range(2)]
    lap = [u.diff(x, 2) + u.diff(y, 2), v.diff(x, 2) + v.diff(y, 2)]
    dv = u.diff(x) + v.diff(y)
    assert z0(divT[0] - (lap[0] + dv.diff(x))) and z0(divT[1] - (lap[1] + dv.diff(y)))  # step 4
    a, b, c, d = sp.symbols("a b c d")
    S = sp.Matrix([[a, b], [b, c]])
    W = sp.Matrix([[0, d], [-d, 0]])
    assert sum(S[i, jj] * W[i, jj] for i in range(2) for jj in range(2)) == 0  # step 8: symmetric : antisymmetric = 0


def test_newton_V3_quadratic_convergence_on_the_cylinder_and_symmetry_N98_R12():  # V3+V7 (C13/N98, R12)
    r = FEM2.cylinder_steady(1.0, "coarse", cache=False)
    res = r["residuals"]
    assert r["converged"] and res[-1] < 1e-10
    assert res[3] < 10 * res[2] ** 2  # quadratic
    assert abs(r["CL"]) < 1e-10 and abs(r["CM"]) < 1e-10  # symmetric mesh ⇒ no lift, no torque
    cyl = read_csv("cylinder_fe_steady.csv")
    assert np.all(np.diff(cyl["CD"]) < 0)  # C_D falls with Re (confined, our geometry: qualitative against unbounded data)
    assert r["CD"] == pytest.approx(cyl["CD"][0], rel=0.02)  # coarse live vs medium cached mesh


def test_fe_time_V3_10_163_orders_one_and_two_N97():  # V3 (C13/N97): α = 1, β = 0 and α = 2, β = 1
    assert abs(FD.ode_scheme_errors("alpha_beta", alpha_t=1.0, beta_t=0.0)["order"] - 1) < 0.1
    assert abs(FD.ode_scheme_errors("alpha_beta", alpha_t=2.0, beta_t=1.0)["order"] - 2) < 0.1
    assert FEM2.time_derivative(3.0, 1.0, 0.5, 0.5, 2.0, 1.0) == pytest.approx(2 * 2 / 0.5 - 0.5)


def test_fe_unsteady_V1_steady_state_is_a_fixed_point_of_march_unsteady():  # V1 (C13/N110): the marching machinery
    m = FEM2.cylinder_channel_mesh(n_theta=16, h_far=1.0, n_r=4)
    bc = FEM2.channel_bcs(m)
    s = FEM2.newton_solve(m, None, 10.0, bc=bc, max_iter=20)
    assert s["converged"]
    r = FEM2.march_unsteady(m, bc, 10.0, 0.25, 3, U0=(s["u"], s["v"], s["p"]))
    f = FEM2.cylinder_forces(m, s, 10.0)
    assert np.allclose(r["CD"], f["CD"], rtol=1e-6) and np.max(np.abs(r["sol"]["u"] - s["u"])) < 1e-7


def test_fe_forces_V7_unsteady_history_strouhal_R13_N110():  # V7 (C13/R13, N110): our cached Re = 100 history
    d = read_csv("cylinder_fe_re100_forces.csv")
    late = d["t"] > d["t"][-1] / 2
    St = ch10.dominant_frequency(d["t"][late], d["CL"][late])
    St2 = ch10.dominant_frequency(d["t"][late], d["CL"][late], "fft")
    assert 0.15 < St < 0.25 and abs(St - St2) < 0.01  # shedding; confined: qualitative against unbounded data (R9)
    assert np.std(d["CL"][late]) > 0.1  # a Hopf bifurcation has happened (symmetric start)
    assert ch10.strouhal_from_period(1 / St) == pytest.approx(St)
    t = np.linspace(0, 50, 5001)
    assert ch10.dominant_frequency(t, np.sin(2 * np.pi * 0.23 * t)) == pytest.approx(0.23, rel=1e-3)
    assert ch10.dominant_frequency(t, np.sin(2 * np.pi * 0.23 * t), "fft") == pytest.approx(0.23, rel=1e-3)


def test_apply_dirichlet_V1_replace_and_penalty_agree_and_nodal_vorticity_N109():  # V1 (C13/N109)
    m = FEM2.structured_square_mesh(3)
    x, y = m.nodes[:, 0], m.nodes[:, 1]
    w = FEM2.nodal_vorticity(m, -y, x)  # rigid rotation: ω = 2
    assert np.allclose(w, 2.0)
    A, f = FEM2.assemble_saddle(m, 1.0, convection=False)
    N, V = m.N, m.V
    Ad = A.toarray()
    assert np.allclose(Ad[2 * N:, :N], Ad[:N, 2 * N:].T) and np.allclose(Ad[2 * N:, 2 * N:], 0)  # [[A, B], [Bᵀ, 0]] (10.137)
    e = FEM2.element_newton_matrices(m, (np.zeros(N), np.zeros(N), np.zeros(V)), 1.0)
    assert e["Auu"].shape == (m.T, 6, 6) and e["Bup"].shape == (m.T, 6, 3)


# ======================================================================================================================
# C14 — the lid-driven cavity benchmark (Ghia 1982, Hou 1995)
# ======================================================================================================================
@needs_ref
def test_ghia_reference_V5_table_reads_and_cross_checked_values():  # V5 (C14): reference/ch10 (SOURCES.md)
    g = ch10.ghia_centreline(100)
    assert g["y"][0] == 0 and g["y"][-1] == 1 and len(g["y"]) == 17 and np.all(np.diff(g["y"]) > 0)
    # three values cross-checked on 2026-09-30 against CMC 36 (2013) Table 1 (ref. 2 = Ghia) — see SOURCES.md
    assert g["u"][list(g["y"]).index(0.4531)] == -0.21090
    assert ch10.ghia_centreline(400)["u"][list(g["y"]).index(0.2813)] == -0.32726
    assert ch10.ghia_centreline(1000)["u"][list(g["y"]).index(0.1719)] == -0.38289
    assert ch10.ghia_vortex_centre(100) == dict(x=0.6172, y=0.7344)
    assert ch10.hou_centres()[100] == (0.6196, 0.7373) and ch10.hou_centres()[400] == (0.5608, 0.6078)


@needs_ref
def test_cavity_mac_V5_centreline_within_one_percent_of_ghia_at_64():  # V5 (C14): live MAC 64², Re = 100
    st = MAC.cavity(100.0, 64, cache=False)
    e = ch10.cavity_error_vs_ghia(st, 100)
    assert e["max_dev"] < 0.01 and e["rel_max"] < 0.01  # ≤ 1 % of the lid speed at every Ghia point
    c = st["centre"]
    gh = ch10.ghia_vortex_centre(100)
    assert abs(c["x"] - gh["x"]) < 1 / 64 and abs(c["y"] - gh["y"]) < 1 / 64
    hx, hy = ch10.hou_centres()[100]
    assert abs(c["x"] - hx) < 1 / 64 and abs(c["y"] - hy) < 1 / 64
    assert np.max(np.abs(MAC.divergence(st["u"], st["v"], st["g"]))) < 1e-11  # V4
    psi = st["psi"]
    assert max(np.max(np.abs(psi[0])), np.max(np.abs(psi[-1])), np.max(np.abs(psi[:, 0])), np.max(np.abs(psi[:, -1]))) < 1e-12
    cl = MAC.cavity_centreline(st)
    assert cl["u"][0] == 0 and cl["u"][-1] == 1.0
    cv = MAC.cavity_centreline_v(st)
    assert cv["v"][0] == 0 and cv["v"][-1] == 0 and cv["v"].max() > 0 > cv["v"].min()


@needs_ref
def test_cavity_mac_V3_deviation_from_ghia_falls_at_second_order():  # V3 (C14/C15): live 16, 24, 32
    ns = (16, 24, 32)
    dev = [ch10.cavity_error_vs_ghia(MAC.cavity(100.0, n, cache=False), 100)["max_dev"] for n in ns]
    assert dev[0] > dev[1] > dev[2]
    assert observed_order([1 / n for n in ns], dev) > 2 - 0.25  # measured 2.36 (pairwise 2.65, 1.99)


def test_cavity_mac_V7_reversed_lid_mirrors_the_flow():  # V7 (C14): u → −u(1 − x), v → v(1 − x) under U → −U
    a = MAC.cavity(100.0, 16, t_end=2.0, tol_steady=0, cache=False, lid=1.0)
    b = MAC.cavity(100.0, 16, t_end=2.0, tol_steady=0, cache=False, lid=-1.0)
    assert np.max(np.abs(b["u"] + a["u"][:, ::-1])) < 1e-13 and np.max(np.abs(b["v"] - a["v"][:, ::-1])) < 1e-13


@needs_ref
def test_cavity_cached_V5_reference_runs_against_ghia_and_hou():  # V5 (C14): our committed 64²/128² tables
    for name, Re, lim in (("cavity_mac_re100_n64.csv", 100, 0.01), ("cavity_mac_re100_n128.csv", 100, 0.01),
                          ("cavity_mac_re400_n128.csv", 400, 0.01)):
        d = read_csv(name)
        g = ch10.ghia_centreline(Re)
        dev = np.max(np.abs(np.interp(g["y"], d["y"], d["u"]) - g["u"]))
        assert dev < lim, (name, dev)
        assert dev == pytest.approx(header_value(name, "max deviation from Ghia 1982"), abs=1e-5)
    hx, hy = ch10.hou_centres()[400]
    assert abs(header_value("cavity_mac_re400_n128.csv", "x") - hx) < 0.01
    assert abs(header_value("cavity_mac_re400_n128.csv", "y") - hy) < 0.01


@needs_ref
def test_cavity_maccormack_V3_V5_reference_runs_converge_to_ghia():  # V3+V5 (C10/C14): our MacCormack tables
    # 32²: 4.14 %, 64²: 1.38 %, 128²: 0.44 % of U (max; tables regenerated with the additive Δt after the review) — converging at order
    # ≈ 1.6 and inside the 1 % benchmark tolerance only at 128² (MAC is at 0.26 % already at 64²)
    ns = (32, 64, 128)
    dev = [header_value(f"cavity_mck_re100_n{n}.csv", "max deviation from Ghia 1982") for n in ns]
    assert dev[0] > dev[1] > dev[2] and dev[2] < 0.01
    assert observed_order([1 / n for n in ns], dev) > 1.4
    g = ch10.ghia_centreline(100)
    for n, dv in zip(ns, dev):
        d = read_csv(f"cavity_mck_re100_n{n}.csv")
        assert np.max(np.abs(np.interp(g["y"], d["y"], d["u"]) - g["u"])) == pytest.approx(dv, abs=1e-4)


# ======================================================================================================================
# C15 — grid convergence, Richardson, GCI (Fig. 10.13, §10.6)
# ======================================================================================================================
def test_gci_V1_exact_on_power_law_models():  # V1 (C15/D23)
    for f0, K, p, r in ((1.0, 1.0, 2.0, 2.0), (-0.3, 5.0, 1.5, 1.5), (2.0, -0.7, 3.0, 2.0)):
        h = 0.01
        f1, f2, f3 = (f0 + K * (h * r ** k) ** p for k in (0, 1, 2))
        g = grid_convergence_index(f1, f2, f3, r)
        assert g["p"] == pytest.approx(p, rel=1e-9) and g["f_ext"] == pytest.approx(f0, rel=1e-9, abs=1e-12)
        assert g["asymptotic_ratio"] == pytest.approx(1.0, rel=0.05) and g["monotone"]
        assert ch10.richardson_three(f1, f2, f3, r) == pytest.approx(dict(p=g["p"], f_ext=g["f_ext"]))
        assert richardson(f2, f1, r, p) == pytest.approx(f0, rel=1e-9, abs=1e-12)
        assert FD.richardson_extrapolate(f2, f1, r, p) == pytest.approx(f0, rel=1e-9, abs=1e-12)
    g = FD.grid_convergence_index(1 + 0.025 ** 2, 1 + 0.05 ** 2, 1 + 0.1 ** 2)
    assert g["p"] == pytest.approx(2.0) and g["f_ext"] == pytest.approx(1.0, abs=1e-12)
    assert not FD.grid_convergence_index(1.0, 1.1, 0.85)["monotone"]  # oscillatory convergence flagged


def test_gci_V7_equal_and_opposite_differences_return_nan_not_an_exception():  # V7 (C15/D23): a degenerate ratio
    # (F3, fixed in loop 1) f = (1.0, 1.1, 1.0) gives |(f3 − f2)/(f2 − f1)| = 1 ⇒ p = 0 ⇒ r^p − 1 = 0 and a
    # ZeroDivisionError escapes from f_ext; an explainer or notebook user choosing such grids crashes instead of seeing NaN.
    g = FD.grid_convergence_index(1.0, 1.1, 1.0)
    assert not g["monotone"] and not np.isfinite(g["f_ext"])


def test_richardson_D23_V2_derivation():  # V2 (C15/D23)
    f0, K, h, r, p = sp.symbols("f0 K h r p", positive=True)
    f1, f2, f3 = f0 + K * h ** p, f0 + K * (r * h) ** p, f0 + K * (r ** 2 * h) ** p  # step 1
    assert z0(sp.powsimp(sp.expand((f3 - f2) - r ** p * (f2 - f1)), force=True))  # steps 2–3: ratio r^p
    ext = f1 + (f1 - f2) / (r ** p - 1)  # step 7
    assert z0(sp.simplify(sp.powsimp(sp.expand(ext - f0), force=True)))


@needs_ref
def test_grid_convergence_V3_cavity_psi_min_second_order_and_richardson():  # V3 (C15/C14): MAC ψ_min, our tables
    f = [header_value(f"cavity_mac_re100_n{n}.csv", "psi_min") for n in (128, 64, 32)]
    g = grid_convergence_index(*f, r=2.0)
    assert abs(g["p"] - 2) < ORDER_TOL and g["monotone"]
    assert g["gci_fine"] < 0.005  # < 0.5 % uncertainty on the 128² value
    live = [MAC.cavity(100.0, n, cache=False)["centre"]["psi_min"] for n in (32, 16)]
    assert live[0] == pytest.approx(f[2], abs=1e-6) and live[1] == pytest.approx(header_value("cavity_mac_re100_n16.csv", "psi_min"), abs=1e-6)


@needs_ref
def test_grid_convergence_V7_block_drag_monotone_lift_to_zero_not_asymptotic():  # V7 (C15/R11): our block tables, Re = 20
    cd, cl = [], []
    for n in (8, 16, 32):
        d = read_csv(f"block_forces_re20_dx{n}.csv")
        late = d["t"] > 20
        cd.append(d["CD"][late].mean())
        cl.append(abs(d["CL"][late].mean()))
    assert cd[0] < cd[1] < cd[2]  # monotone in Δx
    assert cl[0] > cl[1] > cl[2] and cl[2] < 0.02  # the symmetric flow's lift → 0 with refinement
    g = grid_convergence_index(cd[2], cd[1], cd[0], 2.0)
    # the observed order is far below MacCormack's 2 (≈ 0.43 with this time window): the heuristic wall and corner
    # closures (10.151)–(10.154) are not in the asymptotic range on Δx = 1/8 … 1/32 — an Open item, labelled qualitative
    assert 0 < g["p"] < 1.0 and g["gci_fine"] > 0.01


# ======================================================================================================================
# Remaining NOTE items and helpers
# ======================================================================================================================
def test_artificial_compressibility_V1_poiseuille_and_V4_divergence_to_zero_N52():  # V1+V4 (C10/N52): (10.95)
    r = ch10.artificial_compressibility_channel()
    assert r["converged"] and r["max_err"] < 1e-9 and r["div_final"] < 1e-8
    assert r["div_history"][0] > 1e3 * r["div_history"][-1]


def test_book_slips_V1_table_complete_and_evaluators_exist():  # V1 (§8 callouts): R1–R12
    s = ch10.book_slips()
    assert [r["id"] for r in s] == [f"R{k}" for k in range(1, 13)]
    assert all(set(r) >= {"id", "where", "printed", "correct", "evaluator"} for r in s)


def test_advected_gaussian_V2_solves_10_1_and_V4_integral_conserved_N05():  # V2+V4 (C02/N05)
    x, t, u, D, x0, s0 = sp.symbols("x t u D x0 s0", positive=True)
    s2 = s0 ** 2 + 2 * D * t
    T = s0 / sp.sqrt(s2) * sp.exp(-(x - x0 - u * t) ** 2 / (2 * s2))
    assert z0(sp.simplify(T.diff(t) + u * T.diff(x) - D * T.diff(x, 2)))
    xx = np.linspace(-10, 10, 20001)
    for tt in (0.0, 0.5, 2.0):
        v = FD.advected_gaussian(xx, tt, 0.3, 0.05, x0=0.0, s0=0.2, L=None)
        assert np.trapezoid(v, xx) == pytest.approx(0.2 * np.sqrt(2 * np.pi), rel=1e-9)
    xp = np.arange(100) / 100
    assert FD.advected_gaussian(xp, 1.0, 1.0, 0.0) == pytest.approx(FD.advected_gaussian(xp, 0.0, 1.0, 0.0), abs=1e-12)
    rhs = FD.transport_rhs(FD.advected_gaussian(xp, 0.2, 0.5, 0.01), 0.5, 0.01, 0.01)
    dTdt = (FD.advected_gaussian(xp, 0.2 + 1e-6, 0.5, 0.01) - FD.advected_gaussian(xp, 0.2 - 1e-6, 0.5, 0.01)) / 2e-6
    assert np.max(np.abs(rhs - dTdt)) < 0.02 * np.max(np.abs(dTdt))


def test_scripts_drawings_V1_every_helper_returns_an_axes_drawing():  # smoke (Part C.8): scripts/ch10_drawings.py
    import importlib
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(ROOT / "scripts"))
    try:
        dr = importlib.import_module("ch10_drawings")
    finally:
        sys.path.remove(str(ROOT / "scripts"))
    for name in ("spacetime_grid", "stencil_diagram", "characteristic_diagram", "hat_functions", "element_map",
                 "staggered_grid", "p2p1_triangle", "cavity_sketch", "block_channel_sketch", "cylinder_channel_sketch"):
        fig, ax = plt.subplots()
        getattr(dr, name)(ax)
        assert len(ax.lines) + len(ax.patches) + len(ax.texts) + len(ax.collections) > 0, name
        plt.close(fig)


PART_C = """fd_weights taylor_table_sympy stencil_taylor_coefficients test_function stencil_derivative stencil_error fd_derivative
mixed_derivative_onesided transport_rhs advected_gaussian ftcs_coefficients cfl_number diffusion_number cell_peclet
transport_1d_step solve_transport_1d error_norm convergence_study truncation_error_sympy truncation_terms propagate_error
fourier_mode amplification_factor amplification_modulus ftcs_amplification_modulus2 amplification_curve max_amplification
worst_theta is_von_neumann_stable ftcs_stable stability_verdict phase_error numerical_diffusivity ftcs2d_max_amplification
steady_cd_exact cd_layer_thickness steady_cd_fd steady_cd_discrete_exact discrete_root wiggle_indicator advect_periodic
ode_scheme_order lax_demo hat hat_basis interpolate parent_shapes map_to_element map_to_parent shape_slopes
element_matrices_linear element_force element_integrals connectivity assemble_1d assembly_trace solve_steady solve_transport
interior_stencil bilinear_form weak_residual galerkin_equations_sympy MacGrid divergence gradient curl vorticity
convective_terms predictor pressure_poisson_matrix solve_pressure pin_pressure correct project projection_stages step run
dt_limit cavity cavity_centreline cavity_centreline_v streamfunction primary_vortex_centre taylor_green channel_poiseuille
mac_projection_summary maccormack_step maccormack_advection_1d ns_fluxes pressure_isothermal ns_coefficients ns_predictor
ns_corrector maccormack_dt maccormack_dt_asymptotic cavity_density_bc cavity_coefficients weakly_compressible_step
cavity_maccormack block_wall_density block_channel body_forces p2_shape p2_shape_grad p1_shape iso_map jacobian tri_quad_7pt
integrate_element structured_square_mesh p2_node_counts euler_check stokes_cavity infsup_constant infsup_table stokes_p2p1
kovasznay_test assemble_saddle cylinder_channel_mesh channel_bcs element_newton_matrices assemble_newton_system
apply_dirichlet newton_solve time_derivative cylinder_steady weak_form_sympy weak_to_strong_sympy compressible_ns_sympy
maccormack_linear_sympy projection_sympy weak_ns_identity_sympy weak_ns_components_sympy split_linear_system
marchuk_yanenko theta_scheme_linear theta_scheme_amplification_sympy splitting_order collocated_gradient checkerboard
gradient_null_space ghia_centreline ghia_vortex_centre hou_centres cavity_error_vs_ghia strouhal_from_period
dominant_frequency rod_heating_exact artificial_compressibility_channel book_slips cfl_time_step derive_all
richardson_three""".split()


@pytest.mark.parametrize("name", PART_C)
def test_part_c_smoke_every_contract_name_is_exported_and_callable(name):  # smoke (design Part C, 177 rows)
    obj = getattr(ch10, name)
    assert callable(obj)
    assert (obj.__doc__ or "").strip(), name  # every contract function carries a docstring


def test_part_c_V1_direct_calls_of_helpers_otherwise_used_indirectly():  # V1 (Part C 1.5, 3.5, 3.9, 5.9, 5.10, 6.6)
    assert FD.stencil_derivative("sin", 1.0, 0.1) == pytest.approx((np.sin(1.1) - np.sin(0.9)) / 0.2, rel=1e-14)
    assert FD.stencil_derivative(lambda z: z ** 3, 2.0, 0.5, "central2") == pytest.approx(12.0, rel=1e-13)
    g = MAC.MacGrid(9, 7, 1.0, 0.8)
    rng = np.random.default_rng(14)
    us, vs = rng.normal(size=g.u_shape), rng.normal(size=g.v_shape)
    p = rng.normal(size=(7, 9))
    gx, gy = MAC.gradient(p, g)
    cu, cv = MAC.correct(us, vs, p, g, 0.3)
    assert np.allclose(cu, us - 0.3 * gx) and np.allclose(cv, vs - 0.3 * gy)
    u, v, _ = MAC.project(us, vs, g, 0.3)
    psi = MAC.streamfunction(u, v, g)
    assert np.allclose((psi[1:] - psi[:-1]) / g.dy, u) and np.allclose(-(psi[:, 1:] - psi[:, :-1]) / g.dx, v, atol=1e-12)
    assert max(np.max(np.abs(psi[0])), np.max(np.abs(psi[-1])), np.max(np.abs(psi[:, 0])), np.max(np.abs(psi[:, -1]))) < 1e-12
    gg = MAC.MacGrid(20, 20)
    xs = np.arange(21) * gg.dx
    X, Y = np.meshgrid(xs, xs, indexing="xy")
    c = MAC.primary_vortex_centre(-np.sin(np.pi * X) * np.sin(np.pi * Y) * (1 + 0.3 * X), gg)
    assert c["psi_min"] < -1.0 and 0.5 < c["x"] < 0.7 and abs(c["y"] - 0.5) < 1e-9
    m = FEM2.structured_square_mesh(3)
    yv = m.nodes[:, 1]
    dn = m.boundary["all"]
    uex = 0.5 * 10.0 * 8.0 * yv * (1 - yv)
    bc = FEM2.BoundaryConditions(dn, uex[dn], dn, np.zeros(dn.size), pin_pressure=0)
    s = FEM2.stokes_p2p1(m, bc, Re=10.0)  # Stokes: the parabola is exact too (its convective term vanishes)
    assert np.max(np.abs(s["u"] - uex)) < 1e-11
    A, f = FEM2.assemble_newton_system(m, (uex, 0 * uex, np.zeros(m.V)), None, 10.0)
    Ap, _ = FEM2.assemble_newton_system(m, (uex, 0 * uex, np.zeros(m.V)), None, 10.0, printed_10_172=True)
    assert (A != Ap).nnz > 0  # the R6 row differs
    A1, b1 = FEM2.apply_dirichlet(A, f, dn, np.zeros(dn.size), "replace")
    A2, b2 = FEM2.apply_dirichlet(A, f, dn, np.zeros(dn.size), "penalty")
    assert np.allclose(b1[dn], 0) and np.allclose(A1[dn[0]].toarray().ravel()[dn[0]], 1.0)
    assert A2[dn[0], dn[0]] > 1e11
    lam = 1.3  # scalar y' = −λy split as λ/2 + λ/2 (commuting): exact e^{−λt}
    errs = [abs(ch10.marchuk_yanenko([[lam / 2]], [[lam / 2]], [0.0], [0.0], [1.0], dt, int(round(1 / dt)))[0] - np.exp(-lam))
            for dt in (0.02, 0.01, 0.005)]
    assert abs(observed_order([0.02, 0.01, 0.005], errs) - 1) < 0.1
    errs = [abs(ch10.theta_scheme_linear([[lam / 2]], [[lam / 2]], [1.0], dt, int(round(1 / dt)))[0] - np.exp(-lam))
            for dt in (0.02, 0.01, 0.005)]
    assert abs(observed_order([0.02, 0.01, 0.005], errs) - 2) < 0.1


def test_step_and_run_V4_small_cavity_keeps_divergence_zero():  # V4 (C12/N71): MAC.step, MAC.run, predictor
    g = MAC.MacGrid(10, 10, walls=dict(top=1.0, bottom=0.0, left=0.0, right=0.0))
    s = MAC.new_state(g)
    dt = 0.5 * MAC.dt_limit(1.0, 0.0, 100.0, g.dx, safety=1.0)
    divs = []
    r = MAC.run(s, g, 100.0, dt, 30, callback=lambda k, st: divs.append(np.max(np.abs(MAC.divergence(st["u"], st["v"], g)))))
    assert max(divs) < 1e-12 and r["steps"] == 30 and r["state"]["t"] == pytest.approx(30 * dt)
    one = MAC.step(s, g, 100.0, dt)
    us, vs = MAC.predictor(s["u"], s["v"], g, 100.0, dt, lid=1.0)
    assert np.max(np.abs(us[-1, 1:-1])) > 0  # the lid drags the top row (ghost u_g = 2U − u)
    assert np.allclose(one["u"], MAC.project(us, vs, g, dt)[0])


# ======================================================================================================================
# Review fixes (reports/ch10_review.md Must-fixes 1–3 and follow-ups)
# ======================================================================================================================
@pytest.mark.parametrize("n", [20, 33])
def test_primary_vortex_centre_V1_exact_on_a_separable_paraboloid(n):  # V1 (C14/C15, review Must-fix 1)
    # a three-point parabola in x and in y through the discrete minimum is exact for a separable quadratic: the vertex
    # (0.6137, 0.7311) and ψ_min = −0.1 must be returned exactly, and the fitted ψ_min lies below the grid minimum
    g = MAC.MacGrid(n, n)
    xs = np.arange(n + 1) * g.dx
    X, Y = np.meshgrid(xs, xs, indexing="xy")
    psi = -0.1 + 0.8 * (X - 0.6137) ** 2 + 1.7 * (Y - 0.7311) ** 2
    c = MAC.primary_vortex_centre(psi, g)
    assert c["x"] == pytest.approx(0.6137, abs=1e-12) and c["y"] == pytest.approx(0.7311, abs=1e-12)
    assert c["psi_min"] == pytest.approx(-0.1, abs=1e-14)
    assert c["psi_min"] <= psi.min()


def test_primary_vortex_centre_V7_fitted_minimum_below_the_discrete_one_on_cavity_runs():  # V7 (C14, review Must-fix 1)
    for n in (16, 32):
        st = MAC.cavity(100.0, n, cache=False)
        assert st["centre"]["psi_min"] <= st["psi"].min() + 1e-15
    # the regenerated tables: ψ_min 32/64/128² = −0.101363, −0.102986, −0.103392 → p = 2.00, f_ext −0.103528
    f = [header_value(f"cavity_mac_re100_n{n}.csv", "psi_min") for n in (128, 64, 32)]
    assert f == pytest.approx([-0.103392, -0.102986, -0.101363], abs=2e-6)
    g = grid_convergence_index(*f, r=2.0)
    assert abs(g["p"] - 2) < 0.05 and g["f_ext"] == pytest.approx(-0.103528, abs=2e-6) and g["gci_fine"] < 0.0017


def _rod_images(x, t, D=1.0, L=1.0, K=20):
    """Image (erfc) form of the heated-rod solution, walls at 1, rod initially 0 — exact, independent of (10.199)."""
    from scipy.special import erfc
    s = 2 * np.sqrt(D * t)
    x = np.asarray(x, float)
    return sum((-1) ** k * (erfc((k * L + x) / s) + erfc(((k + 1) * L - x) / s)) for k in range(K))


def test_rod_heating_V1_series_agrees_with_the_image_solution_where_sines_vanish():  # V1 (C05/N113, review Must-fix 2)
    # x = L/5 and L/3 make sin((2m − 1)πx/L) = 0 for some m: the old stopping rule quit there (x = 0.2, t = 0.001 gave −0.110)
    assert ch10.rod_heating_exact(0.2, 0.001, 1.0, 1.0) == pytest.approx(7.7442e-06, rel=1e-4)
    assert ch10.rod_heating_exact(0.2, 0.001, 1.0, 1.0) == pytest.approx(float(_rod_images(0.2, 0.001)), abs=1e-11)
    assert ch10.rod_heating_exact(1 / 3, 0.01, 1.0, 1.0) == pytest.approx(0.018425, abs=5e-7)
    assert ch10.rod_heating_exact(1 / 3, 0.01, 1.0, 1.0) == pytest.approx(float(_rod_images(1 / 3, 0.01)), abs=1e-11)


@pytest.mark.parametrize("n", [5, 15, 21, 45])
def test_rod_heating_V7_odd_cell_grids_nonnegative_and_equal_to_images(n):  # V7 (C05/N113, review Must-fix 2)
    x = np.linspace(0, 1, n + 1)
    for t in (0.001, 0.01, 0.1):
        T = ch10.rod_heating_exact(x, t, 1.0, 1.0)
        assert T.min() >= -1e-11  # heating from 0 towards the wall temperature 1: never negative
        assert np.max(np.abs(T - _rod_images(x, t))) < 1e-11


@pytest.mark.parametrize("Re,n", [(5.0, 64), (2.0, 32)])
def test_cavity_maccormack_V7_low_Re_runs_stay_finite_with_the_additive_step(Re, n):  # V7 (C10, review Must-fix 3)
    r = MCK.cavity_maccormack(Re, n=n, t_end=1.0, cache=False, tol_steady=0)
    assert np.all(np.isfinite(r["u"])) and np.max(np.abs(r["u"])) <= 1.0 + 1e-6
    with np.errstate(all="ignore"):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            bad = MCK.cavity_maccormack(Re, n=n, t_end=1.0, cache=False, tol_steady=0, dt_rule="inviscid", check_stability=False)
    assert not np.all(np.isfinite(bad["u"]))  # the inviscid bracket alone blows up here (why the default changed)


def test_cavity_maccormack_V7_too_large_step_raises_and_additive_limits():  # V1+V7 (C10/N58, review Must-fix 3)
    with pytest.raises(ValueError):
        MCK.cavity_maccormack(100.0, n=16, t_end=0.1, cache=False, dt=0.1)
    with pytest.raises(ValueError):
        MCK.cavity_maccormack(100.0, n=16, t_end=0.1, cache=False, dt_rule="nonsense")
    dx = 1 / 32
    assert MCK.maccormack_dt_additive(1, 1, 12.5, dx, dx, 1, 1e-15, 0.8) == pytest.approx(0.8 / (2 / dx + 12.5 * np.sqrt(2) / dx))
    assert MCK.maccormack_dt_additive(0, 0, 0.0, dx, dx, 1, 0.3, 1.0) == pytest.approx(1 / (2 * 0.3 * 2 / dx ** 2))  # Heun limit
    assert MCK.maccormack_dt_additive(1, 1, 12.5, dx, dx, 1, 0.5) < MCK.maccormack_dt_additive(1, 1, 12.5, dx, dx, 1, 0.01)


def test_stability_verdict_V1_reasons_and_flags_with_diffusion_on_a_C_beta_grid():  # V1 (C04/C05/C10)
    # closed forms (ours): upwind with centred diffusion stable ⇔ |C| + 2β ≤ 1; Lax–Wendroff with diffusion ⇔ C² + 2β ≤ 1
    for C in np.linspace(0.05, 1.3, 26):
        for b in np.linspace(0.0, 0.7, 15):
            for scheme, crit in (("upwind", C + 2 * b), ("lax_wendroff", C * C + 2 * b)):
                if abs(crit - 1) < 1e-6:
                    continue
                v = FD.stability_verdict(scheme, C / 2, b)
                assert v["stable"] == (crit < 1), (scheme, C, b)
                if not v["stable"]:
                    if C > 1 + 1e-12:
                        assert v["reason"] == "unstable: C > 1, the characteristic leaves the stencil"
                    elif scheme == "upwind":
                        assert v["reason"] == "unstable: |C| + 2b > 1, the zigzag grows"
                    else:
                        assert v["reason"] == "unstable: C^2 + 2b > 1, the zigzag grows"


def test_numerical_diffusivity_V7_unsteady_schemes_require_the_courant_number():  # V7 (C05/C09)
    with pytest.raises(ValueError):
        FD.numerical_diffusivity(1.0, 0.01, scheme="upwind")
    with pytest.raises(ValueError):
        FD.numerical_diffusivity(1.0, 0.01, scheme="ftcs")
    assert FD.numerical_diffusivity(1.0, 0.01) == pytest.approx(0.005)  # steady (10.94) needs no C


# ======================================================================================================================
# V6 — book values (private json; skipped when absent). The public report states relative differences only.
# ======================================================================================================================
@book_only
def test_book_values_V6_section_10_2_numbers():  # V6 (C01–C05)
    b = book()["sec_10_2_fd"]
    E = FD.truncation_error_sympy("ftcs")
    s = E["symbols"]
    tc = b["truncation_coefficients_10_17"]
    assert float(E["E"].coeff(s["T_tt"]).subs(s["dt"], 1)) == pytest.approx(tc["dt_Ttt"])
    assert float(E["E"].coeff(s["T_xxx"]).subs({s["dx"]: 1, s["u"]: 1})) == pytest.approx(tc["u_dx2_Txxx"])
    assert float(E["E"].coeff(s["T_xxxx"]).subs({s["dx"]: 1, s["D"]: 1})) == pytest.approx(tc["D_dx2_Txxxx"])
    assert FD.is_von_neumann_stable("ftcs", 0.0, b["ftcs_diffusion_limit_beta"]) and not FD.is_von_neumann_stable("ftcs", 0.0, b["ftcs_diffusion_limit_beta"] + 1e-6)
    assert FD.is_von_neumann_stable("upwind", b["cfl_limit"] / 2, 0) and not FD.is_von_neumann_stable("upwind", b["cfl_limit"] / 2 + 1e-6, 0)
    assert FD.discrete_root(b["cell_peclet_limit"] - 1e-9) > 0 > FD.discrete_root(b["cell_peclet_limit"] + 1e-9)
    assert (FD.ode_scheme_order("forward"), FD.convergence_study("ftcs")["order"]) == pytest.approx(
        (b["ftcs_order_time"], b["ftcs_order_space"]), abs=0.15)


@book_only
def test_book_values_V6_section_10_3_to_10_5_numbers():  # V6 (C07, C09, C11, C13, C14)
    bk = book()
    assert np.allclose(FEM1.interior_stencil(0.3, 1, 1)["M_over_h"], bk["sec_10_3_fem"]["consistent_mass_weights_10_63"])
    cd = bk["sec_10_4_convection_dominated"]
    assert abs(100 * np.exp(-1) - cd["T_at_one_layer_thickness_percent"]) <= 0.5  # book rounds to whole percent
    assert abs(100 * np.exp(-2) - cd["T_at_two_layer_thicknesses_percent"]) <= 0.05
    assert FD.numerical_diffusivity(1.0, 1.0) / 1.0 == pytest.approx(cd["upwind_numerical_diffusivity_factor"])
    th = bk["sec_10_4_mac_theta"]
    assert abs((1 - 1 / np.sqrt(2)) - th["theta_scheme_theta"]) < 5e-6  # printed to five digits
    mc = bk["sec_10_4_maccormack"]
    c = MCK.ns_coefficients(1.0, 1.0, 1.0, 1.0)
    assert 1 / c["c5"] == pytest.approx(mc["viscous_cross_coefficient_c5_denominator"])
    cav = bk["sec_10_5_cavity"]
    assert tuple(cav["primary_eddy_hou_re100"]) == ch10.hou_centres()[100]
    for Re, n in ((100, 128), (400, 128)):
        x, y = header_value(f"cavity_mac_re{Re}_n{n}.csv", "x"), header_value(f"cavity_mac_re{Re}_n{n}.csv", "y")
        bx, by = cav[f"primary_eddy_maccormack_re{Re}"]
        assert abs(x - bx) <= cav["primary_eddy_uncertainty"] and abs(y - by) <= cav["primary_eddy_uncertainty"]
    fe = bk["sec_10_5_cylinder_fem"]
    for mesh in ("mesh_coarse", "mesh_fine"):
        V, T = fe[mesh]["pressure_nodes"], fe[mesh]["elements"]
        E = fe[mesh]["velocity_nodes"] - V
        assert FEM2.p2_node_counts(V, E, T)["euler"] == 0  # one hole (the cylinder)
    assert ch10.strouhal_from_period(fe["period_nondimensional"]) == pytest.approx(fe["strouhal"], abs=0.005)
    P, W = FEM2.tri_quad_7pt()
    assert len(W) == fe["quadrature_points"]
    # our confined W = 5d cylinder (the book's geometry is private): St within 5 % of the book's — qualitative (R9)
    d = read_csv("cylinder_fe_re100_forces.csv")
    late = d["t"] > d["t"][-1] / 2
    assert rel(ch10.dominant_frequency(d["t"][late], d["CL"][late]), fe["strouhal"]) < 0.05


# ======================================================================================================================
# slow: live heavy runs and the scripts
# ======================================================================================================================
@slow
def test_scripts_V7_every_ch10_script_runs_headless():  # smoke (Part C.8): 11 scripts with --no-show (cached heavy runs)
    py = sys.executable
    for s in sorted((ROOT / "scripts").glob("ch10_*.py")):
        if s.name == "ch10_common.py":
            continue
        r = subprocess.run([py, str(s), "--no-show"], cwd=str(ROOT), capture_output=True, text=True, timeout=1800)
        assert r.returncode == 0, (s.name, r.stderr[-2000:])


@slow
def test_block_channel_V3_live_coarse_run_reproduces_the_cached_table():  # V3 (C15): Δx = 1/8 live vs reference CSV
    r = MCK.block_channel(20.0, 0.06, 0.125, 30.0, cache=False)
    d = read_csv("block_forces_re20_dx8.csv")  # the table keeps every second force sample, 6 s.f.
    k = np.array([np.argmin(np.abs(r["t"] - t)) for t in d["t"]])  # the same samples
    assert np.allclose(r["t"][k], d["t"], rtol=1e-5)
    assert np.allclose(r["CD"][k], d["CD"], rtol=1e-5, atol=1e-5) and np.allclose(r["CL"][k], d["CL"], rtol=1e-5, atol=1e-5)


@slow
@needs_ref
def test_cavity_maccormack_V5_live_32_against_ghia_and_mac():  # V5 (C10/C14): live MacCormack 32², Re = 100
    r = MCK.cavity_maccormack(100.0, 0.08, 32, 30.0, tol_steady=1e-6, cache=False)
    e = ch10.cavity_error_vs_ghia(r, 100)
    assert e["max_dev"] == pytest.approx(header_value("cavity_mck_re100_n32.csv", "max deviation from Ghia 1982"), abs=2e-4)
    assert abs(r["mass_drift"]) < 0.01
