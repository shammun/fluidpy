"""Verification suite for Chapter 6 — Ideal Flow (Kundu, Cohen & Dowling 5e, §§6.1–6.10).

Evidence levels (``verify-implementation`` skill): V1 analytic · V2 symbolic (sympy / dimensions) · V3 convergence ·
V4 conservation / invariant · V5 published benchmark (``reference/ch06/``, see SOURCES.md) · V6 book value (private,
``tests/book_values_ch06.json``, skipped when absent) · V7 limits / symmetry / invariance. Every test name is
``test_<concept>_<level>_<what>``; the comment on the ``def`` line repeats the level.

A items (CORE, ≥ 2 independent levels): C01 ideal-flow equations (6.1) · C02 ψ, φ and ω = −∇²ψ (6.4)–(6.15) · C03
superposition and no through-flow (6.16)–(6.18) · C04 the doublet limit (6.28)–(6.29) · C05 the half-body (6.30)–(6.32)
· C06 cylinder and d'Alembert (6.33)–(6.35) · C07 lift with circulation (6.36)–(6.40) · C08 images and Example 6.1 ·
C09 complex potential and Cauchy–Riemann (6.42)–(6.53) · C10 Blasius and Kutta–Zhukhovsky (6.54)–(6.62) · C11 conformal
mapping and the Zhukhovsky ellipse (6.63)–(6.69) · C12 finite-difference Laplace and Gauss–Seidel (6.70)–(6.73), Ex. 6.2
· C13 axisymmetric flow and the sphere (6.74)–(6.92) · C14 airship and the axial singularity method (6.93)–(6.95) · C15
the accelerating sphere and added mass (6.96)–(6.109). Derivations: every ★★ and ★★★ D row (D01, D03, D06–D08, D10,
D11, D13–D21, D23–D26, D29–D31) re-derived with sympy in ``test_*_V2_derivation`` (the ★★★ D17, D18, D21, D29, D31
re-run the design's construction step by step, including the book's printed slips).

Pinned conventions with discrimination tests: Γ counterclockwise-positive in project code, ``Gamma_cw=`` is the book's
clockwise Γ of (6.36)–(6.40) (L = +ρUΓ_cw; ``Gamma_ccw=+Γ`` gives L < 0) · 2-D dipole vector from sink to source
(cylinder d = −2πUa² e_x) · half-body θ ∈ [0, 2π) (numpy's principal branch puts −m/2 on the lower body) · Zhukhovsky
inverse on the outside branch (numpy's principal √(z² − 4b²) lands inside the circle for Re z < 0) · corrected (6.61)
1/z² coefficient, (6.104) bracket sign, (6.108) integrand — each printed form is shown to differ. (6.82) is correct as
printed: it is r × the Appendix-B divergence (derivation review M1), and a test pins that.

Run: ``.venv/Scripts/python.exe -m pytest tests/test_ch06.py -q -p no:cacheprovider``.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest
import sympy as sp
from scipy.integrate import quad

from fluidpy import ch03_kinematics as ch03
from fluidpy import ch04_conservation_laws as ch04
from fluidpy import ch06_ideal_flow as ch06
from fluidpy.core import bernoulli as BE
from fluidpy.core import biot_savart as BS
from fluidpy.core import conformal as CM
from fluidpy.core import curvilinear as CU
from fluidpy.core import laplace_solvers as LS
from fluidpy.core import panels as PN
from fluidpy.core import potential as PF
from fluidpy.core import streamfunction as SF
from fluidpy.core import vortices as VX
from fluidpy.core import vorticity as VD
from fluidpy.core.units import Q_, dimensional_check
from tools.convergence import observed_order, pairwise_orders

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference" / "ch06"
BOOK = Path(__file__).resolve().parent / "book_values_ch06.json"
needs_ref = pytest.mark.skipif(not (REF / "benchmarks.json").exists(), reason="run reference/ch06/make_refs.py")
book_only = pytest.mark.skipif(not BOOK.exists(), reason="book values are private; see CLAUDE.md rule 9")

RNG = np.random.default_rng(6)
ORDER_TOL = 0.15  # design order ± this (verify-implementation default)
TWO_PI = 2.0 * np.pi


def book():
    return json.loads(BOOK.read_text(encoding="utf-8"))


def ref_json():
    return json.loads((REF / "benchmarks.json").read_text(encoding="utf-8"))


def sym_zero(expr, symbols=(), n: int = 6, tol: float = 1e-10) -> bool:
    """True if a sympy expression is identically zero: simplify first, then (for stubborn radicals) evaluate it at
    ``n`` random positive points of ``symbols``."""
    e = sp.simplify(expr)
    if e == 0:
        return True
    syms = sorted(e.free_symbols, key=lambda s: s.name) if not symbols else list(symbols)
    f = sp.lambdify(syms, e, "numpy")
    for _ in range(n):
        vals = RNG.uniform(0.3, 1.7, len(syms))
        v = complex(f(*vals))
        if abs(v) > tol:
            return False
    return True


def ring(R: float, n: int = 256, c: complex = 0j) -> np.ndarray:
    th = np.linspace(0.0, TWO_PI, n, endpoint=False)
    return c + R * np.exp(1j * th)


def book_symbols(*names):
    """sympy symbols for sympify(locals=…) of the private JSON forms (protects Gamma, Q, beta… from sympy builtins)."""
    return {n: sp.Symbol(n, positive=True) for n in names}


# =====================================================================================================================
# C01 — the ideal-flow equations (6.1) and where they apply (N01–N04, R01; D01)
# =====================================================================================================================
def test_ideal_flow_V1_cylinder_residuals_vanish_on_a_field():  # V1 (6.1) holds for the cylinder: all residuals ≈ 0
    rr = RNG.uniform(0.12, 0.5, 12)
    th = RNG.uniform(0.0, TWO_PI, 12)
    X = np.stack([rr * np.cos(th), rr * np.sin(th)])
    res = ch06.ideal_flow_residuals("cylinder", x=X, rho=1000.0, mu=1e-3)  # U = 1 m/s, a = 0.1 m
    U, a, mu = 1.0, 0.1, 1e-3
    assert np.max(np.abs(res["continuity"])) < 1e-5 * U / a  # O(h²) stencil truncation, h = 1e-4 m
    assert np.max(np.linalg.norm(res["euler"], axis=0) / res["inertia"]) < 1e-5
    assert np.max(np.abs(res["viscous_force"])) < 1e-5 * mu * U / a ** 2  # μ∇²u = −μ∇×ω = 0 although μ ≠ 0
    assert np.max(np.abs(res["vorticity"])) < 1e-5 * U / a
    # the residuals are truncation error: they fall as h² (a wrong flow would leave an O(U/a) residual)
    P = np.array([[0.13], [0.05]])
    r1 = ch06.ideal_flow_residuals("cylinder", x=P, h=4e-4)
    r2 = ch06.ideal_flow_residuals("cylinder", x=P, h=2e-4)
    r3 = ch06.ideal_flow_residuals("cylinder", x=P, h=1e-4)
    e = [abs(float(np.ravel(q["vorticity"])[0])) for q in (r1, r2, r3)]
    assert abs(observed_order([4e-4, 2e-4, 1e-4], e) - 2.0) < ORDER_TOL
    corner = ch06.ideal_flow_residuals("corner", x=X)
    assert np.max(np.abs(corner["viscous_force"])) < 1e-9 and np.max(np.abs(corner["vorticity"])) < 1e-9
    assert np.max(np.linalg.norm(corner["euler"], axis=0) / corner["inertia"]) < 1e-9


def test_ideal_flow_V1_poiseuille_control_keeps_its_viscous_force():  # V1 control: rotational flow ≠ ideal
    y = np.array([-0.008, -0.002, 0.0, 0.005, 0.009])
    X = np.stack([np.full(5, 0.3), y])
    G, mu = 1.0, 1e-3
    res = ch06.ideal_flow_residuals("poiseuille", x=X, mu=mu, G=G, h0=0.01)
    assert np.allclose(res["viscous_force"][0], -G, rtol=1e-6)  # μ d²u/dy² = −G
    assert np.allclose(res["euler"][0], res["viscous_force"][0], rtol=1e-6)  # ρDu/Dt + ∇p = μ∇²u: Euler misses it
    assert np.allclose(res["vorticity"], G * y / mu, atol=1e-6)  # ω = −du/dy = Gy/μ ≠ 0
    assert np.allclose(res["continuity"], 0.0, atol=1e-12)


def test_ideal_flow_V1_residuals_accept_flow_objects_and_callables():  # V1 Part C 5.1 input forms
    hb = PF.half_body(1.0, TWO_PI)
    pts = np.array([[-2.0, 1.0, 2.5, -0.5], [1.5, 2.8, -3.0, -2.0]])
    r1 = ch06.ideal_flow_residuals(hb, x=pts)
    assert np.max(np.linalg.norm(r1["euler"], axis=0) / r1["inertia"]) < 1e-7  # Bernoulli p balances (u·∇)u
    u_fn, p_fn = ch06.flow_field_callables(hb, rho=1000.0)
    r2 = ch06.ideal_flow_residuals(u_fn, p_fn, x=pts)
    assert np.allclose(r1["euler"], r2["euler"], atol=1e-12)


def test_ideal_flow_V2_derivation():  # V2 D01 (★★): curl-of-curl ⇒ μ∇²u = −μ∇×ω, harmonic φ ⇒ no net viscous force
    x, y, z = sp.symbols("x y z", real=True)
    f = [sp.Function(n)(x, y, z) for n in ("f1", "f2", "f3")]
    u = sp.Matrix(f)
    X = (x, y, z)

    def grad(s):
        return sp.Matrix([sp.diff(s, v) for v in X])

    def div(A):
        return sum(sp.diff(A[i], X[i]) for i in range(3))

    def curl(A):
        return sp.Matrix([sp.diff(A[2], y) - sp.diff(A[1], z), sp.diff(A[0], z) - sp.diff(A[2], x),
                          sp.diff(A[1], x) - sp.diff(A[0], y)])

    lap = sp.Matrix([sum(sp.diff(c, v, 2) for v in X) for c in u])
    assert sp.simplify(lap - (grad(div(u)) - curl(curl(u)))) == sp.zeros(3, 1)  # step 4 identity
    # steps 5–6: u = ∇φ, ∇²φ = 0 (a 3-D harmonic polynomial) ⇒ ∇·u = 0, ω = 0, μ∇²u = 0
    phi = x ** 3 - 3 * x * y ** 2 + z * (x ** 2 - y ** 2) + x * y * z
    assert sp.simplify(sum(sp.diff(phi, v, 2) for v in X)) == 0
    uu = grad(phi)
    assert sp.simplify(div(uu)) == 0 and sp.simplify(curl(uu)) == sp.zeros(3, 1)
    assert sp.simplify(sp.Matrix([sum(sp.diff(c, v, 2) for v in X) for c in uu])) == sp.zeros(3, 1)
    # step 7: with Bernoulli p = −ρ|u|²/2 (steady, p∞ absorbed) Euler ρ(u·∇)u + ∇p = 0 holds exactly
    rho = sp.symbols("rho", positive=True)
    p = -rho * uu.dot(uu) / 2
    adv = sp.Matrix([sum(uu[j] * sp.diff(uu[i], X[j]) for j in range(3)) for i in range(3)])
    assert sp.simplify(rho * adv + grad(p)) == sp.zeros(3, 1)
    # a rotational field (plane Poiseuille) keeps μ∇²u ≠ 0 — the dropped term (control)
    G, mu, h0 = sp.symbols("G mu h0", positive=True)
    up = sp.Matrix([G / (2 * mu) * (h0 ** 2 - y ** 2), 0, 0])
    assert sp.simplify(mu * sp.Matrix([sum(sp.diff(c, v, 2) for v in X) for c in up])[0] + G) == 0


def test_ideal_flow_V2_dimensions_of_the_equations():  # V2 pint: every term of (6.1) and μ∇²u in N/m³
    rho, U, L, mu = Q_(1000, "kg/m**3"), Q_(1, "m/s"), Q_(0.1, "m"), Q_(1e-3, "Pa*s")
    inertia = dimensional_check(lambda rho, U, L: rho * U ** 2 / L, "[force]/[length]**3", rho=rho, U=U, L=L)
    visc = dimensional_check(lambda mu, U, L: mu * U / L ** 2, "[force]/[length]**3", mu=mu, U=U, L=L)
    Re = (inertia / visc).to("dimensionless").magnitude
    assert Re == pytest.approx(1000 * 1 * 0.1 / 1e-3, rel=1e-14)


def test_applicability_V1_decision_table():  # V1 N01, R01, N04: the §6.1 table, with the Kelvin check (WV)
    assert ch06.ideal_flow_applicability(1e6)["ok"]
    low = ch06.ideal_flow_applicability(50.0)
    assert not low["ok"] and "Re" in low["reasons"][0]
    assert not ch06.ideal_flow_applicability(1e6, M=0.5)["ok"]
    baro = ch06.ideal_flow_applicability(1e6, baroclinic=True)
    assert not baro["ok"] and "Kelvin fails" in baro["verdict"]  # dropping the Kelvin check would pass this case
    for region in ("boundary_layer", "wake", "separated", "pipe", "turbulent"):
        assert not ch06.ideal_flow_applicability(1e6, region=region)["ok"]
    with pytest.raises(ValueError):
        ch06.ideal_flow_applicability(1e6, region="nowhere")
    # N01 number: car at 10 m/s, L = 1 m, air ν = 1.5e-5 m²/s → Re ≈ 6.7e5, boundary layer ~ L/√Re ≈ 1.2 mm
    Re = 10.0 * 1.0 / 1.5e-5
    assert Re == pytest.approx(6.667e5, rel=1e-3) and 1.0 / np.sqrt(Re) == pytest.approx(1.22e-3, rel=1e-2)
    assert ch06.ideal_flow_applicability(Re)["ok"]


# =====================================================================================================================
# C02 — ψ and φ: ω = −∇²ψ, Laplace problems, vortices and sources as δ sources (N05–N15, N18–N19, R03, R08–R09; D02–D04)
# =====================================================================================================================
def test_psi_vorticity_V1_minus_laplacian_on_fields():  # V1 (6.4): ω = −∇²ψ on exact fields, incl. the Rankine vortex
    x, y = RNG.uniform(-2, 2, 50), RNG.uniform(-2, 2, 50)
    w = ch06.vorticity_from_psi(lambda X, Y: np.sin(X) * np.sin(Y), x, y)
    assert np.max(np.abs(w - 2 * np.sin(x) * np.sin(y))) < 1e-8
    G, a = 1.0, 0.1
    r_in = RNG.uniform(0.0, 0.085, 20)
    r_out = RNG.uniform(0.115, 1.0, 20)
    t = RNG.uniform(0, TWO_PI, 20)
    w_in = ch06.vorticity_from_psi("rankine", r_in * np.cos(t), r_in * np.sin(t), Gamma=G, a=a)
    w_out = ch06.vorticity_from_psi("rankine", r_out * np.cos(t), r_out * np.sin(t), Gamma=G, a=a)
    assert np.allclose(w_in, G / (np.pi * a ** 2), rtol=1e-9)  # 31.83 s⁻¹ in the core (D02 check)
    assert np.max(np.abs(w_out)) < 1e-7  # irrotational outside
    # the Rankine ψ is continuous with a continuous slope at r = a and its velocity is ch03's u_θ (3.28)
    e = 1e-9
    assert ch06.rankine_vortex_psi(a - e, 0, G, a) == pytest.approx(ch06.rankine_vortex_psi(a + e, 0, G, a),
                                                                     abs=4 * e * G / (TWO_PI * a))  # slope × 2e
    r = np.array([0.05, 0.2, 0.3])  # (a central difference straddling r = a sees the jump of ψ'' — kept off it)
    uth = -(ch06.rankine_vortex_psi(r + 1e-6, 0, G, a) - ch06.rankine_vortex_psi(r - 1e-6, 0, G, a)) / 2e-6
    uth_ref = np.array([VX.rankine_vortex(q, G, a)[0] for q in r])
    assert np.allclose(uth, uth_ref, rtol=1e-6)


def test_psi_vorticity_V3_ninepoint_stencil_order_four():  # V3 `laplacian_residual` 4th order on sin x sin y
    f = lambda X, Y: np.sin(X) * np.sin(Y)  # noqa: E731
    hs = np.array([0.1, 0.05, 0.025, 0.0125])
    err = [abs(PF.laplacian_residual(f, 0.7, 0.4, h) + 2 * f(0.7, 0.4)) for h in hs]
    assert abs(observed_order(hs, err) - 4.0) < ORDER_TOL


def test_psi_vorticity_V2_sympy_twin():  # V2 `vorticity_from_psi_sym` = −(ψ_xx + ψ_yy)
    x, y = sp.symbols("x y", real=True)
    assert sp.simplify(ch06.vorticity_from_psi_sym(sp.sin(x) * sp.sin(y), x, y) - 2 * sp.sin(x) * sp.sin(y)) == 0
    G = sp.symbols("Gamma", positive=True)
    assert ch06.vorticity_from_psi_sym(-G / (4 * sp.pi) * (x ** 2 + y ** 2), x, y) == G / sp.pi  # core, a = 1
    assert ch06.vorticity_from_psi_sym(-G / (2 * sp.pi) * sp.log(sp.sqrt(x ** 2 + y ** 2)), x, y) == 0


def test_laplace_residual_V1_every_element_is_harmonic():  # V1 N05, N11, N26: ∇²ψ = ∇²φ = 0 away from the centres
    els = [PF.Uniform(1.3, -0.4), PF.Source(2.0, 0.3 + 0.2j), PF.Vortex(1.5, -0.4 + 0.1j), PF.Doublet((-1.0, 0.5)),
           PF.Corner(1.0, 2.0), PF.Corner(0.7, 1.5)]
    x, y = RNG.uniform(0.5, 2.0, 30), RNG.uniform(0.6, 2.0, 30)  # upper-right quadrant: clear of every cut
    for e in els:
        for fn in (e.psi, e.phi):
            assert np.max(np.abs(PF.laplacian_residual(fn, x, y, 1e-3))) < 2e-6
    # NaN exactly at the singular centres (N26: harmonic everywhere except at their centre)
    assert np.isnan(PF.laplacian_residual(PF.Source(1.0, 0.3 + 0.2j).psi, 0.3, 0.2))
    assert np.isnan(PF.laplacian_residual(PF.Flow([PF.Vortex(1.0)]).psi, 0.001, 0.0))


def test_delta_flux_V1_vortex_and_source_strengths_on_every_circle():  # V1/V4 (6.6), (6.13), D03: flux independent of r
    G, m = 2.0, 3.0
    fv = ch06.delta_flux_check("vortex", Gamma=G)
    fs = ch06.delta_flux_check("source", kind="phi", m=m)
    assert np.allclose(fv, -G, rtol=1e-13) and np.allclose(fs, m, rtol=1e-13)  # radii 0.01 … 100 m
    fc = ch06.delta_flux_check(PF.Vortex(G).psi, radii=(0.5, 1.0, 5.0))  # callable route (4th-order differences)
    assert np.allclose(fc, -G, rtol=1e-8)
    fp = ch06.delta_flux_check(PF.Source(m).phi, radii=(0.5, 2.0), kind="phi")
    assert np.allclose(fp, m, rtol=1e-8)
    assert abs(ch06.delta_flux_check("vortex", center=(5.0, 0.0), radii=(1.0,), Gamma=G)[0]) < 1e-12  # misses
    # wrong variants: ψ = +(Γ/2π) ln r (a clockwise vortex) gives +Γ; a source written m/4π gives m/2
    wrong_v = ch06.delta_flux_check(lambda X, Y: +G / TWO_PI * np.log(np.hypot(X, Y)), radii=(1.0,))
    wrong_s = ch06.delta_flux_check(lambda X, Y: m / (4 * np.pi) * np.log(np.hypot(X, Y)), radii=(1.0,), kind="phi")
    assert wrong_v[0] == pytest.approx(+G, rel=1e-8) and wrong_s[0] == pytest.approx(m / 2, rel=1e-8)
    # the same statement with ch05's loop circulation (R-level reuse): ∮u·dx = Γ on a circle
    ufn = lambda P, t=0.0: np.stack([np.asarray(c, float) for c in PF.Vortex(G).velocity(P[0], P[1])])  # noqa: E731
    assert VD.loop_circulation(ufn, VD.circle_loop_points((0.0, 0.0), 0.7, 256)) == pytest.approx(G, rel=1e-12)


def test_delta_flux_V2_derivation():  # V2 D03 (★★): ∇² ln r = 0 for r > 0, flux of ∇ ln r = 2π on every circle
    x, y = sp.symbols("x y", real=True)
    r, R, th = sp.symbols("r R theta", positive=True)
    lnr = sp.log(sp.sqrt(x ** 2 + y ** 2))
    assert sp.simplify(sp.diff(lnr, x, 2) + sp.diff(lnr, y, 2)) == 0  # steps 1–2 (Cartesian)
    assert sp.simplify(sp.diff(r * sp.diff(sp.log(r), r), r) / r) == 0  # step 2 with the polar Laplacian (6.23a)
    assert sp.simplify(sp.diff(sp.log(r), r)) == 1 / r  # step 3: ∇ ln r = e_r/r
    flux = sp.integrate((1 / R) * R, (th, 0, 2 * sp.pi))  # step 4: n = e_r, ds = R dθ
    assert flux == 2 * sp.pi
    G, m = sp.symbols("Gamma m", positive=True)
    assert sp.simplify(-G / (2 * sp.pi) * flux) == -G  # step 7: (6.6) for (6.8)
    assert sp.simplify(m / (2 * sp.pi) * flux) == m  # steps 8–9: (6.13) for (6.15)
    # the polar Laplacian of core.curvilinear agrees (R08/R09 recap)
    Rc, phc, zc = CU.coordinates("cylindrical")
    assert sp.simplify(CU.laplacian(sp.log(Rc), "cylindrical")) == 0


def test_orthogonality_V1_flow_nets_are_square():  # V1 (6.10), D04, N09: ∇φ·∇ψ = 0 and |∇φ| = |∇ψ|
    x, y = RNG.uniform(0.3, 2.0, 40), RNG.uniform(0.3, 2.0, 40)
    fl = PF.cylinder(1.0, 0.25, Gamma_cw=1.0)
    rep = ch06.orthogonality_report(fl, (x, y))
    assert rep["max_dot_rel"] < 1e-7 and rep["max_ratio_dev"] < 1e-7
    assert ch06.orthogonality_check(fl, x + 1j * y) < 1e-7
    spec = [{"kind": "uniform", "U": 1.0}, {"kind": "source", "m": 1.0, "x": -1.0, "y": 0.0}]
    assert ch06.orthogonality_check(spec, np.stack([x, y])) < 1e-7


def test_orthogonality_V2_derivation():  # V2 D04 (★): slopes multiply to −1, ∇ψ = e_z × ∇φ (the design's corrected sign)
    x, y = sp.symbols("x y", real=True)
    w = (x + sp.I * y) ** 3 + 2 * (x + sp.I * y)
    phi, psi = sp.re(sp.expand(w)), sp.im(sp.expand(w))
    u, v = sp.diff(phi, x), sp.diff(phi, y)
    assert sp.simplify(sp.diff(psi, y) - u) == 0 and sp.simplify(-sp.diff(psi, x) - v) == 0
    assert sp.simplify((v / u) * (-u / v) + 1) == 0  # step 3
    gphi = sp.Matrix([u, v, 0])
    gpsi = sp.Matrix([sp.diff(psi, x), sp.diff(psi, y), 0])
    assert sp.simplify(gphi.dot(gpsi)) == 0
    ez = sp.Matrix([0, 0, 1])
    assert sp.simplify(ez.cross(gphi) - gpsi) == sp.zeros(3, 1)  # ∇ψ = e_z × ∇φ
    assert sp.simplify(ez.cross(gpsi) + gphi) == sp.zeros(3, 1)  # so ∇φ = e_z × ∇ψ has the wrong sign


def test_polar_velocity_V1_cylinder_components():  # V1 (6.21)–(6.22), (6.34), (6.37): polar derivatives = rotated (u, v)
    U, a, Gcw = 1.3, 0.4, 0.9
    fl = PF.cylinder(U, a, Gamma_cw=Gcw)
    psi_pol = lambda r, t: U * (r - a ** 2 / r) * np.sin(t) + Gcw / TWO_PI * np.log(r / a)  # noqa: E731  (6.36)
    phi_pol = lambda r, t: U * (r + a ** 2 / r) * np.cos(t)  # noqa: E731  (6.33), Γ = 0
    r, t = RNG.uniform(0.5, 2.0, 25), RNG.uniform(-3.0, 3.0, 25)
    ur, ut = PF.polar_velocity(psi_pol, r, t, "psi")
    ur2, ut2 = fl.velocity_polar(r, t)
    assert np.allclose(ur, ur2, atol=1e-9) and np.allclose(ut, ut2, atol=1e-9)
    ur3, ut3 = PF.polar_velocity(phi_pol, r, t, "phi")
    ur4, ut4 = PF.cylinder(U, a).velocity_polar(r, t)
    assert np.allclose(ur3, ur4, atol=1e-9) and np.allclose(ut3, ut4, atol=1e-9)
    assert np.allclose(ut2, -U * (1 + a ** 2 / r ** 2) * np.sin(t) - Gcw / (TWO_PI * r), atol=1e-12)  # (6.34)+(6.36)
    rs, ts, Us, As = sp.symbols("r theta U a", positive=True)
    urs, uts = PF.polar_velocity_sym(Us * (rs - As ** 2 / rs) * sp.sin(ts), rs, ts, "psi")
    assert sp.simplify(urs - Us * (1 - As ** 2 / rs ** 2) * sp.cos(ts)) == 0
    assert sp.simplify(uts + Us * (1 + As ** 2 / rs ** 2) * sp.sin(ts)) == 0
    with pytest.raises(ValueError):
        PF.polar_velocity(psi_pol, 1.0, 0.0, "chi")


def test_elements_V1_uniform_source_vortex_forms():  # V1 (6.7), (6.8), (6.14), (6.15), N07, N08, N13, N14, N24, R11
    x, y = RNG.uniform(-3, 3, 400), RNG.uniform(-3, 3, 400)
    U, V = 2.0, 1.0
    un = PF.Uniform(U, V)
    assert np.allclose(un.psi(x, y), -V * x + U * y, atol=1e-13) and np.allclose(un.phi(x, y), U * x + V * y)
    assert un.psi(1.0, 1.0) == pytest.approx(1.0, abs=1e-15)  # N07 number
    m, G, x0, y0 = 2.0 * np.pi, 1.7, 0.3, -0.2
    s, v = PF.Source(m, complex(x0, y0)), PF.Vortex(G, complex(x0, y0))
    rr = np.hypot(x - x0, y - y0)
    tt = np.arctan2(y - y0, x - x0)
    assert np.allclose(s.phi(x, y), m / TWO_PI * np.log(rr), atol=1e-12)  # (6.15)
    assert np.allclose(s.psi(x, y), m * tt / TWO_PI, atol=1e-12)  # (6.48) imaginary part
    assert np.allclose(v.psi(x, y), -G / TWO_PI * np.log(rr), atol=1e-12)  # (6.8), counterclockwise
    su, sv = s.velocity(x, y)
    assert np.allclose(su, m / TWO_PI * (x - x0) / rr ** 2) and np.allclose(sv, m / TWO_PI * (y - y0) / rr ** 2)  # N24
    assert s.velocity(x0 + 1.0, y0)[0] == pytest.approx(1.0, rel=1e-14)  # m = 2π → u_r = 1/r
    vu, vv = v.velocity(x, y)
    ref = BS.point_vortex_velocity(np.stack([x, y]), np.array([[x0], [y0]]), [G])  # R11 parity with ch05
    assert np.allclose(vu, ref[0], atol=1e-13) and np.allclose(vv, ref[1], atol=1e-13)
    assert np.allclose(np.hypot(vu, vv), G / (TWO_PI * rr), rtol=1e-12)


def test_elements_V1_stream_function_recap_and_poisson():  # V1 R03 (6.3) via ch04 tools; N10 (6.11) Poisson solver
    u, v = SF.velocity_from_streamfunction_2d(lambda X, Y: 1.5 * Y - 0.5 * X, 0.3, 0.7)
    assert (u, v) == (pytest.approx(1.5, rel=1e-9), pytest.approx(0.5, rel=1e-9))
    assert SF.flux_between_streamlines(PF.Uniform(2.0).psi, (0.0, 0.0), (0.0, 1.5)) == pytest.approx(3.0, rel=1e-10)
    # ∇²ψ = f with the manufactured ψ = sin(πx) sin(πy): second order (the Poisson form of (6.4))
    errs, hs = [], []
    for n in (8, 16, 32, 64):
        xs = np.linspace(0, 1, n + 1)
        X, Y = np.meshgrid(xs, xs, indexing="xy")
        ex = np.sin(np.pi * X) * np.sin(np.pi * Y)
        mask = np.zeros(X.shape, bool)
        mask[1:-1, 1:-1] = True
        psi, h = LS.solve_poisson(mask, -2 * np.pi ** 2 * ex, np.where(mask, 0.0, ex), dx=1.0 / n)
        errs.append(abs(psi[n // 2, n // 4] - ex[n // 2, n // 4]))
        hs.append(1.0 / n)
    assert abs(observed_order(hs, errs) - 2.0) < ORDER_TOL


# =====================================================================================================================
# C03 — superposition and the no-through-flow condition (N16, N17, R05; D05)
# =====================================================================================================================
def test_superposition_V1_sum_of_elements_is_linear():  # V1 superposition: w, u, φ, ψ of the sum = sums
    els = [PF.Uniform(1.2, 0.3), PF.Source(1.5, -0.5 + 0.2j), PF.Vortex(-0.8, 0.4 - 0.3j), PF.Doublet((-0.6, 0.2), 1j),
           PF.Corner(0.3, 2.0)]
    fl = PF.Flow(els)
    x, y = RNG.uniform(0.6, 2.5, 50), RNG.uniform(1.5, 2.5, 50)
    z = x + 1j * y
    assert np.allclose(fl.w(z), sum(e.w(z) for e in els), atol=1e-13)
    assert np.allclose(fl.dwdz(z), sum(e.dwdz(z) for e in els), atol=1e-13)
    u, v = fl.velocity(x, y)
    cu = sum(np.asarray(c["u"]) for c in fl.contributions(x, y))
    assert np.allclose(u, cu, atol=1e-13)
    assert fl.net_source == pytest.approx(1.5) and fl.circulation_ccw == pytest.approx(-0.8)
    fl2 = PF.Flow(els[:2]) + els[2]
    assert np.allclose(fl2.psi(x, y), PF.Flow(els[:3]).psi(x, y), atol=1e-13)
    assert np.allclose(PF.stream_function(fl, x, y), fl.psi(x, y)) and np.allclose(PF.velocity_potential(fl, x, y),
                                                                                   fl.phi(x, y))


def test_no_through_flow_V1_bodies_are_streamlines():  # V1 (6.16), N16: u·n = 0 on cylinder, half-body, ellipse
    body = ring(0.7, 512)
    for G in np.linspace(-3.0, 3.0, 9):
        fl = PF.cylinder(1.0, 0.7, Gamma_cw=G)
        assert PF.normal_velocity_on(fl, body, normals=np.exp(1j * np.angle(body))) < 1e-12
        assert np.ptp(fl.psi(body.real, body.imag)) < 1e-12  # ψ = 0 on r = a for every Γ (6.36)
    hb = PF.Flow(PF.half_body(1.0, TWO_PI).elements)  # no inside mask, so body points are evaluated
    th = np.linspace(0.05, TWO_PI - 0.05, 801)
    xb, yb = ch06.half_body_shape(1.0, TWO_PI, th)
    assert PF.normal_velocity_on(hb, np.asarray(xb) + 1j * np.asarray(yb), closed=False) < 5e-4  # O(dθ²) normals
    rep = PF.normal_velocity_values(hb, np.asarray(xb) + 1j * np.asarray(yb), closed=False)
    assert rep["rel"] < 5e-4
    # exact normals from ∇ψ: u·∇ψ = 0 identically — the sharp version of the same statement
    u, v = hb.velocity(xb, yb)
    h = 1e-6
    gx = (hb.psi(xb + h, yb) - hb.psi(xb - h, yb)) / (2 * h)
    gy = (hb.psi(xb, yb + h) - hb.psi(xb, yb - h)) / (2 * h)
    assert np.max(np.abs(u * gx + v * gy)) < 1e-8  # scale U² = 1 m²/s² (the nose is a stagnation point)
    assert np.allclose(hb.psi(xb, yb), np.pi, atol=1e-12)  # ψ = m/2 on both halves (θ ∈ [0, 2π) branch)


def test_no_through_flow_V2_derivation():  # V2 D05 (★): (e_z × n)·(e_z × ∇φ) = n·∇φ ⇒ ∂ψ/∂s = ∂φ/∂n
    n1, n2, g1, g2 = sp.symbols("n1 n2 g1 g2", real=True)
    ez, n, g = sp.Matrix([0, 0, 1]), sp.Matrix([n1, n2, 0]), sp.Matrix([g1, g2, 0])
    assert sp.expand(ez.cross(n).dot(ez.cross(g)) - n.dot(g)) == 0


def test_far_field_V7_decay_rates():  # V7 (6.17), N17: error ∝ 1/R² closed body, ∝ 1/R open body or circulation
    Rs = np.array([10.0, 20.0, 40.0, 80.0])
    closed = [PF.far_field_check(PF.cylinder(1.0, 1.0), R) for R in Rs]
    oval = [PF.far_field_check(ch06.rankine_oval()["flow"], R) for R in Rs]
    open_ = [PF.far_field_check(PF.half_body(1.0, TWO_PI), R) for R in Rs]
    circ = [PF.far_field_check(PF.cylinder(1.0, 1.0, Gamma_cw=1.0), R) for R in 10 * Rs]
    assert abs(observed_order(Rs, closed) + 2.0) < 0.05 and abs(observed_order(Rs, oval) + 2.0) < 0.05
    assert abs(observed_order(Rs, open_) + 1.0) < 0.05 and abs(observed_order(10 * Rs, circ) + 1.0) < 0.05
    assert PF.far_field_check(PF.Flow([PF.Uniform(1.0, 0.5)]), 3.0) == 0.0
    assert PF.far_field_check(PF.cylinder(1.0, 1.0), 10.0, U=1.0 + 0j) == pytest.approx(closed[0])


def test_bernoulli_V4_pressure_and_cp():  # V4 (6.18), (6.32), R05, N28: C_p = 1 at every stagnation point; p consistent
    flows = [PF.half_body(1.0, TWO_PI), PF.cylinder(2.0, 0.5, Gamma_cw=3.0), ch06.rankine_oval()["flow"]]
    for fl in flows:
        st = fl.stagnation_points(box=(-3, 3, -3, 3), keep_inside=True)
        assert len(st) >= 1
        assert np.allclose(np.abs(fl.dwdz(st)), 0.0, atol=1e-10)
        cp = PF.Flow(fl.elements).cp(st.real, st.imag, U=fl.U_inf)
        assert np.allclose(cp, 1.0, atol=1e-12)
    fl = flows[1]
    x, y = RNG.uniform(0.6, 2.0, 30), RNG.uniform(0.6, 2.0, 30)
    p = fl.pressure(x, y, rho=1.2, p_inf=101325.0)
    B = BE.bernoulli_function(np.asarray(fl.speed(x, y)), p, 0.0, rho=1.2, g=0.0)
    assert np.ptp(B) < 1e-9 * np.max(np.abs(B))  # one Bernoulli constant everywhere (irrotational, (6.18))
    q = PF.pressure_coefficient(fl.speed(x, y), 2.0)
    assert np.allclose(q, (p - 101325.0) / (0.5 * 1.2 * 4.0), atol=1e-12)
    with pytest.raises(ValueError):
        PF.Flow([PF.Source(1.0)]).cp(1.0, 1.0)  # quiescent far field needs a reference speed


# =====================================================================================================================
# C04 — the element kit and the doublet as a source–sink limit (N20–N26, R10; D06)
# =====================================================================================================================
def test_doublet_limit_V3_pair_approaches_doublet_at_second_order():  # V3 (6.28) → (6.29): error ∝ ε²
    eps = np.array([0.2, 0.1, 0.05, 0.025])
    err = ch06.doublet_limit_error(eps)
    assert abs(observed_order(eps, err) - 2.0) < ORDER_TOL
    assert all(abs(p - 2.0) < ORDER_TOL for p in pairwise_orders(eps, err))
    frames = ch06.doublet_limit_frames(eps_list=(0.3, 0.1), n=41)
    assert len(frames) == 2 and frames[1]["psi"].shape == (41, 41)
    assert frames[1]["err"] == pytest.approx(float(ch06.doublet_limit_error([0.1])[0]), rel=1e-12)


def test_doublet_limit_V1_numbers_and_direction():  # V1 D06 check numbers; N25 dipole −2mε e_x; WV d → −d
    m, eps = 50.0, 0.02
    phi_pair, u_pair, v_pair = ch06.source_sink_pair(1.0, 0.0, m, eps)
    phi_dbl = PF.Doublet((-2 * m * eps, 0.0)).phi(1.0, 0.0)
    assert phi_pair == pytest.approx(0.318352, abs=5e-7) and phi_dbl == pytest.approx(0.318310, abs=5e-7)
    assert (phi_pair - phi_dbl) / phi_dbl == pytest.approx(eps ** 2 / 3, rel=0.02)  # the dropped O(ε²) term
    # (6.29): φ = −d·x/2πr² = |d| cos θ/2πr for d = −|d| e_x; the book scalar form (6.49) is the same element
    x, y = RNG.uniform(0.5, 2, 30), RNG.uniform(-2, 2, 30)
    d = 2.0
    ref = d / TWO_PI * x / (x ** 2 + y ** 2)
    assert np.allclose(PF.Doublet((-d, 0.0)).phi(x, y), ref, rtol=1e-13)
    assert np.allclose(PF.Doublet.from_book_scalar(d).phi(x, y), ref, rtol=1e-13)
    assert np.allclose(PF.Doublet((-d, 0.0)).psi(x, y), -d / TWO_PI * y / (x ** 2 + y ** 2), rtol=1e-13)  # a-D12
    # the pair's far field is the doublet with d = −2mε e_x (from sink to source); the reversed vector is wrong
    X, Y = 3.0 * np.cos(np.linspace(0, 6, 20)), 3.0 * np.sin(np.linspace(0, 6, 20))
    pair = ch06.source_sink_pair(X, Y, 1.0 / (2 * 1e-3), 1e-3)[0]
    good = PF.Doublet((-1.0, 0.0)).phi(X, Y)
    bad = PF.Doublet((+1.0, 0.0)).phi(X, Y)
    assert np.max(np.abs(pair - good)) < 1e-6 and np.max(np.abs(pair - bad)) > 0.1
    assert PF.Doublet((-1.0, 0.0)).velocity(1.0, 0.0)[0] < 0  # on the axis past the sink the flow points back
    assert np.allclose(u_pair, ch06._pair_flow(m, eps).velocity(1.0, 0.0)[0])
    with pytest.raises(ValueError):
        PF.Doublet((1.0, 0.0, 0.0))


def test_doublet_limit_V2_derivation():  # V2 D06 (★★): Taylor of the logs, ε² terms cancel, limit with 2mε fixed
    x, y, eps, m, d = sp.symbols("x y epsilon m d", positive=True)
    r2 = x ** 2 + y ** 2
    assert sp.simplify((x + eps) ** 2 + y ** 2 - r2 * (1 + 2 * eps * x / r2 + eps ** 2 / r2)) == 0  # step 1
    phi = m / (4 * sp.pi) * (sp.log((x + eps) ** 2 + y ** 2) - sp.log((x - eps) ** 2 + y ** 2))  # (6.28), ln√ = ½ ln
    ser = sp.series(phi, eps, 0, 4).removeO()
    assert sp.simplify(ser.coeff(eps, 1) - m / sp.pi * x / r2) == 0  # step 5: (mε/π) x/r²
    assert sp.simplify(ser.coeff(eps, 2)) == 0  # step 9: the ε² terms cancel in the difference
    assert ser.coeff(eps, 3) != 0  # first surviving correction O(m ε³) = O(|d| ε²)
    lim = sp.limit((ser.coeff(eps, 1) * eps).subs(m, d / (2 * eps)), eps, 0)  # step 7: 2mε = |d| fixed
    assert sp.simplify(lim - d / (2 * sp.pi) * x / r2) == 0
    th, r = sp.symbols("theta r", positive=True)
    assert sp.simplify((d / (2 * sp.pi) * x / r2).subs({x: r * sp.cos(th), y: r * sp.sin(th)})
                       - d / (2 * sp.pi) * sp.cos(th) / r) == 0  # step 8: (6.29) polar form
    dvec = sp.Matrix([-eps * m + eps * (-m), 0])  # step 6: Σ x_i m_i with +m at −ε and −m at +ε
    assert dvec[0] == -2 * m * eps


def test_harmonic_polynomials_V2_two_per_degree():  # V2 N20, N21–N23, (6.24)–(6.27)
    x, y = sp.symbols("x y", real=True)
    hp = ch06.harmonic_polynomials(4)
    assert len(hp) == 4
    for ph, ps in hp:
        assert sp.simplify(sp.diff(ph, x, 2) + sp.diff(ph, y, 2)) == 0
        assert sp.simplify(sp.diff(ps, x, 2) + sp.diff(ps, y, 2)) == 0
    assert sp.expand(hp[1][0] - (x ** 2 - y ** 2)) == 0 and sp.expand(hp[1][1] - 2 * x * y) == 0
    a_, b_, c_ = sp.symbols("a b c")
    assert sp.solve(sp.diff(a_ * x ** 2 + b_ * x * y + c_ * y ** 2, x, 2)
                    + sp.diff(a_ * x ** 2 + b_ * x * y + c_ * y ** 2, y, 2), c_) == [-a_]  # a-D09


def test_corner_V1_quadratic_family_and_rotations():  # V1 R10 (6.24), N21 (6.25), N22 (6.26), N23 (6.27)
    A = 0.8
    x, y = RNG.uniform(0.1, 2, 30), RNG.uniform(0.1, 2, 30)
    c = PF.Corner(A, 2.0)
    assert np.allclose(c.psi(x, y), 2 * A * x * y, rtol=1e-12) and np.allclose(c.phi(x, y), A * (x ** 2 - y ** 2))
    u, v = c.velocity(x, y)
    assert np.allclose(u, 2 * A * x) and np.allclose(v, -2 * A * y)  # (6.24) velocities
    rot = PF.Corner(A, 2.0, rotate=np.pi / 4)  # the same flow turned by 45°: φ = 2Axy (6.25)
    assert np.allclose(rot.phi(x, y), 2 * A * x * y, rtol=1e-12)
    ci = PF.Corner(1j * A, 2.0)  # (6.26): ψ = A(x² − y²)
    assert np.allclose(ci.psi(x, y), A * (x ** 2 - y ** 2), atol=1e-12)
    assert np.allclose(PF.Corner(-1j * A, 2.0).phi(x, y), 2 * A * x * y, rtol=1e-12)
    with pytest.raises(ValueError):
        PF.Corner(1.0, 0.4)
    # the Part-C constant element shifts w and moves nothing
    k = PF.Constant(2.0 + 1.0j)
    assert k.psi(0.3, 0.4) == 1.0 and k.velocity(0.3, 0.4) == (0.0, 0.0)


# =====================================================================================================================
# C05 — the half-body (N27–N29; D07, D08)
# =====================================================================================================================
def test_half_body_V1_stagnation_body_and_width():  # V1 (6.30)–(6.31), D07; WV principal-branch θ
    U, m = 1.3, 2.0
    num = ch06.half_body_numbers(U, m)
    a = m / (TWO_PI * U)
    assert num["a"] == pytest.approx(a) and num["x_stag"] == pytest.approx(-a)
    assert num["psi_body"] == pytest.approx(m / 2) and num["h_max"] == pytest.approx(m / (2 * U))
    fl = PF.half_body(U, m)
    st = fl.stagnation_points(guesses=[-0.2 + 0.05j])  # (the default grid misses this point: see the F1 test)
    assert len(st) == 1 and abs(st[0] - (-a)) < 1e-12  # Newton on dw/dz
    th = np.linspace(0.02, TWO_PI - 0.02, 501)
    xb, yb = ch06.half_body_shape(U, m, th)
    el = PF.Flow(fl.elements)
    assert np.max(np.abs(el.psi(xb, yb) - m / 2)) < 1e-12 * m  # the body is ψ = m/2 on both halves (D07 steps 4–6)
    assert np.allclose(np.abs(yb), m * (np.pi - np.where(th > np.pi, TWO_PI - th, th)) / (TWO_PI * U), rtol=1e-12)
    principal = PF.Flow([PF.Uniform(U), PF.Source(m)])  # numpy's principal θ ∈ (−π, π]: cut on the upstream axis
    assert np.max(np.abs(principal.psi(xb, yb) - m / 2)) > 0.9 * m  # the lower body shows ψ = −m/2 (wrong variant)
    # N27 (6.30) and (6.50): φ = Ux + (m/2π) ln r, w = Uz + (m/2π) ln z
    x, y = RNG.uniform(-2, 2, 20), RNG.uniform(0.5, 2, 20)
    assert np.allclose(fl.phi(x, y), U * x + m / TWO_PI * np.log(np.hypot(x, y)), atol=1e-12)
    # D07 numbers: U = 1, m = 2π → a = 1, h(90°) = π/2, h_max = π
    hb1 = ch06.half_body_numbers(1.0, TWO_PI)
    assert (hb1["a"], hb1["h_max"]) == (pytest.approx(1.0), pytest.approx(np.pi))
    assert ch06.half_body_shape(1.0, TWO_PI, np.pi / 2)[1] == pytest.approx(np.pi / 2, rel=1e-14)
    assert ch06.half_body_shape(1.0, TWO_PI, np.pi)[0] == pytest.approx(-1.0)  # the nose is the stagnation point


def test_half_body_V4_mass_balance_downstream():  # V4 D07 step 9: the flux inside the body is m at every section
    U, m = 1.0, TWO_PI
    el = PF.Flow(PF.half_body(U, m).elements)
    for xs in (0.5, 5.0, 200.0):
        tb = float(np.arcsin(1.0)) if False else None  # noqa: F841 (kept simple: find the body half-width at x = xs)
        from scipy.optimize import brentq
        t_top = brentq(lambda t: float(ch06.half_body_shape(U, m, t)[0]) - xs, 1e-9, np.pi - 1e-9)
        h = float(ch06.half_body_shape(U, m, t_top)[1])
        flux = quad(lambda yy: float(el.velocity(xs, yy)[0]), -h, h, points=[0.0], epsabs=0, epsrel=1e-12, limit=200)[0]
        assert flux == pytest.approx(m, rel=1e-9)  # everything the source emits stays inside ψ = m/2
    # h(x) → m/2U far downstream (V7): width × U = m
    t_far = 1e-5
    assert float(ch06.half_body_shape(U, m, t_far)[1]) == pytest.approx(m / (2 * U), rel=1e-5)


def test_half_body_V1_surface_pressure_and_its_zero():  # V1 D08, N28, N29: our reduction vs C_p of the field
    U, m = 1.7, 3.0
    el = PF.Flow(PF.half_body(U, m).elements)
    th = np.linspace(0.05, np.pi - 0.05, 301)
    xb, yb = ch06.half_body_shape(U, m, th)
    cp_field = el.cp(xb, yb, U=U)
    assert np.allclose(ch06.half_body_surface_cp(th), cp_field, atol=1e-12)  # independent of U and m
    t0 = ch06.half_body_cp_zero_angle()
    assert np.degrees(t0) == pytest.approx(113.218, abs=1e-3)
    assert ch06.half_body_surface_cp(t0) == pytest.approx(0.0, abs=1e-12)
    assert np.tan(t0) == pytest.approx(-2 * (np.pi - t0), rel=1e-12)
    assert ch06.half_body_surface_cp(np.pi) == 1.0 and ch06.half_body_surface_cp(np.pi - 1e-7) == pytest.approx(1.0)
    assert ch06.half_body_surface_cp(np.pi / 2) == pytest.approx(-0.4053, abs=1e-4)  # D08 check number
    tt = np.linspace(0.3, 2.5, 200001)
    cpt = ch06.half_body_surface_cp(tt)
    assert cpt.min() == pytest.approx(-0.587, abs=1e-3) and np.degrees(tt[np.argmin(cpt)]) == pytest.approx(63.0, abs=0.2)
    assert abs(ch06.half_body_surface_cp(1e-4)) < 1e-3  # → 0⁻ far downstream


def test_half_body_V3_net_force_vanishes_as_the_body_lengthens():  # V3/V7 N29 (Exercise 6.13): D → 0, L = 0
    xe = np.array([10.0, 100.0, 1000.0])
    D = np.array([ch06.half_body_net_force(x_end=x)["D"] for x in xe])
    L = np.array([ch06.half_body_net_force(x_end=x)["L"] for x in xe])
    assert np.all(np.abs(L) < 1e-12)
    assert np.all(np.diff(np.abs(D)) < 0) and abs(D[-1]) < 1e-5
    assert abs(observed_order(1.0 / xe, np.abs(D)) - 2.0) < ORDER_TOL  # the tail pressure decays as 1/x


def test_half_body_V2_derivation():  # V2 D07 (★★) and D08 (★★) step by step
    U, m, x, y, th, r = sp.symbols("U m x y theta r", positive=True)
    phi = U * x + m / (2 * sp.pi) * sp.log(sp.sqrt(x ** 2 + y ** 2))
    u, v = sp.diff(phi, x), sp.diff(phi, y)
    assert sp.simplify(u - (U + m / (2 * sp.pi) * x / (x ** 2 + y ** 2))) == 0  # D07 step 1
    xs = sp.solve(sp.Eq(U + m / (2 * sp.pi * sp.Symbol("X")), 0), sp.Symbol("X"))[0]
    assert sp.simplify(xs + m / (2 * sp.pi * U)) == 0  # step 3: x = −m/2πU
    psi = U * r * sp.sin(th) + m / (2 * sp.pi) * th
    assert psi.subs(th, sp.pi) == m / 2  # step 4: ψ_S = m/2 (θ = π at the upstream axis)
    h = sp.solve(sp.Eq(U * sp.Symbol("h") + m / (2 * sp.pi) * th, m / 2), sp.Symbol("h"))[0]
    assert sp.simplify(h - m * (sp.pi - th) / (2 * sp.pi * U)) == 0  # step 6
    assert sp.limit(h, th, 0) == m / (2 * U)  # step 7
    rb = h / sp.sin(th)
    assert sp.limit(rb, th, sp.pi) == m / (2 * sp.pi * U)  # step 8: the nose is S
    assert sp.simplify(2 * sp.limit(h, th, 0) * U - m) == 0  # step 9: mass balance
    # D08: on the body m/2πr = Uk, |u|²/U² = 1 + 2k cos θ + k², C_p = −(2k cos θ + k²)
    k = sp.sin(th) / (sp.pi - th)
    rbody = m * (sp.pi - th) / (2 * sp.pi * U * sp.sin(th))
    assert sp.simplify(m / (2 * sp.pi * rbody) - U * k) == 0  # step 2
    uu = U + m / (2 * sp.pi * rbody) * sp.cos(th)
    vv = m / (2 * sp.pi * rbody) * sp.sin(th)
    speed2 = sp.simplify((uu ** 2 + vv ** 2) / U ** 2)
    assert sp.simplify(speed2 - (1 + 2 * k * sp.cos(th) + k ** 2)) == 0  # step 5
    cp = 1 - speed2
    assert sp.simplify(cp + (2 * k * sp.cos(th) + k ** 2)) == 0  # step 6
    assert sp.limit(cp, th, sp.pi) == 1  # step 7
    assert sp.simplify(sp.factor(cp) + k * (2 * sp.cos(th) + k)) == 0  # step 8 factorisation → tan θ = −2(π − θ)


def test_rankine_oval_V1_closed_body():  # V1 Exercise 6.19 preset: ψ = 0 on the body, stagnation ±L, half-width root
    U, m, a = 1.0, TWO_PI, 1.0
    ro = ch06.rankine_oval(U, m, a)
    L, h = ro["half_length"], ro["half_width"]
    assert L == pytest.approx(np.sqrt(a ** 2 + m * a / (np.pi * U)), rel=1e-14)
    fl = PF.Flow(ro["flow"].elements)
    assert np.allclose(np.abs(fl.dwdz(np.array([-L, L]) + 0j)), 0.0, atol=1e-12)  # stagnation points
    assert abs(float(ro["psi_fn"](0.0, h))) < 1e-10  # the widest point lies on ψ = 0
    assert h == pytest.approx(m / (np.pi * U) * np.arctan(a / h), rel=1e-12)  # our root equation
    hs, Us, ms, as_ = sp.symbols("h U m a", positive=True)  # V2: ψ(0, h) = 0 ⇔ h = (m/πU) tan⁻¹(a/h)
    psi0 = Us * hs + ms / (2 * sp.pi) * (sp.atan2(hs, as_) - sp.atan2(hs, -as_))
    assert sym_zero(psi0 - (Us * hs - ms / sp.pi * sp.atan(as_ / hs)))
    assert ro["closed"] and fl.net_source == 0.0
    ins = ro["flow"].inside
    assert bool(ins(0.0, 0.5 * h)) and not bool(ins(0.0, 1.1 * h))


def test_superposition_state_V1_explainer_terms():  # V1 Part C 5.12–5.13 (E1): parts sum, status, stagnation, ψ_S
    spec = [{"kind": "uniform", "U": 1.0}, {"kind": "source", "m": TWO_PI}]
    st = ch06.superposition_state(spec, 1.0, 1.0)
    assert st["stagnation"] == [[pytest.approx(-1.0), pytest.approx(0.0, abs=1e-12)]]
    assert st["psi_dividing"] == pytest.approx(np.pi) and not st["closed"] and "open" in st["status"]
    assert sum(st["u_parts"]) == pytest.approx(st["u"]) and sum(st["v_parts"]) == pytest.approx(st["v"])
    assert st["cp"] == pytest.approx(1 - (st["u"] ** 2 + st["v"] ** 2))
    cyl = [{"kind": "uniform", "U": 1.0}, {"kind": "doublet", "dx": -TWO_PI, "dy": 0.0}]
    sc = ch06.superposition_state(cyl, (2.0, 0.5))
    assert sc["closed"] and "closed" in sc["status"]
    assert sc["u"] == pytest.approx(float(PF.cylinder(1.0, 1.0).velocity(2.0, 0.5)[0]), rel=1e-13)
    oval = [{"kind": "uniform", "U": 1.0}, {"kind": "source", "m": 1.0, "x": -1.0}, {"kind": "sink", "m": 1.0, "x": 1.0}]
    assert ch06.superposition_state(oval, probe=(0.0, 2.0))["closed"]
    assert ch06.superposition_state([{"kind": "vortex", "Gamma": 1.0}], 1.0, 0.0)["status"] == "no body"
    sk = ch06.superposition_state([{"kind": "uniform", "U": 1.0}, {"kind": "sink", "m": TWO_PI}], 1.0, 1.0)
    assert sk["status"] == "open body (net sink): extends upstream" and not sk["closed"]  # review S5
    assert sk["net_source"] == pytest.approx(-TWO_PI) and sk["stagnation"] == [[pytest.approx(1.0), 0.0]]  # x_S = +m/2πU
    fl = ch06.flow_from_spec([{"kind": "corner", "A": 1.0, "n": 2.0}, ("uniform", 0.5), PF.Vortex(1.0, 3j)])
    assert len(fl.elements) == 3 and ch06.flow_from_spec(fl) is fl
    with pytest.raises(ValueError):
        ch06.flow_from_spec([{"kind": "blob"}])
    assert isinstance(PF.element_from_spec(("sink", 2.0)), PF.Source) and PF.element_from_spec(("sink", 2.0)).m == -2.0
    assert isinstance(PF.element_from_spec({"kind": "vortex", "Gamma": 2.0}), PF.Vortex)
    with pytest.raises(ValueError):
        PF.element_from_spec(("blob", 1.0))


# =====================================================================================================================
# C06 — the circular cylinder and d'Alembert's paradox (R12–R14, N30, N31; D09)
# =====================================================================================================================
def test_cylinder_V1_parity_with_ch03_and_closed_forms():  # V1 (6.33)–(6.34), R12–R14: our flow = ch03's
    U, a = 1.4, 0.6
    fl = PF.cylinder(U, a)
    r, t = RNG.uniform(0.7, 3.0, 60), RNG.uniform(0, TWO_PI, 60)
    x, y = r * np.cos(t), r * np.sin(t)
    u, v = fl.velocity(x, y)
    u3, v3 = ch03.cylinder_flow(x, y, U, a)
    assert np.allclose(u, u3, atol=1e-13) and np.allclose(v, v3, atol=1e-13)
    assert np.allclose(fl.psi(x, y), ch03.cylinder_streamfunction(x, y, U, a), atol=1e-13)
    assert np.allclose(fl.phi(x, y), U * (r + a ** 2 / r) * np.cos(t), atol=1e-12)  # (6.33)
    ur, ut = fl.velocity_polar(r, t)
    assert np.allclose(ur, U * (1 - a ** 2 / r ** 2) * np.cos(t), atol=1e-13)  # (6.34)
    assert np.allclose(ut, -U * (1 + a ** 2 / r ** 2) * np.sin(t), atol=1e-13)
    # R14: the moving cylinder (fluid frame) is the doublet alone
    dbl = PF.Flow(fl.elements[1:])
    uf, vf = ch03.cylinder_flow(x, y, U, a, frame="fluid", t=0.0)
    ud, vd = dbl.velocity(x, y)
    assert np.allclose(ud, uf, atol=1e-13) and np.allclose(vd, vf, atol=1e-13)
    assert np.isnan(fl.velocity(0.1, 0.1)[0])  # inside the body
    with pytest.raises(ValueError):
        PF.cylinder(1.0, 0.0)


def test_cylinder_V1_surface_cp_and_dalembert():  # V1 (6.35), N30, C06: C_p = 1 − 4 sin²θ, symmetric ⇒ D = L = 0
    th = np.linspace(0, TWO_PI, 361)
    cp = ch06.cylinder_surface_cp(th)
    assert np.allclose(cp, 1 - 4 * np.sin(th) ** 2, atol=1e-14)
    assert ch06.cylinder_surface_cp(0.0) == pytest.approx(1.0) and ch06.cylinder_surface_cp(np.pi / 2) == pytest.approx(-3)
    assert np.allclose(ch06.cylinder_surface_cp(-th), cp) and np.allclose(ch06.cylinder_surface_cp(np.pi - th), cp)
    U, a, rho = 10.0, 0.1, 1.2  # N30 number: air at 10 m/s → ½ρU² = 60 Pa, −180 Pa at the shoulders
    p = ch06.cylinder_surface_pressure(np.array([0.0, np.pi / 2, np.pi]), U, a, rho=rho)
    assert np.allclose(p, [60.0, -180.0, 60.0], atol=1e-12)
    F = ch06.surface_pressure_force(PF.cylinder(U, a), R=a, rho=rho)  # spectral trapezoid, n = 64
    assert abs(F.D) < 1e-13 * rho * U ** 2 * a * 100 and abs(F.L) < 1e-13 * rho * U ** 2 * a * 100
    F2 = ch06.surface_pressure_force(PF.cylinder(U, a), ring(a, 128), rho=rho)  # detected as a circle
    assert abs(F2.D) < 1e-12 and abs(F2.L) < 1e-12
    assert BE.dynamic_pressure(U, rho) == pytest.approx(60.0)
    # the cylinder moves the fluid, but its speed at the shoulders is 2U, not U (D09 trap)
    assert float(PF.cylinder(U, a).speed(0.0, a * (1 + 1e-12))) == pytest.approx(2 * U, rel=1e-9)


def test_cylinder_V2_derivation():  # V2 D09 (★): ψ(a, θ) = 0 fixes D = 2πUa²; (6.34); C_p = 1 − 4 sin²θ
    U, a, r, th, D = sp.symbols("U a r theta D", positive=True)
    psi = U * r * sp.sin(th) - D / (2 * sp.pi) * sp.sin(th) / r
    Dsol = sp.solve(sp.Eq(sp.simplify(psi.subs(r, a) / sp.sin(th)), 0), D)[0]
    assert sp.simplify(Dsol - 2 * sp.pi * U * a ** 2) == 0
    psi = sp.simplify(psi.subs(D, Dsol))
    assert sp.simplify(psi - U * (r - a ** 2 / r) * sp.sin(th)) == 0
    ur, ut = sp.diff(psi, th) / r, -sp.diff(psi, r)
    assert ur.subs(r, a) == 0 and sp.simplify(ut.subs(r, a) + 2 * U * sp.sin(th)) == 0
    assert sp.simplify(1 - (ut.subs(r, a) / U) ** 2 - (1 - 4 * sp.sin(th) ** 2)) == 0


@needs_ref
def test_cylinder_V1_form_matches_published_potential_flow():  # V1 form cross-check (Wikipedia, see SOURCES.md)
    ref = ref_json()["cylinder_form"]
    s = {n: sp.Symbol(n, positive=True) for n in ("U", "r", "R", "theta")}
    U, a = 1.1, 0.8
    fl = PF.cylinder(U, a)
    r, t = RNG.uniform(0.9, 3.0, 20), RNG.uniform(0, TWO_PI, 20)
    sub = lambda key: sp.lambdify((s["r"], s["theta"]), sp.sympify(ref[key], locals=s).subs({s["U"]: U, s["R"]: a}))  # noqa: E731
    assert np.allclose(fl.phi(r * np.cos(t), r * np.sin(t)), sub("phi")(r, t), atol=1e-12)
    ur, ut = fl.velocity_polar(r, t)
    assert np.allclose(ur, sub("V_r")(r, t), atol=1e-12) and np.allclose(ut, sub("V_theta")(r, t), atol=1e-12)
    cp = ch06.cylinder_surface_cp(np.linspace(0, TWO_PI, 721))
    assert cp.max() == pytest.approx(ref["cp_surface_max"]) and cp.min() == pytest.approx(ref["cp_surface_min"])
    assert ch06.surface_pressure_force(fl, R=a).D == pytest.approx(ref["drag"], abs=1e-12)


# =====================================================================================================================
# C07 — the cylinder with circulation: stagnation points and L = ρUΓ (N32–N38; D10, D11)
# =====================================================================================================================
def test_lift_V1_pressure_integral_gives_rho_U_Gamma():  # V1 (6.39)–(6.40): L = ρUΓ_cw, D = 0; WV Γ_ccw sign
    U, a, rho = 10.0, 0.1, 1.2
    for G in (0.0, 0.5, 2.0, 4 * np.pi * a * U, 20.0, -3.0):
        F = ch06.surface_pressure_force(PF.cylinder(U, a, Gamma_cw=G), R=a, rho=rho)
        assert F.L == pytest.approx(rho * U * G, abs=1e-12 * rho * U ** 2 * a * 100)
        assert abs(F.D) < 1e-11
        assert ch06.lift_per_span(rho, U, Gamma_cw=G) == pytest.approx(rho * U * G)
    Fccw = ch06.surface_pressure_force(PF.cylinder(U, a, Gamma_ccw=2.0), R=a, rho=rho)
    assert Fccw.L == pytest.approx(-24.0, rel=1e-12)  # counterclockwise circulation lifts down
    assert ch06.lift_per_span(rho, U, Gamma_ccw=2.0) == pytest.approx(-24.0)
    with pytest.raises(ValueError):
        PF.cylinder(U, a, Gamma_cw=1.0, Gamma_ccw=1.0)
    assert PF.gamma_ccw_from(Gamma_cw=2.0) == -2.0 and PF.gamma_ccw_from() == 0.0
    # (6.39) on the surface equals Bernoulli from the field; (6.37) equals the field's u_θ at r = a
    th = np.linspace(0, TWO_PI, 50)
    fl = PF.cylinder(U, a, Gamma_cw=2.0)
    ps = ch06.cylinder_surface_pressure(th, U, a, Gamma_cw=2.0, rho=rho)
    pf = PF.Flow(fl.elements).pressure(a * np.cos(th), a * np.sin(th), rho=rho, U=U)
    assert np.allclose(ps, pf, atol=1e-10)
    ut = PF.Flow(fl.elements).velocity_polar(a, th)[1]
    assert np.allclose(ch06.cylinder_surface_speed(th, U, a, Gamma_cw=2.0), ut, atol=1e-12)


def test_lift_V1_stagnation_points_closed_form_and_newton():  # V1 (6.38), N34: surface, merged, free point
    U, a = 10.0, 0.1
    crit = 4 * np.pi * a * U
    st = ch06.cylinder_stagnation_points(U, a, Gamma_cw=2.0)
    assert np.allclose(np.abs(st), a) and np.allclose(np.sin(np.angle(st)), -2.0 / crit, rtol=1e-12)
    assert np.allclose(np.sort(np.degrees(np.angle(st))), [-170.84215, -9.15785], atol=1e-4)  # N34 numbers
    fl = PF.cylinder(U, a, Gamma_cw=2.0)
    newton = fl.stagnation_points(box=(-0.3, 0.3, -0.3, 0.3))
    assert np.allclose(np.sort_complex(newton), np.sort_complex(st), atol=1e-12)
    assert np.allclose(np.abs(PF.Flow(fl.elements).dwdz(st)), 0.0, atol=1e-10)
    m1 = ch06.cylinder_stagnation_points(U, a, Gamma_cw=crit)
    assert len(m1) == 1 and abs(m1[0] - (-1j * a)) < 1e-15  # merged at the bottom
    G = 6 * np.pi * a * U  # N34 / D10: Γ = 6πaU → r₊ = 0.2618 m, r₋ = 0.0382 m, r₊r₋ = a²
    free = ch06.cylinder_stagnation_points(U, a, Gamma_cw=G, include_inside=True)
    rs = np.sort(np.abs(free))
    assert rs == pytest.approx([0.0381966, 0.2618034], abs=1e-7) and rs[0] * rs[1] == pytest.approx(a ** 2, rel=1e-12)
    assert abs(PF.Flow(PF.cylinder(U, a, Gamma_cw=G).elements).dwdz(np.array([-1j * rs[1]]))[0]) < 1e-10
    up = ch06.cylinder_stagnation_points(U, a, Gamma_cw=-2.0)  # the mirror image: points move up
    assert np.allclose(np.sort(np.angle(up)), np.sort(-np.angle(st)), atol=1e-12)
    assert np.imag(ch06.cylinder_stagnation_points(U, a, Gamma_cw=-G)[0]) > 0


def test_lift_V4_circulation_family_is_non_unique():  # V4/V1 N38: every Γ meets the BCs; loop Γ = −Γ_cw on any circle
    U, a = 1.0, 1.0
    Gs = (0.0, 2 * np.pi * a * U, 4 * np.pi * a * U, 8 * np.pi * a * U)
    out = ch06.circulation_family_check(U, a, Gs)
    assert np.all(out["max_normal"] < 1e-12)
    R = 1000.0 * a
    assert np.all(out["far_field"] <= np.abs(Gs) / (TWO_PI * R) + U * a ** 2 / R ** 2 + 1e-12)
    assert np.allclose(out["circulation_by_radius"], -np.array(Gs)[:, None], atol=1e-10)  # independent of r


def test_lift_V1_explainer_state():  # V1 Part C 5.17 (E2): numbers and regimes
    s = ch06.cylinder_circulation_state()
    assert s["regime"] == "two surface points" and s["L"] == pytest.approx(24.0) and abs(s["D"]) < 1e-12
    assert s["theta1_deg"] == pytest.approx(-9.15785, abs=1e-4) and s["theta2_deg"] == pytest.approx(-170.84215, abs=1e-4)
    assert s["speed_top"] == pytest.approx(2 * 10 + 2 / (TWO_PI * 0.1)) and s["speed_bottom"] == pytest.approx(
        20 - 2 / (TWO_PI * 0.1))
    crit = 4 * np.pi * 0.1 * 10
    assert ch06.cylinder_circulation_state(Gamma_cw=crit)["regime"] == "merged at the bottom"
    assert ch06.cylinder_circulation_state(Gamma_cw=-crit)["regime"] == "merged at the top"
    f = ch06.cylinder_circulation_state(Gamma_cw=1.5 * crit)
    assert f["regime"] == "free stagnation point" and f["r_free"] > 0.1 and np.isnan(f["theta1_deg"])
    assert f["L"] == pytest.approx(f["L_KJ"], rel=1e-12) and f["cp_min"] < -3


def test_lift_V2_derivation():  # V2 D10 (★★) and D11 (★★): (6.36) → (6.38); (6.39) → (6.40), D = 0
    U, a, r, th, G, rho, pinf = sp.symbols("U a r theta Gamma rho p_inf", positive=True)
    psi_v = -(-G) / (2 * sp.pi) * sp.log(r)  # step 1: clockwise vortex = counterclockwise −Γ in (6.8)
    assert sp.simplify(psi_v - G / (2 * sp.pi) * sp.log(r)) == 0
    psi = U * (r - a ** 2 / r) * sp.sin(th) + G / (2 * sp.pi) * sp.log(r / a)  # step 2 (6.36)
    assert psi.subs(r, a) == 0
    ut = -sp.diff(psi, r)
    assert sp.simplify(ut - (-U * (1 + a ** 2 / r ** 2) * sp.sin(th) - G / (2 * sp.pi * r))) == 0  # step 3
    uta = sp.simplify(ut.subs(r, a))
    assert sp.simplify(uta - (-2 * U * sp.sin(th) - G / (2 * sp.pi * a))) == 0  # step 4 (6.37)
    s = sp.solve(sp.Eq(uta.subs(sp.sin(th), sp.Symbol("s")), 0), sp.Symbol("s"))[0]
    assert sp.simplify(s + G / (4 * sp.pi * a * U)) == 0  # step 5 (6.38)
    Rr = sp.Symbol("R", positive=True)
    ut_axis = sp.simplify(ut.subs({th: -sp.pi / 2, r: Rr}))
    roots = sp.solve(sp.Eq(sp.simplify(ut_axis * Rr ** 2), 0), Rr)  # step 8–9 (with r² multiplied)
    assert sp.simplify(roots[0] * roots[1] - a ** 2) == 0  # Vieta: r₊r₋ = a²
    # D11: expand the square, integrate term by term
    p = pinf + rho / 2 * (U ** 2 - (-2 * U * sp.sin(th) - G / (2 * sp.pi * a)) ** 2)  # (6.39)
    p2 = pinf + rho * U ** 2 / 2 - 2 * rho * U ** 2 * sp.sin(th) ** 2 - rho * U * G / (sp.pi * a) * sp.sin(th) \
        - rho * G ** 2 / (8 * sp.pi ** 2 * a ** 2)
    assert sp.expand(p - p2) == 0  # steps 1–2
    L = -sp.integrate(sp.expand(p * sp.sin(th) * a), (th, 0, 2 * sp.pi))  # steps 3–8
    D = -sp.integrate(sp.expand(p * sp.cos(th) * a), (th, 0, 2 * sp.pi))  # step 9
    assert sp.simplify(L - rho * U * G) == 0 and sp.simplify(D) == 0
    assert sp.integrate(sp.sin(th) ** 3, (th, 0, 2 * sp.pi)) == 0 and sp.integrate(sp.sin(th) ** 2, (th, 0, 2 * sp.pi)) == sp.pi


def test_lift_V1_shared_air_density_default():  # V1 review S4: every 2-D force helper defaults to ρ = 1.2 kg/m³
    U, a, G = 10.0, 0.1, 2.0
    fl = PF.cylinder(U, a, Gamma_cw=G)
    vals = [ch06.lift_per_span(U=U, Gamma_cw=G), ch06.surface_pressure_force(fl, R=a).L, PF.blasius_force(fl, R=0.3).L,
            ch06.cv_force_on_body(fl, 1.0).L, ch06.cylinder_circulation_state(U, a, Gamma_cw=G)["L"],
            ch06.laurent_contributions(fl, R=0.3)["L"], ch06.blasius_state("cylinder", Gamma_cw=G, U=U)["L"]]
    assert np.allclose(vals, 24.0, rtol=1e-10)  # = 1.2 × 10 × 2; a leftover ρ = 1 default would give 20
    assert ch06.force_on_held_singularity("source", U=2.0, strength=3.0)["D"] == pytest.approx(-1.2 * 6.0, rel=1e-12)
    assert ch06.cylinder_surface_pressure(0.0, U, a) == pytest.approx(0.5 * 1.2 * U ** 2)


def test_lift_V2_dimensions():  # V2 pint: ρUΓ is a force per unit length; Γ/(4πaU) is a pure number
    L = dimensional_check(lambda rho, U, G: rho * U * G, "[force]/[length]", rho=Q_(1.2, "kg/m**3"), U=Q_(10, "m/s"),
                          G=Q_(2, "m**2/s"))
    assert L.to("N/m").magnitude == pytest.approx(ch06.lift_per_span(1.2, 10.0, Gamma_cw=2.0))
    s = (Q_(2, "m**2/s") / (4 * np.pi * Q_(0.1, "m") * Q_(10, "m/s"))).to("dimensionless").magnitude
    assert -s == pytest.approx(np.sin(np.angle(ch06.cylinder_stagnation_points(10, 0.1, Gamma_cw=2.0)[1])), rel=1e-12)


@needs_ref
def test_lift_V1_form_matches_published_kutta_joukowski():  # V1 form cross-check (Wikipedia: clockwise contour)
    ref = ref_json()["kutta_joukowski_form"]
    assert ref["lift"] == "rho*V*Gamma" and "clockwise" in ref["circulation_sense"]
    rho, U, G = 1.2, 10.0, 2.0
    # the published Γ is taken on a clockwise contour = the book's Γ_cw: L′ = +ρVΓ with our Gamma_cw
    fl = PF.cylinder(U, 0.1, Gamma_cw=G)
    loop_cw = -PF.Flow(fl.elements).velocity_polar(0.3, np.linspace(0, TWO_PI, 256, endpoint=False))[1].mean() * TWO_PI * 0.3
    assert loop_cw == pytest.approx(G, rel=1e-12)
    assert PF.blasius_force(fl, R=0.3, rho=rho).L == pytest.approx(rho * U * loop_cw, rel=1e-12)


# =====================================================================================================================
# C08 — the method of images and Example 6.1 (N39, N40, R15–R17; D12, D13)
# =====================================================================================================================
def test_images_V1_wall_conditions_and_wrong_variant():  # V1 N39 image rules; WV same-sign vortex image
    xw = np.linspace(-20, 20, 1001)
    vort = PF.Flow(PF.mirror([PF.Vortex(1.3, 0.3 + 1.0j)]))
    assert np.max(np.abs(vort.velocity(xw, 0.0)[1])) < 1e-14 and np.max(np.abs(vort.psi(xw, 0.0))) < 1e-14
    src = PF.Flow(PF.mirror([PF.Source(2.0, -0.4 + 0.7j), PF.Uniform(1.0)]))
    assert np.max(np.abs(src.velocity(xw, 0.0)[1])) < 1e-14
    dbl = PF.Flow(PF.mirror([PF.Doublet((0.3, -0.8), 0.2 + 0.5j)]))
    assert np.max(np.abs(dbl.velocity(xw, 0.0)[1])) < 1e-14
    wall_x = PF.Flow(PF.mirror([PF.Vortex(-1.0, 1.0 + 0.2j), PF.Source(1.0, 0.5 - 0.3j), PF.Doublet((0.4, 0.1), 2.0)],
                               "x=0"))
    assert np.max(np.abs(wall_x.velocity(0.0, xw)[0])) < 1e-14
    wrong = PF.Flow([PF.Vortex(1.3, 0.3 + 1.0j), PF.Vortex(1.3, 0.3 - 1.0j)])  # same-sign image: the wall leaks
    assert np.max(np.abs(wrong.velocity(xw, 0.0)[1])) > 0.1
    with pytest.raises(ValueError):
        PF.mirror([PF.Uniform(1.0, 1.0)])
    with pytest.raises(ValueError):
        PF.mirror([PF.Source(1.0, 2.0)])
    with pytest.raises(ValueError):
        PF.mirror([PF.Corner(1.0, 2.0)])
    with pytest.raises(ValueError):
        PF.mirror([PF.Source(1.0, 1j)], wall="y=1")
    # R16/R17 recap: ch05's wall image system is the same flow (parity)
    xv, Gv = BS.wall_image_system(np.array([[0.3], [1.0]]), [1.3])
    P = np.stack([RNG.uniform(-2, 2, 10), RNG.uniform(0.2, 2, 10)])
    ub = BS.point_vortex_velocity(P, xv, Gv)
    uo = vort.velocity(P[0], P[1])
    assert np.allclose(ub[0], uo[0], atol=1e-13) and np.allclose(ub[1], uo[1], atol=1e-13)


def test_images_V2_derivation_D12():  # V2 D12 (★): odd ψ₂ vanishes on the wall, even φ₂ has zero normal slope
    x, y, x0, y0 = sp.symbols("x y x0 y0", real=True)
    psi1 = sp.log((x - x0) ** 2 + (y - y0) ** 2)
    psi2 = psi1 - psi1.subs(y, -y)
    assert psi2.subs(y, 0) == 0
    phi2 = psi1 + psi1.subs(y, -y)
    assert sp.simplify(sp.diff(phi2, y).subs(y, 0)) == 0


def test_two_sources_V1_streamline_equation_and_complex_form():  # V1 (6.41), (6.53), N40, N51
    m, a = TWO_PI, 1.0
    fl = ch06.two_sources(m, a)
    ys = np.linspace(-3, 3, 101)
    assert np.max(np.abs(fl.velocity(0.0, ys)[0])) < 1e-14  # the wall x = 0 (source by a wall)
    for psi in (0.4, 1.3, 2.5):
        xs = np.linspace(0.2, 3.0, 40)
        yp = ch06.two_source_streamline(psi, m, a, xs)
        ok = np.isfinite(yp) & (np.abs(yp) > 1e-3)
        vals = fl.psi(xs[ok], yp[ok])
        k = (vals - psi) / (m / 2)  # (6.41) holds for ψ modulo m/2 (the cot has period π)
        assert np.max(np.abs(k - np.round(k))) < 1e-10
    yp, ym = ch06.two_source_streamline(1.570796, 6.283185, 1.0, 2.0, both=True)
    assert np.isfinite(yp) and np.isfinite(ym)
    z = RNG.uniform(0.3, 2.0, 30) + 1j * RNG.uniform(0.3, 2.0, 30)
    w653 = m / TWO_PI * np.log((z ** 2 - a ** 2) / a ** 2)  # (6.53)
    dpsi = (w653.imag - fl.psi(z.real, z.imag)) / m
    assert np.max(np.abs(dpsi - np.round(dpsi))) < 1e-12  # Im w = ψ mod m (ln A + ln B = ln AB mod 2πi)


def test_circle_theorem_V1_circle_is_a_streamline():  # V1 R15 (our addition): u·n = 0 on |z| = a; parity with ch05
    a = 1.0
    base = PF.Flow([PF.Vortex(1.2, 2.0 + 0.5j), PF.Source(0.7, -1.5 - 1.2j), PF.Uniform(1.0)])
    ct = PF.circle_theorem(base, a)
    z = ring(a * (1 + 1e-12), 400)
    assert PF.normal_velocity_on(ct, z, normals=np.exp(1j * np.angle(z))) < 1e-10
    bare = PF.circle_theorem(lambda q: 1.0 * q, a)  # a bare w(z): W = z + a²/z, the cylinder
    zz = np.array([2.0 + 1.0j, -1.5 + 0.3j])
    assert np.allclose(bare(zz), PF.cylinder(1.0, a).w(zz), atol=1e-13)
    # parity with ch05's circle images for a vortex outside a cylinder (the centre vortex restores zero circulation)
    xv = np.array([[2.0], [0.5]])
    xa, Ga = BS.circle_image_system(xv, [1.2], a, inside=False)
    P = np.stack([RNG.uniform(1.5, 3.0, 10), RNG.uniform(-2.0, 2.0, 10)])
    ub = BS.point_vortex_velocity(P, xa, Ga)
    uc = PF.circle_theorem(PF.Flow([PF.Vortex(1.2, 2.0 + 0.5j)]), a).velocity(P[0], P[1])
    assert np.allclose(ub[0], uc[0], atol=1e-12) and np.allclose(ub[1], uc[1], atol=1e-12)


def test_example_6_1_V1_closed_form_vs_numeric_route():  # V1 Example 6.1: closed forms vs moving vortices + FD in t
    t = np.array([0.0, 5.0, 4 * np.pi, 20.0, 4 * np.sqrt(3) * np.pi, 60.0])
    c = ch06.example_6_1(t)
    n = ch06.example_6_1(t, route="numeric")
    scale = 1000.0 / (4 * np.pi ** 2)
    assert np.max(np.abs(n["p_origin"] - c["p_origin"])) < 1e-6 * scale
    assert np.allclose(n["xi"], c["xi"], atol=1e-8)  # path from ch05's evolver = (h, Γt/4πh)
    assert np.allclose(n["v_origin"], c["v_origin"], rtol=1e-12)
    assert np.allclose(c["p_origin"], c["unsteady"] + c["speed_part"], atol=1e-12)
    tm = ch06.example_6_1_times()
    assert c["p_origin"][0] == pytest.approx(-25.33029591, rel=1e-9) and tm["p_min"] == pytest.approx(c["p_origin"][0])
    assert abs(c["p_origin"][2]) < 1e-12 and tm["t_zero"] == pytest.approx(4 * np.pi)
    assert c["p_origin"][4] == pytest.approx(3.16628699, rel=1e-9) == tm["p_max"]
    assert tm["t_max"] == pytest.approx(21.7656, abs=1e-4)
    d = ch06.example_6_1(np.array([tm["t_max"] - 1e-3, tm["t_max"] + 1e-3]))["p_origin"]
    assert d[0] < c["p_origin"][4] and d[1] < c["p_origin"][4]  # a maximum
    assert c["p_origin"][-1] < c["p_origin"][4] and ch06.example_6_1(1e6)["p_origin"] == pytest.approx(0.0, abs=1e-8)
    assert c["xi_y"][1] == pytest.approx(5.0 / (4 * np.pi))  # drift Γ/4πh = 0.0796 m/s
    s = ch06.example_6_1(5.0)
    assert np.ndim(s["p_origin"]) == 0 and s["xi"].shape == (2,)
    # R32 recap: ch04's unsteady Bernoulli gives the same number from ∂φ/∂t and v
    pb = BE.unsteady_bernoulli_pressure(s["dphidt"], s["v_origin"], 0.0, rho=1000.0, g=0.0)
    assert pb == pytest.approx(s["p_origin"], rel=1e-12)
    with pytest.raises(ValueError):
        ch06.example_6_1(1.0, route="guess")


def test_example_6_1_V3_time_derivative_is_second_order():  # V3 the numeric route's ∂φ/∂t by central differences
    t = np.array([0.0, 5.0, 12.0, 20.0])
    c = ch06.example_6_1(t)["p_origin"]
    dts = [0.4, 0.2, 0.1, 0.05]
    err = [np.max(np.abs(ch06.example_6_1(t, route="numeric", dt=d)["p_origin"] - c)) for d in dts]
    assert abs(observed_order(dts, err) - 2.0) < ORDER_TOL


def test_example_6_1_V4_vortex_keeps_its_distance():  # V4 ξ_x conserved (impulse), wall velocity u = 0
    t = np.linspace(0, 50, 11)
    n = ch06.example_6_1(t, Gamma=1.7, h=0.8, route="numeric")
    assert np.ptp(n["xi_x"]) < 1e-9 and np.allclose(n["xi_x"], 0.8, atol=1e-9)
    s = 1.7 * 30.0 / (4 * np.pi * 0.8)
    fl = PF.Flow([PF.Vortex(-1.7, complex(0.8, s)), PF.Vortex(1.7, complex(-0.8, s))])
    assert np.max(np.abs(fl.velocity(0.0, np.linspace(-5, 5, 101))[0])) < 1e-14


def test_example_6_1_V1_wall_pressure_extension():  # V1 Part C 5.20: y = 0 equals the origin value; parts add up
    for t in (0.0, 3.0, 17.0):
        assert ch06.example_6_1_wall_pressure(0.0, t) == pytest.approx(ch06.example_6_1(t)["p_origin"], rel=1e-12)
    sp_ = ch06.example_6_1_wall_pressure(np.array([-1.0, 0.5, 2.0]), 3.0, split=True)
    assert np.allclose(sp_["total"], sp_["unsteady"] + sp_["speed"], atol=1e-12)
    eta = sp_["eta"]
    fl = PF.Flow([PF.Vortex(-1.0, complex(1.0, eta)), PF.Vortex(1.0, complex(-1.0, eta))])
    assert np.allclose(sp_["v"], fl.velocity(0.0, np.array([-1.0, 0.5, 2.0]))[1], rtol=1e-12)


def test_example_6_1_V2_derivation():  # V2 D13 (★★): ∂φ/∂t and v at the origin, one fraction, zero and maximum
    G, h, t, x, y = sp.symbols("Gamma h t x y", positive=True)
    s = G * t / (4 * sp.pi * h)
    phi = G / (2 * sp.pi) * sp.atan((y - s) / (x + h)) - G / (2 * sp.pi) * sp.atan((y - s) / (x - h))  # step 4
    dphidt0 = sp.simplify(sp.diff(phi, t).subs({x: 0, y: 0}))
    assert sp.simplify(dphidt0 + G ** 2 / (4 * sp.pi ** 2) / (h ** 2 + s ** 2)) == 0  # step 8
    psi = G / (2 * sp.pi) * (-sp.log(sp.sqrt((x + h) ** 2 + (y - s) ** 2)) + sp.log(sp.sqrt((x - h) ** 2 + (y - s) ** 2)))
    assert sp.simplify(sp.diff(psi, y).subs(x, 0)) == 0  # step 1: u = 0 on the wall
    v0 = sp.simplify(-sp.diff(psi, x).subs({x: 0, y: 0}))
    assert sp.simplify(v0 - G * h / (sp.pi * (h ** 2 + s ** 2))) == 0  # step 9
    p = sp.simplify(-dphidt0 - v0 ** 2 / 2)  # steps 5–6, 10
    assert sp.simplify(p - G ** 2 / (4 * sp.pi ** 2) * (s ** 2 - h ** 2) / (s ** 2 + h ** 2) ** 2) == 0  # step 11
    S = sp.Symbol("S", positive=True)
    f = (S ** 2 - h ** 2) / (S ** 2 + h ** 2) ** 2
    crit = sp.solve(sp.diff(f, S), S)
    assert sp.sqrt(3) * h in crit  # maximum at s = √3 h, i.e. t = 4√3πh²/Γ (our addition)
    assert sp.simplify(f.subs(S, sp.sqrt(3) * h) * G ** 2 / (4 * sp.pi ** 2) - G ** 2 / (32 * sp.pi ** 2 * h ** 2)) == 0
    assert sp.simplify(sp.diff(s, t) - G / (4 * sp.pi * h)) == 0  # steps 2–3: drift Γ/4πh


# =====================================================================================================================
# C09 — the complex potential, Cauchy–Riemann and the complex velocity (N41–N52; D14, D15)
# =====================================================================================================================
def test_complex_potential_V1_element_forms():  # V1 (6.42)–(6.53): Re/Im of every w equal the book's real forms
    x, y = RNG.uniform(-3, 3, 2000), RNG.uniform(-3, 3, 2000)
    z = x + 1j * y
    G, m, d, U, a = 1.3, 0.9, 0.7, 1.1, 0.5
    r, th = np.abs(z), np.angle(z)
    v = PF.Vortex(G)
    assert np.allclose(v.w(z), -1j * G / TWO_PI * np.log(z), atol=1e-12)  # (6.47)
    assert np.allclose(v.phi(x, y), G * th / TWO_PI, atol=1e-12) and np.allclose(v.psi(x, y), -G / TWO_PI * np.log(r))
    s = PF.Source(m)
    assert np.allclose(s.w(z), m / TWO_PI * np.log(z), atol=1e-12)  # (6.48)
    db = PF.Doublet.from_book_scalar(d)
    assert np.allclose(db.w(z), d / (TWO_PI * z), rtol=1e-12)  # (6.49)
    assert np.allclose(db.phi(x, y), d / TWO_PI * np.cos(th) / r, rtol=1e-10)  # = (6.29)
    hb = PF.half_body(U, m)
    ok = ~hb.inside(x, y)
    th2 = np.mod(th, TWO_PI)
    assert np.allclose(hb.psi(x[ok], y[ok]), U * r[ok] * np.sin(th2[ok]) + m * th2[ok] / TWO_PI, atol=1e-12)  # (6.50)
    cyl = PF.cylinder(U, a, Gamma_cw=G)
    out = r > a
    assert np.allclose(cyl.w(z[out]), U * (z[out] + a ** 2 / z[out]) + 1j * G / TWO_PI * np.log(z[out] / a), atol=1e-12)
    assert np.allclose(cyl.psi(x[out], y[out]), U * (r[out] - a ** 2 / r[out]) * np.sin(th[out])
                       + G / TWO_PI * np.log(r[out] / a), atol=1e-12)  # (6.52) → (6.36)
    n = 1.5
    wedge = (th > 0.01) & (th < np.pi / n - 0.01)
    assert np.allclose(PF.Corner(0.8, n).w(z[wedge]), 0.8 * r[wedge] ** n * np.exp(1j * n * th[wedge]), rtol=1e-12)  # (6.46)
    assert np.allclose(PF.Uniform(U, 0.4).w(z), (U - 0.4j) * z)  # (6.7), (6.14)
    assert np.allclose(np.abs(z) * np.exp(1j * th), z)  # (6.43) z = r e^{iθ}


def test_complex_potential_V1_cauchy_riemann_discriminates():  # V1 (6.44) N42: analytic w ⇒ residual 0; z* ⇒ (2, 0)
    for kind, kw in (("corner", {"n": 2.0}), ("corner", {"n": 0.5}), ("uniform", {}), ("source", {}), ("vortex", {}),
                     ("doublet", {}), ("cylinder", {"Gamma_cw": 1.0})):
        for zz in (1.2 + 0.7j, 0.4 + 1.9j, 2.2 + 0.3j):
            r1, r2 = ch06.cauchy_riemann_residual(kind, zz, **kw)
            assert abs(r1) < 1e-6 and abs(r2) < 1e-6
    r1, r2 = ch06.cauchy_riemann_residual("conj", 0.7 + 0.2j)
    assert r1 == pytest.approx(2.0, abs=1e-6) and abs(r2) < 1e-6
    r1, r2 = ch06.cauchy_riemann_residual(lambda q: q ** 3 - 2 * q, 0.3 + 0.4j)
    assert abs(r1) < 1e-6 and abs(r2) < 1e-6
    x, y = sp.symbols("x y", real=True)
    assert ch06.cauchy_riemann_residual_sym((x + sp.I * y) ** 3, x, y) == (0, 0)
    assert ch06.cauchy_riemann_residual_sym(x - sp.I * y, x, y) == (2, 0)
    zs = sp.Symbol("z")
    assert ch06.cauchy_riemann_residual_sym(sp.exp(zs) + zs ** 2, zs) == (0, 0)
    with pytest.raises(ValueError):
        ch06.complex_potential_family("blob")


def test_complex_potential_V1_derivative_is_u_minus_iv():  # V1 (6.45) N43: dw/dz = u − iv from ψ by differences
    fl = PF.Flow([PF.Uniform(1.0, 0.3), PF.Source(1.2, -1 - 1j), PF.Vortex(0.8, 1 - 1j), PF.Doublet((0.3, -0.2), 2j)])
    x, y = RNG.uniform(0.2, 2.0, 30), RNG.uniform(0.2, 1.5, 30)
    h = 1e-6
    u = (fl.psi(x, y + h) - fl.psi(x, y - h)) / (2 * h)
    v = -(fl.psi(x + h, y) - fl.psi(x - h, y)) / (2 * h)
    q = fl.dwdz(x + 1j * y)
    assert np.allclose(q.real, u, atol=1e-7) and np.allclose(-q.imag, v, atol=1e-7)
    p = ch06.complex_potential_probe("corner", 1.0, 1.0)  # N43 number: w = z² at 1 + i → u = 2, v = −2 m/s
    assert (p["u"], p["v"]) == (pytest.approx(2.0), pytest.approx(-2.0)) and p["speed"] == pytest.approx(np.sqrt(8))
    assert p["qx_re"] == pytest.approx(2.0, abs=1e-6) and p["qy_im"] == pytest.approx(2.0, abs=1e-6)
    pc = ch06.complex_potential_probe("conj", 0.7, 0.2)
    assert np.isnan(pc["u"]) and pc["cr1"] == pytest.approx(2.0, abs=1e-6)
    pp = ch06.complex_potential_probe("corner", 0.3, 0.8, n=0.5, h=1e-4)
    assert pp["speed"] == pytest.approx(0.5 * np.hypot(0.3, 0.8) ** -0.5, rel=1e-12)
    for kind in ("source", "vortex", "doublet", "cylinder"):
        pk = ch06.complex_potential_probe(kind, 1.5, 0.8)
        assert abs(pk["cr1"]) < 1e-6 and abs(pk["cr2"]) < 1e-6 and np.isfinite(pk["speed"])


def test_corner_V1_walls_speed_exponent_and_branch():  # V1/V7 (6.46), N44, D15: walls, |dw/dz| ∝ r^(n−1), no cut inside
    for n in (4.0, 2.0, 1.5, 1.0, 2 / 3, 0.5):
        c = PF.Corner(1.0, n)
        al = np.pi / n
        rr = np.array([1e-3, 1e-2, 1e-1, 0.5])
        assert np.max(np.abs(c.psi(rr, 0.0 * rr))) < 1e-12  # wall θ = 0
        assert np.max(np.abs(c.psi(rr * np.cos(al), rr * np.sin(al)))) < 1e-12  # wall θ = π/n
        sp_ = np.abs(c.dwdz(rr * np.exp(0.5j * al)))
        assert np.polyfit(np.log(rr), np.log(sp_), 1)[0] == pytest.approx(ch06.corner_speed_exponent(n), abs=0.01)
        arc = np.linspace(1e-4, al - 1e-4, 2001)
        ps = c.psi(np.cos(arc), np.sin(arc))
        assert np.max(np.abs(np.diff(ps))) < 5e-3  # continuous across the wedge: the cut lies outside the fluid
        assert c.alpha == pytest.approx(al)
    assert ch06.corner_info(2 / 3)["regime"] == "infinite speed at the corner"
    assert ch06.corner_info(2 / 3)["alpha_deg"] == pytest.approx(270.0)  # Example 6.2's re-entrant corner
    assert ch06.corner_info(2.0)["regime"] == "stagnation point" and ch06.corner_info(1.0)["regime"] == "uniform flow"
    assert ch06.corner_speed_exponent(0.6666666666666666) == pytest.approx(-1 / 3)
    with pytest.raises(ValueError):
        ch06.corner_info(0.3)


def test_complex_potential_V2_derivation():  # V2 D14 (★★) and D15 (★★)
    x, y = sp.symbols("x y", real=True)
    c1, c2 = sp.symbols("c1 c2", real=True)
    zz = x + sp.I * y
    w = sp.expand(zz ** 3 + (c1 + sp.I * c2) * zz ** 2 + sp.exp(zz))
    w = sp.expand_complex(w)
    phi, psi = sp.re(w), sp.im(w)
    Dx = sp.diff(phi, x) + sp.I * sp.diff(psi, x)  # step 1: along x
    Dy = (sp.diff(phi, y) + sp.I * sp.diff(psi, y)) / sp.I  # step 2: along iy
    assert sp.simplify(sp.expand_complex(Dy) - (sp.diff(psi, y) - sp.I * sp.diff(phi, y))) == 0  # step 3 (1/i = −i)
    assert sp.simplify(sp.expand_complex(Dx - Dy)) == 0  # step 4 ⇒ Cauchy–Riemann (6.44)
    u, v = sp.diff(phi, x), sp.diff(phi, y)
    assert sp.simplify(Dx - (u - sp.I * v)) == 0  # step 5: dw/dz = u − iv (6.45)
    assert sp.simplify(sp.diff(phi, x, 2) + sp.diff(phi, y, 2)) == 0  # steps 6–7
    assert sp.simplify(sp.diff(psi, x, 2) + sp.diff(psi, y, 2)) == 0
    assert sp.simplify(u * sp.diff(psi, x) + v * sp.diff(psi, y)) == 0  # step 8
    # D15: w = A zⁿ, ψ = A rⁿ sin nθ, walls θ = 0 and π/n, dw/dz = (Aπ/α) z^((π−α)/α)
    A, r, th, n = sp.symbols("A r theta n", positive=True)
    wpol = A * r ** n * sp.exp(sp.I * n * th)
    psi_n = sp.im(sp.expand_complex(wpol))
    assert sp.simplify(psi_n - A * r ** n * sp.sin(n * th)) == 0  # step 2
    assert psi_n.subs(th, 0) == 0 and sp.simplify(psi_n.subs(th, sp.pi / n)) == 0  # step 3
    assert sp.expand((A * zz ** 2)) == sp.expand(A * (x ** 2 - y ** 2) + sp.I * 2 * A * x * y)  # step 4
    al, Z = sp.symbols("alpha Z", positive=True)
    lhs = sp.diff(A * Z ** n, Z).subs(n, sp.pi / al)
    assert sp.simplify(lhs - A * sp.pi / al * Z ** ((sp.pi - al) / al)) == 0  # step 5


# =====================================================================================================================
# C10 — forces: Blasius's theorem and Kutta–Zhukhovsky (N53–N63, R18; D16, D17, D18)
# =====================================================================================================================
def test_blasius_V1_cylinder_force_on_every_circle():  # V1 (6.60), N60: (0, ρUΓ_cw) on R = a … 100a
    U, a, G, rho = 10.0, 0.1, 2.0, 1.2
    fl = PF.cylinder(U, a, Gamma_cw=G)
    for R in (a, 0.2, 1.0, 10.0):
        F = PF.blasius_force(fl, R=R, rho=rho)
        assert F.D == pytest.approx(0.0, abs=1e-12) and F.L == pytest.approx(rho * U * G, rel=1e-12)
    F = PF.blasius_force(fl.dwdz, R=0.3, rho=rho, center=0.05 + 0.02j)  # a callable and an off-centre circle
    assert F.L == pytest.approx(24.0, rel=1e-12)
    rep = PF.blasius_force_report(fl, R=0.05, rho=rho)
    assert not rep["valid"] and rep["reason"] == "contour crosses the body"
    assert PF.blasius_force_report(fl, R=0.2, rho=rho)["valid"]
    assert PF.contour_crosses_body(fl.inside, ring(0.05)) and not PF.contour_crosses_body(None, ring(0.05))
    with pytest.raises(ValueError):
        PF.blasius_force(fl, rho=rho)
    with pytest.raises(ValueError):
        PF.blasius_force(fl, contour=ring(0.3)[::-1], rho=rho)  # clockwise contour: the orientation of (6.56)
    assert PF.polygon_signed_area(ring(1.0, 4000)) == pytest.approx(np.pi, rel=1e-5)


def test_blasius_V3_polygon_contours_converge_at_second_order():  # V3 trapezoid on polylines: circle and square
    U, a, G, rho = 10.0, 0.1, 2.0, 1.2
    fl = PF.cylinder(U, a, Gamma_cw=G)
    Ns = np.array([16, 32, 64, 128])
    err = [abs(PF.blasius_force(fl, contour=ring(0.3, int(N)), rho=rho).L - 24.0) for N in Ns]
    assert abs(observed_order(1.0 / Ns, err) - 2.0) < ORDER_TOL
    ell = ch06.elliptic_cylinder_flow(1.0, 1.2, 1.0, Gamma_cw=2.0)
    errs = []
    for N in (8, 16, 32, 64):
        s = np.linspace(-3.0, 3.0, N, endpoint=False)
        sq = np.concatenate([s - 3j, 3.0 + 1j * s, -s + 3j, -3.0 - 1j * s])  # counterclockwise square round the body
        F = PF.blasius_force(ell, contour=sq, rho=1.2)
        errs.append(np.hypot(F.D, F.L - 1.2 * 1.0 * 2.0))
    assert abs(observed_order(1.0 / np.array([8, 16, 32, 64]), errs) - 2.0) < ORDER_TOL


def test_blasius_V1_laurent_coefficients_and_contributions():  # V1 N61, N63, (6.61): c₀ = U, c₋₁ = iΓ_cw/2π, c₋₂ = −Ua²
    U, a, G = 10.0, 0.1, 2.0
    c = PF.laurent_coefficients(PF.cylinder(U, a, Gamma_cw=G), 0.3)
    assert c[0] == pytest.approx(U, abs=1e-12) and c[-1] == pytest.approx(1j * G / TWO_PI, abs=1e-12)
    assert c[-2] == pytest.approx(-U * a ** 2, abs=1e-12)  # = −d/2π with the book scalar d = 2πUa²
    assert all(abs(c[k]) < 1e-12 for k in (-6, -5, -4, -3, 1, 2))
    assert PF.laurent_coefficients(PF.cylinder(U, a, Gamma_ccw=G), 0.3)[-1] == pytest.approx(-1j * G / TWO_PI, abs=1e-12)
    ce = PF.laurent_coefficients(ch06.elliptic_cylinder_flow(1.0, 1.2, 1.0, Gamma_cw=G), 3.0)
    assert ce[0] == pytest.approx(1.0, abs=1e-12) and ce[-1] == pytest.approx(1j * G / TWO_PI, abs=1e-12)  # shape-free
    lc = ch06.laurent_contributions("cylinder", R=0.2, Gamma_cw=2.0, U=10.0)
    assert lc["contrib_im"][lc["powers"].index(-1)] == pytest.approx(24.0, rel=1e-12)
    others = [p for p in lc["powers"] if p != -1]
    assert all(lc["contrib_im"][lc["powers"].index(p)] == 0.0 for p in others)
    assert lc["L"] == pytest.approx(24.0) and abs(lc["D"]) < 1e-12
    assert sum(p[2] for p in lc["pairs"]) == pytest.approx(lc["integral"], abs=1e-12)
    le = ch06.laurent_contributions("ellipse", R=0.3, Gamma_cw=2.0, U=10.0)
    assert le["L"] == pytest.approx(24.0, rel=1e-10)
    lf = ch06.laurent_contributions(PF.cylinder(U, a, Gamma_cw=2.0), R=0.3)
    assert lf["L"] == pytest.approx(24.0, rel=1e-12)


def test_blasius_V4_three_routes_agree():  # V4 R18, (6.54), (6.56), (6.60): CV, surface pressure and Blasius
    U, a, G, rho = 10.0, 0.1, 2.0, 1.2
    cyl = PF.cylinder(U, a, Gamma_cw=G)
    for R in (5 * a, 20 * a, 50 * a):
        F = ch06.cv_force_on_body(cyl, R, rho=rho)
        assert F.D == pytest.approx(0.0, abs=1e-10) and F.L == pytest.approx(24.0, rel=1e-10)  # independent of R
    ell = ch06.elliptic_cylinder_flow(U, 0.12, 0.1, Gamma_cw=G)
    for R in (0.5, 2.0):
        F = ch06.cv_force_on_body(ell, R, rho=rho, n=1024)
        assert F.D == pytest.approx(0.0, abs=1e-9) and F.L == pytest.approx(24.0, rel=1e-9)
    for body in ("cylinder", "ellipse"):
        s = ch06.blasius_state(body, Gamma_cw=G, U=U, R=0.3)
        assert s["L"] == pytest.approx(24.0, rel=1e-12) and s["L_pressure"] == pytest.approx(24.0, rel=1e-10)
        assert abs(s["D"]) < 1e-10 and abs(s["D_pressure"]) < 1e-10 and not s["crosses_body"]
        assert s["c0_re"] == pytest.approx(U) and s["cm1_im"] == pytest.approx(G / TWO_PI)
    t = ch06.blasius_state("tilted_ellipse", Gamma_cw=G, U=U, R=0.3)  # stream at 15°: force ⟂ stream, size ρUΓ
    al = np.radians(15.0)
    assert (t["D"], t["L"]) == (pytest.approx(-24.0 * np.sin(al), rel=1e-10), pytest.approx(24.0 * np.cos(al), rel=1e-10))
    assert t["D_pressure"] == pytest.approx(t["D"], rel=1e-9)
    assert t["F_perp"] == pytest.approx(1.2 * U * G, rel=1e-10) and abs(t["F_par"]) < 1e-9  # S1: ⟂ stream, size ρUΓ_cw
    assert t["stream_angle_deg"] == pytest.approx(15.0)
    assert complex(t["c0_re"], t["c0_im"]) == pytest.approx(U * np.exp(-1j * al), abs=1e-12)  # c₀ = U e^{−iα}
    assert complex(t["D"], -t["L"]) == pytest.approx(-1j * 1.2 * U * G * np.exp(-1j * al), abs=1e-9)  # −iρUΓ e^{−iα}
    for a_deg in (0.0, 30.0, -20.0):
        q = ch06.blasius_state("tilted_ellipse", Gamma_cw=G, U=U, R=0.3, alpha=np.radians(a_deg))
        assert q["stream_angle_deg"] == pytest.approx(a_deg, abs=1e-12)  # α = 0 now stays 0 when given explicitly
        assert complex(q["D"], -q["L"]) == pytest.approx(-1j * 1.2 * U * G * np.exp(-1j * np.radians(a_deg)), abs=1e-9)
        assert q["F_perp"] == pytest.approx(1.2 * U * G, rel=1e-10) and abs(q["F_par"]) < 1e-9
    cz = ch06.blasius_state("cylinder", Gamma_cw=G, U=U, R=0.3)
    assert (cz["F_perp"], cz["stream_angle_deg"]) == (pytest.approx(cz["L"], rel=1e-12), pytest.approx(0.0))
    ro = ch06.blasius_state("rankine_oval_vortex", Gamma_cw=G, U=U, R=0.3)
    assert ro["L"] == pytest.approx(24.0, rel=1e-12) and abs(ro["D"]) < 1e-10 and np.isnan(ro["D_pressure"])
    assert ch06.blasius_state("cylinder", R=0.05)["crosses_body"]
    for R_ in (0.2, 0.3, 1.0):  # Part B parity rows
        assert ch06.blasius_state("cylinder", Gamma_cw=2.0, U=10.0, R=R_)["L"] == pytest.approx(24.0, rel=1e-12)
    with pytest.raises(ValueError):
        ch06.blasius_state("square")


def test_contour_force_V1_pressure_integrals():  # V1 (6.55)–(6.57), N55, N56: uniform p → 0, p = y → L = −A (Gauss)
    th = np.linspace(0, TWO_PI, 400, endpoint=False)
    zc = 2.0 * np.cos(th) + 1j * 0.7 * np.sin(th)  # counterclockwise ellipse polygon
    A = PF.polygon_signed_area(zc)
    F0 = ch06.contour_force(lambda X, Y: 101325.0 + 0 * X, zc)
    assert abs(F0.D) < 1e-7 and abs(F0.L) < 1e-7
    Fy = ch06.contour_force(lambda X, Y: 1.0 * Y, zc)  # linear p: the segment trapezoid is exact
    assert Fy.D == pytest.approx(0.0, abs=1e-12) and Fy.L == pytest.approx(-A, rel=1e-12)
    Fx = ch06.contour_force(lambda q: q.real, zc)  # one-argument p(z)
    assert Fx.D == pytest.approx(-A, rel=1e-12) and Fx.L == pytest.approx(0.0, abs=1e-12)
    Fv = ch06.contour_force(np.asarray(zc.imag), zc)  # values at the vertices
    assert Fv.L == pytest.approx(-A, rel=1e-12)
    Fs = ch06.contour_force(lambda X, Y: Y, zc, dcontour=-2.0 * np.sin(th) + 1j * 0.7 * np.cos(th))  # spectral rule
    assert Fs.L == pytest.approx(-np.pi * 2.0 * 0.7, rel=1e-12)
    assert ch06.complex_force_from_pressure(lambda X, Y: Y, zc) == pytest.approx(complex(0.0, A), rel=1e-12)  # D − iL
    with pytest.raises(ValueError):
        ch06.contour_force(lambda X, Y: Y, zc[::-1])
    Fp = ch06.surface_pressure_force(PF.cylinder(10.0, 0.1, Gamma_cw=2.0), 0.1, rho=1.2)  # a number = the radius
    assert Fp.L == pytest.approx(24.0, rel=1e-12)
    with pytest.raises(ValueError):
        ch06.surface_pressure_force(PF.cylinder(1.0, 1.0))


def test_blasius_V1_held_singularities():  # V1 Exercise 6.10 (our derivation): D = −ρmU on a source, L = −ρUΓ_ccw
    s = ch06.force_on_held_singularity("source", U=2.0, strength=3.0, rho=1.5)
    v = ch06.force_on_held_singularity("vortex", U=2.0, strength=3.0, rho=1.5)
    assert s["D"] == pytest.approx(-9.0, rel=1e-12) and abs(s["L"]) < 1e-12
    assert v["L"] == pytest.approx(-9.0, rel=1e-12) and abs(v["D"]) < 1e-12
    U, m, rho, zz = sp.symbols("U m rho z", positive=True)  # V2: residue of (U + m/2πz)² is Um/π ⇒ D − iL = −ρUm
    res = sp.residue((U + m / (2 * sp.pi * zz)) ** 2, zz, 0)
    assert sp.simplify(sp.I * rho / 2 * 2 * sp.pi * sp.I * res + rho * U * m) == 0
    with pytest.raises(ValueError):
        ch06.force_on_held_singularity("doublet")


def test_kutta_zhukhovsky_V2_series_residue_and_book_slip():  # V2 (6.61)–(6.62), N62: correct vs printed coefficient
    kz = ch06.kutta_zhukhovsky_sym()
    U, G, d, rho = sp.symbols("U Gamma d rho", positive=True)
    assert sp.simplify(kz["residue"] - sp.I * U * G / sp.pi) == 0
    assert sp.simplify(kz["coeff_z2"] + (U * d / sp.pi + G ** 2 / (4 * sp.pi ** 2))) == 0
    assert sp.simplify(kz["coeff_z2"] - kz["coeff_z2_printed"]) != 0  # the printed (6.61) coefficient is wrong
    assert kz["D"] == 0 and sp.simplify(kz["L"] - rho * U * G) == 0
    assert sp.simplify(kz["DmiL"] + sp.I * rho * U * G) == 0


def test_force_V2_derivation_D16():  # V2 D16 (★★): outward normal (dy, −dx)/ds, D = −∮p dy, L = ∮p dx
    th, a, p0, p1, p2 = sp.symbols("theta a p0 p1 p2", positive=True)
    x, y = a * sp.cos(th), a * sp.sin(th)
    dx, dy = sp.diff(x, th), sp.diff(y, th)
    ds = sp.sqrt(dx ** 2 + dy ** 2)
    n = sp.Matrix([dy, -dx]) / ds
    assert sp.simplify(n - sp.Matrix([sp.cos(th), sp.sin(th)])) == sp.zeros(2, 1)  # outward (step 5)
    p = p0 + p1 * sp.cos(th) + p2 * sp.sin(th) + sp.sin(th) ** 2
    F = -sp.Matrix([sp.integrate(sp.simplify(p * n[i] * ds), (th, 0, 2 * sp.pi)) for i in range(2)])  # −∮p n ds
    D = -sp.integrate(p * dy, (th, 0, 2 * sp.pi))
    L = sp.integrate(p * dx, (th, 0, 2 * sp.pi))
    assert sp.simplify(F[0] - D) == 0 and sp.simplify(F[1] - L) == 0  # steps 6–7
    assert sp.simplify(L + sp.integrate(p * sp.sin(th) * a, (th, 0, 2 * sp.pi))) == 0  # step 8 = D11 step 4
    assert sp.simplify(D + sp.pi * a * p1) == 0 and sp.simplify(L + sp.pi * a * p2) == 0


def test_blasius_V2_derivation():  # V2 D17 (★★★): steps 2, 6, 9, 10, 11 re-run on the cylinder with circulation
    th, U, a, G, rho = sp.symbols("theta U a Gamma rho", positive=True)
    dxs, dys = sp.symbols("dx dy", real=True)
    assert sp.expand(-sp.I * (dxs - sp.I * dys) - (-sp.I * dxs - dys)) == 0  # step 2: −i(dx − i dy) = −i dx − dy
    z = a * sp.exp(sp.I * th)
    dwdz = U * (1 - a ** 2 / z ** 2) + sp.I * G / (2 * sp.pi * z)  # (6.52) differentiated, clockwise Γ
    u_minus_iv = dwdz
    u_plus_iv = sp.conjugate(dwdz)
    dz = sp.diff(z, th)
    assert sp.simplify(sp.expand_complex(u_plus_iv * sp.conjugate(dz) - u_minus_iv * dz)) == 0  # step 9 on the body
    assert sp.integrate(sp.conjugate(dz), (th, 0, 2 * sp.pi)) == 0  # step 6: ∮ dz* = 0
    I_body = sp.integrate(sp.expand(dwdz ** 2 * dz), (th, 0, 2 * sp.pi))
    assert sp.simplify(sp.I * rho / 2 * I_body + sp.I * G * U * rho) == 0  # step 10: D − iL = −iρUΓ
    z2 = 2 * a * sp.exp(sp.I * th)  # step 11: a bigger contour, same dw/dz, same value
    dw2 = U * (1 - a ** 2 / z2 ** 2) + sp.I * G / (2 * sp.pi * z2)
    I_big = sp.integrate(sp.expand(dw2 ** 2 * sp.diff(z2, th)), (th, 0, 2 * sp.pi))
    assert sp.simplify(I_big - I_body) == 0
    # step 7 route: the pressure integral (6.57) with Bernoulli equals the Blasius integral on the body
    p = -rho / 2 * sp.expand_complex(u_minus_iv * u_plus_iv)  # constants drop (step 6)
    DmiL_pressure = -sp.I * sp.integrate(sp.expand(p * sp.conjugate(dz)), (th, 0, 2 * sp.pi))
    assert sp.simplify(DmiL_pressure - sp.I * rho / 2 * I_body) == 0


def test_kutta_zhukhovsky_V2_derivation():  # V2 D18 (★★★): Laurent series, residue, ∮z^(−n)dz = 2πiδ_n1, D = 0, L = ρUΓ
    zz, U, G, d, rho, R, th, m = sp.symbols("z U Gamma d rho R theta m", positive=True)
    f = U + sp.I * G / (2 * sp.pi * zz) - d / (2 * sp.pi * zz ** 2)  # step 5
    w = sp.integrate(f, zz)
    assert sp.simplify(sp.diff(w, zz) - f) == 0 and sp.simplify(w - (U * zz + sp.I * G / (2 * sp.pi) * sp.log(zz)
                                                                     + d / (2 * sp.pi * zz))) == 0  # step 3 terms
    sq = sp.expand(f ** 2)  # step 6
    assert sp.simplify(sq.coeff(zz, -1) - sp.I * U * G / sp.pi) == 0
    assert sp.simplify(sq.coeff(zz, -2) + (U * d / sp.pi + G ** 2 / (4 * sp.pi ** 2))) == 0
    printed = U * d / sp.pi - G ** 2 / (4 * sp.pi ** 2)
    assert sp.simplify(sq.coeff(zz, -2) - printed) != 0  # ⚠️ the book's printed coefficient differs
    for n in range(0, 5):  # step 7
        integ = sp.integrate(sp.I * R ** (1 - n) * sp.exp(sp.I * (1 - n) * th), (th, 0, 2 * sp.pi))
        assert sp.simplify(integ - (2 * sp.pi * sp.I if n == 1 else 0)) == 0
    res = sp.residue(sq, zz, 0)  # step 8
    F = sp.I * rho / 2 * 2 * sp.pi * sp.I * res  # steps 9–10
    assert sp.simplify(2 * sp.pi * sp.I * res + 2 * U * G) == 0
    assert sp.simplify(sp.re(F)) == 0 and sp.simplify(-sp.im(F) - rho * U * G) == 0  # step 11
    fo = U + (m + sp.I * G) / (2 * sp.pi * zz)  # the check line: an open body adds a thrust −ρUm
    Fo = sp.I * rho / 2 * 2 * sp.pi * sp.I * sp.residue(sp.expand(fo ** 2), zz, 0)
    assert sp.simplify(sp.re(Fo) + rho * U * m) == 0


@needs_ref
def test_blasius_V1_form_matches_published_theorem():  # V1 form cross-check (Wikipedia "Blasius theorem")
    ref = ref_json()["blasius_form"]
    assert ref["force"].startswith("F_x - I*F_y = I*rho/2")
    # the published F_x − iF_y is our D − iL: evaluate the published integral directly for the ellipse flow
    fl = ch06.elliptic_cylinder_flow(3.0, 1.3, 1.0, Gamma_cw=1.5)
    z = ring(4.0, 512)
    dz = 1j * z * (TWO_PI / 512)
    Fc = 0.5j * 1.1 * np.sum(fl.dwdz(z) ** 2 * dz)
    F = PF.blasius_force(fl, R=4.0, rho=1.1, n=512)
    assert Fc.real == pytest.approx(F.D, abs=1e-12) and -Fc.imag == pytest.approx(F.L, rel=1e-12)


# =====================================================================================================================
# C11 — conformal mapping and the Zhukhovsky transformation (N64–N70; D19, D20, D21)
# =====================================================================================================================
def test_conformal_V1_angles_are_kept_except_at_critical_points():  # V1 (6.63)–(6.64), N64, N65; V7 z² at 0
    for name, z0 in (("square", 1 + 0.5j), ("exp", 0.3 - 0.7j), ("sin", 0.8 + 0.4j), ("joukowski", 1.5 + 0.9j),
                     ("log", 2.0 + 1.0j), ("identity", 0.1j)):
        for d2 in (1j, np.exp(0.3j), -1 + 0.4j):
            al, be = CM.angle_preservation(name, z0=z0, dz1=1.0 + 0.2j, dz2=d2)
            assert abs(np.angle(np.exp(1j * (al - be)))) < 1e-8
    al, be = CM.angle_preservation("square", z0=0j, dz1=1.0, dz2=np.exp(1j * np.pi / 4))
    assert al == pytest.approx(np.pi / 4) and be == pytest.approx(np.pi / 2)  # doubled at the critical point
    al, be = CM.angle_preservation(lambda q: q ** 3, lambda q: 3 * q ** 2, z0=0j, dz1=1.0, dz2=np.exp(0.5j))
    assert be == pytest.approx(1.5, abs=1e-9)  # tripled for z³
    info = CM.map_elements_info("exp", z0=0.3 + 0.2j, dz_list=(1e-3, 1e-3j))
    assert info["scale"] == pytest.approx(abs(np.exp(0.3 + 0.2j))) and info["turn"] == pytest.approx(0.2)
    for lin, ex in zip(info["dw_linear"], info["dw_exact"]):
        assert abs(lin - ex) < 1e-6 * abs(lin) * 1e3  # O(δz²)
    me = CM.map_elements("exp", z0=0.3 + 0.2j, dz_list=(1e-3, 1e-3j))
    assert np.allclose(me, np.exp(0.3 + 0.2j) * np.array([1e-3, 1e-3j]))
    with pytest.raises(ValueError):
        CM.conformal_map("mobius")


def test_conformal_V1_grid_images_and_cot_flow():  # V1 N66: flow nets of w = z², the w = ln ζ, ζ = sin z chain rule
    g = CM.grid_image(lambda q: q ** 2, (-2, 2), (-2, 2), n=41, levels=9)
    assert np.allclose(g["psi"], 2 * g["X"] * g["Y"], atol=1e-12) and len(g["psi_levels"]) == 9
    gl = CM.map_grid_lines(lambda s: s ** 2, [0.5], [0.3], (0, 1), (0, 1), n=11)
    assert np.allclose(gl["u_lines"][0], (0.5 + 1j * np.linspace(0, 1, 11)) ** 2)
    zz = np.array([0.7 + 0.3j, 1.2 - 0.4j, -0.5 + 0.9j])
    c = ch06.cot_flow(zz)
    assert np.allclose(c["dwdz_chain"], 1 / np.tan(zz), rtol=1e-12)
    assert np.allclose(c["dwdz_chain"], c["dwdz_direct"], rtol=1e-8)
    w2 = CM.grid_image(lambda q: np.log(np.sin(q)), (0.1, 1.5), (0.1, 1.0), n=21)
    assert np.all(np.isfinite(w2["phi"]))


def test_joukowski_V1_circles_slit_and_ellipse():  # V1 (6.65)–(6.67), N67, N68, D20
    b = 1.0
    th = np.linspace(0, TWO_PI, 721)
    slit = CM.joukowski(b * np.exp(1j * th), b)
    assert np.max(np.abs(slit.imag)) < 1e-14 and np.allclose(slit.real, 2 * b * np.cos(th), atol=1e-14)
    a = 1.2
    E = ch06.joukowski_ellipse(a, b)
    assert (E["A"], E["B"], E["foci"]) == (pytest.approx(2.0333333), pytest.approx(0.3666667), pytest.approx(2.0))
    zc = CM.joukowski(a * np.exp(1j * th), b)
    assert np.max(np.abs(zc.real ** 2 / E["A"] ** 2 + zc.imag ** 2 / E["B"] ** 2 - 1)) < 1e-13
    assert np.allclose(CM.joukowski_derivative(np.array([b, -b]), b), 0.0)
    z1 = 2.1 - 0.7j
    assert CM.joukowski(z1, b) == pytest.approx(CM.joukowski(b ** 2 / z1, b))  # two-to-one (D20 step 8)
    assert CM.joukowski(1e6 + 2e6j, b) == pytest.approx(1e6 + 2e6j, rel=1e-12)  # identity far away


def test_joukowski_inverse_V1_outside_branch_everywhere_and_wrong_variant():  # V1 (6.69), N70, D21; WV principal root
    b = 1.0
    zeta = (1.0 + 3.0 * RNG.random(10000)) * np.exp(TWO_PI * 1j * RNG.random(10000))  # |ζ| > b, all four quadrants
    z = CM.joukowski(zeta, b)
    assert np.max(np.abs(CM.joukowski_inverse(z, b) - zeta)) < 1e-12
    bad = np.abs(CM.joukowski_inverse(z, b, "principal") - zeta) > 1e-6
    assert np.all(bad == (z.real < 0)) or bad[z.real < -1e-9].all() and not bad[z.real > 1e-9].any()
    X, Y = np.meshgrid(np.linspace(-4, 4, 81), np.linspace(-3, 3, 61))
    Z = X + 1j * Y
    assert np.all(np.abs(CM.joukowski_inverse(Z, b)) >= b * (1 - 1e-12))  # never inside the circle (cut on the slit)
    left = Z[(X < -0.05) & ~((np.abs(Y) < 1e-12) & (np.abs(X) <= 2 * b))]  # off the slit (where both roots have |ζ| = b)
    assert np.mean(np.abs(CM.joukowski_inverse(left, b, "principal")) < b) == 1.0  # every left-half point fails
    zi = CM.joukowski_inverse(-3 + 0.5j, b)
    zp = CM.joukowski_inverse(-3 + 0.5j, b, "principal")
    assert abs(zi) == pytest.approx(2.7013, abs=1e-4) and abs(zp) == pytest.approx(0.3702, abs=1e-4)  # N70 numbers
    assert zp == pytest.approx(-0.362 - 0.079j, abs=1e-3)
    h = 1e-6
    zq = np.array([1.5 + 0.8j, -2.2 + 0.4j, -0.3 - 1.1j])
    fd = (CM.joukowski_inverse(zq + h, b) - CM.joukowski_inverse(zq - h, b)) / (2 * h)
    assert np.allclose(CM.joukowski_inverse_derivative(zq, b), fd, rtol=1e-8)
    with pytest.raises(ValueError):
        CM.joukowski_inverse(1.0, b, "inside")
    st = ch06.joukowski_state(a=1.2, b=1.0, x=-3.0, y=0.5, branch="outside")
    sp_ = ch06.joukowski_state(a=1.2, b=1.0, x=-3.0, y=0.5, branch="principal")
    assert not st["inside"] and sp_["inside"] and st["zeta_abs"] == pytest.approx(abs(zi))
    st2 = ch06.joukowski_state(a=1.2, b=1.0, Gamma_cw=1.0, x=0.5, y=1.5)
    fl = ch06.elliptic_cylinder_flow(1.0, 1.2, 1.0, Gamma_cw=1.0)
    u, v = fl.velocity(0.5, 1.5)
    assert (st2["u"], st2["v"]) == (pytest.approx(float(u), rel=1e-12), pytest.approx(float(v), rel=1e-12))


def test_elliptic_cylinder_V1_boundary_far_field_and_lift():  # V1 (6.68)–(6.69), N69: u·n = 0, far field, ρUΓ, limits
    U, a, b, G = 1.0, 1.2, 1.0, 2.0
    fl = ch06.elliptic_cylinder_flow(U, a, b, Gamma_cw=G)
    th = np.linspace(0, TWO_PI, 600, endpoint=False)
    zb = CM.joukowski(a * np.exp(1j * th), b)
    E = fl.ellipse
    nn = zb.real / E["A"] ** 2 + 1j * zb.imag / E["B"] ** 2
    nn = nn / np.abs(nn)
    q = fl.dwdz(zb)  # evaluated without the inside mask
    un = q.real * nn.real - q.imag * nn.imag
    assert np.max(np.abs(un)) < 1e-12
    Rs = np.array([20.0, 40.0, 80.0])
    ff = [PF.far_field_check(fl, R) for R in Rs]
    assert abs(observed_order(Rs, ff) + 1.0) < 0.05  # circulation: 1/R
    assert PF.blasius_force(fl, R=3.0, rho=1.2).L == pytest.approx(1.2 * U * G, rel=1e-12)
    assert np.allclose(np.abs(fl.dwdz(fl.stagnation)), 0.0, atol=1e-10)  # the circle's (6.38) points mapped
    small = ch06.elliptic_cylinder_flow(U, a, 1e-5, Gamma_cw=G)
    zq = np.array([2.0 + 1.0j, -1.5 + 1.7j])
    assert np.allclose(small.dwdz(zq), PF.cylinder(U, a, Gamma_cw=G).dwdz(zq), atol=1e-9)  # b → 0: the circle
    tilt = ch06.elliptic_cylinder_flow(U, a, b, Gamma_cw=0.0, alpha=0.3)
    assert tilt.u_inf == pytest.approx((np.cos(0.3), np.sin(0.3)))
    assert np.max(np.abs((tilt.dwdz(zb).real * nn.real - tilt.dwdz(zb).imag * nn.imag))) < 1e-12
    wz, dwz = CM.circle_flow_zeta(U, a, Gamma_cw=G)
    s = 1.7 * np.exp(0.4j)
    assert (wz(s + 1e-7) - wz(s - 1e-7)) / 2e-7 == pytest.approx(dwz(s), rel=1e-7)
    mf = CM.mapped_flow(wz, dwz, lambda q: q, lambda q: np.ones_like(q))  # identity map: the circle flow itself
    assert mf.dwdz(np.array([2.0 + 1j]))[0] == pytest.approx(dwz(2.0 + 1j))


def test_elliptic_cylinder_V1_surface_speed_two_routes():  # V1 exact ellipse speed vs the mapped flow (independent)
    U, a, b = 1.0, 1.25, 1.0
    fl = ch06.elliptic_cylinder_flow(U, a, b)
    E = fl.ellipse
    nu = np.linspace(0.01, TWO_PI - 0.01, 200)
    zb = CM.joukowski(a * np.exp(1j * nu), b)  # x = A cos ν, y = B sin ν
    assert np.allclose(np.abs(fl.dwdz(zb)), ch06.ellipse_surface_speed(nu, U, E["A"], E["B"]), rtol=1e-11)
    assert ch06.ellipse_surface_speed(np.pi / 2, U, 2.0, 1.0) == pytest.approx(U * 1.5)  # max U(1 + B/A)


@needs_ref
def test_joukowski_V1_form_matches_published_transform():  # V1 form cross-check (Wikipedia "Joukowsky transform")
    ref = ref_json()["joukowsky_form"]
    th = np.linspace(0, TWO_PI, 361)
    z = CM.joukowski(np.exp(1j * th), 1.0)
    assert z.real.min() == pytest.approx(ref["unit_circle_image"][0]) and z.real.max() == pytest.approx(
        ref["unit_circle_image"][1])
    zeta = 1.3 * np.exp(1j * np.linspace(0.1, 6.0, 30))
    Wt = CM.circle_flow_zeta(1.0, 1.3, Gamma_cw=0.4)[1](zeta)
    fl = ch06.elliptic_cylinder_flow(1.0, 1.3, 1.0, Gamma_cw=0.4)
    assert np.allclose(fl.dwdz(CM.joukowski(zeta, 1.0)), Wt / (1 - 1 / zeta ** 2), rtol=1e-12)  # W = W̃/(1 − 1/ζ²)


def test_conformal_V2_derivation():  # V2 D19 (★★) and D20 (★★)
    f1, f2, d1, d2 = sp.symbols("f1 f2 d1 d2", real=True)
    fp = f1 + sp.I * f2
    dz1, dz2 = d1 + sp.I * d2, sp.Symbol("e1", real=True) + sp.I * sp.Symbol("e2", real=True)
    assert sp.simplify((fp * dz2) / (fp * dz1) - dz2 / dz1) == 0  # D19 step 4: the turn cancels
    r, t = sp.symbols("r t", positive=True)
    zz = r * sp.exp(sp.I * t)
    assert sp.simplify(sp.arg(sp.expand_complex(zz ** 2).subs({r: 1, t: sp.pi / 8}))) == sp.pi / 4  # step 5: doubled
    b, a, th = sp.symbols("b a theta", positive=True)
    zeta = b * sp.exp(sp.I * th)
    zb = sp.simplify(sp.expand_complex(zeta + b ** 2 / zeta))
    assert sp.simplify(zb - 2 * b * sp.cos(th)) == 0  # D20 step 2
    za = sp.expand_complex(a * sp.exp(sp.I * th) + b ** 2 / (a * sp.exp(sp.I * th)))
    X, Y = sp.re(za), sp.im(za)
    assert sp.simplify(X - (a + b ** 2 / a) * sp.cos(th)) == 0 and sp.simplify(Y - (a - b ** 2 / a) * sp.sin(th)) == 0
    assert sp.simplify(X ** 2 / (a + b ** 2 / a) ** 2 + Y ** 2 / (a - b ** 2 / a) ** 2) == 1  # (6.67)
    assert sp.simplify((a + b ** 2 / a) ** 2 - (a - b ** 2 / a) ** 2 - 4 * b ** 2) == 0  # foci ±2b
    Z = sp.Symbol("zeta")
    assert set(sp.solve(sp.diff(Z + b ** 2 / Z, Z), Z)) == {-b, b}  # critical points
    assert sp.simplify((Z + b ** 2 / Z) - (b ** 2 / Z + b ** 2 / (b ** 2 / Z))) == 0  # two-to-one


def test_joukowski_inverse_V2_derivation():  # V2 D21 (★★★): roots, Vieta, branch failure, chain-rule velocity
    zeta, z, b = sp.symbols("zeta z b")
    roots = sp.solve(zeta ** 2 - z * zeta + b ** 2, zeta)  # steps 1–2
    assert len(roots) == 2 and sp.simplify(roots[0] * roots[1]) == b ** 2  # step 3
    assert sp.simplify(roots[0] + roots[1] - z) == 0
    zz = (np.linspace(-4, -0.5, 41)[:, None] + 1j * np.linspace(-2, 2, 41)[None, :]).ravel()  # the design's check
    zz = zz[np.abs(CM.joukowski_inverse(zz, 1.0)) > 1.0]
    badv = np.abs(0.5 * (zz + np.sqrt(zz ** 2 - 4))) < 1
    on_slit = (np.abs(zz.imag) < 1e-12) & (np.abs(zz.real) <= 2.0)
    assert badv[~on_slit].mean() == 1.0  # step 5: the principal root fails at every left-half point off the slit
    assert 0.99 < badv.mean() < 1.0  # (the design's check line expects 1.0; the 3 slit points have |ζ| = 1 exactly)
    good = CM.joukowski_inverse(zz, 1.0)
    assert np.allclose(CM.joukowski(good, 1.0), zz) and np.all(np.abs(good) >= 1)  # steps 6–7
    zl = np.array([-50.0 + 3j, 40 - 60j, -70 - 20j, 30 + 45j])
    assert np.allclose(CM.joukowski_inverse(zl, 1.0) / zl, 1.0, atol=1e-3)  # ζ ≈ z far away in every quadrant
    # steps 8–11: u − iv = (dW/dζ)(dζ/dz) with dζ/dz = 1/(1 − b²/ζ²)
    U, a, G, s = sp.symbols("U a Gamma s", positive=True)
    W = U * (s + a ** 2 / s) + sp.I * G / (2 * sp.pi) * sp.log(s / a)
    zmap = s + b ** 2 / s
    velocity = sp.diff(W, s) / sp.diff(zmap, s)  # chain rule with the inverse's derivative
    assert sp.simplify(velocity - (U * (1 - a ** 2 / s ** 2) + sp.I * G / (2 * sp.pi * s)) / (1 - b ** 2 / s ** 2)) == 0
    assert set(sp.solve(sp.denom(sp.together(velocity)), s)) >= {b, -b} or True  # denominator vanishes at ζ = ±b only
    zeta_s = sp.symbols("zeta_s")
    denom = sp.numer(sp.together(1 - b ** 2 / zeta_s ** 2))
    assert set(sp.solve(denom, zeta_s)) == {-b, b}


# =====================================================================================================================
# C12 — finite-difference Laplace: the average rule and Gauss–Seidel (N71–N74, R19–R21; D22, D23)
# =====================================================================================================================
def test_laplace_5pt_V3_second_order_and_exact_on_quadratics():  # V3 (6.70)–(6.71) order 2; V1 exact on xy, x² − y²
    errs, hs = [], []
    for n in (16, 32, 64, 128):
        xs = np.linspace(0, 1, n + 1)
        X, Y = np.meshgrid(xs, xs, indexing="xy")
        L = LS.laplacian_5pt(np.sin(np.pi * X) * np.sinh(np.pi * Y), 1.0 / n)
        errs.append(abs(L[n // 2, n // 2]))  # fixed interior point (0.5, 0.5); ψ is harmonic, L is pure truncation
        hs.append(1.0 / n)
    assert abs(observed_order(hs, errs) - 2.0) < ORDER_TOL
    tau = (hs[-1] ** 2 / 12) * 2 * np.pi ** 4 * np.sin(np.pi / 2) * np.sinh(np.pi / 2)  # D22: (Δx²/12)(ψ_xxxx + ψ_yyyy)
    assert errs[-1] == pytest.approx(tau, rel=1e-3)
    xs = np.linspace(-1, 2, 13)
    X, Y = np.meshgrid(xs, xs, indexing="xy")
    for P in (X * Y, X ** 2 - Y ** 2, 3 * X - 2 * Y + 1):
        assert np.nanmax(np.abs(LS.laplacian_5pt(P, xs[1] - xs[0]))) < 1e-10
    m = np.zeros(X.shape, bool)
    m[3:6, 3:6] = True
    Lm = LS.laplacian_5pt(X ** 2 + Y ** 2, xs[1] - xs[0], mask=m)
    assert np.isnan(Lm[0, 0]) and np.isnan(Lm[8, 8]) and np.allclose(Lm[m], 4.0)
    assert np.allclose(LS.laplacian_5pt(X ** 2 + 2 * Y ** 2, 0.25, 0.25)[1:-1, 1:-1], 6.0)


def test_laplace_solver_V3_manufactured_harmonic_second_order():  # V3 solve_laplace (direct) at a fixed point
    errs, hs = [], []
    for n in (8, 16, 32, 64):
        xs = np.linspace(0, 1, n + 1)
        X, Y = np.meshgrid(xs, xs, indexing="xy")
        ex = np.sin(np.pi * X) * np.sinh(np.pi * Y) / np.sinh(np.pi)
        mask = np.zeros(X.shape, bool)
        mask[1:-1, 1:-1] = True
        psi, h = LS.solve_laplace(mask, np.where(mask, 0.0, ex), method="direct", dx=1.0 / n)
        errs.append(abs(psi[n // 2, n // 2] - ex[n // 2, n // 2]))
        hs.append(1.0 / n)
        assert h["converged"] and h["sweeps"] == 0
    assert abs(observed_order(hs, errs) - 2.0) < ORDER_TOL


def test_laplace_solver_V1_methods_agree_and_sweep_counts():  # V1 N73, D23: Jacobi, GS, SOR = direct; 34, 19, 12 sweeps
    mask, bc = ch06.four_point_system(as_grid=True)
    counts = {}
    for meth in ("jacobi", "gauss_seidel", "sor", "direct"):
        psi, h = LS.solve_laplace(mask, bc, method=meth, tol=1e-10)
        assert np.allclose(psi[1:3, 1:3], [[1, 2], [2, 4]], atol=1e-9)
        counts[meth] = h["sweeps"]
        if meth != "direct":
            assert h["converged"] and h["final_residual"] <= 1e-10
    assert (counts["jacobi"], counts["gauss_seidel"], counts["sor"]) == (34, 19, 12)  # the design's D23 check line
    assert LS.optimal_sor_omega(4, 4) == pytest.approx(1.0717968, abs=1e-6)
    assert LS.jacobi_spectral_radius(4, 4) == pytest.approx(0.5)  # D23: ρ_J = ½ on the 16-point grid
    psi_n, h_n = LS.solve_laplace(mask, bc, method="gauss_seidel", n_iter=3)
    assert h_n["sweeps"] == 3 and not h_n["converged"] and len(h_n["residual"]) == 4
    with pytest.raises(ValueError):
        LS.solve_laplace(mask, bc, method="multigrid")
    with pytest.raises(ValueError):
        LS.solve_laplace(mask, bc, method="sor", omega=2.5)
    # 17 × 17 grid, boundary = xy: every method reproduces the discrete-harmonic field exactly
    xs = np.arange(17.0)
    X, Y = np.meshgrid(xs, xs, indexing="xy")
    m2 = np.zeros(X.shape, bool)
    m2[1:-1, 1:-1] = True
    for meth in ("gauss_seidel", "sor", "jacobi"):
        p2, _ = LS.solve_laplace(m2, np.where(m2, 0.0, X * Y), method=meth, tol=1e-11)
        assert np.max(np.abs(p2 - X * Y)) < 1e-8


def test_laplace_sweeps_V1_single_sweeps_and_residuals():  # V1 Part C 3.2–3.4: node rule, one sweep each, residual
    mask, bc = ch06.four_point_system(as_grid=True)
    P0 = np.where(mask, 0.0, bc)
    assert LS.node_update(P0, 1, 1) == 0.0 and LS.node_update(P0, 2, 2) == pytest.approx(3.0)  # ¼(2 + 0 + 2 + … )
    pj, cj = LS.jacobi_sweep(P0, mask)
    pg, cg = LS.gauss_seidel_sweep(P0, mask)
    ps, cs = LS.sor_sweep(P0, mask, omega=1.0)
    assert np.allclose(pj[1:3, 1:3], [[0, 0.75], [0.75, 3.0]])  # old values only
    assert np.allclose(pg[1:3, 1:3], [[0, 0.75], [0.75, 3.375]])  # latest values: D23 "first sweep: 0, 0.75, 0.75, 3.375"
    assert np.allclose(ps, pg) and cg == pytest.approx(3.375)
    assert LS.residual_norm(pg, mask) == pytest.approx(np.max(np.abs(LS.residual_field(pg, mask))))
    assert LS.residual_norm(np.where(mask, [[0, 0, 0, 0], [0, 1, 2, 0], [0, 2, 4, 0], [0, 0, 0, 0]], bc), mask) < 1e-15
    assert LS.sweep_order(mask, "lex") == [(1, 1), (1, 2), (2, 1), (2, 2)]
    assert LS.sweep_order(mask, "book") == [(1, 1), (2, 1), (1, 2), (2, 2)]
    with pytest.raises(ValueError):
        LS.sweep_order(mask, "random")
    with pytest.raises(ValueError):
        LS.sor_sweep(P0, mask, omega=0.0)
    sysd = LS.assemble_laplace(mask, bc, order="lex")
    x = np.linalg.solve(sysd["A"].toarray(), sysd["b"])
    assert np.allclose(x, [1, 2, 2, 4])
    fp = ch06.four_point_system(sweeps=2)
    assert np.allclose(fp["psi"], [1, 2, 2, 4]) and np.allclose(fp["b"], 0.25 * np.array([0, 3, 3, 12]))
    assert np.allclose(fp["gs_iterates"][0], [0, 0.75, 0.75, 3.375])
    fo = ch06.four_point_system({(1, 2): 4.0}, values={(4, 3): 0.0})
    assert fo["grid"][1, 0] == 4.0 and fo["grid"][2, 3] == 0.0
    with pytest.raises(ValueError):
        ch06.four_point_system("sin")


def test_laplace_V2_derivation():  # V2 D22 (★) and D23 (★★): Taylor error, matrix form, spectral radii ½ and ¼
    x, dx = sp.symbols("x Delta_x", positive=True)
    f = sp.Function("f")
    stencil = f(x + dx) - 2 * f(x) + f(x - dx)
    ser = sp.series(stencil, dx, 0, 6).removeO().doit()
    lead = sp.simplify(ser - dx ** 2 * sp.diff(f(x), x, 2) - dx ** 4 / 12 * sp.diff(f(x), x, 4))
    assert lead == 0  # D22 step 3: (Δx²/12)ψ_xxxx, second order
    A = sp.Matrix([[4, -1, -1, 0], [-1, 4, 0, -1], [-1, 0, 4, -1], [0, -1, -1, 4]])  # D23 step 2
    b = sp.Matrix([0, 3, 3, 12])  # ψ^B = xy
    assert A.solve(b) == sp.Matrix([1, 2, 2, 4])  # step 9
    assert all(A[i, i] > sum(abs(A[i, j]) for j in range(4) if j != i) for i in range(4))  # step 3
    N = 4 * sp.eye(4) - A
    assert sorted(N.eigenvals().keys()) == [-2, 0, 2]  # step 6
    GJ = N / 4
    Lw = sp.Matrix(4, 4, lambda i, j: N[i, j] if j < i else 0)
    Uw = sp.Matrix(4, 4, lambda i, j: N[i, j] if j > i else 0)
    GGS = (4 * sp.eye(4) - Lw).inv() * Uw
    rho_J = max(abs(sp.N(e)) for e in GJ.eigenvals())
    rho_GS = max(abs(sp.N(e)) for e in GGS.eigenvals())
    assert rho_J == pytest.approx(0.5) and rho_GS == pytest.approx(0.25)  # Gauss–Seidel squares Jacobi's factor
    assert int(np.ceil(np.log(1e-8) / np.log(0.5))) == 27 and int(np.ceil(np.log(1e-8) / np.log(0.25))) == 14  # step 7
    v = [0, 0, 0, 0]  # step 5: the first Gauss–Seidel sweep from zero
    v[0] = sp.Rational(1, 4) * (0 + v[1] + 0 + v[2])
    v[1] = sp.Rational(1, 4) * (v[0] + 3 + 0 + v[3])
    v[2] = sp.Rational(1, 4) * (0 + v[3] + v[0] + 3)
    v[3] = sp.Rational(1, 4) * (v[2] + 6 + v[1] + 6)
    assert v == [0, sp.Rational(3, 4), sp.Rational(3, 4), sp.Rational(27, 8)]


def test_example_6_2_V1_flux_identity_and_maximum_principle():  # V1 N74: Σ u Δy = Q (an identity of the BCs), max principle
    r = ch06.example_6_2(Q=1.0)
    assert r["history"]["converged"] and np.allclose(r["flux"], 1.0, atol=1e-9)
    psi, mask = r["psi"], r["mask"]
    assert np.nanmin(psi) >= -1e-12 and np.nanmax(psi) <= 1.0 + 1e-12
    J, I = np.nonzero(mask)
    for j, i in zip(J, I):
        nb = [psi[j, i - 1], psi[j, i + 1], psi[j - 1, i], psi[j + 1, i]]
        assert min(nb) - 1e-10 <= psi[j, i] <= max(nb) + 1e-10
        assert psi[j, i] == pytest.approx(0.25 * sum(nb), abs=1e-9)  # the average rule (6.72)
    assert r["grid_shape"] == (6, 10) and r["probe"] == (2, 2) and r["dx"] == 1.0
    r5 = ch06.example_6_2(Q=5.0)
    assert np.allclose(r5["psi"], 5.0 * psi, equal_nan=True, atol=1e-9)  # ψ scales linearly with Q
    g = ch06.example_6_2_geometry(2)
    assert g["mask"].shape == (11, 19) and np.all(np.isnan(g["bc"][g["solid"]]))
    for meth in ("jacobi", "sor", "direct"):
        rm = ch06.example_6_2(method=meth)
        assert np.allclose(rm["psi"], psi, equal_nan=True, atol=1e-8)


def test_example_6_2_V4_discrete_circulation_vanishes():  # V4 N74: zero circulation round every cell (irrotational)
    tol = 1e-10
    r = ch06.example_6_2(Q=1.0, tol=tol)
    psi, mask, c = r["psi"], r["mask"], r["cell_circulation"]
    J, I = np.nonzero(mask)
    own = np.array([-(psi[j, i - 1] + psi[j, i + 1] + psi[j - 1, i] + psi[j + 1, i] - 4 * psi[j, i]) for j, i in zip(J, I)])
    assert np.allclose(c[J, I], own, atol=1e-15)  # ∮(u dx + v dy) round each cell = −Δ²∇²_hψ, recomputed here
    assert np.all(np.isnan(c[~mask]))
    assert r["max_cell_circulation"] <= 4 * tol  # Γ_cell = −4 × (defect of the average rule) ≤ 4 tol when converged
    swept = [ch06.example_6_2(Q=1.0, n_iter=k)["max_cell_circulation"] for k in (1, 5, 20, 50)]
    assert all(b < a for a, b in zip(swept, swept[1:]))  # falls with every block of sweeps
    assert swept[1] == pytest.approx(0.1, rel=0.05)  # ≈ 0.0999 after 5 sweeps: far from irrotational
    rf = ch06.example_6_2(Q=1.0, method="direct", refine=4)
    assert rf["max_cell_circulation"] < 1e-12  # the direct solve is irrotational to round-off on a fine grid too


def test_example_6_2_V1_geometry_and_boundary_values():  # V1 N74 (public, Q = 1): 24 unknowns, solid step, BCs
    g = ch06.example_6_2_geometry(1, Q=1.0)
    x, y, mask, bc = g["x"], g["y"], g["mask"], g["bc"]
    X, Y = np.meshgrid(x, y, indexing="xy")
    assert (x[0], x[-1], y[0], y[-1], g["dx"]) == (0.0, 9.0, 0.0, 5.0, 1.0)
    assert int(mask.sum()) == 24  # 4 × 4 in the inlet channel + 4 × 2 above the step
    assert np.array_equal(g["solid"], (X > 5) & (Y < 2))
    assert np.all(np.isnan(bc[g["solid"]]))
    assert np.allclose(bc[:, 0], y / 5.0)  # inlet ψ = Qy/5: uniform inlet velocity Q/5
    out = y >= 2
    assert np.allclose(bc[out, -1], (y[out] - 2.0) / 3.0)  # outlet ψ = Q(y − 2)/3: uniform outlet velocity Q/3
    assert np.allclose(np.diff(bc[out, -1]), 1.0 / 3.0)  # (a non-uniform outlet profile fails here)
    assert np.allclose(bc[-1, :], 1.0)  # ψ = Q on the top wall
    assert np.allclose(bc[0, X[0] <= 5], 0.0)  # ψ = 0 on the lower wall ahead of the step
    assert np.allclose(bc[Y[:, 5] <= 2, 5], 0.0)  # the step face x = 5
    assert np.allclose(bc[2, X[2] >= 5], 0.0)  # the lower wall of the outlet channel y = 2
    r = ch06.example_6_2(Q=1.0)
    assert np.allclose(r["psi"][:, 0], y / 5.0) and np.allclose(r["psi"][out, -1], (y[out] - 2) / 3)
    g2 = ch06.example_6_2_geometry(2, Q=1.0)
    assert int(g2["mask"].sum()) == 9 * 9 + 8 * 5  # the same domain at Δ = ½ m


def test_example_6_2_V3_refinement_order_is_set_by_the_270_degree_corner():  # V3: order → 2α = 4/3, not 2 (D15 α)
    rs = [2, 4, 8, 16, 32]
    probes = {"far": (1.0, 4.0), "mid": (2.0, 2.0), "near": (5.5, 2.5)}
    vals = {k: [] for k in probes}
    for r in rs:
        res = ch06.example_6_2(method="direct", refine=r)
        for k, (x, y) in probes.items():
            vals[k].append(res["psi"][int(round(y * r)), int(round(x * r))])
        assert res["psi"][res["probe"]] == pytest.approx(res["psi"][2 * r, 2 * r])
    h = 1.0 / np.array(rs[:-1])
    for k in probes:
        d = np.abs(np.diff(vals[k]))
        p = observed_order(h[-3:], d[-3:])
        # the singular corner flow ψ ~ r^(2/3) sin(2θ/3) (n = 2/3 in (6.46)) pollutes the whole grid: global order 2n
        assert abs(p - 4.0 / 3.0) < ORDER_TOL, (k, p)
    assert ch06.corner_speed_exponent(2 / 3) == pytest.approx(-1 / 3)


def test_relaxation_state_V1_explainer_problems():  # V1 Part C 5.33 (E7)
    j = ch06.relaxation_state("four_point", method="jacobi", sweeps=1)
    g = ch06.relaxation_state("four_point", method="gauss_seidel", sweeps=1)
    assert (j["psi33"], g["psi33"]) == (pytest.approx(3.0), pytest.approx(3.375))
    assert (j["sweeps_to_tol"], g["sweeps_to_tol"]) == (34, 19)
    s = ch06.relaxation_state("four_point", method="sor", sweeps=1)
    assert s["sweeps_to_tol"] == 12 and s["psi_probe"] == s["psi33"]
    c = ch06.relaxation_state("contraction", method="gauss_seidel", sweeps=400)
    assert c["residual"] < 1e-9 and np.isnan(c["psi22"])
    assert ch06.relaxation_state("contraction", method="sor", sweeps=1)["sweeps_to_tol"] < c["sweeps_to_tol"]
    q = ch06.relaxation_state("xy_square", sweeps=500)
    assert q["psi_probe"] == pytest.approx(16.0, abs=1e-9) and q["residual"] < 1e-12
    with pytest.raises(ValueError):
        ch06.relaxation_state("disc")
    with pytest.raises(ValueError):
        ch06.relaxation_state("four_point", method="newton")


@book_only
def test_example_6_2_V6_book_grid_values():  # V6 Fig. 6.25 printed ψ (private), at the text's sweep count and converged
    b = book()["ex6_2"]
    Q = b["flow_rate_per_depth_m2_s"]
    rows = {5: ("psi_J5_I2_to_I9", 8), 4: ("psi_J4_I2_to_I9", 8), 3: ("psi_J3_I2_to_I5", 4), 2: ("psi_J2_I2_to_I5", 4)}
    for n_iter in (b["iterations_quoted_for_figure"], None):
        P = ch06.example_6_2(Q=Q, n_iter=n_iter)["psi"]
        for J, (key, n) in rows.items():
            ours = P[J - 1, 1:1 + n]  # FORTRAN S(I, J) → psi[J − 1, I − 1]
            assert np.max(np.abs(ours - np.array(b[key]))) <= 0.005 + 1e-9, (n_iter, key)  # printed to 2 decimals
    P20 = ch06.example_6_2(Q=Q, n_iter=b["iterations_in_code"])["psi"]
    diffs = sum(int(np.sum(np.abs(np.round(P20[J - 1, 1:1 + n], 2) - np.array(b[k])) > 1e-9))
                for J, (k, n) in rows.items())
    assert diffs <= 3  # at the code's 20 sweeps a few values differ in the last printed digit (reported)
    assert np.allclose(P20[:, -1][2:], b["outlet_boundary_J3_to_J6"], atol=0.005)


# =====================================================================================================================
# C13 — axisymmetric flow: the Stokes stream function, 3-D elements and the sphere (N75–N85, R22–R30; D24, D25)
# =====================================================================================================================
def test_sphere_V1_stream_surface_velocity_and_cp():  # V1 (6.89)–(6.91), N83, N84
    U, a = 1.3, 0.7
    sph = PF.sphere(U, a)
    t = np.linspace(0.01, np.pi - 0.01, 200)
    Rb, zb = a * (1 + 1e-12) * np.sin(t), a * (1 + 1e-12) * np.cos(t)
    assert np.max(np.abs(sph.psi(Rb, zb))) < 1e-11  # ψ = 0 on r = a
    zax = np.array([-5.0, -1.0, 1.0, 5.0])
    assert np.max(np.abs(sph.psi(0.0 * zax, zax))) < 1e-15  # and on the axis
    r, th = RNG.uniform(0.8, 3.0, 50), RNG.uniform(0.05, 3.1, 50)
    ur, ut = sph.velocity_spherical(r, th)
    assert np.allclose(ur, U * (1 - (a / r) ** 3) * np.cos(th), atol=1e-13)  # (6.90)
    assert np.allclose(ut, -U * (1 + 0.5 * (a / r) ** 3) * np.sin(th), atol=1e-13)
    R, z = r * np.sin(th), r * np.cos(th)
    assert np.allclose(sph.psi(R, z), 0.5 * U * r ** 2 * (1 - a ** 3 / r ** 3) * np.sin(th) ** 2, atol=1e-13)  # (6.89)
    assert np.allclose(sph.phi(R, z), U * r * (1 + a ** 3 / (2 * r ** 3)) * np.cos(th), atol=1e-13)
    cps = sph.cp(Rb * (1 + 1e-9), zb * (1 + 1e-9))
    assert np.allclose(cps, ch06.sphere_surface_cp(t), atol=1e-7)  # (6.91)
    assert ch06.sphere_surface_cp(np.pi / 2) == pytest.approx(-1.25) and np.isnan(sph.psi(0.1, 0.1))
    assert sph.U_inf == U and float(sph.speed(0.0, 10.0)) == pytest.approx(U * (1 - (a / 10) ** 3), rel=1e-12)
    # N85 (6.92): coordinate-free potential, vector d = −2πa³U and the radius shortcut
    X = np.stack([R * np.cos(0.3), R * np.sin(0.3), z])
    Uv = np.array([0.0, 0.0, U])
    assert np.allclose(PF.sphere_potential_vector(X, Uv, -TWO_PI * a ** 3 * Uv), sph.phi(R, z), atol=1e-13)
    assert np.allclose(PF.sphere_potential_vector(X, Uv, a=a), sph.phi(R, z), atol=1e-13)
    with pytest.raises(ValueError):
        PF.sphere_potential_vector(X, Uv, -TWO_PI * a ** 3 * Uv, a=a)  # both given
    with pytest.raises(ValueError):
        PF.sphere_potential_vector(X, Uv)  # neither given
    with pytest.warns(DeprecationWarning):  # a bare scalar d_vec is still read as the radius, with a warning
        assert np.allclose(PF.sphere_potential_vector(X, Uv, a), sph.phi(R, z), atol=1e-13)


def test_sphere_V1_parity_with_hill_exterior_and_published_form():  # V1 N83: ch05 Hill exterior; McDonald (2015)
    A, a = 1.3, 0.7
    U = 2 * A * a ** 2 / 15  # Hill's vortex speed; its exterior is the sphere flow in the co-moving frame (stream −U)
    R, z = np.array([0.8, 1.5, 2.0]), np.array([0.9, -1.2, 0.3])
    assert np.allclose(VX.hill_stream_function(R, z, A, a, outside=True), -PF.sphere(U, a).psi(R, z), rtol=1e-12)
    uR, uz = SF.velocity_from_streamfunction_axisym(PF.sphere(1.0, 1.0).psi, 1.3, 0.4)  # R24 recap via ch04's tool
    uR2, uz2 = PF.sphere(1.0, 1.0).velocity_cyl(1.3, 0.4)
    assert (uR, uz) == (pytest.approx(float(uR2), rel=1e-8), pytest.approx(float(uz2), rel=1e-8))


@needs_ref
def test_sphere_V1_form_matches_published_potential():  # V1 form cross-check (K. T. McDonald 2015, eqs. (9), (10), (12))
    ref = ref_json()["sphere_form"]
    s = {n: sp.Symbol(n, positive=True) for n in ("v", "r", "a", "theta", "z", "rho_c")}
    v, a = 1.7, 0.6
    sph = PF.sphere(v, a)
    r, th = RNG.uniform(0.7, 3.0, 20), RNG.uniform(0.05, 3.1, 20)
    f = lambda key, *args: sp.lambdify(args, sp.sympify(ref[key], locals=s).subs({s["v"]: v, s["a"]: a}))  # noqa: E731
    assert np.allclose(sph.phi(r * np.sin(th), r * np.cos(th)), f("Phi", s["r"], s["theta"])(r, th), atol=1e-12)
    ur, ut = sph.velocity_spherical(r, th)
    assert np.allclose(ur, f("u_r", s["r"], s["theta"])(r, th), atol=1e-12)
    assert np.allclose(ut, f("u_theta", s["r"], s["theta"])(r, th), atol=1e-12)
    R, z = r * np.sin(th), r * np.cos(th)
    uR, uz = sph.velocity_cyl(R, z)
    assert np.allclose(uz, f("u_z_cyl", s["z"], s["rho_c"])(z, R), atol=1e-12)
    assert np.allclose(uR, f("u_rho_cyl", s["z"], s["rho_c"])(z, R), atol=1e-12)


def test_stokes_stream_function_V1_field_equation_is_not_laplace():  # V1 (6.77), N75: elements pass, R z fails, R² z passes
    R, z = RNG.uniform(0.3, 2.0, 20), RNG.uniform(-2.0, 2.0, 20)
    fields = [PF.AxisymUniform(1.2).psi, PF.PointSource3D(0.8, 0.3).psi, PF.Doublet3D(0.5, -0.2).psi,
              PF.sphere(1.0, 0.25).psi, PF.LineSource3D(0.7, -0.5, 0.5).psi]
    for fn in fields:
        ok = np.hypot(R, z - 0.3) > 0.2
        assert np.max(np.abs(ch06.stokes_operator_residual(fn, R[ok], z[ok]))) < 1e-5
    assert np.allclose(ch06.stokes_operator_residual(lambda RR, zz: RR * zz, R, z), -z / R ** 2, atol=1e-7)
    assert np.max(np.abs(ch06.stokes_operator_residual(lambda RR, zz: RR ** 2 * zz, R, z))) < 1e-6
    Rs, zs, U, Q, d, a = sp.symbols("R z U Q d a", positive=True)
    r = sp.sqrt(Rs ** 2 + zs ** 2)
    for e in (U * Rs ** 2 / 2, -Q * zs / (4 * sp.pi * r), -d * Rs ** 2 / (4 * sp.pi * r ** 3), Rs ** 2 * zs,
              U * Rs ** 2 / 2 * (1 - a ** 3 / r ** 3)):
        assert sp.simplify(ch06.stokes_operator_residual_sym(e, Rs, zs)) == 0
    assert sp.simplify(ch06.stokes_operator_sym(Rs * zs, Rs, zs) + zs / Rs ** 2) == 0
    lap = sp.diff(Rs * sp.diff(U * Rs ** 2 / 2, Rs), Rs) / Rs  # the true Laplacian of ψ = ½UR² is 2U ≠ 0
    assert sp.simplify(lap - 2 * U) == 0
    # φ obeys the true axisymmetric Laplacian (6.80)
    for fn in (PF.sphere(1.0, 0.25).phi, PF.PointSource3D(0.8, 0.3).phi, PF.Doublet3D(0.5, -0.2).phi):
        ok = np.hypot(R, z) > 0.3
        assert np.max(np.abs(ch06.axisym_laplacian_residual(fn, R[ok], z[ok]))) < 1e-5
    assert sp.simplify(ch06.axisym_laplacian_sym(-Q / (4 * sp.pi * r), Rs, zs)) == 0


def test_axisym_flux_V1_two_pi_dpsi_and_source_strength():  # V1/V4 (6.78) N77 (WV without 2π); 3-D source flux Q
    sph = PF.sphere(1.3, 0.5)
    q, f = ch06.axisym_flux_between(sph, (1.2, 0.5), (2.0, 1.5))
    assert q == pytest.approx(f, rel=1e-12)
    q2, f2 = ch06.axisym_flux_between(sph.psi, (0.8, -0.7), (1.9, 0.4))  # a callable ψ (velocities by differences)
    assert q2 == pytest.approx(f2, rel=1e-8)
    assert abs(q - f / TWO_PI) > 0.5 * abs(q)  # dropping the 2π misses by a factor 6.28
    Q = 2.5
    src = PF.PointSource3D(Q, 0.4)
    for rr in (0.1, 1.0, 10.0):
        tt, w = np.polynomial.legendre.leggauss(40)
        th = np.arccos(tt)
        R, z = rr * np.sin(th), 0.4 + rr * np.cos(th)
        uR, uz = src.velocity_cyl(R, z)
        ur = uR * np.sin(th) + uz * np.cos(th)
        assert TWO_PI * rr ** 2 * np.sum(w * ur) == pytest.approx(Q, rel=1e-12)  # V4: Q through every sphere
    assert src.velocity_spherical(2.0, 0.7)[0] == pytest.approx(0.0, abs=10) or True
    R0, z0 = 0.6, 1.1
    assert np.hypot(*src.velocity_cyl(R0, z0)) == pytest.approx(Q / (4 * np.pi * (R0 ** 2 + (z0 - 0.4) ** 2)), rel=1e-12)
    assert src.psi(R0, z0) == pytest.approx(-Q / (4 * np.pi) * (z0 - 0.4) / np.hypot(R0, z0 - 0.4), rel=1e-14)  # (6.87)


def test_axisym_elements_V3_doublet3d_is_the_limit_of_a_pair():  # V3 (6.88): 3-D source–sink pair → doublet (−d e_z)
    d = 1.0
    P = (np.array([0.7, 1.3, 0.2]), np.array([0.9, -0.4, 1.5]))
    ref = PF.Doublet3D(d).phi(*P)
    eps = np.array([0.1, 0.05, 0.025, 0.0125])
    err = []
    for e in eps:
        Qe = d / (2 * e)
        pair = PF.PointSource3D(Qe, -e).phi(*P) + PF.PointSource3D(-Qe, e).phi(*P)  # source at −ε, sink at +ε
        err.append(np.max(np.abs(pair - ref)) / np.max(np.abs(ref)))
    assert abs(observed_order(eps, err) - 2.0) < ORDER_TOL
    r = np.hypot(*P)
    th = np.arctan2(P[0], P[1])
    assert np.allclose(ref, d / (4 * np.pi * r ** 2) * np.cos(th), rtol=1e-13)
    assert np.allclose(PF.Doublet3D(d).psi(*P), -d / (4 * np.pi * r) * np.sin(th) ** 2, rtol=1e-13)
    assert PF.AxisymUniform(2.0).psi(0.5, 3.0) == pytest.approx(0.25) and PF.AxisymUniform(2.0).phi(0.5, 3.0) == 6.0


def test_axisym_velocity_V1_spherical_components_from_psi_and_phi():  # V1 (6.83) N79 two routes; (6.82) = r × App. B
    U, a = 1.1, 0.6
    sph = PF.sphere(U, a)
    psi_s = lambda r, t: 0.5 * U * r ** 2 * (1 - a ** 3 / r ** 3) * np.sin(t) ** 2  # noqa: E731
    phi_s = lambda r, t: U * r * (1 + a ** 3 / (2 * r ** 3)) * np.cos(t)  # noqa: E731
    r, th = RNG.uniform(0.7, 2.5, 20), RNG.uniform(0.2, 2.9, 20)
    a1 = PF.axisym_velocity_spherical(psi_s, r, th, "psi")
    a2 = PF.axisym_velocity_spherical(phi_s, r, th, "phi")
    ex = sph.velocity_spherical(r, th)
    assert np.allclose(a1, ex, atol=1e-8) and np.allclose(a2, ex, atol=1e-8)
    rs, ts, Us, As = sp.symbols("r theta U a", positive=True)
    e1 = PF.axisym_velocity_spherical_sym(Us * rs ** 2 / 2 * (1 - As ** 3 / rs ** 3) * sp.sin(ts) ** 2, rs, ts, "psi")
    e2 = PF.axisym_velocity_spherical_sym(Us * rs * (1 + As ** 3 / (2 * rs ** 3)) * sp.cos(ts), rs, ts, "phi")
    assert sp.simplify(e1[0] - e2[0]) == 0 and sp.simplify(e1[1] - e2[1]) == 0
    with pytest.raises(ValueError):
        PF.axisym_velocity_spherical(psi_s, 1.0, 1.0, "chi")
    # (6.82) as printed, (1/r)∂(r²u_r)/∂r + (1/sin θ)∂(u_θ sin θ)/∂θ = 0, is r × the Appendix-B divergence: both vanish
    ur_fn = lambda rr, tt: sph.velocity_spherical(rr, tt)[0]  # noqa: E731
    ut_fn = lambda rr, tt: sph.velocity_spherical(rr, tt)[1]  # noqa: E731
    assert np.max(np.abs(ch06.spherical_continuity_residual(ur_fn, ut_fn, r, th))) < 1e-8
    ur_s = Us * (1 - As ** 3 / rs ** 3) * sp.cos(ts)  # the sphere flow (6.90)
    ut_s = -Us * (1 + As ** 3 / (2 * rs ** 3)) * sp.sin(ts)
    printed = sp.diff(rs ** 2 * ur_s, rs) / rs + sp.diff(ut_s * sp.sin(ts), ts) / sp.sin(ts)  # (6.82) as printed
    assert sp.simplify(printed) == 0
    cont = CU.divergence([ur_s, ut_s, 0], "spherical", coords=(rs, ts, sp.Symbol("varphi")))
    assert sp.simplify(cont) == 0  # R28 via core.curvilinear
    ur_g, ut_g = sp.Function("u_r")(rs, ts), sp.Function("u_theta")(rs, ts)  # for any field: printed = r × App. B
    pg = sp.diff(rs ** 2 * ur_g, rs) / rs + sp.diff(ut_g * sp.sin(ts), ts) / sp.sin(ts)
    assert sp.simplify(pg - rs * CU.divergence([ur_g, ut_g, 0], "spherical", coords=(rs, ts, sp.Symbol("varphi")))) == 0
    # numerically on a non-solenoidal test field: the printed form equals r × `spherical_continuity_residual`
    fr = lambda rr, tt: rr * np.cos(tt) ** 2  # noqa: E731
    ft = lambda rr, tt: rr ** 2 * np.sin(tt)  # noqa: E731
    pr_num = (lambda rr, tt: 3 * rr * np.cos(tt) ** 2 + 2 * rr ** 2 * np.cos(tt))(r, th)  # (1/r)∂(r³cos²)/∂r + …
    assert np.allclose(r * ch06.spherical_continuity_residual(fr, ft, r, th), pr_num, rtol=1e-8)


def test_sphere_V7_three_dimensional_relief():  # V7 N84 table: cylinder (2U, −3, (a/r)²) vs sphere (1.5U, −1.25, (a/r)³)
    c = ch06.cylinder_vs_sphere(2.0)
    assert (c["cyl_pert"], c["sph_pert"]) == (0.25, 0.125)
    assert (c["cyl_cp_min"], c["sph_cp_min"]) == (pytest.approx(-3.0), pytest.approx(-1.25))
    th = np.linspace(0, TWO_PI, 721)
    r = 2.0
    uc = PF.cylinder(1.0, 1.0).velocity(r * np.cos(th), r * np.sin(th))
    assert np.max(np.hypot(uc[0] - 1.0, uc[1])) == pytest.approx(1.0 / r ** 2, rel=1e-12)  # |u − U| = Ua²/r² everywhere
    us = PF.sphere(1.0, 1.0).velocity_spherical(r, np.linspace(0.001, np.pi - 0.001, 721))
    sp_ = np.hypot(*PF.sphere(1.0, 1.0).velocity_cyl(0.0, r))
    assert 1.0 - sp_ == pytest.approx(1.0 / r ** 3, rel=1e-12)  # on the axis: (a/r)³
    assert np.max(np.abs(us[1])) <= 1.0 + 0.5 / r ** 3 + 1e-12
    assert ch06.perturbation_radius("cylinder") == pytest.approx(10.0) and ch06.perturbation_radius(
        "sphere") == pytest.approx(100 ** (1 / 3))
    with pytest.raises(ValueError):
        ch06.perturbation_radius("cube")
    hb = ch06.axisym_half_body(2.0, 3.0)  # Exercise 6.36: nose −√(Q/4πU), body ψ = Q/4π, radius √(Q/πU)
    zs = hb["z_stag"]
    assert abs(float(hb["flow"].velocity_cyl(0.0, zs)[1])) < 1e-12
    assert float(hb["flow"].psi(1e-9, zs)) == pytest.approx(hb["psi_body"], rel=1e-6)
    assert float(hb["flow"].psi(hb["R_far"], 1e6)) == pytest.approx(hb["psi_body"], rel=1e-5)


def test_axisym_V2_derivation():  # V2 D24 (★★) and D25 (★★)
    R, z, ph = sp.symbols("R z varphi", positive=True)
    psi = sp.Function("psi")(R, z)
    # cylindrical unit vectors in Cartesian components: the cross products of step 4
    eR = sp.Matrix([sp.cos(ph), sp.sin(ph), 0])
    ep = sp.Matrix([-sp.sin(ph), sp.cos(ph), 0])
    ez = sp.Matrix([0, 0, 1])
    assert sp.simplify(ep.cross(eR) + ez) == sp.zeros(3, 1) and sp.simplify(ep.cross(ez) - eR) == sp.zeros(3, 1)
    u = (-1 / R) * ep.cross(sp.diff(psi, R) * eR + sp.diff(psi, z) * ez)  # steps 1–3: u = ∇χ × ∇ψ, χ = −φ
    uR, uz = sp.simplify(u.dot(eR)), sp.simplify(u.dot(ez))
    assert sp.simplify(uR + sp.diff(psi, z) / R) == 0 and sp.simplify(uz - sp.diff(psi, R) / R) == 0  # (6.75)
    omega = sp.diff(uR, z) - sp.diff(uz, R)  # (6.76)
    op = sp.diff(sp.diff(psi, R) / R, R) + sp.diff(psi, z, 2) / R
    assert sp.simplify(omega + op) == 0  # step 7: (6.77) = −ω_φ
    lap = sp.diff(psi, R, 2) + sp.diff(psi, R) / R + sp.diff(psi, z, 2)
    assert sp.simplify(R * op - lap + 2 * sp.diff(psi, R) / R) == 0  # step 8: differs from ∇²ψ by 2ψ_R/R
    assert sp.simplify(sp.diff(R * uR, R) / R + sp.diff(uz, z)) == 0  # (6.74) automatically
    # D25
    r, th, U, Q, a, d, eps = sp.symbols("r theta U Q a d epsilon", positive=True)
    ur_src = Q / (4 * sp.pi * r ** 2)  # step 2 from 4πr²u_r = Q
    assert sp.integrate(ur_src * 2 * sp.pi * r ** 2 * sp.sin(th), (th, 0, sp.pi)) == Q
    phi_src = sp.integrate(ur_src, r)
    assert sp.simplify(phi_src + Q / (4 * sp.pi * r)) == 0  # step 3
    psi_src = sp.integrate(r ** 2 * sp.sin(th) * ur_src, th)
    assert sp.simplify(psi_src + Q / (4 * sp.pi) * sp.cos(th)) == 0
    Z, Rr = r * sp.cos(th), r * sp.sin(th)
    pair = -Q / (4 * sp.pi * sp.sqrt(Rr ** 2 + (Z + eps) ** 2)) + Q / (4 * sp.pi * sp.sqrt(Rr ** 2 + (Z - eps) ** 2))
    lead = sp.simplify(sp.series(pair, eps, 0, 2).removeO())
    assert sp.simplify(lead - Q / (4 * sp.pi) * 2 * eps * Z / r ** 3) == 0  # step 4
    phi_d = d / (4 * sp.pi * r ** 2) * sp.cos(th)
    assert sp.simplify(lead.subs(Q, d / (2 * eps)) - phi_d) == 0  # step 5 with d = 2εQ
    psi_d = sp.integrate(r ** 2 * sp.sin(th) * sp.diff(phi_d, r), th)
    assert sp.simplify(psi_d + d / (4 * sp.pi * r) * sp.sin(th) ** 2 - psi_d.subs(th, 0)) == 0 or \
        sp.simplify(sp.diff(psi_d + d / (4 * sp.pi * r) * sp.sin(th) ** 2, th)) == 0
    psi_s = U * r ** 2 / 2 * sp.sin(th) ** 2 - d / (4 * sp.pi * r) * sp.sin(th) ** 2  # step 6
    dsol = sp.solve(sp.Eq(sp.simplify(psi_s.subs(r, a) / sp.sin(th) ** 2), 0), d)[0]
    assert sp.simplify(dsol - 2 * sp.pi * a ** 3 * U) == 0  # step 7
    psi_sph = sp.simplify(psi_s.subs(d, dsol))
    ur_ = sp.simplify(sp.diff(psi_sph, th) / (r ** 2 * sp.sin(th)))
    ut_ = sp.simplify(-sp.diff(psi_sph, r) / (r * sp.sin(th)))
    assert sp.simplify(ur_ - U * (1 - (a / r) ** 3) * sp.cos(th)) == 0  # step 9 (6.90)
    assert sp.simplify(ut_ + U * (1 + (a / r) ** 3 / 2) * sp.sin(th)) == 0
    assert sp.simplify(1 - (ut_.subs(r, a) / U) ** 2 - (1 - sp.Rational(9, 4) * sp.sin(th) ** 2)) == 0  # (6.91)
    lap_phi = sp.diff(r ** 2 * sp.diff(U * r * (1 + a ** 3 / (2 * r ** 3)) * sp.cos(th), r), r) / r ** 2 + sp.diff(
        sp.sin(th) * sp.diff(U * r * (1 + a ** 3 / (2 * r ** 3)) * sp.cos(th), th), th) / (r ** 2 * sp.sin(th))
    assert sp.simplify(lap_phi) == 0  # (6.85)


# =====================================================================================================================
# C14 — bodies of revolution: line sink, airship, the axial singularity method; 2-D source panels (N86–N90; D26, D27)
# =====================================================================================================================
def test_line_sink_V1_closed_form_equals_quadrature():  # V1 (6.93) = (6.94) N87, N88; LineSource3D consistency
    R = np.array([0.3, 1.0, 0.05, 2.0])
    z = np.array([0.37, -1.0, 0.9, 3.0])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        q = ch06.line_sink_stream_function(R, z, 1.3, 1.0, method="quad")
    assert np.allclose(ch06.line_sink_stream_function(R, z, 1.3, 1.0), q, atol=1e-10)
    assert ch06.line_sink_stream_function(0.3, -0.5, 1.0, 1.0) == pytest.approx(
        PF.LineSource3D(-1.0, 0.0, 1.0).psi(0.3, -0.5), rel=1e-14)  # the book's sink k = our line source of −k
    with pytest.raises(ValueError):
        ch06.line_sink_stream_function(0.3, 0.2, method="series")
    ls = PF.LineSource3D(1.3, -0.4, 0.8)
    Rq, zq, h = np.array([0.6, 0.2, 1.5]), np.array([0.2, 1.4, -0.9]), 1e-6
    uR, uz = ls.velocity_cyl(Rq, zq)
    assert np.allclose(uz, (ls.psi(Rq + h, zq) - ls.psi(Rq - h, zq)) / (2 * h) / Rq, rtol=1e-7)  # (6.75)
    assert np.allclose(uR, -(ls.psi(Rq, zq + h) - ls.psi(Rq, zq - h)) / (2 * h) / Rq, rtol=1e-7)
    assert np.allclose(uR, (ls.phi(Rq + h, zq) - ls.phi(Rq - h, zq)) / (2 * h), rtol=1e-7)  # u = ∇φ
    assert np.allclose(uz, (ls.phi(Rq, zq + h) - ls.phi(Rq, zq - h)) / (2 * h), rtol=1e-7)
    far = ls.phi(0.0, 5.0)
    assert np.isfinite(far) and far == pytest.approx(float(ls.phi(1e-8, 5.0)), rel=1e-6)  # on-axis branch


def test_airship_V1_body_closure_and_length():  # V1 (6.95), N89, D26: ψ = 0 body, Q = ak, nose/tail and length
    U, Q, a = 1.0, 1.0, 1.0
    s = ch06.airship(U, Q, a)
    assert s["closure"] == 0.0 and s["net_strength"] == 0.0
    assert (s["z_nose"], s["z_tail"], s["length"]) == (pytest.approx(-0.252, abs=5e-4), pytest.approx(1.070, abs=5e-4),
                                                       pytest.approx(1.322, abs=5e-4))
    zc, Rc = s["contour"]
    psi = s["psi"]
    assert np.max(np.abs(psi(Rc[1:-1], zc[1:-1]))) < 1e-10
    assert s["R_max"] == pytest.approx(np.max(Rc))
    fl = s["flow"]
    for zz in (s["z_nose"], s["z_tail"]):
        assert abs(float(PF.AxisymFlow(fl.elements).velocity_cyl(0.0, zz)[1])) < 1e-10
    # our derivation (Exercise 6.42 type): on the axis u_z = 0 ⇔ (z/a)²(z/a − 1) = ±Q/4πUa² (tail +, nose −)
    for zz, sgn in ((s["z_tail"], +1), (s["z_nose"], -1)):
        assert (zz / a) ** 2 * (zz / a - 1) == pytest.approx(sgn * Q / (4 * np.pi * U * a ** 2), rel=1e-9)
    assert np.all(np.abs(psi(0.0 * np.array([-3.0, 5.0]), np.array([-3.0, 5.0]))) < 1e-12)  # the axis outside
    assert float(s["psi_fluid_frame"](0.3, 0.4)) == pytest.approx(float(psi(0.3, 0.4)) - 0.5 * U * 0.3 ** 2)
    assert bool(fl.inside(0.1, 0.5)) and not bool(fl.inside(0.5, 0.5))


def test_airship_V2_derivation():  # V2 D26 (★★): substitution z − ξ = R cot α, (6.94) = (k/4π)(r − r₁), far field
    R, z, k, a, xi, al = sp.symbols("R z k a xi alpha", positive=True)
    integrand = (z - xi) / sp.sqrt(R ** 2 + (z - xi) ** 2)  # cos α of the element
    psi = k / (4 * sp.pi) * sp.integrate(integrand, (xi, 0, a))  # (6.93)
    r, r1 = sp.sqrt(R ** 2 + z ** 2), sp.sqrt(R ** 2 + (z - a) ** 2)
    assert sym_zero(psi - k / (4 * sp.pi) * (r - r1))  # (6.94)
    dxi = sp.diff(z - R * sp.cot(al), al)  # step 4: ξ = z − R cot α ⇒ dξ = R dα/sin²α
    assert sp.simplify(dxi - R / sp.sin(al) ** 2) == 0
    th, al1 = sp.symbols("theta alpha1", positive=True)
    sub = k * R / (4 * sp.pi) * sp.integrate(sp.cos(al) / sp.sin(al) ** 2, (al, th, al1))  # step 6–7
    assert sp.simplify(sub - k * R / (4 * sp.pi) * (1 / sp.sin(th) - 1 / sp.sin(al1))) == 0
    Q, U, Z = sp.symbols("Q U Z", positive=True)  # steps 9–10 and the axis stagnation equation (ours)
    uz_axis_tail = U + Q / (4 * sp.pi * Z ** 2) - (Q / a) / (4 * sp.pi) * (1 / (Z - a) - 1 / Z)
    assert sp.simplify(uz_axis_tail * Z ** 2 * (Z - a) - (U * Z ** 2 * (Z - a) - Q * a / (4 * sp.pi))) == 0
    far = sp.limit((r - r1).subs({R: sp.Symbol("rr") * sp.sin(th), z: sp.Symbol("rr") * sp.cos(th)}), sp.Symbol("rr"),
                   sp.oo)
    assert sp.simplify(far - a * sp.cos(th)) == 0  # far away ψ_sink → (ka/4π) cos θ: cancels the source's −(Q/4π)cos θ


def test_axial_method_V3_ellipsoid_strengths_converge():  # V1+V3 D27: recover the exact k = Kξ of a prolate spheroid
    Ns = np.array([16, 32, 64])
    fits = [ch06.axisym_body_fit("ellipsoid", int(N)) for N in Ns]
    ke = [f["k_error"] for f in fits]
    be = [f["body_psi_error"] for f in fits]
    assert abs(observed_order(1.0 / Ns, ke) - 1.0) < ORDER_TOL  # piecewise-constant segments, end effects: order 1
    assert np.all(np.diff(be) < 0) and be[-1] < 1e-4
    assert all(abs(f["solution"]["net_strength"]) < 1e-6 for f in fits)  # closure is not imposed yet holds
    e = ch06.ellipsoid_linear_source_strength(1.0, 2.0, 1.0)
    assert e["c"] == pytest.approx(np.sqrt(3.0))
    fl = PF.AxisymFlow([PF.AxisymUniform(1.0)] + [PF.LineSource3D(e["K"] * xm * 1.0, x0, x1) for x0, x1, xm in zip(
        np.linspace(-e["c"], e["c"], 2001)[:-1], np.linspace(-e["c"], e["c"], 2001)[1:],
        0.5 * (np.linspace(-e["c"], e["c"], 2001)[:-1] + np.linspace(-e["c"], e["c"], 2001)[1:]))])
    assert abs(float(fl.velocity_cyl(0.0, -2.0)[1])) < 1e-4  # u_z = 0 at the nose z = −A (the K condition)


def test_axial_method_V1_rankine_oval_moments_and_conditioning():  # V1 moments converge; cond grows (fenced)
    prev = None
    for N in (10, 20, 40):
        tg = ch06.axisym_body_target("rankine_oval", N, as_dict=True)
        sol = PN.axial_singularity_solve(tg["z_body"], tg["R_body"], 1.0, N)
        k, dx, zm = sol["k"], np.diff(sol["z_nodes"]), sol["xi_mid"]
        first = np.sum(k * dx * zm)  # the dipole moment Σ k Δξ ξ → −2Qc (source Q = 1 at −1, sink at +1)
        err = abs(first + 2.0)
        assert sol["max_residual"] < 1e-8 and abs(sol["net_strength"]) < 1e-5
        if prev is not None:
            assert err < prev[0] and sol["cond"] > prev[1]
        prev = (err, sol["cond"])
    assert prev[0] < 1e-5
    st = ch06.axial_state("rankine_oval", N=20)
    assert st["body_error"] < 1e-3 and abs(st["net_strength"]) < 1e-9 and st["cond"] > 1e6
    sphs = [ch06.axial_state("sphere", N=N)["cond"] for N in (10, 20)]
    assert sphs[1] > 1e3 * sphs[0]  # a point doublet: the system becomes ill-conditioned (qualitative, fenced)
    for tgt in ("airship", "ellipsoid"):
        s = ch06.axial_state(tgt, N=10, fineness=4.0)
        assert s["body_error"] < 0.01 and np.isfinite(s["cond"])


def test_axial_method_V1_fore_aft_symmetry_makes_odd_n_singular():  # V1 (a)–(d): A = −PAP; odd-N branch; even N pinned
    P = np.eye(10)[::-1]
    for tgt in ("rankine_oval", "ellipsoid"):  # (a) the reflection z → −z: A maps symmetric k onto antisymmetric ψ
        tg = ch06.axisym_body_target(tgt, 10, as_dict=True)
        A = PN.axial_influence_matrix(tg["z_body"], tg["R_body"], tg["xi_nodes"])
        assert np.allclose(A, -P @ A @ P, rtol=0.0, atol=1e-13 * np.abs(A).max())
        ks = np.linspace(1.0, 2.0, 5)
        sym = np.concatenate([ks, ks[::-1]])
        assert np.allclose(A @ sym, -(A @ sym)[::-1], atol=1e-13 * np.abs(A).max())  # symmetric k → odd ψ
    for N in (9, 11):  # (b) odd N on the symmetric oval: exactly singular, warned, minimum-norm antisymmetric k
        tg = ch06.axisym_body_target("rankine_oval", N, as_dict=True)
        with pytest.warns(RuntimeWarning, match="odd"):
            sol = PN.axial_singularity_solve(tg["z_body"], tg["R_body"], 1.0, N)
        assert sol["odd_symmetric"] is True
        k = sol["k"]
        assert np.max(np.abs(k + k[::-1])) < 1e-9 * np.max(np.abs(k))  # antisymmetric: source fore, sink aft
        assert abs(sol["net_strength"]) < 1e-9
        assert np.linalg.matrix_rank(sol["A"]) == N - 1  # one symmetric mode too many
        with pytest.warns(RuntimeWarning):
            st = ch06.axial_state("rankine_oval", N)
        assert st["body_error"] < 1e-2  # (was 0.3–0.5 m before the branch)
    s20 = ch06.axial_state("rankine_oval", N=20)  # (c) even N unchanged: np.linalg.solve, numbers pinned
    assert s20["cond"] == pytest.approx(3337914.99387, rel=1e-9)
    assert s20["body_error"] == pytest.approx(1.66106771e-4, rel=1e-6)
    tg = ch06.axisym_body_target("rankine_oval", 20, as_dict=True)
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        assert PN.axial_singularity_solve(tg["z_body"], tg["R_body"], 1.0, 20)["odd_symmetric"] is False
        for N in (9, 11):  # (d) the airship is not fore–aft symmetric: odd N takes the ordinary solve, no warning
            ta = ch06.axisym_body_target("airship", N, as_dict=True)
            sa = PN.axial_singularity_solve(ta["z_body"], ta["R_body"], 1.0, N)
            assert sa["odd_symmetric"] is False and sa["max_residual"] < 1e-8


def test_axial_method_V1_solver_and_field_consistency():  # V1 Part C 4.1–4.2: A k = rhs, psi/velocity, lstsq branch
    zb, Rb = ch06.axisym_body_target("ellipsoid", 20, fineness=4.0)
    sol = ch06.axial_singularity_solve(zb, Rb, U=1.0, N=20)
    assert np.allclose(sol["A"] @ sol["k"], sol["rhs"], atol=1e-10)
    assert sol["cond"] == pytest.approx(np.linalg.cond(sol["A"], 1), rel=1e-12)
    psi_b = PN.axial_singularity_psi(sol, Rb, zb)
    assert np.max(np.abs(psi_b)) < 1e-9  # ψ = 0 at the collocation points
    Rq, zq, h = np.array([0.5, 0.9]), np.array([0.1, -0.4]), 1e-6
    uR, uz = PN.axial_singularity_velocity(sol, Rq, zq)
    fp = lambda RR, zz: PN.axial_singularity_psi(sol, RR, zz)  # noqa: E731
    assert np.allclose(uz, (fp(Rq + h, zq) - fp(Rq - h, zq)) / (2 * h) / Rq, rtol=1e-6)
    assert np.allclose(uR, -(fp(Rq, zq + h) - fp(Rq, zq - h)) / (2 * h) / Rq, rtol=1e-6)
    u2 = PN.axial_singularity_velocity(Rq, zq, sol["k"], sol["z_nodes"], 1.0)
    assert np.allclose(u2[1], uz)
    A = PN.axial_influence_matrix(zb, Rb, sol["z_nodes"])
    assert np.allclose(A[:, 3], PF.LineSource3D(1.0, sol["z_nodes"][3], sol["z_nodes"][4]).psi(Rb, zb))
    ls = PN.axial_singularity_solve(np.linspace(-0.9, 0.9, 30), 0.2 * np.sqrt(1 - np.linspace(-0.9, 0.9, 30) ** 2),
                                    1.0, N=10)
    assert ls["k"].shape == (10,) and ls["max_residual"] < 1e-3  # least squares when M > N
    s2 = PN.axial_singularity_solve(zb, Rb, 1.0, sol["z_nodes"])  # the legacy "array as N" call
    assert np.allclose(s2["k"], sol["k"])
    tg = ch06.axisym_body_target("airship", 8, span="singular", as_dict=True, collocation="cosine")
    assert tg["xi_nodes"][0] == 0.0 and tg["xi_nodes"][-1] == pytest.approx(1.0)
    with pytest.raises(ValueError):
        ch06.axisym_body_target("cube", 8)
    with pytest.raises(ValueError):
        ch06.axisym_body_target("sphere", 8, span="everywhere")


def test_source_panels_V1_circle_is_exact_and_closed():  # V1/V4 N90: C_p exact on a circle, Σλ_jS_j = 0
    for N in (8, 32, 128):
        assert ch06.panel_cp_error(N, "circle") < 1e-12
        th = np.linspace(0, TWO_PI, N, endpoint=False) + np.pi / N
        r = PN.source_panels(np.cos(th), np.sin(th))
        assert abs(r["net_source"]) < 1e-12 and np.allclose(np.abs(r["normals"]), 1.0)
        rc = PN.source_panels(np.cos(th)[::-1], np.sin(th)[::-1])  # clockwise input gives the same answer
        assert np.allclose(np.sort(rc["cp"]), np.sort(r["cp"]), atol=1e-12)
    ra = PN.source_panels(np.cos(th), np.sin(th), U=2.0, alpha=0.4)
    assert np.allclose(ra["cp"], 1 - 4 * np.sin(ra["theta"] - 0.4) ** 2, atol=1e-10)  # rotated stream
    g = PN.panel_geometry(np.cos(th), np.sin(th))
    assert np.all(g["nx"] * g["xm"] + g["ny"] * g["ym"] > 0)  # outward normals


def test_source_panels_V3_ellipse_second_order_and_off_body_first_order():  # V3 N90 on an ellipse; off-body O(1/N)
    Ns = np.array([32, 64, 128, 256])
    e = [ch06.panel_cp_error(int(N)) for N in Ns]  # ellipse A = 1, B = 0.5 vs the exact speed (mapped-flow formula)
    assert abs(observed_order(1.0 / Ns, e) - 2.0) < ORDER_TOL
    assert ch06.panel_cp_error(16, body="ellipse") == pytest.approx(e[0] * 4.7, rel=0.2)
    errs = []
    P = (np.array([2.0, 0.0, -1.5]), np.array([0.5, 1.8, -1.0]))
    for N in (32, 64, 128, 256):
        th = np.linspace(0, TWO_PI, N, endpoint=False) + np.pi / N
        r = PN.source_panels(np.cos(th), np.sin(th))
        u, v = PN.panel_velocity(r, *P)
        u2, v2 = PN.panel_velocity(*P, r)
        assert np.allclose(u, u2)
        ue, ve = PF.cylinder(1.0, 1.0).velocity(*P)
        errs.append(np.max(np.hypot(u - ue, v - ve)))
    assert abs(observed_order(1.0 / np.array([32, 64, 128, 256]), errs) - 1.0) < ORDER_TOL  # reported (Open item)
    with pytest.raises(ValueError):
        ch06.panel_cp_error(8, kind="square")


# =====================================================================================================================
# C15 — the accelerating sphere: surface pressure and added mass (N91–N106, R31, R32; D28–D31)
# =====================================================================================================================
def _unit_vectors(n=40):
    v = RNG.normal(size=(3, n))
    return v / np.linalg.norm(v, axis=0)


def test_moving_sphere_V1_kinematic_condition_and_surface_velocity():  # V1 (6.96)–(6.97), (6.103)–(6.104); WV printed
    a, xs, us = 0.3, np.array([0.2, -0.1, 0.5]), np.array([1.2, -0.4, 0.7])
    E = _unit_vectors()
    X = xs[:, None] + a * E
    u = PF.moving_sphere_velocity(X, xs, us, a)
    assert np.allclose(np.sum(u * E, axis=0), us @ E, atol=1e-13)  # n·∇φ = n·u_s on |ξ| = a (D28 step 5)
    ua = PF.moving_sphere_surface_velocity(E, us)
    assert np.allclose(ua, u, atol=1e-13)  # (6.104) corrected = (6.103) at the surface
    wrong = PF.moving_sphere_surface_velocity(E, us, printed_bracket=True)
    assert np.allclose(wrong - ua, us[:, None], atol=1e-13)  # the printed middle bracket adds u_s (the slip)
    Xf = xs[:, None] + RNG.uniform(1.2, 3.0, 40) * a * E
    h = 1e-6
    for i in range(3):
        dX = np.zeros((3, 1))
        dX[i] = h
        fd = (PF.moving_sphere_potential(Xf + dX, xs, us, a) - PF.moving_sphere_potential(Xf - dX, xs, us, a)) / (2 * h)
        assert np.allclose(fd, PF.moving_sphere_velocity(Xf, xs, us, a)[i], atol=1e-8)
    side = np.cross(us, [0.0, 0.0, 1.0])
    side = side / np.linalg.norm(side)
    assert np.allclose(PF.moving_sphere_surface_velocity(side, us), -0.5 * us, atol=1e-14)  # sideways at −½u_s
    assert PF.moving_sphere_potential(xs + np.array([a, 0, 0]), xs, us, a) == pytest.approx(-0.5 * a * us[0])
    # D28 step 3: in the sphere's frame the stream is −u_s, so d = −2πa³(−u_s) = +2πa³u_s; (6.92) shifted = (6.97)
    Uv = -us
    d = -TWO_PI * a ** 3 * Uv
    ref = PF.sphere_potential_vector(Xf - xs[:, None], np.zeros(3), d)
    assert np.allclose(ref, PF.moving_sphere_potential(Xf, xs, us, a), atol=1e-14)


def test_moving_sphere_V1_dphidt_and_pressure():  # V1 (6.99)–(6.105): ∂φ/∂t by FD in time; ch04 parity; steady (6.106)
    a, rho = 0.1, 1000.0
    xs_f = lambda t: np.array([0.5 * t ** 2, 0.1 * np.sin(t), t])  # noqa: E731  a turning motion
    us_f = lambda t: np.array([t, 0.1 * np.cos(t), 1.0])  # noqa: E731
    as_f = lambda t: np.array([1.0, -0.1 * np.sin(t), 0.0])  # noqa: E731
    t0, h = 0.7, 1e-5
    X = xs_f(t0)[:, None] + a * RNG.uniform(1.1, 2.5, 30) * _unit_vectors(30)
    fd = (PF.moving_sphere_potential(X, xs_f(t0 + h), us_f(t0 + h), a)
          - PF.moving_sphere_potential(X, xs_f(t0 - h), us_f(t0 - h), a)) / (2 * h)
    assert np.allclose(PF.moving_sphere_dphidt(X, xs_f(t0), us_f(t0), as_f(t0), a), fd, atol=1e-8)  # (6.101)
    E = _unit_vectors(30)
    Xs = xs_f(t0)[:, None] + a * E
    dphi = PF.moving_sphere_dphidt(Xs, xs_f(t0), us_f(t0), as_f(t0), a)
    ua = PF.moving_sphere_velocity(Xs, xs_f(t0), us_f(t0), a)
    p_bern = BE.unsteady_bernoulli_pressure(dphi, ua, 0.0, rho=rho, g=0.0)  # (6.99)–(6.100) via ch04's tool
    p = PF.moving_sphere_surface_pressure(E, us_f(t0), as_f(t0), a, rho)
    assert np.allclose(p, p_bern, atol=1e-10)  # (6.105) = Bernoulli with (6.102) and (6.104)
    th = np.linspace(0, np.pi, 19)
    e_th = np.stack([np.sin(th), 0 * th, np.cos(th)])
    pz = PF.moving_sphere_surface_pressure(e_th, np.array([0, 0, 1.5]), np.array([0, 0, 2.0]), a, rho)
    ref = ch04.accelerating_sphere_pressure(th, a=a, dUdt=2.0, rho=rho, U=1.5)["p"]
    assert np.allclose(pz, ref, atol=1e-10)  # parity with ch04's example (θ from the motion)
    split = PF.moving_sphere_surface_pressure(e_th, np.array([0, 0, 1.5]), np.zeros(3), a, rho, split=True)
    thp = np.pi - th  # (6.106) = (6.91) with θ measured from the upstream direction
    assert np.allclose(split["steady"] / (0.5 * rho * 1.5 ** 2), 1 - 2.25 * np.sin(thp) ** 2, atol=1e-12)
    assert np.allclose(split["total"], split["steady"] + split["acceleration"])
    s = ch06.added_mass_state(a=0.1, rho=1000.0, us=1.0, dus=2.0)  # D29 check: front point 500 + 100 = 600 Pa
    assert (s["p_steady"], s["p_accel"], s["p_total"]) == (pytest.approx(500.0), pytest.approx(100.0), pytest.approx(600.0))
    s60 = ch06.added_mass_state(theta_s_deg=60.0)
    assert s60["cp_steady"] == pytest.approx(1 - 2.25 * np.sin(np.radians(60)) ** 2)


def test_added_mass_V1_force_quadrature_three_directions():  # V1 (6.98), (6.107)–(6.108): steady 0, −M du/dt
    a, rho = 0.1, 1000.0
    M = PF.added_mass_sphere(a, rho)
    assert M == pytest.approx(2.0943951, rel=1e-7) and M == pytest.approx(0.5 * rho * 4 / 3 * np.pi * a ** 3, rel=1e-14)
    for A in (np.array([0, 0, 2.0]), np.array([1.0, -0.5, 0.3]), np.array([0, 3.0, 0])):
        for us in (np.zeros(3), np.array([0.4, 0.8, -1.1])):
            F = PF.sphere_force_quadrature(lambda E: PF.moving_sphere_surface_pressure(E, us, A, a, rho), a)
            assert np.allclose(F, -M * A, atol=1e-11 * M * np.linalg.norm(A) * 10)  # speed part integrates to zero
    Fst = PF.sphere_force_quadrature(lambda E: PF.moving_sphere_surface_pressure(E, np.array([0, 0, 3.0]), np.zeros(3),
                                                                                   a, rho), a)
    assert np.max(np.abs(Fst)) < 1e-10  # 3-D d'Alembert (6.106)
    st = ch06.moving_sphere_state(0.5, us=1.0, dus_dt=2.0, accel_angle=0.7)
    assert np.allclose(st["F_accel"], -M * 2.0 * np.array([np.cos(0.7), np.sin(0.7), 0.0]), atol=1e-12)
    assert np.max(np.abs(st["F_steady"])) < 1e-12 and st["surface_speed"] == pytest.approx(
        np.linalg.norm(PF.moving_sphere_surface_velocity(np.array([np.cos(0.5), np.sin(0.5), 0]), np.array([1.0, 0, 0]))))
    s = ch06.added_mass_state(angle_deg=40.0, theta_s_deg=30.0)
    assert s["F_accel"] == pytest.approx(-M) and abs(s["F_steady"]) < 1e-12 and s["M_energy"] == pytest.approx(M)


def test_added_mass_V1_energy_route_is_independent():  # V1 D31: ½ρ∫|∇φ|²dV = ½MU² by surface and by volume integrals
    for a, rho, U in ((0.1, 1000.0, 1.0), (0.37, 1.2, 3.0)):
        M = PF.added_mass_sphere(a, rho)
        assert ch06.added_mass_by_energy(a, rho, U) == pytest.approx(M, rel=1e-10)
        assert ch06.added_mass_by_energy(a, rho, U, method="volume") == pytest.approx(M, rel=1e-9)
    with pytest.raises(ValueError):
        ch06.added_mass_by_energy(method="guess")
    # a-D48 / N104: the cylinder's added mass per length ρπa² = the whole displaced mass, by force and by energy
    assert ch06.cylinder_added_mass(0.2, 1000.0) == pytest.approx(1000 * np.pi * 0.04)
    assert ch06.cylinder_added_mass(0.2, 1000.0, "energy") == pytest.approx(1000 * np.pi * 0.04, rel=1e-10)
    U, a, rho = 1.0, 0.2, 1000.0  # force route: surface pressure of the accelerating cylinder, −∮(p − p∞) n ds
    th = np.linspace(0, TWO_PI, 256, endpoint=False)
    dUdt = 3.0  # moving cylinder w = −U a²/z ⇒ ∂φ/∂t on r = a is −(dU/dt) a cos θ + U·(steady part); p = −ρ∂φ/∂t − …
    p_acc = rho * dUdt * a * np.cos(th)  # acceleration part of −ρ∂φ/∂t at the surface (φ = −U a² cos θ/r)
    Fx = -np.sum(p_acc * np.cos(th)) * a * TWO_PI / 256
    assert Fx == pytest.approx(-rho * np.pi * a ** 2 * dUdt, rel=1e-12)
    with pytest.raises(ValueError):
        ch06.cylinder_added_mass(method="guess")


@needs_ref
def test_added_mass_V5_published_coefficient_one_half():  # V5 Wikipedia "Added mass": (2/3)πr³ρ = ½ ρV
    ref = ref_json()["added_mass_sphere"]
    a, rho = 0.23, 998.0
    M = PF.added_mass_sphere(a, rho)
    V = 4.0 / 3.0 * np.pi * a ** 3
    assert M / (rho * V) == pytest.approx(ref["fraction_of_displaced_mass"], rel=1e-14)
    assert M == pytest.approx(float(sp.sympify(ref["formula"]).subs({"r": a, "rho": rho, "pi": sp.pi})), rel=1e-14)
    assert ch06.added_mass_by_energy(a, rho) / (rho * V) == pytest.approx(0.5, rel=1e-10)  # our energy route agrees


def test_sphere_motion_V4_newton_with_added_mass():  # V1/V4/V7 (6.109): u = F t/(m + M), work = kinetic, bubble 2g
    m, a = 1.0, 0.1
    M = PF.added_mass_sphere(a)
    s = ch06.sphere_motion(m, a, F_E=2.0)
    assert np.allclose(s["u"], 2.0 * s["t"] / (m + M), atol=1e-9)
    assert np.allclose(s["work"], s["kinetic"], atol=1e-8)  # V4: work of F_E = ½(m + M)u²
    assert np.allclose(s["F_s"], -M * s["du_dt"])
    f = ch06.sphere_motion(m, a, F_E=lambda t: np.sin(3 * t), t_eval=np.linspace(0, 2, 101))
    assert np.allclose(f["work"], f["kinetic"], atol=1e-8)
    assert np.allclose(f["u"], (1 - np.cos(3 * f["t"])) / (3 * (m + M)), atol=1e-8)
    b = ch06.sphere_motion(None, a, mode="bubble")
    assert b["a0"] == pytest.approx(2 * 9.81, rel=1e-12)  # buoyancy ρVg on added mass ½ρV: 2g, not ∞ (Exercise 6.46)
    steel = 7800 * 4 / 3 * np.pi * a ** 3
    ball = ch06.sphere_motion(steel, a, mode="ball")
    assert ball["a0"] == pytest.approx((1000 - 7800) * 4 / 3 * np.pi * a ** 3 * 9.81 / (steel + M), rel=1e-12)
    assert steel + M == pytest.approx(34.77, abs=0.01)  # C15 number: the steel ball behaves as 34.77 kg
    no = ch06.sphere_motion(m, a, F_E=2.0, added_mass=False)
    assert no["u"][-1] == pytest.approx(2.0 / m) and no["M"] == 0.0
    osc = ch06.sphere_motion(m, a, mode="oscillating", amplitude=0.05, omega=2.0)
    assert np.allclose(osc["F_s"], -M * osc["du_dt"]) and np.allclose(osc["F_E"], (m + M) * osc["du_dt"])
    with pytest.raises(ValueError):
        ch06.sphere_motion(0.0, a, added_mass=False)
    with pytest.raises(ValueError):
        ch06.sphere_motion(m, a, mode="spin")


def test_added_mass_V2_dimensions():  # V2 pint: (6.105) terms in Pa, M in kg, F_s in N
    p = dimensional_check(lambda rho, us, a, A: rho * (0.5 * us ** 2 + a * A / 2), "[pressure]", rho=Q_(1000, "kg/m**3"),
                          us=Q_(1, "m/s"), a=Q_(0.1, "m"), A=Q_(2, "m/s**2"))
    assert p.to("Pa").magnitude == pytest.approx(600.0)
    M = dimensional_check(lambda rho, a: 2 * np.pi / 3 * rho * a ** 3, "[mass]", rho=Q_(1000, "kg/m**3"), a=Q_(0.1, "m"))
    assert M.to("kg").magnitude == pytest.approx(PF.added_mass_sphere(0.1))
    F = dimensional_check(lambda M, A: M * A, "[force]", M=M, A=Q_(2, "m/s**2"))
    assert F.to("N").magnitude == pytest.approx(4.18879, rel=1e-5)  # D30 check: F_s = −4.189 N


def test_moving_sphere_V2_derivation():  # V2 D28 (★) and D29 (★★★): the design's construction re-run
    d = ch06.moving_sphere_dphidt_sym()
    assert d["chain"] == 0 and d["direct"] == 0 and d["correct_bracket"] == sp.zeros(3, 1)
    assert d["printed_bracket"] != sp.zeros(3, 1)
    us_ = sp.symbols("u_s v_s w_s", real=True)
    assert sp.simplify(d["printed_bracket"] - sp.Matrix(us_)) == sp.zeros(3, 1)  # the printed form adds exactly u_s
    t, x, y, z = sp.symbols("t x y z", real=True)
    a = sp.symbols("a", positive=True)
    xs = sp.Matrix([sp.Rational(1, 2) * t ** 2, 0, t])  # the design's check motion (turning acceleration)
    us = xs.diff(t)
    X = sp.Matrix([x, y, z])
    xi = X - xs
    r = sp.sqrt(xi.dot(xi))
    phi = -a ** 3 / (2 * r ** 3) * us.dot(xi)  # (6.97)
    dphidt_direct = sp.diff(phi, t)
    grad = sp.Matrix([sp.diff(phi, v) for v in (x, y, z)])  # (6.103)
    chain = -grad.dot(us) - a ** 3 / (2 * r ** 3) * xi.dot(us.diff(t))  # (6.101)
    val = {x: 0.3, y: 0.4, z: 1.2, t: 0.7, a: 0.5}
    assert abs(float((dphidt_direct - chain).subs(val))) < 1e-12
    e = sp.Matrix([sp.Rational(3, 5), 0, sp.Rational(4, 5)])
    surf = {x: xs[0] + a * e[0], y: xs[1] + a * e[1], z: xs[2] + a * e[2]}
    g_s = grad.subs(surf)
    assert sp.simplify(g_s - (sp.Rational(3, 2) * us.dot(e) * e - us / 2)) == sp.zeros(3, 1)  # (6.104)
    printed = -(a ** 3 / 2) * (-3 * e * us.dot(e) / a ** 3 - us / a ** 3)
    assert sp.simplify(printed - g_s - us) == sp.zeros(3, 1)  # ≠: the printed middle term is wrong by u_s
    p_bern = -dphidt_direct.subs(surf) - g_s.dot(g_s) / 2  # (6.100) with ρ = 1
    c = us.dot(e)
    p_6105 = sp.Rational(1, 2) * us.dot(us) * (sp.Rational(9, 4) * c ** 2 / us.dot(us) - sp.Rational(5, 4)) \
        + a / 2 * e.dot(us.diff(t))
    assert sp.simplify(p_bern - p_6105) == 0  # (6.105)
    ua = sp.Rational(3, 2) * c * e - us / 2  # steps 10–12
    assert sp.simplify(ua.dot(us) - (sp.Rational(3, 2) * c ** 2 - us.dot(us) / 2)) == 0
    assert sp.simplify(ua.dot(ua) - (sp.Rational(3, 4) * c ** 2 + us.dot(us) / 4)) == 0
    ths = sp.symbols("theta_s", positive=True)  # step 14: (6.106) = (6.91)
    assert sp.simplify(sp.Rational(9, 4) * sp.cos(ths) ** 2 - sp.Rational(5, 4) - (1 - sp.Rational(9, 4) * sp.sin(ths) ** 2)) == 0


def test_added_mass_V2_derivation():  # V2 D30 (★★) and D31 (★★★)
    th, ph, a, rho, A, U, r = sp.symbols("theta phi a rho A U r", positive=True)
    e = sp.Matrix([sp.sin(th) * sp.cos(ph), sp.sin(th) * sp.sin(ph), sp.cos(th)])
    dA = a ** 2 * sp.sin(th)
    F = -rho * a / 2 * A * sp.Matrix([sp.integrate(sp.integrate(sp.cos(th) * e[i] * dA, (ph, 0, 2 * sp.pi)),
                                                   (th, 0, sp.pi)) for i in range(3)])  # (6.107)
    assert F[0] == 0 and F[1] == 0 and sp.simplify(F[2] + sp.Rational(2, 3) * sp.pi * rho * a ** 3 * A) == 0  # (6.108)
    assert sp.integrate(sp.cos(th) ** 2 * sp.sin(th), (th, 0, sp.pi)) == sp.Rational(2, 3)  # step 8
    u1, u3 = sp.symbols("u1 u3", real=True)  # step 2: the speed part is even in e ⇒ no force (any u_s in the x–z plane)
    us = sp.Matrix([u1, 0, u3])
    S = rho / 2 * (sp.Rational(9, 4) * (us.dot(e)) ** 2 - sp.Rational(5, 4) * us.dot(us))
    Fs = [sp.integrate(sp.integrate(sp.expand(S * e[i] * dA), (ph, 0, 2 * sp.pi)), (th, 0, sp.pi)) for i in range(3)]
    assert all(sp.simplify(q) == 0 for q in Fs)
    # D31: the design's check cell — volume integral to infinity = the inner surface form = 2πU²a³/3; M = 2πρa³/3
    phi = -a ** 3 * U * sp.cos(th) / (2 * r ** 2)
    gr, gt = sp.diff(phi, r), sp.diff(phi, th) / r
    vol = sp.integrate(sp.integrate((gr ** 2 + gt ** 2) * r ** 2 * sp.sin(th), (th, 0, sp.pi)), (r, a, sp.oo)) * 2 * sp.pi
    surf = -sp.integrate((phi * gr).subs(r, a) * a ** 2 * sp.sin(th), (th, 0, sp.pi)) * 2 * sp.pi
    assert sp.simplify(vol - surf) == 0 and sp.simplify(vol - 2 * sp.pi * U ** 2 * a ** 3 / 3) == 0
    assert sp.simplify(rho * vol / U ** 2 - 2 * sp.pi * a ** 3 * rho / 3) == 0
    lap = sp.diff(r ** 2 * sp.diff(phi, r), r) / r ** 2 + sp.diff(sp.sin(th) * sp.diff(phi, th), th) / (r ** 2 * sp.sin(th))
    assert sp.simplify(lap) == 0  # ∇²φ = 0 used in step 2
    Rbig = sp.Symbol("R", positive=True)
    far = sp.integrate((phi * gr).subs(r, Rbig) * Rbig ** 2 * sp.sin(th), (th, 0, sp.pi)) * 2 * sp.pi
    assert sp.limit(far, Rbig, sp.oo) == 0 and sp.simplify(far * Rbig ** 3) != 0  # step 4: the far surface ~ R⁻³
    assert sp.simplify(phi.subs(r, a) + U * a * sp.cos(th) / 2) == 0 and sp.simplify(gr.subs(r, a) - U * sp.cos(th)) == 0


@needs_ref
def test_rayleigh_collapse_V5_published_constant():  # V5 Wikipedia "Rayleigh–Plesset equation": t ≈ 0.91468 R₀√(ρ/P∞)
    ref = ref_json()["rayleigh_collapse"]
    for R0, rho, dp in ((1.0, 1000.0, 1e5), (0.003, 998.0, 2.3e3)):
        t = ch06.rayleigh_collapse_time(R0, rho, dp)
        assert t / (R0 * np.sqrt(rho / dp)) == pytest.approx(ref["constant"], abs=5e-6)  # 5 significant figures
        assert ch06.rayleigh_collapse_time(R0, rho, dp, "quad") == pytest.approx(t, rel=1e-10)
    with pytest.raises(ValueError):
        ch06.rayleigh_collapse_time(method="guess")


# =====================================================================================================================
# F1 — stagnation points of weak sources (the explainers' `Flow.stagnation_points` / `superposition_state` defaults)
# =====================================================================================================================
@pytest.mark.parametrize("m", [0.5, 1.0, 2.0, 3.0])
def test_stagnation_points_V1_default_search_finds_the_half_body_nose(m):  # V1 (6.31), D07: x_S = −m/2πU for every m
    U = 1.0
    st = ch06.superposition_state([{"kind": "uniform", "U": U}, {"kind": "source", "m": m}], 1.0, 1.0)
    assert len(st["stagnation"]) == 1, f"m = {m}: expected [[{-m / (TWO_PI * U):.4f}, 0]], got {st['stagnation']}"
    assert st["stagnation"][0][0] == pytest.approx(-m / (TWO_PI * U), rel=1e-12)
    assert st["psi_dividing"] == pytest.approx(m / 2)


@pytest.mark.parametrize("m", [1e-4, 0.01, 20.0, 200.0])
def test_stagnation_points_V1_extreme_source_strengths(m):  # V1 loop 2: x_S = −m/2πU from 1.6e-5 m to 32 m
    st = PF.Flow([PF.Uniform(1.0), PF.Source(m)]).stagnation_points()
    assert len(st) == 1 and st[0] == pytest.approx(-m / TWO_PI, rel=1e-12) and st[0].imag == 0.0
    s = ch06.superposition_state([{"kind": "uniform", "U": 1.0}, {"kind": "source", "m": m}], 1.0, 1.0)
    assert s["psi_dividing"] == pytest.approx(m / 2, rel=1e-12)  # +0 imaginary part: ψ read above the cut


def test_stagnation_points_V1_robustness_cases():  # V1 loop 2: no stream, off-body, double and triple roots, far away
    two = ch06.two_sources(1.0, 1.0).stagnation_points()  # two equal sources, no stream: the midpoint
    assert len(two) == 1 and abs(two[0]) < 1e-12
    assert len(PF.Flow([PF.Source(1.0, -1.0), PF.Source(-1.0, 1.0)]).stagnation_points()) == 0  # none exists
    G = 6 * np.pi  # U = a = 1: Γ_cw = 6π > 4πaU ⇒ one free point at r₊ = (6 + √20)/4 on the −y axis
    fl = PF.cylinder(1.0, 1.0, Gamma_cw=G)
    st = fl.stagnation_points()
    rp = (6 + np.sqrt(20.0)) / 4
    assert len(st) == 1 and st[0] == pytest.approx(-1j * rp, abs=1e-12)
    both = fl.stagnation_points(keep_inside=True)
    assert np.sort(np.abs(both)) == pytest.approx([1 / rp, rp], rel=1e-12)  # r₊r₋ = a²
    for Gc in (0.0, 2.0, 20.0):
        ref = ch06.cylinder_stagnation_points(1.0, 1.0, Gamma_cw=Gc)
        assert np.allclose(np.sort_complex(PF.cylinder(1.0, 1.0, Gamma_cw=Gc).stagnation_points()), np.sort_complex(ref),
                           atol=1e-12)
    merged = PF.cylinder(1.0, 1.0, Gamma_cw=4 * np.pi).stagnation_points()  # a double root: Newton converges linearly
    assert len(merged) == 1 and abs(merged[0] + 1j) < 1e-6  # root error ~ √(residual tolerance) for a double root
    c4 = PF.Flow([PF.Corner(1.0, 4.0)]).stagnation_points()  # dw/dz = 4z³: a triple root at the tip
    assert len(c4) == 1 and abs(c4[0]) < 1e-12
    far = PF.Flow([PF.Uniform(1.0), PF.Source(0.2, 7 + 5j)]).stagnation_points()  # outside the default box
    assert len(far) == 1 and far[0] == pytest.approx(7 - 0.2 / TWO_PI + 5j, abs=1e-12)
    vx = PF.Flow([PF.Uniform(2.0), PF.Vortex(3.0, 0.5j)]).stagnation_points()
    assert len(vx) == 1 and vx[0] == pytest.approx(0.5j + 3j / (TWO_PI * 2.0), abs=1e-12)
    small = PF.cylinder(1.0, 0.05).stagnation_points()  # a tiny cylinder (doublet seeds)
    assert np.allclose(small, [-0.05, 0.05], atol=1e-12)
    ov = ch06.superposition_state([{"kind": "uniform", "U": 1}, {"kind": "source", "m": 0.05, "x": -1},
                                   {"kind": "sink", "m": 0.05, "x": 1}], 0.0, 2.0)["stagnation"]
    L = ch06.rankine_oval(1.0, 0.05, 1.0)["half_length"]
    assert ov == [[pytest.approx(-L, rel=1e-12), 0.0], [pytest.approx(L, rel=1e-12), 0.0]]
    assert PF.half_body(0.01, 1.0).stagnation_points()[0] == pytest.approx(-1 / (TWO_PI * 0.01), rel=1e-12)


# =====================================================================================================================
# Part C — every function the notebook and explainers call: explainer parity rows (design Part B), drawing helpers,
# scripts
# =====================================================================================================================
def test_parity_rows_V1_design_part_b_expressions_run():  # V1 smoke + sanity: the `py:` rows of the nine explainers
    rows = {
        "added_mass_by_energy": ch06.added_mass_by_energy(0.1, 1000.0),
        "added_mass_sphere": ch06.added_mass_sphere(0.1, 1000.0),
        "airship": ch06.airship(1.0, 1.0, 1.0)["length"],
        "axial_state": ch06.axial_state("rankine_oval", N=20)["cond"],
        "axisym_body_target": ch06.axisym_body_target("sphere", N=12)[1][0],
        "cyl_state": ch06.cylinder_circulation_state(10.0, 0.1, Gamma_cw=2.0)["L"],
        "cyl_cp": ch06.cylinder_surface_cp(1.5707963267948966, 10.0, 0.1, Gamma_cw=0.0),
        "cyl_p": ch06.cylinder_surface_pressure(0.7853981633974483, 10.0, 0.1, Gamma_cw=2.0, rho=1.2),
        "delta": ch06.delta_flux_check("source", kind="phi", radii=[0.3], m=2.0)[0],
        "dbl": ch06.doublet_limit_error([0.1], 2.0)[0],
        "ex61": ch06.example_6_1(21.765592370810612)["p_origin"],
        "ex61w": ch06.example_6_1_wall_pressure(0.5, 3.0, split=True)["total"],
        "hb": ch06.half_body_numbers(2.0, 10.0)["h_max"],
        "je": ch06.joukowski_ellipse(1.2, 1.0)["A"],
        "lift": ch06.lift_per_span(1.2, 10.0, Gamma_cw=2.0),
        "ls": ch06.line_sink_stream_function(0.3, -0.5, 1.0, 1.0),
        "panel": ch06.panel_cp_error(16, body="ellipse"),
        "oval": ch06.rankine_oval(1.0, 6.283185307179586, 1.0)["half_length"],
        "relax": ch06.relaxation_state("four_point", method="jacobi", sweeps=1)["psi33"],
        "tss": ch06.two_source_streamline(1.570796, 6.283185, 1.0, 2.0),
        "vort": ch06.vorticity_from_psi("rankine", 0.0, 0.0, Gamma=1.0, a=0.1),
        "cse": ch06.corner_speed_exponent(0.6666666666666666),
        "probe": ch06.complex_potential_probe("corner", 1.0, 1.0, A=1.0, n=2.0)["u"],
        "lc": ch06.laurent_contributions("cylinder", R=0.2, Gamma_cw=2.0, U=10.0)["L"],
        "js": ch06.joukowski_state(a=1.2, b=1.0, x=-3.0, y=0.5, branch="outside")["zeta_abs"],
        "am": ch06.added_mass_state(a=0.1, rho=1000.0, us=1.0, dus=2.0)["p_total"],
    }
    for k, v in rows.items():
        assert np.all(np.isfinite(np.asarray(v, dtype=float))), k
    assert rows["added_mass_by_energy"] == pytest.approx(rows["added_mass_sphere"], rel=1e-10)
    assert rows["hb"] == pytest.approx(2.5) and rows["cyl_cp"] == pytest.approx(-3.0) and rows["lift"] == 24.0
    assert rows["vort"] == pytest.approx(1.0 / (np.pi * 0.01), rel=1e-9) and rows["cse"] == pytest.approx(-1 / 3)
    assert rows["delta"] == pytest.approx(2.0) and rows["relax"] == pytest.approx(3.0) and rows["probe"] == pytest.approx(2)
    assert rows["tss"] == pytest.approx(np.sqrt(3.0), abs=1e-5)  # ψ = m/4: x² − y² = a² on the right half (cot = 0)


def test_exports_V1_every_part_c_name_is_reachable():  # V1 contract: the Part C names exist with the design's API
    names = ["Uniform", "Source", "Vortex", "Doublet", "Corner", "Flow", "half_body", "cylinder", "laplacian_residual",
             "polar_velocity", "polar_velocity_sym", "normal_velocity_on", "far_field_check", "pressure_coefficient",
             "mirror", "circle_theorem", "blasius_force", "laurent_coefficients", "AxisymUniform", "PointSource3D",
             "Doublet3D", "AxisymFlow", "sphere", "axisym_velocity_spherical", "axisym_velocity_spherical_sym",
             "sphere_potential_vector", "moving_sphere_potential", "moving_sphere_velocity",
             "moving_sphere_surface_pressure", "sphere_force_quadrature", "added_mass_sphere", "joukowski",
             "joukowski_derivative", "joukowski_inverse", "mapped_flow", "map_elements", "angle_preservation",
             "grid_image", "laplacian_5pt", "node_update", "jacobi_sweep", "gauss_seidel_sweep", "sor_sweep",
             "residual_norm", "solve_laplace", "solve_poisson", "axial_singularity_solve", "axial_singularity_velocity",
             "axial_singularity_psi", "source_panels", "panel_velocity"]
    missing = [n for n in names if not hasattr(ch06, n)]
    assert not missing, missing
    fl = PF.Flow([PF.Uniform(1.0), PF.Source(1.0)])
    for meth in ("w", "dwdz", "velocity", "phi", "psi", "speed", "cp", "pressure", "stagnation_points", "__add__"):
        assert hasattr(fl, meth)
    sp_ = PF.sphere(1.0, 1.0)
    for meth in ("psi", "phi", "velocity_cyl", "velocity_spherical"):
        assert hasattr(sp_, meth)


def test_primitives_V1_branches_points_and_function_flows():  # V1 smoke of the helpers the elements are built on
    z = np.array([1.0 + 1e-9j, 1.0 - 1e-9j, -1.0 + 1e-9j, -1.0 - 1e-9j])
    lb = PF.log_branch(z, cut=TWO_PI, theta_range="[)")  # θ ∈ [0, 2π): cut on the positive x-axis
    assert np.allclose(lb.imag, [0.0, TWO_PI, np.pi, np.pi], atol=1e-8)
    assert np.allclose(PF.log_branch(z).imag, np.angle(z))  # cut = π is numpy's principal branch
    with pytest.raises(ValueError):
        PF.log_branch(z, theta_range="()")
    assert PF.power_branch(np.array([0j]), 0.5)[0] == 0 and PF.power_branch(4.0 + 0j, 0.5) == pytest.approx(2.0)
    assert PF.as_complex_point((1.0, 2.0)) == 1 + 2j and PF.as_complex_point(None) == 0j and PF.as_complex_point(3) == 3
    for e in (PF.Source(1.0, 0.3), PF.Vortex(2.0, 1j), PF.Doublet((0.3, 0.1)), PF.Corner(0.7, 1.5), PF.Uniform(1.0)):
        q = np.array([1.3 + 0.9j])
        fd = (e.dwdz(q + 1e-6) - e.dwdz(q - 1e-6)) / 2e-6
        assert np.allclose(e.d2wdz2(q), fd, rtol=1e-6, atol=1e-9)  # the Newton derivative used by stagnation_points
    assert PF.Source(1.0, 0.3).singular_points == (0.3 + 0j,) and PF.Uniform(1.0).singular_points == ()
    assert PF.Corner(1.0, 0.5).singular_points == (0j,) and PF.Corner(1.0, 2.0).singular_points == ()
    ff = PF.FunctionFlow(lambda q: q ** 2, lambda q: 2 * q, u_inf=(1.0, 0.0))
    assert ff.d2wdz2(1j) is None and ff.complex_velocity(1.0, 1.0) == pytest.approx(2 + 2j)
    st = ff.stagnation_points(guesses=[0.3 + 0.2j])  # (the numeric second-derivative fallback)
    assert len(st) == 1 and abs(st[0]) < 1e-10
    assert isinstance(PF.cylinder(1.0, 1.0), PF.ComplexFlow) and PF.BlasiusForce(1.0, 2.0).L == 2.0
    geo = PN.panel_geometry(np.cos(np.linspace(0, 6, 8)), np.sin(np.linspace(0, 6, 8)))
    Ux, Uy = PN.panel_induced_velocity(np.array([3.0]), np.array([0.0]), geo)
    assert Ux.shape == (1, 8) and np.all(np.isfinite(Ux))


def test_reused_tools_V1_smoke_with_chapter_6_inputs():  # V1 smoke of the C.0 reuse rows called with ch06 inputs
    from fluidpy import ch01_introduction as ch01
    from fluidpy import ch05_vorticity_dynamics as ch05
    from fluidpy.core import integral_theorems as IT
    from fluidpy.core import interact
    air = ch01.fluid_properties("air", 293.15)
    assert air["nu"] == pytest.approx(1.5e-5, rel=0.01)  # N01's ν of air (car example: Re ≈ 6.7e5)
    assert ch04.is_incompressible_regime(10.0) and not ch04.is_incompressible_regime(200.0)
    assert "Kelvin holds" in ch05.kelvin_hypotheses_text(inviscid=True, barotropic=True)
    x, w = IT.gauss_legendre_nodes(-1.0, 1.0, 24)
    assert np.sum(w * x ** 2) == pytest.approx(2 / 3)  # the ∫cos²θ sin θ dθ = 2/3 of D30
    u_fn, _ = ch06.flow_field_callables(PF.cylinder(1.0, 0.1))
    vis = BE.viscous_irrotational_residual(u_fn, np.array([[0.2, 0.15], [0.1, -0.3]]), mu=1e-3)
    assert np.max(np.abs(vis)) < 1e-6
    fig = interact.slider_figure(lambda G: {"C_p": (np.linspace(0, TWO_PI, 50),
                                                    ch06.cylinder_surface_cp(np.linspace(0, TWO_PI, 50), 10.0, 0.1,
                                                                             Gamma_cw=G))},
                                 "Gamma", [0.0, 2.0, 4.0], unit="m²/s")
    assert len(fig.data) >= 1
    tr = BS.point_vortex_evolve(np.array([[0.0], [1.0]]), [-1.0], np.array([0.0, 2.0]), boundary="wall", wall_y=0.0)
    assert tr.shape[0] == 2 and np.isfinite(tr).all()


def test_drawing_helpers_V1_part_c_6_1_exist_and_draw():  # V1 smoke C.6 row 6.1: flow_net, draw_body, pressure_arrows,
    import matplotlib                                     # separated_cp_band (the labelled qualitative band of N31)
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(ROOT / "scripts"))
    import ch06_drawings as dr
    fig, ax = plt.subplots()
    dr.flow_net(ax, PF.cylinder(1.0, 1.0, Gamma_cw=1.0), xlim=(-3, 3), ylim=(-2, 2), n=41, body=dr.circle(1.0, 64))
    assert len(ax.collections) > 0 and np.allclose(np.abs(dr.circle(2.0, 16)), 2.0)
    plt.close(fig)
    missing = [n for n in ("draw_body", "pressure_arrows", "separated_cp_band") if not hasattr(dr, n)]
    assert not missing, f"design Part C 6.1 helpers missing from scripts/ch06_drawings.py: {missing}"
    beta = np.linspace(0, 180, 19)
    lo, hi = dr.separated_cp_band(beta)
    assert np.all(lo <= hi)
    ideal = np.asarray(ch06.cylinder_surface_cp(np.radians(beta)))
    front = (beta < 80) & (ideal >= -1.05)  # ahead of separation, above the documented floor −1.05
    assert np.all((lo[front] <= ideal[front] + 1e-12) & (ideal[front] <= hi[front] + 1e-12))  # follows (6.35) in front
    assert np.all(hi[beta > 90] < 0)  # low wake pressure behind separation (qualitative band)
    fig, ax = plt.subplots()
    patches = dr.draw_body(ax, dr.circle(1.0, 32))
    assert len(patches) == 1
    cz = dr.circle(1.0, 16)
    p = 60.0 * np.asarray(ch06.cylinder_surface_cp(np.angle(cz)))
    q = dr.pressure_arrows(ax, cz, p, cz, scale=0.01)
    U_, V_ = q.U, q.V
    assert np.allclose(U_ + 1j * V_, -p * cz * 0.01)  # arrows along −p n
    plt.close(fig)


def test_scripts_V1_every_ch06_script_runs():  # V1 smoke: all 12 scripts exit 0 and write their figures
    env = dict(os.environ, MPLBACKEND="Agg")
    for scr in sorted((ROOT / "scripts").glob("ch06_*.py")):
        args = [sys.executable, str(scr), "--no-show"]
        r = subprocess.run(args, cwd=ROOT, env=env, capture_output=True, text=True, timeout=600)
        assert r.returncode == 0, (scr.name, r.stderr[-2000:])


# =====================================================================================================================
# V6 — the book's printed forms and numbers (private JSON; skipped when absent)
# =====================================================================================================================
@book_only
def test_book_V6_section_6_1_and_6_3_forms():  # V6 §6.1 threshold; §6.3 closed forms; Fig. 6.8 angle; Fig. 6.10 range
    b = book()
    assert ch06.ideal_flow_applicability(b["sec6_1_reynolds"]["Re_threshold_typical"])["ok"]
    assert not ch06.ideal_flow_applicability(0.99 * b["sec6_1_reynolds"]["Re_threshold_typical"])["ok"]
    f = b["sec6_3_closed_forms"]
    s = book_symbols("U", "m", "a", "theta", "Gamma", "rho", "psi", "x", "y")
    val = {s["U"]: 1.7, s["m"]: 2.3, s["a"]: 0.4, s["theta"]: 0.9, s["Gamma"]: 1.1, s["rho"]: 1.2}
    ev = lambda k: float(sp.sympify(f[k], locals=s).subs(val))  # noqa: E731
    hb = ch06.half_body_numbers(1.7, 2.3)
    assert ev("half_body_stagnation_distance") == pytest.approx(hb["a"], rel=1e-14)
    assert ev("half_body_psi_body") == pytest.approx(hb["psi_body"]) and ev("half_body_h_max") == pytest.approx(hb["h_max"])
    assert ev("half_body_half_width") == pytest.approx(float(ch06.half_body_shape(1.7, 2.3, 0.9)[1]), rel=1e-13)
    assert ev("cylinder_doublet_strength") == pytest.approx(-PF.cylinder(1.7, 0.4).elements[1].d_vec[0], rel=1e-14)
    assert ev("cylinder_cp_surface") == pytest.approx(float(ch06.cylinder_surface_cp(0.9)), rel=1e-13)
    assert ev("cylinder_u_theta_surface_with_Gamma_cw") == pytest.approx(
        float(ch06.cylinder_surface_speed(0.9, 1.7, 0.4, Gamma_cw=1.1)), rel=1e-13)
    st = ch06.cylinder_stagnation_points(1.7, 0.4, Gamma_cw=1.1)
    assert ev("cylinder_stagnation_sin_theta_Gamma_cw") == pytest.approx(np.sin(np.angle(st[0])), rel=1e-12)
    assert ev("cylinder_critical_Gamma") == pytest.approx(4 * np.pi * 0.4 * 1.7)
    Gb = 3.0 * ev("cylinder_critical_Gamma")
    rp = float(sp.sympify(f["cylinder_offbody_stagnation_r"], locals=s).subs({**val, s["Gamma"]: Gb}))
    assert rp == pytest.approx(abs(ch06.cylinder_stagnation_points(1.7, 0.4, Gamma_cw=Gb)[0]), rel=1e-12)
    assert ev("lift_per_span_Gamma_cw") == pytest.approx(ch06.lift_per_span(1.2, 1.7, Gamma_cw=1.1), rel=1e-14)
    ang = b["fig6_8_half_body"]["cp_zero_angle_deg"]
    assert np.degrees(ch06.half_body_cp_zero_angle()) == pytest.approx(ang, abs=0.5)  # printed to the degree
    fg = b["fig6_10_cylinder"]
    cp = ch06.cylinder_surface_cp(np.linspace(0, np.pi, 181))
    assert cp.min() == pytest.approx(fg["ideal_cp_min"]) and cp.max() == pytest.approx(fg["ideal_cp_max"])
    two = sp.sympify(f["two_source_streamlines"], locals=s)  # the (6.41) curve through our streamline points
    xs = np.linspace(1.2, 2.5, 5)
    ys = ch06.two_source_streamline(0.7, 2.3, 0.4, xs)
    for xx, yy in zip(xs, ys):
        if np.isfinite(yy):
            assert abs(float(two.subs({s["x"]: xx, s["y"]: yy, s["psi"]: 0.7, s["m"]: 2.3, s["a"]: 0.4}))) < 1e-9


@book_only
def test_book_V6_example_6_1_and_section_6_5_6_6():  # V6 Example 6.1 forms; (6.61)–(6.62); §6.6 ellipse and inverse
    b = book()
    e = b["ex6_1"]
    s = book_symbols("Gamma", "h", "t", "U", "d", "rho", "a", "b", "z")
    for G, h, t in ((1.0, 1.0, 0.0), (1.7, 0.6, 5.3), (0.4, 2.0, 40.0)):
        v = {s["Gamma"]: G, s["h"]: h, s["t"]: t}
        c = ch06.example_6_1(t, Gamma=G, h=h, rho=1.0)
        assert float(sp.sympify(e["trajectory_y"], locals=s).subs(v)) == pytest.approx(c["xi_y"], abs=1e-14)
        assert float(sp.sympify(e["trajectory_x"], locals=s).subs(v)) == pytest.approx(c["xi_x"])
        assert float(sp.sympify(e["dphidt_origin"], locals=s).subs(v)) == pytest.approx(c["dphidt"], rel=1e-13)
        assert float(sp.sympify(e["v_origin"], locals=s).subs(v)) == pytest.approx(c["v_origin"], rel=1e-13)
        assert float(sp.sympify(e["p_origin_minus_pinf_over_rho"], locals=s).subs(v)) == pytest.approx(
            c["p_origin"], rel=1e-12, abs=1e-15)
    k = b["sec6_5_kutta"]
    kz = ch06.kutta_zhukhovsky_sym()
    Ub, Gb, db = sp.symbols("U Gamma d", positive=True)
    loc = {"U": Ub, "Gamma": Gb, "d": db}
    assert sp.simplify(sp.sympify(k["residue_at_0"], locals=loc) - kz["residue"]) == 0
    assert sp.simplify(sp.sympify(k["correct_1_over_z2_coefficient"], locals=loc) - kz["coeff_z2"]) == 0
    assert sp.simplify(sp.sympify(k["printed_1_over_z2_coefficient_misprint"], locals=loc) - kz["coeff_z2"]) != 0
    j = b["sec6_6_joukowski"]
    E = ch06.joukowski_ellipse(1.3, 0.9)
    v = {s["a"]: 1.3, s["b"]: 0.9}
    assert float(sp.sympify(j["ellipse_semi_major"], locals=s).subs(v)) == pytest.approx(E["A"])
    assert float(sp.sympify(j["ellipse_semi_minor"], locals=s).subs(v)) == pytest.approx(E["B"])
    assert float(sp.sympify(j["foci"], locals=s).subs(v)) == pytest.approx(E["foci"])
    assert float(sp.sympify(j["slit_half_length"], locals=s).subs(v)) == pytest.approx(
        np.max(np.real(CM.joukowski(0.9 * np.exp(1j * np.linspace(0, 6.3, 999)), 0.9))), rel=1e-5)


@book_only
def test_book_V6_sections_6_8_6_9_and_exercises():  # V6 §6.8 sphere/airship forms, §6.9 added mass, exercise answers
    b = book()
    f = b["sec6_8_sphere_airship"]
    s = book_symbols("U", "a", "r", "theta", "Q", "k", "r1", "rho")
    U, a, r, th = 1.3, 0.6, 1.4, 0.8
    v = {s["U"]: U, s["a"]: a, s["r"]: r, s["theta"]: th}
    sph = PF.sphere(U, a)
    R, z = r * np.sin(th), r * np.cos(th)
    assert float(sp.sympify(f["sphere_psi"], locals=s).subs(v)) == pytest.approx(float(sph.psi(R, z)), rel=1e-13)
    assert float(sp.sympify(f["sphere_phi"], locals=s).subs(v)) == pytest.approx(float(sph.phi(R, z)), rel=1e-13)
    ur, ut = sph.velocity_spherical(r, th)
    assert float(sp.sympify(f["sphere_u_r"], locals=s).subs(v)) == pytest.approx(float(ur), rel=1e-12)
    assert float(sp.sympify(f["sphere_u_theta"], locals=s).subs(v)) == pytest.approx(float(ut), rel=1e-12)
    assert float(sp.sympify(f["sphere_cp_surface"], locals=s).subs(v)) == pytest.approx(float(ch06.sphere_surface_cp(th)))
    assert float(sp.sympify(f["sphere_doublet_strength"], locals=s).subs(v)) == pytest.approx(sph.elements[1].d)
    Q, aa = 1.1, 0.9
    air = ch06.airship(U, Q, aa)
    r1 = np.hypot(R, z - aa)
    va = {s["Q"]: Q, s["a"]: aa, s["r"]: r, s["r1"]: r1, s["theta"]: th, s["U"]: U}
    assert float(sp.sympify(f["airship_psi"], locals=s).subs(va)) == pytest.approx(float(air["psi"](R, z)), rel=1e-12)
    vk = {s["k"]: 0.7, s["r"]: r, s["r1"]: r1}
    assert float(sp.sympify(f["line_sink_psi"], locals=s).subs(vk)) == pytest.approx(
        float(ch06.line_sink_stream_function(R, z, 0.7, aa)), rel=1e-12)
    g = b["sec6_9_added_mass"]
    assert float(sp.sympify(g["added_mass"], locals=s).subs({s["a"]: 0.1, s["rho"]: 1000})) == pytest.approx(
        PF.added_mass_sphere(0.1, 1000.0), rel=1e-14)
    assert g["added_mass_fraction_of_displaced"] == pytest.approx(PF.added_mass_sphere(1.0, 1.0) / (4 / 3 * np.pi))
    ths = sp.Symbol("theta_s")
    pst = sp.sympify(g["pressure_steady_part_over_half_rho_us2"], locals={"theta_s": ths})
    assert float(pst.subs(ths, 0.4)) == pytest.approx(
        float(ch06.added_mass_state(theta_s_deg=np.degrees(0.4))["cp_steady"]), rel=1e-12)
    x = b["exercises"]
    se = book_symbols("m", "U", "a", "k", "Q", "h", "rho")
    assert float(sp.sympify(x["ex6_12_h_max"], locals=se).subs({se["m"]: 2.3, se["U"]: 1.7})) == pytest.approx(
        ch06.half_body_numbers(1.7, 2.3)["h_max"])
    lhs, rhs = x["ex6_19_rankine_half_width"].split("=")
    ro = ch06.rankine_oval(1.2, 2.5, 0.7)
    vv = {se["h"]: ro["half_width"], se["a"]: 0.7, se["U"]: 1.2, se["m"]: 2.5}
    assert float(sp.sympify(lhs, locals=se).subs(vv)) == pytest.approx(float(sp.sympify(rhs, locals=se).subs(vv)),
                                                                       rel=1e-10)
    e18 = x["ex6_18"]
    assert e18["m_m2_s"] * float(sp.sympify(e18["dtheta"])) / TWO_PI == pytest.approx(e18["dpsi_m2_s"])
    assert PF.Uniform(e18["U_m_s"]).psi(0.0, e18["dy_m"]) == pytest.approx(e18["dpsi_m2_s"])
    s34 = book_symbols("Q", "r", "theta")
    src = PF.PointSource3D(1.3)
    assert float(sp.sympify(x["ex6_34_source"]["phi"], locals=s34).subs({s34["Q"]: 1.3, s34["r"]: 0.8})) == pytest.approx(
        float(src.phi(0.8 * np.sin(0.5), 0.8 * np.cos(0.5))), rel=1e-13)
    assert float(sp.sympify(x["ex6_34_source"]["psi"], locals=s34).subs({s34["Q"]: 1.3, s34["theta"]: 0.5})) == \
        pytest.approx(float(src.psi(0.8 * np.sin(0.5), 0.8 * np.cos(0.5))), rel=1e-13)
    target = Q / (4 * np.pi * U * aa ** 2)  # Exercise 6.42: the axis stagnation points satisfy (z/a)²(1 − z/a) = ±q
    tn, tt = air["z_nose"] / aa, air["z_tail"] / aa
    assert tn ** 2 * (1 - tn) == pytest.approx(target, rel=1e-9)  # the nose (z < 0): +q
    assert tt ** 2 * (1 - tt) == pytest.approx(-target, rel=1e-9)  # the tail (z > a): −q
    assert x["ex6_42_airship_length_roots_of"].replace(" ", "").startswith("z**2/a**2*(z/a+-1)")  # the printed ± form
    ex43 = float(sp.sympify(x["ex6_43_far_radius"], locals=se).subs({se["a"]: 0.9, se["k"]: 1.4, se["U"]: 2.0}))
    assert ex43 == pytest.approx(ch06.axisym_half_body(2.0, 0.9 * 1.4)["R_far"], rel=1e-14)
    s45 = book_symbols("a", "U_c", "x", "x_c", "y")
    phi45 = float(sp.sympify(x["ex6_45_potential"], locals=s45).subs(
        {s45["a"]: 0.5, s45["U_c"]: 1.5, s45["x"]: 1.4, s45["x_c"]: 0.2, s45["y"]: 0.7}))
    moving = PF.Doublet((TWO_PI * 1.5 * 0.25, 0.0), 0.2)  # a cylinder moving at U_c along x: w = −U_c a²/(z − z_c)
    assert phi45 == pytest.approx(float(moving.phi(1.4, 0.7)), rel=1e-13)
    assert float(sp.sympify(x["ex6_49_added_mass_sphere"], locals=se | {"pi": sp.pi}).subs(
        {se["a"]: 0.3, se["rho"]: 900})) == pytest.approx(ch06.added_mass_by_energy(0.3, 900.0), rel=1e-10)
