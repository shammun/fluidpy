"""Verification suite for Chapter 7 — Gravity Waves (Kundu, Cohen & Dowling 5e, §§7.1–7.8).

Evidence levels (``verify-implementation`` skill): V1 analytic · V2 symbolic (sympy / dimensions) · V3 convergence ·
V4 conservation / invariant · V5 published benchmark (``reference/ch07/``, see SOURCES.md) · V6 book value (private,
``tests/book_values_ch07.json``, skipped when absent) · V7 limits / symmetry / invariance. Every test name is
``test_<concept>_<level>_<what>``; the comment on the ``def`` line repeats the level.

A items (CORE, ≥ 2 independent levels): C01 the sinusoid and its vocabulary (7.1)–(7.9) · C02 the linearised
free-surface problem (7.11)–(7.21) · C03 the dispersion relation (7.28) · C04 phase speed and depth limits (7.29),
(7.45)–(7.52) · C05 particle orbits (7.32)–(7.37) · C06 wave energy and flux (7.38)–(7.44) · C07 capillary–gravity waves
(7.53)–(7.60) · C08 standing waves and seiches (7.61)–(7.65) · C09 group velocity (7.66)–(7.71) · C10 kinematic wave
theory and rays (7.72)–(7.79) · C11 hydraulic jump, KdV and the solitary wave (7.80)–(7.81), (7.87)–(7.88) · C12 Stokes
waves and Stokes drift (7.82)–(7.86) · C13 interfacial waves (7.89)–(7.96) · C14 two-layer modes (7.97)–(7.119) · C15
internal-wave dispersion (7.120)–(7.139) · C16 c ⟂ c_g and F = c_g E (7.140)–(7.159). Derivations: every ★★ and ★★★ D row
(D02–D06, D10, D12, D14–D17, D20–D30, D32–D34, D36, D37) re-derived with sympy in ``test_*_V2_derivation``, the ★★★
D20, D29, D30, D33 step by step through the design's Part F lines; the ★ rows D01, D07, D19 of CORE items too.

Pinned conventions and slips (each with a discriminating assertion): z up, still surface z = 0; η is the elevation (the
level-set function is z − η, R1); ψ with u = ∂ψ/∂z; the (7.66) envelope uses ½Δω t (printed x fails); (7.105) with
e^{i(kx − ωt)} (printed kz fails); internal waves sign-safe in k (printed (7.145) fails for k < 0); g′ with ρ₂ below
(7.117) vs ch04's ρ₁; interfacial E_p = ¼Δρga² from direct integration (the printed middle form gives ⅛).

Run: ``.venv/Scripts/python.exe -m pytest tests/test_ch07.py -q -p no:cacheprovider``.
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
import scipy.sparse as spm
import scipy.sparse.linalg as spla
import sympy as sp
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.special import ellipe, ellipk

from fluidpy import ch01_introduction as ch01
from fluidpy import ch04_conservation_laws as ch04
from fluidpy import ch05_vorticity_dynamics as ch05
from fluidpy import ch07_gravity_waves as ch07
from fluidpy.core import bernoulli as BE
from fluidpy.core import interfaces as IF
from fluidpy.core import kinematics as KN
from fluidpy.core import similarity as SIM
from fluidpy.core import stratification as ST
from fluidpy.core import streamfunction as SF
from fluidpy.core import waves as W
from fluidpy.core.units import Q_, dimensional_check
from tools.convergence import observed_order, pairwise_orders

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference" / "ch07"
BOOK = Path(__file__).resolve().parent / "book_values_ch07.json"
needs_ref = pytest.mark.skipif(not (REF / "benchmarks.json").exists(), reason="run reference/ch07/make_refs.py")
book_only = pytest.mark.skipif(not BOOK.exists(), reason="book values are private; see CLAUDE.md rule 9")

RNG = np.random.default_rng(7)
ORDER_TOL = 0.15  # design order ± this (verify-implementation default)
TWO_PI = 2.0 * np.pi
G = 9.81  # the book's g (ch07 default G_BOOK); passed explicitly


def book():
    return json.loads(BOOK.read_text(encoding="utf-8"))


def ref_json():
    return json.loads((REF / "benchmarks.json").read_text(encoding="utf-8"))


def maxrel(num, ref, floor=0.0):
    num, ref = np.asarray(num, dtype=complex), np.asarray(ref, dtype=complex)
    return float(np.max(np.abs(num - ref)) / max(np.max(np.abs(ref)), floor, 1e-300))


def sym_zero(expr, n: int = 6, tol: float = 1e-10, lo: float = 0.3, hi: float = 1.7) -> bool:
    """True if a sympy expression is identically zero: simplify (after rewriting hyperbolics as exponentials); for
    stubborn forms evaluate at ``n`` random points (all free symbols positive in [lo, hi])."""
    e = sp.simplify(sp.expand(sp.sympify(expr).rewrite(sp.exp)))
    if e == 0:
        return True
    syms = sorted(e.free_symbols, key=lambda s: s.name)
    f = sp.lambdify(syms, e, "numpy")
    for _ in range(n):
        v = complex(f(*RNG.uniform(lo, hi, len(syms))))
        if abs(v) > tol:
            return False
    return True


# =====================================================================================================================
# C01 — the sinusoid and its vocabulary (7.1)–(7.9): N02–N11, D01
# =====================================================================================================================
def test_sinusoid_V1_crests_ride_at_omega_over_k():  # V1 (7.2)–(7.4): every crest moves at ω/k; (7.61) the other way
    k, om, a = 0.0628, 0.785, 1.3
    t = np.linspace(0.0, 20.0, 41)[:, None]
    n = np.arange(-3, 4)[None, :]
    xc = W.crest_positions(t, k, om, n)
    assert np.max(np.abs(W.sinusoid(xc, t, a, k, om) - a)) < 1e-12 * a  # η = a on every crest (7.3)
    assert np.max(np.abs(W.sinusoid(xc + np.pi / k, t, a, k, om) + a)) < 1e-12 * a  # troughs half a λ behind
    v = np.diff(xc, axis=0) / np.diff(t, axis=0)
    assert maxrel(v, np.full_like(v, om / k)) < 1e-12  # (7.4) c = ω/k, the same for every n
    xl = -om / k * t + TWO_PI * n / k  # left-going crests of a cos(kx + ωt) (7.61)
    assert np.max(np.abs(W.sinusoid(xl, t, a, k, om, direction=-1) - a)) < 1e-12 * a
    assert np.max(np.abs(W.sinusoid(xl[1:], t[1:], a, k, om, direction=+1) - a)) > 0.1 * a  # wrong direction fails
    lam, c = TWO_PI / k, om / k  # (7.1): a cos[2π(x − ct)/λ] is the same wave
    x = RNG.uniform(-300, 300, 50)
    assert maxrel(W.sinusoid(x, 3.7, a, k, om), a * np.cos(TWO_PI / lam * (x - c * 3.7))) < 1e-12
    with pytest.raises(ValueError):
        W.sinusoid(0.0, 0.0, direction=2)


def test_phase_speed_V2_derivation_crest_condition():  # V2 D01 (★): (7.3) → x_crest(t) → c = ω/k = λν (7.4)
    x, t = sp.symbols("x t", real=True)
    n = sp.Symbol("n", integer=True)  # the crest label
    k, om, lam, nu = sp.symbols("k omega lambda nu", positive=True)
    xc = sp.solve(sp.Eq(k * x - om * t, 2 * n * sp.pi), x)[0]  # step 3
    assert sp.simplify(xc - (om / k * t + 2 * n * sp.pi / k)) == 0
    assert sp.simplify(sp.diff(xc, t) - om / k) == 0  # step 4–5: n drops out
    assert sp.simplify((om / k).subs({om: 2 * sp.pi * nu, k: 2 * sp.pi / lam}) - lam * nu) == 0  # c = λν
    assert sp.simplify(sp.cos(k * xc - om * t) - 1) == 0  # the crest condition itself (steps 1–2)


def test_wave_parameters_V1_round_trips_and_inconsistent_input():  # V1 (7.1), (7.4), (7.7): N03
    p = W.wave_parameters(lam=100.0, T=8.0)  # the curation's worked number (N03)
    assert p["c"] == pytest.approx(12.5, rel=1e-14) and p["k"] == pytest.approx(TWO_PI / 100.0, rel=1e-14)
    assert p["omega"] == pytest.approx(TWO_PI / 8.0, rel=1e-14) and p["nu"] == pytest.approx(0.125, rel=1e-14)
    for kw in ({"k": p["k"], "omega": p["omega"]}, {"lam": p["lam"], "nu": p["nu"]}, {"lam": p["lam"], "c": p["c"]},
               {"k": p["k"], "c": p["c"]}, {"omega": p["omega"], "c": p["c"]}, {"T": p["T"], "c": p["c"]},
               {"k": p["k"], "T": p["T"], "c": p["c"]}):
        q = W.wave_parameters(**kw)
        assert all(q[key] == pytest.approx(p[key], rel=1e-12) for key in p)
    v = W.wave_parameters(k=[1.0, 1.0], omega=1.0)  # (7.6)–(7.7): λ = 2π/K along K
    assert v["k"] == pytest.approx(np.sqrt(2.0), rel=1e-15) and v["lam"] == pytest.approx(TWO_PI / np.sqrt(2), rel=1e-15)
    with pytest.raises(ValueError):
        W.wave_parameters(k=1.0)
    with pytest.raises(ValueError):
        W.wave_parameters(omega=1.0)
    with pytest.raises(ValueError):
        W.wave_parameters(lam=100.0, T=8.0, c=13.0)  # inconsistent third value


def test_plane_wave_V1_crest_spacing_and_trace_velocities():  # V1 (7.5)–(7.8), N06–N10, Fig. 7.1
    K = np.array([1.0, 1.0])
    om, a = 2.0, 0.5
    Km = np.linalg.norm(K)
    eK = K / Km
    s = TWO_PI * np.arange(-4, 5) / Km
    assert np.max(np.abs(W.plane_wave(s[:, None] * eK, K, om, a) - a)) < 1e-12  # crests 2π/K apart along e_K (7.7)
    xs = np.stack([TWO_PI * np.arange(-4, 5) / K[0], np.zeros(9)], axis=-1)
    assert np.max(np.abs(W.plane_wave(xs, K, om, a) - a)) < 1e-12  # … and 2π/k = 6.28 m apart along x (N08)
    assert TWO_PI / Km == pytest.approx(4.4429, abs=1e-4)
    X = RNG.uniform(-5, 5, (30, 2))
    assert np.max(np.abs(W.plane_wave(X, K, om, a, t=0.3) - W.plane_wave(X.T, K, om, a, t=0.3))) < 1e-15  # both layouts
    c = W.phase_velocity_vector(K, om)
    assert np.linalg.norm(c) == pytest.approx(om / Km, rel=1e-14) and abs(c[0] * K[1] - c[1] * K[0]) < 1e-14  # (7.8) c ∥ K
    tr = W.trace_velocities(K, om)
    assert tr == pytest.approx((om / K[0], om / K[1]), rel=1e-15)
    for _ in range(50):  # random 3-D K: every trace speed ≥ c, and 1/c² = Σ 1/c_i² — they are not components
        K3 = RNG.normal(size=3)
        c3 = np.linalg.norm(W.phase_velocity_vector(K3, 1.3))
        tr3 = np.abs(W.trace_velocities(K3, 1.3))
        assert np.all(tr3 >= c3 * (1 - 1e-12))
        assert np.sum(1.0 / tr3 ** 2) == pytest.approx(1.0 / c3 ** 2, rel=1e-12)
        assert np.linalg.norm(tr3) > c3  # the wrong reading "c = (ω/k, ω/l, ω/m)" overstates the speed
    assert W.trace_velocities([2.0, 0.0], 1.0)[1] == np.inf


def test_plane_wave_V7_rotation_invariance():  # V7: η(Rx; RK) = η(x; K); |c| invariant, c rotates with K
    K = np.array([0.7, -1.9])
    X = RNG.uniform(-4, 4, (40, 2))
    for ang in (0.3, 1.1, 2.9):
        R = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
        assert np.max(np.abs(W.plane_wave(X @ R.T, R @ K, 1.7, 1.0, t=0.4) - W.plane_wave(X, K, 1.7, 1.0, t=0.4))) < 1e-13
        assert np.max(np.abs(W.phase_velocity_vector(R @ K, 1.7) - R @ W.phase_velocity_vector(K, 1.7))) < 1e-14


def test_doppler_V1_probe_frequency_of_a_translated_pattern():  # V1 (7.9), N11: a probe sees ω₀ = ω + U·K
    from scipy.signal import hilbert

    k, om, U = 0.35, 0.9, 1.4
    t = np.linspace(0.0, 400.0, 40001)
    probe = W.sinusoid(2.0 - U * t, t, 1.0, k, om)  # the wave written in the frame moving with the current
    ph = np.unwrap(np.angle(hilbert(probe)))
    mid = slice(4000, 36000)
    slope = np.polyfit(t[mid], ph[mid], 1)[0]  # the analytic signal of cos(c − ω₀t) turns at +ω₀
    assert slope == pytest.approx(om + U * k, rel=1e-6)
    assert float(W.doppler_frequency(om, U, k)) == pytest.approx(om + U * k, rel=1e-15)
    assert float(W.doppler_frequency(0.0, [2.0, 0.5], [0.3, 1.0])) == pytest.approx(2.0 * 0.3 + 0.5, rel=1e-15)
    assert float(W.doppler_frequency(0.0, U, k)) == pytest.approx(U * k)  # frozen pattern (ω = 0)


def test_linear_evolve_V1_single_modes_and_superposition():  # V1 N02 (Fourier superposition; Exercise 7.3)
    L, N, H = 64.0, 256, 3.0
    x = np.linspace(0, L, N, endpoint=False)
    k1, k2 = TWO_PI * 4 / L, TWO_PI * 11 / L
    om1, om2 = float(W.omega_gravity(k1, H, G)), float(W.omega_gravity(k2, H, G))
    t = np.array([0.0, 1.3, 7.9])
    right = W.linear_evolve(0.2 * np.cos(k1 * x), x, t, direction=+1, H=H, g=G)
    assert np.max(np.abs(right - 0.2 * np.cos(k1 * x[None] - om1 * t[:, None]))) < 1e-13
    rest = W.linear_evolve(0.2 * np.cos(k1 * x), x, t, H=H, g=G)  # released from rest → standing wave
    assert np.max(np.abs(rest - 0.2 * np.cos(k1 * x[None]) * np.cos(om1 * t[:, None]))) < 1e-13
    with_rate = W.linear_evolve(0.2 * np.cos(k1 * x), x, t, eta_t0=0.2 * om1 * np.sin(k1 * x), H=H, g=G)
    assert np.max(np.abs(with_rate - right)) < 1e-13  # the right initial rate picks the right-going wave
    two = W.linear_evolve(0.2 * np.cos(k1 * x) + 0.05 * np.sin(k2 * x), x, t, direction=+1, H=H, g=G)
    exact = 0.2 * np.cos(k1 * x[None] - om1 * t[:, None]) + 0.05 * np.sin(k2 * x[None] - om2 * t[:, None])
    assert np.max(np.abs(two - exact)) < 1e-13  # each mode with its own ω(k)


def test_real_field_V1_real_parts_and_product_averages():  # V1 N99, P178: ⟨Re(Ae^{iθ})Re(Be^{iθ})⟩ = ½Re(AB*)
    th = np.linspace(0, TWO_PI, 64, endpoint=False)
    assert np.max(np.abs(W.real_field(0.7, th) - 0.7 * np.cos(th))) < 1e-15
    A, B = 0.3 - 1.2j, -0.8 + 0.5j
    assert np.max(np.abs(W.real_field(A, th) - abs(A) * np.cos(th + np.angle(A)))) < 1e-14
    avg = np.mean(W.real_field(A, th) * W.real_field(B, th))
    assert avg == pytest.approx(0.5 * np.real(A * np.conj(B)), rel=1e-13)
    assert abs(avg - 0.5 * np.real(A * B)) > 0.1  # multiplying the complex fields first is wrong


# =====================================================================================================================
# C02 — the linearised free-surface problem (7.11)–(7.21): R01–R06, N12–N17, D02–D04
# =====================================================================================================================
def test_free_surface_V1_linear_conditions_hold_on_a_field():  # V1 (7.11), (7.12), (7.18), (7.21), (7.55)
    x = RNG.uniform(0, TWO_PI, 25)
    a, k = 0.01, 1.0
    for kH in (0.3, 1.0, 3.0, np.inf):
        H = kH / k
        om = float(W.omega_gravity(k, H, G))
        for t, d in ((0.0, 1), (0.37, 1), (1.1, -1)):
            r = ch07.free_surface_residuals(x, t, a, k, H, G, direction=d)
            assert np.max(np.abs(r["kinematic_linear"])) < 1e-14 * a * om  # (7.18)
            assert np.max(np.abs(r["dynamic_linear"])) < 1e-13 * a * G  # (7.21)
            assert np.max(np.abs(r["bottom"])) < 1e-15 * a * om  # (7.12) (deep: at z = −40/k, e^{−40})
            assert np.max(np.abs(r["laplace"])) < 1e-9 * a * om * k  # (7.11) 4th-order stencil, h = 0.01/k
    sig = 0.0727
    r = ch07.free_surface_residuals(x, 0.2, 1e-4, 300.0, 0.02, G, 1000.0, sig)  # a ripple with tension (7.55)
    geff = G + sig * 300.0 ** 2 / 1000.0
    assert np.max(np.abs(r["dynamic_linear"])) < 1e-13 * 1e-4 * geff


def test_free_surface_V1_exact_residuals_equal_independent_closed_forms():  # V1: the exact (7.16) and (7.19) residuals
    a, k, H = 0.03, 1.0, 1.2
    om = float(W.omega_gravity(k, H, G))
    x = np.linspace(0, TWO_PI, 37)
    t = 0.45
    th = k * x - om * t
    eta = a * np.cos(th)
    sh = np.sinh(k * H)
    u = a * om * np.cosh(k * (eta + H)) / sh * np.cos(th)  # (7.27) at z = η, written out here
    w = a * om * np.sinh(k * (eta + H)) / sh * np.sin(th)
    phi_t = -a * om ** 2 / k * np.cosh(k * (eta + H)) / sh * np.cos(th)
    kin = w - a * om * np.sin(th) - (-a * k * np.sin(th)) * u  # (7.16): φ_z − η_t − η_x φ_x
    dyn = phi_t + 0.5 * (u ** 2 + w ** 2) + G * eta  # (4.83) with p = 0 (7.19)
    r = ch07.free_surface_residuals(x, t, a, k, H, G)
    assert np.max(np.abs(r["kinematic_exact"] - kin)) < 1e-15 + 1e-12 * np.max(np.abs(kin))
    assert np.max(np.abs(r["dynamic_exact"] - dyn)) < 1e-15 + 1e-12 * np.max(np.abs(dyn))
    assert np.max(np.abs(r["kinematic_17"] - (w - a * om * np.sin(th)))) < 1e-15


def test_free_surface_V3_neglected_terms_scale_as_ka():  # V3 D03/D04: |exact residual| ∝ (ka)², relative ∝ ka
    ka = np.array([0.005, 0.01, 0.02, 0.04, 0.08])
    for kH in (1.0, np.inf):
        s = ch07.free_surface_residual_scan(ka, kH=kH, g=G)
        assert abs(observed_order(ka, s["kinematic_abs"]) - 2.0) < ORDER_TOL  # absolute, at fixed k
        assert abs(observed_order(ka, s["dynamic_abs"]) - 2.0) < ORDER_TOL
        assert abs(observed_order(ka, s["kinematic_rel"]) - 1.0) < ORDER_TOL  # ÷ aω: the dropped terms are O(ka)
        assert abs(observed_order(ka, s["dynamic_rel"]) - 1.0) < ORDER_TOL  # ÷ ag
        assert abs(s["slope_abs"] - observed_order(ka, s["kinematic_abs"])) < 1e-12  # the function's own slopes
        assert abs(s["slope_rel_dynamic"] - observed_order(ka, s["dynamic_rel"])) < 1e-12
        assert abs(observed_order(ka, s["kinematic_rel"]) - 2.0) > 0.5  # analysis V25's "(ka)² × aω" is not the scaling
        assert np.max(s["kinematic_linear"]) < 1e-15 and np.max(s["dynamic_linear"]) < 1e-15


def test_surface_normal_V1_unit_normal_and_kinematic_parity():  # V1 (7.13)–(7.16), N13, N14, R04 (ch04 (4.91))
    s = RNG.uniform(-2, 2, 40)
    n = ch07.surface_normal(s)
    assert np.max(np.abs(np.hypot(n[0], n[1]) - 1.0)) < 1e-15 and np.all(n[1] > 0)
    assert np.max(np.abs(n[0] * 1.0 - n[1] * (-s))) < 1e-15  # ∥ ∇f = (−η_x, 1)
    eta_t = RNG.uniform(-1, 1, 40)
    Us = ch07.surface_velocity(eta_t)
    assert np.max(np.abs((n[0] * Us[0] + n[1] * Us[1]) * np.sqrt(s ** 2 + 1) - eta_t)) < 1e-14  # (n·U_s)|∇f| = η_t
    # ch04's material-surface residual Dη_imp/Dt with the level-set function η_imp = z − η (R1) equals (7.16)
    a, k, H, t = 0.04, 1.0, 1.5, 0.3
    om = float(W.omega_gravity(k, H, G))
    xs = np.linspace(0.1, 6.0, 9)
    pts = np.stack([xs, a * np.cos(k * xs - om * t)])
    f_imp = lambda X, tt: X[1] - a * np.cos(k * X[0] - om * tt)  # noqa: E731
    uvec = lambda X, tt: np.stack([np.asarray(ch07.wave_fields(X[0], X[1], tt, a, k, H, G)[c]) for c in ("u", "w")])  # noqa: E731
    res = IF.kinematic_bc_residual(f_imp, uvec, pts, t, h=1e-5, ht=1e-5)
    ref = ch07.free_surface_residuals(xs, t, a, k, H, G)["kinematic_exact"]
    assert np.max(np.abs(res - ref)) < 1e-8 * a * om
    wrong = IF.kinematic_bc_residual(lambda X, tt: a * np.cos(k * X[0] - om * tt) + 0 * X[1], uvec, pts, t, h=1e-5,
                                     ht=1e-5)
    assert np.max(np.abs(wrong - ref)) > 0.1 * a * om  # passing the elevation as the level set (ch04's η) is wrong


def test_linear_bernoulli_V1_parity_with_ch04_unsteady_bernoulli():  # V1 (7.20), (7.30), R06
    phit = RNG.uniform(-3, 3, 20)
    z = RNG.uniform(-5, 0.5, 20)
    lb = ch07.linear_bernoulli_pressure(phit, z, 1025.0, G)
    assert maxrel(lb["p"], BE.unsteady_bernoulli_pressure(phit, 0.0, z, 1025.0, g=G)) < 1e-14
    assert maxrel(lb["p_prime"], -1025.0 * phit) < 1e-12
    assert maxrel(ch07.linear_bernoulli_pressure(0.0, z, 1025.0, G)["p"], -1025.0 * G * z) < 1e-15  # at rest
    a, k, H = 0.2, 0.5, 4.0  # with the wave: −ρφ_t of (7.26) is the p′ of (7.31)
    om = float(W.omega_gravity(k, H, G))
    x, zz, t = RNG.uniform(0, 12, 30), RNG.uniform(-H, 0, 30), 0.8
    phi_t = -a * om ** 2 / k * np.cosh(k * (zz + H)) / np.sinh(k * H) * np.cos(k * x - om * t)
    assert maxrel(ch07.linear_bernoulli_pressure(phi_t, zz, 1000.0, G)["p_prime"],
                  ch07.wave_fields(x, zz, t, a, k, H, G, 1000.0)["p_prime"]) < 1e-12


def test_kinematic_condition_V2_derivation():  # V2 D02 (★★): (7.13) with f = z − η → (7.14), (7.15) → (7.16)
    x, z, t = sp.symbols("x z t", real=True)
    eta = sp.Function("eta")(x, t)
    phi = sp.Function("phi")(x, z, t)
    f = z - eta  # step 1
    gf = sp.Matrix([sp.diff(f, x), sp.diff(f, z)])  # step 2
    assert list(gf) == [-sp.diff(eta, x), 1]
    nrm = sp.sqrt(sp.diff(eta, x) ** 2 + 1)
    n = gf / nrm  # step 3, (7.14)
    Us = sp.Matrix([0, sp.diff(eta, t)])  # step 4, (7.15)
    u = sp.Matrix([sp.diff(phi, x), sp.diff(phi, z)])
    lhs = sp.simplify((n.dot(u) - n.dot(Us)) * nrm)  # step 5
    step6 = -sp.diff(phi, x) * sp.diff(eta, x) + sp.diff(phi, z) - sp.diff(eta, t)
    assert sp.simplify(lhs - step6) == 0
    eq716 = sp.diff(phi, z) - (sp.diff(eta, t) + sp.diff(eta, x) * sp.diff(phi, x))  # step 7
    assert sp.simplify(step6 - eq716) == 0
    material = sp.diff(f, t) + sp.diff(phi, x) * sp.diff(f, x) + sp.diff(phi, z) * sp.diff(f, z)  # step 8, D(z − η)/Dt
    assert sp.simplify(material - eq716) == 0
    sl = 0.37  # the coded (7.14) is this n
    nn = ch07.surface_normal(sl)
    assert np.allclose(nn.ravel(), [-sl / np.sqrt(sl ** 2 + 1), 1 / np.sqrt(sl ** 2 + 1)], atol=1e-15)


def test_linearised_kinematic_V2_derivation():  # V2 D03 (★★): (7.16) → (7.17) → (7.18); both dropped terms O(ka·aω)
    x, z, t, TH = sp.symbols("x z t Theta", real=True)
    a, k, H, om, h = sp.symbols("a k H omega h", positive=True)
    th = k * x - om * t
    eta = a * sp.cos(th)
    phi = a * om / k * sp.cosh(k * (z + H)) / sp.sinh(k * H) * sp.sin(th)  # (7.26)
    F = sp.diff(phi, z) - sp.diff(eta, t) - sp.diff(eta, x) * sp.diff(phi, x)  # (7.16) as a field
    ser = sp.expand(sp.series(F.subs(z, eta), a, 0, 3).removeO())
    assert sp.simplify(ser.coeff(a, 1)) == 0  # (7.18): the O(a) part is φ_z(0) − η_t = 0
    kept = (eta * sp.diff(phi, z, 2) - sp.diff(eta, x) * sp.diff(phi, x)).subs(z, 0)  # steps 4–6: ηφ_zz and η_xφ_x
    assert sym_zero(ser.coeff(a, 2) * a ** 2 - kept)
    size = sp.simplify((kept / (a * om * k * a)).subs(x, (TH + om * t) / k).subs(H, h / k))  # ÷ (ka)(aω)
    assert size.free_symbols <= {TH, h}  # the O(a²) residual is ka·aω times a function of phase and kH only (step 6)


def test_linearised_dynamic_V2_derivation():  # V2 D04 (★★): (7.19) + (4.83) → (7.20) → (7.21)
    x, z, t, TH = sp.symbols("x z t Theta", real=True)
    a, k, H, om, g, h = sp.symbols("a k H omega g h", positive=True)
    th = k * x - om * t
    eta = a * sp.cos(th)
    phi = a * om / k * sp.cosh(k * (z + H)) / sp.sinh(k * H) * sp.sin(th)
    bern = sp.diff(phi, t) + (sp.diff(phi, x) ** 2 + sp.diff(phi, z) ** 2) / 2 + g * z  # (4.83) with p = 0, C = 0
    ser = sp.expand(sp.series(bern.subs(z, eta), a, 0, 3).removeO())
    lin = (sp.diff(phi, t) + g * eta).subs(z, 0)  # steps 5–8: (7.21)
    assert sym_zero(ser.coeff(a, 1) * a - lin)
    assert sym_zero(lin.subs(om, sp.sqrt(g * k * sp.tanh(k * H))))  # it holds exactly with (7.28)
    dropped = (eta * sp.diff(phi, t, z) + (sp.diff(phi, x) ** 2 + sp.diff(phi, z) ** 2) / 2).subs(z, 0)  # steps 3, 7
    assert sym_zero(ser.coeff(a, 2) * a ** 2 - dropped)
    size = sp.simplify((dropped / (k * a * a * om ** 2 / k)).subs(x, (TH + om * t) / k).subs(H, h / k))
    assert size.free_symbols <= {TH, h}  # ka × φ_t: the square term and the transfer are both O(ka)


# =====================================================================================================================
# C03 — the dispersion relation (7.28): N18–N25, D05, D06; V5 Fenton–McKee (1990), Guo (2002)
# =====================================================================================================================
def test_potential_V2_derivation():  # V2 D05 (★★): (7.22) → (7.23) → (7.24) → (7.25) → (7.26), step by step
    x, z, t = sp.symbols("x z t", real=True)
    a, k, H, om = sp.symbols("a k H omega", positive=True)
    A, B = sp.symbols("A B")
    f = sp.Function("f")
    th = k * x - om * t
    trial = f(z) * sp.sin(th)  # step 1 (7.22)
    lap = sp.diff(trial, x, 2) + sp.diff(trial, z, 2)
    assert sp.simplify(lap - (sp.diff(f(z), z, 2) - k ** 2 * f(z)) * sp.sin(th)) == 0  # step 2
    sol = sp.dsolve(sp.diff(f(z), z, 2) - k ** 2 * f(z), f(z)).rhs  # steps 3–4 (7.23)
    assert {sp.exp(k * z), sp.exp(-k * z)} <= sol.atoms(sp.exp)
    fz = A * sp.exp(k * z) + B * sp.exp(-k * z)
    Bsol = sp.solve(sp.diff(fz, z).subs(z, -H), B)[0]  # step 5 (7.24)
    assert sp.simplify(Bsol - A * sp.exp(-2 * k * H)) == 0
    Asol = sp.solve(sp.Eq(k * (A - Bsol), om * a), A)[0]  # steps 6–8 (7.25)
    assert sp.simplify(Asol - a * om / (k * (1 - sp.exp(-2 * k * H)))) == 0
    fz2 = fz.subs(B, Bsol)
    assert sym_zero(fz2 - 2 * A * sp.exp(-k * H) * sp.cosh(k * (z + H)))  # step 9: the cosh regrouping
    assert sym_zero(2 * Asol * sp.exp(-k * H) - a * om / (k * sp.sinh(k * H)))  # step 10
    phi = fz2.subs(A, Asol) * sp.sin(th)
    assert sym_zero(phi - a * om / k * sp.cosh(k * (z + H)) / sp.sinh(k * H) * sp.sin(th))  # step 11 (7.26)
    wrong = sp.diff((A * sp.exp(k * z) + A * sp.exp(2 * k * H) * sp.exp(-k * z)), z).subs(z, -H)
    assert not sym_zero(wrong)  # the wrong-sign bottom constant fails (7.12)
    s = ch07.surface_wave_sympy()  # the coded construction: every residual 0, the wrong one not
    assert all(sp.simplify(v) == 0 for v in s["residuals"].values())
    assert sp.simplify(s["wrong_bottom_residual"]) != 0
    assert sp.simplify(s["B_over_A"] - sp.exp(-2 * sp.Symbol("k", positive=True) * sp.Symbol("H", positive=True))) == 0


def test_dispersion_V2_derivation():  # V2 D06 (★★): (7.21) with (7.26) → ω² = gk tanh kH (7.28) and the T–λ form
    x, z, t = sp.symbols("x z t", real=True)
    a, k, H, om, g, lam, T = sp.symbols("a k H omega g lambda T", positive=True)
    th = k * x - om * t
    phi = a * om / k * sp.cosh(k * (z + H)) / sp.sinh(k * H) * sp.sin(th)
    phit = sp.diff(phi, t)  # step 1
    assert sp.simplify(phit + a * om ** 2 / k * sp.cosh(k * (z + H)) / sp.sinh(k * H) * sp.cos(th)) == 0
    line = sp.simplify(phit.subs(z, 0) / sp.cos(th))  # steps 2–3: −(aω²/k) coth kH
    assert sym_zero(line + a * om ** 2 / k * sp.coth(k * H))
    om2 = sp.solve(sp.Eq(om ** 2 / k * sp.coth(k * H), g), om ** 2)  # steps 4–5
    assert len(om2) == 1 and sym_zero(om2[0] - g * k * sp.tanh(k * H))
    assert sym_zero(sp.sqrt(om2[0]) - sp.sqrt(g * k * sp.tanh(k * H)))  # step 6, the positive root
    Tform = sp.sqrt(2 * sp.pi * lam / g * sp.coth(2 * sp.pi * H / lam))  # step 7
    assert sym_zero((2 * sp.pi / Tform) ** 2 - g * (2 * sp.pi / lam) * sp.tanh(2 * sp.pi * H / lam))
    lam_v, H_v = 37.0, 6.0
    assert float(W.period_from_wavelength(lam_v, H_v, G)) == pytest.approx(
        float(Tform.subs({lam: lam_v, H: H_v, g: G})), rel=1e-14)


def test_dispersion_V1_parity_with_ch04_and_identity():  # V1 (7.28): ω² = gk tanh kH; ch04's test field is this wave
    k = np.logspace(-3, 2, 60)
    for H in (0.5, 7.0, 300.0, np.inf):
        om = W.omega_gravity(k, H, G)
        th = np.ones_like(k) if np.isinf(H) else np.tanh(k * H)
        assert maxrel(om ** 2, G * k * th) < 1e-14
        for kk in (0.05, 0.7, 3.0):
            f7 = ch07.wave_fields(1.3, -0.2, 2.1, 0.05, kk, H, G)
            assert np.isfinite(float(f7["u"])) and np.isfinite(float(f7["w"]))
            if kk * H > 300:  # ch04's plain cosh/sinh overflow there (kH = 900); ch07's forms do not
                continue
            f4 = ch04.linear_wave_surface(1.3, -0.2, 2.1, 0.05, kk, G, H)
            assert f4["omega"] == pytest.approx(float(W.omega_gravity(kk, H, G)), rel=1e-14)
            f7 = ch07.wave_fields(1.3, -0.2, 2.1, 0.05, kk, H, G)
            assert float(f7["u"]) == pytest.approx(float(f4["u"]), rel=1e-12, abs=1e-15)
            assert float(f7["w"]) == pytest.approx(float(f4["w"]), rel=1e-12, abs=1e-15)
    assert float(W.omega_gravity(-0.4, 3.0, G)) == float(W.omega_gravity(0.4, 3.0, G))  # |k|
    assert ch07.G == 9.81 and ch07.G_BOOK == 9.81


def test_dispersion_V7_deep_and_shallow_limits():  # V7 (7.45), (7.49): ω → √(gk), → k√(gH)(1 − (kH)²/6)
    k = 0.8
    assert float(W.omega_gravity(k, 40.0 / k, G)) == pytest.approx(np.sqrt(G * k), rel=1e-12)
    assert float(W.omega_gravity(k, np.inf, G)) == np.sqrt(G * k)
    for kH in (1e-3, 1e-2):
        H = kH / k
        assert float(W.omega_gravity(k, H, G)) == pytest.approx(k * np.sqrt(G * H) * (1 - kH ** 2 / 6), rel=kH ** 4)
    lam = np.logspace(-1, 4, 30)
    assert maxrel(W.phase_speed(TWO_PI / lam, np.inf, G), np.sqrt(G * lam / TWO_PI)) < 1e-14  # (7.45) λ form


def test_wavenumber_from_omega_V1_inverse_is_identity():  # V1 (7.28) inverted by brentq over 7 decades of ω, 7 of H
    om = np.logspace(-4, 3, 36)
    for H in (1e-3, 0.1, 1.0, 10.0, 1e3, 1e4, np.inf):
        k = W.wavenumber_from_omega(om, H, G)
        assert np.max(np.abs(W.omega_gravity(k, H, G) / om - 1)) < 1e-12  # per point, over 7 decades
        if np.isinf(H):
            assert maxrel(k, om ** 2 / G) < 1e-15  # the deep closed form
    k = W.wavenumber_from_omega(om[5:], 0.3, G, 0.0727, 1000.0)  # with surface tension (7.56)
    assert np.max(np.abs(W.omega_capillary_gravity(k, 0.3, 0.0727, 1000.0, G) / om[5:] - 1)) < 1e-12
    assert float(W.wavelength_from_period(12.0, np.inf, G)) == pytest.approx(G * 144.0 / TWO_PI, rel=1e-13)
    lam = W.wavelength_from_period(np.array([4.0, 12.0, 60.0]), 25.0, G)
    assert maxrel(W.period_from_wavelength(lam, 25.0, G), [4.0, 12.0, 60.0]) < 1e-12
    assert float(W.wavenumber_from_omega(0.0, 5.0, G)) == 0.0
    with pytest.raises(ValueError):
        W.wavenumber_from_omega(-1.0, 5.0, G)


def _kd_exact(x_nd):
    """Exact kd for the non-dimensional frequency x = ω√(d/g) (our inverse dispersion, d = 1 m)."""
    om = x_nd * np.sqrt(G / 1.0)
    return W.wavenumber_from_omega(om, 1.0, G) * 1.0


@needs_ref
def test_dispersion_V5_fenton_mckee_1990():  # V5 Fenton & McKee (1990) Eq. (21): max λ error 1.7 %, optimal ν = 1.49
    ref = ref_json()["fenton_mckee_1990"]
    xnd = np.logspace(-2, 2, 4001)
    kd = _kd_exact(xnd)
    om = xnd * np.sqrt(G)
    fm = np.asarray(W.fenton_mckee_kh(om, 1.0, G))
    lam_err = np.abs(kd / fm - 1.0)  # error of the approximate wavelength L = 2π/k
    assert np.max(lam_err) <= ref["max_error_wavelength_percent"] / 100.0  # "always better than 1.7 %"
    assert round(100 * np.max(lam_err), 1) == ref["max_error_wavelength_percent"]  # and it is 1.7 % to the printed digit
    assert lam_err[0] < 1e-4 and lam_err[-1] < 1e-12  # exact in both limits

    def maxerr(nu):
        a = xnd ** 2 * (1.0 / np.tanh(xnd ** nu)) ** (1.0 / nu)  # the paper's one-parameter family
        return np.max(np.abs(kd / a - 1.0))

    nu_opt = minimize_scalar(maxerr, bounds=(1.2, 1.8), method="bounded", options={"xatol": 1e-6}).x
    assert abs(nu_opt / ref["optimal_nu"] - 1.0) < 0.01  # published benchmark tolerance (1 %)
    assert abs(maxerr(1.5) - np.max(lam_err)) < 1e-12  # ν = 3/2 is what fenton_mckee_kh codes


@needs_ref
def test_dispersion_V5_guo_2002():  # V5 Guo (2002) β = 2.4908: 0.75 %; Fenton (2006): 5/2 variant ≈ half the FM error
    ref = ref_json()
    g02 = ref["guo_2002"]
    xnd = np.logspace(-2, 2, 4001)
    kd = _kd_exact(xnd)
    beta = g02["beta"]
    guo_orig = xnd ** 2 * (1.0 - np.exp(-xnd ** beta)) ** (-1.0 / beta)
    e_orig = np.max(np.abs(guo_orig / kd - 1.0))
    assert abs(100 * e_orig - g02["max_rel_error_percent"]) < 0.005  # reproduced to the printed two decimals
    om = xnd * np.sqrt(G)
    e_guo = np.abs(np.asarray(W.guo_kh(om, 1.0, G)) / kd - 1.0)
    e_fm = np.abs(np.asarray(W.fenton_mckee_kh(om, 1.0, G)) / kd - 1.0)
    assert e_guo[0] < 1e-4 and e_guo[-1] < 1e-12  # exact in both limits (Fenton 2006)
    assert 0.4 < np.max(e_guo) / np.max(e_fm) < 0.6  # "about half the error of equation (2)"
    assert 0.007 < np.max(e_guo) < 0.008  # "about 0.7 %" (measured 0.79 %, see the report)


# =====================================================================================================================
# C04 — phase speed (7.29) and the deep/shallow limits (7.45)–(7.52): R07, N26, N37–N46, D07–D09
# =====================================================================================================================
def test_phase_speed_V2_derivation_longer_is_faster():  # V2 D07 (★) and D09 step 3: dc/dλ > 0; c/√(gH) ≈ 1 − (kH)²/6
    lam, H, g, xx, kH = sp.symbols("lambda H g x kH", positive=True)
    c2 = g * lam / (2 * sp.pi) * sp.tanh(2 * sp.pi * H / lam)  # (7.29) λ form, squared
    d = sp.diff(lam * sp.tanh(2 * sp.pi * H / lam), lam)
    xs = 2 * sp.pi * H / lam
    assert sym_zero(d - (sp.sinh(2 * xs) / 2 - xs) / sp.cosh(xs) ** 2)  # step 4
    assert sp.series(sp.sinh(2 * xx) / 2 - xx, xx, 0, 6).removeO().equals(2 * xx ** 3 / 3 + 2 * xx ** 5 / 15)  # > 0
    assert sym_zero(sp.diff(c2, lam) - g / (2 * sp.pi) * d)
    ratio = sp.series(sp.sqrt(sp.tanh(kH) / kH), kH, 0, 4).removeO()  # D09 step 3
    assert sp.simplify(ratio - (1 - kH ** 2 / 6)) == 0


def test_phase_speed_V7_limits_monotonicity_and_bound():  # V7 (7.29), (7.45), (7.49)
    k = np.logspace(-4, 2, 300)
    for H in (0.3, 5.0, 200.0):
        c = W.phase_speed(k, H, G)
        assert np.all(np.diff(c) < 0)  # longer waves (smaller k) are faster: dispersive
        assert np.all(c < np.sqrt(G * H))  # never faster than √(gH)
        assert c[0] == pytest.approx(np.sqrt(G * H), rel=(k[0] * H) ** 2 / 6 * 1.01)
    assert float(W.phase_speed(2.0, 20.0, G)) == pytest.approx(np.sqrt(G / 2.0), rel=1e-15)  # kH = 40
    assert maxrel(W.phase_speed(k, 5.0, G) * k, W.omega_gravity(k, 5.0, G)) < 1e-14  # c = ω/k


def test_depth_regime_V1_book_thresholds_and_errors():  # V1 N38, N43 (R11): errors, not only labels
    d = W.depth_regime(2.0, 1.0)
    assert d["H_over_lambda"] == pytest.approx(1 / np.pi, rel=1e-15)  # kH = 2 ⇔ H/λ = 0.318
    assert d["deep_error"] == pytest.approx(1 - np.sqrt(np.tanh(2.0)), rel=1e-14)
    assert round(100 * d["deep_error"], 2) == 1.82
    s = W.depth_regime(TWO_PI * 0.07, 1.0)
    assert s["shallow_error"] == pytest.approx(1 - np.sqrt(np.tanh(TWO_PI * 0.07) / (TWO_PI * 0.07)), rel=1e-14)
    assert round(100 * s["shallow_error"], 2) == 3.04
    assert W.depth_regime(2.01, 1.0)["regime"] == "deep" and W.depth_regime(1.99, 1.0)["regime"] == "intermediate"
    assert W.depth_regime(TWO_PI * 0.069, 1.0)["regime"] == "shallow"
    assert [ch07.wave_regime_label(v) for v in (3.0, 1.0, 0.2)] == ["deep", "intermediate", "shallow"]
    arr = W.depth_regime(np.array([0.1, 1.0, 5.0]), 1.0)
    assert list(arr["regime"]) == ["shallow", "intermediate", "deep"]


def test_pressure_response_V1_surface_bottom_and_limits():  # V1 (7.31), (7.48), (7.52): N26, N42, N46
    k, H = 0.4, 6.0
    assert float(ch07.pressure_response(k, 0.0, H)) == pytest.approx(1.0, rel=1e-15)
    assert float(ch07.pressure_response(k, -H, H)) == pytest.approx(1 / np.cosh(k * H), rel=1e-14)
    z = np.linspace(-H, 0, 11)
    assert maxrel(ch07.pressure_response(k, z, H), np.cosh(k * (z + H)) / np.cosh(k * H)) < 1e-14
    lam = 30.0
    assert float(ch07.pressure_response(TWO_PI / lam, -lam / 2, np.inf)) == pytest.approx(np.exp(-np.pi), rel=1e-14)
    assert float(ch07.pressure_response(1e-4, -3.0, 5.0)) == pytest.approx(1.0, abs=(1e-4 * 5.0) ** 2)  # hydrostatic,
    # to O((kH)²/2): cosh(2e-4)/cosh(5e-4) = 1 − 1.05e-7
    big = ch07.pressure_response(1.0, np.linspace(-500, 0, 6), 500.0)  # kH = 500: no overflow
    assert np.all(np.isfinite(big)) and maxrel(big, np.exp(np.linspace(-500, 0, 6))) < 1e-12
    for fn, ref in ((W.cosh_over_sinh, lambda k_, z_, H_: np.cosh(k_ * (z_ + H_)) / np.sinh(k_ * H_)),
                    (W.sinh_over_sinh, lambda k_, z_, H_: np.sinh(k_ * (z_ + H_)) / np.sinh(k_ * H_)),
                    (W.cosh_over_cosh, lambda k_, z_, H_: np.cosh(k_ * (z_ + H_)) / np.cosh(k_ * H_))):
        assert maxrel(fn(k, z, H), ref(k, z, H)) < 1e-13
    dp = W.depth_profiles(k, z, H)
    assert set(dp) == {"cosh_sinh", "sinh_sinh", "cosh_cosh"}


def test_wave_fields_V1_closed_forms_and_streamfunction():  # V1 (7.26), (7.27), (7.31), (7.37): N23, N24, N26, N30
    a, k, H, rho = 0.3, 0.45, 5.0, 1025.0
    om = float(W.omega_gravity(k, H, G))
    x, z, t = RNG.uniform(0, 20, 40), RNG.uniform(-H, 0, 40), 1.7
    f = ch07.wave_fields(x, z, t, a, k, H, G, rho)
    th = k * x - om * t
    C, S = np.cosh(k * (z + H)) / np.sinh(k * H), np.sinh(k * (z + H)) / np.sinh(k * H)
    assert maxrel(f["phi"], a * om / k * C * np.sin(th)) < 1e-13
    assert maxrel(f["u"], a * om * C * np.cos(th)) < 1e-13 and maxrel(f["w"], a * om * S * np.sin(th)) < 1e-13
    assert maxrel(f["psi"], a * om / k * S * np.cos(th)) < 1e-13
    assert maxrel(f["p_prime"], rho * G * a * np.cosh(k * (z + H)) / np.cosh(k * H) * np.cos(th)) < 1e-13
    assert maxrel(f["eta"], a * np.cos(th)) < 1e-15 and f["c"] == pytest.approx(om / k)
    psi = lambda X, Z: ch07.wave_fields(X, Z, t, a, k, H, G)["psi"]  # noqa: E731
    u, w = SF.velocity_from_streamfunction_2d(psi, x, z, h=1e-5)  # ch04 convention with y ↦ z: u = ψ_z, w = −ψ_x
    assert maxrel(u, f["u"]) < 1e-8 and maxrel(w, f["w"]) < 1e-8


def test_wave_fields_V2_potential_relations_symbolic():  # V2: ∇²φ = 0, u = φ_x = ψ_z, w = φ_z = −ψ_x, p′ = −ρφ_t
    s = ch07.surface_wave_sympy()
    x, z, t = sp.symbols("x z t", real=True)
    a, k, H, om, g, rho = sp.symbols("a k H omega g rho", positive=True)
    phi, psi, u, w = s["phi"], s["psi"], s["u"], s["w"]
    assert sym_zero(sp.diff(phi, x, 2) + sp.diff(phi, z, 2))
    assert sym_zero(sp.diff(phi, x) - sp.diff(psi, z)) and sym_zero(sp.diff(phi, z) + sp.diff(psi, x))
    assert sym_zero(sp.diff(psi, z) - u) and sym_zero(-sp.diff(psi, x) - w)
    pp = rho * g * a * sp.cosh(k * (z + H)) / sp.cosh(k * H) * sp.cos(k * x - om * t)  # (7.31) last form
    assert sym_zero((-rho * sp.diff(phi, t) - pp).subs(om, sp.sqrt(g * k * sp.tanh(k * H))))


def test_wave_fields_V3_laplacian_residual_is_second_order_truncation():  # V3: 5-point ∇²φ → 0 at O(h²)
    a, k, H, t = 0.2, 0.9, 3.0, 0.4
    X, Z = np.array([0.3, 1.9, 4.4]), np.array([-0.4, -1.5, -2.6])
    phi = lambda x_, z_: ch07.wave_fields(x_, z_, t, a, k, H, G)["phi"]  # noqa: E731
    hs, errs = [0.2, 0.1, 0.05, 0.025], []
    for h in hs:
        lap = (phi(X + h, Z) + phi(X - h, Z) + phi(X, Z + h) + phi(X, Z - h) - 4 * phi(X, Z)) / h ** 2
        errs.append(np.max(np.abs(lap)))
    assert abs(observed_order(hs, errs) - 2.0) < ORDER_TOL  # harmonic: only the stencil's truncation remains


def test_wave_fields_V7_deep_shallow_mirror_and_overflow():  # V7 (7.46)–(7.48), (7.51)–(7.52), (7.61)
    a, k = 0.1, 1.0
    om = np.sqrt(G * k)
    x, z, t = RNG.uniform(0, 6, 20), RNG.uniform(-3, 0, 20), 0.3
    fd = ch07.wave_fields(x, z, t, a, k, np.inf, G)
    assert maxrel(fd["u"], a * om * np.exp(k * z) * np.cos(k * x - om * t)) < 1e-14  # (7.47)
    assert maxrel(fd["p_prime"], 1000 * G * a * np.exp(k * z) * np.cos(k * x - om * t)) < 1e-14  # (7.48)
    kh, H = 1e-3, 1e-3 / k  # shallow (7.51)–(7.52): relative corrections O((kH)²)
    zs = np.linspace(-H, 0, 7)
    fs = ch07.wave_fields(0.2, zs, 0.0, a, k, H, G)
    oms = float(W.omega_gravity(k, H, G))
    assert maxrel(fs["u"], a * oms / kh * np.cos(0.2) * np.ones_like(zs)) < kh ** 2
    assert maxrel(fs["w"], a * oms * (1 + zs / H) * np.sin(0.2)) < kh ** 2
    assert maxrel(fs["p_prime"], 1000 * G * a * np.cos(0.2) * np.ones_like(zs)) < kh ** 2
    fr = ch07.wave_fields(x, z, t, a, k, 2.0, G, direction=+1)  # mirror: (x, u) → (−x, −u) for the left-going wave
    fl = ch07.wave_fields(-x, z, t, a, k, 2.0, G, direction=-1)
    assert maxrel(fl["u"], -fr["u"]) < 1e-14 and maxrel(fl["w"], fr["w"]) < 1e-14 and maxrel(fl["eta"], fr["eta"]) < 1e-15
    fb = ch07.wave_fields(0.0, np.array([-500.0, -250.0, 0.0]), 0.0, a, k, 500.0, G)  # kH = 500
    assert all(np.all(np.isfinite(fb[c])) for c in ("u", "w", "phi", "psi", "p_prime"))


def test_dispersion_state_V1_explainer_numbers():  # V1 Part C 2.8 (E1 state)
    s = ch07.dispersion_state(50.0, 10.0, G)
    k = TWO_PI / 50.0
    assert s["k"] == pytest.approx(k) and s["kH"] == pytest.approx(10 * k)
    assert s["c"] == pytest.approx(float(W.phase_speed(k, 10.0, G)), rel=1e-14)
    assert s["c"] == pytest.approx(8.15, abs=0.005) and s["cg"] / s["c"] == pytest.approx(0.705, abs=5e-4)  # D21 check
    assert s["T"] == pytest.approx(TWO_PI / s["omega"]) and s["c_shallow"] == pytest.approx(np.sqrt(G * 10.0))
    assert s["p_bottom_fraction"] == pytest.approx(1 / np.cosh(10 * k), rel=1e-14)
    assert s["t_cross_1000km_cg_h"] == pytest.approx(1e6 / s["cg"] / 3600)
    assert s["regime"] == "intermediate" and ch07.dispersion_state(50.0)["p_bottom_fraction"] == 0.0


# =====================================================================================================================
# C05 — particle orbits (7.32)–(7.37): R08, N27–N30, N40, N44, D10, D11
# =====================================================================================================================
def test_orbits_V2_derivation():  # V2 D10 (★★): (7.34a, b) integrated → (7.35a, b), zero-mean constants
    x0, z0, t = sp.symbols("x0 z0 t", real=True)
    a, k, H, om = sp.symbols("a k H omega", positive=True)
    C = sp.cosh(k * (z0 + H)) / sp.sinh(k * H)
    S = sp.sinh(k * (z0 + H)) / sp.sinh(k * H)
    ph = k * x0 - om * t
    rhs_x, rhs_z = a * om * C * sp.cos(ph), a * om * S * sp.sin(ph)  # (7.34a, b), steps 4–5
    xi = sp.integrate(rhs_x, t)  # step 6
    ze = sp.integrate(rhs_z, t)  # step 7
    assert sp.simplify(xi - (-a * C * sp.sin(ph))) == 0 and sp.simplify(ze - a * S * sp.cos(ph)) == 0
    T = 2 * sp.pi / om
    assert sp.simplify(sp.integrate(xi, (t, 0, T))) == 0 and sp.simplify(sp.integrate(ze, (t, 0, T))) == 0  # step 8
    assert sym_zero((xi / (a * C)) ** 2 + (ze / (a * S)) ** 2 - 1)  # D11 steps 2–3 (7.36)
    assert sym_zero(sp.sqrt(C ** 2 - S ** 2) * a - a / sp.sinh(k * H))  # D11 step 5: focal half-distance
    assert sym_zero(S / C - sp.tanh(k * (z0 + H)))  # D11 step 4


def test_orbit_V1_ellipses_clockwise_and_constant_foci():  # V1 (7.35)–(7.36), D11, Figs. 7.3–7.4 (N169, N170)
    a, k, H = 0.2, 0.6, 4.0
    om = float(W.omega_gravity(k, H, G))
    t = np.linspace(0, TWO_PI / om, 201)
    foc = []
    for z0 in (-0.1, -1.0, -2.5, -3.9):
        xi, ze = ch07.orbit_linear(1.2, z0, t, a, k, H, G)
        ax = ch07.orbit_semi_axes(z0, a, k, H)
        assert np.max(np.abs(xi ** 2 / ax["A"] ** 2 + ze ** 2 / ax["B"] ** 2 - 1)) < 1e-12  # (7.36)
        area = 0.5 * np.sum(xi[:-1] * ze[1:] - xi[1:] * ze[:-1])
        assert area < 0 and ax["sense"] == "cw"  # clockwise in the x–z plane for a right-going wave
        foc.append(np.sqrt(ax["A"] ** 2 - ax["B"] ** 2))
    assert np.ptp(foc) < 1e-13 and foc[0] == pytest.approx(a / np.sinh(k * H), rel=1e-12)  # same foci at every depth
    assert ch07.orbit_semi_axes(-4.0, a, k, 4.0)["B"] == pytest.approx(0.0, abs=1e-15)  # the bottom parcel slides
    s = ch07.orbit_semi_axes(0.0, 1.0, 3.0, 1.0)  # kH = 3 surface: A = 1.005a, B = a (D11 check)
    assert s["A"] == pytest.approx(1 / np.tanh(3.0), rel=1e-14) and round(float(s["A"]), 3) == 1.005
    assert float(ch07.orbit_semi_axes(0.0, 1.0, 0.3, 1.0)["A"]) == pytest.approx(3.43, abs=0.005)
    d = ch07.orbit_semi_axes(np.array([-0.5, -2.0]), a, k, np.inf)  # deep: circles of radius ae^{kz₀} (7.46)
    assert maxrel(d["A"], a * np.exp(k * np.array([-0.5, -2.0]))) < 1e-14 and maxrel(d["B"], d["A"]) < 1e-15
    kh = 1e-3  # shallow (7.50): A → a/kH, B → a(1 + z/H)
    sh = ch07.orbit_semi_axes(np.array([-0.2e-3, -0.7e-3]), a, 1.0, kh)
    assert maxrel(sh["A"], [a / kh] * 2) < kh ** 2 and maxrel(sh["B"], a * (1 + np.array([-0.2, -0.7]))) < kh ** 2


def test_particle_path_V1_linear_model_and_pathline_parity():  # V1 (7.32)–(7.34): R08, N27, N28
    a, k, H = 0.1, 0.8, 3.0
    om = float(W.omega_gravity(k, H, G))
    te = np.linspace(0, 3 * TWO_PI / om, 61)
    p = ch07.particle_path(0.5, -0.7, te, a, k, H, G, model="linear")
    xi, ze = ch07.orbit_linear(0.5, -0.7, te, a, k, H, G)
    assert np.max(np.abs(p["x"] - (0.5 + xi))) < 1e-9 * a and np.max(np.abs(p["z"] - (-0.7 + ze))) < 1e-9 * a
    ex = ch07.particle_path(0.5, -0.7, te, a, k, H, G, model="exact")
    uvec = lambda r, t: np.array([float(ch07.wave_fields(r[0], r[1], t, a, k, H, G)[c]) for c in ("u", "w")])  # noqa: E731
    pl = KN.pathline(uvec, [ex["x"][0], ex["z"][0]], 0.0, te[1:])
    assert np.max(np.abs(pl[0] - ex["x"][1:])) < 1e-8 * a and np.max(np.abs(pl[1] - ex["z"][1:])) < 1e-8 * a
    with pytest.raises(ValueError):
        ch07.particle_path(0.0, -1.0, te, a, k, H, G, model="nope")


def test_orbit_state_V1_explainer_numbers():  # V1 Part C 2.11 (E2 state)
    s = ch07.orbit_state(-0.5, 0.1, 1.0, 3.0, G)
    ax = ch07.orbit_semi_axes(-0.5, 0.1, 1.0, 3.0)
    assert s["A"] == pytest.approx(float(ax["A"])) and s["B_over_A"] == pytest.approx(np.tanh(2.5), rel=1e-13)
    assert s["drift_per_period"] == pytest.approx(s["drift_speed"] * s["T"]) and s["eulerian_mean"] == 0.0
    assert s["orbital_speed"] == pytest.approx(s["omega"] * s["A"]) and np.isnan(s["drift_numeric"])
    e = ch07.orbit_state(-0.5, 0.02, 1.0, 3.0, G, periods=10, model="exact")
    assert e["drift_numeric"] == pytest.approx(e["drift_per_period"], rel=0.02)  # O(ka) = 2 % at ka = 0.02


# =====================================================================================================================
# C06 — wave energy and energy flux (7.38)–(7.44): N31–N36, D12–D14
# =====================================================================================================================
def test_wave_energy_V1_quadrature_equals_closed_form():  # V1 (7.38)–(7.42) for kH ∈ {0.1, 1, 10, ∞}
    a, k = 0.8, 0.5
    for H in (0.2, 2.0, 20.0, np.inf):
        q = ch07.wave_energy(a, k, H, G, 1025.0, method="quad")
        c = ch07.wave_energy(a, k, H, G, 1025.0)
        assert q["Ek"] == pytest.approx(c["Ek"], rel=1e-8) and q["Ep"] == pytest.approx(c["Ep"], rel=1e-8)
        assert q["Ek"] == pytest.approx(q["Ep"], rel=1e-8)  # equipartition (7.41)
    assert ch07.wave_energy(1.0)["E"] == pytest.approx(0.5 * 1000 * G, rel=1e-15)  # (7.42): 4905 J/m² (D13 check)
    assert float(W.wave_energy_density(0.4, 1025.0, G)) == pytest.approx(0.5 * 1025 * G * 0.16, rel=1e-15)
    assert float(W.wave_energy_density(0.4, drho=3.0, g=G)) == pytest.approx(0.5 * 3.0 * G * 0.16, rel=1e-15)
    with pytest.raises(ValueError):
        ch07.wave_energy(method="nope")


def test_wave_energy_V2_derivation():  # V2 D12 (★★) steps 3–10 and D13: E_k = E_p = ½ρg⟨η²⟩
    z, xx = sp.symbols("z x", real=True)
    a, k, H, g, rho, om, lam = sp.symbols("a k H g rho omega lambda", positive=True)
    xavg = sp.integrate(a ** 2 * sp.cos(2 * sp.pi * xx / lam) ** 2, (xx, 0, lam)) / lam  # step 4
    assert sp.simplify(xavg - a ** 2 / 2) == 0
    Ic = sp.integrate(sp.cosh(k * (z + H)) ** 2, (z, -H, 0))
    Is = sp.integrate(sp.sinh(k * (z + H)) ** 2, (z, -H, 0))
    assert sym_zero(Ic - (H / 2 + sp.sinh(2 * k * H) / (4 * k)))  # step 5
    assert sym_zero(Is - (-H / 2 + sp.sinh(2 * k * H) / (4 * k)))  # step 6
    eb = a ** 2 / 2
    Ek = rho * om ** 2 / (2 * sp.sinh(k * H) ** 2) * (eb * Ic + eb * Is)  # step 3 after the x-averages
    assert sym_zero(Ek - rho * om ** 2 * eb / (2 * sp.sinh(k * H) ** 2) * sp.sinh(2 * k * H) / (2 * k))  # step 7
    assert sym_zero(Ek - rho * om ** 2 * eb * sp.cosh(k * H) / (2 * k * sp.sinh(k * H)))  # step 8
    assert sym_zero(Ek.subs(om, sp.sqrt(g * k * sp.tanh(k * H))) - rho * g * eb / 2)  # steps 9–10 (7.39)
    eta = sp.Symbol("eta", real=True)  # D13 steps 2–3: ∫_{−H}^{η} z dz − ∫_{−H}^{0} z dz = η²/2 (for η < 0 too)
    assert sp.simplify(sp.integrate(z, (z, -H, eta)) - sp.integrate(z, (z, -H, 0)) - eta ** 2 / 2) == 0


def test_energy_flux_V1_quadrature_equals_closed_form():  # V1 (7.43) → (7.44), (7.71): F = E c_g
    a, k = 0.5, 0.35
    for H in (1.0, 8.0, 40.0, np.inf):  # kH = 0.35, 2.8, 14, ∞
        Fq = ch07.energy_flux(a, k, H, G, 1025.0, method="quad")
        Fc = ch07.energy_flux(a, k, H, G, 1025.0)
        assert Fq == pytest.approx(Fc, rel=1e-8)
        assert Fc == pytest.approx(float(W.wave_energy_density(a, 1025.0, G)) * float(W.group_velocity(k, H, G)),
                                   rel=1e-13)  # (7.71)
    om = TWO_PI / 8.0  # a deep-water swell, T = 8 s, a = 1 m: F = ρg²a²T/(8π) (closed form of (7.44) with c_g = g/2ω)
    k8 = om ** 2 / G
    assert ch07.energy_flux(1.0, k8, np.inf, G, 1000.0) == pytest.approx(1000 * G ** 2 * 8.0 / (8 * np.pi), rel=1e-13)


def test_energy_flux_V2_derivation():  # V2 D14 (★★) steps 5–10, and D21 step 10: (7.44) ≡ E c_g
    z, t = sp.symbols("z t", real=True)
    a, k, H, g, rho, om = sp.symbols("a k H g rho omega", positive=True)
    C = sp.cosh(k * (z + H)) / sp.sinh(k * H)
    pu = (rho * a * om ** 2 / k * C * sp.cos(om * t)) * (a * om * C * sp.cos(om * t))  # step 5 at x = 0
    assert sym_zero(pu - rho * a ** 2 * om ** 3 / k * C ** 2 * sp.cos(om * t) ** 2)
    T = 2 * sp.pi / om
    tavg = sp.simplify(sp.integrate(sp.cos(om * t) ** 2, (t, 0, T)) / T)
    assert tavg == sp.Rational(1, 2)  # step 6
    F = rho * a ** 2 * om ** 3 / (2 * k * sp.sinh(k * H) ** 2) * (H / 2 + sp.sinh(2 * k * H) / (4 * k))  # step 7
    assert sym_zero(F - rho * a ** 2 * om ** 3 / (2 * k) * sp.integrate(C ** 2, (z, -H, 0)))
    step8 = rho * a ** 2 * om ** 3 * sp.sinh(2 * k * H) / (8 * k ** 2 * sp.sinh(k * H) ** 2) * (
        1 + 2 * k * H / sp.sinh(2 * k * H))
    assert sym_zero(F - step8)
    Fg = F.subs(om, sp.sqrt(g * k * sp.tanh(k * H)))
    oms = sp.sqrt(g * k * sp.tanh(k * H))
    assert sym_zero(Fg - rho * a ** 2 * g * oms / (4 * k) * (1 + 2 * k * H / sp.sinh(2 * k * H)))  # step 9
    c = oms / k
    assert sym_zero(Fg - (rho * g * a ** 2 / 2) * (c / 2) * (1 + 2 * k * H / sp.sinh(2 * k * H)))  # step 10 (7.44)
    assert sym_zero(Fg - (rho * g * a ** 2 / 2) * sp.diff(oms, k))  # = E dω/dk (7.71)


def test_energy_flux_V4_packet_energy_conserved_and_carried_at_cg():  # V4 ∫η² conserved; V3 its centroid → c_g
    x = np.linspace(-400, 1200, 2 ** 13, endpoint=False)
    k0 = 1.0
    cg = float(W.group_velocity(k0, np.inf, G))
    errs, widths = [], [10.0, 20.0, 40.0]
    for sx in widths:
        eta0 = np.exp(-x ** 2 / (2 * sx ** 2)) * np.cos(k0 * x)
        ev = W.linear_evolve(eta0, x, np.array([0.0, 100.0, 200.0]), direction=+1, g=G)
        E = np.sum(ev ** 2, axis=1)
        assert np.ptp(E) / E[0] < 1e-12  # (7.42)-energy of a one-way linear field is exactly conserved (Parseval)
        env2 = W.envelope(ev) ** 2
        cen = np.sum(x * env2, axis=1) / np.sum(env2, axis=1)
        v = (cen[2] - cen[0]) / 200.0
        errs.append(abs(v / cg - 1))
    assert errs[-1] < 2e-4  # energy moves at c_g ((7.71)) …
    assert abs(observed_order(1.0 / np.array(widths), errs) - 2.0) < ORDER_TOL  # … with an O(δk²) narrow-band error


# =====================================================================================================================
# C07 — capillary–gravity waves (7.53)–(7.60): R09, N49–N57, D15, D16; V5 Wikipedia "Capillary wave"
# =====================================================================================================================
def test_capillary_V2_derivation_tension_condition():  # V2 D15 (★★): (7.53) → (7.54) → (7.55) → (7.56), (7.57)
    xx, z, t = sp.symbols("x z t", real=True)
    a, k, H, g, rho, sig, om = sp.symbols("a k H g rho sigma omega", positive=True)
    eta = sp.Function("eta")(xx)
    kappa = (sp.diff(eta, xx, 2) / (1 + sp.diff(eta, xx) ** 2) ** sp.Rational(3, 2))  # (7.53) step 4
    eps = sp.Symbol("epsilon", positive=True)
    lin = sp.series(kappa.subs(eta, eps * sp.cos(k * xx)).doit(), eps, 0, 2).removeO()
    assert sp.simplify(lin - eps * sp.diff(sp.cos(k * xx), xx, 2)) == 0  # small slope: 1/R ≅ η_xx
    th = k * xx - om * t
    eta_w = a * sp.cos(th)
    assert sp.simplify(sp.diff(eta_w, xx, 2) + k ** 2 * eta_w) == 0  # step 8
    phi = a * om / k * sp.cosh(k * (z + H)) / sp.sinh(k * H) * sp.sin(th)  # step 7: D05 unchanged
    dyn = sp.simplify((sp.diff(phi, t).subs(z, 0) - sig / rho * sp.diff(eta_w, xx, 2) + g * eta_w) / sp.cos(th))  # (7.55)
    om2 = sp.solve(sp.Eq(dyn, 0), om ** 2)
    assert len(om2) == 1 and sym_zero(om2[0] - k * (g + sig * k ** 2 / rho) * sp.tanh(k * H))  # (7.56)
    c756 = sp.sqrt(k * (g + sig * k ** 2 / rho) * sp.tanh(k * H)) / k
    assert sym_zero(c756 - sp.sqrt((g / k + sig * k / rho) * sp.tanh(k * H)))  # (7.57)
    assert sp.simplify(ch07.surface_wave_sympy()["residuals"]["dynamic_7_55"]) == 0
    p = ch07.capillary_surface_pressure(-1e-3 * (TWO_PI / 0.01) ** 2, 0.0727)  # step 3: under a crest p > p_a
    assert float(p) > 0 and float(p) == pytest.approx(0.0727 * 1e-3 * (TWO_PI / 0.01) ** 2, rel=1e-14)
    assert round(float(p), 1) == 28.7  # D15 check (a = 1 mm, λ = 1 cm)


def test_capillary_minimum_V2_derivation():  # V2 D16 (★★): d(c²)/dk = 0 → k_m, λ_m, c_min (7.58); c_g = c there
    k, g, rho, sig = sp.symbols("k g rho sigma", positive=True)
    c2 = g / k + sig * k / rho  # step 1
    km = sp.solve(sp.diff(c2, k), k)  # step 3
    assert len(km) == 1 and sp.simplify(km[0] - sp.sqrt(rho * g / sig)) == 0  # step 4
    assert sp.simplify(g / km[0] - sig * km[0] / rho) == 0  # both restoring terms equal there
    assert sp.simplify(2 * sp.pi / km[0] - 2 * sp.pi * sp.sqrt(sig / (rho * g))) == 0  # step 5 (7.58)
    cmin2 = sp.simplify(c2.subs(k, km[0]))
    assert sp.simplify(cmin2 - 2 * sp.sqrt(g * sig / rho)) == 0  # step 6
    assert sp.simplify(sp.sqrt(cmin2) - (4 * g * sig / rho) ** sp.Rational(1, 4)) == 0  # step 7
    assert sp.simplify(sp.diff(c2, k, 2) - 2 * g / k ** 3) == 0
    om = sp.sqrt(k * (g + sig * k ** 2 / rho))
    assert sp.simplify((sp.diff(om, k) - om / k).subs(k, km[0])) == 0  # tangent through the origin: c_g = c


def test_capillary_minimum_V1_numerical_minimisation():  # V1 (7.58); N54; Ex. 7.10's c_g,min condition
    for sig, rho in ((0.0727, 1000.0), (0.030, 800.0), (0.0485, 13500.0)):
        m = W.capillary_minimum(sig, rho, G)
        r = minimize_scalar(lambda lk: float(W.phase_speed(np.exp(lk), np.inf, G, sig, rho)),
                            bounds=(0.0, 12.0), method="bounded", options={"xatol": 1e-12})
        assert m["k_m"] == pytest.approx(np.exp(r.x), rel=1e-5) and m["c_min"] == pytest.approx(r.fun, rel=1e-12)
        assert float(W.group_velocity(m["k_m"], np.inf, G, sig, rho)) == pytest.approx(m["c_min"], rel=1e-13)
        mg = W.min_group_velocity(sig, rho, G)
        r2 = minimize_scalar(lambda lk: float(W.group_velocity(np.exp(lk), np.inf, G, sig, rho)),
                             bounds=(0.0, 12.0), method="bounded", options={"xatol": 1e-12})
        assert mg["cg_min"] == pytest.approx(r2.fun, rel=1e-12) and mg["k"] == pytest.approx(np.exp(r2.x), rel=1e-5)
        assert sig * mg["k"] ** 2 / (rho * G) == pytest.approx(2 / np.sqrt(3) - 1, rel=1e-14)
    x = sp.Symbol("x", positive=True)  # a-D33: c_g² ∝ x^{−1/2}(1 + 3x)²/(1 + x) is stationary at 3x² + 6x − 1 = 0
    sol = sp.solve(sp.diff(x ** sp.Rational(-1, 2) * (1 + 3 * x) ** 2 / (1 + x), x), x)
    assert any(sp.simplify(s - (2 / sp.sqrt(3) - 1)) == 0 for s in sol)


@needs_ref
def test_capillary_V5_air_water_minimum():  # V5 Wikipedia "Capillary wave": 0.23 m/s at 1.7 cm (IAPWS σ, 20 °C)
    ref = ref_json()["capillary_wave"]
    sig = float(ch01.surface_tension_water(293.15))
    rho = float(ch01.water_density(293.15))  # 998.2 kg/m³ at 20 °C (Kell 1975, ch01)
    m = W.capillary_minimum(sig, rho, G)
    assert round(m["c_min"], 2) == ref["c_min_m_s"] and round(100 * m["lam_m"], 1) == ref["lambda_m_cm"]
    k = np.logspace(-1, 4, 80)  # the interface form with ρ′ = 0 is our (7.56) (deep water)
    wiki = np.sqrt(np.abs(k) * (G + sig / rho * k ** 2))
    assert maxrel(W.omega_capillary_gravity(k, np.inf, sig, rho, G), wiki) < 1e-14


def test_capillary_V7_limits():  # V7 (7.56)–(7.60): σ → 0 gives (7.28); g → 0 gives √(2πσ/ρλ); tension only speeds up
    k = np.logspace(-2, 4, 60)
    for H in (0.05, 2.0, np.inf):
        assert maxrel(W.omega_capillary_gravity(k, H, 0.0, 1000.0, G), W.omega_gravity(k, H, G)) < 1e-15
        assert np.all(W.phase_speed(k, H, G, 0.0727) > W.phase_speed(k, H, G, 0.0))
    lam = np.logspace(-4, -2, 10)
    assert maxrel(W.phase_speed(TWO_PI / lam, np.inf, 0.0, 0.0727, 1000.0), np.sqrt(TWO_PI * 0.0727 / (1000 * lam))) < 1e-14
    cg = W.group_velocity(TWO_PI / lam, np.inf, 0.0, 0.0727, 1000.0)
    assert maxrel(cg, 1.5 * W.phase_speed(TWO_PI / lam, np.inf, 0.0, 0.0727, 1000.0)) < 1e-14  # c_g = 3c/2 (Ex. 7.9)


def test_curvature_V1_circle_and_laplace_jump_parity():  # V1 (7.53), (7.54), R09 (ch04 laplace_jump_from_balance)
    R = 0.7
    xs = np.linspace(-0.5, 0.5, 21)
    eta_x = xs / np.sqrt(R ** 2 - xs ** 2)  # lower half of a circle z = −√(R² − x²): concave up, centre above
    eta_xx = R ** 2 / (R ** 2 - xs ** 2) ** 1.5
    assert maxrel(ch07.curvature(eta_x, eta_xx), np.full_like(xs, 1 / R)) < 1e-14
    assert maxrel(ch07.curvature(0.0, eta_xx, linear=True), eta_xx) < 1e-15
    p = ch07.capillary_surface_pressure(eta_xx, 0.0727, eta_x=eta_x, p_a=100.0)
    dp = IF.laplace_jump_from_balance(0.0727, R, 1e15)
    assert maxrel(100.0 - p, np.full_like(xs, dp)) < 1e-12  # p_a − p = σ/R: higher pressure on the centre's side


def test_capillary_state_V1_explainer_numbers():  # V1 Part C 2.16 (E3 state)
    s = ch07.capillary_state(0.02, 0.0727, 1000.0, np.inf, G)
    k = TWO_PI / 0.02
    assert s["tension_ratio"] == pytest.approx(0.0727 * k ** 2 / (1000 * G), rel=1e-14)
    assert s["c"] == pytest.approx(np.sqrt(s["gravity_term"] + s["tension_term"]), rel=1e-14)
    assert s["regime"] == "crossover" and ch07.capillary_state(0.003)["regime"] == "capillary"
    assert ch07.capillary_state(0.2)["regime"] == "gravity"
    assert s["p_crest_per_a"] == pytest.approx(0.0727 * k ** 2) and s["cg_over_c"] == pytest.approx(s["cg"] / s["c"])
    assert s["c_min"] == pytest.approx(W.capillary_minimum()["c_min"]) and s["cg_min"] < s["c_min"]
    z = ch07.capillary_state(0.02, 0.0, 1000.0)
    assert np.isnan(z["c_min"]) and z["regime"] == "gravity"


# =====================================================================================================================
# C08 — standing waves and seiches (7.61)–(7.65): N58–N63, N177, N178, D17, D18
# =====================================================================================================================
def test_standing_wave_V1_sum_of_opposite_waves():  # V1 (7.61)–(7.63): N58, N59, Fig. 7.11
    a, k, H = 0.15, 0.7, 2.5
    om = float(W.omega_gravity(k, H, G))
    x, z, t = RNG.uniform(0, 15, 40), RNG.uniform(-H, 0, 40), RNG.uniform(0, 10, 40)
    s = ch07.standing_wave_fields(x, z, t, a, k, H, G)
    r, l = ch07.wave_fields(x, z, t, a, k, H, G, direction=+1), ch07.wave_fields(x, z, t, a, k, H, G, direction=-1)
    for c in ("eta", "psi", "u", "w", "phi", "p_prime"):
        assert np.max(np.abs(s[c] - (r[c] + l[c]))) < 1e-13 * max(1.0, np.max(np.abs(s[c])))
    assert maxrel(s["psi"], 2 * a * om / k * np.sinh(k * (z + H)) / np.sinh(k * H) * np.sin(k * x) * np.sin(om * t)) < 1e-13
    assert maxrel(s["u"], 2 * a * om * np.cosh(k * (z + H)) / np.sinh(k * H) * np.sin(k * x) * np.sin(om * t)) < 1e-13
    xn = (np.pi / 2 + np.pi * np.arange(4)) / k  # fixed nodes of η
    assert np.max(np.abs(ch07.standing_wave_fields(xn, 0.0, np.linspace(0, 9, 4)[:, None], a, k, H, G)["eta"])) < 1e-15
    tt = np.linspace(0, TWO_PI / om, 256, endpoint=False)[:, None]  # Ex. 7.8: zero mean energy flux at every x
    ss = ch07.standing_wave_fields(np.array([0.3, 1.1, 2.9]), -0.6, tt, a, k, H, G)
    assert np.max(np.abs(np.mean(ss["p_prime"] * ss["u"], axis=0))) < 1e-10 * 1000 * G * a * a * om


def test_standing_wave_V2_derivation():  # V2 D17 (★★): sum-to-product, ψ₂ with a minus sign, (7.62), (7.63)
    x, z, t = sp.symbols("x z t", real=True)
    a, k, H, om = sp.symbols("a k H omega", positive=True)
    eta = a * sp.cos(k * x - om * t) + a * sp.cos(k * x + om * t)  # steps 1–2
    assert sp.simplify(sp.expand_trig(eta - 2 * a * sp.cos(k * x) * sp.cos(om * t))) == 0  # step 3
    S = sp.sinh(k * (z + H)) / sp.sinh(k * H)
    psi1 = a * om / k * S * sp.cos(k * x - om * t)  # (7.37)
    psi2 = (a * om / k * S * sp.cos(k * x - om * t)).subs(om, -om)  # step 5: ω → −ω flips the prefactor
    assert sp.simplify(psi2 + a * om / k * S * sp.cos(k * x + om * t)) == 0
    psi = psi1 + psi2  # step 6
    assert sp.simplify(sp.expand_trig(psi - 2 * a * om / k * S * sp.sin(k * x) * sp.sin(om * t))) == 0  # step 7 (7.62)
    u = sp.diff(2 * a * om / k * S * sp.sin(k * x) * sp.sin(om * t), z)  # step 8 (7.63)
    assert sp.simplify(u - 2 * a * om * sp.cosh(k * (z + H)) / sp.sinh(k * H) * sp.sin(k * x) * sp.sin(om * t)) == 0
    wrong = psi1 - psi2  # forgetting the minus sign gives cos kx instead of sin kx
    assert sp.simplify(sp.expand_trig(wrong - 2 * a * om / k * S * sp.sin(k * x) * sp.sin(om * t))) != 0


def test_seiche_V1_walls_modes_and_dispersion():  # V1 (7.64)–(7.65), D18, Fig. 7.12
    L, H = 1.5, 0.2
    n = np.arange(5)
    s = ch07.seiche_modes(L, H, n, G)
    k = (n + 1) * np.pi / L
    assert maxrel(s["k"], k) < 1e-15 and maxrel(s["lam"], 2 * L / (n + 1)) < 1e-15
    assert maxrel(s["omega"], W.omega_gravity(k, H, G)) < 1e-14  # (7.65) = (7.28) at k = (n + 1)π/L
    assert round(float(s["T"][0]), 2) == 2.20  # D18 check: bathtub
    for kk in s["k"]:
        u = ch07.standing_wave_fields(np.array([0.0, L]), np.linspace(-H, 0, 5)[:, None], 0.37, 0.01, kk, H, G)["u"]
        assert np.max(np.abs(u)) < 1e-15 + 1e-13 * 0.01  # u = 0 on both walls at every depth
    b = ch07.basin_modes(3.0, 2.0, 0.5, m=1, n=0, g=G)
    assert b["omega"] == pytest.approx(float(ch07.seiche_modes(3.0, 0.5, 0, G)["omega"]), rel=1e-14)
    b2 = ch07.basin_modes(3.0, 2.0, 0.5, m=2, n=3, g=G)
    assert b2["k"] == pytest.approx(np.hypot(2 * np.pi / 3.0, 3 * np.pi / 2.0), rel=1e-15)


def _seiche_steklov_fd(L, H, nx, g=G, nmodes=4):
    """Our independent check of (7.65): the sloshing eigenproblem ∇²φ = 0 in [0, L] × [−H, 0], φ_x = 0 on the walls,
    φ_z = 0 on the bottom, gφ_z = ω²φ on z = 0, by 2nd-order finite differences (ghost-point Neumann conditions) and a
    Schur complement onto the surface nodes. Returns the lowest non-zero ω [rad/s]."""
    nz = max(4, int(round(nx * H / L)))
    dx, dz = L / nx, H / nz
    Nx, Nz = nx + 1, nz + 1
    idx = lambda i, j: j * Nx + i  # noqa: E731
    A = spm.lil_matrix((Nx * Nz, Nx * Nz))
    for j in range(Nz):
        for i in range(Nx):
            p = idx(i, j)
            iw, ie = (i - 1 if i > 0 else 1), (i + 1 if i < nx else nx - 1)
            A[p, idx(iw, j)] += 1 / dx ** 2
            A[p, idx(ie, j)] += 1 / dx ** 2
            A[p, p] += -2 / dx ** 2 - 2 / dz ** 2
            if j == 0:
                A[p, idx(i, 1)] += 2 / dz ** 2
            elif j == nz:
                A[p, idx(i, nz - 1)] += 2 / dz ** 2  # the ghost adds (2/dz)λφ, λ = ω²/g
            else:
                A[p, idx(i, j - 1)] += 1 / dz ** 2
                A[p, idx(i, j + 1)] += 1 / dz ** 2
    A = A.tocsr()
    s = np.array([idx(i, nz) for i in range(Nx)])
    ii = np.setdiff1d(np.arange(Nx * Nz), s)
    Sc = A[s][:, s].toarray() - A[s][:, ii] @ spla.spsolve(A[ii][:, ii].tocsc(), A[ii][:, s].toarray())
    lam = np.sort(np.real(np.linalg.eigvals(-Sc * dz / 2)))
    return np.sqrt(g * lam[lam > 1e-9][:nmodes])


def test_seiche_V3_finite_difference_sloshing_eigenproblem():  # V3 (7.65) at finite depth, 2nd-order FD (ours)
    L, H = 10.0, 2.0
    exact = np.array([float(ch07.seiche_modes(L, H, n, G)["omega"]) for n in range(4)])
    nxs = [40, 80, 160]
    errs = [np.max(np.abs(_seiche_steklov_fd(L, H, nx) / exact - 1)) for nx in nxs]
    assert errs[-1] < 3e-4
    assert abs(observed_order([L / n for n in nxs], errs) - 2.0) < ORDER_TOL


def test_seiche_state_V7_shallow_limit():  # V7 Part C 2.19 (E4): T → 2L/((n + 1)√(gH)) as H/L → 0
    errs = []
    for HL in (0.02, 0.01, 0.005):
        s = ch07.seiche_state(40e3, 40e3 * HL, 0, G)
        assert s["T_min"] == pytest.approx(s["T"] / 60) and s["kH"] == pytest.approx(np.pi * HL)
        assert s["T_shallow"] == pytest.approx(2 * 40e3 / np.sqrt(G * 40e3 * HL), rel=1e-14)
        errs.append(abs(s["shallow_error"]))
    assert abs(observed_order([0.02, 0.01, 0.005], errs) - 2.0) < ORDER_TOL  # error ∝ (kH)²/6
    assert errs[-1] == pytest.approx((np.pi * 0.005) ** 2 / 6, rel=0.01)


# =====================================================================================================================
# C09 — group velocity (7.66)–(7.71): N64–N73, N179–N182, D19–D21
# =====================================================================================================================
def test_beats_V2_derivation():  # V2 D19 (★): sum-to-product with ½Δω t; the printed ½Δω x cannot move
    x, t = sp.symbols("x t", real=True)
    k1, k2, w1, w2, a = sp.symbols("k1 k2 omega1 omega2 a", positive=True)
    eta = a * sp.cos(k1 * x - w1 * t) + a * sp.cos(k2 * x - w2 * t)
    k, w, dk, dw = (k1 + k2) / 2, (w1 + w2) / 2, k2 - k1, w2 - w1
    beats = 2 * a * sp.cos(dk * x / 2 - dw * t / 2) * sp.cos(k * x - w * t)  # (7.66) with t
    assert sp.simplify(sp.expand(sp.expand_trig(sp.expand(eta - beats)))) == 0 or sym_zero(eta - beats)
    printed = 2 * a * sp.cos(dk * x / 2 - dw * x / 2) * sp.cos(k * x - w * t)
    assert not sym_zero(eta - printed)
    assert sp.diff(dk * x / 2 - dw * x / 2, t) == 0  # the printed envelope is frozen in time


def test_beat_wave_V1_equals_sum_printed_fails_and_chord_to_tangent():  # V1 (7.66)–(7.67), N65; V3 chord → tangent
    k1, k2, a, H = 0.9, 1.1, 0.3, 4.0
    x = np.linspace(-60, 60, 801)
    for t in (0.0, 3.3, 17.0):
        b = W.beat_wave(x, t, k1, k2, a, H, G)
        om1, om2 = float(W.omega_gravity(k1, H, G)), float(W.omega_gravity(k2, H, G))
        exact = a * np.cos(k1 * x - om1 * t) + a * np.cos(k2 * x - om2 * t)
        assert np.max(np.abs(b["eta"] - exact)) < 1e-13
        if t > 0:
            bp = W.beat_wave(x, t, k1, k2, a, H, G, printed=True)
            assert np.max(np.abs(bp["eta"] - exact)) > 0.1 * a
    d = W.beat_wave(0.0, 0.0, 0.9, 1.1, 1.0, np.inf, G)  # D19 check: Δω/Δk = 1.568 vs ½√(g/k) = 1.566
    assert round(d["cg_finite"], 3) == 1.568 and round(0.5 * np.sqrt(G), 3) == 1.566
    hs, errs = [0.2, 0.1, 0.05, 0.025], []
    for dk in hs:
        b = W.beat_wave(0.0, 0.0, 1.0 - dk, 1.0 + dk, 1.0, 2.0, G)
        errs.append(abs(b["cg_finite"] - float(W.group_velocity(1.0, 2.0, G))))
    assert abs(observed_order(hs, errs) - 2.0) < ORDER_TOL  # the symmetric chord approaches dω/dk at O(Δk²)


def test_packet_envelope_V2_derivation():  # V2 D20 (★★★): steps 4–12 for a Gaussian spectrum
    kap, xi, x = sp.symbols("kappa xi x", real=True)
    t, s, w2 = sp.symbols("t s omega2", positive=True)
    k0, w0, cg, tt = sp.symbols("k0 omega0 c_g t_", real=True)
    ph = (k0 + kap) * x - (w0 + cg * kap + w2 * kap ** 2 / 2) * tt  # steps 3, 5: k = k0 + κ, ω Taylor-expanded
    assert sp.expand(ph - ((k0 * x - w0 * tt) + kap * (x - cg * tt) - w2 * kap ** 2 * tt / 2)) == 0  # step 6
    A = sp.exp(-kap ** 2 / (2 * s ** 2))  # step 2: a narrow Gaussian spectrum
    a0 = sp.integrate(A * sp.exp(sp.I * kap * x), (kap, -sp.oo, sp.oo))  # step 4: the initial envelope a(x)
    env1 = sp.integrate(A * sp.exp(sp.I * kap * xi), (kap, -sp.oo, sp.oo))  # steps 8–9
    assert sp.simplify(env1 - a0.subs(x, xi)) == 0  # the envelope is a(x − c_g t): (7.68)
    assert sp.simplify(a0 - s * sp.sqrt(2 * sp.pi) * sp.exp(-s ** 2 * x ** 2 / 2)) == 0  # real for symmetric A (step 10)
    # step 11: with ω″ kept the modulus is m2 = 2πs²/√(1 + q²) exp(−s²ξ²/(1 + q²)), q = ω″s²t — the closed form of
    # ∫A e^{iκξ − iω″κ²t/2}dκ (a complex Gaussian integral); check it solves the defining ODE and initial value
    q = w2 * s ** 2 * t
    env2 = s * sp.sqrt(2 * sp.pi) / sp.sqrt(1 + sp.I * q) * sp.exp(-s ** 2 * xi ** 2 / (2 * (1 + sp.I * q)))
    # ∂env/∂t = −(iω″/2)∫κ²A e^{…}dκ = (iω″/2)∂²env/∂ξ² (Schrödinger form of the dropped term)
    assert sp.simplify(sp.diff(env2, t) - sp.I * w2 / 2 * sp.diff(env2, xi, 2)) == 0
    assert sp.simplify(env2.subs(t, 0) - env1.subs(xi, xi)) == 0
    m2 = sp.simplify(sp.expand_complex(env2 * sp.conjugate(env2)))
    assert sp.simplify(sp.diff(m2, xi).subs(xi, 0)) == 0  # the peak stays at ξ = 0, x = c_g t (step 12)
    assert sp.simplify(m2.subs(xi, 0) / m2.subs({xi: 0, t: 0}) - 1 / sp.sqrt(1 + q ** 2)) == 0  # spreading, step 7
    xg = np.linspace(-60, 60, 7)  # the coded closed form (gaussian_packet, order 2) is this envelope
    om_q = lambda kk: 1.0 + 0.8 * (kk - 1.0) + 0.3 * (kk - 1.0) ** 2  # noqa: E731
    gp = W.gaussian_packet(xg, 25.0, 1.0, 1.0, 6.0, omega_fn=om_q, order=2)
    sx = 6.0
    ref = np.array([complex(env2.subs({s: 1 / sx, w2: 0.6, t: 25.0, xi: v - 0.8 * 25.0})) for v in xg]) / (
        np.sqrt(TWO_PI) / sx)
    assert np.max(np.abs(gp["envelope"] - np.abs(ref))) < 1e-12


def test_group_velocity_V2_derivation():  # V2 D21 (★★): dω/dk of (7.28) = (7.69); limits (7.70); F = E c_g (7.71)
    k, H, g, xx = sp.symbols("k H g x", positive=True)
    om = sp.sqrt(g * k * sp.tanh(k * H))
    lhs = sp.diff(om ** 2, k)  # steps 1–2
    assert sym_zero(lhs - (g * sp.tanh(k * H) + g * k * H / sp.cosh(k * H) ** 2))
    cg = lhs / (2 * om)  # step 3
    c = om / k
    assert sym_zero(cg - sp.diff(om, k))
    assert sym_zero(cg / c - sp.Rational(1, 2) * (1 + k * H / (sp.cosh(k * H) ** 2 * sp.tanh(k * H))))  # steps 4–6
    assert sym_zero(1 / (sp.cosh(xx) ** 2 * sp.tanh(xx)) - 2 / sp.sinh(2 * xx))  # step 7 identity
    assert sym_zero(cg - c / 2 * (1 + 2 * k * H / sp.sinh(2 * k * H)))  # (7.69)
    assert sp.limit(2 * xx / sp.sinh(2 * xx), xx, sp.oo) == 0 and sp.limit(2 * xx / sp.sinh(2 * xx), xx, 0) == 1  # (7.70)
    sig, rho = sp.symbols("sigma rho", positive=True)  # the capillary–gravity form coded in group_velocity
    omc = sp.sqrt(k * (g + sig * k ** 2 / rho) * sp.tanh(k * H))
    coded = omc / k / 2 * ((g + 3 * sig * k ** 2 / rho) / (g + sig * k ** 2 / rho) + 2 * k * H / sp.sinh(2 * k * H))
    assert sym_zero(sp.diff(omc, k) - coded)


def test_group_velocity_V1_complex_step_parity():  # V1 (7.67), (7.69): analytic c_g = numerical dω/dk to 1e-12
    k = np.logspace(-2, 3, 50)
    for H, sig in ((np.inf, 0.0), (1.5, 0.0), (0.02, 0.0727), (np.inf, 0.0727)):
        an = W.group_velocity(k, H, G, sig, 1000.0)
        cs = W.group_velocity_numeric(None, k, H=H, g=G, sigma=sig, rho=1000.0)
        assert maxrel(an, cs) < 1e-12
    beta, l = 2e-11, 1e-6  # a Rossby-type relation ω = −βk/(k² + l²) (dispersion given as data)
    kk = np.array([3e-7, 1e-6, 4e-6])
    ros = W.group_velocity_numeric(lambda q: -beta * q / (q ** 2 + l ** 2), kk)
    assert maxrel(ros, -beta * (l ** 2 - kk ** 2) / (kk ** 2 + l ** 2) ** 2) < 1e-12
    m, N = 0.8, 0.01
    iw = W.group_velocity_numeric(lambda q: W.internal_wave_omega(q, m, N), 0.6)
    assert float(iw) == pytest.approx(N * m ** 2 / (0.6 ** 2 + m ** 2) ** 1.5, rel=1e-12)  # (7.145) x-part
    assert float(W.group_velocity_numeric(lambda q: np.abs(q) * 2.0, 0.5, method="central")) == pytest.approx(2.0)
    with pytest.raises(TypeError):
        W.group_velocity_numeric(lambda q: float(np.real(q)), 0.5, method="complex")


def test_group_velocity_numeric_V3_central_differences_fourth_order():  # V3 the central fallback converges at order 4
    ex = float(W.group_velocity(0.7, 2.0, G))
    hs = [0.08, 0.04, 0.02, 0.01]
    errs = [abs(float(W.group_velocity_numeric(lambda q: np.sqrt(G * q * np.tanh(2.0 * q)), 0.7, h=h, method="central"))
                - ex) for h in hs]
    assert abs(observed_order(hs, errs) - 4.0) < ORDER_TOL


def test_group_velocity_V7_limits():  # V7 (7.70): c_g/c → ½ deep, → 1 shallow, 3/2 capillary; c_g = c at k_m
    assert float(W.group_velocity(1.0, 30.0, G) / W.phase_speed(1.0, 30.0, G)) == pytest.approx(0.5, abs=1e-12)
    assert float(W.group_velocity(1.0, np.inf, G) / W.phase_speed(1.0, np.inf, G)) == 0.5
    r = float(W.group_velocity(1e-4, 1.0, G) / W.phase_speed(1e-4, 1.0, G))
    assert r == pytest.approx(1.0, abs=1e-8)
    assert float(W.group_velocity(1e6, np.inf, G, 0.0727) / W.phase_speed(1e6, np.inf, G, 0.0727)) == pytest.approx(
        1.5, abs=1e-6)  # ½(g + 3sk²)/(g + sk²) → 3/2 as sk² ≫ g
    km = W.capillary_minimum(0.0727, 1000.0, G)["k_m"]
    assert float(W.group_velocity(km, np.inf, G, 0.0727)) == pytest.approx(float(W.phase_speed(km, np.inf, G, 0.0727)),
                                                                          rel=1e-14)


def test_gaussian_packet_V1_envelope_rides_at_cg():  # V1 (7.68) order 1 exactly; order 2 = FFT for quadratic ω
    x = np.linspace(-300, 900, 2 ** 13, endpoint=False)
    g1 = W.gaussian_packet(x, 0.0, 1.0, 1.0, 20.0, order=1, g=G)
    for T in (50.0, 200.0):
        gT = W.gaussian_packet(x, T, 1.0, 1.0, 20.0, order=1, g=G)
        assert gT["cg"] == pytest.approx(0.5 * np.sqrt(G), rel=1e-12)
        assert np.max(np.abs(gT["envelope"] - 1.0 * np.exp(-(x - gT["cg"] * T) ** 2 / 800.0))) < 1e-14  # rigid shift
        assert np.max(np.abs(gT["eta"] - gT["envelope"] * np.cos(x - gT["omega0"] * T))) < 1e-12  # carrier at c
    assert np.max(np.abs(g1["eta"] - np.exp(-x ** 2 / 800.0) * np.cos(x))) < 1e-14
    om_q = lambda kk: 1.0 + 0.8 * (kk - 1.0) + 0.3 * (kk - 1.0) ** 2  # noqa: E731
    for T in (50.0, 200.0):
        gp = W.gaussian_packet(x, T, 1.0, 1.0, 20.0, omega_fn=om_q, order=2)
        le = W.linear_evolve(np.exp(-x ** 2 / 800.0) * np.exp(1j * x), x, T, omega_fn=om_q)
        assert np.max(np.abs(gp["eta"] - le)) < 1e-12 and gp["omega2"] == pytest.approx(0.6, rel=1e-8)


def test_group_velocity_V1_negative_k_is_signed_derivative():  # V1 review M1: dω/dk of ω(|k|) is sgn(k)·(7.69)
    # (7.28)/(7.56) depend on |k| only, so dω/dk = sgn(k) ω′(|k|): a left-going wave (k < 0) has c_g < 0, same size.
    k = np.array([1e-4, 0.05, 0.3, 1.0, 7.0, 400.0])
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # kH up to 800: no overflow warning from the complex tanh either
        for H, sig in ((np.inf, 0.0), (2.0, 0.0), (0.05, 0.0), (np.inf, 0.0727), (0.02, 0.0727)):
            pos = W.group_velocity(k, H, G, sig, 1000.0)
            neg = W.group_velocity(-k, H, G, sig, 1000.0)
            assert np.all(pos > 0) and np.array_equal(neg, -pos)  # odd in k
            for kk in (k, -k):
                an = W.group_velocity(kk, H, G, sig, 1000.0)
                cs = W.group_velocity_numeric(None, kk, H=H, g=G, sigma=sig, rho=1000.0)  # complex step
                assert np.all(np.sign(cs) == np.sign(kk)) and maxrel(cs, an) < 1e-12
                ce = W.group_velocity_numeric(lambda q: np.sqrt(np.abs(q) * (G + sig / 1000.0 * q ** 2)
                                                                * (np.tanh(np.abs(q) * H) if np.isfinite(H) else 1.0)),
                                              kk, method="central")  # an independent real ω(|k|), differenced
                assert maxrel(ce, an) < 1e-9
    # closed forms at k = −1: deep water −½√(g/|k|) (7.70); finite depth −(c/2)(1 + 2|k|H/sinh 2|k|H) (7.69)
    assert float(W.group_velocity_numeric(k=-1.0, g=G)) == pytest.approx(-0.5 * np.sqrt(G), rel=1e-13)
    assert float(W.group_velocity(-1.0, np.inf, G)) == pytest.approx(-0.5 * np.sqrt(G), rel=1e-13)
    c2 = np.sqrt(G * np.tanh(2.0))
    ex = -0.5 * c2 * (1.0 + 4.0 / np.sinh(4.0))
    assert float(W.group_velocity_numeric(k=-1.0, H=2.0, g=G)) == pytest.approx(ex, rel=1e-13)
    assert float(W.group_velocity(-1.0, 2.0, G)) == pytest.approx(ex, rel=1e-13)
    # shallow limit (7.70): c_g → −√(gH) for a long left-going wave
    assert float(W.group_velocity_numeric(k=-1e-5, H=1.0, g=G)) == pytest.approx(-np.sqrt(G), rel=1e-9)


def test_gaussian_packet_V1_negative_k0_is_the_mirror_image():  # V1 review M1: k₀ < 0 is the x → −x mirror, c_g < 0
    x = np.linspace(-600, 600, 2 ** 12, endpoint=False)
    for H in (np.inf, 5.0):
        for order in (1, 2):
            gp = W.gaussian_packet(x, 40.0, 1.0, 0.5, 20.0, order=order, H=H, g=G)
            gm = W.gaussian_packet(-x, 40.0, 1.0, -0.5, 20.0, order=order, H=H, g=G)
            # η₀ = a e^{−x²/2σ²} cos k₀x is even in x and in k₀; ω(|k|) is even ⇒ η(x, t; −k₀) = η(−x, t; k₀)
            assert gm["cg"] == pytest.approx(-gp["cg"], rel=1e-13) and gm["cg"] < 0
            assert gm["cg"] == pytest.approx(float(W.group_velocity(-0.5, H, G)), rel=1e-12)
            assert gm["omega0"] == gp["omega0"] and gm["omega2"] == pytest.approx(gp["omega2"], rel=1e-8)
            assert np.max(np.abs(gm["eta"] - gp["eta"])) < 1e-12
    t = 200.0  # deep water, order 1: the envelope rides rigidly to the LEFT at c_g = −½√(g/|k₀|)
    g1 = W.gaussian_packet(x, t, 1.0, -0.5, 20.0, order=1, g=G)
    assert g1["cg"] == pytest.approx(-0.5 * np.sqrt(G / 0.5), rel=1e-12) and abs(g1["cg"]) < 1e3
    assert np.max(np.abs(g1["envelope"] - np.exp(-(x - g1["cg"] * t) ** 2 / 800.0))) < 1e-14
    # order 2 is exact for quadratic dispersion: ω(k) = Ω(|k|), Ω(q) = 1 + 0.8(q − 1) + 0.3(q − 1)², packet at k₀ = −1
    # (its spectrum lies at k < 0, where ω(k) = Ω(−k): c_g = −0.8, ω″ = 0.6) vs the FFT evolution with Ω(|k|)
    om_q = lambda q: 1.0 + 0.8 * (q - 1.0) + 0.3 * (q - 1.0) ** 2  # noqa: E731
    xx = np.linspace(-900, 300, 2 ** 13, endpoint=False)
    for T in (50.0, 200.0):
        gp = W.gaussian_packet(xx, T, 1.0, -1.0, 20.0, omega_fn=lambda q: om_q(-q), order=2)
        le = W.linear_evolve(np.exp(-xx ** 2 / 800.0) * np.exp(-1j * xx), xx, T, omega_fn=om_q)
        assert gp["cg"] == pytest.approx(-0.8, rel=1e-12) and gp["omega2"] == pytest.approx(0.6, rel=1e-8)
        assert np.max(np.abs(gp["eta"] - le)) < 1e-12


def test_linear_evolve_V7_nondispersive_translation():  # V7: ω = ck translates any shape rigidly (shallow limit)
    L, N, c = 50.0, 500, 1.7
    x = np.linspace(0, L, N, endpoint=False)
    eta0 = np.exp(-((x - 10) / 2.0) ** 2) - 0.5 * np.exp(-((x - 20) / 1.0) ** 2)
    dx = L / N
    t = 37 * dx / c
    ev = W.linear_evolve(eta0, x, t, omega_fn=lambda k: c * k, direction=+1)
    assert np.max(np.abs(ev - np.roll(eta0, 37))) < 1e-13


def test_pond_ripples_V3_cosine_sum_converges_to_the_cauchy_poisson_integral():  # V1 t = 0 hump; V3 midpoint rule
    w = 0.01
    xp = np.linspace(-0.2, 0.2, 201)
    assert np.max(np.abs(ch07.pond_ripples(xp, 0.0, width=w, n_modes=512) - np.exp(-xp ** 2 / (2 * w ** 2)))) < 1e-8
    xs, tt = np.array([0.0, 0.03, 0.08]), 0.05

    def integrand(k, x):
        return w * np.sqrt(TWO_PI) * np.exp(-0.5 * (k * w) ** 2) * np.cos(k * x) * np.cos(
            float(W.omega_capillary_gravity(k, np.inf, 0.0727, 1000.0, G)) * tt) / np.pi

    ref = np.array([quad(integrand, 0, 6 / w, args=(v,), limit=2000, epsabs=1e-13, epsrel=1e-12)[0] for v in xs])
    ns = [256, 512, 1024]
    errs = [np.max(np.abs(ch07.pond_ripples(xs, tt, width=w, n_modes=n, g=G) - ref)) for n in ns]
    assert abs(observed_order([1 / n for n in ns], errs) - 2.0) < ORDER_TOL
    two = ch07.pond_ripples(xs, np.array([0.0, tt]), width=w, n_modes=256, g=G)
    assert two.shape == (2, 3) and np.max(np.abs(two[1] - ch07.pond_ripples(xs, tt, width=w, n_modes=256, g=G))) < 1e-15


def test_envelope_V1_hilbert_modulus():  # V1: |η + iH[η]| of a narrow-band packet is its Gaussian envelope
    x = np.linspace(-200, 200, 4096, endpoint=False)
    env = 0.7 * np.exp(-x ** 2 / (2 * 20 ** 2))
    assert np.max(np.abs(W.envelope(env * np.cos(3.0 * x)) - env)) < 1e-12


def test_packet_state_V1_explainer_numbers():  # V1 Part C 2.20 (E5 state)
    s = ch07.packet_state(0.05, np.inf, G, a=1.5, distance=2e6)
    assert s["ratio"] == pytest.approx(0.5, rel=1e-13) and s["regime"] == "cg<c"
    assert s["F"] == pytest.approx(s["E"] * s["cg"]) and s["E"] == pytest.approx(0.5 * 1000 * G * 2.25)
    assert s["arrival_h"] == pytest.approx(2e6 / s["cg"] / 3600) and s["t_crest_cross"] == pytest.approx(
        10 * s["lam"] / (s["c"] - s["cg"]))
    assert ch07.packet_state(1e-4, 1.0, G)["regime"] == "cg=c"
    assert ch07.packet_state(3000.0, np.inf, G, sigma=0.0727)["regime"] == "cg>c"


def test_viscous_decay_V2_derivation():  # V2 a-D68 (demoted, N72): ε = 2νS_ijS_ij of (7.47) → a₀e^{−2νk²t}
    x, z, t = sp.symbols("x z t", real=True)
    a, k, g, nu, rho = sp.symbols("a k g nu rho", positive=True)
    om = sp.sqrt(g * k)
    th = k * x - om * t
    u, w = a * om * sp.exp(k * z) * sp.cos(th), a * om * sp.exp(k * z) * sp.sin(th)  # (7.47)
    Sxx, Szz, Sxz = sp.diff(u, x), sp.diff(w, z), (sp.diff(u, z) + sp.diff(w, x)) / 2
    SS = sp.simplify(Sxx ** 2 + Szz ** 2 + 2 * Sxz ** 2)
    assert sp.simplify(SS - 2 * a ** 2 * om ** 2 * k ** 2 * sp.exp(2 * k * z)) == 0
    D = rho * sp.integrate(2 * nu * SS, (z, -sp.oo, 0))  # dissipation per unit area
    E = rho * g * a ** 2 / 2
    dadt = sp.simplify(-D / sp.diff(E, a))  # dE/dt = −D
    assert sp.simplify(dadt + 2 * nu * k ** 2 * a) == 0
    assert float(W.viscous_decay(0.3, 2.0, 1.3e-6, 1e4)) == pytest.approx(0.3 * np.exp(-2 * 1.3e-6 * 4 * 1e4), rel=1e-15)


# =====================================================================================================================
# C10 — kinematic wave theory and rays (7.72)–(7.79): N47, N48, N74–N80, N174, N175, N184, D22–D24
# =====================================================================================================================
def test_crest_conservation_V2_derivation():  # V2 D22 (★★): (7.73) → (7.74) → (7.75)
    x, t, c = sp.symbols("x t c", real=True)
    theta = sp.Function("theta")(x, t)
    Om = sp.Function("Omega")
    k, om = sp.diff(theta, x), -sp.diff(theta, t)  # (7.73)
    assert sp.simplify(sp.diff(k, t) + sp.diff(om, x)) == 0  # (7.74): mixed partials commute
    wl = Om(k)  # step 5: local dispersion relation
    kk = sp.Symbol("kk")
    assert sp.simplify(sp.diff(wl, x) - sp.diff(Om(kk), kk).subs(kk, k) * sp.diff(k, x)) == 0  # step 6 (chain rule)
    F = sp.Function("F")
    kc = F(x - c * t)
    assert sp.simplify(sp.diff(kc, t) + c * sp.diff(kc, x)) == 0  # step 8: constant along dx/dt = c


def test_frequency_along_rays_V2_derivation():  # V2 D23 (★★): (7.76)–(7.78) → (7.79) and dk/dt = −∂ω/∂x
    x, t, X = sp.symbols("x t X", real=True)
    g = sp.Symbol("g", positive=True)
    kf = sp.Function("k")(x, t)
    Hf = sp.Function("H")
    K = sp.Symbol("K", positive=True)
    Om = sp.sqrt(g * K * sp.tanh(K * Hf(X)))  # the local (7.28) with H = H(x)
    om = Om.subs({K: kf, X: x})  # (7.76)
    cg = sp.diff(Om, K).subs({K: kf, X: x})  # (7.77)
    assert sp.simplify(sp.diff(om, t) - cg * sp.diff(kf, t)) == 0  # steps 3–4 (7.78): at fixed x only k changes
    r774 = sp.diff(kf, t) + sp.diff(om, x)
    assert sp.simplify((sp.diff(om, t) + cg * sp.diff(om, x)) - cg * r774) == 0  # steps 5–6: (7.79) = c_g × (7.74)
    om_x_k = sp.diff(Om, X).subs({K: kf, X: x})  # step 8: (∂ω/∂x)_k
    assert sp.simplify(r774 - (sp.diff(kf, t) + cg * sp.diff(kf, x) + om_x_k)) == 0


def test_snell_V2_derivation():  # V2 D24 (★★): ω(|k|, H(x)) ⇒ k_y conserved; |k| sin α = const; α → 0 at the shore
    kx, ky, xx, yy = sp.symbols("k_x k_y x y", real=True)
    g, om, H, k0, a0 = sp.symbols("g omega H k0 alpha0", positive=True)
    Hf = sp.Function("H")
    w = sp.sqrt(g * sp.sqrt(kx ** 2 + ky ** 2) * sp.tanh(sp.sqrt(kx ** 2 + ky ** 2) * Hf(xx)))
    assert sp.diff(w, yy) == 0  # steps 2–3: dk_y/dt = −∂ω/∂y = 0
    kmag = sp.Symbol("kmag", positive=True)
    alpha = sp.asin(k0 * sp.sin(a0) / kmag)  # steps 6–8
    assert sp.simplify(kmag * sp.sin(alpha) - k0 * sp.sin(a0)) == 0
    ks = om / sp.sqrt(g * H)  # step 9: shallow |k| = ω/√(gH) → ∞ as H → 0
    assert sp.limit(alpha.subs(kmag, ks), H, 0, "+") == 0


def test_ray_trace_V1_homogeneous_straight_at_cg():  # V1: constant depth — straight rays at c_g, k constant
    H = 6.0
    k0 = np.array([0.08, 0.05])
    r = W.ray_trace(None, x0=(0.0, 0.0), k0=k0, t_span=(0.0, 300.0), H_fn=lambda X: H, n_out=31)
    km = np.linalg.norm(k0)
    cg = float(W.group_velocity(km, H, G))
    exp_x = cg * np.outer(k0 / km, r["t"])
    assert np.max(np.abs(r["x"] - exp_x)) < 1e-8 * np.max(np.abs(exp_x))
    assert np.max(np.abs(r["k"] - k0[:, None])) < 1e-12 and r["omega_drift"] < 1e-10  # solve_ivp rtol = 1e-10
    r1 = W.ray_trace(lambda kv, xv: float(W.omega_gravity(np.linalg.norm(kv), H, G)), x0=0.0, k0=0.1,
                     t_span=(0.0, 50.0), n_out=5)
    assert r1["x"][0, -1] == pytest.approx(float(W.group_velocity(0.1, H, G)) * 50.0, rel=1e-8)
    with pytest.raises(ValueError):
        W.ray_trace(None, x0=(0.0, 0.0), k0=(0.1, 0.0))


def test_ray_trace_V4_frequency_conserved_while_k_grows():  # V4 (7.79): ω constant along a ray over a 1:50 beach
    slope, T = 1 / 50.0, 8.0
    om = TWO_PI / T
    x0 = 1000.0
    k0 = float(W.wavenumber_from_omega(om, slope * x0, G))
    al0 = np.radians(30.0)
    r = W.ray_trace(None, x0=(x0, 0.0), k0=(-k0 * np.cos(al0), k0 * np.sin(al0)), t_span=(0.0, 3000.0),
                    H_fn=lambda X: slope * X[0], n_out=400)
    assert r["status"] == 1  # stopped by the shore event
    assert r["omega_drift"] < 1e-8  # (D23 check)
    km = np.hypot(*r["k"])
    assert km[-1] > 10 * km[0]  # the wave shortens as the water shallows
    assert np.ptp(r["k"][1]) < 1e-10 * abs(r["k"][1][0])  # k_y conserved (Snell)
    alpha = np.arcsin(r["k"][1] / km)
    assert np.all(np.diff(alpha) < 1e-12) and alpha[-1] < np.radians(2.0)  # crests turn parallel to the shore


def test_ray_trace_V1_snell_closed_form_parity():  # V1 D24: ray_trace ≡ snell_ray_plane_beach (quadrature of tan α)
    slope, T, x0, al0 = 1 / 40.0, 9.0, 800.0, np.radians(35.0)
    om = TWO_PI / T
    k0 = float(W.wavenumber_from_omega(om, slope * x0, G))
    r = W.ray_trace(None, x0=(x0, 0.0), k0=(-k0 * np.cos(al0), k0 * np.sin(al0)), t_span=(0.0, 400.0),
                    H_fn=lambda X: slope * X[0], n_out=41)
    xs = r["x"][0][1:]
    sn = ch07.snell_ray_plane_beach(xs, al0, x0, T, slope, G)
    assert np.max(np.abs(r["x"][1][1:] - sn["y"])) < 1e-6 * np.max(sn["y"])
    assert np.ptp(sn["snell"]) < 1e-12 * sn["l"] and sn["l"] == pytest.approx(k0 * np.sin(al0), rel=1e-14)
    assert np.max(np.abs(sn["alpha"] - np.arctan2(r["k"][1][1:], -r["k"][0][1:]))) < 1e-8


def test_refraction_state_V1_worked_example():  # V1 Part C 2.24 (E6): T = 8 s, 30° at 20 m → 11.3° at 2 m (D24 check)
    s = ch07.refraction_state(8.0, np.radians(30.0), 20.0, 2.0, G)
    assert round(s["alpha_deg"], 1) == 11.3
    assert s["snell"] == pytest.approx(s["k0"] * 0.5, rel=1e-14) and s["c"] == pytest.approx(s["omega"] / s["k"])
    assert s["lam"] == pytest.approx(TWO_PI / s["k"]) and s["cg"] == pytest.approx(float(W.group_velocity(s["k"], 2.0, G)))


def test_local_wavenumber_V3_fourth_order_and_crest_conservation():  # V3 (7.73) order 4; V1 (7.74), (7.79) residuals
    th = lambda xx, tt: 0.5 * xx + 0.002 * xx ** 2 - 2 * tt + 0.01 * np.sin(0.3 * xx - 0.1 * tt)  # noqa: E731
    kx = lambda xx, tt: 0.5 + 0.004 * xx + 0.003 * np.cos(0.3 * xx - 0.1 * tt)  # noqa: E731
    xs = np.array([1.0, 3.0, 7.0])
    hs, errs = [0.4, 0.2, 0.1, 0.05], []
    for h in hs:
        k, _ = ch07.local_wavenumber_frequency(th, xs, 0.5, h=h)
        errs.append(np.max(np.abs(k - kx(xs, 0.5))))
    assert abs(observed_order(hs, errs) - 4.0) < ORDER_TOL
    r = ch07.crest_conservation_residual(th, xs, 0.5)
    assert np.max(np.abs(r["r774"])) < 1e-9  # (7.74) for any smooth phase (D22 check: 1e-9; nested-FD round-off
    # bound ≈ 2.25·ε·|θ|/(h·h_t) = 1.8e-9 at x = 7 with h = h_t = 1e-3)
    Hx = lambda x_: 12.0 - 0.01 * x_  # noqa: E731  a steady train over a gentle slope: θ = ∫k(x)dx − ωt
    om = TWO_PI / 7.0
    kofx = lambda x_: float(W.wavenumber_from_omega(om, Hx(x_), G))  # noqa: E731

    def theta(xx, tt):
        xx = np.atleast_1d(np.asarray(xx, float))
        vals = np.array([quad(kofx, 0.0, v, epsabs=1e-13, epsrel=1e-13)[0] for v in xx.ravel()]).reshape(xx.shape)
        return vals - om * np.asarray(tt, float)

    rr = ch07.crest_conservation_residual(theta, np.array([100.0, 300.0]), 2.0, H_fn=Hx, h=0.5, ht=0.5)
    assert np.max(np.abs(rr["r_dispersion"])) < 1e-8 * om  # the local dispersion relation holds along the train
    assert np.max(np.abs(rr["r779"])) < 1e-10 * om  # (7.79): ω constant
    assert np.max(np.abs(rr["r775"])) > 1e-6  # (7.75) is for homogeneous media only: k is not carried unchanged
    hom = ch07.crest_conservation_residual(lambda xx, tt: 0.3 * xx - float(W.omega_gravity(0.3, 5.0, G)) * tt,
                                           np.array([1.0, 2.0]), 0.3, omega_fn=lambda k: W.omega_gravity(k, 5.0, G))
    assert np.max(np.abs(hom["r775"])) < 1e-10 and np.max(np.abs(hom["r_dispersion"])) < 1e-10


# =====================================================================================================================
# C11 — hydraulic jump (7.80)–(7.81), nonlinear steepening, KdV (7.87) and the solitary wave (7.88): R10, N81–N85,
# N93–N97, N185, N186, N189, D25, D26
# =====================================================================================================================
def test_hydraulic_jump_V2_derivation_belanger():  # V2 D25 (★★): CV momentum → (7.80) → quadratic → (7.81)
    rho, g, H1, H2, Q, zz, Hh = sp.symbols("rho g H1 H2 Q z H", positive=True)
    Fr, r = sp.symbols("Fr1 r", positive=True)
    face = sp.integrate(rho * g * (Hh - zz), (zz, 0, Hh))  # step 3
    assert sp.simplify(face - rho * g * Hh ** 2 / 2) == 0
    mom = rho * Q * (Q / H2 - Q / H1) - (rho * g * H1 ** 2 / 2 - rho * g * H2 ** 2 / 2)  # step 4 (u = Q/H, step 2)
    e780 = Q ** 2 * (1 / H2 - 1 / H1) - g * (H1 ** 2 - H2 ** 2) / 2  # step 5 (7.80)
    assert sp.simplify(mom / rho - e780) == 0
    assert sp.simplify(e780 - (H1 - H2) * (Q ** 2 / (H1 * H2) - g * (H1 + H2) / 2)) == 0  # step 6: common factor
    Q2 = sp.solve(Q ** 2 / (H1 * H2) - g * (H1 + H2) / 2, Q ** 2)[0]  # step 7
    assert sp.simplify(Q2 - g * H1 * H2 * (H1 + H2) / 2) == 0
    quadr = sp.simplify((Q2 / (g * H1 ** 3)).subs(H2, r * H1))  # step 8: Fr₁² = r(1 + r)/2
    assert sp.simplify(quadr - r * (1 + r) / 2) == 0
    rr = sp.Symbol("rr", real=True)
    roots = sp.solve(rr ** 2 + rr - 2 * Fr ** 2, rr)  # step 9
    assert sp.simplify(roots[0] * roots[1] + 2 * Fr ** 2) == 0  # Vieta: product −2Fr₁² < 0
    pos = [q for q in roots if sp.simplify(q - (-1 + sp.sqrt(1 + 8 * Fr ** 2)) / 2) == 0]
    assert len(roots) == 2 and len(pos) == 1  # step 10 (7.81): the other root is negative
    assert all(sp.simplify(q + (1 + sp.sqrt(1 + 8 * Fr ** 2)) / 2) == 0 for q in roots if q not in pos)


def test_hydraulic_jump_V2_derivation_energy_loss():  # V2 D26 (★★): E₂ − E₁ = −g(H₂ − H₁)³/(4H₁H₂)
    g, H1, H2 = sp.symbols("g H1 H2", positive=True)
    Q2 = g * H1 * H2 * (H1 + H2) / 2
    dE = Q2 / 2 * (1 / H2 ** 2 - 1 / H1 ** 2) + g * (H2 - H1)  # steps 1–2
    assert sp.simplify(1 / H2 ** 2 - 1 / H1 ** 2 - (H1 - H2) * (H1 + H2) / (H1 ** 2 * H2 ** 2)) == 0  # step 3
    assert sp.simplify(Q2 / 2 * (H1 - H2) * (H1 + H2) / (H1 ** 2 * H2 ** 2) - g * (H1 + H2) ** 2 * (H1 - H2) / (
        4 * H1 * H2)) == 0  # step 4
    assert sp.simplify(dE - g * (H2 - H1) / (4 * H1 * H2) * (4 * H1 * H2 - (H1 + H2) ** 2)) == 0  # step 5
    assert sp.expand(4 * H1 * H2 - (H1 + H2) ** 2 + (H2 - H1) ** 2) == 0  # step 6
    assert sp.simplify(dE + g * (H2 - H1) ** 3 / (4 * H1 * H2)) == 0  # step 7


def test_hydraulic_jump_V1_momentum_mass_and_numbers():  # V1 (7.80), (7.81), D25/D26 checks
    for Fr1 in (1.0, 1.4, 3.0, 7.5):
        j = ch07.hydraulic_jump(0.25, Fr1=Fr1, g=G)
        assert float(ch07.jump_momentum_residual(0.25, j["H2"], j["Q"], G)) == pytest.approx(0.0, abs=1e-13)
        assert j["u1"] * 0.25 == pytest.approx(j["u2"] * j["H2"], rel=1e-14)  # mass
        assert j["Fr2"] <= 1.0 + 1e-15 and j["Fr1"] == pytest.approx(Fr1, rel=1e-14)
        assert j["E2"] - j["E1"] == pytest.approx(j["dE"], rel=1e-10, abs=1e-14)  # the closed loss = direct energies
    assert ch07.hydraulic_jump(0.25, Fr1=1.0, g=G)["H2"] == pytest.approx(0.25, rel=1e-15)  # Fr₁ = 1: no jump
    j = ch07.hydraulic_jump(0.1, Fr1=3.0, g=G)  # D25/D26 checks
    assert round(j["H2"], 3) == 0.377 and round(j["dE"], 3) == -1.385 and round(j["head_loss"], 3) == 0.141
    assert ch07.hydraulic_jump(0.1, u1=j["u1"], g=G)["H2"] == pytest.approx(j["H2"], rel=1e-15)
    with pytest.raises(ValueError):
        ch07.hydraulic_jump(0.1, u1=1.0, Fr1=2.0)


def test_hydraulic_jump_V4_second_law():  # V4: mechanical energy never increases across a physical jump
    Fr = np.linspace(1.0, 10.0, 91)
    dE = np.array([ch07.hydraulic_jump(0.3, Fr1=f, g=G)["dE"] for f in Fr])
    assert np.all(dE <= 0.0) and np.all(np.diff(dE) < 0)  # a loss that grows with Fr₁
    sub = [ch07.hydraulic_jump(0.3, Fr1=f, g=G) for f in (0.3, 0.6, 0.9)]  # a "jump down" would create energy
    assert all(s["dE"] > 0 and not s["allowed"] and s["H2"] < 0.3 for s in sub)


@needs_ref
def test_hydraulic_jump_V1_form_belanger_published():  # V1 form (Wikipedia "Hydraulic jumps in rectangular channels")
    ref = ref_json()["belanger"]
    Fr1, y1, y2 = sp.symbols("Fr1 y1 y2", positive=True)
    ratio = sp.lambdify(Fr1, sp.sympify(ref["depth_ratio"], locals={"Fr1": Fr1}))
    loss = sp.lambdify((y1, y2), sp.sympify(ref["head_loss"], locals={"y1": y1, "y2": y2}))
    for f in (1.2, 2.5, 6.0):
        j = ch07.hydraulic_jump(0.4, Fr1=f, g=G)
        assert j["ratio"] == pytest.approx(ratio(f), rel=1e-14)
        assert j["head_loss"] == pytest.approx(loss(0.4, j["H2"]), rel=1e-12)


def test_jump_state_V1_budget_and_moving_bore_parity():  # V1 Part C 2.27 (E7): CV budget, ch04 bore speed
    s = ch07.jump_state(0.2, 2.5, G, 1000.0)
    assert abs(s["residual"]) < 1e-10 * s["p_out"]
    assert s["mom_in"] - s["mom_out"] == pytest.approx(s["p_out"] - s["p_in"], rel=1e-12)
    assert s["power_loss"] == pytest.approx(1000 * G * s["Q"] * s["head_loss"], rel=1e-14)
    assert s["bore_speed"] == pytest.approx(float(ch04.bore_speed(0.2, s["H2"], G)), rel=1e-13)  # Fig. 7.20c frame
    assert s["flow_behind"] == pytest.approx(s["u1"] - s["u2"]) and s["allowed"]
    d = ch07.jump_state(0.1, 3.0, G)  # D25 check numbers: 882.9 − 234.1 = 697.9 − 49.05 N/m
    assert (round(d["mom_in"], 1), round(d["mom_out"], 1), round(d["p_out"], 1), round(d["p_in"], 2)) == (
        882.9, 234.1, 697.9, 49.05)


def test_kdv_V2_solitary_wave_residual():  # V2 Exercise 7.15 (a-D46): (7.88) solves (7.87) with c = c₀(1 + a/2H)
    r = ch07.kdv_residual_sympy()
    assert sp.simplify(r["residual"]) == 0 and all(sp.simplify(c) == 0 for c in r["residual_coefficients"])
    bad = ch07.kdv_residual_sympy(nonlinear_coefficient=3)  # the wrong 3 instead of 3/2
    assert sp.simplify(bad["residual"]) != 0
    xx, tt = sp.symbols("x t", real=True)  # independently, with sech directly (no tanh substitution)
    a, H, g = sp.symbols("a H g", positive=True)
    c0 = sp.sqrt(g * H)
    c = c0 * (1 + a / (2 * H))
    eta = a / sp.cosh(sp.sqrt(3 * a / (4 * H ** 3)) * (xx - c * tt)) ** 2
    kdv = sp.diff(eta, tt) + c0 * sp.diff(eta, xx) + sp.Rational(3, 2) * c0 / H * eta * sp.diff(eta, xx) + \
        c0 * H ** 2 / 6 * sp.diff(eta, xx, 3)
    f = sp.lambdify((xx, tt, a, H, g), kdv, "numpy")
    fx = sp.lambdify((xx, tt, a, H, g), c0 * sp.diff(eta, xx), "numpy")
    pts = RNG.uniform(-3, 3, (20, 2))
    for av, Hv in ((0.1, 1.0), (0.3, 2.0), (0.05, 0.5)):
        res = np.array([f(px, pt, av, Hv, G) for px, pt in pts])
        sc = np.max(np.abs([fx(px, pt, av, Hv, G) for px, pt in pts]))
        assert np.max(np.abs(res)) < 1e-12 * sc  # zero to round-off at 60 random points


def test_kdv_solve_V1_soliton_translates_at_c0_1_plus_a_over_2H():  # V1 (7.88) under our IF-RK4 spectral scheme
    H, a = 1.0, 0.2
    x = np.linspace(-40, 40, 256, endpoint=False)
    eta0 = ch07.solitary_wave(x, 0.0, a, H, G, x0=-20.0)
    c = ch07.solitary_wave_speed(a, H, G)
    assert c == pytest.approx(np.sqrt(G * H) * 1.1, rel=1e-15)
    T = 30.0 / c
    r = ch07.kdv_solve(eta0, x, [T], H, G)
    assert np.max(np.abs(r["eta"][-1] - ch07.solitary_wave(x, T, a, H, G, x0=-20.0))) < 1e-6 * a
    e_c = ch07.solitary_wave(x, 0.0, a, H, G)  # centred: its tails are 3e-14 at the periodic seam
    rhs = ch07.kdv_rhs(e_c, x, H, G)  # a travelling wave: η_t = −c η_x
    ex = np.real(np.fft.ifft(1j * TWO_PI * np.fft.fftfreq(x.size, x[1] - x[0]) * np.fft.fft(e_c)))
    assert np.max(np.abs(rhs + c * ex)) < 1e-9 * c * a
    wrong = ch07.solitary_wave_speed(a, H, G) * 1.0 - np.sqrt(G * H) * a / (2 * H)  # c₀ alone would lag behind
    assert abs(wrong - c) > 0.05


def test_kdv_solve_V3_time_order_four_and_spectral_space():  # V3 IF-RK4 in time (order 4), spectral in x
    H, a = 1.0, 0.2
    x = np.linspace(-40, 40, 256, endpoint=False)
    eta0 = ch07.solitary_wave(x, 0.0, a, H, G, x0=-20.0)
    ref = ch07.kdv_solve(eta0, x, [2.0], H, G, dt=0.0025)["eta"][-1]
    dts = [0.04, 0.02, 0.01]
    errs = [np.max(np.abs(ch07.kdv_solve(eta0, x, [2.0], H, G, dt=d)["eta"][-1] - ref)) for d in dts]
    assert abs(observed_order(dts, errs) - 4.0) < ORDER_TOL
    es = []
    for n in (64, 128):
        xx = np.linspace(-40, 40, n, endpoint=False)
        r = ch07.kdv_solve(ch07.solitary_wave(xx, 0.0, a, H, G, x0=-20.0), xx, [2.0], H, G, dt=0.005)
        es.append(np.max(np.abs(r["eta"][-1] - ch07.solitary_wave(xx, 2.0, a, H, G, x0=-20.0))))
    assert es[1] < es[0] / 300  # spectral: doubling N gains far more than any fixed order (≥ 10^2.5)


def test_kdv_solve_V4_invariants_of_a_splitting_hump():  # V4 mass exactly; momentum and Hamiltonian to O(dt⁴)
    H = 1.0
    x = np.linspace(-100, 100, 512, endpoint=False)
    eta0 = 0.3 * np.exp(-(x / 4.0) ** 2)
    ts = np.linspace(0, 20.0, 11)
    r = ch07.kdv_solve(eta0, x, ts, H, G)
    assert np.ptp(r["mass"]) / abs(r["mass"][0]) < 1e-13
    assert np.ptp(r["momentum"]) / r["momentum"][0] < 1e-8 and np.ptp(r["energy"]) / abs(r["energy"][0]) < 1e-8
    d1 = ch07.kdv_solve(eta0, x, ts, H, G, dt=0.01)
    d2 = ch07.kdv_solve(eta0, x, ts, H, G, dt=0.005)
    q = np.ptp(d1["momentum"]) / np.ptp(d2["momentum"])
    assert np.log2(q) > 3.0  # the drift is the integrator's truncation error (≈ dt⁴), not a leak
    inv = ch07.kdv_invariants(eta0, x)
    assert np.isnan(inv["energy"]) and inv["mass"] == pytest.approx(np.sum(eta0) * (x[1] - x[0]))
    assert r["cached"] is False and r["steps"] > 0


def test_kdv_linear_phase_speed_V1_taylor_of_7_29():  # V1 N94: c₀(1 − k²H²/6) = first two terms of (7.29); error ∝ (kH)⁴
    H = 2.0
    kH = np.array([0.05, 0.1, 0.2, 0.4])
    err = np.abs(ch07.kdv_linear_phase_speed(kH / H, H, G) / W.phase_speed(kH / H, H, G) - 1)
    assert abs(observed_order(kH, err) - 4.0) < ORDER_TOL
    kk, Hs = sp.symbols("k H", positive=True)
    ser = sp.series(sp.sqrt(sp.tanh(kk * Hs) / (kk * Hs)), kk, 0, 4).removeO()
    assert sp.simplify(ser - (1 - kk ** 2 * Hs ** 2 / 6)) == 0
    assert float(ch07.ursell_number(0.1, 10.0, 1.0)) == pytest.approx(10.0)


@needs_ref
def test_cnoidal_V1_form_and_soliton_limit():  # V1 form (Wikipedia "Cnoidal wave"); m → 1 with datum="trough" → (7.88)
    x = np.linspace(-30, 30, 1201)
    Hd, A = 1.0, 0.2
    for m in (0.5, 0.9, 0.999):
        cm = ch07.cnoidal_wave(x, 0.7, Hd, A, m, G)  # mean datum: compare with the published formulas
        K, E = ellipk(m), ellipe(m)
        eta2 = A / m * (1 - m - E / K)
        Delta = Hd * np.sqrt(4 * m * Hd / (3 * A))
        c = np.sqrt(G * Hd) * (1 + A / (m * Hd) * (1 - m / 2 - 1.5 * E / K))
        assert cm["trough"] == pytest.approx(eta2, rel=1e-13) and cm["Delta"] == pytest.approx(Delta, rel=1e-14)
        assert cm["c"] == pytest.approx(c, rel=1e-13) and cm["wavelength"] == pytest.approx(2 * Delta * K, rel=1e-14)
        lam = cm["wavelength"]  # zero mean over one wavelength
        xw = np.linspace(0, lam, 4001)
        e = ch07.cnoidal_wave(xw, 0.0, Hd, A, m, G)["eta"]
        assert abs(np.trapezoid(e, xw) / lam) < 1e-6 * A
    sw = ch07.solitary_wave(x, 0.0, A, Hd, G)
    tr = ch07.cnoidal_wave(x, 0.0, Hd, A, 1 - 1e-8, G, datum="trough")
    assert np.max(np.abs(tr["eta"] - sw)) < 1e-6 * A and tr["c"] == pytest.approx(ch07.solitary_wave_speed(A, Hd, G),
                                                                                  rel=1e-6)
    assert np.max(np.abs(ch07.cnoidal_wave(x, 0.0, Hd, A, 1 - 1e-8, G)["eta"] - sw)) > 1e-3  # mean datum: slow (1/K)
    Delta_s = Hd * np.sqrt(4 * Hd / (3 * A))  # the published solitary width = the book's (3a/4H³)^{−1/2}
    assert Delta_s == pytest.approx(1 / np.sqrt(3 * A / (4 * Hd ** 3)), rel=1e-15)
    with pytest.raises(ValueError):
        ch07.cnoidal_wave(x, 0.0, datum="bottom")


def test_simple_wave_V1_breaking_time_and_first_order_speed():  # V1 N81 (our Riemann simple wave, labelled)
    H, a, k = 2.0, 0.05, 0.2
    x = np.linspace(0, TWO_PI / k, 40001)
    eta0 = a * np.sin(k * x)
    s = ch07.simple_wave_evolve(eta0, x, np.array([0.0, 1.0]), H, G)
    c0 = np.sqrt(G * H)
    dcdx = 1.5 * np.sqrt(G / (H + eta0)) * a * k * np.cos(k * x)  # exact ∂c/∂x of c = 3√(g(H + η)) − 2c₀
    assert s["t_break"] == pytest.approx(1 / np.max(-dcdx), rel=1e-6)
    assert s["t_break"] == pytest.approx(2 * H / (3 * a * k * c0), rel=2 * a / H)  # leading order in a/H
    assert np.max(np.abs(s["x_points"][1] - (x + s["c"]))) < 1e-12 and np.all(s["eta"] == eta0)
    e = np.array([1e-4, 1e-3, 1e-2])
    for model in ("simple", "book"):  # both = c₀(1 + (3/2)η/H) + O(η²) — the KdV nonlinear coefficient
        cm = ch07.nonlinear_wavelet_speed(e, H, G, model)
        err = np.abs(cm - c0 * (1 + 1.5 * e / H))
        assert abs(observed_order(e, err) - 2.0) < ORDER_TOL
    assert float(ch07.nonlinear_wavelet_speed(0.1, H, G)) > float(ch07.nonlinear_wavelet_speed(-0.1, H, G))  # crests win
    with pytest.raises(ValueError):
        ch07.nonlinear_wavelet_speed(0.1, H, G, "nope")


# =====================================================================================================================
# C12 — Stokes waves and Stokes drift (7.82)–(7.86): N86–N92, N187, N188, D27
# =====================================================================================================================
def test_stokes_drift_V2_derivation():  # V2 D27 (★★): ⟨ξu_x + ζu_z⟩ = a²ωke^{2kz₀} (7.85); any depth (7.86)
    x0, z0, t = sp.symbols("x0 z0 t", real=True)
    a, k, om, H = sp.symbols("a k omega H", positive=True)
    phv = k * x0 - om * t
    xi, ze = -a * sp.exp(k * z0) * sp.sin(phv), a * sp.exp(k * z0) * sp.cos(phv)  # (7.46), step 4
    X, Z = sp.symbols("X Z", real=True)
    u = a * om * sp.exp(k * Z) * sp.cos(k * X - om * t)  # (7.47)
    ux, uz = sp.diff(u, X).subs({X: x0, Z: z0}), sp.diff(u, Z).subs({X: x0, Z: z0})
    assert sp.simplify(xi * ux - a ** 2 * om * k * sp.exp(2 * k * z0) * sp.sin(phv) ** 2) == 0  # step 5
    assert sp.simplify(ze * uz - a ** 2 * om * k * sp.exp(2 * k * z0) * sp.cos(phv) ** 2) == 0  # step 6
    assert sp.simplify(xi * ux + ze * uz - a ** 2 * om * k * sp.exp(2 * k * z0)) == 0  # step 7: time-independent
    T = 2 * sp.pi / om
    assert sp.simplify(sp.integrate(u.subs({X: x0, Z: z0}), (t, 0, T))) == 0  # step 3
    C = sp.cosh(k * (Z + H)) / sp.sinh(k * H)
    S = sp.sinh(k * (Z + H)) / sp.sinh(k * H)
    uH = a * om * C * sp.cos(k * X - om * t)
    xiH, zeH = -a * C.subs(Z, z0) * sp.sin(phv), a * S.subs(Z, z0) * sp.cos(phv)
    corr = xiH * sp.diff(uH, X).subs({X: x0, Z: z0}) + zeH * sp.diff(uH, Z).subs({X: x0, Z: z0})
    avg = sp.integrate(sp.expand(corr), (t, 0, T)) / T
    assert sym_zero(avg - a ** 2 * om * k * sp.cosh(2 * k * (z0 + H)) / (2 * sp.sinh(k * H) ** 2))  # (7.86), N91


def test_stokes_drift_V1_closed_form_and_limits():  # V1 (7.85)–(7.86); D27 check number
    a, k = 0.2, 0.5
    z = np.linspace(-6, 0, 13)
    om = float(W.omega_gravity(k, np.inf, G))
    assert maxrel(ch07.stokes_drift(z, a, k, np.inf, G), a ** 2 * om * k * np.exp(2 * k * z)) < 1e-14
    H = 5.0
    omH = float(W.omega_gravity(k, H, G))
    zz = np.linspace(-H, 0, 11)
    assert maxrel(ch07.stokes_drift(zz, a, k, H, G),
                  a ** 2 * omH * k * np.cosh(2 * k * (zz + H)) / (2 * np.sinh(k * H) ** 2)) < 1e-13
    assert float(ch07.stokes_drift(-1.0, a, k, 200.0, G)) == pytest.approx(float(ch07.stokes_drift(-1.0, a, k, np.inf, G)),
                                                                          rel=1e-13)
    om8 = TWO_PI / 8.0  # D27 check: a = 1 m, T = 8 s deep → 4.94 cm/s
    assert round(100 * float(ch07.stokes_drift(0.0, 1.0, om8 ** 2 / G, np.inf, G)), 2) == 4.94
    assert abs(ch07.eulerian_mean_u(-0.8, a, k, H, G)) < 1e-10 * a * omH  # N92: Eulerian mean 0 at a fixed point


def test_stokes_drift_V3_exact_path_lines_converge_at_order_ka():  # V3: numeric drift → (7.86) with error O(ka)
    k, H, z0 = 1.0, 3.0, -0.5
    ka = np.array([0.01, 0.02, 0.04, 0.08])
    rel = np.array([abs(ch07.stokes_drift_numeric(z0, e / k, k, H, G, periods=20) /
                        float(ch07.stokes_drift(z0, e / k, k, H, G)) - 1) for e in ka])
    assert rel[0] < 0.01
    assert abs(observed_order(ka, rel) - 1.0) < ORDER_TOL
    d1 = ch07.stokes_drift_numeric(z0, 0.04, k, H, G, periods=20, rtol=1e-10)
    d2 = ch07.stokes_drift_numeric(z0, 0.04, k, H, G, periods=20, rtol=5e-11)
    assert abs(d1 / d2 - 1) < 1e-7  # integrator-converged
    om = float(W.omega_gravity(k, H, G))
    T = TWO_PI / om
    tl = ch07.particle_path(0.0, z0, [0.0, 20 * T], 0.02, k, H, G, model="taylor1")  # (7.84a, b) model
    assert (tl["x"][-1] - tl["x"][0]) / (20 * T) == pytest.approx(float(ch07.stokes_drift(z0, 0.02, k, H, G)), rel=0.03)


@needs_ref
def test_stokes_drift_V5_deep_water_published():  # V5 Wikipedia "Stokes drift": ωka²e^{2kz}; ≈ 4 % at z = −λ/4
    ref = ref_json()["stokes_drift"]
    k, a = 0.3, 0.4
    lam = TWO_PI / k
    frac = float(ch07.stokes_drift(-lam / 4, a, k, np.inf, G) / ch07.stokes_drift(0.0, a, k, np.inf, G))
    assert round(100 * frac) == ref["quarter_wavelength_fraction_percent"]
    om, kk, aa, zz = sp.symbols("omega k a z")
    form = sp.lambdify((om, kk, aa, zz), sp.sympify(ref["deep_form"], locals={"omega": om, "k": kk, "a": aa, "z": zz}))
    z = np.linspace(-10, 0, 6)
    assert maxrel(ch07.stokes_drift(z, a, k, np.inf, G), form(np.sqrt(G * k), k, a, z)) < 1e-14


def test_stokes_expansion_V2_third_order_coefficients():  # V2 a-D40: (7.82) ½, 3/8 and (7.83) γ = 1; Ex. 7.2 literal γ = 3/8
    r = ch07.stokes_expansion_sympy()
    assert r["alpha"] == sp.Rational(1, 2) and r["gamma"] == 1
    t3 = r["third_order"]
    assert (t3["alpha"], t3["beta"], t3["delta"], t3["gamma"]) == (sp.Rational(1, 2), sp.Rational(-5, 8),
                                                                    sp.Rational(3, 8), 1)
    assert r["exercise_literal"]["gamma"] == sp.Rational(3, 8)  # the truncated set-up (recorded discrepancy)
    # independent expansion (ours): φ = Σ A_n e^{nz} sin nθ, η = ε cos θ + b₂ε² cos 2θ + b₃ε³ cos 3θ, ω = 1 + w₂ε²
    th, e, z = sp.symbols("theta epsilon z", real=True)
    a11, a13, a22, a33, b22, b33, w2, C2 = sp.symbols("a11 a13 a22 a33 b22 b33 w2 C2")
    from sympy.simplify.fu import TR8

    om = 1 + w2 * e ** 2
    phi = (e * a11 + e ** 3 * a13) * sp.exp(z) * sp.sin(th) + e ** 2 * a22 * sp.exp(2 * z) * sp.sin(2 * th) + \
        e ** 3 * a33 * sp.exp(3 * z) * sp.sin(3 * th)
    eta = e * sp.cos(th) + e ** 2 * b22 * sp.cos(2 * th) + e ** 3 * b33 * sp.cos(3 * th)
    kin = (sp.diff(phi, z) + om * sp.diff(eta, th) - sp.diff(eta, th) * sp.diff(phi, th)).subs(z, eta)
    dyn = (-om * sp.diff(phi, th) + (sp.diff(phi, th) ** 2 + sp.diff(phi, z) ** 2) / 2).subs(z, eta) + eta - e ** 2 * C2
    kin = sp.expand(sp.series(kin, e, 0, 4).removeO())
    dyn = sp.expand(sp.series(dyn, e, 0, 4).removeO())
    hm = lambda ex, fn, n: sp.expand(TR8(sp.expand(ex))).coeff(fn(n * th))  # noqa: E731
    sol = {a11: sp.solve(hm(kin.coeff(e, 1), sp.sin, 1), a11)[0]}
    d2 = dyn.coeff(e, 2).subs(sol)
    const = sp.expand(TR8(sp.expand(d2))).subs({sp.cos(2 * th): 0, sp.sin(2 * th): 0, sp.cos(th): 0, sp.sin(th): 0})
    sol.update(sp.solve([hm(kin.coeff(e, 2).subs(sol), sp.sin, 2), hm(d2, sp.cos, 2), const], [a22, b22, C2], dict=True)[0])
    k3, d3 = kin.coeff(e, 3).subs(sol), dyn.coeff(e, 3).subs(sol)
    s3 = sp.solve([hm(k3, sp.sin, 1), hm(d3, sp.cos, 1), hm(k3, sp.sin, 3), hm(d3, sp.cos, 3)], [a13, w2, a33, b33],
                  dict=True)[0]
    assert sol[b22] == sp.Rational(1, 2) and s3[b33] == sp.Rational(3, 8)  # (7.82)
    assert 2 * s3[w2] == 1  # ω² = gk(1 + (ka)²): (7.83)
    assert sol[a22] == 0 and s3[a33] == 0 and s3[a13] == sp.Rational(-1, 8)  # deep water: φ has one harmonic only
    assert s3[w2] + t3["beta"] == s3[a13]  # the coded β (φ = εω(1 + βε²)…) is our a₁₃ = ½ − 5/8 = −1/8


@needs_ref
def test_stokes_wave_V5_speed_and_limiting_steepness():  # V5 HandWiki "Stokes wave": 0.1410633; c = (1 + ½(ka)²)√(g/k)
    ref = ref_json()["stokes_limit"]
    assert ch07.STOKES_LIMIT_STEEPNESS == ref["H_over_lambda"]
    k = 0.4
    ka = np.array([0.02, 0.04, 0.08, 0.16])
    diff = np.abs(ch07.stokes_wave_speed(k, ka / k, G) / np.sqrt(G / k) - (1 + 0.5 * ka ** 2))  # published form
    assert abs(observed_order(ka, diff) - 4.0) < ORDER_TOL  # identical to O((ka)²), differs at O((ka)⁴)
    assert 0.5 * ref["H_over_lambda"] == pytest.approx(0.0705, abs=5e-5)  # max a ≈ half the crest-to-trough height


def test_stokes_wave_profile_V1_harmonics_and_permanence():  # V1 (7.82)–(7.83), Fig. 7.21 (N86, N87, N187)
    a, k = 0.1, 1.0
    x = np.linspace(0, TWO_PI, 721)
    c = float(ch07.stokes_wave_speed(k, a, G))
    assert c == pytest.approx(np.sqrt(G / k * (1 + (k * a) ** 2)), rel=1e-15)
    p3 = ch07.stokes_wave_profile(x, 0.0, a, k, G)
    assert maxrel(p3, a * np.cos(x) + 0.5 * k * a ** 2 * np.cos(2 * x) + 0.375 * k ** 2 * a ** 3 * np.cos(3 * x)) < 1e-14
    assert maxrel(ch07.stokes_wave_profile(x, 0.0, a, k, G, order=1), a * np.cos(x)) < 1e-15
    assert np.max(p3) > -np.min(p3)  # peaked crests, flat troughs
    assert maxrel(ch07.stokes_wave_profile(x + c * 3.3, 3.3, a, k, G), p3) < 1e-12  # permanent form, speed c
    assert abs(np.mean(p3[:-1])) < 1e-15  # the harmonics add no mean level


def test_dyed_line_V1_advances_by_the_stokes_drift():  # V1 Part C 2.12 (Fig. 7.22): each particle gains ū_L T per period
    a, k, H = 0.03, 1.0, 3.0
    z0s = np.array([-0.2, -0.8, -1.6])
    om = float(W.omega_gravity(k, H, G))
    T = TWO_PI / om
    d = ch07.dyed_line(z0s, np.array([0.0, 5 * T, 10 * T]), a, k, H, G)
    for j, z0 in enumerate(z0s):  # the same exact path lines (7.33) as particle_path, started on the line
        pp = ch07.particle_path(0.0, z0, np.array([0.0, 5 * T, 10 * T]), a, k, H, G, model="exact", start="mean")
        assert np.max(np.abs(d["x"][:, j] - pp["x"])) < 1e-8 * a and np.max(np.abs(d["z"][:, j] - pp["z"])) < 1e-8 * a
    adv = (d["x"][2] - d["x"][0]) / 10
    zbar = z0s - np.asarray(ch07.orbit_linear(0.0, z0s, 0.0, a, k, H, G)[1])  # a particle on the line at a crest
    assert maxrel(adv, ch07.stokes_drift(zbar, a, k, H, G) * T) < 3 * k * a  # sits aS above its mean depth; O(ka)
    assert np.all(np.diff(adv) < 0)  # the dyed line leans forward: the top drifts most
    d0 = ch07.dyed_line(z0s, 0.0, a, k, H, G)
    assert np.all(d0["x"] == 0.0) and np.all(d0["z"] == z0s)


# =====================================================================================================================
# C13 — interfacial waves (7.89)–(7.96): R11, N98–N108, N190–N192, D28
# =====================================================================================================================
def test_interface_V2_derivation():  # V2 D28 (★★): (7.89)–(7.94) → A = −B = iωa/k → ω = ε√(gk) (7.95)
    x, z, t = sp.symbols("x z t", real=True)
    a, k, g, om, r1, r2 = sp.symbols("a k g omega rho1 rho2", positive=True)
    A, B = sp.symbols("A B")
    E = sp.exp(sp.I * (k * x - om * t))
    phi1, phi2, zeta = A * sp.exp(-k * z) * E, B * sp.exp(k * z) * E, a * E  # steps 1–2
    assert sp.simplify(sp.diff(phi1, x, 2) + sp.diff(phi1, z, 2)) == 0 and sp.simplify(
        sp.diff(phi2, x, 2) + sp.diff(phi2, z, 2)) == 0  # (7.90)
    As = sp.solve(sp.simplify((sp.diff(phi1, z) - sp.diff(zeta, t)).subs(z, 0) / E), A)[0]  # step 4
    Bs = sp.solve(sp.simplify((sp.diff(phi2, z) - sp.diff(zeta, t)).subs(z, 0) / E), B)[0]  # step 5
    assert sp.simplify(As - sp.I * om * a / k) == 0 and sp.simplify(As + Bs) == 0
    assert sp.simplify(sp.diff(phi1, t).subs({A: As, z: 0}) / E - om ** 2 * a / k) == 0  # step 6
    dyn = (r1 * sp.diff(phi1, t) + r1 * g * zeta - r2 * sp.diff(phi2, t) - r2 * g * zeta).subs({A: As, B: Bs, z: 0}) / E
    sols = sp.solve(sp.expand(sp.simplify(dyn)), om ** 2)  # steps 7–9
    assert len(sols) == 1 and sp.simplify(sols[0] - g * k * (r2 - r1) / (r2 + r1)) == 0
    assert sp.limit(sp.sqrt(g * k * (r2 - r1) / (r2 + r1)), r1, 0) == sp.sqrt(g * k)  # step 10


def test_interface_fields_V1_residuals_vortex_sheet_and_parity():  # V1 (7.90)–(7.95), N107, Fig. 7.24 (N190)
    a, k, r1, r2 = 0.5, 0.2, 1000.0, 1020.0
    x = np.linspace(0, 30, 13)
    r = ch07.interface_residuals(x, 1.7, a, k, r1, r2, G)
    om = float(W.interface_omega(k, r1, r2, G))
    for key in ("kinematic1", "kinematic2"):
        assert np.max(np.abs(r[key])) < 1e-15 * 10 + 1e-14 * a * om
    assert np.max(np.abs(r["dynamic"])) < 1e-12 * r2 * G * a
    assert np.max(np.abs(r["laplace1"])) < 1e-8 * a * om * k and np.max(np.abs(r["laplace2"])) < 1e-8 * a * om * k
    assert np.max(np.abs(r["far1"])) < 1e-15 and np.max(np.abs(r["far2"])) < 1e-15  # (7.91)–(7.92): e^{−40}
    f = ch07.interface_fields(x, 0.0, 1.7, a, k, r1, r2, G)
    assert f["A"] == pytest.approx(1j * om * a / k) and f["B"] == pytest.approx(-f["A"])
    assert np.max(np.abs(f["u1"] + f["u2"])) < 1e-14  # opposite tangential velocities at z = 0: a vortex sheet
    th = k * x - om * 1.7
    assert maxrel(f["gamma_sheet"], 2 * om * a * np.cos(th)) < 1e-14
    assert maxrel(f["gamma_sheet"], ch05.vortex_sheet_strength(f["u1"], f["u2"])) < 1e-13  # u_below − u_above (ccw)
    assert maxrel(f["zeta_interface"], a * np.cos(th)) < 1e-15
    assert maxrel(ch07.interface_fields(x, 0.5, 1.7, a, k, r1, r2, G)["u"], -om * a * np.exp(-0.5 * k) * np.cos(th)) < 1e-13


def test_interface_energy_V1_quarter_from_direct_integration():  # V1 N105 (Ex. 7.18), (7.96): E_k = E_p = ¼Δρga²
    a, k, r1, r2 = 0.1, 1.0, 1000.0, 1020.0
    q = ch07.interface_energy(a, k, r1, r2, G, method="quad")
    quarter = 0.25 * (r2 - r1) * G * a ** 2
    assert q["Ek"] == pytest.approx(quarter, rel=1e-8) and q["Ep"] == pytest.approx(quarter, rel=1e-8)
    assert ch07.interface_energy(a, k, r1, r2, G)["E"] == pytest.approx(float(W.wave_energy_density(a, g=G, drho=r2 - r1)))
    lam = TWO_PI / k  # the printed middle form (g(ρ₂ − ρ₁)/2λ)∫₀^{λ/2}ζ²dx is ⅛, not ¼ (implementer's note 4)
    mid = G * (r2 - r1) / (2 * lam) * quad(lambda xx: (a * np.cos(k * xx)) ** 2, 0, lam / 2)[0]
    assert mid == pytest.approx(quarter / 2, rel=1e-12) and q["Ep"] != pytest.approx(mid, rel=0.1)


def test_interface_omega_V7_limits_and_rayleigh_taylor():  # V7 (7.95): ρ₁ → 0, ρ₁ → ρ₂; ρ₁ > ρ₂ → NaN + warning
    k = np.logspace(-3, 1, 20)
    assert maxrel(W.interface_omega(k, 0.0, 1000.0, G), np.sqrt(G * k)) < 1e-15
    assert np.max(W.interface_omega(k, 1000.0, 1000.0, G)) == 0.0
    assert maxrel(W.interface_omega(k, 998.0, 1002.0, G), np.sqrt(W.eps2_density(998.0, 1002.0)) * np.sqrt(G * k)) < 1e-15
    with pytest.warns(RuntimeWarning):
        assert np.isnan(float(W.interface_omega(1.0, 1020.0, 1000.0, G)))


@needs_ref
def test_interface_V1_form_published_two_fluid_dispersion():  # V1 form: Wikipedia's interface relation with σ = 0
    k = np.logspace(-3, 2, 40)
    r1, r2 = 1.2, 1000.0
    wiki = np.sqrt(np.abs(k) * ((r2 - r1) / (r2 + r1) * G))
    assert maxrel(W.interface_omega(k, r1, r2, G), wiki) < 1e-15
    assert "rho_p" in ref_json()["capillary_wave"]["dispersion"]


# =====================================================================================================================
# C14 — two-layer modes (7.97)–(7.119): N109–N131, N193, N194, D29–D31
# =====================================================================================================================
def _two_layer_symbols():
    return sp.symbols("a k H g omega rho1 rho2", positive=True)


def test_two_layer_constants_V2_derivation():  # V2 D29 (★★★): steps 3–12 line by line; the printed (7.105) fails
    a, k, H, g, om, r1, r2 = _two_layer_symbols()
    x, z, t = sp.symbols("x z t", real=True)
    A, B, C, b = sp.symbols("A B C b")
    I = sp.I
    E = sp.exp(I * (k * x - om * t))
    phi1 = (A * sp.exp(k * z) + B * sp.exp(-k * z)) * E  # (7.104)
    phi2 = C * sp.exp(k * z) * E  # (7.105) with e^{i(kx − ωt)}
    eta, zeta = a * E, b * E  # (7.102), (7.103)
    s3 = sp.simplify((sp.diff(phi1, z) - sp.diff(eta, t)).subs(z, 0) / E)  # step 3: k(A − B) = −iωa
    assert sp.simplify(s3 - (k * (A - B) + I * om * a)) == 0
    s4 = sp.simplify((sp.diff(phi1, t) + g * eta).subs(z, 0) / E)  # step 4: A + B = −iga/ω
    assert sp.simplify(sp.solve(s4, A)[0] + B + I * g * a / om) == 0
    AB = sp.solve([s3, s4], [A, B], dict=True)[0]
    A_b, B_b = -I * a / 2 * (om / k + g / om), I * a / 2 * (om / k - g / om)
    assert sp.simplify(AB[A] - A_b) == 0 and sp.simplify(AB[B] - B_b) == 0  # steps 5–6 (7.106), (7.107)
    s7 = sp.simplify((sp.diff(phi1, z) - sp.diff(phi2, z)).subs(z, -H) / E)  # step 7
    assert sp.simplify(s7 - (k * (A * sp.exp(-k * H) - B * sp.exp(k * H)) - k * C * sp.exp(-k * H))) == 0
    Cs = sp.solve(s7, C)[0]
    assert sp.simplify(Cs - (A - B * sp.exp(2 * k * H))) == 0  # step 8
    C_b = -I * a / 2 * (om / k + g / om) - I * a / 2 * (om / k - g / om) * sp.exp(2 * k * H)
    assert sp.simplify(Cs.subs(AB) - C_b) == 0  # (7.108)
    s9 = sp.simplify((sp.diff(phi1, z) - sp.diff(zeta, t)).subs(z, -H) / E)  # step 9
    bs = sp.solve(s9, b)[0]
    assert sp.simplify(bs - I * k / om * (A * sp.exp(-k * H) - B * sp.exp(k * H))) == 0  # step 10
    assert sp.simplify(I * k / om * A_b * sp.exp(-k * H) - a / 2 * (1 + g * k / om ** 2) * sp.exp(-k * H)) == 0  # step 11
    assert sp.simplify(-I * k / om * B_b * sp.exp(k * H) - a / 2 * (1 - g * k / om ** 2) * sp.exp(k * H)) == 0
    b_b = a / 2 * (1 + g * k / om ** 2) * sp.exp(-k * H) + a / 2 * (1 - g * k / om ** 2) * sp.exp(k * H)
    assert sp.simplify(bs.subs(AB) - b_b) == 0  # step 12 (7.109)
    assert sp.simplify(b_b.subs(H, 0) - a) == 0  # check: H → 0, the interface is the surface
    phi2_printed = C * sp.exp(k * z) * sp.exp(I * (k * z - om * t))  # the printed (7.105)
    r = (sp.diff(phi1, z) - sp.diff(phi2_printed, z)).subs(z, -H).subs(AB).subs(C, Cs.subs(AB))
    assert sp.simplify(sp.diff(r, x)) != 0  # cannot hold for every x
    tl = ch07.two_layer_sympy()  # the coded construction agrees
    assert all(sp.simplify(v) == 0 for v in tl["residuals"].values()) and all(sp.simplify(v) == 0 for v in
                                                                                 tl["bc_residuals"])
    assert sp.simplify(tl["printed_7105_residual"]) != 0


def test_two_layer_dispersion_V2_derivation():  # V2 D30 (★★★): (7.101) with (7.106)–(7.109) → (7.110), step by step
    a, k, H, g, om, r1, r2 = _two_layer_symbols()
    s = sp.Symbol("s", positive=True)
    I = sp.I
    eH, emH = sp.exp(k * H), sp.exp(-k * H)
    A = -I * a / 2 * (om / k + g / om)
    B = I * a / 2 * (om / k - g / om)
    CemH = A * emH - B * eH  # step 4 (D29 step 7)
    b = a / 2 * (1 + g * k / om ** 2) * emH + a / 2 * (1 - g * k / om ** 2) * eH
    lhs3 = -I * om * (r1 * (A * emH + B * eH) - r2 * CemH)  # step 3
    rhs3 = (r2 - r1) * g * b
    lhs5 = -I * om * ((r1 - r2) * A * emH + (r1 + r2) * B * eH)  # step 5
    assert sp.simplify(sp.expand(lhs3 - lhs5)) == 0
    assert sp.simplify(-I * om * A + a / 2 * (om ** 2 / k + g)) == 0  # step 6
    assert sp.simplify(-I * om * B - a / 2 * (om ** 2 / k - g)) == 0
    sub = {om: sp.sqrt(s * g * k)}  # step 7: s = ω²/gk
    left8 = (r2 - r1) * (s + 1) * emH + (r1 + r2) * (s - 1) * eH  # step 8
    right8 = (r2 - r1) * ((1 + 1 / s) * emH + (1 - 1 / s) * eH)
    assert sp.simplify(sp.expand(lhs5.subs(sub) * 2 / (a * g) - left8)) == 0
    assert sp.simplify(sp.expand(rhs3.subs(sub) * 2 / (a * g) - right8)) == 0
    step9 = (s - 1) / s * ((r2 - r1) * (s + 1) * emH + ((r1 + r2) * s - (r2 - r1)) * eH)
    assert sp.simplify(sp.expand(left8 - right8 - step9)) == 0  # step 9: the common factor (s − 1)/s
    brace = (r2 - r1) * (s + 1) * emH + ((r1 + r2) * s - (r2 - r1)) * eH
    step11 = 2 * s * (r1 * sp.sinh(k * H) + r2 * sp.cosh(k * H)) - 2 * (r2 - r1) * sp.sinh(k * H)
    assert sym_zero(brace - step11)  # step 11
    e7110 = (s - 1) * (s * (r1 * sp.sinh(k * H) + r2 * sp.cosh(k * H)) - (r2 - r1) * sp.sinh(k * H))  # (7.110) in s
    assert sym_zero(step9 - 2 / s * e7110)  # step 12: divide by the nonzero 2/s
    residual = (lhs3 - rhs3).subs(sub)  # the full residual of (7.101) = (ag²k/ω²) × (7.110)
    assert sym_zero(residual - a * g / s * e7110)
    tl = ch07.two_layer_sympy()
    assert sym_zero(tl["pressure_residual"] - a * g ** 2 * k / om ** 2 * tl["dispersion_7110"])
    assert sym_zero(tl["common_factor"] - a * g ** 2 * k / om ** 2)


def test_two_layer_modes_V2_derivation():  # V2 D31 (★★): roots (7.111), (7.113); b (7.112); η/ζ (7.114); long waves
    a, k, H, g, om, r1, r2 = _two_layer_symbols()
    kh = sp.Symbol("kh", positive=True)
    b_of = lambda w2: a / 2 * (1 + g * k / w2) * sp.exp(-k * H) + a / 2 * (1 - g * k / w2) * sp.exp(k * H)  # noqa: E731
    assert sp.simplify(b_of(g * k) - a * sp.exp(-k * H)) == 0  # steps 1–2 (7.112)
    w2bc = g * k * (r2 - r1) * sp.sinh(k * H) / (r2 * sp.cosh(k * H) + r1 * sp.sinh(k * H))  # step 4 (7.113)
    s = w2bc / (g * k)
    D = (r2 - r1) * sp.sinh(k * H)
    assert sym_zero(1 + 1 / s - r2 * sp.exp(k * H) / D)  # step 5
    assert sym_zero(1 - 1 / s + (r2 * sp.exp(-k * H) + 2 * r1 * sp.sinh(k * H)) / D)
    assert sym_zero(b_of(w2bc) + a * r1 * sp.exp(k * H) / (r2 - r1))  # step 6
    assert sym_zero(a / b_of(w2bc) + (r2 - r1) / r1 * sp.exp(-k * H))  # step 7 (7.114)
    assert sp.limit(w2bc.subs(H, kh / k).rewrite(sp.exp), kh, sp.oo) == g * k * (r2 - r1) / (r2 + r1)  # step 8 → (7.95)
    ser = sp.series(w2bc.subs(H, kh / k), kh, 0, 2).removeO()  # step 9 (7.115)
    assert sp.simplify(ser - g * k * (r2 - r1) / r2 * kh) == 0
    gp = g * (r2 - r1) / r2
    assert sp.simplify(g * k * (r2 - r1) / r2 * (k * H) / k ** 2 - gp * H) == 0  # step 10: c² = g′H (7.116)–(7.117)
    AplusB = -sp.I * g * a / om  # step 11 (7.119): p′ = iρ₁ω(A + B) = ρ₁ga
    assert sp.simplify(sp.I * r1 * om * AplusB - r1 * g * a) == 0


def test_two_layer_residuals_V1_numeric_and_printed_variant():  # V1 (7.97)–(7.101); T3 slip fails
    x = np.linspace(0, 60, 9)
    for mode in ("barotropic", "baroclinic"):
        r = ch07.two_layer_residuals(x, 3.3, 0.1, 40.0, 1000.0, 1003.0, G, 1.0, mode)
        md = ch07.two_layer_modes(0.1, 40.0, 1000.0, 1003.0, G, 1.0, mode)
        scale = max(abs(md["A"]), abs(md["B"])) * 0.1
        for key in ("laplace2", "kin_surface", "kin_interface_1", "kin_interface_2"):
            assert np.max(np.abs(r[key])) < 1e-12 * max(scale, 1.0), (mode, key)
        assert np.max(np.abs(r["dyn_surface"])) < 1e-12 * G and np.max(np.abs(r["dyn_interface"])) < 1e-9 * 1003 * G
        assert np.max(np.abs(r["decay2"])) < 1e-15 * max(1.0, abs(md["C"]))
    bad = ch07.two_layer_residuals(x, 3.3, 0.1, 40.0, 1000.0, 1003.0, G, 1.0, "baroclinic", printed_7_105=True)
    assert np.max(np.abs(bad["laplace2"])) > 1e-6 and np.max(np.abs(bad["kin_interface_2"])) > 1e-6


def test_two_layer_modes_V1_amplitude_ratios_and_pressure():  # V1 (7.112), (7.114), (7.118), (7.119); Figs. 7.27–7.28
    k, H, r1, r2 = 0.05, 40.0, 1000.0, 1004.0
    bt = ch07.two_layer_modes(k, H, r1, r2, G, 1.0, "barotropic")
    assert bt["omega"] == pytest.approx(np.sqrt(G * k), rel=1e-15)
    assert bt["b"] == pytest.approx(bt["b_book"], rel=1e-12) and np.imag(bt["b"]) == 0  # in phase, e^{−kH}
    bc = ch07.two_layer_modes(k, H, r1, r2, G, 1.0, "baroclinic")
    assert bc["eta_over_zeta"] == pytest.approx(bc["eta_over_zeta_book"], rel=1e-9) and bc["eta_over_zeta"] < 0
    assert bc["p_prime_check"]["top_over_hydrostatic"] == pytest.approx(1.0, rel=1e-12)  # p′ = ρ₁gη at z = 0
    lw = ch07.two_layer_modes(1e-5, H, r1, r2, G, 1.0, "baroclinic", long_wave=True)
    assert lw["c"] == pytest.approx(np.sqrt(G * (r2 - r1) / r2 * H), rel=1e-15)
    assert lw["omega"] / 1e-5 == pytest.approx(lw["c"], rel=1e-12)
    assert lw["omega_exact"] / 1e-5 == pytest.approx(lw["c"], rel=1e-3)  # (7.115) is the kH → 0 limit
    assert lw["eta_over_zeta_exact"] == pytest.approx(-(r2 - r1) / r1, rel=1e-3)  # (7.118)
    assert lw["p_prime_check"]["interface_over_top"] == pytest.approx(1.0, abs=1e-3)  # hydrostatic upper layer
    with pytest.raises(ValueError):
        ch07.two_layer_modes(k, H, r1, r2, G, 1.0, "barotropic", long_wave=True)
    with pytest.raises(ValueError):
        ch07.two_layer_modes(k, H, r1, r2, G, 1.0, "sideways")


def test_two_layer_V7_limits_and_reduced_gravity():  # V7 (7.110)–(7.117); R7 g′ conventions; Ex. 7.20 rigid lid
    r1, r2 = 1000.0, 1005.0
    k = 0.3
    wbt, wbc = W.two_layer_free_surface_omega(k, 200.0, r1, r2, G)  # kH = 60 → (7.95)
    assert float(wbc) == pytest.approx(float(W.interface_omega(k, r1, r2, G)), rel=1e-12)
    assert float(wbt) == pytest.approx(np.sqrt(G * k), rel=1e-15)
    # both roots satisfy (7.110)
    for wv in (wbt, wbc):
        s_ = float(wv) ** 2 / (G * k)
        th = k * 200.0
        assert (s_ - 1) * (s_ * (r1 * np.sinh(th) + r2 * np.cosh(th)) - (r2 - r1) * np.sinh(th)) == pytest.approx(
            0.0, abs=1e-9 * r2 * np.cosh(th))
    kk = np.array([1e-7, 1e-6, 1e-5])  # kH ≤ 3e-4: (7.115) holds to O(ρ₁kH/2ρ₂)
    c_bc = W.two_layer_free_surface_omega(kk, 30.0, r1, r2, G)[1] / kk
    assert maxrel(c_bc, np.full(3, float(W.two_layer_long_wave_speed(30.0, r1, r2, G)))) < 1e-3
    assert float(W.two_layer_free_surface_omega(0.1, 10.0, r1, r1, G)[1]) == 0.0  # ρ₁ = ρ₂: no internal mode
    gl = W.reduced_gravity_book(r1, r2, G)
    assert gl == pytest.approx(G * (r2 - r1) / r2) and W.reduced_gravity_book(r1, r2, G, "upper") == pytest.approx(
        float(SIM.reduced_gravity(r1, r2, G)))  # ch04's ρ₁ form
    assert gl / W.reduced_gravity_book(r1, r2, G, "upper") == pytest.approx(r1 / r2, rel=1e-15)
    assert float(SIM.reduced_gravity(r1, r2, G, ref="lower")) == pytest.approx(gl)
    with pytest.raises(ValueError):
        W.reduced_gravity_book(r1, r2, G, ref="middle")
    big = ch07.two_layer_rigid_lid_omega(0.2, 1e3, 1e3, r1, r2, G)  # h₁, h₂ → ∞ → (7.95)
    assert float(big) == pytest.approx(float(W.interface_omega(0.2, r1, r2, G)), rel=1e-12)
    h1, h2 = 20.0, 60.0
    c_rl = float(ch07.two_layer_rigid_lid_omega(1e-6, h1, h2, r1, r2, G)) / 1e-6
    assert c_rl ** 2 == pytest.approx(G * (r2 - r1) * h1 * h2 / (r1 * h2 + r2 * h1), rel=1e-9)


def test_two_layer_state_V1_explainer_numbers():  # V1 Part C 2.35 (E8)
    s = ch07.two_layer_state(0.05, 40.0, 1000.0, 1002.0, G)
    assert s["omega"] == pytest.approx(s["omega_bc"]) and s["c_bt"] == pytest.approx(s["omega_bt"] / 0.05)
    assert s["g_prime_lower"] == pytest.approx(G * 2 / 1002) and s["g_prime_upper"] == pytest.approx(G * 2 / 1000)
    assert s["c_long"] == pytest.approx(np.sqrt(G * 2 / 1002 * 40.0)) and s["kH"] == pytest.approx(2.0)
    assert s["eps2_density"] == pytest.approx(2 / 2002) and s["T"] == pytest.approx(TWO_PI / s["omega"])
    assert ch07.two_layer_state(0.05, 40.0, mode="barotropic")["eta_over_zeta"] > 0


# =====================================================================================================================
# C15 — internal-wave dispersion (7.120)–(7.139): R12–R20, N132–N147, N196, D32–D34
# =====================================================================================================================
def test_boussinesq_linear_V2_derivation():  # V2 D32 (★★): (7.124) in (7.120)–(7.125) → (7.126), (7.128)–(7.131)
    x, y, z, t = sp.symbols("x y z t", real=True)
    g, rho0, eps = sp.symbols("g rho0 epsilon", positive=True)
    u, v, w, rp = [sp.Function(n)(x, y, z, t) for n in ("u", "v", "w", "rp")]
    rbar = sp.Function("rhobar")(z)
    full = sp.diff(rbar + eps * rp, t) + sum(eps * c * sp.diff(rbar + eps * rp, s) for c, s in ((u, x), (v, y), (w, z)))
    assert sp.diff(rbar, t) == 0 and sp.diff(rbar, x) == 0 and sp.diff(rbar, y) == 0  # step 3
    lin = sp.expand(full).coeff(eps, 1)  # step 4: keep O(ε)
    assert sp.simplify(lin - (sp.diff(rp, t) + w * sp.diff(rbar, z))) == 0  # (7.126): w dρ̄/dz kept
    quadratic = sp.expand(full).coeff(eps, 2)
    assert sp.simplify(quadratic - (u * sp.diff(rp, x) + v * sp.diff(rp, y) + w * sp.diff(rp, z))) == 0  # dropped
    N2 = -g / rho0 * sp.diff(rbar, z)  # (7.127), step 5
    assert sp.simplify(lin - (sp.diff(rp, t) - N2 * rho0 / g * w)) == 0  # (7.131)
    res = ch07.boussinesq_linear_sympy()
    assert all(sp.simplify(res[k]) == 0 for k in ("r7126", "r7128", "r7129", "r7130", "r7131"))


def test_w_equation_V2_derivation():  # V2 D33 (★★★): (4.10), (7.128)–(7.131) → (7.132), (7.133) → (7.134), N(z) allowed
    x, y, z, t = sp.symbols("x y z t", real=True)
    rho0, g = sp.symbols("rho0 g", positive=True)
    N = sp.Function("N")(z)
    u, v, w, p, r = [sp.Function(n)(x, y, z, t) for n in ("u", "v", "w", "p", "r")]
    lapH = lambda f: sp.diff(f, x, 2) + sp.diff(f, y, 2)  # noqa: E731
    R1, R2 = sp.diff(u, t) + sp.diff(p, x) / rho0, sp.diff(v, t) + sp.diff(p, y) / rho0  # (7.128), (7.129)
    R3 = sp.diff(w, t) + sp.diff(p, z) / rho0 + r * g / rho0  # (7.130)
    R4 = sp.diff(r, t) - N ** 2 * rho0 / g * w  # (7.131)
    R5 = sp.diff(u, x) + sp.diff(v, y) + sp.diff(w, z)  # (4.10)
    s1 = sp.diff(R5, t)  # step 1: ∂/∂t of continuity
    s2 = sp.diff(sp.diff(u, t), x) + sp.diff(sp.diff(v, t), y) + sp.diff(w, z, t)  # step 2: swap the derivatives
    assert sp.simplify(s1 - s2) == 0
    s3 = s2 - sp.diff(R1, x) - sp.diff(R2, y)  # step 3: insert ∂u/∂t = −p′_x/ρ₀ − … i.e. subtract the residuals
    assert sp.simplify(sp.expand(s3 - (-lapH(p) / rho0 + sp.diff(w, z, t)))) == 0
    eq132 = lapH(p) / rho0 - sp.diff(w, z, t)  # step 4: (7.132) = −(∂_tR₅ − ∂_xR₁ − ∂_yR₂)
    assert sp.simplify(sp.expand(eq132 + (sp.diff(R5, t) - sp.diff(R1, x) - sp.diff(R2, y)))) == 0
    s5 = sp.diff(R3, t)  # step 5
    s6 = g / rho0 * sp.diff(r, t) - N ** 2 * w  # step 6: (g/ρ₀)∂ρ′/∂t = N²w when (7.131) holds
    assert sp.simplify(sp.expand(s6 - g / rho0 * R4)) == 0
    eq133 = sp.diff(p, t, z) / rho0 + sp.diff(w, t, 2) + N ** 2 * w  # step 7: (7.133) = ∂_tR₃ − (g/ρ₀)R₄
    assert sp.simplify(sp.expand(eq133 - (s5 - g / rho0 * R4))) == 0
    s8 = lapH(eq133)  # step 8 (∇_H² commutes with ∂_t, ∂_z and N(z))
    assert sp.simplify(sp.expand(s8 - (sp.diff(lapH(p), t, z) / rho0 + lapH(sp.diff(w, t, 2)) + N ** 2 * lapH(w)))) == 0
    s9 = sp.diff(eq132, t, z)  # step 9
    w134 = sp.diff(w, t, 2, z, 2) + lapH(sp.diff(w, t, 2)) + N ** 2 * lapH(w)  # steps 10–12 (7.134)
    assert sp.simplify(sp.expand(s9 - s8 + w134)) == 0  # p′ is eliminated exactly, for any N(z)
    res = ch07.boussinesq_linear_sympy()
    assert all(sp.simplify(res[k]) == 0 for k in ("r7132", "r7133", "r7134", "r7137"))


def test_internal_dispersion_V2_derivation():  # V2 D34 (★★): (7.136) in (7.134) → (7.137) → ω = N cos θ (7.139)
    x, y, z, t = sp.symbols("x y z t", real=True)
    k, l, m, om, N, w0 = sp.symbols("k l m omega N w0", positive=True)
    kr = sp.Symbol("k_r", real=True, nonzero=True)
    wp = w0 * sp.exp(sp.I * (k * x + l * y + m * z - om * t))  # step 2
    lap = sp.diff(wp, x, 2) + sp.diff(wp, y, 2) + sp.diff(wp, z, 2)
    assert sp.simplify(lap / wp + (k ** 2 + l ** 2 + m ** 2)) == 0  # step 3
    eq = sp.simplify((sp.diff(lap, t, 2) + N ** 2 * (sp.diff(wp, x, 2) + sp.diff(wp, y, 2))) / wp)  # step 4
    assert sp.simplify(eq - (om ** 2 * (k ** 2 + l ** 2 + m ** 2) - N ** 2 * (k ** 2 + l ** 2))) == 0
    sol = sp.solve(eq, om ** 2)  # step 5 (7.137)
    assert len(sol) == 1 and sp.simplify(sol[0] - N ** 2 * (k ** 2 + l ** 2) / (k ** 2 + l ** 2 + m ** 2)) == 0
    Km = sp.sqrt(kr ** 2 + m ** 2)
    om138 = N * sp.sqrt(kr ** 2 / Km ** 2)  # step 6, l = 0: N|k|/K for either sign of k
    assert sp.simplify(om138 - N * sp.Abs(kr) / Km) == 0
    th = sp.Symbol("theta", positive=True)  # steps 7–8: cos θ = |k|/K
    assert sp.simplify((N * sp.cos(th)).subs(th, sp.acos(sp.Abs(kr) / Km)) - N * sp.Abs(kr) / Km) == 0


def test_internal_wave_omega_V7_limits_and_direction_only():  # V7 (7.137)–(7.139), N145: θ = 0 ⇒ ω = N (ch01 parcel)
    N = 0.012
    assert float(W.internal_wave_omega(1.0, 0.0, N)) == pytest.approx(N, rel=1e-15)
    t = np.linspace(0, 3 * TWO_PI / N, 7)
    assert maxrel(ST.parcel_displacement(t, 1.0, N ** 2), np.cos(float(W.internal_wave_omega(0.3, 0.0, N)) * t)) < 1e-12
    k, m = 0.02, 0.05
    w1 = float(W.internal_wave_omega(k, m, N))
    for s in (0.1, 3.0, 1e4):  # frequency depends on direction only
        assert float(W.internal_wave_omega(s * k, s * m, N)) == pytest.approx(w1, rel=1e-14)
    ang = RNG.uniform(0, TWO_PI, 200)
    ws = W.internal_wave_omega(np.cos(ang), np.sin(ang), N)
    assert np.all(ws <= N * (1 + 1e-15)) and maxrel(ws, N * np.abs(np.cos(ang))) < 1e-14  # ω = N|cos θ| ≤ N
    assert float(W.internal_wave_omega(k, m, N, l=0.03)) == pytest.approx(
        N * np.sqrt((k ** 2 + 0.03 ** 2) / (k ** 2 + 0.03 ** 2 + m ** 2)), rel=1e-15)
    assert float(W.beam_angle(0.6, 1.0)) == pytest.approx(np.arccos(0.6), rel=1e-15)
    assert np.isnan(float(W.beam_angle(1.2, 1.0))) and float(W.beam_angle(1.0, 1.0)) == 0.0


def test_internal_wave_fields_V3_residuals_second_order():  # V3 (4.10), (7.128), (7.130), (7.131): FD residuals → 0 at O(h²)
    X, Z = np.linspace(0, 3, 5), np.linspace(-1, 1, 5)
    hs = [1e-2, 5e-3, 2.5e-3]
    out = {k: [] for k in ("continuity", "momentum_x", "momentum_z", "density")}
    for h in hs:
        f = ch07.internal_wave_fields(X, Z, 0.4, 1.0, 1.5, 1.0, 0.1, residuals=True, h=h, ht=h)
        for key in out:
            out[key].append(np.max(np.abs(f["residuals"][key])))
    for key, e in out.items():
        assert abs(observed_order(hs, e) - 2.0) < ORDER_TOL, key


def test_internal_wave_fields_V1_polarization_and_transversality():  # V1 (7.141), (7.149), (7.153): N144, D35
    k, m, N, w0, r0 = -0.7, 1.3, 0.9, 0.05, 1025.0
    X, Z, T = RNG.uniform(-5, 5, 30), RNG.uniform(-5, 5, 30), 2.3
    f = ch07.internal_wave_fields(X, Z, T, k, m, N, w0, r0, G)
    om = float(W.internal_wave_omega(k, m, N))
    th = k * X + m * Z - om * T
    assert maxrel(f["w"], w0 * np.cos(th)) < 1e-14 and maxrel(f["u"], -m / k * w0 * np.cos(th)) < 1e-14
    assert np.max(np.abs(f["K_dot_u"])) < 1e-15 * 10  # (7.141) K·u = 0
    assert maxrel(f["p_prime"], -om * m * r0 / k ** 2 * w0 * np.cos(th)) < 1e-13
    assert maxrel(f["rho_prime"], N ** 2 * r0 / G * f["zeta_particle"]) < 1e-13  # (7.149) ρ′ = N²ρ₀ζ/g
    assert maxrel(f["zeta_particle"], -w0 / om * np.sin(th)) < 1e-14  # w = ∂ζ/∂t
    ph = np.linspace(0, TWO_PI, 64, endpoint=False)
    prod = np.mean(W.real_field(f["amplitudes"]["rho"], ph) * W.real_field(f["amplitudes"]["w"], ph))
    assert abs(prod) < 1e-14 * abs(f["amplitudes"]["rho"]) * w0  # ρ′ is 90° out of phase with w: no mean product


def test_w_equation_residual_V1_plane_wave_right_and_wrong_frequency():  # V1/V3 (7.134) checked with a callable w
    k, m, N = 1.0, 1.5, 1.0
    om = float(W.internal_wave_omega(k, m, N))
    X, Z = np.linspace(0, 3, 5), np.linspace(-1, 1, 5)
    hs, errs = [0.08, 0.04, 0.02], []
    for h in hs:
        errs.append(np.max(np.abs(ch07.w_equation_residual(lambda a, b, c: 0.1 * np.cos(k * a + m * b - om * c), X, Z,
                                                           0.3, N, h=h, ht=h))))
    assert abs(observed_order(hs, errs) - 2.0) < ORDER_TOL
    bad = ch07.w_equation_residual(lambda a, b, c: 0.1 * np.cos(k * a + m * b - 1.1 * om * c), X, Z, 0.3, N)
    assert np.max(np.abs(bad)) > 100 * errs[-1]


def test_layered_flow_V1_horizontal_nondivergent_layers():  # V1 (7.142), N146
    u = lambda x, y, z: np.sin(x) * np.cos(y) * (1 + z ** 2)  # noqa: E731
    v = lambda x, y, z: -np.cos(x) * np.sin(y) * (1 + z ** 2)  # noqa: E731
    r = ch07.layered_flow_check(u, v, RNG.uniform(-2, 2, 10), RNG.uniform(-2, 2, 10), z=0.7)
    assert np.max(np.abs(r["divergence"])) < 1e-9 and np.all(r["density"] == 0) and np.all(r["momentum_z"] == 0)
    rb = ch07.layered_flow_check(lambda x, y, z: x, lambda x, y, z: 0 * y, 0.3, 0.1)
    assert float(rb["divergence"]) == pytest.approx(1.0)  # a divergent layer is flagged


# =====================================================================================================================
# C16 — c ⟂ c_g, beams, energy and F = c_g E (7.140)–(7.159): N148–N166, N195, N197–N199, D35–D37
# =====================================================================================================================
def test_internal_group_velocity_V2_derivation():  # V2 D36 (★★): ∇_K ω → (7.145); c (7.144); c·c_g = 0 (7.146)
    k, m, N = sp.symbols("k m N", positive=True)
    K = sp.sqrt(k ** 2 + m ** 2)
    om = N * k / K  # step 2
    dk, dm = sp.diff(om, k), sp.diff(om, m)
    assert sp.simplify(dk - N * m ** 2 / K ** 3) == 0 and sp.simplify(dm + N * k * m / K ** 3) == 0  # steps 3–4
    cg = sp.Matrix([dk, dm])
    assert sp.simplify(cg - N * m / K ** 3 * sp.Matrix([m, -k])) == sp.zeros(2, 1)  # step 5 (7.145)
    c = om / K ** 2 * sp.Matrix([k, m])  # step 6 (7.144)
    assert sp.simplify(c.dot(cg)) == 0  # step 8 (7.146)
    assert sp.simplify(c[1] + cg[1]) == 0  # step 9: c_z = −c_g,z
    kr, q = sp.Symbol("k_r", real=True), sp.Symbol("q", positive=True)  # step 7: sign-safe form ω = N|k|/K
    Kr = sp.sqrt(kr ** 2 + m ** 2)
    omr = N * sp.sqrt(kr ** 2) / Kr
    cgr = sp.Matrix([sp.diff(omr, kr), sp.diff(omr, m)])
    cr = omr / Kr ** 2 * sp.Matrix([kr, m])
    for sgn in (1, -1):  # k = ±q: ∂ω/∂k = sgn(k)Nm²/K³, ∂ω/∂m = −N|k|m/K³, and c·c_g = 0 for either sign
        Kq = sp.sqrt(q ** 2 + m ** 2)
        assert sp.simplify(cgr[0].subs(kr, sgn * q) - sgn * N * m ** 2 / Kq ** 3) == 0
        assert sp.simplify(cgr[1].subs(kr, sgn * q) + N * q * m / Kq ** 3) == 0
        assert sp.simplify(cr.dot(cgr).subs(kr, sgn * q)) == 0


def test_polarization_and_flux_V2_derivation():  # V2 D37 (★★): (7.153), (7.158), F = c_g E (7.159)
    k, m, N, w0, rho0, g, om = sp.symbols("k m N w0 rho0 g omega", positive=True)
    P, Rh, U = sp.symbols("P R U")
    p_hat = sp.solve(sp.Eq(-k ** 2 / rho0 * P, (sp.I * m) * (-sp.I * om) * w0), P)[0]  # step 2 from (7.132)
    assert sp.simplify(p_hat + om * m * rho0 / k ** 2 * w0) == 0
    r_hat = sp.solve(sp.Eq(-sp.I * om * Rh, N ** 2 * rho0 / g * w0), Rh)[0]  # step 3 from (7.131)
    assert sp.simplify(r_hat - sp.I * N ** 2 * rho0 / (om * g) * w0) == 0
    u_hat = sp.solve(sp.Eq(-sp.I * om * U, -sp.I * k / rho0 * p_hat), U)[0]  # step 4 from (7.128)
    assert sp.simplify(u_hat + m / k * w0) == 0 and sp.simplify(k * u_hat + m * w0) == 0  # K·u = 0
    Fx = sp.re(p_hat * sp.conjugate(u_hat)) / 2  # steps 7–8
    Fz = sp.re(p_hat * sp.conjugate(w0)) / 2  # step 9
    assert sp.simplify(Fx - rho0 * om * m ** 2 * w0 ** 2 / (2 * k ** 3)) == 0
    assert sp.simplify(Fz + rho0 * om * m * w0 ** 2 / (2 * k ** 2)) == 0
    K = sp.sqrt(k ** 2 + m ** 2)
    E = rho0 / 2 * (m ** 2 / k ** 2 + 1) * w0 ** 2  # (7.157)
    cgE = N * m / K ** 3 * sp.Matrix([m, -k]) * E  # step 11
    assert sp.simplify(cgE - rho0 * N * m * w0 ** 2 / (2 * K * k ** 2) * sp.Matrix([m, -k])) == sp.zeros(2, 1)
    cgE_w = cgE.subs(N, om * K / k)  # step 12 with ω = kN/K
    assert sp.simplify(cgE_w - sp.Matrix([Fx, Fz])) == sp.zeros(2, 1)  # (7.159)
    Ek = rho0 / 4 * (m ** 2 / k ** 2 + 1) * w0 ** 2  # (7.154)
    Ep = (N ** 2 * rho0 / (4 * om ** 2) * w0 ** 2).subs(om, N * k / K)  # (7.155)
    assert sp.simplify(Ek - Ep) == 0  # (7.156)


def test_internal_wave_velocities_V1_gradient_parity_both_signs():  # V1 (7.143)–(7.146), Fig. 7.29 geometry (N195)
    N = 0.8
    for k, m in ((0.6, 0.9), (-0.6, 0.9), (0.6, -0.9), (-1.4, -0.3)):
        v = W.internal_wave_velocities(k, m, N)
        num = W.group_velocity_vector(lambda K: W.internal_wave_omega(K[0], K[1], N), [k, m])
        assert maxrel(v["cg"], num) < 1e-12
        assert abs(v["dot"]) < 1e-15 and v["c"][1] == pytest.approx(-v["cg"][1], rel=1e-14)
        assert np.sign(v["c"][0]) == np.sign(v["cg"][0])  # horizontal parts share a sign
        assert v["theta_K"] == pytest.approx(np.arccos(abs(k) / np.hypot(k, m)), rel=1e-15)
    bad = W.internal_wave_velocities(-0.6, 0.9, N, printed=True)  # printed (7.145) for k < 0: wrong direction
    good = W.internal_wave_velocities(-0.6, 0.9, N)
    assert np.sign(bad["cg"][0]) != np.sign(good["cg"][0])
    s = ch07.internal_wave_state(0.6, 1.0, 1.0, k_sign=-1.0, m_sign=1.0)  # K up-left ⇒ c_g down-left (Fig. 7.29)
    assert s["k"] < 0 and s["m"] > 0 and s["cgx"] < 0 and s["cgz"] < 0 and abs(s["dot"]) < 1e-15
    assert s["beam_from_vertical_deg"] == pytest.approx(np.degrees(np.arccos(0.6)))
    assert s["beam_from_horizontal_deg"] == pytest.approx(90 - s["theta_K_deg"]) and s["T"] == pytest.approx(TWO_PI / 0.6)
    assert s["cg"] == pytest.approx(1.0 * np.sin(np.arccos(0.6)) / 1.0, rel=1e-14)  # |c_g| = N sin θ/K (D36 check)
    d = ch07.internal_wave_state(1 / np.sqrt(2), 1.0, 1.0)  # D36 check: 45°, c = (0.5, 0.5), c_g = (0.5, −0.5)
    assert (d["cx"], d["cz"], d["cgx"], d["cgz"]) == pytest.approx((0.5, 0.5, 0.5, -0.5), abs=1e-14)


def test_internal_wave_energy_V1_closed_forms_equal_averages():  # V1 (7.154)–(7.159) vs period averages of the fields
    k, m, N, w0, r0 = 0.5, -1.2, 0.02, 0.01, 1025.0
    e = ch07.internal_wave_energy(k, m, N, w0, r0, G)
    om = e["omega"]
    tt = np.linspace(0, TWO_PI / om, 400, endpoint=False)
    f = ch07.internal_wave_fields(0.3, -0.2, tt, k, m, N, w0, r0, G)
    Ek = np.mean(0.5 * r0 * (f["u"] ** 2 + f["w"] ** 2))
    Ep = np.mean(0.5 * N ** 2 * r0 * f["zeta_particle"] ** 2)  # (7.150)
    assert Ek == pytest.approx(e["Ek"], rel=1e-12) and Ep == pytest.approx(e["Ep"], rel=1e-12)
    assert e["Ek"] == pytest.approx(e["Ep"], rel=1e-13) and e["E"] == pytest.approx(e["Ek"] + e["Ep"], rel=1e-14)
    F = np.array([np.mean(f["p_prime"] * f["u"]), np.mean(f["p_prime"] * f["w"])])
    assert maxrel(F, e["F"]) < 1e-12 and maxrel(e["F"], e["cgE"]) < 1e-13  # (7.158), (7.159)
    assert maxrel(e["cg"], W.internal_wave_velocities(k, m, N)["cg"]) < 1e-15


def test_internal_energy_budget_V4_pointwise_residual():  # V4 (7.147): ∂E_k/∂t + gρ′w + ∇·(p′u) = 0, FD order 2
    X, Z = np.linspace(0, 3, 5), np.linspace(-1, 1, 5)
    hs, errs = [1e-2, 5e-3, 2.5e-3], []
    for h in hs:
        r = ch07.internal_energy_budget_residual(X, Z, 0.4, 1.0, 1.5, 1.0, 0.1, h=h, ht=h, terms=True)
        errs.append(np.max(np.abs(r["residual"])))
        assert np.max(np.abs(r["dKE_dt"])) > 1.0  # the terms themselves are O(1): the balance is not trivial
    assert abs(observed_order(hs, errs) - 2.0) < ORDER_TOL and errs[-1] < 1e-4
    tt = np.linspace(0, TWO_PI / float(W.internal_wave_omega(1.0, 1.5, 1.0)), 64, endpoint=False)
    conv = ch07.internal_energy_budget_residual(0.3, 0.2, tt, 1.0, 1.5, 1.0, 0.1, terms=True)["conversion"]
    assert abs(np.mean(conv)) < 1e-12 * np.max(np.abs(conv))  # ⟨gρ′w⟩ = 0: no net conversion over a period


def test_internal_pe_interface_limit_V3_epsilon_order_one():  # V3 (7.150)–(7.152), N157: smoothed jump → ¼Δρga²
    eps = [4.0, 2.0, 1.0, 0.5, 0.25]
    errs = [abs(ch07.internal_pe_interface_limit(1.0, 1000.0, 1020.0, G, eps=e, detail=True)["rel_error"]) for e in eps]
    assert abs(observed_order(eps, errs) - 1.0) < ORDER_TOL
    d = ch07.internal_pe_interface_limit(1.0, 1000.0, 1020.0, G, eps=0.05, detail=True)
    assert d["target"] == pytest.approx(0.25 * 20 * G) and abs(d["rel_error"]) < 1e-3


def test_st_andrews_cross_V7_beam_geometry():  # V7 N150, Fig. 7.33 (our illustration): beams at arccos(ω/N), c ⟂ c_g
    x = z = np.linspace(-12, 12, 241)
    for r in (0.3, 0.6, 0.9):
        d = ch07.st_andrews_cross(x, z, 0.0, r, 1.0, 1.0, detail=True)
        assert d["theta"] == pytest.approx(np.arccos(r), rel=1e-15)
        for b in d["beams"]:
            eb = b["e_beam"]
            assert np.degrees(np.arctan2(abs(eb[0]), abs(eb[1]))) == pytest.approx(np.degrees(np.arccos(r)), abs=1e-10)
            assert abs(np.dot(b["c"], b["cg"])) < 1e-15 and np.dot(b["cg"], eb) > 0  # energy away from the source
            assert abs(np.dot(b["K"], eb)) < 1e-12  # phase lines along the beam
            assert float(W.internal_wave_omega(b["K"][0], b["K"][1], 1.0)) == pytest.approx(r, rel=1e-12)
        F = d["field"]
        X, Z = np.meshgrid(x, z)
        R = np.hypot(X, Z)
        ring = (R > 6) & (R < 10)
        ang = np.degrees(np.arctan2(np.abs(X[ring]), np.abs(Z[ring])))  # angle from the vertical
        wts = F[ring] ** 2
        hist, edges = np.histogram(ang, bins=90, range=(0, 90), weights=wts)
        assert abs(0.5 * (edges[np.argmax(hist)] + edges[np.argmax(hist) + 1]) - np.degrees(np.arccos(r))) < 1.0
    with pytest.raises(ValueError):
        ch07.st_andrews_cross(x, z, 0.0, 1.5, 1.0, 1.0)


def test_linear_evolve_2d_V1_plane_wave_and_packet_at_cg():  # V1 (7.144)–(7.145), Fig. 7.32 (N198)
    xg = np.linspace(-60, 60, 128, endpoint=False)
    zg = np.linspace(-60, 60, 128, endpoint=False)
    XX, ZZ = np.meshgrid(xg, zg)
    N = 1.0
    wfun = lambda kk, mm: W.internal_wave_omega(kk, mm, N)  # noqa: E731
    k1, m1 = TWO_PI * 5 / 120, TWO_PI * 8 / 120
    om = float(W.internal_wave_omega(k1, m1, N))
    pw = W.linear_evolve_2d(np.exp(1j * (k1 * XX + m1 * ZZ)), xg, zg, 7.3, wfun)
    assert np.max(np.abs(pw - np.cos(k1 * XX + m1 * ZZ - om * 7.3))) < 1e-12
    xg2 = np.linspace(-60, 60, 256, endpoint=False)
    XX2, ZZ2 = np.meshgrid(xg2, xg2)
    f0 = np.exp(-(XX2 ** 2 + ZZ2 ** 2) / (2 * 8.0 ** 2)) * np.exp(1j * (XX2 + ZZ2))
    fT = W.linear_evolve_2d(f0, xg2, xg2, [0.0, 40.0], wfun)
    E = [W.envelope(ff, axis=-1) ** 2 for ff in fT]
    cen = [(np.sum(XX2 * e) / np.sum(e), np.sum(ZZ2 * e) / np.sum(e)) for e in E]
    v = np.array([(cen[1][0] - cen[0][0]) / 40.0, (cen[1][1] - cen[0][1]) / 40.0])
    cg = W.internal_wave_velocities(1.0, 1.0, N)["cg"]
    assert maxrel(v, cg) < 0.01  # the packet slides along its crests at c_g (narrow-band error O(δk²))
    assert np.dot(v, W.internal_wave_velocities(1.0, 1.0, N)["c"]) < 0.01 * np.linalg.norm(v)  # ⟂ to c


# =====================================================================================================================
# Cross-cutting: dimensions (V2), change of units (V7), scripts, Part C coverage
# =====================================================================================================================
def test_dispersion_family_V2_dimensional_homogeneity():  # V2 pint: every term of the implemented relations
    g, k, H = Q_(9.81, "m/s**2"), Q_(0.3, "1/m"), Q_(4.0, "m")
    rho, sig, a = Q_(1000.0, "kg/m**3"), Q_(0.0727, "N/m"), Q_(0.5, "m")
    dimensional_check(lambda g, k, H: (g * k * np.tanh(float((k * H).to("").magnitude))) ** 0.5, "1/[time]", g=g, k=k, H=H)
    dimensional_check(lambda k, g, sig, rho: (k * (g + sig * k ** 2 / rho)) ** 0.5, "1/[time]", k=k, g=g, sig=sig, rho=rho)
    dimensional_check(lambda g, k, sig, rho: (4 * g * sig / rho) ** 0.25, "velocity", g=g, k=k, sig=sig, rho=rho)
    dimensional_check(lambda sig, rho, g: 2 * np.pi * (sig / (rho * g)) ** 0.5, "[length]", sig=sig, rho=rho, g=g)
    dimensional_check(lambda rho, g, a: rho * g * a ** 2 / 2, "[mass]/[time]**2", rho=rho, g=g, a=a)  # E [J/m²]
    dimensional_check(lambda rho, g, a, k: rho * g * a ** 2 / 2 * (g / k) ** 0.5 / 2, "[mass]*[length]/[time]**3",
                      rho=rho, g=g, a=a, k=k)  # F [W/m]
    dimensional_check(lambda a, g, k: a ** 2 * (g * k) ** 0.5 * k, "velocity", a=a, g=g, k=k)  # (7.85)
    dimensional_check(lambda g, rho: g * (Q_(1002.0, "kg/m**3") - rho) / Q_(1002.0, "kg/m**3"), "acceleration",
                      g=g, rho=rho)  # (7.117)
    N, w0 = Q_(0.01, "1/s"), Q_(0.02, "m/s")
    dimensional_check(lambda rho, N, w0, k: rho * (N * 1.0) * Q_(1.2, "1/m") * w0 ** 2 / (2 * k ** 2),
                      "[mass]/[time]**3", rho=rho, N=N, w0=w0, k=k)  # (7.158) [W/m²]
    dimensional_check(lambda sig, k, a: sig * a * k ** 2, "pressure", sig=sig, k=k, a=a)  # (7.54)


def test_wave_functions_V7_change_of_units():  # V7: rescaling length, time and mass rescales every output consistently
    Lf, Tf, Mf = 100.0, 10.0, 1000.0  # m → cm, s → ds, kg → g
    gS, rS, sS = G, 1000.0, 0.0727
    gP, rP, sP = G * Lf / Tf ** 2, rS * Mf / Lf ** 3, sS * Mf / Tf ** 2
    k, H, a = 0.7, 3.0, 0.2
    kP, HP, aP = k / Lf, H * Lf, a * Lf
    V, Om = Lf / Tf, 1 / Tf
    assert float(W.omega_gravity(kP, HP, gP)) == pytest.approx(float(W.omega_gravity(k, H, gS)) * Om, rel=1e-13)
    assert float(W.phase_speed(kP * 300, HP, gP, sP, rP)) == pytest.approx(float(W.phase_speed(k * 300, H, gS, sS, rS)) * V,
                                                                           rel=1e-13)
    assert float(W.group_velocity(kP * 300, HP, gP, sP, rP)) == pytest.approx(
        float(W.group_velocity(k * 300, H, gS, sS, rS)) * V, rel=1e-13)
    mS, mP = W.capillary_minimum(sS, rS, gS), W.capillary_minimum(sP, rP, gP)
    assert mP["c_min"] == pytest.approx(mS["c_min"] * V, rel=1e-13) and mP["lam_m"] == pytest.approx(mS["lam_m"] * Lf)
    assert ch07.wave_energy(aP, kP, HP, gP, rP)["E"] == pytest.approx(ch07.wave_energy(a, k, H, gS, rS)["E"] * Mf / Tf ** 2,
                                                                      rel=1e-13)
    assert ch07.energy_flux(aP, kP, HP, gP, rP) == pytest.approx(ch07.energy_flux(a, k, H, gS, rS) * Mf * Lf / Tf ** 3,
                                                                 rel=1e-13)
    assert float(ch07.stokes_drift(-0.4 * Lf, aP, kP, HP, gP)) == pytest.approx(
        float(ch07.stokes_drift(-0.4, a, k, H, gS)) * V, rel=1e-13)
    jS, jP = ch07.hydraulic_jump(0.3, Fr1=2.2, g=gS), ch07.hydraulic_jump(0.3 * Lf, Fr1=2.2, g=gP)
    assert jP["H2"] == pytest.approx(jS["H2"] * Lf, rel=1e-14) and jP["dE"] == pytest.approx(jS["dE"] * V ** 2, rel=1e-13)
    assert float(W.interface_omega(kP, rP, 1.02 * rP, gP)) == pytest.approx(float(W.interface_omega(k, rS, 1.02 * rS, gS))
                                                                            * Om, rel=1e-13)
    assert float(W.two_layer_long_wave_speed(HP, rP, 1.01 * rP, gP)) == pytest.approx(
        float(W.two_layer_long_wave_speed(H, rS, 1.01 * rS, gS)) * V, rel=1e-13)
    eS = ch07.internal_wave_energy(k, 1.3 * k, 0.02, 0.01, rS, gS)
    eP = ch07.internal_wave_energy(kP, 1.3 * kP, 0.02 / Tf, 0.01 * V, rP, gP)
    assert eP["E"] == pytest.approx(eS["E"] * Mf / (Lf * Tf ** 2), rel=1e-13)
    assert maxrel(eP["F"], eS["F"] * Mf / Tf ** 3) < 1e-13
    assert float(ch07.seiche_modes(20 * Lf, HP, 1, gP)["T"]) == pytest.approx(float(ch07.seiche_modes(20, H, 1, gS)["T"]) * Tf,
                                                                              rel=1e-13)
    assert float(W.viscous_decay(aP, kP, 1.3e-6 * Lf ** 2 / Tf, 50.0 * Tf)) == pytest.approx(
        float(W.viscous_decay(a, k, 1.3e-6, 50.0)) * Lf, rel=1e-13)
    assert ch07.solitary_wave_speed(aP, HP, gP) == pytest.approx(ch07.solitary_wave_speed(a, H, gS) * V, rel=1e-14)


@pytest.mark.slow
def test_scripts_V1_every_ch07_script_runs():  # V1 smoke (I51): all 11 scripts exit 0 headless
    env = dict(os.environ, MPLBACKEND="Agg")
    scripts = [s for s in sorted((ROOT / "scripts").glob("ch07_*.py")) if s.name != "ch07_drawings.py"]
    assert len(scripts) == 11
    for scr in scripts:
        r = subprocess.run([sys.executable, str(scr), "--no-show"], cwd=ROOT, env=env, capture_output=True, text=True,
                           timeout=600)
        assert r.returncode == 0, (scr.name, r.stderr[-2000:])


def test_part_c_V1_every_contract_function_exists_and_is_scalar_callable():  # V1 smoke: design Part C C.1–C.2
    names = ["sinusoid", "wave_parameters", "crest_positions", "plane_wave", "phase_velocity_vector", "trace_velocities",
             "doppler_frequency", "omega_gravity", "omega_capillary_gravity", "phase_speed", "period_from_wavelength",
             "wavelength_from_period", "wavenumber_from_omega", "fenton_mckee_kh", "guo_kh", "group_velocity",
             "group_velocity_numeric", "group_velocity_vector", "depth_regime", "capillary_minimum", "min_group_velocity",
             "beat_wave", "gaussian_packet", "linear_evolve", "envelope", "linear_evolve_2d", "ray_trace",
             "wave_energy_density", "viscous_decay", "interface_omega", "two_layer_free_surface_omega",
             "two_layer_long_wave_speed", "reduced_gravity_book", "internal_wave_omega", "beam_angle",
             "internal_wave_velocities", "real_field",
             "surface_normal", "surface_velocity", "linear_bernoulli_pressure", "surface_wave_sympy", "wave_fields",
             "free_surface_residuals", "free_surface_residual_scan", "pressure_response", "dispersion_state",
             "particle_path", "orbit_linear", "orbit_semi_axes", "orbit_state", "dyed_line", "wave_energy", "energy_flux",
             "curvature", "capillary_surface_pressure", "capillary_state", "standing_wave_fields", "seiche_modes",
             "basin_modes", "seiche_state", "packet_state", "pond_ripples", "local_wavenumber_frequency",
             "crest_conservation_residual", "snell_ray_plane_beach", "refraction_state", "nonlinear_wavelet_speed",
             "simple_wave_evolve", "hydraulic_jump", "jump_momentum_residual", "jump_state", "stokes_wave_profile",
             "stokes_wave_speed", "stokes_expansion_sympy", "stokes_drift", "stokes_drift_numeric", "eulerian_mean_u",
             "kdv_rhs", "kdv_solve", "kdv_invariants", "kdv_linear_phase_speed", "ursell_number", "solitary_wave",
             "cnoidal_wave", "kdv_residual_sympy", "interface_fields", "interface_residuals", "interface_energy",
             "two_layer_modes", "two_layer_residuals", "two_layer_sympy", "two_layer_rigid_lid_omega",
             "two_layer_state", "boussinesq_linear_sympy", "internal_wave_fields", "w_equation_residual",
             "layered_flow_check", "internal_wave_energy", "internal_energy_budget_residual",
             "internal_pe_interface_limit", "st_andrews_cross", "internal_wave_state", "wave_regime_label"]
    missing = [n for n in names if not callable(getattr(ch07, n, None))]
    assert not missing, missing
    assert ch07.STOKES_LIMIT_STEEPNESS == 0.1410633 and ch07.G_BOOK == 9.81
    for fn, args in ((ch07.dispersion_state, (40.0, 8.0)), (ch07.capillary_state, (0.02,)),
                     (ch07.seiche_state, (100.0, 3.0)), (ch07.packet_state, (0.1,)),
                     (ch07.refraction_state, (9.0, 0.4, 25.0, 3.0)), (ch07.jump_state, (0.2, 2.0)),
                     (ch07.two_layer_state, (0.05, 30.0)), (ch07.internal_wave_state, (0.5,)),
                     (ch07.orbit_state, (-0.3,))):
        out = fn(*args)
        assert all(np.ndim(v) == 0 for v in out.values()), fn.__name__  # scalar-callable (explainer parity rows)
        assert any(isinstance(v, float) and np.isfinite(v) for v in out.values())
    for fn, args in ((W.omega_gravity, (1.0, 2.0)), (W.phase_speed, (1.0,)), (W.group_velocity, (1.0, 2.0)),
                     (ch07.stokes_drift, (-0.2,)), (ch07.pressure_response, (1.0, -0.5, 2.0)),
                     (W.internal_wave_omega, (1.0, 1.0, 1.0)), (W.beam_angle, (0.5, 1.0))):
        assert isinstance(fn(*args), float)


# =====================================================================================================================
# V6 — the book's printed forms and numbers (private JSON; skipped when absent)
# =====================================================================================================================
def _bsyms():
    names = ["a", "k", "m", "l", "K", "omega", "g", "H", "rho", "rho1", "rho2", "rho0", "sigma", "N", "x", "z", "z0", "x0",
             "t", "L", "n", "Fr1", "H1", "H2", "lam", "c", "c0", "w0"]
    return {n: sp.Symbol(n, real=True) for n in names}


@book_only
def test_book_V6_section_7_2_forms_and_numbers():  # V6 §7.2: (7.26)–(7.52) as printed; accuracy claims; ocean numbers
    b = book()
    f = b["sec7_2_closed_forms"]
    S = _bsyms()
    a, k, H, x, z, t, x0, z0 = 0.3, 0.45, 5.0, 2.2, -1.7, 0.9, 1.1, -2.3
    om = float(W.omega_gravity(k, H, G))
    val = {S["a"]: a, S["k"]: k, S["H"]: H, S["g"]: G, S["omega"]: om, S["x"]: x, S["z"]: z, S["t"]: t, S["rho"]: 1000.0,
           S["x0"]: x0, S["z0"]: z0, S["lam"]: TWO_PI / k, S["c"]: om / k}
    ev = lambda key: float(sp.sympify(f[key], locals=S).subs(val))  # noqa: E731
    wf = ch07.wave_fields(x, z, t, a, k, H, G)
    for key, ours in (("phi", wf["phi"]), ("u", wf["u"]), ("w", wf["w"]), ("psi", wf["psi"]), ("p_prime", wf["p_prime"]),
                      ("dispersion_omega", om), ("period_from_wavelength", TWO_PI / om), ("phase_speed", om / k),
                      ("E_total", ch07.wave_energy(a, k, H, G)["E"]), ("energy_flux", ch07.energy_flux(a, k, H, G))):
        assert ev(key) == pytest.approx(float(ours), rel=1e-12), key
    xi, ze = ch07.orbit_linear(x0, z0, t, a, k, H, G)
    assert ev("xi") == pytest.approx(float(xi), rel=1e-12) and ev("zeta") == pytest.approx(float(ze), rel=1e-12)
    n = b["sec7_2_numbers"]
    assert np.tanh(2.0) == pytest.approx(n["tanh_2"], abs=5e-6)
    assert 100 * W.depth_regime(n["deep_threshold_kH"], 1.0)["deep_error"] < n["deep_accuracy_percent"]
    assert 1 / np.pi == pytest.approx(n["deep_threshold_H_over_lambda"], rel=0.006)  # ours 1/π = 0.3183; the book rounds
    se = 100 * W.depth_regime(TWO_PI * n["shallow_threshold_H_over_lambda"], 1.0)["shallow_error"]
    assert round(se) == n["shallow_accuracy_percent"]  # ours 3.04 %; the book rounds
    assert 1 / n["shallow_threshold_H_over_lambda"] == pytest.approx(n["shallow_length_ratio"], rel=0.03)
    lam10 = float(W.wavelength_from_period(n["wind_wave_period_s"], np.inf, G))
    assert 0.0 < lam10 / n["wind_wave_wavelength_m"] - 1 < 0.05  # the book rounds the deep-water λ down (analysis T9)
    for Hd in (n["shelf_depth_m"], n["open_ocean_depth_m"]):  # wind waves are deep water there
        kk = float(W.wavenumber_from_omega(TWO_PI / n["wind_wave_period_s"], Hd, G))
        assert W.depth_regime(kk, Hd)["regime"] == "deep"
    assert np.exp(-np.pi) == pytest.approx(n["deep_pressure_fraction_at_half_wavelength"], abs=0.005)
    lim = {S["a"]: a, S["k"]: k, S["g"]: G, S["H"]: 1e-4, S["omega"]: float(W.omega_gravity(k, 1e-4, G)), S["x"]: x,
           S["z"]: -5e-5, S["t"]: t, S["rho"]: 1000.0, S["x0"]: x0, S["z0"]: -5e-5}
    assert float(sp.sympify(f["shallow_xi"], locals=S).subs(lim)) == pytest.approx(
        float(ch07.orbit_linear(x0, -5e-5, t, a, k, 1e-4, G)[0]), rel=1e-7)
    assert float(sp.sympify(f["shallow_zeta"], locals=S).subs(lim)) == pytest.approx(
        float(ch07.orbit_linear(x0, -5e-5, t, a, k, 1e-4, G)[1]), rel=1e-7)
    assert float(sp.sympify(f["deep_phase_speed"], locals=S).subs(val)) == pytest.approx(np.sqrt(G / k))
    assert float(sp.sympify(f["shallow_phase_speed"], locals=S).subs(val)) == pytest.approx(np.sqrt(G * H))
    dp = {**val, S["omega"]: np.sqrt(G * k)}
    assert float(sp.sympify(f["deep_p_prime"], locals=S).subs(dp)) == pytest.approx(
        float(ch07.wave_fields(x, z, t, a, k, np.inf, G)["p_prime"]), rel=1e-12)


@book_only
def test_book_V6_capillary_standing_and_group_numbers():  # V6 §7.3–7.5: (7.56)–(7.60), (7.62)–(7.65), (7.69), Ex. 7.10
    b = book()
    S = _bsyms()
    c3 = b["sec7_3_capillary"]
    sig, rho, k, H = 0.07, 1000.0, 250.0, 0.03
    val = {S["sigma"]: sig, S["rho"]: rho, S["k"]: k, S["H"]: H, S["g"]: G, S["lam"]: TWO_PI / k}
    ev = lambda d, key, v=val: float(sp.sympify(d[key], locals=S).subs(v))  # noqa: E731
    assert ev(c3, "omega") == pytest.approx(float(W.omega_capillary_gravity(k, H, sig, rho, G)), rel=1e-13)
    m = W.capillary_minimum(sig, rho, G)
    assert ev(c3, "c_min") == pytest.approx(m["c_min"], rel=1e-13) and ev(c3, "lambda_m") == pytest.approx(m["lam_m"])
    assert ev(c3, "pure_capillary_c") == pytest.approx(float(W.phase_speed(k, np.inf, 0.0, sig, rho)), rel=1e-13)
    mb = W.capillary_minimum(c3["sigma_20C_N_per_m"], 1000.0, G)
    assert 100 * mb["c_min"] == pytest.approx(c3["c_min_cm_per_s"], rel=0.005)  # (7.59): ours 23.13 cm/s
    assert 100 * mb["lam_m"] == pytest.approx(c3["lambda_m_cm"], rel=0.005)  # ours 1.714 cm
    s4 = b["sec7_4_standing"]
    a, x, z, t, L = 0.1, 0.8, -0.2, 1.3, 3.0
    om = float(W.omega_gravity(1.1, 0.7, G))
    v4 = {S["a"]: a, S["k"]: 1.1, S["H"]: 0.7, S["omega"]: om, S["x"]: x, S["z"]: z, S["t"]: t, S["g"]: G, S["L"]: L}
    sw = ch07.standing_wave_fields(x, z, t, a, 1.1, 0.7, G)
    assert ev(s4, "psi_standing", v4) == pytest.approx(float(sw["psi"]), rel=1e-12)
    assert ev(s4, "u_standing", v4) == pytest.approx(float(sw["u"]), rel=1e-12)
    for nn in range(3):
        v4[S["n"]] = nn
        sm = ch07.seiche_modes(L, 0.7, nn, G)
        assert ev(s4, "allowed_wavelengths", v4) == pytest.approx(float(sm["lam"]))
        assert ev(s4, "natural_frequencies", v4) == pytest.approx(float(sm["omega"]), rel=1e-13)
    s5 = b["sec7_5_group"]
    v5 = {S["k"]: 0.4, S["H"]: 3.0, S["c"]: float(W.phase_speed(0.4, 3.0, G))}
    assert ev(s5, "cg_surface", v5) == pytest.approx(float(W.group_velocity(0.4, 3.0, G)), rel=1e-13)
    assert float(W.group_velocity(1.0, np.inf, G) / W.phase_speed(1.0, np.inf, G)) == s5["cg_deep_over_c"]
    assert float(W.group_velocity(1e-5, 1.0, G) / W.phase_speed(1e-5, 1.0, G)) == pytest.approx(s5["cg_shallow_over_c"])
    assert float(W.group_velocity(1e4, np.inf, 0.0, 0.07) / W.phase_speed(1e4, np.inf, 0.0, 0.07)) == pytest.approx(
        s5["cg_capillary_over_c"], rel=1e-12)
    lm = sp.Symbol("lam_max", real=True)
    cgmax = float(sp.sympify(s5["cg_max_stone"], locals={**S, "lam_max": lm}).subs({S["g"]: G, lm: 0.5}))
    assert cgmax == pytest.approx(float(W.group_velocity(TWO_PI / 0.5, np.inf, G)), rel=1e-14)  # longest waves lead
    e10 = b["exercises"]["ex7_10"]
    cgm = W.min_group_velocity(e10["sigma"], e10["rho"], G)
    assert 100 * cgm["cg_min"] == pytest.approx(e10["cg_min_cm_per_s"], rel=0.005)  # ours 17.83 cm/s
    assert b["exercises"]["ex7_9_cg_over_c_capillary"] == 1.5
    bw = W.beat_wave(0.0, 0.0, 0.9, 1.1)
    dks = sp.Symbol("dk", real=True)
    mw = float(sp.sympify(s5["modulation_wavelength"], locals={**S, "dk": dks}).subs(dks, bw["dk"]))
    assert mw == pytest.approx(4 * np.pi / 0.2, rel=1e-12)  # two beats per envelope wavelength


@book_only
def test_book_V6_nonlinear_interface_and_internal():  # V6 §7.6–7.8 forms; Figs. 7.33; Exercises 7.6, 7.16, 7.17, 7.20
    b = book()
    S = _bsyms()
    s6 = b["sec7_6_nonlinear"]
    Fr1, H1 = 2.7, 0.3
    j = ch07.hydraulic_jump(H1, Fr1=Fr1, g=G)
    v = {S["Fr1"]: Fr1, S["H1"]: H1, S["H2"]: j["H2"], S["g"]: G}
    ev = lambda d, key, vv: float(sp.sympify(d[key], locals=S).subs(vv))  # noqa: E731
    assert ev(s6, "jump_depth_ratio", v) == pytest.approx(j["ratio"], rel=1e-14)
    assert ev(s6, "jump_energy_change", v) == pytest.approx(j["dE"], rel=1e-13)
    a, k, x, t = 0.08, 1.2, 0.4, 0.7
    c = float(ch07.stokes_wave_speed(k, a, G))
    vs = {S["a"]: a, S["k"]: k, S["x"]: x, S["t"]: t, S["c"]: c, S["g"]: G}
    assert ev(s6, "stokes_profile", vs) == pytest.approx(float(ch07.stokes_wave_profile(x, t, a, k, G)), rel=1e-13)
    assert ev(s6, "stokes_speed", vs) == pytest.approx(c, rel=1e-14)
    assert 0.5 * ch07.STOKES_LIMIT_STEEPNESS == pytest.approx(s6["stokes_max_amplitude_over_lambda"], rel=0.01)
    vd = {S["a"]: a, S["k"]: k, S["omega"]: float(W.omega_gravity(k, 2.0, G)), S["z0"]: -0.5, S["H"]: 2.0}
    assert ev(s6, "stokes_drift_general", vd) == pytest.approx(float(ch07.stokes_drift(-0.5, a, k, 2.0, G)), rel=1e-13)
    vd[S["omega"]] = np.sqrt(G * k)
    assert ev(s6, "stokes_drift_deep", vd) == pytest.approx(float(ch07.stokes_drift(-0.5, a, k, np.inf, G)), rel=1e-13)
    Hd = 1.3
    vk = {S["a"]: 0.2, S["H"]: Hd, S["x"]: 0.9, S["t"]: 0.4, S["c0"]: np.sqrt(G * Hd), S["g"]: G, S["k"]: 0.6,
          S["c"]: ch07.solitary_wave_speed(0.2, Hd, G)}
    assert ev(s6, "soliton", vk) == pytest.approx(float(ch07.solitary_wave(0.9, 0.4, 0.2, Hd, G)), rel=1e-13)
    assert ev(s6, "soliton_speed", vk) == pytest.approx(ch07.solitary_wave_speed(0.2, Hd, G), rel=1e-14)
    assert ev(s6, "kdv_linear_phase_speed", vk) == pytest.approx(float(ch07.kdv_linear_phase_speed(0.6, Hd, G)))
    s7 = b["sec7_7_interface"]
    r1, r2, k, a, H = 1000.0, 1003.0, 0.2, 0.4, 30.0
    vi = {S["rho1"]: r1, S["rho2"]: r2, S["k"]: k, S["g"]: G, S["a"]: a, S["H"]: H}
    assert ev(s7, "omega_interface", vi) == pytest.approx(float(W.interface_omega(k, r1, r2, G)), rel=1e-14)
    assert ev(s7, "Ek_interface", vi) == pytest.approx(ch07.interface_energy(a, k, r1, r2, G)["Ek"], rel=1e-14)
    assert ev(s7, "E_interface", vi) == pytest.approx(ch07.interface_energy(a, k, r1, r2, G)["E"], rel=1e-14)
    wbt, wbc = W.two_layer_free_surface_omega(k, H, r1, r2, G)
    assert ev(s7, "barotropic_root", vi) == pytest.approx(float(wbt) ** 2) and ev(s7, "baroclinic_root", vi) == pytest.approx(
        float(wbc) ** 2, rel=1e-12)
    assert ev(s7, "barotropic_b", vi) == pytest.approx(np.real(ch07.two_layer_modes(k, H, r1, r2, G, a, "barotropic")["b"]))
    assert ev(s7, "baroclinic_eta_over_zeta", vi) == pytest.approx(ch07.two_layer_modes(k, H, r1, r2, G, a)["eta_over_zeta"],
                                                                   rel=1e-9)
    assert ev(s7, "long_wave_speed", vi) == pytest.approx(float(W.two_layer_long_wave_speed(H, r1, r2, G)))
    d20, d30 = float(ch01.water_density(293.15)), float(ch01.water_density(303.15))
    assert round(100 * (d20 - d30) / d20, 1) == s7["density_decrease_percent_per_10C"]
    s8 = b["sec7_8_internal"]
    k, m, l, N, w0, r0 = 0.3, 0.7, 0.2, 0.01, 0.02, 1025.0
    K = np.sqrt(k ** 2 + m ** 2)
    om = float(W.internal_wave_omega(k, m, N))
    vn = {S["k"]: k, S["m"]: m, S["l"]: l, S["N"]: N, S["w0"]: w0, S["rho0"]: r0, S["omega"]: om, S["K"]: K}
    assert ev(s8, "dispersion", vn) == pytest.approx(float(W.internal_wave_omega(k, m, N, l=l)) ** 2, rel=1e-13)
    e = ch07.internal_wave_energy(k, m, N, w0, r0, G)
    assert ev(s8, "Ek", vn) == pytest.approx(e["Ek"], rel=1e-13) and ev(s8, "Ep", vn) == pytest.approx(e["Ep"], rel=1e-13)
    assert float(np.degrees(W.beam_angle(s8["fig7_33_omega_rad_s"], s8["fig7_33_N_rad_s"]))) == pytest.approx(
        s8["fig7_33_angle_deg"], rel=0.01)  # ours 44.77°: the figure label is rounded
    ex = b["exercises"]
    e6 = ex["ex7_6_lake"]
    T6 = float(ch07.basin_modes(e6["L_m"], e6["b_m"], e6["H_m"], *e6["mode"], g=G)["T"]) / 60
    assert 0.0 < T6 / e6["period_min"] - 1 < 0.01  # T8: (7.65) is verified independently; the printed answer is low
    e16 = ex["ex7_16"]
    assert np.sqrt(e16["g_prime"] * e16["thermocline_depth_m"]) == pytest.approx(e16["c_m_s"], rel=0.005)
    gp = W.reduced_gravity_book(float(ch01.water_density(273.15 + e16["T_upper_C"])),
                                float(ch01.water_density(273.15 + e16["T_lower_C"])), G)
    assert gp == pytest.approx(e16["g_prime"], rel=0.02)  # the book's density table vs Kell (1975)
    e17 = ex["ex7_17"]
    assert np.degrees(float(W.beam_angle(e17["omega"], e17["N"]))) == pytest.approx(60.0, rel=1e-13)
    h1, h2 = 20.0, 70.0
    vr = {S["g"]: G, S["k"]: 0.1, S["rho1"]: r1, S["rho2"]: r2, sp.Symbol("h1", real=True): h1,
          sp.Symbol("h2", real=True): h2}
    form = sp.sympify(ex["ex7_20_rigid_lid_dispersion"], locals={**S, "h1": sp.Symbol("h1", real=True),
                                                                 "h2": sp.Symbol("h2", real=True)})
    assert float(form.subs(vr)) == pytest.approx(float(ch07.two_layer_rigid_lid_omega(0.1, h1, h2, r1, r2, G)) ** 2,
                                                 rel=1e-13)
