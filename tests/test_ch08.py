"""Verification suite for Chapter 8 — Laminar Flow (Kundu, Cohen & Dowling 5e, §§8.1–8.7).

Evidence levels (``verify-implementation`` skill): V1 analytic · V2 symbolic (sympy / dimensions) · V3 convergence ·
V4 conservation / invariant · V5 published benchmark (``reference/ch08/``, see SOURCES.md) · V6 book value (private,
``tests/book_values_ch08.json``, skipped when absent) · V7 limits / symmetry / invariance. Every test name is
``test_<concept>_<level>_<what>``; the comment on the ``def`` line repeats the level.

A items (CORE, ≥ 2 independent levels): C01 laminar vs turbulent, Re = Ud/ν, ν = μ/ρ · C02 Couette–Poiseuille (8.5) ·
C03 pipe Poiseuille (8.6)–(8.8) · C04 circular Couette (8.9)–(8.12) · C05 lubrication balance (8.14)–(8.17) · C06
lubrication profile (8.18)–(8.19), gap flux, Reynolds equation, Hele-Shaw · C07 slider bearing (Example 8.1) · C08 thin
film (Example 8.3) · C09 Stokes' first problem (8.20)–(8.31) · C10 the similarity ansatz (8.32), Examples 8.4–8.7 ·
C11 Stokes' second problem (8.33)–(8.38) · C12 Stokes equations (8.39)–(8.43) · C13 Stokes' stream function
(8.44)–(8.49) · C14 Stokes drag (8.50)–(8.52), settling, Millikan · C15 far-field breakdown and Oseen (8.53).
Derivations: every ★★ and ★★★ D row (D04, D06–D08, D10–D15, D17–D19, D21–D24, D26–D33) is re-derived with sympy in a
``test_*_V2_derivation`` test; the seven ★★★ rows (D10, D13, D14, D22, D28, D30, D31) step by step through the design's
Part F lines.

Pinned conventions and slips (each with a discriminating assertion): walls y = 0 (fixed) and h (moving); ``dpdx`` is the
book's dp/dx while ch04's G = −dp/dx; η = y/√(νt); θ from the downstream axis in §8.6; Re = 2aU/ν for the sphere.
Printed slips kept as wrong variants that must FAIL: (8.13b) with ∂p/∂x; (8.17a) without ν; (8.19) with U₀ added;
V = Q/h without 1/h; Example 8.1's (1 − αx/L) integrands and first-power final denominator; ∫ω dy = −U; η₉₅ = ±2.76; the
rear pressure minimum +3μU/2a; Oseen's equation with +∂p/∂x_i; ch05's ω = Ω₁ (factor 2) and γ = +2U; the biharmonic
∇⁴ψ in place of (E²)²ψ; Oseen's 3/8 applied to the diameter Re.

Run: ``.venv/Scripts/python.exe -m pytest tests/test_ch08.py -q -p no:cacheprovider``  (``-m "not slow"`` skips the
script smoke test).
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
from scipy.optimize import brentq
from scipy.special import erf, erfc, erfcinv, erfinv, gamma as gamma_fn

from fluidpy import ch01_introduction as ch01
from fluidpy import ch03_kinematics as ch03
from fluidpy import ch04_conservation_laws as ch04
from fluidpy import ch05_vorticity_dynamics as ch05
from fluidpy import ch08_laminar_flow as ch08
from fluidpy.core import creeping as CRP
from fluidpy.core import curvilinear as CL
from fluidpy.core import diffusion as DIF
from fluidpy.core import laminar as LAM
from fluidpy.core import lubrication as LUB
from fluidpy.core import navier_stokes as NS
from fluidpy.core import potential as PF
from fluidpy.core import similarity as SIM
from fluidpy.core import vortices as VX
from fluidpy.core.units import Q_
from tools.convergence import observed_order, pairwise_orders

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference" / "ch08"
BOOK = Path(__file__).resolve().parent / "book_values_ch08.json"
needs_ref = pytest.mark.skipif(not (REF / "benchmarks.json").exists(), reason="run reference/ch08/make_refs.py")
book_only = pytest.mark.skipif(not BOOK.exists(), reason="book values are private; see CLAUDE.md rule 9")

RNG = np.random.default_rng(8)
ORDER_TOL = 0.15  # design order ± this (verify-implementation default)
EPS = np.finfo(float).eps


def book():
    return json.loads(BOOK.read_text(encoding="utf-8"))


def ref_json():
    return json.loads((REF / "benchmarks.json").read_text(encoding="utf-8"))


def maxrel(num, ref, floor=0.0):
    num, ref = np.asarray(num, dtype=float), np.asarray(ref, dtype=float)
    return float(np.max(np.abs(num - ref)) / max(np.max(np.abs(ref)), floor, 1e-300))


def z0(expr) -> bool:
    return sp.simplify(expr) == 0


# =====================================================================================================================
# C01 — laminar vs turbulent: Re = Ud/ν, the momentum diffusivity ν = μ/ρ, t ~ L²/ν (§8.1)
# =====================================================================================================================
def test_pipe_flow_regime_V1_definition_and_band_edges():  # V1 (C01): Re = Ud/ν on a field; labels switch at 2000/3000
    U = np.linspace(0.01, 0.5, 60)
    d, nu = 0.02, 1.0e-6
    Re, lab = ch08.pipe_flow_regime(U, d, nu)
    assert maxrel(Re, U * d / nu) < 1e-15
    exp = np.where(U * d / nu < 2000, "laminar", np.where(U * d / nu <= 3000, "transitional", "turbulent"))
    assert np.all(lab == exp)
    assert np.all(np.diff(Re) > 0)  # monotone in U
    # exactly at the band edges (the book's 2000 and 3000)
    assert ch08.pipe_flow_regime(0.1, d, nu)[1] == "transitional"  # Re = 2000 exactly
    assert ch08.pipe_flow_regime(0.1 * (1 - 1e-12), d, nu)[1] == "laminar"
    assert ch08.pipe_flow_regime(0.15, d, nu)[1] == "transitional"  # Re = 3000
    assert ch08.pipe_flow_regime(0.15 * (1 + 1e-12), d, nu)[1] == "turbulent"
    # wrong variant: the radius instead of the diameter halves Re and changes the regime of a 2500 flow
    assert ch08.pipe_flow_regime(0.125, d / 2, nu)[1] == "laminar" != ch08.pipe_flow_regime(0.125, d, nu)[1]


def test_pipe_flow_regime_V7_scale_invariance():  # V7 (C01): Re unchanged under U→kU, d→d/k and under a change of units
    for k in (0.1, 3.0, 17.0):
        assert ch08.pipe_flow_regime(0.3 * k, 0.01 / k, 1e-6)[0] == pytest.approx(ch08.pipe_flow_regime(0.3, 0.01, 1e-6)[0],
                                                                                   rel=1e-14)
    # the same flow in cm and s (U [cm/s], d [cm], ν [cm²/s])
    assert ch08.pipe_flow_regime(30.0, 1.0, 1e-2)[0] == pytest.approx(ch08.pipe_flow_regime(0.3, 0.01, 1e-6)[0], rel=1e-14)
    s = ch08.inertia_viscous_scales(2.0, 0.5, 1e-5)
    assert s["ratio"] == pytest.approx(s["inertia"] / s["viscous"], rel=1e-14)
    assert s["ratio"] == pytest.approx(float(SIM.reynolds_number(2.0, 0.5, 1e-5)), rel=1e-14)  # reuse ch04 C15


def test_momentum_diffusivity_V5_air_water_from_published_property_laws():  # V5 (C01/N01): USSA Sutherland + ch01 water
    consts = json.loads((ROOT / "reference" / "ch01" / "ussa1976_constants.json").read_text(encoding="utf-8"))
    beta, S = consts["sutherland_beta_kg_per_m_s_sqrtK"], consts["sutherland_S_K"]  # USSA-1976 (V5 in ch01)
    T = np.array([273.15, 293.15, 313.15])
    for Ti in T:
        mu_air = beta * Ti ** 1.5 / (Ti + S)  # independent Sutherland evaluation from the cited constants
        rho_air = ch08.P_ATM / (ch08.R_AIR * Ti)
        assert ch08.momentum_diffusivity("air", Ti) == pytest.approx(mu_air / rho_air, rel=1e-12)
        assert ch08.momentum_diffusivity("water", Ti) == pytest.approx(float(ch01.water_viscosity(Ti)) /
                                                                       float(ch01.water_density(Ti)), rel=1e-12)
    ratio = ch08.momentum_diffusivity("air", 293.15) / ch08.momentum_diffusivity("water", 293.15)
    assert 14.0 < ratio < 17.0  # analysis V2 band; D01's "air ≈ 15× water"
    assert ch08.momentum_diffusivity("air", 293.15) == pytest.approx(1.5e-5, rel=0.02)
    assert ch08.momentum_diffusivity("water", 293.15) == pytest.approx(1.0e-6, rel=0.02)
    # wrong variant: comparing μ instead of ν inverts the verdict (water is ~55× more viscous)
    assert float(ch01.water_viscosity(293.15)) / float(ch01.sutherland_viscosity(293.15)) > 50
    with pytest.raises(ValueError):
        ch08.momentum_diffusivity("oil")


def test_diffusion_time_V1_D01_numbers():  # V1 (C01/D01): t = L²/ν; 1 cm in water 100 s, in air ≈ 6.7 s
    L = np.array([1e-3, 1e-2, 0.1])
    assert maxrel(ch08.diffusion_time(L, 1e-6), L ** 2 / 1e-6) < 1e-15
    assert ch08.diffusion_time(0.01, 1e-6) == pytest.approx(100.0, rel=1e-14)
    assert ch08.diffusion_time(0.01, 1.5e-5) == pytest.approx(6.667, rel=1e-3)
    s = ch08.inertia_viscous_scales(0.2, 0.05, 1e-6)
    assert s["inertia"] == pytest.approx(0.2 ** 2 / 0.05) and s["viscous"] == pytest.approx(1e-6 * 0.2 / 0.05 ** 2)


def test_wall_bc_residuals_V1_no_through_flow_and_no_slip():  # V1 (R04, R05): (8.2), (8.3)
    n, t = np.array([0.0, 1.0, 0.0]), np.array([1.0, 0.0, 0.0])
    Us = np.array([0.3, 0.0, 0.0])
    assert ch08.wall_bc_residuals(Us, Us, n) == (0.0, 0.0)
    assert ch08.wall_bc_residuals(Us + [0.0, 0.2, 0.0], Us, n) == pytest.approx((0.2, 0.0))  # through-flow only
    assert ch08.wall_bc_residuals(Us + [0.1, 0.0, 0.0], Us, n, t) == pytest.approx((0.0, 0.1))  # slip only
    # tangential part = magnitude of the in-plane slip when no tangent is given (wrong variant: t·u ignored)
    assert ch08.wall_bc_residuals(Us + [0.3, 0.0, 0.4], Us, n)[1] == pytest.approx(0.5)
    # a tilted wall with an unnormalised normal
    nn = np.array([1.0, 1.0, 0.0])
    assert ch08.wall_bc_residuals([1.0, 1.0, 0.0], [0.0, 0.0, 0.0], nn)[0] == pytest.approx(np.sqrt(2.0))


# =====================================================================================================================
# C02 — Couette–Poiseuille flow (8.4)–(8.5), flow rate, shear stress, backflow
# =====================================================================================================================
def test_channel_flow_V1_navier_stokes_residual_and_walls():  # V1 (C02): the (4.39b) residual of (8.5) is round-off only
    h, U, dpdx, mu, rho = 0.01, 0.1, 4.0, 1e-3, 1000.0
    uf = lambda X, T: np.stack([np.asarray(ch08.channel_flow(X[1], h, U, dpdx, mu)), 0.0 * X[1]])  # noqa: E731
    pf = lambda X, T: 1e5 + dpdx * X[0]  # noqa: E731  p = p0 + (dp/dx)x
    X = np.stack([RNG.uniform(0, 0.1, 40), RNG.uniform(0, h, 40)])
    step = h / 10
    terms = NS.ns_incompressible_terms(uf, pf, X, 0.0, rho, mu, (0.0, 0.0), h=step)
    bound = 10 * 4 * EPS * (U + dpdx * h ** 2 / (8 * mu)) / step ** 2 * mu / rho + 10 * EPS * 1e5 / (rho * step)
    assert np.max(np.abs(terms.residual)) < bound  # the stencil is exact for a quadratic: only round-off remains
    assert np.max(np.abs(terms.advective)) < bound  # u·∇u = 0 exactly (D02)
    assert maxrel(terms.pressure[0], -dpdx / rho) < 10 * EPS * 1e5 / (step * abs(dpdx))  # round-off of p0 = 1e5 Pa
    assert maxrel(terms.viscous[0], dpdx / rho) < 1e-5
    assert ch08.channel_flow(0.0, h, U, dpdx, mu) == 0.0 and ch08.channel_flow(h, h, U, dpdx, mu) == pytest.approx(U, rel=1e-15)
    with pytest.raises(ValueError):
        ch08.channel_flow(0.5 * h, h, U, dpdx=1.0, G=1.0)  # both conventions at once


def test_channel_flow_V1_parity_with_ch04_and_the_sign_convention():  # V1 (C02, R07, R08): G = −dp/dx
    h, U, mu = 0.02, 0.05, 2e-3
    y = np.linspace(0, h, 101)
    for dpdx in (-3.0, 0.0, 2.5):
        X = np.stack([np.zeros_like(y), y])
        u4 = NS.exact_solution("couette", X, U=U, h=h, G=-dpdx, mu=mu)[0][0]
        assert maxrel(ch08.channel_flow(y, h, U, dpdx, mu), u4) < 1e-14
        assert maxrel(ch08.channel_flow(y, h, U, mu=mu, G=-dpdx), u4) < 1e-14
    assert maxrel(ch08.channel_flow(y, h, 0.0, -3.0, mu), ch04.plane_poiseuille(y, G=3.0, h=h, mu=mu)) < 1e-14
    # wrong variant: passing ch04's G as dp/dx reverses the parabola
    assert maxrel(ch08.channel_flow(y, h, 0.0, 3.0, mu), ch04.plane_poiseuille(y, G=3.0, h=h, mu=mu)) > 1.0


def test_channel_flow_V2_sympy_residual_and_derivation_engine():  # V2 (C02, N04–N08): parallel_flow_sympy("channel")
    d = ch08.parallel_flow_sympy("channel")
    assert d["residual"] == 0 and d["checks"]["all_zero"]
    assert all(z0(v) for v in d["residuals"].values())
    x, y, h, U, mu, rho, G = d["symbols"]
    assert z0(d["profile"] - (U / h * y - G / (2 * mu) * y * (h - y)))  # (8.5)
    assert z0(mu * sp.diff(d["profile"], y, 2) - G)
    # D03: a non-constant pressure gradient leaves an x-dependence the y-side cannot match
    f = sp.Function("f")
    assert d["nonconstant_requires"] != 0 and z0(d["nonconstant_requires"].subs(f(x), 7 * x).doit() + 7 / rho)


def test_couette_poiseuille_V2_derivation():  # V2 derivation D04 ★★: integrate μu″ = dp/dx twice, no slip at both walls
    y, h, U, mu = sp.symbols("y h U mu", positive=True)
    G, c1, c2 = sp.symbols("G c1 c2", real=True)  # G = dp/dx (the book's constant)
    du = (G * y + c1) / mu  # step 1: μ du/dy = (dp/dx) y + c1
    u = sp.integrate(du, y) + c2 / mu  # step 2: μu = (dp/dx) y²/2 + c1 y + c2
    assert z0(mu * u - (G * y ** 2 / 2 + c1 * y + c2))
    A, B = -c1, -c2  # step 3: the book's A, B are minus our constants
    assert z0((-y ** 2 / 2 * G + mu * u + A * y + B))
    c2v = sp.solve(u.subs(y, 0), c2)[0]  # step 4
    assert c2v == 0
    c1v = sp.solve(sp.Eq(mu * u.subs({c2: 0, y: h}), mu * U), c1)[0]  # steps 5–6
    assert z0(c1v - (mu * U / h - h / 2 * G))
    prof = sp.simplify(u.subs({c1: c1v, c2: 0}))  # step 7
    assert z0(prof - (U / h * y + G / (2 * mu) * (y ** 2 - h * y)))  # step 8
    assert z0(prof - (U / h * y - G / (2 * mu) * y * (h - y)))  # step 9: (8.5)
    # check line: h = 1 cm, U = 0.1 m/s, dp/dx = 4 Pa/m, μ = 1e-3 at y = h/4 → −0.0125 m/s
    val = float(prof.subs({y: 0.0025, h: 0.01, U: 0.1, G: 4.0, mu: 1e-3}))
    assert val == pytest.approx(-0.0125, rel=1e-12)
    assert float(ch08.channel_flow(0.0025, 0.01, 0.1, 4.0, 1e-3)) == pytest.approx(-0.0125, rel=1e-12)


def test_channel_flow_rate_V1_quadrature_and_printed_V_units():  # V1 (N10) + V2 dimensions: Q = ∫u dy, V = Q/h (R9)
    for h, U, dpdx, mu in ((0.01, 0.1, -5.0, 1e-3), (0.003, 0.0, 2.0, 0.05), (0.02, -0.3, 7.0, 1e-3)):
        Q, V = ch08.channel_flow_rate(h, U, dpdx, mu)
        Qq = quad(lambda s: float(ch08.channel_flow(s, h, U, dpdx, mu)), 0, h, epsabs=0, epsrel=1e-13)[0]
        assert Q == pytest.approx(Qq, rel=1e-12) and V == pytest.approx(Qq / h, rel=1e-12)
        if U != 0:
            assert Q == pytest.approx(U * h / 2 * (1 - h ** 2 / (6 * mu * U) * dpdx), rel=1e-12)  # the book's form
    # units: Q [m²/s] vs V [m/s] — the printed middle form "V = ∫u dy" has the dimensions of Q (wrong variant)
    u, y, hq = Q_(1.0, "m/s"), Q_(1.0, "m"), Q_(1.0, "m")
    integral = u * y
    assert integral.dimensionality != u.dimensionality
    assert (integral / hq).dimensionality == u.dimensionality


def test_channel_shear_stress_V1_derivative_and_wall_values():  # V1 (C02, D05, R08): τ = μ du/dy, linear in y
    h, U, dpdx, mu = 0.01, 0.1, -3.0, 1e-3
    y = np.linspace(0.0, h, 51)
    dy = 1e-7
    fd = mu * (np.asarray(ch08.channel_flow(y + dy, h, U, dpdx, mu)) - np.asarray(ch08.channel_flow(y - dy, h, U, dpdx, mu))) / (2 * dy)
    assert maxrel(ch08.channel_shear_stress(y, h, U, dpdx, mu), fd) < 1e-7
    assert np.allclose(np.diff(ch08.channel_shear_stress(y, h, U, dpdx, mu), 2), 0.0, atol=1e-15)  # linear
    # plane Poiseuille: |τ_w| = (h/2)|dp/dx| at both walls, equal and opposite
    t0, th = ch08.channel_shear_stress(0.0, h, 0.0, dpdx, mu), ch08.channel_shear_stress(h, h, 0.0, dpdx, mu)
    assert abs(t0) == pytest.approx(h / 2 * abs(dpdx)) and th == pytest.approx(-t0)
    assert ch08.channel_shear_stress(0.3 * h, h, U, 0.0, mu) == pytest.approx(mu * U / h)  # Couette: uniform


def test_channel_backflow_V1_threshold_by_bisection_and_state():  # V1 (N09, D05): backflow iff dp/dx > 2μU/h²
    h, U, mu = 0.01, 0.1, 1e-3
    thr = ch08.channel_backflow_threshold(U, h, mu)
    assert thr == pytest.approx(2.0, rel=1e-14)  # D05 numbers: 2 Pa/m
    yy = np.linspace(0, h, 4001)[1:-1]
    fmin = lambda g: float(np.min(ch08.channel_flow(yy, h, U, g, mu))) / U + 1e-9  # noqa: E731
    g_star = brentq(fmin, 0.5, 10.0, xtol=1e-12)  # first gradient at which u < 0 somewhere inside
    assert g_star == pytest.approx(thr, rel=2e-3)  # limited by the y grid (reversal starts at y → 0)
    for g, back in ((1.9, False), (2.1, True), (-5.0, False), (6.5, True)):
        st = ch08.couette_poiseuille_state(h, U, g, mu)
        assert st["backflow"] is back
    st = ch08.couette_poiseuille_state(h, U, 4.0, mu)
    assert st["y_reversal"] == pytest.approx(0.005, rel=1e-12)  # D05: 5 mm at 4 Pa/m
    assert float(ch08.channel_flow(st["y_reversal"], h, U, 4.0, mu)) == pytest.approx(0.0, abs=1e-15)
    assert st["zero_flow_dpdx"] == pytest.approx(6 * mu * U / h ** 2)
    assert ch08.channel_flow_rate(h, U, st["zero_flow_dpdx"], mu)[0] == pytest.approx(0.0, abs=1e-16)
    assert st["Q"] == pytest.approx(st["Q_couette"] + st["Q_poiseuille"], rel=1e-14)
    assert st["tau_bottom"] == pytest.approx(float(ch08.channel_shear_stress(0.0, h, U, 4.0, mu)))
    # Poiseuille u_max = 1.5 V (D05)
    sp_ = ch08.couette_poiseuille_state(h, 0.0, -3.0, mu)
    assert sp_["u_max"] == pytest.approx(1.5 * sp_["V"], rel=1e-12) and sp_["y_umax"] == pytest.approx(h / 2)


def test_channel_flow_V7_superposition_and_mirror_symmetry():  # V7 (C02): linearity; Poiseuille symmetric about h/2
    h, mu = 0.01, 1e-3
    y = np.linspace(0, h, 41)
    both = ch08.channel_flow(y, h, 0.2, -3.0, mu)
    assert maxrel(both, np.asarray(ch08.channel_flow(y, h, 0.2, 0.0, mu)) + np.asarray(ch08.channel_flow(y, h, 0.0, -3.0, mu))) < 1e-14
    pois = np.asarray(ch08.channel_flow(y, h, 0.0, -3.0, mu))
    assert maxrel(pois, pois[::-1]) < 1e-12
    # reversing everything reverses u
    assert maxrel(ch08.channel_flow(y, h, -0.2, 3.0, mu), -np.asarray(both)) < 1e-15


# =====================================================================================================================
# C03 — pipe Poiseuille flow (8.6)–(8.8), Hagen–Poiseuille, f = 64/Re
# =====================================================================================================================
def test_pipe_poiseuille_V1_parity_ns_residual_and_ch03():  # V1 (C03): ch04 exact_solution, ch03 pipe_profile, (4.39b)
    a, dpdz, mu, rho = 0.005, -40.0, 1e-3, 1000.0
    R = np.linspace(0, a, 81)
    X = np.stack([R, np.zeros_like(R), np.zeros_like(R)])
    w4 = NS.exact_solution("pipe_poiseuille", X, G=-dpdz, R=a, mu=mu, g=0.0)[0][2]
    assert maxrel(ch08.pipe_poiseuille(R, a, dpdz, mu), w4) < 1e-14
    assert maxrel(ch08.pipe_poiseuille(R, a, mu=mu, G=-dpdz), w4) < 1e-14
    V = ch08.pipe_flow_rate(a, dpdz, mu)[1]
    assert maxrel(ch08.pipe_poiseuille(R, a, dpdz, mu), ch03.pipe_profile(R, z=1e3, U_mean=V, R=a, n0=20.0, L_e=1.0)) < 1e-12
    uf = lambda P, T: np.stack([0 * P[0], 0 * P[0], np.asarray(ch08.pipe_poiseuille(np.hypot(P[0], P[1]), a, dpdz, mu))])  # noqa: E731
    pf = lambda P, T: 1e5 + dpdz * P[2]  # noqa: E731
    th = RNG.uniform(0, 2 * np.pi, 30)
    rr = RNG.uniform(0, a, 30)
    P = np.stack([rr * np.cos(th), rr * np.sin(th), RNG.uniform(0, 1, 30)])
    step = a / 10
    res = NS.ns_incompressible_terms(uf, pf, P, 0.0, rho, mu, (0.0, 0.0, 0.0), h=step).residual
    assert np.max(np.abs(res)) < 100 * EPS * (V * 2 * mu / step ** 2 / rho + 1e5 / (rho * step))
    assert ch08.pipe_poiseuille(a, a, dpdz, mu) == 0.0


def test_pipe_poiseuille_V2_derivation_engine_and_D06():  # V2 derivation D06 ★★: (1/R)(Ru′)′ = G/μ, bounded, u(a) = 0
    d = ch08.parallel_flow_sympy("pipe")
    assert d["residual"] == 0 and d["checks"]["all_zero"]
    assert d["ln_term_unbounded"] == -sp.oo  # wrong variant: keeping A ln R is unbounded on the axis
    R, a, mu = sp.symbols("R a mu", positive=True)
    G, A, B = sp.symbols("G A B", real=True)
    w = sp.Function("w")
    lap = sp.diff(R * sp.diff(w(R), R), R) / R  # step 4: axisymmetric Laplacian
    first = sp.integrate(R * G / mu, R) + A  # step 6: R u′ = R²G/(2μ) + A
    assert z0(first - (R ** 2 / (2 * mu) * G + A))
    gen = sp.integrate(first / R, R) + B  # step 7
    assert z0(gen - (R ** 2 / (4 * mu) * G + A * sp.log(R) + B))
    assert z0(lap.subs(w(R), gen).doit() - G / mu)
    assert sp.limit(gen.subs({A: 1, B: 0, G: 0}), R, 0, "+") == -sp.oo  # step 8: A must vanish
    Bv = sp.solve(gen.subs({A: 0, R: a}), B)[0]  # step 9
    assert z0(Bv + a ** 2 / (4 * mu) * G)
    prof = sp.factor(gen.subs({A: 0, B: Bv}))  # step 10: (8.6)
    assert z0(prof - (R ** 2 - a ** 2) / (4 * mu) * G)
    # check line: a = 1 mm, dp/dz = −1000 Pa/m → u(0) = 0.25 m/s
    assert float(prof.subs({R: 0, a: 1e-3, G: -1000, mu: 1e-3})) == pytest.approx(0.25)
    # wrong variant 1/(2μ): its Laplacian gives 2G/μ, not G/μ
    assert not z0(lap.subs(w(R), (R ** 2 - a ** 2) / (2 * mu) * G).doit() - G / mu)


def test_pipe_flow_rate_V1_quadrature_umax_and_friction_factor():  # V1 (N14–N16): Q by quad, u_max = 2V, f·Re = 64
    for a, dpdz, mu, rho in ((1e-3, -1000.0, 1e-3, 1000.0), (0.02, -3.0, 1.8e-5, 1.2)):
        Q, V, umax = ch08.pipe_flow_rate(a, dpdz, mu)
        Qq = quad(lambda s: float(ch08.pipe_poiseuille(s, a, dpdz, mu)) * 2 * np.pi * s, 0, a, epsabs=0, epsrel=1e-13)[0]
        assert Q == pytest.approx(Qq, rel=1e-12) and V == pytest.approx(Q / (np.pi * a ** 2), rel=1e-14)
        assert umax == pytest.approx(float(ch08.pipe_poiseuille(0.0, a, dpdz, mu)), rel=1e-14) and umax == pytest.approx(2 * V)
        tau0 = float(ch08.pipe_wall_stress(a, dpdz))
        f = 8 * abs(tau0) / (rho * V ** 2)
        Re = V * 2 * a * rho / mu
        assert f == pytest.approx(64.0 / Re, rel=1e-12) and ch08.pipe_friction_factor(Re) == pytest.approx(f, rel=1e-12)
    # D07 check line: a = 1 mm, dp/dz = −1000: V = 0.125, Q = 3.93e-7, τ0 = −0.5 Pa, Re = 250, f = 0.256
    Q, V, _ = ch08.pipe_flow_rate(1e-3, -1000.0, 1e-3)
    assert (V, Q) == (pytest.approx(0.125), pytest.approx(3.927e-7, rel=1e-3))
    assert float(ch08.pipe_wall_stress(1e-3, -1000.0)) == pytest.approx(-0.5)
    assert ch08.pipe_friction_factor(250.0) == pytest.approx(0.256)


def test_pipe_wall_stress_V4_control_volume_force_balance():  # V4 (N15, D07 steps 4–5): πa²Δp + 2πaLτ0 = 0
    a, dpdz, mu, L, rho = 0.004, -25.0, 1e-3, 0.7, 1000.0
    p0, pL = 2e5, 2e5 + dpdz * L
    tau0 = float(ch08.pipe_wall_stress(a, dpdz))
    assert np.pi * a ** 2 * (p0 - pL) + 2 * np.pi * a * L * tau0 == pytest.approx(0.0, abs=1e-15)
    # momentum flux in = out (fully developed): ∫ρw² dA identical at both ends, so the force balance is exact
    mflux = quad(lambda s: rho * float(ch08.pipe_poiseuille(s, a, dpdz, mu)) ** 2 * 2 * np.pi * s, 0, a, epsrel=1e-13)[0]
    assert mflux > 0
    # (8.7) at R = a is (8.8), and τ = μ du/dR by finite differences
    Rr = np.linspace(0, a, 21)
    dR = 1e-8
    fd = mu * (np.asarray(ch08.pipe_poiseuille(Rr + dR, a, dpdz, mu)) - np.asarray(ch08.pipe_poiseuille(Rr - dR, a, dpdz, mu))) / (2 * dR)
    assert maxrel(ch08.pipe_shear_stress(Rr, dpdz), fd) < 1e-7
    assert float(ch08.pipe_shear_stress(a, dpdz)) == pytest.approx(tau0, rel=1e-15)
    # the cylindrical strain rate (R09) of (8.6) gives the same τ_zR
    Rs, ph, zs = CL.coordinates("cylindrical")
    Gs, mus, as_ = sp.symbols("G mu a", real=True)
    S = CL.strain_rate([0, 0, (Rs ** 2 - as_ ** 2) / (4 * mus) * Gs], "cylindrical")
    assert z0(2 * mus * S[0, 2] - Rs / 2 * Gs)


def test_pipe_V2_derivation_D07_flux_friction_factor():  # V2 derivation D07 ★★: τ, τ0, Q, V, u_max, f = 64/Re
    R, a, mu, rho, Lp = sp.symbols("R a mu rho L", positive=True)
    G = sp.Symbol("G", negative=True)
    u = (R ** 2 - a ** 2) / (4 * mu) * G
    tau = mu * sp.diff(u, R)
    assert z0(tau - R / 2 * G)  # step 2 (8.7)
    tau0 = tau.subs(R, a)
    assert z0(tau0 - a / 2 * G)  # step 3 (8.8)
    assert z0(sp.pi * a ** 2 * (-G * Lp) + 2 * sp.pi * a * Lp * tau0)  # steps 4–5: p(0) − p(L) = −G L
    assert z0(sp.integrate((R ** 2 - a ** 2) * R, (R, 0, a)) + a ** 4 / 4)  # step 7
    Q = sp.integrate(u * 2 * sp.pi * R, (R, 0, a))
    assert z0(Q + sp.pi * a ** 4 / (8 * mu) * G)
    V = Q / (sp.pi * a ** 2)
    assert z0(V + a ** 2 / (8 * mu) * G) and z0(u.subs(R, 0) - 2 * V)  # step 8
    assert z0(-tau0 - 4 * mu * V / a)  # step 9
    f = 8 * (-tau0) / (rho * V ** 2)
    Re = rho * V * 2 * a / mu
    assert z0(f - 64 / Re)  # step 10


# =====================================================================================================================
# C04 — circular Couette flow (8.9)–(8.12), pressure, stress, torque, power
# =====================================================================================================================
def test_circular_couette_V1_walls_ode_vorticity():  # V1 (C04): no slip, (1/R)(Ru)′ = 2A, Euler ODE residual
    for R1, R2, O1, O2 in ((0.01, 0.02, 1.0, 0.0), (0.05, 0.06, -2.0, 3.0), (0.3, 1.0, 0.5, 0.5)):
        u, A, B = ch08.circular_couette(np.array([R1, R2]), R1, R2, O1, O2, return_coeffs=True)
        assert maxrel(u, [O1 * R1, O2 * R2]) < 1e-13
        R = np.linspace(R1, R2, 41)
        dR = 1e-4 * R2  # second differences: round-off ~ eps·u/dR², truncation ~ dR²·u⁗
        f = lambda r: np.asarray(ch08.circular_couette(r, R1, R2, O1, O2))  # noqa: E731
        vort = ((R + dR) * f(R + dR) - (R - dR) * f(R - dR)) / (2 * dR) / R  # ω_z = (1/R) d(Ru)/dR
        assert np.max(np.abs(vort - 2 * A)) < 1e-6 * max(abs(A), abs(B) / R1 ** 2)
        upp = (f(R + dR) - 2 * f(R) + f(R - dR)) / dR ** 2
        up = (f(R + dR) - f(R - dR)) / (2 * dR)
        ode = upp + up / R - f(R) / R ** 2  # φ-momentum: u″ + u′/R − u/R² = 0
        assert np.max(np.abs(ode)) < 1e-5 * (abs(A) / R1 + abs(B) / R1 ** 3)
    # Ω₁ = Ω₂: rigid rotation, B = 0 (D08 check)
    _, A, B = ch08.circular_couette(0.5, 0.3, 1.0, 0.7, 0.7, return_coeffs=True)
    assert A == pytest.approx(0.7) and B == pytest.approx(0.0, abs=1e-16)
    # D08 numbers: R1 = 1 cm, R2 = 2 cm, Ω1 = 1, Ω2 = 0 → A = −1/3, B = 1.333e-4, u(1.5 cm) = 3.9 mm/s
    u, A, B = ch08.circular_couette(0.015, 0.01, 0.02, 1.0, 0.0, return_coeffs=True)
    assert (A, B) == (pytest.approx(-1 / 3), pytest.approx(4e-4 / 3)) and u == pytest.approx(3.889e-3, rel=1e-3)


def test_circular_couette_V2_engine_and_D08_derivation():  # V2 derivation D08 ★★ (+ D09 ★ limits)
    d = ch08.parallel_flow_sympy("circular_couette")
    assert d["residual"] == 0 and d["checks"]["all_zero"]
    R, R1, R2 = sp.symbols("R R_1 R_2", positive=True)
    O1, O2, A, B = sp.symbols("Omega_1 Omega_2 A B", real=True)
    v = sp.Function("v")
    first = sp.Eq(sp.diff(R * v(R), R) / R, 2 * A)  # step 5
    sol = sp.dsolve(sp.Eq(sp.diff(R * v(R), R), 2 * A * R), v(R)).rhs  # step 6–7
    C1 = [s for s in sol.free_symbols if s.name == "C1"][0]
    assert z0(sol.subs(C1, B) / 1 - (A * R + B / R)) or z0(sol - (A * R + C1 / R))
    u = A * R + B / R
    assert z0(first.lhs.subs(v(R), u).doit() - 2 * A)
    # the φ-component of the cylindrical vector Laplacian (P186) annihilates u = AR + B/R
    Rc, ph, zc = CL.coordinates("cylindrical")
    assert z0(CL.vector_laplacian([0, A * Rc + B / Rc, 0], "cylindrical")[1])
    cons = sp.solve([sp.Eq(O1 * R1 ** 2, A * R1 ** 2 + B), sp.Eq(O2 * R2 ** 2, A * R2 ** 2 + B)], [A, B], dict=True)[0]  # 8
    assert z0(cons[A] - (O2 * R2 ** 2 - O1 * R1 ** 2) / (R2 ** 2 - R1 ** 2))  # step 9
    assert z0(cons[B] + (O2 - O1) * R1 ** 2 * R2 ** 2 / (R2 ** 2 - R1 ** 2))  # step 10
    prof = u.subs(cons)
    book = ((O2 * R2 ** 2 - O1 * R1 ** 2) * R - (O2 - O1) * R1 ** 2 * R2 ** 2 / R) / (R2 ** 2 - R1 ** 2)  # (8.10)
    assert z0(prof - book)
    # D09: divide by R2², then R2 → ∞ with Ω2 = 0 gives (8.11); R1, Ω1 → 0 gives (8.12)
    assert z0(sp.limit(prof.subs(O2, 0), R2, sp.oo) - O1 * R1 ** 2 / R)
    assert z0(sp.limit(prof.subs(O1, 0), R1, 0) - O2 * R)
    assert z0(2 * sp.pi * R * (O1 * R1 ** 2 / R) - 2 * sp.pi * O1 * R1 ** 2)  # Γ independent of R
    # Wikipedia Taylor–Couette form (V1 cross-check): A = Ω1(μ − η²)/(1 − η²), B = Ω1R1²(1 − μ)/(1 − η²)
    mu_r, eta = O2 / O1, R1 / R2
    assert z0(cons[A] - O1 * (mu_r - eta ** 2) / (1 - eta ** 2))
    assert z0(cons[B] - O1 * R1 ** 2 * (1 - mu_r) / (1 - eta ** 2))


def test_circular_couette_V7_limits_and_ch05_ch03_parity():  # V7 (N21, R12): R2 → ∞, R1 → 0; ch05 ω = 2Ω1 trap
    R1, O1 = 0.02, 3.0
    R = np.linspace(R1, 0.3, 50)
    u_inf = ch08.circular_couette(R, R1, np.inf, O1, 0.0)
    assert maxrel(u_inf, O1 * R1 ** 2 / R) < 1e-15
    assert maxrel(ch08.circular_couette(R, R1, 1e6 * R1, O1, 0.0), u_inf) < 1e-10
    uc5 = ch05.rotating_cylinder_flow(R, R1, 2 * O1)[0]  # ch05 takes the cylinder's vorticity ω = 2Ω1
    assert maxrel(u_inf, uc5) < 1e-14
    assert maxrel(u_inf, ch05.rotating_cylinder_flow(R, R1, O1)[0]) > 0.4  # wrong variant: ω = Ω1 is off by 2
    assert maxrel(u_inf, VX.line_vortex(R, O1 * R1 ** 2)) < 1e-14  # Γ = 2πΩ1R1², B = Γ/2π
    Rr = np.linspace(0, 1.0, 30)
    assert maxrel(ch08.circular_couette(Rr, 0.0, 1.0, 0.0, 0.8), VX.solid_body_rotation(Rr, 0.8)) < 1e-15
    assert maxrel(ch08.circular_couette(Rr[1:], 1e-9, 1.0, 0.0, 0.8), VX.solid_body_rotation(Rr[1:], 0.8)) < 1e-12


def test_circular_couette_V1_navier_stokes_residual_with_pressure():  # V1 (C04, N20): Cartesian (4.39b) residual → 0
    R1, R2, O1, O2, rho, mu = 0.02, 0.05, 2.0, -1.0, 1000.0, 1e-3

    def uf(X, T):
        r = np.hypot(X[0], X[1])
        ut = np.asarray(ch08.circular_couette(r, R1, R2, O1, O2))
        return np.stack([-ut * X[1] / r, ut * X[0] / r])

    pf = lambda X, T: np.asarray(ch08.circular_couette_pressure(np.hypot(X[0], X[1]), R1, R2, O1, O2, rho, 1e5))  # noqa: E731
    th = RNG.uniform(0, 2 * np.pi, 25)
    rr = RNG.uniform(0.025, 0.045, 25)
    P = np.stack([rr * np.cos(th), rr * np.sin(th)])
    errs, steps = [], [4e-4, 2e-4, 1e-4]
    for s in steps:
        t = NS.ns_incompressible_terms(uf, pf, P, 0.0, rho, mu, (0.0, 0.0), h=s)
        errs.append(np.max(np.abs(t.residual)) / np.max(np.abs(t.advective)))
    assert errs[-1] < 1e-4  # ~ (h/R)² at h = 1e-4 m, R ≈ 0.03 m
    assert abs(observed_order(steps, errs) - 2.0) < ORDER_TOL  # the residual is pure O(h²) truncation → 0
    # the advective term is the centripetal −u²/R e_R, balanced by the pressure (N22)
    t = NS.ns_incompressible_terms(uf, pf, P, 0.0, rho, mu, (0.0, 0.0), h=1e-4)
    ut = np.asarray(ch08.circular_couette(rr, R1, R2, O1, O2))
    radial = (t.advective[0] * np.cos(th) + t.advective[1] * np.sin(th))
    assert maxrel(radial, -ut ** 2 / rr) < 1e-4  # O((h/R)²) stencil truncation at h = 1e-4 m
    for case, exp in (("channel", 0), ("pipe", 0)):
        assert ch08.advective_acceleration_check(case)["residual"] == 0
    assert ch08.advective_acceleration_check("circular_couette")["residual"] == 0


def test_circular_couette_pressure_V1_radial_balance():  # V1 (N20): dp/dR = ρu²/R by finite differences, p(R1) = p1
    R1, R2, O1, O2, rho = 0.1, 0.25, 4.0, 1.0, 998.0
    R = np.linspace(R1 + 1e-3, R2 - 1e-3, 30)
    dR = 1e-6
    f = lambda r: np.asarray(ch08.circular_couette_pressure(r, R1, R2, O1, O2, rho, 7.0))  # noqa: E731
    fd = (f(R + dR) - f(R - dR)) / (2 * dR)
    assert maxrel(fd, rho * np.asarray(ch08.circular_couette(R, R1, R2, O1, O2)) ** 2 / R) < 1e-7
    assert f(np.array([R1]))[0] == pytest.approx(7.0)
    # solid-body limit: p = p1 + ρΩ²R²/2 (paraboloid)
    assert ch08.circular_couette_pressure(0.5, 0.0, 1.0, 0.0, 2.0, 1000.0, 0.0) == pytest.approx(1000 * 4 * 0.25 / 2)
    st = ch08.circular_couette_state(R1, R2, O1, O2, rho=rho)
    assert st["dp_gap"] == pytest.approx(float(f(np.array([R2]))[0]) - 7.0, rel=1e-12)


def test_circular_couette_V4_torque_and_power_equal_dissipation():  # V4 (R10, R11): torque constant; power in = ∫ε dA
    mu = 1e-3
    for R1, R2, O1, O2 in ((0.01, 0.02, 1.0, 0.0), (0.05, 0.08, -1.0, 2.5), (0.02, np.inf, 3.0, 0.0)):
        pw = ch08.circular_couette_power(R1, R2, O1, O2, mu)
        Rs = np.array([R1, 1.3 * R1, 1.7 * R1]) if np.isinf(R2) else np.linspace(R1, R2, 7)
        torque = 2 * np.pi * Rs ** 2 * np.asarray(ch08.circular_couette_shear_stress(Rs, R1, R2, O1, O2, mu))
        assert maxrel(torque, np.full_like(Rs, torque[0])) < 1e-13  # angular momentum balance: torque independent of R
        assert pw["torque_inner"] + pw["torque_outer"] == pytest.approx(0.0, abs=1e-18)
        assert pw["power_in"] == pytest.approx(pw["dissipation"], rel=1e-12)
        assert pw["dissipation_quad"] == pytest.approx(pw["dissipation"], rel=1e-10)
        assert pw["power_in"] >= 0
        # σ_Rφ = μR d(u/R)/dR by finite differences
        R = Rs[1]
        dR = 1e-7 * R
        g = lambda r: float(ch08.circular_couette(r, R1, R2, O1, O2)) / r  # noqa: E731
        assert float(ch08.circular_couette_shear_stress(R, R1, R2, O1, O2, mu)) == pytest.approx(
            mu * R * (g(R + dR) - g(R - dR)) / (2 * dR), rel=1e-6)
    # R2 = ∞, Ω2 = 0: power = 4πμΩ1²R1² = ch05's dissipation outside a cylinder (Γ = 2πΩ1R1²)
    R1, O1 = 0.02, 3.0
    pw = ch08.circular_couette_power(R1, np.inf, O1, 0.0, mu)
    assert pw["power_in"] == pytest.approx(4 * np.pi * mu * O1 ** 2 * R1 ** 2, rel=1e-13)
    Rout = 1.0
    c5 = ch05.dissipation_outside_cylinder(R1, Rout, 2 * np.pi * O1 * R1 ** 2, mu)
    assert c5["dissipated"] == pytest.approx(pw["dissipation"] * (1 - R1 ** 2 / Rout ** 2), rel=1e-10)
    assert c5["power_in"] == pytest.approx(pw["power_in"], rel=1e-13)
    assert float(ch08.circular_couette_shear_stress(R1, R1, np.inf, O1, 0.0, mu)) == pytest.approx(-2 * mu * O1, rel=1e-14)
    # wrong variant (R15): the printed (2πR1)σ_Rφu_φ is negative — the power into the fluid is minus it
    printed = 2 * np.pi * R1 * float(ch08.circular_couette_shear_stress(R1, R1, np.inf, O1, 0.0, mu)) * O1 * R1
    assert printed < 0 and -printed == pytest.approx(pw["power_in"], rel=1e-13)


def test_circular_couette_R10_V1_zero_net_viscous_force():  # V1 (R10): the viscous vortex exerts no net viscous force
    R1, O1 = 0.02, 3.0

    def uf(X, T):
        r = np.hypot(X[0], X[1])
        ut = np.asarray(ch08.circular_couette(r, R1, np.inf, O1, 0.0))
        return np.stack([-ut * X[1] / r, ut * X[0] / r])

    P = np.stack([np.array([0.03, 0.0, -0.05]), np.array([0.01, 0.04, 0.02])])
    scale = 1e-3 * O1 * R1 ** 2 / 0.03 ** 3  # μ × (typical |∇²u| of a non-vortex flow of this size)
    hs = [4e-4, 2e-4, 1e-4]
    mags = []
    for hh in hs:
        lap, div2S, mcurl = NS.viscous_force_forms(uf, P, 0.0, 1e-3, hh)
        mags.append(max(np.max(np.abs(f_)) for f_ in (lap, div2S, mcurl)))
    assert mags[-1] < 1e-3 * scale and abs(observed_order(hs, mags) - 2.0) < 0.25  # → 0: pure stencil truncation
    # contrast: the Lamb–Oseen vortex of the same circulation has an O(1) net viscous force at the same points
    sig = 0.02
    ulo = lambda X, T: (lambda rr: np.stack([-np.asarray(VX.gaussian_vortex(rr, 2 * np.pi * O1 * R1 ** 2, sig)[0]) * X[1] / rr,  # noqa: E731
                                             np.asarray(VX.gaussian_vortex(rr, 2 * np.pi * O1 * R1 ** 2, sig)[0]) * X[0] / rr]))(np.hypot(X[0], X[1]))
    lap_lo = NS.viscous_force_forms(ulo, P, 0.0, 1e-3, 1e-4)[0]
    assert np.max(np.abs(lap_lo)) > 0.1 * scale
    assert ch08.circular_couette(0.03, R1, np.inf, O1, 0.0) != 0  # but the flow is not at rest


def test_circular_couette_state_V1_rayleigh_criterion():  # V1 (B1/IF3): Rayleigh stable iff μ = Ω2/Ω1 > η² (co-rotation)
    R1, R2 = 0.5, 1.0
    for O2, stable in ((0.0, False), (0.2, False), (0.3, True), (1.0, True), (2.0, True)):
        st = ch08.circular_couette_state(R1, R2, 1.0, O2)
        assert st["rayleigh_stable"] is stable  # η² = 0.25
        assert st["vorticity"] == pytest.approx(2 * st["A"])
    st = ch08.circular_couette_state(0.01, 0.02, 1.0, 0.0)
    assert st["power_in"] == pytest.approx(st["dissipation"], rel=1e-12)


# =====================================================================================================================
# C05 — lubrication scaling and balance (8.13)–(8.17)
# =====================================================================================================================
def test_lubrication_nondim_V2_coefficient_sets_and_printed_8_13b():  # V2 (C05, N24, N26–N29): (8.15), (8.16a,b)
    d = ch08.lubrication_nondim_sympy()
    L, h, U, rho, mu, Pa, eps, Re, Lam = d["symbols"]
    assert d["matches"] == {"x": True, "y": True, "continuity": True}
    assert d["x_coeffs"] == {"inertia": eps ** 2 * Re, "pressure": 1 / Lam, "diff_along": eps ** 2, "diff_across": 1}
    assert d["y_coeffs"] == {"inertia": eps ** 4 * Re, "pressure": 1 / Lam, "diff_along": eps ** 4, "diff_across": eps ** 2}
    # wrong variant (R6): the printed ∂p/∂x in (8.13b) gives pressure coefficient ε/Λ — not (8.16b)
    assert z0(d["printed_13b_y_coeffs"]["pressure"] - eps / Lam)
    assert d["printed_13b_y_coeffs"] != d["y_coeffs"]
    assert ch08.lubrication_nondim_sympy(printed_8_13b=True)["matches"]["y"] is False
    assert d["profile_8_19_matches"] and z0(d["book_form_top_wall"] - (sp.Symbol("U_h", real=True) + sp.Symbol("U_0", real=True)))


def test_lubrication_scaling_V2_derivation():  # V2 derivation D10 ★★★: every line of Part F steps 1–14
    L, eps, U, rho, mu, Pa = sp.symbols("L epsilon U rho mu P_a", positive=True)
    x, y, t = sp.symbols("x y t", real=True)
    h = eps * L
    Re_L, Lam = rho * U * L / mu, mu * U * L / (Pa * h ** 2)
    # concrete O(1) starred fields make every coefficient a pure ratio (steps 1–2: the substitutions)
    xs, ys, ts = x / L, y / h, U * t / L
    us = sp.exp(xs + 2 * ys + 3 * ts)
    vs = sp.exp(2 * xs - ys + ts)
    ps = sp.exp(-xs + 5 * ys - 2 * ts)
    u, v, p = U * us, eps * U * vs, Pa * ps
    X, Y, T = sp.symbols("X Y T")
    usf, vsf, psf = (f.subs({x: L * X, y: h * Y, t: L * T / U}) for f in (us, vs, ps))
    dstar = lambda f, var, n=1: sp.diff(f, var, n).subs({X: xs, Y: ys, T: ts})  # noqa: E731  starred derivative
    ratio = lambda a, b: sp.simplify(sp.expand(a / b))  # noqa: E731
    # step 3–4: continuity (8.15): both terms carry U/L; dividing leaves no coefficient
    assert z0(ratio(sp.diff(u, x), dstar(usf, X)) - U / L) and z0(ratio(sp.diff(v, y), dstar(vsf, Y)) - U / L)
    # step 5: each x-inertia term is U²/L × its starred twin
    assert z0(ratio(sp.diff(u, t), dstar(usf, T)) - U ** 2 / L)
    assert z0(ratio(u * sp.diff(u, x), usf.subs({X: xs, Y: ys, T: ts}) * dstar(usf, X)) - U ** 2 / L)
    assert z0(ratio(v * sp.diff(u, y), vsf.subs({X: xs, Y: ys, T: ts}) * dstar(usf, Y)) - U ** 2 / L)
    # step 6: −(1/ρ)∂p/∂x = −(P_a/ρL)∂p*/∂x*
    assert z0(ratio(-sp.diff(p, x) / rho, -dstar(psf, X)) - Pa / (rho * L))
    # step 7: μ/ρ ∂²u/∂x² → μU/(ρL²), μ/ρ ∂²u/∂y² → μU/(ρε²L²)
    assert z0(ratio(mu / rho * sp.diff(u, x, 2), dstar(usf, X, 2)) - mu * U / (rho * L ** 2))
    assert z0(ratio(mu / rho * sp.diff(u, y, 2), dstar(usf, Y, 2)) - mu * U / (rho * eps ** 2 * L ** 2))
    # steps 8–10: × ρε²L²/(μU) → {ε²Re_L, 1/Λ, ε², 1} (8.16a)
    fx = rho * eps ** 2 * L ** 2 / (mu * U)
    assert z0(U ** 2 / L * fx - eps ** 2 * Re_L)
    assert z0(Pa / (rho * L) * fx - 1 / Lam)
    assert z0(mu * U / (rho * L ** 2) * fx - eps ** 2) and z0(mu * U / (rho * eps ** 2 * L ** 2) * fx - 1)
    # step 11: y-inertia = εU²/L × starred
    assert z0(ratio(sp.diff(v, t), dstar(vsf, T)) - eps * U ** 2 / L)
    assert z0(ratio(u * sp.diff(v, x), usf.subs({X: xs, Y: ys, T: ts}) * dstar(vsf, X)) - eps * U ** 2 / L)
    assert z0(ratio(v * sp.diff(v, y), vsf.subs({X: xs, Y: ys, T: ts}) * dstar(vsf, Y)) - eps * U ** 2 / L)
    # step 12: pressure P_a/(ρεL), friction μεU/(ρL²) and μU/(ρεL²)
    assert z0(ratio(-sp.diff(p, y) / rho, -dstar(psf, Y)) - Pa / (rho * eps * L))
    assert z0(ratio(mu / rho * sp.diff(v, x, 2), dstar(vsf, X, 2)) - mu * eps * U / (rho * L ** 2))
    assert z0(ratio(mu / rho * sp.diff(v, y, 2), dstar(vsf, Y, 2)) - mu * U / (rho * eps * L ** 2))
    # steps 13–14: × ρε³L²/(μU) → {ε⁴Re_L, 1/Λ, ε⁴, ε²} (8.16b)
    fy = rho * eps ** 3 * L ** 2 / (mu * U)
    assert z0(eps * U ** 2 / L * fy - eps ** 4 * Re_L) and z0(Pa / (rho * eps * L) * fy - 1 / Lam)
    assert z0(mu * eps * U / (rho * L ** 2) * fy - eps ** 4) and z0(mu * U / (rho * eps * L ** 2) * fy - eps ** 2)
    # the printed (8.13b) with ∂p/∂x: coefficient (P_a/ρL)·fy = ε/Λ ≠ 1/Λ (wrong variant)
    assert not z0(ratio(-sp.diff(p, x) / rho, -dstar(psf, X)) * fy - 1 / Lam)
    assert z0(ratio(-sp.diff(p, x) / rho, -dstar(psf, X)) * fy - eps / Lam)


def test_lubrication_balance_V2_derivation_and_units():  # V2 derivation D11 ★★ + dimensions of (8.17a) (R7)
    eps, Re, Lam = sp.symbols("epsilon Re_L Lambda", positive=True)
    px, uyy, uxx, inert = sp.symbols("p_x u_yy u_xx I", real=True)
    x16a = sp.Eq(eps ** 2 * Re * inert, -px / Lam + eps ** 2 * uxx + uyy)
    lim = sp.limit(x16a.rhs - x16a.lhs, eps, 0)  # steps 2–3: inertia and along-gap friction drop
    assert z0(lim - (-px / Lam + uyy))
    L, U, rho, mu, Pa = sp.symbols("L U rho mu P_a", positive=True)
    h = eps * L
    # step 5: × μU/(ρε²L²) with 1/Λ = P_a h²/(μUL): −(P_a/(ρL))p*_x + (μU/(ρh²))u*_yy = −(1/ρ)∂p/∂x + ν∂²u/∂y²
    back = sp.simplify((-px * Pa * h ** 2 / (mu * U * L) + uyy) * mu * U / (rho * eps ** 2 * L ** 2))
    assert z0(back - (-(Pa / L) * px / rho + mu / rho * U / h ** 2 * uyy))
    # units: (1/ρ)∂p/∂x and ν∂²u/∂y² are both m/s²; the printed form without ν is not (R7)
    rho_q, dpdx_q, nu_q, u_q, y_q = Q_(870, "kg/m**3"), Q_(1, "Pa/m"), Q_(1e-4, "m**2/s"), Q_(1, "m/s"), Q_(1, "m")
    a1, a2, printed = dpdx_q / rho_q, nu_q * u_q / y_q ** 2, u_q / y_q ** 2
    assert a1.dimensionality == a2.dimensionality == Q_(1, "m/s**2").dimensionality
    assert printed.dimensionality != a1.dimensionality
    # the engine-film numbers of the design: ε = 1e-3, Re_L = 4350, ε²Re_L = 4.35e-3, Λ(P_a) = 49.3
    s = ch08.lubrication_scales(0.05, 50e-6, 5.0, 870.0, 0.05)
    assert (s["eps"], s["Re_L"], s["eps2_Re_L"]) == (pytest.approx(1e-3), pytest.approx(4350.0), pytest.approx(4.35e-3))
    assert s["Lambda"] == pytest.approx(49.35, rel=1e-3) and s["p_visc"] == pytest.approx(5e6)


def test_lubrication_scales_V1_definitions_and_term_magnitudes():  # V1 (N26, N31) + V2 parity with the sympy coefficients
    L, h, U, rho, mu = 0.05, 50e-6, 5.0, 870.0, 0.05
    s = ch08.lubrication_scales(L, h, U, rho, mu)
    assert s["eps"] == pytest.approx(h / L) and s["Re_L"] == pytest.approx(rho * U * L / mu)
    assert s["Lambda"] == pytest.approx(mu * U * L / (ch08.P_ATM * h ** 2)) and s["p_visc"] == pytest.approx(mu * U * L / h ** 2)
    # wrong variant: Re with the gap instead of the length (a factor ε)
    assert s["Re_L"] != pytest.approx(rho * U * h / mu)
    d = ch08.lubrication_nondim_sympy()
    Ls, hs, Us, rhos, mus, Pas, eps, Re, Lam = d["symbols"]
    vals = {eps: s["eps"], Re: s["Re_L"], Lam: s["Lambda"]}
    m = ch08.lubrication_term_magnitudes(L, h, U, rho, mu, p_scale="atm")
    for k_num, grp, k_sym in (("x_inertia", "x_coeffs", "inertia"), ("x_pressure", "x_coeffs", "pressure"),
                              ("x_diff_along", "x_coeffs", "diff_along"), ("x_diff_across", "x_coeffs", "diff_across"),
                              ("y_inertia", "y_coeffs", "inertia"), ("y_pressure", "y_coeffs", "pressure"),
                              ("y_diff_along", "y_coeffs", "diff_along"), ("y_diff_across", "y_coeffs", "diff_across")):
        assert m[k_num] == pytest.approx(float(sp.sympify(d[grp][k_sym]).subs(vals)), rel=1e-12)
    mv = ch08.lubrication_term_magnitudes(L, h, U, rho, mu, p_scale="viscous")
    assert mv["x_pressure"] == pytest.approx(1.0) and mv["y_pressure"] == pytest.approx(1.0)
    md = ch08.lubrication_term_magnitudes(L, h, U, rho, mu, p_scale="dynamic")
    assert md["x_pressure"] == pytest.approx(s["eps2_Re_L"])
    with pytest.raises(ValueError):
        ch08.lubrication_term_magnitudes(L, h, U, rho, mu, p_scale="bogus")


# =====================================================================================================================
# C06 — lubrication profile (8.18)–(8.19), gap flux, 1-D Reynolds equation, Hele-Shaw
# =====================================================================================================================
def test_lubrication_velocity_V1_walls_flux_and_printed_form():  # V1 (C06, N33, N34): u(0) = U0, u(h) = U_h
    h, dpdx, Uh, U0, mu = 2e-4, -3e3, 1.5, 0.4, 0.02
    assert ch08.lubrication_velocity(0.0, h, dpdx, Uh, U0, mu) == pytest.approx(U0)
    assert ch08.lubrication_velocity(h, h, dpdx, Uh, U0, mu) == pytest.approx(Uh)
    # wrong variant (R8b): the printed form gives u(h) = U_h + U_0
    assert ch08.lubrication_velocity(h, h, dpdx, Uh, U0, mu, form="book") == pytest.approx(Uh + U0)
    assert ch08.lubrication_velocity(h, h, dpdx, Uh, 0.0, mu, form="book") == pytest.approx(Uh)  # U0 = 0 agrees
    q = quad(lambda s: float(ch08.lubrication_velocity(s, h, dpdx, Uh, U0, mu)), 0, h, epsabs=0, epsrel=1e-13)[0]
    assert float(ch08.lubrication_flux(h, dpdx, Uh, U0, mu)) == pytest.approx(q, rel=1e-12)
    # U0 = 0, constant gap: the channel profile (8.5)
    y = np.linspace(0, h, 31)
    assert maxrel(ch08.lubrication_velocity(y, h, dpdx, Uh, 0.0, mu), ch08.channel_flow(y, h, Uh, dpdx, mu)) < 1e-13
    # the profile satisfies (8.17a) μu″ = ∂p/∂x (quadratic → exact second difference)
    dy = h / 50
    upp = (np.asarray(ch08.lubrication_velocity(y + dy, h, dpdx, Uh, U0, mu)) - 2 * np.asarray(ch08.lubrication_velocity(y, h, dpdx, Uh, U0, mu))
           + np.asarray(ch08.lubrication_velocity(y - dy, h, dpdx, Uh, U0, mu))) / dy ** 2
    assert maxrel(mu * upp, np.full_like(y, dpdx)) < 1e-6
    with pytest.raises(ValueError):
        ch08.lubrication_velocity(0.0, h, dpdx, form="other")


def test_lubrication_profile_V2_derivation():  # V2 derivation D12 ★★: (8.18) → the consistent (8.19)
    y, h, mu = sp.symbols("y h mu", positive=True)
    px, Uh, U0, A, B = sp.symbols("p_x U_h U_0 A B", real=True)
    u18 = sp.integrate(sp.integrate(px / mu, y), y) + A * y + B  # steps 2–3
    assert z0(u18 - (px / mu * y ** 2 / 2 + A * y + B))
    Bv = sp.solve(u18.subs(y, 0) - U0, B)[0]  # step 4
    assert Bv == U0
    Av = sp.solve(u18.subs({B: U0, y: h}) - Uh, A)[0]  # step 5
    assert z0(Av - ((Uh - U0) / h - h / (2 * mu) * px))
    u = u18.subs({A: Av, B: Bv})
    assert z0(u - (px / (2 * mu) * (y ** 2 - h * y) + (Uh - U0) * y / h + U0))  # step 6
    cons = -h ** 2 / (2 * mu) * px * (y / h) * (1 - y / h) + Uh * y / h + U0 * (1 - y / h)  # step 7
    assert z0(u - cons)
    book = -h ** 2 / (2 * mu) * px * (y / h) * (1 - y / h) + Uh * y / h + U0  # step 8 (printed)
    assert z0(book.subs(U0, 0) - cons.subs(U0, 0)) and not z0(book.subs(y, h) - Uh)
    d = ch08.lubrication_nondim_sympy()
    assert d["profile_8_19_matches"]


def test_reynolds_equation_V2_derivation():  # V2 derivation D13 ★★★: Leibniz, kinematic walls, flux, Reynolds equation
    d = ch08.reynolds_equation_sympy()
    assert d["leibniz_residual"] == 0 and d["kinematic_cancellation"] == 0 and d["reynolds_residual"] == 0
    # independent rebuild with generic h(x, t), p(x, t)
    x, y, t = sp.symbols("x y t", real=True)
    mu = sp.Symbol("mu", positive=True)
    Uh, U0 = sp.symbols("U_h U_0", real=True)
    h, p = sp.Function("h")(x, t), sp.Function("p")(x, t)
    px = sp.diff(p, x)
    u = -h ** 2 / (2 * mu) * px * (y / h) * (1 - y / h) + Uh * y / h + U0 * (1 - y / h)
    q = sp.integrate(u, (y, 0, h))
    # steps 9–11: q = −h³p_x/(12μ) + (U0 + U_h)h/2 from ∫(y/h)(1 − y/h) = h/6 and ∫y/h = ∫(1 − y/h) = h/2
    hs = sp.Symbol("h", positive=True)
    ys = sp.Symbol("y", positive=True)
    assert z0(sp.integrate(ys / hs * (1 - ys / hs), (ys, 0, hs)) - hs / 6)
    assert z0(sp.integrate(ys / hs, (ys, 0, hs)) - hs / 2) and z0(sp.integrate(1 - ys / hs, (ys, 0, hs)) - hs / 2)
    assert z0(q - (-h ** 3 / (12 * mu) * px + (U0 + Uh) * h / 2))
    # step 2: Leibniz ∫₀ʰ ∂u/∂x dy = ∂q/∂x − u(h) ∂h/∂x
    lhs = sp.integrate(sp.diff(u, x), (y, 0, h))
    assert z0(lhs - (sp.diff(q, x) - u.subs(y, h) * sp.diff(h, x)))
    # steps 1, 3–7: v(0) = 0 and continuity give v(h) = −∫u_x dy; the kinematic wall v(h) = h_t + u(h)h_x ⇒ h_t + q_x = 0
    v_h = -lhs
    ht_from_kinematic = v_h - u.subs(y, h) * sp.diff(h, x)
    assert z0(ht_from_kinematic + sp.diff(q, x))
    # step 12: the pressure form ∂/∂x(h³p_x/12μ) = h_t + ((U0 + U_h)/2)h_x
    hts = sp.Symbol("h_t")
    reyn = sp.diff(h ** 3 / (12 * mu) * px, x) - (hts + (U0 + Uh) / 2 * sp.diff(h, x))
    assert z0(reyn.subs(hts, -sp.diff(q, x)))
    # check line: constant h, ∂p/∂x = 0 → Couette flux (U0 + U_h)h/2
    assert float(ch08.lubrication_flux(1e-4, 0.0, 2.0, 1.0)) == pytest.approx(1.5e-4)


def test_reynolds_pressure_1d_V1_exact_slider_and_V3_order_two():  # V1 + V3 (N32, I15): quad route exact, trapezoid order 2
    L, h0, al, U, mu = 0.05, 50e-6, 0.5, 5.0, 0.05
    pex = lambda x: np.asarray(ch08.slider_bearing(x, h0, al, L, U, mu))  # noqa: E731
    x = np.linspace(0, L, 41)
    p, q = ch08.reynolds_pressure_1d(x, lambda s: h0 * (1 + al * s / L), U_0=-U, U_h=0.0, mu=mu)
    assert maxrel(p, pex(x)) < 1e-10  # pad frame: floor −U, pad at rest
    assert q == pytest.approx(-(1 + al) / (2 + al) * U * h0, rel=1e-12)  # = Example 8.1's C1
    hs, es = [], []
    for N in (21, 41, 81, 161, 321):
        xx = np.linspace(0, L, N)
        pp, _ = ch08.reynolds_pressure_1d(xx, h0 * (1 + al * xx / L), U_0=-U, U_h=0.0, mu=mu)
        es.append(maxrel(pp, pex(xx)))
        hs.append(xx[1] - xx[0])
    assert abs(observed_order(hs, es) - 2.0) < ORDER_TOL
    # end pressures honoured; a uniform gap with no wall motion carries a pure Poiseuille flux
    xx = np.linspace(0, 1.0, 11)
    pp, qq = ch08.reynolds_pressure_1d(xx, np.full(11, 1e-3), 0.0, 0.0, 1e-3, p_left=12.0, p_right=0.0)
    assert pp[0] == pytest.approx(12.0) and pp[-1] == pytest.approx(0.0, abs=1e-12)
    assert qq == pytest.approx(-(1e-3) ** 3 / (12e-3) * (-12.0), rel=1e-12)


def test_hele_shaw_V1_velocity_potential_and_mean():  # V1 (N36): u = ∇φ, no slip at the plates, ū = −(h²/12μ)∇p
    h, mu = 1e-3, 1e-3
    gp = (np.array([-50.0, 20.0]), np.array([10.0, -30.0]))
    z = np.linspace(0, h, 21)
    for gx, gy in zip(*gp):
        u, v = ch08.hele_shaw_velocity(z, h, (gx, gy), mu)
        assert abs(u[0]) + abs(u[-1]) + abs(v[0]) + abs(v[-1]) == 0.0
        # φ = −z(h − z)p/2μ ⇒ ∂φ/∂x = −z(h − z)(∂p/∂x)/2μ = u
        dx = 1e-3
        dphidx = (np.asarray(ch08.hele_shaw_potential(gx * dx, z, h, mu)) - np.asarray(ch08.hele_shaw_potential(-gx * dx, z, h, mu))) / (2 * dx)
        assert maxrel(dphidx, u) < 1e-12
        ub = quad(lambda s: float(ch08.hele_shaw_velocity(s, h, (gx, gy), mu)[0]), 0, h, epsrel=1e-13)[0] / h
        assert ch08.hele_shaw_mean_velocity((gx, gy), h, mu)[0] == pytest.approx(ub, rel=1e-12)
    # wrong variant: φ = +z(h − z)p/2μ gives −u
    assert maxrel(-np.asarray(ch08.hele_shaw_potential(1.0, z, h, mu)), z * (h - z) / (2 * mu)) < 1e-15


def test_hele_shaw_V2_harmonic_potential():  # V2 (N36): ∇²p = 0 ⇒ ∇²φ = 0 in every plane; ideal-flow streamlines
    xs, ys = sp.symbols("x y", real=True)
    zz, h, mu = sp.symbols("z h mu", positive=True)
    for p in (xs * (1 + 1 / (xs ** 2 + ys ** 2)), sp.exp(xs) * sp.cos(ys), xs ** 2 - ys ** 2):
        assert z0(sp.diff(p, xs, 2) + sp.diff(p, ys, 2))
        phi = -zz * (h - zz) / (2 * mu) * p
        assert z0(sp.diff(phi, xs, 2) + sp.diff(phi, ys, 2))
        u = -1 / (2 * mu) * sp.diff(p, xs) * zz * (h - zz)
        assert z0(u - sp.diff(phi, xs))


def test_hele_shaw_cylinder_V1_equals_ideal_cylinder():  # V1 (N36): gap average = the ch06/ch04 ideal cylinder
    U, a, h, mu = 0.02, 0.01, 1e-3, 1e-3
    th = RNG.uniform(0, 2 * np.pi, 40)
    rr = RNG.uniform(1.05 * a, 5 * a, 40)
    X, Y = rr * np.cos(th), rr * np.sin(th)
    hs = ch08.hele_shaw_cylinder(X, Y, None, U, a, h, mu)
    vel = NS.exact_solution("cylinder", np.stack([X, Y]), U=U, a=a)[0]
    assert maxrel(hs["u_mean"], vel[0]) < 1e-12 and maxrel(hs["v_mean"], vel[1]) < 1e-12
    # at height z: parabolic factor 6z(h − z)/h², unit mean across the gap; φ = −z(h − z)p/2μ
    for zf in (0.1, 0.5, 0.9):
        hz = ch08.hele_shaw_cylinder(X, Y, zf * h, U, a, h, mu)
        assert maxrel(hz["u"], 6 * zf * (1 - zf) * vel[0]) < 1e-12
        assert maxrel(hz["phi"], -zf * h * (h - zf * h) / (2 * mu) * np.asarray(hs["p"])) < 1e-12
    # the pressure drives the mean flow: ū = −(h²/12μ)∂p/∂x (finite differences)
    d = 1e-6
    px = (np.asarray(ch08.hele_shaw_cylinder(X + d, Y, None, U, a, h, mu)["p"]) - np.asarray(ch08.hele_shaw_cylinder(X - d, Y, None, U, a, h, mu)["p"])) / (2 * d)
    assert maxrel(-h ** 2 / (12 * mu) * px, hs["u_mean"]) < 1e-6
    assert np.isnan(ch08.hele_shaw_cylinder(0.0, 0.0, None, U, a, h, mu)["p"])


def test_hele_shaw_grid_V3_staircase_first_order():  # V3 (N36): masked five-point Laplace solve converges (staircase ⇒ order 1)
    hs, es = [], []
    for n in (33, 65, 129, 257):
        g = ch08.hele_shaw_streamfunction_grid(n)
        hs.append(g["h"])
        es.append(g["err_far"])
    p = observed_order(hs, es)
    assert abs(p - 1.0) < ORDER_TOL  # the staircase disk boundary is first order (documented in the docstring)
    assert es[-1] < 0.015 and all(np.diff(es) < 0)


# =====================================================================================================================
# C07 — slider bearing (Example 8.1), exact load, optimum taper
# =====================================================================================================================
def test_slider_bearing_V2_engine_and_printed_slips():  # V2 (C07): sympy engine residuals; the two printed slips fail
    d = ch08.slider_bearing_sympy()
    assert all(v == 0 for v in d["residuals"].values())
    assert d["ode_residual_exact"] == 0 and d["bc_residuals"] == (0, 0)
    assert d["ode_residual_book"] != 0 and d["printed_final_satisfies_ode"] is False  # R8 (ii)
    assert d["printed_intermediate_satisfies_ends"] is False  # R8 (i): (1 − αx/L) integrands


def test_slider_bearing_V2_derivation():  # V2 derivation D14 ★★★: Part F steps 3–15, independent of the engine
    x, y = sp.symbols("x y", real=True)
    L, h0, mu, U, al = sp.symbols("L h_0 mu U alpha", positive=True)
    C1, C2, pe, px = sp.symbols("C_1 C_2 p_e p_x", real=True)
    hs = sp.Symbol("h", positive=True)
    u = -hs ** 2 / (2 * mu) * px * (y / hs) * (1 - y / hs) + U * y / hs  # (8.19), U0 = 0, U_h = U
    flux = sp.integrate(u, (y, 0, hs))
    assert z0(flux - (-hs ** 3 / (12 * mu) * px + U * hs / 2))  # step 3
    c1eq = flux - U * hs  # step 4: ∫(u − U)dy
    assert z0(c1eq - (-hs ** 3 / (12 * mu) * px - U * hs / 2))
    dpdx = sp.solve(sp.Eq(C1, c1eq), px)[0]  # step 5
    assert z0(dpdx - (-12 * mu * C1 / hs ** 3 - 6 * mu * U / hs ** 2))
    hx = h0 * (1 + al * x / L)
    dpdx_x = dpdx.subs(hs, hx)  # step 6
    P = 6 * mu * L * C1 / (al * h0 ** 3) * (1 + al * x / L) ** -2 + 6 * mu * U * L / (al * h0 ** 2) * (1 + al * x / L) ** -1 + C2
    assert z0(sp.diff(P, x) - dpdx_x)  # step 7: antiderivatives of (1 + αx/L)^(−3), ^(−2)
    e0 = sp.Eq(P.subs(x, 0), pe)  # step 8
    assert z0(e0.lhs - (6 * mu * L / (al * h0 ** 2) * (C1 / h0 + U) + C2))
    eL = sp.Eq(P.subs(x, L), pe)  # step 9
    assert z0(eL.lhs - (6 * mu * L / (al * h0 ** 2) * (C1 / (h0 * (1 + al) ** 2) + U / (1 + al)) + C2))
    sol = sp.solve([e0, eL], [C1, C2], dict=True)[0]
    assert z0(sol[C1] + (1 + al) / (2 + al) * U * h0)  # step 10
    assert z0(sol[C2] - (pe - 6 * mu * L * U / (al * h0 ** 2) / (2 + al)))  # step 11
    s = sp.Symbol("s", positive=True)
    assert z0(sp.expand(-(1 + al) + (2 + al) * (1 + s) - (1 + s) ** 2 - s * (al - s)))  # step 12's numerator
    p_ex = sp.simplify(P.subs(sol))
    step12 = pe + 6 * mu * L * U / (al * h0 ** 2) * (s * (al - s)) / ((2 + al) * (1 + s) ** 2)
    assert z0(p_ex - step12.subs(s, al * x / L))
    final = pe + 6 * mu * L * U / h0 ** 2 * al * (x / L) * (1 - x / L) / ((2 + al) * (1 + al * x / L) ** 2)  # step 13
    assert z0(p_ex - final)
    printed = pe + 6 * mu * L * U / h0 ** 2 * al * (x / L) * (1 - x / L) / ((2 + al) * (1 + al * x / L))
    assert not z0(sp.diff(printed, x) - dpdx_x.subs(C1, sol[C1]))  # the first power fails step 5's ODE
    lin = sp.series(final - pe, al, 0, 2).removeO()  # step 14
    assert z0(lin - 3 * al * mu * L * U / h0 ** 2 * (x / L) * (1 - x / L))
    assert z0(sp.integrate(lin, (x, 0, L)) - al * mu * L ** 2 * U / (2 * h0 ** 2))  # step 15
    # exact load (N35, ours) and its small-α series
    W = sp.integrate(sp.expand(final - pe), (x, 0, L))
    Wex = 6 * mu * U * L ** 2 / (h0 ** 2 * al ** 2) * (sp.log(1 + al) - 2 * al / (2 + al))
    assert abs(float((W - Wex).subs({mu: 0.05, U: 5, L: 0.05, h0: 5e-5, al: sp.Rational(3, 10)}))) < 1e-8 * float(
        Wex.subs({mu: 0.05, U: 5, L: 0.05, h0: 5e-5, al: sp.Rational(3, 10)}))
    ser = sp.series(Wex * h0 ** 2 / (mu * U * L ** 2), al, 0, 6).removeO()
    assert z0(ser - (al / 2 - 3 * al ** 2 / 4 + sp.Rational(33, 40) * al ** 3 - sp.Rational(13, 16) * al ** 4 + sp.Rational(171, 224) * al ** 5))


def test_slider_bearing_V1_quadrature_load_and_numbers():  # V1 (C07, N35): ∫p dx = W, linear → exact as α → 0; D14 numbers
    L, h0, U, mu = 0.05, 50e-6, 5.0, 0.05
    for al in (1e-4, 5e-4, 2e-3, 0.1, 0.5, 1.2, 3.0):
        Wq = quad(lambda s: float(ch08.slider_bearing(s, h0, al, L, U, mu)), 0, L, epsabs=0, epsrel=1e-13)[0]
        assert ch08.slider_bearing_load(h0, al, L, U, mu, "exact") == pytest.approx(Wq, rel=1e-9)
        Wl = quad(lambda s: float(ch08.slider_bearing(s, h0, al, L, U, mu, model="linear")), 0, L, epsrel=1e-13)[0]
        assert ch08.slider_bearing_load(h0, al, L, U, mu, "linear") == pytest.approx(Wl, rel=1e-12)
    # the series branch (|α| < 1e-3) joins the closed form continuously
    a_ = 1e-3
    assert ch08.slider_bearing_load(h0, a_ * (1 - 1e-9), L, U, mu) == pytest.approx(ch08.slider_bearing_load(h0, a_ * (1 + 1e-9), L, U, mu), rel=1e-8)
    assert ch08.slider_bearing_load(h0, 1e-6, L, U, mu) / ch08.slider_bearing_load(h0, 1e-6, L, U, mu, "linear") == pytest.approx(1.0, abs=2e-6)
    # D14 / N35 numbers: α = 0.1 → 12.5 (linear) vs 10.8 kN/m (exact); α = 0.5 → 62.5 vs 32.8 kN/m
    assert ch08.slider_bearing_load(h0, 0.1, L, U, mu, "linear") == pytest.approx(12500.0)
    assert ch08.slider_bearing_load(h0, 0.1, L, U, mu) == pytest.approx(10.8e3, rel=5e-3)
    assert ch08.slider_bearing_load(h0, 0.5, L, U, mu, "linear") == pytest.approx(62500.0)
    assert ch08.slider_bearing_load(h0, 0.5, L, U, mu) == pytest.approx(32.8e3, rel=5e-3)
    # the printed first-power pressure carries a different (wrong) load
    assert ch08.slider_bearing_load(h0, 0.5, L, U, mu, "book") > 1.05 * ch08.slider_bearing_load(h0, 0.5, L, U, mu)
    with pytest.raises(ValueError):
        ch08.slider_bearing_load(h0, -1.0, L, U, mu)


def test_slider_bearing_V7_reversal_and_ends():  # V7 (C07): p(0) = p(L) = p_e; αU < 0 sucks the pad down; p ∝ U
    L, h0, mu = 0.05, 50e-6, 0.05
    x = np.linspace(0, L, 51)
    for model in ("exact", "linear", "book"):
        p = np.asarray(ch08.slider_bearing(x, h0, 0.4, L, 3.0, mu, p_e=1e5, model=model))
        assert p[0] == pytest.approx(1e5) and p[-1] == pytest.approx(1e5)
    pp = np.asarray(ch08.slider_bearing(x, h0, 0.4, L, 3.0, mu))
    pm = np.asarray(ch08.slider_bearing(x, h0, 0.4, L, -3.0, mu))
    assert maxrel(pm, -pp) < 1e-15 and ch08.slider_bearing_load(h0, 0.4, L, -3.0, mu) < 0
    assert ch08.slider_bearing_load(h0, -0.3, L, 3.0, mu) < 0  # a gap narrowing the wrong way
    # load ∝ 1/h0² (stability of the bearing)
    assert ch08.slider_bearing_load(h0 / 2, 0.4, L, 3.0, mu) == pytest.approx(4 * ch08.slider_bearing_load(h0, 0.4, L, 3.0, mu))


@needs_ref
def test_slider_optimum_taper_V5_san_andres():  # V5 (N35): K_opt 2.1889, W* 0.0267; P_max at K = 2.414, 0.043
    ref = ref_json()["san_andres_slider"]
    opt = ch08.slider_optimum_taper()
    assert opt["K_opt"] == pytest.approx(ref["K_opt"], rel=1e-2)
    assert opt["W_star"] == pytest.approx(ref["W_opt"], rel=1e-2)
    assert abs(opt["K_opt"] - ref["K_opt"]) < 5e-4  # also to the printed 5 figures
    # San Andrés' W(K) equals our exact load normalised by 6μUL²/h_exit² with K = 1 + α (identity on a sweep)
    W_SA = lambda K: (np.log(K) + 2 * (1 - K) / (1 + K)) / (1 - K) ** 2  # noqa: E731
    for al in (0.05, 0.4, 1.1887, 2.0, 5.0):
        assert ch08.slider_bearing_load(1.0, al, 1.0, 1.0, 1.0) / 6.0 == pytest.approx(W_SA(1 + al), rel=1e-12)
    # peak pressure P_max(K) = (K − 1)/(4K(1 + K)) and its maximum 0.043 at K = 2.414
    Pm = lambda al: ch08.slider_bearing_state(1.0, al, 1.0, 1.0, 1.0)["dp_max"] / 6.0  # noqa: E731
    for al in (0.2, 1.0, 1.414, 3.0):
        K = 1 + al
        assert Pm(al) == pytest.approx((K - 1) / (4 * K * (1 + K)), rel=1e-12)
    from scipy.optimize import minimize_scalar
    r = minimize_scalar(lambda a: -Pm(a), bounds=(0.2, 5.0), method="bounded", options={"xatol": 1e-10})
    assert 1 + r.x == pytest.approx(ref["K_Pmax"], rel=1e-3) and -r.fun == pytest.approx(ref["Pmax_max"], rel=1e-2)


def test_slider_bearing_state_V1_explainer_numbers():  # V1 (E3/IF4): peak location, C1, flux at every x, inlet backflow
    L, h0, U, mu = 0.05, 50e-6, 5.0, 0.05
    for al in (0.1, 0.8, 1.5):
        st = ch08.slider_bearing_state(h0, al, L, U, mu)
        x = np.linspace(0, L, 20001)
        p = np.asarray(ch08.slider_bearing(x, h0, al, L, U, mu))
        assert st["dp_max"] == pytest.approx(p.max(), rel=1e-7) and st["x_pmax"] == pytest.approx(x[np.argmax(p)], abs=L / 20000)
        assert st["C1"] == pytest.approx(-(1 + al) / (2 + al) * U * h0)
        assert st["W_exact"] == pytest.approx(ch08.slider_bearing_load(h0, al, L, U, mu))
        assert st["err_linear"] == pytest.approx(100 * (st["W_linear"] / st["W_exact"] - 1))
        assert st["inlet_backflow"] is (al > 1)
        # pad-frame flux ∫(u − U)dy = C1 at every station (slider_gap_velocity)
        for xs in (0.0, 0.3 * L, L):
            hx = h0 * (1 + al * xs / L)
            q = quad(lambda s: float(ch08.slider_gap_velocity(xs, s, h0, al, L, U, mu, frame="pad")), 0, hx, epsrel=1e-12)[0]
            assert q == pytest.approx(st["C1"], rel=1e-9)
    assert np.isnan(ch08.slider_gap_velocity(0.01, 1.0, h0, 0.3, L, U, mu))  # outside the gap
    assert ch08.slider_gap_velocity(0.01, 0.0, h0, 0.3, L, U, mu) == 0.0  # the floor
    assert ch08.slider_gap_velocity(0.0, h0, h0, 0.3, L, U, mu) == pytest.approx(U)  # the pad


# =====================================================================================================================
# C08 — thin film (Example 8.3): profile, flux, the thin-film equation, the spreading bead
# =====================================================================================================================
def test_thin_film_V1_profile_walls_and_flux():  # V1 (C08): u(0) = 0, ∂u/∂y(h) = 0, q = −(ρg/3μ)h³h_x
    h, hx, rho, g, mu = 3e-3, -0.02, 1260.0, 9.81, 1.4
    assert ch08.thin_film_velocity(0.0, h, hx, rho, g, mu) == 0.0
    dy = 1e-9
    top = (float(ch08.thin_film_velocity(h + dy, h, hx, rho, g, mu)) - float(ch08.thin_film_velocity(h - dy, h, hx, rho, g, mu))) / (2 * dy)
    assert abs(top) < 1e-9 * abs(rho * g / mu * hx * h)
    q = quad(lambda s: float(ch08.thin_film_velocity(s, h, hx, rho, g, mu)), 0, h, epsrel=1e-13)[0]
    assert float(ch08.thin_film_flux(h, hx, rho, g, mu)) == pytest.approx(q, rel=1e-12)
    assert float(ch08.thin_film_flux(h, hx, rho, g, mu)) > 0  # flows down the slope (h_x < 0 → +x)
    # wrong variant h³/(2μ): 1.5 times the flux
    assert rho * g / (2 * mu) * h ** 3 * abs(hx) == pytest.approx(1.5 * q)


def test_thin_film_equation_V2_derivation():  # V2 derivation D15 ★★: hydrostatics + (8.18) + stress-free top + CV mass
    x, y, t = sp.symbols("x y t", real=True)
    rho, g, mu, pa = sp.symbols("rho g mu p_a", positive=True)
    h = sp.Function("h")(x, t)
    A, B = sp.symbols("A B", real=True)
    p = pa + rho * g * (h - y)  # step 3
    px = sp.diff(p, x)
    assert z0(px - rho * g * sp.diff(h, x))  # step 4
    u = px / mu * y ** 2 / 2 + A * y + B  # step 5 ((8.18))
    Bv = sp.solve(u.subs(y, 0), B)[0]  # step 6
    Av = sp.solve(sp.diff(u, y).subs({y: h, B: Bv}), A)[0]  # step 7: no stress at the free surface
    assert z0(Av + rho * g / mu * h * sp.diff(h, x))
    uu = u.subs({A: Av, B: Bv})
    assert z0(uu + rho * g / (2 * mu) * sp.diff(h, x) * y * (2 * h - y))  # step 8
    q = sp.integrate(uu, (y, 0, h))
    assert z0(q + rho * g / (3 * mu) * h ** 3 * sp.diff(h, x))  # step 9
    ht = -sp.diff(q, x)  # steps 1–2, 10: h_t + q_x = 0
    assert z0(ht - rho * g / (3 * mu) * sp.diff(h ** 3 * sp.diff(h, x), x))
    # units: ρg/(3μ) is 1/(m s); ∂(h³h_x)/∂x is m²
    beta = Q_(1000, "kg/m**3") * Q_(9.81, "m/s**2") / Q_(1.0, "Pa*s")
    assert (beta * Q_(1, "m**2")).dimensionality == Q_(1, "m/s").dimensionality


def test_viscous_current_similarity_V1_pde_volume_front():  # V1 (N63): Huppert's form solves the PDE and holds ∫h = A
    x, t, A, beta = sp.symbols("x t A beta", positive=True)
    eN = sp.Symbol("eta_N", positive=True)
    xN = eN * (beta * A ** 3 * t) ** sp.Rational(1, 5)
    hc = sp.Rational(3, 10) ** sp.Rational(1, 3) * eN ** sp.Rational(2, 3) * (A ** 2 / beta) ** sp.Rational(1, 5) * t ** sp.Rational(-1, 5)
    h = hc * (1 - x ** 2 / xN ** 2) ** sp.Rational(1, 3)
    assert z0((sp.diff(h, t) - beta * sp.diff(h ** 3 * sp.diff(h, x), x)) / h)
    rho, g, mu, area = 1000.0, 9.81, 1.0, 1e-4
    for tt in (1.0, 100.0):
        hh, xn = ch08.viscous_current_similarity(0.0, tt, area, rho, g, mu, return_front=True)
        v = quad(lambda s: float(ch08.viscous_current_similarity(s, tt, area, rho, g, mu)), 0, float(xn), limit=200)[0]
        assert v == pytest.approx(area, rel=1e-8)
        assert float(ch08.viscous_current_similarity(1.001 * float(xn), tt, area, rho, g, mu)) == 0.0
        st = ch08.thin_film_state(tt, area, rho, g, mu)
        assert st["x_N"] == pytest.approx(float(xn)) and st["h_centre"] == pytest.approx(float(hh))
        assert st["front_speed"] == pytest.approx(float(xn) / (5 * tt)) and st["beta"] == pytest.approx(rho * g / (3 * mu))
    # t^{1/5} and t^{−1/5} exactly
    r1 = ch08.viscous_current_similarity(0.0, 32.0, area, rho, g, mu, return_front=True)
    r0 = ch08.viscous_current_similarity(0.0, 1.0, area, rho, g, mu, return_front=True)
    assert float(r1[1]) / float(r0[1]) == pytest.approx(2.0, rel=1e-13) and float(r1[0]) / float(r0[0]) == pytest.approx(0.5, rel=1e-13)


@needs_ref
def test_viscous_current_eta_N_V5_huppert():  # V5 (N63): η_N = 1.411 (Huppert 1982, Ball & Huppert Appendix a)
    ref = ref_json()["huppert_planar_current"]
    eN = ch08.viscous_current_eta_N()
    assert eN == pytest.approx(ref["eta_N"], rel=1e-3)
    ind = ((1 / 5) * (3 / 10) ** (1 / 3) * np.sqrt(np.pi) * gamma_fn(1 / 3) / gamma_fn(5 / 6)) ** (-3 / 5)
    assert eN == pytest.approx(ind, rel=1e-14)
    # the volume condition fixes η_N: ∫₀¹(1 − y²)^{1/3} dy · (3/10)^{1/3} η_N^{5/3} = 1
    I = quad(lambda s: (1 - s * s) ** (1 / 3), 0, 1, epsrel=1e-13)[0]
    assert I * 0.3 ** (1 / 3) * eN ** (5 / 3) == pytest.approx(1.0, rel=1e-10)


def test_thin_film_spread_V4_volume_conserved():  # V4 (C08): Σh Δx is conserved to the Newton tolerance
    rho, g, mu = 1000.0, 9.81, 1.0
    N, X = 200, 0.3
    dx = 2 * X / N
    x = -X + dx * (np.arange(N) + 0.5)
    h0 = np.where(np.abs(x) < 0.02, 0.004, 0.0) + 0.001 * np.exp(-((x - 0.05) / 0.01) ** 2)
    res = ch08.thin_film_spread(h0, x, [1.0, 10.0, 100.0], rho, g, mu)
    assert np.max(np.abs(res["volume"] / (np.sum(np.maximum(h0, 1e-6)) * dx) - 1)) < 1e-12
    assert np.all(res["h"] > 0) and np.all(np.diff(res["x_front"]) > 0)  # positivity, the front advances


@needs_ref
def test_thin_film_spread_V5_huppert_shape_and_t_one_fifth():  # V5 (C08/N63): the run approaches Huppert's solution
    rho, g, mu = 1000.0, 9.81, 1.0
    # (a) from a box: the front grows as t^{1/5} once the initial shape is forgotten
    N, X = 400, 0.5
    dx = 2 * X / N
    x = -X + dx * (np.arange(N) + 0.5)
    h0 = np.where(np.abs(x) < 0.02, 0.005, 0.0)
    area = float(np.sum(h0) * dx / 2)
    ts = np.logspace(0, 3, 13)
    res = ch08.thin_film_spread(h0, x, ts, rho, g, mu)
    slope = np.polyfit(np.log(ts[-5:]), np.log(res["x_front"][-5:]), 1)[0]
    assert abs(slope - 0.2) < 0.005
    xN = float(ch08.viscous_current_similarity(0.0, ts[-1], area, rho, g, mu, return_front=True)[1])
    assert res["x_front"][-1] / xN == pytest.approx(1.0, abs=0.03)
    # (b) started on the similarity profile at t0 = 10 s, after 70 s (the solver clock starts at 0) the L2 shape error < 2 %
    area = 1e-4
    t0, t1 = 10.0, 80.0
    xN1 = float(ch08.viscous_current_similarity(0.0, t1, area, rho, g, mu, return_front=True)[1])
    for N in (200, 400, 800):
        dx = 2 * 1.6 * xN1 / N
        x = -1.6 * xN1 + dx * (np.arange(N) + 0.5)
        run = ch08.thin_film_spread(ch08.viscous_current_similarity(x, t0, area, rho, g, mu), x, [t1 - t0], rho, g, mu, h_min=1e-7)
        hex_ = np.asarray(ch08.viscous_current_similarity(x, t1, area, rho, g, mu))
        assert np.linalg.norm(run["h"][-1] - hex_) / np.linalg.norm(hex_) < 0.02
        assert run["x_front"][-1] / xN1 == pytest.approx(1.0, abs=0.02)


# =====================================================================================================================
# C09 — Stokes' first problem (8.20)–(8.31)
# =====================================================================================================================
def test_stokes_first_V1_boundaries_collapse_and_parity():  # V1 (C09, N52, N55): BCs, η-collapse, ch04 parity
    y = np.linspace(0, 0.05, 201)
    assert ch08.stokes_first_problem(0.0, 3.0, 2.0, 1e-6) == pytest.approx(2.0)
    assert ch08.stokes_first_problem(1.0, 3.0, 2.0, 1e-6) < 1e-100  # u → 0 far away (erfc, no cancellation)
    assert np.all(np.asarray(ch08.stokes_first_problem(y, 0.0, 1.0, 1e-6)) == 0.0)  # (8.21): at rest for t ≤ 0
    eta = np.linspace(0, 6, 61)
    profiles = []
    for U, nu, t in ((1.0, 1e-6, 10.0), (0.3, 1.5e-5, 2.0), (4.0, 1e-4, 1000.0)):
        profiles.append(np.asarray(ch08.stokes_first_problem(eta * np.sqrt(nu * t), t, U, nu)) / U)
    assert np.max(np.ptp(np.array(profiles), axis=0)) < 1e-14  # the collapse (N55)
    assert maxrel(profiles[0], erfc(eta / 2)) < 1e-15
    assert ch04.stokes_first_problem is ch08.stokes_first_problem  # promoted: ch04 re-exports it
    X = np.stack([np.zeros_like(y), y])
    assert maxrel(ch08.stokes_first_problem(y, 5.0, 0.7, 1e-6), NS.exact_solution("stokes_first", X, 5.0, U=0.7, nu=1e-6)[0][0]) < 1e-14
    # η = y/√(νt) (8.25) and the figures' axis η/2 (R22); the missing-2 variant is rejected in the V2 test below
    assert ch08.similarity_variable(0.02, 100.0, 1e-6) == pytest.approx(2.0) and ch08.similarity_variable(0.02, 100.0, 1e-6, half=True) == pytest.approx(1.0)


def test_stokes_first_V2_sympy_residual_and_wrong_variant():  # V2 (C09): (8.30) solves (8.20); erfc(y/√νt) does not
    d = ch08.similarity_reduce_sympy("stokes1")
    assert d["residual"] == 0 and d["pde_residual"] == 0 and d["matches_8_26"]
    y, t, nu, U = sp.symbols("y t nu U", positive=True)
    good = U * sp.erfc(y / (2 * sp.sqrt(nu * t)))
    bad = U * sp.erfc(y / sp.sqrt(nu * t))
    assert z0(sp.diff(good, t) - nu * sp.diff(good, y, 2))
    assert not z0(sp.diff(bad, t) - nu * sp.diff(bad, y, 2))


def test_stokes_first_V2_derivation_D17_D18_D19():  # V2 derivations D17, D18, D19 ★★: Π groups → η, (8.26), erf
    # D17: five variables, two dimensions → three groups; linearity removes y/Ut
    g = ch08.stokes_first_pi_groups()
    assert len(g) == 3
    dims = {"u": (1, -1), "U": (1, -1), "y": (1, 0), "t": (0, 1), "nu": (2, -1)}
    for grp in g:
        tot = np.sum([np.array(dims[k]) * float(v) for k, v in grp.items()], axis=0)
        assert np.allclose(tot, 0.0)
    y, t, nu, U = sp.symbols("y t nu U", positive=True)
    uU = sp.erfc(y / (2 * sp.sqrt(nu * t)))  # u/U does not depend on U (step 6)
    assert sp.diff(uU, U) == 0
    # D18: ∂η/∂t = −η/(2t); chain rule; t drops out
    eta = y / sp.sqrt(nu * t)
    assert z0(sp.diff(eta, t) + eta / (2 * t))
    assert z0(sp.diff(eta, y) - 1 / sp.sqrt(nu * t))
    F = sp.Function("F")
    e = sp.Symbol("eta", positive=True)
    u = U * F(eta)
    lhs = sp.diff(u, t)
    rhs = nu * sp.diff(u, y, 2)
    red = sp.simplify(((lhs - rhs) * t / U).subs(y, e * sp.sqrt(nu * t)).doit())
    target = -e / 2 * sp.Derivative(F(e), e) - sp.Derivative(F(e), (e, 2))
    assert z0(red - target)
    # D19: F′ = A e^{−η²/4}, F = A∫₀^η e^{−ξ²/4}dξ + B; B = 1, A = −1/√π; F = erfc(η/2)
    Aa, Bb, xi = sp.symbols("A B xi", real=True)
    G = Aa * sp.exp(-e ** 2 / 4)
    assert z0(sp.diff(G, e) + e / 2 * G)  # steps 1–4
    assert z0(sp.integrate(sp.exp(-xi ** 2 / 4), (xi, 0, sp.oo)) - sp.sqrt(sp.pi))  # steps 7–9: 2·√π/2
    Av = sp.solve(Aa * sp.sqrt(sp.pi) + 1, Aa)[0]
    assert z0(Av + 1 / sp.sqrt(sp.pi))
    Fsol = Av * sp.integrate(sp.exp(-xi ** 2 / 4), (xi, 0, e)) + 1
    assert z0((Fsol - sp.erfc(e / 2)).rewrite(sp.erf))  # steps 10–11: (8.30)
    ode = sp.dsolve(sp.Eq(-e / 2 * F(e).diff(e), F(e).diff(e, 2)), F(e)).rhs
    C1, C2 = sorted(ode.free_symbols - {e}, key=str)
    csol = sp.solve([ode.subs(e, 0) - 1, sp.limit(ode, e, sp.oo)], [C1, C2], dict=True)[0]
    assert z0((ode.subs(csol) - sp.erfc(e / 2)).rewrite(sp.erf))
    # check line: y = 1 cm, t = 100 s, ν = 1e-6 → η = 1, u/U = erfc(0.5) = 0.4795
    assert ch08.stokes_first_problem(0.01, 100.0, 1.0, 1e-6) == pytest.approx(0.4795, abs=1e-4)


def test_similarity_ode_solve_V3_bvp_matches_closed_forms():  # V3 (N46, N61): solve_bvp ≤ 1e-7, insensitive to η_max
    for case in ("stokes1", "line_vortex"):
        errs = []
        for em in (10.0, 12.0, 14.0):
            d = ch08.similarity_ode_solve(case, eta_max=em, full=True)
            assert d["success"] and d["max_err"] < 1e-7
            errs.append(d["max_err"])
        from scipy.interpolate import CubicSpline
        e1, F1 = ch08.similarity_ode_solve(case, eta_max=10.0)
        e2, F2 = ch08.similarity_ode_solve(case, eta_max=14.0)
        common = np.linspace(max(e1[0], e2[0]), 6.0, 50)
        assert np.max(np.abs(CubicSpline(e1, F1)(common) - CubicSpline(e2, F2)(common))) < 2e-7
    # line vortex: F(η0 = 1e-4) = 0 truncation is O(η0²)
    d = ch08.similarity_ode_solve("line_vortex", full=True)
    assert d["max_err"] < 1e-8


def test_crank_nicolson_V3_second_order_space_and_time():  # V3 (N39, I25): CN order 2; backward Euler order 1 (wrong variant)
    nu, U, t0, t1 = 1e-6, 1.0, 100.0, 400.0
    Ly = 12 * np.sqrt(nu * t1)  # ≥ 9.5√(νt): the far Dirichlet 0 is exact to 1e-16
    y = np.linspace(0, Ly, 801)
    u0 = ch08.stokes_first_problem(y, t0, U, nu)
    for theta, design in ((0.5, 2.0), (1.0, 1.0)):
        uref = DIF.crank_nicolson_1d(u0, y, (t1 - t0) / 5120, 5120, nu, U, 0.0, startup_be=0, t0=t0, theta=theta)
        dts, errs = [], []
        for n in (10, 20, 40, 80):
            u = DIF.crank_nicolson_1d(u0, y, (t1 - t0) / n, n, nu, U, 0.0, startup_be=0, t0=t0, theta=theta)
            errs.append(np.max(np.abs(u - uref)))
            dts.append((t1 - t0) / n)
        assert abs(observed_order(dts, errs) - design) < ORDER_TOL, (theta, pairwise_orders(dts, errs))
    hs, errs = [], []
    for ny in (41, 81, 161, 321):
        yy = np.linspace(0, Ly, ny)
        u = DIF.crank_nicolson_1d(ch08.stokes_first_problem(yy, t0, U, nu), yy, (t1 - t0) / 2000, 2000, nu, U, 0.0, startup_be=0, t0=t0)
        errs.append(np.max(np.abs(u - ch08.stokes_first_problem(yy, t1, U, nu))))
        hs.append(yy[1] - yy[0])
    assert abs(observed_order(hs, errs) - 2.0) < ORDER_TOL


def test_crank_nicolson_V3_impulsive_start_and_stability():  # V3 (C09): Rannacher start keeps order 2; stable at 50× FTCS
    nu, U, T = 1e-6, 1.0, 100.0
    Ly = 12 * np.sqrt(nu * T)
    dts, e_be, e_cn = [], [], []
    for n in (50, 100, 200, 400):
        y = np.linspace(0, Ly, 4 * n + 1)
        ex = ch08.stokes_first_problem(y, T, U, nu)
        e_be.append(np.max(np.abs(DIF.crank_nicolson_1d(np.zeros_like(y), y, T / n, n, nu, U, 0.0, startup_be=2) - ex)))
        e_cn.append(np.max(np.abs(DIF.crank_nicolson_1d(np.zeros_like(y), y, T / n, n, nu, U, 0.0, startup_be=0) - ex)))
        dts.append(T / n)
    assert abs(observed_order(dts, e_be) - 2.0) < ORDER_TOL
    assert pairwise_orders(dts, e_cn)[-1] < 1.0  # pure CN stalls on the impulsive jump (the reason for startup_be)
    y = np.linspace(0, Ly, 201)
    dt = 50 * (y[1] - y[0]) ** 2 / (2 * nu)
    u = DIF.crank_nicolson_1d(np.zeros_like(y), y, dt, int(T / dt) + 1, nu, U, 0.0)
    assert np.all(np.isfinite(u)) and np.max(np.abs(u)) <= U * (1 + 1e-12)
    # FTCS twin (reuse, ch01) below its limit agrees with (8.30)
    dy = y[1] - y[0]
    dtf = DIF.stable_time_step(nu, dy, 0.4)
    nst = int(round(T / dtf))
    f0 = np.zeros_like(y)
    f0[0] = U
    uf = DIF.ftcs_diffusion_1d(f0, nu, dy, T / nst, nst, values=(U, 0.0), save_every=nst)[-1]
    assert np.max(np.abs(np.asarray(uf) - ch08.stokes_first_problem(y, T, U, nu))) < 2e-2


def test_vorticity_content_V1_plus_U_and_vorticity_field():  # V1 (N53): ∫ω dy = +U for all t (printed −U rejected)
    for t in (0.1, 10.0, 1e4):
        val, err = ch08.vorticity_content(t, 0.7, 1e-6, return_error=True)
        assert val == pytest.approx(0.7, rel=1e-10) and err < 1e-10
        assert val != pytest.approx(-0.7)  # R10
    y = np.linspace(1e-4, 0.02, 60)
    dy = 1e-8
    fd = -(np.asarray(ch08.stokes_first_problem(y + dy, 50.0, 0.7, 1e-6)) - np.asarray(ch08.stokes_first_problem(y - dy, 50.0, 0.7, 1e-6))) / (2 * dy)
    assert maxrel(ch08.stokes_first_vorticity(y, 50.0, 0.7, 1e-6), fd) < 1e-6
    st = ch08.stokes_first_state(100.0, 1.0, 1e-6)
    assert st["vorticity_content"] == 1.0 and st["omega_wall"] == pytest.approx(float(ch08.stokes_first_vorticity(0.0, 100.0, 1.0, 1e-6)))
    assert st["tau_w"] == pytest.approx(1000 * 1e-6 * st["omega_wall"])


def test_diffusion_thickness_V1_level_inversion_and_D20_numbers():  # V1 (N54, D20): erfc(δ/2√νt) = level; 3.643
    for lev in (0.01, 0.05, 0.5):
        d = ch08.diffusion_thickness(100.0, 1e-6, lev)
        assert erfc(d / (2 * np.sqrt(1e-4))) == pytest.approx(lev, rel=1e-12)
        eta = brentq(lambda e: erfc(e / 2) - lev, 0, 20, xtol=1e-14)
        assert d / 1e-2 == pytest.approx(eta, rel=1e-10)
    assert ch08.diffusion_thickness(100.0, 1e-6) == pytest.approx(0.036428, rel=1e-4)  # 3.64 cm
    assert ch08.diffusion_thickness(3600.0, 1e-6) == pytest.approx(0.2186, rel=1e-3)  # 21.9 cm
    assert 2 * erfcinv(0.05) == pytest.approx(2.7718, rel=1e-4)  # the 95 % convention of D20's check
    with pytest.raises(ValueError):
        ch08.diffusion_thickness(1.0, 1e-6, 1.5)
    st = ch08.stokes_first_state(100.0)
    assert st["eta_edge"] == pytest.approx(3.6428, rel=1e-4) and st["delta"] == pytest.approx(st["eta_edge"] * st["sqrt_nut"])


def test_stokes_first_stopped_V1_superposition():  # V1 (N56): equal to (8.30) for t ≤ T; u(0, t > T) = 0; diffusion
    y = np.linspace(0, 0.02, 81)
    T, nu = 50.0, 1e-6
    assert maxrel(ch08.stokes_first_stopped(y, 30.0, T, 1.0, nu), ch08.stokes_first_problem(y, 30.0, 1.0, nu)) < 1e-15
    assert ch08.stokes_first_stopped(0.0, 80.0, T, 1.0, nu) == pytest.approx(0.0, abs=1e-15)
    # diffusion residual (finite differences) for t > T
    t, dt = 120.0, 1e-3
    f = lambda yy, tt: np.asarray(ch08.stokes_first_stopped(yy, tt, T, 1.0, nu))  # noqa: E731
    yi = y[5:-5]
    ut = np.max(np.abs((f(yi, t + dt) - f(yi, t - dt)) / (2 * dt)))
    rs = []
    for dy in (4e-4, 2e-4, 1e-4):  # the residual is pure O(Δy²) truncation of the second difference
        res = (f(yi, t + dt) - f(yi, t - dt)) / (2 * dt) - nu * (f(yi + dy, t) - 2 * f(yi, t) + f(yi - dy, t)) / dy ** 2
        rs.append(np.max(np.abs(res)))
    assert abs(observed_order([4e-4, 2e-4, 1e-4], rs) - 2.0) < ORDER_TOL and rs[-1] < 1e-4 * ut
    # the momentum ∫u dy grows as 2U√(νt/π) before T
    m = quad(lambda s: float(ch08.stokes_first_stopped(s, 30.0, T, 1.0, nu)), 0, 0.1, epsrel=1e-12)[0]
    assert m == pytest.approx(2 * np.sqrt(nu * 30.0 / np.pi), rel=1e-9)


# =====================================================================================================================
# C10 — the similarity ansatz (8.32): Examples 8.4–8.7, vortex sheet, temporal boundary layer, line vortex
# =====================================================================================================================
def test_similarity_reduce_V2_all_cases():  # V2 (C10): reduced ODEs, exponents and closed forms
    d = ch08.similarity_reduce_sympy("stokes1_delta")
    assert d["delta_matches"] and d["residual"] == 0
    d = ch08.similarity_reduce_sympy("vortex_sheet")
    assert d["n"] == sp.Rational(1, 2) and d["residual"] == 0 and d["omega_matches"] and d["pde_residual"] == 0
    d = ch08.similarity_reduce_sympy("line_vortex")
    assert d["residual"] == 0 and d["pde_residual"] == 0 and d["spinup_residual"] == 0 and d["bcs"] == (0, 1)
    d = ch08.similarity_reduce_sympy("spreading")
    assert (d["n"], d["m"]) == (sp.Rational(1, 5), sp.Rational(1, 5)) and d["t_free"] and d["residual"] == 0
    with pytest.raises(ValueError):
        ch08.similarity_reduce_sympy("nope")


def test_similarity_example_8_4_V2_derivation():  # V2 derivation D21 ★★: brackets proportional ⇒ δ = √(2C1νt)
    y, t, nu, U, C1 = sp.symbols("y t nu U C_1", positive=True)
    e = sp.Symbol("eta", positive=True)
    d = sp.Function("delta")(t)
    F = sp.Function("F")
    u = U * F(y / d)
    dut = sp.diff(u, t).subs(y, e * d).doit()  # step 2: −Uη(δ′/δ)F′
    assert z0(dut + U * e * sp.diff(d, t) / d * sp.diff(F(e), e))
    duyy = sp.diff(u, y, 2).subs(y, e * d).doit()  # step 3: (U/δ²)F″
    assert z0(duyy * d ** 2 / U - sp.diff(F(e), e, 2))
    # step 5: dividing (8.20) by U gives the brackets [δ′/δ] and [ν/δ²]
    red = sp.expand((dut - nu * duyy) / U)
    assert z0(red.coeff(sp.diff(F(e), e)) + e * sp.diff(d, t) / d) and z0(red.coeff(sp.diff(F(e), e, 2)) + nu / d ** 2)
    ds = sp.dsolve(sp.Eq(d * sp.diff(d, t), C1 * nu), d)  # step 7
    ds = ds if isinstance(ds, list) else [ds]
    found = False
    for s_ in ds:
        cc = [c for c in s_.rhs.free_symbols if c.name == "C1"]
        e0 = s_.rhs.subs(cc[0], 0) if cc else s_.rhs
        if z0(e0 - sp.sqrt(2 * C1 * nu * t)):
            found = True
    assert found  # step 8: δ(0) = 0
    assert z0(sp.sqrt(2 * C1 * nu * t).subs(C1, sp.Rational(1, 2)) - sp.sqrt(nu * t))
    assert ch08.similarity_reduce_sympy("stokes1_delta")["delta_matches"]


def test_vortex_sheet_V2_derivation():  # V2 derivation D22 ★★★: Part F steps 1–14
    y, t, nu, U, A, D = sp.symbols("y t nu U A D", positive=True)
    n = sp.Symbol("n", real=True)
    e = sp.Symbol("eta", real=True)
    F = sp.Function("F")
    # step 1: ω = −∂u/∂y obeys the same diffusion equation (derivatives commute)
    uf = sp.Function("u")(y, t)
    assert z0(sp.diff(-sp.diff(uf, y), t) - (-sp.diff(sp.diff(uf, t), y)))
    # steps 2–6: ansatz with δ = √(νt) ⇒ −nF − ½ηF′ = F″
    delta = sp.sqrt(nu * t)
    w = A * t ** (-n) * F(y / delta)
    red = sp.simplify(((sp.diff(w, t) - nu * sp.diff(w, y, 2)) * t ** (n + 1) / A).subs(y, e * delta).doit())
    target = -n * F(e) - e / 2 * sp.diff(F(e), e) - sp.diff(F(e), e, 2)
    assert z0(red - target)
    # steps 7–9: −∫ω dy = −A t^{−n} δ ∫F dη must be constant ⇒ t^{−n+1/2} constant ⇒ n = ½
    power = sp.powsimp(t ** (-n) * sp.sqrt(t), force=True)
    assert z0(sp.log(power).expand(force=True).coeff(sp.log(t)) - (sp.Rational(1, 2) - n))
    assert sp.solve(sp.Rational(1, 2) - n, n) == [sp.Rational(1, 2)]
    # step 10: F + ηF′ = (ηF)′
    assert z0(F(e) + e * sp.diff(F(e), e) - sp.diff(e * F(e), e))
    # steps 11–12: F′ + ηF/2 = 0 ⇒ F = De^{−η²/4}, which solves the n = ½ ODE; ηF → 0
    Fs = D * sp.exp(-e ** 2 / 4)
    assert z0(sp.diff(Fs, e) + e / 2 * Fs)
    assert z0(target.subs(n, sp.Rational(1, 2)).subs(F(e), Fs).doit())
    assert sp.limit(e * Fs, e, sp.oo) == 0
    # step 13: ∫e^{−η²/4} dη = 2√π ⇒ AD = −U/√(πν)
    I = sp.integrate(sp.exp(-e ** 2 / 4), (e, -sp.oo, sp.oo))
    assert z0(I - 2 * sp.sqrt(sp.pi))
    AD = sp.Symbol("AD", real=True)
    ADv = sp.solve(sp.Eq(-AD * sp.sqrt(nu) * I, 2 * U), AD)[0]
    assert z0(ADv + U / sp.sqrt(sp.pi * nu))
    # step 14: ω and u; −∫ω dy = 2U for every t; ω = −∂u/∂y; u(0) = 0
    yy = sp.Symbol("y", real=True)
    om = -U / sp.sqrt(sp.pi * nu * t) * sp.exp(-yy ** 2 / (4 * nu * t))
    u = U * sp.erf(yy / (2 * sp.sqrt(nu * t)))
    assert z0(om + sp.diff(u, yy)) and u.subs(yy, 0) == 0
    assert z0(-sp.integrate(om, (yy, -sp.oo, sp.oo)) - 2 * U)
    assert z0(sp.diff(om, t) - nu * sp.diff(om, yy, 2))
    # check line: ±0.95U at η = ±2.772 (not the printed 2.76), width 5.544√(νt); U = 1 cm/s, t = 1 s: peak ω = −5.64 s⁻¹
    assert 2 * erfinv(0.95) == pytest.approx(2.7718, rel=1e-4)
    assert float(erf(2.76 / 2)) < 0.9495  # the printed 2.76 is not the 95 % point (R11)
    ww = ch08.vortex_sheet_diffusion(0.0, 1.0, 0.01, 1e-6)[1]
    assert ww == pytest.approx(-5.642, rel=1e-3)
    assert ch08.transition_width(1.0, 1e-6) == pytest.approx(5.544e-3, rel=1e-3)


def test_vortex_sheet_V4_conserved_jump_and_ch05_parity():  # V4 (N59): −∫ω dy = 2U for all t; ch05 with γ = −2U
    U, nu = 0.3, 1e-6
    for t in (0.5, 50.0, 5e3):
        s = np.sqrt(nu * t)
        val = quad(lambda yy: float(ch08.vortex_sheet_diffusion(yy, t, U, nu)[1]), -40 * s, 40 * s, epsrel=1e-12, limit=200)[0]
        assert -val == pytest.approx(2 * U, rel=1e-10)
    y = np.linspace(-0.01, 0.01, 81)
    u8, w8 = ch08.vortex_sheet_diffusion(y, 20.0, U, nu)
    u5, w5 = ch05.diffusing_vortex_sheet(y, 20.0, -2 * U, nu)
    assert maxrel(u8, u5) < 1e-14 and maxrel(w8, w5) < 1e-14
    assert maxrel(u8, ch05.diffusing_vortex_sheet(y, 20.0, 2 * U, nu)[0]) > 1.0  # wrong variant γ = +2U
    dy = 1e-8
    fd = -(np.asarray(ch08.vortex_sheet_diffusion(y + dy, 20.0, U, nu)[0]) - np.asarray(ch08.vortex_sheet_diffusion(y - dy, 20.0, U, nu)[0])) / (2 * dy)
    assert maxrel(w8, fd) < 1e-6
    wdt = ch08.transition_width(20.0, nu)
    assert float(ch08.vortex_sheet_diffusion(wdt / 2, 20.0, U, nu)[0]) == pytest.approx(0.95 * U, rel=1e-12)


def test_temporal_bl_V1_wall_stress_and_cf():  # V1 (N60): τ_w = μ∂u/∂y(0) by FD; C_f = (2/√π)Re_x^{−1/2}
    U, nu, rho = 0.2, 1.5e-5, 1.2
    for t in (0.1, 3.0):
        d = ch08.temporal_bl_wall_stress(t, U, nu, rho)
        dy = 1e-7 * np.sqrt(nu * t)
        fd = (float(ch08.vortex_sheet_diffusion(dy, t, U, nu)[0]) - float(ch08.vortex_sheet_diffusion(-dy, t, U, nu)[0])) / (2 * dy)
        assert d["tau_w"] == pytest.approx(rho * nu * fd, rel=1e-7)
        assert d["Cf"] == pytest.approx(d["tau_w"] / (0.5 * rho * U ** 2), rel=1e-14)
        assert d["Cf"] == pytest.approx(1.1284 / np.sqrt(d["Re_x"]), rel=1e-4)
        assert d["Cf_coefficient"] == pytest.approx(2 / np.sqrt(np.pi))
    # the upper half is Stokes' first problem seen from the plate (Galilean shift): U − u = stokes_first(2U)
    y = np.linspace(0, 0.01, 21)
    assert maxrel(U - np.asarray(ch08.vortex_sheet_diffusion(y, 2.0, U, nu)[0]), ch08.stokes_first_problem(y, 2.0, U, nu)) < 1e-12


def test_line_vortex_V1_parity_circulation_and_axis():  # V1 (N61, N62): Lamb–Oseen = Gaussian vortex σ² = 4νt
    G, nu = 0.05, 1e-6
    r = np.linspace(0, 0.02, 101)
    for t in (10.0, 1000.0):
        u = ch08.line_vortex_decay(r, t, G, nu)
        assert maxrel(u, VX.gaussian_vortex(r, G, 2 * np.sqrt(nu * t))[0]) < 1e-14
        assert maxrel(u, VX.gaussian_vortex(r, G, np.sqrt(2 * nu * t))[0]) > 0.05  # wrong variant σ² = 2νt
        X = np.stack([r[1:], np.zeros_like(r[1:])])
        ulo = NS.exact_solution("lamb_oseen", X, t, Gamma=G, nu=nu)[0][1]
        assert maxrel(np.asarray(u)[1:], ulo) < 1e-12
    assert ch08.line_vortex_decay(0.0, 10.0, G, nu) == 0.0
    small = 1e-9
    assert ch08.line_vortex_decay(small, 10.0, G, nu) == pytest.approx(G * small / (8 * np.pi * nu * 10.0), rel=1e-8)
    assert 2 * np.pi * 50.0 * ch08.line_vortex_decay(50.0, 10.0, G, nu) == pytest.approx(G, rel=1e-14)  # Γ far away
    # spin-up (Exercise 8.26): circulation at fixed r → Γ as t → ∞, → 0 as t → 0
    assert 2 * np.pi * 0.01 * ch08.line_vortex_spinup(0.01, 1e9, G, nu) == pytest.approx(G, rel=1e-4)
    assert ch08.line_vortex_spinup(0.01, 1e-3, G, nu) < 1e-100


def test_line_vortex_V2_residuals():  # V2 (N61, N62): both solve ∂u/∂t = ν∂/∂r[(1/r)∂(ru)/∂r]
    r, t, nu, G = sp.symbols("r t nu Gamma", positive=True)
    op = lambda q: sp.simplify(sp.diff(q, t) - nu * sp.diff(sp.diff(r * q, r) / r, r))  # noqa: E731
    dec = G / (2 * sp.pi * r) * (1 - sp.exp(-r ** 2 / (4 * nu * t)))
    spin = G / (2 * sp.pi * r) * sp.exp(-r ** 2 / (4 * nu * t))
    assert op(dec) == 0 and op(spin) == 0
    assert op(G / (2 * sp.pi * r) * (1 - sp.exp(-r ** 2 / (2 * nu * t)))) != 0  # wrong σ²


def test_similarity_example_8_7_V2_derivation():  # V2 derivation D23 ★★: n = m = 1/5 from the PDE and the volume
    x, t, A, D, beta = sp.symbols("x t A D beta", positive=True)
    n, m = sp.symbols("n m", real=True)
    F = sp.Function("F")
    e = sp.Symbol("eta", positive=True)
    h = A * t ** (-n) * F(x / (D * t ** m))
    lhs = sp.diff(h, t)
    rhs = beta * sp.diff(h ** 3 * sp.diff(h, x), x)
    # step 4: (h³h_x)_x = A⁴t^{−4n}/δ² (3F²F′² + F³F″)
    flux_x = sp.simplify((rhs / beta).subs(x, e * D * t ** m).doit())
    Fe = F(e)
    target = A ** 4 * t ** (-4 * n) / (D * t ** m) ** 2 * (3 * Fe ** 2 * sp.diff(Fe, e) ** 2 + Fe ** 3 * sp.diff(Fe, e, 2))
    assert z0(flux_x - target)
    # steps 7–9: t^{−n−1} ∝ t^{−4n−2m}, volume t^{−n+m} constant ⇒ n = m = 1/5
    sol = sp.solve([sp.Eq(-n - 1, -4 * n - 2 * m), sp.Eq(-n + m, 0)], [n, m], dict=True)[0]
    assert sol == {n: sp.Rational(1, 5), m: sp.Rational(1, 5)}
    red = sp.simplify(((lhs - rhs) * t ** (sol[n] + 1) / A).subs({n: sol[n], m: sol[m]}).subs(x, e * D * t ** sp.Rational(1, 5)).doit())
    assert not red.has(t)  # the reduced equation is an ODE in η
    # wrong exponents leave t in the equation
    red_bad = sp.simplify(((lhs - rhs) * t ** (sp.Rational(1, 4) + 1) / A).subs({n: sp.Rational(1, 4), m: sp.Rational(1, 4)}).subs(x, e * D * t ** sp.Rational(1, 4)).doit())
    assert red_bad.has(t)


def test_similarity_collapse_V1_right_exponents_only():  # V1 (C10, E6): spread ≈ 0 at the right (n, m), > 1e-2 elsewhere
    for case, n, m in (("stokes1", 0.0, 0.5), ("vortex_sheet", 0.5, 0.5), ("line_vortex", 1.0, 0.5), ("spreading", 0.2, 0.2)):
        assert ch08.similarity_collapse_error(case, n, m) < 1e-12
        assert ch08.similarity_collapse_error(case, n + 0.1, m) > 1e-2
        assert ch08.similarity_collapse_error(case, n, m + 0.1) > 1e-2
    assert ch08.similarity_collapse_error("vortex_sheet", 0.25, 0.5) > 1e-2  # analysis wrong variant n = ¼


# =====================================================================================================================
# C11 — Stokes' second problem (8.33)–(8.38)
# =====================================================================================================================
def test_stokes_second_V2_engine_and_D24_derivation():  # V2 derivation D24 ★★: (8.35) → (8.36) → (8.37) → (8.38)
    d = ch08.stokes_second_sympy()
    assert d["residual"] == 0 and d["roots_match"]
    y, t, w, nu, U = sp.symbols("y t omega nu U", positive=True)
    f = sp.Function("f")
    ut = sp.exp(sp.I * w * t) * f(y)
    ode = sp.simplify((sp.diff(ut, t) - nu * sp.diff(ut, y, 2)) / sp.exp(sp.I * w * t))  # steps 1–3
    assert z0(ode - (sp.I * w * f(y) - nu * sp.diff(f(y), y, 2)))
    assert z0(sp.expand(((1 + sp.I) / sp.sqrt(2)) ** 2) - sp.I)  # step 5: √i = (1 + i)/√2
    k = (1 + sp.I) * sp.sqrt(w / (2 * nu))
    assert z0(sp.expand(nu * k ** 2) - sp.I * w)
    assert sp.re(sp.expand_complex(k)).is_positive  # step 7: e^{+ky} grows ⇒ B = 0
    uu = sp.re(sp.expand_complex(U * sp.exp(sp.I * w * t) * sp.exp(-k * y)))  # steps 8–10
    book = U * sp.exp(-y * sp.sqrt(w / (2 * nu))) * sp.cos(w * t - y * sp.sqrt(w / (2 * nu)))
    assert z0(sp.expand_trig(uu - book))
    assert z0(sp.diff(book, t) - nu * sp.diff(book, y, 2)) and z0(book.subs(y, 0) - U * sp.cos(w * t))
    # wrong variant: the growing root (+k) is unbounded
    assert sp.limit(sp.exp(y * sp.sqrt(w / (2 * nu))), y, sp.oo) == sp.oo


def test_stokes_second_V1_wall_envelope_phase_speed_and_D25():  # V1 (C11, N72, D25): zero-crossing tracking = √(2νω)
    nu, om, U = 1e-6, 2 * np.pi, 0.02
    t = np.linspace(0, 2 * np.pi / om, 50)
    assert maxrel(ch08.stokes_second_problem(0.0, t, U, om, nu), U * np.cos(om * t)) < 1e-15
    s = ch08.stokes_layer(nu, om)
    y = np.linspace(0, 6 * s["delta_e"], 200)
    tt = np.linspace(0, 2 * np.pi / om, 2001)
    amp = np.max(np.abs(np.asarray(ch08.stokes_second_problem(y[:, None], tt[None, :], U, om, nu))), axis=1)
    assert maxrel(amp, ch08.stokes_layer_envelope(y, U, om, nu)) < 2e-5  # the envelope Ue^{−y/δ_e}
    # track the zero crossing ωt − y/δ_e = π/2 in time → dy/dt = ωδ_e = √(2νω)
    ycross = []
    for tc in (0.3, 0.35, 0.4):
        yy = np.linspace(0, 3 * s["delta_e"], 3001)
        uu = np.asarray(ch08.stokes_second_problem(yy, tc, 1.0, om, nu))
        i = np.nonzero(np.sign(uu[:-1]) != np.sign(uu[1:]))[0][0]
        ycross.append(brentq(lambda q: float(ch08.stokes_second_problem(q, tc, 1.0, om, nu)), yy[i], yy[i + 1], xtol=1e-15))
    speed = np.polyfit([0.3, 0.35, 0.4], ycross, 1)[0]
    assert speed == pytest.approx(np.sqrt(2 * nu * om), rel=1e-6) and s["phase_speed"] == pytest.approx(np.sqrt(2 * nu * om))
    # D25 numbers: e^{−2√2} = 0.0591 at 4√(ν/ω); δ_e = 0.56 mm; 4√(ν/ω) = 1.6 mm; crest 3.5 mm/s; ratio 2√2
    assert s["amplitude_at_delta_book"] == pytest.approx(0.05911, rel=1e-3)
    assert float(ch08.stokes_layer_envelope(s["delta_book"], 1.0, om, nu)) == pytest.approx(s["amplitude_at_delta_book"], rel=1e-14)
    assert (s["delta_e"], s["delta_book"], s["phase_speed"]) == (pytest.approx(5.64e-4, rel=1e-3), pytest.approx(1.596e-3, rel=1e-3), pytest.approx(3.545e-3, rel=1e-3))
    assert s["ratio_book_to_e"] == pytest.approx(2 * np.sqrt(2)) and s["wavelength"] == pytest.approx(2 * np.pi * s["delta_e"])
    st = ch08.stokes_layer_state(nu, om, 2 * s["delta_e"], U)
    assert st["amplitude"] == pytest.approx(U * np.exp(-2)) and st["phase_lag"] == pytest.approx(2.0)
    assert st["time_lag"] == pytest.approx(2.0 / om) and st["period"] == pytest.approx(1.0)
    # N73: three groups u/U, ωt, y√(ω/ν) — the Π count (not self-similar)
    g = ch01.pi_groups({"u": "m/s", "U": "m/s", "y": "m", "t": "s", "nu": "m^2/s", "omega": "1/s"}, solution="u")
    assert len(g) == 4  # 6 variables, 2 dimensions → 4 groups; with U fixed by linearity: 3 (the book's)


@needs_ref
def test_stokes_second_V1_form_wikipedia():  # V1 form cross-check (reference/ch08): the published u(y, t) and δ = √(2ν/ω)
    ref = ref_json()["stokes_second_problem"]
    y, t, omega, nu, U = sp.symbols("y t omega nu U", positive=True)
    pub = sp.sympify(ref["u_form"], locals=dict(y=y, t=t, omega=omega, nu=nu, U=U))
    assert z0(pub - ch08.stokes_second_sympy()["u_book"].subs({sp.Symbol("omega", positive=True): omega}))
    assert ch08.stokes_layer(1e-6, 3.0)["delta_e"] == pytest.approx(float(sp.sympify(ref["penetration_depth"]).subs({"nu": 1e-6, "omega": 3.0})))


def test_stokes_second_V3_crank_nicolson_orders_and_transients():  # V3 (C11): CN twin order 2 in y and t; from rest → (8.38)
    nu, om, U = 1e-6, 2 * np.pi, 1.0
    de = np.sqrt(2 * nu / om)
    Ly = 30 * de
    Tp = 2 * np.pi / om
    wall = lambda t: U * np.cos(om * t)  # noqa: E731
    hs, es = [], []
    for ny in (61, 121, 241, 481):
        y = np.linspace(0, Ly, ny)
        u = DIF.crank_nicolson_1d(ch08.stokes_second_problem(y, 0.0, U, om, nu), y, Tp / 4000, 4000, nu, wall, 0.0, startup_be=0)
        es.append(np.max(np.abs(u - ch08.stokes_second_problem(y, Tp, U, om, nu))))
        hs.append(y[1] - y[0])
    assert abs(observed_order(hs, es) - 2.0) < ORDER_TOL
    y = np.linspace(0, Ly, 1601)
    u0 = ch08.stokes_second_problem(y, 0.0, U, om, nu)
    uref = DIF.crank_nicolson_1d(u0, y, Tp / 8000, 8000, nu, wall, 0.0, startup_be=0)
    dts, es = [], []
    for n in (25, 50, 100, 200):
        es.append(np.max(np.abs(DIF.crank_nicolson_1d(u0, y, Tp / n, n, nu, wall, 0.0, startup_be=0) - uref)))
        dts.append(Tp / n)
    assert abs(observed_order(dts, es) - 2.0) < ORDER_TOL
    # from rest: after 10 periods the transient is gone to 2e-3 U within 6δ_e (design Part F D24 check)
    y = np.linspace(0, 20 * de, 401)
    tt, UU = DIF.crank_nicolson_1d(np.zeros_like(y), y, 10 * Tp / 2000, 2000, nu, wall, 0.0, return_all=True, save_every=200)
    m = y < 6 * de
    err = [np.max(np.abs(UU[k] - ch08.stokes_second_problem(y, tt[k], U, om, nu))[m]) for k in range(len(tt))]
    assert err[-1] < 2e-3 * U and err[-1] < err[1]


# =====================================================================================================================
# C12 — the Stokes equations (8.39)–(8.43)
# =====================================================================================================================
def test_low_re_scaling_V2_engine_and_D26_derivation():  # V2 derivation D26 ★★: dynamic vs viscous pressure scale
    d = ch08.low_re_scaling_sympy()
    Re = sp.Symbol("Re", positive=True)
    assert d["dynamic"]["inertia"] == 1 and d["dynamic"]["pressure"] == 1 and z0(d["dynamic"]["viscous"] - 1 / Re)
    assert d["dynamic_limit"]["pressure"] == 0  # (8.41) → the wrong equation 0 = μ∇²u
    assert d["viscous_times_Re"]["pressure"] == 1 and d["viscous_limit"] == {"advective": 0, "pressure": 1, "viscous": 1, "inertia": 0}
    # own derivation: size each term of (8.39) and divide
    rho, U, L, mu = sp.symbols("rho U L mu", positive=True)
    Re_s = rho * U * L / mu
    inertia, visc = rho * U ** 2 / L, mu * U / L ** 2
    for P, want in ((rho * U ** 2, (1, 1, 1 / Re_s)), (mu * U / L, (1, 1 / Re_s, 1 / Re_s))):
        pressure = P / L
        got = (inertia / inertia, pressure / inertia, visc / inertia)
        assert all(z0(g_ - w_) for g_, w_ in zip(got, want))
    # (8.42): dividing by μU/L² gives (Re, 1, 1); Re → 0 leaves ∇p = μ∇²u
    got = (inertia / visc, (mu * U / L) / L / visc, visc / visc)
    assert z0(got[0] - Re_s) and got[1] == 1 and got[2] == 1
    # droplet numbers (D26 check): μU/a = 0.022 Pa vs ρU² = 1.7e-4 Pa
    st = ch08.settling_state(10e-6, 1000.0, 1.2, 1.8e-5)
    assert 1.8e-5 * st["U_t"] / 10e-6 == pytest.approx(0.0218, rel=2e-2) and 1.2 * st["U_t"] ** 2 == pytest.approx(1.75e-4, rel=2e-2)


def test_stokes_residual_V1_sphere_field_and_ideal_flow_fails():  # V1 (C12): ∇p − μ∇²u ≈ 0 for Stokes; ≠ 0 for ideal flow
    U, a, mu = 1.0, 1.0, 1.0
    uf = lambda X: np.stack(CRP.stokes_sphere_velocity_xyz(X[0], X[1], X[2], U, a))  # noqa: E731
    pf = lambda X: CRP.stokes_sphere_pressure(np.sqrt(X[0] ** 2 + X[1] ** 2 + X[2] ** 2), np.arctan2(np.hypot(X[1], X[2]), X[0]), U, a, mu)  # noqa: E731
    P = RNG.normal(size=(3, 50))
    P = P / np.linalg.norm(P, axis=0) * RNG.uniform(1.5, 5.0, 50)
    scale = 3 * mu * U * a / 1.5 ** 3  # |∇p| ~ 3μUa/r³
    res = CRP.stokes_residual(uf, pf, P, mu, 1e-4)
    assert np.max(np.abs(res)) < 1e-6 * scale
    errs = [np.max(np.abs(CRP.stokes_residual(uf, pf, P, mu, h))) for h in (1e-2, 5e-3, 2.5e-3)]
    assert abs(observed_order([1e-2, 5e-3, 2.5e-3], errs) - 2.0) < ORDER_TOL  # pure O(h²) truncation
    sph = PF.sphere(U, a)

    def uid(X):
        rr = np.sqrt(X[0] ** 2 + X[1] ** 2 + X[2] ** 2)
        th = np.arctan2(np.hypot(X[1], X[2]), X[0])
        ur, ut = sph.velocity_spherical(rr, th)
        c, s_ = np.cos(th), np.sin(th)
        rad = np.asarray(ur) * s_ + np.asarray(ut) * c
        rc = np.hypot(X[1], X[2])
        return np.stack([np.asarray(ur) * c - np.asarray(ut) * s_, rad * X[1] / rc, rad * X[2] / rc])

    # the ideal-flow sphere has ∇²u = 0, so its viscous Stokes pressure balance fails with the Stokes pressure
    res_id = CRP.stokes_residual(uid, pf, P, mu, 1e-4)
    assert np.max(np.abs(res_id)) > 0.1 * scale


# =====================================================================================================================
# C13 — Stokes' stream function (8.44)–(8.49)
# =====================================================================================================================
def test_stokes_sphere_sympy_V2_engine():  # V2 (C13, C14): E⁴ψ = 0, roots, constants, p, divergence, drag parts
    d = ch08.stokes_sphere_sympy()
    r, th, a, U, mu = d["symbols"]
    assert d["roots"] == [-1, 1, 2, 4]
    assert d["constants"] == {sp.Symbol("A"): 0, sp.Symbol("B"): U / 2, sp.Symbol("C"): -3 * U * a / 4, sp.Symbol("D"): U * a ** 3 / 4}
    assert d["E4psi"] == 0 and d["divergence"] == 0 and d["curlcurl_identity"] == 0
    assert d["p_check_r"] == 0 and d["p_check_theta"] == 0 and d["sigma_rr_viscous_a"] == 0
    assert z0(d["p"] + 3 * mu * a * U * sp.cos(th) / (2 * r ** 2))
    assert z0(d["D_pressure"] - 2 * sp.pi * mu * a * U) and z0(d["D_friction"] - 4 * sp.pi * mu * a * U)
    assert d["biharmonic_psi"] != 0  # wrong variant (R17): the scalar biharmonic of ψ does not vanish
    assert z0(CRP.E4_residual(d["psi"], r, th))


def test_stokes_stream_function_V2_derivation_D28():  # V2 derivation D28 ★★★: ∇×∇×(A e_φ) = −E²(r sinθ A)/(r sinθ) e_φ
    r, th, ph = CL.coordinates("spherical")
    A = sp.Function("A")(r, th)
    cc = CL.curl(CL.curl([0, 0, A], "spherical"), "spherical")  # step 7–9, built with core.curvilinear
    Psi = r * sp.sin(th) * A
    assert z0(cc[0]) and z0(cc[1])
    assert z0(cc[2] + CRP.E2(Psi, r, th) / (r * sp.sin(th)))  # step 9
    c1 = CL.curl([0, 0, A], "spherical")
    assert z0(c1[0] - sp.diff(Psi, th) / (r ** 2 * sp.sin(th))) and z0(c1[1] + sp.diff(Psi, r) / (r * sp.sin(th)))  # step 8
    # steps 2–6: for any ψ, the (6.83) velocity has ω_φ = −E²ψ/(r sinθ)
    psi = sp.Function("psi")(r, th)
    ur = sp.diff(psi, th) / (r ** 2 * sp.sin(th))
    ut = -sp.diff(psi, r) / (r * sp.sin(th))
    om = CL.curl([ur, ut, 0], "spherical")
    assert z0(om[0]) and z0(om[1])
    assert z0(sp.diff(r * ut, r) + sp.diff(psi, r, 2) / sp.sin(th))  # step 3
    assert z0(om[2] + CRP.E2(psi, r, th) / (r * sp.sin(th)))  # step 6
    # steps 10–12: −∇×∇×ω = −E²(E²ψ)/(r sinθ) e_φ; for (8.48) E⁴ψ = 0 and hence −∇×∇×ω = 0
    a, U = sp.symbols("a U", positive=True)
    psi48 = U * r ** 2 * sp.sin(th) ** 2 * (sp.Rational(1, 2) - 3 * a / (4 * r) + a ** 3 / (4 * r ** 3))
    om48 = -CRP.E2(psi48, r, th) / (r * sp.sin(th))
    ccc = CL.curl(CL.curl([0, 0, om48], "spherical"), "spherical")
    assert all(z0(c) for c in ccc)
    assert z0(CRP.E2(CRP.E2(psi48, r, th), r, th))
    lap = lambda q: sp.diff(r ** 2 * sp.diff(q, r), r) / r ** 2 + sp.diff(sp.sin(th) * sp.diff(q, th), th) / (r ** 2 * sp.sin(th))  # noqa: E731
    assert not z0(lap(lap(psi48)))  # (8.44) is (E²)², not ∇⁴


def test_stokes_stream_function_V2_derivation_D29_D27():  # V2 derivations D29 ★★ (separable ψ) and D27 ★★ (curl kills p)
    r, th = sp.symbols("r theta", positive=True)
    f = sp.Function("f")
    E2f = sp.simplify(CRP.E2(f(r) * sp.sin(th) ** 2, r, th) / sp.sin(th) ** 2)
    assert z0(E2f - (sp.diff(f(r), r, 2) - 2 * f(r) / r ** 2))  # step 3
    g = sp.diff(f(r), r, 2) - 2 * f(r) / r ** 2
    ode = sp.expand(sp.diff(g, r, 2) - 2 * g / r ** 2)  # step 4
    target = sp.diff(f(r), r, 4) - 4 * sp.diff(f(r), r, 2) / r ** 2 + 8 * sp.diff(f(r), r) / r ** 3 - 8 * f(r) / r ** 4
    assert z0(ode - target)
    lam = sp.Symbol("lambda")
    char = sp.expand(sp.simplify(target.subs(f(r), r ** lam).doit() / r ** (lam - 4)))
    assert z0(char - sp.expand(lam * (lam - 1) * (lam - 2) * (lam - 3) - 4 * lam * (lam - 1) + 8 * lam - 8))  # step 5
    assert z0(sp.factor(char) - (lam - 4) * (lam - 2) * (lam - 1) * (lam + 1))  # step 6
    a, U = sp.symbols("a U", positive=True)
    C, D = sp.symbols("C D")
    fg = U / 2 * r ** 2 + C * r + D / r  # steps 7–8: A = 0, B = U/2
    sol = sp.solve([fg.subs(r, a), sp.diff(fg, r).subs(r, a)], [C, D], dict=True)[0]  # steps 9–10
    assert sol == {C: -3 * U * a / 4, D: U * a ** 3 / 4}
    # D27: in Cartesian components the curl kills any pressure and commutes with the Laplacian
    x, y, z = sp.symbols("x y z", real=True)
    p = sp.Function("p")(x, y, z)
    assert all(z0(c) for c in CL.curl(CL.gradient(p, "cartesian"), "cartesian"))
    uu = [sp.Function(n_)(x, y, z) for n_ in ("u", "v", "w")]
    lapv = [sp.diff(c, x, 2) + sp.diff(c, y, 2) + sp.diff(c, z, 2) for c in uu]
    curl_lap = CL.curl(lapv, "cartesian", simplify=False)
    lap_curl = [sp.diff(c, x, 2) + sp.diff(c, y, 2) + sp.diff(c, z, 2) for c in CL.curl(uu, "cartesian", simplify=False)]
    assert all(z0(p_ - q_) for p_, q_ in zip(curl_lap, lap_curl))


def test_stokes_sphere_velocity_V1_walls_far_field_divergence():  # V1 (N87): no slip, → U e_x, ∇·u = 0, (6.83) parity
    U, a = 0.7, 2.0
    th = np.linspace(0.01, np.pi - 0.01, 40)
    ur, ut = ch08.stokes_sphere_velocity(a, th, U, a)
    assert np.max(np.abs(ur)) < 1e-15 * U and np.max(np.abs(ut)) < 1e-15 * U
    ur, ut = ch08.stokes_sphere_velocity(1e6 * a, th, U, a)
    assert np.max(np.abs(np.asarray(ur) - U * np.cos(th))) < 1e-5 * U and np.max(np.abs(np.asarray(ut) + U * np.sin(th))) < 1e-5 * U
    # velocity from the stream function by central differences (6.83)
    r = np.linspace(1.2 * a, 6 * a, 30)
    R, T = np.meshgrid(r, th)
    d = 1e-6
    psi = lambda rr, tt: np.asarray(ch08.stokes_sphere_streamfunction(rr, tt, U, a))  # noqa: E731
    urf = (psi(R, T + d) - psi(R, T - d)) / (2 * d) / (R ** 2 * np.sin(T))
    utf = -(psi(R + d, T) - psi(R - d, T)) / (2 * d) / (R * np.sin(T))
    ur, ut = ch08.stokes_sphere_velocity(R, T, U, a)
    assert maxrel(ur, urf) < 1e-7 and maxrel(ut, utf) < 1e-7
    # ∇·u = 0 by finite differences of the Cartesian field
    P = np.stack([RNG.uniform(-5, 5, 40), RNG.uniform(-5, 5, 40), RNG.uniform(-5, 5, 40)])
    P = P[:, np.linalg.norm(P, axis=0) > 1.5 * a]
    h = 1e-5
    div = sum((np.asarray(ch08.stokes_sphere_velocity_xyz(*(P + h * np.eye(3)[i][:, None]), U, a)[i])
               - np.asarray(ch08.stokes_sphere_velocity_xyz(*(P - h * np.eye(3)[i][:, None]), U, a)[i])) / (2 * h) for i in range(3))
    assert np.max(np.abs(div)) < 1e-8
    assert np.isnan(ch08.stokes_sphere_streamfunction(0.5 * a, 1.0, U, a)) and np.isnan(ch08.stokes_sphere_velocity(0.5 * a, 1.0, U, a)[0])
    # wrong variant: flipping the sign of u_θ breaks no slip only if … it breaks (6.83) consistency
    assert maxrel(-np.asarray(ut), utf) > 1.0


def test_stokes_sphere_V7_frames_and_fore_aft_symmetry():  # V7 (N93): fluid frame = body − stream; |u|(θ) = |u|(π − θ)
    U, a = 1.0, 1.0
    r = np.linspace(1.1, 8, 25)[:, None]
    th = np.linspace(0.05, np.pi / 2, 20)[None, :]
    urb, utb = ch08.stokes_sphere_velocity(r, th, U, a)
    urf, utf = ch08.stokes_sphere_velocity(r, th, U, a, "fluid")
    assert maxrel(np.asarray(urb) - np.asarray(urf), U * np.cos(th) * np.ones_like(r)) < 1e-14
    assert maxrel(np.asarray(utb) - np.asarray(utf), -U * np.sin(th) * np.ones_like(r)) < 1e-14
    s1 = np.hypot(*ch08.stokes_sphere_velocity(r, th, U, a, "fluid"))
    s2 = np.hypot(*ch08.stokes_sphere_velocity(r, np.pi - th, U, a, "fluid"))
    assert maxrel(s1, s2) < 1e-13  # no wake: reversibility of creeping flow
    assert maxrel(ch08.stokes_sphere_streamfunction(r, th, U, a) - np.asarray(ch08.stokes_sphere_streamfunction(r, th, U, a, "fluid")),
                  0.5 * U * r ** 2 * np.sin(th) ** 2) < 1e-14
    # reversing U reverses u and p (linearity)
    assert maxrel(ch08.stokes_sphere_velocity(r, th, -U, a)[0], -np.asarray(urb)) < 1e-15
    assert maxrel(ch08.stokes_sphere_pressure(r, th, -U, a), -np.asarray(ch08.stokes_sphere_pressure(r, th, U, a))) < 1e-15
    # side line: Stokes 0 at the wall, ideal sphere 1.5U; the ideal speed matches core.potential.sphere
    assert ch08.side_line_speed(a, U, a, "stokes") == pytest.approx(0.0, abs=1e-15)
    assert ch08.side_line_speed(a, U, a, "ideal") == pytest.approx(1.5 * U)
    ur_i, ut_i = PF.sphere(U, a).velocity_spherical(np.array([1.3, 2.0]), np.pi / 2)
    assert maxrel(ch08.side_line_speed(np.array([1.3, 2.0]), U, a, "ideal"), np.hypot(ur_i, ut_i)) < 1e-12
    assert ch08.side_line_speed(3.0, U, a, "oseen", Re=1e-9) == pytest.approx(float(ch08.side_line_speed(3.0, U, a, "stokes")), rel=1e-8)


# =====================================================================================================================
# C14 — Stokes drag (8.50)–(8.52), settling, Millikan
# =====================================================================================================================
def test_stokes_pressure_V2_derivation_D30():  # V2 derivation D30 ★★★: ∇p = μ∇²u = −μ∇×ω integrated and checked
    r, th, ph = CL.coordinates("spherical")
    a, U, mu = sp.symbols("a U mu", positive=True)
    f = U / 2 * r ** 2 - 3 * U * a / 4 * r + U * a ** 3 / (4 * r)
    psi = f * sp.sin(th) ** 2
    ur = sp.simplify(sp.diff(psi, th) / (r ** 2 * sp.sin(th)))
    ut = sp.simplify(-sp.diff(psi, r) / (r * sp.sin(th)))
    assert z0(sp.diff(f, r, 2) - 2 * f / r ** 2 - 3 * U * a / (2 * r))  # step 2: E²ψ = (3Ua/2r) sin²θ
    om = -CRP.E2(psi, r, th) / (r * sp.sin(th))
    assert z0(om + 3 * U * a / (2 * r ** 2) * sp.sin(th))  # step 3
    Psi = r * sp.sin(th) * om
    assert z0(Psi + 3 * U * a / (2 * r) * sp.sin(th) ** 2)  # step 4
    cw = CL.curl([0, 0, om], "spherical")
    assert z0(cw[0] + 3 * U * a * sp.cos(th) / r ** 3)  # step 5
    assert z0(cw[1] + 3 * U * a * sp.sin(th) / (2 * r ** 3))  # step 6
    # step 1: μ∇²u (core.curvilinear vector Laplacian, independent of the engine) equals −μ∇×ω
    lap = CL.vector_laplacian([ur, ut, 0], "spherical")
    assert z0(lap[0] + cw[0]) and z0(lap[1] + cw[1]) and z0(lap[2])
    gth = sp.Function("g")(th)
    p = sp.integrate(-mu * cw[0], r) + gth  # steps 7–8
    assert z0(p - (-3 * mu * U * a * sp.cos(th) / (2 * r ** 2) + gth))
    gprime = sp.solve(sp.Eq(sp.diff(p, th) / r, -mu * cw[1]), sp.diff(gth, th))  # step 9
    assert gprime == [0]
    p50 = -3 * mu * a * U * sp.cos(th) / (2 * r ** 2)  # step 10 (8.50)
    assert z0(sp.diff(p50, r) - mu * lap[0]) and z0(sp.diff(p50, th) / r - mu * lap[1])
    # step 11: +3μU/2a at the front (θ = π), −3μU/2a at the rear (θ = 0) — the printed + for the minimum is wrong (R12)
    assert z0(p50.subs({r: a, th: sp.pi}) - 3 * mu * U / (2 * a)) and z0(p50.subs({r: a, th: 0}) + 3 * mu * U / (2 * a))
    # wrong variant: p = +3μaU cosθ/(2r²) fails the momentum balance
    assert not z0(sp.diff(-p50, r) - mu * lap[0])


def test_stokes_drag_V2_derivation_D31():  # V2 derivation D31 ★★★: σ_rr, σ_rθ on r = a, x-traction, ⅓ + ⅔ = 6πμaU
    r, th, ph = CL.coordinates("spherical")
    a, U, mu, pinf = sp.symbols("a U mu p_inf", positive=True)
    ur = U * sp.cos(th) * (1 - 3 * a / (2 * r) + a ** 3 / (2 * r ** 3))
    ut = -U * sp.sin(th) * (1 - 3 * a / (4 * r) - a ** 3 / (4 * r ** 3))
    p = pinf - 3 * mu * a * U * sp.cos(th) / (2 * r ** 2)
    S = CL.strain_rate([ur, ut, 0], "spherical")
    assert z0(S[0, 0].subs(r, a))  # step 4: ∂u_r/∂r = 0 on the sphere
    assert z0(sp.diff(ur, th).subs(r, a))  # step 6
    s_rr = sp.simplify((-p + 2 * mu * S[0, 0]).subs(r, a))
    s_rt = sp.simplify((2 * mu * S[0, 1]).subs(r, a))
    assert z0(s_rt + 3 * mu * U / (2 * a) * sp.sin(th))  # step 7
    assert z0(s_rr - (-pinf + 3 * mu * U / (2 * a) * sp.cos(th)))  # step 8
    tx = sp.simplify(s_rr * sp.cos(th) - s_rt * sp.sin(th))  # steps 2, 9
    assert z0(tx - (-pinf * sp.cos(th) + 3 * mu * U / (2 * a)))  # uniform apart from p∞
    dA = 2 * sp.pi * a ** 2 * sp.sin(th)  # step 10
    assert sp.integrate(-pinf * sp.cos(th) * dA, (th, 0, sp.pi)) == 0  # step 11
    assert sp.integrate(sp.cos(th) ** 2 * sp.sin(th), (th, 0, sp.pi)) == sp.Rational(2, 3)
    assert sp.integrate(sp.sin(th) ** 3, (th, 0, sp.pi)) == sp.Rational(4, 3)
    Dp = sp.integrate(3 * mu * U / (2 * a) * sp.cos(th) ** 2 * dA, (th, 0, sp.pi))  # step 12
    Df = sp.integrate(3 * mu * U / (2 * a) * sp.sin(th) ** 2 * dA, (th, 0, sp.pi))  # step 13
    assert z0(Dp - 2 * sp.pi * mu * a * U) and z0(Df - 4 * sp.pi * mu * a * U) and z0(Dp + Df - 6 * sp.pi * mu * a * U)
    # running integrals of stokes_drag_running (design 3.8)
    c = sp.Symbol("c")
    run_p = sp.integrate(3 * mu * U / (2 * a) * sp.cos(th) ** 2 * dA, (th, 0, sp.Symbol("T", positive=True)))
    assert z0(sp.simplify(run_p - sp.pi * mu * a * U * (1 - sp.cos(sp.Symbol("T", positive=True)) ** 3)))


def test_stokes_drag_V1_quadrature_parts_and_slip_sphere():  # V1 (C14): Gauss–Legendre of the tractions = 6πμaU
    mu, a, U = 1.8e-5, 10e-6, 0.012
    D = ch08.stokes_drag(mu, a, U)
    assert D == pytest.approx(6 * np.pi * mu * a * U, rel=1e-15)
    for n in (1, 2, 4, 16, 64):
        Dq = ch08.sphere_drag_quadrature(lambda th: ch08.stokes_sphere_surface_stresses(th, U, a, mu, p_inf=1e5)[2], a, n)
        assert Dq == pytest.approx(D, rel=1e-12)  # exact already at n = 1: t_x is linear in cos θ
    parts = ch08.stokes_drag(mu, a, U, parts=True)
    assert parts["pressure"] == pytest.approx(D / 3) and parts["friction"] == pytest.approx(2 * D / 3)
    Dp = ch08.sphere_drag_quadrature(lambda th: -np.asarray(ch08.stokes_sphere_pressure(a, th, U, a, mu)) * np.cos(th), a, 16)
    assert Dp == pytest.approx(parts["pressure"], rel=1e-12)
    Df = ch08.sphere_drag_quadrature(lambda th: -np.asarray(ch08.stokes_sphere_surface_stresses(th, U, a, mu)[1]) * np.sin(th), a, 16)
    assert Df == pytest.approx(parts["friction"], rel=1e-12)
    run = ch08.stokes_drag_running(np.pi, mu, a, U)
    assert run["pressure"] == pytest.approx(parts["pressure"]) and run["friction"] == pytest.approx(parts["friction"])
    th = np.linspace(0.2, 3.0, 9)
    dth = 1e-6
    dR = (np.asarray(ch08.stokes_drag_running(th + dth, mu, a, U)["total"]) - np.asarray(ch08.stokes_drag_running(th - dth, mu, a, U)["total"])) / (2 * dth)
    integrand = np.asarray(ch08.stokes_sphere_surface_stresses(th, U, a, mu)[2]) * 2 * np.pi * a ** 2 * np.sin(th)
    assert maxrel(dR, integrand) < 1e-7
    # wrong variant: a slip sphere (bubble, no shear) would have D = 4πμaU — not Stokes' rigid sphere
    assert D != pytest.approx(4 * np.pi * mu * a * U, rel=0.1)
    # surface stresses by finite differences of the velocity field
    h = 1e-6 * a
    thv = np.linspace(0.3, 2.8, 7)
    g = lambda rr: np.asarray(ch08.stokes_sphere_velocity(rr, thv, U, a)[1]) / rr  # noqa: E731
    srt = mu * a * (g(a + h) - g(a)) / h  # one-sided at the wall (u_r = 0 there, so ∂u_r/∂θ = 0)
    assert maxrel(srt, ch08.stokes_sphere_surface_stresses(thv, U, a, mu)[1]) < 1e-4


def test_stokes_pressure_V1_extremes_and_printed_minimum():  # V1 (N88): ±3μU/2a; the printed + minimum is rejected (R12)
    mu, a, U = 1e-3, 0.002, 0.01
    th = np.linspace(0, np.pi, 181)
    p = np.asarray(ch08.stokes_sphere_pressure(a, th, U, a, mu))
    assert p.max() == pytest.approx(1.5 * mu * U / a) and th[np.argmax(p)] == pytest.approx(np.pi)
    assert p.min() == pytest.approx(-1.5 * mu * U / a) and th[np.argmin(p)] == 0.0
    assert p.min() != pytest.approx(+1.5 * mu * U / a)
    srr, srt, tx = ch08.stokes_sphere_surface_stresses(th, U, a, mu)
    assert maxrel(srr, -p) < 1e-14 and np.ptp(tx) < 1e-15 * abs(tx[0]) + 1e-18  # t_x uniform over the sphere
    # droplet (D30 check): 3μU/2a = 0.033 Pa for the 10 µm droplet in air
    st = ch08.settling_state(10e-6, 1000.0, 1.2, 1.8e-5)
    assert 1.5 * 1.8e-5 * st["U_t"] / 10e-6 == pytest.approx(0.033, rel=0.02)


@needs_ref
def test_stokes_law_V1_form_and_drag_coefficient():  # V1 form (reference: Stokes' law) + V1 C_D = 24/Re parity (N92)
    ref = ref_json()["stokes_law"]
    mu, R, v, rho_p, rho_f, g = 1.8e-5, 20e-6, 0.03, 2500.0, 1.2, 9.80665
    assert ch08.stokes_drag(mu, R, v) == pytest.approx(float(sp.sympify(ref["drag_form"].split("=")[1]).subs({"mu": mu, "R": R, "v": v, "pi": np.pi})))
    vs = float(sp.sympify(ref["settling_form"].split("=")[1]).subs({"rho_p": rho_p, "rho_f": rho_f, "mu": mu, "g": g, "R": R}))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        assert ch08.terminal_velocity(R, rho_p, rho_f, mu, g) == pytest.approx(vs, rel=1e-14)
    Re = np.logspace(-3, 0, 20)
    assert maxrel(ch08.stokes_drag_coefficient(Re), SIM.sphere_drag_coefficient(Re, "stokes")) < 1e-15
    nu, a = 1.5e-5, 5e-6
    U = Re * nu / (2 * a)
    CD = 6 * np.pi * 1.2 * nu * a * U / (0.5 * 1.2 * U ** 2 * np.pi * a ** 2)
    assert maxrel(ch08.stokes_drag_coefficient(Re), CD) < 1e-13
    # R20: D = f(μ, U, a) has a single Π group D/(μUa)
    grp = ch01.pi_groups({"D": "N", "mu": "Pa*s", "U": "m/s", "a": "m"}, solution="D")
    assert len(grp) == 1


def test_terminal_velocity_V1_force_balance_round_trip_and_warning():  # V1 (N90): drag = effective weight; radius back
    rho_p, rho, mu = 1000.0, 1.2, 1.8e-5
    a = np.array([2e-6, 5e-6, 10e-6])
    U = ch08.terminal_velocity(a, rho_p, rho, mu)
    W = 4 / 3 * np.pi * a ** 3 * ch08.G0 * (rho_p - rho)
    assert maxrel(6 * np.pi * mu * a * U, W) < 1e-14
    assert maxrel(ch08.radius_from_terminal_velocity(U, rho_p, rho, mu), a) < 1e-14
    st = ch08.settling_state(10e-6, rho_p, rho, mu)
    assert st["U_t"] == pytest.approx(0.0121, rel=1e-2) and st["Re"] == pytest.approx(0.016, rel=2e-2)
    assert st["D"] == pytest.approx(4.1e-11, rel=1e-2) and st["D"] == pytest.approx(st["weight_eff"], rel=1e-14) and st["valid"]
    assert st["C_D"] == pytest.approx(24 / st["Re"])
    with pytest.warns(RuntimeWarning):
        ch08.terminal_velocity(200e-6, rho_p, rho, mu)  # Re > 0.1


@needs_ref
def test_millikan_V5_synthetic_experiment_recovers_e():  # V5 (N91): CODATA e recovered within 1 % (seed 0, 40 drops)
    e = ref_json()["elementary_charge"]["e_C"]
    assert ch08.E_CHARGE == e
    r = ch08.synthetic_millikan()
    assert abs(r["e_est"] / e - 1) < 0.01 and np.all(r["n"] == r["n_true"])
    # millikan_charge is the balance neE = 6πμa(U_rise + U_fall) exactly (noise-free)
    a, rho_p, rho, mu, E = 0.6e-6, 886.0, 1.2, 1.83e-5, 3.2e5
    Uf = float(ch08.terminal_velocity(a, rho_p, rho, mu))
    Ur = (3 * e * E - 4 / 3 * np.pi * a ** 3 * ch08.G0 * (rho_p - rho)) / (6 * np.pi * mu * a)
    assert float(ch08.millikan_charge(Uf, Ur, rho_p, rho, mu, E)) == pytest.approx(3 * e, rel=1e-12)
    assert float(ch08.millikan_charge(Uf, Ur, rho_p, rho, mu, -E)) < 0  # wrong variant: the field sign
    q = np.array([2, 5, 3, 1, 4]) * e * (1 + 1e-4 * RNG.standard_normal(5))
    assert ch08.estimate_elementary_charge(q)["e"] == pytest.approx(e, rel=1e-3)


# =====================================================================================================================
# C15 — far-field breakdown, Oseen (8.53), drag laws
# =====================================================================================================================
def test_inertia_viscous_ratio_V7_linear_growth_in_r_and_Re():  # V7 (C15, D32): log–log slope 1 ± 0.05; ∝ Ua/ν
    r = np.array([50.0, 100.0, 200.0, 500.0])
    for th in (np.pi / 3, np.pi / 2, 2.5):
        rat = ch08.inertia_viscous_ratio(r, th, 1.0, 1.0, 1.0)
        assert abs(np.polyfit(np.log(r), np.log(rat), 1)[0] - 1.0) < 0.05
    r1 = ch08.inertia_viscous_ratio(100.0, 1.0, 1.0, 1.0, 1.0)
    assert ch08.inertia_viscous_ratio(100.0, 1.0, 1.0, 1.0, 0.1) == pytest.approx(10 * r1, rel=1e-12)  # ∝ 1/ν (exact)
    # ∝ U up to the round-off of the nested finite differences: eps·r³/(h_rel²·r²·a) ≈ 2e-6 at r = 100a
    assert ch08.inertia_viscous_ratio(100.0, 1.0, 3.0, 1.0, 1.0) == pytest.approx(3 * r1, rel=2e-5)
    # near the sphere inertia is weak for small Re_a, far away it wins (crossover r ~ a/Re_a)
    Re_a = 0.01
    assert ch08.inertia_viscous_ratio(2.0, 1.0, 1.0, 1.0, 1.0 / Re_a) < 0.1
    # the O(1) prefactor is ≈ 0.05–0.13 (sympy D32: radial ratio = (Re_a r/a)(8cos²θ − 4sin²θ)/(16 cosθ)), so the
    # crossover sits at r ≈ 10–20 a/Re_a: still ∝ 1/Re_a, which is the claim
    ratio_far = ch08.inertia_viscous_ratio(1000.0, 1.0, 1.0, 1.0, 1.0 / Re_a) / (Re_a * 1000.0)
    assert 0.01 < ratio_far < 1.0
    assert ch08.inertia_viscous_ratio(5e4, 1.0, 1.0, 1.0, 1.0 / Re_a) > 1.0


def test_far_field_V2_derivation_D32():  # V2 derivation D32 ★★: u − U ~ Ua/r, ρu·∇u ~ ρU²a/r², μ∇²u ~ μUa/r³
    r, th, ph = CL.coordinates("spherical")
    a, U, nu = sp.symbols("a U nu", positive=True)
    ur = U * sp.cos(th) * (1 - 3 * a / (2 * r) + a ** 3 / (2 * r ** 3))
    ut = -U * sp.sin(th) * (1 - 3 * a / (4 * r) - a ** 3 / (4 * r ** 3))
    dist = (ur - U * sp.cos(th), ut + U * sp.sin(th))
    for c in dist:  # step 1: the disturbance decays as 1/r
        lim = sp.limit(c * r / (U * a), r, sp.oo)
        assert lim != 0 and lim.is_finite is not False
    adv = CL.advective_acceleration([ur, ut, 0], "spherical")
    lap = CL.vector_laplacian([ur, ut, 0], "spherical")
    th0 = sp.pi / 3
    # (the radial component carries the O(1/r²) inertia; the θ-component's leading term is O(a²/r³))
    A2 = sp.limit((adv[0] * r ** 2 / (U ** 2 * a)).subs(th, th0), r, sp.oo)  # steps 2–3: inertia ~ U²a/r²
    V3 = sp.limit((nu * lap[0] * r ** 3 / (nu * U * a)).subs(th, th0), r, sp.oo)  # step 4: viscous ~ νUa/r³
    assert A2 != 0 and A2.is_finite and V3 != 0 and V3.is_finite
    ratio = sp.limit((adv[0] / (nu * lap[0]) / r).subs(th, th0), r, sp.oo)  # steps 5–7: ratio ~ (U/ν) r
    assert z0(ratio - U / nu * A2 / V3)
    assert sp.simplify(ratio * nu / U).is_number
    assert sp.limit((adv[1] * r ** 2 / (U ** 2 * a)).subs(th, th0), r, sp.oo) == 0


def test_oseen_V2_limit_linearisation_and_D33():  # V2 derivation D33 ★★ + N95: Re → 0 of (8.53) is (8.48); the sign of ∇p
    d = ch08.oseen_limit_sympy()
    assert d["difference"] == 0
    r, th, a, U = sp.symbols("r theta a U", positive=True)
    Re = sp.Symbol("Re", positive=True)
    s = Re / 4 * r / a * (1 - sp.cos(th))  # step 1
    assert z0(sp.series(1 - sp.exp(-sp.Symbol("s")), sp.Symbol("s"), 0, 2).removeO() - sp.Symbol("s"))  # step 3
    term = 3 / Re * (1 + sp.cos(th)) * s  # step 4
    assert z0(term - sp.Rational(3, 4) * r / a * sp.sin(th) ** 2)  # step 5
    stokes = U * a ** 2 * (r ** 2 / (2 * a ** 2) + a / (4 * r) - 3 * r / (4 * a)) * sp.sin(th) ** 2  # steps 6–7
    assert z0(stokes - U * r ** 2 * sp.sin(th) ** 2 * (sp.Rational(1, 2) - 3 * a / (4 * r) + a ** 3 / (4 * r ** 3)))
    lin = ch08.oseen_linearisation_sympy()
    assert lin["expansion_residual"] == 0
    # wrong variant (R13): the printed +∂p/∂x differs from the (8.43)-consistent −∂p/∂x
    assert lin["oseen_equation"] != lin["oseen_equation_printed"]
    xs, ys, zs = sp.symbols("x y z", real=True)
    pfun = sp.Function("p")(xs, ys, zs)
    assert z0(lin["oseen_equation_printed"].rhs - lin["oseen_equation"].rhs - 2 * sp.diff(pfun, xs))
    assert z0(lin["oseen_equation"].rhs.coeff(sp.diff(pfun, xs)) + 1)  # −∂p/∂x, consistent with (8.43)


def test_oseen_V1_velocity_axis_wake_and_no_slip_order():  # V1 (N97, N100): u = (6.83) of ψ; ψ = 0 on the axis; O(Re) slip
    U, a, Re = 1.0, 1.0, 1.0
    r = np.linspace(1.2, 6, 20)
    th = np.linspace(0.1, np.pi - 0.1, 25)
    R, T = np.meshgrid(r, th)
    d = 1e-6
    psi = lambda rr, tt: np.asarray(ch08.oseen_streamfunction(rr, tt, U, a, Re))  # noqa: E731
    urf = (psi(R, T + d) - psi(R, T - d)) / (2 * d) / (R ** 2 * np.sin(T))
    utf = -(psi(R + d, T) - psi(R - d, T)) / (2 * d) / (R * np.sin(T))
    ur, ut = ch08.oseen_velocity(R, T, U, a, Re)
    assert maxrel(ur, urf) < 1e-6 and maxrel(ut, utf) < 1e-6
    assert np.max(np.abs(ch08.oseen_streamfunction(r, 0.0, U, a, Re))) < 1e-14
    assert np.max(np.abs(ch08.oseen_streamfunction(r, np.pi, U, a, Re))) < 1e-12
    # Re → 0: Stokes exactly (algebraically stable form)
    assert maxrel(ch08.oseen_streamfunction(R, T, U, a, 1e-10), ch08.stokes_sphere_streamfunction(R, T, U, a)) < 1e-9
    assert maxrel(ch08.oseen_velocity(R, T, U, a, 0.0)[1], ch08.stokes_sphere_velocity(R, T, U, a)[1]) < 1e-14
    # fore–aft asymmetry (a wake) in the fluid frame at Re = 1
    down = np.abs(ch08.oseen_streamfunction(3.0, 0.4, U, a, 1.0, "fluid"))
    up = np.abs(ch08.oseen_streamfunction(3.0, np.pi - 0.4, U, a, 1.0, "fluid"))
    assert abs(down - up) > 1e-2 * max(down, up)
    # no slip only to O(Re): the wall speed scales linearly with Re (R18)
    w = [np.max(np.hypot(*ch08.oseen_velocity(a, np.linspace(0.1, 3.0, 30), U, a, R_))) for R_ in (1e-3, 1e-2, 1e-1)]
    assert abs(observed_order([1e-3, 1e-2, 1e-1], w) - 1.0) < 0.05


@needs_ref
def test_drag_laws_V5_oseen_proudman_pearson_morrison():  # V5 (N99, N101): radius vs diameter Re; Morrison between
    ref = ref_json()["oseen_drag"]
    Re = np.logspace(-2, np.log10(5), 30)
    Ra = Re / 2
    k = ref["oseen_coeff_radius"]
    assert maxrel(ch08.oseen_drag_coefficient(Re), 12 / Ra * (1 + k * Ra)) < 1e-14  # Wikipedia's radius form
    assert maxrel(ch08.proudman_pearson_drag_coefficient(Re), 12 / Ra * (1 + k * Ra + ref["pp_coeff_radius"] * Ra ** 2 * np.log(Ra))) < 1e-14
    # wrong variant: Oseen's 3/8 applied to the diameter Re
    assert np.max(np.abs(np.asarray(ch08.oseen_drag_coefficient(Re)) / (24 / Re * (1 + 3 / 8 * Re)) - 1)) > 0.05
    # all → 24/Re as Re → 0
    small = np.array([1e-6, 1e-5])
    for fn in (ch08.oseen_drag_coefficient, ch08.proudman_pearson_drag_coefficient):
        assert maxrel(np.asarray(fn(small)) * small / 24, np.ones(2)) < 1e-5
    # the book's claim: experiments (Morrison correlation, ch04 V5) lie between Stokes and Oseen for 0.1 ≤ Re ≤ 5
    Rb = np.logspace(-1, np.log10(5), 40)
    m = np.asarray(SIM.sphere_drag_coefficient(Rb, "morrison"))
    assert np.all(m >= np.asarray(ch08.stokes_drag_coefficient(Rb))) and np.all(m <= np.asarray(ch08.oseen_drag_coefficient(Rb)))


# =====================================================================================================================
# Reused functions called with ch08 parameters (design Part C C.0) and the contract (C.1–C.4)
# =====================================================================================================================
def test_reuse_V1_couette_startup_ftcs_and_navier_stokes_presets():  # V1 smoke (R14, C.0): ch01/ch04 reuse sanity
    h, U, nu = 0.01, 0.1, 1e-6
    y = np.linspace(0, h, 21)
    late = ch08.couette_startup_profile(y, 1e6, U, h, nu)
    assert maxrel(late, ch08.channel_flow(y, h, U, 0.0)) < 1e-10  # Exercise 8.31 → plane Couette
    early = ch08.couette_startup_profile(y, 1.0, U, h, nu)  # short times: Stokes' first problem from the top plate
    assert maxrel(early, ch08.stokes_first_problem(h - y, 1.0, U, nu), floor=U) < 1e-6
    # the ch04 presets that ch08 now derives satisfy (4.39b) at ch08-sized points (residual ≪ the viscous term)
    for name, scale_key in (("couette", "h"), ("poiseuille", "h"), ("stokes_first", None), ("lamb_oseen", None), ("solid_body", None)):
        uf, pf, q = NS.exact_solution_fields(name)
        X = np.stack([np.array([0.001, 0.002, 0.003]), np.array([0.002, 0.004, 0.006])])
        t = 1.0
        terms = NS.ns_incompressible_terms(uf, pf, X, t, q["rho"], q["mu"], (0.0, 0.0), h=1e-4, ht=1e-4)
        ref = np.max(np.abs(terms.pressure)) + np.max(np.abs(terms.viscous)) + np.max(np.abs(terms.advective))
        # stencil h = 1e-4 m against √(νt) = 1 mm for the unsteady presets: truncation ≲ (h/√(νt))²/12 ≈ 1e-3
        assert np.max(np.abs(terms.residual)) < 2e-3 * ref + 1e-12, name  # (Couette with G = 0 has no forces at all)


def test_part_c_V1_every_contract_function_exists_and_is_scalar_callable():  # V1 smoke: design Part C C.1–C.4
    names = ["channel_flow", "channel_flow_rate", "channel_shear_stress", "channel_backflow_threshold",
             "couette_poiseuille_state", "pipe_poiseuille", "pipe_shear_stress", "pipe_wall_stress", "pipe_flow_rate",
             "pipe_friction_factor", "circular_couette", "circular_couette_pressure", "circular_couette_shear_stress",
             "circular_couette_power", "circular_couette_state", "similarity_variable", "stokes_first_problem",
             "stokes_first_vorticity", "stokes_first_stopped", "stokes_first_state", "diffusion_thickness",
             "transition_width", "vortex_sheet_diffusion", "temporal_bl_wall_stress", "line_vortex_decay",
             "line_vortex_spinup", "stokes_second_problem", "stokes_layer", "stokes_layer_state", "crank_nicolson_1d",
             "lubrication_scales", "lubrication_term_magnitudes", "lubrication_velocity", "lubrication_flux",
             "reynolds_pressure_1d", "slider_bearing", "slider_bearing_load", "slider_optimum_taper",
             "slider_bearing_state", "slider_gap_velocity", "hele_shaw_velocity", "hele_shaw_potential", "thin_film_flux",
             "thin_film_spread", "viscous_current_similarity", "thin_film_state", "stokes_residual", "E2", "E4_residual",
             "stokes_sphere_streamfunction", "stokes_sphere_velocity", "stokes_sphere_velocity_xyz",
             "stokes_sphere_pressure", "stokes_sphere_surface_stresses", "stokes_drag", "stokes_drag_running",
             "sphere_drag_quadrature", "stokes_drag_coefficient", "oseen_drag_coefficient",
             "proudman_pearson_drag_coefficient", "terminal_velocity", "radius_from_terminal_velocity", "millikan_charge",
             "settling_state", "inertia_viscous_ratio", "oseen_streamfunction", "oseen_velocity", "side_line_speed",
             "stokes_sphere_sympy", "pipe_flow_regime", "inertia_viscous_scales", "momentum_diffusivity",
             "diffusion_time", "wall_bc_residuals", "parallel_flow_sympy", "advective_acceleration_check",
             "lubrication_nondim_sympy", "reynolds_equation_sympy", "slider_bearing_sympy", "hele_shaw_cylinder",
             "similarity_reduce_sympy", "similarity_ode_solve", "vorticity_content", "similarity_collapse_error",
             "stokes_second_sympy", "low_re_scaling_sympy", "oseen_linearisation_sympy", "oseen_limit_sympy",
             "synthetic_millikan", "couette_startup_profile", "ftcs_diffusion_1d", "hele_shaw_mean_velocity",
             "hele_shaw_streamfunction_grid", "thin_film_velocity", "viscous_current_eta_N", "stokes_layer_envelope",
             "estimate_elementary_charge", "stokes_first_pi_groups", "reynolds_number"]
    missing = [n for n in names if not callable(getattr(ch08, n, None))]
    assert not missing, missing
    assert ch08.G_BOOK == 9.81 and ch08.G0 == 9.80665 and ch08.P_ATM == 101325.0
    for fn, args in ((ch08.couette_poiseuille_state, (0.01, 0.1, 4.0)), (ch08.circular_couette_state, (0.01, 0.02, 1.0, 0.0)),
                     (ch08.stokes_first_state, (100.0,)), (ch08.stokes_layer_state, (1e-6, 6.0, 1e-3)),
                     (ch08.slider_bearing_state, (5e-5, 0.3, 0.05, 5.0, 0.05)), (ch08.thin_film_state, (10.0, 1e-4)),
                     (ch08.settling_state, (1e-5, 1000.0, 1.2, 1.8e-5)), (ch08.lubrication_scales, (0.05, 5e-5, 5.0, 870.0, 0.05)),
                     (ch08.stokes_layer, (1e-6, 6.0)), (ch08.temporal_bl_wall_stress, (1.0, 0.1))):
        out = fn(*args)
        assert all(np.ndim(v) == 0 for v in out.values()), fn.__name__  # scalar-callable (explainer parity rows)
    for fn, args in ((ch08.channel_flow, (0.003, 0.01, 0.1, 4.0)), (ch08.pipe_poiseuille, (0.001, 0.002, -10.0)),
                     (ch08.circular_couette, (0.015, 0.01, 0.02, 1.0, 0.0)), (ch08.stokes_first_problem, (0.01, 100.0)),
                     (ch08.stokes_second_problem, (1e-3, 0.3)), (ch08.slider_bearing, (0.02, 5e-5, 0.3, 0.05, 5.0)),
                     (ch08.line_vortex_decay, (0.01, 100.0, 0.05)), (ch08.stokes_drag_coefficient, (0.1,)),
                     (ch08.oseen_drag_coefficient, (0.5,)), (ch08.terminal_velocity, (5e-6, 1000.0, 1.2, 1.8e-5)),
                     (ch08.diffusion_thickness, (100.0,)), (ch08.stokes_sphere_pressure, (2.0, 1.0))):
        assert np.ndim(fn(*args)) == 0 and np.isfinite(float(fn(*args))), fn.__name__


@pytest.mark.slow
def test_scripts_V1_every_ch08_script_runs():  # V1 smoke (I41): all 12 scripts exit 0 headless
    env = dict(os.environ, MPLBACKEND="Agg")
    scripts = [s for s in sorted((ROOT / "scripts").glob("ch08_*.py")) if s.name not in ("ch08_drawings.py", "ch08_common.py")]
    assert len(scripts) == 12
    for scr in scripts:
        r = subprocess.run([sys.executable, str(scr), "--no-show"], cwd=ROOT, env=env, capture_output=True, text=True, timeout=600)
        assert r.returncode == 0, (scr.name, r.stderr[-2000:])


# =====================================================================================================================
# V6 — the book's printed forms and numbers (private JSON; skipped when absent)
# =====================================================================================================================
def _bsyms():
    names = ["y", "h", "U", "mu", "nu", "rho", "dpdx", "dpdz", "R", "a", "R1", "R2", "Omega1", "Omega2", "x", "L", "h0",
             "alpha", "p_e", "t", "eta", "omega", "r", "theta", "Gamma", "Re", "g", "U_h", "U_0", "z", "h_x", "C1"]
    return {n: sp.Symbol(n, real=True) for n in names}


def _ev(expr: str, vals: dict) -> float:
    S = _bsyms()
    return float(sp.sympify(expr, locals=S).subs({S[k]: v for k, v in vals.items()}).evalf())


@book_only
def test_book_V6_section_8_1_and_8_2_forms_and_numbers():  # V6 §8.1–§8.2: Re band, ν ratio, (8.5)–(8.12) as printed
    b = book()
    v = b["sec8_1_values"]
    assert ch08.pipe_flow_regime(v["Re_transition_low"] * 1e-6 / 0.01, 0.01, 1e-6)[1] == "transitional"
    ratio = ch08.momentum_diffusivity("air", 293.15) / ch08.momentum_diffusivity("water", 293.15)
    assert ratio == pytest.approx(v["air_vs_water_diffusivity_ratio"], rel=5e-3)
    assert ch08.momentum_diffusivity("air", 293.15) == pytest.approx(v["nu_air_m2_s"], rel=0.01)
    assert ch08.momentum_diffusivity("water", 293.15) == pytest.approx(v["nu_water_m2_s"], rel=0.01)
    f = b["sec8_2_closed_forms"]
    vals = dict(y=0.0031, h=0.01, U=0.12, mu=1.3e-3, dpdx=-4.5, R=0.0021, a=0.004, dpdz=-30.0, R1=0.01, R2=0.025,
                Omega1=1.7, Omega2=-0.6)
    assert _ev(f["channel_u_8_5"], vals) == pytest.approx(float(ch08.channel_flow(0.0031, 0.01, 0.12, -4.5, 1.3e-3)), rel=1e-12)
    Q, V = ch08.channel_flow_rate(0.01, 0.12, -4.5, 1.3e-3)
    assert _ev(f["channel_Q"], vals) == pytest.approx(Q, rel=1e-12) and _ev(f["channel_V"], vals) == pytest.approx(V, rel=1e-12)
    assert _ev(f["poiseuille_tau"], vals) == pytest.approx(float(ch08.channel_shear_stress(0.0031, 0.01, 0.0, -4.5, 1.3e-3)), rel=1e-12)
    assert _ev(f["pipe_u_8_6"], vals) == pytest.approx(float(ch08.pipe_poiseuille(0.0021, 0.004, -30.0, 1.3e-3)), rel=1e-12)
    assert _ev(f["pipe_tau0_8_8"], vals) == pytest.approx(float(ch08.pipe_wall_stress(0.004, -30.0)), rel=1e-12)
    Qp, Vp, _ = ch08.pipe_flow_rate(0.004, -30.0, 1.3e-3)
    assert _ev(f["pipe_Q"], vals) == pytest.approx(Qp, rel=1e-12) and _ev(f["pipe_V"], vals) == pytest.approx(Vp, rel=1e-12)
    vals_c = dict(vals, R=0.017)
    u, A, B = ch08.circular_couette(0.017, 0.01, 0.025, 1.7, -0.6, return_coeffs=True)
    assert _ev(f["circ_couette_u_8_10"], vals_c) == pytest.approx(u, rel=1e-12)
    assert _ev(f["circ_couette_A"], vals_c) == pytest.approx(A, rel=1e-12) and _ev(f["circ_couette_B"], vals_c) == pytest.approx(B, rel=1e-12)
    assert _ev(f["rotating_cylinder_8_11"], vals_c) == pytest.approx(float(ch08.circular_couette(0.017, 0.01, np.inf, 1.7, 0.0)), rel=1e-12)
    assert _ev(f["sigma_Rphi"], vals_c) == pytest.approx(float(ch08.circular_couette_shear_stress(0.017, 0.01, np.inf, 1.7, 0.0, 1.3e-3)), rel=1e-12)
    assert _ev(f["solid_body_8_12"], vals_c) == pytest.approx(float(ch08.circular_couette(0.017, 0.0, 0.025, 0.0, -0.6)), rel=1e-12)


@book_only
def test_book_V6_section_8_3_lubrication_forms_and_numbers():  # V6 §8.3: ε²Re_L of the oil example; Example 8.1 forms
    b = book()
    s3 = b["sec8_3_lubrication"]
    sc = ch08.lubrication_scales(s3["example_length_m"], s3["example_gap_m"], s3["example_speed_m_s"], 1.0, s3["example_nu_oil_m2_s"])
    assert sc["eps2_Re_L"] == pytest.approx(s3["example_eps2_ReL"], rel=5e-3)
    vals = dict(x=0.013, L=0.05, h0=5e-5, alpha=0.37, U=4.0, mu=0.03, p_e=1e5, y=2e-5, h=6e-5, dpdx=-2e6, U_h=4.0, U_0=0.0)
    assert _ev(s3["ex8_1_p_minus_pe_exact"], vals) + 1e5 == pytest.approx(float(ch08.slider_bearing(0.013, 5e-5, 0.37, 0.05, 4.0, 0.03, 1e5)), rel=1e-12)
    assert _ev(s3["ex8_1_p_minus_pe_as_printed"], vals) + 1e5 == pytest.approx(float(ch08.slider_bearing(0.013, 5e-5, 0.37, 0.05, 4.0, 0.03, 1e5, "book")), rel=1e-12)
    assert _ev(s3["ex8_1_p_minus_pe_as_printed"], vals) != pytest.approx(_ev(s3["ex8_1_p_minus_pe_exact"], vals), rel=1e-3)
    assert _ev(s3["ex8_1_p_linear_alpha"], vals) + 1e5 == pytest.approx(float(ch08.slider_bearing(0.013, 5e-5, 0.37, 0.05, 4.0, 0.03, 1e5, "linear")), rel=1e-12)
    assert _ev(s3["ex8_1_W"], vals) == pytest.approx(ch08.slider_bearing_load(5e-5, 0.37, 0.05, 4.0, 0.03, "linear"), rel=1e-12)
    assert _ev(s3["ex8_1_C1"], vals) == pytest.approx(ch08.slider_bearing_state(5e-5, 0.37, 0.05, 4.0, 0.03)["C1"], rel=1e-12)
    assert _ev(s3["u_8_19"], vals) == pytest.approx(float(ch08.lubrication_velocity(2e-5, 6e-5, -2e6, 4.0, 0.0, 0.03, "book")), rel=1e-12)
    assert _ev(s3["ex8_3_flux"], dict(h=3e-3, h_x=-0.02, rho=1260.0, g=9.81, mu=1.4)) == pytest.approx(float(ch08.thin_film_flux(3e-3, -0.02, 1260.0, 9.81, 1.4)), rel=1e-12)


@book_only
def test_book_V6_section_8_4_to_8_6_numbers_and_slips():  # V6 §8.4–§8.6: 3.64, 2.76 (slip), 5.54, 0.06, 24/Re, drag parts
    b = book()
    s4, s5, s6 = b["sec8_4_similarity"], b["sec8_5_oscillating_plate"], b["sec8_6_sphere"]
    eta99 = ch08.stokes_first_state(1.0)["eta_edge"]
    assert eta99 == pytest.approx(s4["eta_99_8_31"], rel=5e-3) and round(eta99, 2) == s4["eta_99_8_31"]
    w = ch08.transition_width(1.0, 1.0)
    assert w == pytest.approx(s4["ex8_5_width_coeff"], rel=5e-3) and round(w, 2) == s4["ex8_5_width_coeff"]
    # R11: the printed ±2.76 is not the 95 % point (0.43 % off; the book's own width 5.54 = 2 × 2.772)
    assert abs(w / 2 / s4["ex8_5_eta_95_as_printed"] - 1) > 3e-3
    assert ch08.vorticity_content(10.0, 1.0) == pytest.approx(1.0) and s4["vorticity_integral_as_printed"] == "-U"
    assert z0(sp.sympify(s4["A_const"]) - ch08.similarity_reduce_sympy("stokes1")["A"])
    assert (float(ch08.similarity_reduce_sympy("vortex_sheet")["n"]), float(ch08.similarity_reduce_sympy("spreading")["n"])) == (s4["ex8_5_n"], s4["ex8_7_n"])
    amp = ch08.stokes_layer(1e-6, 1.0)["amplitude_at_delta_book"]
    assert round(amp, 2) == s5["amplitude_ratio_at_4"]  # reproduces the printed ≈ 0.06 to its one significant figure
    assert ch08.stokes_layer(1e-6, 1.0)["delta_book"] == pytest.approx(s5["depth_coeff_4"] * np.sqrt(1e-6 / 1.0))
    mu, a, U = 1e-3, 0.002, 0.01
    pmax = float(ch08.stokes_sphere_pressure(a, np.pi, U, a, mu)) / (mu * U / a)
    assert pmax == pytest.approx(s6["p_max_over_muU_a"], rel=1e-12)
    assert float(ch08.stokes_sphere_pressure(a, 0.0, U, a, mu)) == pytest.approx(-_ev("3*mu*U/(2*a)", dict(mu=mu, U=U, a=a)))
    assert float(ch08.stokes_sphere_pressure(a, 0.0, U, a, mu)) != pytest.approx(_ev(s6["p_min_as_printed"], dict(mu=mu, U=U, a=a)))
    parts = ch08.stokes_drag(mu, a, U, parts=True)
    assert parts["pressure"] / parts["total"] == pytest.approx(float(sp.Rational(s6["pressure_fraction"])))
    assert parts["friction"] / parts["total"] == pytest.approx(float(sp.Rational(s6["friction_fraction"])))
    assert _ev(s6["drag_8_51"], dict(mu=mu, a=a, U=U)) == pytest.approx(parts["total"], rel=1e-12)
    assert _ev(s6["CD_8_52"], dict(Re=0.3)) == pytest.approx(float(ch08.stokes_drag_coefficient(0.3)))
    assert _ev(s6["CD_oseen"], dict(Re=0.3)) == pytest.approx(float(ch08.oseen_drag_coefficient(0.3)))
    vals = dict(r=2.3, theta=1.1, U=0.7, a=1.0, Re=0.4, mu=mu)
    assert _ev(s6["psi_8_48"], vals) == pytest.approx(float(ch08.stokes_sphere_streamfunction(2.3, 1.1, 0.7, 1.0)), rel=1e-12)
    assert _ev(s6["psi_fluid_frame"], vals) == pytest.approx(float(ch08.stokes_sphere_streamfunction(2.3, 1.1, 0.7, 1.0, "fluid")), rel=1e-12)
    assert _ev(s6["u_r_8_49"], vals) == pytest.approx(float(ch08.stokes_sphere_velocity(2.3, 1.1, 0.7, 1.0)[0]), rel=1e-12)
    assert _ev(s6["u_theta_8_49"], vals) == pytest.approx(float(ch08.stokes_sphere_velocity(2.3, 1.1, 0.7, 1.0)[1]), rel=1e-12)
    assert _ev(s6["psi_oseen_8_53"], vals) == pytest.approx(float(ch08.oseen_streamfunction(2.3, 1.1, 0.7, 1.0, 0.4)), rel=1e-12)
    for key, fn in (("ex8_5_u", lambda: ch08.vortex_sheet_diffusion(0.002, 3.0, 0.4, 1e-6)[0]),
                    ("ex8_5_omega_z", lambda: ch08.vortex_sheet_diffusion(0.002, 3.0, 0.4, 1e-6)[1]),
                    ("ex8_6_u_theta", lambda: ch08.line_vortex_decay(0.002, 3.0, 0.05, 1e-6)),
                    ("ex8_26_u_theta", lambda: ch08.line_vortex_spinup(0.002, 3.0, 0.05, 1e-6)),
                    ("u_8_30", lambda: ch08.stokes_first_problem(0.002, 3.0, 0.4, 1e-6))):
        assert _ev(s4[key], dict(y=0.002, r=0.002, t=3.0, U=0.4, nu=1e-6, Gamma=0.05)) == pytest.approx(float(fn()), rel=1e-10)
    assert _ev(s5["u_8_38"], dict(y=4e-4, t=0.3, U=0.4, omega=6.0, nu=1e-6)) == pytest.approx(float(ch08.stokes_second_problem(4e-4, 0.3, 0.4, 6.0, 1e-6)), rel=1e-12)
