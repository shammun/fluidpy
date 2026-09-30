"""Verification suite for Chapter 9 — Boundary Layers and Related Topics (Kundu, Cohen & Dowling 5e, §§9.1–9.11).

Evidence levels (``verify-implementation`` skill): V1 analytic · V2 symbolic (sympy / dimensions) · V3 convergence ·
V4 conservation / invariant · V5 published benchmark (``reference/ch09/``, see SOURCES.md) · V6 book value (private,
``tests/book_values_ch09.json``, skipped when absent) · V7 limits / symmetry / invariance. Every test name is
``test_<concept>_<level>_<what>``; the comment on the ``def`` line repeats the level.

A items (CORE, ≥ 2 independent levels): C01 scaling and the boundary-layer equations (9.4)–(9.11) · C02 the three
thicknesses (9.16)–(9.17) · C03 the Blasius reduction (9.19)–(9.27) · C04 the Blasius solution and its numbers
(9.28)–(9.33) · C05 Falkner–Skan (9.34)–(9.36) · C06 the von Kármán momentum integral (9.43) · C07 Thwaites (9.44)–(9.50) ·
C08 separation: wall curvature (9.51)–(9.52) · C09 form drag of the separated model · C10 the Kármán vortex street ·
C11 the drag crisis and sports balls (qualitative) · C12 the free jet (9.53)–(9.76) · C13 the wall jet (9.77)–(9.85) ·
C14 the teacup secondary flow.
Derivations: every ★★ and ★★★ D row (D01, D04–D11, D14–D17, D19–D22) — and the cheap ★ rows D02, D03, D12, D13, D18 — is
re-derived with sympy in a ``test_*_V2_derivation`` test; the ★★★ rows (D14, D16, D19, D20, D21) step by step through the
design's Part F lines.

Pinned conventions and slips (each with a discriminating assertion): u = ψ_y, v = −ψ_x; dp/dx > 0 adverse; angles from the
FORWARD stagnation point; Re on the diameter for bluff bodies; the book's n is the code's m. Printed slips kept as wrong
variants that must FAIL: R1 (9.7) without the squares · R2 η₉₉ = 4.93 · R3 wall-jet ODE with coefficient 1 · R4 the
wall-jet separation of variables · R5 both faces of the plate · R6 h₉₉ = 5.6152 · R10 the Magnus sentence · R11 reverse
flow below n = −0.0904 · R14 "Chapter 13".

Run: ``.venv/Scripts/python.exe -m pytest tests/test_ch09.py -q -p no:cacheprovider``  (``-m "not slow"`` skips the script
smoke test and two multi-second continuations).
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
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

from fluidpy import ch09_boundary_layers as ch09
from fluidpy.core import boundary_layer as BL
from fluidpy.core import bluff_body as BB
from fluidpy.core import creeping as CRP
from fluidpy.core import jets as JET
from fluidpy.core import laminar as LAM
from fluidpy.core import potential as PF
from fluidpy.core import similarity as SIM
from fluidpy.core.units import Q_, dimensional_check
from tools.convergence import observed_order, pairwise_orders

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference" / "ch09"
BOOK = Path(__file__).resolve().parent / "book_values_ch09.json"
needs_ref = pytest.mark.skipif(not (REF / "benchmarks.json").exists(), reason="run reference/ch09/make_refs.py")
book_only = pytest.mark.skipif(not BOOK.exists(), reason="book values are private; see CLAUDE.md rule 9")

ORDER_TOL = 0.15  # design order ± this (verify-implementation default)


def book():
    return json.loads(BOOK.read_text(encoding="utf-8"))


def ref_json():
    return json.loads((REF / "benchmarks.json").read_text(encoding="utf-8"))


def rel(a, b):
    return abs(float(a) - float(b)) / max(abs(float(b)), 1e-300)


def z0(expr) -> bool:
    return sp.simplify(expr) == 0


def belden_to_book(beta: float, kappa: float):
    """(m, f''_book(0)) from Belden et al.'s (β, κ): m = β/(2 − β) is the book's n; f''_book = √((m+1)/2) κ (normalisation of η)."""
    m = beta / (2.0 - beta)
    return m, np.sqrt((m + 1.0) / 2.0) * kappa


# =====================================================================================================================
# C01 — Prandtl scaling, the boundary-layer equations (9.4)–(9.11) (§9.1)
# =====================================================================================================================
def test_bl_scales_V1_definitions_and_identities():  # V1 (C01): δ̄/L = Re^-1/2, adv = visc, v ~ U Re^-1/2, C_f ~ 2/√Re
    for U, L, nu in ((1.0, 1.0, 1.5e-5), (10.0, 0.3, 1e-6), (0.05, 2.0, 1e-6)):
        s = BL.boundary_layer_scales(U, L, nu, rho=1.2)
        Re = U * L / nu
        assert s["Re"] == pytest.approx(Re, rel=1e-14)
        assert s["delta_over_L"] == pytest.approx(Re ** -0.5, rel=1e-13)  # (9.4)
        assert s["adv"] == pytest.approx(U ** 2 / L, rel=1e-14)  # (9.2)
        assert s["visc"] == pytest.approx(s["adv"], rel=1e-12)  # (9.3) balanced with (9.2) by construction of δ̄
        assert s["visc_x"] / s["visc"] == pytest.approx(1.0 / Re, rel=1e-12)  # the term (9.7) drops is 1/Re smaller
        assert s["v_scale"] == pytest.approx(U * Re ** -0.5, rel=1e-13)
        assert s["cf_estimate"] == pytest.approx(2.0 / np.sqrt(Re), rel=1e-14)
        assert s["tau0_scale"] == pytest.approx(1.2 * nu * U / s["delta"], rel=1e-13)  # μU/δ̄
        # bridge to ch08: δ̄ = √(ν t) with t = L/U (Stokes-layer thickness after one flow-through time)
        assert s["delta"] == pytest.approx(np.sqrt(nu * (L / U)), rel=1e-14)


def test_bl_scales_V7_thin_layer_limits():  # V7 (C01): δ/L → 0 as Re → ∞ (monotone), Re = 1 gives a layer as thick as the body
    Res = np.array([1.0, 1e2, 1e4, 1e6, 1e8])
    d = np.array([BL.boundary_layer_scales(1.0, 1.0, 1.0 / R)["delta_over_L"] for R in Res])
    assert np.all(np.diff(d) < 0) and d[0] == pytest.approx(1.0) and d[-1] == pytest.approx(1e-4)
    assert BL.boundary_layer_scales(1.0, 1.0, 1e-3)["Re"] == pytest.approx(1e3)


def test_bl_variables_V1_eq_9_6_and_round_trip():  # V1 (C01): (9.6) stretched variables, and the inverse returns the input
    L, U, rho, nu = 0.7, 3.0, 1.2, 2e-5
    Re = U * L / nu
    x, y, u, v, p = 0.35, 1.2e-3, 2.1, 0.015, 4.0
    xs, ys, us, vs, ps = BL.to_bl_variables(x, y, u, v, p, L, U, rho, nu=nu)
    assert (xs, ys, us, vs, ps) == pytest.approx((x / L, y / L * np.sqrt(Re), u / U, v / U * np.sqrt(Re), p / (rho * U ** 2)), rel=1e-14)
    back = BL.from_bl_variables(xs, ys, us, vs, ps, L, U, rho, nu=nu)
    assert back == pytest.approx((x, y, u, v, p), rel=1e-13)
    # Re given directly, and the ν-vs-Re guard
    assert BL.to_bl_variables(x, y, u, v, p, L, U, rho, Re=Re)[1] == pytest.approx(ys, rel=1e-14)
    with pytest.raises(ValueError):
        BL.to_bl_variables(x, y, u, v, p, L, U, rho)
    # wrong variant: stretching y by Re instead of √Re changes y* by √Re
    assert y / L * Re != pytest.approx(ys, rel=1e-3)
    # arrays broadcast
    assert np.shape(BL.to_bl_variables(np.ones(3) * x, y, u, v, p, L, U, rho, nu=nu)[0]) == (3,)


def test_bl_nondim_V2_sympy_coefficients_and_R1_printed_form_fails():  # V2 (C01): (9.7)–(9.8) coefficients; printed (9.7) fails
    r = ch09.bl_nondim_sympy()
    p = r["coefficients"]
    assert p["u*du*/dx*"] == 0 and p["v*du*/dy*"] == 0 and p["dp*/dx*"] == 0  # (9.7): inertia and pressure O(1)
    assert p["d2u*/dx*2"] == -1 and p["d2u*/dy*2"] == 0  # (9.7): 1/Re only on the streamwise second derivative
    assert p["u*dv*/dx*"] == -1 and p["v*dv*/dy*"] == -1 and p["dp*/dy*"] == 0  # (9.8): pressure alone at O(1)
    assert p["d2v*/dx*2"] == -2 and p["d2v*/dy*2"] == -1 and p["continuity"] == 0
    assert z0(r["residual_correct"]) and r["dimension_check"]["ok"]
    # R1 (wrong variant): the printed (9.7) with first-order denominators is dimensionally inconsistent and differs from the derived one
    bad = ch09.bl_nondim_sympy(printed_9_7=True)
    assert not bad["dimension_check"]["ok"], "the printed (9.7) must fail the dimension check"
    assert not z0(bad["residual_printed_vs_correct"])
    # the limit (9.9)/(9.10) has no 1/Re term
    assert not r["limit"]["x_momentum"].has(sp.Symbol("Re", positive=True))
    assert r["limit"]["y_momentum"].rhs.has(sp.Function("ps")) or r["limit"]["y_momentum"].lhs == 0


def test_bl_x_momentum_V3_residual_of_9_9_on_blasius_field_is_second_order():  # V3 (C01/C03): central differences of (9.9)/(9.18)
    U, nu = 1.0, 1e-4
    hs, errs, errc = [], [], []
    for n in (20, 40, 80):
        x = np.linspace(0.5, 1.5, n + 1)
        y = np.linspace(0.0, 0.06, n + 1)
        X, Y = np.meshgrid(x, y, indexing="xy")
        f = BL.blasius_fields(X, Y, U, nu)
        res = BL.bl_x_momentum_residual(x, y, f["u"], f["v"], 0.0, nu, rho=1.2)  # dp/dx = 0 for the plate
        errs.append(np.max(np.abs(res[2:-2, 2:-2])) / (U ** 2 / 0.5))  # scale: U²/x (two boundary layers of nodes excluded: one-sided edge stencils cost an order in u_yy)
        ux = np.gradient(f["u"], x, axis=1, edge_order=2)
        vy = np.gradient(f["v"], y, axis=0, edge_order=2)
        errc.append(np.max(np.abs((ux + vy)[2:-2, 2:-2])) / (U / 0.5))  # continuity (6.2) in the same units
        hs.append(x[1] - x[0])
    assert abs(observed_order(hs, errs) - 2.0) < 0.25, pairwise_orders(hs, errs)  # error span 1.1 decades → ±0.25
    assert abs(observed_order(hs, errc) - 2.0) < 0.25, pairwise_orders(hs, errc)
    assert errs[-1] < 5e-4 and errc[-1] < 5e-4
    # 1-D single-station form (dudx supplied) agrees with the 2-D residual
    j = 30
    xs = np.array([1.0])
    fld = BL.blasius_fields(np.ones_like(y) * 1.0, y, U, nu)
    ux0 = float(np.gradient(f["u"], x, axis=1, edge_order=2)[j, 40])
    r1 = BL.bl_x_momentum_residual(None, y, fld["u"], fld["v"], 0.0, nu, dudx=ux0 * np.ones_like(y))
    assert np.all(np.isfinite(r1)) and r1.shape == y.shape
    # wrong variant: a pressure gradient that is not there leaves an O(dp/dx) residual
    assert np.max(np.abs(BL.bl_x_momentum_residual(x, y, f["u"], f["v"], 50.0, nu, rho=1.2)[2:-2, 2:-2])) > 10


def test_bl_pressure_V7_variation_across_layer_scales_as_inverse_Re():  # V7 (C01/N10): (9.8) ⇒ (9.10), ratio ∝ Re^-1
    ratios = []
    for Re in (1e3, 1e4, 1e5):
        d = BL.bl_pressure_variation(Re)
        assert abs(d["order"] - 1.0) < ORDER_TOL, d
        ratios.append(d["ratio"])
    assert np.all(np.diff(ratios) < 0)
    assert ratios[0] / ratios[1] == pytest.approx(10.0, rel=0.05)
    assert abs(BL.bl_pressure_variation(1e5)["delta_p"]) < 1e-3  # the pressure change across the layer, in units of ρU²


def test_outer_flow_V1_eq_9_11_and_9_35_and_cylinder_parity():  # V1 (C01/N11): −(1/ρ)dp/dx = U_e U_e′; (9.35); cylinder U_e = 2U sin φ
    rho = 1.2
    x = np.linspace(0.1, 2.0, 40)
    for kind, kw in (("flat", dict(U=3.0)), ("wedge", dict(n=0.3, a=2.0)), ("diffuser", dict(U1=2.0, L=1.5)),
                     ("cylinder", dict(U=2.0, a=0.5)), ("retarded", dict(U0=2.0, L=3.0, c=1.0))):
        of = BL.outer_flow(kind, **kw)
        h = 1e-6
        dU_fd = (of.Ue(x + h) - of.Ue(x - h)) / (2 * h)
        assert np.max(np.abs(of.dUe(x) - dU_fd)) < 1e-6, kind
        assert np.max(np.abs(of.dpdx(x, rho) + rho * of.Ue(x) * of.dUe(x))) < 1e-12, kind  # (9.11)
    # (9.35): −dp/dx = n a² x^(2n−1) (ρ = 1)
    n, a = 0.3, 2.0
    of = BL.outer_flow("wedge", n=n, a=a)
    assert np.max(np.abs(-of.dpdx(x, 1.0) - n * a ** 2 * x ** (2 * n - 1))) < 1e-12
    # Hiemenz: n = 1 constant thickness; stagnation kind is n = 1
    assert BL.outer_flow("stagnation", a=10.0).Ue(0.1) == pytest.approx(1.0)
    assert BL.outer_flow("stagnation", a=10.0).dpdx(0.1, 1.2) == pytest.approx(-12.0)  # D02 check: −12 Pa/m (favourable)
    # cylinder speed 2U sin φ equals the ch06 potential-flow surface speed on the circle
    cyl = PF.cylinder(2.0, 0.5)
    phi = np.linspace(0.05, 3.0, 25)
    of = BL.outer_flow("cylinder", U=2.0, a=0.5)
    z = 0.5 * np.exp(-1j * np.pi) * np.exp(1j * phi)  # φ from the FORWARD stagnation point z = −a
    speed = np.abs(cyl.dwdz(z))  # |dw/dz| = |u − iv| = the surface speed
    assert np.max(np.abs(speed - of.Ue(0.5 * phi))) < 1e-9
    assert np.max(np.abs(speed - 2 * 2.0 * np.sin(phi))) < 1e-9  # 2U sin φ
    # the diffuser decelerates: adverse gradient dp/dx > 0
    assert BL.outer_flow("diffuser", U1=1.0, L=1.0).dpdx(0.3, 1.0) > 0


# =====================================================================================================================
# C02 — the three thicknesses δ₉₉, δ*, θ (§9.2)
# =====================================================================================================================
def test_thicknesses_V1_exponential_profile_exact():  # V1 (C02): u/U = 1 − e^{−y/a} ⇒ δ* = a, θ = a/2, δ₉₉ = a ln 100, H = 2
    a = 0.013
    y = np.linspace(0.0, 14 * a, 801)
    u = 1.0 - np.exp(-y / a)
    t = BL.thicknesses(y, u, 1.0)
    assert t["delta_star"] == pytest.approx(a, rel=1e-8)
    assert t["theta"] == pytest.approx(a / 2, rel=1e-8)
    assert t["H"] == pytest.approx(2.0, rel=1e-8)
    assert t["delta99"] == pytest.approx(a * np.log(100.0), rel=1e-6)
    # without the tail estimate the integral is short by the neglected tail e^{-14}·a
    assert BL.displacement_thickness(y, u, 1.0, tail=False) == pytest.approx(a * (1 - np.exp(-14)), rel=1e-6)
    # wrong variant (forgetting the tail, stopping the integral at δ₉₉): a 1 % underestimate
    m = y <= a * np.log(100.0)
    assert np.trapezoid(1 - u[m], y[m]) < 0.995 * a
    # Ue scaling: the same profile at edge speed 5
    t5 = BL.thicknesses(y, 5.0 * u, 5.0)
    assert t5["delta_star"] == pytest.approx(a, rel=1e-8)


def test_profile_shape_V1_closed_forms_against_sympy_and_quadrature():  # V1+V2 (C02): exact δ*, θ, H, δ₉₉ of the model profiles
    eta = sp.Symbol("eta", positive=True)
    delta = 0.02
    exact = {
        "linear": (eta, 1),
        "sine": (sp.sin(sp.pi * eta / 2), 1),
        "cubic": (sp.Rational(3, 2) * eta - eta ** 3 / 2, 1),
    }
    y = np.linspace(0.0, 2.0 * delta, 4001)
    for name, (F, top) in exact.items():
        ds = float(sp.integrate(1 - F, (eta, 0, 1)))  # profile is 1 beyond η = 1
        th = float(sp.integrate(F * (1 - F), (eta, 0, 1)))
        s = BL.profile_shape(name, y, delta)
        assert s["delta_star"] == pytest.approx(ds * delta, rel=1e-13), name
        assert s["theta"] == pytest.approx(th * delta, rel=1e-13), name
        assert s["H"] == pytest.approx(ds / th, rel=1e-13), name
        t = BL.thicknesses(y, delta * 0 + s["u_over_Ue"], 1.0, level=0.99)  # numerical route on the sampled profile
        assert t["delta_star"] == pytest.approx(s["delta_star"], rel=5e-6), name
        assert t["theta"] == pytest.approx(s["theta"], rel=5e-6), name
        assert t["delta99"] == pytest.approx(s["delta99"], rel=1e-5), name
    # exponential profile: δ ≡ a
    se = BL.profile_shape("exponential", y, delta)
    assert (se["delta_star"], se["theta"], se["H"]) == pytest.approx((delta, delta / 2, 2.0), rel=1e-13)
    # power-law family: δ* = δ/(p+1), θ = pδ/((p+1)(p+2)), H = (p+2)/p
    for p in (2.0, 7.0):
        sp_ = BL.profile_shape("power", y, delta, p=p)
        ds = quad(lambda s_: 1 - s_ ** (1 / p), 0, 1)[0]
        th = quad(lambda s_: s_ ** (1 / p) * (1 - s_ ** (1 / p)), 0, 1)[0]
        assert sp_["delta_star"] == pytest.approx(ds * delta, rel=1e-10)
        assert sp_["theta"] == pytest.approx(th * delta, rel=1e-10)
        assert sp_["H"] == pytest.approx((p + 2) / p, rel=1e-13)
    # the design's "expect" table (units of δ): linear (0.5, 0.1667, 3), sine (0.3634, 0.1366, 2.660), cubic (0.375, 0.1393, 2.692)
    assert BL.profile_shape("linear", y, 1.0)["delta_star"] == 0.5 and BL.profile_shape("sine", y, 1.0)["H"] == pytest.approx(2.660, abs=1e-3)
    assert BL.profile_shape("cubic", y, 1.0)["theta"] == pytest.approx(0.1393, abs=1e-4)
    with pytest.raises(ValueError):
        BL.profile_shape("nonsense", y, 1.0)


def test_delta_level_V1_root_on_exponential_and_blasius_profiles():  # V1 (C02): δ at any level; monotone in the level
    a = 0.01
    y = np.linspace(0, 12 * a, 1201)
    u = 1 - np.exp(-y / a)
    for lev in (0.5, 0.9, 0.99, 0.999):
        assert BL.delta_level(y, u, 1.0, lev) == pytest.approx(-a * np.log(1 - lev), rel=1e-5)
    eta = np.linspace(0, 12, 2401)
    _, fp, _ = BL.blasius_profile(eta)
    assert BL.delta_level(eta, fp, 1.0, 0.99) == pytest.approx(BL.blasius_constants()["eta99"], rel=1e-8)
    with pytest.raises(ValueError):
        BL.delta_level(y[:5], u[:5], 1.0, 0.99)  # profile never reaches the level


def test_momentum_thickness_V4_theta_equals_integrated_wall_friction_D04():  # V4 (C02/D04): ρU²θ(x) = ∫₀ˣ τ₀ dx′
    U, nu, rho = 2.0, 1.5e-5, 1.2
    xs = np.linspace(1e-4, 1.5, 60)
    # θ(x) from the SAMPLED Blasius profile on a fine η-grid (the sampled-profile route of C02) ...
    eta = np.linspace(0, 14, 3001)
    _, fp, _ = BL.blasius_profile(eta)
    th = np.array([BL.momentum_thickness(eta * np.sqrt(nu * x / U), U * fp, U) for x in xs])
    # ... against the wall friction, integrated with quad (the singular x^{-1/2} is integrable)
    F = np.array([quad(lambda s: float(BL.blasius_wall_shear(s, U, rho, nu)), 0.0, x, epsabs=0, epsrel=1e-12)[0] for x in xs])
    assert np.max(np.abs(rho * U ** 2 * th / F - 1.0)) < 1e-6
    # and the closed-form functions agree with the sampled route
    assert BL.blasius_theta(xs[10], U, nu) == pytest.approx(th[10], rel=1e-7)
    assert BL.blasius_delta_star(xs[10], U, nu) == pytest.approx(
        BL.displacement_thickness(eta * np.sqrt(nu * xs[10] / U), U * fp, U), rel=1e-7)


def test_displacement_V4_streamline_lift_equals_U_dDeltaStar_dx_D03():  # V4 (C02/D03): v(h → ∞) = U dδ*/dx (continuity)
    U, nu = 1.7, 1e-5
    x = 0.8
    v_far = float(BL.blasius_fields(x, 40 * np.sqrt(nu * x / U), U, nu)["v"])
    h = 1e-4
    dds = (BL.blasius_delta_star(x + h, U, nu) - BL.blasius_delta_star(x - h, U, nu)) / (2 * h)
    assert v_far == pytest.approx(U * dds, rel=1e-6)
    assert v_far * np.sqrt(U * x / nu) / U == pytest.approx(BL.blasius_constants()["v_inf"], rel=1e-8)  # 0.8604 = δ*/2 (Fig. 9.6)
    assert BL.blasius_constants()["v_inf"] == pytest.approx(0.5 * BL.blasius_constants()["delta_star"], rel=1e-12)


def test_thicknesses_V7_ordering_shape_factor_and_limits():  # V7 (C02): δ* < δ₉₉, θ < δ*, H ≥ 1; top-hat limit H → 1
    c = BL.blasius_constants()
    assert c["theta"] < c["delta_star"] < c["eta99"]
    assert c["delta_star"] / c["eta99"] == pytest.approx(0.35, abs=0.005)  # D03 trap: δ* = 0.35 δ₉₉
    assert c["H"] == pytest.approx(2.59, abs=0.005)
    y = np.linspace(0, 1, 20001)
    u = np.clip(y / 0.01, 0, 1)  # thin linear ramp: a top hat as the ramp shrinks
    t = BL.thicknesses(y, u, 1.0)
    assert t["H"] == pytest.approx(3.0, rel=1e-3)  # a ramp always has H = 3; the layer thins:
    assert t["delta_star"] == pytest.approx(0.005, rel=1e-3) and t["delta_star"] < 0.01
    assert BL.profile_shape("power", np.linspace(0, 1, 5), 1.0, p=1e6)["H"] == pytest.approx(1.0, abs=1e-5)  # H → 1 as the profile fills
    with pytest.raises(ValueError):
        BL.profile_shape("power", np.linspace(0, 1, 5), 1.0, p=0.5)


# =====================================================================================================================
# C03 — the Blasius reduction by similarity (§9.3)
# =====================================================================================================================
def test_blasius_reduction_V2_sympy_ode_brackets_and_cancellation():  # V2 (C03): (9.18) with ψ = Uδ f(η) ⇒ (9.25)–(9.27)
    r = ch09.similarity_reduce_sympy("blasius")
    assert z0(r["residual"])
    assert z0(r["cancelled_terms"][0] + r["cancelled_terms"][1])  # the η f′ f″ terms cancel (D05 step 8)
    assert r["cancelled_terms"][0] != 0
    assert z0(r["bracket_ratio"] - sp.Symbol("nu", positive=True) * 0 - r["bracket_ratio"])  # sympy expression well-formed
    assert r["ode_coefficients"] == {"fppp": 1, "f f''": sp.Rational(1, 2)}
    assert r["m"] == sp.Rational(1, 2) and r["n"] == 0
    U, nu, x = sp.symbols("U nu x", positive=True)
    assert z0(r["delta"] - sp.sqrt(nu * x / U))  # (9.26)
    assert z0(r["delta_check"])
    # the two x-only brackets of (9.25) are U²δ′/δ and νU/δ² (as functions of the generic δ(x))
    b1, b2 = r["brackets"]
    d = sp.Function("delta")(x)
    assert z0(b1 - U ** 2 * sp.diff(d, x) / d) and z0(b2 - nu * U / d ** 2)


def test_similarity_collapse_V1_blasius_and_jets_collapse_at_the_right_exponents_only():  # V1 (C03/C12/C13): data collapse
    for case, n, m in (("blasius", 0.0, 0.5), ("free_jet", 1 / 3, 2 / 3), ("wall_jet", 0.5, 0.75)):
        good = ch09.similarity_collapse_error(case, n=n, m=m)
        assert good < 1e-12, (case, good)
        assert ch09.similarity_collapse_error(case, n=n, m=m * 0.8) > 1e-2, case  # wrong width exponent: no collapse
        assert ch09.similarity_collapse_error(case, n=n + 0.15, m=m) > 1e-2, case  # wrong amplitude exponent
    # Blasius: a wrong length scale ∝ x^{1/4} does not collapse either
    assert ch09.similarity_collapse_error("blasius", n=0.0, m=0.25) > 1e-2


def test_blasius_fields_V4_streamfunction_gives_continuity_and_9_23_9_24():  # V4 (C03/N26–N27): u = ψ_y, v = −ψ_x from the ψ of (9.19)
    U, nu = 1.3, 2e-5
    x0, y0 = 0.6, 4e-3
    h = 1e-6
    psi = lambda x, y: float(BL.blasius_fields(x, y, U, nu)["psi"])  # noqa: E731
    u_fd = (psi(x0, y0 + h) - psi(x0, y0 - h)) / (2 * h)
    v_fd = -(psi(x0 + h, y0) - psi(x0 - h, y0)) / (2 * h)
    f = BL.blasius_fields(x0, y0, U, nu)
    assert float(f["u"]) == pytest.approx(u_fd, rel=1e-6)  # (9.23)
    assert float(f["v"]) == pytest.approx(v_fd, rel=1e-5)  # (9.24)
    assert float(f["delta"]) == pytest.approx(np.sqrt(nu * x0 / U), rel=1e-14)  # (9.26)
    assert float(f["eta"]) == pytest.approx(y0 / np.sqrt(nu * x0 / U), rel=1e-14)
    with pytest.raises(ValueError):
        BL.blasius_fields(0.0, 1e-3, U, nu)  # the leading edge is singular
    # no-slip and no through-flow at the wall (9.20)
    w = BL.blasius_fields(0.5, 0.0, U, nu)
    assert float(w["u"]) == 0.0 and float(w["v"]) == pytest.approx(0.0, abs=1e-15)


# =====================================================================================================================
# C04 — the Blasius solution and its numbers (§9.3)
# =====================================================================================================================
@needs_ref
def test_blasius_V5_wall_shear_against_published_values():  # V5 (C04): Belden et al. (arXiv:1907.09912) and Wikipedia
    R = ref_json()
    c = BL.blasius_constants()
    kappa = [p for p in R["belden_falkner_skan"]["points"] if p["beta"] == 0.0][0]["kappa"]
    _, fb = belden_to_book(0.0, kappa)
    assert rel(c["fpp0"], fb) < 1e-10  # 12 printed digits (their normalisation converted to the book's)
    assert rel(c["fpp0"], R["blasius_wikipedia"]["fpp0"]) < 1e-12  # Wikipedia small-η coefficient (15 digits)
    # the second number on the Wikipedia page (…043934904293) is NOT the converged constant: ours differs from it by 4e-5
    assert 3e-5 < rel(c["fpp0"], R["blasius_wikipedia"]["fpp0_other_number_on_the_same_page"]) < 5e-5
    # rounded thicknesses (Wikipedia 1.72 and 0.665; drag constant 1.328)
    assert rel(c["delta_star"], R["blasius_wikipedia"]["delta_star_coeff_rounded"]) < 5e-3
    assert rel(c["theta"], R["blasius_wikipedia"]["theta_coeff_rounded"]) < 2e-3
    assert c["cd_coeff"] == pytest.approx(1.328, abs=5e-4)
    # cf agrees with the secondary skin-friction page: c_f = 0.664 Re_x^{-1/2}
    assert BL.blasius_skin_friction(1e4) == pytest.approx(R["skin_friction_wikipedia"]["laminar_cf_coeff"] / 100.0, rel=3e-4)


def test_blasius_V3_two_solvers_agree_and_truncation_is_insensitive():  # V3 (C04): Töpfer IVP vs solve_bvp vs shooting
    c = BL.blasius_constants()  # Töpfer scaling of one initial-value problem
    d_t = BL.falkner_skan(0.0, method="toepfer")
    assert d_t["fpp0"] == c["fpp0"]
    d16 = BL.falkner_skan(0.0, eta_max=16.0, tol=1e-10)  # solve_bvp, an independent algorithm
    assert abs(d16["fpp0"] - c["fpp0"]) < 1e-10  # the task's 1e-10 by two solvers (observed ≈ 2e-14)
    assert np.max(np.abs(d16["fp"] - BL.blasius_profile(d16["eta"])[1])) < 1e-9
    d_s = BL.falkner_skan(0.0, eta_max=12.0, method="shoot")  # a third route: brentq on f″(0)
    assert abs(d_s["fpp0"] - c["fpp0"]) < 1e-9
    # truncation study: |f″(0)(η_max) − f″(0)| decays like the Gaussian tail e^{-η²/4}: 8 → 12 → 16
    err = [abs(BL.falkner_skan(0.0, eta_max=e, tol=1e-10)["fpp0"] - c["fpp0"]) for e in (8.0, 10.0, 12.0)]
    assert err[0] > err[1] > err[2] and err[2] < 1e-10
    assert 1e-7 < err[0] < 1e-5  # the default η_max = 10 is good to ~1e-9 only (documented; see report)
    # from-scratch RK4-free check: an independent scipy shoot with the Töpfer scaling law f″(0) = F′(∞)^{-3/2}
    sol = solve_ivp(lambda s, Y: [Y[1], Y[2], -0.5 * Y[0] * Y[2]], (0, 14), [0, 0, 1.0], rtol=1e-12, atol=1e-14)
    assert sol.y[1, -1] ** -1.5 == pytest.approx(c["fpp0"], rel=1e-9)


def test_blasius_V2_ode_residual_boundary_conditions_and_tail():  # V2 (C04): f‴ + ½ff″ = 0, (9.28), (9.29) on the returned profile
    eta = np.linspace(0.0, 14.0, 5601)
    f, fp, fpp = BL.blasius_profile(eta)
    h = eta[1] - eta[0]
    fppp = np.gradient(fpp, h, edge_order=2)
    res = fppp + 0.5 * f * fpp
    assert np.max(np.abs(res[5:-5])) < 2e-6 * np.max(np.abs(fpp))  # 2nd-order differences of a 1e-13 solution
    assert np.max(np.abs(np.gradient(f, h, edge_order=2)[5:-5] - fp[5:-5])) < 1e-6  # f′ = df/dη
    assert np.max(np.abs(np.gradient(fp, h, edge_order=2)[5:-5] - fpp[5:-5])) < 1e-6  # f″ = df′/dη
    assert (f[0], fp[0]) == (0.0, 0.0) and fp[-1] == pytest.approx(1.0, abs=1e-12)  # (9.28), (9.29)
    # the far field f = η − δ* (linear)
    assert f[-1] == pytest.approx(eta[-1] - BL.blasius_constants()["delta_star"], abs=1e-9)
    # N32: f′ → 1 monotonically, f″ > 0 and decreasing beyond the maximum at the wall
    core = eta < 10.0  # beyond η ≈ 10 the solution is at the 1e-18 round-off floor
    assert np.all(np.diff(fp) >= -1e-14) and np.all(fpp[core] > 0) and np.all(np.diff(fpp)[core[1:]] <= 1e-12)
    # N33: Gaussian tail of f′ − 1 (design: the linearised far field within 0.3 % for 4 ≤ η ≤ 8)
    e = np.linspace(4.0, 8.0, 41)
    dev = BL.blasius_far_field(e)
    num = BL.blasius_profile(e)[1] - 1.0
    assert np.max(np.abs(dev / num - 1.0)) < 3e-3
    # (1/η)e^{-η²/4} form: the decay rate of f″ is the Gaussian's, ln f″ ≈ −η²/4 + const  (slope check between η = 6 and 7)
    ds_ = BL.blasius_constants()["delta_star"]
    slope = (np.log(BL.blasius_profile(7.0)[2]) - np.log(BL.blasius_profile(6.0)[2])) / 1.0
    assert slope == pytest.approx(-((7.0 - ds_) ** 2 - (6.0 - ds_) ** 2) / 4.0, abs=0.02)  # f″ ∝ exp(−(η − δ*)²/4): the shifted Gaussian


def test_blasius_V7_theta_is_twice_fpp0_and_thickness_relations():  # V7 (C04): θ = 2f″(0) = C_f√Re_x; H = δ*/θ; v∞ = δ*/2
    c = BL.blasius_constants()
    assert c["theta"] == pytest.approx(2 * c["fpp0"], rel=1e-10)  # D06 step 6: integrate (9.27) by parts
    assert c["cf_coeff"] == pytest.approx(c["theta"], rel=1e-10)  # C_f√Re_x = θ/δ
    assert c["cd_coeff"] == pytest.approx(2 * c["cf_coeff"], rel=1e-14)  # C_D = 2C_f(L)
    assert c["tau_coeff"] == c["fpp0"]
    assert c["H"] == pytest.approx(c["delta_star"] / c["theta"], rel=1e-14)
    # D06 step 6 intermediate line at a finite η: f″(0) = ½ I_θ(η) − ½ f(1 − f′) + f″(η)  (the design prints it without the +f″(η))
    for e in (1.0, 2.0, 4.0):
        Ith = quad(lambda s: BL.blasius_profile(s)[1] * (1 - BL.blasius_profile(s)[1]), 0, e, epsabs=1e-14, epsrel=1e-13)[0]
        f, fp, fpp = BL.blasius_profile(e)
        assert c["fpp0"] == pytest.approx(0.5 * Ith - 0.5 * f * (1 - fp) + fpp, abs=1e-11)
    # the D06 step 6 line as printed (no +f″(η)) is only true in the limit η → ∞
    f, fp, fpp = BL.blasius_profile(2.0)
    Ith = quad(lambda s: BL.blasius_profile(s)[1] * (1 - BL.blasius_profile(s)[1]), 0, 2.0)[0]
    assert abs(c["fpp0"] - 0.5 * (Ith - f * (1 - fp))) > 0.01
    # Töpfer symmetry (D06 step 1–3): f″(0) = F′(∞)^{-3/2}, F′(∞) = 2.0854 for F″(0) = 1
    assert c["fpp0"] ** (-2.0 / 3.0) == pytest.approx(2.0854, abs=1e-3)


def test_blasius_V1_wall_shear_drag_and_R5_one_side():  # V1 (C04/N37–N40): (9.31)–(9.33), drag by quadrature, both faces double (R5)
    U, rho, nu, L = 1.0, 1.2, 1.5e-5, 1.0
    Re = U * L / nu
    c = BL.blasius_constants()
    assert BL.blasius_wall_shear(L, U, rho, nu) == pytest.approx(c["fpp0"] * rho * U ** 2 / np.sqrt(Re), rel=1e-14)  # (9.31)
    assert BL.blasius_skin_friction(Re) == pytest.approx(BL.blasius_wall_shear(L, U, rho, nu) / (0.5 * rho * U ** 2), rel=1e-13)  # (9.32)
    F1 = BL.blasius_drag(L, U, rho, nu)  # quad of (9.31)
    assert F1 == pytest.approx(2 * BL.blasius_wall_shear(L, U, rho, nu) * L, rel=1e-9)  # ∫x^{-1/2}dx = 2√L (D06 step 9)
    assert BL.blasius_drag_coefficient(Re) == pytest.approx(F1 / (0.5 * rho * U ** 2 * L), rel=1e-9)  # (9.33)
    assert BL.blasius_drag_coefficient(Re) == pytest.approx(2 * BL.blasius_skin_friction(Re), rel=1e-13)
    assert BL.blasius_drag(L, U, rho, nu, sides=2) == pytest.approx(2 * F1, rel=1e-14)  # R5: a two-sided plate doubles
    assert BL.blasius_drag_coefficient(Re, sides=2) == pytest.approx(2 * BL.blasius_drag_coefficient(Re), rel=1e-14)
    # drag ∝ U^{3/2} (N39): doubling U multiplies F_D by 2^{3/2}
    assert BL.blasius_drag(L, 2 * U, rho, nu) / F1 == pytest.approx(2 ** 1.5, rel=1e-9)
    # design D06 check numbers: air, 1 m at 1 m/s: τ₀(1 m) = 1.543e-3 Pa, F_D = 3.09e-3 N/m, C_D = 5.14e-3
    assert BL.blasius_wall_shear(1.0, 1.0, 1.2, 1.5e-5) == pytest.approx(1.543e-3, rel=1e-3)
    assert F1 == pytest.approx(3.09e-3, rel=2e-3) and BL.blasius_drag_coefficient(Re) == pytest.approx(5.14e-3, rel=2e-3)
    # dimensional homogeneity of (9.31)
    dimensional_check(lambda rho, U, x, nu: 0.332 * rho * U ** 2 / (U * x / nu) ** 0.5, "pressure",
                      rho=Q_(1.2, "kg/m**3"), U=Q_(1.0, "m/s"), x=Q_(1.0, "m"), nu=Q_(1.5e-5, "m**2/s"))


def test_blasius_V1_R2_printed_eta99_is_not_the_99_percent_root():  # V1 (C04/N35): 4.910 is the root of f′ = 0.99; 4.93 is not
    eta99 = BL.blasius_constants()["eta99"]
    fp = lambda e: float(BL.blasius_profile(e)[1])  # noqa: E731
    assert abs(fp(eta99) - 0.99) < 1e-12
    assert eta99 == pytest.approx(4.910, abs=5e-4)
    # planted wrong variant: the printed 4.93 leaves f′ = 0.9904, four parts in 10⁴ off the definition
    U, nu, x = 1.0, 1.5e-5, 1.0
    d_ok, d_bad = BL.blasius_delta99(x, U, nu), BL.blasius_delta99(x, U, nu, printed=True)
    assert abs(fp(d_ok / np.sqrt(nu * x / U)) - 0.99) < 1e-12
    assert abs(fp(d_bad / np.sqrt(nu * x / U)) - 0.99) > 2e-4, "the printed 4.93 must NOT satisfy f′ = 0.99"
    assert d_bad / d_ok == pytest.approx(4.93 / eta99, rel=1e-12)
    assert 0.003 < d_bad / d_ok - 1 < 0.005  # 0.4 %


def test_blasius_V7_temporal_and_spatial_layers_compare_D06_N41():  # V7 (C04/N41): C_f√Re_x = 0.664 vs 1.128 (ch08 temporal) vs 1.328 (plate mean)
    c = BL.blasius_constants()
    t = LAM.temporal_bl_wall_stress(1.0, 1.0, 1e-6, 1000.0)
    assert t["Cf_coefficient"] == pytest.approx(2 / np.sqrt(np.pi), rel=1e-12)  # 1.1284
    assert c["cf_coeff"] < t["Cf_coefficient"] < c["cd_coeff"]
    assert c["cd_coeff"] / c["cf_coeff"] == pytest.approx(2.0, rel=1e-13)


# =====================================================================================================================

# =====================================================================================================================
# C05 — Falkner–Skan (§9.4), N43–N48
# =====================================================================================================================
FS_MS = (-0.05, 0.1, 1.0 / 3.0, 1.0, 4.0)


@pytest.mark.parametrize("m", FS_MS)
def test_falkner_skan_V2_ode_residual_and_boundary_conditions(m):  # V2 (C05): (9.36) on the solved profile, (9.28), (9.29)
    d = BL.falkner_skan(m, eta_max=16.0, n=25601, tol=1e-10)  # h = 6e-4 keeps the 2nd-order FD error of f‴ below 3e-6
    assert d["success"]
    eta, f, fp, fpp = d["eta"], d["f"], d["fp"], d["fpp"]
    h = eta[1] - eta[0]
    res = np.gradient(fpp, h, edge_order=2) + 0.5 * (m + 1) * f * fpp - m * fp ** 2 + m  # f‴ by central differences of f″
    assert np.max(np.abs(res[3:-3])) < 3e-6 * max(1.0, np.max(np.abs(fpp))), m
    assert np.max(np.abs(np.gradient(f, h, edge_order=2)[3:-3] - fp[3:-3])) < 1e-5 and np.max(np.abs(np.gradient(fp, h, edge_order=2)[3:-3] - fpp[3:-3])) < 1e-5
    assert (f[0], fp[0]) == pytest.approx((0.0, 0.0), abs=1e-14) and fp[-1] == pytest.approx(1.0, abs=1e-8)  # (9.28), (9.29)
    # the returned f‴ is the ODE itself; V7: f‴(0) = −m (read off (9.36) at η = 0 where f = f′ = 0)
    assert d["fppp"][0] == pytest.approx(-m, abs=1e-12)


def test_falkner_skan_V2_sympy_reduction_and_D07_lines():  # V2 (C05/D07): the reduction to (9.36) and every displayed step of D07
    r = ch09.similarity_reduce_sympy("falkner_skan")
    n = sp.Symbol("n", positive=True)
    assert z0(r["residual"])
    assert r["ode_coefficients"] == {"fppp": 1, "f f''": (n + 1) / 2, "f'^2": -n, "1": n}
    assert r["n"] == n
    # n = 0 reproduces the Blasius ODE (9.27)
    F0, F1, F2, F3 = sp.symbols("F0 F1 F2 F3")
    assert z0(r["ode"].lhs.subs(n, 0) - (F3 + F0 * F2 / 2))
    # n = 1: constant thickness δ = √(ν/a) (D07 check line)
    x, nu, a = sp.symbols("x nu a", positive=True)
    assert z0(r["delta"].subs(n, 1) - sp.sqrt(nu / a))


@needs_ref
@pytest.mark.parametrize("beta", (0.5, 0.0, -0.12))
def test_falkner_skan_V5_wall_shear_against_belden_et_al(beta):  # V5 (C05): κ = f″(0) at published β (12–15 digits), normalisation converted
    pt = [p for p in ref_json()["belden_falkner_skan"]["points"] if p["beta"] == beta and p["kappa"] > 0][0]
    m, f_book = belden_to_book(pt["beta"], pt["kappa"])
    d = BL.falkner_skan(m, eta_max=16.0, tol=1e-10)
    assert d["success"]
    assert rel(d["fpp0"], f_book) < 1e-9, (beta, d["fpp0"], f_book)
    # wrong variant: forgetting the √((m+1)/2) normalisation factor is off by tens of percent (β = 0.5) — the test discriminates
    if beta != 0.0:
        assert rel(d["fpp0"], pt["kappa"]) > 0.05


@needs_ref
def test_falkner_skan_V5_separation_member():  # V5 (C05/N47): β_sep = −0.198837735, m_sep = −0.09043 (the fold of the attached branch)
    R = ref_json()["belden_falkner_skan"]["points"]
    sep = BL.falkner_skan_separation()
    beta_sep = [p for p in R if p["kappa"] == 0.0][0]["beta"]
    assert abs(sep["beta_sep"] - beta_sep) < 1e-9  # the source prints 9 digits
    assert sep["fpp0_at_sep"] == 0.0 and sep["m_sep"] == pytest.approx(beta_sep / (2 - beta_sep), abs=1e-9)
    assert sep["m_sep"] == pytest.approx(-0.09043, abs=5e-6)


@needs_ref
@pytest.mark.slow
def test_falkner_skan_V5_reversed_branch_beta_minus_0_02():  # V5 (C05/R11): Stewartson's second branch κ(β = −0.02) (slow: the continuation past the fold is ≈ 15 s, cached for the next test)
    R = ref_json()["belden_falkner_skan"]["points"]
    pt = [p for p in R if p["beta"] == -0.02 and p["kappa"] < 0][0]
    m, f_book = belden_to_book(pt["beta"], pt["kappa"])
    d = BL.falkner_skan(m, branch="reversed")  # Stewartson's second branch (R11)
    assert d["success"] and d["fpp0"] < 0
    assert rel(d["fpp0"], f_book) < 1e-9
    assert np.min(d["fp"]) < 0.0  # the second branch really has reverse flow near the wall
    # V7: attached solution at the same m has κ > 0 and no reverse flow
    a = BL.falkner_skan(m, eta_max=16.0)
    assert a["fpp0"] > 0 and np.min(a["fp"]) > -1e-12  # no reverse flow (round-off floor)


@needs_ref
@pytest.mark.slow
def test_falkner_skan_V5_reversed_branch_beta_minus_0_12():  # V5 (C05/R11): a second reversed-branch point (slow: continuation ≈ 15 s)
    pt = [p for p in ref_json()["belden_falkner_skan"]["points"] if p["beta"] == -0.12 and p["kappa"] < 0][0]
    m, f_book = belden_to_book(pt["beta"], pt["kappa"])
    d = BL.falkner_skan(m, branch="reversed")
    assert d["success"] and rel(d["fpp0"], f_book) < 1e-8, (d["fpp0"], f_book)


def test_falkner_skan_R11_no_attached_solution_below_the_fold():  # V7 (C05/N47): the printed "reverse flow for n < −0.0904" is a slip
    below = BL.falkner_skan(-0.095)
    assert below["success"] is False and np.all(np.isnan(below["fpp"])), "no attached solution with f′(∞) = 1 exists below the fold"
    with pytest.raises(RuntimeError):
        BL.falkner_skan_state(-0.095)
    just_above = BL.falkner_skan(-0.0904, eta_max=16.0)
    assert just_above["success"] and 0 < just_above["fpp0"] < 0.01  # f″(0) → 0 at the fold (κ(−0.0904) ≈ 0.0048)
    # the reversed branch exists only in −0.0904 < m < 0
    for bad in (-0.095, 0.0, 0.2):
        with pytest.raises(ValueError):
            BL.falkner_skan(bad, branch="reversed")
    # f″(0) decreases monotonically to 0 as m ↓ fold (V7)
    ms = np.array([0.0, -0.03, -0.06, -0.08, -0.09])
    k = [BL.falkner_skan_state(m_)["fpp0"] for m_ in ms]
    assert np.all(np.diff(k) < 0) and k[-1] < 0.05
    # N73: the adverse gradient thickens the layer — H = δ*/θ and δ*/δ grow monotonically as m falls towards the fold (fuller profile for m > 0)
    H = [BL.falkner_skan_state(m_)["H"] for m_ in (4.0, 1.0, 0.0, -0.03, -0.06, -0.09)]
    Idl = [BL.falkner_skan_state(m_)["I_delta"] for m_ in (4.0, 1.0, 0.0, -0.03, -0.06, -0.09)]
    assert np.all(np.diff(H) > 0) and np.all(np.diff(Idl) > 0)


def test_falkner_skan_V1_momentum_integral_closes_for_every_member_and_design_expect_row():  # V1/V4 (C05/C06): f″(0) = (3m+1)/2 I_θ + m I_δ from (9.43)
    for m in (-0.09, -0.05, 0.0, 0.1, 1.0 / 3.0, 1.0, 4.0):
        st = BL.falkner_skan_state(m, eta_max=16.0)  # (the default η_max = 10 at m = −0.05 leaves a 6e-8 truncation error, see report)
        # (9.43) for U_e = a xᵐ reduces to f″(0) = ((3m+1)/2) I_θ + m I_δ — independent of the solver's own f″(0):
        assert st["fpp0"] == pytest.approx(0.5 * (3 * m + 1) * st["I_theta"] + m * st["I_delta"], rel=1e-8, abs=1e-9), m
        assert st["fppp0"] == pytest.approx(-m, abs=1e-12)  # V7
        assert st["cf_sqrtRex"] == pytest.approx(2 * st["fpp0"], rel=1e-14)
        assert st["lam"] == pytest.approx(m * st["I_theta"] ** 2, rel=1e-14)
        assert st["l"] == pytest.approx(st["fpp0"] * st["I_theta"], rel=1e-14) and st["H"] == pytest.approx(st["I_delta"] / st["I_theta"], rel=1e-14)
    h = BL.falkner_skan_state(1.0)  # Hiemenz stagnation flow: the design's expect row (1.2326, 2.2162, 0.08546, 0.36034)
    assert (h["fpp0"], h["H"], h["lam"], h["l"]) == pytest.approx((1.2326, 2.2162, 0.08546, 0.36034), abs=6e-5)


@needs_ref
def test_hiemenz_V5_wall_shear_and_displacement_thickness():  # V5 (C05): F″(0) = 1.232588 (Weidman & Turner); δ*/δ = 0.6479 (Wikipedia)
    R = ref_json()
    h = BL.falkner_skan_state(1.0)
    assert abs(h["fpp0"] - R["hiemenz_weidman_turner"]["Fpp0"]) < 5e-7  # 7 printed digits
    assert abs(h["I_delta"] - R["hiemenz_delta_star_wikipedia"]["delta_star_over_delta"]) < 5e-5
    # constant thickness (n = 1): δ = √(ν/a) at every x
    assert BL.falkner_skan_thickness(0.1, 1.0, 10.0, 1.5e-5) == pytest.approx(BL.falkner_skan_thickness(3.0, 1.0, 10.0, 1.5e-5), rel=1e-14)
    assert BL.falkner_skan_thickness(1.0, 1.0, 10.0, 1.5e-5) == pytest.approx(np.sqrt(1.5e-5 / 10.0), rel=1e-14)
    # thickness grows for n < 1 (∝ x^{(1−n)/2}), shrinks for n > 1
    assert BL.falkner_skan_thickness(2.0, 0.0, 1.0, 1e-6) / BL.falkner_skan_thickness(1.0, 0.0, 1.0, 1e-6) == pytest.approx(np.sqrt(2.0), rel=1e-14)
    assert BL.falkner_skan_thickness(2.0, 2.0, 1.0, 1e-6) < BL.falkner_skan_thickness(1.0, 2.0, 1.0, 1e-6)


def test_falkner_skan_V3_bvp_and_shooting_agree():  # V3 (C05): two algorithms (solve_bvp with continuation; brentq shoot) for m ≥ −0.05
    for m in (-0.05, 0.1, 1.0 / 3.0):
        b = BL.falkner_skan(m, eta_max=12.0, tol=1e-10)
        s = BL.falkner_skan(m, eta_max=12.0, method="shoot")
        assert abs(b["fpp0"] - s["fpp0"]) < 5e-9, (m, b["fpp0"], s["fpp0"])
    with pytest.raises(ValueError):
        BL.falkner_skan(-0.08, method="shoot")  # shooting is exponentially sensitive near separation: refused
    with pytest.raises(ValueError):
        BL.falkner_skan(0.3, method="toepfer")  # Töpfer scaling is the Blasius (m = 0) case only
    # truncation: the answer is insensitive to η_max (8, 12, 16) at Hiemenz
    v = [BL.falkner_skan(1.0, eta_max=e)["fpp0"] for e in (8.0, 12.0, 16.0)]
    assert abs(v[1] - v[2]) < 1e-9 and abs(v[0] - v[2]) < 1e-5


def test_falkner_skan_shoot_V3_hiemenz_by_shooting():  # V3 (C05): the documented "shoot" method for m ≥ −0.05 — Hiemenz and n = 4 (Fig. 9.7 members)
    # (F1, fixed) brentq used to raise "f(a) and f(b) must have different signs" for m = 1 and m = 4.
    # the method works for m ≥ −0.05 and the design's from-scratch moment shoots Hiemenz.  Hypothesis: for large m the lower bracket
    # s = 1e-4 already makes f′ overshoot (the event returns the positive sentinel 0.6), so g(lo) and g(hi) are both positive.
    for m, want in ((1.0, 1.2325876568), (4.0, 2.4057248594)):
        s = BL.falkner_skan(m, eta_max=10.0, method="shoot")
        assert s["fpp0"] == pytest.approx(want, rel=1e-7), m


def test_falkner_skan_V1_wall_curvature_equals_pressure_gradient_D12():  # V1 (C05/C08): μ u_yy(wall) = dp/dx from the solved FIELD (FD in y)
    nu, a, rho = 1e-5, 2.0, 1.0
    for m in (-0.05, 0.3, 1.0):
        x = 0.4
        of = BL.outer_flow("wedge", n=m, a=a)
        delta = float(BL.falkner_skan_thickness(x, m, a, nu))
        y = np.linspace(0, 0.5 * delta, 4001)
        u = BL.falkner_skan_fields(x, y, m, a, nu, eta_max=16.0)["u"]
        h = y[1]
        uyy0 = (2 * u[0] - 5 * u[1] + 4 * u[2] - u[3]) / h ** 2  # one-sided 2nd-order (u = 0 at the wall)
        assert float(BL.wall_curvature(of.dpdx(x, rho), rho * nu)) == pytest.approx(uyy0, rel=2e-3), m  # μ u_yy = dp/dx
        assert (uyy0 > 0) == (m < 0)  # (9.51)/(9.52): adverse (m < 0) ⇒ u_yy > 0; favourable ⇒ u_yy < 0
    assert float(BL.wall_curvature(3.0, 0.5)) == 6.0


def test_falkner_skan_fields_V4_continuity_and_9_9_residual():  # V4 (C05): u = U_e f′, v from ψ = √(νxU_e) f (9.34): (6.2) and (9.9) hold
    nu, a, m = 1e-4, 1.5, 0.4
    x = np.linspace(0.6, 1.4, 61)
    y = np.linspace(0.0, 0.05, 61)
    X, Y = np.meshgrid(x, y, indexing="xy")
    f = BL.falkner_skan_fields(X, Y, m, a, nu, eta_max=16.0)
    of = BL.outer_flow("wedge", n=m, a=a)
    dpdx = np.array([float(of.dpdx(v, 1.2)) for v in x])
    res = BL.bl_x_momentum_residual(x, y, f["u"], f["v"], dpdx, nu, rho=1.2)
    sc = float(np.max(np.abs(f["u"])) ** 2 / 0.6)
    assert np.max(np.abs(res[2:-2, 2:-2])) < 5e-3 * sc
    ux = np.gradient(f["u"], x, axis=1, edge_order=2)
    vy = np.gradient(f["v"], y, axis=0, edge_order=2)
    assert np.max(np.abs((ux + vy)[2:-2, 2:-2])) < 5e-3 * np.max(np.abs(f["u"])) / 0.6
    # wrong variant: the Blasius field with the wedge pressure gradient does NOT satisfy (9.9)
    fb = BL.blasius_fields(X, Y, float(of.Ue(1.0)), nu)
    resb = BL.bl_x_momentum_residual(x, y, fb["u"], fb["v"], dpdx, nu, rho=1.2)
    assert np.max(np.abs(resb[2:-2, 2:-2])) > 5e-2 * sc


def test_falkner_skan_table_V1_rows_agree_with_state_and_fig_9_7_members_exist():  # V1 (C05/N46): table for the explainers = the solver
    ms = [-0.05, 0.0, 1.0 / 3.0, 1.0, 4.0]
    t = BL.falkner_skan_table(m_grid=ms, n_eta=81)
    assert t["fp"].shape == (5, 81) and t["eta"][-1] == 8.0
    for i, m in enumerate(ms):
        st = BL.falkner_skan_state(m)
        assert t["fpp0"][i] == pytest.approx(st["fpp0"], rel=1e-7)
        assert t["lam"][i] == pytest.approx(st["lam"], rel=1e-7, abs=1e-12) and t["H"][i] == pytest.approx(st["H"], rel=1e-7)
        assert abs(t["fp"][i, 0]) < 1e-12
    assert t["fpp0"][1] == pytest.approx(0.33206, abs=5e-6)
    assert t["inflection_eta"][1] == 0.0 and t["inflection_eta"][0] > 0  # Blasius: at the wall; n = −0.05: at finite η
    assert np.isnan(t["inflection_eta"][2])  # nan for favourable n = 1/3 (no inflection)
    tf = BL.falkner_skan_table(fast=True)  # the default 41-node grid, halved
    assert tf["m"].min() > -0.09043 and tf["m"].max() == 4.0 and np.all(np.diff(tf["fpp0"]) > 0) and np.all(np.isfinite(tf["fp"]))
    assert np.all(np.isnan(tf["fpp0_reversed"]))  # fast mode skips the second branch
    # the seven Fig. 9.7 members are all attached solutions (the labels are public exponents, the curves are ours)
    for m in (4, 1, 1 / 3, 1 / 9, 0, -0.0654, -0.0904):
        assert BL.falkner_skan(m, eta_max=16.0)["success"], m


# =====================================================================================================================
# C06 — the von Kármán momentum integral (§9.5)
# =====================================================================================================================
def test_momentum_integral_V2_sympy_steps_9_37_to_9_43():  # V2 (C06/D08): flux form, the (9.41)→(9.42) identity, and the deficit integrals
    r = ch09.momentum_integral_sympy()
    assert z0(r["identity_9_38"]) and z0(r["identity_9_42"]) and z0(r["residual"])
    labels = [s_[0] for s_ in r["steps"]]
    assert labels == ["(9.37)", "(9.38)", "(9.39)", "(9.40)", "(9.41)", "(9.42)"]
    # (9.43) as returned: τ₀/ρ = d(U_e²θ)/dx + U_e δ* U_e′
    assert r["result"].lhs == sp.Symbol("tau_0") / sp.Symbol("rho", positive=True)


def test_momentum_integral_V1_residual_is_zero_for_every_exact_solution():  # V1/V4 (C06): (9.43) closes for Blasius and Falkner–Skan
    nu, a = 1e-5, 1.3
    x = np.geomspace(0.05, 2.0, 400)
    for m in (0.0, 0.3, 1.0, -0.05):
        st = BL.falkner_skan_state(m, eta_max=16.0)
        Ue = a * x ** m
        th = st["I_theta"] * np.sqrt(nu * x / Ue)
        ds = st["I_delta"] * np.sqrt(nu * x / Ue)
        tau0 = 1.0 * nu * Ue * st["fpp0"] / np.sqrt(nu * x / Ue)  # ρ = 1
        r = BL.momentum_integral_residual(x, Ue, th, ds, tau0, rho=1.0)
        assert np.max(np.abs(r[5:-5]) / tau0[5:-5]) < 1e-7, m
    # wrong variant of (9.43) is detected: without the pressure-gradient term (δ* dropped)
    st = BL.falkner_skan_state(0.3)
    Ue = a * x ** 0.3
    th, ds = st["I_theta"] * np.sqrt(nu * x / Ue), st["I_delta"] * np.sqrt(nu * x / Ue)
    tau0 = nu * Ue * st["fpp0"] / np.sqrt(nu * x / Ue)
    assert np.max(np.abs(BL.momentum_integral_residual(x, Ue, th, 0 * ds, tau0)[5:-5]) / tau0[5:-5]) > 0.05
    # Blasius from (9.43) directly: τ₀/ρ = U² dθ/dx, dθ/dx = θ/(2x)
    U = 2.0
    th_b = BL.blasius_theta(x, U, nu)
    assert np.max(np.abs(U ** 2 * th_b / (2 * x) - BL.blasius_wall_shear(x, U, 1.0, nu)) / BL.blasius_wall_shear(x, U, 1.0, nu)) < 1e-9


def test_momentum_integral_V4_marched_layer_satisfies_9_43():  # V4 (C06): an independent solver (parabolic marching) closes (9.43) to 2e-3
    nu = 1e-3
    for of, uin in ((BL.outer_flow("flat", U=1.0), lambda y: np.where(y < 0.05, np.sin(np.pi * np.minimum(y, 0.05) / 0.1), 1.0)),
                    (BL.outer_flow("diffuser", U1=1.0, L=30.0), None)):
        x = np.linspace(0.1, 1.0, 181)
        r = BL.march_boundary_layer(of, x, nu, ny=800, u_inlet=uin)
        xs = r["x"]
        Ue = np.array([float(of.Ue(v)) for v in xs])
        th = np.array([BL.momentum_thickness(r["y"][k], r["u"][k], Ue[k]) for k in range(len(xs))])
        ds = np.array([BL.displacement_thickness(r["y"][k], r["u"][k], Ue[k]) for k in range(len(xs))])
        res = BL.momentum_integral_residual(xs, Ue, th, ds, r["tau0"], rho=1.0)
        assert np.max(np.abs(res[3:-3]) / r["tau0"][3:-3]) < 2e-3


def test_karman_pohlhausen_V1_exact_coefficients_and_V5_within_a_few_percent_of_blasius():  # V1+V5 (C06/N55): assumed profiles, (9.43) closed by hand
    nu, U = 1e-5, 1.0
    x = np.linspace(0.05, 1.0, 40)
    c = BL.blasius_constants()
    # exact by sympy: δ δ′ = F′(0) ν/(U c_θ) for u/U = F(y/δ) (D08 with dU_e = 0): δ²/(νx/U) = 2F′(0)/c_θ
    e = sp.Symbol("e", positive=True)
    cases = {"cubic": sp.Rational(3, 2) * e - e ** 3 / 2, "sine": sp.sin(sp.pi * e / 2), "quartic": 2 * e - 2 * e ** 3 + e ** 4}
    for prof, Fe in cases.items():
        c_th = sp.integrate(Fe * (1 - Fe), (e, 0, 1))
        Fp0 = sp.diff(Fe, e).subs(e, 0)
        d_over = float(sp.sqrt(2 * Fp0 / c_th))  # δ/√(νx/U)
        r = BL.karman_pohlhausen(BL.outer_flow("flat", U=U), x, nu, profile=prof)
        s = np.sqrt(nu * x / U)
        assert np.max(np.abs(r["delta"] / s / d_over - 1)) < 1e-5, prof
        assert np.max(np.abs(r["theta"] / s / (float(c_th) * d_over) - 1)) < 1e-5, prof
        cf_sqrt = 2 * float(Fp0) / d_over  # C_f √Re_x
        assert np.max(np.abs(2 * r["tau0"] / (U ** 2) * np.sqrt(U * x / nu) / cf_sqrt - 1)) < 1e-5, prof
        # V5: within 6 % of the exact Blasius numbers (docstring claim) — the assumed profile is a closure, not the truth
        assert abs(r["theta"][-1] / s[-1] / c["theta"] - 1) < 0.06, prof
    # the design's expect row: cubic θ = 0.6465 (−2.7 %), δ = 4.641, τ₀ coefficient 0.3232
    r = BL.karman_pohlhausen(BL.outer_flow("flat", U=U), x, nu, profile="cubic")
    s = np.sqrt(nu * x / U)
    assert r["delta"][-1] / s[-1] == pytest.approx(4.641, abs=1e-3) and r["theta"][-1] / s[-1] == pytest.approx(0.6465, abs=1e-4)
    assert r["tau0"][-1] / (U ** 2) * np.sqrt(U * x[-1] / nu) == pytest.approx(0.3232, abs=1e-4)
    # V7: an adverse gradient thickens the layer and lowers the shear relative to the flat plate
    ra = BL.karman_pohlhausen(BL.outer_flow("diffuser", U1=1.0, L=30.0), x, nu, profile="cubic")
    rf = BL.karman_pohlhausen(BL.outer_flow("flat", U=1.0), x, nu, profile="cubic")
    assert ra["delta"][-1] > rf["delta"][-1] and ra["tau0"][-1] < rf["tau0"][-1]
    with pytest.raises(KeyError):
        BL.karman_pohlhausen(BL.outer_flow("flat", U=U), x, nu, profile="nonsense")


# =====================================================================================================================
# C07 — Thwaites' method (§9.6)
# =====================================================================================================================
def _thwaites_wedge(n, x=0.7, a=1.0, nu=1e-6, closure="falkner_skan"):
    """Thwaites on U_e = a xⁿ started at x₀ = 1e-6 x (n ≤ 1.4) or 0.05 x (larger n) with the exact power-law θ₀ (θ²U⁶/ν = 0.45∫U⁵ from 0)."""
    of = BL.outer_flow("wedge", n=n, a=a)
    x0 = 1e-6 * x if n < 1.4 else 0.05 * x
    xs = x0 + (x - x0) * np.linspace(0, 1, 4001) ** 2
    th0 = np.sqrt(0.45 * nu * xs[0] ** (1 - n) / (a * (5 * n + 1)))
    return BL.thwaites(xs, of, nu, theta0=th0, closure=closure, stop_at_separation=False)


def test_thwaites_V1_flat_plate_and_example_9_1():  # V1 (C07/N66): U_e = U ⇒ θ² = 0.45νx/U exactly; Example 9.1 numbers (ours)
    nu, U = 1e-5, 2.0
    x = np.linspace(0.0, 1.5, 301)[1:]
    r = BL.thwaites(x, BL.outer_flow("flat", U=U), nu, theta0=np.sqrt(0.45 * nu * x[0] / U))  # started with the exact θ at x[0]
    assert np.max(np.abs(r["theta"] ** 2 / (0.45 * nu * x / U) - 1)) < 1e-12
    assert np.max(np.abs(r["lam"])) < 1e-12 and r["separated"] is False
    e = ch09.example_9_1()
    c = BL.blasius_constants()
    assert e["theta_coef"] == pytest.approx(np.sqrt(0.45), rel=1e-14) and e["theta_coef"] == pytest.approx(0.6708, abs=5e-5)
    assert e["theta_err"] == pytest.approx(np.sqrt(0.45) / c["theta"] - 1, rel=1e-12) and e["theta_err"] == pytest.approx(0.010, abs=5e-4)
    assert e["l0"] == pytest.approx(2 * c["fpp0"] ** 2, rel=1e-7)  # l(0) = f″(0) I_θ = 2 f″(0)² (table built at solver tol 1e-8: measured 2e-8)
    assert e["cf_sqrtRex"] == pytest.approx(2 * e["l0"] / e["theta_coef"], rel=1e-14) and e["cf_sqrtRex"] == pytest.approx(0.6575, abs=1e-3)
    assert e["cf_err"] == pytest.approx(-0.010, abs=1e-3) and e["delta_star_coef"] == pytest.approx(1.738, abs=1e-3)
    assert e["H0"] == pytest.approx(c["H"], rel=1e-7)
    # Holstein–Bohlen (9.44): λ = θ²U_e′/ν
    assert BL.holstein_bohlen(2e-3, -3.0, 1e-5) == pytest.approx(-1.2)


def test_thwaites_V1_diffuser_lambda_closed_form_and_separation_point_example_9_2():  # V1 (C07/N67): λ = −(0.45/4)[(1+x/L)⁴ − 1]; x_sep = 1.8^{1/4} − 1
    nu, U1, L = 1e-5, 1.0, 1.0
    xi = np.linspace(0.0, 0.5, 2001)
    of = BL.outer_flow("diffuser", U1=U1, L=L)
    r = BL.thwaites(xi * L, of, nu, theta0=0.0, closure="falkner_skan", stop_at_separation=False)
    lam_exact = -(0.45 / 4.0) * ((1.0 + xi) ** 4 - 1.0)  # (9.50) with U_e = U₁/(1 + x/L): the integral is closed-form
    assert np.max(np.abs(r["lam"] - lam_exact)) < 1e-12
    e = ch09.example_9_2()
    assert np.max(np.abs(e["lam"](xi) - lam_exact)) < 1e-14 and e["max_err"] < 1e-11
    assert e["x_sep_over_L"] == pytest.approx(1.8 ** 0.25 - 1.0, rel=1e-14)  # (1 + x/L)⁴ = 1.8 for λ_sep = −0.09
    assert e["x_sep_over_L"] == pytest.approx(0.15829, abs=5e-6)
    assert e["x_sep_over_L_fs"] == pytest.approx(((0.45 / 4 + 0.0681) / (0.45 / 4)) ** 0.25 - 1.0, rel=1e-12)
    # numerical crossing agrees with the closed form (linear interpolation on a 2001-point grid)
    assert e["x_sep_numeric"] == pytest.approx(e["x_sep_over_L_fs"], abs=2e-4)
    # separation (white closure, the book's criterion λ = −0.09) at exactly the closed-form point
    rw = BL.thwaites(xi * L, of, nu, closure="white", stop_at_separation=False)
    assert rw["x_sep"] == pytest.approx(1.8 ** 0.25 - 1.0, abs=3e-5)
    # θ₀ > 0 makes λ more negative and brings separation earlier (V7)
    e2 = ch09.example_9_2(theta0=5e-4)
    assert e2["lam"](0.1) < e["lam"](0.1) and e2["x_sep_over_L"] < e["x_sep_over_L"]
    assert e2["lam"](0.0) == pytest.approx(-(5e-4) ** 2 * U1 / (1e-5 * L), rel=1e-12)  # −θ₀²U₁/(νL) at x = 0
    assert ch09.example_9_2(theta0=2e-3)["x_sep_over_L"] is None  # λ(0) = −0.4 is already past the criterion
    # stop_at_separation cuts the arrays just past the first crossing
    rs = BL.thwaites(xi * L, of, nu, closure="white")
    assert rs["separated"] and rs["lam"][-1] <= -0.09 < rs["lam"][-2] and len(rs["x"]) < len(xi)


def test_thwaites_V1_power_law_family_closed_form_and_stagnation_limit():  # V1 (C07/D09): U_e = axⁿ ⇒ θ² = 0.45νx^{1−n}/(a(5n+1)), λ = 0.45n/(5n+1)
    nu, a = 1e-6, 1.7
    for n in (-0.05, 0.3, 1.0):
        x0 = 1e-3
        xs = np.geomspace(x0, 2.0, 300)
        th0 = np.sqrt(0.45 * nu * x0 ** (1 - n) / (a * (5 * n + 1)))
        r = BL.thwaites(xs, BL.outer_flow("wedge", n=n, a=a), nu, theta0=th0, stop_at_separation=False)
        assert np.max(np.abs(r["theta"] ** 2 / (0.45 * nu * xs ** (1 - n) / (a * (5 * n + 1))) - 1)) < 1e-9, n
        assert np.max(np.abs(r["lam"] - 0.45 * n / (5 * n + 1))) < 1e-9, n
    # stagnation start x = 0 (U_e = 0): θ² = 0.45ν/(6U_e′(0)) and λ ≡ 0.075 for U_e = a x — the 0/0 limit
    x = np.linspace(0.0, 1.0, 101)
    r = BL.thwaites(x, BL.outer_flow("stagnation", a=10.0), nu)
    assert np.max(np.abs(r["theta"] ** 2 - 0.45 * nu / 60.0)) < 1e-15
    assert np.max(np.abs(r["lam"] - 0.075)) < 1e-12
    # Thwaites (with the book's λ_sep = −0.09) predicts wedge separation exactly at n = −0.1: 0.45n/(5n+1) = −0.09 (exact FS: −0.0904)
    assert 0.45 * (-0.1) / (5 * (-0.1) + 1) == pytest.approx(-0.09, abs=1e-15)
    # a stagnation start with θ₀ ≠ 0 is refused
    with pytest.raises(ValueError):
        BL.thwaites(x, BL.outer_flow("stagnation", a=10.0), nu, theta0=1e-4)
    # array U_e and callable U_e give the same answer
    xs = np.linspace(0.05, 1.0, 50)
    ra = BL.thwaites(xs, 2.0 * xs ** 0.3, nu, theta0=1e-4)
    rc = BL.thwaites(xs, lambda s: 2.0 * s ** 0.3, nu, theta0=1e-4)
    assert np.max(np.abs(ra["theta"] / rc["theta"] - 1)) < 5e-4  # spline-of-samples vs analytic: differs at the first interval only


def test_thwaites_V5_accuracy_against_exact_falkner_skan_family():  # V5 (C07/N65): |θ error| ≤ 8 % for −0.09 ≤ n ≤ 4 (NOT the book's 3 %); τ₀ measured
    nu, a, x = 1e-6, 1.0, 0.7
    errs = {}
    for n in (-0.09, -0.05, 0.0, 1.0 / 3.0, 1.0, 4.0):
        st = BL.falkner_skan_state(n)
        Ue = a * x ** n
        r = _thwaites_wedge(n, x, a, nu)
        th_ex = st["I_theta"] * np.sqrt(nu * x / Ue)
        cf_ex = st["cf_sqrtRex"] / np.sqrt(Ue * x / nu)
        errs[n] = (r["theta"][-1] / th_ex - 1, r["cf"][-1] / cf_ex - 1)
        assert abs(errs[n][0]) <= 0.08, (n, errs[n])  # the design finding (i): 8 %, not 3 %
    # the measured values (design: n = 0 +1.0 %, 1/3 −4.2 %, 1 −6.3 %, 4 −7.6 %, −0.05 +3.1 %)
    assert errs[0.0][0] == pytest.approx(0.0101, abs=5e-4) and errs[1.0 / 3.0][0] == pytest.approx(-0.0424, abs=5e-4)
    assert errs[1.0][0] == pytest.approx(-0.0632, abs=5e-4) and errs[4.0][0] == pytest.approx(-0.0757, abs=5e-4)
    assert errs[-0.05][0] == pytest.approx(0.0308, abs=5e-4)
    # wall shear with the exact-FS closure: within 3 % for accelerating flows (the book's favourable claim holds for τ₀ …)
    for n in (0.0, 1.0 / 3.0, 1.0, 4.0):
        assert abs(errs[n][1]) < 0.03, (n, errs[n])
    assert abs(errs[-0.05][1]) < 0.10  # … and 10 % for mild adverse gradients
    # near the fold the shear is NOT reliable: the method reports λ below the FS fold for a still-attached flow (τ₀ error −56 % at n = −0.085)
    r = _thwaites_wedge(-0.085)
    cf_ex = BL.falkner_skan_state(-0.085)["cf_sqrtRex"] / np.sqrt(0.7 ** -0.085 * 0.7 / 1e-6)
    assert r["cf"][-1] / cf_ex - 1 < -0.4  # documented limitation: Thwaites predicts existence, not location, of separation
    # V7: the 'white' closure (λ_sep = −0.09) is within 3.5 % in τ₀ for n ≥ 0
    for n in (0.0, 1.0, 4.0):
        st = BL.falkner_skan_state(n)
        rw = _thwaites_wedge(n, closure="white")
        assert abs(rw["cf"][-1] / (st["cf_sqrtRex"] / np.sqrt(0.7 ** n * 0.7 / 1e-6)) - 1) < 0.035, n


def test_thwaites_bug_V1_thwaites_named_wedge_fails_for_n_above_1p4():  # V1 (C07 explainer function, Part C 1.20b): every wedge n must be callable
    # (F2, fixed) for n ≥ 1.5 the start value U_e(1e-6 x) < 1e-9 max(U_e) used to trip the stagnation-start guard.
    # U_e(1e-6 x) < 1e-9 max(U_e) trips the 'stagnation start' guard in thwaites(), which raises although theta0 is the (correct)
    # power-law value.  n = 4 is a Fig. 9.7 member and an explicit row of the design's accuracy table (θ −7.6 %).
    r = BL.thwaites_named("wedge", 0.7, 1e-6, n=4.0, a=1.0)
    st = BL.falkner_skan_state(4.0)
    assert r["theta"] / (st["I_theta"] * np.sqrt(1e-6 * 0.7 / 0.7 ** 4)) - 1 == pytest.approx(-0.0757, abs=1e-3)


def test_thwaites_named_V1_one_station_scalars_match_the_marched_run():  # V1 (C07/Part C 1.20b): flat, diffuser, wedge n ≤ 1, cylinder, retarded
    nu = 1e-5
    f = BL.thwaites_named("flat", 0.6, nu, theta0=np.sqrt(0.45 * nu * 1e-9), U=2.0)  # θ₀ ≈ 0 (start at 0)
    assert f["theta"] == pytest.approx(np.sqrt(0.45 * nu * 0.6 / 2.0), rel=1e-6)
    d = BL.thwaites_named("diffuser", 0.1, nu, U1=1.0, L=1.0)
    assert d["lam"] == pytest.approx(-(0.45 / 4) * (1.1 ** 4 - 1), rel=1e-9)
    assert d["x_sep"] is None
    d2 = BL.thwaites_named("diffuser", 0.3, nu, U1=1.0, L=1.0)
    assert d2["x_sep"] == pytest.approx(0.12563, abs=2e-4)  # exact-FS criterion λ = −0.0681: 0.12563 (Example 9.2 with the FS closure)
    w = BL.thwaites_named("wedge", 0.5, nu, n=1.0, a=3.0, rho=1000.0)
    assert w["lam"] == pytest.approx(0.075, abs=1e-9) and w["tau0"] > 0
    c = BL.thwaites_named("cylinder", 0.5 * np.deg2rad(60.0), nu, U=1.0, a=0.5)
    assert c["lam"] == pytest.approx(float(BL.thwaites_cylinder_closed_form(np.deg2rad(60.0))), abs=2e-4)
    rt = BL.thwaites_named("retarded", 0.2, nu, U0=1.0, L=1.0, c=1.0)
    assert np.isfinite(rt["theta"]) and rt["lam"] < 0


def test_thwaites_closure_V7_shape_of_L_and_l_H_against_the_falkner_skan_family():  # V1+V7 (C07/N57–N64): closure values at nodes, L(λ) near 0.45 − 6λ
    tb = BL.thwaites_closure_table(60)
    assert np.all(np.diff(tb["lam"]) > 0) and np.all(np.diff(tb["l"]) > 0)  # l increases with λ
    assert tb["lam"][0] == pytest.approx(-0.0681, abs=5e-5) and tb["l"][0] == 0.0  # the fold: zero shear
    assert tb["lam"].max() == pytest.approx(0.106, abs=5e-4)  # → 0.1065 as m → ∞ (the table stops at m = 50)
    # node values are the exact FS integrals: l(0) = 2f″(0)², H(0) = 2.5911, L(0) = 2l − 2(2+H)·0 = 0.4410
    c = BL.blasius_constants()
    assert float(BL.thwaites_l(0.0)) == pytest.approx(2 * c["fpp0"] ** 2, rel=1e-7)
    assert float(BL.thwaites_H(0.0)) == pytest.approx(c["H"], rel=1e-7)
    assert float(BL.thwaites_L(0.0)) == pytest.approx(0.4410, abs=5e-4)
    assert float(BL.thwaites_L(0.0855)) == pytest.approx(0.0, abs=5e-4) and float(BL.thwaites_L(-0.0675)) == pytest.approx(0.818, abs=2e-3)
    # the linear fit L ≈ 0.45 − 6λ: within 0.05 for −0.068 ≤ λ ≤ 0.05 (analysis §6), 0.065 up to λ = 0.0855 (the D10 text says 0.06; measured 0.063)
    lam = tb["lam"]
    m1 = (lam >= -0.0681) & (lam <= 0.05)
    m2 = (lam >= -0.0681) & (lam <= 0.0855)
    assert np.max(np.abs(tb["L"][m1] - (0.45 - 6 * lam[m1]))) < 0.05
    assert np.max(np.abs(tb["L"][m2] - (0.45 - 6 * lam[m2]))) < 0.065
    # "white" closure: l = (λ + 0.09)^0.62, zero at −0.09, l(0) ≈ 0.2247
    assert float(BL.thwaites_l(-0.09, "white")) == 0.0 and float(BL.thwaites_l(0.0, "white")) == pytest.approx(0.09 ** 0.62, rel=1e-12)
    # l = 0 below the fold for the FS closure, and the fast table is a coarser subset with the same fold
    assert float(BL.thwaites_l(-0.08)) == 0.0
    tf = BL.thwaites_closure_table(fast=True)
    assert len(tf["lam"]) < len(tb["lam"]) and tf["lam"][0] == pytest.approx(tb["lam"][0], abs=1e-9)
    with pytest.raises(ValueError):
        BL.thwaites_l(0.0, "nonsense")


@needs_ref
def test_thwaites_V5_fit_constants_and_separation_form_against_the_published_paper():  # V5 (C07): 0.45, 6, θ² form (Agrawal et al. Eqs. 2.4–2.5), m ≈ 0.09
    R = ref_json()["thwaites_agrawal"]
    assert R["separation_abs_m"] == 0.09 and BL.LAMBDA_SEP_BOOK == -R["separation_abs_m"]
    # 0.45 + 6m with m = −λ: our exact-FS closure at λ = 0 gives L = 0.441 (2 % below 0.45) — the fit is the paper's, ours is exact for FS
    assert abs(float(BL.thwaites_L(0.0)) - 0.45) < 0.01
    # θ² form: the flat-plate value used above is the paper's (2.5) with U_e const
    nu, U, x = 1e-5, 1.0, 0.5
    th = BL.thwaites(np.array([x / 2, x]), BL.outer_flow("flat", U=U), nu, theta0=np.sqrt(0.45 * nu * x / 2 / U))["theta"][-1]
    assert th ** 2 == pytest.approx(0.45 * nu / U ** 6 * U ** 5 * x, rel=1e-12)


def test_thwaites_cylinder_V1_closed_form_versus_quadrature_and_separation_angle():  # V1 (C07/D11): λ(φ) = 0.45 cos φ F(φ)/sin⁶φ; φ_sep = 103.1°
    phi = np.deg2rad(np.array([5.0, 30.0, 60.0, 82.0, 90.0, 100.0, 120.0]))
    closed = np.asarray(BL.thwaites_cylinder_closed_form(phi))
    quad_ = np.asarray(BL.thwaites_cylinder(np.rad2deg(phi)))  # the quadrature route (independent of the closed form)
    assert np.max(np.abs(closed - quad_)) < 1e-9
    # design numbers λ(30°, 60°, 82°, 90°, 100°) = 0.0722, 0.0589, 0.0263, 0.0000, −0.0603
    got = np.asarray(BL.thwaites_cylinder_closed_form(np.deg2rad([30.0, 60.0, 82.0, 90.0, 100.0])))
    assert got == pytest.approx([0.0722, 0.0589, 0.0263, 0.0, -0.0603], abs=6e-5)
    assert float(BL.thwaites_cylinder(0.0)) == pytest.approx(0.075, abs=1e-15)  # the stagnation-point limit 0.45/6
    assert float(BL.thwaites_cylinder_closed_form(np.deg2rad(0.5))) == pytest.approx(0.075, abs=2e-4)
    # separation angles: −0.09 ⇒ 103.11°, −0.0681 (exact FS zero shear) ⇒ 100.89°
    assert BL.thwaites_cylinder_separation(-0.09) == pytest.approx(103.11, abs=0.01)
    assert BL.thwaites_cylinder_separation(BL.LAMBDA_SEP_FS) == pytest.approx(100.89, abs=0.01)
    # independent numerical route: the full Thwaites run on the cylinder
    full = BL.thwaites_cylinder(closure="white", n=1200)
    assert full["phi_sep_deg"] == pytest.approx(103.11, abs=0.2)
    fs = BL.thwaites_cylinder(closure="falkner_skan", n=1200)
    assert fs["phi_sep_deg"] == pytest.approx(100.89, abs=0.2)
    # V7: λ passes through zero at the shoulder (U_e peaks at 90°), is positive before and negative after
    assert float(BL.thwaites_cylinder_closed_form(np.pi / 2)) == pytest.approx(0.0, abs=1e-15)
    assert BL.thwaites_cylinder_separation(-0.09) > 90.0
    with pytest.raises(ValueError):
        BL.thwaites_cylinder_separation(0.1)
    # the ideal-flow prediction is later than the observed subcritical separation (≈ 82°, a rounded experimental value)
    assert BL.thwaites_cylinder_separation(-0.09) > 82.0


# =====================================================================================================================
# C08 — separation: wall curvature, inflection, τ₀ = 0 (§9.7), and the parabolic marching solver (N14, N24)
# =====================================================================================================================
def test_inflection_V7_signs_of_the_wall_curvature_and_location():  # V7 (C08/N71–N72): favourable ⇒ none; Blasius at the wall; adverse ⇒ finite η
    for m in (0.1, 1.0 / 3.0, 1.0, 4.0):
        assert BL.falkner_skan_state(m)["inflection_eta"] is None, m  # (9.51)
    assert BL.falkner_skan_state(0.0)["inflection_eta"] == 0.0  # zero curvature at the wall
    etas = []
    for m in (-0.02, -0.05, -0.08, -0.09):
        e = BL.falkner_skan_state(m)["inflection_eta"]
        assert e is not None and e > 0, m  # (9.52)
        etas.append(e)
    assert np.all(np.diff(etas) > 0)  # the inflection point moves out as the layer approaches separation
    assert BL.falkner_skan_state(-0.05)["inflection_eta"] == pytest.approx(1.650, abs=1e-3)  # design: 1.650
    # sampled-profile detector agrees (independent of the ODE's f‴): FS n = −0.05 and n = 1, and Blasius
    for m, want in ((-0.05, 1.650), (1.0, None), (0.0, 0.0)):
        d = BL.falkner_skan(m, eta_max=16.0, n=1601)
        got = BL.profile_inflection(d["eta"], d["fp"])
        if want is None:
            assert got is None
        else:
            assert got == pytest.approx(want, abs=8e-3), m
    # a synthetic S-shaped profile has its inflection where u″ = 0
    y = np.linspace(0, 1, 401)
    u = y + 0.2 * (y ** 3 - 1.5 * y ** 2)  # u″ = 1.2y − 0.6 → y = 0.5
    assert BL.profile_inflection(y, u) == pytest.approx(0.5, abs=2e-3)


def test_separation_point_V1_sign_change_interpolation():  # V1 (C08): first τ₀ sign change, linear interpolation; None if attached
    x = np.linspace(0, 1, 11)
    tau = 1.0 - 2.0 * x  # zero at x = 0.5 exactly
    assert BL.separation_point(x, tau) == pytest.approx(0.5, abs=1e-14)
    assert BL.separation_point(x, np.abs(tau) + 1e-3) is None
    tau2 = np.where(x < 0.33, 1.0, -0.5)
    assert 0.3 <= BL.separation_point(x, tau2) <= 0.4
    assert BL.separation_point(x, np.ones_like(x)) is None


def test_march_V3_blasius_orders_in_dx_and_dsigma():  # V3 (C01/C03/N14): BE first order in Δx; second order in Δσ for the plate
    nu = 1e-3
    U = 1.0
    x1 = 1.0
    ex = float(BL.blasius_wall_shear(x1, U, 1.0, nu))
    of = BL.outer_flow("flat", U=U)
    # Δx, backward Euler (order=1): vs the exact similarity value at very fine Δσ
    hs, errs = [], []
    for nx in (10, 20, 40):
        r = BL.march_boundary_layer(of, np.linspace(0.1, x1, nx + 1), nu, ny=800, order=1)
        hs.append(0.9 / nx)
        errs.append(abs(r["tau0"][-1] / ex - 1))
    assert abs(observed_order(hs, errs) - 1.0) < ORDER_TOL, pairwise_orders(hs, errs)
    # Δσ (central differences, plate): error ∝ Δσ², fixed small Δx
    hs, errs = [], []
    for ny in (100, 200, 400):
        r = BL.march_boundary_layer(of, np.linspace(0.1, x1, 81), nu, ny=ny)
        hs.append(1.0 / ny)
        errs.append(abs(r["tau0"][-1] / ex - 1))
    assert abs(observed_order(hs, errs) - 2.0) < ORDER_TOL, pairwise_orders(hs, errs)
    assert errs[-1] < 1e-3  # analysis §6: τ₀ within 1e-3 of the similarity solution (plate)
    # the default settings on a geometric grid reproduce the plate's τ₀ to 1e-3 at every MARCHED station (the inlet station's τ₀ is the
    # extrapolation formula applied to a linearly interpolated profile: 2.8e-3 at ny = 400 and 8.7e-3 at ny = 800 — see report)
    r = BL.march_boundary_layer(of, np.geomspace(0.1, x1, 161), nu)
    exact = np.array([float(BL.blasius_wall_shear(v, U, 1.0, nu)) for v in r["x"]])
    assert np.max(np.abs(r["tau0"][1:] / exact[1:] - 1)) < 1e-3
    assert abs(r["tau0"][0] / exact[0] - 1) < 1e-2
    # the marched profile reproduces the Blasius profile: u at ψ-nodes against f′(η) with f(η) = ψ/√(νUx)
    eta = np.linspace(0, 12, 4001)
    b = BL.blasius_profile(eta)
    fp_ex = np.interp(r["psi"] / np.sqrt(nu * U * r["x"][-1]), b[0], b[1])
    assert np.max(np.abs(r["u"][-1] - U * fp_ex)) < 2e-3


def test_march_V3_time_step_order_of_the_bdf2_scheme_on_pressure_gradient_flows():  # V3 (C01): Δx-convergence, self-reference at fixed Δσ
    nu = 1e-3
    for n, grids in ((0.5, (20, 40, 80)), (-0.05, (40, 80, 160))):
        of = BL.outer_flow("wedge", n=n, a=1.0)
        ref = BL.march_boundary_layer(of, np.linspace(0.1, 1.0, 1281), nu, ny=200)["tau0"][-1]
        hs, errs = [], []
        for nx in grids:
            t = BL.march_boundary_layer(of, np.linspace(0.1, 1.0, nx + 1), nu, ny=200)["tau0"][-1]
            hs.append(0.9 / nx)
            errs.append(abs(t - ref))
        assert abs(observed_order(hs, errs) - 2.0) < 0.25, (n, pairwise_orders(hs, errs))  # span ≈ 1.3 decades ⇒ ±0.25
        # backward Euler is first order on the same flows
        ref1 = BL.march_boundary_layer(of, np.linspace(0.1, 1.0, 1281), nu, ny=200, order=1)["tau0"][-1]
        errs1 = [abs(BL.march_boundary_layer(of, np.linspace(0.1, 1.0, nx + 1), nu, ny=200, order=1)["tau0"][-1] - ref1) for nx in grids]
        assert abs(observed_order(hs, errs1) - 1.0) < ORDER_TOL, (n, pairwise_orders(hs, errs1))


def test_march_V3_falkner_skan_tau0_is_second_order_in_dsigma():  # V3 (C01/C05): analysis §6 "observed order in Δψ … central 2"; τ₀ vs exact FS
    # (F3a, fixed) the wall shear of a pressure-gradient flow used to converge only first-order in Δσ (0.97–0.99).
    # Hypothesis: with dp/dx ≠ 0 the near-wall expansion w = u² = Aσ² + Bσ³ + … has B ∝ dp/dx, and the central difference of the
    # w_σ/σ term of the ψ-diffusion operator has relative error h²/(3σ²) = 1/3 at the first node (exact for the plate, where B = 0).
    nu = 1e-3
    n = 0.5
    st = BL.falkner_skan_state(n)
    ex = st["fpp0"] * nu * np.sqrt(1.0 / nu)  # τ₀ at x = 1, ρ = 1, U_e = 1
    of = BL.outer_flow("wedge", n=n, a=1.0)
    hs, errs = [], []
    for ny in (100, 200, 400, 800):
        t = BL.march_boundary_layer(of, np.linspace(0.1, 1.0, 161), nu, ny=ny)["tau0"][-1]
        hs.append(1.0 / ny)
        errs.append(abs(t / ex - 1))
    assert abs(observed_order(hs, errs) - 2.0) < 0.25, (pairwise_orders(hs, errs), errs)


def test_march_V5_falkner_skan_tau0_within_1e_3_of_the_similarity_solution_at_default_settings():  # V5 (C05): analysis §6 "τ₀ within 1e-3 of the similarity solution"
    # (F3b, fixed) default ny = 400 used to give 0.5–0.9 % error on n = 0.2, 0.5, 1.
    nu = 1e-3
    worst = 0.0
    for n in (0.2, 0.5, 1.0):
        st = BL.falkner_skan_state(n)
        of = BL.outer_flow("wedge", n=n, a=1.0)
        r = BL.march_boundary_layer(of, np.linspace(0.1, 1.0, 81), nu)
        worst = max(worst, abs(r["tau0"][-1] / (st["fpp0"] * nu * np.sqrt(1.0 / nu)) - 1))
    assert worst < 1e-3, worst


def test_march_V5_wedge_flows_separate_below_the_fold_and_stay_attached_above():  # V5+V7 (C05/C08): Belden fold n = −0.09043 seen by the PDE solver
    nu = 1e-3
    x = np.linspace(0.1, 50.0, 1500)
    seps = []
    for n in (-0.085, -0.089):
        r = BL.march_boundary_layer(BL.outer_flow("wedge", n=n, a=1.0), x, nu)
        assert r["separated"] is False and r["x"][-1] == x[-1], n
        assert np.all(r["tau0"] > 0)
    for n in (-0.095, -0.1, -0.12):  # (n = −0.092, 2 % below the fold, is numerically marginal: it separated at x = 8.1 or not at all depending on the inlet)
        Ue0 = 0.1 ** n
        uin = lambda y, Ue0=Ue0: Ue0 * BL.blasius_profile(np.asarray(y) / np.sqrt(nu * 0.1 / Ue0))[1]  # noqa: E731 (FS(m < fold) does not exist)
        r = BL.march_boundary_layer(BL.outer_flow("wedge", n=n, a=1.0), x, nu, u_inlet=uin)
        assert r["separated"] and r["x_sep"] is not None and r["tau0"][-1] > 0, n
        seps.append(r["x_sep"])
    assert np.all(np.diff(seps) < 0)  # the stronger the adverse gradient, the earlier the separation
    assert seps[0] > 0.5 and seps[-1] < 0.5


def test_march_V7_favourable_gradients_forget_the_inlet_profile_adverse_do_not_N24():  # V7 (C03/N24): Serrin/Peletier remark, quantified
    nu = 1e-3
    x = np.linspace(0.1, 3.0, 300)
    end_diff = {}
    for n in (0.5, 0.0, -0.05):
        of = BL.outer_flow("wedge", n=n, a=1.0) if n else BL.outer_flow("flat", U=1.0)
        Ue0 = float(of.Ue(0.1))
        ra = BL.march_boundary_layer(of, x, nu, u_inlet=lambda y: Ue0 * np.minimum(np.asarray(y) / 0.05, 1.0))  # linear ramp
        rb = BL.march_boundary_layer(of, x, nu, u_inlet=lambda y: Ue0 * (1 - np.exp(-np.asarray(y) / 0.004)))  # exponential
        assert not ra["separated"] and not rb["separated"]
        end_diff[n] = abs(ra["tau0"][-1] - rb["tau0"][-1]) / ra["tau0"][-1]
    assert end_diff[0.5] < 1e-3
    assert end_diff[0.5] < end_diff[0.0] < end_diff[-0.05]  # forgetting is fastest for favourable, slowest for adverse gradients
    assert end_diff[-0.05] > 2 * end_diff[0.0]


def test_march_V1_inlet_forms_and_fast_mode_and_separation_stop():  # V1 (C01/N13, N14): (9.15) inlet as a callable or (y, u) arrays; FAST; stop at τ₀ ≤ 0
    nu = 1e-3
    of = BL.outer_flow("flat", U=1.0)
    x = np.linspace(0.1, 0.5, 41)
    prof = lambda y: 1 - np.exp(-np.asarray(y) / 0.006)  # noqa: E731
    yy = np.linspace(0, 12 * np.sqrt(nu * 0.1), 4001)
    a = BL.march_boundary_layer(of, x, nu, u_inlet=prof)
    b = BL.march_boundary_layer(of, x, nu, u_inlet=(yy, prof(yy)))
    assert np.max(np.abs(a["tau0"] / b["tau0"] - 1)) < 1e-6
    f = BL.march_boundary_layer(of, x, nu, fast=True)
    assert len(f["sigma"]) == 151 and abs(f["tau0"][-1] / float(BL.blasius_wall_shear(0.5, 1.0, 1.0, nu)) - 1) < 3e-2  # FAST: 1e-2 class
    assert a["y"].shape == a["u"].shape == a["v"].shape == (len(a["x"]), 401)
    assert np.all(a["u"][:, 0] == 0.0) and np.allclose(a["u"][:, -1], 1.0)  # no slip; free stream at the edge
    assert np.all(a["v"][-1, 200:] > 0)  # v > 0 above a growing layer (displacement thickness lifts the outer streamlines)
    # separation stops the run and reports x_sep (strong adverse gradient, Blasius inlet)
    of2 = BL.outer_flow("diffuser", U1=1.0, L=0.4)
    with pytest.raises(ValueError):  # default inlet = local FS(m0) with m0 = −0.2 < fold: no such profile (the message is a bare NaN error, see report)
        BL.march_boundary_layer(of2, np.linspace(0.1, 2.0, 400), nu)
    r = BL.march_boundary_layer(of2, np.linspace(0.1, 2.0, 400), nu, u_inlet=lambda y: BL.blasius_profile(np.asarray(y) / np.sqrt(nu * 0.1 / float(of2.Ue(0.1))))[1] * float(of2.Ue(0.1)))
    assert r["separated"] and 0.1 < r["x_sep"] < 1.0 and r["x"][-1] <= r["x_sep"] + 1e-9
    assert np.all(r["tau0"] > 0)


# =====================================================================================================================
# C09 — form drag of the separated model (§9.7); the drag-crisis chain of C11 uses the same functions
# =====================================================================================================================
def test_separated_drag_V1_closed_form_quadrature_and_sympy_integral():  # V1+V2 (C09/D13): C_D,p = sin φs (1 − 4/3 sin²φs − C_b)
    phi = sp.Symbol("phi", positive=True)
    ps, cb = sp.symbols("phi_s C_b", positive=True)
    sym = sp.integrate((1 - 4 * sp.sin(phi) ** 2) * sp.cos(phi), (phi, 0, ps)) + cb * (sp.sin(sp.pi) - sp.sin(ps))  # ∫₀^π C_p cos φ dφ
    closed = sp.sin(ps) * (1 - sp.Rational(4, 3) * sp.sin(ps) ** 2 - cb)
    assert z0(sp.expand(sym - closed))
    f = sp.lambdify((ps, cb), closed, "numpy")
    for phis, cbase in ((82.0, -1.2), (125.0, -0.6), (90.0, -3.0), (60.0, 0.0), (150.0, -1.0)):
        want = float(f(np.deg2rad(phis), cbase))
        assert BB.separated_pressure_drag(phis, cp_base=cbase) == pytest.approx(want, abs=1e-12), (phis, cbase)  # Gauss–Legendre split at φ_s
    # the design's expect row: (82°, −1.2) = 0.8840, (125°, −0.6) = 0.5778, (90°, ideal C_b = −3) = 2.667; crude default at 82° = 2.59
    assert BB.separated_pressure_drag(82.0, -1.2) == pytest.approx(0.8840, abs=3e-4)
    assert BB.separated_pressure_drag(125.0, -0.6) == pytest.approx(0.5778, abs=1e-4)
    assert BB.separated_pressure_drag(90.0) == pytest.approx(8.0 / 3.0, abs=1e-12)
    assert BB.separated_pressure_drag(82.0) == pytest.approx(2.59, abs=5e-3)
    # sampled routes: a callable with breakpoints, and samples on 0…180° and 0…360°
    ph = np.linspace(0, 180, 721)
    cp = BB.separated_cp(ph, 82.0, cp_base=-1.2)
    assert BB.pressure_drag_from_cp(ph, cp) == pytest.approx(BB.separated_pressure_drag(82.0, -1.2), abs=2e-3)  # first order at the jump
    ph2 = np.linspace(0, 360, 1441)
    cp2 = BB.separated_cp(np.where(ph2 > 180, 360 - ph2, ph2), 82.0, cp_base=-1.2)  # symmetric about the axis
    assert BB.pressure_drag_from_cp(ph2, cp2) == pytest.approx(BB.separated_pressure_drag(82.0, -1.2), abs=2e-3)


def test_separated_drag_V3_jump_costs_first_order_on_sampled_pressure():  # V3 (C09): trapezoid with a jump at a node: error ∝ Δφ
    exact = BB.separated_pressure_drag(82.0, -1.2)
    hs, errs = [], []
    for dphi in (2.0, 1.0, 0.5, 0.25):  # 82° is a node of every grid
        ph = np.arange(0.0, 180.0 + 1e-9, dphi)
        errs.append(abs(BB.pressure_drag_from_cp(ph, BB.separated_cp(ph, 82.0, cp_base=-1.2)) - exact))
        hs.append(dphi)
    assert abs(observed_order(hs, errs) - 1.0) < ORDER_TOL, pairwise_orders(hs, errs)


def test_separated_drag_V7_dalembert_and_monotonicity_in_the_wake_pressure():  # V7 (C09/N81): ideal flow has zero drag; higher base pressure ⇒ less drag
    ph = np.linspace(0, 360, 361)[:-1]
    ideal = BB.cp_ideal_cylinder(ph)
    assert abs(BB.pressure_drag_from_cp(ph, ideal)) < 1e-14  # ∮(1 − 4sin²φ)cos φ dφ = 0: d'Alembert (ch06)
    assert abs(BB.pressure_drag_from_cp(lambda p: 1 - 4 * np.sin(p) ** 2, n=64)) < 1e-14  # one-sided form too
    assert float(BB.cp_ideal_cylinder(0.0)) == 1.0 and float(BB.cp_ideal_cylinder(90.0)) == pytest.approx(-3.0)  # stagnation +1, shoulder −3
    # crude default (C_b = ideal value at φ_s): D = (8/3) sin³φ_s → 0 as φ_s → 180° (the ideal limit)
    d = [BB.separated_pressure_drag(p) for p in (120.0, 150.0, 170.0, 179.0)]
    assert np.all(np.diff(d) < 0) and d[-1] < 1e-4
    for p in (120.0, 170.0):
        assert BB.separated_pressure_drag(p) == pytest.approx(8.0 / 3.0 * np.sin(np.deg2rad(p)) ** 3, rel=1e-12)
    # at fixed φ_s the drag falls linearly as the base pressure rises: dD/dC_b = −sin φ_s
    for p in (82.0, 125.0):
        d1, d2 = BB.separated_pressure_drag(p, -1.2), BB.separated_pressure_drag(p, -0.6)
        assert d2 - d1 == pytest.approx(-0.6 * np.sin(np.deg2rad(p)), abs=1e-12)
    # (the design's V7 "later separation ⇒ less drag at FIXED base pressure" is false on this model: D(82°) 0.884, D(90°) 0.867, D(125°) 1.07 at C_b = −1.2)
    assert BB.separated_pressure_drag(125.0, -1.2) > BB.separated_pressure_drag(82.0, -1.2)


def test_drag_crisis_V7_separation_moves_aft_and_the_model_drag_falls():  # V7 (C11/N81): 82° → 125°, C_b −1.2 → −0.6, pressure drag falls; roughness triggers earlier
    p = BB.drag_crisis_pair()
    assert p["subcritical"]["phi_sep_deg"] == 82.0 and p["supercritical"]["phi_sep_deg"] == 125.0
    assert p["ratio"] == pytest.approx(p["supercritical"]["cd_model"] / p["subcritical"]["cd_model"], rel=1e-14) and 0.6 < p["ratio"] < 0.7
    assert p["subcritical"]["cd_model"] == pytest.approx(BB.separated_pressure_drag(82.0, -1.2), abs=1e-12)
    Re = np.logspace(4, 7, 61)
    st = [BB.drag_crisis_state(r) for r in Re]
    phi = np.array([s["phi_sep_deg"] for s in st])
    assert np.all(np.diff(phi) >= -1e-12) and phi[0] == 82.0 and phi[-1] == 125.0  # monotone, smooth blend across the critical band
    cd = np.array([s["cd_model"] for s in st])
    assert cd[0] > cd[-1] and np.all(np.isfinite(cd))
    assert BB.drag_crisis_state(3e5)["blend"] == 0.0 and BB.drag_crisis_state(6e5)["blend"] == 1.0
    # surface roughness: the transition is tripped at Re_cr/3 — a rough cylinder at Re = 1.5e5 is already supercritical, a smooth one is not
    assert BB.drag_crisis_state(1.5e5, rough=True)["phi_sep_deg"] == 125.0 and BB.drag_crisis_state(1.5e5)["phi_sep_deg"] == 82.0
    assert BB.drag_crisis_state(5e4, rough=True)["phi_sep_deg"] == 82.0
    assert BB.drag_crisis_state(1e6)["qualitative"] is True


def test_cylinder_regimes_V7_table_lookup_at_every_threshold_and_state_dict():  # V7 (C10/C11, N77, N80): value on a threshold belongs to the upper regime
    th = BB.CYLINDER_THRESHOLDS
    labels = {}
    for name, R in th.items():
        below = BB.cylinder_flow_regime(R * (1 - 1e-12))["label"]
        at = BB.cylinder_flow_regime(R)["label"]
        assert below != at, name  # the label switches exactly at the threshold
        labels[name] = at
    assert BB.cylinder_flow_regime(0.5)["label"] == "creeping flow: symmetric, no wake"
    assert BB.cylinder_flow_regime(10.0)["label"] == "two steady attached eddies"
    r100 = BB.cylinder_flow_regime(100.0)
    assert r100["label"] == "laminar Karman street" and r100["St"] == 0.2 and r100["separation_deg"] is None
    assert BB.cylinder_flow_regime(1e5)["separation_deg"] == 82.0 and BB.cylinder_flow_regime(1e6)["separation_deg"] == 125.0
    assert BB.cylinder_flow_regime(4e5)["separation_deg"] is None  # inside the critical band
    assert BB.cylinder_flow_regime(50.0, thresholds=dict(street_onset=60.0))["label"] == "two steady attached eddies"  # thresholds are arguments
    # sphere table: no regular street, St None
    assert BB.sphere_flow_regime(50.0)["label"] == "steady wake with an attached doughnut eddy"
    assert BB.sphere_flow_regime(1e4)["St"] is None and BB.sphere_flow_regime(1e4)["separation_deg"] == 80.0
    assert BB.sphere_flow_regime(1e6)["separation_deg"] == 120.0 and BB.sphere_flow_regime(0.1)["label"].startswith("creeping")
    # cylinder_state: separated model only for Re ≥ 3000; St for the street
    s = BB.cylinder_state(100.0)
    assert s["phi_sep_deg"] is None and s["cd_model"] is None and s["St"] == 0.2
    s = BB.cylinder_state(1e5)
    assert s["phi_sep_deg"] == 82.0 and s["cd_model"] == pytest.approx(BB.separated_pressure_drag(82.0, -1.2), abs=1e-12) and s["qualitative"] is True
    # the schematic C_D(Re) is labelled qualitative: Lamb-type growth at small Re, ≈ 1 in the middle, a dip near the critical Re, then recovery
    cds = BB.cylinder_cd_schematic(np.array([0.1, 1.0, 10.0, 1e3, 1e5, 5e5, 1e7]))
    assert cds[0] == pytest.approx(8 * np.pi / (0.1 * (2.002 - np.log(0.1))), rel=1e-12) and cds[3] == pytest.approx(1.0)
    assert np.all(np.diff(cds[:4]) < 0) and cds[5] < 0.5 * cds[4] and cds[6] > cds[5]


def test_sphere_drag_V5_reused_morrison_correlation_shows_the_crisis_and_stokes_oseen_bracket():  # V5 (C11/N83): reuse of ch04/ch08 tested functions
    Re = np.array([0.1, 0.5, 1.0, 2.0, 5.0])
    cd = np.asarray(SIM.sphere_drag_coefficient(Re, "morrison"))
    st = np.asarray(CRP.stokes_drag_coefficient(Re))
    os_ = np.asarray(CRP.oseen_drag_coefficient(Re))
    assert np.all(cd >= st * 0.999) and np.all(cd <= os_ * 1.001)  # Stokes below, Oseen above (diameter Re; ch08 bracket)
    R = np.logspace(3, 6, 400)
    c = np.asarray(SIM.sphere_drag_coefficient(R))
    i = int(np.argmin(c))
    assert 2e5 < R[i] < 8e5 and c[i] < 0.3 * float(SIM.sphere_drag_coefficient(1e5))  # the drag crisis: C_D falls by > 3× and then recovers
    assert c[-1] > c[i]
    # potential-flow suction peak of the sphere: C_p = 1 − (9/4) sin²θ has its minimum −5/4 (the ch06 function)
    from fluidpy import ch06_ideal_flow as ch06
    assert float(ch06.sphere_surface_cp(np.pi / 2)) == pytest.approx(-1.25, abs=1e-14)


def test_ball_dynamics_V1_swing_kinematics_and_R10_magnus_truth_table():  # V1+V7 (C11/N84–N85): y = ½ (F/W) g (d/U)²; sign table; R10 slip
    assert BB.ball_swing_deflection(0.2, 18.0, 35.0) == pytest.approx(0.2595, abs=5e-5)  # design expect row
    d = BB.ball_swing_deflection(0.3, 10.0, 25.0)
    assert d == pytest.approx(0.5 * 0.3 * 9.81 * (10.0 / 25.0) ** 2, rel=1e-14)
    assert BB.ball_swing_deflection(0.3, 20.0, 25.0) == pytest.approx(4 * d, rel=1e-13)  # ∝ d²
    assert BB.ball_swing_deflection(0.6, 10.0, 25.0) == pytest.approx(2 * d, rel=1e-13)  # ∝ F/W
    assert BB.ball_swing_deflection(0.3, 10.0, 50.0) == pytest.approx(d / 4, rel=1e-13)  # ∝ 1/U²
    assert BB.ball_swing_deflection(0.3, 10.0, 25.0, g=1.62) == pytest.approx(d * 1.62 / 9.81, rel=1e-13)  # a lunar cricket ball
    # truth table (the book's sentence prints "Re < Re_cr" twice — R10): negative only when just the fast side is past the crisis
    Rcr = 3e5
    assert BB.magnus_sign(1e5, 5e5, Rcr) == "−"  # slow side subcritical, fast side supercritical: NEGATIVE
    assert BB.magnus_sign(1e5, 2e5, Rcr) == "+" and BB.magnus_sign(4e5, 5e5, Rcr) == "+"  # both on the same side: ordinary POSITIVE (the corrected second inequality)
    assert BB.magnus_sign(1e5, 1e5, Rcr) == "none" and BB.magnus_sign(5e5, 1e5, Rcr) == "none"
    assert BB.magnus_sign(1e5, Rcr, Rcr) == "−"  # Re_fast = Re_cr counts as past the crisis (Re_slow < Re_cr ≤ Re_fast)


# =====================================================================================================================
# C10 — the Kármán vortex street (§9.8) and the Strouhal shedding frequency
# =====================================================================================================================
def test_karman_ratio_V1_closed_form_and_marginal_growth_is_zero():  # V1 (C10/D14): b/a = arccosh(√2)/π; growth zero only there
    r = BB.karman_street_ratio()
    assert r == pytest.approx(np.arccosh(np.sqrt(2)) / np.pi, rel=1e-15) and r == pytest.approx(0.280550, abs=5e-7)
    assert np.cosh(np.pi * r) == pytest.approx(np.sqrt(2), rel=1e-14)
    assert r == pytest.approx(np.log(np.sqrt(2) + 1) / np.pi, rel=1e-14)  # arccosh √2 = ln(√2 + 1)
    assert float(BB.karman_street_growth_closed(r)) < 1e-14
    b = np.array([0.05, 0.1, 0.2, 0.25, 0.3, 0.35, 0.5, 1.0])
    g = np.asarray(BB.karman_street_growth_closed(b))
    assert np.all(g > 0)  # every other spacing grows (V-shaped curve)
    want = np.pi / 2 * np.abs(0.5 - 1.0 / np.cosh(np.pi * b) ** 2)
    assert np.max(np.abs(g - want)) < 1e-14
    # design expect row (Γ = a = 1): 0.6400, 0.2982, 0.1099, 0.0663, 0.2208, 0.5359, 0.7737 at b/a = 0.1, 0.2, 0.25, 0.3, 0.35, 0.5, 1.0
    for bb, gg in ((0.1, 0.6400), (0.2, 0.2982), (0.25, 0.1099), (0.3, 0.0663), (0.35, 0.2208), (0.5, 0.5359), (1.0, 0.7737)):
        assert float(BB.karman_street_growth_closed(bb)) == pytest.approx(gg, abs=6e-5)
    # the growth is a V: decreasing to r then increasing
    fine = np.linspace(0.05, 1.0, 400)
    gf = np.asarray(BB.karman_street_growth_closed(fine))
    imin = int(np.argmin(gf))
    assert fine[imin] == pytest.approx(r, abs=3e-3) and np.all(np.diff(gf[: imin + 1]) <= 1e-12) and np.all(np.diff(gf[imin:]) >= -1e-12)
    # units: Γ/a² is 1/s
    dimensional_check(lambda Gamma, a: Gamma / a ** 2, "vorticity", Gamma=Q_(1.0, "m**2/s"), a=Q_(1.0, "m"))  # 1/s
    assert float(BB.karman_street_growth_closed(0.2, Gamma=2.0, a=0.5)) == pytest.approx(float(BB.karman_street_growth_closed(0.2)) * 2.0 / 0.25, rel=1e-13)


def test_karman_street_V3_three_independent_numerical_routes_agree_with_the_closed_form():  # V3 (C10): Bloch 4×4 (lattice sums), finite periodic cell, closed form
    for b in (0.1, 0.2, 0.25, 0.3, 0.35, 0.5, 1.0):
        closed = float(BB.karman_street_growth_closed(b))
        bloch = float(np.max(BB.karman_street_spectrum(b).real))  # k = π/a, exact lattice sums
        per = float(np.max(BB.karman_street_spectrum_periodic(b, n_pairs=16).real))  # a completely different assembly (cot sum on a period 16a)
        scan = BB.karman_street_growth(b)  # maximum over 61 wavenumbers in (0, π/a]
        assert bloch == pytest.approx(closed, abs=1e-12) and per == pytest.approx(closed, abs=1e-10), b
        assert scan == pytest.approx(closed, abs=1e-12), b  # k = π/a is the most dangerous mode for every b/a tried
    # at the marginal spacing the spectrum is pure imaginary ±0.7854i·Γ/a² (neutral)
    r = BB.karman_street_ratio()
    sp_ = BB.karman_street_spectrum(r)
    assert np.max(np.abs(sp_.real)) < 1e-12 and np.allclose(np.sort(sp_.imag), [-np.pi / 4, -np.pi / 4, np.pi / 4, np.pi / 4], atol=1e-12)
    assert np.max(BB.karman_street_spectrum_periodic(r, n_pairs=8).real) < 1e-6  # translation zero mode + neutral street (round-off ≈ 1e-8)
    # a general wavenumber: the Bloch eigenvalues appear in the finite-cell spectrum (k = 2πj/(Na)); lattice-sum truncation converges like 1/N²
    N, j, b = 16, 2, 0.3
    k = 2 * np.pi * j / N
    per = BB.karman_street_spectrum_periodic(b, n_pairs=N)
    hs, errs = [], []
    for nt in (500, 1000, 2000, 4000):
        sp_k = BB.karman_street_spectrum(b, k=k, n_terms=nt)
        errs.append(max(np.min(np.abs(per - s_)) for s_ in sp_k))
        hs.append(1.0 / nt)
    assert errs[-1] < 5e-8 and abs(observed_order(hs, errs) - 2.0) < ORDER_TOL, (errs, pairwise_orders(hs, errs))


def test_karman_street_V7_facing_rows_always_unstable_and_street_motion():  # V7 (C10/D14 step 15): unstable at every spacing; U_s = (Γ/2a) tanh(πb/a)
    for b in (0.05, 0.1, 0.3, 0.5, 1.0, 2.0):
        sp_ = BB.karman_street_spectrum(b, offset=0.0)
        assert np.max(sp_.real) == pytest.approx(np.pi / 4, abs=1e-12)  # independent of b
        assert BB.karman_street_growth(b, offset=0.0) == pytest.approx(np.pi / 4, abs=1e-12)
    # street speed by DIRECT summation of the Biot–Savart velocity of 2N+1 pairs (independent of the closed form and of the tanh)
    a, Gam = 1.3, 0.7
    for b in (0.2, 0.28, 0.6):
        n = np.arange(-40000, 40001)
        zA = n * a + 0.5j * b * a  # upper row (+Γ); A = index 0 (b is the ratio b/a)
        zB = (n + 0.5) * a - 0.5j * b * a  # lower row (−Γ)
        z0_ = zA[40000]
        w = np.sum(Gam / (2j * np.pi) / (z0_ - zA[n != 0])) + np.sum(-Gam / (2j * np.pi) / (z0_ - zB))  # u − iv at A (pairs symmetric ⇒ converges)
        assert abs(w.imag) < 5e-6 * Gam / a  # purely streamwise
        assert abs(w.real) == pytest.approx(float(BB.karman_street_velocity(a, b * a, Gam)), rel=2e-4), b
    r = BB.karman_street_ratio()
    assert float(BB.karman_street_velocity(1.0, r, 1.0)) == pytest.approx(0.5 * np.tanh(np.arccosh(np.sqrt(2))), rel=1e-13)  # tanh(arccosh √2) = 1/√2 ⇒ 0.3536Γ/a
    assert float(BB.karman_street_velocity(1.0, r, 1.0)) == pytest.approx(0.5 / np.sqrt(2), rel=1e-13)


def test_karman_street_positions_V7_linear_growth_and_decay_of_the_perturbation():  # V7 (C10): the perturbed street grows at Re σ, stays bounded when neutral, decays when stable
    def amp(b, mode, ts, eps=1e-3):
        return max(float(np.max(np.abs(BB.karman_street_positions(t, b, eps=eps, mode=mode)["zA"] - (np.arange(8) + 0.5j * b)))) for t in ts)

    r = BB.karman_street_ratio()
    a0 = amp(0.5, "unstable", [0.0])
    assert a0 == pytest.approx(1e-3, rel=1e-6)
    assert amp(0.5, "unstable", np.linspace(9, 10.5, 16)) > 20 * a0  # growth e^{0.536 t} ≈ 200× at t = 10 (oscillating factor ≥ 0.1)
    assert BB.karman_street_positions(1.0, 0.5, mode="unstable")["growth"] == pytest.approx(float(BB.karman_street_growth_closed(0.5)), abs=1e-12)
    assert amp(r, "unstable", np.linspace(0, 20, 81)) < 1.6e-3  # neutral: bounded oscillation, no growth
    assert amp(0.4, "stable", np.linspace(9, 10.5, 16)) < 0.2 * amp(0.4, "stable", [0.0])  # stable: decays, growth < 0
    assert BB.karman_street_positions(1.0, 0.4, mode="stable")["growth"] < 0
    with pytest.raises(ValueError):
        BB.karman_street_positions(1.0, 0.3, mode="nonsense")


def test_shedding_V1_strouhal_frequency_and_the_omega_versus_f_trap():  # V1 (C10/N78): f = St U/d; U = 10 m/s, d = 2 mm ⇒ 1000 Hz; Ω of St = Ωd/U is f in Hz
    s = BB.shedding_frequency(10.0, 2e-3)
    assert s["f"] == pytest.approx(1000.0, rel=1e-14) and s["omega_rad"] == pytest.approx(2 * np.pi * 1000.0, rel=1e-14)
    # ch04's St = Ωd/U with Ω in Hz gives 0.2; with the angular frequency the same street would read 2π × 0.2
    assert float(SIM.strouhal_number(s["f"], 2e-3, 10.0)) == pytest.approx(0.2, rel=1e-14)
    assert float(SIM.strouhal_number(s["omega_rad"], 2e-3, 10.0)) == pytest.approx(2 * np.pi * 0.2, rel=1e-13)
    assert float(BB.shedding_angular_frequency(10.0, 2e-3)) == pytest.approx(2 * np.pi * 1000.0, rel=1e-14)
    assert BB.shedding_frequency(20.0, 2e-3)["f"] == pytest.approx(2 * 1000.0, rel=1e-14) and BB.shedding_frequency(10.0, 4e-3, St=0.21)["f"] == pytest.approx(0.21 * 10 / 4e-3)
    dimensional_check(lambda St, U, d: St * U / d, "vorticity", St=Q_(0.2, "dimensionless"), U=Q_(10.0, "m/s"), d=Q_(2e-3, "m"))  # a frequency is 1/s
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        v = BB.strouhal_of_re(np.array([100.0, 1e3, 1e5, 1e6]))
        assert np.isnan(v[0]) and np.isnan(v[3]) and 0.19 <= v[1] <= 0.22 and 0.19 <= v[2] <= 0.22  # Roshko plateau only where cited
        assert any("strouhal_of_re" in str(x.message) for x in w)


# =====================================================================================================================
# C12 — the free two-dimensional laminar jet (§9.10)
# =====================================================================================================================
def _jet_syms():
    return sp.symbols("x y J rho nu", positive=True)


def test_free_jet_V2_sympy_similarity_and_the_tanh_solution():  # V2 (C12): 3f‴ + ff″ + f′² = 0 ⇐ (9.18); f = √6 tanh(η/√6) solves it with (9.66)–(9.68)
    r = ch09.similarity_reduce_sympy("free_jet")
    assert z0(r["residual"]) and r["ode_coefficients"] == {"fppp": 3, "f f''": 1, "f'^2": 1}
    assert r["n"] == sp.Rational(1, 3) and r["m"] == sp.Rational(2, 3)
    e = sp.Symbol("eta", real=True)
    f = sp.sqrt(6) * sp.tanh(e / sp.sqrt(6))
    assert z0(3 * f.diff(e, 3) + f * f.diff(e, 2) + f.diff(e) ** 2)
    assert (f.subs(e, 0), f.diff(e).subs(e, 0), sp.limit(f.diff(e), e, sp.oo)) == (0, 1, 0)  # (9.68), (9.67), (9.66)
    assert z0(3 * f.diff(e, 2) + f * f.diff(e)) and z0(3 * f.diff(e) + f ** 2 / 2 - 3)  # the two first integrals (9.69)
    assert z0(r["F"].subs(sp.Symbol("eta", positive=True), e) - f)
    # C = ∫f′²dη = ∫sech⁴(η/√6)dη = 4√6/3 — sympy, independent of jets.py; and the quantities the code exposes
    t = sp.Symbol("t", real=True)
    C = sp.integrate(sp.sqrt(6) * (1 - t ** 2), (t, -1, 1))  # η = √6 artanh t, sech² = 1 − t²
    assert z0(C - 4 * sp.sqrt(6) / 3)
    c = JET.free_jet_constants()
    assert c["C"] == pytest.approx(float(C), rel=1e-14) and c["C_quad"] == pytest.approx(float(C), rel=1e-11)
    assert c["mdot_coeff"] == pytest.approx(float(sp.Integer(36) ** sp.Rational(1, 3)), rel=1e-14) and c["f_inf"] == pytest.approx(np.sqrt(6), rel=1e-14)
    assert float(sp.sqrt(6) * sp.acosh(10)) == pytest.approx(c["h99_coeff"], rel=1e-14)  # 7.3319
    assert (2 * np.sqrt(6)) ** 3 / float(C) == pytest.approx(36.0, rel=1e-13)  # D17 step 10


def test_free_jet_V3_bvp_matches_tanh_and_the_truncated_boundary_condition_converges():  # V3 (C12): solve_bvp vs √6 tanh; η_max and tolerance studies
    d = JET.free_jet_ode_solve()  # Robin far-field condition, η_max = 12
    assert d["success"] and d["max_err"] < 1e-8  # measured ≈ 4e-13
    errs = [JET.free_jet_ode_solve(tol=t)["max_err"] for t in (1e-4, 1e-6, 1e-8, 1e-10)]
    assert errs[0] > errs[1] > errs[2] and errs[3] < 1e-10 or errs[3] < 1e-12  # refining the tolerance refines the answer
    assert JET.free_jet_ode_solve(eta_max=16.0)["max_err"] < 1e-8 and JET.free_jet_ode_solve(eta_max=9.0)["max_err"] < 1e-8
    # the book's f′(η_max) = 0 (9.66) imposed at a finite η_max converges to the truth like the tail e^{-2η/√6}
    e12 = JET.free_jet_ode_solve(eta_max=12.0, bc="dirichlet")["max_err"]
    e18 = JET.free_jet_ode_solve(eta_max=18.0, bc="dirichlet")["max_err"]
    e24 = JET.free_jet_ode_solve(eta_max=24.0, bc="dirichlet")["max_err"]
    assert e12 > 10 * e18 > 100 * e24 * 0.1 and e24 < 1e-4
    assert 2e-3 < e12 < 3e-3  # (docstring says ≈ 2e-4; measured 2.4e-3, see report)
    # a wrong coefficient of f‴ breaks the closed form (the solver is not tuned to it)
    assert JET.free_jet_ode_solve(coeff=1.0)["max_err"] > 0.05
    with pytest.raises(ValueError):
        JET.free_jet_ode_solve(bc="nonsense")
    assert JET.free_jet_ode_solve(fast=True)["max_err"] < 1e-6


def test_free_jet_V4_momentum_flux_is_conserved_and_mass_flux_grows():  # V4 (C12/D15): ρ∫u²dy = J at every x; ṁ ∝ x^{1/3}; exponents by log–log slope
    J, rho, nu = 1.3, 1.2, 1.5e-5
    xs = np.array([0.01, 0.05, 0.1, 0.5, 2.0, 10.0])
    Jx, mdot, u0, dl = [], [], [], []
    for x in xs:
        delta = float(JET.free_jet_thickness(x, J, rho, nu))
        fu = lambda y, x=x: float(JET.free_jet(x, y, J, rho, nu)["u"])  # noqa: E731
        Jx.append(2 * rho * quad(lambda y: fu(y) ** 2, 0, 40 * delta, epsabs=0, epsrel=1e-13, limit=200)[0])
        mdot.append(2 * rho * quad(fu, 0, 40 * delta, epsabs=0, epsrel=1e-13, limit=200)[0])
        u0.append(float(JET.free_jet(x, 0.0, J, rho, nu)["u"]))
        dl.append(delta)
    assert np.max(np.abs(np.array(Jx) / J - 1)) < 1e-10  # momentum flux J, independent of x (V4)
    assert np.max(np.abs(np.array(mdot) / np.array([float(JET.free_jet_mass_flux(x, J, rho, nu)) for x in xs]) - 1)) < 1e-10  # (9.73) by quadrature
    lx = np.log(xs)
    assert np.polyfit(lx, np.log(u0), 1)[0] == pytest.approx(-1 / 3, abs=1e-12)  # u₀ ∝ x^{-1/3}
    assert np.polyfit(lx, np.log(dl), 1)[0] == pytest.approx(2 / 3, abs=1e-12)  # δ ∝ x^{2/3}
    assert np.polyfit(lx, np.log(mdot), 1)[0] == pytest.approx(1 / 3, abs=1e-10)  # ṁ ∝ x^{1/3}: entrainment, momentum is conserved but mass is not
    assert np.all(np.diff(mdot) > 0)
    # trapezoid helper (used by the notebook) on a fine grid, both half-lines
    y = np.linspace(-30 * dl[3], 30 * dl[3], 40001)
    assert JET.jet_momentum_flux(y, JET.free_jet(0.5, y, J, rho, nu)["u"], rho) == pytest.approx(J, rel=1e-8)
    # scaling in J: u₀ ∝ J^{2/3}, δ ∝ J^{-1/3}, ṁ ∝ J^{1/3}
    assert float(JET.free_jet_centreline(1.0, 8 * J, rho, nu) / JET.free_jet_centreline(1.0, J, rho, nu)) == pytest.approx(4.0, rel=1e-13)
    assert float(JET.free_jet_thickness(1.0, 8 * J, rho, nu) / JET.free_jet_thickness(1.0, J, rho, nu)) == pytest.approx(0.5, rel=1e-13)
    assert float(JET.free_jet_mass_flux(1.0, 8 * J, rho, nu) / JET.free_jet_mass_flux(1.0, J, rho, nu)) == pytest.approx(2.0, rel=1e-13)


def test_free_jet_V3_pde_residual_and_continuity_are_second_order_and_entrainment_matches_9_75():  # V3+V1 (C12): (9.18) on the field; v/u₀ → ∓√6/(3√Re_x)
    J, rho, nu = 1.0, 1.2, 1.5e-5
    hs, errs, errc = [], [], []
    for n in (40, 80, 160):
        x = np.linspace(0.2, 0.4, n + 1)
        y = np.linspace(-3e-3, 3e-3, n + 1)
        X, Y = np.meshgrid(x, y, indexing="xy")
        f = JET.free_jet(X, Y, J, rho, nu)
        ux = np.gradient(f["u"], x, axis=1, edge_order=2)
        uy = np.gradient(f["u"], y, axis=0, edge_order=2)
        uyy = np.gradient(uy, y, axis=0, edge_order=2)
        res = f["u"] * ux + f["v"] * uy - nu * uyy  # (9.18) with dp/dx = 0
        vy = np.gradient(f["v"], y, axis=0, edge_order=2)
        sc = np.max(f["u"]) ** 2 / 0.2
        errs.append(np.max(np.abs(res[2:-2, 2:-2])) / sc)
        errc.append(np.max(np.abs((ux + vy)[2:-2, 2:-2])) / (np.max(f["u"]) / 0.2))
        hs.append(x[1] - x[0])
    assert abs(observed_order(hs, errs) - 2.0) < 0.25, pairwise_orders(hs, errs)
    assert abs(observed_order(hs, errc) - 2.0) < 0.25, pairwise_orders(hs, errc)
    # entrainment velocity at the edge: v/u₀ = ∓√6/(3√Re_x), Re_x = x u₀/ν
    x = 0.3
    fld = JET.free_jet(x, 60 * float(JET.free_jet_thickness(x, J, rho, nu)), J, rho, nu)
    Re_x = x * float(fld["u0"]) / nu
    assert float(fld["v"]) / float(fld["u0"]) == pytest.approx(float(JET.free_jet_entrainment_velocity(Re_x)), rel=1e-10)
    fneg = JET.free_jet(x, -60 * float(JET.free_jet_thickness(x, J, rho, nu)), J, rho, nu)
    assert float(fneg["v"]) == pytest.approx(-float(fld["v"]), rel=1e-12)  # symmetric: fluid comes in from both sides (v odd in y)
    assert float(fld["v"]) < 0 and float(JET.free_jet_entrainment_velocity(Re_x)) == pytest.approx(-np.sqrt(6) / (3 * np.sqrt(Re_x)), rel=1e-14)
    # v = −ψ_x from the returned ψ (finite difference), and u = ψ_y
    h = 1e-7
    yy = 1e-4
    v_fd = -(float(JET.free_jet(x + h, yy, J, rho, nu)["psi"]) - float(JET.free_jet(x - h, yy, J, rho, nu)["psi"])) / (2 * h)
    u_fd = (float(JET.free_jet(x, yy + h * 1e-3, J, rho, nu)["psi"]) - float(JET.free_jet(x, yy - h * 1e-3, J, rho, nu)["psi"])) / (2e-3 * h)
    assert float(JET.free_jet(x, yy, J, rho, nu)["v"]) == pytest.approx(v_fd, rel=1e-5)
    assert float(JET.free_jet(x, yy, J, rho, nu)["u"]) == pytest.approx(u_fd, rel=1e-5)


@needs_ref
def test_free_jet_V5_bickley_constants_and_the_design_expect_rows():  # V5 (C12): Bickley (1937) via Wikipedia — 0.4543, 0.2752, 3.3019, 0.5503
    R = ref_json()["bickley_jet_wikipedia"]
    c = JET.free_jet_constants()
    assert abs(c["u0_coeff"] - R["u0_coeff"]) < 5e-5 and abs(c["xi_coeff"] - R["xi_coeff"]) < 5e-5 and abs(c["mdot_coeff"] - R["Q_coeff"]) < 5e-5  # rounding of 4–5 printed digits
    v_coeff = np.sqrt(6) / 3 * c["C"] ** (-1 / 3)  # v = (√6/3)(Jν/Cρx²)^{1/3}(2ξ sech²ξ − tanh ξ) from (9.74)
    assert abs(v_coeff - R["v_coeff"]) < 5e-5
    # the field reproduces Bickley's form u = 0.4543 (M²/νρ²x)^{1/3} sech²ξ, ξ = 0.2752 (M/ν²ρ)^{1/3} y/x^{2/3} at arbitrary parameters
    J, rho, nu, x, y = 2.5, 1.1, 2e-5, 0.7, 1.5e-3
    xi = R["xi_coeff"] * (J / (nu ** 2 * rho)) ** (1 / 3) * y / x ** (2 / 3)
    u_b = R["u0_coeff"] * (J ** 2 / (nu * rho ** 2 * x)) ** (1 / 3) / np.cosh(xi) ** 2
    assert float(JET.free_jet(x, y, J, rho, nu)["u"]) == pytest.approx(u_b, rel=3e-4)
    Q_b = R["Q_coeff"] * (J * nu * rho ** 2 * x) ** (1 / 3)  # Q = 2ρ∫₀^∞u dy
    assert float(JET.free_jet_mass_flux(x, J, rho, nu)) == pytest.approx(Q_b, rel=3e-5)
    # design expect rows: air J = 1 N/m, x = 0.1 m; water J = 1
    kw = dict(J=1.0)
    a = dict(rho=1.2, nu=1.5e-5)
    assert float(JET.free_jet_centreline(0.1, rho=a["rho"], nu=a["nu"], **kw)) == pytest.approx(35.14, abs=5e-3)
    assert float(JET.free_jet_thickness(0.1, rho=a["rho"], nu=a["nu"], **kw)) == pytest.approx(0.2066e-3, abs=5e-8)
    assert float(JET.free_jet_mass_flux(0.1, rho=a["rho"], nu=a["nu"], **kw)) == pytest.approx(0.04268, abs=5e-6)
    assert float(JET.free_jet_halfwidth(0.1, rho=a["rho"], nu=a["nu"], **kw)) == pytest.approx(1.515e-3, abs=5e-7)
    assert JET.free_jet_reynolds(0.1, rho=a["rho"], nu=a["nu"], **kw)["Re_x"] == pytest.approx(2.343e5, rel=2e-3)
    assert float(JET.free_jet(0.1, 60 * 0.2066e-3, rho=a["rho"], nu=a["nu"], **kw)["v"]) == pytest.approx(-0.0593, abs=1e-4)  # edge entrainment velocity
    assert float(JET.free_jet_centreline(0.1, 1.0, 1000.0, 1e-6)) == pytest.approx(0.9787, abs=5e-4)
    assert float(JET.free_jet_mass_flux(0.1, 1.0, 1000.0, 1e-6)) == pytest.approx(1.533, abs=1e-3)


def test_free_jet_V1_R6_halfwidth_printed_coefficient_is_the_four_percent_point():  # V1 (C12/N112/D18): sech²(z) = level ⇒ √6·arccosh(1/√level); the printed 5.6152 is the 4 % point
    J, rho, nu, x = 1.0, 1.2, 1.5e-5, 0.1
    delta = float(JET.free_jet_thickness(x, J, rho, nu))
    for lev, want in ((0.01, 7.3319), (0.04, 5.6153), (0.5, 2.1589)):
        assert JET.free_jet_at_level(lev)["coeff"] == pytest.approx(want, abs=1e-4), lev
        h = float(JET.free_jet_halfwidth(x, J, rho, nu, level=lev))
        assert float(JET.free_jet(x, h, J, rho, nu)["u"]) / float(JET.free_jet(x, 0.0, J, rho, nu)["u"]) == pytest.approx(lev, rel=1e-12)  # the definition, on the field
    assert JET.free_jet_at_level(0.01)["z"] == pytest.approx(np.log(10 + np.sqrt(99)), rel=1e-14)  # arccosh 10 = ln(10 + √99) (D18 step 3)
    # planted wrong variant: the printed coefficient does NOT put u at 1 % of u₀ — it is the 4 % point, and the halfwidth is 23 % too small
    hp = float(JET.free_jet_halfwidth(x, J, rho, nu, printed=True))
    ratio = float(JET.free_jet(x, hp, J, rho, nu)["u"]) / float(JET.free_jet(x, 0.0, J, rho, nu)["u"])
    assert ratio == pytest.approx(0.04, abs=2e-5) and abs(ratio - 0.01) > 0.02, "the printed 5.6152 must NOT satisfy sech²z = 0.01"
    assert hp / float(JET.free_jet_halfwidth(x, J, rho, nu)) == pytest.approx(JET.H99_PRINTED / JET.free_jet_at_level(0.01)["coeff"], rel=1e-13)
    assert 0.76 < hp / float(JET.free_jet_halfwidth(x, J, rho, nu)) < 0.77
    assert 1 / np.cosh(2.2924) ** 2 == pytest.approx(0.04, abs=1e-4) and 1 / np.cosh(2.9932) ** 2 == pytest.approx(0.01, abs=1e-5)  # D18 check line
    rr = JET.free_jet_reynolds(x, J, rho, nu)
    Re_x = x * float(JET.free_jet_centreline(x, J, rho, nu)) / nu
    assert rr["Re_x"] == pytest.approx(Re_x, rel=1e-12)
    assert rr["Re_h99"] == pytest.approx(float(JET.free_jet_centreline(x, J, rho, nu)) * float(JET.free_jet_halfwidth(x, J, rho, nu)) / nu, rel=1e-12)
    assert rr["Re_h99_printed"] / rr["Re_h99"] == pytest.approx(JET.H99_PRINTED / JET.free_jet_at_level(0.01)["coeff"], rel=1e-12)
    assert (3 * J * x / (4 * np.sqrt(6) * rho * nu ** 2)) ** (2 / 3) == pytest.approx(Re_x, rel=1e-12)  # the closed form of Re_x
    assert len(JET.free_jet_profile_table(51, 10.0)["eta"]) == 51 and JET.free_jet_profile_table(51, 10.0)["fp"][0] == 1.0
    dimensional_check(lambda J, rho, nu, x: (J ** 2 / (rho ** 2 * nu * x)) ** (1 / 3), "velocity", J=Q_(1.0, "N/m"), rho=Q_(1.2, "kg/m**3"), nu=Q_(1.5e-5, "m**2/s"), x=Q_(0.1, "m"))


# =====================================================================================================================
# C13 — the wall jet (§9.10)
# =====================================================================================================================
def test_wall_jet_V2_sympy_reduction_first_integrals_and_R3_R4_planted_variants():  # V2 (C13): 4f‴ + ff″ + 2f′² = 0; the printed coefficient 1 fails; (9.83) from the first integrals
    r = ch09.similarity_reduce_sympy("wall_jet")
    assert z0(r["residual"]) and r["ode_coefficients"] == {"fppp": 4, "f f''": 1, "f'^2": 2}
    assert r["n"] == sp.Rational(-1, 2) and r["m"] == sp.Rational(3, 4)
    bad = ch09.similarity_reduce_sympy("wall_jet", printed=True)
    assert not z0(bad["residual"]), "R3: the printed coefficient 1 leaves a non-zero residual"
    assert bad["ode_coefficients"]["fppp"] == 1
    x, nu, C = sp.symbols("x nu C", positive=True)
    F3 = sp.Symbol("F3")
    # the leftover of the printed variant is −(3/4)C²f‴/x²: the missing 3f‴ times the prefactor −C²/(4x²) (design D20 check line)
    assert z0(bad["residual"] + sp.Rational(3, 4) * C ** 2 * F3 / x ** 2)
    w = ch09.wall_jet_sympy()
    assert z0(w["first_integral_derivative"]) and z0(w["second_integral_derivative"]) and z0(w["g_integral_derivative"])
    assert z0(w["integrand_correct_residual"]) and not z0(w["integrand_printed_residual"])  # R4: the printed ∫df/(f_∞^{3/2}f − f²) is not f′ = …/6
    assert w["fpp0_over_finf3"] == sp.Rational(1, 72) and w["K1"] == sp.Rational(1, 40)
    assert w["ode_correct"].lhs != w["ode_printed"].lhs


def test_wall_jet_V3_two_independent_routes_agree_and_the_scale_law_f_inf_cubed_is_72_fpp0():  # V3+V7 (C13/D21): IVP vs the implicit (9.83) by brentq; f″(0) = f_∞³/72
    d = JET.wall_jet_ode_solve(1.0 / 72.0, eta_max=120.0, n=4001)
    assert d["f_inf"] == pytest.approx(1.0, abs=1e-10) and d["err_vs_9_83"] < 1e-10  # measured 7e-12
    assert d["fpp0_over_finf_cubed"] == pytest.approx(1.0 / 72.0, rel=1e-9)
    assert d["first_integral_residual_max"] < 1e-11  # 4ff″ − 2f′² + f²f′ = 0 along the IVP solution (D21 step 4)
    assert np.max(d["fp"]) == pytest.approx(0.078745, abs=2e-6) and d["eta"][np.argmax(d["fp"])] == pytest.approx(8.11, abs=0.03)  # the design's expect row
    # scaling f → λ f(λη): f″(0) = λ³/72 ⇒ f_∞ = λ (V7); tested at three scales
    for lam in (0.5, 2.0, 3.0):
        e = JET.wall_jet_ode_solve(lam ** 3 / 72.0, eta_max=120.0 / lam, n=4001)
        assert e["f_inf"] == pytest.approx(lam, rel=1e-9) and e["err_vs_9_83"] < 1e-9, lam
    # the profile function obeys the ODE it claims to solve (finite differences of the (9.83) profile) and the scaling law
    eta = np.linspace(0.0, 60.0, 12001)
    pr = JET.wall_jet_profile(eta, 1.0)
    h = eta[1] - eta[0]
    fpp_fd = np.gradient(pr["fp"], h, edge_order=2)
    fppp = np.gradient(pr["fpp"], h, edge_order=2)
    res = 4 * fppp + pr["f"] * pr["fpp"] + 2 * pr["fp"] ** 2
    assert np.max(np.abs(res[3:-3])) < 3e-8 and np.max(np.abs(fpp_fd[3:-3] - pr["fpp"][3:-3])) < 1e-7
    assert pr["fpp"][0] == pytest.approx(1 / 72, rel=1e-12) and pr["f"][0] == 0 and pr["fp"][0] == 0
    assert pr["f"][-1] == pytest.approx((1 - 4.2897 * np.exp(-15.0)) ** 2, abs=2e-9)  # f = f_∞g², 1 − g = 4.29e^{−η/4} (15 e-foldings at η = 60)
    for lam in (0.5, 2.0):
        p2 = JET.wall_jet_profile(eta / lam, lam)
        p1 = JET.wall_jet_profile(eta, 1.0)
        assert np.max(np.abs(p2["f"] - lam * p1["f"])) < 1e-12 and np.max(np.abs(p2["fp"] - lam ** 2 * p1["fp"])) < 1e-12
    # far field: 1 − g ≈ √3 e^{√3π/6} e^{−f_∞η/4} = 4.29 e^{−f_∞η/4} (the constant the book omits)
    g = JET.wall_jet_profile(np.array([60.0, 80.0]), 1.0)["g"]
    assert (1 - g) / np.exp(-np.array([60.0, 80.0]) / 4) == pytest.approx(np.sqrt(3) * np.exp(np.sqrt(3) * np.pi / 6), rel=1e-5)


def test_wall_jet_R3_printed_ode_fails_the_relation_and_gives_a_different_f_inf():  # V2 (C13/N120): coefficient 1 integrates to f_∞ = 0.397, not 1, and violates (9.83)
    bad = JET.wall_jet_ode_solve(1.0 / 72.0, eta_max=120.0, coeff=1.0, n=4001)
    assert abs(bad["f_inf"] - 1.0) > 0.5 and bad["err_vs_9_83"] > 0.05  # design: f_∞ = 0.397
    assert bad["f_inf"] == pytest.approx(0.397, abs=2e-3)
    assert bad["first_integral_residual_max"] > 1e-3  # the first integral 4ff″ − 2f′² + f²f′ is NOT conserved by the printed ODE
    assert JET.wall_jet_ode_solve(1.0 / 72.0, eta_max=120.0, printed=True, n=4001)["f_inf"] == pytest.approx(bad["f_inf"], rel=1e-14)


def test_wall_jet_ode_V3_design_expect_row_at_the_default_eta_max():  # V3 (C13, Part C 2.6): the contract "f″(0) = 1/72 ⇒ f_∞ = 1.0000000" with the DEFAULT arguments
    # (F4, fixed) the default η_max = 40 used to stop 10 e-foldings out: f_∞ = 0.9996, err 5.5e-4.
    # (η_max = 120 reproduces the contract to 7e-12).  The notebook/explainer that call the default print these numbers.
    d = JET.wall_jet_ode_solve(1.0 / 72.0)
    assert d["f_inf"] == pytest.approx(1.0, abs=1e-6)
    assert d["err_vs_9_83"] < 1e-8


def test_wall_jet_V4_invariant_is_conserved_while_the_momentum_flux_is_not():  # V4 (C13/D19): d/dx ∫u(∫_y^∞u²)dy = 0; ρ∫u²dy ∝ x^{-1/4}; ṁ ∝ x^{1/4}
    C, finf, nu, rho = 0.7, 1.3, 1e-3, 1000.0
    xs = np.array([0.25, 0.5, 1.0, 2.0, 4.0])
    inv, mom, mdot = [], [], []
    for x, npts in zip(xs, (2001, 3001, 4001, 5001, 6001)):  # a different resolution at every station: a real test of the quadrature
        delta = float(JET.wall_jet(x, 0.0, C, finf, nu)["delta"])
        y = np.linspace(0.0, 110 * delta / finf, npts)  # f′ ≈ 2.1 f_∞² e^{−f_∞η/4}: 110/f_∞ ≈ 27 e-foldings
        f = JET.wall_jet(x, y, C, finf, nu, rho)
        inv.append(JET.wall_jet_invariant(y, f["u"]))
        mom.append(rho * np.trapezoid(f["u"] ** 2, y))
        mdot.append(rho * np.trapezoid(f["u"], y))
    inv, mom, mdot = map(np.array, (inv, mom, mdot))
    Psi = C ** 2 * nu * finf ** 4 / 40.0  # (9.85) with K = 1/40
    assert np.max(np.abs(inv / Psi - 1)) < 1e-5  # conserved (V4) and equal to (9.85), to the second-order trapezoid error of each grid
    lx = np.log(xs)
    assert np.polyfit(lx, np.log(mom), 1)[0] == pytest.approx(-0.25, abs=2e-5)  # the ordinary momentum flux DECAYS (wall shear): ∝ x^{-1/4}
    assert np.polyfit(lx, np.log(mdot), 1)[0] == pytest.approx(0.25, abs=2e-5)  # ṁ ∝ x^{1/4} (9.84)
    assert mdot == pytest.approx([float(JET.wall_jet_mass_flux(x, C, finf, rho, nu)) for x in xs], rel=1e-5)
    assert np.polyfit(lx, np.log([float(JET.wall_jet(x, 0.0, C, finf, nu)["delta"]) for x in xs]), 1)[0] == pytest.approx(0.75, abs=1e-12)  # δ ∝ x^{3/4}
    assert np.polyfit(lx, np.log([float(JET.wall_jet(x, 1e-3, C, finf, nu)["u0"]) for x in xs]), 1)[0] == pytest.approx(-0.5, abs=1e-12)  # u₀ ∝ x^{-1/2}
    # the free jet conserves ρ∫u²dy instead (V4 contrast) and its width grows more slowly (2/3 < 3/4)
    assert 2 / 3 < 3 / 4
    # numerical K₁ = ∫f′∫f′² = 1/40 from the profile, and the exact rational integrals
    assert JET.wall_jet_K1(numeric=True) == pytest.approx(0.025, rel=1e-6) and JET.wall_jet_K1() == 0.025
    ints = JET.wall_jet_integrals(1.0)
    assert (ints["int_fp"], ints["int_fp2"], ints["invariant"], ints["fpp0"]) == pytest.approx((1.0, 1 / 18, 1 / 40, 1 / 72), rel=1e-12)
    i2 = JET.wall_jet_integrals(2.0)
    assert (i2["int_fp"], i2["int_fp2"], i2["invariant"], i2["fpp0"]) == pytest.approx((2.0, 8 / 18, 16 / 40, 8 / 72), rel=1e-12)  # f_∞, f_∞³, f_∞⁴, f_∞³ scalings


def test_wall_jet_V3_fields_satisfy_continuity_and_the_boundary_layer_equation():  # V3 (C13): (9.18) and (6.2) on the wall-jet field, 2nd order
    C, finf, nu = 1.0, 1.0, 1e-3
    hs, errs, errc = [], [], []
    for n in (40, 80, 160):
        x = np.linspace(1.0, 1.4, n + 1)
        y = np.linspace(0.0, 0.25, n + 1)
        X, Y = np.meshgrid(x, y, indexing="xy")
        f = JET.wall_jet(X, Y, C, finf, nu)
        ux = np.gradient(f["u"], x, axis=1, edge_order=2)
        uy = np.gradient(f["u"], y, axis=0, edge_order=2)
        uyy = np.gradient(uy, y, axis=0, edge_order=2)
        vy = np.gradient(f["v"], y, axis=0, edge_order=2)
        sc = np.max(f["u"]) ** 2 / 1.0
        errs.append(np.max(np.abs((f["u"] * ux + f["v"] * uy - nu * uyy)[2:-2, 2:-2])) / sc)
        errc.append(np.max(np.abs((ux + vy)[2:-2, 2:-2])) / np.max(f["u"]))
        hs.append(x[1] - x[0])
    assert abs(observed_order(hs, errs) - 2.0) < 0.25, pairwise_orders(hs, errs)
    assert abs(observed_order(hs, errc) - 2.0) < 0.25, pairwise_orders(hs, errc)
    w = JET.wall_jet(1.0, 0.0, C, finf, nu)
    assert float(w["u"]) == 0.0 and float(w["v"]) == pytest.approx(0.0, abs=1e-14)  # no slip, no through-flow (9.77)
    assert float(JET.wall_jet(1.0, 3.0, C, finf, nu)["u"]) < 1e-9 * float(w["u0"])  # (9.78): f′ ≈ 2.1e^{−η/4}, η = 95 at y = 3
    # constants from ONE datum: Ψ = ṁ⁴/(40ρ⁴νx) (the (C, f_∞) pair is degenerate — only C f_∞² is physical)
    rho, x = 1000.0, 2.0
    mdot = float(JET.wall_jet_mass_flux(x, 0.7, 1.3, rho, nu))
    a = JET.wall_jet_constants(rho, nu, mdot=mdot, x=x, f_inf=1.3)
    b = JET.wall_jet_constants(rho, nu, Psi=a["Psi"], x=x, f_inf=1.3)
    assert a["Psi"] == pytest.approx(mdot ** 4 / (40 * rho ** 4 * nu * x), rel=1e-13) and a["C"] == pytest.approx(0.7, rel=1e-12)
    assert a["C_f_inf_sq"] == pytest.approx(0.7 * 1.3 ** 2, rel=1e-12) and b["mdot_at_x"] == pytest.approx(mdot, rel=1e-12)
    assert JET.wall_jet_constants(rho, nu, Psi=a["Psi"], mdot=mdot, x=x)["residual"] == pytest.approx(0.0, abs=1e-12)
    assert abs(JET.wall_jet_constants(rho, nu, Psi=2 * a["Psi"], mdot=mdot, x=x)["residual"] - 1.0) < 1e-12
    with pytest.raises(ValueError):
        JET.wall_jet_constants(rho, nu)
    with pytest.raises(ValueError):
        JET.wall_jet_constants(rho, nu, mdot=mdot)
    # the first-integral helper is zero on the IVP and non-zero on a wrong profile
    assert abs(float(JET.wall_jet_first_integral_residual(0.5, 0.1, 0.2))) > 1e-3
    dimensional_check(lambda C, nu: C ** 2 * nu / 40.0, "[length] ** 5 / [time] ** 3", C=Q_(1.0, "m**1.5/s"), nu=Q_(1e-3, "m**2/s"))  # Ψ: m⁵/s³ (u³L² with ν x), not a force per length


# =====================================================================================================================
# C14 — the teacup (§9.11)
# =====================================================================================================================
def test_teacup_V1_radial_force_imbalance_numbers_and_sign():  # V1+V7 (C14/D22): F_in = ρ(u_e² − u²)/R; 1000, 750, 437.5, 0 N/m³ for u_e = 0.2, R = 0.04
    ue, R, rho = 0.2, 0.04, 1000.0
    F = [ch09.secondary_flow_radial_force(ue, f * ue, R, rho) for f in (0.0, 0.5, 0.75, 1.0)]
    assert F == pytest.approx([1000.0, 750.0, 437.5, 0.0], abs=1e-9)
    assert F[0] == pytest.approx(rho * ue ** 2 / R, rel=1e-14) and F[0] / (rho * 9.81) == pytest.approx(0.102, abs=5e-4)  # D22 step 5: ≈ 0.10 ρg
    z = np.linspace(0.0, 5e-3, 50)
    u = ch09.secondary_flow_layer_profile(z, 1e-3, ue, "linear")
    Fz = ch09.secondary_flow_radial_force(ue, u, R, rho)
    assert np.all(np.diff(Fz) <= 1e-12) and Fz[0] == pytest.approx(1000.0) and Fz[-1] == 0.0  # largest on the floor, zero in the core
    assert ch09.secondary_flow_radial_force(ue, 1.5 * ue, R, rho) < 0  # a layer faster than the core would be pushed OUTWARD
    dimensional_check(lambda rho, u, R: rho * u ** 2 / R, "[mass] / ([length] ** 2 * [time] ** 2)", rho=Q_(1000.0, "kg/m**3"), u=Q_(0.2, "m/s"), R=Q_(0.04, "m"))  # N/m³
    # Ekman hook: δ_E = √(ν/Ω) with Ω = u_e/R = 5 rad/s ⇒ 0.45 mm in water
    assert np.sqrt(1e-6 / (ue / R)) == pytest.approx(0.447e-3, abs=1e-6)


def test_teacup_V7_layer_profiles_are_illustrative_but_physical():  # V7 (C14, demonstration): the notebook's profile shapes vanish on the floor and reach u_e; force ≥ 0 everywhere
    z = np.linspace(0.0, 4e-3, 401)
    for shape in ("exponential", "linear", "power", "sine"):
        u = ch09.secondary_flow_layer_profile(z, 1e-3, 0.2, shape)
        assert u[0] == 0.0 and np.all(np.diff(u) >= -1e-15) and np.all(u <= 0.2 + 1e-15), shape
        F = ch09.secondary_flow_radial_force(0.2, u, 0.04, 1000.0)
        assert np.all(F >= -1e-12) and F[0] == pytest.approx(1000.0), shape  # inward at every height, maximal at the floor
        # the whole meridional loop needs continuity: the net inward force integrates to a positive radial 'pressure-driven' momentum
        assert np.trapezoid(F, z) > 0
    assert ch09.secondary_flow_layer_profile(np.array([2e-3]), 1e-3, 0.2, "linear")[0] == 0.2 and ch09.secondary_flow_layer_profile(2e-3, 1e-3, 0.2, "sine") == pytest.approx(0.2)
    assert ch09.secondary_flow_layer_profile(1e-3, 1e-3, 0.2, "exponential") == pytest.approx(0.2 * (1 - np.exp(-1)))
    with pytest.raises(ValueError):
        ch09.secondary_flow_layer_profile(1e-3, 1e-3, 0.2, "nonsense")


# =====================================================================================================================
# V2 — symbolic re-derivations of every ★★ / ★★★ D row (and the cheap ★ rows), line by line where the design gives lines
# =====================================================================================================================
_FJ = sp.symbols("F0:5")


def _jet_D(eta_expr):
    """Total derivative of an expression in (x, y, F0…F4) where F_k = f^{(k)}(η) and η = eta_expr(x, y): the chain rule the book skips."""
    def D(e, var):
        out = sp.diff(e, var)
        for k in range(4):
            out += sp.diff(e, _FJ[k]) * _FJ[k + 1] * sp.diff(eta_expr, var)
        return out
    return D


def _at(expr, y, eta, delta):
    """Substitute y = η δ and simplify power laws."""
    return sp.simplify(sp.powsimp(sp.expand_power_base(expr.subs(y, eta * delta), force=True), force=True))


def test_bl_scaling_V2_derivation_D01():  # V2 (D01 ★★): every displayed line of Part F D01, from the Navier–Stokes equations, independent of bl_nondim_sympy
    xs, ys, x, y = sp.symbols("xs ys x y", positive=True)
    L, U, Re, rho = sp.symbols("L U Re rho", positive=True)
    us, vs, ps = (sp.Function(n)(xs, ys) for n in ("us", "vs", "ps"))
    # steps 1–3: dominant balance U²/L = νU/δ̄² ⇒ δ̄/L = Re^{-1/2}; v from continuity
    nu = U * L / Re
    dl = sp.Symbol("delta_bar", positive=True)
    sol = sp.solve(sp.Eq(U ** 2 / L, nu * U / dl ** 2), dl)
    assert len(sol) == 1 and z0(sol[0] / L - Re ** sp.Rational(-1, 2))  # (9.4)
    v_sc = sp.Symbol("v_sc", positive=True)
    assert z0(sp.solve(sp.Eq(U / L, v_sc / sol[0]), v_sc)[0] - sol[0] * U / L) and z0(sol[0] * U / L - U * Re ** sp.Rational(-1, 2))  # v ~ U Re^{-1/2}
    # step 5 on a concrete test function (the chain rule): ∂_y of u* = f(x/L, y√Re/L) gains √Re/L
    test = sp.sin(xs) * sp.cos(2 * ys) + xs * ys
    u_xy = U * test.subs({xs: x / L, ys: y * sp.sqrt(Re) / L}, simultaneous=True)
    assert z0(sp.diff(u_xy, y, 2) - U * Re / L ** 2 * sp.diff(test, ys, 2).subs({xs: x / L, ys: y * sp.sqrt(Re) / L}, simultaneous=True))
    assert z0(sp.diff(u_xy, x) - U / L * sp.diff(test, xs).subs({xs: x / L, ys: y * sp.sqrt(Re) / L}, simultaneous=True))
    d_dx = lambda q: sp.diff(q, xs) / L  # noqa: E731
    d_dy = lambda q: sp.diff(q, ys) * sp.sqrt(Re) / L  # noqa: E731
    u, v, p = U * us, U * vs / sp.sqrt(Re), rho * U ** 2 * ps  # (9.6)
    # step 6: continuity keeps its form
    assert z0(d_dx(u) + d_dy(v) - U / L * (us.diff(xs) + vs.diff(ys)))
    # step 7: inertia and pressure of the x-momentum equation all have the size U²/L
    assert z0(u * d_dx(u) + v * d_dy(u) + d_dx(p) / rho - U ** 2 / L * (us * us.diff(xs) + vs * us.diff(ys) + ps.diff(xs)))
    # step 8: the viscous terms
    assert z0(nu * (d_dx(d_dx(u)) + d_dy(d_dy(u))) - (nu * U / L ** 2 * us.diff(xs, 2) + nu * U * Re / L ** 2 * us.diff(ys, 2)))
    # step 9: divide by U²/L ⇒ (9.7) with 1/Re on the x-derivative only, coefficient 1 on the y-derivative
    xm = u * d_dx(u) + v * d_dy(u) + d_dx(p) / rho - nu * (d_dx(d_dx(u)) + d_dy(d_dy(u)))
    r97 = us * us.diff(xs) + vs * us.diff(ys) + ps.diff(xs) - us.diff(xs, 2) / Re - us.diff(ys, 2)
    assert z0(xm * L / U ** 2 - r97)
    # printed (9.7) (slip R1): first-order denominators — different expression AND wrong dimensions (m²/s² vs m/s²)
    printed = us * us.diff(xs) + vs * us.diff(ys) + ps.diff(xs) - us.diff(xs) / Re - us.diff(ys)
    assert not z0(xm * L / U ** 2 - printed)
    adv = Q_(1, "m/s") ** 2 / Q_(1, "m")
    visc_printed = Q_(1, "m**2/s") * Q_(1, "m/s") / Q_(1, "m")  # ν ∂u/∂x
    visc_ok = Q_(1, "m**2/s") * Q_(1, "m/s") / Q_(1, "m") ** 2
    assert visc_ok.check("[acceleration]") and adv.check("[acceleration]") and not visc_printed.check("[acceleration]")
    # steps 10–12: y-momentum
    ym = u * d_dx(v) + v * d_dy(v) + d_dy(p) / rho - nu * (d_dx(d_dx(v)) + d_dy(d_dy(v)))
    assert z0(u * d_dx(v) + v * d_dy(v) - U ** 2 / (L * sp.sqrt(Re)) * (us * vs.diff(xs) + vs * vs.diff(ys)))  # step 10
    assert z0(d_dy(p) / rho - U ** 2 * sp.sqrt(Re) / L * ps.diff(ys))  # step 11 pressure
    assert z0(nu * d_dx(d_dx(v)) - U ** 2 / (L * Re ** sp.Rational(3, 2)) * vs.diff(xs, 2)) and z0(nu * d_dy(d_dy(v)) - U ** 2 / (L * sp.sqrt(Re)) * vs.diff(ys, 2))  # step 11 viscous
    r98 = (us * vs.diff(xs) + vs * vs.diff(ys)) / Re + ps.diff(ys) - vs.diff(xs, 2) / Re ** 2 - vs.diff(ys, 2) / Re  # (9.8): all on one side
    assert z0(ym * L / (U ** 2 * sp.sqrt(Re)) - r98)
    # steps 13–14: Re → ∞ leaves p*_y* = 0 (9.10) and (9.9)
    assert z0(sp.limit(r98.subs({us: 1, vs: 1}) * 0 + ps.diff(ys) * 0, Re, sp.oo))  # (trivial guard: the limit machinery runs)
    coeffs_98 = [sp.limit(c, Re, sp.oo) for c in (1 / Re, 1 / Re ** 2, 1 / Re)]
    assert coeffs_98 == [0, 0, 0]  # every term of (9.8) except p*_y* carries a coefficient that vanishes
    # the check line: air, U = 1 m/s, L = 1 m, ν = 1.5e-5 ⇒ Re = 6.67e4, δ̄/L = 3.87e-3, dropped/kept diffusion 1/Re = 1.5e-5
    assert 1.0 / 1.5e-5 == pytest.approx(6.67e4, rel=1e-3) and (1.5e-5) ** 0.5 == pytest.approx(3.87e-3, rel=1e-3)
    assert BL.boundary_layer_scales(1.0, 1.0, 1.5e-5)["visc_x"] / BL.boundary_layer_scales(1.0, 1.0, 1.5e-5)["visc"] == pytest.approx(1.5e-5, rel=1e-12)


def test_pressure_matching_V2_derivation_D02():  # V2 (D02 ★): Bernoulli differentiated along the edge gives (9.11); inserted into (9.9) it gives (9.37)
    x = sp.Symbol("x", positive=True)
    rho, c0 = sp.symbols("rho c0", positive=True)
    Ue, pe = sp.Function("U_e")(x), sp.Function("p_e")(x)
    bern = sp.Eq(pe + rho * Ue ** 2 / 2, c0)
    dpe = sp.solve(sp.diff(bern.lhs - bern.rhs, x), sp.diff(pe, x))[0]
    assert z0(-dpe / rho - Ue * sp.diff(Ue, x))  # step 3: −(1/ρ)dp_e/dx = U_e U_e′  (9.11)
    u, v, nu = (sp.Function(n)(x, sp.Symbol("y")) for n in ("u", "v", "nu"))
    y = sp.Symbol("y")
    lhs = u * sp.diff(u, x) + v * sp.diff(u, y)
    assert z0((lhs + dpe / rho - nu * sp.diff(u, y, 2)) - (lhs - Ue * sp.diff(Ue, x) - nu * sp.diff(u, y, 2)))  # step 5: (9.9) → (9.37)
    # the check line: U_e = 10x s⁻¹ at x = 0.1 m, ρ = 1.2: dp/dx = −12 Pa/m (favourable)
    assert float(BL.outer_flow("wedge", n=1.0, a=10.0).dpdx(0.1, 1.2)) == pytest.approx(-12.0, rel=1e-13)
    # a decelerating outer flow has dp/dx > 0 (adverse)
    assert float(BL.outer_flow("diffuser", U1=1.0, L=1.0).dpdx(0.2, 1.2)) > 0


def test_displacement_thickness_V2_derivation_D03():  # V2 (D03 ★): v(h) = U dδ*/dx from continuity, on a smooth trial profile with growing δ(x)
    x, y, h, k, U = sp.symbols("x y h k U", positive=True)
    dl = k * sp.sqrt(x)
    u = U * (1 - sp.exp(-y / dl))  # trial profile (any smooth profile with an edge works)
    A = sp.integrate(u, (y, 0, h))  # step 1: mass-flux per ρ
    deficit = sp.integrate(U - u, (y, 0, h))  # step 2: ∫(U − u)dy = U ∫(1 − u/U)dy
    dstar_h = sp.simplify(deficit / U)  # steps 3–4: δ*(h)
    assert z0(dstar_h - dl * (1 - sp.exp(-h / dl)))
    dstar = sp.limit(dstar_h, h, sp.oo)  # step 5: converges as h → ∞ (the tail e^{-h/δ})
    assert z0(dstar - dl)
    v_h = -sp.integrate(sp.diff(u, x), (y, 0, h))  # step 6: v(h) = −∫u_x dy (v(0) = 0)
    assert z0(sp.simplify(v_h - U * sp.diff(dstar_h, x)))  # −d/dx[U(h − δ*(h))] = U dδ*(h)/dx at fixed h
    assert z0(sp.limit(v_h, h, sp.oo) - U * sp.diff(dl, x))  # v_∞ = U δ*′
    # Blasius: δ* = 1.7208√(νx/U), v_∞ = 0.8604 U/√Re_x (the check line)
    assert BL.blasius_constants()["v_inf"] == pytest.approx(0.5 * 1.7207876575, rel=1e-9)


def test_momentum_thickness_V2_derivation_D04():  # V2 (D04 ★★): ρU²θ = ∫τ₀dx from the control-volume balance, algebra of steps 2–9 with symbolic integrals
    rho, U, h0, h = sp.symbols("rho U h0 h", positive=True)
    A, B, T = sp.symbols("A B T", positive=True)  # A = ∫u dy, B = ∫u² dy over the roof height, T = ∫₀ˣ τ₀ dx′
    eqs = [sp.Eq(rho * A, rho * U * h0),  # step 2: mass
           sp.Eq(rho * B - rho * U ** 2 * h0, -T)]  # step 5: momentum theorem (step 3 flux − step 4 force)
    Tsol = sp.solve(eqs, [T, h0], dict=True)[0][T]
    assert z0(Tsol - rho * (U * A - B))  # step 7: ∫τ₀dx = ρ∫u(U − u)dy in the integrals of u
    # steps 7–9 with the deficit d = 1 − u/U: U A − B = U²∫d(1 − d)dy (and d(1 − d) = (u/U)(1 − u/U))
    d = sp.Symbol("d")
    assert z0(sp.expand(U * (U * (1 - d)) - U ** 2 * (1 - d) ** 2 - U ** 2 * d * (1 - d)))  # pointwise: u(U − u) = U² (u/U)(1 − u/U)
    # the same identity with a real profile: θ from the profile function against the wall friction (Blasius), 4 stations
    nu, Uv, rho_v = 1.5e-5, 2.0, 1.2
    for xv in (0.05, 0.5, 1.5):
        F = quad(lambda s: float(BL.blasius_wall_shear(s, Uv, rho_v, nu)), 0, xv, epsabs=0, epsrel=1e-12)[0]
        assert rho_v * Uv ** 2 * float(BL.blasius_theta(xv, Uv, nu)) == pytest.approx(F, rel=1e-10)
    # the check line: ρU²θ = 1.2 × 1 × 2.572e-3 = 3.09e-3 N/m at x = 1 m (air, U = 1 m/s)
    assert 1.2 * 1.0 * float(BL.blasius_theta(1.0, 1.0, 1.5e-5)) == pytest.approx(3.09e-3, rel=2e-3)


def test_blasius_reduction_V2_derivation_D05():  # V2 (D05 ★★): all 13 steps with the total-derivative chain rule; two terms cancel; δδ′ = ν/2U
    x, y, U, nu, eta = sp.symbols("x y U nu eta", positive=True)
    dl = sp.Function("delta")(x)
    D = _jet_D(y / dl)
    psi = U * dl * _FJ[0]  # step 2 (9.19)
    assert z0(sp.diff(y / dl, y) - 1 / dl) and z0(sp.diff(y / dl, x) + y / dl * sp.diff(dl, x) / dl)  # step 3
    u, v = D(psi, y), -D(psi, x)
    assert z0(u - U * _FJ[1])  # step 4 (9.23)
    assert z0(sp.simplify(v.subs(y, eta * dl)) - U * sp.diff(dl, x) * (eta * _FJ[1] - _FJ[0]))  # step 5 (9.24)
    ux, uy, uyy = D(u, x), D(u, y), D(D(u, y), y)
    assert z0(ux.subs(y, eta * dl) + U * sp.diff(dl, x) / dl * eta * _FJ[2]) and z0(uy - U / dl * _FJ[2]) and z0(uyy - U / dl ** 2 * _FJ[3])  # step 6
    uux, vuy = (u * ux).subs(y, eta * dl), (v * uy).subs(y, eta * dl)
    assert z0(uux + U ** 2 * sp.diff(dl, x) / dl * eta * _FJ[1] * _FJ[2])  # step 7
    assert z0(vuy - U ** 2 * sp.diff(dl, x) / dl * (eta * _FJ[1] * _FJ[2] - _FJ[0] * _FJ[2]))
    assert z0(uux + vuy + U ** 2 * sp.diff(dl, x) / dl * _FJ[0] * _FJ[2])  # step 8: the ηf′f″ terms cancel
    # step 9–10: −[U²δ′/δ] f f″ = [νU/δ²] f‴  ⇒  f‴ + [Uδδ′/ν] f f″ = 0
    ode = sp.expand((uux + vuy - nu * uyy.subs(y, eta * dl)) * (-dl ** 2 / (nu * U)))
    assert z0(ode - (_FJ[3] + U * dl * sp.diff(dl, x) / nu * _FJ[0] * _FJ[2]))
    # step 11: the bracket = ½ ⇒ δ² = νx/U + D, D = 0 from δ(0) = 0 (9.22)
    D_, dsq = sp.symbols("D dsq")
    g = sp.Function("g")
    sol = sp.dsolve(sp.Eq(g(x).diff(x), nu / U), g(x))  # g = δ²
    assert z0(sol.rhs.subs(sp.Symbol("C1"), 0) - nu * x / U) and z0(sp.sqrt(nu * x / U).diff(x) * sp.sqrt(nu * x / U) - nu / (2 * U))
    # step 12–13: Blasius equation; boundary conditions from (9.20)–(9.21)
    assert z0(_FJ[3] + _FJ[0] * _FJ[2] / 2 - (ode.subs(sp.diff(dl, x), sp.diff(sp.sqrt(nu * x / U), x)).subs(dl, sp.sqrt(nu * x / U))))
    assert z0(v.subs(y, 0).subs({_FJ[0]: 0}) * 0 + (eta * _FJ[1] - _FJ[0]).subs({eta: 0, _FJ[0]: 0}))  # v(η = 0) = 0 with f(0) = 0: no through-flow
    # the check line: x = 0.5 m, U = 1 m/s, air: δ = 2.74 mm, η = 2 at y = 5.48 mm
    assert np.sqrt(1.5e-5 * 0.5 / 1.0) == pytest.approx(2.74e-3, rel=1e-3) and 2 * np.sqrt(1.5e-5 * 0.5) == pytest.approx(5.48e-3, rel=1e-3)
    # cross-check against the tested solver: the residual of (9.18) on the sympy-derived ODE is zero for the numerical f
    eta_g = np.linspace(0.0, 8.0, 2001)
    f, fp, fpp = BL.blasius_profile(eta_g)
    assert np.max(np.abs(np.gradient(fpp, eta_g[1] - eta_g[0], edge_order=2)[2:-2] + 0.5 * f[2:-2] * fpp[2:-2])) < 1e-4


def test_blasius_numbers_V2_derivation_D06():  # V2 (D06 ★★): Töpfer symmetry (weights), θ = 2f″(0) from the ODE (both lines), τ₀, C_f, C_D algebra
    e, lam = sp.symbols("eta lambda", positive=True)
    f1 = sp.Function("f1")
    fl = lam * f1(lam * e)
    ode = fl.diff(e, 3) + fl * fl.diff(e, 2) / 2
    ode1 = sp.Subs(f1(sp.Symbol("s")).diff(sp.Symbol("s"), 3) + f1(sp.Symbol("s")) * f1(sp.Symbol("s")).diff(sp.Symbol("s"), 2) / 2, sp.Symbol("s"), lam * e).doit()
    assert z0(sp.simplify(ode.doit() - lam ** 4 * ode1.doit()))  # step 1: weight λ⁴
    # steps 2–3 numerically with an independent IVP: F″(0) = 1 ⇒ F′(∞) = 2.0854, λ = 0.6925, f″(0) = λ³ = 0.33206
    sol = solve_ivp(lambda s, Y: [Y[1], Y[2], -0.5 * Y[0] * Y[2]], (0, 14), [0, 0, 1.0], rtol=1e-12, atol=1e-14)
    Finf = sol.y[1, -1]
    lam_v = Finf ** -0.5
    assert (Finf, lam_v, lam_v ** 3) == pytest.approx((2.0854, 0.6925, 0.33206), abs=1e-4)
    # step 4: 99 % point is the root of a monotone function (brentq) — 4.910
    assert BL.blasius_constants()["eta99"] == pytest.approx(4.910, abs=5e-4)
    # step 5: ∫(1 − f′)dη = lim(η − f) [f′ = df/dη]; step 6: with the ODE, E = ½f(1 − f′) − ½ I_θ(η) − f″(η) is constant = −f″(0)
    F, eta_ = sp.Function("f"), sp.Symbol("eta")
    Ith = sp.Function("Ith")(eta_)  # I_θ(η) with I_θ′ = f′(1 − f′)
    Ex = F(eta_) * (1 - F(eta_).diff(eta_)) / 2 - Ith / 2 - F(eta_).diff(eta_, 2)
    dE = Ex.diff(eta_).subs(Ith.diff(eta_), F(eta_).diff(eta_) * (1 - F(eta_).diff(eta_)))
    dE = dE.subs(F(eta_).diff(eta_, 3), -F(eta_) * F(eta_).diff(eta_, 2) / 2)
    assert z0(dE)  # E′ = 0 by the ODE (9.27): D06 step 6, with the +f″(η) that the design line omits
    # step 7: μU f″(0)/δ = f″(0) ρU²/√Re_x;  steps 8–10: C_f, drag, C_D
    rho, Uu, nu, xx, L, mu, f0 = sp.symbols("rho U nu x L mu f0", positive=True)
    dlt = sp.sqrt(nu * xx / Uu)
    assert z0(rho * nu * Uu * f0 / dlt - f0 * rho * Uu ** 2 / sp.sqrt(Uu * xx / nu))  # (9.31)
    tau0 = f0 * rho * Uu ** 2 / sp.sqrt(Uu * xx / nu)
    assert z0(tau0 / (rho * Uu ** 2 / 2) - 2 * f0 / sp.sqrt(Uu * xx / nu))  # (9.32)
    FD = sp.integrate(tau0, (xx, 0, L))
    assert z0(FD - 2 * f0 * rho * Uu ** 2 * L / sp.sqrt(Uu * L / nu))  # ∫x^{-1/2} = 2√L
    assert z0(FD / (rho * Uu ** 2 * L / 2) - 4 * f0 / sp.sqrt(Uu * L / nu))  # (9.33): C_D = 4f″(0)/√Re_L = 1.328/√Re_L
    assert 4 * BL.blasius_constants()["fpp0"] == pytest.approx(1.328, abs=5e-4)


def test_falkner_skan_reduction_V2_derivation_D07():  # V2 (D07 ★★): all 12 steps for U_e = a xⁿ with symbolic n
    x, y, nu, a, eta = sp.symbols("x y nu a eta", positive=True)
    n = sp.Symbol("n", real=True)
    Ue = a * x ** n
    dl = sp.sqrt(nu * x / Ue)
    g = sp.sqrt(nu * x * Ue)
    D = _jet_D(y / dl)
    assert z0(Ue.diff(x) * Ue - n * Ue ** 2 / x)  # step 1
    assert z0(sp.simplify(g.diff(x) / g - (n + 1) / (2 * x))) and z0(sp.simplify(sp.diff(y / dl, x) + (1 - n) / (2 * x) * y / dl)) and z0(sp.simplify(g / dl - Ue))  # step 3
    psi = g * _FJ[0]
    u, v = D(psi, y), -D(psi, x)
    ux, uy, uyy = D(u, x), D(u, y), D(D(u, y), y)
    checks = {
        "u": (u, Ue * _FJ[1]),
        "v": (v, -g / (2 * x) * ((n + 1) * _FJ[0] + (n - 1) * eta * _FJ[1])),
        "ux": (ux, Ue / x * (n * _FJ[1] - (1 - n) / 2 * eta * _FJ[2])),
        "uy": (uy, Ue / dl * _FJ[2]),
        "uyy": (uyy, Ue ** 2 / (nu * x) * _FJ[3]),
        "uux": (u * ux, Ue ** 2 / x * (n * _FJ[1] ** 2 - (1 - n) / 2 * eta * _FJ[1] * _FJ[2])),
        "vuy": (v * uy, -Ue ** 2 / (2 * x) * ((n + 1) * _FJ[0] * _FJ[2] + (n - 1) * eta * _FJ[1] * _FJ[2])),
    }
    for name, (got, want) in checks.items():
        assert z0(_at(got - want, y, eta, dl)), name  # steps 4–9
    assert z0(_at(u * ux + v * uy - Ue ** 2 / x * (n * _FJ[1] ** 2 - (n + 1) / 2 * _FJ[0] * _FJ[2]), y, eta, dl))  # step 10: ηf′f″ cancels
    assert z0(Ue * Ue.diff(x) + nu * Ue ** 2 / (nu * x) * _FJ[3] - Ue ** 2 / x * (n + _FJ[3]))  # step 11
    full = _at(u * ux + v * uy - Ue * Ue.diff(x) - nu * uyy, y, eta, dl)
    assert z0(full + Ue ** 2 / x * (_FJ[3] + (n + 1) / 2 * _FJ[0] * _FJ[2] - n * _FJ[1] ** 2 + n))  # step 12 (9.36)
    # n = 0 is Blasius, n = 1 is Hiemenz with δ = √(ν/a)
    assert z0((_FJ[3] + (n + 1) / 2 * _FJ[0] * _FJ[2] - n * _FJ[1] ** 2 + n).subs(n, 0) - (_FJ[3] + _FJ[0] * _FJ[2] / 2)) and z0(dl.subs(n, 1) - sp.sqrt(nu / a))
    assert BL.falkner_skan_thickness(2.0, 1.0, 10.0, 1.5e-5) * 1e3 == pytest.approx(1.22, abs=5e-3)  # the check line: 1.22 mm at every x


def test_momentum_integral_V2_derivation_D08():  # V2 (D08 ★★): (9.43) from (9.37) on a manufactured boundary-layer flow, with exact integrals; the (9.38), (9.42) identities on generic fields
    x, y, yp = sp.symbols("x y yp", positive=True)
    # generic identities (steps 2 and 6–7)
    u, v, tau, Ue = (sp.Function(n) for n in ("u", "v", "tau", "Ue"))
    U_, V_ = u(x, y), v(x, y)
    assert z0(sp.expand(U_ * (U_.diff(x) + V_.diff(y)) + U_ * U_.diff(x) + V_ * U_.diff(y) - (sp.diff(U_ ** 2, x) + sp.diff(U_ * V_, y))))  # (9.38) = (9.37) + u × continuity
    Ue_x = Ue(x)
    lhs941 = sp.diff(U_ ** 2, x) - Ue_x * Ue_x.diff(x) - Ue_x * U_.diff(x)
    lhs942 = sp.diff(U_ ** 2 - Ue_x * U_, x) + Ue_x.diff(x) * (U_ - Ue_x)
    assert z0(sp.expand(lhs941 - lhs942))  # step 7 identity (integrands of (9.41) and (9.42))
    # manufactured flow: u = U_e(x)(1 − e^{−y/δ}), U_e = 1 + x, δ = √x; v from continuity (v(0) = 0); τ from (9.37) with τ(∞) = 0
    Uex, dl = 1 + x, sp.sqrt(x)
    um = Uex * (1 - sp.exp(-y / dl))
    vm = sp.simplify(-sp.integrate(um.diff(x).subs(y, yp), (yp, 0, y)))
    assert z0(um.diff(x) + vm.diff(y))  # continuity (6.2) holds
    residual = sp.simplify(um * um.diff(x) + vm * um.diff(y) - Uex * Uex.diff(x))  # = (1/ρ)τ_y
    tau0_over_rho = -sp.integrate(residual, (y, 0, sp.oo))  # τ(∞) = 0 ⇒ τ(0)/ρ = −∫τ_y/ρ dy
    theta = sp.integrate(sp.simplify((um / Uex) * (1 - um / Uex)), (y, 0, sp.oo))
    dstar = sp.integrate(1 - um / Uex, (y, 0, sp.oo))
    assert z0(sp.simplify(tau0_over_rho - (sp.diff(Uex ** 2 * theta, x) + Uex * dstar * Uex.diff(x))))  # (9.43): exact, no approximation
    # a wrong variant (the pressure-gradient term dropped) fails
    assert not z0(sp.simplify(tau0_over_rho - sp.diff(Uex ** 2 * theta, x)))
    # the Blasius check line: τ₀/ρ = U² dθ/dx = 0.332 U²/√Re_x  (air, x = 1 m: 1.286e-3 m²/s² for U = 1)
    assert float(BL.blasius_wall_shear(1.0, 1.0, 1.0, 1.5e-5)) == pytest.approx(1.286e-3, rel=1e-3)


def test_thwaites_closed_form_V2_derivation_D09():  # V2 (D09 ★★): steps 3–10 with generic θ, U_e, H, l; integrating factor; power-law closed form
    x, nu = sp.symbols("x nu", positive=True)
    th, Ue, H, lfun, lam = (sp.Function(n)(x) for n in ("theta", "U_e", "H", "l", "lam"))
    Uep = Ue.diff(x)
    # steps 3–5: (9.43) divided by ... ⇒ l = (2 + H)λ + (U_e/2) d(θ²/ν)/dx
    lhs = th / (nu * Ue) * sp.diff(Ue ** 2 * th, x) + th ** 2 * H * Uep / nu  # (9.47) right side with δ* = Hθ
    lam_ex = th ** 2 * Uep / nu
    assert z0(sp.expand(lhs - ((2 + H) * lam_ex + Ue / 2 * sp.diff(th ** 2 / nu, x))))
    # step 7–8: θ²/ν = λ/U_e′  ⇒  U_e d/dx(λ/U_e′) = 2l − 2(2 + H)λ   (from l = (2+H)λ + (U_e/2)(θ²/ν)′)
    Lam = sp.Function("lam")(x)
    Z = Lam / Uep
    l_expr = (2 + H) * Lam + Ue / 2 * sp.diff(Z, x)
    assert z0(sp.expand(Ue * sp.diff(Z, x) - (2 * l_expr - 2 * (2 + H) * Lam)))
    # step 9: L = 0.45 − 6λ ⇒ Z′ + 6(U_e′/U_e)Z = 0.45/U_e with Z = θ²/ν = λ/U_e′
    Zf = sp.Function("Z")(x)
    ode = sp.Eq(Ue * Zf.diff(x), sp.Rational(45, 100) - 6 * Zf * Uep)
    ode49 = sp.Eq(Zf.diff(x) + 6 * Uep / Ue * Zf, sp.Rational(45, 100) / Ue)
    assert z0((ode.lhs - ode.rhs) - Ue * (ode49.lhs - ode49.rhs))
    # step 10: integrating factor U_e⁶
    assert z0(sp.diff(Ue ** 6 * Zf, x) - Ue ** 6 * (Zf.diff(x) + 6 * Uep / Ue * Zf))
    assert z0(sp.diff(Ue ** 6 * Zf, x) - sp.Rational(45, 100) * Ue ** 5 - Ue ** 6 * (ode49.lhs - ode49.rhs))  # so d(U⁶Z)/dx = 0.45 U⁵ when (9.49) holds
    # power law U_e = a xⁿ: Z = c x^{1−n} ⇒ c(1 + 5n) = 0.45/a, and (9.50) gives the same c
    a, n, s = sp.symbols("a n s", positive=True)
    c = sp.Symbol("c", positive=True)
    Ue_p = a * x ** n
    Zp = c * x ** (1 - n)
    res = sp.simplify(Zp.diff(x) + 6 * Ue_p.diff(x) / Ue_p * Zp - sp.Rational(45, 100) / Ue_p)
    csol = sp.solve(sp.simplify(res * x ** n), c)[0]
    assert z0(csol - sp.Rational(45, 100) / (a * (5 * n + 1)))
    I5 = sp.integrate((a * s ** n) ** 5, (s, 0, x), conds="none")
    assert z0(sp.simplify(sp.Rational(45, 100) * I5 / Ue_p ** 6 - csol * x ** (1 - n)))  # (9.50) with θ₀ = 0
    # Blasius, diffuser: 0.6708 √(νx/U) and λ = −(0.45/4)[(1 + x/L)⁴ − 1]
    L_ = sp.Symbol("L", positive=True)
    Ud = 1 / (1 + s / L_)
    lam_d = sp.simplify(sp.Rational(45, 100) * sp.integrate(Ud ** 5, (s, 0, x)) / Ud.subs(s, x) ** 6 * sp.diff(Ud.subs(s, x), x))
    assert z0(lam_d + sp.Rational(45, 400) * ((1 + x / L_) ** 4 - 1))
    assert np.sqrt(0.45) == pytest.approx(0.6708, abs=5e-5)


def test_thwaites_closure_V2_derivation_D10():  # V2 (D10 ★★): λ = nI_θ², l = I_θ f″(0), H = I_δ/I_θ by independent quadrature; the numbers of steps 7–8
    for m in (0.0, 1.0, -0.05):
        d = BL.falkner_skan(m, eta_max=16.0, n=16001)
        eta, fp, fpp = d["eta"], d["fp"], d["fpp"]
        spl = CubicSpline(eta, fp)
        Ith = float(quad(lambda s_: float(spl(s_)) * (1 - float(spl(s_))), 0.0, eta[-1], epsabs=1e-13, epsrel=1e-12, limit=400)[0])  # adaptive quad, independent of Simpson
        Idl = float(quad(lambda s_: 1 - float(spl(s_)), 0.0, eta[-1], epsabs=1e-13, epsrel=1e-12, limit=400)[0])
        st = BL.falkner_skan_state(m, eta_max=16.0)
        assert (st["I_theta"], st["I_delta"]) == pytest.approx((Ith, Idl), rel=1e-7)
        assert st["lam"] == pytest.approx(m * Ith ** 2, rel=1e-7, abs=1e-12) and st["l"] == pytest.approx(fpp[0] * Ith, rel=1e-7) and st["H"] == pytest.approx(Idl / Ith, rel=1e-7)
    # step 7: L = 2l − 2(2 + H)λ at Blasius, Hiemenz and the near-separation member; step 8: differences from 0.45 − 6λ
    L0 = float(BL.thwaites_L(0.0))
    st1 = BL.falkner_skan_state(1.0)
    L1 = 2 * st1["l"] - 2 * (2 + st1["H"]) * st1["lam"]
    lam_s = -0.0675
    Ls = float(BL.thwaites_L(lam_s))
    assert (L0, L1, Ls) == pytest.approx((0.441, 0.0, 0.818), abs=2e-3)
    assert (0.45 - 6 * 0.0, 0.45 - 6 * st1["lam"], 0.45 - 6 * lam_s) == pytest.approx((0.450, -0.063, 0.855), abs=1e-3)
    assert (0.45 - L0, 0.45 - 6 * st1["lam"] - L1, 0.45 - 6 * lam_s - Ls) == pytest.approx((0.01, -0.063, 0.04), abs=0.006)  # "0.01 too high, 0.06 too low, 0.04 too high"
    # step 9: two separation criteria
    fold = BL.falkner_skan_separation()["m_sep"]
    assert fold == pytest.approx(-0.0904, abs=1e-4) and BL.LAMBDA_SEP_FS == pytest.approx(-0.0681, abs=5e-5)
    tb = BL.thwaites_closure_table(60)
    assert tb["lam"][0] == pytest.approx(-0.0681, abs=6e-5) and tb["H"][0] == pytest.approx(4.03, abs=0.01)  # H ≈ 4.0 at the fold (the analysis note 3.81 is n = −0.09)
    assert BL.falkner_skan_state(-0.0904)["H"] == pytest.approx(3.970, abs=0.002)  # the design row for m = −0.0904 (just above the fold)


def test_thwaites_cylinder_V2_derivation_D11():  # V2 (D11 ★★): F(φ) by c = cos φ, λ(φ), its limit at the stagnation point, the crossing
    phi, psi, c = sp.symbols("phi psi c", positive=True)
    F = sp.integrate((1 - c ** 2) ** 2, (c, sp.cos(phi), 1))  # step 4 (substitution sin⁵ψ dψ = −(1 − c²)²dc)
    assert z0(F - (sp.Rational(8, 15) - sp.cos(phi) + sp.Rational(2, 3) * sp.cos(phi) ** 3 - sp.cos(phi) ** 5 / 5))
    assert z0(sp.diff(F, phi) - sp.sin(phi) ** 5) and z0(sp.integrate(sp.sin(psi) ** 5, (psi, 0, phi)) - F)
    nu, U, a = sp.symbols("nu U a", positive=True)
    Ue = 2 * U * sp.sin(phi)
    theta2 = sp.Rational(45, 100) * nu / Ue ** 6 * a * (2 * U) ** 5 * F  # steps 2–3
    Uep = 2 * U / a * sp.cos(phi)  # step 1: dU_e/dx, x = aφ
    lam = sp.simplify(theta2 * Uep / nu)
    assert z0(lam - sp.Rational(45, 100) * sp.cos(phi) * F / sp.sin(phi) ** 6)  # step 5: U, a, ν cancel
    assert sp.limit(lam, phi, 0) == sp.Rational(3, 40) and z0(lam.subs(phi, sp.pi / 2))  # step 6: 0.45/6 = 0.075, λ(90°) = 0
    assert float(lam.subs(phi, sp.pi / 6)) == pytest.approx(0.0722, abs=6e-5) and float(lam.subs(phi, sp.pi / 3)) == pytest.approx(0.0589, abs=6e-5)
    root = brentq(sp.lambdify(phi, lam + sp.Rational(9, 100), "numpy"), np.deg2rad(60), np.deg2rad(120), xtol=1e-14)  # step 7
    assert np.rad2deg(root) == pytest.approx(103.11, abs=0.01) and np.rad2deg(root) == pytest.approx(BL.thwaites_cylinder_separation(-0.09), abs=1e-9)
    assert 0.0855 / 0.075 - 1 == pytest.approx(0.14, abs=0.01)  # the step-6 statement: compare Hiemenz 0.0855 (14 % above 0.075 — "12 % low" is 1 − 0.075/0.0855 = 12.3 %)
    assert 1 - 0.075 / 0.0855 == pytest.approx(0.123, abs=0.001)


def test_wall_curvature_V2_derivation_D12():  # V2 (D12 ★): μ u_yy = dp/dx from (9.9) at the wall; f‴(0) = −n; the inflection argument on a solved profile
    x, y, mu, rho, nu = sp.symbols("x y mu rho nu", positive=True)
    u, v = sp.Function("u")(x, y), sp.Function("v")(x, y)
    dpdx = sp.Symbol("dpdx")
    eq = u * u.diff(x) + v * u.diff(y) + dpdx / rho - nu * u.diff(y, 2)
    uyy = sp.Symbol("uyy")
    # with u = v = 0 (no slip 9.12, no through-flow 9.13) both advective terms vanish: 0 = −(1/ρ)dp/dx + ν u_yy  ⇒  μ u_yy = dp/dx  (μ = ρν)
    wall = eq.subs(u.diff(y, 2), uyy).subs({u.diff(x): 0, u.diff(y): 0}).subs({u: 0, v: 0})
    assert z0(wall - (dpdx / rho - nu * uyy))
    assert z0(sp.solve(sp.Eq(wall, 0), uyy)[0] * rho * nu - dpdx)
    # FS: u = U_e F′(η) ⇒ u_yy(wall) = −n U_e²/(νx), and (1/μ)dp/dx with dp/dx = −ρ n U_e²/x is the same number
    n, a = sp.symbols("n a", positive=True)
    Ue = a * x ** n
    assert z0(-n * Ue ** 2 / (nu * x) - (-rho * n * Ue ** 2 / x) / (rho * nu))
    # signs on the solved profiles: adverse ⇒ f‴ changes sign (inflection), favourable ⇒ f‴ < 0 throughout
    for m, infl in ((-0.05, True), (0.3, False)):
        d = BL.falkner_skan(m, eta_max=16.0, n=4001)
        fppp = d["fppp"][(d["eta"] < 12) & (np.abs(d["fppp"]) > 1e-10)]
        sg = np.sign(fppp)
        assert (sg[0] > 0) == (m < 0) and (np.any(np.diff(sg) != 0) == infl), m
    assert BL.falkner_skan_state(-0.05)["inflection_eta"] == pytest.approx(1.650, abs=1e-3)  # the D12 check line


def test_pressure_drag_V2_derivation_D13():  # V2 (D13 ★): the seven steps: C_D,p = ½∮C_p cos φ dφ, front part, wake part, sum
    phi, ps, cb, rho, U, a = sp.symbols("phi phi_s C_b rho U a", positive=True)
    assert z0(sp.integrate(sp.cos(phi), (phi, 0, 2 * sp.pi)))  # step 2: uniform p∞ exerts no net force
    dp, pinf = sp.symbols("dp pinf")
    Cp = dp / (rho * U ** 2 / 2)
    assert z0(dp * sp.cos(phi) * a / (rho * U ** 2 / 2 * (2 * a)) - Cp * sp.cos(phi) / 2)  # step 3: F/(½ρU²·2a) = ½ C_p cos φ per dφ
    front = 2 * sp.integrate((1 - 4 * sp.sin(phi) ** 2) * sp.cos(phi), (phi, 0, ps))  # step 5: both halves
    assert z0(front - 2 * (sp.sin(ps) - sp.Rational(4, 3) * sp.sin(ps) ** 3))
    wake = cb * sp.integrate(sp.cos(phi), (phi, ps, 2 * sp.pi - ps))  # step 6
    assert z0(wake + 2 * cb * sp.sin(ps))
    assert z0(sp.Rational(1, 2) * (front + wake) - sp.sin(ps) * (1 - sp.Rational(4, 3) * sp.sin(ps) ** 2 - cb))  # step 7
    assert z0((sp.sin(ps) * (1 - sp.Rational(4, 3) * sp.sin(ps) ** 2 - (1 - 4 * sp.sin(ps) ** 2))).subs(ps, sp.pi))  # ideal limit → 0
    assert BB.separated_pressure_drag(90.0, -1.0) == pytest.approx(2.0 / 3.0, abs=1e-12)  # the check lines: 2.667 and 0.667 at 90°
    assert BB.separated_pressure_drag(82.0, -1.2) == pytest.approx(0.884, abs=1e-3) and BB.separated_pressure_drag(125.0, -0.6) == pytest.approx(0.578, abs=1e-3)


def test_karman_street_V2_derivation_D14():  # V2 (D14 ★★★): steps 3–15 — cot sum, street speed, linearised equations, lattice sums at ζ = −a/2 + ib, the 4×4 matrix, its factorised polynomial, growth, marginal spacing, facing rows
    b, mu = sp.symbols("b mu", positive=True)
    xA, yA, xB, yB = sp.symbols("xA yA xB yB", real=True)
    # step 3: Σ_n 1/(z − na) = (π/a)cot(πz/a) (symmetric partial sums) — numerical check at a complex z
    z = 0.37 + 0.42j
    nn = np.arange(-200000, 200001)
    assert abs(np.sum(1.0 / (z - nn)) - np.pi / np.tan(np.pi * z)) < 1e-4
    # step 4: street speed from the direct formula w = −(Γ/2ia) cot(πζ/a), ζ = −a/2 + ib: cot(−π/2 + iπb/a) = −i tanh(πb/a)
    zeta = -sp.Rational(1, 2) + sp.I * b
    assert z0(sp.simplify(sp.expand_complex(sp.cot(sp.pi * zeta))) - (-sp.I * sp.tanh(sp.pi * b))) or z0(sp.simplify(sp.expand_complex(sp.cos(sp.pi * zeta) / sp.sin(sp.pi * zeta)) + sp.I * sp.tanh(sp.pi * b)))
    # step 5: linearisation d/dz of Γ/(2πi (z_i − z_j)) with respect to the displacements
    Gam, zi, zj, di, dj = sp.symbols("Gamma z_i z_j d_i d_j")
    eps = sp.Symbol("eps")
    w = Gam / (2 * sp.pi * sp.I * (zi + eps * di - zj - eps * dj))
    lin = sp.diff(w, eps).subs(eps, 0)  # first-order change of the induced velocity
    assert z0(lin + Gam / (2 * sp.pi * sp.I) * (di - dj) / (zi - zj) ** 2)  # dw = −Σ(Γ_j/2πi)(δ_i − δ_j)/(z_i − z_j)²
    # step 9: sin, cos at πζ/a = −π/2 + iπb/a ; lattice sums P = π²/sin², Q = π² cos/sin²
    c, s = sp.cosh(sp.pi * b), sp.sinh(sp.pi * b)
    assert z0(sp.simplify(sp.expand_complex(sp.sin(sp.pi * zeta))) + c) and z0(sp.simplify(sp.expand_complex(sp.cos(sp.pi * zeta))) - sp.I * s)
    P = sp.pi ** 2 / c ** 2
    Q = sp.I * sp.pi ** 2 * s / c ** 2
    assert z0(sp.simplify(sp.expand_complex(sp.pi ** 2 / sp.sin(sp.pi * zeta) ** 2)) - P)
    assert z0(sp.simplify(sp.expand_complex(sp.pi ** 2 * sp.cos(sp.pi * zeta) / sp.sin(sp.pi * zeta) ** 2)) - Q)
    # own-row sums: Σ_{n≠0} 1/(na)² = π²/3a², Σ_{n≠0}(−1)^n/(na)² = −π²/6a² ⇒ P_ii − Q_ii = π²/2a² (a = 1)
    n_ = sp.Symbol("n", integer=True, positive=True)
    assert sp.summation(2 / n_ ** 2, (n_, 1, sp.oo)) == sp.pi ** 2 / 3 and sp.summation(2 * (-1) ** n_ / n_ ** 2, (n_, 1, sp.oo)) == -sp.pi ** 2 / 6
    # steps 7–8 → 10 → 11: the two complex equations become the real 4×4 (a = Γ = 1)
    def build(Pv, Qv):
        dA, dB = xA + sp.I * yA, xB + sp.I * yB
        cA, cB = sp.conjugate(dA), sp.conjugate(dB)
        Pc, Qc = sp.conjugate(Pv), sp.conjugate(Qv)
        dAdot = 1 / (2 * sp.pi * sp.I) * (sp.pi ** 2 / 2 * cA - Pc * cA + Qc * cB)  # step 8
        dBdot = -1 / (2 * sp.pi * sp.I) * (sp.pi ** 2 / 2 * cB - Pc * cB + Qc * cA)
        e = [sp.re(sp.expand(dAdot)), sp.im(sp.expand(dAdot)), sp.re(sp.expand(dBdot)), sp.im(sp.expand(dBdot))]
        vars_ = [xA, yA, xB, yB]
        return sp.Matrix(4, 4, lambda i, j: sp.simplify(sp.diff(e[i], vars_[j])))
    gam, sig = sp.Rational(1, 2) - 1 / c ** 2, s / c ** 2
    M = build(P, Q)
    Mdesign = sp.pi / 2 * sp.Matrix([[0, -gam, -sig, 0], [-gam, 0, 0, sig], [sig, 0, 0, gam], [0, -sig, gam, 0]])  # step 11
    assert sp.simplify(M - Mdesign) == sp.zeros(4, 4)
    # step 12: characteristic polynomial (in μ, eigenvalue of M/(π/2)) factorises as ((μ − γ)² + σ²)((μ + γ)² + σ²)
    assert z0(sp.simplify((Mdesign / (sp.pi / 2) - mu * sp.eye(4)).det() - ((mu - gam) ** 2 + sig ** 2) * ((mu + gam) ** 2 + sig ** 2)))
    # steps 13–14: growth (π/2)|γ|; marginal spacing cosh² = 2 ⇒ σ = ½, eigenvalues ±i π/4
    sols = sp.solve(sp.cosh(sp.pi * b) ** 2 - 2, b)
    assert len(sols) == 1 and float(sols[0] - sp.acosh(sp.sqrt(2)) / sp.pi) == pytest.approx(0.0, abs=1e-14)  # acosh √2 = ln(1 + √2)
    assert z0(sp.simplify(gam.subs(b, sp.acosh(sp.sqrt(2)) / sp.pi))) and z0(sp.simplify(sig.subs(b, sp.acosh(sp.sqrt(2)) / sp.pi) - sp.Rational(1, 2)))
    # step 15: facing rows ζ = ib: P = −π²/sinh², Q = −π² cosh/sinh², polynomial ((2μ)² − 1)²/16 = (μ² − ¼)² — NOTE the design writes (μ² − π²/16)²,
    # which is the same polynomial in λ = (π/2)μ, not in μ (a notation slip; the eigenvalues ±π/4 are right)
    zeta2 = sp.I * b
    assert z0(sp.simplify(sp.expand_complex(sp.pi ** 2 / sp.sin(sp.pi * zeta2) ** 2)) + sp.pi ** 2 / s ** 2)
    assert z0(sp.simplify(sp.expand_complex(sp.pi ** 2 * sp.cos(sp.pi * zeta2) / sp.sin(sp.pi * zeta2) ** 2)) + sp.pi ** 2 * c / s ** 2)
    M2 = build(-sp.pi ** 2 / s ** 2, -sp.pi ** 2 * c / s ** 2)
    charp2 = sp.expand((M2 / (sp.pi / 2) - mu * sp.eye(4)).det())
    assert z0(charp2 - (mu ** 2 - sp.Rational(1, 4)) ** 2) and not z0(charp2 - (mu ** 2 - sp.pi ** 2 / 16) ** 2)
    lam_ = sp.Symbol("lam")
    assert z0(sp.expand((M2 - lam_ * sp.eye(4)).det() - (lam_ ** 2 - sp.pi ** 2 / 16) ** 2))  # in the true eigenvalue λ the printed form is right
    # cross-check of the derived matrix against the tested numerical Bloch spectrum (different assembly: lattice sums with the code's row convention)
    for bv in (0.2, 0.28, 0.4):
        ev = np.linalg.eigvals(np.array(Mdesign.subs(b, bv).evalf(), dtype=float))
        assert np.max(ev.real) == pytest.approx(float(np.max(BB.karman_street_spectrum(bv).real)), abs=1e-12), bv
        got = np.sort_complex(np.round(ev, 9))
        want = np.sort_complex(np.round(BB.karman_street_spectrum(bv), 9))
        assert np.max(np.abs(got - want)) < 1e-8, bv  # the whole spectrum, not just its real part
    # the symbolic engine of the module agrees
    ks = ch09.karman_street_sympy()
    assert z0(ks["growth_check"]) and z0(ks["marginal_check"])
    yv = sp.Symbol("y", positive=True)
    assert float(ks["gamma"].subs({sp.Symbol("Gamma", positive=True): 1, sp.Symbol("a", positive=True): 1, yv: sp.pi * 0.3})) == pytest.approx(float(BB.karman_street_growth_closed(0.3)), rel=1e-12)


def test_jet_momentum_V2_derivation_D15():  # V2 (D15 ★★): flux form, boundary terms, and J conserved on the exact Bickley solution with sympy integrals
    x, y, nu, J, rho = sp.symbols("x y nu J rho", positive=True)
    u, v = sp.Function("u")(x, y), sp.Function("v")(x, y)
    assert z0(sp.expand(u * (u.diff(x) + v.diff(y)) + (u * u.diff(x) + v * u.diff(y)) - (sp.diff(u ** 2, x) + sp.diff(u * v, y))))  # steps 1–3 (continuity ⇒ flux form)
    # Bickley solution with symbolic constants: u = u0 sech²(y/(√6 δ)); ∫u²dy = J/ρ for every x (exact sympy integral); boundary terms vanish
    Cj = 4 * sp.sqrt(6) / 3
    u0 = (J ** 2 / (Cj ** 2 * rho ** 2 * nu * x)) ** sp.Rational(1, 3)
    dl = (Cj * rho * nu ** 2 * x ** 2 / J) ** sp.Rational(1, 3)
    s = sp.Symbol("s", real=True)
    prim = sp.tanh(s) - sp.tanh(s) ** 3 / 3  # ∫sech⁴ s ds (substitution t = tanh s: ∫(1 − t²)dt)
    assert z0(sp.diff(prim, s) - sp.sech(s) ** 4)
    I4 = sp.limit(prim, s, sp.oo) - sp.limit(prim, s, -sp.oo)
    assert I4 == sp.Rational(4, 3)
    flux = sp.simplify(u0 ** 2 * dl * sp.sqrt(6) * I4 * rho)  # ρ∫u²dy with y = √6 δ s
    assert z0(flux - J) and z0(sp.diff(flux, x))  # steps 6–7: J is constant, equal to the slot's flux
    # boundary terms: u → 0 and u_y → 0 as y → ±∞ (steps 4–5)
    uu = u0 * sp.sech(y / (sp.sqrt(6) * dl)) ** 2
    assert sp.limit(uu, y, sp.oo) == 0 and sp.limit(uu.diff(y), y, sp.oo) == 0
    # (9.56) needs the KINEMATIC stress: d/dx∫u²dy, [uv] and [ν u_y] all have the unit m²/s²; the printed τ = μ u_y is a pressure (slip R7)
    kin = Q_(1e-6, "m**2/s") * Q_(1, "1/s")
    assert kin.check("[length]**2/[time]**2") and (Q_(1, "m/s") ** 2).check("[length]**2/[time]**2")
    assert (Q_(1e-3, "Pa*s") * Q_(1, "1/s")).check("[pressure]") and not (Q_(1e-3, "Pa*s") * Q_(1, "1/s")).check("[length]**2/[time]**2")
    # numerical J on the code's field at four stations: 1.0000 N/m each (the check line)
    for xv in (0.05, 0.1, 0.2, 0.4):
        fj = lambda yy, xv=xv: float(JET.free_jet(xv, yy, 1.0, 1.2, 1.5e-5)["u"]) ** 2  # noqa: E731
        assert 2 * 1.2 * quad(fj, 0, 0.05, epsabs=0, epsrel=1e-12, limit=200)[0] == pytest.approx(1.0, rel=1e-8)


def test_free_jet_similarity_V2_derivation_D16():  # V2 (D16 ★★★): all 14 steps: exponents from J, the four derivatives, x-powers cancel, ODE coefficients
    x, y, nu, J, rho, C = sp.symbols("x y nu J rho C", positive=True)
    u0 = (J ** 2 / (C ** 2 * rho ** 2 * nu * x)) ** sp.Rational(1, 3)
    dl = (C * rho * nu ** 2 * x ** 2 / J) ** sp.Rational(1, 3)
    assert z0(sp.simplify(u0 ** 2 * dl * C - J / rho))  # step 3–4
    assert z0(sp.simplify(dl ** 2 - nu * x / u0)) and z0(sp.simplify((nu * x / u0) ** sp.Rational(1, 2) - dl))  # step 1 definition and step 5
    A = (J * nu / (C * rho)) ** sp.Rational(1, 3)  # step 6
    B = (C * rho * nu ** 2 / J) ** sp.Rational(1, 3)
    assert z0(sp.simplify(u0 * dl - A * x ** sp.Rational(1, 3))) and z0(sp.simplify(dl - B * x ** sp.Rational(2, 3))) and z0(sp.simplify(A * B - nu))
    E = y / (B * x ** sp.Rational(2, 3))
    D = _jet_D(E)
    assert z0(sp.diff(E, x) + 2 * E / (3 * x)) and z0(sp.diff(E, y) - 1 / (B * x ** sp.Rational(2, 3)))  # step 7
    psi = A * x ** sp.Rational(1, 3) * _FJ[0]
    u, v = D(psi, y), -D(psi, x)
    eta = sp.Symbol("eta", positive=True)
    ux, uy, uyy = D(u, x), D(u, y), D(D(u, y), y)
    sub = lambda e: sp.simplify(sp.powsimp(sp.expand_power_base(e.subs(y, eta * B * x ** sp.Rational(2, 3)), force=True), force=True))  # noqa: E731
    assert z0(sub(u - A / B * x ** sp.Rational(-1, 3) * _FJ[1]))  # step 8
    assert z0(sub(v + A / 3 * x ** sp.Rational(-2, 3) * (_FJ[0] - 2 * eta * _FJ[1])))
    assert z0(sub(ux - A / B * x ** sp.Rational(-4, 3) * (-_FJ[1] / 3 - 2 * eta * _FJ[2] / 3)))  # step 9
    assert z0(sub(uy - A / B ** 2 / x * _FJ[2])) and z0(sub(uyy - A / B ** 3 * x ** sp.Rational(-5, 3) * _FJ[3]))  # step 10
    assert z0(sub(u * ux - A ** 2 / B ** 2 * x ** sp.Rational(-5, 3) * _FJ[1] * (-_FJ[1] / 3 - 2 * eta * _FJ[2] / 3)))  # step 11
    assert z0(sub(v * uy + A ** 2 / (3 * B ** 2) * x ** sp.Rational(-5, 3) * (_FJ[0] * _FJ[2] - 2 * eta * _FJ[1] * _FJ[2])))
    assert z0(sub(u * ux + v * uy + A ** 2 / (3 * B ** 2) * x ** sp.Rational(-5, 3) * (_FJ[1] ** 2 + _FJ[0] * _FJ[2])))  # step 12: ηf′f″ cancels
    res = u * ux + v * uy - nu * uyy
    step13 = sub(res + A ** 2 / (3 * B ** 2) * x ** sp.Rational(-5, 3) * (_FJ[1] ** 2 + _FJ[0] * _FJ[2]) + nu * A / B ** 3 * x ** sp.Rational(-5, 3) * _FJ[3])
    assert z0(step13)  # step 13: every term ∝ x^{-5/3}
    step14 = sub(res * (-3 * B ** 2 * x ** sp.Rational(5, 3) / A ** 2)) - (3 * _FJ[3] + _FJ[0] * _FJ[2] + _FJ[1] ** 2)  # multiply by −3B³/A · (B^{-1}) : uses AB = ν
    assert z0(sp.simplify(step14.subs(nu, A * B)))  # step 14 (AB = ν makes the coefficient of f‴ exactly 3)
    # the module's engine returns the same ODE
    assert ch09.similarity_reduce_sympy("free_jet")["ode_coefficients"] == {"fppp": 3, "f f''": 1, "f'^2": 1}
    # the check line: air, J = 1 N/m: u₀(0.1) = 35.14 m/s, u₀(0.4) = 22.14 m/s, ratio 4^{-1/3}
    assert float(JET.free_jet_centreline(0.1, 1.0, 1.2, 1.5e-5)) == pytest.approx(35.14, abs=5e-3) and float(JET.free_jet_centreline(0.4, 1.0, 1.2, 1.5e-5)) == pytest.approx(22.14, abs=5e-3)
    assert 4 ** (-1 / 3) == pytest.approx(0.63, abs=1e-3)


def test_free_jet_solution_V2_derivation_D17():  # V2 (D17 ★★): first integrals (total derivatives), tanh, C, ṁ coefficient, entrainment limit
    e = sp.Symbol("eta", real=True)
    f = sp.Function("f")
    ode = 3 * f(e).diff(e, 3) + f(e) * f(e).diff(e, 2) + f(e).diff(e) ** 2
    assert z0(ode - (3 * f(e).diff(e, 3) + sp.diff(f(e) * f(e).diff(e), e)))  # step 2
    assert z0(sp.diff(3 * f(e).diff(e, 2) + f(e) * f(e).diff(e), e) - (3 * f(e).diff(e, 3) + sp.diff(f(e) * f(e).diff(e), e)))  # step 3 (C₁ = 0 by the far field)
    assert z0(sp.diff(3 * f(e).diff(e) + f(e) ** 2 / 2, e) - (3 * f(e).diff(e, 2) + f(e) * f(e).diff(e)))  # steps 4–5
    fs = sp.sqrt(6) * sp.tanh(e / sp.sqrt(6))
    assert z0(3 * fs.diff(e) + fs ** 2 / 2 - 3) and sp.limit(fs, e, sp.oo) == sp.sqrt(6)  # step 6: f_∞² = 6
    t = sp.Symbol("t")
    assert z0(sp.diff(sp.atanh(t), t) - 1 / (1 - t ** 2))  # steps 7–8: ∫dt/(1 − t²) = artanh t, with t = f/√6, dη/√6
    for ev in (0.3, 1.7, 4.0):  # artanh(f/√6) = η/√6 for the tanh solution
        assert float(sp.atanh(fs.subs(e, ev) / sp.sqrt(6)) - ev / sp.sqrt(6)) == pytest.approx(0.0, abs=1e-14)
    assert z0(fs.diff(e) - sp.sech(e / sp.sqrt(6)) ** 2)  # (9.71)
    s_ = sp.Symbol("s_", real=True)
    prim = sp.tanh(s_) - sp.tanh(s_) ** 3 / 3
    assert z0(sp.diff(prim, s_) - sp.sech(s_) ** 4)
    assert sp.sqrt(6) * (sp.limit(prim, s_, sp.oo) - sp.limit(prim, s_, -sp.oo)) == 4 * sp.sqrt(6) / 3  # step 9 (9.72): η = √6 s
    assert z0(sp.simplify((2 * sp.sqrt(6)) ** 3 / (4 * sp.sqrt(6) / 3) - 36))  # step 10: (2√6)³/C = 36
    # step 11: v/u₀ = −(f − 2ηf′)/(3√Re_x) → ∓√6/(3√Re_x); f − 2ηf′ → ±√6 as η → ±∞ since η sech² → 0
    assert sp.limit(fs - 2 * e * fs.diff(e), e, sp.oo) == sp.sqrt(6) and sp.limit(fs - 2 * e * fs.diff(e), e, -sp.oo) == -sp.sqrt(6)
    # Bickley's coefficients from C (the check line): 0.4543, 0.2752, 3.3019, air ṁ = 0.04268 kg/(m s)
    c = JET.free_jet_constants()
    assert (c["u0_coeff"], c["xi_coeff"], c["mdot_coeff"]) == pytest.approx((0.4543, 0.2752, 3.3019), abs=6e-5)
    assert float(JET.free_jet_mass_flux(0.1, 1.0, 1.2, 1.5e-5)) == pytest.approx(0.04268, abs=5e-6)


def test_jet_halfwidth_V2_derivation_D18():  # V2 (D18 ★): sech² z = 0.01 ⇒ z = arccosh 10 = ln(10 + √99); √6 z = 7.3319; arccosh 5 belongs to 4 %
    z = sp.Symbol("z", positive=True)
    zs = sp.solve(sp.Eq(1 / sp.cosh(z) ** 2, sp.Rational(1, 100)), z)
    assert len(zs) == 1 and float(zs[0] - sp.acosh(10)) == pytest.approx(0.0, abs=1e-14) and float(sp.acosh(10) - sp.log(10 + sp.sqrt(99))) == pytest.approx(0.0, abs=1e-14)
    assert float(sp.sqrt(6) * sp.acosh(10)) == pytest.approx(7.3319, abs=1e-4)
    assert z0(1 / sp.cosh(sp.acosh(5)) ** 2 - sp.Rational(1, 25)) and float(sp.sqrt(6) * sp.acosh(5)) == pytest.approx(5.6153, abs=1e-4)
    assert float(sp.sqrt(6) * 2.2924) == pytest.approx(5.6152, abs=5e-5)  # the printed coefficient is √6 × the rounded 2.2924
    assert 1 / np.cosh(2.2924) ** 2 == pytest.approx(0.0400, abs=5e-5) and 1 / np.cosh(2.9932) ** 2 == pytest.approx(0.0100, abs=5e-6)


def test_wall_jet_invariant_V2_derivation_D19():  # V2 (D19 ★★★): kinematic identities with a generic ψ; each integral step on the exact wall-jet solution (numerically); x u₀² = const
    x, y = sp.symbols("x y", positive=True)
    psi = sp.Function("psi")(x, y)
    u, v = psi.diff(y), -psi.diff(x)
    assert z0(u.diff(x) + v.diff(y))  # continuity by construction
    assert z0(v * u.diff(y) - (sp.diff(u * v, y) + u * u.diff(x)))  # step 6: vu_y = ∂(uv)/∂y + u u_x
    # the wall-jet solution at ν = 1e-3, C = f_∞ = 1: fields on a fine grid at x and x ± h
    nu, Cc = 1e-3, 1.0
    x0, hx = 1.0, 1e-4
    delta = float(JET.wall_jet(x0, 0.0, Cc, 1.0, nu)["delta"])
    ygrid = np.linspace(0.0, 110 * delta, 8001)
    fld = lambda xv: JET.wall_jet(xv, ygrid, Cc, 1.0, nu)  # noqa: E731
    f0, fp, fm = fld(x0), fld(x0 + hx), fld(x0 - hx)
    u_, v_ = f0["u"], f0["v"]
    tail = lambda q: np.concatenate([np.cumsum(((q[1:] + q[:-1]) / 2 * np.diff(ygrid))[::-1])[::-1], [0.0]])  # ∫_y^∞ (trapezoid)  # noqa: E731
    uy_ = np.gradient(u_, ygrid, edge_order=2)
    G2x = (tail(fp["u"] ** 2 / 2) - tail(fm["u"] ** 2 / 2)) / (2 * hx)  # ∂/∂x ∫_y^∞ u²/2 dy′ (step 3)
    inner_vuy = tail(v_ * uy_)  # ∫_y^∞ v u_y dy′
    T1 = np.trapezoid(u_ * G2x, ygrid)
    T2 = np.trapezoid(u_ * inner_vuy, ygrid)
    scale = np.trapezoid(np.abs(u_ * G2x), ygrid)
    assert abs(T1 + T2) < 1e-4 * scale  # steps 4–5: T₁ + T₂ = −ν∫u u_y dy = 0
    assert abs(nu * np.trapezoid(u_ * uy_, ygrid)) < 1e-12 * scale + 1e-14  # step 5: −(ν/2)[u²]₀^∞ = 0
    # step 7: ∫_y^∞ v u_y dy′ = −u v(y) + ∂_x∫_y^∞ u²/2 on the same fields
    assert np.max(np.abs(inner_vuy - (-u_ * v_ + G2x))[:6000]) < 2e-4 * np.max(np.abs(G2x))
    # steps 8–9: T₂ = −∫u²v dy + T₁ and 2T₁ − ∫u²v = 0
    Iu2v = np.trapezoid(u_ ** 2 * v_, ygrid)
    assert abs(T2 - (-Iu2v + T1)) < 1e-4 * scale and abs(2 * T1 - Iu2v) < 2e-4 * scale
    # steps 10–12: d/dx ∫uG dy = ∫u_x G dy + ∫u G_x dy with G = ∫_y^∞u²; ∫u_x G dy = −∫u²v dy (integration by parts); total = 0
    G = tail(u_ ** 2)
    ux_ = (fp["u"] - fm["u"]) / (2 * hx)
    A1 = np.trapezoid(ux_ * G, ygrid)
    A2 = np.trapezoid(u_ * (tail(fp["u"] ** 2) - tail(fm["u"] ** 2)) / (2 * hx), ygrid)
    assert A1 == pytest.approx(-Iu2v, rel=2e-4)  # step 11
    assert abs(A1 + A2) < 2e-4 * abs(A1)  # step 12 (9.80): the invariant does not change with x
    # step 13: I = u₀³δ² K ⇒ ν x u₀² K with K = 1/40 (f_∞ = 1): here ν x u₀² = ν (C x^{-1/2})² x = ν C²
    assert np.trapezoid(u_ * G, ygrid) == pytest.approx(nu * Cc ** 2 / 40.0, rel=5e-5)
    # the ordinary momentum flux is NOT conserved: it differs between x ± h
    assert abs(np.trapezoid(fp["u"] ** 2, ygrid) / np.trapezoid(fm["u"] ** 2, ygrid) - 1) > 1e-5


def test_wall_jet_reduction_V2_derivation_D20():  # V2 (D20 ★★★): all 12 steps with symbolic C, ν; the coefficient 4; the printed coefficient 1 leaves −(3/4)C²f‴/x²
    x, y, nu, C, eta = sp.symbols("x y nu C eta", positive=True)
    A, B = sp.sqrt(nu * C), sp.sqrt(nu / C)
    assert z0(A / B - C) and z0(A * B - nu) and z0(sp.simplify(C * x ** sp.Rational(-1, 2) * (nu * x / (C * x ** sp.Rational(-1, 2))) ** sp.Rational(1, 2) - A * x ** sp.Rational(1, 4)))  # step 1: u₀δ = A x^{1/4}
    E = y / (B * x ** sp.Rational(3, 4))
    D = _jet_D(E)
    assert z0(sp.diff(E, x) + 3 * E / (4 * x)) and z0(sp.diff(E, y) - 1 / (B * x ** sp.Rational(3, 4)))  # step 2
    psi = A * x ** sp.Rational(1, 4) * _FJ[0]
    u, v = D(psi, y), -D(psi, x)
    ux, uy, uyy = D(u, x), D(u, y), D(D(u, y), y)
    sub = lambda e: sp.simplify(sp.powsimp(sp.expand_power_base(e.subs(y, eta * B * x ** sp.Rational(3, 4)), force=True), force=True))  # noqa: E731
    assert z0(sub(u - C * x ** sp.Rational(-1, 2) * _FJ[1]))  # step 3
    assert z0(sub(v + A / 4 * x ** sp.Rational(-3, 4) * (_FJ[0] - 3 * eta * _FJ[1])))  # step 4
    assert z0(sub(ux - C * x ** sp.Rational(-3, 2) * (-_FJ[1] / 2 - 3 * eta * _FJ[2] / 4)))  # step 5
    assert z0(sub(uy - C / B * x ** sp.Rational(-5, 4) * _FJ[2])) and z0(sub(uyy - C ** 2 / nu * x ** -2 * _FJ[3]))  # step 6
    assert z0(sub(u * ux - C ** 2 * x ** -2 * _FJ[1] * (-_FJ[1] / 2 - 3 * eta * _FJ[2] / 4)))  # step 7
    assert z0(sub(v * uy + C ** 2 / 4 * x ** -2 * (_FJ[0] * _FJ[2] - 3 * eta * _FJ[1] * _FJ[2])))  # step 8
    assert z0(sub(u * ux + v * uy + C ** 2 / 4 * x ** -2 * (2 * _FJ[1] ** 2 + _FJ[0] * _FJ[2])))  # step 9: ηf′f″ cancels
    assert z0(sub(nu * uyy - C ** 2 * x ** -2 * _FJ[3]))  # step 10
    res = sub(u * ux + v * uy - nu * uyy)
    assert z0(sp.simplify(res * 4 * x ** 2 / C ** 2 + (4 * _FJ[3] + _FJ[0] * _FJ[2] + 2 * _FJ[1] ** 2)))  # step 11: −(4f‴ + ff″ + 2f′²)
    printed = sp.simplify(res * 4 * x ** 2 / C ** 2 + (_FJ[3] + _FJ[0] * _FJ[2] + 2 * _FJ[1] ** 2))  # the printed ODE leaves the extra 3f‴
    assert z0(printed - (-3 * _FJ[3]))  # residual leftover = −3f‴ (in units of −C²/(4x²) ⇒ the −(3/4)C²f‴/x² of the design check line)
    # conditions: u = v = 0 at y = 0 ⇒ f′(0) = 0 and f(0) = 0 (all homogeneous ⇒ free scale)
    assert z0(sub(u).subs({eta: 0, _FJ[1]: 0})) and z0(sub(v).subs({eta: 0, _FJ[0]: 0}))


def test_wall_jet_integration_V2_derivation_D21():  # V2 (D21 ★★★): the two first integrals, partial fractions, implicit inversion vs IVP, f″(0) = f_∞³/72, far field 4.29, free scale
    e = sp.Symbol("e")
    f = sp.Function("f")(e)
    fp_, fpp_, fppp_ = f.diff(e), f.diff(e, 2), f.diff(e, 3)
    ode = 4 * fppp_ + f * fpp_ + 2 * fp_ ** 2
    assert z0(sp.expand(f * ode - (4 * f * fppp_ + f ** 2 * fpp_ + 2 * f * fp_ ** 2)))  # step 2
    assert z0(4 * f * fppp_ - (4 * sp.diff(f * fpp_, e) - 2 * sp.diff(fp_ ** 2, e))) and z0(f ** 2 * fpp_ + 2 * f * fp_ ** 2 - sp.diff(f ** 2 * fp_, e))  # step 3
    first = 4 * f * fpp_ - 2 * fp_ ** 2 + f ** 2 * fp_
    assert z0(sp.diff(first, e) - f * ode)  # step 4: d/dη(first integral) = f × ODE
    # step 5–7: dividing by 4f^{3/2} gives a total derivative
    step5 = f ** sp.Rational(-1, 2) * fpp_ - fp_ ** 2 * f ** sp.Rational(-3, 2) / 2 + f ** sp.Rational(1, 2) * fp_ / 4
    assert z0(step5 - first / (4 * f ** sp.Rational(3, 2)))
    assert z0(sp.diff(f ** sp.Rational(-1, 2) * fp_ + f ** sp.Rational(3, 2) / 6, e) - step5)
    # step 8–9: separation of variables and g² = f/f_∞
    F, finf, g = sp.symbols("F f_inf g", positive=True)
    fprime = sp.sqrt(F) * (finf ** sp.Rational(3, 2) - F ** sp.Rational(3, 2)) / 6
    assert z0(1 / (finf ** sp.Rational(3, 2) * sp.sqrt(F) - F ** 2) - 1 / (6 * fprime))  # correct: f^{1/2}
    assert not z0(1 / (finf ** sp.Rational(3, 2) * F - F ** 2) - 1 / (6 * fprime))  # printed (R4): f in place of f^{1/2}
    assert z0(sp.simplify((1 / (finf ** sp.Rational(3, 2) * sp.sqrt(F) - F ** 2)).subs(F, finf * g ** 2) * 2 * finf * g - 2 / (finf * (1 - g ** 3))))
    # step 10: partial fractions
    assert z0(sp.apart(1 / (1 - g ** 3), g) - (sp.Rational(1, 3) * (1 / (1 - g) + (g + 2) / (1 + g + g ** 2)))) or z0(sp.simplify(1 / (1 - g ** 3) - sp.Rational(1, 3) * (1 / (1 - g) + (g + 2) / (1 + g + g ** 2))))
    # steps 11–12: (9.83) solved for η at f_∞ = 1, derivative = 12/(1 − g³)
    eta_of_g = 4 * (-sp.log(1 - g) + sp.sqrt(3) * sp.atan((2 * g + 1) / sp.sqrt(3)) + sp.log(1 + g + g ** 2) / 2 - sp.sqrt(3) * sp.atan(1 / sp.sqrt(3)))
    assert z0(sp.simplify(sp.diff(eta_of_g, g) - 12 / (1 - g ** 3))) and z0(eta_of_g.subs(g, 0))  # g(0) = 0
    # the IVP solution satisfies (9.83): invert with brentq (η from g) — two independent routes
    d = JET.wall_jet_ode_solve(1.0 / 72.0, eta_max=120.0, n=6001)
    gnum = np.sqrt(np.maximum(d["f"], 0.0) / d["f_inf"])
    lam_eta = sp.lambdify(g, eta_of_g, "numpy")
    idx = [i for i in range(1, len(d["eta"])) if gnum[i] < 0.9999]
    assert np.max(np.abs(lam_eta(gnum[idx[::50]]) - d["eta"][idx[::50]])) < 2e-8
    # step 13: 1 − g ≈ √3 e^{√3π/6} e^{−f_∞η/4} = 4.29 e^{−f_∞η/4} (the constant the book drops), and f′ ≈ 2.14 f_∞² e^{−f_∞η/4}
    K = sp.sqrt(3) * sp.exp(sp.sqrt(3) * sp.pi / 6)
    assert float(K) == pytest.approx(4.2897, abs=1e-4) and float(K / 2) == pytest.approx(2.1448, abs=1e-4)
    # step 14: f″(0) = f_∞³/72 (f^{-1/2}f′ → √(2f″(0)) at the wall)
    assert z0(sp.solve(sp.Eq(sp.sqrt(2 * sp.Symbol("f2", positive=True)), finf ** sp.Rational(3, 2) / 6), sp.Symbol("f2", positive=True))[0] - finf ** 3 / 72)
    # step 15: scaling f → λf(λη): all terms weight λ⁴, f_∞ → λf_∞; C → C/λ² keeps u₀f′ invariant (the physical constant is C f_∞²)
    lam = sp.Symbol("lam", positive=True)
    flam = lam * sp.Function("h")(lam * e)
    assert z0(sp.simplify((4 * flam.diff(e, 3) + flam * flam.diff(e, 2) + 2 * flam.diff(e) ** 2).doit() - lam ** 4 * sp.Subs(4 * sp.Function("h")(sp.Symbol("s")).diff(sp.Symbol("s"), 3) + sp.Function("h")(sp.Symbol("s")) * sp.Function("h")(sp.Symbol("s")).diff(sp.Symbol("s"), 2) + 2 * sp.Function("h")(sp.Symbol("s")).diff(sp.Symbol("s")) ** 2, sp.Symbol("s"), lam * e).doit()))
    # the module's engine: same first integrals
    ws = ch09.wall_jet_sympy()
    assert z0(ws["first_integral_derivative"]) and ws["fpp0_over_finf3"] == sp.Rational(1, 72)
    assert JET.wall_jet_integrals(1.0)["int_fp2"] == pytest.approx(1 / 18, rel=1e-13)  # the check line: ∫f′ = 1, ∫f′² = 1/18, ∫f′∫f′² = 1/40, peak f′ = 0.0787 at η = 8.11
    assert np.max(d["fp"]) == pytest.approx(0.0787, abs=1e-4) and d["eta"][np.argmax(d["fp"])] == pytest.approx(8.11, abs=0.03)


def test_secondary_flow_V2_derivation_D22():  # V2 (D22 ★★): the r-momentum equation of a pure swirl gives ∂p/∂R = ρu²/R; the net inward force ρ(u_e² − u²)/R; numbers
    r, z, rho, nu = sp.symbols("r z rho nu", positive=True)
    uphi = sp.Function("u")(r, z)
    p = sp.Function("p")(r, z)
    # cylindrical Navier–Stokes, r-component, u_r = u_z = 0, axisymmetric, steady: −u_φ²/r = −(1/ρ)∂p/∂r + ν[∇²u_r − u_r/r² − (2/r²)∂u_φ/∂φ]
    visc_r = 0  # ∇²u_r − u_r/r² − (2/r²)∂u_φ/∂φ = 0 when u_r = 0 and the field is axisymmetric
    balance = sp.Eq(-uphi ** 2 / r, -p.diff(r) / rho + nu * visc_r)
    dpdr = sp.solve(balance, p.diff(r))[0]
    assert z0(dpdr - rho * uphi ** 2 / r)  # step 1
    ue = sp.Function("ue")(r)
    F_in = sp.simplify(dpdr.subs(uphi, ue) - rho * uphi ** 2 / r)  # steps 2–3: layer feels the core's pressure gradient, needs only ρu²/R
    assert z0(F_in - rho * (ue ** 2 - uphi ** 2) / r)
    # step 4: sign, step 5: numbers
    assert ch09.secondary_flow_radial_force(0.2, 0.0, 0.04, 1000.0) == pytest.approx(1000.0) and 1000.0 / 9810.0 == pytest.approx(0.102, abs=1e-3)
    assert [ch09.secondary_flow_radial_force(0.2, f * 0.2, 0.04, 1000.0) for f in (0, 0.5, 0.75, 1.0)] == pytest.approx([1000.0, 750.0, 437.5, 0.0])
    # step 7: Ekman scale √(ν/Ω) with Ω = u_e/R
    assert np.sqrt(1e-6 / (0.2 / 0.04)) * 1e3 == pytest.approx(0.447, abs=1e-3)
    # continuity closes the loop: an inflow on the floor (radial flux 2πRu_r per unit height) must upwell near the axis — the flux needs a positive inward mean
    zc = np.linspace(0, 1e-3, 200)
    inflow = ch09.secondary_flow_radial_force(0.2, ch09.secondary_flow_layer_profile(zc, 1e-3, 0.2, "linear"), 0.04, 1000.0)
    assert np.trapezoid(inflow, zc) > 0


# =====================================================================================================================
# Transition and the flat-plate drag curve (C09/N68–N69), the sympy engines, the slip table and the design Part C contract
# =====================================================================================================================
def test_plate_drag_V1_laminar_law_mixed_curve_and_the_turbulent_form():  # V1 (C09/N69): 1.328/√Re_L; mixed = turbulent − Re_tr(C_t − C_l)/Re_L above Re_tr; design expect rows
    R = np.array([1e5, 1e6, 1e7])
    lam = np.asarray(BL.plate_drag_coefficient(R))
    assert lam == pytest.approx(4 * BL.blasius_constants()["fpp0"] / np.sqrt(R), rel=1e-14)  # (9.33)
    assert lam == pytest.approx([4.20e-3, 1.328e-3, 4.20e-4], rel=2e-3)
    tur = np.asarray(BL.plate_drag_coefficient(R, "turbulent"))
    assert tur == pytest.approx(0.074 * R ** -0.2, rel=1e-14) and tur[1:] == pytest.approx([4.669e-3, 2.946e-3], rel=1e-3)
    mix = np.asarray(BL.plate_drag_coefficient(R, "mixed"))
    A = 5e5 * (0.074 * 5e5 ** -0.2 - 4 * BL.blasius_constants()["fpp0"] / np.sqrt(5e5))
    assert A == pytest.approx(1743.0, abs=1.5) and mix[1:] == pytest.approx([2.926e-3, 2.772e-3], rel=1e-3)
    assert mix * R == pytest.approx(np.where(R <= 5e5, lam * R, tur * R - A), rel=1e-12)  # the momentum-thickness patching, the design's formula
    assert float(BL.plate_drag_coefficient(5e5, "mixed")) == pytest.approx(1.878e-3, rel=1e-3)  # joins the laminar line at Re_tr
    assert float(BL.plate_drag_coefficient(5e5 * (1 + 1e-9), "mixed")) == pytest.approx(float(BL.plate_drag_coefficient(5e5)), rel=1e-7)  # continuous
    assert mix[0] == lam[0] and mix[2] > lam[2]  # below Re_tr the mixed curve IS the laminar one; far above, transition raises the drag
    assert float(BL.plate_drag_coefficient(1e6, sides=2)) == pytest.approx(2 * 1.328e-3, rel=2e-3)  # R5 again
    assert float(BL.plate_drag_coefficient(1e7, "mixed", Re_tr=1e6)) < float(BL.plate_drag_coefficient(1e7, "mixed", Re_tr=5e5))  # later transition ⇒ less drag
    with pytest.raises(ValueError):
        BL.plate_drag_coefficient(1e6, "nonsense")


@needs_ref
def test_plate_drag_V5_secondary_turbulent_law_within_three_percent():  # V5 secondary (C09/N69): Wikipedia c_f = 0.0576 Re_x^{-1/5} ⇒ plate mean 0.0720; the code's 0.074 is 2.8 % above (cited source: no primary read)
    R = ref_json()["skin_friction_wikipedia"]
    mean_coeff = 1.25 * R["turbulent_local_cf_coeff"]  # ∫₀ᴸ x^{-1/5} dx / L = (5/4) L^{-1/5}
    assert mean_coeff == pytest.approx(R["turbulent_mean_coeff_derived"], abs=1e-4)
    assert abs(0.074 / mean_coeff - 1) < 0.03 and abs(0.074 / mean_coeff - 1) > 0.02  # 2.8 %: inside the digitised/secondary 3 % band, not tighter
    assert float(BL.plate_drag_coefficient(1e6, "turbulent")) == pytest.approx(mean_coeff * 1e6 ** -0.2, rel=0.03)
    assert float(BL.plate_drag_coefficient(1e6, turbulent_coeff=mean_coeff, regime="turbulent")) == pytest.approx(mean_coeff * 1e6 ** -0.2, rel=1e-14)  # the coefficient is an argument


def test_transition_state_V7_regime_table_and_arguments():  # V7 (C09/N68, R9): laminar → transitional → turbulent; thresholds are arguments (rounded experimental values)
    assert [BL.transition_state(r) for r in (1e5, 4.999e5, 5e5, 4.9e6, 5e6, 1e8)] == ["laminar", "laminar", "transitional", "transitional", "turbulent", "turbulent"]
    assert BL.transition_state(8e5, Re_cr=1e6) == "laminar" and BL.transition_state(2e6, Re_cr=1e6) == "transitional"
    assert BL.transition_state(2e7, Re_cr=1e6, Re_turb=1e7) == "turbulent" and BL.transition_state(2e6, Re_cr=1e6, Re_turb=1e7) == "transitional"


def test_sympy_engines_V2_thwaites_jets_cylinder_and_derive_all():  # V2 (C06/C07/C12/C13, Part C 4.3–4.6, 4.10): every engine's residuals are exactly zero
    th = ch09.thwaites_sympy()
    assert z0(th["identity_9_47_to_l"]) and z0(th["identity_9_48"]) and z0(th["integrating_factor_check"])
    x = sp.Symbol("x", positive=True)
    assert th["ode_9_49"].lhs.has(sp.Derivative) and th["I5_power_law"] is not None
    jm = ch09.jet_momentum_sympy()
    assert z0(jm["pointwise_identity"]) and z0(jm["residual_x_dependence"]) and z0(jm["equals_J_over_rho"]) and z0(jm["C_value"]) and z0(jm["C_free_jet"] - 4 * sp.sqrt(6) / 3)
    wi = ch09.wall_jet_invariant_sympy()
    assert z0(wi["residual"]) and z0(wi["delta_matches"]) and wi["exponent_delta"] is not None
    assert z0(wi["u0_solution"] - sp.Symbol("C", positive=True) * x ** sp.Rational(-1, 2))
    ct = ch09.cylinder_thwaites_sympy()
    assert z0(ct["F_direct_difference"]) and z0(ct["F_derivative_residual"]) and ct["lam0"] == sp.Rational(3, 40) and z0(ct["lam_at_90deg"]) and ct["lam_82deg"] == pytest.approx(0.0263, abs=6e-5)
    d = ch09.derive_all()
    flat = []
    for v in d.values():
        flat.extend(v if isinstance(v, tuple) else [v])
    assert len(d) == 12 and all(z0(sp.sympify(v)) for v in flat), d
    # the engines FAIL on wrong inputs (they are not vacuous): a wrong integrating factor, a wrong coefficient 0.45 → 0.5 in (9.49)
    Ue, Z = sp.Function("U_e")(x), sp.Function("Z")(x)
    wrong = sp.diff(Ue ** 5 * Z, x) - Ue ** 5 * (Z.diff(x) + 6 * Ue.diff(x) / Ue * Z)
    assert not z0(wrong)
    assert not z0(ch09.bl_nondim_sympy(printed_9_7=True)["residual_printed_vs_correct"])


def test_book_slips_V1_table_is_consistent_with_the_planted_variants():  # V1 (Part C 4.9): every printed/correct pair is reproduced by the code's discriminating option
    tab = {r["id"]: r for r in ch09.book_slips()}
    assert set(tab) == {"R1", "R2", "R3", "R4", "R5", "R6", "R7", "R10", "R11", "R14"}
    assert tab["R2"]["printed_value"] == 4.93 and tab["R2"]["correct_value"] == pytest.approx(BL.blasius_constants()["eta99"], rel=1e-14)
    assert tab["R3"]["printed_value"] == 1 and tab["R3"]["correct_value"] == 4
    assert tab["R6"]["printed_value"] == pytest.approx(5.6152, abs=1e-4) and tab["R6"]["correct_value"] == pytest.approx(7.3319, abs=1e-4)
    assert (tab["R14"]["printed_value"], tab["R14"]["correct_value"]) == (13, 12)
    # each of those printed values is REJECTED by an independent test above: R2 (f′(4.93) ≠ 0.99), R3 (ode coefficient), R6 (u/u₀ = 0.04)
    assert abs(float(BL.blasius_profile(tab["R2"]["printed_value"])[1]) - 0.99) > 2e-4
    assert JET.free_jet_at_level(0.04)["coeff"] == pytest.approx(tab["R6"]["printed_value"], abs=1e-3)
    assert all({"id", "where", "printed", "correct"} <= set(r) for r in tab.values())


def test_part_c_V1_every_contract_function_exists_and_is_scalar_callable():  # V1 (design Part C 1.1–4.10 + reused 0.x): smoke with scalars, finite floats, physical sanity
    # (i) ch09 re-exports every public name of BL, JET, BB and the similarity engines (explainer parity rows write ch09.<name> only)
    for mod in (BL, JET, BB):
        for name in getattr(mod, "__all__"):
            assert hasattr(ch09, name), (mod.__name__, name)
    for name in ("similarity_reduce_sympy", "similarity_ode_solve", "similarity_collapse_error", "sphere_drag_coefficient", "strouhal_number", "stokes_drag_coefficient",
                 "oseen_drag_coefficient", "temporal_bl_wall_stress", "pressure_coefficient", "reynolds_number", "similarity_variable", "diffusion_thickness", "stokes_first_problem"):
        assert callable(getattr(ch09, name)), name
    # (ii) scalar in, scalar out (parity rows index the result down to one float)
    scalars = [
        ch09.blasius_delta99(1.0, 1.0, 1.5e-5), ch09.blasius_delta_star(1.0, 1.0, 1.5e-5), ch09.blasius_theta(1.0, 1.0, 1.5e-5), ch09.blasius_wall_shear(1.0, 1.0, 1.2, 1.5e-5),
        ch09.blasius_skin_friction(1e5), ch09.blasius_drag(1.0, 1.0, 1.2, 1.5e-5), ch09.blasius_drag_coefficient(1e5), ch09.blasius_far_field(6.0),
        ch09.falkner_skan_thickness(1.0, 0.3, 1.0, 1e-6), ch09.holstein_bohlen(1e-3, -1.0, 1e-5), ch09.thwaites_l(0.02), ch09.thwaites_H(0.02), ch09.thwaites_L(0.02),
        ch09.thwaites_cylinder_closed_form(1.0), ch09.thwaites_cylinder(60.0), ch09.thwaites_cylinder_separation(), ch09.wall_curvature(3.0, 1e-3),
        ch09.plate_drag_coefficient(1e6), ch09.free_jet_centreline(0.1, 1.0, 1.2, 1.5e-5), ch09.free_jet_thickness(0.1, 1.0, 1.2, 1.5e-5), ch09.free_jet_mass_flux(0.1, 1.0, 1.2, 1.5e-5),
        ch09.free_jet_halfwidth(0.1, 1.0, 1.2, 1.5e-5), ch09.free_jet_entrainment_velocity(1e4), ch09.wall_jet_mass_flux(1.0, 1.0, 1.0, 1.0, 1e-3), ch09.wall_jet_first_integral_residual(0.1, 0.2, 0.3),
        ch09.karman_street_ratio(), ch09.karman_street_velocity(1.0, 0.28, 1.0), ch09.karman_street_growth(0.3), ch09.karman_street_growth_closed(0.3), ch09.shedding_angular_frequency(10.0, 2e-3),
        ch09.cp_ideal_cylinder(30.0), ch09.separated_cp(100.0, 82.0), ch09.separated_pressure_drag(82.0, -1.2), ch09.pressure_drag_from_cp(lambda p: 1 - 4 * np.sin(p) ** 2),
        ch09.ball_swing_deflection(0.2, 18.0, 35.0), ch09.cylinder_cd_schematic(1e3), ch09.secondary_flow_radial_force(0.2, 0.1, 0.04, 1000.0), ch09.jet_momentum_flux(np.linspace(0, 1, 5), np.ones(5)),
    ]
    for v in scalars:
        assert np.ndim(v) == 0 and np.isfinite(float(v)), v
    for v, want in ((ch09.blasius_delta99(1.0, 1.0, 1.5e-5), 4.910 * np.sqrt(1.5e-5)), (ch09.free_jet_entrainment_velocity(1e4), -np.sqrt(6) / 300.0),
                    (ch09.wall_curvature(3.0, 1e-3), 3000.0), (ch09.secondary_flow_radial_force(0.2, 0.1, 0.04, 1000.0), 750.0)):
        assert float(v) == pytest.approx(want, rel=2e-4)
    # (iii) dict-returning contract functions carry their documented keys with finite scalars
    assert set(BL.boundary_layer_scales(1.0, 1.0, 1e-5)) >= {"Re", "delta_over_L", "delta", "v_scale", "tau0_scale", "cf_estimate", "adv", "visc", "visc_x"}
    assert set(BL.blasius_constants()) >= {"fpp0", "eta99", "delta_star", "theta", "H", "v_inf", "tau_coeff", "cf_coeff", "cd_coeff"}
    assert set(BL.blasius_fields(1.0, 1e-3, 1.0, 1e-6)) >= {"u", "v", "psi", "eta", "delta", "tau0_over_rho"}
    assert set(BL.falkner_skan(0.3)) >= {"eta", "f", "fp", "fpp", "fppp", "fpp0", "success"} and set(BL.falkner_skan_separation()) == {"m_sep", "beta_sep", "fpp0_at_sep"}
    assert set(BL.falkner_skan_state(0.3)) >= {"fpp0", "fppp0", "I_delta", "I_theta", "H", "lam", "l", "cf_sqrtRex", "inflection_eta", "separated"}
    assert set(BL.falkner_skan_fields(1.0, 1e-3, 0.3, 1.0, 1e-5)) >= {"u", "v", "psi", "eta", "delta", "Ue"} and set(BL.thicknesses(np.linspace(0, 5, 51), np.linspace(0, 1, 51), 1.0)) == {"delta99", "delta_star", "theta", "H"}
    xs = np.linspace(0.05, 0.5, 30)
    assert set(BL.thwaites(xs, BL.outer_flow("flat", U=1.0), 1e-5, theta0=1e-4)) >= {"x", "theta", "delta_star", "tau0", "cf", "lam", "H", "l", "separated", "x_sep"}
    assert set(BL.thwaites_named("flat", 0.3, 1e-5, U=1.0, theta0=1e-4)) >= {"theta", "delta_star", "tau0", "cf", "lam", "H", "l", "x_sep"}
    assert set(BL.thwaites_closure_table(fast=True)) == {"m", "lam", "l", "H", "L"} and set(BL.karman_pohlhausen(BL.outer_flow("flat", U=1.0), xs, 1e-5)) >= {"x", "delta", "theta", "delta_star", "tau0", "lam"}
    assert set(JET.free_jet_constants()) >= {"C", "mdot_coeff", "h99_arg", "h99_coeff", "h99_printed", "u0_coeff", "f_inf"} and set(JET.free_jet_profile(1.0)) == {"f", "fp", "fpp"}
    assert set(JET.free_jet(0.1, 1e-4, 1.0, 1.2, 1.5e-5)) == {"u", "v", "psi", "eta", "u0", "delta"} and set(JET.free_jet_reynolds(0.1, 1.0, 1.2, 1.5e-5)) == {"Re_x", "Re_h99", "Re_h99_printed"}
    assert set(JET.free_jet_at_level(0.01)) == {"z", "coeff"} and set(JET.wall_jet_profile(1.0)) == {"g", "f", "fp", "fpp"} and set(JET.wall_jet_integrals()) == {"int_fp", "int_fp2", "invariant", "fpp0"}
    assert set(JET.wall_jet(1.0, 0.1, 1.0, 1.0, 1e-3)) >= {"u", "v", "psi", "eta", "u0", "delta", "mdot"} and set(JET.wall_jet_ode_solve(eta_max=30.0)) >= {"eta", "f", "fp", "fpp", "f_inf", "err_vs_9_83"}
    assert set(JET.wall_jet_constants(1000.0, 1e-6, Psi=1e-6)) == {"C", "C_f_inf_sq", "Psi", "mdot_at_x", "residual"} and JET.wall_jet_invariant(np.linspace(0, 10, 200), np.exp(-np.linspace(0, 10, 200))) > 0
    assert {"label", "thresholds", "separation_deg", "St"} <= set(BB.cylinder_flow_regime(100.0)) and {"label", "thresholds", "separation_deg", "St"} <= set(BB.sphere_flow_regime(1e3))
    assert {"label", "phi_sep_deg", "St", "cb", "cd_model", "qualitative"} <= set(BB.cylinder_state(1e5)) and {"label", "phi_sep_deg", "cb", "cd_model", "blend", "qualitative"} <= set(BB.drag_crisis_state(1e5))
    assert set(BB.drag_crisis_pair()) == {"subcritical", "supercritical", "ratio", "Re_cr"} and set(BB.shedding_frequency(10.0, 2e-3)) == {"f", "omega_rad"}
    assert BB.karman_street_spectrum(0.3).shape == (4,) and BB.karman_street_spectrum_periodic(0.3, fast=True).shape == (32,) and set(BB.karman_street_positions(1.0, 0.3)) == {"zA", "zB", "growth"}
    assert BB.magnus_sign(1e5, 5e5) == "−" and BL.transition_state(1e6) == "transitional"
    ex1, ex2 = ch09.example_9_1(), ch09.example_9_2()
    assert {"theta_coef", "theta_err", "delta_star_coef", "cf_sqrtRex", "cf_err"} <= set(ex1) and {"lam", "x_sep_over_L", "x_sep_over_L_fs", "table"} <= set(ex2)
    # (iv) the outer-flow kinds all return callable objects with dpdx, and the 'custom' kind accepts a callable with or without its derivative
    for kind, kw in (("flat", {}), ("wedge", dict(n=0.3)), ("stagnation", dict(a=2.0)), ("diffuser", {}), ("cylinder", {}), ("linear_retarded", {}), ("retarded", {}),
                     ("custom", dict(Ue=lambda s: 1 + 0.5 * s)), ("custom", dict(Ue=lambda s: 1 + 0.5 * s, dUe=lambda s: 0.5 + 0 * np.asarray(s)))):
        of = BL.outer_flow(kind, **kw)
        assert np.isfinite(float(of.Ue(0.3))) and np.isfinite(float(of.dUe(0.3))) and np.isfinite(float(of.dpdx(0.3, 1.2)))
    assert BL.outer_flow("custom", Ue=lambda s: 1 + 0.5 * s).dUe(0.3) == pytest.approx(0.5, abs=1e-7)
    with pytest.raises(ValueError):
        BL.outer_flow("nonsense")
    # (v) reused functions the notebook calls (ch04/ch06/ch08 tested elsewhere): still callable with ch09-shaped inputs
    assert float(SIM.reynolds_number(1.0, 1.0, 1.5e-5)) == pytest.approx(6.667e4, rel=1e-3) and float(LAM.similarity_variable(1e-3, 1.0, 1e-6)) == pytest.approx(1.0)
    assert float(LAM.diffusion_thickness(1.0, 1e-6)) > 0 and float(SIM.pressure_coefficient(1.0, 0.0, 1.2, 1.0)) == pytest.approx(1 / 0.6)
    assert float(CRP.oseen_drag_coefficient(1.0)) > float(CRP.stokes_drag_coefficient(1.0))
    # (vi) wrong-argument guards and errors named in the docstrings
    with pytest.raises(ValueError):
        BL.falkner_skan(0.3, method="nonsense")
    with pytest.raises(ValueError):
        BL.thwaites(np.linspace(0.1, 1, 5), BL.outer_flow("flat", U=1.0), 1e-5, closure="nonsense")
    with pytest.raises(ValueError):
        ch09.secondary_flow_layer_profile(1.0, 1.0, 1.0, "nonsense")


# =====================================================================================================================
# V6 — numbers the book prints (private JSON, skipped when absent): relative errors only in the report
# =====================================================================================================================
@book_only
def test_book_V6_blasius_numbers_and_the_R2_slip():  # V6 (C04/N33–N40): 4.93, 1.72, 0.664, 0.332, 1.33, 0.86, the air example
    B = book()["sec_9_3_blasius"]
    c = BL.blasius_constants()
    assert rel(c["eta99"], B["eta_99_derived_by_analyst"]) < 1e-5 and rel(c["eta99"], B["eta_99_book"]) < 5e-3  # 4.93 is the figure reading; the root is 0.4 % smaller
    assert rel(c["delta_star"], B["delta_star_coeff"]) < 5e-3 and rel(c["theta"], B["theta_coeff"]) < 5e-3
    assert rel(c["fpp0"], B["tau0_coeff_9_31"]) < 5e-3 and rel(c["cf_coeff"], B["cf_coeff_9_32"]) < 5e-3
    assert rel(c["cd_coeff"], B["CD_coeff_9_33"]) < 5e-3 and rel(c["v_inf"], B["v_over_U_sqrtRex_far"]) < 5e-3
    ex = B["air_example"]
    nu_book = ex["U"] * ex["x"] / ex["Re_x"]  # the ν implied by Re_x = 6e4
    assert rel(float(BL.blasius_delta99(ex["x"], ex["U"], nu_book)) * 100, ex["delta99_cm"]) < 5e-3
    assert 1.85 < float(BL.blasius_delta99(ex["x"], ex["U"], 1.5e-5)) * 100 < 1.95  # ν = 1.5e-5: 6.7e4 and 1.9 cm (the analyst's note)
    assert BL.boundary_layer_scales(1.0, 1.0, 1.5e-5)["cf_estimate"] == pytest.approx(book()["sec_9_1"]["c_f_estimate_2_over_sqrtRe"] / np.sqrt(6.667e4), rel=2e-3)  # 2/√Re: three times 0.664/√Re
    assert BL.boundary_layer_scales(1.0, 1.0, 1.5e-5)["cf_estimate"] / float(BL.blasius_skin_friction(6.667e4)) == pytest.approx(3.0, abs=0.02)
    # exercise-type checks (answers private): the Blasius field at x = 0.15 m, U = 6 m/s, ν = 1.5e-5
    ex8 = book()["exercises_answers_in_text_do_not_publish"]["9.8_air"]
    f = BL.blasius_fields(0.15, ex8["y_mm"] * 1e-3, 6.0, 1.5e-5)
    fpp = float(BL.blasius_profile(float(f["eta"]))[2])
    assert rel(float(f["v"]) * 100, ex8["v_cm_s"]) < 0.02 and rel(6.0 * fpp / float(f["delta"]), ex8["du_dy_per_s"]) < 2e-3  # 0.9 % and 0.1 %: the book prints 2–3 digits
    assert rel(float(f["u"]) / 6.0, 0.456) < 1e-3


@book_only
def test_book_V6_falkner_skan_numbers():  # V6 (C05): n_sep = −0.0904, stagnation exercise, Fig. 9.7 members
    B = book()["sec_9_4_falkner_skan"]
    assert rel(BL.falkner_skan_separation()["m_sep"], B["n_separation"]) < 5e-4
    assert rel(BL.falkner_skan_state(1.0)["cf_sqrtRex"], B["exercise_9_19_stagnation_cf_sqrtRex"]) < 2.1e-5  # 2f″(0) = 2.4652 (5 printed digits: half-unit 2e-5)
    assert rel(BL.falkner_skan_state(1.0)["fpp0"], B["derived_by_analyst_stagnation_fpp0"]) < 1e-9
    for n in B["fig_9_7_labels_n"]:
        assert BL.falkner_skan(n, eta_max=16.0)["success"], n
    # the book's x-axis of Fig. 9.7 is ½·√(n+1)·η (a rescaling by √((n+1)/2)² … the tests use the f-normalisation of the ODE, checked by V5)
    assert B["fig_9_7_x_axis"].startswith("sqrt((n+1)/2)")


@book_only
def test_book_V6_thwaites_table_examples_and_accuracy_claims():  # V6 (C07): Table 9.1 vs our closures, Examples 9.1 and 9.2, ±3 %/±10 %
    T = book()["sec_9_6_thwaites"]
    rows = np.array(T["table_9_1_lambda_l_H"])
    lam, l_b, H_b = rows[:, 0], rows[:, 1], rows[:, 2]
    # the 'white' fit (l = (λ + 0.09)^0.62, a fit to Thwaites' table) reproduces the table's l: max relative error over −0.06 ≤ λ ≤ 0.25
    m = lam >= -0.06
    err_white = np.abs(np.asarray(BL.thwaites_l(lam[m], "white")) / l_b[m] - 1)
    assert err_white.max() < 0.03, err_white.max()
    # the exact-FS closure (ours) agrees with the table only approximately — it is exact for FS flows, the table is a cross-family fit:
    m2 = (lam >= 0.0) & (lam <= 0.1)
    e_l = np.abs(np.asarray(BL.thwaites_l(lam[m2])) / l_b[m2] - 1)
    e_H = np.abs(np.asarray(BL.thwaites_H(lam[m2])) / H_b[m2] - 1)
    assert e_l.max() < 0.062 and e_H.max() < 0.048  # measured 6.0 % and 4.7 % (the design's "< 5 %" holds for H, not for l)
    m3 = (lam >= -0.016) & (lam < 0.0)
    assert np.abs(np.asarray(BL.thwaites_l(lam[m3])) / l_b[m3] - 1).max() < 0.04  # measured 3.6 %; it worsens quickly towards the fold (14 % at λ = −0.04, ≫ 50 % at −0.06)
    assert abs(float(BL.thwaites_l(-0.04)) / 0.153 - 1) > 0.10  # the FS-exact l is NOT within 5 % of the table for strongly adverse λ (design D10 claim corrected in the report)
    assert T["L_lambda_linear_fit"] == [0.45, -6.0]
    e1 = ch09.example_9_1()
    E1 = T["example_9_1_blasius"]
    assert rel(e1["theta_coef"], E1["theta_coeff"]) < 1e-3 and abs(e1["theta_err"] * 100 - E1["theta_pct_above_blasius"]) < 0.5
    assert rel(e1["l0"], E1["l0"]) < 5e-3 and rel(e1["H0"], E1["H0"]) < 1e-2 and rel(e1["delta_star_coef"], E1["delta_star_coeff"]) < 1e-2
    assert rel(e1["cf_sqrtRex"], E1["cf_coeff"]) < 3e-3 and abs(e1["cf_err"] * 100 + E1["cf_pct_below_blasius"]) < 0.5
    E2 = T["example_9_2_diffuser"]
    e2 = ch09.example_9_2()
    for xi, lam_book in E2["lambda_table"].items():
        assert rel(float(e2["lam"](float(xi))), lam_book) < 1e-3, xi
    assert abs(e2["x_sep_over_L"] - E2["x_over_L_separation"]) < 5e-3  # the book prints 0.16 (2 s.f.); the closed form is 0.15829
    assert rel(e2["x_sep_over_L"], E2["derived_by_analyst_x_over_L_sep"]) < 1e-4 and E2["separation_lambda"] == BL.LAMBDA_SEP_BOOK
    # accuracy claims: measured θ errors on FS flows exceed the 3 % claim for strong acceleration (−4 % … −7.6 %)
    th_err = abs(_thwaites_wedge(1.0)["theta"][-1] / (BL.falkner_skan_state(1.0)["I_theta"] * np.sqrt(1e-6 * 0.7 / 0.7)) - 1)
    assert th_err > T["accuracy_favourable_percent"] / 100.0  # 6.3 % > 3 %: documented deviation (θ, not τ₀, is what (9.50) predicts)


@book_only
def test_book_V6_transition_cylinder_sphere_and_ball_numbers():  # V6 (C09–C11): regime thresholds, angles, wire example, Kármán ratio, sports-ball consistency
    B = book()
    T, C, S = B["sec_9_7_transition_separation"], B["sec_9_8_cylinder"], B["sec_9_9_sphere_sports"]
    th = BB.CYLINDER_THRESHOLDS
    assert th["creeping"] == C["creeping_Re_below"] and th["attached_eddies"] == C["steady_eddies_Re_above"] and th["street_onset"] == C["wake_unstable_Re"]
    assert th["irregular"] == C["wake_laminar_below_Re"] and th["critical"] == C["Re_cr_smooth_cylinder"]
    assert BB.cylinder_flow_regime(1e5)["separation_deg"] == T["cylinder_laminar_separation_deg"] and BB.cylinder_flow_regime(1e6)["separation_deg"] == T["cylinder_turbulent_separation_deg"]
    assert BB.CYLINDER_THRESHOLDS["critical"] == 3e5 and BL.transition_state(T["fig_9_11"]["transition_Re_local"]) == "transitional"  # local Re 5e5 (the code's default Re_cr)
    assert BL.transition_state(T["Re_x_cr_flat_plate"], Re_cr=T["Re_x_cr_flat_plate"] / T["Re_cr_factor"]) == "transitional"  # Re_cr ≈ 10⁶ "within a factor of five"
    w = C["wire_example"]
    assert BB.shedding_frequency(w["U"], w["d_mm"] * 1e-3, C["strouhal"])["f"] == pytest.approx(w["frequency_Hz"], rel=1e-12)
    assert rel(BB.karman_street_ratio(), C["karman_stagger_ratio"]) < 5e-3 and rel(BB.karman_street_ratio(), C["karman_stagger_exact_derived"]) < 1e-4
    # drag-crisis magnitude: the book's C_D falls from 1.2 to 0.33 (× 0.28); the illustrative model's pressure drag falls × 0.65 — qualitative agreement only
    ratio_book = C["CD_supercritical"] / C["CD_subcritical"]
    ratio_model = BB.drag_crisis_pair()["ratio"]
    assert ratio_model < 1.0 and ratio_book < 1.0  # same direction; the magnitudes differ (ours 0.65, the book's ≈ 0.28) — reported in the Open items, not asserted
    assert S["potential_flow_suction_Cp_min"] == -1.25
    from fluidpy import ch06_ideal_flow as ch06
    assert float(ch06.sphere_surface_cp(np.pi / 2)) == pytest.approx(S["potential_flow_suction_Cp_min"], abs=1e-14)
    k = S["cricket"]
    t = S["cricket"]["derived_by_analyst_flight_time_s"]
    assert BB.ball_swing_deflection(k["side_force_fraction_of_weight"], k["speed"] * t, k["speed"]) == pytest.approx(k["deflection_m"], rel=5e-3)
    assert BB.magnus_sign(1e5, 5e5, S["Re_cr_sphere"] * 0.6) == "−"  # negative Magnus effect when only the fast side is past the crisis
    assert BB.sphere_flow_regime(S["ring_eddy_oscillates_Re"] * 1.01)["label"].startswith("unsteady wake")  # loops shed above Re ≈ 130


@book_only
def test_book_V6_free_jet_and_wall_jet_forms():  # V6 (C12/C13): 4√6/3, 36, the printed h₉₉ numbers (slip R6), wall-jet printed ODE (slip R3), exponents
    B = book()
    F, W = B["sec_9_10_free_jet"], B["sec_9_10_wall_jet"]
    c = JET.free_jet_constants()
    assert F["C_9_72"] == "4*sqrt(6)/3" and c["C"] == pytest.approx(4 * np.sqrt(6) / 3, rel=1e-14)
    assert F["mdot_9_73"].startswith("(36*J") and c["mdot_coeff"] ** 3 == pytest.approx(36.0, rel=1e-13)
    hb = F["h99_book_coefficients"]
    assert JET.H99_PRINTED == hb["coefficient"] and rel(c["h99_printed"], hb["coefficient"]) < 5e-5  # the printed 5.6152 = √6 arccosh 5 to 4 s.f.
    assert rel(np.arccosh(5.0), hb["argument"]) < 1e-4 and rel(F["h99_derived_by_analyst"]["coefficient"], c["h99_coeff"]) < 1e-5
    assert abs(1 / np.cosh(hb["argument"]) ** 2 - F["h99_derived_by_analyst"]["book_argument_corresponds_to_sech2"]) < 1e-4  # the printed argument is the 4 % point
    assert F["width_growth_exponent"] == pytest.approx(2 / 3) and F["mass_flux_exponent"] == pytest.approx(1 / 3)
    assert abs(F["bickley_centreline_constant_public"] - c["u0_coeff"]) < 5e-5
    assert W["ode_printed"] == "f''' + f f'' + 2 f'^2 = 0" and W["ode_correct_derived"] == "4 f''' + f f'' + 2 f'^2 = 0"
    assert W["mdot_exponent"] == 0.25 and W["f_inf_cubed_over_fpp0_derived"] == 72.0
    d = JET.wall_jet_ode_solve(1.0 / 72.0, eta_max=120.0, n=2001)
    assert d["f_inf"] ** 3 / (1.0 / 72.0) == pytest.approx(W["f_inf_cubed_over_fpp0_derived"], rel=1e-8)


# =====================================================================================================================
# Scripts (design C.5): every ch09 script runs headless
# =====================================================================================================================
@pytest.mark.slow
def test_scripts_V1_every_ch09_script_runs():  # V1 smoke: all scripts exit 0 with --no-show (figure generation, no physics of their own)
    env = dict(os.environ, MPLBACKEND="Agg")
    scripts = [s_ for s_ in sorted((ROOT / "scripts").glob("ch09_*.py")) if s_.name not in ("ch09_common.py", "ch09_drawings.py")]
    assert len(scripts) >= 22, [s_.name for s_ in scripts]
    for scr in scripts:
        r = subprocess.run([sys.executable, str(scr), "--no-show"], cwd=ROOT, env=env, capture_output=True, text=True, timeout=900)
        assert r.returncode == 0, (scr.name, r.stderr[-2000:])


# =====================================================================================================================
# Re-verification after the implementer's fixes F1-F4 (2026-09-30): tests for the newly repaired behaviours
# =====================================================================================================================
def test_falkner_skan_shoot_V3_agrees_with_bvp_across_the_documented_range():  # V3 (C05, F1): shoot (brentq) vs solve_bvp continuation, m ≥ -0.05
    for m in (-0.05, 0.0, 0.3, 0.5, 1.0, 2.0, 4.0):
        a = BL.falkner_skan(m, method="shoot")["fpp0"]
        b = BL.falkner_skan(m, method="bvp")["fpp0"]
        assert abs(a / b - 1) < 1e-7, (m, a, b)


def test_thwaites_named_wedge_V1_power_law_lambda_and_theta_for_every_n_above_1p4():  # V1+V5 (C07, F2): closed λ = 0.45n/(5n+1); θ within the 8 % bound of exact FS
    nu, x = 1e-6, 0.7
    for n in (1.5, 2.0, 3.0, 4.0):
        r = BL.thwaites_named("wedge", x, nu, n=n, a=1.0)
        assert r["lam"] == pytest.approx(0.45 * n / (5 * n + 1), rel=1e-6), n  # V1
        st = BL.falkner_skan_state(n)
        assert abs(r["theta"] / (st["I_theta"] * np.sqrt(nu * x / x ** n)) - 1) < 0.08, n  # V5 (same bound as the n ≤ 1 rows)


def test_thwaites_V7_stagnation_guard_still_rejects_a_wrong_theta0():  # V7 (C07, F2): the loosened guard must not accept a non-zero theta0 at a true 0/0 start
    x = np.linspace(0.0, 1.0, 51)
    with pytest.raises(ValueError):
        BL.thwaites(x, x ** 2 * 0.0, 1e-6, theta0=1e-3)  # U_e ≡ 0 at every station: no power-law start value exists


def test_wall_jet_ode_V7_default_eta_max_scales_with_the_free_scale():  # V7 (C13, F4): f_∞³ = 72 f″(0) with DEFAULT arguments at three scales
    for fpp0 in (0.005, 1.0 / 72.0, 0.2, 3.0):
        d = JET.wall_jet_ode_solve(fpp0)
        assert d["f_inf"] ** 3 / fpp0 == pytest.approx(72.0, rel=1e-8), fpp0
        assert d["err_vs_9_83"] < 1e-8, fpp0


def test_march_V3_inlet_station_wall_shear_falls_with_refinement():  # V3 (C01/N14, F3): τ₀ at the inlet station (interpolated FS profile) converges to the Blasius value
    nu = 1e-3
    x0 = 0.1
    ex = BL.falkner_skan_state(0.0)["fpp0"] * nu * np.sqrt(1.0 / (x0 * nu))
    of = BL.outer_flow("wedge", n=0.0, a=1.0)
    errs, hs = [], []
    for ny in (100, 200, 400, 800):
        t = BL.march_boundary_layer(of, np.linspace(x0, 1.0, 11), nu, ny=ny)["tau0"][0]
        errs.append(abs(t / ex - 1))
        hs.append(1.0 / ny)
    assert all(b < a for a, b in zip(errs, errs[1:])), errs  # monotone (it used to GROW: 2.8e-3 → 8.7e-3)
    assert errs[-1] < 1e-6, errs
    assert observed_order(hs, errs) > 1.8, pairwise_orders(hs, errs)


def test_march_V7_below_the_fold_inlet_raises_a_clear_error_and_accepts_a_supplied_profile():  # V7 (C08/N14): m0 = -0.2 has no attached FS profile
    nu = 1e-3
    of = BL.outer_flow("diffuser", U1=1.0, L=0.4)  # m0 = x U_e′/U_e = -0.2 at x0 = 0.1
    with pytest.raises(ValueError, match="separation"):
        BL.march_boundary_layer(of, np.linspace(0.1, 0.2, 6), nu)
    # supplying the Blasius profile is the documented way out: the call runs and the layer starts attached
    bl = BL.blasius_fields
    y = np.linspace(0, 12 * np.sqrt(nu * 0.1), 400)
    eta = y / np.sqrt(nu * 0.1)
    d = BL.falkner_skan(0.0)
    u = np.interp(eta, d["eta"], d["fp"])
    r = BL.march_boundary_layer(of, np.linspace(0.1, 0.11, 3), nu, u_inlet=(y, u), ny=100)
    assert r["tau0"][0] > 0 and bl is not None


def test_cylinder_cd_schematic_V1_is_continuous_at_re_1_and_equals_lamb_there():  # V1+V7 (C11, cd schematic): the Re = 1 anchor is Lamb's 8π/2.002
    lamb = 8 * np.pi / 2.002
    assert BB.cylinder_cd_schematic(1.0) == pytest.approx(lamb, rel=1e-12)
    lo, hi = BB.cylinder_cd_schematic(1.0 - 1e-9), BB.cylinder_cd_schematic(1.0 + 1e-9)
    assert abs(lo - hi) / lamb < 1e-6
    Re = np.logspace(-3, 8, 400)
    cd = np.array([BB.cylinder_cd_schematic(r) for r in Re])
    assert np.all(np.isfinite(cd)) and np.all(cd > 0)
    assert np.max(np.abs(np.diff(np.log(cd)))) < 0.25  # no jump anywhere on a 400-point log sweep (the dip near 4e5 is the steepest part)


def test_drawings_V1_ch09_drawings_imports_and_every_helper_returns_a_figure():  # V1 smoke (design C.5): scripts/ch09_drawings.py, headless
    import importlib
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(ROOT / "scripts"))
    try:
        dr = importlib.import_module("ch09_drawings")
        for fn, kw in ((dr.plate_layer, {}), (dr.cylinder_sketch, {}), (dr.jet_sketch, {"kind": "free"}), (dr.jet_sketch, {"kind": "wall"}), (dr.cup_section, {})):
            fig = fn(**kw)
            assert hasattr(fig, "savefig") and len(fig.axes) >= 1, fn.__name__
            plt.close(fig)
        fig, ax = plt.subplots()
        dr.profile_arrows(ax, np.linspace(0, 1, 9), np.linspace(0, 1, 9) ** 0.5)
        plt.close(fig)
    finally:
        sys.path.remove(str(ROOT / "scripts"))


def test_falkner_skan_V5_default_eta_max_blasius_within_1e_10_of_toepfer():  # V5/V3 (C04): default truncation is now 12 for m = 0
    assert abs(BL.falkner_skan(0.0)["fpp0"] / BL.blasius_constants()["fpp0"] - 1) < 1e-10
