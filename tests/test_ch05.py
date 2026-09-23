"""Verification suite for Chapter 5 — Vorticity Dynamics (Kundu, Cohen & Dowling 5e, §§5.1–5.8).

Evidence levels (``verify-implementation`` skill): V1 analytic · V2 symbolic (sympy / dimensions) · V3 convergence ·
V4 conservation / invariant · V5 published benchmark (``reference/ch05/``, see SOURCES.md) · V6 book value (private,
``tests/book_values_ch05.json``, skipped when absent) · V7 limits / symmetry / invariance. Every test name is
``test_<concept>_<level>_<what>``; the comment on the ``def`` line repeats the level.

A items (CORE, ≥ 2 independent levels): C01 vortex tubes cannot end (5.4) · C02 pressure and stress of the two basic
vortices (5.6)–(5.7) · C03 Kelvin's theorem (5.8)–(5.11) · C04 baroclinic pressure torque (Fig. 5.6, (5.28)) · C05
Helmholtz's theorems · C06 the vorticity equation (5.13) · C07 Biot–Savart (5.14)–(5.16) · C08 filament law (5.17) ·
C09 rotating/baroclinic vorticity equation (5.30) · C10 stretching and tilting (5.31)–(5.32) · C11 absolute
circulation (5.33) and the fluid column · C12 point vortices · C13 images and rings · C14 vortex sheet. Derivations
D02–D16, D18–D20 (every ★★ and ★★★) re-derived with sympy in ``test_*_V2_derivation`` (the ★★★ D09, D10, D11, D15
re-run the design's construction step by step).

Pinned conventions with discrimination tests: ω is the vorticity in (5.1) (tank turns at ω/2) · (5.14) sign +1/(4π)
(the printed −1/(4π) reverses the swirl) · ∇ρ × ∇p (not ∇p × ∇ρ) · sheet strength γ = u_below − u_above (the text;
the caption's u₁ − u₂ is the clockwise sense) · ``ellipk``/``ellipe`` take the parameter m = k² · book e_n points away
from the centre of curvature · (5.30) stretching acts on ω + 2Ω · counterclockwise-positive plane vorticity.

Run: ``.venv/Scripts/python.exe -m pytest tests/test_ch05.py -q -p no:cacheprovider``.
"""
from __future__ import annotations

import inspect
import json
import re
from pathlib import Path

import numpy as np
import pytest
import sympy as sp
from scipy.integrate import dblquad, quad
from scipy.special import ellipe, ellipk, erf

from fluidpy import ch05_vorticity_dynamics as ch05
from fluidpy.core import _stencil as st
from fluidpy.core import biot_savart as BS
from fluidpy.core import curvilinear as CU
from fluidpy.core import integral_theorems as IT
from fluidpy.core import kinematics as K
from fluidpy.core import navier_stokes as NS
from fluidpy.core import vortices as VX
from fluidpy.core import vorticity as VD
from fluidpy.core.units import Q_, dimensional_check
from tools.convergence import observed_order, pairwise_orders

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference" / "ch05"
BOOK = Path(__file__).resolve().parent / "book_values_ch05.json"
needs_ref = pytest.mark.skipif(not (REF / "benchmarks.json").exists(), reason="run reference/ch05/make_refs.py")
book_only = pytest.mark.skipif(not BOOK.exists(), reason="book values are private; see CLAUDE.md rule 9")

RNG = np.random.default_rng(5)
ORDER_TOL = 0.15  # design order ± this (verify-implementation default)
G = 9.81          # the book's (and the chapter's default) g


def book():
    return json.loads(BOOK.read_text(encoding="utf-8"))


def ref_json():
    return json.loads((REF / "benchmarks.json").read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------------------------------------------------
# helpers: sympy expressions → field callables in the core.kinematics layout (components on axis 0, points last)
# ---------------------------------------------------------------------------------------------------------------------
X1, X2, X3, TT = sp.symbols("x1 x2 x3 t", real=True)
XYZ = (X1, X2, X3)


def vec_field(exprs, coords=XYZ):
    """Vector field callable u(x, t) from sympy expressions (x of shape (d,) or (d, N))."""
    f = sp.lambdify((*coords, TT), list(exprs), "numpy")

    def u(x, t=0.0):
        x_ = np.asarray(x, dtype=float)
        shape = x_.shape[1:]
        return np.stack([np.broadcast_to(np.asarray(v, dtype=float), shape) for v in f(*x_, t)])
    return u


def sca_field(expr, coords=XYZ):
    """Scalar field callable F(x, t) from a sympy expression."""
    f = sp.lambdify((*coords, TT), expr, "numpy")

    def F(x, t=0.0):
        x_ = np.asarray(x, dtype=float)
        return np.broadcast_to(np.asarray(f(*x_, t), dtype=float), x_.shape[1:]).copy()
    return F


def s_curl(F, X=XYZ):
    return sp.Matrix([sp.diff(F[2], X[1]) - sp.diff(F[1], X[2]), sp.diff(F[0], X[2]) - sp.diff(F[2], X[0]),
                      sp.diff(F[1], X[0]) - sp.diff(F[0], X[1])])


def s_grad(f, X=XYZ):
    return sp.Matrix([sp.diff(f, v) for v in X])


def s_div(F, X=XYZ):
    return sum(sp.diff(F[i], X[i]) for i in range(3))


def s_dir(a, F, X=XYZ):
    """(a·∇)F."""
    return sp.Matrix([sum(a[j] * sp.diff(F[i], X[j]) for j in range(3)) for i in range(3)])


def s_lap(F, X=XYZ):
    return sp.Matrix([sum(sp.diff(F[i], v, 2) for v in X) for i in range(3)])


def is_zero(M):
    """True when every entry of a sympy Matrix/expression simplifies to 0."""
    M = sp.Matrix(M) if not isinstance(M, sp.Basic) or isinstance(M, sp.MatrixBase) else sp.Matrix([M])
    return all(sp.simplify(sp.expand(e)) == 0 for e in M)


def rel(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.max(np.abs(a - b)) / max(float(np.max(np.abs(b))), 1e-300))


def swirl_field(u_theta):
    """Plane velocity u = u_θ(r) e_θ from a profile callable (r → u_θ)."""
    def u(x, t=0.0):
        X = np.asarray(x, dtype=float)
        r = np.hypot(X[0], X[1])
        ut = np.asarray(u_theta(r), dtype=float)
        return np.stack([-ut * X[1] / r, ut * X[0] / r])
    return u


def lamb_oseen_omega(Gamma, sigma):
    """Exact plane vorticity of the Gaussian (Lamb–Oseen) core, shape (1, …)."""
    return lambda x, t=0.0: np.array([Gamma / (np.pi * sigma ** 2)
                                      * np.exp(-(np.asarray(x)[0] ** 2 + np.asarray(x)[1] ** 2) / sigma ** 2)])


# =====================================================================================================================
# C01 — vortex lines, vortex tubes and their strength; (5.3)–(5.4) tubes cannot end; D01
# =====================================================================================================================
def test_vorticity_field_V1_linear_rotation_gives_twice_b():  # V1 (u = b × x ⇒ ω = 2b, the definition (3.15))
    b = np.array([0.3, -1.2, 0.7])
    u = lambda x, t=0.0: np.cross(b, np.asarray(x, float), axis=0) if np.ndim(x) > 1 else np.cross(b, x)  # noqa
    W = ch05.vorticity_field(u, 1e-3)
    P = RNG.uniform(-2, 2, (3, 20))
    assert np.allclose(W(P), 2 * b[:, None], atol=1e-10)
    # plane fields give ω_z with shape (1, …): the solid body u_θ = ωr/2 of (5.1) has ω_z = ω
    ws = ch05.vorticity_field(swirl_field(lambda r: ch05.solid_body_from_vorticity(r, 3.0)[0]), 1e-4)
    w2 = ws(RNG.uniform(0.2, 1.0, (2, 7)))
    assert w2.shape == (1, 7) and np.allclose(w2, 3.0, rtol=1e-8)


def test_vorticity_field_V3_second_order_on_abc_flow():  # V3 (ABC: ω = u exactly)
    u = ch05.abc_flow(1.0, 0.7, 0.4)
    P = RNG.uniform(-3, 3, (3, 30))
    hs = [0.1, 0.05, 0.025, 0.0125]
    errs = [np.max(np.abs(ch05.vorticity_field(u, h)(P) - u(P))) for h in hs]
    assert abs(observed_order(hs, errs) - 2.0) < ORDER_TOL


def test_vortex_line_V1_abc_vortex_lines_are_streamlines():  # V1 (ω = u ⇒ the curve of (5.3) is the streamline (3.7))
    u = ch05.abc_flow(1.0, 1.0, 1.0)
    W = ch05.vorticity_field(u, 1e-4)  # stencil curl — an independent route to ω
    s, Xw = ch05.vortex_line(W, [0.2, 0.4, -0.1], s_max=2.0, n=201)
    Xs = K.streamline(u, np.array([0.2, 0.4, -0.1]), 0.0, 2.0, True, 201)
    assert Xw.shape == Xs.shape == (3, 201)
    assert np.max(np.abs(Xw - Xs)) < 1e-6
    # the traced curve is tangent to ω: |ω × dX/ds| / |ω| ≈ 0 along it
    dX = np.gradient(Xw, s, axis=1)
    wv = W(Xw)
    tang = np.linalg.norm(np.cross(wv.T, dX.T), axis=1) / np.linalg.norm(wv, axis=0)
    assert np.max(tang[2:-2]) < 1e-3  # np.gradient is only an O(Δs²) check here
    assert s[np.argmin(np.linalg.norm(Xw - np.array([[0.2], [0.4], [-0.1]]), axis=0))] == pytest.approx(0.0, abs=1e-12)


def test_vortex_line_V1_helical_swirl_keeps_zR2_constant():  # V1 (u_φ = aRz: ω = (−aR, 0, 2az), lines zR² = const)
    a = 0.8
    u = ch05.helical_swirl_field(a)
    W = ch05.vorticity_field(u, 1e-4)
    P = RNG.uniform(0.2, 1.0, (3, 10))
    R = np.hypot(P[0], P[1])
    w_exact = np.stack([-a * R * P[0] / R, -a * R * P[1] / R, 2 * a * P[2]])  # ω_R e_R + ω_z e_z in Cartesian
    assert np.allclose(W(P), w_exact, atol=1e-7)
    for x0 in ([0.5, 0.0, 1.0], [0.3, 0.4, 0.6]):
        _, X = ch05.vortex_line(W, x0, s_max=0.6, n=301)
        q = X[2] * (X[0] ** 2 + X[1] ** 2)
        assert np.ptp(q) < 1e-8 * abs(q[0]) + 1e-12  # our integration of (5.3)


def test_vortex_tube_strength_V1_two_routes_on_lamb_oseen():  # V1 (Γ = ∮u·dx = ∫ω·n dA = Γ₀(1 − e⁻¹) at r = σ)
    Gam, nu, t0 = 1.0, 1e-3, 2.5  # σ² = 4ν t₀ = 0.01, σ = 0.1
    u = ch05.lamb_oseen_field(Gam, nu, t0)
    exact = Gam * (1 - np.exp(-1.0))
    lp = IT.planar_loop((0.0, 0.0), None, radius=0.1, n=128)
    errs, hs = [], []
    for nr in (32, 64, 128, 256):
        sf = IT.planar_disc((0.0, 0.0), None, radius=0.1, nr=nr, ntheta=64)
        ts = ch05.vortex_tube_strength(u, lp, sf, omega=lamb_oseen_omega(Gam, 0.1))
        assert isinstance(ts, ch05.TubeStrength)
        assert ts.circulation == pytest.approx(exact, rel=1e-12)  # spectral on a circle
        errs.append(abs(ts.flux - exact))
        hs.append(0.1 / nr)
        assert ts.diff == pytest.approx(ts.circulation - ts.flux, abs=1e-15)
    assert abs(observed_order(hs, errs) - 2.0) < ORDER_TOL  # V3: the midpoint disc of ch02 is second order
    assert errs[-1] < 2e-6
    # flux route from the velocity alone (stencil curl) agrees with the exact-ω route
    sf = IT.planar_disc((0.0, 0.0), None, radius=0.1, nr=256, ntheta=64)
    assert ch05.vortex_tube_strength(u, None, sf).flux == pytest.approx(exact, rel=5e-6)
    # a missing route is NaN; orientation flip (loop clockwise ↔ normal −e_z) changes the sign of both routes
    only = ch05.vortex_tube_strength(u, lp, None)
    assert np.isnan(only.flux) and only.circulation == pytest.approx(exact, rel=1e-12)
    lp3 = IT.planar_loop((0.0, 0.0, 0.0), (0.0, 0.0, -1.0), radius=0.1, n=128)
    sf3 = IT.planar_disc((0.0, 0.0, 0.0), (0.0, 0.0, -1.0), radius=0.1, nr=256, ntheta=64)
    u3 = ch05.planar_field_3d(u)
    om3 = lambda x, t=0.0: np.concatenate([0 * lamb_oseen_omega(Gam, 0.1)(x)] * 2 + [lamb_oseen_omega(Gam, 0.1)(x)])  # noqa
    t3 = ch05.vortex_tube_strength(u3, lp3, sf3, omega=om3)
    assert t3.circulation == pytest.approx(-exact, rel=1e-12) and t3.flux == pytest.approx(-exact, rel=5e-6)


def test_gaussian_tube_V2_flux_function_construction_is_solenoidal():  # V2 (C.3 row 3.4: ∇·ω ≡ 0 by construction)
    R, z, phi = sp.symbols("R z phi", positive=True)
    Gam, a0, L, tw = sp.symbols("Gamma a_0 L tau", positive=True)
    a2 = a0 ** 2 * sp.exp(-z / L)
    chi = Gam / (2 * sp.pi) * (1 - sp.exp(-R ** 2 / a2))
    wz = sp.diff(chi, R) / R
    wR = -sp.diff(chi, z) / R
    wphi = tw * R * wz
    div = sp.diff(R * wR, R) / R + sp.diff(wphi, phi) / R + sp.diff(wz, z)
    assert sp.simplify(div) == 0
    assert sp.simplify(wz - Gam / (sp.pi * a2) * sp.exp(-R ** 2 / a2)) == 0
    assert sp.simplify(wR + Gam / (2 * sp.pi) * R / (a2 * L) * sp.exp(-R ** 2 / a2)) == 0
    # the tube surfaces R/a(z) = const are flux surfaces: ω·∇(R/a) = 0 (ω tangent to the wall)
    s = R / sp.sqrt(a2)
    assert sp.simplify(wR * sp.diff(s, R) + wz * sp.diff(s, z)) == 0
    # the coded values equal these formulas (random points, twisted)
    f = sp.lambdify((R, z, Gam, a0, L, tw), [wR, wphi, wz], "numpy")
    for _ in range(5):
        Rv, zv = RNG.uniform(0.0, 0.3), RNG.uniform(-1, 1)
        got = ch05.gaussian_tube_vorticity(Rv, zv, 1.3, 0.1, 0.8, 4.0)
        assert np.allclose(got, f(Rv, zv, 1.3, 0.1, 0.8, 4.0), rtol=1e-12, atol=1e-14)


def test_gaussian_tube_V1_cartesian_field_and_sections():  # V1 (C.3 rows 3.4–3.5; (5.4) flux independent of z)
    wf = ch05.gaussian_tube_field(1.0, 0.1, 1.0, twist=3.0)
    P = RNG.uniform(-0.2, 0.2, (3, 12))
    R = np.hypot(P[0], P[1])
    wR, wphi, wz = ch05.gaussian_tube_vorticity(R, P[2], 1.0, 0.1, 1.0, 3.0)
    cart = np.stack([wR * P[0] / R - wphi * P[1] / R, wR * P[1] / R + wphi * P[0] / R, wz])
    assert np.allclose(wf(P), cart, rtol=1e-12)
    # the Cartesian field is solenoidal (stencil divergence at round-off × h²), the "broken" one is not
    divs = np.abs(st.div(wf, P, 0.0, 1e-4)) / np.max(np.abs(wf(P)))
    assert np.max(divs) < 1e-6
    bro = ch05.broken_tube_field(1.0, 0.1, 1.0)
    assert np.max(np.abs(st.div(bro, P, 0.0, 1e-4))) > 1.0
    s0, s1 = ch05.gaussian_tube_section(0.0), ch05.gaussian_tube_section(1.0)
    assert (s0["radius"], s0["area"], s0["flux"]) == pytest.approx((0.1, 0.031415927, 0.63212056), rel=1e-7)
    assert (s0["mean_omega"], s0["peak_omega"]) == pytest.approx((20.121022, 31.830989), rel=1e-7)
    assert (s1["radius"], s1["area"], s1["mean_omega"], s1["peak_omega"]) == pytest.approx(
        (0.060653066, 0.011557273, 54.694609, 86.525598), rel=1e-7)
    for z in np.linspace(-1, 3, 9):  # flux constant along the tube; mean vorticity ∝ 1/area
        s = ch05.gaussian_tube_section(z)
        assert s["flux"] == pytest.approx(s0["flux"], rel=1e-14)
        assert s["mean_omega"] * s["area"] == pytest.approx(s0["flux"], rel=1e-12)
        # independent: the flux is ∫₀^{R_t} ω_z 2πR dR by quad
        Rt = s["radius"]
        v, _ = quad(lambda r: ch05.gaussian_tube_vorticity(r, z)[2] * 2 * np.pi * r, 0, Rt, epsabs=0, epsrel=1e-13)
        assert v == pytest.approx(s["flux"], rel=1e-11)
    assert ch05.tube_core_radius(1.0) == pytest.approx(0.1 * np.exp(-0.5), rel=1e-14)


def test_tube_flux_budget_V4_gauss_on_tube_pieces():  # V4 ((5.4): −Γ_lower + side + Γ_upper = ∫∇·ω dV = 0)
    b = ch05.tube_flux_budget("gaussian_tube", 0.0, 1.0, 0.1)
    assert isinstance(b, ch05.TubeFlux)
    assert (b.lower, b.upper) == pytest.approx((-0.6321205588, 0.6321205588), rel=1e-9)
    assert abs(b.side) < 1e-9 and abs(b.total) < 1e-8
    assert b.total == pytest.approx(b.lower + b.side + b.upper, abs=1e-15)
    # other pieces and parameters, twisted too (twist does not change the flux or the wall)
    for (za, zb, R0, kw) in [(-0.5, 2.0, 0.05, {}), (0.0, 1.5, 0.2, dict(twist=6.0)), (1.0, 1.2, 0.1, dict(L=0.4))]:
        bb = ch05.tube_flux_budget("gaussian_tube", za, zb, R0, **kw)
        assert abs(bb.total) < 1e-8 * abs(bb.upper) and abs(bb.side) < 1e-8 * abs(bb.upper)
        assert -bb.lower == pytest.approx(bb.upper, rel=1e-9)
    lo = ch05.tube_flux_budget("lamb_oseen", 0.0, 1.0, 0.1, Gamma=1.0, sigma=0.1)
    assert lo.upper == pytest.approx(1 - np.exp(-1), rel=1e-12) and abs(lo.total) < 1e-14
    # a callable axisymmetric field with its own wall: the same answer as the preset
    cb = ch05.tube_flux_budget(ch05.gaussian_tube_field(), 0.0, 1.0, 0.1, R_of_z=lambda z: 0.1 * np.exp(-np.asarray(z) / 2))
    assert np.allclose(cb, b, atol=1e-12)
    with pytest.raises(ValueError):
        ch05.tube_flux_budget(ch05.gaussian_tube_field(), 0.0, 1.0, 0.1)
    with pytest.raises(ValueError):
        ch05.tube_flux_budget("nope")


def test_tube_flux_budget_V1_broken_field_fails_exactly_as_computed():  # V1 (∇·ω ≠ 0 ⇒ (5.4) fails; our number)
    b = ch05.tube_flux_budget("broken", 0.0, 1.0, 0.1)
    flux0 = 1 - np.exp(-1.0)
    assert b.lower == pytest.approx(-flux0, rel=1e-9) and b.upper == pytest.approx(flux0 * np.exp(-1.0), rel=1e-9)
    assert b.total == pytest.approx(-0.399577, rel=1e-5) and b.side == 0.0
    # = ∫_V ∇·ω dV = ∫_0^1 ∂ω_z/∂z(integrated over the disc) dz = flux0 (e^{−1} − 1)  (independent)
    assert b.total == pytest.approx(flux0 * (np.exp(-1.0) - 1.0), rel=1e-9)


def test_vortex_tube_V1_traced_lines_stay_on_the_tube_wall():  # V1 (lines lie on R/a(z) = const; Fig. 5.1)
    tube = ch05.vortex_tube("gaussian_tube", (0, 0, 0), (0, 0, 1), 0.1, n_lines=12, s_max=0.8, n=41, twist=5.0)
    assert isinstance(tube, ch05.Tube) and tube.array.shape == (12, 3, 41) and tube.s.shape == (41,)
    L = tube.array
    R = np.hypot(L[:, 0], L[:, 1])
    ratio = R / ch05.tube_core_radius(L[:, 2])
    assert np.ptp(ratio) < 1e-8 and ratio.mean() == pytest.approx(1.0, rel=1e-8)
    # twist makes helices: the azimuth advances along each line; the tube narrows downstream (z up)
    ang = np.unwrap(np.arctan2(L[0, 1], L[0, 0]))
    assert np.ptp(ang) > 0.1
    assert R[:, -1].mean() < R[:, 0].mean()
    with pytest.raises(ValueError):  # a seed circle reaching irrotational fluid stops a line early
        ch05.vortex_tube(lambda x, t=0.0: 0.0 * np.asarray(x), (0, 0, 0), (0, 0, 1), 0.1, n_lines=4)


def test_tube_flux_budget_traced_V3_closed_polyhedron_converges():  # V3/V4 (non-axisymmetric route of (5.4))
    sides, totals, uppers = [], [], []
    for nl in (8, 16, 32):
        tube = ch05.vortex_tube("gaussian_tube", (0, 0, 0), (0, 0, 1), 0.1, n_lines=nl, s_max=0.8, n=41, twist=5.0)
        d = ch05.tube_flux_budget_traced("gaussian_tube", tube, twist=5.0)
        f = d["flux"]
        assert abs(f.total) < 5e-9  # a closed surface: Gauss for a solenoidal field
        totals.append(abs(f.total))
        uppers.append(f.upper)
        assert d["mean_omega_upper"] > d["mean_omega_lower"] and d["area_upper"] < d["area_lower"]
    # the end fluxes approach the exact tube strength Γ(1 − e⁻¹) as the polygon → circle (order 2 in 1/n_lines)
    err = [abs(v - (1 - np.exp(-1))) for v in uppers]
    assert abs(observed_order([1 / 8, 1 / 16, 1 / 32], err) - 2.0) < 0.25
    bro = ch05.vortex_tube("gaussian_tube", (0, 0, 0), (0, 0, 1), 0.1, n_lines=16, s_max=0.8, n=41)
    assert abs(ch05.tube_flux_budget_traced("broken", bro)["flux"].total) > 0.1  # the broken field fails


def test_vortex_ring_vorticity_V1_meridional_flux_is_Gamma():  # V1 (a closed tube: the ring's strength)
    v, _ = dblquad(lambda z, R: ch05.vortex_ring_vorticity(R, z, 1.7, 0.5, 0.05), 0.2, 0.8, -0.3, 0.3,
                   epsabs=0, epsrel=1e-11)
    assert v == pytest.approx(1.7, rel=1e-9)


# =====================================================================================================================
# C02 — pressure and viscous stress in the two basic vortices (5.1)–(5.7); R02–R08, N05–N09; D02, D03
# =====================================================================================================================
def test_solid_body_V1_omega_is_the_vorticity_not_the_rate():  # V1 ((5.1) u_θ = ωr/2; WV ω passed as rate)
    r = np.linspace(0.0, 2.0, 9)
    u, wz = ch05.solid_body_from_vorticity(r, 4.0)
    assert np.allclose(u, 2.0 * r, rtol=1e-15) and np.allclose(wz, 4.0)
    W = ch05.vorticity_field(swirl_field(lambda q: ch05.solid_body_from_vorticity(q, 4.0)[0]), 1e-4)
    assert np.allclose(W(RNG.uniform(0.3, 1.0, (2, 6))), 4.0, rtol=1e-9)
    # wrong variant: using ω as ch03's rotation rate doubles the vorticity
    Wbad = ch05.vorticity_field(swirl_field(lambda q: VX.solid_body_rotation(q, 4.0)), 1e-4)
    assert np.allclose(Wbad(RNG.uniform(0.3, 1.0, (2, 6))), 8.0, rtol=1e-9)


def test_line_vortex_V1_circulation_on_every_circle_and_irrotational():  # V1 ((5.2); Γ = 2πB)
    r = np.array([0.0, 0.1, 1.0, 7.0])
    u, wz = ch05.line_vortex_gamma(r, 2.5)
    assert np.isnan(wz[0]) and np.all(wz[1:] == 0.0)
    assert np.allclose(u[1:], 2.5 / (2 * np.pi * r[1:]), rtol=1e-15)
    uf = swirl_field(lambda q: ch05.line_vortex_gamma(q, 2.5)[0])
    for rr, c in ((0.3, (0.0, 0.0)), (2.0, (0.5, -0.3))):
        pts = ch05.circle_loop_points(c, rr, 256)
        assert ch05.loop_circulation(uf, pts) == pytest.approx(2.5, rel=1e-12)  # both enclose the axis
    assert abs(ch05.loop_circulation(uf, ch05.circle_loop_points((2.0, 0.0), 0.5, 256))) < 1e-12  # does not
    assert np.allclose(ch05.vorticity_field(uf, 1e-4)(RNG.uniform(0.5, 1.5, (2, 5))), 0.0, atol=1e-7)


def test_solid_body_pressure_V2_derivation():  # V2 — D02 re-derived: cylindrical Euler → (5.5a, b) → (5.6)
    Rs, Ph, Zs = CU.coordinates("cylindrical")
    om, rho, g, po, zz = sp.symbols("omega rho g p_o z_c", positive=True)
    u = [0, om * Rs / 2, 0]  # (5.1): ω is the vorticity
    acc = CU.advective_acceleration(u, "cylindrical")  # steady: Du/Dt = (u·∇)u in physical components
    assert sp.simplify(acc[0] + (om * Rs / 2) ** 2 / Rs) == 0 and acc[1] == 0 and acc[2] == 0  # steps 3–4: −u_θ²/r e_r
    # steps 5–9: ∂p/∂r = ρu_θ²/r = ρω²r/4, integrate in r with an unknown f(z), fix it by ∂p/∂z = −ρg
    f = sp.Function("f")
    p_trial = sp.integrate(rho * om ** 2 * Rs / 4, Rs) + f(Zs)
    fsol = sp.dsolve(sp.Eq(sp.diff(p_trial, Zs), -rho * g), f(Zs), ics={f(0): po}).rhs
    p = p_trial.subs(f(Zs), fsol)
    assert sp.simplify(p - (po + rho * om ** 2 * Rs ** 2 / 8 - rho * g * Zs)) == 0  # (5.6)
    # the result satisfies both balances (5.5a), (5.5b) — and the trap (ω as the rate) is 4× too big
    assert sp.simplify(-rho * (om * Rs / 2) ** 2 / Rs + sp.diff(p, Rs)) == 0
    assert sp.simplify(sp.diff(p, Zs) + rho * g) == 0
    assert sp.simplify(rho * (om * Rs) ** 2 / Rs / (rho * (om * Rs / 2) ** 2 / Rs)) == 4
    # step 10: an isobar p = p_o + Δp is z = ω²r²/8g − Δp/ρg
    dp = sp.Symbol("Delta_p")
    zis = sp.solve(sp.Eq(p, po + dp), Zs)[0]
    assert sp.simplify(zis - (om ** 2 * Rs ** 2 / (8 * g) - dp / (rho * g))) == 0
    # line vortex (5.7) by the same moves (the demoted a-D07): ∫ρΓ²/(4π²r³) dr → −ρΓ²/(8π²r²), p → p_∞ as r → ∞
    Gm, pinf = sp.symbols("Gamma p_inf", positive=True)
    ut = Gm / (2 * sp.pi * Rs)
    pl = pinf - sp.integrate(rho * ut ** 2 / Rs, (Rs, Rs, sp.oo)) - rho * g * Zs
    assert sp.simplify(pl - (pinf - rho * Gm ** 2 / (8 * sp.pi ** 2 * Rs ** 2) - rho * g * Zs)) == 0
    # the Bernoulli function: grows as ω²r²/4 across the tank (R07), uniform for the line vortex
    assert sp.simplify(((om * Rs / 2) ** 2 / 2 + g * Zs + p / rho) - (po / rho + om ** 2 * Rs ** 2 / 4)) == 0
    assert sp.simplify((ut ** 2 / 2 + g * Zs + pl / rho) - pinf / rho) == 0
    # the code equals the derived forms
    fp = sp.lambdify((Rs, Zs, om, rho, g, po), p)
    fl = sp.lambdify((Rs, Zs, Gm, rho, g, pinf), pl)
    for _ in range(5):
        rv, zv = RNG.uniform(0.01, 1.0), RNG.uniform(-1, 1)
        assert ch05.solid_body_pressure(rv, zv, 7.0, 998.0, 9.8, 1e5) == pytest.approx(fp(rv, zv, 7.0, 998.0, 9.8, 1e5),
                                                                                         rel=1e-14)
        assert ch05.line_vortex_pressure(rv, zv, 1.3, 998.0, 9.8, 1e5) == pytest.approx(fl(rv, zv, 1.3, 998.0, 9.8, 1e5),
                                                                                          rel=1e-14)


def test_solid_body_pressure_V1_gradients_isobars_contract():  # V1 ((5.5a, b), (5.6), Fig. 5.2)
    assert ch05.solid_body_pressure(0.1, 0.0, 10.0) == pytest.approx(125.0, rel=1e-14)
    r, z = RNG.uniform(0.05, 0.5, 8), RNG.uniform(-0.3, 0.3, 8)
    dpdr, dpdz = ch05.solid_body_pressure_gradients(r, z, 6.0, 1000.0, G)
    h = 1e-6
    assert np.allclose(dpdr, (ch05.solid_body_pressure(r + h, z, 6.0) - ch05.solid_body_pressure(r - h, z, 6.0)) / (2 * h),
                       rtol=1e-7)
    assert np.allclose(dpdz, -1000.0 * G) and np.allclose(dpdr, 1000.0 * 36.0 * r / 4, rtol=1e-14)
    # isobars: p(r, z_iso(r)) = p_o + Δp at every r; rim–centre height ω²R²/8g = 12.742 mm (ω = 10 s⁻¹, R = 0.1 m)
    for dpg in (0.0, 0.05, -0.02):
        zi = ch05.isobar_height(r, 10.0, dpg, "solid")
        assert np.allclose(ch05.solid_body_pressure(r, zi, 10.0) / (1000.0 * G), dpg, atol=1e-13)
    assert ch05.isobar_height(0.1, 10.0) - ch05.isobar_height(0.0, 10.0) == pytest.approx(0.0127421, rel=1e-5)
    assert ch05.isobar_height(0.1, 1.0, kind="line") == pytest.approx(-0.129104, rel=1e-5)
    assert ch05.line_vortex_pressure(0.1, 0.0, 1.0) == pytest.approx(-1266.5148, rel=1e-7)
    with pytest.raises(ValueError):
        ch05.isobar_height(0.1, 1.0, kind="rankine")
    with pytest.raises(ValueError):
        ch05.isobar_height(0.1, 1.0, kind="nope")


def test_rotating_tank_V1_volume_conserving_paraboloid():  # V1 (Fig. 5.2 free surface; independent quad volume)
    t = ch05.rotating_tank_free_surface(0.1, 0.1, 5.0)
    assert (t["z_vertex"], t["z_rim"]) == pytest.approx((0.0936290, 0.1063710), rel=1e-6)
    rise = 5.0 ** 2 * 0.1 ** 2 / (4 * G)  # Ω²R²/4g, both ways
    assert t["rim_rise"] == pytest.approx(rise, rel=1e-10) and t["centre_drop"] == pytest.approx(rise, rel=1e-10)
    assert t["r_dry"] == 0.0 and not t["spills"] and np.isnan(t["r_top"])
    # the surface is the isobar of (5.6) with ω = 2Ω_t: height difference = ω²R²/8g
    assert t["z_rim"] - t["z_vertex"] == pytest.approx(ch05.isobar_height(0.1, 10.0) - ch05.isobar_height(0.0, 10.0),
                                                       rel=1e-12)
    # dry bottom: fast spin — the water volume ∫ max(z_s, 0) 2πr dr is still πR²·depth (our quad, not the code's)
    for (R, d, Om, H) in [(0.1, 0.05, 20.0, None), (1.0, 3.0, 40.0, None), (0.3, 0.2, 3.0, 0.21)]:
        s = ch05.rotating_tank_free_surface(R, d, Om, H=H)
        zs = lambda r: s["z_vertex"] + Om ** 2 * r ** 2 / (2 * G)  # noqa: E731
        vol, _ = quad(lambda r: max(zs(r), 0.0) * 2 * np.pi * r, 0, R, epsabs=0, epsrel=1e-12, limit=200,
                      points=[s["r_dry"]] if 0 < s["r_dry"] < R else None)
        assert vol == pytest.approx(np.pi * R ** 2 * d, rel=1e-9)
        assert abs(s["volume_residual"]) < 1e-12
        if s["r_dry"] > 0:
            assert zs(s["r_dry"]) == pytest.approx(0.0, abs=1e-10) and s["uncovered_area"] == pytest.approx(
                np.pi * s["r_dry"] ** 2, rel=1e-14)
        if H is not None:
            assert s["spills"] == (s["z_rim"] > H)
    # a closed tank: the lid clips the surface; volume still conserved (independent)
    c = ch05.rotating_tank_free_surface(1.0, 3.0, 40.0, H=4.0, closed=True)
    zsc = lambda r: min(max(c["z_vertex"] + 1600.0 * r ** 2 / (2 * G), 0.0), 4.0)  # noqa: E731
    vol, _ = quad(lambda r: zsc(r) * 2 * np.pi * r, 0, 1.0, points=[c["r_dry"], c["r_top"]], epsabs=0, epsrel=1e-12,
                  limit=200)
    assert vol == pytest.approx(3.0 * np.pi, rel=1e-9)
    assert c["r_dry"] < c["r_top"] < 1.0


def test_line_vortex_pressure_V1_radial_balance_and_funnel():  # V1 ((5.7) satisfies (5.5a) with (5.2))
    r = np.linspace(0.05, 2.0, 30)
    h = 1e-7
    dpdr = (ch05.line_vortex_pressure(r + h, 0.0, 2.0) - ch05.line_vortex_pressure(r - h, 0.0, 2.0)) / (2 * h)
    ut = ch05.line_vortex_gamma(r, 2.0)[0]
    assert np.allclose(dpdr, 1000.0 * ut ** 2 / r, rtol=1e-6)
    # the funnel: (c − z) r² = const (a cubic surface, not a quadric) on every isobar
    for dpg in (-0.1, -0.5):
        zi = ch05.isobar_height(r, 2.0, dpg, "line")
        assert np.ptp((-dpg - zi) * r ** 2) < 1e-12
    assert ch05.bernoulli_across_vortex("line", r) == pytest.approx(0.0, abs=1e-12)


def test_rankine_pressure_V1_matched_core_and_tornado_contract():  # V1 (R08, Exercise 5.2 geometry; ch04 form)
    Gam, a = 1.0, 0.1
    assert ch05.rankine_pressure(0.1, 0.0, Gam, a) == pytest.approx(-1266.5148, rel=1e-7)
    assert ch05.rankine_pressure(0.0, 0.0, Gam, a) == pytest.approx(-2533.0296, rel=1e-7)
    eps = 1e-9
    pin, pout = ch05.rankine_pressure(a - eps, 0.0, Gam, a), ch05.rankine_pressure(a + eps, 0.0, Gam, a)
    assert pin == pytest.approx(pout, rel=1e-7)  # continuous
    h = 1e-6
    sl_in = (ch05.rankine_pressure(a - 2 * h, 0, Gam, a) - ch05.rankine_pressure(a - 4 * h, 0, Gam, a)) / (2 * h)
    sl_out = (ch05.rankine_pressure(a + 4 * h, 0, Gam, a) - ch05.rankine_pressure(a + 2 * h, 0, Gam, a)) / (2 * h)
    assert sl_in == pytest.approx(sl_out, rel=1e-3)  # continuous slope (both ≈ ρu_θ²/r at a)
    r = np.linspace(0.0, 0.4, 41)
    from fluidpy.core.bernoulli import rankine_vortex_pressure  # ch04's independent implementation
    assert np.allclose(ch05.rankine_pressure(r, 0.0, Gam, a), rankine_vortex_pressure(r, Gam, a), rtol=1e-12)
    # radial balance everywhere (inside: solid body with ω_c = Γ/πa²)
    rr = np.array([0.02, 0.07, 0.15, 0.3])
    dp = (ch05.rankine_pressure(rr + h, 0, Gam, a) - ch05.rankine_pressure(rr - h, 0, Gam, a)) / (2 * h)
    ut = VX.rankine_vortex(rr, Gam, a)[0]
    assert np.allclose(dp, 1000.0 * ut ** 2 / rr, rtol=1e-6)
    # tornado preset (air ρ = 1.2, Γ = 1e4 m²/s, a = 50 m): u_max, core-edge and axis deficits
    tor = ch05.vortex_pressure_scenario("rankine", 50.0, rho=1.2, Gamma=1e4, a=50.0)
    assert tor["u_theta"] == pytest.approx(31.830989, rel=1e-7) and tor["p"] == pytest.approx(-607.92710, rel=1e-7)
    assert ch05.rankine_pressure(0.0, 0.0, 1e4, 50.0, rho=1.2) == pytest.approx(-1215.8542, rel=1e-7)
    # the inverse helpers: Γ from the edge deficit, radius from a pressure
    assert ch05.tornado_circulation(-607.92710185, 50.0, rho=1.2) == pytest.approx(1e4, rel=1e-9)
    for pg in (-100.0, -500.0, -900.0):
        rq = ch05.rankine_radius_at_pressure(pg, 1e4, 50.0, rho=1.2)
        assert ch05.rankine_pressure(rq, 0.0, 1e4, 50.0, rho=1.2) == pytest.approx(pg, rel=1e-10)
    assert np.isnan(ch05.rankine_radius_at_pressure(-5000.0, 1e4, 50.0, rho=1.2))  # below the axis deficit
    # the rankine isobar kind equals the pressure/ρg route
    assert ch05.isobar_height(0.05, Gam, 0.0, "rankine", a=a) == pytest.approx(
        ch05.rankine_pressure(0.05, 0.0, Gam, a) / (1000.0 * G), rel=1e-13)


def test_bernoulli_across_vortex_V1_rotational_vs_irrotational():  # V1 (R07: B grows as ω²r²/4; uniform for (5.2))
    r = np.linspace(0.0, 1.0, 11)
    assert np.allclose(ch05.bernoulli_across_vortex("solid", r, param=6.0), 36.0 * r ** 2 / 4, rtol=1e-12, atol=1e-13)
    assert ch05.bernoulli_across_vortex("solid", 0.1, param=10.0) == pytest.approx(0.25, rel=1e-13)
    assert ch05.bernoulli_across_vortex("solid", 0.5, z=0.3, param=6.0, r_ref=0.2) == pytest.approx(
        36 * (0.25 - 0.04) / 4, rel=1e-12)
    assert np.allclose(ch05.bernoulli_across_vortex("line", r[1:], param=3.0, r_ref=0.4), 0.0, atol=1e-11)
    # Rankine (ch04 value): B(0) − B(∞) = −Γ²/(4π²a²) = −1 for Γ = 2π, a = 1; flat outside the core
    assert ch05.bernoulli_across_vortex("rankine", 0.0, param=2 * np.pi, a=1.0, r_ref=np.inf) == pytest.approx(-1.0,
                                                                                                               rel=1e-12)
    out = ch05.bernoulli_across_vortex("rankine", np.linspace(1.0, 5.0, 9), param=2 * np.pi, a=1.0, r_ref=np.inf)
    assert np.allclose(out, 0.0, atol=1e-12)
    with pytest.raises(ValueError):
        ch05.bernoulli_across_vortex("rankine", 0.5)
    with pytest.raises(ValueError):
        ch05.bernoulli_across_vortex("nope", 0.5)


def test_viscous_stress_V2_derivation():  # V2 — D03: σ_rθ = −μΓ/πr² ≠ 0, net viscous force 0 by three routes
    r, th, mu, Gm = sp.symbols("r theta mu Gamma", positive=True)
    d = ch05.polar_viscous_stress(0, Gm / (2 * sp.pi * r), r, th, mu)
    assert sp.simplify(d["sigma_rtheta"] + mu * Gm / (sp.pi * r ** 2)) == 0  # steps 1–4
    assert d["sigma_rr"] == 0 and d["sigma_thetatheta"] == 0  # step 5
    sig = d["sigma_rtheta"]
    assert sp.simplify(sp.diff(r ** 2 * sig, r) / r ** 2) == 0  # steps 6–7 (metric-correct divergence)
    assert sp.simplify(sp.diff(sig, r)) != 0  # the trap: the naive dσ/dr is not zero
    ut = Gm / (2 * sp.pi * r)
    assert sp.simplify(sp.diff(ut, r, 2) + sp.diff(ut, r) / r - ut / r ** 2) == 0  # step 8
    Rs, Ph, Zs = CU.coordinates("cylindrical")
    vl = CU.vector_laplacian([0, Gm / (2 * sp.pi * Rs), 0], "cylindrical")
    assert all(sp.simplify(c) == 0 for c in vl)  # step 9 (Appendix B vector Laplacian, off the axis)
    # contrast: the Gaussian (Lamb–Oseen) vortex has stress AND a net force μ(∇²u)_θ = μ ∂ω_z/∂r  ((4.40))
    s = sp.Symbol("sigma_c", positive=True)
    ug = Gm / (2 * sp.pi * r) * (1 - sp.exp(-r ** 2 / s ** 2))
    wg = sp.diff(r * ug, r) / r
    dg = ch05.polar_viscous_stress(0, ug, r, th, mu)
    fdiv = sp.diff(r ** 2 * dg["sigma_rtheta"], r) / r ** 2
    assert sp.simplify(fdiv - mu * sp.diff(wg, r)) == 0 and sp.simplify(fdiv) != 0
    # the code: every route of vortex_stress_force agrees with these closed forms
    for rv in (0.03, 0.1, 0.25):
        line = ch05.vortex_stress_force("line", rv, mu=1e-3, Gamma=1.0)
        assert line["sigma_rtheta"] == pytest.approx(-1e-3 / (np.pi * rv ** 2), rel=1e-7)
        assert max(abs(line[k]) for k in ("force_divergence", "force_laplacian", "force_curl")) < 1e-6 * abs(
            line["sigma_rtheta"]) / rv
        gau = ch05.vortex_stress_force("gaussian", rv, mu=1e-3, Gamma=1.0, sigma_c=0.1)
        fg = float(sp.lambdify((r, mu, Gm, s), fdiv)(rv, 1e-3, 1.0, 0.1))
        for k in ("force_divergence", "force_laplacian", "force_curl"):
            assert gau[k] == pytest.approx(fg, rel=1e-5, abs=1e-12)
        sol = ch05.vortex_stress_force("solid", rv, mu=1e-3, omega=2.0)
        assert abs(sol["sigma_rtheta"]) < 1e-12 and max(abs(sol[k]) for k in ("force_divergence", "force_laplacian",
                                                                                 "force_curl")) < 1e-9
    with pytest.raises(ValueError):
        ch05.vortex_stress_force("nope", 0.1)


def test_viscous_force_V1_naive_divergence_misses_the_metric():  # V1 (discrimination of D03's trap)
    d = ch05.polar_net_viscous_force(lambda q: 1.0 / (2 * np.pi * np.asarray(q)), 0.1, mu=1e-3)
    scale = abs(d["sigma_rtheta"]) / 0.1  # the force scale σ/r of the element [N/m³]
    assert abs(d["force_metric"]) < 1e-6 * scale and abs(d["force_laplacian"]) < 1e-6 * scale
    # dσ/dr = 2μΓ/πr³ ≠ 0 — dropping the r² of the polar divergence invents a force
    assert d["force_naive"] == pytest.approx(2e-3 / (np.pi * 0.1 ** 3), rel=1e-6)


def test_stress_vs_net_force_table_V1_three_rows():  # V1 (N09's table)
    rows = ch05.stress_vs_net_force_table()
    names = [r["name"] for r in rows]
    assert names == ["solid body", "line vortex", "Lamb–Oseen"]
    sol, line, lo = rows
    assert abs(sol["sigma_rtheta"]) < 1e-15 and abs(sol["net_force"]) < 1e-12
    assert line["sigma_rtheta"] == pytest.approx(-1e-3 * 0.01 / (np.pi * 0.05 ** 2), rel=1e-7)
    assert abs(line["net_force"]) < 1e-6 * abs(line["sigma_rtheta"]) / 0.05  # nested-difference round-off
    assert abs(lo["sigma_rtheta"]) > 0 and abs(lo["net_force"]) > 1e-3 * abs(lo["sigma_rtheta"]) / 0.05
    # the Lamb–Oseen row's force equals μ∂ω_z/∂r (4.40) at r = 0.05, σ_c² = 4νt
    s2 = 4 * 1e-6 * 100.0
    dw = 0.01 / (np.pi * s2) * np.exp(-0.05 ** 2 / s2) * (-2 * 0.05 / s2)
    assert lo["net_force"] == pytest.approx(1e-3 * dw, rel=1e-5)
    assert lo["omega_z"] == pytest.approx(0.01 / (np.pi * s2) * np.exp(-0.05 ** 2 / s2), rel=1e-12)


def test_rotating_cylinder_V1_no_slip_torque_and_energy():  # V1/V4 (R08, N08: torque −2μΓ at every r; energy)
    a, om = 0.1, 2.0
    u, wz = ch05.rotating_cylinder_flow(np.array([0.05, a, 0.3]), a, om)
    assert u[1] == pytest.approx(om * a / 2, rel=1e-14)  # no slip on the cylinder turning at ω/2
    assert u[2] == pytest.approx(om * a ** 2 / (2 * 0.3), rel=1e-14) and wz[0] == om and wz[2] == 0.0
    Gam = np.pi * a ** 2 * om
    assert ch05.loop_circulation(swirl_field(lambda q: ch05.rotating_cylinder_flow(q, a, om)[0]),
                                 ch05.circle_loop_points((0, 0), 0.4, 256)) == pytest.approx(Gam, rel=1e-12)
    rr = np.geomspace(0.1, 100.0, 7)
    assert np.allclose(ch05.torque_per_length(rr, 1.0, 1e-3), -2e-3, rtol=1e-13)  # independent of r
    d = ch05.dissipation_outside_cylinder(0.1, 1.0, 1.0, 1e-3)
    assert (d["dissipated"], d["power_in"], d["power_out"]) == pytest.approx((0.0315127, 0.0318310, 3.18310e-4), rel=1e-5)
    # independent: ρε = μ(r d(u_θ/r)/dr)² in polar form for u_θ = Γ/2πr, integrated by our own quad
    phi = lambda r: 1e-3 * (-1.0 / (np.pi * r ** 2)) ** 2  # noqa: E731
    v, _ = quad(lambda r: phi(r) * 2 * np.pi * r, 0.1, 1.0, epsabs=0, epsrel=1e-13)
    assert d["dissipated"] == pytest.approx(v, rel=1e-10)
    assert d["dissipated"] == pytest.approx(d["power_in"] - d["power_out"], rel=1e-10) and abs(d["residual"]) < 1e-15
    # power = |torque| × angular speed u_θ/r; as R_out → ∞ everything put in is dissipated
    assert d["power_in"] == pytest.approx(2e-3 * (1.0 / (2 * np.pi * 0.1)) / 0.1, rel=1e-13)
    far = ch05.dissipation_outside_cylinder(0.1, 1e3, 1.0, 1e-3)
    assert far["dissipated"] == pytest.approx(far["power_in"], rel=1e-7)


@pytest.mark.parametrize("kind,kw", [("solid", dict(omega=2.0)), ("line", dict(Gamma=1.0)),
                                     ("rankine", dict(Gamma=1.0, a=0.1)), ("cylinder", dict(omega=2.0, a=0.1))])
def test_vortex_pressure_scenario_V1_full_balance_every_kind(kind, kw):  # V1 (E2: every field, every radius)
    mu, rho = 1e-3, 1000.0
    a = kw.get("a", 0.0)
    Gam = kw.get("Gamma", np.pi * a ** 2 * kw.get("omega", 0.0))
    for r in (0.03, 0.07, 0.15, 0.4, 1.2):
        d = ch05.vortex_pressure_scenario(kind, r, rho=rho, mu=mu, **kw)
        # (5.5a): dp/dr/ρ = u_θ²/r
        assert abs(d["radial_residual"]) <= 1e-6 * max(d["centripetal"], 1e-12)
        assert d["centripetal"] == pytest.approx(d["u_theta"] ** 2 / r, rel=1e-14)
        rigid = kind == "solid" or r < a
        if kind == "solid":
            assert d["u_theta"] == pytest.approx(kw["omega"] * r / 2, rel=1e-14)
            assert d["p"] == pytest.approx(ch05.solid_body_pressure(r, 0.0, kw["omega"], rho), rel=1e-13)
            assert d["B"] == pytest.approx(kw["omega"] ** 2 * r ** 2 / 4, rel=1e-12)
            assert d["omega_z"] == kw["omega"]
        elif r >= a:
            assert d["u_theta"] == pytest.approx(Gam / (2 * np.pi * r), rel=1e-13)
            assert d["p"] == pytest.approx(ch05.line_vortex_pressure(r, 0.0, Gam, rho), rel=1e-12)
            assert abs(d["B"]) < 1e-9 and d["omega_z"] == 0.0
        if rigid:
            assert d["sigma_rtheta"] == 0.0 and d["net_viscous_force"] == 0.0 and d["torque"] == 0.0
        else:
            assert d["sigma_rtheta"] == pytest.approx(-mu * Gam / (np.pi * r ** 2), rel=1e-6)
            assert abs(d["net_viscous_force"]) < 1e-6 * abs(d["sigma_rtheta"]) / r
            assert d["torque"] == pytest.approx(-2 * mu * Gam, rel=1e-6)  # N08: the same at every radius
            assert "no net viscous force" in d["status"]
        assert d["in_fluid"] == (not (kind == "cylinder" and r < a))
        # the surface height through the reference point is the isobar of (5.6)/(5.7)
        zs = ch05.isobar_height(r, kw["omega"] if kind == "solid" else Gam, 0.0, "solid" if kind == "solid" else (
            "line" if kind == "line" else "rankine"), a=a if a > 0 else None)
        assert d["surface_z"] == pytest.approx(zs, rel=1e-13)
    with pytest.raises(ValueError):
        ch05.vortex_pressure_scenario("nope", 0.1)


@pytest.mark.parametrize("kind,kw", [("rankine", dict(Gamma=1.0, a=0.1)), ("cylinder", dict(omega=2.0, a=0.1))])
def test_vortex_pressure_scenario_V1_core_edge_is_consistent(kind, kw):  # V1 (E2 at r = a: core edge / cylinder wall)
    """At r = a the vorticity jumps (Rankine core edge; the cylinder's surface). Just outside, the fluid's stress is
    σ_rθ = −μΓ/πa² and the torque −2μΓ (N08, ``torque_per_length``); the net viscous force is concentrated on the edge
    and has no finite pointwise value. The scenario must not report a step-dependent finite force together with the
    status "no net viscous force", and it must report the fluid's wall stress and torque."""
    mu, a = 1e-3, 0.1
    Gam = kw.get("Gamma", np.pi * a ** 2 * kw.get("omega", 0.0))
    d = ch05.vortex_pressure_scenario(kind, a, mu=mu, **kw)
    just_out = ch05.vortex_pressure_scenario(kind, a * (1 + 1e-3), mu=mu, **kw)
    assert just_out["torque"] == pytest.approx(-2 * mu * Gam, rel=1e-5)  # passes: 0.1 % outside the edge
    assert d["torque"] == pytest.approx(ch05.torque_per_length(a, Gam, mu), rel=1e-3), \
        f"torque at r = a: {d['torque']:.6g} vs −2μΓ = {-2 * mu * Gam:.6g} (σ_rθ = {d['sigma_rtheta']:.6g})"
    if "no net viscous force" in d["status"]:
        assert abs(d["net_viscous_force"]) < 1e-6 * abs(just_out["sigma_rtheta"]) / a, \
            f"net force {d['net_viscous_force']:.6g} N/m³ reported with status {d['status']!r}"
    # the stress on each side, computed independently (σ_rθ = μ r d(u_θ/r)/dr of the core and of the outer profile)
    Omc = Gam / (2 * np.pi * a ** 2)  # core rotation rate u_θ/r (rigid)
    sig_in = mu * a * 0.0 * Omc  # u_θ/r = const inside ⇒ σ_rθ(a⁻) = 0
    sig_out = mu * a * (-2 * Gam / (2 * np.pi * a ** 3))  # d(Γ/2πr²)/dr at a⁺ ⇒ −μΓ/πa²
    assert d["sigma_rtheta"] == pytest.approx(sig_out, rel=1e-7)
    assert d["edge_line_force"] == pytest.approx(sig_out - sig_in, rel=1e-7)  # the jump σ(a⁺) − σ(a⁻)
    assert np.isnan(d["net_viscous_force"]) and "no net viscous force" not in d["status"]
    inside = ch05.vortex_pressure_scenario(kind, a * (1 - 1e-3), mu=mu, **kw)
    assert inside["sigma_rtheta"] == 0.0 and inside["edge_line_force"] == 0.0
    assert just_out["edge_line_force"] == 0.0 and abs(just_out["net_viscous_force"]) < 1e-6 * abs(sig_out) / a
    # the edge force per unit length 2πa × jump balances the torque carried outward: 2πa·(jump)·a = −2μΓ
    assert 2 * np.pi * a * d["edge_line_force"] * a == pytest.approx(-2 * mu * Gam, rel=1e-7)
    # the reported jump is independent of the step (it was 1/h-dependent before the fix)
    assert ch05.vortex_pressure_scenario(kind, a, mu=mu, **kw)["edge_line_force"] == d["edge_line_force"]
    assert ch05.vortex_pressure_scenario("line", 0.1, mu=mu, Gamma=1.0)["edge_line_force"] == 0.0


def test_vortex_pressure_V2_dimensions():  # V2 (pint: (5.6), (5.7), σ_rθ, torque, dissipation power)
    p6 = dimensional_check(lambda rho, omega, r, g, z: rho * omega ** 2 * r ** 2 / 8 - rho * g * z, "[pressure]",
                           rho=Q_(1000, "kg/m**3"), omega=Q_(10, "1/s"), r=Q_(0.1, "m"), g=Q_(G, "m/s**2"),
                           z=Q_(0.0, "m"))
    assert p6.to("Pa").magnitude == pytest.approx(ch05.solid_body_pressure(0.1, 0.0, 10.0), rel=1e-14)
    p7 = dimensional_check(lambda rho, Gamma, r: -rho * Gamma ** 2 / (8 * np.pi ** 2 * r ** 2), "[pressure]",
                           rho=Q_(1000, "kg/m**3"), Gamma=Q_(1, "m**2/s"), r=Q_(0.1, "m"))
    assert p7.to("Pa").magnitude == pytest.approx(ch05.line_vortex_pressure(0.1, 0.0, 1.0), rel=1e-14)
    s = dimensional_check(lambda mu, Gamma, r: -mu * Gamma / (np.pi * r ** 2), "[pressure]", mu=Q_(1e-3, "Pa*s"),
                          Gamma=Q_(1, "m**2/s"), r=Q_(0.1, "m"))
    assert s.to("Pa").magnitude == pytest.approx(ch05.line_vortex_viscous_stress(0.1, 1.0, 1e-3), rel=1e-14)
    tq = dimensional_check(lambda r, sig: 2 * np.pi * r ** 2 * sig, "[force]", r=Q_(0.1, "m"), sig=s)  # N m per m
    assert tq.to("N").magnitude == pytest.approx(-2e-3, rel=1e-12)
    pw = dimensional_check(lambda mu, Gamma, a: mu * Gamma ** 2 / (np.pi * a ** 2), "[power]/[length]",
                           mu=Q_(1e-3, "Pa*s"), Gamma=Q_(1, "m**2/s"), a=Q_(0.1, "m"))
    assert pw.to("W/m").magnitude == pytest.approx(ch05.dissipation_outside_cylinder(0.1, 1.0, 1.0)["power_in"],
                                                   rel=1e-13)


# =====================================================================================================================
# C03 — Kelvin's circulation theorem (5.8)–(5.11); N10–N14, N16, R09; D04, D05
# =====================================================================================================================
def test_loop_helpers_V1_spectral_quadrature_on_closed_loops():  # V1 (the ∮ machinery of (5.9)–(5.11), (5.33))
    pts = ch05.circle_loop_points((0.3, -0.2), 0.7, 64)
    assert pts.shape == (2, 64) and np.allclose(np.hypot(pts[0] - 0.3, pts[1] + 0.2), 0.7)
    assert ch05.loop_length(pts) == pytest.approx(2 * np.pi * 0.7, rel=1e-13)
    assert ch05.loop_vector_area(pts) == pytest.approx([0.0, 0.0, np.pi * 0.49], rel=1e-13, abs=1e-15)
    T = ch05.loop_tangent(pts)  # dx/ds = 2π·0.7·(−sin, cos)
    th = 2 * np.pi * np.arange(64) / 64
    assert np.allclose(T, 2 * np.pi * 0.7 * np.stack([-np.sin(th), np.cos(th)]), atol=1e-11)
    # ∮ u·dx of solid-body rotation b e_z × x = 2b·area; of a gradient field = 0 (exact differential)
    u = lambda x, t=0.0: np.stack([-0.8 * np.asarray(x)[1], 0.8 * np.asarray(x)[0]])  # noqa: E731
    assert ch05.loop_circulation(u, pts) == pytest.approx(2 * 0.8 * np.pi * 0.49, rel=1e-13)
    grad = lambda x, t=0.0: np.stack([np.cos(np.asarray(x)[0]) * np.asarray(x)[1], np.sin(np.asarray(x)[0])])  # noqa
    assert abs(ch05.loop_circulation(grad, pts)) < 1e-13
    # the smoothed square: area side², length 4 side, still spectral
    sq = ch05.square_loop_points((1.0, 2.0), 0.4, 512)
    assert ch05.loop_vector_area(sq)[2] == pytest.approx(0.16, rel=1e-10)
    assert ch05.loop_circulation(u, sq) == pytest.approx(2 * 0.8 * 0.16, rel=1e-10)
    # the length integrand |dx/ds| is only C² at the smoothed corners: algebraic, not spectral, convergence
    ns = [128, 256, 512, 1024]
    le = [abs(ch05.loop_length(ch05.square_loop_points((1.0, 2.0), 0.4, n)) - 1.6) for n in ns]
    assert le[2] < 1e-6 * 1.6 and observed_order([1 / n for n in ns], le) > 2.5
    # a tilted 3-D circle: vector area = πr² n (n its normal), the right-hand sense
    n = np.array([1.0, 2.0, 2.0]) / 3.0
    c3 = ch05.circle_loop_points((0.1, 0.2, 0.3), 0.5, 128, normal=n)
    assert np.allclose(ch05.loop_vector_area(c3), np.pi * 0.25 * n, atol=1e-13)
    assert ch05.loop_line_integral(u(c3[:2]), c3[:2]) == pytest.approx(2 * 0.8 * np.pi * 0.25 * n[2], rel=1e-12)


def test_material_circulation_V4_cellular_loop_stretched_gamma_constant():  # V4 (Kelvin (5.8), N13 E3 "cellular")
    s = ch05.kelvin_scenario("cellular")
    assert s["t_end"] == pytest.approx(6 * 2 * np.pi, rel=1e-14) and s["broken"] is None
    G0 = ch05.kelvin_scenario_circulation("cellular", 0.0)
    assert G0 == pytest.approx(1.155130, rel=1e-6)
    # flux route (independent): ∫∫ 2 sin x sin y dA over the disc r = 0.5 about (π/2, 0.9)
    fl, _ = dblquad(lambda th, r: 2 * np.sin(np.pi / 2 + r * np.cos(th)) * np.sin(0.9 + r * np.sin(th)) * r,
                    0, 0.5, 0, 2 * np.pi, epsabs=0, epsrel=1e-12)
    assert G0 == pytest.approx(fl, rel=1e-11)
    t = np.linspace(0.0, s["t_end"], 7)
    Gt, loops = ch05.material_circulation(s["u"], s["pts0"], t, return_loops=True)
    Ls = [ch05.loop_length(lp) for lp in loops]
    assert Ls[-1] / Ls[0] > 6.5  # the loop is stretched ×6.7 …
    assert np.max(np.abs(Gt - G0)) < 1e-9 * G0  # … and its circulation does not move (contract 1e-12)
    # the loop is material: each point lies on its own pathline (independent integration of one particle)
    k = 137
    pl = K.pathline(s["u"], s["pts0"][:, k], 0.0, t, 1e-12, 1e-14)
    assert np.max(np.abs(pl - loops[:, :, k].T)) < 1e-8
    # the stream function is materially conserved on every loop point (steady flow: ψ constant on pathlines)
    psi = lambda P: np.sin(P[0]) * np.sin(P[1])  # noqa: E731
    assert np.max(np.abs(psi(loops[-1]) - psi(loops[0]))) < 1e-8


@pytest.mark.parametrize("name", ["rankine_straddle", "gaussian"])
def test_material_circulation_V4_loop_straddling_a_core(name):  # V4 ((5.8): steady vortex, loop across the core)
    s = ch05.kelvin_scenario(name)
    t = np.linspace(0.0, 6.0, 4)
    Gt = ch05.material_circulation(s["u"], s["pts0"], t)
    if name == "rankine_straddle":
        # Γ(0) = 2 × lens area (ω = 2 in the core, 0 outside): circle r = 1 ∩ circle r = 0.5 about (0.8, 0)
        d, r1, r2 = 0.8, 1.0, 0.5
        lens = (r1 ** 2 * np.arccos((d * d + r1 * r1 - r2 * r2) / (2 * d * r1))
                + r2 ** 2 * np.arccos((d * d + r2 * r2 - r1 * r1) / (2 * d * r2))
                - 0.5 * np.sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2)))
        assert Gt[0] == pytest.approx(2 * lens, rel=2e-5)  # the kink at the core edge: algebraic convergence
        assert lens == pytest.approx(0.549106, rel=1e-5)
        assert np.ptp(Gt) < 2e-5 * Gt[0]
    else:
        assert np.ptp(Gt) < 1e-9 * abs(Gt[0])


def test_material_circulation_V1_lamb_oseen_decays_as_closed_form():  # V1 ((5.11) at work: N13)
    Gr, dG = ch05.lamb_oseen_circulation(0.005, 10.0, 0.01, 1e-6)
    assert (Gr, dG) == pytest.approx((0.00464739, -3.34538e-4), rel=1e-5)
    s = ch05.kelvin_scenario("lamb_oseen")
    t = np.array([0.0, 10.0, 30.0, 60.0])
    Gt = ch05.material_circulation(s["u"], s["pts0"], t)
    exact = 0.01 * (1 - np.exp(-0.005 ** 2 / (4e-6 * (t + 10.0))))
    assert np.allclose(Gt, exact, rtol=1e-11)  # the circle is material (u_r = 0) and Γ drains outward
    assert np.allclose(ch05.lamb_oseen_circulation(0.005, t, 0.01, 1e-6, t0=10.0)[0], exact, rtol=1e-14)
    # ∂Γ/∂t by a central difference of the closed form = the analytic rate
    h = 1e-3
    fd = (ch05.lamb_oseen_circulation(0.005, 10 + h, 0.01, 1e-6)[0] - ch05.lamb_oseen_circulation(0.005, 10 - h, 0.01,
                                                                                                 1e-6)[0]) / (2 * h)
    assert fd == pytest.approx(dG, rel=1e-7)
    assert ch05.lamb_oseen_circulation(10.0, 1.0, 0.01, 1e-6)[0] == pytest.approx(0.01, rel=1e-14)  # r → ∞: all of Γ
    # the viscous line integral of (5.11), computed from the field, equals ∂Γ/∂t (N13, D05 step 7)
    for (r, tt) in ((0.005, 10.0), (0.003, 5.0), (0.01, 40.0)):
        assert ch05.lamb_oseen_viscous_loop_integral(r, tt, 0.01, 1e-6) == pytest.approx(
            ch05.lamb_oseen_circulation(r, tt, 0.01, 1e-6)[1], rel=1e-6)


def _manufactured_u(x, t=0.0):
    """An unsteady plane field that is not an Euler solution (Γ of material loops changes)."""
    X = np.asarray(x, dtype=float)
    return np.stack([np.sin(X[1]) * (1 + 0.3 * t) + 0.2 * X[0], 0.5 * np.sin(X[0]) - 0.2 * X[1]])


def test_kelvin_rate_V3_total_equals_dGamma_dt_of_material_loop():  # V3 ((5.9) is an identity of kinematics)
    pts0 = ch05.circle_loop_points((0.4, 0.1), 0.6, 256)
    t0 = 0.7
    loop_t0 = ch05.material_loop(_manufactured_u, pts0, [0.0, t0])[-1]
    kr = ch05.kelvin_rate_terms(_manufactured_u, loop_t0, t0, h=1e-5, ht=1e-5)
    assert isinstance(kr, ch05.KelvinRate) and kr.total == pytest.approx(kr.acceleration + kr.contour, abs=1e-15)
    assert abs(kr.contour) < 1e-8 * abs(kr.acceleration)  # ∮u·du = ∮d(½u²) = 0 (N11)
    assert abs(kr.acceleration) > 0.05
    dts, errs = [0.2, 0.1, 0.05, 0.025], []
    for dt in dts:
        G = ch05.material_circulation(_manufactured_u, pts0, [0.0, t0 - dt, t0, t0 + dt])
        errs.append(abs((G[3] - G[1]) / (2 * dt) - kr.total))
    assert abs(observed_order(dts, errs) - 2.0) < ORDER_TOL
    assert errs[-1] < 1e-4 * abs(kr.total)


def test_kelvin_rate_V2_derivation():  # V2 — D04: dΓ/dt of a material loop = ∮(Du/Dt)·dx + ∮d(½|u|²)
    s = sp.Symbol("s", real=True)
    Xf, Yf = sp.Function("X")(s, TT), sp.Function("Y")(s, TT)
    x, y = sp.symbols("x y", real=True)
    uu = x ** 2 * y + TT * y + sp.sin(x)  # a generic (not divergence-free, not Euler) unsteady field
    vv = -x * y ** 2 + x * TT + sp.cos(y)
    on = {x: Xf, y: Yf}
    # steps 1–2: Γ(t) = ∫ u(X(s,t), t)·∂X/∂s ds; step 3: d/dt at fixed labels, with the particles moving at u
    integrand = uu.subs(on) * sp.diff(Xf, s) + vv.subs(on) * sp.diff(Yf, s)
    ddt = sp.diff(integrand, TT)
    material = {sp.Derivative(Xf, TT): uu.subs(on), sp.Derivative(Yf, TT): vv.subs(on)}
    ddt = ddt.subs({sp.Derivative(Xf, s, TT): sp.diff(uu.subs(on), s), sp.Derivative(Yf, s, TT): sp.diff(vv.subs(on), s)})
    ddt = ddt.subs(material)
    # (5.9): acceleration term (Du/Dt along the loop) + contour term u_i D(dx_i)/Dt with D(dx)/Dt = du (steps 5–6)
    Du = (sp.diff(uu, TT) + uu * sp.diff(uu, x) + vv * sp.diff(uu, y)).subs(on)
    Dv = (sp.diff(vv, TT) + uu * sp.diff(vv, x) + vv * sp.diff(vv, y)).subs(on)
    accel = Du * sp.diff(Xf, s) + Dv * sp.diff(Yf, s)
    contour = uu.subs(on) * sp.diff(uu.subs(on), s) + vv.subs(on) * sp.diff(vv.subs(on), s)
    assert sp.simplify(sp.expand(ddt - accel - contour)) == 0
    # step 7: the contour term is an exact s-derivative, ∂/∂s(½|u|²) — its integral over a closed loop vanishes (step 8)
    ke = (uu ** 2 + vv ** 2).subs(on) / 2
    assert sp.simplify(sp.expand(contour - sp.diff(ke, s))) == 0
    # the trap: a multi-valued function (the polar angle) is not single-valued round the origin: ∮dθ = 2π
    th = sp.Symbol("theta")
    assert sp.integrate(sp.diff(th, th), (th, 0, 2 * sp.pi)) == 2 * sp.pi


def test_kelvin_force_terms_V1_barotropic_zero_baroclinic_stokes():  # V1 ((5.10), R09; −∮dp/ρ = ∫(∇ρ×∇p/ρ²)·n dA)
    sq = ch05.square_loop_points((0.0, 0.5), 0.2, 1024)
    lf = ch05.lock_exchange_fields(1000.0, 1025.0, 0.1)
    gp = lambda x, t=0.0: lf["grad_p"](np.asarray(x)[0], np.asarray(x)[1])  # noqa: E731
    kf = ch05.kelvin_force_terms(sq, gp, lf["rho_fn"])
    assert isinstance(kf, ch05.KelvinForces) and kf.total == pytest.approx(kf.pressure + kf.body + kf.viscous)
    # independent Stokes route: ∫∫ (∇ρ×∇p)_z/ρ² over the square (dblquad of the closed-form gradients)
    def bz(y, x):  # (∇ρ × ∇p)_z/ρ² from the closed-form gradients
        gr, gp = lf["grad_rho"](x, y), lf["grad_p"](x, y)
        return float((gr[0] * gp[1] - gr[1] * gp[0]) / lf["rho"](x, y) ** 2)
    stokes, _ = dblquad(bz, -0.1, 0.1, 0.4, 0.6, epsabs=0, epsrel=1e-11)
    assert kf.pressure == pytest.approx(stokes, rel=1e-8) and kf.pressure > 0  # counterclockwise (heavy on the left)
    assert ch05.kelvin_scenario("baroclinic")["rate"] == pytest.approx(stokes, rel=1e-7)
    # wrong variant: the reference density ρ₀ instead of ρ(x) (the "environment" density) removes the source
    assert abs(ch05.kelvin_force_terms(sq, gp, 1012.5).pressure) < 1e-12
    # a barotropic field: ρ = ρ(p) (isentropic power law) with an arbitrary smooth p → ∮dp/ρ = 0
    pf = lambda x, t=0.0: 1e5 * (1 + 0.1 * np.sin(np.asarray(x)[0]) * np.cos(2 * np.asarray(x)[1]))  # noqa: E731
    rf = lambda x, t=0.0: 1.2 * (pf(x) / 1e5) ** (1 / 1.4)  # noqa: E731
    circ = ch05.circle_loop_points((0.2, 0.3), 0.8, 512)
    kb = ch05.kelvin_force_terms(circ, ch05.stencil_gradient_fn(pf, 1e-5), rf)
    scale = np.max(np.abs(ch05.stencil_gradient_fn(pf, 1e-5)(circ))) / 1.2 * 0.8
    assert abs(kb.pressure) < 1e-9 * scale
    # body forces: conservative g = −∇Φ → 0; a non-conservative g = c(−y, x) → 2c·area
    gPhi = lambda x, t=0.0: np.stack([np.cos(np.asarray(x)[0]), 2 * np.asarray(x)[1]])  # noqa: E731
    assert abs(ch05.kelvin_force_terms(circ, None, 1.0, grad_Phi=gPhi).body) < 1e-12
    gnc = lambda x, t=0.0: 0.3 * np.stack([-np.asarray(x)[1], np.asarray(x)[0]])  # noqa: E731
    assert ch05.kelvin_force_terms(circ, None, 1.0, g=gnc).body == pytest.approx(0.6 * np.pi * 0.64, rel=1e-12)


def test_kelvin_V2_derivation():  # V2 — D05: (5.10) → (5.11) → (5.8) and the three sources
    x, y = sp.symbols("x y", real=True)
    p = sp.Function("p")(x, y)
    P = sp.Symbol("P")
    rho_of_p = sp.Function("rho")
    # step 6: barotropic ρ = ρ(p): (1/ρ)∇p is a gradient (of the pressure function), so its curl vanishes
    F = [sp.diff(p, x) / rho_of_p(p), sp.diff(p, y) / rho_of_p(p)]
    assert sp.simplify(sp.diff(F[1], x) - sp.diff(F[0], y)) == 0
    # step 9 (baroclinic): curl(−∇p/ρ)_z = (∇ρ × ∇p)_z/ρ² for independent ρ(x, y), p(x, y) — by Stokes the loop term
    r_ = sp.Function("rho")(x, y)
    Fb = [-sp.diff(p, x) / r_, -sp.diff(p, y) / r_]
    curl = sp.diff(Fb[1], x) - sp.diff(Fb[0], y)
    assert sp.simplify(curl - (sp.diff(r_, x) * sp.diff(p, y) - sp.diff(r_, y) * sp.diff(p, x)) / r_ ** 2) == 0
    # step 5: a conservative body force does no work round a loop: curl(−∇Φ) = 0
    Phi = sp.Function("Phi")(x, y)
    assert sp.simplify(sp.diff(-sp.diff(Phi, y), x) - sp.diff(-sp.diff(Phi, x), y)) == 0
    # step 7 with Lamb–Oseen: ∂Γ/∂t on the circle r equals ∮ν∇²u·dx = 2πr ν(∇²u)_θ  ((5.11))
    r, t, nu, G0 = sp.symbols("r t nu Gamma_0", positive=True)
    ut = G0 / (2 * sp.pi * r) * (1 - sp.exp(-r ** 2 / (4 * nu * t)))
    lap_t = sp.diff(ut, r, 2) + sp.diff(ut, r) / r - ut / r ** 2
    Gam = G0 * (1 - sp.exp(-r ** 2 / (4 * nu * t)))
    assert sp.simplify(sp.diff(Gam, t) - 2 * sp.pi * r * nu * lap_t) == 0
    assert sp.simplify(2 * sp.pi * r * ut - Gam) == 0  # Γ = 2πr u_θ
    # the incompressible, constant-μ caveat: μ∇²u = −μ∇×ω only when ∇·u = 0 (compressible adds ∇(∇·u))
    u3 = [x ** 2 * y, x * y, 0]
    X = (x, y, sp.Symbol("z"))
    lapu = s_lap(sp.Matrix(u3), X)
    curlw = s_curl(s_curl(sp.Matrix(u3), X), X)
    assert sp.simplify(lapu + curlw - s_grad(s_div(sp.Matrix(u3), X), X)) == sp.zeros(3, 1)
    assert s_div(sp.Matrix(u3), X) != 0 and sp.simplify(lapu + curlw) != sp.zeros(3, 1)


def test_kelvin_hypotheses_V1_all_sixteen_combinations():  # V1 (N16 decision table ↔ the surviving terms of (5.10))
    import itertools
    names = ["viscous", "baroclinic", "body", "coriolis"]
    for flags in itertools.product([True, False], repeat=4):
        d = ch05.kelvin_hypotheses(*flags)
        broken = [n for n, ok in zip(names, flags) if not ok]
        assert d["surviving_terms"] == broken and d["holds"] == (not broken)
        txt = ch05.kelvin_hypotheses_text(*flags)
        assert d["verdict"] == txt
        if not broken:
            assert txt.startswith("Kelvin holds") and "(5.8)" in d["text"]
        else:
            assert txt.startswith("Kelvin fails") and all("$" in d["text"] for _ in broken)
            if "viscous" in broken:
                assert "(5.11)" in d["text"]
            if "coriolis" in broken:
                assert "(5.33)" in d["text"]
    assert ch05.kelvin_hypotheses(True, False, True, True)["surviving_terms"] == ["baroclinic"]


def _scenario_scale(s, t):
    P = s["pts0"] if t == 0 else ch05.material_loop(s["u"], s["pts0"], [0.0, t])[-1]
    U = np.asarray(s["u"](P, t))
    return float(np.max(np.linalg.norm(U, axis=0)) ** 2 + 1e-30), P


@pytest.mark.parametrize("name,t", [("cellular", 0.0), ("cellular", 5.0), ("rankine_straddle", 0.0),
                                    ("gaussian", 2.0), ("lamb_oseen", 0.0), ("lamb_oseen", 20.0), ("baroclinic", 0.0),
                                    ("rotating", 0.0), ("rotating", 3.0), ("helmholtz_abc", 1.0),
                                    ("helmholtz", 0.5)])
def test_kelvin_scenario_rate_V1_every_scenario_full_split(name, t):  # V1 (E3: every term of (5.9)/(5.10) per scene)
    s = ch05.kelvin_scenario(name)
    r = ch05.kelvin_scenario_rate(name, t)
    assert set(r) == {"acceleration", "contour", "pressure", "body", "viscous", "coriolis", "total"}
    assert r["total"] == pytest.approx(r["acceleration"] + r["contour"], abs=1e-15)
    assert r["body"] == 0.0
    u2, P = _scenario_scale(s, t)
    tol = 1e-7 * u2  # stencil level, relative to |u|² (the m²/s² scale of the loop terms per unit length)
    assert abs(r["contour"]) < tol
    if name in ("cellular", "rankine_straddle", "gaussian", "helmholtz_abc", "helmholtz"):
        # inviscid, barotropic, inertial: every term 0 and DΓ/Dt = 0
        tol2 = (2e-4 if name == "rankine_straddle" else 1e-7) * u2  # the Rankine kink: stencils straddle the edge
        for k in ("acceleration", "pressure", "viscous", "coriolis", "total"):
            assert abs(r[k]) < tol2, (k, r[k])
    if name == "cellular":  # the p of the scenario is the steady Euler pressure: Du/Dt = −∇p/ρ at the loop points
        a = K.acceleration(s["u"], P, t, 1e-5).a
        gp = s["grad_p"](P, t) / s["rho"]
        assert np.max(np.abs(a + gp)) < 1e-6 * np.max(np.abs(a))
    if name == "lamb_oseen":
        dG = ch05.lamb_oseen_circulation(5e-3, t, 0.01, 1e-6, t0=10.0)[1]
        assert r["viscous"] == pytest.approx(dG, rel=1e-6) and r["acceleration"] == pytest.approx(dG, rel=1e-5)
        assert r["pressure"] == 0.0 and r["coriolis"] == 0.0
    if name == "baroclinic":
        assert r["acceleration"] == r["pressure"] == pytest.approx(s["rate"], rel=1e-14)
        assert r["viscous"] == 0.0 and r["coriolis"] == 0.0 and r["pressure"] > 0
        assert ch05.kelvin_scenario_circulation("baroclinic", 0.5) == pytest.approx(0.5 * s["rate"], rel=1e-14)
    if name == "rotating":  # DΓ/Dt = −∮(2Ω×u)·dx; closed form Γ(t) = 2ΩA₀(1 − e^{−αt})
        dG = 2 * 0.5 * np.pi * 0.2 * np.exp(-0.2 * t)
        assert r["coriolis"] == pytest.approx(dG, rel=1e-9) and r["acceleration"] == pytest.approx(dG, rel=1e-6)
        assert r["pressure"] == 0.0 and r["viscous"] == 0.0


def test_kelvin_scenario_V1_contract_and_errors():  # V1 (C.1 row 1.10–1.11; kelvin_scenario_gamma alias)
    for name in ch05.KELVIN_SCENARIOS:
        s = ch05.kelvin_scenario(name)
        for k in ("u", "pts0", "t_end", "label", "broken", "Omega", "rho", "nu", "grad_p", "visc_force", "p"):
            assert k in s, (name, k)
        assert s["pts0"].shape[0] == (3 if name in ("helmholtz", "helmholtz_abc") else 2)
    assert ch05.kelvin_scenario("lamb_oseen")["broken"] == "viscous"
    assert ch05.kelvin_scenario("baroclinic")["broken"] == "baroclinic"
    assert ch05.kelvin_scenario("rotating")["broken"] == "inertial"
    assert ch05.kelvin_scenario_circulation("rankine_straddle", 0.0) == pytest.approx(1.098212, rel=1e-6)
    assert ch05.kelvin_scenario_circulation("rotating", 5.0) == pytest.approx(1.985865, rel=1e-6)
    assert ch05.kelvin_scenario_gamma("cellular", 0.0) == ch05.kelvin_scenario_circulation("cellular", 0.0)
    assert ch05.kelvin_scenario_circulation("cellular", 3.0) == pytest.approx(1.155130, rel=1e-6)  # material route
    with pytest.raises(ValueError):
        ch05.kelvin_scenario("nope")


# =====================================================================================================================
# C04 — barotropic vs baroclinic element: the pressure torque (Fig. 5.6), N15 lock exchange; D06, D07
# =====================================================================================================================
def test_pressure_torque_V2_derivation():  # V2 — D06: torque about G / I_G, doubled, is ∇ρ×∇p/ρ₀² (× the exact factor)
    r, th, R, rho0, p0 = sp.symbols("r theta R rho_0 p_0", positive=True)
    a, b, c, d = sp.symbols("a b c d", real=True)  # ∇ρ = (a, b), ∇p = (c, d)
    x, y = r * sp.cos(th), r * sp.sin(th)
    rho = rho0 + a * x + b * y
    p = p0 + c * x + d * y
    # steps 1–3: force −∮p n ds = −πR²∇p, its moment about the centre O vanishes (x = R n on a circle)
    nx, ny = sp.cos(th), sp.sin(th)
    prim = p.subs(r, R)
    F = [sp.integrate(-prim * nx * R, (th, 0, 2 * sp.pi)), sp.integrate(-prim * ny * R, (th, 0, 2 * sp.pi))]
    assert sp.simplify(F[0] + sp.pi * R ** 2 * c) == 0 and sp.simplify(F[1] + sp.pi * R ** 2 * d) == 0
    MO = sp.integrate((R * nx * (-prim * ny) - R * ny * (-prim * nx)) * R, (th, 0, 2 * sp.pi))
    assert sp.simplify(MO) == 0
    # steps 4–5: mass and centre of mass (x_G = R²∇ρ/4ρ₀)
    area = lambda f: sp.integrate(sp.integrate(f * r, (r, 0, R)), (th, 0, 2 * sp.pi))  # noqa: E731
    M = area(rho)
    assert sp.simplify(M - rho0 * sp.pi * R ** 2) == 0
    xG, yG = sp.simplify(area(rho * x) / M), sp.simplify(area(rho * y) / M)
    assert sp.simplify(xG - R ** 2 * a / (4 * rho0)) == 0 and sp.simplify(yG - R ** 2 * b / (4 * rho0)) == 0
    # steps 6–7: moment about G = −x_G × F = (πR⁴/4ρ₀)(∇ρ × ∇p)
    MG = sp.simplify(-(xG * F[1] - yG * F[0]))
    assert sp.simplify(MG - sp.pi * R ** 4 * (a * d - b * c) / (4 * rho0)) == 0
    # step 8: I_G exactly (parallel axes); the book-style πρ₀R⁴/2 is its leading term
    IG = sp.simplify(area(rho * ((x - xG) ** 2 + (y - yG) ** 2)))
    assert sp.simplify(IG - sp.pi * rho0 * R ** 4 / 2 * (1 - (a ** 2 + b ** 2) * R ** 2 / (8 * rho0 ** 2))) == 0
    # steps 9–10: Dω/Dt = 2 M_G / I_G = (∇ρ×∇p)/ρ₀² · 1/(1 − |∇ρ|²R²/8ρ₀²)
    spin = sp.simplify(2 * MG / IG)
    exact = (a * d - b * c) / rho0 ** 2 / (1 - (a ** 2 + b ** 2) * R ** 2 / (8 * rho0 ** 2))
    assert sp.simplify(spin - exact) == 0
    assert sp.simplify(sp.limit(spin, R, 0) - (a * d - b * c) / rho0 ** 2) == 0
    # the trap: torque about O (not G) is zero (MO above) — it would predict no spin-up
    # the code equals the exact result (defaults: R = 1 cm, ∇ρ = (10, 0), ∇p = (0, −9810), ρ₀ = 1000)
    t = ch05.pressure_torque_on_element()
    vals = {a: 10.0, b: 0.0, c: 0.0, d: -9810.0, rho0: 1000.0, R: 0.01, p0: 1e5}
    assert t["torque"] == pytest.approx(float(MG.subs(vals)), rel=1e-9)
    assert t["I_G"] == pytest.approx(float(IG.subs(vals)), rel=1e-9)
    assert t["x_G"] == pytest.approx([float(xG.subs(vals)), 0.0], rel=1e-9, abs=1e-18)
    assert t["force"] == pytest.approx([0.0, float(F[1].subs(vals))], rel=1e-9, abs=1e-12)
    assert t["spin_up"] == pytest.approx(float(exact.subs(vals)), rel=1e-9)
    assert t["ratio"] == pytest.approx(1 / (1 - 100 * 1e-4 / (8 * 1e6)), abs=1e-11)  # 1 + 1.25e-9


def test_pressure_torque_V1_contract_barotropic_and_callables():  # V1 (C.1 row 1.12; Fig. 5.6 barotropic element)
    t = ch05.pressure_torque_on_element()
    assert t["x_G"][0] == pytest.approx(2.5e-7, rel=1e-9) and t["torque"] == pytest.approx(-7.7047e-7, rel=1e-4)
    assert t["I_G"] == pytest.approx(1.57079e-5, rel=1e-5) and t["spin_up"] == pytest.approx(-0.0981, rel=1e-7)
    assert t["baroclinic"] == pytest.approx(-0.0981, rel=1e-14) and t["mass"] == pytest.approx(np.pi * 0.1, rel=1e-10)
    assert np.allclose(t["offset"], t["x_G"])  # centre at the origin
    # barotropic: isopycnals parallel to isobars (both horizontal) → no torque, ratio undefined
    bt = ch05.pressure_torque_on_element(grad_rho=(0.0, -3.0), grad_p=(0.0, -9810.0))
    assert abs(bt["torque"]) < 1e-18 and bt["baroclinic"] == 0.0 and np.isnan(bt["ratio"])
    # tilted, parallel gradients anywhere: still no torque
    g = np.array([0.6, -0.8])
    bt2 = ch05.pressure_torque_on_element(grad_rho=5 * g, grad_p=9000 * g, center=(0.3, -0.2))
    assert abs(bt2["spin_up"]) < 1e-9 * (5 * 9000 / 1000.0 ** 2)  # round-off of the rim sums vs |∇ρ||∇p|/ρ₀²
    # callables (f(x, y) form and field form) give the same numbers as the "linear" presets; stencil point value too
    c = np.array([0.3, -0.2])
    pf = lambda x, y: 1e5 + 0.0 * (x - c[0]) - 9810.0 * (y - c[1])  # noqa: E731
    rf = lambda x, t=0.0: 1000.0 + 10.0 * (np.asarray(x)[0] - c[0])  # noqa: E731
    tc = ch05.pressure_torque_on_element(pf, rf, center=c)
    tl = ch05.pressure_torque_on_element(center=c)
    for k in ("torque", "I_G", "spin_up", "mass"):
        assert tc[k] == pytest.approx(tl[k], rel=1e-9)
    assert tc["baroclinic"] == pytest.approx(-0.0981, rel=1e-8)  # through baroclinic_term's stencil
    # a constant density: no torque whatever the pressure field
    assert abs(ch05.pressure_torque_on_element(pf, 1000.0, center=c)["torque"]) < 1e-15


def test_pressure_torque_V3_ratio_converges_order_two_in_R():  # V3 (the torque route → (5.28) as R → 0)
    Rs = [0.4, 0.2, 0.1, 0.05]
    errs = [abs(ch05.pressure_torque_on_element(radius=R, grad_rho=(200.0, 50.0), grad_p=(300.0, -9810.0))["ratio"] - 1)
            for R in Rs]
    assert abs(observed_order(Rs, errs) - 2.0) < ORDER_TOL
    ex = [(200 ** 2 + 50 ** 2) * R ** 2 / 8e6 / (1 - (200 ** 2 + 50 ** 2) * R ** 2 / 8e6) for R in Rs]
    assert np.allclose(errs, ex, rtol=1e-6)


def test_lock_exchange_V2_derivation():  # V2 — D07: t = 0⁺ baroclinic rate at the interface
    x, y, H, g, d = sp.symbols("x y H g delta", positive=True)
    r1, r2 = sp.symbols("rho_1 rho_2", positive=True)
    rb = (r1 + r2) / 2
    rho = rb - (r2 - r1) / 2 * sp.tanh(2 * x / d)  # heavy ρ₂ on the left (x < 0)
    p = rb * g * (H - y)  # hydrostatic with the mean density (fluid still at rest)
    bz = (sp.diff(rho, x) * sp.diff(p, y) - sp.diff(rho, y) * sp.diff(p, x)) / rho ** 2  # (∇ρ × ∇p)_z/ρ²  (5.28)
    at0 = sp.simplify(bz.subs(x, 0))
    assert sp.simplify(at0 - 2 * (r2 - r1) * g / ((r2 + r1) * d)) == 0  # steps 2–6
    assert sp.simplify(sp.diff(rho, x).subs(x, 0) - (r1 - r2) / d) == 0  # step 2: the slope at the centre
    # step 7: positive = counterclockwise for heavy-left; swapping the sides flips it; the circulation per unit height
    assert float(at0.subs({r1: 1000, r2: 1025, g: 9.81, d: 0.1})) > 0
    assert float(at0.subs({r1: 1025, r2: 1000, g: 9.81, d: 0.1})) < 0
    # the source integrated across the interface is independent of δ: ∫ρ′/ρ² dx = [−1/ρ] = 1/ρ₂ − 1/ρ₁, so the column
    # rate is ρ̄gΔρ/(ρ₁ρ₂) per unit height — finite as δ → 0 (the interface becomes a vortex sheet growing at this rate)
    col = sp.diff(p, y) * (1 / r2 - 1 / r1)
    assert sp.simplify(col - rb * g * (r2 - r1) / (r1 * r2)) == 0
    fnum = sp.lambdify(x, bz.subs({r1: 1000, r2: 1025, g: 9.81, d: 0.1, y: 0.5, H: 1}), "numpy")
    num, _ = quad(fnum, -3.0, 3.0, points=[0.0], epsabs=0, epsrel=1e-12, limit=200)
    assert num == pytest.approx(float(col.subs({r1: 1000, r2: 1025, g: 9.81})), rel=1e-9)
    assert ch05.lock_exchange_fields(nx=2001)["circulation_rate"] == pytest.approx(num * 1.0, rel=1e-5)  # H = 1 m
    assert ch05.lock_exchange_initial_vorticity_rate(1000.0, 1025.0, 0.1) == pytest.approx(
        float(at0.subs({r1: 1000, r2: 1025, g: 9.81, d: 0.1})), rel=1e-14)


def test_lock_exchange_V1_field_route_sense_and_contract():  # V1 (N15: 2.4222 s⁻²; the field route; WV ∇p×∇ρ)
    rate = ch05.lock_exchange_initial_vorticity_rate(1000.0, 1025.0, 0.1)
    assert rate == pytest.approx(2.422222, rel=1e-6)
    lf = ch05.lock_exchange_fields(1000.0, 1025.0, 0.1, nx=129)  # odd nx: the grid contains the interface x = 0
    bt = ch05.baroclinic_term(lf["rho_fn"], lf["p_fn"], np.array([0.0, 0.5]), 0.0, 1e-5)
    assert bt[2] == pytest.approx(rate, rel=1e-7) and bt[0] == bt[1] == 0.0
    assert lf["rate"] == rate and np.max(lf["baroclinic_z"]) == pytest.approx(rate, rel=2e-3)
    assert lf["rho"](-1.0, 0.5) == pytest.approx(1025.0, rel=1e-6) and lf["rho"](1.0, 0.5) == pytest.approx(1000.0,
                                                                                                             rel=1e-6)
    # the wrong order ∇p × ∇ρ gives the opposite sense (clockwise)
    gr, gp = lf["grad_rho"](0.0, 0.5), lf["grad_p"](0.0, 0.5)
    assert (gp[0] * gr[1] - gp[1] * gr[0]) / lf["rho"](0.0, 0.5) ** 2 == pytest.approx(-rate, rel=1e-12)
    assert ch05.baroclinic_rate_2d(gr, gp, lf["rho"](0.0, 0.5)) == pytest.approx(rate, rel=1e-12)
    assert ch05.baroclinic_rate_2d([10, 0], [0, -9810], 1000) == pytest.approx(-0.0981, rel=1e-14)
    # the tank-integrated source = ∫ rate × δ-wide profile: finite as δ → 0 (the interface becomes a vortex sheet)
    c1 = ch05.lock_exchange_fields(delta=0.1, nx=801)["circulation_rate"]
    c2 = ch05.lock_exchange_fields(delta=0.05, nx=801)["circulation_rate"]
    assert c1 == pytest.approx(c2, rel=1e-3) and c1 > 0
    # swapping the sides reverses the spin; equal densities give no source
    assert ch05.lock_exchange_initial_vorticity_rate(1025.0, 1000.0, 0.1) == pytest.approx(-rate, rel=1e-14)
    assert ch05.lock_exchange_initial_vorticity_rate(1000.0, 1000.0, 0.1) == 0.0


def test_baroclinic_term_V1_barotropic_zero_V3_order_V2_symbolic():  # V1/V3/V2 ((5.28))
    x, y, z = XYZ
    rho_e = 1000 + 3 * x + sp.sin(y) * z + x * y
    p_e = 1e5 - 9810 * z + 50 * sp.cos(x) * y
    Bs = ch05.baroclinic_term_sym(rho_e, p_e, XYZ)
    assert is_zero(Bs - s_grad(rho_e).cross(s_grad(p_e)) / rho_e ** 2)
    rf, pf = sca_field(rho_e), sca_field(p_e)
    P = RNG.uniform(-1, 1, (3, 6))
    exact = np.array(sp.lambdify(XYZ, list(Bs), "numpy")(*P), dtype=float)
    hs = [0.1, 0.05, 0.025, 0.0125]
    errs = [np.max(np.abs(ch05.baroclinic_term(rf, pf, P, 0.0, h) - exact)) for h in hs]
    assert abs(observed_order(hs, errs) - 2.0) < ORDER_TOL
    # barotropic ρ = ρ(p) with an arbitrary p: zero (V1)
    pb = sca_field(1e5 + 300 * sp.sin(x + 2 * y) * z)
    rb = lambda X, t=0.0: 1.2 * (pb(X) / 1e5) ** (1 / 1.4)  # noqa: E731
    scale = np.max(np.linalg.norm(st.grad(rb, P, 0.0, 1e-4), axis=0) * np.linalg.norm(st.grad(pb, P, 0.0, 1e-4), axis=0)
                   / rb(P) ** 2)
    assert np.max(np.abs(ch05.baroclinic_term(rb, pb, P, 0.0, 1e-4))) < 1e-8 * scale
    # a density varying with x under hydrostatic p(z): (∇ρ × ∇p) along −y × … = e_y ρ_x p_z / ρ² sign check
    bt = ch05.baroclinic_term(lambda X, t=0.0: 1000 + 2 * np.asarray(X)[0], lambda X, t=0.0: -9810 * np.asarray(X)[2],
                              np.array([0.0, 0.0, 0.0]), 0.0, 1e-4)
    assert bt == pytest.approx([0.0, 2 * 9810 / 1e6, 0.0], rel=1e-8, abs=1e-15)  # (2,0,0)×(0,0,−9810) = (0, 2·9810, 0)


def test_baroclinic_element_scenario_V1_closed_form_equals_numeric():  # V1 (E4 readouts)
    for tilt in (0.0, 0.3, 0.5, 1.2, -0.4):
        e = ch05.baroclinic_element_scenario(tilt, 5.0, 0.05, numeric=True)
        assert e["baroclinic"] == pytest.approx(-5.0 * G * np.sin(tilt) / 1000.0, abs=1e-15)
        assert e["spin_up"] == pytest.approx(e["baroclinic"], rel=1e-12, abs=1e-18)
        n = e["numeric"]
        if tilt == 0.0:
            assert e["status"].startswith("barotropic") and e["sense"] == "none" and abs(n["torque"]) < 1e-15
        else:
            assert e["status"].startswith("baroclinic")
            assert n["torque"] == pytest.approx(e["torque"], rel=1e-8)
            assert n["spin_up"] / e["spin_up"] == pytest.approx(1 / (1 - 25 * 0.0025 / 8e6), rel=1e-9)
            assert e["sense"] == ("clockwise" if tilt > 0 else "counterclockwise")
        assert np.allclose(e["offset"], 0.05 ** 2 * e["grad_rho"] / 4000.0)


def test_baroclinic_V2_dimensions():  # V2 (pint: (5.28), the lock rate, the element spin-up)
    b = dimensional_check(lambda gr, gp, rho: gr * gp / rho ** 2, "1/[time]**2", gr=Q_(10, "kg/m**4"),
                          gp=Q_(-9810, "Pa/m"), rho=Q_(1000, "kg/m**3"))
    assert b.to("1/s**2").magnitude == pytest.approx(ch05.baroclinic_rate_2d([10, 0], [0, -9810], 1000), rel=1e-14)
    lr = dimensional_check(lambda dr, g, rs, d: 2 * dr * g / (rs * d), "1/[time]**2", dr=Q_(25, "kg/m**3"),
                           g=Q_(G, "m/s**2"), rs=Q_(2025, "kg/m**3"), d=Q_(0.1, "m"))
    assert lr.to("1/s**2").magnitude == pytest.approx(ch05.lock_exchange_initial_vorticity_rate(1000, 1025, 0.1),
                                                      rel=1e-14)
    tq = dimensional_check(lambda R, gr, gp, r0: np.pi * R ** 4 * gr * gp / (4 * r0), "[force]",  # N m per m
                           R=Q_(0.01, "m"), gr=Q_(10, "kg/m**4"), gp=Q_(-9810, "Pa/m"), r0=Q_(1000, "kg/m**3"))
    assert tq.to("N").magnitude == pytest.approx(ch05.pressure_torque_on_element()["torque"], rel=1e-8)


# =====================================================================================================================
# C05 — Helmholtz's vortex theorems; N17, N51; D08
# =====================================================================================================================
def test_frozen_in_V1_inviscid_vortex_lines_are_material():  # V1/V4 (Helmholtz 1 by the field equation)
    for u, x0 in ((ch05.abc_flow(1.0, 1.0, 1.0), [0.3, 0.2, 0.1]), (ch05.abc_flow(1.0, 0.6, 0.3), [1.0, -0.5, 2.0]),
                  (ch05.stretched_gaussian_vortex_field(1.0, 1.0, 0.5, 0.0), [0.5, 0.2, 0.3])):
        d = ch05.frozen_in_check(u, x0, t_span=(0.0, 2.0))
        assert np.max(d["angle"]) < 1e-8  # δx ∥ ω for all t
        assert np.max(np.abs(d["ratio"] - 1)) < 1e-8  # |ω|/|δx| constant (integrated ω)
        assert np.max(np.abs(d["ratio_field"] - 1)) < 1e-7  # … and with the Eulerian ω(x(t)) of the field itself
        assert np.allclose(d["omega"], d["omega_field"], rtol=1e-7, atol=1e-9)
        assert d["x"].shape == (3, 21) and d["t"][-1] == 2.0
    # the stretched tube: |ω| at the particle grows as e^{αt} × the Gaussian factor — ratio route checked above;
    # the element δx started along ω stretches (|δx| grows) exactly like ω
    d = ch05.frozen_in_check(ch05.stretched_gaussian_vortex_field(1.0, 1.0, 0.5, 0.0), [0.0, 0.0, 0.3],
                             t_span=(0.0, 2.0))
    assert np.linalg.norm(d["delta"][:, -1]) / np.linalg.norm(d["delta"][:, 0]) == pytest.approx(np.exp(1.0), rel=1e-8)


def test_frozen_in_V1_viscous_burgers_breaks_it():  # V1 (wrong-hypothesis check: ν > 0 ⇒ lines are not material)
    ub = ch05.burgers_vortex_field(1e-3, 1.0, 1e-6)
    d = ch05.frozen_in_check(ub, [1e-3, 0.0, 1e-3], t_span=(0.0, 2.0), nu=1e-6)
    assert d["ratio"][-1] < 0.5 and d["ratio_field"][-1] == pytest.approx(d["ratio"][-1], rel=1e-3)
    # ω of the steady Burgers field at the moving particle falls, δx along the axis stretches as e^{αt}
    assert np.linalg.norm(d["delta"][:, -1]) / np.linalg.norm(d["delta"][:, 0]) == pytest.approx(np.exp(2.0), rel=1e-6)


def test_helmholtz_V4_patch_on_a_tube_wall_keeps_zero_flux():  # V4 (N17 / Fig. 5.7; E3's Helmholtz modes)
    s = ch05.kelvin_scenario("helmholtz")
    t = np.linspace(0.0, s["t_end"], 5)
    Gt = ch05.material_circulation(s["u"], s["pts0"], t)
    assert np.max(np.abs(Gt)) < 1e-10  # zero flux through the patch, and it stays zero
    sa = ch05.kelvin_scenario("helmholtz_abc")
    G, loops = ch05.material_circulation(sa["u"], sa["pts0"], np.linspace(0, 2, 5), return_loops=True)
    assert G[0] == pytest.approx(8.016e-10, rel=1e-3) and np.ptp(G) < 1e-6 * G[0]
    for lp in loops:  # the carried patch stays on a tube wall: its normal ⟂ ω (contract ≤ 2e-5)
        A = ch05.loop_vector_area(lp)
        w = sa["u"](lp.mean(axis=1))
        assert abs(A @ w) / (np.linalg.norm(A) * np.linalg.norm(w)) < 2e-5
    # contrast: a patch facing ω (normal ∥ ω) carries flux |ω|·A ≈ 5e-4 ≫ 8e-10
    w0 = sa["u"](np.zeros(3))
    face = ch05.circle_loop_points((0, 0, 0), 0.01, 256, normal=w0)
    assert ch05.loop_circulation(sa["u"], face) == pytest.approx(np.linalg.norm(w0) * np.pi * 1e-4, rel=1e-4)


def test_helmholtz_V2_derivation():  # V2 — D08: vortex surfaces are material; the ABC lines are streamlines
    R, t, al, s0, Gm, nu = sp.symbols("R t alpha sigma_0 Gamma nu", positive=True)
    # inviscid stretched Gaussian vortex: σ² = σ₀²e^{−αt}; flux function χ(R, t) through the circle R
    s2 = s0 ** 2 * sp.exp(-al * t)
    chi = Gm / (2 * sp.pi) * (1 - sp.exp(-R ** 2 / s2))
    uR = -al * R / 2
    assert sp.simplify(sp.diff(chi, t) + uR * sp.diff(chi, R)) == 0  # Dχ/Dt = 0: tube walls χ = const are material
    # with viscosity (σ² = σ₀²e^{−αt} + (4ν/α)(1 − e^{−αt})) the tube walls are not material (Helmholtz fails)
    s2v = s0 ** 2 * sp.exp(-al * t) + 4 * nu / al * (1 - sp.exp(-al * t))
    chiv = Gm / (2 * sp.pi) * (1 - sp.exp(-R ** 2 / s2v))
    assert sp.simplify(sp.diff(chiv, t) + uR * sp.diff(chiv, R)) != 0
    # the ABC flow is Beltrami (ω = u): its vortex lines are its streamlines, which in a steady flow are pathlines
    A, B, C = sp.symbols("A B C")
    x, y, z = XYZ
    u = sp.Matrix([A * sp.sin(z) + C * sp.cos(y), B * sp.sin(x) + A * sp.cos(z), C * sp.sin(y) + B * sp.cos(x)])
    assert is_zero(s_curl(u) - u)
    # and it is a steady Euler solution: (u·∇)u = ∇(½|u|²) − u × ω = ∇(½|u|²) (the Lamb vector vanishes)
    assert is_zero(s_dir(u, u) - s_grad(u.dot(u) / 2))
    # steps 3–5: Stokes on a patch lying in a tube wall — zero flux ⇒ zero circulation (numerically in the V4 test)


# =====================================================================================================================
# C06 — the vorticity equation (5.12)–(5.13); N18–N20, N22, N23, R10, R11; D09 ★★★
# =====================================================================================================================
def test_vorticity_equation_V2_derivation():  # V2 — D09 ★★★ re-run step by step for generic fields (not via fluidpy)
    x, y, z = XYZ
    nu, rho = sp.symbols("nu rho", positive=True)
    A = sp.Matrix([sp.Function(f"A{i}")(x, y, z, TT) for i in range(3)])
    u = s_curl(A)  # a completely generic divergence-free velocity (u = ∇ × A)
    assert sp.expand(s_div(u)) == 0
    w = s_curl(u)
    p, Phi = sp.Function("p")(x, y, z, TT), sp.Function("Phi")(x, y, z)
    # steps 2–3: the curls of −∇p/ρ (ρ constant) and of g = −∇Φ vanish identically
    assert is_zero(s_curl(-s_grad(p) / rho)) and is_zero(s_curl(-s_grad(Phi)))
    # step 4: ∇ × ∂u/∂t = ∂ω/∂t (Schwarz)
    assert is_zero(s_curl(sp.diff(u, TT)) - sp.diff(w, TT))
    # steps 5–6: (u·∇)u = ∇(½u·u) + ω × u (Lamb), so ∇ × [(u·∇)u] = ∇ × (ω × u); the book's ∇(u·u) is also curl-free
    assert is_zero(s_dir(u, u) - s_grad(u.dot(u) / 2) - w.cross(u))
    assert is_zero(s_curl(s_grad(u.dot(u))))
    # step 7: ∇ × (ν∇²u) = ν∇²ω
    assert is_zero(s_curl(nu * s_lap(u)) - nu * s_lap(w))
    # step 9: (B.3.10) for two generic, unrelated vector fields a, b (general form with both divergences)
    a = sp.Matrix([sp.Function(f"a{i}")(x, y, z) for i in range(3)])
    b = sp.Matrix([sp.Function(f"b{i}")(x, y, z) for i in range(3)])
    assert is_zero(s_curl(a.cross(b)) - (s_dir(b, a) - s_dir(a, b) + a * s_div(b) - b * s_div(a)))
    # steps 10–12: curl of the NS residual ≡ the (5.13) residual Dω/Dt − (ω·∇)u − ν∇²ω for this divergence-free u
    ns_res = sp.diff(u, TT) + s_dir(u, u) + s_grad(p) / rho + s_grad(Phi) - nu * s_lap(u)
    r513 = sp.diff(w, TT) + s_dir(u, w) - s_dir(w, u) - nu * s_lap(w)
    assert is_zero(s_curl(ns_res) - r513)
    # without ∇·u = 0 the two differ by ω(∇·u) (the term step 10 drops): a compressible polynomial field
    uc = sp.Matrix([x ** 2 * y, y * z + x, z ** 2 * x * TT])
    wc = s_curl(uc)
    diff_c = s_curl(sp.diff(uc, TT) + s_dir(uc, uc) - nu * s_lap(uc)) - (
        sp.diff(wc, TT) + s_dir(uc, wc) - s_dir(wc, uc) - nu * s_lap(wc))
    assert is_zero(diff_c - wc * s_div(uc)) and not is_zero(wc * s_div(uc))


def test_vorticity_equation_sym_V2_function_outputs():  # V2 (C.1 row 1.15 on generic and on test fields)
    x, y, z = XYZ
    nu = sp.Symbol("nu", positive=True)
    P, Q, R = [sp.Function(n)(x, y, z, TT) for n in "PQR"]
    d = ch05.vorticity_equation_sym([P, Q, R], XYZ, TT, nu, p_expr=sp.Function("p")(x, y, z, TT),
                                    Phi_expr=sp.Function("Phi")(x, y, z))
    assert is_zero(d["curl_pressure"]) and is_zero(d["curl_gravity"]) and is_zero(d["identity_B310"])
    assert is_zero(d["curl_local"] - sp.diff(d["omega"], TT))
    # generic (not divergence-free) field: residual_513 = ω(∇·u), not 0
    assert is_zero(d["residual_513"] - d["omega"] * d["div_u"])
    # the design's check field: u = (yz², xz, −xy) is divergence-free → residual 0; Hill's vortex likewise
    for uf in ([y * z ** 2, x * z, -x * y],
               [x * z / 5, y * z / 5, (1 - 2 * (x ** 2 + y ** 2) - z ** 2) / 5]):  # Hill inside, A = a = 1
        dd = ch05.vorticity_equation_sym(uf, XYZ, TT, nu)
        assert dd["div_u"] == 0 and is_zero(dd["residual_513"]) and is_zero(dd["div_omega"])
    # Hill: ω = A R e_φ = A(−y, x, 0) and (ω·∇)u = (u·∇)ω — steady, inviscid, and ∇²ω = 0
    dh = ch05.vorticity_equation_sym([x * z / 5, y * z / 5, (1 - 2 * (x ** 2 + y ** 2) - z ** 2) / 5], XYZ, TT, nu)
    assert is_zero(dh["omega"] - sp.Matrix([-y, x, 0]))
    assert is_zero(dh["lhs_513"] - dh["rhs_513"])


@pytest.mark.parametrize("name", ["lamb_oseen", "taylor_green", "burgers", "hill", "abc", "stretched_viscous"])
def test_vorticity_terms_V1_exact_solutions_balance(name):  # V1 ((5.13) residual on exact NS/Euler solutions)
    wfn = None
    if name == "lamb_oseen":  # exact ω supplied: ω from u by a nested stencil would sit at its round-off floor here
        u = ch05.planar_field_3d(ch05.lamb_oseen_field(0.01, 1e-6, 10.0))
        P, nu, h, ht = np.array([[4e-3, 1e-3, -2e-3], [1e-3, 6e-3, 3e-3], [0.0, 0.0, 0.5]]), 1e-6, 5e-6, 1e-3

        def wfn(x, t=0.0):
            X = np.asarray(x)
            s2 = 4e-6 * (t + 10.0)
            wz = 0.01 / (np.pi * s2) * np.exp(-(X[0] ** 2 + X[1] ** 2) / s2)
            return np.stack([0 * wz, 0 * wz, wz])
    elif name == "taylor_green":
        def u(x, t=0.0):
            X = np.asarray(x)
            e = np.exp(-2 * 0.01 * t)
            return np.stack([np.cos(X[0]) * np.sin(X[1]) * e, -np.sin(X[0]) * np.cos(X[1]) * e, 0 * X[2]])
        P, nu, h, ht = RNG.uniform(-1, 1, (3, 4)), 0.01, 1e-3, 1e-3
    elif name == "burgers":
        u = ch05.burgers_vortex_field(1e-3, 1.0, 1e-6)
        P, nu, h, ht = np.array([[1e-3, 5e-4, -2e-3], [0.0, 1e-3, 1e-3], [0.0, 2e-3, -1e-3]]), 1e-6, 2e-6, None
    elif name == "hill":
        u = ch05.hill_spherical_vortex_field(1.0, 1.0)
        P, nu, h, ht = np.array([[0.5, 0.2, -0.3], [0.0, 0.3, 0.1], [0.3, -0.4, 0.2]]), 0.0, 1e-3, None
    elif name == "abc":
        u = ch05.abc_flow(1.0, 0.8, 0.5)
        P, nu, h, ht = RNG.uniform(-2, 2, (3, 4)), 0.0, 1e-3, None
    else:  # the stretched Gaussian vortex with viscosity: an exact unsteady NS solution (every term nonzero)
        u = ch05.stretched_gaussian_vortex_field(1.0, 1.0, 0.5, 0.02)
        P, nu, h, ht = np.array([[0.5, 0.2, 1.0], [0.1, 0.8, -0.5], [0.3, 0.2, 0.1]]), 0.02, 1e-3, 1e-4
    b = ch05.vorticity_terms(u, P, 0.7, nu, h=h, ht=ht, omega_fn=wfn)
    assert isinstance(b, ch05.VorticityTerms)
    big = max(np.max(np.abs(getattr(b, k))) for k in ("local", "advective", "stretching_tilting", "diffusion"))
    if name == "lamb_oseen":  # (u·∇)ω = 0 only by cancellation (u ⟂ ∇ω): its stencil truncation scales with |u||ω|/σ
        sig = np.sqrt(4e-6 * 10.7)
        big = max(big, float(np.max(np.linalg.norm(u(P, 0.7), axis=0) * np.abs(wfn(P, 0.7)[2]))) / sig)
    assert np.max(np.abs(b.residual)) < 1e-5 * big
    if name in ("lamb_oseen", "taylor_green"):  # plane: no stretching or tilting; local = diffusion
        assert np.max(np.abs(b.stretching_tilting)) < 1e-9 * big
        assert np.allclose(b.local + b.advective, b.diffusion, rtol=1e-5, atol=1e-5 * big)
    if name == "burgers":  # steady: advective = stretching + diffusion, each nonzero (D18)
        assert np.max(np.abs(b.local)) == 0.0
        assert np.min(np.abs(b.stretching_tilting[2])) > 1.0 and np.min(np.abs(b.diffusion[2])) > 1.0
    if name in ("hill", "abc"):  # inviscid steady: advection = stretching/tilting
        assert np.allclose(b.advective, b.stretching_tilting, atol=1e-5 * big) and np.max(np.abs(b.diffusion)) == 0.0
    if name == "stretched_viscous":
        assert min(np.min(np.abs(b.local[2])), np.min(np.abs(b.stretching_tilting[2])), np.min(np.abs(b.diffusion[2]))) > 1e-3


def test_vorticity_terms_V3_residual_second_order():  # V3 (nested stencils of (5.13): order 2 in h)
    u = ch05.burgers_vortex_field(1e-3, 1.0, 1e-6)
    P = np.array([1e-3, 5e-4, 2e-4])
    hs = [4e-5, 2e-5, 1e-5, 5e-6]
    errs = [float(np.max(np.abs(ch05.vorticity_terms(u, P, 0.0, 1e-6, h=h).residual))) for h in hs]
    assert abs(observed_order(hs, errs) - 2.0) < ORDER_TOL
    # the fully nested route (ω from u by stencils, then ∇ω, ∇²ω) on the unsteady Lamb–Oseen vortex: order 2 as well
    ulo = ch05.planar_field_3d(ch05.lamb_oseen_field(0.01, 1e-6, 10.0))
    hs2 = [8e-4, 4e-4, 2e-4, 1e-4]
    e2 = [float(np.max(np.abs(ch05.vorticity_terms(ulo, np.array([4e-3, 1e-3, 0.0]), 0.0, 1e-6, h=h, ht=1e-3).residual)))
          for h in hs2]
    assert abs(observed_order(hs2, e2) - 2.0) < ORDER_TOL
    # wrong variant: the right Burgers field with the wrong viscosity is not a solution (residual O(the terms))
    wrong = ch05.vorticity_terms(u, P, 0.0, 2e-6, h=5e-6)
    assert np.max(np.abs(wrong.residual)) > 0.3 * np.max(np.abs(wrong.stretching_tilting))


def test_vorticity_budget_preset_V1_contract_and_residuals():  # V1 (C.1 row 1.18: E5/E7 scenes)
    out = {n: ch05.vorticity_budget_preset(n) for n in ch05.VORTICITY_BUDGET_PRESETS}
    for n, d in out.items():
        big = max(abs(v) for k, v in d.items() if k != "residual")
        assert abs(d["residual"]) <= 1e-6 * big, n
        assert all(isinstance(v, float) for v in d.values())
    b = out["burgers"]
    assert b["advective"] == pytest.approx(b["stretching_tilting"] + b["diffusion"], rel=1e-5)
    bb = ch05.burgers_balance(1e-3, 1e-3, 1.0, 1e-6)  # closed forms at R = 1 mm (D18)
    assert b["stretching_tilting"] == pytest.approx(bb["stretching"], rel=1e-6)
    assert b["advective"] == pytest.approx(bb["advective"], rel=1e-6) and b["diffusion"] == pytest.approx(
        bb["diffusion"], rel=1e-5)
    assert out["lamb_oseen"]["local"] == pytest.approx(out["lamb_oseen"]["diffusion"], rel=1e-6)
    assert out["lamb_oseen"]["stretching_tilting"] == 0.0
    rc = out["rotating_column"]
    assert rc["planetary"] == pytest.approx(2 * 0.5 * 0.2, rel=1e-12) and rc["local"] == pytest.approx(0.2, rel=1e-9)
    assert rc["stretching_tilting"] == 0.0 and rc["baroclinic"] == 0.0
    le = out["lock_exchange"]
    assert le["baroclinic"] == pytest.approx(2.4222222, rel=1e-7) and le["local"] == le["baroclinic"]
    hc = ch05.vorticity_budget_preset("hill", component=1)
    assert hc["advective"] == pytest.approx(hc["stretching_tilting"], rel=1e-9) and abs(hc["advective"]) > 0.01
    assert ch05.vorticity_budget_preset("burgers", x=2e-3)["stretching_tilting"] == pytest.approx(
        ch05.burgers_balance(2e-3)["stretching"], rel=1e-6)
    with pytest.raises(ValueError):
        ch05.vorticity_budget_preset("nope")


def test_diffusing_vortex_sheet_V2_symbolic_V1_invariants():  # V2/V1 (N22: 1-D diffusion of a sheet)
    y, t, g, nu = sp.symbols("y t gamma nu", positive=True)
    w = g / (2 * sp.sqrt(sp.pi * nu * t)) * sp.exp(-y ** 2 / (4 * nu * t))
    uu = -g / 2 * sp.erf(y / (2 * sp.sqrt(nu * t)))
    assert sp.simplify(sp.diff(w, t) - nu * sp.diff(w, y, 2)) == 0  # (5.13) → ∂ω/∂t = ν∂²ω/∂y²
    assert sp.simplify(w + sp.diff(uu, y)) == 0  # ω_z = −∂u/∂y (u odd)
    assert sp.simplify(sp.integrate(w, (y, -sp.oo, sp.oo)) - g) == 0  # ∫ω dy = γ at every t
    assert sp.limit(uu, y, sp.oo) == -g / 2  # u(±∞) = ∓γ/2: the jump γ is kept
    f = sp.lambdify((y, t, g, nu), [uu, w], ["numpy", "scipy"])
    Y = np.linspace(-5e-3, 5e-3, 11)
    for tt in (0.1, 1.0, 10.0):
        got = ch05.diffusing_vortex_sheet(Y, tt, 1.3, 1e-6)
        assert np.allclose(got, f(Y, tt, 1.3, 1e-6), rtol=1e-13, atol=1e-15)
    assert ch05.diffusing_vortex_sheet(0.0, 1.0, 1.0, 1e-6)[1] == pytest.approx(282.09479, rel=1e-7)
    assert ch05.diffusing_vortex_sheet(1e-3, 1.0, 1.0, 1e-6)[0] == pytest.approx(-0.26024994, rel=1e-7)


def test_diffusing_vortex_sheet_V3_ftcs_agreement_second_order():  # V3 (``core.diffusion`` FTCS vs the closed form)
    from fluidpy.core.diffusion import ftcs_diffusion_1d
    nu, gam, t0, t1, Lh = 1e-6, 1.0, 1.0, 3.0, 0.02
    errs, hs = [], []
    for n in (101, 201, 401):
        yv = np.linspace(-Lh, Lh, n)
        dy = yv[1] - yv[0]
        dt = 0.2 * dy ** 2 / nu
        steps = int(round((t1 - t0) / dt))
        dt = (t1 - t0) / steps
        w0 = ch05.diffusing_vortex_sheet(yv, t0, gam, nu)[1]
        F = ftcs_diffusion_1d(w0, nu, dy, dt, steps, values=(0.0, 0.0), save_every=steps)
        errs.append(np.max(np.abs(F[-1] - ch05.diffusing_vortex_sheet(yv, t1, gam, nu)[1])))
        hs.append(dy)
    assert abs(observed_order(hs, errs) - 2.0) < ORDER_TOL


def test_hill_vortex_V1_wikipedia_form_and_continuity():  # V1 form cross-check (reference/ch05) + ψ continuity
    A, a = 1.3, 0.7
    U = VX.hill_translation_speed(A, a)
    assert U == pytest.approx(2 * A * a ** 2 / 15, rel=1e-15)
    # our ψ (co-moving frame, vortex moving +z) = Wikipedia's with U_w = −U (their vortex moves −z): inside and outside
    for r, th in ((0.3, 0.4), (0.6, 1.2), (1.0, 2.0), (2.5, 0.7)):
        R, z = r * np.sin(th), r * np.cos(th)
        if r <= a:
            w_psi = -(3 * (-U) / 4) * (1 - r ** 2 / a ** 2) * r ** 2 * np.sin(th) ** 2
        else:
            w_psi = ((-U) / 2) * (1 - a ** 3 / r ** 3) * r ** 2 * np.sin(th) ** 2
        assert VX.hill_stream_function(R, z, A, a) == pytest.approx(w_psi, rel=1e-12, abs=1e-15)
        uR, uz, wphi = ch05.hill_spherical_vortex(R, z, A, a)
        if r <= a:
            assert wphi == pytest.approx((15 * U / (2 * a ** 2)) * r * np.sin(th), rel=1e-12)  # −(15U_w/2a²) r sinθ
    # ψ = 0 on the sphere from both sides
    for th in np.linspace(0.2, 2.9, 7):
        assert abs(VX.hill_stream_function(a * np.sin(th), a * np.cos(th), A, a)) < 1e-15
    # inside, the coded velocity is the stream function's: u_R = −R⁻¹∂ψ/∂z, u_z = R⁻¹∂ψ/∂R (central differences)
    h = 1e-6
    for (R, z) in ((0.2, 0.1), (0.4, -0.3), (0.1, 0.5)):
        psi = lambda RR, zz: VX.hill_stream_function(RR, zz, A, a)  # noqa: E731
        uR, uz, _ = ch05.hill_spherical_vortex(R, z, A, a)
        assert uR == pytest.approx(-(psi(R, z + h) - psi(R, z - h)) / (2 * h) / R, rel=1e-7)
        assert uz == pytest.approx((psi(R + h, z) - psi(R - h, z)) / (2 * h) / R, rel=1e-7)
    assert np.isnan(ch05.hill_spherical_vortex(2.0, 0.0, A, a, outside=False)[0])


def test_hill_vortex_V1_outside_velocity_is_the_stream_functions():  # V1 (N23 "our addition": the exterior flow)
    """Outside the sphere the coded velocity must be derived from the coded (correct) exterior ψ = −(U/2)R²(1 − a³/r³):
    u_R = −R⁻¹∂ψ/∂z = +(3/2)Ua³Rz/r⁵. The sphere is a streamline (normal velocity 0 from both sides) and the tangential
    velocity is continuous there ((Aa²/5) sin θ both sides)."""
    A, a = 1.3, 0.7
    h = 1e-6
    psi = lambda RR, zz: VX.hill_stream_function(RR, zz, A, a)  # noqa: E731
    for (R, z) in ((0.8, 0.9), (1.5, -0.4), (0.3, 1.2)):
        uR, uz, _ = ch05.hill_spherical_vortex(R, z, A, a)
        uR_psi = -(psi(R, z + h) - psi(R, z - h)) / (2 * h) / R
        assert uz == pytest.approx((psi(R + h, z) - psi(R - h, z)) / (2 * h) / R, rel=1e-7)
        assert uR == pytest.approx(uR_psi, rel=1e-6), f"u_R({R}, {z}) = {uR:.6g}, −R⁻¹∂ψ/∂z = {uR_psi:.6g}"
    for th in np.linspace(0.2, 2.9, 7):
        Ro, zo = 1.000001 * a * np.sin(th), 1.000001 * a * np.cos(th)
        uo = np.array(ch05.hill_spherical_vortex(Ro, zo, A, a)[:2])
        U = VX.hill_translation_speed(A, a)
        assert abs(uo @ np.array([np.sin(th), np.cos(th)])) < 1e-5 * U, f"normal velocity at θ = {th:.2f}"
        assert uo @ np.array([np.cos(th), -np.sin(th)]) == pytest.approx(A * a ** 2 / 5 * np.sin(th), rel=1e-4)


def test_burgers_vortex_V1_wikipedia_form_pressure_and_ns():  # V1 form cross-check + NS residual with its pressure
    Gam, al, nu = 1e-3, 1.0, 1e-6
    R = np.linspace(0.0, 8e-3, 17)
    uR, uphi, uz, wz = ch05.burgers_vortex(R, 0.3, Gam, al, nu)
    aw = al / 2  # Wikipedia's α (their v_z = 2α_w z)
    with np.errstate(divide="ignore", invalid="ignore"):
        w_uphi = np.where(R > 0, Gam / (2 * np.pi * R) * (1 - np.exp(-aw * R ** 2 / (2 * nu))), 0.0)
    assert np.allclose(uR, -aw * R) and np.allclose(uz, 2 * aw * 0.3) and np.allclose(uphi, w_uphi, rtol=1e-12)
    assert np.allclose(wz, aw * Gam / (2 * np.pi * nu) * np.exp(-aw * R ** 2 / (2 * nu)), rtol=1e-13)
    assert ch05.burgers_core_radius(1.0, 1e-6) == pytest.approx(0.002, rel=1e-15)
    assert wz[0] == pytest.approx(79.577472, rel=1e-7)
    # ω_z is the curl of the field (stencil) and the full NS equation holds with burgers_pressure (V1, 3-D)
    uf = ch05.burgers_vortex_field(Gam, al, nu)
    P = np.array([[1e-3, -5e-4, 2e-3], [5e-4, 1e-3, -1e-3], [0.1, -0.2, 0.05]])
    assert np.allclose(ch05.vorticity_field(uf, 1e-6)(P)[2], ch05.burgers_vortex(np.hypot(P[0], P[1]), P[2], Gam, al,
                                                                                   nu)[3], rtol=1e-6)
    pf = lambda X, t=0.0: ch05.VX.burgers_pressure(np.hypot(np.asarray(X)[0], np.asarray(X)[1]), np.asarray(X)[2],  # noqa
                                                   Gam, al, nu, 1000.0)
    ns = NS.ns_incompressible_terms(uf, pf, P, 0.0, rho=1000.0, mu=1e-3, g=(0.0, 0.0, 0.0), h=2e-6)
    big = np.max(np.abs(ns.advective))
    assert np.max(np.abs(ns.residual)) < 1e-5 * big
    # the pressure's radial balance dp/dR = ρ(u_φ²/R − α²R/4) and its axial part −ρα²z (D31) by central differences
    Rv, h = np.array([5e-4, 2e-3, 6e-3]), 1e-8
    dpdR = (VX.burgers_pressure(Rv + h, 0.0, Gam, al, nu) - VX.burgers_pressure(Rv - h, 0.0, Gam, al, nu)) / (2 * h)
    up = ch05.burgers_vortex(Rv, 0.0, Gam, al, nu)[1]
    assert np.allclose(dpdR, 1000.0 * (up ** 2 / Rv - al ** 2 * Rv / 4), rtol=1e-5)
    bb = ch05.burgers_balance(np.array([0.0, 1e-3, 3e-3]))
    assert np.allclose(bb["residual"], 0.0, atol=1e-12) and bb["core_radius"] == 0.002
    assert bb["omega_z"][0] == pytest.approx(79.577472, rel=1e-7)


def test_lamb_oseen_field_V1_equals_navier_stokes_exact_solution():  # V1 (C.3 row 3.3 = ch04's exact field)
    u = ch05.lamb_oseen_field(0.02, 1e-5, 3.0)
    P = RNG.uniform(-0.03, 0.03, (2, 12))
    for t in (0.0, 2.0):
        ref, _ = NS.exact_solution("lamb_oseen", P, t, Gamma=0.02, mu=1e-2, rho=1000.0, t0=3.0)
        assert np.allclose(u(P, t), ref, rtol=1e-12, atol=1e-16)
    assert np.allclose(u(np.zeros(2)), 0.0)  # finite (zero) on the axis


# =====================================================================================================================
# C07 — velocity from vorticity: Poisson, Green's function (5.14) (sign corrected), (5.15), Biot–Savart (5.16); D10, D11
# =====================================================================================================================
def test_poisson_green_V2_derivation():  # V2 — D10 ★★★: ∇×ω = −∇²u; G = −1/(4πr); the sign of (5.14)
    x, y, z, eps = sp.symbols("x y z epsilon", positive=True)
    u = sp.Matrix([sp.Function(f"u{i}")(x, y, z) for i in range(3)])
    X = (x, y, z)
    # steps 1–3: ∇×(∇×u) = ∇(∇·u) − ∇²u (generic), so ∇×ω = −∇²u when ∇·u = 0
    assert is_zero(s_curl(s_curl(u, X), X) - (s_grad(s_div(u, X), X) - s_lap(u, X)))
    # step 6: 1/r is harmonic off the source
    r = sp.sqrt(x ** 2 + y ** 2 + z ** 2)
    assert sp.simplify(sum(sp.diff(1 / r, v, 2) for v in X)) == 0
    # step 7: flux of ∇(1/r) through a sphere of any radius ε is −4π (∂(1/r)/∂r = −1/ε², area element ε² sinθ)
    th, ph = sp.symbols("theta phi")
    flux = sp.integrate(sp.integrate(-1 / eps ** 2 * eps ** 2 * sp.sin(th), (th, 0, sp.pi)), (ph, 0, 2 * sp.pi))
    assert flux == -4 * sp.pi
    # step 8: G = −1/(4πr) has unit source strength (flux of ∇G = +1); the trap G = +1/(4πr) has −1
    Gs = -1 / (4 * sp.pi * eps)
    assert sp.simplify(sp.diff(Gs, eps) * 4 * sp.pi * eps ** 2) == 1
    # steps 9–11 constructively: φ = ∫G q d³x′ solves ∇²φ = q for a radial source (shell form of the same integral)
    rr, s_ = sp.symbols("r s", positive=True)
    q = sp.exp(-s_ ** 2)
    phi = -(1 / rr) * sp.integrate(q * s_ ** 2, (s_, 0, rr)) - sp.integrate(q * s_, (s_, rr, sp.oo))
    lap = sp.diff(rr ** 2 * sp.diff(phi, rr), rr) / rr ** 2
    assert sp.simplify(lap - sp.exp(-rr ** 2)) == 0  # +q: the −1/(4π) kernel is the right Green's function
    # hence u = ∫G(−∇′×ω) = +(1/4π)∫(∇′×ω)/|x − x′|: two minus signs; the printed −1/(4π) solves ∇²u = +∇×ω instead
    assert ch05.poisson_green_3d(np.array([1.0, 2.0, 2.0]), np.zeros(3)) == pytest.approx(-1 / (12 * np.pi), rel=1e-15)
    pts = RNG.uniform(-1, 1, (3, 5))
    assert np.allclose(ch05.poisson_green_3d(pts, np.zeros((3, 1))), -1 / (4 * np.pi * np.linalg.norm(pts, axis=0)))


def test_velocity_from_curl_omega_V1_corrected_sign_and_printed_sign():  # V1 (N25: +1/(4π); WV −1/(4π))
    F = ch05.gaussian_tube_fields()
    b = F["bounds"]
    ref = F["u_theta_reference"](0.5)
    assert ref == pytest.approx(0.308838, rel=1e-5)
    nodes, w = ch05.cylinder_quadrature(b["radius"], b["z0"], b["z1"], 40, 64, 40)
    x = np.array([[0.5, 0.0, -0.3 * 0.5], [0.0, 0.5, 0.3 * 0.5], [0.0, 0.0, 0.0]])  # r = 0.5 in three directions
    up = ch05.velocity_from_curl_omega(F["curl_omega"], x, nodes, w)  # default sign +1
    um = ch05.velocity_from_curl_omega(F["curl_omega"], x, nodes, w, sign=-1.0)
    for k in range(2):
        r = np.hypot(x[0, k], x[1, k])
        uth = (-up[0, k] * x[1, k] + up[1, k] * x[0, k]) / r
        assert uth == pytest.approx(ref, rel=5e-6)  # counterclockwise for Γ > 0 along +z
        assert (-um[0, k] * x[1, k] + um[1, k] * x[0, k]) / r == pytest.approx(-ref, rel=5e-6)  # the printed sign
        assert abs(up[2, k]) < 1e-14
    # the (5.14) route equals the Biot–Savart (5.16) route on the same nodes (the surface term of (5.15) is 0 here)
    ub = ch05.biot_savart_volume(F["omega"], x, nodes, w)
    assert np.allclose(up, ub, rtol=1e-9, atol=1e-14)


def test_biot_savart_volume_V3_spectral_convergence_and_long_tube_limit():  # V3/V7 (quadrature; L → ∞ gives (5.2))
    F = ch05.gaussian_tube_fields()
    ref = F["u_theta_reference"](0.5)
    b = F["bounds"]
    errs = []
    for n in ((16, 32, 16), (24, 32, 24), (32, 48, 32), (40, 64, 40)):
        nodes, w = ch05.cylinder_quadrature(b["radius"], b["z0"], b["z1"], *n)
        errs.append(abs(ch05.biot_savart_volume(F["omega"], [0.5, 0, 0], nodes, w)[1] / ref - 1))
    assert all(e2 < e1 / 5 for e1, e2 in zip(errs, errs[1:])) and errs[-1] < 5e-6  # faster than any power
    # the reference (1-D quad, z integrated in closed form) → the infinite line value as L grows (V7)
    uinf = F["u_theta_infinite"](0.5)
    assert uinf == pytest.approx(0.3183099, rel=1e-7)
    for L in (40.0, 400.0):
        assert ch05.gaussian_tube_fields(L=L)["u_theta_reference"](0.5) == pytest.approx(uinf, rel=2.5 * (0.5 / L) ** 2)
    # and the finite tube equals the segment law (D13) for the enclosed circulation, up to the core-size correction
    seg = ch05.segment_speed(0.5, np.arctan2(0.5, 2.0), np.pi - np.arctan2(0.5, 2.0), 1.0)
    assert ref == pytest.approx(seg, rel=1e-3)
    # the induced field is solenoidal (stencil divergence of the quadrature sum, points outside the core)
    nodes, w = ch05.cylinder_quadrature(b["radius"], b["z0"], b["z1"], 24, 32, 24)
    ufield = lambda X, t=0.0: ch05.biot_savart_volume(F["omega"], X, nodes, w)  # noqa: E731
    Pp = np.array([[0.7, 0.2, 0.3], [0.1, -0.8, 0.4], [0.5, 1.5, -2.5]])
    dv = [np.max(np.abs(st.div(ufield, Pp, 0.0, h))) for h in (4e-3, 2e-3, 1e-3)]
    assert abs(observed_order([4e-3, 2e-3, 1e-3], dv) - 2.0) < ORDER_TOL  # only the stencil's truncation: ∇·u = 0


def test_velocity_from_vorticity_fft_V1_taylor_green_abc_V3_box():  # V1/V3 (N24: ∇²u = −∇×ω on a periodic box)
    n = 32
    x = 2 * np.pi * np.arange(n) / n
    X, Y = np.meshgrid(x, x, indexing="xy")
    w = 2 * np.sin(X) * np.sin(Y)  # ω_z of u = (sin x cos y, −cos x sin y)
    uv = ch05.velocity_from_vorticity_fft(w, 2 * np.pi)
    assert np.allclose(uv[0], np.sin(X) * np.cos(Y), atol=1e-13) and np.allclose(uv[1], -np.cos(X) * np.sin(Y),
                                                                                 atol=1e-13)
    Z3, Y3, X3 = np.meshgrid(x, x, x, indexing="ij")
    abc = ch05.abc_flow(1.0, 0.7, 0.3)(np.stack([X3, Y3, Z3]))  # (3, nz, ny, nx), ω = u (Beltrami)
    assert np.allclose(ch05.velocity_from_vorticity_fft(abc, 2 * np.pi), abc, atol=1e-13)
    # a Gaussian vortex in a periodic box: the periodic images' error falls as the box grows
    errs = []
    for Lb in (1.0, 2.0, 4.0):
        m = int(64 * Lb)
        xs = (np.arange(m) - m // 2) * Lb / m
        Xg, Yg = np.meshgrid(xs, xs, indexing="xy")
        wg = 1.0 / (np.pi * 0.05 ** 2) * np.exp(-(Xg ** 2 + Yg ** 2) / 0.05 ** 2)
        wg -= wg.mean()  # the k = 0 mode is dropped by the solver
        uvg = ch05.velocity_from_vorticity_fft(wg, Lb)
        j = m // 2 + int(round(0.2 * m / Lb))  # the grid point nearest x = 0.2, y = 0
        r0 = xs[j]
        errs.append(abs(uvg[1][m // 2, j] - 1.0 / (2 * np.pi * r0) * (1 - np.exp(-(r0 / 0.05) ** 2))))
    uref = 1.0 / (2 * np.pi * 0.2)
    assert errs[0] > errs[1] > errs[2] and errs[2] < 0.01 * uref  # the periodic images fade as the box grows
    with pytest.raises(ValueError):
        ch05.velocity_from_vorticity_fft(np.zeros(5), 1.0)


def test_biot_savart_2d_V1_outside_spectral_V3_dual_grid_interior():  # V1/V3 (plane kernel; point vortices' law)
    Gam, sig = 1.0, 0.1
    out_err, in_err, hs = [], [], []
    for n in (100, 200, 400):
        x = (np.arange(n) + 0.5) * 2.0 / n - 1.0  # cell centres on [−1, 1]
        X, Y = np.meshgrid(x, x, indexing="xy")
        w = Gam / (np.pi * sig ** 2) * np.exp(-(X ** 2 + Y ** 2) / sig ** 2)
        pts = np.array([[0.3, 0.5, 0.05, 0.1], [0.0, 0.0, 0.0, 0.0]])  # cell corners (the dual grid)
        u = ch05.biot_savart_2d(w, X, Y, pts)
        r = np.hypot(*pts)
        ut = Gam / (2 * np.pi * r) * (1 - np.exp(-r ** 2 / sig ** 2))
        uth = (-u[0] * pts[1] + u[1] * pts[0]) / r
        e = np.abs(uth - ut) / ut
        out_err.append(e[:2].max())
        in_err.append(e[2:].max())
        hs.append(2.0 / n)
    assert max(out_err) < 1e-6 and out_err[-1] < 1e-8  # outside the vorticity the rectangle rule is spectral
    assert observed_order(hs, in_err) > 3.5 and in_err[-1] < 1e-6  # inside: points on the dual grid (symmetric)
    # from-scratch double loop (the notebook's transparent version) equals the vectorised kernel
    n = 40
    x = (np.arange(n) + 0.5) * 2.0 / n - 1.0
    X, Y = np.meshgrid(x, x, indexing="xy")
    w = Gam / (np.pi * sig ** 2) * np.exp(-(X ** 2 + Y ** 2) / sig ** 2)
    p = np.array([0.33, -0.21])
    uu = vv = 0.0
    dA = (2.0 / n) ** 2
    for j in range(n):
        for i in range(n):
            rx, ry = p[0] - X[j, i], p[1] - Y[j, i]
            r2 = rx * rx + ry * ry
            uu += -w[j, i] * dA * ry / (2 * np.pi * r2)
            vv += w[j, i] * dA * rx / (2 * np.pi * r2)
    assert np.allclose(ch05.biot_savart_2d(w, X, Y, p), [uu, vv], rtol=1e-12)
    # a node exactly at the field point contributes nothing (self term)
    assert np.all(np.isfinite(ch05.biot_savart_2d(w, X, Y, np.array([x[3], x[5]]))))


def test_curl_theorem_box_V1_polynomial_exact_V3_midpoint_order():  # V1/V3 ((5.15), N27–N28)
    x, y, z = XYZ
    Fp = [x ** 2 * y + z ** 3, y * z ** 2 - x, x * y * z + y ** 3]
    fF = sp.lambdify(XYZ, Fp, "numpy")
    fC = sp.lambdify(XYZ, list(s_curl(sp.Matrix(Fp))), "numpy")
    F_fn = lambda X, Y, Z: np.array([np.broadcast_to(v, np.shape(X)) for v in fF(X, Y, Z)])  # noqa: E731
    C_fn = lambda X, Y, Z: np.array([np.broadcast_to(v, np.shape(X)) for v in fC(X, Y, Z)])  # noqa: E731
    bnd = ((-0.3, 0.7), (0.1, 1.2), (-1.0, 0.5))
    c = ch05.curl_theorem_box(F_fn, bnd, n=6, curl_fn=C_fn)
    assert isinstance(c, ch05.CurlCheck) and np.max(np.abs(c.diff)) < 1e-12 * np.max(np.abs(c.volume))
    exact = [float(sp.integrate(ci, (x, *bnd[0]), (y, *bnd[1]), (z, *bnd[2]))) for ci in s_curl(sp.Matrix(Fp))]
    assert np.allclose(c.volume, exact, rtol=1e-12) and np.allclose(c.surface, exact, rtol=1e-12)
    cfd = ch05.curl_theorem_box(F_fn, bnd, n=6)  # curl from 4th-order differences
    assert np.allclose(cfd.volume, exact, rtol=1e-9)
    # F = b × x on the unit cube: ∫∇×F dV = 2b
    bvec = np.array([0.2, -0.5, 1.0])
    bx = lambda X, Y, Z: np.array([bvec[1] * Z - bvec[2] * Y, bvec[2] * X - bvec[0] * Z, bvec[0] * Y - bvec[1] * X])  # noqa
    assert np.allclose(ch05.curl_theorem_box(bx, (0.0, 1.0), n=4).surface, 2 * bvec, rtol=1e-12)
    # midpoint rule on a smooth non-polynomial field: both sides converge at order 2 to the same vector
    Fs = lambda X, Y, Z: np.array([np.sin(Y) * np.exp(Z), np.cos(X * Z), np.sin(X + Y)])  # noqa: E731
    ref = ch05.curl_theorem_box(Fs, (0.0, 1.0), n=12, rule="gauss").surface
    hs, e = [], []
    for n in (8, 16, 32):
        cc = ch05.curl_theorem_box(Fs, (0.0, 1.0), n=n, rule="midpoint")
        e.append(np.max(np.abs(cc.surface - ref)))
        hs.append(1.0 / n)
    assert abs(observed_order(hs, e) - 2.0) < ORDER_TOL


def test_biot_savart_V2_derivation():  # V2 — D11 ★★★: product rule, ∇′(1/r), the corrected step 5 and the book's slip
    xs = sp.symbols("x y z", real=True)
    xp = sp.symbols("xp yp zp", real=True)
    dvec = sp.Matrix(xs) - sp.Matrix(xp)
    r = sp.sqrt(dvec.dot(dvec))
    W = sp.Matrix([sp.Function(f"w{i}")(*xp) for i in range(3)])  # generic ω(x′)
    phi = 1 / r
    curlp = lambda A: s_curl(A, xp)  # noqa: E731  (∇′×)
    gradp = s_grad(phi, xp)
    assert is_zero(curlp(phi * W) - phi * curlp(W) - gradp.cross(W))  # step 2 (product rule for a curl)
    assert is_zero(gradp - dvec / r ** 3)  # step 4: ∇′(1/|x − x′|) = +(x − x′)/|x − x′|³
    assert is_zero(-gradp.cross(W) - W.cross(dvec) / r ** 3)  # step 5 (correct): −∇′φ × ω = ω × (x − x′)/r³
    book = (dvec / r ** 3).cross(W)  # the book prints +((x − x′)/r³) × ω here
    assert is_zero((-gradp.cross(W) - book) - 2 * W.cross(dvec) / r ** 3)  # … it differs by 2ω×(x − x′)/r³
    # steps 7–9: the k-th component of ∇′×(ω/r) is a divergence of F_i = ε_kij ω_j/r (Gauss per component), so
    # ∫_V′ ∇′×(ω/r) = ∮ n×(ω/r) — checked numerically by curl_theorem_box; steps 10–11: n×ω = 0 on ends ∥ ω, ω = 0 on
    # the side — the tube quadrature of the (5.14) and (5.16) routes agree (test above)
    for k in range(3):
        Fk = sp.Matrix([sum(sp.LeviCivita(k, i, j) * W[j] * phi for j in range(3)) for i in range(3)])
        assert is_zero(s_div(Fk, xp) - curlp(phi * W)[k])


# =====================================================================================================================
# C08 — the filament law (5.17); N29 (segment, infinite line = (5.2)); D12, D13
# =====================================================================================================================
def test_segment_induced_velocity_V1_equals_quadrature_of_5_17():  # V1 (closed form = quad of the filament law)
    for _ in range(6):
        xa, xb = RNG.uniform(-1, 1, 3), RNG.uniform(-1, 1, 3)
        x = RNG.uniform(-1.5, 1.5, 3)
        Gam = RNG.uniform(0.5, 2.0)
        L = np.linalg.norm(xb - xa)
        e = (xb - xa) / L

        def integrand(l, k):
            xp = xa + l * e
            d = x - xp
            return (Gam / (4 * np.pi) * np.cross(e, d) / np.linalg.norm(d) ** 3)[k]
        ref = np.array([quad(integrand, 0, L, args=(k,), epsabs=0, epsrel=1e-12, limit=200)[0] for k in range(3)])
        assert np.allclose(ch05.segment_induced_velocity(x, xa, xb, Gam), ref, rtol=1e-9, atol=1e-13)
    # vectorised over field points and segments; on the segment's line the value is 0 (kernel undefined)
    P = RNG.uniform(-1, 1, (3, 4))
    A, B = RNG.uniform(-1, 1, (3, 3)), RNG.uniform(-1, 1, (3, 3))
    tot = ch05.segment_induced_velocity(P, A, B, 1.0)
    one = sum(ch05.segment_induced_velocity(P, A[:, j], B[:, j], 1.0) for j in range(3))
    assert np.allclose(tot, one, rtol=1e-13)
    assert np.allclose(ch05.segment_induced_velocity([0, 0, 3.0], [0, 0, -1.0], [0, 0, 1.0], 1.0), 0.0)
    assert ch05.segment_induced_velocity([1, 0, 0], [0, 0, -1], [0, 0, 1], 1.0) == pytest.approx([0, 0.1125395, 0],
                                                                                                abs=1e-7)


def test_segment_V1_infinite_line_recovers_5_2_and_wikipedia_form():  # V1 (D13's closing of the loop; form check)
    for d in (0.3, 1.0, 2.5):
        long = ch05.segment_induced_velocity([d, 0, 0], [0, 0, -1e6], [0, 0, 1e6], 1.0)
        assert long[1] == pytest.approx(1 / (2 * np.pi * d), rel=1e-10)  # Γ/2πd = (5.2), along e_z × e_d = e_y
        half = ch05.segment_induced_velocity([d, 0, 0], [0, 0, 0], [0, 0, 1e6], 1.0)
        assert half[1] == pytest.approx(1 / (4 * np.pi * d), rel=1e-10)  # semi-infinite: half
    assert ch05.segment_speed(1.0, 0.0, np.pi, 1.0) == pytest.approx(1 / (2 * np.pi), rel=1e-14)
    assert ch05.segment_speed(1.0, 0.0, np.pi / 2, 1.0) == pytest.approx(0.0795775, rel=1e-6)
    # Wikipedia (aerodynamics): v = Γ/(4πr)(cos A − cos B) with the angles at the ends; same numbers as the vector form
    th_a, th_b = np.arctan2(1.0, 1.0), np.pi - np.arctan2(1.0, 1.0)
    assert ch05.segment_speed(1.0, th_a, th_b, 1.0) == pytest.approx(1 / (4 * np.pi) * (np.cos(th_a) - np.cos(th_b)),
                                                                     rel=1e-15)
    assert ch05.segment_speed(1.0, th_a, th_b, 1.0) == pytest.approx(0.1125395, rel=1e-6)
    # a reversed segment (Γ along −e_z) reverses the velocity; smoothing eps removes the axis singularity
    assert np.allclose(ch05.segment_induced_velocity([1, 0, 0], [0, 0, 1], [0, 0, -1], 1.0), [0, -0.1125395, 0],
                       atol=1e-7)
    assert np.all(np.isfinite(ch05.segment_induced_velocity([1e-9, 0, 0], [0, 0, -1], [0, 0, 1], 1.0, eps=0.01)))


def test_segment_V2_derivation():  # V2 — D13: the segment integral, the substitution, the infinite limit
    d, l, la, lb, th = sp.symbols("d l l_a l_b theta", real=True)
    dp = sp.Symbol("d", positive=True)
    # step 2: e_z × (d e_d − l e_z) = d e_φ (right-hand rule) with e_d = e_x
    assert sp.Matrix([0, 0, 1]).cross(sp.Matrix([dp, 0, -l])) == sp.Matrix([0, dp, 0])
    integ = dp / (dp ** 2 + l ** 2) ** sp.Rational(3, 2)
    F = sp.integrate(integ, l)
    # step 7 in terms of the element angles cos θ = −l/√(d² + l²): ∫_{la}^{lb} = (cos θ_a − cos θ_b)/d
    cos = lambda ll: -ll / sp.sqrt(dp ** 2 + ll ** 2)  # noqa: E731
    assert sp.simplify(F.subs(l, lb) - F.subs(l, la) - (cos(la) - cos(lb)) / dp) == 0
    # steps 5–6: l = −d cot θ ⇒ d·dl/(d² + l²)^{3/2} = sin θ dθ/d
    ths = sp.Symbol("theta", positive=True)
    sub = integ.subs(l, -dp * sp.cot(ths)) * sp.diff(-dp * sp.cot(ths), ths)
    sub = sp.trigsimp(sub).subs(sp.Abs(sp.sin(ths)), sp.sin(ths))  # 0 < θ < π on the segment: sin θ > 0
    assert sp.simplify(sub - sp.sin(ths) / dp) == 0
    # step 8: θ_a → 0, θ_b → π gives 2/d, i.e. Γ/2πd; semi-infinite (θ_b = π/2) gives half
    assert sp.limit(F.subs(l, lb), lb, sp.oo) - sp.limit(F.subs(l, la), la, -sp.oo) == 2 / dp
    assert sp.limit(F.subs(l, lb), lb, sp.oo) - F.subs(l, 0) == 1 / dp


def test_filament_V3_polygon_ring_converges_order_two_and_contracts():  # V3/V1 (C08 figure; contract numbers)
    Ms = [8, 16, 32, 64, 128]
    vals = [ch05.filament_velocity_preset("ring", 0, 0, 0, M=M, component=2) for M in Ms]
    assert vals[0] == pytest.approx(1.054786, rel=1e-6) and vals[1] == pytest.approx(1.013052, rel=1e-6)
    assert vals[3] == pytest.approx(1.000804, rel=1e-6)
    err = [abs(v - 1.0) for v in vals]  # Γ/2R = 1 m/s for R = 0.5 m
    assert abs(observed_order([1 / M for M in Ms], err) - 2.0) < ORDER_TOL
    assert ch05.ring_axis_velocity(0.0, 0.5, 1.0) == 1.0 and ch05.ring_axis_velocity(0.5, 0.5, 1.0) == pytest.approx(
        0.35355339, rel=1e-8)
    # off the ring plane on the axis too
    for z in (0.2, 1.0):
        e = [abs(ch05.filament_velocity_preset("ring", 0, 0, z, M=M, component=2) - ch05.ring_axis_velocity(z, 0.5, 1.0))
             for M in (32, 64, 128)]
        assert abs(observed_order([1 / 32, 1 / 64, 1 / 128], e) - 2.0) < ORDER_TOL
    # square loop of side 1 at its centre: 2√2Γ/πL (four segments seen at 45°–135° from d = 0.5)
    assert ch05.filament_velocity_preset("square", 0, 0, 0, M=64, component=2) == pytest.approx(2 * np.sqrt(2) / np.pi,
                                                                                               rel=1e-12)
    assert ch05.filament_velocity_preset("straight", 1, 0, 0, component=1) == pytest.approx(0.1125395, rel=1e-6)
    assert ch05.filament_velocity_preset("semi_infinite", 1, 0, 0, component=1) == pytest.approx(1 / (4 * np.pi), rel=1e-3)
    uh = ch05.filament_velocity_preset("helix", 0.0, 0.0, 0.0)
    assert np.all(np.isfinite(uh)) and uh[2] > 0  # a helix winding counterclockwise about +z pushes fluid along +z
    for name in ch05.FILAMENT_PRESETS:
        P = ch05.filament_preset(name, 64)
        assert P.shape[0] == 3 and P.shape[1] in (64, 65) and name in ch05.FILAMENT_CLOSED
    with pytest.raises(ValueError):
        ch05.filament_preset("nope")


def test_filament_contributions_V1_sum_equals_filament_velocity():  # V1 (E6 per-segment arrows)
    P = ch05.filament_preset("helix", 40)
    x = np.array([0.5, -0.2, 0.1])
    parts = ch05.filament_contributions(x, P, 1.3, closed=False)
    assert parts.shape == (3, 40)
    assert np.allclose(parts.sum(axis=1), ch05.filament_velocity(x, P, 1.3, closed=False), rtol=1e-13)
    Pr = ch05.filament_preset("ring", 24)
    assert np.allclose(ch05.filament_contributions(x, Pr, 1.0).sum(axis=1), ch05.filament_velocity(x, Pr, 1.0), rtol=1e-13)


def test_filament_V2_derivation():  # V2 — D12: summing (5.17) round a circle gives the ring's axial velocity
    R, z, G, ph = sp.symbols("R z Gamma phi", positive=True)
    xp = sp.Matrix([R * sp.cos(ph), R * sp.sin(ph), 0])  # element position
    e_w = sp.Matrix([-sp.sin(ph), sp.cos(ph), 0])  # e_ω = e_φ (counterclockwise about +z)
    x = sp.Matrix([0, 0, z])
    dvec = x - xp
    du = G / (4 * sp.pi) * e_w.cross(dvec) / (R ** 2 + z ** 2) ** sp.Rational(3, 2) * R  # dl = R dφ
    u = [sp.simplify(sp.integrate(sp.simplify(c), (ph, 0, 2 * sp.pi))) for c in du]
    assert u[0] == 0 and u[1] == 0
    assert sp.simplify(u[2] - G * R ** 2 / (2 * (R ** 2 + z ** 2) ** sp.Rational(3, 2))) == 0
    assert sp.simplify(u[2].subs(z, 0) - G / (2 * R)) == 0  # the ring centre: Γ/2R
    # step 3's "frozen kernel": moving the source point across a core of radius a changes 1/|x−x′|³ by ~3a/|x − x′|
    a_, D = sp.symbols("a D", positive=True)
    change = sp.series((D - a_) ** -3 / D ** -3 - 1, a_, 0, 2).removeO()
    assert sp.simplify(change - 3 * a_ / D) == 0
    f = sp.lambdify((z, R, G), u[2])
    for zv in (0.0, 0.3, 1.0):
        assert ch05.ring_axis_velocity(zv, 0.5, 1.2) == pytest.approx(float(f(zv, 0.5, 1.2)), rel=1e-14)


def test_ring_ring_velocity_V3_elliptic_form_equals_polygon():  # V3/V1 (N42: m = k² parameter; on the axis)
    pts = np.array([[0.2, 0.0, 0.3], [0.8, 0.0, -0.2], [0.45, 0.0, 0.1], [0.5, 0.0, 0.4]]).T
    uR, uz = ch05.ring_ring_velocity(pts[0], pts[2], 0.5, 0.0, 1.0)
    errs = []
    for M in (32, 64, 128, 256):
        uf = ch05.filament_velocity(pts, ch05.filament_preset("ring", M, R=0.5), 1.0)
        errs.append(max(np.max(np.abs(uf[0] - uR)), np.max(np.abs(uf[2] - uz))))
        assert np.max(np.abs(uf[1])) < 1e-12  # no swirl
    assert abs(observed_order([1 / 32, 1 / 64, 1 / 128, 1 / 256], errs) - 2.0) < ORDER_TOL
    z = np.array([-0.7, 0.0, 0.3, 2.0])
    uR0, uz0 = ch05.ring_ring_velocity(0.0, z, 0.5, 0.1, 1.3)
    assert np.allclose(uR0, 0.0) and np.allclose(uz0, ch05.ring_axis_velocity(z - 0.1, 0.5, 1.3), rtol=1e-12)
    # symmetry: u_R is odd and u_z even about the ring's plane
    a1, b1 = ch05.ring_ring_velocity(0.7, 0.35, 0.5, 0.1, 1.0)
    a2, b2 = ch05.ring_ring_velocity(0.7, -0.15, 0.5, 0.1, 1.0)
    assert a1 == pytest.approx(-a2, rel=1e-12) and b1 == pytest.approx(b2, rel=1e-12)
    # wrong variant (the modulus k passed where scipy expects m = k²) misses the polygon result by O(1)
    R_, z_, R0 = 0.8, -0.2, 0.5
    Dp2 = (R_ + R0) ** 2 + z_ ** 2
    k = np.sqrt(4 * R_ * R0 / Dp2)
    Dm2 = (R_ - R0) ** 2 + z_ ** 2
    uz_wrong = 1.0 / (2 * np.pi * np.sqrt(Dp2)) * (ellipk(k) + (R0 ** 2 - R_ ** 2 - z_ ** 2) / Dm2 * ellipe(k))
    uz_right = 1.0 / (2 * np.pi * np.sqrt(Dp2)) * (ellipk(k ** 2) + (R0 ** 2 - R_ ** 2 - z_ ** 2) / Dm2 * ellipe(k ** 2))
    assert uz_right == pytest.approx(ch05.ring_ring_velocity(R_, z_, R0, 0.0, 1.0)[1], rel=1e-14)
    assert abs(uz_wrong - uz_right) > 0.05 * abs(uz_right)


@needs_ref
def test_ring_self_velocity_V5_kelvin_thin_ring():  # V5 (Kelvin's formula, Wikipedia "Vortex ring")
    f = sp.sympify(ref_json()["kelvin_ring_speed"]["formula"])
    Gm, R, a = sp.symbols("Gamma R a")
    for (Rv, av, Gv) in ((1.0, 0.1, 1.0), (0.3, 0.01, 2.5), (2.0, 0.05, -0.7)):
        ref = float(f.subs({Gm: Gv, R: Rv, a: av}))
        assert ch05.ring_self_velocity(Rv, av, Gv) == pytest.approx(ref, rel=1e-14)
    assert ch05.ring_self_velocity(1.0, 0.1, 1.0) == pytest.approx(0.328816, rel=1e-6)
    assert ch05.ring_self_velocity(1.0, 0.1, 1.0, core="hollow") == pytest.approx(
        1 / (4 * np.pi) * (np.log(80) - 0.5), rel=1e-14)
    assert ch05.RING_CORES["uniform"] == 0.25


# =====================================================================================================================
# C09 — the vorticity equation in a rotating frame with baroclinic generation (5.18)–(5.30); R12–R17, N30–N37; D14, D15
# =====================================================================================================================
def _generic_divfree():
    x, y, z = XYZ
    A = sp.Matrix([sp.Function(f"A{i}")(x, y, z, TT) for i in range(3)])
    return s_curl(A)


def test_rotating_lamb_form_V2_derivation():  # V2 — D14: (5.20) → (5.21)–(5.24) → (5.25), generic fields
    x, y, z = XYZ
    nu = sp.Symbol("nu", positive=True)
    Om = sp.Matrix(sp.symbols("Omega0:3", real=True))
    ug = sp.Matrix([sp.Function(f"g{i}")(x, y, z, TT) for i in range(3)])  # generic (not divergence-free)
    wg = s_curl(ug)
    eps = sp.LeviCivita
    grad_u = lambda u, i, j: sp.diff(u[i], XYZ[j])  # noqa: E731  u_{i,j}
    # (5.21): u_j u_{i,j} = −(u × ω)_i + ½(u_j²)_{,i}  (any field)
    lhs = sp.Matrix([sum(ug[j] * grad_u(ug, i, j) for j in range(3)) for i in range(3)])
    assert is_zero(lhs - (-ug.cross(wg) + s_grad(ug.dot(ug) / 2)))
    # (5.22): ε_ijk ω_k = u_{j,i} − u_{i,j}  (any field)
    for i in range(3):
        for j in range(3):
            assert sp.expand(sum(eps(i, j, k) * wg[k] for k in range(3)) - (grad_u(ug, j, i) - grad_u(ug, i, j))) == 0
    # (5.23): ν u_{i,jj} = −ν ε_ijk ω_{k,j} + ν(u_{j,j})_{,i}; the last term vanishes only for ∇·u = 0
    visc = sp.Matrix([-nu * sum(eps(i, j, k) * sp.diff(wg[k], XYZ[j]) for j in range(3) for k in range(3))
                      for i in range(3)])
    assert is_zero(nu * s_lap(ug) - visc - nu * s_grad(s_div(ug)))
    # (5.24): 2ε_ijkΩ_j u_k = −2ε_ijkΩ_k u_j (relabelled dummies)
    for i in range(3):
        a = sum(2 * eps(i, j, k) * Om[j] * ug[k] for j in range(3) for k in range(3))
        b = sum(-2 * eps(i, j, k) * Om[k] * ug[j] for j in range(3) for k in range(3))
        assert sp.expand(a - b) == 0
    # (5.25) residual ≡ (5.20) residual for a divergence-free generic u with generic p, ρ, Φ (g = −∇Φ)
    u = _generic_divfree()
    w = s_curl(u)
    p, rho, Phi = (sp.Function(n)(x, y, z, TT) for n in ("p", "rho", "Phi"))
    R20 = sp.diff(u, TT) + s_dir(u, u) + 2 * Om.cross(u) + s_grad(p) / rho + s_grad(Phi) - nu * s_lap(u)
    R25 = (sp.diff(u, TT) + s_grad(u.dot(u) / 2 + Phi) - u.cross(w + 2 * Om) + s_grad(p) / rho
           + nu * s_curl(w))
    assert is_zero(R20 - R25)


def test_rotating_vorticity_equation_V2_derivation():  # V2 — D15 ★★★ with sp.LeviCivita, as the design's check cell
    X = XYZ
    nu = sp.Symbol("nu", positive=True)
    Om = sp.symbols("Omega0:3", real=True)
    eps = sp.LeviCivita
    d = lambda f, j: sp.diff(f, X[j])  # noqa: E731  (the comma)
    u = [sp.Function(f"u{i}")(*X, TT) for i in range(3)]  # generic velocity (nothing dropped yet)
    w = [sum(eps(n, q, i) * d(u[i], q) for q in range(3) for i in range(3)) for n in range(3)]  # ω_n = ε_nqi u_{i,q}
    R3 = range(3)
    # step 1's Lamb term (third term of (5.26)) vs steps 5–8 with nothing dropped
    lamb = [-sum(eps(n, q, i) * eps(i, j, k) * d(u[j] * (w[k] + 2 * Om[k]), q)
                 for q in R3 for i in R3 for j in R3 for k in R3) for n in R3]
    book = [sum(-d(u[n], j) * (w[j] + 2 * Om[j]) + u[j] * d(w[n], j) + d(u[j], j) * (w[n] + 2 * Om[n])
                - u[n] * d(w[j], j) for j in R3) for n in R3]
    assert all(sp.expand(a - b) == 0 for a, b in zip(lamb, book))
    # step 3: the gradient term vanishes (ε antisymmetric × symmetric second derivative)
    Phi = sp.Function("Phi")(*X)
    assert all(sp.expand(sum(eps(n, q, i) * d(d(sum(ui ** 2 for ui in u) / 2 + Phi, i), q) for q in R3 for i in R3))
               == 0 for n in R3)
    # steps 10–11: the pressure term is the baroclinic vector (5.28)
    p, rho = sp.Function("p")(*X), sp.Function("rho")(*X)
    pres = [-sum(eps(n, q, i) * d(d(p, i) / rho, q) for q in R3 for i in R3) for n in R3]
    baro = s_grad(rho, X).cross(s_grad(p, X)) / rho ** 2
    assert all(sp.simplify(pres[n] - baro[n]) == 0 for n in R3)
    # steps 12–13: the viscous term is νω_{n,jj} + ν(stuff × ∇·ω) — with ω = ∇ × u the extra term vanishes
    visc = [-nu * sum(eps(n, q, i) * eps(i, j, k) * d(d(w[k], j), q) for q in R3 for i in R3 for j in R3 for k in R3)
            for n in R3]
    assert all(sp.expand(visc[n] - nu * sum(d(d(w[n], j), j) for j in R3)) == 0 for n in R3)
    # steps 7, 9: the two dropped terms vanish for ∇·ω ≡ 0 (always) and ∇·u = 0 (the book's silent drop in (5.27))
    assert all(sp.expand(sum(d(w[j], j) for j in R3)) == 0 for _ in [0])  # (5.18) for any u
    # a divergence-free polynomial u: the full curl of (5.25) equals the (5.30) form
    x, y, z = X
    up = sp.Matrix([y * z ** 2 + TT * x ** 2 * y, x * z - TT * x * y ** 2, -x * y + TT * x ** 3])
    assert sp.expand(s_div(up)) == 0
    wp = s_curl(up)
    Omm = sp.Matrix(Om)
    rp = 1000 + x * y + z * TT
    pp = 1e5 - 9810 * z + x ** 2 * y * TT
    Phip = 9.81 * z - (Om[2] ** 2) * (x ** 2 + y ** 2) / 2
    R25 = (sp.diff(up, TT) + s_grad(up.dot(up) / 2 + Phip) - up.cross(wp + 2 * Omm) + s_grad(pp) / rp
           + nu * s_curl(wp))
    R30 = (sp.diff(wp, TT) + s_dir(up, wp) - s_dir(wp + 2 * Omm, up) - s_grad(rp).cross(s_grad(pp)) / rp ** 2
           - nu * s_lap(wp))
    assert is_zero(s_curl(R25) - R30)
    # the fluidpy symbolic budget: its residual is this same curl, and (5.30) reduces to (5.13) when Ω = 0, ρ = const
    bs = ch05.vorticity_budget_sym(list(up), X, TT, nu, tuple(Om), rp, pp)
    assert is_zero(bs["residual_530"] - s_curl(R25)) and bs["reduces_to_513"] is True
    assert is_zero(bs["planetary"] - s_dir(2 * Omm, up)) and is_zero(bs["baroclinic"] - s_grad(rp).cross(s_grad(pp))
                                                                     / rp ** 2)
    t513 = ch05.vorticity_terms_sym(list(up), X, TT, nu)
    assert is_zero(t513["residual"] - (sp.diff(wp, TT) + s_dir(up, wp) - s_dir(wp, up) - nu * s_lap(wp)))


def test_rotating_ns_V1_corotating_rest_and_inertial_oscillation():  # V1 (R14: (5.20) residual)
    Om, g, rho = 0.5, 9.81, 1000.0
    P = RNG.uniform(-1, 1, (3, 8))
    # solid-body rotation seen from the co-rotating frame: u = 0, p = ρ(Ω²R²/2 − gz), g_eff = −g e_z + Ω²(x, y, 0)
    pf = lambda X, t=0.0: rho * (Om ** 2 * (np.asarray(X)[0] ** 2 + np.asarray(X)[1] ** 2) / 2 - g * np.asarray(X)[2])  # noqa
    gf = lambda X, t=0.0: np.stack([Om ** 2 * np.asarray(X)[0], Om ** 2 * np.asarray(X)[1], -g + 0 * np.asarray(X)[2]])  # noqa
    zero = lambda X, t=0.0: 0.0 * np.asarray(X)  # noqa: E731
    assert np.max(np.abs(ch05.rotating_ns_residual(zero, pf, rho, P, 0.0, 1e-6, (0, 0, Om), gf))) < 1e-8
    # wrong variant: the plain gravity (no centrifugal part) leaves Ω²R
    bad = ch05.rotating_ns_residual(zero, pf, rho, P, 0.0, 1e-6, (0, 0, Om), (0.0, 0.0, -g))
    assert np.allclose(bad[:2], Om ** 2 * P[:2], rtol=1e-6)
    # inertial oscillation u = U(cos ft, −sin ft, 0), f = 2Ω, uniform p: ∂u/∂t + 2Ω × u = 0
    U = 0.3
    uo = lambda X, t=0.0: np.stack([U * np.cos(2 * Om * t) + 0 * np.asarray(X)[0],  # noqa: E731
                                    -U * np.sin(2 * Om * t) + 0 * np.asarray(X)[0], 0 * np.asarray(X)[0]])
    pc = lambda X, t=0.0: -rho * g * np.asarray(X)[2]  # noqa: E731
    for t in (0.0, 1.3, 4.0):
        assert np.max(np.abs(ch05.rotating_ns_residual(uo, pc, rho, P, t, 0.0, (0, 0, Om), (0, 0, -g), ht=1e-5))) < 1e-8
    # wrong variant: the Coriolis term without its factor 2 (as a frame with Ω/2) does not balance
    assert np.max(np.abs(ch05.rotating_ns_residual(uo, pc, rho, P, 1.3, 0.0, (0, 0, Om / 2), (0, 0, -g), ht=1e-5))) > 0.1


def test_rotating_lamb_form_V1_equals_5_20_residual_on_random_fields():  # V1 (identity (5.21)–(5.24), numerically)
    x, y, z = XYZ
    psi = sp.Matrix([sp.sin(y) * z * (1 + TT), x * sp.cos(z) + TT * y, sp.exp(x / 3) * y])
    ue = s_curl(psi)  # divergence-free
    pe = 1e5 + 30 * sp.sin(x * y) - 9810 * z + TT * x
    re = 1000 + 2 * x - y * z
    Phe = 9.81 * z - 0.02 * (x ** 2 + y ** 2)
    u, p, rho, Phi = vec_field(ue), sca_field(pe), sca_field(re), sca_field(Phe)
    gf = vec_field(-s_grad(Phe))
    P = RNG.uniform(-0.8, 0.8, (3, 6))
    Om = (0.1, -0.2, 0.3)
    r20 = ch05.rotating_ns_residual(u, p, rho, P, 0.4, 1e-2, Om, gf, h=1e-4, ht=1e-5)
    lf = ch05.rotating_lamb_form_terms(u, p, rho, Phi, P, 0.4, 1e-2, Om, h=1e-4, ht=1e-5)
    assert set(lf) == {"local", "grad_B", "lamb_abs", "pressure", "viscous_curl", "residual"}
    assert np.allclose(lf["residual"], r20, atol=1e-5 * np.max(np.abs(lf["pressure"])))
    assert np.allclose(lf["residual"], lf["local"] + lf["grad_B"] + lf["lamb_abs"] - lf["pressure"] - lf["viscous_curl"])
    # the exact (5.20) residual from sympy agrees with both (the fields are not a solution: residual ≠ 0)
    Omm = sp.Matrix(Om)
    R20 = sp.diff(ue, TT) + s_dir(ue, ue) + 2 * Omm.cross(ue) + s_grad(pe) / re + s_grad(Phe) - 1e-2 * s_lap(ue)
    ex = vec_field(R20)(P, 0.4)
    assert np.allclose(r20, ex, rtol=1e-6, atol=1e-6 * np.max(np.abs(ex)))


def test_vorticity_budget_V1_rotating_column_and_absolute_vorticity():  # V1 ((5.30) numerically; R18; N39)
    u = VD._rotating_column_field(0.5, 0.2)  # a column at rest in the frame at t = 0, stretched at α
    P = RNG.uniform(-1, 1, (3, 5))
    for t in (0.0, 2.0):
        b = ch05.vorticity_budget(u, P, t, 0.0, (0, 0, 0.5))
        assert isinstance(b, ch05.VorticityBudget)
        # ζ(t) = 2Ω(e^{αt} − 1): Dζ/Dt = (ζ + 2Ω)α — stretching of the absolute vorticity
        zeta = 2 * 0.5 * np.expm1(0.2 * t)
        assert np.allclose(b.local[2], (zeta + 1.0) * 0.2, rtol=1e-7)
        assert np.allclose(b.planetary[2], 2 * 0.5 * 0.2, rtol=1e-9) and np.allclose(b.stretching_tilting[2],
                                                                                   zeta * 0.2, rtol=1e-7, atol=1e-12)
        assert np.max(np.abs(b.residual)) < 1e-7
    # wrong variant: the relative vorticity alone (no planetary term) leaves 2Ωα unexplained
    bb = ch05.vorticity_budget(u, P, 0.0, 0.0, (0, 0, 0.0))
    assert np.allclose(bb.residual[2], 0.2, rtol=1e-7)
    # planetary *tilting* in the budget (N39): u = (s z, 0, 0) under Ω e_z gives Dω_x/Dt ∋ 2Ω ∂u/∂z = 2Ωs, not 2Ω ∂w/∂x
    s_, Om_ = 0.3, 0.7
    ush = lambda X, t=0.0: np.stack([s_ * np.asarray(X)[2], 0 * np.asarray(X)[0], 0 * np.asarray(X)[0]])  # noqa: E731
    bsh = ch05.vorticity_budget(ush, P, 0.0, 0.0, (0, 0, Om_))
    assert np.allclose(bsh.planetary, np.array([2 * Om_ * s_, 0, 0])[:, None], atol=1e-10)
    # and for a generic 3-D field the budget's planetary term is 2GΩ of the stencil gradient (`planetary_vorticity_terms`)
    ug = ch05.abc_flow(1.0, 0.6, 0.3)
    Omv3 = np.array([0.2, -0.4, 0.9])
    bg = ch05.vorticity_budget(ug, P, 0.0, 0.0, Omv3, h=1e-4)
    for k in range(P.shape[1]):
        Gk = K.velocity_gradient_at(ug, P[:, k], 0.0, 1e-4)
        assert np.allclose(bg.planetary[:, k], ch05.planetary_vorticity_terms(Gk, Omv3), rtol=1e-9, atol=1e-12)
    # planetary stretching and tilting 2(Ω·∇)u = 2GΩ = 2Ω ∂u/∂z for Ω = Ω e_z (N39)
    G_ = RNG.normal(size=(3, 3))
    assert np.allclose(ch05.planetary_vorticity_terms(G_, [0, 0, 0.7]), 1.4 * G_[:, 2])
    assert ch05.planetary_vorticity_terms(np.diag([0, 0, 1e-5]), [0, 0, 7.292115e-5])[2] == pytest.approx(1.458423e-9,
                                                                                                         rel=1e-6)
    # absolute vorticity ω + 2Ω round-trips ch03's relative vorticity ω′ = ω − 2Ω
    w = RNG.normal(size=3)
    Omv = np.array([0.1, 0.0, 0.5])
    assert np.allclose(K.vorticity_in_rotating_frame(ch05.absolute_vorticity(w, Omv), Omv), w, atol=1e-15)
    assert ch05.absolute_vorticity(0.0, 0.5) == 1.0  # fluid at rest in the frame: ω_a = 2Ω
    # the lock exchange's baroclinic term by vorticity_budget with ρ, p callables
    lf = ch05.lock_exchange_fields()
    bl = ch05.vorticity_budget(lambda X, t=0.0: 0.0 * np.asarray(X), np.array([0.0, 0.5, 0.0]), 0.0, 0.0, (0, 0, 0),
                               lambda X, t=0.0: lf["rho_fn"](X), lambda X, t=0.0: lf["p_fn"](X), h=1e-5)
    assert bl.baroclinic[2] == pytest.approx(2.4222222, rel=1e-6)


def test_vorticity_divergence_V1_solenoidal_to_round_off():  # V1 ((5.18): ∇·ω = 0 for any flow)
    x, y, z = XYZ
    ue = sp.Matrix([sp.sin(x * y) * z, sp.exp(z / 2) * x ** 2, sp.cos(x + y * z)])  # compressible, arbitrary
    u = vec_field(ue)
    P = RNG.uniform(-1, 1, (3, 10))
    wscale = np.max(np.abs(ch05.vorticity_field(u, 1e-3)(P)))
    assert np.max(np.abs(ch05.vorticity_divergence(u, P, 0.0, 1e-3))) < 1e-8 * wscale / 1e-3
    assert isinstance(ch05.vorticity_divergence(u, P[:, 0], 0.0, 1e-3), float)
    assert is_zero(s_div(s_curl(ue)))  # the identity itself
    # contrast: the "broken tube" (not a curl) has ∇·ω ≠ 0
    assert abs(float(st.div(ch05.broken_tube_field(), np.array([0.01, 0.0, 0.2]), 0.0, 1e-4))) > 1.0


# =====================================================================================================================
# C10 — stretching and tilting (5.31)–(5.32); N21, N38, N53; D16, D17, D18
# =====================================================================================================================
def test_stretching_tilting_split_V1_decomposition_V7_invariance():  # V1/V7 ((5.31): stretching ∥ ω, tilting ⟂ ω)
    for _ in range(10):
        w, G_ = RNG.normal(size=3), RNG.normal(size=(3, 3))
        d = ch05.stretching_tilting_split(w, G_)
        assert np.allclose(d["stretching"] + d["tilting"], G_ @ w, atol=1e-13) and np.allclose(d["total"], G_ @ w)
        assert np.linalg.norm(np.cross(d["stretching"], w)) < 1e-12 * np.linalg.norm(w) ** 2 * np.linalg.norm(G_)
        assert abs(d["tilting"] @ w) < 1e-12 * np.linalg.norm(w) ** 2 * np.linalg.norm(G_)
        es = w / np.linalg.norm(w)
        assert d["rate"] == pytest.approx(es @ G_ @ es, rel=1e-13) and np.allclose(d["e_s"], es)
        # rotate ω and G together: magnitudes unchanged, vectors rotated (V7)
        Q, _ = np.linalg.qr(RNG.normal(size=(3, 3)))
        dr = ch05.stretching_tilting_split(Q @ w, Q @ G_ @ Q.T)
        assert dr["rate"] == pytest.approx(d["rate"], rel=1e-12) and np.allclose(dr["tilting"], Q @ d["tilting"],
                                                                                atol=1e-12)
    # 2-D: ω ∥ e_z with a planar G → neither stretching nor tilting (D17 step 5)
    Gp = np.zeros((3, 3))
    Gp[:2, :2] = RNG.normal(size=(2, 2))
    d2 = ch05.stretching_tilting_split([0, 0, 3.0], Gp)
    assert np.allclose(d2["stretching"], 0) and np.allclose(d2["tilting"], 0) and d2["tilt_rate"] == 0.0
    d0 = ch05.stretching_tilting_split([0, 0, 0], Gp)
    assert d0["rate"] == 0.0 and np.allclose(d0["total"], 0)


def test_uniform_strain_V1_presets_V2_expm_solves_the_equation():  # V1/V2 ((5.32) exact in linear flows; D17 numbers)
    assert ch05.uniform_strain_vorticity([0, 0, 1], "axial_stretch", 2.0)[2] == pytest.approx(np.exp(2.0), rel=1e-13)
    assert ch05.uniform_strain_vorticity([1, 0, 0], "shear_tilt", 2.0) == pytest.approx([1.0, 0.0, 2.0], abs=1e-13)
    assert ch05.uniform_strain_vorticity([0, 0, 1], "planar", 2.0) == pytest.approx([0, 0, 1.0], abs=1e-15)
    assert ch05.uniform_strain_vorticity([0, 0, 1], "axial_compression", 1.0)[2] == pytest.approx(np.exp(-1.0))
    ts = np.array([0.0, 0.5, 1.0])
    assert ch05.uniform_strain_vorticity([0, 0, 2.0], "burgers", ts, rate=0.3).shape == (3, 3)
    # V2: for every preset, d/dt(e^{Gt}ω₀) = G e^{Gt}ω₀ and the flow u = Gx is a steady Euler flow (G² symmetric)
    t = sp.Symbol("t")
    for name in ch05.STRAIN_PRESETS:
        Gm = sp.Matrix(ch05.strain_preset(name, 1.0)).applyfunc(sp.nsimplify)
        E = (Gm * t).exp()
        w0 = sp.Matrix(sp.symbols("w0:3"))
        assert is_zero(sp.diff(E * w0, t) - Gm * E * w0)
        assert is_zero(Gm * Gm - (Gm * Gm).T)  # (u·∇)u = G²x is a gradient ⇒ a pressure exists
    with pytest.raises(ValueError):
        ch05.strain_preset("nope")
    # independent: the frozen-in integration of a particle in u = Gx carries ω exactly as e^{Gt}ω₀
    Gs = ch05.strain_preset("shear_tilt", 0.8) + ch05.strain_preset("axial_stretch", 0.3)
    u = lambda X, t=0.0: np.einsum("ij,j...->i...", Gs, np.asarray(X, float))  # noqa: E731
    w0 = np.array([0.3, 0.2, 1.0])
    d = ch05.frozen_in_check(u, [0.1, 0.1, 0.1], delta0=1e-3 * w0, t_span=(0, 2.0))
    # ω of the linear flow itself is 0-ish (weak carried vorticity): compare the element δx (same equation as ω)
    assert np.allclose(d["delta"][:, -1] / 1e-3, ch05.uniform_strain_vorticity(w0, Gs, 2.0), rtol=1e-8)


def test_stretching_tilting_scenario_V1_statuses_and_numbers():  # V1 (E5 state)
    s = ch05.stretching_tilting_scenario("axial_stretch", 2.0, 1.0, 0.0)
    assert s["growth"] == pytest.approx(np.exp(2.0), rel=1e-12) and s["status"].startswith("stretching")
    assert s["angle"] == pytest.approx(0.0, abs=1e-12) and s["stretch_rate"] == pytest.approx(1.0, rel=1e-12)
    c = ch05.stretching_tilting_scenario("axial_compression", 1.0)
    assert c["status"].startswith("compression") and c["growth"] == pytest.approx(np.exp(-1.0))
    tl0 = ch05.stretching_tilting_scenario("shear_tilt", 0.0, 1.0, np.pi / 2)  # ω ∥ e_x, w = x: Gω ⟂ ω
    assert tl0["status"].startswith("tilting") and abs(tl0["stretch_rate"]) < 1e-15 and tl0["tilt_rate"] == 1.0
    tl = ch05.stretching_tilting_scenario("shear_tilt", 2.0, 1.0, np.pi / 2)
    assert tl["omega"] == pytest.approx([1.0, 0.0, 2.0], abs=1e-12)
    assert tl["stretch_rate"] == pytest.approx(0.4, rel=1e-12) and tl["tilt_rate"] == pytest.approx(0.2, rel=1e-12)
    pl = ch05.stretching_tilting_scenario("planar", 1.0)
    assert pl["status"].startswith("2-D") and pl["growth"] == pytest.approx(1.0)
    st_ = ch05.stretched_tube(2.0, 1.0, 10.0, 1e-4)
    assert (st_["omega"], st_["A"], st_["Gamma"]) == pytest.approx((20.0, 5e-5, 1e-3), rel=1e-14)
    L = np.array([1.0, 1.5, 3.0])
    sL = ch05.stretched_tube(L, 1.0, 10.0, 1e-4)
    assert np.allclose(sL["omega"] * sL["A"], 1e-3) and np.allclose(sL["A"] * L, 1e-4)  # ωA and AL conserved
    assert np.allclose(sL["radius_ratio"], np.sqrt(1 / L))


def test_natural_coordinates_V2_derivation():  # V2 — D16 (and D17's projection): (ω·∇)u = ω ∂u/∂s; the helix frame
    Gs = sp.Matrix(3, 3, sp.symbols("g0:9", real=True))
    th, ph, wm = sp.symbols("theta phi omega", positive=True)
    es = sp.Matrix([sp.sin(th) * sp.cos(ph), sp.sin(th) * sp.sin(ph), sp.cos(th)])  # e_s = ω/|ω|
    w = wm * es
    # steps 1–4: (ω·∇)u = Gω = |ω| G e_s = ω ∂u/∂s (∂u/∂s = G e_s, the derivative along the line)
    assert is_zero(Gs * w - wm * (Gs * es))
    # steps 5–6: the along-line component is |ω|(e_s·G e_s): the stretching rate of a material element along e_s
    en = sp.Matrix([sp.cos(th) * sp.cos(ph), sp.cos(th) * sp.sin(ph), -sp.sin(th)])
    em = es.cross(en)
    assert sp.simplify(es.dot(en)) == 0 and sp.simplify(em.norm() - 1) == 0
    comps = [sp.simplify(v.dot(Gs * w) - wm * v.dot(Gs * es)) for v in (es, en, em)]
    assert comps == [0, 0, 0]
    # the helix of Fig. 5.9: arc-length tangent, curvature a/(a² + c²), torsion c/(a² + c²), book e_n = −N (outward)
    a, c, s = sp.symbols("a c s", positive=True)
    k = sp.sqrt(a ** 2 + c ** 2)
    Xh = sp.Matrix([a * sp.cos(s / k), a * sp.sin(s / k), c * s / k])
    T = sp.diff(Xh, s)
    assert sp.simplify(T.norm() - 1) == 0
    Tp = sp.diff(T, s)
    kappa = sp.simplify(Tp.norm())
    assert sp.simplify(kappa - a / k ** 2) == 0
    N = sp.simplify(Tp / kappa)  # Frenet normal: toward the axis
    B = T.cross(N)
    tau = sp.simplify(-sp.diff(B, s).dot(N))
    assert sp.simplify(tau - c / k ** 2) == 0
    fX, fT, fN = (sp.lambdify((s, a, c), M, "numpy") for M in (Xh, T, N))
    for sv in (0.0, 1.3, 4.0):
        fr = ch05.helix_frame(sv, 1.0, 0.3)
        assert np.allclose(ch05.helix(sv, 1.0, 0.3), np.ravel(fX(sv, 1.0, 0.3)), atol=1e-14)
        assert np.allclose(fr["e_s"], np.ravel(fT(sv, 1.0, 0.3)), atol=1e-14)
        assert np.allclose(fr["e_n"], -np.ravel(fN(sv, 1.0, 0.3)), atol=1e-14)  # the book's convention
        assert np.allclose(fr["e_m"], np.cross(fr["e_s"], fr["e_n"]), atol=1e-15)
    f0 = ch05.helix_frame(0.0)
    assert (f0["curvature"], f0["torsion"]) == pytest.approx((1 / 1.09, 0.3 / 1.09), rel=1e-14)


def test_frenet_frame_V1_numeric_frame_on_the_helix():  # V1 (C.5 row 5.14 numerical frame vs closed form)
    s = np.linspace(0.0, 6.0, 601)
    X = ch05.helix(s, 1.0, 0.3)
    fr = ch05.frenet_frame(X, s)
    assert np.allclose(fr["kappa"][5:-5], 1.0 / 1.09, rtol=1e-4)
    i = 300
    hf = ch05.helix_frame(s[i], 1.0, 0.3)
    assert np.allclose(fr["e_s"][:, i], hf["e_s"], atol=1e-5) and np.allclose(fr["e_n"][:, i], hf["e_n"], atol=1e-4)
    assert np.allclose(ch05.helical_vortex_line(1.0, 0.3, s), X)


def test_burgers_balance_V2_derivation():  # V2 — D18: ω ∝ L (ideal tube); the steady Burgers balance
    L, V, G0 = sp.symbols("L V Gamma", positive=True)
    A = V / L  # step 2: AL = V
    om = G0 / A  # step 1: ωA = Γ
    assert sp.simplify(om - G0 * L / V) == 0 and sp.simplify(om.subs(L, 2 * L) / om) == 2  # step 3: ω ∝ L
    R, al, nu = sp.symbols("R alpha nu", positive=True)
    wz = al * G0 / (4 * sp.pi * nu) * sp.exp(-al * R ** 2 / (4 * nu))
    lhs = -al * R / 2 * sp.diff(wz, R)  # steps 5–6: u_R dω_z/dR
    rhs = al * wz + nu / R * sp.diff(R * sp.diff(wz, R), R)  # stretching + cylindrical diffusion
    assert sp.simplify(lhs - rhs) == 0
    # step 7: the balance is a total derivative; step 8: the bracket vanishes (regular at R = 0)
    f = sp.Function("f")(R)
    tot = sp.diff(nu * R * sp.diff(f, R) + al / 2 * R ** 2 * f, R) / R
    assert sp.simplify(tot - (nu / R * sp.diff(R * sp.diff(f, R), R) + al * f + al * R / 2 * sp.diff(f, R))) == 0
    assert sp.simplify(nu * R * sp.diff(wz, R) + al / 2 * R ** 2 * wz) == 0
    # step 9: the amplitude from ∫ω_z 2πR dR = Γ, and the core radius √(4ν/α)
    assert sp.simplify(sp.integrate(wz * 2 * sp.pi * R, (R, 0, sp.oo)) - G0) == 0
    assert sp.simplify(wz.subs(R, sp.sqrt(4 * nu / al)) / wz.subs(R, 0) - sp.exp(-1)) == 0
    # the curl of the book's velocity field is this ω_z (Exercise 5.12 ↔ our (c))
    up = G0 / (2 * sp.pi * R) * (1 - sp.exp(-al * R ** 2 / (4 * nu)))
    assert sp.simplify(sp.diff(R * up, R) / R - wz) == 0
    bb = ch05.burgers_balance(np.array([5e-4, 2e-3]), 1e-3, 1.0, 1e-6)
    pars = {G0: 1e-3, al: 1.0, nu: 1e-6}
    fl = sp.lambdify(R, [e.subs(pars) for e in (lhs, al * wz, nu / R * sp.diff(R * sp.diff(wz, R), R))])
    for i, Rv in enumerate((5e-4, 2e-3)):
        ad, stt, di = (float(v) for v in fl(Rv))
        assert (bb["advective"][i], bb["stretching"][i], bb["diffusion"][i]) == pytest.approx((ad, stt, di), rel=1e-12)


# =====================================================================================================================
# C11 — Kelvin in a rotating frame (5.33): the absolute circulation; the fluid column; R18, N39, N54; D19, D20
# =====================================================================================================================
def test_absolute_circulation_V4_material_loop_in_rotating_frame():  # V4/V1 ((5.33) with a numerically advected loop)
    s = ch05.kelvin_scenario("rotating")
    t = np.linspace(0.0, 10.0, 6)
    Gt, loops = ch05.material_circulation(s["u"], s["pts0"], t, return_loops=True)
    Ga = np.array([ch05.absolute_circulation(s["u"], lp, s["Omega"], tt)[1] for lp, tt in zip(loops, t)])
    assert np.allclose(Gt, 2 * 0.5 * np.pi * (-np.expm1(-0.2 * t)), rtol=1e-9, atol=1e-12)  # Γ changes …
    assert np.max(np.abs(Ga - np.pi)) < 1e-9  # … Γ_a = 2Ω·A_vec + Γ stays 2ΩA₀ = π
    assert ch05.kelvin_scenario_circulation("rotating", 10.0) == pytest.approx(2.716424, rel=1e-6)
    # fluid at rest in the inertial frame, seen from the rotating frame: u′ = −Ω × x ⇒ Γ = −2ΩA, Γ_a = 0
    Om = 0.3
    u = lambda X, t=0.0: np.stack([Om * np.asarray(X)[1], -Om * np.asarray(X)[0]])  # noqa: E731
    pts = ch05.circle_loop_points((0.4, -0.2), 0.7, 256)
    G, Gav = ch05.absolute_circulation(u, pts, Om)  # scalar Ω = Ω e_z
    assert G == pytest.approx(-2 * Om * np.pi * 0.49, rel=1e-12) and abs(Gav) < 1e-12
    # a tilted loop: Γ_a − Γ = 2Ω·A_vec (not 2ΩA)
    n = np.array([0.0, 0.6, 0.8])
    c3 = ch05.circle_loop_points((0.0, 0.0, 0.0), 0.5, 256, normal=n)
    u3 = lambda X, t=0.0: 0.0 * np.asarray(X)  # noqa: E731
    G3, Ga3 = ch05.absolute_circulation(u3, c3, (0.0, 0.0, 0.5))
    assert G3 == 0.0 and Ga3 == pytest.approx(2 * 0.5 * 0.8 * np.pi * 0.25, rel=1e-12)


def test_absolute_circulation_V2_derivation():  # V2 — D19: DA_vec/Dt = ∮u × dx, so D(Γ + 2Ω·A_vec)/Dt = 0
    s = sp.Symbol("s", real=True)
    Xs = sp.Matrix([sp.Function(f"X{i}")(s, TT) for i in range(3)])  # a generic material loop x = X(s, t)
    U = sp.Matrix([sp.Function(f"U{i}")(s, TT) for i in range(3)])  # u at the loop points, u = ∂X/∂t
    Om = sp.Matrix(sp.symbols("Omega0:3", real=True))
    # step 4: triple product (Ω × u)·dx = Ω·(u × dx)
    dX = sp.diff(Xs, s)
    assert sp.expand(Om.cross(U).dot(dX) - Om.dot(U.cross(dX))) == 0
    # steps 6–7: d/dt(½ x × x_s) = u × x_s + ½ ∂_s(x × u), with ∂X/∂t = U (material) — the last term integrates to 0
    dA = sp.diff(Xs.cross(dX) / 2, TT).subs({sp.Derivative(Xs[i], TT): U[i] for i in range(3)})
    dA = dA.subs({sp.Derivative(Xs[i], s, TT): sp.diff(U[i], s) for i in range(3)})
    assert is_zero(dA - (U.cross(dX) + sp.diff(Xs.cross(U), s) / 2))
    # step 8: the Coriolis loop term −2∮(Ω×u)·dx = −2Ω·∮u × dx = −2Ω·DA_vec/Dt  ⇒ D(Γ + 2Ω·A_vec)/Dt = 0
    # steps 9–10 (Stokes for the uniform field Ω × x): ∮(Ω × x)·dx = 2Ω·A_vec on any smooth closed loop
    th = sp.Symbol("theta")
    loop = sp.Matrix([sp.cos(th) + sp.Rational(1, 5) * sp.cos(2 * th), sp.sin(th) - sp.Rational(1, 3) * sp.sin(3 * th),
                      sp.Rational(1, 2) * sp.sin(th) * sp.cos(th)])
    lhs = sp.integrate(sp.expand(Om.cross(loop).dot(sp.diff(loop, th))), (th, 0, 2 * sp.pi))
    Av = [sp.integrate(sp.expand(c), (th, 0, 2 * sp.pi)) / 2 for c in loop.cross(sp.diff(loop, th))]
    assert sp.simplify(lhs - 2 * Om.dot(sp.Matrix(Av))) == 0
    # and ch05.loop_vector_area reproduces A_vec of that loop numerically
    f = sp.lambdify(th, list(loop), "numpy")
    thv = 2 * np.pi * np.arange(256) / 256
    pts = np.array([np.broadcast_to(np.asarray(v, float), thv.shape) for v in f(thv)])
    assert np.allclose(ch05.loop_vector_area(pts), [float(a) for a in Av], atol=1e-13)


def test_fluid_column_V2_derivation():  # V2 — D20: (ζ + f)/h = const and its small-change limit 2Ω ∂w/∂z
    zeta, f, A, h, C1, C2 = sp.symbols("zeta f A h C1 C2", real=True)
    h0, z0 = sp.symbols("h0 zeta0", real=True)
    sol = sp.solve([sp.Eq((zeta + f) * A, C1), sp.Eq(A * h, C2)], [zeta, A], dict=True)[0]
    q = sp.simplify((sol[zeta] + f) / h)
    assert sp.simplify(q - C1 / C2) == 0  # steps 1–4: h-independent
    zeta_h = (z0 + f) * h / h0 - f  # step 5
    assert sp.simplify((zeta_h + f) / h - (z0 + f) / h0) == 0
    # step 6: Dζ/Dt = (ζ + f)(1/h)Dh/Dt; for ζ ≪ f and (1/h)Dh/Dt = ∂w/∂z this is 2Ω∂w/∂z (f = 2Ω)
    t = sp.Symbol("t")
    hf = sp.Function("h")(t)
    zt = (z0 + f) * hf / h0 - f
    assert sp.simplify(sp.diff(zt, t) - (zt + f) / hf * sp.diff(hf, t)) == 0
    for (hh, hh0, z00, ff) in ((1100.0, 1000.0, 0.0, 1e-4), (800.0, 1000.0, 2e-5, 1.2e-4)):
        assert ch05.column_relative_vorticity(hh, hh0, z00, ff) == pytest.approx(float(zeta_h.subs(
            {h: hh, h0: hh0, z0: z00, f: ff})), rel=1e-12)
    assert ch05.column_relative_vorticity(1100.0, 1000.0, 0.0, 1e-4) == pytest.approx(1.0e-5, rel=1e-12)


def test_fluid_column_V1_ridge_trough_and_ring_of_air():  # V1/V4 (N54; C11 numbers)
    x = np.linspace(-3e5, 3e5, 61)
    for depth, sgn in (("ridge", -1), ("trough", +1)):
        c = ch05.column_over_slope(x, depth=depth, lat_deg=45.0)
        assert np.ptp(c["ratio"]) < 1e-12 * np.max(np.abs(c["ratio"]))  # (ζ + f)/h conserved
        assert np.sign(c["zeta"][30]) == sgn  # squashed → anticyclonic; stretched → cyclonic
        assert c["f"] == pytest.approx(1.03126e-4, rel=1e-5)
    sl = ch05.column_over_slope(x, depth="slope")
    assert np.all(np.diff(sl["zeta"]) > 0)
    cb = ch05.column_over_slope(x, depth=lambda q: 1000 + 0 * q)
    assert np.allclose(cb["zeta"], 0.0, atol=1e-18)
    with pytest.raises(ValueError):
        ch05.column_over_slope(x, depth="nope")
    f30, f45, f60 = (float(ch05.coriolis_parameter(np.deg2rad(d))) for d in (30, 45, 60))
    assert (f30, f45, f60) == pytest.approx((7.29212e-5, 1.03126e-4, 1.26303e-4), rel=1e-5)
    A0 = np.pi * 5e5 ** 2
    G1 = ch05.relative_circulation_after_move(0.0, A0, 30.0, A0, 60.0)
    assert G1 == pytest.approx(-4.19261e7, rel=1e-5) and G1 / A0 == pytest.approx(-5.33820e-5, rel=1e-5)
    assert G1 + f60 * A0 == pytest.approx(0.0 + f30 * A0, rel=1e-12)  # Γ_a = Γ + fA conserved
    # the ring also shrinking: a smaller area at the same latitude spins it up cyclonically
    assert ch05.relative_circulation_after_move(0.0, A0, 45.0, 0.5 * A0, 45.0) == pytest.approx(f45 * 0.5 * A0,
                                                                                               rel=1e-12)


# =====================================================================================================================
# C12 — point vortices: pairs and the centre of vorticity (§5.7); N40, N41, N55, N56; D21
# =====================================================================================================================
def test_point_vortex_velocity_V1_single_vortex_and_from_scratch():  # V1 ((5.2) superposed; self term excluded)
    xv = np.array([[0.3], [-0.2]])
    P = RNG.uniform(-1, 1, (2, 8))
    u = ch05.point_vortex_velocity(P, xv, [2.0])
    d = P - xv
    r2 = np.sum(d * d, axis=0)
    assert np.allclose(u, 2.0 / (2 * np.pi) * np.stack([-d[1], d[0]]) / r2, rtol=1e-13)  # counterclockwise
    assert np.allclose(ch05.point_vortex_velocity(xv[:, 0], xv, [2.0]), 0.0)  # a vortex does not move itself
    assert np.all(np.isfinite(ch05.point_vortex_velocity(xv[:, 0] + 1e-9, xv, [2.0], eps=1e-3)))
    # from-scratch double loop for dx_k/dt (the notebook's transparent version) = point_vortex_rhs
    V = RNG.uniform(-1, 1, (2, 5))
    Gm = RNG.uniform(-1, 1, 5)
    rhs = np.zeros((2, 5))
    for k in range(5):
        for j in range(5):
            if j != k:
                rx, ry = V[0, k] - V[0, j], V[1, k] - V[1, j]
                r2 = rx * rx + ry * ry
                rhs[0, k] += -Gm[j] / (2 * np.pi) * ry / r2
                rhs[1, k] += Gm[j] / (2 * np.pi) * rx / r2
    assert np.allclose(ch05.point_vortex_rhs(V, Gm), rhs, rtol=1e-13)


def test_point_vortex_evolve_V4_invariants_and_V1_pair_orbits():  # V4/V1 (Figs. 5.11–5.12; D21)
    pr = ch05.point_vortex_preset("three_vortices")
    t = np.linspace(0.0, pr["t_end"], 201)
    X = ch05.point_vortex_evolve(pr["xv"], pr["Gamma"], t)
    inv = ch05.point_vortex_invariants(X, pr["Gamma"])
    assert np.max(np.abs(inv["P"] - inv["P"][0])) < 1e-10 and np.ptp(inv["I"]) < 1e-10 and np.ptp(inv["H"]) < 1e-10
    # equal pair: period 2π/((Γ₁ + Γ₂)/2πh²) = 19.7392 s; after one period both are back
    vp = ch05.vortex_pair(1.0, 1.0, 1.0)
    assert (vp["V1"], vp["centre_from_1"], vp["rotation_rate"], vp["period"]) == pytest.approx(
        (0.159155, 0.5, 0.318310, 19.7392), rel=1e-5)
    eq = ch05.point_vortex_preset("equal_pair")
    Xe = ch05.point_vortex_evolve(eq["xv"], eq["Gamma"], np.linspace(0, 10 * vp["period"], 11))
    assert np.max(np.abs(Xe[-1] - eq["xv"])) < 1e-8  # ten orbits: invariants and the rate hold
    Xq = ch05.point_vortex_evolve(eq["xv"], eq["Gamma"], [0.0, vp["period"] / 4])
    assert np.allclose(Xq[-1], [[0.0, 0.0], [-0.5, 0.5]], atol=1e-9)  # a quarter turn counterclockwise
    # unequal pair Γ₂ = 3Γ₁: G fixed ¾h from vortex 1; rate (Γ₁ + Γ₂)/2πh²
    un = ch05.point_vortex_preset("unequal_pair")
    vu = ch05.vortex_pair(1.0, 3.0, 1.0)
    assert vu["centre_from_1"] == 0.75 and vu["rotation_rate"] == pytest.approx(4 / (2 * np.pi), rel=1e-14)
    Xu = ch05.point_vortex_evolve(un["xv"], un["Gamma"], np.linspace(0, un["t_end"], 50))
    cg = np.array([ch05.centre_of_vorticity(Xk, un["Gamma"]) for Xk in Xu])
    assert np.allclose(cg, [0.25, 0.0], atol=1e-10)  # −0.5 + 0.75
    ang = np.unwrap(np.arctan2(Xu[:, 1, 1] - Xu[:, 1, 0], Xu[:, 0, 1] - Xu[:, 0, 0]))
    assert (ang[-1] - ang[0]) / un["t_end"] == pytest.approx(vu["rotation_rate"], rel=1e-9)
    # opposite pair: translates at Γ/2πh in +y (both), no rotation; G at infinity
    op = ch05.point_vortex_preset("opposite_pair")
    Xo = ch05.point_vortex_evolve(op["xv"], op["Gamma"], [0.0, 10.0])
    assert np.allclose(Xo[-1] - op["xv"], [[0, 0], [1.5915494, 1.5915494]], atol=1e-7)
    assert ch05.vortex_pair(1.0, -1.0, 1.0)["translation_speed"] == pytest.approx(0.1591549, rel=1e-6)
    assert np.all(np.isnan(ch05.centre_of_vorticity(op["xv"], op["Gamma"])))
    assert ch05.point_vortex_invariants(np.array([[0.0, 0.0], [1.0, 0.0]]), [1, 1])["H"] == 0.0
    # the caption's slip: G is not a stagnation point unless Γ₁ = Γ₂ (fluid at G moves at −1.70 m/s for Γ₂ = 3Γ₁)
    uG = ch05.point_vortex_velocity(np.array([0.25, 0.0]), un["xv"], un["Gamma"])
    assert uG[1] == pytest.approx(1 / (2 * np.pi * 0.75) - 3 / (2 * np.pi * 0.25), rel=1e-12) and uG[1] < -1.6
    assert np.allclose(ch05.point_vortex_velocity(np.zeros(2), eq["xv"], eq["Gamma"]), 0.0, atol=1e-16)
    for name in ch05.POINT_VORTEX_PRESETS:
        pp = ch05.point_vortex_preset(name)
        assert pp["xv"].shape == (2, len(pp["Gamma"])) and pp["t_end"] > 0 and pp["label"]
    with pytest.raises(ValueError):
        ch05.point_vortex_preset("nope")


@needs_ref
def test_point_vortex_pairs_V5_garcia_haziot_rates():  # V5 (García & Haziot 2023, §2.1)
    r = ref_json()["point_vortex_pairs"]
    l = 0.37
    ls = sp.Symbol("l", positive=True)
    Om0 = float(sp.sympify(r["corotating_Omega0"], locals={"l": ls}).subs(ls, l))
    V0c = complex(sp.sympify(r["counter_rotating_V0"], locals={"l": ls, "i": sp.I}).subs(ls, l))
    assert V0c.real == 0.0 and V0c.imag < 0  # along −y for +1 at +l, −1 at −l
    V0 = abs(V0c)
    xv = np.array([[l, -l], [0.0, 0.0]])  # z1 = l, z2 = −l
    assert ch05.vortex_pair(1.0, 1.0, 2 * l)["rotation_rate"] == pytest.approx(Om0, rel=1e-14)
    X = ch05.point_vortex_evolve(xv, [1.0, 1.0], [0.0, 1.0])
    assert np.arctan2(X[-1, 1, 0], X[-1, 0, 0]) == pytest.approx(Om0 * 1.0, rel=1e-9)  # counterclockwise
    Xc = ch05.point_vortex_evolve(xv, [1.0, -1.0], [0.0, 2.0])
    assert np.allclose(Xc[-1] - xv, [[0, 0], [-2 * V0, -2 * V0]], atol=1e-9)  # V₀ = −i/(4πl): along −y
    assert ch05.vortex_pair(1.0, -1.0, 2 * l)["translation_speed"] == pytest.approx(V0, rel=1e-14)


# =====================================================================================================================
# C13 — method of images: a vortex near a wall; circle and channel images; rings near a wall, leap-frogging (§5.7); D22
# =====================================================================================================================
def test_wall_image_V1_no_penetration_and_drift():  # V1 (Fig. 5.14; D22: Γ/4πh along the wall)
    xv, Gm = np.array([[0.2], [0.5]]), np.array([1.0])
    allv, allg = ch05.wall_image_system(xv, Gm, 0.0)
    assert np.allclose(allv[:, 1], [0.2, -0.5]) and allg[1] == -1.0
    xw = np.stack([np.linspace(-5, 5, 200), np.zeros(200)])
    uw = ch05.point_vortex_velocity(xw, allv, allg)
    assert np.max(np.abs(uw[1])) < 1e-15  # u·n = 0 on the wall
    # the wall velocity Γh/π(x² + h²) (twice the single vortex's tangential part, our derivation)
    assert np.allclose(uw[0], 1.0 * 0.5 / (np.pi * ((xw[0] - 0.2) ** 2 + 0.25)), rtol=1e-12)
    assert ch05.vortex_near_wall_speed(1.0, 0.5) == pytest.approx(1 / (4 * np.pi * 0.5), rel=1e-15)
    nw = ch05.point_vortex_preset("near_wall")
    X = ch05.point_vortex_evolve(nw["xv"], nw["Gamma"], [0.0, 10.0], boundary=nw["boundary"], **nw["bp"])
    assert X[-1, 0, 0] == pytest.approx(10 * 0.1591549, rel=1e-6) and X[-1, 1, 0] == pytest.approx(0.5, abs=1e-10)
    assert ch05.image_vortices(xv, Gm, "wall", wall_y=0.0)[1][0] == -1.0
    assert ch05.image_vortices(xv, Gm, None)[0].shape == (2, 0)
    # wrong variant: a same-sign image makes fluid cross the wall
    ubad = ch05.point_vortex_velocity(xw, allv, np.abs(allg))
    assert np.max(np.abs(ubad[1])) > 0.1


def test_circle_images_V1_circle_is_a_streamline_and_knife_pair():  # V1 (Fig. 5.13; Exercise 5.14 geometry)
    c, a = np.array([0.3, 0.0]), 1.0
    xv, Gm = np.array([[0.0, 0.1], [0.15, -0.4]]), np.array([-1.0, 1.5])
    th = 2 * np.pi * np.arange(400) / 400
    ring = c[:, None] + a * np.stack([np.cos(th), np.sin(th)])
    for inside in (True, False):
        allv, allg = ch05.circle_image_system(xv if inside else 3 * xv + 2.0, Gm, a, c, inside=inside)
        un = np.sum(ch05.point_vortex_velocity(ring, allv, allg) * (ring - c[:, None]) / a, axis=0)
        assert np.max(np.abs(un)) < 1e-13  # the circle is a streamline
    # outside: net circulation round the cylinder equals cylinder_circulation (0 by default)
    ov, og = ch05.circle_image_system(3 * xv + 2.0, Gm, a, c, inside=False, cylinder_circulation=0.7)
    uu = swirl_field(lambda r: 0 * r)  # placeholder avoided: compute circulation on a loop hugging the cylinder
    loop = ch05.circle_loop_points(c, 1.001 * a, 512)
    u_img = lambda X, t=0.0: ch05.point_vortex_velocity(X, ov[:, 2:], og[2:])  # noqa: E731  images + centre only
    assert ch05.loop_circulation(u_img, loop) == pytest.approx(0.7, rel=1e-9) and uu is not None
    # the knife-blade pair inside the bucket: moves in −x first, stays inside, then separates along the wall
    kb = ch05.point_vortex_preset("knife_bucket")
    t = np.linspace(0.0, kb["t_end"], 241)
    X = ch05.point_vortex_evolve(kb["xv"], kb["Gamma"], t, boundary=kb["boundary"], **kb["bp"])
    assert X[1, 0, 0] < X[0, 0, 0] and X[1, 0, 1] < X[0, 0, 1]
    rr = np.linalg.norm(X - np.array(kb["bp"]["center"])[None, :, None], axis=1)
    assert np.max(rr) < kb["bp"]["a"]
    sep = np.linalg.norm(X[:, :, 0] - X[:, :, 1], axis=1)
    assert sep[-1] > 3 * sep[0]  # separated along the wall


def test_channel_images_V1_series_closed_form_and_rhs():  # V1 (Exercise 5.13 series; our cot closed form)
    for (h, H) in ((0.25, 1.0), (0.1, 1.0), (0.4, 2.0), (0.9, 1.0)):
        s = ch05.channel_image_velocity(h, H, 1.0)
        assert s == pytest.approx(1.0 / (4 * H) / np.tan(np.pi * h / H), rel=1e-9, abs=1e-12)
        assert ch05.channel_image_velocity_exact(h, H, 1.0) == pytest.approx(s, rel=1e-9, abs=1e-12)
    assert abs(ch05.channel_image_velocity(0.5, 1.0, 1.0)) < 1e-12  # mid-channel: 0
    assert ch05.channel_image_velocity(1e-3, 1.0, 1.0) == pytest.approx(1 / (4 * np.pi * 1e-3), rel=1e-5)  # → wall
    # direct sum over ±20 000 image pairs (independent): vortex +Γ at y = h, images at ±2nH + h (+Γ), ±2nH − h (−Γ)
    h, H = 0.25, 1.0
    n = np.arange(1, 20001)
    ys = np.concatenate([2 * n * H + h, -2 * n * H + h])
    ym = np.concatenate([[-h], 2 * n * H - h, -2 * n * H - h])
    u = -np.sum(1.0 / (2 * np.pi) * (h - ys) / (h - ys) ** 2) + np.sum(1.0 / (2 * np.pi) * (h - ym) / (h - ym) ** 2)
    assert u == pytest.approx(0.25, rel=1e-4)
    # the closed-form channel kernel inside point_vortex_rhs moves a vortex at (Γ/4H)cot(πh/H)
    r = ch05.point_vortex_rhs(np.array([[0.0], [0.25]]), [1.0], "channel", H=1.0)
    assert r[0, 0] == pytest.approx(0.25, rel=1e-12) and abs(r[1, 0]) < 1e-14
    X = ch05.point_vortex_evolve(np.array([[0.0], [0.3]]), [1.0], [0.0, 2.0], boundary="channel", H=1.0)
    assert X[-1, 0, 0] == pytest.approx(2.0 / 4.0 / np.tan(0.3 * np.pi), rel=1e-9)
    with pytest.raises(ValueError):
        ch05.image_vortices(np.zeros((2, 1)), [1.0], "channel", H=1.0)


def test_ring_dynamics_V4_leapfrog_and_ring_toward_a_wall():  # V4/V1 (N42–N43, Fig. 5.15)
    t = np.linspace(0.0, 30.0, 601)
    d = ch05.ring_dynamics([dict(R=1.0, z=0.0, Gamma=1.0, a=0.1), dict(R=1.0, z=0.5, Gamma=1.0, a=0.1)], t)
    assert np.ptp(d["impulse"]) < 1e-7 * d["impulse"][0]  # ΣΓπR² conserved (no wall)
    gap = d["z"][:, 1] - d["z"][:, 0]
    assert np.sum(np.diff(np.sign(gap)) != 0) >= 4  # the rings pass through each other repeatedly
    assert np.allclose(d["a"] ** 2 * d["R"], 0.01, rtol=1e-12)  # core volume a²R conserved
    # a single ring alone moves at Kelvin's speed (no wall, no partner)
    d1 = ch05.ring_dynamics([dict(R=1.0, z=0.0, Gamma=1.0, a=0.1)], [0.0, 2.0])
    assert d1["z"][-1, 0] == pytest.approx(2.0 * ch05.ring_self_velocity(1.0, 0.1, 1.0), rel=1e-9)
    # toward a wall at z = 3: the ring widens monotonically and slows down (Fig. 5.15)
    tw = np.linspace(0.0, 60.0, 301)
    w = ch05.ring_dynamics([dict(R=1.0, z=0.0, Gamma=1.0, a=0.1)], tw, wall_z=3.0)
    R, z = w["R"][:, 0], w["z"][:, 0]
    assert np.all(np.diff(R) > 0) and np.all(z < 3.0)
    vz = np.diff(z) / np.diff(tw)
    assert np.all(np.diff(vz) < 1e-12) and vz[-1] < 0.1 * vz[0]


# =====================================================================================================================
# C14 — the vortex sheet: strength = jump in tangential velocity (§5.8); N44, N60; D23
# =====================================================================================================================
def test_vortex_sheet_strength_V1_text_convention_and_caption_variant():  # V1 (D23; the caption's sign WV)
    assert ch05.vortex_sheet_strength(-1.0, 1.0) == 2.0
    assert ch05.vortex_sheet_strength(-1.0, 1.0, convention="caption") == -2.0
    assert ch05.vortex_sheet_strength(-1.0, 1.0, convention="cw") == -2.0
    with pytest.raises(ValueError):
        ch05.vortex_sheet_strength(0.0, 1.0, convention="nope")
    # a row of counterclockwise filaments (γ > 0): u = −γ/2 above, +γ/2 below — the text's u₂ − u₁ returns +γ
    ua = ch05.discrete_sheet_u(0.0, 0.02, 2.0, 2000, L=10.0)  # y = 4 filament spacings above the row
    ub = ch05.discrete_sheet_u(0.0, -0.02, 2.0, 2000, L=10.0)
    uc = ch05.continuous_sheet_velocity(0.0, 0.02, 2.0, L=10.0)[0]
    assert ua == pytest.approx(uc, rel=1e-6) and ub == pytest.approx(-uc, rel=1e-6) and uc == pytest.approx(-1.0,
                                                                                                           rel=3e-3)
    assert ch05.vortex_sheet_strength(ua, ub) == pytest.approx(2.0, rel=3e-3)
    # the circulation of a thin counterclockwise box round the sheet equals γ ds (D23 by direct loop integration)
    xs = np.linspace(-5, 5, 2000, endpoint=False) + 2.5e-3
    u = lambda X, t=0.0: ch05.vortex_sheet_velocity(X, np.stack([xs, 0 * xs]), 2.0, 10.0 / 2000)  # noqa: E731
    box = ch05.square_loop_points((0.0, 0.0), 0.2, 2048) * np.array([[1.0], [0.25]])  # 0.2 × 0.05 rectangle
    assert ch05.loop_circulation(u, box) == pytest.approx(2.0 * 0.2, rel=5e-3)


def test_discrete_sheet_V3_first_order_in_N_and_continuous_form():  # V3/V1 (N60 convergence; the continuous sheet)
    d = ch05.discrete_sheet_convergence(2.0)
    assert d["l1_error"] == pytest.approx([0.043968, 0.0044131, 0.00044128], rel=1e-4)
    assert abs(observed_order([1 / n for n in d["N"]], d["l1_error"]) - 1.0) < ORDER_TOL
    assert np.allclose(d["u_above"], -d["u_below"])
    # the continuous finite sheet: our closed form equals the quad of the plane kernel along it (independent)
    for (x, y) in ((0.0, 0.05), (0.3, -0.1), (0.7, 0.2)):
        pk = [x] if -0.5 < x < 0.5 else None  # the kernel peaks under the field point
        uref = quad(lambda s: -2.0 / (2 * np.pi) * y / ((x - s) ** 2 + y ** 2), -0.5, 0.5, epsabs=0, epsrel=1e-12,
                    points=pk, limit=200)[0]
        vref = quad(lambda s: 2.0 / (2 * np.pi) * (x - s) / ((x - s) ** 2 + y ** 2), -0.5, 0.5, epsabs=1e-13,
                    epsrel=1e-11, points=pk, limit=200)[0]  # v = 0 at x = 0 by symmetry: an absolute floor
        assert ch05.continuous_sheet_velocity(x, y, 2.0) == pytest.approx((uref, vref), rel=1e-10, abs=1e-13)
    assert ch05.continuous_sheet_velocity(0.0, 0.05, 2.0)[0] == pytest.approx(-0.936549, rel=1e-6)
    assert ch05.discrete_sheet_u(0.0, 0.05, 2.0, 100) == pytest.approx(-0.936551, rel=1e-6)
    ys = np.array([0.05, -0.05])
    assert np.allclose(ch05.discrete_sheet_u(0.0, ys, 2.0, 100), [-0.936551, 0.936551], rtol=1e-6)
    S = np.stack([np.linspace(-0.5, 0.5, 11), np.zeros(11)])
    P = RNG.uniform(-1, 1, (2, 4)) + np.array([[0.0], [0.3]])
    assert np.allclose(ch05.vortex_sheet_velocity(P, S, 2.0, 0.1), ch05.point_vortex_velocity(P, S, 0.2 * np.ones(11)))


def test_sheet_rollup_V1_periodic_kernel_V4_centroid():  # V1/V4 (the Krasny kernel = the periodic image sum)
    N, L = 16, 1.0
    d = ch05.sheet_rollup(N=N, gamma=1.0, amplitude=0.05, delta=0.0, t_eval=[0.0, 1e-7], L=L)
    vel = np.stack([(d["x"][1] - d["x"][0]) / 1e-7, (d["y"][1] - d["y"][0]) / 1e-7])
    X0 = np.stack([d["x"][0], d["y"][0]])
    # (a) the complex form the docstring names: u − iv = Σ (Γ_j/2iL) cot(π(z − z_j)/L), self term excluded
    zz = X0[0] + 1j * X0[1]
    dz = zz[:, None] - zz[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        cotm = np.where(np.eye(N, dtype=bool), 0.0, 1.0 / np.tan(np.pi * dz / L))
    w = np.sum((L / N) / (2j * L) * cotm, axis=1)
    assert np.allclose(vel, np.stack([w.real, -w.imag]), rtol=1e-6, atol=1e-6 * np.max(np.abs(w)))  # Δt = 1e-7 s
    # (b) the direct sum over ±M periodic copies (its symmetric tail is O(1/M): 1e-3 at M = 3000)
    ref = np.zeros((2, N))
    for m in range(-3000, 3001):
        ref += ch05.point_vortex_velocity(X0, X0 + np.array([[m * L], [0.0]]), np.full(N, L / N))
    assert np.allclose(vel, ref, rtol=1e-2, atol=2e-4 * np.max(np.abs(ref)))
    r = ch05.sheet_rollup(N=100)
    assert np.ptp(r["mean_y"]) < 1e-15 and r["x"].shape == (9, 100) and r["Gamma_each"] == pytest.approx(0.01)
    # the perturbation grows (Kelvin–Helmholtz): the displacement at t = 4L/γ exceeds the initial amplitude tenfold
    assert np.max(np.abs(r["y"][-1])) > 10 * 0.01


def test_cylinder_as_vortex_sheet_V1_ring_of_point_vortices():  # V1 (Exercise 5.19 geometry: sheet on a circle)
    Om, a, N = 0.7, 0.4, 2000
    th = 2 * np.pi * (np.arange(N) + 0.5) / N
    S = a * np.stack([np.cos(th), np.sin(th)])
    gam = Om * a  # the sheet strength equal to the rim speed of the rigid patch
    ds = 2 * np.pi * a / N
    for r in (0.6, 1.0, 2.0):
        u = ch05.vortex_sheet_velocity(np.array([r, 0.0]), S, gam, ds)
        assert u[1] == pytest.approx(Om * a ** 2 / r, rel=1e-9)  # outside: the patch's irrotational field Ωa²/x
    assert np.linalg.norm(ch05.vortex_sheet_velocity(np.array([0.1, 0.05]), S, gam, ds)) < 1e-9  # inside: at rest
    # the equivalent vorticity patch: uniform 2Ω inside (a Rankine core of Γ = 2πΩa²)
    u_r, w_r = ch05.rotating_cylinder_flow(np.array([0.2, 1.0]), a, 2 * Om)
    assert w_r[0] == pytest.approx(2 * Om) and u_r[1] == pytest.approx(Om * a ** 2 / 1.0, rel=1e-14)


# =====================================================================================================================
# The contract with the notebook and explainers (design Part C), scalar parity, drawing helpers
# =====================================================================================================================
PART_C = """
vorticity_field vortex_line TubeStrength vortex_tube_strength Tube vortex_tube TubeFlux tube_flux_budget material_loop
loop_circulation loop_length material_circulation KelvinRate kelvin_rate_terms KelvinForces kelvin_force_terms
kelvin_scenario kelvin_scenario_circulation kelvin_scenario_rate pressure_torque_on_element baroclinic_term
baroclinic_term_sym frozen_in_check vorticity_equation_sym VorticityTerms vorticity_terms VorticityBudget
vorticity_budget vorticity_budget_sym vorticity_budget_preset vorticity_divergence rotating_lamb_form_terms
stretching_tilting_split planetary_vorticity_terms loop_vector_area absolute_circulation
poisson_green_3d velocity_from_curl_omega gaussian_tube_fields velocity_from_vorticity_fft biot_savart_volume
biot_savart_2d segment_induced_velocity segment_speed filament_velocity filament_contributions filament_preset
filament_velocity_preset ring_axis_velocity ring_ring_velocity ring_self_velocity point_vortex_velocity
point_vortex_rhs point_vortex_evolve point_vortex_invariants centre_of_vorticity point_vortex_preset
wall_image_system circle_image_system channel_image_velocity vortex_sheet_velocity continuous_sheet_velocity
burgers_vortex burgers_vortex_field burgers_core_radius hill_spherical_vortex hill_spherical_vortex_field
lamb_oseen_field gaussian_tube_vorticity gaussian_tube_field gaussian_tube_section broken_tube_field
vortex_ring_vorticity CurlCheck curl_theorem_box
solid_body_from_vorticity line_vortex_gamma solid_body_pressure_gradients solid_body_pressure line_vortex_pressure
rankine_pressure isobar_height rotating_tank_free_surface bernoulli_across_vortex polar_viscous_stress
line_vortex_viscous_stress vortex_stress_force rotating_cylinder_flow torque_per_length dissipation_outside_cylinder
vortex_pressure_scenario lamb_oseen_circulation lamb_oseen_viscous_loop_integral kelvin_hypotheses
kelvin_hypotheses_text lock_exchange_initial_vorticity_rate lock_exchange_fields baroclinic_rate_2d
rotating_ns_residual absolute_vorticity diffusing_vortex_sheet strain_preset uniform_strain_vorticity stretched_tube
helix helix_frame column_relative_vorticity column_over_slope relative_circulation_after_move helical_swirl_field
cellular_flow abc_flow vortex_pair vortex_near_wall_speed ring_dynamics vortex_sheet_strength discrete_sheet_u
discrete_sheet_convergence sheet_rollup
""".split()
# public additions beyond Part C that the scripts/notebook use (also exercised here)
EXTRA = """
tube_flux_budget_traced planar_field_3d stencil_gradient_fn circle_loop_points square_loop_points loop_tangent
loop_line_integral vorticity_terms_sym KELVIN_SCENARIOS VORTICITY_BUDGET_PRESETS cylinder_quadrature
channel_image_velocity_exact image_vortices RING_CORES FILAMENT_PRESETS FILAMENT_CLOSED POINT_VORTEX_PRESETS
stretched_gaussian_vortex_field tube_core_radius polar_net_viscous_force stress_vs_net_force_table tornado_circulation
rankine_radius_at_pressure kelvin_scenario_gamma burgers_balance stretching_tilting_scenario frenet_frame
helical_vortex_line baroclinic_element_scenario STRAIN_PRESETS
""".split()


def test_contract_V1_every_part_c_name_exists_and_is_reexported():  # V1 (design Part C)
    missing = [n for n in PART_C + EXTRA if not hasattr(ch05, n)]
    assert not missing, missing
    assert len(PART_C) == 119
    # re-exports are the core objects, not copies
    assert ch05.vorticity_terms is VD.vorticity_terms and ch05.ring_ring_velocity is BS.ring_ring_velocity
    assert ch05.burgers_vortex is VX.burgers_vortex and ch05.curl_theorem_box is IT.curl_theorem_box
    assert ch05.abc_flow(1, 1, 1)(np.ones(3)) == pytest.approx(VX.abc_flow_field(1, 1, 1)(np.ones(3)))
    # contracted leading parameters (design Part C signatures)
    sig = lambda n, k=3: list(inspect.signature(getattr(ch05, n)).parameters)[:k]  # noqa: E731
    assert sig("tube_flux_budget", 4) == ["field", "z_a", "z_b", "R0"]
    assert sig("pressure_torque_on_element", 4) == ["p", "rho", "center", "radius"]
    assert sig("vortex_pressure_scenario", 2) == ["kind", "r"]
    assert sig("velocity_from_curl_omega", 5) == ["curl_omega_fn", "x", "nodes", "weights", "sign"]
    assert sig("ring_ring_velocity", 5) == ["R", "z", "R0", "z0", "Gamma0"]
    assert sig("column_relative_vorticity", 4) == ["h", "h0", "zeta0", "f"]
    assert sig("relative_circulation_after_move", 5) == ["Gamma0", "A0", "lat0_deg", "A1", "lat1_deg"]
    assert sig("kelvin_hypotheses", 4) == ["inviscid", "barotropic", "conservative", "inertial"]
    assert sig("vortex_sheet_strength", 3) == ["u_above", "u_below", "convention"]


def test_contract_V1_every_part_c_name_is_exercised_in_this_file():  # V1
    full = Path(__file__).read_text(encoding="utf-8")
    start = full.index('PART_C = """')
    end = full.index('def test_contract_V1_every_part_c_name_exists_and_is_reexported')
    src = full[:start] + full[end:]  # everything except the name lists themselves
    unused = [n for n in PART_C + EXTRA if not re.search(rf"ch05\.{n}\b", src)]
    assert not unused, unused


def test_result_types_V1_named_fields():  # V1 (the NamedTuples the notebook unpacks)
    assert ch05.TubeStrength._fields == ("circulation", "flux", "diff")
    assert ch05.TubeFlux._fields == ("lower", "side", "upper", "total")
    assert ch05.KelvinRate._fields == ("acceleration", "contour", "total")
    assert ch05.KelvinForces._fields == ("pressure", "body", "viscous")
    assert ch05.VorticityTerms._fields == ("local", "advective", "stretching_tilting", "diffusion", "residual")
    assert ch05.VorticityBudget._fields == ("local", "advective", "stretching_tilting", "planetary", "baroclinic",
                                            "diffusion", "residual")
    assert ch05.CurlCheck._fields == ("volume", "surface", "diff")
    tube = ch05.Tube([np.zeros((3, 2))], np.zeros(3), np.array([0, 0, 1.0]), 0.1)
    assert tube.array.shape == (1, 3, 2) and tube.s.size == 0
    kf = ch05.KelvinForces(1.0, 2.0, 3.0)
    assert kf.total == 6.0
    # the cellular flow function equals the scenario's field; ω = 2ψ/ℓ² (a steady Euler flow)
    s = ch05.kelvin_scenario("cellular")
    P = RNG.uniform(0, 3, (2, 7))
    assert np.allclose(ch05.cellular_flow(P), s["u"](P))
    wz = ch05.vorticity_field(lambda X, t=0.0: ch05.cellular_flow(X), 1e-4)(P)[0]
    assert np.allclose(wz, 2 * np.sin(P[0]) * np.sin(P[1]), rtol=1e-7)


def test_scalar_callable_V1_explainer_parity_functions_return_floats():  # V1 (selftest parity rows need floats)
    vals = [ch05.gaussian_tube_section(0.5)["flux"], ch05.solid_body_pressure(0.1, 0.0, 10.0),
            ch05.line_vortex_pressure(0.1, 0.0, 1.0), ch05.rankine_pressure(0.05, 0.0, 1.0, 0.1),
            ch05.isobar_height(0.1, 10.0), ch05.bernoulli_across_vortex("solid", 0.1, param=10.0),
            ch05.line_vortex_viscous_stress(0.1, 1.0), ch05.torque_per_length(0.3, 1.0),
            ch05.lamb_oseen_circulation(0.005, 10.0, 0.01, 1e-6)[0], ch05.kelvin_scenario_circulation("rotating", 5.0),
            ch05.lock_exchange_initial_vorticity_rate(1000.0, 1025.0, 0.1), ch05.baroclinic_rate_2d([10, 0], [0, -9810],
                                                                                                    1000),
            ch05.burgers_core_radius(1.0, 1e-6), ch05.burgers_vortex(1e-3, 0.0, 1e-3, 1.0, 1e-6)[1],
            ch05.ring_axis_velocity(0.2, 0.5, 1.0), ch05.ring_self_velocity(1.0, 0.1, 1.0),
            ch05.segment_speed(1.0, 0.3, 2.0, 1.0), ch05.filament_velocity_preset("ring", 0, 0, 0, component=2),
            ch05.column_relative_vorticity(1100.0, 1000.0), ch05.relative_circulation_after_move(0, 1e10, 30, 1e10, 60),
            ch05.vortex_near_wall_speed(1.0, 0.5), ch05.channel_image_velocity(0.25, 1.0, 1.0),
            ch05.discrete_sheet_u(0.0, 0.05, 2.0, 100), ch05.vortex_sheet_strength(-1.0, 1.0),
            ch05.vortex_pair(1, 3, 1)["centre_from_1"], ch05.ring_ring_velocity(0.3, 0.2, 0.5, 0.0, 1.0)[1],
            ch05.kelvin_scenario_rate("baroclinic")["pressure"], ch05.vorticity_budget_preset("burgers")["advective"],
            ch05.vortex_pressure_scenario("line", 0.2)["p"], ch05.continuous_sheet_velocity(0.0, 0.1, 1.0)[0],
            ch05.lamb_oseen_viscous_loop_integral(0.005, 10.0, 0.01, 1e-6), ch05.tornado_circulation(-2000.0, 15.0),
            ch05.poisson_green_3d(np.ones(3), np.zeros(3)), ch05.channel_image_velocity_exact(0.25, 1.0, 1.0)]
    assert all(isinstance(v, float) for v in vals), [type(v) for v in vals]
    assert isinstance(ch05.kelvin_hypotheses_text(True, False), str)


def test_scripts_V1_drawing_helpers_run():  # V1 smoke (C.6)
    import sys

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(ROOT / "scripts"))
    from scripts.ch05_drawings import (biot_savart_geometry, column_sketch, disc_element, helix_with_frame,  # noqa
                                       loop_with_element, sheet_circuit, tank_section, tube_mesh, vortex_marks)
    fig = plt.figure()
    ax = fig.add_subplot(121)
    ax3 = fig.add_subplot(122, projection="3d")
    tube = ch05.vortex_tube("gaussian_tube", (0, 0, 0), (0, 0, 1), 0.1, n_lines=6, s_max=0.5, n=21)
    tube_mesh(ax3, tube)
    biot_savart_geometry(ax3, ch05.filament_preset("ring", 16), np.array([0.2, 0.1, 0.4]), closed=True)
    tank_section(ax, 0.1, 0.2)
    disc_element(ax, (0.0, 0.0), 0.1)
    loop_with_element(ax, ch05.circle_loop_points((0, 0), 0.5, 64), 3, u=ch05.cellular_flow)
    x = np.linspace(-1e5, 1e5, 11)
    c = ch05.column_over_slope(x)
    column_sketch(ax, c["x"], c["h"], c["zeta"])
    vortex_marks(ax, np.array([[0.0, 1.0], [0.0, 0.0]]), [1.0, -1.0])
    sheet_circuit(ax)
    pf = helix_with_frame()
    assert pf is not None and len(ax.lines) + len(ax.patches) + len(ax.collections) > 0
    plt.close("all")


def test_scripts_V1_every_ch05_script_runs(tmp_path):  # V1 smoke (all nine scripts exit 0; figures written)
    import os
    import subprocess
    import sys
    env = dict(os.environ, MPLBACKEND="Agg")
    for scr in sorted((ROOT / "scripts").glob("ch05_*.py")):
        args = [sys.executable, str(scr)] + (["--no-show"] if scr.name == "ch05_drawings.py" else [])
        r = subprocess.run(args, cwd=ROOT, env=env, capture_output=True, text=True, timeout=600)
        assert r.returncode == 0, (scr.name, r.stderr[-2000:])


# =====================================================================================================================
# V6 — the book's printed forms and exercise numbers (private JSON; skipped when absent)
# =====================================================================================================================
@book_only
def test_book_V6_printed_forms_of_section_5_1_and_5_7():  # V6
    b = book()
    s = {n: sp.Symbol(n, positive=True) for n in ("omega", "a", "r", "rho", "g", "z", "Gamma", "mu", "h", "Gamma1",
                                                   "Gamma2")}
    cyl = b["sec5_1_rotating_cylinder"]
    for r_ in (0.2, 0.7):
        ref = float(sp.sympify(cyl["u_theta_outside"], locals=s).subs({s["omega"]: 3.0, s["a"]: 0.15, s["r"]: r_}))
        assert ch05.rotating_cylinder_flow(r_, 0.15, 3.0)[0] == pytest.approx(ref, rel=1e-14)
    assert float(sp.sympify(cyl["rankine_Gamma"], locals=s).subs({s["a"]: 0.15, s["omega"]: 3.0})) == pytest.approx(
        np.pi * 0.15 ** 2 * 3.0, rel=1e-15)
    pf = b["sec5_1_pressure_forms"]
    vals = {s["rho"]: 998.0, s["omega"]: 7.0, s["r"]: 0.3, s["g"]: 9.81, s["z"]: -0.2, s["Gamma"]: 1.3, s["mu"]: 1e-3}
    assert ch05.solid_body_pressure(0.3, -0.2, 7.0, 998.0, 9.81) == pytest.approx(
        float(sp.sympify(pf["solid_body"], locals=s).subs(vals)), rel=1e-14)
    assert ch05.line_vortex_pressure(0.3, -0.2, 1.3, 998.0, 9.81) == pytest.approx(
        float(sp.sympify(pf["line_vortex"], locals=s).subs(vals)), rel=1e-14)
    assert ch05.line_vortex_viscous_stress(0.3, 1.3, 1e-3) == pytest.approx(
        float(sp.sympify(pf["sigma_r_theta_line_vortex"], locals=s).subs(vals)), rel=1e-14)
    pw = b["sec5_7_pairs_and_wall"]
    pv = {s["Gamma1"]: 1.2, s["Gamma2"]: 0.4, s["h"]: 0.8, s["Gamma"]: 1.2}
    vp = ch05.vortex_pair(1.2, 0.4, 0.8)
    assert vp["V1"] == pytest.approx(float(sp.sympify(pw["V1"], locals=s).subs(pv)), rel=1e-15)
    assert vp["V2"] == pytest.approx(float(sp.sympify(pw["V2"], locals=s).subs(pv)), rel=1e-15)
    assert ch05.vortex_pair(1.2, -1.2, 0.8)["translation_speed"] == pytest.approx(
        float(sp.sympify(pw["opposite_pair_speed"], locals=s).subs(pv)), rel=1e-15)
    assert ch05.vortex_near_wall_speed(1.2, 0.8) == pytest.approx(
        float(sp.sympify(pw["wall_image_speed"], locals=s).subs(pv)), rel=1e-15)
    sh = b["sec5_8_sheet"]
    u1, u2 = sp.symbols("u1 u2")
    assert float(sp.sympify(sh["text_strength"]).subs({u1: -0.3, u2: 0.9})) == ch05.vortex_sheet_strength(-0.3, 0.9)
    assert float(sp.sympify(sh["caption_strength"]).subs({u1: -0.3, u2: 0.9})) == ch05.vortex_sheet_strength(
        -0.3, 0.9, convention="caption")


@book_only
def test_book_V6_exercise_5_1_rotating_tank_uncovered_area():  # V6 (reproduces the exercise's area to < 0.5 %)
    e = book()["ex5_1_rotating_tank"]
    t = ch05.rotating_tank_free_surface(e["tank_diameter_m"] / 2, e["water_depth_m"], e["omega0_rad_s"],
                                        H=e["tank_height_m"], closed=True)
    assert t["uncovered_area"] == pytest.approx(e["uncovered_area_m2"], rel=5e-3)
    # the exercise's hint height ω₀²r²/2g is the (5.6) isobar with the vorticity ω = 2ω₀
    om0, r, g = sp.symbols("omega0 r g")
    hint = float(sp.sympify(e["hint_height"]).subs({om0: 40.0, r: 0.3, g: 9.81}))
    assert hint == pytest.approx(ch05.isobar_height(0.3, 80.0) - ch05.isobar_height(0.0, 80.0), rel=1e-14)


@book_only
def test_book_V6_exercise_5_2_tornado():  # V6 (Γ from the core-edge gauge pressure with ρ(25 °C) = 1.184 kg/m³)
    e = book()["ex5_2_tornado"]
    Gam = ch05.tornado_circulation(e["gauge_pressure_at_core_edge_Pa"], e["radius_m"], rho=1.184)
    assert Gam == pytest.approx(e["circulation_m2_s"], rel=5e-3)
    assert ch05.rankine_pressure(e["radius_m"], 0.0, Gam, e["radius_m"], rho=1.184) == pytest.approx(
        e["gauge_pressure_at_core_edge_Pa"], rel=1e-12)
    r1 = ch05.rankine_radius_at_pressure(e["gauge_from_Pa"], Gam, e["radius_m"])
    r2 = ch05.rankine_radius_at_pressure(e["gauge_to_Pa"], Gam, e["radius_m"])
    assert r2 == pytest.approx(e["radius_m"], rel=1e-9) and r1 > r2
    assert (r1 - r2) / e["translation_speed_m_s"] > 0  # the time the barometer takes (ours)


@book_only
def test_book_V6_exercise_closed_forms():  # V6 (Exercises 5.3, 5.5, 5.6, 5.11, 5.12, 5.13, 5.19, 5.20)
    b = book()
    names = ("a", "R", "z", "rho1", "rho2", "g", "delta", "gamma", "nu", "t", "y", "A", "alpha", "Gamma", "h", "H",
             "x", "Omega_z", "psi")
    s = {n: sp.Symbol(n, positive=True) for n in names}
    # 5.3: helical flow, its vorticity and the vortex-line invariant
    e3 = b["ex5_3_helical_flow"]
    vel = sp.sympify(e3["velocity"], locals=s)
    vort = sp.sympify(e3["vorticity"], locals=s)
    P = np.array([0.3, 0.4, 0.7])
    Rv = np.hypot(P[0], P[1])
    uu = ch05.helical_swirl_field(0.8)(P)
    assert np.hypot(uu[0], uu[1]) == pytest.approx(float(vel[1].subs({s["a"]: 0.8, s["R"]: Rv, s["z"]: 0.7})), rel=1e-14)
    wv = ch05.vorticity_field(ch05.helical_swirl_field(0.8), 1e-5)(P)
    assert wv[2] == pytest.approx(float(vort[2].subs({s["a"]: 0.8, s["z"]: 0.7})), rel=1e-7)
    assert np.hypot(wv[0], wv[1]) == pytest.approx(abs(float(vort[0].subs({s["a"]: 0.8, s["R"]: Rv}))), rel=1e-7)
    # 5.5: lock-exchange rate
    f55 = sp.sympify(b["ex5_5_lock_exchange"]["rate"], locals=s)
    assert ch05.lock_exchange_initial_vorticity_rate(1000.0, 1025.0, 0.1, 9.81) == pytest.approx(
        float(f55.subs({s["rho1"]: 1000, s["rho2"]: 1025, s["g"]: 9.81, s["delta"]: 0.1})), rel=1e-14)
    # 5.6: diffusing sheet
    f56 = sp.sympify(b["ex5_6_diffusing_sheet"]["omega_z"], locals=s)
    assert ch05.diffusing_vortex_sheet(1e-3, 2.0, 1.5, 1e-6)[1] == pytest.approx(
        float(f56.subs({s["gamma"]: 1.5, s["nu"]: 1e-6, s["t"]: 2.0, s["y"]: 1e-3})), rel=1e-13)
    # 5.11: Hill's vortex (inside)
    e11 = b["ex5_11_hill"]
    psi = sp.sympify(e11["psi"], locals=s)
    for (Rv_, zv) in ((0.3, 0.2), (0.5, -0.4)):
        sub = {s["A"]: 1.3, s["a"]: 0.9, s["R"]: Rv_, s["z"]: zv}
        assert VX.hill_stream_function(Rv_, zv, 1.3, 0.9) == pytest.approx(float(psi.subs(sub)), rel=1e-13)
        uR = float((-(1 / s["R"]) * sp.diff(psi, s["z"])).subs(sub))
        uz = float(((1 / s["R"]) * sp.diff(psi, s["R"])).subs(sub))
        got = ch05.hill_spherical_vortex(Rv_, zv, 1.3, 0.9)
        assert (got[0], got[1]) == pytest.approx((uR, uz), rel=1e-13)
        assert got[2] == pytest.approx(float(sp.sympify(e11["omega_phi"], locals=s).subs(sub)), rel=1e-14)
    # 5.12: Burgers' vortex
    e12 = b["ex5_12_burgers"]
    sub = {s["alpha"]: 1.0, s["nu"]: 1e-6, s["Gamma"]: 1e-3, s["R"]: 1.5e-3, s["z"]: 0.2}
    got = ch05.burgers_vortex(1.5e-3, 0.2, 1e-3, 1.0, 1e-6)
    assert got[0] == pytest.approx(float(sp.sympify(e12["u_R"], locals=s).subs(sub)), rel=1e-14)
    assert got[1] == pytest.approx(float(sp.sympify(e12["u_phi"], locals=s).subs(sub)), rel=1e-13)
    assert got[2] == pytest.approx(float(sp.sympify(e12["u_z"], locals=s).subs(sub)), rel=1e-14)
    # 5.13: the channel series (partial sums converge to our value) and the mid-plane answer
    e13 = b["ex5_13_channel"]
    n = np.arange(1, 200001)
    ser = 1.0 / (4 * np.pi * 0.25) * (1 - 2 * np.sum(1.0 / ((n * 1.0 / 0.25) ** 2 - 1)))
    assert ser == pytest.approx(ch05.channel_image_velocity(0.25, 1.0, 1.0), rel=1e-6)
    assert abs(ch05.channel_image_velocity(0.5, 1.0, 1.0) - e13["midplane_answer"]) < 1e-12
    # 5.19: a rotating patch as a sheet on its rim
    e19 = b["ex5_19_cylinder_sheet"]
    Om, a = 0.7, 0.4
    assert ch05.rotating_cylinder_flow(1.3, a, 2 * Om)[0] == pytest.approx(
        float(sp.sympify(e19["u_on_x_axis"], locals=s).subs({s["Omega_z"]: Om, s["a"]: a, s["x"]: 1.3})), rel=1e-14)
    assert float(sp.sympify(e19["patch_vorticity"], locals=s).subs(s["Omega_z"], Om)) == ch05.rotating_cylinder_flow(
        0.1, a, 2 * Om)[1]
    assert float(sp.sympify(e19["sheet_strength"], locals=s).subs({s["Omega_z"]: Om, s["a"]: a})) == pytest.approx(
        ch05.rotating_cylinder_flow(a, a, 2 * Om)[0], rel=1e-14)
    # 5.20: velocity on the wall under a vortex
    e20 = b["ex5_20_wall_image"]
    allv, allg = ch05.wall_image_system(np.array([[0.0], [0.6]]), [1.1])
    xw = np.array([[0.4], [0.0]])
    assert ch05.point_vortex_velocity(xw, allv, allg)[0, 0] == pytest.approx(
        float(sp.sympify(e20["wall_velocity_x"], locals=s).subs({s["Gamma"]: 1.1, s["h"]: 0.6, s["x"]: 0.4})), rel=1e-13)
