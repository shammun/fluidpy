"""Verification suite for Chapter 14 — Aerodynamics (Kundu, Cohen & Dowling 5e, §§14.1–14.9, Eqs. (14.1)–(14.27)).

Evidence levels (``verify-implementation`` skill): V1 analytic · V2 symbolic (sympy / dimensions) · V3 convergence ·
V4 conservation / invariant / integral identity · V5 published benchmark (``reference/ch14/benchmarks.json``, see
SOURCES.md) · V6 book value (private, ``tests/book_values_ch14.json``, skipped when absent) · V7 limits / symmetry /
invariance.  Every test name is ``test_<concept>_<level>_<what>``; the comment on the ``def`` line repeats the level and
the curation ID.

A items (CORE, >= 2 independent levels): C01 aspect ratio (14.1)–(14.2) and level flight · C02 the Kutta condition ·
C03 the starting vortex (Kelvin) · C04 the Zhukhovsky airfoil (14.3)–(14.8) · C05 the flow round it (14.9) · C06 the
Kutta circulation and the lift (14.11)–(14.12) · C07 surface speed and pressure · C08 the shed filament · C09 the
downwash (14.13) · C10 induced drag (14.14)–(14.15) · C11 the sine series and Glauert's integral (14.18)–(14.19) ·
C12 the lifting-line equation (14.16)–(14.20) · C13 lift and induced drag from the coefficients (14.21)–(14.22) ·
C14 the elliptic loading (14.23)–(14.27) · C15 the drag polar · C16 thrust from lift (fish tail, sail, keel).

Derivations: every two- and three-star D row (D03–D06, D08–D10, D12–D22) is re-derived in a ``test_*_V2_derivation``
test with sympy, **independently of the chapter's own sympy engine** (``circular_arc_sym`` is tested separately); the
intermediate lines of the four three-star rows (D03, D04, D15, D18) are checked one by one.

Printed slips kept as wrong variants that must FAIL an independent check (the pages were re-read as images first):
#1 the denominator of the first two members of (14.19) · #2 "epsilon = tan(w/U)" · #4 the thickness coefficient (the
printed 1.3 is 3*sqrt(3)/4 rounded, and the rule is first order in e) · #8 "c = 4a" (the exact lift slope exceeds 2 pi).

All inputs are ours (``ch14.illustrative_inputs()`` or visibly odd numbers); nothing here is quoted from the book.

Run: ``.venv/Scripts/python.exe -m pytest tests/test_ch14.py -q -p no:cacheprovider``  (``-m "not slow"`` skips the
script runs).
"""
from __future__ import annotations

import inspect
import json
import math
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import sympy as sp
from scipy.integrate import dblquad, quad
from scipy.optimize import brentq, minimize_scalar

from fluidpy import ch06_ideal_flow as ch06
from fluidpy import ch14_aerodynamics as ch14
from fluidpy.core import airfoils as AF
from fluidpy.core import biot_savart as BS
from fluidpy.core import conformal as CONF
from fluidpy.core import lifting_line as LL
from fluidpy.core import potential as POT
from fluidpy.core import similarity as SIM
from fluidpy.core.units import Q_
from tools.convergence import observed_order, pairwise_orders

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference" / "ch14"
BOOK = Path(__file__).resolve().parent / "book_values_ch14.json"
book_only = pytest.mark.skipif(not BOOK.exists(), reason="book values are private; see CLAUDE.md rule 9")
slow = pytest.mark.slow
PI = math.pi
I = ch14.illustrative_inputs()
B, E, BETA, ALPHA, US, RHO = I["b"], I["e"], I["beta"], I["alpha"], I["U_section"], I["rho"]   # our section (ours)
SEC = (US, B, E, BETA, ALPHA)                    # positional (U, b, e, beta, alpha) of the section functions
GN = np.array([2.3, 0.0, 0.31, 0.0, -0.12])      # a symmetric three-term loading of ours [m^2/s]
GN_ASYM = np.array([1.9, 0.27, -0.4, 0.11])      # one with even terms (not symmetric about mid-span)
SPAN = 7.0                                       # an odd span of ours [m]


def book():
    return json.loads(BOOK.read_text(encoding="utf-8"))


def bench():
    return json.loads((REF / "benchmarks.json").read_text(encoding="utf-8"))


def rel(a, b):
    return abs(float(a) - float(b)) / max(abs(float(b)), 1e-300)


def close(a, b, rtol=1e-10, atol=0.0):
    """``numpy.allclose`` with NO absolute tolerance unless one is written (numpy's default atol = 1e-8 hides small numbers)."""
    return bool(np.allclose(a, b, rtol=rtol, atol=atol))


def z0(expr):
    """True if a sympy expression is identically zero (several simplifiers tried)."""
    for f in (sp.simplify, lambda q: sp.simplify(sp.expand_trig(q)), lambda q: sp.simplify(sp.expand(q, complex=True)),
              lambda q: sp.simplify(sp.trigsimp(sp.expand(q)))):
        if f(expr) == 0:
            return True
    return False


def pv_symmetric(f, lo, hi, x0, eps=1e-6):
    """Cauchy principal value by symmetric exclusion (independent of the singularity-subtraction used in ``LL``)."""
    a = quad(f, lo, x0 - eps, limit=400, epsabs=1e-13, epsrel=1e-12)[0]
    b = quad(f, x0 + eps, hi, limit=400, epsabs=1e-13, epsrel=1e-12)[0]
    return a + b


# =====================================================================================================================
# 0. Contract: every name of design Part C exists, is re-exported, has the documented shape, and is used in this file
# =====================================================================================================================
PART_C = """
THICKNESS_RULE_COEFFICIENT zhukhovsky_circle zhukhovsky_circle_from_centre zhukhovsky_airfoil zhukhovsky_inverse
section_geometry trailing_edge_angle zhukhovsky_shape circular_arc zero_lift_angle camber_from_beta resolve_lift_drag
section_coefficient circulation_from_cl cl_from_circulation normal_force_coefficient lift_from_cp kutta_circulation
kutta_circulation_root stagnation_points rear_stagnation_point trailing_edge_speed zhukhovsky_flow flat_plate_flow
zhukhovsky_lift_coefficient zhukhovsky_surface surface_circulation zhukhovsky_force flat_plate_forces thin_airfoil
planform aspect_ratio mean_chord span_angle span_position circulation_series trailing_sheet_strength
elliptic_circulation elliptic_downwash downwash downwash_series glauert_integral induced_angle effective_angle
strip_forces lift_and_induced_drag lift_from_series induced_drag_from_series span_efficiency section_lift_coefficient
lifting_line_modes lifting_line_matrix lifting_line_rhs lifting_line_residual solve_lifting_line elliptic_induced_drag
induced_drag_coefficient finite_wing_lift_slope elliptic_wing elliptic_twist_rectangular twist_for_elliptic_loading
horseshoe_velocity horseshoe_downwash horseshoe_wing downwash_along_x ground_effect_factor staircase_from_loading
staircase_downwash rolled_up_pair wake_kinetic_energy_per_length
illustrative_inputs book_slips traps conventions_table compressibility_check level_flight stall_speed surface_point
edge_speed map_point edge_flow_suction trailing_edge_wedge circular_arc_sym suction_side_separation chord_reynolds
stall_regime impulsive_start_lumped_vortex wagner_function wake_rollup tip_vortex_estimate
vortex_pair_energy_per_length wake_vortex_rolling_moment generic_lift_curve generic_drag_curve drag_polar
drag_breakdown drag_share drag_vs_speed speed_min_drag best_glide profile_drag_estimate oscillating_foil apparent_wind
sail_forces sail_drive_curve velocity_made_good sailboat_course_limit sailboat_equilibrium projectile_range_drag
joukowski joukowski_derivative joukowski_inverse circle_flow_zeta lift_per_span cylinder_surface_speed
cylinder_stagnation_points segment_speed segment_induced_velocity point_vortex_velocity point_vortex_evolve
point_vortex_invariants blasius_force pressure_coefficient lift_coefficient drag_coefficient mach_number
reynolds_number
""".split()
FLOW_METHODS = ("w_zeta dwdzeta zeta_of w dwdz_from_zeta dwdz velocity speed stream phi psi inside speed_near_edge "
                "cp").split()
DRAWINGS = ("draw_planform_labelled draw_three_axes draw_section_labelled draw_axes_conventions draw_kelvin_circuit "
            "draw_triangle_OQP draw_circle_angles draw_tip_flow_rear_view draw_horseshoe_staircase draw_force_triangle "
            "draw_span_semicircle draw_stalled_section draw_velocity_triangle draw_sail_polygon fig_four_circles "
            "fig_three_stages fig_pressure_two_ways fig_planforms_loadings").split()


def test_contract_V7_every_part_c_name_exists_and_is_exercised_here():  # V7 — design Part C (C.0–C.3): 127 names
    src = Path(__file__).read_text(encoding="utf-8")
    body = src.split('""".split()', 1)[1]                                # everything after the list itself
    assert len(PART_C) == len(set(PART_C)) == 127
    missing = [n for n in PART_C if not hasattr(ch14, n)]
    assert not missing, missing
    unused = [n for n in PART_C if not re.search(rf"\b(?:ch14|AF|LL)\.{n}\b", body)]
    assert not unused, unused                                            # each name is called by at least one test
    F = ch14.zhukhovsky_flow(*SEC)
    assert all(callable(getattr(F, m)) for m in FLOW_METHODS)
    assert all(re.search(rf"\.{m}\(", body) for m in FLOW_METHODS)
    for attr in ("U", "b", "e", "beta", "alpha", "mu", "a", "Gamma_kutta", "Gamma_cw", "is_kutta"):
        assert hasattr(F, attr), attr
    assert (set(AF.__all__) | set(LL.__all__)) - {"ZhukhovskyFlow"} <= set(PART_C) and set(PART_C) <= set(ch14.__all__)
    assert isinstance(F, ch14.ZhukhovskyFlow) and ch14.ZhukhovskyFlow is AF.ZhukhovskyFlow


def test_contract_V7_reexports_are_the_same_objects():  # V7 — row I32: ch14.<name> reaches the tested originals
    for n in AF.__all__:
        assert getattr(ch14, n) is getattr(AF, n), n
    for n in LL.__all__:
        assert getattr(ch14, n) is getattr(LL, n), n
    assert ch14.joukowski is CONF.joukowski and ch14.joukowski_derivative is CONF.joukowski_derivative
    assert ch14.joukowski_inverse is CONF.joukowski_inverse and ch14.circle_flow_zeta is CONF.circle_flow_zeta
    assert ch14.lift_per_span is ch06.lift_per_span and ch14.cylinder_surface_speed is ch06.cylinder_surface_speed
    assert ch14.cylinder_stagnation_points is ch06.cylinder_stagnation_points
    assert ch14.segment_speed is BS.segment_speed and ch14.segment_induced_velocity is BS.segment_induced_velocity
    assert ch14.point_vortex_velocity is BS.point_vortex_velocity and ch14.point_vortex_evolve is BS.point_vortex_evolve
    assert ch14.point_vortex_invariants is BS.point_vortex_invariants
    assert ch14.blasius_force is POT.blasius_force and ch14.pressure_coefficient is POT.pressure_coefficient
    assert ch14.lift_coefficient is SIM.lift_coefficient and ch14.drag_coefficient is SIM.drag_coefficient
    assert ch14.mach_number is SIM.mach_number and ch14.reynolds_number is SIM.reynolds_number
    assert ch14.AF is AF and ch14.LL is LL


def test_contract_V7_signatures_match_design_part_c():  # V7 — argument order, defaults and keyword-only markers
    def sig(f):
        return str(inspect.signature(f)).replace(" ", "")

    def params(f):
        return [(p.name, p.kind is p.KEYWORD_ONLY) for p in inspect.signature(f).parameters.values()]

    want = {
        ch14.zhukhovsky_circle: [("b", False), ("e", False), ("beta", False), ("a", True)],
        ch14.zhukhovsky_airfoil: [("b", False), ("e", False), ("beta", False), ("n", False), ("a", True)],
        ch14.zhukhovsky_shape: [("b", False), ("e", False), ("beta", False), ("exact", False), ("a", True), ("n", True)],
        ch14.kutta_circulation: [("U", False), ("a", False), ("alpha", False), ("beta", False)],
        ch14.kutta_circulation_root: [("U", False), ("b", False), ("e", False), ("beta", False), ("alpha", False), ("xtol", False)],
        ch14.rear_stagnation_point: [("U", False), ("a", False), ("alpha", False), ("beta", False), ("Gamma_cw", True)],
        ch14.zhukhovsky_flow: [("U", False), ("b", False), ("e", False), ("beta", False), ("alpha", False), ("Gamma_cw", True), ("a", True)],
        ch14.zhukhovsky_surface: [("U", False), ("b", False), ("e", False), ("beta", False), ("alpha", False), ("n", False), ("Gamma_cw", True)],
        ch14.zhukhovsky_force: [("U", False), ("b", False), ("e", False), ("beta", False), ("alpha", False), ("rho", False),
                                ("method", False), ("Gamma_cw", True), ("n", True)],
        ch14.surface_point: [("U", False), ("b", False), ("e", False), ("beta", False), ("alpha", False), ("t", False), ("Gamma_cw", True)],
        ch14.edge_speed: [("U", False), ("b", False), ("e", False), ("beta", False), ("alpha", False), ("d", False), ("Gamma_cw", True)],
        ch14.map_point: [("zeta_re", False), ("zeta_im", False), ("b", False)],
        ch14.downwash_series: [("Gn", False), ("gamma1", False), ("s", False), ("printed_middle", False)],
        ch14.section_lift_coefficient: [("alpha", False), ("w", False), ("U", False), ("beta", False), ("K", True)],
        ch14.solve_lifting_line: [("c", False), ("s", False), ("U", False), ("alpha", False), ("K", True), ("beta", True), ("twist", True),
                                  ("N", True), ("symmetric", True), ("rho", True)],
        ch14.lifting_line_matrix: [("c", False), ("s", False), ("K", False), ("N", False), ("symmetric", True)],
        ch14.elliptic_wing: [("alpha", False), ("AR", False), ("K", True), ("beta", True)],
        ch14.drag_polar: [("CL", False), ("CD0", True), ("AR", True), ("e_span", True)],
        ch14.sail_forces: [("U", False), ("V", False), ("theta", False), ("CL", True), ("CD", True), ("rho", True), ("A", True)],
        ch14.stall_regime: [("Re_c", False), ("Re_laminar", True), ("Re_turbulent", True)],
        ch14.generic_lift_curve: [("alpha", False), ("K", True), ("beta", True), ("alpha_stall", True), ("CL_post", True), ("width", True)],
        ch14.impulsive_start_lumped_vortex: [("U", False), ("c", False), ("alpha", False), ("t_end", False), ("dt", False), ("wake", True),
                                             ("shed_fraction", True), ("eps", True), ("rho", True), ("history", True)],
    }
    for f, w in want.items():
        assert params(f) == w, (f.__name__, params(f))
    # K has no default anywhere (trap T12): calling without it is a TypeError
    for call in (lambda: ch14.solve_lifting_line(1.0, 6.0, 1.0, 0.1), lambda: ch14.elliptic_wing(0.1, 6.0),
                 lambda: ch14.section_lift_coefficient(0.1, 0.0, 1.0), lambda: ch14.lifting_line_matrix(1.0, 6.0, N=4)):
        with pytest.raises(TypeError):
            call()
    pz = inspect.signature(CONF.circle_flow_zeta).parameters
    assert pz["center"].default == 0j and pz["center"].kind is pz["center"].KEYWORD_ONLY and pz["alpha"].default == 0.0
    assert "printed_middle" in sig(ch14.downwash_series) and inspect.signature(ch14.downwash_series).parameters["printed_middle"].default is False
    assert inspect.signature(ch14.zhukhovsky_shape).parameters["n"].default == 20001
    assert inspect.signature(ch14.solve_lifting_line).parameters["N"].default == 20


def test_contract_V7_scalars_in_floats_out():  # V7 — the parity-row convention (design Part C header)
    scalars = [ch14.kutta_circulation(US, 0.27, ALPHA, BETA), ch14.trailing_edge_speed(US, 0.27, B, ALPHA, BETA),
               ch14.zhukhovsky_lift_coefficient(ALPHA, BETA), ch14.zhukhovsky_lift_coefficient(ALPHA, BETA, E, exact=True),
               ch14.aspect_ratio(12.0, 15.0), ch14.mean_chord(15.0, 12.0), ch14.span_angle(1.0, SPAN), ch14.span_position(1.0, SPAN),
               ch14.circulation_series(GN, 1.0), ch14.downwash_series(GN, 1.0, SPAN), ch14.elliptic_downwash(2.0, SPAN),
               ch14.glauert_integral(3, 1.0), ch14.downwash(0.5, lambda y: ch14.elliptic_circulation(y, SPAN, 2.0), SPAN),
               ch14.induced_drag_coefficient(0.5, 9.6), ch14.finite_wing_lift_slope(2 * PI, 6.0), ch14.ground_effect_factor(1.0, 4.0),
               ch14.horseshoe_downwash(0.0, 0.0, 1.0, 2.0), ch14.stall_speed(9e3, 1.2, 15.0, 1.4), ch14.chord_reynolds(30.0, 1.0, 1.5e-5),
               ch14.wagner_function(2.0), ch14.drag_polar(0.5, CD0=0.02, AR=6.0), ch14.drag_share(0.5, CD0=0.02, AR=6.0),
               ch14.velocity_made_good(3.0, 0.7), ch14.sailboat_course_limit(6.0, 8.0), ch14.edge_speed(*SEC, 1e-3),
               ch14.zero_lift_angle(0.02), ch14.camber_from_beta(0.07, 1.0), ch14.section_coefficient(10.0, 1.2, 3.0, 1.0),
               ch14.circulation_from_cl(0.5, 30.0, 1.0), ch14.cl_from_circulation(7.5, 30.0, 1.0), ch14.induced_angle(1.0, 50.0),
               ch14.effective_angle(0.1, 1.0, 50.0), ch14.elliptic_induced_drag(9e3, 1.2, 50.0, 12.0),
               ch14.vortex_pair_energy_per_length(100.0, 9.0, 0.5, 1.2), ch14.speed_min_drag(9e3, 1.2, 15.0, 9.6, 0.025),
               ch14.projectile_range_drag(0.15, 1.2, 0.004, 0.4, 30.0, 0.1), ch14.profile_drag_estimate(1e6, 0.12, regime="laminar", form_factor=1.0)]
    assert all(type(v) is float for v in scalars), [type(v) for v in scalars if type(v) is not float]
    mu, a = ch14.zhukhovsky_circle(B, E, BETA)
    assert type(mu) is complex and type(a) is float
    for d in (ch14.surface_point(*SEC, 1.0), ch14.map_point(0.3, 0.4, 0.25), ch14.tip_vortex_estimate(9e3, 1.2, 50.0, 12.0),
              ch14.best_glide(0.025, 9.6, 0.9), ch14.rolled_up_pair(GN, SPAN), ch14.flat_plate_forces(1.2, 30.0, 1.0, 0.1)):
        assert all(type(v) is float for v in d.values()), d
    assert ch14.downwash_series(GN, np.array([0.5, 1.0]), SPAN).shape == (2,)
    assert ch14.circulation_series(GN, np.linspace(0, PI, 5)).shape == (5,)


def test_tables_V7_inputs_conventions_slips_traps():  # V7 — N-rows of the tables; rule 9 (inputs are ours)
    keys = ("b e beta alpha U_section W s A AR U rho nu CD0 e_span K CL_max CL_max_flap taper alpha_stall CL_post stall_width "
            "CD_sep Re_laminar Re_turbulent U_wind V_over_U CL_sail CD_sail A_sail LD_hull CDA_hull rho_water theta U_fish "
            "V_tail L_tail LD_tail Gamma_v ell").split()
    assert set(keys) <= set(I) and ch14.illustrative_inputs() == I and ch14.illustrative_inputs() is not I       # a fresh dict every call
    # the design's input set (header convention 11) — the implementer switched to it; the notebook's numbers rely on it
    assert (I["b"], I["e"], I["U_section"], I["W"], I["s"], I["A"], I["CD0"], I["e_span"]) == (0.25, 0.08, 30.0, 9000.0, 12.0, 15.0, 0.025, 0.9)
    assert close(I["beta"], np.deg2rad(4.0)) and close(I["alpha"], np.deg2rad(5.0)) and close(I["K"], 2 * PI)
    assert close(I["AR"], ch14.aspect_ratio(I["s"], I["A"])) and close(I["c_mean"], ch14.mean_chord(I["A"], I["s"]))
    assert close(I["V_over_U"], I["V_boat"] / I["U_wind"]) and close(I["LD_tail"], I["L_tail"] / I["D_tail"])
    assert (I["wd_s"], I["wd_c"], I["wd_U"], I["wd_rho"]) == (6.0, 1.0, 1.0, 1.0)
    slips = ch14.book_slips()
    assert list(slips) == [str(k) for k in range(1, 11)]
    assert [slips[k]["kind"] for k in slips] == ["slip", "slip", "loose", "loose", "slip", "loose", "loose", "loose", "loose", "loose"]
    for row in slips.values():
        assert {"id", "kind", "where", "printed", "correct", "corrected", "how_to_tell", "coded_in", "test"} <= set(row)
        assert row["correct"] == row["corrected"]
    assert "printed_middle=True" in slips["1"]["test"] and "arctan" in slips["2"]["test"] and "14.8" in slips["4"]["where"]
    tr = ch14.traps()
    assert [t["id"] for t in tr] == [f"T{k}" for k in range(1, 19)] and all({"what", "where", "code"} <= set(t) for t in tr)
    tab = ch14.conventions_table()
    assert list(tab.columns) == ["symbol", "meanings", "ours", "code"] and len(tab) >= 12
    assert {"Γ", "b", "a", "β", "e", "Λ", "γ", "θ", "K", "A", "w"} <= set(tab["symbol"])


def test_public_V7_no_private_book_number_is_typed_in_this_file():  # V7 — rule 9
    src = Path(__file__).read_text(encoding="utf-8")
    assert "book_values_ch14.json" in src
    if BOOK.exists():
        bv = book()
        for s in bv["_forbidden_public"]:
            assert s not in src.replace("_forbidden_public", ""), s
        for item in bv["_forbidden_public_regex"]:
            assert not re.search(item["regex"], src), item["why"]
        public = src.split("# 99. Private book values", 1)[0]
        ex = bv["exercises"]
        for token in (str(ex["ex_14_13"]["weight_N"]), str(ex["ex_14_14"]["derived_by_analyst"]["Gamma1"]),
                      str(ex["ex_14_16"]["derived_by_analyst"]["M_N_m"]), str(ex["ex_14_13"]["derived_by_analyst"]["Di_N"])):
            assert token not in public, token


# =====================================================================================================================
# 1. The one changed line of core.conformal: circle_flow_zeta(center=0j) must reproduce chapter 6
# =====================================================================================================================
def test_circle_flow_zeta_V7_center_default_reproduces_ch06():  # V7 + V1 — (6.68); the default path is unchanged
    U, a, G = 1.7, 0.6, 2.9
    r, th = np.meshgrid(np.linspace(0.6, 3.0, 13), np.linspace(-3.0, 3.0, 29))
    zeta = r * np.exp(1j * th)
    w0, dw0 = CONF.circle_flow_zeta(U, a, Gamma_cw=G)
    w1, dw1 = ch14.circle_flow_zeta(U, a, Gamma_cw=G, center=0j)
    assert np.array_equal(w0(zeta), w1(zeta)) and np.array_equal(dw0(zeta), dw1(zeta))          # bit for bit
    assert close(w0(zeta), U * (zeta + a ** 2 / zeta) + 1j * G / (2 * PI) * np.log(zeta / a), 1e-13)   # (6.68) as printed
    assert close(dw0(zeta), U * (1 - a ** 2 / zeta ** 2) + 1j * G / (2 * PI * zeta), 1e-13)
    wc, dwc = CONF.circle_flow_zeta(U, a, Gamma_ccw=-G)                                         # the project convention
    assert np.array_equal(wc(zeta), w0(zeta))
    # with a centre and an angle: the same function of (zeta - mu) e^{-i alpha}  (shift, then rotate the stream)
    mu, al = -0.11 + 0.07j, 0.23
    w2, dw2 = CONF.circle_flow_zeta(U, a, Gamma_cw=G, alpha=al, center=mu)
    s = zeta * np.exp(1j * al)                                                                  # points with (zeta' - mu) e^{-i al} = zeta
    lhs = w2(mu + s)
    ref = U * (zeta + a ** 2 / zeta) + 1j * G / (2 * PI) * np.log(s / a)
    assert close(lhs.real, ref.real, 1e-12, 1e-12)
    assert close(dw2(mu + s) * np.exp(1j * al), dw0(zeta), 1e-12)                                # velocity turned by alpha
    # surface speed of the centred circle equals (6.37) = (14.10) of chapter 6
    t = np.linspace(0.0, 2 * PI, 60, endpoint=False)
    ut = (dw0(a * np.exp(1j * t)) * 1j * np.exp(1j * t)).real
    assert close(ut, ch14.cylinder_surface_speed(t, U, a, Gamma_cw=G), 1e-12, 1e-13)
    assert close(ut, -2 * U * np.sin(t) - G / (2 * PI * a), 1e-12, 1e-13)


def test_circle_flow_zeta_V1_offcentre_circle_is_a_streamline_with_the_right_far_field_and_circulation():  # V1 + V4 — N28, D06
    U, G = 3.1, 4.4
    mu, a = ch14.zhukhovsky_circle(0.5, 0.07, 0.12)
    for al in (-0.2, 0.0, 0.31):
        w, dw = ch14.circle_flow_zeta(U, a, Gamma_cw=G, alpha=al, center=mu)
        t = np.linspace(-3.0, 3.0, 301)                                   # avoid the cut of the logarithm at t = pi
        psi = w(mu + a * np.exp(1j * t)).imag
        assert np.ptp(psi) < 1e-12 * U * a                                # the circle is a streamline
        far = dw(mu + 1e7 * np.exp(1j * np.linspace(0, 6, 7)))
        assert close(far, U * np.exp(-1j * al), 1e-6)                     # u - iv -> U e^{-i alpha}
        tt = np.linspace(0.0, 2 * PI, 400, endpoint=False)
        for R, c in ((a, mu), (2.5 * a, mu + 0.3 * a), (9.0, 0.4 + 0.1j)):
            zc = c + R * np.exp(1j * tt)
            loop = np.sum(dw(zc) * 1j * R * np.exp(1j * tt)) * (2 * PI / tt.size)     # circulation + i flux
            assert abs(loop.real + G) < 1e-10 * G and abs(loop.imag) < 1e-10 * G      # counterclockwise circulation = -Gamma_cw
        un = (dw(mu + a * np.exp(1j * tt)) * np.exp(1j * tt)).real                    # radial velocity on the circle
        assert np.max(np.abs(un)) < 1e-12 * U


# =====================================================================================================================
# 2. C01 — the wing in numbers: aspect ratio (14.1), mean chord (14.2), level flight; N01–N08, N83
# =====================================================================================================================
def test_aspect_ratio_V1_planforms_closed_forms():  # V1 — C01 (14.1), N03 (14.2), N02
    s, c0 = 11.0, 1.7
    r = ch14.planform("rectangular", s, c0)
    assert close(r["A"], s * c0) and close(r["AR"], s / c0) and close(r["c_mean"], c0) and close(r["A_trapezoid"], s * c0, 1e-13)
    e = ch14.planform("elliptic", s, c0)
    assert close(e["A"], PI * c0 * s / 4) and close(e["AR"], 4 * s / (PI * c0)) and close(e["c_mean"], PI * c0 / 4)
    t = ch14.planform("tapered", s, c0, taper=0.4)
    assert close(t["A"], 0.5 * c0 * (1 + 0.4) * s) and close(t["A_trapezoid"], t["A"], 1e-13)   # piecewise linear: exact
    assert close(t["chord"](np.array([0.0, s / 2, -s / 2, s / 4])), [c0, 0.4 * c0, 0.4 * c0, 0.7 * c0])
    for p in (r, e, t):
        assert close(p["AR"], ch14.aspect_ratio(p["s"], p["A"])) and close(p["AR"], p["s"] / p["c_mean"])     # both members of (14.1)
        assert close(p["c_mean"], ch14.mean_chord(p["A"], p["s"]))
        area = quad(lambda y: float(p["chord"](y)), -s / 2, s / 2, points=[0.0], epsabs=1e-13, epsrel=1e-13)[0]
        assert close(area, p["A"], 1e-9)                                                                    # A = integral of c(y) dy
        assert close(p["chord"](p["y"]), p["chord"](-p["y"]))                                               # symmetric about mid-span
    assert close(ch14.aspect_ratio(np.array([10.0, 20.0]), 20.0), [5.0, 20.0])
    with pytest.raises(ValueError):
        ch14.planform("delta", s, c0)


def test_aspect_ratio_V2_two_forms_and_dimensions():  # V2 — D01 (sympy) and dimensional homogeneity
    s, A, W, rho, U, CLm = sp.symbols("s A W rho U C_Lmax", positive=True)
    cbar = A / s                                                          # (14.2)
    assert sp.simplify(s ** 2 / A - s / cbar) == 0                        # (14.1): the two members agree
    CL = sp.symbols("C_L", positive=True)
    sol = sp.solve(sp.Eq(CL * sp.Rational(1, 2) * rho * U ** 2 * A, W), CL)[0]      # L = W with (4.108)
    assert sp.simplify(sol - W / (sp.Rational(1, 2) * rho * U ** 2 * A)) == 0
    Umin = sp.solve(sp.Eq(sol.subs(U, sp.Symbol("Um", positive=True)), CLm), sp.Symbol("Um", positive=True))[0]
    assert sp.simplify(Umin - sp.sqrt(2 * W / (rho * A * CLm))) == 0
    num = {s: 12.0, A: 15.0, W: 9000.0, rho: 1.225, U: 30.0, CLm: 1.4}
    assert close(ch14.level_flight(9000.0, 1.225, 30.0, 15.0)["CL"], float(sol.subs(num)))
    assert close(ch14.stall_speed(9000.0, 1.225, 15.0, 1.4), float(Umin.subs(num)))
    # dimensions with pint
    Wq, rq, Uq, Aq, sq = Q_(9000.0, "N"), Q_(1.225, "kg/m**3"), Q_(30.0, "m/s"), Q_(15.0, "m**2"), Q_(12.0, "m")
    assert (sq ** 2 / Aq).check("[]") and (Aq / sq).check("[length]")
    assert (Wq / (0.5 * rq * Uq ** 2 * Aq)).check("[]") and ((2 * Wq / (rq * Aq)) ** 0.5).check("[length]/[time]")
    assert (rq * Uq * Q_(1.0, "m**2/s")).check("[force]/[length]")       # L = rho U Gamma is a force per unit span (P336)
    assert (rq * Q_(1.0, "m**2/s") ** 2).check("[force]")                # rho Gamma^2 is a force: (14.22), (14.25)
    assert (Q_(1.0, "m**2/s") / Q_(1.0, "m")).check("[length]/[time]")   # w = Gamma_1 / 2s is a velocity (14.24)


def test_planform_V3_trapezoid_area_of_the_ellipse_converges_at_order_three_halves():  # V3 — N02 (square-root ends)
    ns, errs = [41, 81, 161, 321, 641], []
    for n in ns:
        p = ch14.planform("elliptic", 9.0, 1.3, n=n)
        errs.append(abs(p["A_trapezoid"] - p["A"]) / p["A"])
    order = observed_order([9.0 / (n - 1) for n in ns], errs)
    assert abs(order - 1.5) < 0.15, (order, errs)
    assert errs[-1] < 1e-4


def test_level_flight_V1_lift_equals_weight_and_stall_speed_inverts():  # V1 — N05, N73
    W, rho, A = I["W"], I["rho"], I["A"]
    U = np.array([I["U_slow"], I["U_cruise"], 80.0])
    lf = ch14.level_flight(W, rho, U, A, CD=0.03)
    assert close(ch14.lift_coefficient(W, rho, U, A), lf["CL"])                            # (4.108) with L = W
    assert close(lf["CL"] * lf["q"] * A, W, 1e-13)                                         # the lift it implies is the weight
    assert close(lf["D"], ch14.drag_coefficient(1.0, rho, U, A) ** -1 * 0.03, 1e-13)       # D = C_D (1/2 rho U^2 A)
    assert close(lf["power"], lf["D"] * U)
    assert close(lf["CL"][0] / lf["CL"][1], (U[1] / U[0]) ** 2)                            # C_L falls as 1/U^2
    none = ch14.level_flight(W, rho, 30.0, A)
    assert none["D"] is None and none["power"] is None and type(none["CL"]) is float
    for CLm in (I["CL_max"], I["CL_max_flap"]):
        Um = ch14.stall_speed(W, rho, A, CLm)
        assert close(ch14.level_flight(W, rho, Um, A)["CL"], CLm, 1e-13)
    assert ch14.stall_speed(W, rho, A, I["CL_max_flap"]) < ch14.stall_speed(W, rho, A, I["CL_max"])   # flaps: slower landing
    assert close(ch14.chord_reynolds(30.0, 1.25, 1.5e-5), ch14.reynolds_number(30.0, 1.25, 1.5e-5))


def test_compressibility_V1_half_mach_squared_against_the_isentropic_density_ratio():  # V1 + V7 — N01
    gam, R, T = 1.4, 287.05, 288.15
    U = np.array([10.0, 50.0, 102.0, 250.0])
    c = ch14.compressibility_check(U, T)
    Ma = U / math.sqrt(gam * 287.0 * T)
    assert close(c["Ma"], ch14.mach_number(U, T=T)) and close(c["Ma"], Ma, 2e-4)            # c = sqrt(gamma R T), R of dry air
    exact = (1 + 0.5 * (gam - 1) * c["Ma"] ** 2) ** (1 / (gam - 1)) - 1                     # rho_0/rho - 1, isentropic
    assert close(c["density_change"], 0.5 * c["Ma"] ** 2)
    err = np.abs(c["density_change"] - exact)
    assert err[2] < 1e-3 and close(c["Ma"][2], 0.3, 2e-3)                                   # the docstring's "1e-3 at Ma = 0.3"
    assert close(err / c["Ma"] ** 4, (2 - gam) / 8, 0.08)                                   # next term of the series: (2 - gamma)/8 Ma^4
    assert list(c["incompressible_ok"]) == [True, True, True, False]
    edge = ch14.compressibility_check(np.array([0.29, 0.31]) * math.sqrt(gam * 287.0 * T), T)                # the limit itself: Ma < 0.3
    assert list(edge["incompressible_ok"]) == [True, False]
    one = ch14.compressibility_check(I["U_cruise"])
    assert one["incompressible_ok"] is True and type(one["Ma"]) is float


def _exact_symmetric(e, b=1.0):
    """Chord and maximum thickness of the symmetric Zhukhovsky section by our own route (closed form + optimiser)."""
    a = b * (1 + e)
    chord = 2 * b + b * (1 + 2 * e) + b / (1 + 2 * e)                     # z(b) - z(-b(1 + 2e))
    y = lambda t: -((-b * e + a * np.exp(1j * t)) + b ** 2 / (-b * e + a * np.exp(1j * t))).imag
    r = minimize_scalar(y, bounds=(1.0, 3.0), method="bounded", options={"xatol": 1e-13})
    zt = (-b * e + a * np.exp(1j * r.x)) + b ** 2 / (-b * e + a * np.exp(1j * r.x))
    return chord, -2 * r.fun, zt.real, r.x


def test_section_geometry_V1_plate_arc_and_symmetric_sections():  # V1 — N07, N24, N26
    g = ch14.zhukhovsky_shape(1.0, 0.0, 0.0)
    assert close(g["chord"], 4.0, 1e-10) and abs(g["camber"]) < 1e-12 and abs(g["t_max"]) < 1e-12        # the flat plate
    for beta in (0.05, 0.2, 0.35):
        g = ch14.zhukhovsky_shape(1.3, 0.0, beta)                                                        # the circular arc
        assert close(g["chord"], 4 * 1.3, 1e-8) and close(g["camber"], 2 * 1.3 * math.tan(beta), 1e-7)
        assert abs(g["t_max"]) < 1e-7 and abs(g["x_camber_over_c"] - 0.5) < 1e-3      # position of a flat maximum
    for e in (0.03, 0.1, 0.2):
        chord, tmax, xt, _ = _exact_symmetric(e, 0.7)
        g = ch14.zhukhovsky_shape(0.7, e, 0.0)
        assert close(g["chord"], chord, 1e-9) and close(g["t_max"], tmax, 1e-7) and abs(g["camber"]) < 1e-9
        assert abs(g["x_tmax"] - xt) < 2e-4 * chord and g["exact"] is True
    # an ellipse given clockwise: chord 2A, thickness 2B at mid-chord, no camber — exact when the tip is a vertex
    t = np.linspace(0.0, -2 * PI, 3000, endpoint=False)
    ge = ch14.section_geometry(1.7 * np.cos(t), 0.2 * np.sin(t))
    assert close(ge["chord"], 3.4, 1e-12) and close(ge["t_max"], 0.4, 1e-6) and abs(ge["camber"]) < 1e-12 and abs(ge["x_tmax"] - 0.5) < 1e-6
    # ... and only first order in the spacing when a ROUNDED trailing edge falls between two vertices (the trailing
    # edge is the vertex of largest x, not refined; every section of this chapter has its sharp edge at index 0) — the docstring of section_geometry now says so
    for n in (1000, 3000, 6000):
        t = np.linspace(0.3, 0.3 - 2 * PI, n, endpoint=False)
        gr = ch14.section_geometry(1.7 * np.cos(t), 0.2 * np.sin(t))
        assert close(gr["chord"], 3.4, 1e-5) and close(gr["t_max"], 0.4, 1e-5) and abs(gr["camber"]) < 3 * 0.2 * PI / n
    # rotated section: the chord angle is reported and the shape numbers do not change (V7)
    A = ch14.zhukhovsky_airfoil(1.0, 0.09, 0.1, 4000)
    g0 = ch14.section_geometry(A["x"], A["y"], te_index=0)
    zr = A["z"] * np.exp(-0.3j)
    g1 = ch14.section_geometry(zr.real, zr.imag, te_index=0)
    assert close([g1["chord"], g1["camber"], g1["t_max"]], [g0["chord"], g0["camber"], g0["t_max"]], 1e-9)
    assert abs((g1["chord_angle"] - g0["chord_angle"]) + 0.3) < 1e-12
    assert set(g0) >= {"chord", "camber", "x_camber", "t_max", "x_tmax", "camber_line", "camber_ratio", "thickness_ratio", "le", "te"}


def test_section_geometry_V3_second_order_in_the_number_of_points():  # V3 — N07 (docstring: second order)
    ns = [250, 500, 1000, 2000, 4000]
    chord, tmax, _, _ = _exact_symmetric(0.1)
    et, ec, ecam = [], [], []
    for n in ns:
        A = ch14.zhukhovsky_airfoil(1.0, 0.1, 0.0, n)
        g = ch14.section_geometry(A["x"], A["y"], te_index=0)
        et.append(abs(g["t_max"] - tmax))
        ec.append(abs(g["chord"] - chord))
        Aa = ch14.zhukhovsky_airfoil(1.0, 0.0, 0.2, n)
        ecam.append(abs(ch14.section_geometry(Aa["x"], Aa["y"], te_index=0)["camber"] - 2 * math.tan(0.2)))
    h = [1.0 / n for n in ns]
    assert observed_order(h, et) > 1.75 and observed_order(h, ecam) > 1.75, (observed_order(h, et), observed_order(h, ecam), et, ecam)
    assert et[-1] < 1e-6 and ecam[-1] < 1e-6 and max(ec) < 1e-5


def test_resolve_lift_drag_V1_rotation_and_section_coefficients():  # V1 — N08, R01, N47
    Fx, Fy = np.array([3.0, -1.2, 0.0]), np.array([4.0, 2.5, 7.0])
    for al in (-0.3, 0.0, 0.2, 1.1):
        L, D = ch14.resolve_lift_drag(Fx, Fy, al)
        s_hat, n_hat = np.array([math.cos(al), math.sin(al)]), np.array([-math.sin(al), math.cos(al)])
        assert close(D, Fx * s_hat[0] + Fy * s_hat[1], 1e-13, 1e-15) and close(L, Fx * n_hat[0] + Fy * n_hat[1], 1e-13, 1e-15)
        assert close(L ** 2 + D ** 2, Fx ** 2 + Fy ** 2, 1e-13)
    L, D = ch14.resolve_lift_drag(3.0, 4.0, 0.0)
    assert (L, D) == (4.0, 3.0)                                            # alpha = 0: lift is F_y
    L, D = ch14.resolve_lift_drag(-math.sin(0.4), math.cos(0.4), 0.4)      # a force normal to the stream is all lift
    assert close(L, 1.0) and abs(D) < 1e-15
    rho, U, c, G = 1.21, 33.0, 0.9, 7.7
    Lp = ch14.lift_per_span(rho, U, Gamma_cw=G)                            # (6.62)
    assert close(Lp, rho * U * G) and Lp > 0
    CL = ch14.section_coefficient(Lp, rho, U, c)
    assert close(CL, ch14.cl_from_circulation(G, U, c)) and close(ch14.circulation_from_cl(CL, U, c), G) and close(CL, 2 * G / (U * c))
    assert close(ch14.section_coefficient(Lp, rho, U, c), ch14.lift_coefficient(Lp, rho, U, c * 1.0))   # chord x unit span as the area


def test_projectile_V7_vacuum_limit_and_V3_shallow_trajectories_approach_the_formula():  # V7 + V3 — N83
    m, rho, A, V0, g = 0.151, 1.19, 0.0043, 31.0, 9.80665
    for th in (0.05, 0.2):
        vac = 2 * V0 ** 2 * math.cos(th) * math.sin(th) / g
        assert close(ch14.projectile_range_drag(m, rho, A, 1e-9, V0, th), vac, 1e-6)                    # C_D -> 0
        assert close(ch14.projectile_range_drag(m, rho, A, 0.0, V0, th), vac, 1e-13)
        k = rho * A * 0.45
        assert close(ch14.projectile_range_drag(m, rho, A, 0.45, V0, th),
                     2 * m / k * math.log(1 + k * V0 ** 2 * math.cos(th) * math.sin(th) / (m * g)), 1e-13)
    gaps = []
    for th in (0.4, 0.2, 0.1, 0.05):
        r = ch14.projectile_range_drag(m, rho, A, 0.45, V0, th, numeric=True)
        assert r["x_numeric"] < r["x_formula"] < r["x_vacuum"]                                           # drag shortens the range
        gaps.append(abs(r["x_numeric"] / r["x_formula"] - 1))
    assert gaps == sorted(gaps, reverse=True) and gaps[-1] < 0.015                                       # -> the formula as theta -> 0
    assert 0.8 < observed_order([0.4, 0.2, 0.1, 0.05], gaps) < 1.2
    r0 = ch14.projectile_range_drag(m, rho, A, 0.0, V0, 0.3, numeric=True)                               # no drag: the parabola exactly
    assert close(r0["x_numeric"], r0["x_vacuum"], 1e-8) and close(r0["t_flight"], 2 * V0 * math.sin(0.3) / g, 1e-8)


# =====================================================================================================================
# 3. C02 — the Kutta condition; N10, N11, N30, R04, R09
# =====================================================================================================================
SECTIONS = [(0.0, 0.0), (0.05, 0.0), (0.0, 0.12), (0.08, 0.07), (0.15, 0.2), (0.1, -0.05)]     # (e, beta) pairs of ours
ANGLES = (-0.15, 0.03, 0.2)


def test_kutta_condition_V1_root_finder_equals_14_11_and_the_edge_point_stagnates():  # V1 — C02, N30 (14.11)
    U, b = 2.7, 0.6
    for e, beta in SECTIONS:
        mu, a = ch14.zhukhovsky_circle(b, e, beta)
        for al in ANGLES:
            GK = ch14.kutta_circulation(U, a, al, beta)
            scale = 4 * PI * U * a
            assert abs(ch14.kutta_circulation_root(U, b, e, beta, al) - GK) < 1e-12 * scale       # independent route: brentq on u_theta(B)
            _, dw = ch14.circle_flow_zeta(U, a, Gamma_cw=GK, alpha=al, center=mu)
            assert abs(complex(dw(b))) < 1e-13 * U                                                # dw/dzeta = 0 at zeta = b
            F = ch14.zhukhovsky_flow(U, b, e, beta, al)
            assert F.is_kutta and close(F.Gamma_cw, GK, 1e-14, 1e-300) and close(F.Gamma_kutta, GK, 1e-14, 1e-300)
            assert np.sign(GK) == np.sign(round(al + beta, 12))                                   # trap T1: clockwise, positive for alpha + beta > 0
            # planted variants (analysis T11): sin(alpha - beta), radius b instead of a, counterclockwise sign
            for bad in (4 * PI * U * a * math.sin(al - beta), 4 * PI * U * b * math.sin(al + beta), -GK):
                if abs(bad - GK) > 1e-3 * scale:
                    _, dwb = ch14.circle_flow_zeta(U, a, Gamma_cw=bad, alpha=al, center=mu)
                    assert abs(complex(dwb(b))) > 1e-4 * U, (e, beta, al)
                    assert not ch14.zhukhovsky_flow(U, b, e, beta, al, Gamma_cw=bad).is_kutta
    assert ch14.kutta_circulation(U, 0.7, -0.12, 0.12) == 0.0                                     # V7: no lift at alpha = -beta
    al = np.linspace(-0.2, 0.3, 11)
    assert close(ch14.kutta_circulation(U, 0.7, al, 0.1), 4 * PI * U * 0.7 * np.sin(al + 0.1), 1e-14, 1e-300)


def test_kutta_condition_V5_published_form_of_the_circulation():  # V5 — C02/N30 against the cited form (b = 1)
    form = bench()["joukowsky_kutta_circulation"]
    assert "asin(mu_y / R)" in form["form"] and form["verified"] == "2026-10-12"
    V = 3.3
    for e, beta in SECTIONS:
        mu, a = ch14.zhukhovsky_circle(1.0, e, beta)                       # the source's map is z = zeta + 1/zeta: b = 1
        R = math.hypot(1 - mu.real, mu.imag)                               # the source's radius: the circle passes through zeta = 1
        assert close(R, a, 1e-13)
        for al in ANGLES:
            cited = 4 * PI * V * R * math.sin(al + math.asin(mu.imag / R))
            assert close(ch14.kutta_circulation(V, a, al, beta), cited, 1e-12, 1e-13)
            assert close(ch14.kutta_circulation_root(V, 1.0, e, beta, al), cited, 1e-11, 1e-12)


def test_kutta_circulation_V2_derivation():  # V2 — D09 (two-star): (14.10) at theta = -(alpha + beta) -> (14.11) -> (14.12)
    U, a, b, G = sp.symbols("U a b Gamma", positive=True)
    al, be, th = sp.symbols("alpha beta theta", real=True)
    zeta = sp.symbols("zeta")
    mu = b - a * sp.exp(-sp.I * be)                                        # the circle through zeta = b (C04)
    w = U * ((zeta - mu) * sp.exp(-sp.I * al) + a ** 2 * sp.exp(sp.I * al) / (zeta - mu)) + sp.I * G / (2 * sp.pi) * sp.log((zeta - mu) / a)
    dw = sp.diff(w, zeta)
    # the tangential speed on the circle, theta from the stream-aligned diameter through Q: (14.10)
    on = dw.subs(zeta, mu + a * sp.exp(sp.I * (th + al)))
    u_th = sp.re(sp.expand(sp.simplify(on * sp.I * sp.exp(sp.I * (th + al))), complex=True))
    assert sp.simplify(u_th - (-2 * U * sp.sin(th) - G / (2 * sp.pi * a))) == 0
    u_r = sp.im(sp.expand(sp.simplify(on * sp.I * sp.exp(sp.I * (th + al))), complex=True))       # -(radial speed): zero on the circle
    assert sp.simplify(u_r) == 0
    # B = zeta = b is at theta = -(alpha + beta)
    assert sp.simplify((mu + a * sp.exp(sp.I * (th + al))).subs(th, -(al + be)) - b) == 0
    GK = sp.solve(sp.Eq((-2 * U * sp.sin(th) - G / (2 * sp.pi * a)).subs(th, -(al + be)), 0), G)[0]
    assert sp.simplify(GK - 4 * sp.pi * U * a * sp.sin(al + be)) == 0                             # (14.11)
    # the same from the Kutta condition written as dw/dzeta = 0 at zeta = b
    G2 = sp.solve(sp.Eq(dw.subs(zeta, b), 0), G)[0]
    assert sp.simplify(sp.expand((G2 - 4 * sp.pi * U * a * sp.sin(al + be)).rewrite(sp.exp))) == 0
    # (14.12): L = rho U Gamma (6.62), c = 4a, sin x = x
    rho, c = sp.symbols("rho c", positive=True)
    CL = (rho * U * GK / (sp.Rational(1, 2) * rho * U ** 2 * c)).subs(c, 4 * a)
    assert sp.simplify(CL - 2 * sp.pi * sp.sin(al + be)) == 0
    x = sp.symbols("x")
    assert sp.series(2 * sp.pi * sp.sin(x), x, 0, 3).removeO() == 2 * sp.pi * x                   # (14.12), with an error of order x^3
    # the finite speed at the cusp (N29): l'Hopital on (14.9), dw/dz -> (d2w/dzeta2) zeta / 2 at zeta = b
    lim = (sp.diff(w, zeta, 2) * zeta / 2).subs(zeta, b).subs(G, 4 * sp.pi * U * a * sp.sin(al + be))
    assert sp.simplify(sp.expand((lim - U * b / a * sp.cos(al + be) * sp.exp(2 * sp.I * be)).rewrite(sp.exp))) == 0
    num = {U: 30.0, a: 0.27, b: 0.25, al: 0.09, be: 0.07}
    assert close(ch14.trailing_edge_speed(30.0, 0.27, 0.25, 0.09, 0.07), abs(complex(lim.subs(num))), 1e-12)


def test_stagnation_points_V1_both_points_for_any_circulation():  # V1 — C02, N14, R09 (14.10) = 0
    U, b, e, beta, al = 2.2, 0.5, 0.09, 0.1, 0.17
    mu, a = ch14.zhukhovsky_circle(b, e, beta)
    GK = ch14.kutta_circulation(U, a, al, beta)
    Gmax = 4 * PI * U * a
    for G in (0.0, 0.5 * GK, GK, 1.7 * GK, -0.6 * GK, 0.999 * Gmax, 1.4 * Gmax, -1.2 * Gmax):
        st = ch14.stagnation_points(U, b, e, beta, al, Gamma_cw=G)
        _, dw = ch14.circle_flow_zeta(U, a, Gamma_cw=G, alpha=al, center=mu)
        assert abs(complex(dw(st["zeta_rear"]))) < 1e-11 * U and abs(complex(dw(st["zeta_front"]))) < 1e-11 * U
        assert close(st["z_rear"], ch14.joukowski(st["zeta_rear"], b), 1e-13) and close(st["z_front"], ch14.joukowski(st["zeta_front"], b), 1e-13)
        assert st["Gamma_cw"] == G and close(st["Gamma_kutta"], GK, 1e-14)
        if abs(G) <= Gmax:
            assert st["on_surface"] and close(abs(st["zeta_rear"] - mu), a, 1e-12) and close(abs(st["zeta_front"] - mu), a, 1e-12)
            k = G / Gmax
            assert close(st["theta_rear"], -math.asin(k), 1e-13, 1e-15) and close(st["theta_front"], PI + math.asin(k), 1e-13)   # (6.38)
            assert close(st["zeta_rear"], mu + a * np.exp(1j * (st["theta_rear"] + al)), 1e-12)        # theta is from the stream-aligned diameter
            assert close(ch14.cylinder_surface_speed(st["theta_rear"], U, a, Gamma_cw=G), 0.0, 0, 1e-12 * U)
            pts6 = ch14.cylinder_stagnation_points(U, a, Gamma_cw=G)                                   # chapter 6, centred, alpha = 0
            mine = np.sort_complex(np.array([(st["zeta_rear"] - mu) * np.exp(-1j * al), (st["zeta_front"] - mu) * np.exp(-1j * al)]))
            assert close(np.sort_complex(np.asarray(pts6)), mine, 1e-9, 1e-9)
            rs = ch14.rear_stagnation_point(U, a, al, beta, Gamma_cw=G)
            assert close(rs["theta"], st["theta_rear"], 1e-14, 1e-16) and close(rs["phi"], rs["theta"] + al, 1e-14, 1e-16)
            assert close(rs["offset"], rs["theta"] + al + beta, 1e-13, 1e-15) and close(rs["arc_offset"], a * rs["offset"], 1e-13, 1e-15)
            assert close(mu + a * np.exp(1j * (rs["offset"] - beta)), st["zeta_rear"], 1e-12)          # offset is the angle ahead of B
        else:
            assert not st["on_surface"] and abs(st["zeta_rear"] - mu) > a > abs(st["zeta_front"] - mu) and math.isnan(st["theta_rear"])
            assert not ch14.rear_stagnation_point(U, a, al, beta, Gamma_cw=G)["on_surface"]
    k = ch14.stagnation_points(U, b, e, beta, al)                         # default: the Kutta value
    assert k["rear_at_edge"] and abs(k["zeta_rear"] - b) < 1e-13 and abs(k["z_rear"] - 2 * b) < 1e-13 and k["Gamma_cw"] == k["Gamma_kutta"]
    r0 = ch14.rear_stagnation_point(U, a, al, beta, Gamma_cw=0.0)         # Fig. 14.9a: no circulation -> on the suction side
    assert close(r0["offset"], al + beta) and r0["offset"] > 0 and not r0["at_edge"]
    rk = ch14.rear_stagnation_point(U, a, al, beta, Gamma_cw=GK)
    assert rk["at_edge"] and abs(rk["offset"]) < 1e-14
    offs = [ch14.rear_stagnation_point(U, a, al, beta, Gamma_cw=f * GK)["offset"] for f in np.linspace(0, 1.5, 16)]
    assert np.all(np.diff(offs) < 0) and offs[-1] < 0                     # more circulation moves it back, past the edge to the lower side


# =====================================================================================================================
# 4. C03 — the starting vortex (Kelvin); N14–N16
# =====================================================================================================================
def test_starting_vortex_V4_kelvin_total_circulation_is_zero_at_every_step():  # V4 — C03, D02: Gamma_foil + Gamma_starting = 0
    U, c, al = 3.0, 0.8, 0.14
    Gs = PI * U * c * math.sin(al)
    for kw in (dict(wake="rigid"), dict(wake="free"), dict(wake="rigid", shed_fraction=0.4), dict(wake="free", eps=0.02)):
        r = ch14.impulsive_start_lumped_vortex(U, c, al, 6 * c / U, 0.05 * c / U, **kw)
        assert r["t"].size == 120 and np.max(np.abs(r["total"])) < 1e-12 * Gs
        assert close(r["Gamma_bound_cw"] + r["Gamma_wake_cw"], r["total"], 0, 1e-15 * Gs)
        assert close(np.sum(r["wake_Gamma_cw"]), r["Gamma_wake_cw"][-1], 1e-12) and close(r["Gamma_steady"], Gs, 1e-13)
        assert np.all(r["Gamma_bound_cw"] > 0) and np.all(np.diff(r["Gamma_bound_cw"]) > 0)       # clockwise bound circulation grows
        assert np.all(r["wake_Gamma_cw"] < 0)                                                     # every shed vortex is counterclockwise
        assert close(r["s_chords"], U * r["t"] / c, 1e-12) and close(r["lift_circulatory"], 1.0 * U * r["Gamma_bound_cw"], 1e-13)
        assert close(r["ratio"], r["Gamma_bound_cw"] / Gs, 1e-13) and r["model"].startswith("lumped vortex")
    # a new vortex at every change of speed or incidence; the total stays zero through all of them
    Uf = lambda t: 3.0 if t < 1.0 else (4.5 if t < 2.0 else 2.0)          # speed up, then slow down
    r = ch14.impulsive_start_lumped_vortex(Uf, c, al, 3.0, 0.01, history=True)
    assert np.max(np.abs(r["total"])) < 1e-12 * Gs
    k1, k2 = np.searchsorted(r["t"], 1.0 - 1e-9), np.searchsorted(r["t"], 2.0 - 1e-9)
    assert r["wake_Gamma_cw"][k1] < 3 * r["wake_Gamma_cw"][k1 - 1] < 0    # a strong new starting vortex at the speed-up
    assert r["wake_Gamma_cw"][k2] > 0                                     # a stopping vortex (clockwise) at the slow-down
    assert r["wake_x_hist"].shape == (300, 300) and np.isnan(r["wake_x_hist"][0, 1]) and np.isfinite(r["wake_x_hist"][-1]).all()
    assert np.all(np.diff(r["wake_x"]) < 0) and r["wake_x"].min() > c     # older vortices are further downstream, all behind the edge
    ra = ch14.impulsive_start_lumped_vortex(U, c, lambda t: 0.1 if t < 0.5 else 0.2, 1.0, 0.01)   # a step in incidence
    assert np.max(np.abs(ra["total"])) < 1e-12 * Gs and close(ra["Gamma_steady"][-1], PI * U * c * math.sin(0.2))
    with pytest.raises(ValueError):
        ch14.impulsive_start_lumped_vortex(U, c, al, 1.0, 0.1, wake="frozen")


def test_starting_vortex_V1_steady_limit_is_the_kutta_value_and_first_step_is_one_third():  # V1 + V7 — C03 limits
    U, c, al = 1.3, 1.1, 0.09
    r = ch14.impulsive_start_lumped_vortex(U, c, al, 400 * c / U, 0.5 * c / U)
    gap, s = 1 - r["ratio"], r["s_chords"]
    assert np.all(np.diff(r["ratio"]) > 0) and 0 < gap[-1] < 1.5e-3       # -> pi U c sin(alpha) = (14.11) with a = b = c/4
    assert close(r["Gamma_steady"][-1], ch14.kutta_circulation(U, c / 4, al, 0.0), 1e-13)
    big = s >= 50
    assert np.all(np.abs(gap[big] * s[big] - 0.5) < 0.06)                 # the deficit decays algebraically, as 1/(2 chords), not exponentially
    first = [float(ch14.impulsive_start_lumped_vortex(U, c, al, 2 * dt, dt)["ratio"][0]) for dt in (1e-2, 1e-3, 1e-4)]
    assert abs(first[-1] - 1 / 3) < 2e-4 and abs(first[0] - 1 / 3) > abs(first[1] - 1 / 3) > abs(first[2] - 1 / 3)   # the model's Gamma_bound starts at 1/3 (exact circulation: 0; Wagner's 1/2 is a lift ratio)
    # the steady algebra by hand: one bound vortex at c/4, no flow through the plate at 3c/4, wake infinitely far
    Gb = U * math.sin(al) * 2 * PI * (0.75 * c - 0.25 * c)
    assert close(Gb, PI * U * c * math.sin(al)) and close(ch14.flat_plate_forces(1.0, U, c, al)["Gamma_cw"], Gb)


def test_starting_vortex_V3_first_order_in_the_time_step():  # V3 — C03 (docstring: explicit, first order in dt)
    U, c, al = 1.0, 1.0, 0.1
    dts = [0.2, 0.1, 0.05, 0.025, 0.0125]
    for t_obs in (1.0, 4.0):
        val = []
        for dt in dts:
            r = ch14.impulsive_start_lumped_vortex(U, c, al, t_obs, dt)
            assert abs(r["t"][-1] - t_obs) < 1e-12
            val.append(float(r["ratio"][-1]))
        d = np.abs(np.diff(val))                                          # successive differences fall like dt^p
        orders = pairwise_orders(dts[:-1], d)
        assert all(abs(p - 1.0) < 0.15 for p in orders), (t_obs, orders, val)
    rf = [float(ch14.impulsive_start_lumped_vortex(U, c, al, 2.0, dt, wake="free")["ratio"][-1]) for dt in (0.1, 0.05, 0.025, 0.0125)]
    of = pairwise_orders([0.1, 0.05, 0.025], np.abs(np.diff(rf)))
    assert all(abs(p - 1.0) < 0.25 for p in of), (of, rf)


def test_starting_vortex_V5_circulatory_lift_against_the_cited_fits_of_wagners_function():  # V5 + V4 — C03 vs Jones' fit (first-hand, two sources)
    """Wagner's function phi(s) is a LIFT ratio: circulatory lift / steady lift after a step in incidence.

    The quantity compared with it is therefore the model's circulatory lift over the steady lift,
    ``lift_circulatory / (rho U Gamma_steady)`` — not a circulation ratio.  The bound circulation is a different
    quantity: by Kelvin's theorem it equals minus what has been shed, so the exact one starts from 0 (review M1), the
    model's from 1/3, and neither from Wagner's phi(0) = 1/2.
    """
    j, g = bench()["wagner_jones_fit"], bench()["wagner_garrick_fit"]
    s = np.linspace(0.0, 120.0, 2401)
    phi = ch14.wagner_function(s)
    assert close(phi, 1 - j["C1"] * np.exp(-j["eps1"] * s) - j["C2"] * np.exp(-j["eps2"] * s), 1e-15)      # the cited coefficients
    assert abs(ch14.wagner_function(0.0) - j["phi_at_zero"]) < 1e-15 and j["phi_at_zero"] == 0.5 and abs(ch14.wagner_function(1e4) - 1) < 1e-12
    assert np.all(np.diff(phi) > 0) and close(j["C1"] + j["C2"], 0.5, 1e-15)
    garrick = 1 - 2 / (4 + s)                                             # an independent published approximation
    assert np.max(np.abs(phi - garrick)) < j["max_abs_error_stated"] + g["max_abs_error_stated"]     # each within its stated error of the exact one
    # the lumped-vortex model's CIRCULATORY LIFT / STEADY LIFT against the fit (s in semi-chords = 2 x chords travelled);
    # rho is not 1, so a lift built without the density, or from the steady circulation, cannot pass
    U, c, al, rho = I["U_section"], 1.0, I["alpha"], I["rho"]
    r = ch14.impulsive_start_lumped_vortex(U, c, al, 20 * c / U, 0.05 * c / U, rho=rho)
    lift_steady = rho * U * r["Gamma_steady"]                             # Kutta–Zhukhovsky lift of the steady flow [N/m]
    assert close(lift_steady, rho * U * PI * U * c * math.sin(al), 1e-13) and rho != 1.0
    lift_ratio = r["lift_circulatory"] / lift_steady                      # the quantity Wagner's function describes
    assert close(r["lift_circulatory"], rho * U * r["Gamma_bound_cw"], 1e-14)                                # rho U Gamma_bound, by hand
    assert [round(float(lift_ratio[k]), 4) for k in (19, 59, 99, 199)] == [0.6463, 0.8156, 0.8792, 0.9388]   # design C.6 item 4
    # 90 % of the steady lift: 6.13 chords at dt = 0.05 c/U with the rigid wake (linear interpolation between steps;
    # 6.22, 6.13, 6.09, 6.06 at dt = 0.1, 0.05, 0.025, 0.0125 c/U — first order in dt, toward about 6.04)
    assert abs(float(np.interp(0.9, lift_ratio, r["s_chords"])) - 6.13) < 0.01
    w = (r["s_chords"] >= 2.0)
    dev = np.abs(lift_ratio - ch14.wagner_function(2 * r["s_chords"]))
    assert np.max(dev[w]) < 0.012, np.max(dev[w])                         # measured 0.0060 over 2-20 chords at dt = 0.05 c/U (was 0.03)
    assert 0.15 < np.max(dev[~w]) < 0.17, np.max(dev[~w])                 # measured 0.166: before two chords the model is NOT Wagner's lift
    # --- what a circulation is not: Wagner's 1/2 at the start.  The model's bound circulation starts well below it ...
    assert 0.5 - float(r["ratio"][0]) > 0.15 and close(r["ratio"][0], r["Gamma_bound_cw"][0] / (PI * U * c * math.sin(al)), 1e-13)   # 0.344 at this dt
    first = [float(ch14.impulsive_start_lumped_vortex(U, c, al, 2 * q * c / U, q * c / U)["ratio"][0]) for q in (1e-2, 1e-3, 1e-4)]
    assert abs(first[-1] - 1 / 3) < 1e-3 and abs(first[0] - 1 / 3) > abs(first[1] - 1 / 3) > abs(first[2] - 1 / 3)   # 0.3938, 0.3399, 0.33400 -> 1/3
    assert all(0.5 - f > 0.10 for f in first) and abs(float(ch14.wagner_function(0.0)) - first[-1] - 1 / 6) < 1e-3   # 1/2 - 1/3: a different quantity
    # ... because Kelvin ties it to what has been shed: Gamma_bound + sum(Gamma_wake) = 0 to round-off at EVERY step,
    # summed here from the individual wake vortices of a run stopped after each of the first steps (not the stored total)
    assert np.max(np.abs(r["total"])) < 1e-13 * r["Gamma_steady"][-1]                                        # measured 2.2e-16
    assert np.max(np.abs(r["Gamma_bound_cw"] + r["Gamma_wake_cw"])) < 1e-13 * r["Gamma_steady"][-1]
    for m in (1, 2, 7, 40):
        q = ch14.impulsive_start_lumped_vortex(U, c, al, m * 0.05 * c / U, 0.05 * c / U, rho=rho)
        assert len(q["wake_Gamma_cw"]) == m and abs(q["Gamma_bound_cw"][-1] + q["wake_Gamma_cw"].sum()) < 1e-13 * q["Gamma_steady"][-1]
        assert close(q["Gamma_bound_cw"][-1], r["Gamma_bound_cw"][m - 1], 1e-12)                             # the same march
    assert r["wake_Gamma_cw"][0] < 0 < r["Gamma_bound_cw"][0] and close(r["wake_Gamma_cw"][0], -r["Gamma_bound_cw"][0], 1e-13)   # the starting vortex is the opposite of the first bound value


# =====================================================================================================================
# 5. C04 — the Zhukhovsky airfoil (14.3)–(14.8); N19–N27, R05–R07
# =====================================================================================================================
def test_zhukhovsky_circle_V1_reductions_and_round_trip():  # V1 + V7 — N25 (14.7), trap T7
    b = 0.8
    for e, beta in SECTIONS:
        mu, a = ch14.zhukhovsky_circle(b, e, beta)
        assert close(abs(b - mu), a, 1e-14) and close(a, b * (1 + e) / math.cos(beta), 1e-14)     # through zeta = b
        assert abs(-b - mu) <= a * (1 + 1e-14)                                                    # zeta = -b is not outside the circle
        e2, be2, a2 = ch14.zhukhovsky_circle_from_centre(mu, b)
        assert close([e2, be2, a2], [e, beta, a], 1e-12, 1e-14)
        mu3, a3 = ch14.zhukhovsky_circle(b, 99.0, beta, a=a)                                      # (a, b, beta) given: e is ignored
        assert close(mu3, mu, 1e-14) and a3 == a
    mu, a = ch14.zhukhovsky_circle(b, 0.07, 0.0)
    assert close(a, b * 1.07) and close(mu, -0.07 * b) and mu.imag == 0                           # (14.7): a = b(1 + e), centre -be
    mu, a = ch14.zhukhovsky_circle(b, 0.0, 0.3)
    assert close(a, b / math.cos(0.3)) and close(mu, 1j * b * math.tan(0.3), 1e-13, 1e-15)        # the arc: centre i b tan(beta)
    # planted variant (slip #4 row): a = b(1 + e) cos(beta) puts zeta = -b outside the circle: no round nose
    a_bad = b * 1.05 * math.cos(0.3)
    assert abs(-b - (b - a_bad * np.exp(-0.3j))) > a_bad


def test_zhukhovsky_airfoil_V1_slit_arc_symmetric_section_and_14_4():  # V1 — C04, N19–N24, R05, R06
    b, n = 0.6, 720
    P = ch14.zhukhovsky_airfoil(b, 0.0, 0.0, n)
    assert np.max(np.abs(P["y"])) < 1e-15 and close(P["x"], 2 * b * np.cos(P["t"]), 1e-13, 1e-15)          # the slit, (6.66)
    assert close(P["x"].max(), 2 * b) and close(P["x"].min(), -2 * b)
    for e, beta in SECTIONS:
        S = ch14.zhukhovsky_airfoil(b, e, beta, n)
        assert close(S["z"], ch14.joukowski(S["zeta"], b), 1e-13, 1e-15)                                   # (14.4) = (14.3)
        assert close(S["zeta"], S["R"] * np.exp(1j * S["theta"]), 1e-14) and close(S["zeta"], S["mu"] + S["a"] * np.exp(1j * S["t"]), 1e-14)
        assert close(S["t"][0], -beta, 1e-15, 1e-16) and abs(S["z"][0] - 2 * b) < 1e-14                    # index 0: the trailing edge
        assert close(S["e"], e, 1e-12, 1e-14) and S["beta"] == beta
        area = 0.5 * np.sum(S["x"] * np.roll(S["y"], -1) - np.roll(S["x"], -1) * S["y"])
        assert area > -1e-15                                                                               # counterclockwise
        k = ch14.map_point(S["zeta"][37].real, S["zeta"][37].imag, b)
        assert close([k["x"], k["y"], k["R"], k["theta"]], [S["x"][37], S["y"][37], S["R"][37], S["theta"][37]], 1e-13, 1e-15)
    for beta in (0.05, 0.2, 0.4):                                                                          # the circular arc (N24)
        S, arc = ch14.zhukhovsky_airfoil(b, 0.0, beta, n), ch14.circular_arc(b, beta)
        assert np.max(np.abs(np.abs(S["z"] - arc["centre"]) - arc["radius"])) < 1e-12 * b                  # every point on the circle
        assert close(arc["centre"], -2j * b / math.tan(2 * beta)) and close(arc["radius"], 2 * b / math.sin(2 * beta))
        assert close(S["y"].max(), arc["height"], 1e-5) and close(arc["height"], 2 * b * math.tan(beta)) and S["y"].min() > -1e-14
        assert arc["chord"] == 4 * b and arc["end_angle"] == 2 * beta and close(arc["a"], S["a"])
        assert close(abs(2 * b - arc["centre"]), arc["radius"], 1e-13)                                     # the arc passes through z = +-2b
        tang = np.angle((S["z"][1] - S["z"][0]) * np.exp(-1j * PI))                                        # direction of the surface leaving the tail
        assert abs(tang + 2 * beta) < 0.02                                                                 # end tangent at 2 beta to the chord
        # the arc is traced twice: the lower half of the circle lands on the same arc
        up, lo = S["z"][S["theta"] > 0], S["z"][S["theta"] < 0]
        assert lo.size > 10 and np.max(np.abs(np.abs(lo - arc["centre"]) - arc["radius"])) < 1e-12 * b and lo.imag.min() > -1e-14
    for e in (0.04, 0.12):                                                                                 # the symmetric section (N25)
        S = ch14.zhukhovsky_airfoil(b, e, 0.0, n)
        assert close(S["z"][1:], np.conj(S["z"][1:][::-1]), 1e-12, 1e-13)                                  # mirror image in the x-axis
        assert S["y"][5] > 0 and close(S["x"].min(), -(b * (1 + 2 * e) + b / (1 + 2 * e)), 1e-13)          # upper surface first; the nose
    Sa = ch14.zhukhovsky_airfoil(1.0, 0.0, 0.1, 100, a=1.1)                                                # (a, b, beta) parametrisation
    assert Sa["a"] == 1.1 and close(Sa["e"], 1.1 * math.cos(0.1) - 1) and abs(Sa["z"][0] - 2.0) < 1e-14


def test_zhukhovsky_shape_V7_rules_of_14_8_are_first_order_and_slip_4():  # V7 + V1 — N26 (14.8), slip #4
    assert ch14.THICKNESS_RULE_COEFFICIENT == 3 * math.sqrt(3) / 4 and abs(AF.THICKNESS_RULE_COEFFICIENT - 1.3) < 1e-3   # printed, rounded
    rule = ch14.zhukhovsky_shape(0.5, 0.06, 0.08, exact=False)
    assert rule["exact"] is False and close(rule["chord"], 4 * 0.5) and close(rule["camber"], 0.5 * 0.08 * 2.0)
    assert close(rule["thickness_ratio"], 3 * math.sqrt(3) / 4 * 0.06) and close(rule["t_max"], rule["thickness_ratio"] * 2.0)
    assert rule["x_tmax"] == -0.5 and rule["x_tmax_over_c"] == 0.25 and close(rule["camber_ratio"], 0.04)
    # the measured contour tends to the rules as e -> 0: coefficient 3 sqrt(3)/4, position x = -b, chord 4b (1 + e^2)
    es = np.array([0.1, 0.05, 0.02, 0.01, 0.005])
    coef, xt, dc = [], [], []
    for e in es:
        m = ch14.zhukhovsky_shape(1.0, e, 0.0)
        chord, tmax, x_t, t_star = _exact_symmetric(e)
        assert close(m["chord"], chord, 1e-9) and close(m["t_max"], tmax, 2e-7)
        coef.append(m["thickness_ratio"] / e)
        xt.append(m["x_tmax"])
        dc.append(m["chord"] / 4 - 1)
    coef, xt, dc = np.array(coef), np.array(xt), np.array(dc)
    assert np.all(coef < 3 * math.sqrt(3) / 4) and np.all(np.diff(coef) > 0)             # the rule over-predicts, less and less
    assert abs(observed_order(es, 3 * math.sqrt(3) / 4 - coef) - 1.0) < 0.15             # ... with a first-order error in e
    assert abs(coef[-1] - 3 * math.sqrt(3) / 4) < 0.01 * 1.3                             # within 1 % at e = 0.005
    assert 0.085 < 1 - coef[0] / (3 * math.sqrt(3) / 4) < 0.10                           # slip #4: about 9 % low at e = 0.1
    assert abs(1.3 - 3 * math.sqrt(3) / 4) < 0.05 * (3 * math.sqrt(3) / 4 - coef[2])     # 1.3 vs 1.299: 25 times below the rule's own error at e = 0.02
    assert abs(observed_order(es, dc) - 2.0) < 0.15 and close(dc, es ** 2 * (1 - 2 * es + 4 * es ** 2), 0.02)   # chord correction: second order
    assert np.all(np.abs(xt + 1.0) < 0.2 * es) and abs(xt[-1] + 1.0) < 1e-3              # thickest near x = -b: the quarter chord
    lim = ch14.zhukhovsky_shape(1.0, 1e-3, 0.0)
    assert abs(lim["thickness_ratio"] / 1e-3 - 3 * math.sqrt(3) / 4) < 2e-3 and abs(lim["x_tmax_over_c"] - 0.25) < 1e-3
    # camber: exactly 2b tan(beta) for the arc; the rule camber = beta c / 2 is first order in beta; thickness changes it little
    for beta in (0.02, 0.05, 0.1):
        m = ch14.zhukhovsky_shape(1.0, 0.0, beta)
        assert close(m["camber_ratio"], math.tan(beta) / 2, 1e-7) and abs(m["camber_ratio"] / (beta / 2) - 1) < 0.4 * beta ** 2
        assert close(ch14.camber_from_beta(beta, m["chord"]), 0.5 * beta * m["chord"]) and close(ch14.zero_lift_angle(beta / 2), -beta)
        mt = ch14.zhukhovsky_shape(1.0, 0.05, beta)
        assert abs(mt["camber_ratio"] / m["camber_ratio"] - 1) < 0.01 and close(mt["e"], 0.05, 1e-12)
    # camber changes the thickness coefficient only at second order in beta
    c0 = ch14.zhukhovsky_shape(1.0, 0.01, 0.0)["thickness_ratio"]
    rise = [ch14.zhukhovsky_shape(1.0, 0.01, be)["thickness_ratio"] / c0 - 1 for be in (0.05, 0.1, 0.2)]
    assert abs(observed_order([0.05, 0.1, 0.2], rise) - 2.0) < 0.2 and rise[-1] < 0.04


def test_trailing_edge_V7_cusp_angle_tends_to_zero_and_a_wedge_keeps_its_angle():  # V7 + V3 — N27, D05, R07
    ns, ang = [100, 200, 400, 800, 1600], []
    for n in ns:
        S = ch14.zhukhovsky_airfoil(1.0, 0.08, 0.07, n)
        ang.append(ch14.trailing_edge_angle(S["x"], S["y"], te_index=0))
    assert abs(observed_order([1.0 / n for n in ns], ang) - 1.0) < 0.05 and ang[-1] < 1e-3       # -> 0, first order in the spacing
    assert close(ang[-1], ch14.trailing_edge_angle(S["x"], S["y"]), 1e-12)                       # default: the point of largest x
    tri = np.array([1.0 + 0j, -1.0 + 0.3j, -1.0 - 0.3j])                                         # a polygon with a known vertex angle
    assert close(ch14.trailing_edge_angle(tri.real, tri.imag), 2 * math.atan(0.3 / 2.0), 1e-13)
    # the corner-flow model of a finite edge angle (R07)
    for n in (1.0, 1.05, 1.5, 2.0):
        wd = ch14.trailing_edge_wedge(n, U=2.0, d=0.5)
        assert close(wd["included_angle"], 2 * PI * (1 - 1 / n), 1e-13, 1e-15) and close(wd["half_fluid_angle"], PI / n) and wd["n"] == n
        r = np.array([0.01, 0.1, 0.5, 2.0])
        for th in (PI / n, -PI / n):
            assert np.max(np.abs(wd["velocity"](r, th)[1])) < 1e-13                               # the two surfaces are streamlines
        assert close(wd["speed"](r), 2.0 * (r / 0.5) ** (n - 1)) and close(wd["Cp"](r), 1 - (r / 0.5) ** (2 * (n - 1)))
        ur, ut = wd["velocity"](r, 0.3)
        assert close(np.hypot(ur, ut), wd["speed"](r), 1e-13)
        phi = lambda rr, tt: 2.0 * 0.5 / n * (rr / 0.5) ** n * np.cos(n * tt)                     # phi = (U d / n)(r/d)^n cos(n theta)
        h = 1e-6
        assert close(ur, (phi(r + h, 0.3) - phi(r - h, 0.3)) / (2 * h), 1e-7) and close(ut, (phi(r, 0.3 + h) - phi(r, 0.3 - h)) / (2 * h * r), 1e-7)
        assert wd["Cp_edge"] == (0.0 if n == 1.0 else 1.0) and close(wd["Cp"](0.0), wd["Cp_edge"], 0, 1e-15)
    assert close(ch14.trailing_edge_wedge(1.0, U=2.0)["velocity"](0.7, 0.4), (2.0 * math.cos(0.4), -2.0 * math.sin(0.4)))   # n = 1: uniform stream


def test_circular_arc_V2_derivation():  # V2 — D03 (three-star): (14.4) -> (14.5) -> cosine rule -> (14.6) -> the arc; each line of Part F
    b, R, th, be = sp.symbols("b R theta beta", positive=True)
    zeta = R * sp.exp(sp.I * th)
    z = zeta + b ** 2 / zeta                                                             # step 1: (14.3) with zeta = R e^{i theta}
    x, y = sp.re(sp.expand(z, complex=True)), sp.im(sp.expand(z, complex=True))
    assert sp.simplify(x - (R + b ** 2 / R) * sp.cos(th)) == 0 and sp.simplify(y - (R - b ** 2 / R) * sp.sin(th)) == 0      # step 2: (14.4)
    assert sp.simplify(x ** 2 / sp.cos(th) ** 2 - y ** 2 / sp.sin(th) ** 2 - 4 * b ** 2) == 0                               # step 4
    assert sp.simplify(x ** 2 * sp.sin(th) ** 2 - y ** 2 * sp.cos(th) ** 2 - 4 * b ** 2 * sp.sin(th) ** 2 * sp.cos(th) ** 2) == 0   # step 5: (14.5)
    # steps 6-7: the circle |zeta - i b tan(beta)| = b / cos(beta), written out, IS the cosine rule of triangle OQP
    circle = sp.expand((zeta - sp.I * b * sp.tan(be)) * sp.conjugate(zeta - sp.I * b * sp.tan(be)), complex=True) - b ** 2 / sp.cos(be) ** 2
    cosrule = R ** 2 + b ** 2 * sp.tan(be) ** 2 - 2 * R * b * sp.tan(be) * sp.sin(th) - b ** 2 / sp.cos(be) ** 2
    assert sp.simplify(circle - cosrule) == 0                                            # the angle at O is 90 deg - theta: cos -> sin(theta)
    assert sp.simplify(cosrule - (R ** 2 - b ** 2 - 2 * R * b * sp.tan(be) * sp.sin(th))) == 0                             # step 8
    assert sp.simplify(cosrule / R - ((R - b ** 2 / R) - 2 * b * sp.tan(be) * sp.sin(th))) == 0                            # step 9: (14.6)
    Y = 2 * b * sp.tan(be) * sp.sin(th) ** 2                                             # y = 2 b tan(beta) sin^2(theta), from (14.6) and (14.4)
    Rsol = b * sp.tan(be) * sp.sin(th) + b * sp.sqrt(sp.tan(be) ** 2 * sp.sin(th) ** 2 + 1)                               # the circle, R(theta) > 0
    assert sp.simplify(cosrule.subs(R, Rsol)) == 0 and sp.simplify((y - Y).subs(R, Rsol)) == 0
    # steps 10-12: eliminate theta
    X2 = 4 * b ** 2 * sp.cos(th) ** 2 + Y ** 2 * sp.cos(th) ** 2 / sp.sin(th) ** 2       # step 11: (14.5) divided by sin^2
    assert sp.simplify((x ** 2 - X2).subs(R, Rsol)) == 0
    step12 = X2 + Y ** 2 + 2 * b * Y * (sp.cot(be) - sp.tan(be)) - 4 * b ** 2
    assert sp.simplify(step12) == 0
    assert sp.simplify(sp.expand_trig(sp.cot(be) - sp.tan(be) - 2 * sp.cot(2 * be))) == 0     # step 13: the double-angle identity
    arc = X2 + (Y + 2 * b * sp.cot(2 * be)) ** 2 - (2 * b / sp.sin(2 * be)) ** 2         # step 14: the circle of the arc
    assert sp.simplify(arc) == 0
    assert sp.simplify(Y.subs(th, sp.pi / 2) - 2 * b * sp.tan(be)) == 0                  # step 15: height 2 b tan(beta)
    assert sp.series(2 * b * sp.tan(be) / (4 * b), be, 0, 3).removeO() == be / 2         # camber/c = beta/2 + O(beta^3): the rule of (14.8)
    # the design's worked number (b = 1, tan(beta) = 3/4: the top zeta = 2i maps to 1.5i) and the coded circle
    assert complex(ch14.joukowski(2j, 1.0)) == 1.5j
    arcn = ch14.circular_arc(1.0, math.atan(0.75))
    assert close(arcn["height"], 1.5) and close(arcn["centre"].imag + arcn["radius"], 1.5) and close(arcn["centre"].imag, -7 / 12, 1e-12)
    # the chapter's own engine agrees (tested separately from the derivation above)
    eng = ch14.circular_arc_sym()
    assert eng["ok"] is True and eng["all_zero"] is True and eng["res_14_5"] == 0 and eng["res_14_6"] == 0 and eng["res_arc"] == 0


def test_thickness_rule_V2_derivation():  # V2 — D04 (three-star): a = b(1 + e) (14.7) -> c = 4b, t_max/c = (3 sqrt(3)/4) e at x = -b
    b, e = sp.symbols("b e", positive=True)
    th = sp.symbols("theta", real=True)
    zeta = -b * e + b * (1 + e) * sp.exp(sp.I * th)                                      # step 1
    eps = e * (1 - sp.exp(-sp.I * th))
    assert sp.simplify(zeta - b * sp.exp(sp.I * th) * (1 + eps)) == 0                    # step 2: zeta = b e^{i theta}(1 + eps)
    assert sp.simplify(zeta.subs(th, 0) - b) == 0                                        # the circle passes through zeta = b
    z = zeta + b ** 2 / zeta
    z1 = sp.series(z, e, 0, 2).removeO()                                                 # first order in e
    step4 = 2 * b * sp.cos(th) + 2 * sp.I * b * eps * sp.sin(th)
    assert sp.simplify(sp.expand((z1 - step4).rewrite(sp.exp))) == 0                     # steps 3-4
    x1 = sp.simplify(sp.re(sp.expand(z1.rewrite(sp.cos), complex=True)))
    y1 = sp.simplify(sp.im(sp.expand(z1.rewrite(sp.cos), complex=True)))
    assert sp.simplify(x1 - (2 * b * sp.cos(th) - 2 * b * e * sp.sin(th) ** 2)) == 0     # step 5
    assert sp.simplify(y1 - 2 * b * e * sp.sin(th) * (1 - sp.cos(th))) == 0
    assert x1.subs(th, 0) == 2 * b and sp.simplify(x1.subs(th, sp.pi) + 2 * b) == 0      # step 6: c = 4b to first order
    assert sp.simplify(y1.subs(th, -th) + y1) == 0 and sp.simplify(x1.subs(th, -th) - x1) == 0      # step 7: symmetric, t = 2y
    dy = sp.diff(y1, th)
    assert sp.simplify(dy - 2 * b * e * (sp.cos(th) - sp.cos(2 * th))) == 0              # step 8
    assert sp.simplify(dy.subs(th, 2 * sp.pi / 3)) == 0 and sp.diff(y1, th, 2).subs(th, 2 * sp.pi / 3) < 0   # step 9: a maximum at 120 deg
    tmax = 2 * y1.subs(th, 2 * sp.pi / 3)
    assert sp.simplify(tmax - 3 * sp.sqrt(3) * b * e) == 0 and sp.simplify(tmax / (4 * b) - 3 * sp.sqrt(3) / 4 * e) == 0   # step 10
    assert sp.simplify(x1.subs(th, 2 * sp.pi / 3).subs(e, 0) + b) == 0                   # step 11: x = -b, the quarter chord
    # the exact chord and its expansion ("c = 4b(1 + e^2 + ...)" in step 6)
    chord = z.subs(th, 0) - z.subs(th, sp.pi)
    assert sp.simplify(chord - (2 * b + b * (1 + 2 * e) + b / (1 + 2 * e))) == 0
    assert sp.expand(sp.series(chord / (4 * b), e, 0, 3).removeO()) == 1 + e ** 2
    # the tail: y ~ b e theta^3 and 2b - x ~ b theta^2 (the cusp of D05)
    assert sp.series(y1, th, 0, 4).removeO() == b * e * th ** 3 and sp.series(2 * b - x1, th, 0, 3).removeO().subs(e, 0) == b * th ** 2
    # numbers of the Check: e = 0.1 -> rule 0.1299, measured 0.1179, chord 4.033 b
    m = ch14.zhukhovsky_shape(1.0, 0.1, 0.0)
    assert round(float(3 * sp.sqrt(3) / 4 * 0.1), 4) == 0.1299 and round(m["thickness_ratio"], 4) == 0.1179 and round(m["chord"], 3) == 4.033
    assert close(ch14.zhukhovsky_shape(1.0, 0.1, 0.0, exact=False)["thickness_ratio"], float(3 * sp.sqrt(3) / 4) * 0.1)


def test_cusp_V2_derivation():  # V2 — D05 (two-star): dz/dzeta = 0 at zeta = b, z - 2b = (zeta - b)^2 / b, angles double
    b, d = sp.symbols("b delta", positive=True)
    zeta = sp.symbols("zeta")
    z = zeta + b ** 2 / zeta
    assert sp.diff(z, zeta).subs(zeta, b) == 0 and sp.simplify(sp.diff(z, zeta) - (1 - b ** 2 / zeta ** 2)) == 0
    assert sp.diff(z, zeta, 2).subs(zeta, b) == 2 / b and sp.diff(z, zeta).subs(zeta, -b) == 0             # the two critical points
    ser = sp.series(z.subs(zeta, b + d), d, 0, 3).removeO()
    assert sp.simplify(ser - (2 * b + d ** 2 / b)) == 0                                                    # z - 2b = delta^2 / b
    phi = sp.symbols("phi", real=True)
    arg_out = sp.arg(sp.exp(sp.I * phi) ** 2)                                                              # a ray at phi leaves at 2 phi
    assert sp.simplify(sp.exp(sp.I * phi) ** 2 - sp.exp(2 * sp.I * phi)) == 0 and arg_out is not None
    # so the two tangents of the circle at B (pi apart) leave the image 2 pi apart: zero included angle (numerically, on our section)
    S = ch14.zhukhovsky_airfoil(1.0, 0.08, 0.07, 20000)
    v_up, v_lo = S["z"][1] - S["z"][0], S["z"][-1] - S["z"][0]
    assert abs(np.angle(v_up / v_lo)) < 1e-3 and abs(np.angle(-v_up) + 2 * 0.07) < 1e-3                    # a cusp pointing along e^{-2 i beta}
    # a finite edge angle: w = A z^n, fluid in |theta| < pi/n, body angle 2 pi (1 - 1/n), speed ~ r^(n-1)
    n, r, A = sp.symbols("n r A", positive=True)
    th = sp.symbols("theta", real=True)
    psi = A * r ** n * sp.im(sp.exp(sp.I * n * th))                                                        # Im of w = A (r e^{i theta})^n
    assert psi == A * r ** n * sp.sin(n * th)
    assert sp.simplify(psi.subs(th, sp.pi / n)) == 0 and sp.simplify(psi.subs(th, -sp.pi / n)) == 0        # the walls are streamlines
    speed = sp.simplify(sp.Abs(sp.diff(A * sp.Symbol("zz", positive=True) ** n, sp.Symbol("zz", positive=True))))
    assert sp.simplify(speed - A * n * sp.Symbol("zz", positive=True) ** (n - 1)) == 0                     # zero at the edge for n > 1
    assert sp.simplify((2 * sp.pi - 2 * sp.pi / n) - 2 * sp.pi * (1 - 1 / n)) == 0


# =====================================================================================================================
# 6. C05 — the flow round the airfoil (14.9); N28, R08
# =====================================================================================================================
def _outside_points(F, n=9):
    """A field of z-plane points outside the section (rays from points of the contour outward in the zeta-plane)."""
    t = np.linspace(0.2, 6.0, 4 * n)
    rho = np.array([1.02, 1.1, 1.5, 3.0, 9.0])
    zeta = F.mu + F.a * rho[:, None] * np.exp(1j * t[None, :])
    return np.asarray(ch14.joukowski(zeta, F.b)), zeta


def test_zhukhovsky_flow_V1_section_is_a_streamline_far_field_and_complex_derivative():  # V1 — C05 (14.9), N28
    for e, beta in SECTIONS[1:]:
        for al, frac in ((0.12, 1.0), (-0.07, 1.0), (0.2, 0.4)):
            F0 = ch14.zhukhovsky_flow(2.1, 0.5, e, beta, al)
            F = ch14.zhukhovsky_flow(2.1, 0.5, e, beta, al, Gamma_cw=frac * F0.Gamma_kutta)
            t = np.linspace(-beta + 0.05, -beta + 2 * PI - 0.05, 400)
            t = t[np.abs((t + PI) % (2 * PI) - PI) < 3.0]                 # stay off the cut of the logarithm
            zs = F.mu + F.a * np.exp(1j * t)
            assert np.ptp(np.asarray(F.w_zeta(zs)).imag) < 1e-12 * F.U * F.a                       # psi constant on the circle = the section
            zsec = np.asarray(ch14.joukowski(zs * (1 + 1e-9), F.b))
            assert np.ptp(F.stream(zsec.real, zsec.imag)) < 1e-6 * F.U * F.a                       # and through the inverse map
            z, zeta = _outside_points(F)
            assert close(F.zeta_of(z), zeta, 1e-10)
            q = np.asarray(F.dwdz(z))
            assert close(q, np.asarray(F.dwdzeta(zeta)) * zeta ** 2 / (zeta ** 2 - F.b ** 2), 1e-10)   # (14.9) as printed
            assert close(q, np.asarray(F.dwdz_from_zeta(zeta)), 1e-10)
            h = 1e-6 * F.b                                                # dw/dz against centred differences of w(z), two directions
            keep = np.abs(np.angle(zeta - F.mu)) < 2.9
            for dz in (h, 1j * h):
                fd = (np.asarray(F.w(z + dz)) - np.asarray(F.w(z - dz))) / (2 * dz)
                assert close(fd[keep], q[keep], 2e-6, 1e-7 * F.U)         # analytic: the same derivative in every direction (Cauchy–Riemann)
            u, v = F.velocity(z.real, z.imag)
            assert close(u, q.real, 1e-14, 1e-300) and close(v, -q.imag, 1e-14, 1e-300)            # dw/dz = u - iv (trap T2)
            assert close(F.speed(z), np.abs(q), 1e-14) and close(F.cp(z), 1 - np.abs(q) ** 2 / F.U ** 2, 1e-12, 1e-13)
            assert close(F.cp(z), ch14.pressure_coefficient(F.speed(z), F.U), 1e-12, 1e-13)
            assert close(F.phi(z.real, z.imag), np.asarray(F.w(z)).real, 1e-14, 1e-300) and close(F.psi(z.real, z.imag), F.stream(z.real, z.imag))
            far = 1e6 * np.exp(1j * np.linspace(0.1, 6.0, 7))
            assert close(F.dwdz(far), F.U * np.exp(-1j * al), 1e-5) and close(F.u_inf, (F.U * math.cos(al), F.U * math.sin(al)))
            assert not np.any(F.inside(z.real, z.imag))
    F = ch14.zhukhovsky_flow(*SEC)
    assert bool(F.inside(0.0, 0.03)) and bool(F.inside(-0.3, 0.02)) and not bool(F.inside(0.0, 0.3)) and not bool(F.inside(0.6, 0.0))
    assert np.isnan(F.speed(0.0 + 0.03j)) and F.label == "Zhukhovsky airfoil"
    assert close(F.dwdz(1e7), complex(29.89, -2.61), 1e-3)                # design D06: far field 29.89 - 2.61 i for U = 30, alpha = 5 deg


def test_zhukhovsky_flow_V4_circulation_is_the_same_in_both_planes():  # V4 — C05, D08 (trap T6)
    tt = np.linspace(0.0, 2 * PI, 2000, endpoint=False)
    for e, beta in ((0.08, 0.07), (0.0, 0.15), (0.12, 0.0)):
        F0 = ch14.zhukhovsky_flow(2.0, 0.5, e, beta, 0.11)
        for G in (F0.Gamma_kutta, 0.0, 10.0, -3.0):
            F = ch14.zhukhovsky_flow(2.0, 0.5, e, beta, 0.11, Gamma_cw=G)
            for R, c in ((1.3, 0.0), (2.5, 0.2 + 0.1j), (40.0, -3.0j)):                           # circles in the z-plane round the section
                zc = c + R * np.exp(1j * tt)
                loop = np.sum(np.asarray(F.dwdz(zc)) * 1j * R * np.exp(1j * tt)) * (2 * PI / tt.size)
                assert abs(loop.real + G) < 1e-9 * max(abs(G), 1.0) and abs(loop.imag) < 1e-9     # circulation -Gamma_cw, no net source
            for rho_ in (1.05, 1.6, 4.0):                                                         # partner contours: a zeta-circle and its image
                zeta = F.mu + rho_ * F.a * np.exp(1j * tt)
                dzeta = 1j * rho_ * F.a * np.exp(1j * tt) * (2 * PI / tt.size)
                dz = dzeta * np.asarray(ch14.joukowski_derivative(zeta, F.b))
                ring_z = np.sum(np.asarray(F.dwdz_from_zeta(zeta)) * dz)
                ring_zeta = np.sum(np.asarray(F.dwdzeta(zeta)) * dzeta)
                assert abs(ring_z - ring_zeta) < 1e-11 * max(abs(G), 1.0) and abs(ring_zeta.real + G) < 1e-10 * max(abs(G), 1.0)
                assert np.max(np.abs(np.asarray(F.dwdz_from_zeta(zeta)) - np.asarray(F.dwdzeta(zeta)))) > 1e-3   # the velocity is NOT the same


def test_zhukhovsky_flow_V7_plate_limit_reproduces_the_chapter_6_route():  # V7 — C05/N34: e = beta = 0 through (6.69)
    U, b, al = 1.9, 0.7, 0.21
    F, Fp = ch14.zhukhovsky_flow(U, b, 0.0, 0.0, al), ch14.flat_plate_flow(U, b, al)
    assert F.mu == 0 and F.a == b and close(F.Gamma_kutta, 4 * PI * U * b * math.sin(al)) and Fp.Gamma_cw == F.Gamma_cw   # N34
    w6, dw6 = ch14.circle_flow_zeta(U, b, Gamma_cw=F.Gamma_cw, alpha=al)
    X, Y = np.meshgrid(np.linspace(-3, 3, 13), np.linspace(-2, 2, 9) + 0.11)
    z = X + 1j * Y
    zeta6 = np.asarray(ch14.joukowski_inverse(z, b, branch="outside"))
    assert close(F.zeta_of(z), zeta6, 1e-12) and close(Fp.dwdz(z), np.asarray(dw6(zeta6)) / np.asarray(ch14.joukowski_derivative(zeta6, b)), 1e-11)
    assert close(np.asarray(Fp.w(z)).imag, np.asarray(w6(zeta6)).imag, 1e-12, 1e-12)
    for G in (0.0, 1.0):
        assert ch14.flat_plate_flow(U, b, al, Gamma_cw=G).Gamma_cw == G
    # mirror symmetry of a symmetric section: alpha -> -alpha conjugates the flow
    Fa, Fb = ch14.zhukhovsky_flow(U, b, 0.1, 0.0, al), ch14.zhukhovsky_flow(U, b, 0.1, 0.0, -al)
    zz, _ = _outside_points(Fa)
    assert close(Fa.dwdz(zz), np.conj(np.asarray(Fb.dwdz(np.conj(zz)))), 1e-10) and close(Fa.Gamma_cw, -Fb.Gamma_cw)


def test_zhukhovsky_inverse_V1_round_trip_and_the_root_outside_the_shifted_circle():  # V1 — R05, trap T8
    for e, beta in SECTIONS:
        F = ch14.zhukhovsky_flow(1.0, 0.5, e, beta, 0.1)
        z, zeta = _outside_points(F)
        back = np.asarray(ch14.zhukhovsky_inverse(z, F.b, F.mu, F.a))
        assert close(back, zeta, 1e-10) and close(ch14.joukowski(back, F.b), z, 1e-12) and np.all(np.abs(back - F.mu) > F.a)
        assert close(ch14.zhukhovsky_inverse(z, F.b, F.mu), back)                                 # without the body test: same root
    X, Y = np.meshgrid(np.linspace(-3, 3, 25), np.linspace(-2, 2, 17) + 0.07)
    z = X + 1j * Y
    assert close(ch14.zhukhovsky_inverse(z, 0.6), ch14.joukowski_inverse(z, 0.6, branch="outside"), 1e-13)   # mu = 0: chapter 6
    F = ch14.zhukhovsky_flow(1.0, 0.5, 0.0, 0.3, 0.0)                                             # a raised circle: part of it lies inside |zeta| = b
    zeta_low = F.mu + 1.01 * F.a * np.exp(1j * np.linspace(-2.0, -1.2, 9))                        # just outside its lower side
    assert np.all(np.abs(zeta_low) < F.b)                                                         # ... yet inside the centred circle
    z_low = np.asarray(ch14.joukowski(zeta_low, F.b))
    assert close(ch14.zhukhovsky_inverse(z_low, F.b, F.mu, F.a), zeta_low, 1e-11)
    wrong = np.asarray(ch14.joukowski_inverse(z_low, F.b, branch="outside"))                      # planted: choose by |zeta| >= b
    assert np.all(np.abs(wrong - F.mu) < F.a) and np.max(np.abs(wrong - zeta_low)) > 0.1          # it lands inside the body
    inner = np.asarray(ch14.zhukhovsky_inverse(np.array([0.0 + 0.05j, 0.3 + 0.2j]), 0.5, *ch14.zhukhovsky_circle(0.5, 0.2, 0.0)))
    assert np.isnan(inner[0].real) and not np.isnan(inner[1].real)                                # NaN inside the section only
    assert type(ch14.zhukhovsky_inverse(3.0 + 1.0j, 0.5)) is complex


def test_offcentre_circle_flow_V2_derivation():  # V2 — D06 (two-star): shift and rotate (6.68)
    U, a, G = sp.symbols("U a Gamma", positive=True)
    al, t = sp.symbols("alpha t", real=True)
    mur, mui = sp.symbols("mu_r mu_i", real=True)
    mu = mur + sp.I * mui
    Z, zeta = sp.symbols("Z zeta")
    w0 = U * (Z + a ** 2 / Z) + sp.I * G / (2 * sp.pi) * sp.log(Z / a)                   # (6.68): centred, aligned
    w = w0.subs(Z, (zeta - mu) * sp.exp(-sp.I * al))                                     # replace Z by (zeta - mu) e^{-i alpha}
    target = U * ((zeta - mu) * sp.exp(-sp.I * al) + a ** 2 * sp.exp(sp.I * al) / (zeta - mu)) + sp.I * G / (2 * sp.pi) * sp.log((zeta - mu) / a)
    diff = sp.diff(w - target, zeta)
    assert sp.simplify(diff) == 0                                                        # same velocity: they differ by a constant only
    const = sp.simplify((w - target).subs({zeta: mu + 2 * a, al: sp.Rational(1, 5)}))
    assert sp.simplify(const - sp.I * G / (2 * sp.pi) * (-sp.I * sp.Rational(1, 5))) == 0          # the constant Gamma alpha / 2 pi (real)
    dw = sp.diff(target, zeta)
    assert sp.limit(dw.subs(zeta, 1 / sp.Symbol("q", positive=True)), sp.Symbol("q", positive=True), 0) == U * sp.exp(-sp.I * al)   # far field
    on = target.subs(zeta, mu + a * sp.exp(sp.I * t))                                    # on the circle
    psi = sp.im(sp.expand(sp.simplify(on.subs(sp.log(sp.exp(sp.I * t)), sp.I * t)), complex=True))
    assert sp.simplify(psi) == 0                                                         # Im w = 0: the circle is a streamline
    res = sp.residue(dw, zeta, mu)                                                       # circulation: 2 pi i x residue = -Gamma
    assert sp.simplify(2 * sp.pi * sp.I * res + G) == 0
    assert sp.simplify(target.subs({mur: 0, mui: 0, al: 0}) - w0.subs(Z, zeta)) == 0     # mu = 0, alpha = 0 returns (6.68)
    assert sp.simplify(sp.diff(target, zeta) * zeta ** 2 / (zeta ** 2 - sp.Symbol("b") ** 2)
                       - sp.diff(target, zeta) / sp.diff(zeta + sp.Symbol("b") ** 2 / zeta, zeta)) == 0     # D07: (14.9) is the chain rule


def test_circulation_invariance_V2_derivation():  # V2 — D08 (two-star): (dw/dz) dz = (dw/dzeta) dzeta
    b, G, U, a = sp.symbols("b Gamma U a", positive=True)
    zeta, mu, al = sp.symbols("zeta mu alpha")
    w = sp.Function("w")(zeta)
    z = zeta + b ** 2 / zeta
    dwdz = sp.diff(w, zeta) / sp.diff(z, zeta)                                           # chain rule (R08)
    assert sp.simplify(dwdz * sp.diff(z, zeta) - sp.diff(w, zeta)) == 0                  # the integrand of the loop integral is unchanged
    assert sp.simplify(dwdz - sp.diff(w, zeta) * zeta ** 2 / (zeta ** 2 - b ** 2)) == 0  # ... although the velocity is rescaled (14.9)
    wz = U * ((zeta - mu) * sp.exp(-sp.I * al) + a ** 2 * sp.exp(sp.I * al) / (zeta - mu)) + sp.I * G / (2 * sp.pi) * sp.log((zeta - mu) / a)
    assert sp.simplify(2 * sp.pi * sp.I * sp.residue(sp.diff(wz, zeta), zeta, mu) + G) == 0        # loop integral of dw = -Gamma in either plane
    assert sp.limit(sp.diff(z, zeta), zeta, sp.oo) == 1                                  # z = zeta far away: the same far stream


# =====================================================================================================================
# 7. C06 — the lift of a Zhukhovsky airfoil (14.11)–(14.12); N29, N31, N32, N34, N18, R10
# =====================================================================================================================
def _slope_factor(e):
    """Exact lift slope of the symmetric Zhukhovsky section over 2 pi: 8 pi (a/c) / 2 pi with our closed-form chord."""
    return 4 * (1 + e) / (3 + 2 * e + 1 / (1 + 2 * e))


def test_zhukhovsky_lift_V1_blasius_force_is_rho_U_gamma_normal_to_the_stream():  # V1 — C06, R04, R10 (6.60), (6.62)
    U, b, rho = 2.4, 0.5, 1.21
    for e, beta in SECTIONS:
        for al in ANGLES:
            F = ch14.zhukhovsky_flow(U, b, e, beta, al)
            for G in (None, 0.0, 0.6 * F.Gamma_kutta + 0.3, -1.1):
                f = ch14.zhukhovsky_force(U, b, e, beta, al, rho, Gamma_cw=G)
                Gam = F.Gamma_kutta if G is None else G
                scale = rho * U * max(abs(Gam), U * b)
                assert abs(f["L"] - rho * U * Gam) < 1e-10 * scale and abs(f["D"]) < 1e-10 * scale          # L = rho U Gamma, no drag
                assert close(f["L"], ch14.lift_per_span(rho, U, Gamma_cw=Gam), 0, 1e-10 * scale)
                assert abs(f["Fx"] + f["L"] * math.sin(al)) < 1e-10 * scale and abs(f["Fy"] - f["L"] * math.cos(al)) < 1e-10 * scale   # lift is not F_y (T5)
                assert f["Gamma_cw"] == Gam and f["method"] == "blasius"
                assert close(f["CL"], 2 * Gam / (U * f["chord"]), 1e-9, 1e-10) and abs(f["CD"]) < 1e-9
                # an independent contour: our own Blasius integral on a different circle
                tt = np.linspace(0, 2 * PI, 720, endpoint=False)
                zc = 0.3 + 0.2j + 5.0 * np.exp(1j * tt)
                Fl = ch14.zhukhovsky_flow(U, b, e, beta, al, Gamma_cw=Gam)
                DiL = 0.5j * rho * np.sum(np.asarray(Fl.dwdz(zc)) ** 2 * 1j * 5.0 * np.exp(1j * tt)) * (2 * PI / tt.size)   # D_x - i L_y
                assert abs(DiL.real - f["Fx"]) < 1e-9 * scale and abs(-DiL.imag - f["Fy"]) < 1e-9 * scale
            fk = ch14.zhukhovsky_force(U, b, e, beta, al, rho)
            if abs(al + beta) > 1e-9:
                assert np.sign(fk["L"]) == np.sign(al + beta)                                               # clockwise Gamma > 0 lifts (T1)
            assert close(fk["CL"], ch14.zhukhovsky_lift_coefficient(al, beta, e, exact=True), 1e-8, 1e-10)  # 8 pi (a/c) sin(alpha + beta)
            assert close(fk["chord"], ch14.zhukhovsky_shape(b, e, beta, n=4001)["chord"], 1e-13)
    for al in (0.02, 0.1, 0.35, -0.2):                                                                      # N34: the flat plate, exactly
        fp = ch14.zhukhovsky_force(U, b, 0.0, 0.0, al, rho)
        assert close(fp["CL"], 2 * PI * math.sin(al), 1e-10) and close(ch14.zhukhovsky_lift_coefficient(al, 0.0, 0.0, exact=True), 2 * PI * math.sin(al), 1e-10)
        assert close(fp["Gamma_cw"], 4 * PI * U * b * math.sin(al), 1e-13) and close(fp["chord"], 4 * b, 1e-10)
    with pytest.raises(ValueError):
        ch14.zhukhovsky_force(U, b, 0.05, 0.05, 0.1, rho, method="panel")


def test_zhukhovsky_lift_V3_pressure_route_converges_to_the_blasius_force_at_second_order():  # V3 + V4 — C06/C07, R10 (two routes)
    rho = RHO
    fb = ch14.zhukhovsky_force(*SEC, rho)
    ns, eL, eD, eM = [250, 500, 1000, 2000, 4000], [], [], []
    for n in ns:
        fp = ch14.zhukhovsky_force(*SEC, rho, method="pressure", n=n)
        eL.append(abs(fp["L"] - fb["L"]) / fb["L"])
        eD.append(abs(fp["D"]) / fb["L"])
        eM.append(abs(fp["M_le"] - fb["M_le"]) / abs(fb["M_le"]))
        assert fp["method"] == "pressure"
    h = [1.0 / n for n in ns]
    for err in (eL, eD, eM):
        assert abs(observed_order(h, err) - 2.0) < 0.15, (observed_order(h, err), err)
    assert eL[-1] < 1e-4 and eD[-1] < 1e-4 and eM[-1] < 1e-4                                 # analysis T7: 1e-4 at n = 4000
    assert close(fb["L"], rho * US * ch14.kutta_circulation(US, ch14.zhukhovsky_circle(B, E, BETA)[1], ALPHA, BETA), 1e-10)
    # the moment: flat plate -> the line of action passes through the quarter chord, for every incidence
    for al in (0.03, 0.15, -0.1):
        f = ch14.zhukhovsky_force(1.7, 0.3, 0.0, 0.0, al, 1.0)
        assert abs(f["xcp_over_c"] - 0.25) < 1e-10 and close(f["x_cp"], -2 * 0.3 + 0.25 * 4 * 0.3, 1e-10)
        assert close(f["M_le"], f["Fy"] * 0.25 * f["chord"], 1e-9)       # counterclockwise +: an upward force behind the nose turns it nose-down
        assert close(f["M0"], f["x_cp"] * f["Fy"], 1e-12)
    # a cambered section: moment about the leading edge from the two routes, and M0 - M_le consistent with the force
    sh = ch14.zhukhovsky_shape(B, E, BETA, n=4001)
    A = ch14.zhukhovsky_airfoil(B, E, BETA, 4001)
    le = ch14.section_geometry(A["x"], A["y"], te_index=0)["le"]
    assert close(fb["M_le"], fb["M0"] - (le.real * fb["Fy"] - le.imag * fb["Fx"]), 1e-12) and close(fb["chord"], sh["chord"])


def test_zhukhovsky_lift_V5_slope_two_pi_per_radian_and_the_quarter_chord():  # V5 — C06 (14.12) against the cited thin-airfoil results
    bm = bench()
    slope = bm["thin_airfoil_lift_slope"]["slope_over_pi"] * PI                               # "2 pi units per radian"
    al = np.linspace(-0.1, 0.2, 7)
    rule = ch14.zhukhovsky_lift_coefficient(al, 0.05)
    assert close(np.diff(rule) / np.diff(al), slope, 1e-12)                                   # (14.12): slope 2 pi
    assert close(rule, slope * (al + 0.05), 1e-13, 1e-15) and close(ch14.zhukhovsky_lift_coefficient(1.0, 0.0), slope)   # parity pin (C.6 item 3)
    cl0 = ch14.zhukhovsky_lift_coefficient(0.0, 0.05)
    assert close(rule, cl0 + slope * al, 1e-13, 1e-15)                                        # the source's cambered form c_l = c_l0 + 2 pi alpha
    assert abs(ch14.zhukhovsky_lift_coefficient(1e-4, 0.0, 0.0, exact=True) / 1e-4 - slope) < 1e-6 * slope      # exact plate slope at small alpha
    assert close(ch14.thin_airfoil(lambda x: 0.0 * x, 0.07)["CL"], slope * 0.07, 1e-12) and ch14.thin_airfoil(lambda x: 0.0 * x, 0.07)["slope"] == slope
    xcp = bm["thin_airfoil_centre_of_pressure_symmetric"]["value"]                            # "one quarter of the chord behind the leading edge"
    assert abs(ch14.zhukhovsky_force(1.0, 0.25, 0.0, 0.0, 0.08, 1.0)["xcp_over_c"] - xcp) < 1e-10
    assert abs(ch14.thin_airfoil(lambda x: 0.0 * x, 0.07)["CM_c4"]) < 1e-14                   # no moment about the quarter chord: the same statement


def test_zhukhovsky_lift_V7_small_angle_rule_and_slip_8_thickness_raises_the_slope():  # V7 — C06, N31, slip #8
    for e in (0.0, 0.02, 0.05, 0.1, 0.2):
        exact = ch14.zhukhovsky_lift_coefficient(1e-3, 0.0, e, exact=True)
        assert close(exact / (2 * PI * 1e-3), _slope_factor(e), 1e-6)                         # our closed form of 8 pi a / c
    assert _slope_factor(0.0) == 1.0 and 1.085 < _slope_factor(0.1) < 1.095                   # about 9 % above 2 pi at e = 0.1
    assert abs(observed_order([0.1, 0.05, 0.02], [_slope_factor(e) - 1 for e in (0.1, 0.05, 0.02)]) - 1.0) < 0.1   # c = 4a is a first-order statement
    for beta in (0.0, 0.05, 0.12):
        assert ch14.zhukhovsky_lift_coefficient(-beta, beta) == 0.0 and abs(ch14.zhukhovsky_lift_coefficient(-beta, beta, 0.08, exact=True)) < 1e-15
        assert close(ch14.zero_lift_angle(math.tan(beta) / 2), -math.tan(beta)) and abs(ch14.zero_lift_angle(beta / 2) + beta) < 1e-15   # N31
        for al in (0.01, 0.05):                                                               # the arc: exact 2 pi sin(alpha + beta) / cos(beta)
            assert close(ch14.zhukhovsky_lift_coefficient(al, beta, 0.0, exact=True), 2 * PI * math.sin(al + beta) / math.cos(beta), 1e-7)
    x = np.array([0.2, 0.1, 0.05, 0.025])                                                     # rule/exact -> 1 as alpha + beta -> 0 (third-order sine error)
    gap = [abs(ch14.zhukhovsky_lift_coefficient(v, 0.0) / ch14.zhukhovsky_lift_coefficient(v, 0.0, 0.0, exact=True) - 1) for v in x]
    assert abs(observed_order(x, gap) - 2.0) < 0.05 and close(gap, x ** 2 / 6, 0.02)
    CL1 = ch14.zhukhovsky_lift_coefficient(np.deg2rad(2.0), 0.03)                             # N35: C_L(a2) = C_L(a1) + 2 pi (a2 - a1)
    assert close(ch14.zhukhovsky_lift_coefficient(np.deg2rad(6.0), 0.03), CL1 + 2 * PI * np.deg2rad(4.0))
    with pytest.raises(ValueError):
        ch14.zhukhovsky_lift_coefficient(0.1, np.array([0.0, 0.1]), 0.05, exact=True)


def test_trailing_edge_speed_V1_finite_limit_of_14_9_with_kutta_and_blow_up_without():  # V1 + V7 — N29 (trap T8)
    F = ch14.zhukhovsky_flow(*SEC)
    lim = ch14.trailing_edge_speed(F.U, F.a, F.b, F.alpha, F.beta)
    assert close(lim, F.U * F.b / F.a * math.cos(F.alpha + F.beta)) and 0 < lim < F.U
    d = np.array([1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8])
    err = np.array([abs(ch14.edge_speed(*SEC, float(q)) / lim - 1) for q in d])
    assert abs(observed_order(d, err) - 0.5) < 0.02                       # the speed approaches its limit like the square root of the distance
    assert err[-1] < 1e-4 and err[3] < 3e-4                               # 2.2e-4 at 1e-6 m, 2.2e-5 at 1e-8 m (the docstring states these)
    assert close(ch14.edge_speed(*SEC, 1e-5), F.speed_near_edge(1e-5), 1e-14) and close(F.speed_near_edge(d), [ch14.edge_speed(*SEC, float(q)) for q in d])
    q_edge = complex(F.dwdz_from_zeta(F.b))
    assert close(q_edge, lim * np.exp(2j * F.beta), 1e-13)                # u - iv = q e^{2 i beta}: the stream leaves along e^{-2 i beta}
    assert close(complex(F.dwdz(2 * F.b + 1e-9 * np.exp(-2j * F.beta))), q_edge, 1e-4)
    for frac in (0.8, 1.2, 0.0):                                          # any other circulation: |dw/dz| ~ d^(-1/2)
        G = frac * F.Gamma_kutta
        sp_ = np.array([ch14.edge_speed(*SEC, float(q), Gamma_cw=G) for q in d[2:]])
        assert abs(observed_order(d[2:], sp_) + 0.5) < 0.01 and sp_[-1] > 100 * F.U
        Fn = ch14.zhukhovsky_flow(*SEC, Gamma_cw=G)
        assert np.isinf(complex(Fn.dwdz_from_zeta(Fn.b)).real) and not Fn.is_kutta
        assert ch14.surface_point(*SEC, -BETA, Gamma_cw=G)["speed"] == float("inf")
    p = ch14.surface_point(*SEC, -BETA)                                    # the scalar wrapper at the edge: the limit, and inf map factor
    assert close(p["speed"], lim, 1e-13) and p["map_factor_abs"] == float("inf") and abs(p["dwdzeta_abs"]) < 1e-12 * F.U
    assert close([p["x"], p["y"], p["zeta_re"], p["zeta_im"]], [2 * B, 0.0, B, 0.0], 1e-14, 1e-15) and close(p["cp"], 1 - (lim / US) ** 2)


def test_thin_airfoil_V1_plate_parabolic_arc_and_exact_coefficients_of_a_cubic_camber_line():  # V1 + V2 — N18 (ours, standard)
    for al in (0.0, 0.04, -0.1):
        r = ch14.thin_airfoil(lambda x: 0.0 * x, al)
        assert close(r["CL"], 2 * PI * al, 1e-13, 1e-15) and abs(r["alpha_L0"]) < 1e-15 and close(r["A0"], al, 1e-13, 1e-15)
    h = 0.03                                                              # parabolic arc y_c = 4 h x (1 - x): slope 4 h (1 - 2x)
    for al in (0.0, 0.05):
        r = ch14.thin_airfoil(lambda x: 4 * h * (1 - 2 * x), al)
        assert close(r["alpha_L0"], -2 * h, 1e-12) and close(r["alpha_L0"], ch14.zero_lift_angle(h), 1e-12)     # = -beta of the Zhukhovsky arc (N31)
        assert close(r["CL"], 2 * PI * (al + 2 * h), 1e-12) and close(r["A1"], 4 * h, 1e-12) and abs(r["A2"]) < 1e-14
        assert close(r["CM_c4"], -PI * h, 1e-12) and close(r["A0"], al, 1e-12, 1e-14)
        assert close(r["CL"], ch14.zhukhovsky_lift_coefficient(al, 2 * h), 1e-12)                               # (14.12) with beta = 2 camber/c
    # exact Fourier coefficients of a cubic camber line with sympy (the theory's own integrals, not the code's quadrature)
    x, th = sp.symbols("x theta", real=True)
    yc = sp.Rational(1, 10) * x * (1 - x) * (sp.Rational(3, 5) - x)       # a reflexed (S-shaped) line of ours
    sl = sp.diff(yc, x).subs(x, (1 - sp.cos(th)) / 2)
    al = sp.Rational(1, 25)
    A0 = al - sp.integrate(sl, (th, 0, sp.pi)) / sp.pi
    A1 = 2 / sp.pi * sp.integrate(sl * sp.cos(th), (th, 0, sp.pi))
    A2 = 2 / sp.pi * sp.integrate(sl * sp.cos(2 * th), (th, 0, sp.pi))
    f = sp.lambdify(x, sp.diff(yc, x), "numpy")
    r = ch14.thin_airfoil(f, float(al))
    assert close([r["A0"], r["A1"], r["A2"]], [float(A0), float(A1), float(A2)], 1e-12, 1e-14)
    assert close(r["CL"], float(2 * sp.pi * (A0 + A1 / 2)), 1e-12) and close(r["CM_c4"], float(sp.pi / 4 * (A2 - A1)), 1e-12)
    aL0 = -sp.integrate(sl * (sp.cos(th) - 1), (th, 0, sp.pi)) / sp.pi
    assert close(r["alpha_L0"], float(aL0), 1e-12) and close(r["CL"], 2 * PI * (float(al) - r["alpha_L0"]), 1e-12)   # C_L = 2 pi (alpha - alpha_L0)
    assert close(ch14.thin_airfoil(f, 0.2)["CM_c4"], r["CM_c4"], 1e-12)   # the quarter chord is the aerodynamic centre: C_M there ignores alpha
    assert close(ch14.thin_airfoil(f, float(al), n=20)["CL"], r["CL"], 1e-12)                 # Gauss–Chebyshev: exact for polynomial slopes


def test_thin_airfoil_V3_discrete_vortex_method_gives_the_same_lift_and_moment():  # V3 — N18 by an independent numerical method
    # N lumped vortices on the chord (quarter-panel points), no flow through the camber line at the three-quarter points:
    # sum_j Gamma_j / (2 pi (x_i - x_j)) = U (alpha - dy_c/dx_i);  C_L = 2 sum Gamma / (U c);  C_M,c/4 = -2 sum Gamma (x - c/4) / (U c^2)
    slope = lambda x: 0.1 * ((1 - x) * (0.6 - x) - x * (0.6 - x) - x * (1 - x))
    al, ref = 0.04, None
    ref = ch14.thin_airfoil(slope, al)
    Ns, eL, eM = [25, 50, 100, 200, 400], [], []
    for N in Ns:
        xv = (np.arange(N) + 0.25) / N
        xc = (np.arange(N) + 0.75) / N
        M = 1.0 / (2 * PI * (xc[:, None] - xv[None, :]))
        G = np.linalg.solve(M, al - slope(xc))
        eL.append(abs(2 * G.sum() - ref["CL"]))
        eM.append(abs(-2 * np.sum(G * (xv - 0.25)) - ref["CM_c4"]))
    assert eL[-1] < 2e-4 * abs(ref["CL"]) and eM[-1] < 2e-3 * abs(ref["CM_c4"]), (eL, eM, ref)
    assert observed_order([1 / N for N in Ns], eL) > 0.9 and observed_order([1 / N for N in Ns], eM) > 0.9
    xv, xc = (np.arange(40) + 0.25) / 40, (np.arange(40) + 0.75) / 40                         # flat plate: the discrete method is exact
    G = np.linalg.solve(1.0 / (2 * PI * (xc[:, None] - xv[None, :])), np.full(40, 0.07))
    assert close(2 * G.sum(), ch14.thin_airfoil(lambda x: 0.0 * x, 0.07)["CL"], 1e-12)


# =====================================================================================================================
# 8. C07 — surface speed and pressure; N09, N10, N12, N13, R03
# =====================================================================================================================
def test_surface_pressure_V1_stagnation_points_bernoulli_and_the_edge_value():  # V1 — C07, R03
    n = 2000
    for e, beta in SECTIONS[1:]:
        for al in (0.1, -0.06):
            S = ch14.zhukhovsky_surface(2.0, 0.5, e, beta, al, n)
            F = ch14.zhukhovsky_flow(2.0, 0.5, e, beta, al)
            assert S["is_kutta"] and close(S["Gamma_cw"], F.Gamma_kutta, 1e-14, 1e-300) and S["z"].size == n and S["t"][0] == -beta
            assert close(S["u_theta"], ch14.cylinder_surface_speed(S["t"] - al, 2.0, F.a, Gamma_cw=F.Gamma_kutta), 1e-12, 1e-13)   # (14.10) = (6.37)
            q = np.abs(np.asarray(F.dwdz_from_zeta(S["zeta"])))                               # the complex chain rule, node by node
            assert close(S["speed"][1:], q[1:], 1e-9, 1e-12) and close(S["speed"], np.abs(S["u_s"]), 1e-15, 0)
            assert close(S["cp"], 1 - (S["speed"] / 2.0) ** 2, 1e-14, 1e-15) and close(S["cp"], ch14.pressure_coefficient(S["speed"], 2.0), 1e-13, 1e-14)
            assert S["cp"].max() <= 1.0 + 1e-14 and close(S["cp_min"], S["cp"].min()) and S["cp_min"] < 0
            lim = ch14.trailing_edge_speed(2.0, F.a, 0.5, al, beta)
            assert close(S["speed"][0], lim, 1e-14) and S["u_s"][0] < 0 and close(S["speed"][1], lim, 0.05)   # finite and continuous at the cusp
            assert len(S["stagnation"]) == 1                              # the rear stagnation point of the circle maps to the cusp: not a stagnation point
            zf = S["stagnation"][0]
            tf = float(np.angle(ch14.stagnation_points(2.0, 0.5, e, beta, al)["zeta_front"] - F.mu))
            pf = ch14.surface_point(2.0, 0.5, e, beta, al, tf)
            assert abs(pf["speed"]) < 1e-10 and abs(pf["cp"] - 1.0) < 1e-12 and abs(complex(pf["x"], pf["y"]) - zf) < 1e-12   # C_p = 1 exactly
            assert close(S["ds_dt"], F.a * np.abs(np.asarray(ch14.joukowski_derivative(S["zeta"], 0.5))), 1e-12, 1e-14)
            assert close(S["s"][-1] + abs(S["z"][0] - S["z"][-1]), np.sum(np.abs(np.diff(np.append(S["z"], S["z"][0])))), 1e-12)
            for side in ("upper", "lower"):
                d = S[side]
                assert np.all(np.diff(d["xc"]) > 0) and d["xc"][0] == 0.0 and d["xc"][-1] == 1.0 and d["index"][-1] == 0
                assert close(d["cp"], S["cp"][d["index"]]) and close(d["x"], S["x"][d["index"]])
            assert S["upper"]["y"].mean() > S["lower"]["y"].mean() and close(S["chord_x"], S["x"].max() - S["x"].min())
            k = 311                                                       # the scalar wrapper at a node
            p = ch14.surface_point(2.0, 0.5, e, beta, al, float(S["t"][k]))
            assert close([p["speed"], p["cp"], p["x"], p["y"]], [S["speed"][k], S["cp"][k], S["x"][k], S["y"][k]], 1e-9, 1e-12)
            assert close(p["speed"], p["dwdzeta_abs"] * p["map_factor_abs"], 1e-12) and close(p["dwdzeta_abs"], abs(S["u_theta"][k]), 1e-10)
            assert close(complex(p["zeta_re"], p["zeta_im"]), S["zeta"][k], 1e-14)
    # without the Kutta circulation: two stagnation points with C_p = 1, an infinite speed at the edge
    F = ch14.zhukhovsky_flow(*SEC)
    S0 = ch14.zhukhovsky_surface(*SEC, n, Gamma_cw=0.5 * F.Gamma_kutta)
    assert not S0["is_kutta"] and np.isinf(S0["speed"][0]) and len(S0["stagnation"]) == 2 and S0["cp"][0] == -np.inf
    for zs in S0["stagnation"]:
        assert abs(ch14.zhukhovsky_flow(*SEC, Gamma_cw=0.5 * F.Gamma_kutta).speed(zs * (1 + 1e-9))) < 1e-3 * US
    # suction on the upper surface near the nose at positive incidence; mirror image of a symmetric section
    S = ch14.zhukhovsky_surface(*SEC, n)
    assert S["xc_cp_min"] < 0.05 and S["y"][np.argmin(S["cp"])] > S["y"][np.argmin(S["x"])] - 1e-3
    assert np.trapezoid(S["lower"]["cp"], S["lower"]["xc"]) > np.trapezoid(S["upper"]["cp"], S["upper"]["xc"])           # more suction above
    Sp, Sm = ch14.zhukhovsky_surface(US, B, 0.1, 0.0, 0.12, 400), ch14.zhukhovsky_surface(US, B, 0.1, 0.0, -0.12, 400)
    assert close(Sp["cp"][1:], Sm["cp"][1:][::-1], 1e-9, 1e-10)
    Sz = ch14.zhukhovsky_surface(US, B, 0.1, 0.0, 0.0, 400)               # symmetric at zero incidence: identical curves, no loop
    assert abs(ch14.normal_force_coefficient(Sz["upper"]["xc"], Sz["upper"]["cp"], Sz["lower"]["xc"], Sz["lower"]["cp"])) < 1e-12


def test_surface_pressure_V4_loop_area_is_the_normal_force_and_the_speed_integral_is_the_circulation():  # V4 + V3 — C07, N09, N10, D10
    rho = RHO
    fb = ch14.zhukhovsky_force(*SEC, rho)
    q = 0.5 * rho * US ** 2
    ns, eC, eG = [500, 1000, 2000, 4000], [], []
    for n in ns:
        S = ch14.zhukhovsky_surface(*SEC, n)
        Cy = ch14.normal_force_coefficient(S["upper"]["xc"], S["upper"]["cp"], S["lower"]["xc"], S["lower"]["cp"])
        eC.append(abs(Cy - fb["Fy"] / (q * S["chord_x"])))                # the loop area against the Blasius F_y
        eG.append(abs(ch14.surface_circulation(S) - S["Gamma_cw"]) / S["Gamma_cw"])
        assert abs(ch14.surface_circulation(S, "circle") - S["Gamma_cw"]) < 1e-12 * S["Gamma_cw"]        # exact in the zeta-plane
    h = [1.0 / n for n in ns]
    assert abs(observed_order(h, eC) - 2.0) < 0.15 and abs(observed_order(h, eG) - 2.0) < 0.15, (eC, eG)
    assert eC[-1] < 1e-4 and eG[2] < 1e-6 and eG[-1] < 2.5e-7             # analysis T6: 1e-4 vs the pressure integral; docstring: 1e-6 at n = 2000
    CL_ideal = ch14.lift_from_cp(S["upper"]["xc"], S["upper"]["cp"], S["lower"]["xc"], S["lower"]["cp"], ALPHA)
    CL_none = ch14.lift_from_cp(S["upper"]["xc"], S["upper"]["cp"], S["lower"]["xc"], S["lower"]["cp"], ALPHA, suction="none")
    assert close(CL_ideal, Cy / math.cos(ALPHA)) and close(CL_none, Cy * math.cos(ALPHA))
    assert close(CL_ideal, 2 * S["Gamma_cw"] / (US * S["chord_x"]), 1e-5)                                # = 2 Gamma / (U c): the lift (T5)
    assert close(CL_ideal - CL_none, CL_ideal * math.sin(ALPHA) ** 2, 1e-12)                             # they differ by O(alpha^2)
    assert round(Cy, 3) == 1.054 and round(CL_ideal, 3) == 1.058                                         # design D10's worked number
    # with any other circulation the circle route still returns it (the polygon route needs a finite edge speed)
    Sg = ch14.zhukhovsky_surface(*SEC, 800, Gamma_cw=3.3)
    assert close(ch14.surface_circulation(Sg, "circle"), 3.3, 1e-12)
    with pytest.raises(ValueError):
        ch14.surface_circulation(Sg, "simpson")
    # N10: the circulation as the difference of the two surface-speed integrals, Gamma = int (U_upper - U_lower) ds
    S = ch14.zhukhovsky_surface(*SEC, 4000)
    il = int(np.argmin(np.abs(S["u_s"][1:]))) + 1                         # the front stagnation point splits the two layers
    seg = np.abs(np.diff(np.append(S["z"], S["z"][0])))
    mid = 0.5 * (S["speed"] + np.roll(S["speed"], -1))
    upper, lower = np.sum((mid * seg)[:il]), np.sum((mid * seg)[il:])
    assert close(upper - lower, S["Gamma_cw"], 2e-5)
    # exactness for piecewise-linear curves (a toy loop of ours)
    xs = np.array([0.0, 0.25, 1.0])
    assert close(ch14.normal_force_coefficient(xs, [-0.6, -0.6, -0.6], xs, [0.4, 0.4, 0.4]), 1.0, 1e-15)
    assert close(ch14.normal_force_coefficient(xs, [-1.0, -0.4, 0.0], xs, [0.5, 0.2, 0.0]), (0.25 * 0.35 + 0.75 * 0.1) + (0.25 * 0.7 + 0.75 * 0.2), 1e-14)
    with pytest.raises(ValueError):
        ch14.lift_from_cp(xs, xs, xs, xs, 0.1, suction="half")


def test_pressure_loop_V2_derivation():  # V2 — D10 (two-star): F_y = -int p_u dx + int p_l dx = closed-loop integral; C_y = loop area
    t, x = sp.symbols("t x", real=True)
    kk, p0, p1, p2, p3, pinf = sp.symbols("k p_0 p_1 p_2 p_3 p_inf", positive=True)
    X, Y = sp.cos(t), kk * sp.sin(t)                                      # a closed contour (ellipse), counterclockwise
    p = lambda xx, yy: p0 + p1 * xx + p2 * yy + p3 * xx * yy ** 2         # any smooth pressure field
    # force on the body: dF = -p n ds with n ds = (dy, -dx)  =>  F_y = closed integral of p dx
    Fy_loop = sp.integrate(p(X, Y) * sp.diff(X, t), (t, 0, 2 * sp.pi))
    yu, yl = kk * sp.sqrt(1 - x ** 2), -kk * sp.sqrt(1 - x ** 2)
    Fy_two = -sp.integrate(p(x, yu), (x, -1, 1)) + sp.integrate(p(x, yl), (x, -1, 1))                    # the book's two integrals
    assert sp.simplify(Fy_loop - Fy_two) == 0
    # Gauss: the same force is minus the volume integral of dp/dy (an independent route)
    r = sp.symbols("r", positive=True)
    dpdy = sp.diff(p(x, sp.Symbol("yy")), sp.Symbol("yy"))
    gauss = -sp.integrate(sp.integrate(dpdy.subs({x: r * sp.cos(t), sp.Symbol("yy"): kk * r * sp.sin(t)}) * kk * r, (r, 0, 1)), (t, 0, 2 * sp.pi))
    assert sp.simplify(Fy_loop - gauss) == 0
    assert sp.integrate(pinf * sp.diff(X, t), (t, 0, 2 * sp.pi)) == 0     # a uniform pressure pushes on nothing: p - p_inf may be used
    Fx_loop = -sp.integrate(p(X, Y) * sp.diff(Y, t), (t, 0, 2 * sp.pi))   # and F_x = -closed integral of p dy
    gx = -sp.integrate(sp.integrate(sp.diff(p(x, sp.Symbol("yy")), x).subs({x: r * sp.cos(t), sp.Symbol("yy"): kk * r * sp.sin(t)}) * kk * r, (r, 0, 1)), (t, 0, 2 * sp.pi))
    assert sp.simplify(Fx_loop - gx) == 0
    rho, U, c = sp.symbols("rho U c", positive=True)
    cu, cl = sp.symbols("C_pu C_pl", real=True)                            # constant C_p on each surface over 0 < x < c
    Fy = -sp.integrate(pinf + cu * rho * U ** 2 / 2, (x, 0, c)) + sp.integrate(pinf + cl * rho * U ** 2 / 2, (x, 0, c))
    assert sp.simplify(Fy / (rho * U ** 2 * c / 2) - (cl - cu)) == 0       # C_y = int (C_pl - C_pu) d(x/c)
    al, L = sp.symbols("alpha L", positive=True)                           # the lift is normal to the stream: F_y = L cos(alpha), F_x = -L sin(alpha)
    Lr, Dr = ch14.resolve_lift_drag(-math.sin(0.1) * 2.0, math.cos(0.1) * 2.0, 0.1)
    assert close(Lr, 2.0) and abs(Dr) < 1e-15 and sp.simplify(L * sp.cos(al) / sp.cos(al) - L) == 0


def test_flat_plate_V1_leading_edge_suction_closes_the_force_triangle():  # V1 + V4 — N12 (ours)
    rho, U, c = 1.19, 31.0, 0.9
    for al in (0.02, 0.1, 0.3):
        f = ch14.flat_plate_forces(rho, U, c, al)
        assert close(f["L"], PI * rho * U ** 2 * c * math.sin(al)) and close(f["CL"], 2 * PI * math.sin(al)) and f["D"] == 0.0
        assert close(f["Gamma_cw"], ch14.kutta_circulation(U, c / 4, al, 0.0)) and close(f["L"], ch14.lift_per_span(rho, U, Gamma_cw=f["Gamma_cw"]))
        assert close(f["N"] * math.cos(al) + f["S"] * math.sin(al), f["L"]) and abs(f["N"] * math.sin(al) - f["S"] * math.cos(al)) < 1e-12 * f["L"]
        assert close(f["S"], PI * rho * U ** 2 * c * math.sin(al) ** 2) and close(f["N"] ** 2 + f["S"] ** 2, f["L"] ** 2)
        fb = ch14.zhukhovsky_force(U, c / 4, 0.0, 0.0, al, rho)           # the Blasius force: chordwise part = -S, normal part = N
        assert close(fb["Fy"], f["N"], 1e-9) and close(fb["Fx"], -f["S"], 1e-9)
        # the edge flow of the plate: |dw/dz|^2 r -> U_o^2 a = U^2 c sin^2(alpha) at the leading edge, in every direction
        F = ch14.flat_plate_flow(U, c / 4, al)
        for ang in (PI, 2.0, -2.5):
            r = np.array([1e-8, 1e-10]) * c                               # |dw/dz|^2 r = K (1 + C sqrt(r) + ...): remove the sqrt(r) term
            k2 = np.asarray(F.speed(-c / 2 + r * np.exp(1j * ang))) ** 2 * r
            K = (k2[1] * math.sqrt(r[0]) - k2[0] * math.sqrt(r[1])) / (math.sqrt(r[0]) - math.sqrt(r[1]))
            assert close(K, U ** 2 * c * math.sin(al) ** 2, 1e-4), (al, ang, K)
        s = ch14.edge_flow_suction(rho, U * math.sin(al), c)              # control volume round the edge with U_o^2 a = U^2 c sin^2(alpha)
        assert close(s["closed_form"], f["S"]) and close(-s["suction"][0], f["S"], 1e-12)
    for eps in (1e-5, 1e-3, 0.3, 5.0):                                    # V4: the same force on every control circle
        s = ch14.edge_flow_suction(1.2, 3.0, 0.4, eps=eps, n=400)
        assert close(s["F_hold"][0], PI * 1.2 * 9.0 * 0.4, 1e-12) and abs(s["F_hold"][1]) < 1e-12 * s["closed_form"]
        assert close(s["momentum_flux"][0], s["closed_form"], 1e-12) and abs(s["pressure_force"][0]) < 1e-12 * s["closed_form"]
        assert s["suction"][0] < 0 and close(s["suction"][0], -s["F_hold"][0])                    # the fluid pulls the plate out of its own edge
    # our own quadrature of the same control volume (momentum flux and pressure), phi = 2 U_o sqrt(a r) cos(theta/2)
    Uo, a, eps = 3.0, 0.4, 0.07
    ur = lambda th: Uo * math.sqrt(a / eps) * math.cos(th / 2)
    ut = lambda th: -Uo * math.sqrt(a / eps) * math.sin(th / 2)
    mx = quad(lambda th: 1.2 * (ur(th) * math.cos(th) - ut(th) * math.sin(th)) * ur(th) * eps, 0, 2 * PI)[0]
    px = quad(lambda th: 0.5 * 1.2 * (ur(th) ** 2 + ut(th) ** 2) * math.cos(th) * eps, 0, 2 * PI)[0]
    assert close(mx, PI * 1.2 * Uo ** 2 * a, 1e-10) and abs(px) < 1e-10 * mx


def test_suction_side_separation_V7_qualitative_trends_only():  # V7 — N13 (label: qualitative; trap T16)
    al = np.deg2rad([0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
    out = [ch14.suction_side_separation(float(a), E, BETA, 1e6) for a in al]
    assert all(o["label"] == "qualitative" and o["separated"] and o["side"] == "upper" for o in out)
    xs = np.array([o["x_sep_over_c"] for o in out])
    assert np.all(np.diff(xs) < 0) and 0 < xs[-1] < 0.02 and 0.3 < xs[0] < 0.7        # the point moves forward as alpha grows
    assert all(o["xc_peak"] < o["x_sep_over_c"] for o in out)                         # separation lies behind the suction peak
    assert np.all(np.diff([o["xc_peak"] for o in out]) < 0)                           # the peak moves to the nose
    up, lo = ch14.suction_side_separation(0.0, 0.1, 0.0, 1e6, side="upper"), ch14.suction_side_separation(0.0, 0.1, 0.0, 1e6, side="lower")
    assert close(up["x_sep_over_c"], lo["x_sep_over_c"], 1e-9) and lo["side"] == "lower"            # symmetric section, zero incidence
    a_hi, a_lo = ch14.suction_side_separation(0.07, E, BETA, 1e7), ch14.suction_side_separation(0.07, E, BETA, 1e5)
    assert close(a_hi["x_sep_over_c"], a_lo["x_sep_over_c"], 1e-9)                    # Thwaites: the position ignores Re ...
    assert close(a_lo["theta_sep"] / a_hi["theta_sep"], 10.0, 1e-6)                   # ... the momentum thickness goes as Re^(-1/2)
    pr = ch14.suction_side_separation(0.1, E, BETA, 1e6, side="pressure")
    assert pr["side"] == "lower" and (not pr["separated"] or pr["x_sep_over_c"] > out[3]["x_sep_over_c"])
    assert ch14.suction_side_separation(-0.2, E, BETA, 1e6)["side"] == "lower"        # negative lift: the suction side is underneath
    assert out[2]["Ue"][0] == 0.0 and close(out[2]["Ue"][-1], ch14.trailing_edge_speed(1.0, ch14.zhukhovsky_circle(1.0, E, BETA)[1], 1.0, al[2], BETA))
    with pytest.raises(ValueError):
        ch14.suction_side_separation(0.1, 0.0, 0.0, 1e6)


# =====================================================================================================================
# 9. C08 — bound vortex and trailing filaments; N36, N37, N39, N40, N49, N50, R11
# =====================================================================================================================
def _gamma_fn(Gn, s):
    return lambda y: ch14.circulation_series(Gn, ch14.span_angle(y, s))


def test_shed_filament_V1_sheet_strength_is_minus_dGamma_dy_and_V4_helmholtz_totals():  # V1 + V4 — C08, D11, N40, N49, N50
    s = SPAN
    ga = np.linspace(0.0, PI, 181)
    y = ch14.span_position(ga, s)
    assert y[0] == -s / 2 and close(y[-1], s / 2) and close(ch14.span_angle(y, s), ga, 1e-12, 1e-7)       # gamma = 0 port tip, pi starboard tip
    assert close(y, -(s / 2) * np.cos(ga))
    for Gn in (GN, GN_ASYM, [1.7]):
        G = ch14.circulation_series(Gn, ga)
        assert abs(G[0]) < 1e-14 and abs(G[-1]) < 1e-14                   # every term vanishes at both tips: Gamma(+-s/2) = 0
        assert close(G, sum(g * np.sin((k + 1) * ga) for k, g in enumerate(Gn)), 1e-14, 1e-15)                # (14.18)
        gi = ga[20:-20]
        fn, h = _gamma_fn(Gn, s), 1e-5
        yi = ch14.span_position(gi, s)
        fd = -(fn(yi + h) - fn(yi - h)) / (2 * h)                         # -dGamma/dy by centred differences
        assert close(ch14.trailing_sheet_strength(Gn, gi, s), fd, 1e-7, 1e-8)
        # Helmholtz: what leaves the starboard half is the mid-span circulation; the whole sheet carries no net circulation
        dG = lambda g: float(np.sum(np.arange(1, len(Gn) + 1) * np.asarray(Gn) * np.cos(np.arange(1, len(Gn) + 1) * g)))   # dGamma/dgamma
        star = quad(lambda g: -dG(g), PI / 2, PI, epsabs=1e-13)[0]
        port = quad(lambda g: -dG(g), 0.0, PI / 2, epsabs=1e-13)[0]
        assert close(star, ch14.circulation_series(Gn, PI / 2), 1e-12) and abs(star + port) < 1e-12
        tot = quad(lambda g: float(ch14.trailing_sheet_strength(Gn, g, s)) * (s / 2) * math.sin(g), 0, PI, epsabs=1e-12)[0]
        assert abs(tot) < 1e-10                                           # the coded strength integrates to zero over the span
    assert close(ch14.circulation_series([1.7], ga), ch14.elliptic_circulation(y, s, 1.7), 1e-12, 1e-8)       # one term: the ellipse (14.23)
    assert close(ch14.elliptic_circulation(0.0, s, 1.7), 1.7) and ch14.elliptic_circulation(s / 2, s, 1.7) == 0.0
    for n in (1, 2, 3, 4):                                                # odd n symmetric about mid-span, even n antisymmetric
        e_n = np.zeros(n)
        e_n[-1] = 1.0
        assert close(ch14.circulation_series(e_n, PI - ga), (-1) ** (n + 1) * ch14.circulation_series(e_n, ga), 1e-12, 1e-13)
    st = ch14.trailing_sheet_strength(GN, ga[1:-1], s)                     # symmetric loading: the two half-sheets rotate in opposite senses
    assert close(st, -st[::-1], 1e-9, 1e-9) and st[-5] > 0 > st[4]
    assert np.isinf(ch14.trailing_sheet_strength([1.0], 0.0, s)) or abs(ch14.trailing_sheet_strength([1.0], 1e-9, s)) > 1e6   # infinite at the tips
    # a staircase 5, 8, 10, 10, 8, 5 sheds 5, 3, 2 on each side: total 10 = the root value (design D11)
    edges, steps = np.arange(-3.0, 4.0), np.array([5.0, 8.0, 10.0, 10.0, 8.0, 5.0])
    w, parts = ch14.staircase_downwash(0.37, edges, steps, parts=True)
    jumps = parts * 4 * PI * (edges - 0.37)                               # Gamma_left - Gamma_right at each edge
    assert close(jumps, [-5, -3, -2, 0, 2, 3, 5], 1e-12, 1e-12) and abs(jumps.sum()) < 1e-12 and close(jumps[4:].sum(), 10.0) and close(w, parts.sum())
    ye, Gs = ch14.staircase_from_loading(_gamma_fn([2.0], s), s, 8)
    assert close(np.diff(ye), s / 8) and close(Gs, ch14.elliptic_circulation(0.5 * (ye[:-1] + ye[1:]), s, 2.0), 1e-12)
    yc, _ = ch14.staircase_from_loading(_gamma_fn([2.0], s), s, 8, spacing="cosine")
    assert close(ch14.span_angle(yc, s), np.linspace(0, PI, 9), 1e-9, 1e-7)
    with pytest.raises(ValueError):
        ch14.staircase_from_loading(_gamma_fn([2.0], s), s, 8, spacing="chebyshev")


def test_tip_vortices_V1_rolled_up_pair_spacing_is_the_conserved_impulse():  # V1 + V4 — N37 (ours)
    s = SPAN
    p = ch14.rolled_up_pair([2.0], s)
    assert close(p["Gamma0"], 2.0) and close(p["spacing"], PI * s / 4) and close(p["descent_speed"], 2.0 / (2 * PI * (PI * s / 4)))
    for Gn in (GN, [2.0, 0.0, -0.3]):
        p = ch14.rolled_up_pair(Gn, s)
        fn = _gamma_fn(Gn, s)
        G0 = float(fn(0.0))
        area = quad(lambda y: float(fn(y)), -s / 2, s / 2, epsabs=1e-12)[0]
        assert close(p["Gamma0"], G0, 1e-13) and close(p["spacing"], area / G0, 1e-9)                         # b' = (1/Gamma_0) int Gamma dy
        n = np.arange(1, len(Gn) + 1)
        cen = quad(lambda g: -(s / 2) * math.cos(g) * -float(np.sum(n * np.asarray(Gn) * np.cos(n * g))), PI / 2, PI, epsabs=1e-13)[0] / G0
        assert close(2 * cen, p["spacing"], 1e-10)                        # = twice the centroid of the starboard half-sheet
        # a pair of line vortices +-Gamma_0 that far apart sinks at Gamma_0 / (2 pi b')  (chapter 5)
        u = ch14.point_vortex_velocity(np.array([[p["spacing"] / 2], [0.0]]), np.array([[-p["spacing"] / 2], [0.0]]), [-G0])
        assert close(-float(np.asarray(u)[1, 0]), p["descent_speed"], 1e-12)
    W, rho, U = I["W"], I["rho"], I["U_cruise"]
    t = ch14.tip_vortex_estimate(W, rho, U, I["s"])
    assert close(ch14.lift_from_series([t["Gamma0"]], I["s"], rho, U), W, 1e-13)                               # L = W through (14.21)
    assert close(t["spacing"], PI * I["s"] / 4) and close(t["descent_speed"], t["Gamma0"] / (2 * PI * t["spacing"]))
    assert close(t["downwash"], ch14.elliptic_downwash(t["Gamma0"], I["s"])) and close(t["induced_drag"], ch14.elliptic_induced_drag(W, rho, U, I["s"]))
    assert close(t["induced_drag"], rho * t["downwash"] * t["Gamma0"] * PI * I["s"] / 4, 1e-13)                # D_i = rho w int Gamma dy
    pr = ch14.rolled_up_pair([t["Gamma0"]], I["s"])
    assert close([pr["spacing"], pr["descent_speed"]], [t["spacing"], t["descent_speed"]])


def test_wake_rollup_V4_invariants_are_conserved_the_spiral_is_qualitative():  # V4 — N37 (label: conserved invariants; qualitative spiral)
    Gn, s = [2.0, 0.0, 0.2], SPAN
    r = ch14.wake_rollup(Gn, s, fast=True)
    assert "qualitative" in r["label"] and r["y"].shape == r["z"].shape == (41, 80) and r["sigma"].shape == (80,)
    assert abs(r["sigma"].sum()) < 1e-13 and close(r["half_circulation"], r["Gamma0"], 1e-13) and close(r["Gamma0"], 1.8)    # each half sheds Gamma_0
    assert np.ptp(r["impulse"]) < 1e-9 * abs(r["impulse"][0])                                                  # vertical impulse sum sigma y
    assert close(r["impulse"][0], -np.sum(r["sigma"] * -r["y"][0]), 1e-13)
    assert np.ptp(r["energy"]) < 1e-5 * abs(r["energy"][0])                                                    # the blob Hamiltonian
    assert np.ptp(r["centroid_y"]) < 1e-8 * s and close(r["centroid_y"][0], r["spacing_theory"] / 2, 2e-3)     # centroids stay b'/2 from the centre
    assert close(r["spacing_theory"], ch14.rolled_up_pair(Gn, s)["spacing"]) and close(r["descent_speed_theory"], ch14.rolled_up_pair(Gn, s)["descent_speed"])
    assert np.all(np.diff(r["centroid_z"]) < 0) and r["centroid_z"][0] == 0.0                                  # the wake sinks
    assert np.max(np.abs(r["y"] + r["y"][:, ::-1])) < 1e-6 * s and np.max(np.abs(r["z"] - r["z"][:, ::-1])) < 1e-6 * s   # mirror symmetry kept
    X = np.stack([r["y"], r["z"]], axis=1)                                 # (T, 2, M): the ch05 invariants of point vortices
    inv = ch14.point_vortex_invariants(X, r["sigma"])
    assert np.ptp(inv["P"][:, 0]) < 1e-9 * abs(inv["P"][0, 0]) and np.max(np.abs(inv["P"][:, 1])) < 1e-9 * abs(inv["P"][0, 0])
    assert close(r["t"][-1], 0.6 * s ** 2 / 1.8) and r["eps"] == 0.05 * s
    # the roll-up itself: the sheet's edge curls up and inward — the tip vortex sits inboard of the tip, the sheet no longer flat
    assert np.ptp(r["z"][-1]) > 0.05 * s and np.ptp(r["z"][0]) == 0.0
    # the centroid of each half sinks at the speed of the rolled-up pair, Gamma_0 / (2 pi b'), from the very start — five
    # times slower than the flat sheet's own 2w = Gamma_1/s (the tip singularity carries the weight of the average)
    mean_sink = -r["centroid_z"][-1] / r["t"][-1]
    assert close(mean_sink, r["descent_speed_theory"], 0.02) and mean_sink < 0.3 * 2.0 / s
    sink = -np.gradient(r["centroid_z"], r["t"])
    assert np.all(np.abs(sink / r["descent_speed_theory"] - 1) < 0.06)


# =====================================================================================================================
# 10. C09 — the downwash (14.13); N41, N60, N61, N65, N66, R12
# =====================================================================================================================
def test_downwash_V1_quadrature_of_14_13_against_closed_forms_and_the_series():  # V1 — C09 (14.13), N56 (14.24), N61
    s = SPAN
    ys = np.array([-3.3, -2.0, -0.6, 0.0, 0.5, 1.9, 3.2, 3.45])
    # (a) parabolic loading Gamma_0 (1 - (2y/s)^2): our own closed form, w = (2 Gamma_0 / (pi s^2)) [s + y1 ln|(s/2 - y1)/(s/2 + y1)|]
    G0 = 3.0
    w = ch14.downwash(ys, lambda y: G0 * (1 - (2 * y / s) ** 2), s)
    exact = 2 * G0 / (PI * s ** 2) * (s + ys * np.log(np.abs((s / 2 - ys) / (s / 2 + ys))))
    assert close(w, exact, 1e-9, 1e-11) and w[3] > 0 and close(w[3], 2 * G0 / (PI * s))                        # downward behind a lifting wing
    assert w[-1] < 0                                                      # ... and an upwash right at the tips of this loading
    # (b) elliptic loading: the same downwash Gamma_1 / 2s at every station (14.24)
    we = ch14.downwash(ys, lambda y: ch14.elliptic_circulation(y, s, 2.0), s)
    assert close(we, 2.0 / (2 * s), 1e-9) and close(ch14.elliptic_downwash(2.0, s), 2.0 / (2 * s)) and type(ch14.downwash(0.3, _gamma_fn([2.0], s), s)) is float
    assert close(w[3] / we[3], (2 * G0 / (PI * s)) / (2.0 / (2 * s)))     # N61: mid-span downwash of the parabola against the ellipse
    # (c) series loadings: quadrature of (14.13) = the last member of (14.19)
    for Gn in ([2.0], [2.0, 0.3], GN, GN_ASYM):
        wq = ch14.downwash(ys, _gamma_fn(Gn, s), s)
        assert close(wq, ch14.downwash_series(Gn, ch14.span_angle(ys, s), s), 1e-8, 1e-10), Gn
    # (d) the principal value by symmetric exclusion in y (no singularity subtraction, no span angle)
    dG = lambda y: -8 * G0 * y / s ** 2
    for y1 in (-2.0, 0.5, 3.0):
        pv = pv_symmetric(lambda y: dG(y) / (y1 - y), -s / 2, s / 2, y1) / (4 * PI)
        assert close(ch14.downwash(y1, lambda y: G0 * (1 - (2 * y / s) ** 2), s), pv, 1e-6, 2e-7)    # the excluded gap leaves 8e-8
    with pytest.raises(ValueError):
        ch14.downwash(0.0, _gamma_fn([2.0], s), s, method="panel")


def test_downwash_V3_a_staircase_of_horseshoes_converges_to_14_13():  # V3 — C09 from-scratch route (Fig. 14.21), C08
    s, G1 = SPAN, 2.0
    fn = lambda y: ch14.elliptic_circulation(y, s, G1)
    Ns, eu, ec, em = [20, 40, 80, 160, 320], [], [], []
    for N in Ns:
        ye, Gs = ch14.staircase_from_loading(fn, s, N)
        ym = 0.5 * (ye[:-1] + ye[1:])
        w = ch14.staircase_downwash(ym, ye, Gs)
        eu.append(abs(w[N // 2] - G1 / (2 * s)))                          # mid-span, uniform panels
        em.append(np.max(np.abs(w[N // 4: 3 * N // 4] - G1 / (2 * s))))   # inner half of the span
        yc, Gc = ch14.staircase_from_loading(fn, s, N, spacing="cosine")
        wc = ch14.staircase_downwash(0.5 * (yc[:-1] + yc[1:]), yc, Gc)
        ec.append(abs(wc[N // 2] - G1 / (2 * s)))                         # mid-span, cosine panels
    h = [1.0 / N for N in Ns]
    pu, pm, pc = observed_order(h, eu), observed_order(h, em), observed_order(h, ec)
    # measured orders for the elliptic loading (20-320 panels): 1.5 with uniform panels (1.517 at mid-span, 1.473 over the
    # inner half; successive pairs 1.54 -> 1.50 and 1.44 -> 1.49), 2.0 with cosine panels (2.023; pairs 2.05 -> 2.01)
    assert abs(pu - 1.5) < 0.05 and abs(pm - 1.5) < 0.05 and abs(pc - 2.0) < 0.05, (pu, pm, pc)
    assert eu[-1] < 1e-5 and ec[-1] < 5e-6 and ec[-1] < eu[-1]            # 3.9e-6 and 2.3e-6 m/s at 320 panels: cosine panels are the better ones
    # one panel = one horseshoe: two half-lines at the distance s/2 -> Gamma / (pi s) at the centre (design D12: 2 m/s for 4 pi, 2 m)
    assert close(ch14.staircase_downwash(0.0, [-1.0, 1.0], [4 * PI]), 2.0) and close(ch14.horseshoe_downwash(0.0, 0.0, 4 * PI, 2.0), 2.0, 1e-12)
    # the same sum built from Biot–Savart horseshoes (bound segments induce nothing on their own line)
    ye, Gs = ch14.staircase_from_loading(fn, s, 12)
    y1 = np.array([-2.1, 0.13, 1.6])
    hs = sum(np.asarray(ch14.horseshoe_downwash(0.0, y1 - 0.5 * (ye[k] + ye[k + 1]), float(Gs[k]), float(ye[k + 1] - ye[k]))) for k in range(12))
    assert close(hs, ch14.staircase_downwash(y1, ye, Gs), 1e-11)
    # R12: one half-line gives half the velocity of an infinite line vortex
    assert close(ch14.segment_speed(0.7, PI / 2, PI, 3.0), 3.0 / (4 * PI * 0.7)) and close(ch14.segment_speed(0.7, 0.0, PI, 3.0), 3.0 / (2 * PI * 0.7))
    assert close(ch14.segment_speed(0.7, 0.0, PI / 2, 3.0), 3.0 / (4 * PI * 0.7))     # the same leg described from its far end (design D12 step 1)


def _biot_savart_numeric(P, a, d, Gamma, length=np.inf):
    """Our own quadrature of (5.17) for a straight filament from a along the unit vector d (not the closed segment formula)."""
    P, a, d = np.asarray(P, float), np.asarray(a, float), np.asarray(d, float)

    def comp(k):
        def f(t):
            r = P - (a + t * d)
            return Gamma / (4 * PI) * np.cross(d, r)[k] / np.linalg.norm(r) ** 3
        return quad(f, 0.0, length, limit=400, epsabs=1e-13, epsrel=1e-12)[0]

    return np.array([comp(0), comp(1), comp(2)])


def test_horseshoe_V1_biot_savart_field_centre_far_wake_and_upwash_outboard():  # V1 + V4 — N60, N65, R12 (trap T10)
    G, s = 3.7, 2.6
    h = s / 2
    assert close(ch14.horseshoe_downwash(0.0, 0.0, G, s), G / (PI * s), 1e-13)                                 # centre of the bound vortex
    assert close(ch14.horseshoe_downwash(1e7, 0.0, G, s), 2 * G / (PI * s), 1e-6) and abs(ch14.horseshoe_downwash(-1e7, 0.0, G, s)) < 1e-6 * G / s
    y = np.array([-3 * s, -0.9 * s, -0.3 * s, 0.1 * s, 0.45 * s, 0.7 * s, 2 * s])
    plane = G * h / (2 * PI * (h ** 2 - y ** 2))                           # two half-lines, added by hand, on the lifting line
    assert close(ch14.horseshoe_downwash(0.0, y, G, s), plane, 1e-12)
    assert np.all(plane[np.abs(y) > h] < 0) and np.all(plane[np.abs(y) < h] > 0)                               # upwash outboard of the tips
    assert ch14.horseshoe_downwash(0.0, 0.4999 * s, G, s) > 1000 * G / (PI * s)                                # unbounded toward the tips
    # every component at generic points: our quadrature of (5.17) on the three filaments
    for P in ((0.7, 0.3, 0.4), (-1.1, 1.9, -0.6), (2.5, -0.2, 0.05)):
        ref = (_biot_savart_numeric(P, (0, -h, 0), (0, 1, 0), G, s) + _biot_savart_numeric(P, (0, h, 0), (1, 0, 0), G)
               - _biot_savart_numeric(P, (0, -h, 0), (1, 0, 0), G))
        got = np.array(ch14.horseshoe_velocity(*P, G, s))
        assert close(got, ref, 1e-8, 1e-10), (P, got, ref)
        assert close(ch14.horseshoe_downwash(P[0], P[1], G, s, z=P[2]), -got[2], 1e-14)                        # w is positive downward (T5)
        seg = (np.asarray(ch14.segment_induced_velocity(np.array(P)[:, None], np.array([0, -h, 0.0]), np.array([0, h, 0.0]), G), float).ravel()
               + np.asarray(ch14.segment_induced_velocity(np.array(P)[:, None], np.array([0, h, 0.0]), np.array([1e7, h, 0.0]), G), float).ravel()
               + np.asarray(ch14.segment_induced_velocity(np.array(P)[:, None], np.array([1e7, -h, 0.0]), np.array([0, -h, 0.0]), G), float).ravel())
        assert close(got, seg, 1e-6, 1e-9)                                # three closed segments with a very long leg
        assert close(np.array(ch14.horseshoe_velocity(*P, G, s, x_end=1e7)), got, 1e-6, 1e-9)
        d = 1e-4                                                          # away from the filaments the field has no divergence and no curl
        J = np.array([(np.array(ch14.horseshoe_velocity(*(np.array(P) + d * e), G, s)) - np.array(ch14.horseshoe_velocity(*(np.array(P) - d * e), G, s))) / (2 * d)
                      for e in np.eye(3)]).T                              # J[i, j] = du_i/dx_j
        scale = np.max(np.abs(J))
        assert abs(np.trace(J)) < 1e-6 * scale and np.max(np.abs(J - J.T)) < 1e-6 * scale
    ux, uy, uz = ch14.horseshoe_velocity(np.array([0.5, 1.0]), np.array([0.2, -0.4]), 0.3, G, s)
    assert ux.shape == (2,) and type(ch14.horseshoe_velocity(0.5, 0.2, 0.3, G, s)[2]) is float
    # the "simplest wing" (Exercise-type result, ours): centre downwash for the whole span gives half the elliptic drag
    rho, U, c = 1.2, 40.0, 0.5
    hw = ch14.horseshoe_wing(G, s, rho, U)
    assert close(hw["L"], rho * U * G * s) and close(hw["w0"], G / (PI * s)) and close(hw["Di_centre"], rho * G ** 2 / PI) and close(hw["ratio"], 0.5)
    q = 0.5 * rho * U ** 2 * s * c
    assert close(hw["Di_centre"] / q, (hw["L"] / q) ** 2 / (2 * PI * (s / c)))                                 # C_Di = C_L^2 / (2 pi Lambda)
    assert close(hw["Di_elliptic_same_lift"], ch14.elliptic_induced_drag(hw["L"], rho, U, s)) and "tips" in hw["note"]


def test_ground_effect_V1_image_horseshoe_cancels_the_normal_velocity_on_the_ground():  # V1 + V7 — N66 (ours)
    G, s = 2.9, 3.4
    for hgt in (0.1 * s, 0.25 * s, s, 5 * s):
        own = ch14.horseshoe_velocity(0.0, 0.0, 0.0, G, s)[2]
        img = ch14.horseshoe_velocity(0.0, 0.0, 2 * hgt, -G, s)[2]         # the wing's centre is 2h above the image horseshoe
        assert close(-(own + img) / (G / (PI * s)), ch14.ground_effect_factor(hgt, s), 1e-12)
        assert close(ch14.ground_effect_factor(hgt, s), 1 - 1 / (1 + (4 * hgt / s) ** 2))
        for x, y in ((0.4, 0.3), (-1.0, 2.5), (6.0, -0.2)):               # on the ground plane (z = -h) nothing goes through
            wall = ch14.horseshoe_velocity(x, y, -hgt, G, s)[2] + ch14.horseshoe_velocity(x, y, hgt, -G, s)[2]
            assert abs(wall) < 1e-13 * G / s
    assert ch14.ground_effect_factor(0.0, s) == 0.0 and abs(ch14.ground_effect_factor(1e6 * s, s) - 1) < 1e-12
    f = ch14.ground_effect_factor(np.linspace(0.01, 3, 40) * s, s)
    assert np.all(np.diff(f) > 0) and close(ch14.ground_effect_factor(s / 4, s), 0.5)                          # half the downwash at h = s/4


def test_downwash_along_x_V7_zero_far_ahead_w_at_the_line_twice_w_far_behind():  # V7 + V1 — N41 (Fig. 14.22)
    s = SPAN
    for Gn, y1 in (([2.0], 0.0), (GN, 0.0), (GN_ASYM, 1.1)):
        w0 = ch14.downwash_series(Gn, ch14.span_angle(y1, s), s)
        x = np.array([-300 * s, -3 * s, -0.5 * s, -0.05 * s, 0.0, 0.05 * s, 0.5 * s, 3 * s, 300 * s])
        w = ch14.downwash_along_x(x, Gn, s, y1=y1)
        assert abs(w[0]) < 2e-3 * abs(w0) and close(w[4], w0, 1e-13) and close(w[-1], 2 * w0, 2e-3)            # 0, w, 2w
        assert close(w + w[::-1], 2 * w0, 1e-9, 1e-11)                    # the sheet's field is odd about the value at the line
        if len(Gn) == 1:
            assert np.all(np.diff(w) > 0)                                 # elliptic loading: a monotonic rise (Fig. 14.22)
        # (a loading with a dip at mid-span has reversed filaments there and is not monotonic near the line: the slope
        #  dw/dx at x = 0 is -(Gamma''/4 pi) int dy/|y - y1| — logarithmically infinite, of the sign of -Gamma'')
    # against a direct Biot–Savart sum over a fine staircase of horseshoes (legs only, and legs + bound vortex)
    # (the station must be a panel mid-point: downstream it lies in the sheet, and only the symmetric sum is the principal value)
    Gn = GN
    ye, Gs = ch14.staircase_from_loading(_gamma_fn(Gn, s), s, 1400)
    yc, wd = 0.5 * (ye[:-1] + ye[1:]), np.diff(ye)
    y1 = float(yc[780])
    for x in (-0.8 * s, -0.2 * s, 0.3 * s, 1.5 * s):
        full = legs = 0.0
        tot = -np.asarray(ch14.horseshoe_velocity(x, y1 - yc, 0.0, 1.0, 1.0)[2]) * 0     # shape only
        for k in range(yc.size):
            uz = ch14.horseshoe_velocity(x, y1 - yc[k], 0.0, float(Gs[k]), float(wd[k]))[2]
            ub = np.asarray(ch14.segment_induced_velocity(np.array([[x], [y1 - yc[k]], [0.0]]), np.array([0, -wd[k] / 2, 0.0]),
                                                          np.array([0, wd[k] / 2, 0.0]), float(Gs[k])), float).ravel()[2]
            full, legs = full - uz, legs - (uz - ub)
        assert close(ch14.downwash_along_x(x, Gn, s, y1=y1, include_bound=True), full, 2e-3), x
        assert close(ch14.downwash_along_x(x, Gn, s, y1=y1), legs, 2e-3), x
        assert tot.shape == yc.shape
    b = ch14.downwash_along_x(np.array([-s, s]), Gn, s, include_bound=True) - ch14.downwash_along_x(np.array([-s, s]), Gn, s)
    assert b[0] < 0 < b[1] and close(b[0], -b[1], 1e-9)                   # the bound vortex: upwash ahead, downwash behind


# =====================================================================================================================
# 11. C10 — induced drag (14.14)–(14.15); N38, N42–N45, N67
# =====================================================================================================================
def test_induced_drag_V1_strip_forces_effective_angle_and_slip_2():  # V1 + V7 — C10, N42 (slip #2), N43 (14.14), N44, N45
    U, w, G, rho, dy = 50.0, np.array([0.0, 1.0, 7.0, 25.0]), 20.0, 1.2, 1.0
    for exact in (False, True):
        f = ch14.strip_forces(G, w, U, rho, dy, exact=exact)
        assert close(f["dL"], rho * U * G * dy) and close(f["dDi"], rho * w * G * dy, 1e-14, 1e-300) and close(f["Ue"], np.hypot(U, w))
        assert close(f["dL"] ** 2 + f["dDi"] ** 2, f["dLe"] ** 2, 1e-13) and close(f["dLe"], rho * f["Ue"] * G * dy)
        assert f["dDi"][0] == 0.0                                          # no downwash, no induced drag
    eps = ch14.strip_forces(G, w, U, rho, dy, exact=True)["eps"]
    tri = np.arctan2(w, U)                                                # the angle of the velocity triangle of Fig. 14.23
    assert close(eps, tri, 1e-14, 1e-300) and close(ch14.induced_angle(w, U, small=False), tri, 1e-14, 1e-300)
    f = ch14.strip_forces(G, w, U, rho, dy, exact=True)
    assert close(f["dLe"] * np.cos(eps), f["dL"], 1e-13) and close(f["dLe"] * np.sin(eps), f["dDi"], 1e-13, 1e-12)   # exact: U_e cos = U, U_e sin = w
    assert close(ch14.induced_angle(w, U), w / U, 1e-15, 0) and close(ch14.strip_forces(G, w, U, rho, dy)["eps"], w / U, 1e-15, 0)
    small = np.array([0.4, 0.2, 0.1, 0.05])
    assert abs(observed_order(small, np.abs(np.arctan(small) - small)) - 3.0) < 0.1                            # eps = w/U + O((w/U)^3)
    printed = np.tan(w[1:] / U)                                           # slip #2: "epsilon = tan(w/U)" as printed
    assert np.all(np.abs(printed - tri[1:]) > 100 * np.abs(tri[1:] - tri[1:])) and np.all(printed > tri[1:])
    assert close((printed - tri[1:]) / (w[1:] / U) ** 3, 2 / 3, 0.25) and abs(printed[-1] - tri[-1]) > 0.08    # wrong by (2/3)(w/U)^3: 5 deg at w/U = 1/2
    assert close(ch14.effective_angle(0.1, w, U), 0.1 - w / U) and ch14.effective_angle(0.1, 0.0, U) == 0.1    # (14.14)
    d = ch14.strip_forces(20.0, 1.0, 50.0, 1.2, 1.0)                       # design D13's worked number
    assert close([d["dL"], d["dDi"], d["eps"]], [1200.0, 24.0, 0.02])


def test_induced_drag_V1_quadrature_of_14_15_equals_14_21_and_14_22():  # V1 + V3 — C10 (14.15), N53 (14.21), C13 (14.22)
    s, rho, U = SPAN, 1.21, 37.0
    for Gn in (GN, GN_ASYM, [1.4]):
        Lx, Dx = ch14.lift_from_series(Gn, s, rho, U), ch14.induced_drag_from_series(Gn, rho)
        eL, eD, ns = [], [], [201, 401, 801, 1601]
        for n in ns:
            ga = np.linspace(0.0, PI, n)
            y, G, w = ch14.span_position(ga, s), ch14.circulation_series(Gn, ga), ch14.downwash_series(Gn, ga, s)
            L, Di = ch14.lift_and_induced_drag(y, G, w, rho, U)
            eL.append(abs(L - Lx) / abs(Lx))
            eD.append(abs(Di - Dx) / Dx)
        if len(Gn) > 1:
            assert abs(observed_order([1 / n for n in ns], eL) - 2.0) < 0.15 and abs(observed_order([1 / n for n in ns], eD) - 2.0) < 0.15
        assert eL[-1] < 1e-5 and eD[-1] < 1e-5
    # midpoint rule in the span angle (exact for these trigonometric polynomials) with the downwash from the direct
    # quadrature of (14.13): no series formula for w anywhere
    Gn = GN_ASYM
    ga = (np.arange(16) + 0.5) * PI / 16
    y, G = ch14.span_position(ga, s), ch14.circulation_series(Gn, ga)
    wq = ch14.downwash(y, _gamma_fn(Gn, s), s)
    jac = (PI / 16) * (s / 2) * np.sin(ga)
    assert close(rho * U * np.sum(G * jac), ch14.lift_from_series(Gn, s, rho, U), 1e-12)
    assert close(rho * np.sum(wq * G * jac), ch14.induced_drag_from_series(Gn, rho), 1e-8)
    assert close(ch14.lift_from_series([12.0, 0.0, 5.0], 12.0, 1.2, 50.0), 6786.0, 1e-4)                       # design D17: whatever Gamma_3 is
    assert round(ch14.induced_drag_from_series([12.0, 0.0, 1.0], 1.2), 1) == 69.3                              # design D18


def test_downwash_integral_V2_derivation():  # V2 — D12 (two-star): the half-line result (R12) summed over the sheet -> (14.13)
    G, Gf, d, s, G0 = sp.symbols("Gamma Gamma_f d s Gamma_0", positive=True)
    x, y, y1, dG = sp.symbols("x y y_1 dGamma_dy", real=True)
    # Biot–Savart (5.17) along a straight filament on the x-axis, seen from a point at the distance d abreast of its start
    kernel = Gf * d / (4 * sp.pi * (x ** 2 + d ** 2) ** sp.Rational(3, 2))
    half = sp.integrate(kernel, (x, 0, sp.oo))
    full = sp.integrate(kernel, (x, -sp.oo, sp.oo))
    assert sp.simplify(half - Gf / (4 * sp.pi * d)) == 0 and sp.simplify(full - 2 * half) == 0      # half of an infinite line vortex
    # direction: circulation Gf about +x at spanwise position y; the point (0, y1, 0) on the lifting line
    xh, r = sp.Matrix([1, 0, 0]), sp.Matrix([0, y1 - y, 0])
    u = Gf / (4 * sp.pi) * xh.cross(r) / r.dot(r)
    assert u[0] == 0 and u[1] == 0 and sp.simplify(u[2] - Gf / (4 * sp.pi * (y1 - y))) == 0         # upward for y1 > y
    dw = -u[2].subs(Gf, -dG)                                              # the shed strength is -(dGamma/dy) dy; w is positive DOWNWARD (T5)
    assert sp.simplify(dw - dG / (4 * sp.pi * (y1 - y))) == 0             # the integrand of (14.13)
    assert sp.simplify(dw - (-dG) / (4 * sp.pi * (y - y1))) == 0          # = the form written before (14.13), -(dGamma/dy) dy / (4 pi (y - y1))
    # one horseshoe: Gamma drops by Gamma at y = +s/2 and rises by Gamma at y = -s/2  ->  w(0) = Gamma / (pi s)
    w0 = (-G) / (4 * sp.pi * (0 - s / 2)) + (+G) / (4 * sp.pi * (0 - (-s / 2)))
    assert sp.simplify(w0 - G / (sp.pi * s)) == 0
    # the parabolic loading in closed form (the principal value by a symmetric gap), as used in the V1 test
    eps = sp.symbols("epsilon", positive=True)
    Gam = G0 * (1 - (2 * y / s) ** 2)
    integrand = sp.diff(Gam, y) / (4 * sp.pi * (y1 - y))
    assert sp.simplify(integrand - 2 * G0 / (sp.pi * s ** 2) * (1 + y1 / (y - y1))) == 0            # y/(y - y1) = 1 + y1/(y - y1)
    t = sp.symbols("t", positive=True)                                    # distance from the singular point on either side
    left = -sp.integrate(1 / t, (t, eps, y1 + s / 2))                     # int dy/(y - y1) over (-s/2, y1 - eps), with y = y1 - t
    right = sp.integrate(1 / t, (t, eps, s / 2 - y1))                     # ... over (y1 + eps, s/2), with y = y1 + t
    pv = sp.simplify(sp.logcombine(left + right, force=True))             # the two ln(eps) cancel
    assert not pv.has(eps)
    w_par = 2 * G0 / (sp.pi * s ** 2) * (s + y1 * pv)
    for yv in (-2.0, 0.0, 0.5, 3.0):
        num = float(w_par.subs({G0: 3.0, s: 7.0, y1: yv}))
        assert close(num, 2 * 3.0 / (PI * 49.0) * (7.0 + yv * math.log(abs((3.5 - yv) / (3.5 + yv)))), 1e-12, 1e-14)
        assert close(ch14.downwash(yv, lambda q: 3.0 * (1 - (2 * q / 7.0) ** 2), 7.0), num, 1e-9, 1e-11)


def test_strip_forces_V2_derivation():  # V2 — D13 (two-star): (14.14), dL = rho U Gamma dy, dD_i = rho w Gamma dy, (14.15); slip #2
    U, w, G, rho, dy, al = sp.symbols("U w Gamma rho dy alpha", positive=True)
    eps = sp.atan(w / U)                                                  # tan(eps) = w / U  (Fig. 14.23)
    Ue = sp.sqrt(U ** 2 + w ** 2)
    assert sp.simplify(Ue * sp.cos(eps) - U) == 0 and sp.simplify(Ue * sp.sin(eps) - w) == 0        # exact, not small-angle
    dLe = rho * Ue * G * dy                                               # Kutta–Zhukhovsky force normal to U_e
    assert sp.simplify(dLe * sp.cos(eps) - rho * U * G * dy) == 0 and sp.simplify(dLe * sp.sin(eps) - rho * w * G * dy) == 0
    x = sp.symbols("x", positive=True)
    assert sp.series(sp.atan(x), x, 0, 4).removeO() == x - x ** 3 / 3     # eps = w/U to within (w/U)^3 / 3
    assert sp.series(sp.tan(x) - sp.atan(x), x, 0, 4).removeO() == 2 * x ** 3 / 3                   # the printed "tan(w/U)" is off by twice that
    assert sp.simplify((al - x) - (al - sp.series(sp.atan(x), x, 0, 2).removeO())) == 0             # (14.14): alpha_e = alpha - w/U
    assert (rho * w * G * dy).subs(w, 0) == 0                             # no downwash (infinite span): no induced drag
    # the two integrals of (14.15) for the elliptic loading and its uniform downwash: (14.21) and (14.25)
    y, s, G1 = sp.symbols("y s Gamma_1", positive=True)
    Gam = G1 * sp.sqrt(1 - (2 * y / s) ** 2)
    L = sp.integrate(rho * U * Gam, (y, -s / 2, s / 2))
    Di = sp.integrate(rho * (G1 / (2 * s)) * Gam, (y, -s / 2, s / 2))
    assert sp.simplify(L - sp.pi * s / 4 * rho * U * G1) == 0 and sp.simplify(Di - sp.pi / 8 * rho * G1 ** 2) == 0
    num = ch14.strip_forces(7.0, 3.0, 40.0, 1.2, 0.5, exact=True)
    assert close(num["eps"], float(eps.subs({w: 3.0, U: 40.0}))) and close(num["dDi"], 1.2 * 3.0 * 7.0 * 0.5)


def _kirchhoff_energy_of_sheet(Gn, s, rho, M=1500):
    """Energy per length of a plane vortex sheet from the log kernel: T = -(rho / 4 pi) iint sigma sigma' ln|y - y'| dy dy'."""
    n = np.arange(1, len(Gn) + 1)
    g1, g2 = (np.arange(M) + 0.5) * PI / M, (np.arange(M) + 0.0) * PI / M + PI / (2 * M) * 0 + PI / M * 0.25      # staggered nodes
    dg = lambda g: -np.sum(n[:, None] * np.asarray(Gn)[:, None] * np.cos(n[:, None] * g[None, :]), axis=0)     # sigma dy = -dGamma
    Y1, Y2 = -(s / 2) * np.cos(g1), -(s / 2) * np.cos(g2)
    K = np.log(np.abs(Y1[:, None] - Y2[None, :]))
    return -rho / (4 * PI) * (dg(g1)[:, None] * dg(g2)[None, :] * K).sum() * (PI / M) ** 2


def test_induced_drag_V4_equals_the_kinetic_energy_left_in_the_wake():  # V4 — N38: D_i U = rate of wake kinetic energy
    s, rho = SPAN, 1.21
    for Gn in (GN, GN_ASYM, [1.4], [1.0, 0.0, 0.0, 0.0, 0.0, 0.3]):
        Di = ch14.induced_drag_from_series(Gn, rho)
        assert close(ch14.wake_kinetic_energy_per_length(Gn, s, rho), Di, 1e-10)                               # (1/2) rho int Gamma w_inf dy, w_inf = 2w
        assert close(ch14.wake_kinetic_energy_per_length(Gn, s, rho, m=64), Di, 1e-10)                         # midpoint rule: spectral
        assert close(ch14.wake_kinetic_energy_per_length(Gn, 3 * s, rho), Di, 1e-10)                           # the span drops out of (14.22)
        e1, e2 = _kirchhoff_energy_of_sheet(Gn, s, rho, 750), _kirchhoff_energy_of_sheet(Gn, s, rho, 1500)     # an independent energy integral
        assert abs(e2 - Di) < 1e-3 * Di and abs(e2 - Di) < 0.6 * abs(e1 - Di), (e1, e2, Di)                    # 2.7e-4, halving with the step
    # elliptic loading: the far wake is the flow of a flat plate of width s moving at w_inf = Gamma_1/s; its energy is
    # (1/2) (rho pi (s/2)^2) w_inf^2 with the plate's added mass equal to that of the circumscribed cylinder (chapter 6)
    G1 = 1.4
    assert close(0.5 * ch06.cylinder_added_mass(s / 2, rho) * (G1 / s) ** 2, ch14.induced_drag_from_series([G1], rho), 1e-12)


def test_vortex_pair_energy_V1_rankine_cores_and_the_field_outside_them():  # V1 + V2 — N67 (ours): (rho Gamma^2 / 2 pi)[ln(s/sigma) + 1/4]
    r, sig, G, rho = sp.symbols("r sigma Gamma rho", positive=True)
    core = sp.integrate(sp.Rational(1, 2) * rho * (G * r / (2 * sp.pi * sig ** 2)) ** 2 * 2 * sp.pi * r, (r, 0, sig))   # solid-body core (3.28)
    assert sp.simplify(core - rho * G ** 2 / (16 * sp.pi)) == 0 and sp.simplify(2 * core - rho * G ** 2 / (2 * sp.pi) / 4) == 0
    Gam, s, sigma, rh = 3.0, 2.0, 0.04, 1.2
    # outside the cores: |u| = Gamma s / (2 pi r1 r2) for the pair; integrate (1/2) rho |u|^2 over the half-plane x > 0 and double it
    f = lambda rr, th: 0.5 * rh * (Gam * s / (2 * PI)) ** 2 / (rr * (rr ** 2 + s ** 2 + 2 * rr * s * math.cos(th)))    # ... r dr dtheta / (r^2 r2^2)
    rmax = lambda th: np.inf if abs(th) <= PI / 2 else s / (2 * abs(math.cos(th)))
    out = 2 * dblquad(f, -PI, PI, lambda th: sigma, rmax, epsabs=1e-10, epsrel=1e-9)[0]
    assert close(out, rh * Gam ** 2 / (2 * PI) * math.log(s / sigma), 1e-3)                                    # leading order in sigma/s
    E = ch14.vortex_pair_energy_per_length(Gam, s, sigma, rh)
    assert close(E, out + 2 * float(core.subs({rho: rh, G: Gam})), 1e-3) and close(E, rh * Gam ** 2 / (2 * PI) * (math.log(s / sigma) + 0.25))
    sg = np.array([0.02, 0.04, 0.08, 0.16])
    Es = ch14.vortex_pair_energy_per_length(Gam, s, sg, rh)
    assert np.all(np.diff(Es) < 0) and close(np.diff(Es), -rh * Gam ** 2 / (2 * PI) * math.log(2), 1e-12)      # a fatter core holds less energy
    assert (Q_(1.2, "kg/m**3") * Q_(1.0, "m**2/s") ** 2).check("[energy]/[length]")


# =====================================================================================================================
# 12. C11 — the sine series and Glauert's integral (14.18)–(14.19); N51, slip #1
# =====================================================================================================================
def test_glauert_integral_V1_quadrature_closed_form_and_recursion():  # V1 — N51, D15 (ours, standard)
    g1s = (0.2, 0.7, PI / 3, PI / 2, 2.0, 2.9)
    for g1 in g1s:
        I_ = [ch14.glauert_integral(n, g1) for n in range(8)]
        for n in range(7):
            assert close(ch14.glauert_integral(n, g1, "quad"), I_[n], 1e-8, 1e-9), (n, g1)                    # singular part subtracted
            assert close(I_[n], PI * math.sin(n * g1) / math.sin(g1), 1e-13, 1e-15)
            if n >= 1:
                assert abs(I_[n + 1] + I_[n - 1] - 2 * math.cos(g1) * I_[n]) < 1e-12                           # the three-term recursion
        for n in (0, 1, 2, 5):                                            # the principal value by symmetric exclusion (a third route)
            pv = pv_symmetric(lambda g: math.cos(n * g) / (math.cos(g) - math.cos(g1)), 0.0, PI, g1, eps=1e-7)
            assert close(pv, I_[n], 2e-5, 2e-5), (n, g1, pv)              # the gap leaves eps cos(g1)/sin^2(g1) <= 2.5e-6
    assert close([ch14.glauert_integral(n, PI / 3) for n in range(5)], [0, PI, PI, 0, -PI], 1e-13, 1e-13)      # design D15's worked numbers
    assert ch14.glauert_integral(0, 1.1) == 0.0 and close(ch14.glauert_integral(1, 0.4), PI)                   # I_0 = 0, I_1 = pi
    for n in (1, 2, 3, 4):                                                # the tips: n and n (-1)^(n+1)  (l'Hopital)
        assert close(ch14.glauert_integral(n, 0.0), PI * n) and close(ch14.glauert_integral(n, PI), PI * n * (-1) ** (n + 1))
        assert close(ch14.glauert_integral(n, 1e-5), PI * n, 1e-8) and close(ch14.glauert_integral(n, PI - 1e-5), PI * n * (-1) ** (n + 1), 1e-8)
    with pytest.raises(ValueError):
        ch14.glauert_integral(1, 1.0, "simpson")


def test_downwash_series_V1_last_member_of_14_19_and_slip_1_the_printed_middle_member_has_the_wrong_sign():  # V1 — C11 (14.19), slip #1
    s = SPAN
    g1 = np.array([0.3, 0.9, PI / 2, 2.2, 2.8])
    for Gn in ([2.0], GN, GN_ASYM):
        n = np.arange(1, len(Gn) + 1)
        last = ch14.downwash_series(Gn, g1, s)
        assert close(last, np.sum(n * np.asarray(Gn) * np.sin(n * g1[:, None]), axis=1) / (2 * s * np.sin(g1)), 1e-13, 1e-15)   # as printed (last member)
        num = lambda g: float(np.sum(n * np.asarray(Gn) * np.cos(n * g)))
        corrected = np.array([pv_symmetric(lambda g: num(g) / (math.cos(g) - math.cos(q)), 0, PI, q) for q in g1]) / (2 * PI * s)
        printed = np.array([pv_symmetric(lambda g: num(g) / (math.cos(q) - math.cos(g)), 0, PI, q) for q in g1]) / (2 * PI * s)
        assert close(corrected, last, 2e-5, 2e-6)                         # denominator cos(gamma) - cos(gamma_1): equals the last member
        assert close(printed, -last, 2e-5, 2e-6)                          # as printed, cos(gamma_1) - cos(gamma): minus the last member
        pm = ch14.downwash_series(Gn, g1, s, printed_middle=True)
        assert close(pm, printed, 2e-5, 2e-6) and close(pm, -last, 1e-8, 1e-10)
        direct = ch14.downwash(ch14.span_position(g1, s), _gamma_fn(Gn, s), s)                                 # (14.13) itself decides
        assert close(direct, last, 1e-8, 1e-10) and not close(direct, pm, 1e-2, 1e-6)
    # the sign that physics demands: behind a lifting wing (Gamma > 0) the air goes down
    assert ch14.downwash_series([2.0], PI / 2, s) > 0 and ch14.horseshoe_downwash(0.0, 0.0, 2.0, s) > 0
    assert close(ch14.downwash_series([2.0], g1, s), 2.0 / (2 * s)) and close(ch14.downwash_series([2.0], g1, s, printed_middle=True), -2.0 / (2 * s), 1e-9)
    assert ch14.book_slips()["1"]["kind"] == "slip" and "cos γ − cos γ₁" in ch14.book_slips()["1"]["correct"]
    for n in (1, 2, 3):                                                   # tip values of each term: n and n (-1)^(n+1), times n Gamma_n / 2s
        e_n = np.zeros(n)
        e_n[-1] = 1.0
        assert close(ch14.downwash_series(e_n, 0.0, s), n * n / (2 * s)) and close(ch14.downwash_series(e_n, PI, s), n * n * (-1) ** (n + 1) / (2 * s))
        assert close(ch14.downwash_series(e_n, 1e-6, s), n * n / (2 * s), 1e-9)


def test_series_in_downwash_integral_V2_derivation():  # V2 — D14 (two-star): (14.13) with (14.18) -> first two members of (14.19), corrected
    s, g, g1 = sp.symbols("s gamma gamma_1", positive=True)
    G = sp.symbols("Gamma1:4", real=True)
    y, y1 = -(s / 2) * sp.cos(g), -(s / 2) * sp.cos(g1)
    assert sp.simplify((y1 - y) - (s / 2) * (sp.cos(g) - sp.cos(g1))) == 0                # the order of the cosines (slip #1)
    assert sp.simplify(sp.diff(y, g) - (s / 2) * sp.sin(g)) == 0 and y.subs(g, 0) == -s / 2 and y.subs(g, sp.pi) == s / 2
    Gam = sum(G[k] * sp.sin((k + 1) * g) for k in range(3))                               # (14.18)
    dGam = sp.diff(Gam, g)                                                                # dGamma = (dGamma/dy) dy
    assert sp.simplify(dGam - sum((k + 1) * G[k] * sp.cos((k + 1) * g) for k in range(3))) == 0
    integrand_y = (dGam / sp.diff(y, g)) / (y1 - y) * sp.diff(y, g) / (4 * sp.pi)         # (1/4 pi)(dGamma/dy) dy / (y1 - y), in gamma
    corrected = sum((k + 1) * G[k] * sp.cos((k + 1) * g) for k in range(3)) / (sp.cos(g) - sp.cos(g1)) / (2 * sp.pi * s)
    printed = sum((k + 1) * G[k] * sp.cos((k + 1) * g) for k in range(3)) / (sp.cos(g1) - sp.cos(g)) / (2 * sp.pi * s)
    assert sp.simplify(integrand_y - corrected) == 0 and sp.simplify(integrand_y + printed) == 0             # the printed integrand is its negative


def test_glauert_integral_V2_derivation():  # V2 — D15 (three-star): I_0 = 0, I_1 = pi, the recursion, the closed form, the tips
    g, g1 = sp.symbols("gamma gamma_1", positive=True)
    n = sp.symbols("n", integer=True, positive=True)
    # step 1: the antiderivative for n = 0
    F0 = sp.log(sp.sin((g1 + g) / 2) / sp.sin((g1 - g) / 2)) / sp.sin(g1)
    assert sp.simplify(sp.diff(F0, g) - 1 / (sp.cos(g) - sp.cos(g1))) == 0
    Fr = sp.log(sp.sin((g + g1) / 2) / sp.sin((g - g1) / 2)) / sp.sin(g1)                 # the same antiderivative to the right of gamma_1
    assert sp.simplify(sp.diff(Fr, g) - 1 / (sp.cos(g) - sp.cos(g1))) == 0
    # steps 2-3: symmetric gap -> I_0 = 0 (the end points give ln 1; the two logarithms of sin(eps/2) cancel)
    eps = sp.symbols("epsilon", positive=True)
    left = F0.subs(g, g1 - eps) - F0.subs(g, 0)
    right = Fr.subs(g, sp.pi) - Fr.subs(g, g1 + eps)
    step3 = sp.log(sp.sin(g1 - eps / 2) / sp.sin(g1 + eps / 2)) / sp.sin(g1)             # the line written in step 3
    assert F0.subs(g, 0) == 0 and sp.simplify(Fr.subs(g, sp.pi)) == 0                    # the end points give ln 1 = 0
    for gv, ev in ((0.9, 0.01), (2.3, 1e-3), (0.4, 0.05)):                               # the two logarithms of sin(eps/2) cancel
        assert abs(complex((left + right - step3).subs({g1: gv, eps: ev}).evalf())) < 1e-12
    assert sp.limit(step3, eps, 0) == 0
    # step 4: I_1 = pi + cos(gamma_1) I_0 = pi
    assert sp.simplify(sp.cos(g) / (sp.cos(g) - sp.cos(g1)) - (1 + sp.cos(g1) / (sp.cos(g) - sp.cos(g1)))) == 0
    # steps 5-6: the product-to-sum identity and the recursion I_{n+1} + I_{n-1} = 2 cos(gamma_1) I_n
    assert sp.simplify(sp.expand_trig(sp.cos((n + 1) * g) + sp.cos((n - 1) * g) - 2 * sp.cos(g) * sp.cos(n * g))) == 0
    assert sp.simplify(2 * sp.cos(g) * sp.cos(n * g) / (sp.cos(g) - sp.cos(g1))
                       - (2 * sp.cos(n * g) + 2 * sp.cos(g1) * sp.cos(n * g) / (sp.cos(g) - sp.cos(g1)))) == 0
    assert sp.integrate(sp.cos(n * g), (g, 0, sp.pi)) == 0                                # int cos(n gamma) over (0, pi) vanishes, n >= 1
    # steps 7-9: the candidate obeys the same recursion with the same two starting values
    S = lambda k: sp.pi * sp.sin(k * g1) / sp.sin(g1)
    assert S(0) == 0 and sp.simplify(S(1) - sp.pi) == 0
    assert sp.simplify(sp.expand_trig(S(n + 1) + S(n - 1) - 2 * sp.cos(g1) * S(n))) == 0
    # direct evaluation for n = 1..6: cos(k gamma) = T_k(cos gamma); (T_k(X) - T_k(X1))/(X - X1) is a polynomial; the rest is T_k(X1) I_0 = 0
    X, X1 = sp.symbols("X X_1")
    for k in range(1, 7):
        Tk = sp.chebyshevt(k, X)
        q = sp.cancel((Tk - Tk.subs(X, X1)) / (X - X1))
        Ik = sp.integrate(sp.expand(q.subs(X, sp.cos(g))), (g, 0, sp.pi))
        assert sp.simplify(sp.expand(Ik - sp.pi * sp.chebyshevu(k - 1, X1))) == 0         # = pi sin(k gamma_1) / sin(gamma_1)
        assert close(float((sp.pi * sp.chebyshevu(k - 1, X1)).subs(X1, math.cos(0.8))), ch14.glauert_integral(k, 0.8), 1e-13)
    # step 10: the last member of (14.19); step 11: the tips
    s_ = sp.symbols("s", positive=True)
    assert sp.simplify(S(n) / (2 * sp.pi * s_) - sp.sin(n * g1) / (2 * s_ * sp.sin(g1))) == 0
    for k in (1, 2, 3, 4):
        assert sp.limit(sp.sin(k * g1) / sp.sin(g1), g1, 0) == k and sp.limit(sp.sin(k * g1) / sp.sin(g1), g1, sp.pi) == (-1) ** (k + 1) * k


# =====================================================================================================================
# 13. C12 — the lifting-line equation (14.16)–(14.20) and its solver; N46–N48, N52, N63, N64, N68
# =====================================================================================================================
def test_lifting_line_matrix_V1_entry_by_entry_and_the_right_hand_side():  # V1 — C12 (14.20), N52
    s, K, U, al, N = 6.0, 5.7, 2.0, 0.07, 5
    cfun = lambda y: 1.0 - 0.3 * np.abs(2 * y / s)
    tw = lambda y: 0.01 * y ** 2
    A, ga = ch14.lifting_line_matrix(cfun, s, K, N)
    n, g2 = ch14.lifting_line_modes(N)
    assert close(ga, np.arange(1, N + 1) * PI / (N + 1)) and close(g2, ga) and close(n, np.arange(1, N + 1)) and A.shape == (N, N)
    r = ch14.lifting_line_rhs(cfun, s, U, al, K, N, beta=0.02, twist=tw)
    for j in range(N):
        yj = -(s / 2) * math.cos(ga[j])
        for k in range(N):
            assert close(A[j, k], (1 + (k + 1) * K * cfun(yj) / (4 * s * math.sin(ga[j]))) * math.sin((k + 1) * ga[j]), 1e-13, 1e-15)   # (14.20)
        assert close(r[j], 0.5 * K * U * cfun(yj) * (al + tw(yj) + 0.02), 1e-13)
    ns, gs = ch14.lifting_line_modes(N, symmetric=True)                    # symmetric loading: odd modes, stations on the port half
    assert close(ns, [1, 3, 5, 7, 9]) and close(gs, np.arange(1, N + 1) * PI / (2 * N)) and close(gs[-1], PI / 2)
    As, _ = ch14.lifting_line_matrix(cfun, s, K, N, symmetric=True)
    assert close(As[2, 3], (1 + 7 * K * cfun(-(s / 2) * math.cos(gs[2])) / (4 * s * math.sin(gs[2]))) * math.sin(7 * gs[2]), 1e-13)
    assert 0 < ga.min() and ga.max() < PI and 0 < gs.min()                # the tips (sin gamma = 0) are never stations (trap T10)
    # the same matrix from every accepted form of the chord and of K
    yst = -(s / 2) * np.cos(ga)
    A1, _ = ch14.lifting_line_matrix(cfun(yst), s, np.full(N, K), N)
    A2, _ = ch14.lifting_line_matrix(cfun, s, lambda y: K + 0 * y, N)
    assert np.array_equal(A, A1) and np.array_equal(A, A2)
    P = ch14.planform("tapered", s, 1.0, taper=0.7)
    assert close(ch14.lifting_line_matrix(P, s, K, N)[0], A, 1e-14)
    assert close(ch14.lifting_line_matrix(1.3, s, K, N)[0], ch14.lifting_line_matrix(lambda y: 1.3 + 0 * y, s, K, N)[0], 1e-15)
    with pytest.raises(ValueError):
        ch14.lifting_line_matrix(np.ones(N + 2), s, K, N)
    sol = ch14.solve_lifting_line(cfun, s, U, al, K=K, beta=0.02, twist=tw, N=N, rho=1.3)
    assert np.array_equal(sol["matrix"], A) and np.array_equal(sol["rhs"], r) and close(sol["Gn"], np.linalg.solve(A, r), 1e-12, 1e-14)
    assert close(sol["cond"], np.linalg.cond(A)) and close(sol["gamma"], ga) and close(sol["y"], yst) and list(sol["n"]) == [1, 2, 3, 4, 5]
    assert close(sol["L"], ch14.lift_from_series(sol["Gn"], s, 1.3, U)) and close(sol["Di"], ch14.induced_drag_from_series(sol["Gn"], 1.3))
    assert sol["e"] == sol["e_span"] and sol["Gamma1"] == sol["Gn"][0]
    # the area of a chord given as a callable is found by a Gauss rule in the span angle applied to each half-span
    # separately, so the kink of this tapered chord at mid-span sits on an end point and the area is exact: measured
    # 5.10000000000001 against the closed form 5.1 (2.0e-15 relative; it was 2.2e-6 with one rule across the kink)
    assert close(sol["A"], P["A"], 1e-12) and close(sol["AR"], P["AR"], 1e-12)
    for s2, c2, tp in ((12.0, 1.6, 0.45), (7.0, 2.1, 0.2), (5.0, 3.0, 0.0), (9.0, 0.8, 1.0)):       # other tapers, a pointed tip, no taper: at most 1.5e-14
        P2 = ch14.planform("tapered", s2, c2, taper=tp)
        assert close(ch14.solve_lifting_line(P2["chord"], s2, U, al, K=K, N=N)["A"], P2["A"], 1e-12), (s2, c2, tp)
    solP = ch14.solve_lifting_line(P, s, U, al, K=K, beta=0.02, twist=tw, N=N, rho=1.3)
    assert solP["A"] == P["A"] and close(solP["AR"], P["AR"], 1e-14) and close(solP["Gn"], sol["Gn"], 1e-13, 1e-15)
    # the two routes (callable chord, planform dict) now give the same coefficients: measured 2.0e-15 relative in C_L and
    # C_Di (they differed by the area error, 2.2e-6, before); the forces never depended on the area
    assert close(solP["CL"], sol["CL"], 1e-12) and close(solP["CDi"], sol["CDi"], 1e-12) and close(solP["L"], sol["L"], 1e-14) and close(solP["Di"], sol["Di"], 1e-14)
    smooth = ch14.solve_lifting_line(lambda y: 1.2 * np.sqrt(np.clip(1 - (2 * y / s) ** 2, 0, None)), s, U, al, K=K, N=N)
    assert close(smooth["A"], PI * 1.2 * s / 4, 1e-12) and close(ch14.solve_lifting_line(1.3, s, U, al, K=K, N=N)["A"], 1.3 * s, 1e-15)
    assert close(sol["CL"], sol["L"] / (0.5 * 1.3 * U ** 2 * sol["A"])) and close(sol["CDi"], sol["Di"] / (0.5 * 1.3 * U ** 2 * sol["A"]))
    assert close(sol["CL"], ch14.solve_lifting_line(cfun, s, U, al, K=K, beta=0.02, twist=tw, N=N)["CL"], 1e-13)   # coefficients ignore rho


def test_lifting_line_V1_elliptic_planform_untwisted_gives_the_elliptic_loading_exactly():  # V1 — C12/C14, N58 (14.27), N59, D20
    for s, c0, K, al, be in ((12.0, 1.6, 2 * PI, 0.05, 0.0), (7.0, 2.1, 5.5, 0.03, 0.02), (30.0, 1.0, 6.1, -0.02, 0.05)):
        P = ch14.planform("elliptic", s, c0)
        for N, sym in ((4, False), (12, False), (9, True)):
            r = ch14.solve_lifting_line(P, s, 3.0, al, K=K, beta=be, N=N, symmetric=sym, rho=1.1)
            assert np.max(np.abs(r["Gn"][1:])) < 1e-12 * abs(r["Gn"][0])                                       # only Gamma_1: the ellipse (14.23)
            slope = K / (1 + K / (PI * P["AR"]))
            assert close(r["CL"], slope * (al + be), 1e-12) and close(r["slope"], slope, 1e-12) and close(r["e_span"], 1.0, 1e-12)
            assert close(ch14.finite_wing_lift_slope(K, P["AR"]), slope) and close(r["CDi"], r["CL"] ** 2 / (PI * P["AR"]), 1e-12)   # (14.26)
            assert close(r["w"], r["Gn"][0] / (2 * s), 1e-11)                                                  # uniform downwash (14.24)
            assert close(r["Gamma"], ch14.elliptic_circulation(r["y"], s, r["Gn"][0]), 1e-11)
            ew = ch14.elliptic_wing(al, P["AR"], K=K, beta=be)
            assert close([ew["CL"], ew["CDi"], ew["slope"]], [r["CL"], r["CDi"], r["slope"]], 1e-12)
            assert close(ew["w_over_U"], r["w"][0] / 3.0, 1e-11) and close(ew["alpha_e"], al - ew["w_over_U"]) and close(ew["w_over_U"], ew["CL"] / (PI * P["AR"]))
            # (14.16), (14.17) at the stations: C_L = K[alpha - w/U + beta] and Gamma = (1/2) U c C_L
            cl = ch14.section_lift_coefficient(al, r["w"], 3.0, be, K=K)
            assert close(ch14.circulation_from_cl(cl, 3.0, P["chord"](r["y"])), r["Gamma"], 1e-11)
            assert close(cl, r["CL"], 1e-11)                               # every section of the untwisted elliptic wing flies at the wing's C_L
    assert close(ch14.finite_wing_lift_slope(2 * PI, 6.0), 2 * PI / (1 + 2 / 6.0)) and round(ch14.finite_wing_lift_slope(2 * PI, 6.0), 2) == 4.71   # design D20
    assert round(ch14.finite_wing_lift_slope(2 * PI, 9.6), 2) == 5.20 and close(ch14.finite_wing_lift_slope(2 * PI, 1e12), 2 * PI, 1e-10)       # Lambda -> inf
    AR = np.array([1e-3, 1e-2])
    assert close(ch14.finite_wing_lift_slope(5.0, AR), PI * AR, 2e-2)      # Lambda -> 0: slope -> pi Lambda
    assert ch14.section_lift_coefficient(0.1, 0.0, 1.0, 0.03, K=2 * PI) == ch14.zhukhovsky_lift_coefficient(0.1, 0.03)    # w = 0, K = 2 pi: (14.12)


def test_lifting_line_V3_rectangular_wing_converges_and_14_17_holds_between_the_stations():  # V3 — C12, N48, N52 (the demo wing)
    s, c, U, al, K = I["wd_s"], I["wd_c"], I["wd_U"], I["wd_alpha"], I["wd_K"]
    Ns = [4, 8, 16, 32, 64, 128]
    sol = [ch14.solve_lifting_line(c, s, U, al, K=K, N=N) for N in Ns]
    CL, CDi, es = [q["CL"] for q in sol], [q["CDi"] for q in sol], [q["e_span"] for q in sol]
    dCL, dCD = np.abs(np.diff(CL)), np.abs(np.diff(CDi))
    oL, oD = observed_order([1 / N for N in Ns[:-1]], dCL), observed_order([1 / N for N in Ns[:-1]], dCD)
    assert oL > 3.0 and oD > 3.0, (oL, oD, CL, CDi)                       # no design order is stated; measured about 3.7 and 3.8
    assert dCL[-1] < 1e-7 and dCD[-1] < 1e-8 and abs(es[-1] - es[-2]) < 1e-6
    assert round(es[-1], 4) == 0.9539 and round(sol[-1]["slope"], 4) == 4.5304                                 # the analyst's converged values
    assert sol[-1]["slope"] < ch14.finite_wing_lift_slope(K, s / c) and es[-1] < 1                             # worse than the elliptic wing
    assert close(sol[-1]["CDi"], ch14.induced_drag_coefficient(sol[-1]["CL"], s / c, es[-1]), 1e-12)
    for q in sol:
        assert np.max(np.abs(q["Gn"][1::2])) < 1e-12 * q["Gn"][0]         # symmetric wing: even coefficients vanish
        res = ch14.lifting_line_residual(q["Gn"], q["gamma"], c, s, U, al, K=K)
        assert np.max(np.abs(res)) < 1e-11 * q["Gn"][0]                   # (14.17) holds at the collocation stations
        assert np.all(q["Gamma"] > 0) and np.all(np.diff(q["w"][: q["w"].size // 2]) < 0)                      # more downwash toward the tips
    gg = np.linspace(0.15, PI - 0.15, 600)                                 # between the stations the residual falls as N grows
    rmax = [np.max(np.abs(ch14.lifting_line_residual(q["Gn"], gg, c, s, U, al, K=K))) / q["Gn"][0] for q in sol]
    assert all(b < a for a, b in zip(rmax, rmax[1:])) and rmax[-1] < 1e-4, rmax
    # independent of the series formula: the downwash of the solved loading by direct quadrature of (14.13)
    q = sol[3]
    yq = np.array([-2.2, -0.9, 0.0, 1.4])
    wq = ch14.downwash(yq, _gamma_fn(q["Gn"], s), s)
    lhs = ch14.circulation_series(q["Gn"], ch14.span_angle(yq, s))
    assert close(lhs, 0.5 * K * U * c * (al - wq / U), 1e-4)               # Gamma = (K/2) U c [alpha - w/U]  (14.17): 2.2e-5 at N = 32
    # the symmetric (odd-mode) solver gives the same wing
    qs = ch14.solve_lifting_line(c, s, U, al, K=K, N=64, symmetric=True)
    assert close(qs["CL"], CL[-1], 1e-7) and close(qs["e_span"], es[-1], 1e-6) and qs["Gn"].size == 127 and np.all(qs["Gn"][1::2] == 0)
    # V7: a very slender rectangular wing tends to the two-dimensional result K alpha
    slender = [ch14.solve_lifting_line(c, AR * c, U, al, K=K, N=48, symmetric=True)["CL"] / (K * al) for AR in (6, 20, 100, 1000)]
    assert all(b > a for a, b in zip(slender, slender[1:])) and slender[-1] > 0.99 and slender[-1] < 1.0
    # an asymmetric wing (aileron-like twist): even coefficients appear, the lift still comes from Gamma_1 alone
    qa = ch14.solve_lifting_line(c, s, U, al, K=K, twist=lambda y: 0.02 * np.sign(y), N=40, rho=1.2)
    assert np.max(np.abs(qa["Gn"][1::2])) > 1e-3 * qa["Gn"][0] and close(qa["L"], PI * s / 4 * 1.2 * U * qa["Gn"][0], 1e-13)


def test_lifting_line_V1_twist_gives_an_elliptic_loading_on_a_non_elliptic_planform():  # V1 — N63, N58 (trap T11)
    c, s, K, U, ar = 1.3, 8.0, 5.6, 40.0, np.deg2rad(3.0)
    t = ch14.elliptic_twist_rectangular(ar, c, s, K=K)
    assert close(t["alpha_t"], K * c * ar / (4 * s)) and close(t["Gamma1_over_U"], 0.5 * K * c * ar) and close(t["AR"], s / c)
    r = ch14.solve_lifting_line(c, s, U, t["alpha_t"], K=K, twist=lambda y: ar * np.sqrt(np.clip(1 - (2 * y / s) ** 2, 0, None)), N=14)
    assert np.max(np.abs(r["Gn"][1:])) < 1e-12 * r["Gn"][0] and close(r["Gn"][0], 0.5 * K * U * c * ar, 1e-12)  # only Gamma_1 = (1/2) K U c alpha_r
    assert close(r["CL"], t["CL"], 1e-12) and close(r["CDi"], t["CDi"], 1e-12) and close(r["e_span"], 1.0, 1e-12)
    assert close(r["w"], r["Gn"][0] / (2 * s), 1e-11) and close(t["alpha_t"], r["w"][0] / U, 1e-11)             # alpha_t is the induced angle
    r_no = ch14.solve_lifting_line(c, s, U, t["alpha_t"] + ar * PI / 4, K=K, N=14)                             # the same wing without twist
    assert r_no["e_span"] < 0.999
    # any planform: the twist that gives a chosen lift coefficient with an elliptic loading
    P = ch14.planform("tapered", s, 1.6, taper=0.4)
    for CLt in (0.3, 0.8):
        tw = lambda y: ch14.twist_for_elliptic_loading(P["chord"], y, s, K=K, CL_target=CLt, A=P["A"])["alpha_local"]
        rt = ch14.solve_lifting_line(P, s, U, 0.0, K=K, twist=tw, N=16)
        assert np.max(np.abs(rt["Gn"][1:])) < 1e-11 * rt["Gn"][0] and close(rt["CL"], CLt, 1e-12) and close(rt["e_span"], 1.0, 1e-11)
        d = ch14.twist_for_elliptic_loading(P["chord"], np.array([0.0, 1.0]), s, K=K, CL_target=CLt)            # area from the chord function
        assert close(d["Gamma1_over_U"], 2 * CLt * P["A"] / (PI * s), 1e-5) and close(d["alpha_i"], d["Gamma1_over_U"] / (2 * s))   # kinked chord: 5e-6
        dA = ch14.twist_for_elliptic_loading(P["chord"], np.array([0.0, 1.0]), s, K=K, CL_target=CLt, A=P["A"])
        assert close(dA["Gamma1_over_U"], 2 * CLt * P["A"] / (PI * s), 1e-14) and close(dA["alpha_i"], CLt / (PI * P["AR"]), 1e-14)   # C_L / (pi Lambda)
    # elliptic planform: the required "twist" is no twist at all
    Pe = ch14.planform("elliptic", s, 1.4)
    te = ch14.twist_for_elliptic_loading(Pe["chord"], np.linspace(-3.5, 3.5, 9), s, K=K, CL_target=0.5, A=Pe["A"])["alpha_local"]
    assert np.ptp(te) < 1e-12 and close(te[0], 0.5 / ch14.finite_wing_lift_slope(K, Pe["AR"]), 1e-12)


def test_lifting_line_equation_V2_derivation():  # V2 — D16 (two-star): (14.16) + Gamma = U c C_L / 2 -> (14.17); with (14.18), (14.19) -> (14.20)
    K, U, c, s, al, be, g = sp.symbols("K U c s alpha beta gamma", positive=True)
    w, CL, Gam, rho = sp.symbols("w C_L Gamma rho", positive=True)
    G = sp.symbols("Gamma1:5", real=True)
    cl = K * (al - w / U + be)                                             # (14.16)
    circ = sp.solve(sp.Eq(rho * U * Gam, CL * rho * U ** 2 * c / 2), Gam)[0]              # L' = rho U Gamma = C_L (1/2) rho U^2 c
    assert sp.simplify(circ - U * c * CL / 2) == 0
    eq17 = circ.subs(CL, cl)                                               # (14.17): Gamma = (K/2) U c [alpha - w/U + beta]
    assert sp.simplify(eq17 - K / 2 * U * c * (al - w / U + be)) == 0
    series = sum(G[k] * sp.sin((k + 1) * g) for k in range(4))             # (14.18)
    wser = sum((k + 1) * G[k] * sp.sin((k + 1) * g) for k in range(4)) / (2 * s * sp.sin(g))     # (14.19), last member
    lhs = sum((1 + (k + 1) * K * c / (4 * s * sp.sin(g))) * G[k] * sp.sin((k + 1) * g) for k in range(4))   # right member of (14.20)
    assert sp.simplify((series - eq17.subs(w, wser)) - (lhs - K / 2 * U * c * (al + be))) == 0
    # limits of the Check: s -> infinity gives the two-dimensional Gamma; an elliptic chord cancels sin(gamma) and one term solves it
    assert sp.limit(lhs, s, sp.oo) == series
    c0 = sp.symbols("c_0", positive=True)
    one = sp.solve(sp.Eq((1 + K * c0 * sp.sin(g) / (4 * s * sp.sin(g))) * G[0] * sp.sin(g), K / 2 * U * c0 * sp.sin(g) * (al + be)), G[0])[0]
    assert sp.simplify(one - K * U * c0 * (al + be) / (2 * (1 + K * c0 / (4 * s)))) == 0 and one.free_symbols.isdisjoint({g})


def test_rolling_moment_V1_closed_form_quadrature_and_limits():  # V1 + V7 — N64 (ours; w upward in that exercise, trap T5)
    rho, U, K, s, c, Gv, ell = I["rho"], I["U_cruise"], I["K"], I["s"], I["c_mean"], I["Gamma_v"], I["ell"]
    m = ch14.wake_vortex_rolling_moment(rho, U, K, s, c, Gv, ell)
    wv = lambda y: Gv / (2 * PI * y) * (1 - math.exp(-abs(y) / ell))       # the vortex's vertical velocity
    mine = quad(lambda y: rho * U * y * (0.5 * K * U * c * wv(y) / U), -s / 2, s / 2, points=[0.0], epsabs=1e-10, epsrel=1e-13)[0]   # int rho U y Gamma dy
    assert close(m["M"], mine, 1e-11) and close(m["M_quad"], m["M"], 1e-10)
    assert close(m["M"], rho * U * K * c * Gv / (2 * PI) * (s / 2 - ell * (1 - math.exp(-s / (2 * ell)))), 1e-13)
    assert close(ch14.wake_vortex_rolling_moment(rho, U, K, s, c, -Gv, ell)["M"], -m["M"])                     # odd in the vortex's sign
    tiny = ch14.wake_vortex_rolling_moment(rho, U, K, s, c, Gv, 1e-9)
    assert close(tiny["M"], rho * U * K * c * Gv * s / (4 * PI), 1e-8) and close(tiny["ratio"], 1.0, 1e-8) and 0 < m["ratio"] < 1
    assert ch14.wake_vortex_rolling_moment(rho, U, K, s, c, Gv, 2 * ell)["M"] < m["M"]                         # a fatter core rolls the follower less


# =====================================================================================================================
# 14. C13 — lift and induced drag from the coefficients (14.21)–(14.22); N53, N54
# =====================================================================================================================
def test_series_forces_V7_every_coefficient_adds_drag_and_span_efficiency_is_at_most_one():  # V7 + V1 — C13 (14.22), N54
    rho, U, s = 1.2, 50.0, 12.0
    rng = np.random.default_rng(14)
    base = ch14.induced_drag_from_series([3.0], rho)
    assert ch14.span_efficiency([3.0]) == 1.0 and close(base, PI / 8 * rho * 9.0)
    for _ in range(40):
        Gn = np.concatenate([[3.0], rng.normal(0, 0.6, 6)])
        e = ch14.span_efficiency(Gn)
        n = np.arange(1, 8)
        assert 0 < e < 1 and close(e, Gn[0] ** 2 / np.sum(n * Gn ** 2))
        assert ch14.induced_drag_from_series(Gn, rho) > base and close(ch14.induced_drag_from_series(Gn, rho), base / e, 1e-13)   # same lift, more drag
        assert close(ch14.lift_from_series(Gn, s, rho, U), ch14.lift_from_series([3.0], s, rho, U), 1e-15)     # the lift sees Gamma_1 only
        A = 15.0
        q = 0.5 * rho * U ** 2 * A
        CL, CDi = ch14.lift_from_series(Gn, s, rho, U) / q, ch14.induced_drag_from_series(Gn, rho) / q
        assert close(CDi, ch14.induced_drag_coefficient(CL, s ** 2 / A, e), 1e-13)                             # C_Di = C_L^2 / (pi Lambda e_span)
    assert close(ch14.induced_drag_from_series([0.0, 0.0, 2.0], rho), PI / 8 * rho * 3 * 4.0) and math.isnan(ch14.span_efficiency([0.0, 0.0]))
    assert ch14.lift_from_series([0.0, 5.0, 7.0], s, rho, U) == 0.0        # no Gamma_1, no lift — but still drag
    assert close(ch14.induced_drag_from_series(2 * GN, rho), 4 * ch14.induced_drag_from_series(GN, rho))       # quadratic in the loading
    assert close(ch14.span_efficiency(GN), ch14.span_efficiency(-3.0 * GN), 1e-14)


def test_lift_from_series_V2_derivation():  # V2 — D17 (two-star): (14.15) with (14.18) -> L = (pi s / 4) rho U Gamma_1
    g, s, rho, U = sp.symbols("gamma s rho U", positive=True)
    G = sp.symbols("Gamma1:6", real=True)
    Gam = sum(G[k] * sp.sin((k + 1) * g) for k in range(5))
    dy = (s / 2) * sp.sin(g)                                               # dy = (s/2) sin(gamma) d gamma
    L = sp.integrate(sp.expand(rho * U * Gam * dy), (g, 0, sp.pi))
    assert sp.simplify(L - sp.pi * s / 4 * rho * U * G[0]) == 0            # only n = 1 survives
    n, m = sp.symbols("n m", integer=True, positive=True)
    assert sp.integrate(sp.sin(g) * sp.sin(g), (g, 0, sp.pi)) == sp.pi / 2 and all(sp.integrate(sp.sin(k * g) * sp.sin(g), (g, 0, sp.pi)) == 0 for k in range(2, 7))
    y = sp.symbols("y", real=True)                                         # the Check: rho U x the area of a half-ellipse
    assert sp.simplify(sp.integrate(G[0] * sp.sqrt(1 - (2 * y / s) ** 2), (y, -s / 2, s / 2)) - sp.pi * s * G[0] / 4) == 0


def test_induced_drag_from_series_V2_derivation():  # V2 — D18 (three-star): (14.15), (14.18), (14.19) -> D_i = (pi/8) rho sum n Gamma_n^2
    g, s, rho = sp.symbols("gamma s rho", positive=True)
    G = sp.symbols("Gamma1:5", real=True)
    Gam = sum(G[m] * sp.sin((m + 1) * g) for m in range(4))                # (14.18), index m
    wsin = sum((n + 1) * G[n] * sp.sin((n + 1) * g) for n in range(4)) / (2 * s)          # step 2: w sin(gamma), index n
    w = wsin / sp.sin(g)                                                   # (14.19), last member
    start = rho * w * Gam * (s / 2) * sp.sin(g)                            # step 1: rho w Gamma dy with dy = (s/2) sin(gamma) d gamma
    step3 = rho / 4 * sum((n + 1) * G[n] * sp.sin((n + 1) * g) for n in range(4)) * Gam   # step 3: the span has cancelled, sin(gamma) once
    assert sp.simplify(start - step3) == 0
    terms = [rho / 4 * (n + 1) * G[n] * G[m] * sp.integrate(sp.sin((n + 1) * g) * sp.sin((m + 1) * g), (g, 0, sp.pi))
             for n in range(4) for m in range(4)]                          # steps 4-6: sixteen pairings
    for n in range(4):                                                     # step 7: orthogonality on the half range
        for m in range(4):
            assert sp.integrate(sp.sin((n + 1) * g) * sp.sin((m + 1) * g), (g, 0, sp.pi)) == (sp.pi / 2 if n == m else 0)
    Di = sp.expand(sum(terms))
    target = sp.pi / 8 * rho * sum((n + 1) * G[n] ** 2 for n in range(4))  # steps 8-9: (14.22)
    assert sp.simplify(Di - target) == 0 and sp.simplify(sp.integrate(sp.expand(step3), (g, 0, sp.pi)) - target) == 0
    assert not Di.has(s)                                                   # the span does not appear
    one_term = sp.integrate(sp.sin(g) * sp.sin(g) / sp.sin(g), (g, 0, sp.pi))             # trap: cancelling sin(gamma) twice (one term is enough)
    assert one_term == 2 and one_term != sp.pi / 2
    same_index = sp.pi / 8 * rho * sum(G[n] ** 2 for n in range(4))        # trap: losing the factor n
    assert sp.simplify(same_index - target) != 0
    num = {G[0]: 12.0, G[1]: 0.0, G[2]: 1.0, G[3]: 0.0, rho: 1.2}
    assert close(float(target.subs(num)), ch14.induced_drag_from_series([12.0, 0.0, 1.0, 0.0], 1.2), 1e-13)


# =====================================================================================================================
# 15. C14 — the elliptic loading (14.23)–(14.27); N55–N59, N62
# =====================================================================================================================
def test_elliptic_loading_V1_14_23_to_14_26_and_the_cited_forms_V5():  # V1 + V5 — C14 (14.26), N55–N57
    rho, U, s, A, W = I["rho"], I["U_cruise"], I["s"], I["A"], I["W"]
    G1 = 4 * W / (PI * s * rho * U)                                        # (14.21) with L = W
    assert close(ch14.lift_from_series([G1], s, rho, U), W, 1e-13)
    Di = ch14.induced_drag_from_series([G1], rho)
    assert close(Di, PI / 8 * rho * G1 ** 2) and close(ch14.elliptic_induced_drag(W, rho, U, s), Di, 1e-13)    # both members of (14.25)
    assert close(Di, rho * ch14.elliptic_downwash(G1, s) * G1 * PI * s / 4, 1e-13)                             # = rho w int Gamma dy with w = Gamma_1/2s
    q = 0.5 * rho * U ** 2 * A
    CL, CDi, AR = W / q, Di / q, ch14.aspect_ratio(s, A)
    assert close(CDi, CL ** 2 / (PI * AR), 1e-13) and close(ch14.induced_drag_coefficient(CL, AR), CDi, 1e-13)  # (14.26)
    assert close(CDi, CL ** 2 / (PI * (s ** 2 / A)), 1e-13)                # the middle member of (14.26)
    bm = bench()
    assert "C_L^2 / (pi AR e)" in bm["induced_drag_coefficient"]["form"] and bm["induced_drag_coefficient"]["e_elliptic"] == 1.0
    for e in (1.0, 0.9, 0.75):                                             # the cited form C_D,i = C_L^2 / (pi AR e)
        assert close(ch14.induced_drag_coefficient(CL, AR, e), CL ** 2 / (PI * AR * e), 1e-14)
    assert ch14.induced_drag_coefficient(CL, AR, bm["induced_drag_coefficient"]["e_elliptic"]) == ch14.induced_drag_coefficient(CL, AR)
    assert "L^2 / (1/2 rho V^2 pi b^2)" in bm["induced_drag_force"]["form"]
    assert close(ch14.elliptic_induced_drag(W, rho, U, s), W ** 2 / (0.5 * rho * U ** 2 * PI * s ** 2), 1e-14)  # the cited D_i = L^2 / (1/2 rho V^2 pi b^2)
    Ls = np.array([1e3, 2e3, 4e3])
    assert close(ch14.elliptic_induced_drag(Ls, rho, U, s) / Ls ** 2, 2 / (PI * rho * U ** 2 * s ** 2), 1e-14)  # quadratic in the lift
    assert close(ch14.elliptic_induced_drag(W, rho, U, 2 * s), Di / 4)     # span, not area, is what counts
    assert close(ch14.elliptic_induced_drag(W, rho, 2 * U, s), Di / 4) and ch14.induced_drag_coefficient(0.0, AR) == 0.0
    assert round(ch14.elliptic_induced_drag(1e4, 1.25, 40.0, 10.0)) == 318 and round(ch14.induced_drag_coefficient(0.5, 5.0), 4) == 0.0159      # design D19
    # the loading and its uniform downwash; a quarter of the way along each wing the circulation is Gamma_1 sqrt(3)/2
    y = np.linspace(-s / 2, s / 2, 41)
    assert close(ch14.elliptic_circulation(y, s, G1), G1 * np.sin(ch14.span_angle(y, s)), 1e-12, 1e-7 * G1)
    assert close(ch14.elliptic_circulation(s / 4, s, G1), G1 * math.sqrt(3) / 2)
    assert close(ch14.downwash(y[5:-5:6], lambda yy: ch14.elliptic_circulation(yy, s, G1), s), G1 / (2 * s), 1e-9)
    # N58 (14.27): dL = rho U Gamma dy = C_L (1/2) rho U^2 c dy with uniform w -> c proportional to Gamma: an elliptic planform
    P = ch14.planform("elliptic", s, 4 * A / (PI * s))
    assert close(P["A"], A, 1e-13)
    cl_sec = rho * U * ch14.elliptic_circulation(y[1:-1], s, G1) / (0.5 * rho * U ** 2 * P["chord"](y[1:-1]))
    assert np.ptp(cl_sec) < 1e-10 * cl_sec[0] and close(cl_sec[0], CL, 1e-12)


def test_elliptic_optimum_V2_derivation():  # V2 — D19 (two-star): fixed L fixes Gamma_1; (14.23)–(14.26)
    L, rho, U, s, A, g = sp.symbols("L rho U s A gamma", positive=True)
    G = sp.symbols("Gamma1:5", real=True)
    G1 = sp.solve(sp.Eq(sp.pi * s / 4 * rho * U * G[0], L), G[0])[0]       # (14.21) inverted
    assert sp.simplify(G1 - 4 * L / (sp.pi * s * rho * U)) == 0
    Di_all = sp.pi / 8 * rho * sum((n + 1) * G[n] ** 2 for n in range(4))  # (14.22)
    excess = sp.expand(Di_all - sp.pi / 8 * rho * G[0] ** 2)
    assert excess == sp.expand(sp.pi / 8 * rho * (2 * G[1] ** 2 + 3 * G[2] ** 2 + 4 * G[3] ** 2))      # a sum of non-negative terms
    assert all(sp.diff(Di_all, G[n]).subs({G[1]: 0, G[2]: 0, G[3]: 0}) == 0 for n in (1, 2, 3)) and all(sp.diff(Di_all, G[n], 2) > 0 for n in (1, 2, 3))
    y = -(s / 2) * sp.cos(g)
    assert sp.simplify(sp.sqrt(1 - (2 * y / s) ** 2) - sp.Abs(sp.sin(g))) == 0            # (14.23): Gamma_1 sin(gamma) = Gamma_1 sqrt(1 - (2y/s)^2)
    w = (1 * G[0] * sp.sin(g) / sp.sin(g)) / (2 * s)                       # (14.19) with one term
    assert sp.simplify(w - G[0] / (2 * s)) == 0                            # (14.24)
    Di = (sp.pi / 8 * rho * G[0] ** 2).subs(G[0], G1)
    assert sp.simplify(Di - 2 * L ** 2 / (sp.pi * rho * U ** 2 * s ** 2)) == 0            # (14.25)
    q = rho * U ** 2 * A / 2
    assert sp.simplify(Di / q - (L / q) ** 2 / (sp.pi * (s ** 2 / A))) == 0               # (14.26): C_Di = C_L^2 / (pi Lambda)
    assert sp.limit(Di, s, sp.oo) == 0                                     # infinite span: no induced drag


def test_finite_wing_lift_slope_V2_derivation():  # V2 — D20 (two-star): (14.27), then (14.17) with (14.24) -> dC_L/d alpha = K / (1 + K / (pi Lambda))
    K, U, s, c0, al, be, rho, g = sp.symbols("K U s c_0 alpha beta rho gamma", positive=True)
    G1 = sp.symbols("Gamma_1", positive=True)
    c = c0 * sp.sin(g)                                                     # elliptic planform, c = c_0 sqrt(1 - (2y/s)^2)
    eq = sp.Eq(G1 * sp.sin(g), K / 2 * U * c * (al - (G1 / (2 * s)) / U + be))            # (14.17) with (14.23), (14.24)
    sol = sp.solve(eq, G1)[0]
    assert sol.free_symbols.isdisjoint({g})                                # one constant satisfies it at every station: c proportional to Gamma
    A = sp.pi * c0 * s / 4                                                 # area of the ellipse
    Lam = s ** 2 / A
    CL = (sp.pi * s / 4 * rho * U * sol) / (rho * U ** 2 * A / 2)          # (14.21) over (1/2) rho U^2 A
    assert sp.simplify(CL - K * (al + be) / (1 + K / (sp.pi * Lam))) == 0
    assert sp.simplify(CL.subs(K, 2 * sp.pi) - 2 * sp.pi * (al + be) / (1 + 2 / Lam)) == 0     # C_L,3D = C_L,2D / (1 + 2 / Lambda) for K = 2 pi
    slope = sp.diff(CL, al)
    AR = sp.symbols("Lambda", positive=True)
    sl = K / (1 + K / (sp.pi * AR))
    assert sp.limit(sl, AR, sp.oo) == K and sp.limit(sl / AR, AR, 0) == sp.pi and sp.simplify(slope - sl.subs(AR, Lam)) == 0
    # (14.27): rho U Gamma dy = C_L (1/2) rho U^2 c dy  =>  the section C_L is the same at every station of this wing
    cl_sec = rho * U * sol * sp.sin(g) / (rho * U ** 2 * c / 2)
    assert sp.simplify(cl_sec - CL) == 0


# =====================================================================================================================
# 16. C15 — the drag breakdown and the polar; N70, N72–N75
# =====================================================================================================================
def test_drag_polar_V1_parts_sum_best_glide_and_the_u_shaped_drag_curve():  # V1 — C15, N74, N75, D21 (trap T9)
    CD0, AR, es, W, rho, A = I["CD0"], I["AR"], I["e_span"], I["W"], I["rho"], I["A"]
    CL = np.linspace(0.0, 1.4, 15)
    pol = ch14.drag_polar(CL, CD0=CD0, AR=AR, e_span=es)
    assert close(pol, CD0 + ch14.induced_drag_coefficient(CL, AR, es), 1e-14) and close(ch14.drag_polar(CL, CD0=CD0, AR=AR), CD0 + CL ** 2 / (PI * AR))
    bd = ch14.drag_breakdown(CL, CD0=CD0, AR=AR, e_span=es, friction_fraction=0.6)
    assert close(bd["profile"] + bd["induced"], bd["total"], 1e-14) and close(bd["total"], pol, 1e-14)
    assert close(bd["friction"] + bd["form"], bd["profile"], 1e-14) and close(bd["friction"] + bd["form"] + bd["induced"], bd["total"], 1e-14)
    assert close(bd["pressure"], bd["form"] + bd["induced"], 1e-14) and close(bd["friction"] + bd["pressure"], bd["total"], 1e-14)
    assert close(bd["induced_share"], ch14.drag_share(CL, CD0=CD0, AR=AR, e_span=es), 1e-14, 1e-300) and bd["induced_share"][0] == 0.0
    assert np.all(np.diff(bd["induced_share"]) > 0) and "friction" not in ch14.drag_breakdown(0.5, CD0=CD0, AR=AR)
    # the best glide: minimise C_D / C_L numerically
    bg = ch14.best_glide(CD0, AR, es)
    opt = minimize_scalar(lambda c: float(ch14.drag_polar(c, CD0=CD0, AR=AR, e_span=es)) / c, bounds=(0.05, 3.0), method="bounded", options={"xatol": 1e-12})
    assert close(bg["CL"], opt.x, 1e-6) and close(bg["LD_max"], 1 / opt.fun, 1e-10) and close(bg["CD"], 2 * CD0) and bg["induced_share"] == 0.5
    assert close(ch14.drag_share(bg["CL"], CD0=CD0, AR=AR, e_span=es), 0.5, 1e-13) and close(bg["glide_angle"], math.atan(1 / bg["LD_max"]))
    assert close(bg["LD_max"], 0.5 * math.sqrt(PI * AR * es / CD0)) and close(bg["CL"], math.sqrt(PI * AR * es * CD0))
    b5 = ch14.best_glide(0.02, 5.0)
    assert round(b5["CL"], 2) == 0.56 and round(b5["LD_max"]) == 14 and round(bg["LD_max"], 1) == 16.5         # design D21's worked numbers
    # level flight: profile drag rises as U^2, induced drag falls as 1/U^2; the least total where the two are equal
    Us = ch14.speed_min_drag(W, rho, A, AR, CD0, es)
    opt = minimize_scalar(lambda u: float(ch14.drag_vs_speed(W, rho, A, AR, CD0, u, es)["total"]), bounds=(10.0, 120.0), method="bounded", options={"xatol": 1e-10})
    assert close(Us, opt.x, 1e-6)
    d = ch14.drag_vs_speed(W, rho, A, AR, CD0, Us, es)
    assert close(d["profile"], d["induced"], 1e-12) and close(d["induced_share"], 0.5, 1e-12) and close(d["LD"], bg["LD_max"], 1e-12) and close(d["CL"], bg["CL"], 1e-12)
    U = np.array([25.0, 30.0, Us, 55.0, 80.0])
    dv = ch14.drag_vs_speed(W, rho, A, AR, CD0, U, es)
    assert close(dv["profile"] * U ** -2, 0.5 * rho * A * CD0) and close(dv["induced"] * U ** 2, 2 * W ** 2 / (PI * rho * I["s"] ** 2 * es))
    assert close(dv["induced"], ch14.elliptic_induced_drag(W, rho, U, I["s"]) / es, 1e-13)                      # (14.25) / e_span, s^2 = Lambda A
    assert close(dv["CL"], ch14.level_flight(W, rho, U, A)["CL"]) and close(dv["total"], dv["profile"] + dv["induced"]) and close(dv["LD"], W / dv["total"])
    assert close(dv["total"], ch14.level_flight(W, rho, U, A, CD=ch14.drag_polar(dv["CL"], CD0=CD0, AR=AR, e_span=es))["D"], 1e-12)
    assert np.all(np.diff(dv["induced_share"]) < 0)
    # trap T9: "the induced drag is the largest part" holds below the best-glide speed only; at our cruise it is a small share
    slow_, cruise = ch14.drag_vs_speed(W, rho, A, AR, CD0, I["U_slow"], es), ch14.drag_vs_speed(W, rho, A, AR, CD0, I["U_cruise"], es)
    assert slow_["induced_share"] > 0.5 > cruise["induced_share"] and I["U_slow"] < Us < I["U_cruise"]
    assert bench()["induced_share_wide_body_cruise"]["value"] < 0.5        # a cited airliner figure agrees in kind (not a test of our code)
    assert ch14.stall_speed(W, rho, A, I["CL_max"]) < Us                   # the best glide is flyable: its C_L is below the maximum
    assert (Q_(2.0, "N") ** 2 / (Q_(1.2, "kg/m**3") * Q_(1.0, "m/s") ** 2 * Q_(1.0, "m") ** 2)).check("[force]")   # V2: 2 W^2 / (pi rho U^2 s^2) is a force


def test_best_glide_V2_derivation():  # V2 — D21 (two-star, ours): minimum of C_D / C_L and of D(U) in level flight
    CL, CD0, k, W, rho, U, A, s, e = sp.symbols("C_L C_D0 k W rho U A s e", positive=True)      # k = pi Lambda e_span
    ratio = (CD0 + CL ** 2 / k) / CL                                       # C_D / C_L = D / L
    CLs = sp.solve(sp.diff(ratio, CL), CL)
    assert len(CLs) == 1 and sp.simplify(CLs[0] - sp.sqrt(k * CD0)) == 0 and sp.diff(ratio, CL, 2).subs(CL, CLs[0]) > 0
    assert sp.simplify((CLs[0] ** 2 / k) - CD0) == 0                       # induced = profile at the best glide
    assert sp.simplify(1 / ratio.subs(CL, CLs[0]) - sp.sqrt(k / CD0) / 2) == 0            # (L/D)_max = (1/2) sqrt(pi Lambda e / C_D0)
    D = rho * U ** 2 * A * CD0 / 2 + 2 * W ** 2 / (sp.pi * rho * U ** 2 * s ** 2 * e)     # level flight: L = W
    Us = sp.solve(sp.diff(D, U), U)
    Us = [u for u in Us if u.is_positive]
    assert len(Us) == 1
    assert sp.simplify((rho * Us[0] ** 2 * A * CD0 / 2) - 2 * W ** 2 / (sp.pi * rho * Us[0] ** 2 * s ** 2 * e)) == 0   # the two parts are equal there
    CL_at = W / (rho * Us[0] ** 2 * A / 2)
    assert sp.simplify(CL_at ** 2 - sp.pi * (s ** 2 / A) * e * CD0) == 0   # ... and C_L is the best-glide value with k = pi (s^2/A) e
    # the induced part of D is (14.25) divided by e, i.e. C_Di = C_L^2 / (pi Lambda e) with L = W
    q = rho * U ** 2 * A / 2
    assert sp.simplify(2 * W ** 2 / (sp.pi * rho * U ** 2 * s ** 2 * e) - (W / q) ** 2 / (sp.pi * (s ** 2 / A) * e) * q) == 0


def test_stall_models_V7_generic_curves_and_regimes_are_qualitative_and_say_so():  # V7 — N70, N72, N32 (label: qualitative; trap T16)
    kw = dict(K=I["K"], beta=I["beta"], alpha_stall=I["alpha_stall"], CL_post=I["CL_post"], width=I["stall_width"])
    for f in (ch14.generic_lift_curve, ch14.generic_drag_curve, ch14.stall_regime, ch14.suction_side_separation, ch14.wake_rollup):
        assert "qualitative" in f.__doc__.lower(), f.__name__
    assert "approximate" in ch14.profile_drag_estimate.__doc__.lower() and "model" in ch14.impulsive_start_lumped_vortex.__doc__.lower()
    al = np.linspace(-0.15, kw["alpha_stall"], 40)
    assert close(ch14.generic_lift_curve(al, **kw), kw["K"] * (al + kw["beta"]), 1e-15, 0)                     # the line (14.12) below stall, exactly
    assert close(ch14.generic_lift_curve(al, **kw), ch14.zhukhovsky_lift_coefficient(al, kw["beta"]), 1e-13, 1e-15)
    post = np.linspace(kw["alpha_stall"] + kw["width"], 0.5, 9)
    assert close(ch14.generic_lift_curve(post, **kw), kw["CL_post"]) and ch14.generic_lift_curve(-kw["beta"], **kw) == 0.0
    a = np.linspace(kw["alpha_stall"] - 0.01, kw["alpha_stall"] + kw["width"] + 0.01, 20001)
    cl = ch14.generic_lift_curve(a, **kw)
    slope = np.gradient(cl, a)
    assert np.max(np.abs(np.diff(cl))) < 60 * (a[1] - a[0]) and np.max(np.abs(np.diff(slope))) < 0.05 * kw["K"]   # continuous, with a continuous slope
    assert close(slope[5], kw["K"], 1e-6) and abs(slope[-5]) < 1e-6
    k = int(np.argmax(cl))
    assert kw["alpha_stall"] < a[k] < kw["alpha_stall"] + kw["width"] and cl[k] > kw["K"] * (kw["alpha_stall"] + kw["beta"])   # a maximum, then a drop
    assert cl[k] > kw["CL_post"]
    dk = dict(CD0=I["CD0"], AR=I["AR"], e_span=I["e_span"], CD_sep=I["CD_sep"], **kw)
    aa = np.linspace(-0.2, 0.45, 2601)
    cd = ch14.generic_drag_curve(aa, **dk)
    assert ch14.generic_drag_curve(-kw["beta"], **dk) == I["CD0"] and cd.min() >= I["CD0"] and close(cd.min(), I["CD0"], 1e-6)
    assert abs(aa[np.argmin(cd)] + kw["beta"]) < 3e-4                      # least drag at zero lift
    below = aa <= kw["alpha_stall"]
    assert close(cd[below], ch14.drag_polar(kw["K"] * (aa[below] + kw["beta"]), CD0=I["CD0"], AR=I["AR"], e_span=I["e_span"]), 1e-13)   # the polar below stall
    assert close(ch14.generic_drag_curve(0.45, **dk), I["CD0"] + I["CL_post"] ** 2 / (PI * I["AR"] * I["e_span"]) + I["CD_sep"])
    assert np.all(np.diff(cd[aa > -kw["beta"]][: np.sum(below & (aa > -kw["beta"]))]) > 0)
    with pytest.raises(TypeError):
        ch14.generic_lift_curve(0.1)                                       # every constant is a required keyword: no book number built in
    # the two-threshold classifier
    lam, tur = I["Re_laminar"], I["Re_turbulent"]
    assert ch14.stall_regime(0.99 * lam, Re_laminar=lam, Re_turbulent=tur) == "laminar separation first"
    assert ch14.stall_regime(1.01 * lam, Re_laminar=lam, Re_turbulent=tur) == "mixed" == ch14.stall_regime(0.99 * tur, Re_laminar=lam, Re_turbulent=tur)
    assert ch14.stall_regime(1.01 * tur, Re_laminar=lam, Re_turbulent=tur) == "transition first"
    with pytest.raises(TypeError):
        ch14.stall_regime(1e6)
    Re = ch14.chord_reynolds(I["U_cruise"], I["c_mean"], I["nu"])
    assert close(Re, I["U_cruise"] * I["c_mean"] / I["nu"]) and ch14.stall_regime(Re, Re_laminar=lam, Re_turbulent=tur) == "transition first"


def test_profile_drag_estimate_V7_twice_the_plate_friction_times_a_stated_form_factor():  # V7 — C15 (label: approximate)
    from fluidpy.core.boundary_layer import plate_drag_coefficient

    for Re, regime in ((2e5, "laminar"), (5e6, "turbulent"), (2e6, "mixed")):
        base = 2 * float(plate_drag_coefficient(Re, regime=regime, sides=1))
        assert close(ch14.profile_drag_estimate(Re, 0.12, regime=regime, form_factor=1.0), base, 1e-14)        # both faces of a flat plate
        assert close(ch14.profile_drag_estimate(Re, 0.12, regime=regime, form_factor=1.4), 1.4 * base, 1e-14)
        assert close(ch14.profile_drag_estimate(Re, 0.12, regime=regime, form_factor=lambda t: 1 + 2 * t), 1.24 * base, 1e-14)
    lam = [ch14.profile_drag_estimate(Re, 0.1, regime="laminar", form_factor=1.0) for Re in (1e5, 4e5)]
    assert close(lam[0] / lam[1], 2.0, 1e-12) and close(lam[0] * math.sqrt(1e5) / 2, 1.328, 1e-3)              # (9.33): 1.328 / sqrt(Re) per face
    assert ch14.profile_drag_estimate(5e6, 0.1, regime="turbulent", form_factor=1.0) > ch14.profile_drag_estimate(5e6, 0.1, regime="laminar", form_factor=1.0)
    with pytest.raises(TypeError):
        ch14.profile_drag_estimate(1e6, 0.12, regime="laminar")            # the form factor is an input, never a remembered correlation


# =====================================================================================================================
# 17. C16 — thrust from lift: fish tail, bird stroke, sail and keel; N77, N79, N81, N82
# =====================================================================================================================
def _sail_vectors(U, V, theta, CL, CD, rho, A):
    """Our own vector construction: apparent wind U_r = U - V, drag along it, lift normal to it (toward the bow)."""
    Uvec = -U * np.array([math.cos(theta), math.sin(theta)])               # the true wind blows FROM the direction theta off the bow
    Ur = Uvec - np.array([V, 0.0])                                         # the boat moves along +x
    q = 0.5 * rho * (Ur @ Ur) * A
    d = Ur / np.linalg.norm(Ur)
    n = np.array([-d[1], d[0]]) if (-d[1]) > 0 else np.array([d[1], -d[0]])     # the normal with a forward component
    F = q * CD * d + q * CL * n
    return Ur, F, q


def test_apparent_wind_V1_vector_subtraction_and_sail_forces_resolved_along_the_course():  # V1 — C16 (ours), D22
    U, V, rho, A, CL, CD = I["U_wind"], I["V_boat"], I["rho"], I["A_sail"], I["CL_sail"], I["CD_sail"]
    for th in np.deg2rad([20.0, 45.0, 90.0, 135.0, 170.0]):
        Ur, F, q = _sail_vectors(U, V, th, CL, CD, rho, A)
        aw = ch14.apparent_wind(U, V, th)
        assert close(aw["vector"], Ur, 1e-13, 1e-13) and close(aw["Ur"], np.linalg.norm(Ur)) and close(aw["Ur"], math.sqrt(U ** 2 + V ** 2 + 2 * U * V * math.cos(th)))
        assert close(aw["beta_a"], math.atan2(U * math.sin(th), U * math.cos(th) + V)) and 0 < aw["beta_a"] < th   # swung toward the bow, and faster
        sf = ch14.sail_forces(U, V, th, CL=CL, CD=CD, rho=rho, A=A)
        assert close(sf["L"], q * CL) and close(sf["D"], q * CD) and close(sf["drive"], F[0], 1e-12, 1e-10) and close(sf["lateral"], -F[1], 1e-12)   # leeward +
        assert close(sf["drive"], sf["L"] * math.sin(aw["beta_a"]) - sf["D"] * math.cos(aw["beta_a"]), 1e-13, 1e-12)
        assert close(sf["drive"] ** 2 + sf["lateral"] ** 2, sf["L"] ** 2 + sf["D"] ** 2, 1e-13) and close(sf["F"], math.hypot(sf["L"], sf["D"]))
        assert close(sf["drive_over_L"], sf["drive"] / sf["L"]) and close(sf["beta_a"], aw["beta_a"]) and close(sf["Ur"], aw["Ur"])
        m = ch14.sail_forces(U, V, -th, CL=CL, CD=CD, rho=rho, A=A)        # the other tack: mirror image
        assert close(m["drive"], sf["drive"], 1e-13, 1e-12) and close(m["lateral"], sf["lateral"]) and close(m["beta_a"], -sf["beta_a"])
        assert close(ch14.apparent_wind(U, V, -th)["Ur"], aw["Ur"]) and close(ch14.velocity_made_good(V, th), V * math.cos(th), 1e-14, 1e-15)
    at_rest = ch14.apparent_wind(U, 0.0, 0.8)
    assert close(at_rest["Ur"], U) and close(at_rest["beta_a"], 0.8)       # V = 0: apparent = true wind
    run = ch14.apparent_wind(U, V, PI)
    assert close(run["Ur"], U - V) and close(abs(run["beta_a"]), PI)       # running before the wind: the boat's speed is subtracted
    # the drive vanishes exactly where tan(beta_a) = D / L, and is negative closer to the wind
    th0 = brentq(lambda t: ch14.apparent_wind(U, V, t)["beta_a"] - math.atan2(CD, CL), 0.01, 1.5, xtol=1e-15)
    assert abs(ch14.sail_forces(U, V, th0, CL=CL, CD=CD, rho=rho, A=A)["drive"]) < 1e-10 * 0.5 * rho * U ** 2 * A
    assert ch14.sail_forces(U, V, 0.8 * th0, CL=CL, CD=CD, rho=rho, A=A)["drive"] < 0 < ch14.sail_forces(U, V, 1.2 * th0, CL=CL, CD=CD, rho=rho, A=A)["drive"]
    ths = np.linspace(th0, np.deg2rad(60.0), 50)
    curve = ch14.sail_drive_curve(ths, U=U, V=V, CL=CL, CD=CD, rho=rho, A=A)
    assert close(curve, [ch14.sail_forces(U, V, t, CL=CL, CD=CD, rho=rho, A=A)["drive"] for t in ths], 1e-13, 1e-10)
    assert np.all(np.diff(curve) > 0)                                      # "the thrust decreases as the angle decreases and vanishes"
    assert close(ch14.sail_forces(U, V, 0.7, CL=CL, CD=0.0, rho=rho, A=A)["drive_over_L"], math.sin(ch14.apparent_wind(U, V, 0.7)["beta_a"]))


def test_sailboat_V1_force_polygon_closes_and_the_course_limit_is_the_sum_of_two_drag_angles():  # V1 + V4 — N81, N82 (ours)
    kw = dict(CL=I["CL_sail"], CD=I["CD_sail"], rho_air=I["rho"], A_sail=I["A_sail"], rho_water=I["rho_water"], CDA_hull=I["CDA_hull"])
    lim = ch14.sailboat_course_limit(I["CL_sail"] / I["CD_sail"], I["LD_hull"])
    assert close(lim, math.atan(I["CD_sail"] / I["CL_sail"]) + math.atan(1 / I["LD_hull"]))
    assert close(ch14.sailboat_course_limit(np.array([6.0, 1e12]), np.array([8.0, 1e12])), [math.atan(1 / 6) + math.atan(1 / 8), 0.0], 1e-12, 1e-11)
    for th in np.deg2rad([25.0, 35.0, 45.0, 70.0, 110.0, 150.0]):
        r = ch14.sailboat_equilibrium(I["U_wind"], th, LD_hull=I["LD_hull"], **kw)
        assert r["sails"] and r["V"] > 0
        hd = 0.5 * I["rho_water"] * r["V"] ** 2 * I["CDA_hull"]
        sf = ch14.sail_forces(I["U_wind"], r["V"], th, CL=kw["CL"], CD=kw["CD"], rho=kw["rho_air"], A=kw["A_sail"])
        assert close(r["hull_drag"], hd) and close(r["drive"], hd, 1e-10) and close(r["drive"], sf["drive"]) and close(r["keel_lift"], r["lateral"])
        _, F, _ = _sail_vectors(I["U_wind"], r["V"], th, kw["CL"], kw["CD"], kw["rho_air"], kw["A_sail"])
        hull = np.array([-hd, r["keel_lift"]])                             # hull drag against the motion, keel lift to windward
        assert np.max(np.abs(F + hull)) < 1e-9 * np.linalg.norm(F)         # the two forces close: steady sailing
        assert abs(r["closure"]) < 1e-10 and close(r["beta_a"], r["eps_A"] + r["eps_H"], 1e-10)                # beta_a = eps_A + eps_H
        assert close(r["eps_A"], math.atan2(kw["CD"], kw["CL"])) and close(math.tan(r["eps_H"]), 1 / r["LD_hull_required"], 1e-10)
        assert close(r["vmg"], r["V"] * math.cos(th)) and close(r["LD_hull_required"], r["lateral"] / hd, 1e-10)
        assert r["keel_holds"] == (r["beta_a"] >= lim - 1e-12)             # the keel can hold the course only outside the limit
        V2 = brentq(lambda v: _sail_vectors(I["U_wind"], v, th, kw["CL"], kw["CD"], kw["rho_air"], kw["A_sail"])[1][0]
                    - 0.5 * I["rho_water"] * v ** 2 * I["CDA_hull"], 1e-6, 50.0, xtol=1e-13)
        assert close(r["V"], V2, 1e-9)                                     # our own root of drive = hull drag
    irons = ch14.sailboat_equilibrium(I["U_wind"], 0.5 * math.atan2(kw["CD"], kw["CL"]), **kw)                 # pointing too high: no drive at rest
    assert irons["sails"] is False and irons["V"] == 0.0 and irons["vmg"] == 0.0 and irons["drive"] < 0
    free = ch14.sailboat_equilibrium(I["U_wind"], 0.6, **kw)               # no keel limit given: not checked
    assert free["keel_holds"] is True
    d45 = ch14.sailboat_equilibrium(I["U_wind"], I["theta"], LD_hull=I["LD_hull"], **kw)
    assert d45["beta_a"] < I["theta"] and d45["beta_a"] > lim and d45["keel_holds"]
    vmg = [ch14.sailboat_equilibrium(I["U_wind"], t, **kw)["vmg"] for t in np.deg2rad(np.arange(15.0, 90.0, 2.5))]
    k = int(np.argmax(vmg))
    assert 0 < k < len(vmg) - 1 and vmg[k] > 0                             # tacking: a best course angle strictly between the limit and a beam reach


def test_oscillating_foil_V1_thrust_from_a_sideways_stroke_and_the_power_identity():  # V1 + V7 — N77, N79 (trap T18)
    U, V, L, D = I["U_fish"], I["V_tail"], I["L_tail"], I["D_tail"]
    r0 = ch14.oscillating_foil(U, V, L)
    th = math.atan2(V, U)
    assert close(r0["theta"], th) and close(r0["Ur"], math.hypot(U, V)) and close(r0["thrust"], L * math.sin(th)) and close(r0["side_force"], L * math.cos(th))
    assert close(r0["power_out"], U * L * math.sin(th)) and close(r0["power_in"], V * L * math.cos(th)) and close(r0["efficiency"], 1.0, 1e-14)
    for Vv in (0.1, 0.5, 2.0, 7.0):                                        # the identity of the triangle: any stroke speed
        assert close(ch14.oscillating_foil(U, Vv, L)["efficiency"], 1.0, 1e-13)
    # vectors: water relative to the tail (-U, -V); drag along it, lift normal to it with a forward component
    rel = np.array([-U, -V]) / math.hypot(U, V)
    F = D * rel + L * np.array([-rel[1], rel[0]]) * (1 if -rel[1] > 0 else -1)
    r = ch14.oscillating_foil(U, V, L, D)
    assert close(r["thrust"], F[0], 1e-13) and close(r["side_force"], -F[1], 1e-13)                            # the side force opposes the stroke
    assert close(r["efficiency"], (F[0] * U) / (-F[1] * V), 1e-13) and 0 < r["efficiency"] < 1
    assert close(r["thrust"], L * math.sin(th) - D * math.cos(th)) and close(r["side_force"], L * math.cos(th) + D * math.sin(th))
    stop = ch14.oscillating_foil(U, U * D / L, L, D)                       # tan(theta) = D / L: no thrust left
    assert abs(stop["thrust"]) < 1e-14 * L and abs(stop["efficiency"]) < 1e-13
    effs = ch14.oscillating_foil(U, np.array([0.2, 0.75, 3.0]), L, D)["efficiency"]
    assert effs.shape == (3,) and np.all(effs < 1) and effs[1] > effs[0]
    # a bird's two strokes (N79): a fast downstroke gives thrust, a slow feathered upstroke (little lift) costs little
    down, up = ch14.oscillating_foil(10.0, 4.0, 6.0, 0.5), ch14.oscillating_foil(10.0, 1.0, 0.5, 0.3)
    assert down["thrust"] > 0 and 0.5 * (down["thrust"] + up["thrust"]) > 0 and abs(up["thrust"]) < 0.3 * down["thrust"]


def test_sailing_V2_derivation():  # V2 — D22 (two-star, ours): U_r = U - V, drive = L sin(beta_a) - D cos(beta_a), beta_a = eps_A + eps_H
    U, V, L, D, Lh, Dh = sp.symbols("U V L D L_h D_h", positive=True)
    th, ba = sp.symbols("theta beta_a", positive=True)
    Uvec = sp.Matrix([-U * sp.cos(th), -U * sp.sin(th)])                   # wind FROM the angle theta off the bow
    Ur = Uvec - sp.Matrix([V, 0])
    assert sp.simplify(Ur.dot(Ur) - (U ** 2 + V ** 2 + 2 * U * V * sp.cos(th))) == 0
    tan_ba = (-Ur[1]) / (-Ur[0])                                           # the apparent wind comes from (cos beta_a, sin beta_a)
    assert sp.simplify(tan_ba - U * sp.sin(th) / (U * sp.cos(th) + V)) == 0
    assert sp.simplify(tan_ba.subs(V, 0) - sp.tan(th)) == 0                # V = 0: beta_a = theta
    assert sp.simplify(sp.diff(tan_ba, V)) != 0 and sp.simplify(sp.diff(tan_ba, V) + U * sp.sin(th) / (U * sp.cos(th) + V) ** 2) == 0   # falls with V
    d = sp.Matrix([-sp.cos(ba), -sp.sin(ba)])                              # unit vector along the apparent wind
    n = sp.Matrix([sp.sin(ba), -sp.cos(ba)])                               # normal to it, toward the bow
    assert d.dot(n) == 0 and sp.simplify(n.dot(n) - 1) == 0
    F = D * d + L * n
    drive, lateral = F[0], -F[1]
    assert sp.simplify(drive - (L * sp.sin(ba) - D * sp.cos(ba))) == 0 and sp.simplify(lateral - (L * sp.cos(ba) + D * sp.sin(ba))) == 0
    assert sp.simplify(drive.subs(ba, sp.atan(D / L))) == 0                # zero at tan(beta_a) = D / L
    assert sp.simplify(drive ** 2 + lateral ** 2 - L ** 2 - D ** 2) == 0
    # steady sailing: hull drag = drive, keel lift = lateral; then tan(eps_H) = D_h / L_h = tan(beta_a - eps_A)
    epsA = sp.atan(D / L)
    lhs = sp.tan(ba - epsA)
    assert sp.simplify(sp.expand_trig(lhs) - drive / lateral) == 0         # so beta_a = eps_A + eps_H
    # the fish tail: U sin(theta) = V cos(theta) when V/U = tan(theta): the two powers are equal
    t2 = sp.atan(V / U)
    assert sp.simplify(U * L * sp.sin(t2) - V * L * sp.cos(t2)) == 0
    eta = ((L * sp.sin(t2) - D * sp.cos(t2)) * U) / ((L * sp.cos(t2) + D * sp.sin(t2)) * V)
    assert sp.simplify(eta.subs(D, 0) - 1) == 0 and sp.simplify(eta.subs(V, U * D / L)) == 0


# =====================================================================================================================
# 18. The design's worked numbers (analysis/ch14_design.md, Part F "Check" lines and C.6) — all computed from OUR inputs
# =====================================================================================================================
def test_design_numbers_V7_worked_numbers_of_part_f_are_reproduced():  # V7 — the notebook will print these; they must be right
    assert ch14.aspect_ratio(10.0, 20.0) == 5.0 and close(ch14.level_flight(1e4, 1.25, 40.0, 20.0)["CL"], 0.5)          # D01
    assert round(ch14.flat_plate_forces(1.0, 10.0, 1.0, np.deg2rad(30.0))["Gamma_cw"], 1) == 15.7                       # D02
    mu, a = ch14.zhukhovsky_circle(B, E, BETA)
    GK = ch14.kutta_circulation(US, a, ALPHA, BETA)
    exact, rule = ch14.zhukhovsky_lift_coefficient(ALPHA, BETA, E, exact=True), ch14.zhukhovsky_lift_coefficient(ALPHA, BETA)
    assert round(GK, 2) == 15.96 and round(exact, 3) == 1.058 and round(rule, 3) == 0.987                               # D09
    f = ch14.zhukhovsky_force(*SEC, RHO)
    assert close(f["CL"], exact, 1e-9) and close(f["L"], RHO * US * GK, 1e-10)
    assert close(I["AR"], 9.6) and close(I["c_mean"], 1.25)                                                             # the aircraft of convention 11
    lf_slow, lf_cruise = ch14.level_flight(I["W"], I["rho"], I["U_slow"], I["A"]), ch14.level_flight(I["W"], I["rho"], I["U_cruise"], I["A"])
    assert lf_slow["CL"] < I["CL_max"] and lf_cruise["CL"] < lf_slow["CL"] and close(lf_slow["CL"] / lf_cruise["CL"], (55.0 / 30.0) ** 2)
    assert ch14.compressibility_check(I["U_cruise"])["incompressible_ok"] and not ch14.compressibility_check(250.0)["incompressible_ok"]   # N01
    demo = ch14.solve_lifting_line(I["wd_c"], I["wd_s"], I["wd_U"], I["wd_alpha"], K=I["wd_K"], N=40, symmetric=True, rho=I["wd_rho"])
    assert round(demo["e_span"], 3) == 0.954 and round(demo["CL"], 3) == 0.395                                          # the demo wing, Lambda = 6
    b = ch14.sailboat_equilibrium(I["U_wind"], I["theta"], CL=I["CL_sail"], CD=I["CD_sail"], rho_air=I["rho"], A_sail=I["A_sail"],
                                  rho_water=I["rho_water"], CDA_hull=I["CDA_hull"], LD_hull=I["LD_hull"])
    assert b["sails"] and b["keel_holds"] and 0 < b["V"] < I["U_wind"]
    fish = ch14.oscillating_foil(I["U_fish"], I["V_tail"], I["L_tail"], I["D_tail"])
    assert close(fish["theta"], math.atan(0.75)) and close(fish["thrust"], 10 * 0.6 - 1 * 0.8) and close(fish["side_force"], 10 * 0.8 + 1 * 0.6)   # a 3-4-5 triangle


# =====================================================================================================================
# 19. Figure helpers and scripts (design C.4)
# =====================================================================================================================
def test_drawings_V7_every_sketch_and_composite_returns_a_figure():  # V7 — design C.4: scripts/ch14_drawings.py (18 helpers)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure

    sys.path.insert(0, str(ROOT / "scripts"))
    try:
        import ch14_drawings as DR
    finally:
        sys.path.remove(str(ROOT / "scripts"))
    assert set(DRAWINGS) <= set(DR.__all__)
    S = ch14.zhukhovsky_surface(*SEC, 600)
    A = ch14.zhukhovsky_airfoil(B, E, BETA, 400)
    special = {"fig_three_stages": lambda: DR.fig_three_stages(*SEC), "fig_pressure_two_ways": lambda: DR.fig_pressure_two_ways(S),
               "draw_section_labelled": lambda: DR.draw_section_labelled(A["x"], A["y"], ch14.section_geometry(A["x"], A["y"], te_index=0)),
               "draw_horseshoe_staircase": lambda: DR.draw_horseshoe_staircase(5)}
    for name in DRAWINGS:
        fig = special[name]() if name in special else getattr(DR, name)()
        assert isinstance(fig, Figure) and len(fig.axes) >= 1 and any(len(ax.lines) + len(ax.patches) + len(ax.collections) > 0 for ax in fig.axes), name
        plt.close(fig)
    for kind in ("fish", "bird_down", "bird_up", "sail"):
        fig = DR.draw_velocity_triangle(kind)
        assert isinstance(fig, Figure)
        plt.close(fig)
    fig, ax = plt.subplots()
    assert DR.draw_span_semicircle(ax=ax) is fig                           # drawing into a given Axes returns its Figure
    plt.close("all")


@slow
def test_scripts_V7_every_chapter_script_runs_headless(tmp_path):  # V7 — the nine demos + the drawings + make_refs exit 0
    ran = 0
    for s in sorted((ROOT / "scripts").glob("ch14_*.py")):
        if s.name == "ch14_common.py":
            continue
        r = subprocess.run([sys.executable, str(s), "--no-show", "--fast", "--out", str(tmp_path)], cwd=str(ROOT), capture_output=True,
                           text=True, encoding="utf-8", errors="replace", timeout=900)
        assert r.returncode == 0, (s.name, r.stderr[-2000:])
        ran += 1
    assert ran == 9 and len(list(tmp_path.glob("*.png"))) >= 9
    before = (REF / "benchmarks.json").read_bytes()
    r = subprocess.run([sys.executable, str(REF / "make_refs.py")], cwd=str(ROOT), capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=120)
    assert r.returncode == 0 and (REF / "benchmarks.json").read_bytes() == before      # idempotent


def test_reference_V5_sources_file_cites_every_entry_and_no_book_table():  # V5 — data-and-benchmarks policy
    src = (REF / "SOURCES.md").read_text(encoding="utf-8")
    bm = bench()
    for key, row in bm.items():
        if key.startswith("_"):
            continue
        assert key in src and row["verified"] == "2026-10-12" and row["url"].startswith("https://") and row["url"] in src, key
        assert "book" not in row["source"].lower() and "kundu" not in row["source"].lower()
    assert "arXiv:2104.15122" in src and "flow.byu.edu" in src and "not read" in src and "Read and deliberately not used" in src
    assert len([k for k in bm if not k.startswith("_")]) == 8
    doc = ch14.wagner_function.__doc__
    assert "0.165" in doc and "0.0455" in doc and "0.335" in doc and "semi-chords" in doc     # the docstring states the same fit


# =====================================================================================================================
# 99. Private book values (V6) — read from tests/book_values_ch14.json; skipped when the file is absent
# =====================================================================================================================
@book_only
def test_book_V6_section_shape_rules_and_lift_slope():  # V6 — (14.8), slips #4 and #8, §14.5, §14.7
    bv = book()
    conf = bv["sec_14_4_conformal"]
    assert rel(ch14.THICKNESS_RULE_COEFFICIENT, conf["eq_14_8_thickness_coefficient"]) < 1e-3          # the printed coefficient is ours rounded
    d = conf["derived_by_analyst"]
    assert rel(ch14.THICKNESS_RULE_COEFFICIENT, d["thickness_coefficient_exact_first_order"]) < 1e-14
    for e, key in ((0.02, "t_over_c_over_e_at_e_0.02"), (0.05, "t_over_c_over_e_at_e_0.05"), (0.10, "t_over_c_over_e_at_e_0.10")):
        assert rel(ch14.zhukhovsky_shape(1.0, e, 0.0)["thickness_ratio"] / e, d[key]) < 5e-5, key
    assert rel(ch14.zhukhovsky_shape(1.0, 0.1, 0.0)["chord"] / 4, d["chord_over_4b_at_e_0.10"]) < 1e-6
    lift = bv["sec_14_5_lift"]["derived_by_analyst"]
    assert rel(ch14.zhukhovsky_lift_coefficient(0.1, 0.0, 0.1, exact=True) / (2 * PI * 0.1), lift["CL_exact_over_2pi_alpha_at_e_0.1_alpha_0.1_beta_0"]) < 5e-5
    ch = bv["sec_14_7_characteristics"]
    K_deg = ch["K_per_degree_approx"]
    assert rel(K_deg * 180 / PI, ch["derived_by_analyst"]["K_per_rad_from_0.1_per_deg"]) < 1e-12        # per degree -> per radian (trap T12)
    assert abs(K_deg * 180 / PI / (2 * PI) - 1) < 0.10                                                  # the quoted slope is within 10 % of 2 pi
    assert abs(bv["sec_14_6_lifting_line"]["K_most_airfoils_approx"] / (2 * PI) - 1) < 0.05
    lim = bv["sec_14_1_introduction"]["incompressible_mach_limit"]
    c = 340.29
    assert ch14.compressibility_check(0.99 * lim * c)["incompressible_ok"] and not ch14.compressibility_check(1.01 * lim * c)["incompressible_ok"]
    # rule 9: our illustrative inputs are visibly not the book's exercise inputs
    ex = bv["exercises"]
    assert I["W"] != ex["ex_14_13"]["weight_N"] and I["s"] != ex["ex_14_13"]["span_m"] and I["K"] != ex["ex_14_14"]["K_per_rad"]
    assert I["Gamma_v"] != ex["ex_14_16"]["Gamma_prime"] and I["alpha_stall"] != np.deg2rad(ch["alpha_max_lift_deg_approx"])
    assert not (ch["post_stall_CL"][0] <= I["CL_post"] <= ch["post_stall_CL"][1])
    # slip #5: the printed birth year is not the cited one
    assert ch14.book_slips()["5"]["correct"] == "1868" and not bv["sec_14_3_airfoil_sections"]["historical"]["lanchester_dates_printed"].startswith("1868")


@book_only
def test_book_V6_exercises_on_sections():  # V6 — exercises with stated answers (section part)
    ex = book()["exercises"]
    e4 = ex["ex_14_4"]
    xb = e4["break_x_over_c"]
    xs = np.array([0.0, xb, 1.0])
    Cy = ch14.normal_force_coefficient(xs, [e4["Cp_upper_const"], e4["Cp_upper_const"], e4["Cp_trailing_edge"]],
                                       xs, [e4["Cp_lower_const"], e4["Cp_lower_const"], e4["Cp_trailing_edge"]])
    CL = Cy * math.cos(np.deg2rad(e4["alpha_deg"]))                        # pressure only: the normal force projected on the lift direction
    assert rel(Cy, e4["derived_by_analyst"]) < 1e-12 and rel(CL, e4["answer_CL"]) < 0.005
    assert rel(ch14.lift_from_cp(xs, [e4["Cp_upper_const"]] * 2 + [e4["Cp_trailing_edge"]], xs, [e4["Cp_lower_const"]] * 2 + [e4["Cp_trailing_edge"]],
                                 np.deg2rad(e4["alpha_deg"]), suction="none"), e4["answer_CL"]) < 0.005
    U, b = 1.0, 1.0                                                         # the flat plate: velocities at a named point in both planes
    for al in (0.05, 0.2):
        F = ch14.flat_plate_flow(U, b, al)
        qz = complex(F.dwdzeta(-2 * b + 0j))
        assert close([qz.real, -qz.imag], [0.75 * U * math.cos(al), 2.25 * U * math.sin(al)], 1e-12)
        u, v = F.velocity(-2.5 * b, 0.0)
        assert close([u, v], [U * math.cos(al), 3 * U * math.sin(al)], 1e-12) and close(complex(ch14.joukowski(-2 * b, b)), -2.5 * b)
        assert close(F.Gamma_cw, 4 * PI * U * b * math.sin(al))
    e8 = ex["ex_14_8"]
    S = ch14.zhukhovsky_airfoil(e8["b"], beta=e8["beta"], n=400, a=e8["a"])
    assert abs(S["z"][0] - 2 * e8["b"]) < 1e-13 and close(S["e"], e8["a"] * math.cos(e8["beta"]) / e8["b"] - 1) and S["e"] > 0
    F8 = ch14.zhukhovsky_flow(1.0, e8["b"], 0.0, e8["beta"], np.deg2rad(e8["alpha_deg"]), a=e8["a"])
    assert close(F8.Gamma_cw, 4 * PI * e8["a"] * math.sin(np.deg2rad(e8["alpha_deg"]) + e8["beta"])) and F8.a == e8["a"]
    assert np.ptp(np.asarray(F8.w_zeta(F8.mu + F8.a * np.exp(1j * np.linspace(-3, 3, 200)))).imag) < 1e-12
    e9 = ex["ex_14_9"]
    cl = e9["CL_at_zero_incidence"] + 2 * PI * np.deg2rad(e9["alpha_deg"])
    assert rel(cl, e9["derived_by_analyst_CL"]) < 1e-4
    beta9 = e9["CL_at_zero_incidence"] / (2 * PI)
    assert rel(ch14.zhukhovsky_lift_coefficient(np.deg2rad(e9["alpha_deg"]), beta9), e9["derived_by_analyst_CL"]) < 1e-4
    for n in ex["ex_14_2_n_values"]:
        wd = ch14.trailing_edge_wedge(n)
        assert close(wd["included_angle"], 2 * PI * (1 - 1 / n)) and wd["Cp_edge"] == 1.0 and 0 < wd["included_angle"] < PI


@book_only
def test_book_V6_exercises_on_wings():  # V6 — exercises with checkable results (wing part)
    ex = book()["exercises"]
    e = ex["ex_14_13"]
    d = e["derived_by_analyst"]
    U = e["speed_km_h"] / 3.6
    t = ch14.tip_vortex_estimate(e["weight_N"], d["rho_assumed"], U, e["span_m"])
    assert rel(ch14.elliptic_induced_drag(e["weight_N"], d["rho_assumed"], U, e["span_m"]), d["Di_N"]) < 1e-4 and rel(t["induced_drag"], d["Di_N"]) < 1e-4
    assert rel(t["Gamma0"], d["Gamma1_m2_s"]) < 2e-4
    assert rel(ch14.elliptic_circulation(e["span_m"] / 4, e["span_m"], t["Gamma0"]), d["Gamma_halfway_m2_s"]) < 2e-4
    e = ex["ex_14_14"]
    d = e["derived_by_analyst"]
    ar = np.deg2rad(e["alpha_r_deg"])
    tw = ch14.elliptic_twist_rectangular(ar, e["c_m"], e["span_m_printed_as_b"], K=e["K_per_rad"])
    G1 = tw["Gamma1_over_U"] * e["U_m_s"]
    assert rel(np.rad2deg(tw["alpha_t"]), d["alpha_t_deg"]) < 1e-4 and rel(G1, d["Gamma1"]) < 1e-5
    assert rel(ch14.lift_from_series([G1], e["span_m_printed_as_b"], e["rho"], e["U_m_s"]), d["L_N"]) < 1e-5
    assert rel(ch14.induced_drag_from_series([G1], e["rho"]), d["Di_N"]) < 1e-4
    r = ch14.solve_lifting_line(e["c_m"], e["span_m_printed_as_b"], e["U_m_s"], tw["alpha_t"], K=e["K_per_rad"], rho=e["rho"], N=12,
                                twist=lambda y: ar * np.sqrt(np.clip(1 - (2 * y / e["span_m_printed_as_b"]) ** 2, 0, None)))
    assert rel(r["L"], d["L_N"]) < 1e-5 and rel(r["Di"], d["Di_N"]) < 1e-4 and np.max(np.abs(r["Gn"][1:])) < 1e-10 * r["Gn"][0]
    e = ex["ex_14_16"]
    ell = e["s"] / (2 * e["s_over_2ell"])
    m = ch14.wake_vortex_rolling_moment(e["rho"], e["U"], e["K"], e["s"], e["c"], e["Gamma_prime"], ell)
    assert rel(m["M"], e["derived_by_analyst"]["M_N_m"]) < 1e-5
    e = ex["ex_14_19"]
    E1 = ch14.vortex_pair_energy_per_length(1.0, e["span_m"], e["core_m"], 1.0)
    E2 = ch14.vortex_pair_energy_per_length(1.0, e["span_m"], 2 * e["core_m"], 1.0)
    assert abs(100 * (E2 / E1 - 1) - e["derived_by_analyst"]["change_percent_approx"]) < 1.5            # a doubled core: about that much less drag
    hw = ch14.horseshoe_wing(2.0, 3.0, 1.2, 10.0)                          # Exercise-type reference results: L = rho U Gamma s, D_i = rho Gamma^2 / pi
    assert close(hw["L"], 1.2 * 10.0 * 2.0 * 3.0) and close(hw["Di_centre"], 1.2 * 4.0 / PI)
    w0 = ch14.downwash(0.0, lambda y: 1.0 - (2 * y / 3.0) ** 2, 3.0)       # the parabolic loading: 4/pi times the elliptic mid-span downwash
    assert close(w0 / ch14.elliptic_downwash(1.0, 3.0), 4 / PI, 1e-9)
