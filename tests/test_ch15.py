"""Verification suite for Chapter 15 — Compressible Flow (Kundu, Cohen & Dowling 5e, §§15.1–15.10, Eqs. (15.1)–(15.59)).

Evidence levels (``verify-implementation`` skill): V1 analytic · V2 symbolic (sympy / dimensions) · V3 convergence ·
V4 conservation / invariant · V5 published benchmark (``reference/ch15/benchmarks.json``, see SOURCES.md: NACA Report
1135 read first-hand as page images, plus four secondary pages) · V6 book value (private,
``tests/book_values_ch15.json``, skipped when absent) · V7 limits / symmetry / guards.  Every test name is
``test_<concept>_<level>_<what>``; the comment on the ``def`` line repeats the level and the curation ID.

A items (CORE, >= 2 independent levels): C01 the wave equation (15.12) · C02 d'Alembert (15.15) · C03 the duct energy
law (15.22) · C04 stagnation ratios (15.28) · C05 the area ratio (15.31) · C06 the area–velocity relation (15.34) ·
C07 the normal-shock jump (15.39) · C08 the entropy jump (15.43) · C09 the convergent nozzle · C10 the Laval nozzle ·
C11 the duct with friction and heat (15.46) · C12 the friction and heating theorems · C13 the Mach angle (15.49) ·
C14 the oblique shock (15.53) · C15 the Prandtl–Meyer function (15.56) · C16 the supersonic flat plate (15.59).

Derivations: every two- and three-star D row (D01–D03, D06, D07, D10, D11, D13–D15, D17–D24, D26–D30) is re-derived in
a ``test_*_V2_derivation`` test with sympy, **independently of the chapter's own sympy engines** (those are tested
separately); the intermediate lines of the four three-star rows (D13, D17, D21, D30) are checked one by one.

Printed slips kept as wrong variants that must FAIL an independent check: (15.32) third member
(``printed_15_32``), (15.47) first numerator (``printed_15_47``), (15.48) second relation (``printed_15_48``).

All inputs are ours (``ch15.illustrative_inputs()`` or visibly odd numbers); nothing here is quoted from the book.

Run: ``.venv/Scripts/python.exe -m pytest tests/test_ch15.py -q -p no:cacheprovider``  (``-m "not slow"`` skips the
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
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq, minimize_scalar

from fluidpy import ch07_gravity_waves as ch07
from fluidpy import ch15_compressible_flow as ch15
from fluidpy.core import acoustics as AC
from fluidpy.core import gasdynamics as GD
from fluidpy.core import similarity as SIM
from fluidpy.core import thermo as TH
from fluidpy.core.units import Q_, assert_dimensions
from tools.convergence import observed_order, pairwise_orders

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference" / "ch15"
BOOK = Path(__file__).resolve().parent / "book_values_ch15.json"
book_only = pytest.mark.skipif(not BOOK.exists(), reason="book values are private; see CLAUDE.md rule 9")
slow = pytest.mark.slow
PI = math.pi
G = 1.4
I = ch15.illustrative_inputs()
GAMMAS = (1.1, 1.25, 1.4, 5.0 / 3.0)
M_SUB = np.array([0.07, 0.23, 0.5, 0.81, 0.96])                # our Mach numbers, subsonic
M_SUP = np.array([1.04, 1.37, 2.0, 3.3, 6.5])                  # and supersonic
AIR = GD.PerfectGas(1.4, 287.0)


def book():
    return json.loads(BOOK.read_text(encoding="utf-8"))


def bench():
    return json.loads((REF / "benchmarks.json").read_text(encoding="utf-8"))


def rel(a, b):
    return abs(float(a) - float(b)) / max(abs(float(b)), 1e-300)


def close(a, b, rtol=1e-10, atol=0.0):
    """``numpy.allclose`` with NO absolute tolerance unless one is written."""
    return bool(np.allclose(a, b, rtol=rtol, atol=atol))


def z0(expr):
    """True if a sympy expression is identically zero (several simplifiers tried)."""
    for f in (sp.simplify, lambda q: sp.simplify(sp.expand_trig(q)), lambda q: sp.simplify(sp.together(sp.expand(q))),
              lambda q: sp.simplify(sp.trigsimp(sp.expand(q))), lambda q: sp.simplify(sp.powsimp(sp.expand_log(q, force=True), force=True))):
        if f(expr) == 0:
            return True
    return False


def state(M, p=1.0, T=1.0, gamma=G, R=1.0):
    """A state dict (rho, u, p) of a perfect gas at Mach number M — written here, independent of ``GD``."""
    return dict(rho=p / (R * T), u=M * math.sqrt(gamma * R * T), p=p)


def fluxes(s, gamma=G):
    """Mass, momentum and energy fluxes per unit area of a state dict (the three laws, written here)."""
    h = gamma / (gamma - 1.0) * s["p"] / s["rho"]
    return np.array([s["rho"] * s["u"], s["p"] + s["rho"] * s["u"] ** 2, h + 0.5 * s["u"] ** 2])


def t0_ratio(M, g=G):
    return 1.0 + 0.5 * (g - 1.0) * np.asarray(M, float) ** 2


def tol_digits(x, digits):
    """Half a unit of the last printed decimal."""
    return 0.5 * 10.0 ** (-digits)


# =====================================================================================================================
# 0. Contract: every public name exists, is re-exported, has the documented shape, and is used in this file
# =====================================================================================================================
def test_contract_V7_every_public_name_is_exercised_here():  # V7 — design Part C (C.0–C.3): every name the builders call
    src = Path(__file__).read_text(encoding="utf-8")
    body = src.split("# 0. Contract", 1)[1]
    assert (len(GD.__all__), len(AC.__all__), len(ch15._OWN), len(ch15._REUSED)) == (82, 30, 37, 16)
    names = (set(GD.__all__) | set(AC.__all__) | set(ch15._OWN) | set(ch15._REUSED)) - {"AC", "GD"}
    assert names | {"AC", "GD"} == set(ch15.__all__)
    missing = [n for n in sorted(names) if not hasattr(ch15, n)]
    assert not missing, missing
    unused = [n for n in sorted(names) if not re.search(rf"\b(?:ch15|GD|AC)\.{n}\b", body)]
    assert not unused, unused                                            # each name is called by at least one test


def test_contract_V7_reexports_are_the_same_objects():  # V7 — ch15.<name> reaches the tested originals
    for n in GD.__all__:
        assert getattr(ch15, n) is getattr(GD, n), n
    for n in AC.__all__:
        assert getattr(ch15, n) is getattr(AC, n), n
    assert ch15.AC is AC and ch15.GD is GD
    assert ch15.mach_number is SIM.mach_number and ch15.perfect_gas_sound_speed is TH.perfect_gas_sound_speed
    assert ch15.shallow_simple_wave_evolve is ch07.simple_wave_evolve and ch15.hydraulic_jump is ch07.hydraulic_jump
    assert ch15.simple_wave_evolve is AC.simple_wave_evolve              # the gas version keeps the plain name


def test_contract_V7_signatures_carry_their_conventions():  # V7 — convention-carrying arguments are keyword-only
    def kwonly(f):
        return {k for k, p in inspect.signature(f).parameters.items() if p.kind is p.KEYWORD_ONLY}
    assert "branch" in kwonly(GD.mach_from_area_ratio) and inspect.signature(GD.mach_from_area_ratio).parameters["branch"].default is inspect.Parameter.empty
    assert "method" in kwonly(ch15.flat_plate_supersonic) and inspect.signature(ch15.flat_plate_supersonic).parameters["method"].default is inspect.Parameter.empty
    assert {"printed_15_32", "printed_15_47", "printed_15_48", "dA_over_A", "df", "dq"} <= kwonly(GD.influence_coefficients)
    assert {"delta", "sigma", "branch", "gamma"} <= kwonly(GD.oblique_shock) and "allow_expansion" in kwonly(GD.normal_shock)
    assert all(inspect.signature(GD.influence_coefficients).parameters[k].default is False for k in ("printed_15_32", "printed_15_47", "printed_15_48"))
    with pytest.raises(TypeError):
        GD.mach_from_area_ratio(2.0)
    with pytest.raises(TypeError):
        ch15.flat_plate_supersonic(0.05, 2.0)
    for f in (GD.normal_shock, GD.duct_outlet_mach, GD.oblique_shock, GD.prandtl_meyer, AC.dalembert, ch15.flat_plate_supersonic, GD.influence_coefficients):
        d = f.__doc__
        assert "Book:" in d and "Assumptions" in d and "Validation" in d, f.__name__     # CLAUDE.md rule 5 (the label text is reported, not asserted)


def test_contract_V7_scalars_in_floats_out():  # V7 — the parity-row convention (design Part C header)
    for v in (GD.area_ratio(2.0), GD.normal_shock_mach(2.0), GD.prandtl_meyer(2.0), GD.mach_angle(2.0), GD.duct_outlet_mach(0.5, 0.1),
              GD.stagnation_pressure_ratio(0.5), GD.mach_from_area_ratio(2.0, branch="supersonic"), GD.deflection_angle(2.0, 0.8),
              AC.condensation(1.0, 1.2, 340.0), AC.simple_wave_speed(5.0, 340.0), ch15.ackeret_cp(0.03, 2.0), GD.weak_wave_pressure(2.0, 0.03),
              GD.normal_shock(2.0)["p_ratio"], GD.oblique_shock(2.0, delta=0.1)["sigma"], GD.influence_coefficients(0.5, dA_over_A=0.01)["du_over_u"]):
        assert isinstance(v, float), type(v)
    arr = GD.area_ratio(np.array([0.5, 2.0]))
    assert isinstance(arr, np.ndarray) and arr.shape == (2,)


PARITY_QUOTED = {  # numbers the design prints beside its parity rows (ours; a consistency check of design against code)
    'ch15.area_ratio(2.0, 1.4)': 1.6875, 'ch15.mach_from_area_ratio(2.0, 1.4, branch="subsonic")': 0.30590,
    'ch15.mach_from_area_ratio(2.0, 1.4, branch="supersonic")': 2.19720, 'ch15.mass_flux_function(1.0, 1.4)': 0.68473,
    'ch15.laval_area(0.2)': 1.375, 'ch15.normal_shock_entropy(2.0, 1.4)': 0.130916,
    'ch15.normal_shock_stagnation_pressure_ratio(2.0, 1.4)': 0.720874, 'ch15.moving_shock(101325.0, 288.15, p2_over_p1=2.0)["M_s"]': 1.36277,
    'ch15.laval_critical_pressures(2.0, 1.4)["p_b"]': 0.93716, 'ch15.laval_critical_pressures(2.0, 1.4)["p_d"]': 0.51340,
    'ch15.laval_critical_pressures(2.0, 1.4)["p_f"]': 0.09393, 'ch15.laval_shock_location(0.7, 2.0, 1.4)["M1"]': 1.86271,
    'ch15.laval_shock_location(0.7, 2.0, 1.4)["A_shock_over_At"]': 1.51009, 'ch15.laval_exit_state(0.7, 2.0, 1.4)["M_exit"]': 0.40669,
    'ch15.mass_flow_vs_back_pressure(0.8, 1.4)["mdot"]': 0.56066, 'ch15.duct_parameter_A(0.5, 0.1, 0.0, 1.4)': 0.168,
    'ch15.duct_outlet_mach(0.5, 0.1, 0.0, 1.4)': 0.58867, 'ch15.duct_outlet_mach(2.0, 0.1, 0.0, 1.4)': 1.88215,
    'ch15.duct_choking_limits(0.5, 1.4)["f_max"]': 0.22750, 'ch15.mach_cone_delay(1.0e4, 1.5, 295.0)': 25.2663,
    'ch15.oblique_shock(2.0, delta=np.deg2rad(10.0))["M2"]': 1.64052, 'ch15.oblique_shock(2.0, delta=np.deg2rad(10.0))["p_ratio"]': 1.70658,
    'ch15.deflection_angle(2.0, np.pi/4, 1.4)': 0.257324, 'ch15.weak_wave_pressure(2.0, np.deg2rad(10.0), 1.4)': 0.564293,
    'ch15.prandtl_meyer(2.0, 1.4)': 0.460414, 'ch15.prandtl_meyer_max(1.4)': 2.276853,
    'ch15.prandtl_meyer_turn(2.0, np.deg2rad(10.0), 1.4)["M2"]': 2.38489, 'ch15.prandtl_meyer_turn(2.0, np.deg2rad(10.0), 1.4)["p_ratio"]': 0.54797,
    'ch15.ackeret_pressure(np.deg2rad(2.0), 2.0, 1.4)': 0.112859, 'ch15.ackeret_cp(np.deg2rad(2.0), 2.0)': 0.040307,
    'ch15.simple_wave_speed(5.0, 340.0, 1.4)': 346.0,
}


def test_contract_V7_design_parity_expressions_evaluate():  # V7 — every `py:` expression of design Part B runs and is finite
    txt = (ROOT / "analysis" / "ch15_design.md").read_text(encoding="utf-8")
    exprs = sorted({" ".join(e.split()) for e in re.findall(r"py:\s*'([^']+)'", txt)})
    assert len(exprs) >= 60, len(exprs)
    ns = {"ch15": ch15, "np": np}
    values = {}
    for e in exprs:
        v = eval(e, ns)                                                   # noqa: S307 — our own committed design file
        assert np.all(np.isfinite(np.asarray(v, dtype=float))), e
        values[e] = v
    hit = 0
    for e, quoted in PARITY_QUOTED.items():
        if e in values:
            hit += 1
            assert rel(values[e], quoted) < 6e-5, (e, values[e], quoted)   # the design prints 5–6 significant digits
    assert hit >= 25, hit


# =====================================================================================================================
# 1. §15.1 — the gas, the Mach number, the regimes  (N02–N09, R01–R08; T1, T2)
# =====================================================================================================================
def test_perfect_gas_V1_specific_heats_and_sound_speed():  # V1 — (15.1) specific-heat and sound-speed lines; R01, R02, R07
    for g in GAMMAS:
        gas = GD.PerfectGas(g, 311.0)
        assert close(gas.cp - gas.cv, 311.0, 1e-13) and close(gas.cp / gas.cv, g, 1e-13)
        T = np.array([211.0, 288.15, 933.0])
        assert close(gas.c(T), ch15.perfect_gas_sound_speed(T, g, 311.0), 1e-13) and isinstance(gas.c(300.0), float)
        p0, r0 = 0.9e5, 1.1                                                  # c^2 = (dp/drho)_s on p = p0 (rho/rho0)^gamma  (1.19)
        c_num = ch15.sound_speed_from_eos(lambda r, s: ch15.isentropic_pressure(r, p0, r0, g), r0)
        assert rel(c_num, math.sqrt(g * p0 / r0)) < 1e-6
        assert close(GD.mach_from_speed(250.0, 263.0, gas), 250.0 / math.sqrt(g * 311.0 * 263.0), 1e-13)
    assert GD.AIR_BOOK.gamma == 1.4 and GD.AIR_BOOK.R == 287.0 and close(GD.AIR_BOOK.cp, 1004.5, 1e-13) and close(GD.AIR_BOOK.cv, 717.5, 1e-13)
    assert close(GD.mach_from_speed(170.0, 288.15), ch15.mach_number(170.0, T=288.15), 1e-13)   # default gas = project air
    for bad in ((1.0, 287.0), (0.9, 287.0), (1.4, 0.0)):
        with pytest.raises(ValueError):
            GD.PerfectGas(*bad)


def test_perfect_gas_V2_dimensions_of_the_dimensional_wrappers():  # V2 — units of c, mass flow, wall force, enthalpy, choked flow
    gm, R, T, p0, A = 1.4, Q_(287.0, "J/(kg*K)"), Q_(500.0, "K"), Q_(5.0e5, "Pa"), Q_(0.01, "m**2")
    c = (gm * R * T) ** 0.5
    assert_dimensions(c, "velocity", "c")
    assert close(c.to("m/s").magnitude, AIR.c(500.0), 1e-13)
    mdot = p0 * A * (gm / (R * T)) ** 0.5 * (2 / (gm + 1)) ** ((gm + 1) / (2 * (gm - 1)))
    assert_dimensions(mdot, "mass_flow", "choked flow")
    assert close(mdot.to("kg/s").magnitude, GD.choked_mass_flow(5.0e5, 500.0, 0.01, AIR), 1e-12)
    rho, u = Q_(1.2, "kg/m**3"), Q_(300.0, "m/s")
    assert_dimensions(rho * u * A, "mass_flow")
    assert_dimensions(rho * u * A * u + p0 * A, "force")
    assert close(GD.mass_flow(1.2, 300.0, 0.01), 3.6, 1e-14)
    assert_dimensions(Q_(1004.5, "J/(kg*K)") * T + 0.5 * u ** 2, "[length]**2/[time]**2", "h0")
    assert close(GD.total_enthalpy(500.0, 300.0, 1004.5), 1004.5 * 500.0 + 45000.0, 1e-14)
    assert close(GD.total_enthalpy(500.0, 300.0, 1004.5), ch15.stagnation_enthalpy(1004.5 * 500.0, 300.0), 1e-14)
    assert_dimensions(rho * c, "[mass]/([length]**2*[time])", "impedance")
    assert close(AC.acoustic_impedance(1.2, 340.0), 408.0, 1e-14)
    # homogeneity of the functions themselves: a change of units is a rescaling (V7)
    assert close(GD.choked_mass_flow(5.0, 500.0, 0.01, AIR), 1e-5 * GD.choked_mass_flow(5.0e5, 500.0, 0.01, AIR), 1e-13)
    assert close(GD.choked_mass_flow(5.0e5, 4 * 500.0, 0.01, AIR), 0.5 * GD.choked_mass_flow(5.0e5, 500.0, 0.01, AIR), 1e-13)   # trap T10: ~ p0/sqrt(T0)


def test_mach_regime_V1_both_sides_of_every_threshold():  # V1 (demonstration) — N02, trap T3
    assert ch15.mach_regime(0.0) == "incompressible" and ch15.mach_regime(1e-9) == "subsonic" and ch15.mach_regime(0.7999) == "subsonic"
    assert ch15.mach_regime(0.8) == "transonic" and ch15.mach_regime(1.0) == "transonic" and ch15.mach_regime(1.2) == "transonic"
    assert ch15.mach_regime(1.2001) == "supersonic" and ch15.mach_regime(4.999) == "supersonic" and ch15.mach_regime(5.0) == "hypersonic"
    assert ch15.mach_regime(0.85, transonic=(0.9, 1.1)) == "subsonic" and ch15.mach_regime(3.5, hypersonic=3.2) == "hypersonic"
    for bad in (-0.1, float("nan"), float("inf")):
        with pytest.raises(ValueError):
            ch15.mach_regime(bad)
    cc = ch15.compressibility_check(0.3 * 340.0)                             # R03: the incompressible limit of ch14
    assert close(cc["Ma"], GD.mach_from_speed(0.3 * 340.0, 288.15), 1e-6) and abs(cc["density_change"] - (GD.stagnation_density_ratio(cc["Ma"]) - 1)) < 0.25 * cc["Ma"] ** 4


def test_entropy_line_V1_two_members_agree_and_isentrope():  # V1 — (15.1) entropy line, (15.2); R04, R05
    rng = np.random.default_rng(15)
    for _ in range(20):
        T1, T2 = rng.uniform(200, 900, 2)
        p1, p2 = rng.uniform(2e4, 9e5, 2)
        R, cp = 287.0, 1004.5
        v1, v2 = R * T1 / p1, R * T2 / p2
        a = ch15.perfect_gas_entropy_change_p(T1, p1, T2, p2, cp=cp, R=R)
        b = ch15.perfect_gas_entropy_change(T1, v1, T2, v2, cv=cp - R, R=R)
        assert abs(a - b) < 1e-9 * max(1.0, abs(a))
    Tr, rr = ch15.isentropic_ratios(0.37, 1.4)
    assert close(Tr, 0.37 ** (0.4 / 1.4), 1e-13) and close(rr, 0.37 ** (1 / 1.4), 1e-13) and close(Tr * rr, 0.37, 1e-13)     # p = rho R T
    assert abs(ch15.perfect_gas_entropy_change_p(300.0, 1e5, 300.0 * Tr, 0.37e5, cp=1004.5, R=287.0)) < 1e-9             # an isentrope keeps S


# =====================================================================================================================
# 2. §15.2 — acoustics: C01 the wave equation (15.12), C02 d'Alembert (15.15)–(15.17); steepening (ours)
# =====================================================================================================================
def test_acoustic_relation_V1_linearisation_error_is_the_quadratic_term():  # V1 — (15.5), the Taylor expansion; N10–N13 (T3)
    rho0, p0 = 1.21, 1.0e5
    for g in GAMMAS:
        amp = rho0 * np.array([1e-1, 1e-2, 1e-3, 1e-4])
        ex = AC.isentropic_pressure_perturbation(amp, rho0, p0, g, "exact")
        e1 = ex - AC.isentropic_pressure_perturbation(amp, rho0, p0, g, 1)
        e2 = ex - AC.isentropic_pressure_perturbation(amp, rho0, p0, g, 2)
        assert abs(observed_order(amp[1:], np.abs(e1[1:])) - 2.0) < 0.05          # the dropped term is quadratic …
        assert close(e1[2:] / amp[2:] ** 2, 0.5 * g * (g - 1) * p0 / rho0 ** 2, 2e-3)   # … with the coefficient gamma(gamma-1)p0/(2 rho0^2)
        assert abs(observed_order(amp[:3], np.abs(e2[:3])) - 3.0) < 0.1            # order 2 leaves a cubic remainder
        c = math.sqrt(g * p0 / rho0)
        assert close(AC.isentropic_pressure_perturbation(amp, rho0, p0, g), c ** 2 * amp, 1e-13)        # (15.5) with (15.6)
        assert close(AC.acoustic_density(c ** 2 * amp, c), amp, 1e-13)
        pp = c ** 2 * amp
        assert close(AC.linearisation_error(pp, rho0, c, g), (0.5 * g * (g - 1) * p0 / rho0 ** 2 * amp ** 2) / pp, 1e-12)   # dropped / kept
        assert close(AC.condensation(pp, rho0, c), amp / rho0, 1e-13)
        Tex = 300.0 * ((1 + pp / p0) ** ((g - 1) / g) - 1)                        # exact isentrope (15.2)
        assert close(AC.acoustic_temperature(pp[2:], p0, 300.0, g), Tex[2:], 2e-3) and abs(AC.acoustic_temperature(pp[3], p0, 300.0, g) / Tex[3] - 1) < 1e-4
        assert close(AC.newton_sound_speed(288.15, 287.0) / ch15.perfect_gas_sound_speed(288.15, g, 287.0), 1 / math.sqrt(g), 1e-13)   # Newton vs Laplace
    with pytest.raises(ValueError):
        AC.isentropic_pressure_perturbation(0.01, rho0, p0, 1.4, 3)
    assert close(AC.sound_pressure_level(1.0), 20 * math.log10(1 / 2e-5), 1e-14) and close(AC.sound_pressure_level(2e-5), 0.0, 1e-14, 1e-12)
    assert close(AC.sound_pressure_level(10.0) - AC.sound_pressure_level(1.0), 20.0, 1e-12)


def _lin3d():
    """Index-form bookkeeping for D01/D02 (three space dimensions, uniform stream U_j, constants rho0, p0, c)."""
    t, eps = sp.symbols("t epsilon", real=True)
    X = sp.symbols("x1:4", real=True)
    U = sp.symbols("U1:4", real=True)
    rho0, p0, c = sp.symbols("rho_0 p_0 c", positive=True)
    r = sp.Function("r")(*X, t)
    p = sp.Function("p")(*X, t)
    u = [sp.Function(f"u{i}")(*X, t) for i in (1, 2, 3)]
    D = lambda q: sp.diff(q, t) + sum(U[j] * sp.diff(q, X[j]) for j in range(3))     # noqa: E731
    return t, eps, X, U, rho0, p0, c, r, p, u, D


def test_linearised_laws_V2_derivation():  # V2 — D01 (two stars): (15.4) in continuity and Euler → (15.8), (15.9) → (15.10), (15.11); C01
    t, eps, X, U, rho0, p0, c, r, p, u, D = _lin3d()
    g_ = sp.Symbol("gamma", positive=True)
    # step 1–2: Taylor expansion of p(rho) at constant entropy for p = p0 (rho/rho0)^gamma: p' = c^2 rho' + ..., c^2 = gamma p0/rho0  (15.5), (15.6)
    rp = sp.Symbol("rhop")
    ser = sp.series(p0 * ((rho0 + rp) / rho0) ** g_ - p0, rp, 0, 3).removeO()
    assert z0(ser.coeff(rp, 1) - g_ * p0 / rho0) and z0(ser.coeff(rp, 2) - g_ * (g_ - 1) * p0 / (2 * rho0 ** 2))
    # steps 3–6: the split (15.4) with a bookkeeping epsilon
    rho, pres, vel = rho0 + eps * r, p0 + eps * p, [U[i] + eps * u[i] for i in range(3)]
    cont = sp.expand(sp.diff(rho, t) + sum(sp.diff(rho * vel[j], X[j]) for j in range(3)))
    assert cont.coeff(eps, 0) == 0                                                   # the zeroth order vanishes: U_j, rho0 are constants
    assert z0(cont.coeff(eps, 1) - (D(r) + rho0 * sum(sp.diff(u[j], X[j]) for j in range(3))))          # (15.8), last member
    assert z0(cont.coeff(eps, 2) - sum(sp.diff(r * u[j], X[j]) for j in range(3)))                       # what is dropped is second order
    for i in range(3):
        mom = sp.diff(vel[i], t) + sum(vel[j] * sp.diff(vel[i], X[j]) for j in range(3)) + sp.diff(pres, X[i]) / rho
        m0, m1 = mom.subs(eps, 0).doit(), sp.diff(mom, eps).subs(eps, 0).doit()
        assert z0(m0) and z0(m1 - (D(u[i]) + sp.diff(p, X[i]) / rho0))                                    # (15.9), last member
        m2 = sp.diff(mom, eps, 2).subs(eps, 0).doit() / 2
        assert z0(m2 - (sum(u[j] * sp.diff(u[i], X[j]) for j in range(3)) - r * sp.diff(p, X[i]) / rho0 ** 2))
    # steps 7–10: p' = c^2 rho' (15.5) removes rho' from (15.8)  →  (15.10); (15.9) is (15.11)
    e10 = (D(r) + rho0 * sum(sp.diff(u[j], X[j]) for j in range(3))).subs(r, p / c ** 2).doit()
    assert z0(e10 - (D(p) / c ** 2 + rho0 * sum(sp.diff(u[j], X[j]) for j in range(3))))


def test_wave_equation_V2_derivation():  # V2 — D02 (two stars): D(15.10) − rho0 ∂_i(15.11) = (15.12); C01
    t, eps, X, U, rho0, p0, c, r, p, u, D = _lin3d()
    e10 = D(p) / c ** 2 + rho0 * sum(sp.diff(u[j], X[j]) for j in range(3))
    e11 = [D(u[i]) + sp.diff(p, X[i]) / rho0 for i in range(3)]
    wave = D(D(p)) / c ** 2 - sum(sp.diff(p, X[j], 2) for j in range(3))             # (15.12)
    assert z0(sp.expand(D(e10) - rho0 * sum(sp.diff(e11[i], X[i]) for i in range(3)) - wave))
    # the two operations commute only because U_j is constant: with U_1 = U_1(x_1) the elimination leaves a remainder
    Uv = sp.Function("V")(X[0])
    Dv = lambda q: sp.diff(q, t) + Uv * sp.diff(q, X[0])                              # noqa: E731
    rem = sp.expand(Dv(Dv(p) / c ** 2 + rho0 * sp.diff(u[0], X[0])) - rho0 * sp.diff(Dv(u[0]) + sp.diff(p, X[0]) / rho0, X[0])
                    - (Dv(Dv(p)) / c ** 2 - sp.diff(p, X[0], 2)))
    assert not z0(rem)
    # a plane wave in any direction n (|n| = 1) with speed c + U·n solves (15.12)
    n = sp.symbols("n1:4", real=True)
    f = sp.Function("f")
    arg = sum(n[j] * X[j] for j in range(3)) - (c + sum(n[j] * U[j] for j in range(3))) * t
    res = (D(D(f(arg))) / c ** 2 - sum(sp.diff(f(arg), X[j], 2) for j in range(3))).doit()
    assert z0(res.subs(n[2], sp.sqrt(1 - n[0] ** 2 - n[1] ** 2)))


def test_wave_equation_V2_chapter_engine_agrees():  # V2 — ch15.linearised_acoustics_sympy (the notebook's check cell); C01
    e = ch15.linearised_acoustics_sympy()
    s = e["symbols"]
    x, t, U, c, rho0 = s["x"], s["t"], s["U"], s["c"], s["rho0"]
    r, u, p = s["rho1"], s["u1"], s["p1"]
    assert z0(e["continuity_linear"] - (sp.diff(r, t) + U * sp.diff(r, x) + rho0 * sp.diff(u, x)))
    assert z0(e["continuity_dropped"] - sp.diff(r * u, x)) and z0(e["momentum_linear"] - (sp.diff(u, t) + U * sp.diff(u, x) + sp.diff(p, x) / rho0))
    assert z0(e["momentum_dropped"] - (u * sp.diff(u, x) - r * sp.diff(p, x) / rho0 ** 2))
    assert e["cross_residual"] == 0 and e["dalembert_residual"] == 0
    assert z0(e["wave_operator"] - ((sp.diff(p, t, 2) + 2 * U * sp.diff(p, x, t) + U ** 2 * sp.diff(p, x, 2)) / c ** 2 - sp.diff(p, x, 2)))
    assert z0(e["eq_15_10"] - ((sp.diff(p, t) + U * sp.diff(p, x)) / c ** 2 + rho0 * sp.diff(u, x))) and z0(e["eq_15_11"] - e["momentum_linear"])


def test_wave_equation_V3_residual_of_a_travelling_pulse_vanishes_at_second_order():  # V3 — AC.convected_wave_residual on (15.17); C01
    c, U = 340.0, 170.0
    f, g_ = AC.pulse("gaussian", 0.5, 20.0), AC.pulse("packet", 0.7, 8.0, x0=1.0)
    x = np.linspace(-2.0, 3.0, 41)
    sol = lambda xx, tt: AC.dalembert_in_stream(xx, tt, c, U, f, g_)                  # noqa: E731
    hs = np.array([0.08, 0.04, 0.02, 0.01])
    res = [np.max(np.abs(AC.convected_wave_residual(sol, x, 1.3e-3, c, U, h))) for h in hs]
    scale = np.max(np.abs((f(x + 1e-3) - 2 * f(x) + f(x - 1e-3)) / 1e-6))             # size of p_xx
    assert abs(observed_order(hs, res) - 2.0) < 0.15 and res[-1] < 2e-3 * scale
    wrong = lambda xx, tt: AC.dalembert_in_stream(xx, tt, c, -U, f, g_)               # noqa: E731 — the stream reversed: not a solution
    assert np.max(np.abs(AC.convected_wave_residual(wrong, x, 1.3e-3, c, U, 0.01))) > 0.5 * scale
    still = lambda xx, tt: AC.dalembert(xx, tt, c, f, g_)                             # noqa: E731
    r0 = [np.max(np.abs(AC.wave_equation_residual(still, x, 1.3e-3, c, h))) for h in hs]
    assert max(r0) < 1e-9 * scale        # U → 0 gives (15.13); with the stencil's own time step h/c the differences are exact for f(x ∓ ct): round-off only
    assert close(AC.wave_equation_residual(still, x, 1.3e-3, c, 0.02), AC.convected_wave_residual(still, x, 1.3e-3, c, 0.0, 0.02), 1e-14, 1e-14)


def test_dalembert_V2_derivation():  # V2 — D03 (two stars): characteristic coordinates turn (15.13) into p_xi_zeta = 0 → (15.15); C02
    x, t = sp.symbols("x t", real=True)
    c = sp.Symbol("c", positive=True)
    xi, ze = sp.symbols("xi zeta", real=True)
    P = sp.Function("P")
    sub = {xi: x - c * t, ze: x + c * t}
    p_xt = P(xi, ze)
    # chain rule, one derivative at a time (steps 2–5)
    px = sp.diff(p_xt, xi) + sp.diff(p_xt, ze)
    pt = -c * sp.diff(p_xt, xi) + c * sp.diff(p_xt, ze)
    pxx = sp.diff(px, xi) + sp.diff(px, ze)
    ptt = -c * sp.diff(pt, xi) + c * sp.diff(pt, ze)
    assert z0(sp.expand(ptt / c ** 2 - pxx) + 4 * sp.diff(p_xt, xi, ze))              # (15.13) ≡ −4 P_xi_zeta: only the mixed term survives
    # steps 6–8: P_xi_zeta = 0 integrates to f(xi) + g(zeta); check the result in the original variables for arbitrary shapes
    f, g_ = sp.Function("f"), sp.Function("g")
    sol = f(x - c * t) + g_(x + c * t)
    assert z0((sp.diff(sol, t, 2) / c ** 2 - sp.diff(sol, x, 2)).doit())
    assert sp.diff(f(xi) + g_(ze), xi, ze) == 0
    # and the particle velocity (15.16) satisfies the linear momentum law (15.11) and continuity (15.10) with U = 0  (D04)
    rho0 = sp.Symbol("rho_0", positive=True)
    uu = (f(x - c * t) - g_(x + c * t)) / (rho0 * c)
    assert z0((sp.diff(uu, t) + sp.diff(sol, x) / rho0).doit()) and z0((sp.diff(sol, t) / c ** 2 + rho0 * sp.diff(uu, x)).doit())
    wrong = (f(x - c * t) + g_(x + c * t)) / (rho0 * c)                               # planted: the sign of g lost
    assert not z0((sp.diff(wrong, t) + sp.diff(sol, x) / rho0).doit())


def test_dalembert_V1_shapes_travel_unchanged():  # V1 — (15.15): f moves right at c, g left at c; AC.pulse shapes; C02
    c = 340.0
    x = np.linspace(-3, 3, 601)
    for kind in AC.PULSE_KINDS:
        f = AC.pulse(kind, 0.5, 20.0, x0=0.3)
        tt = 2.0e-3
        assert close(AC.dalembert(x, tt, c, f=f), f(x - c * tt), 1e-14, 1e-14) and close(AC.dalembert(x, tt, c, g=f), f(x + c * tt), 1e-14, 1e-14)
        assert close(AC.dalembert(x + c * tt, tt, c, f=f), AC.dalembert(x, 0.0, c, f=f), 1e-12, 1e-12)
    assert close(AC.pulse("gaussian", 0.5, 20.0, x0=0.3)(0.3), 20.0, 1e-14) and close(AC.pulse("gaussian", 0.5, 20.0)(0.5), 20.0 / math.e, 1e-13)
    assert abs(AC.pulse("tophat", 1.0, 3.0)(0.0) - 3.0) < 1e-3 and abs(AC.pulse("tophat", 1.0, 3.0)(2.0)) < 1e-10 and AC.pulse("sine", 1.0, 3.0)(0.6) == 0.0
    assert close(AC.pulse("sine", 2.0, 3.0)(0.5), 3.0, 1e-13) and close(AC.pulse("packet", 1.0, 3.0)(0.0), 3.0, 1e-14)
    with pytest.raises(ValueError):
        AC.pulse("square")
    assert AC.dalembert(0.1, 0.2, c) == 0.0 and isinstance(AC.dalembert(0.1, 0.0, c, f=AC.pulse()), float)


def test_particle_velocity_V1_equals_the_time_integral_of_the_pressure_gradient():  # V1 — (15.16) against the quadrature of (15.14); C02, T19
    c, rho0 = 340.0, 1.2
    f, g_ = AC.pulse("gaussian", 0.4, 20.0, x0=-3.0), AC.pulse("gaussian", 0.6, 12.0, x0=4.0)     # both far from x at t = 0
    p_fn = lambda xx, tt: AC.dalembert(xx, tt, c, f, g_)                                # noqa: E731
    for xo, tt in ((0.2, 9.0e-3), (-0.5, 7.5e-3), (1.0, 10.5e-3), (0.0, 12.0e-3)):
        u_quad = AC.velocity_from_pressure(p_fn, xo, tt, rho0, t0=0.0, h=1e-4, limit=400)
        u_cf = AC.plane_wave_velocity(xo, tt, c, rho0, f, g_)
        assert abs(u_quad - u_cf) < 1e-7 * 20.0 / (rho0 * c)
    x = np.linspace(-6, 6, 25)
    assert np.all(AC.plane_wave_velocity(x, 0.0, c, rho0, f=f) >= 0) and np.all(AC.plane_wave_velocity(x, 0.0, c, rho0, g=g_) <= 0)   # compression: fluid moves with the wave
    assert close(AC.plane_wave_velocity(x, 0.0, c, rho0, f=f) * AC.acoustic_impedance(rho0, c), f(x), 1e-13, 1e-300)             # p' = rho0 c u'
    assert close(AC.plane_wave_velocity(x, 0.0, c, rho0, f=f) / c, AC.condensation(f(x), rho0, c), 1e-13, 1e-300)                 # u'/c = condensation


def test_acoustic_energy_V4_is_constant_while_two_pulses_cross():  # V4 — invariant of (15.10)–(15.11); independent of how (15.16) was derived; C02
    c, rho0 = 340.0, 1.2
    f, g_ = AC.pulse("gaussian", 0.4, 20.0, x0=-1.0), AC.pulse("gaussian", 0.4, 14.0, x0=1.2)
    x = np.linspace(-12, 12, 24001)
    E, Ebad = [], []
    for tt in np.linspace(0.0, 8e-3, 9):
        p = AC.dalembert(x, tt, c, f, g_)
        E.append(AC.acoustic_energy(p, AC.plane_wave_velocity(x, tt, c, rho0, f, g_), x, rho0, c))
        Ebad.append(AC.acoustic_energy(p, (f(x - c * tt) + g_(x + c * tt)) / (rho0 * c), x, rho0, c))     # planted: (f + g)/(rho0 c)
    assert np.ptp(E) < 1e-10 * E[0] and np.ptp(Ebad) > 0.05 * E[0]
    exact = quad(lambda s: f(s) ** 2, -8, 8, epsabs=1e-14, epsrel=1e-13, limit=400)[0] + quad(lambda s: g_(s) ** 2, -8, 8, epsabs=1e-14, epsrel=1e-13, limit=400)[0]
    assert rel(E[0], exact / (rho0 * c ** 2)) < 1e-8                                  # each part carries ∫ shape^2 dx / (rho0 c^2)


def test_dalembert_initial_value_V1_reproduces_both_initial_conditions():  # V1 + V2 — Exercise-type formulas f, g from F, G; N17, N18
    c = 3.0
    F = lambda s: np.exp(-np.asarray(s, float) ** 2)                                   # noqa: E731
    Gf = lambda s: 2.0 * np.asarray(s, float) * np.exp(-np.asarray(s, float) ** 2)     # noqa: E731
    x = np.linspace(-2, 2, 17)
    assert close(AC.dalembert_initial_value(x, 0.0, c, F, Gf), F(x), 1e-13)
    dt = 1e-5
    pt = (AC.dalembert_initial_value(x, dt, c, F, Gf) - AC.dalembert_initial_value(x, -dt, c, F, Gf)) / (2 * dt)
    assert close(pt, Gf(x), 1e-6, 1e-7)
    f, g_ = AC.dalembert_fg(F, Gf, c, x_ref=-0.7)
    for tt in (0.0, 0.11, 0.4):
        assert close(AC.dalembert(x, tt, c, f, g_), AC.dalembert_initial_value(x, tt, c, F, Gf), 1e-10, 1e-12)
    assert close(AC.dalembert_initial_value(x, 0.3, c, F), 0.5 * (F(x - 0.9) + F(x + 0.9)), 1e-13)     # released from rest: two equal halves
    res = AC.wave_equation_residual(lambda xx, tt: AC.dalembert_initial_value(xx, tt, c, F, Gf), x, 0.2, c, 0.01)
    assert np.max(np.abs(res)) < 5e-3
    xs, ts, cs = sp.symbols("x t c", positive=True)                                    # sympy, a polynomial pair (V2)
    Fs, Gs = lambda s: 1 + s ** 2, lambda s: 3 * s ** 2                                # noqa: E731
    z = sp.Symbol("z")
    sol = (Fs(xs - cs * ts) + Fs(xs + cs * ts)) / 2 + sp.integrate(Gs(z), (z, xs - cs * ts, xs + cs * ts)) / (2 * cs)
    assert z0(sp.diff(sol, ts, 2) / cs ** 2 - sp.diff(sol, xs, 2)) and z0(sol.subs(ts, 0) - Fs(xs)) and z0(sp.diff(sol, ts).subs(ts, 0) - Gs(xs))
    num = AC.dalembert_initial_value(0.4, 0.3, 2.0, lambda s: 1 + np.asarray(s) ** 2, lambda s: 3 * np.asarray(s) ** 2)
    assert rel(num, sol.subs({xs: 0.4, ts: 0.3, cs: 2.0})) < 1e-12


def test_waves_in_a_stream_V1_galilean_shift_of_dalembert():  # V1 + V7 — (15.17) and the two speeds U ± c; D05, N19, N20
    c = 340.0
    f, g_ = AC.pulse("gaussian", 0.5, 20.0), AC.pulse("sine", 1.5, 9.0, x0=0.4)
    x = np.linspace(-4, 4, 161)
    for U in (0.0, 170.0, 510.0, -90.0):
        for tt in (0.0, 1.0e-3, 4.0e-3):
            assert close(AC.dalembert_in_stream(x, tt, c, U, f, g_), AC.dalembert(x - U * tt, tt, c, f, g_), 1e-13, 1e-13)
        sp_, sm = AC.wave_speeds(U, c)
        assert sp_ == U + c and sm == U - c and AC.upstream_signal_possible(U, c) == (abs(U) < c)
        ch = AC.xt_characteristics(U, c, [-1.0, 0.5], np.linspace(0, 5e-3, 6))
        assert ch["x_plus"].shape == (2, 6) and ch["speed_plus"] == U + c and ch["speed_minus"] == U - c
        assert close(AC.dalembert_in_stream(ch["x_plus"], ch["t"][None, :], c, U, f=f), f(np.array([-1.0, 0.5]))[:, None] * np.ones((1, 6)), 1e-12, 1e-12)   # f constant on its lines
        assert close(AC.dalembert_in_stream(ch["x_minus"], ch["t"][None, :], c, U, g=g_), g_(np.array([-1.0, 0.5]))[:, None] * np.ones((1, 6)), 1e-12, 1e-12)
    assert AC.wave_speeds(I["U_stream_sub"], c)[1] < 0 < AC.wave_speeds(I["U_stream_sup"], c)[1]            # U − c changes sign at U = c
    assert AC.wave_speeds(c, c)[1] == 0.0 and not AC.upstream_signal_possible(c, c)
    assert list(AC.upstream_signal_possible(np.array([100.0, 400.0]), c)) == [True, False]


def test_leapfrog_V3_second_order_and_exact_at_courant_one():  # V3 — AC.wave_equation_step (our scheme, DEVIATION); N114
    shape = lambda s: np.exp(-((((np.asarray(s) + 0.5) % 1.0) - 0.5) / 0.08) ** 2)     # noqa: E731 — periodic Gaussian
    orders = {}
    for C in (0.5, 0.8):
        hs, errs = [], []
        for n in (100, 200, 400, 800):
            x = np.arange(n) / n
            dt = C / n
            ns = int(round(0.4 / dt))
            p_prev, p_now = shape(x), shape(x - dt)                                     # exact start at two levels (c = 1)
            for _ in range(ns - 1):
                p_now, p_prev = AC.wave_equation_step(p_now, p_prev, C), p_now
            hs.append(1.0 / n)
            errs.append(np.max(np.abs(p_now - shape(x - ns * dt))))
        orders[C] = observed_order(hs, errs)
        assert abs(orders[C] - 2.0) < 0.15, (C, errs)
    n = 64                                                                             # C = 1: exact translation
    x = np.arange(n) / n
    p_prev, p_now = shape(x), shape(x - 1.0 / n)
    for _ in range(150):
        p_now, p_prev = AC.wave_equation_step(p_now, p_prev, 1.0), p_now
    assert np.max(np.abs(p_now - shape(x - 151.0 / n))) < 1e-11
    with pytest.raises(ValueError):
        AC.wave_equation_step(p_now, p_prev, 1.0001)                                   # the stability limit is asserted
    with pytest.raises(ValueError):
        AC.wave_equation_step(p_now, p_prev, 0.5, bc="open")
    q = AC.wave_equation_step(np.sin(PI * x), np.sin(PI * x), 0.5, bc="fixed")
    assert q[0] == 0.0 and q[-1] == 0.0
    # just above the limit the scheme blows up (the roll-based step itself, without the guard): amplification of the sawtooth mode
    C = 1.02
    saw = (-1.0) ** np.arange(n)
    amp = np.max(np.abs(2 * saw - saw + C ** 2 * (np.roll(saw, -1) - 2 * saw + np.roll(saw, 1))))
    assert amp > 3.0 and abs(2 - 1 + 1.0 ** 2 * (-4)) == 3.0                           # |1 − 4C^2| > 3 for C > 1, = 3 at C = 1 (neutral with p_prev = p_now)


def test_simple_wave_V4_riemann_invariant_and_isentrope():  # V4 — the invariant u − 2c/(gamma − 1) across the wave (ours); N21, N22
    c0 = 340.0
    x = np.linspace(0.0, 1.0, 801)
    for g in (1.4, 5.0 / 3.0, 2.0):
        u0 = I["steep_u"] * np.sin(2 * PI * x)
        s = AC.simple_wave_evolve(u0, x, 1.0e-4, c0, g)
        assert np.ptp(s["u"] - 2 * s["c"] / (g - 1)) < 1e-10 * c0 and close(s["u"] - 2 * s["c"] / (g - 1), -2 * c0 / (g - 1), 1e-12)
        assert close(s["p_over_p0"], s["rho_over_rho0"] ** g, 1e-12) and close(s["c"] ** 2 / c0 ** 2, s["p_over_p0"] / s["rho_over_rho0"], 1e-12)
        assert close(s["speed"], s["u"] + s["c"], 1e-13) and close(s["speed"], AC.simple_wave_speed(u0, c0, g), 1e-14)
        assert close(s["x_points"], x + (c0 + 0.5 * (g + 1) * u0) * 1.0e-4, 1e-14) and s["x"] is s["x_points"] and s["valid"] is True
        tb = 2.0 / ((g + 1) * I["steep_u"] * 2 * PI)                                   # sine: t_b = 2/((gamma + 1) a k)
        assert rel(s["t_break"], tb) < 1e-4 and close(AC.breaking_time(None, None, g, max_compression_rate=I["steep_u"] * 2 * PI), tb, 1e-14)
        many = AC.simple_wave_evolve(u0, x, np.array([0.5 * tb, 2.0 * tb]), c0, g)
        assert many["x_points"].shape == (2, 801) and list(many["valid"]) == [True, False]
        assert np.all(np.diff(many["x_points"][0]) > 0) and np.any(np.diff(many["x_points"][1]) < 0)      # single-valued before t_b, folded after
    assert AC.breaking_time(np.linspace(0, 1, 9), np.linspace(0, 1, 9)) == math.inf                       # a pure expansion never breaks
    assert close(AC.simple_wave_speed(1e-9, c0), c0, 1e-10)                                                # small amplitude travels at c0


def test_breaking_time_V1_first_crossing_of_neighbouring_characteristics():  # V1 — t_b against a direct search for the first crossing; N21
    c0, g = 340.0, 1.4
    x = np.linspace(-4.0, 4.0, 4001)
    u0 = 12.0 * np.exp(-x ** 2) + 3.0 * np.exp(-(x - 1.5) ** 2 / 0.3)
    speed = AC.simple_wave_speed(u0, c0, g)
    dspeed, dx = np.diff(speed), np.diff(x)
    t_cross = np.min(np.where(dspeed < 0, -dx / np.where(dspeed < 0, dspeed, -1.0), np.inf))    # first time two neighbours meet
    tb = AC.breaking_time(u0, x, g)
    assert rel(tb, t_cross) < 2e-4
    xf = np.linspace(-4.0, 4.0, 400001)
    slope = np.max(-np.gradient(12.0 * np.exp(-xf ** 2) + 3.0 * np.exp(-(xf - 1.5) ** 2 / 0.3), xf))
    assert rel(tb, 2.0 / ((g + 1) * slope)) < 1e-5


def test_simple_wave_V1_gamma_two_is_the_shallow_water_wave_of_chapter_7():  # V1 — the hydraulic analogy (gamma = 2, c0 = sqrt(gH)); N22, R16
    H, grav = 1.3, 9.81
    x = np.linspace(-6.0, 6.0, 1201)
    eta0 = 0.11 * np.exp(-x ** 2)
    sw = ch15.shallow_simple_wave_evolve(eta0, x, 0.35, H, grav)
    c0 = math.sqrt(grav * H)
    u0 = 2.0 * (np.sqrt(grav * (H + eta0)) - c0)
    gas = AC.simple_wave_evolve(u0, x, 0.35, c0, 2.0)
    assert close(gas["x_points"], sw["x_points"], 1e-12) and close(gas["c"], np.sqrt(grav * (H + eta0)), 1e-12)
    assert rel(gas["t_break"], sw["t_break"]) < 1e-6
    assert close(gas["rho_over_rho0"], (H + eta0) / H, 1e-12)                          # the depth plays the density
    j = ch15.hydraulic_jump(1.0, Fr1=2.0)                                              # T21: the jump loses energy, the shock stagnation pressure
    assert j["dE"] < 0 and j["H2"] > 1.0 and GD.normal_shock(2.0)["p0_ratio"] < 1.0 and GD.normal_shock(2.0)["T0_ratio"] == 1.0


def test_mach_cone_V1_envelope_of_the_wave_fronts():  # V1 — D25, C13: the common tangent of the circles makes sin(mu) = 1/M with the path
    c = 340.0
    for M in (1.2, 1.6, 2.0, 3.7):
        fr = AC.source_wavefronts(M * c, c, 3.0, n_fronts=200, dt=0.01)
        d = fr["x_source"] - fr["x_centres"]                                           # distance from the source back to each centre
        assert np.all(d > 0) and close(np.arcsin(fr["radii"] / d), GD.mach_angle(M), 1e-12)     # the tangent from the source to every circle
        assert close(fr["mu"], math.asin(1 / M), 1e-14) and np.all(fr["ahead"] < 0) and fr["M"] == pytest.approx(M)
        assert close(fr["t_emit"], 3.0 - 0.01 * np.arange(1, 201), 1e-13) and close(fr["radii"], c * (3.0 - fr["t_emit"]), 1e-12) and close(fr["x_centres"], M * c * fr["t_emit"], 1e-13)
        # zone of action: a point is reached by some front  <=>  it lies inside the cone (x measured downstream from the source)
        rng = np.random.default_rng(3)
        xr, yr = rng.uniform(0.02 * c, 1.9 * c, 400), rng.uniform(-1.5 * c, 1.5 * c, 400)       # behind the source by xr, off the path by yr
        dense = AC.source_wavefronts(M * c, c, 3.0, n_fronts=3000, dt=1e-3)
        px = dense["x_source"] - xr
        inside_any = np.any((px[:, None] - dense["x_centres"][None, :]) ** 2 + yr[:, None] ** 2 <= dense["radii"][None, :] ** 2, axis=1)
        zone = AC.in_zone_of_action(xr, yr, M)
        edge = np.abs(xr - np.abs(yr) * math.sqrt(M * M - 1)) < 0.02 * c               # skip points within a front spacing of the cone
        far = xr * xr + yr * yr > (2.9 * c) ** 2
        assert np.array_equal(inside_any[~edge & ~far], zone[~edge & ~far])
        assert close(GD.mach_from_angle(GD.mach_angle(M)), M, 1e-13)
        h = I["h_flight"]                                                              # delay: the cone reaches the ground point after the overflight
        assert close(AC.mach_cone_delay(h, M, c), h / (M * c * math.tan(math.asin(1 / M))), 1e-13)
        assert AC.in_zone_of_action(M * c * AC.mach_cone_delay(h, M, c) * 1.001, h, M) and not AC.in_zone_of_action(M * c * AC.mach_cone_delay(h, M, c) * 0.999, h, M)
    sub = AC.source_wavefronts(0.5 * c, c, 3.0, n_fronts=5, dt=0.2)
    assert np.all(sub["ahead"] > 0) and math.isnan(sub["mu"])                          # a subsonic source is always inside its own fronts
    assert math.isnan(GD.mach_angle(0.99)) and close(GD.mach_angle(1.0), PI / 2, 1e-15) and math.isnan(AC.mach_cone_delay(1e4, 0.8, c))
    assert math.isnan(GD.mach_from_angle(0.0)) and math.isnan(GD.mach_from_angle(1.7)) and AC.in_zone_of_action(-5.0, 2.0, 0.6)
    # Doppler: fronts bunch ahead (spacing (c − U) dt) and spread behind ((c + U) dt)
    ahead_spacing = np.abs(np.diff(sub["x_centres"] + sub["radii"]))
    behind_spacing = np.abs(np.diff(sub["x_centres"] - sub["radii"]))
    assert close(ahead_spacing, (c - 0.5 * c) * 0.2, 1e-12) and close(behind_spacing, (c + 0.5 * c) * 0.2, 1e-12)
    assert close(AC.doppler_factor(0.5 * c, c, 0.0), c * 0.2 / ahead_spacing[0], 1e-12) and close(AC.doppler_factor(0.5 * c, c, PI), c * 0.2 / behind_spacing[0], 1e-12)


# =====================================================================================================================
# 3. §15.3–§15.5 — C03 the duct laws (15.18)–(15.26), C04 stagnation state (15.27)–(15.30), C05 area ratio (15.31),
#    C06 area–velocity (15.32)–(15.34) and the influence coefficients (with the three corrected slips)
# =====================================================================================================================
def test_duct_energy_equation_V2_derivation():  # V2 — D06 (two stars): (15.20) on the duct's control volume → (15.21) → (15.22); C03
    r1, r2, u1, u2, A1, p1, p2, e1, e2, Q = sp.symbols("rho_1 rho_2 u_1 u_2 A_1 p_1 p_2 e_1 e_2 Q", positive=True)
    mdot = r1 * u1 * A1                                                   # (15.18)
    A2 = mdot / (r2 * u2)
    # flux of (e + u^2/2) out of the control volume: the normal points against the flow on the inlet, with it on the outlet; walls pass nothing
    lhs = -r1 * (e1 + u1 ** 2 / 2) * u1 * A1 + r2 * (e2 + u2 ** 2 / 2) * u2 * A2
    assert z0(lhs - (-(e1 + u1 ** 2 / 2) * mdot + (e2 + u2 ** 2 / 2) * mdot))          # left side of (15.21)
    work = u1 * p1 * A1 - u2 * p2 * A2                                    # pressure work on the two end faces (tau_ij n_j = -p n_i there)
    eq21 = lhs - work - mdot * Q                                          # (15.21): ... = (upA)_1 - (upA)_2 + mdot Q
    h1, h2 = e1 + p1 / r1, e2 + p2 / r2
    assert z0(eq21 / mdot - (h2 + u2 ** 2 / 2 - h1 - u1 ** 2 / 2 - Q))                 # (15.22)
    assert not z0(eq21 / mdot - (e2 + u2 ** 2 / 2 - e1 - u1 ** 2 / 2 - Q))             # without the flow work it is not the energy law


def test_duct_mass_and_momentum_V2_derivation():  # V2 — D07 (two stars): (15.18) → (15.19); (15.24) → (15.25) → (15.26); C03
    x = sp.Symbol("x", real=True)
    rho, u, A, p = (sp.Function(n, positive=True)(x) for n in ("rho", "u", "A", "p"))
    Ff = sp.Function("F_f")(x)
    mdot = sp.Symbol("mdot", positive=True)
    assert z0(sp.diff(sp.log(rho * u * A), x) - (sp.diff(rho, x) / rho + sp.diff(u, x) / u + sp.diff(A, x) / A))      # (15.19)
    # (15.24) between x and x + dx: d(mdot u) = -d(pA) + dF, with dF = p dA (the sloping wall) - F_f dx (friction)
    eq24 = mdot * sp.diff(u, x) + sp.diff(p * A, x) - (p * sp.diff(A, x) - Ff)
    assert z0(eq24 - (mdot * sp.diff(u, x) + A * sp.diff(p, x) + Ff))                   # (15.25): mdot du/dx = -A dp/dx - F_f
    forgot = mdot * sp.diff(u, x) + sp.diff(p * A, x) + Ff                              # the wall's pressure force forgotten: a spurious p dA/dx
    assert z0(forgot - (mdot * sp.diff(u, x) + A * sp.diff(p, x) + Ff) - p * sp.diff(A, x))
    eq26 = (eq24.subs(Ff, 0).subs(mdot, rho * u * A)) / (rho * A)
    assert z0(eq26 - (u * sp.diff(u, x) + sp.diff(p, x) / rho))                         # (15.26)
    e = ch15.duct_equations_sympy()                                                     # the chapter's engine agrees, and shows the slip
    assert e["eq_15_19"] == 0 and e["eq_15_25"] == 0 and e["eq_15_26"] == 0 and e["eq_15_34"] == 0 and e["eq_15_33_printed"] != 0


def _dAdx(x, At=1.0, Ain=2.5, Ae=2.0, xt=0.4, L=1.0):
    return np.where(x <= xt, -2 * (Ain - At) * (xt - x) / xt ** 2, 2 * (Ae - At) * (x - xt) / (L - xt) ** 2)


def test_duct_laws_V4_mass_energy_entropy_along_an_isentropic_nozzle():  # V4 — (15.18), (15.22), (15.2) along GD.nozzle_isentropic; C03, C06
    x = np.linspace(0.0, 1.0, 501)                                        # the throat x = 0.4 is a grid point
    for g in (1.25, 1.4):
        for branch, A_star in (("subsonic", 0.8), ("subsonic", None), ("transonic", None), ("supersonic", None)):
            n = GD.nozzle_isentropic(x, GD.laval_area, g, branch=branch, A_star=A_star)
            rho, T, p, M, A = n["rho_over_rho0"], n["T_over_T0"], n["p_over_p0"], n["M"], n["A"]
            u = M * np.sqrt(T)                                            # in units of c0
            mdot = rho * u * A
            assert np.ptp(mdot) < 1e-9 * mdot[0], (branch, np.ptp(mdot))                 # (15.18)
            assert np.ptp(T + 0.5 * (g - 1) * u ** 2) < 1e-12                            # (15.22) with Q = 0: h + u^2/2 in units of h0
            assert np.ptp(np.log(p) - g * np.log(rho)) < 1e-11 and close(p, rho * T, 1e-12)      # entropy constant; p = rho R T
            assert n["x_throat"] == pytest.approx(0.4) and n["i_throat"] == 200 and close(n["A"], GD.laval_area(x), 1e-15)
            if branch == "transonic":
                assert M[200] == 1.0 and np.all(M[:200] < 1) and np.all(M[201:] > 1) and np.all(np.diff(M) > 0)    # accelerates all the way
                assert close(mdot[0], GD.mass_flux_function(1.0, g) / math.sqrt(g) * n["A_star"], 1e-10)          # sonic throat = choked flux
                slope = (M[201] - M[199]) / (x[201] - x[199])
                assert 1.0 < slope < 10.0                                                # passes M = 1 with a finite slope
            if branch == "subsonic" and A_star == 0.8:
                assert np.argmax(M) == 200 and M[200] < 1.0 and close(n["A_star"], 0.8)   # a venturi: fastest at the throat, not sonic
                s1 = dict(rho=rho[0], u=u[0], p=p[0] / g, h=T[0] / (g - 1))               # units: rho0, c0, rho0 c0^2 (p0 = rho0 c0^2 / gamma)
                s2 = dict(rho=rho[-1], u=u[-1], p=p[-1] / g, h=T[-1] / (g - 1))
                assert abs(GD.energy_residual(s1, s2, 0.0, g)) < 1e-12
                assert abs(GD.energy_residual({k: s1[k] for k in ("rho", "u", "p")}, {k: s2[k] for k in ("rho", "u", "p")}, 0.0, g)) < 1e-12
                assert abs(GD.energy_residual(s1, s2, 0.3, g) + 0.3) < 1e-12              # Q enters with the sign of (15.22)
                F = GD.wall_force(s1, s2, A[0], A[-1])                                    # (15.24) against the pressure integral over the wall
                pfun = lambda xx: float(GD.nozzle_isentropic(np.array([xx, 0.4]), GD.laval_area, g, branch="subsonic", A_star=0.8)["p_over_p0"][0]) / g   # noqa: E731
                integral = quad(lambda xx: pfun(xx) * float(_dAdx(xx)), 0.0, 1.0, points=[0.4], epsabs=1e-12, epsrel=1e-11, limit=200)[0]
                assert abs(F - integral) < 1e-8 * abs(p[0] / g * A[0]), (F, integral)
    with pytest.raises(ValueError):
        GD.nozzle_isentropic(x, GD.laval_area, branch="subsonic", A_star=1.01)           # choked upstream
    with pytest.raises(ValueError):
        GD.nozzle_isentropic(x, GD.laval_area, branch="sonic")
    sym = GD.nozzle_isentropic(np.linspace(0, 1, 201), lambda s: 1.0 + 3.0 * (s - 0.5) ** 2, branch="subsonic", A_star=0.7)
    assert close(sym["M"], sym["M"][::-1], 1e-12)                                        # mirror symmetry of a symmetric venturi
    assert GD.laval_area(0.0) == 2.5 and GD.laval_area(0.4) == 1.0 and GD.laval_area(1.0) == 2.0 and close(GD.laval_area(0.7), 1.25, 1e-14)
    assert abs((GD.laval_area(0.4 + 1e-6) - GD.laval_area(0.4 - 1e-6)) / 2e-6) < 1e-5     # zero slope at the throat
    assert close(GD.laval_area(0.2, A_throat=0.01, A_inlet=0.025, A_exit=0.02), 0.01375, 1e-13)


def test_stagnation_ratios_V1_from_the_energy_law_and_the_isentrope():  # V1 — D08, C04: (15.27) → (15.28); (15.2) → (15.29), (15.30)
    M = np.concatenate([M_SUB, M_SUP])
    for g in GAMMAS:
        gas = GD.PerfectGas(g, 301.0)
        T = 263.0
        u = M * gas.c(T)
        T0 = ch15.stagnation_temperature(T, u, cp=gas.cp)                 # (15.27), the ch04 function
        assert close(GD.stagnation_temperature_ratio(M, g), T0 / T, 1e-13)
        assert close(GD.stagnation_pressure_ratio(M, g), (T0 / T) ** (g / (g - 1)), 1e-13)           # (15.2) between static and stagnation
        assert close(GD.stagnation_density_ratio(M, g), (T0 / T) ** (1 / (g - 1)), 1e-13)
        assert close(GD.stagnation_pressure_ratio(M, g), GD.stagnation_density_ratio(M, g) * GD.stagnation_temperature_ratio(M, g), 1e-13)   # p0 = rho0 R T0
        p = 0.7e5
        for Mi, T0i in zip(M, T0):                                         # the stagnation state has the entropy of the static state
            dS = ch15.perfect_gas_entropy_change_p(T, p, T0i, p * GD.stagnation_pressure_ratio(Mi, g), cp=gas.cp, R=gas.R)
            assert abs(dS) < 1e-9 * gas.cp
        assert close(GD.mach_from_pressure_ratio(1 / GD.stagnation_pressure_ratio(M, g), g), M, 1e-12)
        s = GD.isentropic_state(M, g)
        assert close(s["T0_over_T"], T0 / T, 1e-13) and close(s["p_over_p0"] * s["p0_over_p"], 1.0, 1e-14) and close(s["rho_over_rho0"] * s["rho0_over_rho"], 1.0, 1e-14)
        assert close(s["A_over_Astar"], GD.area_ratio(M, g), 1e-15) and close(s["T_over_T0"], T / T0, 1e-13) and close(s["M"], M, 0)
        assert close(GD.total_enthalpy(T, u, gas.cp), gas.cp * T0, 1e-13)
        assert close(GD.max_speed(T0[0], gas.cp), math.sqrt(2 * gas.cp * T0[0]), 1e-14)
        Mbig = 1e4                                                         # u → u_max as M → infinity at fixed T0
        assert rel(Mbig * gas.c(500.0 / GD.stagnation_temperature_ratio(Mbig, g)), GD.max_speed(500.0, gas.cp)) < 1e-7
    assert math.isnan(GD.mach_from_pressure_ratio(1.2)) and math.isnan(GD.mach_from_pressure_ratio(0.0)) and GD.mach_from_pressure_ratio(1.0) == 0.0
    tab = GD.isentropic_table(M, 1.4)
    assert list(tab.columns) == ["M", "p/p0", "rho/rho0", "T/T0", "A/A*"] and close(tab["p/p0"].to_numpy(), 1 / GD.stagnation_pressure_ratio(M), 1e-14)
    # limits (V7): gamma → 1 gives exp(M^2/2); M → 0 gives the incompressible Bernoulli value p + rho u^2/2
    assert close(GD.stagnation_pressure_ratio(M_SUB, 1 + 1e-7), np.exp(M_SUB ** 2 / 2), 1e-6)
    small = 1e-3
    pb = ch15.stagnation_pressure(1.0, small * math.sqrt(1.4), 1.0)        # p = rho = 1, c^2 = gamma
    assert abs(GD.stagnation_pressure_ratio(small) - pb) < 1e-12 and close(pb, 1 + 0.7e-6, 1e-14)


def test_stagnation_ratios_V5_naca_report_1135_rows():  # V5 — NACA Report 1135, Tables I and II (gamma = 7/5), read first-hand; C04, C05
    b = bench()
    t1, t2 = b["naca1135_table1_subsonic"], b["naca1135_table2_supersonic"]
    assert t1["columns"][:8] == ["M", "p/pt", "rho/rhot", "T/Tt", "beta", "q/pt", "A/A*", "V/a*"] and len(t1["rows"]) >= 3 and len(t2["rows"]) >= 6
    misses, n = [], 0
    for tab in (t1, t2):
        for row in tab["rows"]:
            M = row[0]
            s = GD.isentropic_state(M, 1.4)
            ours = [M, s["p_over_p0"], s["rho_over_rho0"], s["T_over_T0"], math.sqrt(abs(M * M - 1)), 0.5 * 1.4 * M * M * s["p_over_p0"],
                    s["A_over_Astar"], M * math.sqrt(s["T_over_T0"] * 1.2)]
            for k in range(1, 8):
                n += 1
                half = tol_digits(row[k], tab["decimals"][k])
                assert abs(ours[k] - row[k]) <= 2 * half * 1.001, (M, tab["columns"][k], ours[k], row[k])      # one unit of the last printed digit
                if abs(ours[k] - row[k]) > half * 1.001:
                    misses.append((M, tab["columns"][k], ours[k], row[k]))
    # 62 of 63 entries agree to half a unit of the last digit; the one that does not sits 2.4e-6 from the rounding boundary
    assert n == 63 and [(m[0], m[1]) for m in misses] == [(0.55, "A/A*")] and abs(misses[0][2] - 1.25495) < 3e-6, misses
    for tab in (t1, t2):
        for row in tab["rows"]:
            M = row[0]
            if M < 1:
                assert rel(GD.mach_from_area_ratio(row[6], branch="subsonic"), M) < 2e-4
                assert rel(GD.mach_from_pressure_ratio(row[1]), M) < 2e-4
            else:
                assert rel(GD.mach_from_area_ratio(row[6], branch="supersonic"), M) < 5e-4


def test_pitot_compressibility_V2_series_and_incompressible_limit():  # V2 — (15.29) expanded: p0 − p = (rho u^2/2)(1 + M^2/4 + (2 − gamma) M^4/24 + …); N32, N33
    g, M = sp.symbols("gamma M", positive=True)
    exact = ((1 + (g - 1) / 2 * M ** 2) ** (g / (g - 1)) - 1) / (g * M ** 2 / 2)
    ser = sp.series(exact, M, 0, 6).removeO()
    assert z0(ser - (1 + M ** 2 / 4 + (2 - g) * M ** 4 / 24))
    dens = sp.series((1 + (g - 1) / 2 * M ** 2) ** (1 / (g - 1)), M, 0, 4).removeO()
    assert z0(dens - (1 + M ** 2 / 2))                                     # (15.30): the density changes by M^2/2 to leading order
    for gv in GAMMAS:
        Ms = np.array([1e-4, 1e-2, 0.1])
        assert np.all(np.abs(GD.pitot_compressibility_factor(Ms, gv) - (1 + Ms ** 2 / 4 + (2 - gv) * Ms ** 4 / 24)) < 0.1 * Ms ** 6 + 1e-15)   # next term is O(M^6)
        Mb = np.array([0.3, 0.7, 0.95])
        direct = (GD.stagnation_pressure_ratio(Mb, gv) - 1) / (0.5 * gv * Mb ** 2)
        assert close(GD.pitot_compressibility_factor(Mb, gv), direct, 1e-12)
        fr = np.array([1e-6, 0.01, 0.05, 0.4])
        Mf = GD.mach_for_density_change(fr, gv)
        assert close(GD.stagnation_density_ratio(Mf, gv) - 1, fr, 1e-9) and rel(Mf[0], math.sqrt(2e-6)) < 1e-5
    assert GD.pitot_compressibility_factor(0.0) == 1.0 and abs(GD.pitot_compressibility_factor(0.3) - 1.0227) < 5e-5
    assert abs(GD.pitot_compressibility_factor(1e-7) - 1.0) < 1e-13         # no cancellation at tiny M


def test_entropy_from_stagnation_states_V1_matches_the_state_function():  # V1 — D09, N34: S2 − S1 = −R ln(p02/p01) + Cp ln(T02/T01); T17
    gas = GD.PerfectGas(1.33, 296.0)
    for pr, Tr in ((0.72, 1.0), (0.95, 1.3), (1.4, 0.8)):
        ref = ch15.perfect_gas_entropy_change_p(410.0, 3.0e5, 410.0 * Tr, 3.0e5 * pr, cp=gas.cp, R=gas.R)
        assert close(GD.entropy_change_from_stagnation(pr, Tr, gas), ref, 1e-12)
    dS = GD.entropy_change_from_stagnation(np.array([0.5, 0.9, 1.0]), 1.0, gas)
    assert np.all(dS >= 0) and close(GD.stagnation_pressure_loss(dS, gas), [0.5, 0.9, 1.0], 1e-13)      # T0 kept, S up  =>  p0 down
    for M1 in M_SUP:                                                       # across a shock it is (15.43) times Cv  (V4: two routes to one state function)
        s = GD.normal_shock(M1, gas.gamma)
        assert close(GD.entropy_change_from_stagnation(s["p0_ratio"], 1.0, gas), s["dS_over_cv"] * gas.cv, 1e-9, 1e-13)
        assert close(s["dS_over_R"], s["dS_over_cv"] / (gas.gamma - 1), 1e-13)


def test_area_ratio_V1_built_from_mass_conservation_and_its_two_inverses():  # V1 — C05, D10 numerically: A/A* = rho* c* / (rho u); trap T8
    for g in GAMMAS:
        M = np.concatenate([M_SUB, M_SUP])
        rho, T = 1 / GD.stagnation_density_ratio(M, g), 1 / GD.stagnation_temperature_ratio(M, g)
        rs, Ts = 1 / GD.stagnation_density_ratio(1.0, g), 1 / GD.stagnation_temperature_ratio(1.0, g)
        assert close(GD.area_ratio(M, g), rs * math.sqrt(Ts) / (rho * M * np.sqrt(T)), 1e-13)
        assert GD.area_ratio(1.0, g) == pytest.approx(1.0, abs=1e-15) and np.all(GD.area_ratio(M, g) > 1.0)
        assert close(GD.mach_from_area_ratio(GD.area_ratio(M_SUB, g), g, branch="subsonic"), M_SUB, 1e-9)
        assert close(GD.mach_from_area_ratio(GD.area_ratio(M_SUP, g), g, branch="supersonic"), M_SUP, 1e-9)
        assert close(GD.mass_flux_function(M, g) * GD.area_ratio(M, g), GD.mass_flux_function(1.0, g), 1e-13)     # same mdot through A and A*
        sr = GD.sonic_ratios(g)
        st = GD.isentropic_state(1.0, g)
        assert close([sr["T"], sr["p"], sr["rho"]], [st["T_over_T0"], st["p_over_p0"], st["rho_over_rho0"]], 1e-14) and close(sr["c"], math.sqrt(sr["T"]), 1e-15)
        assert close(GD.critical_pressure_ratio(g), sr["p"], 1e-15)
    ar = 1.7
    Ms, Mp = GD.mach_from_area_ratio(ar, branch="subsonic"), GD.mach_from_area_ratio(ar, branch="supersonic")
    assert Ms < 1 < Mp and close(GD.area_ratio(Ms), ar, 1e-12) and close(GD.area_ratio(Mp), ar, 1e-12)             # two Mach numbers share one area
    assert GD.mach_from_area_ratio(1.0, branch="subsonic") == 1.0 and GD.mach_from_area_ratio(1.0, branch="supersonic") == 1.0
    eps = 1e-10                                                            # the minimum is quadratic: M − 1 ≈ ± sqrt(2 eps/(d2A/dM2)) with d2/dM2 = 4/(gamma + 1)
    assert rel(GD.mach_from_area_ratio(1 + eps, branch="supersonic") - 1, math.sqrt(eps * 1.2)) < 1e-3
    assert rel(1 - GD.mach_from_area_ratio(1 + eps, branch="subsonic"), math.sqrt(eps * 1.2)) < 1e-3
    assert rel(GD.mach_from_area_ratio(1e6, branch="subsonic"), GD.mass_flux_function(1.0) / math.sqrt(1.4) / 1e6) < 1e-6      # M → flux/(sqrt(gamma) A/A*) as M → 0
    for bad in (0.99, 0.0, -1.0):
        with pytest.raises(ValueError):
            GD.mach_from_area_ratio(bad, branch="subsonic")
    with pytest.raises(ValueError):
        GD.mach_from_area_ratio(2.0, branch="transonic")
    assert math.isinf(GD.area_ratio(0.0))


def test_area_ratio_V2_derivation():  # V2 — D10 (two stars): the chain (rho*/rho0)(rho0/rho)(c*/c)(1/M) → (15.31), minimum 1 at M = 1; C05
    g, M = sp.symbols("gamma M", positive=True)
    Tf = 1 + (g - 1) / 2 * M ** 2                                          # T0/T (15.28)
    chain = (2 / (g + 1)) ** (1 / (g - 1)) * Tf ** (1 / (g - 1)) * sp.sqrt(2 / (g + 1) * Tf) / M      # rho*/rho0 · rho0/rho · sqrt(T*/T) / M
    printed = (2 / (g + 1) * Tf) ** (sp.Rational(1, 2) * (g + 1) / (g - 1)) / M                        # (15.31)
    assert sp.simplify(1 / (g - 1) + sp.Rational(1, 2) - (g + 1) / (2 * (g - 1))) == 0                 # the exponents add
    for gv in (sp.Rational(7, 5), sp.Rational(5, 3), sp.Rational(9, 8)):
        for mv in (sp.Rational(1, 3), sp.Rational(3, 2), 4):
            assert abs(sp.N((chain / printed).subs({g: gv, M: mv}), 30) - 1) < sp.Float("1e-25")
    dlog = sp.simplify(sp.diff(sp.log(printed), M))
    assert z0(dlog - (M ** 2 - 1) / (M * Tf))                              # d ln(A/A*)/dM = (M^2 − 1)/(M (1 + (gamma−1)M^2/2)): zero only at M = 1
    assert sp.simplify(printed.subs(M, 1)) == 1 and sp.simplify(sp.diff(dlog, M).subs(M, 1) - 4 / (g + 1)) == 0    # a minimum (curvature 4/(gamma+1) > 0)
    # the same sign change is the area–velocity relation (15.34): d ln A / d ln u at constant p0, T0
    dlnu = sp.diff(sp.log(M / sp.sqrt(Tf)), M)
    assert z0(dlog / dlnu + (1 - M ** 2))
    fn = sp.lambdify((M, g), printed, "numpy")                               # the code is this expression
    for gv in GAMMAS:
        Mv = np.concatenate([M_SUB, [1.0], M_SUP])
        assert close(GD.area_ratio(Mv, gv), fn(Mv, gv), 1e-13)


def test_sonic_ratios_and_choked_flow_V2_derivation():  # V2 — D11 (two stars): (15.28)–(15.30) at M = 1 and mdot_max = rho* c* A*; C05, C09
    g, R, T0, p0, As, M = sp.symbols("gamma R T_0 p_0 A_s M", positive=True)
    Tf = 1 + (g - 1) / 2 * M ** 2
    assert sp.simplify(1 / Tf.subs(M, 1) - 2 / (g + 1)) == 0
    rho_s = p0 / (R * T0) * (2 / (g + 1)) ** (1 / (g - 1))                 # rho* = rho0 (2/(gamma+1))^{1/(gamma−1)}, rho0 = p0/(R T0)
    c_s = sp.sqrt(g * R * T0 * 2 / (g + 1))                                # c* = sqrt(gamma R T*)
    formula = p0 * As * sp.sqrt(g / (R * T0)) * (2 / (g + 1)) ** ((g + 1) / (2 * (g - 1)))
    for gv in (sp.Rational(7, 5), sp.Rational(5, 3), sp.Rational(5, 4)):
        assert sp.simplify((rho_s * c_s * As / formula).subs(g, gv)) == 1
    assert sp.simplify(((g + 1) / (2 * (g - 1))).subs(g, sp.Rational(7, 5))) == 3      # the exponent is 3 for gamma = 1.4
    flux = sp.sqrt(g) * M * Tf ** (-(g + 1) / (2 * (g - 1)))               # mdot sqrt(R T0)/(p0 A)
    d = sp.simplify(sp.diff(sp.log(flux), M))
    assert z0(d - (1 - M ** 2) / (M * Tf))                                 # the flux is largest at M = 1
    for gv in GAMMAS:
        gas = GD.PerfectGas(gv, 310.0)
        want = 4.0e5 * 0.013 * math.sqrt(gv / (310.0 * 450.0)) * (2 / (gv + 1)) ** ((gv + 1) / (2 * (gv - 1)))
        assert close(GD.choked_mass_flow(4.0e5, 450.0, 0.013, gas), want, 1e-13)
        assert close(GD.mass_flow_from_mach(1.0, 4.0e5, 450.0, 0.013, gas), want, 1e-13)
        Mg = np.linspace(0.05, 4, 400)
        assert abs(Mg[np.argmax(GD.mass_flux_function(Mg, gv))] - 1.0) < 0.011 and np.all(GD.mass_flux_function(Mg, gv) <= GD.mass_flux_function(1.0, gv) * (1 + 1e-15))
        s = GD.isentropic_state(0.63, gv)                                  # mass_flow_from_mach = rho u A from the static state
        rho, T = 4.0e5 / (310.0 * 450.0) * s["rho_over_rho0"], 450.0 * s["T_over_T0"]
        assert close(GD.mass_flow_from_mach(0.63, 4.0e5, 450.0, 0.013, gas), GD.mass_flow(rho, 0.63 * gas.c(T), 0.013), 1e-13)
    assert abs(GD.mass_flux_function(1.0, 1.4) - 0.6847) < 5e-5


def test_choked_flow_V5_cited_critical_ratio_and_gas_table():  # V5 — Wikipedia "Choked flow" (secondary): 0.528 and the table of minimum pressure ratios; C09
    b = bench()["choked_flow"]
    assert abs(GD.critical_pressure_ratio(1.4) - b["critical_ratio_air"]) < 5e-4
    for name, gam, ratio in b["min_upstream_over_downstream"]:
        assert abs(1 / GD.critical_pressure_ratio(gam) - ratio) < 6e-4 * ratio, (name, 1 / GD.critical_pressure_ratio(gam), ratio)
    # the page's mass-flow form with C_d = 1: mdot = A sqrt(gamma rho0 P0 (2/(gamma+1))^((gamma+1)/(gamma−1)))
    gas = GD.PerfectGas(1.31, 518.0)
    P0, T0, A = 2.2e5, 330.0, 3.0e-4
    rho0 = P0 / (gas.R * T0)
    assert close(GD.choked_mass_flow(P0, T0, A, gas), A * math.sqrt(1.31 * rho0 * P0 * (2 / 2.31) ** (2.31 / 0.31)), 1e-12)


def _linear_duct_solution():
    """The four linearised laws solved by sympy (independent of the coded 4 x 4 solve)."""
    g, M, a, f, q = sp.symbols("gamma M a f q", real=True)
    dr, du, dp, dT = sp.symbols("dr du dp dT", real=True)
    sol = sp.solve([dr + du + a, dp + g * M ** 2 * du + f, dT + (g - 1) * M ** 2 * du - q, dp - dr - dT], [dr, du, dp, dT], dict=True)[0]
    Tf = 1 + (g - 1) / 2 * M ** 2
    dM = sol[du] - sol[dT] / 2
    out = dict(drho_over_rho=sol[dr], du_over_u=sol[du], dp_over_p=sol[dp], dT_over_T=sol[dT], dM_over_M=dM,
               ds_over_cp=sol[dT] - (g - 1) / g * sol[dp], u_du_over_h=(g - 1) * M ** 2 * sol[du],
               dT0_over_T0=sol[dT] + (g - 1) * M ** 2 * dM / Tf, dp0_over_p0=sol[dp] + g * M ** 2 * dM / Tf)
    return (g, M, a, f, q), out


def test_influence_coefficients_V2_sympy_solution_of_the_four_linear_laws():  # V2 — (15.19), (15.33), (15.34), (15.47), (15.48) corrected; C06, C12
    (g, M, a, f, q), sol = _linear_duct_solution()
    assert z0(sol["du_over_u"] - (-a + f + q) / (1 - M ** 2))                              # (15.34) + second relation of (15.47) + corrected (15.48)
    assert z0(sol["drho_over_rho"].subs({f: 0, q: 0}) + M ** 2 * sol["du_over_u"].subs({f: 0, q: 0}))    # (15.33): d rho/rho = −M^2 du/u
    assert z0(-sol["dp_over_p"].subs({a: 0, q: 0}) - (1 + (g - 1) * M ** 2) / (1 - M ** 2) * f)          # (15.47) first relation, numerator corrected
    assert z0(sol["dT_over_T"].subs({a: 0, f: 0}) - (1 - g * M ** 2) / (1 - M ** 2) * q)                 # (15.48) first relation (as printed)
    assert z0(sol["u_du_over_h"].subs({a: 0, f: 0}) - (g - 1) * M ** 2 / (1 - M ** 2) * q)               # (15.48) second relation, sign corrected
    assert z0(sol["ds_over_cp"].subs({f: 0, q: 0})) and z0(sol["dT0_over_T0"].subs(q, 0)) and z0(sol["dp0_over_p0"].subs({f: 0, q: 0}))   # area change is isentropic
    assert z0(sol["ds_over_cp"] - ((g - 1) / g * f + q))                                   # T ds = dq + friction heat: entropy rises with f > 0
    assert z0(sol["dp0_over_p0"].subs({a: 0, q: 0}) + f) and z0(sol["dT0_over_T0"] - q / (1 + (g - 1) / 2 * M ** 2))   # friction costs p0 one for one
    fn = {k: sp.lambdify((g, M, a, f, q), v, "numpy") for k, v in sol.items()}
    for gv in GAMMAS:
        for Mv in (0.0, 0.3, 0.8, 1.2, 2.5):
            for kw in (dict(dA_over_A=0.013), dict(df=0.021), dict(dq=-0.017), dict(dA_over_A=-0.02, df=0.01, dq=0.03)):
                c = GD.influence_coefficients(Mv, gv, **kw)
                for k in fn:
                    want = float(fn[k](gv, Mv, kw.get("dA_over_A", 0.0), kw.get("df", 0.0), kw.get("dq", 0.0)))
                    assert abs(c[k] - want) < 1e-12 * max(1.0, abs(want)), (gv, Mv, kw, k, c[k], want)
    # each printed form fails the same symbolic check
    for Mv in (0.5, 0.9, 2.0):
        ok = GD.influence_coefficients(Mv, G, dA_over_A=1.0, df=1.0, dq=1.0)
        p32 = GD.influence_coefficients(Mv, G, dA_over_A=1.0, printed_15_32=True)
        p47 = GD.influence_coefficients(Mv, G, df=1.0, printed_15_47=True)
        p48 = GD.influence_coefficients(Mv, G, dq=1.0, printed_15_48=True)
        assert abs(p32["du_over_u"] - float(fn["du_over_u"](G, Mv, 1.0, 0.0, 0.0))) > 0.1 * abs(p32["du_over_u"])
        assert p32["drho_over_rho"] / p32["du_over_u"] == pytest.approx(Mv ** 2)            # the printed sign: d rho/rho = +M^2 du/u
        assert abs(p47["dp_over_p"] - float(fn["dp_over_p"](G, Mv, 0.0, 1.0, 0.0))) > 0.15 * abs(float(fn["dp_over_p"](G, Mv, 0.0, 1.0, 0.0)))
        assert p48["u_du_over_h"] == pytest.approx(-float(fn["u_du_over_h"](G, Mv, 0.0, 0.0, 1.0)))     # the printed sign of (15.48)
        assert all(np.isfinite(v) for v in ok.values())
    assert close(GD.area_velocity(0.5, 0.01), GD.influence_coefficients(0.5, dA_over_A=0.01)["du_over_u"], 1e-15)


def _fd_area(M, g, key, h=1e-5):
    br = "subsonic" if M < 1 else "supersonic"
    A0 = GD.area_ratio(M, g)

    def logs(scale):
        Mx = GD.mach_from_area_ratio(A0 * scale, g, branch=br)
        s = GD.isentropic_state(Mx, g)
        return dict(drho_over_rho=math.log(s["rho_over_rho0"]), du_over_u=math.log(Mx * math.sqrt(s["T_over_T0"])), dp_over_p=math.log(s["p_over_p0"]),
                    dT_over_T=math.log(s["T_over_T0"]), dM_over_M=math.log(Mx))[key]
    return (logs(1 + h) - logs(1 - h)) / (math.log(1 + h) - math.log(1 - h))


def test_influence_area_V1_finite_differences_of_15_31_and_printed_15_32_fails():  # V1 — C06 against (15.31); slip #4 (printed (15.32)) must fail
    for g in (1.4, 5.0 / 3.0):
        for M in (0.3, 0.5, 0.9, 1.3, 2.0, 3.5):
            c = GD.influence_coefficients(M, g, dA_over_A=1.0)
            pr = GD.influence_coefficients(M, g, dA_over_A=1.0, printed_15_32=True)
            for key in ("drho_over_rho", "du_over_u", "dp_over_p", "dT_over_T", "dM_over_M"):
                fd = _fd_area(M, g, key)
                assert abs(c[key] - fd) < 1e-6 * max(1.0, abs(fd)), (g, M, key, c[key], fd)
            fd_u, fd_r = _fd_area(M, g, "du_over_u"), _fd_area(M, g, "drho_over_rho")
            assert abs(pr["du_over_u"] - fd_u) > 0.15 * abs(fd_u), (M, pr["du_over_u"], fd_u)            # the printed form misses by 18 % and more
            assert abs(pr["drho_over_rho"] - fd_r) > 0.15 * abs(fd_r)
            assert close(GD.area_velocity(M, 1.0), -1 / (1 - M * M), 1e-14) and close(c["du_over_u"], GD.area_velocity(M, 1.0), 1e-14)
    assert np.sign(GD.influence_coefficients(2.0, dA_over_A=1.0, printed_15_32=True)["du_over_u"]) == -1 and GD.area_velocity(2.0, 1.0) > 0   # printed: a widening duct would slow a supersonic stream


def _fd_duct(M, g, key, which, h=2e-4):
    def one(hh):
        vals = []
        for sgn in (+1, -1):
            kw = {which: sgn * hh}
            M2 = GD.duct_outlet_mach(M, gamma=g, **kw)
            s = GD.duct_state_ratios(M, M2, gamma=g, **kw)
            vals.append(dict(drho_over_rho=math.log(s["rho_ratio"]), du_over_u=math.log(s["u_ratio"]), dp_over_p=math.log(s["p_ratio"]),
                             dT_over_T=math.log(s["T_ratio"]), dM_over_M=math.log(M2 / M), dp0_over_p0=math.log(s["p0_ratio"]),
                             dT0_over_T0=math.log(s["T0_ratio"]), ds_over_cp=s["dS_over_cp"])[key])
        return (vals[0] - vals[1]) / (2 * hh)
    return (4 * one(h / 2) - one(h)) / 3                                   # Richardson: O(h^4)


def test_influence_friction_heat_V1_finite_differences_of_15_46_and_printed_15_47_15_48_fail():  # V1 — C12 against (15.46); slips #11, #12 must fail
    keys = ("drho_over_rho", "du_over_u", "dp_over_p", "dT_over_T", "dM_over_M", "dp0_over_p0", "dT0_over_T0", "ds_over_cp")
    for g in (1.4, 1.25):
        for M in (0.3, 0.5, 0.9, 1.5, 2.0, 3.0):
            for which, kw in (("f", dict(df=1.0)), ("q", dict(dq=1.0))):
                c = GD.influence_coefficients(M, g, **kw)
                for key in keys:
                    fd = _fd_duct(M, g, key, which)
                    assert abs(c[key] - fd) < 1e-6 * max(1.0, abs(fd)), (g, M, which, key, c[key], fd)
            fd_p = _fd_duct(M, g, "dp_over_p", "f")
            p47 = GD.influence_coefficients(M, g, df=1.0, printed_15_47=True)["dp_over_p"]
            miss = 2 * (g - 1) * M * M / (1 + (g - 1) * M * M)                                          # printed (15.47): wrong numerator — it misses by exactly this fraction
            assert abs(abs(p47 - fd_p) / abs(fd_p) - miss) < 1e-5 and miss > 0.03, (M, p47, fd_p, miss)
            fd_u = _fd_duct(M, g, "du_over_u", "q")
            p48 = GD.influence_coefficients(M, g, dq=1.0, printed_15_48=True)
            assert np.sign(p48["du_over_u"]) == -np.sign(fd_u) and abs(p48["u_du_over_h"] + (g - 1) * M * M * fd_u) < 1e-5 * abs(fd_u) + 1e-7   # printed (15.48): wrong sign
    fd, pr = _fd_duct(2.0, G, "dp_over_p", "f"), GD.influence_coefficients(2.0, G, df=1.0, printed_15_47=True)["dp_over_p"]
    assert fd > 0 > pr                                                     # at M = 2 the printed form even has the wrong sign
    assert _fd_duct(0.5, G, "du_over_u", "q") > 0                           # heating accelerates a subsonic stream (the book's own text)


def test_friction_differentials_V2_derivation():  # V2 — D22 (two stars): differentials of (15.45) with q = 0 → (15.47), numerator corrected; C12
    g, M, eps, F = sp.symbols("gamma M epsilon F", positive=True)
    a, b, c = sp.symbols("a b c", real=True)                              # p2 = p1(1 + eps a), rho2 = rho1(1 + eps b), u2 = u1(1 + eps c)
    p1, r1 = sp.symbols("p_1 rho_1", positive=True)
    u1 = M * sp.sqrt(g * p1 / r1)
    p2, r2, u2 = p1 * (1 + eps * a), r1 * (1 + eps * b), u1 * (1 + eps * c)
    h = lambda p, r: g / (g - 1) * p / r                                   # noqa: E731
    laws = [r2 * u2 - r1 * u1, p2 + r2 * u2 ** 2 + p1 * eps * F - p1 - r1 * u1 ** 2, h(p2, r2) + u2 ** 2 / 2 - h(p1, r1) - u1 ** 2 / 2]     # (15.45), f = eps F, q = 0
    lin = [sp.simplify(sp.diff(e, eps).subs(eps, 0)) for e in laws]
    sol = sp.solve(lin, [a, b, c], dict=True)[0]
    assert z0(sol[c] - F / (1 - M ** 2)) and z0(sol[b] + sol[c])                           # du/u = −d rho/rho = df/(1 − M^2)   (second relation)
    assert z0(-sol[a] - (1 + (g - 1) * M ** 2) / (1 - M ** 2) * F)                         # −dp/p1 = [1 + (gamma−1)M^2]/(1 − M^2) df   (corrected)
    assert not z0(-sol[a] - (1 - (g - 1) * M ** 2) / (1 - M ** 2) * F)                     # the printed numerator 1 − (gamma−1)M^2 does not follow
    # the intermediate lines of Part F: energy gives dT/T = −(gamma−1)M^2 du/u; the gas law then dp/p = −[1 + (gamma−1)M^2] du/u
    dT = sol[a] - sol[b]
    assert z0(dT + (g - 1) * M ** 2 * sol[c]) and z0(sol[a] + (1 + (g - 1) * M ** 2) * sol[c])


def test_heating_differentials_V2_derivation():  # V2 — D23 (two stars): differentials of (15.45) with f = 0 → (15.48), second relation corrected; C12
    g, M, eps, Qh = sp.symbols("gamma M epsilon Q_h", positive=True)
    a, b, c = sp.symbols("a b c", real=True)
    p1, r1 = sp.symbols("p_1 rho_1", positive=True)
    u1 = M * sp.sqrt(g * p1 / r1)
    p2, r2, u2 = p1 * (1 + eps * a), r1 * (1 + eps * b), u1 * (1 + eps * c)
    h = lambda p, r: g / (g - 1) * p / r                                   # noqa: E731
    laws = [r2 * u2 - r1 * u1, p2 + r2 * u2 ** 2 - p1 - r1 * u1 ** 2, h(p2, r2) + u2 ** 2 / 2 - h(p1, r1) * (1 + eps * Qh) - u1 ** 2 / 2]    # (15.45), q = eps Q_h
    sol = sp.solve([sp.simplify(sp.diff(e, eps).subs(eps, 0)) for e in laws], [a, b, c], dict=True)[0]
    dT = sol[a] - sol[b]
    assert z0(sol[a] + g * M ** 2 * sol[c]) and z0(dT - (1 - g * M ** 2) * sol[c])         # momentum: dp/p = −gamma M^2 du/u, so dT/T = (1 − gamma M^2) du/u
    assert z0(sol[c] - Qh / (1 - M ** 2))                                                  # energy: du/u = dq/(1 − M^2)
    assert z0(dT - (1 - g * M ** 2) / (1 - M ** 2) * Qh)                                   # (15.48) first relation
    udu_over_h = (g - 1) * M ** 2 * sol[c]                                                 # u^2 = (gamma − 1) M^2 h
    assert z0(udu_over_h - (g - 1) * M ** 2 / (1 - M ** 2) * Qh) and not z0(udu_over_h - (1 - g) * M ** 2 / (1 - M ** 2) * Qh)   # corrected sign; printed sign fails
    assert sp.solve(sp.Eq((1 - g * M ** 2), 0), M) == [1 / sp.sqrt(g)]                     # dT changes sign at M^2 = 1/gamma, not at M = 1 (T18)


def test_influence_coefficients_V7_limits_and_sign_pattern():  # V7 — the four consequences of (15.34) and the two "paradoxes"; C06, C12, T18
    c0 = GD.influence_coefficients(0.0, dA_over_A=0.02)
    assert c0["du_over_u"] == pytest.approx(-0.02) and c0["drho_over_rho"] == 0.0 and c0["dp_over_p"] == 0.0     # M → 0: du/u = −dA/A, constant density (slip #5)
    for M in (0.2, 0.7, 0.95):
        assert GD.area_velocity(M, -1.0) > 0 > GD.area_velocity(M, 1.0)                    # subsonic: converging speeds up
        assert GD.duct_role(M, -1) == "nozzle" and GD.duct_role(M, +1) == "diffuser"
        assert GD.influence_coefficients(M, df=1.0)["dM_over_M"] > 0 > GD.influence_coefficients(M, df=1.0)["dp_over_p"]     # friction speeds a subsonic flow up
        assert GD.influence_coefficients(M, dq=1.0)["dM_over_M"] > 0
    for M in (1.05, 1.8, 4.0):
        assert GD.area_velocity(M, 1.0) > 0 > GD.area_velocity(M, -1.0)                    # supersonic: diverging speeds up
        assert GD.duct_role(M, +1) == "nozzle" and GD.duct_role(M, -1) == "diffuser"
        assert GD.influence_coefficients(M, df=1.0)["dM_over_M"] < 0 < GD.influence_coefficients(M, df=1.0)["dp_over_p"]
        assert GD.influence_coefficients(M, dq=1.0)["dM_over_M"] < 0
    assert GD.duct_role(0.5, 0) == "none" and GD.duct_role(1.0, 1) == "sonic"
    for g in GAMMAS:
        lo, hi = 1 / math.sqrt(g), 1.0                                     # heating cools the gas only for 1/gamma < M^2 < 1
        for M, sign in ((0.9 * lo, 1), (0.5 * (lo + hi), -1), (1.3, 1)):
            assert np.sign(GD.influence_coefficients(M, g, dq=1.0)["dT_over_T"]) == sign and GD.rayleigh_temperature_sign(M, g) == sign
        assert GD.rayleigh_temperature_sign(lo, g) == 0.0 or abs(1 - g * lo * lo) < 1e-15
        assert math.isnan(GD.rayleigh_temperature_sign(1.0, g))
    for M in (0.4, 2.2):
        cf, cq = GD.influence_coefficients(M, df=0.01), GD.influence_coefficients(M, dq=0.01)
        assert cf["ds_over_cp"] > 0 and cf["dp0_over_p0"] < 0 and abs(cf["dT0_over_T0"]) < 1e-15 and cq["ds_over_cp"] > 0 and cq["dT0_over_T0"] > 0 and cq["dp0_over_p0"] < 0
    lin = GD.influence_coefficients(np.array([0.5, 2.0]), dA_over_A=np.array([0.01, 0.02]))
    assert lin["du_over_u"].shape == (2,) and close(lin["du_over_u"], [-0.01 / 0.75, 0.02 / 3.0], 1e-13)


def test_influence_coefficients_V7_sonic_point_is_infinite_as_documented():  # V7 — docstring: "±inf at M = 1 unless −dA/A + df + dq vanishes (then NaN)"
    assert math.isinf(GD.area_velocity(1.0, 0.1))                           # (15.34) itself: infinite
    assert math.isnan(GD.influence_coefficients(1.0, G)["du_over_u"])       # 0/0: nothing drives the flow
    for kw in (dict(dA_over_A=0.1), dict(df=0.1), dict(dq=-0.1), dict(dA_over_A=0.1, df=0.03)):
        v = GD.influence_coefficients(1.0, G, **kw)["du_over_u"]
        assert math.isinf(v), (kw, v)                                       # one non-zero driver alone must give ±inf, not NaN


SINGULAR_KEYS = ("drho_over_rho", "du_over_u", "dp_over_p", "dT_over_T", "dM_over_M", "u_du_over_h")      # carry 1/(1 − M^2)
REGULAR_KEYS = ("ds_over_cp", "dT0_over_T0", "dp0_over_p0")                                               # do not


def test_influence_entropy_and_stagnation_changes_V2_derivation():  # V2 — C06, C12: ds/Cp, dT0/T0, dp0/p0 from the finite laws of a short duct ((15.18), (15.24), (15.45)); no 1/(1 − M^2) in them
    g, M, eps = sp.symbols("gamma M epsilon", positive=True)
    al, F, Qh = sp.symbols("alpha F Q_h", real=True)                         # A2 = A1(1 + eps alpha), friction eps F p1 A1, heat eps Q_h h1
    a, b, c = sp.symbols("a b c", real=True)
    p1, r1, A1 = sp.symbols("p_1 rho_1 A_1", positive=True)
    u1 = M * sp.sqrt(g * p1 / r1)
    p2, r2, u2, A2 = p1 * (1 + eps * a), r1 * (1 + eps * b), u1 * (1 + eps * c), A1 * (1 + eps * al)
    h = lambda p, r: g / (g - 1) * p / r                                     # noqa: E731
    laws = [r2 * u2 * A2 - r1 * u1 * A1,                                                                           # mass (15.18)
            (p2 + r2 * u2 ** 2) * A2 - (p1 + r1 * u1 ** 2) * A1 - (p1 + p2) / 2 * (A2 - A1) + p1 * A1 * eps * F,   # momentum (15.24) with the wall's pressure force and the friction force
            h(p2, r2) + u2 ** 2 / 2 - h(p1, r1) * (1 + eps * Qh) - u1 ** 2 / 2]                                    # energy (15.22)
    sol = sp.solve([sp.simplify(sp.diff(e, eps).subs(eps, 0)) for e in laws], [a, b, c], dict=True)[0]
    assert z0(sol[c] - (-al + F + Qh) / (1 - M ** 2))                        # the driver −dA/A + df + dq over 1 − M^2: this key has the pole
    assert z0(sp.limit((1 - M ** 2) * sol[c], M, 1) - (-al + F + Qh))
    # the three state functions of the outlet state, written from their definitions (not from the linear system of the code)
    Msq2 = u2 ** 2 * r2 / (g * p2)
    Tfac2, Tf = 1 + (g - 1) / 2 * Msq2, 1 + (g - 1) / 2 * M ** 2
    lin = lambda expr: sp.simplify(sp.diff(expr, eps).subs(eps, 0).subs(sol))    # noqa: E731
    dS = lin(sp.log(p2 / r2) - (g - 1) / g * sp.log(p2))                     # s/Cp = ln T − [(gamma − 1)/gamma] ln p
    dT0 = lin(sp.log(p2 / r2 * Tfac2))                                       # T0 = T (1 + (gamma − 1) M^2/2)   (15.28)
    dp0 = lin(sp.log(p2) + g / (g - 1) * sp.log(Tfac2))                      # p0 = p (T0/T)^{gamma/(gamma − 1)}  (15.29)
    assert z0(dS - (Qh + (g - 1) / g * F))                                   # no alpha, no M: T ds = heat added + friction heat
    assert z0(dT0 - Qh / Tf)                                                 # only heat changes T0
    assert z0(dp0 - (-F - g * M ** 2 / 2 * Qh / Tf))                         # friction costs p0 one for one; heating costs gamma M^2/2 per dT0/T0
    assert z0(dp0 - g / (g - 1) * (dT0 - dS))                                # Gibbs: ds/Cp = dT0/T0 − [(gamma − 1)/gamma] dp0/p0
    for e in (dS, dT0, dp0):
        assert sp.simplify(sp.denom(sp.together(e)).subs(M, 1)) != 0 and sp.limit(e, M, 1).is_finite                 # regular at the sonic point
    # the code: the three closed forms on a field of Mach numbers that contains M = 1 exactly, for zero and non-zero increments
    Mv = np.array([0.0, 0.3, 0.9, 1.0 - 1e-9, 1.0, 1.0 + 1e-9, 1.4, 3.0])
    for gv in GAMMAS:
        tf = 1 + 0.5 * (gv - 1) * Mv ** 2
        for av, fv, qv in ((0.0, 0.0, 0.0), (0.07, 0.0, 0.0), (0.0, 0.021, 0.0), (0.0, 0.0, -0.017), (0.1, 0.1, 0.0), (-0.02, 0.01, 0.03), (0.2, 0.03, 0.05)):
            cc = GD.influence_coefficients(Mv, gv, dA_over_A=av, df=fv, dq=qv)
            assert close(cc["ds_over_cp"], np.full(Mv.shape, qv + (gv - 1) / gv * fv), 1e-14, 1e-17), (gv, av, fv, qv, cc["ds_over_cp"])
            assert close(cc["dT0_over_T0"], qv / tf, 1e-14, 1e-17) and close(cc["dp0_over_p0"], -fv - 0.5 * gv * Mv ** 2 * qv / tf, 1e-14, 1e-17)
            one = GD.influence_coefficients(1.0, gv, dA_over_A=av, df=fv, dq=qv)                                     # the scalar call at the sonic point
            assert all(isinstance(one[k], float) and math.isfinite(one[k]) for k in REGULAR_KEYS)
            assert close([one[k] for k in REGULAR_KEYS], [qv + (gv - 1) / gv * fv, 2 * qv / (gv + 1), -fv - gv * qv / (gv + 1)], 1e-14, 1e-17)


def test_influence_coefficients_V7_sonic_point_signs_zero_drivers_and_printed_forms():  # V7 — the documented contract at M = 1, key by key (fix of F1)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error")                                       # the function owns its 1/0 and 0/0: no numpy warning may leak
        for gv in (1.25, 1.4):
            for kw in (dict(dA_over_A=0.1), dict(dA_over_A=-0.1), dict(df=0.1), dict(dq=0.1), dict(dq=-0.1), dict(dA_over_A=0.1, df=0.03),
                       dict(dA_over_A=0.2, df=0.03, dq=0.05), dict(df=0.02, dq=-0.3)):
                at, below = GD.influence_coefficients(1.0, gv, **kw), GD.influence_coefficients(1.0 - 1e-6, gv, **kw)
                for k in SINGULAR_KEYS:                                      # ±inf with the sign of the subsonic side
                    assert math.isinf(at[k]) and np.sign(at[k]) == np.sign(below[k]) and abs(below[k]) > 1e3, (gv, kw, k, at[k], below[k])
                assert all(math.isfinite(at[k]) for k in REGULAR_KEYS) and close([at[k] for k in REGULAR_KEYS], [below[k] for k in REGULAR_KEYS], 1e-5, 1e-17)    # continuous through M = 1 (they change by O(1e-6) over 1e-6 in M)
                assert at["du_over_u"] == GD.area_velocity(1.0, -(-kw.get("dA_over_A", 0.0) + kw.get("df", 0.0) + kw.get("dq", 0.0)))    # (15.34) with the whole driver
            for kw in (dict(), dict(dA_over_A=0.1, df=0.1), dict(df=0.25, dq=-0.25), dict(dA_over_A=0.5, df=0.25, dq=0.25)):     # the driver is exactly zero: 0/0
                at = GD.influence_coefficients(1.0, gv, **kw)
                assert all(math.isnan(at[k]) for k in SINGULAR_KEYS) and all(math.isfinite(at[k]) for k in REGULAR_KEYS), (kw, at)
        arr = GD.influence_coefficients(np.array([0.5, 1.0, 1.0, 1.0, 2.0]), G, dA_over_A=np.array([0.1, 0.1, -0.1, 0.0, 0.1]))    # an array that passes through M = 1
        assert close(arr["du_over_u"][[0, 4]], [-0.1 / 0.75, 0.1 / 3.0], 1e-14) and arr["du_over_u"][1] == -math.inf and arr["du_over_u"][2] == math.inf and math.isnan(arr["du_over_u"][3])
        assert all(arr[k].shape == (5,) for k in arr) and all(np.all(arr[k] == 0.0) for k in REGULAR_KEYS)
        # the printed (wrong) forms keep their own behaviour at M = 1
        p32 = GD.influence_coefficients(1.0, G, dA_over_A=0.1, printed_15_32=True)       # printed (15.32): the area part is regular, −(dA/A)/(1 + M^2)
        assert close([p32["du_over_u"], p32["drho_over_rho"], p32["dp_over_p"], p32["ds_over_cp"]], [-0.05, -0.05, 0.07, 0.0], 1e-14, 1e-17) and all(math.isfinite(v) for v in p32.values())
        assert all(v == 0.0 for v in GD.influence_coefficients(1.0, G, printed_15_32=True).values())
        assert GD.influence_coefficients(1.0, G, dA_over_A=0.1, df=0.1, printed_15_32=True)["du_over_u"] == math.inf           # its driver is df + dq alone
        p47 = GD.influence_coefficients(1.0, G, df=0.1, printed_15_47=True)              # printed (15.47): a pole in the entropy and in p0 that the correct form does not have
        assert p47["du_over_u"] == math.inf and math.isinf(p47["ds_over_cp"]) and math.isinf(p47["dp0_over_p0"]) and p47["dT0_over_T0"] == 0.0
        p48, ok = GD.influence_coefficients(1.0, G, dq=0.1, printed_15_48=True), GD.influence_coefficients(1.0, G, dq=0.1)
        assert p48["du_over_u"] == -math.inf and ok["du_over_u"] == math.inf and math.isinf(p48["dT0_over_T0"]) and math.isfinite(ok["dT0_over_T0"])    # printed (15.48): heating would decelerate
        for opt in ("printed_15_47", "printed_15_48"):                                   # with no friction / no heat a printed option changes nothing
            assert GD.influence_coefficients(1.0, G, dA_over_A=0.1, **{opt: True}) == GD.influence_coefficients(1.0, G, dA_over_A=0.1)


# =====================================================================================================================
# 4. §15.6 — C07 the normal shock (15.35)–(15.42), C08 the entropy jump (15.43)–(15.44); moving shock; shock structure
# =====================================================================================================================
def test_normal_shock_V4_the_three_jump_laws_hold_for_the_returned_state():  # V4 — (15.35)–(15.37) with fluxes written here; C07
    for g in GAMMAS:
        for M1 in np.concatenate([M_SUP, [1.0, 1.0001, 40.0]]):
            s = GD.normal_shock(M1, g)
            s1, s2 = state(M1, 1.0, 1.0, g), state(s["M2"], s["p_ratio"], s["T_ratio"], g)
            f1, f2 = fluxes(s1, g), fluxes(s2, g)
            assert np.max(np.abs(f2 / f1 - 1)) < 5e-13, (g, M1, f2 / f1 - 1)
            assert close(s["rho_ratio"], s["p_ratio"] / s["T_ratio"], 1e-13) and close(s["u_ratio"] * s["rho_ratio"], 1.0, 1e-14) and close(s2["rho"], s["rho_ratio"], 1e-13)
            res = GD.normal_shock_residuals(s1, s2, g)
            assert max(abs(res[k]) for k in ("mass", "momentum", "energy")) < 5e-13
            fl = GD.normal_shock_fluxes(M1, g)
            assert close([fl["mass"], fl["momentum"], fl["energy"]], 1.0, 1e-12)
            assert close([fl["p1_part"] + fl["ke1_part"], fl["p2_part"] + fl["ke2_part"], fl["h1_part"] + fl["k1_part"], fl["h2_part"] + fl["k2_part"]], 1.0, 1e-12)
            assert close(fl["ke1_part"] / fl["p1_part"], g * M1 ** 2, 1e-12) and close(fl["k2_part"] / fl["h2_part"], 0.5 * (g - 1) * s["M2"] ** 2, 1e-12)
    bad = GD.normal_shock_residuals(state(2.0), state(0.6, 4.5, 1.6875))                 # M2 slightly wrong: the residuals see it
    assert abs(bad["mass"]) > 1e-2
    assert abs(GD.normal_shock_residuals(state(2.0), state(2.0))["energy"]) < 1e-15      # "nothing happens" also satisfies the laws


def test_normal_shock_pressure_jump_V2_derivation():  # V2 — D13 (three stars), every line of Part F: (15.35)–(15.38) → (15.39); C07
    g, a, r, M = sp.symbols("gamma a r M_1", positive=True)                # units rho1 = u1 = 1: a = p1/(rho1 u1^2) = 1/(gamma M1^2); r = u2/u1
    p1 = a
    rho2 = 1 / r                                                          # step 2: mass (15.35)
    p2 = p1 + (1 - r)                                                     # step 3: momentum (15.36): p2 = p1 + rho1 u1^2 (1 − r)
    assert z0((p1 - p2) - (-1 + rho2 * r ** 2))                           # … which is (15.36) itself
    p2_over_rho2 = sp.expand(p2 / rho2)
    assert z0(p2_over_rho2 - r * (a + (1 - r)))                           # step 4
    energy = g / (g - 1) * p1 + sp.Rational(1, 2) - g / (g - 1) * p2_over_rho2 - r ** 2 / 2      # step 5: (15.38), side 1 − side 2
    quad_ = sp.expand(2 * (g - 1) * energy)                                # steps 6–8: clear the fractions: a quadratic in r
    assert sp.degree(quad_, r) == 2 and z0(quad_.subs(r, 1))               # r = 1 ("nothing happens") is always a root
    assert z0(quad_ - (r - 1) * ((g + 1) * r - (2 * g * a + g - 1)))       # step 9: the factorised form of Part F
    roots = sp.solve(quad_, r)
    other = [q for q in roots if sp.simplify(q - 1) != 0]
    assert len(roots) == 2 and len(other) == 1
    r_s = other[0]
    assert z0(r_s - ((g - 1) / (g + 1) + 2 * g * a / (g + 1)))             # step 10
    assert z0(1 + r_s - 2 * g * (a + 1) / (g + 1))                         # Vieta: the roots sum to 2 gamma (a + 1)/(gamma + 1)
    assert z0((1 - r_s) - 2 * (1 - g * a) / (g + 1))                       # step 11
    pr = sp.simplify((1 + (1 - r_s) / a).subs(a, 1 / (g * M ** 2)))        # step 12–13: p2/p1 = 1 + (1 − r)/a with a = 1/(gamma M1^2)
    assert z0(pr - (1 + 2 * g / (g + 1) * (M ** 2 - 1)))                   # (15.39)
    assert z0(r_s.subs(a, 1 / (g * M ** 2)) - ((g - 1) * M ** 2 + 2) / ((g + 1) * M ** 2))      # the root printed by the design's check cell
    assert pr.subs({g: sp.Rational(7, 5), M: 2}) == sp.Rational(9, 2) and z0(pr.subs(M, 1) - 1)  # the worked number and the no-shock limit
    assert sp.limit(pr / M ** 2, M, sp.oo) == 2 * g / (g + 1)
    for gv in GAMMAS:
        assert close(GD.normal_shock_pressure_ratio(M_SUP, gv), 1 + 2 * gv / (gv + 1) * (M_SUP ** 2 - 1), 1e-14)


def test_normal_shock_downstream_mach_V2_derivation():  # V2 — D14 (two stars): p1(1 + gamma M1^2) = p2(1 + gamma M2^2) with (15.39) → (15.40); C07
    g, M1, X = sp.symbols("gamma M_1 X", positive=True)                    # X = M2^2
    pr = 1 + 2 * g / (g + 1) * (M1 ** 2 - 1)
    sol = sp.solve(sp.Eq(1 + g * M1 ** 2, pr * (1 + g * X)), X)
    assert len(sol) == 1 and z0(sol[0] - ((g - 1) * M1 ** 2 + 2) / (2 * g * M1 ** 2 + 1 - g))       # (15.40)
    f = lambda x: ((g - 1) * x + 2) / (2 * g * x + 1 - g)                  # noqa: E731 — the map M1^2 → M2^2
    Y = sp.Symbol("Y", positive=True)
    assert z0(f(f(Y)) - Y)                                                 # its own inverse (trap T6)
    assert sp.limit(f(Y), Y, sp.oo) == (g - 1) / (2 * g) and f(1) == 1
    for gv in GAMMAS:
        M2 = GD.normal_shock_mach(M_SUP, gv)
        assert close(GD.normal_shock_mach(M2, gv), M_SUP, 1e-12) and np.all(M2 < 1)
        assert math.isnan(GD.normal_shock_mach(0.99 * math.sqrt((gv - 1) / (2 * gv)), gv))             # no real root below sqrt((gamma−1)/2gamma)


def test_normal_shock_density_and_temperature_V2_derivation():  # V2 — D15 (two stars): (15.41) from mass + momentum + (15.39); (15.42) = (15.39)/(15.41); C07
    g, M = sp.symbols("gamma M_1", positive=True)
    r = sp.Symbol("r", positive=True)                                      # u2/u1 = rho1/rho2
    pr = 1 + 2 * g / (g + 1) * (M ** 2 - 1)
    sol = sp.solve(sp.Eq(pr, 1 + g * M ** 2 * (1 - r)), r)                 # momentum in units of p1: p2/p1 = 1 + gamma M1^2 (1 − u2/u1)
    assert len(sol) == 1 and z0(1 / sol[0] - (g + 1) * M ** 2 / ((g - 1) * M ** 2 + 2))             # (15.41)
    T = pr * sol[0]
    assert z0(T - (1 + 2 * (g - 1) / (g + 1) ** 2 * (g * M ** 2 + 1) / M ** 2 * (M ** 2 - 1)))      # (15.42), rearranged
    assert sp.limit(1 / sol[0], M, sp.oo) == (g + 1) / (g - 1)
    e = ch15.normal_shock_sympy()                                          # the chapter's engine reaches the same four relations
    assert all(e[k] == 0 for k in ("residual_15_39", "residual_15_40", "residual_15_41", "residual_15_42")) and any(sp.simplify(q - 1) == 0 for q in e["roots"])
    Ms, gs = sorted(e["rho_ratio"].free_symbols, key=str)                  # M_1, gamma
    assert rel(e["rho_ratio"].subs({gs: 1.4, Ms: 2.0}), 8 / 3) < 1e-13 and rel(e["T_ratio"].subs({gs: 1.4, Ms: 2.0}), 1.6875) < 1e-13
    assert rel(e["p_ratio"].subs({gs: 1.4, Ms: 2.0}), 4.5) < 1e-13 and rel(e["M2_squared"].subs({gs: 1.4, Ms: 2.0}), 1 / 3) < 1e-13
    for gv in GAMMAS:
        assert close(GD.normal_shock_temperature_ratio(M_SUP, gv), GD.normal_shock_pressure_ratio(M_SUP, gv) / GD.normal_shock_density_ratio(M_SUP, gv), 1e-13)


def test_normal_shock_V5_naca_1135_and_cited_table_rows():  # V5 — NACA Report 1135 Table II (first-hand) and the secondary table; C07, C08
    b = bench()
    t2 = b["naca1135_table2_supersonic"]
    cols = t2["columns"]
    misses, n = [], 0
    for row in t2["rows"]:
        M1 = row[0]
        s = GD.normal_shock(M1, 1.4)
        ours = {"M2": s["M2"], "p2/p1": s["p_ratio"], "rho2/rho1": s["rho_ratio"], "T2/T1": s["T_ratio"], "pt2/pt1": s["p0_ratio"],
                "p1/pt2": 1 / (GD.stagnation_pressure_ratio(M1) * s["p0_ratio"])}
        for name, val in ours.items():
            k = cols.index(name)
            n += 1
            half = tol_digits(row[k], t2["decimals"][k])
            assert abs(val - row[k]) <= 2 * half * 1.001, (M1, name, val, row[k])
            if abs(val - row[k]) > half * 1.001:
                misses.append((M1, name, val, row[k]))
    assert n == 36 and not misses, misses                                  # all 36 entries agree to half a unit of the last printed digit
    w = b["normal_shock_tables"]
    for row in w["rows"]:
        s = GD.normal_shock(row[0], 1.4)
        assert close([s["M2"], s["p_ratio"], s["rho_ratio"], s["T_ratio"], s["p0_ratio"]], row[1:], 0, 5.1e-5), row
    tab = GD.normal_shock_table([r[0] for r in w["rows"]])
    assert list(tab.columns) == ["M1", "M2", "p2/p1", "rho2/rho1", "T2/T1", "p02/p01"] and close(tab["p02/p01"].to_numpy(), [r[5] for r in w["rows"]], 0, 5.1e-5)
    # the report's forms in xi = p2/p1 (its equations (105), (106), (111), (127)) at our own shock strengths
    for g in GAMMAS:
        xi = GD.normal_shock_pressure_ratio(M_SUP, g)
        rr = ((g + 1) * xi + (g - 1)) / ((g - 1) * xi + (g + 1))
        assert close(GD.normal_shock_density_ratio(M_SUP, g), rr, 1e-13) and close(GD.hugoniot_pressure_ratio(rr, g), xi, 1e-12)        # (106) and its inverse
        assert close(M_SUP ** 2, ((g + 1) * xi + (g - 1)) / (2 * g), 1e-13)                                                           # (105)
        assert close(GD.normal_shock_stagnation_pressure_ratio(M_SUP, g), xi ** (-1 / (g - 1)) * rr ** (g / (g - 1)), 1e-11)            # (111)
        assert close((xi - 1) / (rr - 1), g * (xi + 1) / (rr + 1), 1e-11)                                                             # (127), with p1 = rho1 = 1


def test_normal_shock_V7_limits_guard_and_hugoniot():  # V7 — M1 → 1, M1 → infinity, the upstream guard (trap T6), the Hugoniot beside the isentrope; C07
    for g in GAMMAS:
        s1 = GD.normal_shock(1.0, g)
        assert all(s1[k] == 1.0 for k in ("M2", "p_ratio", "rho_ratio", "u_ratio", "T_ratio", "p0_ratio", "T0_ratio")) and s1["dS_over_cv"] == 0.0 and s1["admissible"] is True
        lim, big = GD.strong_shock_limits(g), GD.normal_shock(1e6, g)      # the corrections are O(1/M1^2)
        assert rel(big["rho_ratio"], lim["rho_ratio"]) < 1e-9 and rel(big["M2"], lim["M2"]) < 1e-9 and rel(big["p_ratio"] / 1e12, lim["p_ratio_over_M1sq"]) < 1e-9
        assert rel(big["T_ratio"] / big["p_ratio"], lim["T_ratio_over_p_ratio"]) < 1e-9
        assert rel(g * 1e12 / big["p_ratio"], lim["u1sq_rho1_over_p2"]) < 1e-9               # rho1 u1^2/p2 = gamma M1^2 p1/p2
        r = np.linspace(1.0, 0.98 * (g + 1) / (g - 1), 50)
        assert np.all(GD.hugoniot_pressure_ratio(r, g) >= r ** g * (1 - 1e-14))              # the Hugoniot lies above the isentrope for compressions
        d = np.array([1e-2, 1e-3])
        gap = GD.hugoniot_pressure_ratio(1 + d, g) - (1 + d) ** g
        assert abs(math.log(gap[0] / gap[1]) / math.log(10) - 3.0) < 0.03                    # … and touches it to second order (gap ~ d^3)
        assert math.isnan(GD.hugoniot_pressure_ratio(1.01 * (g + 1) / (g - 1), g))
    with pytest.raises(ValueError, match="UPSTREAM"):
        GD.normal_shock(0.8)
    with pytest.raises(ValueError):
        GD.normal_shock(np.array([1.5, 0.99]))
    ex = GD.normal_shock(0.8, allow_expansion=True)                         # the forbidden root is the mirror image of an admissible shock
    back = GD.normal_shock(ex["M2"])
    assert ex["admissible"] is False and ex["dS_over_cv"] < 0 and ex["p_ratio"] < 1 and ex["p0_ratio"] > 1 and ex["M2"] > 1
    assert close(back["M2"], 0.8, 1e-13) and close(back["p_ratio"] * ex["p_ratio"], 1.0, 1e-13) and close(back["dS_over_cv"], -ex["dS_over_cv"], 1e-11)
    arr = GD.normal_shock(np.array([1.2, 2.0]))
    assert arr["p_ratio"].shape == (2,) and list(arr["admissible"]) == [True, True]
    assert close(GD.normal_shock_pressure_ratio(0.8), ex["p_ratio"], 1e-15)                 # the bare formula does not refuse (documented)


def test_shock_entropy_V1_equals_the_entropy_state_function_of_the_jump_state():  # V1 — (15.43) against core.thermo's entropy; C08
    import mpmath as mp
    for g in GAMMAS:
        cv, R = 1.0 / (g - 1), 1.0
        for M1 in np.concatenate([M_SUP, [0.7, 0.9, 1.5]]):
            s = GD.normal_shock(M1, g, allow_expansion=True)
            ref = ch15.perfect_gas_entropy_change(1.0, 1.0, s["T_ratio"], 1.0 / s["rho_ratio"], cv=cv, R=R) / cv
            assert abs(GD.normal_shock_entropy(M1, g) - ref) < 1e-12 + 1e-10 * abs(ref), (g, M1)
            assert (GD.normal_shock_entropy(M1, g) > 0) == (M1 > 1) == GD.shock_admissible(M1)       # the second law picks M1 > 1
            assert close(s["p0_ratio"], math.exp(-GD.normal_shock_entropy(M1, g) / (g - 1)), 1e-13)
        mp.mp.dps = 40                                                     # near M1 = 1 the logs cancel to third order: compare with 40-digit arithmetic
        gm = mp.mpf(g)                                                     # all of the reference in 40 digits: the first two orders cancel only if a, b, gamma are consistent
        a, b = 2 * gm / (gm + 1), (gm - 1) / (gm + 1)
        for m in (1e-6, 1e-3, 0.0199, 0.0201, -0.0199, -0.0201, 0.3):
            M1d = math.sqrt(1 + m)
            md = mp.mpf(M1d) ** 2 - 1                                       # the exact m of this double M1
            exact = mp.log1p(a * md) - gm * mp.log1p(md) + gm * mp.log1p(b * md)
            got = GD.normal_shock_entropy(M1d, g)
            # the code forms m = (M1 − 1)(M1 + 1), which keeps its relative precision as M1 → 1: no allowance for the cancellation is needed
            assert abs(got - float(exact)) < 1e-10 * abs(float(exact)), (g, m, got, float(exact))
        Mg = np.linspace(1.0, 8.0, 400)
        assert np.all(np.diff(GD.normal_shock_entropy(Mg, g)) > 0)
    # (15.43) gives S2 − S1 = 0 at M1 = 1 exactly (no jump: a sound wave); the second law asks only S2 − S1 >= 0, so M1 = 1 is
    # admissible — the same rule as normal_shock's "admissible" flag (review Should-fix 5)
    assert list(GD.shock_admissible(np.array([0.5, 1.0, 1.5]))) == [False, True, True]
    Mgrid = np.concatenate([np.linspace(0.5, 3.0, 26), [1.0, np.nextafter(1.0, 0.0), np.nextafter(1.0, 2.0)]])
    assert 1.0 in Mgrid and np.array_equal(GD.shock_admissible(Mgrid), Mgrid >= 1.0)
    for g in GAMMAS:
        for M1 in Mgrid:
            s = GD.normal_shock(M1, g, allow_expansion=True)
            assert s["admissible"] == GD.shock_admissible(M1) == (s["dS_over_cv"] >= 0.0), (g, M1)   # the two rules and the entropy sign agree
        assert GD.normal_shock(1.0, g)["dS_over_cv"] == 0.0 and GD.shock_admissible(1.0) is True
    assert math.isnan(GD.normal_shock_entropy(0.2))                         # p2/p1 <= 0: no state at all


def test_shock_entropy_V1_full_relative_precision_next_to_mach_one():  # V1 + V2 — (15.43) against 50-digit arithmetic for |M1 − 1| down to 1e-12, and against its own series from sympy; C08
    import mpmath as mp
    gs, es = sp.symbols("gamma epsilon", positive=True)
    ser = sp.expand(sp.series(sp.log(1 + 2 * gs / (gs + 1) * es) + gs * sp.log(1 + (gs - 1) / (gs + 1) * es) - gs * sp.log(1 + es), es, 0, 6).removeO())
    assert sp.simplify(ser.coeff(es, 1)) == 0 and sp.simplify(ser.coeff(es, 2)) == 0
    c3, c4, c5 = (sp.lambdify(gs, sp.simplify(ser.coeff(es, k)), "math") for k in (3, 4, 5))
    old = mp.mp.dps
    mp.mp.dps = 50
    try:
        for g in GAMMAS:
            gm = mp.mpf(g)
            a, b = 2 * gm / (gm + 1), (gm - 1) / (gm + 1)
            for d in (1e-12, -1e-12, 1e-9, -1e-9, 1e-7, -1e-7, 1e-5, -1e-5, 1e-3, -1e-3, 9e-3, -9e-3):      # |M1^2 − 1| < 0.02
                M1 = 1.0 + d
                md = mp.mpf(M1) ** 2 - 1                                    # the exact M1^2 − 1 of this double
                exact = mp.log1p(a * md) - gm * mp.log1p(md) + gm * mp.log1p(b * md)
                got = GD.normal_shock_entropy(M1, g)
                assert abs(mp.mpf(got) - exact) < mp.mpf("1e-13") * abs(exact), (g, d, got, float(exact))    # M1**2 − 1 in doubles would lose 2.2e-16/|m|: 1e-9 at d = 1e-7
                assert (got > 0) == (d > 0)
    finally:
        mp.mp.dps = old
    for g in GAMMAS:                                                        # the series sympy derives from (15.43): cubic + quartic, remainder of fifth order
        for m in (1e-3, 3e-3, 1e-2, -1e-3, -1e-2):
            got = GD.normal_shock_entropy(math.sqrt(1 + m), g)
            assert abs(got - (c3(g) * m ** 3 + c4(g) * m ** 4)) < 1.5 * abs(c5(g) * m ** 5) + 1e-13 * abs(got), (g, m)
    up, down = GD.normal_shock_entropy(np.nextafter(1.0, 2.0)), GD.normal_shock_entropy(np.nextafter(1.0, 0.0))
    assert GD.normal_shock_entropy(1.0) == 0.0 and 0.0 < up < 1e-44 and -1e-44 < down < 0.0                 # no jump at M1 = 1; the sign is right one ulp to either side


def test_weak_shock_entropy_V2_derivation():  # V2 — D17 (three stars), every line of Part F: (15.43) → (15.44a) → (15.44b); C08
    g, eps, M, D = sp.symbols("gamma epsilon M_1 Delta", positive=True)
    a, b = 2 * g / (g + 1), (g - 1) / (g + 1)
    sub = {eps: M ** 2 - 1}                                                 # step 1
    assert z0((1 + a * eps).subs(sub) - (1 + 2 * g / (g + 1) * (M ** 2 - 1)))                        # step 2: p2/p1 = 1 + a eps
    assert z0(((1 + b * eps) / (1 + eps)).subs(sub) - ((g - 1) * M ** 2 + 2) / ((g + 1) * M ** 2))   # step 3: rho1/rho2 = (1 + b eps)/(1 + eps)
    # steps 4–8: ln(1 + x) = x − x^2/2 + x^3/3 − x^4/4 applied to the three brackets
    ln = lambda x, n: sum((-1) ** (k + 1) * x ** k / k for k in range(1, n + 1))                     # noqa: E731
    S = sp.expand(ln(a * eps, 4) + g * ln(b * eps, 4) - g * ln(eps, 4))
    c = [sp.simplify(S.coeff(eps, k)) for k in (1, 2, 3, 4)]
    assert c[0] == 0 and c[1] == 0                                          # steps 6–7: the first two orders cancel identically
    assert z0(c[2] - 2 * g * (g - 1) / (3 * (g + 1) ** 2))                  # step 8–9: (15.44a)
    assert z0(c[3] + 2 * g ** 2 * (g - 1) / (g + 1) ** 3)                   # the next term quoted in Part F's check
    direct = sp.series(sp.expand_log(sp.log((1 + a * eps) * ((1 + b * eps) / (1 + eps)) ** g), force=True), eps, 0, 5).removeO()   # sympy expands (15.43) on its own
    assert z0(sp.expand(direct) - (c[2] * eps ** 3 + c[3] * eps ** 4))
    eps_of_D = (g + 1) / (2 * g) * D                                        # step 10: (15.39) inverted
    assert z0((a * eps).subs(eps, eps_of_D) - D)
    assert z0((c[2] * eps ** 3).subs(eps, eps_of_D) - (g ** 2 - 1) / (12 * g ** 2) * D ** 3)         # step 11: (15.44b)
    # the same two coefficients are printed in NACA Report 1135 (its (147), divided by gamma − 1, and its weak-shock series in xi)
    assert z0(c[2] / (g - 1) - 2 * g / (3 * (g + 1) ** 2)) and z0(c[3] / (g - 1) + 2 * g ** 2 / (g + 1) ** 3)
    f = bench()["naca1135_forms"]["eq_146_gamma_7_5"]
    assert (c[2] / (g - 1)).subs(g, sp.Rational(7, 5)) == sp.Rational(*f["cubic"]) and (-c[3] / (g - 1)).subs(g, sp.Rational(7, 5)) == sp.Rational(*f["quartic"])
    e = ch15.weak_shock_series(4)                                           # the chapter's engine
    assert e["coefficients"][0] == 0 and e["coefficients"][1] == 0 and e["residual_15_44a"] == 0 and e["residual_15_44b"] == 0
    assert z0(e["coefficients"][3] - c[3].subs(g, sp.Symbol("gamma", positive=True))) and z0(e["pressure_coefficient"] - (sp.Symbol("gamma", positive=True) ** 2 - 1) / (12 * sp.Symbol("gamma", positive=True) ** 2))
    for gv in GAMMAS:                                                       # the coded forms
        Mv = np.array([1.001, 1.05, 1.1])
        assert close(GD.weak_shock_entropy(Mv, gv), 2 * gv * (gv - 1) / (3 * (gv + 1) ** 2) * (Mv ** 2 - 1) ** 3, 1e-13)
        dp = GD.normal_shock_pressure_ratio(Mv, gv) - 1
        assert close(GD.weak_shock_entropy(dp, gv, form="pressure"), GD.weak_shock_entropy(Mv, gv, form="mach"), 1e-12)
    with pytest.raises(ValueError):
        GD.weak_shock_entropy(1.1, form="density")
    # the numbers of Part F's check: at M1 = 1.1 the cubic gives 6.0e-4 against the exact 4.29e-4 (Part F prints 4.2e-4 — it should read 4.3e-4;
    # reported to the designer); at M1 = 1.05, 7.0e-5 against 5.9e-5
    assert abs(GD.weak_shock_entropy(1.1) - 6.0e-4) < 0.05e-4 and abs(GD.normal_shock_entropy(1.1) - 4.29e-4) < 0.005e-4
    assert abs(GD.weak_shock_entropy(1.05) - 7.0e-5) < 0.05e-5 and abs(GD.normal_shock_entropy(1.05) - 5.9e-5) < 0.05e-5
    two_terms = 2 * 1.4 * 0.4 / (3 * 2.4 ** 2) * 0.21 ** 3 - 2 * 1.4 ** 2 * 0.4 / 2.4 ** 3 * 0.21 ** 4       # cubic + quartic: 3.8e-4, on the other side of the exact value
    assert two_terms < GD.normal_shock_entropy(1.1) < GD.weak_shock_entropy(1.1)


def test_weak_shock_entropy_V1_cubic_law_and_its_first_correction():  # V1 — (15.43) against (15.44a) as M1 → 1: slope 3, ratio 1 − 3 gamma m/(gamma + 1); C08
    for g in GAMMAS:
        m = np.array([1e-2, 3e-3, 1e-3, 3e-4, 1e-4])
        M1 = np.sqrt(1 + m)
        ex = GD.normal_shock_entropy(M1, g)
        assert abs(observed_order(m, ex) - 3.0) < 0.05
        ratio = ex / GD.weak_shock_entropy(M1, g)
        assert close(ratio - 1, -3 * g / (g + 1) * m, 0.03, 1e-9)           # the quartic term of the series
        assert np.all(GD.normal_shock_entropy(np.sqrt(1 - m), g) < 0) and close(GD.normal_shock_entropy(np.sqrt(1 - m[2:]), g), -ex[2:], 5e-3)     # odd in M1^2 − 1
        xi = GD.normal_shock_pressure_ratio(M1, g)                          # NACA's weak-shock series in xi: p02/p01 = 1 − (gamma+1)/(12 gamma^2)(xi−1)^3 + (gamma+1)/(8 gamma^2)(xi−1)^4
        series = 1 - (g + 1) / (12 * g ** 2) * (xi - 1) ** 3 + (g + 1) / (8 * g ** 2) * (xi - 1) ** 4
        assert np.all(np.abs(GD.normal_shock_stagnation_pressure_ratio(M1, g) - series) < 2.0 * (xi - 1) ** 5 + 3e-16)


def test_shock_stagnation_pressure_ratio_V2_derivation():  # V2 — D18 (two stars): three routes to p02/p01 (chain of ratios, exp(−dS/R), A1*/A2*); C08
    g, M = sp.symbols("gamma M_1", positive=True)
    M2sq = ((g - 1) * M ** 2 + 2) / (2 * g * M ** 2 + 1 - g)
    pr = 1 + 2 * g / (g + 1) * (M ** 2 - 1)
    rr = (g + 1) * M ** 2 / ((g - 1) * M ** 2 + 2)
    p0 = lambda x: (1 + (g - 1) / 2 * x) ** (g / (g - 1))                   # noqa: E731 — (15.29) as a function of M^2
    chain = pr * p0(M2sq) / p0(M ** 2)                                      # p02/p01 = (p2/p1)(p0/p)_2/(p0/p)_1
    closed = rr ** (g / (g - 1)) * ((g + 1) / (2 * g * M ** 2 - (g - 1))) ** (1 / (g - 1))      # the form "ours" (NACA's (142) with sin(theta) = 1)
    entropy = sp.exp(-sp.log(pr * rr ** (-g)) / (g - 1))                    # exp(−(S2 − S1)/R), R = (gamma − 1) Cv, with (15.43)
    assert z0(sp.expand_log(sp.log(closed) - sp.log(entropy), force=True))  # two of the three routes agree identically
    assert z0(1 + (g - 1) / 2 * M2sq - (1 + (g - 1) / 2 * M ** 2) / (pr / rr))      # T0 unchanged: (T0/T)_2 = (T0/T)_1 / (T2/T1) — the step that turns the chain into the closed form
    for gv in (sp.Rational(7, 5), sp.Rational(5, 3), sp.Rational(9, 8)):
        for mv in (sp.Rational(11, 10), 2, 7):
            vals = [sp.N(e.subs({g: gv, M: mv}), 30) for e in (chain, closed, entropy)]
            assert abs(vals[0] / vals[1] - 1) < sp.Float("1e-25") and abs(vals[2] / vals[1] - 1) < sp.Float("1e-25")
    for gv in GAMMAS:                                                       # third route: the sonic area grows by the same factor (mass conservation at fixed T0)
        M2 = GD.normal_shock_mach(M_SUP, gv)
        k = GD.normal_shock_stagnation_pressure_ratio(M_SUP, gv)
        assert close(k, GD.area_ratio(M2, gv) / GD.area_ratio(M_SUP, gv), 1e-11)                 # A1*/A2*
        assert close(k, GD.mass_flux_function(M_SUP, gv) / GD.mass_flux_function(M2, gv), 1e-11)
        assert close(k, GD.normal_shock_pressure_ratio(M_SUP, gv) * GD.stagnation_pressure_ratio(M2, gv) / GD.stagnation_pressure_ratio(M_SUP, gv), 1e-11)
        assert np.all(k < 1) and GD.normal_shock_stagnation_pressure_ratio(0.8, gv) > 1


def test_moving_shock_V4_conservation_in_the_laboratory_frame():  # V4 — N62–N64: fluxes through a surface moving with the front; trap T7
    for gas in (GD.PerfectGas(1.4, 287.0), GD.PerfectGas(5.0 / 3.0, 2077.0)):
        g, cv = gas.gamma, gas.cv
        for Ms in (1.0005, I["M_s"], 2.7, 9.0):
            m = GD.moving_shock(0.93e5, 281.0, M_s=Ms, gas=gas)
            W, ug, r1, r2, p1, p2 = m["u_shock"], m["u_gas"], m["rho1"], m["rho2"], 0.93e5, m["p2"]
            assert rel(r2 * (W - ug), r1 * W) < 1e-12                       # mass through the moving front
            assert rel(p2 - p1, r1 * W * ug) < 1e-9 * max(1.0, p1 / (p2 - p1))     # momentum given to the gas set in motion
            assert rel(p2 * ug, r1 * W * (cv * (m["T2"] - 281.0) + 0.5 * ug ** 2)) < 1e-9 * max(1.0, p1 / (p2 - p1))    # work of the piston = energy gained
            assert close(m["c1"], gas.c(281.0), 1e-14) and close(W, Ms * m["c1"], 1e-14) and close(m["M_gas"], ug / m["c2"], 1e-14) and W >= m["c1"]
            assert close(m["T0_behind_lab"], m["T2"] + ug ** 2 / (2 * gas.cp), 1e-12) and m["T0_behind_lab"] > 281.0     # T0 depends on the frame
            assert m["p0_behind_lab"] >= p2 and close(r1, p1 / (gas.R * 281.0), 1e-14)
            again = GD.moving_shock(0.93e5, 281.0, p2_over_p1=p2 / p1, gas=gas)                 # NACA (105) inverted
            assert rel(again["M_s"], Ms) < 1e-9 and rel(again["u_gas"], ug) < 1e-6
        weak = GD.moving_shock(0.93e5, 281.0, p2_over_p1=1 + 1e-6, gas=gas)   # weak limit: a sound wave, u' = p'/(rho0 c)  (15.16)
        assert rel(weak["u_shock"], weak["c1"]) < 1e-6 and rel(weak["u_gas"], 1e-6 * 0.93e5 / AC.acoustic_impedance(weak["rho1"], weak["c1"])) < 1e-5
    for kw in ({}, dict(M_s=1.5, p2_over_p1=2.0)):
        with pytest.raises(ValueError):
            GD.moving_shock(1e5, 300.0, **kw)
    with pytest.raises(ValueError):
        GD.moving_shock(1e5, 300.0, M_s=0.9)


def test_blast_constant_V1_crude_model_integrates_to_the_closed_form():  # V1 (+ cited comparison, qualitative) — R17, N64; trap T22
    for g in (1.4, 5.0 / 3.0):
        for geom, vol in (("hemisphere", 2 * PI / 3), ("sphere", 4 * PI / 3)):
            K = ch15.blast_constant_crude(g, geom)
            C2 = 0.5 * (g + 1) * (g - 1) / vol                              # (dR/dt)^2 R^3 rho1/E from u1^2 = (gamma+1) p2/(2 rho1), p2 = (gamma−1) E/Volume
            E, rho, t0, t1 = 3.0e9, 1.1, 1e-3, 0.05
            sol = solve_ivp(lambda t, R: math.sqrt(C2 * E / rho) * R ** -1.5, (t0, t1), [K * (E * t0 ** 2 / rho) ** 0.2], rtol=1e-11, atol=1e-12)
            assert rel(sol.y[0, -1], K * (E * t1 ** 2 / rho) ** 0.2) < 1e-8
            assert rel(ch15.blast_radius(E, t1, rho, K ** -5), K * (E * t1 ** 2 / rho) ** 0.2) < 1e-13    # ch01 writes D = (E t^2/(K rho))^(1/5): its K is ours to the power −5 (T22)
            assert close(K, (75 * (g * g - 1) / (16 * PI) * (1.0 if geom == "hemisphere" else 0.5)) ** 0.2, 1e-13)
        assert close(ch15.blast_constant_crude(g, "hemisphere") / ch15.blast_constant_crude(g, "sphere"), 2 ** 0.2, 1e-13)
    with pytest.raises(ValueError):
        ch15.blast_constant_crude(1.4, "cylinder")
    S = json.loads((ROOT / "reference" / "ch01" / "benchmarks.json").read_text(encoding="utf-8"))["taylor_blast"]["S_gamma"]    # cited similarity constant (sphere, gamma = 1.4)
    low = 1 - ch15.blast_constant_crude(1.4, "sphere") / S
    assert 0.08 < low < 0.11                                                # the crude constant is 9.4 % low: "qualitative" (Open item) — not "a few per cent"
    lim = GD.strong_shock_limits(1.4)
    assert close(lim["u1sq_rho1_over_p2"], 1.2, 1e-14) and close(lim["rho_ratio"], 6.0, 1e-13)


def test_shock_structure_V2_first_integrals_and_roots():  # V2 + V4 — R19, N65–N69: the unnumbered chain of §15.6
    e = ch15.shock_structure_sympy()
    assert all(e[k] == 0 for k in ("momentum_integral", "energy_integral", "becker_form", "structure_equation"))
    # independent: with H = I and p = rho R T the momentum integral, scaled with V, is  U' = −[(gamma+1)/(2 gamma)] (U1 − U)(U − U2)/U
    g, U, U1, U2, th = sp.symbols("gamma U U_1 U_2 theta", positive=True)       # th = 2 I (gamma−1)/((gamma+1) V^2) = U1 U2 (Vieta)
    quad_ = -U ** 2 + 2 * g / (g + 1) * U - th                                  # right side of the printed structure equation over V^2
    r = sp.solve(quad_, U)
    assert z0(r[0] + r[1] - 2 * g / (g + 1)) and z0(r[0] * r[1] - th)
    assert z0(quad_.subs(th, U1 * U2).subs(U1, 2 * g / (g + 1) - U2) - (2 * g / (g + 1) - U2 - U) * (U - U2))   # ≡ (U1 − U)(U − U2): slip #7 is only the missing V^2
    for gas in (AIR, GD.PerfectGas(5.0 / 3.0, 208.0)):
        gm = gas.gamma
        for M1 in (1.2, I["M1_structure"], 3.5):
            c = ch15.shock_structure_constants(M1, 281.0, 0.93e5, gas)
            s = GD.normal_shock(M1, gm)
            T2, p2, u2 = 281.0 * s["T_ratio"], 0.93e5 * s["p_ratio"], c["u2"]
            assert rel(c["m"] * u2 + p2, c["m"] * c["V"]) < 1e-12 and rel(gas.cp * T2 + 0.5 * u2 ** 2, c["I"]) < 1e-12     # the same m, V, I downstream (V4)
            assert rel(c["m"], p2 / (gas.R * T2) * u2) < 1e-12
            disc = math.sqrt(1 - 2 * (gm ** 2 - 1) * c["I"] / (gm ** 2 * c["V"] ** 2))       # the book's root formula
            U1b, U2b = gm / (gm + 1) * (1 + disc), gm / (gm + 1) * (1 - disc)
            R1, R2 = ch15.becker_shock_roots(M1, gm)
            assert close([R1, R2], [U1b, U2b], 1e-11) and close([c["U1"], c["U2"]], [R1, R2], 1e-12) and close(R1 + R2, 2 * gm / (gm + 1), 1e-14)
            assert close(R1 / R2, s["rho_ratio"], 1e-12)                        # U1/U2 = u1/u2 = rho2/rho1
    assert close(ch15.becker_shock_roots(2.0, 1.4), [28 / 33, 7 / 22], 1e-14)   # slip #8: these plateaus belong to M1 = 2
    R1, R2 = ch15.becker_shock_roots(np.array([1.0, 2.0]))
    assert R1.shape == (2,) and R1[0] == R2[0]


def test_becker_shock_profile_V1_satisfies_the_implicit_relation_and_the_ode():  # V1 + V3 + V4 + V7 — N66–N69, Fig. 15.11 (Becker, Pr'' = 1)
    for g, M1 in ((1.4, 1.3), (1.4, 2.0), (5.0 / 3.0, 3.0)):
        U1, U2 = ch15.becker_shock_roots(M1, g)
        kk = (g + 1) / (2 * g)
        eta = np.linspace(-9.0, 7.0, 161)
        U = ch15.becker_shock_profile(eta, M1, g)
        lhs = lambda u: (U1 * np.log(U1 - u) - U2 * np.log(u - U2)) / (U1 - U2)          # noqa: E731 — the printed implicit relation
        mid = 0.5 * (U1 + U2)
        ok = (U1 - U > 1e-6) & (U - U2 > 1e-6)
        assert np.max(np.abs(lhs(U[ok]) - lhs(mid) - kk * eta[ok])) < 1e-8 and ok.sum() > 60
        assert close(ch15.becker_shock_profile(0.0, M1, g), mid, 1e-12) and np.all(np.diff(U) < 0) and np.all((U < U1) & (U > U2))
        h = 1e-4                                                              # the ODE: dU/d eta by central differences
        dU = (ch15.becker_shock_profile(eta[ok] + h, M1, g) - ch15.becker_shock_profile(eta[ok] - h, M1, g)) / (2 * h)
        assert np.max(np.abs(dU - ch15.becker_shock_ode(eta[ok], U[ok], M1, g))) < 1e-7
        sol = solve_ivp(lambda e_, u_: ch15.becker_shock_ode(e_, u_, M1, g), (-3.0, 6.0), [float(ch15.becker_shock_profile(-3.0, M1, g))], rtol=1e-11, atol=1e-13, dense_output=True)
        ee = np.linspace(-3.0, 6.0, 37)
        assert np.max(np.abs(sol.sol(ee)[0] - ch15.becker_shock_profile(ee, M1, g))) < 1e-8      # V3: the independent numerical route
        assert abs(ch15.becker_shock_profile(-80.0, M1, g) - U1) < 1e-9 and abs(ch15.becker_shock_profile(80.0, M1, g) - U2) < 1e-9      # plateaus
        assert close(ch15.becker_shock_profile(eta + 1.5, M1, g, eta0=1.5), U, 1e-10)             # the origin of eta is arbitrary (T20)
        # V4: momentum flux with the viscous stress and total enthalpy are constant through the profile (dimensional constants of our own state)
        gas = GD.PerfectGas(g, 287.0)
        c = ch15.shock_structure_constants(M1, 281.0, 0.93e5, gas)
        u = U * c["V"]
        T = (c["I"] - 0.5 * u ** 2) / gas.cp                                  # Becker: H = I everywhere
        visc = c["m"] * c["V"] * ch15.becker_shock_ode(eta, U, M1, g)         # mu'' du/dx = m V dU/d eta
        assert np.max(np.abs(c["m"] * u + c["m"] / u * gas.R * T - visc - c["m"] * c["V"])) < 1e-10 * c["m"] * c["V"]
        fine = np.linspace(-12, 12, 48001)                                    # thickness against the numerical maximum slope
        slope = np.max(np.abs(ch15.becker_shock_ode(fine, ch15.becker_shock_profile(fine, M1, g), M1, g)))
        assert rel(ch15.becker_shock_thickness(M1, g), (U1 - U2) / slope) < 1e-6
    th = ch15.becker_shock_thickness(np.array([1.001, 1.01, 1.1, 2.0, 5.0]))
    assert np.all(np.diff(th) < 0) and th[0] > 1000 and math.isinf(ch15.becker_shock_thickness(1.0))     # weak shocks are thick
    assert rel(th[0] * (1.001 ** 2 - 1), 4 * 1.4) < 5e-3                      # V7: thickness → 4 gamma/(M1^2 − 1) as M1 → 1 (from U1 − U2 ≈ 2 gamma m/(gamma+1)^2)
    nu, du = 1.5e-5, 300.0                                                    # order of magnitude: delta ~ nu/(u1 − u2) is of the order of the mean free path
    d = ch15.shock_thickness_estimate(du, nu)
    assert close(d, nu / du, 1e-15) and 0.1 < d / ch15.mean_free_path_air(288.15) < 10
    assert_dimensions(Q_(nu, "m**2/s") / Q_(du, "m/s"), "[length]", "thickness")


# =====================================================================================================================
# 5. §15.5, §15.7 — C09 the convergent nozzle, C10 the Laval nozzle; the worked nozzle chains with our inputs
# =====================================================================================================================
def test_convergent_nozzle_V1_piecewise_mass_flow_against_back_pressure():  # V1 + V4 — C09, D19 (ours): rises, then the choked plateau; trap T10
    for g in GAMMAS:
        crit = GD.critical_pressure_ratio(g)
        r = np.linspace(0.0, 1.0, 2001)
        s = GD.mass_flow_vs_back_pressure(r, g)
        assert np.all(np.diff(s["mdot"]) <= 1e-15) and s["mdot"][-1] == 0.0 and close(s["mdot"][r <= crit], GD.mass_flux_function(1.0, g), 1e-14)
        assert np.array_equal(s["choked"], r <= crit) and close(s["p_exit_over_p0"], np.maximum(r, crit), 1e-15) and np.all(s["M_exit"] <= 1.0)
        assert close(s["M_exit"][r > crit], GD.mach_from_pressure_ratio(r[r > crit], g), 1e-14) and close(s["mdot_over_choked"], s["mdot"] / GD.mass_flux_function(1.0, g), 1e-14)
        h = np.array([1e-3, 1e-4, 1e-5])
        drop = 1 - GD.mass_flow_vs_back_pressure(crit + h, g)["mdot_over_choked"]
        assert np.all(drop > 0) and abs(observed_order(h, drop) - 2.0) < 0.05       # value AND slope continuous at the critical ratio: the curve leaves the plateau quadratically
    gas = GD.PerfectGas(1.4, 287.0)
    p0, T0, Ae = I["p0"], I["T0"], I["A_conv_exit"]
    for pB in (0.95 * p0, 0.7 * p0, 0.53 * p0, 0.5 * p0, 0.1 * p0):
        n = GD.convergent_nozzle(p0, T0, pB, Ae, gas)
        rho_e = n["p_exit"] / (gas.R * n["T_exit"])
        assert rel(n["mdot"], GD.mass_flow(rho_e, n["u_exit"], Ae)) < 1e-12 and close(n["u_exit"], n["M_exit"] * gas.c(n["T_exit"]), 1e-14)     # mdot = rho u A at the exit
        assert close(n["mdot_max"], GD.choked_mass_flow(p0, T0, Ae, gas), 1e-14) and n["choked"] == (pB / p0 <= GD.critical_pressure_ratio(1.4))
        assert n["matched"] == (not n["choked"]) and (close(n["p_exit"], pB, 1e-14) if not n["choked"] else n["p_exit"] > pB)                  # exit pressure = back pressure only when subsonic
        assert rel(gas.cp * n["T_exit"] + 0.5 * n["u_exit"] ** 2, gas.cp * T0) < 1e-13
    assert GD.convergent_nozzle(p0, T0, p0, Ae, gas)["mdot"] == 0.0
    assert close(GD.convergent_nozzle(2 * p0, T0, 0.2 * p0, Ae, gas)["mdot"], 2 * GD.convergent_nozzle(p0, T0, 0.1 * p0, Ae, gas)["mdot"], 1e-13)   # choked flow still rises with p0
    dim = GD.mass_flow_vs_back_pressure(0.8, 1.4, dimensionless=False, p0=p0, T0=T0, A_exit=Ae, gas=gas)
    assert close(dim["mdot"], GD.mass_flow_vs_back_pressure(0.8)["mdot"] * p0 * Ae / math.sqrt(287.0 * T0), 1e-14)
    for bad in (1.01, -0.01):
        with pytest.raises(ValueError):
            GD.mass_flow_vs_back_pressure(bad)
    with pytest.raises(ValueError):
        GD.mass_flow_vs_back_pressure(0.8, dimensionless=False)


def test_convergent_nozzle_mass_flow_V1_one_gamma_per_call():  # V1 — C09, review Must-fix 1: gamma and R from the same gas (the old code mixed gamma = 1.4 with gas.R)
    g = GD.PerfectGas(1.25, 320.0)
    p0, T0, Ae = 3e6, 3000.0, 0.01
    gm, R = g.gamma, g.R

    def by_hand(pe_over_p0):                                               # independent single-gamma evaluation of rho u A at the exit
        Me = math.sqrt(2 / (gm - 1) * ((1 / pe_over_p0) ** ((gm - 1) / gm) - 1))                     # (15.29) inverted
        return p0 * Ae / math.sqrt(R * T0) * math.sqrt(gm) * Me * (1 + 0.5 * (gm - 1) * Me ** 2) ** (-(gm + 1) / (2 * (gm - 1)))

    star = p0 * Ae * math.sqrt(gm / (R * T0)) * (2 / (gm + 1)) ** ((gm + 1) / (2 * (gm - 1)))           # choked flow, closed form
    s = GD.mass_flow_vs_back_pressure(0.3, dimensionless=False, p0=p0, T0=T0, A_exit=Ae, gas=g)
    assert s["choked"] is True and rel(s["mdot"], GD.choked_mass_flow(p0, T0, Ae, g)) < 1e-14 and rel(s["mdot"], star) < 1e-13
    assert abs(s["mdot"] - 20.149037) < 1e-5 and abs(s["mdot"] - 20.966) > 0.8                       # not the old mixed-gamma value
    crit = (2 / (gm + 1)) ** (gm / (gm - 1))
    assert rel(s["p_exit_over_p0"], crit) < 1e-14 and s["M_exit"] == 1.0 and rel(by_hand(crit), star) < 1e-13
    r = np.linspace(crit + 1e-3, 0.999, 60)                                 # off the plateau: subsonic exit, p_exit = p_B
    sub = GD.mass_flow_vs_back_pressure(r, dimensionless=False, p0=p0, T0=T0, A_exit=Ae, gas=g)
    assert not np.any(sub["choked"]) and close(sub["mdot"], [by_hand(x) for x in r], 1e-12)
    assert close(sub["mdot_over_choked"], sub["mdot"] / star, 1e-12) and np.all(np.diff(sub["mdot"]) < 0)
    for x in (0.3, 0.8):                                                    # an explicit gamma that agrees gives the same answer
        a = GD.mass_flow_vs_back_pressure(x, dimensionless=False, p0=p0, T0=T0, A_exit=Ae, gas=g)
        b = GD.mass_flow_vs_back_pressure(x, 1.25, dimensionless=False, p0=p0, T0=T0, A_exit=Ae, gas=g)
        c = GD.mass_flow_vs_back_pressure(x, gamma=1.25, dimensionless=False, p0=p0, T0=T0, A_exit=Ae, gas=g)
        assert a["mdot"] == b["mdot"] == c["mdot"] and a["M_exit"] == b["M_exit"] == c["M_exit"]
    nd = GD.mass_flow_vs_back_pressure(0.8, gas=g)                          # dimensionless with a gas: gamma = gas.gamma too
    assert close(nd["mdot"], GD.mass_flow_vs_back_pressure(0.8, 1.25)["mdot"], 1e-15)
    with pytest.raises(ValueError, match="gamma"):
        GD.mass_flow_vs_back_pressure(0.3, 1.4, dimensionless=False, p0=p0, T0=T0, A_exit=Ae, gas=g)
    with pytest.raises(ValueError, match="gamma"):
        GD.mass_flow_vs_back_pressure(0.3, gamma=1.4, gas=g)
    # air / default calls unchanged: gamma=None means 1.4; the parity value of the explainers still holds
    rr = np.linspace(0.0, 1.0, 101)
    for k in ("mdot", "M_exit", "p_exit_over_p0", "mdot_over_choked"):
        assert np.array_equal(GD.mass_flow_vs_back_pressure(rr)[k], GD.mass_flow_vs_back_pressure(rr, 1.4)[k])
        assert np.array_equal(GD.mass_flow_vs_back_pressure(rr)[k], GD.mass_flow_vs_back_pressure(rr, gas=AIR)[k])
    assert abs(GD.mass_flow_vs_back_pressure(0.8, 1.4)["mdot"] - 0.56066) < 5e-6
    gm, R = 1.4, 287.0                                                      # dimensional air call against the same hand formula
    da = GD.mass_flow_vs_back_pressure(0.8, dimensionless=False, p0=p0, T0=T0, A_exit=Ae, gas=AIR)
    assert rel(da["mdot"], by_hand(0.8)) < 1e-12
    assert rel(GD.convergent_nozzle(p0, T0, 0.8 * p0, Ae, g)["mdot"], GD.mass_flow_vs_back_pressure(0.8, dimensionless=False, p0=p0, T0=T0, A_exit=Ae, gas=g)["mdot"]) < 1e-15


def test_convergent_nozzle_mass_flow_V2_derivation():  # V2 — D19 (two stars, ours): rho u A with (15.28)–(15.30) → the flux function; zero slope at the plateau; C09
    g, M, R, T0, p0, A = sp.symbols("gamma M R T_0 p_0 A", positive=True)
    Tf = 1 + (g - 1) / 2 * M ** 2
    T, p = T0 / Tf, p0 * Tf ** (-g / (g - 1))
    mdot = p / (R * T) * M * sp.sqrt(g * R * T) * A                          # rho u A, rho = p/(R T), u = M sqrt(gamma R T)
    target = sp.sqrt(g) * M * Tf ** (-(g + 1) / (2 * (g - 1)))
    assert sp.simplify(-g / (g - 1) + 1 - sp.Rational(1, 2) + (g + 1) / (2 * (g - 1))) == 0     # the exponents: −gamma/(gamma−1) + 1 − 1/2 = −(gamma+1)/(2(gamma−1))
    for gv in (sp.Rational(7, 5), sp.Rational(5, 3), sp.Rational(5, 4)):
        assert sp.simplify((mdot * sp.sqrt(R * T0) / (p0 * A) / target).subs(g, gv)) == 1
    dflux = sp.simplify(sp.diff(sp.log(target), M))
    dp = sp.simplify(sp.diff(sp.log(p), M))
    assert z0(dflux - (1 - M ** 2) / (M * Tf)) and z0(dp + g * M / Tf)        # d(flux)/dM = 0 at M = 1 while dp/dM ≠ 0 there: d(mdot)/d(p_B) = 0 at the critical ratio
    assert sp.simplify((p / p0).subs(M, 1) - (2 / (g + 1)) ** (g / (g - 1))) == 0


def test_laval_nozzle_V1_three_thresholds_and_seven_regimes():  # V1 — C10, Fig. 15.13 (formulas ours): p_f < p_d < p_b < p0 and the names on both sides
    for g in (1.25, 1.4, 5.0 / 3.0):
        for ar in (1.02, 1.5, I["Ae_over_At"], 5.0, 25.0):
            c = GD.laval_critical_pressures(ar, g)
            assert c["p_f"] < c["p_d"] < c["p_b"] < 1.0 and c["M_sub"] < 1 < c["M_sup"] and c["p_star"] == GD.critical_pressure_ratio(g)
            assert close(c["p_b"], 1 / GD.stagnation_pressure_ratio(GD.mach_from_area_ratio(ar, g, branch="subsonic"), g), 1e-14)
            assert close(c["p_d"], c["p_f"] * GD.normal_shock_pressure_ratio(c["M_sup"], g), 1e-14) and close(GD.area_ratio(c["M_sup"], g), ar, 1e-11)
            e = 1e-6
            names = [GD.laval_regime(x, ar, g) for x in (1.0, 0.5 * (1 + c["p_b"]), c["p_b"] * (1 + e), c["p_b"], c["p_b"] * (1 - e), c["p_d"] * (1 + e), c["p_d"],
                                                         c["p_d"] * (1 - e), c["p_f"] * (1 + e), c["p_f"], c["p_f"] * (1 - e), 0.0)]
            assert names == ["no flow", "subsonic", "subsonic", "sonic throat", "shock in nozzle", "shock in nozzle", "shock at exit", "overexpanded",
                             "overexpanded", "design", "underexpanded", "underexpanded"], (g, ar, names)
    one = GD.laval_critical_pressures(1.0)
    assert one["p_b"] == one["p_d"] == one["p_f"] == one["p_star"]            # no diverging part: the three thresholds collapse on p*
    arr = GD.laval_critical_pressures(np.array([1.5, 3.0]))
    assert arr["p_d"].shape == (2,)
    for bad in (1.1, -0.1):
        with pytest.raises(ValueError):
            GD.laval_regime(bad, 2.0)


def _shock_in_nozzle_by_mass_conservation(r, ar, g):
    """Independent route to the shock in a Laval nozzle: chain of ratios for p02/p01, mass conservation for the exit Mach number."""
    def exit_pressure(M1):
        M2 = math.sqrt(((g - 1) * M1 ** 2 + 2) / (2 * g * M1 ** 2 + 1 - g))
        k = (1 + 2 * g / (g + 1) * (M1 ** 2 - 1)) * (t0_ratio(M2, g) / t0_ratio(M1, g)) ** (g / (g - 1))
        flux = lambda M: M * t0_ratio(M, g) ** (-(g + 1) / (2 * (g - 1)))    # noqa: E731
        Me = brentq(lambda M: flux(M) * k * ar - flux(1.0), 1e-9, 1.0, xtol=1e-15)        # same mdot at the throat (p01, A_t) and at the exit (p02, A_e)
        return k * t0_ratio(Me, g) ** (-g / (g - 1)), k, Me
    Msup = brentq(lambda M: t0_ratio(M, g) ** ((g + 1) / (2 * (g - 1))) / M / t0_ratio(1.0, g) ** ((g + 1) / (2 * (g - 1))) - ar, 1.0, 50.0, xtol=1e-15)
    M1 = brentq(lambda M: exit_pressure(M)[0] - r, 1.0 + 1e-12, Msup, xtol=1e-14)
    return (M1,) + exit_pressure(M1)[1:]


def test_laval_nozzle_shock_position_V2_derivation():  # V2 + V4 — D20 (two stars, ours): p0 A* is conserved, so A2* = A_t p01/p02; the shock found by mass conservation; C10
    g, R, T0, p01, p02, A1s, A2s = sp.symbols("gamma R T_0 p_01 p_02 A_1s A_2s", positive=True)
    K = sp.sqrt(g / (R * T0)) * (2 / (g + 1)) ** ((g + 1) / (2 * (g - 1)))    # mdot = K p0 A*  (D11) with T0 the same on both sides
    assert sp.solve(sp.Eq(K * p01 * A1s, K * p02 * A2s), A2s) == [A1s * p01 / p02]
    for g_ in (1.4, 1.25):
        for ar in (1.3, I["Ae_over_At"], 6.0):
            c = GD.laval_critical_pressures(ar, g_)
            for w in (0.03, 0.4, I["pB_over_p0_shock"] if (ar == 2.0 and g_ == 1.4) else 0.7, 0.97):
                r = c["p_d"] + w * (c["p_b"] - c["p_d"]) if w != I["pB_over_p0_shock"] else w
                s = GD.laval_shock_location(r, ar, g_)
                M1, k, Me = _shock_in_nozzle_by_mass_conservation(r, ar, g_)
                assert rel(s["M1"], M1) < 1e-8 and rel(s["p02_over_p01"], k) < 1e-9 and rel(s["M_exit"], Me) < 1e-9, (g_, ar, r)
                assert close(s["M2"], GD.normal_shock_mach(s["M1"], g_), 1e-14) and close(s["A_shock_over_At"], GD.area_ratio(s["M1"], g_), 1e-14)
                assert close(s["p_exit_over_p0"], r, 1e-9) and 1.0 <= s["A_shock_over_At"] <= ar * (1 + 1e-12) and s["M_exit"] < 1


def test_laval_nozzle_V4_consistent_with_shock_and_isentropic_relations():  # V4 + V7 — C10: mass flow along the nozzle, the two ends of the shock range, continuity
    ar, g = I["Ae_over_At"], 1.4
    c = GD.laval_critical_pressures(ar, g)
    top, bottom = GD.laval_shock_location(c["p_b"], ar, g), GD.laval_shock_location(c["p_d"], ar, g)
    assert abs(top["M1"] - 1) < 1e-6 and abs(top["A_shock_over_At"] - 1) < 1e-10 and rel(top["M_exit"], c["M_sub"]) < 1e-6        # shock at the throat as p_B → p_b
    assert rel(bottom["M1"], c["M_sup"]) < 1e-7 and rel(bottom["A_shock_over_At"], ar) < 1e-7                                    # … and in the exit plane at p_d
    rs = np.linspace(c["p_d"], c["p_b"], 41)
    M1s = np.array([GD.laval_shock_location(r, ar, g)["M1"] for r in rs])
    assert np.all(np.diff(M1s) < 0)                                           # a lower back pressure pulls the shock downstream
    for bad in (c["p_d"] * 0.99, c["p_b"] * 1.01):
        with pytest.raises(ValueError):
            GD.laval_shock_location(bad, ar, g)
    e = 1e-7                                                                   # every output of laval_exit_state is continuous across p_b
    hi, at, lo = GD.laval_exit_state(c["p_b"] * (1 + e), ar, g), GD.laval_exit_state(c["p_b"], ar, g), GD.laval_exit_state(c["p_b"] * (1 - e), ar, g)
    assert hi["regime"] == "subsonic" and not hi["choked"] and at["regime"] == "sonic throat" and lo["regime"] == "shock in nozzle" and lo["choked"]
    assert abs(hi["M_exit"] - lo["M_exit"]) < 1e-5 and abs(hi["mdot_over_choked"] - 1) < 1e-6 and abs(hi["M_throat"] - 1) < 3e-3 and at["M_exit"] == c["M_sub"]
    ex = GD.laval_exit_state(c["p_d"], ar, g)
    assert ex["regime"] == "shock at exit" and ex["M_exit"] == c["M_sup"] and ex["A_shock_over_At"] == ar
    for r, name in ((0.3, "overexpanded"), (c["p_f"], "design"), (0.02, "underexpanded")):
        st = GD.laval_exit_state(r, ar, g)
        assert st["regime"] == name and st["M_exit"] == c["M_sup"] and st["p_exit_over_p0"] == c["p_f"] and st["mdot_over_choked"] == 1.0 and math.isnan(st["M1_shock"])
    nf = GD.laval_exit_state(1.0, ar, g)
    assert nf["regime"] == "no flow" and nf["M_exit"] == 0.0 and nf["mdot_over_choked"] == 0.0 and not nf["choked"]
    x = np.linspace(0.0, 1.0, 801)                                             # throat x = 0.4 on the grid
    flux1 = GD.mass_flux_function(1.0, g)
    for r in (0.99, 0.95, c["p_b"], 0.85, I["pB_over_p0_shock"], 0.55, 0.3, c["p_f"], 0.02):
        n = GD.laval_nozzle(x, GD.laval_area, r, g)
        p0loc = n["p_over_p0"] * GD.stagnation_pressure_ratio(n["M"], g)       # local stagnation pressure over p01
        mdot = p0loc * GD.mass_flux_function(n["M"], g) * n["A"]               # ~ rho u A (T0 constant)
        assert np.ptp(mdot) < 1e-8 * flux1, (r, n["regime"], np.ptp(mdot))     # the same mass flow at every station, across the shock too
        assert rel(mdot[0], n["mdot_over_choked"] * flux1 * 1.0) < 1e-8 and n["x_throat"] == pytest.approx(0.4) and n["regime"] == GD.laval_regime(r, ar, g)
        if n["regime"] in ("subsonic", "sonic throat", "shock in nozzle"):
            assert abs(n["p_over_p0"][-1] - r) < 1e-8 and abs(n["p_exit_over_p0"] - r) < 1e-8 and n["M"][-1] < 1      # a subsonic exit takes the back pressure
        else:
            assert abs(n["p_over_p0"][-1] - c["p_f"]) < 1e-10 and n["M"][-1] > 1 and math.isnan(n["x_shock"])           # a supersonic exit does not feel it
        if n["regime"] == "shock in nozzle":
            xs = n["x_shock"]
            assert 0.4 < xs < 1.0 and rel(GD.laval_area(xs), n["A_shock_over_At"]) < 1e-10
            before, after = n["M"][(x > 0.4) & (x < xs)], n["M"][x > xs]
            assert np.all(before > 1) and np.all(after < 1) and close(np.unique(np.round(p0loc[x > xs], 12)), n["p02_over_p01"], 1e-10)
            arr = GD.laval_nozzle(x, GD.laval_area(x), r, g)                    # an array of areas instead of a callable: x_shock by interpolation
            assert abs(arr["x_shock"] - xs) < 2e-3 and close(arr["M_exit"], n["M_exit"], 1e-12)
            s = GD.laval_shock_location(r, ar, g)                               # and the dimensional chain run forwards (the worked example with our inputs)
            w = ch15.nozzle_with_shock(s["M1"], I["A_throat"], I["A_exit"], I["p0"], I["T0"], AIR)
            assert rel(w["p3"] / I["p0"], r) < 1e-8 and rel(w["M3"], s["M_exit"]) < 1e-9 and rel(w["A_shock"], s["A_shock_over_At"] * I["A_throat"]) < 1e-12
            assert rel(w["mdot"], w["mdot_exit"]) < 1e-11 and rel(w["A_star2"], I["A_throat"] / s["p02_over_p01"]) < 1e-9 and rel(w["p02"], I["p0"] * s["p02_over_p01"]) < 1e-9
            assert rel(w["T3"], I["T0"] / GD.stagnation_temperature_ratio(w["M3"])) < 1e-13 and close(w["M2"], s["M2"], 1e-12)
    with pytest.raises(ValueError):
        ch15.nozzle_with_shock(2.5, I["A_throat"], I["A_exit"], I["p0"], I["T0"])       # this shock would stand outside the nozzle


def test_rocket_nozzle_V1_thrust_and_throat_from_two_routes():  # V1 — N46, N47 with our own inputs; trap T9
    gas = GD.PerfectGas(I["gamma_rocket"], I["R_rocket"])
    n = ch15.rocket_nozzle(I["thrust"], I["p0_rocket"], I["T0_rocket"], I["p_exit_rocket"], gas)
    g = gas.gamma
    assert rel(n["mdot"] * n["u_e"], I["thrust"]) < 1e-13 and rel(g * I["p_exit_rocket"] * n["M_e"] ** 2 * n["A_e"], I["thrust"]) < 1e-12     # thrust = mdot u_e = gamma p_e M_e^2 A_e
    assert rel(n["mdot"], n["mdot_choked"]) < 1e-12 and rel(n["Ae_over_At"], GD.area_ratio(n["M_e"], g)) < 1e-14 and rel(n["A_throat"] * n["Ae_over_At"], n["A_e"]) < 1e-14
    assert rel(I["p0_rocket"] / I["p_exit_rocket"], GD.stagnation_pressure_ratio(n["M_e"], g)) < 1e-12 and rel(n["T_e"] * GD.stagnation_temperature_ratio(n["M_e"], g), I["T0_rocket"]) < 1e-14
    assert rel(n["rho_e"], I["p_exit_rocket"] / (gas.R * n["T_e"])) < 1e-14 and n["u_e"] < GD.max_speed(I["T0_rocket"], gas.cp) and n["M_e"] > 1
    assert rel(ch15.nozzle_thrust(n["mdot"], n["u_e"], I["p_exit_rocket"], I["p_exit_rocket"], n["A_e"]), I["thrust"]) < 1e-13                # matched: no pressure term
    vac = ch15.nozzle_thrust(n["mdot"], n["u_e"], I["p_exit_rocket"], 0.0, n["A_e"])
    assert rel(vac - I["thrust"], I["p_exit_rocket"] * n["A_e"]) < 1e-12                # off design the pressure term appears


def test_rocket_nozzle_V7_warns_without_a_gas():  # V7 — review Should-fix 3: gas=None falls back to air and says so (trap T9)
    args = (I["thrust"], I["p0_rocket"], I["T0_rocket"], I["p_exit_rocket"])
    with pytest.warns(UserWarning, match="no gas given"):
        air = ch15.rocket_nozzle(*args)
    with pytest.warns(UserWarning):
        ch15.rocket_nozzle(*args, gas=None)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error")                                     # with a gas: silent
        own = ch15.rocket_nozzle(*args, GD.PerfectGas(I["gamma_rocket"], I["R_rocket"]))
        ch15.rocket_nozzle(*args, gas=GD.PerfectGas())
    assert rel(air["A_throat"], ch15.rocket_nozzle(*args, gas=GD.PerfectGas())["A_throat"]) < 1e-15       # the fallback is the project air
    assert rel(air["A_throat"], own["A_throat"]) > 1e-3                     # and it is not the exhaust answer


def test_small_helpers_V1_against_independent_closed_forms():  # V1 — label audit after review Should-fix 8: duct_role, rayleigh_temperature_sign, strong_shock_limits, laval_area against closed forms they do not use
    gv = np.array(GAMMAS)
    lim = [GD.strong_shock_limits(g) for g in GAMMAS]                       # M1 → infinity of the jump laws (15.38)–(15.41), written out by hand
    assert close([d["rho_ratio"] for d in lim], (gv + 1) / (gv - 1), 1e-15) and close([d["M2"] for d in lim], np.sqrt((gv - 1) / (2 * gv)), 1e-15)
    assert close([d["p_ratio_over_M1sq"] for d in lim], 2 * gv / (gv + 1), 1e-15) and close([d["T_ratio_over_p_ratio"] for d in lim], (gv - 1) / (gv + 1), 1e-15)
    assert close([d["u1sq_rho1_over_p2"] for d in lim], (gv + 1) / 2, 1e-15)
    x = np.linspace(0.0, 1.0, 201)                                          # two parabolas joined at the throat (0.4, 1.0)
    exact = np.where(x <= 0.4, 1.0 + 1.5 * ((0.4 - x) / 0.4) ** 2, 1.0 + 1.0 * ((x - 0.4) / 0.6) ** 2)
    assert close(GD.laval_area(x), exact, 1e-15)
    for g in GAMMAS:                                                        # duct_role from the isentropic A/A*(M) curve, not from (15.34)
        for M in np.concatenate([np.linspace(0.05, 0.95, 19), np.linspace(1.05, 4.0, 20)]):
            h = 1e-6
            dA_dM = (GD.area_ratio(M + h, g) - GD.area_ratio(M - h, g)) / (2 * h)
            for s in (-1.0, 1.0):
                assert GD.duct_role(M, s) == ("nozzle" if s / dA_dM > 0 else "diffuser"), (g, M, s)   # nozzle = M (hence u) rises
            # Rayleigh line (constant area, heat only): T/T* = [(1 + g) M/(1 + g M^2)]^2, T0/T0* rises with heating
            TT = lambda m: ((1 + g) * m / (1 + g * m * m)) ** 2                                           # noqa: E731
            T0T0 = lambda m: (1 + g) * m * m * (2 + (g - 1) * m * m) / (1 + g * m * m) ** 2               # noqa: E731
            dT = (TT(M + h) - TT(M - h)) / (T0T0(M + h) - T0T0(M - h))      # dT/dT0 along the Rayleigh line = sign of dT under heating
            if abs(1 - g * M * M) > 1e-3:
                assert GD.rayleigh_temperature_sign(M, g) == np.sign(dT), (g, M)


def test_stagnation_pressure_V1_incompressible_wrapper_is_not_15_29():  # V1 — review Should-fix 1: the trap pinned — ch15.stagnation_pressure is ch04's p + rho U^2/2, ~5 % below the isentropic p0 at M = 0.8
    from fluidpy.core import bernoulli as BER
    assert ch15.stagnation_pressure_incompressible is ch15.stagnation_pressure and "stagnation_pressure_incompressible" not in ch15.__all__
    p, rho, g = 1.0e5, 1.2, 1.4
    M = np.array([0.0, 0.1, 0.3, 0.5, 0.8, 0.95])
    U = M * math.sqrt(g * p / rho)                                          # c^2 = gamma p / rho
    inc = ch15.stagnation_pressure(p, U, rho)
    assert close(inc, ch15.stagnation_pressure_incompressible(p, U, rho), 0) and close(inc, BER.stagnation_pressure(p, U, rho), 0)
    assert close(inc, p + 0.5 * rho * U ** 2, 1e-15)                        # p + rho U^2 / 2, i.e. p (1 + gamma M^2 / 2)
    p0 = p * GD.stagnation_pressure_ratio(M, g)                             # (15.29): p0/p = [1 + (gamma − 1) M^2 / 2]^(gamma/(gamma − 1))
    low = 1 - inc / p0
    assert close(low, 1 - (1 + 0.5 * g * M ** 2) / (1 + 0.5 * (g - 1) * M ** 2) ** (g / (g - 1)), 1e-12, 1e-16)
    assert low[0] == 0.0 and np.all(np.diff(low) > 0)                       # agree at rest, the gap grows with M
    assert 0.049 < low[4] < 0.051 and abs(low[4] - 0.050081) < 1e-6         # about 5 % low at M = 0.8 (1.448 vs 1.5243): the trap being pinned
    assert low[2] < 0.0015                                                # below 1 % at M = 0.3 (the usual incompressibility rule)
    assert close((p0 - p)[1:] / (inc - p)[1:], GD.pitot_compressibility_factor(M[1:], g), 1e-12)   # the ratio of the two excesses


# =====================================================================================================================
# 6. §15.8 — C11 the duct with friction and heat (15.45)–(15.46), C12 the two theorems, the choking limits
# =====================================================================================================================
DUCT_CASES = [(g, M1, wf, wq) for g in (1.4, 1.25) for M1 in (0.3, 0.6, 0.9, 1.4, 2.5) for wf, wq in ((0.6, 0.0), (0.0, 0.6), (0.3, 0.3), (0.2, -0.5), (0.0, -1.0))]


def test_duct_outlet_mach_V4_the_three_laws_hold_on_both_roots():  # V4 — (15.45) with fluxes written here, both roots of (15.46); C11
    n_roots = 0
    for g, M1, wf, wq in DUCT_CASES:
        lim = GD.duct_choking_limits(M1, g)
        f, q = wf * lim["f_max"], wq * lim["q_max"]
        A = GD.duct_parameter_A(M1, f, q, g)
        assert close(A, M1 ** 2 * (1 + 0.5 * (g - 1) * M1 ** 2 + q) / (1 + g * M1 ** 2 - f) ** 2, 1e-14)
        for branch in ("continuous", "shocked"):
            M2 = GD.duct_outlet_mach(M1, f, q, g, branch=branch)
            if math.isnan(M2):
                continue
            n_roots += 1
            r = GD.duct_state_ratios(M1, M2, f, q, g)
            s1, s2 = state(M1, 1.0, 1.0, g), state(M2, r["p_ratio"], r["T_ratio"], g)
            F1, F2, h1 = fluxes(s1, g), fluxes(s2, g), g / (g - 1)
            assert abs(F2[0] / F1[0] - 1) < 1e-11, (g, M1, f, q, branch)                       # rho1 u1 = rho2 u2
            assert abs(F2[1] + f * 1.0 - F1[1]) < 1e-11 * F1[1]                                # p1 + rho1 u1^2 = p2 + rho2 u2^2 + p1 f
            assert abs(F2[2] - F1[2] - h1 * q) < 1e-11 * F1[2]                                 # h1 + u1^2/2 + h1 q = h2 + u2^2/2
            res = GD.duct_residuals(s1, s2, f, q, g)
            assert max(abs(res[k]) for k in res) < 1e-11
            assert close(M2 ** 2 * (1 + 0.5 * (g - 1) * M2 ** 2) / (1 + g * M2 ** 2) ** 2, A, 1e-11)       # the outlet side of the definition of A
            assert close(r["rho_ratio"] * r["u_ratio"], 1.0, 1e-13) and close(r["T0_ratio"], 1 + q / t0_ratio(M1, g), 1e-13)
            assert close(r["p0_ratio"], r["p_ratio"] * GD.stagnation_pressure_ratio(M2, g) / GD.stagnation_pressure_ratio(M1, g), 1e-12)
            dS = ch15.perfect_gas_entropy_change(1.0, 1.0, r["T_ratio"], 1.0 / r["rho_ratio"], cv=1 / (g - 1), R=1.0) * (g - 1)
            assert abs(r["dS_over_cv"] - dS) < 1e-11 and close(r["dS_over_cp"], r["dS_over_cv"] / g, 1e-14) and close(GD.duct_entropy_change(M1, M2, f, q, g), r["dS_over_cv"], 0, 1e-15)
            if q == 0.0 and f > 0:
                assert (r["dS_over_cv"] > 0) == (branch == "continuous" or M1 > 1) and r["p0_ratio"] < 1 or branch == "shocked"   # friction: entropy up, p0 down
    assert n_roots >= 80
    bad = GD.duct_residuals(state(0.5), state(0.6), 0.1, 0.0)
    assert abs(bad["momentum"]) > 1e-3


def test_duct_outlet_mach_V2_derivation():  # V2 — D21 (three stars), every line of Part F: (15.45) → A → the quadratic → (15.46); C11
    g, M1, M2, f, q, R = sp.symbols("gamma M_1 M_2 f q R", positive=True)
    p1, p2, T1, T2 = sp.symbols("p_1 p_2 T_1 T_2", positive=True)
    rho = lambda p, T: p / (R * T)                                             # noqa: E731
    u = lambda M, T: M * sp.sqrt(g * R * T)                                    # noqa: E731
    h = lambda T: g * R / (g - 1) * T                                          # noqa: E731
    mom = p1 + rho(p1, T1) * u(M1, T1) ** 2 - (p2 + rho(p2, T2) * u(M2, T2) ** 2 + p1 * f)
    assert z0(mom - (p1 * (1 + g * M1 ** 2 - f) - p2 * (1 + g * M2 ** 2)))     # step 1: rho u^2 = gamma p M^2
    en = h(T1) + u(M1, T1) ** 2 / 2 + h(T1) * q - (h(T2) + u(M2, T2) ** 2 / 2)
    assert z0(en - (h(T1) * (1 + (g - 1) / 2 * M1 ** 2 + q) - h(T2) * (1 + (g - 1) / 2 * M2 ** 2)))      # step 2
    mass = rho(p1, T1) * u(M1, T1) - rho(p2, T2) * u(M2, T2)
    assert z0(mass - sp.sqrt(g / R) * (p1 * M1 / sp.sqrt(T1) - p2 * M2 / sp.sqrt(T2)))                   # step 3: the mass flux is p M/sqrt(T)
    pr = (1 + g * M1 ** 2 - f) / (1 + g * M2 ** 2)                             # p2/p1 (step 5 inverted)
    Tr = (1 + (g - 1) / 2 * M1 ** 2 + q) / (1 + (g - 1) / 2 * M2 ** 2)         # T2/T1 (step 6)
    lhs = (M2 / M1) ** 2 - (1 / pr) ** 2 * Tr                                  # step 4 squared: (M2/M1)^2 = (p1/p2)^2 T2/T1
    A_in = M1 ** 2 * (1 + (g - 1) / 2 * M1 ** 2 + q) / (1 + g * M1 ** 2 - f) ** 2
    A_out = M2 ** 2 * (1 + (g - 1) / 2 * M2 ** 2) / (1 + g * M2 ** 2) ** 2
    assert z0(sp.together(lhs) * M1 ** 2 * (1 + (g - 1) / 2 * M2 ** 2) / (1 + g * M2 ** 2) ** 2 - (A_out - A_in))        # step 7: the two sides of A
    A, X = sp.symbols("A X", positive=True)
    quad_ = sp.expand(X * (1 + (g - 1) / 2 * X) - A * (1 + g * X) ** 2)        # step 8: X = M2^2
    assert z0(quad_ - (((g - 1) / 2 - A * g ** 2) * X ** 2 + (1 - 2 * A * g) * X - A))                   # step 9
    a2, a1, a0 = (g - 1) / 2 - A * g ** 2, 1 - 2 * A * g, -A
    disc = sp.expand(a1 ** 2 - 4 * a2 * a0)
    assert z0(disc - (1 - 2 * A * (g + 1)))                                    # step 11: the discriminant (printed by the design's check cell)
    for sign in (+1, -1):                                                      # steps 10, 12: (15.46) solves the quadratic for either sign
        root = (-(1 - 2 * A * g) + sign * sp.sqrt(1 - 2 * A * (g + 1))) / ((g - 1) - 2 * A * g ** 2)
        assert z0(quad_.subs(X, root))
    A0 = A_in.subs({f: 0, q: 0})                                               # step 13: f = q = 0
    assert z0(disc.subs(A, A0) - (M1 ** 2 - 1) ** 2 / (g * M1 ** 2 + 1) ** 2)
    assert z0(quad_.subs({A: A0, X: M1 ** 2})) and z0(quad_.subs({A: A0, X: ((g - 1) * M1 ** 2 + 2) / (2 * g * M1 ** 2 + 1 - g)}))      # the roots: M1^2 and (15.40)
    assert z0(disc.subs(A, 1 / (2 * (g + 1)))) and z0(quad_.subs({A: 1 / (2 * (g + 1)), X: 1}))          # one root, M2 = 1, when A = 1/(2(gamma+1))
    assert abs(GD.duct_parameter_A(0.5, 0.1, 0.0) - 0.168) < 5e-4 and abs(1 - 2 * GD.duct_parameter_A(0.5, 0.1, 0.0) * 2.4 - 0.1936) < 5e-4       # Part F's numbers
    assert abs(GD.duct_outlet_mach(0.5, 0.1) - 0.589) < 5e-4 and abs(GD.duct_outlet_mach(0.5, 0.1, branch="shocked") - 1.936) < 5e-4


def test_duct_outlet_mach_V1_without_friction_or_heat_returns_the_inlet_or_the_shock():  # V1 + V7 — C11, trap T12: which root continues M2 = M1
    for g in GAMMAS:
        M = np.concatenate([M_SUB, M_SUP])
        c, s = GD.duct_outlet_mach(M, 0.0, 0.0, g, branch="both")
        assert close(c, M, 1e-12)                                              # the continuous root: nothing happens
        assert close(s[M > 1], GD.normal_shock_mach(M[M > 1], g), 1e-12)       # the other root for M1 > 1: behind a normal shock (15.40)
        ok = M[M < 1] > math.sqrt((g - 1) / (2 * g))
        assert close(s[M < 1][ok], GD.normal_shock_mach(M[M < 1][ok], g), 1e-11) and np.all(np.isnan(s[M < 1][~ok]))
        dS = GD.duct_entropy_change(M, s, 0.0, 0.0, g)
        assert np.all(dS[M > 1] > 0) and np.all(dS[M < 1][ok] < 0) and close(dS[M > 1], GD.normal_shock_entropy(M[M > 1], g), 1e-9)     # the excluded root loses entropy
        assert np.max(np.abs(GD.duct_entropy_change(M, c, 0.0, 0.0, g))) < 1e-12
        plus, minus = GD.duct_outlet_mach(M, 0.0, 0.0, g, branch="plus"), GD.duct_outlet_mach(M, 0.0, 0.0, g, branch="minus")
        assert close(plus[M < 1], M[M < 1], 1e-12) and close(minus[M > 1], M[M > 1], 1e-12)      # "+" continues the subsonic inlet, "−" the supersonic one
        assert np.max(np.abs(plus[M > 1] - M[M > 1])) > 0.03                   # planted: one sign for all M1 is not continuous with M2 = M1
        assert np.all(np.nan_to_num(np.abs(minus[M < 1] - M[M < 1]), nan=1.0) > 0.03)
    assert GD.duct_outlet_mach(1.0, branch="both") == (1.0, 1.0)
    cur = GD.duct_curve(np.linspace(0.2, 3.5, 34), f=0.1)
    assert set(cur) == {"M1", "M2_continuous", "M2_shocked", "choked", "dS_continuous", "dS_shocked"} and cur["choked"].any() and not cur["choked"].all()
    assert np.array_equal(cur["choked"], GD.duct_choking_limits(cur["M1"])["f_max"] < 0.1)             # choked exactly where f exceeds f_max(M1)
    sub = (cur["M1"] < 1) & ~cur["choked"]
    assert np.all(cur["M2_continuous"][sub] > cur["M1"][sub]) and np.all(cur["M2_continuous"][(cur["M1"] > 1) & ~cur["choked"]] < cur["M1"][(cur["M1"] > 1) & ~cur["choked"]])
    with pytest.raises(ValueError):
        GD.duct_outlet_mach(0.5, branch="upper")
    assert math.isnan(GD.duct_parameter_A(0.5, f=1.5)) and math.isnan(GD.duct_outlet_mach(0.5, f=1.5))     # the wall force exceeds the inlet impulse


def test_duct_outlet_mach_V7_removable_singularity_of_the_printed_form():  # V7 — (15.46) is 0/0 where 2 A gamma^2 = gamma − 1; the code takes the limit
    g = 1.4
    for f in (0.0, 0.1):
        M0 = brentq(lambda M: GD.duct_parameter_A(M, f, 0.0, g) - (g - 1) / (2 * g * g), 0.05, 0.95, xtol=1e-15)
        A = GD.duct_parameter_A(M0, f, 0.0, g)
        b, D, den = 1 - 2 * A * g, 1 - 2 * A * (g + 1), (g - 1) - 2 * A * g * g
        assert abs(den) < 1e-14 and abs(-b + math.sqrt(D)) < 1e-14             # numerator and denominator of the "+" root vanish together
        M2 = GD.duct_outlet_mach(M0, f, 0.0, g)
        assert math.isfinite(M2) and close(M2 ** 2 * (1 + 0.2 * M2 ** 2) / (1 + 1.4 * M2 ** 2) ** 2, A, 1e-12)        # a true root there
        h = 1e-6
        left, right = GD.duct_outlet_mach(M0 - h, f, 0.0, g), GD.duct_outlet_mach(M0 + h, f, 0.0, g)
        assert abs(left + right - 2 * M2) < 1e-9                               # and it joins its neighbours smoothly
        if f == 0.0:
            assert close(M0, math.sqrt((g - 1) / (2 * g)), 1e-12) and close(M2, M0, 1e-12)      # where the shock partner (15.40) runs off to infinity
            assert math.isnan(GD.duct_outlet_mach(M0 * (1 - 1e-9), branch="shocked")) and GD.duct_outlet_mach(M0 * (1 + 1e-6), branch="shocked") > 500


def test_duct_choking_limits_V1_outlet_is_sonic_at_the_limit_and_no_root_beyond():  # V1 + V7 — C12, D24 (ours): f_max, q_max; behaviour exactly at the limit
    worst = 0.0
    for g in (1.4, 1.25, 5.0 / 3.0):
        bound = 0.5 * (g + 1) * math.sqrt(4 * 2.220446049250313e-16) * 1.01    # |M2 − 1| = ((gamma+1)/2) sqrt(D) and round-off leaves |D| <= a few ulp
        for M1 in (0.2, 0.35, 0.5, 0.7, 0.9, 0.97, 1.03, 1.3, 2.0, 2.7, 3.5):
            lim = GD.duct_choking_limits(M1, g)
            assert lim["f_max"] > 0 and lim["q_max"] > 0 and lim["f_max"] < 1 + g * M1 ** 2
            for kw in (dict(f=lim["f_max"]), dict(q=lim["q_max"])):
                A = GD.duct_parameter_A(M1, gamma=g, **kw)
                assert abs(2 * A * (g + 1) - 1) < 9e-16                         # the discriminant of (15.46) vanishes (to 4 ulp)
                for M2 in GD.duct_outlet_mach(M1, gamma=g, branch="both", **kw):
                    worst = max(worst, abs(M2 - 1))
                    assert abs(M2 - 1) < bound, (g, M1, kw, M2 - 1)             # 1 to within the square root of round-off (the fold is a branch point)
                    assert abs(M2 ** 2 * (1 + 0.5 * (g - 1) * M2 ** 2) / (1 + g * M2 ** 2) ** 2 * 2 * (g + 1) - 1) < 1e-14       # backward error: it solves A(M2) = A to round-off
                key = "f" if "f" in kw else "q"
                over = {key: kw[key] * (1 + 1e-9)}
                assert all(math.isnan(v) for v in GD.duct_outlet_mach(M1, gamma=g, branch="both", **over))      # just beyond: no real root (choked)
                under = {key: kw[key] * (1 - 1e-6)}
                c, s = GD.duct_outlet_mach(M1, gamma=g, branch="both", **under)
                assert (c < 1) == (M1 < 1) and (s > 1) == (M1 < 1) and abs(c - 1) < 5e-3                         # just inside: two roots, one on each side of 1
        one = GD.duct_choking_limits(1.0, g)
        assert abs(one["f_max"]) < 1e-15 and abs(one["q_max"]) < 1e-15 and all(math.isnan(v) for v in GD.duct_outlet_mach(1.0, 1e-6, 0.0, g, branch="both"))
    assert 1e-9 < worst < 4e-8                                                  # measured: 1.3e-8 — honest, and what the docstring of duct_outlet_mach implies
    loc = GD.duct_choking_locus(1.4)
    assert loc["M1"].shape == (400,) and close(loc["f_max"], GD.duct_choking_limits(loc["M1"])["f_max"], 0, 0) and close(GD.duct_choking_locus(1.4, [0.5, 2.0])["q_max"], [0.46875, 0.46875], 1e-12)


def test_duct_choking_limits_V2_derivation():  # V2 — D24 (two stars, ours): the vanishing discriminant of (15.46) solved for f and for q; C12
    g, M, f, q = sp.symbols("gamma M_1 f q", positive=True)
    A = M ** 2 * (1 + (g - 1) / 2 * M ** 2 + q) / (1 + g * M ** 2 - f) ** 2
    fs = sp.solve(sp.Eq(2 * A.subs(q, 0) * (g + 1), 1), f)
    target_f = 1 + g * M ** 2 - M * sp.sqrt(2 * (g + 1) * (1 + (g - 1) / 2 * M ** 2))
    assert len(fs) == 2 and any(z0(r - target_f) for r in fs)
    other = [r for r in fs if not z0(r - target_f)][0]
    assert z0(other - (1 + g * M ** 2) - (1 + g * M ** 2 - target_f))           # the other root has f > 1 + gamma M1^2: excluded
    qs = sp.solve(sp.Eq(2 * A.subs(f, 0) * (g + 1), 1), q)
    assert len(qs) == 1 and z0(qs[0] - ((1 + g * M ** 2) ** 2 / (2 * (g + 1) * M ** 2) - 1 - (g - 1) / 2 * M ** 2))
    assert z0(target_f.subs(M, 1)) and z0(qs[0].subs(M, 1))                     # both limits vanish at M1 = 1
    assert z0(sp.factor(qs[0]) - (M ** 2 - 1) ** 2 / (2 * (g + 1) * M ** 2))    # q_max = (M1^2 − 1)^2/(2(gamma+1)M1^2) >= 0
    Ac = 1 / (2 * (g + 1))                                                      # at the limit M2 = 1: X = −(1 − 2 A gamma)/((gamma−1) − 2 A gamma^2)
    assert z0(-(1 - 2 * Ac * g) / ((g - 1) - 2 * Ac * g ** 2) - 1)


def test_fanno_rayleigh_V5_cited_reference_state_forms():  # V5 (forms) + V1 — N86–N89 (ours): the usual starred ratios; trap T2
    b = bench()
    assert "4 f L*/D_h" in b["fanno_flow"]["length"] and "T0/T0*" in b["rayleigh_flow"]["ratios"]
    for g in GAMMAS:
        M = np.array([0.2, 0.5, 0.9, 1.0, 1.4, 2.5, 4.0])
        B = 2 / (g + 1) * (1 + 0.5 * (g - 1) * M ** 2)
        fa = GD.fanno_ratios(M, g)
        assert close(fa["p"], 1 / (M * np.sqrt(B)), 1e-13) and close(fa["rho"], np.sqrt(B) / M, 1e-13) and close(fa["T"], 1 / B, 1e-13) and close(fa["u"], M / np.sqrt(B), 1e-13)
        assert close(fa["p0"], B ** ((g + 1) / (2 * (g - 1))) / M, 1e-13)
        assert close(fa["fL"], (1 - M ** 2) / (g * M ** 2) + (g + 1) / (2 * g) * np.log(M ** 2 / B), 1e-12, 1e-15) and np.all(fa["fL"] >= -1e-15)
        ln = GD.fanno_line_Ts(M, g)
        assert close(ln["ds_over_cp"], np.log(M ** ((g - 1) / g) * B ** (-(g + 1) / (2 * g))), 1e-12, 1e-15) and close(ln["T_over_Tstar"], fa["T"], 0) and np.all(ln["ds_over_cp"] <= 1e-15)
        ra = GD.rayleigh_ratios(M, g)
        assert close(ra["p"], (g + 1) / (1 + g * M ** 2), 1e-13) and close(ra["rho"], (1 + g * M ** 2) / ((g + 1) * M ** 2), 1e-13) and close(ra["T"], (g + 1) ** 2 * M ** 2 / (1 + g * M ** 2) ** 2, 1e-13)
        assert close(ra["u"], (g + 1) * M ** 2 / (1 + g * M ** 2), 1e-13) and close(ra["T0"], 2 * (g + 1) * M ** 2 * (1 + 0.5 * (g - 1) * M ** 2) / (1 + g * M ** 2) ** 2, 1e-13)
        assert close(ra["p0"], (g + 1) / (1 + g * M ** 2) * B ** (g / (g - 1)), 1e-13)
        for d in (fa, ra):
            assert all(abs(float(np.asarray(v)[3]) - (0.0 if k == "fL" else 1.0)) < 1e-14 for k, v in d.items())      # every ratio is 1 (and fL = 0) at M = 1
        rl = GD.rayleigh_line_Ts(np.linspace(0.3, 3.0, 2701), g)
        assert abs(rl["M"][np.argmax(rl["T_over_Tstar"])] - 1 / math.sqrt(g)) < 1.5e-3 and abs(rl["M"][np.argmax(rl["ds_over_cp"])] - 1.0) < 1.5e-3 and np.all(rl["ds_over_cp"] <= 1e-15)
        # the cited differentials against the influence coefficients: dM^2/M^2 per dT0/T0 (Rayleigh) and per 4f dx/D (Fanno, df = (gamma M^2/2) d(4fx/D))
        for Mv in (0.4, 0.8, 1.7):
            Tf = t0_ratio(Mv, g)
            assert close(2 * GD.influence_coefficients(Mv, g, dq=1.0)["dM_over_M"] * Tf, (1 + g * Mv ** 2) / (1 - Mv ** 2) * Tf, 1e-12)
            assert close(2 * GD.influence_coefficients(Mv, g, df=1.0)["dM_over_M"] * 0.5 * g * Mv ** 2, g * Mv ** 2 / (1 - Mv ** 2) * Tf, 1e-12)
            num = quad(lambda m: 1.0 / (m * GD.influence_coefficients(m, g, df=1.0)["dM_over_M"] * 0.5 * g * m * m), Mv, 1.0, epsabs=1e-13, epsrel=1e-12)[0]
            assert abs(num - GD.fanno_ratios(Mv, g)["fL"]) < 1e-9               # the friction length is the integral of the corrected (15.47)
    assert abs(1 / math.sqrt(1.4) - b["rayleigh_flow"]["mach_of_maximum_static_enthalpy"]) < 5e-4      # the cited 0.845 (gamma = 1.4)


def test_fanno_and_rayleigh_V1_bridge_to_the_book_parameters_f_and_q():  # V1 — N88, trap T11: friction lengths add, the lumped f does not
    for g in (1.4, 1.3):
        for M1, Mm, M2 in ((0.4, 0.55, 0.7), (2.5, 2.0, 1.6), (0.2, 0.6, 0.95)):
            b = GD.fanno_length_to_f(M1, M2, g)
            assert rel(GD.duct_outlet_mach(M1, b["f"], 0.0, g), M2) < 1e-10 and b["f"] > 0 and b["friction_length"] > 0
            b1, b2 = GD.fanno_length_to_f(M1, Mm, g), GD.fanno_length_to_f(Mm, M2, g)
            assert close(b1["friction_length"] + b2["friction_length"], b["friction_length"], 1e-12) and abs(b1["f"] + b2["f"] - b["f"]) > 1e-3 * b["f"]
            r = GD.duct_state_ratios(M1, M2, b["f"], 0.0, g)
            fa1, fa2 = GD.fanno_ratios(M1, g), GD.fanno_ratios(M2, g)
            assert close([r["p_ratio"], r["T_ratio"], r["rho_ratio"], r["p0_ratio"]], [fa2["p"] / fa1["p"], fa2["T"] / fa1["T"], fa2["rho"] / fa1["rho"], fa2["p0"] / fa1["p0"]], 1e-10)
            ra1, ra2 = GD.rayleigh_ratios(M1, g), GD.rayleigh_ratios(M2, g)      # heating from M1 to M2: q from the T0 ratio
            q = (ra2["T0"] / ra1["T0"] - 1) * t0_ratio(M1, g)
            assert rel(GD.duct_outlet_mach(M1, 0.0, q, g), M2) < 1e-10 and q > 0
            rq = GD.duct_state_ratios(M1, M2, 0.0, q, g)
            assert close([rq["p_ratio"], rq["T_ratio"], rq["u_ratio"], rq["p0_ratio"], rq["T0_ratio"]], [ra2["p"] / ra1["p"], ra2["T"] / ra1["T"], ra2["u"] / ra1["u"], ra2["p0"] / ra1["p0"], ra2["T0"] / ra1["T0"]], 1e-10)
    assert close(GD.duct_choking_limits(0.5)["q_max"], (1 / GD.rayleigh_ratios(0.5)["T0"] - 1) * t0_ratio(0.5), 1e-12)      # q_max is the heat that reaches T0*


def test_friction_and_heating_theorems_V1_drive_the_mach_number_towards_one():  # V1 — C12: the two "surprises" as finite changes of (15.46); T18
    for g in (1.4, 5.0 / 3.0):
        for M1 in (0.2, 0.5, 0.8, 0.95, 1.1, 1.6, 3.0):
            lim = GD.duct_choking_limits(M1, g)
            fs, qs = np.linspace(0, 0.98, 15) * lim["f_max"], np.linspace(0, 0.98, 15) * lim["q_max"]
            Mf, Mq = GD.duct_outlet_mach(M1, fs, 0.0, g), GD.duct_outlet_mach(M1, 0.0, qs, g)
            sign = 1.0 if M1 < 1 else -1.0
            assert np.all(sign * np.diff(Mf) > 0) and np.all(sign * np.diff(Mq) > 0) and np.all(sign * (1 - Mf) > 0) and np.all(sign * (1 - Mq) > 0)     # towards 1, never across
            rf, rq = GD.duct_state_ratios(M1, Mf, fs, 0.0, g), GD.duct_state_ratios(M1, Mq, 0.0, qs, g)
            assert np.all(np.diff(rf["p0_ratio"]) < 0) and np.all(sign * np.diff(rf["p_ratio"]) < 0) and np.all(np.diff(rf["dS_over_cv"]) > 0)             # friction: p0 always falls
            assert np.all(np.diff(rq["p0_ratio"]) < 0) and np.all(np.diff(rq["T0_ratio"]) > 0)
            Mc = GD.duct_outlet_mach(M1, 0.0, -0.2, g)                           # cooling drives the Mach number away from 1
            assert sign * (Mc - M1) < 0
        lo = 1 / math.sqrt(g)
        Mb = 0.5 * (lo + 1)
        qb = 0.2 * GD.duct_choking_limits(Mb, g)["q_max"]
        assert GD.duct_state_ratios(Mb, GD.duct_outlet_mach(Mb, 0.0, qb, g), 0.0, qb, g)["T_ratio"] < 1          # heat added, temperature down (1/gamma < M^2 < 1)
        qa = 0.2 * GD.duct_choking_limits(0.5 * lo, g)["q_max"]
        assert GD.duct_state_ratios(0.5 * lo, GD.duct_outlet_mach(0.5 * lo, 0.0, qa, g), 0.0, qa, g)["T_ratio"] > 1


# =====================================================================================================================
# 7. §15.9 — C13 the Mach angle (15.49), C14 the oblique shock (15.50)–(15.55), C15 the Prandtl–Meyer function (15.56)
# =====================================================================================================================
def test_mach_angle_and_prandtl_meyer_V5_naca_report_1135_rows():  # V5 — NACA Report 1135 Table II, columns nu and mu (degrees), read first-hand; C13, C15
    t2 = bench()["naca1135_table2_supersonic"]
    k_nu, k_mu = t2["columns"].index("nu_deg"), t2["columns"].index("mu_deg")
    half_misses = []
    for row in t2["rows"]:
        d_nu = abs(math.degrees(GD.prandtl_meyer(row[0], 1.4)) - row[k_nu])
        assert d_nu <= 0.001 * 1.001, (row[0], math.degrees(GD.prandtl_meyer(row[0])), row[k_nu])          # one unit of the last printed digit
        if d_nu > 0.0005 * 1.001:
            half_misses.append(row[0])
        assert abs(math.degrees(GD.mach_angle(row[0])) - row[k_mu]) <= 0.005 * 1.001, (row[0], math.degrees(GD.mach_angle(row[0])), row[k_mu])
        assert rel(GD.mach_from_prandtl_meyer(math.radians(row[k_nu]), 1.4), row[0]) < 2e-5 and rel(GD.mach_from_angle(math.radians(row[k_mu])), row[0]) < 4e-4
    assert half_misses == [1.6]                                              # 5 of 6 to half a unit; at M = 1.60 the report prints 14.861 where (15.56) gives 14.86035
    assert abs(math.degrees(math.asin(1 / 1.4)) - bench()["naca1135_forms"]["wedge_limit_deg_gamma_7_5"]) < 0.05


def _oblique_2d(M1, o, g):
    """Mass, normal-momentum, tangential-velocity and energy mismatches of an oblique shock, from velocity vectors in two dimensions."""
    sig, dlt = o["sigma"], o["delta"]
    V1, V2 = M1 * math.sqrt(g), o["M2"] * math.sqrt(g * o["T_ratio"])           # p1 = rho1 = T1 = 1, R = 1
    n, t = np.array([math.sin(sig), -math.cos(sig)]), np.array([math.cos(sig), math.sin(sig)])
    v1, v2 = np.array([V1, 0.0]), V2 * np.array([math.cos(dlt), math.sin(dlt)])  # the flow is turned by delta towards the shock
    un1, un2, ut1, ut2 = v1 @ n, v2 @ n, v1 @ t, v2 @ t
    h1, h2 = g / (g - 1), g / (g - 1) * o["T_ratio"]
    return (o["rho_ratio"] * un2 / un1 - 1, (o["p_ratio"] + o["rho_ratio"] * un2 ** 2) / (1 + un1 ** 2) - 1, ut2 / ut1 - 1 if abs(ut1) > 1e-9 else ut2 - ut1,
            (h2 + 0.5 * V2 ** 2) / (h1 + 0.5 * V1 ** 2) - 1)


def test_oblique_shock_V4_jump_laws_in_two_dimensions():  # V4 — mass, normal momentum, tangential momentum, energy with velocity vectors; C14, D26
    n = 0
    for g in (1.4, 1.25, 5.0 / 3.0):
        for M1 in (1.15, 2.0, I["M_wedge"], 5.0, 12.0):
            mu = math.asin(1 / M1)
            for sig in np.linspace(mu * 1.001, 0.5 * PI * 0.999, 9):
                o = GD.oblique_shock(M1, sigma=float(sig), gamma=g)
                assert max(abs(v) for v in _oblique_2d(M1, o, g)) < 2e-11, (g, M1, sig, _oblique_2d(M1, o, g))
                c = GD.oblique_shock_components(M1, float(sig), g)
                assert close([c["Mn1"], c["Mt1"], c["u1_over_V1"], c["v_over_V1"]], [M1 * math.sin(sig), M1 * math.cos(sig), math.sin(sig), math.cos(sig)], 1e-14)
                assert close(c["Mn2"], GD.normal_shock_mach(c["Mn1"], g), 1e-14) and close(c["M2"], c["Mn2"] / math.sin(sig - c["delta"]), 1e-13) and c["Mn2"] < 1 < c["Mn1"]
                assert close(c["u2_over_V1"], math.sin(sig) / o["rho_ratio"], 1e-13) and close(c["V2_over_V1"], o["M2"] * math.sqrt(o["T_ratio"]) / M1, 1e-11) and close(c["c2_over_c1"], math.sqrt(o["T_ratio"]), 1e-14)
                assert close(o["p_ratio"], 1 + 2 * g / (g + 1) * (M1 ** 2 * math.sin(sig) ** 2 - 1), 1e-13) and close(o["rho_ratio"], math.tan(sig) / math.tan(sig - o["delta"]), 1e-10)   # (15.50), (15.51)
                assert o["branch"] == "given" and o["attached"] is True and o["dS_over_cv"] > 0 and o["p0_ratio"] < 1
                n += 1
    assert n == 135
    for M1, dl in ((2.0, 0.1), (I["M_wedge"], I["wedge_half_angle"]), (4.0, 0.5)):
        for br in ("weak", "strong"):
            assert max(abs(v) for v in _oblique_2d(M1, GD.oblique_shock(M1, delta=dl, branch=br), 1.4)) < 2e-10


def test_oblique_shock_V5_naca_report_1135_forms():  # V5 (forms) — the report's equations (124), (132), (137), (138), (139b), (142), (145), (167), (168), (170); C14
    forms = bench()["naca1135_forms"]
    assert "cot(delta)" in forms["eq_138"] and "Prandtl" in forms["eq_124"]
    for g in (1.4, 1.3, 5.0 / 3.0):
        for M1 in (1.3, 2.0, 3.7, 9.0):
            mu = math.asin(1 / M1)
            th = np.linspace(mu * 1.002, 0.5 * PI * 0.998, 21)                   # the report's theta is our sigma
            c = GD.oblique_shock_components(M1, th, g)
            dl = GD.deflection_angle(M1, th, g)
            m = M1 ** 2 * np.sin(th) ** 2
            assert close(1 / np.tan(dl), np.tan(th) * ((g + 1) * M1 ** 2 / (2 * (m - 1)) - 1), 1e-10)                                           # (138)
            assert close(np.tan(dl), (M1 ** 2 * np.sin(2 * th) - 2 / np.tan(th)) / (2 + M1 ** 2 * (g + np.cos(2 * th))), 1e-11)                 # (139b)
            M2sq = ((g + 1) ** 2 * M1 ** 4 * np.sin(th) ** 2 - 4 * (m - 1) * (g * m + 1)) / ((2 * g * m - (g - 1)) * ((g - 1) * m + 2))
            assert close(c["M2"] ** 2, M2sq, 1e-10)                                                                                             # (132): M2 by another route
            assert close(c["V2_over_V1"] ** 2, 1 - 4 * (m - 1) * (g * m + 1) / ((g + 1) ** 2 * M1 ** 4 * np.sin(th) ** 2), 1e-11)               # (137)
            a_star_sq = 2 / (g + 1) * (1 + 0.5 * (g - 1) * M1 ** 2) / M1 ** 2    # a*^2/V1^2
            assert close(c["u1_over_V1"] * c["u2_over_V1"], a_star_sq - (g - 1) / (g + 1) * c["v_over_V1"] ** 2, 1e-12)                         # (124): Prandtl's relation
            o = GD.oblique_shock(M1, sigma=th, gamma=g)
            assert close(o["p0_ratio"], ((g + 1) * m / ((g - 1) * m + 2)) ** (g / (g - 1)) * ((g + 1) / (2 * g * m - (g - 1))) ** (1 / (g - 1)), 1e-10)    # (142)
            assert close((o["p_ratio"] - 1) / (0.5 * g * M1 ** 2), 4 * (m - 1) / ((g + 1) * M1 ** 2), 1e-12)                                    # (145)
            s167 = ((g + 1) * M1 ** 2 - (3 - g) + math.sqrt((g + 1) * ((g + 1) * M1 ** 4 - 2 * (3 - g) * M1 ** 2 + (g + 9)))) / (4 * g * M1 ** 2)
            s168 = ((g + 1) * M1 ** 2 - 4 + math.sqrt((g + 1) * ((g + 1) * M1 ** 4 + 8 * (g - 1) * M1 ** 2 + 16))) / (4 * g * M1 ** 2)
            assert close(math.sin(GD.sonic_deflection(M1, g)["sigma_sonic"]) ** 2, s167, 1e-13) and close(math.sin(GD.max_deflection(M1, g)["sigma_max"]) ** 2, s168, 1e-13)
    for M1, tol in ((1.02, 0.05), (1.005, 0.015)):                              # (169), (170): near-sonic approximations, radians, gamma = 7/5
        k = (M1 ** 2 - 1) ** 1.5 / M1 ** 2
        assert rel(GD.max_deflection(M1)["delta_max"], forms["eq_170"]["coefficient_gamma_7_5"] * k) < tol
        assert rel(GD.sonic_deflection(M1)["delta_sonic"], forms["eq_169"]["coefficient_gamma_7_5"] * k) < tol
    assert abs(4 / (3 * math.sqrt(3) * 2.4) - forms["eq_170"]["coefficient_gamma_7_5"]) < 5e-5 and abs(1 / (math.sqrt(2) * 2.4) - forms["eq_169"]["coefficient_gamma_7_5"]) < 5e-5


def test_oblique_shock_normal_component_V2_derivation():  # V2 — D26 (two stars): tangential momentum → v continuous → a normal shock on M1 sin(sigma); C14
    r1, r2, u1, u2, v1, v2, g, p1, p2 = sp.symbols("rho_1 rho_2 u_1 u_2 v_1 v_2 gamma p_1 p_2", positive=True)
    sol = sp.solve([sp.Eq(r1 * u1, r2 * u2), sp.Eq(r1 * u1 * v1, r2 * u2 * v2)], [r2, v2], dict=True)[0]      # mass and tangential momentum
    assert sol[v2] == v1
    en = lambda p, r, u, v: g / (g - 1) * p / r + (u ** 2 + v ** 2) / 2          # noqa: E731
    assert z0((en(p1, r1, u1, v1) - en(p2, r2, u2, v2)).subs(v2, v1) - (g / (g - 1) * p1 / r1 + u1 ** 2 / 2 - g / (g - 1) * p2 / r2 - u2 ** 2 / 2))   # v^2/2 cancels from the energy law
    e = ch15.tangential_momentum_sympy()
    assert e["residual"] == 0 and str(e["v2"]) == "v_1"
    sig, dl, V1, V2 = sp.symbols("sigma delta V_1 V_2", positive=True)           # geometry of Fig. 15.18: (u, v) = V (sin, cos) on each side
    u1g, vg, u2g, v2g = V1 * sp.sin(sig), V1 * sp.cos(sig), V2 * sp.sin(sig - dl), V2 * sp.cos(sig - dl)
    V2s = sp.solve(sp.Eq(vg, v2g), V2)[0]
    assert z0(u1g / u2g.subs(V2, V2s) - sp.tan(sig) / sp.tan(sig - dl))          # rho2/rho1 = u1/u2 = tan(sigma)/tan(sigma − delta): second member of (15.51)
    M1 = sp.Symbol("M_1", positive=True)
    c1 = V1 / M1
    assert z0(u1g / c1 - M1 * sp.sin(sig))                                       # the normal Mach number is M1 sin(sigma): (15.39)–(15.41) apply with it


def test_deflection_relation_V2_derivation():  # V2 — D27 (two stars): tan(sigma)/tan(sigma − delta) = rho2/rho1 solved for tan(delta) → (15.53); C14
    g, M, t, r, td = sp.symbols("gamma M_1 t r t_d", positive=True)              # t = tan(sigma), t_d = tan(delta)
    sol = sp.solve(sp.Eq((t - td) / (1 + t * td), t / r), td)                    # tangent subtraction formula; tan(sigma − delta) = tan(sigma)/r
    assert len(sol) == 1 and z0(sol[0] - t * (r - 1) / (r + t ** 2))             # tan(delta) is collected before solving
    s2 = t ** 2 / (1 + t ** 2)                                                   # sin^2(sigma)
    rr = (g + 1) * M ** 2 * s2 / ((g - 1) * M ** 2 * s2 + 2)                     # (15.51), first member
    cos2 = (1 - t ** 2) / (1 + t ** 2)
    printed = 2 / t * (M ** 2 * s2 - 1) / (M ** 2 * (g + cos2) + 2)              # (15.53)
    assert z0(sol[0].subs(r, rr) - printed)
    assert z0(printed.subs(t, 1 / sp.sqrt(M ** 2 - 1))) and sp.limit(printed, t, sp.oo) == 0      # delta = 0 at sigma = mu (tan mu = 1/sqrt(M^2 − 1)) and at sigma = pi/2
    for g_ in GAMMAS:
        for M1 in (1.2, 2.0, 6.0):
            sig = np.linspace(math.asin(1 / M1), 0.5 * PI - 1e-6, 15)
            rho = GD.normal_shock_density_ratio(M1 * np.sin(sig), g_)
            assert close(np.tan(sig - GD.deflection_angle(M1, sig, g_)) * rho, np.tan(sig), 1e-9, 1e-12)
            assert abs(GD.deflection_angle(M1, math.asin(1 / M1), g_)) < 1e-12 and abs(GD.deflection_angle(M1, 0.5 * PI, g_)) < 1e-12 and GD.deflection_angle(M1, 0.5 * math.asin(1 / M1), g_) < 0


def test_weak_wave_rule_V2_derivation():  # V2 — D28 (two stars): (15.50) in (15.53) → (15.54); the limit sigma → mu → (15.55); C14
    g, M = sp.symbols("gamma M_1", positive=True)
    s = sp.Symbol("sigma", positive=True)
    Delta = 2 * g / (g + 1) * (M ** 2 * sp.sin(s) ** 2 - 1)                      # (p2 − p1)/p1 from (15.50)
    t53 = 2 * sp.cot(s) * (M ** 2 * sp.sin(s) ** 2 - 1) / (M ** 2 * (g + sp.cos(2 * s)) + 2)
    t54 = 2 * sp.cot(s) * (g + 1) / (2 * g) * Delta / (M ** 2 * (g + 1 - 2 * sp.sin(s) ** 2) + 2)
    assert z0(t53 - t54)                                                         # (15.54) is (15.53) with cos 2 sigma = 1 − 2 sin^2 sigma
    # the limit, without the three replacements of the text: d Delta/d delta at sigma = mu by the chain rule
    dD, dd = sp.diff(Delta, s), sp.diff(sp.atan(t53), s)
    for gv in (sp.Rational(7, 5), sp.Rational(5, 3)):
        for mv in (sp.Rational(3, 2), 2, 5):
            mu = sp.asin(1 / sp.Integer(1) / mv)
            slope = sp.N((dD / dd).subs({g: gv, M: mv, s: mu}), 30)
            assert abs(slope - sp.N(gv * mv ** 2 / sp.sqrt(mv ** 2 - 1), 30)) < sp.Float("1e-20")       # (15.55)
    # the text's route: tan(delta) ≈ delta, cot(sigma) → sqrt(M^2 − 1), sin^2(sigma) → 1/M^2
    D = sp.Symbol("D", positive=True)
    delta_lin = 2 * sp.sqrt(M ** 2 - 1) * (g + 1) / (2 * g) * D / (M ** 2 * (g + 1 - 2 / M ** 2) + 2)
    assert z0(D / delta_lin - g * M ** 2 / sp.sqrt(M ** 2 - 1))
    e = ch15.weak_wave_sympy()
    assert e["residual_15_54"] == 0 and e["residual_15_55"] == 0
    for gv in GAMMAS:
        assert close(GD.weak_wave_pressure(M_SUP, 0.013, gv), gv * M_SUP ** 2 * 0.013 / np.sqrt(M_SUP ** 2 - 1), 1e-14)


def test_oblique_shock_V1_two_roots_maximum_deflection_and_detachment():  # V1 + V7 — C14, N96–N99, trap T13
    for g in (1.4, 5.0 / 3.0):
        for M1 in (1.2, 2.0, I["M_wedge"], 4.0, 8.0):
            mu = math.asin(1 / M1)
            md, sd = GD.max_deflection(M1, g), GD.sonic_deflection(M1, g)
            opt = minimize_scalar(lambda s: -GD.deflection_angle(M1, s, g), bounds=(mu, 0.5 * PI), method="bounded", options=dict(xatol=1e-13))
            assert rel(md["delta_max"], -opt.fun) < 1e-10 and abs(md["sigma_max"] - opt.x) < 1e-5            # the closed form is the maximum of (15.53)
            assert abs(GD.oblique_shock_components(M1, sd["sigma_sonic"], g)["M2"] - 1) < 1e-11 and sd["sigma_sonic"] < md["sigma_max"] and sd["delta_sonic"] < md["delta_max"]
            for w in (0.05, 0.5, 0.9, 0.9999):
                dl = w * md["delta_max"]
                wk, st = GD.oblique_shock(M1, delta=dl, branch="weak", gamma=g), GD.oblique_shock(M1, delta=dl, branch="strong", gamma=g)
                assert abs(GD.deflection_angle(M1, wk["sigma"], g) - dl) < 1e-10 and abs(GD.deflection_angle(M1, st["sigma"], g) - dl) < 1e-10
                assert mu < wk["sigma"] < md["sigma_max"] < st["sigma"] < 0.5 * PI and st["M2"] < 1 and st["p_ratio"] > wk["p_ratio"] and wk["branch"] == "weak"
                assert (wk["M2"] > 1) == (dl < sd["delta_sonic"])                # a weak root between the sonic and the maximum deflection is subsonic behind
                again = GD.oblique_shock(M1, sigma=wk["sigma"], gamma=g)
                assert close(again["delta"], dl, 1e-9) and close(again["M2"], wk["M2"], 1e-12)
            w0, s0 = GD.oblique_shock(M1, delta=0.0, gamma=g), GD.oblique_shock(M1, delta=0.0, branch="strong", gamma=g)
            assert close(w0["sigma"], mu, 1e-14) and abs(w0["p_ratio"] - 1) < 1e-12 and close(w0["M2"], M1, 1e-9)                # weak → a Mach wave
            assert close(s0["sigma"], 0.5 * PI, 1e-15) and close(s0["M2"], GD.normal_shock_mach(M1, g), 1e-12) and close(s0["p_ratio"], GD.normal_shock_pressure_ratio(M1, g), 1e-12)   # strong → a normal shock
            top = GD.oblique_shock(M1, delta=md["delta_max"], gamma=g)
            assert close(top["sigma"], md["sigma_max"], 1e-12) and close(GD.oblique_shock(M1, delta=md["delta_max"], branch="strong", gamma=g)["sigma"], md["sigma_max"], 1e-12)
            with pytest.raises(ValueError, match="detached"):
                GD.oblique_shock(M1, delta=1.001 * md["delta_max"], gamma=g)
            arr = GD.oblique_shock(M1, delta=np.array([0.5, 1.5]) * md["delta_max"], gamma=g)
            assert list(arr["attached"]) == [True, False] and math.isnan(arr["sigma"][1]) and math.isnan(arr["p_ratio"][1]) and math.isnan(arr["delta"][1])
            assert rel(GD.detachment_mach(md["delta_max"], g), M1) < 1e-8          # the smallest Mach number that keeps this wedge's shock attached
            cur = GD.oblique_shock_curve(M1, g, n=801)
            assert abs(cur["delta"][0]) < 1e-15 and abs(cur["delta"][-1]) < 1e-12 and cur["M2"][0] == M1 and close(cur["mu"], mu, 1e-15) and cur["sigma"].shape == (801,)
            assert rel(cur["delta"].max(), md["delta_max"]) < 1e-5 and abs(cur["sigma"][cur["i_max"]] - md["sigma_max"]) < 2e-3 and abs(cur["M2"][cur["i_sonic"]] - 1) < 5e-3
            assert cur["delta_max"] == md["delta_max"] and cur["sigma_sonic"] == sd["sigma_sonic"] and cur["delta_sonic"] == sd["delta_sonic"] and cur["sigma_max"] == md["sigma_max"]
        assert rel(GD.max_deflection(1e7, g)["delta_max"], math.asin(1 / g)) < 1e-9 and GD.max_deflection(1.0, g)["delta_max"] == 0.0     # limits
        assert math.isinf(GD.detachment_mach(math.asin(1 / g) * 1.0001, g)) and GD.detachment_mach(0.0, g) == 1.0
    w = ch15.wedge_flow(I["M_wedge"], I["wedge_half_angle"])                     # the wedge of our inputs: attached, weak root
    o = GD.oblique_shock(I["M_wedge"], delta=I["wedge_half_angle"])
    assert w["attached"] and not w["subsonic_behind"] and all(close(w[k], o[k], 1e-14) for k in ("sigma", "M2", "p_ratio", "T_ratio", "rho_ratio", "p0_ratio"))
    assert close(w["delta_max"], GD.max_deflection(I["M_wedge"])["delta_max"], 0) and close(w["M_detach"], GD.detachment_mach(I["wedge_half_angle"]), 0) and w["M_detach"] < I["M_wedge"]
    d = ch15.wedge_flow(1.3, I["wedge_half_angle"])                              # the same wedge at a Mach number below M_detach: detached
    assert not d["attached"] and math.isnan(d["sigma"]) and d["subsonic_behind"] and close(d["axis_p_ratio"], GD.normal_shock_pressure_ratio(1.3), 0) and d["M_detach"] > 1.3
    for kw in (dict(), dict(delta=0.1, sigma=0.87), dict(delta=-0.05), dict(delta=0.1, branch="middle")):
        with pytest.raises(ValueError):
            GD.oblique_shock(2.0, **kw)
    with pytest.raises(ValueError):
        GD.oblique_shock(0.9, delta=0.05)


def test_oblique_shock_V7_shock_angle_below_the_mach_angle_is_refused():  # V7 — module convention "outside its range a function returns NaN or raises"; trap T6 for oblique shocks
    M1, mu = 2.0, math.asin(0.5)
    o = None
    try:
        o = GD.oblique_shock(M1, sigma=0.75 * mu)                               # M1 sin(sigma) < 1: an "expansion shock" (entropy would fall)
    except ValueError:
        return                                                                  # refused: fine
    bad = (o["dS_over_cv"], o["p_ratio"], o["attached"])
    assert math.isnan(o["p_ratio"]) or o["attached"] is False, f"oblique_shock(sigma < mu) silently returned an entropy-decreasing state: dS/Cv, p2/p1, attached = {bad}"


OBLIQUE_NUMERIC = ("sigma", "delta", "M2", "Mn1", "Mn2", "p_ratio", "rho_ratio", "T_ratio", "p0_ratio", "dS_over_cv")


def test_oblique_shock_V7_given_shock_angle_outside_its_range_scalars_raise_arrays_flag():  # V7 — the documented guard (fix of F2): M1 sin(sigma) >= 1 and sigma <= pi/2, to 1e-12; trap T6 for oblique shocks
    for g in (1.4, 5.0 / 3.0):
        for M1 in (1.2, 2.0, 5.0):
            mu = math.asin(1 / M1)
            for sig in (0.75 * mu, mu * (1 - 1e-9), 0.0, -0.3, -mu, 0.5 * PI + 1e-9, 2.0, PI - 0.5 * mu, math.asin((1 - 2e-12) / M1), 0.5 * PI + 2e-12):
                with pytest.raises(ValueError, match="not a shock"):         # scalars: refused, with the reason
                    GD.oblique_shock(M1, sigma=sig, gamma=g)
                assert all(math.isnan(v) for v in GD.oblique_shock_components(M1, sig, g).values()), (M1, sig)      # the components: NaN in every key
            for sig in (mu * (1 - 1e-13), math.asin((1 - 0.5e-12) / M1)):   # just below the Mach angle, within the allowance: the Mach wave, no jump at all
                assert 1 - 1e-12 <= M1 * math.sin(sig) < 1
                o, c = GD.oblique_shock(M1, sigma=sig, gamma=g), GD.oblique_shock_components(M1, sig, g)
                assert o["Mn1"] == 1.0 and o["delta"] == 0.0 and o["p_ratio"] == 1.0 and o["rho_ratio"] == 1.0 and o["dS_over_cv"] == 0.0 and o["p0_ratio"] == 1.0 and o["attached"] is True
                assert close(o["M2"], M1, 1e-9) and o["sigma"] == sig and c["Mn1"] == 1.0 and c["delta"] == 0.0 and close(c["V2_over_V1"], 1.0, 1e-12)
            o = GD.oblique_shock(M1, sigma=mu, gamma=g)                      # sigma = asin(1/M1) as a double: M1 sin(sigma) is 1 within an ulp on either side
            assert o["attached"] is True and abs(o["Mn1"] - 1) < 3e-16 and abs(o["delta"]) < 1e-15 and 1.0 <= o["p_ratio"] < 1 + 1e-15 and 0.0 <= o["dS_over_cv"] < 1e-40 and close(o["M2"], M1, 1e-9)
            for sig in (0.5 * PI, 0.5 * PI + 0.5e-12):                       # sigma = pi/2 to round-off: the normal shock
                o = GD.oblique_shock(M1, sigma=sig, gamma=g)
                assert o["attached"] is True and abs(o["delta"]) < 1e-11 and close(o["p_ratio"], GD.normal_shock_pressure_ratio(M1, g), 1e-12) and close(o["M2"], GD.normal_shock_mach(M1, g), 1e-12)
            # arrays: NaN in every numerical key and attached=False exactly at the refused entries; nothing returned lowers the entropy
            sg = np.linspace(-0.5, 3.6, 4101)
            o, c = GD.oblique_shock(M1, sigma=sg, gamma=g), GD.oblique_shock_components(M1, sg, g)
            inside = (M1 * np.sin(sg) >= 1 - 1e-12) & (sg <= 0.5 * PI + 1e-12)
            assert 100 < inside.sum() < sg.size and np.array_equal(o["attached"], inside) and o["attached"].dtype == bool and o["branch"] == "given"
            for k in OBLIQUE_NUMERIC:
                assert o[k].shape == sg.shape and np.array_equal(np.isfinite(o[k]), inside) and np.all(np.isnan(o[k][~inside])), k
            for k, v in c.items():
                assert v.shape == sg.shape and np.array_equal(np.isfinite(v), inside), k
            assert np.all(o["Mn1"][inside] >= 1) and np.all(o["p_ratio"][inside] >= 1) and np.all(o["dS_over_cv"][inside] >= 0) and np.all(o["p0_ratio"][inside] <= 1) and np.all(o["delta"][inside] > -1e-12)
            i = int(np.flatnonzero(inside)[inside.sum() // 2])              # an entry inside equals the scalar call
            one = GD.oblique_shock(M1, sigma=float(sg[i]), gamma=g)
            assert all(close(o[k][i], one[k], 1e-15, 1e-300) for k in OBLIQUE_NUMERIC)
    # the upstream Mach number: subsonic is refused; an array flags it; NaN propagates quietly (no exception, attached=False)
    for M1 in (0.9, 0.3):
        with pytest.raises(ValueError, match="not a shock"):
            GD.oblique_shock(M1, sigma=1.2)
    am = GD.oblique_shock(np.array([0.8, 2.0, 3.0]), sigma=0.6)
    assert list(am["attached"]) == [False, True, True] and math.isnan(am["p_ratio"][0]) and math.isnan(am["sigma"][0]) and close(am["p_ratio"][1:], [GD.oblique_shock(2.0, sigma=0.6)["p_ratio"], GD.oblique_shock(3.0, sigma=0.6)["p_ratio"]], 1e-15)
    for kw in (dict(M1=math.nan, sigma=1.2), dict(M1=2.0, sigma=math.nan)):
        o = GD.oblique_shock(kw["M1"], sigma=kw["sigma"])
        assert o["attached"] is False and all(math.isnan(o[k]) for k in OBLIQUE_NUMERIC if k != "sigma")
    s1 = GD.oblique_shock(1.0, sigma=0.5 * PI)                               # M1 = 1: the only "shock" is the Mach wave at 90 degrees
    assert s1["p_ratio"] == 1.0 and s1["M2"] == 1.0 and s1["delta"] == 0.0 and s1["attached"] is True
    # sigma = asin(1/M1) evaluated in doubles lands one ulp below M1 sin(sigma) = 1 for many M1: those must be Mach waves, not refusals
    Ms = np.linspace(1.01, 9.0, 800)
    low = Ms[Ms * np.sin(np.arcsin(1 / Ms)) < 1.0]
    assert low.size > 20
    for M1 in low[:25]:
        o = GD.oblique_shock(float(M1), sigma=math.asin(1 / M1))
        assert o["Mn1"] == 1.0 and o["p_ratio"] == 1.0 and o["dS_over_cv"] == 0.0 and o["attached"] is True
        w = GD.oblique_shock(float(M1), delta=0.0)                           # and the weak root of zero deflection is that same wave
        assert w["p_ratio"] == 1.0 and w["attached"] is True and close(w["M2"], M1, 1e-9)
    # the chart: no curve below M1 = 1; every point of a curve is a compression
    for bad in (0.9, 0.0, math.nan):
        with pytest.raises(ValueError, match="supersonic"):
            GD.oblique_shock_curve(bad)
    for M1 in (1.0, 1.05, 2.0, 7.0):
        cur = GD.oblique_shock_curve(M1, n=301)
        assert np.all(np.isfinite(cur["delta"])) and np.all(cur["delta"] >= 0) and np.all(np.isfinite(cur["M2"])) and cur["M2"][0] == M1 and cur["delta"][0] == 0.0
        assert np.all(M1 * np.sin(cur["sigma"]) >= 1 - 1e-12) and cur["sigma"][-1] == 0.5 * PI


def test_deflection_angle_V7_bare_formula_outside_the_range_and_its_guarded_callers():  # V7 — (15.53) is left unguarded on purpose (the root finder needs it); what it returns there, pinned
    for g in (1.4, 5.0 / 3.0):
        for M1 in (1.3, 2.0, 5.0):
            mu = math.asin(1 / M1)
            inside = np.linspace(mu, 0.5 * PI, 201)
            assert np.all(GD.deflection_angle(M1, inside, g) > -1e-15)       # in range: a compression turn
            below = np.linspace(0.05 * mu, 0.95 * mu, 50)
            d = GD.deflection_angle(M1, below, g)
            assert np.all(d < 0)                                             # 0 < sigma < mu: negative — the mark of the forbidden expansion root
            assert np.all(GD.deflection_angle(M1, np.linspace(0.5 * PI + 0.01, PI - mu - 0.01, 50), g) < 0)       # pi/2 < sigma < pi − mu: negative as well
            assert close(GD.deflection_angle(M1, -inside[1:-1], g), -GD.deflection_angle(M1, inside[1:-1], g), 1e-12)        # odd in sigma (the mirrored shock)
            mirrored = GD.deflection_angle(M1, PI - below, g)                # pi − mu < sigma < pi: M1 sin(sigma) < 1 and yet a POSITIVE angle — no mark here
            assert close(mirrored, -d, 1e-9) and np.all(mirrored > 0) and np.all(M1 * np.sin(PI - below) < 1)
            for s in (below, PI - below, -below):                            # which is why the guarded callers decide, not the sign of delta
                o = GD.oblique_shock(M1, sigma=s, gamma=g)
                assert not o["attached"].any() and np.all(np.isnan(o["delta"])) and np.all(np.isnan(GD.oblique_shock_components(M1, s, g)["delta"]))
    assert GD.deflection_angle(0.8, 0.9) < 0                                 # a subsonic stream: no shock at any angle


def test_weak_wave_V1_exact_shock_and_expansion_approach_the_linear_rule_from_both_sides():  # V1 + V5 (form) — (15.55), N100–N101, trap T16
    for g in (1.4, 5.0 / 3.0):
        for M in (1.5, 2.0, 3.0):
            b2 = M * M - 1
            c1, c2 = g * M * M / math.sqrt(b2), g * M * M * ((g + 1) * M ** 4 - 4 * b2) / (4 * b2 ** 2)       # NACA's expansion series: first and second coefficients
            d = np.radians([0.25, 0.5, 1.0, 2.0])
            sh = np.array([GD.weak_wave_error(M, float(x), g)["ratio"] for x in d])
            ex = np.array([GD.weak_wave_error(M, -float(x), g)["ratio"] for x in d])
            assert np.all(sh > 1) and np.all(ex < 1)                             # the shock is stronger, the expansion weaker, than the linear rule
            assert abs(observed_order(d, sh - 1) - 1.0) < 0.1 and abs(observed_order(d, 1 - ex) - 1.0) < 0.1        # the miss is first order in the angle
            assert rel((sh[0] - 1) / d[0], c2 / c1) < 0.02 and rel((1 - ex[0]) / d[0], c2 / c1) < 0.02               # with the second-order (Busemann) coefficient on both sides
            turn = np.radians([0.5, 1.0, 2.0])
            pm = GD.prandtl_meyer_turn(M, turn, g)["p_ratio"]
            rem = np.abs(pm - (1 - c1 * turn + c2 * turn ** 2))
            assert abs(observed_order(turn, rem) - 3.0) < 0.15                   # the cited series reproduces the exact expansion to second order
            assert close(GD.weak_wave_pressure(M, -0.01, g), -GD.weak_wave_pressure(M, 0.01, g), 1e-15) and close(ch15.ackeret_pressure(0.01, M, g), GD.weak_wave_pressure(M, 0.01, g), 1e-15)
    e = GD.weak_wave_error(2.0, math.radians(5.0))
    assert e["kind"] == "shock" and close(e["ratio"], e["exact"] / e["linear"], 1e-14) and close(e["rel_error"], e["linear"] / e["exact"] - 1, 1e-13) and abs(e["ratio"] - 1.118) < 5e-4
    assert GD.weak_wave_error(2.0, -0.05)["kind"] == "expansion" and GD.weak_wave_error(2.0, 0.0) == dict(linear=0.0, exact=0.0, ratio=1.0, rel_error=0.0, kind="none")
    assert math.isnan(GD.weak_wave_error(2.0, math.radians(30.0))["exact"]) and math.isnan(GD.weak_wave_pressure(1.0, 0.01)) and math.isnan(GD.weak_wave_pressure(0.7, 0.01))


def test_prandtl_meyer_V1_closed_form_equals_the_quadrature_of_its_integrand():  # V1 — (15.56) integral form against the closed form; C15
    for g in GAMMAS:
        for M in (1.0001, 1.3, 2.0, 4.5, 20.0):
            val, err = quad(lambda t: 2 * t * t / ((1 + t * t) * (2 + (g - 1) * (1 + t * t))), 0.0, math.sqrt(M * M - 1), epsabs=1e-14, epsrel=1e-13, limit=200)     # t = sqrt(M^2 − 1): smooth integrand
            assert abs(GD.prandtl_meyer(M, g) - val) < 1e-11 + 10 * err, (g, M)
            direct = quad(lambda m: math.sqrt(m * m - 1) / (1 + 0.5 * (g - 1) * m * m) / m, 1.0, M, epsabs=1e-13, epsrel=1e-12, limit=400)[0]                         # and in M itself
            assert abs(GD.prandtl_meyer(M, g) - direct) < 1e-9
        Mg = np.linspace(1.0, 30.0, 500)
        nu = GD.prandtl_meyer(Mg, g)
        assert nu[0] == 0.0 and np.all(np.diff(nu) > 0) and np.all(nu < GD.prandtl_meyer_max(g))
        assert close(GD.prandtl_meyer_max(g), 0.5 * PI * (math.sqrt((g + 1) / (g - 1)) - 1), 1e-15) and rel(GD.prandtl_meyer(1e9, g), GD.prandtl_meyer_max(g)) < 1e-8
        assert close(GD.mach_from_prandtl_meyer(GD.prandtl_meyer(M_SUP, g), g), M_SUP, 1e-10)
    assert math.isnan(GD.prandtl_meyer(0.9)) and math.isnan(GD.mach_from_prandtl_meyer(-0.01)) and math.isnan(GD.mach_from_prandtl_meyer(GD.prandtl_meyer_max())) and GD.mach_from_prandtl_meyer(0.0) == 1.0
    assert abs(math.degrees(GD.prandtl_meyer_max(1.4)) - 130.45) < 0.005 and abs(math.degrees(GD.prandtl_meyer(2.0)) - 26.38) < 0.005


def test_prandtl_meyer_differential_relation_V2_derivation():  # V2 — D29 (two stars): tangential speed kept across a weak wave + T0 constant → d delta; C15
    V, dV, mu, dd, M, g = sp.symbols("V dV mu ddelta M gamma", positive=True)
    eps = sp.Symbol("epsilon", positive=True)
    # a wave at the Mach angle mu to the stream; behind it the stream is turned by d delta towards the wave and has the speed V + dV
    tang = (V + eps * dV) * sp.cos(mu - eps * dd) - V * sp.cos(mu)
    first = sp.diff(tang, eps).subs(eps, 0)
    sol = sp.solve(first, dd)[0]
    assert z0(sol + dV / V * sp.cot(mu)) and z0(sp.cot(sp.asin(1 / M)) - sp.sqrt(M ** 2 - 1))        # d delta = −(dV/V) sqrt(M^2 − 1)
    c0 = sp.Symbol("c_0", positive=True)
    Tf = 1 + (g - 1) / 2 * M ** 2
    Vm = M * c0 / sp.sqrt(Tf)                                                    # V = M c, c = c0/sqrt(T0/T) at constant T0
    assert z0(sp.diff(sp.log(Vm), M) - 1 / (M * Tf))                             # dV/V = dM/M + dc/c = dM/(M (1 + (gamma−1)M^2/2))
    ddelta_dM = -sp.sqrt(M ** 2 - 1) * sp.diff(sp.log(Vm), M)
    assert z0(ddelta_dM + sp.sqrt(M ** 2 - 1) / M / Tf)                          # the unnumbered relation before (15.56)
    # the same step gives the weak-wave rule: dp = −rho V dV (15.26) and rho V^2 = gamma p M^2  →  dp/p = gamma M^2 d delta / sqrt(M^2 − 1)   (15.55)
    dp_over_p = -g * M ** 2 * (dV / V)
    assert z0(dp_over_p.subs(dV, sp.solve(sp.Eq(dd, -dV / V * sp.sqrt(M ** 2 - 1)), dV)[0]) - g * M ** 2 * dd / sp.sqrt(M ** 2 - 1))
    for gv in GAMMAS:                                                            # numerically: d nu/dM by central differences is the integrand
        Mv = np.array([1.2, 2.0, 5.0])
        fd = (GD.prandtl_meyer(Mv + 1e-6, gv) - GD.prandtl_meyer(Mv - 1e-6, gv)) / 2e-6
        assert close(fd, np.sqrt(Mv ** 2 - 1) / Mv / (1 + 0.5 * (gv - 1) * Mv ** 2), 1e-8)


def test_prandtl_meyer_closed_form_V2_derivation():  # V2 — D30 (three stars), every line of Part F: the integral → (15.56); C15
    g, t, M = sp.symbols("gamma t M", positive=True)
    lam = sp.sqrt((g - 1) / (g + 1))
    Mt = sp.sqrt(1 + t ** 2)                                                     # step 1: t = sqrt(M^2 − 1)
    assert z0(sp.diff(Mt, t) / Mt - t / (1 + t ** 2))                            # step 2: dM/M = t dt/(1 + t^2)
    assert z0(1 + (g - 1) / 2 * (1 + t ** 2) - (g + 1) / 2 * (1 + lam ** 2 * t ** 2))      # step 3: the bracket
    integrand_t = t / (1 + (g - 1) / 2 * (1 + t ** 2)) * t / (1 + t ** 2)        # step 4: sqrt(M^2−1)/(1 + (gamma−1)M^2/2) dM/M in t
    assert z0(integrand_t - 2 / (g + 1) * t ** 2 / ((1 + t ** 2) * (1 + lam ** 2 * t ** 2)))
    assert z0(t ** 2 / ((1 + t ** 2) * (1 + lam ** 2 * t ** 2)) - (1 / (1 + lam ** 2 * t ** 2) - 1 / (1 + t ** 2)) / (1 - lam ** 2))     # step 5: partial fractions
    assert z0(2 / (g + 1) / (1 - lam ** 2) - 1)                                  # step 6: the prefactor is exactly 1
    assert z0(integrand_t - (1 / (1 + lam ** 2 * t ** 2) - 1 / (1 + t ** 2)))
    F1, F2 = sp.atan(lam * t) / lam, sp.atan(t)                                  # steps 7–8: the two arctangent integrals
    assert z0(sp.diff(F1, t) - 1 / (1 + lam ** 2 * t ** 2)) and z0(sp.diff(F2, t) - 1 / (1 + t ** 2))
    nu = F1 - F2                                                                 # step 9: both vanish at t = 0
    assert nu.subs(t, 0) == 0
    closed = sp.sqrt((g + 1) / (g - 1)) * sp.atan(sp.sqrt((g - 1) / (g + 1) * (M ** 2 - 1))) - sp.atan(sp.sqrt(M ** 2 - 1))       # step 10: (15.56)
    for gv in (sp.Rational(7, 5), sp.Rational(5, 3)):
        for mv in (sp.Rational(11, 10), 2, 5):
            assert abs(sp.N((nu.subs(t, sp.sqrt(M ** 2 - 1)) - closed).subs({g: gv, M: mv}), 30)) < sp.Float("1e-25")
            assert abs(sp.N((sp.diff(closed, M) - sp.sqrt(M ** 2 - 1) / (1 + (g - 1) / 2 * M ** 2) / M).subs({g: gv, M: mv}), 30)) < sp.Float("1e-25")
    assert z0(sp.limit(nu.subs(g, sp.Rational(7, 5)), t, sp.oo) - (sp.sqrt(6) - 1) * sp.pi / 2)             # M → infinity: (sqrt((gamma+1)/(gamma−1)) − 1) pi/2
    num2 = sp.N(closed.subs({g: sp.Rational(7, 5), M: 2}) * 180 / sp.pi, 9)
    assert abs(num2 - sp.Float("26.379761")) < sp.Float("1e-5") and abs(sp.N((sp.sqrt(6) - 1) * 90, 9) - sp.Float("130.45408")) < sp.Float("1e-4")     # the numbers printed by the design's check cell
    assert abs(sp.N(sp.sqrt(6) * sp.atan(1 / sp.sqrt(2)) * 180 / sp.pi, 6) - sp.Float("86.38")) < sp.Float("0.01")       # Part F's worked number: 2.449 x 35.26 deg − 60 deg
    e = ch15.prandtl_meyer_sympy()                                               # the chapter's engine
    assert e["derivative_residual"] == 0 and e["residual_numeric"] < 1e-25 and e["value_at_1"] == 0
    for gv in GAMMAS:
        assert close(GD.prandtl_meyer(M_SUP, gv), [float(closed.subs({g: gv, M: float(m)})) for m in M_SUP], 1e-13)


def test_prandtl_meyer_turn_V4_isentropic_and_composes():  # V4 + V7 — N104, trap T14: delta + nu(M) = const; entropy and total enthalpy unchanged
    for g in (1.4, 1.25):
        for M1 in (1.0, 1.5, I["M_fan"], 4.0):
            for turn in np.radians([2.0, 15.0, 40.0]):
                s = GD.prandtl_meyer_turn(M1, float(turn), g)
                assert close(s["nu2"], s["nu1"] + turn, 1e-14) and close(GD.prandtl_meyer(s["M2"], g), s["nu2"], 1e-10) and s["M2"] > M1 and close(s["delta_book"], -turn, 0)
                assert close(s["p_ratio"], s["T_ratio"] ** (g / (g - 1)), 1e-13) and close(s["rho_ratio"], s["T_ratio"] ** (1 / (g - 1)), 1e-13)                    # entropy unchanged
                assert close(s["T_ratio"] * t0_ratio(s["M2"], g), t0_ratio(M1, g), 1e-13) and s["p_ratio"] < 1                                                  # total enthalpy unchanged
                assert close(s["mu1"], GD.mach_angle(M1), 0) and close(s["mu2"], GD.mach_angle(s["M2"]), 0)
                half = GD.prandtl_meyer_turn(GD.prandtl_meyer_turn(M1, 0.5 * turn, g)["M2"], 0.5 * turn, g)                                                     # two half-turns compose
                assert rel(half["M2"], s["M2"]) < 1e-9
                back = GD.prandtl_meyer_turn(s["M2"], -float(turn), g)       # an isentropic compression undoes it — also back to M1 = 1 exactly (nu = 0 ± round-off)
                assert abs(back["M2"] - M1) < 1e-8 and close(back["p_ratio"] * s["p_ratio"], 1.0, 1e-8)
    out = GD.prandtl_meyer_turn(2.0, math.radians(120.0))
    assert math.isnan(out["M2"]) and math.isnan(out["nu2"]) and math.isnan(GD.prandtl_meyer_turn(2.0, -math.radians(30.0))["M2"])                               # outside [0, nu_max)
    arr = GD.prandtl_meyer_turn(2.0, np.radians([0.0, 10.0]))
    assert arr["M2"].shape == (2,) and close(arr["M2"][0], 2.0, 1e-12)


def test_prandtl_meyer_turn_V7_return_to_mach_one_and_the_round_off_allowance():  # V7 — nu(M1) + turn in (−1e-12, 0) is the sonic state, not NaN; nothing else is forgiven
    n = 0
    for g in GAMMAS:
        for turn in np.radians(np.linspace(0.1, 60.0, 300)):                # out from M = 1 and back by the same angle: never NaN
            out = GD.prandtl_meyer_turn(1.0, float(turn), g)
            back = GD.prandtl_meyer_turn(out["M2"], -float(turn), g)
            assert abs(back["M2"] - 1.0) < 1e-8 and close(back["p_ratio"] * out["p_ratio"], 1.0, 1e-8) and back["nu2"] >= 0.0, (g, turn, back["M2"])
            n += back["M2"] == 1.0
    assert n > 100                                                          # many of these land a few ulp below zero and are set to the sonic state exactly
    for g in (1.4, 5.0 / 3.0):
        for M1 in (1.5, 2.0, 4.0):
            nu1 = GD.prandtl_meyer(M1, g)
            for over in (0.0, 1e-15, 0.5e-12):                              # inside the allowance: M2 = 1 exactly, with the sonic ratios
                s = GD.prandtl_meyer_turn(M1, -(nu1 + over), g)
                assert s["M2"] == 1.0 and s["nu2"] == 0.0 and s["mu2"] == 0.5 * PI and close(s["T_ratio"], t0_ratio(M1, g) / t0_ratio(1.0, g), 1e-14) and close(s["p_ratio"], s["T_ratio"] ** (g / (g - 1)), 1e-13)
            for over in (2e-12, 1e-9, 1e-3):                                # beyond it: still no state
                s = GD.prandtl_meyer_turn(M1, -(nu1 + over), g)
                assert all(math.isnan(s[k]) for k in ("M2", "nu2", "mu2", "p_ratio", "T_ratio", "rho_ratio")), (M1, over, s)
            arr = GD.prandtl_meyer_turn(M1, np.array([-(nu1 + 1e-14), -(nu1 + 1e-9), 0.1]), g)
            assert arr["M2"][0] == 1.0 and math.isnan(arr["M2"][1]) and arr["M2"][2] > M1 and arr["nu2"][0] == 0.0 and math.isnan(arr["nu2"][1])
        tiny = GD.prandtl_meyer_turn(1.0, 1e-13, g)                          # the allowance is one-sided: a small positive nu is a real (slightly supersonic) state
        assert tiny["M2"] > 1.0 and tiny["nu2"] == 1e-13 and close(GD.prandtl_meyer(tiny["M2"], g), 1e-13, 1e-6)
        assert math.isnan(GD.prandtl_meyer_turn(2.0, GD.prandtl_meyer_max(g) - GD.prandtl_meyer(2.0, g) + 1e-9, g)["M2"])      # the other end is unchanged: nu >= nu_max has no state


def test_expansion_fan_V4_mass_flux_between_rays_and_streamline_direction():  # V4 + V1 + V3 (the staircase) — N104–N105, Fig. 15.23 (geometry ours)
    g, M1, turn = 1.4, I["M_fan"], I["fan_turn"]
    f = ch15.expansion_fan(M1, turn, g, n_rays=2001, r_start=0.7)
    s = GD.prandtl_meyer_turn(M1, turn, g)
    assert close(f["ray_angle"][0], GD.mach_angle(M1), 1e-14) and close(f["ray_angle"][-1], s["mu2"] - turn, 1e-12) and close(f["M2"], s["M2"], 1e-12) and close(f["p2_over_p1"], s["p_ratio"], 1e-12)
    assert f["wall_angle2"] == -turn and f["streamline_r"][0] == pytest.approx(0.7) and close(f["theta"], np.linspace(0, turn, 2001), 1e-15) and close(f["mu"], GD.mach_angle(f["M"]), 0)
    rho, c = f["T_over_T1"] ** (1 / (g - 1)), np.sqrt(f["T_over_T1"])            # the velocity normal to a ray is the local sound speed: r rho c is the mass flux below the streamline
    assert np.ptp(f["streamline_r"] * rho * c) < 1e-12 * 0.7 and close(f["p_over_p1"], f["T_over_T1"] ** (g / (g - 1)), 1e-13)
    x, y = f["streamline_x"], f["streamline_y"]
    slope = np.diff(y) / np.diff(x)                                              # the streamline is tangent to the local flow direction −theta (independent of the mass argument)
    mid = 0.5 * (f["theta"][1:] + f["theta"][:-1])
    assert np.max(np.abs(np.arctan(slope) + mid)) < 1e-6
    assert close(x, f["streamline_r"] * np.cos(f["ray_angle"]), 0) and np.all(np.diff(f["M"]) > 0) and np.all(np.diff(f["ray_angle"]) < 0)       # the rays fan out clockwise
    with pytest.raises(ValueError):
        ch15.expansion_fan(2.0, 0.0)
    orders = {}
    Ns = np.array([20, 40, 80, 160, 320])                                        # the staircase of weak waves converges at first order (explicit Euler)
    st = [ch15.fan_staircase(M1, turn, int(N), g) for N in Ns]
    for key in ("M_error", "p_error", "p_weak_error"):
        orders[key] = observed_order(1.0 / Ns, [abs(s_[key]) for s_ in st])
        assert abs(orders[key] - 1.0) < 0.15, (key, orders[key])
    one = st[0]
    assert one["M"].shape == (21,) and one["M"][0] == M1 and one["p_over_p1"][0] == 1.0 and one["p_over_p1_weak"][0] == 1.0 and close(one["theta"], np.linspace(0, turn, 21), 1e-15)
    assert close(one["M2_exact"], s["M2"], 0) and close(one["p2_exact"], s["p_ratio"], 0) and close(one["p_over_p1"], GD.stagnation_pressure_ratio(M1, g) / GD.stagnation_pressure_ratio(one["M"], g), 1e-14)
    with pytest.raises(ValueError):
        ch15.fan_staircase(1.0, 0.1)


def test_curved_wall_V1_mach_lines_first_cross_on_the_envelope():  # V1 + V7 — N102–N103, Fig. 15.22: compression lines converge, expansion lines diverge
    M1, th_max = 2.0, math.radians(8.0)
    ramp = lambda x: th_max * np.clip(np.asarray(x, float), 0.0, 1.0) ** 2 * (3 - 2 * np.clip(np.asarray(x, float), 0.0, 1.0))     # noqa: E731 — smooth turn into the flow
    x = np.linspace(0.0, 1.0, 1601)
    w = ch15.curved_wall_mach_lines(M1, ramp, x)
    assert w["compression"] and w["M"][0] == pytest.approx(M1) and np.all(np.diff(w["M"]) <= 1e-12) and close(w["theta"], ramp(x), 0) and close(w["line_angle"], w["theta"] + w["mu"], 1e-15)
    assert close(GD.prandtl_meyer(w["M"]), GD.prandtl_meyer(M1) - w["theta"], 1e-10, 1e-12)             # nu(M) = nu(M1) − theta_w along the wall
    assert rel(w["y_wall"][-1], quad(lambda s: math.tan(float(ramp(s))), 0.0, 1.0, epsabs=1e-13)[0]) < 1e-6
    phi = w["line_angle"]                                                        # envelope of the family of lines: distance a = sin(mu)/(cos(theta_w) d phi/dx) along each line
    dphi = np.gradient(phi, x)
    a = np.where(dphi > 1e-9, np.sin(w["mu"]) / (np.cos(w["theta"]) * np.where(dphi > 1e-9, dphi, 1.0)), np.inf)
    i = int(np.argmin(a))
    env = np.array([x[i] + a[i] * math.cos(phi[i]), w["y_wall"][i] + a[i] * math.sin(phi[i])])
    assert np.linalg.norm(np.array(w["first_crossing"]) - env) < 0.01 * a[i] and w["crossings"].shape[1] == 2 and len(w["crossings"]) > 100
    e = ch15.curved_wall_mach_lines(M1, lambda s: -ramp(s), x)                   # the same wall turning away: an expansion
    assert not e["compression"] and e["first_crossing"] is None and e["crossings"].shape == (0, 2) and np.all(np.diff(e["M"]) >= -1e-12) and np.all(np.diff(e["line_angle"]) <= 1e-12)


# =====================================================================================================================
# 8. §15.10 — C16 the flat plate in supersonic flow (15.57)–(15.59); thin sections of any shape (ours)
# =====================================================================================================================
def test_flat_plate_V1_linear_theory_formulas():  # V1 — D31, (15.57)–(15.59); trap T15
    for g in (1.4, 5.0 / 3.0):
        for M in (1.3, I["M_inf"], 3.5):
            beta = math.sqrt(M * M - 1)
            for al in np.radians([0.5, 2.0, I["alpha"] * 180 / PI, 9.0]):
                a = ch15.ackeret_flat_plate(float(al), M)
                assert close([a["CL"], a["CD"], a["L_over_D"], a["lift_slope"]], [4 * al / beta, 4 * al ** 2 / beta, 1 / al, 4 / beta], 1e-13)            # (15.59)
                p = ch15.flat_plate_supersonic(float(al), M, g, method="linear")
                assert close([p["p_upper"], p["p_lower"]], [1 - g * M * M * al / beta, 1 + g * M * M * al / beta], 1e-13) and p["p_upper"] < 1 < p["p_lower"]   # expansion above, compression below
                assert close([p["CL"], p["CD"], p["CN"]], [a["CL"], a["CD"], a["CL"]], 1e-12) and close(p["CD"] / p["CL"], al, 1e-12)                   # D = L alpha
                assert close([p["L_over_pb"], p["D_over_pb"]], [2 * al * g * M * M / beta, 2 * al ** 2 * g * M * M / beta], 1e-12) and p["method"] == "linear"   # (15.58)
                assert close(ch15.ackeret_pressure(float(al), M, g), g * M * M * al / beta, 1e-14) and close(ch15.ackeret_cp(float(al), M), ch15.ackeret_pressure(float(al), M, g) / (0.5 * g * M * M), 1e-13)
                n = ch15.flat_plate_supersonic(-float(al), M, g, method="linear")
                assert close(n["CL"], -p["CL"], 1e-13) and close(n["CD"], p["CD"], 1e-13) and close([n["p_upper"], n["p_lower"]], [p["p_lower"], p["p_upper"]], 1e-14)
            assert close(ch15.lift_slope_vs_mach(M), 4 / beta, 1e-14)
    assert all(math.isnan(v) for v in ch15.ackeret_flat_plate(0.05, 0.9).values()) and math.isnan(ch15.ackeret_cp(0.05, 1.0)) and math.isnan(ch15.ackeret_pressure(0.05, 0.6))
    assert math.isnan(ch15.flat_plate_supersonic(0.05, 0.9, method="linear")["CL"])
    assert close(ch15.lift_slope_vs_mach(np.array([0.0, 0.3, 0.79])), 2 * PI, 1e-15) and math.isnan(ch15.lift_slope_vs_mach(1.0)) and math.isnan(ch15.lift_slope_vs_mach(0.8))
    assert close(ch15.lift_slope_vs_mach(0.6, model="prandtl_glauert"), 2 * PI / 0.8, 1e-14) and close(ch15.lift_slope_vs_mach(0.9, band=(0.95, 1.05)), 2 * PI, 1e-15)
    assert rel(ch15.zhukhovsky_lift_coefficient(1e-4, 0.0) / 1e-4, ch15.lift_slope_vs_mach(0.2)) < 1e-6        # the 2 pi of (14.12)
    with pytest.raises(ValueError):
        ch15.lift_slope_vs_mach(0.5, model="karman_tsien")
    with pytest.raises(ValueError):
        ch15.flat_plate_supersonic(0.05, 2.0, method="newtonian")


def test_flat_plate_V7_linear_theory_is_the_small_angle_limit_of_shock_expansion():  # V7 — C16: exact weak shock + exact fan → (15.59); the order of the limit measured (two closed forms, no discretisation)
    orders = {}
    for g in (1.4, 5.0 / 3.0):
        for M in (1.6, I["M_inf"], 3.0):
            al = np.radians([0.25, 0.5, 1.0, 2.0, 4.0])
            ex = [ch15.flat_plate_supersonic(float(a), M, g, method="shock_expansion") for a in al]
            li = [ch15.flat_plate_supersonic(float(a), M, g, method="linear") for a in al]
            eL = np.array([abs(e["CL"] / l["CL"] - 1) for e, l in zip(ex, li)])
            eD = np.array([abs(e["CD"] / l["CD"] - 1) for e, l in zip(ex, li)])
            oL, oD = observed_order(al, eL), observed_order(al, eD)
            orders[(g, M)] = (oL, oD)
            # MEASURED: second order. The second-order (Busemann) pressure terms are equal on the two surfaces and cancel in p3 − p2,
            # so the plate's lift and drag coefficients miss by O(alpha^2) — one order better than a single surface (and than the docstring says).
            assert abs(oL - 2.0) < 0.15 and abs(oD - 2.0) < 0.15, (g, M, oL, oD)
            assert eL[0] < 1e-4 and eD[0] < 1e-4
            b2 = M * M - 1
            c1, c2 = g * M * M / math.sqrt(b2), g * M * M * ((g + 1) * M ** 4 - 4 * b2) / (4 * b2 ** 2)
            one = np.array([(e["p_lower"] - 1) / (c1 * a) - 1 for e, a in zip(ex, al)])                    # a single surface is only first-order accurate
            assert abs(observed_order(al, np.abs(one)) - 1.0) < 0.1 and rel(one[0] / al[0], c2 / c1) < 0.02
            for e, a in zip(ex, al):
                assert close(e["CD"] / e["CL"], math.tan(a), 1e-12) and e["p_upper"] < 1 < e["p_lower"] and e["method"] == "shock_expansion"            # the force is normal to the plate
                assert close(e["p_lower"], GD.oblique_shock(M, delta=float(a), gamma=g)["p_ratio"], 0) and close(e["p_upper"], GD.prandtl_meyer_turn(M, float(a), g)["p_ratio"], 0)
                assert close(e["CN"], (e["p_lower"] - e["p_upper"]) / (0.5 * g * M * M), 1e-13) and close(e["L_over_pb"], e["CL"] * 0.5 * g * M * M, 1e-13)
    assert all(abs(o - 2.0) < 0.15 for pair in orders.values() for o in pair)
    z = ch15.flat_plate_supersonic(0.0, 2.0, method="shock_expansion")
    assert z["CL"] == 0.0 and z["CD"] == 0.0 and z["p_upper"] == 1.0
    with pytest.raises(ValueError, match="detached"):
        ch15.flat_plate_supersonic(math.radians(30.0), 2.0, method="shock_expansion")
    al = I["alpha"]
    w = ch15.flat_plate_wave_pattern(al, I["M_inf"])
    e = ch15.flat_plate_supersonic(al, I["M_inf"], method="shock_expansion")
    assert close([w["p2"], w["p3"]], [e["p_upper"], e["p_lower"]], 1e-14) and w["M2"] > I["M_inf"] > w["M3"] and w["plate_angle"] == -al
    assert close(w["le_fan_first"], GD.mach_angle(I["M_inf"]), 0) and close(w["le_fan_last"], GD.mach_angle(w["M2"]) - al, 1e-14) and close(w["le_shock"], -GD.oblique_shock(I["M_inf"], delta=al)["sigma"], 0)
    assert w["le_fan_first"] > w["le_fan_last"] > -al and w["le_shock"] < -al and w["te_shock"] > -al and w["te_fan_first"] < w["te_fan_last"] < 0
    assert abs(w["p4_upper"] - 1) < 5 * al ** 3 and abs(w["p4_lower"] - 1) < 5 * al ** 3 and abs(w["M4_upper"] - I["M_inf"]) < 5 * al ** 3       # both streams return (nearly) to the free stream
    tiny = ch15.flat_plate_wave_pattern(1e-5, I["M_inf"])                        # alpha → 0: every wave is a Mach line
    mu = GD.mach_angle(I["M_inf"])
    assert all(abs(abs(tiny[k]) - mu) < 1e-4 for k in ("le_fan_first", "le_fan_last", "le_shock", "te_shock", "te_fan_first", "te_fan_last"))
    with pytest.raises(ValueError):
        ch15.flat_plate_wave_pattern(0.0, 2.0)


def test_ackeret_airfoil_V1_thickness_camber_and_incidence_drag():  # V1 + V3 — N112 (ours, Exercise-type): C_L independent of shape; C_D = (4/beta)[alpha^2 + <yc'^2> + <yt'^2>]
    M, beta = I["M_inf"], math.sqrt(I["M_inf"] ** 2 - 1)
    tau, al, hc = I["thickness_ratio"], I["alpha"], 0.02
    zero = lambda x: 0.0 * np.asarray(x, float)                                  # noqa: E731
    flat = ch15.ackeret_airfoil(zero, zero, al, M)
    assert close([flat["CL"], flat["CD"], flat["CD_incidence"]], [4 * al / beta, 4 * al ** 2 / beta, 4 * al ** 2 / beta], 1e-12) and abs(flat["CD_thickness"]) < 1e-15 and abs(flat["CD_camber"]) < 1e-15
    assert close(flat["cp_lower"] - flat["cp_upper"], 4 * al / beta, 1e-12) and flat["x"].shape == (2001,)
    yt = lambda x: 0.5 * tau * np.sin(PI * np.asarray(x, float))                 # noqa: E731 — half-sine half-thickness
    sine = ch15.ackeret_airfoil(yt, lambda x: -yt(x), al, M)
    assert rel(sine["CD_thickness"], PI ** 2 * tau ** 2 / (2 * beta)) < 1e-5 and rel(sine["CL"], 4 * al / beta) < 1e-10 and abs(sine["CD_camber"]) < 1e-15       # the coefficient pi^2/2
    wedge = lambda x: tau * np.minimum(np.asarray(x, float), 1 - np.asarray(x, float))      # noqa: E731 — double wedge of thickness tau: slopes ±tau
    dw = ch15.ackeret_airfoil(wedge, lambda x: -wedge(x), 0.0, M)
    assert rel(dw["CD"], 4 * tau ** 2 / beta) < 2e-3 and abs(dw["CL"]) < 1e-12   # the classical double-wedge wave drag 4 tau^2/beta
    yc = lambda x: 4 * hc * np.asarray(x, float) * (1 - np.asarray(x, float))    # noqa: E731 — parabolic camber line
    cam = ch15.ackeret_airfoil(lambda x: yc(x) + yt(x), lambda x: yc(x) - yt(x), al, M)
    assert rel(cam["CL"], 4 * al / beta) < 1e-10 and rel(cam["CD_camber"], 4 / beta * 16 * hc ** 2 / 3) < 1e-5 and abs(cam["CD_cross"]) < 1e-12     # camber gives drag but no lift
    assert rel(cam["CD"], cam["CD_incidence"] + cam["CD_camber"] + cam["CD_thickness"] + cam["CD_cross"]) < 1e-12
    neg = ch15.ackeret_airfoil(yt, lambda x: -yt(x), -al, M)
    assert close(neg["CD"], sine["CD"], 1e-12) and close(neg["CL"], -sine["CL"], 1e-12)                # a symmetric section: C_D even, C_L odd in alpha
    errs = [abs(ch15.ackeret_airfoil(yt, lambda x: -yt(x), 0.0, M, n=n)["CD"] / (PI ** 2 * tau ** 2 / (2 * beta)) - 1) for n in (51, 101, 201, 401)]
    assert abs(observed_order(1 / np.array([50.0, 100.0, 200.0, 400.0]), errs) - 2.0) < 0.15            # V3: slopes by second-order differences


# =====================================================================================================================
# 9. Optional numerics (ours): the exact Riemann solver and the one-dimensional Euler run (N114)
# =====================================================================================================================
RIEMANN = {"sod": ((1.0, 0.0, 1.0), (0.125, 0.0, 0.1)), "strong": ((1.0, 0.0, 300.0), (1.0, 0.0, 0.3)), "collide": ((1.0, 1.7, 1.0), (1.0, -1.7, 1.0)),
           "moving": ((2.0, 0.4, 3.0), (1.1, -0.2, 0.9))}


def _rh(rho_a, u_a, p_a, rho_b, u_b, p_b, S, g):
    fa = np.array([rho_a * (u_a - S), p_a + rho_a * (u_a - S) ** 2, g / (g - 1) * p_a / rho_a + 0.5 * (u_a - S) ** 2])
    fb = np.array([rho_b * (u_b - S), p_b + rho_b * (u_b - S) ** 2, g / (g - 1) * p_b / rho_b + 0.5 * (u_b - S) ** 2])
    return np.max(np.abs(fb / fa - 1))


def test_riemann_exact_V4_shocks_satisfy_the_jump_laws_and_fans_are_isentropic():  # V4 — GD.riemann_exact against (15.35)–(15.37) in the shock frame
    for g in (1.4, 5.0 / 3.0):
        for name, (L, R) in RIEMANN.items():
            s = GD.riemann_exact(L, R, g, 0.0)
            ps, us = s["p_star"], s["u_star"]
            for side, (r0, u0, p0), rs, S, wave, sgn in (("L", L, s["rho_star_left"], s["shock_speed_left"], s["left_wave"], -1), ("R", R, s["rho_star_right"], s["shock_speed_right"], s["right_wave"], +1)):
                c0 = math.sqrt(g * p0 / r0)
                if wave == "shock":
                    assert ps > p0 and _rh(r0, u0, p0, rs, us, ps, S, g) < 1e-10, (name, side)
                    assert sgn * (S - u0) > c0                                    # the front outruns sound in the gas ahead
                else:
                    assert ps <= p0 and math.isnan(S) and rel(ps / rs ** g, p0 / r0 ** g) < 1e-12                        # isentropic
                    cs = math.sqrt(g * ps / rs)
                    assert abs((us - sgn * 2 * cs / (g - 1)) - (u0 - sgn * 2 * c0 / (g - 1))) < 1e-10 * c0             # the Riemann invariant through the fan
            a = 1.2 * max(abs(L[1]) + math.sqrt(g * L[2] / L[0]), abs(R[1]) + math.sqrt(g * R[2] / R[0]), abs(np.nan_to_num(s["shock_speed_left"])), abs(np.nan_to_num(s["shock_speed_right"])))
            # global balance: d/dt of the integral over [−a t, a t] equals flux in − flux out; integrate piece by piece between the wave edges (Gauss–Legendre)
            cL, cR = math.sqrt(g * L[2] / L[0]), math.sqrt(g * R[2] / R[0])
            left_edges = [s["shock_speed_left"]] if s["left_wave"] == "shock" else [L[1] - cL, us - math.sqrt(g * ps / s["rho_star_left"])]
            right_edges = [s["shock_speed_right"]] if s["right_wave"] == "shock" else [us + math.sqrt(g * ps / s["rho_star_right"]), R[1] + cR]
            edges = [-a] + left_edges + [us] + right_edges + [a]
            assert all(e2 > e1 for e1, e2 in zip(edges[:-1], edges[1:]))          # the waves are ordered: left wave, contact, right wave
            cons = lambda r, u, p: np.array([r, r * u, p / (g - 1) + 0.5 * r * u * u])                # noqa: E731
            flx = lambda r, u, p: np.array([r * u, r * u * u + p, u * (p / (g - 1) + 0.5 * r * u * u + p)])   # noqa: E731
            nodes, weights = np.polynomial.legendre.leggauss(48)
            integral = np.zeros(3)
            for e1, e2 in zip(edges[:-1], edges[1:]):
                xi = 0.5 * (e1 + e2) + 0.5 * (e2 - e1) * nodes
                w = GD.riemann_exact(L, R, g, xi)
                integral += 0.5 * (e2 - e1) * (cons(w["rho"], w["u"], w["p"]) @ weights)
                assert np.all(w["p"] > 0) and np.all(w["rho"] > 0)
            want = a * (cons(*L) + cons(*R)) + flx(*L) - flx(*R)
            scale = a * (np.abs(cons(*L)) + np.abs(cons(*R))) + np.abs(flx(*L)) + np.abs(flx(*R))
            assert np.max(np.abs(integral - want) / scale) < 1e-11, (name, integral, want)
            xi = np.linspace(-a, a, 101)
            w = GD.riemann_exact(L, R, g, xi)
            assert w["rho"][0] == L[0] and w["p"][-1] == R[2]
            m = GD.riemann_exact((R[0], -R[1], R[2]), (L[0], -L[1], L[2]), g, -xi)                           # mirror symmetry
            assert close(m["rho"], w["rho"], 1e-12) and close(m["u"], -w["u"], 1e-12, 1e-12) and close(m["p_star"], ps, 1e-12)
    sod = GD.riemann_exact(*RIEMANN["sod"], 1.4, 0.0)                              # the right shock of Sod's problem is a moving shock into still gas
    ms = GD.moving_shock(0.1, 0.1 / 0.125, p2_over_p1=sod["p_star"] / 0.1, gas=GD.PerfectGas(1.4, 1.0))
    assert rel(ms["u_shock"], sod["shock_speed_right"]) < 1e-12 and rel(ms["u_gas"], sod["u_star"]) < 1e-10 and rel(ms["rho2"], sod["rho_star_right"]) < 1e-12
    assert sod["left_wave"] == "rarefaction" and sod["right_wave"] == "shock"
    fan = GD.riemann_exact(*RIEMANN["sod"], 1.4, np.linspace(-1.1, -0.2, 7))       # inside the left fan u − c = x/t
    inside = (fan["rho"] < 1.0) & (fan["rho"] > sod["rho_star_left"])
    assert inside.sum() >= 4 and close((fan["u"] - np.sqrt(1.4 * fan["p"] / fan["rho"]))[inside], np.linspace(-1.1, -0.2, 7)[inside], 1e-12)
    tr = GD.riemann_exact((1.3, 0.3, 2.0), (1.3, 0.3, 2.0), 1.4, np.array([-5.0, 0.0, 5.0]))
    assert close(tr["rho"], 1.3, 1e-14) and close(tr["u"], 0.3, 1e-13) and close(tr["p_star"], 2.0, 1e-12)
    with pytest.raises(ValueError, match="vacuum"):
        GD.riemann_exact((1.0, -10.0, 1.0), (1.0, 10.0, 1.0))
    col = GD.riemann_exact(*RIEMANN["collide"], 1.4, 0.0)
    assert col["left_wave"] == col["right_wave"] == "shock" and abs(col["u_star"]) < 1e-12 and close(col["shock_speed_left"], -col["shock_speed_right"], 1e-12)


def _simple_wave_state(x, amp, c0, g, rho00):
    u0 = amp * np.sin(2 * PI * x)
    cc = c0 + 0.5 * (g - 1) * u0
    return rho00 * (cc / c0) ** (2 / (g - 1)), u0, rho00 * c0 ** 2 / g * (cc / c0) ** (2 * g / (g - 1))


def test_euler1d_V3_second_order_for_a_smooth_wave_and_first_order_at_a_captured_shock():  # V3 + V4 — ch15.euler1d_run (MacCormack + our smoothing; DEVIATION)
    g, c0, rho00, amp = 1.4, 340.0, 1.2, 5.0
    tb = AC.breaking_time(None, None, g, max_compression_rate=amp * 2 * PI)
    T = 0.5 * tb
    xf = np.arange(40000) / 40000
    sw = AC.simple_wave_evolve(amp * np.sin(2 * PI * xf), xf, T, c0, g)             # the exact simple wave before breaking
    xp = np.mod(sw["x_points"], 1.0)
    o = np.argsort(xp)
    hs, errs = [], []
    for n in (100, 200, 400, 800):
        x = np.arange(n) / n
        rho, u0, p = _simple_wave_state(x, amp, c0, g, rho00)
        r = ch15.euler1d_run(rho, u0, p, x, T, g, eps_av=0.0)
        errs.append(np.max(np.abs(r["u"] - np.interp(x, xp[o], sw["u"][o], period=1.0))))
        hs.append(1.0 / n)
        assert abs(r["mass"] / r["mass0"] - 1) < 1e-13 and abs(r["energy"] / r["energy0"] - 1) < 1e-13 and abs(r["momentum"] - r["momentum0"]) < 1e-11 * rho00 * c0       # V4: conservative to round-off
        assert r["t"] == pytest.approx(T) and r["n_steps"] > 0 and r["rho"].shape == (n,)
    po = pairwise_orders(hs, errs)
    assert abs(observed_order(hs[1:], errs[1:]) - 2.0) < 0.15 and abs(po[-1] - 2.0) < 0.05, (errs, po)
    x = np.arange(400) / 400                                                         # past breaking, with the smoothing on: a captured shock, still conservative
    rho, u0, p = _simple_wave_state(x, 30.0, c0, g, rho00)
    r = ch15.euler1d_run(rho, u0, p, x, 3.0 * AC.breaking_time(None, None, g, max_compression_rate=30.0 * 2 * PI), g)
    assert abs(r["mass"] / r["mass0"] - 1) < 1e-13 and abs(r["energy"] / r["energy0"] - 1) < 1e-13 and abs(r["momentum"] - r["momentum0"]) < 1e-10 * rho00 * c0
    assert np.max(-np.diff(r["u"])) > 5 * 30.0 * 2 * PI / 400 and np.all(np.isfinite(r["p"])) and np.all(r["p"] > 0)       # a jump a few cells thick has formed
    L, R = RIEMANN["sod"]
    hs, e_rho, e_u, e_p = [], [], [], []
    for n in (100, 200, 400, 800):
        x = (np.arange(n) + 0.5) / n
        r = ch15.euler1d_run(np.where(x < 0.5, L[0], R[0]), np.zeros(n), np.where(x < 0.5, L[2], R[2]), x, 0.2, g, bc="transmissive")
        ex = GD.riemann_exact(L, R, g, (x - 0.5) / 0.2)
        hs.append(1.0 / n)
        e_rho.append(np.mean(np.abs(r["rho"] - ex["rho"])))
        e_u.append(np.mean(np.abs(r["u"] - ex["u"])))
        e_p.append(np.mean(np.abs(r["p"] - ex["p"])))
    o_rho, o_u, o_p = observed_order(hs, e_rho), observed_order(hs, e_u), observed_order(hs, e_p)
    assert abs(o_p - 1.0) < 0.15 and abs(o_u - 1.0) < 0.15 and abs(o_rho - 1.0) < 0.25, (o_rho, o_u, o_p)      # L1 order 1 at the shock; the smeared contact holds the density back
    S = GD.riemann_exact(L, R, g, 0.0)["shock_speed_right"]
    i_sh = int(np.argmax(-np.diff(r["p"])))                                           # the steepest pressure drop on the finest grid sits at the exact shock position
    assert abs(x[i_sh] - (0.5 + S * 0.2)) < 4.0 / 800 and close(r["p"][0], L[2], 1e-12) and close(r["rho"][-1], R[0], 1e-12)
    xs = np.linspace(0, 1, 16)
    with pytest.raises(ValueError):
        ch15.euler1d_run(np.ones(16), 0.0, np.ones(16), xs ** 2, 0.1)                 # non-uniform grid
    with pytest.raises(ValueError):
        ch15.euler1d_run(np.ones(16), 0.0, np.ones(16), xs, 0.1, cfl=1.1)
    with pytest.raises(ValueError):
        ch15.euler1d_run(np.ones(16), 0.0, np.ones(16), xs, 0.1, bc="reflecting")
    rest = ch15.euler1d_run(1.2, 0.0, 1.0e5, xs, 1e-3)                                # a uniform state stays uniform
    assert close(rest["rho"], 1.2, 1e-14) and np.max(np.abs(rest["u"])) < 1e-10 and close(rest["p"], 1.0e5, 1e-13)


# =====================================================================================================================
# 10. Tables, sketches, scripts, reference data
# =====================================================================================================================
def test_tables_V7_inputs_conventions_slips_and_traps():  # V7 — I31: the chapter's tables
    assert ch15.illustrative_inputs() == I and len(I) == 43 and all(isinstance(v, (int, float)) for v in I.values())
    c = GD.laval_critical_pressures(I["Ae_over_At"], I["gamma"])
    assert close(I["Ae_over_At"], I["A_exit"] / I["A_throat"], 1e-14) and c["p_d"] < I["pB_over_p0_shock"] < c["p_b"]                 # the "shock in the nozzle" input really gives one
    assert I["wedge_half_angle"] < GD.max_deflection(I["M_wedge"])["delta_max"] and I["alpha"] < GD.max_deflection(I["M_inf"])["delta_max"]
    lim_sub, lim_sup = GD.duct_choking_limits(I["M1_duct_sub"]), GD.duct_choking_limits(I["M1_duct_sup"])
    assert I["f"] < min(lim_sub["f_max"], lim_sup["f_max"]) and I["q"] < min(lim_sub["q_max"], lim_sup["q_max"])                      # none of the duct inputs is choked
    assert GD.prandtl_meyer(I["M_fan"]) + I["fan_turn"] < GD.prandtl_meyer_max() and I["U_stream_sub"] < AIR.c(I["T_air"]) < I["U_stream_sup"] and I["M_s"] > 1
    assert AC.linearisation_error(I["pulse_amplitude"], I["p_air"] / (I["R"] * I["T_air"]), AIR.c(I["T_air"])) < 1e-4                 # the sound pulse is acoustic
    conv, slips, traps = ch15.conventions_table(), ch15.book_slips(), ch15.traps()
    assert conv.shape == (9, 4) and list(conv.columns) == ["item", "here", "elsewhere", "code"]
    assert slips.shape == (17, 7) and list(slips["n"]) == list(range(1, 18)) and list(slips.loc[slips["changes_physics"], "n"]) == [4, 11, 12]
    assert set(slips["kind"]) == {"slip", "loose"} and int((slips["kind"] == "slip").sum()) == 9
    for n, opt in ((4, "printed_15_32"), (11, "printed_15_47"), (12, "printed_15_48")):
        assert opt in slips.loc[slips["n"] == n, "code"].iloc[0]
    assert traps.shape == (22, 3) and list(traps["id"]) == [f"T{k}" for k in range(1, 23)]


def test_drawings_V7_every_sketch_returns_a_figure():  # V7 — design C.4: scripts/ch15_drawings.py (12 helpers)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure

    sys.path.insert(0, str(ROOT / "scripts"))
    try:
        import ch15_drawings as DR
    finally:
        sys.path.remove(str(ROOT / "scripts"))
    names = ("draw_duct_control_volume draw_reservoir_nozzle_tank draw_four_ducts draw_shock_two_frames draw_shock_triangles draw_wedge_streamlines "
             "draw_mach_triangle draw_plate_waves draw_bow_shock_sketch draw_jet_sketch draw_fan draw_wavefront_panels").split()
    assert set(names) == set(DR.__all__)
    wf = ch15.wedge_flow(I["M_wedge"], I["wedge_half_angle"])
    special = {"draw_shock_two_frames": lambda: DR.draw_shock_two_frames(GD.moving_shock(I["p_air"], I["T_air"], M_s=I["M_s"])),
               "draw_shock_triangles": lambda: DR.draw_shock_triangles(I["M_wedge"], wf["sigma"], I["wedge_half_angle"]),
               "draw_wedge_streamlines": lambda: DR.draw_wedge_streamlines(I["M_wedge"], I["wedge_half_angle"], wf["sigma"]),
               "draw_mach_triangle": lambda: DR.draw_mach_triangle(I["M_flight"]),
               "draw_plate_waves": lambda: DR.draw_plate_waves(ch15.flat_plate_wave_pattern(I["alpha"], I["M_inf"])),
               "draw_fan": lambda: DR.draw_fan(ch15.expansion_fan(I["M_fan"], I["fan_turn"])),
               "draw_wavefront_panels": lambda: DR.draw_wavefront_panels([AC.source_wavefronts(U, 340.0, 4.0, 4, 1.0) for U in (0.0, 170.0, 680.0)])}
    for name in names:
        for call in ([special[name]] if name in special else []) + [getattr(DR, name)]:
            fig = call()
            assert isinstance(fig, Figure) and len(fig.axes) >= 1 and any(len(ax.lines) + len(ax.patches) + len(ax.collections) + len(ax.texts) > 0 for ax in fig.axes), name
            plt.close(fig)
    for kind in ("overexpanded", "underexpanded"):
        fig = DR.draw_jet_sketch(kind)
        assert isinstance(fig, Figure)
        plt.close(fig)
    assert isinstance(DR.draw_wedge_streamlines(1.3, math.radians(20.0)), Figure)      # the detached case is drawn, not refused
    fig, ax = plt.subplots()
    assert DR.draw_mach_triangle(2.0, ax=ax) is fig                                   # drawing into a given Axes returns its Figure
    plt.close("all")


@slow
def test_scripts_V7_every_chapter_script_runs_headless(tmp_path):  # V7 — the demos + the drawings + make_refs exit 0
    ran = 0
    for s in sorted((ROOT / "scripts").glob("ch15_*.py")):
        if s.name == "ch15_common.py":
            continue
        r = subprocess.run([sys.executable, str(s), "--no-show", "--fast", "--out", str(tmp_path)], cwd=str(ROOT), capture_output=True,
                           text=True, encoding="utf-8", errors="replace", timeout=900)
        assert r.returncode == 0, (s.name, r.stderr[-2000:])
        ran += 1
    assert ran == 8 and len(list(tmp_path.glob("*.png"))) >= 8
    before = (REF / "benchmarks.json").read_bytes()
    r = subprocess.run([sys.executable, str(REF / "make_refs.py")], cwd=str(ROOT), capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=120)
    assert r.returncode == 0 and (REF / "benchmarks.json").read_bytes() == before      # idempotent


def test_reference_V5_sources_file_cites_every_entry_and_no_book_table():  # V5 — data-and-benchmarks policy
    src = (REF / "SOURCES.md").read_text(encoding="utf-8")
    bm = bench()
    keys = [k for k in bm if not k.startswith("_")]
    assert len(keys) == 7
    for key in keys:
        row = bm[key]
        assert key in src and row["verified"] == "2026-10-13" and row["url"].startswith("https://") and row["url"] in src, key
        assert "kundu" not in row["source"].lower() and "textbook" not in row["source"].lower()
    flat = " ".join(src.split())
    assert "NACA Report 1135" in flat and "public domain" in flat and "secondary" in flat and "Read and deliberately not used" in flat
    assert sum(len(bm[k]["rows"]) for k in ("naca1135_table1_subsonic", "naca1135_table2_supersonic")) == 9


# =====================================================================================================================
# 99. Private book values (V6) — read from tests/book_values_ch15.json; skipped when the file is absent
# =====================================================================================================================
def _dec(s):
    return len(s.split(".")[1]) if "." in s else 0


@book_only
def test_book_V6_tables_15_1_and_15_2_within_half_a_unit_except_the_misprints():  # V6 — slips #9, #10
    bv = book()
    t1 = bv["table_15_1_isentropic_gamma_1.4"]
    miss = []
    for row in t1["rows"]:
        M = float(row[0])
        if M == 0.0:
            assert row[4] == "inf" and math.isinf(GD.area_ratio(0.0))
            continue
        s = GD.isentropic_state(M, 1.4)
        ours = [M, s["p_over_p0"], s["rho_over_rho0"], s["T_over_T0"], s["A_over_Astar"]]
        for k in range(1, 5):
            if abs(ours[k] - float(row[k])) > 0.5 * 10.0 ** (-_dec(row[k])) * 1.02:
                miss.append((row[0], t1["columns"][k]))
    assert len(t1["rows"]) == 251 and sorted(miss) == sorted((k.split("=")[1], v["column"]) for k, v in t1["misprints"].items()), miss
    for k, v in t1["misprints"].items():
        M = float(k.split("=")[1])
        ours = GD.isentropic_state(M, 1.4)["rho_over_rho0" if v["column"] == "rho/rho0" else "A_over_Astar"]
        assert abs(ours - v["formula"]) < 5e-5 and abs(ours - float(v["printed"])) > 0.5
    t2 = bv["table_15_2_normal_shock_gamma_1.4"]
    miss, truncated = [], 0
    for row in t2["rows"]:
        M1 = float(row[0])
        s = GD.normal_shock(M1, 1.4)
        ours = [M1, s["M2"], s["p_ratio"], s["T_ratio"], s["p0_ratio"]]
        for k in range(1, 5):
            units = abs(ours[k] - float(row[k])) / 10.0 ** (-_dec(row[k]))
            if units > 0.6:
                miss.append((row[0], t2["columns"][k]))
            elif units > 0.5 * 1.001:
                truncated += 1                                                  # printed digit truncated rather than rounded (0.51–0.56 of a unit)
    assert len(t2["rows"]) == 101 and sorted(miss) == sorted((k.split("=")[1], v["column"]) for k, v in t2["misprints"].items()), miss
    assert truncated == 6


@book_only
def test_book_V6_worked_examples_15_1_and_15_2():  # V6 — Examples 15.1 and 15.2 (the book reads tables to three digits: 1 %)
    bv = book()
    e = bv["example_15_1_rocket_nozzle"]
    i, pr, d = e["inputs"], e["printed"], e["derived_by_analyst_exact_functions"]
    n = ch15.rocket_nozzle(i["thrust_N"], 1e3 * i["p0_kPa"], i["T0_K"], i["p_ambient_Pa_from_appendix"], GD.PerfectGas(i["gamma"], i["R"]))
    for ours, key in ((n["M_e"], "Me"), (n["Ae_over_At"], "Ae_over_Astar"), (n["T_e"], "Te_K"), (n["rho_e"], "rho_e"), (n["u_e"], "ue_m_per_s"), (n["A_e"], "Ae_m2"), (n["A_throat"], "Astar_m2")):
        assert rel(ours, pr[key]) < 0.01, (key, ours, pr[key])
    for ours, key in ((n["M_e"], "Me"), (n["Ae_over_At"], "Ae_over_Astar"), (n["T_e"], "Te"), (n["rho_e"], "rho_e"), (n["u_e"], "ue"), (n["A_e"], "Ae"), (n["A_throat"], "Astar")):
        assert rel(ours, d[key]) < 1e-9, key
    assert rel(i["p_ambient_Pa_from_appendix"] / (1e3 * i["p0_kPa"]), pr["pe_over_p0"]) < 1e-3
    e = bv["example_15_2_nozzle_with_shock"]
    i, pr, d = e["inputs"], e["printed"], e["derived_by_analyst_exact_functions"]
    w = ch15.nozzle_with_shock(i["M1"], i["A_throat_m2"], i["A_exit_m2"], 1e3 * i["p0_kPa"], i["T0_K"], GD.AIR_BOOK)
    for ours, key in ((w["A_shock"], "A1_m2"), (w["M2"], "M2"), (w["p02"] / 1e3, "p02_kPa"), (w["A_star2"], "A2star_m2"), (w["T3"], "T3_K"), (w["p02"] / (1e3 * i["p0_kPa"]), "p02_over_p01")):
        assert rel(ours, pr[key]) < 0.01, (key, ours, pr[key])
    for ours, key in ((w["A_shock"], "A1"), (w["M2"], "M2"), (w["p02"] / 1e3, "p02_kPa"), (w["A_star2"], "A2star"), (w["M3"], "M3"), (w["T3"], "T3")):
        assert rel(ours, d[key]) < 1e-9, key
    back = GD.laval_shock_location(w["p3"] / (1e3 * i["p0_kPa"]), i["A_exit_m2"] / i["A_throat_m2"])          # the same case found from its back pressure
    assert rel(back["M1"], i["M1"]) < 1e-7


@book_only
def test_book_V6_numbers_quoted_in_the_text_and_figures():  # V6 — §15.1, §15.2, §15.6, §15.8, §15.9; slips #8, #15; trap T3
    bv = book()
    s1 = bv["sec_15_1_introduction"]
    assert GD.AIR_BOOK.R == s1["air_constants_eq_15_3"]["R"] and GD.AIR_BOOK.gamma == s1["air_constants_eq_15_3"]["gamma"]
    assert rel(s1["air_constants_eq_15_3"]["Cp"] / s1["air_constants_eq_15_3"]["Cv"], 1.4003) < 1e-4                    # the rounded constants are not mutually exact (T4)
    lo, hi = s1["transonic_range"]
    assert ch15.mach_regime(lo) == "transonic" and ch15.mach_regime(hi) == "transonic" and ch15.mach_regime(0.999 * lo) == "subsonic"
    Mh = 0.5 * (s1["hypersonic_threshold"] + 5.0)                                     # between the book's threshold and the customary one: the keyword decides (T3)
    assert ch15.mach_regime(Mh) == "supersonic" and ch15.mach_regime(Mh, hypersonic=s1["hypersonic_threshold"]) == "hypersonic"
    assert ch15.compressibility_check(0.99 * s1["incompressible_mach_limit"] * 340.29)["incompressible_ok"]
    s2 = bv["sec_15_2_acoustics"]
    assert rel(1.4 * TH.P_ATM, s2["rho0_c2_Pa_approx"]) < 0.02 and rel(ch15.perfect_gas_sound_speed(288.15), s2["sound_speed_air_15C_m_per_s"]) < 0.005
    assert AC.condensation(s2["ordinary_sound_pressure_Pa_order"], 1.225, 340.3) < s2["condensation_typically_less_than"]
    s6 = bv["sec_15_6_normal_shock"]
    assert abs(GD.strong_shock_limits(1.4)["M2"] - s6["fig_15_8"]["asymptote_printed"]) < 5e-4
    assert rel(GD.normal_shock_entropy(5.0), s6["fig_15_9"]["value_at_M1_5_approx"]) < 0.02
    th = s6["shock_thickness"]
    assert rel(ch15.shock_thickness_estimate(th["velocity_jump_m_per_s"], th["nu_m2_per_s"]), th["thickness_m"]) < 1e-12
    fig = s6["fig_15_11"]
    U1, U2 = ch15.becker_shock_roots(2.0, 1.4)
    assert abs(U1 - fig["U1"]) < 1e-6 and abs(U2 - fig["U2"]) < 1e-5                    # the printed plateaus are those of M1 = 2 …
    Uc1, Uc2 = ch15.becker_shock_roots(fig["M1_in_caption"], 1.4)
    assert abs(Uc1 - fig["U1"]) > 0.02 and abs(Uc2 - fig["U2"]) > 0.02                  # … not of the captioned Mach number (slip #8)
    s8 = bv["sec_15_8_duct"]["derived_by_analyst"]
    for M1, v in s8["choking_limits"].items():
        lim = GD.duct_choking_limits(float(M1))
        assert rel(lim["f_max"], v["f_max"]) < 1e-10 and rel(lim["q_max"], v["q_max"]) < 1e-10
    for M1, v in s8["differentials_vs_finite_differences"].items():
        M = float(M1)
        cf, cq = GD.influence_coefficients(M, df=1.0), GD.influence_coefficients(M, dq=1.0)
        assert rel(cf["dp_over_p"], v["analyst_form_1_plus"]) < 1e-12 and rel(cf["du_over_u"], v["printed_1_over_1_minus_M2"]) < 1e-12
        assert rel(cq["dT_over_T"], v["printed_(1-gM2)/(1-M2)"]) < 1e-12 and rel(cq["u_du_over_h"], v["analyst_(g-1)M2/(1-M2)"]) < 1e-12
        assert rel(GD.influence_coefficients(M, df=1.0, printed_15_47=True)["dp_over_p"], v["printed_form_1_minus"]) < 1e-12
        assert rel(GD.influence_coefficients(M, dq=1.0, printed_15_48=True)["u_du_over_h"], v["printed_(1-g)M2/(1-M2)"]) < 1e-12
        assert rel(cf["dp_over_p"], v["dp_over_p1_per_df_numeric"]) < 1e-4 and rel(cq["u_du_over_h"], v["u_du_over_h1_per_dq_numeric"]) < 1e-4
    f16 = bv["sec_15_8_duct"]["fig_15_16"]                                              # slip #15: the supersonic heating labels are the decimal values
    turns = [brentq(lambda M: GD.duct_choking_limits(M)["q_max"] - qd, 1.0001, 6.0) for qd in f16["q_values_subsonic_labels"]]
    assert abs(turns[0] - 1.4) < 0.05 and max(turns) < 2.0                              # with the decimal values the curves turn where the figure shows them
    whole = [brentq(lambda M: GD.duct_choking_limits(M)["q_max"] - float(lab.split("=")[1]), 1.0001, 9.0) for lab in f16["supersonic_positive_labels_as_printed"]]
    assert min(whole) > 2.5 and max(whole) > 4.0                                        # with the printed whole numbers they would turn far to the right, the last beyond the axis
    s9 = bv["sec_15_9_planar_waves"]
    assert abs(math.degrees(GD.max_deflection(2.0)["delta_max"]) - s9["delta_max_for_M1_2_deg_printed"]) < 0.5
    d9 = s9["derived_by_analyst"]
    for M1, (dm, sm) in d9["delta_max_deg_and_sigma_deg"].items():
        md = GD.max_deflection(float(M1))
        assert abs(math.degrees(md["delta_max"]) - dm) < 1e-8 and abs(math.degrees(md["sigma_max"]) - sm) < 5e-6      # the analyst located sigma with a minimiser (a flat maximum)
    for M1, nu in d9["prandtl_meyer_deg"].items():
        assert abs(math.degrees(GD.prandtl_meyer(float(M1))) - nu) < 1e-10
    for deg, ratio in d9["ratio_exact_over_linear_M1_2"].items():
        assert rel(GD.weak_wave_error(2.0, math.radians(float(deg)))["ratio"], ratio) < 1e-10


@book_only
def test_book_V6_exercises_with_checkable_answers():  # V6 — Exercises 15.3, 15.7–15.9, 15.12–15.15, 15.17; slip #17
    ex = book()["exercises"]
    e = ex["15.3"]
    assert rel(GD.mach_for_density_change(e["criterion_percent"] / 100.0), e["derived_by_analyst"]["b_M_exact"]) < 1e-10
    assert rel(math.sqrt(2 * e["criterion_percent"] / 100.0), e["derived_by_analyst"]["b_M_leading_term"]) < 1e-12
    e = ex["15.7"]
    i, d = e["inputs"], e["derived_by_analyst"]
    T0 = ch15.stagnation_temperature(i["T1_K"], i["u1"], cp=1004.0)
    T2 = T0 - i["u2"] ** 2 / (2 * 1004.0)
    p01, p02 = i["p1_kPa"] * (T0 / i["T1_K"]) ** 3.5, i["p2_kPa"] * (T0 / T2) ** 3.5
    assert rel(T0, d["T0"]) < 1e-10 and rel(T2, d["T2"]) < 1e-10 and rel(p01, d["p01_kPa"]) < 1e-10 and rel(p02, d["p02_kPa"]) < 1e-10 and rel(p01 - p02, d["loss_kPa"]) < 1e-9
    M1 = GD.mach_from_speed(i["u1"], i["T1_K"], GD.AIR_BOOK)
    assert rel(i["p1_kPa"] * GD.stagnation_pressure_ratio(M1), p01) < 2e-3            # (15.29) with gamma = 1.4 and the rounded Cp differ in the fourth digit (T4)
    assert rel(GD.entropy_change_from_stagnation(p02 / p01, 1.0, GD.AIR_BOOK), d["dS_J_per_kgK"]) < 1e-9
    assert 0.03 < 1 - (p01 - p02) / e["printed_loss_kPa_nearly"] < 0.04                 # slip #17: the quoted loss is about 3.5 % above ours — reported, not resolved
    e = ex["15.8"]
    d = e["derived_by_analyst"]
    ms = GD.moving_shock(TH.P_ATM, 288.15, p2_over_p1=1e3 * e["inputs"]["p2_kPa"] / TH.P_ATM, gas=GD.AIR_BOOK)
    assert rel(ms["M_s"], d["M1"]) < 1e-9 and rel(ms["c1"], d["c1"]) < 1e-9 and rel(ms["u_shock"], d["shock_speed"]) < 1e-9 and rel(ms["u_gas"], d["gas_speed_behind"]) < 1e-9
    e = ex["15.9"]
    K = ch15.blast_constant_crude(1.4, "hemisphere")
    assert rel(K, e["derived_by_analyst"]["K"]) < 1e-12 and abs(100 * (K / e["K_similarity_gamma_1.4_printed"] - 1) - e["derived_by_analyst"]["error_vs_1.033_percent"]) < 1e-6
    assert rel(math.degrees(GD.max_deflection(1e6)["delta_max"]), ex["15.12"]["derived_by_analyst_delta_max_Minf_deg"]) < 1e-9
    assert all(math.radians(a) != I["wedge_half_angle"] for a in ex["15.12"]["half_angles_deg"])
    e = ex["15.13"]
    i, d = e["inputs"], e["derived_by_analyst"]
    dl = math.radians(i["turn_deg"])
    for br in ("weak", "strong"):
        o = GD.oblique_shock(i["M1"], delta=dl, branch=br)
        assert close([math.degrees(o["sigma"]), o["M2"], o["p_ratio"]], d[br], 1e-9)
    we = GD.weak_wave_error(i["M1"], dl)
    assert rel(we["linear"], d["linear_dp_over_p"]) < 1e-10 and rel(we["exact"], d["exact_dp_over_p"]) < 1e-10 and rel(we["rel_error"], d["linear_rel_err"]) < 1e-9
    assert rel(math.degrees(GD.mach_angle(i["M1"])), d["mach_angle_deg"]) < 1e-12
    e = ex["15.14"]
    i, d = e["inputs"], e["derived_by_analyst"]
    al = math.radians(i["alpha_deg"])
    p = ch15.flat_plate_supersonic(al, i["M_inf"], method="linear")
    Lf = p["L_over_pb"] * TH.P_ATM * i["chord_m"]
    assert rel(Lf, d["L_N_per_m"]) < 1e-9 and rel(p["D_over_pb"] * TH.P_ATM * i["chord_m"], d["D_N_per_m"]) < 1e-9 and rel(p["CL"], d["CL"]) < 1e-10 and rel(p["CD"], d["CD"]) < 1e-10
    assert rel(Lf * math.cos(al), d["L_with_cos"]) < 1e-9 and rel(Lf * math.sin(al), d["D_with_sin"]) < 1e-9
    assert rel(GD.mass_flux_function(1.0), ex["15.15"]["derived_by_analyst_mdot_sqrt(RT0)_over_p0A"]) < 1e-12
    d = ex["15.17"]["derived_by_analyst"]
    tau, M = 0.05, 2.0
    r = ch15.ackeret_airfoil(lambda x: tau * np.sin(PI * np.asarray(x, float)), lambda x: 0.0 * np.asarray(x, float), 0.0, M)       # the exercise's shape: sine above, flat below
    beta = math.sqrt(M * M - 1)
    assert rel(r["CD"] * beta / (2 * tau ** 2), d["pi^2/2"]) < 1e-5 and abs(r["CL"] * beta / (2 * tau)) < 1e-10


@book_only
def test_book_V6_no_book_value_in_the_public_files_of_this_chapter():  # V6 — rule 9: the private file lists what must not appear in public files
    bv = book()
    public = [Path(__file__), REF / "make_refs.py", REF / "SOURCES.md", REF / "benchmarks.json", ROOT / "tests" / "ch15_verify_figures.py"]
    text = "\n".join(p.read_text(encoding="utf-8") for p in public if p.exists())
    hits = [s for s in bv["_forbidden_public"] if s in text]
    assert not hits, hits
    for rule in bv["_forbidden_public_regex"]:
        assert not re.search(rule["regex"], text), rule["why"]
    e1, e2 = bv["example_15_1_rocket_nozzle"]["inputs"], bv["example_15_2_nozzle_with_shock"]["inputs"]
    assert I["thrust"] != e1["thrust_N"] and I["T0_rocket"] != e1["T0_K"] and I["p0"] != 1e3 * e2["p0_kPa"] and I["T0"] != e2["T0_K"] and I["A_throat"] != e2["A_throat_m2"]
