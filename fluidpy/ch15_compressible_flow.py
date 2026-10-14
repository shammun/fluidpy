"""Chapter 15 — Compressible Flow (Kundu, Cohen & Dowling 5e, §15.1–§15.10, Eqs. (15.1)–(15.59)).

Three conservation laws and the gas law, written between two stations in terms of one number, the Mach number M — and
the factor 1 − M² that appears in every answer changes sign at the speed of sound.

**This module re-exports every public name of the two new core modules**

* ``fluidpy.core.acoustics`` (``AC``) — sound as the linear limit (15.4)–(15.7), the wave equation (15.12)–(15.13),
  d'Alembert's solution (15.14)–(15.17), the finite-amplitude simple wave (ours), wave fronts of a moving source;
* ``fluidpy.core.gasdynamics`` (``GD``) — the perfect gas (15.1)–(15.3), the duct laws (15.18)–(15.26), stagnation and
  sonic states (15.27)–(15.31), the area rule (15.32)–(15.34), normal shocks (15.35)–(15.44), nozzles at different back
  pressures (§15.7), ducts with friction and heat (15.45)–(15.48), Mach angle (15.49), oblique shocks (15.50)–(15.55),
  the Prandtl–Meyer function (15.56);

and the earlier-chapter functions the chapter restates (``mach_number`` (4.111), ``perfect_gas_sound_speed`` (1.27),
``sound_speed_from_eos`` (1.19), ``isentropic_ratios`` (1.25)–(1.26), the two entropy forms of (15.1),
``stagnation_temperature`` and the *incompressible* ``stagnation_pressure`` of ch04, ``hydraulic_jump`` and the
shallow-water ``simple_wave_evolve`` of ch07 (as ``shallow_simple_wave_evolve``), ``blast_radius`` and
``mean_free_path_air`` of ch01, ``compressibility_check`` and the 2π lift slope of ch14), so that ``ch15.<name>``
reaches everything a notebook cell or an explainer parity row needs.

What is added here: the chapter's tables (:func:`illustrative_inputs`, :func:`conventions_table`, :func:`book_slips`,
:func:`traps`), the regime classifier, the sympy engines of the derivations, the worked nozzle examples with our own
inputs, the internal structure of a shock (Becker's solution), wedge flow, the expansion fan and its staircase, curved
walls, Ackeret's thin-airfoil theory (15.57)–(15.59) beside the exact shock–expansion result, and an optional
one-dimensional Euler run.

Conventions: see the two core modules.  Station 1 upstream; subscript 0 = stagnation (but the *undisturbed* state in
§15.2); γ always an argument; angles in radians; σ, δ = shock angle, deflection (β, θ elsewhere);
``prandtl_meyer_turn(turn=)`` counts an expansion positive (the book's δ = −turn); the book's f = F/(p₁A) and q = Q/h₁
are primary, the usual Fanno/Rayleigh forms are labelled "ours".  **Three printed slips are coded corrected** with a
``printed_*`` option each (:func:`fluidpy.core.gasdynamics.influence_coefficients`).  Inputs used in examples are ours
(:func:`illustrative_inputs`), never the book's.
"""
from __future__ import annotations

from typing import Callable

import numpy as np
from scipy.optimize import brentq

from .core import acoustics as AC
from .core import gasdynamics as GD
from .core.acoustics import *  # noqa: F401,F403
from .core.gasdynamics import *  # noqa: F401,F403
from .core._util import as_scalar_if_0d
from .core.airfoils import zhukhovsky_lift_coefficient  # noqa: F401
from .core.bernoulli import stagnation_enthalpy, stagnation_pressure, stagnation_temperature  # noqa: F401
from .core.gasdynamics import PerfectGas
from .core.similarity import mach_number  # noqa: F401
from .core.thermo import (GAMMA_AIR, P_ATM, R_AIR, isentropic_pressure, isentropic_ratios,  # noqa: F401
                          perfect_gas_entropy_change, perfect_gas_entropy_change_p, perfect_gas_sound_speed,
                          sound_speed_from_eos)
from .ch01_introduction import blast_radius, mean_free_path_air  # noqa: F401
from .ch07_gravity_waves import hydraulic_jump  # noqa: F401
from .ch07_gravity_waves import simple_wave_evolve as shallow_simple_wave_evolve  # noqa: F401
from .ch14_aerodynamics import compressibility_check  # noqa: F401

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d

_REUSED = ["mach_number", "perfect_gas_sound_speed", "sound_speed_from_eos", "isentropic_ratios", "isentropic_pressure",
           "perfect_gas_entropy_change", "perfect_gas_entropy_change_p", "stagnation_enthalpy",
           "stagnation_temperature", "stagnation_pressure", "hydraulic_jump", "shallow_simple_wave_evolve",
           "blast_radius", "mean_free_path_air", "compressibility_check", "zhukhovsky_lift_coefficient"]
_OWN = ["illustrative_inputs", "conventions_table", "book_slips", "traps", "mach_regime",
        "linearised_acoustics_sympy", "duct_equations_sympy", "normal_shock_sympy", "weak_shock_series",
        "shock_structure_sympy", "tangential_momentum_sympy", "weak_wave_sympy", "prandtl_meyer_sympy",
        "rocket_nozzle", "nozzle_thrust", "nozzle_with_shock", "blast_constant_crude", "shock_thickness_estimate",
        "shock_structure_constants", "becker_shock_roots", "becker_shock_ode", "becker_shock_profile",
        "becker_shock_thickness", "wedge_flow", "expansion_fan", "fan_staircase", "curved_wall_mach_lines",
        "ackeret_pressure", "ackeret_cp", "flat_plate_supersonic", "ackeret_flat_plate", "ackeret_airfoil",
        "flat_plate_wave_pattern", "lift_slope_vs_mach", "euler1d_run", "AC", "GD"]
__all__ = sorted(set(AC.__all__) | set(GD.__all__) | set(_REUSED) | set(_OWN))


# ======================================================================================================================
# tables
# ======================================================================================================================
def illustrative_inputs() -> dict:
    """The one set of inputs every worked number, script, notebook cell and explainer of the chapter uses — ours,
    visibly not the book's.

    Book: none (the book's example and exercise values stay in the private test file).
    Returns a flat dict.
    **Gas**: ``gamma`` [–], ``R`` [J/(kg K)] (project air), ``T_air`` [K], ``p_air`` [Pa].
    **Sound** (§15.2): ``pulse_amplitude`` [Pa], ``pulse_width`` [m], ``U_stream_sub``, ``U_stream_sup`` [m/s],
    ``steep_u`` [m/s] (velocity amplitude of the finite wave), ``steep_wavelength`` [m].
    **Nozzle** (§15.4–§15.7): ``p0`` [Pa], ``T0`` [K], ``A_throat``, ``A_inlet``, ``A_exit`` [m²], ``x_throat``,
    ``L`` [m] (arguments of ``laval_area``), ``Ae_over_At`` [–], ``pB_over_p0_shock`` [–] (a back pressure that puts
    a shock in the diverging part), ``A_conv_exit`` [m²] (convergent nozzle).
    **Rocket** (§15.5): ``thrust`` [N], ``p0_rocket`` [Pa], ``T0_rocket`` [K], ``p_exit_rocket`` [Pa],
    ``gamma_rocket`` [–], ``R_rocket`` [J/(kg K)].
    **Shock** (§15.6): ``M1`` [–], ``M_s`` [–] (moving shock), ``M1_structure`` [–].
    **Duct** (§15.8): ``M1_duct_sub``, ``M1_duct_sup`` [–], ``f`` [–], ``q`` [–].
    **Planar waves** (§15.9): ``M_flight`` [–], ``h_flight`` [m], ``M_wedge`` [–], ``wedge_half_angle`` [rad],
    ``M_fan`` [–], ``fan_turn`` [rad].
    **Plate** (§15.10): ``M_inf`` [–], ``alpha`` [rad], ``chord`` [m], ``thickness_ratio`` [–].
    Assumptions: none.
    Validation: V7 — tests/test_ch15.py: test_tables_V7_inputs_conventions_slips_and_traps.
    Label: analytic.
    """
    d = np.deg2rad
    return dict(gamma=1.4, R=R_AIR, T_air=288.15, p_air=P_ATM,
                pulse_amplitude=20.0, pulse_width=0.5, U_stream_sub=170.0, U_stream_sup=510.0, steep_u=30.0,
                steep_wavelength=1.0,
                p0=5.0e5, T0=500.0, A_throat=0.010, A_inlet=0.025, A_exit=0.020, x_throat=0.4, L=1.0, Ae_over_At=2.0,
                pB_over_p0_shock=0.75, A_conv_exit=0.010,
                thrust=5.0e4, p0_rocket=3.0e6, T0_rocket=3000.0, p_exit_rocket=P_ATM, gamma_rocket=1.25, R_rocket=320.0,
                M1=2.0, M_s=1.5, M1_structure=2.0,
                M1_duct_sub=0.4, M1_duct_sup=2.5, f=0.1, q=0.2,
                M_flight=1.6, h_flight=10000.0, M_wedge=2.5, wedge_half_angle=float(d(12.0)), M_fan=2.0,
                fan_turn=float(d(15.0)),
                M_inf=2.0, alpha=float(d(5.0)), chord=1.0, thickness_ratio=0.06)


def conventions_table():
    """The conventions of the chapter beside the rival standard of most other sources ("here / elsewhere").

    Book: §15.1–§15.10 (traps T1, T2, T5, T6, T11, T13, T14 of the analysis; decision 8 of the curation).
    Returns a ``pandas.DataFrame`` with columns "item", "here", "elsewhere", "code" (nine rows).
    Assumptions: none.
    Validation: V7 — tests/test_ch15.py: test_tables_V7_inputs_conventions_slips_and_traps.
    Label: analytic.
    """
    import pandas as pd

    rows = [
        ("subscript 0", "three meanings: reference scale (§15.1), undisturbed state of the acoustic split (§15.2), "
         "stagnation state (§15.4 on)", "p_t, T_t ('total') for the stagnation state in most tables",
         "p0, T0 = stagnation in GD; rho0, p0 = undisturbed in AC"),
        ("asterisk", "A* = control surface in (15.20), (15.23); sonic throat area from (15.31) on; A in (15.46) is a "
         "pure number", "starred Fanno/Rayleigh states = sonic state of the same frictional/heated process",
         "A_over_Astar; duct_parameter_A; fanno_ratios / rayleigh_ratios say which star"),
        ("stations", "1 upstream → 2 downstream (3 = exit in the nozzle-with-shock example)",
         "0 → 1 on some reference pages", "M1, M2; normal_shock refuses M1 < 1"),
        ("oblique-shock angles", "σ shock angle, δ deflection, both from the upstream flow direction",
         "β (wave angle) and θ (deflection)", "sigma=, delta="),
        ("turning angle", "δ > 0 when the wall turns into the flow (compression); an expansion is δ < 0 and raises ν",
         "most tables count the expansion turn positive", "prandtl_meyer_turn(turn=) expansion +; delta_book = −turn"),
        ("friction and heat", "lumped f = F/(p₁A) and q = Q/h₁ referred to the inlet state",
         "4f_F L/D (Fanno) and T₀/T₀* (Rayleigh)", "f=, q=; fanno_ratios, rayleigh_ratios, fanno_length_to_f (ours)"),
        ("wall force F", "(15.24): force of the walls on the fluid, +x; (15.45): total retarding force (f > 0)",
         "thrust = force of the fluid on the walls", "wall_force (on the fluid); f (retarding)"),
        ("entropy", "S per unit mass here (s in ch01)", "s", "dS, dS_over_cv, ds_over_cp"),
        ("overloaded letters", "f, g = d'Alembert's shapes (§15.2) but f = friction parameter (§15.8); μ = Mach angle "
         "and viscosity; ν = Prandtl–Meyer function and kinematic viscosity; δ = deflection and shock thickness; "
         "b = chord; c = sound speed; U = u/V inside the shock", "c = chord, b = span in ch14",
         "mu (angle), nu (function), delta (angle), chord"),
    ]
    return pd.DataFrame(rows, columns=["item", "here", "elsewhere", "code"])


def book_slips():
    """Suspected slips of the chapter and what the code does about each (no book number is reproduced).

    Book: §15.2–§15.10 (analysis §9a, 17 rows).  ``kind``: "slip" = an error on the page; "loose" = an inconsistency
    or imprecision.  ``changes_physics`` marks the three that are coded corrected with a ``printed_*`` option.
    Returns a ``pandas.DataFrame`` with columns "n", "kind", "where", "printed", "correct", "code", "changes_physics".
    Assumptions: none.
    Validation: V7 — tests/test_ch15.py: test_tables_V7_inputs_conventions_slips_and_traps.
    Label: analytic.
    """
    import pandas as pd

    rows = [
        (1, "slip", "sentence after (15.6)", "c is larger in 'monotonic' gases", "monatomic (the argument is about γ)", "—", False),
        (2, "loose", "(15.9), (15.11), (15.14), (15.16)", "P + p′ and ρ_o", "p₀ + p′ and ρ₀ as in (15.4)", "rho0", False),
        (3, "loose", "sentence before (15.10)", "'using (15.6) to eliminate ρ′'", "using (15.5), p′ ≅ c²ρ′", "linearised_acoustics_sympy", False),
        (4, "slip", "(15.32), third member", "u du = −dp/ρ = c² dρ/ρ", "u du = −dp/ρ = −c² dρ/ρ", "influence_coefficients(printed_15_32=True)", True),
        (5, "loose", "§15.5, constant-density remark", "dA/dx + du/dx = 0", "(1/A)dA/dx + (1/u)du/dx = 0", "influence_coefficients(M=0)", False),
        (6, "slip", "Examples 15.1, 15.2", "lower-case m for the exit Mach number; a ratio written as an equality of temperatures", "M_e; T₃/T₀₃", "—", False),
        (7, "loose", "shock-structure equation", "dimensional quadratic set ≡ (U₁ − U)(U − U₂) before dividing by V²", "≡ V²(U₁ − U)(U − U₂)", "shock_structure_sympy", False),
        (8, "slip", "caption of Fig. 15.11", "the upstream Mach number quoted for the two plateau speeds", "the plateaus drawn belong to M₁ = 2 for γ = 1.4 (U₁ = 28/33, U₂ = 7/22)", "becker_shock_roots(2.0)", False),
        (9, "slip", "Table 15.1, two entries", "a wrong leading digit in one density and one area entry", "(15.30), (15.31) evaluated", "isentropic_table", False),
        (10, "slip", "Table 15.2, one entry", "last digit of one temperature ratio", "(15.42) evaluated", "normal_shock_table", False),
        (11, "slip", "(15.47), first relation", "−dp/p₁ = [1 − (γ − 1)M²]/(1 − M²) df", "−dp/p₁ = [1 + (γ − 1)M²]/(1 − M²) df", "influence_coefficients(printed_15_47=True)", True),
        (12, "slip", "(15.48), second relation", "u du/h₁ = (1 − γ)M²/(1 − M²) dq", "u du/h₁ = (γ − 1)M²/(1 − M²) dq", "influence_coefficients(printed_15_48=True)", True),
        (13, "slip", "§15.10, before the plate pressures", "'found from (15.47)'", "from (15.57)", "ackeret_pressure", False),
        (14, "loose", "§15.8", "'once M₂ is calculated from (15.45)'", "from (15.46)", "duct_outlet_mach", False),
        (15, "loose", "Fig. 15.16, supersonic heating curves", "parameter labels without their decimal point", "the same decimal values as on the subsonic side", "duct_choking_limits", False),
        (16, "loose", "Exercise 15.1(b)", "'the simplified equation in part b)'", "part a)", "—", False),
        (17, "loose", "Exercise 15.7", "a quoted stagnation-pressure loss that the chapter's own constants do not reproduce", "reported as a discrepancy, not resolved", "stagnation_pressure_ratio", False),
    ]
    return pd.DataFrame(rows, columns=["n", "kind", "where", "printed", "correct", "code", "changes_physics"])


def traps():
    """The chapter's traps: places where the page is right and the reader goes wrong (analysis §9b, T1–T22).

    Book: §15.1–§15.10.
    Returns a ``pandas.DataFrame`` with columns "id", "trap", "guard" (what the code or the notebook does).
    Assumptions: none.
    Validation: V7 — tests/test_ch15.py: test_tables_V7_inputs_conventions_slips_and_traps.
    Label: analytic.
    """
    import pandas as pd

    rows = [
        ("T1", "subscript 0 means reference scale, undisturbed state or stagnation state", "AC uses rho0, p0 for the undisturbed state; GD uses p0, T0 for stagnation"),
        ("T2", "three different stars and an A that is not an area", "A_over_Astar vs duct_parameter_A; Fanno/Rayleigh docstrings name their star"),
        ("T3", "regime boundaries are nominal", "mach_regime takes the thresholds as keywords"),
        ("T4", "γ and the gas constants: rounded constants are not mutually exact", "PerfectGas(gamma, R); AIR_BOOK for book numbers"),
        ("T5", "letters reused (f, g, F, Q, q, U, V, μ, ν, δ, σ, b, c, m)", "conventions_table()"),
        ("T6", "which Mach number is upstream: (15.40) is its own inverse", "normal_shock raises for M1 < 1 unless allow_expansion=True"),
        ("T7", "stagnation quantities depend on the frame", "moving_shock returns T0_behind_lab, p0_behind_lab"),
        ("T8", "the area ratio has two inverses, and A* is not conserved across a shock", "mach_from_area_ratio(branch=); laval_shock_location"),
        ("T9", "the matched-nozzle thrust formula silently drops (p_e − p_a)A_e", "nozzle_thrust keeps the pressure term"),
        ("T10", "'choked' does not mean the mass flow cannot be raised", "choked_mass_flow ∝ p0/√T0"),
        ("T11", "sign and meaning of the wall force; f is not 4fL/D and does not add", "wall_force (on the fluid); fanno_length_to_f"),
        ("T12", "which root of (15.46)", "duct_outlet_mach(branch='continuous' | 'shocked' | 'plus' | 'minus')"),
        ("T13", "oblique-shock conventions; 'weak' is a root, not a small wave", "oblique_shock(delta=, sigma=, branch=)"),
        ("T14", "the sign of the turning angle", "prandtl_meyer_turn(turn=) expansion +; weak_wave_pressure δ + into the flow"),
        ("T15", "§15.10 geometry: b is the chord; upper surface expands", "flat_plate_supersonic returns p_upper < 1 < p_lower"),
        ("T16", "how small is 'small'?", "weak_wave_error; flat_plate_supersonic(method='shock_expansion')"),
        ("T17", "adiabatic is not isentropic", "T0_ratio = 1 but p0_ratio < 1 in normal_shock, duct_state_ratios"),
        ("T18", "the two 'paradoxes' are regime statements", "rayleigh_temperature_sign; influence_coefficients"),
        ("T19", "acoustics: (15.14) is an indefinite integral; u′ is not the wave speed", "velocity_from_pressure(t0=); condensation"),
        ("T20", "shock structure: Pr″ is built on μ″ = (4/3)μ + μ_v; η has an arbitrary origin", "becker_shock_profile(eta0=)"),
        ("T21", "the hydraulic-jump analogy is partial", "hydraulic_jump (energy lost) beside normal_shock (p0 lost)"),
        ("T22", "the blast wave of the exercise is hemispherical", "blast_constant_crude(geometry='hemisphere')"),
    ]
    return pd.DataFrame(rows, columns=["id", "trap", "guard"])


# ======================================================================================================================
# §15.1 — regimes
# ======================================================================================================================
def mach_regime(M: float, *, transonic: tuple[float, float] = (0.8, 1.2), hypersonic: float = 5.0) -> str:
    """Name the flow regime that goes with a Mach number.

    Book: §15.1 (the five nominal classes: incompressible M = 0, subsonic 0 < M < 1, transonic M ≈ 1, supersonic
    M > 1, hypersonic M ≫ 1).  The band and the hypersonic threshold are conventions (trap T3) and therefore
    keywords; the default hypersonic threshold 5 is the customary one — pass the value a source uses.
    Parameters: M [–] ≥ 0; transonic (lower, upper) band [–]; hypersonic threshold [–].
    Returns "incompressible", "subsonic", "transonic", "supersonic" or "hypersonic".  Rule: M = 0 → incompressible;
    inside the closed transonic band → transonic; below it → subsonic; at or above the hypersonic threshold →
    hypersonic; otherwise supersonic.  Raises ``ValueError`` for NaN or negative input.
    Assumptions: a classifier, not a verdict — "incompressible" in practice means M < 0.3
    (:func:`compressibility_check`).
    Validation: V1 — tests/test_ch15.py: test_mach_regime_V1_both_sides_of_every_threshold.
    Label: analytic.
    """
    M = float(M)
    if not np.isfinite(M) or M < 0.0:
        raise ValueError("Mach number must be a finite number >= 0")
    if M == 0.0:
        return "incompressible"
    if transonic[0] <= M <= transonic[1]:
        return "transonic"
    if M < transonic[0]:
        return "subsonic"
    return "hypersonic" if M >= hypersonic else "supersonic"


# ======================================================================================================================
# sympy engines of the derivations
# ======================================================================================================================
def linearised_acoustics_sympy() -> dict:
    """Linearise continuity and the Euler equation about a uniform stream and combine them into the wave equation
    (one space dimension, a bookkeeping parameter ε on every primed quantity).

    Book: §15.2, Eqs. (15.4) (the split), (15.8), (15.9) (linearised laws and what is dropped), (15.10), (15.11)
    (with p′ = c²ρ′ (15.5)), (15.12) and its general solution (15.17).
    Returns dict of sympy objects: ``symbols``; ``continuity_linear`` (coefficient of ε: ρ′_t + Uρ′_x + ρ₀u′_x,
    (15.8)); ``continuity_dropped`` (coefficient of ε²: (ρ′u′)_x); ``momentum_linear`` (u′_t + Uu′_x + p′_x/ρ₀,
    (15.9)); ``momentum_dropped`` (ε²: u′u′_x − ρ′p′_x/ρ₀²); ``eq_15_10``, ``eq_15_11``; ``wave_operator`` ((15.12)
    in one dimension applied to p′); ``cross_residual`` (D(15.10) − ρ₀∂ₓ(15.11) − wave operator, = 0);
    ``dalembert_residual`` ((15.12) applied to f(x − (c + U)t) + g(x + (c − U)t), = 0).
    Assumptions: uniform ρ₀, p₀, U, c; no gravity; inviscid.
    Validation: V2 — tests/test_ch15.py: test_wave_equation_V2_chapter_engine_agrees.
    Label: symbolic.
    """
    import sympy as sp

    x, t, eps = sp.symbols("x t epsilon", real=True)
    rho0, p0, U, c = sp.symbols("rho_0 p_0 U c", positive=True)
    r, u, p = (sp.Function(n)(x, t) for n in ("rho1", "u1", "p1"))
    rho, vel, pres = rho0 + eps * r, U + eps * u, p0 + eps * p  # Eq. (15.4)
    cont = sp.expand(sp.diff(rho, t) + sp.diff(rho * vel, x))  # (4.7)
    mom = sp.diff(vel, t) + vel * sp.diff(vel, x) + sp.diff(pres, x) / rho  # (4.41), g = 0
    # Taylor coefficients in ε by differentiation (1/(ρ₀ + ερ′) ≈ 1/ρ₀ − ερ′/ρ₀² + …)
    mom_ser = sp.expand(sum(sp.diff(mom, eps, k).subs(eps, 0).doit() * eps ** k / sp.factorial(k) for k in range(3)))
    D = lambda q: sp.diff(q, t) + U * sp.diff(q, x)  # noqa: E731
    e10 = D(p) / c ** 2 + rho0 * sp.diff(u, x)  # Eq. (15.10)
    e11 = D(u) + sp.diff(p, x) / rho0  # Eq. (15.11)
    wave = sp.expand(D(D(p)) / c ** 2 - sp.diff(p, x, 2))  # Eq. (15.12)
    f, g = sp.Function("f"), sp.Function("g")
    sol = f(x - (c + U) * t) + g(x + (c - U) * t)  # Eq. (15.17)
    dres = sp.simplify((D(D(sol)) / c ** 2 - sp.diff(sol, x, 2)).doit())
    return dict(symbols=dict(x=x, t=t, eps=eps, rho0=rho0, p0=p0, U=U, c=c, rho1=r, u1=u, p1=p),
                continuity_linear=cont.coeff(eps, 1), continuity_dropped=sp.factor(cont.coeff(eps, 2)),
                momentum_linear=mom_ser.coeff(eps, 1), momentum_dropped=mom_ser.coeff(eps, 2), eq_15_10=e10,
                eq_15_11=e11, wave_operator=wave,
                cross_residual=sp.simplify(sp.expand(D(e10) - rho0 * sp.diff(e11, x)) - wave), dalembert_residual=dres)


def duct_equations_sympy() -> dict:
    """The differential forms of the duct laws, checked symbolically.

    Book: §15.3, Eqs. (15.19) (logarithmic differential of (15.18)), (15.25) (the two members of the momentum law),
    (15.26); §15.5, (15.33) → (15.34).
    Returns dict of sympy expressions that must all be zero: ``eq_15_19`` (d ln(ρuA)/dx − (ρ′/ρ + u′/u + A′/A)),
    ``eq_15_25`` (−d(pA)/dx + p dA/dx − (−A dp/dx)), ``eq_15_26`` ((ρuA u′ + A p′)/(ρA) − (u u′ + p′/ρ)),
    ``eq_15_34`` ((15.19) with dρ/ρ = −M² du/u solved for du/u, minus −(dA/A)/(1 − M²)), and ``eq_15_33_printed``
    (the same with the printed sign of (15.32), dρ/ρ = +M² du/u: **not** zero — it is the size of the slip).
    Assumptions: steady quasi-one-dimensional flow.
    Validation: V2 — tests/test_ch15.py: test_duct_mass_and_momentum_V2_derivation.
    Label: symbolic.
    """
    import sympy as sp

    x, M = sp.symbols("x M", positive=True)
    rho, u, A, p = (sp.Function(n, positive=True)(x) for n in ("rho", "u", "A", "p"))
    e19 = sp.simplify(sp.diff(sp.log(rho * u * A), x) - (sp.diff(rho, x) / rho + sp.diff(u, x) / u + sp.diff(A, x) / A))
    e25 = sp.simplify(-sp.diff(p * A, x) + p * sp.diff(A, x) + A * sp.diff(p, x))
    e26 = sp.simplify((rho * u * A * sp.diff(u, x) + A * sp.diff(p, x)) / (rho * A) - (u * sp.diff(u, x) + sp.diff(p, x) / rho))
    du, dA = sp.symbols("du_over_u dA_over_A", real=True)
    sol = sp.solve(sp.Eq(-M ** 2 * du + du + dA, 0), du)[0]  # (15.19) with (15.33)
    sol_printed = sp.solve(sp.Eq(M ** 2 * du + du + dA, 0), du)[0]
    target = -dA / (1 - M ** 2)  # Eq. (15.34)
    return dict(eq_15_19=e19, eq_15_25=e25, eq_15_26=e26, eq_15_34=sp.simplify(sol - target),
                eq_15_33_printed=sp.simplify(sol_printed - target))


def normal_shock_sympy() -> dict:
    """Solve the three jump laws for the non-trivial root and compare it with the printed relations.

    Book: §15.6, Eqs. (15.35), (15.36), (15.38) → (15.39)–(15.42) ("after some algebra").
    Method: with r = ρ₂/ρ₁ (so u₂ = u₁/r by (15.35)) and ρ₁u₁² = γp₁M₁², momentum gives p₂/p₁ = 1 + γM₁²(1 − 1/r) and
    the energy law becomes a quadratic in r whose roots are r = 1 (nothing happens) and the shock.
    Returns dict: ``rho_ratio`` (the non-trivial root), ``p_ratio``, ``T_ratio``, ``M2_squared`` (sympy expressions
    in γ, M₁); ``roots`` (both roots of the quadratic); ``residual_15_39``, ``residual_15_40``, ``residual_15_41``,
    ``residual_15_42`` (derived − printed, simplified: all zero).
    Assumptions: perfect gas, constant γ.
    Validation: V2 — tests/test_ch15.py: test_normal_shock_density_and_temperature_V2_derivation.
    Label: symbolic.
    """
    import sympy as sp

    g, M, r = sp.symbols("gamma M_1 r", positive=True)
    pr = 1 + g * M ** 2 * (1 - 1 / r)  # Eq. (15.36) with (15.35), in units of p₁
    # Eq. (15.38) in units of p₁/ρ₁: γ/(γ−1) + γM²/2 = γ/(γ−1)·(p₂/p₁)/r + γM²/(2r²)
    energy = sp.together(g / (g - 1) + g * M ** 2 / 2 - g / (g - 1) * pr / r - g * M ** 2 / (2 * r ** 2))
    roots = sp.solve(sp.numer(energy), r)
    shock = [sp.simplify(q) for q in roots if sp.simplify(q - 1) != 0][0]
    p_ratio = sp.simplify(pr.subs(r, shock))
    T_ratio = sp.simplify(p_ratio / shock)
    M2sq = sp.simplify(M ** 2 / (shock * p_ratio))  # M₂² = ρ₂u₂²/(γp₂)
    b39 = 1 + 2 * g / (g + 1) * (M ** 2 - 1)
    b40 = ((g - 1) * M ** 2 + 2) / (2 * g * M ** 2 + 1 - g)
    b41 = (g + 1) * M ** 2 / ((g - 1) * M ** 2 + 2)
    b42 = 1 + 2 * (g - 1) / (g + 1) ** 2 * (g * M ** 2 + 1) / M ** 2 * (M ** 2 - 1)
    return dict(rho_ratio=shock, p_ratio=p_ratio, T_ratio=T_ratio, M2_squared=M2sq, roots=roots,
                residual_15_39=sp.simplify(p_ratio - b39), residual_15_40=sp.simplify(M2sq - b40),
                residual_15_41=sp.simplify(shock - b41), residual_15_42=sp.simplify(T_ratio - b42))


def weak_shock_series(order: int = 4) -> dict:
    """Series of the entropy jump in m = M₁² − 1: the first two coefficients vanish, the third is (15.44a).

    Book: §15.6, Eq. (15.43) expanded (Exercise 15.5) → (15.44a), and (15.44b) with (15.39).
    Parameters: order — highest power of m kept (default 4).
    Returns dict: ``m`` (symbol), ``series`` (sympy polynomial in m), ``coefficients`` (list c₁ … c_order,
    simplified), ``residual_15_44a`` (c₃ − 2γ(γ − 1)/(3(γ + 1)²), = 0), ``pressure_coefficient`` (coefficient of
    ((p₂ − p₁)/p₁)³ after m = (γ + 1)Δ/(2γ)), ``residual_15_44b`` (that minus (γ² − 1)/(12γ²), = 0).
    Assumptions: perfect gas; |m| ≪ 1.
    Validation: V2 — tests/test_ch15.py: test_weak_shock_entropy_V2_derivation.
    Label: symbolic.
    """
    import sympy as sp

    g = sp.Symbol("gamma", positive=True)
    m = sp.Symbol("m", real=True)
    # Eq. (15.43) with M₁² = 1 + m; ln of a product written as a sum (all factors positive near m = 0)
    expr = sp.log(1 + 2 * g / (g + 1) * m) + g * sp.log(((g - 1) * (1 + m) + 2) / ((g + 1) * (1 + m)))
    ser = sp.series(expr, m, 0, order + 1).removeO()
    coeffs = [sp.simplify(ser.coeff(m, k)) for k in range(1, order + 1)]
    c3 = coeffs[2] if order >= 3 else None
    pc = sp.simplify(c3 * ((g + 1) / (2 * g)) ** 3) if c3 is not None else None
    return dict(m=m, series=sum(c * m ** (k + 1) for k, c in enumerate(coeffs)), coefficients=coeffs,
                residual_15_44a=sp.simplify(c3 - 2 * g * (g - 1) / (3 * (g + 1) ** 2)) if c3 is not None else None,
                pressure_coefficient=pc,
                residual_15_44b=sp.simplify(pc - (g ** 2 - 1) / (12 * g ** 2)) if pc is not None else None)


def shock_structure_sympy() -> dict:
    """The steady one-dimensional Navier–Stokes and energy equations inside a shock, their first integrals and
    Becker's reduction, checked symbolically.

    Book: §15.6 "Normal shock structure" (unnumbered chain): d(ρu)/dx = 0, ρu u′ + p′ = (μ″u′)′,
    ρu h′ − u p′ = μ″u′² + (kT′)′ with μ″ = (4/3)μ + μ_v;  integrals ρu = m, mu + p = μ″u′ + mV,
    m(h + ½u²) = μ″uu′ + kT′ + mI;  for Pr″ = μ″C_p/k = 1 the energy integral reads m(H − I) = (μ″)H′ with
    H = C_pT + ½u², whose only bounded solution is H = I (Becker);  then
    −[2γ/(γ + 1)](μ″/m)u u′ = −u² + [2γ/(γ + 1)]uV − 2I(γ − 1)/(γ + 1).
    Returns dict of sympy expressions that must be zero: ``momentum_integral`` (d/dx of the momentum integral minus
    the momentum equation), ``energy_integral`` (d/dx of the energy integral minus [energy equation + u × momentum
    equation]), ``becker_form`` (energy integral with k = μ″C_p minus [m(H − I) − μ″H′]), ``structure_equation``
    (the momentum integral with H = I and p = ρRT, times u/m, minus the printed quadratic form).
    Assumptions: perfect gas; μ″(x), k(x) arbitrary in the first two checks; Pr″ = 1 in the last two.
    Validation: V2 — tests/test_ch15.py: test_shock_structure_V2_first_integrals_and_roots.
    Label: symbolic.
    """
    import sympy as sp

    x = sp.Symbol("x", real=True)
    m, V, I, R, cp, g = sp.symbols("m V I R C_p gamma", positive=True)
    u, T, mu2, k = (sp.Function(n, positive=True)(x) for n in ("u", "T", "mu2", "k"))
    rho = m / u  # mass integral
    p, h = rho * R * T, cp * T
    mom_eq = rho * u * sp.diff(u, x) + sp.diff(p, x) - sp.diff(mu2 * sp.diff(u, x), x)
    en_eq = rho * u * sp.diff(h, x) - u * sp.diff(p, x) - mu2 * sp.diff(u, x) ** 2 - sp.diff(k * sp.diff(T, x), x)
    mom_int = m * u + p - mu2 * sp.diff(u, x) - m * V
    en_int = m * (h + u ** 2 / 2) - mu2 * u * sp.diff(u, x) - k * sp.diff(T, x) - m * I
    H = cp * T + u ** 2 / 2
    becker = sp.simplify(en_int.subs(k, mu2 * cp).doit() - (m * (H - I) - mu2 * sp.diff(H, x)))
    # structure equation: T from H = I, R/C_p = (γ − 1)/γ
    Tb = (I - u ** 2 / 2) / cp
    mom_b = (m * u + (m / u) * R * Tb - mu2 * sp.diff(u, x) - m * V).subs(R, cp * (g - 1) / g)
    lhs = -2 * g / (g + 1) * (mu2 / m) * u * sp.diff(u, x)
    rhs = -u ** 2 + 2 * g / (g + 1) * u * V - 2 * I * (g - 1) / (g + 1)
    struct = sp.simplify((lhs - rhs) - (2 * g / (g + 1)) * (u / m) * mom_b)
    return dict(momentum_integral=sp.simplify(sp.diff(mom_int, x) - mom_eq),
                energy_integral=sp.simplify(sp.diff(en_int, x) - (en_eq + u * mom_eq)), becker_form=becker,
                structure_equation=struct)


def tangential_momentum_sympy() -> dict:
    """Why the shock-parallel velocity does not change across an oblique shock.

    Book: §15.9 (Fig. 15.18b; Exercise 15.18): the control volume of the normal shock with a velocity (u, v);
    mass ρ₁u₁ = ρ₂u₂ and tangential momentum ρ₁u₁v₁ = ρ₂u₂v₂ (no pressure force along the shock).
    Returns dict: ``v2`` (sympy solution, = v₁), ``residual`` (v₂ − v₁ = 0).
    Assumptions: plane shock, uniform states, no shear stress on the control surface.
    Validation: V2 — tests/test_ch15.py: test_oblique_shock_normal_component_V2_derivation.
    Label: symbolic.
    """
    import sympy as sp

    r1, r2, u1, u2, v1, v2 = sp.symbols("rho_1 rho_2 u_1 u_2 v_1 v_2", positive=True)
    sol = sp.solve([sp.Eq(r1 * u1, r2 * u2), sp.Eq(r1 * u1 * v1, r2 * u2 * v2)], [r2, v2], dict=True)[0]
    return dict(v2=sol[v2], residual=sp.simplify(sol[v2] - v1))


def weak_wave_sympy() -> dict:
    """From the exact deflection relation to the weak-wave rule.

    Book: §15.9, Eq. (15.54) (Eq. (15.53) rewritten with the pressure jump of (15.50)) and its small-angle limit
    Eq. (15.55).
    Returns dict: ``residual_15_54`` ((15.53) − (15.54) with Δ = (p₂ − p₁)/p₁ from (15.50); = 0),
    ``slope`` (dΔ/dδ at σ = μ₁ from (15.54) with tan δ → δ, sin σ → 1/M₁), ``residual_15_55`` (slope −
    γM₁²/√(M₁² − 1); = 0).
    Assumptions: δ → 0 on the weak branch (σ → μ₁).
    Validation: V2 — tests/test_ch15.py: test_weak_wave_rule_V2_derivation.
    Label: symbolic.
    """
    import sympy as sp

    g, M = sp.symbols("gamma M_1", positive=True)
    s = sp.Symbol("sigma", positive=True)
    Delta = 2 * g / (g + 1) * (M ** 2 * sp.sin(s) ** 2 - 1)  # Eq. (15.50)
    t53 = 2 * sp.cot(s) * (M ** 2 * sp.sin(s) ** 2 - 1) / (M ** 2 * (g + sp.cos(2 * s)) + 2)  # Eq. (15.53)
    t54 = 2 * sp.cot(s) * (g + 1) / (2 * g) * Delta / (M ** 2 * (g + 1 - 2 * sp.sin(s) ** 2) + 2)  # Eq. (15.54)
    res54 = sp.simplify(sp.expand_trig(t53 - t54))
    D = sp.Symbol("Delta", real=True)
    # at σ = μ₁: cot μ = √(M² − 1), sin²σ = 1/M²; tan δ ≈ δ
    delta_lin = 2 * sp.sqrt(M ** 2 - 1) * (g + 1) / (2 * g) * D / (M ** 2 * (g + 1 - 2 / M ** 2) + 2)
    slope = sp.simplify(D / delta_lin)
    return dict(residual_15_54=res54, slope=slope,
                residual_15_55=sp.simplify(slope - g * M ** 2 / sp.sqrt(M ** 2 - 1)))


def prandtl_meyer_sympy() -> dict:
    """The Prandtl–Meyer integrand and its closed form: the derivative of the closed form is the integrand.

    Book: §15.9, the differential dδ = −[√(M² − 1)/M] dM/(1 + ½(γ − 1)M²) (unnumbered) and Eq. (15.56).
    Returns dict: ``integrand``, ``closed_form`` (sympy, in γ and M > 1), ``derivative_residual`` (d/dM closed form −
    integrand, evaluated exactly in the variable b = √(M² − 1); = 0), ``residual_numeric`` (largest |residual| over a few (γ, M) pairs, as a float),
    ``value_at_1`` (ν(1) = 0).
    Assumptions: perfect gas; M > 1.
    Validation: V2 — tests/test_ch15.py: test_prandtl_meyer_closed_form_V2_derivation.
    Label: symbolic.
    """
    import sympy as sp

    g, M = sp.symbols("gamma M", positive=True)
    integrand = sp.sqrt(M ** 2 - 1) / (1 + (g - 1) / 2 * M ** 2) / M
    k = (g + 1) / (g - 1)
    closed = sp.sqrt(k) * sp.atan(sp.sqrt((M ** 2 - 1) / k)) - sp.atan(sp.sqrt(M ** 2 - 1))  # Eq. (15.56)
    raw = sp.diff(closed, M) - integrand
    num = max(abs(complex(raw.subs({g: gv, M: mv}).evalf(30))) for gv in (sp.Rational(7, 5), sp.Rational(5, 3))
              for mv in (sp.Rational(11, 10), 2, 5))
    # exact check in the variable b = √(M² − 1) > 0 (M = √(1 + b²), dM/db = b/M): everything becomes rational in b
    b = sp.Symbol("b", positive=True)
    closed_b = sp.sqrt(k) * sp.atan(b / sp.sqrt(k)) - sp.atan(b)
    integrand_b = b / (1 + (g - 1) / 2 * (1 + b ** 2)) / (1 + b ** 2) * b  # integrand × dM/db
    res = sp.simplify(sp.together(sp.diff(closed_b, b) - integrand_b))
    return dict(integrand=integrand, closed_form=closed, derivative_residual=res, residual_numeric=float(num),
                value_at_1=closed.subs(M, 1))


# ======================================================================================================================
# §15.5, §15.7 — worked nozzle chains (our inputs)
# ======================================================================================================================
def nozzle_thrust(mdot, u_e, p_e, p_ambient, A_e):
    """Thrust of a nozzle: momentum flux plus the pressure term, ṁu_e + (p_e − p_a)A_e.

    Book: §15.5, Example 15.1 uses Thrust = ρ_e u_e² A_e (= ṁu_e) for a pressure-matched nozzle — the pressure term is
    silently zero there (trap T9); from the momentum law (15.24) on a control volume round the engine.
    Parameters: mdot [kg/s]; u_e exit speed [m/s]; p_e exit-plane pressure [Pa]; p_ambient [Pa]; A_e exit area [m²].
    Returns the thrust [N].
    Assumptions: steady flow; uniform exit; negligible inlet momentum (rocket).
    Validation: V1 — tests/test_ch15.py: test_rocket_nozzle_V1_thrust_and_throat_from_two_routes.
    Label: analytic.
    """
    return _S(_F(mdot) * _F(u_e) + (_F(p_e) - _F(p_ambient)) * _F(A_e))


def rocket_nozzle(thrust: float, p0: float, T0: float, p_exit: float, gas: PerfectGas | None = None) -> dict:
    """Size a pressure-matched convergent–divergent nozzle for a given thrust (the chain of Example 15.1, exact).

    Book: §15.5, Example 15.1: p_e/p₀ → M_e (15.29) → T_e/T₀ (15.28), A_e/A* (15.31) → ρ_e = p_e/RT_e,
    u_e = M_e√(γRT_e) → A_e = Thrust/(ρ_e u_e²) → A* = A_e/(A_e/A*).  The book reads a table; here the relations are
    inverted exactly.
    Parameters: thrust [N]; p0 chamber (stagnation) pressure [Pa]; T0 chamber temperature [K]; p_exit exit pressure
    [Pa] (= ambient: matched); gas (γ and R of the combustion products).
    Returns dict: ``M_e`` [–], ``T_e`` [K], ``rho_e`` [kg/m³], ``u_e`` [m/s], ``A_e``, ``A_throat`` [m²],
    ``Ae_over_At`` [–], ``mdot`` [kg/s], ``mdot_choked`` [kg/s] (independent route through the throat).
    Assumptions: isentropic flow of a perfect gas with constant γ; chamber at rest; exit pressure = ambient.
    Validation: V1 — tests/test_ch15.py: test_rocket_nozzle_V1_thrust_and_throat_from_two_routes.
    Label: analytic.
    """
    g = GD._gas(gas)
    Me = float(GD.mach_from_pressure_ratio(p_exit / p0, g.gamma))
    Te = T0 / float(GD.stagnation_temperature_ratio(Me, g.gamma))
    rho_e = p_exit / (g.R * Te)
    ue = Me * float(g.c(Te))
    Ae = thrust / (rho_e * ue ** 2)
    ar = float(GD.area_ratio(Me, g.gamma))
    return dict(M_e=Me, T_e=Te, rho_e=rho_e, u_e=ue, A_e=Ae, A_throat=Ae / ar, Ae_over_At=ar, mdot=rho_e * ue * Ae,
                mdot_choked=float(GD.choked_mass_flow(p0, T0, Ae / ar, g)))


def nozzle_with_shock(M1: float, A_throat: float, A_exit: float, p0: float, T0: float,
                      gas: PerfectGas | None = None) -> dict:
    """A convergent–divergent nozzle with a normal shock at a known Mach number in the diverging part (the chain of
    Example 15.2, exact).

    Book: §15.7, Example 15.2: M₁ → A₁/A₁* (area at the shock) → M₂, p₀₂/p₀₁ → A₂* = A_t p₀₁/p₀₂ → A₃/A₃* →
    M₃ on the **subsonic** branch → T₃/T₀, p₃ (station 3 = exit).  A* is not conserved across the shock (trap T8).
    Parameters: M1 Mach number just ahead of the shock [–] (1 ≤ M₁ ≤ M_sup(A_exit/A_throat)); A_throat, A_exit [m²];
    p0 [Pa], T0 [K] upstream stagnation state; gas.
    Returns dict: ``A_shock`` [m²], ``M2`` [–], ``p02`` [Pa], ``A_star2`` [m²], ``M3`` [–], ``T3`` [K], ``p3`` [Pa]
    (= the back pressure that holds this shock), ``mdot`` [kg/s], ``mdot_exit`` [kg/s] (from the exit state).
    Assumptions: isentropic except across the shock; adiabatic (T₀ constant throughout).
    Validation: V4 — tests/test_ch15.py: test_laval_nozzle_V4_consistent_with_shock_and_isentropic_relations.
    Label: conserved.
    """
    g = GD._gas(gas)
    A1 = A_throat * float(GD.area_ratio(M1, g.gamma))
    if A1 > A_exit * (1.0 + 1e-12):
        raise ValueError("the shock would stand outside the nozzle: A/A* at M1 exceeds A_exit/A_throat")
    s = GD.normal_shock(M1, g.gamma)
    Astar2 = A_throat / s["p0_ratio"]
    M3 = float(GD.mach_from_area_ratio(A_exit / Astar2, g.gamma, branch="subsonic"))
    T3 = T0 / float(GD.stagnation_temperature_ratio(M3, g.gamma))
    p02 = p0 * s["p0_ratio"]
    p3 = p02 / float(GD.stagnation_pressure_ratio(M3, g.gamma))
    return dict(A_shock=A1, M2=s["M2"], p02=p02, A_star2=Astar2, M3=M3, T3=T3, p3=p3,
                mdot=float(GD.choked_mass_flow(p0, T0, A_throat, g)),
                mdot_exit=float(GD.mass_flow_from_mach(M3, p02, T0, A_exit, g)))


# ======================================================================================================================
# §15.6 — blast wave (crude), shock thickness, shock structure
# ======================================================================================================================
def blast_constant_crude(gamma: float = 1.4, geometry: str = "hemisphere") -> float:
    """The constant K of R(t) = K (E/ρ₁)^{1/5} t^{2/5} from the crude strong-shock model of Exercise 15.9 (ours).

    Book: Exercise 15.9 (with Example 1.4): put the energy E uniformly behind the shock, ρ₂e₂ = E/Volume with
    e = p/((γ − 1)ρ), so p₂ = (γ − 1)E/Volume; the strong-shock limit of (15.39)–(15.41) gives
    dR/dt = u₁ = √(½(γ + 1)p₂/ρ₁); integrating from R(0) = 0:  K = [75(γ² − 1)/(16π)]^{1/5} for a hemisphere
    (Volume = (2/3)πR³) and [75(γ² − 1)/(32π)]^{1/5} for a sphere (Volume = (4/3)πR³).
    Parameters: gamma [–]; geometry "hemisphere" or "sphere".   Returns K [–].
    Crude: the true similarity solution has a non-uniform interior; the sphere constant of this model is about 9 %
    (9.4 %) below the cited similarity value (``reference/ch01/SOURCES.md``).  It is in the normalisation of this
    exercise, not that of ``reference/ch01`` (trap T22): ``ch01.blast_radius`` writes R = (E t²/(Kρ))^{1/5}, so its K
    is this function's K to the power −5.
    Assumptions: strong shock (p₂ ≫ p₁), perfect gas, uniform pressure behind the front, all the energy thermal.
    Validation: V1 — the closed form equals an integration of its own ODE (1e-8); against the cited similarity constant
    the sphere value is 9.4 % low, so the model itself is only qualitative (verification report, Open item 4) —
    tests/test_ch15.py: test_blast_constant_V1_crude_model_integrates_to_the_closed_form.
    Label: qualitative (the closed form of the crude model: analytic).
    """
    if geometry == "hemisphere":
        vol = 2.0 / 3.0 * np.pi
    elif geometry == "sphere":
        vol = 4.0 / 3.0 * np.pi
    else:
        raise ValueError('geometry must be "hemisphere" or "sphere"')
    C2 = 0.5 * (gamma + 1.0) * (gamma - 1.0) / vol  # (dR/dt)² R³ ρ₁/E
    return float((2.5 * np.sqrt(C2)) ** 0.4)


def shock_thickness_estimate(du, nu):
    """Order-of-magnitude thickness of a shock: δ ~ ν/(u₁ − u₂).

    Book: §15.6 "Normal shock structure":  (u₁ − u₂)δ/ν ~ 1 (a Reynolds number of order one on the velocity jump).
    Parameters: du velocity jump u₁ − u₂ [m/s]; nu kinematic viscosity [m²/s].   Returns δ [m].
    Compare with ``mean_free_path_air``: a few mean free paths — the continuum model is at its limit.
    Assumptions: an estimate, not a result.
    Validation: V1 — an order-of-magnitude estimate, compared with becker_shock_thickness — tests/test_ch15.py:
    test_becker_shock_profile_V1_satisfies_the_implicit_relation_and_the_ode.
    Label: analytic (the formula; as a thickness it is an order-of-magnitude estimate).
    """
    return _S(_F(nu) / _F(du))


def shock_structure_constants(M1: float, T1: float, p1: float, gas: PerfectGas | None = None) -> dict:
    """The three constants of integration of the shock-structure equations for a given upstream state.

    Book: §15.6 (unnumbered): ρu = m, mu + p = μ″u′ + mV, m(h + ½u²) = μ″uu′ + kT′ + mI, evaluated far upstream
    where the gradients vanish:  m = ρ₁u₁,  V = u₁ + p₁/m,  I = C_pT₁ + ½u₁².
    Parameters: M1 [–] > 1; T1 [K]; p1 [Pa]; gas.
    Returns dict: ``m`` [kg/(m² s)], ``V`` [m/s], ``I`` [J/kg], ``u1``, ``u2`` [m/s], ``U1``, ``U2`` [–] (= u/V).
    Assumptions: perfect gas.
    Validation: V2, V1 — tests/test_ch15.py: test_shock_structure_V2_first_integrals_and_roots,
    test_becker_shock_profile_V1_satisfies_the_implicit_relation_and_the_ode.
    Label: symbolic, analytic.
    """
    g = GD._gas(gas)
    u1 = M1 * float(g.c(T1))
    rho1 = p1 / (g.R * T1)
    m = rho1 * u1
    V = u1 + p1 / m
    s = GD.normal_shock(M1, g.gamma)
    return dict(m=m, V=V, I=g.cp * T1 + 0.5 * u1 ** 2, u1=u1, u2=u1 * s["u_ratio"], U1=u1 / V,
                U2=u1 * s["u_ratio"] / V)


def becker_shock_roots(M1, gamma: float = 1.4):
    """Scaled speeds U = u/V far upstream and far downstream of a shock: U₁ = γM₁²/(γM₁² + 1), U₂ the same with M₂.

    Book: §15.6 (unnumbered): U₁,₂ = [γ/(γ + 1)]{1 ± [1 − 2(γ² − 1)I/(γ²V²)]^{1/2}}.  The closed forms in the Mach
    numbers are ours (V = u + p/(ρu) = u(γM² + 1)/(γM²) on either side); U₁ + U₂ = 2γ/(γ + 1).
    Parameters: M1 [–] ≥ 1; gamma [–].   Returns (U₁, U₂) [–].
    Assumptions: perfect gas.
    Validation: V2, V1 — tests/test_ch15.py: test_shock_structure_V2_first_integrals_and_roots,
    test_becker_shock_profile_V1_satisfies_the_implicit_relation_and_the_ode.
    Label: symbolic, analytic.
    """
    M1 = _F(M1)
    M2 = _F(GD.normal_shock_mach(M1, gamma))
    return (_S(gamma * M1 ** 2 / (gamma * M1 ** 2 + 1.0)), _S(gamma * M2 ** 2 / (gamma * M2 ** 2 + 1.0)))


def becker_shock_ode(eta, U, M1: float, gamma: float = 1.4):
    """Right-hand side dU/dη of the shock-structure equation for Pr″ = 1.

    Book: §15.6 (unnumbered): −U(U₁ − U)⁻¹(U − U₂)⁻¹ dU = [(γ + 1)/(2γ)](m/μ″)dx, with η ≡ m∫dx/μ″.
    Parameters: eta [–] (unused: the equation is autonomous; kept for ``solve_ivp``); U = u/V [–]; M1 [–]; gamma.
    Returns dU/dη [–] (negative between the plateaus: the gas slows down).
    Assumptions: perfect gas, Pr″ = μ″C_p/k = 1 (Becker), steady one-dimensional continuum flow.
    Validation: V1 — tests/test_ch15.py: test_becker_shock_profile_V1_satisfies_the_implicit_relation_and_the_ode.
    Label: analytic.
    """
    U1, U2 = becker_shock_roots(M1, gamma)
    U = _F(U)
    return _S(-(gamma + 1.0) / (2.0 * gamma) * (U1 - U) * (U - U2) / U)


def becker_shock_profile(eta, M1: float, gamma: float = 1.4, *, eta0: float = 0.0):
    """Velocity profile U(η) = u/V through a normal shock (Becker's exact solution for Pr″ = 1).

    Book: §15.6 (unnumbered implicit result); Fig. 15.11:
    [U₁ ln(U₁ − U) − U₂ ln(U − U₂)]/(U₁ − U₂) = [(γ + 1)/(2γ)] η,  η = m∫dx/μ″ (roughly a count of mean free paths).
    Parameters: eta stretched coordinate [–] (scalar or array; upstream is η → −∞); M1 [–] > 1; gamma [–];
    eta0 — where the profile passes through the mid-value ½(U₁ + U₂) (the additive constant of the integral; the
    origin of η is arbitrary — our choice).
    Returns U [–], between U₂ (downstream) and U₁ (upstream).
    Method: with U = U₂ + (U₁ − U₂)/(1 + e^z) the implicit relation becomes U₁z/(U₁ − U₂) − ln(1 + e^z) + ln 2 =
    [(γ + 1)/(2γ)](η − η₀), monotonic in z; solved by ``brentq`` per point with a bracket from the two asymptotes.
    Assumptions: as :func:`becker_shock_ode`; the book restricts the result to moderate Mach numbers.
    Validation: V1 — tests/test_ch15.py: test_becker_shock_profile_V1_satisfies_the_implicit_relation_and_the_ode.
    Label: analytic.
    """
    U1, U2 = becker_shock_roots(float(M1), gamma)
    dU = U1 - U2
    if dU <= 0.0:
        return _S(np.full_like(_F(eta), U1))
    kk = (gamma + 1.0) / (2.0 * gamma)

    def one(e: float) -> float:
        tgt = kk * (e - eta0) - np.log(2.0)
        fun = lambda z: U1 * z / dU - np.logaddexp(0.0, z) - tgt  # noqa: E731
        za, zb = tgt * dU / U1, tgt * dU / U2
        lo, hi = min(za, zb) - 2.0, max(za, zb) + 2.0
        z = brentq(fun, lo, hi, xtol=1e-13, rtol=1e-14)
        return U2 + dU / (1.0 + np.exp(z)) if z < 700.0 else U2

    return _S(np.vectorize(one, otypes=[float])(_F(eta)))


def becker_shock_thickness(M1, gamma: float = 1.4):
    """Maximum-slope thickness of the Becker profile in the stretched coordinate: Δη = (U₁ − U₂)/max|dU/dη| (ours).

    Book: §15.6, Fig. 15.11 (thickness read off the plot).  From the structure equation the slope is largest at
    U = √(U₁U₂), where (U₁ − U)(U − U₂)/U = (√U₁ − √U₂)², so Δη = [2γ/(γ + 1)](√U₁ + √U₂)/(√U₁ − √U₂).
    Parameters: M1 [–] > 1; gamma [–].   Returns Δη [–] (→ ∞ as M₁ → 1: weak shocks are thick).
    Assumptions: as :func:`becker_shock_ode`.
    Validation: V1 — tests/test_ch15.py: test_becker_shock_profile_V1_satisfies_the_implicit_relation_and_the_ode.
    Label: analytic.
    """
    U1, U2 = becker_shock_roots(M1, gamma)
    a, b = np.sqrt(_F(U1)), np.sqrt(_F(U2))
    with np.errstate(divide="ignore"):
        return _S(2.0 * gamma / (gamma + 1.0) * (a + b) / (a - b))


# ======================================================================================================================
# §15.9 — wedge, expansion fan, curved wall
# ======================================================================================================================
def wedge_flow(M1: float, half_angle: float, gamma: float = 1.4) -> dict:
    """Supersonic flow past a wedge (or into a compression corner): attached weak shock, or detached.

    Book: §15.9, Figs. 15.20, 15.21 (the attached shock is the weak root; beyond δ_max a curved detached shock stands
    ahead of the body — normal on the axis, strong near it, weak far out).
    Parameters: M1 [–] > 1; half_angle wedge half-angle = flow deflection δ [rad] ≥ 0; gamma [–].
    Returns dict: ``attached`` (bool), ``sigma`` [rad], ``M2`` [–], ``p_ratio`` (surface pressure p₂/p₁), ``T_ratio``,
    ``rho_ratio``, ``p0_ratio`` [–], ``delta_max`` [rad], ``M_detach`` (smallest M₁ with an attached shock at this
    angle), ``subsonic_behind`` (bool).  Detached: ``sigma``, ``M2`` … are NaN and ``axis_p_ratio`` is the
    normal-shock pressure ratio on the centre streamline (the only number the one-dimensional theory gives; the
    stand-off distance and the shock shape are not computed).
    Assumptions: two-dimensional wedge (not a cone), perfect gas, inviscid.
    Validation: V1 — attached case = the weak root of oblique_shock on both sides of δ_max; the detached case is a
    sketch (only the centre-line normal shock is verified — verification report, Open item 6) — tests/test_ch15.py:
    test_oblique_shock_V1_two_roots_maximum_deflection_and_detachment.
    Label: analytic (attached case); qualitative (detached case).
    """
    dm = float(GD.max_deflection(M1, gamma)["delta_max"])
    out = dict(delta_max=dm, M_detach=float(GD.detachment_mach(half_angle, gamma)),
               axis_p_ratio=float(GD.normal_shock_pressure_ratio(M1, gamma)))
    if half_angle > dm:
        out.update(attached=False, sigma=np.nan, M2=np.nan, p_ratio=np.nan, T_ratio=np.nan, rho_ratio=np.nan,
                   p0_ratio=np.nan, subsonic_behind=True)
        return out
    o = GD.oblique_shock(M1, delta=half_angle, branch="weak", gamma=gamma)
    out.update(attached=True, sigma=o["sigma"], M2=o["M2"], p_ratio=o["p_ratio"], T_ratio=o["T_ratio"],
               rho_ratio=o["rho_ratio"], p0_ratio=o["p0_ratio"], subsonic_behind=bool(o["M2"] < 1.0))
    return out


def expansion_fan(M1: float, turn: float, gamma: float = 1.4, n_rays: int = 7, *, r_start: float = 1.0) -> dict:
    """The centred Prandtl–Meyer fan at a sharp convex corner: its rays, the state on each and one streamline.

    Book: §15.9, Fig. 15.23 and Eq. (15.56) with δ + ν(M) = const.  Geometry (ours): the stream arrives along +x over
    a wall on y = 0, x < 0; at the origin the wall turns **down** (away from the flow) by ``turn``.  On the ray where
    the flow has turned by θ (0 … turn) the Mach number is ν⁻¹(ν(M₁) + θ), the flow direction is −θ and the ray makes
    the angle φ = μ(M) − θ with +x: first ray at μ₁ to the upstream wall, last at μ₂ to the downstream wall.
    A streamline crosses the ray of Mach number M at the radius r = r_start [(1 + ½(γ − 1)M²)/(1 + ½(γ − 1)M₁²)]^{
    (γ+1)/(2(γ−1))} (ours: from mass conservation between two rays, r ρ c = const, because the velocity component
    normal to a ray is the local sound speed).
    Parameters: M1 [–] ≥ 1; turn [rad] > 0 (expansion; the book's δ = −turn); gamma [–]; n_rays ≥ 2; r_start radius
    at which the streamline meets the first ray [m].
    Returns dict: ``theta`` (n,) [rad], ``M`` (n,), ``mu`` (n,) [rad], ``ray_angle`` (n,) [rad] (from +x,
    counter-clockwise), ``p_over_p1``, ``T_over_T1`` (n,), ``M2``, ``p2_over_p1``, ``streamline_r`` (n,) [m],
    ``streamline_x``, ``streamline_y`` (n,) [m], ``wall_angle2`` (= −turn) [rad].
    Assumptions: perfect gas; isentropic; uniform flow upstream; no wave arrives from above.
    Validation: V1, V4 — tests/test_ch15.py: test_expansion_fan_V4_mass_flux_between_rays_and_streamline_direction.
    Label: analytic, conserved.
    """
    if turn <= 0.0:
        raise ValueError("expansion_fan needs turn > 0 (a turn away from the flow)")
    th = np.linspace(0.0, turn, int(n_rays))
    s = GD.prandtl_meyer_turn(M1, th, gamma)
    M, mu = np.atleast_1d(s["M2"]), np.atleast_1d(s["mu2"])
    phi = mu - th
    t1 = 1.0 + 0.5 * (gamma - 1.0) * M1 ** 2
    r = r_start * ((1.0 + 0.5 * (gamma - 1.0) * M ** 2) / t1) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    return dict(theta=th, M=M, mu=mu, ray_angle=phi, p_over_p1=np.atleast_1d(s["p_ratio"]),
                T_over_T1=np.atleast_1d(s["T_ratio"]), M2=float(M[-1]), p2_over_p1=float(np.atleast_1d(s["p_ratio"])[-1]),
                streamline_r=r, streamline_x=r * np.cos(phi), streamline_y=r * np.sin(phi), wall_angle2=-float(turn))


def fan_staircase(M1: float, turn: float, N: int = 10, gamma: float = 1.4) -> dict:
    """An expansion built from N equal small turns, each a weak wave — the from-scratch route to the fan.

    Book: §15.9: dδ = −[√(M² − 1)/M] dM/(1 + ½(γ − 1)M²) (unnumbered, before (15.56)) for the Mach number and the
    weak-wave rule (15.55) (p₂ − p₁)/p₁ = γM²δ/√(M² − 1) for the pressure, applied step by step with δ = −turn/N
    (explicit Euler: each step uses the state before it).  ``p_over_p1`` is the isentropic pressure (15.29) at the
    staircase Mach number with p₀ constant (the fan is isentropic); ``p_over_p1_weak`` is the running product of the
    weak-wave steps (15.55).
    Parameters: M1 [–] > 1; turn total expansion angle [rad] (> 0; a negative value is a staircase of weak
    compressions); N number of steps; gamma [–].
    Returns dict: ``theta`` (N + 1,) cumulative turn [rad], ``M`` (N + 1,), ``p_over_p1`` (N + 1,),
    ``p_over_p1_weak`` (N + 1,), ``M2_exact``, ``p2_exact`` (from the exact function), ``M_error``, ``p_error``,
    ``p_weak_error`` (relative errors of the last step; they fall like 1/N).
    Assumptions: as (15.55), (15.56).
    Validation: V4, V3 — first-order convergence to prandtl_meyer_turn (observed order 1.00 in M and p, N = 20 … 320) —
    tests/test_ch15.py: test_expansion_fan_V4_mass_flux_between_rays_and_streamline_direction.
    Label: conserved, converged.
    """
    if M1 <= 1.0:
        raise ValueError("fan_staircase needs M1 > 1 (the differential relation is singular at M = 1)")
    N = int(N)
    dth = turn / N
    M, p = np.empty(N + 1), np.empty(N + 1)
    M[0], p[0] = M1, 1.0
    for i in range(N):
        b = np.sqrt(M[i] ** 2 - 1.0)
        M[i + 1] = M[i] + dth * M[i] * (1.0 + 0.5 * (gamma - 1.0) * M[i] ** 2) / b  # dδ = −dθ, text before (15.56)
        p[i + 1] = p[i] * (1.0 - gamma * M[i] ** 2 * dth / b)  # Eq. (15.55) with δ = −dθ
    ex = GD.prandtl_meyer_turn(M1, turn, gamma)
    p_isen = _F(GD.stagnation_pressure_ratio(M1, gamma)) / _F(GD.stagnation_pressure_ratio(M, gamma))  # Eq. (15.29)
    return dict(theta=np.linspace(0.0, turn, N + 1), M=M, p_over_p1=p_isen, p_over_p1_weak=p, M2_exact=ex["M2"],
                p2_exact=ex["p_ratio"], M_error=float(M[-1] / ex["M2"] - 1.0),
                p_error=float(p_isen[-1] / ex["p_ratio"] - 1.0), p_weak_error=float(p[-1] / ex["p_ratio"] - 1.0))


def curved_wall_mach_lines(M1: float, wall_angle_fn: Callable, x, gamma: float = 1.4) -> dict:
    """Mach lines sent out by a smoothly curved wall: they converge on a compression and diverge on an expansion.

    Book: §15.9, Fig. 15.22 (simple wave: flow direction and Mach number are constant along each Mach line when no
    wave arrives from the other side).  The stream arrives along +x at M₁ above the wall; ``wall_angle_fn(x)`` is the
    wall inclination θ_w [rad], **positive when the wall turns into the flow** (the book's δ > 0, compression).  At
    each station ν(M) = ν(M₁) − θ_w, and the Mach line leaves the wall at the angle θ_w + μ(M) from +x.
    Parameters: M1 [–] > 1; wall_angle_fn callable x [m] → θ_w [rad] (θ_w = 0 upstream); x stations [m] (increasing,
    n ≥ 2); gamma [–].
    Returns dict: ``x_wall``, ``y_wall`` (n,) [m] (y by trapezoidal integration of tan θ_w), ``theta`` (n,) [rad],
    ``M`` (n,), ``mu`` (n,) [rad], ``line_angle`` (n,) [rad], ``crossings`` (k, 2) [m] (intersections of neighbouring
    lines above the wall), ``first_crossing`` ((x, y) nearest the wall, or ``None`` — an expansion: the lines
    diverge), ``compression`` (bool).
    Assumptions: isentropic simple wave up to the first crossing; beyond it a shock forms (not computed).
    Validation: V1 — tests/test_ch15.py: test_curved_wall_V1_mach_lines_first_cross_on_the_envelope.
    Label: analytic.
    """
    x = np.atleast_1d(_F(x))
    th = np.atleast_1d(_F(wall_angle_fn(x)))
    s = GD.prandtl_meyer_turn(M1, -th, gamma)
    M, mu = np.atleast_1d(s["M2"]), np.atleast_1d(s["mu2"])
    yw = np.concatenate([[0.0], np.cumsum(0.5 * (np.tan(th[1:]) + np.tan(th[:-1])) * np.diff(x))])
    ang = th + mu
    cx, cy, dist = [], [], []
    for i in range(x.size - 1):
        d1 = np.array([np.cos(ang[i]), np.sin(ang[i])])
        d2 = np.array([np.cos(ang[i + 1]), np.sin(ang[i + 1])])
        det = d1[0] * (-d2[1]) + d2[0] * d1[1]
        if abs(det) < 1e-14:
            continue
        rhs = np.array([x[i + 1] - x[i], yw[i + 1] - yw[i]])
        a = (rhs[0] * (-d2[1]) + d2[0] * rhs[1]) / det
        b = (d1[0] * rhs[1] - d1[1] * rhs[0]) / det
        if a > 0.0 and b > 0.0:
            cx.append(x[i] + a * d1[0])
            cy.append(yw[i] + a * d1[1])
            dist.append(a)
    cr = np.column_stack([cx, cy]) if cx else np.empty((0, 2))
    first = tuple(cr[int(np.argmin(dist))]) if cx else None
    return dict(x_wall=x, y_wall=yw, theta=th, M=M, mu=mu, line_angle=ang, crossings=cr, first_crossing=first,
                compression=bool(cx))


# ======================================================================================================================
# §15.10 — thin airfoil theory in supersonic flow
# ======================================================================================================================
def _beta(M_inf):
    M = _F(M_inf)
    with np.errstate(invalid="ignore"):
        return np.sqrt(np.where(M > 1.0, M ** 2 - 1.0, np.nan))


def ackeret_pressure(delta, M_inf, gamma: float = 1.4):
    """Pressure on a surface inclined at a small angle δ to a supersonic stream: (p − p_∞)/p_∞ = γM_∞²δ/√(M_∞² − 1).

    Book: §15.10, Eq. (15.57) (Ackeret; (15.55) with the free-stream values).  δ > 0 where the surface turns the
    flow into itself (compression), δ < 0 where it turns away.
    Parameters: delta [rad]; M_inf [–] > 1; gamma [–].   Returns (p − p_∞)/p_∞ [–] (NaN for M_∞ ≤ 1).
    Assumptions: thin body, small angles, M_∞ not close to 1 and not hypersonic.
    Validation: V1 — tests/test_ch15.py:
    test_weak_wave_V1_exact_shock_and_expansion_approach_the_linear_rule_from_both_sides.
    Label: analytic.
    """
    return _S(gamma * _F(M_inf) ** 2 * _F(delta) / _beta(M_inf))  # Eq. (15.57)


def ackeret_cp(delta, M_inf):
    """Pressure coefficient of linearised supersonic flow, C_p = (p − p_∞)/(½ρ_∞U_∞²) = 2δ/√(M_∞² − 1) (ours).

    Book: §15.10, Eq. (15.57) divided by ½ρ_∞U_∞²/p_∞ = ½γM_∞² (the book's "ρU² = γpM²").
    Parameters: delta [rad]; M_inf [–] > 1.   Returns C_p [–] (independent of γ; NaN for M_∞ ≤ 1).
    Assumptions: as :func:`ackeret_pressure`.
    Validation: V1 — tests/test_ch15.py: test_flat_plate_V1_linear_theory_formulas.
    Label: analytic.
    """
    return _S(2.0 * _F(delta) / _beta(M_inf))


def ackeret_flat_plate(alpha, M_inf) -> dict:
    """Lift and wave-drag coefficients of a flat plate at a small angle of attack in supersonic flow.

    Book: §15.10, Eq. (15.59)  C_L ≅ 4α/√(M_∞² − 1),  C_D ≅ 4α²/√(M_∞² − 1) (per unit span, on the chord b).
    Parameters: alpha [rad]; M_inf [–] > 1.
    Returns dict [–]: ``CL``, ``CD`` (wave drag — exists in inviscid flow), ``L_over_D`` (= 1/α), ``lift_slope``
    (4/√(M_∞² − 1) per radian, against 2π in incompressible flow).  NaN for M_∞ ≤ 1.
    Assumptions: linearised (Ackeret) theory; fails as M_∞ → 1.
    Validation: V1 — the formulas of (15.59) and L/D = 1/α; that the linear plate is the small-α limit of
    shock–expansion (error second order in α) is shown through flat_plate_supersonic in
    test_flat_plate_V7_linear_theory_is_the_small_angle_limit_of_shock_expansion — tests/test_ch15.py:
    test_flat_plate_V1_linear_theory_formulas.
    Label: analytic.
    """
    a, b = _F(alpha), _beta(M_inf)
    with np.errstate(divide="ignore", invalid="ignore"):
        return dict(CL=_S(4.0 * a / b), CD=_S(4.0 * a ** 2 / b), L_over_D=_S(1.0 / a + 0.0 * b),  # Eq. (15.59)
                    lift_slope=_S(4.0 / b + 0.0 * a))


def flat_plate_supersonic(alpha: float, M_inf: float, gamma: float = 1.4, *, method: str) -> dict:
    """Pressures and forces on a flat plate at incidence in a supersonic stream — linear theory or exact
    shock–expansion.

    Book: §15.10, Fig. 15.24, the plate pressures (p₂ − p_∞)/p_∞ = −γM_∞²α/√(M_∞² − 1) (upper, expansion),
    (p₃ − p_∞)/p_∞ = +γM_∞²α/√(M_∞² − 1) (lower, compression), Eq. (15.58) L = (p₃ − p₂)b cos α, D = (p₃ − p₂)b sin α
    and Eq. (15.59).
    Parameters: alpha angle of attack [rad] (nose up positive); M_inf [–] > 1; gamma [–];
    method — a **required keyword** — "linear" (the book's: (15.57) on both sides, cos α ≈ 1, sin α ≈ α) or "shock_expansion" (ours: the exact
    Prandtl–Meyer expansion on the suction side and the exact weak oblique shock on the pressure side, exact
    resolution of the normal force).
    Returns dict [–]: ``p_upper`` (p₂/p_∞), ``p_lower`` (p₃/p_∞), ``CN`` (normal-force coefficient
    (p₃ − p₂)/(½γM_∞²p_∞)), ``CL``, ``CD``, ``L_over_pb`` (L/(p_∞b)), ``D_over_pb``, ``method``.
    ``shock_expansion`` raises ``ValueError("detached …")`` when |α| exceeds δ_max(M_∞).
    Assumptions: two-dimensional, inviscid, perfect gas; waves from the two sides do not interact with the plate.
    Validation: V1, V7 — the linear result is the small-α limit of shock–expansion with an error of **second** order in
    α for C_L and C_D (measured 2.02, 2.00, 1.99 at M_∞ = 1.6, 2, 3: the second-order (Busemann) pressure terms are
    equal on the two surfaces and cancel in p₃ − p₂); a single surface pressure is first order (measured 1.02); C_D/C_L
    = α (linear), tan α (exact) — tests/test_ch15.py: test_flat_plate_V1_linear_theory_formulas,
    test_flat_plate_V7_linear_theory_is_the_small_angle_limit_of_shock_expansion.
    Label: analytic.
    """
    a = abs(float(alpha))
    sgn = 1.0 if alpha >= 0 else -1.0
    q = 0.5 * gamma * M_inf ** 2  # ½ρ_∞U_∞²/p_∞
    if method == "linear":
        dp = float(ackeret_pressure(a, M_inf, gamma))
        p_suction, p_pressure = 1.0 - dp, 1.0 + dp
        CN = (p_pressure - p_suction) / q
        CL, CD = CN, CN * a  # Eq. (15.58): cos α ≈ 1, sin α ≈ α
    elif method == "shock_expansion":
        if a == 0.0:
            p_suction = p_pressure = 1.0
        else:
            p_suction = float(GD.prandtl_meyer_turn(M_inf, a, gamma)["p_ratio"])
            p_pressure = float(GD.oblique_shock(M_inf, delta=a, branch="weak", gamma=gamma)["p_ratio"])
        CN = (p_pressure - p_suction) / q
        CL, CD = CN * np.cos(a), CN * np.sin(a)  # Eq. (15.58), exact resolution
    else:
        raise ValueError('method must be "linear" or "shock_expansion"')
    p_up, p_lo = (p_suction, p_pressure) if sgn > 0 else (p_pressure, p_suction)
    return dict(p_upper=p_up, p_lower=p_lo, CN=sgn * CN, CL=sgn * CL, CD=float(CD), L_over_pb=sgn * CL * q,
                D_over_pb=float(CD) * q, method=method)


def ackeret_airfoil(yu_fn: Callable, yl_fn: Callable, alpha: float, M_inf: float, n: int = 2001) -> dict:
    """Linearised supersonic flow past a thin section of any shape (ours; Exercise 15.17).

    Book: §15.10, Eq. (15.57) applied to each surface: with the chord on 0 ≤ x ≤ 1 and the section tilted nose-up by
    α, the upper surface is inclined to the stream by dy_u/dx − α and the lower by α − dy_l/dx (each positive when it
    turns the flow into itself), so C_p,u = 2(y_u′ − α)/√(M_∞² − 1), C_p,l = 2(α − y_l′)/√(M_∞² − 1);
    C_L = ∫(C_p,l − C_p,u)dx = 4α/√(M_∞² − 1) for a closed section (independent of thickness and camber),
    C_D = ∫[C_p,u(y_u′ − α) + C_p,l(α − y_l′)]dx = [4/√(M_∞² − 1)][α² + mean(y_c′²) + mean(y_t′²)], with the camber
    line y_c = ½(y_u + y_l) and the half-thickness y_t = ½(y_u − y_l).
    Parameters: yu_fn, yl_fn callables x/chord ∈ [0, 1] → y/chord of the upper and lower surface; alpha [rad];
    M_inf [–] > 1; n number of chordwise points (odd; default 2001).
    Returns dict: ``x`` (n,), ``cp_upper``, ``cp_lower`` (n,), ``CL``, ``CD``, ``CD_incidence`` (4α²/β),
    ``CD_camber``, ``CD_thickness``, ``CD_cross`` (−8α·mean(y_c′)/β: zero when the chord joins the two edges).
    Method: slopes by explicit second-order differences (one-sided at the edges); composite Simpson quadrature.
    Assumptions: thin section, small angles, sharp edges, attached waves; inviscid.
    Validation: V1, V3 — flat plate → (15.59); double wedge 4τ²/β and half-sine π²τ²/2β closed forms; thickness integral
    converges at order 1.97; C_L independent of shape — tests/test_ch15.py:
    test_ackeret_airfoil_V1_thickness_camber_and_incidence_drag.
    Label: analytic, converged.
    """
    from scipy.integrate import simpson

    n = int(n) | 1
    x = np.linspace(0.0, 1.0, n)
    yu, yl = _F(yu_fn(x)) + 0.0 * x, _F(yl_fn(x)) + 0.0 * x
    su, sl = AC._ddx(yu, x), AC._ddx(yl, x)  # surface slopes dy/dx
    b = float(_beta(M_inf))
    cpu, cpl = 2.0 * (su - alpha) / b, 2.0 * (alpha - sl) / b  # Eq. (15.57) as a pressure coefficient
    CL = float(simpson(cpl - cpu, x=x))
    CD = float(simpson(cpu * (su - alpha) + cpl * (alpha - sl), x=x))
    sc, st = 0.5 * (su + sl), 0.5 * (su - sl)  # slopes of the camber line and of the half-thickness
    return dict(x=x, cp_upper=cpu, cp_lower=cpl, CL=CL, CD=CD, CD_incidence=4.0 * alpha ** 2 / b,
                CD_camber=float(4.0 / b * simpson(sc ** 2, x=x)), CD_thickness=float(4.0 / b * simpson(st ** 2, x=x)),
                CD_cross=float(-8.0 * alpha / b * simpson(sc, x=x)))


def flat_plate_wave_pattern(alpha: float, M_inf: float, gamma: float = 1.4) -> dict:
    """The four waves of a flat plate at incidence in a supersonic stream and the states between them.

    Book: §15.10, Fig. 15.24: at the leading edge an expansion fan above and an oblique shock below; at the trailing
    edge the opposite waves bring both streams back (nearly) to the free-stream direction.
    Geometry (ours): free stream along +x; plate from the origin to (cos α, −sin α) (nose up, α > 0); angles from +x,
    counter-clockwise positive.
    Parameters: alpha [rad] > 0; M_inf [–] > 1; gamma [–].
    Returns dict: ``le_fan_first``, ``le_fan_last`` [rad] (rays of the upper leading-edge fan: μ_∞ and μ₂ − α);
    ``le_shock`` [rad] (lower leading-edge shock, = −σ); ``te_shock`` [rad] (upper trailing-edge shock, = σ_te − α);
    ``te_fan_first``, ``te_fan_last`` [rad] (lower trailing-edge fan: −α − μ₃ and −μ₄); ``M2`` (above the plate),
    ``M3`` (below), ``M4_upper``, ``M4_lower`` (behind the trailing-edge waves); ``p2``, ``p3``, ``p4_upper``,
    ``p4_lower`` (over p_∞); ``plate_angle`` (= −α).
    The trailing-edge waves are taken to turn each stream by exactly α (back to the free-stream direction, as the
    book says); the two exit pressures then differ slightly — in the exact solution a slip line leaves the trailing
    edge at a small angle, which is not computed.
    Assumptions: two-dimensional, inviscid, perfect gas; attached waves (α < δ_max).
    Validation: V7 — leading-edge waves and plate pressures against flat_plate_supersonic; the trailing-edge closure is
    a sketch (p₄ returns to p_∞ to O(α³) — verification report, Open item 6) — tests/test_ch15.py:
    test_flat_plate_V7_linear_theory_is_the_small_angle_limit_of_shock_expansion.
    Label: analytic (leading-edge waves, pressures); qualitative (trailing-edge closure).
    """
    if alpha <= 0.0:
        raise ValueError("flat_plate_wave_pattern expects alpha > 0 (nose up)")
    up = GD.prandtl_meyer_turn(M_inf, alpha, gamma)
    lo = GD.oblique_shock(M_inf, delta=alpha, branch="weak", gamma=gamma)
    te_u = GD.oblique_shock(up["M2"], delta=alpha, branch="weak", gamma=gamma)
    te_l = GD.prandtl_meyer_turn(lo["M2"], alpha, gamma)
    return dict(le_fan_first=float(GD.mach_angle(M_inf)), le_fan_last=up["mu2"] - alpha, le_shock=-lo["sigma"],
                te_shock=te_u["sigma"] - alpha, te_fan_first=-alpha - te_l["mu1"], te_fan_last=-te_l["mu2"],
                M2=up["M2"], M3=lo["M2"], M4_upper=te_u["M2"], M4_lower=te_l["M2"], p2=up["p_ratio"],
                p3=lo["p_ratio"], p4_upper=up["p_ratio"] * te_u["p_ratio"], p4_lower=lo["p_ratio"] * te_l["p_ratio"],
                plate_angle=-float(alpha))


def lift_slope_vs_mach(M, *, model: str = "book", band: tuple[float, float] = (0.8, 1.2)):
    """Lift-curve slope dC_L/dα of a thin flat section against Mach number: 2π below the speed of sound, 4/√(M² − 1)
    above, and nothing in between.

    Book: §15.10 (text after (15.59): compare C_L ≅ 4α/√(M_∞² − 1) with the incompressible C_L = 2πα of the preceding
    chapter (14.12); "these expressions do not hold at transonic speeds").
    Parameters: M [–] ≥ 0; model "book" (2π for every M below the band) or "prandtl_glauert" (2π/√(1 − M²) below the
    band — **not in the book**, the usual subsonic compressibility correction, shown for context only); band —
    (lower, upper) Mach numbers between which NaN is returned (our choice, the transonic band of
    :func:`mach_regime`).
    Returns dC_L/dα [1/rad].
    Assumptions: thin flat plate, small α, inviscid; linear theories on both sides.
    Validation: V1 — tests/test_ch15.py: test_flat_plate_V1_linear_theory_formulas.
    Label: analytic.
    """
    M = _F(M)
    with np.errstate(invalid="ignore", divide="ignore"):
        if model == "book":
            sub = np.full_like(M, 2.0 * np.pi)
        elif model == "prandtl_glauert":
            sub = 2.0 * np.pi / np.sqrt(np.where(M < 1.0, 1.0 - M ** 2, np.nan))
        else:
            raise ValueError('model must be "book" or "prandtl_glauert"')
        sup = 4.0 / np.sqrt(np.where(M > 1.0, M ** 2 - 1.0, np.nan))  # Eq. (15.59)
        return _S(np.where(M < band[0], sub, np.where(M > band[1], sup, np.nan)))


# ======================================================================================================================
# optional (ours): one-dimensional unsteady Euler equations
# ======================================================================================================================
def euler1d_run(rho0, u0, p0, x, t_end: float, gamma: float = 1.4, *, cfl: float = 0.5, eps_av: float = 2.0,
                bc: str = "periodic") -> dict:
    """March the one-dimensional Euler equations with MacCormack's scheme (**not in the book — ours**; optional
    numerical companion of the steepening argument of §15.2 and the moving shock of §15.6).

    ∂/∂t (ρ, ρu, ρE) + ∂/∂x (ρu, ρu² + p, u(ρE + p)) = 0,  ρE = p/(γ − 1) + ½ρu².
    Book: none for the scheme here; the two-step scheme is ch10's (10.100)–(10.102)
    (``core.maccormack.maccormack_step``), alternating forward/backward and backward/forward predictor–corrector.
    DEVIATION: the scheme and the artificial viscosity are ours — the book solves no unsteady compressible problem
    in this chapter.  After each step a conservative pressure-switched smoothing is added,
    U_i += d_{i+½} − d_{i−½}, d_{i+½} = min(0.45, eps_av·max(s_i, s_{i+1}))(U_{i+1} − U_i),
    s_i = |p_{i+1} − 2p_i + p_{i−1}|/(p_{i+1} + 2p_i + p_{i−1}); it is third order small where the flow is smooth
    and first order at a captured shock (which is then a few cells thick).  ``eps_av=0`` switches it off.
    Parameters: rho0 [kg/m³], u0 [m/s], p0 [Pa] initial fields at the nodes x (n,) [m] (uniform); t_end [s];
    gamma [–]; cfl Courant number on |u| + c (≤ 1; default 0.5); eps_av smoothing coefficient [–]; bc "periodic" or
    "transmissive" (zero-gradient ends: valid until a wave reaches an end).
    Returns dict: ``x``, ``rho``, ``u``, ``p`` at ``t`` (= t_end), ``n_steps``, ``mass``, ``momentum``, ``energy``
    (Σ over nodes × Δx at the end), ``mass0``, ``momentum0``, ``energy0`` (at the start).
    Cost: about 1 s for 800 nodes to a time of one domain crossing.
    Assumptions: perfect gas, inviscid; uniform grid.
    Validation: V4, V3 — tests/test_ch15.py:
    test_euler1d_V3_second_order_for_a_smooth_wave_and_first_order_at_a_captured_shock.
    Label: conserved, converged.
    """
    from .core.maccormack import maccormack_step

    x = _F(x)
    dx = float(x[1] - x[0])
    if not np.allclose(np.diff(x), dx, rtol=1e-9, atol=0.0):
        raise ValueError("euler1d_run needs a uniform grid")
    if cfl > 1.0:
        raise ValueError("cfl must be <= 1")
    if bc not in ("periodic", "transmissive"):
        raise ValueError('bc must be "periodic" or "transmissive"')
    rho = _F(rho0) + 0.0 * x
    u = _F(u0) + 0.0 * x
    p = _F(p0) + 0.0 * x
    U = np.array([rho, rho * u, p / (gamma - 1.0) + 0.5 * rho * u ** 2])

    def prim(W):
        r = W[0]
        v = W[1] / r
        return r, v, (gamma - 1.0) * (W[2] - 0.5 * r * v ** 2)

    def flux(W):
        r, v, pr = prim(W)
        return np.array([W[1], W[1] * v + pr, v * (W[2] + pr)])

    tot0 = U.sum(axis=1) * dx
    t, n = 0.0, 0
    while t < t_end * (1.0 - 1e-14):
        r, v, pr = prim(U)
        if np.any(pr <= 0.0) or np.any(r <= 0.0):
            raise FloatingPointError("euler1d_run: negative pressure or density (reduce cfl or raise eps_av)")
        dt = min(cfl * dx / float(np.max(np.abs(v) + np.sqrt(gamma * pr / r))), t_end - t)
        U = maccormack_step(U, flux, None, dt, dx, arrangement="F/B" if n % 2 == 0 else "B/F")
        if bc == "transmissive":  # the periodic step contaminates only the two end nodes: reset them first
            U[:, 0], U[:, -1] = U[:, 1], U[:, -2]
        # DEVIATION: pressure-switched artificial viscosity (ours) — the book gives no scheme for this problem;
        # MacCormack's scheme alone oscillates at a captured shock.
        if eps_av > 0.0:
            pr = prim(U)[2]
            s = np.abs(np.roll(pr, -1) - 2.0 * pr + np.roll(pr, 1)) / (np.roll(pr, -1) + 2.0 * pr + np.roll(pr, 1))
            nu_ = np.minimum(0.45, eps_av * np.maximum(s, np.roll(s, -1)))
            d = nu_[None, :] * (np.roll(U, -1, axis=1) - U)
            if bc == "transmissive":
                d[:, -1] = 0.0
            U = U + d - np.roll(d, 1, axis=1)
        if bc == "transmissive":
            U[:, 0], U[:, -1] = U[:, 1], U[:, -2]
        t += dt
        n += 1
    r, v, pr = prim(U)
    tot = U.sum(axis=1) * dx
    return dict(x=x, rho=r, u=v, p=pr, t=t, n_steps=n, mass=float(tot[0]), momentum=float(tot[1]), energy=float(tot[2]),
                mass0=float(tot0[0]), momentum0=float(tot0[1]), energy0=float(tot0[2]))
