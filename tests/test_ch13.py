"""Verification suite for Chapter 13 — Geophysical Fluid Dynamics (Kundu, Cohen & Dowling 5e, §§13.1–13.18).

Evidence levels (``verify-implementation`` skill): V1 analytic · V2 symbolic (sympy / dimensions) · V3 convergence ·
V4 conservation / invariant / integral identity · V5 published benchmark (``reference/ch13/benchmarks.json``, see
SOURCES.md) · V6 book value (private, ``tests/book_values_ch13.json``, skipped when absent) · V7 limits / symmetry /
invariance.  Every test name is ``test_<concept>_<level>_<what>``; the comment on the ``def`` line repeats the level and
the curation ID.

A items (CORE, >= 2 independent levels): C01 thin-shell equations (13.2)–(13.10) · C02 geostrophic balance
(13.11)–(13.13) · C03 thermal wind and Taylor–Proudman (13.15), (13.21) · C04 surface Ekman layer (13.22)–(13.29) ·
C05 Ekman transport and pumping (13.30) · C06 bottom Ekman layer (13.41) · C07 shallow-water set (13.44)–(13.45) ·
C08 vertical normal modes (13.56), (13.64), (13.65), (13.69) · C09 the cubic (13.76) · C10 Poincaré waves (13.80),
(13.82) · C11 Kelvin wave (13.87) · C12 Rossby radius and geostrophic adjustment · C13 potential vorticity (13.94) ·
C14 inertia–gravity waves (13.96), (13.99), (13.112) · C15 Rossby waves (13.117)–(13.120), Rayleigh–Kuo (13.124) ·
C16 Eady problem (13.136), (13.141) · C17 two-dimensional turbulence (13.143)–(13.145).

Derivations: every two- and three-star D row (D01 D02 D04–D27 D29) is re-derived in a ``test_*_V2_derivation`` test
**independently of the chapter's own sympy engines** (which are tested separately).

Printed slips kept as wrong variants that must FAIL the independent check (each page was re-read as an image first):
#1 the sign of the pressure gradient in (13.2) · #2 the friction force without 1/rho · #3 the sign in (13.17) ·
#4 the far-field direction of (13.26) · #6 "the first root at NH/c = 1" · #7 "negligible for omega >> f".
``ekman_surface_printed`` is true as printed for f > 0 (a trap, not a slip) and must refuse f <= 0.

All inputs are ours (``ch13.illustrative_inputs()`` or visibly odd numbers); nothing here is quoted from the book.

Run: ``.venv/Scripts/python.exe -m pytest tests/test_ch13.py -q -p no:cacheprovider``  (``-m "not slow"`` skips the
script runs and the cache regeneration).
"""
from __future__ import annotations

import inspect
import json
import math
import re
import subprocess
import sys
import warnings
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest
import sympy as sp
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq, minimize_scalar

from fluidpy import ch13_geophysical_fluid_dynamics as ch13
from fluidpy.core import gfd as GFD
from fluidpy.core import rotating as ROT
from fluidpy.core import shallow_water as SW
from fluidpy.core import stability as ST
from fluidpy.core import stratification as STRAT
from fluidpy.core import vertical_modes as VM
from fluidpy.core import waves as WAV
from fluidpy.core.units import Q_
from tools.convergence import observed_order

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference" / "ch13"
BOOK = Path(__file__).resolve().parent / "book_values_ch13.json"
book_only = pytest.mark.skipif(not BOOK.exists(), reason="book values are private; see CLAUDE.md rule 9")
slow = pytest.mark.slow
G0 = ch13.G0
PI = math.pi
I = ch13.illustrative_inputs()
LAT = I["lat"]                                   # 35 degrees north, in radians (ours)
F_N = float(ch13.coriolis_parameter(LAT))        # f > 0
F_S = -F_N                                       # the same latitude south
HEMI = (F_N, F_S)


def book():
    return json.loads(BOOK.read_text(encoding="utf-8"))


def rel(a, b):
    return abs(float(a) - float(b)) / max(abs(float(b)), 1e-300)


def bench():
    return json.loads((REF / "benchmarks.json").read_text(encoding="utf-8"))


def close(a, b, rtol=1e-5, atol=0.0):
    """``numpy.allclose`` with NO absolute tolerance unless one is given: numpy's default atol = 1e-8 would pass any two
    numbers of the size of f (1e-4 1/s) or beta (1e-11 1/(m s)) — a mutant with sin for cos in beta survived it."""
    return bool(np.allclose(a, b, rtol=rtol, atol=atol))


# =====================================================================================================================
# 0. Contract: every name of design Part C exists, is re-exported, has the documented shape, and is used in this file
# =====================================================================================================================
PART_C = """
OMEGA_EARTH OMEGA_SOLAR_DAY EARTH_RADIUS_MEAN SIDEREAL_DAY coriolis_parameter hemisphere earth_rotation_local
inertial_period f_plane beta_parameter beta_plane beta_plane_error vertical_velocity_scale coriolis_acceleration_local
thin_layer_terms eddy_friction scale_height buoyancy_frequency_sq geostrophic_velocity geostrophic_from_field
geostrophic_from_height geostrophic_streamfunction rossby_number ekman_number thermal_wind_shear
thermal_wind_from_temperature thermal_wind_integrate taylor_proudman_residual parcel_adjust ekman_depth ekman_surface
ekman_transport ekman_transport_partial ekman_residual ekman_vorticity_balance ekman_solve ekman_pumping
ekman_pumping_bottom ekman_finite_depth ekman_bottom ekman_bottom_transport ekman_force_balance
eddy_viscosity_from_depth long_wave_speed equivalent_depth baroclinic_mode_speed shallow_water_omega
shallow_water_branches shallow_water_discriminant shallow_water_regime dispersion_term_sizes poincare_omega
poincare_group_velocity poincare_amplitudes poincare_fields poincare_orbit inertial_oscillation inertial_radius
kelvin_wave kelvin_omega kelvin_decay_side kelvin_residuals rossby_radius rossby_radius_internal
rossby_radius_two_layer geostrophic_adjustment_1d adjustment_energy inertia_gravity_m2 inertia_gravity_band
inertia_gravity_omega inertia_gravity_regime inertia_gravity_group_velocity wkb_vertical_structure
inertia_gravity_fields inertia_gravity_hodograph lee_wave_m rossby_omega rossby_phase_speed rossby_group_velocity
rossby_max_frequency rossby_omega_circle rossby_long_wave_speed stationary_rossby_wavelength rossby_packet_spectrum
potential_vorticity step_vorticity absolute_vorticity_gradient rayleigh_kuo_criterion eady_alpha eady_factors
eady_phase_speed eady_growth_rate eady_critical eady_fastest eady_max_growth_rate eady_time_scale eady_mode
eady_vertical_velocity eady_fluxes fjortoft_transfer enstrophy_spectrum two_d_cascade_spectrum rhines_length
Modes vertical_modes vertical_modes_shooting modes_uniform_N rigid_lid_error w_structure rho_structure project
reconstruct orthogonality_matrix modal_amplitudes wkb_mode_speed
ShallowWater gaussian_bump geostrophic_state continuity_tendency momentum_tendencies relative_vorticity
sw_potential_vorticity energy dt_limit step run advect_particles linear_1d_step linear_1d_run qg_linear_evolve
qg_linear_evolve_1d barotropic_vorticity_rhs barotropic_run barotropic_invariants spectral_centroids
illustrative_inputs conventions_table book_slips lapse_rate_table wind_from_to ocean_density_gradient_budget
atmosphere_layers ocean_profile_idealized thermocline_N2 anisotropic_eddy_stress term_table_thin_layer
pressure_centre jet_section viscous_layer_thicknesses ekman_surface_printed coastal_upwelling flow_over_step
vertical_structure_solve wkb_error w_equation_rotating_residual lee_wave_field basin_crossing_time
rayleigh_kuo_eigs eady_basic_state eady_matrix eady_numeric_eigs load_reference_run
boussinesq_rotating_sympy perturbation_form_sympy eddy_friction_force_sympy thin_layer_sympy taylor_proudman_sympy
hydrostatic_linear_set_sympy v_equation_sympy rotating_internal_wave_set_sympy w_equation_rotating_sympy
pv_conservation_sympy qg_vorticity_sympy eady_qg_sympy
fig_ekman_surface fig_ekman_bottom fig_mode_roots fig_vertical_modes fig_poincare_kelvin_dispersion
fig_kelvin_sections fig_inertia_gravity_orbit fig_rossby_dispersion
""".split()
# Names beyond the contract that the scripts and caches use; each is exercised below as well.
EXTRA = """coriolis_parameter_deg thin_layer_residuals ekman_pumping_from_curl sverdrup_transport eady_wavelengths
make_state kelvin_state volume interpolate particle_potential_vorticity spectral_wavenumbers barotropic_velocity
barotropic_time_step barotropic_spectrum zonal_energy_fraction random_vorticity traps reexport_audit
friction_force_dimensions sympy_engines reference_run write_reference_runs explainer_constants""".split()


def test_contract_V7_every_part_c_name_exists_and_is_exercised_here():  # contract self-audit (ch12 lesson 5)
    src = Path(__file__).read_text(encoding="utf-8")
    body = src.split("def test_contract_V7_every_part_c_name_exists_and_is_exercised_here", 1)[1]
    missing = [n for n in PART_C + EXTRA if not hasattr(ch13, n)]
    assert missing == []
    unused = [n for n in PART_C + EXTRA if not re.search(r"\b(?:ch13|GFD|VM|SW)\." + re.escape(n) + r"\b", body)]
    assert unused == [], unused
    assert len(PART_C) == len(set(PART_C)) and len(PART_C) >= 175


def test_contract_V7_reexports_are_the_same_objects():  # C.4
    a = ch13.reexport_audit()
    assert a["missing"] == [] and a["aliased"] == {"SW.potential_vorticity": "sw_potential_vorticity"}
    assert a["same_object"] + len(a["aliased"]) == a["total"]
    assert ch13.potential_vorticity is GFD.potential_vorticity and ch13.sw_potential_vorticity is SW.potential_vorticity
    assert ch13.GFD is GFD and ch13.VM is VM and ch13.SW is SW
    assert ch13.coriolis_parameter is ROT.coriolis_parameter and ch13.OMEGA_EARTH == ROT.OMEGA_EARTH


SIGNATURES = {
    "thin_layer_terms": "U L H lat_rad drho_over_rho0 nu_H nu_v g Omega",
    "geostrophic_velocity": "dpdx dpdy f rho0",
    "thermal_wind_from_temperature": "dTdx dTdy f alpha g",
    "ekman_surface": "z tau_x tau_y rho nu_v f U_g V_g as_complex",
    "ekman_solve": "z K f tau rho V_g bottom top",
    "ekman_bottom": "z U_g V_g nu_v f as_complex",
    "shallow_water_omega": "k l c f0 beta",
    "kelvin_wave": "x y t eta0 k H f g direction",
    "geostrophic_adjustment_1d": "x eta0 H f g",
    "adjustment_energy": "eta0 H f g L rho",
    "inertia_gravity_omega": "k m N f l approx",
    "rossby_omega": "k l beta f0 c U",
    "eady_growth_rate": "k l N f H U0",
    "vertical_modes_shooting": "z N2_fn g n lid",
    "modes_uniform_N": "N H g n_modes lid nz",
    "linear_1d_run": "eta0 dx dt n_steps H f g save_every",
    "barotropic_run": "n L zeta0 t_end dt beta nu dealias save_every",
    "jet_section": "y z dT width lat_rad alpha T0 lapse u_surface",
    "flow_over_step": "x U beta f0 h0 h1",
    "eady_numeric_eigs": "k l N f H U0 n",
    "lapse_rate_table": "dT_dz Gamma_a",
}


def test_contract_V7_signatures_match_design_part_c():  # C.1–C.4
    for name, want in SIGNATURES.items():
        have = list(inspect.signature(getattr(ch13, name)).parameters)
        assert have[:len(want.split())] == want.split(), (name, have)
    p = inspect.signature(ch13.vertical_modes).parameters
    assert [p[k].default for k in ("g", "n_modes", "lid", "method")] == [G0, 4, "free", "fd"]
    for name, kwonly in (("thermal_wind_from_temperature", "alpha"), ("lapse_rate_table", "Gamma_a"),
                         ("barotropic_run", "t_end"), ("jet_section", "dT"), ("pressure_centre", "dp"),
                         ("eady_basic_state", "N")):
        assert inspect.signature(getattr(ch13, name)).parameters[kwonly].kind is inspect.Parameter.KEYWORD_ONLY


def test_contract_V7_scalars_in_floats_out_and_f_zero_raises():  # convention 2 of Part C
    scalar_calls = [
        (ch13.beta_parameter, (LAT,)), (ch13.beta_plane, (1e5, LAT)), (ch13.rossby_number, (14.0, F_N, 1.4e6)),
        (ch13.ekman_number, (7.0, F_S, 9e3)), (ch13.ekman_depth, (0.03, F_S)), (ch13.long_wave_speed, (4200.0,)),
        (ch13.poincare_omega, (1e-5, F_N, 200.0)), (ch13.rossby_radius, (3.1, F_S)), (ch13.rossby_omega, (-1e-5, 2e-6, 2e-11)),
        (ch13.potential_vorticity, (1e-5, F_N, 300.0)), (ch13.eady_alpha, (1e-6, 0.0, 1.1e-2, F_S)),
        (ch13.rhines_length, (13.0, 1.6e-11)), (ch13.inertial_period, (F_S,)), (ch13.kelvin_omega, (1e-5, 3.0)),
        (ch13.eady_growth_rate, (2e-6, 0.0, 1.1e-2, F_N, 9e3, 27.0)), (ch13.lee_wave_m, (10.0, 1e-2, 2e-4)),
    ]
    for fn, args in scalar_calls:
        out = fn(*args)
        assert isinstance(out, float) and math.isfinite(out), fn.__name__
    for out in ch13.geostrophic_velocity(1e-3, -2e-3, F_S, 1.2) + ch13.ekman_transport(0.07, 0.0, 1027.0, F_S):
        assert isinstance(out, float)
    assert isinstance(ch13.eady_phase_speed(1.3, 27.0), complex)
    for fn, args in ((ch13.geostrophic_velocity, (1.0, 1.0, 0.0, 1.2)), (ch13.rossby_number, (1.0, 0.0, 1.0)),
                     (ch13.ekman_depth, (1.0, 0.0)), (ch13.ekman_transport, (0.1, 0.0, 1000.0, 0.0)),
                     (ch13.rossby_radius, (3.0, 0.0)), (ch13.thermal_wind_shear, (1e-6, 0.0, 0.0, 1000.0)),
                     (ch13.geostrophic_adjustment_1d, (np.array([1.0]), 0.1, 10.0, 0.0)),
                     (ch13.kelvin_wave, (0.0, 0.0, 0.0, 0.1, 1e-5, 10.0, 0.0))):
        with pytest.raises(ValueError):
            fn(*args)
    assert ch13.inertial_period(0.0) == np.inf                       # the one documented inf


# =====================================================================================================================
# 1. §13.2–§13.4 — stratification, the equations of motion, the thin shell (C01)
# =====================================================================================================================
def test_planetary_V1_coriolis_parameter_beta_and_rotation_vector():  # V1, C01 (13.8), (13.10)
    lat = np.deg2rad(np.array([-73.0, -35.0, -4.0, 0.0, 11.0, 47.0, 88.0]))
    Om = ch13.OMEGA_EARTH
    assert close(ch13.coriolis_parameter(lat), 2 * Om * np.sin(lat), rtol=1e-15, atol=0)
    assert close(ch13.coriolis_parameter_deg(np.rad2deg(lat)), ch13.coriolis_parameter(lat), rtol=1e-14, atol=1e-20)
    assert close(ch13.beta_parameter(lat), 2 * Om * np.cos(lat) / ch13.EARTH_RADIUS_MEAN, rtol=1e-15)
    for la in lat:
        v = ch13.earth_rotation_local(la)
        assert abs(np.linalg.norm(v) - Om) < 1e-19 and v[0] == 0.0 and abs(2 * v[2] - ch13.coriolis_parameter(la)) < 1e-19
        assert ch13.f_plane(la) == ch13.coriolis_parameter(la)
        h = ch13.hemisphere(la)
        assert h["sign"] == np.sign(np.sin(la))
        assert (h["name"], h["turns"], h["cyclonic"]) == {1.0: ("northern", "right", "counter-clockwise"),
                                                          -1.0: ("southern", "left", "clockwise"),
                                                          0.0: ("equator", "none", "none")}[h["sign"]]
    assert abs(ch13.SIDEREAL_DAY - 86164.1) < 0.1                    # one turn against the stars
    assert abs(ch13.OMEGA_EARTH / ch13.OMEGA_SOLAR_DAY - (1 + 1 / 365.2422)) < 2e-6   # V1: one extra turn a year (trap T1)
    assert ch13.inertial_period(F_S) == ch13.inertial_period(F_N) == 2 * PI / F_N


def test_beta_plane_V1_error_is_second_order_in_y():  # V1 + V7, N17
    y = np.array([1e5, 2e5, 4e5, 8e5])
    err = np.abs(ch13.beta_plane_error(y, LAT))
    exact = (ch13.beta_plane(y, LAT) - 2 * ch13.OMEGA_EARTH * np.sin(LAT + y / ch13.EARTH_RADIUS_MEAN)) / (
        2 * ch13.OMEGA_EARTH * np.sin(LAT + y / ch13.EARTH_RADIUS_MEAN))
    assert close(ch13.beta_plane_error(y, LAT), exact, rtol=1e-12)
    assert abs(observed_order(y, err) - 2.0) < 0.1                   # Taylor remainder: (y/R)^2/2 to leading order
    assert close(err[0], 0.5 * (y[0] / ch13.EARTH_RADIUS_MEAN) ** 2, rtol=0.05)
    assert ch13.beta_plane(0.0, -LAT) == ch13.coriolis_parameter(-LAT) < 0


def test_coriolis_components_V2_derivation():  # V2, D02 (13.7)–(13.9)
    Om, th, u, v, w = sp.symbols("Omega theta u v w", real=True)
    cor = (2 * sp.Matrix([0, Om * sp.cos(th), Om * sp.sin(th)])).cross(sp.Matrix([u, v, w]))
    f = 2 * Om * sp.sin(th)
    assert sp.simplify(cor[0] - (2 * Om * (w * sp.cos(th) - v * sp.sin(th)))) == 0       # step: full x-component
    assert sp.simplify(cor[1] - f * u) == 0 and sp.simplify(cor[2] + 2 * Om * u * sp.cos(th)) == 0
    assert sp.simplify(cor[0].subs(w, 0) + f * v) == 0                                    # thin shell: drop w cos(theta)
    # the dropped term against the kept one is (W/U) cot(theta) = (H/L) cot(theta): small away from the equator
    ratio = sp.simplify((cor[0] - cor[0].subs(w, 0)) / (f * v))
    assert sp.simplify(ratio - w * sp.cos(th) / (v * sp.sin(th))) == 0
    rng = np.random.default_rng(13)
    for la in np.deg2rad([-61.0, -8.0, 23.0, 77.0]):
        uu = rng.standard_normal(3)
        full = np.cross(2 * ch13.earth_rotation_local(la), uu)
        assert close(ch13.coriolis_acceleration_local(*uu, la, thin=False), full, rtol=1e-13, atol=1e-18)
        thin = np.array(ch13.coriolis_acceleration_local(*uu, la))
        assert close(thin[1:], full[1:], rtol=1e-13) and abs(thin[0] - full[0] + 2 * ch13.OMEGA_EARTH * uu[2] * np.cos(la)) < 1e-18
    assert math.isclose(ch13.vertical_velocity_scale(14.0, 9e3, 1.4e6), 14.0 * 9e3 / 1.4e6)


def test_friction_force_V2_derivation():  # V2, D01 (13.5) -> (13.6); slip #2
    x, y, z = sp.symbols("x y z", real=True)
    nH, nv, rho = sp.symbols("nu_H nu_v rho", positive=True)
    # a solenoidal polynomial field (curl of a vector potential) and one that is not
    Apot = sp.Matrix([x * y ** 2 * z, x ** 2 * z ** 2 + y ** 3, x * y * z ** 2])
    sol = sp.Matrix([Apot[2].diff(y) - Apot[1].diff(z), Apot[0].diff(z) - Apot[2].diff(x), Apot[1].diff(x) - Apot[0].diff(y)])
    comp = sol + sp.Matrix([x ** 3, 0, 0])
    X = (x, y, z)

    def force(U):
        u, v, w = U
        t = sp.zeros(3, 3)
        t[0, 2] = t[2, 0] = rho * nv * u.diff(z) + rho * nH * w.diff(x)          # (13.5)
        t[1, 2] = t[2, 1] = rho * nv * v.diff(z) + rho * nH * w.diff(y)
        t[0, 1] = t[1, 0] = rho * nH * (u.diff(y) + v.diff(x))
        t[0, 0], t[1, 1], t[2, 2] = 2 * rho * nH * u.diff(x), 2 * rho * nH * v.diff(y), 2 * rho * nv * w.diff(z)
        F = [sum(t[i, j].diff(X[j]) for j in range(3)) / rho for i in range(3)]   # per unit mass: 1/rho (slip #2)
        eq6 = [nH * (q.diff(x, 2) + q.diff(y, 2)) + nv * q.diff(z, 2) for q in U]  # (13.6)
        return [sp.expand(a - b) for a, b in zip(F, eq6)], t

    assert sp.simplify(sum(sol[i].diff(X[i]) for i in range(3))) == 0
    d_sol, _ = force(sol)
    assert d_sol == [0, 0, 0]                                         # (13.6) holds once continuity is used
    d_comp, _ = force(comp)
    div = sum(comp[i].diff(X[i]) for i in range(3))
    assert sp.expand(d_comp[0] - nH * div.diff(x)) == 0 and d_comp[0] != 0        # the cross terms are nu grad(div u)
    assert sp.expand(d_comp[2] - nv * div.diff(z)) == 0
    eng = ch13.eddy_friction_force_sympy()
    assert eng["ok"] and eng["residual"] == 0
    # numeric tensor against (13.5) on a random velocity gradient; symmetric by construction
    G = np.random.default_rng(5).standard_normal((3, 3))
    tau = ch13.anisotropic_eddy_stress(G, 310.0, 0.02, 1027.0)
    assert close(tau, tau.T) and math.isclose(tau[0, 2], 1027.0 * (0.02 * G[0, 2] + 310.0 * G[2, 0]))
    assert math.isclose(tau[2, 2], 2 * 1027.0 * 0.02 * G[2, 2]) and math.isclose(tau[0, 1], 1027.0 * 310.0 * (G[0, 1] + G[1, 0]))
    assert math.isclose(ch13.eddy_friction(2e-9, -3e-4, 310.0, 0.02), 310.0 * 2e-9 + 0.02 * -3e-4)


def test_friction_force_V2_dimensions_printed_form_fails():  # V2 dimensions, slip #2 (printed form must fail)
    stress_gradient = Q_(1.0, "Pa") / Q_(1.0, "m")
    accel = Q_(310.0, "m**2/s") * Q_(1.0, "m/s") / Q_(1.0, "m") ** 2          # right-hand side of (13.6)
    assert not stress_gradient.check("[length]/[time]**2")                      # printed: d(tau)/dx is N/m^3
    assert (stress_gradient / Q_(1027.0, "kg/m**3")).check("[length]/[time]**2") and accel.check("[length]/[time]**2")
    assert ch13.friction_force_dimensions(printed=True)["is_acceleration"] is False
    good = ch13.friction_force_dimensions()
    assert good["is_acceleration"] is True and (good["M"], good["L"], good["T"]) == (0, 1, -2)


def test_boussinesq_V2_rest_state_printed_sign_fails():  # V2, slip #1 in (13.2)
    z, g, rho0 = sp.symbols("z g rho_0", positive=True)
    rb = sp.Function("rhobar")(z)
    pb = sp.Function("pbar")(z)
    hydro = {pb.diff(z): -rb * g}                                     # (13.3)
    corrected = (-pb.diff(z) / rho0 - g * rb / rho0).subs(hydro)      # vertical component of (13.2) at rest
    printed = (+pb.diff(z) / rho0 - g * rb / rho0).subs(hydro)
    assert sp.simplify(corrected) == 0
    assert sp.simplify(printed + 2 * g * rb / rho0) == 0 and sp.simplify(printed) != 0   # a resting fluid would accelerate at 2g rho/rho0
    assert ch13.boussinesq_rotating_sympy()["ok"] is True and ch13.boussinesq_rotating_sympy()["coriolis_work"] == 0
    bad = ch13.boussinesq_rotating_sympy(printed=True)
    assert bad["ok"] is False and bad["residual"] != 0
    assert ch13.perturbation_form_sympy()["ok"] and ch13.thin_layer_sympy()["ok"]


def test_thin_layer_V1_term_sizes_and_residuals():  # V1, C01 (N12, N14)
    U, L, H = I["U_syn"], I["L_syn"], I["atm_H"]
    for lat in (LAT, -LAT):
        t = ch13.thin_layer_terms(U, L, H, lat, nu_H=1e4, nu_v=5.0)
        f = float(ch13.coriolis_parameter(lat))
        assert t["f"] == f and math.isclose(t["W"], U * H / L) and math.isclose(t["Ro"], U / (abs(f) * L))
        assert math.isclose(t["Ro"], ch13.rossby_number(U, f, L)) and math.isclose(t["aspect"], H / L)
        assert math.isclose(t["x"]["acceleration"] / t["x"]["coriolis"], t["Ro"])
        assert math.isclose(t["x"]["coriolis_w"] / t["x"]["coriolis"], (H / L) / abs(np.tan(lat)))   # the dropped term (D02)
        assert math.isclose(t["x"]["friction_v"], 5.0 * U / H ** 2) and math.isclose(t["x"]["friction_H"], 1e4 * U / L ** 2)
        assert t["x"]["pressure"] == max(v for k, v in t["x"].items() if k != "pressure")
        assert t["z"]["buoyancy"] > 1e3 * t["z"]["acceleration"]                 # hydrostatic to a part in a thousand
        tab = ch13.term_table_thin_layer(U, L, H, lat, nu_H=1e4, nu_v=5.0)
        assert math.isclose(tab.loc["acceleration", "x_over_coriolis"], t["Ro"]) and tab.loc["coriolis", "x_over_coriolis"] == 1.0
        assert math.isclose(tab.loc["buoyancy", "z"], t["z"]["buoyancy"]) and np.isnan(tab.loc["buoyancy", "x"])
    with pytest.raises(ValueError):
        ch13.thin_layer_terms(U, L, H, 0.0)
    # a geostrophic, hydrostatic state leaves no residual in the thin-shell equations (13.9), either hemisphere
    for f in HEMI:
        u, v = ch13.geostrophic_velocity(2.1e-3, -0.7e-3, f, 1.2)
        r = ch13.thin_layer_residuals(u=u, v=v, f=f, dpdx=2.1e-3, dpdy=-0.7e-3, dpdz=-0.05 * G0, rho=0.05, rho0=1.2)
        assert all(abs(r[c]["residual"]) < 1e-15 for c in "xyz")
        r2 = ch13.thin_layer_residuals(u=u, v=v, f=-f, dpdx=2.1e-3, dpdy=-0.7e-3, rho0=1.2)   # wrong hemisphere: not balanced
        assert abs(r2["x"]["residual"]) > 1e-3
    assert math.isclose(ch13.scale_height(331.0), 331.0 ** 2 / G0)


def test_stratification_V1_density_budget_profiles_and_atmosphere():  # V1, N03, N04, N06 (13.1)
    b = ch13.ocean_density_gradient_budget(-4.9e-3, 1031.0, 1507.0)
    assert math.isclose(b["adiabatic"], -G0 * 1031.0 / 1507.0 ** 2) and math.isclose(b["potential"], -4.9e-3 - b["adiabatic"])
    assert math.isclose(b["compression_share"], b["adiabatic"] / -4.9e-3) and b["potential"] < 0 < b["compression_share"] < 1
    assert math.isclose(b["potential"], float(STRAT.ocean_potential_density_gradient(-4.9e-3, 1031.0, 1507.0)), rel_tol=1e-12)
    z = np.linspace(-3000.0, 0.0, 1201)
    o = ch13.ocean_profile_idealized(z)
    N2_num = ch13.buoyancy_frequency_sq(z, o["rho_theta"], 1027.0)
    below = z < -60.0
    assert np.max(np.abs(N2_num[below] - o["N2"][below])) < 2e-4 * np.max(o["N2"])      # N^2 = -(g/rho0) d(rho)/dz
    assert np.all(o["N2"][z > -50.0] == 0.0) and np.all(np.diff(o["T"]) >= 0) and close(o["N"] ** 2, o["N2"])
    hs, es = [], []                                                    # V3: the stencil is second order up to the ends
    for n in (41, 81, 161, 321):
        zz = np.linspace(-1000.0, 0.0, n)
        rho = 1027.0 * (1 + 2e-3 * np.exp(zz / 300.0))
        es.append(np.max(np.abs(ch13.buoyancy_frequency_sq(zz, rho, 1027.0) + G0 * 2e-3 / 300.0 * np.exp(zz / 300.0))))
        hs.append(zz[1] - zz[0])
    assert abs(observed_order(hs, es) - 2.0) < 0.15
    th = ch13.thermocline_N2(z)
    assert math.isclose(th.max(), 8.0e-3 ** 2, rel_tol=1e-3) and abs(z[np.argmax(th)] + 300.0) < 3.0 and np.all(th > 0)
    a = ch13.atmosphere_layers(np.array([3000.0, 16000.0, 28000.0]))
    assert list(a["layer"]) == ["troposphere", "stratosphere", "stratosphere"] and a["dT_dz"][0] == -6.5e-3
    assert np.all(a["N2"] > 0) and a["N2"][1] > 2 * a["N2"][0]            # the stratosphere is far more stable
    assert close(a["N2"], STRAT.brunt_vaisala_sq_from_lapse(a["T"], a["dT_dz"]))


def test_lapse_rate_V7_both_conventions_same_verdict():  # V7 invariance under the sign convention (project rule)
    Ga = float(STRAT.adiabatic_lapse_rate())                              # Kundu sign, about -9.8 K/km
    assert Ga < 0 and math.isclose(Ga, -G0 / STRAT.CP_AIR)
    with pytest.raises(TypeError):
        ch13.lapse_rate_table(-6.5e-3)                                    # Gamma_a is a required keyword
    for dT, verdict in ((-6.5e-3, "stable"), (-11.7e-3, "unstable"), (Ga, "neutral"), (2.3e-3, "stable")):
        t = ch13.lapse_rate_table(dT, Gamma_a=Ga)
        assert list(t.index) == ["Kundu", "meteorology"] and list(t["verdict"]) == [verdict, verdict]
        assert t.loc["meteorology", "value_K_per_km"] == -t.loc["Kundu", "value_K_per_km"]
        assert t.loc["meteorology", "adiabatic_K_per_km"] == -t.loc["Kundu", "adiabatic_K_per_km"] > 0
        k = STRAT.lapse_rate_stability(dT, Ga, convention="kundu")
        m = STRAT.lapse_rate_stability(dT, Ga, convention="meteorology")
        assert k.verdict == m.verdict == verdict and k.code == m.code and m.Gamma == -k.Gamma
        sign = np.sign(float(STRAT.brunt_vaisala_sq_from_lapse(260.0, dT)))
        assert verdict == "neutral" or sign == k.code                     # the verdict is the sign of N^2
    t = ch13.lapse_rate_table(-6.5e-3, Gamma_a=Ga)
    assert ">" in t.loc["Kundu", "criterion"] and "<" in t.loc["meteorology", "criterion"]


def test_tables_V7_inputs_conventions_slips_traps_and_names():  # smoke, C.4 tables
    assert ch13.illustrative_inputs() == I and ch13.illustrative_inputs() is not I             # a fresh dict each call
    assert math.isclose(I["lat"], np.deg2rad(35.0)) and I["wavelengths"] == (7.3e5, 3.1e6) and len(I) == 21
    ct = ch13.conventions_table()
    assert list(ct.columns) == ["symbol", "meanings", "ours", "code"] and len(ct) >= 15 and "f" in set(ct["symbol"])
    slips = ch13.book_slips()
    assert list(slips) == [str(i) for i in range(1, 14)]
    for row in slips.values():
        assert {"where", "printed", "correct", "taught_in"} <= set(row) and not re.search(r"\d\.\d{2,}", row["printed"])
    tr = ch13.traps()
    assert len(tr) == 17 and [t["id"] for t in tr] == [f"T{i}" for i in range(1, 18)]
    w = ch13.wind_from_to(9.0, 0.0)
    assert (w["wind_name"], w["current_name"], w["from_deg"], w["to_deg"]) == ("westerly", "eastward", 270.0, 90.0)
    assert ch13.wind_from_to(0.0, -4.0)["wind_name"] == "northerly" and ch13.wind_from_to(0.0, 0.0)["wind_name"] == "calm"
    v = ch13.viscous_layer_thicknesses(1.3e-6, t=77.0, x=0.9, U=0.4, f=F_S)
    assert math.isclose(v["diffusive"], math.sqrt(1.3e-6 * 77.0)) and math.isclose(v["boundary_layer"], math.sqrt(1.3e-6 * 0.9 / 0.4))
    assert math.isclose(v["ekman"], math.sqrt(2 * 1.3e-6 / F_N)) and ch13.viscous_layer_thicknesses(1e-6)["ekman"] is None
    assert set(ch13.sympy_engines()) == {"boussinesq_rotating", "perturbation_form", "eddy_friction_force", "thin_layer",
                                         "taylor_proudman", "hydrostatic_linear_set", "v_equation",
                                         "rotating_internal_wave_set", "w_equation_rotating", "pv_conservation",
                                         "qg_vorticity", "eady_qg"}


# =====================================================================================================================
# 2. §13.5 — geostrophic balance, thermal wind, Taylor–Proudman (C02, C03)
# =====================================================================================================================
def test_geostrophic_V1_flow_along_isobars_both_hemispheres():  # V1 + V7, C02 (13.11)–(13.13)
    x = np.linspace(-1.5e6, 1.5e6, 61)
    for lat, kind in ((LAT, "low"), (-LAT, "low"), (LAT, "high")):
        pc = ch13.pressure_centre(x, x, dp=1900.0, R_c=4.1e5, lat_rad=lat, kind=kind)
        X, Y = np.meshgrid(x, x)
        f = pc["f"]
        assert f == float(ch13.coriolis_parameter(lat))
        px, py = -pc["p"] * X / 4.1e5 ** 2, -pc["p"] * Y / 4.1e5 ** 2
        assert np.max(np.abs(pc["u"] * px + pc["v"] * py)) < 1e-12 * np.max(np.abs(px)) * np.max(np.abs(pc["u"]))   # u . grad p = 0
        assert close(-f * pc["v"], -px / 1.2, rtol=1e-12, atol=1e-18) and close(f * pc["u"], -py / 1.2, rtol=1e-12, atol=1e-18)
        circ = np.sign((X * pc["v"] - Y * pc["u"])[np.hypot(X, Y) > 1e5])  # sense of rotation about the centre
        want = (1.0 if kind == "low" else -1.0) * np.sign(f)               # a low is cyclonic: counter-clockwise for f > 0
        assert np.all(circ == want)
        psi = ch13.geostrophic_streamfunction(pc["p"], f, 1.2)             # u = -dpsi/dy, v = dpsi/dx (trap T14)
        dx = x[1] - x[0]
        assert close(pc["u"][1:-1, :], -(psi[2:, :] - psi[:-2, :]) / (2 * dx), atol=2e-2 * np.max(np.abs(pc["u"])))
    u, v = ch13.geostrophic_velocity(1.7e-3, 0.0, F_N, 1.2)
    us, vs = ch13.geostrophic_velocity(1.7e-3, 0.0, F_S, 1.2)
    assert u == 0.0 and v > 0 and vs == -v                                 # low pressure on the left (north), on the right (south)
    assert math.isclose(ch13.rossby_number(I["U_syn"], F_S, I["L_syn"]), I["U_syn"] / (F_N * I["L_syn"]))
    assert math.isclose(ch13.ekman_number(7.0, F_S, 9e3), 7.0 / (F_N * 9e3 ** 2))


def test_geostrophic_V3_gridded_gradient_is_second_order():  # V3, C02 (13.11)–(13.12), (13.116)
    hs, eu, eh = [], [], []
    for n in (24, 48, 96, 192):
        L = 2.0e6
        x = np.linspace(0.0, L, n + 1)
        X, Y = np.meshgrid(x, x)
        p = 800.0 * np.sin(2 * PI * X / L + 0.4) * np.cos(1.5 * PI * Y / L)
        u, v = ch13.geostrophic_from_field(p, x[1], x[1], F_S, 1.2)
        ue = 800.0 * 1.5 * PI / L * np.sin(2 * PI * X / L + 0.4) * np.sin(1.5 * PI * Y / L) / (1.2 * F_S)
        ve = 800.0 * 2 * PI / L * np.cos(2 * PI * X / L + 0.4) * np.cos(1.5 * PI * Y / L) / (1.2 * F_S)
        eu.append(max(np.max(np.abs(u - ue)), np.max(np.abs(v - ve))))
        uh, vh = ch13.geostrophic_from_height(p / (1.2 * G0), x[1], x[1], F_S)
        eh.append(max(np.max(np.abs(uh - ue * 1.2 / 1.2)), np.max(np.abs(vh - ve))))
        hs.append(x[1])
    assert abs(observed_order(hs, eu) - 2.0) < 0.15 and abs(observed_order(hs, eh) - 2.0) < 0.15   # edges included
    assert eu[-1] < 2e-3 * 800.0 * 2 * PI / 2.0e6 / (1.2 * F_N)


def test_parcel_adjust_V2_ode_and_geostrophic_limit():  # V2 + V7, ours (curation note 1)
    t, r, f = sp.symbols("t r f", real=True)
    Gs, V0 = sp.symbols("G V0")
    lam = r + sp.I * f
    V = -Gs / lam + (V0 + Gs / lam) * sp.exp(-lam * t)
    assert sp.simplify(V.diff(t) + lam * V + Gs) == 0 and sp.simplify(V.subs(t, 0) - V0) == 0   # dV/dt + (r + if)V = -G
    tt = np.linspace(0.0, 6e5, 2001)
    for fs in HEMI:
        Gc = (1.3e-3 + 0.4e-3j) / 1.2
        V = ch13.parcel_adjust(tt, Gc, fs, r=2e-5)
        num = solve_ivp(lambda s, y: -(2e-5 + 1j * fs) * y - Gc, (0, tt[-1]), [0j], t_eval=tt, rtol=1e-10, atol=1e-13).y[0]
        assert np.max(np.abs(V - num)) < 1e-6 * np.max(np.abs(V))
        ug, vg = ch13.geostrophic_velocity(1.3e-3, 0.4e-3, fs, 1.2)
        Vinf = ch13.parcel_adjust(2e9, Gc, fs, r=2e-8)                      # vanishing drag: settles on the geostrophic wind
        assert abs(Vinf - (ug + 1j * vg)) < 1.5 * (2e-8 / F_N) * abs(ug + 1j * vg)   # the offset is r/|f| exactly
        osc = ch13.parcel_adjust(tt, Gc, fs)                                # no drag at all: an inertial circle about it
        assert close(np.abs(osc - (ug + 1j * vg)), abs(ug + 1j * vg), rtol=1e-12)
        Vd = ch13.parcel_adjust(1e9, Gc, fs, r=3e-5)                        # with drag: a component down the pressure gradient
        assert np.real(Vd * np.conj(Gc)) < 0
    assert ch13.parcel_adjust(10.0, 2.0 + 0j, 0.0) == -20.0                 # f = r = 0: free fall down the gradient


def test_thermal_wind_V2_derivation():  # V2, D04 (13.15) and its temperature form
    x, y, z = sp.symbols("x y z", real=True)
    f, g, rho0, al = sp.symbols("f g rho_0 alpha", nonzero=True)
    P = sp.Function("P")(x, y, z)
    v, u = P.diff(x) / (rho0 * f), -P.diff(y) / (rho0 * f)                  # (13.11), (13.12)
    rho = -P.diff(z) / g                                                    # (13.14)
    assert sp.simplify(v.diff(z) + g / (rho0 * f) * rho.diff(x)) == 0       # (13.15a)
    assert sp.simplify(u.diff(z) - g / (rho0 * f) * rho.diff(y)) == 0       # (13.15b)
    T = sp.Function("T")(x, y, z)
    rho_T = -rho0 * al * T                                                  # linear equation of state, perturbation form
    assert sp.simplify(g / (rho0 * f) * rho_T.diff(y) + g * al / f * T.diff(y)) == 0    # ours: du/dz = -(g alpha/f) dT/dy
    assert sp.simplify(-g / (rho0 * f) * rho_T.diff(x) - g * al / f * T.diff(x)) == 0   # ours: dv/dz = +(g alpha/f) dT/dx
    for fs in HEMI:
        du, dv = ch13.thermal_wind_shear(3e-7, -5e-7, fs, 1.2)
        assert math.isclose(du, G0 / (1.2 * fs) * -5e-7) and math.isclose(dv, -G0 / (1.2 * fs) * 3e-7)
        a = 1 / 280.0
        duT, dvT = ch13.thermal_wind_from_temperature(3e-7 / (-1.2 * a), -5e-7 / (-1.2 * a), fs, alpha=a)
        assert math.isclose(duT, du, rel_tol=1e-12) and math.isclose(dvT, dv, rel_tol=1e-12)   # same physics, rho' = -rho0 alpha T'
    du_n, _ = ch13.thermal_wind_shear(0.0, 4e-7, F_N, 1.2)                  # denser (colder) air toward the pole, north
    du_s, _ = ch13.thermal_wind_shear(0.0, -4e-7, F_S, 1.2)                 # denser air toward the pole, south
    assert du_n > 0 and du_s > 0                                            # westerlies increase with height in BOTH hemispheres
    with pytest.raises(TypeError):
        ch13.thermal_wind_from_temperature(1e-6, 0.0, F_N, 1 / 280.0)       # alpha is keyword-only


def test_thermal_wind_V1_jet_section_in_exact_balance():  # V1 + V4, C03 (13.15b), jet_section, thermal_wind_integrate
    y = np.linspace(-1.5e6, 1.5e6, 121)
    z = np.linspace(0.0, 1.1e4, 45)
    a = 1 / 270.0
    for lat in (LAT, -LAT):
        f = float(ch13.coriolis_parameter(lat))
        js = ch13.jet_section(y, z, dT=17.0, width=3.3e5, lat_rad=lat, alpha=a, u_surface=2.0)
        assert js["T"].shape == js["U"].shape == (45, 121)
        dTdy_num = np.gradient(js["T"], y, axis=1)[:, 2:-2]
        assert np.max(np.abs(dTdy_num - js["dTdy"][None, 2:-2])) < 4e-3 * np.max(np.abs(js["dTdy"]))
        dUdz = np.gradient(js["U"], z, axis=0)
        assert close(dUdz, (-G0 * a / f) * js["dTdy"][None, :] * np.ones_like(dUdz), rtol=1e-10, atol=1e-16)
        assert np.all(js["U"][-1] - 2.0 >= 0) and js["U"][-1, 60] == js["U"].max()      # a westerly jet in both hemispheres
        assert np.all(np.diff(js["T"][0] * np.sign(lat)) <= 0)                           # warm toward the equator
        drho_dy = -1.2 * a * js["dTdy"][60]
        ui = ch13.thermal_wind_integrate(z, drho_dy, f, 1.2, u_ref=2.0)
        assert close(ui, js["U"][:, 60], rtol=1e-12)
    zz = np.linspace(0.0, 9e3, 91)                                         # depth-dependent gradient: trapezoid, order 2
    ex = G0 / (1.2 * F_N) * 4e-7 * 3e3 * (1 - np.exp(-zz / 3e3))
    assert np.max(np.abs(ch13.thermal_wind_integrate(zz, 4e-7 * np.exp(-zz / 3e3), F_N, 1.2) - ex)) < 1e-4 * ex[-1]


def test_taylor_proudman_V2_derivation():  # V2, D05 (13.16)–(13.21); slip #3 (printed sign fails)
    x, y, z = sp.symbols("x y z", real=True)
    Om, rho, g = sp.symbols("Omega rho g", positive=True)
    q = sp.Function("q")(x, y)
    p = q - rho * g * z                                                    # (13.14) with uniform density
    v = p.diff(x) / (2 * Om * rho)                                         # (13.16): -2 Omega v = -(1/rho) dp/dx
    u_ok = -p.diff(y) / (2 * Om * rho)                                     # corrected (13.17): +2 Omega u = -(1/rho) dp/dy
    u_pr = +p.diff(y) / (2 * Om * rho)                                     # as printed:        -2 Omega u = -(1/rho) dp/dy
    assert sp.simplify(u_ok.diff(x) + v.diff(y)) == 0                      # cross-differentiation -> zero divergence -> dw/dz = 0
    assert u_ok.diff(z) == 0 and v.diff(z) == 0                            # (13.20)
    div_printed = sp.simplify(u_pr.diff(x) + v.diff(y))
    assert sp.simplify(div_printed - q.diff(x, y) / (Om * rho)) == 0 and div_printed != 0   # printed sign: divergence survives
    assert sp.simplify(u_pr.diff(x) - v.diff(y)) == 0                      # what the printed pair gives instead
    # consistency with the chapter's own (13.12) at f = 2 Omega decides which sign is the slip
    f = sp.Symbol("f")
    assert sp.simplify((-p.diff(y) / (rho * f)).subs(f, 2 * Om) - u_ok) == 0
    good, bad = ch13.taylor_proudman_sympy(), ch13.taylor_proudman_sympy(printed=True)
    assert good["ok"] and good["du_dz"] == 0 and good["dv_dz"] == 0 and bad["ok"] is False and bad["printed_identity"] == 0
    col = ch13.taylor_proudman_residual(lambda a, b, c: (-0.3 * b, 0.3 * a, 0.05), (1.0, 2.0, 3.0))
    assert all(abs(col[k]) < 1e-12 for k in ("du_dz", "dv_dz", "dw_dz"))
    sheared = ch13.taylor_proudman_residual(lambda a, b, c: (0.7 * c, 0.0, 0.0), (1.0, 2.0, 3.0))
    assert math.isclose(sheared["du_dz"], 0.7, rel_tol=1e-9)


# =====================================================================================================================
# 3. §13.6–§13.7 — Ekman layers (C04, C05, C06)
# =====================================================================================================================
EK = dict(tau=I["tau"], rho=I["rho_ocean"], nu=I["nu_v_ocean"])


def test_ekman_surface_V2_derivation():  # V2, D06 (13.22)–(13.29); slip #4 (the root bounded upward fails)
    z = sp.Symbol("z", real=True)
    f, nu, rho, tau = sp.symbols("f nu_v rho tau", positive=True)
    u, v = sp.Function("u")(z), sp.Function("v")(z)
    V = u + sp.I * v
    e22, e23 = -f * v - nu * u.diff(z, 2), f * u - nu * v.diff(z, 2)          # (13.22), (13.23)
    assert sp.expand(e22 + sp.I * e23 - (sp.I * f * V - nu * V.diff(z, 2))) == 0   # i x (13.23) + (13.22) = (13.27)
    d = sp.sqrt(2 * nu / f)                                                    # (13.29)
    lam = (1 + sp.I) / d
    assert sp.simplify(lam ** 2 - sp.I * f / nu) == 0                          # both exponents of (13.28) solve (13.27)
    A = tau * d * (1 - sp.I) / (2 * rho * nu)
    Vs = A * sp.exp(lam * z)
    assert sp.simplify(rho * nu * Vs.diff(z).subs(z, 0) - tau) == 0            # (13.24)–(13.25): rho nu dV/dz = tau, real
    amp = tau / (rho * sp.sqrt(f * nu)) * sp.exp(z / d)
    u_pr = amp * sp.cos(-z / d + sp.pi / 4)                                    # the components as printed (f > 0)
    v_pr = -amp * sp.sin(-z / d + sp.pi / 4)
    assert sp.simplify(sp.expand_complex(Vs) - (u_pr + sp.I * v_pr)).rewrite(sp.exp).simplify() == 0
    assert sp.limit(sp.Abs(sp.exp(z / d)), z, -sp.oo) == 0                     # kept root decays into the ocean (z -> -oo)
    assert sp.limit(sp.exp(-z / d), z, -sp.oo) == sp.oo                        # the other root: bounded only as z -> +oo
    # numeric: the coded spiral against the printed components (f > 0), the governing equations and both conditions
    zz = np.linspace(-9.0, 0.0, 721) * float(ch13.ekman_depth(EK["nu"], F_N))
    uc, vc = ch13.ekman_surface(zz, EK["tau"], 0.0, EK["rho"], EK["nu"], F_N)
    up, vp = ch13.ekman_surface_printed(zz, EK["tau"], EK["rho"], EK["nu"], F_N)
    assert close(uc, up, rtol=0, atol=1e-15) and close(vc, vp, rtol=0, atol=1e-15)
    with pytest.raises(ValueError):
        ch13.ekman_surface_printed(zz, EK["tau"], EK["rho"], EK["nu"], F_S)   # true as printed for f > 0 only: a trap, not a slip
    assert np.all(np.isnan(ch13.ekman_surface(np.array([1.0, 30.0]), EK["tau"], 0.0, EK["rho"], EK["nu"], F_N)[0]))


def test_ekman_surface_V1_spiral_residual_stress_and_hemispheres():  # V1 + V7, C04
    for f in HEMI:
        s = np.sign(f)
        d = float(ch13.ekman_depth(EK["nu"], f))
        assert math.isclose(d, math.sqrt(2 * EK["nu"] / abs(f))) and math.isclose(ch13.ekman_depth(EK["nu"], f, "pi"), PI * d)
        zz = np.linspace(-10 * d, 0.0, 4001)
        for tx, ty in ((EK["tau"], 0.0), (-0.03, 0.05)):
            V = ch13.ekman_surface(zz, tx, ty, EK["rho"], EK["nu"], f, as_complex=True)
            rx, ry = ch13.ekman_residual(zz, V.real, V.imag, EK["nu"], f)
            assert max(np.max(np.abs(rx)), np.max(np.abs(ry))) < 2e-5 * abs(f) * np.max(np.abs(V))    # (13.22)–(13.23), O(dz^2)
            dV0 = (3 * V[-1] - 4 * V[-2] + V[-3]) / (2 * (zz[1] - zz[0]))
            assert abs(EK["rho"] * EK["nu"] * dV0 - (tx + 1j * ty)) < 2e-5 * abs(tx + 1j * ty)         # surface stress
            ang = np.angle(V[-1] / (tx + 1j * ty))
            assert math.isclose(ang, -s * PI / 4, abs_tol=1e-12)                # 45 degrees to the right (f > 0), left (f < 0)
            assert math.isclose(abs(V[-1]), abs(tx + 1j * ty) / (EK["rho"] * math.sqrt(abs(f) * EK["nu"])), rel_tol=1e-12)
            assert close(np.abs(V), abs(V[-1]) * np.exp(zz / d), rtol=1e-12)                    # e-folding scale delta
            turn = np.unwrap(np.angle(V))                                        # rotates clockwise with depth for f > 0
            assert close(np.gradient(turn, zz)[5:-5], s / d, rtol=1e-6)
            Vg = ch13.ekman_surface(zz, tx, ty, EK["rho"], EK["nu"], f, U_g=0.2, V_g=-0.1, as_complex=True)
            assert close(Vg - V, 0.2 - 0.1j)
    zb = np.linspace(-200.0, 0.0, 4001)                                          # (13.31): tilting = diffusion of vorticity
    for f in HEMI:
        b = ch13.ekman_vorticity_balance(zb, EK["tau"], EK["rho"], EK["nu"], f)
        ub, vb = ch13.ekman_surface(zb, EK["tau"], 0.0, EK["rho"], EK["nu"], f)
        dz = zb[1] - zb[0]
        d1 = lambda a: (a[2:] - a[:-2]) / (2 * dz)  # noqa: E731
        d2 = lambda a: (a[2:] - 2 * a[1:-1] + a[:-2]) / dz ** 2  # noqa: E731
        assert close(b["omega_y"][1:-1], d1(ub), atol=1e-5 * np.max(np.abs(b["omega_y"])))      # omega_y = du/dz
        assert close(b["omega_x"][1:-1], -d1(vb), atol=1e-5 * np.max(np.abs(b["omega_x"])))     # omega_x = -dv/dz
        assert close(b["tilt_x"], b["diff_x"], rtol=1e-12, atol=1e-22) and close(b["tilt_y"], b["diff_y"], rtol=1e-12, atol=1e-22)
        assert close(b["tilt_x"][1:-1], -f * d1(vb), atol=1e-5 * np.max(np.abs(b["tilt_x"])))   # -f dv/dz
        assert close(b["diff_x"][1:-1], EK["nu"] * d2(b["omega_y"]), atol=1e-4 * np.max(np.abs(b["diff_x"])))


def test_ekman_transport_V2_derivation():  # V2, D07 (13.30) two ways
    z = sp.Symbol("z", real=True)
    f, nu, rho, tau = sp.symbols("f nu_v rho tau", positive=True)
    d = sp.sqrt(2 * nu / f)
    lam = (1 + sp.I) / d
    Vs = tau * d * (1 - sp.I) / (2 * rho * nu) * sp.exp(lam * z)
    M = sp.simplify(sp.integrate(Vs, (z, -sp.oo, 0), conds="none"))            # route 1: integrate the spiral
    assert sp.simplify(sp.re(M)) == 0 and sp.simplify(sp.im(M) + tau / (rho * f)) == 0
    # route 2: integrate -rho f v = d(tau_x)/dz and rho f u = d(tau_y)/dz once; only "stress -> 0 at depth" is needed
    tx0, ty0 = sp.symbols("tau_x0 tau_y0", real=True)
    My, Mx = -tx0 / (rho * f), ty0 / (rho * f)
    assert sp.simplify(My.subs(tx0, tau) - sp.im(M)) == 0 and Mx.subs(ty0, 0) == sp.re(M)
    assert (Q_(0.07, "Pa") / (Q_(1027.0, "kg/m**3") * Q_(8e-5, "1/s"))).check("[length]**2/[time]")   # m^2/s per unit width


def test_ekman_transport_V1_quadrature_partial_and_hemispheres():  # V1 + V7, C05 (13.30), D07
    for f in HEMI:
        d = float(ch13.ekman_depth(EK["nu"], f))
        for tx, ty in ((EK["tau"], 0.0), (0.02, -0.06)):
            Mx, My = ch13.ekman_transport(tx, ty, EK["rho"], f)
            qx = quad(lambda s: ch13.ekman_surface(s, tx, ty, EK["rho"], EK["nu"], f)[0], -40 * d, 0, epsabs=1e-13, limit=400)[0]
            qy = quad(lambda s: ch13.ekman_surface(s, tx, ty, EK["rho"], EK["nu"], f)[1], -40 * d, 0, epsabs=1e-13, limit=400)[0]
            assert abs(qx - Mx) < 1e-9 * math.hypot(Mx, My) and abs(qy - My) < 1e-9 * math.hypot(Mx, My)
            assert abs(Mx * tx + My * ty) < 1e-15                                # at right angles to the stress
            assert np.sign(tx * My - ty * Mx) == -np.sign(f)                     # to the right for f > 0, to the left for f < 0
            zz = np.linspace(-6 * d, 0.0, 13)
            px, py = ch13.ekman_transport_partial(zz, tx, ty, EK["rho"], EK["nu"], f)
            for j in (2, 7):
                q = quad(lambda s: ch13.ekman_surface(s, tx, ty, EK["rho"], EK["nu"], f)[1], zz[j], 0, epsabs=1e-14)[0]
                assert abs(q - py[j]) < 1e-10 * math.hypot(Mx, My)
            assert px[-1] == 0.0 and py[-1] == 0.0
            fx, fy = ch13.ekman_transport_partial(-60 * d, tx, ty, EK["rho"], EK["nu"], f)
            assert math.isclose(fx, Mx, abs_tol=1e-12) and math.isclose(fy, My, abs_tol=1e-12)
    Mx, My = ch13.ekman_transport(I["tau"], 0.0, I["rho_ocean"], F_N)
    assert Mx == 0.0 and math.isclose(My, -I["tau"] / (I["rho_ocean"] * F_N))    # independent of nu_v
    up = ch13.coastal_upwelling(0.09, "east", LAT)                               # northward wind, land to the east, north
    dn = ch13.coastal_upwelling(-0.09, "east", LAT)
    so = ch13.coastal_upwelling(0.09, "east", -LAT)
    assert up["upwelling"] is False and dn["upwelling"] is True and so["upwelling"] is True
    assert math.isclose(dn["transport_offshore"], 0.09 / (1027.0 * F_N)) and up["kelvin_direction"] == "poleward"
    assert ch13.coastal_upwelling(0.09, "west", LAT)["upwelling"] is True


def test_ekman_pumping_V2_derivation():  # V2, D09 (ours): w_E = (1/rho) curl_z(tau/f); w = sign(f) (delta/2) zeta_g
    x, y = sp.symbols("x y", real=True)
    rho, nu = sp.symbols("rho nu_v", positive=True)
    f = sp.Function("f")(y)
    tx, ty = sp.Function("tau_x")(x, y), sp.Function("tau_y")(x, y)
    Mx, My = ty / (rho * f), -tx / (rho * f)                                    # (13.30) at each point
    wE = Mx.diff(x) + My.diff(y)                                               # continuity integrated over the layer
    assert sp.simplify(wE - ((ty / f).diff(x) - (tx / f).diff(y)) / rho) == 0
    curl = ty.diff(x) - tx.diff(y)
    assert sp.simplify(wE - (curl / (rho * f) + f.diff(y) * tx / (rho * f ** 2))) == 0   # the beta term of the general form
    psi = sp.Function("psi", real=True)(x, y)                                  # non-divergent geostrophic flow above the bottom
    Ug, Vg = -psi.diff(y), psi.diff(x)
    zeta = Vg.diff(x) - Ug.diff(y)
    d = sp.Symbol("delta", positive=True)
    for s in (1, -1):
        assert sp.simplify(1 / (1 + sp.I * s) - (1 - sp.I * s) / 2) == 0
        Mx, My = -d / 2 * (Ug + s * Vg), -d / 2 * (Vg - s * Ug)                # real and imaginary parts of -(U + iV) delta/(1 + is)
        assert sp.expand(Mx + sp.I * My + (Ug + sp.I * Vg) * d * (1 - sp.I * s) / 2) == 0
        w_top = -(Mx.diff(x) + My.diff(y))                                     # continuity integrated over the layer
        assert sp.simplify(w_top - s * d / 2 * zeta) == 0
    # numeric, gridded: a cyclonic (f > 0) stress curl gives upwelling; mirrored in the south
    hs, es = [], []
    for n in (33, 65, 129):
        xs = np.linspace(-1.0e6, 1.0e6, n)
        X, Y = np.meshgrid(xs, xs)
        r2 = (X ** 2 + Y ** 2) / 3.0e5 ** 2
        txn, tyn = -0.1 * Y / 3.0e5 * np.exp(-r2), 0.1 * X / 3.0e5 * np.exp(-r2)            # counter-clockwise stress
        w = ch13.ekman_pumping(txn, tyn, xs[1] - xs[0], xs[1] - xs[0], 1027.0, F_N)
        exact = 0.1 / 3.0e5 * (2 - 2 * r2) * np.exp(-r2) / (1027.0 * F_N)
        es.append(np.max(np.abs(w - exact)))
        hs.append(xs[1] - xs[0])
        ws = ch13.ekman_pumping(txn, tyn, xs[1] - xs[0], xs[1] - xs[0], 1027.0, F_S)
        assert close(ws, -w)
    assert abs(observed_order(hs, es) - 2.0) < 0.15 and w[n // 2, n // 2] > 0
    assert math.isclose(ch13.ekman_pumping_from_curl(2e-7, 1027.0, F_N), 2e-7 / (1027.0 * F_N))
    assert math.isclose(ch13.ekman_pumping_from_curl(0.0, 1027.0, F_N, beta=2e-11, tau_x=0.1), 2e-11 * 0.1 / (1027.0 * F_N ** 2))
    assert math.isclose(ch13.sverdrup_transport(-1e-7, 1027.0, 2e-11), -1e-7 / (1027.0 * 2e-11))
    dlt = float(ch13.ekman_depth(I["nu_v_atm"], F_N))
    assert math.isclose(ch13.ekman_pumping_bottom(1.5e-5, I["nu_v_atm"], F_N), 0.5 * dlt * 1.5e-5)
    assert math.isclose(ch13.ekman_pumping_bottom(-1.5e-5, I["nu_v_atm"], F_S), 0.5 * dlt * 1.5e-5)   # cyclonic in the south is zeta < 0


def test_ekman_bottom_V2_derivation():  # V2, D08 (13.33)–(13.41), transport, 45 degrees, the overshoot
    z = sp.Symbol("z", positive=True)
    f, nu, U = sp.symbols("f nu_v U", positive=True)
    d = sp.sqrt(2 * nu / f)
    u = U * (1 - sp.exp(-z / d) * sp.cos(z / d))                                # (13.41a)
    v = U * sp.exp(-z / d) * sp.sin(z / d)                                      # (13.41b)
    assert sp.simplify(-f * v - nu * u.diff(z, 2)) == 0                         # (13.33)
    assert sp.simplify(f * u - nu * v.diff(z, 2) - f * U) == 0                  # (13.34)
    assert u.subs(z, 0) == 0 and v.subs(z, 0) == 0 and sp.limit(u, z, sp.oo) == U and sp.limit(v, z, sp.oo) == 0
    assert sp.simplify(sp.integrate(v, (z, 0, sp.oo)) - U * d / 2) == 0         # transport to the left of the interior flow
    assert sp.limit(v / u, z, 0) == 1                                           # 45 degrees at the surface
    s = sp.Symbol("s", positive=True)                                           # s = z/delta
    us = 1 - sp.exp(-s) * sp.cos(s)
    assert sp.simplify(us.diff(s).subs(s, 3 * sp.pi / 4)) == 0 and us.diff(s, 2).subs(s, 3 * sp.pi / 4) < 0
    umax = us.subs(s, 3 * sp.pi / 4)
    assert sp.simplify(umax - (1 + sp.exp(-3 * sp.pi / 4) / sp.sqrt(2))) == 0
    assert abs(float(umax) - 1.0670) < 5e-5 and abs(float(us.subs(s, sp.pi)) - (1 + math.exp(-PI))) < 1e-15
    assert float(umax) > float(us.subs(s, sp.pi))                               # the largest overshoot is below z = pi delta


def test_ekman_bottom_V1_overshoot_transport_angle_forces_hemispheres():  # V1 + V4 + V7, C06, N47, N49
    for f in HEMI:
        s = np.sign(f)
        nu, Ug = I["nu_v_atm"], I["U_g"]
        d = float(ch13.ekman_depth(nu, f))
        zz = np.linspace(0.0, 10 * d, 200001)
        u, v = ch13.ekman_bottom(zz, Ug, 0.0, nu, f)
        j = int(np.argmax(u))
        assert abs(zz[j] / d - 3 * PI / 4) < 1e-4 and abs(u[j] / Ug - (1 + math.exp(-3 * PI / 4) / math.sqrt(2))) < 1e-9
        assert abs(u[j] / Ug - 1.0670) < 5e-5 and u[0] == 0.0 and v[0] == 0.0
        assert math.isclose(ch13.ekman_bottom(PI * d, Ug, 0.0, nu, f)[0] / Ug, 1 + math.exp(-PI), rel_tol=1e-12)
        assert np.all(s * v[1:2000] > 0)                                         # toward low pressure: left of U for f > 0
        rx, ry = ch13.ekman_residual(zz[::100], u[::100], v[::100], nu, f, U_g=Ug)
        assert max(np.max(np.abs(rx)), np.max(np.abs(ry))) < 1e-5 * abs(f) * Ug
        Mx, My = ch13.ekman_bottom_transport(Ug, 0.0, nu, f)
        q = quad(lambda t: ch13.ekman_bottom(t, Ug, 0.0, nu, f)[1], 0, 40 * d, epsabs=1e-12, limit=400)[0]
        qd = quad(lambda t: ch13.ekman_bottom(t, Ug, 0.0, nu, f)[0] - Ug, 0, 40 * d, epsabs=1e-12, limit=400)[0]
        assert math.isclose(My, s * 0.5 * Ug * d, rel_tol=1e-12) and math.isclose(Mx, -0.5 * Ug * d, rel_tol=1e-12)
        assert abs(q - My) < 1e-8 * abs(My) and abs(qd - Mx) < 1e-8 * abs(Mx)
        fb = ch13.ekman_force_balance(zz[::2000], Ug, nu, f, rho=1.2)
        assert np.max(np.abs(fb["sum"][0])) < 1e-18 and np.max(np.abs(fb["sum"][1])) < 1e-15   # three forces close
        assert math.isclose(fb["angle_to_isobars"][0], PI / 4) and fb["angle_to_isobars"][-1] < 1e-3
        assert close(fb["pressure"][1], 1.2 * f * Ug) and close(fb["coriolis"][1], -1.2 * f * fb["u"])
        d2u = np.gradient(np.gradient(u, zz), zz)[2000:-2000:2000]
        assert close(fb["friction"][0][1:-1], 1.2 * nu * d2u, atol=2e-4 * 1.2 * abs(f) * Ug)   # friction = rho nu u''
        assert close(fb["speed"], np.hypot(fb["u"], fb["v"]))
        Vt = ch13.ekman_bottom(zz[:50], 3.0, -4.0, nu, f, as_complex=True)       # any direction of the interior flow
        assert close(Vt, (3.0 - 4.0j) / Ug * ch13.ekman_bottom(zz[:50], Ug, 0.0, nu, f, as_complex=True))
    assert math.isclose(ch13.eddy_viscosity_from_depth(730.0, F_S), 0.5 * F_N * 730.0 ** 2)
    assert math.isclose(ch13.ekman_depth(ch13.eddy_viscosity_from_depth(730.0, F_S), F_S), 730.0)


def test_ekman_solve_V3_second_order_and_limits():  # V3 + V1, N37 (ours): d/dz(K dV/dz) = i f (V - V_g)
    f, rho, tau, nu = F_S, 1027.0, 0.07 + 0.02j, 0.03
    d = float(ch13.ekman_depth(nu, f))
    hs, e_const, e_var = [], [], []
    zr = np.linspace(-12 * d, 0.0, 5121)
    Kfn = lambda s: nu * (1 + 0.5 * np.exp(s / d))  # noqa: E731
    ref = ch13.ekman_solve(zr, Kfn, f, tau=tau, rho=rho)
    for n in (81, 161, 321, 641):
        zg = np.linspace(-12 * d, 0.0, n)
        V = ch13.ekman_solve(zg, nu, f, tau=tau, rho=rho)
        ex = ch13.ekman_surface(zg, tau.real, tau.imag, rho, nu, f, as_complex=True)
        e_const.append(np.max(np.abs(V - ex)) / np.max(np.abs(ex)))
        Vv = ch13.ekman_solve(zg, Kfn, f, tau=tau, rho=rho)
        e_var.append(np.max(np.abs(Vv - ref[::(5120 // (n - 1))])) / np.max(np.abs(ref)))
        hs.append(zg[1] - zg[0])
    assert abs(observed_order(hs, e_const) - 2.0) < 0.15 and e_const[-1] < 1e-4
    assert abs(observed_order(hs[:3], e_var[:3]) - 2.0) < 0.15
    # the transport is set by the stress alone, whatever K(z) is (first integral of the equation)
    w = np.zeros_like(zr)
    w[:-1] += 0.5 * np.diff(zr)
    w[1:] += 0.5 * np.diff(zr)
    M = np.sum(ref * w)
    assert abs(M - (-1j * tau / (rho * f))) < 2e-4 * abs(tau / (rho * f))
    # bottom layer and finite depth
    zb = np.linspace(0.0, 14 * d, 1401)
    Vb = ch13.ekman_solve(zb, nu, f, V_g=0.3 + 0.1j, bottom="noslip", top="geostrophic")
    assert np.max(np.abs(Vb - ch13.ekman_bottom(zb, 0.3, 0.1, nu, f, as_complex=True))) < 2e-4 * abs(0.3 + 0.1j)
    H = 5 * d
    zf = np.linspace(-H, 0.0, 1001)
    Vf = ch13.ekman_finite_depth(zf, tau.real, tau.imag, 0.05, 0.0, H, nu, rho, f)
    Vn = ch13.ekman_solve(zf, nu, f, tau=tau, rho=rho, V_g=0.05, bottom="noslip")
    assert np.max(np.abs(Vf - Vn)) < 1e-4 * np.max(np.abs(Vf)) and abs(Vf[0]) < 1e-15
    assert abs(rho * nu * (3 * Vf[-1] - 4 * Vf[-2] + Vf[-3]) / (2 * (zf[1] - zf[0])) - tau) < 1e-4 * abs(tau)
    deep = ch13.ekman_finite_depth(np.array([-2 * d, 0.0]), tau.real, tau.imag, 0.0, 0.0, 40 * d, nu, rho, f)
    assert close(deep, ch13.ekman_surface(np.array([-2 * d, 0.0]), tau.real, tau.imag, rho, nu, f, as_complex=True), rtol=1e-12)
    for bad in (dict(K=np.full(81, nu)), dict(K=-nu)):
        with pytest.raises(ValueError):
            ch13.ekman_solve(np.linspace(-1.0, 0.0, 81), bad["K"], f, tau=tau, rho=rho)
    with pytest.raises(ValueError):
        ch13.ekman_solve(np.linspace(-1.0, 0.0, 81), nu, f)                      # a stress condition needs tau and rho


# =====================================================================================================================
# 4. §13.8–§13.9 — shallow-water equations and vertical normal modes (C07, C08)
# =====================================================================================================================
def test_shallow_water_set_V2_derivation():  # V2, D10 (13.44)–(13.45)
    x, y, z, t, eps = sp.symbols("x y z t epsilon", real=True)
    H, g, f = sp.symbols("H g f", positive=True)
    u, v, eta = (sp.Function(n)(x, y, t) for n in ("u", "v", "eta"))
    w = -(u.diff(x) + v.diff(y)) * z                                           # continuity with u, v independent of z, w(0) = 0
    assert sp.simplify(u.diff(x) + v.diff(y) + w.diff(z)) == 0
    kin = eta.diff(t) + u * eta.diff(x) + v * eta.diff(y) - w.subs(z, H + eta)   # w = D(eta)/Dt at the surface
    flux = eta.diff(t) + (u * (H + eta)).diff(x) + (v * (H + eta)).diff(y)       # (13.44)
    assert sp.expand(kin - flux) == 0
    p = sp.Symbol("rho") * g * (H + eta - z)                                    # hydrostatic pressure
    assert sp.simplify(-p.diff(x) / sp.Symbol("rho") + g * eta.diff(x)) == 0 and p.diff(x).diff(z) == 0
    sub = {u: eps * u, v: eps * v, eta: eps * eta}
    lin = sp.expand(flux.subs(sub, simultaneous=True).doit()).coeff(eps, 1)
    assert sp.expand(lin - (eta.diff(t) + H * (u.diff(x) + v.diff(y)))) == 0    # (13.45a)
    mom = u.diff(t) + u * u.diff(x) + v * u.diff(y) - f * v + g * eta.diff(x)
    assert sp.expand(sp.expand(mom.subs(sub, simultaneous=True).doit()).coeff(eps, 1) - (u.diff(t) - f * v + g * eta.diff(x))) == 0
    assert math.isclose(ch13.long_wave_speed(I["ocean_H"]), math.sqrt(G0 * I["ocean_H"]))
    assert math.isclose(ch13.equivalent_depth(ch13.long_wave_speed(37.0)), 37.0)


def test_vertical_modes_V2_derivation():  # V2, D11 (13.52)–(13.69); orthogonality with weight 1 for BOTH lids (claim 2)
    z = sp.Symbol("z", real=True)
    N, H, g, c, cm, cn = sp.symbols("N H g c c_m c_n", positive=True)
    # (a) uniform N with a free surface: (13.67)–(13.68) solve (13.56) and (13.65); (13.64) gives (13.69)
    psi = sp.cos(N * z / c) - c * N / g * sp.sin(N * z / c)
    assert sp.simplify((psi.diff(z) / N ** 2).diff(z) + psi / c ** 2) == 0      # (13.56)
    assert sp.simplify((psi.diff(z) + N ** 2 / g * psi).subs(z, 0)) == 0        # (13.65)
    bottom = sp.simplify(psi.diff(z).subs(z, -H) * c / (N * sp.cos(N * H / c)))
    assert sp.simplify(bottom - (sp.tan(N * H / c) - c * N / g)) == 0           # (13.64) -> tan(NH/c) = cN/g (13.69)
    W = sp.integrate(psi, (z, -H, z))                                           # (13.53): w's shape is the integral of psi
    assert W.subs(z, -H) == 0 and sp.simplify(W.diff(z) - psi) == 0             # w = 0 at the bottom
    assert sp.simplify((psi.diff(z) / N ** 2 + W / c ** 2).diff(z)) == 0        # (13.55): psi'/N^2 + W/c^2 is constant in z,
    assert sp.simplify((psi.diff(z) / N ** 2 + W / c ** 2).subs(z, -H) - psi.diff(z).subs(z, -H) / N ** 2) == 0   # = 0 by (13.64)
    # (b) Lagrange identity for arbitrary N^2(z): d/dz[(psi_n psi_m' - psi_m psi_n')/N^2] = (1/c_n^2 - 1/c_m^2) psi_m psi_n
    N2 = sp.Function("N2")(z)
    pm, pn = sp.Function("psi_m")(z), sp.Function("psi_n")(z)
    br = (pn * pm.diff(z) - pm * pn.diff(z)) / N2
    ode = {pm.diff(z, 2): N2.diff(z) / N2 * pm.diff(z) - N2 / cm ** 2 * pm, pn.diff(z, 2): N2.diff(z) / N2 * pn.diff(z) - N2 / cn ** 2 * pn}
    assert sp.simplify(br.diff(z).subs(ode) - (1 / cn ** 2 - 1 / cm ** 2) * pm * pn) == 0
    # (c) the bracket vanishes at the bottom (13.64) and at the top for a lid AND for the free surface (13.65)
    a, b, n2 = sp.symbols("a b n2", positive=True)                              # psi_m(0), psi_n(0), N^2(0)
    bracket = lambda dm, dn: (b * dm - a * dn) / n2  # noqa: E731
    assert bracket(0, 0) == 0                                                   # rigid lid or bottom: psi' = 0
    assert sp.simplify(bracket(-n2 / g * a, -n2 / g * b)) == 0                  # free surface: psi' = -(N^2/g) psi, both modes
    # so int psi_m psi_n dz = 0 with weight 1; the surface term psi_m(0) psi_n(0)/g belongs to the energy relation only
    # energy relation: psi_m (psi_n'/N^2)' = (psi_m psi_n'/N^2)' - psi_m' psi_n'/N^2, whose top value is -psi_m psi_n/g
    assert sp.simplify((pm * pn.diff(z) / N2).diff(z) - pm.diff(z) * pn.diff(z) / N2 - pm * (pn.diff(z) / N2).diff(z)) == 0
    assert sp.simplify(a * (-n2 / g * b) / n2 + a * b / g) == 0


def test_vertical_modes_V1_orthogonality_weight_one_exact_quadrature():  # V1 + V4, C08 (N60): claim 2, numerically
    N, H = I["ocean_N"], I["ocean_H"]
    ex = VM.modes_uniform_N(N, H, n_modes=5, lid="free", nz=401)
    X = N * H / ex.c
    assert np.max(np.abs(np.tan(X) - ex.c * N / G0)) < 1e-12                    # (13.69)
    assert np.all((X > np.arange(5) * PI) & (X < np.arange(5) * PI + PI / 2))   # one root per branch, slightly above n pi
    psi = lambda n, s: np.cos(N * s / ex.c[n]) - ex.c[n] * N / G0 * np.sin(N * s / ex.c[n])  # noqa: E731
    for m, n in ((0, 1), (0, 3), (1, 2), (2, 4)):
        ov = quad(lambda s: psi(m, s) * psi(n, s), -H, 0, epsabs=1e-13, epsrel=1e-13, limit=400)[0]
        norm = math.sqrt(quad(lambda s: psi(m, s) ** 2, -H, 0)[0] * quad(lambda s: psi(n, s) ** 2, -H, 0)[0])
        assert abs(ov) / norm < 1e-12                                           # weight 1, free surface: orthogonal
        with_surface_term = ov + psi(m, 0.0) * psi(n, 0.0) * ex.c[m] * ex.c[n] / G0   # a wrongly added surface term (scaled c_m c_n/g)
        assert abs(with_surface_term) / norm > 1e-6                             # ... would break it: the alternative is refuted
        lhs = quad(lambda s: psi(m, s) * psi(m, s), -H, 0)[0] / ex.c[m] ** 2    # energy relation: here the surface term belongs
        dps = lambda s: (-N / ex.c[m] * np.sin(N * s / ex.c[m]) - N * N / G0 * np.cos(N * s / ex.c[m]))  # noqa: E731
        rhs = quad(lambda s: dps(s) ** 2 / N ** 2, -H, 0)[0] + psi(m, 0.0) ** 2 / G0
        assert abs(lhs - rhs) < 1e-10 * abs(lhs)
    z = np.linspace(-H, 0.0, 601)
    for N2 in (N ** 2, ch13.thermocline_N2(z)):
        for lid in ("free", "rigid"):
            md = ch13.vertical_modes(z, N2, n_modes=5, lid=lid)
            assert np.max(np.abs(ch13.orthogonality_matrix(md) - np.eye(5))) < 1e-12
            En = ch13.orthogonality_matrix(md, kind="energy", N2=N2 * np.ones_like(z), lid=lid)
            assert np.max(np.abs(En - np.eye(5))) < 1e-8
            assert close(md.psi[:, -1], 1.0) and np.all(np.diff(md.c) < 0) and close(md.He, md.c ** 2 / G0)
            assert isinstance(md, ch13.Modes) and md[0][1] == md.c[1] and np.isclose(md.weights.sum(), H)
    wrong = ch13.orthogonality_matrix(ch13.vertical_modes(z, N ** 2, n_modes=4, lid="free"), kind="energy", N2=N ** 2, lid="rigid")
    assert np.max(np.abs(wrong - np.eye(4))) > 1e-6                             # dropping the surface term from the energy form fails


def test_vertical_modes_V3_finite_volume_second_order_and_cross_methods():  # V3 + V1, C08 (13.56), (13.64), (13.65)
    N, H = I["ocean_N"], I["ocean_H"]
    orders = {}
    for lid in ("free", "rigid"):
        hs, errs, epsi = [], [], []
        for nz in (41, 81, 161, 321):
            z = np.linspace(-H, 0.0, nz)
            md = ch13.vertical_modes(z, N ** 2, n_modes=4, lid=lid)
            ex = ch13.modes_uniform_N(N, H, n_modes=4, lid=lid, nz=nz)
            errs.append(np.abs(md.c / ex.c - 1))
            epsi.append(np.max(np.abs(md.psi - ex.psi)))
            hs.append(H / (nz - 1))
        errs = np.array(errs)
        first_bc = 1 if lid == "free" else 0
        orders[lid] = [observed_order(hs, errs[:, j]) for j in range(first_bc, 4)]
        assert all(abs(o - 2.0) < 0.15 for o in orders[lid]), orders
        assert max(epsi) < 5e-6 and (epsi[-1] < 1e-10 or abs(observed_order(hs, epsi) - 2.0) < 0.15)   # shapes: exact (lid) or order 2
        if lid == "free":
            assert errs[-1, 0] < 1e-7                                           # barotropic speed: already at 1e-8
    ex = ch13.modes_uniform_N(N, H, n_modes=4, lid="rigid")
    assert close(ex.c, [ch13.baroclinic_mode_speed(N, H, n) for n in (1, 2, 3, 4)], rtol=1e-14)   # (13.71)
    assert close(ex.psi[1], np.cos(2 * PI * ex.z / H))
    fr = ch13.modes_uniform_N(N, H, n_modes=3)
    assert abs(fr.c[0] / math.sqrt(G0 * H) - 1) < 1e-3 and abs(fr.c[1] / ch13.baroclinic_mode_speed(N, H, 1) - 1) < 1e-3
    rl = ch13.rigid_lid_error(N, H)
    assert math.isclose(rl, fr.c[1] / ex.c[0] - 1, rel_tol=1e-9) and math.isclose(rl, -N * N * H / (G0 * PI ** 2), rel_tol=2e-3)
    # three independent solvers on a thermocline profile: finite volume, Chebyshev collocation, RK4 shooting
    z = np.linspace(-H, 0.0, 1601)
    Nfun = lambda s: ch13.thermocline_N2(s)  # noqa: E731
    fd = ch13.vertical_modes(z, Nfun(z), n_modes=3)
    cb = ch13.vertical_modes(z, Nfun, n_modes=3, method="cheb", n_cheb=96)
    for n in (0, 1, 2):
        cs, ps = ch13.vertical_modes_shooting(z[::4], lambda s: float(Nfun(s)), n=n)
        assert abs(fd.c[n] / cs - 1) < 2.5e-4                                   # second-order truncation at 1601 nodes
        if n >= 1:
            assert abs(cb.c[n] / cs - 1) < 1e-6 and np.max(np.abs(np.interp(z[::4], z, cb.psi[n]) - ps)) < 2e-4
    c0s, _ = ch13.vertical_modes_shooting(z[::4], lambda s: float(Nfun(s)), n=0)
    assert abs(cb.c[0] / c0s - 1) < 2e-4          # measured 8.5e-5: the non-symmetric Chebyshev solve loses the barotropic digits
    cr, pr = ch13.vertical_modes_shooting(np.linspace(-H, 0.0, 401), lambda s: N ** 2, n=2, lid="rigid")
    assert abs(cr / ch13.baroclinic_mode_speed(N, H, 2) - 1) < 1e-8 and abs(pr[0] - 1.0) < 1e-7
    with pytest.raises(ValueError):
        ch13.vertical_modes_shooting(z, lambda s: N ** 2, n=0, lid="rigid")
    with pytest.raises(ValueError):
        ch13.vertical_modes(z, -Nfun(z))


def test_vertical_modes_V1_structures_projection_and_amplitudes():  # V1, C08 (13.52)–(13.54), (13.60)–(13.62)
    N, H = I["ocean_N"], I["ocean_H"]
    md = ch13.modes_uniform_N(N, H, n_modes=4, lid="rigid", nz=801)
    for n in range(4):
        k = (n + 1) * PI / H
        assert np.max(np.abs(ch13.w_structure(md, n) - np.sin(k * md.z) / k)) < 2e-5 * H           # integral of psi, zero at both ends
        assert np.max(np.abs(ch13.rho_structure(md, n) + k * np.sin(k * md.z))) < 2e-4 * k          # derivative of psi
    z = np.linspace(-H, 0.0, 801)
    mt = ch13.vertical_modes(z, ch13.thermocline_N2(z), n_modes=4)
    coeffs = np.array([0.3, -1.1, 0.0, 0.45])
    prof = ch13.reconstruct(mt, coeffs)
    assert close(ch13.project(mt, prof), coeffs, atol=1e-12) and close(prof, coeffs @ mt.psi)
    one = ch13.project(mt, mt.psi[2])
    assert close(one, [0, 0, 1, 0], atol=1e-12)
    am = ch13.modal_amplitudes(2.7, 4e-6, p_n=0.31, rho0=1027.0)
    assert math.isclose(am["w_n"], 4e-6 / 2.7 ** 2) and math.isclose(am["rho_n"], -1027.0 * 0.31 / G0)
    assert (Q_(1.0, "m**2/s**3") / Q_(1.0, "m/s") ** 2).check("1/[time]")       # w_n = (1/c^2) dp_n/dt is 1/s (trap T12)
    assert math.isclose(ch13.wkb_mode_speed(z, np.full_like(z, N), 2), ch13.baroclinic_mode_speed(N, H, 2), rel_tol=1e-12)
    wk = [abs(ch13.wkb_mode_speed(z, np.sqrt(ch13.thermocline_N2(z)), n) / ch13.vertical_modes(z, ch13.thermocline_N2(z), n_modes=7, lid="rigid").c[n - 1] - 1) for n in (1, 3, 6)]
    assert wk[2] < wk[0] and wk[2] < 0.1                                        # WKB improves with the mode number
    assert ch13.hydrostatic_linear_set_sympy()["ok"]


def test_vertical_modes_V7_first_root_is_small_printed_statement_fails():  # V7, slip #6 (after (13.69))
    N, H = I["ocean_N"], I["ocean_H"]
    X0 = N * H / ch13.modes_uniform_N(N, H, n_modes=1).c[0]
    assert X0 < 0.1 and abs(X0 - N * math.sqrt(H / G0)) < 1e-3 * X0             # the root is N sqrt(H/g): a few hundredths
    resid = lambda X: math.tan(X) - (N * N * H / G0) / X  # noqa: E731         (13.69) written in X = NH/c
    assert abs(resid(X0)) < 1e-10 and abs(resid(1.0)) > 1.0                     # X = 1 is not a root: the printed "= 1" fails
    assert abs(math.tan(X0) / X0 - 1) < 2e-3 and abs(math.tan(1.0) / 1.0 - 1) > 0.5   # tan x ~ x holds at X0, not at 1


# =====================================================================================================================
# 5. §13.10–§13.12 — the dispersion cubic, Poincaré and Kelvin waves, the Rossby radius, adjustment (C09–C12)
# =====================================================================================================================
def test_v_equation_V2_derivation():  # V2, D12 (13.72)–(13.75), the book's own route with f -> f0 except in df/dy
    x, y, t = sp.symbols("x y t", real=True)
    g, H, f0, beta = sp.symbols("g H f_0 beta", positive=True)
    u, v, eta = (sp.Function(n)(x, y, t) for n in ("u", "v", "eta"))
    D = u.diff(x) + v.diff(y)
    E1, E2, E3 = u.diff(t) - f0 * v + g * eta.diff(x), v.diff(t) + f0 * u + g * eta.diff(y), eta.diff(t) + H * D   # (13.45)
    e72 = u.diff(t, 2) - f0 * v.diff(t) - g * H * D.diff(x)                    # (13.72)
    e73 = v.diff(t, 2) + f0 * u.diff(t) - g * H * D.diff(y)                    # (13.73)
    assert sp.expand(e72 - (E1.diff(t) - g * E3.diff(x))) == 0 and sp.expand(e73 - (E2.diff(t) - g * E3.diff(y))) == 0
    e74 = v.diff(t, 3) + f0 * (f0 * v.diff(t) + g * H * D.diff(x)) - g * H * D.diff(y, t)   # (13.74)
    assert sp.expand(e74 - (e73.diff(t) - f0 * e72)) == 0
    vort = (u.diff(y) - v.diff(x)).diff(t) - f0 * D - beta * v                 # the vorticity equation of the page
    fy = f0 + beta * y
    exact_vort = (u.diff(t) - fy * v + g * eta.diff(x)).diff(y) - (v.diff(t) + fy * u + g * eta.diff(y)).diff(x)
    assert sp.expand(exact_vort - vort.subs(f0, fy)) == 0                      # where beta v comes from: d(f u)/dy
    e75 = v.diff(t, 3) - g * H * (v.diff(x, 2) + v.diff(y, 2)).diff(t) + f0 ** 2 * v.diff(t) - g * H * beta * v.diff(x)
    assert sp.expand(e75 - (e74 + g * H * vort.diff(x))) == 0                  # (13.75)
    eng = ch13.v_equation_sympy()
    assert eng["ok"] and sp.expand(eng["neglected_relative"] - ((f0 + beta * y) ** 2 / f0 ** 2 - 1)) == 0


def test_dispersion_cubic_V2_derivation():  # V2, D13 (13.76): plane wave in (13.75), discriminant, regimes
    x, y, t = sp.symbols("x y t", real=True)
    k, l, om = sp.symbols("k l omega", real=True)
    c, f0, beta = sp.symbols("c f_0 beta", positive=True)
    v = sp.exp(sp.I * (k * x + l * y - om * t))
    e75 = v.diff(t, 3) - c ** 2 * (v.diff(x, 2) + v.diff(y, 2)).diff(t) + f0 ** 2 * v.diff(t) - c ** 2 * beta * v.diff(x)
    cubic = om ** 3 - c ** 2 * om * (k ** 2 + l ** 2) - f0 ** 2 * om - c ** 2 * beta * k    # (13.76)
    assert sp.simplify(e75 / v - sp.I * cubic) == 0
    P, q = sp.symbols("P q", positive=True)
    assert sp.expand(sp.discriminant(om ** 3 - P * om - q, om) - (4 * P ** 3 - 27 * q ** 2)) == 0
    r1, r2, r3 = sp.symbols("r1 r2 r3")
    poly = sp.Poly(sp.expand((om - r1) * (om - r2) * (om - r3)), om)
    assert poly.coeff_monomial(om ** 2) == -(r1 + r2 + r3) and poly.coeff_monomial(1) == -r1 * r2 * r3   # sum = 0, product = +q
    # with the wrong sign of the beta term the slow wave would travel east
    slow = sp.solve(cubic.subs(om ** 3, 0), om)[0]
    assert sp.simplify(slow + beta * k / (k ** 2 + l ** 2 + f0 ** 2 / c ** 2)) == 0      # the low-frequency root is (13.118)
    fast = sp.solve((cubic + c ** 2 * beta * k), om)
    assert any(sp.simplify(s ** 2 - (f0 ** 2 + c ** 2 * (k ** 2 + l ** 2))) == 0 for s in fast if s != 0)   # (13.82)


def test_shallow_water_omega_V1_roots_against_numpy_and_invariants():  # V1 + V4, C09 (13.76): claim 1
    worst_res, worst_sum, worst_np, n_real, n_nan = 0.0, 0.0, 0.0, 0, 0
    for lat in np.deg2rad([4.0, 12.0, 35.0, 60.0, -35.0, -76.0]):
        f0, beta = float(ch13.coriolis_parameter(lat)), float(ch13.beta_parameter(lat))
        for c in (203.0, 2.9, 0.41):
            for lam in np.geomspace(3e4, 4e7, 9):
                for ang in (0.0, 0.7, 2.5, PI):
                    k, l = 2 * PI / lam * math.cos(ang), 2 * PI / lam * math.sin(ang)
                    w = np.array(ch13.shallow_water_omega(k, l, c, f0, beta))
                    P, q = c ** 2 * (k * k + l * l) + f0 ** 2, c ** 2 * beta * k
                    disc = ch13.shallow_water_discriminant(k, l, c, f0, beta)
                    assert math.isclose(disc, 4 * P ** 3 - 27 * q ** 2, rel_tol=1e-12)
                    if disc < 0:
                        assert np.all(np.isnan(w))
                        n_nan += 1
                        continue
                    assert np.all(np.isfinite(w)) and w[0] <= w[1] <= w[2]
                    n_real += 1
                    worst_res = max(worst_res, float(np.max(np.abs(w ** 3 - P * w - q) / (np.abs(w) ** 3 + P * np.abs(w) + abs(q)))))
                    worst_sum = max(worst_sum, abs(w.sum()) / np.abs(w).max())
                    assert abs(np.prod(w) - q) <= 1e-9 * abs(q) + 1e-300        # product of the roots = +c^2 beta k
                    assert abs(w[0] * w[1] + w[0] * w[2] + w[1] * w[2] + P) < 1e-9 * P
                    r = np.sort(np.roots([1.0, 0.0, -P, -q]).real)
                    worst_np = max(worst_np, float(np.max(np.abs(w[[0, 2]] - r[[0, 2]]) / np.abs(r[[0, 2]]))))
                    assert abs(w[1] - r[1]) <= 1e-6 * abs(r[1]) + 1e-9 * abs(r[2])   # numpy's slow root carries the round-off
                    br = ch13.shallow_water_branches(k, l, c, f0, beta)
                    assert (br["poincare_minus"], br["rossby"], br["poincare_plus"]) == tuple(w) and br["kelvin"] == c * k
                    ts = ch13.dispersion_term_sizes(k, l, c, f0, beta, w[1])
                    assert abs(ts["sum"]) < 1e-9 * (abs(ts["beta"]) + abs(ts["rotation"]) + 1e-300)
    assert n_real >= 600 and worst_res < 5e-15 and worst_sum < 1e-12 and worst_np < 1e-9
    # NaN exactly where the discriminant is negative (planetary scales near the equator, outside beta-plane scaling)
    kk = np.geomspace(1e-9, 1e-4, 400)
    disc = ch13.shallow_water_discriminant(kk, 0.0, 300.0, 1e-6, 2.2e-11)
    w = np.array(ch13.shallow_water_omega(kk, 0.0, 300.0, 1e-6, 2.2e-11))
    assert np.any(disc < 0) and np.any(disc > 0) and np.array_equal(np.isnan(w[0]), disc < 0) and np.array_equal(np.isnan(w[1]), disc < 0)
    for kneg in kk[disc < 0][::25]:
        assert np.sum(np.abs(np.roots([1.0, 0.0, -(300.0 ** 2 * kneg ** 2 + 1e-12), -(300.0 ** 2 * 2.2e-11 * kneg)]).imag) > 0) == 2


def test_shallow_water_omega_V7_limits_poincare_kelvin_rossby_and_slip7():  # V7, C09 (N83); slip #7
    f0, beta = F_N, float(ch13.beta_parameter(LAT))
    c = float(ch13.baroclinic_mode_speed(I["ocean_N"], I["ocean_H"], 1))
    Lam = c / f0
    K = np.geomspace(0.05, 300.0, 60) / Lam
    wm, wr, wp = (np.array(a) for a in ch13.shallow_water_omega(-K, 0.0, c, f0, beta))
    assert close(wp, ch13.poincare_omega(K, f0, c), rtol=2e-3) and close(wm, -ch13.poincare_omega(K, f0, c), rtol=2e-3)
    assert close(wr, ch13.rossby_omega(-K, 0.0, beta, f0, c), rtol=2e-3) and np.all(wr > 0)       # westward: k < 0, omega > 0
    short = K * Lam > 50
    assert close(wp[short], ch13.kelvin_omega(K[short], c), rtol=3e-4)                            # omega -> cK: non-rotating
    assert np.all(np.abs(wp) > f0) and np.all(np.abs(wm) > f0) and np.all(np.abs(wr) < 0.05 * f0)       # two super-, one sub-inertial
    w_b0 = ch13.shallow_water_omega(K, 0.3 * K, c, f0, 0.0)
    assert close(w_b0[1], 0.0, atol=1e-18) and close(w_b0[2], ch13.poincare_omega(np.hypot(K, 0.3 * K), f0, c), rtol=1e-13)
    w_f0 = ch13.shallow_water_omega(K, 0.0, c, 0.0, 0.0)
    assert close(w_f0[2], c * K, rtol=1e-13)                                                     # f = beta = 0: ch07 long waves
    assert close(w_f0[2] / K, float(WAV.phase_speed(1e-9, H=c * c / WAV.G_BOOK, g=WAV.G_BOOK)), rtol=1e-6)
    for fs in HEMI:                                                                                   # the roots depend on f0^2 only
        assert close(ch13.shallow_water_omega(-K, 0.0, c, fs, beta), (wm, wr, wp), rtol=1e-14)
    # slip #7: on the slow root (omega << f) omega^3 is the SMALLEST term; at omega = 3f it is the LARGEST
    k1 = -1.0 / Lam
    slow = ch13.shallow_water_omega(k1, 0.0, c, f0, beta)[1]
    ts = {n: abs(v) for n, v in ch13.dispersion_term_sizes(k1, 0.0, c, f0, beta, slow).items() if n != "sum"}
    assert abs(slow) < 0.02 * f0 and min(ts, key=ts.get) == "omega3" and ts["omega3"] < 1e-3 * ts["beta"]
    k3 = math.sqrt(8.0) / Lam                                           # where omega = 3 f on the fast branch
    fastw = ch13.shallow_water_omega(k3, 0.0, c, f0, beta)[2]
    tf = {n: abs(v) for n, v in ch13.dispersion_term_sizes(k3, 0.0, c, f0, beta, fastw).items() if n != "sum"}
    assert math.isclose(fastw / f0, 3.0, rel_tol=1e-3) and max(tf, key=tf.get) == "omega3" and tf["omega3"] > 500 * tf["beta"]
    assert ch13.shallow_water_regime(slow, f0) == dict(ratio=abs(slow) / f0, regime="low", neglect="omega^3")
    assert ch13.shallow_water_regime(3.5 * f0, -f0)["neglect"] == "f0^2 omega" and ch13.shallow_water_regime(2 * f0, f0)["regime"] == "near-inertial"


def test_poincare_V2_derivation():  # V2, D14 (13.77)–(13.82) and the group velocity
    k, l, om, f, g, H = sp.symbols("k l omega f g H", positive=True)
    uh, vh, eh = sp.symbols("uhat vhat etahat")
    sol = sp.solve([-sp.I * om * uh - f * vh + sp.I * k * g * eh, -sp.I * om * vh + f * uh + sp.I * l * g * eh], [uh, vh])   # (13.77), (13.78)
    assert sp.simplify(sol[uh] - g * eh / (om ** 2 - f ** 2) * (om * k + sp.I * f * l)) == 0             # (13.80a)
    assert sp.simplify(sol[vh] - g * eh / (om ** 2 - f ** 2) * (-sp.I * f * k + om * l)) == 0            # (13.80b)
    e79 = sp.simplify((-sp.I * om * eh + sp.I * H * (k * sol[uh] + l * sol[vh])) * (om ** 2 - f ** 2) / (sp.I * eh))
    assert sp.expand(e79 + om * (om ** 2 - f ** 2 - g * H * (k ** 2 + l ** 2))) == 0                     # (13.82) (and omega = 0)
    w = sp.sqrt(f ** 2 + g * H * (k ** 2 + l ** 2))
    assert sp.simplify(w.diff(k) - g * H * k / w) == 0 and sp.simplify(w.diff(l) - g * H * l / w) == 0   # c_g = c^2 K/omega
    assert sp.simplify(sp.sqrt(w.diff(k) ** 2 + w.diff(l) ** 2) * w / sp.sqrt(k ** 2 + l ** 2) - g * H) == 0   # c_p c_g = c^2


def test_poincare_V1_fields_satisfy_the_linear_set_and_orbits():  # V1 + V7, C10 (13.45), (13.80), (13.83)
    H = I["ocean_H"]
    c = math.sqrt(G0 * H)
    for f in HEMI:
        k, l = 2 * PI / 2.3e6, -2 * PI / 4.1e6
        w = float(ch13.poincare_omega(math.hypot(k, l), f, c))
        h, ht = 50.0, 2.0
        x0, y0, t0 = 3.1e5, -2.2e5, 1234.0
        F = lambda x, y, t: ch13.poincare_fields(x, y, t, k, l, 0.4, H, f)  # noqa: E731
        d = lambda i, var: ((F(x0 + h, y0, t0)[i] - F(x0 - h, y0, t0)[i]) / (2 * h) if var == "x" else  # noqa: E731
                            (F(x0, y0 + h, t0)[i] - F(x0, y0 - h, t0)[i]) / (2 * h) if var == "y" else
                            (F(x0, y0, t0 + ht)[i] - F(x0, y0, t0 - ht)[i]) / (2 * ht))
        e, u, v = F(x0, y0, t0)
        amp_u = G0 * 0.4 * math.hypot(k, l) / w
        fd = 2 * (w * ht) ** 2 / 6                                                           # truncation of the central difference
        assert abs(d(0, "t") + H * (d(1, "x") + d(2, "y"))) < fd * 0.4 * w                   # (13.45a)
        assert abs(d(1, "t") - f * v + G0 * d(0, "x")) < fd * amp_u * w                      # (13.45b)
        assert abs(d(2, "t") + f * u + G0 * d(0, "y")) < fd * amp_u * w                      # (13.45c)
        cg = ch13.poincare_group_velocity(k, l, f, c)
        dk = 1e-12
        fdx = (ch13.poincare_omega(math.hypot(k + dk, l), f, c) - ch13.poincare_omega(math.hypot(k - dk, l), f, c)) / (2 * dk)
        assert math.isclose(cg[0], fdx, rel_tol=1e-5) and math.isclose(math.hypot(*cg) * w / math.hypot(k, l), c * c, rel_tol=1e-12)
        assert math.hypot(*cg) < c < w / math.hypot(k, l)
        uh, vh = ch13.poincare_amplitudes(k, l, w, f, G0, 0.4)
        assert abs(-1j * w * 0.4 + 1j * H * (k * uh + l * vh)) < 1e-12 * w                   # (13.79)
        tt = np.linspace(0.0, 2 * PI / w, 400, endpoint=False)
        ob = ch13.poincare_orbit(tt, k, 0.4, H, f)
        w1 = ob["omega"]
        assert math.isclose(ob["axis_ratio"], w1 / abs(f)) and math.isclose(np.ptp(ob["u"]) / np.ptp(ob["v"]), w1 / abs(f), rel_tol=1e-4)
        sense = np.sign(np.mean(ob["u"] * np.gradient(ob["v"], tt) - ob["v"] * np.gradient(ob["u"], tt)))
        assert sense == -np.sign(f) and ob["sense"] == ("clockwise" if f > 0 else "counter-clockwise")
        assert close(np.gradient(ob["x"], tt)[2:-2], ob["u"][2:-2], rtol=1e-3, atol=1e-6 * np.max(np.abs(ob["u"])))
        ti = np.linspace(0.0, 2 * PI / abs(f), 2001)
        ui, vi, xi, yi = ch13.inertial_oscillation(ti, 0.2, -0.1, f)                         # K -> 0: inertial circle
        assert close(np.hypot(ui, vi), math.hypot(0.2, 0.1)) and math.isclose(np.ptp(xi) / 2, ch13.inertial_radius(math.hypot(0.2, 0.1), f), rel_tol=1e-5)
        assert np.sign(np.mean(ui * np.gradient(vi, ti) - vi * np.gradient(ui, ti))) == -np.sign(f)
        assert close(np.gradient(xi, ti)[2:-2], ui[2:-2], atol=1e-5) and abs(xi[-1]) + abs(yi[-1]) < 1e-9   # closed circle
    assert math.isclose(ch13.poincare_omega(0.0, F_S, c), F_N)                               # omega -> |f|, not 0, as K -> 0
    assert math.isclose(ch13.inertial_radius(I["q_inertial"], F_S), I["q_inertial"] / F_N)
    with pytest.raises(ValueError):
        ch13.poincare_fields(0.0, 0.0, 0.0, 0.0, 0.0, 0.1, H, F_N)


def test_kelvin_V2_derivation():  # V2, D15 (13.84)–(13.87): which way a trapped wave travels
    x, y, t = sp.symbols("x y t", real=True)
    g, H, k, eta0 = sp.symbols("g H k eta_0", positive=True)
    f = sp.Symbol("f", real=True, nonzero=True)
    c = sp.sqrt(g * H)
    for d in (1, -1):
        eta = eta0 * sp.exp(-d * f * y / c) * sp.cos(k * (x - d * c * t))                    # (13.87) for travel toward d x
        u = d * sp.sqrt(g / H) * eta
        assert sp.simplify(eta.diff(t) + H * u.diff(x)) == 0                                 # (13.84a)
        assert sp.simplify(u.diff(t) + g * eta.diff(x)) == 0                                 # (13.84b)
        assert sp.simplify(f * u + g * eta.diff(y)) == 0                                     # (13.84c): exact geostrophy across
        # the envelope exp(-d f y/c) decays into y > 0 only when d f > 0: coast on the right for f > 0, on the left for f < 0
    F = sp.Function("F")(y)
    om = sp.Symbol("omega", positive=True)
    e, uu = F * sp.cos(k * x - om * t), sp.Function("G")(y) * sp.cos(k * x - om * t)
    sol = sp.solve([sp.simplify((e.diff(t) + H * uu.diff(x)) / sp.sin(k * x - om * t)), sp.simplify((uu.diff(t) + g * e.diff(x)) / sp.sin(k * x - om * t))], [sp.Function("G")(y), om], dict=True)
    assert any(sp.simplify(s[om] - k * c) == 0 for s in sol)                                 # (13.86): omega = k sqrt(gH)


def test_kelvin_V1_residuals_trapping_and_hemispheres():  # V1 + V7, C11 (13.84), (13.87)
    H, eta0, k = 61.0, 0.3, 2 * PI / 5.2e5
    c = math.sqrt(G0 * H)
    for f in HEMI:
        d = 1 if f > 0 else -1
        res = ch13.kelvin_residuals(3.3e4, 4.7e4, 111.0, eta0, k, H, f, h=10.0, ht=1.0)
        assert all(abs(v) < 1e-7 for v in res.values())
        side = ch13.kelvin_decay_side(f, d)
        assert side["trapped"] is True and side["coast_on"] == ("right" if f > 0 else "left")
        assert ch13.kelvin_decay_side(f, -d)["trapped"] is False
        with pytest.raises(ValueError):
            ch13.kelvin_wave(0.0, 1.0, 0.0, eta0, k, H, f, direction=-d)                     # would grow offshore
        y = np.linspace(0.0, 6 * c / abs(f), 61)
        e, u = ch13.kelvin_wave(0.0, y, 0.0, eta0, k, H, f)
        assert close(e, eta0 * np.exp(-y * abs(f) / c), rtol=1e-13) and close(u, d * math.sqrt(G0 / H) * e)
        assert math.isclose(e[10] / e[0], math.exp(-y[10] / ch13.rossby_radius(c, f)), rel_tol=1e-12)
        e1, _ = ch13.kelvin_wave(d * c * 500.0, y, 500.0, eta0, k, H, f)                     # the pattern moves at c toward d x
        assert close(e1, e, rtol=1e-12)
        assert np.all(np.isnan(ch13.kelvin_wave(0.0, np.array([-1.0]), 0.0, eta0, k, H, f)[0]))   # land
    assert math.isclose(ch13.kelvin_omega(k, c), k * c)
    assert math.isclose(ch13.rossby_radius(c, F_S), c / F_N)
    assert math.isclose(ch13.rossby_radius_internal(1.1e-2, 9e3, F_S, with_pi=False), 1.1e-2 * 9e3 / F_N)
    assert math.isclose(ch13.rossby_radius_internal(1.1e-2, 9e3, F_S, n=2), 1.1e-2 * 9e3 / (2 * PI * F_N))
    gp = G0 * I["drho"] / (I["rho_ocean"] + I["drho"])
    assert math.isclose(ch13.rossby_radius_two_layer(I["H1"], I["rho_ocean"], I["rho_ocean"] + I["drho"], F_S), math.sqrt(gp * I["H1"]) / F_N)
    assert math.isclose(gp, WAV.reduced_gravity_book(I["rho_ocean"], I["rho_ocean"] + I["drho"], G0), rel_tol=1e-12)   # ch07's g'
    with pytest.raises(ValueError):
        ch13.rossby_radius_two_layer(100.0, 1028.0, 1027.0, F_N)


def test_geostrophic_adjustment_V2_derivation():  # V2, D16 (ours): end state, PV, 1/3 of the released energy (claim 5)
    x = sp.Symbol("x", real=True)
    xp = sp.Symbol("x_p", positive=True)
    eta0, g, H = sp.symbols("eta_0 g H", positive=True)
    for s in (1, -1):                                                          # both hemispheres: f = s |f|
        fa = sp.Symbol("fabs", positive=True)
        f = s * fa
        Lam = sp.sqrt(g * H) / fa
        # solve the ODE of step 6 from scratch on x > 0 and impose boundedness and oddness
        E = sp.Function("E")
        gen = sp.dsolve(sp.Eq(E(xp).diff(xp, 2) - E(xp) / Lam ** 2, -eta0 / Lam ** 2), E(xp)).rhs
        consts = sorted(gen.free_symbols - {xp, eta0, g, H, fa}, key=str)
        cand = [gen.subs(dict(zip(consts, vals))) for vals in ((0, sp.Symbol("a")), (sp.Symbol("a"), 0))]
        bounded = [e for e in cand if sp.limit(e.subs({eta0: 1, g: 1, H: 1, fa: 1, sp.Symbol("a"): 1}), xp, sp.oo).is_finite][0]
        a = sp.solve(bounded.subs(xp, 0), sp.Symbol("a"))[0]
        eta_r = sp.simplify(bounded.subs(sp.Symbol("a"), a))
        assert sp.simplify(eta_r - eta0 * (1 - sp.exp(-xp / Lam))) == 0
        v_r = g / f * eta_r.diff(xp)                                           # steady: f v = g d(eta)/dx
        assert sp.simplify(v_r - s * g * eta0 / sp.sqrt(g * H) * sp.exp(-xp / Lam)) == 0     # the jet flips with sign(f)
        assert sp.simplify(v_r.diff(xp) - f * eta_r / H + f * eta0 / H) == 0   # linear PV: v_x - f eta/H keeps its initial value
        eta_l = -eta_r.subs(xp, -x)                                            # odd extension to x < 0
        assert sp.simplify(eta_l.subs(x, 0)) == 0 and sp.simplify(eta_l.diff(x).subs(x, 0) - eta_r.diff(xp).subs(xp, 0)) == 0
        PE = 2 * sp.integrate(sp.Rational(1, 2) * g * (eta0 ** 2 - eta_r ** 2), (xp, 0, sp.oo))
        KE = 2 * sp.integrate(sp.Rational(1, 2) * H * v_r ** 2, (xp, 0, sp.oo))
        assert sp.simplify(KE / PE - sp.Rational(1, 3)) == 0
        assert sp.simplify(PE - 3 * g * eta0 ** 2 * Lam / 2) == 0 and sp.simplify(KE - g * eta0 ** 2 * Lam / 2) == 0
    xs, ts = sp.symbols("xs ts", real=True)                                    # steps 1–3: the conserved quantity
    u, v, e = (sp.Function(n)(xs, ts) for n in ("u", "v", "e"))
    fs = sp.Symbol("f", real=True)
    pv_t = (v.diff(xs) - fs * e / H).diff(ts)
    assert sp.expand(pv_t - ((v.diff(ts) + fs * u).diff(xs) - fs / H * (e.diff(ts) + H * u.diff(xs)))) == 0


def test_geostrophic_adjustment_V1_end_state_energy_and_jet_sign():  # V1 + V4 + V7, C12 (N19): claim 5
    H, eta0 = 1.3, 0.05
    c = math.sqrt(G0 * H)
    for f in HEMI:
        Lam = c / abs(f)
        x = np.linspace(-9 * Lam, 9 * Lam, 3601)
        e, v = ch13.geostrophic_adjustment_1d(x, eta0, H, f)
        assert close(e, -e[::-1], atol=1e-18) and close(v, v[::-1])
        assert np.all(np.sign(v) == np.sign(f)) and math.isclose(np.abs(v).max(), G0 * eta0 / c, rel_tol=1e-12)
        dx = x[1] - x[0]
        inner = slice(5, -5)
        geo = f * v - G0 * np.gradient(e, dx)                                  # f v = g d(eta)/dx away from the kink at x = 0
        mask = np.abs(x) > 3 * dx
        assert np.max(np.abs(geo[mask])) < 1e-5 * G0 * eta0 / Lam
        pv = np.gradient(v, dx) - f * e / H                                    # the linear PV equals its initial value everywhere
        assert np.max(np.abs((pv + f * eta0 * np.sign(x) / H)[mask][inner])) < 1e-5 * abs(f) * eta0 / H
        en = ch13.adjustment_energy(eta0, H, f, rho=1027.0)
        pe = 2 * quad(lambda s: 0.5 * 1027.0 * G0 * (eta0 ** 2 - (eta0 * (1 - math.exp(-s / Lam))) ** 2), 0, 60 * Lam, epsrel=1e-12, limit=300)[0]
        ke = 2 * quad(lambda s: 0.5 * 1027.0 * H * (G0 * eta0 / c * math.exp(-s / Lam)) ** 2, 0, 60 * Lam, epsrel=1e-12)[0]
        assert math.isclose(en["pe_released"], pe, rel_tol=1e-10) and math.isclose(en["ke_jet"], ke, rel_tol=1e-10)
        assert abs(en["ratio"] - 1 / 3) < 1e-14 and math.isclose(en["radiated"], pe - ke, rel_tol=1e-9)
        fin = ch13.adjustment_energy(eta0, H, f, L=2.0 * Lam, rho=1027.0)
        pe2 = 2 * quad(lambda s: 0.5 * 1027.0 * G0 * (eta0 ** 2 - (eta0 * (1 - math.exp(-s / Lam))) ** 2), 0, 2 * Lam, epsrel=1e-12)[0]
        ke2 = 2 * quad(lambda s: 0.5 * 1027.0 * H * (G0 * eta0 / c * math.exp(-s / Lam)) ** 2, 0, 2 * Lam, epsrel=1e-12)[0]
        assert math.isclose(fin["pe_released"], pe2, rel_tol=1e-10) and math.isclose(fin["ke_jet"], ke2, rel_tol=1e-10)
    assert math.isclose(json.loads((REF / "explainer_constants.json").read_text(encoding="utf-8"))["adjustment_ratio"], 1 / 3, rel_tol=1e-14)


# =====================================================================================================================
# 6. §13.13 — potential vorticity, flow over a step (C13)
# =====================================================================================================================
def test_potential_vorticity_V2_derivation():  # V2, D17 (13.88)–(13.94), every intermediate line of the page
    x, y, t = sp.symbols("x y t", real=True)
    g, f0, beta = sp.symbols("g f_0 beta", positive=True)
    u, v, h = (sp.Function(n)(x, y, t) for n in ("u", "v", "h"))
    b = sp.Function("b")(x, y)                                                 # an uneven bottom: eta = h + b
    f = f0 + beta * y
    D = lambda q: q.diff(t) + u * q.diff(x) + v * q.diff(y)  # noqa: E731
    e88 = D(u) - f * v + g * (h + b).diff(x)
    e89 = D(v) + f * u + g * (h + b).diff(y)
    e90 = h.diff(t) + (u * h).diff(x) + (v * h).diff(y)
    zeta, div = v.diff(x) - u.diff(y), u.diff(x) + v.diff(y)
    e91 = zeta.diff(t) + (u * v.diff(x) + v * v.diff(y)).diff(x) - (u * u.diff(x) + v * u.diff(y)).diff(y) + f * div + beta * v
    assert sp.expand(e91 - (e89.diff(x) - e88.diff(y))) == 0                   # (13.91): cross-differentiate; beta v = d(f u)/dy part
    e92 = D(zeta) + (zeta + f) * div + beta * v
    assert sp.expand(e92 - e91) == 0                                           # (13.92): the four nonlinear terms regroup
    assert sp.expand(e90 - (D(h) + h * div)) == 0 and sp.simplify(D(f) - beta * v) == 0
    e93 = D(zeta + f) - (zeta + f) / h * D(h)
    assert sp.simplify(e93 - (e92 - (zeta + f) / h * e90)) == 0                # (13.93): eliminate the divergence with continuity
    e94 = D((zeta + f) / h)
    assert sp.simplify(e94 - e93 / h) == 0                                     # (13.94): quotient rule
    assert ch13.pv_conservation_sympy()["ok"]
    assert math.isclose(ch13.potential_vorticity(-2e-5, F_S, 310.0), (-2e-5 + F_S) / 310.0)
    for fs in HEMI:                                                            # stretching makes (zeta + f) larger in magnitude
        z1 = ch13.step_vorticity(fs, 4000.0, 4400.0)
        assert math.isclose(ch13.potential_vorticity(z1, fs, 4400.0), ch13.potential_vorticity(0.0, fs, 4000.0), rel_tol=1e-14)
        assert np.sign(z1) == np.sign(fs) and math.isclose(z1, fs * 0.1)


def test_flow_over_step_V2_derivation():  # V2, D18 (ours): streamline ODE both ways, and the far-downstream latitude (claim 6)
    x = sp.Symbol("x", real=True)
    beta, f0, h0, h1, Uabs = sp.symbols("beta f_0 h_0 h_1 U", positive=True)
    Yp = f0 * (h1 - h0) / (beta * h0)
    kap = sp.sqrt(beta / Uabs)
    pv = lambda Y, U, h: (U * Y.diff(x, 2) + f0 + beta * Y) / h - f0 / h0  # noqa: E731   linearised (zeta + f)/h - f0/h0
    Ye = Yp * (1 - sp.cos(kap * x))                                            # eastward flow, downstream (x > 0)
    assert sp.simplify(pv(Ye, Uabs, h1)) == 0 and Ye.subs(x, 0) == 0 and Ye.diff(x).subs(x, 0) == 0
    assert sp.simplify(Uabs * Ye.diff(x, 2).subs(x, 0) - f0 * (h1 - h0) / h0) == 0      # the jump of the page: zeta = f (h1 - h0)/h0
    assert sp.simplify(2 * sp.pi / kap - 2 * sp.pi * sp.sqrt(Uabs / beta)) == 0
    Yup, Ydn = Yp / 2 * sp.exp(-kap * x), Yp * (1 - sp.exp(kap * x) / 2)       # westward flow: arrives from x > 0
    assert sp.simplify(pv(Yup, -Uabs, h0)) == 0 and sp.simplify(pv(Ydn, -Uabs, h1)) == 0
    assert sp.simplify(Yup.subs(x, 0) - Ydn.subs(x, 0)) == 0 and sp.simplify((Yup.diff(x) - Ydn.diff(x)).subs(x, 0)) == 0
    assert sp.limit(Yup, x, sp.oo) == 0 and sp.simplify(sp.limit(Ydn, x, -sp.oo) - Yp) == 0
    jump = -Uabs * (Ydn.diff(x, 2) - Yup.diff(x, 2)).subs(x, 0)
    assert sp.simplify(jump - f0 * (h1 - h0) / h0) == 0                         # the same vorticity jump across the step
    # exact statement, no linearisation: a uniform stream far downstream has zeta = 0, so f/h1 = f0/h0 fixes its latitude
    Yinf = sp.Symbol("Y_inf", real=True)
    sol = sp.solve(sp.Eq((f0 + beta * Yinf) / h1, f0 / h0), Yinf)[0]
    assert sp.simplify(sol - Yp) == 0
    back_at_origin = sp.simplify(((0 + f0 + beta * 0) / h1 - f0 / h0))          # "back at its original latitude" with zeta = 0
    assert sp.simplify(back_at_origin - f0 * (h0 - h1) / (h0 * h1)) == 0 and back_at_origin != 0   # violates (13.94) unless h1 = h0


def test_flow_over_step_V1_coded_streamlines_and_permanent_shift():  # V1 + V4, N102, N103: claim 6, numerically
    beta, f0, h0, h1, U = float(ch13.beta_parameter(LAT)), F_N, 4000.0, 3800.0, I["U_mean"]
    Yp = f0 * (h1 - h0) / (beta * h0)
    x = np.linspace(-9e6, 9e6, 36001)
    dx = x[1] - x[0]
    for Us in (U, -U):
        r = ch13.flow_over_step(x, Us, beta, f0, h0, h1)
        Y = r["Y"]
        upstream = x < 0 if Us > 0 else x > 0
        h = np.where(upstream, h0, h1)
        Ypp = (Y[2:] - 2 * Y[1:-1] + Y[:-2]) / dx ** 2
        q = (Us * Ypp + f0 + beta * Y[1:-1]) / h[1:-1]                          # potential vorticity along the streamline
        keep = np.abs(x[1:-1]) > 2 * dx
        assert np.max(np.abs(q[keep] * h0 / f0 - 1)) < 2e-7                     # conserved on both sides of the step
        assert close(r["zeta"][1:-1][keep], Us * Ypp[keep], atol=2e-4 * abs(f0) * 0.05)
        assert math.isclose(r["decay_length"], math.sqrt(U / beta))
        if Us > 0:
            assert np.all(Y[x <= 0] == 0.0) and math.isclose(r["wavelength"], ch13.stationary_rossby_wavelength(U, beta))
            lam = r["wavelength"]
            down = (x > 0) & (x < lam)
            assert math.isclose(np.mean(Y[down]), Yp, rel_tol=2e-3)             # the meander is about the NEW latitude Y_p
            assert math.isclose(Y[down].min(), 2 * Yp, rel_tol=1e-6) and Yp < 0   # shallower: southward for f > 0
        else:
            assert r["wavelength"] is None and np.all(np.diff(Y) > 0) and Y[x > 0][0] < 0   # monotonic; turns before the step
            assert math.isclose(Y[np.argmin(np.abs(x))], Yp / 2, rel_tol=1e-9)
            assert math.isclose(Y[0], Yp, rel_tol=1e-4) and abs(Y[-1]) < 1e-4 * abs(Yp)       # settles at Y_p, not back at 0
            q_back = (0.0 + f0 + beta * 0.0) / h1                               # the "returns to its latitude" alternative
            assert abs(q_back * h0 / f0 - 1) > 0.05                             # misses PV conservation by (h0 - h1)/h1
    south = ch13.flow_over_step(np.array([1e7]), -U, beta, F_S, h0, h1)["Y"]
    assert south > 0                                                            # f < 0: |f| must fall, which is northward in the south
    with pytest.raises(ValueError):
        ch13.flow_over_step(x, 0.0, beta, f0, h0, h1)


# =====================================================================================================================
# 7. §13.14 — internal waves with rotation (C14)
# =====================================================================================================================
def test_w_equation_V2_derivation():  # V2, D19 (13.95) -> (13.96), N may depend on z and is never differentiated
    x, y, z, t = sp.symbols("x y z t", real=True)
    f, g, rho0 = sp.symbols("f g rho_0", positive=True)
    N2 = sp.Function("N2")(z)
    u, v, w, p, rho = (sp.Function(n)(x, y, z, t) for n in ("u", "v", "w", "p", "rho"))
    R1 = u.diff(x) + v.diff(y) + w.diff(z)                                     # the five members of (13.95), all on the left
    R2 = u.diff(t) - f * v + p.diff(x) / rho0
    R3 = v.diff(t) + f * u + p.diff(y) / rho0
    R4 = w.diff(t) + p.diff(z) / rho0 + rho * g / rho0
    R5 = rho.diff(t) - rho0 * N2 / g * w
    lapH = lambda q: q.diff(x, 2) + q.diff(y, 2)  # noqa: E731
    divR, curlR = R2.diff(x) + R3.diff(y), R3.diff(x) - R2.diff(y)
    A = divR.diff(t) + f * curlR                                               # steps 1–6: divergence and curl remove zeta
    assert sp.expand(A - ((u.diff(x) + v.diff(y)).diff(t, 2) + f ** 2 * (u.diff(x) + v.diff(y)) + lapH(p).diff(t) / rho0)) == 0
    B = R4.diff(t) - g / rho0 * R5                                             # steps 7–8: vertical + density remove rho
    assert sp.expand(B - (w.diff(t, 2) + p.diff(z, t) / rho0 + N2 * w)) == 0
    weq = (lapH(w) + w.diff(z, 2)).diff(t, 2) + N2 * lapH(w) + f ** 2 * w.diff(z, 2)   # (13.96)
    combo = lapH(B) - A.diff(z) + (R1.diff(t, 2) + f ** 2 * R1).diff(z)        # steps 9–11: eliminate p, use continuity
    assert sp.expand(combo - weq) == 0                                         # every solution of (13.95) satisfies (13.96)
    assert ch13.rotating_internal_wave_set_sympy()["ok"] and ch13.w_equation_rotating_sympy()["ok"]
    lim = ch13.w_equation_rotating_sympy()["f0_limit"]
    sy = {str(q): q for q in lim.free_symbols}
    assert sp.simplify(lim - sy["N"] ** 2 * (sy["k"] ** 2 + sy["l"] ** 2) / (sy["k"] ** 2 + sy["l"] ** 2 + sy["m"] ** 2)) == 0


def test_vertical_structure_and_dispersion_V2_derivation():  # V2, D20 (13.97)–(13.99), (13.112) as printed
    x, y, z, t = sp.symbols("x y z t", real=True)
    k, l, m, om, N, f, th = sp.symbols("k l m omega N f theta", positive=True)
    Wh = sp.Function("What")(z)
    w = Wh * sp.exp(sp.I * (k * x + l * y - om * t))
    weq = (w.diff(x, 2) + w.diff(y, 2) + w.diff(z, 2)).diff(t, 2) + N ** 2 * (w.diff(x, 2) + w.diff(y, 2)) + f ** 2 * w.diff(z, 2)
    m2 = (k ** 2 + l ** 2) * (N ** 2 - om ** 2) / (om ** 2 - f ** 2)           # (13.99)
    assert sp.simplify(weq / sp.exp(sp.I * (k * x + l * y - om * t)) - (f ** 2 - om ** 2) * (Wh.diff(z, 2) + m2 * Wh)) == 0
    om2 = (N ** 2 * k ** 2 + f ** 2 * m ** 2) / (k ** 2 + m ** 2)
    assert sp.simplify((om2 - f ** 2) - k ** 2 / m ** 2 * (N ** 2 - om2)) == 0  # (13.112): omega^2 - f^2 = (k^2/m^2)(N^2 - omega^2)
    assert sp.simplify(om2.subs(m, k * sp.tan(th)) - (f ** 2 * sp.sin(th) ** 2 + N ** 2 * sp.cos(th) ** 2)) == 0   # tan(theta) = m/k
    assert sp.simplify(om2.subs(f, 0) - N ** 2 * k ** 2 / (k ** 2 + m ** 2)) == 0            # f -> 0: ch07's omega = N cos(theta)


def test_group_velocity_inertia_gravity_V2_derivation():  # V2, D21 (ours; the book leaves it as an exercise)
    k, m, N, f = sp.symbols("k m N f", positive=True)
    om = sp.sqrt((N ** 2 * k ** 2 + f ** 2 * m ** 2) / (k ** 2 + m ** 2))
    fac = (N ** 2 - f ** 2) * k * m / ((m ** 2 + k ** 2) ** sp.Rational(3, 2) * sp.sqrt(m ** 2 * f ** 2 + k ** 2 * N ** 2))
    assert sp.simplify(om.diff(k) - fac * m) == 0 and sp.simplify(om.diff(m) + fac * k) == 0
    assert sp.simplify(om / (k ** 2 + m ** 2) * (k * om.diff(k) + m * om.diff(m))) == 0      # phase velocity is normal to c_g
    assert sp.simplify(om.diff(k).subs(f, 0) - N * m ** 2 / (k ** 2 + m ** 2) ** sp.Rational(3, 2)) == 0   # ch07's result


def test_inertia_gravity_V1_dispersion_band_regimes_group_velocity():  # V1 + V7, C14 (13.99), (13.112); claim 10
    N, f = 5.3e-3, F_S
    k = np.geomspace(1e-5, 1e-2, 40)
    m = np.geomspace(3e-4, 3e-2, 40)[::-1]
    w = ch13.inertia_gravity_omega(k, m, N, f)
    assert np.all((w > abs(f)) & (w < N))
    assert close(w ** 2 - f ** 2, k ** 2 / m ** 2 * (N ** 2 - w ** 2), rtol=1e-10)     # (13.112) exactly as printed
    assert close(ch13.inertia_gravity_m2(k, 0.0, w, N, f), m ** 2, rtol=1e-9)          # (13.99) inverts it
    th = np.arctan2(m, k)
    assert close(w ** 2, f ** 2 * np.sin(th) ** 2 + N ** 2 * np.cos(th) ** 2, rtol=1e-12)
    assert close(ch13.inertia_gravity_omega(k, m, N, -f), w)                           # f enters squared
    # f -> 0 reduces to chapter 7's internal waves, frequency and group velocity alike
    assert close(ch13.inertia_gravity_omega(k, m, N, 0.0), WAV.internal_wave_omega(k, m, N), rtol=1e-14)
    assert close(ch13.inertia_gravity_omega(0.7 * k, m, N, 0.0, l=0.4 * k), WAV.internal_wave_omega(0.7 * k, m, N, 0.4 * k), rtol=1e-14)
    cg0 = ch13.inertia_gravity_group_velocity(1.1e-3, 2.3e-3, N, 0.0)
    ch7 = WAV.group_velocity_vector(lambda K: WAV.internal_wave_omega(K[0], K[1], N), np.array([1.1e-3, 2.3e-3]))
    assert close(cg0, ch7, rtol=1e-6)
    for kk, mm in ((1.1e-3, 2.3e-3), (4e-5, -6e-3), (-3e-4, 9e-4)):
        cg = ch13.inertia_gravity_group_velocity(kk, mm, N, f)
        h = 1e-9
        fd = ((ch13.inertia_gravity_omega(kk + h, mm, N, f) - ch13.inertia_gravity_omega(kk - h, mm, N, f)) / (2 * h),
              (ch13.inertia_gravity_omega(kk, mm + h, N, f) - ch13.inertia_gravity_omega(kk, mm - h, N, f)) / (2 * h))
        assert close(cg, fd, rtol=1e-5) and abs(cg[0] * kk + cg[1] * mm) < 1e-12 * math.hypot(*cg) * math.hypot(kk, mm)
        assert np.sign(cg[1]) == -np.sign(mm) * np.sign(kk) * np.sign(kk)                    # vertical phase and group: opposite
    # approximations and the regime with its stated errors in m^2
    kk, mm = 2e-5, 4e-3
    full = ch13.inertia_gravity_omega(kk, mm, N, f)
    assert math.isclose(ch13.inertia_gravity_omega(kk, mm, N, f, approx="hydrostatic"), math.sqrt(f * f + (N * kk / mm) ** 2))
    assert math.isclose(ch13.inertia_gravity_omega(kk, mm, N, f, approx="midrange"), N * kk / mm)
    assert math.isclose(ch13.inertia_gravity_omega(kk, mm, N, f, approx="nonrotating"), N * kk / math.hypot(kk, mm))
    reg = ch13.inertia_gravity_regime(full, N, f)
    m2x = ch13.inertia_gravity_m2(kk, 0.0, full, N, f)
    assert math.isclose(reg["err_hydrostatic"], (kk ** 2 * N ** 2 / (full ** 2 - f ** 2)) / m2x - 1, rel_tol=1e-9)
    assert math.isclose(reg["err_nonrotating"], (kk ** 2 * (N ** 2 - full ** 2) / full ** 2) / m2x - 1, rel_tol=1e-9)
    assert math.isclose(reg["err_midrange"], (kk ** 2 * N ** 2 / full ** 2) / m2x - 1, rel_tol=1e-9)
    assert ch13.inertia_gravity_regime(1.02 * abs(f), N, f)["regime"] == "near-inertial"
    assert ch13.inertia_gravity_regime(0.98 * N, N, f)["regime"] == "near-buoyancy"
    assert ch13.inertia_gravity_regime(math.sqrt(abs(f) * N), N, f)["regime"] == "mid"
    assert [ch13.inertia_gravity_band(v, N, f)["where"] for v in (0.5 * abs(f), 3 * abs(f), 2 * N)] == ["below f", "in band", "above N"]
    assert ch13.inertia_gravity_m2(kk, 0.0, 0.5 * abs(f), N, f) < 0 and ch13.inertia_gravity_m2(kk, 0.0, 2 * N, N, f) < 0
    with pytest.raises(ValueError):
        ch13.inertia_gravity_m2(kk, 0.0, abs(f), N, f)
    with pytest.raises(ValueError):
        ch13.inertia_gravity_regime(2 * N, N, f)


def test_inertia_gravity_V1_fields_hodograph_and_w_equation_residual():  # V1 + V7, C14 (13.96), (13.108), (13.110)
    N = 5.3e-3
    for f in HEMI:
        k, om = 2 * PI / 8.0e3, 2.9 * abs(f)
        m = math.sqrt(ch13.inertia_gravity_m2(k, 0.0, om, N, f))
        wave = lambda x, y, z, t: math.cos(k * x + m * z - om * t)  # noqa: E731
        r_on = ch13.w_equation_rotating_residual(wave, 10.0, 0.0, -300.0, 50.0, N, f, h=0.5, ht=2.0)
        scale = om ** 2 * (k ** 2 + m ** 2)
        assert abs(r_on) < 1e-5 * scale
        off = lambda x, y, z, t: math.cos(k * x + m * z - 1.3 * om * t)  # noqa: E731
        assert abs(ch13.w_equation_rotating_residual(off, 10.0, 0.0, -300.0, 50.0, N, f, h=0.5, ht=2.0)) > 0.3 * scale
        # the WKB fields with constant N are an exact plane wave: continuity and the y-momentum equation by differences
        x = np.linspace(0.0, 1.6e4, 9)
        z = np.linspace(-800.0, 0.0, 1601)
        Nf = lambda s: N + 0.0 * s  # noqa: E731
        F = lambda t, dxs=0.0: ch13.inertia_gravity_fields(x + dxs, z, t, k, om, Nf, f)  # noqa: E731
        u0, v0, w0 = F(0.0)
        hx, ht = 0.5, 1.0
        ux = (F(0.0, hx)[0] - F(0.0, -hx)[0]) / (2 * hx)
        wz = np.gradient(w0, z, axis=0)
        assert np.max(np.abs(ux + wz)[5:-5]) < 1e-4 * k * np.max(np.abs(u0))                 # du/dx + dw/dz = 0
        vt = (F(ht)[1] - F(-ht)[1]) / (2 * ht)
        assert np.max(np.abs(vt + f * u0)) < 1e-6 * abs(f) * np.max(np.abs(u0))              # dv/dt + f u = 0
        assert math.isclose(np.max(np.abs(v0)) / np.max(np.abs(u0)), abs(f) / om, rel_tol=1e-3)   # axis ratio |f|/omega
        tt = np.linspace(0.0, 2 * PI / om, 400, endpoint=False)
        hu, hv = ch13.inertia_gravity_hodograph(tt, om, f)
        assert math.isclose(np.ptp(hv) / np.ptp(hu), abs(f) / om, rel_tol=1e-4)
        assert np.sign(np.mean(hu * np.gradient(hv, tt) - hv * np.gradient(hu, tt))) == -np.sign(f)   # clockwise for f > 0
    zz = np.linspace(-500.0, 0.0, 501)
    wk = ch13.wkb_vertical_structure(zz, 0.02, A0=1.5)
    assert close(wk, 1.5 / math.sqrt(0.02) * np.exp(1j * 0.02 * (zz - zz[0])), rtol=1e-12)   # constant m: exact
    assert close(ch13.wkb_vertical_structure(zz, 0.02, sign=-1), np.conj(ch13.wkb_vertical_structure(zz, 0.02)))
    num = ch13.vertical_structure_solve(zz, lambda s: N, 2 * PI / 8.0e3, 3e-4, 1e-4, w0=0.0, dw0=1.0)
    m_c = math.sqrt(ch13.inertia_gravity_m2(2 * PI / 8.0e3, 0.0, 3e-4, N, 1e-4))
    assert np.max(np.abs(num - np.sin(m_c * (zz - zz[0])) / m_c)) < 1e-7 / m_c                 # (13.100) with constant m


def test_wkb_V3_error_falls_as_one_over_Hm():  # V3-style scaling, N114 (13.104): claim 7 (wkb_error)
    f, om, Hs, N0 = 1e-4, 3e-4, 1000.0, 5e-3
    Nf = lambda s: N0 * (1 + 0.5 * np.tanh(s / Hs))  # noqa: E731
    z = np.linspace(-3 * Hs, 3 * Hs, 3001)
    Hm, err = [], []
    for eps in (0.4, 0.2, 0.1, 0.05):
        k = (1 / (eps * Hs)) * math.sqrt((om * om - f * f) / (N0 * N0 - om * om))
        r = ch13.wkb_error(z, Nf, k, om, f)
        Hm.append(r["Hm"])
        err.append(r["max_rel_error"])
    assert abs(observed_order(1 / np.array(Hm[1:]), err[1:]) - 1.0) < 0.1                     # error proportional to 1/(H m)
    prod = np.array(err) * np.array(Hm)
    assert abs(prod[-1] / prod[-2] - 1) < 0.01 and 0.1 < prod[-1] < 0.5 and err[-1] < 0.01
    with pytest.raises(ValueError):
        ch13.wkb_error(z, Nf, 1e-3, 2 * N0, f)                                               # evanescent somewhere


def test_lee_waves_V1_vertical_wavenumber_field_and_tilt():  # V1 + V7, N126 (13.113)
    U, N, h0 = 11.0, 1.1e-2, 120.0
    k = 2 * PI / 2.3e4
    m = ch13.lee_wave_m(U, N, k)
    assert math.isclose(m, math.sqrt(N * N / U ** 2 - k * k))
    assert math.isclose(ch13.inertia_gravity_omega(k, m, N, 0.0), k * U, rel_tol=1e-12)       # intrinsic frequency = kU: stationary
    x, z = np.linspace(0.0, 4.6e4, 401), np.linspace(0.0, 6e3, 241)
    lw = ch13.lee_wave_field(x, z, U, N, k, h0)
    assert lw["tilt"] == "upstream" and math.isclose(lw["m"], m)
    disp = z[:, None] - lw["psi"] / U                                                         # streamline displacement
    assert close(disp[0], h0 * np.cos(k * x), atol=1e-9)                                # follows the ground
    one = slice(0, 200)                                                                       # exactly one wavelength in x
    for j in (0, 20, 40, 77):
        ph = np.angle(np.sum(disp[j, one] * np.exp(-1j * k * x[one])))                        # phase of the crest at this height
        assert abs(np.exp(1j * ph) - np.exp(1j * m * z[j])) < 1e-9                            # crest at x = -m z/k: leans upstream
    assert close(lw["w"][:, 1:-1], U * np.gradient(disp, x, axis=1)[:, 1:-1], atol=1e-3 * U * h0 * k)   # w = U d(delta)/dx
    ev = ch13.lee_wave_field(x, z, U, N, 3 * N / U, h0)
    assert ev["m"] is None and ev["tilt"].startswith("none")
    assert np.max(np.abs(ev["psi"][-1] / U - z[-1])) < 1e-3 * h0                              # decays with height
    with pytest.raises(ValueError, match="evanescent"):
        ch13.lee_wave_m(U, N, 3 * N / U)


# =====================================================================================================================
# 8. §13.15–§13.16 — Rossby waves and barotropic instability (C15)
# =====================================================================================================================
def test_qg_vorticity_V2_derivation():  # V2, D22 (13.94) -> (13.114) -> (13.115) -> (13.117)
    x, y, t, eps = sp.symbols("x y t epsilon", real=True)
    g, H, f0, beta = sp.symbols("g H f_0 beta", positive=True)
    u, v, eta = (sp.Function(n)(x, y, t) for n in ("u", "v", "eta"))
    zeta = v.diff(x) - u.diff(y)
    D = lambda q: q.diff(t) + u * q.diff(x) + v * q.diff(y)  # noqa: E731
    f = f0 + beta * y
    e114 = (H + eta) * (D(zeta) + beta * v) - (zeta + f) * D(eta)              # (13.114) before f -> f0
    assert sp.simplify(D((zeta + f) / (H + eta)) * (H + eta) ** 2 - e114) == 0
    sub = {u: eps * u, v: eps * v, eta: eps * eta}
    lin = sp.expand(e114.subs(sub, simultaneous=True).doit()).coeff(eps, 1)
    assert sp.expand(lin - (H * zeta.diff(t) + H * beta * v - f * eta.diff(t))) == 0          # (13.115) with f -> f0 (y-scale small)
    e115 = H * zeta.diff(t) + H * beta * v - f0 * eta.diff(t)
    geo = {u: -g / f0 * eta.diff(y), v: g / f0 * eta.diff(x)}                  # (13.116)
    q = e115.subs(geo, simultaneous=True).doit()
    assert sp.simplify(zeta.subs(geo, simultaneous=True).doit() - g / f0 * (eta.diff(x, 2) + eta.diff(y, 2))) == 0
    e117 = (eta.diff(x, 2) + eta.diff(y, 2) - f0 ** 2 / (g * H) * eta).diff(t) + beta * eta.diff(x)   # (13.117), c^2 = gH
    assert sp.simplify(q * f0 / (g * H) - e117) == 0
    div_geo = (u.diff(x) + v.diff(y)).subs(geo, simultaneous=True).doit()
    assert sp.simplify(div_geo) == 0          # why the divergence must NOT be taken from (13.116): it would lose the stretching term
    assert ch13.qg_vorticity_sympy()["ok"]


def test_rossby_dispersion_V2_derivation():  # V2, D23 (13.118)–(13.120), circles, group velocity, the maximum frequency
    x, y, t = sp.symbols("x y t", real=True)
    k, l = sp.symbols("k l", real=True)
    beta, f0, c, U, om0 = sp.symbols("beta f_0 c U omega_0", positive=True)
    F = f0 ** 2 / c ** 2
    om = -beta * k / (k ** 2 + l ** 2 + F)                                     # (13.118)
    eta = sp.exp(sp.I * (k * x + l * y - om * t))
    assert sp.simplify(((eta.diff(x, 2) + eta.diff(y, 2) - F * eta).diff(t) + beta * eta.diff(x)) / eta) == 0
    assert sp.simplify(om / k + beta / (k ** 2 + l ** 2 + F)) == 0             # (13.119): westward for every wavelength
    circle = (k + beta / (2 * om0)) ** 2 + l ** 2 - ((beta / (2 * om0)) ** 2 - F)
    assert sp.simplify(circle - (k ** 2 + l ** 2 + F + beta * k / om0)) == 0   # completing the square of omega = omega_0
    D = k ** 2 + l ** 2 + F
    assert sp.simplify(om.diff(k) - beta * (k ** 2 - l ** 2 - F) / D ** 2) == 0 and sp.simplify(om.diff(l) - 2 * beta * k * l / D ** 2) == 0
    kmax = sp.solve(om.subs(l, 0).diff(k), k)
    assert set(sp.simplify(s) for s in kmax) == {f0 / c, -f0 / c}
    assert sp.simplify(om.subs({l: 0, k: -f0 / c}) - beta * c / (2 * f0)) == 0  # omega_max = beta c/(2 f0) at k c/f0 = -1
    assert sp.simplify(sp.limit(om / k, k, 0).subs(l, 0) + beta * c ** 2 / f0 ** 2) == 0     # long waves: -beta Lambda^2
    cx = U - beta / (k ** 2 + l ** 2)                                          # (13.120), barotropic: stationary where it vanishes
    ks = sp.solve(cx.subs(l, 0), k)
    assert any(sp.simplify(2 * sp.pi / s - 2 * sp.pi * sp.sqrt(U / beta)) == 0 for s in ks)


def test_rossby_V1_frequency_group_velocity_circle_and_limits():  # V1 + V7, C15 (13.118)–(13.120)
    lat = I["lat_rossby"]
    c = float(ch13.baroclinic_mode_speed(I["ocean_N"], I["ocean_H"], 1))
    beta = float(ch13.beta_parameter(lat))
    for f0 in (float(ch13.coriolis_parameter(lat)), -float(ch13.coriolis_parameter(lat))):
        Lam = c / abs(f0)
        k = -np.geomspace(0.02, 40.0, 300) / Lam
        w = ch13.rossby_omega(k, 0.0, beta, f0, c)
        assert np.all(w > 0) and close(w, -beta * k / (k ** 2 + 1 / Lam ** 2))          # k < 0 gives omega > 0 (trap T8)
        mx = ch13.rossby_max_frequency(beta, f0, c)
        assert math.isclose(mx["omega_max"], beta * c / (2 * abs(f0))) and math.isclose(mx["k"], -1 / Lam)
        assert w.max() <= mx["omega_max"] * (1 + 1e-12) and math.isclose(w.max(), mx["omega_max"], rel_tol=1e-3)
        assert math.isclose(ch13.rossby_omega(mx["k"], 0.0, beta, f0, c), mx["omega_max"], rel_tol=1e-13)
        cgx, cgy = ch13.rossby_group_velocity(k, 0.0, beta, f0, c)
        assert np.all(cgx[np.abs(k) * Lam < 0.99] < 0) and np.all(cgx[np.abs(k) * Lam > 1.01] > 0) and np.all(cgy == 0)
        assert math.isclose(ch13.rossby_phase_speed(-1e-3 / Lam, 0.0, beta, f0, c), ch13.rossby_long_wave_speed(beta, f0, c), rel_tol=1e-5)
        assert math.isclose(ch13.rossby_group_velocity(-1e-3 / Lam, 0.0, beta, f0, c)[0], -beta * Lam ** 2, rel_tol=1e-5)   # non-dispersive
        assert math.isclose(ch13.basin_crossing_time(5.1e6, lat, c), 5.1e6 / (beta * Lam ** 2), rel_tol=1e-12)
        for kk, ll, Um in ((-2.1 / Lam, 0.9 / Lam, 0.0), (-0.4 / Lam, -1.7 / Lam, 0.03), (3.0 / Lam, 0.2 / Lam, -0.02)):
            h = 1e-6 / Lam
            fd = ((ch13.rossby_omega(kk + h, ll, beta, f0, c, Um) - ch13.rossby_omega(kk - h, ll, beta, f0, c, Um)) / (2 * h),
                  (ch13.rossby_omega(kk, ll + h, beta, f0, c, Um) - ch13.rossby_omega(kk, ll - h, beta, f0, c, Um)) / (2 * h))
            assert close(ch13.rossby_group_velocity(kk, ll, beta, f0, c, Um), fd, rtol=1e-6)
            assert math.isclose(ch13.rossby_phase_speed(kk, ll, beta, f0, c, Um), ch13.rossby_omega(kk, ll, beta, f0, c, Um) / kk, rel_tol=1e-12)
            assert math.isclose(ch13.rossby_omega(kk, ll, beta, f0, c, Um) - ch13.rossby_omega(kk, ll, beta, f0, c), Um * kk, rel_tol=1e-9)   # Doppler
        w0 = 0.6 * mx["omega_max"]
        cir = ch13.rossby_omega_circle(w0, beta, f0, c)
        ang = np.linspace(0.0, 2 * PI, 37)
        kc, lc = cir["center_k"] + cir["radius"] * np.cos(ang), cir["radius"] * np.sin(ang)
        assert cir["exists"] and close(ch13.rossby_omega(kc, lc, beta, f0, c), w0, rtol=1e-10)
        cg = np.array(ch13.rossby_group_velocity(kc, lc, beta, f0, c))
        inward = (cir["center_k"] - kc) * cg[0] + (0.0 - lc) * cg[1]
        assert np.all(inward > 0)                                                             # group velocity points to the centre
        assert ch13.rossby_omega_circle(1.2 * mx["omega_max"], beta, f0, c)["exists"] is False
    U = I["U_mean"]
    lam_s = ch13.stationary_rossby_wavelength(U, beta)
    assert math.isclose(lam_s, 2 * PI * math.sqrt(U / beta)) and abs(ch13.rossby_phase_speed(2 * PI / lam_s, 0.0, beta, U=U)) < 1e-12 * U
    assert ch13.rossby_group_velocity(2 * PI / lam_s, 0.0, beta, U=U)[0] == pytest.approx(2 * U)   # the wake lies downstream
    kp, ap = ch13.rossby_packet_spectrum(-2e-5, 2e-6, 41)
    assert math.isclose(ap.sum(), 1.0) and math.isclose(kp[20], -2e-5) and math.isclose(kp[-1] - kp[0], 12e-6) and np.argmax(ap) == 20
    with pytest.raises(ValueError):
        ch13.stationary_rossby_wavelength(-3.0, beta)                                         # no stationary wave in easterlies


def test_qg_linear_evolve_V1_plane_wave_and_packet_group_velocity():  # V1, C15 (13.117) and curation note 8
    f0, beta = float(ch13.coriolis_parameter(I["lat_rossby"])), float(ch13.beta_parameter(I["lat_rossby"]))
    c, L, n = 2.9, 2.0e6, 64
    x = np.arange(n) * L / n
    X, Y = np.meshgrid(x, x)
    k, l = -3 * 2 * PI / L, 2 * 2 * PI / L
    om = ch13.rossby_omega(k, l, beta, f0, c)
    t = 0.37 * 2 * PI / abs(om)
    e = ch13.qg_linear_evolve(np.cos(k * X + l * Y) + 0.3 * np.sin(2 * k * X), x, x, t, beta, f0, c)
    om2 = ch13.rossby_omega(2 * k, 0.0, beta, f0, c)
    assert np.max(np.abs(e - np.cos(k * X + l * Y - om * t) - 0.3 * np.sin(2 * k * X - om2 * t))) < 1e-12
    assert np.max(np.abs(ch13.qg_linear_evolve(np.cos(l * Y), x, x, t, beta, f0, c) - np.cos(l * Y))) < 1e-13   # zonal flow is steady
    kk, aa = ch13.rossby_packet_spectrum(-2e-5, 2e-6, 41)
    xx = np.linspace(-4e6, 4e6, 16001)
    env = lambda tt: np.abs(np.exp(1j * np.multiply.outer(xx, kk)) @ ch13.qg_linear_evolve_1d(aa, kk, 0.0, tt, beta, f0, c))  # noqa: E731
    t1 = 3e6
    speed = (xx[np.argmax(env(t1))] - xx[np.argmax(env(0.0))]) / t1
    assert math.isclose(speed, ch13.rossby_group_velocity(-2e-5, 0.0, beta, f0, c)[0], rel_tol=0.03)   # packet moves at c_g
    assert np.sign(speed) != np.sign(ch13.rossby_phase_speed(-2e-5, 0.0, beta, f0, c))                 # short waves: energy goes east
    field = ch13.qg_linear_evolve_1d(aa, kk, 0.0, 0.0, beta, f0, c, x=np.array([0.0, 1e5]))
    assert math.isclose(field[0], 1.0, rel_tol=1e-12) and close(ch13.qg_linear_evolve_1d(aa, kk, 0.0, 0.0, beta, f0, c), aa)


def test_rayleigh_kuo_V2_derivation():  # V2, D24 (13.122) -> (13.123) -> the normal-mode equation and (13.124)
    x, y, t, eps = sp.symbols("x y t epsilon", real=True)
    beta, k = sp.symbols("beta k", positive=True)
    c = sp.Symbol("c")
    Psi, phi = sp.Function("Psi")(y), sp.Function("phi")(x, y, t)
    psi = Psi + eps * phi                                                      # u = -dpsi/dy, v = dpsi/dx (trap T14)
    u, v = -psi.diff(y), psi.diff(x)
    zeta = psi.diff(x, 2) + psi.diff(y, 2)
    full = zeta.diff(t) + u * zeta.diff(x) + v * (zeta.diff(y) + beta)         # (13.122) with f = f0 + beta y
    U = -Psi.diff(y)
    lap = phi.diff(x, 2) + phi.diff(y, 2)
    e123 = lap.diff(t) + U * lap.diff(x) + (beta - U.diff(y, 2)) * phi.diff(x)  # (13.123)
    assert sp.expand(sp.expand(full).coeff(eps, 1) - e123) == 0 and sp.expand(full).coeff(eps, 0) == 0
    ph = sp.Function("phihat")(y)
    mode = ph * sp.exp(sp.I * k * (x - c * t))
    nm = sp.simplify(e123.subs(phi, mode).doit() / (sp.I * k * sp.exp(sp.I * k * (x - c * t))))
    assert sp.simplify(nm - ((U - c) * (ph.diff(y, 2) - k ** 2 * ph) + (beta - U.diff(y, 2)) * ph)) == 0
    assert sp.simplify((-U.diff(y)).diff(y) + beta - (beta - U.diff(y, 2))) == 0   # d(zeta_bar + f)/dy = beta - U'' (13.124)


SECH2 = dict(U=lambda y: 1 / np.cosh(y) ** 2, Up=lambda y: -2 * np.tanh(y) / np.cosh(y) ** 2,
             Upp=lambda y: (4 * np.tanh(y) ** 2 - 2 / np.cosh(y) ** 2) / np.cosh(y) ** 2)


def _rk_shoot(k, beta, c0, Y=14.0, want_mode=False):
    """Independent route (ours, written for this test): integrate (U - c)(phi'' - k^2 phi) + (beta - U'') phi = 0 from
    y = -Y, where phi ~ exp(kappa y) with kappa^2 = k^2 + beta/c (Re kappa > 0), to the axis and require phi'(0) = 0
    (the sinuous mode of a symmetric jet); secant iteration on complex c."""
    def march(c, dense=False):
        kap = np.sqrt(k * k + beta / c + 0j)
        kap = -kap if kap.real < 0 else kap
        rhs = lambda y, s: [s[1], (k * k - (beta - SECH2["Upp"](y)) / (SECH2["U"](y) - c)) * s[0]]  # noqa: E731
        return solve_ivp(rhs, (-Y, 0.0), [1 + 0j, kap], method="DOP853", rtol=1e-11, atol=1e-13, dense_output=dense), kap
    F = lambda c: (lambda s: s.y[1, -1] / s.y[0, -1])(march(c)[0])  # noqa: E731
    c1, c2 = c0, c0 * (1 + 1e-3)
    f1, f2 = F(c1), F(c2)
    for _ in range(40):
        c3 = c2 - f2 * (c2 - c1) / (f2 - f1)
        c1, f1, c2 = c2, f2, c3
        f2 = F(c2)
        if abs(c2 - c1) < 1e-12:
            break
    if want_mode:
        return c2, march(c2, dense=True)
    return c2, abs(f2)


def test_rayleigh_kuo_V3_growth_rates_against_independent_shooting():  # V3 + V1, C15 (13.123)–(13.124): claim 9
    y = np.linspace(-6.0, 6.0, 4001)
    assert math.isclose(SECH2["Upp"](y).max(), 2 / 3, rel_tol=1e-6)            # largest U'' of the sech^2 jet
    found = {}
    for beta in (0.0, 0.3, 0.6):
        r = ch13.rayleigh_kuo_eigs(1.4, SECH2["U"], SECH2["Up"], SECH2["Upp"], beta, bc="decay", y_max=16.0, parity="even")
        cs, res = _rk_shoot(1.4, beta, r["c"])
        assert res < 1e-10 and abs(r["c"] - cs) < 1e-7 and r["stable"] is False
        assert math.isclose(r["growth_rate"], 1.4 * r["c"].imag)
        found[beta] = r["growth_rate"]
        crit = ch13.rayleigh_kuo_criterion(y, SECH2["U"](y), beta)
        assert crit["changes_sign"] is True and len(crit["y_zero"]) == (2 if beta == 0.0 else 4)
    assert [round(found[b], 4) for b in (0.0, 0.3, 0.6)] == [0.1191, 0.0744, 0.0151]         # ours, computed twice
    assert found[0.0] > found[0.3] > found[0.6] > 0                                          # beta stabilises
    gone = ch13.rayleigh_kuo_eigs(1.4, SECH2["U"], SECH2["Up"], SECH2["Upp"], 0.7, bc="decay", y_max=16.0, parity="even")
    assert gone["stable"] is True and gone["growth_rate"] == 0.0                             # beta > max U'': (13.124) forbids growth
    crit = ch13.rayleigh_kuo_criterion(y, SECH2["U"](y), 0.7)
    assert crit["changes_sign"] is False and crit["min"] > 0 and crit["y_zero"] == []
    G = ch13.absolute_vorticity_gradient(y, SECH2["U"](y), 0.3)
    assert np.max(np.abs(G - (0.3 - SECH2["Upp"](y)))) < 2e-5
    assert np.max(np.abs(ch13.absolute_vorticity_gradient(y, SECH2["U"], 0.3) - (0.3 - SECH2["Upp"](y)))) < 2e-5   # callable form
    # V4: the integral behind the criterion — c_i * int (beta - U'') |phi|^2/|U - c|^2 dy = 0 — on the shooting mode
    cs, (sol, kap) = _rk_shoot(1.4, 0.3, 0.458 + 0.053j, want_mode=True)
    yy = np.linspace(-14.0, 0.0, 28001)
    ph = sol.sol(yy)[0]
    w = (0.3 - SECH2["Upp"](yy)) * np.abs(ph) ** 2 / np.abs(SECH2["U"](yy) - cs) ** 2
    tail = lambda sgn: (0.3 * abs(ph[0]) ** 2 / abs(cs) ** 2) / (2 * kap.real)  # noqa: E731   analytic far-field part
    integral = np.sum(0.5 * (w[1:] + w[:-1]) * np.diff(yy)) + tail(1)
    assert abs(integral) < 2e-6 * np.sum(0.5 * (np.abs(w[1:]) + np.abs(w[:-1])) * np.diff(yy))


def test_rayleigh_contour_V7_beta_default_unchanged_and_radiating_mode_limit():  # V7 + V3, ST.rayleigh_eigs_contour(beta=)
    assert inspect.signature(ST.rayleigh_eigs_contour).parameters["beta"].default == 0.0
    kw = dict(bc="decay", y_max=16.0, parity="even")
    a = ST.rayleigh_eigs_contour(0.9, SECH2["U"], SECH2["Up"], SECH2["Upp"], **kw)
    b = ST.rayleigh_eigs_contour(0.9, SECH2["U"], SECH2["Up"], SECH2["Upp"], beta=0.0, **kw)
    assert np.array_equal(a, b) and len(a) == 1                                # the ch11 path is untouched
    cs, res = _rk_shoot(0.9, 0.0, a[0])
    assert res < 1e-10 and abs(a[0] - cs) < 1e-7                               # and still right: independent shooting
    for beta in (0.1, 0.2, 0.3):
        r = ch13.rayleigh_kuo_eigs(0.9, SECH2["U"], SECH2["Up"], SECH2["Upp"], beta, **kw)
        assert abs(r["c"] - _rk_shoot(0.9, beta, r["c"])[0]) < 1e-7
    assert ch13.rayleigh_kuo_eigs(0.9, SECH2["U"], SECH2["Up"], SECH2["Upp"], 0.45, **kw)["stable"] is True
    # k = 0.9, beta = -1: the far field is a weakly damped radiating wave, kappa = sqrt(k^2 + beta/c), Re kappa = 0.07,
    # so a box of half-width 16 with phi = 0 at its ends is not converged; the shooting value (radiation condition) is
    c14, r14 = _rk_shoot(0.9, -1.0, 0.88 + 0.064j, Y=14.0)
    c24, r24 = _rk_shoot(0.9, -1.0, 0.88 + 0.064j, Y=24.0)
    assert r14 < 1e-10 and abs(c14 - c24) < 1e-10                              # independent of where the shooting starts
    kap = np.sqrt(0.81 - 1.0 / c14)
    assert abs(kap.real) < 0.1                                                 # e-folding length 14: longer than the box
    box = ch13.rayleigh_kuo_eigs(0.9, SECH2["U"], SECH2["Up"], SECH2["Upp"], -1.0, **kw)
    assert abs(box["growth_rate"] / (0.9 * c14.imag) - 1) < 0.08               # measured: 0.0546 against 0.0575, 5 % low
    assert abs(0.9 * c14.imag - 0.057498) < 5e-6


# =====================================================================================================================
# 9. §13.17 — the Eady problem (C16)
# =====================================================================================================================
def test_eady_equation_V2_derivation():  # V2, D25 (13.128), (13.130)–(13.136), (13.139)–(13.140)
    x, y, z, t = sp.symbols("x y z t", real=True)
    f, N, rho0, U0, H, g, k, l = sp.symbols("f N rho_0 U_0 H g k l", positive=True)
    c = sp.Symbol("c")
    U = U0 * z / H
    rhobar = rho0 * (1 - N ** 2 * z / g) + f * rho0 * U0 / (g * H) * y          # uniform N, uniform shear
    assert sp.simplify(U.diff(z) - g / (f * rho0) * rhobar.diff(y)) == 0        # (13.128): thermal wind of the basic state
    assert sp.simplify(-g / rho0 * rhobar.diff(z) - N ** 2) == 0
    p = sp.Function("p")(x, y, z, t)
    adv = lambda q: q.diff(t) + U * q.diff(x)  # noqa: E731
    rho_p, v_p = -p.diff(z) / g, p.diff(x) / (rho0 * f)                         # (13.134) hydrostatic, (13.131) geostrophic
    w = sp.Symbol("w")
    w_sol = sp.solve(adv(rho_p) + v_p * rhobar.diff(y) + w * rhobar.diff(z), w)[0]   # (13.133): linearised density equation
    e135 = -(adv(p.diff(z)) - U.diff(z) * p.diff(x)) / (rho0 * N ** 2)          # (13.135)
    assert sp.simplify(w_sol - e135) == 0
    zeta = (p.diff(x, 2) + p.diff(y, 2)) / (rho0 * f)                           # (13.132)
    e130 = adv(zeta) - f * e135.diff(z)                                         # (13.130)
    e136 = adv(p.diff(x, 2) + p.diff(y, 2) + f ** 2 / N ** 2 * p.diff(z, 2)) / (rho0 * f)   # (13.136)
    assert sp.simplify(sp.expand(e130 - e136)) == 0                             # the two shear terms cancel under d/dz
    al = N * sp.sqrt(k ** 2 + l ** 2) / f                                       # (13.139)
    A, B = sp.symbols("A B")
    ph = A * sp.cosh(al * (z - H / 2)) + B * sp.sinh(al * (z - H / 2))          # (13.140)
    mode = ph * sp.exp(sp.I * (k * x + l * y - k * c * t))
    assert sp.simplify(e136.subs(p, mode).doit()) == 0
    assert ch13.eady_qg_sympy()["ok"]
    bs = ch13.eady_basic_state(np.linspace(-1e6, 1e6, 5), np.linspace(0.0, 9e3, 7), N=1.1e-2, f=F_S, H=9e3, U0=27.0, rho0=1.2)
    assert close(bs["U"], 27.0 * np.linspace(0, 1, 7)) and math.isclose(bs["drho_dy"], F_S * 1.2 * 27.0 / (G0 * 9e3))
    assert math.isclose(bs["slope"], F_S * 27.0 / (1.1e-2 ** 2 * 9e3))          # isopycnal slope = -rho_y/rho_z
    assert math.isclose(-(bs["rho"][0, 4] - bs["rho"][0, 0]) / 2e6 / ((bs["rho"][6, 0] - bs["rho"][0, 0]) / 9e3), bs["slope"], rel_tol=1e-9)
    assert math.isclose(ch13.eady_alpha(3e-6, 4e-6, 1.1e-2, F_S), 1.1e-2 * 5e-6 / F_N)


def test_eady_phase_speed_V2_derivation():  # V2, D26: w' = 0 on both lids -> the 2 x 2 system -> (13.141)
    z = sp.Symbol("z", real=True)
    al, U0, H = sp.symbols("alpha U_0 H", positive=True)
    c, A, B = sp.symbols("c A B")
    ph = A * sp.cosh(al * (z - H / 2)) + B * sp.sinh(al * (z - H / 2))
    # (13.135) with p' = phat e^{ik(x - ct)}: w' = 0 becomes (U - c) phat' - (U0/H) phat = 0 on each lid
    bc = lambda zz: ((U0 * zz / H - c) * ph.diff(z) - U0 / H * ph).subs(z, zz)  # noqa: E731
    q = al * H / 2
    row0 = [sp.expand(bc(0)).coeff(s) for s in (A, B)]
    rowH = [sp.expand(bc(H)).coeff(s) for s in (A, B)]
    page0 = [al * c * sp.sinh(q) - U0 / H * sp.cosh(q), -al * c * sp.cosh(q) + U0 / H * sp.sinh(q)]          # as printed
    pageH = [al * (U0 - c) * sp.sinh(q) - U0 / H * sp.cosh(q), al * (U0 - c) * sp.cosh(q) - U0 / H * sp.sinh(q)]
    assert all(sp.simplify(a - b) == 0 for a, b in zip(row0, page0)) and all(sp.simplify(a - b) == 0 for a, b in zip(rowH, pageH))
    det = sp.expand(sp.Matrix([page0, pageH]).det())
    quad_c = sp.Poly(det, c)
    assert quad_c.degree() == 2
    # the quadratic's roots against (13.141), to 30 digits at exact rational points on both sides of the cut-off
    for aH in (sp.Rational(1, 3), sp.Rational(8, 5), sp.Rational(23, 10), sp.Rational(5, 2), sp.Integer(4)):
        vals = {al: aH, H: 1, U0: sp.Rational(7, 3)}
        qa, qb, qc = (sp.N(co.subs(vals), 40) for co in quad_c.all_coeffs())    # the quadratic a c^2 + b c + c0 = 0, 40 digits
        disc = sp.sqrt(sp.N(qb ** 2 - 4 * qa * qc, 40))
        roots = [(-qb + s_ * disc) / (2 * qa) for s_ in (1, -1)]
        x = aH / 2
        rad = sp.sqrt(sp.N((x - sp.tanh(x)) * (x - sp.coth(x)), 40))
        target = [sp.Rational(7, 6) + s * sp.Rational(7, 3) / aH * rad for s in (1, -1)]
        for tgt in target:
            assert min(abs(sp.N(r - tgt, 30)) for r in roots) < sp.Float("1e-24")
    # the coded matrix is the page's, and it is singular exactly on the coded phase speed
    for aH in (0.4, 1.6, 2.39, 2.6, 3.7):
        for root in (1, -1):
            cc = ch13.eady_phase_speed(aH, 27.0, root)
            M = ch13.eady_matrix(cc, aH, 27.0, 9e3)
            num = np.array([[complex(e.subs({al: aH / 9e3, H: 9e3, U0: 27.0, c: cc})) for e in page0],
                            [complex(e.subs({al: aH / 9e3, H: 9e3, U0: 27.0, c: cc})) for e in pageH]])
            assert close(M, num, rtol=1e-12) and abs(np.linalg.det(M)) < 1e-12 * np.linalg.norm(M) ** 2
        assert abs(np.linalg.det(ch13.eady_matrix(ch13.eady_phase_speed(aH, 27.0) * 1.05, aH, 27.0, 9e3))) > 1e-4 * np.linalg.norm(M) ** 2


def test_eady_cutoff_and_fastest_wave_V2_derivation():  # V2 + V1, D27 (13.142): claim 4, our own root-finding
    xc = brentq(lambda s: s - 1 / math.tanh(s), 0.5, 3.0, xtol=1e-15)          # alpha_c H / 2 = coth(alpha_c H / 2)
    assert abs(ch13.eady_critical() - 2 * xc) < 1e-12 and abs(2 * xc - 2.39936) < 5e-6
    sig = lambda a: -a * math.sqrt(max(-(a / 2 - math.tanh(a / 2)) * (a / 2 - 1 / math.tanh(a / 2)), 0.0)) / a  # noqa: E731
    r = minimize_scalar(sig, bounds=(0.05, 2 * xc), method="bounded", options=dict(xatol=1e-13))
    fast = ch13.eady_fastest()
    assert abs(fast["alphaH"] - r.x) < 1e-6 and abs(fast["sigma_nd"] + r.fun) < 1e-12
    assert abs(fast["alphaH"] - 1.60612) < 5e-6 and abs(fast["sigma_nd"] - 0.30982) < 5e-6 and fast["cr_over_U0"] == 0.5
    assert math.isclose(fast["ci_over_U0"], fast["sigma_nd"] / fast["alphaH"], rel_tol=1e-12)
    fac = ch13.eady_factors(np.array([0.5, 2 * xc, 3.0]))
    assert np.all(fac["tanh_factor"] > 0) and fac["coth_factor"][0] < 0 < fac["coth_factor"][2] and abs(fac["coth_factor"][1]) < 1e-12
    assert close(fac["product"], fac["tanh_factor"] * fac["coth_factor"])
    small = ch13.eady_factors(1e-4)                                             # series branch joins the direct form
    assert math.isclose(small["product"], -(0.5e-4) ** 2 / 3, rel_tol=1e-6)
    assert math.isclose(ch13.eady_phase_speed(0.0, 1.0).imag, 1 / (2 * math.sqrt(3)), rel_tol=1e-12)   # long-wave limit of c_i/U0
    wl = ch13.eady_wavelengths()
    assert math.isclose(wl["cutoff_over_Lambda"], PI / xc) and math.isclose(wl["fastest_over_Lambda"], 2 * PI / fast["alphaH"])
    # the unstable band of (13.142) with l = 0: growth exactly for N H k/|f| below the cut-off
    N, H, U0 = I["atm_N"], I["atm_H"], I["atm_U0"]
    for f in HEMI:
        kc = 2 * xc * abs(f) / (N * H)
        ks = kc * np.array([0.2, 0.67, 0.99, 1.01, 1.6])
        g_ = ch13.eady_growth_rate(ks, 0.0, N, f, H, U0)
        assert np.all(g_[:3] > 0) and np.all(g_[3:] == 0)
        kf = fast["alphaH"] * abs(f) / (N * H)
        assert math.isclose(ch13.eady_growth_rate(kf, 0.0, N, f, H, U0), fast["sigma_nd"] * abs(f) * U0 / (N * H), rel_tol=1e-9)
        assert math.isclose(ch13.eady_max_growth_rate(f, N, U0 / H), fast["sigma_nd"] * abs(f) * U0 / (N * H), rel_tol=1e-12)
        assert math.isclose(ch13.eady_time_scale(f, N, U0 / H), 1 / ch13.eady_max_growth_rate(f, N, U0 / H))
        assert ch13.eady_growth_rate(kf, 0.6 * kf, N, f, H, U0) < ch13.eady_growth_rate(kf, 0.0, N, f, H, U0)   # l = 0 grows fastest
    const = json.loads((REF / "explainer_constants.json").read_text(encoding="utf-8"))
    assert abs(const["eady_critical"] - 2 * xc) < 1e-12 and abs(const["eady_fastest"]["sigma_nd"] + r.fun) < 1e-12


def test_eady_V5_published_growth_rate_and_wavenumbers():  # V5 benchmark (Emanuel, MIT OCW 12.803 Lecture 19; NCL)
    bm = bench()
    fast = ch13.eady_fastest()
    for key, ours in (("eady_max_growth_coefficient", fast["sigma_nd"]), ("eady_fastest_wavenumber", fast["alphaH"]),
                      ("eady_cutoff_wavenumber_rounded", ch13.eady_critical())):
        b = bm[key]
        assert b["verified"] == "2026-10-07" and b["url"].startswith("https://")
        half_unit = 0.5 * 10.0 ** (-(b["digits"] - 1) + math.floor(math.log10(b["value"])))
        assert abs(ours - b["value"]) <= half_unit, (key, ours)                 # agreement to every digit the source prints


def test_eady_V3_chebyshev_eigenvalues_against_the_closed_form():  # V3 + V1, independent numerical route (claim 4)
    N, H, U0 = I["atm_N"], I["atm_H"], I["atm_U0"]
    for f in HEMI:
        for aH in (0.4, 1.0, 1.6061, 2.2, 2.39, 2.6, 3.5):
            k = aH * abs(f) / (N * H)
            ne = ch13.eady_numeric_eigs(k, 0.0, N, f, H, U0, n=48)
            ca = ch13.eady_phase_speed(aH, U0)
            assert abs(ne["c"] - ca) < 2e-8 * U0
            assert math.isclose(ne["growth_rate"], k * max(ca.imag, 0.0), rel_tol=1e-6, abs_tol=1e-12 * abs(f))
    k = 1.6061 * F_N / (N * H)
    errs = [abs(ch13.eady_numeric_eigs(k, 0.0, N, F_N, H, U0, n=n)["c"] - ch13.eady_phase_speed(1.6061, U0)) for n in (4, 6, 8, 12)]
    assert errs[0] > 30 * errs[1] > 900 * errs[2] and errs[3] < 1e-8 * U0       # spectral: faster than any power
    lo, hi = ch13.eady_phase_speed(3.0, U0, -1), ch13.eady_phase_speed(3.0, U0, 1)
    assert lo.imag == 0 and hi.imag == 0 and 0 < lo.real < U0 / 2 < hi.real < U0   # two neutral edge waves beyond the cut-off
    assert math.isclose(lo.real + hi.real, U0)


def test_eady_V1_mode_structure_vertical_velocity_and_fluxes():  # V1 + V7, N160, N163
    N, H, U0, rho0 = I["atm_N"], I["atm_H"], I["atm_U0"], 1.2
    aH = ch13.eady_fastest()["alphaH"]
    z = np.linspace(0.0, H, 201)
    for f in HEMI:
        k = aH * abs(f) / (N * H)
        md = ch13.eady_mode(z, aH, U0, H)
        assert math.isclose(np.abs(md["p_hat"]).max(), 1.0, rel_tol=1e-9) and close(md["amplitude"], np.abs(md["p_hat"]))
        assert np.all(np.diff(md["phase"]) > 0) and math.isclose(md["phase"][-1] - md["phase"][0], PI / 2, rel_tol=1e-3)   # westward tilt
        assert close(md["amplitude"], md["amplitude"][::-1], rtol=1e-9)   # symmetric about mid-depth
        al = aH / H
        pzz = np.gradient(np.gradient(md["p_hat"], z), z)
        assert np.max(np.abs(pzz - al ** 2 * md["p_hat"])[3:-3]) < 2e-3 * al ** 2   # phat'' = alpha^2 phat
        w = ch13.eady_vertical_velocity(z, aH, U0, H, k, N, f, rho0)
        assert abs(w[0]) < 1e-12 * np.abs(w).max() and abs(w[-1]) < 1e-12 * np.abs(w).max()    # w' = 0 on both lids
        decay = ch13.eady_vertical_velocity(z, aH, U0, H, k, N, f, rho0, root=-1)
        assert abs(decay[0]) < 1e-12 * np.abs(decay).max()                      # the decaying partner satisfies them too
        fl = ch13.eady_fluxes(z, aH, U0, H, N, f, rho0)
        assert fl["v_T_sign"] == 1.0 and math.isclose(fl["phase_tilt"], PI / 2, rel_tol=1e-3)
        assert np.all(fl["w_rho"][1:-1] < 0)                                    # light fluid rises: potential energy is released
        assert np.all(np.sign(fl["v_rho"]) == -np.sign(f))                      # heat goes poleward in both hemispheres
        assert close(fl["v_rho"], fl["v_rho"][0], rtol=1e-9)              # independent of height (no interior PV flux)
        neutral = ch13.eady_fluxes(z, 3.0, U0, H, N, f, rho0)
        assert neutral["v_T_sign"] == 0.0 and np.max(np.abs(neutral["v_rho"])) < 1e-12 * np.max(np.abs(fl["v_rho"]))
    with pytest.raises(ValueError):
        ch13.eady_mode(z, 0.0, U0, H)


# =====================================================================================================================
# 10. §13.18 — two-dimensional turbulence (C17)
# =====================================================================================================================
def test_fjortoft_V2_two_constraints_and_ratios():  # V2 + V4, D28 (13.145)
    K0, K1, K2, S0, S1, S2 = sp.symbols("K0 K1 K2 S0 S1 S2", positive=True)
    sol = sp.solve([sp.Eq(S0, S1 + S2), sp.Eq(K0 ** 2 * S0, K1 ** 2 * S1 + K2 ** 2 * S2)], [S1, S2], dict=True)[0]
    assert sp.simplify(sol[S1] / sol[S2] - (K2 - K0) * (K2 + K0) / ((K0 - K1) * (K1 + K0))) == 0          # (13.145a)
    assert sp.simplify(K1 ** 2 * sol[S1] / (K2 ** 2 * sol[S2]) - K1 ** 2 / K2 ** 2 * (K2 ** 2 - K0 ** 2) / (K0 ** 2 - K1 ** 2)) == 0
    for K0n, K1n, K2n in ((1.0, 0.4, 3.0), (7.0, 5.5, 7.9), (2.0, 1.0, 2.2)):
        r = ch13.fjortoft_transfer(K0n, K1n, K2n, S0=2.5)
        assert math.isclose(r["S1"] + r["S2"], 2.5) and math.isclose(K1n ** 2 * r["S1"] + K2n ** 2 * r["S2"], K0n ** 2 * 2.5)
        assert math.isclose(r["energy_ratio"], float(sol[S1].subs({K0: K0n, K1: K1n, K2: K2n, S0: 1}) / sol[S2].subs({K0: K0n, K1: K1n, K2: K2n, S0: 1})))
        assert r["S1"] > 0 and r["S2"] > 0
    octave = ch13.fjortoft_transfer(1.0, 1 / 3, 3.0)                            # a symmetric (in log K) triad, ours
    assert octave["energy_ratio"] > 1 > octave["enstrophy_ratio"]               # energy to large scales, enstrophy to small
    with pytest.raises(ValueError):
        ch13.fjortoft_transfer(1.0, 2.0, 3.0)


def test_cascade_spectra_V2_derivation():  # V2, D29: K^2 S(K), the -3 and -5/3 ranges by dimensions, the Rhines length
    # exponents of [length, time] for S [L^3 T^-2], K [L^-1], enstrophy flux [T^-3], energy flux [L^2 T^-3]
    a, b = sp.symbols("a b")
    ens = sp.solve([sp.Eq(3, -b), sp.Eq(-2, -3 * a)], [a, b])                   # S = C alpha^a K^b
    assert (ens[a], ens[b]) == (sp.Rational(2, 3), -3)
    en = sp.solve([sp.Eq(3, 2 * a - b), sp.Eq(-2, -3 * a)], [a, b])             # S = C eps^a K^b
    assert (en[a], en[b]) == (sp.Rational(2, 3), sp.Rational(-5, 3))
    K = np.geomspace(0.5, 200.0, 400)
    S = ch13.two_d_cascade_spectrum(K, 9.0, 2.0, 5.0)
    slope = np.gradient(np.log(S), np.log(K))
    assert close(slope[K < 8.0][2:-2], -5 / 3, atol=1e-9) and close(slope[K > 10.0][2:-2], -3.0, atol=1e-9)
    assert math.isclose(ch13.two_d_cascade_spectrum(4.0, 9.0, 2.0, 5.0, C_E=1.7), 1.7 * 2.0 ** (2 / 3) * 4.0 ** (-5 / 3))
    assert math.isclose(ch13.two_d_cascade_spectrum(20.0, 9.0, 2.0, 5.0, C_Z=0.6), 0.6 * 5.0 ** (2 / 3) * 20.0 ** (-3))
    assert close(ch13.enstrophy_spectrum(K, S), K ** 2 * S)               # vorticity = curl u: a factor K^2 in spectra
    assert (Q_(1.0, "m**3/s**2") * Q_(1.0, "1/m") ** 2 * Q_(1.0, "1/m")).check("1/[time]**2")   # int K^2 S dK is an enstrophy
    u, beta = I["u_rms"][0], float(ch13.beta_parameter(LAT))
    Lr = ch13.rhines_length(u, beta)
    assert math.isclose(Lr, math.sqrt(u / beta)) and (Q_(1.0, "m/s") / Q_(1.0, "1/(m*s)")).check("[length]**2")
    assert math.isclose(u / Lr, beta * Lr, rel_tol=1e-12)                       # eddy turnover rate = Rossby-wave frequency there
    assert math.isclose(abs(ch13.rossby_omega(1 / Lr, 0.0, beta)), u / Lr, rel_tol=1e-12)


# =====================================================================================================================
# 11. Numerical models (our choices of scheme): C-grid shallow water, the 1-D line, barotropic vorticity (claim 7)
# =====================================================================================================================
def _poincare_state(md, t, k, l, H, f, amp=0.3):
    g = md.grids()
    e = ch13.poincare_fields(g["center"][0], g["center"][1], t, k, l, amp, H, f)[0]
    u = ch13.poincare_fields(g["u"][0], g["u"][1], t, k, l, amp, H, f)[1]
    v = ch13.poincare_fields(g["v"][0], g["v"][1], t, k, l, amp, H, f)[2]
    return SW.make_state(md, e, u, v)


def test_cgrid_V3_poincare_wave_second_order_in_space_third_in_time():  # V3, C07/C10 (13.45): SW.run, SW.step
    H, L = I["ocean_H"], 4.0e6
    k, l = 2 * PI / L, 2 * 2 * PI / L
    for f in HEMI:
        T = 2 * PI / ch13.poincare_omega(math.hypot(k, l), f, math.sqrt(G0 * H))
        hs, errs, drift = [], [], []
        for n in (16, 32, 64):
            md = ch13.ShallowWater(n, n, L, L, H, f)
            s0 = _poincare_state(md, 0.0, k, l, H, f)
            r = ch13.run(md, s0, T, dt=ch13.dt_limit(md, 0.5), save_every=10 ** 6)
            errs.append(np.max(np.abs(r["eta"][-1] - _poincare_state(md, T, k, l, H, f)["eta"])) / 0.3)
            drift.append(r["energy"][-1] / r["energy"][0] - 1)
            hs.append(md.dx)
            assert r["t"][-1] == pytest.approx(T) and r["eta"].shape == (2, n, n)
            assert abs(SW.volume(md, dict(eta=r["eta"][-1], u=r["u"][-1], v=r["v"][-1])) - SW.volume(md, s0)) < 1e-12 * 0.3 * L * L
        assert abs(observed_order(hs, errs) - 2.0) < 0.15 and errs[-1] < 0.01            # measured 1.98; 0.9 % on 64^2
        assert all(d < 0 for d in drift) and abs(drift[-1]) < 2e-4                        # energy drift of the time stepping
    md = ch13.ShallowWater(16, 16, L, L, H, F_N)
    s0 = _poincare_state(md, 0.0, k, l, H, F_N)
    dtm = ch13.dt_limit(md, 1.0)
    assert math.isclose(dtm, math.sqrt(3) / math.sqrt(F_N ** 2 + 4 * G0 * H * (2 / md.dx ** 2))) and math.isclose(ch13.dt_limit(md), 0.5 * dtm)
    ref = ch13.run(md, s0, 40 * dtm, dt=dtm / 32, save_every=10 ** 6)["eta"][-1]
    dts, e_t, e_E = [], [], []
    for div in (1, 2, 4, 8):
        rr = ch13.run(md, s0, 40 * dtm, dt=dtm / div, save_every=10 ** 6)
        e_t.append(np.max(np.abs(rr["eta"][-1] - ref)))
        e_E.append(abs(rr["energy"][-1] / rr["energy"][0] - 1))
        dts.append(dtm / div)
        assert np.isfinite(rr["energy"][-1]) and rr["energy"][-1] <= rr["energy"][0] * (1 + 1e-12)   # stable up to the limit
    assert abs(observed_order(dts, e_t) - 3.0) < 0.15 and abs(observed_order(dts, e_E) - 3.0) < 0.15   # SSP-RK3
    with pytest.raises(ValueError, match="stability limit"):
        ch13.step(md, s0, 1.01 * dtm)
    one = ch13.step(md, s0, 0.5 * dtm)
    assert set(one) == {"eta", "u", "v"} and one["eta"].shape == (16, 16)


def test_cgrid_V4_volume_and_energy_conserved_by_the_space_discretisation():  # V4, C07 (13.44), (13.88)–(13.90)
    rng = np.random.default_rng(1)
    for bc in ("periodic", "channel", "closed"):
        md = ch13.ShallowWater(12, 10, 1.0e6, 8.0e5, 500.0, 1.0e-4, beta=(2e-11 if bc != "periodic" else 0.0), bc=bc)
        st = SW.make_state(md, 0.5 * rng.standard_normal((10, 12)), 0.1 * rng.standard_normal((10, 12)), 0.1 * rng.standard_normal((10, 12)))
        if bc != "periodic":
            assert np.all(st["v"][0] == 0.0)                                    # wall-normal velocity removed
        if bc == "closed":
            assert np.all(st["u"][:, 0] == 0.0)
        for linear in (True, False):
            de = ch13.continuity_tendency(md, st, linear)
            du, dv = ch13.momentum_tendencies(md, st, linear)
            assert abs(np.sum(de)) < 1e-13 * np.sum(np.abs(de))                 # flux form: volume to round-off
            a = 1e-3
            E = lambda s: ch13.energy(md, dict(eta=st["eta"] + s * de, u=st["u"] + s * du, v=st["v"] + s * dv), linear)["total"]  # noqa: E731
            dE = (E(a) - E(-a)) / (2 * a)                                       # dE/dt along the tendencies
            scale = math.sqrt(np.sum(de ** 2) * np.sum(st["eta"] ** 2)) * md.g * md.dx * md.dy
            assert abs(dE) < 1e-8 * scale                                       # energy-conserving Coriolis/vorticity average
        en = ch13.energy(md, st)
        assert math.isclose(en["total"], en["kinetic"] + en["potential"]) and math.isclose(en["potential"], 0.5 * md.g * np.sum(st["eta"] ** 2) * md.dx * md.dy)
    with pytest.raises(ValueError):
        ch13.ShallowWater(8, 8, 1e6, 1e6, 100.0, 1e-4, beta=2e-11)              # a beta-plane cannot be periodic in y
    with pytest.raises(ValueError):
        ch13.ShallowWater(8, 8, 1e6, 1e6, 100.0, 1e-4, bottom=150.0)


def test_cgrid_V1_balanced_states_vorticity_and_kelvin_wave_both_hemispheres():  # V1 + V3 + V7, C02, C11, C13
    hs, es = [], []
    for n in (32, 64, 128):
        md = ch13.ShallowWater(n, n, 4.0e6, 4.0e6, 300.0, F_S)
        bump = ch13.gaussian_bump(md, 1.0, radius=4.0e5)
        st = ch13.geostrophic_state(md, bump)
        du, dv = ch13.momentum_tendencies(md, st)
        es.append(max(np.max(np.abs(du)), np.max(np.abs(dv))))
        hs.append(md.dx)
    assert abs(observed_order(hs, es) - 2.0) < 0.15                             # geostrophic balance holds to O(dx^2)
    assert math.isclose(bump.max(), np.exp(-(md.dx / 2) ** 2 / 4.0e5 ** 2), rel_tol=1e-12) and bump.shape == (128, 128)
    zeta = ch13.relative_vorticity(md, st)
    Xq, Yq = md.grids()["corner"]
    r2 = ((Xq - 2.0e6) ** 2 + (Yq - 2.0e6) ** 2) / 4.0e5 ** 2
    exact = G0 / F_S * (r2 - 2) / 4.0e5 ** 2 * np.exp(-r2 / 2)                  # (g/f) laplacian of the bump (13.116)
    assert np.max(np.abs(zeta - exact)) < 0.01 * np.max(np.abs(exact))
    assert np.sign(zeta[64, 64]) == np.sign(-F_S)                               # a high is anticyclonic: zeta opposite to f
    q = ch13.sw_potential_vorticity(md, st)
    hq = 300.0 + 0.25 * (st["eta"] + np.roll(st["eta"], 1, 1) + np.roll(st["eta"], 1, 0) + np.roll(np.roll(st["eta"], 1, 0), 1, 1))
    assert close(q, (zeta + F_S) / hq, rtol=1e-13)                        # (13.94) on the grid
    ql = ch13.sw_potential_vorticity(md, st, linear=True)
    assert np.max(np.abs(ql - q)) < 2e-4 * np.max(np.abs(q))                    # linear form: O(eta/H)^2 apart
    # a Kelvin wave leaning on the south wall of a channel returns to itself after one period, in each hemisphere
    for f in HEMI:
        ch = ch13.ShallowWater(64, 64, 2.0e6, 1.2e6, 20.0, f, bc="channel")
        ks = SW.kelvin_state(ch, 0.2)
        c = math.sqrt(G0 * 20.0)
        r = ch13.run(ch, ks, 2.0e6 / c, save_every=10 ** 6)
        assert np.max(np.abs(r["eta"][-1] - ks["eta"])) < 4e-3 * 0.2            # measured 0.24 %
        assert np.max(np.abs(r["v"][-1])) < 1e-3 * np.max(np.abs(r["u"][-1]))   # no flow across the channel
        quarter = ch13.run(ch, ks, 0.25 * 2.0e6 / c, save_every=10 ** 6)["eta"][-1]
        shift = int(round(np.sign(f) * 16))                                     # a quarter wavelength toward sign(f) x
        assert np.max(np.abs(quarter - np.roll(ks["eta"], shift, axis=1))) < 6e-3 * 0.2
        assert np.argmax(np.abs(ks["eta"]).max(axis=1)) == 0                    # largest at the wall, decays as exp(-y/Lambda)
    with pytest.raises(ValueError):
        SW.kelvin_state(ch13.ShallowWater(8, 8, 1e6, 1e6, 10.0, 1e-4), 0.1)     # needs a wall


def test_particles_V1_interpolation_advection_and_pv_on_paths():  # V1 + V4, C13 (13.94): advect_particles
    md = ch13.ShallowWater(16, 12, 1.6e6, 1.2e6, 100.0, F_N)
    g = md.grids()
    for loc in ("center", "u", "v", "corner"):
        X, Y = g[loc]
        a = 2.0 + 3e-6 * X - 1e-6 * Y + 4e-12 * X * Y                          # bilinear: reproduced exactly
        xp, yp = np.array([4.3e5, 9.9e5]), np.array([3.1e5, 7.7e5])
        assert close(SW.interpolate(md, a, xp, yp, loc), 2.0 + 3e-6 * xp - 1e-6 * yp + 4e-12 * xp * yp, rtol=1e-13)
    t = np.linspace(0.0, 4.0e4, 9)
    hist = dict(t=t, u=np.full((9, 12, 16), 1.5), v=np.full((9, 12, 16), -0.5))
    xs, ys = ch13.advect_particles(md, hist, np.array([3.0e5, 8.0e5]), np.array([6.0e5, 4.0e5]))
    assert close(xs[-1], np.array([3.0e5, 8.0e5]) + 1.5 * 4.0e4) and close(ys[-1], np.array([6.0e5, 4.0e5]) - 0.5 * 4.0e4)
    assert xs.shape == (9, 2)
    r = ch13.reference_run("pv_particles", fast=True)                           # nonlinear vortex, 32^2, one inertial period
    assert np.max(np.abs(r["q"] / r["q"][0] - 1)) < 0.02 and np.ptp(r["q"][0]) / np.mean(r["q"][0]) > 0.3   # measured 1.5 % of a 43 % spread
    full = ch13.load_reference_run("pv_particles")                              # the cached 48^2 run: tighter
    assert full is not None and np.max(np.abs(full["q"] / full["q"][0] - 1)) < 0.008
    assert np.max(np.hypot(full["xp"][-1] - full["xp"][0], full["yp"][-1] - full["yp"][0])) > 1.0e5   # the particles really moved
    st = SW.make_state(md)
    qq = SW.particle_potential_vorticity(md, dict(t=[0.0], eta=[st["eta"]], u=[st["u"]], v=[st["v"]]), np.array([[5e5]]), np.array([[5e5]]))
    assert math.isclose(qq[0, 0], F_N / 100.0, rel_tol=1e-13)                   # at rest: q = f/H


def test_linear_1d_V3_first_order_in_time_second_in_space():  # V3, curation note 7: linear_1d_step / linear_1d_run
    H, f, a = 60.0, 1.0e-4, 0.4
    c = math.sqrt(G0 * H)
    L = 6 * c / f
    kap = 2 * PI / L
    om = math.sqrt(f * f + c * c * kap * kap)
    T = 1.7 * 2 * PI / om

    def final(n, dt_target):                                                    # eta at time T on n cells
        dx = L / n
        nst = int(round(T / dt_target))
        x = (np.arange(n) + 0.5) * dx
        r = ch13.linear_1d_run(a * np.cos(kap * x), dx, T / nst, nst, H, f, save_every=10 ** 7)
        exact = a * (1 - (c * kap / om) ** 2 * (1 - math.cos(om * T))) * np.cos(kap * x)   # exact solution from rest (ours)
        return r["eta"][-1], exact

    dt0 = 0.8 * (L / 64) / c
    e_dt = [np.max(np.abs(final(64, dt0 / d)[0] - final(64, dt0 / (2 * d))[0])) for d in (1, 2, 4, 8)]
    assert abs(observed_order([dt0 / d for d in (1, 2, 4, 8)], e_dt) - 1.0) < 0.15        # forward–backward: first order
    hs, e_dx = [], []
    for n in (16, 32, 64):
        dts = 0.8 * (L / 64) / c / 4
        e1, ex = final(n, dts)
        e2, _ = final(n, dts / 2)
        e_dx.append(np.max(np.abs(2 * e2 - e1 - ex)) / a)                        # Richardson in dt removes the time error
        hs.append(L / n)
    assert abs(observed_order(hs, e_dx) - 2.0) < 0.15
    with pytest.raises(ValueError):
        ch13.linear_1d_step(np.zeros(8), np.zeros(9), np.zeros(8), 1.0, 1.01 / c, H, 0.0)   # CFL above 1 refused
    with pytest.raises(ValueError):
        ch13.linear_1d_step(np.zeros(8), np.zeros(8), np.zeros(8), 1.0, 0.1 / c, H, 0.0)
    # f = 0: a pulse splits in two that travel at sqrt(gH) (the non-rotating limit, ch07)
    n, dx = 400, 100.0
    x = (np.arange(n) + 0.5) * dx
    r = ch13.linear_1d_run(np.exp(-((x - 2.0e4) / 800.0) ** 2), dx, 0.5 * dx / c, 300, H, 0.0, save_every=10 ** 7)
    right = x > 2.0e4
    assert math.isclose(x[right][np.argmax(r["eta"][-1][right])] - 2.0e4, c * 300 * 0.5 * dx / c, rel_tol=0.02)
    assert math.isclose(r["eta"][-1].max(), 0.5, rel_tol=0.02) and np.all(r["v"] == 0.0)


def test_linear_1d_V4_adjustment_keeps_pv_and_reaches_the_end_state():  # V4 + V1, C12: claims 5 and 7
    H, eta0 = 1.3, 0.05
    c = math.sqrt(G0 * H)
    for f in HEMI:
        Lam = c / abs(f)
        n = 600
        dx = 120 * Lam / n
        x = (np.arange(n) + 0.5) * dx - 60 * Lam
        dt = 0.5 * dx / c
        nper = int(round(2 * PI / abs(f) / dt))
        r = ch13.linear_1d_run(eta0 * np.sign(x), dx, dt, 6 * nper, H, f)
        assert r["pv"].shape == (2, n - 1) and r["eta"].shape == (6 * nper + 1, n) and r["u"].shape[1] == n + 1
        assert np.max(np.abs(r["pv"][1] - r["pv"][0])) < 1e-12 * np.max(np.abs(r["pv"][0]))      # discrete linear PV: round-off
        assert abs(r["eta"][-1].sum() - r["eta"][0].sum()) < 1e-12 * eta0 * n                    # volume
        assert np.all(r["u"][:, 0] == 0.0) and np.all(r["u"][:, -1] == 0.0)                      # walls
        ee, ve = ch13.geostrophic_adjustment_1d(x, eta0, H, f)
        near = np.abs(x) < 3 * Lam
        err = []
        for j in (2, 4, 6):                                                      # averages over successive inertial periods
            em = r["eta"][(j - 1) * nper + 1:j * nper + 1].mean(axis=0)
            vm = r["v"][(j - 1) * nper + 1:j * nper + 1].mean(axis=0)
            err.append(max(np.max(np.abs(em - ee)[near]) / eta0, np.max(np.abs(vm - ve)[near]) / np.max(np.abs(ve))))
        assert err[0] > err[1] > err[2] and err[2] < 0.008                       # the transient radiates away (measured 0.6 %)
        assert np.sign(vm[n // 2]) == np.sign(f)                                 # jet along the step, direction set by sign(f)
        wide = np.abs(x) < 20 * Lam
        ratio = np.sum(0.5 * H * vm[wide] ** 2) / np.sum(0.5 * G0 * (eta0 ** 2 - em[wide] ** 2))
        assert abs(ratio - 1 / 3) < 0.01                                         # a third of the released energy stays (measured 0.335)


def test_barotropic_V4_energy_and_enstrophy_and_strict_dealiasing():  # V4 + V3, C17 (13.143)–(13.144): claim 7
    n, L = 64, 2 * PI
    z0 = SW.random_vorticity(n, L, seed=1)
    assert math.isclose(2 * ch13.barotropic_invariants(z0, L)["energy"], 1.0, rel_tol=1e-12) and abs(z0.mean()) < 1e-14
    assert np.array_equal(z0, SW.random_vorticity(n, L, seed=1)) and not np.array_equal(z0, SW.random_vorticity(n, L, seed=2))
    dt = SW.barotropic_time_step(z0, L)
    u0, v0 = SW.barotropic_velocity(z0, L)
    assert math.isclose(dt, 0.4 * (L / n) / np.max(np.hypot(u0, v0)))
    drift = []
    for div in (1, 2, 4):
        r = ch13.barotropic_run(n, L, z0, t_end=1.0, dt=dt / div, save_every=10 ** 6)
        drift.append((dt / div, abs(r["energy"][-1] / r["energy"][0] - 1), abs(r["enstrophy"][-1] / r["enstrophy"][0] - 1)))
    assert drift[0][1] < 5e-6 and drift[0][2] < 3e-5                            # inviscid, de-aliased: invariants (measured 2e-6, 1e-5)
    oE = observed_order([d[0] for d in drift], [d[1] for d in drift])
    oZ = observed_order([d[0] for d in drift], [d[2] for d in drift])
    assert oE > 3.85 and oZ > 3.85                                              # the drift is the RK4 error: at least dt^4 (measured 5.0, 4.9)
    # the truncated equations conserve both exactly: <psi, rhs> = 0 and <zeta, rhs> = 0 at any instant, with beta too
    KX, KY = SW.spectral_wavenumbers(n, L)
    zh = np.fft.rfft2(z0)
    rhs = np.fft.irfft2(ch13.barotropic_vorticity_rhs(zh, KX, KY, beta=3.0), s=(n, n))
    K2 = KX ** 2 + KY ** 2
    psi = np.fft.irfft2(-zh * np.where(K2 > 0, 1 / np.where(K2 > 0, K2, 1.0), 0.0), s=(n, n))
    assert abs(np.mean(z0 * rhs)) < 1e-13 * math.sqrt(np.mean(z0 ** 2) * np.mean(rhs ** 2))
    assert abs(np.mean(psi * rhs)) < 1e-13 * math.sqrt(np.mean(psi ** 2) * np.mean(rhs ** 2))
    with np.errstate(all="ignore"), warnings.catch_warnings():                 # the same run without the 2/3 rule
        warnings.simplefilter("ignore")
        raw = ch13.barotropic_run(n, L, z0, t_end=1.0, dt=dt, dealias=False, save_every=10 ** 6)
    bad = abs(raw["enstrophy"][-1] / raw["enstrophy"][0] - 1)
    assert not np.isfinite(bad) or bad > 100 * drift[0][2]                      # aliasing destroys the invariants (measured: blow-up)
    # strict '<' at the cut: modes (8, 0) and (8, 3) on 48^2 make a product at k_x = 16 = exactly 2/3 of Nyquist
    m = 48
    x = np.arange(m) * L / m
    X, Y = np.meshgrid(x, x)
    KX2, KY2 = SW.spectral_wavenumbers(m, L)
    zt = np.cos(8 * X) + np.cos(8 * X + 3 * Y + 0.3)
    kept = np.abs(ch13.barotropic_vorticity_rhs(np.fft.rfft2(zt), KX2, KY2))
    full = np.abs(ch13.barotropic_vorticity_rhs(np.fft.rfft2(zt), KX2, KY2, dealias=False))
    at_cut = np.isclose(np.abs(KX2), 16.0) & np.isclose(np.abs(KY2), 3.0)
    assert full[at_cut].max() > 0.1 * full.max() and kept[at_cut].max() == 0.0   # present in the product, removed by the filter
    assert kept[np.isclose(np.abs(KX2), 0.0) & np.isclose(np.abs(KY2), 3.0)].max() > 0.1 * full.max()   # the resolved partner stays


def test_barotropic_V1_single_modes_rossby_wave_viscous_decay_and_diagnostics():  # V1, C15/C17 (13.122)
    n, L = 32, 2 * PI
    x = np.arange(n) * L / n
    X, Y = np.meshgrid(x, x)
    kk, ll, beta = 3.0, 2.0, 5.0
    z0 = 1e-4 * np.cos(kk * X + ll * Y)
    r = ch13.barotropic_run(n, L, z0, t_end=2.0, dt=0.01, beta=beta, save_every=50)
    om = ch13.rossby_omega(kk, ll, beta)                                        # a single wave is an exact nonlinear solution
    assert np.max(np.abs(r["zeta"][-1] - 1e-4 * np.cos(kk * X + ll * Y - om * 2.0))) < 1e-8 * 1e-4
    assert r["t"].tolist() == pytest.approx([0.0, 0.5, 1.0, 1.5, 2.0]) and r["zeta"].shape == (5, n, n)
    rv = ch13.barotropic_run(n, L, z0, t_end=1.0, dt=0.01, nu=0.02, save_every=10 ** 6)
    assert math.isclose(rv["enstrophy"][-1] / rv["enstrophy"][0], math.exp(-2 * 0.02 * 13.0 * 1.0), rel_tol=1e-8)   # exp(-2 nu K^2 t)
    inv = ch13.barotropic_invariants(z0, L)
    assert math.isclose(inv["enstrophy"], 0.5e-8, rel_tol=1e-12) and math.isclose(inv["energy"], 0.25e-8 / 13.0, rel_tol=1e-12)
    cen = ch13.spectral_centroids(z0, L)
    assert math.isclose(cen["K_E"], math.sqrt(13.0), rel_tol=1e-12) and math.isclose(cen["K_Z"], math.sqrt(13.0), rel_tol=1e-12)
    u, v = SW.barotropic_velocity(z0, L)                                        # u = -dpsi/dy, v = dpsi/dx, laplacian(psi) = zeta
    assert close(u, -1e-4 * ll / 13.0 * np.sin(kk * X + ll * Y), atol=1e-18) and close(v, 1e-4 * kk / 13.0 * np.sin(kk * X + ll * Y), atol=1e-18)
    assert SW.zonal_energy_fraction(np.cos(4 * Y), L) == pytest.approx(1.0) and SW.zonal_energy_fraction(np.cos(4 * X), L) == pytest.approx(0.0, abs=1e-20)
    zr = SW.random_vorticity(64, L, seed=4)
    Kb, Eb = SW.barotropic_spectrum(zr, L)[:2]
    assert math.isclose(np.sum(Eb) * (Kb[1] - Kb[0]), ch13.barotropic_invariants(zr, L)["energy"], rel_tol=1e-10)   # sum S dK = energy
    c2 = ch13.spectral_centroids(zr, L)
    assert c2["K_Z"] > c2["K_E"]                                                # enstrophy sits at higher wavenumber than energy
    with pytest.raises(ValueError):
        ch13.barotropic_run(n, L, z0 * 1e6, t_end=1.0, dt=1.0)                  # above the advective limit
    with pytest.raises(ValueError):
        ch13.barotropic_run(31, L, z0[:31, :31], t_end=1.0, dt=0.01)


def test_cached_runs_V4_committed_files_show_the_physics():  # V4 + V7, reference/ch13/*.npz (our own output)
    kb = ch13.load_reference_run("kelvin_basin")
    c, Lam = float(kb["c"]), float(kb["Lambda"])
    assert math.isclose(c, math.sqrt(G0 * float(kb["H"]))) and math.isclose(Lam, c / abs(float(kb["f"]))) and float(kb["f"]) > 0
    dxg = kb["x"][1] - kb["x"][0]
    pos = [np.unravel_index(np.argmax(kb["eta"][i]), kb["eta"][i].shape) for i in (0, 2)]
    assert pos[0][0] == 0 and pos[1][0] == 0                                    # the pulse leans on the south wall ...
    speed = (kb["x"][pos[1][1]] - kb["x"][pos[0][1]]) / (kb["t"][2] - kb["t"][0])
    assert abs(speed - c) < 1.5 * dxg / (kb["t"][2] - kb["t"][0])               # ... and runs toward +x at sqrt(g H): coast on its right
    j5 = np.unravel_index(np.argmax(kb["eta"][5]), kb["eta"][5].shape)
    assert j5[1] == len(kb["x"]) - 1 and j5[0] > 5                              # then turns the corner and climbs the east wall
    across = kb["eta"][0][:, pos[0][1]]
    fit = np.polyfit(kb["y"][:12], np.log(across[:12]), 1)[0]
    assert math.isclose(-1 / fit, Lam, rel_tol=0.05)                            # offshore decay scale = Rossby radius
    assert abs(np.sum(kb["eta"][-1]) - np.sum(kb["eta"][0])) < 2e-3 * np.sum(np.abs(kb["eta"][0]))   # volume, to float16 storage
    tf, tb = ch13.load_reference_run("turbulence_f"), ch13.load_reference_run("turbulence_beta")
    for run_ in (tf, tb):
        assert np.all(np.diff(run_["energy"]) < 0) and np.all(np.diff(run_["enstrophy"]) < 0)      # viscous: both decay,
        assert run_["enstrophy"][-1] / run_["enstrophy"][0] < 0.5 * run_["energy"][-1] / run_["energy"][0]   # enstrophy much faster
        assert run_["K_E"][-1] < 0.5 * run_["K_E"][0]                           # energy moves to larger scales (inverse cascade)
        assert math.isclose(ch13.barotropic_invariants(run_["zeta"][0], float(run_["L"]))["energy"], run_["energy"][0], rel_tol=2e-3)
    zf = SW.zonal_energy_fraction(tf["zeta"][-1], float(tf["L"]))
    zb = SW.zonal_energy_fraction(tb["zeta"][-1], float(tb["L"]))
    assert float(tf["beta"]) == 0.0 and float(tb["beta"]) > 0 and zb > 2 * zf   # beta organises the flow into zonal bands
    assert ch13.load_reference_run("kelvin_basin", folder=ROOT / "no_such_folder") is None
    with pytest.raises(ValueError):
        ch13.load_reference_run("not_a_run")
    assert ch13.REFERENCE_RUNS == ("kelvin_basin", "turbulence_f", "turbulence_beta", "pv_particles")


def test_cached_runs_V7_two_files_regenerate_within_float16_storage():  # claim 8 (fast half; the rest is the slow test)
    for name in ("kelvin_basin", "pv_particles"):
        fresh, stored = ch13.reference_run(name), ch13.load_reference_run(name)
        assert set(fresh) == set(stored)
        for key in fresh:
            a, b = np.asarray(fresh[key], float), np.asarray(stored[key], float)
            assert a.shape == b.shape and np.max(np.abs(a - b)) <= 2e-3 * max(np.max(np.abs(a)), 1e-300), (name, key)
    const = ch13.explainer_constants()
    disk = json.loads((REF / "explainer_constants.json").read_text(encoding="utf-8"))
    assert disk["eady_fastest"] == const["eady_fastest"] and disk["inputs"]["wavelengths"] == list(I["wavelengths"])
    assert disk["omega_earth"] == ch13.OMEGA_EARTH and disk["earth_radius"] == ch13.EARTH_RADIUS_MEAN


@slow
def test_cached_runs_V7_all_files_regenerate_from_the_script(tmp_path):  # claim 8: scripts/ch13_make_caches.py
    r = subprocess.run([sys.executable, str(ROOT / "scripts" / "ch13_make_caches.py"), "--no-show", "--which", "all", "--dest", str(tmp_path),
                        "--out", str(tmp_path / "figs")], cwd=str(ROOT), capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=900)
    assert r.returncode == 0, r.stderr[-2000:]
    for name in ch13.REFERENCE_RUNS:
        new, old = ch13.load_reference_run(name, folder=tmp_path), ch13.load_reference_run(name)
        assert new is not None and set(new) == set(old)
        for key in new:
            a, b = np.asarray(new[key], float), np.asarray(old[key], float)
            assert a.shape == b.shape and np.max(np.abs(a - b)) <= 2e-3 * max(np.max(np.abs(a)), 1e-300), (name, key)
    assert json.loads((tmp_path / "explainer_constants.json").read_text(encoding="utf-8")) == json.loads(
        (REF / "explainer_constants.json").read_text(encoding="utf-8"))
    assert sum(f.stat().st_size for f in tmp_path.glob("*.npz")) < 2 * 1024 * 1024
    wr = ch13.write_reference_runs(["pv_particles"], tmp_path / "again")
    assert wr["pv_particles"]["path"].exists() and "explainer_constants" in wr


# =====================================================================================================================
# 12. Figure helpers, sketches, scripts
# =====================================================================================================================
def test_figures_V7_helpers_return_drawn_figures():  # smoke, C.4 figure helpers
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    made = [ch13.fig_ekman_surface(), ch13.fig_ekman_bottom(), ch13.fig_mode_roots(), ch13.fig_vertical_modes(),
            ch13.fig_poincare_kelvin_dispersion(), ch13.fig_kelvin_sections(), ch13.fig_inertia_gravity_orbit(),
            ch13.fig_rossby_dispersion()]
    for fig in made:
        assert isinstance(fig, matplotlib.figure.Figure) and len(fig.axes) >= 1
        assert sum(len(ax.lines) + len(ax.collections) for ax in fig.axes) >= 2
        assert all(ax.get_xlabel() != "" for ax in fig.axes if ax.lines)
    hodo = made[0].axes[0].lines[0]
    assert np.all(np.isfinite(hodo.get_xdata())) and len(hodo.get_xdata()) > 50
    hodo = ch13.fig_ekman_surface(lat_rad=LAT).axes[0].lines[0]
    south = ch13.fig_ekman_surface(lat_rad=-LAT).axes[0].lines[0]              # the southern spiral is the mirror image
    assert close(south.get_xdata(), hodo.get_xdata()) and close(south.get_ydata(), -np.asarray(hodo.get_ydata()))
    plt.close("all")


def test_drawings_V7_sketch_helpers_return_figures():  # smoke, C.4b scripts/ch13_drawings.py
    import importlib
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(ROOT / "scripts"))
    try:
        dr = importlib.import_module("ch13_drawings")
    finally:
        sys.path.remove(str(ROOT / "scripts"))
    names = ["draw_tangent_plane", "draw_parcel_balance", "draw_thermal_wind_wedge", "draw_taylor_column",
             "draw_ekman_layer_sketch", "draw_shallow_layer", "draw_mode_stack", "draw_kelvin_sections", "draw_pv_column",
             "draw_step_flow_sketch", "draw_eady_wedge", "draw_cascade_arrows"]
    for nm in names:
        fig = getattr(dr, nm)()
        assert isinstance(fig, matplotlib.figure.Figure) and len(fig.axes) >= 1, nm
    assert isinstance(dr.draw_parcel_balance(hemisphere=-1), matplotlib.figure.Figure)
    assert isinstance(dr.draw_ekman_layer_sketch(kind="bottom"), matplotlib.figure.Figure)
    plt.close("all")


@slow
def test_scripts_V7_run_clean(tmp_path):  # smoke: every scripts/ch13_*.py exits 0 (--no-show --fast)
    ran = 0
    for s in sorted((ROOT / "scripts").glob("ch13_*.py")):
        if s.name in ("ch13_common.py", "ch13_make_caches.py"):                 # the cache writer has its own test (with --dest)
            continue
        r = subprocess.run([sys.executable, str(s), "--no-show", "--fast", "--out", str(tmp_path)], cwd=str(ROOT), capture_output=True,
                           text=True, encoding="utf-8", errors="replace", timeout=900)
        assert r.returncode == 0, (s.name, r.stderr[-2000:])
        ran += 1
    assert ran == 13
    r = subprocess.run([sys.executable, str(REF / "make_refs.py")], cwd=str(ROOT), capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=120)
    assert r.returncode == 0, r.stderr[-2000:]


def test_reference_V5_sources_file_cites_every_file_and_no_book_table():  # data-and-benchmarks policy
    src = (REF / "SOURCES.md").read_text(encoding="utf-8")
    for name in ("kelvin_basin.npz", "turbulence_f.npz", "turbulence_beta.npz", "pv_particles.npz", "explainer_constants.json",
                 "benchmarks.json", "make_refs.py"):
        assert name in src and (REF / name).exists(), name
    assert "Emanuel" in src and "12.803" in src and "ncl.ucar.edu" in src and "2026-10-07" in src
    before = (REF / "benchmarks.json").read_bytes()
    bm = bench()
    assert bm["eady_max_growth_coefficient"]["value"] == 0.3098 and bm["eady_fastest_wavenumber"]["value"] == 1.606
    assert all("verified" in v and "source" in v for k, v in bm.items() if not k.startswith("_"))
    assert (REF / "benchmarks.json").read_bytes() == before


# =====================================================================================================================
# 13. Numbers the book prints (V6, private: skipped when tests/book_values_ch13.json is absent)
# =====================================================================================================================
@book_only
def test_book_V6_planetary_numbers_and_scales():  # V6 §13.2–§13.5
    b = book()
    s4 = b["sec_13_4_thin_layer"]
    assert rel(ch13.OMEGA_SOLAR_DAY, s4["Omega_rad_s"]) < 5e-3                   # two-digit rounding of one turn per solar day
    assert rel(2 * ch13.OMEGA_SOLAR_DAY, s4["f_at_poles_abs_s"]) < 5e-3
    assert rel(ch13.EARTH_RADIUS_MEAN / 1e3, s4["earth_radius_km"]) < 1e-12
    assert rel(ch13.OMEGA_EARTH / ch13.OMEGA_SOLAR_DAY - 1, s4["derived_by_analyst"]["relative_difference"]) < 5e-3
    ro = b["sec_13_5_geostrophic"]["rossby_number_example"]
    assert rel(ch13.rossby_number(ro["U_m_s"], ro["f_s"], ro["L_km"] * 1e3), ro["Ro"]) < 1e-12
    s2, s3 = b["sec_13_2_stratification"], b["sec_13_3_equations"]
    ad = ch13.ocean_density_gradient_budget(-1.0, 1025.0, s2["ocean_sound_speed_m_s"], g=9.81)["adiabatic"]
    assert rel(-ad, s2["derived_by_analyst"]["g_rho_over_c2_with_rho_1025"]) < 1e-3
    assert rel(-ad, s2["ocean_adiabatic_density_gradient_kg_m4"]) < 0.10         # the page rounds to one digit
    assert rel(ch13.scale_height(s2["ocean_sound_speed_m_s"], g=9.81) / 1e3, s3["derived_by_analyst"]["c2_over_g_km_c1520"]) < 1e-3
    assert rel(ch13.scale_height(s2["ocean_sound_speed_m_s"], g=9.81) / 1e3, s3["scale_height_ocean_km"]) < 0.20   # "about": order of magnitude
    assert rel(ch13.scale_height(340.0, g=9.81) / 1e3, s3["scale_height_atmosphere_km"]) < 0.20
    Ga = float(STRAT.adiabatic_lapse_rate())
    assert rel(-Ga * 1e3, s2["dry_adiabatic_decrease_C_per_km_approx"]) < 0.03   # the round value of the page
    assert ch13.lapse_rate_table(-s2["troposphere_lapse_C_per_km"] * 1e-3, Gamma_a=Ga)["verdict"].tolist() == ["stable", "stable"]
    assert rel(2 * PI / s2["N_max_thermocline_s"] / 60.0, s2["N_max_period_min"]) < 0.06


@book_only
def test_book_V6_ekman_layers_and_exercises():  # V6 §13.6–§13.7, exercises (slip #10 recorded)
    b = book()
    s7 = b["sec_13_7_ekman_bottom"]
    d_formula = ch13.ekman_depth(b["sec_13_3_equations"]["nu_air_m2_s"], s7["f_s"])
    assert rel(d_formula, s7["derived_by_analyst"]["delta_from_formula_nu_1p5e-5_f_1e-4_m"]) < 1e-3
    assert rel(s7["laminar_air_delta_m"] / d_formula, s7["derived_by_analyst"]["printed_over_formula"]) < 5e-3
    assert rel(d_formula, s7["laminar_air_delta_m"]) > 0.2                       # slip #10: the printed thickness is not the formula's
    assert rel(ch13.eddy_viscosity_from_depth(s7["observed_boundary_layer_km"] * 1e3, s7["f_s"]), s7["implied_eddy_viscosity_m2_s"]) < 1e-12
    ex = b["exercises"]
    f45 = float(ch13.coriolis_parameter(np.deg2rad(45.0)))
    slope = f45 * ex["13_1"]["U_m_s"] / 9.81                                     # geostrophy: g d(eta)/dy = f U
    u, _ = ch13.geostrophic_velocity(0.0, -1000.0 * 9.81 * slope, f45, 1000.0)
    assert rel(u, ex["13_1"]["U_m_s"]) < 1e-12 and rel(slope * 1e5, ex["13_1"]["answer_slope_cm_per_km"]) < 5e-3
    Mx, My = ch13.ekman_bottom_transport(ex["13_3"]["U_m_s"], 0.0, ex["13_3"]["nu_v_m2_s"], f45)
    assert rel(My, ex["13_3"]["answer_transport_m2_s"]) < 1e-3
    f_tank = 2 * (2 * PI * ex["13_2"]["rpm"] / 60.0)
    assert rel(ch13.ekman_depth(ex["13_2"]["nu_m2_s"], f_tank), ex["13_2"]["derived_by_analyst_delta_m_with_f_2Omega"]) < 2e-3
    e5 = ex["13_5"]
    fS = float(ch13.coriolis_parameter(np.deg2rad(-e5["latitude_deg_S"])))
    c5 = math.sqrt(9.81 * e5["density_jump_kg_m3"] / 1025.0 * e5["thermocline_depth_m"])
    assert rel(c5, e5["derived_by_analyst"]["c_m_s"]) < 5e-3 and rel(ch13.rossby_radius(c5, fS) / 1e3, e5["derived_by_analyst"]["decay_scale_km"]) < 5e-3
    side = ch13.kelvin_decay_side(fS, -1)                                        # southward along an eastern-boundary coast
    assert side["coast_on"] == "left" and "left" in e5["derived_by_analyst"]["direction"]
    e4 = ex["13_4"]["derived_by_analyst_12p42h"]
    f45, H4 = f45, ex["13_4"]["depth_km"] * 1e3
    om = 2 * PI / (12.42 * 3600.0)
    c4 = math.sqrt(9.81 * H4)
    K = math.sqrt(om ** 2 - f45 ** 2) / c4
    assert rel(ch13.poincare_omega(K, f45, c4), om) < 1e-12 and rel(2 * PI / K / 1e3, e4["wavelength_km"]) < 2e-3
    assert rel(om / f45, e4["axis_ratio"]) < 2e-3 and rel(om / K, e4["c_m_s"]) < 2e-3
    assert rel(ch13.poincare_group_velocity(K, 0.0, f45, c4)[0], e4["cg_m_s"]) < 2e-3


@book_only
def test_book_V6_modes_waves_and_regimes():  # V6 §13.9–§13.12 (slips #6, #11 recorded)
    b = book()
    s9 = b["sec_13_9_normal_modes"]
    N, H = s9["N_typical_s"], s9["H_km"] * 1e3
    md = ch13.modes_uniform_N(N, H, g=9.81, n_modes=4)
    assert np.max(np.abs(md.c / np.array(s9["derived_by_analyst"]["exact_roots_g9p81_N1e-3_H5000"]) - 1)) < 1e-6
    assert rel(N * H / md.c[0], s9["derived_by_analyst"]["NH_over_c0"]) < 1e-3 and N * H / md.c[0] < 0.1      # slip #6: not 1
    assert rel(ch13.baroclinic_mode_speed(N, H, 1), s9["derived_by_analyst"]["c1_NH_over_pi"]) < 1e-4
    assert rel(ch13.baroclinic_mode_speed(N, H, 1), s9["c1_m_s"]) < 0.25         # "about": the page rounds 1.6 up
    assert rel(ch13.equivalent_depth(s9["c1_m_s"], g=9.81), s9["He_m"]) < 0.03   # consistent with its own rounded speed
    s10 = b["sec_13_10_regimes"]
    hi = s10["high_frequency_example"]
    assert rel(hi["beta"] / (hi["omega_s"] * 2 * PI / (hi["wavelength_km"] * 1e3)), s10["derived_by_analyst"]["ratio_high"]) < 1e-3
    lo = s10["low_frequency_example"]
    kk = 2 * PI / (lo["wavelength_km"] * 1e3)
    for cc, key in ((lo["c_barotropic_m_s"], "ratio_barotropic"), (lo["c_baroclinic_m_s"], "ratio_baroclinic")):
        ts = ch13.dispersion_term_sizes(kk, 0.0, cc, 1e-4, lo["beta"], lo["omega_s"])
        assert rel(abs(ts["omega3"] / ts["beta"]), lo[key]) < 0.01
    s11, s12 = b["sec_13_11_poincare"]["inertial_example"], b["sec_13_12_kelvin"]
    assert rel(ch13.inertial_radius(s11["q_m_s"], 1e-4) / 1e3, s11["radius_km"]) < 1e-12
    dk = s12["deep_sea_example"]
    c = ch13.long_wave_speed(dk["H_km"] * 1e3, g=9.81)
    assert rel(c, dk["c_m_s"]) < 0.01 and rel(ch13.rossby_radius(c, dk["f_s"]) / 1e3, dk["Lambda_km"]) < 0.01   # two-digit rounding
    Li = ch13.rossby_radius_internal(N, H, dk["f_s"]) / 1e3
    assert rel(Li, s12["derived_by_analyst"]["NH_over_pi_f_with_sec_13_9_values_km"]) < 5e-3
    assert s12["internal_rossby_radius_ocean_km"] / Li > 2.5                     # slip #11: inconsistent inputs, not the formula


@book_only
def test_book_V6_rossby_eady_and_turbulence():  # V6 §13.15–§13.18 (slip #12 recorded)
    b = book()
    s15 = b["sec_13_15_rossby"]
    inp = s15["inputs"]
    for cc, want in zip((inp["c_barotropic_m_s"], inp["c_baroclinic_m_s"]), s15["derived_by_analyst"]["omega_max_over_f0"]):
        mx = ch13.rossby_max_frequency(inp["beta"], inp["f0_s"], cc)
        assert rel(mx["omega_max"] / inp["f0_s"], want) < 1e-9
        assert rel(mx["omega_max"] * inp["f0_s"] / (inp["beta"] * cc), s15["omega_max_f0_over_beta_c"]) < 1e-12
        assert rel(mx["k"] * cc / inp["f0_s"], s15["kc_over_f0_at_max"]) < 1e-12
    assert rel(-ch13.rossby_long_wave_speed(inp["beta"], inp["f0_s"], inp["c_baroclinic_m_s"]), s15["derived_by_analyst"]["cx_long"]) < 1e-9
    shortest = 2 * PI / ch13.rossby_max_frequency(inp["beta"], inp["f0_s"], inp["c_baroclinic_m_s"])["omega_max"]
    assert rel(shortest, 365.25 * 86400.0) < 0.01                               # a year to 0.5 % with the page's round inputs
    e7 = b["exercises"]["13_7"]
    f45 = float(ch13.coriolis_parameter(np.deg2rad(45.0)))
    Nn = math.sqrt(9.81 / 300.0 * e7["potential_temperature_rise_C"] / (e7["height_km"] * 1e3))
    assert rel(Nn, e7["answer_N_s"]) < 5e-4
    c1 = ch13.baroclinic_mode_speed(Nn, e7["height_km"] * 1e3, 1)
    assert rel(c1, e7["answer_c1_m_s"]) < 5e-4
    assert rel(-ch13.rossby_long_wave_speed(2e-11, f45, c1), e7["answer_cx_m_s"]) < 5e-3            # needs the round beta
    cx_lat = -ch13.rossby_long_wave_speed(float(ch13.beta_parameter(np.deg2rad(45.0))), f45, c1)
    assert rel(cx_lat, e7["derived_by_analyst"]["cx_with_beta_at_45deg"]) < 2e-3 and rel(cx_lat, e7["answer_cx_m_s"]) > 0.15   # slip #12
    s17 = b["sec_13_17_baroclinic"]
    assert rel(ch13.eady_critical(), s17["alpha_c_H"]) < 5e-3 and rel(ch13.eady_critical(), s17["derived_by_analyst"]["alpha_c_H"]) < 1e-9
    wl = ch13.eady_wavelengths()
    assert rel(wl["cutoff_over_Lambda"], s17["cutoff_wavelength_over_Lambda"]) < 0.01 and rel(wl["fastest_over_Lambda"], s17["fastest_wavelength_over_Lambda"]) < 5e-3
    assert rel(ch13.eady_fastest()["sigma_nd"], s17["derived_by_analyst"]["max_growth_times_NH_over_fU0"]) < 1e-8
    fj = b["sec_13_18_geostrophic_turbulence"]["example"]
    r = ch13.fjortoft_transfer(1.0, fj["K1_over_K0"], fj["K2_over_K0"])
    assert rel(r["energy_ratio"], fj["S1_over_S2"]) < 1e-12 and rel(r["enstrophy_ratio"], fj["enstrophy_ratio"]) < 1e-12


def test_public_V7_no_private_book_number_is_typed_in_this_file():  # CLAUDE.md rule 9
    if not BOOK.exists():
        pytest.skip("private list absent")
    b = book()
    src = Path(__file__).read_text(encoding="utf-8")
    for lit in b.get("_forbidden_public", []):
        assert lit not in src, lit
    for ent in b.get("_forbidden_public_regex", []):
        pat = ent if isinstance(ent, str) else ent["regex"]
        assert re.search(pat, src) is None, pat
