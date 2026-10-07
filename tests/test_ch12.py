"""Verification suite for Chapter 12 — Turbulence (Kundu, Cohen & Dowling 5e, §§12.1–12.13).

Evidence levels (``verify-implementation`` skill): V1 analytic · V2 symbolic (sympy / dimensions) · V3 convergence ·
V4 conservation / invariant / integral identity · V5 published benchmark (``reference/ch12/``, see SOURCES.md) · V6 book
value (private, ``tests/book_values_ch12.json``, skipped when absent) · V7 limits / symmetry / invariance.  ``stat`` marks a
sampled check: fixed seed, tolerance 5 standard errors of the estimator.
Every test name is ``test_<concept>_<level>_<what>``; the comment on the ``def`` line repeats the level and the curation ID.

A items (CORE, ≥ 2 independent levels): C01 ensemble average and its rules (12.1)–(12.9) · C02 lag correlation, integral
scale (12.17)–(12.18) · C03 spectrum pair (12.20)–(12.22) · C04 Reynolds-averaged momentum (12.30) · C05 isotropic
dissipation (12.43) · C06 energy budgets (12.46)–(12.47) · C07 Kolmogorov scales (12.50) · C08 −5/3 law (12.54) · C09 plane
jet (12.62)–(12.68) · C10 law of the wall (12.80)–(12.82) · C11 logarithmic law (12.88), (12.93) · C12 eddy viscosity and
mixing length (12.94)–(12.101) · C13 k–ε (12.103)–(12.105) · C14 Richardson numbers (12.107)–(12.109) · C15 Monin–Obukhov
(12.110)–(12.111) · C16 Taylor dispersion (12.119)–(12.129).

Derivations: every ★★ / ★★★ D row (D02 D04 D05 D06 D07 D08 D09 D10 D12 D13 D14 D15 D16 D17 D19 D21 D23 D26) is re-derived in a
``test_*_V2_derivation`` test **independently of the chapter's own sympy engines** (which are tested separately).  D06, D09,
D10 and the temperature-variance budget use an explicit finite ensemble: polynomial solenoidal fields with random signs,
averaged exactly over the four sign combinations, so that "average" is a finite sum and each budget is an identity
``mean(multiplier × Navier–Stokes residual) = left − right of the book's equation`` that sympy expands to zero.

Printed slips kept as wrong variants that must FAIL: #1 ``inertial_spectrum_1d(printed=True)`` · #3 the scalar flux "at
y = 0" · #5 ``general_similarity_check`` on the exponential family · #6 ``rans_eddy_viscosity_residual(printed=True)`` · #12
``eddy_diffusivity_asymptote(printed=True)`` and ``dispersion_regime`` · #13 ``smoke_plume_width`` slopes · #15
``temperature_variance_sympy()["printed_check"]``.

Run: ``.venv/Scripts/python.exe -m pytest tests/test_ch12.py -q -p no:cacheprovider``  (``-m "not slow"`` skips the script runs).
"""
from __future__ import annotations

import inspect
import itertools
import json
import math
import re
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest
import sympy as sp
from scipy import stats as sstats
from scipy.integrate import quad, solve_ivp
from scipy.signal import welch
from scipy.special import erf

from fluidpy import ch12_turbulence as ch12
from fluidpy.core import boundary_layer as BL
from fluidpy.core import dimensional as DIM
from fluidpy.core import jets as JET
from fluidpy.core import laminar as LAM
from fluidpy.core import stratification as STRAT
from fluidpy.core import turbstats as TS
from fluidpy.core import wall_turbulence as WT
from fluidpy.core.units import Q_, dimensional_check
from tools.convergence import observed_order

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference" / "ch12"
BOOK = Path(__file__).resolve().parent / "book_values_ch12.json"
book_only = pytest.mark.skipif(not BOOK.exists(), reason="book values are private; see CLAUDE.md rule 9")
needs_dns = pytest.mark.skipif(not (REF / "lee_moser_2015_channel_mean.csv").exists(),
                               reason="run reference/ch12/make_refs.py")
slow = pytest.mark.slow
G0 = ch12.G0
PI = math.pi
KAP, BLOG = 0.41, 5.0          # the labelled "classical" textbook pair (WT.LOG_LAW_CONSTANTS["classical"])


def book():
    return json.loads(BOOK.read_text(encoding="utf-8"))


def rel(a, b):
    return abs(float(a) - float(b)) / max(abs(float(b)), 1e-300)


def z0(expr) -> bool:
    return sp.simplify(sp.expand(expr)) == 0


def dns(case: int) -> dict:
    """One Lee & Moser (2015) mean profile from reference/ch12 (subset; see SOURCES.md)."""
    rows = [ln.split(",") for ln in (REF / "lee_moser_2015_channel_mean.csv").read_text(encoding="utf-8").splitlines()
            if ln and not ln.startswith("#") and not ln.startswith("case")]
    a = np.array([[float(v) for v in r] for r in rows if int(r[0]) == case])
    return dict(Re_tau=float(a[0, 1]), y_over_delta=a[:, 2], yplus=a[:, 3], Uplus=a[:, 4], dUdy_plus=a[:, 5])


# ======================================================================================================================
# Contract: every name of design Part C exists, is re-exported by ch12 and is exercised in this file
# ======================================================================================================================
CONTRACT = (
    "FREE_SHEAR_CONSTANTS", "K_EPSILON_CONSTANTS", "LOG_LAW_CONSTANTS", "anisotropy_tensor", "autocorrelation",
    "batchelor_scale", "book_slips", "boundary_layer_stress_from_profile", "buoyant_production", "cascade_tiers",
    "central_moment", "channel_energy_budget", "channel_energy_budget_at", "channel_mixing_length",
    "channel_momentum_residual", "channel_pressure_gradient", "channel_total_stress", "check_averaging_rules",
    "closure_count", "coles_wake", "composite_profile", "composite_profile_plus", "convective_eddy_diffusivity",
    "convective_velocity_scale", "conventions", "correlated_pair", "correlation", "correlation_coefficient",
    "correlation_from_spectrum", "correlation_spectrum_pair", "correlation_time", "cross_correlation",
    "defect_law_groups", "diffusivity_from_variance", "dimensionless_shear", "dispersion_local_slope",
    "dispersion_rate_from_particles", "dispersion_regime", "displaced_parcel_uv", "dissipation_isotropic",
    "dissipation_outer_scaling", "dissipation_rate", "divergence_rms", "dns_grid_points", "drag_coefficient_neutral",
    "eddy_diffusivity_estimate", "eddy_diffusivity_exponential", "eddy_diffusivity_taylor",
    "eddy_viscosity_from_data", "eddy_viscosity_stress", "effective_samples", "ensemble_average",
    "fit_inertial_range", "fit_log_law", "flux_from_gradient_richardson", "flux_richardson",
    "flux_richardson_surface_layer", "free_shear_centerline", "free_shear_exponents", "free_shear_flow",
    "free_shear_profile", "frequency_to_wavenumber_spectrum", "friction_law_from_overlap", "friction_reynolds_number",
    "friction_velocity", "friction_velocity_from_wind", "from_wall_units", "frozen_field_probe", "gaussian_profile",
    "general_similarity_check", "gradient_diffusion_flux", "gradient_moments_isotropic_sympy",
    "gradient_richardson_surface_layer", "gradient_richardson_thermal", "inertial_range_decades",
    "inertial_spectrum_1d", "inertial_spectrum_3d", "integral_scale", "integral_scale_from_spectrum",
    "isotropic_correlation_tensor", "isotropic_scales", "isotropic_tensor_divergence_sympy", "isotropy_report",
    "jet_momentum_flux_per_span", "jet_tke_budget", "k_epsilon_decay", "k_epsilon_eddy_viscosity",
    "k_epsilon_length_scale", "k_epsilon_loglayer_kappa", "k_epsilon_rhs", "kolmogorov_constants",
    "kolmogorov_normalize_spectrum", "kolmogorov_scales", "langevin_particles", "law_of_the_wall_groups",
    "layer_name", "local_reynolds_number_exponent", "log_law", "log_law_crossing", "log_law_defect",
    "log_law_indicator", "longitudinal_transverse_correlation", "make_ensemble", "mass_fraction_from_volume_fraction",
    "mean_divergence", "mean_energy_budget", "mean_energy_budget_sympy", "mean_heat_flux", "mean_scalar_flux",
    "mean_stress_tensor", "mean_to_turbulent_dissipation_ratio", "mixing_length_eddy_viscosity",
    "mixing_length_intercept", "mixing_length_stress", "mixing_length_wall_profile", "mixing_potential_energy_change",
    "mixture_density", "model_spectrum", "moment", "monin_obukhov_from_fluxes", "monin_obukhov_length",
    "nagib_chauhan_B", "nagib_chauhan_kappa", "one_dimensional_from_3d", "one_equation_closure",
    "overlap_matching_sympy", "parcel_uv_expected", "periodogram", "pipe_bulk_velocity_loglaw",
    "pipe_friction_factor_turbulent", "pipe_pressure_gradient", "plane_jet_cross_velocity",
    "plane_jet_eddy_viscosity_profile", "plane_jet_entrainment_velocity", "plane_jet_mass_fraction",
    "plane_jet_mean_velocity", "plane_jet_reynolds_stress", "plane_jet_similarity_sympy", "plane_jet_stress_profile",
    "plane_jet_volume_flux", "plume_concentration", "product_average_split", "profile_integrals", "random_walk",
    "rans_2d_residual", "rans_eddy_viscosity_residual", "rans_momentum_residual", "rans_sympy", "reynolds_decompose",
    "reynolds_stress", "reynolds_stress_budget_sympy", "reynolds_stress_production", "richardson_diffusivity",
    "rough_wall_log_law", "round_jet_distance_for_mass_fraction", "scalar_flux_per_span", "scalar_spectrum",
    "scale_ordering", "scale_separation", "scale_table", "shear_flow_eddy_viscosity_solve", "shear_production",
    "shell_spectrum", "skin_friction_zpg", "slot_mass_flux", "slot_momentum_flux", "smoke_plume_width",
    "smooth_signal", "spalding_slope", "spalding_uplus", "spalding_yplus", "spatial_correlation",
    "spectrum_from_correlation", "spectrum_variance", "standard_error", "standard_error_of_mean", "statistics",
    "stoichiometric_mass_fraction", "stratified_tke_budget", "stress_partition", "surface_layer_regime",
    "surface_layer_state", "surface_layer_wind", "synthetic_solenoidal_field", "taylor_dispersion",
    "taylor_dispersion_exponential", "taylor_dispersion_gaussian", "taylor_dispersion_rate", "taylor_frozen",
    "taylor_microscale", "taylor_reynolds_number", "temperature_variance_budget", "temperature_variance_sympy",
    "thin_shear_layer_terms", "time_average", "time_average_exp_cos", "tke_budget", "tke_budget_sympy",
    "total_stress", "transverse_from_longitudinal", "turbulence_regime", "turbulent_heat_flux",
    "turbulent_kinetic_energy", "turbulent_prandtl", "velocity_covariance", "velocity_defect", "virtual_origin_fit",
    "viscous_length", "viscous_sublayer", "volume_average", "wall_stress_from_pressure_gradient", "wall_units",
    "white_noise_field", "wrong_exponent_fluxes", "zpg_boundary_layer",
)  # the 215 names of analysis/ch12_design.md Part C (C.1–C.3), extracted 2026-10-07
PARITY = (
    'ch12.time_average_exp_cos(5.0, 1.0, 0.5, 10.0, 2*np.pi, 0.4)[2]',
    'ch12.time_average_exp_cos(5.0, 1.0, 0.5, 10.0, 2*np.pi, 20.0)[1]',
    'ch12.time_average_exp_cos(5.0, 1.0, 0.5, 10.0, 2*np.pi, 0.4)[0]',
    'ch12.standard_error_of_mean(0.3, 8)',
    'ch12.correlation_spectrum_pair ("exponential", 1.0, 0.5)["S"](2.0)',
    'ch12.correlation_spectrum_pair("gaussian", 1.0, 0.5)["Lambda_t"]',
    'ch12.correlation_spectrum_pair("damped_cosine", 1.0, 0.5, 6.0)["S"](6.0)',
    'ch12.integral_scale_from_spectrum(0.15915494309, 1.0)',
    'ch12.parcel_uv_expected(2.0, 0.1, 0.5, 1.0)',
    'ch12.parcel_uv_expected(-3.0, 0.2, 0.4, 0.3)',
    'ch12.eddy_viscosity_from_data(-0.1, 2.0)',
    'ch12.kolmogorov_scales(1e-6, 1.0)[0]',
    'ch12.kolmogorov_scales(1.5e-5, 0.125)[1]',
    'ch12.scale_separation(1e7)["eta_over_L"]',
    'ch12.inertial_spectrum_1d(100.0, 1.0)',
    'ch12.model_spectrum(50.0, 1.0, 1e-6, L=1.0, kind="pope")',
    'ch12.channel_energy_budget_at(12.0, 1000.0, 0.41, 26.0)["production"]',
    'ch12.channel_energy_budget_at(100.0, 1000.0, 0.41, 26.0)["mean_dissipation"]',
    'ch12.channel_energy_budget_at(30.0, 1000.0, 0.41, 26.0)["Uplus"]',
    'ch12.plane_jet_mean_velocity(4.0, 0.0, 1.0, 1.0, C5="from_invariant", xi_half=0.1)',
    'ch12.plane_jet_volume_flux(1.0, 1.0, 1.0, C5="from_invariant", xi_half=0.1)',
    'ch12.wrong_exponent_fluxes(4.0, -0.4, 1.0)[0]',
    'ch12.free_shear_exponents("round_plume")["velocity"]*1.0',
    'ch12.log_law(100.0, kappa=0.41, B=5.0)',
    'ch12.spalding_uplus(12.0, kappa=0.41, B=5.0)',
    'ch12.log_law_crossing(kappa=0.41, B=5.0)',
    'ch12.stress_partition(30.0, 1000.0, kappa=0.41, B=5.0)["viscous"]',
    'ch12.layer_name(12.0, 0.012)',
    'ch12.mixing_length_wall_profile(10.0, 0.41) ["slope"]',
    'ch12.mixing_length_wall_profile (100.0, 0.41, damping="van_driest", A_plus=26.0)["Uplus"]',
    'ch12.mixing_length_intercept(0.41)',
    'ch12.mixing_length_intercept(0.41, A_plus=26.0)',
    'ch12.monin_obukhov_from_fluxes(0.108, -30.0, 1.2, 1005.0, 300.0, kappa=0.4)',
    'ch12.surface_layer_wind(10.0, 0.3, 0.03, 82.98165, kappa=0.4)',
    'ch12.surface_layer_wind(10.0, 0.3, 0.03, -16.6, kappa=0.4, unstable="businger_dyer")',
    'ch12.flux_richardson_surface_layer(10.0, 82.98165)',
    'ch12.turbulence_regime(0.12)',
    'ch12.gradient_richardson_thermal(0.010, 0.1, 1/300, Gamma_a=-0.0098)["verdict_kundu"]',
    'ch12.taylor_dispersion_exponential(10.0, 1.0, 10.0)',
    'ch12.taylor_dispersion_exponential(0.01, 1.0, 10.0)',
    'ch12.taylor_dispersion_gaussian(10.0, 1.0, 10.0)',
    'ch12.eddy_diffusivity_exponential(10.0, 1.0, 10.0)',
    'ch12.dispersion_local_slope(10.0, 10.0)',
    'ch12.smoke_plume_width(1000.0, 5.0, 0.5, 20.0)',
    'ch12.dispersion_regime(3.0, 10.0)',
    'ch12.k_epsilon_decay(1.0, 1.0, 10.0, C_eps2=1.92)[0]',
    'ch12.k_epsilon_loglayer_kappa(0.09, 1.44, 1.92, 1.3)',
)  # the `py:` expressions of the explainer storyboards (design Part B), verbatim
# END-OF-CONTRACT-LIST


def test_contract_V1_every_part_c_name_exists_and_is_reexported():  # V1 · design Part C / convention 1
    missing = [n for n in CONTRACT if not hasattr(ch12, n)]
    assert missing == []
    for mod in (TS, WT):
        public = [n for n, o in vars(mod).items() if inspect.isfunction(o) and o.__module__ == mod.__name__
                  and not n.startswith("_")]
        assert all(getattr(ch12, n) is getattr(mod, n) for n in public), mod.__name__
        assert len(public) >= 35


def test_contract_V1_every_part_c_name_is_called_by_a_test():  # self-audit of this file (Part C: ≥ a smoke test each)
    src = Path(__file__).read_text(encoding="utf-8")
    body = src.split("# END-OF-CONTRACT-LIST", 1)[1]
    not_used = [n for n in CONTRACT if not re.search(r"(?<![A-Za-z_0-9])" + re.escape(n) + r"(?![A-Za-z_0-9])", body)]
    assert not_used == [], not_used


def test_contract_V1_docstrings_cite_book_and_validation():  # CLAUDE.md rule 5
    bad = []
    for n in CONTRACT:
        o = getattr(ch12, n)
        if not inspect.isfunction(o):
            continue
        doc = inspect.getdoc(o) or ""
        if "Book:" not in doc or "Validation" not in doc or "Assum" not in doc:
            bad.append(n)
    assert bad == [], bad


def test_parity_V1_every_design_parity_expression_evaluates():  # design Part B: the `py:` rows of E1–E10 and B1
    for expr in PARITY:
        v = eval(expr, {"ch12": ch12, "np": np, "math": math})  # noqa: S307 - our own fixed strings
        assert isinstance(v, str) or math.isfinite(float(v)), expr
    scope = {"ch12": ch12, "np": np}
    # pinned values (closed forms worked by hand here, not copied from the functions)
    assert rel(eval(PARITY[3], scope), 0.3 / math.sqrt(8)) < 1e-14
    assert rel(ch12.parcel_uv_expected(2.0, 0.1, 0.5, 1.0), -0.1) < 1e-14
    assert rel(ch12.kolmogorov_scales(1e-6, 1.0)[0], 1e-6 ** 0.75) < 1e-14
    assert rel(ch12.scale_separation(1e7)["eta_over_L"], 1e7 ** -0.75) < 1e-14
    assert rel(ch12.log_law(100.0, kappa=0.41, B=5.0), math.log(100.0) / 0.41 + 5.0) < 1e-14
    assert rel(ch12.taylor_dispersion_exponential(10.0, 1.0, 10.0), 200.0 * math.exp(-1.0)) < 1e-13
    assert rel(ch12.eddy_diffusivity_exponential(10.0, 1.0, 10.0), 10.0 * (1 - math.exp(-1.0))) < 1e-14
    assert rel(ch12.dispersion_local_slope(10.0, 10.0), (1 - math.exp(-1.0)) / math.exp(-1.0)) < 1e-13
    assert rel(ch12.k_epsilon_loglayer_kappa(0.09, 1.44, 1.92, 1.3), math.sqrt(0.3 * 0.48 * 1.3)) < 1e-14
    assert ch12.layer_name(12.0, 0.012) == "buffer layer" and ch12.turbulence_regime(0.12) == "shear-driven"
    assert ch12.dispersion_regime(3.0, 10.0) == "ballistic"      # the E10 row that sits exactly on t = 0.3 Λ_t


# ======================================================================================================================
# Book slips (design convention 8): each planted wrong variant fails its check
# ======================================================================================================================
def test_book_slips_V1_table_lists_fifteen_and_names_existing_switches():  # N-rows of §9
    slips = ch12.book_slips()
    assert [d["id"] for d in slips] == list(range(1, 16))
    assert {"id", "where", "printed", "corrected", "how_to_tell", "coded_in", "test"} <= set(slips[0])
    for fn in ("inertial_spectrum_1d", "rans_eddy_viscosity_residual", "eddy_diffusivity_asymptote"):
        assert "printed" in inspect.signature(getattr(ch12, fn)).parameters
    conv = ch12.conventions()
    assert {"symbol", "meaning", "code_name", "unit", "note"} <= set(conv.columns) and len(conv) >= 10
    assert hasattr(ch12, "ensemble_average") and "ensemble_average" in set(conv["code_name"])


def test_slip01_inertial_spectrum_V2_printed_exponent_fails_units_and_slope():  # V2 · C08 · slip #1
    eps, k1 = Q_(0.5, "m**2/s**3"), Q_(40.0, "1/m")
    good = eps ** Fraction(2, 3) * k1 ** Fraction(-5, 3)
    bad = eps ** Fraction(2, 3) * k1 ** Fraction(5, 3)
    assert good.check("[length]**3/[time]**2") and not bad.check("[length]**3/[time]**2")
    k = np.geomspace(1.0, 1e3, 40)
    slope = lambda S: np.polyfit(np.log(k), np.log(S), 1)[0]  # noqa: E731
    assert abs(slope(ch12.inertial_spectrum_1d(k, 0.5)) + 5 / 3) < 1e-12
    assert abs(slope(ch12.inertial_spectrum_1d(k, 0.5, printed=True)) - 5 / 3) < 1e-12   # the page: rises with k
    assert ch12.inertial_spectrum_1d(100.0, 0.5, printed=True) > 1e5 * ch12.inertial_spectrum_1d(100.0, 0.5)


def test_slip03_scalar_flux_V4_integral_at_fixed_x_is_the_invariant_not_a_value_at_y0():  # V4 · C09 (N118) · slip #3
    rho, Js, Ms = 1.2, 3.0, 0.4
    kw = dict(C5="from_invariant", xi_half=0.1)
    flux, at_y0 = [], []
    for x in (0.5, 1.0, 2.0, 4.0):
        y = np.linspace(-1.2 * x, 1.2 * x, 4001)
        U = ch12.plane_jet_mean_velocity(x, y, Js, rho, **kw)
        Y = ch12.plane_jet_mass_fraction(x, y, Ms, Js, rho, C6="from_invariant", xi_half_Y=0.15, **kw)
        flux.append(ch12.scalar_flux_per_span(y, U, Y, rho))
        at_y0.append(rho * float(np.interp(0.0, y, U * Y)))     # what "[…]_{y=0}" would be: one point of the integrand
    assert np.max(np.abs(np.array(flux) / Ms - 1)) < 1e-9      # (12.70): ρ∫ȲU dy = Ṁ_s at every x
    assert rel(at_y0[0] / at_y0[-1], 8.0) < 1e-9               # the y = 0 reading decays as 1/x — it is not conserved


def test_slip05_general_similarity_V2_exponential_family_fails_power_family_passes():  # V2 · C09 (N125) · slip #5
    x, a, m, n = sp.symbols("x a m n", positive=True)

    def coeffs(delta, Ucl, Psi):
        c1 = delta * sp.diff(Ucl, x) / Ucl
        return sp.simplify(c1), sp.simplify(c1 + sp.diff(delta, x)), sp.simplify(Psi / Ucl ** 2)

    c1, c2, c3 = coeffs(x ** m, x ** n, x ** (2 * n + m - 1))
    assert z0(c1 / c3 - n) and z0(c2 / c3 - (m + n))            # proportional: one common factor x^{m−1}
    e1, e2, e3 = coeffs(sp.exp(a * x), sp.exp(-a * x), sp.exp(-a * x))
    assert e2 == 0 and e1 != 0 and e3 != 0                      # the page's family: middle coefficient ≡ 0
    assert sp.simplify(sp.diff(sp.exp(-a * x) ** 2 * sp.exp(a * x), x)) != 0   # and U_CL²δ is not constant
    g = ch12.general_similarity_check(lambda s: np.exp(0.7 * s), lambda s: np.exp(-0.7 * s), lambda s: np.exp(-0.7 * s), 2.0)
    assert abs(g["c2"]) < 1e-8 and abs(g["c1"]) > 1.0 and abs(g["momentum_flux_exponent"] + 1.4) < 1e-6
    for (mm, nn) in ((1.0, -0.5), (1.0, -0.4), (0.5, -0.25)):
        g = ch12.general_similarity_check(lambda s: s ** mm, lambda s: s ** nn, lambda s: s ** (2 * nn + mm - 1), 2.0)
        fac = 2.0 ** (mm - 1)
        assert max(abs(g["c1"] - nn * fac), abs(g["c2"] - (mm + nn) * fac), abs(g["c3"] - fac)) < 1e-7
        assert abs(g["momentum_flux_exponent"] - (2 * nn + mm)) < 1e-7    # flat only when m + 2n = 0


def _manufactured_eddy_closure():
    U = lambda p: np.array([np.sin(p[0]) * np.cos(p[1]), -np.cos(p[0]) * np.sin(p[1])])   # noqa: E731 solenoidal
    nuT = lambda p: 0.3 + 0.1 * np.sin(p[0] + 2 * p[1])                                     # noqa: E731
    e = lambda p: 0.5 + 0.2 * np.cos(p[0] - p[1])                                           # noqa: E731
    P = lambda p: 0.7 * p[0] ** 2 - 1.3 * p[0] * p[1] + 2.0 * p[1]                          # noqa: E731
    return U, P, nuT, e


def test_slip06_rans_eddy_viscosity_V1_corrected_matches_12_30_and_printed_differs():  # V1 · C12 (N165) · slip #6
    U, P, nuT, e = _manufactured_eddy_closure()
    nu, rho = 0.05, 1.3
    for pt in (np.array([0.3, 0.7]), np.array([1.1, -0.4])):
        r97 = ch12.rans_eddy_viscosity_residual(U, P, nuT, e, pt, nu=nu, rho=rho, h=1e-3)
        fields = {"U": lambda p, t: U(p), "P": lambda p, t: P(p),
                  "uu": lambda p, t: ch12.eddy_viscosity_stress(_grad(U, p), float(nuT(p)), float(e(p)))}
        r30 = ch12.rans_momentum_residual(fields, pt, nu=nu, rho0=rho, h=1e-3)
        assert np.max(np.abs(r97 - r30)) < 2e-5 * max(1.0, np.max(np.abs(r30)))     # (12.94) into (12.30) = (12.97)
        bad = ch12.rans_eddy_viscosity_residual(U, P, nuT, e, pt, nu=nu, rho=rho, h=1e-3, printed=True)
        gP = np.array([1.4 * pt[0] - 1.3 * pt[1], -1.3 * pt[0] + 2.0])
        assert np.max(np.abs((bad - r97) - (gP.sum() - gP) / rho)) < 1e-5            # the page's index: off by O(1)
        assert np.max(np.abs(bad - r97)) > 0.1


def _grad(U, p, h=1e-5):
    p = np.asarray(p, float)
    cols = []
    for j in range(p.size):
        dp = np.zeros_like(p)
        dp[j] = h
        cols.append((np.asarray(U(p + dp)) - np.asarray(U(p - dp))) / (2 * h))
    return np.stack(cols, axis=1)


def test_slip12_eddy_diffusivity_V7_long_time_constant_and_printed_condition_fails():  # V7 · C16 (N215) · slip #12
    u2, Lam = 0.36, 20.0
    exact = lambda t: ch12.eddy_diffusivity_exponential(t, u2, Lam)  # noqa: E731
    assert rel(ch12.eddy_diffusivity_asymptote(200.0, u2, Lam, which="long"), exact(200.0)) < 1e-4
    assert rel(ch12.eddy_diffusivity_asymptote(0.2, u2, Lam, which="short"), exact(0.2)) < 6e-3
    assert math.isnan(ch12.eddy_diffusivity_asymptote(0.2, u2, Lam, which="long"))          # condition t ≫ Λ_t not met
    assert math.isnan(ch12.eddy_diffusivity_asymptote(200.0, u2, Lam, which="short"))
    bad = ch12.eddy_diffusivity_asymptote(0.2, u2, Lam, which="long", printed=True)         # (12.129) as printed
    assert 90 < bad / exact(0.2) < 110
    t = np.array([0.01, 0.3, 1.0, 3.0, 100.0]) * Lam
    assert list(ch12.dispersion_regime(t, Lam)) == ["ballistic", "ballistic", "transition", "diffusive", "diffusive"]


def test_dispersion_regime_V7_inclusive_boundaries_and_both_sides():  # V7 · C16 · flagged item 3
    for Lam in (0.1, 1.0, 7.0, 10.0, 123.4):
        assert ch12.dispersion_regime(0.3 * Lam, Lam) == "ballistic"
        assert ch12.dispersion_regime(0.3 * Lam * (1 + 1e-9), Lam) == "transition"
        assert ch12.dispersion_regime(3.0 * Lam * (1 - 1e-9), Lam) == "transition"
        assert ch12.dispersion_regime(3.0 * Lam, Lam) == "diffusive"
    assert ch12.dispersion_regime(3.0, 10.0) == "ballistic" and ch12.dispersion_regime(30.0, 10.0) == "diffusive"
    assert ch12.dispersion_regime(1.0, 2.0, lo=0.5, hi=0.6) == "ballistic"
    # the regime names follow the local slope of the exact curve: 2 → 1
    assert ch12.dispersion_local_slope(0.3, 1.0) > 1.89 and ch12.dispersion_local_slope(3.0, 1.0) < 1.4


def test_slip13_smoke_plume_V1_linear_near_square_root_far():  # V1 · C16 (N212) · slip #13
    U, w, Lam = 5.0, 0.5, 20.0
    slope = lambda x: (math.log(ch12.smoke_plume_width(x * 1.001, U, w, Lam))  # noqa: E731
                       - math.log(ch12.smoke_plume_width(x / 1.001, U, w, Lam))) / (2 * math.log(1.001))
    assert abs(slope(0.001 * U * Lam) - 1.0) < 1e-3          # near the source: Z_rms ∝ x (12.121)
    assert abs(slope(1000.0 * U * Lam) - 0.5) < 1e-3         # far: ∝ √x (12.123) — the caption has them swapped
    x = 1000.0
    assert rel(ch12.smoke_plume_width(x, U, w, Lam) ** 2, ch12.taylor_dispersion_exponential(x / U, w * w, Lam)) < 1e-13
    assert rel(ch12.smoke_plume_width(1e-3, U, w, Lam), w * 1e-3 / U) < 1e-4
    z = np.linspace(-400.0, 400.0, 4001)
    for xx in (100.0, 1000.0):
        c = ch12.plume_concentration(xx, z, 2.0, U, ch12.smoke_plume_width(xx, U, w, Lam))
        assert rel(np.trapezoid(c, z), 2.0 / U) < 1e-9       # V4: ∫c dz = Q/U at every x


def test_slip15_temperature_variance_V2_engine_flags_the_printed_factor():  # V2 · C15 (N190) · slip #15
    r = ch12.temperature_variance_sympy()
    assert r["check"] is True and r["printed_check"] is False
    assert r["printed_minus_derived"] != 0
    assert "kappa" in {str(s) for s in r["printed_minus_derived"].free_symbols}


# ======================================================================================================================
# §12.3  C01 — ensemble averages, moments, the rules (12.1)–(12.11); Example 12.1
# ======================================================================================================================
def test_ensemble_average_V1_matches_hand_loops_and_moment_identities():  # V1 · C01 (12.1), (12.10), (12.11)
    rng = np.random.default_rng(3)
    s = rng.gamma(2.0, 1.5, size=(500, 7))
    acc1 = np.zeros(7)
    acc2 = np.zeros(7)
    for row in s:                      # the from-scratch loop of the notebook
        acc1 += row
        acc2 += row ** 2
    mean, var = acc1 / 500, acc2 / 500 - (acc1 / 500) ** 2
    assert np.allclose(TS.ensemble_average(s), mean, rtol=1e-13)
    assert np.allclose(TS.moment(s, 1), mean, rtol=1e-13) and np.allclose(TS.moment(s, 2), acc2 / 500, rtol=1e-13)
    assert np.allclose(TS.central_moment(s, 2), var, rtol=1e-11)
    m1, m2, m3 = (TS.moment(s, k) for k in (1, 2, 3))
    assert np.allclose(TS.central_moment(s, 3), m3 - 3 * m1 * m2 + 2 * m1 ** 3, rtol=1e-10)
    assert np.allclose(TS.moment(s, 0), 1.0)                                     # (12.5) with m = 0: mean(A) = A
    st = TS.statistics(s[:, 0])
    stn = TS.statistics(s[:, 0], normalized=True)
    assert rel(st["variance"], var[0]) < 1e-11 and rel(st["std"], math.sqrt(var[0])) < 1e-11
    assert rel(st["skewness"], TS.central_moment(s[:, 0], 3)) < 1e-12          # book: raw central moments
    assert rel(stn["skewness"], sstats.skew(s[:, 0])) < 1e-10
    assert rel(stn["kurtosis"], sstats.kurtosis(s[:, 0], fisher=False)) < 1e-10
    assert rel(st["rms"] ** 2, st["variance"] + st["mean"] ** 2) < 1e-12


def test_ensemble_average_stat_gaussian_moments_within_5_standard_errors():  # stat · C01 (12.1)
    N, sig = 200_000, 1.7
    u = np.random.default_rng(11).normal(0.0, sig, N)
    for m, exact, var_of_term in ((1, 0.0, sig ** 2), (2, sig ** 2, 2 * sig ** 4), (3, 0.0, 15 * sig ** 6),
                                  (4, 3 * sig ** 4, 96 * sig ** 8)):
        assert abs(TS.moment(u, m) - exact) < 5 * math.sqrt(var_of_term / N), m
    stn = TS.statistics(u, normalized=True)
    assert abs(stn["skewness"]) < 5 * math.sqrt(6 / N) and abs(stn["kurtosis"] - 3) < 5 * math.sqrt(24 / N)
    # volume average (12.3)
    f = np.random.default_rng(1).normal(size=(4, 5, 6))
    assert rel(TS.volume_average(f, m=2), np.mean(f ** 2)) < 1e-14
    x = np.linspace(0, 1, 11)
    w = np.gradient(x ** 2)
    assert rel(TS.volume_average(3 + 2 * x, weights=w), np.sum((3 + 2 * x) * w) / np.sum(w)) < 1e-14


def test_make_ensemble_stat_mean_variance_and_lag_correlation():  # stat · C01 (N07) / C02
    t = np.linspace(0.0, 4.0, 81)
    N, sig, tc = 20_000, 0.8, 0.5
    s = TS.make_ensemble(N, t, lambda tt: 2.0 + np.sin(tt), sig, tc, seed=5)
    assert s.shape == (N, 81)
    assert np.max(np.abs(TS.ensemble_average(s) - (2.0 + np.sin(t)))) < 5 * sig / math.sqrt(N)
    assert np.max(np.abs(TS.central_moment(s, 2) - sig ** 2)) < 5 * sig ** 2 * math.sqrt(2 / N)
    mean, fl = TS.reynolds_decompose(s)
    assert np.max(np.abs(fl.mean(axis=0))) < 1e-12 and np.allclose(mean + fl, s)        # (12.24)–(12.26)
    for lagk in (1, 5, 20):
        r = TS.correlation_coefficient(fl[:, 30], fl[:, 30 + lagk])
        rex = math.exp(-(t[lagk] - t[0]) / tc)
        assert abs(r - rex) < 5 * (1 - rex ** 2) / math.sqrt(N), lagk
    assert np.allclose(TS.make_ensemble(3, t, 1.5, 0.0, tc), 1.5)                      # V7 σ = 0
    assert np.array_equal(TS.make_ensemble(4, t, 0.0, 1.0, tc, seed=9), TS.make_ensemble(4, t, 0.0, 1.0, tc, seed=9))
    # standard error (N20): scatter of N-member means ∝ N^(−1/2)
    assert rel(TS.standard_error_of_mean(0.3, 8), 0.3 / math.sqrt(8)) < 1e-15
    assert rel(TS.standard_error(s[:400, 10]), np.std(s[:400, 10], ddof=1) / 20.0) < 1e-12
    Ns = np.array([25, 100, 400, 1600])
    scat = [np.std(s[: 12 * n, 40].reshape(12, n).mean(axis=1), ddof=1) for n in Ns]
    assert abs(np.polyfit(np.log(Ns), np.log(scat), 1)[0] + 0.5) < 0.25               # 12 groups: ±0.25 is ~2 s.e.


def test_averaging_rules_V1_linear_rules_round_off_product_rule_fails():  # V1 · C01 (12.4)–(12.9) · D01
    t = np.linspace(0.0, 2.0, 41)
    su = TS.make_ensemble(300, t, lambda tt: 1.0 + tt, 0.5, 0.3, seed=1)
    sv = 0.6 * su + TS.make_ensemble(300, t, lambda tt: np.cos(tt), 0.4, 0.2, seed=2)
    r = TS.check_averaging_rules(su, sv, t, A=2.0)
    linear = ("sum", "constant", "d/dt", "integral dt", "d/dx", "integral dx", "average of average")
    assert set(linear) | {"product"} <= set(r)
    assert max(r[k] for k in linear) < 1e-11
    mp, pm, cov = TS.product_average_split(su, sv)
    assert np.max(np.abs(mp - pm - cov)) < 1e-12
    assert rel(r["product"], np.max(np.abs(cov))) < 1e-10 and r["product"] > 0.05     # the one rule that fails
    assert np.allclose(cov[5], np.cov(su[:, 5], sv[:, 5], bias=True)[0, 1], rtol=1e-10)
    assert np.allclose(TS.correlation(*[x - x.mean(axis=0) for x in (su, sv)]), cov, atol=1e-12)


def test_averaging_rules_V2_derivation_D01_finite_ensemble_identities():  # V2 · D01 (★; cheap, kept)
    a1, a2, a3, b1, b2, b3, A = sp.symbols("a1 a2 a3 b1 b2 b3 A")
    t = sp.symbols("t")
    us = [a1 * sp.sin(t), a2 * t ** 2, a3 * sp.exp(t)]
    vs = [b1 * t, b2 * sp.cos(t), b3]
    avg = lambda xs: sum(xs) / len(xs)  # noqa: E731
    assert z0(avg([u + v for u, v in zip(us, vs)]) - avg(us) - avg(vs))                 # (12.4)
    assert z0(avg([A * u for u in us]) - A * avg(us))                                   # (12.5)
    assert z0(avg([sp.diff(u, t) for u in us]) - sp.diff(avg(us), t))                   # (12.6)
    assert z0(avg([sp.integrate(u, (t, 0, 1)) for u in us]) - sp.integrate(avg(us), (t, 0, 1)))   # (12.7)
    U, V = avg(us), avg(vs)
    cov = avg([(u - U) * (v - V) for u, v in zip(us, vs)])
    assert z0(avg([u * v for u, v in zip(us, vs)]) - U * V - cov)                       # mean(ũṽ) = ŪV̄ + mean(uv)
    assert not z0(cov)


def test_time_average_V2_derivation_example_12_1_window_integral():  # V2 · C01 (N23, a-D2) (12.2)
    t, s, A, B, tau, om, D = sp.symbols("t s A B tau omega Delta", positive=True)
    avg = sp.integrate(A * sp.exp(-s / tau) + B * sp.cos(om * s), (s, t - D / 2, t + D / 2)) / D
    x1, x2 = D / (2 * tau), om * D / 2
    form = sp.sinh(x1) / x1 * A * sp.exp(-t / tau) + sp.sin(x2) / x2 * B * sp.cos(om * t)
    assert sp.simplify((avg - form).rewrite(sp.exp)) == 0
    vals = dict(t=5.0, A=1.0, B=0.5, tau=10.0, omega=2 * PI, window=0.4)
    got = ch12.time_average_exp_cos(**vals)
    num = float(form.subs({t: 5.0, A: 1.0, B: 0.5, tau: 10.0, om: 2 * PI, D: 0.4}))
    assert rel(got[0], num) < 1e-13
    assert rel(got[1], math.sinh(0.02) / 0.02) < 1e-13 and rel(got[2], math.sin(0.4 * PI) / (0.4 * PI)) < 1e-13


def test_time_average_V1_sliding_window_matches_closed_form_and_limits():  # V1/V7 · C01 (12.2)
    t = np.linspace(0.0, 20.0, 20001)
    A, B, tau, om = 1.0, 0.5, 10.0, 2 * PI
    u = A * np.exp(-t / tau) + B * np.cos(om * t)
    for win in (0.4, 1.0, 2.5):
        num = TS.time_average(t, u, win)
        ex = ch12.time_average_exp_cos(t, A, B, tau, om, win)[0]
        ok = np.isfinite(num)
        assert np.all(np.isnan(num[t < win / 2 - 1e-9])) and ok.sum() > 15000
        assert np.max(np.abs(num[ok] - ex[ok])) < 2e-6
    assert abs(ch12.time_average_exp_cos(3.0, A, B, tau, om, 1.0)[2]) < 1e-15        # one whole period kills the wave
    assert rel(ch12.time_average_exp_cos(3.0, A, B, tau, om, 1e-6)[0], A * math.exp(-0.3) + B * math.cos(om * 3.0)) < 1e-10
    assert rel(np.nanmean(TS.time_average(t, u, 1.0, m=2)[10000:10010]),             # mean square over a period
               np.mean((A * np.exp(-t[10000:10010] / tau)) ** 2) + B * B / 2) < 2e-3


def test_time_average_V3_second_order_in_sample_spacing():  # V3 · C01 (12.2)
    errs, hs = [], []
    for n in (201, 401, 801, 1601):
        t = np.linspace(0.0, 4.0, n)
        num = TS.time_average(t, np.exp(-t / 2.0) + 0.3 * np.cos(5.0 * t), 1.0)
        ex = ch12.time_average_exp_cos(t, 1.0, 0.3, 2.0, 5.0, 1.0)[0]
        ok = np.isfinite(num)
        errs.append(np.max(np.abs(num[ok] - ex[ok])))
        hs.append(t[1] - t[0])
    assert abs(observed_order(hs, errs) - 2.0) < 0.15


# ======================================================================================================================
# §12.4  C02 — correlations, integral scale, Taylor microscale (12.12)–(12.19)
# ======================================================================================================================
def test_schwartz_V2_derivation_D02_discriminant_bound():  # V2 · D02 (12.16)
    lam, a, b, c = sp.symbols("lambda a b c", real=True)       # a = mean(v²) > 0, b = mean(uv), c = mean(u²)
    q = c + 2 * lam * b + lam ** 2 * a                          # mean((u + λv)²) ≥ 0 for every real λ
    lam_min = sp.solve(sp.diff(q, lam), lam)[0]
    assert z0(q.subs(lam, lam_min) - (c - b ** 2 / a))          # minimum ≥ 0  ⇔  b² ≤ a c  ⇔  |mean(uv)| ≤ u_rms v_rms
    assert z0(sp.discriminant(q, lam) - 4 * (b ** 2 - a * c))   # at most one real root ⇒ discriminant ≤ 0
    rng = np.random.default_rng(0)
    for _ in range(50):
        u, v = rng.normal(size=(2, 40))
        assert abs(TS.correlation(u, v)) <= math.sqrt(TS.correlation(u, u) * TS.correlation(v, v)) * (1 + 1e-14)
        assert -1.0 <= TS.correlation_coefficient(u, v) <= 1.0
    u = rng.normal(size=100)
    assert rel(TS.correlation_coefficient(u, 3 * u), 1.0) < 1e-14 and rel(TS.correlation_coefficient(u, -0.2 * u), -1.0) < 1e-14
    assert rel(TS.correlation_coefficient(5 * u, 7 * (u + 1e-1 * rng.normal(size=100))),
               TS.correlation_coefficient(u, u + 1e-1 * rng.normal(size=100))) < 0.05   # scale invariance (loose: new noise)
    w = u - u.mean()
    v = rng.normal(size=100)
    v -= v.mean()
    assert rel(TS.correlation_coefficient(w, v), np.corrcoef(w, v)[0, 1]) < 1e-12


def test_correlated_pair_stat_requested_coefficient():  # stat · C04/C05 (N25, N63)
    n = 100_000
    for r in (-0.8, 0.0, 0.45):
        u, v = TS.correlated_pair(n, r, 2.0, 0.5, seed=4)
        assert abs(TS.correlation_coefficient(u, v) - r) < 5 * (1 - r * r) / math.sqrt(n)
        assert abs(np.std(u) - 2.0) < 5 * 2.0 / math.sqrt(2 * n) and abs(np.std(v) - 0.5) < 5 * 0.5 / math.sqrt(2 * n)
    u, v = TS.correlated_pair(1000, 1.0, seed=1)
    assert rel(TS.correlation_coefficient(u, v), 1.0) < 1e-12


def test_autocorrelation_V1_fft_equals_direct_and_cosine_closed_form():  # V1 · C02 (12.17)
    rng = np.random.default_rng(2)
    u = rng.normal(size=600)
    for unb in (True, False):
        l1, R1 = TS.autocorrelation(u, 0.1, max_lag=100, method="fft", unbiased=unb)
        l2, R2 = TS.autocorrelation(u, 0.1, max_lag=100, method="direct", unbiased=unb)
        assert np.allclose(l1, l2) and np.max(np.abs(R1 - R2)) < 1e-10
    w = u - u.mean()
    assert rel(R1[0], np.mean(w * w)) < 1e-12
    hand = np.array([np.sum(w[: 600 - k] * w[k:]) / (600 - k) for k in range(5)])      # the notebook's loop over lags
    assert np.max(np.abs(TS.autocorrelation(u, 0.1, max_lag=100)[1][:5] - hand)) < 1e-12
    dt, om = 0.01, 2 * PI
    t = np.arange(0, 200.0, dt)                                                        # whole number of periods
    lag, R = TS.autocorrelation(3.0 * np.cos(om * t + 0.3), dt, max_lag=300)
    # R = ½A² cos ωτ; the finite-record remainder ½A² mean(cos(2ωt + ωτ + 2φ)) over n − k samples is bounded by
    # ½A²/((n − k) sin ω dt) (geometric sum) — the tolerance is that bound, not a fitted number
    assert np.max(np.abs(R - 4.5 * np.cos(om * lag))) < 4.5 / ((t.size - 300) * math.sin(om * dt))
    # cross-correlation: delayed copy peaks at the delay; R_uv(τ) = R_vu(−τ) (D03)
    s = TS.smooth_signal(4096, 0.05, TS.correlation_spectrum_pair("gaussian", 1.0, 0.4)["S"], seed=3)
    d = 17
    lg, Ruv = TS.cross_correlation(s, np.roll(s, d), 0.05, max_lag=200)
    assert rel(lg[np.argmax(Ruv)], d * 0.05) < 1e-12
    lg2, Rvu = TS.cross_correlation(np.roll(s, d), s, 0.05, max_lag=200)
    assert np.max(np.abs(Ruv - Rvu[::-1])) < 1e-10
    la, Ra = TS.autocorrelation(s, 0.05, max_lag=200)
    lc, Rc = TS.cross_correlation(s, s, 0.05, max_lag=200)
    assert np.max(np.abs(Rc[lc >= -1e-12] - Ra)) < 1e-10 and np.max(np.abs(Rc - Rc[::-1])) < 1e-10   # R_11 even


def test_autocorrelation_stat_ou_record_error_shrinks_as_inverse_root_length():  # stat/V3 · C02 (12.17), (12.18)
    sig, tc, dt = 1.0, 0.5, 0.05
    errs, Ts = [], []
    nseed = 16
    for nt in (10_000, 40_000, 160_000):
        e = []
        for seed in range(nseed):
            u = TS.make_ensemble(1, np.arange(nt) * dt, 0.0, sig, tc, seed=100 + seed)[0]
            lag, R = TS.autocorrelation(u, dt, max_lag=40)
            e.append(np.mean((R - sig ** 2 * np.exp(-lag / tc)) ** 2))
        errs.append(math.sqrt(np.mean(e)))
        Ts.append(nt * dt)
    # statistical convergence ∝ T^(−1/2).  Each rms error is itself a sample (relative scatter ≈ 0.5/√nseed), so the
    # fitted slope over a factor 16 in T has a standard error ≈ √2·(0.5/√nseed)/ln 16 ≈ 0.064; tolerance 4 s.e.
    assert abs(np.polyfit(np.log(Ts), np.log(errs), 1)[0] + 0.5) < 0.25
    assert errs[-1] < 5 * sig ** 2 * math.sqrt(2 * tc / Ts[-1])            # and of the size 2Λ/T predicts
    lag, R = TS.autocorrelation(u, dt, max_lag=400)
    Lam = TS.integral_scale(lag, R)                                        # to the first zero (noisy tail)
    assert abs(Lam - tc) < 5 * tc * math.sqrt(2 * tc / Ts[-1]) * 4          # Λ error ~ several × relative R error
    assert rel(TS.effective_samples(Ts[-1], 0.5), Ts[-1] / 0.5) < 1e-15


def test_integral_and_micro_scales_V1_closed_form_pairs():  # V1 · C02 (12.18), (12.19) · D04
    tau = np.linspace(0.0, 12.0, 24001)
    ex = TS.correlation_spectrum_pair("exponential", 1.3, 0.7)
    ga = TS.correlation_spectrum_pair("gaussian", 1.3, 0.7)
    dc = TS.correlation_spectrum_pair("damped_cosine", 1.3, 0.7, 6.0)
    assert rel(TS.integral_scale(tau, ex["r"](tau), upto="all"), 0.7) < 1e-6 and ex["Lambda_t"] == 0.7
    assert rel(TS.integral_scale(tau, ga["r"](tau), upto="all"), math.sqrt(PI) / 2 * 0.7) < 1e-8
    assert rel(ga["Lambda_t"], math.sqrt(PI) / 2 * 0.7) < 1e-14 and ga["lambda_t"] == 0.7
    assert rel(TS.integral_scale(tau, dc["r"](tau), upto="all"), 0.7 / (1 + (6.0 * 0.7) ** 2)) < 1e-5
    assert rel(dc["t_c"], PI / 12.0) < 1e-14 and math.isinf(ex["t_c"]) and math.isnan(ex["lambda_t"])
    assert rel(TS.correlation_time(tau, dc["r"](tau)), PI / 12.0) < 1e-6
    assert math.isinf(TS.correlation_time(tau, ex["r"](tau)))
    fine = np.linspace(0.0, 0.2, 201)
    assert rel(TS.taylor_microscale(fine, ga["r"](fine)), 0.7) < 1e-6                  # λ_t = τ_c for the Gaussian
    # an exponential (OU) correlation has a cusp: the "microscale" shrinks with the sample spacing instead of converging
    lam = [TS.taylor_microscale(np.arange(8) * h, ex["r"](np.arange(8) * h)) for h in (0.02, 0.01, 0.005)]
    assert lam[0] > 1.3 * lam[1] > 1.3 ** 2 * lam[2]
    for p in (ex, ga, dc):
        assert rel(p["S0"], p["variance"] * p["Lambda_t"] / PI) < 1e-13                # S_e(0) = σ²Λ_t/π
        assert rel(TS.integral_scale_from_spectrum(p["S"](0.0), p["variance"]), p["Lambda_t"]) < 1e-13
        assert rel(p["R"](0.3), p["variance"] * p["r"](0.3)) < 1e-14 and p["r"](0.0) == 1.0


def test_taylor_microscale_V2_derivation_D04_osculating_parabola_and_derivative_form():  # V2 · D04 (12.19)
    tau, lam = sp.symbols("tau lambda", positive=True)
    a1, a2, a3, w1, w2, w3 = sp.symbols("a1 a2 a3 omega1 omega2 omega3", positive=True)
    # a stationary signal Σ a_n cos(ω_n t + φ_n) with random phases: R(τ) = Σ ½a_n² cos ω_n τ
    R = sum(a ** 2 / 2 * sp.cos(w * tau) for a, w in ((a1, w1), (a2, w2), (a3, w3)))
    u2 = R.subs(tau, 0)
    r = R / u2
    assert sp.diff(r, tau).subs(tau, 0) == 0                                  # even ⇒ no linear term
    lam2 = -2 / sp.diff(r, tau, 2).subs(tau, 0)                                 # (12.19)
    assert z0(sp.series(r, tau, 0, 3).removeO() - (1 - tau ** 2 / lam2))        # the osculating parabola 1 − τ²/λ²
    dudt2 = sum(a ** 2 * w ** 2 / 2 for a, w in ((a1, w1), (a2, w2), (a3, w3)))  # mean((du/dt)²), phase-averaged
    assert z0(lam2 - 2 * u2 / dudt2)                                            # λ_t² = 2 mean(u²)/mean((du/dt)²)
    # numerically, on a smooth random record with a Gaussian correlation (λ_t = τ_c)
    pair = TS.correlation_spectrum_pair("gaussian", 1.0, 0.4)
    u = TS.smooth_signal(2 ** 17, 0.01, pair["S"], seed=7)
    lag, Rn = TS.autocorrelation(u, 0.01, max_lag=60)
    lam_r = TS.taylor_microscale(lag, Rn, fit_points=6)
    k = 2 * PI * np.fft.rfftfreq(u.size, 0.01)
    du = np.fft.irfft(1j * k * np.fft.rfft(u), n=u.size)
    lam_d = math.sqrt(2 * np.mean(u * u) / np.mean(du * du))
    assert rel(lam_r, lam_d) < 5e-3 and rel(lam_d, 0.4) < 0.05                 # the two definitions agree; ≈ τ_c


# ======================================================================================================================
# §12.4  C03 — spectrum ↔ correlation (12.20)–(12.23)
# ======================================================================================================================
def test_spectrum_pair_V2_derivation_D05_fourier_transforms_and_normalisation():  # V2 · D05 (12.20)–(12.22)
    tau, om, tc, sg = sp.symbols("tau omega tau_c sigma", positive=True)
    S_exp = 2 * sp.integrate(sg ** 2 * sp.exp(-tau / tc) * sp.cos(om * tau), (tau, 0, sp.oo)) / (2 * sp.pi)
    assert z0(S_exp - sg ** 2 * tc / (sp.pi * (1 + om ** 2 * tc ** 2)))
    S_g = 2 * sp.integrate(sg ** 2 * sp.exp(-tau ** 2 / tc ** 2) * sp.cos(om * tau), (tau, 0, sp.oo)) / (2 * sp.pi)
    assert z0(S_g - sg ** 2 * tc * sp.exp(-om ** 2 * tc ** 2 / 4) / (2 * sp.sqrt(sp.pi)))
    for S, Lam in ((S_exp, tc), (S_g, sp.sqrt(sp.pi) * tc / 2)):
        assert z0(2 * sp.integrate(S, (om, 0, sp.oo)) - sg ** 2)                # (12.22): ∫S dω = variance
        assert z0(S.subs(om, 0) - sg ** 2 * Lam / sp.pi)                        # S_e(0) = σ²Λ_t/π
    # the functions return these closed forms
    w = np.linspace(0.0, 30.0, 301)
    for kind, expr in (("exponential", S_exp), ("gaussian", S_g)):
        f = sp.lambdify(om, expr.subs({sg: 1.3, tc: 0.7}), "numpy")
        assert np.max(np.abs(TS.correlation_spectrum_pair(kind, 1.3, 0.7)["S"](w) - f(w))) < 1e-14
    dc = TS.correlation_spectrum_pair("damped_cosine", 1.3, 0.7, 6.0)
    num = quad(lambda s: dc["R"](s) * math.cos(4.0 * s), 0, 60, limit=800)[0] / PI
    assert rel(dc["S"](4.0), num) < 1e-8


def test_spectrum_from_correlation_V1_transform_roundtrip_and_variance():  # V1 · C03 (12.20), (12.21), (12.22)
    tau = np.linspace(0.0, 25.0, 5001)
    w = np.linspace(0.0, 40.0, 2001)
    h = tau[1] - tau[0]
    for kind in ("exponential", "gaussian"):
        p = TS.correlation_spectrum_pair(kind, 1.3, 0.7)
        S = TS.spectrum_from_correlation(tau, p["R"](tau), w)
        # Filon rule = exact transform of the piecewise-linear interpolant: |error| ≤ (h²/8)·∫|R″|dτ/π (interpolation
        # bound; the cusp of the exponential adds no error because τ = 0 is a node)
        d2 = np.abs(np.gradient(np.gradient(p["R"](tau), h), h))[2:-2]
        assert np.max(np.abs(S - p["S"](w))) < (h * h / 8) * np.trapezoid(d2, dx=h) / PI
        S2 = TS.spectrum_from_correlation(tau[::2], p["R"](tau[::2]), w[:400])
        ratio = np.max(np.abs(S2 - p["S"](w[:400]))) / np.max(np.abs(S[:400] - p["S"](w[:400])))
        assert abs(math.log2(ratio) - 2.0) < 0.15                                    # V3: second order in the lag step
        St = TS.spectrum_from_correlation(tau, p["R"](tau), w[:200], rule="trapezoid")
        assert np.max(np.abs(St - p["S"](w[:200]))) < 5e-5 * p["S0"]
        both = TS.spectrum_from_correlation(np.concatenate([-tau[:0:-1], tau]), p["R"](np.abs(np.concatenate([-tau[:0:-1], tau]))), w[:50])
        assert np.max(np.abs(both - S[:50])) < 1e-9
    g = TS.correlation_spectrum_pair("gaussian", 1.3, 0.7)
    wg = np.linspace(0.0, 30.0, 6001)
    lags = np.array([0.0, 0.2, 0.7, 1.5])
    Rb = TS.correlation_from_spectrum(wg, g["S"](wg), lags)
    hw = wg[1] - wg[0]
    d2S = np.abs(np.gradient(np.gradient(g["S"](wg), hw), hw))[2:-2]
    assert np.max(np.abs(Rb - g["R"](lags))) < 2 * (hw * hw / 8) * np.trapezoid(d2S, dx=hw)   # (12.21), no 1/2π here    assert rel(TS.spectrum_variance(wg, g["S"](wg)), 1.69) < 1e-8                    # (12.22) two-sided, half given
    assert rel(TS.spectrum_variance(wg, 2 * g["S"](wg), two_sided=False), 1.69) < 1e-8
    wl = np.linspace(-4000.0, 4000.0, 800001)
    e = TS.correlation_spectrum_pair("exponential", 1.3, 0.7)
    assert rel(TS.spectrum_variance(wl, e["S"](wl)), 1.69) < 3e-4                    # Lorentzian tail: slow
    # frozen-turbulence change of variable keeps the variance (N40)
    k1, S11 = TS.frequency_to_wavenumber_spectrum(wg, g["S"](wg), 8.0)
    assert rel(TS.spectrum_variance(k1, S11), TS.spectrum_variance(wg, g["S"](wg))) < 1e-12
    assert np.allclose(k1, wg / 8.0) and np.allclose(S11, 8.0 * g["S"](wg))
    x, uu = TS.taylor_frozen(np.arange(5.0), np.arange(5.0) ** 2, 8.0)
    assert np.allclose(x, 8.0 * np.arange(5.0)) and np.allclose(uu, np.arange(5.0) ** 2)


def test_periodogram_V4_parseval_and_welch_cross_check():  # V4 · C03 (N38) (12.22)
    rng = np.random.default_rng(8)
    u = rng.normal(size=4096) + 0.7
    d = 0.02
    om, S, w = TS.periodogram(u, d, return_weights=True)
    var = np.var(u)
    assert rel(TS.spectrum_variance(om, S, weights=w), var) < 1e-12                  # Parseval, exactly
    assert rel(np.sum(w * S), var) < 1e-12 and rel(om[1], 2 * PI / (4096 * d)) < 1e-12
    om1, S1, w1 = TS.periodogram(u, d, two_sided=False, return_weights=True)
    om2, S2 = TS.periodogram(u, d, one_sided=True)
    assert np.allclose(S1, S2) and rel(np.sum(w1 * S1), var) < 1e-12
    assert np.allclose(S1[1:-1], 2 * S[1:-1], rtol=1e-12)                            # one-sided = 2 × two-sided
    assert abs(TS.spectrum_variance(om, S) / var - 1) < 5e-3                         # trapezoid agrees to O(1/n)
    # from scratch: np.fft.rfft scaled by hand (the notebook's cell)
    uh = np.fft.rfft(u - u.mean())
    S_hand = np.abs(uh) ** 2 * d / (2 * PI * u.size)
    assert np.allclose(S[1:-1], S_hand[1:-1], rtol=1e-10)
    # a sine puts ½A² in its own bin
    t = np.arange(4096) * d
    j = 100
    omj = 2 * PI * j / (4096 * d)
    o3, S3, w3 = TS.periodogram(2.0 * np.sin(omj * t), d, return_weights=True)
    assert np.argmax(S3) == j and rel(w3[j] * S3[j], 2.0) < 1e-10
    # scipy.signal.welch in cyclic frequency, one-sided: S(ω) two-sided = P(f)/(4π)
    f, P = welch(u, fs=1 / d, window="boxcar", nperseg=512, noverlap=0, detrend="constant", scaling="density")
    o4, S4 = TS.periodogram(u, d, segments=8)
    assert np.allclose(o4, 2 * PI * f) and np.max(np.abs(S4[1:-1] / (P[1:-1] / (4 * PI)) - 1)) < 0.02


def test_periodogram_stat_ou_record_gives_lorentzian_and_smooth_signal_gives_its_spectrum():  # stat · C03
    sig, tc, dt, n, segs = 1.0, 0.5, 0.05, 2 ** 18, 256
    u = TS.make_ensemble(1, np.arange(n) * dt, 0.0, sig, tc, seed=21)[0]
    om, S = TS.periodogram(u, dt, segments=segs)
    exact = TS.correlation_spectrum_pair("exponential", sig, tc)["S"](om)
    band = (om > 0.5) & (om < 8.0)                         # resolved, far below Nyquist (aliasing < 1 %)
    zsc = (S[band] / exact[band] - 1) * math.sqrt(segs)
    assert np.max(np.abs(zsc)) < 5 + 0.02 * math.sqrt(segs) and abs(np.mean(zsc)) < 5 / math.sqrt(band.sum()) + 0.3
    pair = TS.correlation_spectrum_pair("gaussian", 1.5, 0.4)
    s = TS.smooth_signal(2 ** 17, 0.02, pair["S"], seed=5)
    assert abs(s.mean()) < 1e-10 and abs(np.var(s) - 2.25) < 5 * 2.25 * math.sqrt(2 * 2 * pair["Lambda_t"] / (2 ** 17 * 0.02))
    lag, R = TS.autocorrelation(s, 0.02, max_lag=40)
    assert np.max(np.abs(R / R[0] - pair["r"](lag))) < 0.05


def test_spatial_correlation_V1_direct_shift_and_synthetic_field_statistics():  # V1/V4 · C03 (12.23), N02
    L, n = 2 * PI, 64
    E = lambda K: K ** 4 * np.exp(-K ** 2 / 18.0)  # noqa: E731
    u, v = TS.synthetic_solenoidal_field(n, L, E, seed=2, dim=2)
    dx = L / n
    r, R = TS.spatial_correlation(u, v, dx, axis=-1)
    direct = np.array([np.mean(u * np.roll(v, -k, axis=-1)) for k in range(len(r))])
    assert np.max(np.abs(R - direct)) < 1e-10 * np.max(np.abs(direct)) and rel(R[0], np.mean(u * v)) < 1e-9
    assert np.allclose(r, np.arange(len(r)) * dx)
    urms = math.sqrt(np.mean(u * u + v * v) / 2)
    assert ch12.divergence_rms(u, v, dx) < 1e-11 * urms / dx                          # V4: solenoidal to round-off
    a, b = TS.white_noise_field(n, seed=2)
    assert ch12.divergence_rms(a, b, dx) > 0.5 / dx                                   # noise is not a flow
    assert 0 < ch12.divergence_rms(u, v, dx, method="central") < 0.2 * urms / dx
    K, Es = TS.shell_spectrum((u, v), L)
    assert rel(np.sum(Es) * (2 * PI / L), 0.5 * np.mean(u * u + v * v)) < 0.02       # Σ E ΔK = kinetic energy
    assert abs(K[np.argmax(Es)] - 6.0) <= 2.0                                         # peak of K⁴e^{−K²/18} at K = 6
    ens = [TS.synthetic_solenoidal_field(32, L, E, seed=s, dim=2) for s in range(6)]
    su, sv = np.stack([e[0] for e in ens]), np.stack([e[1] for e in ens])
    rep = ch12.mean_divergence(su, sv, L / 32, L / 32, report=True)                   # (12.27), (12.28)
    scale = np.sqrt(np.mean(su ** 2)) / (L / 32)
    assert rep["mean_rms"] < 1e-11 * scale and rep["fluctuation_rms"] < 1e-11 * scale
    assert rep["commutation_residual"] < 1e-11 * scale
    assert ch12.mean_divergence(su, sv, L / 32, L / 32).shape == (32, 32)
    rc = ch12.mean_divergence(su, sv, L / 32, L / 32, method="central", report=True)
    assert rc["commutation_residual"] < 1e-11 * scale < rc["member_rms"]              # (12.8) holds for any linear operator


def test_frozen_field_probe_V7_error_vanishes_and_grows_linearly():  # V7 · C03 (N40)
    F = np.sin(2 * PI * 4 * np.arange(512) / 512) + 0.5 * np.cos(2 * PI * 9 * np.arange(512) / 512)
    assert ch12.frozen_field_probe(F, 10.0, 0.0, 0.01)["reconstruction_error"] < 1e-12
    ratios = np.array([0.002, 0.004, 0.008])
    errs = [ch12.frozen_field_probe(F, 10.0, 10.0 * q, 0.01, n_realizations=256, seed=1)["reconstruction_error"] for q in ratios]
    assert abs(np.polyfit(np.log(ratios), np.log(errs), 1)[0] - 1.0) < 0.2
    p = ch12.frozen_field_probe(F, 10.0, 0.05, 0.01)
    assert p["record"].shape == p["t"].shape and rel(p["ratio"], 0.005) < 1e-12 and np.allclose(p["x_frozen"], 10.0 * p["t"])


# ======================================================================================================================
# An explicit finite ensemble for the averaged equations (independent of the chapter's sympy engines)
# ======================================================================================================================
_X = sp.symbols("x1 x2 x3", real=True)
_T = sp.Symbol("t", real=True)
_NU, _RHO0, _GR, _AL, _T0, _KTH = sp.symbols("nu rho0 g alpha T0 kappa_th", positive=True)
_ENS: dict = {}


def _curl(A):
    x1, x2, x3 = _X
    return [sp.diff(A[2], x2) - sp.diff(A[1], x3), sp.diff(A[0], x3) - sp.diff(A[2], x1), sp.diff(A[1], x1) - sp.diff(A[0], x2)]


def _ensemble() -> dict:
    """Four-member ensemble: mean fields + s1·a + s2·b + s1·s2·c with signs s1, s2 = ±1 (exact averages by summing)."""
    if _ENS:
        return _ENS
    x1, x2, x3 = _X
    t = _T
    U = _curl((x2 * x3 ** 2 / 2 + t * x2, x1 ** 2 * x3 / 3, x1 * x2 ** 2 / 2 - x3 * x1))
    a = _curl((x2 ** 2 * x3, x3 * x1 + t * x1 * x2, x1 * x2))
    b = _curl((x3 ** 2 + x1 * x2 * x3, x1 ** 2 * x2, x2 * x3))
    c = _curl((x1 * x2, x2 * x3 ** 2, t * x3 * x1 + x1 ** 2 * x2))
    P = x1 ** 2 * x2 - x3 * x2 + t * x1
    q = (x1 * x2, x3 ** 2 + x1, x2 * x3 * x1)
    Tm = 2 * x3 + x1 * x2 + t * x3 ** 2
    th = (x1 + x3 ** 2, x2 * x3, x1 * x3 + t * x2)
    members = []
    for s1, s2 in itertools.product((1, -1), repeat=2):
        u = [s1 * a[i] + s2 * b[i] + s1 * s2 * c[i] for i in range(3)]
        members.append(dict(u=u, p=s1 * q[0] + s2 * q[1] + s1 * s2 * q[2], T=s1 * th[0] + s2 * th[1] + s1 * s2 * th[2]))

    def avg(fn):
        return sp.expand(sum(fn(m) for m in members) / 4)

    def D(f, j):
        return sp.diff(f, _X[j])

    def lap(f):
        return sum(sp.diff(f, xx, 2) for xx in _X)

    def ns_residual(m, i):
        """Boussinesq momentum equation of one member, everything on the left: 0 for a real flow; here a forcing."""
        ut = [U[k] + m["u"][k] for k in range(3)]
        out = sp.diff(ut[i], t) + sum(ut[j] * D(ut[i], j) for j in range(3)) + D(P + m["p"], i) / _RHO0 - _NU * lap(ut[i])
        if i == 2:
            out += _GR * (1 - _AL * (Tm + m["T"] - _T0))
        return out

    def heat_residual(m):
        ut = [U[k] + m["u"][k] for k in range(3)]
        Tt = Tm + m["T"]
        return sp.diff(Tt, t) + sum(ut[j] * D(Tt, j) for j in range(3)) - _KTH * lap(Tt)

    _ENS.update(U=U, P=P, Tm=Tm, members=members, avg=avg, D=D, lap=lap, ns=ns_residual, heat=heat_residual)
    return _ENS


def test_rans_V2_derivation_D06_reynolds_averaged_momentum_from_a_finite_ensemble():  # V2 · C04 · D06 (12.27)–(12.30)
    E = _ensemble()
    U, P, Tm, avg, D = E["U"], E["P"], E["Tm"], E["avg"], E["D"]
    assert sp.expand(sum(D(U[i], i) for i in range(3))) == 0                              # (12.27)
    for m in E["members"]:
        assert sp.expand(sum(D(m["u"][i], i) for i in range(3))) == 0                     # (12.28)
    assert all(avg(lambda m, i=i: m["u"][i]) == 0 for i in range(3))                      # (12.26)
    uu = [[avg(lambda m, i=i, j=j: m["u"][i] * m["u"][j]) for j in range(3)] for i in range(3)]
    assert uu[0][1] != 0 and uu[0][0] != 0
    kd = lambda i, j: 1 if i == j else 0  # noqa: E731
    tau = [[-P * kd(i, j) + _RHO0 * _NU * (D(U[i], j) + D(U[j], i)) - _RHO0 * uu[i][j] for j in range(3)] for i in range(3)]
    for i in range(3):
        lhs = sp.diff(U[i], _T) + sum(U[j] * D(U[i], j) for j in range(3))
        rhs = -_GR * (1 - _AL * (Tm - _T0)) * kd(i, 2) + sum(D(tau[i][j], j) for j in range(3)) / _RHO0
        averaged = avg(lambda m, i=i: E["ns"](m, i))
        assert sp.expand(averaged - (lhs - rhs)) == 0                                     # (12.30)
        no_stress = rhs + sum(D(uu[i][j], j) for j in range(3))                           # "mean of a product = product of means"
        assert sp.expand(averaged - (lhs - no_stress)) != 0                               # … leaves the Reynolds stress out


def test_mean_energy_V2_derivation_D09_budget_of_the_mean_flow():  # V2 · C06 · D09 (12.46)
    E = _ensemble()
    U, P, Tm, avg, D = E["U"], E["P"], E["Tm"], E["avg"], E["D"]
    uu = [[avg(lambda m, i=i, j=j: m["u"][i] * m["u"][j]) for j in range(3)] for i in range(3)]
    S = [[(D(U[i], j) + D(U[j], i)) / 2 for j in range(3)] for i in range(3)]
    Ebar = sum(Ui ** 2 for Ui in U) / 2
    lhs = sp.diff(Ebar, _T) + sum(U[j] * D(Ebar, j) for j in range(3))
    flux = [-U[j] * P / _RHO0 + 2 * _NU * sum(U[i] * S[i][j] for i in range(3)) - sum(uu[i][j] * U[i] for i in range(3))
            for j in range(3)]
    exchange = sum(uu[i][j] * D(U[i], j) for i in range(3) for j in range(3))
    rho_bar = _RHO0 * (1 - _AL * (Tm - _T0))

    def rhs(sign_exchange=1, diss=2):
        return (sum(D(flux[j], j) for j in range(3)) - diss * _NU * sum(S[i][j] ** 2 for i in range(3) for j in range(3))
                + sign_exchange * exchange - _GR / _RHO0 * rho_bar * U[2])

    derived = sum(U[i] * avg(lambda m, i=i: E["ns"](m, i)) for i in range(3))             # U_i × (12.30)
    assert sp.expand(derived - (lhs - rhs())) == 0
    assert sp.expand(derived - (lhs - rhs(sign_exchange=-1))) != 0                        # sign of the exchange term matters
    assert sp.expand(derived - (lhs - rhs(diss=1))) != 0                                  # and the factor 2 of 2νS̄S̄
    r = ch12.mean_energy_budget_sympy()
    assert r["check"] is True and r["exchange"] != 0


def test_tke_budget_V2_derivation_D10_from_the_fluctuation_equation():  # V2 · C06 · D10 (12.47)
    E = _ensemble()
    U, avg, D = E["U"], E["avg"], E["D"]
    R3 = range(3)
    uu = [[avg(lambda m, i=i, j=j: m["u"][i] * m["u"][j]) for j in R3] for i in R3]
    ebar = sum(uu[i][i] for i in R3) / 2
    lhs = sp.diff(ebar, _T) + sum(U[j] * D(ebar, j) for j in R3)
    Sp = lambda m, i, j: (D(m["u"][i], j) + D(m["u"][j], i)) / 2  # noqa: E731

    def rhs(triple=sp.Rational(1, 2), prod_sign=-1, diss=2):
        flux = [-avg(lambda m, j=j: m["p"] * m["u"][j]) / _RHO0
                + 2 * _NU * avg(lambda m, j=j: sum(m["u"][i] * Sp(m, i, j) for i in R3))
                - triple * avg(lambda m, j=j: sum(m["u"][i] ** 2 for i in R3) * m["u"][j]) for j in R3]
        return (sum(D(flux[j], j) for j in R3)
                - diss * _NU * avg(lambda m: sum(Sp(m, i, j) ** 2 for i in R3 for j in R3))
                + prod_sign * sum(uu[i][j] * D(U[i], j) for i in R3 for j in R3)
                + _GR * _AL * avg(lambda m: m["u"][2] * m["T"]))

    derived = avg(lambda m: sum(m["u"][i] * E["ns"](m, i) for i in R3))        # mean(u_i × member equation)
    assert sp.expand(derived - (lhs - rhs())) == 0                             # (12.47)
    assert avg(lambda m: sum(m["u"][i] ** 2 for i in R3) * m["u"][0]) != 0      # the triple correlation is not trivially 0
    assert sp.expand(derived - (lhs - rhs(triple=1))) != 0                     # ½ mean(u_i²u_j), not mean(u_i²u_j) (slip #10)
    assert sp.expand(derived - (lhs - rhs(prod_sign=1))) != 0                  # production enters with the opposite sign of (12.46)
    assert sp.expand(derived - (lhs - rhs(diss=1))) != 0
    # intermediate line of Part F: the mean equation removes nothing from mean(u_i N_i) because mean(u_i) = 0
    for i in R3:
        mean_eq = avg(lambda m, i=i: E["ns"](m, i))
        assert avg(lambda m, i=i: m["u"][i] * mean_eq) == 0
    r = ch12.tke_budget_sympy()
    assert r["check"] is True and r["viscous_identity"] is True
    rs = ch12.reynolds_stress_budget_sympy()
    assert rs["all_ok"] is True and all(rs["checks"].values()) and len(rs["checks"]) == 6
    assert sp.simplify(rs["trace_half"] - r["derived"]) == 0 or z0(rs["trace_half"] - r["derived"])   # ½ trace of (12.35) = (12.47)


def test_temperature_variance_V2_derivation_budget_and_printed_factor():  # V2 · C15 (N190, a-D41) (12.112) · slip #15
    E = _ensemble()
    U, Tm, avg, D, lap = E["U"], E["Tm"], E["avg"], E["D"], E["lap"]
    R3 = range(3)
    half_var = avg(lambda m: m["T"] ** 2) / 2
    lhs = sp.diff(half_var, _T) + sum(U[j] * D(half_var, j) for j in R3)

    def rhs(mol=sp.Rational(1, 2)):
        return (-sum(avg(lambda m, j=j: m["u"][j] * m["T"]) * D(Tm, j) for j in R3)
                - sum(D(avg(lambda m, j=j: m["T"] ** 2 * m["u"][j]) / 2, j) for j in R3)
                + _KTH * lap(2 * mol * half_var)
                - _KTH * avg(lambda m: sum(D(m["T"], j) ** 2 for j in R3)))

    derived = avg(lambda m: m["T"] * E["heat"](m))
    assert sp.expand(derived - (lhs - rhs())) == 0                 # κ ∂²(½ mean(T′²)) — the derived form
    assert sp.expand(derived - (lhs - rhs(mol=1))) != 0            # the page's κ ∂ mean(T′²)/∂z (no ½) does not close
    z = np.linspace(0.0, 50.0, 26)
    b = ch12.temperature_variance_budget(z, 300.0 - 0.02 * z, 0.05 * np.ones_like(z), 8e-4 * np.ones_like(z))
    assert np.allclose(b["production"], 0.05 * 0.02) and np.all(b["dissipation"] < 0) and b["transport_inferred"]
    assert np.max(np.abs(b["production"] + b["dissipation"] + b["transport"])) < 1e-15   # inferred transport closes it
    assert np.all(ch12.temperature_variance_budget(z, 300.0 + 0.02 * z, 0.05 * np.ones_like(z), 0 * z)["production"] < 0)


# ======================================================================================================================
# §12.5  C04 — Reynolds stress and the averaged equations (12.24)–(12.35)
# ======================================================================================================================
def test_rans_engine_V2_checks_and_term_structure():  # V2 · C04 (12.27)–(12.34)
    r = ch12.rans_sympy()
    assert r["checks"] == {"eq_12_27": True, "eq_12_28": True, "eq_12_30": True, "eq_12_31": True, "eq_12_34": True}
    assert len(r["mean_momentum"]) == 3 and len(r["fluctuation_momentum"]) == 3 and len(r["steps"]) >= 4
    s = ch12.sympy_summary("overlap_matching", cache=False)
    assert s["name"] == "overlap_matching" and all(s["checks"].values()) and isinstance(s["latex"], dict)
    assert set(ch12.SYMPY_ENGINES) >= {"rans", "tke_budget", "plane_jet_similarity"}


def _mms_rans_3d():
    x, y, z, t = sp.symbols("x y z t")
    nu, rho0, g, al, T0 = 0.07, 1.3, 9.0, 0.02, 1.0
    U = [sp.sin(x) * sp.cos(y) * (1 + t), -sp.cos(x) * sp.sin(y) * (1 + t) + 0.2 * sp.sin(z), 0.3 * sp.cos(x + y)]
    P = sp.cos(x) * sp.sin(z) + 0.5 * y ** 2
    T = 1.0 + 0.3 * sp.sin(x + z)
    uu = sp.Matrix(3, 3, lambda i, j: 0)
    base = [[0.5 + 0.1 * sp.sin(x + y), 0.05 * sp.cos(y * z), 0.02 * sp.sin(x)],
            [0, 0.4 + 0.1 * sp.cos(z), 0.03 * sp.cos(x - z)],
            [0, 0, 0.3 + 0.05 * sp.sin(y)]]
    for i in range(3):
        for j in range(3):
            uu[i, j] = base[min(i, j)][max(i, j)]
    X = (x, y, z)
    res = []
    for i in range(3):
        lhs = sp.diff(U[i], t) + sum(U[j] * sp.diff(U[i], X[j]) for j in range(3))
        tau = [-P * (1 if i == j else 0) + rho0 * nu * (sp.diff(U[i], X[j]) + sp.diff(U[j], X[i])) - rho0 * uu[i, j]
               for j in range(3)]
        rhs = -g * (1 - al * (T - T0)) * (1 if i == 2 else 0) + sum(sp.diff(tau[j], X[j]) for j in range(3)) / rho0
        res.append(lhs - rhs)
    f = lambda e: sp.lambdify((x, y, z, t), e, "numpy")  # noqa: E731
    fU, fP, fT, fuu, fres = f(U), f(P), f(T), f(uu), f(res)
    fields = {"U": lambda p, tt: np.array(fU(*p, tt), float), "P": lambda p, tt: float(fP(*p, tt)),
              "T": lambda p, tt: float(fT(*p, tt)), "uu": lambda p, tt: np.array(fuu(*p, tt), float)}
    return fields, (lambda p, tt: np.array(fres(*p, tt), float)), dict(nu=nu, rho0=rho0, g=g, alpha=al, T0=T0)


def test_rans_momentum_residual_V3_manufactured_fields_second_order():  # V1/V3 · C04 (12.30) — contract-only name
    fields, exact, kw = _mms_rans_3d()
    pt, t0 = np.array([0.4, 0.9, -0.3]), 0.2
    ex = exact(pt, t0)
    assert np.max(np.abs(ex)) > 0.1
    hs = [0.08, 0.04, 0.02, 0.01]
    errs = [np.max(np.abs(ch12.rans_momentum_residual(fields, pt, h=h, ht=h, t=t0, **kw) - ex)) for h in hs]
    assert abs(observed_order(hs, errs) - 2.0) < 0.15 and errs[-1] < 1e-4
    no_g = ch12.rans_momentum_residual(fields, pt, h=0.01, ht=0.01, t=t0, nu=kw["nu"], rho0=kw["rho0"])
    Tb = fields["T"](pt, t0)
    assert abs((ex - no_g)[2] - kw["g"] * (1 - kw["alpha"] * (Tb - kw["T0"]))) < 1e-4      # buoyancy on the last axis only
    assert np.max(np.abs((ex - no_g)[:2])) < 1e-4


def test_rans_2d_residual_V3_manufactured_fields_second_order():  # V1/V3 · C09 (12.58)–(12.60)
    x, y = sp.symbols("x y")
    nu, rho = 0.03, 1.2
    psi = sp.sin(x) * sp.cos(2 * y) + 0.3 * x * y
    U, V = sp.diff(psi, y), -sp.diff(psi, x)
    P = sp.cos(x) * y + x ** 2
    uu, uv, vv = 0.4 + 0.1 * sp.sin(x * y), 0.05 * sp.cos(x - y), 0.3 + 0.1 * sp.cos(y)
    rx = U * sp.diff(U, x) + V * sp.diff(U, y) + sp.diff(P, x) / rho - nu * (sp.diff(U, x, 2) + sp.diff(U, y, 2)) \
        + sp.diff(uu, x) + sp.diff(uv, y)
    ry = U * sp.diff(V, x) + V * sp.diff(V, y) + sp.diff(P, y) / rho - nu * (sp.diff(V, x, 2) + sp.diff(V, y, 2)) \
        + sp.diff(uv, x) + sp.diff(vv, y)
    L = lambda e: sp.lambdify((x, y), e, "numpy")  # noqa: E731
    fs = [L(e) for e in (U, V, P, uu, uv, vv)]
    ex = np.array([L(rx)(0.7, 0.4), L(ry)(0.7, 0.4)])
    hs = [0.08, 0.04, 0.02, 0.01]
    errs = []
    for h in hs:
        r = ch12.rans_2d_residual(*fs, 0.7, 0.4, nu, rho, h=h, return_continuity=True)
        errs.append(max(abs(r[0] - ex[0]), abs(r[1] - ex[1])))
        assert abs(r[2]) < 5 * h * h                                                  # (12.58): solenoidal by construction
    assert abs(observed_order(hs, errs) - 2.0) < 0.15 and errs[-1] < 1e-3
    assert len(ch12.rans_2d_residual(*fs, 0.7, 0.4, nu, rho)) == 2


def test_reynolds_stress_V1_covariance_sign_rotation_and_mean_stress():  # V1 · C04 (12.30), N49
    rng = np.random.default_rng(6)
    A = np.array([[1.0, 0.0, 0.0], [-0.6, 0.8, 0.0], [0.2, 0.1, 0.5]])
    s = rng.normal(size=(5000, 3)) @ A.T + np.array([3.0, 0.0, -1.0])
    cov = TS.velocity_covariance(s)
    fl = s - s.mean(axis=0)
    hand = np.zeros((3, 3))
    for row in fl:                                    # the notebook's explicit sums
        hand += np.outer(row, row)
    assert np.allclose(cov, hand / 5000, rtol=1e-12) and np.allclose(cov, np.cov(s.T, bias=True), rtol=1e-12)
    R = TS.reynolds_stress(s, rho0=1.2)
    assert np.allclose(R, -1.2 * cov) and np.allclose(R, R.T)
    assert cov[0, 1] < 0 and R[0, 1] > 0                                    # mean(uv) < 0 ⇒ positive shear stress
    assert rel(TS.turbulent_kinetic_energy(s), 0.5 * np.trace(cov)) < 1e-13
    th = 0.7
    C = np.array([[math.cos(th), -math.sin(th), 0], [math.sin(th), math.cos(th), 0], [0, 0, 1]])
    assert np.allclose(TS.velocity_covariance(s @ C), C.T @ cov @ C, atol=1e-12)        # V7 tensor transformation
    b = TS.anisotropy_tensor(cov)
    assert abs(np.trace(b)) < 1e-14 and np.allclose(b, b.T)
    assert np.max(np.abs(TS.anisotropy_tensor(2.5 * np.eye(3)))) < 1e-15
    gradU = np.zeros((3, 3))
    gradU[0, 1] = 4.0
    tau = ch12.mean_stress_tensor(101.0, gradU, 1.8e-5, 1.2, cov)
    assert rel(tau[0, 1], 1.8e-5 * 4.0 - 1.2 * cov[0, 1]) < 1e-13 and np.allclose(tau, tau.T)
    iso = ch12.mean_stress_tensor(101.0, gradU, 1.8e-5, 1.2, 0.5 * np.eye(3))
    assert rel(iso[1, 1], -(101.0 + 0.6)) < 1e-14 and rel(iso[0, 1], 7.2e-5) < 1e-13     # isotropic part acts as a pressure
    assert rel(ch12.shear_production(cov, gradU), -cov[0, 1] * 4.0) < 1e-13
    Pij = ch12.reynolds_stress_production(cov, gradU)
    assert rel(0.5 * np.trace(Pij), ch12.shear_production(cov, gradU)) < 1e-13
    assert rel(Pij[0, 0], -2 * cov[0, 1] * 4.0) < 1e-13 and abs(Pij[1, 1]) + abs(Pij[2, 2]) < 1e-15
    W = np.array([[0.3, 1.0, 0.0], [-2.0, 0.1, 0.5], [0.7, 0.0, -0.4]])                   # trace 0
    assert abs(ch12.shear_production(0.9 * np.eye(3), W)) < 1e-15                         # isotropic stress: no production


def test_displaced_parcel_stat_sign_and_size_of_uv():  # stat · C04 (N50)
    for dUdy, c in ((2.0, 1.0), (-3.0, 0.3), (2.0, 0.0)):
        r = ch12.displaced_parcel_uv(dUdy, 0.1, 0.5, n=100_000, seed=3, correlation=c)
        exp = ch12.parcel_uv_expected(dUdy, 0.1, 0.5, c)
        assert rel(exp, -c * 0.5 * 0.1 * dUdy) < 1e-14 if c else exp == 0
        assert abs(r["uv"] - exp) < 5 * r["stderr"]
        assert rel(r["uv"], r["estimate"]) < 1e-10 or abs(r["uv"] - r["estimate"]) < 1e-14
        u, v, ell = r["samples"]
        assert rel(np.mean(u * v), r["uv"]) < 1e-10 or abs(r["uv"]) < 1e-3
        assert np.allclose(u, -ell * dUdy)
        if c:
            assert np.sign(r["uv"]) == -np.sign(dUdy)
            assert rel(r["nu_T"], -r["uv"] / dUdy) < 1e-12 and r["nu_T"] > 0
            assert rel(ch12.eddy_viscosity_from_data(r["uv"], dUdy), r["nu_T"]) < 1e-12
    r = ch12.displaced_parcel_uv(2.0, 0.1, 0.5, n=50_000, seed=1, u_extra_rms=0.3)
    assert abs(r["uv"] + 0.1) < 5 * r["stderr"] and abs(r["r_uv"]) < 0.9


def test_mean_fluxes_V1_heat_scalar_and_mixture_relations():  # V1 · C04 (12.32), (12.34), N54, N55
    Q = ch12.mean_heat_flux(np.array([0.0, 0.0, -0.5]), np.array([0.0, 0.0, 0.05]), 0.026, 1.2, 1005.0)
    assert rel(Q[2], 0.026 * 0.5 + 1.2 * 1005.0 * 0.05) < 1e-13 and Q[0] == 0
    H = ch12.turbulent_heat_flux(0.5, 0.2, 0.5, 1.2, 1005.0)
    assert rel(H, 1.2 * 1005.0 * 0.5 * 0.5 * 0.2) < 1e-14 and ch12.turbulent_heat_flux(0.5, 0.2, -0.5, 1.2, 1005.0) == -H
    assert rel(ch12.mean_heat_flux(np.zeros(3), np.array([0, 0, 0.5 * 0.5 * 0.2]), 0.026, 1.2, 1005.0)[2], H) < 1e-14
    Fy = ch12.mean_scalar_flux(np.array([-2.0, 0.0]), np.array([0.01, 0.0]), 2e-5)
    assert rel(Fy[0], 2e-5 * 2.0 + 0.01) < 1e-14
    qq = dimensional_check(lambda k, dTdz, rho, cp, wT: -k * dTdz + rho * cp * wT, "[mass]/[time]**3",
                           k=Q_(0.026, "W/(m*K)"), dTdz=Q_(-0.5, "K/m"), rho=Q_(1.2, "kg/m**3"), cp=Q_(1005.0, "J/(kg*K)"),
                           wT=Q_(0.05, "K*m/s"))
    assert rel(qq.to("W/m**2").magnitude, Q[2]) < 1e-12
    assert ch12.mixture_density(0.0, 0.7, 1.2) == 1.2 and ch12.mixture_density(1.0, 0.7, 1.2) == 0.7
    assert rel(ch12.mixture_density(0.25, 0.7, 1.2), 0.25 * 0.7 + 0.75 * 1.2) < 1e-15
    assert ch12.mass_fraction_from_volume_fraction(0.0, 0.7, 1.2) == 0 and ch12.mass_fraction_from_volume_fraction(1.0, 0.7, 1.2) == 1
    assert rel(ch12.mass_fraction_from_volume_fraction(0.3, 1.2, 1.2), 0.3) < 1e-15
    assert rel(ch12.mass_fraction_from_volume_fraction(0.25, 0.7, 1.2), 0.25 * 0.7 / (0.25 * 0.7 + 0.75 * 1.2)) < 1e-15


def test_closure_count_V1_component_formula():  # V1 · C04 (N60)
    comp = lambda m: (m + 1) * (m + 2) // 2  # noqa: E731 independent components of a symmetric tensor of order m in 3-D
    rows = ch12.closure_count(4)
    eq, unk = 4, 4
    for n, row in enumerate(rows, start=1):
        if n > 1:
            eq += comp(n)
        unk += comp(n + 1)
        assert (row["level"], row["equations"], row["unknowns"]) == (n, eq, unk) and row["closed"] is False
        assert row["unknowns"] > row["equations"]
    assert rows[0]["equations"] == 4 and rows[0]["unknowns"] == 10 and rows[1]["unknowns_with_extra"] == 32


# ======================================================================================================================
# §12.6  C05 — isotropic turbulence (12.36)–(12.45)
# ======================================================================================================================
def test_isotropic_tensor_V2_derivation_D07_form_incompressibility_and_scales():  # V2 · C05 · D07 (12.40), (12.41)
    x, y, z = sp.symbols("x y z", positive=True)
    X = (x, y, z)
    rr = sp.sqrt(x ** 2 + y ** 2 + z ** 2)
    r, L, lam, a4, u2 = sp.symbols("r L lambda a4 u2", positive=True)
    kd = lambda i, j: 1 if i == j else 0  # noqa: E731
    # step: divergence of the general isotropic form (12.40) is r_i (r F′ + 4F + G′/r) — two explicit (F, G) pairs
    for Fr, Gr in ((sp.exp(-r ** 2), 1 / (1 + r ** 2)), (r ** 2 * sp.exp(-r), sp.cos(r) * sp.exp(-r ** 2))):
        R = [[Fr.subs(r, rr) * X[i] * X[j] + Gr.subs(r, rr) * kd(i, j) for j in range(3)] for i in range(3)]
        claim = (r * sp.diff(Fr, r) + 4 * Fr + sp.diff(Gr, r) / r).subs(r, rr)
        for i in range(3):
            assert sp.simplify(sum(sp.diff(R[i][j], X[j]) for j in range(3)) - X[i] * claim) == 0
    # step: with F = u2 (f − g)/r², G = u2 g the condition r F′ + 4F + G′/r = 0 gives g = f + (r/2) f′
    f, g = sp.Function("f"), sp.Function("g")
    F = u2 * (f(r) - g(r)) / r ** 2
    G = u2 * g(r)
    cond = sp.simplify(r * sp.diff(F, r) + 4 * F + sp.diff(G, r) / r)
    gsol = sp.solve(cond, g(r))
    assert len(gsol) == 1 and z0(gsol[0] - (f(r) + r * sp.diff(f(r), r) / 2))
    # step: (12.41) is (12.40) with that g
    ri, rj, dij = sp.symbols("r_i r_j delta_ij")
    general = (F * ri * rj + G * dij).subs(g(r), gsol[0])
    assert z0(general - u2 * (f(r) * dij + r / 2 * sp.diff(f(r), r) * (dij - ri * rj / r ** 2)))
    # step: Λ_g = Λ_f/2 (integration by parts) for three decaying f
    for fr in (sp.exp(-r ** 2 / L ** 2), sp.exp(-r / L), (1 + r ** 2 / L ** 2) ** -2):
        gr = fr + r * sp.diff(fr, r) / 2
        assert z0(sp.integrate(gr, (r, 0, sp.oo)) - sp.integrate(fr, (r, 0, sp.oo)) / 2)
    # step: λ_g = λ_f/√2 from the Taylor series f = 1 − r²/λ_f² + a4 r⁴
    fs = 1 - r ** 2 / lam ** 2 + a4 * r ** 4
    gs = sp.expand(fs + r * sp.diff(fs, r) / 2)
    assert z0(-2 / sp.diff(gs, r, 2).subs(r, 0) - lam ** 2 / 2)
    # the chapter's engine agrees, and the general form with an independent g is not solenoidal
    rs = sp.Symbol("r", positive=True)       # the engine's own symbol
    for trial in (None, sp.exp(-rs), 1 / (1 + rs ** 2)):
        e = ch12.isotropic_tensor_divergence_sympy(trial) if trial is not None else ch12.isotropic_tensor_divergence_sympy()
        assert e["check"] is True and e["general_fails"] is True


def test_isotropic_dissipation_V2_derivation_D08_gradient_moments_2_4_minus1_and_30():  # V2 · C05 · D08 (12.42)–(12.43)
    x, y, z = sp.symbols("x y z")
    X = (x, y, z)
    lam, a4, u2, nu = sp.symbols("lambda a4 u2 nu", positive=True)
    r2 = x ** 2 + y ** 2 + z ** 2
    kd = lambda i, j: 1 if i == j else 0  # noqa: E731
    # (12.41) with f = 1 − r²/λ² + a4 r⁴ is a polynomial in r_i:  (r/2) f′ = −r²/λ² + 2 a4 r⁴
    fpoly = 1 - r2 / lam ** 2 + a4 * r2 ** 2
    half_rfp = -r2 / lam ** 2 + 2 * a4 * r2 ** 2
    R = [[u2 * (fpoly * kd(i, j) + half_rfp * kd(i, j) - (-1 / lam ** 2 + 2 * a4 * r2) * X[i] * X[j]) for j in range(3)]
         for i in range(3)]
    at0 = {x: 0, y: 0, z: 0}
    mom = lambda i, k, j, l: sp.expand(-sp.diff(R[i][j], X[k], X[l]).subs(at0))  # noqa: E731,E741 mean(∂u_i/∂x_k ∂u_j/∂x_l)
    unit = u2 / lam ** 2
    assert z0(mom(0, 0, 0, 0) - 2 * unit)          # mean((∂u1/∂x1)²)
    assert z0(mom(0, 1, 0, 1) - 4 * unit)          # mean((∂u1/∂x2)²)
    assert z0(mom(0, 1, 1, 0) + unit)              # mean((∂u1/∂x2)(∂u2/∂x1))
    # (12.42) summed in full over i and j (not through the book's shortcut 6{…})
    eps = nu / 2 * sum(mom(i, j, i, j) + 2 * mom(i, j, j, i) + mom(j, i, j, i) for i in range(3) for j in range(3))
    assert z0(eps - 30 * nu * unit)                                       # ε̄ = 30 ν u2/λ_f²
    assert z0(eps - 15 * nu * mom(0, 0, 0, 0))                            # = 15 ν mean((∂u1/∂x1)²)
    assert z0(eps - 6 * nu * (mom(0, 0, 0, 0) + mom(0, 1, 0, 1) + mom(0, 1, 1, 0)))     # the book's first form
    assert z0(sum(mom(i, i, j, j) for i in range(3) for j in range(3)))   # incompressible: mean((∂u_i/∂x_i)²) = 0
    e = ch12.gradient_moments_isotropic_sympy()
    assert (e["m11"], e["m12"], e["cross"], e["eps_factor"], e["eps_factor_g"], e["eps_from_m11"]) == (2, 4, -1, 30, 15, 15)
    assert e["ratio_lambda"] == sp.Rational(1, 2)


def test_isotropic_relations_V1_gaussian_f_three_routes_to_the_dissipation():  # V1 · C05 (12.39)–(12.44)
    Lc, u2, nu = 0.05, 0.36, 1.5e-5
    r = np.linspace(0.0, 12 * Lc, 24001)
    f = np.exp(-r ** 2 / Lc ** 2)
    g_exact = (1 - r ** 2 / Lc ** 2) * f
    fc = lambda s: np.exp(-s ** 2 / Lc ** 2)  # noqa: E731
    assert np.max(np.abs(ch12.transverse_from_longitudinal(r, f) - g_exact)) < 1e-5
    assert np.max(np.abs(ch12.transverse_from_longitudinal(r, fc, fprime=lambda s: -2 * s / Lc ** 2 * fc(s)) - g_exact)) < 1e-14
    assert np.max(np.abs(ch12.transverse_from_longitudinal(r, fc, dim=2, fprime=lambda s: -2 * s / Lc ** 2 * fc(s))
                         - (1 - 2 * r ** 2 / Lc ** 2) * f)) < 1e-14                  # 2-D: g = d(rf)/dr
    sc = ch12.isotropic_scales(r, f)
    assert rel(sc["Lambda_f"], math.sqrt(PI) / 2 * Lc) < 1e-7 and rel(sc["Lambda_ratio"], 0.5) < 1e-5
    assert rel(sc["lambda_f"], Lc) < 1e-5 and rel(sc["lambda_ratio"], 1 / math.sqrt(2)) < 1e-5
    e_f = ch12.dissipation_isotropic(nu, u2, lambda_f=sc["lambda_f"])
    e_g = ch12.dissipation_isotropic(nu, u2, lambda_g=sc["lambda_g"])
    fpp0 = (f[1] - 2 * f[0] + f[1]) / (r[1] - r[0]) ** 2                             # f″(0) by differences (f even)
    e_fd = -15 * nu * u2 * fpp0
    e_dx = ch12.dissipation_isotropic(nu, dudx_sq=2 * u2 / Lc ** 2)
    exact = 30 * nu * u2 / Lc ** 2
    assert max(rel(e_f, exact), rel(e_g, exact), rel(e_fd, exact)) < 2e-5 and rel(e_dx, exact) < 1e-14
    assert rel(ch12.dissipation_isotropic(nu, u2, lambda_g=Lc) / ch12.dissipation_isotropic(nu, u2, lambda_f=Lc), 0.5) < 1e-14
    with pytest.raises((ValueError, TypeError)):
        ch12.dissipation_isotropic(nu, u2, lambda_f=Lc, lambda_g=Lc)
    assert rel(ch12.taylor_reynolds_number(u2, Lc, nu), Lc * 0.6 / nu) < 1e-14
    # the tensor itself
    rv = np.array([0.03, 0.0, 0.0])
    Rt = ch12.isotropic_correlation_tensor(rv, fc, u2=u2, fprime=lambda s: -2 * s / Lc ** 2 * fc(s))
    assert rel(Rt[0, 0], u2 * fc(0.03)) < 1e-13 and rel(Rt[1, 1], u2 * (1 - 0.03 ** 2 / Lc ** 2) * fc(0.03)) < 1e-12
    assert rel(Rt[2, 2], Rt[1, 1]) < 1e-13 and abs(Rt[0, 1]) < 1e-16
    Rg = ch12.isotropic_correlation_tensor(rv, fc, g=lambda s: (1 - s ** 2 / Lc ** 2) * fc(s), u2=u2, incompressible=False)
    assert np.allclose(Rg, Rt, atol=1e-13)
    q = np.array([0.02, -0.01, 0.015])
    Rq = ch12.isotropic_correlation_tensor(q, fc, u2=u2, fprime=lambda s: -2 * s / Lc ** 2 * fc(s))
    n = q / np.linalg.norm(q)
    assert rel(n @ Rq @ n, u2 * fc(np.linalg.norm(q))) < 1e-12 and np.allclose(Rq, Rq.T)
    R0 = ch12.isotropic_correlation_tensor(np.array([1e-9, 0, 0]), fc, u2=u2, fprime=lambda s: -2 * s / Lc ** 2 * fc(s))
    assert rel(np.trace(R0), 3 * u2) < 1e-9                                           # R_ii(0) = 2ē
    div = sum((ch12.isotropic_correlation_tensor(q + 1e-5 * np.eye(3)[j], fc, u2=u2, fprime=lambda s: -2 * s / Lc ** 2 * fc(s))
               - ch12.isotropic_correlation_tensor(q - 1e-5 * np.eye(3)[j], fc, u2=u2, fprime=lambda s: -2 * s / Lc ** 2 * fc(s)))[:, j]
              for j in range(3)) / 2e-5
    assert np.max(np.abs(div)) < 1e-5 * u2 / Lc                                       # ∂R_ij/∂r_j = 0


def test_dissipation_rate_V1_strain_rotation_and_V4_spectral_identity_on_a_3d_field():  # V1/V4 · C05 (12.42)
    nu = 2e-3
    strain = np.array([[[1.0, 0, 0], [0, -1.0, 0], [0, 0, 0]]])
    rot = np.array([[[0, -1.0, 0], [1.0, 0, 0], [0, 0, 0]]])
    shear = np.array([[[0, 3.0, 0], [0, 0, 0], [0, 0, 0]]])
    assert rel(ch12.dissipation_rate(strain, nu), 4 * nu) < 1e-14 and ch12.dissipation_rate(rot, nu) == 0
    assert rel(ch12.dissipation_rate(shear, nu), nu * 9.0) < 1e-14
    L, n = 2 * PI, 32
    E = lambda K: K ** 4 * np.exp(-K ** 2 / 8.0)  # noqa: E731
    u, v, w = TS.synthetic_solenoidal_field(n, L, E, seed=1, dim=3)
    dx = L / n
    assert ch12.divergence_rms(u, v, dx, w=w) < 1e-11 * np.sqrt(np.mean(u * u)) / dx
    gs = ch12.velocity_gradient_samples(u, v, w, dx)
    assert gs.shape == (n ** 3, 3, 3) and np.max(np.abs(np.trace(gs, axis1=1, axis2=2))) < 1e-9 * np.max(np.abs(gs))
    k1 = 2 * PI * np.fft.fftfreq(n, dx)
    KZ, KY, KX = np.meshgrid(k1, k1, k1, indexing="ij")
    K2 = KX ** 2 + KY ** 2 + KZ ** 2
    spec = sum(np.abs(np.fft.fftn(c)) ** 2 for c in (u, v, w)) / n ** 6
    assert rel(ch12.dissipation_rate(gs, nu), nu * np.sum(K2 * spec)) < 1e-10        # ε̄ = ν Σ K²|û|² (solenoidal, periodic)
    # one Fourier mode: exact derivative
    xx = np.arange(n) * dx
    mode = np.broadcast_to(np.sin(3 * xx), (n, n, n))
    g1 = ch12.velocity_gradient_samples(mode, 0 * mode, 0 * mode, dx).reshape(n, n, n, 3, 3)
    assert np.max(np.abs(g1[0, 0, :, 0, 0] - 3 * np.cos(3 * xx))) < 1e-11


@slow
def test_isotropic_field_stat_ratios_two_to_one_and_g_from_f():  # stat · C05 (12.36), (12.37), (12.41), (12.43)
    L, n = 2 * PI, 32
    E = lambda K: K ** 4 * np.exp(-K ** 2 / 8.0)  # noqa: E731
    ratios, e15, gerr, shear = [], [], [], []
    for seed in range(8):
        u, v, w = TS.synthetic_solenoidal_field(n, L, E, seed=seed, dim=3)
        rep = ch12.isotropy_report(u, v, w, dx=L / n)
        ratios.append(rep["transverse_over_longitudinal"])
        shear.append(rep["shear_coefficients"]["r_12"])
        gs = ch12.velocity_gradient_samples(u, v, w, L / n)
        e15.append(15 * np.mean(rep["longitudinal"]) / ch12.dissipation_rate(gs, 1.0))
        rr, f, g = TS.longitudinal_transverse_correlation(u, v, L / n, w=w)
        assert rel(f[0], 1.0) < 1e-12 and rel(g[0], 1.0) < 1e-12
        gerr.append((g - ch12.transverse_from_longitudinal(rr, f))[1:6])
    for arr, target in ((ratios, 2.0), (e15, 1.0), (shear, 0.0)):
        se = np.std(arr, ddof=1) / math.sqrt(len(arr))
        assert abs(np.mean(arr) - target) < 5 * se, (target, np.mean(arr), se)
    gerr = np.array(gerr)
    assert np.all(np.abs(gerr.mean(axis=0)) < 5 * gerr.std(axis=0, ddof=1) / math.sqrt(8) + 0.01)  # + O(Δr²) of np.gradient
    r2d = [ch12.isotropy_report(*TS.synthetic_solenoidal_field(64, L, E, seed=s, dim=2), dx=L / 64)["transverse_over_longitudinal"]
           for s in range(12)]
    se = np.std(r2d, ddof=1) / math.sqrt(12)
    assert abs(np.mean(r2d) - 3.0) < 5 * se and abs(np.mean(r2d) - 2.0) > 5 * se        # 2-D solenoidal field: 3, not 2


# ======================================================================================================================
# §12.7  C06 — the two energy budgets (12.46), (12.47)
# ======================================================================================================================
def test_mean_energy_budget_V4_laminar_poiseuille_work_equals_dissipation():  # V4 · C06 (12.46), laminar limit
    h, nu, rho, dPdx = 0.02, 1e-5, 1.2, -3.0
    y = np.linspace(0.0, h, 801)
    U = -dPdx / (2 * rho * nu) * y * (h - y)
    b = ch12.mean_energy_budget(y, U, 0 * y, nu, dPdx=dPdx, rho=rho)
    I = b["integrals"]
    exact = dPdx ** 2 * h ** 3 / (12 * rho ** 2 * nu)                       # ∫ν U′² dy for the parabola
    assert rel(I["pressure_work"], exact) < 1e-5 and rel(-I["viscous_dissipation"], exact) < 1e-5
    assert np.max(np.abs(b["loss_to_turbulence"])) == 0
    vt = [abs(ch12.mean_energy_budget(yy, -dPdx / (2 * rho * nu) * yy * (h - yy), 0 * yy, nu, dPdx=dPdx, rho=rho)
              ["integrals"]["viscous_transport"]) for yy in (np.linspace(0.0, h, 401), y, np.linspace(0.0, h, 1601))]
    assert vt[1] < 1e-4 * exact and abs(math.log2(vt[0] / vt[2]) / 2 - 2.0) < 0.25      # transport integrates to 0 at O(Δy²)
    # the only discretisation error is that of d/dy(νUU′) (a cubic): central differences are off by (Δy²/6)·ν(UU′)‴,
    # a constant, = 2νa²Δy² with U = a y(h − y) — the residual must be exactly that
    a_ = -dPdx / (2 * rho * nu)
    assert np.max(np.abs(b["residual"][2:-2] / (2 * nu * a_ ** 2 * (y[1] - y[0]) ** 2) - 1)) < 1e-3
    b2 = ch12.mean_energy_budget(y, U, 0 * y, nu, rho=rho)                  # dP/dx inferred from the momentum balance
    assert rel(b2["integrals"]["pressure_work"], exact) < 1e-4


def test_channel_energy_budget_V4_pointwise_closure_and_integral_identity():  # V4 · C06 (12.46), (12.47) · item 6
    for Re in (180.0, 1000.0, 5200.0):
        b = ch12.channel_energy_budget(Re, KAP, 26.0)
        flat = b["pressure_work"] + b["transport"] - b["mean_dissipation"] - b["production"]
        assert np.max(np.abs(flat)) < 1e-12                                                   # (12.46) at every height
        assert np.all(b["mean_dissipation"] >= 0) and np.all(b["production"] >= -1e-15) and np.all(b["pressure_work"] >= 0)
        m, tb = b["mean"], b["turbulence"]
        assert np.max(np.abs(m["pressure_work"] + m["transport"] + m["viscous_dissipation"] + m["loss_to_turbulence"])) < 1e-12
        assert np.array_equal(m["viscous_dissipation"], -b["mean_dissipation"]) and np.array_equal(m["loss_to_turbulence"], -b["production"])
        assert np.array_equal(tb["production"], b["production"]) and np.array_equal(b["turb_sink"], b["production"])
        assert np.max(np.abs(tb["production"] + tb["dissipation_plus_transport"])) == 0       # (12.47) closes by residual
        I = b["integrals"]
        assert abs(b["identity_residual"]) < 1e-3 and rel(b["direct_fraction"], I["dissipation"] / I["work"]) < 1e-13
        assert abs(I["transport"]) < 1e-3 * I["work"]                                         # transport only moves energy
        assert rel(I["work"], b["U_bulk_plus"]) < 1e-9                                        # ∫U⁺/Re_τ dy⁺ = U_bulk⁺
        # independent quadrature of the model (my own closed forms, not the function's arrays)
        lT = lambda yp: min(KAP * yp, 0.09 * Re) * (1 - math.exp(-yp / 26.0))  # noqa: E731
        s = lambda yp: 2 * (1 - yp / Re) / (1 + math.sqrt(1 + 4 * lT(yp) ** 2 * (1 - yp / Re)))  # noqa: E731
        pts = [1.0, 5.0, 12.0, 30.0, 0.09 * Re / KAP, 0.5 * Re]
        diss = quad(lambda yp: s(yp) ** 2, 0, Re, points=pts, limit=400, epsabs=1e-12, epsrel=1e-12)[0]
        prod = quad(lambda yp: (1 - yp / Re - s(yp)) * s(yp), 0, Re, points=pts, limit=400, epsabs=1e-12, epsrel=1e-12)[0]
        Ub = quad(lambda yp: (1 - yp / Re) * s(yp), 0, Re, points=pts, limit=400, epsabs=1e-12, epsrel=1e-12)[0]
        assert rel(diss + prod, Ub) < 1e-10                                                   # the identity, exactly
        assert rel(I["dissipation"], diss) < 2e-3 and rel(I["production"], prod) < 2e-3 and rel(b["U_bulk_plus"], Ub) < 2e-3
        assert np.max(b["production"]) <= 0.25 + 1e-12 and 8 < b["yplus_peak_production"] < 14  # P ≤ τ⁺²/4 near y⁺ ≈ 10–12
    hi, lo = ch12.channel_energy_budget(5200.0, KAP, 26.0), ch12.channel_energy_budget(180.0, KAP, 26.0)
    assert hi["direct_fraction"] < lo["direct_fraction"] < 1                                  # more goes via turbulence at high Re


def test_channel_energy_budget_V3_second_order_in_the_grid():  # V3 · C06 / C12 (channel_mixing_length)
    Re = 1000.0
    lT = lambda yp: min(KAP * yp, 0.09 * Re) * (1 - math.exp(-yp / 26.0))  # noqa: E731
    s = lambda yp: 2 * (1 - yp / Re) / (1 + math.sqrt(1 + 4 * lT(yp) ** 2 * (1 - yp / Re)))  # noqa: E731
    Ucl = quad(s, 0, Re, points=[1, 5, 12, 30, 0.09 * Re / KAP], limit=400, epsabs=1e-13, epsrel=1e-13)[0]
    ns = [200, 400, 800, 1600]
    errs = [abs(ch12.channel_mixing_length(Re, KAP, n=n)["U_cl_plus"] - Ucl) for n in ns]
    assert abs(observed_order([1.0 / n for n in ns], errs) - 2.0) < 0.15 and errs[-1] < errs[0] / 50
    res = [abs(ch12.channel_energy_budget(Re, KAP, 26.0, n=n)["identity_residual"]) for n in ns]
    assert res[-1] < res[0] / 20


def test_channel_energy_budget_at_V1_closed_forms_and_agreement_with_the_arrays():  # V1 · C06 (E5 parity) · item 6
    Re = 1000.0
    b = ch12.channel_energy_budget(Re, KAP, 26.0, n=1600)
    for yp in (0.5, 5.0, 12.0, 30.0, 100.0, 400.0, 900.0):
        a = ch12.channel_energy_budget_at(yp, Re, KAP, 26.0)
        tau = 1 - yp / Re
        lT = min(KAP * yp, 0.09 * Re) * (1 - math.exp(-yp / 26.0))
        s = 2 * tau / (1 + math.sqrt(1 + 4 * lT * lT * tau))
        assert rel(a["total"], tau) < 1e-14 and rel(a["slope"], s) < 1e-12 and rel(a["lT_plus"], lT) < 1e-12
        assert abs(s + lT * lT * s * s - tau) < 1e-14                                         # the model's quadratic
        assert rel(a["uv_plus"], tau - s) < 1e-11 and rel(a["mean_dissipation"], s * s) < 1e-12
        assert rel(a["production"], (tau - s) * s) < 1e-11 and a["turb_sink"] == a["production"]
        assert a["production"] <= tau ** 2 / 4 + 1e-15
        assert abs(a["pressure_work"] + a["transport"] - a["mean_dissipation"] - a["production"]) < 1e-13
        assert abs(a["viscous_transport"] + a["reynolds_transport"] - a["transport"]) < 1e-12
        assert rel(a["dissipation_over_production"], s / (tau - s)) < 1e-9
        Uq = quad(lambda q: 2 * (1 - q / Re) / (1 + math.sqrt(1 + 4 * (min(KAP * q, 0.09 * Re) * (1 - math.exp(-q / 26.0))) ** 2
                                                              * (1 - q / Re))), 0, yp, limit=400, epsabs=1e-13, epsrel=1e-13)[0]
        assert rel(a["Uplus"], Uq) < 1e-9 and rel(a["pressure_work"], Uq / Re) < 1e-9
        for key, arr in (("Uplus", b["Uplus"]), ("production", b["production"]), ("mean_dissipation", b["mean_dissipation"]),
                         ("transport", b["transport"])):
            assert abs(np.interp(yp, b["yplus"], arr) - a[key]) < 2e-3 * max(abs(a[key]), 1e-2), (key, yp)
    eq = ch12.channel_energy_budget_at(12.0, Re, KAP, None)       # no damping: more Reynolds stress at the same height
    assert eq["uv_plus"] > ch12.channel_energy_budget_at(12.0, Re, KAP, 26.0)["uv_plus"]


def test_channel_energy_budget_V1_generic_budget_functions_agree_with_the_model_channel():  # V1 · C06
    Re = 550.0
    c = ch12.channel_mixing_length(Re, KAP, n=1600)
    y, U, uv = c["yplus"], c["Uplus"], -c["uv_plus"]
    me = ch12.mean_energy_budget(y, U, uv, 1.0, dPdx=-1.0 / Re, rho=1.0)        # wall units: ν = 1, dP/dx = −1/Re_τ
    assert np.max(np.abs(me["loss_to_turbulence"] + c["production"])) < 2e-3
    assert np.max(np.abs(me["residual"][3:-3])) < 5e-3
    tk = ch12.tke_budget(y, U, uv, eps=c["production"])
    assert np.max(np.abs(tk["production"] + me["loss_to_turbulence"])) < 1e-12   # the same term, opposite sign
    assert tk["transport_inferred"] is True and np.max(np.abs(tk["production"] + tk["dissipation"] + tk["transport"])) < 1e-12
    tk2 = ch12.tke_budget(y, U, uv, eps=0.8 * c["production"], transport=0 * y, buoyancy=0.05 * c["production"])
    assert np.max(np.abs(tk2["residual"] - 0.25 * c["production"])) < 2e-3
    st = WT.total_stress(y, U, uv, 1.0, 1.0)
    assert np.max(np.abs(st["total"] - (1 - y / Re))[2:-2]) < 5e-3 and np.allclose(st["total"], st["viscous"] + st["reynolds"])
    assert np.max(np.abs(c["total"] - (1 - y / Re))) < 1e-13 and np.allclose(c["uv_plus"] + c["dUdy_plus"], c["total"])
    assert rel(c["Re_bulk"], 2 * Re * c["U_bulk_plus"]) < 1e-12 and rel(c["Cf"], 2 / c["U_bulk_plus"] ** 2) < 1e-12
    assert rel(ch12.mean_to_turbulent_dissipation_ratio(1e4), 1e-4) < 1e-14
    assert rel(ch12.mean_to_turbulent_dissipation_ratio(1e4, urms_over_U=0.1), 1e-2) < 1e-12
    assert rel(ch12.buoyant_production(0.05, 1 / 300.0), G0 * 0.05 / 300.0) < 1e-14 and ch12.buoyant_production(-0.05, 1 / 300.0) < 0


def test_explainer_tables_V1_json_reproduces_the_functions():  # V1 · C06/C12 · item 6 (cache ≠ evidence: recompute)
    T = json.loads((REF / "explainer_tables.json").read_text(encoding="utf-8"))
    assert "ours" in T["note"] and T["kappa"] == KAP and T["A_plus"] == 26.0 and set(T["channel"]) == {"180", "550", "1000", "5200"}
    sig6 = lambda a, b: abs(a - b) <= 6e-6 * max(abs(b), 1e-300) + 1e-12  # noqa: E731
    for key, row in T["channel"].items():
        Re = float(key)
        yp = np.array(row["yplus"])
        assert yp[0] == 0 and rel(yp[-1], Re) < 1e-5 and np.all(np.diff(yp) > 0)
        for i in (1, 20, 45, 70, len(yp) - 2):
            # the stored y⁺ is itself rounded to 6 digits (the table was evaluated at the unrounded height), so each
            # stored value must lie between the function's values at y⁺(1 ∓ 6e-6), widened by its own 6-digit rounding
            lo, a, hi = (ch12.channel_energy_budget_at(float(yp[i]) * f, Re, KAP, 26.0) for f in (1 - 6e-6, 1.0, 1 + 6e-6))
            for col, src in (("Uplus", "Uplus"), ("uv_plus", "uv_plus"), ("production", "production"),
                             ("pressure_work", "pressure_work"), ("mean_dissipation", "mean_dissipation"),
                             ("turb_sink", "turb_sink"), ("transport", "transport"), ("dUdy_plus", "slope")):
                v, pad = row[col][i], 6e-6 * abs(a[src]) + 1e-12
                assert min(lo[src], hi[src]) - pad <= v <= max(lo[src], hi[src]) + pad, (key, col, i, v, a[src])
            assert sig6(row["Uplus"][i], a["Uplus"])
            assert abs(row["pressure_work"][i] + row["transport"][i] - row["mean_dissipation"][i] - row["production"][i]) < 2e-5
        assert rel(row["production_peak"], max(row["production"])) < 2e-2 and row["production_peak"] <= 0.25 + 1e-6
    vd = T["van_driest"]
    for i in (0, 30, 60, 111):
        y = vd["yplus"][i]
        assert sig6(vd["Uplus"][i], ch12.mixing_length_wall_profile(y, KAP, damping="van_driest", A_plus=26.0)["Uplus"])
        assert sig6(vd["Uplus_undamped"][i], ch12.mixing_length_wall_profile(y, KAP)["Uplus"])
    assert sig6(vd["B"], ch12.mixing_length_intercept(KAP, 26.0)) and sig6(vd["B_undamped"], (math.log(4 * KAP) - 1) / KAP)


@slow
def test_explainer_tables_V1_fresh_call_and_writer_equal_the_file(tmp_path):  # V1 · design C.3b
    on_disk = (REF / "explainer_tables.json").read_text(encoding="utf-8")
    out = ch12.write_reference_tables(dest=tmp_path)
    assert Path(out["explainer_tables"]).read_text(encoding="utf-8") == on_disk
    fresh = ch12.explainer_tables()
    assert json.loads(json.dumps(fresh)) == json.loads(on_disk)


def test_mixing_potential_energy_V1_linear_profile_closed_form():  # V1 · C06 (N82)
    z = np.linspace(0.0, 80.0, 1601)
    al, rho0, G = 2e-4, 1000.0, 0.03
    st = ch12.mixing_potential_energy_change(z, 290.0 + G * z, al, rho0)
    assert rel(st["dPE"], rho0 * G0 * al * G * 80.0 ** 3 / 12) < 1e-6 and st["stable_initially"] and st["heat_flux_sign"] == -1
    un = ch12.mixing_potential_energy_change(z, 290.0 - G * z, al, rho0)
    assert rel(un["dPE"], -st["dPE"]) < 1e-9 and not un["stable_initially"] and un["heat_flux_sign"] == 1
    assert np.allclose(st["T_final"], 290.0 + G * 40.0, rtol=1e-9)


# ======================================================================================================================
# §12.7  C07 — Kolmogorov scales and the cascade (12.48)–(12.52)
# ======================================================================================================================
def test_kolmogorov_scales_V2_exponents_by_dimensional_analysis_D11():  # V2 · C07 (12.50) · D11
    # ν [L² T⁻¹], ε̄ [L² T⁻³]: solve a·(2, −1) + b·(2, −3) = target dimension (L, T) with exact rationals
    M = sp.Matrix([[2, 2], [-1, -3]])
    sol = lambda L, T: tuple(M.solve(sp.Matrix([L, T])))  # noqa: E731
    assert sol(1, 0) == (sp.Rational(3, 4), sp.Rational(-1, 4))        # η
    assert sol(1, -1) == (sp.Rational(1, 4), sp.Rational(1, 4))        # u_K
    assert sol(0, 1) == (sp.Rational(1, 2), sp.Rational(-1, 2))        # τ_η
    e = ch12.kolmogorov_exponents()
    assert e == {"eta": (Fraction(3, 4), Fraction(-1, 4)), "u_K": (Fraction(1, 4), Fraction(1, 4)),
                 "tau_eta": (Fraction(1, 2), Fraction(-1, 2))}
    grp = DIM.solve_exponents("eta", ["nu", "eps"], {"eta": "m", "nu": "m**2/s", "eps": "m**2/s**3"})
    assert (grp["nu"], grp["eps"]) == (Fraction(-3, 4), Fraction(1, 4))   # η ν^(−3/4) ε^(1/4) is the dimensionless group
    nu, eps = Q_(1.5e-5, "m**2/s"), Q_(0.125, "m**2/s**3")
    eta = dimensional_check(lambda nu, eps: (nu ** 3 / eps) ** 0.25, "[length]", nu=nu, eps=eps)
    uK = dimensional_check(lambda nu, eps: (nu * eps) ** 0.25, "velocity", nu=nu, eps=eps)
    got = ch12.kolmogorov_scales(1.5e-5, 0.125)
    assert rel(got[0], eta.to("m").magnitude) < 1e-13 and rel(got[1], uK.to("m/s").magnitude) < 1e-13
    assert not ((nu ** 3 * eps) ** 0.25).check("[length]")               # the ν³·ε̄ mutant has the wrong units


def test_kolmogorov_scales_V1_unit_reynolds_number_and_scale_separation():  # V1 · C07 (12.48)–(12.52)
    nu = np.array([1e-6, 1.5e-5, 1e-4])
    eps = np.array([1e-4, 0.125, 30.0])
    eta, uK, tau = ch12.kolmogorov_scales(nu, eps)
    assert np.max(np.abs(eta * uK / nu - 1)) < 1e-14 and np.max(np.abs(tau * uK / eta - 1)) < 1e-14
    assert np.max(np.abs(nu / tau ** 2 / eps - 1)) < 1e-14                # ε̄ = ν/τ_η²
    assert isinstance(ch12.kolmogorov_scales(1e-6, 1.0)[0], float)
    for dU, L, v, c in ((5.0, 1000.0, 1.5e-5, 1.0), (0.3, 0.2, 1e-6, 0.7)):
        e = ch12.dissipation_outer_scaling(dU, L, c)
        assert rel(e, c * dU ** 3 / L) < 1e-15
        Re = dU * L / v
        s = ch12.scale_separation(Re, c_eps=c)
        et, uk, tt = ch12.kolmogorov_scales(v, e)
        assert rel(s["eta_over_L"], et / L) < 1e-12 and rel(s["uK_over_dU"], uk / dU) < 1e-12
        assert rel(s["tau_eta_over_T"], tt * dU / L) < 1e-12
        lam_f = math.sqrt(30 * v * dU ** 2 / e)                             # (12.43) with u2 = ΔU²
        assert rel(s["lambdaT_over_L"], lam_f / L) < 1e-12 and rel(s["R_lambda"], lam_f * dU / v) < 1e-12
        sg = ch12.scale_separation(Re, c_eps=c, microscale="g")
        assert rel(sg["lambdaT_over_L"] / s["lambdaT_over_L"], 1 / math.sqrt(2)) < 1e-12
        assert rel(ch12.inertial_range_decades(Re, c_eps=c), math.log10(L / et)) < 1e-12
        assert rel(ch12.dns_grid_points(Re, c_eps=c), (L / et) ** 3) < 1e-10
    Res = np.geomspace(1e3, 1e9, 13)
    sl = lambda key: np.polyfit(np.log(Res), np.log([ch12.scale_separation(R)[key] for R in Res]), 1)[0]  # noqa: E731
    assert abs(sl("eta_over_L") + 0.75) < 1e-12 and abs(sl("lambdaT_over_L") + 0.5) < 1e-12 and abs(sl("R_lambda") - 0.5) < 1e-12
    assert rel(ch12.dns_grid_points(1e6), 1e6 ** 2.25) < 1e-12 and rel(ch12.inertial_range_decades(1e7), 5.25) < 1e-13
    o = ch12.scale_ordering(1e4)
    assert o["ordered"] and o["eta_over_L"] < o["lambda_over_L"] < 1 and rel(o["lambda_over_L"], math.sqrt(15 / 1e4)) < 1e-12
    assert not ch12.scale_ordering(10.0)["ordered"]                         # λ_g > L below Re_L = 15
    assert rel(ch12.richardson_diffusivity(10.0, 1e-4) / ch12.richardson_diffusivity(1.0, 1e-4), 10 ** (4 / 3)) < 1e-13
    assert (Q_(1e-4, "m**2/s**3") ** Fraction(1, 3) * Q_(10.0, "m") ** Fraction(4, 3)).check("[length]**2/[time]")


def test_cascade_and_scale_table_V1_composition_of_the_scale_functions():  # V1/V7 · C07 (N87, N90)
    c = ch12.cascade_tiers(1.0, 2.0, 1e-5, ratio=2.0)
    assert np.allclose(c["velocity"] ** 3 / c["size"], c["eps"]) and rel(c["eps"], 8.0) < 1e-14
    assert np.allclose(c["Re"], c["velocity"] * c["size"] / 1e-5) and np.allclose(c["turnover"], c["size"] / c["velocity"])
    assert np.all(np.diff(c["size"]) < 0) and c["size"][-1] >= c["eta"] and 1.0 <= c["Re"][-1] <= 2 ** (4 / 3) + 1e-9
    assert rel(c["Re"][0], c["Re_L"]) < 1e-12 and c["n_tiers"] == len(c["size"])
    n1, n2 = ch12.cascade_tiers(1.0, 2.0, 1e-5)["n_tiers"], ch12.cascade_tiers(1.0, 2.0, 1e-7)["n_tiers"]
    assert abs((n2 - n1) - 0.75 * math.log(100) / math.log(2)) <= 1.0
    rows = ch12.scale_table()
    assert {r["name"] for r in rows} == {"kitchen_mixer", "wind_tunnel", "atmospheric_boundary_layer", "ocean_thermocline"}
    for r in rows:
        assert r["name"] in ch12.SCALE_CASES
        e = ch12.dissipation_outer_scaling(r["dU"], r["L"])
        et, uk, tt = ch12.kolmogorov_scales(r["nu"], e)
        assert rel(r["eps"], e) < 1e-13 and rel(r["eta"], et) < 1e-13 and rel(r["u_K"], uk) < 1e-13 and rel(r["tau_eta"], tt) < 1e-13
        assert rel(r["Re_L"], r["dU"] * r["L"] / r["nu"]) < 1e-13 and rel(r["decades"], math.log10(r["L"] / et)) < 1e-12
        assert rel(r["lambda_g"], math.sqrt(15 * r["nu"] * r["dU"] ** 2 / e)) < 1e-12
        assert rel(r["R_lambda"], ch12.taylor_reynolds_number(r["dU"] ** 2, r["lambda_g"], r["nu"])) < 1e-12
        assert rel(r["grid_points"], (r["L"] / et) ** 3) < 1e-10
    abl = ch12.scale_table("atmospheric_boundary_layer")
    assert 1e-4 < abl["eta"] < 1e-3                                        # a fraction of a millimetre
    assert rel(ch12.scale_table({"L": 2.0, "dU": 3.0, "nu": 1e-5})["eps"], 13.5) < 1e-13


# ======================================================================================================================
# §12.7  C08 — Kolmogorov's spectrum (12.53)–(12.55)
# ======================================================================================================================
def test_inertial_spectrum_V2_derivation_D12_one_group_and_the_18_over_55():  # V2 · C08 · D12 (12.53)–(12.55)
    a, b = sp.symbols("a b")
    # S11 [L³ T⁻²] = ε̄^a [L^{2a} T^{−3a}] · k1^b [L^{−b}]
    sol = sp.solve([sp.Eq(2 * a - b, 3), sp.Eq(-3 * a, -2)], [a, b])
    assert sol[a] == sp.Rational(2, 3) and sol[b] == sp.Rational(-5, 3)          # −5/3, not the printed +5/3
    # with ν kept there are two groups: S11/(ν^{5/4} ε̄^{1/4}) and k1 ν^{3/4} ε̄^{−1/4}  (12.53)
    p, q = sp.symbols("p q")
    s2 = sp.solve([sp.Eq(2 * p + 2 * q, 3), sp.Eq(-p - 3 * q, -2)], [p, q])
    assert s2[p] == sp.Rational(5, 4) and s2[q] == sp.Rational(1, 4)
    # isotropic link E11(k1) = ∫_{k1}^∞ (E/K)(1 − k1²/K²) dK for E = C ε̄^{2/3} K^{−5/3}
    K, k1, C, eps = sp.symbols("K k1 C epsilon", positive=True)
    E11 = sp.integrate(C * eps ** sp.Rational(2, 3) * K ** sp.Rational(-5, 3) / K * (1 - k1 ** 2 / K ** 2), (K, k1, sp.oo))
    assert z0(E11 / (C * eps ** sp.Rational(2, 3) * k1 ** sp.Rational(-5, 3)) - sp.Rational(18, 55))
    kc = ch12.kolmogorov_constants(1.5)
    assert rel(kc["C1_one_sided"], 18 * 1.5 / 55) < 1e-14 and rel(kc["C1_two_sided"], 9 * 1.5 / 55) < 1e-14 and kc["C"] == 1.5
    assert ch12.KOLMOGOROV_C == 1.5
    num = ch12.one_dimensional_from_3d(7.0, lambda Kk: ch12.inertial_spectrum_3d(Kk, 0.4))
    assert rel(num, 18 / 55 * 1.5 * 0.4 ** (2 / 3) * 7.0 ** (-5 / 3)) < 1e-8
    assert rel(ch12.one_dimensional_from_3d(7.0, lambda Kk: ch12.inertial_spectrum_3d(Kk, 0.4), one_sided=False), num / 2) < 1e-12
    assert rel(ch12.inertial_spectrum_1d(7.0, 0.4, two_sided=False), num) < 1e-8          # (12.54), one-sided
    assert rel(ch12.inertial_spectrum_1d(7.0, 0.4), num / 2) < 1e-8                       # (12.55): two-sided = half
    assert rel(ch12.inertial_spectrum_1d(7.0, 0.4, C1=0.5, two_sided=False), 0.5 * 0.4 ** (2 / 3) * 7.0 ** (-5 / 3)) < 1e-14


def test_inertial_spectrum_V5_one_dimensional_constant_against_sreenivasan():  # V5 · C08 (N95)
    bm = json.loads((REF / "benchmarks.json").read_text(encoding="utf-8"))["kolmogorov_C1_one_sided"]
    c1 = ch12.kolmogorov_constants()["C1_one_sided"]
    assert abs(c1 - bm["value"]) < bm["std"]            # 0.491 against 0.53 ± 0.055 (one standard deviation of the survey)
    assert abs(c1 / bm["value"] - 1) < 0.10             # 7 %: inside the survey's scatter, outside its 95 % interval of the mean


def test_model_spectrum_V4_dissipation_and_energy_integrals():  # V4 · C08 (N97) · item 8 (c_L)
    eps = 0.7
    for nu in (1e-4, 1e-6):
        eta = (nu ** 3 / eps) ** 0.25
        d = quad(lambda lk: 2 * nu * math.exp(3 * lk) * ch12.model_spectrum(math.exp(lk), eps, nu), math.log(1e-9),
                 math.log(60 / eta), limit=400, epsabs=0, epsrel=1e-11)[0]
        assert rel(d, eps) < 1e-8                                                  # 2ν∫K²E dK = ε̄ exactly (Pao, no L)
    nu, L = 1e-9, 1.0
    eta = (nu ** 3 / eps) ** 0.25
    Ek = lambda lk, **kw: ch12.model_spectrum(math.exp(lk), eps, nu, L=L, **kw) * math.exp(lk)  # noqa: E731
    k = quad(Ek, math.log(1e-7), math.log(60 / eta), limit=600, epsabs=0, epsrel=1e-10)[0]
    assert rel(k, (eps * L) ** (2 / 3)) < 1e-3              # c_L = 6.78 is the value for which L = ē^{3/2}/ε̄ (C = 1.5, p0 = 2)
    k5 = quad(lambda lk: Ek(lk, c_L=5.0), math.log(1e-7), math.log(60 / eta), limit=600, epsabs=0, epsrel=1e-10)[0]
    assert rel(k5, (eps * L) ** (2 / 3)) > 0.05             # another c_L does not have that property
    d = quad(lambda lk: 2 * nu * math.exp(3 * lk) * ch12.model_spectrum(math.exp(lk), eps, nu, L=L), math.log(1e-7),
             math.log(60 / eta), limit=600, epsabs=0, epsrel=1e-10)[0]
    assert rel(d, eps) < 1e-2                               # with the energy range: within 1 % when L/η ≳ 10³
    Kb = np.geomspace(1e3, 1e5, 40)                         # well inside the inertial range (L/η ≈ 6·10⁶)
    slope, const = ch12.fit_inertial_range(Kb, ch12.model_spectrum(Kb, eps, nu, L=L), (1e3, 1e5), eps=eps)
    assert abs(slope + 5 / 3) < 0.02 and abs(const - 1.5) < 0.02
    assert ch12.model_spectrum(50.0, 1.0, 1e-6, L=1.0, kind="pope") == ch12.model_spectrum(50.0, 1.0, 1e-6, L=1.0, kind="pao")
    lowK = np.array([1e-3, 2e-3])
    El = ch12.model_spectrum(lowK, eps, nu, L=L)
    assert abs(math.log(El[1] / El[0]) / math.log(2.0) - 2.0) < 1e-3     # E ∝ K^{p0} at the largest scales


def test_inertial_spectrum_V1_fit_and_kolmogorov_normalisation_collapse():  # V1 · C08 (12.53), (12.54)
    K = np.geomspace(2.0, 500.0, 60)
    slope, a = ch12.fit_inertial_range(K, ch12.inertial_spectrum_3d(K, 0.3, C=1.5), (2.0, 500.0))
    assert abs(slope + 5 / 3) < 1e-12 and rel(a, 1.5 * 0.3 ** (2 / 3)) < 1e-11
    s2, c2 = ch12.fit_inertial_range(K, ch12.inertial_spectrum_3d(K, 0.3, C=1.5), (5.0, 100.0), eps=0.3)
    assert abs(s2 + 5 / 3) < 1e-12 and rel(c2, 1.5) < 1e-12
    hand = np.polyfit(np.log(K), np.log(ch12.inertial_spectrum_1d(K, 0.3)), 1)        # the notebook's from-scratch fit
    assert abs(hand[0] + 5 / 3) < 1e-12 and rel(math.exp(hand[1]), 9 * 1.5 / 55 * 0.3 ** (2 / 3)) < 1e-11
    curves = []
    for nu, eps in ((1e-6, 1e-3), (1.5e-5, 2.0), (1e-4, 50.0)):
        eta = (nu ** 3 / eps) ** 0.25
        kk = np.array([1e-3, 1e-2, 1e-1]) / eta
        x, Phi = ch12.kolmogorov_normalize_spectrum(kk, ch12.inertial_spectrum_1d(kk, eps), nu, eps)
        assert np.allclose(x, [1e-3, 1e-2, 1e-1])
        curves.append(Phi)
    assert np.allclose(curves[0], curves[1], rtol=1e-12) and np.allclose(curves[0], curves[2], rtol=1e-12)
    assert np.allclose(curves[0], 9 * 1.5 / 55 * np.array([1e-3, 1e-2, 1e-1]) ** (-5 / 3), rtol=1e-12)   # Φ = C1 (k1η)^(−5/3)
    S = dimensional_check(lambda eps, k1: 0.25 * eps ** Fraction(2, 3) * k1 ** Fraction(-5, 3), "[length]**3/[time]**2",
                          eps=Q_(0.3, "m**2/s**3"), k1=Q_(10.0, "1/m"))
    assert S.magnitude > 0


# ======================================================================================================================
# §12.8  C09 — free shear flows: the plane jet (12.56)–(12.75), Table 12.1 exponents
# ======================================================================================================================
JET_KW = dict(C5="from_invariant", xi_half=0.1)       # ξ½ = 0.10: a labelled illustrative half-width (design convention 10)


def test_plane_jet_V2_derivation_D13_thin_layer_flux_form_and_invariant():  # V2 · C09 · D13 (12.61)–(12.62)
    x, y = sp.symbols("x y")
    U, V, uv = sp.Function("U")(x, y), sp.Function("V")(x, y), sp.Function("uv")(x, y)
    cont = sp.diff(U, x) + sp.diff(V, y)                                           # (12.58)
    thin = U * sp.diff(U, x) + V * sp.diff(U, y) + sp.diff(uv, y)                  # (12.61), everything on the left
    flux = sp.diff(U ** 2, x) + sp.diff(U * V + uv, y)                             # ∂(U²)/∂x + ∂(UV + mean(uv))/∂y
    assert z0(thin + U * cont - flux)        # so ∫ dy gives d/dx ∫U² dy = −[UV + mean(uv)]_{−∞}^{+∞} = 0  ⇒ (12.62)
    # the cross-stream equation (12.60) reduced to ∂(P + ρ mean(v²))/∂y = 0 is what removes ∂P/∂x (stated; no algebra)


def test_plane_jet_V2_derivation_D14_similarity_equation():  # V2 · C09 · D14 (12.63) (★★★)
    x, y, xi = sp.symbols("x y xi", positive=True)
    delta, Ucl, Psi = sp.Function("delta")(x), sp.Function("U_CL")(x), sp.Function("Psi")(x)
    c1 = delta * sp.diff(Ucl, x) / Ucl
    c2 = c1 + sp.diff(delta, x)
    c3 = Psi / Ucl ** 2
    # two explicit profile pairs: I(ξ) = ∫_0^ξ F, F = I′, and an odd G
    for I_xi, G_xi in ((sp.tanh(xi), xi * sp.exp(-xi ** 2)), (xi / (1 + xi ** 2), sp.sin(xi) * sp.exp(-xi ** 2))):
        F_xi = sp.diff(I_xi, xi)
        at = lambda e: e.subs(xi, y / delta)  # noqa: E731
        U = Ucl * at(F_xi)                                                         # (12.56)
        V = -sp.diff(Ucl * delta * at(I_xi), x)                                    # V = −∫_0^y ∂U/∂x dy, V(0) = 0
        assert sp.simplify(sp.diff(V, y) + sp.diff(U, x)) == 0                     # step: continuity (12.58)
        assert sp.simplify(V.subs(y, 0)) == 0
        # step: ∂U/∂x = U′_CL F − U_CL F′ ξ δ′/δ (the similarity variable moves with x)
        assert sp.simplify(sp.diff(U, x) - (sp.diff(Ucl, x) * at(F_xi)
                                            - Ucl * at(sp.diff(F_xi, xi)) * (y / delta) * sp.diff(delta, x) / delta)) == 0
        # step: V = −U′_CL δ I − U_CL δ′ (I − ξF)
        assert sp.simplify(V - (-sp.diff(Ucl, x) * delta * at(I_xi)
                                - Ucl * sp.diff(delta, x) * (at(I_xi) - y / delta * at(F_xi)))) == 0
        residual = U * sp.diff(U, x) + V * sp.diff(U, y) - sp.diff(Psi * at(G_xi), y)        # (12.61) with −mean(uv) = ΨG
        book = c1 * F_xi ** 2 - c2 * sp.diff(F_xi, xi) * I_xi - c3 * sp.diff(G_xi, xi)        # (12.63), left − right
        assert sp.simplify((residual * delta / Ucl ** 2).subs(y, xi * delta) - book) == 0
        wrong = c1 * F_xi ** 2 + c2 * sp.diff(F_xi, xi) * I_xi - c3 * sp.diff(G_xi, xi)       # sign of the middle term
        assert sp.simplify((residual * delta / Ucl ** 2).subs(y, xi * delta) - wrong) != 0
    e = ch12.plane_jet_similarity_sympy()
    assert e["check"] is True and e["stress_integral_check"] is True
    byname = lambda ex: ex.subs({s_: sp.Symbol(s_.name) for s_ in ex.free_symbols})  # noqa: E731 (engine's own symbols)
    m, n = sp.symbols("m n")
    assert e["exponential_family"][1] == 0 and z0(byname(e["momentum_flux_exponent"]) - (m + 2 * n))
    pf = [byname(sp.sympify(c)) for c in e["power_family"]]
    assert z0(sp.simplify(pf[0] / pf[2]) - n) and z0(sp.simplify(pf[1] / pf[2]) - (m + n))


def test_plane_jet_V2_derivation_D15_exponents_stress_profile_and_volume_flux():  # V2 · C09 · D15 (12.64)–(12.68) · item 8
    x, C4, rho, Js, I2, I1 = sp.symbols("x C4 rho J_s I2 I1", positive=True)
    gam = sp.Symbol("gamma", real=True)
    xi = sp.Symbol("xi", real=True)
    Ucl = C4 * x ** gam                                              # first of (12.64) with δ = x
    J = rho * Ucl ** 2 * x * I2                                      # (12.65)
    g_sol = sp.solve(sp.diff(sp.log(J), x), gam)
    assert g_sol == [sp.Rational(-1, 2)]                             # 2γ + 1 = 0
    c1 = sp.simplify(x * sp.diff(Ucl, x) / Ucl).subs(gam, g_sol[0])
    c2 = c1 + 1
    assert (c1, c2) == (sp.Rational(-1, 2), sp.Rational(1, 2))
    C4_sol = sp.solve(sp.Eq(J.subs(gam, g_sol[0]), Js), C4)[0]
    assert z0(C4_sol - sp.sqrt(Js / (rho * I2)))                     # C5 = C4 (ρ/J_s)^{1/2} = I2^{−1/2}
    # (12.63) with those values integrates once: C3 G = −½ F ∫_0^ξ F  (any F = I′)
    Ifun = sp.Function("I")(xi)
    F = sp.diff(Ifun, xi)
    lhs = c1 * F ** 2 - c2 * sp.diff(F, xi) * Ifun
    assert z0(lhs - sp.diff(-sp.Rational(1, 2) * F * Ifun, xi))      # ⇒ C3 G′ = d/dξ(−½ F I), G(0) = 0 since I(0) = 0
    assert not z0(lhs - sp.diff(+sp.Rational(1, 2) * F * Ifun, xi))  # the sign of analysis row 121 does not satisfy (12.63)
    Vdot = Ucl.subs(gam, g_sol[0]) * x * I1                          # (12.68): U_CL δ ∫F dξ
    assert sp.simplify(sp.diff(sp.log(Vdot), x) * x) == sp.Rational(1, 2)
    # the coded profile: Gaussian F, closed-form integral
    xs = np.linspace(-0.4, 0.4, 81)
    a = math.log(2) / 0.1 ** 2
    Fg = np.exp(-a * xs ** 2)
    Ig = 0.5 * math.sqrt(PI / a) * erf(math.sqrt(a) * xs)
    assert np.max(np.abs(ch12.plane_jet_stress_profile(xs, xi_half=0.1) - (-0.5 * Fg * Ig))) < 1e-14
    assert np.max(np.abs(ch12.plane_jet_stress_profile(xs, F=lambda s: np.exp(-a * s * s)) - (-0.5 * Fg * Ig))) < 1e-10
    assert np.allclose(ch12.plane_jet_stress_profile(xs, xi_half=0.1, C3=4.0) * 4.0, -0.5 * Fg * Ig, atol=1e-14)
    assert np.all(ch12.plane_jet_stress_profile(xs[xs > 0], xi_half=0.1) < 0)       # −mean(uv) < 0 where dU/dy < 0


def test_plane_jet_stress_V1_sign_from_the_momentum_equation_by_quadrature():  # V1 · C09 (12.61), (12.67) · item 8
    """Independent of every similarity formula: integrate U U_x + V U_y = −∂mean(uv)/∂y outward from the axis."""
    rho, Js, x0 = 1.2, 3.0, 2.0
    y = np.linspace(0.0, 0.8, 4001)
    h = 1e-4
    Uf = lambda xx: ch12.plane_jet_mean_velocity(xx, y, Js, rho, **JET_KW)  # noqa: E731
    U = Uf(x0)
    Ux = (Uf(x0 + h) - Uf(x0 - h)) / (2 * h)
    V = -np.concatenate([[0.0], np.cumsum(0.5 * (Ux[1:] + Ux[:-1]) * np.diff(y))])   # continuity, V(0) = 0
    Uy = np.gradient(U, y, edge_order=2)
    rhs = U * Ux + V * Uy
    minus_uv = np.concatenate([[0.0], np.cumsum(0.5 * (rhs[1:] + rhs[:-1]) * np.diff(y))])   # −mean(uv)(y), zero on the axis
    coded = ch12.plane_jet_reynolds_stress(x0, y, Js, rho, **JET_KW)
    scale = np.max(np.abs(coded))
    assert np.max(np.abs(minus_uv - coded)) < 2e-5 * scale
    assert np.all(coded[1:2000] < 0) and scale > 1e-3                       # negative above the axis: same sign as dU/dy
    assert np.max(np.abs(ch12.plane_jet_cross_velocity(x0, y, Js, rho, **JET_KW) - V)) < 1e-6 * np.max(np.abs(V))
    nuT = ch12.eddy_viscosity_from_data(-coded[5:1500], Uy[5:1500])
    assert np.all(nuT > 0)                                                  # an eddy viscosity of the right sign exists
    assert np.allclose(ch12.plane_jet_reynolds_stress(x0, -y, Js, rho, **JET_KW), -coded, atol=1e-15)   # odd in y
    assert np.allclose(ch12.plane_jet_reynolds_stress(x0, y, Js, rho, C3=7.0, **JET_KW), coded, rtol=1e-12)
    assert rel(ch12.plane_jet_reynolds_stress(2 * x0, 2 * 0.1, Js, rho, **JET_KW) / ch12.plane_jet_reynolds_stress(x0, 0.1, Js, rho, **JET_KW), 0.5) < 1e-12
    t = ch12.thin_shear_layer_terms(x0, y[::40], Js, rho, 1.5e-5, **JET_KW)
    big = np.max(np.abs(t["stress_gradient"]))
    assert np.max(np.abs(t["residual"])) < 1e-6 * big and t["viscous_over_stress"] < 1e-3 and t["ratio_V_over_U"] < 0.2
    assert np.max(np.abs(t["advection_x"] + t["advection_y"] - t["stress_gradient"])) < 1e-6 * big


def test_plane_jet_V4_momentum_flux_invariant_and_growing_volume_flux():  # V4 · C09 (12.62), (12.66), (12.68)
    rho, Js = 1.2, 3.0
    Jn, Vn, Vf = [], [], []
    xs = (0.5, 1.0, 2.0, 8.0)
    for x in xs:
        y = np.linspace(-1.5 * x, 1.5 * x, 6001)
        U = ch12.plane_jet_mean_velocity(x, y, Js, rho, **JET_KW)
        Jn.append(ch12.jet_momentum_flux_per_span(y, U, rho))
        assert rel(Jn[-1], JET.jet_momentum_flux(y, U, rho)) < 1e-14 and rel(Jn[-1], rho * np.trapezoid(U * U, y)) < 1e-13
        Vn.append(np.trapezoid(U, y))
        Vf.append(ch12.plane_jet_volume_flux(x, Js, rho, **JET_KW))
        assert rel(ch12.plane_jet_centerline_velocity(x, Js, rho, **JET_KW), U[3000]) < 1e-13
        assert rel(np.interp(0.1 * x, y, U), 0.5 * U[3000]) < 1e-6          # half-width y½ = ξ½ x: the jet is a wedge
    assert np.max(np.abs(np.array(Jn) / Js - 1)) < 1e-10                    # the invariant, equal to J_s
    assert np.max(np.abs(np.array(Vn) / np.array(Vf) - 1)) < 1e-10
    assert rel(Vf[3] / Vf[1], math.sqrt(8.0)) < 1e-13 and rel(Vf[2] / Vf[0], 2.0) < 1e-13   # ∝ x^{1/2}: entrainment
    I = ch12.profile_integrals(0.1, 0.15)
    a = math.log(2) / 0.01
    assert rel(I["I1"], quad(lambda s: math.exp(-a * s * s), -np.inf, np.inf)[0]) < 1e-12
    assert rel(I["I2"], quad(lambda s: math.exp(-2 * a * s * s), -np.inf, np.inf)[0]) < 1e-12 and rel(I["I2"], I["I1"] / math.sqrt(2)) < 1e-14
    assert rel(I["IHF"], quad(lambda s: math.exp(-a * s * s - math.log(2) / 0.15 ** 2 * s * s), -np.inf, np.inf)[0]) < 1e-12
    assert ch12.profile_integrals(0.1)["IHF"] is None
    assert rel(ch12.plane_jet_centerline_velocity(1.0, Js, rho, **JET_KW), math.sqrt(Js / rho / I["I2"])) < 1e-13   # C5 = I2^(−1/2)
    assert rel(ch12.plane_jet_centerline_velocity(4.0, Js, rho, C5=2.0, xi_half=0.1), 2.0 * math.sqrt(Js / rho) / 2.0) < 1e-14
    assert rel(ch12.plane_jet_mean_velocity(4.5, 0.0, Js, rho, x0=0.5, **JET_KW), ch12.plane_jet_mean_velocity(4.0, 0.0, Js, rho, **JET_KW)) < 1e-14
    # entrainment: v_e = ½ dV̇/dx = −V(+∞)
    ve = ch12.plane_jet_entrainment_velocity(2.0, Js, rho, **JET_KW)
    dV = (ch12.plane_jet_volume_flux(2.0 + 1e-5, Js, rho, **JET_KW) - ch12.plane_jet_volume_flux(2.0 - 1e-5, Js, rho, **JET_KW)) / 2e-5
    assert rel(ve, 0.5 * dV) < 1e-8 and rel(ve, -ch12.plane_jet_cross_velocity(2.0, 50.0, Js, rho, **JET_KW)) < 1e-10
    assert rel(ve / ch12.plane_jet_centerline_velocity(2.0, Js, rho, **JET_KW), I["I1"] / 4) < 1e-12
    assert ch12.plane_jet_cross_velocity(2.0, 0.02, Js, rho, **JET_KW) > 0 > ch12.plane_jet_cross_velocity(2.0, 0.6, Js, rho, **JET_KW)
    # the wrong-exponent demonstration (E6) and the fit of a virtual origin
    assert ch12.wrong_exponent_fluxes(4.0, -0.5, 1.0) == (1.0, 2.0)
    mom, vol = ch12.wrong_exponent_fluxes(4.0, -0.4, 1.0)
    assert rel(mom, 4.0 ** 0.2) < 1e-14 and rel(vol, 4.0 ** 0.6) < 1e-14
    S, x0 = ch12.virtual_origin_fit(np.array([1.0, 2.0, 3.0, 5.0]), 0.2 * (np.array([1.0, 2.0, 3.0, 5.0]) - 0.07))
    assert rel(S, 0.2) < 1e-12 and rel(x0, 0.07) < 1e-10
    assert rel(ch12.slot_momentum_flux(1.1, 20.0, 0.01), 1.1 * 400 * 0.01) < 1e-15 and rel(ch12.slot_mass_flux(1.1, 20.0, 0.01), 0.22) < 1e-15
    assert (Q_(1.1, "kg/m**3") * Q_(20.0, "m/s") ** 2 * Q_(0.01, "m")).check("[force]/[length]")
    assert ch12.gaussian_profile(0.0, 0.1) == 1.0 and rel(ch12.gaussian_profile(0.1, 0.1), 0.5) < 1e-15
    assert ch12.gaussian_profile(-0.07, 0.1) == ch12.gaussian_profile(0.07, 0.1)


def test_plane_jet_scalar_V4_mass_flux_invariant_and_centerline_decay():  # V4 · C09 (12.69)–(12.71)
    rho, Js, Ms = 1.2, 3.0, 0.4
    Y1 = ch12.plane_jet_mass_fraction(1.0, 0.0, Ms, Js, rho, C6="from_invariant", xi_half_Y=0.15, **JET_KW)
    Y4 = ch12.plane_jet_mass_fraction(4.0, 0.0, Ms, Js, rho, C6="from_invariant", xi_half_Y=0.15, **JET_KW)
    assert rel(Y1 / Y4, 2.0) < 1e-13                                         # Y_CL ∝ x^{−1/2}
    I = ch12.profile_integrals(0.1, 0.15)
    C5 = I["I2"] ** -0.5
    assert rel(Y1, Ms / math.sqrt(rho * Js) / (C5 * I["IHF"])) < 1e-12      # C6 = 1/(C5 ∫HF dξ)
    Yc = ch12.plane_jet_mass_fraction(4.0, 0.0, Ms, Js, rho, C6=1.7, xi_half_Y=0.15)
    assert rel(Yc, 1.7 * Ms / math.sqrt(rho * Js) / 2.0) < 1e-14            # (12.71) with a given constant
    assert rel(ch12.plane_jet_mass_fraction(4.0, 0.60, Ms, Js, rho, C6=1.7, xi_half_Y=0.15), 0.5 * Yc) < 1e-13
    assert (Q_(Ms, "kg/(m*s)") / (Q_(rho, "kg/m**3") * Q_(Js, "N/m")) ** 0.5 * Q_(1.0, "m") ** -0.5).check("[]")


def test_free_shear_exponents_V2_derivation_D16_each_flow_from_its_invariant():  # V2 · C09 · D16 (Table 12.1 exponents)
    m, n, s = sp.symbols("m n s")
    systems = {
        "plane_jet": ([sp.Eq(m, 1), sp.Eq(2 * n + m, 0), sp.Eq(n + s + m, 0)], True),           # U²δ; U Y δ
        "round_jet": ([sp.Eq(m, 1), sp.Eq(2 * n + 2 * m, 0), sp.Eq(n + s + 2 * m, 0)], True),     # U²δ²; U Y δ²
        "plane_wake": ([sp.Eq(n + m, 0), sp.Eq(m - 1, n), sp.Eq(s + m, 0)], True),               # U∞ΔUδ; dδ/dx ~ ΔU/U∞; U∞Yδ
        "round_wake": ([sp.Eq(n + 2 * m, 0), sp.Eq(m - 1, n), sp.Eq(s + 2 * m, 0)], True),
        "plane_plume": ([sp.Eq(m, 1), sp.Eq(n + s + m, 0), sp.Eq(2 * n + m - 1, s + m)], True),   # buoyancy flux; d(U²δ)/dx ~ Yδ
        "round_plume": ([sp.Eq(m, 1), sp.Eq(n + s + 2 * m, 0), sp.Eq(2 * n + 2 * m - 1, s + 2 * m)], True),
        "shear_layer": ([sp.Eq(m, 1), sp.Eq(n, 0)], False),
    }
    expected_reynolds = {"plane_jet": "1/2", "round_jet": "0", "plane_wake": "0", "round_wake": "-1/3",
                         "plane_plume": "1", "round_plume": "2/3", "shear_layer": "1"}
    assert set(systems) == set(ch12.FREE_SHEAR_FLOWS)
    for flow, (eqs, has_s) in systems.items():
        sol = sp.solve(eqs, [m, n, s] if has_s else [m, n], dict=True)[0]
        ex = ch12.free_shear_exponents(flow, return_equations=True)
        assert isinstance(ex["width"], Fraction) and len(ex["equations"]) == 2 and ex["invariant"]
        assert sp.Rational(ex["width"].numerator, ex["width"].denominator) == sol[m], flow
        assert sp.Rational(ex["velocity"].numerator, ex["velocity"].denominator) == sol[n], flow
        if has_s:
            assert sp.Rational(ex["scalar"].numerator, ex["scalar"].denominator) == sol[s], flow
        else:
            assert ex["scalar"] is None
        assert ex["reynolds"] == ex["width"] + ex["velocity"] == ch12.local_reynolds_number_exponent(flow)
        assert ch12.local_reynolds_number_exponent(flow) == Fraction(expected_reynolds[flow])
    assert "equations" not in ch12.free_shear_exponents("plane_jet")


def test_free_shear_flows_V4_invariants_do_not_vary_with_x():  # V4/V7 · C09 (N122–N124), illustrative constants only
    assert ch12.FREE_SHEAR_CONSTANTS == {}                    # no default set: no book table in a public file
    cst = {"C_U": 2.0, "C_Y": 1.5, "xi_half_U": 0.1, "xi_half_Y": 0.15}       # labelled illustrative
    com = dict(d=0.01, U0=10.0, rho_s=1.0, rho=1.2, Y0=1.0)
    inv = {k: [] for k in ("plane_jet", "round_jet", "plane_wake", "round_wake", "plane_plume", "round_plume")}
    for x in (1.0, 2.0, 4.0):
        c = np.linspace(0.0, 1.5 * x, 6001)
        pj = ch12.free_shear_flow("plane_jet", x, c, constants=cst, **com)
        inv["plane_jet"].append((np.trapezoid(pj["U"] ** 2, c), np.trapezoid(pj["U"] * pj["Y"], c)))
        rj = ch12.free_shear_flow("round_jet", x, c, constants=cst, **com)
        inv["round_jet"].append((np.trapezoid(rj["U"] ** 2 * c, c), np.trapezoid(rj["U"] * rj["Y"] * c, c)))
        pp = ch12.free_shear_flow("plane_plume", x, c, constants=cst, **com)
        inv["plane_plume"].append((np.trapezoid(pp["U"] * pp["Y"], c),))
        rp = ch12.free_shear_flow("round_plume", x, c, constants=cst, **com)
        inv["round_plume"].append((np.trapezoid(rp["U"] * rp["Y"] * c, c),))
        cw = np.linspace(0.0, 12.0 * 0.1 * math.sqrt(0.01 * x), 6001)             # 12 half-widths of the plane wake
        pw = ch12.free_shear_flow("plane_wake", x, cw, constants=cst, theta=0.01, U_inf=10.0)
        inv["plane_wake"].append((np.trapezoid(10.0 - pw["U"], cw),))
        cr = np.linspace(0.0, 12.0 * 0.1 * (1e-4 * x) ** (1 / 3), 6001)            # 12 half-widths of the round wake
        rw = ch12.free_shear_flow("round_wake", x, cr, constants=cst, theta=0.01, U_inf=10.0)
        inv["round_wake"].append((np.trapezoid((10.0 - rw["U"]) * cr, cr),))
        for res, flow in ((pj, "plane_jet"), (rj, "round_jet"), (pp, "plane_plume"), (rp, "round_plume")):
            assert rel(res["width"], 0.1 * x) < 1e-13 and res["exponents"]["flow"] == flow
            assert rel(np.interp(res["width"], c, res["U"]), 0.5 * res["centreline"]) < 1e-6
        assert rel(pw["width"], 0.1 * math.sqrt(0.01 * x)) < 1e-13 and rel(rw["width"], 0.1 * (1e-4 * x) ** (1 / 3)) < 1e-13
    for flow, vals in inv.items():
        v = np.array(vals)
        assert np.max(np.abs(v / v[0] - 1)) < 1e-9, flow
    for flow in inv:
        kw = dict(theta=0.01, U_inf=10.0) if "wake" in flow else com
        c1 = ch12.free_shear_centerline(flow, 1.0, constants=cst, **kw)
        c8 = ch12.free_shear_centerline(flow, 8.0, constants=cst, **kw)
        ex = ch12.free_shear_exponents(flow)
        assert abs(math.log(c8["U_CL"] / c1["U_CL"]) / math.log(8.0) - float(ex["velocity"])) < 1e-12, flow
        if c1["Y_CL"] is not None:
            assert abs(math.log(c8["Y_CL"] / c1["Y_CL"]) / math.log(8.0) - float(ex["scalar"])) < 1e-12, flow
    assert rel(ch12.free_shear_centerline("plane_jet", 0.04, constants=cst, **com)["U_CL"], 2.0 * 10.0 * math.sqrt(1 / 1.2) / 2.0) < 1e-13
    sl = ch12.free_shear_flow("shear_layer", 2.0, np.array([-10.0, 0.0, 10.0]), constants={"dxi80_coeff": 0.1}, U1=10.0, U2=5.0)
    assert np.allclose(sl["U"], [5.0, 7.5, 10.0]) and rel(sl["dxi80"], 0.1 * 5.0 / 7.5) < 1e-14
    edge = ch12.free_shear_flow("shear_layer", 2.0, np.array([-1.0, 1.0]) * sl["width"], constants={"dxi80_coeff": 0.1}, U1=10.0, U2=5.0)
    assert np.allclose(edge["U"], [5.5, 9.5], rtol=1e-12)                    # the central 80 % of the velocity difference
    with pytest.raises((KeyError, ValueError, TypeError)):
        ch12.free_shear_flow("plane_jet", 1.0, 0.0, constants={}, **com)
    assert ch12.free_shear_profile("round_jet", 0.0, xi_half=0.1) == 1.0 and rel(ch12.free_shear_profile("plane_wake", 0.3, xi_half=0.3), 0.5) < 1e-15
    assert rel(ch12.free_shear_profile("shear_layer", 0.0, xi_half=0.1), 0.5) < 1e-15
    assert ch12.free_shear_profile("shear_layer", -5.0, xi_half=0.1) < 1e-12 and rel(ch12.free_shear_profile("shear_layer", 5.0, xi_half=0.1), 1.0) < 1e-12
    xr = ch12.round_jet_distance_for_mass_fraction(0.05, 0.01, 1.0, 1.2, 1.0, 1.5)
    assert rel(ch12.free_shear_centerline("round_jet", xr, constants=cst, **com)["Y_CL"], 0.05) < 1e-13
    st = ch12.stoichiometric_mass_fraction(2.016, 28.96, 0.5, 0.21)          # hydrogen in air, public molar masses
    assert rel(st["v_fuel"] + st["v_oxidiser"], 1.0) < 1e-15 and rel(st["v_fuel"], 0.42 / 1.42) < 1e-14
    assert rel(st["Y_fuel"], 0.42 * 2.016 / (0.42 * 2.016 + 28.96)) < 1e-13 and rel(st["M_mixture"], (0.42 * 2.016 + 28.96) / 1.42) < 1e-13


def test_plane_jet_closures_V2_sech2_profile_and_qualitative_energy_budget():  # V2 / qualitative · C09 (N128, N130)
    xi, a = sp.symbols("xi a", positive=True)
    F = sp.sech(a * xi) ** 2
    Iint = sp.tanh(a * xi) / a
    nu_hat = 1 / (4 * a ** 2)
    assert sp.simplify((-F * Iint / 2 - nu_hat * sp.diff(F, xi)).rewrite(sp.exp)) == 0     # −½F∫F = ν̂F′
    xs = np.linspace(0.0, 0.3, 31)
    p = ch12.plane_jet_eddy_viscosity_profile(xs, 0.1)
    assert rel(p["a"], math.acosh(math.sqrt(2)) / 0.1) < 1e-14 and rel(p["nu_hat"], 1 / (4 * p["a"] ** 2)) < 1e-14
    assert np.allclose(p["F"], 1 / np.cosh(p["a"] * xs) ** 2) and rel(np.interp(0.1, xs, p["F"]), 0.5) < 1e-12
    assert p["F"][-1] > ch12.gaussian_profile(0.3, 0.1)                   # fatter tails than the Gaussian fit
    b = ch12.jet_tke_budget(np.linspace(0.0, 0.3, 61), 0.1)
    assert np.max(np.abs(b["advection"] + b["production"] + b["dissipation"] + b["transport"])) < 1e-14
    assert abs(b["production"][0]) < 1e-15 and np.all(b["production"] >= 0) and np.all(b["dissipation"] <= 0)
    assert 0.3 < np.linspace(0.0, 0.3, 61)[np.argmax(b["production"])] / 0.1 < 1.5      # peaks near the largest shear
    assert "qualitative" in str(b.get("label", "qualitative"))


# ======================================================================================================================
# §12.9  C10 — wall units, the linear stress, the viscous sublayer (12.76)–(12.82), (12.90)–(12.91)
# ======================================================================================================================
def test_wall_units_V2_pi_groups_and_dimensions_D18():  # V2 · C10 · D18 (12.79)–(12.81)
    g = WT.law_of_the_wall_groups()
    by = {next(k for k in d if d[k] == 1 and k in ("U", "y")): d for d in g}
    assert by["U"] == {"U": 1, "rho": Fraction(1, 2), "tau0": Fraction(-1, 2), "nu": 0}          # U/u_*
    assert by["y"] == {"y": 1, "rho": Fraction(-1, 2), "tau0": Fraction(1, 2), "nu": -1}         # y u_*/ν
    d = WT.defect_law_groups()
    byd = {next(k for k in e if e[k] == 1 and k in ("U", "y")): e for e in d}
    assert byd["U"]["delta"] == 0 and byd["y"] == {"y": 1, "rho": 0, "tau0": 0, "delta": -1}      # (12.83): U/u_*, y/δ
    tau0, rho, nu, y, U = Q_(0.3, "Pa"), Q_(1.2, "kg/m**3"), Q_(1.5e-5, "m**2/s"), Q_(2e-3, "m"), Q_(4.0, "m/s")
    us = dimensional_check(lambda tau0, rho: (tau0 / rho) ** 0.5, "velocity", tau0=tau0, rho=rho)
    assert (y * us / nu).check("[]") and (U / us).check("[]") and not (tau0 / rho).check("[length]/[time]")
    yp, Up, ust, lnu = WT.wall_units(2e-3, 4.0, 0.3, 1.2, 1.5e-5)
    assert rel(ust, us.to("m/s").magnitude) < 1e-14 and rel(ust, WT.friction_velocity(0.3, 1.2)) < 1e-15
    assert rel(lnu, 1.5e-5 / ust) < 1e-15 and rel(WT.viscous_length(1.5e-5, ust), lnu) < 1e-15
    assert rel(yp, 2e-3 * ust / 1.5e-5) < 1e-14 and rel(Up, 4.0 / ust) < 1e-14
    yb, Ub = WT.from_wall_units(yp, Up, ust, 1.5e-5)
    assert rel(yb, 2e-3) < 1e-14 and rel(Ub, 4.0) < 1e-14
    assert rel(WT.friction_reynolds_number(0.5, 0.03, 1.5e-5), 1000.0) < 1e-13
    assert WT.viscous_sublayer(3.2) == 3.2 and np.allclose(WT.viscous_sublayer(np.array([0.0, 1.0])), [0.0, 1.0])


def test_channel_stress_V2_derivation_D17_linear_total_stress_and_pressure_gradients():  # V2 · C10 · D17 (12.76)–(12.77), (12.90)–(12.91)
    y, h, tau0, d, G, c0 = sp.symbols("y h tau0 d G c0")
    tau = c0 + G * y                                         # dτ̄/dy = dP/dx = G, the same at every y (12.76)–(12.77)
    sol = sp.solve([sp.Eq(tau.subs(y, 0), tau0), sp.Eq(tau.subs(y, h), -tau0)], [c0, G])   # symmetry: τ̄(h) = −τ̄(0)
    assert z0(sol[G] + 2 * tau0 / h) and z0(tau.subs(sol) - tau0 * (1 - 2 * y / h))          # (12.90)
    Gp = sp.solve(sp.Eq(-G * sp.pi * d ** 2 / 4, tau0 * sp.pi * d), G)[0]                    # plug: pressure force = wall force
    assert z0(Gp + 4 * tau0 / d)                                                             # (12.91)
    yy = np.linspace(0.0, 0.1, 51)
    t = WT.channel_total_stress(yy, 0.1, 0.3)
    assert np.allclose(t, 0.3 * (1 - 2 * yy / 0.1)) and t[0] == 0.3 and abs(t[25]) < 1e-16 and rel(t[-1], -0.3) < 1e-14
    dP = WT.channel_pressure_gradient(0.3, 0.1)
    assert rel(dP, -6.0) < 1e-14 and np.max(np.abs(WT.channel_momentum_residual(yy, t, dP))) < 1e-10
    assert np.max(np.abs(WT.channel_momentum_residual(yy, 0.3 * (1 - yy / 0.1), dP))) > 2.9      # the half-height mutant
    assert rel(WT.pipe_pressure_gradient(0.3, 0.1), -12.0) < 1e-14
    assert rel(WT.wall_stress_from_pressure_gradient(-12.0, 0.1), 0.3) < 1e-14
    assert rel(abs(LAM.pipe_wall_stress(0.05, dpdz=-12.0)), WT.wall_stress_from_pressure_gradient(-12.0, 0.1)) < 1e-13   # ch08 (8.8)
    assert (Q_(0.3, "Pa") / Q_(0.1, "m")).check("[pressure]/[length]")


@needs_dns
def test_viscous_sublayer_V5_lee_moser_dns_follows_u_plus_equals_y_plus():  # V5 · C10 (12.82) · item 7
    for case in (180, 550, 1000, 2000, 5200):
        d = dns(case)
        m = (d["yplus"] > 1e-6) & (d["yplus"] < 1.0)                 # the first row is the wall itself
        assert m.sum() >= 3
        yv, Uv = d["yplus"][m], d["Uplus"][m]
        # U⁺ = y⁺ (12.82) with the first correction of the linear stress (12.77), dU⁺/dy⁺ = 1 − y⁺/Re_τ:
        assert np.max(np.abs(Uv / WT.viscous_sublayer(yv) - 1)) < 1 / (2 * d["Re_tau"]) + 5e-4        # within y⁺/(2Re_τ) + 0.05 %
        assert np.max(np.abs(Uv - (yv - yv ** 2 / (2 * d["Re_tau"]))) / yv) < 5e-4
        m5 = (d["yplus"] > 1e-6) & (d["yplus"] <= 5.0)
        dev5 = np.max(np.abs(d["Uplus"][m5] / d["yplus"][m5] - 1))
        assert 0.01 < dev5 < 0.08                               # a few per cent by y⁺ = 5: why the sublayer "ends" there
        assert abs(d["dUdy_plus"][0] - 1.0) < 1e-6              # wall slope: τ0 = μ dU/dy
        # the total stress is linear (12.77): viscous part dU⁺/dy⁺ never exceeds 1 − y⁺/Re_τ
        assert np.all(d["dUdy_plus"] <= 1 - d["yplus"] / d["Re_tau"] + 1e-3)
    d = dns(5200)
    sp_ = WT.stress_partition(d["yplus"][1:], d["Re_tau"], kappa=KAP, B=BLOG)
    assert np.max(np.abs(sp_["total"] - (1 - d["yplus"][1:] / d["Re_tau"]))) < 1e-13
    assert np.max(np.abs(sp_["viscous"] + sp_["reynolds"] - sp_["total"])) < 1e-13 and np.all(sp_["reynolds"] >= 0)
    inner = d["y_over_delta"][1:] < 0.2
    assert np.max(np.abs(sp_["viscous"] - d["dUdy_plus"][1:])[inner]) < 0.08      # model partition: "approximate" (measured 0.073)
    w = WT.stress_partition(0.0, 1000.0, kappa=KAP, B=BLOG)
    assert w["viscous"] == 1.0 and w["total"] == 1.0 and w["reynolds"] == 0.0
    from scipy.optimize import brentq
    for Re in (1e3, 1e5):
        yeq = brentq(lambda q: WT.stress_partition(q, Re, kappa=KAP, B=BLOG)["viscous"]
                     - WT.stress_partition(q, Re, kappa=KAP, B=BLOG)["reynolds"], 3.0, 30.0)
        assert 9.0 < yeq < 11.0                                  # equal shares in the buffer layer (9.9 and 9.8 here)


def test_layer_name_V7_boundaries_and_both_sides():  # V7 · C10/C11 (N139, N149)
    assert WT.layer_name(4.999, 0.001) == "viscous sublayer" and WT.layer_name(5.0, 0.001) == "buffer layer"
    assert WT.layer_name(29.999, 0.01) == "buffer layer" and WT.layer_name(30.0, 0.01) == "logarithmic layer"
    assert WT.layer_name(500.0, 0.15) == "logarithmic layer" and WT.layer_name(500.0, 0.1501) == "wake region"
    assert WT.layer_name(3.0, 0.5) == "viscous sublayer" and WT.layer_name(12.0, 0.5) == "buffer layer"    # y⁺ decides first
    assert list(WT.layer_name(np.array([1.0, 10.0, 100.0, 4000.0]), np.array([2e-4, 2e-3, 2e-2, 0.8]))) == [
        "viscous sublayer", "buffer layer", "logarithmic layer", "wake region"]
    assert WT.layer_name(40.0, 40.0 / 180.0) == "wake region"                 # δ⁺ < 200: no logarithmic layer left
    assert WT.layer_name(8.0, 0.001, sublayer=10.0) == "viscous sublayer"


def test_boundary_layer_stress_V1_blasius_profile_returns_the_viscous_stress():  # V1/V3 · C10 (12.78)
    Uinf, nu, rho = 2.0, 1e-5, 1.2

    def run(nx, ny):
        x = np.linspace(0.4, 0.6, nx)
        y = np.linspace(0.0, 0.02, ny)
        eta = y[:, None] * np.sqrt(Uinf / (nu * x[None, :]))
        f, fp, fpp = BL.blasius_profile(eta)
        U = Uinf * fp
        V = 0.5 * np.sqrt(nu * Uinf / x[None, :]) * (eta * fp - f)
        got = WT.boundary_layer_stress_from_profile(x, y, U, V, 0.0, rho)
        exact = rho * nu * Uinf * np.sqrt(Uinf / (nu * x[None, :])) * (fpp - 0.332057336)       # μ ∂U/∂y − τ0
        j = nx // 2
        return np.max(np.abs(got[:, j] - exact[:, j])), abs(got[-1, j]), abs(exact[0, j] - 0.0), rho * nu * Uinf * math.sqrt(Uinf / (nu * x[j])) * 0.332057336

    e1, edge1, _, tau0 = run(21, 201)
    e2, edge2, _, _ = run(41, 401)
    assert e2 < 2e-3 * tau0 and 3.0 < e1 / e2 < 5.0                  # second order in the grid
    assert rel(edge2, tau0) < 2e-3                                   # at the edge τ̄ → 0: the integral returns τ0


# ======================================================================================================================
# §12.9  C11 — the logarithmic law (12.83)–(12.93)
# ======================================================================================================================
def test_log_law_V2_derivation_D19_overlap_matching():  # V2 · C11 · D19 (12.84)–(12.89)
    y, ustar, nu, delta, kap, A, B, Uinf = sp.symbols("y u_* nu delta kappa A B U_inf", positive=True)
    yp, xi = sp.symbols("y_plus xi", positive=True)
    f, Fd = sp.Function("f"), sp.Function("F")
    U_in = ustar * f(y * ustar / nu)                                 # (12.80)
    U_out = Uinf - ustar * Fd(y / delta)                             # (12.84)
    dU_in = sp.diff(U_in, y).doit()
    dU_out = sp.diff(U_out, y).doit()
    # (12.85), (12.86): the two expressions for dU/dy; multiplied by y/u_* each depends on one variable only (12.87)
    lhs = sp.simplify((dU_in * y / ustar).subs(y, yp * nu / ustar))
    rhs = sp.simplify((dU_out * y / ustar).subs(y, xi * delta))
    assert lhs.free_symbols <= {yp} and rhs.free_symbols <= {xi}
    assert z0(lhs - yp * sp.diff(f(yp), yp)) and z0(rhs + xi * sp.diff(Fd(xi), xi))
    f_sol = sp.dsolve(sp.Eq(yp * sp.diff(f(yp), yp), 1 / kap), f(yp)).rhs          # a function of y⁺ = a function of ξ = 1/κ
    F_sol = sp.dsolve(sp.Eq(-xi * sp.diff(Fd(xi), xi), 1 / kap), Fd(xi)).rhs
    assert z0(sp.diff(f_sol, yp) - 1 / (kap * yp)) and z0(sp.diff(F_sol, xi) + 1 / (kap * xi))
    log_in = sp.log(yp) / kap + B                                    # (12.88)
    log_out = -sp.log(xi) / kap + A                                  # (12.89)
    total = sp.expand_log((log_in + log_out).subs(xi, yp / sp.Symbol("delta_plus", positive=True)), force=True)
    assert z0(total - (sp.log(sp.Symbol("delta_plus", positive=True)) / kap + A + B))   # friction law: y eliminated
    e = ch12.overlap_matching_sympy()
    assert e["checks"] == {"gradients_match": True, "f_satisfies": True, "F_satisfies": True, "friction_law": True}


def test_log_law_V1_slope_defect_form_friction_law_and_crossing():  # V1 · C11 (12.87)–(12.89)
    yp = np.geomspace(30.0, 3e4, 200)
    Up = WT.log_law(yp, kappa=KAP, B=BLOG)
    assert np.max(np.abs(WT.log_law_indicator(yp, Up) * KAP - 1)) < 1e-4            # y⁺ dU⁺/dy⁺ = 1/κ (natural log)
    assert rel(WT.log_law(math.e, kappa=KAP, B=BLOG) - WT.log_law(1.0, kappa=KAP, B=BLOG), 1 / KAP) < 1e-13
    assert np.max(np.abs(WT.log_law_indicator(yp, yp) / yp - 1)) < 1e-3             # on U⁺ = y⁺ the indicator is y⁺
    xi = np.array([0.01, 0.1, 1.0])
    Fd = WT.log_law_defect(xi, kappa=KAP, A=1.0)
    assert np.allclose(Fd, -np.log(xi) / KAP + 1.0) and Fd[-1] == 1.0
    for Re in (1e3, 1e5):
        fl = WT.friction_law_from_overlap(Re, kappa=KAP, A=1.0, B=BLOG)
        for y_ in (40.0, 0.1 * Re):
            assert rel(WT.log_law(y_, kappa=KAP, B=BLOG) + WT.log_law_defect(y_ / Re, kappa=KAP, A=1.0), fl) < 1e-13
    yc = WT.log_law_crossing(kappa=KAP, B=BLOG)
    assert abs(yc - WT.log_law(yc, kappa=KAP, B=BLOG)) < 1e-9 and 10.7 < yc < 10.9
    assert WT.log_law_crossing(kappa=KAP, B=5.5) > yc and WT.layer_name(yc, 0.001) == "buffer layer"
    xi_, de = WT.velocity_defect(np.array([0.01, 0.05]), np.array([15.0, 18.0]), 20.0, 0.8, 0.1)
    assert np.allclose(xi_, [0.1, 0.5]) and np.allclose(de, [6.25, 2.5])
    c = WT.LOG_LAW_CONSTANTS
    assert c["classical"]["kappa"] == 0.41 and c["classical"]["B"] == 5.0 and all("citation" in v for v in c.values())
    with pytest.raises(TypeError):
        WT.log_law(100.0)                                                           # κ and B have no silent defaults


def test_fit_log_law_V1_exact_recovery_and_the_bias_of_a_composite_profile():  # V1 · C11 (12.88) · flagged item 1
    yp = np.geomspace(1.0, 5000.0, 400)
    fit = WT.fit_log_law(yp, WT.log_law(yp, kappa=0.39, B=4.3), window=(30.0, 0.15), Re_tau=5000.0)
    assert abs(fit["kappa"] - 0.39) < 1e-10 and abs(fit["B"] - 4.3) < 1e-9 and fit["rms_residual"] < 1e-12
    assert fit["yplus_min"] >= 30.0 and fit["yplus_max"] <= 750.0 and fit["n_points"] > 100
    m = (yp >= 30.0) & (yp <= 750.0)
    slope, icpt = np.polyfit(np.log(yp[m]), WT.log_law(yp[m], kappa=0.39, B=4.3), 1)     # the notebook's from-scratch fit
    assert abs(1 / slope - fit["kappa"]) < 1e-12 and abs(icpt - fit["B"]) < 1e-10
    assert WT.fit_log_law(yp, WT.log_law(yp, kappa=0.39, B=4.3), window=(100.0, 400.0))["yplus_max"] <= 400.0
    # Spalding's curve approaches the log law from BELOW (the defect at y⁺ = 30 is ≈ 0.66), so a fit that starts at
    # y⁺ = 30 is biased: κ and B come out low.  The function is right; the window must start further out.
    assert -0.75 < WT.spalding_uplus(30.0, kappa=KAP, B=BLOG) - WT.log_law(30.0, kappa=KAP, B=BLOG) < -0.55
    assert abs(WT.spalding_uplus(1000.0, kappa=KAP, B=BLOG) - WT.log_law(1000.0, kappa=KAP, B=BLOG)) < 1e-3
    comp = WT.composite_profile_plus(yp, 5000.0, kappa=KAP, B=BLOG, Pi=0.0)
    f30 = WT.fit_log_law(yp, comp, window=(30.0, 0.15), Re_tau=5000.0)
    assert f30["kappa"] < 0.39 and f30["B"] < 4.2                                     # the design's (0.41 ± 0.01, 5.0 ± 0.2) is not met
    ypL = np.geomspace(1.0, 1e5, 2000)
    compL = WT.composite_profile_plus(ypL, 1e5, kappa=KAP, B=BLOG, Pi=0.0)
    kap_seq = [WT.fit_log_law(ypL, compL, window=(lo, 0.15), Re_tau=1e5) for lo in (30.0, 100.0, 300.0, 1000.0)]
    err = [abs(k["kappa"] - KAP) for k in kap_seq]
    assert err[0] > err[1] > err[2] > err[3] and err[2] < 0.002 and abs(kap_seq[2]["B"] - BLOG) < 0.06
    assert abs(kap_seq[3]["kappa"] - KAP) < 1e-4 and abs(kap_seq[3]["B"] - BLOG) < 2e-3


@needs_dns
def test_log_law_V5_lee_moser_von_karman_constant_from_the_dns_profile():  # V5 · C11 (12.87), (12.88) · items 1, 7
    bm = json.loads((REF / "benchmarks.json").read_text(encoding="utf-8"))["lee_moser_kappa"]
    d = dns(5200)
    assert abs(d["Re_tau"] - 5185.897) < 1e-3
    for lo in (3 * math.sqrt(d["Re_tau"]), 350.0):               # two usual starts of the logarithmic region
        fit = WT.fit_log_law(d["yplus"], d["Uplus"], window=(lo, 0.15), Re_tau=d["Re_tau"])
        assert abs(fit["kappa"] - bm["value"]) < bm["uncertainty"], (lo, fit)      # κ = 0.384 ± 0.004
        assert fit["n_points"] >= 20 and fit["rms_residual"] < 5e-3
    band = (d["yplus"] > 350.0) & (d["y_over_delta"] < 0.15)
    assert abs(1 / np.mean(d["yplus"][band] * d["dUdy_plus"][band]) - bm["value"]) < bm["uncertainty"]   # the authors' own dU/dy
    ind = WT.log_law_indicator(d["yplus"][1:], d["Uplus"][1:])
    assert abs(1 / np.mean(ind[band[1:]]) - bm["value"]) < bm["uncertainty"]                             # our indicator function
    f30 = WT.fit_log_law(d["yplus"], d["Uplus"], window=(30.0, 0.15), Re_tau=d["Re_tau"])
    assert f30["kappa"] > bm["value"] + bm["uncertainty"]       # starting at y⁺ = 30 is too early for real data as well
    fit = WT.fit_log_law(d["yplus"], d["Uplus"], window=(350.0, 0.15), Re_tau=d["Re_tau"])
    logm = (d["yplus"] > 350.0) & (d["y_over_delta"] < 0.15)
    assert np.max(np.abs(WT.log_law(d["yplus"][logm], kappa=fit["kappa"], B=fit["B"]) - d["Uplus"][logm])) < 0.01
    # the classical textbook pair is a different, coarser description of the same data: within 3 % over 30 < y⁺ < 0.15 δ⁺
    wide = (d["yplus"] > 30.0) & (d["y_over_delta"] < 0.15)
    assert np.max(np.abs(WT.log_law(d["yplus"][wide], kappa=KAP, B=BLOG) / d["Uplus"][wide] - 1)) < 0.03


@needs_dns
def test_spalding_profile_V5_against_lee_moser_inner_layer_approximate():  # V5 (approximate) · C11 (N153) · item 7
    worst = 0.0
    for case in (1000, 2000, 5200):
        d = dns(case)
        inner = (d["yplus"] > 0) & (d["y_over_delta"] < 0.15)
        dev = WT.spalding_uplus(d["yplus"][inner], kappa=KAP, B=BLOG) - d["Uplus"][inner]
        worst = max(worst, float(np.max(np.abs(dev))))
        sub = d["yplus"][inner] < 3.0
        assert np.max(np.abs(dev[sub])) < 0.03
    assert worst < 1.0                       # within one wall unit everywhere (measured 0.77, in the buffer layer): "approximate"


def test_spalding_and_composite_profiles_V1_limits_inverse_and_slope():  # V1 · C11 (N153, N154)
    U = np.linspace(0.0, 26.0, 131)
    y = WT.spalding_yplus(U, kappa=KAP, B=BLOG)
    kU = KAP * U
    assert np.allclose(y, U + math.exp(-KAP * BLOG) * (np.exp(kU) - 1 - kU - kU ** 2 / 2 - kU ** 3 / 6), rtol=1e-10, atol=1e-13)
    small = 1e-3
    assert rel(WT.spalding_yplus(small, kappa=KAP, B=BLOG) - small, math.exp(-KAP * BLOG) * (KAP * small) ** 4 / 24) < 1e-3
    back = WT.spalding_uplus(y[1:], kappa=KAP, B=BLOG)
    assert np.max(np.abs(back - U[1:])) < 1e-8 and WT.spalding_uplus(0.0, kappa=KAP, B=BLOG) == 0.0
    assert rel(WT.spalding_uplus(0.01, kappa=KAP, B=BLOG), 0.01) < 1e-8                         # → y⁺ (12.82)
    # → log law (12.88), slowly: rearranging the formula, U⁺ − [ln(y⁺)/κ + B] = ln(1 + c/y⁺)/κ exactly, with
    # c = −U⁺ + e^{−κB}[1 + κU⁺ + (κU⁺)²/2 + (κU⁺)³/6] = O(U⁺³)
    for ybig in (1e3, 1e4, 1e6):
        Ub_ = WT.spalding_uplus(ybig, kappa=KAP, B=BLOG)
        kU_ = KAP * Ub_
        cc = -Ub_ + math.exp(-KAP * BLOG) * (1 + kU_ + kU_ ** 2 / 2 + kU_ ** 3 / 6)
        assert abs(Ub_ - WT.log_law(ybig, kappa=KAP, B=BLOG) - math.log1p(cc / ybig) / KAP) < 1e-9
    assert abs(WT.spalding_uplus(1e6, kappa=KAP, B=BLOG) - WT.log_law(1e6, kappa=KAP, B=BLOG)) < 1e-3
    Um = 0.5 * (U[1:] + U[:-1])
    num = np.diff(U) / np.diff(y)
    assert np.max(np.abs(WT.spalding_slope(Um, kappa=KAP, B=BLOG) / num - 1)) < 2e-3
    assert WT.spalding_slope(0.0, kappa=KAP, B=BLOG) == 1.0
    Ubig = WT.spalding_uplus(5e3, kappa=KAP, B=BLOG)
    assert rel(WT.spalding_slope(Ubig, kappa=KAP, B=BLOG), 1 / (KAP * 5e3)) < 2e-3              # → 1/(κy⁺)
    xi = np.linspace(0.0, 1.0, 101)
    for kind in ("cubic", "sin2"):
        W = WT.coles_wake(xi, kind=kind)
        assert W[0] == 0 and rel(W[-1], 1.0) < 1e-15 and abs(W[1] - W[0]) < 1e-3 and abs(W[-1] - W[-2]) < 1e-3
        assert np.all(np.diff(W) >= 0) and rel(WT.coles_wake(0.5, kind=kind), 0.5) < 1e-14
    assert np.max(np.abs(WT.coles_wake(xi) - WT.coles_wake(xi, kind="sin2"))) < 0.02
    assert np.allclose(WT.coles_wake(xi), 3 * xi ** 2 - 2 * xi ** 3) and WT.coles_wake(1.7) == 1.0
    yp = np.geomspace(0.5, 4000.0, 60)
    assert np.allclose(WT.composite_profile_plus(yp, 4000.0, kappa=KAP, B=BLOG, Pi=0.0), WT.spalding_uplus(yp, kappa=KAP, B=BLOG))
    cp = WT.composite_profile_plus(yp, 4000.0, kappa=KAP, B=BLOG, Pi=0.5)
    assert np.allclose(cp - WT.spalding_uplus(yp, kappa=KAP, B=BLOG), 2 * 0.5 / KAP * WT.coles_wake(yp / 4000.0))
    assert abs(cp[-1] - (WT.spalding_uplus(4000.0, kappa=KAP, B=BLOG) + 2 * 0.5 / KAP)) < 1e-12  # edge value
    assert abs(cp[-1] - (math.log(4000.0) / KAP + BLOG + 2 * 0.5 / KAP)) < 0.01                  # ≈ log law + 2Π/κ
    assert WT.composite_profile_plus(9000.0, 4000.0, kappa=KAP, B=BLOG, Pi=0.5) == cp[-1]       # held beyond the edge
    us, nu, dl = 0.9, 1.5e-5, 4000.0 * 1.5e-5 / 0.9
    Udim = WT.composite_profile(yp * nu / us, dl, us, nu, kappa=KAP, B=BLOG, Pi=0.5)
    assert np.allclose(Udim, us * cp, rtol=1e-10)


def test_rough_wall_V1_derivation_D20_log_law_with_roughness_length():  # V1/V2 · C11 · D20 (12.93)
    yp, y0p, kap = sp.symbols("y_plus y0_plus kappa", positive=True)
    Bsym = sp.solve(sp.Eq(sp.log(y0p) / kap + sp.Symbol("B"), 0), sp.Symbol("B"))[0]           # U = 0 at y = y0
    assert z0(sp.expand_log(sp.log(yp) / kap + Bsym - sp.log(yp / y0p) / kap, force=True)) and z0(Bsym + sp.log(y0p) / kap)
    us, y0, nu = 0.4, 0.03, 1.5e-5
    z = np.array([0.03, 1.0, 10.0, 100.0])
    U = WT.rough_wall_log_law(z, us, y0, kappa=0.4)
    assert U[0] == 0 and np.allclose(U, us / 0.4 * np.log(z / y0))
    assert np.allclose(U, us * WT.log_law(z * us / nu, kappa=0.4, B=-math.log(y0 * us / nu) / 0.4), atol=1e-12)
    assert math.isnan(WT.rough_wall_log_law(0.01, us, y0, kappa=0.4))
    assert rel(WT.friction_velocity_from_wind(U[2], 10.0, y0, kappa=0.4), us) < 1e-14
    CD = WT.drag_coefficient_neutral(10.0, y0, kappa=0.4)
    assert rel(CD, (us / U[2]) ** 2) < 1e-13 and WT.drag_coefficient_neutral(10.0, 1e-4, kappa=0.4) < CD
    assert rel(WT.rough_wall_log_law(10.0, us, y0, kappa=0.4, ) , WT.rough_wall_log_law(10.0, us, y0, kappa=0.4)) == 0
    # the wind does not depend on viscosity at all (fully rough): same U for any ν — nothing to pass
    assert "nu" not in inspect.signature(WT.rough_wall_log_law).parameters


def test_wall_correlations_V1_nagib_chauhan_zpg_fits_and_pipe_friction():  # V1/V7 · C11 (12.92), N150, N151, N160
    for B in (-3.0, 1.0, 4.2, 5.0, 11.0):
        k = WT.nagib_chauhan_kappa(B)
        assert rel(k * B, 1.6 * (math.exp(0.1663 * B) - 1)) < 1e-13                   # (12.92)
        assert rel(WT.nagib_chauhan_B(k), B) < 1e-8
    assert rel(WT.nagib_chauhan_kappa(1e-9), 1.6 * 0.1663) < 1e-6 and math.isnan(WT.nagib_chauhan_kappa(13.0))
    assert WT.nagib_chauhan_kappa(5.0) > WT.nagib_chauhan_kappa(4.0)
    with pytest.raises(ValueError):
        WT.nagib_chauhan_B(0.05)
    prev = None
    for Rex in (1e6, 1e7, 1e8, 1e9):
        b = WT.zpg_boundary_layer(Rex * 1.5e-5 / 30.0, 30.0, 1.5e-5, kappa=0.384)
        x = Rex * 1.5e-5 / 30.0
        assert rel(b["Re_x"], Rex) < 1e-12 and rel(b["theta"], 0.016 * x * Rex ** -0.15) < 1e-13
        assert rel(b["delta_star"], b["theta"] * math.exp(7.11 * 0.384 / math.log(b["Re_theta"]))) < 1e-13
        br = math.log(b["Re_delta_star"]) / 0.384 + 3.30
        assert rel(b["delta99"], 0.2 * b["delta_star"] * br) < 1e-13 and rel(b["Cf"], 2.0 / br ** 2) < 1e-13
        assert rel(b["Cf"], 2 / b["U_inf_plus"] ** 2) < 1e-13 and rel(b["u_star"], 30.0 / b["U_inf_plus"]) < 1e-13   # C_f = 2/(U∞⁺)²
        assert b["theta"] < b["delta_star"] < b["delta99"] and 1.2 < b["H"] < 1.5 and rel(b["H"], b["delta_star"] / b["theta"]) < 1e-14
        assert rel(b["Re_tau"], b["delta99"] * b["u_star"] / 1.5e-5) < 1e-12 and b["valid"] == (Rex > 1e6)
        cf = [WT.skin_friction_zpg(Rex, law=law, kappa=0.384) for law in ("monkewitz", "schultz_grunow", "white")]
        assert rel(cf[0], b["Cf"]) < 1e-13 and min(cf) > 0.664 / math.sqrt(Rex)
        # the three fits agree within 10 % up to Re_x = 10⁸; at 10⁹ White's and Schultz-Grunow's differ by 12 %
        assert (max(cf) / min(cf) < 1.10) if Rex <= 1e8 else (1.10 < max(cf) / min(cf) < 1.13)
        assert rel(cf[1], 0.370 * math.log10(Rex) ** -2.584) < 1e-13 and rel(cf[2], 0.455 / math.log(0.06 * Rex) ** 2) < 1e-13
        assert prev is None or cf[0] < prev
        prev = cf[0]
    assert rel(WT.skin_friction_zpg(1e7, law="power_fifth"), 0.8 * 0.074 * 1e7 ** -0.2) < 1e-13
    # pipe: area mean of the log law (integration by parts) and the implicit friction law
    a, us, nu = 0.05, 0.3, 1e-6
    quadr = quad(lambda y: 2 / a ** 2 * us * (math.log(y * us / nu) / KAP + BLOG) * (a - y), 0, a, epsabs=0, epsrel=1e-13)[0]
    assert rel(WT.pipe_bulk_velocity_loglaw(a, us, nu, kappa=KAP, B=BLOG), quadr) < 1e-10
    assert rel(WT.pipe_bulk_velocity_loglaw(a, us, nu, kappa=KAP, B=BLOG), us * (math.log(a * us / nu) / KAP + BLOG - 1.5 / KAP)) < 1e-14
    prev_f = None
    slope = math.log(10) / (KAP * math.sqrt(8))                 # 1/√f = slope·log10(Re√f) + icpt from the log law
    icpt = (BLOG - 1.5 / KAP - math.log(2 * math.sqrt(8)) / KAP) / math.sqrt(8)
    for Red in (1e4, 1e5, 1e6, 1e7):
        fP = WT.pipe_friction_factor_turbulent(Red)
        assert abs(1 / math.sqrt(fP) - (2.0 * math.log10(Red * math.sqrt(fP)) - 0.8)) < 1e-10        # Prandtl's law
        fL = WT.pipe_friction_factor_turbulent(Red, kappa=KAP, B=BLOG)
        Uav_over_us = math.sqrt(8 / fL)
        assert abs(Uav_over_us - (math.log(Red / 2 / Uav_over_us) / KAP + BLOG - 1.5 / KAP)) < 1e-9  # f = 8(u_*/U_av)², a⁺ = Re_d/(2U_av⁺)
        assert abs(1 / math.sqrt(fL) - (slope * math.log10(Red * math.sqrt(fL)) + icpt)) < 1e-9      # the same law in Prandtl's form
        assert fP > 64 / Red and (prev_f is None or fP < prev_f)
        assert 1.0 < fL / fP < 1.10       # the log law integrated over the whole pipe over-predicts f by 5–9 % (why −1.02 was tuned to −0.8)
        prev_f = fP
    bm = json.loads((REF / "benchmarks.json").read_text(encoding="utf-8"))["prandtl_pipe_law"]
    assert abs(slope - bm["slope_log10"]) < 0.02 and abs(icpt - bm["intercept"]) < 0.25              # 1.99 and −1.02 vs 2.0 and −0.8


# ======================================================================================================================
# §12.10  C12 — eddy viscosity and the mixing length (12.94)–(12.102)
# ======================================================================================================================
def test_mixing_length_V2_derivation_D21_root_limits_closed_form_and_intercept():  # V2 · C12 · D21 (12.100)–(12.101)
    s, yp, kap = sp.symbols("s y_plus kappa", positive=True)
    roots = sp.solve(sp.Eq(s + kap ** 2 * yp ** 2 * s ** 2, 1), s)                 # ν dU/dy + κ²y²(dU/dy)² = u_*² in wall units
    pos = [r for r in roots if r.subs({kap: 0.41, yp: 10.0}) > 0]
    assert len(pos) == 1
    slope = 2 / (1 + sp.sqrt(1 + 4 * kap ** 2 * yp ** 2))
    assert sp.simplify(pos[0] - slope) == 0                                        # the positive root, stable form
    assert sp.limit(slope, yp, 0) == 1                                             # sublayer (12.82)
    assert sp.limit(slope * kap * yp, yp, sp.oo) == 1                              # → 1/(κy⁺): log layer (12.101)
    Ucf = (sp.asinh(2 * kap * yp) - (sp.sqrt(1 + 4 * kap ** 2 * yp ** 2) - 1) / (2 * kap * yp)) / kap
    assert sp.simplify(sp.diff(Ucf, yp) - slope) == 0 and sp.limit(Ucf, yp, 0) == 0   # the closed-form first integral
    Bund = sp.limit(Ucf - sp.log(yp) / kap, yp, sp.oo)
    assert sp.simplify(Bund - (sp.log(4 * kap) - 1) / kap) == 0                    # intercept without damping
    assert rel(ch12.mixing_length_intercept(KAP), (math.log(4 * KAP) - 1) / KAP) < 1e-13
    for y_ in (0.3, 10.0, 250.0):
        p = ch12.mixing_length_wall_profile(y_, KAP)
        assert rel(p["slope"], float(slope.subs({kap: KAP, yp: y_}))) < 1e-13 and rel(p["Uplus"], float(Ucf.subs({kap: KAP, yp: y_}))) < 1e-11
        assert abs(p["slope"] + p["lT_plus"] ** 2 * p["slope"] ** 2 - 1) < 1e-14 and rel(p["lT_plus"], KAP * y_) < 1e-15
        assert rel(p["uv_plus"], 1 - p["slope"]) < 1e-12 and rel(p["nuT_over_nu"], p["lT_plus"] ** 2 * p["slope"]) < 1e-12
        assert isinstance(p["Uplus"], float) and rel(p["B"], (math.log(4 * KAP) - 1) / KAP) < 1e-12


def test_mixing_length_V1_van_driest_profile_by_independent_quadrature():  # V1 · C12 (12.100) with damping
    A = 26.0
    s = lambda q: 2 / (1 + math.sqrt(1 + 4 * (KAP * q * (1 - math.exp(-q / A))) ** 2))  # noqa: E731
    for y_ in (1.0, 12.0, 100.0, 3000.0):
        Uq = quad(s, 0, y_, points=[1, 5, 12, 30, 100][: 1 + int(y_ > 1) + int(y_ > 12) * 2 + int(y_ > 100)], limit=400,
                  epsabs=1e-13, epsrel=1e-13)[0]
        p = ch12.mixing_length_wall_profile(y_, KAP, damping="van_driest", A_plus=A)
        assert rel(p["Uplus"], Uq) < 1e-8 and rel(p["slope"], s(y_)) < 1e-12
        assert p["Uplus"] > ch12.mixing_length_wall_profile(y_, KAP)["Uplus"]          # damping leaves more viscous sublayer
    Bq = quad(s, 0, 2e4, points=[1, 5, 12, 30, 100, 1000], limit=600, epsabs=1e-12, epsrel=1e-12)[0] - math.log(2e4) / KAP
    Bd = ch12.mixing_length_intercept(KAP, A_plus=A)
    assert abs(Bd - Bq) < 2e-4 and 5.0 < Bd < 5.5 and ch12.mixing_length_intercept(KAP, A_plus=0) == ch12.mixing_length_intercept(KAP)
    assert rel(ch12.mixing_length_wall_profile(2e4, KAP, damping="van_driest", A_plus=A)["B"], Bd) < 1e-10
    arr = ch12.mixing_length_wall_profile(np.array([1.0, 12.0, 100.0]), KAP, damping="van_driest")
    assert arr["Uplus"].shape == (3,) and rel(arr["Uplus"][1], ch12.mixing_length_wall_profile(12.0, KAP, damping="van_driest")["Uplus"]) < 1e-12


def test_mixing_length_V3_trapezoid_profile_converges_at_second_order():  # V3 · C12 (12.100)
    errs, hs = [], []
    for n in (201, 401, 801, 1601):
        yp = np.linspace(0.0, 60.0, n)
        num = ch12.mixing_length_wall_profile(yp, KAP, method="trapezoid")["Uplus"]
        ex = ch12.mixing_length_wall_profile(yp, KAP)["Uplus"]
        errs.append(np.max(np.abs(num - ex)))
        hs.append(yp[1])
    assert abs(observed_order(hs, errs) - 2.0) < 0.15
    # the notebook's from-scratch version: root point by point, then cumulative trapezoid
    yp = np.linspace(0.0, 60.0, 6001)
    root = (-1 + np.sqrt(1 + 4 * (KAP * yp) ** 2)) / np.where(yp > 0, 2 * (KAP * yp) ** 2, 1.0)
    root[0] = 1.0
    Uc = np.concatenate([[0.0], np.cumsum(0.5 * (root[1:] + root[:-1]) * np.diff(yp))])
    assert np.max(np.abs(Uc - ch12.mixing_length_wall_profile(yp, KAP)["Uplus"])) < 1e-5


def test_eddy_viscosity_V1_hypothesis_and_mixing_length_relations():  # V1 · C12 (12.94)–(12.98), N169
    G = np.array([[0.2, 3.0, 0.0], [0.5, -0.5, 1.0], [0.0, 0.4, 0.3]])              # divergence-free: trace 0
    uu = ch12.eddy_viscosity_stress(G, 0.02, 1.5)
    assert rel(np.trace(uu), 2 * 1.5) < 1e-14 and np.allclose(uu, uu.T)
    assert rel(uu[0, 1], -0.02 * 3.5) < 1e-14 and rel(uu[0, 0], 1.0 - 0.02 * 0.4) < 1e-14
    shear = np.zeros((3, 3))
    shear[0, 1] = 4.0
    uvs = ch12.eddy_viscosity_stress(shear, 0.02, 1.5)[0, 1]
    assert rel(uvs, -0.08) < 1e-14 and rel(ch12.eddy_viscosity_from_data(uvs, 4.0), 0.02) < 1e-13
    assert math.isnan(ch12.eddy_viscosity_from_data(0.0, 0.0))
    assert ch12.shear_production(ch12.eddy_viscosity_stress(shear, 0.02, 1.5), shear) > 0     # ν_T > 0 ⇒ production > 0
    assert ch12.gradient_diffusion_flux(-0.01, 5.0) == 0.05 and ch12.gradient_diffusion_flux(0.01, 5.0) == -0.05
    assert np.allclose(ch12.gradient_diffusion_flux(np.array([1.0, -2.0]), 0.5), [-0.5, 1.0])
    assert rel(ch12.eddy_diffusivity_estimate(0.2, 1.5, c=0.5), 0.15) < 1e-15
    assert (Q_(0.2, "m") * Q_(1.5, "m/s")).check("[length]**2/[time]")
    for dU in (3.0, -3.0):
        st = ch12.mixing_length_stress(dU, 0.1)
        nt = ch12.mixing_length_eddy_viscosity(dU, 0.1)
        assert rel(st, 0.01 * 9.0 * np.sign(dU)) < 1e-14 and rel(nt, 0.03) < 1e-14 and rel(st, nt * dU) < 1e-14
    assert ch12.mixing_length_stress(-3.0, 0.1) == -ch12.mixing_length_stress(3.0, 0.1)        # odd in the shear
    w = ch12.convective_velocity_scale(1000.0, 2.0, 300.0)
    assert rel(w, math.sqrt(G0 * 1000.0 * 2.0 / 300.0)) < 1e-14 and rel(ch12.convective_eddy_diffusivity(1000.0, 2.0, 300.0), w * 1000.0) < 1e-14
    assert ((Q_(G0, "m/s**2") * Q_(1000.0, "m") * Q_(2.0, "K") / Q_(300.0, "K")) ** 0.5).check("[length]/[time]")
    oe = ch12.one_equation_closure(1.44, 0.1, 0.55, 0.17, 1.0)
    assert rel(oe["u_T"], 0.55 * 1.2) < 1e-14 and rel(oe["nu_T"], 0.1 * 0.66) < 1e-14 and rel(oe["eps"], 0.17 * 1.2 ** 3 / 0.1) < 1e-13
    assert rel(oe["transport_diffusivity"], oe["nu_T"]) < 1e-15


def test_shear_flow_solver_V3_poiseuille_exact_and_second_order_with_an_eddy_viscosity():  # V1/V3 · C12 (12.99)
    h, nu, rho, dPdx = 1.0, 0.02, 1.3, -0.7
    y = np.sort(np.concatenate([[0.0, h], np.random.default_rng(0).uniform(0.02, 0.98, 37)]))     # any spacing
    lam = ch12.shear_flow_eddy_viscosity_solve(y, lambda yf, s: 0.0 * yf, dPdx, rho, nu)
    assert lam["converged"] and np.max(np.abs(lam["U"] + dPdx / (2 * rho * nu) * y * (h - y))) < 1e-10   # ν_T = 0: the parabola
    nuT = lambda yy: 0.05 * (1 + np.sin(PI * yy / h))  # noqa: E731 smooth, symmetric
    exact = lambda yy: quad(lambda q: dPdx / rho * (q - h / 2) / (nu + nuT(q)), 0, yy, epsabs=1e-13)[0]  # noqa: E731
    errs, hs = [], []
    for n in (21, 41, 81, 161):
        yg = np.linspace(0.0, h, n)
        sol = ch12.shear_flow_eddy_viscosity_solve(yg, lambda yf, s: nuT(yf), dPdx, rho, nu)
        errs.append(max(abs(sol["U"][i] - exact(yg[i])) for i in range(0, n, (n - 1) // 10)))
        hs.append(yg[1])
    assert abs(observed_order(hs, errs) - 2.0) < 0.15
    cou = ch12.shear_flow_eddy_viscosity_solve(np.linspace(0, h, 21), lambda yf, s: 0.0 * yf, 0.0, rho, nu, bc=(0.0, 2.0))
    assert np.max(np.abs(cou["U"] - 2.0 * np.linspace(0, h, 21))) < 1e-12                              # Couette


def test_shear_flow_solver_V1_mixing_length_rule_reproduces_the_model_channel():  # V1 · C12 (12.99) + (12.100): two solvers
    Re = 550.0
    c = ch12.channel_mixing_length(Re, KAP, n=400)
    y = np.concatenate([c["yplus"], (2 * Re - c["yplus"])[-2::-1]])                # both halves, wall units (ν = 1, ρ = 1)

    def nuT(yf, s):
        yw = np.minimum(yf, 2 * Re - yf)
        l = np.minimum(KAP * yw, 0.09 * Re) * (1 - np.exp(-yw / 26.0))  # noqa: E741
        return l ** 2 * np.abs(s)

    sol = ch12.shear_flow_eddy_viscosity_solve(y, nuT, -1.0 / Re, 1.0, 1.0, tol=1e-11, max_iter=4000)
    assert sol["converged"]
    half = sol["U"][: c["yplus"].size]
    assert np.max(np.abs(half - c["Uplus"])) < 5e-3 * c["U_cl_plus"]              # iterative BVP vs the closed-form root
    assert np.max(np.abs(sol["U"] - sol["U"][::-1])) < 1e-6                       # symmetric about the centreline


@needs_dns
def test_channel_mixing_length_V5_against_lee_moser_dns_approximate():  # V5 (approximate) · C12 (N181) · item 7
    out = {}
    for case in (550, 1000, 2000, 5200):
        d = dns(case)
        c = ch12.channel_mixing_length(d["Re_tau"], KAP, n=800)
        dU = np.interp(d["yplus"], c["yplus"], c["Uplus"]) - d["Uplus"]
        Ub_dns = np.trapezoid(d["Uplus"], d["y_over_delta"]) / d["y_over_delta"][-1]
        out[case] = (float(np.max(np.abs(dU))), c["Cf"] / (2 / Ub_dns ** 2) - 1)
        assert np.max(np.abs(dU)) < 1.0                              # within one wall unit everywhere
        assert abs(out[case][1]) < 0.05                              # friction coefficient within 5 % (measured ≤ 1.2 %)
        assert 8 < c["yplus_peak_production"] < 14
    d = dns(180)                                                     # at Re_τ = 180 there is no log layer: the model is poorer
    c = ch12.channel_mixing_length(d["Re_tau"], KAP, n=800)
    assert 0.5 < np.max(np.abs(np.interp(d["yplus"], c["yplus"], c["Uplus"]) - d["Uplus"])) < 2.0


# ======================================================================================================================
# §12.10  C13 — the k–ε model (12.103)–(12.105)
# ======================================================================================================================
def test_k_epsilon_V2_derivation_D23_decay_exponent_and_log_layer_constant():  # V2 · C13 · D23 (and D22 units)
    t, t0, n, e0, C2, C1, Cmu, sig, us, kap, y = sp.symbols("t t0 n e0 C_eps2 C_eps1 C_mu sigma_eps u_* kappa y", positive=True)
    e = e0 * (1 + t / t0) ** (-n)
    eps = -sp.diff(e, t)                                              # dē/dt = −ε̄
    resid = sp.simplify((sp.diff(eps, t) + C2 * eps ** 2 / e) / (eps ** 2 / e))      # dε̄/dt = −C_ε2 ε̄²/ē
    assert sp.solve(resid, n) == [1 / (C2 - 1)]                       # n = 1/(C_ε2 − 1)
    assert z0(eps.subs(t, 0) - n * e0 / t0)                           # t0 = n e0/ε0
    # log layer: P = ε̄ = u_*³/(κy), ē = u_*²/√C_μ, ν_T = C_μ ē²/ε̄ in (12.105) with its diffusion term
    eps_l = us ** 3 / (kap * y)
    e_l = us ** 2 / sp.sqrt(Cmu)
    nuT = Cmu * e_l ** 2 / eps_l
    assert z0(nuT - kap * us * y)                                     # the mixing-length value ν_T = κ u_* y
    eq = sp.diff(nuT / sig * sp.diff(eps_l, y), y) + C1 * eps_l * eps_l / e_l - C2 * eps_l ** 2 / e_l
    ksol = sp.solve(sp.simplify(eq * y ** 2 / us ** 4), kap)
    assert len(ksol) == 1 and z0(ksol[0] ** 2 - sp.sqrt(Cmu) * (C2 - C1) * sig)       # κ² = √C_μ (C_ε2 − C_ε1) σ_ε
    K = ch12.K_EPSILON_CONSTANTS
    assert rel(ch12.k_epsilon_loglayer_kappa(), math.sqrt(math.sqrt(K["C_mu"]) * (K["C_eps2"] - K["C_eps1"]) * K["sigma_eps"])) < 1e-14
    assert 0.38 < ch12.k_epsilon_loglayer_kappa() < 0.45
    nut = dimensional_check(lambda e, eps: 0.09 * e ** 2 / eps, "kinematic_viscosity", e=Q_(1.5, "m**2/s**2"), eps=Q_(0.3, "m**2/s**3"))
    assert rel(ch12.k_epsilon_eddy_viscosity(1.5, 0.3), nut.magnitude) < 1e-14
    assert (Q_(1.5, "m**2/s**2") ** 1.5 / Q_(0.3, "m**2/s**3")).check("[length]") and rel(ch12.k_epsilon_length_scale(1.5, 0.3), 1.5 ** 1.5 / 0.3) < 1e-14


def test_k_epsilon_decay_V1_closed_form_ivp_and_rk4_from_scratch():  # V1 · C13 (12.103), (12.105)
    e0, eps0, C2 = 1.0, 2.0, 1.92
    t = np.linspace(0.0, 10.0, 41)
    e, eps, n, t0 = ch12.k_epsilon_decay(e0, eps0, t, C_eps2=C2)
    assert rel(n, 1 / 0.92) < 1e-14 and rel(t0, n * e0 / eps0) < 1e-14
    assert np.allclose(e, e0 * (1 + t / t0) ** -n, rtol=1e-13) and np.allclose(eps, eps0 * (1 + t / t0) ** -(n + 1), rtol=1e-13)
    ei, epsi, _, _ = ch12.k_epsilon_decay(e0, eps0, t, C_eps2=C2, method="ivp")
    assert np.max(np.abs(ei / e - 1)) < 1e-7 and np.max(np.abs(epsi / eps - 1)) < 1e-7
    sol = solve_ivp(lambda tt, q: ch12.k_epsilon_rhs(q[0], q[1], 0.0), (0, 10), [e0, eps0], rtol=1e-10, atol=1e-13, t_eval=t)
    assert np.max(np.abs(sol.y[0] / e - 1)) < 1e-7                              # the model's source terms give the same decay
    q, hstep = np.array([e0, eps0]), 0.005                                      # RK4 loop (the notebook's from-scratch cell)
    f = lambda q: np.array([-q[1], -C2 * q[1] ** 2 / q[0]])  # noqa: E731
    for _ in range(2000):
        k1 = f(q); k2 = f(q + hstep / 2 * k1); k3 = f(q + hstep / 2 * k2); k4 = f(q + hstep * k3)  # noqa: E702
        q = q + hstep / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    assert rel(q[0], e[-1]) < 1e-8 and rel(q[1], eps[-1]) < 1e-8
    de, deps = ch12.k_epsilon_rhs(1.5, 0.3, 0.3)
    assert de == 0.0 and rel(deps, (1.44 - 1.92) * 0.3 ** 2 / 1.5) < 1e-13      # P = ε̄: ē steady, ε̄ not (diffusion balances it)
    de, deps = ch12.k_epsilon_rhs(1.5, 0.3, 0.5)
    assert rel(de, 0.2) < 1e-13 and rel(deps, 1.44 * 0.5 * 0.3 / 1.5 - 1.92 * 0.09 / 1.5) < 1e-13
    assert isinstance(ch12.k_epsilon_decay(1.0, 1.0, 10.0)[0], float)
    lt = ch12.k_epsilon_length_scale(e, eps)
    assert np.all(np.diff(lt) > 0) and np.all(np.diff(ch12.k_epsilon_eddy_viscosity(e, eps)) < 0)   # eddies grow, ν_T decays (n > 1)


def test_k_epsilon_V5_constants_are_the_published_standard_set():  # V5 · C13 (N178)
    bm = json.loads((REF / "benchmarks.json").read_text(encoding="utf-8"))["k_epsilon_constants"]
    K = ch12.K_EPSILON_CONSTANTS
    assert (K["C_mu"], K["C_eps1"], K["C_eps2"], K["sigma_e"], K["sigma_eps"]) == (
        bm["C_mu"], bm["C_eps1"], bm["C_eps2"], bm["sigma_k"], bm["sigma_eps"])
    assert len(K) == 5                                                          # five constants (slip #9: the page lists six)


def test_k_epsilon_channel_V7_optional_demo_converges_qualitative():  # V7 · C13 (optional; not in the contract) — qualitative
    c = ch12.k_epsilon_channel(1000.0, B=BLOG)
    assert c["converged"] and rel(c["kappa_model"], ch12.k_epsilon_loglayer_kappa()) < 1e-14
    assert np.all(np.diff(c["Uplus"]) > 0) and np.all(c["e_plus"] > 0) and np.all(c["eps_plus"] > 0)
    assert rel(c["e_plus"][0], 1 / math.sqrt(0.09)) < 1e-10 and 0.003 < c["Cf"] < 0.008


# ======================================================================================================================
# §12.11  C14 — Richardson numbers (12.106)–(12.109)
# ======================================================================================================================
GAMMA_A = STRAT.adiabatic_lapse_rate()            # Kundu's sign: −g/C_p ≈ −9.76 K/km (never typed rounded)


def test_flux_richardson_V1_hand_value_signs_and_degenerate_inputs():  # V1/V7 · C14 (12.107)
    al = 1 / 300.0
    Rf = ch12.flux_richardson(-0.02, -0.09, 0.075, al)
    assert rel(Rf, (G0 * al * 0.02) / (0.09 * 0.075)) < 1e-14 and Rf > 0              # downward heat flux: stable
    assert rel(ch12.flux_richardson(+0.02, -0.09, 0.075, al), -Rf) < 1e-14            # upward heat flux: Rf < 0
    assert rel(ch12.flux_richardson(-0.02, +0.09, -0.075, al), Rf) < 1e-14            # mirrored shear: same production
    assert rel(Rf, -ch12.buoyant_production(-0.02, al) / ch12.shear_production(np.array([[0, -0.09], [-0.09, 0]]),
                                                                                 np.array([[0, 0.075], [0, 0]]))) < 1e-13
    assert math.isnan(ch12.flux_richardson(0.0, -0.09, 0.0, al)) and math.isinf(ch12.flux_richardson(-0.02, -0.09, 0.0, al))
    arr = ch12.flux_richardson(np.array([-0.02, 0.0, 0.02]), -0.09, 0.075, al)
    assert arr[0] > 0 and arr[1] == 0 and arr[2] < 0
    assert (Q_(G0, "m/s**2") / Q_(300.0, "K") * Q_(0.02, "K*m/s") / (Q_(0.09, "m**2/s**2") * Q_(0.075, "1/s"))).check("[]")
    # regimes, both sides of each boundary and the boundaries themselves
    assert [ch12.turbulence_regime(v) for v in (-2.0, -1e-12, 0.0, 0.2499, 0.25, 3.0)] == [
        "convective", "convective", "shear-driven", "shear-driven", "decaying", "decaying"]
    assert ch12.turbulence_regime(float("nan")) == "undefined" and ch12.turbulence_regime(0.3, Rf_cr=0.5) == "shear-driven"
    assert list(ch12.turbulence_regime(np.array([-1.0, 0.1, 1.0]))) == ["convective", "shear-driven", "decaying"]


def test_richardson_V2_derivation_D24_reduced_budget_and_ri_equals_prandtl_times_rf():  # V2 · C14 · D24 (12.106)–(12.109)
    g, al, nuT, kT, dUdz, dthdz = sp.symbols("g alpha nu_T kappa_T dUdz dthetadz", positive=True)
    uw = -nuT * dUdz                                                  # (12.94) for U(z)
    wT = -kT * dthdz                                                  # (12.95)
    Rf = (-g * al * wT) / (-uw * dUdz)                                # (12.107)
    Ri = g * al * dthdz / dUdz ** 2                                   # (12.108), N² = gα dθ/dz
    assert z0(Ri - (nuT / kT) * Rf)                                   # (12.109)
    eps = sp.Symbol("epsilon", positive=True)
    P = -uw * dUdz
    assert z0((P + g * al * wT - eps) - (P * (1 - Rf) - eps))         # steady (12.106) without transport: ε̄ = P(1 − Rf)
    assert rel(ch12.turbulent_prandtl(0.6, 0.4), 1.5) < 1e-15
    assert rel(ch12.flux_from_gradient_richardson(0.3, Pr_T=1.5), 0.2) < 1e-15 and ch12.flux_from_gradient_richardson(0.3) == 0.3
    # numbers: build the fluxes from eddy coefficients and recover Ri = Pr_T Rf through the two functions
    nu_T, k_T, dU, dth, a = 0.6, 0.4, 0.05, 0.01, 1 / 290.0
    Rf_n = ch12.flux_richardson(-k_T * dth, -nu_T * dU, dU, a)
    Ri_n = ch12.gradient_richardson_thermal(dth, dU, a)["Ri"]
    assert rel(Ri_n, ch12.turbulent_prandtl(nu_T, k_T) * Rf_n) < 1e-13 and rel(ch12.flux_from_gradient_richardson(Ri_n, 1.5), Rf_n) < 1e-13


def test_stratified_tke_budget_V1_terms_and_optional_shear_argument():  # V1 · C14 (12.106) · flagged item 10
    al = 1 / 300.0
    one = ch12.stratified_tke_budget(10.0, 5.0, -0.09, -0.02, 0.002, al, dUdz=0.075)        # one height, shear given
    assert rel(one["shear_production"], 0.09 * 0.075) < 1e-14 and rel(one["buoyancy"], -G0 * al * 0.02) < 1e-14
    assert one["dissipation"] == -0.002 and rel(one["residual"], -(one["shear_production"] + one["buoyancy"] - 0.002)) < 1e-13
    assert rel(one["Rf"], ch12.flux_richardson(-0.02, -0.09, 0.075, al)) < 1e-14 and one["regime"] == "shear-driven"
    assert isinstance(one["Rf"], float) and isinstance(one["regime"], str)
    z = np.linspace(2.0, 50.0, 97)
    us, kap = 0.3, 0.4
    U = us / kap * np.log(z / 0.03)
    uw = -us ** 2 * np.ones_like(z)
    wT = -0.02 * np.ones_like(z)
    eps = us ** 3 / (kap * z)
    prof = ch12.stratified_tke_budget(z, U, uw, wT, eps, al)                                 # shear by differencing U
    given = ch12.stratified_tke_budget(z, U, uw, wT, eps, al, dUdz=us / (kap * z))           # shear supplied (the new option)
    # differencing a logarithmic profile: central differences of ln z are off by (Δz/z)²/3 (relative), no more
    dz = z[1] - z[0]
    assert np.all(np.abs(prof["shear_production"][1:-1] / given["shear_production"][1:-1] - 1) < 0.4 * (dz / z[1:-1]) ** 2)
    assert np.allclose(given["shear_production"], us ** 3 / (kap * z), rtol=1e-13)
    assert np.allclose(given["Rf"], ch12.flux_richardson(wT, uw, us / (kap * z), al), rtol=1e-13)
    assert np.allclose(given["residual"], -given["buoyancy"], rtol=1e-10)                    # P = ε̄ here, so the rest is buoyancy
    L_M = ch12.monin_obukhov_length(us, -0.02, alpha=al, kappa=kap)
    assert np.allclose(given["Rf"], z / L_M, rtol=1e-12)                                     # (12.111) from (12.107)
    assert set(np.unique(given["regime"])) <= {"shear-driven", "decaying"} and given["regime"][0] == "shear-driven"
    with pytest.raises(Exception):
        ch12.stratified_tke_budget(10.0, 5.0, -0.09, -0.02, 0.002, al)                       # one height needs dUdz


def _comparator(text: str) -> str:
    return re.search(r":\s*\S+\s*([<>=])\s*\S+\s*K/km", text).group(1)


def _numbers(text: str) -> list:
    tail = text.split(":", 1)[1].replace("−", "-")
    return [float(v) for v in re.findall(r"-?\d+\.?\d*", tail)][:2]


def test_gradient_richardson_V1_in_situ_and_potential_routes_and_both_conventions():  # V1 · C14 (12.108) · flagged item 9
    T, dU = 288.0, 0.01
    al = 1 / T
    assert -9.80e-3 < GAMMA_A < -9.70e-3                                 # Kundu's sign: negative
    for dTdz in (-6.5e-3, 0.0, +10e-3, -12e-3, GAMMA_A):
        r = ch12.gradient_richardson_thermal(dTdz, dU, al, Gamma_a=GAMMA_A)
        via_theta = ch12.gradient_richardson_thermal(dTdz - GAMMA_A, dU, al)           # Γ_a = 0: already a θ-gradient
        assert abs(r["Ri"] - via_theta["Ri"]) < 1e-12 and abs(r["N2"] - G0 * al * (dTdz - GAMMA_A)) < 1e-18
        assert abs(r["dthetadz"] - (dTdz - GAMMA_A)) < 1e-18 and r["dtheta_dz"] == r["dthetadz"]
        assert abs(r["N2"] - STRAT.brunt_vaisala_sq_from_lapse(T, dTdz)) < 1e-12        # ch01's N² from the lapse rate
        word = "neutral" if dTdz == GAMMA_A else ("stable" if dTdz > GAMMA_A else "unstable")
        assert r["verdict"] == word == via_theta["verdict"]
        # each text states the criterion ("stable ⇔ …") and then the numbers with the comparison that actually holds;
        # read that comparison back: Kundu's sign: > stable, < unstable; meteorological sign: the reverse
        assert r["verdict_kundu"].startswith("stable ⇔ dT/dz > Γa:") and r["verdict_met"].startswith("stable ⇔ Γ < Γa:")
        ck, cm = _comparator(r["verdict_kundu"]), _comparator(r["verdict_met"])
        assert {">": "stable", "<": "unstable", "=": "neutral"}[ck] == word          # both conventions: the same verdict
        assert {"<": "stable", ">": "unstable", "=": "neutral"}[cm] == word
        nk, nm = _numbers(r["verdict_kundu"]), _numbers(r["verdict_met"])
        assert np.allclose(nk, [-v for v in nm], atol=1e-9) and abs(nk[0] - dTdz * 1e3) < 0.051 and abs(nk[1] - GAMMA_A * 1e3) < 0.051
        k = STRAT.lapse_rate_stability(dTdz, GAMMA_A, convention="kundu")
        m = STRAT.lapse_rate_stability(dTdz, GAMMA_A, convention="meteorology")
        assert r["verdict_kundu"] == k.text and r["verdict_met"] == m.text and k.verdict == m.verdict == word
        assert "dT/dz" in r["verdict_kundu"] and "Γ" in r["verdict_met"] and r["verdict_kundu"] != r["verdict_met"]
        rm = ch12.gradient_richardson_thermal(dTdz, dU, al, Gamma_a=GAMMA_A, convention="meteorology")
        assert rm["text"] == r["verdict_met"] and rm["text_other"] == r["verdict_kundu"] and r["text"] == r["verdict_kundu"]
        assert rm["Ri"] == r["Ri"]                                        # the convention changes the words, never the number
    std = ch12.gradient_richardson_thermal(-6.5e-3, dU, al, Gamma_a=GAMMA_A)
    assert std["verdict"] == "stable" and std["Ri"] > 0                    # cooling at 6.5 K/km with height is STABLE
    assert "−6.5" in std["verdict_kundu"] and "6.5" in std["verdict_met"] and "−6.5" not in std["verdict_met"]
    mutant = ch12.gradient_richardson_thermal(-6.5e-3, dU, al)             # forgetting Γ_a calls the standard atmosphere unstable
    assert mutant["verdict"] == "unstable" and mutant["Ri"] < 0
    assert ch12.gradient_richardson_thermal(0.0, dU, al, Gamma_a=GAMMA_A)["Ri"] > std["Ri"]     # isothermal: more stable still
    assert math.isinf(ch12.gradient_richardson_thermal(0.01, 0.0, al)["Ri"]) and ch12.gradient_richardson_thermal(0.01, 0.0, al)["Ri"] > 0
    assert ch12.gradient_richardson_thermal(-0.01, 0.0, al)["Ri"] == -math.inf and math.isnan(ch12.gradient_richardson_thermal(0.0, 0.0, al)["Ri"])
    arr = ch12.gradient_richardson_thermal(np.array([-12e-3, -6.5e-3]), dU, al, Gamma_a=GAMMA_A)
    assert arr["Ri"][0] < 0 < arr["Ri"][1] and arr["verdict_kundu"] is None


# ======================================================================================================================
# §12.11  C15 — Monin–Obukhov length and the stratified surface layer (12.110)–(12.114)
# ======================================================================================================================
def test_monin_obukhov_V1_definition_sign_and_surface_layer_richardson_number_D25():  # V1/V2 · C15 · D25 (12.110), (12.111)
    us, kap, T = 0.3, 0.4, 300.0
    for wT in (-0.025, +0.1):
        L = ch12.monin_obukhov_length(us, wT, T=T, kappa=kap)
        assert rel(L, -us ** 3 / (kap * (1 / T) * G0 * wT)) < 1e-14 and np.sign(L) == -np.sign(wT)
        assert rel(ch12.monin_obukhov_length(us, wT, alpha=1 / T, kappa=kap), L) < 1e-15
        for z in (2.0, 10.0, 60.0):
            Rf_log = ch12.flux_richardson(wT, -us ** 2, us / (kap * z), 1 / T)       # (12.107) with the log-law gradient
            assert rel(ch12.flux_richardson_surface_layer(z, L), Rf_log) < 1e-13 and rel(Rf_log, z / L) < 1e-13
        assert rel(ch12.monin_obukhov_from_fluxes(1.2 * us ** 2, 1.2 * 1005.0 * wT, 1.2, 1005.0, T, kappa=kap), L) < 1e-13
    assert ch12.monin_obukhov_length(us, 0.0, T=T, kappa=kap) == math.inf and ch12.flux_richardson_surface_layer(10.0, math.inf) == 0.0
    assert (Q_(us, "m/s") ** 3 / (Q_(1 / T, "1/K") * Q_(G0, "m/s**2") * Q_(0.025, "K*m/s"))).check("[length]")
    with pytest.raises(TypeError):
        ch12.monin_obukhov_length(us, -0.02, T=T)                                    # κ has no default
    # symbolic: Rf = z/L_M, and integrating φ_m = 1 + βz/L_M gives the log-linear profile
    z, z0s, LM, beta, u, k = sp.symbols("z z0 L_M beta u_* kappa", positive=True)
    zz = sp.Symbol("zeta_", positive=True)
    prof = sp.integrate(u / (k * zz) * (1 + beta * zz / LM), (zz, z0s, z))
    assert z0(sp.expand_log(prof - u / k * (sp.log(z / z0s) + beta * (z - z0s) / LM), force=True))
    g_, a_, w_ = sp.symbols("g alpha wT")
    Lsym = -u ** 3 / (k * a_ * g_ * w_)
    assert z0((-g_ * a_ * w_) / (u ** 2 * u / (k * z)) - z / Lsym)
    assert rel(ch12.gradient_richardson_surface_layer(10.0, 80.0, Pr_T=1.3), 1.3 * 10.0 / 80.0) < 1e-15
    assert rel(ch12.gradient_richardson_surface_layer(10.0, 80.0, Pr_T=1.3, consistent=True), 1.3 * 0.125 / 1.625) < 1e-14
    assert rel(ch12.gradient_richardson_surface_layer(1e9, 80.0, consistent=True), 1 / 5.0) < 1e-6      # saturates at Pr_T/β
    assert math.isnan(ch12.gradient_richardson_surface_layer(30.0, -80.0, consistent=True))             # φ_m ≤ 0


def test_surface_layer_wind_V1_log_linear_profile_limits_and_shear_function():  # V1/V7 · C15 (N188)
    us, z0_, kap = 0.3, 0.03, 0.4
    z = np.array([0.03, 1.0, 10.0, 50.0])
    neutral = WT.surface_layer_wind(z, us, z0_, kappa=kap)
    assert np.allclose(neutral, WT.rough_wall_log_law(z, us, z0_, kappa=kap), atol=1e-14)               # L_M = ∞: (12.93)
    assert np.allclose(WT.surface_layer_wind(z, us, z0_, -math.inf, kappa=kap), neutral)
    st = WT.surface_layer_wind(z, us, z0_, 80.0, kappa=kap)
    un = WT.surface_layer_wind(z, us, z0_, -400.0, kappa=kap)
    assert np.allclose(st, us / kap * (np.log(z / z0_) + 5.0 * z / 80.0)) and np.all(st[1:] > neutral[1:]) and np.all(un[1:] < neutral[1:])
    assert rel(WT.surface_layer_wind(10.0, us, z0_, 80.0, kappa=kap, beta=4.7), us / kap * (math.log(10 / z0_) + 4.7 * 10 / 80.0)) < 1e-14
    assert math.isnan(WT.surface_layer_wind(0.01, us, z0_, 80.0, kappa=kap)) and isinstance(WT.surface_layer_wind(10.0, us, z0_, 80.0, kappa=kap), float)
    # φ_m = (κz/u_*) dU/dz, by differences of the profile, on each branch
    for L, mode in ((80.0, "log_linear"), (-400.0, "log_linear"), (-20.0, "businger_dyer"), (80.0, "businger_dyer")):
        for zz in (2.0, 10.0, 40.0):
            dz = 1e-4 * zz
            dUdz = (WT.surface_layer_wind(zz + dz, us, z0_, L, kappa=kap, unstable=mode)
                    - WT.surface_layer_wind(zz - dz, us, z0_, L, kappa=kap, unstable=mode)) / (2 * dz)
            assert rel(kap * zz / us * dUdz, WT.dimensionless_shear(zz / L, unstable=mode)) < 1e-6, (L, mode, zz)
    assert WT.dimensionless_shear(0.0) == 1.0 and WT.dimensionless_shear(0.0, unstable="businger_dyer") == 1.0
    assert rel(WT.dimensionless_shear(0.2), 2.0) < 1e-15 and rel(WT.dimensionless_shear(-0.5, unstable="businger_dyer"), 9.0 ** -0.25) < 1e-15
    assert math.isnan(WT.dimensionless_shear(-0.2)) and math.isnan(WT.dimensionless_shear(-1.0))        # 1 + βζ ≤ 0: no meaning
    assert WT.surface_layer_wind(10.0, us, z0_, -5.0, kappa=kap) < WT.surface_layer_wind(5.0, us, z0_, -5.0, kappa=kap)  # outside its range


def test_surface_layer_wind_V5_businger_dyer_integral_with_and_without_the_z0_term():  # V1/V5 · C15 (N188) · flagged item 4
    bm = json.loads((REF / "benchmarks.json").read_text(encoding="utf-8"))["businger_dyer_unstable"]
    us, z0_, kap = 0.3, 0.03, 0.4
    phi = lambda zeta: (1 - bm["coefficient"] * zeta) ** bm["exponent"]  # noqa: E731  φ_m = (1 − 16ζ)^(−1/4), ζ < 0
    assert rel(WT.dimensionless_shear(-0.37, unstable="businger_dyer"), phi(-0.37)) < 1e-15

    def psi(zeta):                  # ψ_m(ζ) = ∫_0^ζ (1 − φ_m)/ζ′ dζ′ by quadrature (independent of the closed form)
        return quad(lambda s: (1 - phi(s)) / s, 0.0, zeta, epsabs=1e-13, epsrel=1e-13)[0]

    for L in (-5.0, -20.0, -200.0):
        for z in (1.0, 10.0, 50.0):
            exact = us / kap * quad(lambda s: phi(s / L) / s, z0_, z, epsabs=1e-13, epsrel=1e-13)[0]     # ∫_{z0}^{z} φ_m u_*/(κz) dz
            full = WT.surface_layer_wind(z, us, z0_, L, kappa=kap, unstable="businger_dyer", psi_at_z0=True)
            common = WT.surface_layer_wind(z, us, z0_, L, kappa=kap, unstable="businger_dyer")
            assert rel(full, exact) < 1e-10                                              # with ψ_m(z0/L): the exact integral
            assert rel(common, us / kap * (math.log(z / z0_) - psi(z / L))) < 1e-10      # default: the cited common form
            assert abs((full - common) - us / kap * psi(z0_ / L)) < 1e-11                # they differ by (u_*/κ)ψ_m(z0/L) only
            assert 0 < full - common < 4.2 * us / kap * z0_ / abs(L)                      # … of order z0/|L_M| (ψ ≈ −4ζ, small ζ)
            assert common < us / kap * math.log(z / z0_)                                 # unstable: less wind than neutral
        assert WT.surface_layer_wind(z0_, us, z0_, L, kappa=kap, unstable="businger_dyer", psi_at_z0=True) == 0.0
        u0 = WT.surface_layer_wind(z0_, us, z0_, L, kappa=kap, unstable="businger_dyer")
        assert u0 < 0 and rel(u0, -us / kap * psi(z0_ / L)) < 1e-9                       # default: U(z0) = −(u_*/κ)ψ_m(z0/L) ≠ 0
    # the stable side is the book's log-linear form in both modes
    assert WT.surface_layer_wind(10.0, us, z0_, 80.0, kappa=kap, unstable="businger_dyer") == WT.surface_layer_wind(10.0, us, z0_, 80.0, kappa=kap)


def test_surface_layer_regime_V7_boundaries_both_signs_and_neutral():  # V7 · C15 (N187) · flagged item 5
    R = ch12.surface_layer_regime
    assert R(10.0, math.inf) == "neutral" and R(10.0, -math.inf) == "neutral"
    assert R(10.0, 100.0) == "forced convection" and R(10.0, -100.0) == "forced convection"       # z ≪ |L_M|, either sign
    assert R(9.999, 10.0) == "forced convection" and R(9.999, -10.0) == "forced convection"
    assert R(10.0, 10.0) == "stable" and R(10.0, -10.0) == "free convection"                       # the boundary |z/L_M| = 1
    assert R(50.0, 10.0) == "stable" and R(50.0, -10.0) == "free convection"
    assert R(10.0, 100.0, crossover=0.05) == "stable" and R(10.0, -100.0, crossover=0.05) == "free convection"
    assert list(R(np.array([1.0, 30.0, 30.0]), np.array([-20.0, -20.0, 20.0]))) == ["forced convection", "free convection", "stable"]
    # "forced convection" is not a statement about static stability: the same string for a stable and an unstable layer
    s_st = ch12.surface_layer_state(0.3, -30.0, 288.0, 0.03, 10.0, 1.2, 1005.0, 0.4)
    s_un = ch12.surface_layer_state(0.3, +30.0, 288.0, 0.03, 10.0, 1.2, 1005.0, 0.4)
    assert s_st["layer"] == s_un["layer"] == "forced convection" and s_st["verdict"] == "stable" and s_un["verdict"] == "unstable"


def test_surface_layer_state_V1_composition_signs_and_both_verdicts():  # V1/V7 · C15 (E9 entry point) · flagged item 9
    us, T, z0_, z, rho, cp, kap = 0.3, 288.0, 0.03, 10.0, 1.2, 1005.0, 0.4
    s = ch12.surface_layer_state(us, -30.0, T, z0_, z, rho, cp, kap)
    wT = -30.0 / (rho * cp)
    L = ch12.monin_obukhov_length(us, wT, T=T, kappa=kap)
    assert rel(s["wT"], wT) < 1e-14 and rel(s["L_M"], L) < 1e-13 and rel(s["zeta"], z / L) < 1e-13 and rel(s["Rf"], z / L) < 1e-13
    assert rel(s["Ri"], s["Rf"]) < 1e-14 and rel(s["phi_m"], 1 + 5 * z / L) < 1e-13 and rel(s["Rf_profile"], (z / L) / s["phi_m"]) < 1e-13
    assert rel(s["U"], WT.surface_layer_wind(z, us, z0_, L, kappa=kap)) < 1e-13 and rel(s["U_neutral"], us / kap * math.log(z / z0_)) < 1e-13
    assert s["regime"] == ch12.turbulence_regime(s["Rf"]) and s["layer"] == ch12.surface_layer_regime(z, L) and s["valid"] is True
    assert rel(s["z_crit"], 0.25 * L) < 1e-13 and s["z_Rf_cr"] == s["z_crit"]
    assert rel(s["Gamma_a"], -G0 / cp) < 1e-14 and rel(s["dtheta_dz"], -wT * s["phi_m"] / (kap * us * z)) < 1e-13
    assert rel(s["dTdz"], s["dtheta_dz"] + s["Gamma_a"]) < 1e-13
    # flux–gradient consistency: κ_T = −wT/(dθ/dz) = κ u_* z/φ_m (Pr_T = 1)
    assert rel(-wT / s["dtheta_dz"], kap * us * z / s["phi_m"]) < 1e-13
    k = STRAT.lapse_rate_stability(s["dTdz"], s["Gamma_a"], convention="kundu")
    m = STRAT.lapse_rate_stability(s["dTdz"], s["Gamma_a"], convention="meteorology")
    assert s["verdict_kundu"] == k.text and s["verdict_met"] == m.text and s["verdict"] == k.verdict == m.verdict == "stable"
    assert s["text"] == s["verdict_kundu"] and s["text_other"] == s["verdict_met"]
    sm = ch12.surface_layer_state(us, -30.0, T, z0_, z, rho, cp, kap, convention="meteorology")
    assert sm["text"] == s["verdict_met"] and sm["U"] == s["U"]
    u = ch12.surface_layer_state(us, +30.0, T, z0_, z, rho, cp, kap)                 # daytime: everything flips sign
    assert rel(u["L_M"], -L) < 1e-13 and u["Rf"] < 0 and u["regime"] == "convective" and u["verdict"] == "unstable"
    assert _comparator(u["verdict_kundu"]) == "<" and _comparator(u["verdict_met"]) == ">"      # unstable in both conventions
    assert _comparator(s["verdict_kundu"]) == ">" and _comparator(s["verdict_met"]) == "<"      # stable in both
    assert u["U"] < u["U_neutral"] < s["U"] and math.isnan(u["z_crit"]) and u["dTdz"] < u["Gamma_a"] < 0
    n = ch12.surface_layer_state(us, 0.0, T, z0_, z, rho, cp, kap)                   # neutral limit
    assert math.isinf(n["L_M"]) and n["Rf"] == 0 and rel(n["U"], n["U_neutral"]) < 1e-14 and n["verdict"] == "neutral"
    assert n["layer"] == "neutral" and n["regime"] == "shear-driven"
    strong = ch12.surface_layer_state(us, +300.0, T, z0_, z, rho, cp, kap)           # z/L_M < −1/β: outside the log-linear range
    assert strong["valid"] is False and math.isnan(strong["phi_m"])
    bd = ch12.surface_layer_state(us, +300.0, T, z0_, z, rho, cp, kap, unstable="businger_dyer")
    assert bd["valid"] is True and 0 < bd["phi_m"] < 1 and 0 < bd["U"] < bd["U_neutral"]
    assert rel(bd["U"], WT.surface_layer_wind(z, us, z0_, bd["L_M"], kappa=kap, unstable="businger_dyer")) < 1e-12


def test_scalar_spectrum_V2_exponents_slopes_and_batchelor_scale():  # V2/V1 · C15 (12.113), (12.114), N193
    a, b, c = sp.symbols("a b c")
    # S_T [K² m] = ε̄_T^a [K² s⁻¹]^a · ε̄^b [m² s⁻³]^b · K^c [m⁻¹]^c
    sol = sp.solve([sp.Eq(2 * a, 2), sp.Eq(2 * b - c, 1), sp.Eq(-a - 3 * b, 0)], [a, b, c])
    assert (sol[a], sol[b], sol[c]) == (1, sp.Rational(-1, 3), sp.Rational(-5, 3))
    assert ch12.scalar_spectrum_exponents() == {"eps_T": Fraction(1), "eps": Fraction(-1, 3), "K": Fraction(-5, 3)}
    nu, kth, eps, epsT = 1e-6, 1e-9, 1e-6, 1e-7                # ν/κ = 1000 (a salt-like scalar): a decade and a half of K⁻¹
    eta = (nu ** 3 / eps) ** 0.25
    etaT = ch12.batchelor_scale(nu, kth, eps)
    assert rel(etaT, eta * math.sqrt(kth / nu)) < 1e-14 and etaT < eta
    K1 = np.geomspace(1e-3, 1e-1, 20) / eta
    K2 = np.geomspace(1.5, 0.5 * eta / etaT, 20) / eta          # 1/η < K < 1/η_T
    S1, S2 = ch12.scalar_spectrum(K1, eps, epsT, nu, kth), ch12.scalar_spectrum(K2, eps, epsT, nu, kth)
    assert abs(np.polyfit(np.log(K1), np.log(S1), 1)[0] + 5 / 3) < 1e-12 and abs(np.polyfit(np.log(K2), np.log(S2), 1)[0] + 1.0) < 1e-12
    assert rel(S1[0], epsT * eps ** (-1 / 3) * K1[0] ** (-5 / 3)) < 1e-13                    # (12.113) with C_T = 1
    assert rel(S2[0], epsT * math.sqrt(nu / eps) / K2[0]) < 1e-13                           # (12.114) amplitude
    j = 1 / eta
    assert rel(ch12.scalar_spectrum(j * (1 - 1e-9), eps, epsT, nu, kth), ch12.scalar_spectrum(j * (1 + 1e-9), eps, epsT, nu, kth)) < 1e-8
    assert ch12.scalar_spectrum(10 / eta, eps, epsT, nu, kth, cutoff=True) < ch12.scalar_spectrum(10 / eta, eps, epsT, nu, kth)
    gas = ch12.scalar_spectrum(np.array([2.0, 20.0]) / eta, eps, epsT, 1.5e-5, 2.1e-5)        # ν/κ < 1: only the −5/3 branch
    assert rel(gas[0] / gas[1], 10 ** (5 / 3)) < 1e-12
    assert rel(ch12.scalar_spectrum(K1[3], eps, epsT, nu, kth, C_T=0.7), 0.7 * S1[3]) < 1e-13
    assert (Q_(epsT, "K**2/s") * Q_(eps, "m**2/s**3") ** Fraction(-1, 3) * Q_(5.0, "1/m") ** Fraction(-5, 3)).check("[temperature]**2*[length]")


# ======================================================================================================================
# §12.12  C16 — Taylor's theory of turbulent dispersion (12.115)–(12.129)
# ======================================================================================================================
def test_taylor_dispersion_V2_derivation_D26_double_integral_by_parts():  # V2 · C16 · D26 (12.117)–(12.119)
    t, tp, tau, u2, Lam, tc = sp.symbols("t t_p tau u2 Lambda t_c", positive=True)
    r = sp.Function("r")
    single = 2 * u2 * sp.Integral((t - tau) * r(tau), (tau, 0, t))                   # (12.119): 2u2 t∫(1 − τ/t) r dτ
    rate = sp.diff(single, t).doit()
    assert z0(rate - 2 * u2 * sp.Integral(r(tau), (tau, 0, t)))                      # (12.117): d mean(X²)/dt
    assert z0(sp.diff(rate, t).doit() - 2 * u2 * r(t)) and single.subs(t, 0).doit() == 0
    for rr, closed in ((sp.exp(-tau / Lam), 2 * u2 * Lam ** 2 * (t / Lam - 1 + sp.exp(-t / Lam))),
                       (sp.exp(-tau ** 2 / tc ** 2), 2 * u2 * (t * sp.sqrt(sp.pi) * tc / 2 * sp.erf(t / tc)
                                                              - tc ** 2 / 2 * (1 - sp.exp(-t ** 2 / tc ** 2))))):
        one = 2 * u2 * sp.integrate((t - tau) * rr, (tau, 0, t))                                 # (12.119)
        two = 2 * u2 * sp.integrate(sp.integrate(rr, (tau, 0, tp)), (tp, 0, t))                  # (12.118)
        assert sp.simplify(one - two) == 0 and sp.simplify(one - closed) == 0
    # D27: the two limits of the exponential case, with the constant offset of the diffusive branch
    Xe = 2 * u2 * Lam ** 2 * (t / Lam - 1 + sp.exp(-t / Lam))
    assert z0(sp.series(Xe, t, 0, 4).removeO() - (u2 * t ** 2 - u2 * t ** 3 / (3 * Lam)))        # (12.120): u2 t² first
    assert sp.limit(Xe - (2 * u2 * Lam * t - 2 * u2 * Lam ** 2), t, sp.oo) == 0                  # (12.122) − 2u2∫τ r dτ
    assert z0(sp.integrate(tau * sp.exp(-tau / Lam), (tau, 0, sp.oo)) - Lam ** 2)
    # D28: D_T = ½ d mean(X²)/dt (12.127) and its limits (12.128), (12.129, condition corrected)
    DT = sp.diff(Xe, t) / 2
    assert z0(DT - u2 * Lam * (1 - sp.exp(-t / Lam))) and sp.limit(DT, t, sp.oo) == u2 * Lam and z0(sp.series(DT, t, 0, 2).removeO() - u2 * t)


def test_taylor_dispersion_V1_closed_forms_two_integral_forms_and_limits():  # V1/V7 · C16 (12.117)–(12.123), (12.127)
    u2, Lam, tc = 0.36, 12.0, 9.0
    r_e = lambda s: np.exp(-np.abs(s) / Lam)  # noqa: E731
    r_g = lambda s: np.exp(-(s / tc) ** 2)    # noqa: E731
    for t in (0.5, 6.0, 12.0, 90.0):
        Xe = ch12.taylor_dispersion_exponential(t, u2, Lam)
        assert rel(Xe, 2 * u2 * Lam ** 2 * (t / Lam - 1 + math.exp(-t / Lam))) < 1e-12
        assert rel(ch12.taylor_dispersion(t, r_e, u2), Xe) < 1e-9 and rel(ch12.taylor_dispersion(t, r_e, u2, form="double"), Xe) < 1e-8
        mine = 2 * u2 * quad(lambda s: (t - s) * math.exp(-s / Lam), 0, t, epsabs=0, epsrel=1e-13)[0]
        assert rel(Xe, mine) < 1e-11
        Xg = ch12.taylor_dispersion_gaussian(t, u2, tc)
        assert rel(ch12.taylor_dispersion(t, r_g, u2), Xg) < 1e-9
        assert rel(Xg, 2 * u2 * quad(lambda s: (t - s) * math.exp(-(s / tc) ** 2), 0, t, epsabs=0, epsrel=1e-13)[0]) < 1e-11
        assert rel(ch12.taylor_dispersion_rate(t, r_e, u2), 2 * u2 * Lam * (1 - math.exp(-t / Lam))) < 1e-10      # (12.117)
        De = ch12.eddy_diffusivity_exponential(t, u2, Lam)
        assert rel(ch12.eddy_diffusivity_taylor(t, r_e, u2), De) < 1e-10 and rel(De, 0.5 * ch12.taylor_dispersion_rate(t, r_e, u2)) < 1e-10
        h = 1e-4 * t
        num = (ch12.taylor_dispersion_exponential(t + h, u2, Lam) - ch12.taylor_dispersion_exponential(t - h, u2, Lam)) / (4 * h)
        assert rel(De, num) < 1e-7                                                   # D_T = ½ d mean(X²)/dt (12.127)
        sl = t / Xe * 2 * num                                        # d ln mean(X²)/d ln t = t (2 D_T)/mean(X²)
        assert rel(ch12.dispersion_local_slope(t, Lam), sl) < 1e-6
    for small in (1e-6, 1e-3, 1e-2):                                                 # ballistic (12.120)–(12.121)
        t = small * Lam
        assert abs(ch12.taylor_dispersion_exponential(t, u2, Lam) / (u2 * t * t) - (1 - small / 3)) < small ** 2
        assert abs(ch12.taylor_dispersion_gaussian(t, u2, tc) / (u2 * t * t) - 1) < 2 * (t / tc) ** 2
        assert abs(ch12.eddy_diffusivity_exponential(t, u2, Lam) / (u2 * t) - 1) < small                 # (12.128)
    t = 500 * Lam                                                                    # diffusive (12.122)–(12.123), with the offset
    assert rel(ch12.taylor_dispersion_exponential(t, u2, Lam), 2 * u2 * Lam * t - 2 * u2 * Lam ** 2) < 1e-13
    assert abs(ch12.taylor_dispersion_exponential(t, u2, Lam) / (2 * u2 * Lam * t) - 1) < 1 / 499
    Lg = math.sqrt(PI) / 2 * tc
    assert abs(ch12.taylor_dispersion_gaussian(500 * tc, u2, tc) - (2 * u2 * Lg * 500 * tc - u2 * tc ** 2)) < 1e-8
    assert rel(ch12.eddy_diffusivity_exponential(t, u2, Lam), u2 * Lam) < 1e-12 and rel(ch12.eddy_diffusivity_taylor(900.0, r_g, u2), u2 * Lg) < 1e-9
    assert abs(ch12.dispersion_local_slope(1e-4, 1.0) - 2.0) < 1e-4 and abs(ch12.dispersion_local_slope(1e5, 1.0) - 1.0) < 1e-4
    assert abs(ch12.dispersion_local_slope(0.03, 1.0) - (2 - 0.01)) < 1e-4            # 2 − x/3 for small x
    arr = ch12.taylor_dispersion_exponential(np.array([0.0, 6.0]), u2, Lam)
    assert arr[0] == 0 and arr.shape == (2,) and isinstance(ch12.taylor_dispersion_exponential(6.0, u2, Lam), float)
    tt = np.linspace(0.0, 30.0, 301)
    assert np.allclose(ch12.diffusivity_from_variance(tt, 2 * 0.7 * tt), 0.7, rtol=1e-12)          # σ² = 2Dt (12.126)
    Dn = ch12.diffusivity_from_variance(tt, ch12.taylor_dispersion_exponential(tt, u2, Lam))
    # central differences of ½ mean(X²): truncation (Δt²/6)|D_T″| ≤ (Δt²/6) u2/Λ_t
    assert np.max(np.abs(Dn[2:-2] - ch12.eddy_diffusivity_exponential(tt[2:-2], u2, Lam))) < 1.01 * (0.1 ** 2 / 6) * u2 / Lam
    assert (Q_(u2, "m**2/s**2") * Q_(Lam, "s")).check("[length]**2/[time]") and (Q_(u2, "m**2/s**2") * Q_(Lam, "s") * Q_(6.0, "s")).check("[length]**2")


def test_langevin_particles_stat_dispersion_matches_taylor_formula_within_5_standard_errors():  # stat/V4 · C16 (12.115)–(12.119)
    n, urms, Lam = 20_000, 0.6, 10.0
    t = np.concatenate([[0.0], np.geomspace(0.05, 100.0, 40)])                      # uneven steps: the update is exact
    X, u = ch12.langevin_particles(n, t, urms, Lam, seed=7)
    assert X.shape == (n, t.size) == u.shape and np.all(X[:, 0] == 0)
    X2 = np.mean(X ** 2, axis=0)
    se = np.std(X ** 2, axis=0, ddof=1) / math.sqrt(n)
    exact = ch12.taylor_dispersion_exponential(t, urms ** 2, Lam)
    assert np.all(np.abs(X2[1:] - exact[1:]) < 5 * se[1:])                           # (12.119) at every output time
    assert abs(np.mean(u ** 2, axis=0) - urms ** 2).max() < 5 * urms ** 2 * math.sqrt(2 / n) * 1.3   # stationary velocity
    for j in (10, 25, 35):                                                           # Lagrangian correlation e^{−τ/Λ}
        rr = TS.correlation_coefficient(u[:, 5], u[:, j])
        rex = math.exp(-(t[j] - t[5]) / Lam)
        assert abs(rr - rex) < 5 * (1 - rex ** 2) / math.sqrt(n) + 5 * abs(np.mean(u[:, 5])) / urms / math.sqrt(n)
    Xu = np.mean(X * u, axis=0)
    seXu = np.std(X * u, axis=0, ddof=1) / math.sqrt(n)
    DT = ch12.eddy_diffusivity_exponential(t, urms ** 2, Lam)
    assert np.all(np.abs(Xu[1:] - DT[1:]) < 5 * seXu[1:])                            # mean(Xu) = D_T (12.116), (12.127)
    tu = np.linspace(0.0, 60.0, 601)
    Xb, ub = ch12.langevin_particles(n, tu, urms, Lam, seed=8)
    lhs, rhs = ch12.dispersion_rate_from_particles(tu, Xb, ub)
    assert np.allclose(rhs, 2 * np.mean(Xb * ub, axis=0)) and lhs.shape == rhs.shape == tu.shape
    assert np.max(np.abs(lhs[5:-5] - rhs[5:-5])) < 0.05 * np.max(rhs)                # (12.116): the same curve, noise ~ N^(−1/2)
    assert abs(np.mean(lhs[100:-5] - rhs[100:-5])) < 0.005 * np.max(rhs)             # and no systematic difference
    X3, u3 = ch12.langevin_particles(2000, tu[:50], urms, Lam, seed=1, dim=3)
    assert X3.shape == (2000, 50, 3) and abs(np.mean(X3[:, -1, 0] * X3[:, -1, 1])) < 5 * np.std(X3[:, -1, 0] * X3[:, -1, 1]) / math.sqrt(2000)
    assert np.array_equal(ch12.langevin_particles(10, tu[:5], urms, Lam, seed=3)[0], ch12.langevin_particles(10, tu[:5], urms, Lam, seed=3)[0])


def test_random_walk_stat_mean_square_distance_grows_as_n_steps():  # stat · C16 (12.124), (12.125) (a-D45)
    nw, ns, L = 20_000, 200, 0.5
    for dim in (1, 2, 3):
        R = ch12.random_walk(ns, nw, L=L, dim=dim, seed=dim)
        assert R.shape == (nw, ns + 1, dim) and np.all(R[:, 0] == 0)
        steps = np.diff(R, axis=1)
        assert np.max(np.abs(np.sqrt((steps ** 2).sum(axis=2)) - L)) < 1e-12        # every step has length L
        R2 = (R ** 2).sum(axis=2)
        for k in (1, 10, 50, 200):
            se = np.std(R2[:, k], ddof=1) / math.sqrt(nw)
            assert abs(R2[:, k].mean() - k * L * L) < 5 * se + 1e-12, (dim, k)      # (12.125): mean(R_n²) = nL²
        assert abs(np.mean((R[:, 50] * steps[:, 50]).sum(axis=1))) < 5 * np.std((R[:, 50] * steps[:, 50]).sum(axis=1)) / math.sqrt(nw)  # (12.124) cross term
    from_scratch = np.cumsum(np.random.default_rng(0).choice([-1.0, 1.0], size=(nw, ns)), axis=1)   # the notebook's loop
    assert abs(np.mean(from_scratch[:, -1] ** 2) - ns) < 5 * math.sqrt(2.0) * ns / math.sqrt(nw)
    p = 0.6
    Rp = ch12.random_walk(2000, 4000, L=1.0, dim=2, persistence=p, seed=4)
    R2p = (Rp ** 2).sum(axis=2)[:, -1]
    target = 2000 * (1 + p) / (1 - p)
    assert abs(R2p.mean() - target) < 5 * np.std(R2p, ddof=1) / math.sqrt(4000) + 2 * p / (1 - p) ** 2 * 1.5   # large-n law (+ O(1) offset)
    assert R2p.mean() > 3 * 2000                                                    # memory makes the walk spread faster


# ======================================================================================================================
# Drawing helpers, scripts, reference data and documentation
# ======================================================================================================================
def test_drawings_V7_helpers_draw_without_error():  # smoke · design C.3b (pure drawing, no physics)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sys.path.insert(0, str(ROOT / "scripts"))
    import ch12_drawings as DR

    names = ("scales_sketch", "parcel_sketch", "stress_element_sketch", "fg_geometry_sketch", "free_shear_sketch",
             "wall_layers_sketch", "surface_layer_sketch", "plume_sketch")
    for nm in names:
        fig, ax = plt.subplots()
        getattr(DR, nm)(ax)
        assert len(ax.lines) + len(ax.patches) + len(ax.collections) + len(ax.texts) > 0, nm
        plt.close(fig)


@slow
def test_scripts_and_make_refs_V7_run_clean(tmp_path):  # V7 smoke (scripts/ch12_*.py --no-show --fast)
    py = sys.executable
    for s in sorted((ROOT / "scripts").glob("ch12_*.py")):
        if s.name in ("ch12_common.py", "ch12_tables.py"):       # ch12_tables rewrites reference/ch12 (covered by the writer test)
            continue
        r = subprocess.run([py, str(s), "--no-show", "--fast", "--out", str(tmp_path)], cwd=str(ROOT), capture_output=True,
                           text=True, encoding="utf-8", errors="replace", timeout=900)
        assert r.returncode == 0, (s.name, r.stderr[-2000:])
    before = {f.name: f.read_bytes() for f in REF.glob("*") if f.suffix in (".csv", ".json")}
    r = subprocess.run([py, str(REF / "make_refs.py")], cwd=str(ROOT), capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=300)
    assert r.returncode == 0, r.stderr[-2000:]
    assert {f.name: f.read_bytes() for f in REF.glob("*") if f.suffix in (".csv", ".json")} == before     # idempotent


def test_reference_V5_sources_file_cites_every_data_file_and_no_book_table():  # data-and-benchmarks policy
    src = (REF / "SOURCES.md").read_text(encoding="utf-8")
    for name in ("lee_moser_2015_channel_mean.csv", "benchmarks.json", "explainer_tables.json"):
        assert name in src and (REF / name).exists()
    assert "Lee" in src and "774" in src and "turbulence.oden.utexas.edu" in src and "2026-10-07" in src
    bm = json.loads((REF / "benchmarks.json").read_text(encoding="utf-8"))
    assert bm["lee_moser_kappa"]["value"] == 0.384 and bm["lee_moser_cases"]["5200"]["Re_tau"] == 5185.897
    assert WT.LOG_LAW_CONSTANTS["channel_dns_lee_moser_2015"]["kappa"] == bm["lee_moser_kappa"]["value"]
    d = dns(180)
    assert abs(d["Re_tau"] - 182.088) < 1e-9 and d["yplus"][-1] < d["Re_tau"] and np.all(np.diff(d["yplus"]) > 0)


def test_docstrings_V1_no_stale_validation_phrases():  # ch01 lesson: stale "pending" text
    stale = re.compile(r"not yet pinned|benchmark pending|TODO|to be confirmed by the verifier", re.I)
    hits = [n for n in CONTRACT if inspect.isfunction(getattr(ch12, n)) and stale.search(inspect.getdoc(getattr(ch12, n)) or "")]
    assert hits == [], hits


# ======================================================================================================================
# V6 — numbers printed by the book (private; skipped when the JSON is absent; nothing is copied into this file)
# ======================================================================================================================
@book_only
def test_book_V6_cascade_free_shear_and_examples():  # V6 (C07, C08, C09, N126)
    b = book()
    mx = b["sec_12_7_cascade"]["mixer_example"]
    assert rel(ch12.kolmogorov_scales(mx["nu_m2_s"], mx["eps_m2_s3"])[0], mx["eta_m"]) < 5e-3
    kc = ch12.kolmogorov_constants(b["sec_12_7_cascade"]["kolmogorov_constant_3D"])
    assert abs(kc["C1_two_sided"] / b["sec_12_7_cascade"]["kolmogorov_constant_S11_double_sided"] - 1) < 0.02
    t = b["sec_12_8_free_shear"]["table_12_1"]
    for key, flow in (("planar_jet", "plane_jet"), ("planar_plume", "plane_plume"), ("round_jet", "round_jet"),
                      ("round_plume", "round_plume")):
        ex = ch12.free_shear_exponents(flow)
        assert abs(float(ex["velocity"]) - t[key]["U_exponent_x"]) < 1e-12 and abs(float(ex["scalar"]) - t[key]["Y_exponent_x"]) < 1e-12
    for key, flow in (("planar_wake", "plane_wake"), ("round_wake", "round_wake")):
        assert abs(float(ch12.free_shear_exponents(flow)["velocity"]) - t[key]["dU_exponent_x"]) < 1e-12
    e2 = b["sec_12_8_free_shear"]["example_12_2"]
    i, bk = e2["inputs"], e2["book"]
    st = ch12.stoichiometric_mass_fraction(i["M_CH4"], i["M_air"], 2.0, i["O2_volume_fraction"])       # CH4 + 2 O2
    assert abs(st["v_fuel"] - bk["v_CH4"]) < 5e-4 and abs(st["Y_fuel"] / bk["Y_CH4"] - 1) < 5e-3
    rj = t["round_jet"]
    x = ch12.round_jet_distance_for_mass_fraction(st["Y_fuel"], i["d_m"], i["M_CH4"], i["M_air"], 1.0, rj["C_Y"])
    assert abs(x / bk["x_m"] - 1) < 5e-3
    cl = ch12.free_shear_centerline("round_jet", x, constants={"C_U": rj["C_U"], "C_Y": rj["C_Y"]}, d=i["d_m"], U0=i["U0_m_s"],
                                    rho_s=i["M_CH4"], rho=i["M_air"])
    assert abs(cl["U_CL"] / bk["U_CL_m_s"] - 1) < 2e-2          # the book rounds to two digits
    assert abs(cl["Y_CL"] / st["Y_fuel"] - 1) < 1e-12


@book_only
def test_book_V6_wall_modelling_and_stratified_numbers():  # V6 (C10, C11, C12, C13, C15)
    b = book()
    w = b["sec_12_9_wall"]
    e3 = w["example_12_3"]
    i, bk = e3["inputs"], e3["book"]
    z = WT.zpg_boundary_layer(i["x_m"], i["U_m_s"], i["nu_air_m2_s"], kappa=i["kappa"])
    assert abs(z["theta"] / bk["theta_m"] - 1) < 5e-3 and abs(z["delta_star"] / bk["delta_star_m"] - 1) < 5e-3
    assert abs(z["delta99"] / bk["delta99_m"] - 1) < 5e-3 and abs(z["Re_x"] / bk["Re_x"] - 1) < 1e-2
    for flow in ("channel", "pipe", "zpg_boundary_layer"):
        pair = w["log_law_constants_by_flow"][flow]
        assert abs(WT.nagib_chauhan_kappa(pair["B"]) / pair["kappa"] - 1) < 0.02            # (12.92) holds for the book's pairs
    assert WT.layer_name(w["viscous_sublayer_limit_yplus"] - 1e-9, 1e-4) == "viscous sublayer"
    assert WT.layer_name(w["buffer_layer_yplus"][1], 1e-3) == "logarithmic layer"
    x34 = b["exercises"]["12_34"]
    slope = math.log(10) / (x34["kappa"] * math.sqrt(8))
    icpt = (x34["B"] - 1.5 / x34["kappa"] - math.log(2 * math.sqrt(8)) / x34["kappa"]) / math.sqrt(8)
    assert abs(slope - x34["derived_by_analyst"]["slope"]) < 1e-5 and abs(icpt - x34["derived_by_analyst"]["intercept"]) < 1e-5
    x33 = b["exercises"]["12_33"]
    tau0 = WT.wall_stress_from_pressure_gradient(x33["dpdx_Pa_m"], x33["d_m"])
    assert abs(tau0 / x33["hint_tau0_Pa"] - 1) < 5e-3
    cv = b["sec_12_10_modeling"]["convection_example"]
    assert abs(ch12.convective_velocity_scale(cv["L_m"], cv["dT_K"], 300.0) / cv["w_m_s"] - 1) < 0.06       # one-digit book value
    assert abs(ch12.convective_eddy_diffusivity(cv["L_m"], cv["dT_K"], 300.0) / cv["kappa_T_m2_s"] - 1) < 0.06
    kb = b["sec_12_10_modeling"]["k_epsilon_constants"]
    assert all(ch12.K_EPSILON_CONSTANTS[k] == kb[k] for k in ("C_mu", "C_eps1", "C_eps2", "sigma_e", "sigma_eps"))
    s = b["sec_12_11_stratified"]
    assert inspect.signature(ch12.turbulence_regime).parameters["Rf_cr"].default == s["Rf_critical"]
    assert inspect.signature(WT.surface_layer_wind).parameters["beta"].default == s["log_linear_coefficient"]
    x37 = b["exercises"]["12_37"]
    L = ch12.monin_obukhov_from_fluxes(x37["tau_Pa"], x37["heat_flux_W_m2"], 1.2, 1005.0, 300.0, kappa=0.41)
    assert abs(L / x37["derived_by_analyst_L_M_m_rho1p2_T300_kappa0p41"] - 1) < 5e-3 and L < 0
    x38 = b["exercises"]["12_38"]
    Lam = TS.correlation_spectrum_pair("gaussian", x38["u_rms_m_s"], x38["t_c_s"])["Lambda_t"]
    assert abs(Lam - x38["derived_by_analyst"]["Lambda_t_s"]) < 1e-6
    assert abs(ch12.eddy_diffusivity_taylor(50.0, lambda q: np.exp(-(q / x38["t_c_s"]) ** 2), x38["u_rms_m_s"] ** 2)
               - x38["derived_by_analyst"]["D_T_m2_s"]) < 1e-6
