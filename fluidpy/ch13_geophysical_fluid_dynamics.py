"""Chapter 13 — Geophysical Fluid Dynamics (Kundu, Cohen & Dowling 5e, §13.1–§13.18, Eqs. (13.1)–(13.145)).

Rotation and stratification on a thin shell: geostrophic and thermal-wind balance, Ekman layers, shallow-water and
internal waves with rotation (Poincaré, Kelvin, inertia–gravity, Rossby), vertical normal modes, potential vorticity,
barotropic and baroclinic (Eady) instability, and two-dimensional turbulence.

**This module re-exports every public name of the three core modules**

* ``fluidpy.core.gfd`` (``GFD``) — closed forms: planetary parameters, balances, Ekman layers, dispersion relations,
  Rossby radii, potential vorticity, the Eady problem, Fjørtoft and Rhines;
* ``fluidpy.core.vertical_modes`` (``VM``) — the Sturm–Liouville modes of N²(z);
* ``fluidpy.core.shallow_water`` (``SW``) — the numerical models (C-grid shallow water, the 1-D march, exact linear
  QG evolution, the pseudo-spectral barotropic model);

so that ``ch13.<name>`` reaches everything (explainer parity rows write ``ch13.<name>`` only).  One name exists in two
of them: ``ch13.potential_vorticity`` is the closed form ``GFD.potential_vorticity(zeta, f, h)``; the gridded
``SW.potential_vorticity(model, state)`` is ``ch13.sw_potential_vorticity``.  :func:`reexport_audit` checks this.

What is added here: the tables of the chapter (:func:`illustrative_inputs`, :func:`conventions_table`,
:func:`book_slips`, :func:`traps`, :func:`lapse_rate_table`), small helpers for §13.1–§13.5 and §13.13–§13.17, twelve
sympy engines that re-derive the chapter's results (two of them with a ``printed=True`` switch that reproduces a
printed slip and fails), figure helpers and the cached model runs of ``reference/ch13``.

Conventions: see ``fluidpy.core.gfd`` (x east, y north, z up; f signed, every direction-dependent function takes the
sign of f into account; SI units; latitudes in radians).  **Lapse rates**: computed in Kundu's convention Γ ≡ dT/dz
(Γ_a = −g/C_p ≈ −9.8 K/km, stable when dT/dz > Γ_a) and always shown beside the meteorological Γ_met ≡ −dT/dz
(Γ_d ≈ +9.8 K/km, stable when Γ_met < Γ_d); ``Gamma_a`` is a required keyword wherever a verdict is returned.
**Book slips** are implemented in corrected form; :func:`book_slips` lists them with the planted wrong variant a test
can fail.  Inputs used in examples are ours (:func:`illustrative_inputs`), never the book's.
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable

import numpy as np

from .core import gfd as GFD
from .core import shallow_water as SW
from .core import vertical_modes as VM
from .core.shallow_water import *  # noqa: F401,F403
from .core.vertical_modes import *  # noqa: F401,F403
from .core.gfd import *  # noqa: F401,F403  (last: GFD owns the bare name potential_vorticity)
from .core.shallow_water import potential_vorticity as sw_potential_vorticity  # noqa: F401
from .core._util import as_scalar_if_0d
from .core.gfd import (EARTH_RADIUS_MEAN, OMEGA_EARTH, beta_parameter, coriolis_parameter, ekman_depth,
                       ekman_surface, ekman_transport, geostrophic_velocity, inertia_gravity_m2, lee_wave_m,
                       rossby_long_wave_speed, thin_layer_terms)
from .core.stratification import adiabatic_lapse_rate, brunt_vaisala_sq_from_lapse, lapse_rate_stability
from .core.thermo import G0

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d

_OWN = ["illustrative_inputs", "conventions_table", "book_slips", "traps", "lapse_rate_table", "reexport_audit",
        "wind_from_to", "ocean_density_gradient_budget", "atmosphere_layers", "ocean_profile_idealized",
        "thermocline_N2", "anisotropic_eddy_stress", "friction_force_dimensions", "term_table_thin_layer",
        "pressure_centre", "jet_section", "viscous_layer_thicknesses", "ekman_surface_printed", "coastal_upwelling",
        "flow_over_step", "vertical_structure_solve", "wkb_error", "w_equation_rotating_residual", "lee_wave_field",
        "basin_crossing_time", "rayleigh_kuo_eigs", "eady_basic_state", "eady_matrix", "eady_numeric_eigs",
        "boussinesq_rotating_sympy", "perturbation_form_sympy", "eddy_friction_force_sympy", "thin_layer_sympy",
        "taylor_proudman_sympy", "hydrostatic_linear_set_sympy", "v_equation_sympy",
        "rotating_internal_wave_set_sympy", "w_equation_rotating_sympy", "pv_conservation_sympy", "qg_vorticity_sympy",
        "eady_qg_sympy", "sympy_engines", "fig_ekman_surface", "fig_ekman_bottom", "fig_mode_roots",
        "fig_vertical_modes", "fig_poincare_kelvin_dispersion", "fig_kelvin_sections", "fig_inertia_gravity_orbit",
        "fig_rossby_dispersion", "REFERENCE_RUNS", "reference_run", "write_reference_runs", "load_reference_run",
        "explainer_constants", "GFD", "VM", "SW"]
__all__ = sorted(set(GFD.__all__) | set(VM.__all__) | set(SW.__all__) | set(_OWN))


def _root() -> Path:
    from .core.project import repo_root

    return repo_root()


# =====================================================================================================================
# 0. Tables: inputs, conventions, slips, traps, the two-convention lapse-rate table
# =====================================================================================================================

def illustrative_inputs() -> dict:
    """The inputs every worked number, script and explainer of the chapter uses — ours, visibly not the book's.

    Book: none (the book's own example values stay in the private test file).
    Returns dict: ``lat`` (35° in rad), ``lat_ekman`` (60°), ``lat_rossby`` (12°) [rad]; ``ocean_H`` [m], ``ocean_N``
    [rad/s]; ``atm_H`` [m], ``atm_N`` [rad/s], ``atm_U0`` [m/s]; ``tau`` [N/m²], ``nu_v_ocean`` [m²/s]; ``U_g`` [m/s],
    ``nu_v_atm`` [m²/s]; ``U_syn`` [m/s], ``L_syn`` [m]; ``H1`` [m], ``drho`` [kg/m³], ``rho_ocean`` [kg/m³];
    ``U_mean`` [m/s]; ``wavelengths`` (two, [m]); ``q_inertial`` [m/s]; ``u_rms`` (atmosphere, ocean) [m/s].
    Use the southern-hemisphere twin of any latitude by negating it.
    Assumptions: none.   Validation: V1 a table.  Label: analytic.
    """
    return dict(lat=float(np.deg2rad(35.0)), lat_ekman=float(np.deg2rad(60.0)), lat_rossby=float(np.deg2rad(12.0)),
                ocean_H=4200.0, ocean_N=2.7e-3, atm_H=9000.0, atm_N=1.1e-2, atm_U0=27.0, tau=0.07, nu_v_ocean=0.03,
                U_g=12.0, nu_v_atm=7.0, U_syn=14.0, L_syn=1.4e6, H1=120.0, drho=3.1, rho_ocean=1027.0, U_mean=17.0,
                wavelengths=(7.3e5, 3.1e6), q_inertial=0.23, u_rms=(13.0, 0.08))


def conventions_table():
    """The overloaded letters of the chapter and how the notebook, the explainers and the code tell them apart.

    Book: §13.1–§13.18 (trap T13: the same letter means different things in different sections).
    Returns a ``pandas.DataFrame`` with columns "symbol", "meanings", "ours", "code".
    Assumptions: none.   Validation: V1 a table.  Label: analytic.
    """
    import pandas as pd

    rows = [
        ("f", "Coriolis parameter f = 2Ω sin θ, signed (positive north, negative south)",
         "f, and f0 at the central latitude; a frequency is always ω", "f, f0"),
        ("Ω, ω", "Ω rotation rate of the Earth; ω wave frequency; ω_x, ω_y horizontal vorticity (Ekman layer)",
         "Ω (sidereal value), ω, vorticity always with its subscript", "Omega, omega, omega_x"),
        ("β", "df/dy = 2Ω cos θ0 / R; in §13.3 also the haline contraction coefficient", "β = df/dy only; β_S for salt",
         "beta"),
        ("N", "buoyancy frequency from potential density", "N; counts are n, nx, n_modes", "N, N2"),
        ("H, h, η", "layer depth, ocean depth, WKB scale of N, lid separation, scale height c²/g; h total depth; "
         "η surface displacement", "H depth; H_N; H_s; H_e equivalent depth; h; η", "H, He, h, eta"),
        ("c", "sound speed (§13.2–13.3); long-wave speed sqrt(gH) (§13.8 on); modal speed c_n; complex phase speed "
         "(§13.16–13.17)", "c_s sound; c long-wave; c_n; c = c_r + i c_i said each time", "c, c_n"),
        ("k, l, m, K", "eastward, northward, vertical wavenumbers; K magnitude; K also the perturbation kinetic energy "
         "(§13.17) and K0, K1, K2 of Fjørtoft's triad", "k, l, m, K; KE for energy; K_v(z) eddy viscosity profile",
         "k, l, m, K"),
        ("ζ", "relative vorticity ∂v/∂x − ∂u/∂y", "ζ; ζ_g geostrophic", "zeta"),
        ("θ", "latitude (§13.4); angle of the wavevector above the horizontal (§13.14); potential temperature",
         "θ latitude; θ_K wavevector angle; θ_p potential temperature", "lat_rad"),
        ("δ", "Ekman e-folding thickness sqrt(2ν_v/|f|); the oceanographic Ekman depth is πδ", "δ; D_E = πδ",
         "ekman_depth(convention=)"),
        ("Λ, λ", "Rossby radius: c/|f| external, sqrt(g'H)/|f| or NH/(nπ|f|) internal, NH/|f| (no π) in §13.17; "
         "λ wavelength", "Λ with the wave speed named; Λ_E = NH/|f|", "rossby_radius*"),
        ("ψ", "vertical mode ψ_n(z) (§13.9); stream function with u = −∂ψ/∂y, v = ∂ψ/∂x (§13.5, §13.16)",
         "ψ_n modes; ψ stream function with its sign stated", "psi"),
        ("E, Ro", "Ekman number ν/(|f|L²); Rossby number U/(|f|L) (ch04 used U/(2ΩL))", "E, Ro",
         "ekman_number, rossby_number"),
        ("α", "thermal expansion (§13.3); Eady wavenumber α = N K/|f| (§13.17); enstrophy flux (§13.18)",
         "α_T; α (Eady); α_Z", "alpha (keyword), alphaH, alpha_ens"),
        ("U, V, τ", "velocity scale, geostrophic interior velocity, mean current, U0 z/H; V complex velocity u + iv; "
         "τ wind stress", "U named each time; V; V_g; τ", "U, U_g, V_g, tau_x, tau_y"),
        ("z = 0", "sea surface with the ocean below (§13.6, §13.9, §13.14); solid surface (§13.7); flat bottom "
         "(§13.8); lower lid (§13.17)", "stated at the top of each block", "per function docstring"),
        ("Γ", "the book quotes lapse rates in words (rates of decrease)", "Kundu Γ = dT/dz computed; meteorological "
         "Γ_met = −dT/dz shown beside it", "dT_dz, Gamma_a (keyword)"),
        ("S(K)", "energy spectrum with mean(u²) = ∫S dK: one-sided, no factor ½ (ch12 used two-sided spectra)",
         "stated beside every spectral statement", "enstrophy_spectrum, barotropic_spectrum"),
    ]
    return pd.DataFrame(rows, columns=["symbol", "meanings", "ours", "code"])


def book_slips() -> dict:
    """The printed slips of chapter 13 that fluidpy corrects, each with the way to tell the two versions apart.

    Book: §13.3–§13.17 (found by reading the page images against dimensions, limits and the book's own neighbouring
    equations; slips #1–#8 were re-read on the rendered pages when this module was written).
    Returns a dict keyed "1" … "13"; each value is a dict with ``id``, ``where``, ``printed``, ``correct`` (alias
    ``corrected``), ``how_to_tell``, ``taught_in``, ``coded_in`` (the function that carries the corrected form) and
    ``test`` (the planted wrong variant a test must fail; "—" for wording-only slips).  Descriptions are ours and
    contain no number from the book.  ``pd.DataFrame(book_slips()).T`` gives the notebook's table.
    Slip #13 (a historical attribution) is of medium confidence and is not asserted anywhere.
    Assumptions: none.   Validation: V1 every ``printed=True`` switch named here exists and fails its check.
    """
    rows = [
        ("(13.2), §13.3", r"momentum equation with $+\frac{1}{\rho_0}\nabla p$",
         r"$\frac{D\mathbf u}{Dt}+2\boldsymbol\Omega\times\mathbf u=-\frac{1}{\rho_0}\nabla p-\frac{g\rho}{\rho_0}\mathbf e_z+\mathbf F$",
         "a fluid at rest must satisfy the equation; the next two displays on the page carry the minus",
         "C01 (R04)", "boussinesq_rotating_sympy", "boussinesq_rotating_sympy(printed=True): the rest state leaves a residual"),
        ("§13.3 and (13.6)", r"$F_i=\partial\tau_{ij}/\partial x_j$ called a force per unit mass",
         r"$F_i=\frac{1}{\rho}\,\partial\tau_{ij}/\partial x_j$",
         "dimensions: stress gradient is N/m³, the right side of (13.6) is m/s²", "C01 (R10), D01",
         "eddy_friction, eddy_friction_force_sympy", "friction_force_dimensions(printed=True) is not an acceleration"),
        ("(13.17), §13.5", r"$-2\Omega u=-\frac1\rho\frac{\partial p}{\partial y}$",
         r"$2\Omega u=-\frac1\rho\frac{\partial p}{\partial y}$",
         "(13.12) with f = 2Ω; with the printed sign cross-differentiation gives ∂u/∂x − ∂v/∂y = 0, not the divergence",
         "C03 (N21), D05", "taylor_proudman_sympy", "taylor_proudman_sympy(printed=True) does not reach ∂w/∂z = 0"),
        ("(13.26), §13.6", r"$u,v\to0$ as $z\to\infty$", r"$u,v\to0$ as $z\to-\infty$",
         "z = 0 is the sea surface and the ocean lies below; the text after (13.29) says −∞", "C04 (N28), D06",
         "ekman_surface", "ekman_surface is NaN for z > 0; the root bounded upward violates the stress condition"),
        ("§13.4, after the inertial period", "the subscript of T_i is said to refer to a vector component",
         "it does not: i is a label for 'inertial'", "T_i is a scalar", "C01 (N15)", "inertial_period", "—"),
        ("§13.9, after (13.69)", r"the first root is said to occur for $NH/c_n=1$", r"$NH/c_0\ll1$",
         "the next step uses tan x ≈ x; the root is N sqrt(H/g), a few hundredths", "C08 (R22)", "modes_uniform_N",
         "N*H/modes_uniform_N(...).c[0] is below 0.1 for ocean-like inputs"),
        ("§13.10, end of the low-frequency paragraph",
         r"$\omega\gg f$ given as the range where the first term of the cubic is negligible", r"$\omega\ll f$",
         "the paragraph is about very slow waves; for ω ≫ f the ω³ term is one of the two dominant ones",
         "C09 (N83), D13", "dispersion_term_sizes, shallow_water_regime",
         "dispersion_term_sizes: ω³ is the largest term at ω = 3f and the smallest on the slow root"),
        ("§13.14, opening", "N is said to be taken depth independent", "depth dependent, N(z)",
         "(13.99) writes N²(z) and the WKB section needs it", "C14 (N105)", "inertia_gravity_m2, wkb_vertical_structure",
         "—"),
        ("cross-references", "sections and an equation of earlier chapters cited by a wrong number",
         "§4.9 (Boussinesq, buoyancy), §8.4 (impulsively started plate), §5.6 (vortex stretching), (7.127) (N²)",
         "section titles of the book's own table of contents", "conventions block", "—", "—"),
        ("§13.7, laminar Ekman thickness of air", "a thickness about a quarter below sqrt(2ν/f) for its own inputs",
         "the value of the formula", "recompute", "C06 (N48)", "ekman_depth", "private test asserts the formula"),
        ("§13.12, typical internal Rossby radius", "about three times NH/(πf) with the chapter's own typical N, H, f",
         "an inconsistency of inputs, not of the formula", "recompute", "C12 (R24)", "rossby_radius_internal",
         "private test records both"),
        ("exercise on the long Rossby-wave speed", "an answer that follows only from a round β",
         "with β at the stated latitude the speed is about a fifth smaller", "recompute both ways", "C15 (N135)",
         "rossby_long_wave_speed", "private test"),
        ("§13.17, opening", "the theory of baroclinic instability attributed to one author (first name misspelt)",
         "not asserted here (medium confidence): the problem solved in §13.17 is known as the Eady problem",
         "history of the subject — unread first-hand", "C16 opening", "—", "—"),
    ]
    out = {}
    for i, (where, printed, correct, how, taught, coded, test) in enumerate(rows, start=1):
        out[str(i)] = dict(id=i, where=where, printed=printed, correct=correct, corrected=correct, how_to_tell=how,
                           taught_in=taught, coded_in=coded, test=test)
    return out


def traps() -> list:
    """Things that are true as printed but easy to misread — kept apart from :func:`book_slips`.

    Book: §13.1–§13.18.   Returns a list of 17 dicts (``id`` "T1" … "T17", ``what``, ``where`` — the block that
    carries the callout, ``code`` — how the code guards against it).
    Assumptions: none.   Validation: V1 a table.  Label: analytic.
    """
    rows = [
        ("the book's Ω is one turn per solar day; the sidereal value is 0.27 % larger", "C01",
         "OMEGA_EARTH (sidereal) by default; OMEGA_SOLAR_DAY for comparison"),
        ("'close to neutral' compares the troposphere with the moist adiabat; it is stable to dry displacements; N² "
         "is from potential density", "§13.2 recap", "lapse_rate_table(dT_dz, Gamma_a=...) in both conventions"),
        ("from the thin-shell equations on, p and ρ are perturbations with the primes dropped", "C01, C03",
         "docstrings say 'perturbation'"),
        ("where z = 0 is changes from section to section", "C04, C06, C07, C08, C16", "stated per function"),
        ("every formula with sqrt(f), 'to the right', 'clockwise', exp(−fy/c) assumes f > 0", "C04, C11",
         "abs(f) and sign(f) everywhere; *_printed forms refuse f <= 0"),
        ("δ is an e-folding scale; the 'Ekman depth' is πδ", "C04", "ekman_depth(convention='efold' | 'pi')"),
        ("the 'orbit' figure is a velocity hodograph; the axis ratio is quoted as ω/f in §13.11 and f/ω in §13.14",
         "C10, C14", "poincare_orbit returns both the hodograph and the path"),
        ("with ω > 0 a Rossby wave has k < 0; the 'maximum' phase speed is a maximum of magnitude", "C15",
         "rossby_omega returns signed ω for signed k"),
        ("'all roots real, two superinertial' holds under β-plane scaling; the fast roots have opposite signs", "C09",
         "shallow_water_omega returns NaN where the discriminant is negative"),
        ("f is replaced by f0 except where it is differentiated; PV conservation itself needs no such step",
         "C09, C13", "pv_conservation_sympy uses an arbitrary f(y)"),
        ("three Rossby radii: external, internal (with π), and NH/f without π in the Eady problem", "C12, C16",
         "rossby_radius_internal(with_pi=)"),
        ("the modal amplitudes have different units: u_n, v_n [m/s], p_n [m²/s²], w_n [1/s], ρ_n [kg/m²]", "C08",
         "modal_amplitudes docstring"),
        ("overloaded letters (θ, α, β, H, c, K, ζ, η, ψ, l, U, τ, i)", "conventions block", "conventions_table()"),
        ("here u = −∂ψ/∂y, v = ∂ψ/∂x; ch04 and ch11 use the opposite sign", "C02, C15",
         "geostrophic_streamfunction, barotropic_velocity"),
        ("approximations that enter silently (continuity and uniform coefficients in the friction force; uniform "
         "density in the no-shear step; no vertical advection of momentum in the Eady set; sqrt(m) frozen in the WKB "
         "velocity)", "D01, D05, D25, C14", "docstrings list them under Assumptions"),
        ("§13.18's spectrum is one-sided with no factor ½; ch12's was two-sided", "C17",
         "barotropic_spectrum states its normalisation"),
        ("two Rossby numbers: U/(fL) here, U/(2ΩL) in ch04 — they differ by sin θ", "C02", "rossby_number docstring"),
    ]
    return [dict(id=f"T{i}", what=w, where=wh, code=c) for i, (w, wh, c) in enumerate(rows, start=1)]


def lapse_rate_table(dT_dz, *, Gamma_a):
    """The stability verdict of a temperature gradient in **both** lapse-rate conventions, as a two-row table.

    Book: §13.2 (lapse rates are quoted there in words, as rates of decrease); ch01 §1.10 for the criterion.
    Parameters
    ----------
    dT_dz : environment temperature gradient in Kundu's sign [K/m] (a standard troposphere is −6.5e-3).
    Gamma_a : adiabatic temperature gradient dT_a/dz in Kundu's sign [K/m] (negative: −g/C_p for dry air) —
        **required keyword, no default**; obtain it from ``core.stratification.adiabatic_lapse_rate()``.
    Returns
    -------
    ``pandas.DataFrame`` with index ("Kundu", "meteorology") and columns "symbol", "value_K_per_km",
    "adiabatic_K_per_km", "criterion" (our own string, e.g. "dT/dz = −6.5 > −9.8 K/km" and
    "Γ_met = 6.5 < Γ_d = 9.8 K/km"), "verdict" ("stable" / "neutral" / "unstable" — the same in both rows).
    Assumptions: dry (unsaturated) displacements; the verdict comes from
    ``core.stratification.lapse_rate_stability``.
    Validation: V7 both rows give the same verdict; neutral exactly at dT_dz = Gamma_a.  Label: analytic.
    """
    import pandas as pd

    x, ga = float(dT_dz), float(Gamma_a)
    res = lapse_rate_stability(x, ga)
    op_k = {1: ">", 0: "=", -1: "<"}[res.code]
    op_m = {1: "<", 0: "=", -1: ">"}[res.code]
    digits = 1
    while res.code != 0 and digits < 6 and f"{x * 1e3:.{digits}f}" == f"{ga * 1e3:.{digits}f}":
        digits += 1
    fm = lambda v: f"{v * 1e3:.{digits}f}".replace("-", "−")  # noqa: E731
    rows = [dict(symbol="dT/dz", value_K_per_km=x * 1e3, adiabatic_K_per_km=ga * 1e3,
                 criterion=f"dT/dz = {fm(x)} {op_k} {fm(ga)} K/km", verdict=res.verdict),
            dict(symbol="Γ_met = −dT/dz", value_K_per_km=-x * 1e3, adiabatic_K_per_km=-ga * 1e3,
                 criterion=f"Γ_met = {fm(-x)} {op_m} Γ_d = {fm(-ga)} K/km", verdict=res.verdict)]
    return pd.DataFrame(rows, index=["Kundu", "meteorology"])


def reexport_audit() -> dict:
    """Check that every public name of the three core modules is reachable from this module.

    Book: none (machinery; the chapter 12 lesson — a forgotten re-export left many names unreachable).
    Returns dict: ``total`` (public names in ``GFD``, ``VM``, ``SW``), ``same_object`` (reachable as the same
    object under the same name), ``aliased`` (dict name → alias for names that exist in two modules), ``missing``
    (list — must be empty).
    Assumptions: none.   Validation: V1 ``missing == []``.  Label: analytic.
    """
    import sys

    me = sys.modules[__name__]
    total, same, aliased, missing = 0, 0, {}, []
    for mod in (GFD, VM, SW):
        for name in mod.__all__:
            total += 1
            obj = getattr(mod, name)
            if getattr(me, name, None) is obj:
                same += 1
            elif mod is SW and name == "potential_vorticity" and getattr(me, "sw_potential_vorticity", None) is obj:
                aliased[f"SW.{name}"] = "sw_potential_vorticity"
            else:
                missing.append(f"{mod.__name__.split('.')[-1]}.{name}")
    return dict(total=total, same_object=same, aliased=aliased, missing=missing)


# =====================================================================================================================
# 1. §13.1–§13.5 helpers
# =====================================================================================================================

_COMPASS = ["north", "north-east", "east", "south-east", "south", "south-west", "west", "north-west"]


def wind_from_to(u, v) -> dict:
    """Name a horizontal flow both ways: as a wind (where it comes from) and as a current (where it goes).

    Book: §13.1 (an "easterly wind" blows from the east; the same flow in the ocean is a "westward current").
    Parameters: u eastward, v northward velocity components [m/s] (scalars).
    Returns dict: ``wind_name`` (e.g. "westerly"), ``current_name`` (e.g. "eastward"), ``from_deg`` and ``to_deg``
    (compass bearings, degrees clockwise from north, in [0, 360)); "calm" and NaN for zero velocity.
    Assumptions: 8-point compass.   Validation: V1 the four cardinal directions.  Label: analytic.
    """
    u, v = float(u), float(v)
    if u == 0.0 and v == 0.0:
        return dict(wind_name="calm", current_name="calm", from_deg=float("nan"), to_deg=float("nan"))
    to = float(np.rad2deg(np.arctan2(u, v)) % 360.0)
    frm = (to + 180.0) % 360.0
    idx = lambda a: int(np.floor((a + 22.5) / 45.0)) % 8  # noqa: E731
    return dict(wind_name=_COMPASS[idx(frm)] + "erly", current_name=_COMPASS[idx(to)] + "ward", from_deg=frm, to_deg=to)


def ocean_density_gradient_budget(drho_dz, rho, c, g: float = G0) -> dict:
    """Split an in-situ density gradient into its adiabatic (compression) part and the potential-density gradient.

    Book: §13.2, Eq. (13.1): dρ_θ/dz = dρ/dz + gρ/c².
    Parameters: drho_dz in-situ gradient [kg/m⁴] (negative: density decreases upward); rho [kg/m³]; c speed of sound
    [m/s]; g [m/s²].
    Returns dict: ``in_situ`` (= drho_dz), ``adiabatic`` (= −gρ/c², the gradient of a neutrally stratified column),
    ``potential`` (= in_situ − adiabatic, Eq. (13.1); negative ⇔ statically stable), ``compression_share``
    (= adiabatic/in_situ: the fraction of the observed increase with depth that is mere compression).
    Assumptions: hydrostatic; c² = (∂p/∂ρ) at constant entropy and salinity.
    Validation: V1 identity with ``core.stratification.ocean_potential_density_gradient``.  Label: analytic.
    """
    a = -g * _F(rho) / _F(c) ** 2
    d = _F(drho_dz)
    return dict(in_situ=_S(d), adiabatic=_S(a), potential=_S(d - a), compression_share=_S(a / d))  # Eq. (13.1)


def atmosphere_layers(z) -> dict:
    """Temperature, lapse rate, static stability and layer name of the standard atmosphere at height z.

    Book: §13.2 (troposphere, tropopause, stratosphere, stratopause); the profile is the public U.S. Standard
    Atmosphere 1976 of ``core.statics`` (ch01), not the book's figure.
    Parameters: z geopotential altitude [m], 0 … 84 852 (scalar or array).
    Returns dict: ``T`` [K]; ``dT_dz`` [K/m] in Kundu's sign (negative where temperature falls with height; the
    meteorological lapse rate is its negative); ``N2`` [1/s²] = (g/T)(dT/dz + g/C_p) for dry air; ``layer`` — one of
    "troposphere", "stratosphere", "mesosphere" (a string, or an array of strings).
    Assumptions: dry perfect gas; the isothermal layer above 11 km and everything up to the stratopause (47 km and
    the isothermal layer to 51 km) is counted as stratosphere.
    Validation: V1 layer table of USSA-1976; N² > 0 at every height.  Label: analytic.
    """
    from .core import statics as _st

    zz = _F(z)
    T = _F(_st.standard_atmosphere(zz)[0])
    idx = np.clip(np.searchsorted(_st.USSA_BASES, zz, side="right") - 1, 0, len(_st.USSA_LAPSE) - 1)
    dT = np.asarray(_st.USSA_LAPSE)[idx]
    names = np.where(zz < 11000.0, "troposphere", np.where(zz < 51000.0, "stratosphere", "mesosphere"))
    return dict(T=_S(T), dT_dz=_S(dT), N2=_S(_F(brunt_vaisala_sq_from_lapse(T, dT))),
                layer=str(names) if names.ndim == 0 else names)


def ocean_profile_idealized(z, *, T_surface: float = 18.0, T_deep: float = 2.0, h_mixed: float = 50.0,
                            h_thermocline: float = 400.0, alpha: float = 2.0e-4, rho0: float = 1027.0, g: float = G0) -> dict:
    """An idealised ocean: a mixed layer over an exponential thermocline, with its analytic buoyancy frequency.

    Book: §13.2 describes the mixed layer, the thermocline and the maximum of N in words and a figure; this profile
    is **ours** (illustrative; not data).
    T = T_surface for z ≥ −h_mixed, T = T_deep + (T_surface − T_deep) exp[(z + h_mixed)/h_thermocline] below;
    ρ_θ = ρ₀ [1 − α (T − T_deep)];  N² = −(g/ρ₀) dρ_θ/dz = g α dT/dz.
    Parameters: z height [m], ≤ 0 (sea surface at 0); T_surface, T_deep [°C]; h_mixed, h_thermocline [m]; alpha
    thermal expansion [1/K]; rho0 density of the deep water [kg/m³]; g [m/s²] — all keyword-only with defaults.
    Returns dict: ``T`` [°C], ``rho_theta`` [kg/m³], ``N2`` [1/s²] (analytic; zero in the mixed layer, largest just
    below it), ``N`` [rad/s].
    Assumptions: density set by temperature alone, linear equation of state.
    Validation: V1 ``buoyancy_frequency_sq`` of ``rho_theta`` reproduces ``N2`` (order 2).  Label: analytic.
    """
    zz = _F(z)
    below = zz < -h_mixed
    e = np.exp(np.minimum(zz + h_mixed, 0.0) / h_thermocline)
    T = np.where(below, T_deep + (T_surface - T_deep) * e, T_surface)
    N2 = np.where(below, g * alpha * (T_surface - T_deep) / h_thermocline * e, 0.0)
    return dict(T=_S(T), rho_theta=_S(rho0 * (1.0 - alpha * (T - T_deep))), N2=_S(N2), N=_S(np.sqrt(N2)))


def thermocline_N2(z, *, N_deep: float = 5.0e-4, N_peak: float = 8.0e-3, z_t: float = -300.0, width: float = 150.0):
    """A smooth thermocline: N²(z) = (N_deep + (N_peak − N_deep) exp(−((z − z_t)/width)²))².

    Book: §13.2 (the buoyancy frequency has a maximum in the thermocline) — the profile is **ours**, for the
    vertical-mode examples of §13.9.
    Parameters: z [m] (≤ 0); N_deep, N_peak [rad/s]; z_t depth of the maximum [m]; width [m] — keyword-only.
    Returns N² [1/s²] (> 0 everywhere, so the mode solvers accept it).
    Assumptions: none (a shape).
    Validation (measured, defaults, H = 4200 m): ``vertical_modes`` on 1601 nodes gives c₁ = 1.7715, c₂ = 0.6179,
    c₃ = 0.4659 m/s (c₁ is 0.17 % lower on 401 nodes; shooting agrees to 1e-4); the WKB estimate
    ``wkb_mode_speed`` gives 1.302, 0.651, 0.434 m/s — off by −27 %, +5 % and −7 %: WKB is not reliable for the
    gravest mode of a sharp thermocline.  Label: analytic.
    """
    zz = _F(z)
    return _S((N_deep + (N_peak - N_deep) * np.exp(-((zz - z_t) / width) ** 2)) ** 2)


def anisotropic_eddy_stress(G, nu_H, nu_v, rho):
    """Turbulent stress tensor of the chapter's anisotropic eddy-viscosity model.

    Book: §13.3, Eq. (13.5): τ_xz = ρν_v ∂u/∂z + ρν_H ∂w/∂x, τ_yz = ρν_v ∂v/∂z + ρν_H ∂w/∂y,
    τ_xy = ρν_H (∂u/∂y + ∂v/∂x), τ_xx = 2ρν_H ∂u/∂x, τ_yy = 2ρν_H ∂v/∂y, τ_zz = 2ρν_v ∂w/∂z.
    Parameters: G velocity-gradient tensor, G[i, j] = ∂u_i/∂x_j [1/s] (3 × 3); nu_H, nu_v eddy viscosities [m²/s];
    rho [kg/m³].
    Returns τ (3 × 3, symmetric) [Pa].
    Assumptions: the book's model, which is **not frame-indifferent**: a rigid rotation in a vertical plane
    (u = ωz, w = −ωx) gives τ_xz = ρω(ν_v − ν_H) ≠ 0.  With ν_H = ν_v it reduces to the Newtonian 2ρνS_ij.
    Validation: V7 both statements above.  Label: analytic.
    """
    G = _F(G)
    if G.shape != (3, 3):
        raise ValueError("anisotropic_eddy_stress: G must be 3 x 3 with G[i, j] = du_i/dx_j")
    t = np.zeros((3, 3))
    t[0, 2] = t[2, 0] = rho * nu_v * G[0, 2] + rho * nu_H * G[2, 0]   # Eq. (13.5), tau_xz
    t[1, 2] = t[2, 1] = rho * nu_v * G[1, 2] + rho * nu_H * G[2, 1]   # tau_yz
    t[0, 1] = t[1, 0] = rho * nu_H * (G[0, 1] + G[1, 0])              # tau_xy
    t[0, 0] = 2.0 * rho * nu_H * G[0, 0]
    t[1, 1] = 2.0 * rho * nu_H * G[1, 1]
    t[2, 2] = 2.0 * rho * nu_v * G[2, 2]
    return t


def friction_force_dimensions(printed: bool = False) -> dict:
    """Dimensions (exponents of mass, length, time) of the friction force built from the stress.

    Book: §13.3 — the text calls F_i = ∂τ_ij/∂x_j a force per unit mass (slip #2); the corrected form divides by ρ.
    Parameters: printed — if True, return the dimensions of the printed form (no 1/ρ).
    Returns dict: ``M``, ``L``, ``T`` exponents of F_i and ``is_acceleration`` (True for the corrected form:
    L T⁻²; False for the printed one: M L⁻² T⁻², a force per unit volume).
    Assumptions: stress in Pa = M L⁻¹ T⁻².   Validation: V2 this is the planted variant of slip #2.
    Label: analytic.
    """
    M, L, T = 1, -1 - 1, -2          # d(tau)/dx
    if not printed:
        M, L = M - 1, L + 3          # divided by rho = M L^-3
    return dict(M=M, L=L, T=T, is_acceleration=(M, L, T) == (0, 1, -2))


def term_table_thin_layer(U, L, H, lat_rad, **kw):
    """Sizes of every term of the thin-shell momentum equations for a choice of scales, as a table.

    Book: §13.4, Eq. (13.9); built on :func:`fluidpy.core.gfd.thin_layer_terms` (same arguments and keywords).
    Returns a ``pandas.DataFrame``: rows "acceleration", "coriolis", "coriolis_w", "pressure", "buoyancy",
    "friction_H", "friction_v"; columns "x" and "z" (sizes [m/s²]; NaN where the term does not appear) and
    "x_over_coriolis", "z_over_coriolis" (ratios to the horizontal Coriolis term |f|U — the first is the Rossby
    number in the "acceleration" row).
    Assumptions: as ``thin_layer_terms``.   Validation: V1 the ratio in the acceleration row equals Ro.
    Label: analytic.
    """
    import pandas as pd

    t = thin_layer_terms(U, L, H, lat_rad, **kw)
    names = ["acceleration", "coriolis", "coriolis_w", "pressure", "buoyancy", "friction_H", "friction_v"]
    cor = t["x"]["coriolis"]
    x = [t["x"].get(n, float("nan")) for n in names]
    z = [t["z"].get(n, float("nan")) for n in names]
    return pd.DataFrame(dict(x=x, z=z, x_over_coriolis=[v / cor for v in x], z_over_coriolis=[v / cor for v in z]),
                        index=names)


def pressure_centre(x, y, *, dp, R_c, lat_rad, rho0: float = 1.2, kind: str = "low") -> dict:
    """A Gaussian pressure centre and the geostrophic wind round it.

    Book: §13.5, Eqs. (13.11)–(13.12) (and the figure of circular flow round a low and a high).
    p′ = ∓dp exp(−r²/2R_c²) (− for a low), u = −(1/ρ₀f) ∂p/∂y, v = (1/ρ₀f) ∂p/∂x (analytic derivatives).
    Parameters: x, y coordinates [m] — 1-D arrays (the fields are then ``[j, i]`` = (y, x)) or scalars; dp central
    pressure anomaly [Pa] > 0; R_c radius [m]; lat_rad latitude [rad] (not the equator); rho0 [kg/m³]; kind "low" or
    "high" — all physical parameters keyword-only.
    Returns dict: ``p`` [Pa], ``u``, ``v`` [m/s], ``f`` [1/s].  The flow is counter-clockwise round a low for f > 0
    and clockwise for f < 0; the wind speed peaks at r = R_c.
    Assumptions: geostrophic (Ro ≪ 1), no friction.   Validation: V1 u·∇p = 0; V7 the hemisphere flips the sense.
    Label: analytic.
    """
    if kind not in ("low", "high"):
        raise ValueError('kind must be "low" or "high"')
    xx, yy = _F(x), _F(y)
    if xx.ndim == 1 and yy.ndim == 1:
        X, Y = np.meshgrid(xx, yy)
    else:
        X, Y = np.broadcast_arrays(xx, yy)
    f = float(coriolis_parameter(lat_rad))
    p = (-1.0 if kind == "low" else 1.0) * dp * np.exp(-(X ** 2 + Y ** 2) / (2.0 * R_c ** 2))
    u, v = geostrophic_velocity(-p * X / R_c ** 2, -p * Y / R_c ** 2, f, rho0)
    return dict(p=_S(p), u=_S(u), v=_S(v), f=f)


def jet_section(y, z, *, dT, width, lat_rad, alpha, T0: float = 280.0, lapse: float = -6.5e-3, u_surface: float = 0.0,
                g: float = G0) -> dict:
    """A frontal zone and the jet that balances it: temperature and zonal wind in exact thermal-wind balance.

    Book: §13.5, Eq. (13.15) in its temperature form (**ours — not in the book** as a formula):
    T(y, z) = T0 + lapse·z − (dT/2) tanh(y/width)·sign(latitude),   ∂U/∂z = −(g α/f) ∂T/∂y.
    The sign factor puts the cold air on the poleward side in either hemisphere, so the jet is westerly in both.
    Parameters: y northward coordinate [m] (1-D), z height [m] (1-D); dT temperature contrast across the front [K];
    width [m]; lat_rad [rad]; alpha thermal expansion [1/K] (required: 1/T0 for a perfect gas); T0 [K]; lapse dT/dz
    in Kundu's sign [K/m]; u_surface wind at z = 0 [m/s] — physical parameters keyword-only.
    Returns dict: ``T`` (z, y) [K], ``U`` (z, y) [m/s], ``dTdy`` (y,) [K/m].
    Assumptions: geostrophic and hydrostatic; ∂T/∂y independent of height (so U is linear in z); Boussinesq.
    Validation: V1 ∂U/∂z equals ``thermal_wind_from_temperature``; V7 westerly in both hemispheres.  Label: analytic.
    """
    yy, zz = _F(y), _F(z)
    s = float(np.sign(np.sin(float(lat_rad))))
    f = float(coriolis_parameter(lat_rad))
    if f == 0.0:
        raise ValueError("jet_section: no thermal-wind balance at the equator (f = 0)")
    dTdy = -0.5 * dT * s / width / np.cosh(yy / width) ** 2
    T = T0 + lapse * zz[:, None] - 0.5 * dT * s * np.tanh(yy / width)[None, :]
    U = u_surface + (-g * float(alpha) / f) * dTdy[None, :] * zz[:, None]
    return dict(T=T, U=U, dTdy=dTdy)


def viscous_layer_thicknesses(nu, *, t=None, x=None, U=None, f=None) -> dict:
    """Three ways a viscous layer gets its thickness: growing in time, growing downstream, or held by rotation.

    Book: §13.6 (opening): a balance of friction with ∂u/∂t (ch08), with u·∇u (ch09) or with the Coriolis force
    (the Ekman layer, Eq. (13.29)).
    Parameters: nu viscosity [m²/s]; t time [s]; x distance downstream [m] with U free-stream speed [m/s]; f Coriolis
    parameter [1/s] — keyword-only, each optional.
    Returns dict: ``diffusive`` = sqrt(νt), ``boundary_layer`` = sqrt(νx/U), ``ekman`` = sqrt(2ν/|f|) [m]; ``None``
    for those whose inputs are missing.  (Scales, without the order-one factors of the exact solutions.)
    Assumptions: laminar or constant eddy viscosity.   Validation: V1 identities.  Label: analytic.
    """
    return dict(diffusive=None if t is None else float(np.sqrt(nu * t)),
                boundary_layer=None if (x is None or U is None) else float(np.sqrt(nu * x / U)),
                ekman=None if f is None else float(ekman_depth(nu, f)))


def ekman_surface_printed(z, tau, rho, nu_v, f):
    """The surface Ekman spiral exactly as the book prints it (northern hemisphere, stress along x).

    Book: §13.6, the unnumbered solution after Eq. (13.29):
    u = (τ/ρ)/sqrt(f ν_v) · e^{z/δ} cos(−z/δ + π/4),   v = −(τ/ρ)/sqrt(f ν_v) · e^{z/δ} sin(−z/δ + π/4).
    Parameters: z ≤ 0 [m]; tau stress along x [N/m²]; rho [kg/m³]; nu_v [m²/s]; f [1/s] **> 0 only**.
    Returns (u, v) [m/s].   Raises ValueError for f ≤ 0 (the printed form takes sqrt(f); use ``ekman_surface``).
    Assumptions: as ``ekman_surface``.   Validation: V1 equals ``ekman_surface`` for f > 0.  Label: analytic.
    """
    if float(f) <= 0:
        raise ValueError("ekman_surface_printed: the printed cos/sin form holds for f > 0 only; use ekman_surface")
    zz = _F(z)
    d = float(ekman_depth(nu_v, f))
    a = tau / rho / np.sqrt(f * nu_v) * np.exp(zz / d)
    return _S(a * np.cos(-zz / d + np.pi / 4.0)), _S(-a * np.sin(-zz / d + np.pi / 4.0))


def coastal_upwelling(tau_alongshore, coast_side: str, lat_rad, rho: float = 1027.0) -> dict:
    """Does an alongshore wind drive surface water away from a north–south coast (upwelling)?

    Book: §13.12 (closing paragraph: the Ekman transport of Eq. (13.30) moves surface water offshore, deeper water
    rises, and the disturbance of the thermocline travels along the coast as an internal Kelvin wave).
    Parameters: tau_alongshore northward wind stress [N/m²] (negative = southward); coast_side "east" or "west" —
    the side of the water on which the land lies ("east" for the west coast of a continent); lat_rad [rad];
    rho [kg/m³].
    Returns dict: ``transport_offshore`` [m²/s] (Ekman transport per unit length of coast, positive away from the
    land), ``upwelling`` (bool: transport_offshore > 0), ``kelvin_direction`` ("poleward" when the land is on the
    east, "equatorward" when it is on the west — the same in both hemispheres, because the trapped wave keeps the
    coast on its right for f > 0 and on its left for f < 0).
    Assumptions: steady Ekman layer, straight meridional coast, deep water offshore.
    Validation: V7 an equatorward wind on an eastern-boundary coast upwells in both hemispheres.  Label: analytic.
    """
    if coast_side not in ("east", "west"):
        raise ValueError('coast_side must be "east" or "west"')
    Mx, _ = ekman_transport(0.0, tau_alongshore, rho, coriolis_parameter(lat_rad))
    off = -Mx if coast_side == "east" else Mx
    return dict(transport_offshore=float(off), upwelling=bool(off > 0),
                kelvin_direction="poleward" if coast_side == "east" else "equatorward")


# =====================================================================================================================
# 2. §13.13–§13.17 helpers
# =====================================================================================================================

def flow_over_step(x, U, beta, f0, h0, h1) -> dict:
    """Streamline of a uniform zonal flow crossing a step in depth at x = 0 on a β-plane (linearised).

    Book: §13.13 describes both cases in words and sketches; the quantitative solution is **ours — not in the book**.
    A column that crosses from depth h₀ to h₁ keeps (ζ + f)/h (Eq. (13.94)).  With the streamline displacement Y(x),
    v = U Y′, ζ = U Y″ and f = f₀ + βY, conservation gives upstream of the step U Y″ + βY = 0 and downstream
    U Y″ + βY = f₀ (h₁ − h₀)/h₀, with Y and Y′ continuous at x = 0.  Writing Y_p = f₀ (h₁ − h₀)/(β h₀), κ = sqrt(β/|U|):
    * **eastward, U > 0** (upstream x < 0, depth h₀): Y = 0 for x < 0, Y = Y_p (1 − cos κx) for x > 0 — a standing
      Rossby wave of wavelength 2π sqrt(U/β);
    * **westward, U < 0** (upstream x > 0, depth h₀): Y = ½Y_p e^{−κx} for x > 0, Y = Y_p (1 − ½ e^{κx}) for x < 0 —
      no oscillation, and the flow feels the step before it reaches it.
    Parameters: x [m]; U zonal speed [m/s] (sign = direction, non-zero); beta [1/(m s)] > 0; f0 [1/s]; h0 upstream
    depth, h1 downstream depth [m].
    Returns dict: ``Y`` [m] (northward displacement of the streamline that is at y = 0 far upstream), ``zeta`` = U Y″
    [1/s], ``wavelength`` [m] (2π sqrt(U/β) for U > 0, None for U < 0), ``decay_length`` = sqrt(|U|/β) [m].
    Far downstream ζ → 0 requires f₀ + βY = f₀ h₁/h₀: the streamline ends (U < 0) or oscillates about (U > 0) the
    latitude shifted by Y_p.
    Assumptions: |βY| ≪ |f₀|, |ζ| ≪ |f₀| (small step), steady, speed U unchanged across the step, uniform in y.
    Validation: V1 ζ just past the step equals ``step_vorticity``; V1 wavelength equals
    ``stationary_rossby_wavelength``; V7 westward case monotonic.  Label: analytic (ours).
    """
    xx = _F(x)
    U, beta = float(U), float(beta)
    if U == 0 or beta <= 0:
        raise ValueError("flow_over_step: needs U != 0 and beta > 0")
    Yp = f0 * (h1 - h0) / (beta * h0)
    kap = np.sqrt(beta / abs(U))
    if U > 0:
        down = xx > 0
        Y = np.where(down, Yp * (1.0 - np.cos(kap * xx)), 0.0)
        Ypp = np.where(down, Yp * kap ** 2 * np.cos(kap * xx), 0.0)
        lam = float(2.0 * np.pi / kap)
    else:
        up = xx > 0
        e = np.exp(-kap * np.abs(xx))
        Y = np.where(up, 0.5 * Yp * e, Yp * (1.0 - 0.5 * e))
        Ypp = np.where(up, 0.5 * Yp * kap ** 2 * e, -0.5 * Yp * kap ** 2 * e)
        lam = None
    return dict(Y=_S(Y), zeta=_S(U * Ypp), wavelength=lam, decay_length=float(1.0 / kap))


def vertical_structure_solve(z, N_fn: Callable, k, omega, f, w0: complex = 0.0, dw0: complex = 1.0):
    """Numerical solution of the vertical-structure equation of internal waves, d²ŵ/dz² + m²(z) ŵ = 0.

    Book: §13.14, Eq. (13.100) with m²(z) of Eq. (13.99) (l = 0).
    Parameters: z levels [m], increasing (output points; integration starts at z[0]); N_fn callable z → N [rad/s];
    k [rad/m]; omega [rad/s]; f [1/s]; w0, dw0 — ŵ and dŵ/dz at z[0] (real or complex).
    Returns ŵ at the levels (complex array if the initial values are complex, else real).
    Numerics: ``scipy.integrate.solve_ivp`` (DOP853, rtol 1e-10, atol 1e-12 of the solution scale).
    Assumptions: linear, Boussinesq, f-plane.   Validation: V1 uniform N gives a sinusoid of wavenumber m.
    Label: converged.
    """
    from scipy.integrate import solve_ivp

    zz = _F(z)
    cplx = np.iscomplexobj(w0) or np.iscomplexobj(dw0)
    y0 = np.array([w0, dw0], dtype=complex if cplx else float)

    def rhs(s, y):
        return [y[1], -float(inertia_gravity_m2(k, 0.0, omega, float(N_fn(s)), f)) * y[0]]

    sol = solve_ivp(rhs, (zz[0], zz[-1]), y0, t_eval=zz, method="DOP853", rtol=1e-10,
                    atol=1e-12 * max(abs(w0), abs(dw0) * abs(zz[-1] - zz[0]), 1e-300))
    if not sol.success:
        raise RuntimeError("vertical_structure_solve: integration failed: " + sol.message)
    return sol.y[0]


def wkb_error(z, N_fn: Callable, k, omega, f) -> dict:
    """How good is the WKB solution for a given N(z)?  Compare it with the numerical solution started from the same
    values.

    Book: §13.14, Eqs. (13.101)–(13.104) and the condition H m ≫ 1.
    Parameters: z levels [m], increasing, fine enough to resolve the wave (≥ 20 per vertical wavelength); N_fn
    callable z → N [rad/s]; k [rad/m]; omega [rad/s] with |f| < ω < N everywhere on z; f [1/s].
    Returns dict: ``max_rel_error`` — the largest |ŵ_numerical − ŵ_WKB| over z divided by the largest |ŵ_WKB|, the
    numerical solution of Eq. (13.100) being started with the WKB value and slope at z[0] (upper sign); ``Hm`` — the
    smallest value on z of m·N/|dN/dz| (the book's H m with H = N/|dN/dz|; ``inf`` for uniform N).
    Assumptions: as ``wkb_vertical_structure``.
    Validation (measured; N = N₀ (1 + 0.5 z/D) on −D ≤ z ≤ 0 with N₀ = 2.7e-3 rad/s, k = 1e-3 rad/m,
    ω = 4e-4 rad/s, f at 35°; 4001 levels): D = 0.5, 1, 2, 4 km give Hm = 1.6, 3.3, 6.6, 13.2 and errors 0.104,
    0.045, 0.020, 0.0094 — the error falls in proportion to 1/Hm (about 0.15/Hm here); uniform N gives 2e-10
    (the WKB form is then exact).  Label: approximate (WKB), error measured.
    """
    zz = _F(z)
    Nz = np.array([float(N_fn(v)) for v in zz])
    m2 = _F(inertia_gravity_m2(k, 0.0, omega, Nz, f))
    if np.any(m2 <= 0):
        raise ValueError("wkb_error: the wave does not propagate everywhere on z (m^2 <= 0 somewhere)")
    m = np.sqrt(m2)
    wk = GFD.wkb_vertical_structure(zz, m, 1.0, +1)
    dm0 = (-3.0 * m[0] + 4.0 * m[1] - m[2]) / (zz[2] - zz[0])
    num = vertical_structure_solve(zz, N_fn, k, omega, f, w0=wk[0], dw0=(1j * m[0] - 0.5 * dm0 / m[0]) * wk[0])
    dN = np.gradient(Nz, zz)
    with np.errstate(divide="ignore"):
        Hm = np.where(dN != 0, m * Nz / np.abs(np.where(dN != 0, dN, 1.0)), np.inf)
    return dict(max_rel_error=float(np.max(np.abs(num - wk)) / np.max(np.abs(wk))), Hm=float(np.min(Hm)))


def w_equation_rotating_residual(w_fn: Callable, x, y, z, t, N, f, h: float = 1.0, ht: float = 1.0) -> float:
    """Residual of the w-equation of rotating internal waves for any trial field w(x, y, z, t).

    Book: §13.14, Eq. (13.96): ∂²/∂t² ∇²w + N² ∇_H² w + f² ∂²w/∂z² = 0.
    Parameters: w_fn callable (x, y, z, t) → w; the point x, y, z [m], t [s]; N [rad/s]; f [1/s]; h difference step
    in space [m]; ht in time [s].
    Returns the residual [1/(m s³)] by nested second-order central differences (second differences in space inside a
    second difference in time).  Compare it with the size of one term, e.g. N²K_h²|w|: choose h and ht a few
    hundredths of the wavelength and period (truncation ∝ h², ht²; round-off ∝ 1e-16/(h² ht²)).
    Assumptions: uniform N.   Validation: V1 zero (to truncation) for a plane wave obeying Eq. (13.112), non-zero for
    a wrong frequency.  Label: converged.
    """
    def d2(fn, i, p, step):
        a, b = list(p), list(p)
        a[i] += step
        b[i] -= step
        return (fn(*a) - 2.0 * fn(*p) + fn(*b)) / step ** 2

    def lap(*p):
        return sum(d2(w_fn, i, p, h) for i in (0, 1, 2))

    p0 = (float(x), float(y), float(z), float(t))
    return float(d2(lap, 3, p0, ht) + N ** 2 * (d2(w_fn, 0, p0, h) + d2(w_fn, 1, p0, h)) + f ** 2 * d2(w_fn, 2, p0, h))


def lee_wave_field(x, z, U, N, k, h0) -> dict:
    """Linear flow of a uniform stratified stream over sinusoidal topography h = h₀ cos kx (lee waves).

    Book: §13.14 "Lee Wave" (Eq. (13.113) with ω = kU; the phase lines of a stationary wave tilt upstream with
    height because its energy must travel upward from the ground).  The field itself is ours.
    Streamline displacement δ = h₀ cos(kx + mz) for k < N/|U| (m = sqrt(N²/U² − k²)), δ = h₀ cos(kx) e^{−μz} with
    μ = sqrt(k² − N²/U²) otherwise; stream function of the total flow ψ = U (z − δ) (u = ∂ψ/∂z, w = −∂ψ/∂x — the
    ch04 sign convention), so w = U ∂δ/∂x.
    Parameters: x [m] (1-D or scalar), z height above the mean ground [m] (1-D or scalar); U > 0 mean wind [m/s];
    N [rad/s]; k [rad/m]; h0 [m].
    Returns dict: ``psi`` (z, x) [m²/s], ``w`` (z, x) [m/s], ``m`` [rad/m] (None when evanescent), ``tilt``
    ("upstream" — lines of constant phase x = −(m/k) z + const lean against the wind — or "none (evanescent)").
    Assumptions: linear (k h₀, m h₀ ≪ 1), uniform U and N, no rotation, radiation condition aloft.
    Validation: V1 w = U ∂h/∂x at z = 0; V7 the phase tilt has the sign of −U.  Label: analytic.
    """
    xx, zz = _F(x), _F(z)
    X, Z = (np.meshgrid(xx, zz) if (xx.ndim == 1 and zz.ndim == 1) else np.broadcast_arrays(xx, zz))
    try:
        m = float(lee_wave_m(U, N, k))
        delta = h0 * np.cos(k * X + m * Z)
        ddx = -h0 * k * np.sin(k * X + m * Z)
        tilt = "upstream"
    except ValueError:
        m = None
        mu = float(np.sqrt(k ** 2 - (N / U) ** 2))
        delta = h0 * np.cos(k * X) * np.exp(-mu * Z)
        ddx = -h0 * k * np.sin(k * X) * np.exp(-mu * Z)
        tilt = "none (evanescent)"
    return dict(psi=_S(U * (Z - delta)), w=_S(U * ddx), m=m, tilt=tilt)


def basin_crossing_time(L, lat_rad, c):
    """Time a long Rossby wave takes to cross a basin of width L: L / (β c²/f₀²).

    Book: §13.15 (the long-wave speed c_x ≈ −βc²/f₀²; "years to cross the ocean at mid-latitudes").
    Parameters: L width [m]; lat_rad latitude [rad] (not the equator; the formula fails within a few degrees of it);
    c long-wave speed of the mode [m/s].
    Returns the time [s].
    Assumptions: long (non-dispersive) waves, β and f₀ evaluated at one latitude.
    Validation: V1 identity with ``rossby_long_wave_speed``.  Label: analytic.
    """
    return _S(_F(L) / np.abs(_F(rossby_long_wave_speed(beta_parameter(lat_rad), coriolis_parameter(lat_rad), c))))


def rayleigh_kuo_eigs(k, U: Callable, Up: Callable, Upp: Callable, beta, domain=(-1.0, 1.0), N: int = 120,
                      delta: float = 0.2, **kw) -> dict:
    """Leading eigenvalue of the barotropic stability problem on a β-plane.

    Book: §13.16, the normal-mode form of Eq. (13.123): (U − c)(ψ̂″ − k²ψ̂) + (β − U″) ψ̂ = 0 — Rayleigh's equation
    (11.81) with U″ replaced by U″ − β.  Solved by ``core.stability.rayleigh_eigs_contour`` (ch11; Chebyshev
    collocation on a complex path) through its ``beta`` keyword.  Here ψ is defined by u′ = −∂ψ/∂y, v′ = ∂ψ/∂x; the
    eigenvalue c does not depend on that sign choice.
    Parameters: k wavenumber (non-dimensional, scaled by 1/L); U, Up, Upp callables y ↦ U, U′, U″ (U and Upp
    analytic, evaluated at complex y); beta (scaled by U₀/L²); domain, N, delta and any further keyword of
    ``rayleigh_eigs_contour`` (``bc="decay"``, ``y_max``, ``parity`` …).
    Returns dict: ``c`` (complex leading eigenvalue; ``nan`` when no growing mode is found — meaning "stable or
    unresolved", see the limits of the solver), ``growth_rate`` = k·Im c (0.0 then), ``stable`` (bool).
    Assumptions: inviscid, parallel, barotropic; the eigen-solver contract of ch11 (leading mode only; check box and
    N doubling for anything quoted).
    Validation (measured; jet U = sech²y, even mode, ``bc="decay"``, y_max = 16, N = 120): growth rate k·c_i at
    k = 0.9 is 0.1608, 0.1373, 0.1064, 0.0670 for β = 0, 0.1, 0.2, 0.3 and no growing mode at 0.45; at k = 1.4 it
    falls from 0.1191 (β = 0) to 0.0151 (β = 0.6) and no growing mode is found at β = 0.7 — beyond max U″ = 2/3,
    as the Rayleigh–Kuo condition requires; the values change by less than 1e-9 with y_max = 24, N = 160.  Not
    converged at default resolution for strongly negative β at small k (k = 0.9, β = −1: 0.055 against 0.060), where
    the critical level approaches the jet maximum (limit (iii) of the solver).  V7 β = 0 reproduces ch11.
    Label: converged (inherits ``rayleigh_eigs_contour``; ranges above).
    """
    from .core.stability import rayleigh_eigs_contour

    # beta goes in as its own argument: folding it into Upp would corrupt the complex path (which uses U'' itself)
    c = rayleigh_eigs_contour(k, U, Up, Upp, domain=domain, N=N, delta=delta, beta=float(beta), **kw)
    if len(c) == 0:
        return dict(c=complex("nan"), growth_rate=0.0, stable=True)
    return dict(c=complex(c[0]), growth_rate=float(k * c[0].imag), stable=False)


def eady_basic_state(y, z, *, N, f, H, U0, rho0, g: float = G0) -> dict:
    """Basic state of the Eady problem: uniform stratification, uniform shear, and the sloping density surfaces that
    the thermal wind requires.

    Book: §13.17, Eqs. (13.127)–(13.128): dU/dz = (g/fρ₀) ∂ρ̄/∂y with U = U₀ z/H.
    ρ̄(y, z) = ρ₀ (1 − N² z/g) + (f ρ₀ U₀/(g H)) y.
    Parameters: y [m] (1-D), z height above the lower lid [m] (1-D); N [rad/s]; f [1/s] non-zero; H [m]; U0 [m/s];
    rho0 [kg/m³]; g — physical parameters keyword-only.
    Returns dict: ``rho`` (z, y) [kg/m³], ``U`` (z,) [m/s], ``slope`` = f U₀/(N² H) (dz/dy of a density surface;
    positive = rising toward +y for f U₀ > 0), ``drho_dy`` [kg/m⁴].
    Assumptions: Boussinesq, geostrophic and hydrostatic basic state.
    Validation: V1 ``thermal_wind_shear`` of ``drho_dy`` gives U₀/H.  Label: analytic.
    """
    yy, zz = _F(y), _F(z)
    fa = float(f)
    if fa == 0:
        raise ValueError("eady_basic_state: f = 0")
    drho_dy = fa * rho0 * U0 / (g * H)   # Eq. (13.128)
    rho = rho0 * (1.0 - N ** 2 * zz[:, None] / g) + drho_dy * yy[None, :]
    return dict(rho=rho, U=U0 * zz / H, slope=fa * U0 / (N ** 2 * H), drho_dy=drho_dy)


def eady_matrix(c, alphaH, U0, H):
    """The 2 × 2 matrix of the Eady problem: coefficients of (A, B) in the two lid conditions.

    Book: §13.17, the pair of homogeneous equations before Eq. (13.141), with q = αH/2 and α = alphaH/H:
    row 1 (z = 0): A[αc sinh q − (U₀/H) cosh q] + B[−αc cosh q + (U₀/H) sinh q] = 0;
    row 2 (z = H): A[α(U₀ − c) sinh q − (U₀/H) cosh q] + B[α(U₀ − c) cosh q − (U₀/H) sinh q] = 0.
    Parameters: c trial phase speed [m/s] (complex allowed); alphaH (dimensionless); U0 [m/s]; H [m].
    Returns a complex array (2, 2) [m/s·1/m].  Its determinant vanishes at the two roots of Eq. (13.141).
    Assumptions: Eady problem.   Validation: V2 |det| at ``eady_phase_speed`` is round-off.  Label: analytic.
    """
    a, q = alphaH / H, 0.5 * alphaH
    c = complex(c)
    return np.array([[a * c * np.sinh(q) - U0 / H * np.cosh(q), -a * c * np.cosh(q) + U0 / H * np.sinh(q)],
                     [a * (U0 - c) * np.sinh(q) - U0 / H * np.cosh(q), a * (U0 - c) * np.cosh(q) - U0 / H * np.sinh(q)]],
                    dtype=complex)


def eady_numeric_eigs(k, l, N, f, H, U0, n: int = 64) -> dict:
    """Eady eigenvalue by Chebyshev collocation — an independent numerical route to Eq. (13.141).

    Book: §13.17, Eq. (13.136) with normal modes (13.137): (U − c)(p̂″ − α²p̂) = 0 on 0 < z < H, with the lid
    conditions c p̂′ + (U₀/H) p̂ = 0 at z = 0 and (U₀ − c) p̂′ − (U₀/H) p̂ = 0 at z = H (from Eq. (13.135)).
    Written as the generalised eigenproblem A p = c B p with interior rows U (D² − α²), (D² − α²) and the two
    boundary rows; solved with ``scipy.linalg.eig``.  The interior rows carry a continuous spectrum c = U(z_j)
    (singular modes); the two Eady modes are the eigenvectors with p̂″ − α²p̂ = 0, picked by that residual.
    Parameters: k, l [rad/m]; N [rad/s]; f [1/s]; H [m]; U0 [m/s]; n polynomial degree.
    Returns dict: ``c`` (complex [m/s]: the one of the two Eady modes with the larger imaginary part — for a neutral
    pair the faster) and ``growth_rate`` = |k|·Im c [1/s] (≥ 0).
    Assumptions: Eady problem.
    Validation (measured, our atmosphere inputs, n = 64): agrees with ``eady_phase_speed`` to 6e-11 of U₀ at
    αH = 0.5, 1.6, 2.2 and 3.0.  Label: converged.
    """
    from scipy.linalg import eig

    from .core.stability import cheb

    al = float(GFD.eady_alpha(k, l, N, f))
    D, z = cheb(int(n), (0.0, H))                    # z descending: z[0] = H, z[-1] = 0
    Lop = D @ D - al ** 2 * np.eye(n + 1)
    A = np.diag(U0 * z / H) @ Lop
    B = Lop.copy()
    A[-1], B[-1] = U0 / H * np.eye(n + 1)[-1], -D[-1]                 # z = 0: (U0/H) p = -c p'
    A[0], B[0] = U0 * D[0] - U0 / H * np.eye(n + 1)[0], D[0]          # z = H: U0 p' - (U0/H) p = c p'
    w, V = eig(A, B)
    ok = np.isfinite(w)
    w, V = w[ok], V[:, ok]
    res = np.linalg.norm(Lop[1:-1] @ V, axis=0) / np.maximum(np.linalg.norm(V, axis=0), 1e-300)
    pick = np.argsort(res)[:2]                        # the two modes with p'' - alpha^2 p = 0
    c = w[pick]
    c = c[np.lexsort((c.real, c.imag))][-1]           # larger Im; for a neutral pair the larger Re
    return dict(c=complex(c), growth_rate=float(abs(k) * max(c.imag, 0.0)))


# =====================================================================================================================
# 3. Sympy engines (each returns at least "residual" — 0 when the result holds — and "ok")
# =====================================================================================================================

def boussinesq_rotating_sympy(printed: bool = False) -> dict:
    """Does a fluid at rest satisfy the rotating Boussinesq momentum equation?

    Book: §13.3, Eq. (13.2) with the hydrostatic state of rest (13.3).  As printed, the pressure term of (13.2) has a
    plus sign (slip #1); the corrected equation is Du/Dt + 2Ω × u = −(1/ρ₀)∇p − (gρ/ρ₀) e_z + F.
    Parameters: printed — build the printed (+) form instead.
    Returns dict: ``residual`` — the vertical component of the right-hand side for u = 0, p = p̄(z), ρ = ρ̄(z) with
    dp̄/dz = −ρ̄g (0 for the corrected form; −2gρ̄/ρ₀ for the printed one); ``ok``; ``coriolis_work`` — u·(2Ω × u)
    (identically 0: the Coriolis force does no work); ``equation`` (text).
    Assumptions: Boussinesq.   Validation: V2 this is the check; ``printed=True`` must fail.  Label: symbolic.
    """
    import sympy as sp

    z, g, rho0 = sp.symbols("z g rho_0", positive=True)
    rhobar = sp.Function("rhobar")(z)
    dpbar = -rhobar * g                                   # Eq. (13.3)
    sign = 1 if printed else -1
    res = sp.simplify(sign * dpbar / rho0 - g * rhobar / rho0)   # Eq. (13.2b), vertical component, u = 0
    u = sp.Matrix(sp.symbols("u v w", real=True))
    Om = sp.Matrix(sp.symbols("Omega_x Omega_y Omega_z", real=True))
    work = sp.simplify(u.dot((2 * Om).cross(u)))
    return dict(residual=res, ok=bool(res == 0), coriolis_work=work,
                equation="Du/Dt + 2 Omega x u = " + ("+" if printed else "-") + "(1/rho0) grad p - (g rho/rho0) e_z + F")


def perturbation_form_sympy() -> dict:
    """The pressure and gravity terms keep their form when p and ρ are replaced by perturbations from rest.

    Book: §13.3, Eq. (13.4) and the unnumbered identity after it:
    −(1/ρ₀)∇p − (gρ/ρ₀) e_z = −(1/ρ₀)∇p′ − (gρ′/ρ₀) e_z, using dp̄/dz = −ρ̄g (13.3).
    Returns dict: ``residual`` (sum of squares of the three components of left − right; 0), ``ok``.
    Assumptions: Boussinesq.   Validation: V2.  Label: symbolic.
    """
    import sympy as sp

    x, y, z, t = sp.symbols("x y z t", real=True)
    g, rho0 = sp.symbols("g rho_0", positive=True)
    rhobar = sp.Function("rhobar")(z)
    pbar = sp.Function("pbar")(z)
    pp, rp = sp.Function("pp")(x, y, z, t), sp.Function("rhop")(x, y, z, t)
    p, rho = pbar + pp, rhobar + rp                      # Eq. (13.4)
    left = [-p.diff(x) / rho0, -p.diff(y) / rho0, -p.diff(z) / rho0 - g * rho / rho0]
    right = [-pp.diff(x) / rho0, -pp.diff(y) / rho0, -pp.diff(z) / rho0 - g * rp / rho0]
    diff = [sp.simplify((a - b).subs(pbar.diff(z), -rhobar * g)) for a, b in zip(left, right)]
    res = sp.simplify(sum(d ** 2 for d in diff))
    return dict(residual=res, ok=bool(res == 0), components=tuple(diff))


def eddy_friction_force_sympy() -> dict:
    """The friction force follows from the anisotropic stress — once continuity is used.

    Book: §13.3: (1/ρ) ∂τ_ij/∂x_j with the stresses of Eq. (13.5) equals Eq. (13.6) plus a term ν ∂(∇·u)/∂x_i that
    vanishes for ∇·u = 0 (ν_H for the x and y components, ν_v for z) — the step the book takes silently.
    Returns dict: ``residual`` (sum of squares of the three differences F_i − [(13.6)_i + ν ∂_i(∇·u)]; 0), ``ok``,
    ``extra_terms`` (the three terms proportional to derivatives of ∇·u).
    Assumptions: uniform ρ, ν_H, ν_v.   Validation: V2.  Label: symbolic.
    """
    import sympy as sp

    x, y, z = sp.symbols("x y z", real=True)
    nH, nv, rho = sp.symbols("nu_H nu_v rho", positive=True)
    u, v, w = (sp.Function(n)(x, y, z) for n in "uvw")
    txz = rho * nv * u.diff(z) + rho * nH * w.diff(x)
    tyz = rho * nv * v.diff(z) + rho * nH * w.diff(y)
    txy = rho * nH * (u.diff(y) + v.diff(x))
    txx, tyy, tzz = 2 * rho * nH * u.diff(x), 2 * rho * nH * v.diff(y), 2 * rho * nv * w.diff(z)   # Eq. (13.5)
    F = [(txx.diff(x) + txy.diff(y) + txz.diff(z)) / rho, (txy.diff(x) + tyy.diff(y) + tyz.diff(z)) / rho,
         (txz.diff(x) + tyz.diff(y) + tzz.diff(z)) / rho]
    div = u.diff(x) + v.diff(y) + w.diff(z)
    eq6 = [nH * (q.diff(x, 2) + q.diff(y, 2)) + nv * q.diff(z, 2) for q in (u, v, w)]              # Eq. (13.6)
    extra = [nH * div.diff(x), nH * div.diff(y), nv * div.diff(z)]
    diffs = [sp.simplify(sp.expand(a - b - c)) for a, b, c in zip(F, eq6, extra)]
    res = sp.simplify(sum(d ** 2 for d in diffs))
    return dict(residual=res, ok=bool(res == 0), extra_terms=tuple(extra))


def thin_layer_sympy() -> dict:
    """A geostrophic, hydrostatic state built from any pressure field satisfies the steady inviscid thin-shell
    equations, and its shear is the thermal wind.

    Book: §13.4, Eq. (13.9) (accelerations and friction dropped); §13.5, Eqs. (13.11)–(13.15).
    Returns dict: ``residual`` (sum of squares of the three momentum residuals and the two thermal-wind differences;
    0), ``ok``, ``thermal_wind`` (the two differences).
    Assumptions: constant f, Ro ≪ 1, E ≪ 1.   Validation: V2.  Label: symbolic.
    """
    import sympy as sp

    x, y, z = sp.symbols("x y z", real=True)
    f, g, rho0 = sp.symbols("f g rho_0", nonzero=True)
    P = sp.Function("P")(x, y, z)
    rho = -P.diff(z) / g                                          # Eq. (13.14)
    u, v = -P.diff(y) / (rho0 * f), P.diff(x) / (rho0 * f)        # Eqs. (13.11), (13.12)
    rx = -f * v + P.diff(x) / rho0                                # Eq. (13.9a)
    ry = f * u + P.diff(y) / rho0                                 # Eq. (13.9b)
    rz = P.diff(z) / rho0 + g * rho / rho0                        # Eq. (13.9c)
    tw = (sp.simplify(v.diff(z) + g / (rho0 * f) * rho.diff(x)), sp.simplify(u.diff(z) - g / (rho0 * f) * rho.diff(y)))
    res = sp.simplify(rx ** 2 + ry ** 2 + rz ** 2 + tw[0] ** 2 + tw[1] ** 2)
    return dict(residual=res, ok=bool(res == 0), thermal_wind=tw)


def taylor_proudman_sympy(printed: bool = False) -> dict:
    """The Taylor–Proudman theorem from the geostrophic balance of a homogeneous fluid.

    Book: §13.5, Eqs. (13.16)–(13.21).  The corrected pair is −2Ωv = −(1/ρ)∂p/∂x and 2Ωu = −(1/ρ)∂p/∂y; the book
    prints the second with a minus on the left (slip #3).
    Parameters: printed — use the printed (13.17).
    Returns dict: ``residual`` — the horizontal divergence ∂u/∂x + ∂v/∂y of the velocity given by the two balances
    (0 for the corrected pair, so continuity gives ∂w/∂z = 0, Eq. (13.19); p_xy/(Ωρ) for the printed pair); ``ok``;
    ``printed_identity`` — ∂u/∂x − ∂v/∂y (0 for the printed pair: the relation it yields instead); ``du_dz``,
    ``dv_dz`` (both 0 with hydrostatic balance and uniform density, Eq. (13.20)).
    Assumptions: Ro ≪ 1, E ≪ 1, steady, uniform density, rotation about the vertical.
    Validation: V2 this is the check; ``printed=True`` must fail.  Label: symbolic.
    """
    import sympy as sp

    x, y, z = sp.symbols("x y z", real=True)
    Om, rho, g = sp.symbols("Omega rho g", positive=True)
    q = sp.Function("q")(x, y)
    p = q - rho * g * z                                   # hydrostatic, uniform density: dp/dz = -rho g
    s = -1 if printed else 1
    v = p.diff(x) / (2 * Om * rho)                        # Eq. (13.16)
    u = -s * p.diff(y) / (2 * Om * rho)                   # Eq. (13.17): s*2*Omega*u = -(1/rho) dp/dy
    res = sp.simplify(u.diff(x) + v.diff(y))
    return dict(residual=res, ok=bool(res == 0), printed_identity=sp.simplify(u.diff(x) - v.diff(y)),
                du_dz=sp.simplify(u.diff(z)), dv_dz=sp.simplify(v.diff(z)))


def hydrostatic_linear_set_sympy() -> dict:
    """A baroclinic Poincaré mode of a uniformly stratified layer satisfies the five hydrostatic linear equations.

    Book: §13.9, Eqs. (13.47)–(13.51).  Trial solution (l = 0, rigid-lid mode ψ = cos mz, c = N/m):
    p/ρ₀ = p̂ cos(mz) e^{i(kx − ωt)}, with ω² = f² + c²k² (Eq. (13.82) with gH → c²).
    Returns dict: ``residual`` (sum of squared moduli of the five residuals after substituting ω; 0), ``ok``,
    ``residuals`` (the five expressions).
    Assumptions: linear, hydrostatic, uniform N.   Validation: V2.  Label: symbolic.
    """
    import sympy as sp

    x, z, t = sp.symbols("x z t", real=True)
    k, m, N, f, g, rho0, ph = sp.symbols("k m N f g rho_0 phat", positive=True)
    om = sp.sqrt(f ** 2 + N ** 2 * k ** 2 / m ** 2)
    E = sp.exp(sp.I * (k * x - om * t))
    p = rho0 * ph * sp.cos(m * z) * E
    uh = k * ph * om / (om ** 2 - f ** 2)
    u = uh * sp.cos(m * z) * E
    v = -sp.I * f / om * u
    rho = -p.diff(z) / g                                  # Eq. (13.50)
    w = g / (rho0 * N ** 2) * rho.diff(t)                 # Eq. (13.51)
    res = [u.diff(x) + w.diff(z),                         # Eq. (13.47), d/dy = 0
           u.diff(t) - f * v + p.diff(x) / rho0,          # Eq. (13.48)
           v.diff(t) + f * u,                             # Eq. (13.49)
           -p.diff(z) - g * rho,                          # Eq. (13.50)
           rho.diff(t) - rho0 * N ** 2 / g * w]           # Eq. (13.51)
    res = [sp.simplify(r / E) for r in res]
    tot = sp.simplify(sum(sp.Abs(r) ** 2 for r in res))
    return dict(residual=tot, ok=bool(tot == 0), residuals=tuple(res))


def v_equation_sympy() -> dict:
    """One equation for v from the linear shallow-water set on a β-plane — exactly, with f = f₀ + βy.

    Book: §13.10, Eqs. (13.72)–(13.75).  With E1, E2, E3 the residuals of the x-momentum, y-momentum and continuity
    members of Eq. (13.45), the combination ∂_t(∂_tE2 − g∂_yE3) − f(∂_tE1 − g∂_xE3) + gH ∂_x(∂_yE1 − ∂_xE2) is
    identically ∂³v/∂t³ − gH ∂_t∇²v + f² ∂v/∂t − gHβ ∂v/∂x: Eq. (13.75) with f² in place of f₀².  The only
    approximation in (13.75) is therefore f² → f₀².
    Returns dict: ``residual`` (combination − exact v-equation; 0), ``ok``, ``exact`` (the exact equation's
    left-hand side), ``neglected_relative`` (f²/f₀² − 1 = 2βy/f₀ + (βy/f₀)², the relative size of what (13.75) drops
    in its third term).
    Assumptions: linear, hydrostatic, β-plane.   Validation: V2.  Label: symbolic.
    """
    import sympy as sp

    x, y, t = sp.symbols("x y t", real=True)
    g, H, f0, beta = sp.symbols("g H f_0 beta", positive=True)
    f = f0 + beta * y
    u, v, eta = (sp.Function(n)(x, y, t) for n in ("u", "v", "eta"))
    E1 = u.diff(t) - f * v + g * eta.diff(x)              # Eq. (13.45b)
    E2 = v.diff(t) + f * u + g * eta.diff(y)              # Eq. (13.45c)
    E3 = eta.diff(t) + H * (u.diff(x) + v.diff(y))        # Eq. (13.45a)
    comb = (E2.diff(t) - g * E3.diff(y)).diff(t) - f * (E1.diff(t) - g * E3.diff(x)) + g * H * (E1.diff(y) - E2.diff(x)).diff(x)
    exact = v.diff(t, 3) - g * H * (v.diff(x, 2) + v.diff(y, 2)).diff(t) + f ** 2 * v.diff(t) - g * H * beta * v.diff(x)
    res = sp.simplify(sp.expand(comb - exact))
    return dict(residual=res, ok=bool(res == 0), exact=exact, neglected_relative=sp.expand(f ** 2 / f0 ** 2 - 1))


def rotating_internal_wave_set_sympy() -> dict:
    """Eliminate u, v, p, ρ from the rotating, non-hydrostatic linear set: the operator acting on w.

    Book: §13.14, Eqs. (13.95) and (13.96).  Writing ∂_t, ∂_x, ∂_y, ∂_z as commuting symbols T, X, Y, Z (uniform N),
    the determinant of the 5 × 5 operator matrix of (13.95) is (1/ρ₀)·T·[T²(X² + Y² + Z²) + N²(X² + Y²) + f²Z²]:
    one time derivative times the operator of Eq. (13.96).
    Returns dict: ``residual`` (determinant − coefficient·T·operator, expanded; 0), ``ok``, ``operator`` (the
    polynomial of (13.96)), ``coefficient``.
    Assumptions: linear, Boussinesq, f-plane, uniform N.   Validation: V2; V7 f = 0 gives ch07's operator.
    Label: symbolic.
    """
    import sympy as sp

    T, X, Y, Z = sp.symbols("T X Y Z")
    N, f, g, rho0 = sp.symbols("N f g rho_0", positive=True)
    # unknowns (u, v, w, p, rho); rows: continuity, x, y, z momentum, density — Eq. (13.95)
    M = sp.Matrix([[X, Y, Z, 0, 0],
                   [T, -f, 0, X / rho0, 0],
                   [f, T, 0, Y / rho0, 0],
                   [0, 0, T, Z / rho0, g / rho0],
                   [0, 0, -rho0 * N ** 2 / g, 0, T]])
    det = sp.expand(M.det())
    op = T ** 2 * (X ** 2 + Y ** 2 + Z ** 2) + N ** 2 * (X ** 2 + Y ** 2) + f ** 2 * Z ** 2     # Eq. (13.96)
    coeff = sp.simplify(det / (T * op))
    res = sp.expand(det - coeff * T * op)
    ok = bool(res == 0) and not (coeff.free_symbols & {T, X, Y, Z})
    return dict(residual=res, ok=ok, operator=op, coefficient=coeff)


def w_equation_rotating_sympy() -> dict:
    """A plane wave in the w-equation gives the inertia–gravity dispersion relation.

    Book: §13.14, Eq. (13.96) with w ∝ e^{i(kx + ly + mz − ωt)} (uniform N) ⇒ ω²(k² + l² + m²) = N²(k² + l²) + f²m²,
    i.e. Eq. (13.112) for l = 0; with f = 0 it is ch07's ω = N cos θ.
    Returns dict: ``residual`` (Eq. (13.96) applied to the plane wave with ω from the relation, divided by the wave;
    0), ``ok``, ``omega_sq`` (the relation), ``f0_limit`` (ω² at f = 0).
    Assumptions: uniform N, f-plane.   Validation: V2.  Label: symbolic.
    """
    import sympy as sp

    x, y, z, t = sp.symbols("x y z t", real=True)
    k, l, m, N, f = sp.symbols("k l m N f", positive=True)
    om2 = (N ** 2 * (k ** 2 + l ** 2) + f ** 2 * m ** 2) / (k ** 2 + l ** 2 + m ** 2)
    om = sp.sqrt(om2)
    w = sp.exp(sp.I * (k * x + l * y + m * z - om * t))
    lap = w.diff(x, 2) + w.diff(y, 2) + w.diff(z, 2)
    res = sp.simplify((lap.diff(t, 2) + N ** 2 * (w.diff(x, 2) + w.diff(y, 2)) + f ** 2 * w.diff(z, 2)) / w)
    return dict(residual=res, ok=bool(res == 0), omega_sq=om2, f0_limit=sp.simplify(om2.subs(f, 0)))


def pv_conservation_sympy() -> dict:
    """Potential vorticity is conserved following the motion in shallow water — exactly, for any f(y) and any bottom.

    Book: §13.13, Eqs. (13.88)–(13.94).  The book reaches D/Dt[(ζ + f)/h] = 0 through the β-plane replacement
    f → f₀ in the stretching term and then restores f; here u_t, v_t, h_t are taken from (13.88)–(13.90) with an
    arbitrary f(y) and η = h + b(x, y), and D/Dt[(v_x − u_y + f)/h] is expanded.
    Returns dict: ``residual`` (h²·Dq/Dt after substitution; 0), ``ok``.
    Assumptions: one homogeneous hydrostatic layer, inviscid.   Validation: V2.  Label: symbolic.
    """
    import sympy as sp

    x, y, t = sp.symbols("x y t", real=True)
    g = sp.symbols("g", positive=True)
    u, v, h = (sp.Function(n)(x, y, t) for n in ("u", "v", "h"))
    f = sp.Function("f")(y)
    b = sp.Function("b")(x, y)
    eta = h + b
    ut = -u * u.diff(x) - v * u.diff(y) + f * v - g * eta.diff(x)      # Eq. (13.88)
    vt = -u * v.diff(x) - v * v.diff(y) - f * u - g * eta.diff(y)      # Eq. (13.89)
    ht = -(u * h).diff(x) - (v * h).diff(y)                            # Eq. (13.90)
    zeta = v.diff(x) - u.diff(y)
    q = (zeta + f) / h                                                 # Eq. (13.94)
    qt = (vt.diff(x) - ut.diff(y)) / h - (zeta + f) * ht / h ** 2
    res = sp.simplify(sp.expand((qt + u * q.diff(x) + v * q.diff(y)) * h ** 2))
    return dict(residual=res, ok=bool(res == 0))


def qg_vorticity_sympy() -> dict:
    """A plane wave in the quasi-geostrophic vorticity equation gives the Rossby dispersion relation.

    Book: §13.15, Eq. (13.117) with η ∝ e^{i(kx + ly − ωt)} ⇒ ω = −βk/(k² + l² + f₀²/c²), Eq. (13.118).
    Returns dict: ``residual`` (Eq. (13.117) applied to the wave with that ω, divided by the wave; 0), ``ok``,
    ``omega``.
    Assumptions: quasi-geostrophic, linear, β-plane.   Validation: V2.  Label: symbolic.
    """
    import sympy as sp

    x, y, t = sp.symbols("x y t", real=True)
    k, l = sp.symbols("k l", real=True)
    beta, f0, c = sp.symbols("beta f_0 c", positive=True)
    om = -beta * k / (k ** 2 + l ** 2 + f0 ** 2 / c ** 2)              # Eq. (13.118)
    eta = sp.exp(sp.I * (k * x + l * y - om * t))
    res = sp.simplify(((eta.diff(x, 2) + eta.diff(y, 2) - f0 ** 2 / c ** 2 * eta).diff(t) + beta * eta.diff(x)) / eta)
    return dict(residual=res, ok=bool(res == 0), omega=om)


def eady_qg_sympy() -> dict:
    """Two steps of the Eady problem: the interior equation, and the phase speed from the 2 × 2 determinant.

    Book: §13.17.  (a) Substituting ζ′ = ∇_H²p′/(ρ₀f) (13.132) and w′ of (13.135) into the perturbation vorticity
    equation (13.130) with U = U₀z/H gives (∂_t + U∂_x)[∇_H²p′ + (f²/N²)∂²p′/∂z²] = 0, Eq. (13.136).
    (b) The determinant of the two lid conditions (``eady_matrix``) is a quadratic in c whose roots are Eq. (13.141).
    Returns dict: ``residual`` (the interior residual of (a); 0), ``ok`` (both steps hold), ``determinant_residual``
    (the determinant divided by its leading coefficient, minus (c − U₀/2)² + (U₀/αH)²·(product of the two brackets
    … ) rewritten in exponentials and simplified; 0), ``determinant_check`` (the same difference evaluated with 30
    digits at αH = 8/5, U₀ = 1, c = 3/10 + i/7 — a guard in case the symbolic simplification is incomplete).
    Assumptions: Eady problem (quasi-geostrophic, uniform N and shear, f-plane).   Validation: V2.  Label: symbolic.
    """
    import sympy as sp

    x, y, z, t = sp.symbols("x y z t", real=True)
    f, N, rho0, U0, H = sp.symbols("f N rho_0 U_0 H", positive=True)
    p = sp.Function("p")(x, y, z, t)
    U = U0 * z / H
    adv = lambda q: q.diff(t) + U * q.diff(x)  # noqa: E731
    zeta = (p.diff(x, 2) + p.diff(y, 2)) / (rho0 * f)                                 # Eq. (13.132)
    w = -(adv(p.diff(z)) - U.diff(z) * p.diff(x)) / (rho0 * N ** 2)                   # Eq. (13.135)
    lhs = adv(zeta) - f * w.diff(z)                                                   # Eq. (13.130)
    target = adv(p.diff(x, 2) + p.diff(y, 2) + f ** 2 / N ** 2 * p.diff(z, 2)) / (rho0 * f)   # Eq. (13.136)
    res_a = sp.simplify(sp.expand(lhs - target))
    a, c = sp.symbols("a c")                    # a = alpha*H (dimensionless); lengths scaled by H
    q = a / 2
    M = sp.Matrix([[a * c * sp.sinh(q) - U0 * sp.cosh(q), -a * c * sp.cosh(q) + U0 * sp.sinh(q)],
                   [a * (U0 - c) * sp.sinh(q) - U0 * sp.cosh(q), a * (U0 - c) * sp.cosh(q) - U0 * sp.sinh(q)]])
    det = sp.expand(M.det())
    lead = det.coeff(c, 2)
    target_b = (c - U0 / 2) ** 2 - (U0 / a) ** 2 * (q - sp.tanh(q)) * (q - sp.coth(q))   # Eq. (13.141) squared
    diff_b = det / lead - target_b
    res_b = sp.simplify(sp.expand(diff_b.rewrite(sp.exp)))
    chk = abs(sp.N(diff_b.subs({a: sp.Rational(8, 5), U0: 1, c: sp.Rational(3, 10) + sp.I / 7}), 30))
    ok_b = bool(res_b == 0) or bool(chk < sp.Float("1e-25"))
    return dict(residual=res_a, ok=bool(res_a == 0) and ok_b, determinant_residual=res_b, determinant_check=chk)


def sympy_engines() -> dict:
    """The twelve symbolic engines of the chapter by name (for the script and the notebook's "the engine agrees" cells).

    Book: §13.3–§13.17 (see each function).   Returns dict name → callable.
    Assumptions: none.   Validation: V2 each returns ``ok`` True (and False with ``printed=True`` where offered).
    """
    return dict(boussinesq_rotating=boussinesq_rotating_sympy, perturbation_form=perturbation_form_sympy,
                eddy_friction_force=eddy_friction_force_sympy, thin_layer=thin_layer_sympy,
                taylor_proudman=taylor_proudman_sympy, hydrostatic_linear_set=hydrostatic_linear_set_sympy,
                v_equation=v_equation_sympy, rotating_internal_wave_set=rotating_internal_wave_set_sympy,
                w_equation_rotating=w_equation_rotating_sympy, pv_conservation=pv_conservation_sympy,
                qg_vorticity=qg_vorticity_sympy, eady_qg=eady_qg_sympy)


# =====================================================================================================================
# 4. Figure helpers (our analogues of the book's figures; every keyword defaults to illustrative_inputs())
# =====================================================================================================================

def _colors():
    from .core.style import COLORS

    return COLORS


def fig_ekman_surface(tau=None, nu_v=None, lat_rad=None, rho=None, n: int = 400):
    """Hodograph and profiles of the surface Ekman spiral (our analogue of the book's two-panel figure in §13.6).

    Book: §13.6, the solution after Eq. (13.29).   Parameters: tau [N/m²], nu_v [m²/s], lat_rad [rad], rho [kg/m³]
    (defaults: ``illustrative_inputs()`` at the Ekman latitude); n points.   Returns a matplotlib Figure.
    Assumptions: as ``ekman_surface``.   Validation: qualitative (figure).
    """
    import matplotlib.pyplot as plt

    I, C = illustrative_inputs(), _colors()
    tau = I["tau"] if tau is None else tau
    nu_v = I["nu_v_ocean"] if nu_v is None else nu_v
    lat = I["lat_ekman"] if lat_rad is None else lat_rad
    rho = I["rho_ocean"] if rho is None else rho
    f = float(coriolis_parameter(lat))
    d = float(ekman_depth(nu_v, f))
    z = np.linspace(-1.5 * np.pi * d, 0.0, n)
    u, v = ekman_surface(z, tau, 0.0, rho, nu_v, f)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    ax[0].plot(u, v, color=C["accent"])
    for s in (0.0, np.pi / 4, np.pi / 2, np.pi):
        us, vs = ekman_surface(-s * d, tau, 0.0, rho, nu_v, f)
        ax[0].annotate("", (us, vs), (0, 0), arrowprops=dict(arrowstyle="->", color=C["teal"], lw=1.2))
        ax[0].text(us, vs, f" {s / np.pi:.2g}π", fontsize=8, color=C["muted"])
    ax[0].annotate("", (0.6 * np.max(np.abs(u)), 0.0), (0, 0), arrowprops=dict(arrowstyle="->", color=C["orange"], lw=2))
    ax[0].text(0.6 * np.max(np.abs(u)), 0.0, " wind stress", color=C["orange"], va="bottom", fontsize=9)
    ax[0].set(xlabel="u [m/s]", ylabel="v [m/s]", title="hodograph (labels: −z/δ)", aspect="equal")
    ax[1].plot(u, -z / d, color=C["accent"], label="u")
    ax[1].plot(v, -z / d, color=C["teal"], label="v")
    ax[1].invert_yaxis()
    ax[1].set(xlabel="velocity [m/s]", ylabel="depth −z/δ", title=f"profiles, δ = {d:.1f} m")
    ax[1].legend()
    return fig


def fig_ekman_bottom(U_g=None, nu_v=None, lat_rad=None, n: int = 400):
    """Hodograph and profiles of the Ekman layer on a rigid surface (our analogue of the figure in §13.7).

    Book: §13.7, Eq. (13.41).   Parameters: U_g [m/s], nu_v [m²/s], lat_rad [rad] (defaults from
    ``illustrative_inputs()``); n points.   Returns a matplotlib Figure.
    Assumptions: as ``ekman_bottom``.   Validation: qualitative (figure).
    """
    import matplotlib.pyplot as plt

    I, C = illustrative_inputs(), _colors()
    U = I["U_g"] if U_g is None else U_g
    nu_v = I["nu_v_atm"] if nu_v is None else nu_v
    lat = I["lat_ekman"] if lat_rad is None else lat_rad
    f = float(coriolis_parameter(lat))
    d = float(ekman_depth(nu_v, f))
    z = np.linspace(0.0, 2.2 * np.pi * d, n)
    u, v = GFD.ekman_bottom(z, U, 0.0, nu_v, f)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    ax[0].plot(u, v, color=C["accent"])
    for s in (np.pi / 4, np.pi / 2, np.pi):
        us, vs = GFD.ekman_bottom(s * d, U, 0.0, nu_v, f)
        ax[0].plot([us], [vs], "o", color=C["teal"], ms=4)
        ax[0].text(us, vs, f" {s / np.pi:.2g}π", fontsize=8, color=C["muted"])
    ax[0].axvline(U, color=C["muted"], lw=0.8, ls="--")
    ax[0].set(xlabel="u [m/s]", ylabel="v [m/s] (toward low pressure)", title="hodograph (labels: z/δ)", aspect="equal")
    ax[1].plot(u, z / d, color=C["accent"], label="u")
    ax[1].plot(v, z / d, color=C["teal"], label="v")
    ax[1].axvline(U, color=C["muted"], lw=0.8, ls="--")
    ax[1].set(xlabel="velocity [m/s]", ylabel="height z/δ", title=f"profiles, δ = {d:.0f} m")
    ax[1].legend()
    return fig


def fig_mode_roots(N=None, H=None, g: float = G0, n_roots: int = 3, exaggerate: float = 200.0):
    """The two curves whose crossings give the mode speeds for uniform N: tan X against ε/X, X = NH/c.

    Book: §13.9, Eq. (13.69) (our analogue of the book's sketch).  ε = N²H/g is so small for an ocean that the
    crossings would be invisible at true scale: the ε/X curve is drawn multiplied by ``exaggerate`` (stated in the
    legend) and the true roots are marked.
    Parameters: N [rad/s], H [m] (defaults from ``illustrative_inputs()``), g, n_roots, exaggerate.
    Returns a matplotlib Figure.   Assumptions: uniform N.   Validation: qualitative (figure).
    """
    import matplotlib.pyplot as plt

    I, C = illustrative_inputs(), _colors()
    N = I["ocean_N"] if N is None else N
    H = I["ocean_H"] if H is None else H
    eps = N * N * H / g
    fig, ax = plt.subplots(figsize=(8, 4.2))
    for n in range(n_roots):
        X = np.linspace(n * np.pi - np.pi / 2 + 0.05, n * np.pi + np.pi / 2 - 0.05, 300)
        X = X[X > 0]
        ax.plot(X, np.tan(X), color=C["accent"], label="tan X" if n == 0 else None)
    Xc = np.linspace(0.02, (n_roots - 0.5) * np.pi, 400)
    ax.plot(Xc, exaggerate * eps / Xc, color=C["orange"], label=f"{exaggerate:g} × (N²H/g)/X")
    modes = VM.modes_uniform_N(N, H, g, n_roots)
    Xr = N * H / modes.c
    ax.plot(Xr, np.zeros_like(Xr), "o", color=C["teal"], label="roots X_n = NH/c_n (true scale)")
    ax.set(ylim=(-4, 4), xlabel="X = N H / c", ylabel="", title=f"mode speeds for uniform N: X₀ = {Xr[0]:.4f}, X₁ = {Xr[1]:.4f}")
    ax.legend(fontsize=8)
    return fig


def fig_vertical_modes(z=None, N2=None, n_modes: int = 4, lid: str = "free", g: float = G0):
    """Mode shapes ψ_n(z) of a stratification beside its N(z) profile.

    Book: §13.9, Eq. (13.56) (our analogue of the book's sketch of the first modes).
    Parameters: z nodes [m] and N2 at them [1/s²] (defaults: 401 nodes over the illustrative ocean depth with
    ``thermocline_N2``); n_modes; lid; g.   Returns a matplotlib Figure.
    Assumptions: as ``vertical_modes``.   Validation: qualitative (figure).
    """
    import matplotlib.pyplot as plt

    I, C = illustrative_inputs(), _colors()
    if z is None:
        z = np.linspace(-I["ocean_H"], 0.0, 401)
    if N2 is None:
        N2 = thermocline_N2(z)
    modes = VM.vertical_modes(z, N2, g, n_modes, lid)
    fig, ax = plt.subplots(1, 2, figsize=(9, 4.4), sharey=True, gridspec_kw=dict(width_ratios=[1, 2]))
    ax[0].plot(np.sqrt(_F(N2) * np.ones_like(z)) * 1e3, z, color=C["blue"])
    ax[0].set(xlabel="N [10⁻³ rad/s]", ylabel="z [m]", title="stratification")
    first = 0 if lid == "free" else 1
    for j in range(modes.psi.shape[0]):
        ax[1].plot(modes.psi[j], z, label=f"n = {j + first}: c = {modes.c[j]:.3g} m/s")
    ax[1].axvline(0, color=C["muted"], lw=0.8)
    ax[1].set(xlabel="ψ_n (ψ = 1 at the surface)", title=f"vertical modes ({lid} surface)")
    ax[1].legend(fontsize=8)
    return fig


def fig_poincare_kelvin_dispersion(H=None, lat_rad=None, g: float = G0, with_rossby: bool = True):
    """Frequency against wavenumber for Poincaré, Kelvin and (optionally) Rossby waves on logarithmic axes.

    Book: §13.11–§13.12 (our analogue of the Poincaré–Kelvin sketch; the logarithmic axes and the third, slow branch
    are our additions).   Parameters: H [m], lat_rad [rad] (defaults from ``illustrative_inputs()``); g; with_rossby.
    Returns a matplotlib Figure.   Assumptions: linear shallow water.   Validation: qualitative (figure).
    """
    import matplotlib.pyplot as plt

    I, C = illustrative_inputs(), _colors()
    H = I["ocean_H"] if H is None else H
    lat = I["lat"] if lat_rad is None else lat_rad
    f, beta = float(coriolis_parameter(lat)), float(beta_parameter(lat))
    c = float(np.sqrt(g * H))
    Lam = c / abs(f)
    K = np.geomspace(0.05, 50.0, 300) / Lam
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    ax.loglog(K * Lam, GFD.poincare_omega(K, f, c) / abs(f), color=C["accent"], label="Poincaré: ω² = f² + c²K²")
    ax.loglog(K * Lam, c * K / abs(f), color=C["teal"], ls="--", label="Kelvin: ω = cK")
    if with_rossby:
        ax.loglog(K * Lam, np.abs(GFD.rossby_omega(-K, 0.0, beta, f, c)) / abs(f), color=C["amber"],
                  label="Rossby: ω = −βk/(k² + f²/c²)")
    ax.axhline(1.0, color=C["muted"], lw=0.8)
    ax.set(xlabel="K Λ  (Λ = c/|f|)", ylabel="ω / |f|", title=f"shallow-water waves, Λ = {Lam / 1e3:.0f} km")
    ax.legend(fontsize=8)
    return fig


def fig_kelvin_sections(H=None, lat_rad=None, eta0: float = 0.5, wavelength: float | None = None, g: float = G0):
    """A Kelvin wave seen two ways: the surface across the shelf under a crest and a trough, and in plan.

    Book: §13.12, Eq. (13.87) (our analogue of the book's transverse sections).   Parameters: H [m], lat_rad [rad]
    (defaults from ``illustrative_inputs()``), eta0 [m], wavelength [m] (default 6 Λ), g.
    Returns a matplotlib Figure.   Assumptions: as ``kelvin_wave``.   Validation: qualitative (figure).
    """
    import matplotlib.pyplot as plt

    I, C = illustrative_inputs(), _colors()
    H = I["ocean_H"] if H is None else H
    lat = I["lat"] if lat_rad is None else lat_rad
    f = float(coriolis_parameter(lat))
    Lam = float(np.sqrt(g * H)) / abs(f)
    lam = 6.0 * Lam if wavelength is None else wavelength
    k = 2.0 * np.pi / lam
    y = np.linspace(0.0, 4.0 * Lam, 200)
    x = np.linspace(0.0, 2.0 * lam, 240)
    d = 1 if f > 0 else -1
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    crest, _ = GFD.kelvin_wave(0.0, y, 0.0, eta0, k, H, f, g)
    ax[0].plot(y / Lam, crest, color=C["accent"], label="under a crest")
    ax[0].plot(y / Lam, -crest, color=C["rose"], label="under a trough")
    ax[0].plot(y / Lam, eta0 * np.exp(-y / Lam), color=C["muted"], ls=":", label="η₀ e^{−y/Λ}")
    ax[0].set(xlabel="distance from the coast y/Λ", ylabel="η [m]", title=f"cross-shore sections, Λ = {Lam / 1e3:.0f} km")
    ax[0].legend(fontsize=8)
    X, Y = np.meshgrid(x, y)
    eta, u = GFD.kelvin_wave(X, Y, 0.0, eta0, k, H, f, g)
    im = ax[1].pcolormesh(X / lam, Y / Lam, eta, cmap="RdBu_r", shading="auto", vmin=-eta0, vmax=eta0)
    sl = (slice(None, None, 25), slice(None, None, 20))
    ax[1].quiver(X[sl] / lam, Y[sl] / Lam, u[sl], 0 * u[sl], color=C["ink"], width=0.004)
    ax[1].set(xlabel="x / wavelength", ylabel="y/Λ", title="plan view: travels toward " + ("+x" if d > 0 else "−x"))
    fig.colorbar(im, ax=ax[1], label="η [m]")
    return fig


def fig_inertia_gravity_orbit(omega_over_f: float = 2.5, lat_rad=None):
    """Horizontal velocity hodograph of an inertia–gravity wave: an ellipse with axis ratio |f|/ω.

    Book: §13.14, Eq. (13.110) (our analogue of the book's hodograph; the same ellipse as the Poincaré wave of
    §13.11).   Parameters: omega_over_f (> 1); lat_rad (its sign decides the sense; default ``illustrative_inputs()``).
    Returns a matplotlib Figure.   Assumptions: as ``inertia_gravity_hodograph``.   Validation: qualitative (figure).
    """
    import matplotlib.pyplot as plt

    I, C = illustrative_inputs(), _colors()
    lat = I["lat"] if lat_rad is None else lat_rad
    fig, ax = plt.subplots(1, 2, figsize=(9, 4.0))
    for a, la in zip(ax, (abs(lat), -abs(lat))):
        f = float(coriolis_parameter(la))
        om = omega_over_f * abs(f)
        t = np.linspace(0.0, 2.0 * np.pi / om, 200)
        u, v = GFD.inertia_gravity_hodograph(t, om, f)
        a.plot(u, v, color=C["accent"])
        for frac in (0.0, 0.25, 0.5):
            ui, vi = GFD.inertia_gravity_hodograph(frac * 2.0 * np.pi / om, om, f)
            a.plot([ui], [vi], "o", color=C["teal"], ms=5)
            a.text(ui, vi, f" ωt = {2 * frac:g}π", fontsize=8)
        a.set(xlabel="u / û", ylabel="v / û", aspect="equal",
              title=("f > 0: clockwise" if f > 0 else "f < 0: counter-clockwise") + f", axis ratio {1 / omega_over_f:.2f}")
    return fig


def fig_rossby_dispersion(lat_rad=None, c=None):
    """Rossby-wave frequency against zonal wavenumber, and the circles of constant frequency in the wavenumber plane.

    Book: §13.15, Eq. (13.118) and the circle form after it (our analogue of the book's two-panel figure, with our
    own frequencies).   Parameters: lat_rad [rad], c [m/s] (defaults: the illustrative latitude and the first
    baroclinic speed of the uniform-N ocean).   Returns a matplotlib Figure.
    Assumptions: as ``rossby_omega``.   Validation: qualitative (figure).
    """
    import matplotlib.pyplot as plt

    I, C = illustrative_inputs(), _colors()
    lat = I["lat_rossby"] if lat_rad is None else lat_rad
    c = float(GFD.baroclinic_mode_speed(I["ocean_N"], I["ocean_H"], 1)) if c is None else c
    f, beta = float(coriolis_parameter(lat)), float(beta_parameter(lat))
    Lam = c / abs(f)
    wmax = GFD.rossby_max_frequency(beta, f, c)["omega_max"]
    k = -np.linspace(0.0, 4.0, 300) / Lam
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    ax[0].plot(k * Lam, GFD.rossby_omega(k, 0.0, beta, f, c) / wmax, color=C["accent"])
    ax[0].plot([-1.0], [1.0], "o", color=C["orange"])
    ax[0].axvline(-1.0, color=C["muted"], ls="--", lw=0.8)
    ax[0].text(-0.95, 0.2, "c_gx < 0\n(long waves)", fontsize=8)
    ax[0].text(-3.2, 0.2, "c_gx > 0 (short waves)", fontsize=8)
    ax[0].set(xlabel="k Λ", ylabel="ω / ω_max", title=f"l = 0; ω_max = βc/2|f| → period {2 * np.pi / wmax / 86400:.0f} days")
    th = np.linspace(0, 2 * np.pi, 200)
    for frac, col in zip((0.45, 0.7, 0.9), (C["accent"], C["teal"], C["orange"])):
        cir = GFD.rossby_omega_circle(frac * wmax, beta, f, c)
        ax[1].plot((cir["center_k"] + cir["radius"] * np.cos(th)) * Lam, cir["radius"] * np.sin(th) * Lam, color=col,
                   label=f"ω = {frac:g} ω_max")
        kk, ll = cir["center_k"] + cir["radius"] * np.cos(2.2), cir["radius"] * np.sin(2.2)
        cg = np.array(GFD.rossby_group_velocity(kk, ll, beta, f, c))
        cg = cg / np.hypot(*cg) * 0.5 / Lam
        ax[1].annotate("", ((kk + cg[0]) * Lam, (ll + cg[1]) * Lam), (kk * Lam, ll * Lam),
                       arrowprops=dict(arrowstyle="->", color=col))
    ax[1].set(xlabel="k Λ", ylabel="l Λ", aspect="equal", title="circles of constant ω; arrows: group velocity")
    ax[1].legend(fontsize=8)
    return fig


# =====================================================================================================================
# 5. Cached model runs (reference/ch13) — our own output only
# =====================================================================================================================

REFERENCE_RUNS = ("kelvin_basin", "turbulence_f", "turbulence_beta", "pv_particles")


def reference_run(name: str, fast: bool = False) -> dict:
    """Compute one of the chapter's model runs (the content of ``reference/ch13/<name>.npz``).

    Book: none of these runs is in the book — they are **our numerical illustrations** of §13.12 (a Kelvin wave
    going round a basin), §13.18 (decaying two-dimensional turbulence without and with β) and §13.13 (potential
    vorticity following particles).  All are deterministic (seeded).
    Parameters: name — one of ``REFERENCE_RUNS``; fast — a smaller, shorter version (for a live fallback in the
    notebook; the cached files are made with ``fast=False``).
    Returns a dict of numpy arrays (frames as float32 here; the cache stores float16) and 0-d parameter arrays:
    * "kelvin_basin": ``x``, ``y`` [m], ``t`` [s], ``eta`` (nt, ny, nx) [m], ``Lambda``, ``c``, ``f``, ``H``,
      ``eta0`` — a first-baroclinic-like layer (equivalent depth 1.3 m) at 35° N in a closed square basin 10 Λ wide,
      started with a Kelvin pulse on the south wall;
    * "turbulence_f", "turbulence_beta": ``x`` [box units], ``t``, ``zeta`` (nt, n, n), ``energy``, ``enstrophy``,
      ``K_E``, ``K_Z`` per saved frame, ``beta``, ``nu``, ``L`` — non-dimensional (box 2π, unit rms speed), seed 3;
    * "pv_particles": ``x``, ``y``, ``t``, ``eta`` (nt, ny, nx), ``xp``, ``yp``, ``q`` (nt, n_particles), ``f``,
      ``H`` — a nonlinear geostrophic vortex in a doubly periodic box.
    Assumptions: as the models in ``core.shallow_water``.   Validation: reproduced exactly by re-running (the
    script ``scripts/ch13_make_caches.py`` asserts the stored file against a fresh run).  Label: qualitative
    (illustrations), with the conservation numbers quoted in ``core.shallow_water``.
    """
    I = illustrative_inputs()
    f = float(coriolis_parameter(I["lat"]))
    if name == "kelvin_basin":
        He, n = 1.3, (32 if fast else 64)
        c = float(np.sqrt(G0 * He))
        Lam = c / abs(f)
        Lb = 10.0 * Lam
        md = SW.ShallowWater(n, n, Lb, Lb, He, f, bc="closed")
        Xc, Yc = md.grids()["center"]
        Xu, Yu = md.grids()["u"]
        eta0, sig, x0 = 0.5, 1.5 * Lam, 0.3 * Lb
        shape = lambda X, Y: eta0 * np.exp(-Y / Lam) * np.exp(-((X - x0) ** 2) / (2 * sig ** 2))  # noqa: E731
        st = SW.make_state(md, shape(Xc, Yc), np.sqrt(G0 / He) * shape(Xu, Yu), None)
        t_end = (1.5 if fast else 4.4) * Lb / c       # 4 Lb/c is one circuit of the basin
        dt = SW.dt_limit(md, 0.6)
        nst = int(np.ceil(t_end / dt))
        r = SW.run(md, st, t_end, dt=t_end / nst, save_every=max(nst // (8 if fast else 22), 1))
        return dict(x=md.x_c, y=md.y_c, t=r["t"], eta=r["eta"].astype(np.float32), Lambda=np.array(Lam), c=np.array(c),
                    f=np.array(f), H=np.array(He), eta0=np.array(eta0))
    if name in ("turbulence_f", "turbulence_beta"):
        n, L = (32 if fast else 64), 2.0 * np.pi
        beta = 0.0 if name == "turbulence_f" else 8.0
        nu = 4e-4 if fast else 2e-4
        z0 = SW.random_vorticity(n, L, k_peak=8.0 if not fast else 5.0, u_rms=1.0, seed=3)
        dt = SW.barotropic_time_step(z0, L, beta=beta, nu=nu)
        t_end = 6.0 if fast else 16.0
        nst = int(np.ceil(t_end / dt))
        r = SW.barotropic_run(n, L, z0, t_end=t_end, dt=t_end / nst, beta=beta, nu=nu, save_every=max(nst // 8, 1))
        return dict(x=np.arange(n) * L / n, t=r["t"], zeta=r["zeta"].astype(np.float32), energy=r["energy"],
                    enstrophy=r["enstrophy"], K_E=r["K_E"], K_Z=r["K_Z"], beta=np.array(beta), nu=np.array(nu),
                    L=np.array(L))
    if name == "pv_particles":
        n, Lb, Hs = (32 if fast else 48), 4.0e6, 300.0
        md = SW.ShallowWater(n, n, Lb, Lb, Hs, f)
        st = SW.geostrophic_state(md, SW.gaussian_bump(md, 30.0, radius=4.0e5))
        rng = np.random.default_rng(0)
        npart = 24
        xp0 = Lb * (0.3 + 0.4 * rng.random(npart))
        yp0 = Lb * (0.3 + 0.4 * rng.random(npart))
        t_end = (1.0 if fast else 3.0) * 2.0 * np.pi / abs(f)
        r = SW.run(md, st, t_end, linear=False, save_every=2)
        xp, yp = SW.advect_particles(md, r, xp0, yp0)
        q = SW.particle_potential_vorticity(md, r, xp, yp)
        keep = np.unique(np.linspace(0, len(r["t"]) - 1, 13).astype(int))
        return dict(x=md.x_c, y=md.y_c, t=r["t"][keep], eta=r["eta"][keep].astype(np.float32), xp=xp[keep],
                    yp=yp[keep], q=q[keep], f=np.array(f), H=np.array(Hs))
    raise ValueError(f"unknown reference run {name!r}; choose from {REFERENCE_RUNS}")


_FRAME_KEYS = {"kelvin_basin": ("eta",), "turbulence_f": ("zeta",), "turbulence_beta": ("zeta",), "pv_particles": ("eta",)}


def write_reference_runs(which="all", dest=None) -> dict:
    """Compute the chapter's model runs and write them to ``reference/ch13/<name>.npz`` (frames as float16), together
    with ``explainer_constants.json``.

    Book: none (our own model output; see :func:`reference_run`).
    Parameters: which — "all" or an iterable of names from ``REFERENCE_RUNS``; dest folder (default
    ``reference/ch13``).
    Returns dict name → dict(path, bytes, fresh) where ``fresh`` is the just-computed run (so the caller can compare
    it with what :func:`load_reference_run` reads back).
    Assumptions: none.   Validation: V1 the round trip agrees to the float16 quantisation.  Label: analytic.
    """
    import json

    folder = Path(dest) if dest is not None else _root() / "reference" / "ch13"
    folder.mkdir(parents=True, exist_ok=True)
    names = REFERENCE_RUNS if which == "all" else tuple(which)
    out = {}
    for name in names:
        run_ = reference_run(name)
        store = {k: (v.astype(np.float16) if k in _FRAME_KEYS[name] else np.asarray(v)) for k, v in run_.items()}
        path = folder / f"{name}.npz"
        np.savez_compressed(path, **store)
        out[name] = dict(path=path, bytes=path.stat().st_size, fresh=run_)
    cpath = folder / "explainer_constants.json"
    cpath.write_text(json.dumps(explainer_constants(), indent=1), encoding="utf-8")
    out["explainer_constants"] = dict(path=cpath, bytes=cpath.stat().st_size, fresh=explainer_constants())
    return out


def load_reference_run(name: str, folder=None):
    """Load a cached model run from ``reference/ch13/<name>.npz``.

    Book: none (our own model output; see :func:`reference_run` for the keys of each run).
    Parameters: name — "kelvin_basin", "turbulence_f", "turbulence_beta" or "pv_particles"; folder (default
    ``reference/ch13``).
    Returns a dict of arrays (frames converted back to float64), or ``None`` when the file is absent — the caller
    then falls back to ``reference_run(name, fast=True)``.
    Assumptions: none.   Validation: V1 round trip (``scripts/ch13_make_caches.py``).  Label: analytic.
    """
    if name not in REFERENCE_RUNS:
        raise ValueError(f"unknown reference run {name!r}; choose from {REFERENCE_RUNS}")
    path = (Path(folder) if folder is not None else _root() / "reference" / "ch13") / f"{name}.npz"
    if not path.exists():
        return None
    with np.load(path) as data:
        return {k: (data[k].astype(float) if data[k].dtype == np.float16 else data[k]) for k in data.files}


def explainer_constants() -> dict:
    """The constants an explainer could compute itself, computed here once — to pin its parity rows, never to type.

    Book: §13.17 (the Eady numbers are our computed values of Eq. (13.141); the book prints rounded ones).
    Returns dict: ``eady_critical`` (α_cH), ``eady_fastest`` (dict of ``eady_fastest()``), ``eady_wavelengths``,
    ``adjustment_ratio`` (kinetic energy of the adjusted jet over the potential energy released), ``omega_earth``,
    ``earth_radius`` and ``inputs`` (``illustrative_inputs()`` with tuples as lists).
    Assumptions: none.   Validation: V1 equal to the functions named (the cache script asserts it).  Label: analytic.
    """
    I = illustrative_inputs()
    f = float(coriolis_parameter(I["lat"]))
    return dict(eady_critical=GFD.eady_critical(), eady_fastest=GFD.eady_fastest(), eady_wavelengths=GFD.eady_wavelengths(),
                adjustment_ratio=GFD.adjustment_energy(0.5, I["ocean_H"], f)["ratio"], omega_earth=OMEGA_EARTH,
                earth_radius=EARTH_RADIUS_MEAN, inputs={k: (list(v) if isinstance(v, tuple) else v) for k, v in I.items()})
