"""Chapter 14 — Aerodynamics (Kundu, Cohen & Dowling 5e, §14.1–§14.9, Eqs. (14.1)–(14.27)).

Lift is circulation; the sharp trailing edge chooses it (Kutta condition); the wing tips make you pay for it (trailing
vortices, downwash, induced drag) — and the same lift, normal to the *relative* stream, drives fish, birds and sails.

**This module re-exports every public name of the two new core modules**

* ``fluidpy.core.airfoils`` (``AF``) — Zhukhovsky shapes (14.3)–(14.8), the flow at incidence (14.9), the Kutta
  circulation (14.10)–(14.11), the lift coefficient (14.12), surface pressure, section forces, thin-airfoil theory;
* ``fluidpy.core.lifting_line`` (``LL``) — planforms (14.1)–(14.2), downwash (14.13)–(14.14), strip forces (14.15),
  the lifting-line equation (14.16)–(14.20), series results (14.21)–(14.22), the elliptic loading (14.23)–(14.27),
  horseshoe vortices;

and the earlier-chapter functions the chapter restates (``lift_per_span`` (6.62), ``cylinder_surface_speed`` (6.37),
``joukowski`` (6.65) = (14.3), ``segment_speed`` / ``segment_induced_velocity`` ((5.17) integrated), ``blasius_force``
(6.60), ``lift_coefficient`` / ``drag_coefficient`` (4.107)–(4.108), ``mach_number``, ``reynolds_number``), so that
``ch14.<name>`` reaches everything an explainer parity row needs.

What is added here: the chapter's tables (:func:`illustrative_inputs`, :func:`conventions_table`, :func:`book_slips`,
:func:`traps`), level flight and the drag polar (§14.2, §14.7), the leading-edge suction and trailing-edge wedge of
the exercises, the lumped-vortex impulsive start (the starting vortex, §14.3), wake roll-up and tip-vortex estimates
(§14.6), the oscillating foil (§14.8), sailing (§14.9), and a sympy engine for the circular arc of §14.4.

Conventions: see the two core modules.  **Γ is clockwise** (keyword ``Gamma_cw=``); x is the mapped chord axis in
§14.4–§14.5 (lift is not F_y), y is the span in §14.6 and w is the downwash, positive downward; angles in radians;
the lift-curve slope K is per radian and always passed explicitly; the span efficiency is ``e_span`` and the
thickness parameter ``e``.  **Qualitative models of ours** (not book formulas) say so in their label:
:func:`suction_side_separation`, :func:`wake_rollup` past the first turn, :func:`generic_lift_curve`,
:func:`generic_drag_curve`, :func:`stall_regime`; :func:`profile_drag_estimate` is approximate;
:func:`impulsive_start_lumped_vortex` is a model.  Inputs used in examples are ours (:func:`illustrative_inputs`),
never the book's.
"""
from __future__ import annotations

from typing import Callable

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

from .core import airfoils as AF
from .core import lifting_line as LL
from .core.airfoils import *  # noqa: F401,F403
from .core.lifting_line import *  # noqa: F401,F403
from .core._util import as_scalar_if_0d
from .core.biot_savart import (point_vortex_evolve, point_vortex_invariants, point_vortex_velocity,  # noqa: F401
                               segment_induced_velocity, segment_speed)
from .core.boundary_layer import plate_drag_coefficient, thwaites
from .core.conformal import (circle_flow_zeta, joukowski, joukowski_derivative, joukowski_inverse,  # noqa: F401
                             mapped_flow)
from .core.potential import blasius_force, pressure_coefficient  # noqa: F401
from .core.similarity import drag_coefficient, lift_coefficient, mach_number, reynolds_number  # noqa: F401
from .core.thermo import G0
from .ch06_ideal_flow import cylinder_stagnation_points, cylinder_surface_speed, lift_per_span  # noqa: F401

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d

_REUSED = ["lift_per_span", "cylinder_surface_speed", "cylinder_stagnation_points", "joukowski",
           "joukowski_derivative", "joukowski_inverse", "circle_flow_zeta", "mapped_flow", "blasius_force",
           "pressure_coefficient", "segment_speed", "segment_induced_velocity", "point_vortex_velocity",
           "point_vortex_evolve", "point_vortex_invariants", "lift_coefficient", "drag_coefficient", "mach_number",
           "reynolds_number"]
_OWN = ["illustrative_inputs", "conventions_table", "book_slips", "traps", "compressibility_check", "level_flight",
        "stall_speed", "surface_point", "edge_speed", "map_point", "edge_flow_suction", "trailing_edge_wedge", "suction_side_separation", "chord_reynolds",
        "stall_regime", "impulsive_start_lumped_vortex", "wagner_function", "wake_rollup", "tip_vortex_estimate",
        "vortex_pair_energy_per_length", "wake_vortex_rolling_moment", "generic_lift_curve", "generic_drag_curve",
        "drag_polar", "drag_breakdown", "drag_share", "profile_drag_estimate", "best_glide", "drag_vs_speed",
        "speed_min_drag", "oscillating_foil", "apparent_wind", "sail_forces", "sail_drive_curve",
        "velocity_made_good", "sailboat_course_limit", "sailboat_equilibrium", "projectile_range_drag",
        "circular_arc_sym", "AF", "LL"]
__all__ = sorted(set(AF.__all__) | set(LL.__all__) | set(_REUSED) | set(_OWN))


# ======================================================================================================================
# tables
# ======================================================================================================================
def illustrative_inputs() -> dict:
    """The one set of inputs every worked number, script, notebook cell and explainer of the chapter uses — ours,
    visibly not the book's (it is the input set of ``analysis/ch14_design.md``, header convention 11).

    Book: none (the book's example and exercise values stay in the private test file).
    Returns a flat dict.
    **Section**: ``b`` [m] (chord ≈ 1 m), ``e`` [–], ``beta`` [rad], ``alpha`` [rad], ``U_section`` [m/s].
    **Aircraft**: ``W`` [N], ``s`` [m], ``A`` [m²], ``AR``, ``c_mean`` [m], ``U_slow``, ``U_cruise`` [m/s], ``U``
    (= ``U_cruise``), ``rho`` [kg/m³] (sea level, USSA-1976), ``nu`` [m²/s], ``CD0``, ``e_span``, ``K`` [1/rad],
    ``CL_max``, ``CL_max_flap``, ``taper``.
    **Demo wing** (rectangular, Λ = 6, non-dimensional): ``wd_s`` [m], ``wd_c`` [m], ``wd_alpha`` [rad], ``wd_U``
    [m/s], ``wd_rho`` [kg/m³], ``wd_K`` [1/rad].
    **Qualitative models**: ``alpha_stall`` [rad], ``CL_post``, ``stall_width`` [rad], ``CD_sep``, ``Re_laminar``,
    ``Re_turbulent``.
    **Boat**: ``U_wind``, ``V_boat`` [m/s], ``V_over_U``, ``theta`` [rad], ``A_sail`` [m²], ``CL_sail``, ``CD_sail``,
    ``LD_hull``, ``CDA_hull`` [m²], ``rho_water`` [kg/m³].
    **Fish**: ``U_fish``, ``V_tail`` [m/s], ``L_tail``, ``D_tail`` [N], ``LD_tail``.
    **Wake encounter**: ``Gamma_v`` [m²/s], ``ell`` [m].
    Assumptions: none.   Validation: V1 a table.  Label: analytic.
    """
    s, A = 12.0, 15.0
    return dict(b=0.25, e=0.08, beta=float(np.deg2rad(4.0)), alpha=float(np.deg2rad(5.0)), U_section=30.0,
                W=9000.0, s=s, A=A, AR=s ** 2 / A, c_mean=A / s, U_slow=30.0, U_cruise=55.0, U=55.0, rho=1.225,
                nu=1.46e-5, CD0=0.025, e_span=0.9, K=float(2.0 * np.pi), CL_max=1.4, CL_max_flap=2.0, taper=0.45,
                wd_s=6.0, wd_c=1.0, wd_alpha=float(np.deg2rad(5.0)), wd_U=1.0, wd_rho=1.0, wd_K=float(2.0 * np.pi),
                alpha_stall=float(np.deg2rad(14.0)), CL_post=0.9, stall_width=float(np.deg2rad(2.0)), CD_sep=0.15,
                Re_laminar=5.0e5, Re_turbulent=3.0e6,
                U_wind=6.0, V_boat=3.0, V_over_U=0.5, theta=float(np.deg2rad(45.0)), A_sail=20.0, CL_sail=1.2,
                CD_sail=0.2, LD_hull=8.0, CDA_hull=0.05, rho_water=1000.0,
                U_fish=1.0, V_tail=0.75, L_tail=10.0, D_tail=1.0, LD_tail=10.0, Gamma_v=310.0, ell=1.7)


def conventions_table():
    """The overloaded letters of the chapter and how the notebook, the explainers and the code tell them apart.

    Book: §14.2–§14.9 (traps T1–T5, T12, T13: Γ clockwise; b, a, β, c, e, Λ, γ, θ change meaning between sections).
    Returns a ``pandas.DataFrame`` with columns "symbol", "meanings", "ours", "code".
    Assumptions: none.   Validation: V1 a table.  Label: analytic.
    """
    import pandas as pd

    rows = [
        ("Γ", "clockwise circulation round the section (book: (6.62), (14.10), (14.11)); bound circulation Γ(y) in §14.6",
         "always labelled 'Γ (clockwise)'; lift ρUΓ > 0 for Γ > 0", "Gamma_cw= (keyword); Gamma_ccw = −Gamma_cw"),
        ("dw/dz", "u − iv (ch06 and here); ch13's Ekman velocity was u + iv", "complex velocity, conjugate of the vector",
         "ZhukhovskyFlow.dwdz, .velocity"),
        ("b", "constant of the map z = ζ + b²/ζ (chord ≃ 4b); the span in two exercises", "map constant only; span is s", "b"),
        ("a", "radius of the ζ-plane circle; a length scale in the leading-edge exercise; cylinder radius in ch06",
         "circle radius, a = b(1 + e)/cos β", "a"),
        ("β", "camber angle of the circle (§14.4–§14.5) = minus the zero-lift angle (§14.6–§14.7); df/dy in ch13",
         "camber angle / −(zero-lift angle), radians", "beta"),
        ("e", "thickness parameter, a ≡ b(1 + e) (14.7)", "thickness parameter; span efficiency is e_span", "e, e_span"),
        ("c", "chord; wave speed in ch07 and ch13", "chord [m]", "c"),
        ("Λ", "aspect ratio s²/A (14.1); the Rossby radius in ch13", "aspect ratio", "AR"),
        ("γ", "span angle, y = −(s/2) cos γ (§14.6); sheet strength in ch05; ratio of specific heats", "span angle", "gamma"),
        ("θ", "polar angle of ζ in (14.4); angle from the stream-aligned diameter in (14.10); stroke / course angle in §14.8–§14.9",
         "three different angles, each drawn where it is used", "theta (local to each function), t (contour parameter)"),
        ("x, y", "§14.3–§14.5: section plane, x the mapped chord, y up; §14.6: y along the span", "axes drawn at the top of each part",
         "(Fx, Fy) and (L, D) both returned"),
        ("w", "downwash, positive downward (§14.6); the complex potential (§14.5); upward in one exercise", "downwash (down +)",
         "w; horseshoe_velocity returns uz (up +)"),
        ("K", "section lift-curve slope, per radian in (14.16)–(14.20), per degree in the prose of §14.7", "per radian", "K= (required)"),
        ("A", "planform area in (4.107)–(4.108) as used here and in (14.26); frontal area for bluff bodies and the projectile exercise",
         "planform area unless stated", "A"),
        ("U, V, U_r", "true stream / wind U, speed of the foil or boat V, relative (apparent) velocity U_r = U − V", "three arrows",
         "apparent_wind, oscillating_foil"),
    ]
    return pd.DataFrame(rows, columns=["symbol", "meanings", "ours", "code"])


def book_slips() -> dict:
    """The printed slips and loose statements of chapter 14, each with the way to tell the two versions apart.

    Book: §14.3–§14.6 (found by reading the page images against the book's own neighbouring equations; slips #1, #2,
    #4 and #8 were re-read on the rendered pages p732, p734, p739, p741 when this module was written).
    Returns a dict keyed "1" … "10"; each value has ``id``, ``kind`` ("slip" = false as printed, "loose" = an
    inconsistency), ``where``, ``printed``, ``correct`` (alias ``corrected``), ``how_to_tell``, ``coded_in`` and
    ``test`` (the planted wrong variant a test can fail; "—" for wording-only items).
    Assumptions: none.   Validation: V1 a table; the computable rows are checked by the named variants.
    Label: analytic.
    """
    rows = [
        ("slip", "(14.19), first two members", "denominator cos γ₁ − cos γ with the prefactor +1/(2πs)",
         "cos γ − cos γ₁ (or a minus sign in front); the last member and (14.20)–(14.26) are right",
         "from (14.13) with y = −(s/2) cos γ: y₁ − y = (s/2)(cos γ − cos γ₁); the printed middle member gives −Γ₁/2s for the "
         "elliptic loading, direct quadrature of (14.13) gives +Γ₁/2s",
         "LL.downwash_series, LL.glauert_integral, LL.downwash", "LL.downwash_series(printed_middle=True)"),
        ("slip", "unnumbered line before (14.14)", "ε = tan(w/U)", "ε = tan⁻¹(w/U) ≃ w/U",
         "the velocity triangle of Fig. 14.23: tan ε = w/U", "LL.induced_angle", "np.tan(w/U) instead of np.arctan(w/U)"),
        ("loose", "sentence after dw(y₁), §14.6", "'at location y (< y₁) along the lifting line'",
         "the velocity is induced at y₁; it is downward for a filament at y > y₁ on the starboard half (as the caption of Fig. 14.21 says)",
         "the formula −(dΓ/dy)dy/(4π(y − y₁)) and (14.13) are consistent with the caption", "LL.staircase_downwash", "—"),
        ("loose", "(14.8)", "'camber = ≅ ½βc' (two relation signs); coefficient 1.3",
         "camber ≃ ½βc; the first-order thickness coefficient is 3√3/4 = 1.299, and the rule is first order in e",
         "measure the contour: at e = 0.1 the ratio t_max/(c e) is about 9 % below 1.3",
         "AF.zhukhovsky_shape(exact=True/False), AF.THICKNESS_RULE_COEFFICIENT", "a = b(1 + e) cos β instead of /cos β"),
        ("slip", "historical notes, §14.3", "Lanchester's year of birth", "1868",
         "a cited biography; the printed year would make him a teenager when his car company was founded", "—", "—"),
        ("loose", "text after (14.26)", "(4.107) and (4.108) 'in Section 4.3'", "§4.11 (dynamic similarity)",
         "the dimensionless coefficients are defined in the similarity section of chapter 4", "core.similarity", "—"),
        ("loose", "(14.22) against (14.18)–(14.20)", "upper limit N in (14.22), ∞ elsewhere", "either, consistently",
         "(14.21)–(14.22) are written for the truncated solution", "LL.induced_drag_from_series", "—"),
        ("loose", "text after (14.11)", "'chord length c (≃ 4a)' and 'camber/chord ratio β/2' beside c ≃ 4b of (14.8)",
         "both are first-order statements; exactly C_L = 8π(a/c) sin(α + β)",
         "at finite thickness the lift slope exceeds 2π (about 9 % at e = 0.1)",
         "AF.zhukhovsky_lift_coefficient(exact=True)", "—"),
        ("loose", "Exercises 14.14 and 14.18", "span written b", "s (b is the map constant)", "notation drift inside the exercises",
         "LL (span is always s)", "—"),
        ("loose", "end of §14.6", "the 'Λ-pattern' of migrating birds, printed with the aspect-ratio symbol", "a V formation",
         "same glyph as the aspect ratio Λ of (14.1)", "LL.horseshoe_velocity (upwash outboard of the tips)", "—"),
    ]
    out = {}
    for i, (kind, where, printed, correct, how, coded, test) in enumerate(rows, start=1):
        out[str(i)] = dict(id=i, kind=kind, where=where, printed=printed, correct=correct, corrected=correct,
                           how_to_tell=how, coded_in=coded, test=test)
    return out


def traps() -> list:
    """Things that are true as printed but easy to misread — kept apart from :func:`book_slips`.

    Book: §14.2–§14.9.   Returns a list of 18 dicts (``id`` "T1" … "T18", ``what``, ``where`` — the section, ``code``
    — how the code guards against it).
    Assumptions: none.   Validation: V1 a table.  Label: analytic.
    """
    rows = [
        ("Γ is the clockwise circulation; positive lift at positive α needs Γ_cw > 0", "§14.3, §14.5",
         "keyword-only Gamma_cw=; returned values named Gamma_cw"),
        ("dw/dz = u − iv here; ch13's Ekman velocity was u + iv", "§14.5", "ZhukhovskyFlow.velocity returns (u, v)"),
        ("b, a, β, c, e, Λ, γ change meaning or are overloaded", "all", "conventions_table(); e_span for the span efficiency"),
        ("three different θ: polar angle of ζ (14.4); from the stream-aligned diameter (14.10); stroke and course angles",
         "§14.4, §14.5, §14.8, §14.9", "contour parameter is t; rear_stagnation_point returns theta, phi and offset"),
        ("axes switch: x is the mapped chord in §14.4–§14.5 (lift is not F_y); y is the span in §14.6 and w is positive downward",
         "§14.4–§14.6", "zhukhovsky_force returns (Fx, Fy) and (L, D); horseshoe_velocity returns uz (up) and horseshoe_downwash −uz"),
        ("'the circulation does not change under the map' and 'equal velocities at infinity' are stated, not shown; the velocity is "
         "not equal anywhere else", "§14.5", "ZhukhovskyFlow.dwdz applies the factor ζ²/(ζ² − b²) of (14.9)"),
        ("which circle parameters: a = b(1 + e) is defined for the symmetric case only, a = b/cos β for the arc only", "§14.4",
         "zhukhovsky_circle: μ = b − a e^{−iβ}, a = b(1 + e)/cos β; also accepts (a, b, β)"),
        ("(14.9) is 0/0 at the trailing edge with the Kutta circulation and infinite otherwise; the inverse map needs the root "
         "outside the shifted circle", "§14.5", "analytic limit trailing_edge_speed; zhukhovsky_inverse picks the larger |ζ − μ|"),
        ("'induced drag is the largest part of the total drag' holds at high lift coefficient (exactly half at best glide), not at cruise",
         "§14.6, §14.7", "drag_share, drag_vs_speed, best_glide"),
        ("(14.13) is a principal value; the downwash of a single horseshoe is infinite at its tips; its 'better than elliptic' drag "
         "is an artefact", "§14.6", "singularity subtracted in downwash; collocation stations avoid the tips; horseshoe_wing note"),
        ("elliptic loading is not elliptic planform: twist can give it on any planform", "§14.6",
         "twist_for_elliptic_loading, elliptic_twist_rectangular"),
        ("K is per radian in (14.16)–(14.20) and quoted per degree in §14.7", "§14.6, §14.7", "K is a required argument, per radian"),
        ("A is the planform area in (4.107)–(4.108) as used here; a section uses the chord per unit span", "§14.3, §14.6",
         "section_coefficient (chord) vs lift_coefficient (area)"),
        ("the Kutta condition is a model of a viscous outcome, not a theorem of ideal flow; it fails at stall and in fast unsteady motion",
         "§14.3", "Gamma_cw can be set to any value; generic_lift_curve is labelled qualitative"),
        ("Kelvin's theorem needs only the circuit itself to stay in inviscid constant-density fluid; the foil and its boundary layers are inside",
         "§14.3", "impulsive_start_lumped_vortex keeps Γ_bound + ΣΓ_wake = 0 at every step"),
        ("the separation estimate, the wake roll-up past the first turn, the generic curves and the Reynolds-number regimes are "
         "qualitative models of ours", "§14.3, §14.6, §14.7", "label 'qualitative' in the docstring and in the returned dict"),
        ("typical values quoted in the prose and every exercise input are private", "all", "illustrative_inputs()"),
        ("the fish and sail formulas are quasi-steady and drag-free in the book; 'ideally perfect' is an identity of the velocity triangle",
         "§14.8, §14.9", "oscillating_foil(D=…) and sail_forces(CD=…) include drag"),
    ]
    return [dict(id=f"T{i}", what=w, where=wh, code=c) for i, (w, wh, c) in enumerate(rows, start=1)]


# ======================================================================================================================
# §14.1–§14.2: incompressibility check, level flight
# ======================================================================================================================
def compressibility_check(U, T: float = 288.15) -> dict:
    """Is the flow slow enough to be treated as incompressible?  Mach number and the size of the density change.

    Book: §14.1 ("incompressible aerodynamics if the flow speeds are low enough (Mach number < 0.3)"), with (4.111)
    Ma = U/c.  Ours: the isentropic density change between the stream and a stagnation point is ≈ ½Ma² to leading
    order (4.5 % at Ma = 0.3).
    Parameters: U speed [m/s]; T air temperature [K] (default 288.15 K, sea level).
    Returns dict: ``Ma`` [–], ``density_change`` ≈ ½Ma² [–], ``incompressible_ok`` (Ma < 0.3; bool or array).
    Assumptions: perfect gas (dry air), isentropic stagnation; leading order in Ma².
    Validation: V1 against the exact isentropic density ratio (1e-3 at Ma = 0.3: series truncation); V7 Ma → 0.
    Label: analytic.
    """
    Ma = _F(mach_number(U, T=T))  # Eq. (4.111)
    ok = Ma < 0.3
    return {"Ma": _S(Ma), "density_change": _S(0.5 * Ma ** 2), "incompressible_ok": bool(ok) if ok.ndim == 0 else ok}


def level_flight(W, rho, U, A, CD=None) -> dict:
    """Steady level flight: the lift coefficient the wing must fly at, and the drag and power if C_D is given.

    Book: §14.2 (lift "must be generated by the wings to counter the weight", thrust is negative drag) with (4.108):
    L = W ⇒ C_L = W/(½ρU²A); T = D = C_D ½ρU²A; power = D U.
    Parameters: W weight [N]; rho [kg/m³]; U [m/s]; A planform area [m²]; CD drag coefficient [–] or None.
    Returns dict: ``CL`` [–], ``q`` dynamic pressure [Pa], ``D`` [N], ``power`` [W] (None without CD).
    Assumptions: steady, level, unaccelerated; thrust along the flight path.
    Validation: V1 L recomputed from C_L equals W (1e-12); V2 dimensions.  Label: analytic.
    """
    q = 0.5 * _F(rho) * _F(U) ** 2
    CL = _F(W) / (q * _F(A))  # Eq. (4.108) with L = W
    D = None if CD is None else _S(_F(CD) * q * _F(A))  # Eq. (4.107)
    return {"CL": _S(CL), "q": _S(q), "D": D, "power": None if CD is None else _S(_F(D) * _F(U))}


def stall_speed(W, rho, A, CL_max):
    """Slowest speed of level flight, U_min = √(2W/(ρ A C_L,max)) [m/s].

    Book: §14.7 (flaps raise the lift coefficient, so the aircraft can land more slowly — stated).  Ours: L = W with
    (4.108) at C_L = C_L,max.
    Parameters: W [N]; rho [kg/m³]; A planform area [m²]; CL_max [–].   Returns U_min [m/s].
    Assumptions: steady level flight.   Validation: V1 inverts :func:`level_flight`.  Label: analytic.
    """
    return _S(np.sqrt(2.0 * _F(W) / (_F(rho) * _F(A) * _F(CL_max))))


# ======================================================================================================================
# scalar wrappers for explainer parity rows (no new physics)
# ======================================================================================================================
def surface_point(U: float, b: float, e: float, beta: float, alpha: float, t: float, *,
                  Gamma_cw: float | None = None) -> dict:
    """One point of the circle and its image on the section: position, the two factors of (14.9), speed and C_p.

    Book: §14.5, Eq. (14.9) dw/dz = (dw/dζ) ζ²/(ζ² − b²), at the circle point ζ = μ + a e^{it} (t = polar angle about
    the centre Q from the ξ direction; t = −β is the trailing-edge point B); C_p = 1 − |dw/dz|²/U² (Bernoulli).
    A scalar wrapper round ``zhukhovsky_flow(...).dwdzeta`` and ``.dwdz_from_zeta`` — no new physics.
    Parameters: U [m/s]; b [m]; e [–]; beta, alpha [rad]; t [rad]; Gamma_cw [m²/s] (None → the Kutta value (14.11)).
    Returns dict of floats: ``zeta_re``, ``zeta_im`` [m]; ``x``, ``y`` [m] (14.3); ``dwdzeta_abs`` [m/s];
    ``map_factor_abs`` = |ζ²/(ζ² − b²)| [–] (inf at the edge); ``speed`` [m/s]; ``cp`` [–].  At t = −β the speed is
    the limit U(b/a)cos(α + β) of ``trailing_edge_speed`` with the Kutta circulation and inf otherwise.
    Assumptions: steady ideal flow, attached.
    Validation: V1 equals ``zhukhovsky_surface`` at its nodes; speed = dwdzeta_abs · map_factor_abs away from the edge.
    Label: analytic.
    """
    F = AF.zhukhovsky_flow(U, b, e, beta, alpha, Gamma_cw=Gamma_cw)
    tt = float(t)
    d = (tt + F.beta + np.pi) % (2.0 * np.pi) - np.pi  # angular distance from the trailing-edge point
    zeta = complex(F.b) if abs(d) < 1e-13 else complex(F.mu + F.a * np.exp(1j * tt))
    z = complex(joukowski(zeta, F.b))
    den = abs(zeta ** 2 - F.b ** 2)
    speed = float(abs(complex(F.dwdz_from_zeta(zeta))))  # Eq. (14.9)
    return {"zeta_re": zeta.real, "zeta_im": zeta.imag, "x": z.real, "y": z.imag,
            "dwdzeta_abs": float(abs(complex(F.dwdzeta(zeta)))),
            "map_factor_abs": float(abs(zeta) ** 2 / den) if den > 0.0 else float("inf"), "speed": speed,
            "cp": float(1.0 - (speed / F.U) ** 2)}


def edge_speed(U: float, b: float, e: float, beta: float, alpha: float, d: float, *,
               Gamma_cw: float | None = None) -> float:
    """Speed at the distance d behind the trailing edge, along the direction e^{−2iβ} in which the cusp points [m/s].

    Book: §14.5, text after (14.9): the velocity at the trailing edge is infinite unless B is a stagnation point of
    the circle flow; then it is finite.  Equals ``zhukhovsky_flow(U, b, e, beta, alpha, Gamma_cw=…).speed_near_edge(d)``
    — a scalar wrapper, no new physics.
    Parameters: U [m/s]; b [m]; e [–]; beta, alpha [rad]; d [m] (> 0); Gamma_cw [m²/s] (None → Kutta).
    Returns the speed [m/s]: → U(b/a)cos(α + β) as d → 0 with the Kutta value, ∝ d^{−1/2} otherwise.
    Assumptions: steady ideal flow.   Validation: V1 limit against ``trailing_edge_speed``.  Label: analytic.
    """
    return float(AF.zhukhovsky_flow(U, b, e, beta, alpha, Gamma_cw=Gamma_cw).speed_near_edge(float(d)))


def map_point(zeta_re: float, zeta_im: float, b: float) -> dict:
    """Image of one ζ-plane point under the Zhukhovsky map, by the real and imaginary parts (14.4).

    Book: §14.4, Eq. (14.4) x = (R + b²/R) cos θ, y = (R − b²/R) sin θ for ζ = R e^{iθ} (from z = ζ + b²/ζ (14.3)).
    Parameters: zeta_re, zeta_im [m] (not both zero); b [m].
    Returns dict of floats: ``R`` [m], ``theta`` [rad] (polar angle of ζ from the origin), ``x``, ``y`` [m].
    Assumptions: none (geometry).   Validation: V1 equals ``joukowski(zeta, b)``.  Label: analytic.
    """
    R = float(np.hypot(zeta_re, zeta_im))
    th = float(np.arctan2(zeta_im, zeta_re))
    return {"R": R, "theta": th, "x": float((R + b ** 2 / R) * np.cos(th)),  # Eq. (14.4)
            "y": float((R - b ** 2 / R) * np.sin(th))}  # Eq. (14.4)


# ======================================================================================================================
# §14.3–§14.4: the sharp edges (Exercises 14.2 and 14.7)
# ======================================================================================================================
def edge_flow_suction(rho: float, Uo: float, a: float, eps: float = 1e-3, n: int = 2000) -> dict:
    """Leading-edge suction on a plate of no thickness, from a control volume round its edge.

    Book: §14.3 ("the leading edge suction that occurs on foils of negligible thickness") and Exercise 14.7:
    φ = 2U_o√(ar) cos(θ/2) is the flow round the edge of a plate lying on the positive x-axis; the control-volume
    momentum equation ∮ρu(u·n)dξ = −∮p n dξ + F on a circle of radius ε gives the force F that holds the plate.
    Result (ours): u_r = U_o√(a/r) cos(θ/2), u_θ = −U_o√(a/r) sin(θ/2), p − p_∞ = −½ρU_o²a/r; the momentum flux is
    (πρU_o²a, 0) and the pressure integral vanishes, for every ε.  So F = +πρU_o²a x̂ must be applied to hold the
    plate: the fluid pulls the plate with πρU_o²a in the −x direction — out of its own edge, a suction.
    Parameters: rho [kg/m³]; Uo velocity scale [m/s]; a length scale [m]; eps radius of the control circle [m];
    n nodes of the periodic trapezoid rule in θ.
    Returns dict: ``F_hold`` (Fx, Fy) [N/m]; ``suction`` (force of the fluid on the plate, (−πρU_o²a, 0)) [N/m];
    ``momentum_flux``, ``pressure_force`` (the two integrals) [N/m]; ``closed_form`` = πρU_o²a [N/m].
    For a flat plate of chord c at incidence α the edge flow has U_o²a = U²c sin²α, which gives S = πρU²c sin²α of
    :func:`fluidpy.core.airfoils.flat_plate_forces`.
    Assumptions: steady ideal flow; plate of zero thickness.
    Validation: V1 quadrature = closed form; V4 independent of ε.  Label: analytic.
    """
    th = (np.arange(int(n)) + 0.5) * 2.0 * np.pi / int(n)
    amp = float(Uo) * np.sqrt(float(a) / float(eps))
    ur, ut = amp * np.cos(0.5 * th), -amp * np.sin(0.5 * th)
    ux = ur * np.cos(th) - ut * np.sin(th)
    uy = ur * np.sin(th) + ut * np.cos(th)
    dxi = float(eps) * 2.0 * np.pi / int(n)
    mom = (float(rho) * np.sum(ux * ur) * dxi, float(rho) * np.sum(uy * ur) * dxi)
    pg = -0.5 * float(rho) * (ur ** 2 + ut ** 2)  # p − p_∞ (Bernoulli; the fluid is at rest far away)
    pf = (-np.sum(pg * np.cos(th)) * dxi, -np.sum(pg * np.sin(th)) * dxi)
    Fh = (float(mom[0] - pf[0]), float(mom[1] - pf[1]))
    return {"F_hold": Fh, "suction": (-Fh[0], -Fh[1]), "momentum_flux": (float(mom[0]), float(mom[1])),
            "pressure_force": (float(pf[0]), float(pf[1])), "closed_form": float(np.pi * rho * Uo ** 2 * a)}


def trailing_edge_wedge(n: float, U: float = 1.0, d: float = 1.0) -> dict:
    """Corner-flow model of a trailing edge of finite included angle: φ = (Ud/n)(r/d)ⁿ cos(nθ).

    Book: §14.4, Fig. 14.15 and Exercise 14.2 (ch06's power-law potentials w = Azⁿ).  The fluid fills |θ| < π/n (the
    rays θ = ±π/n are streamlines: the two surfaces), so the body's included angle is 2π(1 − 1/n); the speed is
    U(r/d)^{n−1}: n = 1 a cusp (uniform flow, no stagnation point), 1 < n ≤ 2 a wedge with a stagnation point at the
    edge.
    Parameters: n exponent (1 ≤ n ≤ 2); U speed at r = d [m/s]; d reference length [m].
    Returns dict: ``included_angle`` [rad]; ``half_fluid_angle`` = π/n [rad]; ``velocity(r, theta)`` → (u_r, u_θ)
    [m/s]; ``speed(r)`` [m/s]; ``Cp(r)`` = 1 − (r/d)^{2(n−1)} (reference: the speed at r = d); ``Cp_edge``.
    Assumptions: steady ideal flow near the edge.
    Validation: V1 n = 1 uniform flow, C_p(edge) = 0; n > 1: C_p(edge) = 1; u_θ = 0 on θ = ±π/n.  Label: analytic.
    """
    n, U, d = float(n), float(U), float(d)

    def velocity(r, theta):
        f = U * (_F(r) / d) ** (n - 1.0)
        return _S(f * np.cos(n * _F(theta))), _S(-f * np.sin(n * _F(theta)))

    def speed(r):
        return _S(U * (_F(r) / d) ** (n - 1.0))

    def Cp(r):
        return _S(1.0 - (_F(r) / d) ** (2.0 * (n - 1.0)))

    return {"included_angle": float(2.0 * np.pi * (1.0 - 1.0 / n)), "half_fluid_angle": float(np.pi / n),
            "velocity": velocity, "speed": speed, "Cp": Cp, "Cp_edge": 0.0 if n == 1.0 else 1.0, "n": n}


# ======================================================================================================================
# §14.3, §14.7: stall (qualitative)
# ======================================================================================================================
def chord_reynolds(U, c, nu):
    """Chord Reynolds number Re_c = Uc/ν [–].

    Book: §14.7 (the maximum lift coefficient depends on the Reynolds number).
    Parameters: U [m/s]; c chord [m]; nu kinematic viscosity [m²/s].   Label: analytic.
    """
    return _S(_F(U) * _F(c) / _F(nu))


def stall_regime(Re_c: float, *, Re_laminar: float, Re_turbulent: float) -> str:
    """Which of the two stall behaviours §14.7 describes a chord Reynolds number falls in — a labelled classifier.

    Book: §14.7 (at moderate Reynolds number the suction-side layer separates while still laminar: early stall, low
    maximum lift; at high Reynolds number it becomes turbulent first and stays attached longer: later stall, higher
    maximum lift — the same physics as the drag crisis of ch09).  The book's ranges are private; the two thresholds
    are required keywords.
    Returns "laminar separation first" (Re_c < Re_laminar), "transition first" (Re_c > Re_turbulent) or "mixed".
    Assumptions: a two-threshold cartoon of a continuous dependence.   Validation: V7 both sides of each threshold.
    Label: qualitative.
    """
    if Re_c < Re_laminar:
        return "laminar separation first"
    if Re_c > Re_turbulent:
        return "transition first"
    return "mixed"


def suction_side_separation(alpha: float, e: float, beta: float, Re_c: float, n: int = 400, side: str = "suction",
                            closure: str = "falkner_skan") -> dict:
    """Where a laminar boundary layer would separate on a Zhukhovsky section — a qualitative estimate of ours.

    Book: §14.3 and §14.7 (stall: the suction-side boundary layer cannot climb the adverse pressure gradient behind
    the suction peak and separates ahead of the trailing edge; Fig. 14.26) — described, not computed.  Ours: Thwaites'
    method of §9.6, (9.44)–(9.50) (``core.boundary_layer.thwaites``), on the ideal surface speed of
    ``zhukhovsky_surface`` from the front stagnation point to the trailing edge.
    **Qualitative**: the layer is taken laminar and the outer flow attached (Kutta); a real section at flight
    Reynolds numbers becomes turbulent and stays attached much further, and after separation the pressure itself
    changes.  The estimate shows the trend (the point moves forward as α grows), not the maximum lift.  In Thwaites'
    method the separation *position* does not depend on Re_c (only the layer's thickness does).
    Parameters
    ----------
    alpha [rad]; e [–] (> 0: a rounded nose); beta [rad]; Re_c chord Reynolds number [–]; n surface points;
    side "suction" (default), "pressure", "upper" or "lower"; closure — Thwaites closure.
    Returns
    -------
    dict: ``x_sep_over_c`` (chordwise position of separation from the leading edge along the chord axis; NaN if the
    layer stays attached to the edge); ``separated`` (bool); ``s_sep`` arc length from the stagnation point [in units
    of b = 1]; ``z_sep`` (complex, b = 1); ``xc_peak`` (position of the suction peak on that side); ``theta_sep``
    momentum thickness at separation / chord; ``label`` = "qualitative"; ``side``; ``xc``, ``Ue`` (the path).
    Assumptions: laminar, steady, attached ideal outer flow; b = 1, U = 1 (non-dimensional).
    Validation: V7 moves forward monotonically with α; symmetric section at α = 0: same on both sides.
    Label: qualitative.
    """
    if e <= 0.0:
        raise ValueError("a section with a rounded nose is needed (e > 0): the plate's leading-edge speed is infinite")
    b, U = 1.0, 1.0
    F = AF.zhukhovsky_flow(U, b, e, beta, alpha)
    a, mu = F.a, F.mu
    t_front = (np.pi + 2.0 * (alpha + beta)) % (2.0 * np.pi)  # u_θ = 0 of (14.10) at θ = π + (α + β)
    lifting_up = (alpha + beta) >= 0.0
    upper = {"suction": lifting_up, "pressure": not lifting_up, "upper": True, "lower": False}[side]
    t = np.linspace(t_front, 0.0 if upper else 2.0 * np.pi, int(n))
    zeta = mu + a * np.exp(1j * (t - beta))
    jac = np.abs(1.0 - b ** 2 / zeta ** 2)
    u_th = -2.0 * U * np.sin(t - beta - alpha) - F.Gamma_cw / (2.0 * np.pi * a)  # Eq. (14.10)
    with np.errstate(divide="ignore", invalid="ignore"):
        Ue = np.abs(u_th) / jac  # Eq. (14.9) on the circle
    Ue[0] = 0.0
    Ue[-1] = float(AF.trailing_edge_speed(U, a, b, alpha, beta))
    z = np.asarray(joukowski(zeta, b))
    ds = a * jac * abs(t[1] - t[0])
    s = np.concatenate([[0.0], np.cumsum(0.5 * (ds[:-1] + ds[1:]))])
    sh = AF.zhukhovsky_shape(b, e, beta, exact=True, n=4001)
    c = sh["chord"]
    x_all = AF.zhukhovsky_airfoil(b, e, beta, 2001)["x"]
    xmin, cx = float(x_all.min()), float(x_all.max() - x_all.min())
    xc = (z.real - xmin) / cx
    r = thwaites(s, Ue, nu=U * c / float(Re_c), closure=closure, stop_at_separation=False)
    ssep = r["x_sep"]
    k = int(np.argmax(Ue))
    out = {"separated": ssep is not None, "s_sep": float("nan"), "x_sep_over_c": float("nan"), "z_sep": complex("nan"),
           "theta_sep": float("nan"), "xc_peak": float(xc[k]), "label": "qualitative", "side": "upper" if upper else "lower",
           "xc": xc, "Ue": Ue, "s": s, "lam": r["lam"]}
    if ssep is not None:
        ts = float(np.interp(ssep, s, t))
        zs = complex(joukowski(mu + a * np.exp(1j * (ts - beta)), b))
        out.update(s_sep=float(ssep), x_sep_over_c=float((zs.real - xmin) / cx), z_sep=zs,
                   theta_sep=float(np.interp(ssep, s, r["theta"]) / c))
    return out


def generic_lift_curve(alpha, *, K: float, beta: float, alpha_stall: float, CL_post: float, width: float):
    """A drawn model of a measured lift curve: the line K(α + β) up to stall, then a smooth drop to a plateau.

    Book: §14.5, Fig. 14.17 and §14.7, Fig. 14.25 (measured C_L against α: a linear range, a maximum, a drop at stall,
    a post-stall plateau) — figures, not formulas; their numbers are private.  Ours:
    C_L = (1 − σ) K(α + β) + σ C_L,post with σ = smoothstep((α − α_stall)/width) (0 below α_stall, 1 above
    α_stall + width).  Every constant is a required keyword.
    Parameters: alpha [rad]; K [1/rad]; beta [rad]; alpha_stall [rad] (where the curve leaves the line);
    CL_post plateau value [–]; width [rad] of the drop.   Returns C_L [–].
    Assumptions: an empirical cartoon; positive stall only.
    Validation: V7 equals K(α + β) exactly for α ≤ α_stall; continuous with a continuous slope.  Label: qualitative.
    """
    al = _F(alpha)
    x = np.clip((al - float(alpha_stall)) / float(width), 0.0, 1.0)
    sig = x * x * (3.0 - 2.0 * x)
    return _S((1.0 - sig) * float(K) * (al + float(beta)) + sig * float(CL_post))


def generic_drag_curve(alpha, *, CD0: float, AR: float, e_span: float, K: float, beta: float, alpha_stall: float,
                       CL_post: float, width: float, CD_sep: float):
    """A drawn model of a measured drag curve: profile drag + induced drag + a rise where the flow separates.

    Book: §14.7, Fig. 14.25 (drag rising roughly as α² below stall, then steeply) — a figure; numbers private.
    Ours: C_D = C_D0 + C_L(α)²/(πΛ e_span) + C_D,sep σ(α), with C_L and σ of :func:`generic_lift_curve`.
    Parameters: alpha [rad]; CD0 [–]; AR [–]; e_span [–]; K, beta, alpha_stall, CL_post, width as in
    :func:`generic_lift_curve`; CD_sep extra drag coefficient of the separated flow [–].   Returns C_D [–].
    Assumptions: an empirical cartoon.   Validation: V7 minimum C_D0 at zero lift (α = −β).  Label: qualitative.
    """
    al = _F(alpha)
    x = np.clip((al - float(alpha_stall)) / float(width), 0.0, 1.0)
    sig = x * x * (3.0 - 2.0 * x)
    CL = _F(generic_lift_curve(al, K=K, beta=beta, alpha_stall=alpha_stall, CL_post=CL_post, width=width))
    return _S(float(CD0) + CL ** 2 / (np.pi * float(AR) * float(e_span)) + float(CD_sep) * sig)


# ======================================================================================================================
# §14.3: the starting vortex (a lumped-vortex model, ours)
# ======================================================================================================================
def impulsive_start_lumped_vortex(U, c: float, alpha, t_end: float, dt: float, *, wake: str = "rigid",
                                  shed_fraction: float = 0.25, eps: float | None = None, rho: float = 1.0,
                                  history: bool = False) -> dict:
    """Growth of the bound circulation of a flat plate after a sudden start, with the shed (starting) vorticity.

    Book: §14.3, Figs. 14.9–14.10: by Kelvin's theorem the circulation round a material circuit that started in
    fluid at rest stays zero, Γ_foil + Γ_starting = 0; a new vortex is shed at every change of speed or incidence.
    The book gives no model for the transient.  Ours — the lumped-vortex model: one bound vortex at the quarter
    chord, no flow through the plate at the three-quarter-chord point, one new wake vortex per time step placed
    ``shed_fraction``·U dt behind the trailing edge, and Γ_bound + ΣΓ_wake = 0 imposed at every step (two linear
    equations per step).  The wake then moves with the stream (``wake="rigid"``) or with the stream plus the velocity
    induced by all vortices (``wake="free"``, desingularised with ``eps``; explicit Euler, first order in dt).
    The steady limit is exactly the Kutta value of a flat plate, Γ = πUc sin α ((14.11) with a = b = c/4).
    Parameters
    ----------
    U : stream speed [m/s] — a float or a callable U(t).   c : chord [m].   alpha : incidence [rad] — a float or a
    callable α(t).   t_end [s]; dt [s] (both open choices of the model).   wake : "rigid" | "free".
    shed_fraction [–]; eps core size [m] (None → 0.05c); rho [kg/m³]; history : also return the wake positions at
    every step (arrays padded with NaN).
    Returns
    -------
    dict: ``t`` [s]; ``s_chords`` = ∫U dt/c; ``Gamma_bound_cw`` [m²/s]; ``Gamma_wake_cw`` (sum of the shed vortices);
    ``total`` (their sum: zero to rounding); ``Gamma_steady`` = πU(t)c sin α(t); ``ratio`` = bound/steady;
    ``lift_circulatory`` = ρUΓ_bound [N/m] (the added-mass part of the unsteady lift is not modelled);
    ``wake_x``, ``wake_y`` [m] (plate frame: leading edge at 0, trailing edge at c), ``wake_Gamma_cw`` at the end;
    with ``history``: ``wake_x_hist``, ``wake_y_hist`` (steps × steps).
    Assumptions: ideal flow outside a thin wake; small plate (one bound vortex); a model, not a simulation of the
    boundary layer — the first value is 1/3 of the steady one where Wagner's exact thin-airfoil result is 1/2.
    Validation: V4 total = 0 at every step (1e-12); V1 Γ_bound → πUc sin α; V3 first order in dt; V5 growth against
    :func:`wagner_function` (a few per cent after the first two chords, model).  Label: conserved, analytic (limit);
    the transient is a model.
    """
    Uf = U if callable(U) else (lambda t: float(U))
    af = alpha if callable(alpha) else (lambda t: float(alpha))
    c = float(c)
    nt = int(round(float(t_end) / float(dt)))
    xb, xc = 0.25 * c, 0.75 * c
    eps = 0.05 * c if eps is None else float(eps)
    xw, yw, gw = np.zeros(nt), np.zeros(nt), np.zeros(nt)  # counterclockwise-positive strengths internally
    tt, Gb, Gwsum, sch, Gst = (np.zeros(nt) for _ in range(5))
    hx = np.full((nt, nt), np.nan) if history else None
    hy = np.full((nt, nt), np.nan) if history else None
    dist = 0.0
    kb = 1.0 / (2.0 * np.pi * (xc - xb))
    for k in range(nt):
        t = (k + 1) * float(dt)
        Uk, ak = float(Uf(t)), float(af(t))
        ux, uy = Uk * np.cos(ak), Uk * np.sin(ak)
        xn, yn = c + shed_fraction * ux * dt, shed_fraction * uy * dt
        kn = (xc - xn) / (2.0 * np.pi * ((xc - xn) ** 2 + yn ** 2))
        old = (gw[:k] * (xc - xw[:k]) / (2.0 * np.pi * ((xc - xw[:k]) ** 2 + yw[:k] ** 2))).sum()
        # no flow through the plate at 3c/4, and Kelvin's theorem:
        M = np.array([[kb, kn], [1.0, 1.0]])
        rhs = np.array([-uy - old, -gw[:k].sum()])
        gb, gnew = np.linalg.solve(M, rhs)
        xw[k], yw[k], gw[k] = xn, yn, gnew
        m = k + 1
        tt[k], Gb[k], Gwsum[k] = t, -gb, -gw[:m].sum()
        dist += Uk * dt
        sch[k] = dist / c
        Gst[k] = np.pi * Uk * c * np.sin(ak)
        vx, vy = np.full(m, ux), np.full(m, uy)
        if wake == "free":
            P = np.vstack([xw[:m], yw[:m]])
            ui = point_vortex_velocity(P, P, gw[:m], eps) + point_vortex_velocity(P, np.array([[xb], [0.0]]), [gb], eps)
            vx, vy = vx + ui[0], vy + ui[1]
        elif wake != "rigid":
            raise ValueError('wake must be "rigid" or "free"')
        xw[:m] += vx * dt
        yw[:m] += vy * dt
        if history:
            hx[k, :m], hy[k, :m] = xw[:m], yw[:m]
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(Gst != 0.0, Gb / Gst, np.nan)
    out = {"t": tt, "s_chords": sch, "Gamma_bound_cw": Gb, "Gamma_wake_cw": Gwsum, "total": Gb + Gwsum,
           "Gamma_steady": Gst, "ratio": ratio, "lift_circulatory": float(rho) * np.array([Uf(q) for q in tt]) * Gb,
           "wake_x": xw, "wake_y": yw, "wake_Gamma_cw": -gw, "model": "lumped vortex (ours)"}
    if history:
        out.update(wake_x_hist=hx, wake_y_hist=hy)
    return out


def wagner_function(s):
    """Wagner's function φ(s): growth of the circulatory lift of a thin airfoil after a step in incidence (a fit).

    Book: none (§14.3 describes the start-up in words: the steady pattern is reached "after the flow has moved a
    chord or two").  R. T. Jones' two-exponential approximation, φ = 1 − 0.165 e^{−0.0455 s} − 0.335 e^{−0.3 s}, with
    s = Ut/(c/2) the distance travelled in **semi-chords**; φ(0) = ½, φ → 1.
    Source: BYU FLOW Lab, Aeroelasticity.jl documentation, "Wagner's Function"
    (https://flow.byu.edu/Aeroelasticity.jl/dev/models/aerodynamics/wagner; coefficients read from the fetched page on
    2026-10-12 by the implementer — the verifier should confirm first-hand before using it as V5).
    Parameters: s distance in semi-chords [–] (≥ 0).   Returns φ [–].
    Assumptions: thin airfoil, small incidence, incompressible; an approximation of the exact function (about 1 %).
    Validation: V7 φ(0) = 0.5, φ(∞) = 1, monotone.  Label: benchmark (cited fit; verifier to confirm the source).
    """
    s = _F(s)
    return _S(1.0 - 0.165 * np.exp(-0.0455 * s) - 0.335 * np.exp(-0.3 * s))


# ======================================================================================================================
# §14.6: the wake
# ======================================================================================================================
def wake_rollup(Gn, s: float, n: int = 200, t_end: float | None = None, eps: float | None = None, n_t: int = 41,
                fast: bool = False, rtol: float = 1e-7, atol: float = 1e-9) -> dict:
    """Roll-up of the trailing vortex sheet behind a wing, as desingularised point vortices in the cross-flow plane.

    Book: §14.6, Figs. 14.19–14.20 (the sheet rolls up into two tip vortices) — described.  Ours: the sheet of
    strength −dΓ/dy (14.18) is cut into ``n`` filaments at equal steps of the span angle, each a plane point vortex
    of circulation −ΔΓ about +x (counterclockwise on the starboard side, seen from behind), evolved by
    ``core.biot_savart.point_vortex_evolve`` with the core size ``eps`` (vortex blobs).  Time is a stand-in for
    distance behind the wing, x = Ut.
    **Qualitative past the first turn**: the turns of the spiral depend on n and eps; what is robust — and returned —
    are the invariants: total circulation of each half Γ₀, the vertical impulse ΣΓ_k y_k (so the spacing of the two
    centroids, b′ = (πs/4)Γ₁/Γ₀) and the descent of the centroid.
    Parameters
    ----------
    Gn [m²/s]; s [m]; n number of vortices (even; ``fast`` → 80); t_end [s] (None → 0.6 s²/Γ₀);
    eps [m] (None → 0.05 s); n_t output times; rtol, atol of the integrator.
    Returns
    -------
    dict: ``t`` (n_t,) [s]; ``y``, ``z`` (n_t, n) [m] (rear view: y to starboard, z up); ``sigma`` (n,) [m²/s];
    ``Gamma0``; ``impulse`` (n_t,) = Σσ_k y_k [m³/s]; ``half_circulation`` (starboard); ``centroid_y`` (n_t,) of the
    starboard half [m] (theory b′/2); ``centroid_z`` (n_t,) [m]; ``spacing_theory`` b′ [m]; ``descent_speed_theory``
    [m/s]; ``energy`` (n_t,) [m⁴/s²] — the Hamiltonian of the blob dynamics, −(1/4π)Σ_{j<k}σ_jσ_k ln(r_jk² + eps²);
    ``eps``; ``label``.
    Assumptions: two-dimensional unsteady analogy of the steady wake (slender wake); inviscid; planar initial sheet.
    Validation: V4 impulse and blob Hamiltonian conserved (to the integrator's tolerance); centroid spacing b′.  Label: conserved
    (invariants); qualitative (the spiral).
    """
    g = np.atleast_1d(_F(Gn))
    s = float(s)
    n = 80 if fast else int(n)
    ge = np.linspace(0.0, np.pi, n + 1)
    Ge = _F(LL.circulation_series(g, ge))
    sigma = Ge[:-1] - Ge[1:]  # −ΔΓ on each panel: circulation about +x
    y0 = -0.5 * s * np.cos(0.5 * (ge[:-1] + ge[1:]))
    pair = LL.rolled_up_pair(g, s)
    G0 = pair["Gamma0"]
    t_end = 0.6 * s ** 2 / abs(G0) if t_end is None else float(t_end)
    eps = 0.05 * s if eps is None else float(eps)
    t = np.linspace(0.0, t_end, int(n_t))
    X = point_vortex_evolve(np.vstack([y0, np.zeros(n)]), sigma, t, eps=eps, rtol=rtol, atol=atol)
    y, z = X[:, 0, :], X[:, 1, :]
    star = y0 > 0.0
    Gs = float(sigma[star].sum())
    # Hamiltonian of the desingularised (blob) dynamics: H_δ = −(1/4π) Σ_{j<k} σ_jσ_k ln(r_jk² + δ²) — conserved exactly
    dy = y[:, :, None] - y[:, None, :]
    dzz = z[:, :, None] - z[:, None, :]
    iu = np.triu_indices(n, 1)
    Hd = -np.sum((sigma[:, None] * sigma[None, :])[iu] * np.log((dy ** 2 + dzz ** 2)[:, iu[0], iu[1]] + eps ** 2), axis=1) / (4.0 * np.pi)
    return {"t": t, "y": y, "z": z, "sigma": sigma, "Gamma0": G0, "impulse": (sigma * y).sum(axis=1),
            "half_circulation": Gs, "centroid_y": (sigma[star] * y[:, star]).sum(axis=1) / Gs,
            "centroid_z": (sigma[star] * z[:, star]).sum(axis=1) / Gs, "spacing_theory": pair["spacing"],
            "descent_speed_theory": pair["descent_speed"], "energy": Hd,
            "eps": eps, "label": "qualitative past the first turn; invariants conserved"}


def tip_vortex_estimate(W: float, rho: float, U: float, s: float) -> dict:
    """Strength, spacing and sink rate of the tip-vortex pair behind an aircraft with an elliptic loading.

    Book: §14.6 (tip vortices with the mid-span circulation Γ₀, Fig. 14.20; the cover photograph) with (14.21),
    (14.24): L = W = (πs/4)ρUΓ₁ ⇒ Γ₀ = Γ₁ = 4W/(πρUs); spacing b′ = πs/4 and descent speed Γ₀/(2πb′) (ours, from
    the conserved impulse); downwash at the wing w = Γ₁/2s.
    Parameters: W weight [N]; rho [kg/m³]; U flight speed [m/s]; s span [m].
    Returns dict: ``Gamma0`` [m²/s], ``spacing`` [m], ``descent_speed`` [m/s], ``downwash`` [m/s],
    ``induced_drag`` [N] (14.25).
    Assumptions: steady level flight, elliptic loading, complete roll-up.
    Validation: V1 consistent with ``rolled_up_pair`` and (14.25).  Label: analytic.
    """
    G0 = 4.0 * float(W) / (np.pi * float(rho) * float(U) * float(s))  # Eq. (14.21) with L = W
    bp = 0.25 * np.pi * float(s)
    return {"Gamma0": float(G0), "spacing": float(bp), "descent_speed": float(G0 / (2.0 * np.pi * bp)),
            "downwash": float(LL.elliptic_downwash(G0, s)), "induced_drag": float(LL.elliptic_induced_drag(W, rho, U, s))}


def vortex_pair_energy_per_length(Gamma, s, sigma, rho):
    """Kinetic energy per unit length of a counter-rotating pair of Rankine vortices, (ρΓ²/2π)[ln(s/σ) + ¼] [J/m].

    Book: §14.6, Exercise 14.19 (tip-vortex core size and induced drag; Rankine cores (3.28)) — result ours: outside
    the cores the pair's energy is (ρΓ²/2π) ln(s/σ); each core in solid-body rotation holds ρΓ²/(16π).  A larger core
    holds less energy per unit length, and that energy per length is the induced drag (D_iU = energy rate).
    Parameters: Gamma circulation of each vortex [m²/s]; s distance between the two centres [m] (the pair's spacing —
    not the wing span); sigma core radius [m]; rho [kg/m³].   Returns the energy per length [J/m = N].
    Assumptions: σ ≪ s (cores do not deform each other); two-dimensional far wake.
    Validation: V7 grows as ln(s/σ); V2 dimensions.  Label: analytic (leading order in σ/s).
    """
    G = _F(Gamma)
    return _S(_F(rho) * G ** 2 / (2.0 * np.pi) * (np.log(_F(s) / _F(sigma)) + 0.25))


def wake_vortex_rolling_moment(rho: float, U: float, K: float, s: float, c: float, Gamma_v: float, ell: float) -> dict:
    """Rolling moment on a small rectangular wing flying along the axis of another aircraft's tip vortex.

    Book: §14.6, Exercise 14.16: M = ∫_{−s/2}^{s/2} ρU y Γ(y) dy with the vertical velocity
    w(y) = (Γ′/2πy)[1 − exp(−|y|/ℓ)] (**upward positive in this exercise — the opposite of the chapter's downwash**),
    the wing's own trailing vortices ignored.  Result (ours): strip theory, Γ(y) = ½KUc·(w/U) from (14.17) (the
    uniform part of the angle of attack gives an even Γ and no moment), so
    M = (ρUKcΓ′/2π)[s/2 − ℓ(1 − e^{−s/2ℓ})].
    Parameters: rho [kg/m³]; U speed of the small aircraft [m/s]; K its section slope [1/rad]; s span [m]; c chord
    [m]; Gamma_v circulation Γ′ of the large aircraft's vortex [m²/s]; ell core length scale ℓ [m].
    Returns dict: ``M`` [N m] (closed form), ``M_quad`` (quadrature), ``ratio`` to the ℓ → 0 limit ρUKcΓ′s/(4π).
    Assumptions: w ≪ U; lifting-strip theory without the small wing's own downwash (over-estimates the moment).
    Validation: V1 closed form = quadrature (1e-10); V7 odd in Γ′; ℓ → 0 limit.  Label: analytic.
    """
    pref = float(rho) * float(U) * float(K) * float(c) * float(Gamma_v) / (2.0 * np.pi)
    M = pref * (0.5 * s - ell * (-np.expm1(-0.5 * s / ell)))

    def f(y):
        return 0.0 if y == 0.0 else -np.expm1(-abs(y) / ell)

    Mq = pref * quad(f, 0.0, 0.5 * s, epsabs=1e-13, epsrel=1e-13)[0]
    return {"M": float(M), "M_quad": float(Mq), "ratio": float(M / (pref * 0.5 * s))}


# ======================================================================================================================
# §14.7: the drag breakdown and the polar
# ======================================================================================================================
def drag_polar(CL, *, CD0: float, AR: float, e_span: float = 1.0):
    """Drag polar C_D = C_D0 + C_L²/(πΛ e_span) [–]: profile drag plus induced drag.

    Book: §14.7 (drag = profile drag (friction + form) + induced drag) with (14.26) C_Di = C_L²/(πΛ).  The constant
    C_D0 and e_span are ours (the book states the breakdown in words).
    Parameters: CL [–]; CD0 profile-drag coefficient [–]; AR = Λ [–]; e_span [–].   Returns C_D [–].
    Assumptions: profile drag independent of lift (a parabolic polar); below stall.
    Validation: V1 parts of :func:`drag_breakdown` sum to it.  Label: analytic.
    """
    return _S(float(CD0) + _F(CL) ** 2 / (np.pi * float(AR) * float(e_span)))  # Eq. (14.26) + profile drag


def drag_share(CL, *, CD0: float, AR: float, e_span: float = 1.0):
    """Fraction of the total drag that is induced drag, C_Di/C_D [–].

    Book: §14.6–§14.7 ("the induced drag generally makes the largest contribution to the total drag" — trap T9: true
    at high lift coefficient; exactly ½ at the best glide; smaller at cruise).
    Parameters: as :func:`drag_polar`.   Returns the fraction in [0, 1).
    Assumptions: parabolic polar.   Validation: V1 = ½ at ``best_glide``'s C_L.  Label: analytic.
    """
    cdi = _F(CL) ** 2 / (np.pi * float(AR) * float(e_span))
    return _S(cdi / (float(CD0) + cdi))


def drag_breakdown(CL, *, CD0: float, AR: float, e_span: float = 1.0, friction_fraction: float | None = None) -> dict:
    """The drag of a wing split into named parts that sum to the total.

    Book: §14.7: total drag = friction drag (shear stress) + pressure drag (normal stress); pressure drag = induced
    drag + form drag; profile drag = friction + form (the drag the section alone would have).
    Parameters: CL [–]; CD0 profile-drag coefficient [–]; AR; e_span; friction_fraction — the share of the profile
    drag that is friction [–] (None → the friction/form split is not reported).
    Returns dict of coefficients: ``profile``, ``induced``, ``total``, ``induced_share``; with friction_fraction also
    ``friction``, ``form``, ``pressure`` (= form + induced).
    Assumptions: parabolic polar; the friction/form split is an input, not computed.
    Validation: V1 friction + form + induced = total; profile + induced = total (1e-12).  Label: analytic.
    """
    cdi = _F(CL) ** 2 / (np.pi * float(AR) * float(e_span))  # Eq. (14.26)
    out = {"profile": float(CD0), "induced": _S(cdi), "total": _S(float(CD0) + cdi), "induced_share": _S(cdi / (float(CD0) + cdi))}
    if friction_fraction is not None:
        fr = float(friction_fraction) * float(CD0)
        out.update(friction=fr, form=float(CD0) - fr, pressure=_S(float(CD0) - fr + cdi))
    return out


def profile_drag_estimate(Re_c, t_over_c: float, *, regime: str, form_factor) -> float:
    """Rough profile-drag coefficient of a section: flat-plate friction on both faces times a form factor.

    Book: §14.7 (profile drag = friction + form drag) — no formula.  Ours: C_D0 ≈ FF · 2 C_F(Re_c) with C_F the
    one-face flat-plate drag coefficient of ``core.boundary_layer.plate_drag_coefficient`` ((9.33) laminar; the
    turbulent one-fifth-power fit, whose coefficient is still secondary-sourced — an open item carried from ch09).
    The form factor is **an input** (a number, or a callable of t/c): no correlation is quoted from memory.
    Parameters: Re_c [–]; t_over_c thickness ratio [–]; regime "laminar" | "turbulent" | "mixed"; form_factor ≥ 1.
    Returns C_D0 [–] (based on the planform area).
    Assumptions: attached flow, zero-incidence plate friction; no pressure-gradient effect on the friction.
    Validation: V7 form_factor = 1 → twice the plate value.  Label: approximate (qualitative).
    """
    ff = float(form_factor(t_over_c)) if callable(form_factor) else float(form_factor)
    return float(ff * 2.0 * float(plate_drag_coefficient(Re_c, regime=regime, sides=1)))


def best_glide(CD0: float, AR: float, e_span: float = 1.0) -> dict:
    """Best lift-to-drag ratio of a parabolic polar: (L/D)_max = ½√(πΛ e_span/C_D0) at C_L = √(πΛ e_span C_D0).

    Book: §14.7 (C_Di ∝ C_L² ∝ α²; the polar is not written).  Ours: d(C_D/C_L)/dC_L = 0 ⇒ C_D0 = C_Di — at the best
    glide the induced drag is exactly half of the total.
    Parameters: CD0 [–]; AR [–]; e_span [–].
    Returns dict: ``CL``, ``CD`` (= 2C_D0), ``LD_max``, ``induced_share`` (= 0.5), ``glide_angle`` [rad].
    Assumptions: parabolic polar; below stall.
    Validation: V1 numerical minimisation of C_D/C_L (1e-6).  Label: analytic.
    """
    k = np.pi * float(AR) * float(e_span)
    CL = float(np.sqrt(k * float(CD0)))
    LD = float(0.5 * np.sqrt(k / float(CD0)))
    return {"CL": CL, "CD": 2.0 * float(CD0), "LD_max": LD, "induced_share": 0.5, "glide_angle": float(np.arctan(1.0 / LD))}


def drag_vs_speed(W: float, rho: float, A: float, AR: float, CD0: float, U, e_span: float = 1.0) -> dict:
    """Drag of an aircraft in level flight against speed: profile drag rising as U², induced drag falling as 1/U².

    Book: §14.6–§14.7 with (14.25): D_i = 2L²/(πρU²s²) (÷ e_span), L = W, s² = ΛA; profile drag ½ρU²A C_D0 (ours).
    Parameters: W [N]; rho [kg/m³]; A [m²]; AR; CD0; U speed(s) [m/s]; e_span.
    Returns dict: ``profile``, ``induced``, ``total`` [N]; ``CL``; ``induced_share``; ``LD``.
    Assumptions: steady level flight; parabolic polar; below stall (check C_L against C_L,max).
    Validation: V1 minimum of ``total`` at :func:`speed_min_drag` where the two parts are equal.  Label: analytic.
    """
    U = _F(U)
    q = 0.5 * float(rho) * U ** 2
    Dp = q * float(A) * float(CD0)
    Di = 2.0 * float(W) ** 2 / (np.pi * float(rho) * U ** 2 * float(AR) * float(A) * float(e_span))  # Eq. (14.25)
    return {"profile": _S(Dp), "induced": _S(Di), "total": _S(Dp + Di), "CL": _S(float(W) / (q * float(A))),
            "induced_share": _S(Di / (Dp + Di)), "LD": _S(float(W) / (Dp + Di))}


def speed_min_drag(W: float, rho: float, A: float, AR: float, CD0: float, e_span: float = 1.0) -> float:
    """Speed of least drag (best glide) in level flight, U* = √(2W/(ρA C_L*)) with C_L* = √(πΛ e_span C_D0) [m/s].

    Book: §14.7 (ours, from :func:`best_glide` and L = W).   Parameters: as :func:`drag_vs_speed`.
    Assumptions: parabolic polar.   Validation: V1 minimiser of ``drag_vs_speed``.  Label: analytic.
    """
    return float(np.sqrt(2.0 * float(W) / (float(rho) * float(A) * best_glide(CD0, AR, e_span)["CL"])))


# ======================================================================================================================
# §14.8: thrust from lift
# ======================================================================================================================
def oscillating_foil(U, V, L, D=0.0) -> dict:
    """A foil moving sideways at V while the animal moves forward at U: thrust, side force and the two powers.

    Book: §14.8, Fig. 14.27: the lift L is normal to the relative velocity U_r and has the forward component
    L sin θ, working at the rate UL sin θ; the tail pushes sideways against L cos θ at the rate VL cos θ; since
    V/U = tan θ the two rates are equal ("the conversion of energy is ideally perfect").  Ours: with a drag D along
    U_r, thrust = L sin θ − D cos θ and side force = L cos θ + D sin θ, so the efficiency
    η = (L sin θ − D cos θ)U/((L cos θ + D sin θ)V) is below 1 and falls to 0 at tan θ = D/L.
    Parameters: U forward speed [m/s]; V sideways speed of the foil [m/s]; L lift [N]; D drag [N].
    Returns dict: ``theta`` [rad]; ``Ur`` [m/s]; ``thrust``, ``side_force`` [N]; ``power_out`` = thrust·U,
    ``power_in`` = side force·V [W]; ``efficiency`` [–].
    Assumptions: quasi-steady; the same triangle on both half-strokes; no wake losses beyond D.
    Validation: V1 D = 0: efficiency 1 for any θ; V7 efficiency → 0 at tan θ = D/L.  Label: analytic.
    """
    U_, V_, L_, D_ = _F(U), _F(V), _F(L), _F(D)
    th = np.arctan2(V_, U_)  # V/U = tan θ
    thrust = L_ * np.sin(th) - D_ * np.cos(th)
    side = L_ * np.cos(th) + D_ * np.sin(th)
    pin, pout = side * V_, thrust * U_
    with np.errstate(divide="ignore", invalid="ignore"):
        eta = np.where(pin != 0.0, pout / pin, np.nan)
    return {"theta": _S(th), "Ur": _S(np.hypot(U_, V_)), "thrust": _S(thrust), "side_force": _S(side),
            "power_out": _S(pout), "power_in": _S(pin), "efficiency": _S(eta)}


# ======================================================================================================================
# §14.9: sailing against the wind
# ======================================================================================================================
def apparent_wind(U, V, theta) -> dict:
    """Apparent wind on a boat: U_r = U − V (true wind minus boat velocity).

    Book: §14.9, Fig. 14.29 (wind speed U, sailing speed V, apparent wind U_r).  Axes (ours): x along the boat's
    course, the true wind blowing *from* the direction at the angle θ off the bow (θ = 0 dead upwind, π/2 a beam
    reach), so U = −U(cos θ, sin θ), V = (V, 0).
    Parameters: U true wind speed [m/s]; V boat speed [m/s]; theta course angle off the true wind [rad].
    Returns dict: ``Ur`` apparent wind speed [m/s] = √(U² + V² + 2UV cos θ); ``beta_a`` apparent-wind angle off the
    bow [rad] (≤ θ: the boat's own motion swings the wind forward); ``vector`` (U_rx, U_ry) [m/s].
    Assumptions: steady motion on a straight course; no leeway in this triangle.
    Validation: V1 V = 0: apparent = true wind; V7 mirror symmetry θ → −θ.  Label: analytic.
    """
    U_, V_, th = _F(U), _F(V), _F(theta)
    urx, ury = -U_ * np.cos(th) - V_, -U_ * np.sin(th)
    return {"Ur": _S(np.hypot(urx, ury)), "beta_a": _S(np.arctan2(-ury, -urx)), "vector": (_S(urx), _S(ury))}


def sail_forces(U, V, theta, *, CL: float, CD: float, rho: float, A: float) -> dict:
    """Lift and drag of a sail in the apparent wind, resolved into the driving and lateral components.

    Book: §14.9, Fig. 14.29: lift perpendicular to U_r and drag parallel to U_r; the resultant F is split into a
    driving component along the boat's motion and a lateral one.  Ours: with β_a the apparent-wind angle off the
    bow, drive = L sin β_a − D cos β_a (zero at tan β_a = D/L), lateral = L cos β_a + D sin β_a (to leeward).
    Parameters: U, V [m/s]; theta [rad] (see :func:`apparent_wind`); CL, CD of the sail [–] (keywords); rho air
    density [kg/m³]; A sail area [m²].
    Returns dict: ``beta_a`` [rad], ``Ur`` [m/s], ``L``, ``D``, ``drive``, ``lateral`` (leeward +), ``F`` [N],
    ``drive_over_L`` [–].
    Assumptions: the sail is trimmed to the same C_L, C_D on every course; forces in the horizontal plane; quasi-steady.
    Validation: V1 drive = 0 exactly at tan β_a = D/L; drive² + lateral² = L² + D²; V7 port/starboard mirror.
    Label: analytic.
    """
    aw = apparent_wind(U, V, theta)
    Ur, ba = _F(aw["Ur"]), np.abs(_F(aw["beta_a"]))
    q = 0.5 * float(rho) * Ur ** 2 * float(A)
    L, D = q * float(CL), q * float(CD)
    drive = L * np.sin(ba) - D * np.cos(ba)
    lat = L * np.cos(ba) + D * np.sin(ba)
    with np.errstate(divide="ignore", invalid="ignore"):
        dol = np.where(L != 0.0, drive / L, np.nan)
    return {"beta_a": aw["beta_a"], "Ur": aw["Ur"], "L": _S(L), "D": _S(D), "drive": _S(drive), "lateral": _S(lat),
            "F": _S(np.hypot(L, D)), "drive_over_L": _S(dol)}


def sail_drive_curve(theta, *, U: float, V: float, CL: float, CD: float, rho: float, A: float):
    """Driving force of the sail against the course angle θ off the true wind [N] (:func:`sail_forces` over θ).

    Book: §14.9 ("the thrust decreases as the angle θ decreases and normally vanishes" below a limiting angle).
    Parameters: theta [rad] (array or scalar); the rest as :func:`sail_forces`.   Returns the drive [N].
    Assumptions: as :func:`sail_forces` (fixed boat speed V along the curve).   Label: analytic.
    """
    return sail_forces(U, V, theta, CL=CL, CD=CD, rho=rho, A=A)["drive"]


def velocity_made_good(V, theta):
    """Speed of progress straight upwind on a tack at the angle θ to the wind, V cos θ [m/s].

    Book: §14.9 (a zig-zag course at ±θ reaches an upwind destination) — ours.
    Parameters: V boat speed [m/s]; theta course angle off the true wind [rad].   Label: analytic.
    """
    return _S(_F(V) * np.cos(_F(theta)))


def sailboat_course_limit(LD_sail, LD_hull):
    """Closest steady course to the apparent wind: β_a = ε_A + ε_H, the two drag angles ε = tan⁻¹(D/L) [rad].

    Book: §14.9 (the keel is a second wing whose lift balances the sail's lateral force; the limiting angle is quoted,
    not derived).  Ours: in steady sailing the aerodynamic force (at ε_A behind the normal to the apparent wind) is
    balanced by the hydrodynamic force (at ε_H behind the normal to the course), which closes only when
    β_a = ε_A + ε_H.
    Parameters: LD_sail lift-to-drag ratio of the sail (and rig) [–]; LD_hull lift-to-drag ratio of hull and keel [–].
    Returns β_a,min [rad] — the angle off the *apparent* wind; the true-wind angle is larger (:func:`apparent_wind`).
    Assumptions: steady, forces in the horizontal plane, no heel.
    Validation: V1 :func:`sailboat_equilibrium` satisfies it.  Label: analytic.
    """
    return _S(np.arctan(1.0 / _F(LD_sail)) + np.arctan(1.0 / _F(LD_hull)))


def sailboat_equilibrium(U: float, theta: float, *, CL: float, CD: float, rho_air: float, A_sail: float,
                         rho_water: float, CDA_hull: float, LD_hull: float | None = None) -> dict:
    """Steady boat speed on a given course: the speed at which the sail's drive equals the hull's drag.

    Book: §14.9 (the driving component "performs work in moving the boat" against the hull's drag; the keel's lift
    opposes the lateral force).  Ours: hull drag ½ρ_w V² (C_D A)_hull; ``brentq`` on drive(V) − hull drag(V); the
    keel then has to supply the lateral force, which fixes the hull's lift-to-drag ratio and closes the force polygon
    (β_a = ε_A + ε_H).
    Parameters: U true wind [m/s]; theta course off the true wind [rad]; CL, CD of the sail; rho_air [kg/m³];
    A_sail [m²]; rho_water [kg/m³]; CDA_hull drag area of hull and keel along the course [m²]; LD_hull — the best
    lift-to-drag ratio the hull and keel can reach [–] (None → not checked).
    Returns dict: ``V`` [m/s]; ``sails`` (False if the drive is not positive at rest: "in irons"); ``beta_a`` [rad];
    ``drive``, ``lateral``, ``hull_drag``, ``keel_lift`` [N]; ``LD_hull_required``; ``eps_A``, ``eps_H`` [rad];
    ``closure`` = β_a − ε_A − ε_H (zero); ``vmg`` [m/s]; ``keel_holds`` (required |L/D| ≤ ``LD_hull``; True when
    ``LD_hull`` is None) — below the course limit of :func:`sailboat_course_limit` the keel cannot hold the boat.
    Assumptions: steady; hull drag quadratic in speed with a constant drag area; no heel, no waves.
    Validation: V1 polygon closes (closure 1e-10); drive = hull drag.  Label: analytic.
    """
    kw = dict(CL=CL, CD=CD, rho=rho_air, A=A_sail)
    hd = lambda V: 0.5 * float(rho_water) * V ** 2 * float(CDA_hull)  # noqa: E731
    f = lambda V: float(sail_forces(U, V, theta, **kw)["drive"]) - hd(V)  # noqa: E731
    if f(0.0) <= 0.0:
        r = sail_forces(U, 0.0, theta, **kw)
        return {"V": 0.0, "sails": False, "beta_a": r["beta_a"], "drive": r["drive"], "lateral": r["lateral"],
                "hull_drag": 0.0, "keel_lift": r["lateral"], "LD_hull_required": float("nan"), "eps_A": float(np.arctan2(CD, CL)),
                "eps_H": float("nan"), "closure": float("nan"), "vmg": 0.0, "keel_holds": True}
    hi = float(U)
    while f(hi) > 0.0:
        hi *= 2.0
    V = float(brentq(f, 0.0, hi, xtol=1e-13, rtol=1e-13))
    r = sail_forces(U, V, theta, **kw)
    eA, eH = float(np.arctan2(CD, CL)), float(np.arctan2(hd(V), r["lateral"]))
    need = float(r["lateral"] / hd(V))
    return {"V": V, "sails": True, "beta_a": r["beta_a"], "drive": r["drive"], "lateral": r["lateral"], "hull_drag": hd(V),
            "keel_lift": r["lateral"], "LD_hull_required": need, "eps_A": eA, "eps_H": eH,
            "closure": float(abs(r["beta_a"]) - eA - eH), "vmg": float(velocity_made_good(V, theta)),
            "keel_holds": True if LD_hull is None else bool(abs(need) <= float(LD_hull))}


# ======================================================================================================================
# Exercise 14.1 and the sympy engine of §14.4
# ======================================================================================================================
def projectile_range_drag(m: float, rho: float, A: float, CD: float, V0: float, theta0: float, g: float = G0,
                          numeric: bool = False):
    """Range of a projectile with a constant drag coefficient on a shallow trajectory.

    Book: Exercise 14.1 (formula given in the exercise): x ≃ (2m/ρAC_D) ln(1 + ρAC_D V_o² cos θ_o sin θ_o/(mg)); as
    C_D → 0, ln(1 + ε) → ε and x → 2V_o² cos θ_o sin θ_o/g, the range without drag.
    Parameters: m mass [kg]; rho air density [kg/m³]; A frontal area [m²]; CD [–]; V0 launch speed [m/s];
    theta0 launch angle [rad]; g [m/s²]; numeric → also integrate the full equations of the exercise
    (``solve_ivp``, RK45, rtol 1e-10, terminal event at y = 0).
    Returns the range [m] (float), or with ``numeric`` a dict: ``x_formula``, ``x_numeric``, ``x_vacuum`` [m],
    ``t_flight`` [s].
    Assumptions: V_x² ≫ V_y², mg ≫ D sin θ (shallow trajectory); level ground; no lift.
    Validation: V7 C_D → 0 limit; V3 the numerical range approaches the formula as θ_o → 0.  Label: analytic.
    """
    k = float(rho) * float(A) * float(CD)
    arg = k * V0 ** 2 * np.cos(theta0) * np.sin(theta0) / (m * g)
    x = 2.0 * m / k * np.log1p(arg) if k > 0 else 2.0 * V0 ** 2 * np.cos(theta0) * np.sin(theta0) / g
    if not numeric:
        return float(x)

    def rhs(t, q):
        vx, vy = q[2], q[3]
        d = 0.5 * k * np.hypot(vx, vy) / m
        return [vx, vy, -d * vx, -g - d * vy]

    def ground(t, q):
        return q[1]

    ground.terminal, ground.direction = True, -1
    sol = solve_ivp(rhs, (0.0, 1e4), [0.0, 0.0, V0 * np.cos(theta0), V0 * np.sin(theta0)], events=ground, rtol=1e-10,
                    atol=1e-12, first_step=1e-6)
    return {"x_formula": float(x), "x_numeric": float(sol.y[0, -1]),
            "x_vacuum": float(2.0 * V0 ** 2 * np.cos(theta0) * np.sin(theta0) / g), "t_flight": float(sol.t[-1])}


def circular_arc_sym() -> dict:
    """Sympy check of §14.4: (14.4) → (14.5), the cosine rule → (14.6), and the circle of the arc.

    Book: §14.4, Eqs. (14.4)–(14.6) and the unnumbered x² + (y + 2b cot 2β)² = (2b csc 2β)².
    Returns dict of sympy expressions: ``x``, ``y`` (14.4); ``res_14_5`` = x² sin²θ − y² cos²θ − 4b² sin²θ cos²θ;
    ``res_14_6`` — the cosine rule b²/cos²β = R² + b² tan²β − 2Rb tan β cos(90° − θ) divided by R, minus
    (R − b²/R − 2b tan β sin θ); ``res_arc`` — the circle equation with y = 2b tan β sin²θ (from (14.6)) and x² from
    (14.5); ``ok`` (bool; alias ``all_zero``).  All three residuals simplify to 0.
    Assumptions: 0 < θ < π, 0 < β < π/4, R, b > 0.   Validation: V2 exact.  Label: symbolic.
    """
    import sympy as sp

    R, b, th, be = sp.symbols("R b theta beta", positive=True)
    x = (R + b ** 2 / R) * sp.cos(th)  # Eq. (14.4)
    y = (R - b ** 2 / R) * sp.sin(th)  # Eq. (14.4)
    r5 = sp.simplify(sp.expand(x ** 2 * sp.sin(th) ** 2 - y ** 2 * sp.cos(th) ** 2 - 4 * b ** 2 * sp.sin(th) ** 2 * sp.cos(th) ** 2))
    cosine = R ** 2 + b ** 2 * sp.tan(be) ** 2 - 2 * R * b * sp.tan(be) * sp.sin(th) - b ** 2 / sp.cos(be) ** 2
    r6 = sp.simplify(sp.expand(cosine / R - (R - b ** 2 / R - 2 * b * sp.tan(be) * sp.sin(th))))
    ya = 2 * b * sp.tan(be) * sp.sin(th) ** 2  # Eq. (14.6): y/sin θ = 2b tan β sin θ
    x2 = ya ** 2 * sp.cos(th) ** 2 / sp.sin(th) ** 2 + 4 * b ** 2 * sp.cos(th) ** 2  # Eq. (14.5) solved for x²
    arc = x2 + (ya + 2 * b * sp.cot(2 * be)) ** 2 - (2 * b / sp.sin(2 * be)) ** 2
    ra = sp.simplify(sp.expand(sp.expand_trig(arc.rewrite(sp.sin).rewrite(sp.cos))))
    if ra != 0:  # second route: everything in tan(β) and sin θ
        tb, sn = sp.symbols("t_b s_n", positive=True)
        sub = arc.subs({sp.cot(2 * be): (1 - tb ** 2) / (2 * tb), sp.sin(2 * be): 2 * tb / (1 + tb ** 2), sp.tan(be): tb,
                        sp.cos(th) ** 2: 1 - sn ** 2, sp.sin(th): sn})
        ra = sp.simplify(sp.expand(sub))
    ok = bool(r5 == 0 and r6 == 0 and ra == 0)
    return {"x": x, "y": y, "res_14_5": r5, "res_14_6": r6, "res_arc": ra, "all_zero": ok, "ok": ok}
