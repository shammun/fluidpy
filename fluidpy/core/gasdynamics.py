"""Perfect-gas flow relations (Kundu, Cohen & Dowling 5e, Ch. 15 — "GD" in the analysis and curation files).

One-dimensional steady flow of a perfect gas with constant specific heats, written with one number, the Mach number M:
stagnation and sonic reference states (15.27)–(15.31), the area rule (15.32)–(15.34), normal shocks (15.35)–(15.44),
nozzles at different back pressures (§15.7), constant-area ducts with friction and heat (15.45)–(15.48), and the planar
waves of §15.9: Mach angle (15.49), oblique shocks (15.50)–(15.55) and the Prandtl–Meyer function (15.56).

Conventions (every function)
----------------------------
* SI units, kelvin, absolute pressures; angles in **radians**; ``gamma`` is always an argument (default 1.4); R and C_p
  enter only through :class:`PerfectGas` (``AIR_BOOK`` carries the rounded constants of (15.3); the project default is
  ``core.thermo.R_AIR`` ≈ 287.05 — pass ``AIR_BOOK`` to reproduce a book number).  **One gas per call**: a function
  that takes both ``gamma`` and ``gas`` (:func:`mass_flow_vs_back_pressure`) takes γ from ``gas`` when ``gas`` is
  given and raises ``ValueError`` if an explicit ``gamma`` disagrees with ``gas.gamma``.
* **Three sets of air constants live side by side** (trap T4) — they differ in the fourth significant figure:

  ==========================================  ==================  =====================================
  where                                        R [J/(kg K)]        C_p [J/(kg K)]
  ==========================================  ==================  =====================================
  ``core.bernoulli.CP_BOOK`` (ch04 default of  —                   1004.5
  ``stagnation_temperature(T, U, cp=…)``)
  ``PerfectGas()`` (project air, γ = 1.4)      287.058 (R_AIR)     γR/(γ − 1) = 1004.70
  ``AIR_BOOK`` (rounded constants of (15.3))   287.0               1004.5
  ==========================================  ==================  =====================================

  Consequence: ``stagnation_temperature(300, 400)`` (C_p = 1004.5) = 379.6416 K but
  300 × ``stagnation_temperature_ratio(mach_from_speed(400, 300))`` (PerfectGas(), C_p = 1004.70) = 379.6255 K — a
  1e-6 agreement check between the two routes **fails unless C_p is passed consistently**
  (``stagnation_temperature(T, U, cp=gas.cp)`` with the same ``gas`` used for M).
* **Station 1 is upstream, station 2 downstream.**  Subscript 0 = stagnation state, asterisk = sonic state reached
  *isentropically* (``A_over_Astar``).  The starred states of :func:`fanno_ratios` / :func:`rayleigh_ratios` are the
  sonic state of the *same frictional / heated process* — a different star, never to be mixed with A/A*.
* State dictionaries have the keys ``"rho"`` [kg/m³], ``"u"`` [m/s], ``"p"`` [Pa].
* The book's duct parameters are ``f = F/(p₁A)`` (F = total *retarding* wall force) and ``q = Q/h₁``; the number A of
  (15.46) is :func:`duct_parameter_A` (it is not an area).
* Oblique shock: ``sigma`` = shock angle, ``delta`` = deflection, both from the upstream flow direction (β and θ in
  most other sources).  :func:`prandtl_meyer_turn` counts an **expansion** turn positive; the book's δ (positive when
  the wall turns into the flow) is ``-turn``.
* Scalars in → Python floats out (explainer parity rows); arrays in → arrays out.  Outside its range a function returns
  NaN or raises ``ValueError`` with a stated meaning ("detached", "choked", M₁ < 1, A/A* < 1).

**Three book slips are coded corrected** (the printed forms stay available so that a test can show they fail):
(15.32) third member needs a minus sign (``printed_15_32``); (15.47) first numerator is 1 + (γ − 1)M²
(``printed_15_47``); (15.48) second relation carries (γ − 1)M² (``printed_15_48``) — see
:func:`influence_coefficients`.

Formulas labelled "ours" are not printed in the book (it gives a table column, a figure or an exercise instead); each
follows from the book equations named in its docstring.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.optimize import brentq

from ._util import as_scalar_if_0d
from .thermo import R_AIR

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d

__all__ = [
    "PerfectGas", "AIR_BOOK", "mach_from_speed",
    "mass_flow", "mass_flow_from_mach", "total_enthalpy", "energy_residual", "wall_force",
    "influence_coefficients", "area_velocity", "duct_role",
    "stagnation_temperature_ratio", "stagnation_pressure_ratio", "stagnation_density_ratio",
    "mach_from_pressure_ratio", "isentropic_state", "isentropic_table", "sonic_ratios", "critical_pressure_ratio",
    "max_speed", "pitot_compressibility_factor", "mach_for_density_change", "entropy_change_from_stagnation",
    "stagnation_pressure_loss", "area_ratio", "mach_from_area_ratio", "mass_flux_function", "choked_mass_flow",
    "laval_area", "nozzle_isentropic",
    "normal_shock", "normal_shock_mach", "normal_shock_pressure_ratio", "normal_shock_density_ratio",
    "normal_shock_temperature_ratio", "normal_shock_entropy", "weak_shock_entropy", "shock_admissible",
    "normal_shock_stagnation_pressure_ratio", "normal_shock_residuals", "normal_shock_fluxes", "normal_shock_table",
    "hugoniot_pressure_ratio", "strong_shock_limits", "moving_shock",
    "convergent_nozzle", "mass_flow_vs_back_pressure", "laval_critical_pressures", "laval_regime",
    "laval_shock_location", "laval_exit_state", "laval_nozzle",
    "duct_parameter_A", "duct_outlet_mach", "duct_state_ratios", "duct_entropy_change", "duct_residuals",
    "duct_choking_limits", "duct_curve", "duct_choking_locus", "rayleigh_temperature_sign", "fanno_ratios",
    "fanno_length_to_f", "rayleigh_ratios", "fanno_line_Ts", "rayleigh_line_Ts",
    "mach_angle", "mach_from_angle", "deflection_angle", "oblique_shock_components", "max_deflection",
    "sonic_deflection", "oblique_shock", "oblique_shock_curve", "detachment_mach", "weak_wave_pressure",
    "weak_wave_error", "prandtl_meyer", "prandtl_meyer_max", "mach_from_prandtl_meyer", "prandtl_meyer_turn",
    "riemann_exact",
]


def _vec(fn: Callable, *args):
    """Apply a scalar function element-wise (root finders) and return a float for scalar input."""
    return _S(np.vectorize(fn, otypes=[float])(*args))


# ======================================================================================================================
# §15.1 — the gas
# ======================================================================================================================
@dataclass(frozen=True)
class PerfectGas:
    """A perfect gas with constant specific heats: p = ρRT, C_p − C_v = R, γ = C_p/C_v.

    Book: §15.1, Eq. (15.1) (specific-heat line: C_v = R/(γ − 1), C_p = γR/(γ − 1); sound-speed line c = √(γRT)).

    Parameters
    ----------
    gamma : ratio of specific heats γ [–], > 1 (default 1.4).
    R : specific gas constant [J/(kg K)] (default ``core.thermo.R_AIR`` ≈ 287.05; the book's rounded 287 is
        :data:`AIR_BOOK`).

    Attributes ``cp``, ``cv`` [J/(kg K)]; method ``c(T)`` [m/s] for T [K].

    **Air constants (trap T4)**: ``PerfectGas()`` has R = 287.058, hence C_p = 1004.70; :data:`AIR_BOOK` has R = 287.0,
    C_p = 1004.5; ch04's ``stagnation_temperature`` defaults to ``CP_BOOK`` = 1004.5.  Mixing them moves results in
    the fourth or fifth significant figure (379.6416 K vs 379.6255 K for T = 300 K, u = 400 m/s), so a 1e-6 check
    between ``stagnation_temperature`` and T·:func:`stagnation_temperature_ratio` fails unless ``cp=gas.cp`` is
    passed with the same ``gas`` used for the Mach number.
    Assumptions: thermally and calorically perfect gas.
    Validation: V1, V2 — tests/test_ch15.py: test_perfect_gas_V1_specific_heats_and_sound_speed,
    test_sonic_ratios_and_choked_flow_V2_derivation.
    Label: analytic, symbolic.
    """

    gamma: float = 1.4
    R: float = R_AIR

    def __post_init__(self):
        if not (self.gamma > 1.0) or not (self.R > 0.0):
            raise ValueError("PerfectGas needs gamma > 1 and R > 0")

    @property
    def cp(self) -> float:
        """C_p = γR/(γ − 1) [J/(kg K)], Eq. (15.1)."""
        return self.gamma * self.R / (self.gamma - 1.0)  # Eq. (15.1)

    @property
    def cv(self) -> float:
        """C_v = R/(γ − 1) [J/(kg K)], Eq. (15.1)."""
        return self.R / (self.gamma - 1.0)  # Eq. (15.1)

    def c(self, T):
        """Speed of sound c = √(γRT) [m/s] at temperature T [K], Eq. (15.1)."""
        return _S(np.sqrt(self.gamma * self.R * _F(T)))  # Eq. (15.1)


#: Air with the rounded constants of (15.3): γ = 1.40, R = 287 m²/(s² K) (then C_p = 1004.5, C_v = 717.5; the book's
#: own rounded C_p = 1004 and C_v = 717 give C_p/C_v = 1.4003 — trap T4).
AIR_BOOK = PerfectGas(1.4, 287.0)


def _gas(gas: PerfectGas | None) -> PerfectGas:
    return PerfectGas() if gas is None else gas


def mach_from_speed(u, T, gas: PerfectGas | None = None):
    """Mach number M = u/c with c = √(γRT).

    Book: §15.1 ((4.111) restated) with the sound-speed line of Eq. (15.1).
    Parameters: u speed [m/s]; T static temperature [K]; gas (default air, γ = 1.4, R = R_AIR).
    Returns M [–].
    Assumptions: perfect gas; T is the *static* temperature at the point.
    Validation: V1 — tests/test_ch15.py: test_perfect_gas_V1_specific_heats_and_sound_speed.
    Label: analytic.
    """
    g = _gas(gas)
    return _S(_F(u) / np.sqrt(g.gamma * g.R * _F(T)))


# ======================================================================================================================
# §15.3 — the three duct laws
# ======================================================================================================================
def mass_flow(rho, u, A):
    """Mass flow rate ṁ = ρuA of a quasi-one-dimensional duct flow.

    Book: §15.3, Eq. (15.18)  ρ₁u₁A₁ = ρ₂u₂A₂, or ρuA = ṁ = const.
    Parameters: rho density [kg/m³]; u speed [m/s]; A cross-section [m²].   Returns ṁ [kg/s].
    Assumptions: steady flow; cross-section averages.
    Validation: V2 — tests/test_ch15.py: test_sonic_ratios_and_choked_flow_V2_derivation.
    Label: symbolic.
    """
    return _S(_F(rho) * _F(u) * _F(A))  # Eq. (15.18)


def mass_flux_function(M, gamma: float = 1.4):
    """Dimensionless mass flux ṁ√(RT₀)/(p₀A) = √γ M [1 + ½(γ − 1)M²]^{−(γ+1)/(2(γ−1))} (ours).

    Book: §15.4 — follows from (15.18) with (15.28)–(15.30): ρu = (p/RT) M √(γRT).  Maximum at M = 1.
    Parameters: M [–] ≥ 0; gamma [–].   Returns the flux function [–] (0.6847 at M = 1 for γ = 1.4).
    Assumptions: perfect gas; p₀, T₀ are the local stagnation values.
    Validation: V1, V2 — tests/test_ch15.py: test_area_ratio_V1_built_from_mass_conservation_and_its_two_inverses,
    test_sonic_ratios_and_choked_flow_V2_derivation.
    Label: analytic, symbolic.
    """
    M = _F(M)
    return _S(np.sqrt(gamma) * M * (1.0 + 0.5 * (gamma - 1.0) * M ** 2) ** (-(gamma + 1.0) / (2.0 * (gamma - 1.0))))


def mass_flow_from_mach(M, p0, T0, A, gas: PerfectGas | None = None):
    """Mass flow ṁ = p₀A/√(RT₀) × :func:`mass_flux_function` (M).

    Book: §15.3–§15.4, Eq. (15.18) with (15.28)–(15.30).
    Parameters: M [–]; p0 stagnation pressure [Pa]; T0 stagnation temperature [K]; A area [m²]; gas.
    Returns ṁ [kg/s].
    Assumptions: perfect gas; p0, T0 local stagnation values at the section.
    Validation: V2 — tests/test_ch15.py: test_sonic_ratios_and_choked_flow_V2_derivation.
    Label: symbolic.
    """
    g = _gas(gas)
    return _S(_F(p0) * _F(A) / np.sqrt(g.R * _F(T0)) * _F(mass_flux_function(M, g.gamma)))


def total_enthalpy(T, u, cp: float):
    """Stagnation (total) enthalpy h₀ ≡ h + ½u² = C_pT + ½u².

    Book: §15.4 (definition of h₀, unnumbered) and Eq. (15.27)  C_pT₀ = C_pT + ½u².
    Parameters: T static temperature [K]; u speed [m/s]; cp [J/(kg K)].   Returns h₀ [J/kg].
    Assumptions: perfect gas (h = C_pT).
    Validation: V2 — tests/test_ch15.py: test_perfect_gas_V2_dimensions_of_the_dimensional_wrappers.
    Label: symbolic.
    """
    return _S(cp * _F(T) + 0.5 * _F(u) ** 2)  # Eq. (15.27)


def _h(state: dict, gamma: float):
    if "h" in state:
        return _F(state["h"])
    return gamma / (gamma - 1.0) * _F(state["p"]) / _F(state["rho"])


def energy_residual(state1: dict, state2: dict, Q=0.0, gamma: float = 1.4):
    """Residual of the duct energy law, h₂ + ½u₂² − h₁ − ½u₁² − Q (zero when the law holds).

    Book: §15.3, Eq. (15.22)  h₂ + ½u₂² − h₁ − ½u₁² = Q (valid with friction and shocks between the stations).
    Parameters: state1, state2 dicts with ``"u"`` [m/s] and either ``"h"`` [J/kg] or ``"p"`` [Pa] and ``"rho"``
    [kg/m³] (then h = γp/((γ − 1)ρ)); Q heat added per unit mass between 1 and 2 [J/kg]; gamma [–].
    Returns the residual [J/kg].
    Assumptions: steady flow, no shaft work, no body force.
    Validation: V4 — tests/test_ch15.py: test_duct_laws_V4_mass_energy_entropy_along_an_isentropic_nozzle.
    Label: conserved.
    """
    return _S(_h(state2, gamma) + 0.5 * _F(state2["u"]) ** 2 - _h(state1, gamma) - 0.5 * _F(state1["u"]) ** 2
              - _F(Q))  # Eq. (15.22)


def wall_force(state1: dict, state2: dict, A1, A2):
    """x-force F exerted **on the fluid** by the duct walls between stations 1 and 2.

    Book: §15.3, Eq. (15.24)  −ṁu₁ + ṁu₂ = (pA)₁ − (pA)₂ + F.
    Parameters: state1, state2 dicts (``rho``, ``u``, ``p``); A1, A2 areas [m²].
    Returns F [N], positive in +x (the flow direction).  The retarding force of (15.45) is −F in a constant-area duct.
    Assumptions: steady flow; ṁ = ρ₁u₁A₁.
    Validation: V4 — tests/test_ch15.py: test_duct_laws_V4_mass_energy_entropy_along_an_isentropic_nozzle.
    Label: conserved.
    """
    mdot = _F(state1["rho"]) * _F(state1["u"]) * _F(A1)
    return _S(mdot * (_F(state2["u"]) - _F(state1["u"])) - _F(state1["p"]) * _F(A1) + _F(state2["p"]) * _F(A2))  # (15.24)


# ======================================================================================================================
# §15.5, §15.8 — differential (influence) coefficients: area, friction, heat
# ======================================================================================================================
def influence_coefficients(M, gamma: float = 1.4, *, dA_over_A=0.0, df=0.0, dq=0.0, printed_15_32: bool = False,
                           printed_15_47: bool = False, printed_15_48: bool = False) -> dict:
    """Fractional changes of ρ, u, p, T, M and entropy over a short length of duct with a small area change dA/A, a
    small friction increment df and a small heat increment dq.

    Book: §15.3 Eq. (15.19) dρ/ρ + du/u + dA/A = 0; §15.5 Eqs. (15.32)–(15.34); §15.8 Eqs. (15.47), (15.48).
    The four linearised laws at a station where p = p₁, T = T₁, h = h₁ (df = dF/(pA), dq = dQ/h):

    * mass (15.19)            dρ/ρ + du/u + dA/A = 0
    * momentum (15.25)/(15.45) dp/p + γM² du/u = −df
    * energy (15.22)/(15.45)   dT/T + (γ − 1)M² du/u = dq
    * state p = ρRT            dp/p = dρ/ρ + dT/T

    give du/u = (−dA/A + df + dq)/(1 − M²)  — (15.34) for the area, the second relation of (15.47) for friction,
    and u du/h₁ = (γ − 1)M² dq/(1 − M²) for heat; then dp/p = −γM² du/u − df, dT/T = dq − (γ − 1)M² du/u.

    **Corrected slips** (each printed form is an option that a test must fail):

    * ``printed_15_32=True`` — the page prints u du = −dp/ρ = **+**c² dρ/ρ; correct is −c² dρ/ρ (dp = c² dρ at
      constant entropy; the book's own next line (15.33) dρ/ρ = −M² du/u needs the minus).  The option evaluates the
      *area* part with dρ/ρ = +M² du/u, hence du/u = −(dA/A)/(1 + M²).  Keys it changes: ``du_over_u``,
      ``drho_over_rho``, ``dp_over_p``, ``dT_over_T``, ``dM_over_M``, ``u_du_over_h`` (not ``ds_over_cp``,
      ``dT0_over_T0``, ``dp0_over_p0``, which do not contain dA/A).
    * ``printed_15_47=True`` — the page prints −dp/p₁ = [1 − (γ − 1)M²]/(1 − M²) df; correct numerator is
      1 + (γ − 1)M².  The option replaces the friction part of ``dp_over_p`` and, through dp/p, changes
      ``ds_over_cp`` and ``dp0_over_p0`` (keys it changes: ``dp_over_p``, ``ds_over_cp``, ``dp0_over_p0``).
    * ``printed_15_48=True`` — the page prints u du/h₁ = (1 − γ)M²/(1 − M²) dq; correct is (γ − 1)M²/(1 − M²) dq
      (heating accelerates a subsonic flow).  The option flips the heat part of du/u and everything built on it —
      keys it changes: ``du_over_u``, ``drho_over_rho``, ``dp_over_p``, ``dM_over_M``, ``u_du_over_h``,
      ``ds_over_cp``, ``dT0_over_T0``, ``dp0_over_p0`` (all but ``dT_over_T``, whose heat part is taken from the
      energy law and never flipped).

    **Inlet-referred increments.**  The book's df and dq are referred to the *inlet* state (df_book = dF/(p₁A),
    dq_book = dQ/h₁), the local ones here to the local state (df = dF/(pA), dq = dQ/h); hence
    df_book = (p/p₁)·df_local and dq_book = (h/h₁)·dq_local.  That is why the second relation of (15.47),
    du/u = (p₁/p) df/(1 − M²), carries the factor p₁/p: written with the local df it is du/u = df/(1 − M²).  At the
    inlet itself (p = p₁, h = h₁) the two coincide.

    Parameters
    ----------
    M : local Mach number [–], ≥ 0.   gamma : [–].
    dA_over_A : fractional area change dA/A [–].   df : friction increment dF/(pA) [–] (≥ 0, retarding).
    dq : heat increment dQ/h [–] (either sign).

    Returns
    -------
    dict of [–] quantities: ``drho_over_rho``, ``du_over_u``, ``dp_over_p``, ``dT_over_T``, ``dM_over_M``,
    ``ds_over_cp``, ``u_du_over_h`` (= (γ − 1)M² du/u), ``dT0_over_T0``, ``dp0_over_p0``.

    **At the sonic point M = 1** (the driver D = −dA/A + df + dq):

    * ``du_over_u``, ``drho_over_rho``, ``dp_over_p``, ``dT_over_T``, ``dM_over_M``, ``u_du_over_h`` are **±inf**
      when D ≠ 0 — also when only one of the three increments is non-zero — and **NaN** when the driver D is exactly
      zero (in floating point) (a 0/0: the sonic
      point then needs the throat condition, which this local rule does not contain).  The sign of the infinity is
      that of the subsonic side (1 − M² evaluates to +0), as in :func:`area_velocity`.
    * ``ds_over_cp`` = dq + [(γ − 1)/γ] df, ``dT0_over_T0`` = dq/(1 + ½(γ − 1)M²) and ``dp0_over_p0`` =
      −df − ½γM² dq/(1 + ½(γ − 1)M²) do not contain 1/(1 − M²) at all: they stay **finite** at M = 1.
    * With a ``printed_*`` option the rule is the same for the driver of that printed form: ``printed_15_32`` makes
      the area part regular (−(dA/A)/(1 + M²)), so the driver is df + dq and a pure area change returns finite
      values at M = 1; ``printed_15_48`` makes the driver −dA/A + df − dq (and gives ``ds_over_cp``,
      ``dT0_over_T0``, ``dp0_over_p0`` a spurious pole when dq ≠ 0); ``printed_15_47`` adds a pole to
      ``ds_over_cp`` and ``dp0_over_p0`` when df ≠ 0 and changes the driver of ``dp_over_p``.

    Method: every key is formed as (driver)·1/(1 − M²) + (regular part) with the driver summed **before** the
    division, so an increment that is exactly zero cannot produce 0·inf = NaN.
    Assumptions: steady quasi-one-dimensional flow, perfect gas, first order in the increments.
    Validation: V2 (sympy solution of the four linearised laws, all nine keys), V1 (finite differences of (15.31)
    and of (15.46), Richardson), V7 (M → 0, sign pattern, the sonic point; the three printed forms fail the same
    checks) — tests/test_ch15.py: test_influence_coefficients_V2_sympy_solution_of_the_four_linear_laws,
    test_influence_area_V1_finite_differences_of_15_31_and_printed_15_32_fails,
    test_influence_friction_heat_V1_finite_differences_of_15_46_and_printed_15_47_15_48_fail,
    test_influence_coefficients_V7_limits_and_sign_pattern,
    test_influence_coefficients_V7_sonic_point_is_infinite_as_documented,
    test_influence_entropy_and_stagnation_changes_V2_derivation (ds/C_p, dT₀/T₀, dp₀/p₀ from the finite laws),
    test_influence_coefficients_V7_sonic_point_signs_zero_drivers_and_printed_forms (the M = 1 contract key by key).
    Label: analytic, symbolic.
    """
    M = _F(M)
    a, f_, q_ = _F(dA_over_A), _F(df), _F(dq)
    M2 = M ** 2
    g1 = gamma - 1.0
    zero = np.zeros(np.broadcast_shapes(M.shape, a.shape, f_.shape, q_.shape))
    with np.errstate(divide="ignore", invalid="ignore"):
        inv = 1.0 / (1.0 - M2)

        def extra(driver):  # a pole that only a printed (wrong) form has: absent when its increment is exactly zero
            driver = driver + zero
            return np.where(driver == 0.0, 0.0, driver * inv)

        tfac = 1.0 + 0.5 * g1 * M2
        q_48 = -q_ if printed_15_48 else q_  # printed (15.48): u du/h₁ = (1 − γ)M² dq/(1 − M²)
        if printed_15_32:
            D_u = f_ + q_48 + zero  # printed (15.32): dρ/ρ = +M² du/u into (15.19) — the area part has no pole
            D_T = f_ + q_ + zero
            du_reg = -a / (1.0 + M2) + zero
            no_pole = (f_ == 0.0) & (q_ == 0.0) & (zero == 0.0)  # nothing but the (regular) printed area part

            def pole(driver):
                return np.where(no_pole, 0.0, driver * inv)
        else:
            D_u = -a + f_ + q_48 + zero  # the driver of du/u
            D_T = -a + f_ + q_ + zero  # the driver that enters the energy law (the heat part of dT/T is never flipped)
            du_reg = zero

            def pole(driver):  # ±inf at M = 1; NaN (0/0) when the driver vanishes there
                return driver * inv
        d47 = 2.0 * g1 * M2 * f_ if printed_15_47 else zero  # printed (15.47): numerator 1 − (γ − 1)M²
        d48 = g1 * M2 * (q_ - q_48)  # (γ − 1)M² × (driver of the energy law − driver of du/u): 0 unless printed_15_48
        e47, e48 = extra(d47), extra(d48)
        du = pole(D_u) + du_reg  # Eqs. (15.34), (15.47) second relation, (15.48) second relation (corrected)
        drho = -a - du  # mass, Eq. (15.19); with dA = 0 this is (15.33) for the area part
        dp = pole(-gamma * M2 * D_u + d47) - gamma * M2 * du_reg - f_  # momentum; Eq. (15.47) first relation (corr.)
        dT = q_ - g1 * M2 * (pole(D_T) + du_reg)  # energy; Eq. (15.48) first relation, (1 − γM²)dq/(1 − M²)
        dM = pole(tfac * D_u + 0.5 * d48) + tfac * du_reg - 0.5 * q_  # dM/M = du/u − ½ dT/T
        ds = q_ + g1 / gamma * f_ - e48 - g1 / gamma * e47 + zero  # ds/C_p = dT/T − [(γ − 1)/γ] dp/p
        dT0 = (q_ - e48) / tfac + zero  # dT₀/T₀ = dT/T + (γ − 1)M² (dM/M)/(1 + ½(γ − 1)M²)
        dp0 = -f_ - 0.5 * gamma * M2 * (q_ - e48) / tfac + e47 + zero  # dp₀/p₀ = dp/p + γM² (dM/M)/(1 + ½(γ − 1)M²)
    out = dict(drho_over_rho=drho, du_over_u=du, dp_over_p=dp, dT_over_T=dT, dM_over_M=dM, ds_over_cp=ds,
               u_du_over_h=g1 * M2 * du, dT0_over_T0=dT0, dp0_over_p0=dp0)
    return {k: _S(v) for k, v in out.items()}


def area_velocity(M, dA_over_A):
    """Fractional speed change for a fractional area change in isentropic duct flow: du/u = −(dA/A)/(1 − M²).

    Book: §15.5, Eq. (15.34)  du/u = −[1/(1 − M²)] dA/A.
    Parameters: M [–]; dA_over_A [–].   Returns du/u [–] (±inf at M = 1 for dA ≠ 0).
    Assumptions: steady, frictionless, adiabatic (isentropic) quasi-one-dimensional flow.
    Validation: V1, V2, V7 — tests/test_ch15.py:
    test_influence_area_V1_finite_differences_of_15_31_and_printed_15_32_fails,
    test_influence_coefficients_V2_sympy_solution_of_the_four_linear_laws,
    test_influence_coefficients_V7_limits_and_sign_pattern.
    Label: analytic, symbolic.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        return _S(-_F(dA_over_A) / (1.0 - _F(M) ** 2))  # Eq. (15.34)


def duct_role(M: float, dA_sign: float) -> str:
    """Does this piece of duct act as a nozzle (speeds the flow up) or a diffuser (slows it down)?

    Book: §15.5, consequences (i), (ii) of Eq. (15.34); Fig. 15.4.
    Parameters: M local Mach number [–]; dA_sign sign of dA/dx (negative = converging, positive = diverging).
    Returns "nozzle", "diffuser", "none" (dA = 0) or "sonic" (M = 1: only possible at a throat).
    Assumptions: isentropic quasi-one-dimensional flow.
    Validation: V1, V7 — tests/test_ch15.py: test_small_helpers_V1_against_independent_closed_forms,
    test_influence_coefficients_V7_limits_and_sign_pattern.
    Label: analytic.
    """
    if dA_sign == 0:
        return "none"
    if M == 1.0:
        return "sonic"
    du = -np.sign(dA_sign) / (1.0 - float(M) ** 2)  # Eq. (15.34)
    return "nozzle" if du > 0 else "diffuser"


# ======================================================================================================================
# §15.4 — stagnation and sonic reference states
# ======================================================================================================================
def stagnation_temperature_ratio(M, gamma: float = 1.4):
    """T₀/T = 1 + ½(γ − 1)M².

    Book: §15.4, Eq. (15.28)  T₀/T = 1 + u²/(2C_pT) = 1 + ½(γ − 1)M².
    Parameters: M [–]; gamma [–].   Returns T₀/T [–] (≥ 1).
    Assumptions: perfect gas; adiabatic deceleration to rest (friction allowed).
    Validation: V1 — tests/test_ch15.py: test_stagnation_ratios_V1_from_the_energy_law_and_the_isentrope.
    Label: analytic.
    """
    return _S(1.0 + 0.5 * (gamma - 1.0) * _F(M) ** 2)  # Eq. (15.28)


def stagnation_pressure_ratio(M, gamma: float = 1.4):
    """p₀/p = [1 + ½(γ − 1)M²]^{γ/(γ−1)}.

    Book: §15.4, Eq. (15.29).
    Parameters: M [–]; gamma [–].   Returns p₀/p [–] (≥ 1).
    Assumptions: perfect gas; *isentropic* deceleration to rest.
    Validation: V1, V5, V2 — tests/test_ch15.py: test_stagnation_ratios_V1_from_the_energy_law_and_the_isentrope,
    test_normal_shock_V5_naca_1135_and_cited_table_rows, test_shock_stagnation_pressure_ratio_V2_derivation.
    Label: analytic, benchmark, symbolic.
    """
    return _S((1.0 + 0.5 * (gamma - 1.0) * _F(M) ** 2) ** (gamma / (gamma - 1.0)))  # Eq. (15.29)


def stagnation_density_ratio(M, gamma: float = 1.4):
    """ρ₀/ρ = [1 + ½(γ − 1)M²]^{1/(γ−1)}.

    Book: §15.4, Eq. (15.30).
    Parameters: M [–]; gamma [–].   Returns ρ₀/ρ [–] (≥ 1).
    Assumptions: perfect gas; isentropic deceleration to rest.
    Validation: V1, V2 — tests/test_ch15.py: test_stagnation_ratios_V1_from_the_energy_law_and_the_isentrope,
    test_pitot_compressibility_V2_series_and_incompressible_limit.
    Label: analytic, symbolic.
    """
    return _S((1.0 + 0.5 * (gamma - 1.0) * _F(M) ** 2) ** (1.0 / (gamma - 1.0)))  # Eq. (15.30)


def mach_from_pressure_ratio(p_over_p0, gamma: float = 1.4):
    """Mach number from the static-to-stagnation pressure ratio (inverse of (15.29)).

    Book: §15.4, Eq. (15.29) inverted: M² = [2/(γ − 1)] [(p/p₀)^{−(γ−1)/γ} − 1].
    Parameters: p_over_p0 [–] in (0, 1]; gamma [–].   Returns M [–] (NaN outside (0, 1]).
    Assumptions: as (15.29).
    Validation: V1, V5 — tests/test_ch15.py: test_stagnation_ratios_V1_from_the_energy_law_and_the_isentrope,
    test_stagnation_ratios_V5_naca_report_1135_rows.
    Label: analytic, benchmark.
    """
    r = _F(p_over_p0)
    with np.errstate(invalid="ignore", divide="ignore"):
        M2 = 2.0 / (gamma - 1.0) * np.expm1(-(gamma - 1.0) / gamma * np.log(r))
        return _S(np.where((r > 0) & (r <= 1.0), np.sqrt(np.maximum(M2, 0.0)), np.nan))


def area_ratio(M, gamma: float = 1.4):
    """Area ratio to the sonic throat, A/A* = (1/M) {[2/(γ + 1)] [1 + ½(γ − 1)M²]}^{(γ+1)/(2(γ−1))}.

    Book: §15.4, Eq. (15.31) (stated; Exercise 15.4).  A* is the area at which the *same isentropic* flow would be
    sonic — it changes across a shock, with friction and with heat.
    Parameters: M [–] > 0; gamma [–].   Returns A/A* [–] (≥ 1, minimum 1 at M = 1; inf at M = 0).
    Assumptions: steady isentropic quasi-one-dimensional flow of a perfect gas.
    Validation: V1, V2 — tests/test_ch15.py: test_area_ratio_V1_built_from_mass_conservation_and_its_two_inverses,
    test_shock_stagnation_pressure_ratio_V2_derivation.
    Label: analytic, symbolic.
    """
    M = _F(M)
    with np.errstate(divide="ignore", invalid="ignore"):
        return _S((2.0 / (gamma + 1.0) * (1.0 + 0.5 * (gamma - 1.0) * M ** 2)) ** (0.5 * (gamma + 1.0) / (gamma - 1.0))
                  / M)  # Eq. (15.31)


def _mach_from_area_scalar(ar: float, gamma: float, branch: str) -> float:
    if not np.isfinite(ar):
        return np.nan
    if ar < 1.0 - 1e-12:
        raise ValueError(f"A/A* = {ar!r} < 1: no isentropic flow passes through an area smaller than the sonic throat")
    if ar <= 1.0 + 1e-15:
        return 1.0
    fun = lambda m: float(area_ratio(m, gamma)) - ar  # noqa: E731
    if branch == "subsonic":
        lo = min(0.5, 0.1 / ar)
        root = brentq(fun, lo, 1.0, xtol=1e-15, rtol=1e-15)
    elif branch == "supersonic":
        hi = 2.0
        while fun(hi) < 0.0:
            hi *= 2.0
        root = brentq(fun, 1.0, hi, xtol=1e-15, rtol=1e-15)
    else:
        raise ValueError('branch must be "subsonic" or "supersonic"')
    assert abs(fun(root)) <= 1e-9 * ar, "area-ratio inversion did not converge"
    return float(root)


def mach_from_area_ratio(A_over_Astar, gamma: float = 1.4, *, branch: str):
    """The Mach number that belongs to an area ratio — on the subsonic or on the supersonic branch (ours).

    Book: §15.4, Eq. (15.31) inverted; the two branches are the two columns a reader finds in Table 15.1 (trap T8).
    Parameters: A_over_Astar [–] ≥ 1; gamma [–]; branch "subsonic" (M ≤ 1) or "supersonic" (M ≥ 1) — a **required
    keyword** (there is no default: the choice is physics, not a convenience).
    Returns M [–]; exactly 1 at A/A* = 1.  Raises ``ValueError`` for A/A* < 1.
    Method: ``brentq`` on each side of M = 1 (bracket from the physics), residual asserted.  Near A/A* = 1 the root is
    conditioned like √(A/A* − 1) (the minimum is quadratic).
    Assumptions: as :func:`area_ratio`.
    Validation: V1, V5 — tests/test_ch15.py: test_area_ratio_V1_built_from_mass_conservation_and_its_two_inverses,
    test_stagnation_ratios_V5_naca_report_1135_rows.
    Label: analytic, benchmark.
    """
    return _vec(lambda a: _mach_from_area_scalar(float(a), gamma, branch), _F(A_over_Astar))


def isentropic_state(M, gamma: float = 1.4) -> dict:
    """All isentropic ratios at one Mach number (one row of our version of the isentropic table).

    Book: §15.4, Eqs. (15.28)–(15.31); Table 15.1 (never reproduced — generated here).
    Parameters: M [–]; gamma [–].
    Returns dict: ``M``, ``T_over_T0``, ``p_over_p0``, ``rho_over_rho0``, ``A_over_Astar``, ``T0_over_T``,
    ``p0_over_p``, ``rho0_over_rho`` [–].
    Assumptions: perfect gas, isentropic.
    Validation: V1, V5 — tests/test_ch15.py: test_stagnation_ratios_V1_from_the_energy_law_and_the_isentrope,
    test_stagnation_ratios_V5_naca_report_1135_rows.
    Label: analytic, benchmark.
    """
    t = _F(stagnation_temperature_ratio(M, gamma))
    return dict(M=_S(_F(M)), T_over_T0=_S(1.0 / t), p_over_p0=_S(t ** (-gamma / (gamma - 1.0))),
                rho_over_rho0=_S(t ** (-1.0 / (gamma - 1.0))), A_over_Astar=area_ratio(M, gamma), T0_over_T=_S(t),
                p0_over_p=_S(t ** (gamma / (gamma - 1.0))), rho0_over_rho=_S(t ** (1.0 / (gamma - 1.0))))


def isentropic_table(M, gamma: float = 1.4):
    """Our isentropic-flow table: columns M, p/p0, rho/rho0, T/T0, A/A*.

    Book: §15.4, Eqs. (15.28)–(15.31) (the layout of Table 15.1; the values are computed, not transcribed).
    Parameters: M array [–]; gamma [–].   Returns a ``pandas.DataFrame``.
    Assumptions: as :func:`isentropic_state`.
    Validation: V1 — tests/test_ch15.py: test_stagnation_ratios_V1_from_the_energy_law_and_the_isentrope.
    Label: analytic.
    """
    import pandas as pd

    M = np.atleast_1d(_F(M))
    s = isentropic_state(M, gamma)
    return pd.DataFrame({"M": M, "p/p0": s["p_over_p0"], "rho/rho0": s["rho_over_rho0"], "T/T0": s["T_over_T0"],
                         "A/A*": s["A_over_Astar"]})


def sonic_ratios(gamma: float = 1.4) -> dict:
    """Sonic (critical) state over stagnation state: T*/T₀ = 2/(γ + 1), p*/p₀, ρ*/ρ₀, c*/c₀ (ours).

    Book: §15.4 (sonic conditions defined in words; the formulas are (15.28)–(15.30) at M = 1).
    Parameters: gamma [–].   Returns dict ``T``, ``p``, ``rho``, ``c`` [–] (0.8333, 0.5283, 0.6339, 0.9129 for 1.4).
    Assumptions: perfect gas, isentropic.
    Validation: V1 — tests/test_ch15.py: test_area_ratio_V1_built_from_mass_conservation_and_its_two_inverses.
    Label: analytic.
    """
    t = 2.0 / (gamma + 1.0)
    return dict(T=t, p=t ** (gamma / (gamma - 1.0)), rho=t ** (1.0 / (gamma - 1.0)), c=float(np.sqrt(t)))


def critical_pressure_ratio(gamma: float = 1.4) -> float:
    """p*/p₀ = [2/(γ + 1)]^{γ/(γ−1)} — the back-pressure ratio at and below which a convergent nozzle is choked (ours).

    Book: §15.7 (choking described in words), from Eq. (15.29) at M = 1.
    Parameters: gamma [–].   Returns p*/p₀ [–] (0.5283 for γ = 1.4).
    Assumptions: perfect gas, isentropic flow up to the exit.
    Validation: V1, V5 — tests/test_ch15.py: test_area_ratio_V1_built_from_mass_conservation_and_its_two_inverses,
    test_choked_flow_V5_cited_critical_ratio_and_gas_table.
    Label: analytic, benchmark.
    """
    return float((2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0)))


def max_speed(T0, cp: float):
    """Largest speed a gas can reach from a reservoir at T₀: u_max = √(2C_pT₀) (T → 0, M → ∞) (ours).

    Book: §15.4, Eq. (15.27) with T = 0 (Exercise 15.6).
    Parameters: T0 stagnation temperature [K]; cp [J/(kg K)].   Returns u_max [m/s].
    Assumptions: adiabatic steady flow of a perfect gas.
    Validation: V1 — tests/test_ch15.py: test_stagnation_ratios_V1_from_the_energy_law_and_the_isentrope.
    Label: analytic.
    """
    return _S(np.sqrt(2.0 * cp * _F(T0)))


def pitot_compressibility_factor(M, gamma: float = 1.4):
    """(p₀ − p)/(½ρu²) — how much a pitot reading exceeds the incompressible p + ½ρu² (exact; 1 + M²/4 + … ).

    Book: §15.4, Eq. (15.29) with ρu² = γpM² (Exercise 15.3).
    Parameters: M [–] ≥ 0; gamma [–].   Returns the factor [–] (1 at M = 0; 1.0227 at M = 0.3 for γ = 1.4).
    Method: ``expm1``/``log1p`` so that small M is free of cancellation.
    Assumptions: perfect gas, isentropic stagnation (subsonic pitot; no bow shock).
    Validation: V2 — tests/test_ch15.py: test_pitot_compressibility_V2_series_and_incompressible_limit.
    Label: symbolic.
    """
    M = _F(M)
    x = 0.5 * (gamma - 1.0) * M ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        fac = np.expm1(gamma / (gamma - 1.0) * np.log1p(x)) / (0.5 * gamma * M ** 2)
    return _S(np.where(M == 0.0, 1.0, fac))


def mach_for_density_change(fraction, gamma: float = 1.4):
    """Mach number at which the stagnation density exceeds the static density by a given fraction (ours).

    Book: §15.4, Eq. (15.30) inverted: ρ₀/ρ = 1 + fraction ⇒ M² = [2/(γ − 1)] [(1 + fraction)^{γ−1} − 1].
    Parameters: fraction (ρ₀ − ρ)/ρ [–] ≥ 0; gamma [–].   Returns M [–] (≈ √(2·fraction) for small values).
    Assumptions: isentropic.
    Validation: V2 — tests/test_ch15.py: test_pitot_compressibility_V2_series_and_incompressible_limit.
    Label: symbolic.
    """
    return _S(np.sqrt(2.0 / (gamma - 1.0) * np.expm1((gamma - 1.0) * np.log1p(_F(fraction)))))


def entropy_change_from_stagnation(p02_over_p01, T02_over_T01=1.0, gas: PerfectGas | None = None):
    """Entropy change between two sections from their stagnation states: S₂ − S₁ = −R ln(p₀₂/p₀₁) + C_p ln(T₀₂/T₀₁).

    Book: §15.4 (unnumbered display on the page of (15.27)–(15.30); from the entropy line of (15.1)).
    Parameters: p02_over_p01 [–]; T02_over_T01 [–] (1 for adiabatic flow); gas.
    Returns S₂ − S₁ [J/(kg K)].
    Assumptions: perfect gas; the stagnation state has the entropy of the static state.
    Validation: V1 — tests/test_ch15.py: test_entropy_from_stagnation_states_V1_matches_the_state_function.
    Label: analytic.
    """
    g = _gas(gas)
    return _S(-g.R * np.log(_F(p02_over_p01)) + g.cp * np.log(_F(T02_over_T01)))


def stagnation_pressure_loss(dS, gas: PerfectGas | None = None):
    """Stagnation-pressure ratio that goes with an entropy rise in adiabatic flow: p₀₂/p₀₁ = exp(−ΔS/R).

    Book: §15.4 (S₂ − S₁ = −R ln(p₀₂/p₀₁) with T₀₂ = T₀₁, hence p₀₂ < p₀₁).
    Parameters: dS entropy rise S₂ − S₁ [J/(kg K)]; gas.   Returns p₀₂/p₀₁ [–] (< 1 for dS > 0).
    Assumptions: adiabatic (T₀ constant), perfect gas.
    Validation: V1 — tests/test_ch15.py: test_entropy_from_stagnation_states_V1_matches_the_state_function.
    Label: analytic.
    """
    return _S(np.exp(-_F(dS) / _gas(gas).R))


def choked_mass_flow(p0, T0, A_star, gas: PerfectGas | None = None):
    """Largest mass flow a throat can pass: ṁ_max = p₀A* √(γ/(RT₀)) [2/(γ + 1)]^{(γ+1)/(2(γ−1))} (ours).

    Book: §15.4/§15.7 (choking; Exercise 15.15) — (15.18) at the sonic throat with (15.28)–(15.30) at M = 1.
    Parameters: p0 [Pa]; T0 [K]; A_star throat area [m²]; gas.   Returns ṁ_max [kg/s].
    Note (trap T10): "choked" fixes ṁ against the *back pressure* only; raising p₀ or lowering T₀ still changes it.
    Assumptions: isentropic flow from the reservoir to a sonic throat, perfect gas.
    Validation: V1, V2, V5 — tests/test_ch15.py: test_convergent_nozzle_V1_piecewise_mass_flow_against_back_pressure,
    test_sonic_ratios_and_choked_flow_V2_derivation, test_choked_flow_V5_cited_critical_ratio_and_gas_table.
    Label: analytic, symbolic, benchmark.
    """
    return mass_flow_from_mach(1.0, p0, T0, A_star, gas)


# ======================================================================================================================
# §15.5 — a duct of varying area
# ======================================================================================================================
def laval_area(x, A_throat: float = 1.0, A_inlet: float = 2.5, A_exit: float = 2.0, x_throat: float = 0.4,
               L: float = 1.0):
    """A smooth convergent–divergent test shape (ours): two parabolas meeting with zero slope at the throat.

    A(x) = A_t + (A_in − A_t)[(x_t − x)/x_t]² for x ≤ x_t,  A_t + (A_e − A_t)[(x − x_t)/(L − x_t)]² for x ≥ x_t.
    Book: §15.5, Figs. 15.5–15.6 (shape only; the numbers are ours).
    Parameters: x axial coordinate [m], 0 ≤ x ≤ L; A_throat, A_inlet, A_exit areas [m²]; x_throat, L [m].
    Returns A(x) [m²].
    Assumptions: none (geometry).
    Validation: V1, V7 — A(0), A(x_t), A(L) and zero slope at the throat (in the test of the duct laws) —
    tests/test_ch15.py: test_laval_nozzle_V4_consistent_with_shock_and_isentropic_relations,
    test_small_helpers_V1_against_independent_closed_forms.
    Label: analytic.
    """
    x = _F(x)
    up = A_throat + (A_inlet - A_throat) * ((x_throat - x) / x_throat) ** 2
    dn = A_throat + (A_exit - A_throat) * ((x - x_throat) / (L - x_throat)) ** 2
    return _S(np.where(x <= x_throat, up, dn))


def nozzle_isentropic(x, A_fn, gamma: float = 1.4, *, branch: str = "transonic", A_star: float | None = None) -> dict:
    """Isentropic solution along a duct of given area: M(x) and the ratios to the stagnation state.

    Book: §15.5 (Figs. 15.5–15.6) with Eq. (15.31) inverted point by point.
    Parameters
    ----------
    x : stations [m] (1-D array).   A_fn : callable A(x) [m²] or an array of areas at ``x``.   gamma : [–].
    branch : "subsonic" (subsonic everywhere: a venturi — the throat need not be sonic), "supersonic" (supersonic
        everywhere) or "transonic" (subsonic upstream of the minimum area, sonic there, supersonic downstream).
    A_star : sonic reference area [m²]; default min A (sonic throat).  Must be ≤ min A; with ``branch="subsonic"``
        and A_star < min A the throat Mach number is below 1.

    Returns
    -------
    dict: ``x``, ``A``, ``M``, ``p_over_p0``, ``T_over_T0``, ``rho_over_rho0``, ``A_star``, ``i_throat``, ``x_throat``.
    **The throat must be a grid point**: the sonic area and the switch of branch are taken from the smallest sampled
    area, so a grid that steps over the true minimum puts the throat (and A*) slightly wrong.
    Assumptions: steady, frictionless, adiabatic quasi-one-dimensional flow; no shock (see :func:`laval_nozzle`).
    Validation: V4 — tests/test_ch15.py: test_duct_laws_V4_mass_energy_entropy_along_an_isentropic_nozzle.
    Label: conserved.
    """
    x = np.atleast_1d(_F(x))
    A = np.atleast_1d(_F(A_fn(x) if callable(A_fn) else A_fn))
    it = int(np.argmin(A))
    As = float(A[it]) if A_star is None else float(A_star)
    if As > A[it] * (1.0 + 1e-12):
        raise ValueError("A_star exceeds the smallest area: the flow would be choked upstream of this duct")
    ar = np.maximum(A / As, 1.0)
    if branch == "subsonic":
        M = np.atleast_1d(mach_from_area_ratio(ar, gamma, branch="subsonic"))
    elif branch == "supersonic":
        M = np.atleast_1d(mach_from_area_ratio(ar, gamma, branch="supersonic"))
    elif branch == "transonic":
        M = np.where(np.arange(A.size) <= it, np.atleast_1d(mach_from_area_ratio(ar, gamma, branch="subsonic")),
                     np.atleast_1d(mach_from_area_ratio(ar, gamma, branch="supersonic")))
    else:
        raise ValueError('branch must be "subsonic", "supersonic" or "transonic"')
    s = isentropic_state(M, gamma)
    return dict(x=x, A=A, M=M, p_over_p0=np.atleast_1d(s["p_over_p0"]), T_over_T0=np.atleast_1d(s["T_over_T0"]),
                rho_over_rho0=np.atleast_1d(s["rho_over_rho0"]), A_star=As, i_throat=it, x_throat=float(x[it]))


# ======================================================================================================================
# §15.6 — normal shock waves
# ======================================================================================================================
def normal_shock_pressure_ratio(M1, gamma: float = 1.4):
    """p₂/p₁ = 1 + [2γ/(γ + 1)](M₁² − 1) across a normal shock.

    Book: §15.6, Eq. (15.39).
    Parameters: M1 upstream Mach number [–]; gamma [–].   Returns p₂/p₁ [–].
    The bare formula: it does not refuse M₁ < 1 (the forbidden expansion root) — :func:`normal_shock` does.
    Assumptions: perfect gas; steady, adiabatic, no wall force, constant area.
    Validation: V2, V5, V7 — tests/test_ch15.py: test_normal_shock_pressure_jump_V2_derivation,
    test_normal_shock_V5_naca_1135_and_cited_table_rows, test_normal_shock_V7_limits_guard_and_hugoniot.
    Label: symbolic, benchmark.
    """
    return _S(1.0 + 2.0 * gamma / (gamma + 1.0) * (_F(M1) ** 2 - 1.0))  # Eq. (15.39)


def normal_shock_mach(M1, gamma: float = 1.4):
    """Downstream Mach number of a normal shock, M₂² = [(γ − 1)M₁² + 2]/[2γM₁² + 1 − γ].

    Book: §15.6, Eq. (15.40); Fig. 15.8.
    Parameters: M1 [–]; gamma [–].   Returns M₂ [–]; NaN where M₁² ≤ (γ − 1)/(2γ) (no real root).
    **Trap T6**: the map is its own inverse — feeding it a subsonic M₁ returns the supersonic partner without any
    error (the expansion shock, forbidden by the second law).  Use :func:`normal_shock` for a guarded call.
    Assumptions: as (15.39).
    Validation: V4, V2 — tests/test_ch15.py: test_oblique_shock_V4_jump_laws_in_two_dimensions,
    test_normal_shock_downstream_mach_V2_derivation.
    Label: symbolic, conserved.
    """
    M1 = _F(M1)
    num = (gamma - 1.0) * M1 ** 2 + 2.0
    den = 2.0 * gamma * M1 ** 2 + 1.0 - gamma
    with np.errstate(divide="ignore", invalid="ignore"):
        return _S(np.where(den > 0.0, np.sqrt(num / np.where(den > 0.0, den, np.nan)), np.nan))  # Eq. (15.40)


def normal_shock_density_ratio(M1, gamma: float = 1.4):
    """ρ₂/ρ₁ = u₁/u₂ = (γ + 1)M₁²/[(γ − 1)M₁² + 2] across a normal shock.

    Book: §15.6, Eq. (15.41).
    Parameters: M1 [–]; gamma [–].   Returns ρ₂/ρ₁ [–] (→ (γ + 1)/(γ − 1) as M₁ → ∞).
    Assumptions: as (15.39).
    Validation: V2, V5 — tests/test_ch15.py: test_normal_shock_density_and_temperature_V2_derivation,
    test_normal_shock_V5_naca_1135_and_cited_table_rows.
    Label: symbolic, benchmark.
    """
    M1 = _F(M1)
    return _S((gamma + 1.0) * M1 ** 2 / ((gamma - 1.0) * M1 ** 2 + 2.0))  # Eq. (15.41)


def normal_shock_temperature_ratio(M1, gamma: float = 1.4):
    """T₂/T₁ = 1 + [2(γ − 1)/(γ + 1)²] [(γM₁² + 1)/M₁²] (M₁² − 1) across a normal shock.

    Book: §15.6, Eq. (15.42).
    Parameters: M1 [–]; gamma [–].   Returns T₂/T₁ [–].
    Assumptions: as (15.39).
    Validation: V2 — tests/test_ch15.py: test_normal_shock_density_and_temperature_V2_derivation.
    Label: symbolic.
    """
    M1 = _F(M1)
    return _S(1.0 + 2.0 * (gamma - 1.0) / (gamma + 1.0) ** 2 * (gamma * M1 ** 2 + 1.0) / M1 ** 2
              * (M1 ** 2 - 1.0))  # Eq. (15.42)


def normal_shock_entropy(M1, gamma: float = 1.4):
    """Entropy jump of a normal shock, (S₂ − S₁)/C_v = ln[(p₂/p₁)(ρ₁/ρ₂)^γ].

    Book: §15.6, Eq. (15.43); Fig. 15.9.  Negative for M₁ < 1 (why expansion shocks do not occur in a perfect gas).
    Parameters: M1 [–]; gamma [–].   Returns (S₂ − S₁)/C_v [–] (NaN where p₂/p₁ ≤ 0).
    Method: with m = M₁² − 1, a = 2γ/(γ + 1), b = (γ − 1)/(γ + 1) the jump is ln(1 + am) − γ ln(1 + m) + γ ln(1 + bm);
    the first two orders in m cancel exactly, so for |m| < 0.02 the function sums the Taylor series from the cubic
    term on (terms (−1)^{k+1}(a^k − γ + γb^k)m^k/k, k = 3…16) instead of subtracting logarithms.  m is formed as
    (M₁ − 1)(M₁ + 1), not M₁² − 1, so that it keeps full relative precision as M₁ → 1.
    Assumptions: as (15.39).
    Validation: V1, V2 — entropy state function of the jump state; 40-digit arithmetic near M₁ = 1; log–log slope 2.997
    against M₁² − 1 — tests/test_ch15.py: test_shock_entropy_V1_equals_the_entropy_state_function_of_the_jump_state,
    test_weak_shock_entropy_V2_derivation; full relative precision for |M₁ − 1| down to 1e-12 against 50-digit
    arithmetic and the sympy series — test_shock_entropy_V1_full_relative_precision_next_to_mach_one.
    Label: analytic, symbolic.
    """
    M1 = _F(M1)
    m = (M1 - 1.0) * (M1 + 1.0)  # M₁² − 1 without the cancellation near M₁ = 1
    a, b = 2.0 * gamma / (gamma + 1.0), (gamma - 1.0) / (gamma + 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        big = np.log1p(a * m) - gamma * np.log1p(m) + gamma * np.log1p(b * m)  # Eq. (15.43)
    small = np.zeros_like(m)
    for k in range(16, 2, -1):
        small = small + (-1.0) ** (k + 1) * (a ** k - gamma + gamma * b ** k) / k * m ** k
    return _S(np.where(np.abs(m) < 0.02, small, big))


def weak_shock_entropy(x, gamma: float = 1.4, *, form: str = "mach"):
    """Leading (cubic) term of the entropy jump of a weak shock.

    Book: §15.6, Eq. (15.44a)  (S₂ − S₁)/C_v ≅ [2γ(γ − 1)/(3(γ + 1)²)] (M₁² − 1)³  (``form="mach"``, x = M₁), and
    Eq. (15.44b)  (S₂ − S₁)/C_v ≅ [(γ² − 1)/(12γ²)] ((p₂ − p₁)/p₁)³  (``form="pressure"``, x = (p₂ − p₁)/p₁).
    Parameters: x — M₁ [–] or the pressure jump [–] according to ``form``; gamma [–].
    Returns the approximate (S₂ − S₁)/C_v [–].
    Assumptions: |M₁² − 1| ≪ 1.
    Validation: V1, V2 — tests/test_ch15.py: test_weak_shock_entropy_V1_cubic_law_and_its_first_correction,
    test_weak_shock_entropy_V2_derivation.
    Label: analytic, symbolic.
    """
    x = _F(x)
    if form == "mach":
        return _S(2.0 * gamma * (gamma - 1.0) / (3.0 * (gamma + 1.0) ** 2) * (x ** 2 - 1.0) ** 3)  # Eq. (15.44a)
    if form == "pressure":
        return _S((gamma ** 2 - 1.0) / (12.0 * gamma ** 2) * x ** 3)  # Eq. (15.44b)
    raise ValueError('form must be "mach" or "pressure"')


def shock_admissible(M1):
    """Second-law test for a normal shock in a perfect gas: S₂ − S₁ ≥ 0 ⇔ M₁ ≥ 1.

    Book: §15.6 (text after (15.44a); Fig. 15.9).
    Parameters: M1 [–].   Returns bool (array): True for M₁ ≥ 1.
    **At M₁ = 1 exactly** S₂ − S₁ = 0 (no jump: the "shock" is a sound wave), which the second law (ΔS ≥ 0)
    allows, so this returns True there — the same rule as the ``admissible`` flag of :func:`normal_shock`.
    Assumptions: perfect gas with constant γ (dense gases near the critical point can behave differently).
    Validation: V1 — tests/test_ch15.py: test_shock_entropy_V1_equals_the_entropy_state_function_of_the_jump_state.
    Label: analytic.
    """
    r = _F(M1) >= 1.0
    return bool(r) if r.ndim == 0 else r


def normal_shock_stagnation_pressure_ratio(M1, gamma: float = 1.4):
    """Stagnation-pressure ratio across a normal shock, p₀₂/p₀₁ = exp(−(S₂ − S₁)/R) = A₁*/A₂* (ours).

    Book: §15.6 (stated: T₀ and h₀ unchanged, p₀ and ρ₀ fall; only the last column of Table 15.2).  From (15.43) and
    the stagnation-entropy relation of §15.4:
    p₀₂/p₀₁ = [(γ + 1)M₁²/((γ − 1)M₁² + 2)]^{γ/(γ−1)} [(γ + 1)/(2γM₁² − (γ − 1))]^{1/(γ−1)}.
    Parameters: M1 [–]; gamma [–].   Returns p₀₂/p₀₁ [–] (≤ 1 for M₁ ≥ 1; > 1 for the forbidden M₁ < 1).
    Method: exp(−(S₂ − S₁)/(C_v(γ − 1))) with :func:`normal_shock_entropy` (accurate near M₁ = 1).
    Assumptions: as (15.39).
    Validation: V1, V2, V5 — tests/test_ch15.py: test_weak_shock_entropy_V1_cubic_law_and_its_first_correction,
    test_shock_stagnation_pressure_ratio_V2_derivation, test_normal_shock_V5_naca_1135_and_cited_table_rows.
    Label: analytic, symbolic, benchmark.
    """
    return _S(np.exp(-_F(normal_shock_entropy(M1, gamma)) / (gamma - 1.0)))


def normal_shock(M1, gamma: float = 1.4, *, allow_expansion: bool = False) -> dict:
    """Everything across a normal shock for a given **upstream** Mach number.

    Book: §15.6, Eqs. (15.39)–(15.43) (Rankine–Hugoniot relations) and the second-law selection M₁ > 1.
    Parameters
    ----------
    M1 : upstream Mach number [–] (station 1 is upstream).   gamma : [–].
    allow_expansion : if False (default) any M₁ < 1 raises ``ValueError`` — because (15.40) is its own inverse, a
        downstream value passed by mistake would silently return an expansion shock (trap T6).  If True the
        conservation-law root is returned for M₁ < 1 as well, flagged ``admissible=False``.

    Returns
    -------
    dict [–]: ``M1``, ``M2``, ``p_ratio`` (p₂/p₁), ``rho_ratio`` (ρ₂/ρ₁ = u₁/u₂), ``u_ratio`` (u₂/u₁), ``T_ratio``,
    ``p0_ratio`` (p₀₂/p₀₁), ``T0_ratio`` (= 1), ``dS_over_cv``, ``dS_over_R``, ``admissible`` (M₁ ≥ 1: ΔS ≥ 0;
    True at M₁ = 1, where ΔS = 0 — the same rule as :func:`shock_admissible`).
    Assumptions: perfect gas, constant γ; the shock is a discontinuity (thin compared with the duct).
    Validation: V4, V2, V5, V7 — tests/test_ch15.py:
    test_normal_shock_V4_the_three_jump_laws_hold_for_the_returned_state,
    test_shock_structure_V2_first_integrals_and_roots, test_normal_shock_V5_naca_1135_and_cited_table_rows,
    test_normal_shock_V7_limits_guard_and_hugoniot.
    Label: symbolic, conserved, benchmark.
    """
    M1a = _F(M1)
    if not allow_expansion and np.any(M1a < 1.0):
        raise ValueError("normal_shock needs the UPSTREAM Mach number M1 >= 1 (station 1 is upstream; a shock from "
                         "M1 < 1 would lower the entropy). Pass allow_expansion=True to see the forbidden root.")
    ds = _F(normal_shock_entropy(M1a, gamma))
    rr = _F(normal_shock_density_ratio(M1a, gamma))
    adm = M1a >= 1.0
    return dict(M1=_S(M1a), M2=normal_shock_mach(M1a, gamma), p_ratio=normal_shock_pressure_ratio(M1a, gamma),
                rho_ratio=_S(rr), u_ratio=_S(1.0 / rr), T_ratio=normal_shock_temperature_ratio(M1a, gamma),
                p0_ratio=_S(np.exp(-ds / (gamma - 1.0))), T0_ratio=_S(np.ones_like(M1a)), dS_over_cv=_S(ds),
                dS_over_R=_S(ds / (gamma - 1.0)), admissible=bool(adm) if adm.ndim == 0 else adm)


def normal_shock_residuals(state1: dict, state2: dict, gamma: float = 1.4) -> dict:
    """Relative residuals of the three jump laws for a pair of states (zero for a true shock pair).

    Book: §15.6, Eqs. (15.35) ρ₁u₁ = ρ₂u₂, (15.36) p₁ − p₂ = −ρ₁u₁² + ρ₂u₂², (15.37)/(15.38)
    [γ/(γ − 1)]p₁/ρ₁ + ½u₁² = [γ/(γ − 1)]p₂/ρ₂ + ½u₂².
    Parameters: state1 (upstream), state2 (downstream) dicts (``rho``, ``u``, ``p``); gamma [–].
    Returns dict ``mass``, ``momentum``, ``energy`` [–], each (side 2 − side 1)/side 1.
    Assumptions: perfect gas; constant area, no wall force, no heat.
    Validation: V4 — tests/test_ch15.py: test_normal_shock_V4_the_three_jump_laws_hold_for_the_returned_state.
    Label: conserved.
    """
    r1, u1, p1 = _F(state1["rho"]), _F(state1["u"]), _F(state1["p"])
    r2, u2, p2 = _F(state2["rho"]), _F(state2["u"]), _F(state2["p"])
    m1, m2 = r1 * u1, r2 * u2  # Eq. (15.35)
    j1, j2 = p1 + r1 * u1 ** 2, p2 + r2 * u2 ** 2  # Eq. (15.36)
    e1 = gamma / (gamma - 1.0) * p1 / r1 + 0.5 * u1 ** 2  # Eq. (15.38)
    e2 = gamma / (gamma - 1.0) * p2 / r2 + 0.5 * u2 ** 2
    return dict(mass=_S((m2 - m1) / m1), momentum=_S((j2 - j1) / j1), energy=_S((e2 - e1) / e1))


def normal_shock_fluxes(M1, gamma: float = 1.4) -> dict:
    """Mass, momentum and energy fluxes on both sides of a normal shock, in units of the upstream values, split into
    their parts (the paired bars of the shock explainer).

    Book: §15.6, Eqs. (15.35)–(15.37) evaluated with (15.39)–(15.42).
    Parameters: M1 [–] (any M₁ with a real M₂; the forbidden side M₁ < 1 is allowed so that it can be shown);
    gamma [–].
    Returns dict [–]: ``mass`` (ρ₂u₂/ρ₁u₁), ``momentum`` ((p₂ + ρ₂u₂²)/(p₁ + ρ₁u₁²)), ``energy``
    ((h₂ + ½u₂²)/(h₁ + ½u₁²)) — all 1 — and the parts as fractions of the upstream total: ``p1_part``, ``ke1_part``
    (pressure and ρu² shares of the momentum flux upstream), ``p2_part``, ``ke2_part``; ``h1_part``, ``k1_part``
    (enthalpy and kinetic shares of the total enthalpy upstream), ``h2_part``, ``k2_part``.
    Assumptions: as (15.39).
    Validation: V4 — tests/test_ch15.py: test_normal_shock_V4_the_three_jump_laws_hold_for_the_returned_state.
    Label: conserved.
    """
    M1 = _F(M1)
    M2 = _F(normal_shock_mach(M1, gamma))
    pr, rr, Tr = (_F(normal_shock_pressure_ratio(M1, gamma)), _F(normal_shock_density_ratio(M1, gamma)),
                  _F(normal_shock_temperature_ratio(M1, gamma)))
    j1 = 1.0 + gamma * M1 ** 2
    e1 = 1.0 + 0.5 * (gamma - 1.0) * M1 ** 2
    p2, k2 = pr / j1, pr * gamma * M2 ** 2 / j1
    h2, kk2 = Tr / e1, Tr * 0.5 * (gamma - 1.0) * M2 ** 2 / e1
    out = dict(mass=rr * (M2 / M1) * np.sqrt(Tr), momentum=p2 + k2, energy=h2 + kk2, p1_part=1.0 / j1,
               ke1_part=gamma * M1 ** 2 / j1, p2_part=p2, ke2_part=k2, h1_part=1.0 / e1,
               k1_part=0.5 * (gamma - 1.0) * M1 ** 2 / e1, h2_part=h2, k2_part=kk2)
    return {k: _S(v) for k, v in out.items()}


def normal_shock_table(M1, gamma: float = 1.4):
    """Our normal-shock table: columns M1, M2, p2/p1, rho2/rho1, T2/T1, p02/p01.

    Book: §15.6, Eqs. (15.39)–(15.42) and the stagnation-pressure ratio (the layout of Table 15.2; computed here).
    Parameters: M1 array [–] ≥ 1; gamma [–].   Returns a ``pandas.DataFrame``.
    Assumptions: as :func:`normal_shock`.
    Validation: V5 — tests/test_ch15.py: test_normal_shock_V5_naca_1135_and_cited_table_rows.
    Label: benchmark.
    """
    import pandas as pd

    M1 = np.atleast_1d(_F(M1))
    s = normal_shock(M1, gamma)
    return pd.DataFrame({"M1": M1, "M2": s["M2"], "p2/p1": s["p_ratio"], "rho2/rho1": s["rho_ratio"],
                         "T2/T1": s["T_ratio"], "p02/p01": s["p0_ratio"]})


def hugoniot_pressure_ratio(rho_ratio, gamma: float = 1.4):
    """The Hugoniot curve: p₂/p₁ = [(γ + 1)ρ₂/ρ₁ − (γ − 1)]/[(γ + 1) − (γ − 1)ρ₂/ρ₁] (ours).

    Book: §15.6 (Rankine–Hugoniot relations named) — M₁ eliminated between (15.39) and (15.41).  Compare with the
    isentrope p₂/p₁ = (ρ₂/ρ₁)^γ (15.2): the two curves touch to second order at ρ₂/ρ₁ = 1.
    Parameters: rho_ratio ρ₂/ρ₁ [–], below (γ + 1)/(γ − 1); gamma [–].   Returns p₂/p₁ [–] (NaN beyond the limit).
    Assumptions: as (15.39).
    Validation: V5, V7 — tests/test_ch15.py: test_normal_shock_V5_naca_1135_and_cited_table_rows,
    test_normal_shock_V7_limits_guard_and_hugoniot.
    Label: benchmark.
    """
    r = _F(rho_ratio)
    den = (gamma + 1.0) - (gamma - 1.0) * r
    with np.errstate(divide="ignore", invalid="ignore"):
        return _S(np.where(den > 0.0, ((gamma + 1.0) * r - (gamma - 1.0)) / den, np.nan))


def strong_shock_limits(gamma: float = 1.4) -> dict:
    """Limits of the normal-shock relations as M₁ → ∞ (ours; Exercise 15.9a).

    Book: §15.6, Eqs. (15.39)–(15.42) with M₁² ≫ 1; Fig. 15.8 (asymptote of M₂).
    Parameters: gamma [–].
    Returns dict [–]: ``rho_ratio`` = (γ + 1)/(γ − 1); ``M2`` = √((γ − 1)/(2γ)); ``p_ratio_over_M1sq`` = 2γ/(γ + 1);
    ``T_ratio_over_p_ratio`` = (γ − 1)/(γ + 1); ``u1sq_rho1_over_p2`` = (γ + 1)/2 (so u₁ → √(½(γ + 1)p₂/ρ₁)).
    Assumptions: perfect gas with constant γ (real gases dissociate first).
    Validation: V1, V7 — tests/test_ch15.py: test_small_helpers_V1_against_independent_closed_forms,
    test_normal_shock_V7_limits_guard_and_hugoniot.
    Label: analytic.
    """
    return dict(rho_ratio=(gamma + 1.0) / (gamma - 1.0), M2=float(np.sqrt((gamma - 1.0) / (2.0 * gamma))),
                p_ratio_over_M1sq=2.0 * gamma / (gamma + 1.0), T_ratio_over_p_ratio=(gamma - 1.0) / (gamma + 1.0),
                u1sq_rho1_over_p2=0.5 * (gamma + 1.0))


def moving_shock(p1, T1, *, M_s=None, p2_over_p1=None, gas: PerfectGas | None = None) -> dict:
    """A normal shock running into gas at rest (laboratory frame), by a Galilean shift of the stationary shock.

    Book: §15.6 "Moving normal shock wave in a stationary medium", Fig. 15.10: subtract u₁ everywhere — the shock moves
    at u₁ = M₁c₁ > c₁ and the gas behind follows it at u₁ − u₂.
    Parameters
    ----------
    p1, T1 : pressure [Pa] and temperature [K] of the still gas ahead.
    M_s : shock Mach number u_shock/c₁ [–] (≥ 1), **or** p2_over_p1 : the pressure ratio [–] (≥ 1) — give one.
    gas : :class:`PerfectGas` (default air).

    Returns
    -------
    dict: ``M_s`` [–], ``u_shock`` [m/s], ``u_gas`` [m/s] (speed of the gas behind, in the direction of the shock),
    ``p2`` [Pa], ``T2`` [K], ``rho1``, ``rho2`` [kg/m³], ``c1``, ``c2`` [m/s], ``M_gas`` = u_gas/c₂ [–],
    ``T0_behind_lab`` [K] and ``p0_behind_lab`` [Pa] (stagnation values of the moving gas in the laboratory frame —
    *higher* than T₁, p₁: stagnation quantities depend on the frame, trap T7).
    Assumptions: perfect gas; plane shock of constant speed.
    Validation: V4 — tests/test_ch15.py: test_moving_shock_V4_conservation_in_the_laboratory_frame.
    Label: conserved.
    """
    g = _gas(gas)
    if (M_s is None) == (p2_over_p1 is None):
        raise ValueError("give exactly one of M_s and p2_over_p1")
    if M_s is None:
        M_s = np.sqrt(1.0 + (g.gamma + 1.0) / (2.0 * g.gamma) * (_F(p2_over_p1) - 1.0))  # Eq. (15.39) inverted
    s = normal_shock(M_s, g.gamma)
    p1, T1 = _F(p1), _F(T1)
    c1 = np.sqrt(g.gamma * g.R * T1)
    u1 = _F(s["M1"]) * c1
    u_gas = u1 * (1.0 - _F(s["u_ratio"]))  # u₁ − u₂
    T2, p2 = T1 * _F(s["T_ratio"]), p1 * _F(s["p_ratio"])
    c2 = np.sqrt(g.gamma * g.R * T2)
    Mg = u_gas / c2
    out = dict(M_s=s["M1"], u_shock=u1, u_gas=u_gas, p2=p2, T2=T2, rho1=p1 / (g.R * T1), rho2=p2 / (g.R * T2), c1=c1,
               c2=c2, M_gas=Mg, T0_behind_lab=T2 * _F(stagnation_temperature_ratio(Mg, g.gamma)),
               p0_behind_lab=p2 * _F(stagnation_pressure_ratio(Mg, g.gamma)))
    return {k: _S(v) for k, v in out.items()}


# ======================================================================================================================
# §15.7 — nozzles at different back pressures (every formula of this section is ours: the book has none)
# ======================================================================================================================
def _one_gamma(gamma: float | None, gas: PerfectGas | None, where: str) -> float:
    """The single γ of a call that accepts both ``gamma`` and ``gas`` (``gamma=None`` = "1.4, or gas.gamma")."""
    if gas is None:
        return 1.4 if gamma is None else float(gamma)
    if gamma is not None and not np.isclose(float(gamma), gas.gamma, rtol=1e-12, atol=0.0):
        raise ValueError(f"{where}: gamma={gamma!r} disagrees with gas.gamma={gas.gamma!r} — one call, one gas "
                         "(omit gamma when gas is given)")
    return float(gas.gamma)


def mass_flow_vs_back_pressure(pB_over_p0, gamma: float | None = None, *, dimensionless: bool = True, p0=None,
                               T0=None, A_exit=None, gas: PerfectGas | None = None) -> dict:
    """Flow through a **convergent** nozzle as the back pressure is lowered: rises, then stays at the choked value.

    Book: §15.7, Fig. 15.12 (described in words).  Ours: for p_B/p₀ > p*/p₀ the exit is subsonic and p_exit = p_B, so
    M_exit follows from (15.29); for p_B/p₀ ≤ p*/p₀ = [2/(γ + 1)]^{γ/(γ−1)} the exit is sonic, p_exit = p* > p_B;
    dimensional ṁ = p₀A_exit/√(RT₀) · √γ M_e [1 + ½(γ − 1)M_e²]^{−(γ+1)/(2(γ−1))} with **one** γ throughout.
    Parameters: pB_over_p0 [–] in [0, 1]; gamma [–] — default ``None`` means 1.4, or ``gas.gamma`` when ``gas`` is
    given (an explicit ``gamma`` that disagrees with ``gas.gamma`` raises ``ValueError``: γ and R must come from the
    same gas); dimensionless — if True ``mdot`` is ṁ√(RT₀)/(p₀A_exit) [–]; if False give p0 [Pa], T0 [K], A_exit [m²]
    (and gas, default air) and ``mdot`` is in kg/s.
    Returns dict: ``mdot``, ``M_exit`` [–], ``p_exit_over_p0`` [–], ``choked`` (bool), ``mdot_over_choked`` [–].
    Assumptions: isentropic flow inside the nozzle; exit = smallest area; reservoir at rest.
    Validation: V1, V4 — tests/test_ch15.py: test_convergent_nozzle_V1_piecewise_mass_flow_against_back_pressure,
    test_convergent_nozzle_mass_flow_V1_one_gamma_per_call.
    Label: analytic, conserved.
    """
    gamma = _one_gamma(gamma, gas, "mass_flow_vs_back_pressure")
    r = _F(pB_over_p0)
    if np.any(r > 1.0) or np.any(r < 0.0):
        raise ValueError("pB_over_p0 must lie in [0, 1] (a back pressure above the supply pressure reverses the flow)")
    crit = critical_pressure_ratio(gamma)
    choked = r <= crit
    Me = np.where(choked, 1.0, _F(mach_from_pressure_ratio(np.where(choked, 1.0, r), gamma)))
    flux = _F(mass_flux_function(Me, gamma))
    mdot = flux
    if not dimensionless:
        if p0 is None or T0 is None or A_exit is None:
            raise ValueError("dimensionless=False needs p0, T0 and A_exit")
        mdot = flux * _F(p0) * _F(A_exit) / np.sqrt(_gas(gas).R * _F(T0))
    return dict(mdot=_S(mdot), M_exit=_S(Me), p_exit_over_p0=_S(np.where(choked, crit, r)),
                choked=bool(choked) if choked.ndim == 0 else choked,
                mdot_over_choked=_S(flux / float(mass_flux_function(1.0, gamma))))


def convergent_nozzle(p0: float, T0: float, pB: float, A_exit: float, gas: PerfectGas | None = None) -> dict:
    """Exit state and mass flow of a convergent nozzle fed from a reservoir (dimensional wrapper).

    Book: §15.7, Fig. 15.12.
    Parameters: p0 supply (stagnation) pressure [Pa]; T0 [K]; pB back pressure [Pa] (≤ p0); A_exit [m²]; gas.
    Returns dict: ``M_exit`` [–], ``p_exit`` [Pa], ``T_exit`` [K], ``u_exit`` [m/s], ``mdot`` [kg/s], ``mdot_max``
    [kg/s], ``choked`` (bool), ``matched`` (p_exit = p_B, true exactly when the exit is subsonic or just sonic).
    Assumptions: as :func:`mass_flow_vs_back_pressure`.
    Validation: V1, V4 — tests/test_ch15.py: test_convergent_nozzle_V1_piecewise_mass_flow_against_back_pressure.
    Label: analytic, conserved.
    """
    g = _gas(gas)
    s = mass_flow_vs_back_pressure(pB / p0, g.gamma, dimensionless=False, p0=p0, T0=T0, A_exit=A_exit, gas=g)
    Te = T0 / float(stagnation_temperature_ratio(s["M_exit"], g.gamma))
    pe = s["p_exit_over_p0"] * p0
    return dict(M_exit=s["M_exit"], p_exit=pe, T_exit=Te, u_exit=s["M_exit"] * float(g.c(Te)), mdot=s["mdot"],
                mdot_max=float(choked_mass_flow(p0, T0, A_exit, g)), choked=s["choked"],
                matched=bool(pB / p0 >= critical_pressure_ratio(g.gamma)))


def laval_critical_pressures(Ae_over_At, gamma: float = 1.4) -> dict:
    """The three back pressures that divide the regimes of a convergent–divergent nozzle (ours).

    Book: §15.7, Fig. 15.13 (the dividing curves are drawn; no formula is given).  With a sonic throat:
    p_b/p₀ = p/p₀ at M_sub(A_e/A_t) (subsonic everywhere, throat just sonic);
    p_f/p₀ = p/p₀ at M_sup(A_e/A_t) (the design, pressure-matched, supersonic exit);
    p_d/p₀ = (p_f/p₀) × p₂/p₁ at M_sup (a normal shock standing in the exit plane).  Always p_f < p_d < p_b < p₀.
    Parameters: Ae_over_At exit-to-throat area ratio [–] ≥ 1; gamma [–].
    Returns dict [–]: ``p_b``, ``p_d``, ``p_f``, ``M_sub``, ``M_sup``, ``p_star`` (throat pressure when choked).
    Assumptions: isentropic flow except across a normal shock; quasi-one-dimensional.
    Validation: V1, V2, V4 — tests/test_ch15.py: test_laval_nozzle_V1_three_thresholds_and_seven_regimes,
    test_laval_nozzle_shock_position_V2_derivation, test_laval_nozzle_V4_consistent_with_shock_and_isentropic_relations.
    Label: analytic, conserved, symbolic.
    """
    ar = _F(Ae_over_At)
    Msub = _F(mach_from_area_ratio(ar, gamma, branch="subsonic"))
    Msup = _F(mach_from_area_ratio(ar, gamma, branch="supersonic"))
    pb = 1.0 / _F(stagnation_pressure_ratio(Msub, gamma))
    pf = 1.0 / _F(stagnation_pressure_ratio(Msup, gamma))
    pd_ = pf * _F(normal_shock_pressure_ratio(Msup, gamma))
    return dict(p_b=_S(pb), p_d=_S(pd_), p_f=_S(pf), M_sub=_S(Msub), M_sup=_S(Msup),
                p_star=critical_pressure_ratio(gamma))


LAVAL_REGIMES = ("no flow", "subsonic", "sonic throat", "shock in nozzle", "shock at exit", "overexpanded", "design",
                 "underexpanded")


def laval_regime(pB_over_p0: float, Ae_over_At: float, gamma: float = 1.4, *, rtol: float = 1e-9) -> str:
    """Name the operating condition of a convergent–divergent nozzle for a given back pressure.

    Book: §15.7, Fig. 15.13 (seven conditions from fully subsonic to under-expanded).
    Parameters: pB_over_p0 [–] in [0, 1]; Ae_over_At [–] ≥ 1; gamma [–]; rtol — relative band within which a back
    pressure counts as *equal* to a dividing value (our numerical choice).
    Returns one of :data:`LAVAL_REGIMES`: "no flow" (p_B = p₀); "subsonic" (p_b < p_B < p₀: a venturi, not choked);
    "sonic throat" (p_B = p_b); "shock in nozzle" (p_d < p_B < p_b); "shock at exit" (p_B = p_d); "overexpanded"
    (p_f < p_B < p_d: oblique shocks outside); "design" (p_B = p_f); "underexpanded" (p_B < p_f: expansion fans
    outside).
    Assumptions: as :func:`laval_critical_pressures`.
    Validation: V1, V4 — tests/test_ch15.py: test_laval_nozzle_V1_three_thresholds_and_seven_regimes,
    test_laval_nozzle_V4_consistent_with_shock_and_isentropic_relations.
    Label: analytic, conserved.
    """
    r = float(pB_over_p0)
    if r > 1.0 or r < 0.0:
        raise ValueError("pB_over_p0 must lie in [0, 1]")
    c = laval_critical_pressures(Ae_over_At, gamma)
    eq = lambda a, b: abs(a - b) <= rtol * b  # noqa: E731
    if eq(r, 1.0):
        return "no flow"
    if eq(r, c["p_b"]):
        return "sonic throat"
    if r > c["p_b"]:
        return "subsonic"
    if eq(r, c["p_d"]):
        return "shock at exit"
    if r > c["p_d"]:
        return "shock in nozzle"
    if eq(r, c["p_f"]):
        return "design"
    return "overexpanded" if r > c["p_f"] else "underexpanded"


def laval_shock_location(pB_over_p0: float, Ae_over_At: float, gamma: float = 1.4) -> dict:
    """Strength and position (as an area) of the normal shock inside the diverging part for a given back pressure
    (ours).

    Book: §15.7, Fig. 15.13 curves between "sonic throat" and "shock at exit"; the method of Example 15.2 run
    backwards.  Find M₁ with  p_B/p₀ = (p₀₂/p₀₁)(M₁) × p/p₀ at M_sub((A_e/A_t)(p₀₂/p₀₁)) — the sonic area grows across
    the shock, A₂* = A_t p₀₁/p₀₂ (trap T8) — then A_shock/A_t = A/A* at M₁.
    Parameters: pB_over_p0 [–] with p_d ≤ p_B/p₀ ≤ p_b; Ae_over_At [–] > 1; gamma [–].
    Returns dict [–]: ``M1`` (Mach number just ahead of the shock), ``M2``, ``A_shock_over_At``, ``p02_over_p01``,
    ``M_exit`` (subsonic), ``p_exit_over_p0`` (= p_B/p₀).
    Method: ``brentq`` on M₁ ∈ [1, M_sup(A_e/A_t)]; residual asserted.  Raises ``ValueError`` outside the range.
    Assumptions: isentropic on both sides of a normal shock; the exit pressure equals the back pressure (subsonic exit).
    Validation: V2, V4 — tests/test_ch15.py: test_laval_nozzle_shock_position_V2_derivation,
    test_laval_nozzle_V4_consistent_with_shock_and_isentropic_relations.
    Label: conserved, symbolic.
    """
    r, ar = float(pB_over_p0), float(Ae_over_At)
    c = laval_critical_pressures(ar, gamma)
    if r > c["p_b"] * (1.0 + 1e-12) or r < c["p_d"] * (1.0 - 1e-12):
        raise ValueError(f"no normal shock inside the nozzle: p_B/p0 = {r:.6g} is outside [p_d, p_b] = "
                         f"[{c['p_d']:.6g}, {c['p_b']:.6g}]")

    def exit_p(M1: float) -> tuple[float, float, float]:
        k = float(normal_shock_stagnation_pressure_ratio(M1, gamma))
        Me = float(mach_from_area_ratio(max(ar * k, 1.0), gamma, branch="subsonic"))
        return k / float(stagnation_pressure_ratio(Me, gamma)), k, Me

    g = lambda M1: exit_p(M1)[0] - r  # noqa: E731
    if g(1.0) <= 0.0:
        M1 = 1.0
    elif g(c["M_sup"]) >= 0.0:
        M1 = float(c["M_sup"])
    else:
        M1 = float(brentq(g, 1.0, c["M_sup"], xtol=1e-14, rtol=1e-14))
    pe, k, Me = exit_p(M1)
    assert abs(pe - r) < 1e-8, "shock-location root did not converge"
    return dict(M1=M1, M2=float(normal_shock_mach(M1, gamma)), A_shock_over_At=float(area_ratio(M1, gamma)),
                p02_over_p01=k, M_exit=Me, p_exit_over_p0=pe)


def laval_exit_state(pB_over_p0: float, Ae_over_At: float, gamma: float = 1.4) -> dict:
    """Everything the back-pressure explainer reads out, without a nozzle shape (ours).

    Book: §15.7, Fig. 15.13.
    Parameters: pB_over_p0 [–] in (0, 1]; Ae_over_At [–] > 1; gamma [–].
    Returns dict: ``regime`` (see :func:`laval_regime`); ``M_exit`` [–] (exit-plane Mach number); ``p_exit_over_p0``
    [–] (exit-*plane* pressure: equals p_B/p₀ only when the exit is subsonic; p_f/p₀ when supersonic);
    ``M1_shock``, ``A_shock_over_At``, ``p02_over_p01`` (NaN, NaN, 1 when there is no shock inside; a shock in the
    exit plane is reported with its upstream state as the exit state and ``p02_over_p01`` = 1 inside the nozzle);
    ``M_throat`` [–]; ``mdot_over_choked`` [–] (1 whenever the throat is sonic); ``choked`` (bool).
    Assumptions: as :func:`laval_shock_location`; outside the exit the jet adjusts through oblique shocks or
    expansion fans that are not computed here.
    Validation: V4 — tests/test_ch15.py: test_laval_nozzle_V4_consistent_with_shock_and_isentropic_relations.
    Label: conserved.
    """
    r, ar = float(pB_over_p0), float(Ae_over_At)
    reg = laval_regime(r, ar, gamma)
    c = laval_critical_pressures(ar, gamma)
    out = dict(regime=reg, M1_shock=np.nan, A_shock_over_At=np.nan, p02_over_p01=1.0, M_throat=1.0,
               mdot_over_choked=1.0, choked=True)
    if reg in ("no flow", "subsonic"):
        Me = float(mach_from_pressure_ratio(r, gamma))
        if Me == 0.0:
            out.update(M_exit=0.0, p_exit_over_p0=1.0, M_throat=0.0, mdot_over_choked=0.0, choked=False)
            return out
        At_over_Astar = float(area_ratio(Me, gamma)) / ar
        Mt = float(mach_from_area_ratio(max(At_over_Astar, 1.0), gamma, branch="subsonic"))
        out.update(M_exit=Me, p_exit_over_p0=r, M_throat=Mt, mdot_over_choked=1.0 / max(At_over_Astar, 1.0),
                   choked=False)
    elif reg == "sonic throat":
        out.update(M_exit=float(c["M_sub"]), p_exit_over_p0=float(c["p_b"]))
    elif reg == "shock in nozzle":
        s = laval_shock_location(r, ar, gamma)
        out.update(M_exit=s["M_exit"], p_exit_over_p0=r, M1_shock=s["M1"], A_shock_over_At=s["A_shock_over_At"],
                   p02_over_p01=s["p02_over_p01"])
    else:  # shock at exit, overexpanded, design, underexpanded: supersonic up to the exit plane
        out.update(M_exit=float(c["M_sup"]), p_exit_over_p0=float(c["p_f"]))
        if reg == "shock at exit":
            out.update(M1_shock=float(c["M_sup"]), A_shock_over_At=ar)
    return out


def laval_nozzle(x, A_fn, pB_over_p0: float, gamma: float = 1.4) -> dict:
    """Mach number and pressure along a convergent–divergent nozzle for a given back pressure, shock included.

    Book: §15.7, Fig. 15.13.
    Parameters: x stations [m] (increasing, inlet to exit); A_fn callable A(x) [m²] or array of areas at ``x``;
    pB_over_p0 [–] in (0, 1]; gamma [–].
    Returns dict: ``x``, ``A``, ``M``, ``p_over_p0`` (p(x)/p₀₁, with the upstream stagnation pressure), ``regime``,
    ``x_shock`` [m] (NaN without a shock inside), ``A_shock_over_At``, ``M1_shock``, ``p02_over_p01``, ``M_exit``,
    ``p_exit_over_p0`` (exit plane), ``i_throat``, ``x_throat``, ``mdot_over_choked``.
    **The throat must be a grid point** (the throat area is the smallest sampled area).
    Stations with A(x) < A_shock downstream of the throat are supersonic (ahead of the shock), the rest subsonic with
    A₂* = A_t p₀₁/p₀₂ and p/p₀₁ = (p₀₂/p₀₁) p/p₀₂; ``x_shock`` is found by ``brentq`` on A(x) = A_shock when ``A_fn``
    is callable and by linear interpolation otherwise.
    Assumptions: as :func:`laval_shock_location`; the diverging part is monotonic.
    Validation: V4 — tests/test_ch15.py: test_laval_nozzle_V4_consistent_with_shock_and_isentropic_relations.
    Label: conserved.
    """
    x = np.atleast_1d(_F(x))
    A = np.atleast_1d(_F(A_fn(x) if callable(A_fn) else A_fn))
    it = int(np.argmin(A))
    At, ar = float(A[it]), float(A[-1] / A[it])
    st = laval_exit_state(pB_over_p0, ar, gamma)
    reg = st["regime"]
    idx = np.arange(A.size)
    x_shock = np.nan
    if reg in ("no flow", "subsonic"):
        if st["M_exit"] == 0.0:
            M = np.zeros_like(A)
        else:
            Astar = A[-1] / float(area_ratio(st["M_exit"], gamma))
            M = np.atleast_1d(mach_from_area_ratio(np.maximum(A / Astar, 1.0), gamma, branch="subsonic"))
        p = 1.0 / np.atleast_1d(_F(stagnation_pressure_ratio(M, gamma)))
    else:
        Msub = np.atleast_1d(mach_from_area_ratio(np.maximum(A / At, 1.0), gamma, branch="subsonic"))
        Msup = np.atleast_1d(mach_from_area_ratio(np.maximum(A / At, 1.0), gamma, branch="supersonic"))
        if reg == "sonic throat":
            M = Msub
            p = 1.0 / np.atleast_1d(_F(stagnation_pressure_ratio(M, gamma)))
        elif reg == "shock in nozzle":
            k, Ash = st["p02_over_p01"], st["A_shock_over_At"] * At
            if callable(A_fn) and it < A.size - 1:
                x_shock = float(brentq(lambda s: float(A_fn(s)) - Ash, x[it], x[-1], xtol=1e-14)) \
                    if (A[it] - Ash) * (A[-1] - Ash) < 0 else float(x[it] if Ash <= A[it] else x[-1])
            else:
                x_shock = float(np.interp(Ash, A[it:], x[it:]))
            behind = (idx > it) & (x > x_shock)
            M2s = np.atleast_1d(mach_from_area_ratio(np.maximum(A * k / At, 1.0), gamma, branch="subsonic"))
            M = np.where(idx <= it, Msub, np.where(behind, M2s, Msup))
            p = 1.0 / np.atleast_1d(_F(stagnation_pressure_ratio(M, gamma)))
            p = np.where(behind, k * p, p)
        else:
            M = np.where(idx <= it, Msub, Msup)
            p = 1.0 / np.atleast_1d(_F(stagnation_pressure_ratio(M, gamma)))
    return dict(x=x, A=A, M=M, p_over_p0=p, regime=reg, x_shock=x_shock, A_shock_over_At=st["A_shock_over_At"],
                M1_shock=st["M1_shock"], p02_over_p01=st["p02_over_p01"], M_exit=st["M_exit"],
                p_exit_over_p0=st["p_exit_over_p0"], i_throat=it, x_throat=float(x[it]),
                mdot_over_choked=st["mdot_over_choked"])


# ======================================================================================================================
# §15.8 — constant-area duct with friction and heat (the book's f and q; "Fanno"/"Rayleigh" forms are ours)
# ======================================================================================================================
def duct_parameter_A(M1, f=0.0, q=0.0, gamma: float = 1.4):
    """The pure number A of the duct solution: A ≡ M₁²(1 + ½(γ − 1)M₁² + q)/(1 + γM₁² − f)².

    Book: §15.8 (unnumbered display before (15.46)).  **A is not an area** (trap T2).
    Parameters: M1 inlet Mach number [–]; f = F/(p₁A_duct) friction parameter [–] (F = total retarding wall force,
    f ≥ 0); q = Q/h₁ heat parameter [–] (either sign); gamma [–].
    Returns A [–] (NaN where 1 + γM₁² − f ≤ 0: the wall force would exceed the inlet impulse).
    Assumptions: steady flow of a perfect gas in a duct of constant area.
    Validation: V4, V2, V1, V7 — tests/test_ch15.py: test_duct_outlet_mach_V4_the_three_laws_hold_on_both_roots,
    test_duct_outlet_mach_V2_derivation,
    test_duct_outlet_mach_V1_without_friction_or_heat_returns_the_inlet_or_the_shock,
    test_duct_outlet_mach_V7_removable_singularity_of_the_printed_form.
    Label: symbolic, conserved, analytic.
    """
    M1, f, q = _F(M1), _F(f), _F(q)
    den = 1.0 + gamma * M1 ** 2 - f
    with np.errstate(divide="ignore", invalid="ignore"):
        return _S(np.where(den > 0.0, M1 ** 2 * (1.0 + 0.5 * (gamma - 1.0) * M1 ** 2 + q) / den ** 2, np.nan))


def duct_outlet_mach(M1, f=0.0, q=0.0, gamma: float = 1.4, *, branch: str = "continuous"):
    """Outlet Mach number of a constant-area duct with wall friction f and heat exchange q.

    Book: §15.8, Eq. (15.46)  M₂² = {−(1 − 2Aγ) ± [1 − 2A(γ + 1)]^{1/2}}/{(γ − 1) − 2Aγ²}, with A from
    :func:`duct_parameter_A`; Figs. 15.15, 15.16.
    Parameters
    ----------
    M1 : inlet Mach number [–] > 0.   f : F/(p₁A_duct) [–].   q : Q/h₁ [–].   gamma : [–].
    branch : "continuous" — the root that continues M₂ = M₁ as f, q → 0 (subsonic stays subsonic, supersonic stays
        supersonic): the "+" sign for M₁ < 1 and the "−" sign for M₁ > 1 (trap T12; at M₁ = 1 the "+" root);
        "shocked" — the other root (for M₁ > 1 the state behind a normal shock somewhere in the duct; for M₁ < 1 the
        supersonic state excluded by the second law);  "plus" / "minus" — the sign of (15.46) as printed;
        "both" — tuple (continuous, shocked).

    Returns
    -------
    M₂ [–]; **NaN where no real root exists** — the discriminant 1 − 2A(γ + 1) < 0 means the duct is choked (the
    assumed f or q cannot be passed at this M₁), or the root M₂² is negative.
    Method: each root is evaluated both as printed and in the conjugate form M₂² = 2A/[(1 − 2Aγ) ± √D] (the two are
    identical because (1 − 2Aγ)² − D = −2A[(γ − 1) − 2Aγ²]); the form with the larger denominator is used, which
    removes the 0/0 of the printed form at 2Aγ² = γ − 1.  A discriminant between −1e-13 and 0 is treated as 0 (the
    choking limit itself, where round-off of the O(1) terms decides the sign); anything more negative is choked.
    Assumptions: perfect gas; f lumps the whole duct into one number referred to the inlet pressure (trap T11).
    Validation: V4, V2, V1, V7 — tests/test_ch15.py: test_duct_outlet_mach_V4_the_three_laws_hold_on_both_roots,
    test_duct_outlet_mach_V2_derivation,
    test_duct_outlet_mach_V1_without_friction_or_heat_returns_the_inlet_or_the_shock,
    test_duct_outlet_mach_V7_removable_singularity_of_the_printed_form.
    Label: symbolic, conserved, analytic.
    """
    M1 = _F(M1)
    A = _F(duct_parameter_A(M1, f, q, gamma))
    b = 1.0 - 2.0 * A * gamma
    D = 1.0 - 2.0 * A * (gamma + 1.0)
    D = np.where((D < 0.0) & (D > -1e-13), 0.0, D)  # exactly at the choking limit D = 0 up to round-off of O(1) terms
    den = (gamma - 1.0) - 2.0 * A * gamma ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        rD = np.sqrt(np.where(D >= 0.0, D, np.nan))

        def root(sign: float):
            d2 = b + sign * rD
            X_printed = (-b + sign * rD) / den  # Eq. (15.46)
            X_conj = 2.0 * A / d2
            X = np.where(np.abs(d2) >= np.abs(den), X_conj, X_printed)
            return np.where(X >= 0.0, np.sqrt(np.where(X >= 0.0, X, np.nan)), np.nan)

        plus, minus = root(+1.0), root(-1.0)
    sub = M1 <= 1.0
    cont, shock = np.where(sub, plus, minus), np.where(sub, minus, plus)
    if branch == "continuous":
        return _S(cont)
    if branch == "shocked":
        return _S(shock)
    if branch == "plus":
        return _S(plus)
    if branch == "minus":
        return _S(minus)
    if branch == "both":
        return _S(cont), _S(shock)
    raise ValueError('branch must be "continuous", "shocked", "plus", "minus" or "both"')


def duct_state_ratios(M1, M2, f=0.0, q=0.0, gamma: float = 1.4) -> dict:
    """Outlet-to-inlet ratios of a constant-area duct once M₂ is known.

    Book: §15.8 (Mach-number forms of the momentum and energy laws, unnumbered, after (15.45)):
    p₁(1 + γM₁² − f) = p₂(1 + γM₂²),  h₁(1 + ½(γ − 1)M₁² + q) = h₂(1 + ½(γ − 1)M₂²).
    Parameters: M1, M2 [–]; f, q [–]; gamma [–].
    Returns dict [–]: ``p_ratio`` (p₂/p₁), ``T_ratio``, ``rho_ratio``, ``u_ratio`` (u₂/u₁), ``p0_ratio``,
    ``T0_ratio``, ``dS_over_cv``, ``dS_over_cp``.
    Assumptions: perfect gas; M₂ from :func:`duct_outlet_mach` with the same f, q.
    Validation: V4, V1 — tests/test_ch15.py: test_duct_outlet_mach_V4_the_three_laws_hold_on_both_roots,
    test_fanno_and_rayleigh_V1_bridge_to_the_book_parameters_f_and_q.
    Label: conserved, analytic.
    """
    M1, M2, f, q = _F(M1), _F(M2), _F(f), _F(q)
    t1, t2 = 1.0 + 0.5 * (gamma - 1.0) * M1 ** 2, 1.0 + 0.5 * (gamma - 1.0) * M2 ** 2
    pr = (1.0 + gamma * M1 ** 2 - f) / (1.0 + gamma * M2 ** 2)
    Tr = (t1 + q) / t2
    rr = pr / Tr
    with np.errstate(divide="ignore", invalid="ignore"):
        ds = np.log(pr) - gamma * np.log(rr)
    out = dict(p_ratio=pr, T_ratio=Tr, rho_ratio=rr, u_ratio=1.0 / rr, p0_ratio=pr * (t2 / t1) ** (gamma / (gamma - 1.0)),
               T0_ratio=(t1 + q) / t1, dS_over_cv=ds, dS_over_cp=ds / gamma)
    return {k: _S(v) for k, v in out.items()}


def duct_entropy_change(M1, M2, f=0.0, q=0.0, gamma: float = 1.4):
    """(S₂ − S₁)/C_v between inlet and outlet of the duct, ln[(p₂/p₁)(ρ₁/ρ₂)^γ].

    Book: §15.8 (the shaded "ΔS < 0" regions of Figs. 15.15, 15.16) with the entropy line of Eq. (15.1).
    Parameters: M1, M2 [–]; f, q [–]; gamma [–].   Returns (S₂ − S₁)/C_v [–].
    With q = 0 a negative value marks a root excluded by the second law; with q ≠ 0 the entropy may fall (cooling).
    Assumptions: perfect gas.
    Validation: V4, V1 — tests/test_ch15.py: test_duct_outlet_mach_V4_the_three_laws_hold_on_both_roots,
    test_duct_outlet_mach_V1_without_friction_or_heat_returns_the_inlet_or_the_shock.
    Label: conserved, analytic.
    """
    return duct_state_ratios(M1, M2, f, q, gamma)["dS_over_cv"]


def duct_residuals(state1: dict, state2: dict, f=0.0, q=0.0, gamma: float = 1.4) -> dict:
    """Relative residuals of the three duct laws with friction and heat (zero for a consistent pair of states).

    Book: §15.8, Eq. (15.45)  ρ₁u₁ = ρ₂u₂,  p₁ + ρ₁u₁² = p₂ + ρ₂u₂² + p₁f,  h₁ + ½u₁² + h₁q = h₂ + ½u₂².
    Parameters: state1, state2 dicts (``rho``, ``u``, ``p``); f = F/(p₁A) [–]; q = Q/h₁ [–]; gamma [–].
    Returns dict ``mass``, ``momentum``, ``energy`` [–] (right side − left side, over the left side).
    Assumptions: perfect gas (h = γp/((γ − 1)ρ)); constant area.
    Validation: V4 — tests/test_ch15.py: test_duct_outlet_mach_V4_the_three_laws_hold_on_both_roots.
    Label: conserved.
    """
    r1, u1, p1 = _F(state1["rho"]), _F(state1["u"]), _F(state1["p"])
    r2, u2, p2 = _F(state2["rho"]), _F(state2["u"]), _F(state2["p"])
    h1, h2 = gamma / (gamma - 1.0) * p1 / r1, gamma / (gamma - 1.0) * p2 / r2
    m1, j1, e1 = r1 * u1, p1 + r1 * u1 ** 2, h1 + 0.5 * u1 ** 2 + h1 * _F(q)
    return dict(mass=_S((r2 * u2 - m1) / m1), momentum=_S((p2 + r2 * u2 ** 2 + p1 * _F(f) - j1) / j1),
                energy=_S((h2 + 0.5 * u2 ** 2 - e1) / e1))  # Eq. (15.45)


def duct_choking_limits(M1, gamma: float = 1.4) -> dict:
    """Largest friction parameter (q = 0) and largest heat parameter (f = 0) a duct can take at a given M₁ (ours).

    Book: §15.8 (choking by friction or heat; Exercise 15.16).  The discriminant of (15.46) vanishes, 2A(γ + 1) = 1,
    and M₂ = 1:  f_max = 1 + γM₁² − M₁√(2(γ + 1)(1 + ½(γ − 1)M₁²)),
    q_max = (1 + γM₁²)²/(2(γ + 1)M₁²) − 1 − ½(γ − 1)M₁².  Both vanish at M₁ = 1.
    Parameters: M1 [–] > 0; gamma [–].   Returns dict ``f_max``, ``q_max`` [–].
    Assumptions: as :func:`duct_outlet_mach`.
    Validation: V1 — at the limit duct_outlet_mach returns 1 to √(round-off) (|M₂ − 1| ≤ 1.3e-8 measured, bound ½(γ +
    1)√(4 ulp)) on both roots, and NaN just beyond — tests/test_ch15.py:
    test_duct_choking_limits_V1_outlet_is_sonic_at_the_limit_and_no_root_beyond.
    Label: analytic.
    """
    M1 = _F(M1)
    t = 1.0 + 0.5 * (gamma - 1.0) * M1 ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        fmax = 1.0 + gamma * M1 ** 2 - M1 * np.sqrt(2.0 * (gamma + 1.0) * t)
        qmax = (1.0 + gamma * M1 ** 2) ** 2 / (2.0 * (gamma + 1.0) * M1 ** 2) - t
    return dict(f_max=_S(fmax), q_max=_S(qmax))


def duct_curve(M1, f=0.0, q=0.0, gamma: float = 1.4) -> dict:
    """One curve of the M₂–M₁ charts: both roots over an array of inlet Mach numbers, NaN where choked.

    Book: §15.8, Eq. (15.46); Figs. 15.15 (f as the parameter) and 15.16 (q as the parameter).
    Parameters: M1 array [–]; f, q [–] (scalars); gamma [–].
    Returns dict: ``M1``, ``M2_continuous``, ``M2_shocked``, ``choked`` (bool array: no real root),
    ``dS_continuous``, ``dS_shocked`` ((S₂ − S₁)/C_v on each root).
    Assumptions: as :func:`duct_outlet_mach`.
    Validation: V1 — tests/test_ch15.py:
    test_duct_outlet_mach_V1_without_friction_or_heat_returns_the_inlet_or_the_shock.
    Label: analytic.
    """
    M1 = np.atleast_1d(_F(M1))
    c, s = duct_outlet_mach(M1, f, q, gamma, branch="both")
    c, s = np.atleast_1d(c), np.atleast_1d(s)
    return dict(M1=M1, M2_continuous=c, M2_shocked=s, choked=np.isnan(c) & np.isnan(s),
                dS_continuous=np.atleast_1d(duct_entropy_change(M1, c, f, q, gamma)),
                dS_shocked=np.atleast_1d(duct_entropy_change(M1, s, f, q, gamma)))


def duct_choking_locus(gamma: float = 1.4, M1=None) -> dict:
    """The choking limits as curves over M₁ (the envelope where the curves of Figs. 15.15, 15.16 turn round).

    Book: §15.8 (limiting curves meet M₂ = 1 at right angles); ours, from :func:`duct_choking_limits`.
    Parameters: gamma [–]; M1 array [–] (default 0.05 … 4, 400 points).
    Returns dict ``M1``, ``f_max``, ``q_max`` [–].
    Assumptions: as :func:`duct_choking_limits`.
    Validation: V1 — tests/test_ch15.py: test_duct_choking_limits_V1_outlet_is_sonic_at_the_limit_and_no_root_beyond.
    Label: analytic.
    """
    M1 = np.linspace(0.05, 4.0, 400) if M1 is None else np.atleast_1d(_F(M1))
    lim = duct_choking_limits(M1, gamma)
    return dict(M1=M1, f_max=np.atleast_1d(lim["f_max"]), q_max=np.atleast_1d(lim["q_max"]))


def rayleigh_temperature_sign(M, gamma: float = 1.4):
    """Sign of dT/dq in a frictionless constant-area duct: −1 for 1/γ < M² < 1 (adding heat *cools* the gas), +1
    elsewhere, 0 at M² = 1/γ.

    Book: §15.8, Eq. (15.48) first relation  dT/T₁ = [(1 − γM²)/(1 − M²)] dq and the text ("1/γ < M₁² < 1").
    Parameters: M [–]; gamma [–].   Returns −1.0, 0.0 or +1.0 (NaN at M = 1).
    Assumptions: perfect gas, no friction.
    Validation: V1, V7 — tests/test_ch15.py: test_small_helpers_V1_against_independent_closed_forms,
    test_influence_coefficients_V7_limits_and_sign_pattern.
    Label: analytic.
    """
    M2 = _F(M) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        return _S(np.where(M2 == 1.0, np.nan, np.sign((1.0 - gamma * M2) / (1.0 - M2))))  # Eq. (15.48)


def fanno_ratios(M, gamma: float = 1.4) -> dict:
    """Adiabatic constant-area flow with friction in the usual reference-state form (**not in the book — ours**).

    Book: §15.8, Eqs. (15.45), (15.47) (the laws these forms integrate; the forms themselves are not printed).
    The book never says "Fanno flow"; it uses the lumped f of (15.45).  Here every ratio is to the **sonic state of
    the same frictional process** (a star that is *not* the isentropic A* — trap T2):
    T/T* = (γ + 1)/(2 + (γ − 1)M²) (T₀ constant), p/p* = (1/M)√(T/T*) (mass), ρ/ρ* = (1/M)√((2 + (γ − 1)M²)/(γ + 1)),
    p₀/p₀* = (1/M)[(2 + (γ − 1)M²)/(γ + 1)]^{(γ+1)/(2(γ−1))}, and the friction length to choking
    4f_F L*/D = (1 − M²)/(γM²) + [(γ + 1)/(2γ)] ln[(γ + 1)M²/(2 + (γ − 1)M²)]  (f_F = Fanning friction factor, D =
    hydraulic diameter; obtained by integrating the corrected (15.47) with df = ½γM² d(4f_F x/D)).
    Parameters: M [–] > 0; gamma [–].
    Returns dict [–]: ``T``, ``p``, ``rho``, ``u`` (u/u*), ``p0``, ``fL`` (4f_F L*/D ≥ 0).
    Assumptions: perfect gas, adiabatic, constant area, constant friction factor.
    Validation: V1, V5 — tests/test_ch15.py: test_fanno_and_rayleigh_V1_bridge_to_the_book_parameters_f_and_q,
    test_fanno_rayleigh_V5_cited_reference_state_forms.
    Label: analytic, benchmark.
    """
    M = _F(M)
    t = (gamma + 1.0) / (2.0 + (gamma - 1.0) * M ** 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = dict(T=t, p=np.sqrt(t) / M, rho=1.0 / (M * np.sqrt(t)), u=M * np.sqrt(t),
                   p0=t ** (-(gamma + 1.0) / (2.0 * (gamma - 1.0))) / M,
                   fL=(1.0 - M ** 2) / (gamma * M ** 2) + (gamma + 1.0) / (2.0 * gamma) * np.log(t * M ** 2))
    return {k: _S(v) for k, v in out.items()}


def fanno_length_to_f(M1, M2, gamma: float = 1.4) -> dict:
    """Bridge between the usual friction length and the book's lumped friction parameter (ours).

    Book: §15.8, Eq. (15.45) (definition f = F/(p₁A)) — the link is not in the book.  For an adiabatic duct that takes
    the flow from M₁ to M₂:  4f_F L/D = fL(M₁) − fL(M₂) (:func:`fanno_ratios`), and the momentum law gives
    f = 1 + γM₁² − (p₂/p₁)(1 + γM₂²) with p₂/p₁ = (M₁/M₂)√(T₂/T₁).
    Parameters: M1, M2 [–] on the same side of 1 with M₂ between M₁ and 1; gamma [–].
    Returns dict [–]: ``friction_length`` (4f_F L/D), ``f`` (the book's parameter for the whole duct).
    Note (trap T11): f refers to the *inlet* pressure, so two ducts in series do not add their f; friction lengths do.
    Assumptions: as :func:`fanno_ratios`.
    Validation: V1 — tests/test_ch15.py: test_fanno_and_rayleigh_V1_bridge_to_the_book_parameters_f_and_q.
    Label: analytic.
    """
    a, b = fanno_ratios(M1, gamma), fanno_ratios(M2, gamma)
    M1, M2 = _F(M1), _F(M2)
    pr = _F(b["p"]) / _F(a["p"])
    return dict(friction_length=_S(_F(a["fL"]) - _F(b["fL"])),
                f=_S(1.0 + gamma * M1 ** 2 - pr * (1.0 + gamma * M2 ** 2)))


def rayleigh_ratios(M, gamma: float = 1.4) -> dict:
    """Frictionless constant-area flow with heat exchange in the usual reference-state form (**not in the book —
    ours**).

    Book: §15.8, Eqs. (15.45), (15.48) (the laws these forms integrate; the forms themselves are not printed).
    The book never says "Rayleigh flow"; it uses q = Q/h₁ in (15.45).  Every ratio is to the **sonic state of the
    same heated process** (not the isentropic star — trap T2):  p/p* = (1 + γ)/(1 + γM²) (momentum),
    T/T* = M²(1 + γ)²/(1 + γM²)² (mass), u/u* = ρ*/ρ = (1 + γ)M²/(1 + γM²),
    T₀/T₀* = 2(γ + 1)M²(1 + ½(γ − 1)M²)/(1 + γM²)², p₀/p₀* = [(1 + γ)/(1 + γM²)][(2 + (γ − 1)M²)/(γ + 1)]^{γ/(γ−1)}.
    Parameters: M [–] ≥ 0; gamma [–].
    Returns dict [–]: ``T``, ``p``, ``rho``, ``u``, ``T0``, ``p0``.  T/T* is largest at M = 1/√γ.
    Assumptions: perfect gas, no friction, constant area.
    Validation: V1, V5 — tests/test_ch15.py: test_fanno_and_rayleigh_V1_bridge_to_the_book_parameters_f_and_q,
    test_fanno_rayleigh_V5_cited_reference_state_forms.
    Label: analytic, benchmark.
    """
    M = _F(M)
    k = (1.0 + gamma) / (1.0 + gamma * M ** 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = dict(T=(M * k) ** 2, p=k, rho=1.0 / (k * M ** 2), u=k * M ** 2,
                   T0=2.0 * (gamma + 1.0) * M ** 2 * (1.0 + 0.5 * (gamma - 1.0) * M ** 2) / (1.0 + gamma * M ** 2) ** 2,
                   p0=k * ((2.0 + (gamma - 1.0) * M ** 2) / (gamma + 1.0)) ** (gamma / (gamma - 1.0)))
    return {k_: _S(v) for k_, v in out.items()}


def _line_Ts(ratios: dict, M, gamma: float) -> dict:
    with np.errstate(divide="ignore", invalid="ignore"):
        ds = np.log(_F(ratios["T"])) - (gamma - 1.0) / gamma * np.log(_F(ratios["p"]))
    return dict(M=_S(_F(M)), ds_over_cp=_S(ds), T_over_Tstar=ratios["T"])


def fanno_line_Ts(M, gamma: float = 1.4) -> dict:
    """The friction (Fanno) line in the T–s plane: (s − s*)/C_p and T/T* along the process (ours).

    Book: §15.8 (friction drives the flow to M = 1; the entropy is largest there) — the T–s picture is not in the
    book.  (s − s*)/C_p = ln(T/T*) − [(γ − 1)/γ] ln(p/p*) with :func:`fanno_ratios`; the star is the sonic state of
    the same frictional process.
    Parameters: M [–] > 0; gamma [–].   Returns dict ``M``, ``ds_over_cp`` (≤ 0, maximum 0 at M = 1), ``T_over_Tstar``.
    Assumptions: as :func:`fanno_ratios`.
    Validation: V1, V5 — tests/test_ch15.py: test_fanno_rayleigh_V5_cited_reference_state_forms.
    Label: analytic, benchmark.
    """
    return _line_Ts(fanno_ratios(M, gamma), M, gamma)


def rayleigh_line_Ts(M, gamma: float = 1.4) -> dict:
    """The heating (Rayleigh) line in the T–s plane: (s − s*)/C_p and T/T* along the process (ours).

    Book: §15.8 (heating drives the flow to M = 1; between M² = 1/γ and 1 the temperature falls while heat is added)
    — the T–s picture is not in the book.  The star is the sonic state of the same heated process.
    Parameters: M [–] > 0; gamma [–].   Returns dict ``M``, ``ds_over_cp`` (maximum 0 at M = 1), ``T_over_Tstar``
    (maximum at M = 1/√γ).
    Assumptions: as :func:`rayleigh_ratios`.
    Validation: V1, V5 — tests/test_ch15.py: test_fanno_rayleigh_V5_cited_reference_state_forms.
    Label: analytic, benchmark.
    """
    return _line_Ts(rayleigh_ratios(M, gamma), M, gamma)


# ======================================================================================================================
# §15.9 — Mach angle, oblique shocks, weak waves, Prandtl–Meyer expansion
# ======================================================================================================================
def mach_angle(M):
    """Mach angle μ = sin⁻¹(1/M): half-angle of the Mach cone (or wedge).

    Book: §15.9, Eq. (15.49)  sin μ = 1/M.
    Parameters: M [–].   Returns μ [rad]; π/2 at M = 1; **NaN for M < 1** (no cone: the source is subsonic).
    Assumptions: weak (acoustic) disturbances in a uniform stream.
    Validation: V1, V5 — tests/test_ch15.py: test_mach_cone_V1_envelope_of_the_wave_fronts,
    test_mach_angle_and_prandtl_meyer_V5_naca_report_1135_rows.
    Label: analytic, benchmark.
    """
    M = _F(M)
    with np.errstate(divide="ignore", invalid="ignore"):
        return _S(np.where(M >= 1.0, np.arcsin(1.0 / np.where(M >= 1.0, M, 1.0)), np.nan))  # Eq. (15.49)


def mach_from_angle(mu):
    """Mach number from the Mach angle, M = 1/sin μ.

    Book: §15.9, Eq. (15.49) inverted.
    Parameters: mu [rad] in (0, π/2].   Returns M [–] (NaN outside the range).
    Assumptions: as :func:`mach_angle`.
    Validation: V1, V5 — tests/test_ch15.py: test_mach_cone_V1_envelope_of_the_wave_fronts,
    test_mach_angle_and_prandtl_meyer_V5_naca_report_1135_rows.
    Label: analytic, benchmark.
    """
    mu = _F(mu)
    with np.errstate(divide="ignore", invalid="ignore"):
        return _S(np.where((mu > 0.0) & (mu <= 0.5 * np.pi + 1e-15), 1.0 / np.sin(mu), np.nan))


_SIGMA_RTOL = 1e-12  # round-off allowance on M₁ sin σ ≥ 1 and on σ ≤ π/2 (oblique-shock functions that take sigma)


def deflection_angle(M1, sigma, gamma: float = 1.4):
    """Flow deflection δ behind an oblique shock of angle σ: tan δ = 2 cot σ (M₁² sin²σ − 1)/(M₁²(γ + cos 2σ) + 2).

    Book: §15.9, Eq. (15.53); Fig. 15.19.  (σ, δ are β, θ in most other sources.)
    Parameters: M1 upstream Mach number [–] > 1; sigma shock angle from the upstream flow direction [rad];
    gamma [–].   Returns δ [rad] (0 at σ = μ₁ and at σ = π/2; negative on 0 < σ < μ₁ — not a shock).
    **The bare formula**: like :func:`normal_shock_mach` it does not refuse an angle outside [μ₁, π/2].  On
    0 < σ < μ₁ the value is negative — the mark of the forbidden expansion root (M₁ sin σ < 1, entropy would fall) —
    but elsewhere outside the range (σ > π/2, σ ≤ 0, M₁ < 1) the sign follows the trigonometry, not the physics:
    **the sign of δ is not a validity test outside [μ₁, π/2]; use** :func:`oblique_shock` **(raises / flags) or**
    :func:`oblique_shock_components` **(NaN)**, which check M₁ sin σ ≥ 1 and σ ≤ π/2.
    Assumptions: perfect gas; plane shock; uniform flow on each side.
    Validation: V2, V5, V1, V7 — tests/test_ch15.py: test_deflection_relation_V2_derivation,
    test_oblique_shock_V5_naca_report_1135_forms, test_oblique_shock_V1_two_roots_maximum_deflection_and_detachment,
    test_deflection_angle_V7_bare_formula_outside_the_range_and_its_guarded_callers.
    Label: symbolic, benchmark, analytic.
    """
    M1, s = _F(M1), _F(sigma)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = 2.0 / np.tan(s) * (M1 ** 2 * np.sin(s) ** 2 - 1.0) / (M1 ** 2 * (gamma + np.cos(2.0 * s)) + 2.0)  # (15.53)
    return _S(np.arctan(t))


def oblique_shock_components(M1, sigma, gamma: float = 1.4) -> dict:
    """An oblique shock resolved along and across the shock: a normal shock on the normal components plus an
    unchanged shock-parallel speed.

    Book: §15.9 (unnumbered displays with Fig. 15.18): (u₁, v) = V₁(sin σ, cos σ), (u₂, v) = V₂(sin(σ − δ),
    cos(σ − δ)); M_n1 = M₁ sin σ > 1, M_n2 = M₂ sin(σ − δ) < 1.
    Parameters: M1 [–] > 1; sigma [rad] in [μ₁, π/2]; gamma [–].
    Returns dict: ``Mn1``, ``Mt1`` (= M₁ cos σ), ``Mn2``, ``M2``, ``delta`` [rad], ``u1_over_V1`` (sin σ),
    ``v_over_V1`` (cos σ), ``u2_over_V1``, ``V2_over_V1``, ``c2_over_c1``.
    **Outside the range** — M₁ sin σ < 1 (σ below the Mach angle, or M₁ < 1: the entropy-lowering "expansion shock",
    trap T6 for oblique shocks) or σ > π/2 — **every key is NaN** (scalars and arrays alike).  Round-off is forgiven:
    1 − 10⁻¹² ≤ M₁ sin σ < 1 is treated as the Mach wave (M_n1 = 1, δ = 0) and σ may exceed π/2 by 10⁻¹² rad
    (``_SIGMA_RTOL``).
    Assumptions: as :func:`deflection_angle`.
    Validation: V4, V5, V1, V7 — tests/test_ch15.py: test_oblique_shock_V4_jump_laws_in_two_dimensions,
    test_oblique_shock_V5_naca_report_1135_forms, test_oblique_shock_V1_two_roots_maximum_deflection_and_detachment,
    test_oblique_shock_V7_given_shock_angle_outside_its_range_scalars_raise_arrays_flag,
    test_deflection_angle_V7_bare_formula_outside_the_range_and_its_guarded_callers.
    Label: conserved, benchmark, analytic.
    """
    M1, s = np.broadcast_arrays(_F(M1), _F(sigma))
    with np.errstate(invalid="ignore"):
        Mn1 = M1 * np.sin(s)
        mach_wave = (Mn1 < 1.0) & (Mn1 >= 1.0 - _SIGMA_RTOL)  # σ = μ₁ to round-off
        ok = (Mn1 >= 1.0 - _SIGMA_RTOL) & (s <= 0.5 * np.pi + _SIGMA_RTOL)
    Mn1 = np.where(mach_wave, 1.0, np.where(ok, Mn1, np.nan))
    s = np.where(ok, s, np.nan)
    rr = _F(normal_shock_density_ratio(Mn1, gamma))
    Mn2 = _F(normal_shock_mach(Mn1, gamma))
    d = np.where(mach_wave, 0.0, _F(deflection_angle(M1, s, gamma)))
    u2 = np.sin(s) / rr
    V2 = np.hypot(u2, np.cos(s))
    with np.errstate(divide="ignore", invalid="ignore"):
        M2 = Mn2 / np.sin(s - d)  # Eq. (15.52)
    out = dict(Mn1=Mn1, Mt1=M1 * np.cos(s), Mn2=Mn2, M2=M2, delta=d, u1_over_V1=np.sin(s), v_over_V1=np.cos(s),
               u2_over_V1=u2, V2_over_V1=V2, c2_over_c1=np.sqrt(_F(normal_shock_temperature_ratio(Mn1, gamma))))
    return {k: _S(v) for k, v in out.items()}


def max_deflection(M1, gamma: float = 1.4) -> dict:
    """Largest deflection an attached oblique shock can produce, and the shock angle at which it occurs (ours).

    Book: §15.9, Fig. 15.19 (δ_max read from the chart; no formula).  Setting dδ/dσ = 0 in (15.53) gives a quadratic
    in sin²σ:  sin²σ_max = {(γ + 1)M₁² − 4 + √((γ + 1)[(γ + 1)M₁⁴ + 8(γ − 1)M₁² + 16])}/(4γM₁²).
    Parameters: M1 [–] ≥ 1; gamma [–].
    Returns dict ``delta_max`` [rad], ``sigma_max`` [rad] (δ_max → sin⁻¹(1/γ) as M₁ → ∞: 45.58° for γ = 1.4).
    Assumptions: as :func:`deflection_angle`.
    Validation: V5, V1, V7 (limits M₁ → 1 and M₁ → ∞) — tests/test_ch15.py:
    test_oblique_shock_V5_naca_report_1135_forms, test_oblique_shock_V1_two_roots_maximum_deflection_and_detachment.
    Label: benchmark, analytic.
    """
    M1 = _F(M1)
    g1 = gamma + 1.0
    with np.errstate(divide="ignore", invalid="ignore"):
        s2 = (g1 * M1 ** 2 - 4.0 + np.sqrt(g1 * (g1 * M1 ** 4 + 8.0 * (gamma - 1.0) * M1 ** 2 + 16.0))) / (
            4.0 * gamma * M1 ** 2)
        sig = np.where(M1 >= 1.0, np.arcsin(np.sqrt(np.clip(s2, 0.0, 1.0))), np.nan)
    return dict(delta_max=_S(np.maximum(_F(deflection_angle(M1, sig, gamma)), 0.0)), sigma_max=_S(sig))


def sonic_deflection(M1, gamma: float = 1.4) -> dict:
    """The point of the δ–σ curve where the flow behind the shock is exactly sonic (the M₂ = 1 locus) (ours).

    Book: §15.9, Fig. 15.19 (locus drawn; no formula).  M₂ = 1 in (15.52) with (15.53) gives
    sin²σ_s = {(γ + 1)M₁² − (3 − γ) + √((γ + 1)[(γ + 1)M₁⁴ − 2(3 − γ)M₁² + γ + 9])}/(4γM₁²).
    σ_s < σ_max: weak shocks with δ between δ_sonic and δ_max are subsonic behind (trap T13).
    Parameters: M1 [–] ≥ 1; gamma [–].   Returns dict ``delta_sonic`` [rad], ``sigma_sonic`` [rad].
    Assumptions: as :func:`deflection_angle`.
    Validation: V5, V1 — tests/test_ch15.py: test_oblique_shock_V5_naca_report_1135_forms,
    test_oblique_shock_V1_two_roots_maximum_deflection_and_detachment.
    Label: benchmark, analytic.
    """
    M1 = _F(M1)
    g1 = gamma + 1.0
    with np.errstate(divide="ignore", invalid="ignore"):
        s2 = (g1 * M1 ** 2 - (3.0 - gamma) + np.sqrt(g1 * (g1 * M1 ** 4 - 2.0 * (3.0 - gamma) * M1 ** 2 + gamma + 9.0))
              ) / (4.0 * gamma * M1 ** 2)
        sig = np.where(M1 >= 1.0, np.arcsin(np.sqrt(np.clip(s2, 0.0, 1.0))), np.nan)
    return dict(delta_sonic=_S(np.maximum(_F(deflection_angle(M1, sig, gamma)), 0.0)), sigma_sonic=_S(sig))


def _sigma_from_delta(M1: float, delta: float, gamma: float, branch: str, raise_detached: bool) -> float:
    if not (np.isfinite(M1) and np.isfinite(delta)):
        return np.nan
    if M1 < 1.0:
        raise ValueError("oblique_shock needs a supersonic upstream Mach number M1 >= 1")
    if delta < 0.0:
        raise ValueError("delta < 0 is a turn away from the flow: an expansion (use prandtl_meyer_turn), not a shock")
    md = max_deflection(M1, gamma)
    mu = float(np.arcsin(1.0 / M1))
    if delta > md["delta_max"]:
        if delta > md["delta_max"] * (1.0 + 1e-12) + 1e-15:
            if raise_detached:
                raise ValueError(f"detached: delta = {np.degrees(delta):.4f} deg exceeds delta_max = "
                                 f"{np.degrees(md['delta_max']):.4f} deg at M1 = {M1:g}; no attached oblique shock")
            return np.nan
        return float(md["sigma_max"])
    fun = lambda s: float(deflection_angle(M1, s, gamma)) - delta  # noqa: E731
    if branch == "weak":
        if delta == 0.0:
            return mu
        lo, hi = mu, float(md["sigma_max"])
    elif branch == "strong":
        if delta == 0.0:
            return 0.5 * np.pi
        lo, hi = float(md["sigma_max"]), 0.5 * np.pi
    else:
        raise ValueError('branch must be "weak" or "strong"')
    if fun(lo) * fun(hi) > 0.0:  # delta within round-off of delta_max
        return float(md["sigma_max"])
    root = float(brentq(fun, lo, hi, xtol=1e-15, rtol=1e-15))
    assert abs(fun(root)) < 1e-10, "oblique-shock root did not converge"
    return root


def oblique_shock(M1, *, delta=None, sigma=None, branch: str = "weak", gamma: float = 1.4) -> dict:
    """An oblique shock from the upstream Mach number and either the deflection δ or the shock angle σ.

    Book: §15.9, Eqs. (15.50) p₂/p₁ = 1 + [2γ/(γ + 1)](M₁² sin²σ − 1), (15.51) ρ₂/ρ₁ = tan σ/tan(σ − δ),
    (15.52) M₂² sin²(σ − δ) = [(γ − 1)M₁² sin²σ + 2]/[2γM₁² sin²σ + 1 − γ], (15.53); Figs. 15.18–15.20.
    Parameters
    ----------
    M1 : upstream Mach number [–] ≥ 1.
    delta : deflection angle [rad] ≥ 0 (θ elsewhere), **or** sigma : shock angle [rad] (β elsewhere) in [μ₁, π/2] —
        give exactly one.
    branch : with ``delta``: "weak" (smaller σ; the attached shock of a wedge) or "strong" (larger σ, always subsonic
        behind).  "Weak"/"strong" name the two roots, not the size of the compression (trap T13).
    gamma : [–].

    Returns
    -------
    dict: ``sigma`` [rad], ``delta`` [rad], ``M2``, ``Mn1``, ``Mn2``, ``p_ratio``, ``rho_ratio``, ``T_ratio``,
    ``p0_ratio``, ``dS_over_cv`` [–], ``attached`` (bool), ``branch``.
    Scalar input with δ > δ_max raises ``ValueError("detached: …")``; array input returns NaN there and
    ``attached=False``.
    **A given ``sigma`` outside [μ₁, π/2] is refused in the same way** (trap T6 for oblique shocks): M₁ sin σ < 1
    would send a subsonic normal Mach number through the jump formulas and return an "expansion shock" with
    p₂/p₁ < 1 and S₂ < S₁.  Scalar input raises ``ValueError("not a shock: …")``; array input returns NaN in every
    numerical key and ``attached=False`` at those entries.  Round-off is forgiven as in
    :func:`oblique_shock_components` (σ = μ₁ to 10⁻¹² in M₁ sin σ is the Mach wave: no jump, δ = 0).
    Method: closed forms through the normal-shock functions at M_n1 = M₁ sin σ; with ``delta`` the angle is the
    ``brentq`` root of (15.53) on (μ₁, σ_max) or (σ_max, π/2), residual asserted.
    Assumptions: as :func:`deflection_angle`.
    Validation: V4, V5, V1, V7 — four jump laws in two dimensions (135 shocks, < 2e-11); both roots, δ_max and
    detachment; guard σ < μ₁ refused — tests/test_ch15.py: test_oblique_shock_V4_jump_laws_in_two_dimensions,
    test_oblique_shock_V5_naca_report_1135_forms, test_oblique_shock_V1_two_roots_maximum_deflection_and_detachment,
    test_oblique_shock_V7_shock_angle_below_the_mach_angle_is_refused,
    test_oblique_shock_V7_given_shock_angle_outside_its_range_scalars_raise_arrays_flag,
    test_deflection_angle_V7_bare_formula_outside_the_range_and_its_guarded_callers.
    Label: conserved, benchmark, analytic.
    """
    if (delta is None) == (sigma is None):
        raise ValueError("give exactly one of delta and sigma")
    M1a = _F(M1)
    if sigma is None:
        scalar = M1a.ndim == 0 and np.ndim(delta) == 0
        sig = _F(np.vectorize(lambda m, d: _sigma_from_delta(float(m), float(d), gamma, branch, scalar),
                              otypes=[float])(M1a, _F(delta)))
    else:
        sig = _F(sigma)
        branch = "given"
        with np.errstate(invalid="ignore"):
            mn = M1a * np.sin(sig)
            bad = (mn < 1.0 - _SIGMA_RTOL) | (sig > 0.5 * np.pi + _SIGMA_RTOL)
        if bad.ndim == 0 and bool(bad):
            raise ValueError(f"not a shock: sigma = {np.degrees(float(sig)):.4f} deg at M1 = {float(M1a):g} gives "
                             f"M1 sin(sigma) = {float(mn):.6g}; an oblique shock needs the shock angle between the "
                             "Mach angle mu1 = asin(1/M1) and 90 deg (M1 sin(sigma) >= 1, else the entropy would "
                             "fall)")
        sig = np.where(bad, np.nan, sig)  # arrays: NaN and attached=False at the refused entries
    comp = oblique_shock_components(M1a, sig, gamma)
    Mn1 = _F(comp["Mn1"])
    att = np.isfinite(_F(comp["Mn1"]) + sig)
    sig = sig + 0.0 * Mn1  # the common (broadcast) shape
    out = dict(sigma=_S(sig), delta=comp["delta"] if sigma is not None else _S(np.where(att, _F(delta), np.nan)),
               M2=comp["M2"], Mn1=comp["Mn1"], Mn2=comp["Mn2"],
               p_ratio=normal_shock_pressure_ratio(Mn1, gamma),  # Eq. (15.50)
               rho_ratio=normal_shock_density_ratio(Mn1, gamma),  # Eq. (15.51)
               T_ratio=normal_shock_temperature_ratio(Mn1, gamma),
               p0_ratio=normal_shock_stagnation_pressure_ratio(Mn1, gamma),
               dS_over_cv=normal_shock_entropy(Mn1, gamma), attached=bool(att) if att.ndim == 0 else att,
               branch=branch)
    return out


def oblique_shock_curve(M1: float, gamma: float = 1.4, n: int = 400) -> dict:
    """One curve of the deflection–shock-angle chart: δ(σ) and M₂(σ) for σ from the Mach angle to π/2.

    Book: §15.9, Eq. (15.53); Fig. 15.19.
    Parameters: M1 [–] > 1; gamma [–]; n number of points (default 400).
    Returns dict: ``sigma`` [rad] (n,), ``delta`` [rad], ``M2`` [–], ``i_max`` and ``i_sonic`` (indices nearest the
    maximum-deflection and M₂ = 1 points), ``delta_max``, ``sigma_max``, ``delta_sonic``, ``sigma_sonic`` [rad],
    ``mu`` [rad].
    The curve stays inside [μ₁, π/2] by construction (its first point is the Mach wave: δ = 0, M₂ = M₁), so no
    entry is an expansion root.  ``M1`` < 1 raises ``ValueError`` (no Mach angle, no curve).
    Assumptions: as :func:`deflection_angle`.
    Validation: V1 — tests/test_ch15.py: test_oblique_shock_V1_two_roots_maximum_deflection_and_detachment.
    Label: analytic.
    """
    if not float(M1) >= 1.0:
        raise ValueError("oblique_shock_curve needs a supersonic upstream Mach number M1 >= 1 (below it there is no "
                         "Mach angle and no oblique shock)")
    mu = float(mach_angle(M1))
    sig = np.linspace(mu, 0.5 * np.pi, int(n))
    c = oblique_shock_components(M1, sig, gamma)
    md, sd = max_deflection(M1, gamma), sonic_deflection(M1, gamma)
    M2 = np.array(c["M2"], dtype=float)
    M2[0] = M1  # at σ = μ the wave is a Mach wave: no change
    d = np.maximum(np.asarray(c["delta"], dtype=float), 0.0)
    return dict(sigma=sig, delta=d, M2=M2, i_max=int(np.argmin(np.abs(sig - md["sigma_max"]))),
                i_sonic=int(np.argmin(np.abs(sig - sd["sigma_sonic"]))), delta_max=md["delta_max"],
                sigma_max=md["sigma_max"], delta_sonic=sd["delta_sonic"], sigma_sonic=sd["sigma_sonic"], mu=mu)


def detachment_mach(delta, gamma: float = 1.4):
    """Smallest upstream Mach number at which a wedge of half-angle δ still carries an attached shock (ours).

    Book: §15.9, Figs. 15.19–15.21 (detachment described; no formula).  Root of δ_max(M₁) = δ; δ_max rises
    monotonically from 0 at M₁ = 1 to sin⁻¹(1/γ) as M₁ → ∞.
    Parameters: delta [rad] ≥ 0; gamma [–].
    Returns M₁ [–]; 1 for δ = 0; ``inf`` for δ ≥ sin⁻¹(1/γ) (never attached).
    Method: ``brentq`` on M₁ with an expanding bracket.
    Assumptions: as :func:`max_deflection`.
    Validation: V1 — tests/test_ch15.py: test_oblique_shock_V1_two_roots_maximum_deflection_and_detachment.
    Label: analytic.
    """
    lim = float(np.arcsin(1.0 / gamma))

    def one(d: float) -> float:
        if d <= 0.0:
            return 1.0
        if d >= lim:
            return np.inf
        fun = lambda m: float(max_deflection(m, gamma)["delta_max"]) - d  # noqa: E731
        hi = 2.0
        while fun(hi) < 0.0:
            hi *= 2.0
            if hi > 1e8:
                return np.inf
        return float(brentq(fun, 1.0, hi, xtol=1e-14, rtol=1e-14))

    return _vec(one, _F(delta))


def weak_wave_pressure(M, delta, gamma: float = 1.4):
    """Pressure change across a weak wave that turns a supersonic stream by a small angle:
    (p₂ − p₁)/p₁ = γM²δ/√(M² − 1).

    **ARGUMENT ORDER: (M, delta)** — the siblings ``ch15.ackeret_pressure(delta, M_inf)`` and
    ``ch15.ackeret_cp(delta, M_inf)`` take the angle FIRST.  Call by keyword when in doubt.

    Book: §15.9, Eq. (15.55) (and (15.57) with M = M_∞).  **Sign**: δ > 0 when the wall turns *into* the flow
    (compression); δ < 0 for a turn away (expansion) — the same formula serves both because the entropy change is
    third order in the pressure change (15.44b).
    Parameters: M [–] > 1; delta [rad], |δ| ≪ 1; gamma [–].   Returns (p₂ − p₁)/p₁ [–] (NaN for M ≤ 1).
    Assumptions: small deflection; fails as M → 1 (the coefficient diverges) and near δ_max (trap T16).
    Validation: V2, V5, V1 — tests/test_ch15.py: test_weak_wave_rule_V2_derivation,
    test_weak_wave_V1_exact_shock_and_expansion_approach_the_linear_rule_from_both_sides.
    Label: symbolic, benchmark, analytic.
    """
    M = _F(M)
    with np.errstate(divide="ignore", invalid="ignore"):
        return _S(np.where(M > 1.0, gamma * M ** 2 * _F(delta) / np.sqrt(np.where(M > 1.0, M ** 2 - 1.0, np.nan)),
                           np.nan))  # Eq. (15.55)


def weak_wave_error(M: float, delta: float, gamma: float = 1.4) -> dict:
    """How good is the weak-wave rule?  The linear (15.55) beside the exact oblique shock (δ > 0) or the exact
    Prandtl–Meyer expansion (δ < 0) (ours).

    Book: §15.9, Eq. (15.55) and the remark that it holds for weak expansions too.
    Parameters: M [–] > 1; delta [rad] (book sign: + into the flow); gamma [–].
    Returns dict [–]: ``linear``, ``exact`` ((p₂ − p₁)/p₁), ``ratio`` (exact/linear), ``rel_error``
    ((linear − exact)/exact), ``kind`` ("shock", "expansion" or "none").  ``exact`` is NaN for a detached shock.
    Assumptions: as the two exact functions.
    Validation: V5, V1 — tests/test_ch15.py:
    test_weak_wave_V1_exact_shock_and_expansion_approach_the_linear_rule_from_both_sides.
    Label: benchmark, analytic.
    """
    lin = float(weak_wave_pressure(M, delta, gamma))
    if delta == 0.0:
        return dict(linear=0.0, exact=0.0, ratio=1.0, rel_error=0.0, kind="none")
    if delta > 0.0:
        if delta > float(max_deflection(M, gamma)["delta_max"]):
            ex = np.nan
        else:
            ex = float(oblique_shock(M, delta=delta, gamma=gamma)["p_ratio"]) - 1.0
        kind = "shock"
    else:
        ex = float(prandtl_meyer_turn(M, -delta, gamma)["p_ratio"]) - 1.0
        kind = "expansion"
    return dict(linear=lin, exact=ex, ratio=ex / lin, rel_error=(lin - ex) / ex, kind=kind)


def prandtl_meyer(M, gamma: float = 1.4):
    """Prandtl–Meyer function ν(M) = √((γ + 1)/(γ − 1)) tan⁻¹√((γ − 1)(M² − 1)/(γ + 1)) − tan⁻¹√(M² − 1).

    Book: §15.9, Eq. (15.56) (the closed form of ∫₁^M √(M² − 1)/(1 + ½(γ − 1)M²) dM/M).  ν is the angle through which
    a sonic stream must be turned (isentropically, away from itself) to reach Mach number M.
    Parameters: M [–]; gamma [–].   Returns ν [rad] (0 at M = 1; **NaN for M < 1**).
    Assumptions: perfect gas; isentropic simple wave (waves of one family only).
    Validation: V1, V2, V5 — tests/test_ch15.py:
    test_prandtl_meyer_V1_closed_form_equals_the_quadrature_of_its_integrand,
    test_prandtl_meyer_differential_relation_V2_derivation, test_mach_angle_and_prandtl_meyer_V5_naca_report_1135_rows
    (also used by test_prandtl_meyer_turn_V4_isentropic_and_composes, whose conservation evidence belongs to
    :func:`prandtl_meyer_turn`).
    Label: analytic, symbolic, benchmark.
    """
    M = _F(M)
    k = (gamma + 1.0) / (gamma - 1.0)
    with np.errstate(invalid="ignore"):
        b = np.sqrt(np.where(M >= 1.0, M ** 2 - 1.0, np.nan))
        return _S(np.sqrt(k) * np.arctan(b / np.sqrt(k)) - np.arctan(b))  # Eq. (15.56)


def prandtl_meyer_max(gamma: float = 1.4) -> float:
    """Largest possible turning angle, ν(M → ∞) = (π/2)(√((γ + 1)/(γ − 1)) − 1) (ours).

    Book: §15.9, Eq. (15.56) in the limit M → ∞.
    Parameters: gamma [–].   Returns ν_max [rad] (130.45° for γ = 1.4: expansion to vacuum).
    Assumptions: as :func:`prandtl_meyer` (a continuum, perfect-gas idealisation).
    Validation: V1 — tests/test_ch15.py: test_prandtl_meyer_V1_closed_form_equals_the_quadrature_of_its_integrand.
    Label: analytic.
    """
    return float(0.5 * np.pi * (np.sqrt((gamma + 1.0) / (gamma - 1.0)) - 1.0))


def mach_from_prandtl_meyer(nu, gamma: float = 1.4):
    """Mach number that belongs to a Prandtl–Meyer angle (inverse of (15.56)).

    Book: §15.9, Eq. (15.56) inverted (the step "would determine M₂" after δ₁ + ν(M₁) = δ₂ + ν(M₂)).
    Parameters: nu [rad] in [0, ν_max); gamma [–].   Returns M [–] ≥ 1 (NaN outside the range).
    Method: ``brentq`` with an expanding bracket; residual asserted.
    Assumptions: as :func:`prandtl_meyer`.
    Validation: V1, V5 — tests/test_ch15.py: test_prandtl_meyer_V1_closed_form_equals_the_quadrature_of_its_integrand,
    test_mach_angle_and_prandtl_meyer_V5_naca_report_1135_rows.
    Label: analytic, benchmark.
    """
    numax = prandtl_meyer_max(gamma)

    def one(v: float) -> float:
        if not np.isfinite(v) or v < 0.0 or v >= numax:
            return np.nan
        if v == 0.0:
            return 1.0
        fun = lambda m: float(prandtl_meyer(m, gamma)) - v  # noqa: E731
        hi = 2.0
        while fun(hi) < 0.0:
            hi *= 2.0
            if hi > 1e12:
                return np.nan
        root = float(brentq(fun, 1.0, hi, xtol=1e-15, rtol=1e-15))
        assert abs(fun(root)) < 1e-10, "Prandtl–Meyer inversion did not converge"
        return root

    return _vec(one, _F(nu))


_NU_TOL = 1e-12  # [rad] round-off allowance below ν = 0 in prandtl_meyer_turn


def prandtl_meyer_turn(M1, turn, gamma: float = 1.4) -> dict:
    """Mach number and isentropic ratios after a supersonic stream is turned through an angle by a simple wave.

    Book: §15.9, δ₁ + ν(M₁) = δ₂ + ν(M₂) (unnumbered, after (15.56)); Fig. 15.23.
    **Sign convention**: ``turn`` > 0 is an **expansion** (the wall turns away from the flow), so ν(M₂) = ν(M₁) +
    turn.  The book's δ is positive when the wall turns *into* the flow: δ₂ − δ₁ = −turn (trap T14).  A negative
    ``turn`` is an isentropic (gradual) compression — valid only until the Mach lines cross.
    Parameters: M1 [–] ≥ 1; turn [rad]; gamma [–].
    Returns dict: ``M2`` [–], ``nu1``, ``nu2`` [rad], ``mu1``, ``mu2`` [rad], ``p_ratio`` (p₂/p₁), ``T_ratio``,
    ``rho_ratio`` [–], ``delta_book`` (= −turn) [rad].  NaN where ν(M₁) + turn is outside [0, ν_max).
    Round-off at the sonic end: a compression that brings the stream exactly back to M = 1 can leave ν(M₁) + turn
    a few units of round-off below zero (−10⁻¹⁷ rad instead of 0); a value in (−``nu_tol``, 0) with ``nu_tol`` =
    10⁻¹² rad (``_NU_TOL``, about 2 × 10⁻⁷ arc-seconds) is set to 0, so M₂ = 1 is returned and not NaN.  A turn
    that overshoots by more than that is still NaN.
    Assumptions: perfect gas; isentropic; waves of one family only.
    Validation: V1, V5, V4 — entropy and total enthalpy unchanged; two half-turns compose; a compression undoes the turn
    — tests/test_ch15.py: test_weak_wave_V1_exact_shock_and_expansion_approach_the_linear_rule_from_both_sides,
    test_prandtl_meyer_turn_V4_isentropic_and_composes; V7 (the return to M = 1 and the round-off allowance) —
    test_prandtl_meyer_turn_V7_return_to_mach_one_and_the_round_off_allowance.
    Label: analytic, benchmark, conserved.
    """
    M1 = _F(M1)
    nu1 = _F(prandtl_meyer(M1, gamma))
    nu2 = nu1 + _F(turn)
    with np.errstate(invalid="ignore"):
        nu2 = np.where((nu2 < 0.0) & (nu2 > -_NU_TOL), 0.0, nu2)  # back to M = 1 exactly, to round-off
    M2 = _F(mach_from_prandtl_meyer(nu2, gamma))
    t = _F(stagnation_temperature_ratio(M1, gamma)) / _F(stagnation_temperature_ratio(M2, gamma))  # T₂/T₁, T₀ const
    out = dict(M2=M2, nu1=nu1, nu2=np.where(np.isfinite(M2), nu2, np.nan), mu1=_F(mach_angle(M1)),
               mu2=_F(mach_angle(M2)), p_ratio=t ** (gamma / (gamma - 1.0)), T_ratio=t,
               rho_ratio=t ** (1.0 / (gamma - 1.0)), delta_book=-_F(turn) + 0.0 * M1)
    return {k: _S(v) for k, v in out.items()}


# ======================================================================================================================
# optional (ours): the exact Riemann solver — a moving shock, a contact and a fan from one initial jump
# ======================================================================================================================
def riemann_exact(left, right, gamma: float = 1.4, xi=0.0) -> dict:
    """Exact solution of the one-dimensional Riemann problem for a perfect gas (**not in the book — ours**; optional
    companion of the moving shock of §15.6 and of the steepening argument of §15.2).

    Two uniform states meet at x = 0 at t = 0; for t > 0 the solution depends on ξ = x/t only: a shock or a
    rarefaction fan on each side of a contact surface.
    Book: none (the building blocks are (15.39)–(15.41) for the shock and the isentropic relations (15.2) with the
    simple-wave relation u ± 2c/(γ − 1) = const for the fan).
    Parameters: left, right — (ρ [kg/m³], u [m/s], p [Pa]) of the two states; gamma [–]; xi = x/t [m/s] (scalar or
    array).
    Returns dict: ``rho``, ``u``, ``p`` at ``xi``; ``p_star`` [Pa], ``u_star`` [m/s] (between the waves);
    ``rho_star_left``, ``rho_star_right``; ``left_wave``, ``right_wave`` ("shock" or "rarefaction");
    ``shock_speed_left``, ``shock_speed_right`` [m/s] (NaN for a rarefaction).
    Method: the standard pressure function f(p) = f_L + f_R + u_R − u_L (shock branch from the jump relations,
    rarefaction branch from the isentrope) solved by ``brentq``; residual asserted.  Raises ``ValueError`` if the two
    states separate so fast that a vacuum forms.
    Assumptions: perfect gas with constant γ; inviscid; one space dimension.
    Validation: V4, V3 — tests/test_ch15.py: test_riemann_exact_V4_shocks_satisfy_the_jump_laws_and_fans_are_isentropic,
    test_euler1d_V3_second_order_for_a_smooth_wave_and_first_order_at_a_captured_shock.
    Label: conserved, converged.
    """
    rL, uL, pL = (float(v) for v in left)
    rR, uR, pR = (float(v) for v in right)
    g = gamma
    cL, cR = np.sqrt(g * pL / rL), np.sqrt(g * pR / rR)
    if 2.0 * (cL + cR) / (g - 1.0) <= uR - uL:
        raise ValueError("vacuum forms between the two states (they separate too fast)")

    def fK(p: float, rK: float, pK: float, cK: float) -> float:
        if p > pK:  # shock
            return (p - pK) * np.sqrt(2.0 / ((g + 1.0) * rK) / (p + (g - 1.0) / (g + 1.0) * pK))
        return 2.0 * cK / (g - 1.0) * ((p / pK) ** ((g - 1.0) / (2.0 * g)) - 1.0)  # rarefaction

    fun = lambda p: fK(p, rL, pL, cL) + fK(p, rR, pR, cR) + uR - uL  # noqa: E731
    lo, hi = 1e-12 * min(pL, pR), max(pL, pR)
    while fun(hi) < 0.0:
        hi *= 2.0
    ps = float(brentq(fun, lo, hi, xtol=1e-14 * hi, rtol=1e-14))
    assert abs(fun(ps)) < 1e-8 * (cL + cR), "Riemann pressure root did not converge"
    us = 0.5 * (uL + uR) + 0.5 * (fK(ps, rR, pR, cR) - fK(ps, rL, pL, cL))
    k = (g - 1.0) / (g + 1.0)
    xi = _F(xi)
    rho, u, p = np.empty_like(xi), np.empty_like(xi), np.empty_like(xi)
    # left side
    if ps > pL:
        rsL = rL * (ps / pL + k) / (k * ps / pL + 1.0)
        SL = uL - cL * np.sqrt((g + 1.0) / (2.0 * g) * ps / pL + (g - 1.0) / (2.0 * g))
        headL = tailL = SL
    else:
        rsL = rL * (ps / pL) ** (1.0 / g)
        SL = np.nan
        headL, tailL = uL - cL, us - cL * (ps / pL) ** ((g - 1.0) / (2.0 * g))
    if ps > pR:
        rsR = rR * (ps / pR + k) / (k * ps / pR + 1.0)
        SR = uR + cR * np.sqrt((g + 1.0) / (2.0 * g) * ps / pR + (g - 1.0) / (2.0 * g))
        headR = tailR = SR
    else:
        rsR = rR * (ps / pR) ** (1.0 / g)
        SR = np.nan
        headR, tailR = uR + cR, us + cR * (ps / pR) ** ((g - 1.0) / (2.0 * g))
    with np.errstate(invalid="ignore"):
        cfL = 2.0 / (g + 1.0) * (cL + 0.5 * (g - 1.0) * (uL - xi))  # inside the left fan
        cfR = 2.0 / (g + 1.0) * (cR - 0.5 * (g - 1.0) * (uR - xi))  # inside the right fan
        conds = [xi <= headL, (xi > headL) & (xi < tailL), (xi >= tailL) & (xi <= us), (xi > us) & (xi <= tailR),
                 (xi > tailR) & (xi < headR), xi >= headR]
        rho = np.select(conds, [rL, rL * (cfL / cL) ** (2.0 / (g - 1.0)), rsL, rsR,
                                rR * (cfR / cR) ** (2.0 / (g - 1.0)), rR])
        u = np.select(conds, [uL, 2.0 / (g + 1.0) * (cL + 0.5 * (g - 1.0) * uL + xi), us, us,
                              2.0 / (g + 1.0) * (-cR + 0.5 * (g - 1.0) * uR + xi), uR])
        p = np.select(conds, [pL, pL * (cfL / cL) ** (2.0 * g / (g - 1.0)), ps, ps,
                              pR * (cfR / cR) ** (2.0 * g / (g - 1.0)), pR])
    return dict(rho=_S(rho), u=_S(u), p=_S(p), p_star=ps, u_star=float(us), rho_star_left=float(rsL),
                rho_star_right=float(rsR), left_wave="shock" if ps > pL else "rarefaction",
                right_wave="shock" if ps > pR else "rarefaction", shock_speed_left=float(SL),
                shock_speed_right=float(SR))
