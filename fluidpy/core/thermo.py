"""Classical thermodynamics of a simple compressible substance and the perfect gas (reused by Ch. 4, 7, 13, 15).

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e, §1.8 (first law, equations of state, enthalpy, specific heats,
entropy, Gibbs relations, speed of sound, thermal expansion) and §1.9 (perfect gas).

Conventions (see ``knowledge/notation.md``)
-------------------------------------------
* Everything per unit mass (lower-case e, h, q, w, s, v) in SI: J/kg, J/(kg K), m^3/kg. Temperatures in kelvin.
* First-law sign: ``q`` is heat added **to** the system and ``w`` is work done **on** the system, so a reversible
  expansion does work ``-p dv`` on the particle (Eq. 1.11).
* Gas constants: ``R_U`` is per **kilomole** (J kmol^-1 K^-1) as in the book; ``M_w`` is in kg/kmol; ``R = R_U/M_w``
  is per kilogram. Mixing mol and kmol is a factor-1000 trap.

Constants and their sources
---------------------------
* ``K_B``, ``N_A``: exact values of the 2019 SI redefinition (CODATA; NIST, https://physics.nist.gov/cuu/Constants/).
* ``M_W_AIR`` = 28.9644 kg/kmol: U.S. Standard Atmosphere 1976 (NASA-TM-X-74335), sea-level mean molecular weight.
* ``GAMMA_AIR`` = 1.4: USSA-1976 convention for air; ``CP_AIR``, ``CV_AIR`` are *derived* from it and ``R_AIR``
  through (1.23)–(1.24), not quoted.
* ``G0`` = 9.80665 m/s^2 (standard gravity, exact by definition, 3rd CGPM 1901); ``P_ATM`` = 101325 Pa (standard
  atmosphere, exact by definition). ``P_REF`` = 1.0e5 Pa is the reference pressure for potential temperature
  (meteorological convention; the book says "nearly equal to 100 kPa").
* Molar masses in ``MOLAR_MASS``: IUPAC standard atomic weights (N2, O2, Ar as tabulated in USSA-1976).
* ``VDW_CO2_MOLAR``: CRC Handbook of Chemistry and Physics (via the Wikipedia van der Waals data page).
* Tait exponent n = 7.15 for water: R. H. Cole, *Underwater Explosions*, Princeton (1948), modified Tait equation.
"""
from __future__ import annotations

from typing import Callable, Mapping

import numpy as np

from ._util import as_scalar_if_0d, require_positive

# --------------------------------------------------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------------------------------------------------
K_B: float = 1.380649e-23          #: Boltzmann constant [J/K] (exact, SI 2019)
N_A: float = 6.02214076e23         #: Avogadro constant [1/mol] (exact, SI 2019)
N_A_KMOL: float = N_A * 1.0e3      #: Avogadro's number per kilomole [1/kmol] (the book's A_o)
R_U: float = K_B * N_A_KMOL        #: universal gas constant [J/(kmol K)] = k_B A_o (§1.9), 8314.462618...
M_W_AIR: float = 28.9644           #: mean molecular weight of dry air [kg/kmol] (USSA-1976)
R_AIR: float = R_U / M_W_AIR       #: gas constant of dry air [J/(kg K)], ≈ 287.05
GAMMA_AIR: float = 1.4             #: ratio of specific heats of air (USSA-1976 convention)
CP_AIR: float = GAMMA_AIR * R_AIR / (GAMMA_AIR - 1.0)   #: derived C_p of air [J/(kg K)], ≈ 1004.7
CV_AIR: float = R_AIR / (GAMMA_AIR - 1.0)               #: derived C_v of air [J/(kg K)], ≈ 717.6
G0: float = 9.80665                #: standard gravitational acceleration [m/s^2] (exact)
G_BOOK: float = 9.81               #: the book's rounded g, the default of the Ch. 4 modules (its examples use 9.81)
P_ATM: float = 101325.0            #: standard atmospheric pressure [Pa] (exact)
P_REF: float = 1.0e5               #: reference pressure p_o for potential temperature/density [Pa]

#: Molar masses [kg/kmol] (IUPAC standard atomic weights; air: USSA-1976).
MOLAR_MASS: dict[str, float] = {
    "He": 4.002602,
    "H2O": 18.01528,
    "N2": 28.0134,
    "O2": 31.9988,
    "Ar": 39.948,
    "CO2": 44.0095,
    "air": M_W_AIR,
}

#: Ideal ratio of specific heats from equipartition for rigid molecules (kinetic theory): 3 translational degrees of
#: freedom (monatomic, 5/3), + 2 rotational (linear/diatomic, 7/5), + 3 rotational (non-linear polyatomic, 4/3).
GAMMA_BY_ATOMICITY: dict[str, float] = {"monatomic": 5.0 / 3.0, "diatomic": 7.0 / 5.0, "triatomic": 4.0 / 3.0}
GAMMA_IDEAL = GAMMA_BY_ATOMICITY  # alias

#: van der Waals constants of CO2 per mole: a [Pa m^6/mol^2], b [m^3/mol] (CRC Handbook: 3.640 L^2 bar/mol^2,
#: 0.04267 L/mol).
VDW_CO2_MOLAR: dict[str, float] = {"a": 0.3640, "b": 4.267e-5}
#: The same constants per unit mass: a [Pa m^6/kg^2] = a_molar/M^2, b [m^3/kg] = b_molar/M with M = 44.0095e-3 kg/mol.
VDW_CO2: dict[str, float] = {"a": VDW_CO2_MOLAR["a"] / (MOLAR_MASS["CO2"] * 1e-3) ** 2,
                             "b": VDW_CO2_MOLAR["b"] / (MOLAR_MASS["CO2"] * 1e-3)}

#: Tait exponent for water (Cole 1948).
TAIT_N_WATER: float = 7.15


# --------------------------------------------------------------------------------------------------------------------
# Perfect gas: molecular and continuum forms
# --------------------------------------------------------------------------------------------------------------------
def molecular_gas_pressure(n, V, T):
    """Pressure of n non-interacting molecules in a volume V at temperature T (molecular perfect-gas law).

    Book: §1.9, Eq. (1.21) ``pV = n k_B T``.

    Parameters
    ----------
    n : float or array_like
        Number of molecules [-] (a count, not moles).
    V : float or array_like
        Container volume [m^3], > 0.
    T : float or array_like
        Absolute temperature [K].

    Returns
    -------
    p : float or ndarray
        Average pressure on the container walls [Pa].

    Notes
    -----
    Assumptions: attractive forces negligible; V/n much larger than the volume of one molecule.

    Validation: V1 n k_B T/V (rel 1e-14) and the D11 chain: equals :func:`perfect_gas_pressure` of the same gas
    (rel 1e-12); it is the reference n k_B T that sampled Maxwellian molecules reproduce within 1 % in the kinetic
    pressure test. Label: analytic.
    """
    require_positive("V", V)
    p = np.asarray(n, dtype=float) * K_B * np.asarray(T, dtype=float) / np.asarray(V, dtype=float)  # Eq. (1.21)
    return as_scalar_if_0d(p)


def gas_constant(M_w):
    """Specific gas constant of a gas with molecular weight M_w.

    Book: §1.9, Eq. (1.22) (``R = R_u / M_w`` with ``R_u = k_B A_o``).

    Parameters
    ----------
    M_w : float or array_like
        Molecular weight [kg/kmol] (e.g. 28.9644 for air). **Per kilomole**, as in the book.

    Returns
    -------
    R : float or ndarray
        Gas constant [J/(kg K)].

    Notes
    -----
    Assumptions: perfect gas; M_w of a mixture is its mole-fraction-weighted mean.

    Validation: V1 ``gas_constant(M_W_AIR)`` = R_AIR (rel 1e-15), R_U = k_B A_o; V5 R_AIR vs USSA-1976 R*/M0
    (rel 1.7e-5, tolerance 5e-5) and p/(rho T) of every USSA Table 1 row (rel 2e-4); CODATA k_B, N_A exact.
    Label: analytic, benchmark.
    """
    require_positive("M_w", M_w)
    return as_scalar_if_0d(R_U / np.asarray(M_w, dtype=float))  # Eq. (1.22): R = R_u / M_w


def molecule_mass(M_w):
    """Mass of one molecule, ``m = M_w / A_o``.

    Book: §1.9 (text before Eq. (1.22): m = M_w/A_o, A_o = Avogadro's number per kg-mole).

    Parameters
    ----------
    M_w : float or array_like
        Molecular weight [kg/kmol].

    Returns
    -------
    m : float or ndarray
        [kg] (≈ 4.81e-26 kg for air).

    Validation: V1 k_B/m = R for air (rel 1e-12, D11 step 4); V5 through :func:`mean_molecular_speed`, which matches
    the USSA-1976 Table 2 mean particle speed (rel 1.9e-5). Label: analytic, benchmark.
    """
    return as_scalar_if_0d(np.asarray(M_w, dtype=float) / N_A_KMOL)  # §1.9: m = M_w / A_o


molecular_mass = molecule_mass  #: alias (design Part C name)


def perfect_gas_pressure(rho, T, R=R_AIR):
    """Pressure of a perfect gas from its density and temperature.

    Book: §1.9, Eq. (1.22) ``p = rho R T``.

    Parameters
    ----------
    rho : float or array_like
        Density [kg/m^3].
    T : float or array_like
        Absolute temperature [K].
    R : float, optional
        Specific gas constant [J/(kg K)]; default dry air ``R_AIR``.

    Returns
    -------
    p : float or ndarray
        Absolute pressure [Pa].

    Notes
    -----
    Assumptions: perfect gas (molecular interactions negligible), continuum (Kn << 1).

    Validation: V1 D11 chain equals :func:`molecular_gas_pressure` (rel 1e-12); V2 pint rho R T is [pressure] and
    matches the code (rel 1e-14); a = b = 0 van der Waals reduces to it (rel 1e-14). Label: analytic, symbolic.
    """
    return as_scalar_if_0d(np.asarray(rho, dtype=float) * R * np.asarray(T, dtype=float))  # Eq. (1.22)


def perfect_gas_density(p, T, R=R_AIR):
    """Density of a perfect gas from pressure and temperature.

    Book: §1.9, Eq. (1.22) rearranged: ``rho = p / (R T)``.

    Parameters
    ----------
    p : float or array_like
        Absolute pressure [Pa].
    T : float or array_like
        Absolute temperature [K], > 0.
    R : float, optional
        Specific gas constant [J/(kg K)]; default ``R_AIR``.

    Returns
    -------
    rho : float or ndarray
        Density [kg/m^3].

    Notes
    -----
    Assumptions: perfect gas.

    Validation: V5 USSA-1976 sea level ``perfect_gas_density(101325, 288.15)`` vs Table 1 1.2250 kg/m^3 (rel −1.8e-5,
    tolerance 1e-4); V1 worked number 101325/(287.058 · 288.15) (rel 1e-5); round trip through
    :func:`perfect_gas_state`. Label: benchmark, analytic.
    """
    require_positive("T", T)
    return as_scalar_if_0d(np.asarray(p, dtype=float) / (R * np.asarray(T, dtype=float)))  # Eq. (1.22)


def perfect_gas_state(p=None, rho=None, T=None, R=R_AIR):
    """Fill in the missing one of (p, rho, T) for a perfect gas: two properties fix the state.

    Book: §1.8, Eq. (1.12) (thermal equation of state ``p = p(v, T)``), specialised with §1.9, Eq. (1.22).

    Parameters
    ----------
    p : float or array_like or None
        Absolute pressure [Pa].
    rho : float or array_like or None
        Density [kg/m^3].
    T : float or array_like or None
        Absolute temperature [K].
    R : float, optional
        Specific gas constant [J/(kg K)].

    Returns
    -------
    (p, rho, T) : tuple
        The complete state [Pa, kg/m^3, K].

    Raises
    ------
    ValueError
        Unless exactly two of the three are given.

    Notes
    -----
    Assumptions: single-component perfect gas (a simple compressible substance; seawater would need salinity too).

    Validation: V1 any two reproduce the third (rel 1e-14); V7 none, one or three arguments raise. Label: analytic.
    """
    given = [x is not None for x in (p, rho, T)]
    if sum(given) != 2:
        raise ValueError("give exactly two of p, rho, T (two properties fix the state of a simple substance)")
    if p is None:
        p = perfect_gas_pressure(rho, T, R)  # Eq. (1.22)
    elif rho is None:
        rho = perfect_gas_density(p, T, R)  # Eq. (1.22)
    else:
        require_positive("rho", rho)
        T = as_scalar_if_0d(np.asarray(p, dtype=float) / (R * np.asarray(rho, dtype=float)))  # Eq. (1.22)
    return p, rho, T


def specific_volume(rho):
    """Specific volume, the volume per unit mass.

    Book: §1.8, ``v = 1/rho`` (defined just before Eq. 1.11).

    Parameters
    ----------
    rho : float or array_like
        Density [kg/m^3], > 0.

    Returns
    -------
    v : float or ndarray
        Specific volume [m^3/kg].

    Validation: V1 rho = 4 gives 0.25 exactly. Label: analytic.
    """
    require_positive("rho", rho)
    return as_scalar_if_0d(1.0 / np.asarray(rho, dtype=float))  # §1.8: v = 1/rho


def enthalpy(e, p, v):
    """Specific enthalpy from internal energy, pressure and specific volume.

    Book: §1.8, Eq. (1.13) ``h ≡ e + p v``.

    Parameters
    ----------
    e : float or array_like
        Specific internal energy [J/kg].
    p : float or array_like
        Absolute pressure [Pa].
    v : float or array_like
        Specific volume [m^3/kg].

    Returns
    -------
    h : float or ndarray
        Specific enthalpy [J/kg].

    Notes
    -----
    Assumptions: none beyond the definition (any substance).

    Validation: V1 e + p v with perfect-gas e and v = R T/p equals :func:`perfect_gas_enthalpy` (rel 1e-12).
    Label: analytic.
    """
    return as_scalar_if_0d(np.asarray(e, dtype=float) + np.asarray(p, dtype=float) * np.asarray(v, dtype=float))  # Eq. (1.13)


def perfect_gas_internal_energy(T, cv=CV_AIR, T_ref=0.0, e_ref=0.0):
    """Specific internal energy of a perfect gas with constant C_v: ``e = e_ref + C_v (T − T_ref)``.

    Book: §1.9 (a perfect gas has e = e(T), claim before Eq. 1.25) with §1.8, Eq. (1.15) integrated at constant C_v.

    Parameters
    ----------
    T : float or array_like
        Absolute temperature [K].
    cv : float, optional
        Specific heat at constant volume [J/(kg K)], constant; default ``CV_AIR``.
    T_ref, e_ref : float, optional
        Reference temperature [K] and the energy assigned to it [J/kg] (only differences of e are physical).

    Returns
    -------
    e : float or ndarray
        Specific internal energy [J/kg].

    Notes
    -----
    Assumptions: perfect gas, constant specific heats (the book notes C_v actually rises with T).

    Validation: V1 h − e = R T at 300 K with :func:`perfect_gas_enthalpy` (rel 1e-12). Label: analytic.
    """
    return as_scalar_if_0d(e_ref + cv * (np.asarray(T, dtype=float) - T_ref))  # e = e(T), de = C_v dT


def perfect_gas_enthalpy(T, cp=CP_AIR, T_ref=0.0, h_ref=0.0):
    """Specific enthalpy of a perfect gas with constant C_p: ``h = h_ref + C_p (T − T_ref)``.

    Book: §1.9 (h = h(T)) with §1.8, Eq. (1.14) integrated at constant C_p.

    Parameters
    ----------
    T : float or array_like
        Absolute temperature [K].
    cp : float, optional
        Specific heat at constant pressure [J/(kg K)], constant; default ``CP_AIR``.
    T_ref, h_ref : float, optional
        Reference temperature [K] and enthalpy [J/kg].

    Returns
    -------
    h : float or ndarray
        Specific enthalpy [J/kg].

    Notes
    -----
    Assumptions: perfect gas, constant C_p. With ``T_ref = 0`` and zero references, h − e = (C_p − C_v) T = R T.

    Validation: V1 :func:`specific_heat_cp` of it returns C_p (rel 1e-6); h − e = R T (rel 1e-12). Label: analytic.
    """
    return as_scalar_if_0d(h_ref + cp * (np.asarray(T, dtype=float) - T_ref))  # h = h(T), dh = C_p dT


# --------------------------------------------------------------------------------------------------------------------
# Partial derivatives with a variable held fixed
# --------------------------------------------------------------------------------------------------------------------
def partial_derivative(f: Callable[..., float], var: str, at: Mapping[str, float], rel_step: float = 1e-6,
                       step: float | None = None):
    """Partial derivative of ``f`` with respect to one keyword argument, all others held fixed (central difference).

    Book: §1.8, the ``( )_p`` / ``( )_v`` / ``( )_s`` notation of Eqs. (1.14), (1.15), (1.19), (1.20).

    Parameters
    ----------
    f : callable
        Function of keyword arguments, e.g. ``f(T=..., p=...)``; must accept numpy arrays or floats.
    var : str
        Name of the argument to vary (e.g. ``"T"``).
    at : mapping
        Values of all arguments of ``f`` at the evaluation point (SI units of f's arguments).
    rel_step : float, optional
        Relative step: ``h = rel_step * max(|x|, 1)``. Default 1e-6 (truncation ~h^2, round-off ~eps/h).
    step : float, optional
        Absolute step h, overriding ``rel_step`` (for convergence studies).

    Returns
    -------
    dfdx : float or ndarray
        ``(f(x + h) − f(x − h)) / (2 h)`` [units of f per unit of var].

    Notes
    -----
    Method: second-order central difference (our choice; the book defines the derivatives only).
    Assumptions: f smooth near the point.

    Validation: V3 observed order 2.000 (2 ± 0.15) on ∂(y sin x)/∂x as ``step`` halves from 0.1 to 0.0125; V1 through
    C_p, C_v, α and c of known equations of state (rel 1e-6), and the van der Waals Gibbs identity (rel 1e-5).
    Label: analytic, converged.
    """
    x = np.asarray(at[var], dtype=float)
    h = step if step is not None else rel_step * np.maximum(np.abs(x), 1.0)
    plus = dict(at)
    minus = dict(at)
    plus[var] = x + h
    minus[var] = x - h
    return as_scalar_if_0d((np.asarray(f(**plus), dtype=float) - np.asarray(f(**minus), dtype=float)) / (2.0 * h))


def specific_heat_cp(h_of_Tp: Callable[[float, float], float], T, p, rel_step: float = 1e-6):
    """Specific heat at constant pressure from an enthalpy function h(T, p).

    Book: §1.8, Eq. (1.14) ``C_p ≡ (∂h/∂T)_p``.

    Parameters
    ----------
    h_of_Tp : callable
        ``h(T, p)`` [J/kg] with T [K] and p [Pa] as positional arguments.
    T : float
        Temperature [K].
    p : float
        Pressure [Pa] (held fixed).
    rel_step : float, optional
        Relative finite-difference step (see :func:`partial_derivative`).

    Returns
    -------
    cp : float
        [J/(kg K)].

    Notes
    -----
    Assumptions: single-component substance; h a smooth state function.

    Validation: V1 perfect-gas h = C_p T returns C_p (rel 1e-6). Label: analytic.
    """
    return partial_derivative(lambda T, p: h_of_Tp(T, p), "T", {"T": T, "p": p}, rel_step)  # Eq. (1.14)


def specific_heat_cv(e_of_Tv: Callable[[float, float], float], T, v, rel_step: float = 1e-6):
    """Specific heat at constant volume from an internal-energy function e(T, v).

    Book: §1.8, Eq. (1.15) ``C_v ≡ (∂e/∂T)_v``.

    Parameters
    ----------
    e_of_Tv : callable
        ``e(T, v)`` [J/kg] with T [K] and v [m^3/kg] as positional arguments.
    T : float
        Temperature [K].
    v : float
        Specific volume [m^3/kg] (held fixed).
    rel_step : float, optional
        Relative finite-difference step.

    Returns
    -------
    cv : float
        [J/(kg K)].

    Validation: V1 van der Waals e = C_v T − a/v (CO2) returns C_v (rel 1e-6). Label: analytic.
    """
    return partial_derivative(lambda T, v: e_of_Tv(T, v), "T", {"T": T, "v": v}, rel_step)  # Eq. (1.15)


def sound_speed_from_eos(p_of_rho_s: Callable[[float, float], float], rho, s=0.0, rel_step: float = 1e-6):
    """Speed of sound from an equation of state p(rho, s), differentiated at constant entropy.

    Book: §1.8, Eq. (1.19) ``c^2 = (∂p/∂rho)_s`` (proof deferred by the book to Ch. 15).

    Parameters
    ----------
    p_of_rho_s : callable
        ``p(rho, s)`` [Pa] with rho [kg/m^3] and s [J/(kg K)] positional.
    rho : float
        Density [kg/m^3].
    s : float, optional
        Specific entropy [J/(kg K)] (held fixed); ignored by EOS that are already isentropes.
    rel_step : float, optional
        Relative finite-difference step.

    Returns
    -------
    c : float
        Speed of sound [m/s].

    Raises
    ------
    ValueError
        If (∂p/∂rho)_s <= 0 (no real sound speed).

    Notes
    -----
    Assumptions: infinitesimal isentropic disturbances. Incompressible limit: ∂rho/∂p → 0 gives c → ∞.

    Validation: V1 p = K rho^γ gives sqrt(γ p/rho) at three densities (rel 1e-8) and 340.29 m/s at sea level
    (rel 1e-4); on the perfect-gas isentrope it equals :func:`perfect_gas_sound_speed` (rel 1e-8); V7 Tait liquid
    c = sqrt(K0/rho0) and c ∝ sqrt(K0) (incompressible limit c → ∞); a non-positive derivative raises.
    Label: analytic.
    """
    dpdrho = partial_derivative(lambda rho, s: p_of_rho_s(rho, s), "rho", {"rho": rho, "s": s}, rel_step)
    if np.any(np.asarray(dpdrho) <= 0):
        raise ValueError(f"(dp/drho)_s = {dpdrho!r} <= 0: no real speed of sound")
    return as_scalar_if_0d(np.sqrt(dpdrho))  # Eq. (1.19): c^2 = (dp/drho)_s


def thermal_expansion_coefficient(rho_of_Tp: Callable[[float, float], float], T, p, rel_step: float = 1e-6):
    """Thermal expansion coefficient from a density function rho(T, p).

    Book: §1.8, Eq. (1.20) ``alpha ≡ −(1/rho) (∂rho/∂T)_p``.

    Parameters
    ----------
    rho_of_Tp : callable
        ``rho(T, p)`` [kg/m^3] with T [K] and p [Pa] positional.
    T : float
        Temperature [K].
    p : float
        Pressure [Pa] (held fixed).
    rel_step : float, optional
        Relative finite-difference step.

    Returns
    -------
    alpha : float
        [1/K]; positive when the fluid expands on heating (negative for water below ~4 °C).

    Validation: V1 perfect gas returns 1/T at 200, 300, 400 K (rel 1e-6); linear EOS rho0(1 − a0 (T − T0)) returns
    a0/(1 − a0 ΔT) (rel 1e-6); water (Kell) α < 0 at 2 °C and > 0 at 8 °C. Label: analytic.
    """
    rho = np.asarray(rho_of_Tp(T, p), dtype=float)
    drho_dT = partial_derivative(lambda T, p: rho_of_Tp(T, p), "T", {"T": T, "p": p}, rel_step)
    return as_scalar_if_0d(-drho_dT / rho)  # Eq. (1.20)


# --------------------------------------------------------------------------------------------------------------------
# Specific heats of a perfect gas
# --------------------------------------------------------------------------------------------------------------------
def cv_from_cp(cp, R=R_AIR):
    """C_v of a perfect gas from C_p.

    Book: §1.9, Eq. (1.23) ``R = C_p − C_v``.

    Parameters
    ----------
    cp : float or array_like
        Specific heat at constant pressure [J/(kg K)].
    R : float, optional
        Gas constant [J/(kg K)].

    Returns
    -------
    cv : float or ndarray
        [J/(kg K)].

    Notes
    -----
    Assumptions: perfect gas (uses e = e(T), h = h(T)).

    Validation: V1 ``cv_from_cp(CP_AIR)`` = CV_AIR (rel 1e-13), C_p − C_v = R. Label: analytic.
    """
    return as_scalar_if_0d(np.asarray(cp, dtype=float) - R)  # Eq. (1.23)


def gamma_from_cp(cp, R=R_AIR):
    """Ratio of specific heats of a perfect gas from C_p and R.

    Book: §1.9, Eq. (1.24) ``gamma ≡ C_p/C_v`` with (1.23).

    Parameters
    ----------
    cp : float or array_like
        [J/(kg K)], > R.
    R : float, optional
        Gas constant [J/(kg K)].

    Returns
    -------
    gamma : float or ndarray
        [-].

    Validation: V1 round trip ``gamma_from_cp(cp_from_gamma(5/3))`` (rel 1e-13); C_p = 1005 gives 1.40 ± 0.002.
    Label: analytic.
    """
    cp = np.asarray(cp, dtype=float)
    return as_scalar_if_0d(cp / (cp - R))  # Eq. (1.24) with (1.23)


def cp_from_gamma(gamma, R=R_AIR):
    """C_p of a perfect gas from gamma and R: ``C_p = gamma R / (gamma − 1)``.

    Book: §1.9, Eqs. (1.23)–(1.24) solved for C_p.

    Parameters
    ----------
    gamma : float or array_like
        Ratio of specific heats [-], > 1.
    R : float, optional
        Gas constant [J/(kg K)].

    Returns
    -------
    cp : float or ndarray
        [J/(kg K)].

    Validation: V1 round trip with :func:`gamma_from_cp` (rel 1e-13); ``cp_from_gamma(1.4) − cv_from_gamma(1.4)`` = R.
    Label: analytic.
    """
    g = np.asarray(gamma, dtype=float)
    return as_scalar_if_0d(g * R / (g - 1.0))  # (1.23)+(1.24)


def cv_from_gamma(gamma, R=R_AIR):
    """C_v of a perfect gas from gamma and R: ``C_v = R / (gamma − 1)``.

    Book: §1.9, Eqs. (1.23)–(1.24) solved for C_v.

    Parameters
    ----------
    gamma : float or array_like
        [-], > 1.
    R : float, optional
        Gas constant [J/(kg K)].

    Returns
    -------
    cv : float or ndarray
        [J/(kg K)].

    Validation: V1 ``cp_from_gamma(1.4) − cv_from_gamma(1.4)`` = R_AIR (rel 1e-13). Label: analytic.
    """
    g = np.asarray(gamma, dtype=float)
    return as_scalar_if_0d(R / (g - 1.0))  # (1.23)+(1.24)


# --------------------------------------------------------------------------------------------------------------------
# Isentropic perfect gas
# --------------------------------------------------------------------------------------------------------------------
def isentropic_pressure(rho, p0, rho0, gamma=GAMMA_AIR):
    """Pressure along an isentrope of a perfect gas with constant specific heats.

    Book: §1.9, Eq. (1.25) ``p/rho^gamma = const``, written as ``p = p0 (rho/rho0)^gamma``.

    Parameters
    ----------
    rho : float or array_like
        Density [kg/m^3].
    p0, rho0 : float
        Pressure [Pa] and density [kg/m^3] of a reference state on the same isentrope.
    gamma : float, optional
        Ratio of specific heats [-].

    Returns
    -------
    p : float or ndarray
        [Pa].

    Notes
    -----
    Assumptions: perfect gas, constant C_p and C_v, adiabatic and frictionless (isentropic).

    Validation: V2 sympy D14 step by step (C_v dT = −p dv, dp/p = −γ dv/v, p v^γ constant); V3 a from-scratch Euler
    march of dp/drho = γ p/rho converges to it at order 0.998; V1 doubling rho gives p0 2^1.4 (rel 1e-14), consistent
    with :func:`isentropic_ratios` (bicycle pump, rel 1e-12). Label: symbolic, converged, analytic.
    """
    return as_scalar_if_0d(p0 * (np.asarray(rho, dtype=float) / rho0) ** gamma)  # Eq. (1.25)


def isentropic_ratios(p_over_p0, gamma=GAMMA_AIR):
    """Temperature and density ratios along an isentrope, as functions of the pressure ratio.

    Book: §1.9, Eq. (1.26) ``T/T0 = (p/p0)^((gamma−1)/gamma)`` and ``rho/rho0 = (p/p0)^(1/gamma)``.

    Parameters
    ----------
    p_over_p0 : float or array_like
        Pressure ratio p/p0 [-], > 0.
    gamma : float, optional
        Ratio of specific heats [-].

    Returns
    -------
    (T_over_T0, rho_over_rho0) : tuple of float or ndarray
        [-].

    Notes
    -----
    Assumptions: as :func:`isentropic_pressure`. gamma → 1 gives T/T0 → 1 (isothermal limit).

    Validation: V1 (T/T0)/(rho/rho0)^(γ−1) = 1 and p/p0 = (rho/rho0)(T/T0) (1e-14); bicycle pump p ratio 2 → 351.3 K,
    rho ratio 1.641; V7 γ → 1 gives T/T0 → 1 (1e-8); p ratio 0 raises. Label: analytic.
    """
    r = np.asarray(p_over_p0, dtype=float)
    require_positive("p_over_p0", r)
    T_ratio = r ** ((gamma - 1.0) / gamma)  # Eq. (1.26), temperature
    rho_ratio = r ** (1.0 / gamma)  # Eq. (1.26), density
    return as_scalar_if_0d(T_ratio), as_scalar_if_0d(rho_ratio)


def perfect_gas_sound_speed(T, gamma=GAMMA_AIR, R=R_AIR):
    """Speed of sound in a perfect gas.

    Book: §1.9, Eq. (1.27) ``c = sqrt(gamma R T)`` (from (1.19) with (1.25) and (1.22)).

    Parameters
    ----------
    T : float or array_like
        Absolute temperature [K].
    gamma : float, optional
        Ratio of specific heats [-].
    R : float, optional
        Gas constant [J/(kg K)].

    Returns
    -------
    c : float or ndarray
        [m/s].

    Notes
    -----
    Assumptions: perfect gas with constant gamma; isentropic (not isothermal — Newton's sqrt(R T) error).

    Validation: V5 USSA-1976 Table 1 sound-speed column at 0–50 km (rel 2e-4; max 2.1e-5); V1 equals
    :func:`sound_speed_from_eos` on the isentrope through sea level (rel 1e-8); c/sqrt(R T) = sqrt(γ) (not Newton's
    isothermal value). Label: benchmark, analytic.
    """
    require_positive("T", T)
    return as_scalar_if_0d(np.sqrt(gamma * R * np.asarray(T, dtype=float)))  # Eq. (1.27)


def perfect_gas_expansion_coefficient(T):
    """Thermal expansion coefficient of a perfect gas.

    Book: §1.9, Eq. (1.28) ``alpha = 1/T``.

    Parameters
    ----------
    T : float or array_like
        Absolute temperature [K], > 0.

    Returns
    -------
    alpha : float or ndarray
        [1/K].

    Notes
    -----
    Assumptions: perfect gas (rho = p/RT at fixed p).

    Validation: V1 equals :func:`thermal_expansion_coefficient` of the perfect-gas density at 200, 300, 400 K
    (rel 1e-6). Label: analytic.
    """
    require_positive("T", T)
    return as_scalar_if_0d(1.0 / np.asarray(T, dtype=float))  # Eq. (1.28)


# --------------------------------------------------------------------------------------------------------------------
# Processes: first law, path functions, entropy
# --------------------------------------------------------------------------------------------------------------------
_SINGLE_KINDS = ("isothermal", "isochoric", "isobaric", "isentropic")
_TWO_LEG_KINDS = ("isochoric-isobaric", "isobaric-isochoric", "isentropic-isochoric", "isothermal-isochoric")


def _check_leg(kind: str, v1: float, T1: float, v2: float, T2: float, gamma: float, rtol: float = 1e-9) -> None:
    """Raise ValueError if (v2, T2) cannot be reached from (v1, T1) by a single leg of this kind."""
    if kind == "isothermal":
        ok = np.isclose(T1, T2, rtol=rtol)
        need = "T2 == T1"
    elif kind == "isochoric":
        ok = np.isclose(v1, v2, rtol=rtol)
        need = "v2 == v1"
    elif kind == "isobaric":
        ok = np.isclose(T1 / v1, T2 / v2, rtol=rtol)
        need = "T2/v2 == T1/v1 (same p)"
    elif kind == "isentropic":
        ok = np.isclose(T1 * v1 ** (gamma - 1.0), T2 * v2 ** (gamma - 1.0), rtol=rtol)
        need = "T v^(gamma-1) unchanged"
    else:
        raise ValueError(f"unknown process kind {kind!r}")
    if not ok:
        raise ValueError(f"a single {kind} leg cannot join ({v1}, {T1}) to ({v2}, {T2}): needs {need}")


def _corner_state(kind: str, v1: float, T1: float, v2: float, T2: float, gamma: float) -> tuple[float, float]:
    """Intermediate state where the two legs of a two-leg path meet."""
    first, _second = kind.split("-")
    if kind == "isochoric-isobaric":
        return v1, T2 * v1 / v2  # v fixed at v1, then p fixed at p2 = R T2/v2
    if kind == "isobaric-isochoric":
        return v2, T1 * v2 / v1  # p fixed at p1 until v2, then v fixed
    if kind == "isentropic-isochoric":
        return v2, T1 * (v1 / v2) ** (gamma - 1.0)  # isentrope to v2, then v fixed
    if kind == "isothermal-isochoric":
        return v2, T1  # isotherm to v2, then v fixed
    raise ValueError(f"unknown two-leg kind {kind!r}")


def _leg_samples(kind: str, v1: float, T1: float, v2: float, T2: float, n: int, gamma: float):
    if kind == "isochoric":
        T = np.linspace(T1, T2, n)
        return np.full(n, float(v1)), T
    v = np.linspace(v1, v2, n)
    if kind == "isothermal":
        return v, np.full(n, float(T1))
    if kind == "isobaric":
        return v, T1 * v / v1  # p = R T/v constant
    return v, T1 * (v1 / v) ** (gamma - 1.0)  # isentropic: T v^(gamma-1) = const, from (1.25) with p v = R T


def process_path(kind: str, state1, state2, n: int = 401, gamma: float = GAMMA_AIR, R: float = R_AIR):
    """Sampled reversible perfect-gas process between two states, as one leg or two legs meeting at a corner.

    Book: §1.8 (heat and work are path functions, internal energy a state function; reversible processes,
    Eqs. (1.10)–(1.11)) and §1.9 (isentrope (1.25)). The path generator itself is ours.

    Parameters
    ----------
    kind : str
        Single legs ``"isothermal"`` (needs T1 = T2), ``"isochoric"`` (v1 = v2), ``"isobaric"`` (p1 = p2),
        ``"isentropic"`` (T v^(γ−1) equal); two legs ``"isochoric-isobaric"``, ``"isobaric-isochoric"``,
        ``"isentropic-isochoric"``, ``"isothermal-isochoric"`` (any two end states).
    state1, state2 : tuple (v, T)
        Specific volume [m^3/kg] and temperature [K] at the start and end.
    n : int, optional
        Samples per leg including both ends (default 401); two-leg paths have 2n − 1 samples.
    gamma : float, optional
        Ratio of specific heats [-].
    R : float, optional
        Gas constant [J/(kg K)] (for the pressure).

    Returns
    -------
    dict
        ``"v"`` [m^3/kg], ``"T"`` [K], ``"p"`` = R T/v [Pa] along the path; ``"corner"``: index of the corner sample
        (None for a single leg).

    Raises
    ------
    ValueError
        If a single-leg kind cannot join the two states.

    Notes
    -----
    Parameterisation: uniform in v (in T for isochoric legs). Assumptions: perfect gas, reversible, constant γ.

    Validation: V1 trapezoid totals over its samples reproduce the closed forms of :func:`path_heat_work_totals` for
    isothermal, isochoric and four two-leg paths (rel 1e-5 to 1e-6 at n = 4001); isochoric samples give p = R T/v
    exactly; corner index n − 1; an unreachable single leg raises. Leg-wise constancy of p (isobaric) and p v^γ
    (isentropic) is not tested separately. Label: analytic.
    """
    v1, T1 = map(float, state1)
    v2, T2 = map(float, state2)
    kind = kind.lower()
    if kind in _SINGLE_KINDS:
        _check_leg(kind, v1, T1, v2, T2, gamma)
        v, T = _leg_samples(kind, v1, T1, v2, T2, n, gamma)
        corner = None
    elif kind in _TWO_LEG_KINDS:
        first, second = kind.split("-")
        vc, Tc = _corner_state(kind, v1, T1, v2, T2, gamma)
        va, Ta = _leg_samples(first, v1, T1, vc, Tc, n, gamma)
        vb, Tb = _leg_samples(second, vc, Tc, v2, T2, n, gamma)
        v, T = np.concatenate((va, vb[1:])), np.concatenate((Ta, Tb[1:]))
        corner = n - 1
    else:
        raise ValueError(f"unknown process kind {kind!r}; use one of {_SINGLE_KINDS + _TWO_LEG_KINDS}")
    return {"v": v, "T": T, "p": R * T / v, "corner": corner}


def _leg_totals(kind: str, v1: float, T1: float, v2: float, T2: float, cv: float, R: float) -> dict[str, float]:
    """Exact heat, work, Δe and Δs for one reversible perfect-gas leg (constant C_v)."""
    cp = cv + R  # Eq. (1.23)
    de = cv * (T2 - T1)  # perfect gas: de = C_v dT
    if kind == "isothermal":
        w = -R * T1 * np.log(v2 / v1)  # −∫ p dv with p = R T/v
        ds = R * np.log(v2 / v1)
    elif kind == "isochoric":
        w = 0.0
        ds = cv * np.log(T2 / T1)
    elif kind == "isobaric":
        w = -R * (T2 - T1)  # −p (v2 − v1) with p v = R T
        ds = cp * np.log(T2 / T1)
    elif kind == "isentropic":
        w = de  # q = 0
        ds = 0.0
    else:
        raise ValueError(kind)
    q = de - w  # Eq. (1.10)
    return {"q": float(q), "w": float(w), "de": float(de), "ds": float(ds), "int_dq_over_T": float(ds)}


def path_heat_work_totals(kind: str, v1: float, T1: float, v2: float, T2: float, gamma: float = GAMMA_AIR,
                          R: float = R_AIR):
    """Exact totals of heat, work, internal-energy and entropy change along a reversible perfect-gas path.

    Book: §1.8, Eq. (1.10) ``δq + δw = Δe``, Eq. (1.11) ``de = dq − p dv``, Eq. (1.16) ``s2 − s1 = ∫ dq_rev/T`` and
    the Gibbs relation (1.18), with the perfect gas (1.22)–(1.24) and constant specific heats. Closed forms leg by leg:
    isothermal w = −R T ln(v2/v1); isochoric w = 0; isobaric w = −R ΔT; isentropic q = 0.

    Parameters
    ----------
    kind : str
        As :func:`process_path` (single or two-leg kinds).
    v1, T1 : float
        Start state [m^3/kg], [K]. Note: argument order differs from :func:`irreversible_process` (``kind, T1, v1``)
        and :func:`perfect_gas_entropy_change` (``T1, v1, T2, v2``) — here v comes before T. Pass keywords when in
        doubt.
    v2, T2 : float
        End state [m^3/kg], [K].
    gamma : float, optional
        Ratio of specific heats [-] (C_v = R/(γ − 1)).
    R : float, optional
        Gas constant [J/(kg K)].

    Returns
    -------
    dict of float
        ``q``, ``w`` (work done on the gas), ``de`` [J/kg]; ``ds``, ``int_dq_over_T`` [J/(kg K)] (equal for these
        reversible paths).

    Notes
    -----
    Scalar-callable closed forms for explainer parity. Assumptions: reversible, perfect gas, constant γ.

    Validation: V1 equals the trapezoid totals of :func:`process_heat_work` on n = 4001 sampled paths (rel 1e-5) and
    their ∫dq/T (rel 1e-5); isobaric w = −R ΔT, isochoric w = 0, isothermal q = −w; V4 Δe and Δs identical for four
    two-leg paths between the same states (ptp < 1e-9 J/kg and < 1e-12 J/(kg K)) while q and w differ by > 1e4 J/kg;
    Δs equals :func:`perfect_gas_entropy_change` (rel 1e-12); unknown kinds raise. Label: analytic, conserved.
    """
    cv = R / (gamma - 1.0)
    kind = kind.lower()
    if kind in _SINGLE_KINDS:
        _check_leg(kind, v1, T1, v2, T2, gamma)
        return _leg_totals(kind, v1, T1, v2, T2, cv, R)
    if kind in _TWO_LEG_KINDS:
        first, second = kind.split("-")
        vc, Tc = _corner_state(kind, v1, T1, v2, T2, gamma)
        a = _leg_totals(first, v1, T1, vc, Tc, cv, R)
        b = _leg_totals(second, vc, Tc, v2, T2, cv, R)
        return {k: a[k] + b[k] for k in a}
    raise ValueError(f"unknown process kind {kind!r}")


def irreversible_process(kind: str, T1: float, v1: float, cv: float = CV_AIR, R: float = R_AIR, T2: float | None = None,
                         v2: float | None = None):
    """Totals for two classic irreversible processes of a perfect gas: stirring at constant volume and free expansion.

    Book: §1.8 — the stirring example (temperature rises without heat transfer, so C_v cannot be *defined* as heat per
    degree) and the second law (ii), Clausius–Duhem, taught with the actual heat: ``s2 − s1 ≥ ∫ δq/T``. Δs is
    evaluated with the Gibbs relation (1.18), valid for any process between equilibrium states.

    Note: argument order differs from :func:`path_heat_work_totals` (``kind, v1, T1, v2, T2``) — here T1 comes
    before v1, as in :func:`perfect_gas_entropy_change`. Pass keywords when in doubt.

    Parameters
    ----------
    kind : {"stirring", "free_expansion"}
        ``"stirring"``: rigid insulated container, v fixed, stirring work raises T from T1 to T2 (needs T2 >= T1).
        ``"free_expansion"``: insulated gas expands into vacuum from v1 to v2 (needs v2 >= v1); T unchanged (e = e(T)).
    T1 : float
        Initial temperature [K].
    v1 : float
        Initial specific volume [m^3/kg].
    cv : float, optional
        Constant C_v [J/(kg K)].
    R : float, optional
        Gas constant [J/(kg K)].
    T2 : float, optional
        Final temperature [K] (stirring).
    v2 : float, optional
        Final specific volume [m^3/kg] (free expansion).

    Returns
    -------
    dict of float
        ``q`` (0), ``w`` (stirring work on the gas; 0 for free expansion), ``de`` [J/kg]; ``ds`` from (1.18) and
        ``int_dq_over_T`` (0) [J/(kg K)] — so ``ds > int_dq_over_T``.

    Validation: V1 stirring 300 → 330 K gives ds = C_v ln 1.1, free expansion v → 2v gives R ln 2 (rel 1e-14), both
    with ds > ∫δq/T = 0; equals :func:`free_expansion` and the end of :func:`stirred_isochoric_process`; cooling by
    stirring, compression "free expansion" and unknown kinds raise. Label: analytic.
    """
    kind = kind.lower()
    if kind == "stirring":
        if T2 is None or T2 < T1:
            raise ValueError("stirring needs T2 >= T1 (stirring work can only heat the gas)")
        de = cv * (T2 - T1)
        return {"q": 0.0, "w": float(de), "de": float(de), "ds": float(cv * np.log(T2 / T1)),  # (1.18) at dv = 0
                "int_dq_over_T": 0.0}
    if kind == "free_expansion":
        if v2 is None or v2 < v1:
            raise ValueError("free expansion needs v2 >= v1")
        return {"q": 0.0, "w": 0.0, "de": 0.0, "ds": float(R * np.log(v2 / v1)),  # (1.18) at dT = 0
                "int_dq_over_T": 0.0}
    raise ValueError("kind must be 'stirring' or 'free_expansion'")


def join_paths(*legs):
    """Concatenate (v, T) legs end to end, dropping the duplicated joint samples.

    Parameters
    ----------
    *legs : tuple of (ndarray, ndarray)
        ``(v, T)`` arrays (e.g. ``(d["v"], d["T"])`` from :func:`process_path`); each leg must start where the previous
        one ended.

    Returns
    -------
    (v, T) : tuple of ndarray
        The joined path [m^3/kg], [K].

    Validation: V1 joins the four legs of the closed test cycle (whose Δe, q + w and Δs close to < 1e-5); legs that do
    not meet raise. Label: analytic.
    """
    vs, Ts = [], []
    for i, (v, T) in enumerate(legs):
        v = np.asarray(v, dtype=float)
        T = np.asarray(T, dtype=float)
        if i > 0:
            if not (np.isclose(v[0], vs[-1][-1]) and np.isclose(T[0], Ts[-1][-1])):
                raise ValueError(f"leg {i} does not start where leg {i - 1} ended")
            v, T = v[1:], T[1:]
        vs.append(v)
        Ts.append(T)
    return np.concatenate(vs), np.concatenate(Ts)


def process_heat_work(v_path, T_path, cv=CV_AIR, R=R_AIR):
    """Cumulative work, internal-energy change and heat along a reversible perfect-gas path.

    Book: §1.8, Eq. (1.10) ``δq + δw = Δe`` and Eq. (1.11) ``de = dq − p dv`` (reversible, p dv work only), with the
    perfect gas (1.22) and constant C_v.

    Parameters
    ----------
    v_path : array_like
        Specific volume along the path [m^3/kg], length N >= 2.
    T_path : array_like
        Temperature along the path [K], length N.
    cv : float, optional
        Constant specific heat at constant volume [J/(kg K)]; default ``CV_AIR``.
    R : float, optional
        Gas constant [J/(kg K)].

    Returns
    -------
    dict
        ``"w"``: cumulative work done **on** the gas, ``−∫ p dv`` [J/kg] (trapezoid rule);
        ``"de"``: cumulative internal-energy change ``C_v (T − T_0)`` [J/kg];
        ``"q"``: cumulative heat added, ``de − w`` [J/kg] (Eq. 1.10);
        ``"p"``: pressure along the path ``R T / v`` [Pa]. All arrays start at 0 (p excepted).

    Notes
    -----
    Method: cumulative trapezoid (second order in the sample spacing) — refine ``n`` in :func:`process_path` for
    convergence studies. Assumptions: reversible (quasi-static, frictionless), perfect gas, constant C_v.

    Validation: V1 isothermal doubling w = −R T ln 2 (rel 1e-6 at n = 4001, −59.7 kJ/kg at 300 K), q = −w, Δe = 0;
    isochoric w = 0, q = C_v ΔT (rel 1e-14); V3 trapezoid work on an isentrope converges at order 1.999 (n = 11…81);
    V4 a closed reversible four-leg cycle returns Δe = 0 and q + w = 0 (1e-9) with net work > 1e3 J/kg; unequal or
    too-short inputs raise. Label: analytic, converged, conserved.
    """
    v = np.asarray(v_path, dtype=float)
    T = np.asarray(T_path, dtype=float)
    if v.shape != T.shape or v.ndim != 1 or v.size < 2:
        raise ValueError("v_path and T_path must be 1-D arrays of the same length >= 2")
    p = R * T / v  # Eq. (1.22) with v = 1/rho
    dw = -0.5 * (p[1:] + p[:-1]) * np.diff(v)  # reversible work on the particle, −p dv (trapezoid)
    w = np.concatenate(([0.0], np.cumsum(dw)))
    de = cv * (T - T[0])  # perfect gas: de = C_v dT
    q = de - w  # Eq. (1.10): δq + δw = Δe  →  q = Δe − w
    return {"q": q, "w": w, "de": de, "p": p}


def entropy_change_reversible(q_path, T_path, cumulative: bool = False):
    """Entropy change along a reversible path as the integral of dq/T.

    Book: §1.8, Eq. (1.16) ``s2 − s1 = ∫ dq_rev / T`` (and its differential form (1.17) ``T ds = dq``).

    Parameters
    ----------
    q_path : array_like
        Cumulative heat added along a *reversible* path [J/kg] (e.g. ``process_heat_work(...)["q"]``).
    T_path : array_like
        Temperature at the same samples [K].
    cumulative : bool, optional
        If True return the running integral (same length as the inputs, starting at 0); else the total.

    Returns
    -------
    ds : float or ndarray
        [J/(kg K)].

    Notes
    -----
    Method: trapezoid rule on (1/T) dq. Assumptions: the heat is reversible. For an irreversible process the actual
    ∫δq/T is only a lower bound on Δs (Clausius–Duhem, taught with the actual δq; see analysis §9 typo 4).

    Validation: V1 on four two-leg paths it equals the closed-form Δs of :func:`path_heat_work_totals` (rel 1e-5); V4 a
    closed reversible four-leg cycle returns Δs = −8.4e-6 J/(kg K) (trapezoid, n = 2001 per leg; bound 1e-5).
    Label: analytic, conserved.
    """
    q = np.asarray(q_path, dtype=float)
    T = np.asarray(T_path, dtype=float)
    require_positive("T_path", T)
    ds = 0.5 * (1.0 / T[1:] + 1.0 / T[:-1]) * np.diff(q)  # Eq. (1.17): ds = dq/T (trapezoid)
    s = np.concatenate(([0.0], np.cumsum(ds)))  # Eq. (1.16)
    return s if cumulative else float(s[-1])


def perfect_gas_entropy_change(T1, v1, T2, v2, cv=CV_AIR, R=R_AIR):
    """Entropy change of a perfect gas between two states given by (T, v): integrated Gibbs relation.

    Book: §1.8, Eq. (1.18) ``T ds = de + p dv`` with de = C_v dT and p = R T / v: ``Δs = C_v ln(T2/T1) + R ln(v2/v1)``.

    Parameters
    ----------
    T1, T2 : float or array_like
        Temperatures [K].
    v1, v2 : float or array_like
        Specific volumes [m^3/kg].
    cv : float, optional
        Constant C_v [J/(kg K)].
    R : float, optional
        Gas constant [J/(kg K)].

    Returns
    -------
    ds : float or ndarray
        s2 − s1 [J/(kg K)].

    Notes
    -----
    Assumptions: perfect gas, constant C_v. Valid for **any** process between the states (reversible or not) because
    (1.18) relates state functions only: free expansion (T fixed, v doubles, q = 0) gives Δs = R ln 2 > 0.

    Validation: V2 sympy: the integrated (T, v) form satisfies T ds = de + p dv coefficient by coefficient (D10);
    V1 equals the (T, p) form for heating at constant p, C_p ln 2 = 696.4 J/(kg K) (rel 1e-12); an isentrope from
    (1.26) gives 0 (1e-12); equals the Δs of :func:`path_heat_work_totals` (rel 1e-12). Label: symbolic, analytic.
    """
    ds = cv * np.log(np.asarray(T2, dtype=float) / np.asarray(T1, dtype=float)) + R * np.log(
        np.asarray(v2, dtype=float) / np.asarray(v1, dtype=float))  # Eq. (1.18) integrated: T ds = de + p dv
    return as_scalar_if_0d(ds)


def perfect_gas_entropy_change_p(T1, p1, T2, p2, cp=CP_AIR, R=R_AIR):
    """Entropy change of a perfect gas between two states given by (T, p): ``Δs = C_p ln(T2/T1) − R ln(p2/p1)``.

    Book: §1.8, Eq. (1.18) second form ``T ds = dh − v dp`` with dh = C_p dT and v = R T / p.

    Parameters
    ----------
    T1, T2 : float or array_like
        Temperatures [K].
    p1, p2 : float or array_like
        Absolute pressures [Pa].
    cp : float, optional
        Constant C_p [J/(kg K)].
    R : float, optional
        Gas constant [J/(kg K)].

    Returns
    -------
    ds : float or ndarray
        [J/(kg K)].

    Validation: V2 sympy: the integrated (T, p) form satisfies T ds = dh − v dp; V1 agrees with
    :func:`perfect_gas_entropy_change` for the same states (rel 1e-12) and equals C_p ln 2 at constant p (1e-14).
    Label: symbolic, analytic.
    """
    ds = cp * np.log(np.asarray(T2, dtype=float) / np.asarray(T1, dtype=float)) - R * np.log(
        np.asarray(p2, dtype=float) / np.asarray(p1, dtype=float))  # Eq. (1.18): T ds = dh − v dp
    return as_scalar_if_0d(ds)


def stirred_isochoric_process(v, T1, T2, n: int = 201, cv=CV_AIR, R=R_AIR):
    """Irreversible stirring of a rigid, insulated gas: temperature rises with no heat added.

    Book: §1.8 (the stirring example after the C_v definition and after Eq. (1.18); Clausius–Duhem inequality ii).

    Parameters
    ----------
    v : float
        Specific volume, fixed [m^3/kg].
    T1, T2 : float
        Initial and final temperature [K], T2 >= T1 (stirring can only heat).
    n : int, optional
        Number of samples.
    cv, R : float, optional
        Constant C_v and gas constant [J/(kg K)].

    Returns
    -------
    dict
        ``"T"`` [K], ``"v"`` [m^3/kg], ``"q"`` (zero, J/kg), ``"w"`` stirring work done on the gas (= Δe) [J/kg],
        ``"de"`` [J/kg], ``"s"`` cumulative entropy change from (1.18) [J/(kg K)], ``"int_dq_over_T"`` (zero).

    Notes
    -----
    Shows (a) q ≠ C_v ΔT when work is irreversible, so C_v cannot be *defined* as heat per degree, and (b)
    s2 − s1 = C_v ln(T2/T1) > ∫δq/T = 0 (Clausius–Duhem with the actual δq). Assumptions: perfect gas, constant C_v,
    uniform state at each instant.

    Validation: V1 q = 0 at every sample; final s equals ``irreversible_process("stirring")`` ds = C_v ln(T2/T1)
    (rel 1e-12); ``stirred_isochoric_path`` is the same object. Label: analytic.
    """
    if T2 < T1:
        raise ValueError("stirring work only heats the gas: need T2 >= T1")
    T = np.linspace(T1, T2, n)
    de = cv * (T - T1)
    return {
        "T": T,
        "v": np.full(n, float(v)),
        "q": np.zeros(n),
        "w": de.copy(),  # all energy enters as irreversible stirring work
        "de": de,
        "s": cv * np.log(T / T1),  # Eq. (1.18) at dv = 0: T ds = de
        "int_dq_over_T": np.zeros(n),
    }


#: Alias with the name used in ``analysis/ch01_curation.md`` §8.
stirred_isochoric_path = stirred_isochoric_process


def free_expansion(v1, v2, T, R=R_AIR):
    """Joule free expansion of a perfect gas into vacuum: q = w = Δe = 0 but entropy rises.

    Book: §1.8 (second law ii, Clausius–Duhem: with no heat exchanged, s2 − s1 ≥ 0), evaluated with Eq. (1.18) for a
    perfect gas.

    Parameters
    ----------
    v1, v2 : float
        Specific volumes before and after [m^3/kg], v2 >= v1.
    T : float
        Temperature [K] (unchanged for a perfect gas, e = e(T)).
    R : float, optional
        Gas constant [J/(kg K)].

    Returns
    -------
    dict
        ``"q"``, ``"w"``, ``"de"`` (all 0.0 J/kg), ``"ds"`` = R ln(v2/v1) [J/(kg K)], ``"int_dq_over_T"`` = 0.0.

    Validation: V1 doubling v gives R ln 2 and equals ``irreversible_process("free_expansion")`` exactly.
    Label: analytic.
    """
    if v2 < v1:
        raise ValueError("a free expansion needs v2 >= v1")
    return {"q": 0.0, "w": 0.0, "de": 0.0, "ds": float(R * np.log(v2 / v1)), "int_dq_over_T": 0.0}


# --------------------------------------------------------------------------------------------------------------------
# Non-ideal equations of state used for contrast
# --------------------------------------------------------------------------------------------------------------------
def tait_pressure(rho, rho0=1000.0, K0=2.2e9, n=TAIT_N_WATER, p0=P_ATM):
    """Isentropic (modified Tait / Murnaghan) equation of state of a stiff liquid.

    Book: used with §1.8, Eq. (1.19) to show the incompressible limit c → ∞ (the EOS itself is not in the book).
    Form: ``p = p0 + (K0/n) [ (rho/rho0)^n − 1 ]`` (Cole 1948), so ``(∂p/∂rho)_s = K0/rho0`` at rho = rho0.

    Parameters
    ----------
    rho : float or array_like
        Density [kg/m^3].
    rho0 : float, optional
        Reference density at p0 [kg/m^3] (default 1000, our round value for water).
    K0 : float, optional
        Isentropic bulk modulus at the reference state [Pa] (default 2.2e9, our round value of the order of water's,
        giving c ≈ 1483 m/s at rho0).
    n : float, optional
        Tait exponent [-]; default 7.15 for water (Cole 1948).
    p0 : float, optional
        Reference pressure [Pa].

    Returns
    -------
    p : float or ndarray
        [Pa].

    Notes
    -----
    Assumptions: along one isentrope (entropy fixed), liquid state. Not in the book — our contrast model.

    Validation: V1 sound speed at rho0 from :func:`sound_speed_from_eos` equals sqrt(K0/rho0) (rel 1e-6; 1483 m/s);
    V7 c ∝ sqrt(K0) over K0 = 2.2e7…2.2e11 Pa; p(rho0) = p0; K0 = 0 raises. Label: analytic.
    """
    require_positive("K0", K0)
    return as_scalar_if_0d(p0 + (K0 / n) * ((np.asarray(rho, dtype=float) / rho0) ** n - 1.0))


def van_der_waals_constants_per_mass(a_molar, b_molar, M_w):
    """Convert van der Waals constants from per-mole to per-unit-mass form.

    Parameters
    ----------
    a_molar : float
        [Pa m^6/mol^2].
    b_molar : float
        [m^3/mol].
    M_w : float
        Molecular weight [kg/kmol] (converted to kg/mol internally: kmol vs mol trap).

    Returns
    -------
    (a, b) : tuple of float
        a [Pa m^6/kg^2], b [m^3/kg].

    Validation: V1 CRC CO2 molar constants convert to ``VDW_CO2`` (rel 1e-12; the kmol → mol factor is exercised).
    No pint dimension check. Label: analytic.
    """
    M = M_w / 1.0e3  # kg/kmol -> kg/mol
    return a_molar / M ** 2, b_molar / M


def van_der_waals_pressure(T, v, a, b, R):
    """Van der Waals thermal equation of state per unit mass: ``p = R T/(v − b) − a/v^2``.

    Book: contrast model for §1.8 Eq. (1.12) and the §1.9 claim that only a perfect gas has e = e(T) (not in the book).

    Parameters
    ----------
    T : float or array_like
        Temperature [K].
    v : float or array_like
        Specific volume [m^3/kg], > b.
    a : float
        Attraction constant [Pa m^6/kg^2].
    b : float
        Co-volume [m^3/kg].
    R : float
        Specific gas constant [J/(kg K)].

    Returns
    -------
    p : float or ndarray
        [Pa].

    Notes
    -----
    Assumptions: classical van der Waals fluid. a = b = 0 recovers the perfect gas (1.22).

    Validation: V1 numerical check of the Gibbs-relation identity (∂e/∂v)_T = T (∂p/∂T)_v − p with
    :func:`van_der_waals_internal_energy` for CO2 at 300 K (central differences, rel 1e-5); V7 a = b = 0 equals
    :func:`perfect_gas_pressure` (rel 1e-14). No sympy check. Label: analytic.
    """
    T = np.asarray(T, dtype=float)
    v = np.asarray(v, dtype=float)
    return as_scalar_if_0d(R * T / (v - b) - a / v ** 2)


def van_der_waals_internal_energy(T, v, a, cv, e_ref=0.0):
    """Internal energy of a van der Waals fluid with constant C_v: ``e = C_v T − a/v + e_ref``.

    Book: contrast for the §1.9 claim (Exercise 1.10) that e = e(T) only for a perfect gas. Derived from Gibbs (1.18):
    (∂e/∂v)_T = T (∂p/∂T)_v − p = a/v^2, integrated at fixed T.

    Parameters
    ----------
    T : float or array_like
        Temperature [K].
    v : float or array_like
        Specific volume [m^3/kg].
    a : float
        Attraction constant [Pa m^6/kg^2].
    cv : float
        Constant C_v [J/(kg K)].
    e_ref : float, optional
        Additive constant [J/kg].

    Returns
    -------
    e : float or ndarray
        [J/kg].

    Validation: V1 numerical: (∂e/∂v)_T equals T (∂p/∂T)_v − p from :func:`van_der_waals_pressure` (rel 1e-5) and
    a/v^2 (rel 1e-6); :func:`specific_heat_cv` of it returns C_v (rel 1e-6). No sympy check. Label: analytic.
    """
    return as_scalar_if_0d(cv * np.asarray(T, dtype=float) - a / np.asarray(v, dtype=float) + e_ref)


def helmholtz_free_energy(e, T, s):
    """Helmholtz free energy per unit mass f = e − Ts [J/kg], Eq. (4.94) (Ch. 4 §4.10, surface tension revisited).

    Book: §4.10, Eqs. (4.94)–(4.95): df = de − T ds − s dT; for a reversible isothermal change (1.18) gives
    df = −p dv, the work done on the system. Parameters: e [J/kg], T [K], s [J/(kg K)].
    Validation: V2 sympy: a perfect gas at constant T has df = −p dv. Label: analytic, symbolic.
    """
    return as_scalar_if_0d(np.asarray(e, dtype=float) - np.asarray(T, dtype=float) * np.asarray(s, dtype=float))


__all__ = [
    "helmholtz_free_energy", "G_BOOK",
    "K_B", "N_A", "N_A_KMOL", "R_U", "M_W_AIR", "R_AIR", "GAMMA_AIR", "CP_AIR", "CV_AIR", "G0", "P_ATM", "P_REF",
    "MOLAR_MASS", "GAMMA_BY_ATOMICITY", "GAMMA_IDEAL", "VDW_CO2_MOLAR", "VDW_CO2", "TAIT_N_WATER",
    "molecular_gas_pressure", "gas_constant", "molecule_mass", "molecular_mass", "perfect_gas_pressure", "perfect_gas_density", "perfect_gas_state",
    "specific_volume", "enthalpy", "perfect_gas_internal_energy", "perfect_gas_enthalpy",
    "partial_derivative", "specific_heat_cp", "specific_heat_cv", "sound_speed_from_eos",
    "thermal_expansion_coefficient", "cv_from_cp", "gamma_from_cp", "cp_from_gamma", "cv_from_gamma",
    "isentropic_pressure", "isentropic_ratios", "perfect_gas_sound_speed", "perfect_gas_expansion_coefficient",
    "process_path", "path_heat_work_totals", "irreversible_process", "join_paths", "process_heat_work", "entropy_change_reversible", "perfect_gas_entropy_change",
    "perfect_gas_entropy_change_p", "stirred_isochoric_process", "stirred_isochoric_path", "free_expansion",
    "tait_pressure", "van_der_waals_constants_per_mass", "van_der_waals_pressure", "van_der_waals_internal_energy",
]
