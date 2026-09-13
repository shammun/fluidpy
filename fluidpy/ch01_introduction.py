"""Chapter 1 — Introduction: continuum, molecular transport, surface tension, statics, thermodynamics, stratification,
dimensional analysis.

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e (2012), Ch. 1, §§1.1–1.11, Eqs. (1.1)–(1.40), Examples 1.1–1.5.

This module holds the chapter-specific physics and **re-exports** every reusable primitive the chapter uses from
``fluidpy.core`` (``thermo``, ``statics``, ``stratification``, ``diffusion``, ``dimensional``, ``units``), so that
notebooks and explainer parity rows can call everything as ``ch01.<name>`` (``tools/shot.py`` sees the chapter
module and the core package only).

Conventions used throughout (recorded in ``knowledge/notation.md`` by the knowledge-keeper)
-----------------------------------------------------------------------------------------
* SI units, temperatures in kelvin inside every function (``celsius_to_kelvin`` at the boundary).
* z positive upward (§1.7, §1.10); p0 = pressure at z = 0; ``p_ref`` = reference pressure of θ and rho_θ.
* Lapse rate Γ ≡ dT/dz (book sign; negative when temperature falls with height). ``adiabatic_lapse_rate`` is negative.
* First law: q = heat added to, w = work done on the system (per unit mass).
* Example 1.1: the book's α is measured from the horizontal, so the vertical tension is σ sin α; α is the complement
  of the usual contact angle θ_c (sin α = cos θ_c). ``alpha_from_contact_angle`` converts.
* Book-quoted numbers are not in this file (they live in the git-ignored ``tests/book_values_ch01.json``); constants
  below cite public sources.
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np

from .core._util import as_scalar_if_0d, require_nonnegative, require_positive
from .core.units import KELVIN_OFFSET, celsius_to_kelvin, kelvin_to_celsius  # noqa: F401  (re-export)
from .core.thermo import (  # noqa: F401  (re-exports for ch01.<name>)
    CP_AIR, CV_AIR, G0, GAMMA_AIR, GAMMA_BY_ATOMICITY, GAMMA_IDEAL, K_B, M_W_AIR, MOLAR_MASS, N_A, N_A_KMOL, P_ATM,
    P_REF, R_AIR, R_U, TAIT_N_WATER, VDW_CO2, VDW_CO2_MOLAR,
    cp_from_gamma, cv_from_cp, cv_from_gamma, enthalpy, entropy_change_reversible, free_expansion, gamma_from_cp,
    gas_constant, irreversible_process, isentropic_pressure, isentropic_ratios, join_paths, molecular_gas_pressure,
    molecule_mass, partial_derivative, path_heat_work_totals, perfect_gas_density, perfect_gas_enthalpy,
    perfect_gas_entropy_change, perfect_gas_entropy_change_p, perfect_gas_expansion_coefficient,
    perfect_gas_internal_energy, perfect_gas_pressure, perfect_gas_sound_speed, perfect_gas_state, process_heat_work,
    process_path, sound_speed_from_eos, specific_heat_cp, specific_heat_cv, specific_volume, stirred_isochoric_path,
    stirred_isochoric_process, tait_pressure, thermal_expansion_coefficient, van_der_waals_constants_per_mass,
    van_der_waals_internal_energy, van_der_waals_pressure,
)
from .core.statics import (  # noqa: F401
    USSA_BASES, USSA_LAPSE, USSA_LAPSE0, USSA_P0, USSA_T0, USSA_Z1, USSA_Z_TOP,
    absolute_pressure, atmosphere_from_temperature, buoyancy_force, gauge_pressure, hydrostatic_pressure_uniform,
    integrate_hydrostatic, isothermal_density, isothermal_pressure, layered_pressure, linear_lapse_pressure,
    net_pressure_force_on_box, scale_height, standard_atmosphere,
)
from .core.stratification import (  # noqa: F401
    adiabatic_lapse_rate, brunt_vaisala_sq, brunt_vaisala_sq_from_lapse, brunt_vaisala_sq_from_theta,
    classify_stability, isentropic_density_gradient, lapse_rate, ocean_potential_density_gradient,
    parcel_acceleration_atmosphere, parcel_displacement, parcel_ode, parcel_ode_atmosphere, parcel_ode_from_gradients,
    parcel_temperature, potential_density, potential_temperature, potential_temperature_gradient, stability_timescale,
    temperature_from_potential,
)
from .core.diffusion import (  # noqa: F401
    couette_startup_profile, derivative_2nd_order, ftcs_diffusion_1d, ftcs_stable_time_step, gaussian_spreading,
    stable_time_step,
)
from .core.dimensional import (  # noqa: F401
    BASIS, BLAST, LATEX_SYMBOLS, PENDULUM, PIPE, PRESET_INFO, PYTHAGORAS, RAYLEIGH, SCALE_HEIGHT, SPHERE_DRAG,
    UNIT_SYSTEMS, dimension_vector, dimensional_matrix, exponent_matrix, group_dimension, group_expression, group_latex,
    group_value, groups_independent, minor_determinant, pi_groups, rank_by_minors, rescale_units, solve_exponents,
)

# ====================================================================================================================
# §1.3 Solids, liquids and gases
# ====================================================================================================================
def shear_deformation_history(t, tau, kind: str = "fluid", G=None, mu=None, t_off=None, tau_y=None):
    """Shear strain angle γ(t) of a material element under a constant shear stress applied at t = 0 (removed at t_off).

    Book: §1.3 (definition of a fluid: it deforms continuously under any shear stress, however small; a solid deflects
    to a fixed shape and springs back; plastic and viscoelastic materials in between; Fig. 1.1).
    Models (ours, not in the book): elastic solid γ = τ/G; Newtonian fluid dγ/dt = τ/μ (from Eq. (1.3) with
    du/dy = dγ/dt); Bingham plastic dγ/dt = (τ − τ_y)/μ above the yield stress τ_y, rigid below; Maxwell viscoelastic
    γ = τ/G + τ t/μ, whose elastic part springs back on unloading.

    Parameters
    ----------
    t : float or array_like
        Time [s] (the load acts for 0 <= t < t_off).
    tau : float
        Applied shear stress [Pa], >= 0.
    kind : {"solid", "fluid", "bingham", "maxwell"}, optional
        Material model.
    G : float, optional
        Elastic shear modulus [Pa] (solid, maxwell).
    mu : float, optional
        Viscosity [Pa s] (fluid, bingham, maxwell); for Bingham the plastic viscosity.
    t_off : float, optional
        Time the load is removed [s]; default never.
    tau_y : float, optional
        Yield stress [Pa] (bingham; None means 0).

    Returns
    -------
    gamma : float or ndarray
        Shear strain angle [rad] (small-strain measure).

    Notes
    -----
    After unloading: solid returns to 0; fluid and Bingham keep their final strain (zero memory); Maxwell recovers
    only its elastic part τ/G (partial memory, "viscoelastic"). Assumptions: small strains, constant stress while loaded.

    Validation (planned): V1 solid constant τ/G then 0; fluid τ t/μ then τ t_off/μ; V7 μ → ∞ gives no flow, τ <= τ_y
    Bingham gives no flow. Label: pending.
    """
    t = np.asarray(t, dtype=float)
    require_nonnegative("tau", tau)
    t_off = np.inf if t_off is None else float(t_off)
    loaded = (t >= 0) & (t < t_off)
    t_flow = np.clip(t, 0.0, t_off)  # time the stress has acted so far
    kind = kind.lower()
    if kind == "solid":
        if G is None:
            raise ValueError("solid needs G")
        gamma = np.where(loaded, tau / G, 0.0)  # Hookean: γ = τ/G, springs back when unloaded
    elif kind == "fluid":
        if mu is None:
            raise ValueError("fluid needs mu")
        require_positive("mu", mu)
        gamma = tau * t_flow / mu  # Eq. (1.3): τ = μ dγ/dt  →  γ = τ t/μ, kept after unloading
    elif kind == "bingham":
        if mu is None:
            raise ValueError("bingham needs mu")
        gamma = max(tau - (tau_y or 0.0), 0.0) * t_flow / mu  # flows only above the yield stress
    elif kind == "maxwell":
        if G is None or mu is None:
            raise ValueError("maxwell needs G and mu")
        gamma = np.where(loaded, tau / G, 0.0) + tau * t_flow / mu  # elastic part recovers, viscous part stays
    else:
        raise ValueError(f"unknown kind {kind!r}")
    return as_scalar_if_0d(np.where(t < 0, 0.0, gamma))


def traction_components(traction, normal):
    """Split a traction (stress vector on a surface) into its normal stress and its shear-stress vector.

    Book: §1.3 (normal vs shear stresses; compression and tension; a fluid at rest carries no shear, §1.7).

    Parameters
    ----------
    traction : array_like, shape (3,) or (..., 3)
        Force per unit area exerted on the face by the material outside it [Pa] (tuples and lists accepted).
    normal : array_like, shape (3,) or (..., 3)
        Outward normal of the face; normalised internally [-].

    Returns
    -------
    (sigma_n, tau_vec, tau_mag) : tuple
        ``sigma_n`` = t·n̂ [Pa], signed: positive = tension (pulling outward), negative = compression;
        ``tau_vec`` = t − sigma_n n̂ [Pa], the tangential (shear) part, orthogonal to n̂;
        ``tau_mag`` = |tau_vec| [Pa].

    Validation (planned): V1 tau_vec · n̂ = 0 and sigma_n n̂ + tau_vec = t; pure pressure gives tau_mag = 0. Label: pending.
    """
    t = np.asarray(traction, dtype=float)
    n = np.asarray(normal, dtype=float)
    n = n / np.linalg.norm(n, axis=-1, keepdims=True)
    sigma_n = np.sum(t * n, axis=-1)  # dot product: normal component
    tau_vec = t - np.asarray(sigma_n)[..., None] * n  # remaining tangential part
    tau_mag = np.linalg.norm(tau_vec, axis=-1)
    return as_scalar_if_0d(sigma_n), tau_vec, as_scalar_if_0d(tau_mag)


molecular_mass = molecule_mass  # alias (the chapter module's earlier name)


def number_density(rho, M_w):
    """Number of molecules per unit volume from density and molecular weight: ``n/V = rho A_o / M_w``.

    Book: §1.3–§1.4 (molecular picture), using rho = m n/V from §1.9.

    Parameters
    ----------
    rho : float or array_like
        Density [kg/m^3].
    M_w : float or array_like
        Molecular weight [kg/kmol].

    Returns
    -------
    n : float or ndarray
        [1/m^3] (air at sea level ≈ 2.5e25).

    Validation (planned): V1 equals p/(k_B T) for a perfect gas. Label: pending.
    """
    return as_scalar_if_0d(np.asarray(rho, dtype=float) * N_A_KMOL / np.asarray(M_w, dtype=float))  # rho = m n/V


def mean_molecular_spacing(rho, M_w):
    """Mean distance between neighbouring molecules, ``(n/V)^(-1/3)``.

    Book: §1.3 (gases vs liquids: gas molecules are much farther apart; C04 in the curation).

    Parameters
    ----------
    rho : float or array_like
        Density [kg/m^3].
    M_w : float or array_like
        Molecular weight [kg/kmol].

    Returns
    -------
    spacing : float or ndarray
        [m] (water ≈ 0.31 nm, sea-level air ≈ 3.4 nm).

    Validation (planned): V1 closed form; ratio air/water ≈ 11. Label: pending.
    """
    return as_scalar_if_0d(number_density(rho, M_w) ** (-1.0 / 3.0))


# ====================================================================================================================
# §1.4 Continuum hypothesis
# ====================================================================================================================
def maxwellian_velocities(n: int, T: float, m: float, rng=None, seed: int = 0):
    """Random molecular velocities from the Maxwell–Boltzmann distribution (Gaussian components).

    Book: §1.4 (pressure as the average of molecular collisions; kinetic picture of temperature). Kinetic theory
    (not printed in the book): each Cartesian component is normal with zero mean and variance k_B T/m.

    Parameters
    ----------
    n : int
        Number of molecules.
    T : float
        Temperature [K].
    m : float
        Molecular mass [kg] (see :func:`molecular_mass`).
    rng : numpy.random.Generator, optional
        Random generator; default ``np.random.default_rng(seed)``.
    seed : int, optional
        Seed used when ``rng`` is None (default 0, reproducible).

    Returns
    -------
    u : ndarray, shape (n, 3)
        Velocities [m/s].

    Validation (planned): V1 component variance k_B T/m within sampling error (3σ). Label: pending.
    """
    rng = np.random.default_rng(seed) if rng is None else rng
    return rng.normal(0.0, np.sqrt(K_B * T / m), size=(int(n), 3))  # σ = sqrt(k_B T / m)


def molecular_pressure(number_density, m, velocities):
    """Kinetic-theory pressure from molecular velocities: ``p = (1/3) (N/V) m <|u|^2>``.

    Book: §1.4 (pressure = statistical average of the collision force per unit area); derivation D34 added by the
    curation (momentum 2 m u_x per wall bounce, <u_x^2> = <|u|^2>/3).

    Parameters
    ----------
    number_density : float
        Number of molecules per unit volume N/V [1/m^3].
    m : float
        Molecular mass [kg].
    velocities : array_like, shape (N, 3)
        Sampled molecular velocities [m/s].

    Returns
    -------
    p : float
        [Pa].

    Notes
    -----
    Assumptions: ideal (non-interacting) molecules, isotropic velocity distribution, elastic wall collisions.

    Validation (planned): V4 with Maxwellian samples (N = 1e6, seed 0) equals n k_B T (Eq. 1.21) within 0.5 %. Label: pending.
    """
    u = np.asarray(velocities, dtype=float)
    mean_sq = float(np.mean(np.sum(u ** 2, axis=-1)))
    return number_density * m * mean_sq / 3.0  # D34: p = (1/3) n m <|u|^2>


def wall_impact_pressure(velocities, m, number_density, dt=1e-12, area=1.0):
    """Pressure as momentum delivered to a wall per unit area and time by the molecules that hit it (D34, step by step).

    Book: §1.4 (pressure as the average collision force per unit area); derivation D34 added by the curation.
    Expected-count form: each sampled velocity stands for a class of n/N molecules per unit volume; those with u_x > 0
    within u_x dt of the wall hit it in dt, i.e. ``(n/N) · area · u_x · dt`` molecules, each delivering ``2 m u_x``.
    Pressure = total momentum / (area · dt) = ``(2 n m / N) Σ_{u_x>0} u_x^2``.

    Parameters
    ----------
    velocities : array_like, shape (N, 3)
        Sampled molecular velocities [m/s] (x = wall-normal direction, wall at +x).
    m : float
        Molecular mass [kg].
    number_density : float
        Molecules per unit volume [1/m^3].
    dt : float, optional
        Counting interval [s] (cancels; kept to show the construction).
    area : float, optional
        Wall area [m^2] (cancels).

    Returns
    -------
    p : float
        [Pa]; its expectation is n m <u_x^2> = n k_B T for a Maxwellian gas.

    Validation (planned): V4 equals :func:`molecular_pressure` and n k_B T within sampling error. Label: pending.
    """
    u = np.asarray(velocities, dtype=float)
    ux = u[:, 0]
    hitting = ux > 0.0  # only molecules moving toward the wall reach it
    count_per_class = number_density / u.shape[0] * area * ux[hitting] * dt  # molecules of each class hitting in dt
    momentum = np.sum(count_per_class * 2.0 * m * ux[hitting])  # each bounce reverses u_x: 2 m u_x
    return float(momentum / (area * dt))  # force per area


def density_expected(L, number_density, m, gradient=0.0, L_flow=1.0):
    """Expected density measured in a cube of side L when the true density varies linearly across the cube.

    Book: §1.4 (continuum: properties are averages over a small region around the point). Model (ours): the cube spans
    [0, L] from the point along x and the density is rho (1 + gradient · x / L_flow), so the box average is
    ``n m (1 + gradient · L / (2 L_flow))``.

    Parameters
    ----------
    L : float or array_like
        Box side [m].
    number_density : float
        Number density at the point [1/m^3].
    m : float
        Molecular mass [kg].
    gradient : float, optional
        Fractional density change per ``L_flow`` [-].
    L_flow : float, optional
        Length over which the fractional change ``gradient`` occurs (the flow's own length scale) [m].

    Returns
    -------
    rho_box : float or ndarray
        [kg/m^3].

    Validation (planned): V1 gradient = 0 gives n m for every L. Label: pending.
    """
    return as_scalar_if_0d(number_density * m * (1.0 + gradient * np.asarray(L, dtype=float) / (2.0 * L_flow)))


def density_noise_expected(L, number_density):
    """Expected relative fluctuation of the molecule count in a cube of side L: ``(n L^3)^(-1/2)``.

    Book: §1.4 (continuum hypothesis); derivation D35 added by the curation (Poisson counting: std/mean = N^(-1/2),
    so the relative noise falls as L^(-3/2)).

    Parameters
    ----------
    L : float or array_like
        Box side [m].
    number_density : float
        Molecules per unit volume [1/m^3].

    Returns
    -------
    rel_noise : float or ndarray
        std(rho)/mean(rho) [-].

    Notes
    -----
    Deterministic companion of :func:`sample_density` for explainer parity (seeded draws cannot be matched in JS).

    Validation (planned): V1 slope −3/2 on log–log; V4 matches sample statistics within sampling error. Label: pending.
    """
    return as_scalar_if_0d((number_density * np.asarray(L, dtype=float) ** 3) ** -0.5)  # D35


def sample_density(L_box, number_density, m, rng=None, n_samples: int = 200, gradient: float = 0.0,
                   L_flow: float = 1.0, seed: int = 0):
    """Measured density δm/δV in repeated random cubes of side L: sample mean and standard deviation.

    Book: §1.4 (continuum hypothesis: rho(x) = δm/δV is meaningful only when δV holds very many molecules yet is small
    compared with the flow's length scale).

    Parameters
    ----------
    L_box : float or array_like
        Cube side(s) [m].
    number_density : float
        Molecules per unit volume at the point [1/m^3].
    m : float
        Molecular mass [kg].
    rng : numpy.random.Generator, optional
        Random generator; default ``np.random.default_rng(seed)``.
    n_samples : int, optional
        Number of independent boxes per size (default 200).
    gradient : float, optional
        Fractional density change per ``L_flow`` across the box (see :func:`density_expected`).
    L_flow : float, optional
        Macroscopic length scale [m] (default 1 m).
    seed : int, optional
        Seed used when ``rng`` is None (default 0).

    Returns
    -------
    (mean, std) : tuple of float or ndarray
        Sample mean and sample standard deviation of δm/δV [kg/m^3], same shape as ``L_box``.

    Notes
    -----
    Method: molecule counts are Poisson with mean λ = n L^3 (1 + gradient L/(2 L_flow)) (ideal-gas positions are
    independent); for λ > 1e7 the normal approximation λ + sqrt(λ) Z is used (numpy's Poisson sampler overflows).
    Assumptions: uncorrelated molecular positions (ideal gas); a liquid's short-range order reduces the noise.

    Validation (planned): V1 relative std ∝ L^(-3/2) (slope −1.5 ± 0.05); V7 mean → rho for large L with no gradient.
    Label: pending.
    """
    rng = np.random.default_rng(seed) if rng is None else rng
    L = np.atleast_1d(np.asarray(L_box, dtype=float))
    means = np.empty_like(L)
    stds = np.empty_like(L)
    for i, Li in enumerate(L):
        lam = number_density * Li ** 3 * max(1.0 + gradient * Li / (2.0 * L_flow), 0.0)
        if lam > 1e7:
            counts = lam + np.sqrt(lam) * rng.standard_normal(n_samples)
        else:
            counts = rng.poisson(lam, n_samples).astype(float)
        rho = counts * m / Li ** 3  # δm/δV
        means[i] = rho.mean()
        stds[i] = rho.std(ddof=1)
    if np.ndim(L_box) == 0:
        return float(means[0]), float(stds[0])
    return means, stds


def knudsen_number(l, L):
    """Knudsen number, the ratio of the molecular mean free path to the length scale of interest.

    Book: §1.4, ``Kn = l/L``; the continuum approximation is valid when Kn << 1.

    Parameters
    ----------
    l : float or array_like
        Mean free path [m].
    L : float or array_like
        Length scale of the problem (body size, pore diameter, …) [m].

    Returns
    -------
    Kn : float or ndarray
        [-].

    Validation (planned): V1 closed form; pint dimensionless. Label: pending.
    """
    require_positive("L", L)
    return as_scalar_if_0d(np.asarray(l, dtype=float) / np.asarray(L, dtype=float))  # §1.4: Kn = l/L


def mean_free_path_jennings(mu, rho, p):
    """Mean free path of gas molecules from viscosity, density and pressure (Jennings 1988 formula).

    Book: §1.4 (the mean free path l in Kn = l/L; the book quotes an order-of-magnitude value for air). Formula (cited
    addition): ``λ = sqrt(π/8) μ / (0.4987445 sqrt(rho p))`` — S. G. Jennings, J. Aerosol Sci. 19, 159 (1988),
    doi:10.1016/0021-8502(88)90219-4, as reported by Tsalikis, Mavrantzas & Pratsinis, Aerosol Sci. Technol. (2024).

    Parameters
    ----------
    mu : float or array_like
        Dynamic viscosity [Pa s].
    rho : float or array_like
        Density [kg/m^3].
    p : float or array_like
        Pressure [Pa].

    Returns
    -------
    l : float or ndarray
        Mean free path [m] (≈ 67 nm for air at 300 K, 1 atm).

    Validation (planned): V5 67.3 nm at 300 K, 1 atm with Sutherland μ; V6 book's order of magnitude (private JSON).
    Label: pending.
    """
    return as_scalar_if_0d(np.sqrt(np.pi / 8.0) * np.asarray(mu, dtype=float)
                           / (0.4987445 * np.sqrt(np.asarray(rho, dtype=float) * np.asarray(p, dtype=float))))


def mean_free_path_air(T, p=P_ATM):
    """Mean free path of air from temperature and pressure (Jennings formula with Sutherland μ and perfect-gas rho).

    Book: §1.4 (mean free path in Kn = l/L). Combines :func:`sutherland_viscosity` (USSA-1976 constants),
    :func:`perfect_gas_density` (Eq. 1.22) and :func:`mean_free_path_jennings` (Jennings 1988).

    Parameters
    ----------
    T : float or array_like
        Temperature [K].
    p : float or array_like, optional
        Pressure [Pa].

    Returns
    -------
    l : float or ndarray
        [m] (≈ 67 nm at 300 K, 1 atm; grows as 1/p with altitude).

    Validation (planned): V5 67.3 nm at 300 K, 101325 Pa (Jennings via Tsalikis et al. 2024). Label: pending.
    """
    mu = sutherland_viscosity(T)
    rho = perfect_gas_density(p, T)
    return mean_free_path_jennings(mu, rho, p)


# ====================================================================================================================
# §1.5 Molecular transport
# ====================================================================================================================
#: Dry-air composition, U.S. Standard Atmosphere 1976 main constituents: name -> (mole fraction [-], molar mass
#: [kg/kmol]). The remaining trace gases (Ne, He, Kr, ...) make up < 3e-5 of the mixture.
DRY_AIR: dict[str, tuple[float, float]] = {"N2": (0.78084, 28.0134), "O2": (0.209476, 31.9988),
                                           "Ar": (0.00934, 39.948), "CO2": (0.000314, 44.0095)}




def mass_fractions(mole_fractions, molar_masses):
    """Mass fractions of a mixture from its mole fractions: ``Y_i = X_i M_i / Σ_j X_j M_j``.

    Book: §1.5 (mass fraction Y of a constituent; partial density rho Y, defined before Eq. (1.1)).

    Parameters
    ----------
    mole_fractions : mapping or array_like
        Mole (number) fractions X_i [-] (need not sum exactly to 1; they are normalised).
    molar_masses : mapping or array_like
        Molar masses M_i [kg/kmol], same keys/order.

    Returns
    -------
    Y : dict or ndarray
        Mass fractions [-] summing to 1 (dict if a mapping was given).

    Validation (planned): V1 sums to 1; mean molar mass of the USSA composition ≈ 28.96. Label: pending.
    """
    if isinstance(mole_fractions, Mapping):
        keys = list(mole_fractions)
        X = np.array([mole_fractions[k] for k in keys], dtype=float)
        M = np.array([molar_masses[k] for k in keys], dtype=float)
        Y = X * M / np.sum(X * M)
        return {k: float(y) for k, y in zip(keys, Y)}
    X = np.asarray(mole_fractions, dtype=float)
    M = np.asarray(molar_masses, dtype=float)
    return X * M / np.sum(X * M)  # Y_i = X_i M_i / Σ X_j M_j


def fick_mass_flux(rho, kappa_m, grad_Y):
    """Diffusive mass flux of a constituent, down its mass-fraction gradient.

    Book: §1.5, Eq. (1.1) ``J_m = −rho κ_m ∇Y``.

    Parameters
    ----------
    rho : float or array_like
        Mixture density [kg/m^3].
    kappa_m : float
        Mass diffusivity [m^2/s], >= 0 (second law iii).
    grad_Y : float or array_like
        Mass-fraction gradient [1/m] (a vector as the last axis, or a scalar component).

    Returns
    -------
    J_m : float or ndarray
        Mass flux [kg/(m^2 s)], antiparallel to ∇Y.

    Notes
    -----
    Assumptions: binary diffusion, linear in the gradient, isothermal and isobaric (no Soret/pressure diffusion).

    Validation (planned): V1 J · ∇Y < 0; pint kg m^-2 s^-1; negative κ_m raises. Label: pending.
    """
    require_nonnegative("kappa_m", kappa_m)
    # numpy broadcasting: for a field of vectors pass rho with a trailing axis, e.g. rho[..., None]
    return as_scalar_if_0d(-np.asarray(rho, dtype=float) * kappa_m * np.asarray(grad_Y, dtype=float))  # Eq. (1.1)


def fourier_heat_flux(k, grad_T):
    """Conductive heat flux, down the temperature gradient.

    Book: §1.5, Eq. (1.2) ``q = −k ∇T``.

    Parameters
    ----------
    k : float
        Thermal conductivity [W/(m K)], >= 0.
    grad_T : float or array_like
        Temperature gradient [K/m].

    Returns
    -------
    q : float or ndarray
        Heat flux [W/m^2] (= J m^-2 s^-1).

    Validation (planned): V1 antiparallel to ∇T; pint W m^-2; negative k raises. Label: pending.
    """
    require_nonnegative("k", k)
    return as_scalar_if_0d(-k * np.asarray(grad_T, dtype=float))  # Eq. (1.2)


def newton_shear_stress(mu, du_dy):
    """Shear stress in a unidirectional shear flow (Newton's law of friction).

    Book: §1.5, Eq. (1.3) ``τ = μ (du/dy)``.

    Parameters
    ----------
    mu : float or array_like
        Dynamic viscosity [Pa s], >= 0.
    du_dy : float or array_like
        Velocity gradient [1/s].

    Returns
    -------
    tau : float or ndarray
        [Pa]; the momentum flux from fast to slow fluid.

    Notes
    -----
    Assumptions: Newtonian fluid, u = u(y) only (the tensor form comes in Ch. 4).

    Validation (planned): V1 linear profile gives uniform τ; pint Pa; negative μ raises. Label: pending.
    """
    require_nonnegative("mu", mu)
    return as_scalar_if_0d(np.asarray(mu, dtype=float) * np.asarray(du_dy, dtype=float))  # Eq. (1.3)


def shear_stress_profile(u, y, mu):
    """Shear stress τ(y) = μ du/dy along a sampled velocity profile.

    Book: §1.5, Eq. (1.3), evaluated on data (Fig. 1.3).

    Parameters
    ----------
    u : array_like, shape (N,)
        Velocity samples [m/s].
    y : array_like, shape (N,)
        Wall-normal coordinate [m], monotonic, N >= 3.
    mu : float
        Dynamic viscosity [Pa s].

    Returns
    -------
    tau : ndarray, shape (N,)
        [Pa].

    Notes
    -----
    Method: second-order central differences in the interior and second-order one-sided stencils at the ends, written
    out explicitly (:func:`fluidpy.core.diffusion.derivative_2nd_order`).

    Validation (planned): V1 parabolic channel profile u = U0 (1 − (y/b)^2) gives τ = −2 μ U0 y/b^2 exactly (quadratic);
    V3 order 2 ± 0.15 on a sine profile. Label: pending.
    """
    require_nonnegative("mu", mu)
    return mu * derivative_2nd_order(u, y)  # Eq. (1.3) with 2nd-order stencils


def wall_shear_history(u_hist, dy: float, mu: float):
    """Wall shear stress μ ∂u/∂y at both walls for every stored profile of a plate-driven flow.

    Book: §1.5, Eq. (1.3) at y = 0 and y = h (the E2 explainer's τ_w(t) panel).

    Parameters
    ----------
    u_hist : array_like, shape (nt, N)
        Velocity profiles on a uniform grid from y = 0 to y = h (e.g. output of ``ftcs_diffusion_1d``) [m/s].
    dy : float
        Grid spacing [m].
    mu : float
        Dynamic viscosity [Pa s].

    Returns
    -------
    (tau_bottom, tau_top) : tuple of ndarray, shape (nt,)
        μ ∂u/∂y at y = 0 and y = h [Pa] (both → μ U/h in steady Couette flow).

    Notes
    -----
    Stencils: (−3u_0 + 4u_1 − u_2)/(2Δy) and (3u_{N−1} − 4u_{N−2} + u_{N−3})/(2Δy), second order.

    Validation (planned): V1 linear profile gives μ U/h at both walls. Label: pending.
    """
    u = np.atleast_2d(np.asarray(u_hist, dtype=float))
    tau_bottom = mu * (-3.0 * u[:, 0] + 4.0 * u[:, 1] - u[:, 2]) / (2.0 * dy)  # Eq. (1.3) at y = 0
    tau_top = mu * (3.0 * u[:, -1] - 4.0 * u[:, -2] + u[:, -3]) / (2.0 * dy)  # Eq. (1.3) at y = h
    return tau_bottom, tau_top


def viscosity_power_law(T, mu_ref, T_ref, n=0.5):
    """Power-law temperature dependence of viscosity, ``μ = μ_ref (T/T_ref)^n``.

    Book: §1.5 (for gases μ varies approximately as T^(1/2), because the random thermal speed does).

    Parameters
    ----------
    T : float or array_like
        Temperature [K].
    mu_ref : float
        Viscosity at T_ref [Pa s].
    T_ref : float
        Reference temperature [K].
    n : float, optional
        Exponent [-] (0.5 = the book's kinetic-theory estimate; ≈ 0.7 fits air better).

    Returns
    -------
    mu : float or ndarray
        [Pa s].

    Validation (planned): V1 μ(T_ref) = μ_ref; V5 relative error vs Sutherland at 216.65 K reported (≈ 9 % for n = 0.5).
    Label: pending.
    """
    return as_scalar_if_0d(mu_ref * (np.asarray(T, dtype=float) / T_ref) ** n)  # §1.5: μ ∝ T^(1/2) for gases


def sutherland_viscosity(T, beta=1.458e-6, S=110.4):
    """Viscosity of air from Sutherland's law with the U.S. Standard Atmosphere 1976 constants.

    Book: §1.5 (gas viscosity rises with temperature). Benchmark model (cited addition): ``μ = β T^(3/2) / (T + S)``
    with β = 1.458e-6 kg/(m s K^(1/2)) and S = 110.4 K (NASA-TM-X-74335).

    Parameters
    ----------
    T : float or array_like
        Temperature [K].
    beta : float, optional
        [kg m^-1 s^-1 K^-1/2].
    S : float, optional
        Sutherland temperature [K].

    Returns
    -------
    mu : float or ndarray
        [Pa s] (1.7894e-5 at 288.15 K).

    Validation (planned): V5 USSA-1976 table at 0 km (1.7894e-5) and 11–20 km. Label: pending.
    """
    T = np.asarray(T, dtype=float)
    require_positive("T", T)
    return as_scalar_if_0d(beta * T ** 1.5 / (T + S))


def water_viscosity(T):
    """Dynamic viscosity of liquid water at atmospheric pressure (Vogel-type fit).

    Book: §1.5 (the viscosity of a liquid decreases with temperature). Fit (cited addition):
    ``μ = A · 10^(B/(T − C))`` with A = 2.414e-5 Pa s, B = 247.8 K, C = 140 K (widely reproduced engineering fit,
    e.g. T. Al-Shemmeri, *Engineering Fluid Mechanics* (2012); accuracy about ±1 % from 0 to 100 °C).

    Parameters
    ----------
    T : float or array_like
        Temperature [K], about 273–373 K.

    Returns
    -------
    mu : float or ndarray
        [Pa s] (≈ 1.00e-3 at 20 °C).

    Validation (planned): V5 against a published table (e.g. IAPWS 2008 at 20 °C: 1.0016e-3) within 1 %. Label: pending.
    """
    T = np.asarray(T, dtype=float)
    if np.any(T <= 140.0):
        raise ValueError("water_viscosity fit needs T > 140 K (use 273–373 K)")
    return as_scalar_if_0d(2.414e-5 * 10.0 ** (247.8 / (T - 140.0)))


def kinematic_viscosity(mu, rho):
    """Kinematic viscosity, the diffusivity of momentum.

    Book: §1.5, Eq. (1.4) ``ν ≡ μ/rho``.

    Parameters
    ----------
    mu : float or array_like
        Dynamic viscosity [Pa s].
    rho : float or array_like
        Density [kg/m^3], > 0.

    Returns
    -------
    nu : float or ndarray
        [m^2/s].

    Validation (planned): V1 closed form; V5 USSA sea level ν = 1.4607e-5 m^2/s. Label: pending.
    """
    require_positive("rho", rho)
    return as_scalar_if_0d(np.asarray(mu, dtype=float) / np.asarray(rho, dtype=float))  # Eq. (1.4)


def thermal_diffusivity(k, rho, cp):
    """Thermal diffusivity ``κ = k/(rho C_p)``, the diffusivity of heat in the model equation ∂T/∂t = κ ∂²T/∂y².

    Book: §1.5 (Fourier's law (1.2) with the energy content rho C_p T; formal energy equation in Ch. 4).

    Parameters
    ----------
    k : float or array_like
        Thermal conductivity [W/(m K)].
    rho : float or array_like
        Density [kg/m^3].
    cp : float or array_like
        Specific heat [J/(kg K)].

    Returns
    -------
    kappa : float or ndarray
        [m^2/s].

    Validation (planned): V1 closed form; pint m^2/s. Label: pending.
    """
    return as_scalar_if_0d(np.asarray(k, dtype=float) / (np.asarray(rho, dtype=float) * np.asarray(cp, dtype=float)))


def diffusion_time(L, D):
    """Time scale for diffusion across a distance L: ``t ~ L^2 / D``.

    Book: §1.5 (diffusive smoothing); scaling from the model diffusion equation (ours).

    Parameters
    ----------
    L : float or array_like
        Distance [m].
    D : float or array_like
        Diffusivity (ν, κ or κ_m) [m^2/s].

    Returns
    -------
    t : float or ndarray
        [s].

    Validation (planned): V1 closed form; Couette start-up reaches steady state for t ≳ 0.5 h^2/ν. Label: pending.
    """
    require_positive("D", D)
    return as_scalar_if_0d(np.asarray(L, dtype=float) ** 2 / np.asarray(D, dtype=float))


# ====================================================================================================================
# §1.6 Surface tension
# ====================================================================================================================
IAPWS_TC: float = 647.096  #: critical temperature of water [K] (IAPWS R1-76(2014))


def surface_tension_water(T):
    """Surface tension of water against its vapour (≈ air), IAPWS release R1-76(2014).

    Book: §1.6 (surface tension σ as force per unit length or energy per unit area; depends on temperature).
    Correlation (cited addition): ``σ = B τ^μ (1 + b τ)``, τ = 1 − T/T_c, B = 235.8 mN/m, b = −0.625, μ = 1.256,
    T_c = 647.096 K (https://iapws.org/public/documents/CH-L9/Surf-H2O-2014.pdf).

    Parameters
    ----------
    T : float or array_like
        Temperature [K], 248.15 K (supercooled extrapolation) … T_c.

    Returns
    -------
    sigma : float or ndarray
        [N/m] (0.07274 at 20 °C).

    Validation (planned): V5 IAPWS Table 1 at 20 °C (72.74 mN/m) and 25 °C (71.97); V7 σ → 0 at T_c, monotone. Label: pending.
    """
    T = np.asarray(T, dtype=float)
    if np.any(T > IAPWS_TC) or np.any(T < 248.15):
        raise ValueError("surface_tension_water is valid for 248.15 K <= T <= 647.096 K")
    tau = 1.0 - T / IAPWS_TC
    return as_scalar_if_0d(235.8e-3 * tau ** 1.256 * (1.0 - 0.625 * tau))


def laplace_pressure_jump(sigma, R1, R2=None):
    """Pressure jump across a curved interface (Laplace pressure).

    Book: §1.6, Eq. (1.5) ``p_i − p_o = σ (1/R1 + 1/R2)``; a sphere (R1 = R2 = R) gives 2σ/R.

    Parameters
    ----------
    sigma : float or array_like
        Surface tension [N/m].
    R1 : float or array_like
        First principal radius of curvature [m]; ``np.inf`` for a flat direction.
    R2 : float or array_like, optional
        Second principal radius [m]; default R1 (sphere). Radii are **signed**: positive when the centre of curvature is
        on the inner (i) side; a saddle has radii of opposite sign.

    Returns
    -------
    dp : float or ndarray
        p_i − p_o [Pa], higher on the concave side.

    Notes
    -----
    Assumptions: static interface, uniform σ.

    Validation (planned): V1 sphere 2σ/R; flat (∞, ∞) gives 0; saddle R2 = −R1 gives 0; pint Pa. Label: pending.
    """
    R1 = np.asarray(R1, dtype=float)
    R2 = R1 if R2 is None else np.asarray(R2, dtype=float)
    with np.errstate(divide="ignore"):
        curv = np.where(np.isinf(R1), 0.0, 1.0 / R1) + np.where(np.isinf(R2), 0.0, 1.0 / R2)
    return as_scalar_if_0d(np.asarray(sigma, dtype=float) * curv)  # Eq. (1.5)


# ====================================================================================================================
# §1.7 Fluid statics
# ====================================================================================================================
def wedge_pressure_difference(rho, dz, theta, g=G0):
    """Pressure differences on the faces of a small triangular fluid wedge at rest.

    Book: §1.7, force balances on the wedge of Fig. 1.5 leading to Eq. (1.6): x-balance gives p1 = p3; z-balance gives
    ``p2 − p1 − (1/2) rho g dz = 0``; as dz → 0 at fixed θ, p1 = p2 = p3.

    Parameters
    ----------
    rho : float
        Density [kg/m^3].
    dz : float or array_like
        Height of the wedge (length of the vertical face) [m].
    theta : float or array_like
        Angle between the sloping face and the horizontal [rad], 0 < θ < π/2.
    g : float, optional
        Gravitational acceleration [m/s^2].

    Returns
    -------
    dict
        ``"p2_minus_p1"`` = rho g dz/2 [Pa] (→ 0 linearly with dz, independent of θ);
        ``"p3_minus_p1"`` = 0 [Pa].

    Notes
    -----
    Assumptions: fluid at rest (no shear stress), gravity the only body force, unit thickness.

    Validation (planned): V1 p3 − p1 = 0; V3 p2 − p1 ∝ dz (order 1 ± 0.05); independent of θ. Label: pending.
    """
    dz = np.asarray(dz, dtype=float)
    theta = np.asarray(theta, dtype=float)
    return {"p2_minus_p1": as_scalar_if_0d(0.5 * rho * g * dz + 0.0 * theta),  # §1.7: p2 − p1 − (1/2) rho g dz = 0
            "p3_minus_p1": as_scalar_if_0d(0.0 * dz * theta)}  # §1.7: (p1 ds) sin θ − p3 dz = 0 with dz = sin θ ds


def wedge_face_forces(p1, rho, dz, theta, g=G0):
    """Forces (per unit thickness) on the three faces of the Fig. 1.5 wedge and the residuals of both balances.

    Book: §1.7, ``(p1 ds) sin θ − p3 dz = 0`` and ``−(p1 ds) cos θ + p2 dx − (1/2) rho g dx dz = 0`` with
    dz = sin θ ds and dx = cos θ ds.

    Parameters
    ----------
    p1 : float
        Pressure on the sloping face [Pa].
    rho : float
        Density [kg/m^3].
    dz : float or array_like
        Vertical face length [m].
    theta : float or array_like
        Slope angle [rad].
    g : float, optional
        [m/s^2].

    Returns
    -------
    dict
        ``ds``, ``dx`` [m]; ``p2``, ``p3`` [Pa] from the balances; force magnitudes per unit thickness [N/m]:
        ``F1 = p1 ds`` (sloping face), ``F2 = p2 dx`` (bottom), ``F3 = p3 dz`` (vertical face),
        ``W = (1/2) rho g dx dz`` (weight); residuals ``res_x``, ``res_z`` [N/m] (zero by construction);
        ``weight_to_face_force`` = W/F2 [-] (→ 0 as dz → 0: the weight vanishes faster than the face forces).

    Validation (planned): V1 residuals vanish; W/F2 ∝ dz. Label: pending.
    """
    dz = np.asarray(dz, dtype=float)
    theta = np.asarray(theta, dtype=float)
    ds = dz / np.sin(theta)  # dz = sin θ ds
    dx = ds * np.cos(theta)  # cos θ ds = dx
    p3 = p1 * ds * np.sin(theta) / dz  # x-balance
    W = 0.5 * rho * g * dx * dz
    p2 = (p1 * ds * np.cos(theta) + W) / dx  # z-balance
    res_x = p1 * ds * np.sin(theta) - p3 * dz
    res_z = -p1 * ds * np.cos(theta) + p2 * dx - W
    out = {"ds": ds, "dx": dx, "p2": p2, "p3": p3, "F1": p1 * ds, "F2": p2 * dx, "F3": p3 * dz, "W": W,
           "res_x": res_x, "res_z": res_z, "weight_to_face_force": W / (p2 * dx)}
    return {k: as_scalar_if_0d(v) for k, v in out.items()}


def alpha_from_contact_angle(theta_c):
    """Convert a conventional contact angle θ_c (measured through the liquid from the wall) to the book's α.

    Book: Example 1.1 / Fig. 1.7 draws α between the meniscus tangent and the horizontal, so the vertical tension
    component is σ sin α; hence α = π/2 − θ_c (analysis §9, convention 3).

    Parameters
    ----------
    theta_c : float or array_like
        Contact angle [rad] (0 = perfect wetting).

    Returns
    -------
    alpha : float or ndarray
        The book's angle [rad] (π/2 = perfect wetting).

    Validation (planned): V1 θ_c = 0 → π/2; round trip. Label: pending.
    """
    return as_scalar_if_0d(np.pi / 2.0 - np.asarray(theta_c, dtype=float))


contact_angle_to_alpha = alpha_from_contact_angle  # alias


def capillary_rise(sigma, alpha, rho, R, g=G0):
    """Height a liquid climbs in a narrow tube.

    Book: §1.7, Example 1.1: ``σ (2πR) sin α = rho g h (πR^2)`` ⇒ ``h = 2 σ sin α / (rho g R)``.

    Parameters
    ----------
    sigma : float or array_like
        Surface tension [N/m].
    alpha : float or array_like
        The book's angle between the meniscus and the horizontal at the wall [rad] (π/2 = full wetting). **Not** the
        usual contact angle θ_c: α = π/2 − θ_c (see :func:`alpha_from_contact_angle`).
    rho : float or array_like
        Liquid density [kg/m^3] (the gas above is neglected).
    R : float or array_like
        Tube radius [m].
    g : float, optional
        [m/s^2].

    Returns
    -------
    h : float or ndarray
        Rise of the meniscus bottom above the outside free surface [m] (negative for a non-wetting liquid, α < 0).

    Notes
    -----
    Assumptions: narrow tube (spherical meniscus), meniscus volume neglected, static, uniform σ. Consistent with
    (1.5): meniscus radius R/sin α gives p_E = p_atm − 2σ sin α/R = p_atm − rho g h.

    Validation (planned): V1 force-balance residual 0; consistency with :func:`laplace_pressure_jump` and Eq. (1.9);
    α = 0 gives 0; h ∝ 1/R. Label: pending.
    """
    return as_scalar_if_0d(2.0 * np.asarray(sigma, dtype=float) * np.sin(np.asarray(alpha, dtype=float))
                           / (np.asarray(rho, dtype=float) * g * np.asarray(R, dtype=float)))  # Example 1.1


def capillary_rise_deg(sigma, alpha_deg, rho, R, g=G0):
    """:func:`capillary_rise` with the book's angle α given in degrees.

    Book: §1.7, Example 1.1.

    Parameters
    ----------
    sigma : float or array_like
        [N/m].
    alpha_deg : float or array_like
        The book's α [degrees] (90 = full wetting).
    rho, R : float or array_like
        [kg/m^3], [m].
    g : float, optional
        [m/s^2].

    Returns
    -------
    h : float or ndarray
        [m].

    Validation (planned): V1 equals capillary_rise(np.radians(alpha_deg)). Label: pending.
    """
    return capillary_rise(sigma, np.radians(alpha_deg), rho, R, g)


# ====================================================================================================================
# §1.8 Thermodynamics helpers (kinetic time scales, property fits)
# ====================================================================================================================
def mean_molecular_speed(T, m):
    """Mean thermal speed of gas molecules, ``c̄ = sqrt(8 k_B T/(π m))`` (Maxwell–Boltzmann mean).

    Book: §1.8 (relaxation to equilibrium is fast because molecular speeds and collision rates are high). Kinetic
    theory formula (not printed in the book).

    Parameters
    ----------
    T : float or array_like
        Temperature [K].
    m : float or array_like
        Molecular mass [kg].

    Returns
    -------
    c_bar : float or ndarray
        [m/s] (≈ 470 m/s for air at 300 K).

    Validation (planned): V4 equals the mean |u| of :func:`maxwellian_velocities` samples. Label: pending.
    """
    return as_scalar_if_0d(np.sqrt(8.0 * K_B * np.asarray(T, dtype=float) / (np.pi * np.asarray(m, dtype=float))))


def collision_time(l, T, m):
    """Mean time between molecular collisions, ``t_c = l / c̄``.

    Book: §1.8 (relaxation time of a fluid particle; C23 in the curation).

    Parameters
    ----------
    l : float or array_like
        Mean free path [m].
    T : float or array_like
        Temperature [K].
    m : float or array_like
        Molecular mass [kg].

    Returns
    -------
    t_c : float or ndarray
        [s] (≈ 1e-10 s for air at room conditions).

    Validation (planned): V1 closed form. Label: pending.
    """
    return as_scalar_if_0d(np.asarray(l, dtype=float) / mean_molecular_speed(T, m))


def water_density(T):
    """Density of air-free liquid water at 1 atm (Kell 1975 correlation).

    Book: §1.8, used with Eq. (1.20) to show that α of water changes sign near 4 °C. Correlation (cited addition):
    G. S. Kell, J. Chem. Eng. Data 20, 97 (1975):
    ``rho = (999.83952 + 16.945176 t − 7.9870401e-3 t^2 − 46.170461e-6 t^3 + 105.56302e-9 t^4 − 280.54253e-12 t^5)
    / (1 + 16.879850e-3 t)`` with t in °C.

    Parameters
    ----------
    T : float or array_like
        Temperature [K], 273.15–373.15 K (t on IPTS-68; the ITS-90 difference < 0.03 K is ignored).

    Returns
    -------
    rho : float or ndarray
        [kg/m^3] (maximum ≈ 999.97 near 277.1 K).

    Validation (planned): V5 maximum near 3.98 °C; value at 20 °C ≈ 998.2 (NIST tables). Label: pending.
    """
    t = np.asarray(T, dtype=float) - KELVIN_OFFSET
    num = (999.83952 + 16.945176 * t - 7.9870401e-3 * t ** 2 - 46.170461e-6 * t ** 3 + 105.56302e-9 * t ** 4
           - 280.54253e-12 * t ** 5)
    return as_scalar_if_0d(num / (1.0 + 16.879850e-3 * t))


def seawater_density_eos80(T, S):
    """Density of seawater at atmospheric pressure (UNESCO one-atmosphere equation of state, EOS-80).

    Book: §1.10 (seawater density depends on temperature, pressure and salinity S; potential density at constant S).
    Correlation (cited addition): Millero & Poisson, Deep-Sea Res. 28A, 625 (1981); UNESCO Tech. Pap. Mar. Sci. 44
    (Fofonoff & Millard 1983), Eq. (13), as implemented in the ``seawater`` package (``dens0``/``smow``).

    Parameters
    ----------
    T : float or array_like
        Temperature [K] (ITS-90; converted to IPTS-68 internally, t68 = 1.00024 t90).
    S : float or array_like
        Salinity [g/kg], used as practical salinity (PSS-78).

    Returns
    -------
    rho : float or ndarray
        [kg/m^3] at p = 0 dbar gauge.

    Notes
    -----
    The book defines S as kg of salt per kg of water; S [g/kg] = 1000 × that. Practical salinity is numerically close
    to g/kg (≈ 0.5 % difference), which is ignored here.

    Validation (planned): V5 UNESCO check values: S = 35, t68 = 0 °C → 1028.10633; S = 0, t68 = 30 °C → 995.65113.
    Label: pending.
    """
    # DEVIATION: salinity in g/kg is fed to EOS-80 as practical salinity (PSS-78); the two differ by about 0.5 %.
    s = np.asarray(S, dtype=float)
    t68 = 1.00024 * (np.asarray(T, dtype=float) - KELVIN_OFFSET)
    a = (999.842594, 6.793952e-2, -9.095290e-3, 1.001685e-4, -1.120083e-6, 6.536332e-9)
    b = (8.24493e-1, -4.0899e-3, 7.6438e-5, -8.2467e-7, 5.3875e-9)
    c = (-5.72466e-3, 1.0227e-4, -1.6546e-6)
    d = 4.8314e-4
    smow = a[0] + (a[1] + (a[2] + (a[3] + (a[4] + a[5] * t68) * t68) * t68) * t68) * t68
    rho = (smow + (b[0] + (b[1] + (b[2] + (b[3] + b[4] * t68) * t68) * t68) * t68) * s
           + (c[0] + (c[1] + c[2] * t68) * t68) * s * np.sqrt(np.maximum(s, 0.0)) + d * s ** 2)
    return as_scalar_if_0d(rho)


def seawater_linear_coefficients(T0=283.15, S0=35.0):
    """Coefficients of a linear seawater equation of state from the linearisation of EOS-80 at (T0, S0).

    Book: §1.10 (seawater; the thermal expansion coefficient (1.20) and its salinity analogue).

    Parameters
    ----------
    T0 : float, optional
        Reference temperature [K] (default 10 °C).
    S0 : float, optional
        Reference salinity [g/kg] (default 35).

    Returns
    -------
    (rho0, alpha_T, beta_S) : tuple of float
        rho0 [kg/m^3]; alpha_T = −(1/rho0) ∂rho/∂T [1/K]; beta_S = (1/rho0) ∂rho/∂S [1/(g/kg)]
        (at the defaults: 1026.95, 1.669e-4, 7.61e-4).

    Notes
    -----
    Method: central differences of :func:`seawater_density_eos80` (steps 1e-3 K and 1e-3 g/kg).

    Validation (planned): V1 the linear EOS matches EOS-80 to second order near (T0, S0). Label: pending.
    """
    rho0 = seawater_density_eos80(T0, S0)
    dT, dS = 1e-3, 1e-3
    drho_dT = (seawater_density_eos80(T0 + dT, S0) - seawater_density_eos80(T0 - dT, S0)) / (2 * dT)
    drho_dS = (seawater_density_eos80(T0, S0 + dS) - seawater_density_eos80(T0, S0 - dS)) / (2 * dS)
    return float(rho0), float(-drho_dT / rho0), float(drho_dS / rho0)


def seawater_density_linear(T, S, rho0=1027.0, alpha_T=1.67e-4, beta_S=7.6e-4, T0=283.15, S0=35.0):
    """Linear equation of state of seawater, ``rho = rho0 [1 − α_T (T − T0) + β_S (S − S0)]``.

    Book: §1.10 (seawater: rho = rho(T, p, S); a parcel keeps its salinity). Default coefficients are the rounded
    linearisation of EOS-80 (Millero & Poisson 1981) at 10 °C and 35 g/kg, see :func:`seawater_linear_coefficients`
    (1026.95 kg/m^3, 1.669e-4 1/K, 7.61e-4 per g/kg).

    Parameters
    ----------
    T : float or array_like
        Temperature [K].
    S : float or array_like
        Salinity [g/kg].
    rho0 : float, optional
        Reference density [kg/m^3].
    alpha_T : float, optional
        Thermal expansion coefficient [1/K].
    beta_S : float, optional
        Haline contraction coefficient [1/(g/kg)].
    T0 : float, optional
        Reference temperature [K].
    S0 : float, optional
        Reference salinity [g/kg].

    Returns
    -------
    rho : float or ndarray
        [kg/m^3] (surface pressure; compressibility enters separately through (1.35)).

    Validation (planned): V1 rho0 at (T0, S0); V5 within 0.5 kg/m^3 of EOS-80 over 5–15 °C, 33–37 g/kg. Label: pending.
    """
    return as_scalar_if_0d(rho0 * (1.0 - alpha_T * (np.asarray(T, dtype=float) - T0)
                                   + beta_S * (np.asarray(S, dtype=float) - S0)))


def fluid_properties(name: str, T: float = 293.15, p: float = P_ATM):
    """Density, dynamic and kinematic viscosity of air or water from the cited correlations in this module.

    Book: §1.5 (μ and ν of gases and liquids; Eq. (1.4)).

    Parameters
    ----------
    name : {"air", "water"}
        Fluid.
    T : float, optional
        Temperature [K].
    p : float, optional
        Pressure [Pa] (air only).

    Returns
    -------
    dict
        ``rho`` [kg/m^3], ``mu`` [Pa s], ``nu`` [m^2/s]; plus ``sigma`` [N/m] for water.

    Notes
    -----
    Air: perfect gas (1.22) with R_AIR and Sutherland μ. Water: Kell density, Vogel viscosity, IAPWS surface tension.

    Validation (planned): V5 air at 288.15 K matches USSA sea level (rho 1.2250, μ 1.7894e-5, ν 1.4607e-5). Label: pending.
    """
    name = name.lower()
    if name == "air":
        rho = perfect_gas_density(p, T)
        mu = sutherland_viscosity(T)
        return {"rho": rho, "mu": mu, "nu": kinematic_viscosity(mu, rho)}
    if name == "water":
        rho = water_density(T)
        mu = water_viscosity(T)
        return {"rho": rho, "mu": mu, "nu": kinematic_viscosity(mu, rho), "sigma": surface_tension_water(T)}
    raise ValueError("fluid_properties knows 'air' and 'water' (other fluids need a cited source)")


# Property rows used to build FLUIDS (T [K]: values). Sources: F. P. Incropera, D. P. DeWitt et al., *Fundamentals of
# Heat and Mass Transfer*, Appendix A — Table A.4 (air), A.5 (engine oil, glycerin), A.6 (saturated water),
# A.8 (binary diffusion coefficients); honey: FAO, *Value-added products from beekeeping*, Ch. 2 (Munro 1943 viscosity,
# White 1975 specific gravity) and the ranges reported in "Comparison of thermal and rheologic properties of Slovak
# mixed flower honey and forest honey" (k 0.43–0.45 W/(m K), c_p 2448–2575 J/(kg K)).
_T20 = 293.15
_AIR_A4 = {250.0: {"cp": 1006.0, "k": 22.3e-3}, 300.0: {"cp": 1007.0, "k": 26.3e-3}}
_WATER_A6 = {290.0: {"cp": 4184.0, "k": 598e-3}, 295.0: {"cp": 4181.0, "k": 606e-3}}
_OIL_A5 = {290.0: {"rho": 890.0, "cp": 1868.0, "mu": 0.999, "k": 0.145},
           300.0: {"rho": 884.1, "cp": 1909.0, "mu": 0.486, "k": 0.145}}
_GLYCERIN_A5 = {290.0: {"rho": 1265.8, "cp": 2367.0, "mu": 1.85, "k": 0.286},
                300.0: {"rho": 1259.9, "cp": 2427.0, "mu": 0.799, "k": 0.286}}
_D_H2O_AIR_298 = 0.26e-4       # Table A.8: H2O–air at 298 K [m^2/s]
_D_GLYCEROL_H2O_298 = 0.94e-9  # Table A.8: glycerol–H2O (dilute) at 298 K [m^2/s]


def _interp_rows(rows: dict, T: float, log_keys=("mu",)) -> dict:
    """Linear interpolation between two table rows (log-linear for viscosity)."""
    (T1, a), (T2, b) = sorted(rows.items())
    w = (T - T1) / (T2 - T1)
    out = {}
    for k in a:
        out[k] = float(np.exp((1 - w) * np.log(a[k]) + w * np.log(b[k]))) if k in log_keys else (1 - w) * a[k] + w * b[k]
    return out


def _build_fluids() -> dict[str, dict[str, float]]:
    air = _interp_rows(_AIR_A4, _T20)
    water = _interp_rows(_WATER_A6, _T20)
    oil = _interp_rows(_OIL_A5, _T20)
    gly = _interp_rows(_GLYCERIN_A5, _T20)
    mu_w = float(water_viscosity(_T20))
    fluids = {
        "air": {"rho": float(perfect_gas_density(P_ATM, _T20)), "mu": float(sutherland_viscosity(_T20)),
                "k": air["k"], "cp": air["cp"], "kappa_m": _D_H2O_AIR_298},
        "water": {"rho": float(water_density(_T20)), "mu": mu_w, "k": water["k"], "cp": water["cp"],
                  "kappa_m": _D_GLYCEROL_H2O_298},
        "glycerine": {"rho": gly["rho"], "mu": gly["mu"], "k": gly["k"], "cp": gly["cp"]},
        "honey": {"rho": 1.4237 * float(water_density(_T20)), "mu": 18.96, "k": 0.44, "cp": 2.51e3},
        "engine_oil": {"rho": oil["rho"], "mu": oil["mu"], "k": oil["k"], "cp": oil["cp"]},
    }
    # DEVIATION: no tabulated solute diffusivity exists for glycerine, honey or oil; estimated with the Stokes–Einstein
    # scaling D ∝ T/μ from glycerol in water (same solute, same temperature) — order of magnitude only.
    for name in ("glycerine", "honey", "engine_oil"):
        fluids[name]["kappa_m"] = _D_GLYCEROL_H2O_298 * mu_w / fluids[name]["mu"]
    for f in fluids.values():
        f["nu"] = f["mu"] / f["rho"]  # Eq. (1.4)
    return fluids


#: Properties at 20 °C (293.15 K), 1 atm: ``rho`` [kg/m^3], ``mu`` [Pa s], ``nu`` [m^2/s], ``k`` [W/(m K)],
#: ``cp`` [J/(kg K)], ``kappa_m`` [m^2/s] (diffusivity of a small solute/vapour: water vapour in air; glycerol in water;
#: Stokes–Einstein estimates for the viscous liquids). Air: perfect gas (1.22) + Sutherland; water: Kell density +
#: Vogel viscosity; k, cp and the other liquids from Incropera Appendix A (interpolated to 293.15 K, μ log-linearly);
#: honey (≈ 17 % water): FAO/Munro/White values. Not the book's numbers.
FLUIDS: dict[str, dict[str, float]] = {}


# ====================================================================================================================
# §1.10 Synthetic atmospheric profile (our data, not digitised from the book)
# ====================================================================================================================
def synthetic_boundary_layer_profile(z, T_surface=288.15, mixed_top=800.0, inversion_depth=200.0,
                                     inversion_dT_dz=0.01, mixed_dT_dz=-0.0098, upper_dT_dz=-0.0045):
    """Piecewise-linear temperature profile of a lower atmosphere: near-neutral mixed layer, inversion, stable layer.

    Book: §1.10 and Fig. 1.9a (a near-neutral layer by the ground with slope ≈ Γ_a, an inversion where T increases with
    height, and a layer above cooling more slowly than Γ_a). The numbers are ours (not digitised from the book).

    Parameters
    ----------
    z : float or array_like
        Height above ground [m], >= 0.
    T_surface : float, optional
        Surface temperature [K].
    mixed_top : float, optional
        Top of the mixed layer [m].
    inversion_depth : float, optional
        Thickness of the inversion layer [m].
    inversion_dT_dz : float, optional
        Γ = dT/dz inside the inversion [K/m] (positive: temperature rises with height).
    mixed_dT_dz : float, optional
        Γ in the mixed layer [K/m] (default −9.8e-3, marginally steeper than Γ_a = −9.76e-3: near-neutral).
    upper_dT_dz : float, optional
        Γ above the inversion [K/m] (default −4.5e-3 > Γ_a: stable).

    Returns
    -------
    T : float or ndarray
        [K], continuous at the layer boundaries.

    Validation (planned): V1 slopes per layer equal the inputs; continuity at the kinks. Label: pending.
    """
    z = np.asarray(z, dtype=float)
    z_edges, slopes, T_edges = _bl_layers(T_surface, mixed_top, inversion_depth, inversion_dT_dz, mixed_dT_dz,
                                          upper_dT_dz)
    k = np.clip(np.searchsorted(z_edges, z, side="right") - 1, 0, 2)
    return as_scalar_if_0d(T_edges[k] + slopes[k] * (z - z_edges[k]))


def _bl_layers(T_surface, mixed_top, inversion_depth, inversion_dT_dz, mixed_dT_dz, upper_dT_dz):
    z_edges = np.array([0.0, mixed_top, mixed_top + inversion_depth])
    slopes = np.array([mixed_dT_dz, inversion_dT_dz, upper_dT_dz])
    T_edges = T_surface + np.concatenate(([0.0], np.cumsum(slopes[:-1] * np.diff(z_edges))))
    return z_edges, slopes, T_edges


def synthetic_boundary_layer_column(z, T_surface=288.15, mixed_top=800.0, inversion_depth=200.0,
                                    inversion_dT_dz=0.01, mixed_dT_dz=-0.0098, upper_dT_dz=-0.0045,
                                    p_surface=P_ATM, p_ref=P_REF, cp=CP_AIR, R=R_AIR, g=G0, gamma=GAMMA_AIR):
    """Complete static column for :func:`synthetic_boundary_layer_profile`: T, Γ, p, rho, θ, rho_θ, N^2 and stability.

    Book: §1.10 — hydrostatics (1.8) with (1.22) layer by layer (closed form, :func:`linear_lapse_pressure`),
    potential temperature (1.31) and (1.32), potential density (1.33), N^2 = (g/θ) dθ/dz, classification after (1.29).

    Parameters
    ----------
    z : array_like
        Heights [m], >= 0.
    T_surface, mixed_top, inversion_depth, inversion_dT_dz, mixed_dT_dz, upper_dT_dz : float, optional
        Profile parameters (see :func:`synthetic_boundary_layer_profile`).
    p_surface : float, optional
        Surface pressure p0 [Pa].
    p_ref : float, optional
        Reference pressure for θ and rho_θ [Pa].
    cp, R, g, gamma : float, optional
        Gas properties [J/(kg K)], [J/(kg K)], gravity [m/s^2], γ (defaults satisfy (γ−1)/γ = R/C_p).

    Returns
    -------
    dict of ndarray
        ``z`` [m], ``T`` [K], ``dT_dz`` [K/m], ``p`` [Pa], ``rho`` [kg/m^3], ``theta`` [K], ``dtheta_dz`` [K/m],
        ``rho_theta`` [kg/m^3], ``N2`` [1/s^2], ``stability`` (labels, tolerance 1e-8 1/s^2).

    Validation (planned): V1 θ rho_θ = p_ref/R; V4 dθ/dz from (1.32) equals finite differences of θ. Label: pending.
    """
    z = np.atleast_1d(np.asarray(z, dtype=float))
    z_edges, slopes, T_edges = _bl_layers(T_surface, mixed_top, inversion_depth, inversion_dT_dz, mixed_dT_dz,
                                          upper_dT_dz)
    p_edges = [p_surface]
    for i in range(2):
        p_edges.append(linear_lapse_pressure(z_edges[i + 1] - z_edges[i], p_edges[-1], T_edges[i], slopes[i], R, g))
    k = np.clip(np.searchsorted(z_edges, z, side="right") - 1, 0, 2)
    T = T_edges[k] + slopes[k] * (z - z_edges[k])
    dT_dz = slopes[k]
    p = np.array([linear_lapse_pressure(zz - z_edges[kk], p_edges[kk], T_edges[kk], slopes[kk], R, g)
                  for zz, kk in zip(z, k)])  # (1.8) + (1.22), layer by layer
    rho = p / (R * T)  # Eq. (1.22)
    theta = potential_temperature(T, p, p_ref, gamma)  # Eq. (1.31)
    dtheta_dz = potential_temperature_gradient(T, dT_dz, theta, cp, g)  # Eq. (1.32)
    N2 = brunt_vaisala_sq_from_theta(theta, dtheta_dz, g)
    return {"z": z, "T": T, "dT_dz": dT_dz, "p": p, "rho": rho, "theta": theta, "dtheta_dz": dtheta_dz,
            "rho_theta": potential_density(rho, p, p_ref, gamma), "N2": N2,
            "stability": np.asarray(classify_stability(N2, tol=1e-8))}


# ====================================================================================================================
# §1.11 Dimensional analysis: examples
# ====================================================================================================================
def poiseuille_pressure_drop(mu, U, dx, d):
    """Laminar pressure drop in a round pipe (Hagen–Poiseuille), used as a data generator for the pipe Π collapse.

    Book: §1.11 pipe example (1.38)–(1.40); the laminar result itself is derived in Ch. 8:
    ``Δp = 32 μ U Δx / d^2`` with U the mean velocity, i.e. Π1 = 32 (Δx/d)(μ/(rho U d)).

    Parameters
    ----------
    mu : float or array_like
        Dynamic viscosity [Pa s].
    U : float or array_like
        Mean velocity [m/s].
    dx : float or array_like
        Distance between pressure taps [m].
    d : float or array_like
        Pipe diameter [m].

    Returns
    -------
    dp : float or ndarray
        Pressure drop [Pa].

    Notes
    -----
    Assumptions: steady, fully developed laminar flow (Re ≲ 2000), smooth pipe (ε irrelevant).

    Validation (planned): V2 dimensionless form Π1 = 32 Π2 Π4. Label: pending.
    """
    return as_scalar_if_0d(32.0 * np.asarray(mu, dtype=float) * np.asarray(U, dtype=float)
                           * np.asarray(dx, dtype=float) / np.asarray(d, dtype=float) ** 2)


def pythagoras_phi(beta):
    """The shape function of Example 1.3: area of a right triangle with unit hypotenuse and acute angle β.

    Book: §1.11, Example 1.3 (a = C^2 φ(β)); explicit form (our check): legs C cos β and C sin β give
    ``φ(β) = (1/2) sin β cos β = (1/4) sin 2β``.

    Parameters
    ----------
    beta : float or array_like
        Acute angle [rad].

    Returns
    -------
    phi : float or ndarray
        [-].

    Validation (planned): V1 A^2 φ + B^2 φ = C^2 φ for random triangles. Label: pending.
    """
    return as_scalar_if_0d(0.25 * np.sin(2.0 * np.asarray(beta, dtype=float)))


#: Taylor's similarity constant for a spherical blast in a gas with γ = 1.4: K = S(1.4)^-5 = 0.856 in E = K rho D^5/t^2
#: (G. I. Taylor, Proc. R. Soc. Lond. A 201, 159 (1950), as reported by J. S. Díaz, arXiv:2009.05674).
TAYLOR_K_GAMMA14: float = 0.856


def blast_energy(D, t, rho, K=1.0):
    """Energy of an intense point blast from the blast-wave radius at a given time.

    Book: §1.11, Example 1.4: the lone group Π1 = E t^2/(rho D^5) is a constant K, so ``E = K rho D^5 / t^2``.

    Parameters
    ----------
    D : float or array_like
        Blast-wave radius [m].
    t : float or array_like
        Time since the release [s].
    rho : float
        Undisturbed air density [kg/m^3].
    K : float, optional
        Dimensionless constant, not given by dimensional analysis (default 1.0; Taylor's full solution gives
        ``TAYLOR_K_GAMMA14`` = 0.856 for a spherical blast in γ = 1.4 gas).

    Returns
    -------
    E : float or ndarray
        [J].

    Notes
    -----
    Assumptions: strong blast (ambient pressure negligible), γ fixed (hidden in K), point release.

    Validation (planned): V1 round trip with :func:`blast_radius`; D ∝ t^(2/5) (log–log slope 0.4). Label: pending.
    """
    return as_scalar_if_0d(K * rho * np.asarray(D, dtype=float) ** 5 / np.asarray(t, dtype=float) ** 2)  # Example 1.4


def blast_radius(E, t, rho, K=1.0):
    """Blast-wave radius at time t for energy E: ``D = (E t^2 / (K rho))^(1/5)``.

    Book: §1.11, Example 1.4 inverted (the testable form D ∝ t^(2/5)).

    Parameters
    ----------
    E : float or array_like
        Released energy [J].
    t : float or array_like
        Time [s].
    rho : float
        Air density [kg/m^3].
    K : float, optional
        Dimensionless constant (see :func:`blast_energy`).

    Returns
    -------
    D : float or ndarray
        [m].

    Validation (planned): V1 blast_energy(blast_radius(E, t)) = E. Label: pending.
    """
    return as_scalar_if_0d((np.asarray(E, dtype=float) * np.asarray(t, dtype=float) ** 2 / (K * rho)) ** 0.2)


def rayleigh_scattering_ratio(V, d, lam, phi3=1.0):
    """Scattered-to-incident intensity ratio for a particle much smaller than the wavelength (Rayleigh law).

    Book: §1.11, Example 1.5: ``S/I = (λ/d)^2 (V/λ^3)^2 φ3(n_s) = V^2 φ3(n_s) / (d^2 λ^4)``.

    Parameters
    ----------
    V : float or array_like
        Particle volume [m^3].
    d : float or array_like
        Distance from the particle to the observer [m].
    lam : float or array_like
        Wavelength [m].
    phi3 : float, optional
        The undetermined function of the refractive index φ3(n_s) [-] (default 1: only ratios are meaningful).

    Returns
    -------
    ratio : float or ndarray
        S/I [-].

    Notes
    -----
    Assumptions: λ >> V^(1/3) (dipole scattering), far field (energy conservation gives S ∝ d^-2), elastic scattering.

    Validation (planned): V1 blue (450 nm)/red (700 nm) = (700/450)^4 ≈ 5.85; S ∝ V^2, ∝ d^-2; group dimensionless.
    Label: pending.
    """
    lam = np.asarray(lam, dtype=float)
    return as_scalar_if_0d(phi3 * np.asarray(V, dtype=float) ** 2
                           / (np.asarray(d, dtype=float) ** 2 * lam ** 4))  # Example 1.5


def wavelength_to_rgb(lam_nm):
    """Approximate display colour of monochromatic light (visual aid only, not physics).

    Piecewise-linear approximation of the visible spectrum (after Dan Bruton's widely used algorithm), with intensity
    fall-off near the ends of the visible range.

    Parameters
    ----------
    lam_nm : float or array_like
        Wavelength [nm], 380–780 (outside → black).

    Returns
    -------
    (r, g, b) : tuple of float or ndarray
        Colour components in [0, 1].

    Validation: visual only (not evidence). Label: qualitative.
    """
    w = np.atleast_1d(np.asarray(lam_nm, dtype=float))
    r = np.zeros_like(w)
    gr = np.zeros_like(w)
    b = np.zeros_like(w)
    m = (w >= 380) & (w < 440)
    r[m], b[m] = -(w[m] - 440) / 60.0, 1.0
    m = (w >= 440) & (w < 490)
    gr[m], b[m] = (w[m] - 440) / 50.0, 1.0
    m = (w >= 490) & (w < 510)
    gr[m], b[m] = 1.0, -(w[m] - 510) / 20.0
    m = (w >= 510) & (w < 580)
    r[m], gr[m] = (w[m] - 510) / 70.0, 1.0
    m = (w >= 580) & (w < 645)
    r[m], gr[m] = 1.0, -(w[m] - 645) / 65.0
    m = (w >= 645) & (w <= 780)
    r[m] = 1.0
    fade = np.where(w < 420, 0.3 + 0.7 * (w - 380) / 40.0, np.where(w > 700, 0.3 + 0.7 * (780 - w) / 80.0, 1.0))
    fade = np.where((w < 380) | (w > 780), 0.0, fade)
    out = [np.clip(c * fade, 0.0, 1.0) for c in (r, gr, b)]
    if np.ndim(lam_nm) == 0:
        return tuple(float(c[0]) for c in out)
    return tuple(c.reshape(np.shape(lam_nm)) for c in out)


FLUIDS.update(_build_fluids())

_NOT_EXPORTED = {"annotations", "Mapping", "Sequence", "np", "as_scalar_if_0d", "require_nonnegative",
                 "require_positive"}
__all__ = sorted(name for name in dir() if not name.startswith("_") and name not in _NOT_EXPORTED)
