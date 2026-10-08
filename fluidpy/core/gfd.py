"""Geophysical fluid dynamics: planetary parameters, geostrophic and thermal-wind balance, Ekman layers, the dispersion
relations of rotating shallow water and of rotating stratified fluid, Rossby radii, potential vorticity, the Eady problem,
and the two-dimensional cascade.

Book: Kundu, Cohen & Dowling 5e, Ch. 13 §13.3–§13.18, Eqs. (13.5)–(13.145) (every equation coded here was read from the
rendered pages chapters/pages/ch13/p653–p715).

Conventions (``knowledge/notation.md``)
---------------------------------------
* x east, y north, z up; (u, v, w).  SI units.  Latitude in radians (``lat_rad``); degrees only in ``*_deg`` names.
* **Sign of f.**  f = 2Ω sin(latitude) is positive in the northern hemisphere and negative in the southern.  The book
  writes its formulas for f > 0 (``sqrt(f)``, "to the right", "clockwise", ``exp(-f y/c)``).  Every function here
  accepts either sign: it uses ``abs(f)`` for scales and ``sign(f)`` for directions.  Functions with f in a denominator
  raise ``ValueError`` at f = 0 (the equator) instead of returning infinity.
* **Which z = 0** is stated per function: the sea surface with the ocean below (surface Ekman layer), the solid surface
  (bottom Ekman layer), the lower lid (Eady problem).
* Complex horizontal velocity V = u + i v (book §13.6); complex stress tau = tau_x + i tau_y.
* The rotation rate is the sidereal ``OMEGA_EARTH`` of ``core.rotating``; the book rounds one turn per solar day
  (``OMEGA_SOLAR_DAY``, 0.27 % smaller).
* Scalars in, Python floats (or tuples / dicts of floats) out, so explainer parity rows can call every closed form.
* Items marked "ours — not in the book" are standard results the book states only in words or not at all (thermal wind
  in temperature form, Ekman pumping, Sverdrup balance, an Ekman layer with K(z), geostrophic adjustment, group
  velocities, the Eady growth rate in physical units).
"""
from __future__ import annotations

from typing import Callable

import numpy as np
from scipy.linalg import solve_banded
from scipy.optimize import brentq, minimize_scalar

from ._util import as_scalar_if_0d
from .rotating import OMEGA_EARTH, coriolis_parameter
from .thermo import G0
from .waves import internal_wave_omega

__all__ = [
    # constants and geometry
    "OMEGA_EARTH", "OMEGA_SOLAR_DAY", "EARTH_RADIUS_MEAN", "SIDEREAL_DAY", "coriolis_parameter",
    "coriolis_parameter_deg", "hemisphere", "earth_rotation_local", "inertial_period", "f_plane", "beta_parameter",
    "beta_plane", "beta_plane_error", "vertical_velocity_scale", "coriolis_acceleration_local", "thin_layer_terms",
    "thin_layer_residuals", "eddy_friction", "scale_height", "buoyancy_frequency_sq",
    # balances
    "geostrophic_velocity", "geostrophic_from_field", "geostrophic_from_height", "geostrophic_streamfunction",
    "rossby_number", "ekman_number", "thermal_wind_shear", "thermal_wind_from_temperature", "thermal_wind_integrate",
    "taylor_proudman_residual", "parcel_adjust",
    # Ekman layers
    "ekman_depth", "ekman_surface", "ekman_transport", "ekman_transport_partial", "ekman_residual",
    "ekman_vorticity_balance", "ekman_solve", "ekman_pumping", "ekman_pumping_from_curl", "sverdrup_transport",
    "ekman_pumping_bottom", "ekman_finite_depth", "ekman_bottom", "ekman_bottom_transport", "ekman_force_balance",
    "eddy_viscosity_from_depth",
    # shallow-water waves
    "long_wave_speed", "equivalent_depth", "baroclinic_mode_speed", "shallow_water_omega", "shallow_water_branches",
    "shallow_water_discriminant", "shallow_water_regime", "dispersion_term_sizes", "poincare_omega",
    "poincare_group_velocity", "poincare_amplitudes", "poincare_fields", "poincare_orbit", "inertial_oscillation",
    "inertial_radius", "kelvin_wave", "kelvin_omega", "kelvin_decay_side", "kelvin_residuals", "rossby_radius",
    "rossby_radius_internal", "rossby_radius_two_layer", "geostrophic_adjustment_1d", "adjustment_energy",
    # internal waves with rotation
    "inertia_gravity_m2", "inertia_gravity_band", "inertia_gravity_omega", "inertia_gravity_regime",
    "inertia_gravity_group_velocity", "wkb_vertical_structure", "inertia_gravity_fields", "inertia_gravity_hodograph",
    "lee_wave_m",
    # Rossby waves
    "rossby_omega", "rossby_group_velocity", "rossby_phase_speed", "rossby_max_frequency", "rossby_omega_circle",
    "rossby_long_wave_speed", "stationary_rossby_wavelength", "rossby_packet_spectrum",
    # vorticity and instability
    "potential_vorticity", "step_vorticity", "absolute_vorticity_gradient", "rayleigh_kuo_criterion", "eady_alpha",
    "eady_phase_speed", "eady_growth_rate", "eady_factors", "eady_critical", "eady_fastest", "eady_wavelengths",
    "eady_max_growth_rate",
    "eady_time_scale", "eady_mode", "eady_vertical_velocity", "eady_fluxes",
    # two-dimensional turbulence
    "fjortoft_transfer", "enstrophy_spectrum", "two_d_cascade_spectrum", "rhines_length",
]

OMEGA_SOLAR_DAY: float = 2.0 * np.pi / 86400.0   #: one turn per solar day [rad/s] — the book's rounding (trap T1)
EARTH_RADIUS_MEAN: float = 6.371e6               #: mean radius of the Earth [m] (the sphere used for beta)
SIDEREAL_DAY: float = 2.0 * np.pi / OMEGA_EARTH  #: one turn relative to the stars [s] (about 86 164 s)

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d


def _nonzero_f(f, who: str):
    """Return f as a float array; raise at f = 0 (the equator), where every f-in-the-denominator formula fails."""
    fa = _F(f)
    if np.any(fa == 0.0):
        raise ValueError(f"{who}: f = 0 (the equator) — this balance divides by the Coriolis parameter and does not "
                         "hold there")
    return fa


def _tup(*vals):
    return tuple(_S(v) for v in vals)


def _d1(z: np.ndarray, y: np.ndarray) -> np.ndarray:
    """First derivative on a non-uniform 1-D grid: 3-point Lagrange stencils, second order at every node (one-sided at
    the two ends).  Explicit stencils — not ``np.gradient``, which is first order at the edges."""
    z = _F(z)
    y = np.asarray(y)
    n = z.size
    if n < 3:
        raise ValueError("need at least 3 nodes for a second-order derivative")
    d = np.empty_like(y, dtype=np.result_type(y.dtype, float))
    h0, h1 = z[1:-1] - z[:-2], z[2:] - z[1:-1]
    d[1:-1] = (-h1 / (h0 * (h0 + h1)) * y[:-2] + (h1 - h0) / (h0 * h1) * y[1:-1] + h0 / (h1 * (h0 + h1)) * y[2:])
    a, b = z[1] - z[0], z[2] - z[0]
    d[0] = -(a + b) / (a * b) * y[0] + b / (a * (b - a)) * y[1] - a / (b * (b - a)) * y[2]
    a, b = z[-1] - z[-2], z[-1] - z[-3]
    d[-1] = (a + b) / (a * b) * y[-1] - b / (a * (b - a)) * y[-2] + a / (b * (b - a)) * y[-3]
    return d


def _d2(z: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Second derivative at the interior nodes of a non-uniform grid (3-point; second order on a uniform or smoothly
    stretched grid).  The two end values are NaN."""
    z = _F(z)
    y = np.asarray(y)
    d = np.full(y.shape, np.nan, dtype=np.result_type(y.dtype, float))
    h0, h1 = z[1:-1] - z[:-2], z[2:] - z[1:-1]
    d[1:-1] = 2.0 * (y[:-2] / (h0 * (h0 + h1)) - y[1:-1] / (h0 * h1) + y[2:] / (h1 * (h0 + h1)))
    return d


def _ddx_grid(a: np.ndarray, d: float, axis: int, periodic: bool = False) -> np.ndarray:
    """∂/∂x on a uniform grid: centred differences inside, second-order one-sided stencils at the two edges (or the
    periodic wrap)."""
    a = np.asarray(a, dtype=float)
    if periodic:
        return (np.roll(a, -1, axis) - np.roll(a, 1, axis)) / (2.0 * d)
    a = np.moveaxis(a, axis, -1)
    if a.shape[-1] < 3:
        raise ValueError("need at least 3 points along each differentiated axis")
    out = np.empty_like(a)
    out[..., 1:-1] = (a[..., 2:] - a[..., :-2]) / (2.0 * d)
    out[..., 0] = (-3.0 * a[..., 0] + 4.0 * a[..., 1] - a[..., 2]) / (2.0 * d)
    out[..., -1] = (3.0 * a[..., -1] - 4.0 * a[..., -2] + a[..., -3]) / (2.0 * d)
    return np.moveaxis(out, -1, axis)


# =====================================================================================================================
# 1. Planetary parameters and the thin-layer equations (§13.3, §13.4)
# =====================================================================================================================

def coriolis_parameter_deg(lat_deg, Omega: float = OMEGA_EARTH):
    """Coriolis parameter f = 2Ω sin(latitude) for a latitude in degrees — interface wrapper of
    :func:`fluidpy.core.rotating.coriolis_parameter` (the single definition of f).

    Book: §13.4, Eq. (13.8).
    Parameters
    ----------
    lat_deg : latitude [degrees], positive north.   Omega : rotation rate [rad/s] (default sidereal).
    Returns
    -------
    f : Coriolis parameter [1/s]; negative in the southern hemisphere.
    Assumptions: spherical Earth, local vertical component of 2Ω only.
    Validation: V1 pole = 2Ω, equator = 0, odd in latitude.  Label: analytic.
    """
    return coriolis_parameter(np.deg2rad(_F(lat_deg)), Omega)


def hemisphere(lat_rad) -> dict:
    """Which way rotation deflects at this latitude: the sign of f and the words that go with it.

    Book: §13.1 (all statements are for the northern hemisphere; the sense reverses in the southern) and §13.4,
    Eq. (13.8).
    Parameters
    ----------
    lat_rad : latitude [rad] (scalar).
    Returns
    -------
    dict: ``sign`` (+1.0, −1.0 or 0.0 — the sign of f), ``name`` ("northern", "southern", "equator"), ``turns``
    ("right", "left", "none" — the side toward which the Coriolis force turns a moving parcel), ``cyclonic``
    ("counter-clockwise", "clockwise", "none" — the sense of flow round a low).
    Assumptions: none.   Validation: V7 flips with the sign of the latitude.  Label: analytic.
    """
    s = float(np.sign(np.sin(float(lat_rad))))
    if s > 0:
        return dict(sign=1.0, name="northern", turns="right", cyclonic="counter-clockwise")
    if s < 0:
        return dict(sign=-1.0, name="southern", turns="left", cyclonic="clockwise")
    return dict(sign=0.0, name="equator", turns="none", cyclonic="none")



def earth_rotation_local(lat_rad, Omega: float = OMEGA_EARTH):
    """Components of the Earth's rotation vector in the local tangent-plane axes (x east, y north, z up).

    Book: §13.4, the unnumbered line before Eq. (13.7): Ω = (0, Ω cos θ, Ω sin θ), θ the latitude.
    Parameters
    ----------
    lat_rad : latitude [rad] (scalar).   Omega : rotation rate [rad/s].
    Returns
    -------
    ndarray (3,) = (Omega_x, Omega_y, Omega_z) [rad/s].
    Assumptions: spherical Earth.   Validation: V1 Omega_z = f/2; magnitude = Omega.  Label: analytic.
    """
    lat = float(lat_rad)
    return np.array([0.0, Omega * np.cos(lat), Omega * np.sin(lat)])



def inertial_period(f):
    """Inertial period T_i = 2π/|f| — the time a parcel takes to go once round an inertial circle.

    Book: §13.4 (unnumbered, after Eq. (13.8)); §13.11 "Inertial Motion".  (The subscript i is a label for "inertial",
    not a vector component — slip #5.)
    Parameters
    ----------
    f : Coriolis parameter [1/s], either sign.
    Returns
    -------
    T_i [s]; ``inf`` at f = 0 (no inertial oscillation at the equator).
    Assumptions: none.   Validation: V1 2π/|f|; even in f.  Label: analytic.
    """
    fa = np.abs(_F(f))
    with np.errstate(divide="ignore"):
        return _S(np.where(fa > 0, 2.0 * np.pi / np.where(fa > 0, fa, 1.0), np.inf))


def f_plane(lat0_rad, Omega: float = OMEGA_EARTH):
    """f-plane: the constant Coriolis parameter f₀ = 2Ω sin θ₀ of a region centred on latitude θ₀.

    Book: §13.4 "f-Plane Model" (unnumbered).
    Parameters: lat0_rad central latitude [rad]; Omega [rad/s].   Returns f0 [1/s].
    Assumptions: region small enough that f hardly varies (time scales below weeks, lengths below thousands of km).
    Validation: V1 equals :func:`coriolis_parameter`.  Label: analytic.
    """
    return coriolis_parameter(lat0_rad, Omega)


def beta_parameter(lat_rad, Omega: float = OMEGA_EARTH, R: float = EARTH_RADIUS_MEAN):
    """Northward gradient of the Coriolis parameter, β = df/dy = 2Ω cos θ₀ / R.

    Book: §13.4, the definition under Eq. (13.10) (dθ/dy = 1/R).
    Parameters
    ----------
    lat_rad : latitude θ₀ [rad].   Omega : rotation rate [rad/s].   R : planetary radius [m].
    Returns
    -------
    beta [1/(m s)]; positive in both hemispheres, largest at the equator, zero at the poles.
    Assumptions: spherical Earth.
    Validation: V1 complex-step derivative of 2Ω sin(θ₀ + y/R); V7 even in latitude.  Label: analytic.
    """
    return _S(2.0 * Omega * np.cos(_F(lat_rad)) / R)  # Eq. (13.10), definition of beta


def beta_plane(y, lat0_rad, Omega: float = OMEGA_EARTH, R: float = EARTH_RADIUS_MEAN):
    """β-plane Coriolis parameter f = f₀ + β y.

    Book: §13.4, Eq. (13.10).
    Parameters
    ----------
    y : northward distance from the central latitude [m].   lat0_rad : central latitude θ₀ [rad].
    Omega [rad/s], R [m] : planet.
    Returns
    -------
    f [1/s].
    Assumptions: first-order Taylor expansion about θ₀ (|y| ≪ R).
    Validation: V3 error against the sphere is second order in y/R (:func:`beta_plane_error`).  Label: analytic.
    """
    return _S(coriolis_parameter(lat0_rad, Omega) + beta_parameter(lat0_rad, Omega, R) * _F(y))  # Eq. (13.10)


def beta_plane_error(y, lat0_rad, Omega: float = OMEGA_EARTH, R: float = EARTH_RADIUS_MEAN):
    """Relative error of the β-plane against the sphere: (f₀ + βy − f_exact)/f_exact, f_exact = 2Ω sin(θ₀ + y/R).

    Book: §13.4, Eq. (13.10) (the truncation of the Taylor series — the size of the error is ours).
    Parameters: as :func:`beta_plane`.
    Returns the relative error (dimensionless); leading term +½ (y/R)² for small y/R (the β-plane overestimates the
    magnitude of f on both sides of θ₀).
    Raises ValueError where f_exact = 0 (the point lies on the equator).
    Assumptions: spherical Earth.   Validation: V3 observed order 2 in y/R.  Label: analytic.
    """
    exact = _nonzero_f(2.0 * Omega * np.sin(float(lat0_rad) + _F(y) / R), "beta_plane_error")
    return _S((_F(beta_plane(y, lat0_rad, Omega, R)) - exact) / exact)



def vertical_velocity_scale(U, H, L):
    """Thin-layer scale of the vertical velocity from continuity: W ~ U H / L.

    Book: §13.4 (unnumbered scaling W/U ~ H/L).
    Parameters: U horizontal velocity scale [m/s]; H depth scale [m]; L horizontal length scale [m].
    Returns W [m/s] (an order of magnitude, not a value).
    Assumptions: incompressible, all three terms of continuity of the same order.
    Validation: V1 identity.  Label: analytic.
    """
    return _S(_F(U) * _F(H) / _F(L))


def coriolis_acceleration_local(u, v, w, lat_rad, Omega: float = OMEGA_EARTH, thin: bool = True):
    """Coriolis acceleration 2Ω × u in local axes: (−f v, f u, −2Ω u cos θ) for a thin layer.

    Book: §13.4, Eq. (13.7) and the determinant above it.  With ``thin=False`` the x-component keeps the term the book
    drops, 2Ω(w cos θ − v sin θ).
    Parameters
    ----------
    u, v, w : velocity components east, north, up [m/s].   lat_rad : latitude [rad].   Omega [rad/s].
    thin : drop w cos θ against v sin θ (the thin-layer approximation) if True.
    Returns
    -------
    (a_x, a_y, a_z) [m/s²] — the *acceleration* terms on the left of the momentum equation (the Coriolis force per unit
    mass is minus this).
    Assumptions: spherical Earth; ``thin=True`` needs w ≪ v.
    Validation: V1 equals ``core.rotating.coriolis_acceleration`` with the full Ω vector when thin=False; the
    difference with thin=True is exactly 2Ω w cos θ.  Label: analytic.
    """
    u, v, w, lat = _F(u), _F(v), _F(w), _F(lat_rad)
    f = 2.0 * Omega * np.sin(lat)
    ax = -f * v if thin else 2.0 * Omega * (w * np.cos(lat) - v * np.sin(lat))
    return _tup(ax, f * u, -2.0 * Omega * u * np.cos(lat))  # Eq. (13.7)


def eddy_friction(lap_h, d2z, nu_H, nu_v):
    """Friction force per unit mass with anisotropic eddy viscosities: F = ν_H ∇_H² q + ν_v ∂²q/∂z² for q = u, v or w.

    Book: §13.3, Eq. (13.6) — in the corrected reading F_i = (1/ρ) ∂τ_ij/∂x_j (the book omits the 1/ρ, slip #2).
    Parameters
    ----------
    lap_h : horizontal Laplacian ∂²q/∂x² + ∂²q/∂y² [1/(m s)].   d2z : ∂²q/∂z² [1/(m s)].
    nu_H, nu_v : horizontal and vertical eddy viscosities [m²/s].
    Returns
    -------
    F [m/s²].
    Assumptions: ∇·u = 0 and uniform ν_H, ν_v (both used silently by the book to cancel the cross terms).
    Validation: V2 sympy: divergence of (13.5)/ρ minus ν ∂(∇·u) (``ch13.eddy_friction_force_sympy``).  Label: symbolic.
    """
    return _S(_F(nu_H) * _F(lap_h) + _F(nu_v) * _F(d2z))  # Eq. (13.6)


def thin_layer_terms(U, L, H, lat_rad, *, drho_over_rho0: float = 1e-3, nu_H: float = 0.0, nu_v: float = 0.0,
                     g: float = G0, Omega: float = OMEGA_EARTH) -> dict:
    """Order-of-magnitude size of every term of the thin-shell momentum equations for a choice of scales.

    Book: §13.4, Eq. (13.9) with the scaling W/U ~ H/L and the full Coriolis acceleration before Eq. (13.7); §13.5,
    Eq. (13.13) for Ro.  (The scale analysis as a table is ours.)
    Parameters
    ----------
    U : horizontal velocity scale [m/s].   L : horizontal length scale [m].   H : vertical scale [m].
    lat_rad : latitude [rad] (not the equator).
    drho_over_rho0 : relative density perturbation ρ′/ρ₀ (dimensionless; an input you must judge — default 1e-3).
    nu_H, nu_v : eddy viscosities [m²/s].   g [m/s²].   Omega [rad/s].
    Returns
    -------
    dict: ``W`` = UH/L [m/s]; ``aspect`` = H/L; ``f`` [1/s] (signed); ``Ro`` = U/(|f|L); ``x`` → dict of sizes [m/s²]:
    "acceleration" U²/L, "coriolis" |f|U, "coriolis_w" 2Ω cos θ·W (the term the thin-layer approximation drops),
    "friction_H" ν_H U/L², "friction_v" ν_v U/H², "pressure" (taken equal to the largest of the others: the pressure
    gradient adjusts to balance it); ``z`` → dict: "acceleration" UW/L, "coriolis" 2ΩU cos θ, "buoyancy"
    g·drho_over_rho0, "friction_H" ν_H W/L², "friction_v" ν_v W/H².
    Raises ValueError at the equator (Ro divides by f).
    Assumptions: one scale per variable; an estimate, not a solution.
    Validation: V1 identities; coriolis_w/coriolis = cot θ · H/L.  Label: analytic.
    """
    U, L, H, lat = float(U), float(L), float(H), float(lat_rad)
    f = float(_nonzero_f(2.0 * Omega * np.sin(lat), "thin_layer_terms"))
    W = U * H / L
    x = dict(acceleration=U * U / L, coriolis=abs(f) * U, coriolis_w=float(2.0 * Omega * abs(np.cos(lat)) * W),
             friction_H=nu_H * U / L ** 2, friction_v=nu_v * U / H ** 2)
    x["pressure"] = max(x.values())
    z = dict(acceleration=U * W / L, coriolis=float(2.0 * Omega * U * abs(np.cos(lat))), buoyancy=g * drho_over_rho0,
             friction_H=nu_H * W / L ** 2, friction_v=nu_v * W / H ** 2)
    return dict(W=W, aspect=H / L, f=f, Ro=U / (abs(f) * L), x=x, z=z)



def thin_layer_residuals(*, u=0.0, v=0.0, f, dpdx=0.0, dpdy=0.0, dpdz=0.0, rho=0.0, rho0, g: float = G0, Du_Dt=0.0,
                         Dv_Dt=0.0, Dw_Dt=0.0, lap_h_u=0.0, lap_h_v=0.0, lap_h_w=0.0, d2u_dz2=0.0, d2v_dz2=0.0,
                         d2w_dz2=0.0, nu_H=0.0, nu_v=0.0) -> dict:
    """Every term of the thin-shell momentum equations at a point, named, with all terms moved to the left so each
    row sums to its residual.

    Book: §13.4, Eq. (13.9) (primes dropped: p and rho are perturbations from the state of rest).
    Parameters (keyword-only)
    -------------------------
    u, v [m/s]; f [1/s]; dpdx, dpdy, dpdz perturbation pressure gradients [Pa/m]; rho perturbation density [kg/m³];
    rho0 reference density [kg/m³]; g [m/s²]; Du_Dt, Dv_Dt, Dw_Dt material accelerations [m/s²]; lap_h_* horizontal
    Laplacians and d2*_dz2 vertical second derivatives of u, v, w [1/(m s)]; nu_H, nu_v eddy viscosities [m²/s].
    Returns
    -------
    dict with keys "x", "y", "z"; each a dict of terms [m/s²] — "acceleration", "coriolis", "pressure", "buoyancy"
    (z only), "friction" — and "residual" (their sum; zero when the equation is satisfied).
    Assumptions: thin layer, Boussinesq, eddy-viscosity friction.
    Validation: V2 zero residuals for a geostrophic + hydrostatic state.  Label: analytic.
    """
    out = {}
    fx = eddy_friction(lap_h_u, d2u_dz2, nu_H, nu_v)
    fy = eddy_friction(lap_h_v, d2v_dz2, nu_H, nu_v)
    fz = eddy_friction(lap_h_w, d2w_dz2, nu_H, nu_v)
    out["x"] = dict(acceleration=_S(Du_Dt), coriolis=_S(-_F(f) * _F(v)), pressure=_S(_F(dpdx) / rho0),
                    friction=_S(-_F(fx)))  # Eq. (13.9a)
    out["y"] = dict(acceleration=_S(Dv_Dt), coriolis=_S(_F(f) * _F(u)), pressure=_S(_F(dpdy) / rho0),
                    friction=_S(-_F(fy)))  # Eq. (13.9b)
    out["z"] = dict(acceleration=_S(Dw_Dt), pressure=_S(_F(dpdz) / rho0), buoyancy=_S(g * _F(rho) / rho0),
                    friction=_S(-_F(fz)))  # Eq. (13.9c)
    for row in out.values():
        row["residual"] = _S(sum(_F(val) for val in row.values()))
    return out



def scale_height(c, g: float = G0):
    """Scale height c²/g of a compressible medium — the depth over which the Boussinesq set stops being accurate.

    Book: §13.3 (text after Eq. (13.2)).  Not the isothermal scale height RT/g of ``core.statics.scale_height``.
    Parameters: c speed of sound [m/s]; g [m/s²].   Returns c²/g [m].
    Assumptions: none (a definition).   Validation: V1 identity.  Label: analytic.
    """
    return _S(_F(c) ** 2 / g)


def buoyancy_frequency_sq(z, rho, rho0: float, g: float = G0):
    """Buoyancy frequency squared from a potential-density profile, N² = −(g/ρ₀) dρ/dz.

    Book: §13.2 (unnumbered definition; ρ is the *potential* density, ρ₀ a constant reference value).
    Parameters
    ----------
    z : heights [m], increasing upward (non-uniform spacing allowed, ≥ 3 nodes).
    rho : potential density at the nodes [kg/m³].   rho0 : reference density [kg/m³].   g [m/s²].
    Returns
    -------
    N2 : N² at the nodes [1/s²]; positive where the profile is statically stable.
    Assumptions: Boussinesq (constant ρ₀ in the denominator).
    Numerics: explicit 3-point stencils, second order at every node including the ends (not ``np.gradient``).
    Validation: V1 exact for a linear profile; V3 order 2 on an exponential profile.  Label: converged.
    """
    return -g / float(rho0) * _d1(_F(z), _F(rho))


# =====================================================================================================================
# 2. Geostrophic balance, thermal wind, Taylor–Proudman (§13.5)
# =====================================================================================================================

def geostrophic_velocity(dpdx, dpdy, f, rho0):
    """Geostrophic velocity from the horizontal pressure gradient: u = −(1/ρ₀f) ∂p/∂y, v = (1/ρ₀f) ∂p/∂x.

    Book: §13.5, Eqs. (13.11)–(13.12).
    Parameters
    ----------
    dpdx, dpdy : horizontal pressure gradients [Pa/m].   f : Coriolis parameter [1/s], either sign, non-zero.
    rho0 : reference density [kg/m³].
    Returns
    -------
    (u, v) [m/s].  The flow is along isobars with low pressure on the left for f > 0 and on the right for f < 0.
    Assumptions: Ro ≪ 1, E ≪ 1, steady.   Raises ValueError at f = 0.
    Validation: V1 u·∇p = 0; V7 the circulation round a low reverses with the sign of f.  Label: analytic.
    """
    fa = _nonzero_f(f, "geostrophic_velocity")
    return _tup(-_F(dpdy) / (rho0 * fa), _F(dpdx) / (rho0 * fa))  # Eqs. (13.11), (13.12)


def geostrophic_from_field(p, dx: float, dy: float, f, rho0: float, periodic: bool = False):
    """Geostrophic velocity of a gridded pressure field.

    Book: §13.5, Eqs. (13.11)–(13.12).
    Parameters
    ----------
    p : pressure at the nodes of a uniform grid, shape [ny, nx] = (y, x) [Pa].   dx, dy : spacings [m].
    f : Coriolis parameter [1/s] — a scalar or an array broadcastable to p (a β-plane column ``f[:, None]``).
    rho0 : reference density [kg/m³].   periodic : wrap the differences (doubly periodic field) if True.
    Returns
    -------
    (u, v) arrays [m/s] at the nodes.
    Numerics: centred differences inside, second-order one-sided stencils at the edges; second order overall.
    Assumptions: as :func:`geostrophic_velocity`.
    Validation: V1 Gaussian high against its analytic wind; V3 order 2; V4 ∇·u = 0 for constant f.  Label: converged.
    """
    p = _F(p)
    return geostrophic_velocity(_ddx_grid(p, dx, -1, periodic), _ddx_grid(p, dy, -2, periodic), f, rho0)


def geostrophic_from_height(eta, dx: float, dy: float, f, g: float = G0, periodic: bool = False):
    """Geostrophic velocity of a gridded surface (or interface) height: u = −(g/f) ∂η/∂y, v = (g/f) ∂η/∂x.

    Book: §13.15, Eq. (13.116) (the shallow-water form of (13.11)–(13.12) with p = ρ g η).
    Parameters: eta [m] on a uniform grid [ny, nx]; dx, dy [m]; f [1/s] scalar or broadcastable; g [m/s²] (pass the
    reduced gravity for an interface); periodic as :func:`geostrophic_from_field`.
    Returns (u, v) [m/s] at the nodes.
    Assumptions: hydrostatic homogeneous layer, Ro ≪ 1.   Validation: V1 against :func:`geostrophic_from_field` with
    p = ρ g η.  Label: converged.
    """
    return geostrophic_from_field(_F(eta) * g, dx, dy, f, 1.0, periodic)  # Eq. (13.116)


def geostrophic_streamfunction(p, f, rho0):
    """Stream function of geostrophic flow on an f-plane, ψ = p/(f ρ₀), with u = −∂ψ/∂y, v = ∂ψ/∂x.

    Book: §13.5 (text after (13.11)–(13.12): "p/fρ₀ can be regarded as a stream function").  Sign convention of this
    chapter (trap T14): ch04 and ch11 use u = +∂ψ/∂y.
    Parameters: p [Pa]; f constant [1/s]; rho0 [kg/m³].   Returns psi [m²/s].
    Assumptions: f constant.   Validation: V1 its curl reproduces :func:`geostrophic_velocity`.  Label: analytic.
    """
    fa = _nonzero_f(f, "geostrophic_streamfunction")
    return _S(_F(p) / (fa * rho0))


def rossby_number(U, f, L):
    """Rossby number Ro = U / (|f| L): nonlinear acceleration over Coriolis force.

    Book: §13.5, Eq. (13.13).  The ch04 sibling ``core.similarity.rossby_number(U, Omega, l)`` is U/(2Ωl); the two
    agree at the pole and differ by sin(latitude) elsewhere (trap T17).
    Parameters: U velocity scale [m/s]; f [1/s] (either sign); L horizontal length scale [m].
    Returns Ro (dimensionless, ≥ 0).   Raises ValueError at f = 0.
    Assumptions: a scale estimate.   Validation: V1 equals the ch04 form at the pole.  Label: analytic.
    """
    fa = _nonzero_f(f, "rossby_number")
    return _S(np.abs(_F(U)) / (np.abs(fa) * _F(L)))  # Eq. (13.13)


def ekman_number(nu, f, L):
    """Ekman number E = ν / (|f| L²): viscous force over Coriolis force.

    Book: §13.5, Eq. (13.18).
    Parameters: nu (eddy or molecular) viscosity [m²/s]; f [1/s]; L length scale [m] (use the depth for a vertical
    Ekman number: E = (δ/L)²/2 with δ the Ekman thickness).
    Returns E (dimensionless).   Raises ValueError at f = 0.
    Assumptions: a scale estimate.   Validation: V1 E = ½(δ/L)² with :func:`ekman_depth`.  Label: analytic.
    """
    fa = _nonzero_f(f, "ekman_number")
    return _S(_F(nu) / (np.abs(fa) * _F(L) ** 2))  # Eq. (13.18)


def thermal_wind_shear(drho_dx, drho_dy, f, rho0, g: float = G0):
    """Thermal wind: vertical shear of the geostrophic velocity from the horizontal density gradient.

    Book: §13.5, Eq. (13.15): ∂v/∂z = −(g/ρ₀f) ∂ρ/∂x and ∂u/∂z = (g/ρ₀f) ∂ρ/∂y.
    Parameters
    ----------
    drho_dx, drho_dy : horizontal density gradients [kg/m⁴].   f [1/s], non-zero.   rho0 [kg/m³].   g [m/s²].
    Returns
    -------
    (du_dz, dv_dz) [1/s] — note the order: the u-shear first.
    Assumptions: geostrophic and hydrostatic balance.   Raises ValueError at f = 0.
    Validation: V2 z-derivative of the geostrophic velocity of a hydrostatic pressure field; V7 denser to the north
    gives eastward shear for f > 0 and westward for f < 0.  Label: symbolic.
    """
    fa = _nonzero_f(f, "thermal_wind_shear")
    return _tup(g / (rho0 * fa) * _F(drho_dy), -g / (rho0 * fa) * _F(drho_dx))  # Eq. (13.15)


def thermal_wind_from_temperature(dTdx, dTdy, f, *, alpha, g: float = G0):
    """Thermal wind in temperature form: ∂u/∂z = −(gα/f) ∂T/∂y, ∂v/∂z = (gα/f) ∂T/∂x.

    Book: §13.5 states this in words only; **ours — not in the book** as a formula.  It is Eq. (13.15) with the linear
    equation of state ρ′/ρ₀ = −α T′ of §13.3.
    Parameters
    ----------
    dTdx, dTdy : horizontal temperature gradients [K/m] (potential temperature for a deep layer of air).
    f [1/s], non-zero.
    alpha : thermal expansion coefficient [1/K] — **required keyword, no default** (1/T₀ for a perfect gas, about
        2e-4 for sea water; passing it silently wrong changes the answer by an order of magnitude).
    g [m/s²].
    Returns
    -------
    (du_dz, dv_dz) [1/s].
    Assumptions: as :func:`thermal_wind_shear`, plus density a function of temperature only.
    Validation: V1 equals :func:`thermal_wind_shear` under ρ′ = −ρ₀αT′.  Label: analytic.
    """
    fa = _nonzero_f(f, "thermal_wind_from_temperature")
    a = float(alpha)
    return _tup(-g * a / fa * _F(dTdy), g * a / fa * _F(dTdx))



def thermal_wind_integrate(z, drho_dy, f, rho0, g: float = G0, u_ref: float = 0.0):
    """Zonal geostrophic wind u(z) obtained by integrating the thermal-wind shear upward from the first level.

    Book: §13.5, Eq. (13.15b), integrated in z (the integration is ours).
    Parameters
    ----------
    z : levels [m], increasing (≥ 2).   drho_dy : northward density gradient at the levels [kg/m⁴] (array or scalar).
    f [1/s], non-zero.   rho0 [kg/m³].   g [m/s²].   u_ref : wind at z[0] [m/s].
    Returns
    -------
    u : wind at the levels [m/s].
    Numerics: cumulative trapezoid rule (second order; exact for a shear linear in z).
    Assumptions: as :func:`thermal_wind_shear`.   Validation: V1 uniform gradient gives a linear profile.
    Label: analytic.
    """
    z = _F(z)
    fa = float(_nonzero_f(f, "thermal_wind_integrate"))
    shear = g / (rho0 * fa) * (_F(drho_dy) * np.ones_like(z))  # Eq. (13.15b)
    return u_ref + np.concatenate([[0.0], np.cumsum(0.5 * (shear[1:] + shear[:-1]) * np.diff(z))])


def taylor_proudman_residual(u_fn: Callable, x, h: float = 1e-4) -> dict:
    """∂u/∂z of a velocity field at a point — zero for a flow that obeys the Taylor–Proudman theorem.

    Book: §13.5, Eq. (13.21) ∂u/∂z = 0 (from (13.19) and (13.20)).
    Parameters
    ----------
    u_fn : callable (x, y, z) → (u, v, w) [m/s].   x : the point (x, y, z) [m].   h : difference step [m].
    Returns
    -------
    dict: ``du_dz``, ``dv_dz``, ``dw_dz`` [1/s] by a central difference (second order in h).
    Assumptions: the theorem itself needs Ro ≪ 1, E ≪ 1, steady flow, uniform density.
    Validation: V1 zero for a z-independent flow, non-zero for a sheared one.  Label: analytic.
    """
    x0, y0, z0 = (float(c) for c in x)
    up = np.asarray(u_fn(x0, y0, z0 + h), dtype=float)
    um = np.asarray(u_fn(x0, y0, z0 - h), dtype=float)
    d = (up - um) / (2.0 * h)
    return dict(du_dz=float(d[0]), dv_dz=float(d[1]), dw_dz=float(d[2]))



def parcel_adjust(t, G, f, r: float = 0.0, V0: complex = 0j):
    """Velocity of a parcel released in a uniform pressure gradient on an f-plane, with optional linear drag.

    Book: the mechanism is described in words in §13.5 ("how is such a motion set up?"); the closed form is **ours —
    not in the book**.  It solves dV/dt + (r + i f) V = −G, the horizontal momentum equations (13.9a, b) for a
    uniform pressure-gradient acceleration, with V = u + i v.
    Parameters
    ----------
    t : time since release [s].
    G : complex pressure-gradient acceleration (1/ρ₀)(∂p/∂x + i ∂p/∂y) [m/s²].
    f : Coriolis parameter [1/s] (either sign; f = 0 allowed if r > 0).
    r : linear drag rate [1/s] (0 = frictionless).   V0 : complex velocity at t = 0 [m/s].
    Returns
    -------
    V : complex velocity u + i v [m/s].  For r = 0 the parcel circles the geostrophic velocity i G / f with period
    2π/|f| (inertial oscillation); for r > 0 it spirals into −G/(r + i f), which has a component down the pressure
    gradient.
    Assumptions: uniform G and f, linear drag, no feedback of the parcel on the pressure field.
    Validation: V2 residual of the ODE; V7 r = 0, t-average = geostrophic velocity.  Label: analytic.
    """
    lam = r + 1j * float(f)
    if lam == 0:
        return _S(V0 - G * _F(t))
    Vs = -G / lam
    return _S(Vs + (V0 - Vs) * np.exp(-lam * _F(t)))


# =====================================================================================================================
# 3. Ekman layers (§13.6, §13.7)
# =====================================================================================================================

def ekman_depth(nu_v, f, convention: str = "efold"):
    """Ekman-layer thickness δ = sqrt(2 ν_v / |f|).

    Book: §13.6, Eq. (13.29) (and §13.7 after (13.40)).
    Parameters
    ----------
    nu_v : vertical eddy viscosity [m²/s].   f : Coriolis parameter [1/s], either sign, non-zero.
    convention : "efold" (default) — the e-folding scale δ of the book; "pi" — the depth πδ at which the current of
        the surface spiral first opposes the surface current (the oceanographer's "Ekman depth"; trap T6).
    Returns
    -------
    thickness [m].
    Assumptions: constant ν_v.   Raises ValueError at f = 0 (no Ekman layer at the equator).
    Validation: V1 4ν_v doubles δ; equals ch08's Stokes-layer thickness with ω → |f|.  Label: analytic.
    """
    fa = np.abs(_nonzero_f(f, "ekman_depth"))
    delta = np.sqrt(2.0 * _F(nu_v) / fa)  # Eq. (13.29)
    if convention == "efold":
        return _S(delta)
    if convention == "pi":
        return _S(np.pi * delta)
    raise ValueError('convention must be "efold" or "pi"')


def _ekman_lambda(nu_v, f):
    """(1 + i s)/δ with s = sign(f): the root of λ² = i f / ν_v that has a positive real part."""
    s = np.sign(float(f))
    return (1.0 + 1j * s) / float(ekman_depth(nu_v, f))


def ekman_surface(z, tau_x, tau_y, rho, nu_v, f, U_g: float = 0.0, V_g: float = 0.0, as_complex: bool = False):
    """Velocity in the Ekman layer below a wind-stressed sea surface (the Ekman spiral).

    Book: §13.6, Eqs. (13.22)–(13.29) and the unnumbered solution after (13.29); the interior flow of Fig. 13.7 is
    added as (U_g, V_g).  Corrected far-field condition: u, v → (U_g, V_g) as z → −∞ (the book prints z → ∞, slip #4).
    Written for either sign of f: with s = sign(f),
    V = V_g + (tau_x + i tau_y)(1 − i s) / (ρ sqrt(2 ν_v |f|)) · exp[(1 + i s) z / δ].
    Parameters
    ----------
    z : height [m], z = 0 at the sea surface, **z ≤ 0 in the water**.
    tau_x, tau_y : wind stress on the surface [N/m²].   rho : water density [kg/m³].
    nu_v : vertical eddy viscosity [m²/s] (constant).   f : Coriolis parameter [1/s], non-zero, either sign.
    U_g, V_g : depth-independent geostrophic velocity under the layer [m/s].   as_complex : return V = u + i v.
    Returns
    -------
    (u, v) [m/s], or the complex V.  NaN for z > 0 (above the water).
    Assumptions: steady, horizontally uniform, constant ν_v, infinitely deep.
    Validation: V2 residual of (13.22)–(13.23) and the stress condition for f > 0 and f < 0; V1 surface current
    |tau|/(ρ sqrt(ν_v |f|)) at 45° to the right (f > 0) or left (f < 0) of the stress.  Label: symbolic.
    """
    z = _F(z)
    lam = _ekman_lambda(nu_v, f)
    A = (tau_x + 1j * tau_y) / (rho * nu_v * lam)   # from rho nu_v dV/dz = tau at z = 0, Eqs. (13.24)-(13.25)
    with np.errstate(over="ignore", invalid="ignore"):
        V = np.where(z <= 0.0, (U_g + 1j * V_g) + A * np.exp(lam * np.minimum(z, 0.0)), np.nan + 0j)  # Eq. (13.28), B = 0
    if as_complex:
        return _S(V)
    return _tup(V.real, V.imag)


def ekman_transport(tau_x, tau_y, rho, f):
    """Ekman volume transport of the surface layer: (M_x, M_y) = (tau_y, −tau_x) / (ρ f).

    Book: §13.6, Eq. (13.30) (stress along x: ∫u dz = 0, ∫v dz = −tau/(ρ f)).
    Parameters: tau_x, tau_y wind stress [N/m²]; rho [kg/m³]; f [1/s], non-zero.
    Returns (M_x, M_y) [m²/s] — the depth-integrated velocity (relative to any interior geostrophic flow), 90° to the
    right of the stress for f > 0 and 90° to the left for f < 0, **independent of the eddy viscosity**.
    Assumptions: steady, horizontally uniform; stress vanishes at depth.  Follows from integrating −ρ f v = dτ/dz, so it
    holds for any K(z).   Raises ValueError at f = 0.
    Validation: V4 quadrature of :func:`ekman_surface` and of :func:`ekman_solve` with three K(z).  Label: analytic.
    """
    fa = _nonzero_f(f, "ekman_transport")
    return _tup(_F(tau_y) / (rho * fa), -_F(tau_x) / (rho * fa))  # Eq. (13.30)


def ekman_transport_partial(z, tau_x, tau_y, rho, nu_v, f):
    """Transport of the surface Ekman layer between depth z and the surface, in closed form.

    Book: §13.6 — the integral of the unnumbered solution after Eq. (13.29) from z to 0 (the partial integral is ours);
    it tends to Eq. (13.30) as z → −∞.
    Parameters: z ≤ 0 [m]; the rest as :func:`ekman_surface` (no interior flow).
    Returns (M_x, M_y) [m²/s] = ∫_z^0 (u, v) dz′.
    Assumptions: as :func:`ekman_surface`.
    Validation: V1 quadrature; V7 → :func:`ekman_transport` at depth.  Label: analytic.
    """
    z = np.minimum(_F(z), 0.0)
    lam = _ekman_lambda(nu_v, f)
    A = (tau_x + 1j * tau_y) / (rho * nu_v * lam)
    M = A / lam * (1.0 - np.exp(lam * z))
    return _tup(M.real, M.imag)


def ekman_residual(z, u, v, nu_v, f, U_g: float = 0.0, V_g: float = 0.0):
    """Residuals of the Ekman balance for sampled profiles u(z), v(z).

    Book: §13.6, Eqs. (13.22)–(13.23), and §13.7, Eqs. (13.33)–(13.34) (the pressure gradient written as f U_g).
    r_x = −f (v − V_g) − ν_v u″,   r_y = f (u − U_g) − ν_v v″.
    Parameters
    ----------
    z : nodes [m], monotonic (≥ 3).   u, v : velocities at the nodes [m/s].   nu_v [m²/s].   f [1/s].
    U_g, V_g : geostrophic velocity [m/s].
    Returns
    -------
    (r_x, r_y) arrays [m/s²] at the **interior** nodes z[1:-1] (length len(z) − 2).
    Numerics: 3-point second differences, second order on uniform and smoothly stretched grids.
    Assumptions: steady, horizontally uniform, constant ν_v.
    Validation: V1 on the closed-form spiral (4001 uniform nodes over 14 δ) the residual is 1.4e-6 of f|V| —
    truncation of the stencil, not a property of the solution.  Label: analytic.
    """
    z, u, v = _F(z), _F(u), _F(v)
    return (-f * (v - V_g) - nu_v * _d2(z, u))[1:-1], (f * (u - U_g) - nu_v * _d2(z, v))[1:-1]



def ekman_vorticity_balance(z, tau_x, rho, nu_v, f) -> dict:
    """Why the Ekman layer does not grow: tilting of planetary vorticity against vertical diffusion of horizontal
    vorticity, evaluated on the surface spiral.

    Book: §13.6, Eq. (13.31): −f dv/dz = ν_v d²ω_y/dz² and −f du/dz = ν_v d²ω_x/dz², with ω_x = −dv/dz, ω_y = du/dz.
    Parameters: z ≤ 0 [m]; tau_x stress along x [N/m²]; rho [kg/m³]; nu_v [m²/s]; f [1/s] non-zero.
    Returns dict: ``tilt_x`` (= −f dv/dz) and ``diff_x`` (= ν_v d²ω_y/dz²) — the two sides of the first member;
    ``tilt_y`` (= −f du/dz) and ``diff_y`` (= ν_v d²ω_x/dz²) — the two sides of the second [1/s²]; ``omega_x``,
    ``omega_y`` [1/s].  tilt = diff to round-off.
    Assumptions: as :func:`ekman_surface`.   Derivatives are analytic (V′ = λV).
    Validation: V1 the two sides agree at round-off for both signs of f.  Label: analytic.
    """
    lam = _ekman_lambda(nu_v, f)
    V = _F(0.0) + ekman_surface(z, tau_x, 0.0, rho, nu_v, f, as_complex=True)
    V1, V3 = lam * V, lam ** 3 * V
    return dict(tilt_x=_S(-f * V1.imag), diff_x=_S(nu_v * V3.real),       # Eq. (13.31a): omega_y = du/dz
                tilt_y=_S(-f * V1.real), diff_y=_S(nu_v * (-V3.imag)),    # Eq. (13.31b): omega_x = -dv/dz
                omega_x=_S(-V1.imag), omega_y=_S(V1.real))



def ekman_solve(z, K, f, *, tau=None, rho=None, V_g: complex = 0.0, bottom: str = "decay", top: str = "stress"):
    """Ekman layer with a depth-dependent eddy viscosity K(z), solved numerically.

    Book: **ours — not in the book** (the book notes that the constant-viscosity spiral is rarely observed).  The
    equation is the complex form (13.27)/(13.37) with K inside the derivative: d/dz (K dV/dz) = i f (V − V_g).
    Parameters
    ----------
    z : nodes [m], strictly increasing (non-uniform spacing allowed; cluster nodes where K is small).
    K : eddy viscosity **at the faces** between nodes [m²/s]: a callable K(z) (evaluated at the face mid-points), a
        scalar, or an array of length ``len(z) − 1``.  An array of length ``len(z)`` (node values) is refused.
    f : Coriolis parameter [1/s], non-zero, either sign.
    tau : surface stress tau_x + i tau_y (complex) or a pair (tau_x, tau_y) [N/m²]; needed for ``top="stress"``.
    rho : density [kg/m³]; needed for ``top="stress"``.
    V_g : complex geostrophic velocity U_g + i V_g [m/s].
    bottom : "decay" (V = V_g at z[0]: the deep end of a surface layer), "noslip" (V = 0 at z[0]: a solid surface) or
        "stressfree" (K dV/dz = 0 at z[0]).
    top : "stress" (ρ K dV/dz = tau at z[-1]), "geostrophic" (V = V_g at z[-1]: the top of a bottom layer) or
        "stressfree".
    Returns
    -------
    V : complex velocity u + i v at the nodes [m/s].
    Numerics (ours): finite-volume discretisation on node-centred control volumes with face fluxes K (V_{j+1} −
    V_j)/h (second order on uniform and smoothly stretched grids); half control volumes carry the flux boundary
    conditions; one complex tridiagonal solve (``scipy.linalg.solve_banded``).
    Assumptions: steady, horizontally uniform, K > 0.
    Validation: V1 constant K against :func:`ekman_surface` / :func:`ekman_bottom`; V3 order 2; V4 transport =
    tau/(ρ f) for any K(z).  Label: converged.
    """
    z = _F(z)
    n = z.size
    if n < 3 or np.any(np.diff(z) <= 0):
        raise ValueError("ekman_solve: z must be strictly increasing with at least 3 nodes")
    fa = float(_nonzero_f(f, "ekman_solve"))
    h = np.diff(z)
    zf = 0.5 * (z[1:] + z[:-1])
    if callable(K):
        Kf = _F(K(zf)) * np.ones(n - 1)
    else:
        Kf = _F(K)
        if Kf.ndim == 0:
            Kf = np.full(n - 1, float(Kf))
        elif Kf.size == n:
            raise ValueError("ekman_solve: K has one value per node; give K at the faces (length len(z) - 1) or a "
                             "callable K(z)")
        elif Kf.size != n - 1:
            raise ValueError("ekman_solve: K must be a scalar, a callable or an array of length len(z) - 1")
    if np.any(Kf <= 0):
        raise ValueError("ekman_solve: K must be positive")
    Vg = complex(V_g)
    c = Kf / h                                   # face conductances
    lo = np.zeros(n, dtype=complex)
    di = np.zeros(n, dtype=complex)
    up = np.zeros(n, dtype=complex)
    rhs = np.zeros(n, dtype=complex)
    vol = 0.5 * (h[1:] + h[:-1])
    lo[1:-1], up[1:-1] = c[:-1], c[1:]
    di[1:-1] = -(c[:-1] + c[1:]) - 1j * fa * vol   # d/dz(K dV/dz) - i f V = -i f V_g
    rhs[1:-1] = -1j * fa * vol * Vg
    if bottom == "decay":
        di[0], rhs[0] = 1.0, Vg
    elif bottom == "noslip":
        di[0], rhs[0] = 1.0, 0.0
    elif bottom == "stressfree":
        di[0], up[0], rhs[0] = -c[0] - 1j * fa * 0.5 * h[0], c[0], -1j * fa * 0.5 * h[0] * Vg
    else:
        raise ValueError('bottom must be "decay", "noslip" or "stressfree"')
    if top in ("stress", "stressfree"):
        if top == "stress":
            if tau is None or rho is None:
                raise ValueError('ekman_solve: top="stress" needs tau and rho')
            t = complex(tau[0], tau[1]) if isinstance(tau, (tuple, list, np.ndarray)) else complex(tau)
            flux = t / float(rho)
        else:
            flux = 0.0
        di[-1], lo[-1] = -c[-1] - 1j * fa * 0.5 * h[-1], c[-1]
        rhs[-1] = -flux - 1j * fa * 0.5 * h[-1] * Vg
    elif top == "geostrophic":
        di[-1], rhs[-1] = 1.0, Vg
    else:
        raise ValueError('top must be "stress", "geostrophic" or "stressfree"')
    ab = np.zeros((3, n), dtype=complex)
    ab[0, 1:], ab[1], ab[2, :-1] = up[:-1], di, lo[1:]
    return solve_banded((1, 1), ab, rhs)


def ekman_pumping(tau_x, tau_y, dx: float, dy: float, rho: float, f, periodic: bool = False):
    """Ekman pumping velocity at the base of the surface layer from a gridded wind stress: w_E = (1/ρ) curl_z(τ/f).

    Book: **ours — not in the book** (the book gives the transport, Eq. (13.30), and the consequence for lows and
    highs in §13.7).  It is the horizontal divergence of the Ekman transport (tau_y, −tau_x)/(ρ f), which continuity
    turns into a vertical velocity at the base of the layer (w = 0 at the surface).
    Parameters
    ----------
    tau_x, tau_y : wind stress at the nodes of a uniform grid [ny, nx] [N/m²].   dx, dy [m].   rho [kg/m³].
    f : Coriolis parameter [1/s], scalar or broadcastable (``f[:, None]`` for a β-plane), non-zero.
    periodic : wrap the differences if True.
    Returns
    -------
    w_E [m/s] at the nodes; positive = upwelling (suction of interior water into the layer).
    Numerics: centred differences, second order (one-sided second order at the edges).
    Assumptions: steady Ekman layer much thinner than the scale of the wind.
    Validation: V1 manufactured stress with a known curl; V3 order 2.  Label: converged.
    """
    fa = _nonzero_f(f, "ekman_pumping")
    return (_ddx_grid(_F(tau_y) / fa * np.ones_like(_F(tau_y)), dx, -1, periodic)
            - _ddx_grid(_F(tau_x) / fa * np.ones_like(_F(tau_x)), dy, -2, periodic)) / rho


def ekman_pumping_from_curl(curl_tau, rho, f, beta: float = 0.0, tau_x: float = 0.0):
    """Ekman pumping from the local curl of the wind stress: w_E = curl_z τ/(ρ f) + β tau_x/(ρ f²).

    Book: **ours — not in the book**; the point form of :func:`ekman_pumping` (the second term comes from the
    northward change of f inside curl_z(τ/f)).
    Parameters: curl_tau = ∂tau_y/∂x − ∂tau_x/∂y [N/m³]; rho [kg/m³]; f [1/s] non-zero; beta [1/(m s)]; tau_x [N/m²].
    Returns w_E [m/s], positive upward.   Raises ValueError at f = 0.
    Assumptions: as :func:`ekman_pumping`.   Validation: V1 against the gridded form.  Label: analytic.
    """
    fa = _nonzero_f(f, "ekman_pumping_from_curl")
    return _S(_F(curl_tau) / (rho * fa) + beta * _F(tau_x) / (rho * fa ** 2))


def sverdrup_transport(curl_tau, rho, beta):
    """Sverdrup balance: depth-integrated northward transport V = curl_z τ / (ρ β).

    Book: **ours — not in the book** (§13.7 names the wind-driven circulation and Stommel's result without presenting
    them).  It follows from the Ekman pumping above and the linear vorticity balance β v = f ∂w/∂z of the interior.
    Parameters: curl_tau [N/m³]; rho [kg/m³]; beta [1/(m s)] > 0.
    Returns V [m²/s] — northward volume transport per unit width, Ekman layer included.
    Assumptions: steady, linear, frictionless interior, no flow at great depth.
    Validation: V1 identity.  Label: analytic.
    """
    if np.any(_F(beta) <= 0):
        raise ValueError("sverdrup_transport: beta must be positive")
    return _S(_F(curl_tau) / (rho * _F(beta)))


def ekman_bottom(z, U_g, V_g, nu_v, f, as_complex: bool = False):
    """Velocity in the Ekman layer above a rigid surface under a geostrophic flow (U_g, V_g).

    Book: §13.7, Eqs. (13.33)–(13.41): u = U[1 − e^{−z/δ} cos(z/δ)], v = U e^{−z/δ} sin(z/δ) for f > 0 and flow
    along x.  Written for either sign of f and any direction of the interior flow: with s = sign(f),
    V = (U_g + i V_g) [1 − exp(−(1 + i s) z/δ)].
    Parameters
    ----------
    z : height above the surface [m], z ≥ 0.   U_g, V_g : geostrophic velocity far above [m/s].
    nu_v : vertical eddy viscosity [m²/s].   f [1/s], non-zero.   as_complex : return V = u + i v.
    Returns
    -------
    (u, v) [m/s] or complex V.  NaN for z < 0.  Near the surface the flow is turned 45° to the left of the
    geostrophic flow for f > 0 (toward low pressure) and 45° to the right for f < 0 (also toward low pressure).
    Assumptions: steady, horizontally uniform, constant ν_v, no slip at z = 0.
    Validation: V2 residual of (13.33)–(13.34) and both boundary conditions; V1 v first vanishes at z = πδ, where
    u = U (1 + e^{−π}) (4.3 % above U); the largest overshoot of u is e^{−3π/4}/sqrt(2) = 6.7 % at z = 3πδ/4.
    Label: symbolic.
    """
    z = _F(z)
    lam = _ekman_lambda(nu_v, f)
    with np.errstate(over="ignore", invalid="ignore"):
        V = np.where(z >= 0.0, (U_g + 1j * V_g) * (1.0 - np.exp(-lam * np.maximum(z, 0.0))), np.nan + 0j)  # Eq. (13.41)
    if as_complex:
        return _S(V)
    return _tup(V.real, V.imag)


def ekman_bottom_transport(U_g, V_g, nu_v, f):
    """Transport of the bottom Ekman layer *relative to the geostrophic flow*: ∫₀^∞ (V − V_gcomplex) dz.

    Book: §13.7, the unnumbered integral after Eq. (13.41): ∫₀^∞ v dz = U sqrt(ν_v/2f) = ½ U δ (flow along x, f > 0).
    In general, with s = sign(f): M_x = −(δ/2)(U_g + s V_g), M_y = (δ/2)(s U_g − V_g).
    Parameters: U_g, V_g geostrophic velocity [m/s]; nu_v [m²/s]; f [1/s] non-zero.
    Returns (M_x, M_y) [m²/s].  For flow along +x: a deficit ½Uδ along the flow and a cross-flow transport ½Uδ toward
    low pressure — to the left for f > 0, to the right for f < 0.
    Assumptions: as :func:`ekman_bottom`.   Validation: V1 quadrature of :func:`ekman_bottom`.  Label: analytic.
    """
    lam = _ekman_lambda(nu_v, f)
    M = -(_F(U_g) + 1j * _F(V_g)) / lam
    return _tup(M.real, M.imag)


def ekman_pumping_bottom(zeta_g, nu_v, f):
    """Vertical velocity at the top of a bottom Ekman layer under geostrophic vorticity ζ_g: w = sign(f) (δ/2) ζ_g.

    Book: **ours — not in the book** as a formula; §13.7 gives the consequence (air converges and rises in a low,
    sinks in a high).  It is minus the divergence of :func:`ekman_bottom_transport` for a non-divergent interior flow.
    Parameters: zeta_g = ∂V_g/∂x − ∂U_g/∂y [1/s]; nu_v [m²/s]; f [1/s] non-zero.
    Returns w [m/s], positive upward: rising over a cyclone (ζ_g of the sign of f) in either hemisphere.
    Assumptions: as :func:`ekman_bottom`, layer thin against the scale of the interior flow, f constant.
    Validation: V1 divergence of the gridded transport of a manufactured vortex.  Label: analytic.
    """
    return _S(np.sign(float(f)) * 0.5 * float(ekman_depth(nu_v, f)) * _F(zeta_g))


def ekman_finite_depth(z, tau_x, tau_y, U_g, V_g, H, nu_v, rho, f):
    """Ekman flow between a wind-stressed surface (z = 0) and a no-slip bottom (z = −H): a surface layer over a
    bottom layer in water of finite depth.

    Book: **ours — not in the book**: the exact solution of d²V/dz² = (i f/ν_v)(V − V_g) (Eqs. (13.27), (13.37)) with
    ρ ν_v dV/dz = tau at z = 0 and V = 0 at z = −H; the book only describes an observed profile of this kind.
    Parameters
    ----------
    z : height [m], −H ≤ z ≤ 0.   tau_x, tau_y : surface stress [N/m²].   U_g, V_g : interior geostrophic velocity
    [m/s].   H : water depth [m].   nu_v [m²/s].   rho [kg/m³].   f [1/s] non-zero.
    Returns
    -------
    V : complex velocity u + i v [m/s] (an array for array z); NaN outside −H ≤ z ≤ 0.
    Numerics: written with decaying exponentials only, so H ≫ δ does not overflow.
    Assumptions: steady, horizontally uniform, constant ν_v.
    Validation: V2 residual and both boundary conditions; V7 H ≫ δ reproduces :func:`ekman_surface` near the surface
    and :func:`ekman_bottom` near the bottom.  Label: analytic.
    """
    z = _F(z)
    lam = _ekman_lambda(nu_v, f)
    Vg = U_g + 1j * V_g
    zz = np.clip(z, -H, 0.0)
    den = 1.0 + np.exp(-2.0 * lam * H)
    cosh_ratio = np.exp(-lam * (zz + H)) * (1.0 + np.exp(2.0 * lam * zz)) / den      # cosh(lam z)/cosh(lam H)
    sinh_ratio = np.exp(lam * zz) * (1.0 - np.exp(-2.0 * lam * (zz + H))) / den      # sinh(lam (z+H))/cosh(lam H)
    V = Vg - Vg * cosh_ratio + (tau_x + 1j * tau_y) / (rho * nu_v * lam) * sinh_ratio
    return _S(np.where((z >= -H) & (z <= 0.0), V, np.nan + 0j))



def ekman_force_balance(z, U_g, nu_v, f, rho: float = 1.0) -> dict:
    """The three forces on a parcel in the bottom Ekman layer and the angle its velocity makes with the isobars.

    Book: §13.7, Eqs. (13.32)–(13.34) and the force triangle of the section's last figure (interior flow along +x).
    Coriolis force = ρ f (v, −u); pressure force = ρ f (0, U_g); friction = ρ ν_v (u″, v″) = ρ f (−v, u − U_g).
    Parameters
    ----------
    z : height above the surface [m], ≥ 0.   U_g : geostrophic speed along +x [m/s].   nu_v [m²/s].
    f [1/s] non-zero.   rho : density [kg/m³] (default 1: forces per unit mass).
    Returns
    -------
    dict: ``coriolis``, ``pressure``, ``friction`` — each a pair (F_x, F_y) [N/m³; m/s² for rho = 1]; ``u``, ``v``,
    ``speed`` [m/s]; ``angle_to_isobars`` [rad, ≥ 0] — the angle between the velocity and the geostrophic direction,
    a magnitude: toward low pressure for 0 < z < πδ (to the left of U_g for f > 0, to the right for f < 0; π/4 in
    the limit z → 0, which is the value returned at z = 0), and slightly toward HIGH pressure for πδ < z < 2πδ, where
    the cross-isobar velocity reverses (under 1°; the sign of ``v`` tells which); ``sum`` — the vector sum of the three forces (zero to round-off).
    Assumptions: as :func:`ekman_bottom`.
    Validation: V1 the three vectors sum to zero at every height; angle π/4 at the ground, → 0 above the layer.
    Label: analytic.
    """
    z = _F(z)
    u, v = (_F(c) for c in ekman_bottom(z, U_g, 0.0, nu_v, f))
    cor = (rho * f * v, -rho * f * u)
    pre = (0.0 * u, rho * f * U_g + 0.0 * u)
    fri = (-rho * f * v, rho * f * (u - U_g))
    if float(U_g) == 0.0:
        ang = np.zeros(z.shape)
    else:   # angle of V relative to the geostrophic direction; the limit pi/4 at the ground, where V = 0
        ang = np.where(z > 0.0, np.abs(np.angle(np.where(z > 0.0, (u + 1j * v) / U_g, 1.0))), np.pi / 4.0)
    return dict(coriolis=_tup(*cor), pressure=_tup(*pre), friction=_tup(*fri), u=_S(u), v=_S(v),
                speed=_S(np.hypot(u, v)), angle_to_isobars=_S(ang),
                sum=_tup(cor[0] + pre[0] + fri[0], cor[1] + pre[1] + fri[1]))



def eddy_viscosity_from_depth(delta, f):
    """Eddy viscosity implied by an observed Ekman-layer thickness: ν_v = |f| δ² / 2.

    Book: §13.7 (the inversion of Eq. (13.29) used to estimate ν_v from the depth of the atmospheric boundary layer).
    Parameters: delta e-folding thickness [m]; f [1/s] non-zero.   Returns nu_v [m²/s].
    Assumptions: constant-viscosity Ekman layer.   Validation: V1 inverse of :func:`ekman_depth`.  Label: analytic.
    """
    fa = np.abs(_nonzero_f(f, "eddy_viscosity_from_depth"))
    return _S(0.5 * fa * _F(delta) ** 2)


# =====================================================================================================================
# 4. Shallow-water waves with rotation (§13.8, §13.10–§13.12)
# =====================================================================================================================

def long_wave_speed(H, g: float = G0):
    """Speed of long gravity waves on a layer of depth H: c = sqrt(g H).

    Book: §13.9, Eq. (13.70) and §13.12, Eq. (13.86) (the shallow-water limit of ch07).
    Parameters: H depth (or equivalent depth) [m] ≥ 0; g gravity or reduced gravity [m/s²].   Returns c [m/s].
    Assumptions: wavelength ≫ depth, small amplitude.
    Validation: V7 shallow limit of ``core.waves.phase_speed``.  Label: analytic.
    """
    return _S(np.sqrt(g * _F(H)))


def equivalent_depth(c, g: float = G0):
    """Equivalent depth H_e = c²/g of a vertical mode with long-wave speed c.

    Book: §13.8, Eq. (13.46) and §13.9, Eq. (13.62).
    Parameters: c modal speed [m/s]; g [m/s²].   Returns H_e [m].
    Assumptions: hydrostatic normal mode (each mode then obeys the shallow-water equations with H → H_e).
    Validation: V1 inverse of :func:`long_wave_speed`.  Label: analytic.
    """
    return _S(_F(c) ** 2 / g)


def baroclinic_mode_speed(N, H, n: int = 1):
    """Long-wave speed of the n-th baroclinic mode of a layer of uniform buoyancy frequency: c_n = N H / (n π).

    Book: §13.9, Eq. (13.71) (n = 1, 2, 3, …; the rigid-lid limit of Eq. (13.69)).
    Parameters: N buoyancy frequency [rad/s]; H depth [m]; n mode number ≥ 1.   Returns c_n [m/s].
    Raises ValueError for n < 1 (n = 0 is the barotropic mode, c₀ = sqrt(gH), :func:`long_wave_speed`).
    Assumptions: uniform N, hydrostatic, rigid lid (relative error against the free surface about N²H/(g n²π²)).
    Validation: V1 against ``core.vertical_modes.modes_uniform_N``.  Label: analytic.
    """
    if int(n) < 1:
        raise ValueError("baroclinic_mode_speed: n must be >= 1 (n = 0 is the barotropic mode, sqrt(g H))")
    return _S(_F(N) * _F(H) / (int(n) * np.pi))  # Eq. (13.71)


def shallow_water_discriminant(k, l, c, f0, beta):
    """Discriminant 4 (c²K² + f₀²)³ − 27 (c² β k)² of the shallow-water cubic; ≥ 0 means three real frequencies.

    Book: §13.10, Eq. (13.76) (the book states that all roots are real; the criterion is ours).
    Parameters: k, l wavenumbers [rad/m]; c long-wave speed [m/s]; f0 [1/s]; beta [1/(m s)].
    Returns the discriminant [1/s⁶].  It is negative only for planetary-scale waves next to the equator, where the
    β-plane cubic is not valid anyway.
    Assumptions: β-plane, f² → f₀².   Validation: V1 sign agrees with ``numpy.roots``.  Label: analytic.
    """
    P = _F(c) ** 2 * (_F(k) ** 2 + _F(l) ** 2) + _F(f0) ** 2
    return _S(4.0 * P ** 3 - 27.0 * (_F(c) ** 2 * beta * _F(k)) ** 2)


def shallow_water_omega(k, l, c, f0, beta):
    """The three frequencies of linear shallow-water waves on a β-plane: roots of ω³ − (c²K² + f₀²) ω − c² β k = 0.

    Book: §13.10, Eq. (13.76) (repeated as Eq. (13.121)), from the v-equation (13.75) with v ∝ exp[i(kx + ly − ωt)].
    Parameters
    ----------
    k, l : eastward and northward wavenumbers [rad/m] (signed).   c : long-wave speed sqrt(gH) or a modal c_n [m/s].
    f0 : Coriolis parameter [1/s].   beta : df/dy [1/(m s)].
    Returns
    -------
    (omega_minus, omega_rossby, omega_plus) [rad/s]: the westward- and eastward-travelling inertia–gravity (Poincaré)
    roots (negative and positive) and the slow root between them, whose sign is opposite to that of β k (westward
    phase).  The labels "westward" for the negative and "eastward" for the positive fast root hold for k > 0 only:
    the phase moves along x at ω/k, so for k < 0 the negative root travels east.  The tuple is always in ascending
    order.  All three are NaN where the discriminant is negative (no three real roots: see
    :func:`shallow_water_discriminant`).
    Numerics: the two fast roots by the trigonometric solution of the depressed cubic; the slow root — smaller by
    one to five orders of magnitude — from the product of the roots (ω₋ ω_R ω₊ = c²βk), so that it does not lose its
    digits to cancellation; then four Newton steps on each root.  Measured on our ocean inputs (external and first
    baroclinic mode, wavelengths 730 and 3100 km, 35°): cubic residual of every root below 1e-15 of its largest
    term, sum of the roots below 3e-16 of ω₊, agreement with ``numpy.roots`` (companion matrix) to 5e-16 for all
    three roots.  NaN (all three) when the discriminant is negative — for the
    external mode that is a band of planetary-scale waves at low latitude, where β c²|k| is no longer small against
    (f₀² + c²K²)^{3/2} and the cubic has one real and two complex roots.
    Assumptions: linear, hydrostatic, β-plane with f² → f₀².
    Validation: V1 cubic residual; independent route ``numpy.roots``; V7 β → 0: ±sqrt(f₀² + c²K²) and 0.
    Label: analytic.
    """
    k, l, c, f0 = np.broadcast_arrays(_F(k), _F(l), _F(c), _F(f0))
    P = c ** 2 * (k ** 2 + l ** 2) + f0 ** 2          # omega^3 - P omega - q = 0
    q = c ** 2 * beta * k
    with np.errstate(invalid="ignore", divide="ignore"):
        r = 2.0 * np.sqrt(P / 3.0)
        arg = np.where(P > 0, 3.0 * q / np.where(P > 0, P * r, 1.0), 0.0)
        ok = np.abs(arg) <= 1.0
        th = np.arccos(np.clip(arg, -1.0, 1.0)) / 3.0
        wp = r * np.cos(th)                             # largest root
        wm = r * np.cos(th + 2.0 * np.pi / 3.0)         # most negative root
        wr = np.where(wp * wm != 0, q / np.where(wp * wm != 0, wp * wm, 1.0), 0.0)   # product of the roots = +q
        for _ in range(4):   # Newton polish of all three (quadratic convergence; the start is already close)
            wm, wr, wp = (w - np.where(3.0 * w ** 2 - P != 0, (w ** 3 - P * w - q) / np.where(3.0 * w ** 2 - P != 0,
                          3.0 * w ** 2 - P, 1.0), 0.0) for w in (wm, wr, wp))
    bad = ~ok
    wm, wr, wp = (np.where(bad, np.nan, w) for w in (wm, wr, wp))
    return _tup(wm, wr, wp)  # Eq. (13.76)


def kelvin_omega(k, c):
    """Frequency of a Kelvin wave, ω = c k (non-dispersive; any frequency is possible).

    Book: §13.12, the line before Eq. (13.86) (ω = ±k sqrt(gH)); the sign that is trapped against the coast is chosen
    by :func:`kelvin_decay_side`.
    Parameters: k alongshore wavenumber [rad/m]; c = sqrt(gH) or a modal speed [m/s].   Returns omega [rad/s].
    Assumptions: as :func:`kelvin_wave`.   Validation: V1 identity.  Label: analytic.
    """
    return _S(_F(c) * _F(k))


def shallow_water_branches(k, l, c, f0, beta) -> dict:
    """The named wave branches at one wavenumber: the three roots of the cubic and the Kelvin line.

    Book: §13.10, Eq. (13.76); §13.11, Eq. (13.82); §13.12, Eq. (13.86); §13.15, Eq. (13.118).
    Parameters: as :func:`shallow_water_omega`.
    Returns dict [rad/s]: ``poincare_plus``, ``poincare_minus``, ``rossby`` (the cubic's roots) and ``kelvin`` (c·k,
    the boundary-trapped companion, which is not a root of the cubic).
    Assumptions: as :func:`shallow_water_omega`.   Validation: V7 β → 0: ±sqrt(f₀² + c²K²) and 0.  Label: analytic.
    """
    wm, wr, wp = shallow_water_omega(k, l, c, f0, beta)
    return dict(poincare_plus=wp, poincare_minus=wm, rossby=wr, kelvin=kelvin_omega(k, c))



def shallow_water_regime(omega, f0) -> dict:
    """Name the regime of a shallow-water wave from its frequency relative to f, and the term of the cubic that may
    be dropped there.

    Book: §13.10, the three ranges discussed after Eq. (13.76).  The closing sentence of the low-frequency paragraph
    prints "ω ≫ f" where ω ≪ f is meant (slip #7); the band edges (|ω/f₀| = 3 and 1) are ours.
    Parameters: omega [rad/s]; f0 [1/s] (non-zero).
    Returns dict: ``ratio`` = |ω/f₀|; ``regime`` — "high" (ratio > 3: gravity waves that do not feel rotation),
    "near-inertial" (1 ≤ ratio ≤ 3: Poincaré waves), "low" (ratio < 1: Rossby waves); ``neglect`` — the term of
    ω³ − c²K²ω − f₀²ω − c²βk = 0 that may be dropped: "f0^2 omega" (high; the β term is negligible too),
    "c^2 beta k" (near-inertial), "omega^3" (low).
    Assumptions: a labelling aid, not a theorem.   Validation: V1 the three bands.  Label: analytic.
    """
    r = abs(float(omega)) / abs(float(_nonzero_f(f0, "shallow_water_regime")))
    if r > 3.0:
        return dict(ratio=r, regime="high", neglect="f0^2 omega")
    if r >= 1.0:
        return dict(ratio=r, regime="near-inertial", neglect="c^2 beta k")
    return dict(ratio=r, regime="low", neglect="omega^3")



def dispersion_term_sizes(k, l, c, f0, beta, omega) -> dict:
    """The four terms of the shallow-water cubic at a given frequency, to see which ones balance.

    Book: §13.10, Eq. (13.76): ω³ − c²K² ω − f₀² ω − c² β k = 0.
    Parameters: k, l [rad/m]; c [m/s]; f0 [1/s]; beta [1/(m s)]; omega [rad/s] (a root, or any trial frequency).
    Returns dict [1/s³]: ``omega3`` (ω³), ``gravity`` (−c²K²ω), ``rotation`` (−f₀²ω), ``beta`` (−c²βk), ``sum`` (zero
    for a root).
    Use: on a fast root at ω = 3f the ω³ term is the largest of the four; on the slow root it is the smallest — the
    corrected reading of the book's closing sentence (slip #7).
    Assumptions: as :func:`shallow_water_omega`.   Validation: V1 sum = 0 at each root.  Label: analytic.
    """
    K2 = float(k) ** 2 + float(l) ** 2
    w = float(omega)
    out = dict(omega3=w ** 3, gravity=-c ** 2 * K2 * w, rotation=-f0 ** 2 * w, beta=-c ** 2 * beta * float(k))
    out["sum"] = out["omega3"] + out["gravity"] + out["rotation"] + out["beta"]
    return out



def poincare_omega(K, f, c):
    """Frequency of rotational gravity (Poincaré, Sverdrup) waves: ω = sqrt(f² + c² K²).

    **Call by keyword** — ``poincare_omega(K=…, f=…, c=…)``: the order (K, f, c) differs from the (k, l, c, f0, …)
    order of ``shallow_water_omega`` and from (k, l, f, c) of ``poincare_group_velocity``.

    Book: §13.11, Eqs. (13.81)–(13.82) with c² = gH.
    Parameters: K horizontal wavenumber magnitude sqrt(k² + l²) [rad/m]; f [1/s] (either sign); c [m/s].
    Returns omega [rad/s] ≥ |f| (positive root).
    Assumptions: linear shallow water on an f-plane, unbounded.
    Validation: V7 f = 0 gives cK (ch07 long waves); K = 0 gives |f| (inertial oscillation).  Label: analytic.
    """
    return _S(np.sqrt(_F(f) ** 2 + (_F(c) * _F(K)) ** 2))  # Eq. (13.82)


def poincare_group_velocity(k, l, f, c):
    """Group velocity of Poincaré waves: c_g = c² (k, l) / ω.

    Book: **ours — not in the book** (the gradient of Eq. (13.82) in wavenumber space).
    Parameters: k, l [rad/m]; f [1/s]; c [m/s].   Returns (c_gx, c_gy) [m/s]; magnitude below c, and
    c_phase · c_group = c².  Zero at K = 0.
    Assumptions: as :func:`poincare_omega`.   Validation: V1 finite-difference ∂ω/∂k.  Label: analytic.
    """
    k, l = _F(k), _F(l)
    w = _F(poincare_omega(np.sqrt(k ** 2 + l ** 2), f, c))
    with np.errstate(invalid="ignore", divide="ignore"):
        return _tup(np.where(w > 0, c ** 2 * k / np.where(w > 0, w, 1.0), 0.0),
                    np.where(w > 0, c ** 2 * l / np.where(w > 0, w, 1.0), 0.0))


def poincare_amplitudes(k, l, omega, f, g, eta_hat):
    """Complex velocity amplitudes of a rotating shallow-water plane wave of surface amplitude η̂.

    Book: §13.11, Eq. (13.80): û = g η̂ (ωk + i f l)/(ω² − f²), v̂ = g η̂ (−i f k + ωl)/(ω² − f²).
    Parameters: k, l [rad/m]; omega [rad/s]; f [1/s]; g [m/s²]; eta_hat [m] (may be complex).
    Returns (u_hat, v_hat) complex [m/s] for fields ∝ exp[i(kx + ly − ωt)].  NaN at ω² = f² (the inertial limit, where
    η̂ → 0 and the ratio is undefined).
    Assumptions: as :func:`poincare_omega`.   Validation: V2 satisfies (13.77)–(13.79).  Label: symbolic.
    """
    d = _F(omega) ** 2 - _F(f) ** 2
    with np.errstate(invalid="ignore", divide="ignore"):
        fac = np.where(d != 0, g * np.asarray(eta_hat) / np.where(d != 0, d, 1.0), np.nan)
    return _S(fac * (omega * _F(k) + 1j * f * _F(l))), _S(fac * (-1j * f * _F(k) + omega * _F(l)))  # Eq. (13.80)


def poincare_fields(x, y, t, k, l, eta_hat, H, f, g: float = G0):
    """Surface height and velocity of a Poincaré plane wave (real parts), with ω = +sqrt(f² + gH K²).

    Book: §13.11, Eqs. (13.80), (13.82), (13.83).
    Parameters: x, y [m]; t [s]; k, l [rad/m] (not both zero); eta_hat real amplitude [m]; H [m]; f [1/s]; g [m/s²].
    Returns (eta, u, v): η [m] and the velocity components [m/s].
    Raises ValueError for k = l = 0 (the inertial oscillation, where η = 0: use :func:`inertial_oscillation`).
    Assumptions: as :func:`poincare_omega`.   Validation: V1 zero residual in (13.45).  Label: analytic.
    """
    K = np.hypot(float(k), float(l))
    if K == 0:
        raise ValueError("poincare_fields: k = l = 0 is the inertial oscillation (eta = 0); use inertial_oscillation")
    w = float(poincare_omega(K, f, np.sqrt(g * H)))
    uh, vh = poincare_amplitudes(k, l, w, f, g, eta_hat)
    ph = np.exp(1j * (k * _F(x) + l * _F(y) - w * _F(t)))
    return _tup((eta_hat * ph).real, (uh * ph).real, (vh * ph).real)



def poincare_orbit(t, k, eta_hat, H, f, g: float = G0) -> dict:
    """Velocity hodograph and particle orbit at a fixed point of a Poincaré wave travelling along +x.

    Book: §13.11, Eq. (13.83) at x = 0: u = (ω η̂/kH) cos ωt, v = −(f η̂/kH) sin ωt.
    Parameters: t [s]; k > 0 [rad/m]; eta_hat [m]; H [m]; f [1/s]; g [m/s²].
    Returns dict: ``u``, ``v`` [m/s] (the velocity at x = 0); ``x``, ``y`` [m] (the particle's displacement about its
    mean position, the time integral of the velocity: x = (η̂/kH) sin ωt, y = (f η̂/ωkH) cos ωt); ``omega`` [rad/s];
    ``axis_ratio`` = ω/|f| (major/minor, the same for hodograph and orbit; ``inf`` for f = 0: a straight line);
    ``sense`` ("clockwise" for f > 0, "counter-clockwise" for f < 0, "rectilinear" for f = 0).
    The long axis lies along the direction of propagation.  The book's figure of this "orbit" shows velocity
    amplitudes: the path has the same shape with both axes divided by ω (trap T7).
    Assumptions: as :func:`poincare_omega`; small displacements (the velocity is evaluated at the mean position).
    Validation: V1 signed area of the hodograph; dx/dt = u.  Label: analytic.
    """
    w = float(poincare_omega(abs(float(k)), f, np.sqrt(g * H)))
    a = eta_hat / (k * H)
    tt = _F(t)
    sense = "clockwise" if f > 0 else ("counter-clockwise" if f < 0 else "rectilinear")
    return dict(u=_S(w * a * np.cos(w * tt)), v=_S(-f * a * np.sin(w * tt)),  # Eq. (13.83)
                x=_S(a * np.sin(w * tt)), y=_S(f * a / w * np.cos(w * tt)), omega=w,
                axis_ratio=w / abs(f) if f != 0 else float("inf"), sense=sense)



def inertial_oscillation(t, u0, v0, f):
    """Inertial motion: the velocity and the path of a parcel subject to the Coriolis force alone.

    Book: §13.11 "Inertial Motion" (unnumbered): ∂u/∂t − f v = 0, ∂v/∂t + f u = 0, with solution u = q cos ft,
    v = −q sin ft for a parcel starting along x.  General start: u + i v = (u₀ + i v₀) exp(−i f t), and the position
    relative to the starting point is x + i y = i (u + i v − u₀ − i v₀)/f.
    Parameters: t [s]; u0, v0 initial velocity [m/s]; f [1/s] (f = 0 gives the straight line x = u₀t, y = v₀t).
    Returns (u, v, x, y): velocity [m/s] and displacement [m].  Constant speed; a circle of radius q/|f| traced
    clockwise for f > 0 and counter-clockwise for f < 0, period 2π/|f|.
    Assumptions: no pressure gradient, no friction, f constant.
    Validation: V4 speed conserved; returns to the start after one inertial period.  Label: analytic.
    """
    tt = _F(t)
    V0 = u0 + 1j * v0
    V = V0 * np.exp(-1j * f * tt)
    X = 1j * (V - V0) / f if f != 0 else V0 * tt
    return _tup(V.real, V.imag, X.real, X.imag)



def inertial_radius(q, f):
    """Radius of an inertial circle, r = q/|f|.

    Book: §13.11 "Inertial Motion" (unnumbered; Coriolis acceleration f q = centripetal acceleration r f²).
    Parameters: q speed [m/s]; f [1/s] non-zero.   Returns r [m].
    Assumptions: as :func:`inertial_oscillation`.   Validation: V1 path of the integrated velocity.  Label: analytic.
    """
    return _S(np.abs(_F(q)) / np.abs(_nonzero_f(f, "inertial_radius")))


def rossby_radius(c, f):
    """Rossby radius of deformation Λ = c/|f|.

    Book: §13.12 (unnumbered definition after Eq. (13.87)).
    Parameters: c long-wave speed [m/s] — sqrt(gH) (external), sqrt(g′H) or a modal c_n (internal); f [1/s] non-zero.
    Returns Lambda [m].
    Assumptions: none (a definition).   Validation: V1 e-folding distance of :func:`kelvin_wave`.  Label: analytic.
    """
    return _S(_F(c) / np.abs(_nonzero_f(f, "rossby_radius")))


def rossby_radius_internal(N, H, f, n: int = 1, with_pi: bool = True):
    """Internal Rossby radius of a uniformly stratified layer: N H/(n π |f|), or N H/|f| without the π.

    Book: §13.12 (Λ = NH/πf for n = 1) and §13.17 (Λ ≡ HN/f, "it is usual to omit the factor π") — two definitions
    under one name (trap T11).
    Parameters: N [rad/s]; H [m]; f [1/s] non-zero; n mode number ≥ 1; with_pi True for c_n/|f|, False for the §13.17
    convention (n is then ignored).
    Returns Lambda [m].
    Assumptions: uniform N, rigid lid.   Validation: V1 :func:`baroclinic_mode_speed` / |f|.  Label: analytic.
    """
    fa = np.abs(_nonzero_f(f, "rossby_radius_internal"))
    if with_pi:
        return _S(_F(baroclinic_mode_speed(N, H, n)) / fa)
    return _S(_F(N) * _F(H) / fa)


def rossby_radius_two_layer(H1, rho1, rho2, f, g: float = G0):
    """Internal Rossby radius of a thin upper layer over a deep lower one: sqrt(g′ H₁)/|f|, g′ = g (ρ₂ − ρ₁)/ρ₂.

    Book: §13.12 (c = sqrt(g′H) with the reduced gravity defined on the following page; ch07 Eq. (7.117)).
    Parameters: H1 upper-layer thickness [m]; rho1 < rho2 layer densities [kg/m³]; f [1/s] non-zero; g [m/s²].
    Returns Lambda [m].
    Assumptions: lower layer much deeper than the upper; long waves.
    Validation: V1 ``core.waves.reduced_gravity_book``.  Label: analytic.
    """
    if np.any(_F(rho2) <= _F(rho1)):
        raise ValueError("rossby_radius_two_layer: needs rho2 > rho1 (stable layering)")
    gp = g * (_F(rho2) - _F(rho1)) / _F(rho2)
    return rossby_radius(np.sqrt(gp * _F(H1)), f)


def kelvin_decay_side(f, direction) -> dict:
    """Is a Kelvin wave travelling in a given direction along the coast y = 0 trapped in the sea y ≥ 0?

    Book: §13.12 (the choice of the decaying solution η̂ = η₀ e^{−fy/c} before Eq. (13.87), made there for f > 0 and
    travel along +x).  With travel direction d = ±1 along x the cross-shore structure is exp(−d f y/c).
    Parameters: f [1/s] non-zero; direction +1 (toward +x) or −1.
    Returns dict: ``trapped`` (bool: d·sign(f) > 0 — the amplitude decays away from the coast into y > 0) and
    ``coast_on`` ("right" or "left": the side of the direction of travel on which a trapped wave keeps the coast in
    this hemisphere — right for f > 0, left for f < 0).
    Assumptions: as :func:`kelvin_wave`.   Validation: V7 both hemispheres and both directions.  Label: analytic.
    """
    fa = float(_nonzero_f(f, "kelvin_decay_side"))
    d = 1 if float(direction) > 0 else -1
    return dict(trapped=bool(d * fa > 0), coast_on="right" if fa > 0 else "left")



def kelvin_wave(x, y, t, eta0, k, H, f, g: float = G0, direction=None):
    """Coastal Kelvin wave: surface height and alongshore current in the sea y ≥ 0 bounded by the coast y = 0.

    Book: §13.12, Eqs. (13.84)–(13.87): η = η₀ e^{−fy/c} cos k(x − ct), u = η₀ sqrt(g/H) e^{−fy/c} cos k(x − ct),
    v = 0, c = sqrt(gH) (written there for f > 0 and travel toward +x).
    General form used here, with d = direction and Λ = c/|f|:
    η = η₀ e^{−y/Λ} cos k(x − d c t),   u = d sqrt(g/H) η.
    Parameters
    ----------
    x : alongshore coordinate [m].   y : distance from the coast [m], ≥ 0.   t [s].
    eta0 : amplitude at the coast [m].   k : alongshore wavenumber [rad/m].   H : depth (or equivalent depth) [m].
    f [1/s] non-zero.   g [m/s²] (reduced gravity for an internal wave).
    direction : +1 (toward +x) or −1; default sign(f) — the only direction that is trapped.
    Returns
    -------
    (eta, u) : surface height [m] and alongshore velocity [m/s]; v = 0 identically.  NaN for y < 0 (land).
    Raises ValueError if direction·sign(f) < 0: that wave would grow offshore and is not a solution.
    Assumptions: linear shallow water, f-plane, straight vertical coast, v ≡ 0.
    Validation: V2 satisfies the three equations (13.84); V1 e-folding distance c/|f|.  Label: symbolic.
    """
    fa = float(_nonzero_f(f, "kelvin_wave"))
    d = (1 if fa > 0 else -1) if direction is None else (1 if float(direction) > 0 else -1)
    if d * fa < 0:
        raise ValueError("kelvin_wave: a wave travelling this way would grow away from the coast (direction*sign(f) < 0)"
                         " — in this hemisphere the trapped wave travels the other way")
    c = float(np.sqrt(g * H))
    x, y, t = np.broadcast_arrays(_F(x), _F(y), _F(t))
    env = np.where(y >= 0.0, np.exp(-np.abs(fa) * np.maximum(y, 0.0) / c), np.nan)
    eta = eta0 * env * np.cos(k * (x - d * c * t))  # Eq. (13.87)
    return _tup(eta, d * np.sqrt(g / H) * eta)



def kelvin_residuals(x, y, t, eta0, k, H, f, g: float = G0, direction=None, h: float = 1.0, ht: float = 1.0) -> dict:
    """Residuals of the three Kelvin-wave equations evaluated on :func:`kelvin_wave` by central differences.

    Book: §13.12, Eq. (13.84): ∂η/∂t + H ∂u/∂x = 0, ∂u/∂t = −g ∂η/∂x, f u = −g ∂η/∂y.
    Parameters: as :func:`kelvin_wave` (a scalar point with y > h); h difference step in x and y [m]; ht in t [s].
    Returns dict (dimensionless): ``continuity``, ``x_momentum``, ``y_geostrophy`` — each residual divided by the
    size of its largest term at this distance from the coast (|η₀| k c e^{−y/Λ}, g|η₀| k e^{−y/Λ},
    |f| sqrt(g/H)|η₀| e^{−y/Λ}).  With the default steps they are of order (k h)² or round-off, far below 1e-6.
    Assumptions: as :func:`kelvin_wave`.   Validation: V2 this is the check.  Label: analytic.
    """
    kw = dict(eta0=eta0, k=k, H=H, f=f, g=g, direction=direction)
    c = float(np.sqrt(g * H))
    fa = abs(float(f))
    amp = abs(eta0) * float(np.exp(-fa * float(y) / c))
    e_tp, u_tp = kelvin_wave(x, y, t + ht, **kw)
    e_tm, u_tm = kelvin_wave(x, y, t - ht, **kw)
    e_xp, u_xp = kelvin_wave(x + h, y, t, **kw)
    e_xm, u_xm = kelvin_wave(x - h, y, t, **kw)
    e_yp, _ = kelvin_wave(x, y + h, t, **kw)
    e_ym, _ = kelvin_wave(x, y - h, t, **kw)
    _, u0 = kelvin_wave(x, y, t, **kw)
    return dict(continuity=float(((e_tp - e_tm) / (2 * ht) + H * (u_xp - u_xm) / (2 * h)) / (amp * abs(k) * c)),
                x_momentum=float(((u_tp - u_tm) / (2 * ht) + g * (e_xp - e_xm) / (2 * h)) / (g * amp * abs(k))),
                y_geostrophy=float((f * u0 + g * (e_yp - e_ym) / (2 * h)) / (fa * np.sqrt(g / H) * amp)))



def geostrophic_adjustment_1d(x, eta0, H, f, g: float = G0):
    """End state of the geostrophic adjustment of a step in surface height.

    Book: **ours — not in the book** (the book describes how a balanced flow is set up, §13.5, and defines the Rossby
    radius, §13.12).  Initial state: fluid at rest with η = η₀ sgn(x).  The linear shallow-water equations (13.45)
    conserve ζ − f η/H at every x; requiring the end state to be geostrophic (f v = g ∂η/∂x) gives
    η″ − η/Λ² = −η₀ sgn(x)/Λ², hence, with Λ = sqrt(gH)/|f|, c = sqrt(gH) and s = sign(f),
    η = η₀ sgn(x) (1 − e^{−|x|/Λ}),   v = s (g η₀/c) e^{−|x|/Λ},   u = 0.
    Parameters: x [m]; eta0 half-height of the step [m]; H [m]; f [1/s] non-zero; g [m/s²].
    Returns (eta, v): the adjusted surface [m] and the jet along the step [m/s] (high water on its right for f > 0,
    on its left for f < 0).
    Raises ValueError at f = 0: without rotation there is no steady end state (everything radiates away).
    Assumptions: linear, inviscid, f-plane, unbounded in x, uniform in y; the gravity waves have left.
    Validation: V1 satisfies geostrophy and linear potential-vorticity conservation; V3 the 1-D shallow-water march
    converges to it.  Label: analytic (ours).
    """
    fa = float(_nonzero_f(f, "geostrophic_adjustment_1d"))
    x = _F(x)
    c = float(np.sqrt(g * H))
    e = np.exp(-np.abs(x) * abs(fa) / c)
    return _tup(eta0 * np.sign(x) * (1.0 - e), np.sign(fa) * g * eta0 / c * e)



def adjustment_energy(eta0, H, f, g: float = G0, L=None, rho: float = 1000.0) -> dict:
    """Energy budget of the adjusted step of :func:`geostrophic_adjustment_1d`, per unit length of the step.

    Book: **ours — not in the book**.
    ΔPE = ½ ρ g ∫ (η_initial² − η²) dx,   KE = ½ ρ H ∫ v² dx,   integrated over all x (or over |x| ≤ L).
    Parameters: eta0 [m]; H [m]; f [1/s] non-zero; g [m/s²]; L half-width of the integration [m] (None: infinite);
    rho density [kg/m³].
    Returns dict [J/m]: ``pe_released`` (= (3/2) ρ g η₀² Λ for L = None), ``ke_jet`` (= (1/2) ρ g η₀² Λ),
    ``radiated`` (their difference, carried away by the gravity waves) and ``ratio`` = ke_jet/pe_released — computed
    from the two integrals (1/3 for L = None: only a third of the potential energy released stays in the jet).
    Assumptions: as :func:`geostrophic_adjustment_1d`.
    Validation: V1 quadrature of the closed-form profiles.  Label: analytic (ours).
    """
    fa = float(_nonzero_f(f, "adjustment_energy"))
    Lam = float(np.sqrt(g * H)) / abs(fa)
    e1 = 1.0 if L is None or np.isinf(L) else -np.expm1(-L / Lam)          # 1 - exp(-L/Lambda)
    e2 = 1.0 if L is None or np.isinf(L) else -np.expm1(-2.0 * L / Lam)    # 1 - exp(-2L/Lambda)
    pe = rho * g * eta0 ** 2 * (2.0 * Lam * e1 - 0.5 * Lam * e2)
    ke = rho * g * eta0 ** 2 * 0.5 * Lam * e2
    return dict(pe_released=float(pe), ke_jet=float(ke), radiated=float(pe - ke),
                ratio=float(ke / pe) if pe > 0 else float("nan"))



# =====================================================================================================================
# 5. Internal waves with rotation (§13.14)
# =====================================================================================================================

def inertia_gravity_m2(k, l, omega, N, f):
    """Vertical wavenumber squared of an internal wave in a rotating stratified fluid.

    Book: §13.14, Eq. (13.99) (and (13.109) for l = 0): m² = (k² + l²)(N² − ω²)/(ω² − f²).
    Parameters: k, l [rad/m]; omega [rad/s]; N buoyancy frequency (local value; may be an array) [rad/s]; f [1/s].
    Returns m² [1/m²]: positive for |f| < ω < N (the wave propagates vertically), negative outside the band (trapped).
    Raises ValueError at ω² = f² (the inertial limit: the vertical wavenumber is unbounded).
    Assumptions: linear, Boussinesq, f-plane, N varying slowly (WKB) or uniform.
    Validation: V1 equivalent to :func:`inertia_gravity_omega`.  Label: analytic.
    """
    d = _F(omega) ** 2 - _F(f) ** 2
    if np.any(d == 0):
        raise ValueError("inertia_gravity_m2: omega^2 = f^2 (the inertial limit) — m^2 is unbounded")
    return _S((_F(k) ** 2 + _F(l) ** 2) * (_F(N) ** 2 - _F(omega) ** 2) / d)  # Eq. (13.99)



def inertia_gravity_band(omega, N, f) -> dict:
    """Is a frequency inside the internal-wave band |f| < ω < N?

    Book: §13.14 (unnumbered inequality f < ω < N after Eq. (13.100), assuming N > f).
    Parameters: omega [rad/s]; N [rad/s]; f [1/s] (N > |f| assumed).
    Returns dict: ``propagating`` (bool: |f| < ω < N) and ``where`` ("below f", "in band" or "above N").
    Assumptions: as :func:`inertia_gravity_m2`.   Validation: V1 sign of m².  Label: analytic.
    """
    w, fa, Na = abs(float(omega)), abs(float(f)), abs(float(N))
    where = "below f" if w <= fa else ("above N" if w >= Na else "in band")
    return dict(propagating=where == "in band", where=where)



def inertia_gravity_omega(k, m, N, f, l=0.0, approx: str | None = None):
    """Frequency of inertia–gravity waves: ω² = (N² K_h² + f² m²)/(K_h² + m²) = f² sin²θ + N² cos²θ.

    Book: §13.14, Eq. (13.112) (ω² − f² = (k²/m²)(N² − ω²)) and the unnumbered form on the following page, θ the
    angle of the wavevector above the horizontal, tan θ = m/k.
    Parameters
    ----------
    k, l : horizontal wavenumbers [rad/m] (K_h² = k² + l²).   m : vertical wavenumber [rad/m].   N [rad/s].   f [1/s].
    approx : None (exact) or one of the three limiting forms listed after the dispersion relation —
        "nonrotating" (ω ~ N: N K_h/sqrt(K_h² + m²)), "hydrostatic" (ω ~ f: sqrt(f² + N²K_h²/m²)),
        "midrange" (f ≪ ω ≪ N: N K_h/|m|).
    Returns
    -------
    omega [rad/s] ≥ 0.  Depends on the direction of the wavevector only.  Complex-safe in k and m (for complex-step
    group velocities).  With f = 0 and approx None it delegates to ``core.waves.internal_wave_omega`` (ch07).
    Assumptions: linear, Boussinesq, f-plane, uniform N.
    Validation: V7 f = 0 equals ch07; ω → |f| as k/m → 0, ω → N as m/k → 0.  Label: analytic.
    """
    cx = np.iscomplexobj(k) or np.iscomplexobj(m) or np.iscomplexobj(l)
    kk, mm, ll = (np.asarray(a) if cx else _F(a) for a in (k, m, l))
    kh2 = kk * kk + ll * ll
    if approx is None:
        if float(f) == 0.0:
            return internal_wave_omega(k, m, float(N), l)
        return _S(np.sqrt((N ** 2 * kh2 + f ** 2 * mm * mm) / (kh2 + mm * mm)))  # Eq. (13.112)
    if approx == "nonrotating":
        return _S(np.sqrt(N ** 2 * kh2 / (kh2 + mm * mm)))
    if approx == "hydrostatic":
        return _S(np.sqrt(f ** 2 + kh2 * N ** 2 / (mm * mm)))
    if approx == "midrange":
        return _S(np.sqrt(kh2 * N ** 2 / (mm * mm)))
    raise ValueError('approx must be None, "nonrotating", "hydrostatic" or "midrange"')



def inertia_gravity_regime(omega, N, f) -> dict:
    """Which limiting form of the internal-wave dispersion relation fits a frequency, and how well each one does.

    Book: §13.14, the three limiting cases listed after Eq. (13.112), applied to m² of Eq. (13.109).
    Parameters: omega, N [rad/s]; f [1/s] (|f| < ω < N).
    Returns dict: ``err_nonrotating`` = −f²/ω² (f² dropped against ω²), ``err_hydrostatic`` = ω²/(N² − ω²) (ω²
    dropped against N²), ``err_midrange`` (both dropped) — the relative error in m² of each approximation at this ω;
    ``regime`` — "mid" when both single approximations are within 10 % (|f| ≪ ω ≪ N), otherwise "near-inertial"
    (the hydrostatic form is the better one) or "near-buoyancy" (the non-rotating form is the better one).  The
    10 % threshold is ours.
    Raises ValueError outside the band (there is no propagating wave to classify).
    Assumptions: as :func:`inertia_gravity_m2`.   Validation: V1 errors against the exact m².  Label: analytic.
    """
    w2, f2, N2 = float(omega) ** 2, float(f) ** 2, float(N) ** 2
    if not f2 < w2 < N2:
        raise ValueError("inertia_gravity_regime: omega is outside the band |f| < omega < N")
    e_nr = -f2 / w2
    e_hy = w2 / (N2 - w2)
    e_mid = N2 * (w2 - f2) / (w2 * (N2 - w2)) - 1.0
    if abs(e_nr) < 0.1 and abs(e_hy) < 0.1:
        regime = "mid"
    else:
        regime = "near-inertial" if abs(e_hy) < abs(e_nr) else "near-buoyancy"
    return dict(regime=regime, err_nonrotating=e_nr, err_hydrostatic=e_hy, err_midrange=e_mid)



def inertia_gravity_group_velocity(k, m, N, f):
    """Group velocity of inertia–gravity waves in the (x, z) plane.

    Book: §13.14 asks for it in an exercise and uses the result (c_g perpendicular to c, opposite vertical
    components); the formula is **ours — not in the book**:
    c_gx = k m² (N² − f²)/(ω K⁴),   c_gz = −m k² (N² − f²)/(ω K⁴),   K² = k² + m².
    Parameters: k, m [rad/m]; N [rad/s]; f [1/s].
    Returns (c_gx, c_gz) [m/s].  k c_gx + m c_gz = 0 (energy moves along the crests); phase up means energy down.
    Zero when the wavevector is exactly horizontal or vertical.
    Assumptions: as :func:`inertia_gravity_omega`, N > |f|.
    Validation: V1 complex-step gradient of :func:`inertia_gravity_omega`.  Label: analytic.
    """
    k, m = _F(k), _F(m)
    K4 = (k ** 2 + m ** 2) ** 2
    w = np.sqrt((N ** 2 * k ** 2 + f ** 2 * m ** 2) / (k ** 2 + m ** 2))
    with np.errstate(invalid="ignore", divide="ignore"):
        fac = np.where(w > 0, (N ** 2 - f ** 2) / np.where(w > 0, w * K4, 1.0), 0.0)
    return _tup(k * m ** 2 * fac, -m * k ** 2 * fac)


def wkb_vertical_structure(z, m, A0: float = 1.0, sign: int = +1):
    """WKB vertical structure of the vertical velocity: ŵ = (A₀/sqrt(m)) exp(± i ∫ m dz).

    Book: §13.14, Eq. (13.104) (from (13.101)–(13.103)).
    Parameters
    ----------
    z : levels [m], increasing.   m : local vertical wavenumber at the levels [rad/m], > 0 (from
    :func:`inertia_gravity_m2`).   A0 : amplitude constant.   sign : +1 upper sign (upward phase propagation for
    k, m, ω > 0), −1 lower sign.
    Returns
    -------
    w_hat : complex array; the phase integral starts at z[0] (cumulative trapezoid rule).  NaN where m is not
    positive (a turning level: the approximation fails there).
    Assumptions: N(z) varies slowly over a vertical wavelength (H m ≫ 1).
    Validation: V1 exact for uniform m; V3 error ∝ 1/(Hm) (measured about 0.15/Hm, ``ch13.wkb_error``).  Label: approximate (WKB).
    """
    z, m = _F(z), _F(m) * np.ones_like(_F(z))
    ok = m > 0
    ms = np.where(ok, m, 1.0)
    phase = np.concatenate([[0.0], np.cumsum(0.5 * (ms[1:] + ms[:-1]) * np.diff(z))])
    return np.where(ok, A0 / np.sqrt(ms) * np.exp(sign * 1j * phase), np.nan + 0j)  # Eq. (13.104)


def inertia_gravity_fields(x, z, t, k, omega, N_fn: Callable, f, A0: float = 1.0, sign: int = +1):
    """Velocity field of an internal wave in slowly varying stratification (WKB), l = 0.

    Book: §13.14, Eq. (13.108) with the dispersion relation (13.109):
    u = ∓(A₀ sqrt(m)/k) cos φ,  v = ∓(A₀ f sqrt(m)/ωk) sin φ,  w = (A₀/sqrt(m)) cos φ,  φ = kx ± ∫ m dz − ωt.
    Parameters
    ----------
    x : horizontal positions [m] (scalar or 1-D).   z : levels [m], increasing, 1-D (the phase integral starts at
    z[0]).   t [s].   k > 0 [rad/m].   omega > 0 [rad/s].   N_fn : callable z → N(z) [rad/s].   f [1/s].
    A0 : amplitude constant [m^{1/2}/s].   sign : +1 (upper signs: phase moves upward) or −1.
    Returns
    -------
    (u, v, w) : real arrays [m/s] on the grid (z, x) — shape (len(z), len(x)), or (len(z),) for scalar x.  NaN at
    levels where m² ≤ 0 (beyond a turning level the wave does not propagate and the WKB form fails).
    Assumptions: N(z) varies slowly over a vertical wavelength (H m ≫ 1).
    Validation: V2 continuity ik û + dŵ/dz = 0 to WKB order; hodograph clockwise for f > 0.
    Label: approximate (WKB).
    """
    z = _F(z)
    xx = _F(x)
    m2 = _F(inertia_gravity_m2(k, 0.0, omega, _F(N_fn(z)) * np.ones_like(z), f))
    ok = m2 > 0
    m = np.sqrt(np.where(ok, m2, 1.0))
    integ = np.concatenate([[0.0], np.cumsum(0.5 * (m[1:] + m[:-1]) * np.diff(z))])
    if xx.ndim:
        m, integ, ok = m[:, None], integ[:, None], ok[:, None] & np.ones(xx.shape, dtype=bool)[None, :]
        xx = xx[None, :]
    ph = k * xx + sign * integ - omega * float(t)
    u = -sign * A0 * np.sqrt(m) / k * np.cos(ph)                       # Eq. (13.108)
    v = -sign * A0 * f * np.sqrt(m) / (omega * k) * np.sin(ph)
    w = A0 / np.sqrt(m) * np.cos(ph)
    return np.where(ok, u, np.nan), np.where(ok, v, np.nan), np.where(ok, w, np.nan)



def inertia_gravity_hodograph(t, omega, f, sign: int = +1):
    """Horizontal velocity hodograph of an internal wave at a fixed point: u = ∓cos ωt, v = ±(f/ω) sin ωt.

    Book: §13.14, Eq. (13.110) (amplitude of u set to one).
    Parameters: t [s]; omega > 0 [rad/s]; f [1/s]; sign +1 (upper signs) or −1.
    Returns (u, v) (dimensionless, scaled by the amplitude of u): an ellipse with axis ratio |f|/ω (minor/major),
    long axis along the direction of propagation, clockwise for f > 0 with either sign choice, anticlockwise for
    f < 0.  (§13.11 quotes the same ellipse as ω/f, major/minor — trap T7.)
    Assumptions: as :func:`inertia_gravity_fields`.   Validation: V1 signed area.  Label: analytic.
    """
    return _tup(-sign * np.cos(omega * _F(t)), sign * f / omega * np.sin(omega * _F(t)))  # Eq. (13.110)


def lee_wave_m(U, N, k):
    """Vertical wavenumber of a stationary lee wave: m = sqrt(N²/U² − k²).

    Book: §13.14 "Lee Wave": Eq. (13.113) with the intrinsic frequency ω = kU of a wave that stands still over the
    ground, i.e. U = N/sqrt(k² + m²) (unnumbered).
    Parameters: U mean wind [m/s], non-zero; N [rad/s]; k horizontal wavenumber of the topography [rad/m] (magnitude).
    Returns m [rad/m] ≥ 0.
    Raises ValueError("evanescent: k > N/U") when the disturbance decays with height instead of propagating.
    Assumptions: uniform U and N, rotation negligible (ω ≫ f), linear.
    Validation: V1 the Doppler-shifted frequency vanishes.  Label: analytic.
    """
    r = (_F(N) / _F(U)) ** 2 - _F(k) ** 2
    if np.any(r < 0):
        raise ValueError("evanescent: k > N/U")
    return _S(np.sqrt(r))



# =====================================================================================================================
# 6. Rossby waves (§13.15)
# =====================================================================================================================

def _rossby_den(k, l, f0, c):
    F = 0.0 if np.isinf(c) else (float(f0) / float(c)) ** 2
    return k * k + l * l + F


def rossby_omega(k, l, beta, f0: float = 0.0, c: float = np.inf, U: float = 0.0):
    """Rossby-wave frequency ω = U k − β k / (k² + l² + f₀²/c²), signed.

    Book: §13.15, Eq. (13.118), with the Doppler shift of Eq. (13.120).  The book "regards ω as positive", which
    forces k < 0; here ω carries its sign: k > 0 gives ω < 0 for U = 0 — either way the phase moves westward (trap T8).
    Parameters
    ----------
    k, l : eastward and northward wavenumbers [rad/m], signed.   beta [1/(m s)].
    f0 [1/s] and c [m/s] : Coriolis parameter and long-wave speed; f₀/c = 1/Λ.  Defaults f0 = 0, c = ∞ give the
    barotropic (non-divergent) wave ω = −βk/K².   U : uniform eastward mean flow [m/s].
    Returns
    -------
    omega [rad/s].  NaN where the denominator vanishes (k = l = 0 with f₀/c = 0).  Complex-safe in k and l.
    Assumptions: quasi-geostrophic, linear, β-plane, ω ≪ f₀; not valid within a few degrees of the equator.
    Validation: V2 plane wave has zero residual in (13.117); V7 the slow root of the cubic (13.76).  Label: symbolic.
    """
    cx = np.iscomplexobj(k) or np.iscomplexobj(l)
    kk, ll = (np.asarray(a) if cx else _F(a) for a in (k, l))
    D = _rossby_den(kk, ll, f0, c)
    with np.errstate(invalid="ignore", divide="ignore"):
        w = U * kk - beta * kk / np.where(D != 0, D, np.nan)  # Eqs. (13.118), (13.120)
    return _S(w)


def rossby_group_velocity(k, l, beta, f0: float = 0.0, c: float = np.inf, U: float = 0.0):
    """Group velocity of Rossby waves: c_gx = U + β (k² − l² − f₀²/c²)/D², c_gy = 2 β k l/D², D = k² + l² + f₀²/c².

    Book: §13.15 defines c_g = (∂ω/∂k, ∂ω/∂l) and shows its direction in a figure; the components are **ours — not
    in the book**.
    Parameters: as :func:`rossby_omega`.
    Returns (c_gx, c_gy) [m/s].  For l = 0 and U = 0: westward for long waves (k² < f₀²/c²), eastward for short ones,
    zero at k = ±f₀/c.
    Assumptions: as :func:`rossby_omega`.   Validation: V1 complex-step gradient of ω.  Label: analytic.
    """
    k, l = _F(k), _F(l)
    D = _rossby_den(k, l, f0, c)
    F = D - k * k - l * l
    with np.errstate(invalid="ignore", divide="ignore"):
        D2 = np.where(D != 0, D * D, np.nan)
        return _tup(U + beta * (k * k - l * l - F) / D2, 2.0 * beta * k * l / D2)


def rossby_phase_speed(k, l, beta, f0: float = 0.0, c: float = np.inf, U: float = 0.0):
    """Eastward phase speed of Rossby waves, c_x = ω/k = U − β/(k² + l² + f₀²/c²).

    Book: §13.15, Eqs. (13.119) and (13.120).
    Parameters: as :func:`rossby_omega`.   Returns c_x [m/s]: negative (westward) for every wave when U = 0.
    Assumptions: as :func:`rossby_omega`.   Validation: V1 ω/k.  Label: analytic.
    """
    D = _rossby_den(_F(k), _F(l), f0, c)
    with np.errstate(invalid="ignore", divide="ignore"):
        return _S(U - beta / np.where(D != 0, D, np.nan))  # Eq. (13.120)


def rossby_max_frequency(beta, f0, c) -> dict:
    """Largest Rossby-wave frequency, ω_max = β c/(2|f₀|), reached at l = 0, k = −|f₀|/c.

    Book: §13.15 (text after the definition of group velocity: ω_max f₀/βc = ½ at kc/f₀ = −1).
    Parameters: beta [1/(m s)]; f0 [1/s] non-zero; c [m/s].
    Returns dict: ``omega_max`` [rad/s] and ``k`` = −|f₀|/c [rad/m] (the wavenumber at which it occurs, with the
    book's convention ω > 0, k < 0; the group velocity vanishes there).
    Assumptions: as :func:`rossby_omega`, U = 0.   Validation: V1 maximum of ω over k.  Label: analytic.
    """
    fa = abs(float(_nonzero_f(f0, "rossby_max_frequency")))
    return dict(omega_max=float(beta) * float(c) / (2.0 * fa), k=-fa / float(c))



def rossby_omega_circle(omega, beta, f0, c) -> dict:
    """The circle of wavenumbers that share one Rossby frequency: (k + β/2ω)² + l² = (β/2ω)² − f₀²/c².

    Book: §13.15 (unnumbered rearrangement of Eq. (13.118), ω > 0 as in the book).
    Parameters: omega [rad/s], non-zero; beta [1/(m s)]; f0 [1/s]; c [m/s].
    Returns dict: ``center_k`` = −β/(2ω), ``center_l`` = 0.0, ``radius`` [rad/m] (the square root of the right-hand
    side; 0.0 when that is negative) and ``exists`` (bool: the right-hand side is positive, i.e. ω does not exceed
    :func:`rossby_max_frequency`).  The group velocity points from the circle toward its centre.
    Assumptions: as :func:`rossby_omega`, U = 0.   Validation: V1 points on the circle return ω.  Label: analytic.
    """
    w = float(omega)
    if w == 0:
        raise ValueError("rossby_omega_circle: omega must be non-zero")
    F = 0.0 if np.isinf(c) else (float(f0) / float(c)) ** 2
    r2 = (beta / (2.0 * w)) ** 2 - F
    return dict(center_k=-beta / (2.0 * w), center_l=0.0, radius=float(np.sqrt(max(r2, 0.0))), exists=bool(r2 > 0))



def rossby_long_wave_speed(beta, f0, c):
    """Phase (and group) speed of long, non-dispersive Rossby waves: c_x ≈ −β c²/f₀² = −β Λ².

    Book: §13.15 (unnumbered limit of Eq. (13.119) for k² + l² → 0).
    Parameters: beta [1/(m s)]; f0 [1/s] non-zero; c [m/s].   Returns c_x [m/s] (negative: westward).
    Assumptions: wavelength ≫ 2πΛ.   Validation: V7 limit of :func:`rossby_phase_speed`.  Label: analytic.
    """
    return _S(-_F(beta) * _F(c) ** 2 / _nonzero_f(f0, "rossby_long_wave_speed") ** 2)


def stationary_rossby_wavelength(U, beta):
    """Wavelength of the barotropic Rossby wave held stationary by an eastward flow: λ = 2π sqrt(U/β).

    Book: §13.15 (unnumbered, from Eq. (13.120) with c_x = 0, l = 0, f₀²/c² negligible) and §13.13 (flow over a step).
    Parameters: U eastward mean flow [m/s] > 0; beta [1/(m s)] > 0.
    Returns lambda [m].
    Raises ValueError for U ≤ 0: a westward flow cannot hold a Rossby wave still (its phase already moves west).
    Assumptions: barotropic, l = 0.   Validation: V1 ω = 0 at k = sqrt(β/U).  Label: analytic.
    """
    if np.any(_F(U) <= 0) or np.any(_F(beta) <= 0):
        raise ValueError("stationary_rossby_wavelength: needs an eastward flow U > 0 and beta > 0")
    return _S(2.0 * np.pi * np.sqrt(_F(U) / _F(beta)))



def rossby_packet_spectrum(k0, sigma, n):
    """Wavenumbers and Gaussian weights of a wave packet centred on k₀ (to build a Rossby packet as a sum of modes).

    Book: not in the book — a helper of ours for the packet figures of §13.15.
    Parameters: k0 central wavenumber [rad/m]; sigma spectral width [rad/m] > 0; n number of modes (odd keeps k₀).
    Returns (k, amplitudes): n wavenumbers evenly spaced over k₀ ± 3σ [rad/m] and weights exp(−(k − k₀)²/2σ²)
    normalised to sum 1, so that the packet's peak height is the sum of the amplitudes.
    Assumptions: none.   Validation: V1 weights sum to 1; envelope width 1/σ.  Label: analytic.
    """
    if sigma <= 0 or int(n) < 1:
        raise ValueError("rossby_packet_spectrum: sigma > 0 and n >= 1 are required")
    k = k0 + sigma * np.linspace(-3.0, 3.0, int(n))
    a = np.exp(-0.5 * ((k - k0) / sigma) ** 2)
    return k, a / a.sum()



# =====================================================================================================================
# 7. Potential vorticity and the instabilities of the very long waves (§13.13, §13.16, §13.17)
# =====================================================================================================================

def potential_vorticity(zeta, f, h):
    """Shallow-water potential vorticity q = (ζ + f)/h, conserved following the motion.

    Book: §13.13, Eq. (13.94).
    Parameters: zeta relative vorticity ∂v/∂x − ∂u/∂y [1/s]; f Coriolis parameter at the column's latitude [1/s];
    h layer thickness [m] > 0.
    Returns q [1/(m s)].
    Assumptions: one homogeneous hydrostatic layer, inviscid (the conservation law is exact for shallow water with
    any f(y); the book reaches it through the f → f₀ replacement and then restores f, trap T10).
    Validation: V2 D q/Dt = 0 from (13.88)–(13.90) by sympy; V4 constant on particles of the numerical model.
    Label: symbolic.
    """
    return _S((_F(zeta) + _F(f)) / _F(h))  # Eq. (13.94)


def step_vorticity(f, h0, h1):
    """Relative vorticity a column acquires on crossing a step from depth h₀ to h₁: ζ = f (h₁ − h₀)/h₀.

    Book: §13.13 (unnumbered, from f/h₀ = (ζ + f)/h₁ for a column with no relative vorticity upstream).
    Parameters: f at the upstream latitude [1/s]; h0, h1 depths before and after [m].
    Returns zeta [1/s]: of sign opposite to f (anticyclonic) when the layer gets shallower.
    Assumptions: as :func:`potential_vorticity`.   Validation: V7 equals ch05 ``column_relative_vorticity``.
    Label: analytic.
    """
    return _S(_F(f) * (_F(h1) - _F(h0)) / _F(h0))


def absolute_vorticity_gradient(y, U, beta, h: float | None = None):
    """Northward gradient of the absolute vorticity of a zonal current, d(ζ̄ + f)/dy = β − d²U/dy².

    Book: §13.16, Eq. (13.124).
    Parameters
    ----------
    y : northward coordinate [m] (or non-dimensional) — the sample points (≥ 4, increasing).
    U : the current — samples U(y_j) on the grid, or a callable U(y) (then differentiated by a central second
        difference of step ``h``, default 1e-4 of the span of y).
    beta : planetary vorticity gradient, in units consistent with U and y.
    Returns
    -------
    array β − U″ at every point of y.
    Numerics (samples): 3-point second differences inside; at the two ends the second derivative of the cubic through
    the four nearest nodes — second order throughout.
    Assumptions: barotropic zonal flow U(y).   Validation: V1 on U = sech²y sampled at 401 points over |y| ≤ 5 the
    largest error against the analytic β − U″ is 8.3e-4.  Label: analytic.
    """
    y = _F(y)
    if callable(U):
        step = float(h) if h is not None else 1e-4 * max(float(np.ptp(y)), 1.0)
        Upp = (_F(U(y + step)) - 2.0 * _F(U(y)) + _F(U(y - step))) / step ** 2
    else:
        Uy = _F(U)
        if y.size < 4:
            raise ValueError("absolute_vorticity_gradient: need at least 4 samples")
        Upp = _d2(y, Uy)
        for end, sl in ((0, slice(0, 4)), (-1, slice(-4, None))):
            s = (y[sl] - y[end]) / (y[sl][-1] - y[sl][0])          # scaled, origin at the end node
            cf = np.polyfit(s, Uy[sl], 3)
            Upp[end] = 2.0 * cf[1] / (y[sl][-1] - y[sl][0]) ** 2
    return beta - Upp  # Eq. (13.124)



def rayleigh_kuo_criterion(y, U, beta, tol: float = 1e-12) -> dict:
    """Rayleigh–Kuo necessary condition for barotropic instability: β − U″ must change sign inside the flow.

    Book: §13.16, Eq. (13.124) and the statement around it (Rayleigh's inflection-point criterion of ch11 with −U″
    replaced by β − U″).
    Parameters: y grid (increasing); U samples or callable (see :func:`absolute_vorticity_gradient`); beta; tol —
    relative threshold (of max|β − U″|) below which a value counts as zero.
    Returns dict: ``changes_sign`` (bool — the necessary condition is met), ``y_zero`` (list of linearly interpolated
    zero crossings), ``min`` and ``max`` of β − U″ on the grid.  **Necessary, not sufficient**: a profile that meets
    it can still be stable; one that fails it is stable.
    Assumptions: inviscid, barotropic, parallel flow, β-plane.
    Validation: V1 analytic jets with β on either side of max U″; V7 β = 0 equals ch11 ``rayleigh_criterion``.
    Label: analytic.
    """
    y = _F(y)
    G = _F(absolute_vorticity_gradient(y, U, beta))
    scale = max(float(np.max(np.abs(G))), 1e-300)
    s = np.sign(np.where(np.abs(G) <= tol * scale, 0.0, G))
    idx = np.nonzero(s[:-1] * s[1:] < 0)[0]
    y0 = y[idx] - G[idx] * (y[idx + 1] - y[idx]) / (G[idx + 1] - G[idx])
    return dict(changes_sign=bool(np.any(s > 0) and np.any(s < 0)), y_zero=[float(v) for v in y0], min=float(G.min()),
                max=float(G.max()))



def eady_alpha(k, l, N, f):
    """Vertical decay rate of quasi-geostrophic perturbations: α = (N/|f|) sqrt(k² + l²).

    Book: §13.17, Eq. (13.139).  (This α is neither the thermal expansion coefficient of §13.3 nor the enstrophy flux
    of §13.18 — trap T13.)
    Parameters: k, l horizontal wavenumbers [rad/m]; N [rad/s]; f [1/s] non-zero.   Returns alpha [1/m] ≥ 0.
    Assumptions: uniform N, f-plane.   Validation: V1 αH = K·(NH/|f|).  Label: analytic.
    """
    return _S(_F(N) / np.abs(_nonzero_f(f, "eady_alpha")) * np.sqrt(_F(k) ** 2 + _F(l) ** 2))  # Eq. (13.139)


def _eady_brackets(x):
    """x − tanh x and x − coth x for x = αH/2 ≥ 0, with series below 1e-3 where the first loses digits."""
    x = _F(x)
    small = x < 1e-3
    xs = np.where(small, 1.0, x)
    t1 = np.where(small, x ** 3 / 3.0 - 2.0 * x ** 5 / 15.0 + 17.0 * x ** 7 / 315.0, xs - np.tanh(xs))
    with np.errstate(divide="ignore", invalid="ignore"):
        t2 = np.where(x > 0, x - 1.0 / np.tanh(np.where(x > 0, x, 1.0)), -np.inf)
        prod = np.where(small, -x ** 2 / 3.0 + 16.0 * x ** 4 / 45.0, t1 * np.where(x > 0, t2, 0.0))
    return t1, t2, prod


def eady_factors(alphaH) -> dict:
    """The two brackets under the square root of the Eady phase speed and their product.

    Book: §13.17, Eq. (13.141): (αH/2 − tanh(αH/2)) and (αH/2 − coth(αH/2)).
    Parameters: alphaH = α H (dimensionless) ≥ 0.
    Returns dict: ``tanh_factor`` = x − tanh x (always ≥ 0), ``coth_factor`` = x − coth x (negative below the
    cut-off; −inf at αH = 0), ``product`` (negative ⇔ unstable; → −(αH)²/12 as αH → 0), with x = αH/2.
    Numerics: series below x = 1e-3 (the first bracket is a difference of nearly equal numbers there).
    Assumptions: Eady problem.   Validation: V1 against direct evaluation in extended precision.  Label: analytic.
    """
    t1, t2, prod = _eady_brackets(0.5 * _F(alphaH))
    return dict(tanh_factor=_S(t1), coth_factor=_S(t2), product=_S(prod))



def eady_phase_speed(alphaH, U0, root: int = +1):
    """Complex eastward phase speed of Eady waves.

    Book: §13.17, Eq. (13.141): c = U₀/2 ± (U₀/αH) sqrt[(αH/2 − tanh(αH/2))(αH/2 − coth(αH/2))].
    Parameters: alphaH = α H (dimensionless) ≥ 0; U0 velocity at the upper lid [m/s] (the basic flow is U₀ z/H);
    root +1 (the growing mode below the cut-off, the faster edge wave above it) or −1.
    Returns c [m/s], complex.  Below the cut-off c = U₀/2 ± i c_i; above it c is real and lies between 0 and U₀.
    At αH = 0 the limit U₀/2 ± i U₀/(2 sqrt 3) is returned.
    Numerics: the radicand is evaluated as a product (never joined across the cut-off by interpolation).
    Assumptions: Eady problem (uniform N and shear, f-plane, rigid lids at z = 0, H, quasi-geostrophic).
    Validation: V2 root of the 2 × 2 determinant; V3 Chebyshev eigenvalues (``ch13.eady_numeric_eigs``).
    Label: symbolic.
    """
    a = _F(alphaH)
    _, _, prod = _eady_brackets(0.5 * a)
    with np.errstate(divide="ignore", invalid="ignore"):
        rad = np.where(a > 0, np.sqrt(prod.astype(complex)) / np.where(a > 0, a, 1.0), 1j / (2.0 * np.sqrt(3.0)))
    return _S(0.5 * U0 + (1 if root >= 0 else -1) * U0 * rad)  # Eq. (13.141)


def eady_growth_rate(k, l, N, f, H, U0):
    """Growth rate of an Eady wave in physical units: σ = |k| c_i [1/s].

    Book: §13.17 (the perturbation grows as exp(k c_i t), with c from Eq. (13.141)); the dimensional evaluation is
    **ours — not in the book**.
    Parameters: k eastward and l northward wavenumber [rad/m]; N [rad/s]; f [1/s] non-zero; H layer depth [m];
    U0 velocity difference across the layer [m/s].
    Returns sigma [1/s] ≥ 0; exactly 0 at and beyond the cut-off (αH ≥ 2.3994) and for k = 0.
    Assumptions: as :func:`eady_phase_speed`.   Validation: V1 maximum equals :func:`eady_max_growth_rate` for l = 0.
    Label: analytic.
    """
    aH = _F(eady_alpha(k, l, N, f)) * H
    ci = np.abs(np.imag(_F(0.0) + eady_phase_speed(aH, abs(U0))))
    return _S(np.abs(_F(k)) * ci)


def eady_critical() -> float:
    """Cut-off of the Eady instability: the root α_c H of αH/2 = coth(αH/2), about 2.3994.

    Book: §13.17 (unnumbered marginal condition; the page prints the root rounded to two figures).
    Returns alpha_c H (dimensionless) — computed with ``brentq`` (xtol 1e-14), never typed.
    Assumptions: Eady problem.   Validation: V1 residual of x = coth x.  Label: analytic.
    """
    return 2.0 * brentq(lambda x: x - 1.0 / np.tanh(x), 1.0, 2.0, xtol=1e-14, rtol=1e-14)


def eady_fastest() -> dict:
    """The fastest-growing Eady wave (l = 0), found numerically from Eq. (13.141).

    Book: §13.17 states the most unstable wavelength as a multiple of Λ = NH/f without deriving it ("it can be
    shown"); the numbers here are **ours, computed**.
    Returns dict (all dimensionless): ``alphaH`` (1.6061), ``sigma_nd`` = σ_max N H/(|f| U₀) = the maximum of
    αH·c_i/U₀ (0.30982), ``ci_over_U0`` at the maximum (0.19290), ``cr_over_U0`` (0.5: the wave moves with the
    mid-level flow).  The wavelengths follow as 2π/αH in units of Λ (:func:`eady_wavelengths`).
    Numerics: bounded Brent maximisation (xatol 1e-12) on (0, α_c H).
    Assumptions: Eady problem, l = 0, Λ ≡ NH/|f| (the π-less Rossby radius of §13.17).
    Validation: V1 dσ/dα = 0 by finite differences; V3 against the Chebyshev eigenvalues.  Label: analytic.
    """
    ac = eady_critical()

    def neg(a):
        return -a * abs(complex(eady_phase_speed(a, 1.0)).imag)

    res = minimize_scalar(neg, bounds=(1e-6, ac - 1e-9), method="bounded", options=dict(xatol=1e-12))
    a = float(res.x)
    c = complex(eady_phase_speed(a, 1.0))
    return dict(alphaH=a, sigma_nd=float(-res.fun), ci_over_U0=float(c.imag), cr_over_U0=float(c.real))



def eady_wavelengths() -> dict:
    """Wavelengths of the Eady problem in units of the Rossby radius Λ = NH/|f| (l = 0): the short-wave cut-off and the
    fastest-growing wave.

    Book: §13.17 (the two unnumbered multiples of Λ after Eq. (13.142), printed there to two figures; here computed).
    Returns dict (dimensionless): ``cutoff_over_Lambda`` = 2π/(α_c H) (2.6187: shorter waves are stable) and
    ``fastest_over_Lambda`` = 2π/αH at the maximum growth rate (3.9120).
    Assumptions: Eady problem, l = 0.   Validation: V1 from :func:`eady_critical` and :func:`eady_fastest`.
    Label: analytic.
    """
    return dict(cutoff_over_Lambda=2.0 * np.pi / eady_critical(), fastest_over_Lambda=2.0 * np.pi / eady_fastest()["alphaH"])



def eady_max_growth_rate(f, N, dUdz):
    """Largest Eady growth rate in physical units: σ_max = 0.30982 |f| |dU/dz| / N.

    **Call by keyword** — ``eady_max_growth_rate(f=…, N=…, dUdz=…)``: the order (f, N, dUdz) differs from the
    (N, f) order of ``eady_alpha`` and ``eady_growth_rate``, and swapping f and N changes the answer by (N/f)².

    Book: **ours — not in the book** (the storm-track "Eady growth rate"); the coefficient is the computed maximum of
    Eq. (13.141), ``eady_fastest()["sigma_nd"]``, not a typed constant.
    Parameters: f [1/s]; N [rad/s] > 0; dUdz vertical shear U₀/H [1/s].   Returns sigma_max [1/s].
    Assumptions: Eady problem.   Validation: V1 maximum of :func:`eady_growth_rate` over k.  Label: analytic.
    """
    return _S(eady_fastest()["sigma_nd"] * np.abs(_F(f)) * np.abs(_F(dUdz)) / _F(N))



def eady_time_scale(f, N, dUdz):
    """e-folding time of the fastest Eady wave, 1/σ_max [s].

    **Call by keyword** — ``eady_time_scale(f=…, N=…, dUdz=…)``: the order is (f, N, dUdz), not (N, f).

    Book: **ours — not in the book** (see :func:`eady_max_growth_rate`).
    Parameters: f [1/s]; N [rad/s]; dUdz [1/s].   Returns the time [s] (inf for zero shear or f = 0).
    Assumptions: Eady problem.   Validation: V1 reciprocal.  Label: analytic.
    """
    s = _F(eady_max_growth_rate(f, N, dUdz))
    with np.errstate(divide="ignore"):
        return _S(np.where(s > 0, 1.0 / np.where(s > 0, s, 1.0), np.inf))


def _eady_AB(alphaH, U0, root=+1):
    """Coefficients (A, B) of p̂ = A cosh α(z − H/2) + B sinh α(z − H/2) from the condition at z = 0, scaled so that
    the largest |p̂| in the layer is 1; also returns c."""
    a = float(alphaH)
    if a <= 0:
        raise ValueError("the Eady mode needs alphaH > 0")
    c = complex(eady_phase_speed(a, U0, root)) / U0
    q = 0.5 * a
    # lower lid: c dp/dz + (U0/H) p = 0  ->  A[-a c sinh q + cosh q] + B[a c cosh q - sinh q] = 0
    B = (a * c * np.sinh(q) - np.cosh(q)) / (a * c * np.cosh(q) - np.sinh(q))
    s = a * (np.linspace(0.0, 1.0, 201) - 0.5)
    peak = float(np.max(np.abs(np.cosh(s) + B * np.sinh(s))))
    return (1.0 + 0j) / peak, B / peak, c * U0



def eady_mode(z, alphaH, U0, H, root: int = +1) -> dict:
    """Vertical structure of the pressure perturbation of an Eady wave.

    Book: §13.17, Eq. (13.140) p̂ = A cosh α(z − H/2) + B sinh α(z − H/2), with B/A from the boundary condition at
    the lower lid (the first of the two unnumbered homogeneous equations before Eq. (13.141)).
    Parameters: z height above the lower lid [m], 0 ≤ z ≤ H; alphaH (dimensionless) > 0; U0 [m/s]; H [m];
    root as :func:`eady_phase_speed` (+1: the growing mode below the cut-off).
    Returns dict: ``p_hat`` (complex, normalised so that the largest |p̂| in the layer is 1), ``amplitude`` |p̂|,
    ``phase`` arg p̂ [rad] (unwrapped for an array), ``A``, ``B`` (complex).  With p′ = Re[p̂ exp ik(x − ct)] a phase
    that increases with height means the pressure pattern tilts westward with height (for k > 0) — the signature of
    the growing mode.
    Assumptions: Eady problem.   Validation: V2 both lid conditions; V1 tilt westward for root = +1.  Label: analytic.
    """
    z = _F(z)
    A, B, _ = _eady_AB(alphaH, U0, root)
    s = alphaH * (z / H - 0.5)
    p = A * np.cosh(s) + B * np.sinh(s)  # Eq. (13.140)
    ph = np.angle(p)
    if ph.ndim:
        ph = np.unwrap(ph)
    return dict(p_hat=_S(p), amplitude=_S(np.abs(p)), phase=_S(ph), A=complex(A), B=complex(B))



def eady_vertical_velocity(z, alphaH, U0, H, k, N, f, rho0, root: int = +1):
    """Complex amplitude of the vertical velocity of an Eady wave of unit pressure amplitude.

    Book: §13.17, Eq. (13.135) with p′ ∝ exp ik(x − ct): ŵ = −(ik/ρ₀N²)[(U − c) dp̂/dz − (U₀/H) p̂], U = U₀ z/H.
    Parameters: z [m]; alphaH; U0 [m/s]; H [m]; k eastward wavenumber [rad/m]; N [rad/s]; f [1/s] (unused by the
    formula itself; kept so that the argument list matches :func:`eady_fluxes`); rho0 [kg/m³]; root as
    :func:`eady_phase_speed`.
    Returns w_hat (complex) [m/s per Pa of pressure amplitude] for the mode of :func:`eady_mode` (max |p̂| = 1).
    Zero at both lids.
    Assumptions: Eady problem.   Validation: V1 ŵ(0) = ŵ(H) = 0.  Label: analytic.
    """
    z = _F(z)
    A, B, c = _eady_AB(alphaH, U0, root)
    al = alphaH / H
    s = al * (z - 0.5 * H)
    p = A * np.cosh(s) + B * np.sinh(s)
    pz = al * (A * np.sinh(s) + B * np.cosh(s))
    return _S(-1j * float(k) / (rho0 * N ** 2) * ((U0 * z / H - c) * pz - U0 / H * p))  # Eq. (13.135)



def eady_fluxes(z, alphaH, U0, H, N, f, rho0, g: float = G0, k=None, root: int = +1) -> dict:
    """Eddy density fluxes of an Eady wave, averaged over a wavelength: mean(w′ρ′) and mean(v′ρ′).

    Book: §13.17 "Energetics": dK/dt = −g ∫ w′ρ′ dV (unnumbered) — growth needs light fluid rising and dense fluid
    sinking; with v′ from Eq. (13.131), ρ′ from Eq. (13.134) and w′ from Eq. (13.135).
    Parameters: z [m]; alphaH; U0 [m/s]; H [m]; N [rad/s]; f [1/s] non-zero; rho0 [kg/m³]; g [m/s²]; k eastward
    wavenumber [rad/m] (default αH·|f|/(N H): the wave with l = 0); root as :func:`eady_phase_speed`.
    Returns dict: ``w_rho`` (mean w′ρ′ at z; negative at every interior height for the growing mode: potential
    energy is released), ``v_rho`` (mean v′ρ′ at z), ``v_T_sign`` (+1.0 when the eddy heat flux, proportional to
    −v′ρ′, is poleward — the case for the growing mode in either hemisphere; −1.0 equatorward; 0.0 for a neutral
    wave), ``phase_tilt`` [rad] (phase of p̂ at the upper lid minus that at the lower lid; positive = westward tilt
    with height).  Fluxes are per unit squared pressure amplitude (the mode has max |p̂| = 1).
    Assumptions: Eady problem; averages of products of two waves, mean(a′b′) = ½ Re(â b̂*).
    Validation: V1 signs; V7 fluxes vanish for the neutral waves beyond the cut-off.  Label: analytic.
    """
    z = _F(z)
    fa = float(_nonzero_f(f, "eady_fluxes"))
    A, B, c = _eady_AB(alphaH, U0, root)
    al = alphaH / H
    kk = alphaH * abs(fa) / (N * H) if k is None else float(k)

    def fields(zz):
        s = al * (zz - 0.5 * H)
        p = A * np.cosh(s) + B * np.sinh(s)
        pz = al * (A * np.sinh(s) + B * np.cosh(s))
        w = -1j * kk / (rho0 * N ** 2) * ((U0 * zz / H - c) * pz - U0 / H * p)
        return p, w, 1j * kk * p / (rho0 * fa), -pz / g      # p, w, v (13.131), rho (13.134)

    p, w, v, r = fields(z)
    _, _, vm, rm = fields(np.array(0.5 * H))
    vr_mid = 0.5 * float(np.real(vm * np.conj(rm)))
    scale = 0.5 * abs(vm) * abs(rm)
    sgn = 0.0 if abs(vr_mid) <= 1e-12 * max(scale, 1e-300) else float(np.sign(-np.sign(fa) * vr_mid))
    p0, pH = fields(np.array(0.0))[0], fields(np.array(float(H)))[0]
    return dict(w_rho=_S(0.5 * np.real(w * np.conj(r))), v_rho=_S(0.5 * np.real(v * np.conj(r))), v_T_sign=sgn,
                phase_tilt=float(np.angle(pH / p0)))



# =====================================================================================================================
# 8. Two-dimensional (geostrophic) turbulence (§13.18)
# =====================================================================================================================

def fjortoft_transfer(K0, K1, K2, S0: float = 1.0) -> dict:
    """Fjørtoft's argument: energy S₀ at wavenumber K₀ moved to K₁ < K₀ < K₂ with energy and enstrophy conserved.

    Book: §13.18, the two conservation sums S₀ = S₁ + S₂, K₀²S₀ = K₁²S₁ + K₂²S₂ and Eq. (13.145).
    Parameters: K0, K1, K2 wavenumbers [rad/m] with K1 < K0 < K2; S0 initial spectral energy (any unit).
    Returns dict: ``S1``, ``S2`` (the two shares, same unit as S0), ``energy_ratio`` = S₁/S₂ =
    (K₂² − K₀²)/(K₀² − K₁²), ``enstrophy_ratio`` = K₁²S₁/(K₂²S₂).  For K₁ and K₂ a factor apart on either side of K₀
    more energy goes to the larger scale and more enstrophy to the smaller.
    Assumptions: inviscid two-dimensional flow; transfer to exactly two wavenumbers (a thought experiment).
    Validation: V1 direct 2 × 2 solve; both sums conserved.  Label: analytic.
    """
    K0, K1, K2 = float(K0), float(K1), float(K2)
    if not K1 < K0 < K2:
        raise ValueError("fjortoft_transfer: needs K1 < K0 < K2")
    S1 = S0 * (K2 ** 2 - K0 ** 2) / (K2 ** 2 - K1 ** 2)
    S2 = S0 * (K0 ** 2 - K1 ** 2) / (K2 ** 2 - K1 ** 2)
    return dict(S1=S1, S2=S2, energy_ratio=S1 / S2, enstrophy_ratio=K1 ** 2 * S1 / (K2 ** 2 * S2))  # Eq. (13.145)


def enstrophy_spectrum(K, S):
    """Enstrophy spectrum K² S(K) of a two-dimensional isotropic field with energy spectrum S(K).

    Book: §13.18 (unnumbered: mean ζ² = ∫₀^∞ K² S(K) dK, with mean u² = ∫₀^∞ S(K) dK — one-sided in K, no factor ½;
    trap T16).
    Parameters: K [rad/m]; S energy spectrum [m³/s²].   Returns K²S [m/s²].
    Assumptions: isotropic, homogeneous, two-dimensional, non-divergent.
    Validation: V1 integral equals the mean-square vorticity of a synthetic field.  Label: analytic.
    """
    return _S(_F(K) ** 2 * _F(S))


def two_d_cascade_spectrum(K, K0, eps, alpha_ens, C_E: float = 1.0, C_Z: float = 1.0):
    """Model energy spectrum of forced two-dimensional turbulence: ε^{2/3} K^{−5/3} below the forcing wavenumber K₀
    (inverse energy cascade) and α^{2/3} K^{−3} above it (forward enstrophy cascade).

    Book: §13.18 (the two inertial ranges, S ∝ ε^{2/3}K^{−5/3} and S ∝ α^{2/3}K^{−3}).  The two prefactors are **not
    given by the book**: C_E = C_Z = 1 are placeholders of order one, not measured constants.
    Parameters: K [rad/m] > 0; K0 forcing wavenumber [rad/m]; eps energy injection rate [m²/s³]; alpha_ens enstrophy
    flux [1/s³]; C_E, C_Z dimensionless.
    Returns S(K) [m³/s²] (discontinuous at K₀ unless α = ε K₀²·(C_E/C_Z)^{3/2}).
    Assumptions: stationary, isotropic, both ranges inertial.
    Validation: V1 exponents from the Π theorem (exact rationals).  Label: analytic (dimensional).
    """
    K = _F(K)
    return _S(np.where(K < K0, C_E * eps ** (2.0 / 3.0) * K ** (-5.0 / 3.0), C_Z * alpha_ens ** (2.0 / 3.0) * K ** (-3.0)))


def rhines_length(u_rms, beta):
    """Rhines length sqrt(u/β): the eddy size at which two-dimensional turbulence turns into Rossby waves.

    Book: §13.18 (unnumbered "Rhines length"; the original paper is unread by us — no citation is implied).
    Parameters: u_rms root-mean-square eddy speed [m/s]; beta [1/(m s)] > 0.   Returns the length [m].
    Assumptions: an order-of-magnitude scale (no numerical factor is implied).
    Validation: V1 identity; qualitative comparison with the jet spacing of the β-plane run.  Label: analytic.
    """
    return _S(np.sqrt(np.abs(_F(u_rms)) / _F(beta)))
