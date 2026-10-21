"""Atmospheric balances on pressure surfaces: geostrophic wind from geopotential height on the sphere, the isothermal
scale height and the hypsometric thickness of a layer.

Vallis: §2.6, Eq. (2.166) (hypsometric equation; scale height RT/g in the text after it); §2.8, Eqs. (2.198),
(2.210) and (2.215) (geostrophic balance in spherical and pressure coordinates, f × u_g = −∇_p Φ).
Kundu: §13.5, Eqs. (K 13.11)–(K 13.12) (geostrophic velocity, ``fluidpy.core.gfd``); §1.10 (scale height,
``fluidpy.core.statics.scale_height``).

Conventions
-----------
* Gridded fields are (..., lat, lon): longitude on the last axis, latitude on the second-to-last. Latitude and longitude
  coordinates (``lat``, ``lon``) are given in **degrees** at this interface (as in every reanalysis file) and
  converted to radians inside. Latitude may run S→N or N→S (PSL files run N→S); the derivative uses the actual
  coordinate values, so either order is correct.
* Geopotential height Z = Φ/g₀ with g₀ = 9.80665 m s⁻² (WMO; the NCEP/ERA5 definition).
* Pressures in Pa inside functions; temperatures in K.
"""
from __future__ import annotations

import numpy as np

from ..core._util import as_scalar_if_0d
from ..core.statics import scale_height as _core_scale_height
from .constants import EARTH_RADIUS, G0, OMEGA, R_D
from .frames import coriolis_parameter

__all__ = ["geostrophic_wind_from_height", "scale_height", "hypsometric_thickness"]

_S = as_scalar_if_0d


def _ddlon_periodic(Z: np.ndarray, dlam: float) -> np.ndarray:
    """Centred difference ∂Z/∂λ on a uniform, globally periodic longitude axis (last axis)."""
    return (np.roll(Z, -1, axis=-1) - np.roll(Z, 1, axis=-1)) / (2.0 * dlam)


def _ddlon_open(Z: np.ndarray, lam: np.ndarray) -> np.ndarray:
    """∂Z/∂λ on a regional (non-periodic) longitude axis: centred inside, second-order one-sided at the ends."""
    return np.gradient(Z, lam, axis=-1, edge_order=2)


def geostrophic_wind_from_height(Z, lat, lon, g0: float = G0, Omega: float = OMEGA, a: float = EARTH_RADIUS,
                                 lat_min: float = 5.0, periodic: bool | None = None):
    """Geostrophic wind on a pressure surface from its geopotential height Z(ϑ, λ) on the sphere:
    u_g = −(g/f) ∂Z/∂y, v_g = (g/f) ∂Z/∂x, with dy = a dϑ and dx = a cos ϑ dλ.

    Vallis: §2.8, Eq. (2.215) f × u_g = −∇_p Φ with Φ = gZ, in the spherical form of Eqs. (2.198)/(2.210a,b):
    f u_g = −(1/a) ∂Φ/∂ϑ, f v_g = (1/(a cos ϑ)) ∂Φ/∂λ. Kundu: §13.5, Eqs. (K 13.11)–(K 13.12) (Cartesian form,
    ``core.gfd.geostrophic_from_height``).

    Parameters
    ----------
    Z : array_like, shape (..., nlat, nlon) — geopotential height [m] on one pressure surface (leading axes, e.g.
        time or level, broadcast).
    lat : array_like, shape (nlat,) — latitudes [degrees], strictly monotonic (either direction).
    lon : array_like, shape (nlon,) — longitudes [degrees], uniformly spaced and increasing.
    g0 : float — the g₀ that defines Z [m s⁻²] (default 9.80665).
    Omega : float — rotation rate [rad/s].   a : float — planetary radius [m].
    lat_min : float — |ϑ| below which the result is NaN [degrees] (default 5°; f → 0 there and geostrophy
        fails). Maps of the full sphere (e.g. the 500 hPa scene) pass ``lat_min=15`` to hide the near-equatorial
        blow-up of 1/f; the library default stays 5°.
    periodic : bool or None — wrap the longitude difference (a global grid). ``None`` (default) decides from the grid:
        periodic when nlon·Δλ = 360° to 1e-6.

    Returns
    -------
    (u_g, v_g) : ndarrays like Z [m/s] — eastward and northward geostrophic wind at the grid nodes. NaN within
    ±``lat_min`` of the equator and at the poles (|ϑ| = 90°, where cos ϑ = 0).

    Numerics (our choice; the book gives no scheme): centred second-order differences in λ (periodic for a global grid,
    one-sided second-order at the ends of a regional grid) and in ϑ (``np.gradient`` with the actual coordinates,
    second-order at the first and last latitude); f = 2Ω sin ϑ at each node (not an f-plane).
    Assumptions: hydrostatic, Ro ≪ 1, thin shell (r ≈ a), traditional approximation.
    Validation: V1 exact (round-off) for fields quadratic in (ϑ, λ), edge nodes included; the solid-body test
    Z = Z₀ − (Ω a u₀/g) sin²ϑ gives u_g = u₀ cos ϑ, v_g = 0 to second order in the grid spacing (≈ 0.025 m s⁻¹ at
    2.5°, u₀ = 20 m s⁻¹); V2 sympy (solid body, wave, f k×u_g = −∇Φ); V3 second-order convergence; V7 cyclonic (anticlockwise seen from above) flow round a NH low, clockwise round an
    SH low; northward wind east of a NH trough. Label: converged.
    """
    Z = np.asarray(Z, dtype=float)
    lat_r = np.deg2rad(np.asarray(lat, dtype=float))
    lon_d = np.asarray(lon, dtype=float)
    if Z.ndim < 2 or Z.shape[-2:] != (lat_r.size, lon_d.size):
        raise ValueError(f"Z must have shape (..., {lat_r.size}, {lon_d.size}); got {Z.shape}")
    if lon_d.size < 3 or lat_r.size < 3:
        raise ValueError("need at least 3 latitudes and 3 longitudes")
    dlon = np.diff(lon_d)
    if not np.allclose(dlon, dlon[0], rtol=0, atol=1e-6) or dlon[0] <= 0:
        raise ValueError("longitudes must be uniformly spaced and increasing")
    if periodic is None:
        periodic = abs(lon_d.size * dlon[0] - 360.0) < 1e-6
    lam = np.deg2rad(lon_d)
    dZ_dlam = _ddlon_periodic(Z, np.deg2rad(dlon[0])) if periodic else _ddlon_open(Z, lam)
    dZ_dlat = np.gradient(Z, lat_r, axis=-2, edge_order=2)
    f = np.asarray(coriolis_parameter(lat_r, Omega))[:, None]  # Eq. (2.80): f = 2Ω sin ϑ, column over latitude
    coslat = np.cos(lat_r)[:, None]
    bad = (np.abs(lat_r) < np.deg2rad(lat_min)) | (np.abs(np.abs(lat_r) - np.pi / 2) < 1e-9)
    with np.errstate(divide="ignore", invalid="ignore"):
        u_g = -(g0 / (f * a)) * dZ_dlat  # Eq. (2.198)/(2.215): f u_g = −(g₀/a) ∂Z/∂ϑ
        v_g = (g0 / (f * a * coslat)) * dZ_dlam  # Eq. (2.198)/(2.215): f v_g = (g₀/(a cos ϑ)) ∂Z/∂λ
    mask = np.broadcast_to(bad[:, None], Z.shape[-2:])
    u_g = np.where(mask, np.nan, u_g)
    v_g = np.where(mask, np.nan, v_g)
    return u_g, v_g


def scale_height(T, R: float = R_D, g: float = G0):
    """Isothermal (pressure) scale height H = R T / g — the height over which pressure falls by a factor e.

    Vallis: §2.6, text after Eq. (2.166) (RT/g ≈ 7 km at 240 K, 8.2 km at 280 K); §1.6 (isothermal atmosphere).
    Kundu: §1.10 — calls ``fluidpy.core.statics.scale_height`` (one implementation).

    Parameters
    ----------
    T : float or array_like — temperature [K].   R : float — gas constant [J kg⁻¹ K⁻¹] (dry air default).
    g : float — gravity [m s⁻²].

    Returns
    -------
    H : float or ndarray [m].

    Assumptions: ideal gas, isothermal layer, hydrostatic.
    Validation: V1 linear in T; ≈ 7.0 km at 240 K and ≈ 8.2 km at 280 K (Vallis's numbers, V6). Label: analytic.
    """
    return _core_scale_height(T, R, g)  # H = R T / g


def hypsometric_thickness(p1, p2, T, p=None, R: float = R_D, g: float = G0):
    """Thickness z(p₂) − z(p₁) of the layer between two pressure levels: the hypsometric equation
    z(p₂) − z(p₁) = −(R/g) ∫_{p₁}^{p₂} T d ln p.

    Vallis: §2.6, Eq. (2.166). Kundu: §1.10 (hydrostatic ideal-gas atmosphere; not numbered there).

    Parameters
    ----------
    p1, p2 : float or array_like — bounding pressures [Pa] (p₁ > p₂ gives a positive thickness: p₂ is the upper level).
    T : float or array_like — either the layer-mean temperature [K] (when ``p`` is None; broadcasts with p1, p2), or a
        temperature profile [K] at the pressure levels ``p`` (1-D).
    p : array_like or None — pressure levels of the profile T [Pa] (1-D, any order, must span [p₂, p₁]).
    R : float — gas constant [J kg⁻¹ K⁻¹].   g : float — gravity [m s⁻²].

    Returns
    -------
    Δz : float or ndarray [m]. For a mean temperature: (R T̄/g) ln(p₁/p₂). For a profile: trapezoidal integral of T in
    ln p, with T linearly interpolated in ln p at p₁ and p₂.

    Assumptions: hydrostatic, ideal gas, constant g (Φ = g z).
    Validation: V1 isothermal closed form; 1000→100 hPa at 240 K ≈ 16 km (Vallis's rule of thumb, V6); V1 a profile
    with T linear in ln p is integrated exactly (the trapezoid in ln p is exact for it, off-grid bounds included);
    V3 second-order convergence when T is nonlinear in ln p (USSA-1976 troposphere); V5 USSA-1976 5 km thickness.
    Label: analytic.
    """
    if p is None:
        Tb = np.asarray(T, dtype=float)
        return _S(R * Tb / g * np.log(np.asarray(p1, dtype=float) / np.asarray(p2, dtype=float)))  # Eq. (2.166)
    pp = np.asarray(p, dtype=float)
    TT = np.asarray(T, dtype=float)
    if pp.ndim != 1 or TT.shape != pp.shape:
        raise ValueError("a profile needs 1-D p and T of the same length")
    p1f, p2f = float(p1), float(p2)
    lo, hi = min(p1f, p2f), max(p1f, p2f)
    if lo < pp.min() * (1 - 1e-12) or hi > pp.max() * (1 + 1e-12):
        raise ValueError("p1 and p2 must lie within the profile's pressure range")
    order = np.argsort(pp)
    lnp, Ts = np.log(pp[order]), TT[order]
    inner = (lnp > np.log(lo)) & (lnp < np.log(hi))
    x = np.concatenate([[np.log(lo)], lnp[inner], [np.log(hi)]])
    y = np.interp(x, lnp, Ts)
    integral = np.trapezoid(y, x) if hasattr(np, "trapezoid") else np.trapz(y, x)  # ∫_{lo}^{hi} T d ln p
    # Eq. (2.166): z(p2) − z(p1) = −(R/g) ∫_{p1}^{p2} T d ln p ; ∫_{p1}^{p2} = +∫_{lo}^{hi} if p1 < p2 else −∫
    sign = 1.0 if p1f < p2f else -1.0
    return float(-(R / g) * sign * integral)
