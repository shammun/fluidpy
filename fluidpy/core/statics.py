"""Fluid statics: gauge pressure, hydrostatic balance, isothermal and layered atmospheres, buoyancy.

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e, §1.7 (Eqs. 1.6–1.9, Example 1.1) and §1.10 (static profiles,
isothermal atmosphere, scale height). Reused by Ch. 4 (hydrostatic base state), Ch. 7 (layered fluids), Ch. 13.

Conventions: **z is positive upward** (§1.7, §1.10); depth below a level is h = −z. Pressures are absolute unless the
name says ``gauge``. ``p0`` is the pressure at z = 0 (Eq. 1.9); it is not the potential-temperature reference
pressure ``p_ref`` of ``core.stratification``.
"""
from __future__ import annotations

from typing import Callable, Sequence

import numpy as np
from scipy.integrate import solve_ivp

from ._util import as_scalar_if_0d, require_positive
from .thermo import G0, P_ATM, R_AIR

# U.S. Standard Atmosphere 1976 defining constants (NASA-TM-X-74335).
USSA_T0: float = 288.15        #: sea-level temperature [K]
USSA_P0: float = 101325.0      #: sea-level pressure [Pa]
USSA_LAPSE0: float = -6.5e-3   #: dT/dz of the 0–11 km layer [K/m] (book sign convention: negative = cooling upward)
#: Base geopotential altitudes of the seven layers [m] and their lapse rates dT/dz [K/m] (book sign).
USSA_BASES = np.array([0.0, 11000.0, 20000.0, 32000.0, 47000.0, 51000.0, 71000.0])
USSA_LAPSE = np.array([-6.5e-3, 0.0, 1.0e-3, 2.8e-3, 0.0, -2.8e-3, -2.0e-3])
USSA_Z_TOP: float = 84852.0    #: top of the seventh layer, geopotential altitude [m]
USSA_Z1: float = 11000.0       #: top of the first layer (tropopause base), geopotential altitude [m]


def gauge_pressure(p, p_atm=P_ATM):
    """Gauge pressure: absolute pressure minus atmospheric pressure.

    Book: §1.7, ``p_gauge = p − p_atm``.

    Parameters
    ----------
    p : float or array_like
        Absolute pressure [Pa].
    p_atm : float, optional
        Atmospheric pressure [Pa]; default 101325 Pa (standard atmosphere).

    Returns
    -------
    p_gauge : float or ndarray
        [Pa]; 0 at atmospheric pressure, −p_atm in a vacuum.

    Validation (planned): V1 p_atm → 0, vacuum → −p_atm. Label: pending.
    """
    return as_scalar_if_0d(np.asarray(p, dtype=float) - p_atm)  # §1.7: p_gauge = p − p_atm


def absolute_pressure(p_gauge, p_atm=P_ATM):
    """Absolute pressure from a gauge reading: ``p = p_gauge + p_atm``.

    Book: §1.7 (inverse of the gauge definition).

    Parameters
    ----------
    p_gauge : float or array_like
        Gauge pressure [Pa].
    p_atm : float, optional
        Atmospheric pressure [Pa].

    Returns
    -------
    p : float or ndarray
        [Pa].

    Validation (planned): V1 round trip with :func:`gauge_pressure`. Label: pending.
    """
    return as_scalar_if_0d(np.asarray(p_gauge, dtype=float) + p_atm)


def hydrostatic_pressure_uniform(z, p0, rho, g=G0):
    """Pressure in a fluid of uniform density at rest.

    Book: §1.7, Eq. (1.9) ``p = p0 − rho g z`` (integral of (1.8) with rho constant).

    Parameters
    ----------
    z : float or array_like
        Height above the reference level, positive upward [m].
    p0 : float
        Pressure at z = 0 [Pa].
    rho : float
        Uniform density [kg/m^3].
    g : float, optional
        Gravitational acceleration [m/s^2].

    Returns
    -------
    p : float or ndarray
        [Pa]; at depth h = −z the rise is rho g h.

    Notes
    -----
    Assumptions: static fluid, uniform rho and g.

    Validation (planned): V2 sympy dp/dz + rho g = 0; V1 p(0) = p0 and depth h gives +rho g h. Label: pending.
    """
    return as_scalar_if_0d(p0 - rho * g * np.asarray(z, dtype=float))  # Eq. (1.9)


def integrate_hydrostatic(z, rho_fn: Callable[[float, float], float], p0: float, g: float = G0,
                          z0: float | None = None, rtol: float = 1e-10, atol: float | None = None,
                          method: str = "DOP853"):
    """Integrate the hydrostatic law for a density that may depend on height and pressure.

    Book: §1.7, Eq. (1.8) ``dp/dz = −rho g``.

    Parameters
    ----------
    z : array_like
        Heights where p is wanted [m], z upward, any order.
    rho_fn : callable
        ``rho(z, p)`` [kg/m^3] with z [m] and p [Pa] (e.g. ``lambda z, p: p/(R*T(z))`` for a perfect gas).
    p0 : float
        Pressure at ``z0`` [Pa].
    g : float, optional
        Gravitational acceleration [m/s^2].
    z0 : float, optional
        Level where p = p0 [m]; default ``z[0]``.
    rtol : float, optional
        Relative tolerance of ``solve_ivp`` (default 1e-10).
    atol : float, optional
        Absolute tolerance [Pa]; default ``1e-12 * max(|p0|, 1)``.
    method : str, optional
        ``solve_ivp`` method (default ``"DOP853"``).

    Returns
    -------
    p : ndarray
        Pressure at ``z`` [Pa].

    Notes
    -----
    Method: adaptive Runge–Kutta integration upward and downward from z0 with dense output (our choice).
    Assumptions: static fluid; gravity uniform and vertical.

    Validation (planned): V1 constant rho equals Eq. (1.9); isothermal gas equals :func:`isothermal_pressure`; linear
    T(z) equals :func:`linear_lapse_pressure`; V3 error falls with rtol. Label: pending.
    """
    z = np.atleast_1d(np.asarray(z, dtype=float))
    z0 = float(z[0]) if z0 is None else float(z0)
    atol = 1e-12 * max(abs(p0), 1.0) if atol is None else atol

    def rhs(zz, p):
        return [-rho_fn(zz, p[0]) * g]  # Eq. (1.8): dp/dz = −rho g

    p = np.empty_like(z)
    up = z >= z0
    down = ~up
    for mask, z_end in ((up, z[up].max() if up.any() else None), (down, z[down].min() if down.any() else None)):
        if z_end is None:
            continue
        if z_end == z0:
            p[mask] = p0
            continue
        sol = solve_ivp(rhs, (z0, z_end), [p0], method=method, rtol=rtol, atol=atol, dense_output=True)
        if not sol.success:
            raise RuntimeError(f"hydrostatic integration failed: {sol.message}")
        p[mask] = sol.sol(z[mask])[0]
    return p


def atmosphere_from_temperature(z, T_fn: Callable[[np.ndarray], np.ndarray], p0: float = P_ATM, R: float = R_AIR,
                                g: float = G0, z0: float | None = None, rtol: float = 1e-10):
    """Static perfect-gas column from a prescribed temperature profile: p(z), rho(z), T(z).

    Book: §1.10 (hydrostatics (1.8) plus the equation of state (1.12)/(1.22) tie p, rho and T together, so one profile
    fixes the other two).

    Parameters
    ----------
    z : array_like
        Heights [m], upward.
    T_fn : callable
        ``T(z)`` [K], vectorised.
    p0 : float
        Pressure at ``z0`` [Pa].
    R : float, optional
        Gas constant [J/(kg K)].
    g : float, optional
        Gravitational acceleration [m/s^2].
    z0 : float, optional
        Level of p0 [m]; default ``z[0]``.
    rtol : float, optional
        Integration tolerance.

    Returns
    -------
    (p, rho, T) : tuple of ndarray
        [Pa], [kg/m^3], [K].

    Notes
    -----
    Assumptions: static, perfect gas, uniform g.

    Validation (planned): V1 isothermal case equals :func:`isothermal_pressure`; V5 USSA-1976 0–11 km layer reproduces
    the tabulated pressures at 5 and 11 km. Label: pending.
    """
    z = np.atleast_1d(np.asarray(z, dtype=float))
    p = integrate_hydrostatic(z, lambda zz, pp: pp / (R * float(T_fn(zz))), p0, g=g, z0=z0, rtol=rtol)  # (1.8)+(1.22)
    T = np.asarray(T_fn(z), dtype=float) * np.ones_like(z)
    rho = p / (R * T)  # Eq. (1.22)
    return p, rho, T


def isothermal_pressure(z, p0, T, R=R_AIR, g=G0):
    """Pressure in an isothermal perfect-gas atmosphere.

    Book: §1.10, "Scale Height of the Atmosphere": ``dp/dz = −rho g = −p g/(R T)`` integrates to
    ``p(z) = p0 exp(−g z/(R T))``.

    Parameters
    ----------
    z : float or array_like
        Height [m], upward.
    p0 : float
        Pressure at z = 0 [Pa].
    T : float
        Uniform temperature [K].
    R : float, optional
        Gas constant [J/(kg K)].
    g : float, optional
        Gravitational acceleration [m/s^2].

    Returns
    -------
    p : float or ndarray
        [Pa].

    Notes
    -----
    Assumptions: isothermal, perfect gas, constant g and R.

    Validation (planned): V2 sympy ODE residual; V1 p(H)/p0 = e^-1 with H = R T/g. Label: pending.
    """
    require_positive("T", T)
    return as_scalar_if_0d(p0 * np.exp(-g * np.asarray(z, dtype=float) / (R * T)))  # §1.10: p = p0 e^{−gz/RT}


def isothermal_density(z, rho0, T, R=R_AIR, g=G0):
    """Density in an isothermal perfect-gas atmosphere: ``rho = rho0 exp(−g z/(R T))`` (follows p by (1.22)).

    Book: §1.10 (isothermal atmosphere; rho = p/(R T) with T fixed).

    Parameters
    ----------
    z : float or array_like
        Height [m].
    rho0 : float
        Density at z = 0 [kg/m^3].
    T : float
        Temperature [K].
    R, g : float, optional
        Gas constant [J/(kg K)], gravity [m/s^2].

    Returns
    -------
    rho : float or ndarray
        [kg/m^3].

    Validation (planned): V1 equals isothermal_pressure/(R T). Label: pending.
    """
    require_positive("T", T)
    return as_scalar_if_0d(rho0 * np.exp(-g * np.asarray(z, dtype=float) / (R * T)))


def scale_height(T, R=R_AIR, g=G0):
    """Scale height of an isothermal atmosphere: the height over which pressure falls by a factor e.

    Book: §1.10, ``H = R T / g``.

    Parameters
    ----------
    T : float or array_like
        Temperature [K].
    R : float, optional
        Gas constant [J/(kg K)].
    g : float, optional
        Gravitational acceleration [m/s^2].

    Returns
    -------
    H : float or ndarray
        [m].

    Notes
    -----
    Assumptions: isothermal perfect gas. The same law follows from dimensional analysis (Example 1.2) with constant 1.

    Validation (planned): V1 H = R T/g and linear in T; V6 book's value at its mean temperature (private JSON).
    Label: pending.
    """
    return as_scalar_if_0d(R * np.asarray(T, dtype=float) / g)  # §1.10: H = RT/g


def linear_lapse_pressure(z, p0, T0, dT_dz, R=R_AIR, g=G0):
    """Pressure in a perfect-gas layer whose temperature changes linearly with height, T = T0 + (dT/dz) z.

    Book: §1.10 (hydrostatics (1.8) with (1.22) for a prescribed lapse rate Γ ≡ dT/dz; closed form derived by us:
    ``p = p0 (T/T0)^(−g/(R Γ))``, reducing to the isothermal law as Γ → 0).

    Parameters
    ----------
    z : float or array_like
        Height above the level of p0 and T0 [m].
    p0 : float
        Pressure at z = 0 [Pa].
    T0 : float
        Temperature at z = 0 [K].
    dT_dz : float
        Lapse rate Γ = dT/dz [K/m] in the **book's sign**: negative when T falls with height (USSA troposphere −6.5e-3).
    R, g : float, optional
        Gas constant [J/(kg K)], gravity [m/s^2].

    Returns
    -------
    p : float or ndarray
        [Pa].

    Notes
    -----
    Assumptions: static perfect gas, constant g; T must stay positive over z.

    Validation (planned): V1 equals :func:`integrate_hydrostatic`; V7 dT_dz → 0 recovers :func:`isothermal_pressure`;
    V5 USSA-1976 table at 5 and 11 km. Label: pending.
    """
    z = np.asarray(z, dtype=float)
    if dT_dz == 0.0:
        return isothermal_pressure(z, p0, T0, R, g)
    T = T0 + dT_dz * z
    require_positive("T along the layer", T)
    return as_scalar_if_0d(p0 * (T / T0) ** (-g / (R * dT_dz)))  # integral of dp/p = −g dz/(R (T0 + Γ z))


def standard_atmosphere(z, g=G0, R=R_AIR):
    """U.S. Standard Atmosphere 1976 below 84.852 km geopotential altitude: temperature, pressure, density.

    Book: §1.10 (a static atmosphere is fixed by its temperature profile through (1.8) and (1.22); the isothermal
    scale-height model is compared with it). Data model: NASA-TM-X-74335 (1976) — seven layers with linear T(z):
    base geopotential altitudes 0, 11, 20, 32, 47, 51, 71 km and lapse rates dT/dz = −6.5, 0, +1.0, +2.8, 0, −2.8,
    −2.0 K/km (book sign), T = 288.15 K and p = 101325 Pa at sea level (layer table as reproduced by J. Hawley,
    "Formulae and code for the U.S. Standard Atmosphere (1976)"). Base temperatures and pressures follow by closed-form
    integration layer by layer.

    Parameters
    ----------
    z : float or array_like
        Geopotential altitude [m], 0 <= z <= 84852.
    g : float, optional
        Standard gravity [m/s^2].
    R : float, optional
        Gas constant [J/(kg K)] (USSA uses R* = 8314.32 J/(kmol K); our CODATA-based R_AIR differs by ~2e-5).

    Returns
    -------
    (T, p, rho) : tuple of float or ndarray
        [K], [Pa], [kg/m^3].

    Notes
    -----
    Assumptions: hydrostatic dry perfect gas with constant molecular weight; geopotential (not geometric) altitude.

    Validation (planned): V5 PDAS USSA-1976 table values at 0, 5, 11, 20, 32, 47 km (1e-3 rel). Label: pending.
    """
    z = np.asarray(z, dtype=float)
    if np.any(z < 0) or np.any(z > USSA_Z_TOP):
        raise ValueError(f"standard_atmosphere covers 0 <= z <= {USSA_Z_TOP:.0f} m (geopotential)")
    Tb = [USSA_T0]
    pb = [USSA_P0]
    for i in range(len(USSA_BASES) - 1):
        dz = USSA_BASES[i + 1] - USSA_BASES[i]
        Tb.append(Tb[i] + USSA_LAPSE[i] * dz)
        pb.append(linear_lapse_pressure(dz, pb[i], Tb[i], USSA_LAPSE[i], R, g))
    Tb, pb = np.array(Tb), np.array(pb)
    k = np.clip(np.searchsorted(USSA_BASES, z, side="right") - 1, 0, len(USSA_BASES) - 1)
    dz = z - USSA_BASES[k]
    L = USSA_LAPSE[k]
    T = Tb[k] + L * dz  # linear temperature in each layer
    with np.errstate(divide="ignore", invalid="ignore"):
        p_lin = pb[k] * (T / Tb[k]) ** (-g / (R * np.where(L == 0, 1.0, L)))  # (1.8)+(1.22), dT/dz ≠ 0
    p_iso = pb[k] * np.exp(-g * dz / (R * Tb[k]))  # (1.8)+(1.22), isothermal layer
    p = np.where(L == 0, p_iso, p_lin)
    rho = p / (R * T)  # Eq. (1.22)
    return as_scalar_if_0d(T), as_scalar_if_0d(p), as_scalar_if_0d(rho)


def layered_pressure(z, thicknesses: Sequence[float], densities: Sequence[float], p_surface: float = P_ATM,
                     g: float = G0):
    """Hydrostatic pressure in a stack of uniform-density liquid layers below a free surface.

    Book: §1.7, Eq. (1.9) applied layer by layer (pressure is continuous across each interface; within a layer
    ``p = p_top − rho_i g (z − z_top)``).

    Parameters
    ----------
    z : float or array_like
        Height [m] with z = 0 at the free surface and z < 0 inside the liquid (depth h = −z).
    thicknesses : sequence of float
        Layer thicknesses [m], listed from the top layer down.
    densities : sequence of float
        Layer densities [kg/m^3], same order; below the last layer the last density continues.
    p_surface : float, optional
        Pressure at the free surface [Pa] (default standard atmosphere).
    g : float, optional
        [m/s^2].

    Returns
    -------
    p : float or ndarray
        Absolute pressure [Pa]; for z > 0 (above the surface) ``p_surface`` is returned (air weight neglected).

    Notes
    -----
    Assumptions: static, immiscible layers of uniform density, uniform g.

    Validation (planned): V1 one layer equals :func:`hydrostatic_pressure_uniform`; slope −rho_i g in each layer;
    continuity at interfaces. Label: pending.
    """
    th = np.asarray(thicknesses, dtype=float)
    rh = np.asarray(densities, dtype=float)
    if th.size != rh.size or th.size == 0:
        raise ValueError("thicknesses and densities must be non-empty and of equal length")
    z_int = -np.concatenate(([0.0], np.cumsum(th)))  # interface heights, top first (0, −h1, −h1−h2, …)
    p_int = p_surface + g * np.concatenate(([0.0], np.cumsum(rh * th)))  # pressure at each interface
    z = np.asarray(z, dtype=float)
    k = np.clip(np.searchsorted(-z_int, -z, side="right") - 1, 0, th.size - 1)  # layer index containing z
    p = p_int[k] - rh[k] * g * (z - z_int[k])  # Eq. (1.9) within layer k
    return as_scalar_if_0d(np.where(z > 0, p_surface, p))


def buoyancy_force(rho_fluid, volume, g=G0):
    """Archimedes' buoyancy: upward net pressure force on a fully submerged body.

    Book: §1.7 (consequence of Eq. (1.9); derivation D37 added by the curation: the pressure difference between the
    bottom and top faces of a body in a fluid at rest equals the weight of the displaced fluid).

    Parameters
    ----------
    rho_fluid : float or array_like
        Density of the surrounding fluid [kg/m^3] (uniform).
    volume : float or array_like
        Submerged volume [m^3].
    g : float, optional
        Gravitational acceleration [m/s^2].

    Returns
    -------
    F_b : float or ndarray
        Upward force [N] (positive = up).

    Notes
    -----
    Assumptions: static uniform fluid; the body's own density is irrelevant to this force.

    Validation (planned): V1 equals the z-component of :func:`net_pressure_force_on_box` for Eq. (1.9). Label: pending.
    """
    return as_scalar_if_0d(np.asarray(rho_fluid, dtype=float) * g * np.asarray(volume, dtype=float))


def net_pressure_force_on_box(p_fn: Callable[[np.ndarray, np.ndarray, np.ndarray], np.ndarray],
                              box: Sequence[float], n: int = 41):
    """Net force exerted by a pressure field on a rectangular box: ``F = −∮ p n dA``, face by face.

    Book: §1.7 (the force balances on a fluid cube behind Eqs. (1.7)–(1.8), Fig. 1.6) and derivation D37 (buoyancy).

    Parameters
    ----------
    p_fn : callable
        ``p(x, y, z)`` [Pa], vectorised over numpy arrays (z upward).
    box : sequence of 6 floats
        ``(x0, x1, y0, y1, z0, z1)`` [m] with x0 < x1, y0 < y1, z0 < z1.
    n : int, optional
        Gauss–Legendre points per direction on each face (exact for polynomial p of degree <= 2n − 1).

    Returns
    -------
    F : ndarray, shape (3,)
        Force on the box [N] (x, y, z components); for Eq. (1.9) it is (0, 0, rho g V).

    Notes
    -----
    Method: tensor-product Gauss–Legendre quadrature on the six faces (our choice).

    Validation (planned): V1 horizontal components vanish (Pascal, Eq. 1.7) and F_z = rho g V for Eq. (1.9). Label: pending.
    """
    x0, x1, y0, y1, z0, z1 = map(float, box)
    xi, wi = np.polynomial.legendre.leggauss(n)

    def nodes(a, b):
        return 0.5 * (b - a) * xi + 0.5 * (a + b), 0.5 * (b - a) * wi

    xs, wx = nodes(x0, x1)
    ys, wy = nodes(y0, y1)
    zs, wz = nodes(z0, z1)

    def face_integral(A, B, wA, wB, make_xyz):
        AA, BB = np.meshgrid(A, B, indexing="ij")
        X, Y, Z = make_xyz(AA, BB)
        return float(np.sum(np.outer(wA, wB) * p_fn(X, Y, Z)))

    Fx = (face_integral(ys, zs, wy, wz, lambda a, b: (np.full_like(a, x0), a, b))
          - face_integral(ys, zs, wy, wz, lambda a, b: (np.full_like(a, x1), a, b)))
    Fy = (face_integral(xs, zs, wx, wz, lambda a, b: (a, np.full_like(a, y0), b))
          - face_integral(xs, zs, wx, wz, lambda a, b: (a, np.full_like(a, y1), b)))
    Fz = (face_integral(xs, ys, wx, wy, lambda a, b: (a, b, np.full_like(a, z0)))  # bottom face pushes up
          - face_integral(xs, ys, wx, wy, lambda a, b: (a, b, np.full_like(a, z1))))  # top face pushes down
    return np.array([Fx, Fy, Fz])


__all__ = [
    "USSA_T0", "USSA_P0", "USSA_LAPSE0", "USSA_Z1", "USSA_BASES", "USSA_LAPSE", "USSA_Z_TOP", "layered_pressure",
    "gauge_pressure", "absolute_pressure", "hydrostatic_pressure_uniform", "integrate_hydrostatic",
    "atmosphere_from_temperature", "isothermal_pressure", "isothermal_density", "scale_height",
    "linear_lapse_pressure", "standard_atmosphere", "buoyancy_force", "net_pressure_force_on_box",
]
