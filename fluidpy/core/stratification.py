"""Static stability of stratified fluids: displaced parcels, Brunt–Väisälä frequency, lapse rates, potential
temperature and density, the ocean criterion.

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e, §1.10, Eqs. (1.29)–(1.35) and the unnumbered parcel equation.
Reused by Ch. 4 (Boussinesq N), Ch. 7 (internal waves, omega <= N), Ch. 11 (Richardson number), Ch. 13 (GFD).

Sign conventions (binding user decision, the book's own)
-------------------------------------------------------
* z is positive upward; ζ is the upward displacement of a parcel from its rest height z_o.
* The **lapse rate is Γ ≡ dT/dz** (book §1.10). It is *negative* where temperature falls with height. Every
  environment lapse-rate argument is named ``dT_dz``. :func:`adiabatic_lapse_rate` returns Γ_a = −g α T / C_p, which
  is **negative** (≈ −9.8e-3 K/m for dry air). A layer is stable when dT/dz > Γ_a (equivalently dθ/dz > 0).
  ⚠️ Most meteorology (including the lapse-rate-feedback literature) uses the opposite, standard convention
  Γ ≡ −dT/dz: there Γ_a ≈ +9.8 K/km and a layer is stable when Γ < Γ_a. Same physics, flipped sign and inequality.
  Every computation here uses Kundu's sign; the meteorology form is produced, never hidden, by
  :func:`lapse_rate_convention` (Γ in either convention) and :func:`lapse_rate_stability` (verdict plus the criterion
  written out with numbers in either convention). See knowledge/notation.md.
* ``p_ref`` (default 1.0e5 Pa) is the reference pressure p_o of (1.31)/(1.33); it is kept separate from the surface
  pressure ``p0`` of ``core.statics``.
"""
from __future__ import annotations

from typing import Callable, NamedTuple

import numpy as np
from scipy.integrate import solve_ivp

from ._util import as_scalar_if_0d, require_positive
from .thermo import CP_AIR, G0, GAMMA_AIR, P_REF, R_AIR


# --------------------------------------------------------------------------------------------------------------------
# Brunt–Väisälä frequency and classification
# --------------------------------------------------------------------------------------------------------------------
def brunt_vaisala_sq(rho0, drho_dz, drho_a_dz, g=G0):
    """Square of the Brunt–Väisälä (buoyancy) frequency from the environment and isentropic density gradients.

    Book: §1.10, Eq. (1.29) ``N^2 = −(g/rho(z_o)) (drho/dz − drho_a/dz)``.

    Parameters
    ----------
    rho0 : float or array_like
        Density of the medium at the rest height z_o [kg/m^3].
    drho_dz : float or array_like
        Equilibrium (environment) density gradient [kg/m^4].
    drho_a_dz : float or array_like
        Isentropic density gradient followed by a displaced parcel [kg/m^4] (0 for an incompressible parcel;
        −rho g/c^2 for a compressible one, see :func:`isentropic_density_gradient`).
    g : float, optional
        Gravitational acceleration [m/s^2].

    Returns
    -------
    N2 : float or ndarray
        [1/s^2]; > 0 stable (oscillation at N), < 0 unstable (growth rate sqrt(−N^2)), = 0 neutral.

    Notes
    -----
    Assumptions: small displacement (first order in ζ), frictionless adiabatic parcel in instant pressure equilibrium
    with its surroundings, static background, no mixing.

    Validation (planned): V2 sympy isothermal perfect gas gives g^2/(c_p T); V4 equals (g/θ) dθ/dz for random smooth
    T(z); V1 constant rho with incompressible parcel gives 0. Label: pending.
    """
    rho0 = np.asarray(rho0, dtype=float)
    require_positive("rho0", rho0)
    N2 = -(g / rho0) * (np.asarray(drho_dz, dtype=float) - np.asarray(drho_a_dz, dtype=float))  # Eq. (1.29)
    return as_scalar_if_0d(N2)


def brunt_vaisala_sq_from_theta(theta, dtheta_dz, g=G0):
    """N^2 of a perfect-gas atmosphere from the potential-temperature gradient: ``N^2 = (g/θ) dθ/dz``.

    Book: §1.10, links Eq. (1.29) with Eqs. (1.31)–(1.34) (derivation D36, added by the curation; the book states only
    that the sign of dθ/dz decides stability).

    Parameters
    ----------
    theta : float or array_like
        Potential temperature [K].
    dtheta_dz : float or array_like
        Its vertical gradient [K/m].
    g : float, optional
        Gravitational acceleration [m/s^2].

    Returns
    -------
    N2 : float or ndarray
        [1/s^2].

    Notes
    -----
    Assumptions: perfect gas with constant γ, hydrostatic background, small displacements.

    Validation (planned): V4 equals :func:`brunt_vaisala_sq` with drho_a_dz = −rho g/c^2 on smooth profiles. Label: pending.
    """
    theta = np.asarray(theta, dtype=float)
    require_positive("theta", theta)
    return as_scalar_if_0d(g / theta * np.asarray(dtheta_dz, dtype=float))  # D36: N^2 = (g/θ) dθ/dz


def brunt_vaisala_sq_from_lapse(T, dT_dz, cp=CP_AIR, g=G0):
    """N^2 of a perfect-gas atmosphere from temperature and lapse rate: ``N^2 = (g/T) (dT/dz + g/C_p)``.

    Book: §1.10, Eq. (1.32) ``(T/θ) dθ/dz = dT/dz + g/C_p = Γ − Γ_a`` combined with N^2 = (g/θ) dθ/dz (D36).

    Parameters
    ----------
    T : float or array_like
        Temperature [K].
    dT_dz : float or array_like
        Environment lapse rate in Kundu's convention Γ = dT/dz [K/m] (negative when T falls with height;
        the meteorology lapse rate is −dT_dz).
    cp : float, optional
        Specific heat at constant pressure [J/(kg K)].
    g : float, optional
        Gravitational acceleration [m/s^2].

    Returns
    -------
    N2 : float or ndarray
        [1/s^2]; zero when dT/dz = Γ_a = −g/C_p.

    Validation (planned): V1 dT/dz = −g/cp gives 0; V2 equals :func:`brunt_vaisala_sq_from_theta` via (1.32). Label: pending.
    """
    T = np.asarray(T, dtype=float)
    require_positive("T", T)
    return as_scalar_if_0d(g / T * (np.asarray(dT_dz, dtype=float) + g / cp))  # (1.32) with N^2 = (g/θ) dθ/dz


def classify_stability(N2, tol: float = 1e-12):
    """Label a static state stable, neutral or unstable from the sign of N^2.

    Book: §1.10 (after Eq. (1.29): N^2 > 0 stable, N^2 < 0 unstable, N^2 = 0 neutral).

    Parameters
    ----------
    N2 : float or array_like
        Squared buoyancy frequency [1/s^2].
    tol : float, optional
        Absolute tolerance [1/s^2]: |N2| <= tol counts as neutral (round-off guard; typical |N^2| ~ 1e-4).

    Returns
    -------
    label : str or ndarray of str
        ``"stable"``, ``"neutral"`` or ``"unstable"`` (array of labels for array input).

    Validation (planned): V1 sign cases; V7 tolerance boundary and array handling. Label: pending.
    """
    arr = np.asarray(N2, dtype=float)
    lab = np.where(arr > tol, "stable", np.where(arr < -tol, "unstable", "neutral"))
    return str(lab) if lab.ndim == 0 else lab


def stability_timescale(N2, tol: float = 1e-12):
    """Time scale of a displaced parcel: oscillation period (stable) or e-folding time (unstable).

    Book: §1.10 (solutions of the parcel equation: cos(N t) for N^2 > 0, cosh(sqrt(−N^2) t) for N^2 < 0).

    Parameters
    ----------
    N2 : float
        Squared buoyancy frequency [1/s^2].
    tol : float, optional
        Neutral tolerance [1/s^2] (as :func:`classify_stability`).

    Returns
    -------
    (kind, seconds) : tuple (str, float)
        ``("period", 2π/N)`` for N^2 > 0, ``("efold", 1/sqrt(−N^2))`` for N^2 < 0, ``("none", inf)`` when neutral.

    Validation (planned): V1 period 2π/N matches the zero crossings of :func:`parcel_displacement`. Label: pending.
    """
    N2 = float(N2)
    label = classify_stability(N2, tol)
    if label == "stable":
        return "period", float(2.0 * np.pi / np.sqrt(N2))  # cos(N t) repeats after 2π/N
    if label == "unstable":
        return "efold", float(1.0 / np.sqrt(-N2))  # cosh(σ t) grows by e every 1/σ
    return "none", float("inf")


# --------------------------------------------------------------------------------------------------------------------
# Displaced parcel
# --------------------------------------------------------------------------------------------------------------------
def parcel_displacement(t, zeta0, N2, w0=0.0):
    """Analytic displacement of a parcel obeying the linear parcel equation ``ζ'' + N^2 ζ = 0``.

    Book: §1.10, unnumbered parcel equation ``d^2ζ/dt^2 − (g/rho(z_o)) (drho/dz − drho_a/dz) ζ = 0`` (derivation D18,
    Exercise 1.13), whose coefficient is −N^2 of Eq. (1.29).

    Parameters
    ----------
    t : float or array_like
        Time since release [s].
    zeta0 : float
        Initial upward displacement [m].
    N2 : float or array_like
        Squared buoyancy frequency [1/s^2].
    w0 : float, optional
        Initial vertical velocity [m/s] (default 0: release from rest, as in the book).

    Returns
    -------
    zeta : float or ndarray
        Displacement [m]: ``ζ0 cos Nt + (w0/N) sin Nt`` (N^2 > 0), ``ζ0 + w0 t`` (N^2 = 0),
        ``ζ0 cosh σt + (w0/σ) sinh σt`` with σ = sqrt(−N^2) (N^2 < 0).

    Notes
    -----
    Assumptions: linearised (small ζ), undamped (the book's remark that viscosity and conduction arrest the oscillation
    is outside this model).

    Validation (planned): V2 sympy residual of ζ'' + N^2 ζ = 0; V1 period 2π/N; V7 N^2 → 0 limit. Label: pending.
    """
    t = np.asarray(t, dtype=float)
    N2 = np.asarray(N2, dtype=float)
    a = np.sqrt(np.abs(N2))
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        stable = zeta0 * np.cos(a * t) + np.where(a > 0, w0 / np.where(a > 0, a, 1.0), 0.0) * np.sin(a * t)
        unstable = zeta0 * np.cosh(a * t) + np.where(a > 0, w0 / np.where(a > 0, a, 1.0), 0.0) * np.sinh(a * t)
    neutral = zeta0 + w0 * t
    zeta = np.where(N2 > 0, stable, np.where(N2 < 0, unstable, neutral))  # solutions of ζ'' + N^2 ζ = 0
    return as_scalar_if_0d(zeta)


def parcel_ode(t_span, zeta0: float, rho_env_fn: Callable[[float], float], rho_parcel_fn: Callable[[float], float],
               z0: float = 0.0, w0: float = 0.0, g: float = G0, t_eval=None, zeta_max: float | None = None,
               rtol: float = 1e-10, atol: float = 1e-12):
    """Nonlinear motion of a displaced parcel under its weight and buoyancy.

    Book: §1.10, Newton's second law "including weight and buoyancy for the displaced element" (Fig. 1.8), before
    linearisation: ``rho_p V ζ'' = −rho_p V g + rho_e V g``, i.e. ``ζ'' = −g (rho_p(ζ) − rho_e(z0 + ζ)) / rho_p(ζ)``.

    Parameters
    ----------
    t_span : (float, float)
        Start and end time [s].
    zeta0 : float
        Initial upward displacement [m].
    rho_env_fn : callable
        Environment density ``rho_e(z)`` [kg/m^3] at absolute height z [m].
    rho_parcel_fn : callable
        Parcel density ``rho_p(ζ)`` [kg/m^3] as a function of its displacement ζ [m] from z0 (it follows the isentrope
        through its rest state).
    z0 : float, optional
        Rest height z_o [m].
    w0 : float, optional
        Initial vertical velocity [m/s].
    g : float, optional
        Gravitational acceleration [m/s^2].
    t_eval : array_like, optional
        Output times [s]; default the solver's own steps.
    zeta_max : float, optional
        Stop when |ζ| reaches this value [m] (terminal event: the cap for N^2 < 0 runaway, which would otherwise grow
        like cosh(sqrt(−N^2) t) until overflow). Default ``1e3 * max(|zeta0|, 1e-3)`` (≈ 7.6 e-folding times from rest).
        A second terminal event stops the run if either density reaches zero.
    rtol, atol : float, optional
        ``solve_ivp`` tolerances (RK45).

    Returns
    -------
    (t, zeta) : tuple of ndarray
        Times [s] and displacement [m] (shorter than ``t_eval`` if the cap was hit).

    Notes
    -----
    Assumptions: frictionless adiabatic parcel, pressure equal to its surroundings at every height, no mixing.
    Buoyancy = weight of displaced fluid (Archimedes, derivation D37).

    Validation (planned): V3 small ζ0 matches :func:`parcel_displacement` with error ∝ ζ0^2; V1 period 2π/N. Label: pending.
    """
    cap = 1e3 * max(abs(zeta0), 1e-3) if zeta_max is None else float(zeta_max)

    def rhs(_t, y):
        rp = rho_parcel_fn(y[0])
        return [y[1], -g * (rp - rho_env_fn(z0 + y[0])) / rp]  # Newton II with weight and buoyancy

    def hit_cap(_t, y):
        return abs(y[0]) - cap

    def density_vanishes(_t, y):  # guard: stop before a linear density model reaches zero (non-physical)
        return min(rho_parcel_fn(y[0]), rho_env_fn(z0 + y[0]))

    hit_cap.terminal = True
    density_vanishes.terminal = True
    sol = solve_ivp(rhs, t_span, [zeta0, w0], t_eval=t_eval, events=(hit_cap, density_vanishes), rtol=rtol,
                    atol=atol)
    if sol.status == -1:
        raise RuntimeError(f"parcel integration failed: {sol.message}")
    return sol.t, sol.y[0]


def _sample_parcel(t, zeta0, rho_env_fn, rho_parcel_fn, g, w0, zeta_max):
    """Run :func:`parcel_ode` from t = 0 and sample at the requested times (NaN after a runaway cap)."""
    t_arr = np.atleast_1d(np.asarray(t, dtype=float))
    t_end = float(t_arr.max())
    if t_end <= 0.0:
        out = np.full_like(t_arr, float(zeta0))
    else:
        ts, zs = parcel_ode((0.0, t_end), zeta0, rho_env_fn, rho_parcel_fn, z0=0.0, w0=w0, g=g,
                            t_eval=np.unique(np.concatenate(([0.0], t_arr))), zeta_max=zeta_max)
        out = np.interp(t_arr, ts, zs) if ts.size > 1 else np.full_like(t_arr, float(zeta0))
        out = np.where(t_arr <= ts[-1] + 1e-12, out, np.nan)
    return as_scalar_if_0d(out[0]) if np.ndim(t) == 0 else out


def parcel_ode_from_gradients(t, zeta0, rho0, drho_dz, drho_a_dz, g=G0, w0=0.0, zeta_max=None):
    """Nonlinear parcel displacement for linear environment and parcel density profiles (ocean / lab-tank mode).

    Book: §1.10 parcel argument with the book's first-order densities kept exact: environment
    ``rho(z_o + ζ) = rho0 + (drho/dz) ζ`` and parcel ``rho0 + (drho_a/dz) ζ``; ``ζ'' = −g (rho_p − rho_e)/rho_p``.

    Parameters
    ----------
    t : float or array_like
        Output time(s) [s], >= 0.
    zeta0 : float
        Initial displacement [m].
    rho0 : float
        Density at the rest height [kg/m^3].
    drho_dz : float
        Environment density gradient [kg/m^4].
    drho_a_dz : float
        Parcel (isentropic) density gradient [kg/m^4].
    g : float, optional
        Gravitational acceleration [m/s^2].
    w0 : float, optional
        Initial velocity [m/s].
    zeta_max : float, optional
        Runaway cap [m] (see :func:`parcel_ode`); times after the cap return NaN.

    Returns
    -------
    zeta : float or ndarray
        Displacement at ``t`` [m].

    Notes
    -----
    Scalar-callable. For |ζ| << rho0/|drho_a_dz| it reduces to :func:`parcel_displacement` with N^2 from Eq. (1.29).

    Validation (planned): V3 convergence to the linear solution as ζ0 → 0 (error ∝ ζ0^2). Label: pending.
    """
    return _sample_parcel(t, zeta0, lambda z: rho0 + drho_dz * z, lambda zeta: rho0 + drho_a_dz * zeta, g, w0,
                          zeta_max)


def parcel_acceleration_atmosphere(zeta, T0, dT_dz, cp=CP_AIR, g=G0):
    """Buoyant acceleration of a dry parcel displaced by ζ in an atmosphere with constant lapse rate.

    Book: §1.10 parcel argument for a perfect gas (our closed model of the buoyancy term): at equal pressures
    rho ∝ 1/T, so ``−g (rho_p − rho_e)/rho_p = g (T_p − T_e)/T_e`` with the parcel on the dry adiabat
    ``T_p = T0 − g ζ/C_p`` (Eq. 1.30 with α = 1/T) and the environment ``T_e = T0 + (dT/dz) ζ``.

    Parameters
    ----------
    zeta : float or array_like
        Upward displacement from the rest height [m].
    T0 : float
        Temperature at the rest height [K].
    dT_dz : float
        Environment lapse rate in Kundu's convention Γ = dT/dz [K/m] (negative when T falls with height;
        the meteorology lapse rate is −dT_dz).
    cp : float, optional
        [J/(kg K)].
    g : float, optional
        [m/s^2].

    Returns
    -------
    a : float or ndarray
        d^2ζ/dt^2 [m/s^2]; for small ζ it equals −N^2 ζ with N^2 = (g/T0)(dT/dz + g/C_p).

    Notes
    -----
    Model choice (ours; the book linearises at once): an isentropic parcel in hydrostatic surroundings actually cools at
    dT_p/dz = −(g/C_p)(T_p/T_e), which equals the constant −g/C_p only where T_p = T_e. Using the constant rate is
    exact to first order in ζ, so the linear limit (and N^2) are unaffected; the difference is O(ζ^2).

    Validation (planned): V1 derivative at ζ = 0 equals −:func:`brunt_vaisala_sq_from_lapse`; zero for dT/dz = −g/C_p.
    Label: pending.
    """
    zeta = np.asarray(zeta, dtype=float)
    T_p = T0 - g * zeta / cp  # parcel follows the dry adiabat, Eq. (1.30)
    T_e = T0 + dT_dz * zeta  # environment profile
    return as_scalar_if_0d(g * (T_p - T_e) / T_e)  # buoyancy at equal pressure: rho ∝ 1/T


def parcel_ode_atmosphere(t, zeta0, T0, dT_dz, cp=CP_AIR, g=G0, w0=0.0, zeta_max=None, rtol=1e-10, atol=1e-12):
    """Nonlinear dry-parcel displacement driven by :func:`parcel_acceleration_atmosphere` (E4 nonlinear path).

    Book: §1.10 parcel argument in a perfect-gas atmosphere (Fig. 1.8); ζ'' = g (T_p − T_e)/T_e.

    Parameters
    ----------
    t : float or array_like
        Output time(s) [s], >= 0.
    zeta0 : float
        Initial upward displacement [m].
    T0 : float
        Environment temperature at the rest height [K].
    dT_dz : float
        Environment lapse rate in Kundu's convention Γ = dT/dz [K/m] (meteorology: −dT_dz).
    cp, g : float, optional
        [J/(kg K)], [m/s^2].
    w0 : float, optional
        Initial velocity [m/s].
    zeta_max : float, optional
        Runaway cap [m] (terminal event for N^2 < 0); later times return NaN. Default ``1e3 * max(|zeta0|, 1e-3)``.
        The run also stops (NaN afterwards) if the parcel or environment temperature would fall below 10 % of T0,
        where the linear temperature profiles stop being physical.
    rtol, atol : float, optional
        ``solve_ivp`` tolerances.

    Returns
    -------
    zeta : float or ndarray
        [m].

    Notes
    -----
    Environment argument ``dT_dz`` is Kundu's Γ ≡ dT/dz (meteorology Γ = −dT_dz).

    Validation (planned): V3 matches :func:`parcel_displacement` for small ζ0; V1 dT/dz = −g/C_p keeps ζ = ζ0.
    Label: pending.
    """
    cap = 1e3 * max(abs(zeta0), 1e-3) if zeta_max is None else float(zeta_max)
    t_arr = np.atleast_1d(np.asarray(t, dtype=float))
    t_end = float(t_arr.max())
    if t_end <= 0.0:
        out = np.full_like(t_arr, float(zeta0))
        return as_scalar_if_0d(out[0]) if np.ndim(t) == 0 else out

    def rhs(_t, y):
        return [y[1], parcel_acceleration_atmosphere(y[0], T0, dT_dz, cp, g)]

    def hit_cap(_t, y):
        return abs(y[0]) - cap

    def too_cold(_t, y):  # guard: linear T profiles reach absolute zero far from the rest height
        return min(T0 - g * y[0] / cp, T0 + dT_dz * y[0]) - 0.1 * T0

    hit_cap.terminal = True
    too_cold.terminal = True
    sol = solve_ivp(rhs, (0.0, t_end), [zeta0, w0], t_eval=np.unique(np.concatenate(([0.0], t_arr))),
                    events=(hit_cap, too_cold), rtol=rtol, atol=atol)
    ts, zs = sol.t, sol.y[0]
    out = np.interp(t_arr, ts, zs) if ts.size > 1 else np.full_like(t_arr, float(zeta0))
    out = np.where(t_arr <= ts[-1] + 1e-12, out, np.nan)
    return as_scalar_if_0d(out[0]) if np.ndim(t) == 0 else out


# --------------------------------------------------------------------------------------------------------------------
# Lapse rates and potential temperature / density
# --------------------------------------------------------------------------------------------------------------------
def lapse_rate(T, z):
    """Lapse rate of a sampled profile in Kundu's convention Γ ≡ dT/dz (negative where T falls with height).

    Meteorology convention Γ ≡ −dT/dz: negate the result and flip every stability inequality
    (see :func:`lapse_rate_convention`, :func:`lapse_rate_stability`).

    Book: §1.10, ``Γ ≡ dT/dz``.

    Parameters
    ----------
    T : array_like
        Temperature samples [K].
    z : array_like
        Heights [m], strictly monotonic (uniform or not).

    Returns
    -------
    dT_dz : ndarray
        [K/m], second-order accurate everywhere (``np.gradient`` with ``edge_order=2``).

    Validation (planned): V1 linear profile returns its slope exactly; V3 order 2 on a smooth profile. Label: pending.
    """
    return np.gradient(np.asarray(T, dtype=float), np.asarray(z, dtype=float), edge_order=2)  # Γ ≡ dT/dz


def adiabatic_lapse_rate(T=None, cp=CP_AIR, alpha=None, g=G0):
    """Convention Γ ≡ dT/dz (Kundu's sign; NEGATIVE, ≈ −9.8 K/km in dry air; stable when dT/dz > Γ_a).

    Meteorology convention Γ ≡ −dT/dz: the same quantity is +9.8 K/km and stability flips to Γ < Γ_a
    (:func:`lapse_rate_convention`, :func:`lapse_rate_stability` write both forms).

    Adiabatic temperature gradient Γ_a = dT_a/dz of a parcel moving isentropically through a hydrostatic fluid.

    Book: §1.10, Eq. (1.30) ``dT_a/dz ≡ Γ_a = −g α T / C_p`` (derivation D19, Exercise 1.14). For a perfect gas
    α = 1/T (1.28) and Γ_a = −g/C_p. Most meteorology uses Γ ≡ −dT/dz instead (Γ_a ≈ +9.8 K/km, stable when
    Γ < Γ_a); to compare with that literature, negate this return value.

    Parameters
    ----------
    T : float or array_like, optional
        Temperature [K]; required when ``alpha`` is given (general fluid), ignored for a perfect gas.
    cp : float or array_like, optional
        Specific heat at constant pressure [J/(kg K)]; default ``CP_AIR``.
    alpha : float or array_like, optional
        Thermal expansion coefficient [1/K]; ``None`` means perfect gas (α = 1/T).
    g : float, optional
        Gravitational acceleration [m/s^2].

    Returns
    -------
    Gamma_a : float or ndarray
        dT/dz [K/m], **negative** (≈ −9.76e-3 K/m for dry air with C_p = 1004.7): the parcel cools as it rises.
        Under the meteorological convention Γ ≡ −dT/dz the same quantity is +9.76e-3 K/m.

    Notes
    -----
    Assumptions: isentropic parcel, hydrostatic surroundings, single-component fluid.

    Validation (planned): V1 perfect gas −g/C_p; general form with α = 1/T agrees; V5 |Γ_a| ≈ 9.8 K/km (AMS Glossary);
    V6 book value (private JSON). Label: pending.
    """
    if alpha is None:
        return as_scalar_if_0d(-g / np.asarray(cp, dtype=float))  # Eq. (1.30) with α = 1/T (1.28)
    if T is None:
        raise ValueError("adiabatic_lapse_rate needs T when alpha is given")
    return as_scalar_if_0d(-g * np.asarray(alpha, dtype=float) * np.asarray(T, dtype=float)
                           / np.asarray(cp, dtype=float))  # Eq. (1.30)


def parcel_temperature(T0, z, z0=0.0, cp=CP_AIR, g=G0):
    """Temperature of a dry parcel lifted adiabatically from (z0, T0): the dry adiabat through a point.

    Book: §1.10, Eq. (1.30) integrated for a perfect gas: ``T_a(z) = T0 + Γ_a (z − z0)`` with Γ_a = −g/C_p.

    Parameters
    ----------
    T0 : float
        Parcel temperature at z0 [K].
    z : float or array_like
        Height [m].
    z0 : float, optional
        Starting height [m].
    cp : float, optional
        [J/(kg K)].
    g : float, optional
        [m/s^2].

    Returns
    -------
    T_a : float or ndarray
        [K].

    Validation (planned): V1 slope equals :func:`adiabatic_lapse_rate`; θ constant along it (1.31). Label: pending.
    """
    return as_scalar_if_0d(T0 + adiabatic_lapse_rate(cp=cp, g=g) * (np.asarray(z, dtype=float) - z0))


def potential_temperature(T, p, p_ref=P_REF, gamma=GAMMA_AIR):
    """Potential temperature: the temperature a parcel reaches when brought adiabatically to the reference pressure.

    Book: §1.10, Eq. (1.31) ``T = θ (p/p_o)^((γ−1)/γ)``, solved for θ: ``θ = T (p_o/p)^((γ−1)/γ)``.

    Parameters
    ----------
    T : float or array_like
        Temperature [K].
    p : float or array_like
        Pressure [Pa], > 0.
    p_ref : float, optional
        Reference pressure p_o [Pa] (default 1.0e5; not the surface pressure p0).
    gamma : float, optional
        Ratio of specific heats [-].

    Returns
    -------
    theta : float or ndarray
        [K].

    Notes
    -----
    Assumptions: perfect gas, constant γ, reversible adiabatic process (D20).

    Validation (planned): V1 T = θ at p = p_ref; V7 θ constant along a dry adiabat built with
    ``atmosphere_from_temperature``. Label: pending.
    """
    p = np.asarray(p, dtype=float)
    require_positive("p", p)
    return as_scalar_if_0d(np.asarray(T, dtype=float) * (p_ref / p) ** ((gamma - 1.0) / gamma))  # Eq. (1.31) inverted


def temperature_from_potential(theta, p, p_ref=P_REF, gamma=GAMMA_AIR):
    """Temperature from potential temperature and pressure.

    Book: §1.10, Eq. (1.31) ``T(z) = θ(z) (p(z)/p_o)^((γ−1)/γ)`` (the printed form).

    Parameters
    ----------
    theta : float or array_like
        Potential temperature [K].
    p : float or array_like
        Pressure [Pa].
    p_ref : float, optional
        Reference pressure p_o [Pa].
    gamma : float, optional
        [-].

    Returns
    -------
    T : float or ndarray
        [K].

    Validation (planned): V1 round trip with :func:`potential_temperature`. Label: pending.
    """
    p = np.asarray(p, dtype=float)
    require_positive("p", p)
    return as_scalar_if_0d(np.asarray(theta, dtype=float) * (p / p_ref) ** ((gamma - 1.0) / gamma))  # Eq. (1.31)


def potential_temperature_gradient(T, dT_dz, theta=None, cp=CP_AIR, g=G0, p=None, p_ref=P_REF,
                                   gamma=GAMMA_AIR):
    """dθ/dz from the lapse rate: (T/θ) dθ/dz = Γ − Γ_a with Kundu's Γ ≡ dT/dz (stable when Γ > Γ_a, dθ/dz > 0).

    Meteorology convention Γ ≡ −dT/dz: the same equation reads (T/θ) dθ/dz = Γ_a − Γ, stable when Γ < Γ_a. The
    ``dT_dz`` argument is always Kundu's sign (convert a meteorology Γ with :func:`lapse_rate_convention`).

    Book: §1.10, Eq. (1.32) ``(T/θ) dθ/dz = dT/dz + g/C_p = d(T − T_a)/dz = Γ − Γ_a``, so
    ``dθ/dz = (θ/T)(dT/dz + g/C_p)``.

    Parameters
    ----------
    T : float or array_like
        Temperature [K].
    dT_dz : float or array_like
        Environment lapse rate in Kundu's convention Γ = dT/dz [K/m] (negative when T falls with height;
        the meteorology lapse rate is −dT_dz).
    theta : float or array_like, optional
        Potential temperature [K]; if omitted it is computed from ``p``.
    cp : float, optional
        [J/(kg K)].
    g : float, optional
        [m/s^2].
    p : float or array_like, optional
        Pressure [Pa] (needed when ``theta`` is None).
    p_ref : float, optional
        Reference pressure [Pa].
    gamma : float, optional
        [-]; must be consistent with cp and R ((γ−1)/γ = R/C_p) for (1.32) to hold exactly.

    Returns
    -------
    dtheta_dz : float or ndarray
        [K/m]; > 0 stable, 0 neutral, < 0 unstable.

    Notes
    -----
    Assumptions: perfect gas, hydrostatic, constant C_p (D21).

    Validation (planned): V2 sympy residual of (1.32) = 0; V1 matches centred differences of θ(z) along a synthetic
    profile. Label: pending.
    """
    T = np.asarray(T, dtype=float)
    if theta is None:
        if p is None:
            raise ValueError("give theta or p")
        theta = potential_temperature(T, p, p_ref, gamma)
    theta = np.asarray(theta, dtype=float)
    return as_scalar_if_0d(theta / T * (np.asarray(dT_dz, dtype=float) + g / cp))  # Eq. (1.32)


# --------------------------------------------------------------------------------------------------------------------
# The two lapse-rate conventions (user decision: compute with Kundu's, always show the meteorology form too)
# --------------------------------------------------------------------------------------------------------------------
#: The two accepted values of ``convention`` in :func:`lapse_rate_convention` and :func:`lapse_rate_stability`.
LAPSE_RATE_CONVENTIONS: tuple[str, ...] = ("kundu", "meteorology")


def _convention(convention: str) -> str:
    key = str(convention).strip().lower()
    if key not in LAPSE_RATE_CONVENTIONS:
        raise ValueError(f"unknown lapse-rate convention {convention!r}; use 'kundu' (Γ ≡ dT/dz) or "
                         "'meteorology' (Γ ≡ −dT/dz)")
    return key


def lapse_rate_convention(dT_dz, convention: str = "kundu"):
    """Convention switch: Γ_met = −Γ_Kundu; convert the number AND flip the stability inequality.

    Returns the lapse rate Γ expressed in the requested convention: ``"kundu"`` gives ``dT_dz`` unchanged
    (Γ ≡ dT/dz), ``"meteorology"`` gives ``−dT_dz`` (Γ ≡ −dT/dz).

    Book: §1.10 defines ``Γ ≡ dT/dz`` (Kundu's convention, used by every fluidpy computation). The standard
    meteorology convention (AMS Glossary, most atmospheric-science texts and the lapse-rate-feedback literature)
    defines ``Γ ≡ −dT/dz``. The physics is identical; the sign of every lapse rate and the direction of every
    stability inequality flip (Kundu: stable when Γ > Γ_a; meteorology: stable when Γ < Γ_a).

    Parameters
    ----------
    dT_dz : float or array_like
        Vertical temperature gradient in Kundu's sign [K/m] (negative when T falls with height). Works for an
        environment profile and for the adiabatic value from :func:`adiabatic_lapse_rate` alike.
    convention : {"kundu", "meteorology"}, optional
        Target convention (case-insensitive); any other string raises ValueError.

    Returns
    -------
    Gamma : float or ndarray
        Γ in the requested convention [K/m]: ``dT_dz`` for "kundu", ``−dT_dz`` for "meteorology".

    Notes
    -----
    The map is its own inverse: ``lapse_rate_convention(Gamma_met, "meteorology")`` turns a meteorology lapse rate
    back into the ``dT_dz`` that fluidpy functions expect. Example: the dry adiabat is −9.76e-3 K/m (Kundu) and
    +9.76e-3 K/m = +g/C_p (meteorology).

    Validation (planned): V1 Kundu identity; meteorology Γ_a = +g/C_p; involution (applying twice returns the input).
    Label: pending.
    """
    x = np.asarray(dT_dz, dtype=float)
    Gamma = x if _convention(convention) == "kundu" else -x  # Γ_Kundu ≡ dT/dz ; Γ_met ≡ −dT/dz
    return as_scalar_if_0d(Gamma)


class LapseStability(NamedTuple):
    """Result of :func:`lapse_rate_stability`; the first two fields are ``(verdict, inequality_text)``.

    Fields: ``verdict`` ("stable" / "neutral" / "unstable"); ``text`` (criterion with numbers); ``code`` (+1 stable,
    0 neutral, −1 unstable); ``Gamma`` and ``Gamma_a`` [K/m] in the requested convention; ``margin`` [K/m], positive
    when stable and identical in both conventions (= (T/θ) dθ/dz of Eq. 1.32 for a perfect gas).
    """

    verdict: str
    text: str
    code: int
    Gamma: float
    Gamma_a: float
    margin: float


def _fmt_number(x: float, digits: int, ascii_minus: bool) -> str:
    s = f"{x:.{digits}f}"
    if float(s) == 0.0:
        s = s.lstrip("-")  # never print "−0.0"
    return s if ascii_minus else s.replace("-", "−")


def lapse_rate_stability(dT_dz, Gamma_a=None, convention: str = "kundu", tol: float = 1e-9, per_km: bool = True,
                         prefix: bool = True, digits: int | None = None, ascii_minus: bool = False,
                         cp=CP_AIR, g=G0) -> LapseStability:
    """Stability verdict from the lapse rate plus the criterion written out with numbers, in either convention.

    Book: §1.10, Eqs. (1.30) and (1.32): ``(T/θ) dθ/dz = dT/dz + g/C_p = Γ − Γ_a`` with Kundu's ``Γ ≡ dT/dz``, so a
    layer is stable when ``dT/dz > Γ_a`` (θ increases upward), neutral when equal, unstable when ``dT/dz < Γ_a``.
    Meteorology convention ``Γ ≡ −dT/dz``: both numbers change sign and the inequality flips, stable when ``Γ < Γ_a``.

    Parameters
    ----------
    dT_dz : float
        Environment temperature gradient in **Kundu's sign** [K/m] (standard troposphere −6.5e-3). Convert a
        meteorology lapse rate first with ``lapse_rate_convention(Gamma_met, "meteorology")``.
    Gamma_a : float, optional
        Adiabatic temperature gradient dT_a/dz in **Kundu's sign** [K/m] (negative for air); default
        ``adiabatic_lapse_rate(cp=cp, g=g)`` = −g/C_p (dry perfect gas).
    convention : {"kundu", "meteorology"}, optional
        Convention of ``text``, ``Gamma`` and ``Gamma_a`` in the result (``verdict``, ``code``, ``margin`` do not
        depend on it).
    tol : float, optional
        Neutral band |margin| <= tol [K/m] (default 1e-9 K/m = 1e-6 K/km, a round-off guard).
    per_km : bool, optional
        Numbers in the text in K/km (default) or K/m.
    prefix : bool, optional
        Start the text with the criterion itself: ``"stable ⇔ dT/dz > Γa: "`` (Kundu) or ``"stable ⇔ Γ < Γa: "``
        (meteorology). False gives the bare evaluated inequality, e.g. ``"−6.5 > −9.8 K/km"``.
    digits : int, optional
        Decimals in the text; default 1 in K/km and 5 in K/m. Raised automatically (up to 6) when the two numbers
        would print alike although the layer is not neutral (e.g. −9.8 vs Γa = −9.76 K/km prints "−9.80 < −9.76").
    ascii_minus : bool, optional
        Use "-" instead of the typographic minus "−" (U+2212, default).
    cp : float, optional
        Specific heat at constant pressure [J/(kg K)] for the default Γ_a.
    g : float, optional
        Gravitational acceleration [m/s^2] for the default Γ_a.

    Returns
    -------
    LapseStability
        NamedTuple ``(verdict, text, code, Gamma, Gamma_a, margin)``. Kundu: ``Gamma = dT_dz``, ``Gamma_a`` negative,
        ``margin = Gamma − Gamma_a``, text ``"stable ⇔ dT/dz > Γa: −6.5 > −9.8 K/km"``. Meteorology:
        ``Gamma = −dT_dz``, ``Gamma_a`` positive, ``margin = Gamma_a − Gamma``, text
        ``"stable ⇔ Γ < Γa: 6.5 < 9.8 K/km"``. The relation printed is the true one (">", "<" or "="), so an unstable
        layer reads e.g. ``"stable ⇔ dT/dz > Γa: −12.0 < −9.8 K/km"``.

    Notes
    -----
    Assumptions: dry (unsaturated) parcel, hydrostatic surroundings; with the default Γ_a also a perfect gas with
    constant C_p. ``code`` equals the sign of N^2 = (g/T)(dT/dz − Γ_a) (D36) for T > 0. Scalar input only (the text is
    one sentence); use :func:`lapse_rate_convention` for arrays.

    Validation (planned): V1 texts for −6.5 K/km in both conventions; V1 meteorology Γ_a = +g/C_p; V4 ``code``
    identical in both conventions over −15…+10 K/km and equal to the sign of :func:`brunt_vaisala_sq_from_lapse`.
    Label: pending.
    """
    if np.ndim(dT_dz) != 0 or (Gamma_a is not None and np.ndim(Gamma_a) != 0):
        raise ValueError("lapse_rate_stability is scalar-only; use lapse_rate_convention for arrays")
    conv = _convention(convention)
    x = float(dT_dz)
    Ga_k = float(adiabatic_lapse_rate(cp=cp, g=g)) if Gamma_a is None else float(Gamma_a)
    margin = x - Ga_k  # Eq. (1.32): (T/θ) dθ/dz = Γ − Γ_a in Kundu's convention; > 0 stable
    code = 1 if margin > tol else (-1 if margin < -tol else 0)
    verdict = {1: "stable", 0: "neutral", -1: "unstable"}[code]
    if conv == "kundu":
        Gamma, Gamma_a_c = x, Ga_k  # stable ⇔ dT/dz > Γ_a
        op = {1: ">", 0: "=", -1: "<"}[code]
        crit = "stable ⇔ dT/dz > Γa: "
    else:
        Gamma, Gamma_a_c = -x, -Ga_k  # Γ_met = −dT/dz, Γ_a,met = −Γ_a: stable ⇔ Γ < Γ_a (inequality flipped)
        op = {1: "<", 0: "=", -1: ">"}[code]
        crit = "stable ⇔ Γ < Γa: "
    scale, unit = (1.0e3, "K/km") if per_km else (1.0, "K/m")
    d = (1 if per_km else 5) if digits is None else int(digits)
    a, b = Gamma * scale, Gamma_a_c * scale
    while code != 0 and d < 6 and _fmt_number(a, d, True) == _fmt_number(b, d, True):
        d += 1
    body = f"{_fmt_number(a, d, ascii_minus)} {op} {_fmt_number(b, d, ascii_minus)} {unit}"
    return LapseStability(verdict, (crit if prefix else "") + body, code, Gamma, Gamma_a_c, margin)


def potential_density(rho, p, p_ref=P_REF, gamma=GAMMA_AIR):
    """Potential density: the density a parcel reaches when brought isentropically to the reference pressure.

    Book: §1.10, Eq. (1.33) ``rho = rho_θ (p/p_o)^(1/γ)``, solved for rho_θ: ``rho_θ = rho (p_o/p)^(1/γ)``.

    Parameters
    ----------
    rho : float or array_like
        Density [kg/m^3].
    p : float or array_like
        Pressure [Pa].
    p_ref : float, optional
        Reference pressure p_o [Pa].
    gamma : float, optional
        [-].

    Returns
    -------
    rho_theta : float or ndarray
        [kg/m^3].

    Notes
    -----
    Assumptions: perfect gas, constant γ. θ rho_θ = p_o/R everywhere (D23), hence (1.34).

    Validation (planned): V4 invariant θ rho_θ = p_ref/R; V1 (1.34) by finite differences. Label: pending.
    """
    p = np.asarray(p, dtype=float)
    require_positive("p", p)
    return as_scalar_if_0d(np.asarray(rho, dtype=float) * (p_ref / p) ** (1.0 / gamma))  # Eq. (1.33) inverted


def isentropic_density_gradient(rho, c, g=G0):
    """Density gradient followed by a parcel displaced isentropically (and at constant salinity) in a static fluid.

    Book: §1.10, before Eq. (1.35): ``drho_a/dz = (∂rho_a/∂p)_{s,S} dp_a/dz = −rho_a g/c^2 ≅ −rho g/c^2``.

    Parameters
    ----------
    rho : float or array_like
        Density [kg/m^3].
    c : float or array_like
        Speed of sound [m/s].
    g : float, optional
        [m/s^2].

    Returns
    -------
    drho_a_dz : float or ndarray
        [kg/m^4] (negative).

    Notes
    -----
    Assumptions: (∂rho/∂p)_{s,S} = 1/c^2 (inverse of (1.19)); the parcel's pressure follows hydrostatics; rho_a ≅ rho.

    Validation (planned): V1 −rho g/c^2; V2 for a perfect gas equals the isentropic parcel gradient from (1.26).
    Label: pending.
    """
    c = np.asarray(c, dtype=float)
    require_positive("c", c)
    return as_scalar_if_0d(-np.asarray(rho, dtype=float) * g / c ** 2)  # §1.10: drho_a/dz = −rho g/c^2


def ocean_potential_density_gradient(drho_dz, rho, c, g=G0):
    """Ocean static-stability indicator: the compressibility-corrected density gradient.

    Book: §1.10, Eq. (1.35) ``drho_θ/dz = drho/dz − drho_a/dz ≅ drho/dz + rho g/c^2`` (negative: stable; zero: neutral;
    positive: unstable).

    Parameters
    ----------
    drho_dz : float or array_like
        Observed (in-situ) density gradient [kg/m^4], z upward.
    rho : float or array_like
        Density [kg/m^3].
    c : float or array_like
        Speed of sound [m/s].
    g : float, optional
        [m/s^2].

    Returns
    -------
    drho_theta_dz : float or ndarray
        [kg/m^4]. Has the **same sign** as the potential-density gradient; the book's equality holds only up to a positive
        factor ≈ 1 near the reference pressure (analysis §9, typo 7). Its sign is the opposite of N^2.

    Notes
    -----
    Assumptions: seawater parcel at fixed salinity, rho_a ≅ rho, small fractional density changes.

    Validation (planned): V1 sign equals −sign(N^2) from :func:`brunt_vaisala_sq`; V7 c → ∞ gives drho/dz. Label: pending.
    """
    c = np.asarray(c, dtype=float)
    require_positive("c", c)
    return as_scalar_if_0d(np.asarray(drho_dz, dtype=float)
                           + np.asarray(rho, dtype=float) * g / c ** 2)  # Eq. (1.35)


__all__ = [
    "brunt_vaisala_sq", "brunt_vaisala_sq_from_theta", "brunt_vaisala_sq_from_lapse", "classify_stability",
    "stability_timescale", "parcel_displacement", "parcel_ode", "parcel_ode_from_gradients",
    "parcel_acceleration_atmosphere", "parcel_ode_atmosphere",
    "lapse_rate", "adiabatic_lapse_rate", "LAPSE_RATE_CONVENTIONS", "LapseStability", "lapse_rate_convention", "lapse_rate_stability",
    "parcel_temperature", "potential_temperature", "temperature_from_potential",
    "potential_temperature_gradient", "potential_density", "isentropic_density_gradient",
    "ocean_potential_density_gradient",
]
