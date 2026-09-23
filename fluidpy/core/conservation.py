"""Integral conservation laws of mass, momentum, energy and angular momentum on any control volume.

Book: Kundu, Cohen & Dowling 5e, Ch. 4 — mass (4.1)–(4.5); momentum (4.13)–(4.17); energy (4.46)–(4.48); angular
momentum (4.64)–(4.65) (rendered pages chapters/pages/ch04/p124–p131, p148–p153).

A control volume is any ``core.transport.ControlVolume`` (fixed or moving box, growing sphere/cylinder/cone, moving
ellipse): its ``volume_nodes`` / ``surface_nodes`` give quadrature points, outward normals n, area elements and the
surface velocity b (Fig. 3.18). Every budget returns each term separately and the **residual**, which vanishes when the
fields obey the law:

* ``storage`` d/dt ∫_{V*} (·) dV — the central difference in t of the volume integral over the *moving* CV (the
  definition (3.31), independent of the RTT). With ``material=True`` (the (4.1)/(4.13)/(4.46) check) the surface moves
  with the fluid, b = u, at the instant: storage is then the RTT (3.35) value ∫∂(·)/∂t dV + ∮(·)u·n dA of the material
  volume coincident with the CV (the coincident-CV step of (4.4), (4.16));
* ``outflux`` ∮ (·)(u − b)·n dA — the flux relative to the moving surface (u and b in the same frame).

Sign conventions: n points out of the CV; forces are those exerted **on the fluid inside** (a drag on a body is −F_D on
the fluid, Example 4.1); z is up and g = (0, 0, −g), g = 9.81 m/s² by default (the book's value). A 2-D control
volume (e.g. ``MovingEllipse2D`` in the x–y plane) keeps only the first two components of g, so the default vertical
gravity is **dropped** (a horizontal plane); pass ``g=(0.0, -9.81)`` to put gravity along −y in a vertical plane.
Book typos handled: (4.15) ends with a spurious "= 0" (it is an identity — (3.35) rearranged); (4.51) writes dA for the
volume integrals (dV) — analysis §9 items 1–2.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from ._stencil import DEFAULT_HT, as_field, ev, gvec
from .thermo import G_BOOK as G0
from .transport import ControlVolume

__all__ = ["MassBudget", "mass_budget", "MomentumBudget", "momentum_budget", "EnergyBudget", "energy_budget",
           "interval_mass_budget", "angular_momentum_flux", "angular_momentum_flux_surface",
           "AngularMomentumBudget", "angular_momentum_budget"]

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731


def _vol_int(F: Callable, cv: ControlVolume, t: float, n: int, lead: tuple = ()) -> np.ndarray:
    X, W = cv.volume_nodes(t, n)
    return np.sum(ev(F, X, t, lead) * W, axis=-1)


def _dt_field(F: Callable, dFdt: Callable | None, X, t, lead):
    if dFdt is not None:
        return ev(dFdt, X, t, lead)
    ht = DEFAULT_HT
    return (ev(F, X, t + ht, lead) - ev(F, X, t - ht, lead)) / (2.0 * ht)


def _storage(F: Callable, dFdt: Callable | None, cv: ControlVolume, t: float, n: int, dt: float,
             material: bool, u: Callable, lead: tuple = ()) -> np.ndarray:
    """d/dt ∫ F dV over the moving CV (central difference, (3.31)) or over the coincident material volume (RTT)."""
    if not material:
        return (_vol_int(F, cv, t + dt, n, lead) - _vol_int(F, cv, t - dt, n, lead)) / (2.0 * dt)  # Eq. (3.31)
    X, W = cv.volume_nodes(t, n)
    vol = np.sum(_dt_field(F, dFdt, X, t, lead) * W, axis=-1)
    Xs, N, dA, _ = cv.surface_nodes(t, n)
    un = np.einsum("ik,ik->k", ev(u, Xs, t, (Xs.shape[0],)), N)
    return vol + np.sum(ev(F, Xs, t, lead) * un * dA, axis=-1)  # Eq. (3.35) with b = u


def _rel_flux(u: Callable, material: bool, cv: ControlVolume, t: float, n: int):
    """Surface nodes and (u − b)·n there (b = the CV's surface velocity, or u for a material surface)."""
    Xs, N, dA, bcv = cv.surface_nodes(t, n)
    U = ev(u, Xs, t, (Xs.shape[0],))
    B = U if material else bcv
    return Xs, N, dA, U, np.einsum("ik,ik->k", U - B, N)


def _traction(Xs, N, t, traction, tau, d):
    if traction is not None:
        return ev(lambda X, T: traction(X, N, T), Xs, t, (d,))
    if tau is not None:
        return np.einsum("ijk,ik->jk", ev(tau, Xs, t, (d, d)), N)  # Eq. (2.15): f_j = n_i τ_ij
    return np.zeros((d, Xs.shape[1]))


# ======================================================================================================================
# mass (4.1)–(4.5)
# ======================================================================================================================
@dataclass(frozen=True)
class MassBudget:
    """Terms of Eq. (4.5) [kg/s]: ``storage`` d/dt∫_{V*}ρ dV + ``outflux`` ∮ρ(u − b)·n dA = ``residual`` (0 when mass is
    conserved); ``local`` = ∫_{V*}∂ρ/∂t dV (NaN unless ∂ρ/∂t is given); ``mass`` = ∫_{V*}ρ dV [kg].
    Book: §4.2, Eq. (4.5). Label: analytic.
    """

    storage: float
    outflux: float
    local: float
    residual: float
    mass: float


def mass_budget(rho, u: Callable, cv: ControlVolume, t: float = 0.0, drho_dt: Callable | None = None,
                material: bool = False, n: int = 24, dt: float = 1e-4) -> MassBudget:
    """Mass budget of an arbitrarily moving control volume, Eq. (4.5).

    Book: §4.2, Eq. (4.5): d/dt ∫_{V*(t)} ρ dV + ∫_{A*(t)} ρ(u − b)·n dA = 0; b = u gives the material volume (4.1),
    b = 0 a fixed CV ((4.2)/(4.6)). Derived from (4.1)–(4.4) by the coincident-CV step (our D01).

    Parameters
    ----------
    rho : density field ρ(x, t) [kg/m³] (callable or constant);  u : velocity field [m/s]
    cv : a ``core.transport.ControlVolume``;  t : time [s]
    drho_dt : ∂ρ/∂t field [kg/(m³ s)] (fills ``local``; else ``local`` is NaN)
    material : True → the surface moves with the fluid (b = u): outflux 0 and storage = the rate of the coincident
        material volume ((4.1) check); False → the CV's own surface velocity b
    n : quadrature nodes per direction;  dt : time step of the storage central difference [s]

    Returns
    -------
    :class:`MassBudget` (plain floats).

    Assumptions: continuum, smooth ρ and u, piecewise smooth CV surface.
    Validation: V4 fixed box in the expanding flow: storage = −outflux (residual O(dt²)); material: outflux 0 and
    storage 0; V3 a moving, growing sphere: residual → 0 at order 2 in dt. Label: conserved, converged.
    """
    rho_f = as_field(rho)
    R = lambda X, T: ev(rho_f, X, T)  # noqa: E731
    storage = float(_storage(R, drho_dt, cv, t, n, dt, material, u))
    Xs, N, dA, U, rel = _rel_flux(u, material, cv, t, n)
    outflux = float(np.sum(R(Xs, t) * rel * dA))  # ∮ ρ(u − b)·n dA
    local = float(_vol_int(drho_dt, cv, t, n)) if drho_dt is not None else float("nan")
    mass = float(_vol_int(R, cv, t, n))
    return MassBudget(storage, outflux, local, storage + outflux, mass)  # Eq. (4.5)


def interval_mass_budget(x0: float, x1: float, t: float = 0.0, dx0dt: float = 0.0, dx1dt: float = 0.0,
                         flow="expanding", dt: float = 1e-5, **p) -> dict:
    """One-dimensional mass budget (4.5) per unit cross-section for the interval [x0, x1] whose ends move at dx0dt,
    dx1dt (fixed CV: 0, 0; material interval: the local fluid speeds).

    Book: §4.2, Eq. (4.5) read in 1-D (with ẋ = u at both ends it is (4.1)); the C01 tiny example and E1's mass mode.

    Parameters
    ----------
    x0, x1 : end positions at t [m];  t [s];  dx0dt, dx1dt : end speeds [m/s]
    flow : "expanding" (keywords a [1/s] = 1.0, rho0 [kg/m³] = 1.0: u = ax/(1 + at), ρ = ρ₀/(1 + at)) or a pair
        (rho_fn, u_fn) of callables of (x, t) (floats or arrays of x)
    dt : central-difference step for the storage of a callable flow [s]

    Returns
    -------
    dict(storage = d/dt∫_{x0(t)}^{x1(t)}ρ dx, flux_left = ρ(u − ẋ₀) at x0 (inflow when > 0), flux_right = ρ(u − ẋ₁) at
    x1, net_outflux = flux_right − flux_left, local = ∫∂ρ/∂t dx, residual = storage + net_outflux) [kg/(m² s)].

    Validation: V1 expanding flow (a = 1, ρ₀ = 1) on [1, 2] at t = 0: storage −1, flux_left 1, flux_right 2, residual 0;
    material ends (ẋ = u): all fluxes and the storage 0. Label: analytic.
    """
    from scipy.integrate import quad
    if isinstance(flow, str):
        if flow != "expanding":
            raise ValueError("flow must be 'expanding' or a pair (rho_fn, u_fn)")
        a, rho0 = float(p.get("a", 1.0)), float(p.get("rho0", 1.0))
        s = 1.0 + a * t
        rho_fn = lambda x, tt: rho0 / (1.0 + a * tt) + 0.0 * np.asarray(x, dtype=float)  # noqa: E731
        u_fn = lambda x, tt: a * np.asarray(x, dtype=float) / (1.0 + a * tt)  # noqa: E731
        local = -a * rho0 / s ** 2 * (x1 - x0)  # ∫∂ρ/∂t dx (ρ uniform in x)
        storage = local + rho0 / s * (dx1dt - dx0dt)  # Leibniz (3.30), exact
    else:
        rho_fn, u_fn = flow
        I = lambda tt: quad(lambda x: float(rho_fn(x, tt)), x0 + dx0dt * (tt - t), x1 + dx1dt * (tt - t),  # noqa: E731
                            epsabs=0.0, epsrel=1e-12)[0]
        storage = (I(t + dt) - I(t - dt)) / (2.0 * dt)
        local = quad(lambda x: (float(rho_fn(x, t + dt)) - float(rho_fn(x, t - dt))) / (2.0 * dt), x0, x1,
                     epsabs=0.0, epsrel=1e-10)[0]
    fl = float(rho_fn(x0, t)) * (float(u_fn(x0, t)) - dx0dt)  # ρ(u − b) at the left end (inflow if > 0)
    fr = float(rho_fn(x1, t)) * (float(u_fn(x1, t)) - dx1dt)  # ρ(u − b) at the right end (outflow if > 0)
    net = fr - fl
    return {"storage": float(storage), "flux_left": fl, "flux_right": fr, "net_outflux": net, "local": float(local),
            "residual": float(storage + net)}  # Eq. (4.5), 1-D


# ======================================================================================================================
# momentum (4.13)–(4.17)
# ======================================================================================================================
@dataclass(frozen=True)
class MomentumBudget:
    """Terms of Eq. (4.17) [N] (arrays (d,)): ``storage`` d/dt∫ρu dV + ``outflux`` ∮ρu(u − b)·n dA = ``body`` ∫ρg dV +
    ``surface`` ∮f dA; ``residual`` = storage + outflux − body − surface. Book: §4.4, Eq. (4.17). Label: analytic.
    """

    storage: np.ndarray
    outflux: np.ndarray
    body: np.ndarray
    surface: np.ndarray
    residual: np.ndarray


def momentum_budget(rho, u: Callable, cv: ControlVolume, t: float = 0.0, g=(0.0, 0.0, -G0),
                    traction: Callable | None = None, n: int = 24, dt: float = 1e-4, tau: Callable | None = None,
                    material: bool = False, d_rho_u_dt: Callable | None = None) -> MomentumBudget:
    """Linear-momentum budget of an arbitrarily moving control volume, Eq. (4.17).

    Book: §4.4, Eq. (4.17): d/dt∫_{V*}ρu dV + ∫_{A*}ρu(u − b)·n dA = ∫_{V*}ρg dV + ∫_{A*}f(n, x, t) dA; b = u gives
    Newton's second law for a material volume (4.13)–(4.14); the coincident-CV equalities (4.16a–d) are our D05.

    Parameters
    ----------
    rho [kg/m³]; u [m/s]; cv : ControlVolume; t [s]; g : body force per mass [m/s²] (first d components used)
    traction : callable f(x, n_hat, t) → (d, K) force per area on the fluid from outside [Pa]; None → no surface force
    n, dt : quadrature nodes, storage time step [s]
    tau : alternatively a stress field τ(x, t) [Pa] (then f_j = n_iτ_ij, (2.15))
    material : b = u (the (4.13) check);  d_rho_u_dt : ∂(ρu)/∂t field for the material storage

    Returns
    -------
    :class:`MomentumBudget` — arrays (d,) [N].

    Validation: V1 uniform flow through a fixed box: outflux 0, residual 0; fluid at rest under gravity with
    τ = −p_s δ: surface = −body (Archimedes, ``core.statics.net_pressure_force_on_box``); V4 material b = u reproduces
    (4.13). Label: analytic, conserved.
    """
    rho_f = as_field(rho)
    X0, _ = cv.volume_nodes(t, n)
    d = X0.shape[0]
    rho_u = lambda X, T: ev(rho_f, X, T) * ev(u, X, T, (d,))  # noqa: E731
    storage = _storage(rho_u, d_rho_u_dt, cv, t, n, dt, material, u, (d,))
    Xs, N, dA, U, rel = _rel_flux(u, material, cv, t, n)
    outflux = np.sum(ev(rho_f, Xs, t) * U * rel * dA, axis=-1)  # ∮ ρu (u − b)·n dA
    Xv, W = cv.volume_nodes(t, n)
    body = np.sum(ev(rho_f, Xv, t) * gvec(g, d, (Xv.shape[1],)) * W, axis=-1)  # ∫ ρg dV
    surface = np.sum(_traction(Xs, N, t, traction, tau, d) * dA, axis=-1)  # ∮ f dA
    return MomentumBudget(storage, outflux, body, surface, storage + outflux - body - surface)  # Eq. (4.17)


# ======================================================================================================================
# energy (4.46)–(4.48)
# ======================================================================================================================
@dataclass(frozen=True)
class EnergyBudget:
    """Terms of Eq. (4.48) [W]: ``storage`` d/dt∫ρ(e + ½|u|²)dV + ``outflux`` ∮ρ(e + ½|u|²)(u − b)·n dA = ``body_work``
    ∫ρg·u dV + ``surface_work`` ∮f·u dA − ``heat_out`` ∮q·n dA; ``residual``. Book: §4.8, Eq. (4.48). Label: analytic.
    """

    storage: float
    outflux: float
    body_work: float
    surface_work: float
    heat_out: float
    residual: float


def energy_budget(rho, u: Callable, e, cv: ControlVolume, t: float = 0.0, g=(0.0, 0.0, -G0),
                  traction: Callable | None = None, q: Callable | None = None, n: int = 24, dt: float = 1e-4,
                  tau: Callable | None = None, material: bool = False) -> EnergyBudget:
    """Total-energy budget (internal + kinetic) of an arbitrarily moving control volume, Eq. (4.48).

    Book: §4.8, Eqs. (4.46)–(4.48): d/dt∫_{V*}ρ(e + ½|u|²)dV + ∫_{A*}(ρe + ½ρ|u|²)(u − b)·n dA = ∫_{V*}ρg·u dV +
    ∫_{A*}f·u dA − ∫_{A*}q·n dA — energy changes by the work of body and surface forces minus the heat that leaves
    (q·n > 0 out).

    Parameters
    ----------
    rho [kg/m³]; u [m/s]; e : internal energy [J/kg] (callable or constant); cv; t [s]; g [m/s²];
    traction : as :func:`momentum_budget`; q : heat-flux field [W/m²] (None → 0); n, dt; tau; material

    Returns :class:`EnergyBudget` (floats) [W].
    Validation: V1 adiabatic inviscid uniform flow: all terms 0; steady Couette with viscous heating (fixed CV spanning
    the gap): shear work in at the moving wall = heat out through the walls; V4 material CV reproduces (4.46).
    Label: analytic, conserved.
    """
    rho_f, e_f = as_field(rho), as_field(e)
    X0, _ = cv.volume_nodes(t, n)
    d = X0.shape[0]
    E = lambda X, T: ev(rho_f, X, T) * (ev(e_f, X, T) + 0.5 * np.sum(ev(u, X, T, (d,)) ** 2, axis=0))  # noqa: E731
    storage = float(_storage(E, None, cv, t, n, dt, material, u))
    Xs, N, dA, U, rel = _rel_flux(u, material, cv, t, n)
    outflux = float(np.sum(E(Xs, t) * rel * dA))
    Xv, W = cv.volume_nodes(t, n)
    body = float(np.sum(ev(rho_f, Xv, t) * np.sum(gvec(g, d, (Xv.shape[1],)) * ev(u, Xv, t, (d,)), axis=0) * W))
    f = _traction(Xs, N, t, traction, tau, d)
    surf = float(np.sum(np.sum(f * U, axis=0) * dA))  # ∮ f·u dA
    heat = 0.0 if q is None else float(np.sum(np.einsum("ik,ik->k", ev(q, Xs, t, (d,)), N) * dA))  # ∮ q·n dA
    return EnergyBudget(storage, outflux, body, surf, heat, storage + outflux - body - surf + heat)  # Eq. (4.48)


# ======================================================================================================================
# angular momentum (4.64)–(4.65)
# ======================================================================================================================
def angular_momentum_flux_surface(rho, u: Callable, X, N, dA, t: float = 0.0, origin=(0.0, 0.0, 0.0)) -> np.ndarray:
    """∫ (r × ρu)(u·n) dA over a sampled surface (points X (3, K), outward normals N (3, K), areas dA (K,)) [N m].

    Book: §4.9, Eq. (4.65) (the flux term), Example 4.6. r is measured from ``origin`` (a point on the chosen axis).
    Label: analytic.
    """
    X_, N_ = _F(X), _F(N)
    rho_f = as_field(rho)
    U = ev(u, X_, t, (3,))
    r = X_ - _F(origin).reshape(3, 1)
    return np.sum(np.cross(r, ev(rho_f, X_, t) * U, axis=0) * np.einsum("ik,ik->k", U, N_) * _F(dA), axis=-1)


def angular_momentum_flux(rho, u: Callable, cv: ControlVolume, t: float = 0.0, origin=(0.0, 0.0, 0.0),
                          n: int = 24) -> np.ndarray:
    """Net outflow of angular momentum ∮_{A_o}(r × ρu)(u·n) dA through a stationary CV, Eq. (4.65) [N m], (3,).

    Book: §4.9, Eq. (4.65) (stationary CV, b = 0). 3-D control volumes only. Label: analytic.
    Validation: V1 two small jet discs of a sprinkler give 2aρAU² cos α (Example 4.6's closed form).
    """
    Xs, N, dA, _ = cv.surface_nodes(t, n)
    return angular_momentum_flux_surface(rho, u, Xs, N, dA, t, origin)


@dataclass(frozen=True)
class AngularMomentumBudget:
    """Terms of Eq. (4.65) [N m] (3-vectors): ``storage`` d/dt∫(r × ρu)dV + ``outflux`` ∮(r × ρu)(u·n)dA =
    ``body_torque`` ∫(r × ρg)dV + ``surface_torque`` ∮(r × f)dA; ``residual``. Book: §4.9, Eq. (4.65). Label: analytic.
    """

    storage: np.ndarray
    outflux: np.ndarray
    body_torque: np.ndarray
    surface_torque: np.ndarray
    residual: np.ndarray


def angular_momentum_budget(rho, u: Callable, cv: ControlVolume, t: float = 0.0, g=(0.0, 0.0, -G0),
                            traction: Callable | None = None, origin=(0.0, 0.0, 0.0), n: int = 24, dt: float = 1e-4,
                            tau: Callable | None = None) -> AngularMomentumBudget:
    """Angular-momentum principle for a stationary control volume, Eq. (4.65).

    Book: §4.9, Eqs. (4.64)–(4.65): d/dt∫_{V_o}(r × ρu)dV + ∫_{A_o}(r × ρu)(u·n)dA = ∫_{V_o}(r × ρg)dV + ∫_{A_o}(r × f)dA.
    Parameters as :func:`momentum_budget` (3-D CV, used as fixed). Returns :class:`AngularMomentumBudget`.
    Validation: V1 rigid rotation in a fixed cylinder with τ = −pδ: storage, outflux and torques about the axis 0.
    Label: analytic.
    """
    rho_f = as_field(rho)
    o = _F(origin).reshape(3, 1)
    H = lambda X, T: np.cross(X - o, ev(rho_f, X, T) * ev(u, X, T, (3,)), axis=0)  # noqa: E731
    storage = (_vol_int(H, cv, t + dt, n, (3,)) - _vol_int(H, cv, t - dt, n, (3,))) / (2.0 * dt)
    Xs, N, dA, _ = cv.surface_nodes(t, n)
    outflux = angular_momentum_flux_surface(rho, u, Xs, N, dA, t, origin)
    Xv, W = cv.volume_nodes(t, n)
    body = np.sum(np.cross(Xv - o, ev(rho_f, Xv, t) * gvec(g, 3, (Xv.shape[1],)), axis=0) * W, axis=-1)
    f = _traction(Xs, N, t, traction, tau, 3)
    surf = np.sum(np.cross(Xs - o, f, axis=0) * dA, axis=-1)
    return AngularMomentumBudget(storage, outflux, body, surf, storage + outflux - body - surf)  # Eq. (4.65)
