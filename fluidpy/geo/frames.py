"""Rotating frames for Vallis Ch. 2 and the Chapter-0B Coriolis explainers: the Earth-frame view of free (inertial)
motion on a turntable and on a sphere, the Coriolis parameter f and its gradient β, the horizontal Coriolis
acceleration −f k̂ × u, and f-plane inertial oscillations.

Vallis: §2.1, Eqs. (2.7) and (2.12) (relative vs inertial velocity; Coriolis and centrifugal terms); §2.3,
Eqs. (2.78)–(2.81) (f-plane and β-plane); §2.9, Eqs. (2.233)–(2.234) (inertial circles of radius U/f, period 2π/f);
§3.8, Eqs. (3.126)–(3.127) (inertial oscillations).
Kundu: §4.7, Eqs. (K 4.42)–(K 4.45) and Fig. 4.8 (the pole projectile, ``fluidpy.core.rotating``); §13.4,
Eqs. (K 13.7)–(K 13.10) and §13.11 "Inertial Motion" (``fluidpy.core.gfd``).

Conventions
-----------
* Latitude ϑ (Vallis's symbol; θ is potential temperature) and longitude λ in **radians** inside every function; the
  ``*_deg`` arguments/returns are the only degree interfaces.
* Tangent plane: x east, y north, z up; k̂ the local vertical. Earth's Ω = 7.292115e-5 s⁻¹ (sidereal, from
  ``fluidpy.geo.constants``), a = 6.371e6 m.
* Rotation sense is always stated with the viewpoint: "clockwise seen from above" means looking down the local
  vertical (northern-hemisphere inertial circles); for f < 0 the sense reverses.
* Physics already in ``fluidpy.core`` is imported, not re-implemented (f, β, inertial oscillation, inertial radius and
  period); no ``core`` signature is changed.
"""
from __future__ import annotations

import numpy as np

from ..core._util import as_scalar_if_0d
from ..core import gfd as _gfd
from ..core import rotating as _rot
from .constants import EARTH_RADIUS, OMEGA
from .vectors import cross, k_cross

__all__ = ["coriolis_parameter", "beta_parameter", "coriolis_acceleration_horizontal", "coriolis_horizontal",
           "turntable_puck", "inertial_circle", "inertial_oscillation", "inertial_oscillation_sphere",
           "to_inertial_lon", "great_circle", "great_circle_earth_frame"]

_S = as_scalar_if_0d
_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731


# ---------------------------------------------------------------------------------------------------------------------
# f, β and the horizontal Coriolis acceleration
# ---------------------------------------------------------------------------------------------------------------------

def coriolis_parameter(lat, Omega: float = OMEGA):
    """Coriolis parameter f = 2Ω sin ϑ — twice the local vertical component of the planet's rotation.

    Vallis: §2.3, Eq. (2.78) (f₀ = 2Ω^z = 2Ω sin ϑ₀) and Eq. (2.80) (f = 2Ω sin ϑ).
    Kundu: §13.4, Eq. (K 13.8); defined once in ``fluidpy.core.rotating.coriolis_parameter`` (called here).

    Parameters
    ----------
    lat : float or array_like — latitude ϑ [rad], positive north.
    Omega : float — planetary rotation rate Ω [rad/s] (default Earth, sidereal).

    Returns
    -------
    f : float or ndarray [s⁻¹]; ≈ 1.03e-4 s⁻¹ at 45°N, 0 at the equator, negative in the southern hemisphere.

    Assumptions: spherical planet; only the local-vertical component of 2Ω (traditional approximation).
    Validation: V1 f(π/2) = 2Ω, f(0) = 0, odd in ϑ; equals ``core.rotating.coriolis_parameter``. Label: analytic.
    """
    return _rot.coriolis_parameter(lat, Omega)  # Eq. (2.80): f = 2Ω sin ϑ


def beta_parameter(lat, Omega: float = OMEGA, a: float = EARTH_RADIUS):
    """β = ∂f/∂y = 2Ω cos ϑ / a — the northward gradient of the Coriolis parameter (β-plane f = f₀ + βy).

    Vallis: §2.3.2, Eqs. (2.80)–(2.81) (β = 2Ω cos ϑ₀ / a); restated with Eq. (2.203).
    Kundu: §13.4, Eq. (K 13.10); calls ``fluidpy.core.gfd.beta_parameter``.

    Parameters
    ----------
    lat : float or array_like — latitude ϑ [rad].
    Omega : float — rotation rate Ω [rad/s].   a : float — planetary radius [m] (default 6.371e6).

    Returns
    -------
    beta : float or ndarray [m⁻¹ s⁻¹]; ≈ 1.6e-11 at 45°, ≥ 0 in both hemispheres, 0 at the poles.

    Assumptions: spherical planet, dy = a dϑ.
    Validation: V1 derivative of 2Ω sin ϑ divided by a (finite/complex step); V7 even in ϑ. Label: analytic.
    """
    return _gfd.beta_parameter(lat, Omega, a)  # Eq. (2.81): β = 2Ω cos ϑ₀ / a


def coriolis_acceleration_horizontal(u, v, f):
    """Horizontal Coriolis acceleration (force per unit mass) −f k̂ × u = (f v, −f u).

    Vallis: §2.3, Eq. (2.79) (Du/Dt + f₀ × u = −∇p/ρ, so the force is −f k × u) and §2.9, Eq. (2.233)
    (−f × u = −fU n: to the right of the motion for f > 0). Kundu: §13.4, Eq. (K 13.7) (the left-hand-side term
    (−f v, f u) is the negative of this); §4.7, Eq. (K 4.45) (−2Ω × u).

    Parameters
    ----------
    u, v : float or array_like — eastward and northward velocity [m/s] (broadcast together).
    f : float or array_like — Coriolis parameter [s⁻¹] (either sign; use :func:`coriolis_parameter`).

    Returns
    -------
    (a_x, a_y) : floats or ndarrays [m/s²]. Perpendicular to u (does no work), magnitude |f||u|, pointing to the
    right of the motion for f > 0 (northern hemisphere) and to the left for f < 0.

    Assumptions: traditional approximation (only the vertical component of Ω), horizontal motion.
    Validation: V1 a·u = 0, |a| = |f||u|; equals the horizontal part of ``core.rotating.coriolis_acceleration``
    with Ω = (0, 0, f/2); V7 sign flips with f. Label: analytic.
    """
    u_, v_, f_ = np.broadcast_arrays(_F(u), _F(v), _F(f))
    kxu = k_cross(np.stack([u_, v_], axis=-1))  # k × u = (−v, u)
    a = -f_[..., None] * kxu  # Eq. (2.79)/(2.233): −f k × u
    return _S(a[..., 0]), _S(a[..., 1])


def coriolis_horizontal(u, v, lat, Omega: float = OMEGA):
    """Horizontal Coriolis acceleration −f k̂ × u = (f v, −f u) at latitude ϑ, with f = 2Ω sin ϑ (the PoC-B name).

    Vallis: §2.2.5, Eq. (2.55) (Du/Dt + f × u = −∇_z p/ρ, f = 2Ω sin ϑ k) and §2.3, Eq. (2.79). Kundu: §13.4,
    Eq. (K 13.7). Calls :func:`coriolis_parameter` and :func:`coriolis_acceleration_horizontal`.

    Parameters
    ----------
    u, v : float or array_like — eastward, northward velocity [m/s].   lat : float or array_like — latitude ϑ [rad].
    Omega : float — rotation rate [rad/s].

    Returns
    -------
    (a_x, a_y) [m/s²] — to the right of u in the NH (ϑ > 0), to the left in the SH, zero on the equator.

    Assumptions: traditional approximation.   Validation: V1 ⟂ u, |a| = |f||u|; V7 mirror in ϑ. Label: analytic.
    """
    return coriolis_acceleration_horizontal(u, v, coriolis_parameter(lat, Omega))  # Eq. (2.55): −f k × u


# ---------------------------------------------------------------------------------------------------------------------
# Free motion on a turntable (flat rotating frame)
# ---------------------------------------------------------------------------------------------------------------------

def _rotate_z(x, y, ang):
    """Coordinates of the point (x, y) in axes turned by +ang (i.e. the point rotated by −ang)."""
    c, s = np.cos(ang), np.sin(ang)
    return x * c + y * s, -x * s + y * c


def turntable_puck(t, x0: float, y0: float, u0: float, v0: float, Omega: float, launch: str = "rotating") -> dict:
    """A frictionless puck on a flat turntable: its straight inertial path and the same path seen by an observer
    riding the table (closed form).

    Vallis: §2.1, Eq. (2.7) v_I = v_R + Ω × r (launch velocity conversion) and Eq. (2.12) (in the rotating frame the
    path obeys dv_R/dt = −2Ω × v_R − Ω × (Ω × r) — the Coriolis and centrifugal terms).
    Kundu: §4.7, Fig. 4.8 and ``fluidpy.core.rotating.projectile_paths`` (the special case x0 = y0 = 0, v0 = 0,
    launch="inertial").

    The puck feels no horizontal force, so in the inertial (laboratory) frame it moves in a straight line at constant
    speed: x_I = x0 + U t, y_I = y0 + V t. The table (and its frame) has turned through Ωt, so the rotating-frame
    coordinates are that point rotated by −Ωt: x_R = x_I cos Ωt + y_I sin Ωt, y_R = −x_I sin Ωt + y_I cos Ωt.

    Parameters
    ----------
    t : float or array_like — time(s) [s], t = 0 at launch (frames coincide at t = 0).
    x0, y0 : float — launch position [m] relative to the rotation axis.
    u0, v0 : float — launch velocity [m/s]: relative to the table if ``launch="rotating"`` (a hand on the table
        pushes it), or inertial (laboratory) if ``launch="inertial"``.
    Omega : float — table rotation rate [rad/s], > 0 anticlockwise seen from above (like the NH Earth), < 0 clockwise.
    launch : "rotating" (default) or "inertial" — the frame of (u0, v0).

    Returns
    -------
    dict with ``t`` [s]; ``x_inertial``, ``y_inertial``, ``x_rotating``, ``y_rotating`` [m] (arrays like t, floats for
    scalar t); ``U``, ``V`` the inertial velocity [m/s]; ``angle`` = Ωt [rad] turned by the table.

    Assumptions: frictionless, flat table (no centripetal support: the puck keeps its inertial straight line).
    Validation: V1 equals the solution of the rotating-frame ODE (2.12) (solve_ivp, DOP853, rtol 1e-10); matches
    ``core.rotating.projectile_paths`` for a launch from the axis; V4 |x_R| = |x_I| (rotation preserves distance from
    the axis); V7 Ωt → 0 gives the inertial path; NH (Ω > 0) deflection to the right. Label: analytic.
    """
    if launch not in ("rotating", "inertial"):
        raise ValueError("launch must be 'rotating' or 'inertial'")
    Om = float(Omega)
    U, V = float(u0), float(v0)
    if launch == "rotating":  # Eq. (2.7): v_I = v_R + Ω × r, Ω × r = Ω(−y, x)
        U, V = U - Om * float(y0), V + Om * float(x0)
    tt = _F(t)
    xi = float(x0) + U * tt  # straight line in the inertial frame
    yi = float(y0) + V * tt
    ang = Om * tt
    xr, yr = _rotate_z(xi, yi, ang)  # rotate the inertial point by −Ωt
    return dict(t=_S(tt), x_inertial=_S(xi), y_inertial=_S(yi), x_rotating=_S(xr), y_rotating=_S(yr),
                U=U, V=V, angle=_S(ang))


# ---------------------------------------------------------------------------------------------------------------------
# f-plane inertial oscillations
# ---------------------------------------------------------------------------------------------------------------------

def inertial_circle(U, f) -> dict:
    """Size, period and sense of an f-plane inertial circle for speed U and Coriolis parameter f.

    Vallis: §2.9.2, after Eq. (2.234) (U/R + f = 0: inertia circles of radius U/f and period 2π/f = 1 day/(2 sin ϑ));
    §3.8.2, Eq. (3.126). Kundu: §13.11 "Inertial Motion"; calls ``core.gfd.inertial_radius`` and
    ``core.gfd.inertial_period``.

    Parameters
    ----------
    U : float — speed |u| [m/s].   f : float — Coriolis parameter [s⁻¹], nonzero (either sign).

    Returns
    -------
    dict: ``radius`` = U/|f| [m]; ``period`` = 2π/|f| [s]; ``omega`` = |f| [rad/s]; ``sense`` —
    "clockwise seen from above" for f > 0 (northern hemisphere), "anticlockwise seen from above" for f < 0;
    ``turns`` — "right" (f > 0) / "left" (f < 0), the side toward which the parcel is deflected.

    Assumptions: f-plane (constant f), no pressure gradient, no friction. Raises ValueError at f = 0.
    Validation: V1 at 45°N and U = 10 m/s, radius ≈ 97 km, period ≈ 16.9 h; V7 sense flips with f. Label: analytic.
    """
    f = float(f)
    r = _gfd.inertial_radius(U, f)  # raises at f = 0
    T = _gfd.inertial_period(f)
    pos = f > 0
    return dict(radius=float(r), period=float(T), omega=abs(f),
                sense="clockwise seen from above" if pos else "anticlockwise seen from above",
                turns="right" if pos else "left")


def inertial_oscillation(t, u0: float, v0: float, f: float, x0: float = 0.0, y0: float = 0.0) -> dict:
    """Velocity and path of a free parcel on the f-plane (Coriolis force only): an inertial circle.

    Vallis: §3.8.2, Eqs. (3.126)–(3.127) (∂u/∂t − f₀v = 0, ∂v/∂t + f₀u = 0; "free particles in a rotating frame");
    §2.9.2 after Eq. (2.234) (radius U/f, period 2π/f). Kundu: §13.11; the solution u + iv = (u₀ + iv₀) e^{−ift} is
    ``core.gfd.inertial_oscillation`` (called here).

    Parameters
    ----------
    t : float or array_like — time(s) [s].   u0, v0 : float — initial velocity east, north [m/s].
    f : float — Coriolis parameter f = 2Ω sin ϑ [s⁻¹] (0 gives the straight line).
    x0, y0 : float — starting position [m] (default origin).

    Returns
    -------
    dict: ``u``, ``v`` [m/s]; ``x``, ``y`` [m] (absolute, start at (x0, y0)); ``xc``, ``yc`` circle centre [m]
    (to the right of the initial motion for f > 0); plus the :func:`inertial_circle` keys when f ≠ 0.

    Assumptions: f-plane, no pressure gradient, no friction (the "inertial" name is historical: the motion is not
    straight in space — Vallis §2.9.2).
    Validation: V1 returns to the start after 2π/|f|; distance from the centre = U/|f|; V4 speed conserved;
    solves (3.127) (finite differences). Label: analytic.
    """
    u, v, x, y = _gfd.inertial_oscillation(t, float(u0), float(v0), float(f))  # Eq. (3.126) solution
    out = dict(u=u, v=v, x=_S(_F(x) + x0), y=_S(_F(y) + y0))
    if f != 0:
        out["xc"] = float(x0) + float(v0) / f  # centre = start − i(u0 + i v0)/f
        out["yc"] = float(y0) - float(u0) / f
        out.update(inertial_circle(np.hypot(u0, v0), f))
    return out


# ---------------------------------------------------------------------------------------------------------------------
# Coriolis-only (inertial) motion on the sphere — the PoC-B globe paths
# ---------------------------------------------------------------------------------------------------------------------

def _wrap_lon(lon):
    """Wrap longitudes [rad] to (−π, π]."""
    w = np.mod(np.asarray(lon, dtype=float) + np.pi, 2.0 * np.pi) - np.pi
    return np.where(w == -np.pi, np.pi, w)


def _sphere_rhs(s, Omega: float, a: float, metric: bool, coriolis: bool):
    lam, th, u, v = s
    m = (2.0 * Omega * np.sin(th) if coriolis else 0.0) + (u * np.tan(th) / a if metric else 0.0)
    return np.array([u / (a * np.cos(th)),  # Eq. (2.28): dλ/dt = u/(a cos ϑ)
                     v / a,  # Eq. (2.28): dϑ/dt = v/a
                     m * v,  # Eq. (2.50a), ∇p = 0: du/dt = (f + u tan ϑ/a) v
                     -m * u])  # Eq. (2.50b), ∇p = 0: dv/dt = −(f + u tan ϑ/a) u


def _rk4(s, h, *args):
    k1 = _sphere_rhs(s, *args)
    k2 = _sphere_rhs(s + 0.5 * h * k1, *args)
    k3 = _sphere_rhs(s + 0.5 * h * k2, *args)
    k4 = _sphere_rhs(s + h * k3, *args)
    return s + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def inertial_oscillation_sphere(lat0: float, lon0: float, u0: float, v0: float, t, dt: float = 60.0,
                                Omega: float = OMEGA, a: float = EARTH_RADIUS, metric: bool = True,
                                coriolis: bool = True, lat_max_deg: float = 89.0):
    """Earth-frame path of a frictionless parcel on the sphere feeling only the Coriolis force (no pressure gradient):
    inertial oscillations with f = 2Ω sin ϑ varying along the path (fixed-step RK4).

    Vallis: §2.2.4, Eqs. (2.50a,b) with ∂p = 0 — du/dt = (f + u tan ϑ/a) v, dv/dt = −(f + u tan ϑ/a) u, f = 2Ω sin ϑ
    (defined after Eq. (2.51)) — and Eq. (2.28) for dλ/dt = u/(a cos ϑ), dϑ/dt = v/a; vector form Eq. (2.55)
    Du/Dt + f × u = 0. f-plane limit: §3.8.2, Eq. (3.126) and :func:`inertial_oscillation`. Kundu: §13.11.

    Parameters
    ----------
    lat0, lon0 : float — start latitude ϑ₀ and longitude λ₀ [rad] (|ϑ₀| < ``lat_max_deg``).
    u0, v0 : float — initial eastward and northward velocity relative to the Earth [m/s].
    t : float or array_like — output time(s) [s]; negative times integrate backward (the "where-from" probe).
    dt : float — RK4 step [s] (default 60 s; the explainer's JS mirror uses the same step and grid, so outputs at
        multiples of dt agree to round-off). Off-grid output times take one partial step from the last grid state.
    Omega : float — rotation rate [rad/s].   a : float — radius [m].
    metric : bool — keep the metric term u tan ϑ / a (default True; with ``coriolis=False`` the path is a great
        circle, :func:`great_circle`).
    coriolis : bool — keep f (default True); False gives the "if the Earth did not spin" ghost.
    lat_max_deg : float — beyond |ϑ| = this [degrees] the (λ, ϑ) form is singular (cos ϑ → 0); later outputs are NaN.

    Returns
    -------
    (lon, lat, u, v) : longitude [rad, wrapped to (−π, π]], latitude [rad], velocity east/north [m/s], each shaped
    like t (floats for scalar t).

    Numerics: classical RK4, fixed step (our choice; deterministic and mirrored in JS). Speed |u| is conserved by the
    equations (the Coriolis and metric terms are ⟂ u); RK4 keeps it to ~1e-9 relative over 2 days at dt = 60 s.
    Assumptions: thin shell, traditional approximation, no friction, no pressure gradient, w = 0. Motion in the
    Earth frame: the space view is :func:`to_inertial_lon`.
    Validation: V7 f-plane limit — on the sphere the loops drift westward by the β-effect, ≈ 15 % of the radius U/|f|
    per period at U = 10 m s⁻¹ and 45° (departure linear in U/(|f|a)); the f-plane circle (radius U/|f|, period
    2π/|f|) is the U/(|f|a) → 0 limit; V1 ``coriolis=False`` reproduces the great circle (closed form) to the integrator error; V3 RK4 error ∝ dt⁴;
    V4 speed conserved; V7 NH right / SH left mirror (v(−ϑ₀) = −v(ϑ₀)). Label: converged.
    """
    tt = _F(t)
    flat = np.atleast_1d(tt).ravel()
    args = (float(Omega), float(a), bool(metric), bool(coriolis))
    lim = np.deg2rad(lat_max_deg)
    s0 = np.array([float(lon0), float(lat0), float(u0), float(v0)])
    out = np.full((flat.size, 4), np.nan)
    h = float(dt)
    if h <= 0:
        raise ValueError("dt must be > 0")
    for sgn in (1.0, -1.0):
        idx = np.where(flat >= 0)[0] if sgn > 0 else np.where(flat < 0)[0]
        if idx.size == 0:
            continue
        order = idx[np.argsort(np.abs(flat[idx]))]
        s, n, dead = s0.copy(), 0, False
        for k in order:
            T = abs(flat[k])
            while not dead and (n + 1) * h <= T * (1 + 1e-12) + 1e-9:
                s = _rk4(s, sgn * h, *args)
                n += 1
                dead = abs(s[1]) > lim
            if dead:
                continue
            rem = T - n * h
            sk = _rk4(s, sgn * rem, *args) if rem > 1e-9 else s
            out[k] = sk if abs(sk[1]) <= lim else np.nan
    out[:, 0] = _wrap_lon(out[:, 0])
    shp = tt.shape
    return tuple(_S(out[:, i].reshape(shp)) for i in range(4))


def to_inertial_lon(lon, t, Omega: float = OMEGA):
    """Longitude of an Earth-frame point as seen from space (the inertial frame that coincided with the Earth frame at
    t = 0): λ_I = λ_E + Ωt; latitude is unchanged.

    Vallis: §2.1, Eq. (2.7) (v_I = v_R + Ω × r: the planet carries every point eastward at Ω a cos ϑ, i.e. dλ/dt
    gains Ω). Kundu: §4.7, Eq. (K 4.42).

    Parameters
    ----------
    lon : float or array_like — Earth-frame longitude λ_E [rad].   t : float or array_like — time [s] (broadcast).
    Omega : float — rotation rate [rad/s].

    Returns
    -------
    λ_I [rad], wrapped to (−π, π].

    Validation: V1 offset (λ_I − λ_E) mod 2π = Ωt mod 2π; inverse with −Ω. Label: analytic.
    """
    return _S(_wrap_lon(_F(lon) + float(Omega) * _F(t)))  # Eq. (2.7): λ_I = λ_E + Ωt


def great_circle(lat0: float, lon0: float, u0: float, v0: float, t, a: float = EARTH_RADIUS):
    """Great-circle path at constant speed on a non-rotating sphere — the PoC-B "f = 0 ghost" (closed form).

    Vallis: §2.2.4, Eqs. (2.50a,b) with f = 0 and ∇p = 0 (only the metric terms remain; their solutions are great
    circles) and Eq. (2.28). Kundu: none (geometry). Calls :func:`great_circle_earth_frame` with Ω = 0.

    Parameters
    ----------
    lat0, lon0 : float — start [rad].   u0, v0 : float — initial velocity east/north [m/s].
    t : float or array_like — time(s) [s] (negative allowed).   a : float — radius [m].

    Returns
    -------
    (lon, lat) [rad] (lon wrapped to (−π, π]), shaped like t; the start point for every t if u0 = v0 = 0.

    Assumptions: perfect sphere, no rotation, no forces along the surface.
    Validation: V1 |r| = a, closes after 2πa/U; equals :func:`inertial_oscillation_sphere` with ``coriolis=False``
    to the RK4 error. Label: analytic.
    """
    U = float(np.hypot(u0, v0))
    tt = _F(t)
    if U == 0.0:
        return _S(np.full(tt.shape, _wrap_lon(lon0))), _S(np.full(tt.shape, float(lat0)))
    g = great_circle_earth_frame(tt, lat0, lon0, float(np.arctan2(u0, v0)), U, Omega=0.0, a=a, launch="inertial")
    return g["lon"], g["lat"]


# ---------------------------------------------------------------------------------------------------------------------
# A great circle (space-straight path) on the sphere seen from the rotating Earth
# ---------------------------------------------------------------------------------------------------------------------

def _unit_vectors(lat: float, lon: float):
    cl, sl, co, so = np.cos(lat), np.sin(lat), np.cos(lon), np.sin(lon)
    r = np.array([cl * co, cl * so, sl])
    east = np.array([-so, co, 0.0])
    north = np.array([-sl * co, -sl * so, cl])
    return r, east, north


def great_circle_earth_frame(t, lat0: float, lon0: float, azimuth: float, speed: float, Omega: float = OMEGA,
                             a: float = EARTH_RADIUS, launch: str = "earth") -> dict:
    """A projectile skimming a frictionless sphere along a great circle (straight "as seen from space") and the same
    path in the rotating Earth frame, as latitude/longitude over time (closed form).

    **Not used by PoC B** (the globe explainer uses :func:`inertial_oscillation_sphere`, :func:`to_inertial_lon` and
    :func:`great_circle`): on a bulge-free frictionless sphere this puck slides toward the equator under the
    unbalanced centrifugal effect, which misleads. Kept for a labelled Like…/Unlike… "trap" card only.

    Vallis: §2.1, Eq. (2.7) (v_I = v_R + Ω × r at launch) and Eq. (2.12) (the Coriolis and centrifugal terms that bend
    the Earth-frame path); §2.2, Eq. (2.28) (λ, ϑ coordinates). Kundu: §4.7, Fig. 4.8 (the same idea at the pole,
    ``core.rotating.projectile_paths``).

    Inertial frame: the motion stays on the sphere |r| = a with no tangential force, so it is a great circle at
    constant speed V: r_I(t) = a [cos(Vt/a) r̂₀ + sin(Vt/a) d̂], with r̂₀ the launch point and d̂ the unit inertial
    launch direction. Earth frame: the planet has turned through Ωt about its axis, so r_E(t) = R_z(−Ωt) r_I(t)
    (same latitude, longitude reduced by Ωt).

    Parameters
    ----------
    t : float or array_like — time(s) since launch [s] (frames coincide at t = 0).
    lat0, lon0 : float — launch latitude ϑ₀ and longitude λ₀ [rad].
    azimuth : float — launch heading [rad], measured clockwise from north (0 north, π/2 east).
    speed : float — launch speed [m/s] (relative to the ground if ``launch="earth"``, inertial if "inertial").
    Omega : float — rotation rate [rad/s] (default Earth).   a : float — sphere radius [m].
    launch : "earth" (default; the ground's eastward speed Ωa cos ϑ₀ is added, Eq. (2.7)) or "inertial".

    Returns
    -------
    dict: ``lat`` [rad] and ``lon`` [rad, wrapped to (−π, π]] of the Earth-frame path; ``lat_inertial``,
    ``lon_inertial`` [rad] of the great circle itself; ``xyz_earth``, ``xyz_inertial`` [m], shape (..., 3);
    ``speed_inertial`` V [m/s]; ``period`` 2πa/V [s] (one inertial lap). Scalar t gives 0-d entries.

    Assumptions: perfect sphere, frictionless, gravity exactly radial, no centrifugal support — a "space-straight"
    teaching idealisation. On the real (oblate, in-equilibrium) Earth the horizontal part of effective gravity cancels
    the centrifugal force and a free puck instead executes inertial oscillations (:func:`inertial_oscillation`; Vallis
    §2.9.2 notes that "inertial" circles are not truly inertial).
    # DEVIATION: the explainer's idealised sphere, not Vallis's thin-shell equations with effective gravity — chosen to
    #   show a path that is straight in space; labelled "idealised" wherever it is drawn.
    Validation: V1 |r| = a, lap closes after 2πa/V, latitude of r_E equals latitude of r_I, lon difference = −Ωt;
    V7 Ω → 0 gives the great circle; a short-time NH launch deflects to the right of the ground track. Label: analytic.
    """
    if launch not in ("earth", "inertial"):
        raise ValueError("launch must be 'earth' or 'inertial'")
    r0, east, north = _unit_vectors(float(lat0), float(lon0))
    d = np.cos(azimuth) * north + np.sin(azimuth) * east  # heading clockwise from north
    vI = float(speed) * d
    if launch == "earth":  # Eq. (2.7): v_I = v_R + Ω × r (= Ω a cos ϑ₀ eastward)
        vI = vI + cross(np.array([0.0, 0.0, float(Omega)]), a * r0)
    V = float(np.linalg.norm(vI))
    if V == 0.0:
        raise ValueError("zero inertial speed: the projectile does not move in space")
    dI = vI / V
    tt = _F(t)
    s = V * tt / a  # arc angle travelled
    xyz_I = a * (np.cos(s)[..., None] * r0 + np.sin(s)[..., None] * dI)
    ang = float(Omega) * tt
    xe, ye = _rotate_z(xyz_I[..., 0], xyz_I[..., 1], ang)  # R_z(−Ωt) r_I
    xyz_E = np.stack([xe, ye, xyz_I[..., 2]], axis=-1)
    lat_i = np.arcsin(np.clip(xyz_I[..., 2] / a, -1.0, 1.0))
    lon_i = np.arctan2(xyz_I[..., 1], xyz_I[..., 0])
    lon_e = np.arctan2(ye, xe)
    return dict(lat=_S(lat_i), lon=_S(lon_e), lat_inertial=_S(lat_i), lon_inertial=_S(lon_i),
                xyz_earth=xyz_E, xyz_inertial=xyz_I, speed_inertial=V, period=2.0 * np.pi * a / V)
