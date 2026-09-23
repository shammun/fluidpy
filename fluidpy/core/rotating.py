"""Noninertial (translating and rotating) frames: velocity and acceleration transformation, apparent body forces,
Coriolis and centrifugal effects, effective gravity.

Book: Kundu, Cohen & Dowling 5e, Ch. 4 §4.7, Eqs. (4.42)–(4.45), Figs. 4.6–4.9, Example 4.5 (rendered pages
chapters/pages/ch04/p143–p148).

Conventions (``knowledge/notation.md``)
---------------------------------------
* Primes denote the **noninertial frame** O′1′2′3′ (position x′, velocity u′, acceleration a′). (Elsewhere primes mean
  rotated axes (ch02), a translating frame (ch03) or perturbations (§4.9) — read the context.)
* O′ moves at U(t) = dX/dt and rotates at Ω(t) relative to the inertial O123; one clock (t = t′).
* **Sign language.** The *acceleration* terms of (4.43)/(4.44) are +2Ω × u′ (Coriolis) and +Ω × (Ω × x′)
  (centripetal); moved to the right-hand side of (4.45) they become the apparent *forces per unit mass*
  −2Ω × u′ (Coriolis force) and −Ω × (Ω × x′) = +Ω²R e_R (centrifugal). The book uses one name for both: (4.43)
  calls +2Ω × u′ "the Coriolis acceleration", while the text after (4.45) calls −2Ω × u "the Coriolis acceleration"
  that deflects a particle to the right in the NH (the sign there is right; only the name is shared — analysis §9).
* Earth: Ω = 7.292115e-5 rad/s (WGS-84), Ω_z > 0 in the northern hemisphere (deflection to the right).
"""
from __future__ import annotations

import numpy as np

from ._util import as_scalar_if_0d
from .tensors import rotation_matrix_3d
from .thermo import G_BOOK as G0

__all__ = ["OMEGA_EARTH", "EARTH_A", "EARTH_B", "EARTH_RADIUS", "rotating_basis", "basis_rate", "basis_rate_exact",
           "inertial_velocity",
           "frame_acceleration_terms", "apparent_body_forces", "coriolis_acceleration", "coriolis_force",
           "centrifugal_acceleration", "centrifugal_potential", "effective_gravity", "earth_oblateness_diameter",
           "coriolis_parameter", "projectile_paths", "coriolis_deflection"]

OMEGA_EARTH: float = 7.292115e-5   #: Earth's rotation rate [rad/s] (WGS-84)
EARTH_A: float = 6378137.0         #: WGS-84 equatorial radius [m]
EARTH_B: float = 6356752.3142      #: WGS-84 polar radius [m]
EARTH_RADIUS: float = EARTH_A      #: radius of the spherical Earth used by :func:`effective_gravity` [m]

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d


def _vec3(v) -> np.ndarray:
    """A 3-vector (or (3, N) array); a scalar Ω means rotation about z."""
    a = _F(v)
    if a.ndim == 0:
        return np.array([0.0, 0.0, float(a)])
    return a


def _bc(a, b):
    a, b = _vec3(a), _vec3(b)
    if a.ndim < b.ndim:
        a = a.reshape(a.shape + (1,) * (b.ndim - a.ndim))
    if b.ndim < a.ndim:
        b = b.reshape(b.shape + (1,) * (a.ndim - b.ndim))
    return a, b


def _cross(a, b) -> np.ndarray:
    a, b = _bc(a, b)
    return np.cross(a, b, axis=0)


def rotating_basis(Omega, t: float, E0=None) -> np.ndarray:
    """Unit vectors e′_i(t) of a frame rotating at constant Ω from E0 at t = 0 (rows of the result).

    Book: §4.7, Fig. 4.7: each e′_i sweeps a cone about Ω with de′_i/dt = Ω × e′_i (our D14).
    Parameters: Omega : 3-vector [rad/s] (scalar → about z); t [s]; E0 : (3, 3) rows at t = 0 (default identity).
    Returns (3, 3) rows e′_1, e′_2, e′_3. Validation: V3 central difference of e′_i = Ω × e′_i (order 2). Label: analytic.
    """
    Om = _vec3(Omega)
    w = float(np.linalg.norm(Om))
    E = np.eye(3) if E0 is None else _F(E0)
    if w == 0.0:
        return E.copy()
    Rm = rotation_matrix_3d(Om / w, w * float(t))  # active rotation by angle |Ω|t about Ω̂
    return (Rm @ E.T).T


def basis_rate(Omega, t: float, h: float = 1e-6) -> np.ndarray:
    """de′_i/dt of the rotating unit vectors by a central difference of :func:`rotating_basis` (rows) [1/s] — equal to
    Ω × e′_i (Fig. 4.7, text before (4.43); our D14 check). Parameters: Omega [rad/s]; t [s]; h : time step [s].
    Validation: V1/V3 equals :func:`basis_rate_exact` to O(h²). Book: §4.7, Fig. 4.7. Label: converged."""
    return (rotating_basis(Omega, t + h) - rotating_basis(Omega, t - h)) / (2.0 * h)


def basis_rate_exact(Omega, E) -> np.ndarray:
    """Ω × e′_i for each row of E — the geometric result de′_i/dt = Ω × e′_i of Fig. 4.7. Book: §4.7, Fig. 4.7. Label: analytic."""
    Om = _vec3(Omega)
    return np.stack([np.cross(Om, e) for e in _F(E)])


def inertial_velocity(U, u_prime, Omega, x_prime) -> np.ndarray:
    """Velocity in the inertial frame from rotating-frame quantities: u = U + u′ + Ω × x′, Eq. (4.42).

    Book: §4.7, Eq. (4.42) (x = X + x′ differentiated; the dx′_i/dt e′_i part is u′, the x′_i de′_i/dt part is
    Ω × x′). Parameters (3-vectors or (3, N)): U [m/s] frame-origin velocity, u′ [m/s], Ω [rad/s], x′ [m].
    Returns u [m/s]. Validation: V1 a point fixed in the rotating frame moves at Ω × x′. Label: analytic.
    """
    Up, up = _bc(U, u_prime)
    return Up + up + _cross(Omega, x_prime)  # Eq. (4.42)


def frame_acceleration_terms(a_prime, u_prime, x_prime, Omega, dOmega_dt=0.0, dU_dt=0.0) -> dict:
    """The five parts of a fluid particle's inertial acceleration seen from a noninertial frame, Eq. (4.43)/(4.44).

    Book: §4.7, Eq. (4.43): a = dU/dt + a′ + 2Ω × u′ + (dΩ/dt) × x′ + Ω × (Ω × x′) — frame acceleration, relative
    acceleration, Coriolis acceleration, angular-acceleration term, centripetal acceleration (our D24 ★★★ derives it:
    one Ω × u′ comes from the turning axes, the other from the moving position).

    Parameters
    ----------
    a_prime [m/s²], u_prime [m/s], x_prime [m] : in the rotating frame (3-vectors or (3, N));
    Omega [rad/s] (scalar → about z); dOmega_dt [rad/s²]; dU_dt [m/s²] (frame-origin acceleration)

    Returns
    -------
    dict(frame=dU/dt, relative=a′, coriolis=2Ω × u′, angular=(dΩ/dt) × x′, centripetal=Ω × (Ω × x′), total=a) [m/s²]

    Validation: V1 a path given in the inertial frame, differentiated in the rotating frame, plus these terms equals
    its inertial acceleration (200 random cases); V2 sympy re-derivation with general U(t), Ω(t); the "no factor 2"
    variant fails. Label: analytic, symbolic.
    """
    ap, up = _bc(a_prime, u_prime)
    fr = _bc(dU_dt, ap)[0] * np.ones_like(ap)
    cor = 2.0 * _cross(Omega, up)
    ang = _cross(dOmega_dt, x_prime)
    cen = _cross(Omega, _cross(Omega, x_prime))
    ang, _ = _bc(ang, ap)
    cen, _ = _bc(cen, ap)
    return {"frame": fr, "relative": ap, "coriolis": cor, "angular": ang, "centripetal": cen,
            "total": fr + ap + cor + ang + cen}  # Eq. (4.43)


def apparent_body_forces(u_prime, x_prime, Omega, dOmega_dt=0.0, dU_dt=0.0, g=(0.0, 0.0, -G0)) -> dict:
    """Body forces per unit mass in the noninertial momentum equation: the bracket of Eq. (4.45).

    Book: §4.7, Eq. (4.45): ρD′u′/Dt = −∇′p + ρ[g − dU/dt − 2Ω × u′ − (dΩ/dt) × x′ − Ω × (Ω × x′)] + μ∇′²u′ (derived
    from the incompressible, constant-μ (4.39b); ∇′²u′ = ∇²u needs U, Ω uniform in space — a rigid frame).

    Returns dict(gravity=g, frame=−dU/dt, coriolis=−2Ω × u′, angular=−(dΩ/dt) × x′, centrifugal=−Ω × (Ω × x′),
    total) [m/s²]. With Ω = (0, 0, Ω): centrifugal = Ω²R e_R (points away from the axis).
    Validation: V1 total = g − (terms of :func:`frame_acceleration_terms` other than a′); inertial frame → g only.
    Label: analytic.
    """
    up = _vec3(u_prime)
    gv, _ = _bc(_F(g), up)
    fr = -_bc(dU_dt, up)[0] * np.ones_like(up)
    cor = -2.0 * _cross(Omega, up)
    ang, _ = _bc(-_cross(dOmega_dt, x_prime), up)
    cen, _ = _bc(-_cross(Omega, _cross(Omega, x_prime)), up)
    gv = gv * np.ones_like(up)
    return {"gravity": gv, "frame": fr, "coriolis": cor, "angular": ang, "centrifugal": cen,
            "total": gv + fr + cor + ang + cen}  # Eq. (4.45), [·]-bracket


def coriolis_acceleration(Omega, u) -> np.ndarray:
    """Apparent Coriolis acceleration felt in the rotating frame — the Coriolis **force per unit mass** −2Ω × u′ of the
    bracket of (4.45) [m/s²] (what the book's text after Fig. 4.7 calls "the Coriolis acceleration −2Ω × u").

    Book: §4.7, Eq. (4.45). The *kinematic* term of (4.43)/(4.44) is +2Ω × u′ (``frame_acceleration_terms["coriolis"]``)
    — the same vector with the opposite sign, moved to the other side of the equation (analysis §9, sign language).
    It depends on the velocity, not the position, and is ⟂ u (does no work). Example: Ω = (0, 0, Ω), u = (u, 0, 0) →
    (0, −2Ωu, 0): to the right of the motion in the NH; cylindrical form −2Ω_z u_R e_φ for flow out of a high.
    Validation: V1 the "+2Ω × u" (no sign flip) and "−Ω × u" (no factor 2) variants fail the projectile test.
    Label: analytic.
    """
    return -2.0 * _cross(Omega, u)  # Eq. (4.45): −2Ω × u′


def coriolis_force(Omega, u) -> np.ndarray:
    """Alias of :func:`coriolis_acceleration`: −2Ω × u′ [m/s² = N/kg] (Eq. (4.45)). Book: §4.7, Eq. (4.45). Label: analytic."""
    return coriolis_acceleration(Omega, u)


def centrifugal_acceleration(Omega, x) -> np.ndarray:
    """Centrifugal apparent acceleration −Ω × (Ω × x′) [m/s²]; for Ω = (0, 0, Ω) it is +Ω²R e_R.

    Book: §4.7, Eq. (4.45) and the text after Fig. 4.9 (the (4.43) term Ω × (Ω × x′) is the *centripetal*
    acceleration; this is its negative). Validation: V2 sympy = −∇(−½Ω²R²). Label: analytic.
    """
    return -_cross(Omega, _cross(Omega, x))


def centrifugal_potential(R, Omega):
    """Potential of the centrifugal acceleration, Φ_c = −½Ω²R² [J/kg], so −∇Φ_c = Ω²R e_R (Exercise 4.43; our D27).

    Book: §4.7, text after Fig. 4.9 ("a body-force potential for the new term can be found"). Label: analytic.
    """
    return _S(-0.5 * _F(Omega) ** 2 * _F(R) ** 2)


def effective_gravity(lat_rad, g_n: float = 9.8, Omega: float = OMEGA_EARTH, a: float = EARTH_A):
    """Effective gravity g_e = g + Ω²R e_R on a spherical Earth at latitude φ: magnitude and deflection from the radial.

    Book: §4.7, Fig. 4.9 and the text after it: g_e = g + Ω²R e_R, R = a cos φ (distance from the axis); the
    equipotentials are ⟂ g_e and sea level is one of them.

    Parameters
    ----------
    lat_rad : latitude φ [rad] (float or array);  g_n : Newtonian gravity toward the centre [m/s²];
    Omega : rotation rate [rad/s];  a : radius [m]

    Returns
    -------
    (g_e, deflection) : |g_e| [m/s²] and the angle between g_e and the radial (inward) direction [rad], ≥ 0 (toward
    the equator); both scalar-callable.

    Assumptions: spherical Earth with radial Newtonian gravity (the bulge itself is ignored).
    Validation: V1 pole: g_e = g_n, deflection 0; equator: g_e = g_n − Ω²a, deflection 0; deflection peaks near 45°;
    WGS-84 Ω²a = 0.0339 m/s². Label: analytic.
    """
    phi = _F(lat_rad)
    c = float(Omega) ** 2 * float(a) * np.cos(phi)  # |Ω²R e_R| with R = a cos φ
    radial_in = float(g_n) - c * np.cos(phi)  # inward-radial component of g_e
    tangential = c * np.sin(phi)  # toward the equator
    return _S(np.hypot(radial_in, tangential)), _S(np.arctan2(tangential, radial_in))


def earth_oblateness_diameter(a: float = EARTH_A, b: float = EARTH_B) -> float:
    """Equatorial minus polar diameter 2(a − b) [m] of the reference ellipsoid (WGS-84 default ≈ 42.77 km; the book
    rounds to 42 km). Book: §4.7, text after Fig. 4.9. Label: benchmark (WGS-84 values)."""
    return 2.0 * (float(a) - float(b))


def coriolis_parameter(lat_rad, Omega: float = OMEGA_EARTH):
    """Coriolis parameter f = 2Ω sin φ [1/s] — the local vertical component of 2Ω (forward pointer to Ch. 13).

    Book: §4.7 (Ω_z > 0 in the NH, < 0 in the SH); f is named in Ch. 13. Label: analytic.
    """
    return _S(2.0 * float(Omega) * np.sin(_F(lat_rad)))


def projectile_paths(u0: float, Omega: float, t):
    """A projectile launched horizontally from the pole (rotation axis) along x at speed u0, frictionless and with
    gravity balanced: its straight inertial path and the same path seen from the rotating frame (Fig. 4.8).

    Book: §4.7, Fig. 4.8: forward distance ut, deflection ≈ Ωut², angular deflection Ωt — "the projectile in fact
    travels in a straight line if observed from outer space".

    Parameters
    ----------
    u0 : launch speed [m/s];  Omega : rotation rate about +z [rad/s] (> 0 NH, < 0 SH);  t : time(s) [s]

    Returns
    -------
    (inertial_xy, rotating_xy) : arrays (2, N) [m] (x′ = u0 t cos Ωt, y′ = −u0 t sin Ωt; NH → y′ < 0: to the right).

    Validation: V1 equals the solution of the rotating-frame ODE a′ = −2Ω × u′ − Ω × (Ω × x′)
    (``ch04.coriolis_projectile``) to 1e-9; |y′| → Ωu0t² as Ωt → 0. Label: analytic.
    """
    t_ = np.atleast_1d(_F(t))
    xi = np.stack([u0 * t_, 0.0 * t_])
    th = float(Omega) * t_
    xr = np.stack([u0 * t_ * np.cos(th), -u0 * t_ * np.sin(th)])  # inertial point seen in axes turned by Ωt
    return xi, xr


def coriolis_deflection(u: float, Omega: float, t, exact: bool = False):
    """Sideways deflection of the pole projectile: Ωut² (the book's small-Ωt estimate) or, with ``exact=True``, the
    true lateral offset ut sin Ωt seen in the rotating frame [m]. Book: §4.7, text before Fig. 4.8. Scalar-callable.
    Validation: V1 exact/estimate = sin(Ωt)/(Ωt) → 1 as Ωt → 0 (relative error (Ωt)²/6). Label: analytic.
    """
    t_ = _F(t)
    if exact:
        return _S(float(u) * t_ * np.sin(abs(float(Omega)) * t_))
    return _S(abs(float(Omega)) * float(u) * t_ ** 2)
