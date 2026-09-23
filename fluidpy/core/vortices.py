"""Axisymmetric plane vortices: solid-body rotation, the irrotational (line) vortex, the Rankine and Gaussian vortices,
their vorticity and circulation, and the polar-coordinate vorticity formula.

Book: Ch. 3 §3.5, Eqs. (3.22)–(3.29), Figs. 3.15–3.16 (rendered pages chapters/pages/ch03/p110–p112). Plane polar
coordinates (r, θ), counter-clockwise positive; u_r radial and u_θ azimuthal velocity [m/s]; ω_z the vorticity normal
to the plane [1/s]; Γ the **circulation** [m²/s] (not the ch01 lapse rate or the ch02 shear rate; ch03 §3.5's shear
rate is γ). Profile callables of r alone are accepted wherever an axisymmetric u_θ(r) is meant; polar fields that
depend on θ are callables f(r, θ).

Reuse: Γ(r) of the Gaussian vortex, the Lamb–Oseen vortex (σ² = 4νt, Ch. 5/8), cyclone and tornado models (Ch. 13),
wing-tip vortices (Ch. 14).
"""
from __future__ import annotations

import inspect
from typing import Callable

import numpy as np
import sympy as sp
from scipy.optimize import brentq
from scipy.special import lambertw

from ._util import as_scalar_if_0d

__all__ = ["solid_body_rotation", "line_vortex", "rankine_vortex", "gaussian_vortex", "gaussian_vortex_max_radius",
           "vortex_profile", "VORTEX_KINDS", "polar_vorticity_z", "polar_vorticity_z_sym", "circulation_circle",
           "mean_vorticity_in_disc", "vortex_velocity_field"]

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d


def solid_body_rotation(r, omega0):
    """Solid-body rotation u_r = 0, u_θ = ω₀ r, Eq. (3.22).

    Book: §3.5, Eq. (3.22), Fig. 3.15 (a tank spun up until the transients die out). ω₀ [rad/s] is the angular
    velocity of every particle about the origin; the vorticity (3.23) is 2ω₀ everywhere and S = 0.
    Parameters: r [m] (≥ 0), ω₀ [rad/s]. Returns u_θ [m/s].
    Validation: V1 Cartesian twin = ch02 ``solid_body_rotation_field``; ω = 2ω₀. Label: analytic.
    """
    return _S(float(omega0) * _F(r))  # Eq. (3.22): u_θ = ω₀ r


def line_vortex(r, B):
    """Irrotational (ideal line) vortex u_r = 0, u_θ = B/r, Eq. (3.25).

    Book: §3.5, Eq. (3.25), Fig. 3.16 (the figure writes C for B; we use B). Γ = 2πB around the origin for every r
    (3.26); ω_z = 0 for r > 0 and a δ-function core at r = 0 (3.27). Parameters: r [m], B [m²/s].
    Returns u_θ [m/s]; NaN at r = 0 (singular).
    Validation: V1 Γ independent of r; (1/r) d(r u_θ)/dr = 0. Label: analytic.
    """
    r_ = _F(r)
    with np.errstate(divide="ignore", invalid="ignore"):
        u = np.where(r_ > 0.0, float(B) / np.where(r_ > 0.0, r_, 1.0), np.nan)  # Eq. (3.25): u_θ = B/r
    return _S(u)


def rankine_vortex(r, Gamma, sigma):
    """Rankine vortex: uniform vorticity core of radius σ, irrotational outside, Eq. (3.28).

    Book: §3.5, Eq. (3.28): ω_z = Γ/πσ² (r ≤ σ), 0 (r > σ); u_θ = (Γ/2πσ²) r (r ≤ σ), Γ/2πr (r > σ). u_θ is
    continuous and maximal (Γ/2πσ) at r = σ; the circulation of any circle r ≥ σ is Γ.

    Parameters
    ----------
    r : radius [m] (float or array, ≥ 0);  Gamma : total circulation Γ [m²/s];  sigma : core radius σ [m] (> 0)

    Returns
    -------
    (u_theta [m/s], omega_z [1/s])

    Validation: V1 continuity at σ, max at σ, Γ(r ≥ σ) = Γ; V2 (3.23) applied to each piece gives ω_z; form
    cross-check against the standard Rankine-vortex profile (not a V5 benchmark). Label: analytic, symbolic.
    """
    r_ = _F(r)
    G, s = float(Gamma), float(sigma)
    inside = r_ <= s
    with np.errstate(divide="ignore", invalid="ignore"):
        u = np.where(inside, G * r_ / (2.0 * np.pi * s ** 2), G / (2.0 * np.pi * np.where(inside, 1.0, r_)))  # (3.28)
    w = np.where(inside, G / (np.pi * s ** 2), 0.0)  # Eq. (3.28): ω_z = Γ/πσ² inside, 0 outside
    return _S(u), _S(w)


def gaussian_vortex(r, Gamma, sigma):
    """Gaussian vortex: smooth core, Eq. (3.29).

    Book: §3.5, Eq. (3.29): ω_z = (Γ/πσ²) exp(−r²/σ²), u_θ = (Γ/2πr)(1 − exp(−r²/σ²)). Evaluated as
    −expm1(−r²/σ²)/r so there is no cancellation near r = 0 (u_θ(0) = 0 exactly). It is the Lamb–Oseen vortex with
    σ² = 4νt (Ch. 5/8 pointer); r ≪ σ → solid body with ω₀ = Γ/2πσ², r ≫ σ → line vortex Γ/2πr; u_θ peaks at
    r ≈ 1.1209σ (:func:`gaussian_vortex_max_radius`).

    Parameters / Returns: as :func:`rankine_vortex`.

    Validation: V1 Γ(r) = Γ(1 − e^{−r²/σ²}); limits; rel. error < 1e-14 at r = 1e-9σ; V2 (3.23) gives ω_z from u_θ;
    form cross-check against the standard Lamb–Oseen profile with σ² = 4νt (not a V5 benchmark).
    Label: analytic, symbolic.
    """
    r_ = _F(r)
    G, s = float(Gamma), float(sigma)
    x = (r_ / s) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        u = np.where(r_ > 0.0, G / (2.0 * np.pi) * (-np.expm1(-x)) / np.where(r_ > 0.0, r_, 1.0), 0.0)  # Eq. (3.29)
    w = G / (np.pi * s ** 2) * np.exp(-x)  # Eq. (3.29): ω_z = (Γ/πσ²) e^{−r²/σ²}
    return _S(u), _S(w)


def gaussian_vortex_max_radius(sigma: float = 1.0, method: str = "brentq") -> float:
    """Radius of the maximum of u_θ in the Gaussian vortex: r_max ≈ 1.1209σ (§3.5, Exercise 3.26; D20).

    Setting du_θ/dr = 0 for (3.29) with x = r²/σ² gives 1 + 2x = eˣ (the root x = 0 is the trivial one); the other
    root is x* = 1.25643… and r_max = σ√x*.

    method
        ``"brentq"`` — root of 1 + 2x − eˣ bracketed in (0.5, 3) (sign change: +0.35 at 0.5, −13.1 at 3);
        ``"lambertw"`` — closed form x* = −W₋₁(−e^{−1/2}/2) − ½ (``scipy.special.lambertw``, branch −1).

    Returns r_max [m] (same unit as σ).

    Validation: V1 residual < 1e-14, the two methods agree to 1e-13; V5 α = r_max²/σ² = 1.256 (Canivete Cuissa &
    Steiner 2022, Eq. 10); V6 book value (private JSON). Label: analytic, benchmark.
    Book: §3.5 after Eq. (3.29) (Exercise 3.26).
    """
    if method == "brentq":
        x = brentq(lambda x: 1.0 + 2.0 * x - np.exp(x), 0.5, 3.0, xtol=1e-15, rtol=4 * np.finfo(float).eps)  # 1 + 2x = eˣ
    elif method == "lambertw":
        x = float(np.real(-lambertw(-0.5 * np.exp(-0.5), k=-1)) - 0.5)  # x* = −W₋₁(−e^{−1/2}/2) − ½
    else:
        raise ValueError('method must be "brentq" or "lambertw"')
    return float(sigma) * float(np.sqrt(x))


VORTEX_KINDS = ("solid", "line", "rankine", "gaussian")


def vortex_profile(kind: str, r, Gamma: float = 2 * np.pi, sigma: float = 1.0):
    """One dispatcher over the four §3.5 profiles, all parametrised by (Γ, σ) so that they can be compared.

    * ``"solid"`` — solid-body rotation (3.22) with ω₀ = Γ/(2πσ²), i.e. the circulation of the circle r = σ is Γ
      (our parametrisation for comparability; the book gives ω₀ directly); ω_z = 2ω₀ = Γ/πσ² everywhere.
    * ``"line"`` — irrotational vortex (3.25) with B = Γ/2π (Γ = 2πB, (3.26)); ω_z = 0 for r > 0 (NaN at r = 0 for u_θ).
    * ``"rankine"`` — (3.28);  ``"gaussian"`` — (3.29).

    Returns (u_theta [m/s], omega_z [1/s]). Units: r, σ [m]; Γ [m²/s]. Scalar-callable (explainer parity).
    Label: analytic.
    Book: §3.5, Eqs. (3.22), (3.25), (3.28), (3.29).
    """
    G, s = float(Gamma), float(sigma)
    if kind == "solid":
        w0 = G / (2.0 * np.pi * s ** 2)
        return solid_body_rotation(r, w0), _S(np.full(np.shape(r), 2.0 * w0))
    if kind == "line":
        return line_vortex(r, G / (2.0 * np.pi)), _S(np.zeros(np.shape(r)))
    if kind == "rankine":
        return rankine_vortex(r, G, s)
    if kind == "gaussian":
        return gaussian_vortex(r, G, s)
    raise ValueError(f"kind must be one of {VORTEX_KINDS}")


def _polar_fn(f: Callable | None) -> Callable:
    """Accept f(r) (axisymmetric) or f(r, θ); None → 0."""
    if f is None:
        return lambda r, th: np.zeros(np.broadcast(np.asarray(r), np.asarray(th)).shape)
    try:
        n_args = len([p for p in inspect.signature(f).parameters.values()
                      if p.default is inspect.Parameter.empty and p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD)])
    except (TypeError, ValueError):
        n_args = 2
    if n_args == 1:
        return lambda r, th: _F(f(r)) + 0.0 * _F(th)
    return f


def _profile_fn(u_theta, **kw) -> Callable:
    """A profile given as a callable u_θ(r) (or f(r, θ)) or as a kind name of :func:`vortex_profile` (keywords
    ``Gamma``, ``sigma``) → callable."""
    if isinstance(u_theta, str):
        kind = u_theta
        return lambda r, *_: vortex_profile(kind, r, **kw)[0]  # noqa: E731
    return u_theta


def polar_vorticity_z(u_r: Callable | None, u_theta, r, theta=0.0, h: float = 1e-5, **kw):
    """Vorticity normal to the plane in polar coordinates, Eq. (3.23): ω_z = (1/r) ∂(r u_θ)/∂r − (1/r) ∂u_r/∂θ.

    Book: §3.5, Eq. (3.23) (from Appendix B; derived in our D17 from the circulation round a small polar sector,
    Exercise 3.20). For (3.22) it gives 2ω₀, for (3.25) zero (r > 0).

    Parameters
    ----------
    u_r, u_theta : callables f(r) or f(r, θ) [m/s] (``u_r=None`` means u_r = 0); ``u_theta`` may also be a kind
        name of :func:`vortex_profile` with keywords ``Gamma``, ``sigma`` (``**kw``)
    r : radius [m] (> h);  theta : angle [rad];  h : step [m] for r and [rad] for θ (explicit 2nd-order central)

    Returns
    -------
    ω_z [1/s]

    Validation: V1 equals the Cartesian curl of the same field (non-axisymmetric u_r, u_θ); V3 order 2;
    V2 ``polar_vorticity_z_sym`` for (3.22) → 2ω₀, (3.25) → 0. Label: analytic, converged, symbolic.
    """
    ur, ut = _polar_fn(u_r), _polar_fn(_profile_fn(u_theta, **kw))
    r_, th = _F(r), _F(theta)
    d_rut = ((r_ + h) * _F(ut(r_ + h, th)) - (r_ - h) * _F(ut(r_ - h, th))) / (2.0 * h)  # ∂(r u_θ)/∂r
    d_ur = (_F(ur(r_, th + h)) - _F(ur(r_, th - h))) / (2.0 * h)  # ∂u_r/∂θ
    return _S(d_rut / r_ - d_ur / r_)  # Eq. (3.23)


def polar_vorticity_z_sym(u_r_expr, u_theta_expr, r: sp.Symbol, theta: sp.Symbol):
    """Symbolic Eq. (3.23): ω_z = (1/r) ∂(r u_θ)/∂r − (1/r) ∂u_r/∂θ (simplified). Label: symbolic.
    Book: §3.5, Eq. (3.23).
    """
    return sp.simplify(sp.diff(r * sp.sympify(u_theta_expr), r) / r - sp.diff(sp.sympify(u_r_expr), theta) / r)


def vortex_velocity_field(profile, **params) -> Callable:
    """Cartesian plane velocity field u(x, t) = u_θ(r) e_θ = u_θ(r)(−y/r, x/r) of an axisymmetric vortex.

    ``profile`` is a kind of :func:`vortex_profile` ("solid", "line", "rankine", "gaussian"; params Gamma, sigma) or a
    callable u_θ(r). The field follows the kinematics convention (x of shape (2,) or (2, N), t ignored — steady); u = 0
    at r = 0 except for the line vortex (NaN). Book: §3.5, Eqs. (3.22), (3.25), (3.28), (3.29).

    Validation: V1 solid body equals ch02 ``solid_body_rotation_field``; its numerical curl is 2ω₀. Label: analytic.
    """
    if callable(profile):
        prof = profile
    else:
        kind = str(profile)
        prof = lambda r: vortex_profile(kind, r, **params)[0]  # noqa: E731

    def u(x, t=0.0):
        x_ = _F(x)
        r = np.hypot(x_[0], x_[1])
        ut = _F(prof(r))
        with np.errstate(divide="ignore", invalid="ignore"):
            f = np.where(r > 0.0, ut / np.where(r > 0.0, r, 1.0), np.where(np.isfinite(ut), 0.0, np.nan))
        return np.stack([-f * x_[1], f * x_[0]])  # u_θ e_θ, e_θ = (−y, x)/r
    return u


def circulation_circle(u_theta, r, center=(0.0, 0.0), n: int = 512, **kw) -> float:
    """Circulation Γ = ∮ u·ds around the circle of radius r about ``center``, Eqs. (3.18), (3.24), (3.26).

    Book: §3.5, Eq. (3.24): Γ = ∫₀^{2π} u_θ r dθ = 2π r u_θ for a circle centred on an axisymmetric vortex
    (2πr²ω₀ for solid-body rotation, 2πB for the line vortex (3.26)). For a centred circle and a profile u_θ(r) this
    formula is used; for an off-centre circle (Exercise 3.23: "any circuit") the loop integral of the Cartesian field
    is computed with ``core.integral_theorems.circulation`` (midpoint rule, n points — spectrally accurate for smooth
    periodic integrands).

    Parameters
    ----------
    u_theta : profile callable u_θ(r) [m/s] or a kind name of :func:`vortex_profile` with keywords ``Gamma``,
        ``sigma`` (``**kw``; defaults Γ = 2π m²/s, σ = 1 m)
    r : circle radius [m];  center : (x, y) [m];  n : loop samples (off-centre only)

    Returns
    -------
    Γ [m²/s]

    Validation: V1 solid body 2πr²ω₀ for centred and off-centre circles; line vortex 2πB for every r enclosing the
    origin and 0 for circles that exclude it. Label: analytic.
    """
    from .integral_theorems import circulation, planar_loop
    from .kinematics import as_coord_field

    prof = _profile_fn(u_theta, **kw)
    c = _F(center)
    if np.allclose(c, 0.0):
        return float(2.0 * np.pi * float(r) * float(_F(prof(float(r)))))  # Eq. (3.24): Γ = 2π r u_θ
    loop = planar_loop(c, radius=float(r), n=n)
    return circulation(as_coord_field(vortex_velocity_field(prof)), loop)  # Eq. (3.18): ∮ u·ds


def mean_vorticity_in_disc(u_theta, r, **kw):
    """Mean vorticity in the disc of radius r about the axis: Γ(r)/(πr²) = 2u_θ(r)/r, Eq. (3.27).

    Book: §3.5, Eq. (3.27): [ω_z]_{r→0} = lim (1/A)∫ω_z dA = lim (1/πr²)∮u·ds = lim 2B/r² for the line vortex —
    infinite at the axis with a finite area integral (a δ-function core). For solid-body rotation it is 2ω₀ at every r.
    Parameters: u_theta callable u_θ(r) [m/s] or a kind name (keywords ``Gamma``, ``sigma``), r [m] (> 0).
    Returns [1/s].
    Validation: V1 line vortex 2B/r² (log–log slope −2); solid body 2ω₀. Label: analytic.
    """
    r_ = _F(r)
    return _S(2.0 * _F(_profile_fn(u_theta, **kw)(r_)) / r_)  # (1/πr²)·2πr u_θ
