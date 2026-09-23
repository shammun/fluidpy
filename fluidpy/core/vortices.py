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
           "mean_vorticity_in_disc", "vortex_velocity_field",
           # Ch. 5 additions (exact 3-D vortex solutions and vorticity test fields)
           "burgers_vortex", "burgers_vortex_field", "burgers_core_radius", "burgers_pressure",
           "hill_spherical_vortex", "hill_stream_function", "hill_spherical_vortex_field", "hill_vortex_field",
           "hill_translation_speed", "lamb_oseen_field", "tube_core_radius", "gaussian_tube_vorticity",
           "gaussian_tube_field", "gaussian_tube_section", "broken_tube_field", "vortex_ring_vorticity",
           "stretched_gaussian_vortex_field", "abc_flow_field"]

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


# ======================================================================================================================
# Ch. 5 additions: exact 3-D vortex solutions and vorticity test fields
# ======================================================================================================================
def _cyl(x):
    """Cartesian points x (3,) or (3, N) → (x, y, z, R) with R = √(x² + y²)."""
    x_ = _F(x)
    return x_[0], x_[1], x_[2], np.hypot(x_[0], x_[1])


def burgers_core_radius(alpha: float, nu: float) -> float:
    """Core radius √(4ν/α) of Burgers' vortex (the e-folding radius of ω_z) [m].

    Book: Exercise 5.12 (§5.4 test field; our D18: stretching αω_z balances radial diffusion inside this radius).
    Parameters: alpha α [1/s] (> 0, axial stretching rate), nu ν [m²/s]. Label: analytic.
    """
    return float(np.sqrt(4.0 * float(nu) / float(alpha)))  # σ_B = √(4ν/α)


def burgers_vortex(R, z, Gamma: float, alpha: float, nu: float):
    """Burgers' vortex: axisymmetric strain plus a swirl whose vorticity is held steady by stretching against diffusion.

    Book: §5.4 / Exercise 5.12 (rendered page p221): u_R = −αR/2, u_φ = (Γ/2πR)[1 − exp(−αR²/4ν)], u_z = αz; the
    vorticity ω_z = (1/R) d(R u_φ)/dR = (αΓ/4πν) exp(−αR²/4ν) is our evaluation (part c). Steady because the
    stretching term αω_z of (5.13) balances the diffusion ν∇²ω_z and the inward advection (D18).

    Parameters
    ----------
    R : cylindrical radius [m] (≥ 0);  z : axial coordinate [m];  Gamma : circulation Γ [m²/s] (as R → ∞)
    alpha : stretching rate α [1/s] (> 0: elongation along z);  nu : kinematic viscosity ν [m²/s] (> 0)

    Returns
    -------
    (u_R, u_phi, u_z, omega_z) [m/s, m/s, m/s, 1/s]; u_φ(0) = 0 (cancellation-free via expm1).

    Assumptions: incompressible, constant ρ and ν, α > 0 for a steady concentrated core.
    Validation: V2 sympy NS residual 0 with :func:`burgers_pressure` (verifier); V1 ω_z by (3.23); V1 cross-check of
    the form against Wikipedia "Burgers vortex" (their α = book α/2). Label: symbolic, analytic.
    """
    R_, z_ = _F(R), _F(z)
    a = float(alpha) / (4.0 * float(nu))
    q = a * R_ ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        uphi = np.where(R_ > 0.0, float(Gamma) / (2.0 * np.pi) * (-np.expm1(-q)) / np.where(R_ > 0.0, R_, 1.0), 0.0)
    uR = -0.5 * float(alpha) * R_  # u_R = −αR/2  (Exercise 5.12)
    uz = float(alpha) * z_ + 0.0 * R_  # u_z = αz
    wz = float(alpha) * float(Gamma) / (4.0 * np.pi * float(nu)) * np.exp(-q)  # ω_z = (αΓ/4πν) e^{−αR²/4ν}
    return _S(uR), _S(uphi), _S(uz), _S(wz)


def burgers_pressure(R, z, Gamma: float, alpha: float, nu: float, rho: float = 1000.0, p0: float = 0.0):
    """Pressure of Burgers' vortex, p(R, z) with p(0, 0) = p₀ (Exercise 5.12 part b; our derivation D31).

    From the radial and axial Navier–Stokes components with u of :func:`burgers_vortex` (the viscous terms of u_R and
    u_z vanish because both are linear): ∂p/∂R = ρ(u_φ²/R − α²R/4), ∂p/∂z = −ρα²z, so
    p = p₀ − ρα²(R²/8 + z²/2) + ρ∫₀^R u_φ²/R′ dR′, and with s = αR²/4ν the swirl integral is closed-form:
    ∫₀^R u_φ²/R′ dR′ = (Γ/2π)²(α/8ν)[2 ln 2 − (1 − e^{−s})²/s − 2(E₁(s) − E₁(2s))] (E₁ = exponential integral).

    Parameters: as :func:`burgers_vortex`, plus rho [kg/m³], p0 [Pa]. Returns p [Pa]. Label: analytic.
    Book: Exercise 5.12 (the book asks for p; it does not print it).
    """
    from scipy.special import exp1

    R_, z_ = _F(R), _F(z)
    a = float(alpha) / (4.0 * float(nu))
    s = a * R_ ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        small = s < 1e-6
        ss = np.where(small, 1.0, s)
        F = 2.0 * np.log(2.0) - (np.expm1(-ss)) ** 2 / ss - 2.0 * (exp1(ss) - exp1(2.0 * ss))
        F = np.where(small, s - 0.0, F)  # series: F(s) = s + O(s²) near the axis
    swirl = (float(Gamma) / (2.0 * np.pi)) ** 2 * 0.5 * a * F  # ∫₀^R u_φ²/R′ dR′ = (Γ/2π)² (a/2) F(s)
    return _S(float(p0) - float(rho) * float(alpha) ** 2 * (R_ ** 2 / 8.0 + z_ ** 2 / 2.0) + float(rho) * swirl)


def burgers_vortex_field(Gamma: float, alpha: float, nu: float) -> Callable:
    """Cartesian velocity field u(x, t) (3-D, steady) of Burgers' vortex, for stencil and trajectory tools.

    u = (u_R x − u_φ y, u_R y + u_φ x)/R, u_z = αz (Exercise 5.12). x has shape (3,) or (3, N); t is ignored.
    Book: Exercise 5.12 (§5.4). Label: analytic.
    """
    a = float(alpha) / (4.0 * float(nu))
    G, al = float(Gamma), float(alpha)

    def u(x, t=0.0):
        X, Y, Z, R = _cyl(x)
        with np.errstate(divide="ignore", invalid="ignore"):
            f = np.where(R > 0.0, G / (2.0 * np.pi) * (-np.expm1(-a * R ** 2)) / np.where(R > 0.0, R, 1.0) ** 2,
                         G * a / (2.0 * np.pi))  # u_φ/R (finite on the axis: Γa/2π)
        return np.stack([-0.5 * al * X - f * Y, -0.5 * al * Y + f * X, al * Z])
    return u


def hill_translation_speed(A: float, a: float) -> float:
    """Speed U = 2Aa²/15 at which Hill's spherical vortex moves through fluid at rest (along +z for A > 0).

    Our matching of the exterior potential flow past the sphere to the interior of Exercise 5.11 (tangential velocity
    continuous at the sphere: −3U/2 = −Aa²/5); cross-check: Wikipedia "Hill's spherical vortex" ω = (15U/2a²) r sin θ.
    Units: A [1/(m s)], a [m] → U [m/s]. Book: Exercise 5.11. Label: analytic.
    """
    return 2.0 * float(A) * float(a) ** 2 / 15.0


def _hill_parts(R, z, A, a, outside):
    R_, z_ = _F(R), _F(z)
    A_, a_ = float(A), float(a)
    r2 = R_ ** 2 + z_ ** 2
    ins = r2 <= a_ ** 2
    psi_in = A_ / 10.0 * R_ ** 2 * (a_ ** 2 - R_ ** 2 - z_ ** 2)  # Exercise 5.11: ψ = (Aa⁴/10)(R²/a²)(1 − R²/a² − z²/a²)
    uR_in = A_ * R_ * z_ / 5.0  # u_R = −R⁻¹ ∂ψ/∂z
    uz_in = A_ / 5.0 * (a_ ** 2 - 2.0 * R_ ** 2 - z_ ** 2)  # u_z = R⁻¹ ∂ψ/∂R
    w_in = A_ * R_  # ω = AR e_φ
    if outside:
        U = hill_translation_speed(A_, a_)
        r = np.sqrt(np.where(ins, a_ ** 2, r2))
        k = a_ ** 3 / r ** 3
        psi_out = -0.5 * U * R_ ** 2 * (1.0 - k)  # uniform stream −U e_z past a sphere
        uR_out = -1.5 * U * a_ ** 3 * R_ * z_ / r ** 5  # u_R = −R⁻¹∂ψ/∂z
        uz_out = -U * (1.0 - k) - 1.5 * U * a_ ** 3 * R_ ** 2 / r ** 5  # u_z = R⁻¹∂ψ/∂R
        w_out = 0.0 * R_
    else:
        psi_out = uR_out = uz_out = w_out = np.full(np.shape(r2), np.nan)
    pick = lambda i, o: _S(np.where(ins, i, o))  # noqa: E731
    return pick(uR_in, uR_out), pick(uz_in, uz_out), pick(w_in, w_out), pick(psi_in, psi_out)


def hill_spherical_vortex(R, z, A: float, a: float, outside: bool = True):
    """Hill's spherical vortex in the frame moving with it (Exercise 5.11): (u_R, u_z, ω_φ).

    Book: Exercise 5.11 (rendered page p221): ψ(R, z) = (Aa⁴/10)(R²/a²)(1 − R²/a² − z²/a²) for R² + z² ≤ a²,
    u_R = −R⁻¹∂ψ/∂z = ARz/5, u_z = R⁻¹∂ψ/∂R = (A/5)(a² − 2R² − z²), ω = AR e_φ (ω_φ = ∂u_R/∂z − ∂u_z/∂R). A steady
    inviscid solution of (5.13) in which the ring-shaped vortex lines are stretched as they move outward (u_R/R > 0)
    exactly as fast as advection carries ω_φ = AR. Outside (our addition, ``outside=True``): the uniform stream −U e_z
    (U = 2Aa²/15) past the sphere, ψ = −(U/2)R²(1 − a³/r³); ψ = 0 on the sphere from both sides and the tangential
    velocity is continuous. ``outside=False`` returns NaN outside the sphere.

    Parameters: R cylindrical radius [m] (≥ 0); z axial coordinate [m] (centre at 0); A [1/(m s)]; a sphere radius [m].
    Returns (u_R, u_z, omega_phi) [m/s, m/s, 1/s] (ψ: :func:`hill_stream_function`).
    Assumptions: inviscid, steady in the co-moving frame, axisymmetric without swirl.
    Validation: V2 sympy: ∇·u = 0, ω = AR e_φ, (5.13) with ν = 0 satisfied; V1 ψ continuous at r = a; U = 2Aa²/15 vs
    the Wikipedia form (V1 cross-check). Label: symbolic, analytic.
    """
    uR, uz, w, _ = _hill_parts(R, z, A, a, outside)
    return uR, uz, w


def hill_stream_function(R, z, A: float, a: float, outside: bool = True):
    """Stokes stream function ψ(R, z) [m³/s] of :func:`hill_spherical_vortex` (inside: Exercise 5.11; outside: the
    potential flow past the sphere in the co-moving frame). Book: Exercise 5.11. Label: analytic."""
    return _hill_parts(R, z, A, a, outside)[3]


def hill_spherical_vortex_field(A: float, a: float, outside: bool = True) -> Callable:
    """Cartesian velocity u(x, t) (3-D, steady, co-moving frame) of :func:`hill_spherical_vortex`.
    Book: Exercise 5.11. Label: analytic."""
    def u(x, t=0.0):
        X, Y, Z, R = _cyl(x)
        uR, uz, _ = hill_spherical_vortex(R, Z, A, a, outside)
        uR, uz = _F(uR), _F(uz)
        with np.errstate(divide="ignore", invalid="ignore"):
            c = np.where(R > 0.0, X / np.where(R > 0.0, R, 1.0), 0.0)
            s = np.where(R > 0.0, Y / np.where(R > 0.0, R, 1.0), 0.0)
        return np.stack([uR * c, uR * s, uz])
    return u


hill_vortex_field = hill_spherical_vortex_field  # alias


def lamb_oseen_field(Gamma: float, nu: float, t0: float = 0.0) -> Callable:
    """Plane velocity u(x, t) of the Lamb–Oseen vortex, u_θ = (Γ/2πr)(1 − e^{−r²/σ²}), σ² = 4ν(t + t₀) (x (2,) or
    (2, N)). The Gaussian vortex (3.29) with the core diffusing; the same field as
    ``core.navier_stokes.exact_solution("lamb_oseen")`` (velocity only). Parameters: Gamma [m²/s]; nu [m²/s]; t0 [s]
    (age at t = 0). Book: §3.5 (3.29), §5.2 (5.11) test field. Label: analytic."""
    G, nu_, t0_ = float(Gamma), float(nu), float(t0)

    def u(x, t=0.0):
        X = _F(x)
        s2 = 4.0 * nu_ * (float(t) + t0_)
        r2 = X[0] ** 2 + X[1] ** 2
        with np.errstate(divide="ignore", invalid="ignore"):
            f = np.where(r2 > 0.0, G / (2.0 * np.pi) * (-np.expm1(-r2 / s2)) / np.where(r2 > 0.0, r2, 1.0),
                         G / (2.0 * np.pi * s2))  # u_θ/r (finite on the axis)
        return np.stack([-f * X[1], f * X[0]])
    return u


def tube_core_radius(z, a0: float = 0.1, L: float = 1.0):
    """Core radius a(z) = a₀e^{−z/2L} of :func:`gaussian_tube_vorticity` [m] (a² = a₀²e^{−z/L}).
    Book: §5.1 (Fig. 5.1, a tube that narrows along its length). Label: analytic."""
    return _S(float(a0) * np.exp(-_F(z) / (2.0 * float(L))))


def gaussian_tube_vorticity(R, z, Gamma: float = 1.0, a0: float = 0.1, L: float = 1.0, twist: float = 0.0):
    """A **narrowing Gaussian vortex tube**, exactly solenoidal: (ω_R, ω_φ, ω_z) at (R, z) [1/s].

    Construction (ours, for (5.4) and E1): flux function χ(R, z) = (Γ/2π)(1 − e^{−R²/a(z)²}) with a(z)² = a₀²e^{−z/L};
    ω_z = (1/R)∂χ/∂R = (Γ/πa²)e^{−R²/a²}, ω_R = −(1/R)∂χ/∂z = −(Γ/2π)(R/(a²L))e^{−R²/a²}, so
    (1/R)∂(Rω_R)/∂R + ∂ω_z/∂z = 0 exactly; ω_φ = τRω_z (τ = ``twist`` [1/m]) makes the lines helices without changing
    ∇·ω. The tube surfaces are R/a(z) = const: the tube narrows downstream (z increasing), the vorticity grows as 1/a²
    and the flux through the tube through (R₀, 0) stays Γ(1 − e^{−R₀²/a₀²}) — Eq. (5.4).
    Parameters: R, z [m]; Gamma Γ [m²/s]; a0 a₀ [m]; L [m]; twist τ [1/m].
    Book: §5.1, Eq. (5.4), Fig. 5.1. Label: analytic.
    """
    R_, z_ = _F(R), _F(z)
    a2 = float(a0) ** 2 * np.exp(-z_ / float(L))
    e = np.exp(-R_ ** 2 / a2)
    wz = float(Gamma) / (np.pi * a2) * e
    wR = -float(Gamma) / (2.0 * np.pi) * R_ / (a2 * float(L)) * e
    return _S(wR), _S(float(twist) * R_ * wz), _S(wz)


def gaussian_tube_field(Gamma: float = 1.0, a0: float = 0.1, L: float = 1.0, twist: float = 0.0) -> Callable:
    """Cartesian vorticity callable ω(x, t) → (3, …) [1/s] of :func:`gaussian_tube_vorticity`.
    Book: §5.1, Eq. (5.4). Label: analytic."""
    G, a0_, L_, tw = float(Gamma), float(a0), float(L), float(twist)

    def omega(x, t=0.0):
        X, Y, Z, R = _cyl(x)
        a2 = a0_ ** 2 * np.exp(-Z / L_)
        e = np.exp(-R ** 2 / a2)
        wz = G / (np.pi * a2) * e
        c1 = -G / (2.0 * np.pi) / (a2 * L_) * e  # ω_R / R
        c2 = tw * wz  # ω_φ / R
        return np.stack([c1 * X - c2 * Y, c1 * Y + c2 * X, wz])
    return omega


def gaussian_tube_section(z, R0: float = 0.1, Gamma: float = 1.0, a0: float = 0.1, L: float = 1.0) -> dict:
    """Cross-section at height z of the tube of :func:`gaussian_tube_vorticity` that passes through (R₀, 0):
    radius R₀e^{−z/2L}, area πR₀²e^{−z/L}, flux Γ(1 − e^{−R₀²/a₀²}) (the same at every z — Eq. (5.4)), mean vorticity
    flux/area, peak vorticity Γ/(πa(z)²) on the axis.
    Returns dict(radius [m], area [m²], flux [m²/s], mean_omega [1/s], peak_omega [1/s]); scalar-callable.
    Book: §5.1, Eq. (5.4) (tube strength independent of where it is measured). Label: analytic."""
    zf = float(z)
    rad = float(R0) * np.exp(-zf / (2.0 * float(L)))
    area = np.pi * rad ** 2
    flux = float(Gamma) * (-np.expm1(-(float(R0) / float(a0)) ** 2))
    a2 = float(a0) ** 2 * np.exp(-zf / float(L))
    return dict(radius=float(rad), area=float(area), flux=float(flux), mean_omega=float(flux / area),
                peak_omega=float(Gamma / (np.pi * a2)))


def broken_tube_field(Gamma: float = 1.0, a0: float = 0.1, L: float = 1.0) -> Callable:
    """⚠️ **Not a vorticity field** (∇·ω = ∂ω_z/∂z ≠ 0): ω = (Γ/πa₀²)e^{−R²/a₀²}e^{−z/L} e_z — a "tube" whose flux
    fades along its length, kept only to show Gauss' budget (5.4) failing (E1's "broken" preset). Returns ω(x, t).
    Book: §5.1, Eq. (5.4) (counterexample). Label: analytic."""
    G, a0_, L_ = float(Gamma), float(a0), float(L)

    def omega(x, t=0.0):
        X, Y, Z, R = _cyl(x)
        wz = G / (np.pi * a0_ ** 2) * np.exp(-R ** 2 / a0_ ** 2) * np.exp(-Z / L_)
        return np.stack([0.0 * wz, 0.0 * wz, wz])
    return omega


def vortex_ring_vorticity(R, z, Gamma: float = 1.0, R0: float = 0.5, a: float = 0.05):
    """Azimuthal vorticity of a Gaussian-core vortex ring of radius R₀ in the plane z = 0,
    ω_φ = (Γ/πa²)exp(−((R − R₀)² + z²)/a²) [1/s] (a ≪ R₀; the closed tube of E1). Its flux through any meridional
    half-plane is Γ (to e^{−R₀²/a²}). Book: §5.1 (tubes close on themselves), §5.7 (rings). Label: analytic."""
    return _S(float(Gamma) / (np.pi * float(a) ** 2) * np.exp(-((_F(R) - float(R0)) ** 2 + _F(z) ** 2) / float(a) ** 2))


def stretched_gaussian_vortex_field(Gamma: float = 1.0, sigma0: float = 1.0, alpha: float = 1.0,
                                    nu: float = 0.0) -> Callable:
    """Axisymmetric strain u_R = −αR/2, u_z = αz plus a Gaussian swirl whose core evolves exactly:
    σ²(t) = σ₀²e^{−αt} + (4ν/α)(1 − e^{−αt}), u_φ = (Γ/2πR)(1 − e^{−R²/σ²(t)}), ω_z = (Γ/πσ²)e^{−R²/σ²}.

    An exact unsteady Navier–Stokes solution (Euler when ν = 0) that contains the chapter's vortices as limits:
    ν = 0 — the inviscid tube shrinks, σ = σ₀e^{−αt/2}, and its peak vorticity grows as e^{αt} (stretching, (5.32));
    t → ∞ — Burgers' vortex (Exercise 5.12), σ² = 4ν/α; α → 0 — Lamb–Oseen, σ² = σ₀² + 4νt. Along a particle path
    the circulation Γ(R, t) = Γ(1 − e^{−R²/σ²}) of the material circle is conserved when ν = 0 (Kelvin (5.8)).

    Parameters: Gamma Γ [m²/s]; sigma0 σ₀ [m]; alpha α [1/s] (≠ 0; use a tiny value for the Lamb–Oseen limit);
    nu ν [m²/s]. Returns u(x, t) → (3, …) [m/s]. Book: §5.4–5.6 test field (our construction). Label: analytic.
    """
    G, s0, al, nu_ = float(Gamma), float(sigma0), float(alpha), float(nu)

    def sigma2(t):
        e = np.exp(-al * t)
        return s0 ** 2 * e + (4.0 * nu_ / al) * (1.0 - e) if al != 0.0 else s0 ** 2 + 4.0 * nu_ * t

    def u(x, t=0.0):
        X, Y, Z, R = _cyl(x)
        s2 = sigma2(float(t))
        with np.errstate(divide="ignore", invalid="ignore"):
            f = np.where(R > 0.0, G / (2.0 * np.pi) * (-np.expm1(-R ** 2 / s2)) / np.where(R > 0.0, R, 1.0) ** 2,
                         G / (2.0 * np.pi * s2))
        return np.stack([-0.5 * al * X - f * Y, -0.5 * al * Y + f * X, al * Z])
    u.sigma2 = sigma2
    return u


def abc_flow_field(A: float = 1.0, B: float = 1.0, C: float = 1.0) -> Callable:
    """Arnold–Beltrami–Childress flow u = (A sin z + C cos y, B sin x + A cos z, C sin y + B cos x): a steady Euler
    solution with ω = ∇ × u = u (Beltrami), so vortex lines are streamlines — a test field for (5.3), (5.13) with
    ν = 0 and Helmholtz's first theorem. Units: A, B, C [m/s], x [m] with unit wavenumber. Label: analytic.
    Book: §5.1, §5.3–5.4 (test field; not in the book)."""
    def u(x, t=0.0):
        X, Y, Z = _F(x)[0], _F(x)[1], _F(x)[2]
        return np.stack([A * np.sin(Z) + C * np.cos(Y), B * np.sin(X) + A * np.cos(Z), C * np.sin(Y) + B * np.cos(X)])
    return u
