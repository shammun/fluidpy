"""Lifting-line theory for wings of finite span: planforms, the trailing vortex sheet and its downwash, the sine-series
solution of Prandtl's equation, lift and induced drag, the elliptic loading, and horseshoe vortices.

Book: Kundu, Cohen & Dowling 5e, Ch. 14 §14.2 and §14.6, Eqs. (14.1)–(14.2), (14.13)–(14.27) (read from the rendered
pages chapters/pages/ch14/p719, p738–p742) and Exercises 14.10–14.18 (results only; no exercise text or numbers).

Conventions (``knowledge/notation.md``)
---------------------------------------
* **y runs along the span** (−s/2 port tip … +s/2 starboard tip), x downstream, z up.  The lifting line is on the
  y-axis; the stream U is along +x.
* **w is the downwash, positive downward** (the book's sign in (14.13)–(14.24)); velocity components returned as
  (ux, uy, uz) are in the axes above, so downwash = −uz.
* Γ(y) is the bound circulation with the book's sense (lift ρUΓ per unit span upward for Γ > 0).  The filament shed
  per unit span has strength −dΓ/dy about the +x direction (counterclockwise seen from behind on the starboard half).
* Span angle γ of (14.18): y = −(s/2) cos γ, γ = 0 at the port tip and π at the starboard tip.  ``Gn`` is the array
  (Γ₁, Γ₂, …) of the sine series Γ = Σ Γ_n sin nγ, so ``Gn[n − 1]`` multiplies sin nγ.
* **K is the section lift-curve slope per radian** and is a required argument everywhere (a default of 2π would
  silently mean "Zhukhovsky section"; trap T12).  β(y) is minus the zero-lift angle of the section [rad].
* Span efficiency is ``e_span`` (never ``e``, which is the thickness parameter of ``core.airfoils``).
* Slip #1: the first two members of (14.19) are printed with the denominator cos γ₁ − cos γ; it must be
  cos γ − cos γ₁.  The last member and everything after it are right.  ``downwash_series(printed_middle=True)``
  evaluates the printed middle member (it has the wrong sign).
* Scalars in, floats out.
"""
from __future__ import annotations

from typing import Callable

import numpy as np
from scipy.integrate import quad

from ._util import as_scalar_if_0d
from .biot_savart import segment_induced_velocity

__all__ = [
    "planform", "aspect_ratio", "mean_chord", "span_angle", "span_position", "circulation_series",
    "trailing_sheet_strength", "elliptic_circulation", "glauert_integral", "downwash", "downwash_series",
    "elliptic_downwash", "induced_angle", "effective_angle", "strip_forces", "lift_and_induced_drag",
    "lift_from_series", "induced_drag_from_series", "span_efficiency", "section_lift_coefficient",
    "lifting_line_modes", "lifting_line_matrix", "lifting_line_rhs", "lifting_line_residual", "solve_lifting_line",
    "elliptic_induced_drag", "induced_drag_coefficient", "finite_wing_lift_slope", "elliptic_wing",
    "twist_for_elliptic_loading", "elliptic_twist_rectangular", "horseshoe_velocity", "horseshoe_downwash",
    "horseshoe_wing", "downwash_along_x", "ground_effect_factor", "staircase_from_loading", "staircase_downwash",
    "rolled_up_pair", "wake_kinetic_energy_per_length",
]

_S = as_scalar_if_0d
_FOUR_PI = 4.0 * np.pi


def _F(a):
    return np.asarray(a, dtype=float)


def _Gn(Gn):
    g = np.atleast_1d(_F(Gn))
    return g, np.arange(1, g.size + 1, dtype=float)


# ======================================================================================================================
# planform (§14.2)
# ======================================================================================================================
def aspect_ratio(s, A):
    """Aspect ratio of a planform, Λ ≡ s²/A = s/c̄ (14.1) [–].

    Book: §14.2, Eq. (14.1).   Parameters: s wingspan [m]; A planform area [m²].   Returns Λ [–].
    Assumptions: none (definition).   Validation: V1 rectangle Λ = s/c; V2 both forms of (14.1) agree with (14.2).
    Label: analytic.
    """
    return _S(_F(s) ** 2 / _F(A))  # Eq. (14.1)


def mean_chord(A, s):
    """Average chord length, c̄ = A/s (14.2) [m].

    Book: §14.2, Eq. (14.2).   Parameters: A planform area [m²]; s wingspan [m].   Returns c̄ [m].
    Assumptions: none (definition).   Validation: V1 rectangle c̄ = c.  Label: analytic.
    """
    return _S(_F(A) / _F(s))  # Eq. (14.2)


def planform(kind: str, s: float, c_root: float, taper: float = 1.0, n: int = 201) -> dict:
    """Chord distribution c(y), area, aspect ratio and mean chord of a simple straight planform.

    Book: §14.2, Fig. 14.2 and Eqs. (14.1)–(14.2); Fig. 14.24 (the elliptic planform).  A = ∫ c(y) dy (implied).
    Parameters
    ----------
    kind : "rectangular" (c = c_root), "tapered" (c = c_root[1 − (1 − taper)|2y/s|], taper = tip chord / root chord)
    or "elliptic" (c = c_root √(1 − (2y/s)²)).   s : span [m].   c_root : chord at mid-span [m].   taper [–].
    n : number of spanwise samples (uniform in y).
    Returns
    -------
    dict: ``y`` [m], ``c`` [m] (arrays), ``chord`` (callable c(y)), ``A`` [m²] (closed form), ``A_trapezoid`` [m²]
    (trapezoid rule on the samples), ``AR`` = Λ (14.1), ``c_mean`` (14.2) [m], ``kind``, ``s``, ``c_root``, ``taper``.
    Assumptions: unswept, symmetric about mid-span.
    Validation: V1 rectangle Λ = s/c; ellipse A = πc_root s/4, Λ = 4s/(πc_root); trapezoid → closed form as n grows
    (order 1.5 for the ellipse: square-root ends).  Label: analytic.
    """
    s, c0, tp = float(s), float(c_root), float(taper)
    if kind == "rectangular":
        chord, A = (lambda y: c0 + 0.0 * _F(y)), c0 * s
    elif kind == "tapered":
        chord, A = (lambda y: c0 * (1.0 - (1.0 - tp) * np.abs(2.0 * _F(y) / s))), 0.5 * c0 * (1.0 + tp) * s
    elif kind == "elliptic":
        chord, A = (lambda y: c0 * np.sqrt(np.clip(1.0 - (2.0 * _F(y) / s) ** 2, 0.0, None))), 0.25 * np.pi * c0 * s
    else:
        raise ValueError('kind must be "rectangular", "tapered" or "elliptic"')
    y = np.linspace(-0.5 * s, 0.5 * s, int(n))
    c = chord(y)
    return {"y": y, "c": c, "chord": chord, "A": float(A), "A_trapezoid": float(np.trapezoid(c, y)),
            "AR": float(aspect_ratio(s, A)), "c_mean": float(mean_chord(A, s)), "kind": kind, "s": s, "c_root": c0,
            "taper": tp if kind == "tapered" else (1.0 if kind == "rectangular" else 0.0)}


# ======================================================================================================================
# the sine series (14.18) and the trailing sheet
# ======================================================================================================================
def span_angle(y, s):
    """Span angle γ of a spanwise position, from y = −(s/2) cos γ [rad] (γ = 0 port tip, π starboard tip).

    Book: §14.6, unnumbered line before (14.18).   Parameters: y [m] (|y| ≤ s/2); s span [m].   Returns γ [rad].
    Assumptions: none.   Validation: V7 round trip with :func:`span_position`.  Label: analytic.
    """
    return _S(np.arccos(np.clip(-2.0 * _F(y) / _F(s), -1.0, 1.0)))


def span_position(gamma, s):
    """Spanwise position of a span angle, y = −(s/2) cos γ [m].

    Book: §14.6, unnumbered line before (14.18).   Parameters: gamma [rad]; s [m].   Returns y [m].
    Assumptions: none.   Validation: V7 round trip.  Label: analytic.
    """
    return _S(-0.5 * _F(s) * np.cos(_F(gamma)))


def circulation_series(Gn, gamma):
    """Bound circulation from its sine series, Γ = Σ Γ_n sin(nγ) (14.18) [m²/s].

    Book: §14.6, Eq. (14.18).   Parameters: Gn coefficients (Γ₁, Γ₂, …) [m²/s]; gamma span angle [rad].
    Returns Γ(γ) [m²/s] (shape of gamma).
    Assumptions: Γ vanishes at both tips (every term does).
    Validation: V1 one term = the ellipse (14.23); V7 zero at γ = 0, π.  Label: analytic.
    """
    g, n = _Gn(Gn)
    ga = _F(gamma)
    return _S(np.sum(g * np.sin(n * ga[..., None]), axis=-1))  # Eq. (14.18)


def trailing_sheet_strength(Gn, gamma, s):
    """Strength of the trailing vortex sheet per unit span, −dΓ/dy [m/s] (circulation about +x per metre of span).

    Book: §14.6, Fig. 14.21 and the text: "the magnitude of the trailing vortex filament of width dy is −(dΓ/dy)dy".
    With (14.18): −dΓ/dy = −Σ nΓ_n cos nγ / ((s/2) sin γ).
    Parameters: Gn [m²/s]; gamma [rad] (0 < γ < π); s [m].   Returns −dΓ/dy [m/s] (infinite at the tips in general).
    Assumptions: filaments leave straight back (assumption (iii) of the theory).
    Validation: V1 finite difference of :func:`circulation_series`; V7 odd in y for a symmetric loading.
    Label: analytic.
    """
    g, n = _Gn(Gn)
    ga = _F(gamma)
    with np.errstate(divide="ignore", invalid="ignore"):
        return _S(-np.sum(n * g * np.cos(n * ga[..., None]), axis=-1) / (0.5 * float(s) * np.sin(ga)))


def elliptic_circulation(y, s, Gamma1):
    """Elliptic distribution of circulation, Γ = Γ₁ sin γ = Γ₁ √(1 − (2y/s)²) (14.23) [m²/s].

    Book: §14.6, Eq. (14.23).   Parameters: y [m]; s [m]; Gamma1 mid-span circulation [m²/s].   Returns Γ(y).
    Assumptions: only the first coefficient of (14.18) is non-zero.
    Validation: V1 equals ``circulation_series([Gamma1], span_angle(y, s))``.  Label: analytic.
    """
    return _S(_F(Gamma1) * np.sqrt(np.clip(1.0 - (2.0 * _F(y) / _F(s)) ** 2, 0.0, None)))  # Eq. (14.23)


def glauert_integral(n: int, gamma1: float, method: str = "closed") -> float:
    """Glauert's integral, PV∫₀^π cos(nγ) dγ/(cos γ − cos γ₁) = π sin(nγ₁)/sin γ₁.

    Book: §14.6, Eq. (14.19) ("the final equality comes from evaluating the integral" — not shown).  Ours (standard):
    the principal value; note the order of the two cosines in the denominator (slip #1).
    ``method="closed"``: the right-hand side (limit n(±1)^{n+1} at the tips); ``method="quad"``: adaptive quadrature of
    the regular integrand [cos nγ − cos nγ₁]/(cos γ − cos γ₁) (allowed because the n = 0 integral is zero).
    Parameters: n integer ≥ 0; gamma1 [rad] in (0, π).   Returns the value [–].
    Assumptions: principal value.   Validation: V1 quad = closed for n = 0…6 (1e-8); recursion
    I_{n+1} + I_{n−1} = 2 cos γ₁ I_n.  Label: analytic.
    """
    g1 = float(gamma1)
    if method == "closed":
        sg = np.sin(g1)
        if abs(sg) < 1e-12:
            return float(np.pi * n * (1.0 if np.cos(g1) > 0 else (-1.0) ** (n + 1)))
        return float(np.pi * np.sin(n * g1) / sg)
    if method != "quad":
        raise ValueError('method must be "closed" or "quad"')
    c1, cn1 = np.cos(g1), np.cos(n * g1)

    def f(g):
        d = np.cos(g) - c1
        if abs(d) < 1e-14:
            return -n * np.sin(n * g1) / (-np.sin(g1)) if abs(np.sin(g1)) > 0 else 0.0
        return (np.cos(n * g) - cn1) / d

    val, _ = quad(f, 0.0, np.pi, points=[g1], limit=200, epsabs=1e-12, epsrel=1e-12)
    return float(val)


def _ddgamma(G: Callable, ga, h: float = 1e-3):
    """dG/dγ by fourth-order differences (one-sided within 2h of the tips, where G may have a kink when reflected)."""
    ga = _F(ga)
    cen = (-G(ga + 2 * h) + 8.0 * G(ga + h) - 8.0 * G(ga - h) + G(ga - 2 * h)) / (12.0 * h)
    fwd = (-25.0 * G(ga) + 48.0 * G(ga + h) - 36.0 * G(ga + 2 * h) + 16.0 * G(ga + 3 * h) - 3.0 * G(ga + 4 * h)) / (12.0 * h)
    bwd = (25.0 * G(ga) - 48.0 * G(ga - h) + 36.0 * G(ga - 2 * h) - 16.0 * G(ga - 3 * h) + 3.0 * G(ga - 4 * h)) / (12.0 * h)
    return np.where(ga < 2.5 * h, fwd, np.where(ga > np.pi - 2.5 * h, bwd, cen))


def downwash(y1, Gamma_fn: Callable, s: float, method: str = "quad", h: float = 1e-3):
    """Downwash on the lifting line from the whole trailing sheet, by direct quadrature of (14.13) [m/s].

    Book: §14.6, Eq. (14.13) w(y₁) = (1/4π)∫_{−s/2}^{s/2} (dΓ/dy) dy/(y₁ − y) — a principal value (the bound vortex
    induces nothing on itself).
    Method (ours): with y = −(s/2) cos γ, w = (1/2πs) PV∫₀^π g(γ) dγ/(cos γ − cos γ₁), g = dΓ/dγ; the singular part is
    removed with PV∫₀^π dγ/(cos γ − cos γ₁) = 0, leaving the regular integrand [g(γ) − g(γ₁)]/(cos γ − cos γ₁), which
    ``scipy.integrate.quad`` integrates with a break point at γ₁.  g by fourth-order differences of Γ(y(γ)) with the
    step ``h`` [rad].
    Parameters
    ----------
    y1 : station(s) [m], |y1| < s/2.   Gamma_fn : callable Γ(y) [m²/s], vectorised, zero at y = ±s/2.   s : span [m].
    method : "quad" (the only one; kept for the explainer contract).   h : difference step in γ [rad].
    Returns
    -------
    w(y₁) [m/s], positive downward.
    Assumptions: lifting-line assumptions (large aspect ratio, straight bound vortex, straight non-interacting
    trailing filaments parallel to the stream).
    Validation: V1 equals :func:`downwash_series` for 1-, 2-, 3-term loadings (1e-8); elliptic loading → Γ₁/2s (14.24)
    at every station.  Label: analytic, converged.
    """
    if method != "quad":
        raise ValueError('method must be "quad"')
    s = float(s)

    def G(ga):
        return _F(Gamma_fn(-0.5 * s * np.cos(ga)))

    def one(yy):
        g1 = float(np.arccos(np.clip(-2.0 * yy / s, -1.0, 1.0)))
        c1 = np.cos(g1)
        gp1 = float(_ddgamma(G, g1, h))

        def f(ga):
            d = np.cos(ga) - c1
            if abs(d) < 1e-9:
                return 0.0  # a removable point of measure zero (quad never needs it at the break point)
            return (float(_ddgamma(G, ga, h)) - gp1) / d

        val, _ = quad(f, 0.0, np.pi, points=[g1], limit=400, epsabs=1e-12, epsrel=1e-11)
        return val / (2.0 * np.pi * s)  # Eq. (14.13) in the span angle

    yv = _F(y1)
    out = np.array([one(float(q)) for q in np.atleast_1d(yv).ravel()]).reshape(np.atleast_1d(yv).shape)
    return _S(out if yv.ndim else out[0])


def downwash_series(Gn, gamma1, s, printed_middle: bool = False):
    """Downwash of the sine-series loading, w(y₁) = (1/2s) Σ nΓ_n sin(nγ₁)/sin γ₁ (14.19, last member) [m/s].

    Book: §14.6, Eq. (14.19).  The last member is coded.  Slip #1: the first two members are printed with the
    denominator cos γ₁ − cos γ; from (14.13) with y = −(s/2) cos γ it must be cos γ − cos γ₁, so the printed middle
    member equals *minus* the last one.  ``printed_middle=True`` evaluates that printed member by quadrature
    (:func:`glauert_integral` with the sign of the printed denominator) — a test must show it has the opposite sign.
    Parameters: Gn [m²/s]; gamma1 [rad]; s [m]; printed_middle (bool).   Returns w [m/s], positive downward.
    At the tips sin nγ/sin γ → n (γ = 0) and n(−1)^{n+1} (γ = π).
    Assumptions: lifting-line assumptions; series differentiable term by term.
    Validation: V1 equals :func:`downwash` (quadrature of (14.13)); one term → Γ₁/2s (14.24).  Label: analytic.
    """
    g, n = _Gn(Gn)
    ga = _F(gamma1)
    s = float(s)
    if printed_middle:
        def one(q):
            return sum(nn * gg * (-glauert_integral(int(nn), q, "quad")) for nn, gg in zip(n, g)) / (2.0 * np.pi * s)
        out = np.array([one(float(q)) for q in np.atleast_1d(ga).ravel()]).reshape(np.atleast_1d(ga).shape)
        return _S(out if ga.ndim else out[0])
    sg = np.sin(ga)[..., None]
    tip = np.abs(sg) < 1e-12
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.sin(n * ga[..., None]) / np.where(tip, 1.0, sg)
    lim = np.where(np.cos(ga)[..., None] > 0, n, n * (-1.0) ** (n + 1))
    ratio = np.where(tip, lim, ratio)
    return _S(np.sum(n * g * ratio, axis=-1) / (2.0 * s))  # Eq. (14.19), last member


def elliptic_downwash(Gamma1, s):
    """Downwash of the elliptic loading, w = Γ₁/2s (14.24) [m/s] — the same at every spanwise station.

    Book: §14.6, Eq. (14.24).   Parameters: Gamma1 [m²/s]; s [m].   Returns w [m/s] (downward).
    Assumptions: elliptic loading (14.23).   Validation: V1 quadrature of (14.13).  Label: analytic.
    """
    return _S(_F(Gamma1) / (2.0 * _F(s)))  # Eq. (14.24)


def induced_angle(w, U, small: bool = True):
    """Angle by which the downwash turns the local stream down, ε = tan⁻¹(w/U) ≃ w/U [rad].

    Book: §14.6, unnumbered line before (14.14) and Fig. 14.23.  Slip #2: the page prints "ε = tan(w/U)"; the figure
    shows tan ε = w/U.   Parameters: w downwash [m/s]; U stream speed [m/s]; small → the approximation w/U.
    Returns ε [rad].   Assumptions (small=True): w ≪ U.
    Validation: V7 the two forms agree to O((w/U)³).  Label: analytic.
    """
    r = _F(w) / _F(U)
    return _S(r if small else np.arctan(r))


def effective_angle(alpha, w, U):
    """Effective angle of attack of a spanwise station, α_e = α − ε = α − w/U (14.14) [rad].

    Book: §14.6, Eq. (14.14).   Parameters: alpha geometric angle of attack [rad]; w downwash [m/s]; U [m/s].
    Returns α_e [rad].   Assumptions: w ≪ U.   Validation: V7 w = 0 → α.  Label: analytic.
    """
    return _S(_F(alpha) - _F(w) / _F(U))  # Eq. (14.14)


# ======================================================================================================================
# forces
# ======================================================================================================================
def strip_forces(Gamma, w, U, rho, dy, exact: bool = False) -> dict:
    """Lift and induced drag of a strip dy of the wing (the Kutta–Zhukhovsky force normal to the local stream).

    Book: §14.6, Fig. 14.23 and the unnumbered lines dL = dL_e cos ε = ρU_eΓ dy cos ε ≃ ρUΓ dy,
    dD_i = dL_e sin ε = ρU_eΓ dy sin ε ≃ ρwΓ dy, with dL_e = ρU_eΓ dy normal to U_e.
    Because U_e cos ε = U and U_e sin ε = w exactly, the two results hold for any w once Γ is given; ``exact`` only
    changes the reported ε (tan⁻¹(w/U) instead of w/U).
    Parameters: Gamma [m²/s]; w downwash [m/s]; U [m/s]; rho [kg/m³]; dy strip width [m].
    Returns dict: ``dL``, ``dDi``, ``dLe`` [N]; ``eps`` [rad]; ``Ue`` [m/s].
    Assumptions: each strip is a two-dimensional section in the resultant stream.
    Validation: V1 dL² + dDi² = dLe²; V7 w = 0 → no induced drag.  Label: analytic.
    """
    G, w_, U_, r, d = _F(Gamma), _F(w), _F(U), _F(rho), _F(dy)
    Ue = np.hypot(U_, w_)
    return {"dL": _S(r * U_ * G * d), "dDi": _S(r * w_ * G * d), "dLe": _S(r * Ue * G * d),
            "eps": induced_angle(w_, U_, small=not exact), "Ue": _S(Ue)}


def lift_and_induced_drag(y, Gamma, w, rho, U) -> tuple:
    """Lift and induced drag of the wing by quadrature, L = ∫ρUΓ dy, D_i = ∫ρwΓ dy (14.15) [N].

    Book: §14.6, Eq. (14.15).   Parameters: y stations [m] (increasing, tip to tip); Gamma [m²/s] and w [m/s] at
    them; rho [kg/m³]; U [m/s].   Returns (L, D_i) [N] (trapezoid rule; use stations uniform in γ for accuracy).
    Assumptions: lifting-line assumptions; w ≪ U.
    Validation: V1 equals (14.21), (14.22) for a series loading.  Label: analytic, converged.
    """
    y, G, w_ = _F(y), _F(Gamma), _F(w)
    L = float(rho) * float(U) * np.trapezoid(G, y)  # Eq. (14.15)
    Di = float(rho) * np.trapezoid(w_ * G, y)  # Eq. (14.15)
    return float(L), float(Di)


def lift_from_series(Gn, s, rho, U):
    """Lift of the wing from the first coefficient only, L = (πs/4) ρUΓ₁ (14.21) [N].

    Book: §14.6, Eq. (14.21).   Parameters: Gn [m²/s]; s [m]; rho [kg/m³]; U [m/s].   Returns L [N].
    Assumptions: lifting-line assumptions.   Validation: V1 quadrature of (14.15).  Label: analytic.
    """
    g, _ = _Gn(Gn)
    return float(0.25 * np.pi * float(s) * float(rho) * float(U) * g[0])  # Eq. (14.21)


def induced_drag_from_series(Gn, rho):
    """Induced drag of the wing, D_i = (π/8) ρ Σ nΓ_n² (14.22) [N] — every coefficient adds a positive amount.

    Book: §14.6, Eq. (14.22) (upper limit N: written for the truncated solution).
    Parameters: Gn [m²/s]; rho [kg/m³].   Returns D_i [N].
    Assumptions: lifting-line assumptions; planar wake.
    Validation: V1 quadrature of (14.15); V4 equals the wake kinetic energy per unit length.  Label: analytic.
    """
    g, n = _Gn(Gn)
    return float(np.pi / 8.0 * float(rho) * np.sum(n * g ** 2))  # Eq. (14.22)


def span_efficiency(Gn) -> float:
    """Span efficiency e_span = Γ₁²/Σ nΓ_n² ≤ 1 (ours), so that C_Di = C_L²/(πΛ e_span).

    Book: §14.6, text after (14.22): the best lift-to-drag ratio has Γ_n = 0 for n > 1.  Ours: the ratio of the
    induced drag of the elliptic loading with the same lift to the actual one, from (14.21), (14.22).
    Parameters: Gn [m²/s].   Returns e_span [–] (1 only for the elliptic loading; NaN if Γ₁ = 0).
    Assumptions: lifting-line assumptions.   Validation: V7 ≤ 1 with equality for one term.  Label: analytic.
    """
    g, n = _Gn(Gn)
    den = np.sum(n * g ** 2)
    return float(g[0] ** 2 / den) if den > 0 else float("nan")


def section_lift_coefficient(alpha, w, U, beta=0.0, *, K):
    """Section lift coefficient on a finite wing, C_L = K[α − w(y)/U + β(y)] (14.16) [–].

    Book: §14.6, Eq. (14.16).   Parameters: alpha [rad]; w downwash [m/s]; U [m/s]; beta = −(zero-lift angle) [rad];
    K lift-curve slope [1/rad] (keyword, required).   Returns C_L [–].
    Assumptions: linear section law (no stall).   Validation: V7 w = 0, K = 2π → (14.12).  Label: analytic.
    """
    return _S(_F(K) * (_F(alpha) - _F(w) / _F(U) + _F(beta)))  # Eq. (14.16)


# ======================================================================================================================
# the lifting-line equation (14.17), (14.20)
# ======================================================================================================================
def _at(q, y, default=0.0):
    """A scalar, a callable of y, a planform dict, or None, evaluated at the stations y."""
    if q is None:
        return np.full_like(_F(y), default)
    if isinstance(q, dict):
        return _F(q["chord"](y)) + 0.0 * _F(y)
    if callable(q):
        return _F(q(y)) + 0.0 * _F(y)
    q = _F(q)
    if q.ndim == 0:
        return np.full_like(_F(y), float(q))
    if q.shape != _F(y).shape:
        raise ValueError("an array input must be given at the collocation stations (same length)")
    return q


def _area(c, s, m: int = 400) -> float:
    if isinstance(c, dict):
        return float(c["A"])
    if not callable(c):
        c = _F(c)
        if c.ndim == 0:
            return float(c) * float(s)
        raise ValueError("give the chord as a scalar, a callable c(y) or a planform dict (the area is needed)")
    xg, wg = np.polynomial.legendre.leggauss(int(m))
    ga = 0.5 * np.pi * (xg + 1.0)
    return float(0.5 * np.pi * np.sum(wg * _F(c(-0.5 * s * np.cos(ga))) * 0.5 * s * np.sin(ga)))


def lifting_line_modes(N: int, symmetric: bool = False) -> tuple:
    """Mode numbers n and collocation angles γ_j of the truncated lifting-line system.

    Book: §14.6, text after (14.20) ("truncating the sum after N terms, and then requiring its validity at N points").
    Ours: general — n = 1…N at γ_j = jπ/(N + 1); ``symmetric`` (loading even in y) — odd n = 1, 3, …, 2N − 1 at
    γ_j = jπ/(2N), j = 1…N (port half, the last station at mid-span).  The tips, where sin γ = 0, are never stations.
    Returns (n, gamma).   Label: analytic.
    """
    N = int(N)
    j = np.arange(1, N + 1, dtype=float)
    if symmetric:
        return 2.0 * j - 1.0, j * np.pi / (2.0 * N)
    return j, j * np.pi / (N + 1.0)


def lifting_line_matrix(c, s: float, K, N: int, *, symmetric: bool = False) -> tuple:
    """Matrix of the collocated lifting-line equation: A[j, k] = (1 + n_k K c_j/(4 s sin γ_j)) sin(n_k γ_j).

    Book: §14.6, Eq. (14.20) (K/2)Uc(α + β) = Σ (1 + nKc/(4s sin γ)) Γ_n sin(nγ) required at N stations.
    Parameters: c chord — scalar [m], callable c(y), planform dict, or array at the stations; s span [m]; K slope
    [1/rad] (scalar, callable K(y) or array at the stations); N number of terms; symmetric (see
    :func:`lifting_line_modes`).
    Returns (A, gamma): the (N, N) matrix [–] and the station angles [rad].
    Assumptions: lifting-line assumptions.   Validation: V1 entry by entry against two nested loops.  Label: analytic.
    """
    n, ga = lifting_line_modes(N, symmetric)
    y = -0.5 * float(s) * np.cos(ga)
    cj, Kj = _at(c, y), _at(K, y)
    A = (1.0 + n[None, :] * (Kj * cj)[:, None] / (4.0 * float(s) * np.sin(ga)[:, None])) * np.sin(n[None, :] * ga[:, None])
    return A, ga  # Eq. (14.20), left-hand side


def lifting_line_rhs(c, s: float, U: float, alpha: float, K, N: int, *, beta=0.0, twist=None,
                     symmetric: bool = False) -> np.ndarray:
    """Right-hand side of the collocated lifting-line equation, r_j = (K/2) U c_j (α + twist_j + β_j) [m²/s].

    Book: §14.6, Eq. (14.20), left member.   Parameters: as :func:`lifting_line_matrix`; U [m/s]; alpha [rad];
    beta (scalar / callable / array) [rad]; twist — geometric twist added to α (same forms) [rad].
    Returns r (N,) [m²/s].   Label: analytic.
    """
    _, ga = lifting_line_modes(N, symmetric)
    y = -0.5 * float(s) * np.cos(ga)
    return 0.5 * _at(K, y) * float(U) * _at(c, y) * (float(alpha) + _at(twist, y) + _at(beta, y))  # Eq. (14.20)


def lifting_line_residual(Gn, gamma, c, s: float, U: float, alpha: float, *, K, beta=0.0, twist=None):
    """Residual of Prandtl's lifting-line equation (14.17) at arbitrary span angles, for a series loading [m²/s].

    Book: §14.6, Eq. (14.17) Γ(y) = (K/2) U c(y)[α − w(y)/U + β(y)] with w from (14.13) (here its series form, the
    last member of (14.19)).  Returns Γ − (K/2)Uc[α + twist − w/U + β]: zero at the collocation stations of
    :func:`solve_lifting_line`, small between them once N is large enough.
    Parameters: Gn [m²/s]; gamma [rad]; c, K, beta, twist as in :func:`lifting_line_matrix` (callables of y or
    scalars); s [m]; U [m/s]; alpha [rad].   Returns the residual [m²/s] (shape of gamma).
    Assumptions: lifting-line assumptions.   Validation: V3 max residual between stations → 0 as N grows.
    Label: converged.
    """
    ga = _F(gamma)
    y = -0.5 * float(s) * np.cos(ga)
    G = _F(circulation_series(Gn, ga))
    w = _F(downwash_series(Gn, ga, s))
    return _S(G - 0.5 * _at(K, y) * float(U) * _at(c, y) * (float(alpha) + _at(twist, y) - w / float(U) + _at(beta, y)))


def solve_lifting_line(c, s: float, U: float, alpha: float, *, K, beta=0.0, twist=None, N: int = 20,
                       symmetric: bool = False, rho: float = 1.0) -> dict:
    """Solve the lifting-line equation (14.20) for the coefficients Γ_n by collocation, and report the wing's forces.

    Book: §14.6, Eqs. (14.17)–(14.22): N terms of (14.18), (14.20) required at N spanwise stations ⇒ N linear
    equations; then L = (πs/4)ρUΓ₁ (14.21) and D_i = (π/8)ρ Σ nΓ_n² (14.22).
    Method: dense solve of :func:`lifting_line_matrix` · Γ = :func:`lifting_line_rhs` (``numpy.linalg.solve``).
    The book leaves N and the stations open — both are arguments here (:func:`lifting_line_modes`).
    Parameters
    ----------
    c : chord — scalar [m], callable c(y) [m] or a :func:`planform` dict.   s : span [m].   U : speed [m/s].
    alpha : geometric angle of attack of the wing [rad].   K : section lift-curve slope [1/rad] (required; scalar or
    callable K(y)).   beta : −(zero-lift angle) [rad] (scalar or callable).   twist : local angle added to α [rad]
    (callable of y or scalar).   N : number of terms/stations (default 20).   symmetric : use odd modes only (inputs
    must be even in y).   rho : density [kg/m³] (default 1 — forces then per unit density).
    Returns
    -------
    dict: ``Gn`` (Γ₁ … Γ_nmax, zeros at the even places when symmetric) [m²/s]; ``n``; ``gamma``, ``y`` (stations);
    ``Gamma``, ``w`` at the stations [m²/s], [m/s]; ``L``, ``Di`` [N]; ``A`` [m²]; ``AR``; ``CL``, ``CDi`` [–];
    ``e_span`` (alias ``e``); ``slope`` = dC_L/dα [1/rad] (solved with unit α, no twist, β = 0); ``cond`` (condition
    number); ``matrix``, ``rhs``; ``Gamma1``.
    Assumptions: lifting-line assumptions (large aspect ratio, unswept, planar straight wake), linear sections.
    Validation: V1 elliptic planform: Γ_n = 0 for n > 1 and C_L = K(α + β)/(1 + K/πΛ); V3 rectangular wing converges
    with N; residual of (14.17) between stations → 0; V7 Λ → ∞: C_L → K(α + β).  Label: analytic, converged.
    """
    s, U = float(s), float(U)
    n, ga = lifting_line_modes(N, symmetric)
    M, _ = lifting_line_matrix(c, s, K, N, symmetric=symmetric)
    r = lifting_line_rhs(c, s, U, alpha, K, N, beta=beta, twist=twist, symmetric=symmetric)
    sol = np.linalg.solve(M, r)  # Eq. (14.20) at N stations
    r1 = lifting_line_rhs(c, s, U, 1.0, K, N, symmetric=symmetric)
    sol1 = np.linalg.solve(M, r1)
    nmax = int(n[-1])
    Gn = np.zeros(nmax)
    Gn[n.astype(int) - 1] = sol
    A = _area(c, s)
    AR = s ** 2 / A
    q = 0.5 * float(rho) * U ** 2 * A
    L = lift_from_series(Gn, s, rho, U)
    Di = induced_drag_from_series(Gn, rho)
    slope = 0.25 * np.pi * s * float(rho) * U * sol1[0] / q
    es = span_efficiency(Gn)
    return {"Gn": Gn, "n": np.arange(1, nmax + 1), "gamma": ga, "y": -0.5 * s * np.cos(ga),
            "Gamma": _F(circulation_series(Gn, ga)), "w": _F(downwash_series(Gn, ga, s)), "L": L, "Di": Di,
            "A": float(A), "AR": float(AR), "CL": float(L / q), "CDi": float(Di / q), "e_span": es, "e": es,
            "slope": float(slope), "cond": float(np.linalg.cond(M)), "matrix": M, "rhs": r, "Gamma1": float(Gn[0])}


# ======================================================================================================================
# the elliptic loading (14.23)–(14.27)
# ======================================================================================================================
def elliptic_induced_drag(L, rho, U, s):
    """Induced drag of the elliptic loading, D_i = (π/8)ρΓ₁² = 2L²/(πρU²s²) (14.25) [N].

    Book: §14.6, Eq. (14.25) (second form, with (14.21)).   Parameters: L lift [N]; rho [kg/m³]; U [m/s]; s span [m].
    Returns D_i [N].   Assumptions: elliptic loading; lifting-line assumptions.
    Validation: V1 both forms of (14.25) agree; V5 cited form L²/(½ρU²πs²).  Label: analytic.
    """
    return _S(2.0 * _F(L) ** 2 / (np.pi * _F(rho) * _F(U) ** 2 * _F(s) ** 2))  # Eq. (14.25)


def induced_drag_coefficient(CL, AR, e_span: float = 1.0):
    """Induced drag coefficient, C_Di = C_L²/(πΛ) (14.26) (divided by e_span for a non-elliptic loading, ours) [–].

    Book: §14.6, Eq. (14.26) C_Di = D_i/(½ρU²A) = C_L²/(π(s²/A)) = C_L²/(πΛ).
    Parameters: CL [–]; AR = Λ aspect ratio [–]; e_span span efficiency [–] (1 = elliptic loading).
    Returns C_Di [–].   Assumptions: lifting-line assumptions; A is the planform area.
    Validation: V1 from (14.25); V5 cited form C_L²/(π AR e).  Label: analytic.
    """
    return _S(_F(CL) ** 2 / (np.pi * _F(AR) * _F(e_span)))  # Eq. (14.26)


def finite_wing_lift_slope(K, AR):
    """Lift-curve slope of an untwisted elliptic wing, dC_L/dα = K/(1 + K/(πΛ)) [1/rad].

    Book: §14.6, Exercise 14.15 (C_L,3D = C_L,2D/(1 + 2/Λ) for K = 2π), from (14.17) with (14.24) and (14.27).
    Parameters: K section slope [1/rad]; AR = Λ [–].   Returns the wing's slope [1/rad].
    Assumptions: elliptic planform, uniform section, no twist; lifting-line assumptions.
    Validation: V1 ``solve_lifting_line`` on an elliptic planform; V7 Λ → ∞ → K.  Label: analytic.
    """
    return _S(_F(K) / (1.0 + _F(K) / (np.pi * _F(AR))))


def elliptic_wing(alpha, AR, *, K, beta=0.0) -> dict:
    """Untwisted elliptic wing: lift and induced-drag coefficients, induced angle.

    Book: §14.6, Eqs. (14.24), (14.26), (14.27) and Exercise 14.15.  C_L = K(α + β)/(1 + K/(πΛ)),
    w/U = Γ₁/(2sU) = C_L/(πΛ), C_Di = C_L²/(πΛ).
    Parameters: alpha [rad]; AR = Λ [–]; K [1/rad] (keyword, required); beta [rad].
    Returns dict: ``CL``, ``CDi``, ``slope`` [1/rad], ``w_over_U`` (= induced angle ε [rad]), ``alpha_e`` [rad].
    Assumptions: as :func:`finite_wing_lift_slope`.   Validation: V1 ``solve_lifting_line``.  Label: analytic.
    """
    sl = finite_wing_lift_slope(K, AR)
    CL = sl * (_F(alpha) + _F(beta))
    eps = CL / (np.pi * _F(AR))
    return {"CL": _S(CL), "CDi": induced_drag_coefficient(CL, AR), "slope": sl, "w_over_U": _S(eps),
            "alpha_e": _S(_F(alpha) - eps)}


def twist_for_elliptic_loading(c, y, s: float, *, K, CL_target: float, A: float | None = None) -> dict:
    """Local angle of attack that gives an elliptic loading on a wing of any planform.

    Book: §14.6, text after (14.27) and Exercise 14.14 (an elliptic loading can be had on a non-elliptic planform by
    twist).  From (14.17) with (14.23), (14.24): α(y) + β(y) = 2Γ₁ sin γ/(K U c(y)) + Γ₁/(2sU), and
    Γ₁/U = 2C_L A/(πs) from (14.21).
    Parameters: c chord — callable c(y), scalar or array at ``y`` [m]; y stations [m]; s [m]; K [1/rad] (keyword);
    CL_target wing lift coefficient [–]; A planform area [m²] (None → from c: scalar or callable).
    Returns dict: ``alpha_local`` = α(y) + β(y) [rad] at ``y``; ``alpha_i`` = Γ₁/(2sU) [rad] (uniform induced angle);
    ``Gamma1_over_U`` [m].
    Assumptions: lifting-line assumptions; linear sections.
    Validation: V1 ``solve_lifting_line(twist=…)`` returns no coefficient beyond the first.  Label: analytic.
    """
    y = _F(y)
    s = float(s)
    if A is None:
        A = _area(c, s)
    G1U = 2.0 * float(CL_target) * float(A) / (np.pi * s)  # Eq. (14.21) with L = C_L ½ρU²A
    sg = np.sqrt(np.clip(1.0 - (2.0 * y / s) ** 2, 0.0, None))
    al = 2.0 * G1U * sg / (_at(K, y) * _at(c, y)) + G1U / (2.0 * s)  # Eq. (14.17) with (14.23), (14.24)
    return {"alpha_local": _S(al), "alpha_i": float(G1U / (2.0 * s)), "Gamma1_over_U": float(G1U)}


def elliptic_twist_rectangular(alpha_r: float, c: float, s: float, *, K: float) -> dict:
    """Rectangular wing with the elliptic twist α_w = α_r √(1 − (2y/s)²): the flight angle that makes its loading elliptic.

    Book: §14.6, Exercise 14.14 (result ours): matching the constant and the sin γ parts of (14.17) with (14.23),
    (14.24) gives Γ₁ = ½KUc α_r and α_t = Γ₁/(2sU) = Kc α_r/(4s).
    Parameters: alpha_r root twist [rad]; c chord [m]; s span [m]; K [1/rad] (keyword).
    Returns dict: ``alpha_t`` [rad]; ``Gamma1_over_U`` [m]; ``CL``; ``CDi``; ``AR``.
    Assumptions: uniform uncambered section; lifting-line assumptions.
    Validation: V1 ``solve_lifting_line`` with this twist and angle has Γ_n = 0 for n > 1.  Label: analytic.
    """
    G1U = 0.5 * float(K) * float(c) * float(alpha_r)
    AR = float(s) / float(c)
    CL = 0.5 * np.pi * G1U / float(c)  # (14.21) / (½ρU² s c)
    return {"alpha_t": float(K) * float(c) * float(alpha_r) / (4.0 * float(s)), "Gamma1_over_U": G1U, "CL": float(CL),
            "CDi": float(CL ** 2 / (np.pi * AR)), "AR": AR}


# ======================================================================================================================
# horseshoe vortices (Exercises 14.10, 14.11, 14.17, 14.18; Fig. 14.22)
# ======================================================================================================================
def _half_line(x, y, z, y0, z0=0.0):
    """Velocity per unit Γ of a straight semi-infinite filament from (0, y0, z0) to x = +∞ (circulation about +x)."""
    dy, dz = y - y0, z - z0
    d2 = dy ** 2 + dz ** 2
    r = np.sqrt(x ** 2 + d2)
    with np.errstate(divide="ignore", invalid="ignore"):
        f = (1.0 + x / r) / (_FOUR_PI * d2)  # (Γ/4πd)(cos θ_a − cos θ_b) with θ_b = π, D13 of ch05
    return 0.0 * f, -dz * f, dy * f


def horseshoe_velocity(x, y, z, Gamma: float, s: float, x_end: float | None = None) -> tuple:
    """Velocity induced by one horseshoe vortex: a bound segment on the y-axis and two trailing legs along +x.

    Book: §14.6, Exercises 14.10–14.11 and 14.17 (the formula of the straight segment is (5.17) integrated, ch05
    D13: (Γ/4πd)(cos θ_a − cos θ_b)); Fig. 14.22.  The bound vortex runs from (0, −s/2, 0) to (0, +s/2, 0) with
    circulation Γ (lift upward for a stream along +x); the starboard leg leaves along +x with +Γ, the port leg
    returns with −Γ.  The semi-infinite legs use the closed angle formula (not a long finite segment); with ``x_end``
    the legs are cut at x = x_end (three finite segments by ``core.biot_savart.segment_induced_velocity``).
    Parameters: x, y, z field point(s) [m] (broadcast); Gamma [m²/s]; s span of the horseshoe [m]; x_end [m] or None.
    Returns (ux, uy, uz) [m/s] in (x downstream, y span, z up) axes — the downwash is −uz.
    Assumptions: straight filaments of zero core size (infinite speed on them).
    Validation: V1 downwash at the centre Γ/(πs); far downstream on the centre line 2× that, far upstream 0; equals
    three finite segments as x_end → ∞ (1e-6).  Label: analytic.
    """
    x, y, z = np.broadcast_arrays(_F(x), _F(y), _F(z))
    shp = x.shape
    P = np.vstack([x.ravel(), y.ravel(), z.ravel()])
    h = 0.5 * float(s)
    G = float(Gamma)
    with np.errstate(divide="ignore", invalid="ignore"):
        ub = segment_induced_velocity(P, np.array([0.0, -h, 0.0]), np.array([0.0, h, 0.0]), G)
        ub = np.asarray(ub, dtype=float).reshape(3, -1)
        if x_end is None:
            a1 = _half_line(P[0], P[1], P[2], +h)
            a2 = _half_line(P[0], P[1], P[2], -h)
            u = [ub[k] + G * (a1[k] - a2[k]) for k in range(3)]
        else:
            xe = float(x_end)
            u1 = np.asarray(segment_induced_velocity(P, np.array([0.0, h, 0.0]), np.array([xe, h, 0.0]), G), float)
            u2 = np.asarray(segment_induced_velocity(P, np.array([xe, -h, 0.0]), np.array([0.0, -h, 0.0]), G), float)
            u = [ub[k] + u1.reshape(3, -1)[k] + u2.reshape(3, -1)[k] for k in range(3)]
    return tuple(_S(np.asarray(q).reshape(shp)) for q in u)


def horseshoe_downwash(x, y, Gamma: float, s: float, z=0.0):
    """Downwash (positive downward) of one horseshoe vortex at (x, y, z) [m/s] = −uz of :func:`horseshoe_velocity`.

    Book: §14.6, Exercises 14.11, 14.17.   Label: analytic.
    """
    return _S(-_F(horseshoe_velocity(x, y, z, Gamma, s)[2]))


def horseshoe_wing(Gamma: float, s: float, rho: float, U: float) -> dict:
    """The simplest wing — one horseshoe vortex of constant Γ: lift, centre downwash and the drag it suggests.

    Book: §14.6, Exercise 14.11 (results ours): L = ρUΓs; w(0) = Γ/(πs) (two half-lines at the distance s/2);
    using that centre value for the whole span, D_i = ρ w(0) Γ s = ρΓ²/π and C_Di = C_L²/(2πΛ) — half the elliptic
    value (14.26).  That "better than elliptic" drag is an artefact: the true downwash of a horseshoe grows without
    bound toward the tips (trap T10), so ∫ρwΓ dy diverges.
    Parameters: Gamma [m²/s]; s [m]; rho [kg/m³]; U [m/s].
    Returns dict: ``L`` [N], ``w0`` [m/s], ``Di_centre`` [N], ``Di_elliptic_same_lift`` [N], ``ratio`` (= ½),
    ``note``.
    Assumptions: constant Γ along the span (unphysical at the tips).
    Validation: V1 w0 from :func:`horseshoe_velocity`.  Label: analytic.
    """
    L = float(rho) * float(U) * float(Gamma) * float(s)
    w0 = float(Gamma) / (np.pi * float(s))
    Di = float(rho) * w0 * float(Gamma) * float(s)
    De = float(elliptic_induced_drag(L, rho, U, s))
    return {"L": L, "w0": w0, "Di_centre": Di, "Di_elliptic_same_lift": De, "ratio": Di / De,
            "note": "centre downwash used for the whole span; the true downwash is unbounded at the tips"}


def ground_effect_factor(h, s):
    """Factor by which the ground reduces the centre downwash of a horseshoe vortex, 1 − 1/(1 + (4h/s)²) [–].

    Book: §14.6, Exercise 14.18 (result ours): an image horseshoe of opposite strength at the depth h below the
    ground (2h below the wing) cancels the normal velocity on the ground; at the wing's centre its two legs, at the
    distance √((s/2)² + (2h)²), give an upwash Γ/(πs)·1/(1 + (4h/s)²).
    Parameters: h height of the wing above the ground [m]; s span [m].   Returns the factor in [0, 1).
    Assumptions: single horseshoe; flat ground; image method.
    Validation: V1 two :func:`horseshoe_velocity` calls; V7 h → ∞ → 1, h → 0 → 0.  Label: analytic.
    """
    return _S(1.0 - 1.0 / (1.0 + (4.0 * _F(h) / _F(s)) ** 2))


def downwash_along_x(x, Gn, s: float, y1: float = 0.0, include_bound: bool = False):
    """Downwash of a series loading at the points (x, y₁, 0) ahead of and behind the lifting line [m/s].

    Book: §14.6, Fig. 14.22 (zero far ahead, w at the lifting line, 2w far behind — "the actual profile can be
    determined from the Biot–Savart law").  Ours: each trailing filament is a half-line, so
    w(x, y₁) = (1/4π)∫ (dΓ/dy)/(y₁ − y)·[1 + x/√(x² + (y₁ − y)²)] dy; the singular part is w(0, y₁)(1 + sgn x) and
    the remainder is a regular integral (``scipy.integrate.quad`` in the span angle).  ``include_bound`` adds the
    bound vortex's own field ∫ xΓ dy/(4π(x² + (y₁ − y)²)^{3/2}) (upwash ahead, downwash behind, singular at x → 0);
    the book's figure shows the sheet only.
    Parameters: x [m] (positive downstream; any shape); Gn [m²/s]; s [m]; y1 [m]; include_bound (bool).
    Returns w [m/s] (downward positive).
    Assumptions: straight rigid wake in the plane z = 0.
    Validation: V7 x → −∞: 0; x = 0: (14.19); x → +∞: 2× (14.19).  Label: analytic.
    """
    g, n = _Gn(Gn)
    s = float(s)
    g1 = float(span_angle(y1, s))
    c1 = np.cos(g1)
    w0 = float(downwash_series(g, g1, s))

    def one(xx):
        if xx == 0.0:
            return w0
        sgn = np.sign(xx)

        def f(ga):
            d = 0.5 * s * (np.cos(ga) - c1)  # y₁ − y
            if abs(d) < 1e-14 * s:
                return 0.0
            return float(np.sum(n * g * np.cos(n * ga))) / d * (xx / np.hypot(xx, d) - sgn)

        val, _ = quad(f, 0.0, np.pi, points=[g1], limit=400, epsabs=1e-12, epsrel=1e-10)
        w = w0 * (1.0 + sgn) + val / _FOUR_PI
        if include_bound:
            def fb(ga):
                d = 0.5 * s * (np.cos(ga) - c1)
                return xx * float(np.sum(g * np.sin(n * ga))) * 0.5 * s * np.sin(ga) / (xx ** 2 + d ** 2) ** 1.5

            vb, _ = quad(fb, 0.0, np.pi, points=[g1], limit=400, epsabs=1e-12, epsrel=1e-10)
            w += vb / _FOUR_PI
        return w

    xv = _F(x)
    out = np.array([one(float(q)) for q in np.atleast_1d(xv).ravel()]).reshape(np.atleast_1d(xv).shape)
    return _S(out if xv.ndim else out[0])


def staircase_from_loading(Gamma_fn: Callable, s: float, N: int, spacing: str = "uniform") -> tuple:
    """Replace a loading Γ(y) by N panels of constant circulation (N nested horseshoe vortices side by side).

    Book: §14.6, Fig. 14.21 (a filament of strength dΓ leaves wherever the bound strength changes).
    Parameters: Gamma_fn callable Γ(y) [m²/s]; s [m]; N panels; spacing "uniform" (equal widths) or "cosine"
    (equal steps of the span angle).
    Returns (y_edges (N + 1,), Gamma_steps (N,)): panel edges [m] and the circulation on each panel [m²/s] (Γ at the
    panel mid-point).   Label: analytic.
    """
    N = int(N)
    if spacing == "uniform":
        ye = np.linspace(-0.5 * s, 0.5 * s, N + 1)
    elif spacing == "cosine":
        ye = -0.5 * float(s) * np.cos(np.linspace(0.0, np.pi, N + 1))
    else:
        raise ValueError('spacing must be "uniform" or "cosine"')
    return ye, _F(Gamma_fn(0.5 * (ye[:-1] + ye[1:])))


def staircase_downwash(y1, y_edges, Gamma_steps, parts: bool = False):
    """Downwash on the lifting line from a staircase loading: a finite sum of half-line filaments [m/s].

    Book: §14.6, the unnumbered dw(y₁) = −(dΓ/dy)dy/(4π(y − y₁)) before (14.13), summed over the filaments that leave
    at the panel edges: the filament at y_k carries Γ_left − Γ_right about +x and induces
    (Γ_left − Γ_right)/(4π(y_k − y₁)) downward at y₁ (half the velocity of an infinite line vortex).  The sum tends to
    the principal value (14.13) as the panels are refined (evaluate at panel mid-points).
    Parameters: y1 station(s) [m] (not on an edge); y_edges (M + 1,) [m]; Gamma_steps (M,) [m²/s]; parts → also
    return the (…, M + 1) array of each filament's contribution.
    Returns w(y₁) [m/s] (and the parts).
    Assumptions: straight semi-infinite filaments in the plane of the wing.
    Validation: V1 one panel at its centre: Γ/(πs); V3 → :func:`downwash` as the number of panels grows (first order).
    Label: analytic, converged.
    """
    ye, G = _F(y_edges), _F(Gamma_steps)
    jump = np.concatenate([[0.0], G]) - np.concatenate([G, [0.0]])  # Γ_left − Γ_right at each edge
    y1a = _F(y1)
    with np.errstate(divide="ignore", invalid="ignore"):
        contrib = jump / (_FOUR_PI * (ye - y1a[..., None]))
    w = np.sum(contrib, axis=-1)
    return (_S(w), contrib) if parts else _S(w)


# ======================================================================================================================
# the far wake
# ======================================================================================================================
def rolled_up_pair(Gn, s: float) -> dict:
    """The pair of tip vortices a loading rolls up into: circulation, spacing and descent speed.

    Book: §14.6, Fig. 14.20 (each half of the sheet rolls up into one vortex whose circulation equals the mid-span
    value Γ₀).  Ours: the impulse per unit length ρ∫Γ dy is conserved, so the spacing of the pair is
    b′ = (1/Γ₀)∫Γ dy = (πs/4)Γ₁/Γ₀ (πs/4 for the elliptic loading); the pair sinks at Γ₀/(2πb′).
    Parameters: Gn [m²/s]; s [m].   Returns dict: ``Gamma0`` [m²/s], ``spacing`` b′ [m], ``descent_speed`` [m/s].
    Assumptions: complete roll-up of each half into one line vortex; symmetric loading.
    Validation: V1 elliptic: spacing πs/4; V4 equals the centroid of the half-sheet.  Label: analytic.
    """
    g, n = _Gn(Gn)
    G0 = float(np.sum(g * np.sin(n * 0.5 * np.pi)))
    bp = 0.25 * np.pi * float(s) * g[0] / G0
    return {"Gamma0": G0, "spacing": float(bp), "descent_speed": float(G0 / (2.0 * np.pi * bp))}


def wake_kinetic_energy_per_length(Gn, s: float, rho: float, m: int = 256) -> float:
    """Kinetic energy of the cross-flow left behind per unit length of wake, ½ρ∫Γ w_∞ dy [J/m = N].

    Book: §14.6 ("the rate of generation of trailing-vortex kinetic energy must equal the rate of work done against
    the induced drag, namely D_iU" — stated).  Ours: far behind the wing the sheet's flow is two-dimensional with the
    potential jump Γ(y) across it and the normal velocity w_∞ = 2w (Fig. 14.22), so by Green's theorem the energy per
    unit length is ½ρ∫Γ w_∞ dy = ρ∫Γ w dy — the induced drag of (14.15).  Evaluated by the midpoint rule in the span
    angle (``m`` nodes), not from (14.22).
    Parameters: Gn [m²/s]; s [m]; rho [kg/m³]; m quadrature nodes.   Returns the energy per length [N].
    Assumptions: straight, rigid, planar wake.   Validation: V4 equals (14.22) to 1e-10.  Label: conserved.
    """
    ga = (np.arange(int(m)) + 0.5) * np.pi / int(m)
    G = _F(circulation_series(Gn, ga))
    w_inf = 2.0 * _F(downwash_series(Gn, ga, s))
    return float(0.5 * float(rho) * np.sum(G * w_inf * 0.5 * float(s) * np.sin(ga)) * np.pi / int(m))
