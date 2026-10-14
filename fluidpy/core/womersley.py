"""Oscillatory (pulsatile) flow in a rigid circular pipe — Womersley's solution — and the complex Bessel functions
it needs.  The cylindrical sibling of Stokes' second problem (``core.laminar.stokes_second_problem``).

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e (2012), §16.3: the Womersley number (16.12); "Pulsatile Flow in a
Rigid Cylindrical Tube" (16.69)–(16.100).  Every equation was transcribed from the rendered page images
(chapters/pages/ch16/p823, p837–p840).  Re-exported by ``fluidpy.ch16_biofluid_mechanics``; the viscous elastic-tube
waves of the same chapter (16.130)–(16.136) call :func:`bessel_ratio_j0_j1`.

Conventions (analysis ch16 §9, traps T1, T4, T25)
* **Time convention e^{+iωt}** (the book's (16.73)–(16.75)): −∂p/∂x = G e^{iωt}, u = U(r)e^{iωt}, the physical
  velocity is the **real part** and belongs to the gradient −∂p/∂x = G cos ωt.  The Bessel argument is i^{3/2}α r/a
  with i^{3/2} = e^{3iπ/4}.  (The elastic-tube analysis uses e^{i(kx − ωt)} — never mix the two.)
* G = −dp/dx is the pressure-gradient amplitude [Pa/m] (ch08 ``G=``).  A harmonic series is
  −∂p/∂x = G₀ + Re Σ_{n≥1} G_n e^{inωt} (G₀ real, G_n complex); harmonic n has α_n = α√n (16.86).
* α = a√(ω/ν) (16.12) uses the **radius** a; Re = 2ūa/ν and St = 2aω/ū are the chapter's diameter forms, so
  α² = Re·St/4 (trap T6).
* Large α: J₀(i^{3/2}α) grows like e^{α/√2}; every ratio of Bessel functions is evaluated with the exponentially
  scaled ``scipy.special.jve`` so nothing overflows (trap T25).
* Wall shear stress is reported as the positive-for-forward-flow magnitude τ_w = −μ ∂u/∂r|_a (the book's (16.92) calls
  it a "shear rate" and writes ∂u/∂r — slip S8).
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache

import numpy as np
from scipy.special import jve

from ._util import as_scalar_if_0d, require_positive

__all__ = ["I12", "I32", "bessel_j0_series", "bessel_j1_series", "bessel_ratio_j0_j1", "womersley_number",
           "stokes_number", "womersley_from_re_st", "womersley_F1", "womersley_F2", "womersley_F2_series",
           "womersley_F2_asymptote", "womersley_profile_complex", "womersley_velocity", "womersley_velocity_nd",
           "womersley_harmonics", "womersley_flow_rate", "womersley_flow_rate_harmonics", "womersley_wall_shear",
           "womersley_XY",
           "womersley_flow_ratio", "womersley_phase_lag", "womersley_small_alpha", "womersley_large_alpha",
           "womersley_ode_residual", "womersley_regime", "womersley_pde_solve"]

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d

I12: complex = complex(np.exp(0.25j * np.pi))   #: i^{1/2} = e^{iπ/4}
I32: complex = complex(np.exp(0.75j * np.pi))   #: i^{3/2} = e^{3iπ/4}


def _S_c(x):
    """Scalar-if-0d for complex results (keeps the complex type)."""
    arr = np.asarray(x)
    return complex(arr) if arr.ndim == 0 else arr


# ======================================================================================================================
# Bessel functions of complex argument
# ======================================================================================================================
def bessel_j0_series(z, nterms: int = 60):
    """J₀(z) by its power series Σ_m (−1)^m (z/2)^{2m}/(m!)² for complex z (the from-scratch version; ours).

    Book: §16.3, the J₀ of (16.77)–(16.81) (the series is the standard definition, not printed in the book).
    Parameters
    ----------
    z : complex argument (|z| ≲ 25 for double-precision accuracy: the terms grow to ~e^{|z|} before they cancel);
    nterms : number of terms (60 is enough for |z| ≤ 25).
    Returns
    -------
    J₀(z) (complex).
    Assumptions: none.
    Validation: V1 against ``scipy.special.jv`` to 1e-12 relative for |z| ≤ 15. Label: analytic.
    """
    z = np.asarray(z, dtype=complex)
    term = np.ones_like(z)
    s = term.copy()
    q = -(z / 2.0) ** 2
    for m in range(1, int(nterms)):
        term = term * q / (m * m)
        s = s + term
    return _S_c(s)


def bessel_j1_series(z, nterms: int = 60):
    """J₁(z) by its power series (z/2) Σ_m (−1)^m (z/2)^{2m}/(m!(m + 1)!) for complex z (ours).

    Book: §16.3, the J₁ of (16.91) (standard series).
    Parameters
    ----------
    z : complex argument;  nterms : number of terms.
    Returns
    -------
    J₁(z) (complex).
    Assumptions: none.  Validation: V1 against ``scipy.special.jv`` (|z| ≤ 15, 1e-12). Label: analytic.
    """
    z = np.asarray(z, dtype=complex)
    term = z / 2.0
    s = term.copy()
    q = -(z / 2.0) ** 2
    for m in range(1, int(nterms)):
        term = term * q / (m * (m + 1))
        s = s + term
    return _S_c(s)


def bessel_ratio_j0_j1(z):
    """J₀(z)/J₁(z) for complex z without overflow (both scaled by e^{−|Im z|} with ``scipy.special.jve``).

    Book: §16.3, the ratio J₀(β̂)/J₁(β̂) of (16.136) and 2J₁/(zJ₀) of (16.91).
    Parameters
    ----------
    z : complex argument (≠ 0).
    Returns
    -------
    ratio (complex).
    Assumptions: none.  Validation: V1 against ``jv`` ratio where ``jv`` does not overflow. Label: analytic.
    """
    z = np.asarray(z, dtype=complex)
    return _S_c(jve(0, z) / jve(1, z))


def _j0_ratio(z_in, z_out):
    """J₀(z_in)/J₀(z_out) without overflow: jve(0, z_in)/jve(0, z_out) · e^{|Im z_in| − |Im z_out|}."""
    z_in = np.asarray(z_in, dtype=complex)
    z_out = np.asarray(z_out, dtype=complex)
    return jve(0, z_in) / jve(0, z_out) * np.exp(np.abs(z_in.imag) - np.abs(z_out.imag))


# ======================================================================================================================
# Dimensionless groups
# ======================================================================================================================
def womersley_number(a, omega, nu):
    """Womersley number α = a√(ω/ν): the tube radius against the oscillatory (Stokes) layer.

    Book: §16.3, Eq. (16.12).  α = √2 a/δ_e with δ_e = √(2ν/ω) of ch08's Stokes layer (ours).
    Parameters
    ----------
    a : tube radius [m];  omega : angular frequency ω [rad/s] (heart rate in beats/min × 2π/60);
    nu : kinematic viscosity [m²/s].
    Returns
    -------
    alpha [–].
    Assumptions: none.  Validation: V1 hand value; V1 = ½√(Re·St) (:func:`womersley_from_re_st`). Label: analytic.
    """
    nu = _F(nu)
    require_positive("nu", nu)
    return _S(_F(a) * np.sqrt(_F(omega) / nu))  # Eq. (16.12)


def stokes_number(a, omega, nu):
    """Stokes number α² = a²ω/ν (the square of the Womersley number, as the book names it).

    Book: §16.3, text after (16.12).
    Parameters
    ----------
    a [m];  omega [rad/s];  nu [m²/s].
    Returns
    -------
    alpha² [–].
    Assumptions: none.  Validation: V1. Label: analytic.
    """
    return _S(_F(womersley_number(a, omega, nu)) ** 2)


def womersley_from_re_st(Re, St):
    """Womersley number from the chapter's Reynolds and Strouhal numbers, α = ½√(Re·St) (ours, derivation a-d45).

    Book: §16.3, text after (16.12): Re = 2au/ν and St = 2aω/u (diameter forms), so Re·St = 4a²ω/ν = 4α².
    (Other texts define α² = 2π Re St with different Re, St — trap T6.)
    Parameters
    ----------
    Re : 2au/ν [–];  St : 2aω/u [–].
    Returns
    -------
    alpha [–].
    Assumptions: the chapter's definitions.  Validation: V1 equals :func:`womersley_number`. Label: analytic.
    """
    return _S(0.5 * np.sqrt(_F(Re) * _F(St)))


def womersley_regime(alpha, *, low: float = 1.0, high: float = 10.0) -> str:
    """Name the regime of oscillatory pipe flow from α: "quasi-steady" (α < low), "transitional", "inertial" (α > high)
    (ours).

    Book: §16.3, text after (16.12) ("for α < 1, viscous effects are dominant … quasi-steady") and the limits (16.87),
    (16.88).  The upper threshold is our choice (the plug + Stokes-layer picture is clear above about 10); both are
    keywords so the explainer can show them.
    Parameters
    ----------
    alpha : Womersley number [–] (scalar);  low, high : thresholds [–].
    Returns
    -------
    str.
    Assumptions: none.  Validation: V1 (three test values). Label: analytic.
    """
    alpha = float(alpha)
    if alpha < low:
        return "quasi-steady"
    if alpha > high:
        return "inertial"
    return "transitional"


# ======================================================================================================================
# The two Womersley functions
# ======================================================================================================================
def womersley_F1(r, a, alpha):
    """F₁ = J₀(i^{3/2}α r/a)/J₀(i^{3/2}α) — the shape function of the oscillatory profile.

    Book: §16.3, Eq. (16.81) (book writes F₁(α); it depends on r/a too).  Evaluated with exponentially scaled Bessel
    functions so it never overflows (trap T25); F₁ → 1 everywhere as α → 0.
    Parameters
    ----------
    r : radial position(s) [m], 0 ≤ r ≤ a;  a : radius [m];  alpha : Womersley number [–] (≥ 0).
    Returns
    -------
    F1 (complex) [–].
    Assumptions: rigid tube, e^{+iωt} convention.
    Validation: V1 against :func:`bessel_j0_series` ratio for α ≤ 10; F₁(a) = 1. Label: analytic.
    """
    s = _F(r) / float(a)
    alpha = float(alpha)
    if alpha == 0.0:
        return _S_c(np.ones_like(s, dtype=complex))
    return _S_c(_j0_ratio(I32 * alpha * s, I32 * alpha))  # Eq. (16.81)


def womersley_F2(alpha):
    """F₂(α) = 2J₁(i^{3/2}α)/(i^{3/2}α J₀(i^{3/2}α)) — the flow-rate and wall-shear function.

    Book: §16.3, Eq. (16.91) (from ∫₀ᵃ J₀(λr) r dr = (a/λ)J₁(λa), derivation D16).  F₂(0) = 1 (limit).
    Parameters
    ----------
    alpha : Womersley number [–] (≥ 0), scalar or array.
    Returns
    -------
    F2 (complex) [–].
    Assumptions: rigid tube, e^{+iωt}.
    Validation: V1 = 1 − (8/a²)∫₀ᵃ ... (quadrature of F₁); V7 small-α series (16.93) error ∝ α⁴ and corrected large-α
    asymptote (16.96) error ∝ α⁻³. Label: analytic.
    """
    al = np.asarray(alpha, dtype=float)
    z = I32 * al
    with np.errstate(invalid="ignore", divide="ignore"):
        out = 2.0 / (z * bessel_ratio_j0_j1(np.where(al == 0.0, 1.0, z)))  # Eq. (16.91), J₁/J₀ via the scaled ratio
    out = np.where(al == 0.0, 1.0 + 0.0j, out)
    # very small α: 2J₁/(zJ₀) loses digits through jve's internal series; use the convergent series below 1e-3
    small = (al > 0.0) & (al < 1e-3)
    if np.any(small):
        out = np.where(small, _F2_series_eval(al, 6), out)
    return _S_c(out)


@lru_cache(maxsize=16)
def _F2_series_coeffs(order: int) -> tuple:
    """Exact coefficients c_k of 2J₁(z)/(zJ₀(z)) = Σ_k c_k x^k, x = z², up to x^order (power-series division)."""
    from math import factorial

    n = order + 1
    num = [Fraction((-1) ** m, 4 ** m * factorial(m) * factorial(m + 1)) for m in range(n)]  # 2J₁(z)/z
    den = [Fraction((-1) ** m, 4 ** m * factorial(m) ** 2) for m in range(n)]                # J₀(z)
    c: list[Fraction] = []
    for k in range(n):
        s = num[k] - sum((den[j] * c[k - j] for j in range(1, k + 1)), Fraction(0))
        c.append(s / den[0])
    return tuple(c)


def _F2_series_eval(alpha, order: int):
    x = (I32 * np.asarray(alpha, dtype=float)) ** 2  # x = z² = −iα²
    c = _F2_series_coeffs(order)
    out = np.zeros_like(x, dtype=complex)
    for k in range(order, -1, -1):
        out = out * x + float(c[k])
    return out


def womersley_F2_series(alpha, order: int = 2):
    """Small-α series of F₂ to the power α^order: order 2 gives the book's F₂ ≈ 1 − iα²/8 (16.93); higher even orders
    come from exact power-series division of 2J₁(z)/(zJ₀(z)) in x = z² = −iα² (order 4 adds −α⁴/48; ours).

    Book: §16.3, Eq. (16.93)  F₂(α) ≈ 1 − iα²/8 − O(α⁴).
    Parameters
    ----------
    alpha : Womersley number [–];  order : highest power of α kept (even, ≥ 2; 2 reproduces (16.93)).
    Returns
    -------
    F2 approximation (complex).
    Assumptions: α ≪ 1.
    Validation: V7 error against :func:`womersley_F2` ∝ α^{order+2} (measured slope). Label: analytic.
    """
    order = int(order)
    if order < 2 or order % 2:
        raise ValueError("order must be an even integer >= 2 (powers of alpha)")
    return _S_c(_F2_series_eval(alpha, order // 2))  # Eq. (16.93) for order = 2


def womersley_F2_asymptote(alpha, *, printed_16_96: bool = False):
    """Large-α asymptote of F₂: corrected F₂ ≈ (2/(i^{1/2}α))(1 + i^{3/2}/(2α)) = 2/(i^{1/2}α) + i/α².

    Book: §16.3, Eq. (16.96), printed F₂ ≈ (2/(i^{1/2}α))(1 + 1/(2α)).
    DEVIATION: slip S10 — the printed correction 1/(2α) should be i^{3/2}/(2α) (from J₁/J₀ → +i[1 − i/(2z) + …] for
    Im z → +∞, z = i^{3/2}α; numerically J₁/J₀ = −0.00712 + 0.99293i at α = 50, against +i[1 − i/(2z)] =
    −0.00707 + 0.99293i).  The printed form is less accurate than the leading term alone (error ∝ α⁻² against
    α⁻³ for the corrected form); ``printed_16_96=True`` returns it for the test that must fail.
    Parameters
    ----------
    alpha : Womersley number [–] (> 0);  printed_16_96 : return the printed form.
    Returns
    -------
    F2 approximation (complex).
    Assumptions: α ≫ 1.
    Validation: V7 error slopes against :func:`womersley_F2` (corrected −3, printed −2). Label: analytic.
    """
    al = np.asarray(alpha, dtype=float)
    require_positive("alpha", al)
    corr = 1.0 / (2.0 * al) if printed_16_96 else I32 / (2.0 * al)
    return _S_c(2.0 / (I12 * al) * (1.0 + corr))  # Eq. (16.96) (corrected unless printed_16_96)


# ======================================================================================================================
# Velocity, flow rate, wall shear
# ======================================================================================================================
def womersley_profile_complex(r, a, alpha):
    """Complex profile shape U(r)/(Ga²/μ) = (1 − F₁)/(iα²) of one harmonic, so that u = Re[G a²/μ · shape · e^{iωt}].

    Book: §16.3, Eqs. (16.80), (16.82): U(r) = −(iG/ωρ)(1 − F₁) = (Ga²/(iμα²))(1 − F₁).  α → 0 limit (1 − r²/a²)/4.
    Parameters
    ----------
    r [m];  a [m];  alpha [–] (≥ 0).
    Returns
    -------
    shape (complex) [–].
    Assumptions: rigid tube, e^{+iωt}.
    Validation: V2 residual of (16.76) with the corrected right side −G/μ (:func:`womersley_ode_residual`).
    Label: symbolic, analytic.
    """
    s = _F(r) / float(a)
    alpha = float(alpha)
    if alpha < 1e-3:  # quasi-steady series: (1 − F₁)/(iα²) = (1 − s²)/4 − iα²(3 − 4s² + s⁴)/64 + …
        return _S_c((1.0 - s ** 2) / 4.0 - 1j * alpha ** 2 * (3.0 - 4.0 * s ** 2 + s ** 4) / 64.0)
    return _S_c((1.0 - _F1_arr(s, alpha)) / (1j * alpha ** 2))  # Eq. (16.82)


def _F1_arr(s, alpha):
    return _j0_ratio(I32 * alpha * np.asarray(s, dtype=float), I32 * alpha)


def womersley_velocity(r, t, a, omega, rho, mu, G):
    """Womersley's oscillatory velocity in a rigid pipe for the gradient −∂p/∂x = G cos ωt.

    Book: §16.3, Eqs. (16.73)–(16.82): u(r, t) = Re[(Ga²/(iμα²))(1 − F₁)e^{iωt}], α = a√(ωρ/μ).
    Parameters
    ----------
    r : radial position(s) [m];  t : time(s) [s] (broadcast with r);  a : radius [m];  omega : ω [rad/s];
    rho : density [kg/m³];  mu : viscosity [Pa s];  G : gradient amplitude −∂p/∂x [Pa/m] (real; a complex G shifts the
    phase: −∂p/∂x = Re G e^{iωt}).
    Returns
    -------
    u [m/s] (real).
    Assumptions: rigid straight pipe, fully developed, axisymmetric, Newtonian, long wave (no convective inertia).
    Validation: V1 F.D. residual of (16.72); V3 :func:`womersley_pde_solve` (Crank–Nicolson, order 2); V7 α → 0
    Poiseuille, α → ∞ Stokes layer. Label: analytic, converged.
    """
    alpha = float(womersley_number(a, omega, float(mu) / float(rho)))
    shape = np.asarray(womersley_profile_complex(r, a, alpha), dtype=complex)
    u = (G * float(a) ** 2 / float(mu)) * shape * np.exp(1j * float(omega) * _F(t))  # Eq. (16.82)
    return _S(np.real(u))


def womersley_velocity_nd(r_over_a, wt, alpha):
    """Dimensionless Womersley velocity u/(Ga²/4μ) = Re[4(1 − F₁)/(iα²) e^{iωt}] (for explainer parity; ours).

    Book: §16.3, Eq. (16.82) divided by the steady centre-line speed Ga²/(4μ) of (16.16), so the α → 0 limit is
    (1 − (r/a)²) cos ωt.
    Parameters
    ----------
    r_over_a : r/a [–];  wt : phase ωt [rad];  alpha : Womersley number [–].
    Returns
    -------
    u* [–].
    Assumptions: as :func:`womersley_velocity`.
    Validation: V1 equals :func:`womersley_velocity` scaled. Label: analytic.
    """
    shape = np.asarray(womersley_profile_complex(r_over_a, 1.0, float(alpha)), dtype=complex)
    return _S(np.real(4.0 * shape * np.exp(1j * _F(wt))))


def _harmonics(G_n):
    """(G₀, [G₁, G₂, …]) from a scalar (pure harmonic n = 1, no mean) or a sequence (G₀ first)."""
    arr = np.atleast_1d(np.asarray(G_n, dtype=complex))
    if np.ndim(G_n) == 0:
        return 0.0, arr
    return float(np.real(arr[0])), arr[1:]


def womersley_harmonics(r, t, a, omega, rho, mu, G_n):
    """Velocity for a periodic gradient written as a Fourier series −∂p/∂x = G₀ + Re Σ_{n≥1} G_n e^{inωt}.

    Book: §16.3, Eqs. (16.83)–(16.86): u = u₀(r) + Σ u_n(r)e^{inωt}, u₀ = G₀a²/(4μ)(1 − r²/a²) (Poiseuille),
    u_n = G_n a²/(iμα_n²)(1 − F₁(α_n)), α_n = α√n.  The book writes the series without "Re" (complex shorthand).
    Parameters
    ----------
    r [m];  t [s] (broadcast);  a [m];  omega : fundamental ω [rad/s];  rho [kg/m³];  mu [Pa s];
    G_n : scalar G (a single harmonic n = 1, no mean) or sequence [G₀, G₁, …] [Pa/m] (G₀ real, G_n complex).
    Returns
    -------
    u [m/s] (real).
    Assumptions: as :func:`womersley_velocity`; linear superposition.
    Validation: V1 superposition of :func:`womersley_velocity` and :func:`fluidpy.core.laminar.pipe_poiseuille`.
    Label: analytic.
    """
    G0, Gs = _harmonics(G_n)
    r_, t_ = np.broadcast_arrays(_F(r), _F(t))
    a, mu = float(a), float(mu)
    alpha = float(womersley_number(a, omega, mu / float(rho)))
    u = G0 * a ** 2 / (4.0 * mu) * (1.0 - (r_ / a) ** 2) + 0.0 * t_  # Eq. (16.85)
    for n, Gn in enumerate(Gs, start=1):
        if Gn == 0:
            continue
        shape = np.asarray(womersley_profile_complex(r_, a, alpha * np.sqrt(n)), dtype=complex)  # α_n = α√n
        u = u + np.real(Gn * a ** 2 / mu * shape * np.exp(1j * n * float(omega) * t_))  # Eq. (16.86)
    return _S(u)


def womersley_flow_rate(t, a, omega, rho, mu, G):
    """Volume flow rate of oscillatory pipe flow for −∂p/∂x = G cos ωt: Q = Re[(πa⁴/(iμα²))G(1 − F₂(α))e^{iωt}]
    (a sequence [G₀, G₁, …] gives the harmonic sum Q₀ + Re Σ πa⁴G_n/(iμα_n²)(1 − F₂(α_n))e^{inωt}).

    Book: §16.3, Eqs. (16.90)–(16.91) (single mode) and (16.89) (harmonics; Q₀ = πG₀a⁴/(8μ), (16.95)); (16.98) is the
    real/imaginary split for one mode.
    Parameters
    ----------
    t : time(s) [s];  a [m];  omega [rad/s];  rho [kg/m³];  mu [Pa s];  G : scalar amplitude [Pa/m] (single harmonic
    n = 1, no mean — the meaning shared with :func:`womersley_harmonics` and :func:`womersley_wall_shear`) or
    sequence [G₀, G₁, …] (G₀ the mean).  :func:`womersley_flow_rate_harmonics` refuses a bare scalar.
    Returns
    -------
    Q [m³/s] (real).
    Assumptions: as :func:`womersley_velocity`.
    Validation: V1 quadrature ∫₀ᵃ 2πr u dr of :func:`womersley_harmonics` to 1e-10. Label: analytic.
    """
    G0, Gs = _harmonics(G)
    t_ = _F(t)
    a, mu = float(a), float(mu)
    alpha = float(womersley_number(a, omega, mu / float(rho)))
    Q = np.pi * G0 * a ** 4 / (8.0 * mu) + 0.0 * t_  # Eq. (16.95) (mean part)
    for n, Gn in enumerate(Gs, start=1):
        if Gn == 0:
            continue
        an = alpha * np.sqrt(n)
        coef = np.pi * a ** 4 / mu * _flow_shape(an)
        Q = Q + np.real(Gn * coef * np.exp(1j * n * float(omega) * t_))  # Eq. (16.90)
    return _S(Q)


def womersley_flow_rate_harmonics(t, a, omega, rho, mu, G_n):
    """Flow rate for a harmonic series of gradients, G_n = [G₀, G₁, …] (G₀ the mean): the sum (16.89) to N harmonics.

    Book: §16.3, Eq. (16.89) (printed to ∞; the sum runs to the N harmonics of (16.83)–(16.84), slip S7, loose).
    Parameters: t [s];  a [m];  omega [rad/s];  rho;  mu;  G_n sequence [G₀, G₁, …] [Pa/m] (a one-element sequence
    [G₀] is a steady mean gradient).  Returns Q [m³/s].
    A bare scalar raises ValueError: :func:`womersley_flow_rate` reads a scalar as the n = 1 amplitude (no mean),
    while a "harmonics" call could be read as the mean G₀ — the two readings differ by orders of magnitude, so the
    ambiguous input is refused (pass [G₀] for a mean, or call :func:`womersley_flow_rate` for one harmonic).
    Validation: V1 quadrature of :func:`womersley_harmonics`. Label: analytic.
    """
    if np.ndim(G_n) == 0:
        raise ValueError("womersley_flow_rate_harmonics needs a sequence [G0, G1, ...]; a scalar is ambiguous "
                         "(use [G0] for a steady mean, or womersley_flow_rate(t, a, omega, rho, mu, G1) for one "
                         "harmonic)")
    G_n = np.atleast_1d(np.asarray(G_n, dtype=complex))
    if G_n.size == 1:
        G_n = np.concatenate([G_n, [0.0]])
    return womersley_flow_rate(t, a, omega, rho, mu, list(G_n))  # Eq. (16.89)


def _flow_shape(alpha):
    """(1 − F₂)/(iα²) = Σ_{k≥1} c_k x^{k−1} (x = −iα², the series of :func:`womersley_F2_series`) for small α —
    limit 1/8 − iα²/48 — and the Bessel form otherwise."""
    alpha = float(alpha)
    if alpha < 1e-2:
        x = -1j * alpha ** 2
        c = _F2_series_coeffs(4)
        return sum(float(c[k]) * x ** (k - 1) for k in range(1, 5))
    return (1.0 - complex(womersley_F2(alpha))) / (1j * alpha ** 2)


def womersley_wall_shear(t, a, omega, rho, mu, G_n):
    """Wall shear stress τ_w = −μ ∂u/∂r|_{r=a} = G₀a/2 + Re Σ (a/2)G_n F₂(α_n)e^{inωt}.

    Book: §16.3, Eq. (16.92), printed as a "shear rate" τ = ∂u/∂r with an unnamed F.
    DEVIATION: slip S8 — the printed values are the wall shear *stress* −μ∂u/∂r, and F is F₂ of (16.91) (from
    dJ₀/dz = −J₁: −μ ∂u_n/∂r|_a = (G_n a/2)·2J₁(z)/(zJ₀(z))).
    Parameters
    ----------
    t [s];  a [m];  omega [rad/s];  rho [kg/m³];  mu [Pa s];  G_n : scalar G or [G₀, G₁, …] [Pa/m].
    Returns
    -------
    tau_w [Pa] (positive when the wall drags against forward flow).
    Assumptions: as :func:`womersley_velocity`.
    Validation: V1 central difference of :func:`womersley_harmonics` at the wall (1e-8). Label: analytic.
    """
    G0, Gs = _harmonics(G_n)
    t_ = _F(t)
    a = float(a)
    alpha = float(womersley_number(a, omega, float(mu) / float(rho)))
    tau = G0 * a / 2.0 + 0.0 * t_
    for n, Gn in enumerate(Gs, start=1):
        if Gn == 0:
            continue
        F2 = complex(womersley_F2(alpha * np.sqrt(n)))
        tau = tau + np.real(a / 2.0 * Gn * F2 * np.exp(1j * n * float(omega) * t_))  # Eq. (16.92), corrected
    return _S(tau)


def womersley_XY(alpha):
    """Real and imaginary parts of Z = 1 − F₂(α) = X + iY.

    Book: §16.3, Eq. (16.97).  Small α: X ≈ α⁴/48, Y ≈ α²/8; large α: Z → 1.
    Parameters
    ----------
    alpha : Womersley number [–].
    Returns
    -------
    (X, Y) [–].
    Assumptions: as :func:`womersley_F2`.
    Validation: V1 from :func:`womersley_F2`. Label: analytic.
    """
    Z = 1.0 - np.asarray(womersley_F2(alpha), dtype=complex)  # Eq. (16.97)
    return _S(np.real(Z)), _S(np.imag(Z))


def womersley_flow_ratio(alpha):
    """Flow amplitude relative to the steady (Poiseuille) value at the same gradient amplitude:
    |Q|/|Q₀| = 8|Z|/α² = 8√(X² + Y²)/α².

    Book: §16.3, Eqs. (16.95), (16.99) (the ratio is ours).  → 1 as α → 0, ≈ 8/α² as α → ∞ (inertia-limited).
    Parameters
    ----------
    alpha : Womersley number [–] (≥ 0).
    Returns
    -------
    ratio [–].
    Assumptions: as :func:`womersley_F2`.
    Validation: V1 |womersley_flow_rate|/(πGa⁴/8μ). Label: analytic.
    """
    al = np.asarray(alpha, dtype=float)
    out = np.empty_like(al)
    flat = al.reshape(-1)
    res = out.reshape(-1)
    for i, v in enumerate(flat):
        res[i] = 8.0 * abs(_flow_shape(v))
    return _S(out)


def womersley_phase_lag(alpha):
    """Phase lag φ of the flow behind the gradient G cos ωt: tan φ = X/Y, φ ∈ [0, π/2) [rad].

    Book: §16.3, Eq. (16.100).  From (16.98), Re Q = |Q| cos(ωt − φ) with cos φ = Y/|Z|, sin φ = X/|Z|, so
    φ = arctan2(X, Y): 0 as α → 0 (in phase, viscous), → π/2 as α → ∞ (a driven mass).
    Parameters
    ----------
    alpha : Womersley number [–].
    Returns
    -------
    phi [rad].
    Assumptions: as :func:`womersley_F2`.
    Validation: V1 the time of maximum of :func:`womersley_flow_rate`. Label: analytic.
    """
    al = np.asarray(alpha, dtype=float)
    out = np.empty_like(al)
    flat = al.reshape(-1)
    res = out.reshape(-1)
    for i, v in enumerate(flat):
        s = _flow_shape(v) * 1j * max(v, 1e-300) ** 2  # = Z = X + iY (with the small-α series)
        res[i] = np.arctan2(s.real, s.imag) if v > 0 else 0.0  # Eq. (16.100)
    return _S(out)


def womersley_small_alpha(r, t, a, omega, rho, mu, G):
    """Quasi-steady limit α → 0: u ≈ (Ga²/4μ)(1 − r²/a²) cos ωt (Poiseuille following the gradient).

    Book: §16.3, Eq. (16.87).  Same signature as :func:`womersley_velocity` (a drop-in ghost curve).
    Parameters
    ----------
    r [m];  t [s];  a [m];  omega [rad/s];  rho [kg/m³] (unused, kept for the common signature);  mu [Pa s];  G [Pa/m].
    Returns
    -------
    u [m/s].
    Assumptions: α ≪ 1.
    Validation: V7 relative difference from :func:`womersley_velocity` ∝ α². Label: analytic.
    """
    a = float(a)
    return _S(np.real(G * a ** 2 / (4.0 * float(mu)) * (1.0 - (_F(r) / a) ** 2) * np.exp(1j * float(omega) * _F(t))))


def womersley_large_alpha(r, t, a, omega, rho, mu, G):
    """Large-α limit: u ≈ Re[(Ga²/(iμα²)){1 − exp[−√(ω/2ν)(1 + i)(a − r)]}e^{iωt}] — a plug ρ∂u/∂t = G cos ωt with a
    Stokes layer at the wall.

    Book: §16.3, Eq. (16.88) ("Stokes' second problem" of ch08, δ_e = √(2ν/ω)).  For harmonic n pass ω → nω (the
    book's exponent omits n — slip S7, loose).
    Parameters
    ----------
    r [m];  t [s];  a [m];  omega [rad/s];  rho [kg/m³];  mu [Pa s];  G [Pa/m].
    Returns
    -------
    u [m/s].
    Assumptions: α ≫ 1.
    Validation: V7 difference from :func:`womersley_velocity` → 0 as α grows (∝ 1/α near the wall). Label: analytic.
    """
    a, mu, rho, omega = float(a), float(mu), float(rho), float(omega)
    nu = mu / rho
    alpha2 = a ** 2 * omega / nu
    k = np.sqrt(omega / (2.0 * nu))
    shape = (1.0 - np.exp(-k * (1.0 + 1j) * (a - _F(r)))) / (1j * alpha2)
    return _S(np.real(G * a ** 2 / mu * shape * np.exp(1j * omega * _F(t))))  # Eq. (16.88)


def womersley_ode_residual(r, a, omega, rho, mu, G, *, h: float | None = None):
    """Left side of (16.76), U″ + U′/r − (iωρ/μ)U, for the coded U(r), by central differences (ours).

    Book: §16.3, Eq. (16.76), printed with right side +G/μ.
    DEVIATION: slip S6 — substituting (16.75) into (16.74) gives right side **−G/μ**; the particular solution
    G/(iωρ) of (16.77) is consistent only with −G/μ.  This function returns the left side so a test can show it equals
    −G/μ (and not +G/μ).
    Parameters
    ----------
    r : radial position(s) [m], strictly inside (0, a);  a [m];  omega [rad/s];  rho [kg/m³];  mu [Pa s];  G [Pa/m];
    h : difference step [m] (default a·1e-4).
    Returns
    -------
    lhs (complex) [1/(m s)] — compare with −G/μ.
    Assumptions: as :func:`womersley_velocity`.
    Validation: V1 (the residual test itself). Label: analytic.
    """
    a, mu, rho, omega = float(a), float(mu), float(rho), float(omega)
    h = a * 1e-4 if h is None else float(h)
    alpha = float(womersley_number(a, omega, mu / rho))

    def U(x):
        return G * a ** 2 / mu * np.asarray(womersley_profile_complex(x, a, alpha), dtype=complex)

    r_ = _F(r)
    Up, U0, Um = U(r_ + h), U(r_), U(r_ - h)
    return _S_c((Up - 2.0 * U0 + Um) / h ** 2 + (Up - Um) / (2.0 * h) / r_ - 1j * omega * rho / mu * U0)


# ======================================================================================================================
# Independent numerical route: Crank–Nicolson in r for (16.72)
# ======================================================================================================================
def womersley_pde_solve(a, omega, rho, mu, G, *, nr: int = 200, steps_per_period: int = 400, periods: int = 8,
                        periodic: bool = True, fast: bool = False) -> dict:
    """Solve ρ∂u/∂t = G cos ωt + μ(∂²u/∂r² + (1/r)∂u/∂r), u(a, t) = 0, by Crank–Nicolson (ours).

    Book: §16.3, Eq. (16.72) with (16.73) (the PDE whose periodic solution is (16.82)).  The scheme is ours:
    second-order central differences on nodes r_j = jΔr, j = 0…nr (u_nr = 0 at the wall), the symmetry row at r = 0
    uses ∇²u → 4(u₁ − u₀)/Δr² (L'Hôpital: u_r/r → u_rr), Crank–Nicolson in time with the forcing at mid-step.
    DEVIATION: with ``periodic=True`` (default) the periodic state is found directly — the map over one period is
    affine, u(T) = M u(0) + b, so u(0) = (I − M)⁻¹b — instead of marching from rest until the start-up transient
    (decay rate ≈ 5.78ν/a²) has died, which takes many periods at large α.  ``periodic=False`` marches ``periods``
    periods from rest.
    Parameters
    ----------
    a [m];  omega [rad/s];  rho [kg/m³];  mu [Pa s];  G [Pa/m];  nr : radial intervals;  steps_per_period : time steps
    per period;  periods : periods to march (only with periodic=False);  periodic : solve for the periodic state;
    fast : coarse grid (nr = 80, 200 steps) for notebooks.
    Returns
    -------
    dict: r (nr + 1,) [m], t (steps_per_period + 1,) [s] over one period starting at a multiple of T, u (nt, nr + 1)
    [m/s], dr, dt.
    Cost: ≈ 0.2 s at the defaults (periodic), ≈ 0.05 s fast.
    Assumptions: as :func:`womersley_velocity`.
    Validation: V3 observed order 2 (Δr and Δt halved together) against :func:`womersley_velocity`. Label: converged.
    """
    from scipy.linalg import solve_banded

    if fast:
        nr, steps_per_period = 80, 200
    a, omega, rho, mu, G = float(a), float(omega), float(rho), float(mu), float(G)
    nu = mu / rho
    nr, ns = int(nr), int(steps_per_period)
    r = np.linspace(0.0, a, nr + 1)
    dr = r[1] - r[0]
    T = 2.0 * np.pi / omega
    dt = T / ns
    n = nr  # unknowns u_0 … u_{nr−1}; u_nr = 0
    # Laplacian matrix L (tridiagonal) on the unknowns
    lo = np.zeros(n)
    di = np.zeros(n)
    up = np.zeros(n)
    di[0], up[0] = -4.0 / dr ** 2, 4.0 / dr ** 2
    j = np.arange(1, n)
    rj = r[1:n]
    lo[1:] = 1.0 / dr ** 2 - 1.0 / (2.0 * dr * rj)
    di[1:] = -2.0 / dr ** 2
    up[1:] = 1.0 / dr ** 2 + 1.0 / (2.0 * dr * rj)
    del j
    c = 0.5 * nu * dt
    ab = np.zeros((3, n))  # (I − c L) in banded form
    ab[0, 1:] = -c * up[:-1]
    ab[1, :] = 1.0 - c * di
    ab[2, :-1] = -c * lo[1:]

    def apply_rhs(U):
        """(I + c L) U for U of shape (n,) or (n, m)."""
        out = (1.0 + c * di)[:, None] * U if U.ndim == 2 else (1.0 + c * di) * U
        if U.ndim == 2:
            out[:-1] += c * up[:-1, None] * U[1:]
            out[1:] += c * lo[1:, None] * U[:-1]
        else:
            out[:-1] += c * up[:-1] * U[1:]
            out[1:] += c * lo[1:] * U[:-1]
        return out

    def step(U, tk, forcing=True):
        rhs = apply_rhs(U)
        if forcing:
            f = G / rho * np.cos(omega * (tk + 0.5 * dt)) * dt  # mid-step forcing (second order)
            rhs = rhs + (f if U.ndim == 1 else f)
        return solve_banded((1, 1), ab, rhs)

    if periodic:
        # period map: M (homogeneous propagation of the identity) and b (forced response from rest)
        Mmat = np.eye(n)
        b = np.zeros(n)
        for k in range(ns):
            Mmat = solve_banded((1, 1), ab, apply_rhs(Mmat))
            b = step(b, k * dt)
        u0 = np.linalg.solve(np.eye(n) - Mmat, b)
    else:
        u0 = np.zeros(n)
        for k in range(ns * int(periods)):
            u0 = step(u0, k * dt)
    U = np.empty((ns + 1, nr + 1))
    U[:, -1] = 0.0
    U[0, :-1] = u0
    cur = u0
    for k in range(ns):
        cur = step(cur, k * dt)
        U[k + 1, :-1] = cur
    return dict(r=r, t=np.arange(ns + 1) * dt, u=U, dr=dr, dt=dt)
