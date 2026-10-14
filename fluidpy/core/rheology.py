"""Generalised-Newtonian constitutive laws and their steady pipe flows: Newtonian, Casson (yield stress), power law,
Herschel–Bulkley.

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e (2012), Ch. 16, §16.2 "Nature of Blood" (Eqs. (16.5)–(16.7): the
yield-stress correlation, the Casson law in the plasma-scaled form, Whitmore's fit) and §16.3 "Laminar Flow of a Casson
Fluid in a Rigid-Walled Tube" (Eqs. (16.207)–(16.220)).  Every equation was transcribed from the rendered page images
(chapters/pages/ch16/p817, p866–p868).  The power-law and Herschel–Bulkley fluxes are ours (standard results, checked
by quadrature); the book names them only in its exercises.

Used by Ch. 16 at three sites (the rheology of whole blood, apparent viscosity, the Casson pipe) and by any later
non-Newtonian example.

Conventions (analysis ch16 §9, traps T5, T14, T21, T27)
* Shear rate γ̇ ≥ 0 [1/s] and shear stress τ ≥ 0 [Pa] are **magnitudes** (the book's τ = −μ du/dr > 0 in a pipe).
* Casson constants: the book writes the law twice — (16.6) √(τ/μ_p) = k_c√γ̇ + √(τ_y/μ_p) with a *dimensionless*
  k_c, and (16.207) τ^{1/2} = τ_y^{1/2} + K_c γ̇^{1/2} with a *dimensional* K_c [(Pa s)^{1/2}].  They are the same law
  with K_c = k_c √μ_p (:func:`casson_constants`).  K_c² [Pa s] is the high-shear (asymptotic) viscosity.
* Apparent viscosity μ_app = τ/γ̇ is the **secant** slope (as in (16.32)); :func:`differential_viscosity` gives the
  tangent dτ/dγ̇ (the book's words "slope at a point" describe the tangent; its formula is the secant — trap T27).
* Pipe: radius a [m], length L [m], pressure drop Δp = p₁ − p₂ > 0 [Pa] drives flow in +x; wall stress
  τ_w = aΔp/(2L); plug radius a_c = 2Lτ_y/Δp; no flow when a_c ≥ a (Δp ≤ 2Lτ_y/a).
* SI throughout; CGS (dyn/cm², poise) only through converters.  Every function is vectorised and scalar-callable.
"""
from __future__ import annotations

import numpy as np

from ._util import as_scalar_if_0d, require_positive

__all__ = ["DYN_CM2", "POISE", "WHITMORE_KC", "WHITMORE_SQRT_TAUY_OVER_MUP", "newtonian_stress", "apparent_viscosity",
           "casson_constants", "casson_stress", "casson_shear_rate", "casson_apparent_viscosity",
           "differential_viscosity", "casson_stress_plasma_form", "whitmore_casson", "yield_stress_correlation",
           "casson_wall_stress", "casson_min_pressure_drop", "casson_plug_radius", "casson_pipe_velocity",
           "casson_core_velocity", "casson_pipe_flow_rate", "casson_poiseuille_flux_ratio", "casson_flux_ratio",
           "power_law_pipe_flow_rate", "power_law_pipe_velocity", "herschel_bulkley_pipe_flow_rate"]

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d

DYN_CM2: float = 0.1     #: 1 dyn/cm² in Pa (exact)
POISE: float = 0.1       #: 1 P = 1 dyn s/cm² in Pa s (exact); 1 cP = 1e-3 Pa s
WHITMORE_KC: float = 1.53                 #: Whitmore's k_c in (16.7) [–]
WHITMORE_SQRT_TAUY_OVER_MUP: float = 2.0  #: Whitmore's √(τ_y/μ_p) in (16.7) [s^{-1/2}]


# ======================================================================================================================
# Constitutive laws
# ======================================================================================================================
def newtonian_stress(gamma_dot, mu):
    """Newtonian shear stress τ = μγ̇.

    Book: §16.2 "Nature of Blood" (first line of p. 790, the definition the chapter starts from; ch01 (1.7)).
    Parameters
    ----------
    gamma_dot : shear rate γ̇ [1/s] (≥ 0).   mu : dynamic viscosity μ [Pa s].
    Returns
    -------
    tau : shear stress [Pa].
    Assumptions: Newtonian fluid.
    Validation: V1 (the Casson law with τ_y = 0 and K_c² = μ reduces to it). Label: analytic.
    """
    return _S(_F(mu) * _F(gamma_dot))


def apparent_viscosity(tau, gamma_dot):
    """Apparent viscosity μ_app = τ/γ̇ — the **secant** slope of the flow curve through the origin (trap T27).

    Book: §16.2 (definition after (16.5)); the same definition by Poiseuille analogy in (16.32).
    Parameters
    ----------
    tau : shear stress [Pa];  gamma_dot : shear rate [1/s] (> 0).
    Returns
    -------
    mu_app [Pa s].
    Assumptions: none.  Validation: V1 (Newtonian ⇒ μ). Label: analytic.
    """
    g = _F(gamma_dot)
    require_positive("gamma_dot", g)
    return _S(_F(tau) / g)


def casson_constants(k_c, sqrt_tau_y_over_mu_p, mu_p):
    """Convert the plasma-scaled Casson constants of (16.6) into the dimensional ones of (16.207).

    Book: §16.2, Eq. (16.6)  √(τ/μ_p) = k_c√γ̇ + √(τ_y/μ_p), and §16.3, Eq. (16.207)  τ^{1/2} = τ_y^{1/2} + K_cγ̇^{1/2}.
    Multiplying (16.6) by √μ_p gives (16.207) with K_c = k_c√μ_p and τ_y = (√(τ_y/μ_p))² μ_p (derivation D02, ours).
    Parameters
    ----------
    k_c : Casson viscosity coefficient [–];  sqrt_tau_y_over_mu_p : √(τ_y/μ_p) [s^{-1/2}];  mu_p : plasma viscosity [Pa s].
    Returns
    -------
    (tau_y, Kc) : yield stress [Pa], Casson constant [(Pa s)^{1/2}].
    Assumptions: homogeneous fluid.  Validation: V2 (identity of the two forms, tested by
    :func:`casson_stress_plasma_form`). Label: symbolic.
    """
    mu_p = _F(mu_p)
    require_positive("mu_p", mu_p)
    return _S(_F(sqrt_tau_y_over_mu_p) ** 2 * mu_p), _S(_F(k_c) * np.sqrt(mu_p))


def casson_stress(gamma_dot, tau_y, Kc):
    """Casson shear stress at a given shear rate: τ = (√τ_y + K_c√γ̇)².

    Book: §16.3, Eq. (16.207) (the same law as (16.6)).  At γ̇ = 0 the value returned is τ_y, the stress at the yield
    surface (a material at rest may carry any stress ≤ τ_y).
    Parameters
    ----------
    gamma_dot : shear rate γ̇ [1/s] (≥ 0);  tau_y : yield stress [Pa] (≥ 0);  Kc : Casson constant [(Pa s)^{1/2}].
    Returns
    -------
    tau [Pa].
    Assumptions: homogeneous Casson fluid, steady simple shear.
    Validation: V1 inverse of :func:`casson_shear_rate`; V7 τ_y → 0 gives Newtonian with μ = K_c². Label: analytic.
    """
    g = _F(gamma_dot)
    if np.any(g < 0):
        raise ValueError("gamma_dot is a magnitude: must be >= 0")
    return _S((np.sqrt(_F(tau_y)) + _F(Kc) * np.sqrt(g)) ** 2)  # Eq. (16.207)


def casson_shear_rate(tau, tau_y, Kc):
    """Shear rate of a Casson fluid at a given stress: γ̇ = (√τ − √τ_y)²/K_c² for τ ≥ τ_y, and γ̇ = 0 below yield.

    Book: §16.3, Eq. (16.207) (piecewise), equivalently (16.211) γ̇ = [τ + τ_y − 2√(ττ_y)]/K_c².
    Parameters
    ----------
    tau : shear stress [Pa] (≥ 0);  tau_y : yield stress [Pa];  Kc : Casson constant [(Pa s)^{1/2}] (> 0).
    Returns
    -------
    gamma_dot [1/s].
    Assumptions: as :func:`casson_stress`.
    Validation: V1 round trip with :func:`casson_stress` to 1e-12; zero below yield. Label: analytic.
    """
    tau, tau_y, Kc = _F(tau), _F(tau_y), _F(Kc)
    require_positive("Kc", Kc)
    d = np.sqrt(np.maximum(tau, 0.0)) - np.sqrt(tau_y)
    return _S(np.where(tau >= tau_y, d ** 2 / Kc ** 2, 0.0))  # Eq. (16.207), (16.211)


def casson_apparent_viscosity(gamma_dot, tau_y, Kc):
    """Casson apparent (secant) viscosity μ_app = τ/γ̇ = (√(τ_y/γ̇) + K_c)².

    Book: §16.2 definition of μ_app with (16.6); Fig. 16.8 (ours: the closed form, derivation D02).  Limits: γ̇ → ∞
    gives K_c² (the asymptotic viscosity); γ̇ → 0 diverges like τ_y/γ̇.
    Parameters
    ----------
    gamma_dot : shear rate [1/s] (> 0);  tau_y : yield stress [Pa];  Kc : [(Pa s)^{1/2}].
    Returns
    -------
    mu_app [Pa s].
    Assumptions: as :func:`casson_stress`.
    Validation: V1 = casson_stress/γ̇; V7 both limits. Label: analytic.
    """
    g = _F(gamma_dot)
    require_positive("gamma_dot", g)
    return _S((np.sqrt(_F(tau_y) / g) + _F(Kc)) ** 2)


def differential_viscosity(gamma_dot, tau_y, Kc):
    """Casson differential (tangent) viscosity dτ/dγ̇ = K_c² + K_c√(τ_y/γ̇) (trap T27; ours).

    Book: §16.2 ("the slope … at any given point on the curve"), with (16.207); the closed form is ours.  Always below
    the secant :func:`casson_apparent_viscosity` (they differ by τ_y/γ̇ + K_c√(τ_y/γ̇)) and equal to it when τ_y = 0.
    Parameters
    ----------
    gamma_dot : [1/s] (> 0);  tau_y : [Pa];  Kc : [(Pa s)^{1/2}].
    Returns
    -------
    dtau_dgamma [Pa s].
    Assumptions: as :func:`casson_stress`.
    Validation: V1 central difference of :func:`casson_stress`. Label: analytic.
    """
    g = _F(gamma_dot)
    require_positive("gamma_dot", g)
    Kc = _F(Kc)
    return _S(Kc ** 2 + Kc * np.sqrt(_F(tau_y) / g))


def casson_stress_plasma_form(gamma_dot, k_c, sqrt_tau_y_over_mu_p, mu_p):
    """Casson stress from the plasma-scaled form (16.6): τ = μ_p (k_c√γ̇ + √(τ_y/μ_p))².

    Book: §16.2, Eq. (16.6).
    Parameters
    ----------
    gamma_dot [1/s];  k_c [–];  sqrt_tau_y_over_mu_p [s^{-1/2}];  mu_p plasma viscosity [Pa s].
    Returns
    -------
    tau [Pa].
    Assumptions: homogeneous fluid.
    Validation: V2 equals :func:`casson_stress` with :func:`casson_constants`. Label: symbolic.
    """
    return _S(_F(mu_p) * (_F(k_c) * np.sqrt(_F(gamma_dot)) + _F(sqrt_tau_y_over_mu_p)) ** 2)  # Eq. (16.6)


def whitmore_casson(gamma_dot, mu_p, *, k_c: float = WHITMORE_KC, s0: float = WHITMORE_SQRT_TAUY_OVER_MUP) -> dict:
    """Whitmore's least-squares Casson fit for whole blood: √(τ/μ_p) = 1.53√γ̇ + 2.0.

    Book: §16.2, Eq. (16.7); Figs. 16.7–16.8.  The constants carry hidden units: 1.53 is k_c (dimensionless) and 2.0 is
    s₀ = √(τ_y/μ_p) in s^{-1/2}, so τ_y = s₀²μ_p = 4μ_p (in Pa with μ_p in Pa s) and K_c² = 1.53²μ_p = 2.34μ_p
    (trap T14).
    Parameters
    ----------
    gamma_dot : shear rate [1/s] (> 0);  mu_p : plasma viscosity [Pa s];  k_c [–];  s0 [s^{-1/2}].
    Returns
    -------
    dict: tau [Pa], mu_app [Pa s] (secant), mu_tangent [Pa s] (dτ/dγ̇), tau_y [Pa], Kc [(Pa s)^{1/2}].
    Assumptions: homogeneous blood, the fit's range of shear rates.
    Validation: V1 against :func:`casson_stress` with :func:`casson_constants` (1.53, 2.0, μ_p). Label: analytic.
    """
    tau_y, Kc = casson_constants(k_c, s0, mu_p)
    tau = casson_stress_plasma_form(gamma_dot, k_c, s0, mu_p)  # Eq. (16.7)
    return dict(tau=tau, mu_app=casson_apparent_viscosity(gamma_dot, tau_y, Kc),
                mu_tangent=differential_viscosity(gamma_dot, tau_y, Kc), tau_y=tau_y, Kc=Kc)


def yield_stress_correlation(H, C_F, *, warn_range: bool = True):
    """Yield stress of blood from hematocrit and fibrinogen: √τ_y = (H − 0.1)(C_F + 0.5), τ_y in dyn/cm² — returned in Pa.

    Book: §16.2, Eq. (16.5): H the hematocrit **as a fraction** (> 0.1), C_F the fibrinogen content in g per 100 mL,
    stated range 0.21 < C_F < 0.46.  The hidden units of the empirical constants make τ_y come out in dyn/cm²; we
    multiply by 0.1 to return Pa (DEVIATION: SI return value; the correlation itself is unchanged).
    Parameters
    ----------
    H : hematocrit [–], 0.1 < H < 1 (a value > 1 is taken as a percent by mistake and raises — trap T13);
    C_F : fibrinogen concentration [g/100 mL];  warn_range : warn if C_F lies outside the stated range.
    Returns
    -------
    tau_y [Pa].
    Assumptions: the correlation's data range.
    Validation: V1 hand value; V7 pint conversion of dyn/cm² to Pa. Label: analytic.
    """
    H, C_F = _F(H), _F(C_F)
    if np.any(H > 1.0):
        raise ValueError("H is the hematocrit as a FRACTION (e.g. 0.45), not a percent (trap T13)")
    if np.any(H <= 0.1):
        raise ValueError("(16.5) holds only for H > 0.1")
    if warn_range and np.any((C_F < 0.21) | (C_F > 0.46)):
        import warnings

        warnings.warn("C_F outside the stated range 0.21 < C_F < 0.46 g/100 mL of (16.5)", stacklevel=2)
    sqrt_tau_dyn = (H - 0.1) * (C_F + 0.5)  # Eq. (16.5), √(dyn/cm²)
    return _S(sqrt_tau_dyn ** 2 * DYN_CM2)


# ======================================================================================================================
# Casson fluid in a rigid pipe (§16.3, (16.207)–(16.220))
# ======================================================================================================================
def casson_wall_stress(dp, L, a):
    """Wall shear stress of fully developed pipe flow, τ_w = aΔp/(2L) (any fluid: a force balance).

    Book: §16.3, Eq. (16.209) at r = a (τ_rx(r) = (p₁ − p₂)r/(2L)), from the cylinder force balance (16.208).  The
    first member of (16.209) is printed (r/2)dp/dx; with dp/dx = −Δp/L it should read −(r/2)dp/dx (slip S25, loose).
    Parameters
    ----------
    dp : pressure drop Δp = p₁ − p₂ [Pa];  L : length [m];  a : radius [m].
    Returns
    -------
    tau_w [Pa] (positive for Δp > 0).
    Assumptions: steady, fully developed, axisymmetric.
    Validation: V4 force balance τ_w 2πaL = πa²Δp. Label: conserved.
    """
    return _S(_F(a) * _F(dp) / (2.0 * _F(L)))  # Eq. (16.209) at r = a


def casson_min_pressure_drop(L, a, tau_y):
    """Smallest pressure drop that moves a yield-stress fluid through a pipe: (p₁ − p₂)_min = 2Lτ_y/a.

    Book: §16.3, text after (16.209) (τ_w = τ_y).  The book's worked example quotes the dyn/cm² value 1000× too small
    (slip S26); the code is the formula.
    Parameters
    ----------
    L : length [m];  a : radius [m];  tau_y : yield stress [Pa].
    Returns
    -------
    dp_min [Pa].
    Assumptions: as :func:`casson_wall_stress`.
    Validation: V1 τ_w(dp_min) = τ_y; :func:`casson_pipe_flow_rate` → 0 there. Label: analytic.
    """
    return _S(2.0 * _F(L) * _F(tau_y) / _F(a))


def casson_plug_radius(dp, L, tau_y):
    """Radius of the unsheared plug, a_c = 2Lτ_y/Δp = a τ_y/τ_w (where τ(r) = τ_y).

    Book: §16.3, Eq. (16.210); Fig. 16.21.  a_c ≥ a means the whole pipe is below yield: no flow (trap T21).
    Parameters
    ----------
    dp : Δp [Pa] (> 0);  L : length [m];  tau_y : yield stress [Pa].
    Returns
    -------
    a_c [m].
    Assumptions: as :func:`casson_wall_stress`.
    Validation: V1 τ(a_c) = τ_y. Label: analytic.
    """
    dp = _F(dp)
    require_positive("dp", dp)
    return _S(2.0 * _F(L) * _F(tau_y) / dp)  # Eq. (16.210)


def casson_core_velocity(a, dp, L, tau_y, Kc):
    """Speed of the Casson plug, u_c = (1/4K_c²)(Δp/L)(√a − √a_c)³(√a + ⅓√a_c); zero if a_c ≥ a.

    Book: §16.3, Eq. (16.218) (from (16.216) at r/a = a_c/a = τ_y/τ_w, (16.217)).
    Parameters
    ----------
    a : radius [m];  dp : Δp [Pa];  L : length [m];  tau_y : [Pa];  Kc : [(Pa s)^{1/2}].
    Returns
    -------
    u_c [m/s].
    Assumptions: steady, fully developed, rigid pipe, Casson fluid, no slip.
    Validation: V1 equals :func:`casson_pipe_velocity` at r = a_c; V7 τ_y → 0 gives the Poiseuille centre-line speed with
    μ = K_c². Label: analytic.
    """
    a, dp, L, Kc = _F(a), _F(dp), _F(L), _F(Kc)
    ac = np.minimum(_F(casson_plug_radius(dp, L, tau_y)), a)
    sa, sc = np.sqrt(a), np.sqrt(ac)
    return _S(dp / L / (4.0 * Kc ** 2) * (sa - sc) ** 3 * (sa + sc / 3.0))  # Eq. (16.218)


def casson_pipe_velocity(r, a, dp, L, tau_y, Kc):
    """Velocity profile of a Casson fluid in a pipe: a flat plug for r < a_c and the sheared annulus (16.215).

    Book: §16.3, Eq. (16.215)  u = (1/4K_c²)(Δp/L)[(a² − r²) − (8/3)√a_c(√a³ − √r³) + 2a_c(a − r)] for a_c ≤ r ≤ a,
    u = u_c (16.218) for r < a_c; (16.216) is the same in τ_w and τ_y.  Zero everywhere if a_c ≥ a.
    Parameters
    ----------
    r : radial position(s) [m], 0 ≤ r ≤ a;  a : radius [m];  dp : Δp [Pa] (> 0);  L : length [m];  tau_y : [Pa];
    Kc : [(Pa s)^{1/2}].
    Returns
    -------
    u [m/s].
    Assumptions: as :func:`casson_core_velocity`.
    Validation: V2 sympy: −du/dr = (√τ − √τ_y)²/K_c² with τ = Δp r/(2L), u(a) = 0, du/dr(a_c) = 0; V1 continuity at a_c.
    Label: symbolic, analytic.
    """
    r = _F(r)
    a, dp, L, Kc = float(a), float(dp), float(L), float(Kc)
    ac = float(casson_plug_radius(dp, L, tau_y))
    if ac >= a:
        return _S(np.zeros_like(r))
    rr = np.clip(r, ac, a)  # inside the plug the profile is flat at its value at a_c
    u = dp / L / (4.0 * Kc ** 2) * ((a ** 2 - rr ** 2) - 8.0 / 3.0 * np.sqrt(ac) * (a ** 1.5 - rr ** 1.5)
                                    + 2.0 * ac * (a - rr))  # Eq. (16.215)
    return _S(u)


def casson_pipe_flow_rate(a, dp, L, tau_y, Kc):
    """Volume flux of a Casson fluid in a pipe (16.220):
    Q = (π/8)(1/K_c²)(Δp/L)a⁴[1 − (16/7)ξ^{1/2} + (4/3)ξ − (1/21)ξ⁴], ξ = a_c/a; zero for ξ ≥ 1.

    Book: §16.3, Eq. (16.219) (plug + annulus) and Eq. (16.220) ("after considerable algebra"; verified by sympy,
    derivation D33).
    Parameters
    ----------
    a : radius [m];  dp : Δp [Pa];  L : length [m];  tau_y : [Pa];  Kc : [(Pa s)^{1/2}].
    Returns
    -------
    Q [m³/s].
    Assumptions: as :func:`casson_core_velocity`.
    Validation: V1 = πa_c²u_c + ∫2πru dr by quadrature; V7 τ_y → 0 ⇒ πa⁴Δp/(8K_c²L); ξ → 1 ⇒ 0. Label: analytic,
    symbolic.
    """
    a, dp, L, Kc = _F(a), _F(dp), _F(L), _F(Kc)
    xi = _F(casson_plug_radius(dp, L, tau_y)) / a
    bracket = 1.0 - 16.0 / 7.0 * np.sqrt(np.minimum(xi, 1.0)) + 4.0 / 3.0 * xi - xi ** 4 / 21.0
    Q = np.pi / 8.0 / Kc ** 2 * dp / L * a ** 4 * bracket  # Eq. (16.220)
    return _S(np.where(xi < 1.0, Q, 0.0))


def casson_poiseuille_flux_ratio(xi):
    """Q_Casson/Q_Poiseuille at the same Δp with μ = K_c²: 1 − (16/7)ξ^{1/2} + (4/3)ξ − (1/21)ξ⁴, ξ = a_c/a = τ_y/τ_w.

    Book: §16.3, the bracket of Eq. (16.220).  Zero for ξ ≥ 1; the leading loss is (16/7)√ξ (curation E1 'aha').
    Parameters
    ----------
    xi : a_c/a [–] (≥ 0).
    Returns
    -------
    ratio [–].
    Assumptions: as :func:`casson_pipe_flow_rate`.
    Validation: V1 (16.220)/(16.17). Label: analytic.
    """
    xi = _F(xi)
    b = 1.0 - 16.0 / 7.0 * np.sqrt(np.minimum(np.maximum(xi, 0.0), 1.0)) + 4.0 / 3.0 * xi - xi ** 4 / 21.0
    return _S(np.where(xi < 1.0, b, 0.0))


casson_flux_ratio = casson_poiseuille_flux_ratio  # design name (C.2)


# ======================================================================================================================
# Other generalised-Newtonian pipe flows (ours; the book names them in Exercises 16.5–16.6 only)
# ======================================================================================================================
def power_law_pipe_velocity(r, a, dp, L, m, n):
    """Velocity of a power-law fluid τ = m γ̇ⁿ in a pipe: u = [n/(n + 1)](τ_w/m)^{1/n} a [1 − (r/a)^{1 + 1/n}] (ours).

    Book: none (the power law is named in Exercise 16.5; standard result, derived from τ = Δp r/(2L) of (16.209)).
    Parameters
    ----------
    r [m];  a [m];  dp : Δp [Pa];  L [m];  m : consistency [Pa sⁿ];  n : flow index [–] (> 0).
    Returns
    -------
    u [m/s].
    Assumptions: steady, fully developed, rigid pipe.
    Validation: V7 n = 1 ⇒ Poiseuille with μ = m. Label: analytic.
    """
    tw = _F(casson_wall_stress(dp, L, a))
    a, n = float(a), float(n)
    return _S(n / (n + 1.0) * (tw / float(m)) ** (1.0 / n) * a * (1.0 - (_F(r) / a) ** (1.0 + 1.0 / n)))


def power_law_pipe_flow_rate(a, dp, L, m, n):
    """Flux of a power-law fluid τ = m γ̇ⁿ in a pipe: Q = [πn/(3n + 1)] a³ (τ_w/m)^{1/n} (ours).

    Book: none (Exercise 16.5 names the problem; the result is standard).
    Parameters
    ----------
    a [m];  dp : Δp [Pa];  L [m];  m [Pa sⁿ];  n [–].
    Returns
    -------
    Q [m³/s].
    Assumptions: as :func:`power_law_pipe_velocity`.
    Validation: V1 quadrature of :func:`power_law_pipe_velocity`; V7 n = 1 ⇒ πa⁴Δp/(8mL). Label: analytic.
    """
    tw = _F(casson_wall_stress(dp, L, a))
    n = _F(n)
    return _S(np.pi * n / (3.0 * n + 1.0) * _F(a) ** 3 * (tw / _F(m)) ** (1.0 / n))


def herschel_bulkley_pipe_flow_rate(a, dp, L, m, n, tau0):
    """Flux of a Herschel–Bulkley fluid τ = τ₀ + m γ̇ⁿ (τ > τ₀) in a pipe (ours):
    Q = πa³(τ_w/m)^{1/n}(1 − φ)^{1+1/n}[(1 − φ)²/(3 + 1/n) + 2φ(1 − φ)/(2 + 1/n) + φ²/(1 + 1/n)], φ = τ₀/τ_w.

    Book: none (Exercise 16.6 names the problem; standard result, here derived by integrating by parts
    Q = (πa³/τ_w³)∫₀^{τ_w} τ² γ̇(τ) dτ).
    Parameters
    ----------
    a [m];  dp : Δp [Pa];  L [m];  m [Pa sⁿ];  n [–];  tau0 : yield stress [Pa].
    Returns
    -------
    Q [m³/s] (0 when τ_w ≤ τ₀).
    Assumptions: as :func:`power_law_pipe_velocity`.
    Validation: V7 τ₀ = 0 ⇒ :func:`power_law_pipe_flow_rate`; V1 quadrature of the Rabinowitsch integral. Label: analytic.
    """
    tw = _F(casson_wall_stress(dp, L, a))
    n, m = _F(n), _F(m)
    phi = np.minimum(_F(tau0) / tw, 1.0)
    k = 1.0 / n
    q = (np.pi * _F(a) ** 3 * (tw / m) ** k * (1.0 - phi) ** (1.0 + k)
         * ((1.0 - phi) ** 2 / (3.0 + k) + 2.0 * phi * (1.0 - phi) / (2.0 + k) + phi ** 2 / (1.0 + k)))
    return _S(np.where(phi < 1.0, q, 0.0))
