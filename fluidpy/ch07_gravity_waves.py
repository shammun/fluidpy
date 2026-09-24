"""Chapter 7 — Gravity waves: the linear free-surface problem and its solution (φ, u, w, ψ, p′, the dispersion relation
derived), particle orbits and Stokes drift, wave energy and energy flux, surface tension, standing waves and seiches,
group velocity and kinematic wave theory (rays, refraction), nonlinear shallow-water steepening, the hydraulic jump,
Stokes waves, KdV and the solitary wave, interfacial and two-layer waves (barotropic and baroclinic modes, reduced
gravity) and internal waves in a continuously stratified fluid (ω = N cos θ, K·u = 0, c ⟂ c_g, beams, F = c_g E) —
plus re-exports of the reusable wave primitives in ``fluidpy.core.waves``.

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e (2012), Ch. 7, §§7.1–7.8, Eqs. (7.1)–(7.159), Figs. 7.1–7.33.
Every equation was transcribed from the rendered page images (chapters/pages/ch07/p282–p331, printed pp. 255–304).

Where the physics lives (all public names of ``core.waves`` are re-exported, so ``ch07.<name>`` reaches every callable)
--------------------------------------------------------------------------------------------------------------------
* ``core.waves`` (W) — wave vocabulary (7.1)–(7.9); ω(k) of every family ((7.28), (7.56), (7.95), (7.110)–(7.113),
  (7.137)); inverse dispersion; c, c_g (analytic and for any ω by complex step, vector form); depth regimes; capillary
  minimum; beats, packets, FFT evolution (1-D, 2-D), rays; wave energy density; reduced gravity (7.117).
* this module — the chapter's fields (surface, standing, interfacial, two-layer, internal), boundary-condition residuals,
  sympy derivations, particle paths and Stokes drift, energy integrals, seiches, kinematic-wave checks, nonlinear waves
  (simple-wave steepening, hydraulic jump, Stokes waves, KdV solver, solitary and cnoidal waves) and the explainer
  state functions (``*_state``).

Conventions (analysis §9 R1–R11, curation §8)
* x horizontal (propagation), **z up**, still surface z = 0, bottom z = −H (Fig. 7.2); ``H = np.inf`` = deep water.
  §7.7 two-layer problem: origin at the mean free surface, interface at z = −H.
* **η is the surface elevation** η(x, t) here (the surface function is f = z − η, (7.14)); in ch04 (4.90) η was the
  level-set function itself — call ``core.interfaces.kinematic_bc_residual`` with z − η, never with the elevation.
* p is **gauge** pressure; p′ ≡ p + ρgz in §7.2 (7.30) but p′ = p − p̄(z) in §7.8 (7.124).
* ψ: u = ∂ψ/∂z, w = −∂ψ/∂x (ch04's planar convention with y → z; GFD often uses the opposite sign — Ch. 13 trap).
* ω = angular frequency (not vorticity); ζ = particle vertical excursion (§7.2, §7.6, §7.8) *and* interface
  displacement (§7.7) — code keys ``zeta`` (particle) / ``zeta`` of an interface are documented per function; θ = local
  phase θ(x, t) (7.72) *and* the angle of K above the horizontal (7.139) (``theta_phase`` / ``theta_K``).
* Vortex-sheet strength γ = u_below − u_above (ch05 counterclockwise convention, ``ch05.vortex_sheet_strength``).
* Energies: surface/interfacial E per unit **horizontal area** and F per unit crest **length**; internal-wave E per
  unit **volume** and F per unit **area** (text after (7.159)).
* g = 9.81 m/s² (``G_BOOK``) in this chapter's functions; ρ = 1000 kg/m³ unless stated.

Book slips handled (never coded as printed; printed variants kept only as labelled options for wrong-variant tests):
T1 (7.66) envelope ½Δω **x** → ½Δω **t** (``beat_wave(printed=True)``); T2 "u from (7.28)" before (7.44) means (7.27);
T3 (7.105) e^{i(k**z** − ωt)} → e^{i(k**x** − ωt)} (``two_layer_residuals(printed_7_105=True)``); T4 (7.98) ∂φ₁/**d**z;
T5 "y = 0" → z = 0 (p. 254); T6 "(7.88)" in the Ursell remark means (7.87); the printed (7.138)/(7.145) assume k > 0
(``internal_wave_velocities(printed=True)``); T10 (p. 288) the interfacial E_p line's middle form
(g(ρ₂ − ρ₁)/2λ)∫₀^{λ/2}ζ²dx equals ⅛(ρ₂ − ρ₁)ga², not ¼ — its first and last forms agree (¼), so the middle one needs
1/λ instead of 1/(2λ); :func:`interface_energy` uses the consistent ¼; T11 (p. 288) the text cites Exercise 7.16 for
the interfacial E_k — it is Exercise 7.18. Also: Exercise 7.2 set up literally (potential amplitude fixed at aω/k)
gives γ = 3/8, not the book's γ = 1 of (7.83) — see :func:`stokes_expansion_sympy`. Book-quoted numbers stay in the
git-ignored ``tests/book_values_ch07.json``.

Reduced gravity — which g′? ``ch07.reduced_gravity`` is the re-exported ch04 function
``core.similarity.reduced_gravity`` whose default ``ref="upper"`` is g(ρ₂ − ρ₁)/ρ₁ (ρ₁ in the denominator, after
(4.105)) — **not** this chapter's (7.117) g′ = g(ρ₂ − ρ₁)/ρ₂. For (7.117) call ``reduced_gravity_book(rho1, rho2)``
(default ``ref="lower"``) or ``reduced_gravity(rho1, rho2, ref="lower")``. The re-export is kept for backward
compatibility (ch04 recap R10).
"""
from __future__ import annotations

import functools
import hashlib
from pathlib import Path
from typing import Callable

import numpy as np
import sympy as sp
from scipy import integrate
from scipy.integrate import quad, solve_ivp
from scipy.special import ellipe, ellipj, ellipk

from .core import waves as W
from .core._util import as_scalar_if_0d
from .core.potential import laplacian_residual
from .core.similarity import froude_number, reduced_gravity  # noqa: F401  (re-exported: R10 recap, g′ conventions)
from .core.stratification import brunt_vaisala_sq, parcel_displacement  # noqa: F401  (re-exported: R18, N145)
from .core.thermo import G_BOOK
from .core.waves import *  # noqa: F401,F403
from .core.waves import (cosh_over_cosh, cosh_over_sinh, group_velocity, group_velocity_numeric, internal_wave_omega,
                         internal_wave_velocities, omega_capillary_gravity, omega_gravity, phase_speed,
                         sinh_over_sinh, wavenumber_from_omega)

G = G_BOOK  #: default g [m/s²] of this chapter's functions (the book's 9.81)
STOKES_LIMIT_STEEPNESS = 0.1410633
"""Height-to-wavelength ratio (crest to trough, H/λ) of the limiting (highest) deep-water Stokes wave, ± 4·10⁻⁷,
about one-seventh — the classical limiting steepness of the strongly-nonlinear-wave literature reviewed by L. W.
Schwartz & J. D. Fenton, "Strongly nonlinear waves", Annu. Rev. Fluid Mech. 14, 39–60 (1982). Digits as quoted by
HandWiki "Physics:Stokes wave" (https://handwiki.org/wiki/Physics:Stokes_wave, read 2026-09-24; it cites Schwartz &
Fenton for the one-seventh statement and footnotes the ± digits to Dyachenko, Lushnikov & Korotkevich 2016, whose
arXiv abstract 1507.02784 gives only the 2π/3 = 120° crest angle, Stokes 1880). The book's a_max ≈ 0.07λ (§7.6, after
(7.83)) is about half of it (a is measured from the mean level). Label: benchmark."""

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d
_TWO_PI = 2.0 * np.pi


def _omega(k, H, g, sigma=0.0, rho=1000.0):
    return float(omega_capillary_gravity(k, H, sigma, rho, g)) if sigma else float(omega_gravity(k, H, g))


def _dz_deep(k, H):
    """A depth at which a deep-water field is negligible (e^{−40}) or the bottom for finite H."""
    return -float(H) if np.isfinite(H) else -40.0 / float(k)


# ======================================================================================================================
# §7.2 the linear free-surface problem
# ======================================================================================================================
def surface_normal(eta_x):
    """Upward unit normal of the free surface f = z − η(x, t) = 0: n = (−(∂η/∂x) e_x + e_z)/√((∂η/∂x)² + 1).

    Book: §7.2, Eq. (7.14). ⚠️ η is the elevation here; ch04 (4.90) used η for the level-set function itself.
    Parameters: eta_x surface slope ∂η/∂x [–] (array allowed). Returns n with components on the first axis, shape
    (2, …) [–]. Validation: V1, V2 — tests/test_ch07.py: test_kinematic_condition_V2_derivation,
    test_surface_normal_V1_unit_normal_and_kinematic_parity. Checks: V1 unit length, parallel to ∇(z − η). Label:
    analytic.
    """
    s = _F(eta_x)
    den = np.sqrt(s ** 2 + 1.0)
    return np.stack([-s / den, 1.0 / den + 0.0 * s])  # Eq. (7.14)


def surface_velocity(eta_t):
    """Velocity of the surface point at fixed x: U_s = (∂η/∂t) e_z [m/s] (only its normal part matters).

    Book: §7.2, Eq. (7.15). Parameters: eta_t ∂η/∂t [m/s]. Returns (2, …) array. Validation: V1 — tests/test_ch07.py:
    test_surface_normal_V1_unit_normal_and_kinematic_parity. Checks: V1 (n·U_s)|∇f| = ∂η/∂t (the right side of (7.16)).
    Label: analytic.
    """
    e = _F(eta_t)
    return np.stack([0.0 * e, e])  # Eq. (7.15)


def linear_bernoulli_pressure(dphi_dt, z, rho: float = 1000.0, g: float = G) -> dict:
    """Gauge pressure from the linearised unsteady Bernoulli equation ∂φ/∂t + p/ρ + gz ≅ 0 and its perturbation part.

    Book: §7.2, Eqs. (7.20) (|∇φ|² of (4.83) dropped; constant fixed on the undisturbed surface far away) and (7.30)
    p′ ≡ p + ρgz = −ρ ∂φ/∂t. Parameters: dphi_dt [m²/s²]; z [m] (up); rho [kg/m³]; g [m/s²].
    Returns dict(p [Pa] gauge, p_prime [Pa]). Validation: V1 — tests/test_ch07.py:
    test_linear_bernoulli_V1_parity_with_ch04_unsteady_bernoulli. Checks: V1 parity with
    ``core.bernoulli.unsteady_bernoulli_pressure(speed=0)``; at rest p = −ρgz. Label: analytic.
    """
    pt = _F(dphi_dt)
    p = -float(rho) * (pt + float(g) * _F(z))  # Eq. (7.20)
    return {"p": _S(p), "p_prime": _S(p + float(rho) * float(g) * _F(z))}  # Eq. (7.30)


def wave_fields(x, z, t, a: float = 0.1, k: float = 1.0, H=np.inf, g: float = G, rho: float = 1000.0,
                direction: int = +1, sigma: float = 0.0) -> dict:
    """The linear progressive surface wave: η, φ, ψ, u, w and p′ at (x, z, t) (broadcast).

    θ = kx − sωt (s = direction ±1):
    η = a cos θ (7.2)/(7.61); φ = s(aω/k)[cosh k(z + H)/sinh kH] sin θ (7.26); u = s aω[cosh k(z + H)/sinh kH] cos θ,
    w = s aω[sinh k(z + H)/sinh kH] sin θ (7.27); ψ = s(aω/k)[sinh k(z + H)/sinh kH] cos θ (7.37);
    p′ = ρ(aω²/k)[cosh k(z + H)/sinh kH] cos θ = ρga[cosh k(z + H)/cosh kH] cos θ (7.31) (with σ: ρa(g + σk²/ρ)[…]).
    Deep water (H = ∞): every ratio → e^{kz} ((7.47), (7.48)). The left-going wave (s = −1) is the mirror image x → −x.
    ω from (7.28) (σ = 0) or (7.56).
    Book: §7.2, Eqs. (7.26), (7.27), (7.31), (7.37), (7.47)–(7.48), (7.51)–(7.52); §7.4 (7.61).
    Parameters: x, z [m] (z up, −H ≤ z ≲ η); t [s]; a amplitude [m]; k [rad/m] > 0; H depth [m]; g [m/s²]; rho
    [kg/m³]; direction ±1; sigma surface tension [N/m] (changes ω only).
    Returns dict(eta [m], phi [m²/s], psi [m²/s], u, w [m/s], p_prime [Pa], omega [rad/s], c [m/s]).
    Assumptions: inviscid, irrotational, constant density, ka ≪ 1, air ignored. Overflow-safe hyperbolic ratios.
    Validation: V1, V2, V3, V7 — tests/test_ch07.py: test_wave_fields_V1_closed_forms_and_streamfunction,
    test_wave_fields_V2_potential_relations_symbolic, test_wave_fields_V3_laplacian_residual_is_second_order_truncation,
    test_wave_fields_V7_deep_shallow_mirror_and_overflow, test_dispersion_V1_parity_with_ch04_and_identity. Checks: V2
    ∇²φ = 0, u = φ_x = ψ_z, w = φ_z = −ψ_x; V1 parity with ``ch04.linear_wave_surface``; V3 FD Laplacian order; V7
    deep/shallow limits, no overflow at kH = 500. Label: analytic, symbolic.
    """
    if direction not in (1, -1):
        raise ValueError("direction must be +1 or -1")
    s = float(direction)
    k = float(k)
    om = _omega(k, H, g, sigma, rho)
    x_, z_, t_ = np.broadcast_arrays(_F(x), _F(z), _F(t))
    th = k * x_ - s * om * t_
    Cz = _F(cosh_over_sinh(k, z_, H))
    Sz = _F(sinh_over_sinh(k, z_, H))
    Pz = _F(cosh_over_cosh(k, z_, H))
    g_eff = float(g) + float(sigma) * k ** 2 / float(rho)
    return {"eta": _S(float(a) * np.cos(th)),  # Eq. (7.2)
            "phi": _S(s * float(a) * om / k * Cz * np.sin(th)),  # Eq. (7.26)
            "psi": _S(s * float(a) * om / k * Sz * np.cos(th)),  # Eq. (7.37)
            "u": _S(s * float(a) * om * Cz * np.cos(th)),  # Eq. (7.27)
            "w": _S(s * float(a) * om * Sz * np.sin(th)),  # Eq. (7.27)
            "p_prime": _S(float(rho) * float(a) * g_eff * Pz * np.cos(th)),  # Eq. (7.31)
            "omega": om, "c": om / k}


def pressure_response(k, z, H=np.inf):
    """Pressure response factor p′/(ρga cos θ) = cosh k(z + H)/cosh kH [–] — the bottom-gauge transfer function.

    Book: §7.2, Eq. (7.31); deep limit e^{kz} (7.48) (4 % at z = −λ/2: a bottom sensor is a low-pass filter); shallow
    limit 1 (7.52, hydrostatic). Parameters: k [rad/m]; z [m] (−H ≤ z ≤ 0); H [m]. Overflow-safe.
    Validation: V1 — tests/test_ch07.py: test_pressure_response_V1_surface_bottom_and_limits. Checks: V1 = 1 at z = 0,
    1/cosh kH at the bottom, e^{−π} = 0.0432 at z = −λ/2 (deep). Label: analytic.
    """
    return cosh_over_cosh(k, z, H)  # Eq. (7.31)


def free_surface_residuals(x, t, a: float = 0.01, k: float = 1.0, H=np.inf, g: float = G, rho: float = 1000.0,
                           sigma: float = 0.0, h: float | None = None, direction: int = +1) -> dict:
    """Residuals of the field equation and of the exact and linearised boundary conditions for the linear wave.

    * laplace: ∇²φ (7.11) at mid-depth by the O(h⁴) stencil of ``core.potential.laplacian_residual`` (y ↦ z) [1/s];
    * bottom: w at z = −H (7.12) (deep: at z = −40/k) [m/s];
    * kinematic_exact: (∂φ/∂z)_{z=η} − ∂η/∂t − (∂η/∂x)(∂φ/∂x)_{z=η} (7.16) [m/s] — O(ka)·aω (the neglected slope and
      Taylor terms);
    * kinematic_17: (∂φ/∂z)_{z=η} − ∂η/∂t (7.17) [m/s]; kinematic_linear: (∂φ/∂z)_{z=0} − ∂η/∂t (7.18) [m/s] (= 0);
    * dynamic_exact: (∂φ/∂t + ½|∇φ|² + gη − (σ/ρ)κ)_{z=η}, the full Bernoulli (4.83) with (7.19)/(7.53), κ the exact
      curvature [m²/s²] — O(ka)·ag; dynamic_linear: (∂φ/∂t)_{z=0} + gη − (σ/ρ)η_xx (7.21)/(7.55) [m²/s²] (= 0).
    Book: §7.2, Eqs. (7.11), (7.12), (7.16)–(7.21); §7.3 (7.53)–(7.55); Exercise 7.2 (this is its quantitative content).
    Parameters: x [m] (array); t [s]; a, k, H, g, rho, sigma, direction as :func:`wave_fields`; h Laplacian step [m]
    (default 0.01/k). Returns dict of arrays + scales dict(aomega, ag, c, g_over_k).
    Validation: V1, V3 — tests/test_ch07.py: test_free_surface_V1_linear_conditions_hold_on_a_field,
    test_free_surface_V1_exact_residuals_equal_independent_closed_forms,
    test_free_surface_V3_neglected_terms_scale_as_ka. Checks: V1 linear residuals 0 to round-off, exact residuals equal
    independent closed forms; V7 asymptotic order in ka (exact residuals ∝ (ka)² when normalised by c and g/k, ∝ ka when
    normalised by aω, ag — an order-of-the-neglected-terms check, not a grid convergence). Label: analytic.
    """
    k = float(k)
    s = float(direction)
    f0 = wave_fields(x, 0.0, t, a, k, H, g, rho, direction, sigma)
    om = f0["omega"]
    x_ = _F(x)
    th = k * x_ - s * om * _F(t)
    eta = float(a) * np.cos(th)
    eta_x = -float(a) * k * np.sin(th)
    eta_xx = -float(a) * k ** 2 * np.cos(th)
    eta_t = s * float(a) * om * np.sin(th)
    fe = wave_fields(x_, eta, t, a, k, H, g, rho, direction, sigma)
    Ce = _F(cosh_over_sinh(k, eta, H))
    C0 = _F(cosh_over_sinh(k, 0.0, H))
    phi_t_eta = -float(a) * om ** 2 / k * Ce * np.cos(th)  # ∂φ/∂t of (7.26) at z = η
    phi_t_0 = -float(a) * om ** 2 / k * C0 * np.cos(th)
    kappa = eta_xx / (1.0 + eta_x ** 2) ** 1.5  # (7.53) exact curvature
    sr = float(sigma) / float(rho)
    z_mid = 0.5 * _dz_deep(k, H) if np.isfinite(H) else -0.5 / k
    hh = 0.01 / k if h is None else float(h)
    lap = laplacian_residual(lambda X, Z: wave_fields(X, Z, t, a, k, H, g, rho, direction, sigma)["phi"],
                             x_, np.full_like(x_, z_mid), hh)
    bottom = wave_fields(x_, _dz_deep(k, H), t, a, k, H, g, rho, direction, sigma)["w"]
    return {"laplace": _S(lap),  # Eq. (7.11)
            "bottom": _S(bottom),  # Eq. (7.12)
            "kinematic_exact": _S(_F(fe["w"]) - eta_t - eta_x * _F(fe["u"])),  # Eq. (7.16)
            "kinematic_17": _S(_F(fe["w"]) - eta_t),  # Eq. (7.17)
            "kinematic_linear": _S(_F(f0["w"]) - eta_t),  # Eq. (7.18)
            "dynamic_exact": _S(phi_t_eta + 0.5 * (_F(fe["u"]) ** 2 + _F(fe["w"]) ** 2) + float(g) * eta
                                - sr * kappa),  # Eqs. (7.19)/(7.53) with (4.83)
            "dynamic_linear": _S(phi_t_0 + float(g) * eta - sr * eta_xx),  # Eq. (7.21)/(7.55)
            "scales": {"aomega": float(a) * om, "ag": float(a) * float(g), "c": om / k, "g_over_k": float(g) / k}}


def free_surface_residual_scan(ka_values, kH=1.0, g: float = G, n_x: int = 64, k: float = 1.0) -> dict:
    """How big the linearisation error is: max over a wavelength of the residuals of the exact kinematic (7.16) and
    dynamic (7.19) conditions for the linear wave, against the steepness ka at fixed k, with log–log slopes.

    kinematic_abs [m/s], dynamic_abs [m²/s²]: absolute, ∝ (ka)² at fixed k (slope 2); kinematic_rel (÷ aω), dynamic_rel
    (÷ ga): ∝ ka (slope 1) — the dropped terms are O(ka) relative to the kept ones. The linearised conditions (7.18),
    (7.21) are satisfied to round-off (kinematic_linear, dynamic_linear). Book: §7.2 (7.16)–(7.21), Exercise 7.2 — our
    diagnostic for the C02 figure and backup explainer B1.
    Parameters: ka_values (array, each ≪ 1); kH (np.inf = deep; default 1); g [m/s²]; n_x points per wavelength;
    k [rad/m] (fixed; a = ka/k). Returns dict(ka, kinematic_abs, dynamic_abs, kinematic_rel, dynamic_rel,
    kinematic_linear, dynamic_linear, slope_abs (kinematic), slope_abs_dynamic, slope_rel (kinematic),
    slope_rel_dynamic). Validation: V3 — tests/test_ch07.py: test_free_surface_V3_neglected_terms_scale_as_ka. Checks:
    V7 asymptotic order in ka — slopes 2 and 1 (± 0.15) (the test is filed under V3; it measures the order of the
    neglected terms, not a discretisation). Label: analytic.
    """
    ka = np.atleast_1d(_F(ka_values))
    H = np.inf if np.isinf(kH) else float(kH) / float(k)
    x = np.linspace(0.0, _TWO_PI / float(k), int(n_x), endpoint=False)
    out = {n: [] for n in ("kinematic_abs", "dynamic_abs", "kinematic_rel", "dynamic_rel", "kinematic_linear",
                           "dynamic_linear")}
    for e in ka:
        r = free_surface_residuals(x, 0.0, e / k, k, H, g)
        sc = r["scales"]
        ke, de = np.max(np.abs(r["kinematic_exact"])), np.max(np.abs(r["dynamic_exact"]))
        out["kinematic_abs"].append(ke)
        out["dynamic_abs"].append(de)
        out["kinematic_rel"].append(ke / sc["aomega"])
        out["dynamic_rel"].append(de / sc["ag"])
        out["kinematic_linear"].append(np.max(np.abs(r["kinematic_linear"])))
        out["dynamic_linear"].append(np.max(np.abs(r["dynamic_linear"])))
    res = {n: np.array(v) for n, v in out.items()}
    res["ka"] = ka
    if ka.size >= 2:
        lk = np.log(ka)
        res["slope_abs"] = float(np.polyfit(lk, np.log(res["kinematic_abs"]), 1)[0])
        res["slope_abs_dynamic"] = float(np.polyfit(lk, np.log(res["dynamic_abs"]), 1)[0])
        res["slope_rel"] = float(np.polyfit(lk, np.log(res["kinematic_rel"]), 1)[0])
        res["slope_rel_dynamic"] = float(np.polyfit(lk, np.log(res["dynamic_rel"]), 1)[0])
    return res


@functools.lru_cache(maxsize=1)
def surface_wave_sympy() -> dict:
    """The §7.2 derivation (curation D05, D06) re-done in sympy, each step with its residual (0 when correct). Cached.

    Steps: trial φ = f(z) sin(kx − ωt) (7.22) → Laplace (7.11) gives f″ − k²f = 0 → f = Ae^{kz} + Be^{−kz} (7.23) →
    bottom (7.12) gives B = Ae^{−2kH} (7.24) → kinematic (7.18) gives k(A − B) = ωa (7.25) → A, B → φ regrouped into
    (7.26) → u, w (7.27) → dynamic (7.21) gives ω² = gk tanh kH (7.28); with surface tension (7.55) → (7.56); ψ (7.37)
    reproduces u, w; p′ (7.31) both forms. A wrong-sign bottom constant (B = Ae^{+2kH}) is kept as a discriminating
    residual that must NOT vanish.
    Book: §7.2, Eqs. (7.22)–(7.28), (7.31), (7.37); §7.3 (7.55)–(7.56). Returns dict of sympy expressions: trial, ode,
    general (= f_general), B_over_A (e^{−2kH}), A, B, phi, u, w, psi, regroup_identity, laplace_residual,
    bottom_residual, kinematic_residual (each 0), dynamic_line (coefficient of cos(kx − ωt) in (7.21) with (7.26):
    g − (ω²/k)coth kH, times a), dispersion (ω² = gk tanh kH), omega2_capillary, ``residuals`` (dict name → simplified
    expression, all 0) and ``wrong_bottom_residual`` (non-zero).
    Validation: V2 — tests/test_ch07.py: test_potential_V2_derivation, test_wave_fields_V2_potential_relations_symbolic,
    test_capillary_V2_derivation_tension_condition. Checks: V2 (every residual is exactly 0). Label: symbolic.
    """
    x, z, t = sp.symbols("x z t", real=True)
    a, k, H, g, om, rho, sig = sp.symbols("a k H g omega rho sigma", positive=True)
    A, B = sp.symbols("A B")
    f = sp.Function("f")
    th = k * x - om * t
    trial = f(z) * sp.sin(th)  # Eq. (7.22)
    lap = sp.diff(trial, x, 2) + sp.diff(trial, z, 2)
    ode = sp.simplify(lap / sp.sin(th))  # f'' − k² f
    fsol = sp.dsolve(sp.Eq(ode, 0), f(z)).rhs
    phi_gen = (A * sp.exp(k * z) + B * sp.exp(-k * z)) * sp.sin(th)  # Eq. (7.23)
    bottom = sp.simplify(sp.diff(phi_gen, z).subs(z, -H) / sp.sin(th))
    B_sol = sp.solve(sp.Eq(bottom, 0), B)[0]  # Eq. (7.24)
    eta = a * sp.cos(th)  # Eq. (7.2)
    kin = sp.simplify((sp.diff(phi_gen, z).subs(z, 0) - sp.diff(eta, t)) / sp.sin(th))  # Eq. (7.25)
    AB = sp.solve([sp.Eq(kin, 0), sp.Eq(B, B_sol)], [A, B], dict=True)[0]
    phi = phi_gen.subs(AB)
    phi_book = a * om / k * sp.cosh(k * (z + H)) / sp.sinh(k * H) * sp.sin(th)  # Eq. (7.26)
    u_book = a * om * sp.cosh(k * (z + H)) / sp.sinh(k * H) * sp.cos(th)  # Eq. (7.27)
    w_book = a * om * sp.sinh(k * (z + H)) / sp.sinh(k * H) * sp.sin(th)
    psi_book = a * om / k * sp.sinh(k * (z + H)) / sp.sinh(k * H) * sp.cos(th)  # Eq. (7.37)
    dyn = sp.simplify((sp.diff(phi_book, t).subs(z, 0) + g * eta) / sp.cos(th))  # Eq. (7.21), coefficient of cos
    om2 = sp.solve(sp.Eq(dyn, 0), om)
    disp = g * k * sp.tanh(k * H)  # Eq. (7.28): ω²
    dyn_sig = sp.simplify((sp.diff(phi_book, t).subs(z, 0) - sig / rho * sp.diff(eta, x, 2) + g * eta)
                          / sp.cos(th))  # Eq. (7.55)
    disp_sig = k * (g + sig * k ** 2 / rho) * sp.tanh(k * H)  # Eq. (7.56)
    pp1 = -rho * sp.diff(phi_book, t)  # Eq. (7.31) first form
    pp2 = rho * g * a * sp.cosh(k * (z + H)) / sp.cosh(k * H) * sp.cos(th)  # Eq. (7.31) last form

    def zero(e):
        return sp.simplify(sp.expand(sp.simplify(e.rewrite(sp.exp))))

    res = {"ode_solution": zero(sp.diff(fsol, z, 2) - k ** 2 * fsol),
           "phi_7_26": zero(phi - phi_book),
           "u_7_27": zero(sp.diff(phi_book, x) - u_book),
           "w_7_27": zero(sp.diff(phi_book, z) - w_book),
           "laplace_7_11": zero(sp.diff(phi_book, x, 2) + sp.diff(phi_book, z, 2)),
           "bottom_7_12": zero(sp.diff(phi_book, z).subs(z, -H)),
           "kinematic_7_18": zero(sp.diff(phi_book, z).subs(z, 0) - sp.diff(eta, t)),
           "dynamic_7_21": zero(dyn.subs(om, sp.sqrt(disp))),
           "dynamic_7_55": zero(dyn_sig.subs(om, sp.sqrt(disp_sig))),
           "psi_u": zero(sp.diff(psi_book, z) - u_book),
           "psi_w": zero(-sp.diff(psi_book, x) - w_book),
           "p_prime_7_31": zero((pp1 - pp2).subs(om, sp.sqrt(disp)))}
    wrong = sp.simplify(sp.diff(phi_gen.subs(B, A * sp.exp(2 * k * H)), z).subs(z, -H) / sp.sin(th))
    return {"trial": trial, "ode": ode, "f_general": fsol, "general": fsol, "B_of_A": B_sol,
            "B_over_A": sp.simplify(B_sol / A), "A": AB[A], "B": AB[B], "phi": phi_book, "u": u_book, "w": w_book,
            "psi": psi_book, "regroup_identity": res["phi_7_26"], "laplace_residual": res["laplace_7_11"],
            "bottom_residual": res["bottom_7_12"], "kinematic_residual": res["kinematic_7_18"], "dynamic_line": dyn,
            "omega_solutions": om2, "omega2": disp, "dispersion": sp.Eq(om ** 2, disp),
            "omega2_capillary": disp_sig, "residuals": res, "wrong_bottom_residual": wrong}


# ======================================================================================================================
# §7.2 particle orbits, streamlines, energy
# ======================================================================================================================
def _vel_fns(a, k, H, g):
    om = float(omega_gravity(k, H, g))

    def uw(xp, zp, t):
        th = k * xp - om * t
        return (a * om * _F(cosh_over_sinh(k, zp, H)) * np.cos(th),  # Eq. (7.33)
                a * om * _F(sinh_over_sinh(k, zp, H)) * np.sin(th))

    def grads(xp, zp, t):
        th = k * xp - om * t
        C = _F(cosh_over_sinh(k, zp, H))
        Sz = _F(sinh_over_sinh(k, zp, H))
        return (-a * om * k * C * np.sin(th), a * om * k * Sz * np.cos(th),  # u_x, u_z
                a * om * k * Sz * np.cos(th), a * om * k * C * np.sin(th))  # w_x, w_z
    return om, uw, grads


def orbit_linear(x0, z0, t, a: float = 0.1, k: float = 1.0, H=np.inf, g: float = G) -> tuple:
    """Linearised particle excursions ξ = −a[cosh k(z₀ + H)/sinh kH] sin(kx₀ − ωt), ζ = a[sinh k(z₀ + H)/sinh kH]
    cos(kx₀ − ωt) — closed ellipses traversed clockwise (right-going wave).

    Book: §7.2, Eqs. (7.35a, b) (integrals of (7.34a, b)); deep (7.46) circles of radius ae^{kz₀}; shallow (7.50).
    Parameters: x0, z0 mean position [m]; t [s]; a [m]; k [rad/m]; H [m]; g. Returns the tuple (ξ, ζ) [m] (ζ = the
    particle's vertical excursion, ``zeta_particle``); the particle is at (x₀ + ξ, z₀ + ζ).
    Validation: V1, V2 — tests/test_ch07.py: test_orbits_V2_derivation,
    test_orbit_V1_ellipses_clockwise_and_constant_foci, test_particle_path_V1_linear_model_and_pathline_parity. Checks:
    V1 the ellipse (7.36) = 1 at all t; clockwise (negative signed area); equals ``particle_path(model="linear")``.
    Label: analytic.
    """
    k = float(k)
    om = float(omega_gravity(k, H, g))
    th = k * _F(x0) - om * _F(t)
    xi = -float(a) * _F(cosh_over_sinh(k, z0, H)) * np.sin(th)  # Eq. (7.35a)
    ze = float(a) * _F(sinh_over_sinh(k, z0, H)) * np.cos(th)  # Eq. (7.35b)
    return _S(xi), _S(ze)


def orbit_semi_axes(z0, a: float = 0.1, k: float = 1.0, H=np.inf) -> dict:
    """Semi-axes of the particle orbit ellipse (7.36): A = a cosh k(z₀ + H)/sinh kH (horizontal, semi-major),
    B = a sinh k(z₀ + H)/sinh kH (vertical); focal half-distance √(A² − B²) = a/sinh kH, the same at every depth
    (0 in deep water: circles); B = 0 at the bottom; sense "clockwise" for a right-going wave.

    Book: §7.2, Eq. (7.36) and the text after it; Fig. 7.4. Parameters: z0 mean depth [m]; a [m]; k [rad/m]; H [m].
    Returns dict(A, B, focal_half, sense = "cw"). Validation: V1 — tests/test_ch07.py:
    test_orbit_V1_ellipses_clockwise_and_constant_foci, test_orbit_state_V1_explainer_numbers. Checks: V1 focal distance
    independent of z₀; deep A = B = ae^{kz₀}; shallow A → a/kH, B → a(1 + z₀/H). Label: analytic.
    """
    k = float(k)
    A = float(a) * _F(cosh_over_sinh(k, z0, H))  # Eq. (7.36)
    B = float(a) * _F(sinh_over_sinh(k, z0, H))
    focal = 0.0 if np.isinf(H) else float(a) / np.sinh(k * float(H)) if k * float(H) < 700 else 0.0
    return {"A": _S(A), "B": _S(B), "focal_half": float(focal), "sense": "cw"}


def particle_path(x0, z0, t_eval, a: float = 0.1, k: float = 1.0, H=np.inf, g: float = G, model: str = "exact",
                  start: str = "orbit", rtol: float = 1e-10, atol: float = 1e-12) -> dict:
    """Path of a fluid particle under the linear wave (7.27) by integrating the path-line equations.

    model = "exact": dx_p/dt = u(x_p, z_p, t), dz_p/dt = w(x_p, z_p, t) (7.32)/(7.33) — the orbits do not close: Stokes
    drift; "linear": right sides frozen at the mean position (x₀, z₀) (7.34a, b) — closed ellipses (7.35); "taylor1":
    velocity Taylor-expanded to first order in ξ = x_p − x₀, ζ = z_p − z₀ (7.84a, b).
    start = "orbit": the particle starts on its linear orbit (x₀ + ξ(0), z₀ + ζ(0)), so its mean position is ≈ (x₀, z₀);
    "mean": it starts at (x₀, z₀) itself.
    Book: §7.2, Eqs. (7.32)–(7.35); §7.6, Eqs. (7.84a, b); Figs. 7.3, 7.22.
    Parameters: x0, z0 [m] (scalars); t_eval (n,) [s] increasing from the start time; a, k, H, g; rtol, atol
    (``solve_ivp`` DOP853 — our choice). Returns dict(t, x, z) [s, m, m].
    Validation: V1, V3 — tests/test_ch07.py: test_particle_path_V1_linear_model_and_pathline_parity,
    test_stokes_drift_V3_exact_path_lines_converge_at_order_ka, test_dyed_line_V1_advances_by_the_stokes_drift. Checks:
    V1 "linear" = :func:`orbit_linear` (1e-10); V7 "exact" drift per period → (7.86) with relative error O(ka)
    (asymptotic order); V3 halving rtol changes the drift < 1e-7; "taylor1" drift = (7.86) within 3 %; V1 the dyed line
    advances by ū_L T per period. Label: analytic, converged (the rtol study only).
    """
    k = float(k)
    a = float(a)
    om, uw, grads = _vel_fns(a, k, H, g)
    t_ = _F(t_eval)
    if start == "orbit":
        xi0, ze0 = orbit_linear(x0, z0, t_[0], a, k, H, g)
        y0 = [float(x0) + float(xi0), float(z0) + float(ze0)]
    elif start == "mean":
        y0 = [float(x0), float(z0)]
    else:
        raise ValueError("start must be 'orbit' or 'mean'")
    X0, Z0 = float(x0), float(z0)

    if model == "exact":
        def rhs(t, y):
            u, w = uw(y[0], y[1], t)
            return [float(u), float(w)]  # Eq. (7.33)
    elif model == "linear":
        def rhs(t, y):
            u, w = uw(X0, Z0, t)
            return [float(u), float(w)]  # Eq. (7.34a, b)
    elif model == "taylor1":
        def rhs(t, y):
            u, w = uw(X0, Z0, t)
            ux, uz, wx, wz = grads(X0, Z0, t)
            xi, ze = y[0] - X0, y[1] - Z0
            return [float(u + xi * ux + ze * uz), float(w + xi * wx + ze * wz)]  # Eq. (7.84a, b)
    else:
        raise ValueError("model must be 'exact', 'linear' or 'taylor1'")
    sol = solve_ivp(rhs, (t_[0], t_[-1]), y0, method="DOP853", t_eval=t_, rtol=rtol, atol=atol)
    return {"t": sol.t, "x": sol.y[0], "z": sol.y[1]}


def stokes_drift(z0, a: float = 0.1, k: float = 1.0, H=np.inf, g: float = G):
    """Stokes drift (mean Lagrangian velocity) ū_L = a²ωk cosh 2k(z₀ + H)/(2 sinh² kH) [m/s]; deep water a²ωk e^{2kz₀}.

    Book: §7.6, Eqs. (7.85) (deep, from the time average of (7.84a) with (7.46)–(7.47)) and (7.86) (any depth,
    Exercise 7.14). Decays twice as fast as the orbits; the vertical drift is zero. Overflow-safe form
    a²ωk e^{2kz₀}(1 + e^{−4k(z₀+H)})/(1 − e^{−2kH})².
    Parameters: z0 mean depth [m]; a [m]; k [rad/m]; H [m]; g. Returns ū_L [m/s].
    Validation: V1, V2, V3, V5 — tests/test_ch07.py: test_stokes_drift_V2_derivation,
    test_stokes_drift_V1_closed_form_and_limits, test_stokes_drift_V3_exact_path_lines_converge_at_order_ka,
    test_stokes_drift_V5_deep_water_published. Checks: V2 sympy time average of (7.84a); V1 deep limit; V3 vs
    :func:`stokes_drift_numeric`; V5 Wikipedia "Stokes drift" deep form. Label: analytic, symbolic.
    """
    k = float(k)
    om = float(omega_gravity(k, H, g))
    z_ = _F(z0)
    H_ = float(H)
    with np.errstate(over="ignore", invalid="ignore"):
        ratio = np.exp(2 * k * z_) * (1.0 + np.exp(-4.0 * k * (z_ + H_))) / np.expm1(-2.0 * k * H_) ** 2
    return _S(float(a) ** 2 * om * k * ratio)  # Eq. (7.86) (→ (7.85) as H → ∞)


def stokes_drift_numeric(z0, a: float = 0.1, k: float = 1.0, H=np.inf, g: float = G, periods: int = 20,
                         x0: float = 0.0, rtol: float = 1e-10, atol: float = 1e-12) -> float:
    """Stokes drift measured from exact path lines: (x_p(NT) − x_p(0))/(NT) [m/s], the particle starting on its linear
    orbit about (x₀, z₀). Book: §7.6, Fig. 7.22 (the open orbit advancing ū_L T per period) — our numerical check of
    (7.86). Parameters: as :func:`particle_path`; periods N. Validation: V3 — tests/test_ch07.py:
    test_stokes_drift_V3_exact_path_lines_converge_at_order_ka. Checks: V7 → (7.86) with relative error O(ka)
    (asymptotic order in ka, not a convergence study); V3 integrator-converged (halving rtol changes the drift < 1e-7).
    Label: analytic, converged (the rtol study only)."""
    om = float(omega_gravity(k, H, g))
    T = _TWO_PI / om
    p = particle_path(x0, z0, [0.0, periods * T], a, k, H, g, "exact", "orbit", rtol, atol)
    return float((p["x"][-1] - p["x"][0]) / (periods * T))


def eulerian_mean_u(z, a: float = 0.1, k: float = 1.0, H=np.inf, g: float = G, x: float = 0.0) -> float:
    """Period-averaged horizontal velocity at a fixed point always under water, (1/T)∫₀ᵀ u dt [m/s] — zero for the periodic
    irrotational wave (the Eulerian mean), although the Lagrangian mean (Stokes drift) is not.

    Book: §7.6, text after (7.86) (u = u|_{z=−H} + ∫∂w/∂x dz by irrotationality). ``quad`` over one period.
    Validation: V1 — tests/test_ch07.py: test_stokes_drift_V1_closed_form_and_limits. Checks: V1 = 0 to 1e-10. Label:
    analytic."""
    om = float(omega_gravity(k, H, g))
    T = _TWO_PI / om
    val, _ = quad(lambda tt: float(wave_fields(x, z, tt, a, k, H, g)["u"]), 0.0, T, epsabs=1e-14, limit=200)
    return val / T


def wave_energy(a: float = 1.0, k: float = 1.0, H=np.inf, g: float = G, rho: float = 1000.0,
                method: str = "closed") -> dict:
    """Wave energy per unit horizontal area: kinetic E_k, potential E_p and total E [J/m²].

    "closed": E_k = ½ρg⟨η²⟩ (7.39) = E_p (7.41) (equipartition) = ¼ρga²; E = ½ρga² (7.42).
    "quad": E_k = (ρ/2λ)∫₀^λ∫_{−H}^0(u² + w²)dz dx (7.38) and E_p = (ρg/λ)∫₀^λ∫₀^η z dz dx (7.40) by ``dblquad`` over one
    wavelength at t = 0 (deep water: z ≥ −40/k).
    Book: §7.2, Eqs. (7.38)–(7.42). Parameters: a [m]; k [rad/m]; H [m]; g; rho [kg/m³]; method.
    Validation: V1, V2, V7 — tests/test_ch07.py: test_wave_energy_V2_derivation,
    test_wave_energy_V1_quadrature_equals_closed_form, test_wave_functions_V7_change_of_units. Checks: V1 quad = closed
    (1e-8) for kH ∈ {0.1, 1, 10, ∞}; V2 sympy z-integral. Label: analytic.
    """
    if method == "closed":
        Ek = 0.25 * float(rho) * float(g) * float(a) ** 2  # Eq. (7.39) with ⟨η²⟩ = a²/2
        Ep = 0.25 * float(rho) * float(g) * float(a) ** 2  # Eq. (7.41)
        return {"Ek": Ek, "Ep": Ep, "E": Ek + Ep}  # Eq. (7.42)
    if method != "quad":
        raise ValueError("method must be 'closed' or 'quad'")
    k = float(k)
    lam = _TWO_PI / k
    zb = _dz_deep(k, H)

    def ke(zz, xx):
        f = wave_fields(xx, zz, 0.0, a, k, H, g, rho)
        return 0.5 * float(rho) * (float(f["u"]) ** 2 + float(f["w"]) ** 2)

    Ek, _ = integrate.dblquad(ke, 0.0, lam, lambda xx: zb, lambda xx: 0.0, epsabs=1e-12, epsrel=1e-11)  # (7.38)
    Ep, _ = integrate.dblquad(lambda zz, xx: float(rho) * float(g) * zz, 0.0, lam, lambda xx: 0.0,
                              lambda xx: float(a) * np.cos(k * xx), epsabs=1e-12, epsrel=1e-11)  # Eq. (7.40)
    return {"Ek": Ek / lam, "Ep": Ep / lam, "E": (Ek + Ep) / lam}


def energy_flux(a: float = 1.0, k: float = 1.0, H=np.inf, g: float = G, rho: float = 1000.0,
                method: str = "closed") -> float:
    """Time-averaged energy flux per unit crest length across a vertical plane, F = ⟨∫_{−H}^0 p′u dz⟩ [W/m].

    "closed": F = [½ρga²][(c/2)(1 + 2kH/sinh 2kH)] (7.44) = E c_g (7.71); "quad": the double integral of p′u from
    :func:`wave_fields` over one period and the depth at x = 0 (7.43) (deep: z ≥ −40/k).
    Book: §7.2, Eqs. (7.43)–(7.44); §7.5 (7.71). (The book's "u from (7.28)" before (7.44) means (7.27).)
    Validation: V1, V2 — tests/test_ch07.py: test_energy_flux_V2_derivation,
    test_energy_flux_V1_quadrature_equals_closed_form. Checks: V1 quad = closed (1e-8); V2 sympy (7.44) ≡ E c_g. Label:
    analytic, symbolic.
    """
    k = float(k)
    if method == "closed":
        c = float(phase_speed(k, H, g))
        kH = np.inf if np.isinf(H) else k * float(H)
        fac = 1.0 + (0.0 if np.isinf(kH) else float(W._x_over_sinh(2 * kH)))
        return 0.5 * float(rho) * float(g) * float(a) ** 2 * 0.5 * c * fac  # Eq. (7.44)
    if method != "quad":
        raise ValueError("method must be 'closed' or 'quad'")
    om = float(omega_gravity(k, H, g))
    T = _TWO_PI / om
    zb = _dz_deep(k, H)

    def pu(zz, tt):
        f = wave_fields(0.0, zz, tt, a, k, H, g, rho)
        return float(f["p_prime"]) * float(f["u"])

    F, _ = integrate.dblquad(pu, 0.0, T, lambda tt: zb, lambda tt: 0.0, epsabs=1e-10, epsrel=1e-11)  # Eq. (7.43)
    return F / T


# ======================================================================================================================
# §7.3 surface tension
# ======================================================================================================================
def curvature(eta_x, eta_xx, linear: bool = False):
    """Curvature 1/R of the surface graph z = η(x): η_xx/(1 + η_x²)^{3/2} [1/m], or its small-slope form η_xx.

    Book: §7.3, Eq. (7.53) (positive when the surface is concave up — centre of curvature above, in the air).
    Parameters: eta_x [–]; eta_xx [1/m]; linear. Validation: V1 — tests/test_ch07.py:
    test_curvature_V1_circle_and_laplace_jump_parity. Checks: V1 circle of radius R → 1/R. Label: analytic.
    """
    if linear:
        return _S(_F(eta_xx))
    return _S(_F(eta_xx) / (1.0 + _F(eta_x) ** 2) ** 1.5)  # Eq. (7.53)


def capillary_surface_pressure(eta_xx, sigma: float = 0.0727, eta_x=None, p_a: float = 0.0):
    """Pressure just inside the liquid at the surface with surface tension, p = p_a − σ/R [Pa] (gauge when p_a = 0):
    p = −σ ∂²η/∂x² for small slopes (7.54); exact curvature when eta_x is given (7.53). Under a crest (η_xx < 0) the
    liquid pressure exceeds p_a (the centre of curvature is in the liquid). Book: §7.3, Eqs. (7.53)–(7.54), (1.5).
    Parameters: eta_xx [1/m]; sigma [N/m]; eta_x [–] or None; p_a [Pa]. Validation: V1, V2 — tests/test_ch07.py:
    test_capillary_V2_derivation_tension_condition, test_curvature_V1_circle_and_laplace_jump_parity. Checks: V1 parity
    with ``core.interfaces.laplace_jump_from_balance`` (one infinite radius). Label: analytic."""
    kap = curvature(0.0 if eta_x is None else eta_x, eta_xx, linear=eta_x is None)
    return _S(float(p_a) - float(sigma) * _F(kap))  # Eq. (7.53)/(7.54)


def capillary_state(lam, sigma: float = 0.0727, rho: float = 1000.0, H=np.inf, g: float = G) -> dict:
    """Explainer E3 state: everything about a capillary–gravity wave of wavelength λ.

    Returns dict(k, omega, c, cg, gravity_term = g/k [m²/s²], tension_term = σk/ρ [m²/s²], tension_ratio = (σk/ρ)/(g/k),
    regime ∈ {"capillary" (ratio > 2), "crossover", "gravity" (ratio < 0.5)}, c_min, lam_m, k_m, cg_min,
    p_crest_per_a = σk² [Pa/m] (liquid pressure excess under a crest per metre of amplitude, (7.54)), cg_over_c).
    Book: §7.3, Eqs. (7.54), (7.56)–(7.60); §7.5 (c_g, c_g,min). Scalar-callable. Validation: V1 — tests/test_ch07.py:
    test_capillary_state_V1_explainer_numbers. Label: analytic.
    """
    k = _TWO_PI / float(lam)
    gt, tt = float(g) / k, float(sigma) * k / float(rho)
    c = float(phase_speed(k, H, g, sigma, rho))
    cg = float(group_velocity(k, H, g, sigma, rho))
    if sigma > 0:
        m = W.capillary_minimum(sigma, rho, g)
        cgm = W.min_group_velocity(sigma, rho, g)["cg_min"] if g > 0 else np.nan
    else:
        m, cgm = {"c_min": np.nan, "lam_m": np.nan, "k_m": np.nan}, np.nan
    ratio = tt / gt if gt > 0 else np.inf
    regime = "capillary" if ratio > 2.0 else ("gravity" if ratio < 0.5 else "crossover")
    return {"k": k, "omega": c * k, "c": c, "cg": cg, "gravity_term": gt, "tension_term": tt, "tension_ratio": ratio,
            "regime": regime, "c_min": m["c_min"], "lam_m": m["lam_m"], "k_m": m["k_m"], "cg_min": cgm,
            "p_crest_per_a": float(sigma) * k ** 2, "cg_over_c": cg / c}


# ======================================================================================================================
# §7.4 standing waves and seiches
# ======================================================================================================================
def standing_wave_fields(x, z, t, a: float = 0.1, k: float = 1.0, H=np.inf, g: float = G,
                         rho: float = 1000.0) -> dict:
    """Standing wave = sum of equal right- and left-going waves: η = 2a cos kx cos ωt, nodes at kx = ±π/2, ±3π/2, …

    ψ = (2aω/k)[sinh k(z + H)/sinh kH] sin kx sin ωt (7.62); u = 2aω[cosh k(z + H)/sinh kH] sin kx sin ωt (7.63);
    w = −∂ψ/∂x = −2aω[sinh k(z + H)/sinh kH] cos kx sin ωt; φ = −(2aω/k)[cosh k(z + H)/sinh kH] cos kx sin ωt;
    p′ = 2ρga[cosh k(z + H)/cosh kH] cos kx cos ωt.
    Book: §7.4, text before (7.62), Eqs. (7.62)–(7.63), Fig. 7.11. Parameters as :func:`wave_fields`.
    Returns dict(eta, psi, u, w, phi, p_prime, omega). Validation: V1, V2 — tests/test_ch07.py:
    test_standing_wave_V2_derivation, test_standing_wave_V1_sum_of_opposite_waves. Checks: V1 equals the sum of
    ``wave_fields(direction=+1)`` and ``(direction=−1)`` (1e-13); nodes fixed. Label: analytic.
    """
    k = float(k)
    om = float(omega_gravity(k, H, g))
    x_, z_, t_ = np.broadcast_arrays(_F(x), _F(z), _F(t))
    C = _F(cosh_over_sinh(k, z_, H))
    Sz = _F(sinh_over_sinh(k, z_, H))
    P = _F(cosh_over_cosh(k, z_, H))
    ckx, skx, cwt, swt = np.cos(k * x_), np.sin(k * x_), np.cos(om * t_), np.sin(om * t_)
    a = float(a)
    return {"eta": _S(2 * a * ckx * cwt),  # η = 2a cos kx cos ωt
            "psi": _S(2 * a * om / k * Sz * skx * swt),  # Eq. (7.62)
            "u": _S(2 * a * om * C * skx * swt),  # Eq. (7.63)
            "w": _S(-2 * a * om * Sz * ckx * swt),
            "phi": _S(-2 * a * om / k * C * ckx * swt),
            "p_prime": _S(2 * float(rho) * float(g) * a * P * ckx * cwt), "omega": om}


def seiche_modes(L: float, H: float, n=0, g: float = G) -> dict:
    """Natural modes of a closed basin (lake) of length L and uniform depth H with vertical walls at x = 0, L.

    u = 0 at both walls ⇒ sin kL = 0 ⇒ kL = (n + 1)π, λ = 2L/(n + 1) (7.64); ω = √((πg(n + 1)/L) tanh((n + 1)πH/L))
    (7.65) = (7.28) at that k. Book: §7.4, Eqs. (7.63)–(7.65), Fig. 7.12.
    Parameters: L [m]; H [m]; n mode number 0, 1, 2, … (array allowed); g. Returns dict(k [rad/m], lam [m], omega
    [rad/s], T [s]). Validation: V1, V3 — tests/test_ch07.py: test_seiche_V1_walls_modes_and_dispersion,
    test_seiche_V3_finite_difference_sloshing_eigenproblem. Checks: V1 u(0) = u(L) = 0; V3 FD shallow-water eigenproblem
    → these frequencies. Label: analytic.
    """
    n_ = _F(n)
    k = (n_ + 1.0) * np.pi / float(L)  # kL = (n + 1)π
    om = np.sqrt(np.pi * float(g) * (n_ + 1.0) / float(L) * np.tanh((n_ + 1.0) * np.pi * float(H) / float(L)))  # (7.65)
    return {"k": _S(k), "lam": _S(2.0 * float(L) / (n_ + 1.0)),  # Eq. (7.64)
            "omega": _S(om), "T": _S(_TWO_PI / om)}


def basin_modes(L: float, b: float, H: float, m=1, n=0, g: float = G) -> dict:
    """Modes of a rectangular basin L × b (depth H, vertical walls): k² = (mπ/L)² + (nπ/b)², ω = √(gk tanh kH)
    (Exercises 7.5–7.6; our statement). η ∝ cos(mπx/L) cos(nπy/b). m, n ≥ 0, not both 0.
    Returns dict(k, lam = 2π/k, omega, T). Book: §7.4 (7.65) generalised. Validation: V1 — tests/test_ch07.py:
    test_seiche_V1_walls_modes_and_dispersion. Label: analytic."""
    m_, n_ = _F(m), _F(n)
    k = np.sqrt((m_ * np.pi / float(L)) ** 2 + (n_ * np.pi / float(b)) ** 2)
    om = _F(omega_gravity(k, H, g))
    return {"k": _S(k), "lam": _S(_TWO_PI / k), "omega": _S(om), "T": _S(_TWO_PI / om)}


def seiche_state(L: float, H: float, n: int = 0, g: float = G) -> dict:
    """Explainer E4 state: seiche mode n with the shallow-water period T_s = 2L/((n + 1)√(gH)) and its error.
    Returns dict(k, lam, kH, omega, T [s], T_shallow [s], shallow_error = T_s/T − 1, T_min = T/60 [min]).
    Book: §7.4, (7.64)–(7.65), (7.49). Scalar-callable. Validation: V7 — tests/test_ch07.py:
    test_seiche_state_V7_shallow_limit. Label: analytic."""
    s = seiche_modes(L, H, n, g)
    Ts = 2.0 * float(L) / ((float(n) + 1.0) * np.sqrt(float(g) * float(H)))
    return {**{kk: float(v) for kk, v in s.items()}, "kH": float(s["k"]) * float(H), "T_shallow": float(Ts),
            "shallow_error": float(Ts / float(s["T"]) - 1.0), "T_min": float(s["T"]) / 60.0}


def packet_state(k0: float, H=np.inf, g: float = G, sigma: float = 0.0, rho: float = 1000.0, a: float = 1.0,
                 distance: float = 1.0e6, group_length: float | None = None) -> dict:
    """Explainer E5 / C09 worked-number state for a wave group of carrier wavenumber k₀.

    Returns dict(omega, lam, c, cg, ratio = c_g/c, E = ½ρga² [J/m²] (7.42), F = E c_g [W/m] (7.71), t_crest_cross =
    group_length/|c − c_g| [s] (time for a crest to run through a group of that length; inf when c_g = c),
    arrival_h = distance/c_g [h] (energy arrival from a storm), regime ∈ {"cg<c", "cg=c", "cg>c"} (±0.5 %)).
    group_length defaults to 10λ. Book: §7.5, Eqs. (7.67)–(7.71), (7.42). Scalar-callable. Validation: V1 —
    tests/test_ch07.py: test_packet_state_V1_explainer_numbers. Label: analytic.
    """
    k = float(k0)
    om = float(omega_capillary_gravity(k, H, sigma, rho, g))
    c = om / k
    cg = float(group_velocity(k, H, g, sigma, rho))
    lam = _TWO_PI / k
    Lg = 10.0 * lam if group_length is None else float(group_length)
    E = float(W.wave_energy_density(a, rho, g))
    r = cg / c
    regime = "cg=c" if abs(r - 1.0) <= 5e-3 else ("cg<c" if r < 1.0 else "cg>c")
    return {"omega": om, "lam": lam, "c": c, "cg": cg, "ratio": r, "E": E, "F": E * cg,  # Eq. (7.71)
            "t_crest_cross": Lg / abs(c - cg) if abs(c - cg) > 0 else np.inf, "arrival_h": float(distance) / cg / 3600.0,
            "regime": regime}


def pond_ripples(x, t, width: float = 0.01, n_modes: int = 256, k_max: float | None = None, H=np.inf,
                 sigma: float = 0.0727, rho: float = 1000.0, g: float = G, a: float = 1.0):
    """Stone-in-a-pond (Cauchy–Poisson) wave train from a released Gaussian hump η₀ = a e^{−x²/(2w²)} (at rest), as a
    direct cosine sum η(x, t) = Σ_j A_j cos(k_j x) cos(ω(k_j)t) with A_j = (1/π)F(k_j)Δk, F(k) = a w√(2π) e^{−k²w²/2} the
    cosine transform of the hump, k_j = (j + ½)Δk, Δk = k_max/n (midpoint rule; k_max defaults to 6/w). ω from (7.56)
    (capillary–gravity: the calm centre inside c_g,min·t, the fastest long waves in front). The sum is periodic in x
    with period 2π/Δk — keep |x| well inside it. An explainer can mirror it term by term (no FFT).
    Book: §7.5, Fig. 7.16 and the stone-in-a-pond paragraph (N72) — our construction (curation "Cauchy–Poisson problem,
    our extension"). Parameters: x [m] (array); t [s] scalar or (M,); width w [m]; n_modes; k_max [rad/m]; H; sigma
    [N/m]; rho; g; a [m]. Returns η [m] (shape of x, or (M, len(x))).
    Validation: V3 — tests/test_ch07.py: test_pond_ripples_V3_cosine_sum_converges_to_the_cauchy_poisson_integral.
    Checks: V1 t = 0 reproduces the hump (midpoint sum → integral); V1 parity with ``linear_evolve``. Label: analytic.
    """
    w_ = float(width)
    km = 6.0 / w_ if k_max is None else float(k_max)
    dk = km / int(n_modes)
    kj = (np.arange(int(n_modes)) + 0.5) * dk
    Aj = float(a) * w_ * np.sqrt(_TWO_PI) * np.exp(-0.5 * (kj * w_) ** 2) * dk / np.pi  # cosine-transform weights
    wj = _F(omega_capillary_gravity(kj, H, sigma, rho, g))
    x_ = _F(x)
    t_ = _F(t)
    cosx = np.cos(np.multiply.outer(x_, kj))  # (…, n)
    if t_.ndim == 0:
        return _S(cosx @ (Aj * np.cos(wj * float(t_))))
    return np.stack([cosx @ (Aj * np.cos(wj * tt)) for tt in t_])


# ======================================================================================================================
# §7.5 group velocity, kinematic wave theory, refraction
# ======================================================================================================================
def _d1(f: Callable, x, h):
    return (-f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h) + f(x - 2 * h)) / (12.0 * h)


def local_wavenumber_frequency(theta_fn: Callable, x, t, h: float = 1e-4, ht: float | None = None) -> tuple:
    """Local wavenumber and frequency of a slowly varying wave train η = a(x, t) cos θ(x, t) (7.72):
    k ≡ ∂θ/∂x, ω ≡ −∂θ/∂t (7.73), by fourth-order central differences.

    Book: §7.5, Eqs. (7.72)–(7.73). Parameters: theta_fn θ(x, t) [rad]; x [m]; t [s]; h [m], ht [s] (default h) steps
    (choose ≪ the wavelength and period). Returns the tuple (k [rad/m], ω [rad/s]). Validation: V3 — tests/test_ch07.py:
    test_local_wavenumber_V3_fourth_order_and_crest_conservation. Checks: V1 θ = kx − ωt exact; V3 order 4. Label:
    analytic, converged."""
    x_, t_ = _F(x), _F(t)
    ht = h if ht is None else ht
    k = _d1(lambda xx: _F(theta_fn(xx, t_)), x_, h)  # Eq. (7.73)
    om = -_d1(lambda tt: _F(theta_fn(x_, tt)), t_, ht)
    return _S(k), _S(om)


def crest_conservation_residual(theta_fn: Callable, x, t, H_fn: Callable | None = None, h: float = 1e-3,
                                ht: float | None = None, omega_fn: Callable | None = None, g: float = G) -> dict:
    """Residuals of kinematic wave theory for a phase field θ(x, t):
    r774 = ∂k/∂t + ∂ω/∂x (7.74, crest conservation — identically 0 for any smooth θ);
    with a dispersion relation (surface gravity waves on depth H_fn(x), or ``omega_fn`` ω(k) for a homogeneous medium):
    r_dispersion = ω_local − ω(k_local, x); r775 = ∂k/∂t + c_g ∂k/∂x (7.75, homogeneous media only);
    r779 = ∂ω/∂t + c_g ∂ω/∂x (7.79, frequency constant along rays, also for H = H(x)).
    Book: §7.5, Eqs. (7.73)–(7.79). Parameters: theta_fn; x [m]; t [s]; H_fn or omega_fn; h [m], ht [s] (default h)
    steps; g. Returns dict(r774 [, r775, r779, r_dispersion], k, omega [, cg]). Validation: V2, V3 — tests/test_ch07.py:
    test_crest_conservation_V2_derivation, test_local_wavenumber_V3_fourth_order_and_crest_conservation. Checks: V1
    plane wave; V3 order 4 in h. Label: converged.
    """
    x_, t_ = _F(x), _F(t)
    ht = h if ht is None else ht

    def k_of(xx, tt):
        return _d1(lambda q: _F(theta_fn(q, tt)), xx, h)

    def w_of(xx, tt):
        return -_d1(lambda q: _F(theta_fn(xx, q)), tt, ht)

    k = k_of(x_, t_)
    om = w_of(x_, t_)
    k_t = _d1(lambda q: k_of(x_, q), t_, ht)
    k_x = _d1(lambda q: k_of(q, t_), x_, h)
    w_t = _d1(lambda q: w_of(x_, q), t_, ht)
    w_x = _d1(lambda q: w_of(q, t_), x_, h)
    out = {"k": _S(k), "omega": _S(om), "r774": _S(k_t + w_x)}  # Eq. (7.74)
    if H_fn is not None:
        Hx = _F(H_fn(x_))
        cg = _F(group_velocity(k, Hx, g))
        wdisp = _F(omega_gravity(k, Hx, g))  # Eq. (7.76)
    elif omega_fn is not None:
        cg = _F(group_velocity_numeric(omega_fn, k))
        wdisp = _F(omega_fn(k))
    else:
        return out
    out.update({"cg": _S(cg), "r_dispersion": _S(om - wdisp),
                "r775": _S(k_t + cg * k_x),  # Eq. (7.75)
                "r779": _S(w_t + cg * w_x)})  # Eq. (7.79)
    return out


def snell_ray_plane_beach(x, alpha0: float, x0: float, T: float, slope: float, g: float = G) -> dict:
    """Closed-form ray over a plane beach H(x) = slope·x (shore at x = 0, x offshore): with the full dispersion relation
    (7.28) and ω = 2π/T fixed along the ray (7.79), the alongshore wavenumber is conserved, k(x) sin α(x) = k₀ sin α₀
    (Snell's law; α = angle between the ray/wavenumber and the offshore-to-shore normal −e_x), and the ray is
    y(x) = ∫_x^{x₀} tan α(x′) dx′.

    Book: §7.2 refraction (Fig. 7.8, words only) and §7.5 (7.79) — the Snell form is ours (curation D24).
    Parameters: x positions ≤ x₀ (array, > 0) [m]; alpha0 incidence angle at x₀ [rad]; x0 start [m]; T period [s];
    slope dH/dx [–]; g. Returns dict(x, y, H, k, alpha, alpha_deg, c, cg, l (alongshore wavenumber), snell = k sin α
    (= l at every point)) along the ray.
    Validation: V1, V2 — tests/test_ch07.py: test_snell_V2_derivation, test_ray_trace_V1_snell_closed_form_parity.
    Checks: V1 parity with :func:`ray_trace` over the same beach; k sin α constant. Label: analytic.
    """
    om = _TWO_PI / float(T)
    x_ = np.atleast_1d(_F(x))
    Hs = float(slope) * x_
    k = np.array([float(wavenumber_from_omega(om, Hi, g)) for Hi in Hs])  # ω fixed along the ray, (7.28) inverted
    k0 = float(wavenumber_from_omega(om, float(slope) * float(x0), g))
    l = k0 * np.sin(float(alpha0))  # Snell: k sin α conserved
    alpha = np.arcsin(np.clip(l / k, -1, 1))

    def tan_a(xx):
        kk = float(wavenumber_from_omega(om, float(slope) * xx, g))
        s = l / kk
        return s / np.sqrt(1.0 - s * s)

    y = np.array([quad(tan_a, xi, float(x0), epsabs=1e-12, epsrel=1e-12, limit=200)[0] for xi in x_])
    return {"x": _S(x_), "y": _S(y), "H": _S(Hs), "k": _S(k), "alpha": _S(alpha), "alpha_deg": _S(np.degrees(alpha)),
            "c": _S(om / k), "cg": _S(_F(group_velocity(k, Hs, g))), "l": float(l),
            "snell": _S(k * np.sin(alpha))}  # Snell: k sin α = const


def refraction_state(T: float, alpha0: float, H0: float, H: float, g: float = G) -> dict:
    """Explainer E6 state: one point of a plane-beach ray. A wave of period T arrives from depth H₀ at angle α₀ to the
    depth-contour normal; at depth H the wavenumber follows from ω = 2π/T fixed (7.28) and the angle from Snell's law
    k sin α = k₀ sin α₀ (curation D24, ours). Returns dict(omega, k0, k, alpha [rad], alpha_deg, c, cg, snell = k sin α,
    lam). Book: §7.2 refraction (Fig. 7.8); §7.5 (7.79). Scalar-callable. Validation: V1 — tests/test_ch07.py:
    test_refraction_state_V1_worked_example. Label: analytic."""
    om = _TWO_PI / float(T)
    k0 = float(wavenumber_from_omega(om, H0, g))
    k = float(wavenumber_from_omega(om, H, g))
    s = k0 * np.sin(float(alpha0)) / k
    alpha = float(np.arcsin(np.clip(s, -1.0, 1.0)))
    return {"omega": om, "k0": k0, "k": k, "alpha": alpha, "alpha_deg": float(np.degrees(alpha)), "c": om / k,
            "cg": float(group_velocity(k, H, g)), "snell": k * np.sin(alpha), "lam": _TWO_PI / k}


def wave_regime_label(kH) -> str:
    """"deep" (kH > 2), "intermediate" or "shallow" (H/λ < 0.07) — the book's thresholds (§7.2, after (7.45) and
    (7.49)); thin wrapper of :func:`~fluidpy.core.waves.depth_regime` with k = 1. Scalar-callable. Validation: V1 —
    tests/test_ch07.py: test_depth_regime_V1_book_thresholds_and_errors. Label: analytic."""
    return W.depth_regime(1.0, float(kH))["regime"]


def dispersion_state(lam: float, H=np.inf, g: float = G, rho: float = 1000.0, a: float = 1.0) -> dict:
    """Explainer E1 state for a wave of wavelength λ on depth H: dict(k, kH, H_over_lambda, omega, T, c, cg, c_deep =
    √(g/k), c_shallow = √(gH), regime, deep_error, shallow_error, p_bottom_fraction = 1/cosh kH, p_surface_amp = ρga [Pa],
    t_cross_1000km_c_h, t_cross_1000km_cg_h = 10⁶ m / c (or c_g) in hours).
    Book: §7.2, (7.28), (7.29), (7.31), (7.45), (7.49); §7.5 (7.69). Scalar-callable. Validation: V1 —
    tests/test_ch07.py: test_dispersion_state_V1_explainer_numbers. Label: analytic."""
    k = _TWO_PI / float(lam)
    om = float(omega_gravity(k, H, g))
    c = om / k
    cg = float(group_velocity(k, H, g))
    reg = W.depth_regime(k, H)
    return {"k": k, "kH": float(reg["kH"]), "H_over_lambda": float(reg["H_over_lambda"]), "omega": om,
            "T": _TWO_PI / om, "c": c, "cg": cg, "c_deep": float(np.sqrt(float(g) / k)),
            "c_shallow": float(np.sqrt(float(g) * float(H))), "regime": reg["regime"],
            "deep_error": float(reg["deep_error"]), "shallow_error": float(reg["shallow_error"]),
            "p_bottom_fraction": 0.0 if np.isinf(H) else float(pressure_response(k, -float(H), H)),
            "p_surface_amp": float(rho) * float(g) * float(a), "t_cross_1000km_c_h": 1e6 / c / 3600.0,
            "t_cross_1000km_cg_h": 1e6 / cg / 3600.0}


# ======================================================================================================================
# §7.6 nonlinear waves: steepening, hydraulic jump, Stokes waves, KdV
# ======================================================================================================================
def nonlinear_wavelet_speed(eta, H: float, g: float = G, model: str = "simple"):
    """Speed of a small wavelet riding on a finite-amplitude shallow-water wave [m/s].

    model = "book" (§7.6, text with Fig. 7.19): c = c′ + u with c′ = √(gH′), H′ = H + η the local depth, u the particle
    speed — here u from linear shallow-water theory u = (c₀/H)η ((7.51) with (7.52)): c = √(g(H + η)) + c₀η/H.
    model = "simple" (Riemann simple wave, **our extension**, exact for the nonlinear shallow-water equations):
    c = 3√(g(H + η)) − 2c₀. Both = c₀(1 + (3/2)η/H) + O(η²) — the KdV coefficient of (7.87). Crests outrun troughs.
    Parameters: eta [m]; H [m]; g; model. Returns c [m/s]. Validation: V1 — tests/test_ch07.py:
    test_simple_wave_V1_breaking_time_and_first_order_speed. Checks: V1 first-order agreement. Label: analytic.
    """
    e = _F(eta)
    c0 = np.sqrt(float(g) * float(H))
    if model == "book":
        return _S(np.sqrt(float(g) * (float(H) + e)) + c0 * e / float(H))  # c = c′ + u
    if model == "simple":
        return _S(3.0 * np.sqrt(float(g) * (float(H) + e)) - 2.0 * c0)
    raise ValueError("model must be 'simple' or 'book'")


def simple_wave_evolve(eta0, x, t, H: float, g: float = G) -> dict:
    """Nonlinear steepening of a right-going shallow-water hump by characteristics: each surface point keeps its η and
    moves at c(η) = 3√(g(H + η)) − 2√(gH) (Riemann simple wave), so x_i(t) = x_i + c(η₀_i) t; the profile becomes
    multivalued (breaks) at t_b = 1/max(−dc/dx) (for a sine a sin kx: t_b = 2H/(3akc₀) to leading order in a/H — the
    steepest point moves off η = 0 at the next order because c′(η) ∝ (H + η)^{−1/2}).

    Book: §7.6, Fig. 7.19 (qualitative: "crests overtake troughs", method of characteristics named) — this simple-wave
    solution is **our extension** (labelled). Parameters: eta0 (N,) [m] at x (N,) [m]; t scalar or (M,) [s]; H; g.
    Returns dict(x_points (M, N) or (N,) positions of the surface points (also under the key ``x``), eta (N,)
    (unchanged along characteristics), t_break [s], c (N,)). Validation: V1 — tests/test_ch07.py:
    test_simple_wave_V1_breaking_time_and_first_order_speed. Checks: V1 breaking time of a sine vs the formula. Label:
    analytic.
    """
    e0 = _F(eta0)
    x_ = _F(x)
    c = 3.0 * np.sqrt(float(g) * (float(H) + e0)) - 2.0 * np.sqrt(float(g) * float(H))
    dcdx = np.gradient(c, x_)
    tb = 1.0 / np.max(-dcdx) if np.max(-dcdx) > 0 else np.inf
    t_ = _F(t)
    xs = x_[None, :] + c[None, :] * np.atleast_1d(t_)[:, None]
    xp = xs[0] if t_.ndim == 0 else xs
    return {"x_points": xp, "x": xp, "eta": e0, "t_break": float(tb), "c": c}


def jump_momentum_residual(H1, H2, Q, g: float = G):
    """Momentum balance of a stationary hydraulic jump per unit width and density,
    Q²(1/H₂ − 1/H₁) − ½g(H₁² − H₂²) [m³/s²] (0 for a physical pair of conjugate depths).
    Book: §7.6, Eq. (7.80) (CV (4.17) with d/dt = 0, b = 0, uniform velocity, hydrostatic faces; Fig. 7.20b).
    Validation: V1, V2 — tests/test_ch07.py: test_hydraulic_jump_V2_derivation_belanger,
    test_hydraulic_jump_V1_momentum_mass_and_numbers. Label: analytic."""
    return _S(_F(Q) ** 2 * (1.0 / _F(H2) - 1.0 / _F(H1)) - 0.5 * float(g) * (_F(H1) ** 2 - _F(H2) ** 2))  # Eq. (7.80)


def hydraulic_jump(H1: float, u1: float | None = None, Fr1: float | None = None, g: float = G) -> dict:
    """Conjugate (Bélanger) depth, speeds and energy loss of a stationary hydraulic jump.

    H₂/H₁ = ½(−1 + √(1 + 8Fr₁²)) (7.81), Fr₁² = Q²/gH₁³ = u₁²/gH₁ (4.104); mass Q = u₁H₁ = u₂H₂; energy of a surface
    particle E = u²/2 + gH, E₂ − E₁ = −(H₂ − H₁)g(H₂ − H₁)²/(4H₁H₂) (text after (7.81)): a loss (head loss
    (H₂ − H₁)³/(4H₁H₂)) only if H₂ > H₁ ⇔ Fr₁ > 1 — Fr₁ < 1 would create mechanical energy (second law forbids).
    Book: §7.6, Eqs. (7.80)–(7.81) and the energy line; Fig. 7.20.
    Parameters: H1 upstream depth [m]; u1 [m/s] or Fr1 [–]; g. Returns dict(H2, u1, u2, Fr1, Fr2, ratio, Q [m²/s],
    E1, E2, dE [J/kg], head_loss [m], allowed (Fr₁ ≥ 1)). Validation: V1, V2, V4 — tests/test_ch07.py:
    test_hydraulic_jump_V2_derivation_belanger, test_hydraulic_jump_V2_derivation_energy_loss,
    test_hydraulic_jump_V1_momentum_mass_and_numbers, test_hydraulic_jump_V4_second_law,
    test_hydraulic_jump_V1_form_belanger_published. Checks: V2 sympy (7.80) ⇒ quadratic and the energy identity; V1
    momentum residual 0; V1 same *form* as Wikipedia's Bélanger equation (a form check, not a numerical benchmark); V4
    second law; parity with ``ch04.bore_speed``. Label: analytic, symbolic.
    """
    H1 = float(H1)
    if (u1 is None) == (Fr1 is None):
        raise ValueError("give exactly one of u1 or Fr1")
    if u1 is None:
        u1 = float(Fr1) * np.sqrt(float(g) * H1)
    u1 = float(u1)
    Fr1v = u1 / np.sqrt(float(g) * H1)  # Eq. (4.104)
    r = 0.5 * (-1.0 + np.sqrt(1.0 + 8.0 * Fr1v ** 2))  # Eq. (7.81)
    H2 = r * H1
    Q = u1 * H1
    u2 = Q / H2
    E1 = 0.5 * u1 ** 2 + float(g) * H1
    E2 = 0.5 * u2 ** 2 + float(g) * H2
    dE = -(H2 - H1) * float(g) * (H2 - H1) ** 2 / (4.0 * H1 * H2)  # E₂ − E₁ (text after (7.81))
    return {"H2": H2, "u1": u1, "u2": u2, "Fr1": Fr1v, "Fr2": u2 / np.sqrt(float(g) * H2), "ratio": r, "Q": Q,
            "E1": E1, "E2": E2, "dE": dE, "head_loss": -dE / float(g), "allowed": bool(Fr1v >= 1.0)}


def jump_state(H1: float, Fr1: float, g: float = G, rho: float = 1000.0, frame: str = "stationary") -> dict:
    """Explainer E7 state: the hydraulic jump with its momentum budget and the moving-bore frame (per metre of width).

    Momentum budget of the CV of Fig. 7.20b (7.80): mom_in = ρQu₁ (inflow of momentum), mom_out = ρQu₂, p_in = ½ρgH₁²
    (pressure force on face 1, pushing +x), p_out = ½ρgH₂² (on face 2, pushing −x); residual = (mom_in − mom_out) +
    (p_in − p_out) = 0 [N/m]. power_loss = ρgQ·head_loss [W/m]. frame = "bore": the same jump moving at
    bore_speed = u₁ = √(gH₂(H₁ + H₂)/2H₁) into still water of depth H₁ (Fig. 7.20c; = ``ch04.bore_speed(H1, H2)``); the
    water behind follows at u₁ − u₂ (flow_behind).
    Book: §7.6, (7.80)–(7.81) and the jump energy line; §4.4 Example 4.3. Returns dict(H2, ratio, u1, u2, Q, Fr2, dE,
    head_loss, power_loss, mom_in, mom_out, p_in, p_out, residual, allowed, bore_speed, flow_behind, frame).
    Scalar-callable. Validation: V1 — tests/test_ch07.py: test_jump_state_V1_budget_and_moving_bore_parity. Label:
    analytic.
    """
    j = hydraulic_jump(H1, Fr1=Fr1, g=g)
    r = float(rho)
    mi, mo = r * j["Q"] * j["u1"], r * j["Q"] * j["u2"]
    pi_, po = 0.5 * r * float(g) * float(H1) ** 2, 0.5 * r * float(g) * j["H2"] ** 2
    return {"H2": j["H2"], "ratio": j["ratio"], "u1": j["u1"], "u2": j["u2"], "Q": j["Q"], "Fr2": j["Fr2"],
            "dE": j["dE"], "head_loss": j["head_loss"], "power_loss": r * float(g) * j["Q"] * j["head_loss"],
            "mom_in": mi, "mom_out": mo, "p_in": pi_, "p_out": po,
            "residual": (mi - mo) + (pi_ - po),  # Eq. (7.80) × ρ: ρQ(u₂ − u₁) = ½ρg(H₁² − H₂²)
            "allowed": j["allowed"], "bore_speed": j["u1"], "flow_behind": j["u1"] - j["u2"], "frame": frame}


def stokes_wave_speed(k, a, g: float = G):
    """Amplitude-dependent speed of a deep-water Stokes wave c = √((g/k)(1 + k²a²)) [m/s] (to O((ka)²)).

    Book: §7.6, Eq. (7.83) (Stokes 1847). The coefficient γ = 1 in c² = (g/k)(1 + γk²a²) comes from the **consistent
    third-order** Stokes expansion (:func:`stokes_expansion_sympy` ``third_order``), in which the potential's first
    harmonic carries a correction −(5/8)(ka)² relative to aω/k (−(1/8)(ka)² relative to a√(g/k)), fixed only by the
    O((ka)³) kinematic balance. The literal Exercise 7.2 ansatz (potential amplitude held at aω/k, kinematic condition
    truncated at (ka)¹) gives γ = 3/8 instead (``exercise_literal``) — it is not a derivation of (7.83).
    Parameters: k [rad/m] > 0; a [m]; g [m/s²]. Returns c [m/s].
    Validation: V1, V2, V5 — tests/test_ch07.py: test_stokes_expansion_V2_third_order_coefficients,
    test_stokes_wave_V5_speed_and_limiting_steepness, test_stokes_wave_profile_V1_harmonics_and_permanence. Checks: V2
    consistent third-order expansion (γ = 1); V5 HandWiki "Stokes wave" c = (1 + ½(ka)²)√(g/k). Label: analytic."""
    return _S(np.sqrt(float(g) / _F(k) * (1.0 + (_F(k) * _F(a)) ** 2)))  # Eq. (7.83)


def stokes_wave_profile(x, t, a: float = 0.1, k: float = 1.0, g: float = G, order: int = 3):
    """Stokes wave to third order η = a cos k(x − ct) + ½ka² cos 2k(x − ct) + (3/8)k²a³ cos 3k(x − ct) [m], c from (7.83):
    peaked crests, flat troughs; all harmonics move at one speed, so the profile is permanent. Asymptotic only — keep
    ka ≲ 0.3; the 120° limiting crest is not a property of the truncated series.
    Book: §7.6, Eq. (7.82), Fig. 7.21. Parameters: x [m]; t [s]; a [m]; k [rad/m]; g; order 1–3.
    Validation: V1, V2 — tests/test_ch07.py: test_stokes_expansion_V2_third_order_coefficients,
    test_stokes_wave_profile_V1_harmonics_and_permanence. Checks: V2 the coefficients ½ and 3/8 from the consistent
    third-order expansion (:func:`stokes_expansion_sympy`); V1 permanence and harmonics. Label: analytic."""
    k = float(k)
    c = float(stokes_wave_speed(k, a, g))
    th = k * (_F(x) - c * _F(t))
    eta = float(a) * np.cos(th)
    if order >= 2:
        eta = eta + 0.5 * k * float(a) ** 2 * np.cos(2 * th)
    if order >= 3:
        eta = eta + 0.375 * k ** 2 * float(a) ** 3 * np.cos(3 * th)  # Eq. (7.82)
    return _S(eta)


@functools.lru_cache(maxsize=1)
def stokes_expansion_sympy() -> dict:
    """The deep-water Stokes expansion (Exercise 7.2; curation's demoted a-D40) in sympy, non-dimensional (k = g = 1,
    ε = ka), on z = η with the exact kinematic condition (7.16) and the full-Bernoulli dynamic condition p = 0 (the
    Bernoulli constant absorbs the mean), each projected on its Fourier harmonics.

    * ``third_order`` — η = ε cos θ + αε² cos 2θ + δε³ cos 3θ, φ = εω(1 + βε²)e^{z} sin θ, ω² = 1 + γε²: the kinematic
      condition at O(ε²) (sin 2θ) gives α = ½; at O(ε³) (sin θ, sin 3θ) β = −5/8, δ = 3/8; the dynamic condition at O(ε³)
      (cos θ) then gives γ = 1 — the book's (7.82) coefficients ½ and 3/8 and (7.83) c² = (g/k)(1 + k²a²).
    * ``exercise_literal`` — the exercise's truncated set-up with the potential amplitude fixed at aω/k (β = 0) and the
      kinematic condition only to (ka)¹: α = ½, but the O(ε³) dynamic cos θ balance gives γ = 3/8, not 1. γ = 1 needs the
      O((ka)²) correction of the potential's first harmonic (β = −5/8), which the kinematic condition fixes only at the
      next order (our finding, confirmed by the derivation review M2: (7.83) is right, the truncated ansatz is not).
    Book: §7.6, Eqs. (7.82)–(7.83); Exercise 7.2. Returns dict(alpha (= 1/2), gamma (= 1) — the consistent third-order
    values, third_order=dict(alpha, beta, delta, gamma), exercise_literal=dict(alpha, gamma)). Cached. Validation: V2 —
    tests/test_ch07.py: test_stokes_expansion_V2_third_order_coefficients. Label: symbolic.
    """
    th, eps, z = sp.symbols("theta epsilon z", real=True)
    al, be, de, ga = sp.symbols("alpha beta delta gamma")

    from sympy.simplify.fu import TR8  # products of sines and cosines → sums of harmonics

    def proj(expr, fn, n):
        return sp.expand(TR8(sp.expand(expr))).coeff(fn(n * th))  # Fourier coefficient of fn(nθ)

    def conditions(eta, amp):
        # travelling wave: every field depends on θ = x − ωt, so ∂/∂x = d/dθ and ∂/∂t = −ω d/dθ
        om = sp.sqrt(1 + ga * eps ** 2)  # Eq. (7.83), k = g = 1
        phi = eps * om * amp * sp.exp(z) * sp.sin(th)
        phi_x, phi_z, phi_t = sp.diff(phi, th), sp.diff(phi, z), -om * sp.diff(phi, th)
        eta_x, eta_t = sp.diff(eta, th), -om * sp.diff(eta, th)
        kin = (phi_z - eta_t - eta_x * phi_x).subs(z, eta)  # Eq. (7.16)
        dyn = (phi_t + (phi_x ** 2 + phi_z ** 2) / 2).subs(z, eta) + eta  # p = 0 at z = η, full Bernoulli
        return (sp.expand(sp.series(kin, eps, 0, 4).removeO()), sp.expand(sp.series(dyn, eps, 0, 4).removeO()))

    kin, dyn = conditions(eps * sp.cos(th) + al * eps ** 2 * sp.cos(2 * th) + de * eps ** 3 * sp.cos(3 * th),
                          1 + be * eps ** 2)
    a_sol = sp.solve(proj(kin.coeff(eps, 2), sp.sin, 2), al)[0]
    k3 = kin.coeff(eps, 3).subs(al, a_sol)
    b_sol = sp.solve(proj(k3, sp.sin, 1), be)[0]
    d_sol = sp.solve(proj(k3, sp.sin, 3), de)[0]
    g_sol = sp.solve(proj(dyn.coeff(eps, 3).subs({al: a_sol, be: b_sol, de: d_sol}), sp.cos, 1), ga)[0]
    kinL, dynL = conditions(eps * sp.cos(th) + al * eps ** 2 * sp.cos(2 * th), 1)
    aL = sp.solve(proj(kinL.coeff(eps, 2), sp.sin, 2), al)[0]
    gL = sp.solve(proj(dynL.coeff(eps, 3).subs(al, aL), sp.cos, 1), ga)[0]
    return {"alpha": a_sol, "gamma": g_sol,
            "third_order": {"alpha": a_sol, "beta": b_sol, "delta": d_sol, "gamma": g_sol},
            "exercise_literal": {"alpha": aL, "gamma": gL}}


def solitary_wave_speed(a: float, H: float, g: float = G) -> float:
    """Solitary-wave speed c = c₀(1 + a/2H), c₀ = √(gH) [m/s] — taller is faster. Book: §7.6, after (7.88).
    Validation: V1, V2 — tests/test_ch07.py: test_kdv_V2_solitary_wave_residual,
    test_kdv_solve_V1_soliton_translates_at_c0_1_plus_a_over_2H, test_cnoidal_V1_form_and_soliton_limit. Label:
    analytic."""
    return float(np.sqrt(float(g) * float(H)) * (1.0 + float(a) / (2.0 * float(H))))


def solitary_wave(x, t, a: float = 0.2, H: float = 1.0, g: float = G, x0: float = 0.0):
    """Solitary wave (soliton) of the KdV equation η = a sech²[(3a/4H³)^{1/2}(x − x₀ − ct)], c = c₀(1 + a/2H) [m].
    Book: §7.6, Eq. (7.88), Fig. 7.23b (Russell 1844; check by substitution, Exercise 7.15 → :func:`kdv_residual_sympy`).
    Parameters: x [m]; t [s]; a crest height [m]; H depth [m]; g; x0 [m]. Validation: V1, V2 — tests/test_ch07.py:
    test_kdv_V2_solitary_wave_residual, test_kdv_solve_V1_soliton_translates_at_c0_1_plus_a_over_2H,
    test_cnoidal_V1_form_and_soliton_limit. Checks: V2 residual 0; V1 KdV solver translates it at c; V5 Wikipedia
    "Cnoidal wave" solitary limit. Label: analytic, symbolic."""
    c = solitary_wave_speed(a, H, g)
    arg = np.sqrt(3.0 * float(a) / (4.0 * float(H) ** 3)) * (_F(x) - float(x0) - c * _F(t))
    return _S(float(a) / np.cosh(arg) ** 2)  # Eq. (7.88)


def cnoidal_wave(x, t, H: float = 1.0, height: float = 0.2, m: float = 0.9, g: float = G,
                 datum: str = "mean") -> dict:
    """Cnoidal wave, the periodic permanent-form KdV solution η = η₂ + A cn²((x − ct)/Δ | m) with A = height (crest to
    trough), Δ = H√(4mH/(3A)), c = c₀[1 + 3η₂/(2H) + A(2m − 1)/(2mH)], wavelength λ = 2K(m)Δ. ``datum="mean"``: trough
    η₂ = −A(E/K − 1 + m)/m (zero mean over a wavelength; approaches the soliton only like 1/K(m), logarithmically);
    ``datum="trough"``: η₂ = 0 (trough at the still-water level) — then m → 1 gives exactly the solitary wave (7.88).
    Book: §7.6, text before (7.88) and Fig. 7.23a (named; not printed) — our derivation from (7.87) with
    (cn²)″ = −6m cn⁴ + 4(2m − 1)cn² + 2(1 − m); cross-checked with Wikipedia "Cnoidal wave". scipy ``ellipj`` (m = k²).
    Parameters: x [m]; t [s]; H depth [m]; height A [m]; m elliptic parameter 0 < m < 1; g.
    Returns dict(eta, c, wavelength, trough, Delta). Validation: V1 — tests/test_ch07.py:
    test_cnoidal_V1_form_and_soliton_limit. Checks: V1 m → 1 limit vs :func:`solitary_wave` (1e-6); V1 KdV residual
    (spectral) ~ 0. Label: analytic.
    """
    A = float(height)
    Hh = float(H)
    m = float(m)
    c0 = np.sqrt(float(g) * Hh)
    Delta = Hh * np.sqrt(4.0 * m * Hh / (3.0 * A))
    Km, Em = ellipk(m), ellipe(m)
    if datum not in ("mean", "trough"):
        raise ValueError("datum must be 'mean' or 'trough'")
    eta2 = -A * (Em / Km - 1.0 + m) / m if datum == "mean" else 0.0
    c = c0 * (1.0 + 1.5 * eta2 / Hh + A * (2.0 * m - 1.0) / (2.0 * m * Hh))
    _, cn, _, _ = ellipj((_F(x) - c * _F(t)) / Delta, m)
    return {"eta": _S(eta2 + A * cn ** 2), "c": float(c), "wavelength": float(2.0 * Km * Delta), "trough": float(eta2),
            "Delta": float(Delta)}


def ursell_number(a, lam, H):
    """Ursell ratio aλ²/H³ [–] of nonlinear to dispersive terms in KdV (7.87): ≳ 16 → steepening to a jump; smaller →
    balance and permanent forms (cnoidal, solitary). Book: §7.6, text after (7.87) (the book says "(7.88)"; the terms are
    those of (7.87)). Validation: V1 — tests/test_ch07.py: test_kdv_linear_phase_speed_V1_taylor_of_7_29. Label:
    analytic."""
    return _S(_F(a) * _F(lam) ** 2 / _F(H) ** 3)


def kdv_linear_phase_speed(k, H: float, g: float = G):
    """Phase speed of the linearised KdV equation c = c₀(1 − k²H²/6) [m/s] — the first two Taylor terms of (7.29).
    Book: §7.6, text after (7.87). Validation: V1 — tests/test_ch07.py: test_kdv_linear_phase_speed_V1_taylor_of_7_29.
    Checks: V1 error vs (7.29) ∝ (kH)⁴. Label: analytic."""
    return _S(np.sqrt(float(g) * float(H)) * (1.0 - _F(k) ** 2 * float(H) ** 2 / 6.0))


def _kdv_setup(x, H, g):
    x_ = _F(x)
    N = x_.size
    dx = x_[1] - x_[0]
    kx = _TWO_PI * np.fft.fftfreq(N, d=dx)
    c0 = np.sqrt(float(g) * float(H))
    return x_, N, dx, kx, c0


def kdv_rhs(eta, x, H: float, g: float = G):
    """Right side of the KdV equation, ∂η/∂t = −c₀η_x − (3/2)(c₀/H)ηη_x − (1/6)c₀H²η_xxx [m/s], by spectral derivatives
    on a periodic grid. Book: §7.6, Eq. (7.87). Parameters: eta (N,) [m]; x (N,) uniform periodic [m]; H; g.
    Validation: V1 — tests/test_ch07.py: test_kdv_solve_V1_soliton_translates_at_c0_1_plus_a_over_2H. Label:
    analytic."""
    x_, N, dx, kx, c0 = _kdv_setup(x, H, g)
    E = np.fft.fft(_F(eta))
    ex = np.real(np.fft.ifft(1j * kx * E))
    exxx = np.real(np.fft.ifft((1j * kx) ** 3 * E))
    e = _F(eta)
    return -c0 * ex - 1.5 * c0 / float(H) * e * ex - c0 * float(H) ** 2 / 6.0 * exxx  # Eq. (7.87)


def kdv_invariants(eta, x, H: float | None = None, g: float = G) -> dict:
    """Conserved integrals of KdV (7.87) on a periodic domain (rectangle rule = spectrally accurate for periodic data):
    mass ∫η dx [m²]; momentum ∫η² dx [m³] (the KdV "momentum", ∝ the wave potential energy ½ρg∫η²); energy = the
    Hamiltonian ∫(αη³/6 − βη_x²/2) dx [m⁴/s] with α = (3/2)c₀/H, β = c₀H²/6 (η_x spectral; NaN unless H is given — the
    c₀η_x term only adds total derivatives). Book: §7.6 (7.87) — our diagnostic. Works on (N,) or (M, N) arrays.
    Validation: V4 — tests/test_ch07.py: test_kdv_solve_V4_invariants_of_a_splitting_hump. Checks: V4 all three
    conserved by ``kdv_solve``. Label: conserved."""
    e = _F(eta)
    x_ = _F(x)
    dx = float(x_[1] - x_[0])
    out = {"mass": _S(np.sum(e, axis=-1) * dx), "momentum": _S(np.sum(e ** 2, axis=-1) * dx)}
    if H is None:
        out["energy"] = _S(np.full(e.shape[:-1], np.nan)) if e.ndim > 1 else float("nan")
        return out
    c0 = np.sqrt(float(g) * float(H))
    al, be = 1.5 * c0 / float(H), c0 * float(H) ** 2 / 6.0
    kx = _TWO_PI * np.fft.fftfreq(x_.size, d=dx)
    ex = np.real(np.fft.ifft(1j * kx * np.fft.fft(e, axis=-1), axis=-1))
    out["energy"] = _S(np.sum(al * e ** 3 / 6.0 - be * ex ** 2 / 2.0, axis=-1) * dx)
    return out


def kdv_solve(eta0, x, t_out, H: float, g: float = G, dt: float | None = None, dealias: bool = True,
              cache: str | Path | None = None) -> dict:
    """Solve the KdV equation (7.87) on a periodic domain: Fourier pseudo-spectral in x, integrating-factor RK4 in t.

    The stiff linear part L(k) = −i(c₀k − c₀H²k³/6) is integrated exactly (integrating factor e^{Lt}); the nonlinear term
    −(3/4)(c₀/H) ∂(η²)/∂x by classical RK4 (Trefethen, *Spectral Methods in MATLAB*, p27 scheme), with 2/3-rule
    dealiasing of the product. **Our scheme** (the book states (7.87) but gives no numerical method).
    Book: §7.6, Eq. (7.87); Figs. 7.23 (solitary waves), the Ursell balance.
    Parameters: eta0 (N,) [m] on a uniform periodic x (N,) [m]; t_out increasing output times [s] (first ≥ 0); H depth
    [m]; g; dt [s] (default 0.1 dx/c_max, c_max = c₀(1 + 1.5 max|η₀|/H) — stricter than the analysis's 0.5 for
    accuracy, see the comment in the code — then shortened to land on each output time);
    dealias; cache: path of an ``.npz`` file (e.g. outputs/ch07/kdv_run.npz) — reused when it was written with the
    same inputs (SHA-1 of η₀, x, t_out, H, g, dt, dealias; any change recomputes and overwrites).
    Returns dict(t, eta (M, N), mass, momentum, energy (M,) — :func:`kdv_invariants`, dt, steps, cached (bool)).
    Validation: V1, V3, V4 — tests/test_ch07.py: test_kdv_solve_V1_soliton_translates_at_c0_1_plus_a_over_2H,
    test_kdv_solve_V3_time_order_four_and_spectral_space, test_kdv_solve_V4_invariants_of_a_splitting_hump. Checks: V1
    the solitary wave (7.88) translates at c₀(1 + a/2H); V3 order 4 in dt, spectral in N; V4 mass, momentum, energy
    conserved. Label: analytic, converged, conserved.
    """
    if cache is not None:
        key = hashlib.sha1()
        for arr in (_F(eta0), _F(x), np.atleast_1d(_F(t_out)), np.array([float(H), float(g), -1.0 if dt is None
                                                                         else float(dt), float(dealias)])):
            key.update(np.ascontiguousarray(arr).tobytes())
        digest = key.hexdigest()
        cpath = Path(cache)
        if cpath.exists():
            with np.load(cpath, allow_pickle=False) as d:
                if str(d["digest"]) == digest:
                    return {"t": d["t"], "eta": d["eta"], "mass": d["mass"], "momentum": d["momentum"],
                            "energy": d["energy"], "dt": float(d["dt"]), "steps": int(d["steps"]), "cached": True}
        res = kdv_solve(eta0, x, t_out, H, g, dt, dealias, cache=None)
        cpath.parent.mkdir(parents=True, exist_ok=True)
        np.savez(cpath, digest=digest, **{kk: v for kk, v in res.items() if kk != "cached"})
        return res
    x_, N, dx, kx, c0 = _kdv_setup(x, H, g)
    Hf = float(H)
    Lop = -1j * (c0 * kx - c0 * Hf ** 2 * kx ** 3 / 6.0)  # linear part of (7.87) in Fourier space
    mask = (np.abs(kx) <= (2.0 / 3.0) * np.max(np.abs(kx))) if dealias else np.ones(N, bool)
    alpha = 1.5 * c0 / Hf

    def nl(V):
        e = np.real(np.fft.ifft(V))
        return -0.5 * alpha * 1j * kx * np.fft.fft(e * e) * mask  # −(3/2)(c₀/H)ηη_x = −(3/4)(c₀/H)(η²)_x

    t_out_ = np.atleast_1d(_F(t_out))
    e0 = _F(eta0)
    cmax = c0 * (1.0 + 1.5 * np.max(np.abs(e0)) / Hf)
    # default step: CFL 0.1 on c_max (the integrating factor makes any dt stable; the RK4 error of the nonlinear part
    # ∝ dt⁴ — with 0.5 dx/c_max the invariants of a splitting hump drift ~2e-4, with 0.1 dx/c_max ~4e-7)
    dtm = 0.1 * dx / cmax if dt is None else float(dt)
    V = np.fft.fft(e0)
    tc = 0.0
    out = []
    steps = 0
    for to in t_out_:
        span = to - tc
        if span < -1e-12:
            raise ValueError("t_out must be increasing and >= 0")
        n = int(np.ceil(span / dtm - 1e-9)) if span > 0 else 0
        if n > 0:
            h = span / n
            E1 = np.exp(Lop * h / 2.0)
            E2 = E1 * E1
            for _ in range(n):
                a_ = h * nl(V)
                b_ = h * nl(E1 * (V + a_ / 2.0))
                c_ = h * nl(E1 * V + b_ / 2.0)
                d_ = h * nl(E2 * V + E1 * c_)
                V = E2 * V + (E2 * a_ + 2.0 * E1 * (b_ + c_) + d_) / 6.0  # IF-RK4 step
            steps += n
        tc = to
        out.append(np.real(np.fft.ifft(V)))
    eta = np.array(out)
    inv = kdv_invariants(eta, x_, H, g)
    return {"t": t_out_, "eta": eta, "mass": np.atleast_1d(inv["mass"]), "momentum": np.atleast_1d(inv["momentum"]),
            "energy": np.atleast_1d(inv["energy"]), "dt": dtm, "steps": steps, "cached": False}


def kdv_residual_sympy(nonlinear_coefficient=sp.Rational(3, 2)) -> dict:
    """Exercise 7.15: substitute the solitary wave (7.88) into KdV (7.87) and simplify — the residual is exactly 0 with
    c = c₀(1 + a/2H). Done in the variable T = tanh(β(x − ct)), β = √(3a/4H³) (d/dx → β(1 − T²)d/dT), so the residual is
    a polynomial in T. ``nonlinear_coefficient`` ≠ 3/2 is the discriminating wrong variant (non-zero residual).
    Book: §7.6, Eqs. (7.87)–(7.88). Returns dict(residual, residual_coefficients, speed). Validation: V2 —
    tests/test_ch07.py: test_kdv_V2_solitary_wave_residual. Label: symbolic."""
    a, H, g, T = sp.symbols("a H g T", positive=True)
    c0 = sp.sqrt(g * H)
    c = c0 * (1 + a / (2 * H))
    beta = sp.sqrt(3 * a / (4 * H ** 3))
    eta = a * (1 - T ** 2)  # a sech² = a(1 − tanh²)

    def dx(f):
        return sp.expand(sp.diff(f, T) * beta * (1 - T ** 2))

    ex = dx(eta)
    res = -c * ex + c0 * ex + nonlinear_coefficient * c0 * eta / H * ex + sp.Rational(1, 6) * c0 * H ** 2 * dx(dx(ex))
    res = sp.expand(sp.simplify(res))
    return {"residual": sp.simplify(res), "residual_coefficients": sp.Poly(res, T).all_coeffs(), "speed": c}


# ======================================================================================================================
# §7.7 waves on a density interface
# ======================================================================================================================
def interface_fields(x, z, t, a: float = 0.1, k: float = 1.0, rho1: float = 1000.0, rho2: float = 1020.0,
                     g: float = G) -> dict:
    """Interfacial wave between two infinitely deep fluids (ρ₁ above, ρ₂ > ρ₁ below), from the complex solution
    ζ = a e^{i(kx−ωt)} (7.89), φ₁ = A e^{−kz} e^{i(kx−ωt)}, φ₂ = B e^{kz} e^{i(kx−ωt)}, A = −B = iωa/k (from (7.93)),
    ω = ε√(gk) (7.95); real parts returned.
    u₁ = −ωa e^{−kz} cos θ, u₂ = ωa e^{kz} cos θ (opposite: a vortex sheet), w₁ = ωa e^{−kz} sin θ, w₂ = ωa e^{kz} sin θ.
    Book: §7.7, Eqs. (7.89)–(7.95), the u₁, u₂ lines after (7.96), Fig. 7.24.
    Parameters: x, z [m] (interface at z = 0); t [s]; a [m]; k [rad/m]; rho1, rho2 [kg/m³]; g.
    Returns dict(zeta_interface (interface displacement; also ``zeta``), phi1, phi2, u1, u2, w1, w2, phi, u, w (fluid 1
    for z ≥ 0, fluid 2 below), A, B (complex), omega, gamma_sheet = u₂ − u₁ at z = 0 (u_below − u_above, ch05 convention)
    = 2ωa cos θ).
    Validation: V1, V2 — tests/test_ch07.py: test_interface_V2_derivation,
    test_interface_fields_V1_residuals_vortex_sheet_and_parity. Checks: V1 residuals of (7.90)–(7.94) = 0; A = −B =
    iωa/k. Label: analytic.
    """
    k = float(k)
    om = float(W.interface_omega(k, rho1, rho2, g))
    x_, z_, t_ = np.broadcast_arrays(_F(x), _F(z), _F(t))
    ph = np.exp(1j * (k * x_ - om * t_))
    A = 1j * om * float(a) / k  # from (7.93)
    B = -A
    e1, e2 = np.exp(-k * z_), np.exp(k * z_)
    phi1, phi2 = A * e1 * ph, B * e2 * ph
    u1, u2 = 1j * k * phi1, 1j * k * phi2  # u = ∂φ/∂x
    w1, w2 = -k * phi1, k * phi2  # w = ∂φ/∂z
    up = z_ >= 0
    th = k * x_ - om * t_
    zi = _S(np.real(float(a) * ph))  # Eq. (7.89)
    return {"zeta_interface": zi, "zeta": zi,
            "phi1": _S(np.real(phi1)), "phi2": _S(np.real(phi2)), "u1": _S(np.real(u1)), "u2": _S(np.real(u2)),
            "w1": _S(np.real(w1)), "w2": _S(np.real(w2)),
            "phi": _S(np.real(np.where(up, phi1, phi2))), "u": _S(np.real(np.where(up, u1, u2))),
            "w": _S(np.real(np.where(up, w1, w2))), "A": complex(A), "B": complex(B), "omega": om,
            "gamma_sheet": _S(2.0 * om * float(a) * np.cos(th))}


def interface_residuals(x, t, a: float = 0.1, k: float = 1.0, rho1: float = 1000.0, rho2: float = 1020.0,
                        g: float = G, h: float | None = None) -> dict:
    """Residuals of the two-fluid problem for :func:`interface_fields`: Laplace in each fluid (7.90) (O(h⁴) stencil at
    z = ±0.5/k), decay (7.91)–(7.92) (|φ| at z = ±40/k), kinematic (7.93) on both sides and pressure continuity (7.94) at
    z = 0 (analytic derivatives of the complex fields). Book: §7.7, Eqs. (7.90)–(7.94). Returns dict of arrays.
    Validation: V1 — tests/test_ch07.py: test_interface_fields_V1_residuals_vortex_sheet_and_parity. Checks: V1 all ~ 0.
    Label: analytic."""
    k = float(k)
    x_ = _F(x)
    hh = 0.01 / k if h is None else float(h)
    f0 = interface_fields(x_, 0.0, t, a, k, rho1, rho2, g)
    om = f0["omega"]
    A, B = f0["A"], f0["B"]
    ph = np.exp(1j * (k * x_ - om * _F(t)))
    zeta = float(a) * ph
    lap1 = laplacian_residual(lambda X, Z: interface_fields(X, Z, t, a, k, rho1, rho2, g)["phi1"], x_,
                              np.full_like(x_, 0.5 / k), hh)
    lap2 = laplacian_residual(lambda X, Z: interface_fields(X, Z, t, a, k, rho1, rho2, g)["phi2"], x_,
                              np.full_like(x_, -0.5 / k), hh)
    far1 = interface_fields(x_, 40.0 / k, t, a, k, rho1, rho2, g)["phi1"]
    far2 = interface_fields(x_, -40.0 / k, t, a, k, rho1, rho2, g)["phi2"]
    dz1, dz2 = -k * A * ph, k * B * ph
    zt = -1j * om * zeta
    dyn = float(rho1) * (-1j * om * A * ph) + float(rho1) * float(g) * zeta \
        - (float(rho2) * (-1j * om * B * ph) + float(rho2) * float(g) * zeta)  # Eq. (7.94)
    return {"laplace1": _S(lap1), "laplace2": _S(lap2), "far1": _S(far1), "far2": _S(far2),  # (7.90)–(7.92)
            "kinematic1": _S(np.real(dz1 - zt)), "kinematic2": _S(np.real(dz2 - zt)),  # Eq. (7.93)
            "dynamic": _S(np.real(dyn))}


def interface_energy(a: float = 0.1, k: float = 1.0, rho1: float = 1000.0, rho2: float = 1020.0, g: float = G,
                     method: str = "closed") -> dict:
    """Energy of an interfacial wave per unit horizontal area: E_k = ¼(ρ₂ − ρ₁)ga² (Exercise 7.18), E_p = ¼(ρ₂ − ρ₁)ga²
    (column swap, Fig. 7.25), E = ½(ρ₂ − ρ₁)ga² (7.96) [J/m²].

    "quad": E_k = (1/λ)∫₀^λ∫_{−∞}^{∞} ½ρ(u² + w²) dz dx and E_p = (1/λ)∫₀^λ∫₀^ζ (ρ₂ − ρ₁)g z dz dx (the energy of lifting
    ρ₂ above z = 0 and lowering ρ₁) by ``dblquad`` (|z| ≤ 40/k).
    Book: §7.7, E_k and E_p lines before and (7.96). (The printed middle form of E_p with ∫₀^{λ/2} and 1/2λ gives ⅛; its
    outer forms give ¼ — see the module notes.) Validation: V1 — tests/test_ch07.py:
    test_interface_energy_V1_quarter_from_direct_integration. Checks: V1 quad = closed (1e-8). Label: analytic.
    """
    drho = float(rho2) - float(rho1)
    if method == "closed":
        Ek = 0.25 * drho * float(g) * float(a) ** 2
        Ep = 0.25 * drho * float(g) * float(a) ** 2
        return {"Ek": Ek, "Ep": Ep, "E": Ek + Ep}  # Eq. (7.96)
    if method != "quad":
        raise ValueError("method must be 'closed' or 'quad'")
    k = float(k)
    lam = _TWO_PI / k
    zmax = 40.0 / k

    def ke(zz, xx):
        f = interface_fields(xx, zz, 0.0, a, k, rho1, rho2, g)
        r = float(rho1) if zz >= 0 else float(rho2)
        return 0.5 * r * (float(f["u"]) ** 2 + float(f["w"]) ** 2)

    Ek1, _ = integrate.dblquad(ke, 0.0, lam, lambda xx: 0.0, lambda xx: zmax, epsabs=1e-12, epsrel=1e-11)
    Ek2, _ = integrate.dblquad(ke, 0.0, lam, lambda xx: -zmax, lambda xx: 0.0, epsabs=1e-12, epsrel=1e-11)
    Ep, _ = integrate.dblquad(lambda zz, xx: drho * float(g) * zz, 0.0, lam, lambda xx: 0.0,
                              lambda xx: float(a) * np.cos(k * xx), epsabs=1e-12, epsrel=1e-11)
    return {"Ek": (Ek1 + Ek2) / lam, "Ep": Ep / lam, "E": (Ek1 + Ek2 + Ep) / lam}


def _two_layer_constants(k, H, om, a, g):
    """(7.106)–(7.109): A, B, C, b for η = a e^{i(kx−ωt)} (a real)."""
    A = -0.5j * a * (om / k + g / om)  # Eq. (7.106)
    B = 0.5j * a * (om / k - g / om)  # Eq. (7.107)
    C = -0.5j * a * (om / k + g / om) - 0.5j * a * (om / k - g / om) * np.exp(2 * k * H)  # Eq. (7.108)
    b = 0.5 * a * (1 + g * k / om ** 2) * np.exp(-k * H) + 0.5 * a * (1 - g * k / om ** 2) * np.exp(k * H)  # (7.109)
    return A, B, C, b


def two_layer_modes(k: float, H: float, rho1: float = 1000.0, rho2: float = 1002.0, g: float = G, a: float = 1.0,
                    mode: str = "baroclinic", long_wave: bool = False) -> dict:
    """A layer of thickness H (ρ₁, free surface on top) over an infinitely deep fluid ρ₂: the barotropic and baroclinic
    modes with their constants and amplitude ratio.

    Surface η = a e^{i(kx−ωt)} (7.102) (a real), interface at z = −H: ζ = b e^{i(kx−ωt)} (7.103); φ₁ = (Ae^{kz} + Be^{−kz})
    e^{i(kx−ωt)} (7.104), φ₂ = Ce^{kz}e^{i(kx−ωt)} (7.105, read e^{i(kx−ωt)}); A, B, C, b from (7.106)–(7.109).
    mode "barotropic": ω² = gk (7.111), b = ae^{−kH} (7.112) (in phase); "baroclinic": (7.113), η/ζ =
    −((ρ₂ − ρ₁)/ρ₁)e^{−kH} (7.114) (antiphase, interface ≫ surface). long_wave=True (kH ≪ 1, baroclinic): ω² =
    kg((ρ₂ − ρ₁)/ρ₂)kH (7.115), c = √(g′H) (7.116)–(7.117), η/ζ = −(ρ₂ − ρ₁)/ρ₁ (7.118) — the exact values are returned
    alongside. p′ in the upper layer: iρ₁ω(Ae^{kz} + Be^{−kz}) e^{i(kx−ωt)} → ρ₁gη (7.119) as kH → 0.
    Book: §7.7, Eqs. (7.97)–(7.119), Figs. 7.27–7.28. Parameters: k [rad/m]; H [m]; rho1 < rho2 [kg/m³]; g; a [m].
    Returns dict(omega, c, A, B, C, b (complex), eta_over_zeta (a/b, real), p_top, p_interface (complex p′ amplitudes at
    z = 0 and −H), p_hydrostatic = ρ₁ga, p_prime_check = dict(top_over_hydrostatic (= 1: (7.119) at z = 0),
    interface_over_top (→ 1 as kH → 0: depth-independent, hydrostatic)), mode, plus omega_exact, eta_over_zeta_exact when
    long_wave). Validation: V1, V2 — tests/test_ch07.py: test_two_layer_modes_V2_derivation,
    test_two_layer_modes_V1_amplitude_ratios_and_pressure, test_two_layer_residuals_V1_numeric_and_printed_variant.
    Checks: V2 :func:`two_layer_sympy`; V1 :func:`two_layer_residuals` = 0; V7 limits. Label: analytic.
    """
    k, H, a = float(k), float(H), float(a)
    w_bt, w_bc = W.two_layer_free_surface_omega(k, H, rho1, rho2, g)
    if mode == "barotropic":
        om = float(w_bt)
    elif mode == "baroclinic":
        om = float(w_bc)
    else:
        raise ValueError("mode must be 'barotropic' or 'baroclinic'")
    A, B, C, b = _two_layer_constants(k, H, om, a, float(g))
    p_top = complex(1j * float(rho1) * om * (A + B))  # Eq. (7.119) at z = 0
    p_int = complex(1j * float(rho1) * om * (A * np.exp(-k * H) + B * np.exp(k * H)))
    out = {"omega": om, "c": om / k, "A": complex(A), "B": complex(B), "C": complex(C), "b": complex(b),
           "eta_over_zeta": float(np.real(a / b)), "p_top": p_top, "p_interface": p_int,
           "p_hydrostatic": float(rho1) * float(g) * a,
           "p_prime_check": {"top_over_hydrostatic": float(np.real(p_top) / (float(rho1) * float(g) * a)),
                             "interface_over_top": float(np.real(p_int / p_top))}, "mode": mode}
    if mode == "barotropic":
        out["b_book"] = a * np.exp(-k * H)  # Eq. (7.112)
    else:
        out["eta_over_zeta_book"] = -(float(rho2) - float(rho1)) / float(rho1) * np.exp(-k * H)  # Eq. (7.114)
    if long_wave:
        if mode != "baroclinic":
            raise ValueError("long_wave applies to the baroclinic mode")
        gp = W.reduced_gravity_book(rho1, rho2, g, "lower")
        out["omega_exact"], out["eta_over_zeta_exact"] = out["omega"], out["eta_over_zeta"]
        out["omega"] = float(np.sqrt(k * float(g) * (float(rho2) - float(rho1)) / float(rho2) * k * H))  # (7.115)
        out["c"] = float(np.sqrt(gp * H))  # Eq. (7.116)
        out["eta_over_zeta"] = -(float(rho2) - float(rho1)) / float(rho1)  # Eq. (7.118)
    return out


def two_layer_residuals(x, t, k: float = 0.1, H: float = 50.0, rho1: float = 1000.0, rho2: float = 1002.0,
                        g: float = G, a: float = 1.0, mode: str = "baroclinic", printed_7_105: bool = False) -> dict:
    """Residuals of the layer-over-deep-fluid problem (7.90), (7.97)–(7.101) for :func:`two_layer_modes`' solution
    (real parts). ``printed_7_105=True`` uses φ₂ = Ce^{kz}e^{i(k**z** − ωt)} as printed (§9 T3) — its Laplace and interface
    residuals are then non-zero (the wrong variant). Book: §7.7, Eqs. (7.97)–(7.105). Returns dict(laplace2, decay2,
    kin_surface (7.98), dyn_surface (7.99), kin_interface_1, kin_interface_2 (7.100), dyn_interface (7.101)).
    Validation: V1 — tests/test_ch07.py: test_two_layer_residuals_V1_numeric_and_printed_variant. Checks: V1 all ~ 0 for
    the correct form. Label: analytic."""
    k = float(k)
    md = two_layer_modes(k, H, rho1, rho2, g, a, mode)
    om, A, B, C, b = md["omega"], md["A"], md["B"], md["C"], md["b"]
    x_ = _F(x)
    ph = np.exp(1j * (k * x_ - om * _F(t)))

    def phi2(z):
        if printed_7_105:
            return C * np.exp(k * z) * np.exp(1j * (k * z - om * _F(t))) + 0.0 * x_
        return C * np.exp(k * z) * ph

    def phi2_z(z):
        return ((k + 1j * k) if printed_7_105 else k) * phi2(z)

    def phi1(z):
        return (A * np.exp(k * z) + B * np.exp(-k * z)) * ph

    def phi1_z(z):
        return k * (A * np.exp(k * z) - B * np.exp(-k * z)) * ph

    zi = -float(H)
    eta = float(a) * ph
    zeta = b * ph
    # φ₂,xx + φ₂,zz: correct form (ik)² + k² = 0; printed form has no x-dependence, ∂²/∂z² = ((1 + i)k)² ≠ 0
    lap_coef = ((1 + 1j) * k) ** 2 if printed_7_105 else (1j * k) ** 2 + k ** 2
    lap2 = lap_coef * phi2(zi - 1.0 / k)  # Eq. (7.90)
    return {"laplace2": _S(np.real(lap2)),
            "decay2": _S(np.real(phi2(zi - 40.0 / k))),  # Eq. (7.97)
            "kin_surface": _S(np.real(phi1_z(0.0) + 1j * om * eta)),  # Eq. (7.98)
            "dyn_surface": _S(np.real(-1j * om * phi1(0.0) + float(g) * eta)),  # Eq. (7.99)
            "kin_interface_1": _S(np.real(phi1_z(zi) + 1j * om * zeta)),  # Eq. (7.100)
            "kin_interface_2": _S(np.real(phi2_z(zi) + 1j * om * zeta)),
            "dyn_interface": _S(np.real(float(rho1) * (-1j * om * phi1(zi)) + float(rho1) * float(g) * zeta
                                        - float(rho2) * (-1j * om * phi2(zi)) - float(rho2) * float(g) * zeta))}  # (7.101)


@functools.lru_cache(maxsize=1)
def two_layer_sympy() -> dict:
    """The §7.7 layer-over-deep-fluid algebra in sympy (curation D29–D30, ★★★): solve (7.98)–(7.100) for A, B, C, b and
    compare with (7.106)–(7.109); substitute into (7.101) and factor — the residual is (ag²k/ω²) × the left side of
    (7.110); check its two roots (7.111), (7.113), the barotropic ratio (7.112) and the baroclinic ratio (7.114).
    The printed (7.105) φ₂ = Ce^{kz}e^{i(kz − ωt)} (slip T3) fails Laplace (7.90): printed_7105_residual ≠ 0.
    Book: §7.7, Eqs. (7.98)–(7.114); Exercise 7.19. Returns dict(A, B, C, b (solved), bc_residuals (the four conditions
    (7.98), (7.99), (7.100) ×2 with the book's (7.106)–(7.109) substituted: 4 zeros), pressure_residual ((7.101) with
    them), common_factor (= ag²k/ω²), dispersion_7110 (its left side), factored (the factorised residual),
    printed_7105_residual, residuals (dict, all 0), ratio_7_101_to_7_110, omega2_baroclinic). Cached. Validation: V2 —
    tests/test_ch07.py: test_two_layer_constants_V2_derivation, test_two_layer_dispersion_V2_derivation. Label:
    symbolic."""
    a, k, H, g, om, r1, r2 = sp.symbols("a k H g omega rho1 rho2", positive=True)
    A, B, C, b = sp.symbols("A B C b")
    I = sp.I
    E, Ei = sp.exp(k * H), sp.exp(-k * H)
    eqs = [sp.Eq(k * (A - B), -I * om * a),  # (7.98): ∂φ₁/∂z = ∂η/∂t at z = 0
           sp.Eq(-I * om * (A + B) + g * a, 0),  # (7.99)
           sp.Eq(k * (A * Ei - B * E), k * C * Ei),  # (7.100) φ₁_z = φ₂_z at z = −H
           sp.Eq(k * C * Ei, -I * om * b)]  # (7.100) φ₂_z = ζ_t
    sol = sp.solve(eqs, [A, B, C, b], dict=True)[0]
    A_b = -I * a / 2 * (om / k + g / om)  # Eq. (7.106)
    B_b = I * a / 2 * (om / k - g / om)  # Eq. (7.107)
    C_b = -I * a / 2 * (om / k + g / om) - I * a / 2 * (om / k - g / om) * sp.exp(2 * k * H)  # Eq. (7.108)
    b_b = a / 2 * (1 + g * k / om ** 2) * Ei + a / 2 * (1 - g * k / om ** 2) * E  # Eq. (7.109)
    dyn = (r1 * (-I * om * (A * Ei + B * E)) + r1 * g * b - r2 * (-I * om * C * Ei) - r2 * g * b)  # (7.101)
    dyn_s = sp.simplify(dyn.subs(sol))
    lhs_710 = (om ** 2 / (g * k) - 1) * (om ** 2 / (g * k) * (r1 * sp.sinh(k * H) + r2 * sp.cosh(k * H))
                                         - (r2 - r1) * sp.sinh(k * H))  # Eq. (7.110)
    ratio = sp.simplify((dyn_s / lhs_710).rewrite(sp.exp))
    om_bc2 = g * k * (r2 - r1) * sp.sinh(k * H) / (r2 * sp.cosh(k * H) + r1 * sp.sinh(k * H))  # Eq. (7.113)

    def z0(e):
        return sp.simplify(sp.expand(e.rewrite(sp.exp)))

    res = {"A_7_106": z0(sol[A] - A_b), "B_7_107": z0(sol[B] - B_b), "C_7_108": z0(sol[C] - C_b),
           "b_7_109": z0(sol[b] - b_b),
           "root_7_111": z0(lhs_710.subs(om, sp.sqrt(g * k))),
           "root_7_113": z0(lhs_710.subs(om, sp.sqrt(om_bc2))),
           "b_7_112": z0(b_b.subs(om, sp.sqrt(g * k)) - a * Ei),
           "eta_over_zeta_7_114": z0(a / b_b.subs(om, sp.sqrt(om_bc2)) + (r2 - r1) / r1 * Ei)}
    book = {A: A_b, B: B_b, C: C_b, b: b_b}
    bc = [z0((e.lhs - e.rhs).subs(book)) for e in eqs]
    pres = z0(dyn.subs(book))
    zz, xx, tt = sp.symbols("z x t", real=True)
    phi2_printed = sp.Symbol("C") * sp.exp(k * zz) * sp.exp(I * (k * zz - om * tt))  # (7.105) as printed
    printed_res = sp.simplify((sp.diff(phi2_printed, xx, 2) + sp.diff(phi2_printed, zz, 2)) / phi2_printed)
    return {"A": sol[A], "B": sol[B], "C": sol[C], "b": sol[b], "bc_residuals": bc, "pressure_residual": pres,
            "common_factor": sp.simplify(ratio), "dispersion_7110": lhs_710,
            "factored": sp.Mul(ratio, lhs_710, evaluate=False), "printed_7105_residual": printed_res,
            "residuals": res, "ratio_7_101_to_7_110": ratio, "omega2_baroclinic": om_bc2}


def two_layer_rigid_lid_omega(k, h1: float, h2: float, rho1: float, rho2: float, g: float = G):
    """Interfacial waves between rigid lids (upper layer h₁ over lower layer h₂; Exercise 7.20, our derivation):
    ω² = gk(ρ₂ − ρ₁)/(ρ₁ coth kh₁ + ρ₂ coth kh₂) [rad/s]; h₁, h₂ → ∞ recovers (7.95); long waves c² → g(ρ₂ − ρ₁)h₁h₂/
    (ρ₁h₂ + ρ₂h₁) ≈ g′h₁h₂/(h₁ + h₂). Book: §7.7 (Exercise 7.20). Validation: V7 — tests/test_ch07.py:
    test_two_layer_V7_limits_and_reduced_gravity. Label: analytic."""
    kk = np.abs(_F(k))
    coth = lambda q: 1.0 / np.tanh(q)  # noqa: E731
    w2 = float(g) * kk * (float(rho2) - float(rho1)) / (float(rho1) * coth(kk * h1) + float(rho2) * coth(kk * h2))
    with np.errstate(invalid="ignore"):
        return _S(np.where(w2 >= 0, np.sqrt(np.maximum(w2, 0)), np.nan))


def two_layer_state(k: float, H: float, rho1: float = 1000.0, rho2: float = 1002.0, g: float = G,
                    mode: str = "baroclinic") -> dict:
    """Explainer E8 state: dict(omega_bt, omega_bc, omega (of ``mode``), c_bt, c_bc, T (period of mode), eta_over_zeta
    (of mode), g_prime_lower (7.117), g_prime_upper (ch04), c_long = √(g′H), eps2_density, omega_two_deep (7.95), kH,
    mode). Book: §7.7, (7.95), (7.110)–(7.118). Scalar-callable. Validation: V1 — tests/test_ch07.py:
    test_two_layer_state_V1_explainer_numbers. Label: analytic."""
    w_bt, w_bc = W.two_layer_free_surface_omega(k, H, rho1, rho2, g)
    md = two_layer_modes(k, H, rho1, rho2, g, 1.0, mode)
    kf = float(k)
    return {"omega_bt": float(w_bt), "omega_bc": float(w_bc), "omega": md["omega"], "c_bt": float(w_bt) / kf,
            "c_bc": float(w_bc) / kf, "T": _TWO_PI / md["omega"], "eta_over_zeta": md["eta_over_zeta"],
            "g_prime_lower": W.reduced_gravity_book(rho1, rho2, g, "lower"),
            "g_prime_upper": W.reduced_gravity_book(rho1, rho2, g, "upper"),
            "c_long": float(W.two_layer_long_wave_speed(H, rho1, rho2, g)),
            "eps2_density": float(W.eps2_density(rho1, rho2)),
            "omega_two_deep": float(W.interface_omega(k, rho1, rho2, g)), "kH": kf * float(H), "mode": mode}


# ======================================================================================================================
# §7.8 internal waves in a continuously stratified fluid
# ======================================================================================================================
@functools.lru_cache(maxsize=1)
def boussinesq_linear_sympy() -> dict:
    """The §7.8 reduction (curation D32–D33, ★★★) in sympy, for an arbitrary background N(z):
    (7.124) in (7.120)–(7.122) with the hydrostatic state (7.123) gives (7.128)–(7.130); (7.125) linearised gives (7.126);
    with (7.127) it is (7.131); and the w-equation (7.134) is the exact combination
    ∂²_t∇²w + N²∇²_H w = ∇²_H(∂_t R₃ − (g/ρ₀)R₄) + ∂_t∂_z(∂_t R₅ − ∂_x R₁ − ∂_y R₂)
    of the equation residuals R₁–R₃ (7.128)–(7.130), R₄ (7.131), R₅ (4.10) — so every solution of the linear set obeys
    (7.134), without constant N. Also: the plane wave (7.136) in (7.134) gives (7.137).
    Book: §7.8, Eqs. (7.120)–(7.137), (4.9)–(4.10). Returns dict(r7126, r7128, r7129, r7130, r7131, r7132 (∂_t R₅ −
    ∂_x R₁ − ∂_y R₂ minus (7.132)'s two sides), r7133 (∂_t R₃ − (g/ρ₀)R₄ minus (7.133)'s), r7134, r7137 — each 0;
    w_equation (the left side of (7.134)), residuals (the same dict), dispersion, N2). Cached. Validation: V2 —
    tests/test_ch07.py: test_boussinesq_linear_V2_derivation, test_w_equation_V2_derivation. Label: symbolic.
    """
    x, y, z, t = sp.symbols("x y z t", real=True)
    g, rho0, eps = sp.symbols("g rho0 epsilon", positive=True)
    u, v, w, pp, rp = [sp.Function(n)(x, y, z, t) for n in ("u", "v", "w", "p1", "r1")]
    rbar = sp.Function("rhobar")(z)
    pbar = sp.Function("pbar")(z)
    N2 = -g / rho0 * sp.diff(rbar, z)  # Eq. (7.127)
    hydro = {sp.diff(pbar, z): -rbar * g}  # Eq. (7.123)
    p_tot, r_tot = pbar + pp, rbar + rp  # Eq. (7.124)
    mom_z_full = sp.diff(w, t) + sp.diff(p_tot, z) / rho0 + r_tot * g / rho0  # (7.122) residual
    mom_z_pert = sp.diff(w, t) + sp.diff(pp, z) / rho0 + rp * g / rho0  # (7.130) residual
    res_130 = sp.simplify(mom_z_full.subs(hydro) - mom_z_pert)
    res_128 = sp.simplify((sp.diff(u, t) + sp.diff(p_tot, x) / rho0) - (sp.diff(u, t) + sp.diff(pp, x) / rho0))
    res_129 = sp.simplify((sp.diff(v, t) + sp.diff(p_tot, y) / rho0) - (sp.diff(v, t) + sp.diff(pp, y) / rho0))
    # (7.125) with small-amplitude scaling: u, v, w, ρ′ → ε·(…); keep O(ε)
    dens = sum(sp.diff(r_tot.subs(rp, eps * rp), s) * c for s, c in ((t, 1), (x, eps * u), (y, eps * v), (z, eps * w)))
    lin = sp.expand(dens).coeff(eps, 1)
    res_126 = sp.simplify(lin - (sp.diff(rp, t) + w * sp.diff(rbar, z)))  # Eq. (7.126)
    res_131 = sp.simplify((sp.diff(rp, t) + w * sp.diff(rbar, z)) - (sp.diff(rp, t) - N2 * rho0 / g * w))  # (7.131)
    R1 = sp.diff(u, t) + sp.diff(pp, x) / rho0
    R2 = sp.diff(v, t) + sp.diff(pp, y) / rho0
    R3 = mom_z_pert
    R4 = sp.diff(rp, t) - N2 * rho0 / g * w
    R5 = sp.diff(u, x) + sp.diff(v, y) + sp.diff(w, z)

    def lapH(f):
        return sp.diff(f, x, 2) + sp.diff(f, y, 2)

    weq = sp.diff(lapH(w) + sp.diff(w, z, 2), t, 2) + N2 * lapH(w)  # Eq. (7.134)
    Q1 = sp.diff(R5, t) - sp.diff(R1, x) - sp.diff(R2, y)
    Q2 = sp.diff(R3, t) - g / rho0 * R4
    res_132 = sp.simplify(sp.expand(Q1 - (sp.diff(w, z, t) - lapH(pp) / rho0)))  # Eq. (7.132)
    res_133 = sp.simplify(sp.expand(Q2 - (sp.diff(pp, t, z) / rho0 + sp.diff(w, t, 2) + N2 * w)))  # Eq. (7.133)
    combo = lapH(Q2) + sp.diff(Q1, t, z)
    res_134 = sp.simplify(sp.expand(weq - combo))
    # plane wave (7.136) with constant N
    k, l, m, om, Nc, w0 = sp.symbols("k l m omega N w0", positive=True)
    wp = w0 * sp.exp(sp.I * (k * x + l * y + m * z - om * t))
    weq_c = sp.diff(sp.diff(wp, x, 2) + sp.diff(wp, y, 2) + sp.diff(wp, z, 2), t, 2) + Nc ** 2 * (
        sp.diff(wp, x, 2) + sp.diff(wp, y, 2))
    disp = sp.solve(sp.Eq(sp.simplify(weq_c / wp), 0), om)
    res_137 = sp.simplify(sp.simplify(weq_c / wp).subs(om, Nc * sp.sqrt((k ** 2 + l ** 2) / (k ** 2 + l ** 2 + m ** 2))))
    res = {"r7126": res_126, "r7128": res_128, "r7129": res_129, "r7130": res_130, "r7131": res_131,
           "r7132": res_132, "r7133": res_133, "r7134": res_134, "r7137": res_137}
    return {**res, "w_equation": weq, "residuals": res, "dispersion": disp, "N2": N2}


def internal_wave_fields(x, z, t, k: float = 1.0, m: float = 1.0, N: float = 1.0, w0: float = 1.0,
                         rho0: float = 1000.0, g: float = G, residuals: bool = False, h: float | None = None,
                         ht: float | None = None) -> dict:
    """Plane internal gravity wave in the x–z plane with every field from the polarization relations (7.153):
    w = ŵ e^{i(kx+mz−ωt)} (7.136), u = −(m/k)ŵ e^{iθ}, p′ = −(ωmρ₀/k²)ŵ e^{iθ}, ρ′ = (iN²ρ₀/ωg)ŵ e^{iθ}, particle vertical
    displacement ζ = (i/ω)ŵ e^{iθ} (w = ∂ζ/∂t, so ρ′ = N²ρ₀ζ/g (7.149)); ω = N|k|/K (7.138, sign-safe); real parts.
    Book: §7.8, Eqs. (7.128)–(7.131), (7.136)–(7.141), (7.149), (7.153).
    Parameters: x, z [m]; t [s]; k ≠ 0, m [rad/m]; N [rad/s]; w0 = ŵ real amplitude [m/s]; rho0 [kg/m³]; g; residuals
    → add second-order finite-difference residuals of (4.10), (7.128), (7.130), (7.131) with steps h [m] (default
    1e-3/K) and ht [s] (default 1e-3/ω).
    Returns dict(u, w, p_prime [Pa], rho_prime [kg/m³], zeta_particle [m] (also ``zeta``), omega, K_dot_u (k u + m w,
    = 0: (7.141)), amplitudes (complex dict), [residuals]). Validation: V1, V2, V3 — tests/test_ch07.py:
    test_polarization_and_flux_V2_derivation, test_internal_wave_fields_V1_polarization_and_transversality,
    test_internal_wave_fields_V3_residuals_second_order. Checks: V2 polarization by sympy; V3 FD residuals order 2; V1
    K·u = 0. Label: analytic, converged.
    """
    k, m = float(k), float(m)
    om = float(internal_wave_omega(k, m, N))
    amps = {"w": complex(w0), "u": complex(-m / k * w0), "p": complex(-om * m * float(rho0) / k ** 2 * w0),  # (7.153)
            "rho": complex(1j * float(N) ** 2 * float(rho0) / (om * float(g)) * w0), "zeta": complex(1j / om * w0)}

    def fields(X, Z, T):
        ph = np.exp(1j * (k * _F(X) + m * _F(Z) - om * _F(T)))
        return {n: np.real(A * ph) for n, A in amps.items()}

    f = fields(x, z, t)
    out = {"u": _S(f["u"]), "w": _S(f["w"]), "p_prime": _S(f["p"]), "rho_prime": _S(f["rho"]),
           "zeta_particle": _S(f["zeta"]), "zeta": _S(f["zeta"]), "omega": om, "K_dot_u": _S(k * f["u"] + m * f["w"]),  # Eq. (7.141)
           "amplitudes": amps}
    if residuals:
        K = np.hypot(k, m)
        hh = 1e-3 / K if h is None else float(h)
        hT = 1e-3 / om if ht is None else float(ht)
        X, Z, T = _F(x), _F(z), _F(t)

        def d(name, dx=0.0, dz=0.0, dt=0.0, step=1.0):
            return (fields(X + dx, Z + dz, T + dt)[name] - fields(X - dx, Z - dz, T - dt)[name]) / (2 * step)

        cont = d("u", dx=hh, step=hh) + d("w", dz=hh, step=hh)  # (4.10)
        momx = d("u", dt=hT, step=hT) + d("p", dx=hh, step=hh) / float(rho0)  # Eq. (7.128)
        momz = d("w", dt=hT, step=hT) + d("p", dz=hh, step=hh) / float(rho0) + f["rho"] * float(g) / float(rho0)  # (7.130)
        dens = d("rho", dt=hT, step=hT) - float(N) ** 2 * float(rho0) / float(g) * f["w"]  # Eq. (7.131)
        out["residuals"] = {"continuity": _S(cont), "momentum_x": _S(momx), "momentum_z": _S(momz),
                            "density": _S(dens)}
    return out


def w_equation_residual(w_fn: Callable, x, z, t, N: float, h: float = 1e-2, ht: float = 1e-2):
    """Residual of the internal-wave equation ∂²/∂t²(∂²w/∂x² + ∂²w/∂z²) + N²∂²w/∂x² (7.134) in the x–z plane for a
    callable w(x, z, t), by nested second-order central differences (steps h [m], ht [s]; truncation O(h² + ht²),
    round-off ~ 1e-16/(h²ht²) — hence the default 1e-2 for unit-scale K and ω; scale the steps with 1/K and 1/ω).
    Book: §7.8, Eq. (7.134). Validation: V1 — tests/test_ch07.py:
    test_w_equation_residual_V1_plane_wave_right_and_wrong_frequency. Checks: V1 ≈ 0 for the plane wave with ω from
    (7.138), ≠ 0 with a wrong ω. Label: converged."""
    X, Z, T = _F(x), _F(z), _F(t)

    def d2(f, var, s):
        if var == "x":
            return lambda a, b, c: (f(a + s, b, c) - 2 * f(a, b, c) + f(a - s, b, c)) / s ** 2
        if var == "z":
            return lambda a, b, c: (f(a, b + s, c) - 2 * f(a, b, c) + f(a, b - s, c)) / s ** 2
        return lambda a, b, c: (f(a, b, c + s) - 2 * f(a, b, c) + f(a, b, c - s)) / s ** 2

    wf = lambda a, b, c: _F(w_fn(a, b, c))  # noqa: E731
    wxx, wzz = d2(wf, "x", h), d2(wf, "z", h)
    lap = lambda a, b, c: wxx(a, b, c) + wzz(a, b, c)  # noqa: E731
    return _S(d2(lap, "t", ht)(X, Z, T) + float(N) ** 2 * wxx(X, Z, T))  # Eq. (7.134)


def layered_flow_check(u_fn: Callable, v_fn: Callable, x, y, z=0.0, h: float = 1e-5) -> dict:
    """The steady solution w = p′ = ρ′ = 0 (§7.8, before Fig. 7.30): any horizontal, horizontally non-divergent (u, v)
    with arbitrary z-dependence, ∂u/∂x + ∂v/∂y = 0 (7.142), satisfies (4.10) and (7.128)–(7.131) — decoupled horizontal
    layers (pancake flow, blocking). Returns dict(divergence (7.142) by central differences, momentum_x = 0,
    momentum_z = 0, density = 0) for callables u(x, y, z), v(x, y, z). Book: §7.8, Eq. (7.142). Validation: V1 —
    tests/test_ch07.py: test_layered_flow_V1_horizontal_nondivergent_layers. Label: analytic."""
    X, Y, Z = _F(x), _F(y), _F(z)
    div = (_F(u_fn(X + h, Y, Z)) - _F(u_fn(X - h, Y, Z))) / (2 * h) \
        + (_F(v_fn(X, Y + h, Z)) - _F(v_fn(X, Y - h, Z))) / (2 * h)  # Eq. (7.142)
    zero = 0.0 * div
    return {"divergence": _S(div), "momentum_x": _S(zero), "momentum_z": _S(zero), "density": _S(zero)}


def internal_wave_energy(k: float = 1.0, m: float = 1.0, N: float = 1.0, w0: float = 1.0, rho0: float = 1000.0,
                         g: float = G) -> dict:
    """Mean energy and energy flux of a plane internal wave, per unit volume and per unit area.

    E_k = ¼ρ₀(m²/k² + 1)ŵ² (7.154); E_p = N²ρ₀ŵ²/(4ω²) (7.155); E_k = E_p (7.156); E = ½ρ₀(m²/k² + 1)ŵ² (7.157);
    F = ⟨p′u⟩ = (ρ₀ωmŵ²/2k²)(m/k e_x − e_z) (7.158); c_g E (unnumbered line) = F (7.159).
    Book: §7.8, Eqs. (7.147)–(7.159). Parameters: k ≠ 0, m [rad/m]; N [rad/s]; w0 [m/s]; rho0 [kg/m³]; g.
    Returns dict(Ek, Ep, E [J/m³], F (2,) [W/m²], cgE (2,), cg (2,) [m/s], omega). Sign-safe for k < 0.
    Validation: V1, V2 — tests/test_ch07.py: test_polarization_and_flux_V2_derivation,
    test_internal_wave_energy_V1_closed_forms_equal_averages. Checks: V1 E_k = E_p, F = c_g E (1e-14); closed forms =
    averages of the real fields. Label: analytic.
    """
    k, m, w0, r0 = float(k), float(m), float(w0), float(rho0)
    v = internal_wave_velocities(k, m, N)
    om = v["omega"]
    Ek = 0.25 * r0 * (m ** 2 / k ** 2 + 1.0) * w0 ** 2  # Eq. (7.154)
    Ep = float(N) ** 2 * r0 / (4.0 * om ** 2) * w0 ** 2  # Eq. (7.155)
    E = 0.5 * r0 * (m ** 2 / k ** 2 + 1.0) * w0 ** 2  # Eq. (7.157)
    F = r0 * om * m * w0 ** 2 / (2.0 * k ** 2) * np.array([m / k, -1.0])  # Eq. (7.158)
    return {"Ek": Ek, "Ep": Ep, "E": E, "F": F, "cgE": v["cg"] * E, "cg": v["cg"], "omega": om}  # Eq. (7.159)


def internal_energy_budget_residual(x, z, t, k: float = 1.0, m: float = 1.0, N: float = 1.0, w0: float = 1.0,
                                    rho0: float = 1000.0, g: float = G, h: float | None = None,
                                    ht: float | None = None, terms: bool = False):
    """Pointwise internal-wave energy equation (7.147) ∂/∂t[½ρ₀(u² + w²)] + gρ′w + ∇·(p′u) = 0 for the plane wave of
    :func:`internal_wave_fields`, by second-order central differences (h [m], default 1e-3/K; ht [s], default 1e-3/ω).
    Book: §7.8, Eqs. (7.147)–(7.148). Returns the residual [W/m³]; terms=True returns dict(dKE_dt, conversion = gρ′w,
    div_flux, residual). Validation: V4 — tests/test_ch07.py: test_internal_energy_budget_V4_pointwise_residual. Checks:
    V4 residual → 0 at order 2. Label: conserved, converged."""
    k_, m_ = float(k), float(m)
    om = float(internal_wave_omega(k_, m_, N))
    K = np.hypot(k_, m_)
    hh = 1e-3 / K if h is None else float(h)
    hT = 1e-3 / om if ht is None else float(ht)
    X, Z, T = _F(x), _F(z), _F(t)

    def f(a, b, c):
        return internal_wave_fields(a, b, c, k_, m_, N, w0, rho0, g)

    ke = lambda a, b, c: 0.5 * float(rho0) * (_F(f(a, b, c)["u"]) ** 2 + _F(f(a, b, c)["w"]) ** 2)  # noqa: E731
    dke = (ke(X, Z, T + hT) - ke(X, Z, T - hT)) / (2 * hT)
    f0 = f(X, Z, T)
    conv = float(g) * _F(f0["rho_prime"]) * _F(f0["w"])
    fx = lambda a, b, c: _F(f(a, b, c)["p_prime"]) * _F(f(a, b, c)["u"])  # noqa: E731
    fz = lambda a, b, c: _F(f(a, b, c)["p_prime"]) * _F(f(a, b, c)["w"])  # noqa: E731
    div = (fx(X + hh, Z, T) - fx(X - hh, Z, T)) / (2 * hh) + (fz(X, Z + hh, T) - fz(X, Z - hh, T)) / (2 * hh)
    res = _S(dke + conv + div)  # Eq. (7.147)
    if terms:
        return {"dKE_dt": _S(dke), "conversion": _S(conv), "div_flux": _S(div), "residual": res}
    return res


def internal_pe_interface_limit(a: float = 1.0, rho1: float = 1000.0, rho2: float = 1020.0, g: float = G,
                                rho0: float = 1000.0, eps: float = 1.0, k: float = 0.01, detail: bool = False):
    """Consistency of the continuous potential energy (7.150) E_p = ½N²ρ₀ζ² with the interface value (7.151)
    ¼(ρ₂ − ρ₁)ga²: smooth the density jump to ρ̄ = ½(ρ₁ + ρ₂) − ½(ρ₂ − ρ₁) tanh(z/ε), so N² = (g/ρ₀)(ρ₂ − ρ₁)
    sech²(z/ε)/(2ε) → (g/ρ₀)(ρ₂ − ρ₁)δ(z) (7.152); with the two-fluid particle displacement ζ = a e^{−k|z|} cos θ the
    column integral averaged over a wavelength, ∫½N²ρ₀⟨ζ²⟩dz, → ¼(ρ₂ − ρ₁)ga² as ε → 0 with error O(kε).
    Book: §7.8, Eqs. (7.150)–(7.152). Parameters: a [m]; rho1, rho2 [kg/m³]; g; rho0 reference density [kg/m³] (it
    cancels); eps profile half-width [m]; k wavenumber of the interfacial wave [rad/m] (the relative error is
    −2 ln 2·kε + O((kε)²) ≈ −1.39 kε, since ∫₀^∞ u sech²u du = ln 2 — the column value falls short; the defaults
    ε = 1 m, k = 0.01 rad/m give −1.37 %, measured). Returns ⟨∫E_p dz⟩ [J/m²]; detail=True → dict(Ep_column, target,
    rel_error). Our construction (the book argues with δ(z) directly). Validation: V3 — tests/test_ch07.py:
    test_internal_pe_interface_limit_V3_epsilon_order_one. Checks: V3 order ≥ 1 in ε. Label: converged."""
    r0 = float(rho0)
    drho = float(rho2) - float(rho1)
    e = float(eps)

    def integrand(zz):
        q = np.exp(-2.0 * abs(zz) / e)
        N2 = float(g) / r0 * drho / (2.0 * e) * 4.0 * q / (1.0 + q) ** 2  # smoothed (7.152): sech² = 4q/(1 + q)²
        return 0.5 * N2 * r0 * 0.5 * float(a) ** 2 * np.exp(-2.0 * float(k) * abs(zz))  # (7.150) averaged, ⟨cos²⟩ = ½

    zmax = 60.0 * e + 40.0 / float(k)
    pts = sorted({min(e, zmax / 2), min(10 * e, zmax / 2)})
    val = 2.0 * quad(integrand, 0.0, zmax, points=pts, epsabs=1e-14, epsrel=1e-12, limit=400)[0]
    target = 0.25 * drho * float(g) * float(a) ** 2  # Eq. (7.151)
    if detail:
        return {"Ep_column": val, "target": target, "rel_error": val / target - 1.0}
    return float(val)


def st_andrews_cross(x, z, t, omega: float = 0.71, N: float = 1.0, width: float = 1.0, amp: float = 1.0,
                     wavelength: float | None = None, detail: bool = False):
    """Illustrative St Andrew's cross: four beams of internal plane waves radiating from an oscillating source at the
    origin, each along a direction θ = arccos(ω/N) from the vertical with a Gaussian cross-beam envelope of half-width
    ``width``; in each beam K ⟂ beam axis (so ω(K) = N cos θ exactly) with the sign that makes c_g point away from the
    source — phase lines move across the beam, energy along it (c ⟂ c_g).
    Book: §7.8, text before and Fig. 7.33 (Mowbray & Rarity 1967) — **our illustration** (a superposition of beams, not
    a solution of the forced problem; labelled).
    Parameters: x (Nx,), z (Nz,) [m] (or 2-D arrays of the same shape); t [s]; omega < N [rad/s]; N; width [m]; amp;
    wavelength of the phase lines [m] (default 2·width). Returns the ρ′-like field (Nz, Nx) [–, scaled by amp];
    detail=True → dict(field, theta (from the vertical), beams: list of dict(e_beam, K, c, cg)).
    Validation: V7 — tests/test_ch07.py: test_st_andrews_cross_V7_beam_geometry. Checks: V7 beam axis at arccos(ω/N); c
    ⟂ c_g in each beam. Label: qualitative (illustration).
    """
    th = float(W.beam_angle(omega, N))
    if not np.isfinite(th):
        raise ValueError("omega must satisfy 0 <= omega <= N")
    lamb = 2.0 * float(width) if wavelength is None else float(wavelength)
    K0 = _TWO_PI / lamb
    X, Z = (np.meshgrid(_F(x), _F(z), indexing="xy") if np.ndim(x) == 1 else (_F(x), _F(z)))
    field = np.zeros_like(X, dtype=float)
    beams = []
    for sx in (1.0, -1.0):
        for sz in (1.0, -1.0):
            eb = np.array([sx * np.sin(th), sz * np.cos(th)])  # beam direction (c_g), θ from the vertical
            Kd = np.array([-eb[1], eb[0]])  # ⟂ beam
            v = internal_wave_velocities(K0 * Kd[0], K0 * Kd[1], N)
            if np.dot(v["cg"], eb) < 0:
                Kd = -Kd
                v = internal_wave_velocities(K0 * Kd[0], K0 * Kd[1], N)
            s_al = X * eb[0] + Z * eb[1]
            s_cr = -X * eb[1] + Z * eb[0]
            env = np.exp(-s_cr ** 2 / (2.0 * float(width) ** 2)) * 0.5 * (1.0 + np.tanh(s_al / float(width)))
            field += float(amp) * env * np.cos(K0 * (Kd[0] * X + Kd[1] * Z) - float(omega) * float(t))
            beams.append({"e_beam": eb, "K": K0 * Kd, "c": v["c"], "cg": v["cg"]})
    if detail:
        return {"field": field, "theta": th, "beams": beams}
    return field


def internal_wave_state(omega_over_N: float, N: float = 1.0, K: float = 1.0, k_sign: float = 1.0,
                        m_sign: float = 1.0) -> dict:
    """Explainer E9 state: an internal wave of frequency ω = (ω/N)·N with |K| = K; k = k_sign·K cos θ, m = m_sign·K sin θ
    (k_sign = −1, m_sign = +1 is the book's Fig. 7.29 geometry: K up-left, c_g down-left).
    Returns dict(theta_K (angle of K above the horizontal, rad), theta_K_deg, beam_from_vertical_deg (= θ_K: analysis
    R10), beam_from_horizontal_deg (= 90° − θ_K), k, m, omega, cx, cz, cgx, cgz, c = |c|, cg = |c_g| [m/s], dot (c·c_g,
    = 0), T = 2π/ω, T_N = 2π/N [s]). Book: §7.8, (7.138)–(7.146), Figs. 7.29, 7.31, 7.33. Scalar-callable.
    Validation: V1 — tests/test_ch07.py: test_internal_wave_velocities_V1_gradient_parity_both_signs. Label:
    analytic."""
    r = float(omega_over_N)
    th = float(np.arccos(np.clip(r, 0.0, 1.0)))
    k = float(np.sign(k_sign) or 1.0) * float(K) * np.cos(th)
    m = float(np.sign(m_sign) or 1.0) * float(K) * np.sin(th)
    v = internal_wave_velocities(k, m, N) if abs(k) > 0 else {"omega": 0.0, "c": np.zeros(2), "cg": np.zeros(2),
                                                            "dot": 0.0}
    om = r * float(N)
    return {"theta_K": th, "theta_K_deg": float(np.degrees(th)), "beam_from_vertical_deg": float(np.degrees(th)),
            "beam_from_horizontal_deg": 90.0 - float(np.degrees(th)), "k": float(k), "m": float(m), "omega": om,
            "cx": float(v["c"][0]), "cz": float(v["c"][1]), "cgx": float(v["cg"][0]), "cgz": float(v["cg"][1]),
            "c": float(np.linalg.norm(v["c"])), "cg": float(np.linalg.norm(v["cg"])), "dot": float(v["dot"]),
            "T": _TWO_PI / om if om > 0 else np.inf, "T_N": _TWO_PI / float(N)}


# ======================================================================================================================
# explainer state for particle orbits and dye lines (E2)
# ======================================================================================================================
def dyed_line(z0s, t, a: float = 0.1, k: float = 1.0, H=np.inf, g: float = G, x0: float = 0.0,
              rtol: float = 1e-9, atol: float = 1e-12) -> dict:
    """Positions at time(s) t of fluid particles that start on the vertical line x = x₀ at depths z₀ (t = 0), following
    the exact path lines (7.32)/(7.33) of the linear wave — the dyed line of Fig. 7.22 leaning forward (Stokes drift).
    The z₀ are **starting** positions, not orbit centres: a particle released at (x₀, z₀) circles about the mean depth
    z̄ = z₀ − ζ(x₀, z₀, 0) (ζ from (7.35b); at a crest, x₀ = 0, it starts a·sinh k(z₀ + H)/sinh kH above z̄), so its
    drift is ū_L(z̄), equal to ū_L(z₀) only to O(ka).
    Book: §7.6, Fig. 7.22 (our computation). Parameters: z0s (n,) [m]; t scalar or (M,) [s] (≥ 0); a, k, H, g; x0 [m];
    rtol, atol (``solve_ivp`` DOP853, all particles integrated together). Returns dict(t, x (M, n) — or (n,) for
    scalar t, z likewise).
    Validation: V1 — tests/test_ch07.py: test_dyed_line_V1_advances_by_the_stokes_drift. Checks: V1 identical to
    ``particle_path(model="exact", start="mean")``; net advance per period of each particle = ū_L(z̄)T within O(ka) (z̄
    the mean depth); the top drifts most. Label: analytic.
    """
    zs = np.atleast_1d(_F(z0s))
    n = zs.size
    k = float(k)
    om, uw, _ = _vel_fns(float(a), k, H, g)
    t_ = np.atleast_1d(_F(t))
    tt = np.unique(np.concatenate([[0.0], t_]))

    def rhs(tc, y):
        u, w = uw(y[:n], y[n:], tc)
        return np.concatenate([_F(u), _F(w)])

    y0 = np.concatenate([np.full(n, float(x0)), zs])
    if float(tt[-1]) == 0.0:  # only t = 0 asked for: the line itself
        Y, st = y0[:, None], np.array([0.0])
    else:
        sol = solve_ivp(rhs, (0.0, float(tt[-1])), y0, method="DOP853", t_eval=tt, rtol=rtol, atol=atol)
        Y, st = sol.y, sol.t
    idx = np.searchsorted(st, t_)
    X = Y[:n, idx].T
    Z = Y[n:, idx].T
    return {"t": t_, "x": X[0] if np.ndim(t) == 0 else X, "z": Z[0] if np.ndim(t) == 0 else Z}


def orbit_state(z0: float, a: float = 0.1, k: float = 1.0, H=np.inf, g: float = G, periods: int = 1,
                model: str = "linear") -> dict:
    """Explainer E2 state for one depth: dict(A, B, focal_half (7.36), B_over_A, omega, T, drift_speed (7.86),
    drift_per_period = ū_L T, eulerian_mean (0, N92), drift_numeric (distance drifted per period from exact path lines
    over ``periods`` periods; NaN unless model = "exact"), orbital_speed = ω·A).
    Book: §7.2 (7.35)–(7.36); §7.6 (7.85)–(7.86). Scalar-callable. Validation: V1 — tests/test_ch07.py:
    test_orbit_state_V1_explainer_numbers. Label: analytic."""
    ax = orbit_semi_axes(z0, a, k, H)
    om = float(omega_gravity(k, H, g))
    T = _TWO_PI / om
    uL = float(stokes_drift(z0, a, k, H, g))
    A_, B_ = float(ax["A"]), float(ax["B"])
    out = {"A": A_, "B": B_, "focal_half": ax["focal_half"], "B_over_A": B_ / A_ if A_ else np.nan, "omega": om,
           "T": T, "drift_speed": uL, "drift_per_period": uL * T, "eulerian_mean": 0.0, "drift_numeric": np.nan,
           "orbital_speed": om * A_}
    if model == "exact":
        out["drift_numeric"] = stokes_drift_numeric(z0, a, k, H, g, periods=periods) * T
    return out
