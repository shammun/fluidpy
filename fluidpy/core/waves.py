"""Linear wave primitives: the vocabulary of waves, dispersion relations of every Ch. 7 wave family (surface gravity,
capillary–gravity, interfacial, two-layer, internal), their inverses, phase and group velocity (scalar and vector, for
*any* dispersion relation), depth regimes, linear spectral evolution of 1-D and 2-D wave fields, beats and packets,
Hamiltonian ray tracing and wave energy.

Book: Kundu, Cohen & Dowling 5e, Ch. 7 §§7.1–7.8, Eqs. (7.1)–(7.9), (7.28)–(7.29), (7.42), (7.45), (7.49), (7.56)–(7.60),
(7.66)–(7.71), (7.73)–(7.79), (7.95)–(7.96), (7.110)–(7.117), (7.137)–(7.146). Every equation was read from the rendered
page images (chapters/pages/ch07/p282–p331). Ch. 11 (Kelvin–Helmholtz base state), Ch. 13 (Poincaré/Kelvin/Rossby waves,
reduced-gravity models, WKB rays) and Ch. 15 (acoustics) reuse this module.

Conventions (analysis §9 R1–R11, curation §8)
* x horizontal along propagation, **z up**, still surface z = 0, flat bottom z = −H; ``H = np.inf`` means deep water.
* ω is the (intrinsic, (7.9)) **angular frequency** [rad/s] — not the vorticity of Ch. 2–6. ω ≥ 0; the sign of k sets
  the direction (cos(kx − ωt) with k < 0 travels to −x).
* g defaults to ``G_BOOK`` = 9.81 m/s² (the book's value; ``core`` elsewhere defaults to G0 = 9.80665 — pass g explicitly
  when it matters).
* Reduced gravity: (7.117) g′ = g(ρ₂ − ρ₁)/ρ₂ (**lower** density in the denominator) — :func:`reduced_gravity_book`;
  ch04's ``core.similarity.reduced_gravity`` uses ρ₁ (``ref="upper"``). g′_lower/g′_upper = ρ₁/ρ₂;
  g′_lower/g′_mean = (ρ₁ + ρ₂)/(2ρ₂) (mean ρ̄ = (ρ₁ + ρ₂)/2).
* Sign of k: every ω(k) here depends on |k|; :func:`phase_speed` returns the speed |c| ≥ 0, while :func:`group_velocity`
  and :func:`group_velocity_numeric` return the **signed** dω/dk = sgn(k)·c_g (negative for a left-going wave).
* ⚠ Sibling argument orders differ (kept for API stability; call by keyword): ``omega_capillary_gravity(k, H, sigma,
  rho, g)`` with σ = 0.0727 N/m by default, but ``phase_speed``/``group_velocity``/``period_from_wavelength``/
  ``wavenumber_from_omega``/``wavelength_from_period(…, H, g, sigma, rho)`` with σ = 0 by default.
* Internal waves: the printed (7.138) ω = kN/K and (7.145) c_g = (Nm/K³)(m e_x − k e_z) assume k > 0; here ω = N|k|/K
  and c_g = ∇_K ω (sign-safe; the printed form is kept as ``printed=True`` for wrong-variant tests). θ is K's angle
  above the horizontal = the beam (c_g, particle motion) angle from the **vertical**.
* Vectors: a wavenumber vector ``K`` is a 1-D array (k, l, m) or (k, m). A point set ``X`` passed to
  :func:`plane_wave` has its components on the **last** axis, shape (…, d) (design contract); components-first arrays
  (d, …) are also accepted when the last axis is not d. Vector *outputs* of ch07 field functions (e.g.
  ``ch07.surface_normal``) carry their components on the first axis.
* Complex notation (§7.7 on): Re{} is dropped during the algebra and restored at the end (:func:`real_field`).

Numerics: overflow-safe hyperbolic ratios (cosh k(z + H)/sinh kH = e^{kz}(1 + e^{−2k(z+H)})/(1 − e^{−2kH}) — exact, no
overflow at kH = 500, and e^{kz} automatically when H = ∞); inverse dispersion by ``brentq`` with a physics bracket and an
asserted residual; group velocity by the complex-step derivative when ω is complex-safe, else fourth-order central
differences.
"""
from __future__ import annotations

import warnings
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from ._util import as_scalar_if_0d
from .thermo import G_BOOK

__all__ = [
    # vocabulary (§7.1)
    "sinusoid", "wave_parameters", "crest_positions", "plane_wave", "phase_velocity_vector", "trace_velocities",
    "doppler_frequency", "real_field",
    # hyperbolic depth profiles (§7.2)
    "depth_profiles", "cosh_over_sinh", "sinh_over_sinh", "cosh_over_cosh",
    # surface dispersion family (§7.2–7.3)
    "omega_gravity", "omega_capillary_gravity", "phase_speed", "group_velocity", "period_from_wavelength",
    "wavenumber_from_omega", "wavelength_from_period", "fenton_mckee_kh", "guo_kh", "depth_regime",
    "capillary_minimum", "min_group_velocity",
    # group velocity for any ω (§7.5)
    "group_velocity_numeric", "group_velocity_vector", "beat_wave", "gaussian_packet",
    # spectral evolution and rays (§7.5)
    "linear_evolve", "envelope", "linear_evolve_2d", "ray_trace",
    # energy
    "wave_energy_density", "viscous_decay",
    # interfaces and layers (§7.7)
    "eps2_density", "interface_omega", "two_layer_free_surface_omega", "two_layer_long_wave_speed",
    "reduced_gravity_book",
    # internal waves (§7.8)
    "internal_wave_omega", "beam_angle", "internal_wave_velocities",
]

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d
_TWO_PI = 2.0 * np.pi


def _isinf(H) -> np.ndarray:
    return np.isinf(np.asarray(H, dtype=float))


# ======================================================================================================================
# §7.1 vocabulary
# ======================================================================================================================
def sinusoid(x, t, a: float = 1.0, k: float = 1.0, omega: float = 1.0, direction: int = +1):
    """Sinusoidal travelling wave η = a cos(kx − ωt) (direction +1, (7.2)) or a cos(kx + ωt) (direction −1, (7.61)).

    Book: §7.1, Eqs. (7.1)–(7.2); §7.4, Eq. (7.61). (7.1) a cos[2π(x − ct)/λ] is the same wave with k = 2π/λ, c = ω/k.
    Parameters: x [m]; t [s]; a amplitude [m]; k wavenumber [rad/m]; omega angular frequency [rad/s]; direction ±1.
    Returns η [m] (broadcast of x and t). Scalar-callable.
    Validation: V1 — tests/test_ch07.py: test_sinusoid_V1_crests_ride_at_omega_over_k. Checks: V1 η(x_crest(t), t) = a
    on the crests of :func:`crest_positions`; V7 direction −1 moves the crests to −x. Label: analytic.
    """
    if direction not in (1, -1):
        raise ValueError("direction must be +1 or -1")
    return _S(float(a) * np.cos(float(k) * _F(x) - direction * float(omega) * _F(t)))  # Eq. (7.2) / (7.61)


def wave_parameters(*, k=None, lam=None, omega=None, T=None, nu=None, c=None) -> dict:
    """All the wave descriptors from any two independent ones: k = 2π/λ, T = 2π/ω = 1/ν, ω = 2πν, c = ω/k = λν.

    Book: §7.1, Eqs. (7.1)–(7.4) (and (7.7) λ = 2π/K when ``k`` is a wavenumber *vector*).
    Parameters (give exactly enough, extra values must be consistent): k [rad/m] (scalar, or a vector K whose magnitude
    is used); lam wavelength [m]; omega [rad/s]; T period [s]; nu cyclic frequency [Hz]; c phase speed [m/s].
    Returns dict(k, lam, omega, T, nu, c) as floats. Raises ValueError when under-determined or inconsistent (rel 1e-9).
    Validation: V1 — tests/test_ch07.py: test_wave_parameters_V1_round_trips_and_inconsistent_input. Checks: V1 round
    trips (k, T) → all → (λ, ν) → all; c = ω/k = λν. Label: analytic.
    """
    given = {"k": k, "lam": lam, "omega": omega, "T": T, "nu": nu, "c": c}
    if k is not None and np.ndim(k) > 0:
        k = float(np.linalg.norm(_F(k)))  # Eq. (7.6): K² = k² + l² + m²
    kk = float(k) if k is not None else (_TWO_PI / float(lam) if lam is not None else None)  # k = 2π/λ
    if omega is not None:
        om = float(omega)
    elif T is not None:
        om = _TWO_PI / float(T)
    elif nu is not None:
        om = _TWO_PI * float(nu)
    else:
        om = None
    if kk is None and om is None:
        raise ValueError("give one spatial (k or lam) or one temporal (omega, T, nu) quantity plus one more")
    if kk is None:
        if c is None:
            raise ValueError("under-determined: need k, lam or c together with the frequency")
        kk = om / float(c)
    if om is None:
        if c is None:
            raise ValueError("under-determined: need omega, T, nu or c together with the wavenumber")
        om = float(c) * kk
    out = {"k": kk, "lam": _TWO_PI / kk, "omega": om, "T": _TWO_PI / om if om != 0 else np.inf,
           "nu": om / _TWO_PI, "c": om / kk}  # Eq. (7.4): c = ω/k = λν
    for name, val in given.items():
        if val is None or name == "k":
            continue
        ref = out[name]
        if not np.isclose(float(val), ref, rtol=1e-9, atol=0.0):
            raise ValueError(f"inconsistent input: {name} = {val} but the other inputs give {ref}")
    return out


def crest_positions(t, k: float, omega: float, n=0):
    """Positions of the crests (phase = 2nπ) of η = a cos(kx − ωt): x_crest = (ω/k)t + 2nπ/k [m].

    Book: §7.1, Eq. (7.3) solved for x_crest (the line after (7.3)); its speed is (7.4). Parameters: t [s]; k [rad/m];
    omega [rad/s]; n integer crest index (array allowed). Validation: V1 — tests/test_ch07.py:
    test_sinusoid_V1_crests_ride_at_omega_over_k. Checks: V1 Δx/Δt = ω/k. Label: analytic.
    """
    return _S(float(omega) / float(k) * _F(t) + _TWO_PI * _F(n) / float(k))  # Eq. (7.3)


def plane_wave(X, K, omega: float, a: float = 1.0, t=0.0):
    """Three-dimensional (or 2-D) plane wave η = a cos(K·x − ωt) = a cos(kx + ly + mz − ωt).

    Book: §7.1, Eqs. (7.5)–(7.7). Parameters: X points, shape (…, d) with the components on the **last** axis (design
    contract; an array whose last axis is not d but whose first axis is, is read components-first) [m]; K wavenumber
    vector (d,) [rad/m]; omega [rad/s]; a [m]; t [s]. Returns η (…) [m]. Crests are the planes K·x − ωt = 2nπ, a
    distance λ = 2π/K apart along e_K (7.7). Validation: V1, V7 — tests/test_ch07.py:
    test_plane_wave_V1_crest_spacing_and_trace_velocities, test_plane_wave_V7_rotation_invariance. Checks: V1 crest
    spacing along e_K = 2π/K, along x = 2π/k. Label: analytic.
    """
    X_ = _F(X)
    K_ = _F(K)
    d = K_.size
    if X_.shape[-1] == d:
        KX = np.tensordot(X_, K_, axes=([-1], [0]))
    elif X_.shape[0] == d:
        KX = np.tensordot(K_, X_, axes=([0], [0]))
    else:
        raise ValueError(f"X must have {d} components on its last (or first) axis")
    return _S(float(a) * np.cos(KX - float(omega) * _F(t)))  # Eq. (7.5)


def phase_velocity_vector(K, omega: float) -> np.ndarray:
    """Phase-velocity vector c = (ω/K) e_K, e_K = K/K [m/s].

    Book: §7.1, Eq. (7.8) (re-displayed with (7.143) in §7.8). Parameters: K (d,) [rad/m]; omega [rad/s]. Returns (d,)
    array. Validation: V1, V7 — tests/test_ch07.py: test_plane_wave_V1_crest_spacing_and_trace_velocities,
    test_plane_wave_V7_rotation_invariance. Checks: V1 |c| = ω/K; rotation invariance. Label: analytic.
    """
    K_ = _F(K)
    Kmag = np.linalg.norm(K_)
    return float(omega) / Kmag * K_ / Kmag  # Eq. (7.8)


def trace_velocities(K, omega: float) -> tuple:
    """Trace velocities c_x = ω/k, c_y = ω/l, c_z = ω/m [m/s] (inf where a component is 0) — the speeds at which the
    crests cut each axis; each ≥ c = ω/K and **not** the components of the vector c.

    Book: §7.1, text after (7.8) and the inset of Fig. 7.1. Parameters: K (d,) [rad/m]; omega [rad/s]. Returns a
    **tuple** of d Python floats (c_x, c_y[, c_z]) [m/s] (design contract), not an array.
    Validation: V1 — tests/test_ch07.py: test_plane_wave_V1_crest_spacing_and_trace_velocities. Checks: V1 c_x ≥ c for
    random K; the reciprocals add as 1/c² = Σ 1/c_i². Label: analytic.
    """
    K_ = _F(K)
    with np.errstate(divide="ignore"):
        tr = np.where(K_ != 0, float(omega) / np.where(K_ != 0, K_, 1.0), np.inf)
    return tuple(float(v) for v in tr)  # a tuple of floats (design contract)


def doppler_frequency(omega, U, K):
    """Observed frequency at a fixed point in a uniform stream U: ω₀ = ω + U·K [rad/s] (Doppler shift).

    Book: §7.1, Eq. (7.9); the frozen pattern ω = 0 gives ω₀ = Uk. ω is the intrinsic frequency (the rest of the chapter
    uses intrinsic ω). Parameters: omega intrinsic [rad/s]; U (d,) [m/s]; K (d,) [rad/m]. Scalar U, K allowed (1-D).
    Validation: V1 — tests/test_ch07.py: test_doppler_V1_probe_frequency_of_a_translated_pattern. Checks: V1 frozen
    pattern → Uk; V1 FFT of a probe signal of a translated pattern. Label: analytic.
    """
    return _S(_F(omega) + float(np.dot(np.atleast_1d(_F(U)), np.atleast_1d(_F(K)))))  # Eq. (7.9)


def real_field(amp, phase):
    """Re{amp e^{i·phase}} — restoring the real field from the complex notation of §7.7–7.8.

    Book: §7.7, text before (7.89): ζ = Re{a exp[i(kx − ωt)]}, Re{} dropped during the algebra. Products (energies,
    fluxes) need real parts first: ⟨Re(Ae^{iθ}) Re(Be^{iθ})⟩ = ½ Re(AB*). Parameters: amp complex amplitude (any units);
    phase [rad]. Validation: V1 — tests/test_ch07.py: test_real_field_V1_real_parts_and_product_averages. Checks: V1
    Re{a e^{iθ}} = a cos θ. Label: analytic.
    """
    return _S(np.real(np.asarray(amp, dtype=complex) * np.exp(1j * _F(phase))))


# ======================================================================================================================
# hyperbolic depth structure (overflow-safe)
# ======================================================================================================================
def cosh_over_sinh(k, z, H=np.inf):
    """cosh k(z + H)/sinh kH [–], the depth structure of φ and u (7.26)–(7.27); e^{kz} when H = ∞.

    Computed as e^{kz}(1 + e^{−2k(z+H)})/(1 − e^{−2kH}) (algebraically identical, no overflow at kH = 500; −expm1 keeps
    small kH accurate). Book: §7.2, Eqs. (7.26), (7.27), (7.35); deep limit before (7.46). k > 0 [rad/m]; z [m] (≥ −H);
    H [m]. Validation: V1, V7 — tests/test_ch07.py: test_pressure_response_V1_surface_bottom_and_limits,
    test_wave_fields_V1_closed_forms_and_streamfunction, test_wave_fields_V7_deep_shallow_mirror_and_overflow. Label:
    analytic.
    """
    k_, z_, H_ = _F(k), _F(z), _F(H)
    with np.errstate(over="ignore", invalid="ignore"):
        return _S(np.exp(k_ * z_) * (1.0 + np.exp(-2.0 * k_ * (z_ + H_))) / (-np.expm1(-2.0 * k_ * H_)))


def sinh_over_sinh(k, z, H=np.inf):
    """sinh k(z + H)/sinh kH [–], the depth structure of w and ψ (7.27), (7.37); e^{kz} when H = ∞; 0 at the bottom.
    Overflow-safe form e^{kz}(−expm1(−2k(z + H)))/(−expm1(−2kH)). Book: §7.2, Eqs. (7.27), (7.35b), (7.37). Validation:
    V1, V7 — tests/test_ch07.py: test_pressure_response_V1_surface_bottom_and_limits,
    test_wave_fields_V1_closed_forms_and_streamfunction, test_wave_fields_V7_deep_shallow_mirror_and_overflow. Label:
    analytic.
    """
    k_, z_, H_ = _F(k), _F(z), _F(H)
    with np.errstate(over="ignore", invalid="ignore"):
        return _S(np.exp(k_ * z_) * (-np.expm1(-2.0 * k_ * (z_ + H_))) / (-np.expm1(-2.0 * k_ * H_)))


def cosh_over_cosh(k, z, H=np.inf):
    """cosh k(z + H)/cosh kH [–], the pressure response factor of (7.31); e^{kz} when H = ∞ (7.48), → 1 when kH → 0
    (7.52, hydrostatic). Overflow-safe form e^{kz}(1 + e^{−2k(z+H)})/(1 + e^{−2kH}). Book: §7.2, Eqs. (7.31), (7.48),
    (7.52). Validation: V1, V7 — tests/test_ch07.py: test_pressure_response_V1_surface_bottom_and_limits,
    test_wave_fields_V7_deep_shallow_mirror_and_overflow. Label: analytic."""
    k_, z_, H_ = _F(k), _F(z), _F(H)
    with np.errstate(over="ignore", invalid="ignore"):
        return _S(np.exp(k_ * z_) * (1.0 + np.exp(-2.0 * k_ * (z_ + H_))) / (1.0 + np.exp(-2.0 * k_ * H_)))


def depth_profiles(k, z, H=np.inf) -> dict:
    """The three hyperbolic depth ratios of §7.2 at once: dict(cosh_sinh, sinh_sinh, cosh_cosh) [–] (see
    :func:`cosh_over_sinh`, :func:`sinh_over_sinh`, :func:`cosh_over_cosh`). Book: §7.2, (7.26)–(7.31). Validation: V1 —
    tests/test_ch07.py: test_pressure_response_V1_surface_bottom_and_limits. Label: analytic."""
    return {"cosh_sinh": cosh_over_sinh(k, z, H), "sinh_sinh": sinh_over_sinh(k, z, H),
            "cosh_cosh": cosh_over_cosh(k, z, H)}


def _tanh_kH(k, H):
    """tanh(kH) that is 1 for H = ∞ and accepts complex k (complex-step derivatives).

    numpy's complex tanh evaluates sinh/cosh internally and raises a spurious "overflow" RuntimeWarning for
    Re(kH) ≳ 355 although the returned value (1 + O(1e-300)i) is correct — the warning is silenced (verification O9).
    """
    H_ = np.asarray(H, dtype=float)
    with np.errstate(over="ignore"):  # complex tanh at large kH: correct value, spurious overflow flag (O9)
        if H_.ndim == 0:
            if np.isinf(H_):
                return np.ones_like(k) if np.ndim(k) else (1.0 + 0.0 * k)
            return np.tanh(k * float(H_))
        inf = np.isinf(H_)
        return np.where(inf, 1.0 + 0.0 * k, np.tanh(k * np.where(inf, 1.0, H_)))


def _abs_k(k):
    """|k| for real input; for complex input (a complex step k + ih) the analytic continuation of |k|, i.e. −k where
    Re k < 0 and k elsewhere — so Im ω(|k + ih|)/h = sgn(k)·ω′(|k|) = dω/dk, the **signed** group velocity (a left-going
    wave has c_g < 0). Returning k unchanged for Re k < 0 (the pre-review code) evaluated √(gk) on the wrong branch and
    gave c_g ≈ 1e20 m/s (review M1)."""
    if np.iscomplexobj(k):
        return np.where(np.real(k) < 0, -k, k)
    return np.abs(_F(k))


def _x_over_sinh(x):
    """x/sinh x for x ≥ 0 (→ 1 at 0, → 0 at ∞), overflow-safe."""
    x_ = np.minimum(_F(x), 1e6)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        big = 2.0 * x_ * np.exp(-x_) / (-np.expm1(-2.0 * x_))
        small = np.where(x_ == 0, 1.0, x_ / np.sinh(np.where(x_ == 0, 1.0, x_)))
    return np.where(x_ > 20.0, big, small)


# ======================================================================================================================
# §7.2–7.3 surface-wave dispersion family
# ======================================================================================================================
def omega_gravity(k, H=np.inf, g: float = G_BOOK):
    """Dispersion relation of linear surface gravity waves ω = √(gk tanh kH) [rad/s].

    Book: §7.2, Eq. (7.28) (derived from the dynamic condition (7.21) applied to (7.26)); deep limit ω = √(gk) (7.45),
    shallow limit ω = k√(gH) (7.49); ω(k, x) with H = H(x) in (7.76).
    Parameters: k wavenumber [rad/m] (|k| is used; complex k accepted for complex-step derivatives); H depth [m] (np.inf
    = deep); g [m/s²] (default G_BOOK = 9.81). Returns ω [rad/s] ≥ 0.
    Assumptions: inviscid, irrotational, constant density, small amplitude (ka ≪ 1), no surface tension, air ignored.
    Validation: V1, V2, V6, V7 — tests/test_ch07.py: test_dispersion_V2_derivation,
    test_dispersion_V1_parity_with_ch04_and_identity, test_dispersion_V7_deep_and_shallow_limits,
    test_wave_functions_V7_change_of_units, test_book_V6_section_7_2_forms_and_numbers. Checks: V2 sympy φ (7.26) with
    this ω satisfies (7.11), (7.12), (7.18), (7.21); V7 deep/shallow limits; V1 parity with ch04.linear_wave_surface's
    ω. Label: analytic, symbolic.
    """
    kk = _abs_k(k)
    return _S(np.sqrt(float(g) * kk * _tanh_kH(kk, H)))  # Eq. (7.28)


def omega_capillary_gravity(k, H=np.inf, sigma: float = 0.0727, rho: float = 1000.0, g: float = G_BOOK):
    """Capillary–gravity dispersion relation ω = √(k(g + σk²/ρ) tanh kH) [rad/s].

    Book: §7.3, Eq. (7.56) (the dynamic condition (7.55) with the Laplace jump (7.54) replaces (7.21); φ keeps the form
    (7.26)). σ = 0 recovers (7.28); g = 0, H = ∞ gives pure capillary waves (7.60).
    Parameters: k [rad/m] (|k|; complex allowed); H [m]; sigma surface tension [N/m] (default 0.0727, clean water near
    20 °C); rho liquid density [kg/m³]; g [m/s²]. Returns ω [rad/s].
    Assumptions: as :func:`omega_gravity` plus a clean interface of constant σ; the air's density neglected.
    ⚠ Argument order (sigma, rho, g) and default σ = 0.0727 differ from :func:`phase_speed`/:func:`group_velocity`
    (g, sigma, rho; σ = 0) — pass σ, ρ, g by keyword when mixing them (review should-fix 6; signatures kept stable).
    Validation: V1, V2, V5, V6, V7 — tests/test_ch07.py: test_capillary_V2_derivation_tension_condition,
    test_capillary_V7_limits, test_capillary_V5_air_water_minimum, test_wavenumber_from_omega_V1_inverse_is_identity,
    test_book_V6_capillary_standing_and_group_numbers. Checks: V2 sympy (7.55) + (7.26) ⇒ (7.56); V7 σ → 0 and g → 0
    limits. Label: analytic, symbolic.
    """
    kk = _abs_k(k)
    return _S(np.sqrt(kk * (float(g) + float(sigma) * kk ** 2 / float(rho)) * _tanh_kH(kk, H)))  # Eq. (7.56)


def phase_speed(k, H=np.inf, g: float = G_BOOK, sigma: float = 0.0, rho: float = 1000.0):
    """Phase speed c = ω/k = √((g/k + σk/ρ) tanh kH) [m/s] of linear gravity (σ = 0) or capillary–gravity waves.

    Book: §7.2, Eq. (7.29); §7.3, Eq. (7.57); limits (7.45) √(g/k) (deep), (7.49) √(gH) (shallow), (7.60) √(2πσ/ρλ)
    (pure capillary, g = 0). Longer gravity waves are faster ⇒ dispersive.
    Parameters: k [rad/m] (|k|); H [m] (np.inf deep); g [m/s²]; sigma [N/m] (default 0); rho [kg/m³].
    Returns c [m/s] — the speed |ω/k| ≥ 0 for either sign of k (unlike :func:`group_velocity`, which is signed).
    Assumptions: as :func:`omega_capillary_gravity`. ⚠ Argument order (H, g, sigma, rho) and σ = 0 default differ from
    :func:`omega_capillary_gravity` (H, sigma, rho, g; σ = 0.0727) — call by keyword when mixing them.
    Validation: V1, V2, V7 — tests/test_ch07.py: test_phase_speed_V2_derivation_longer_is_faster,
    test_phase_speed_V7_limits_monotonicity_and_bound, test_dispersion_V7_deep_and_shallow_limits,
    test_kdv_linear_phase_speed_V1_taylor_of_7_29, test_wave_functions_V7_change_of_units. Checks: V7 kH → ∞ ⇒ √(g/k),
    kH → 0 ⇒ √(gH)(1 − (kH)²/6); V1 c = ω/k. Label: analytic.
    """
    kk = _abs_k(k)
    return _S(np.sqrt((float(g) / kk + float(sigma) * kk / float(rho)) * _tanh_kH(kk, H)))  # Eq. (7.29)/(7.57)


def group_velocity(k, H=np.inf, g: float = G_BOOK, sigma: float = 0.0, rho: float = 1000.0):
    """Group velocity c_g = dω/dk of (7.56) [m/s]; for σ = 0 it is the book's c_g = (c/2)[1 + 2kH/sinh 2kH].

    Book: §7.5, Eqs. (7.67), (7.69), limits (7.70) c_g = c/2 (deep), c (shallow); pure capillary c_g = 3c/2 (Exercise 7.9,
    text after (7.70)); local c_g = ∂ω(k, x)/∂k (7.77). The capillary–gravity form (our differentiation of (7.56)):
    c_g = (c/2)[(g + 3σk²/ρ)/(g + σk²/ρ) + 2kH/sinh 2kH].
    ⚠ Argument order (H, g, sigma, rho; σ = 0) differs from :func:`omega_capillary_gravity` (H, sigma, rho, g;
    σ = 0.0727) — call by keyword when mixing them.
    **Sign convention** (review M1): the returned value is the signed derivative dω/dk of ω(|k|), i.e.
    sgn(k)·(c/2)[… + 2|k|H/sinh 2|k|H] — positive for k > 0 (identical to (7.69)), negative for a left-going wave
    (k < 0), the same number :func:`group_velocity_numeric` gives. k = 0 is treated as k → 0⁺.
    Parameters: k [rad/m] (either sign); H [m]; g [m/s²]; sigma [N/m]; rho [kg/m³]. Returns c_g [m/s] (signed, same
    sign as k).
    Validation: V1, V2, V7 — tests/test_ch07.py: test_group_velocity_V2_derivation (sympy d/dk √(gk tanh kH) − (7.69)
    = 0, capillary form), test_group_velocity_V1_complex_step_parity, test_group_velocity_V1_negative_k_is_signed_derivative,
    test_group_velocity_V7_limits (½, 1, 3/2; c_g = c at k_m), test_capillary_minimum_V1_numerical_minimisation,
    test_wave_functions_V7_change_of_units. Label: analytic, symbolic.
    """
    k_ = _F(k)
    sgn = np.where(k_ < 0, -1.0, 1.0)  # dω/dk = sgn(k) ω′(|k|)
    kk = np.abs(k_)
    s = float(sigma) / float(rho)
    c = _F(phase_speed(kk, H, g, sigma, rho))
    ratio = (float(g) + 3.0 * s * kk ** 2) / (float(g) + s * kk ** 2)
    kH = np.where(_isinf(H), np.inf, kk * np.where(_isinf(H), 1.0, _F(H)))
    return _S(sgn * 0.5 * c * (ratio + _x_over_sinh(2.0 * kH)))  # Eq. (7.69) (σ = 0: ratio = 1), signed


def period_from_wavelength(lam, H=np.inf, g: float = G_BOOK, sigma: float = 0.0, rho: float = 1000.0):
    """Wave period from the wavelength, T = √((2πλ/g) coth(2πH/λ)) [s] (σ = 0), or 2π/ω(2π/λ) with surface tension.

    Book: §7.2, Eq. (7.28), second form. Parameters: lam [m]; H [m]; g [m/s²]; sigma [N/m]; rho [kg/m³].
    Validation: V1, V2 — tests/test_ch07.py: test_dispersion_V2_derivation,
    test_wavenumber_from_omega_V1_inverse_is_identity. Checks: V1 equals 2π/ω(2π/λ); deep T = √(2πλ/g). Label: analytic.
    """
    lam_ = _F(lam)
    if float(sigma) == 0.0:
        th = _F(_tanh_kH(_TWO_PI / lam_, H))
        return _S(np.sqrt(_TWO_PI * lam_ / float(g) / th))  # Eq. (7.28): T = √((2πλ/g) coth(2πH/λ))
    return _S(_TWO_PI / _F(omega_capillary_gravity(_TWO_PI / lam_, H, sigma, rho, g)))


def _k_from_omega_scalar(om: float, H: float, g: float, sigma: float, rho: float, rtol: float) -> float:
    if om == 0.0:
        return 0.0
    if om < 0 or not np.isfinite(om):
        raise ValueError(f"omega must be finite and >= 0; got {om}")
    s = sigma / rho
    kd = om * om / g if g > 0 else np.inf
    if s == 0.0 and np.isinf(H):
        return kd  # deep water, closed form k = ω²/g
    f = (lambda kk: (g * kk + s * kk ** 3) * (np.tanh(kk * H) if np.isfinite(H) else 1.0) - om * om)
    if s == 0.0:
        ks = om / np.sqrt(g * H)
        lo = max(kd, ks)
        hi = kd / np.tanh(kd * H)  # k ≥ kd ⇒ tanh kH ≥ tanh(kd H) ⇒ k ≤ kd/tanh(kd H)
        lo, hi = lo * (1 - 1e-12), hi * (1 + 1e-12)
        if f(lo) > 0:
            lo = lo * 0.5
    else:
        lo = 1e-300
        hi = max(kd if np.isfinite(kd) else 0.0, (om * om / s) ** (1.0 / 3.0),
                 om / np.sqrt(g * H) if (np.isfinite(H) and g > 0) else 0.0, 1e-12)
        while f(hi) <= 0:
            hi *= 2.0
    k = brentq(f, lo, hi, xtol=1e-300, rtol=rtol, maxiter=500)
    res = abs(f(k)) / (om * om)
    if res > 1e-12:
        raise RuntimeError(f"inverse dispersion did not converge (relative residual {res:.2e})")
    return k


def wavenumber_from_omega(omega, H=np.inf, g: float = G_BOOK, sigma: float = 0.0, rho: float = 1000.0,
                          rtol: float = 1e-15):
    """Invert the dispersion relation: the wavenumber k [rad/m] with ω(k) = omega for (7.28) (σ = 0) or (7.56).

    Book: §7.2, Eq. (7.28) (inverse); §7.3, (7.56). Deep water and σ = 0: closed form k = ω²/g. Otherwise ``brentq``
    on (gk + σk³/ρ) tanh kH − ω² with the physics bracket [max(ω²/g, ω/√(gH)), (ω²/g)/tanh(ω²H/g)] (both limits bound the
    true root; ω(k) is monotonic), residual asserted < 1e-12 relative.
    Parameters: omega [rad/s] (≥ 0, array allowed); H [m]; g [m/s²]; sigma [N/m]; rho [kg/m³]; rtol for brentq.
    Returns k [rad/m]. Validation: V1, V6 — tests/test_ch07.py: test_wavenumber_from_omega_V1_inverse_is_identity,
    test_ray_trace_V1_snell_closed_form_parity, test_book_V6_section_7_2_forms_and_numbers. Checks: V1 inverse ∘ forward
    = identity over ω ∈ [1e-4, 1e3], H ∈ [1e-3, 1e4] (the evidence for this function); the explicit approximations of
    Fenton & McKee (1990, max λ error 1.7 %) and Guo (2002, 0.79 % for the coded exponent 5/2) are themselves checked
    *against* it (their V5 tests), so they do not raise its label. Label: analytic.
    """
    om = _F(omega)
    H_ = float(H)
    out = np.vectorize(lambda w: _k_from_omega_scalar(float(w), H_, float(g), float(sigma), float(rho), rtol),
                       otypes=[float])(om)
    return _S(out)


def wavelength_from_period(T, H=np.inf, g: float = G_BOOK, sigma: float = 0.0, rho: float = 1000.0):
    """Wavelength λ = 2π/k [m] of a wave of period T [s] on depth H [m] (inverse of (7.28)); deep water λ = gT²/2π.

    Book: §7.2, Eq. (7.28) and the ocean example after (7.45) (T ≈ 10 s: our 156.1 m in deep water).
    Validation: V1, V6 — tests/test_ch07.py: test_wavenumber_from_omega_V1_inverse_is_identity,
    test_book_V6_section_7_2_forms_and_numbers. Checks: V1 T = 10 s deep ⇒ 156.13 m. Label: analytic.
    """
    return _S(_TWO_PI / _F(wavenumber_from_omega(_TWO_PI / _F(T), H, g, sigma, rho)))


def fenton_mckee_kh(omega, H, g: float = G_BOOK):
    """Explicit approximation kH = (ω²H/g)[coth((ω√(H/g))^{3/2})]^{2/3} [–] (Fenton & McKee 1990), maximum wavelength
    error 1.7 % (the primary Fenton & McKee (1990) bound; ours measures 1.66 % in λ), exact in both limits — an
    approximation checked against :func:`wavenumber_from_omega`, not used by the physics.

    Source: J. D. Fenton, "A note on two approximations to the linear dispersion relation for surface gravity water
    waves" (2006), Eq. (2), https://johndfenton.com/Papers/Dispersion-Relation.pdf (read 2026-09-24).
    Parameters: omega [rad/s]; H [m]; g [m/s²]. Validation: V5 — tests/test_ch07.py:
    test_dispersion_V5_fenton_mckee_1990, test_dispersion_V5_guo_2002. Label: benchmark.
    """
    x = _F(omega) * np.sqrt(_F(H) / float(g))
    return _S(x ** 2 * (1.0 / np.tanh(x ** 1.5)) ** (2.0 / 3.0))


def guo_kh(omega, H, g: float = G_BOOK):
    """Explicit approximation kH = (ω²H/g)(1 − exp(−(ω√(H/g))^{5/2}))^{−2/5} [–] (Guo 2002, in Fenton's 5/2 form),
    maximum error 0.79 % for this exponent (Guo's own fitted β = 2.4908 gives 0.75 %), exact in both limits.
    Source: Fenton (2006) Eq. (3), https://johndfenton.com/Papers/Dispersion-Relation.pdf (read 2026-09-24).
    Parameters: omega [rad/s]; H [m]; g [m/s²]. Validation: V5 — tests/test_ch07.py: test_dispersion_V5_guo_2002. Label:
    benchmark.
    """
    x = _F(omega) * np.sqrt(_F(H) / float(g))
    return _S(x ** 2 * (-np.expm1(-x ** 2.5)) ** (-0.4))


def depth_regime(k, H, deep_kH: float = 2.0, shallow_H_over_lambda: float = 0.07) -> dict:
    """Deep / intermediate / shallow classification with the size of each limit's error.

    Book: §7.2, text after (7.45) (deep if kH > 2, i.e. H > 0.32λ, "within 2 %") and after (7.49) (shallow if
    H < 0.07λ, "better than 3 %"). Errors returned (never only the label, analysis R11):
    deep_error = 1 − √(tanh kH) (relative error of c = √(g/k), (7.45)); shallow_error = 1 − √(tanh kH/kH) (of c = √(gH),
    (7.49)). Parameters: k [rad/m]; H [m]; thresholds (book values by default).
    Returns dict(regime, kH, H_over_lambda, deep_error, shallow_error) (floats/str for scalar input).
    Validation: V1, V6 — tests/test_ch07.py: test_depth_regime_V1_book_thresholds_and_errors,
    test_book_V6_section_7_2_forms_and_numbers. Checks: V1 kH = 2 ⇒ H/λ = 0.3183, deep error 1.82 %; H = 0.07λ ⇒ shallow
    error 3.04 %. Label: analytic.
    """
    kk = np.abs(_F(k))
    H_ = _F(H)
    kH = kk * H_
    th = np.where(np.isinf(kH), 1.0, np.tanh(np.where(np.isinf(kH), 1.0, kH)))
    deep_err = 1.0 - np.sqrt(th)
    with np.errstate(divide="ignore", invalid="ignore"):
        shallow_err = np.where(np.isinf(kH), 1.0, 1.0 - np.sqrt(th / np.where(kH == 0, 1.0, kH)))
    shallow_err = np.where(kH == 0, 0.0, shallow_err)
    HL = kH / _TWO_PI
    regime = np.where(kH > deep_kH, "deep", np.where(HL < shallow_H_over_lambda, "shallow", "intermediate"))
    reg = regime.item() if regime.ndim == 0 else regime
    return {"regime": reg, "kH": _S(kH), "H_over_lambda": _S(HL), "deep_error": _S(deep_err),
            "shallow_error": _S(shallow_err)}


def capillary_minimum(sigma: float = 0.0727, rho: float = 1000.0, g: float = G_BOOK) -> dict:
    """Minimum phase speed of deep-water capillary–gravity waves: c_min = (4gσ/ρ)^{1/4} at λ_m = 2π√(σ/ρg).

    Book: §7.3, Eq. (7.58) (dc/dλ = 0 in (7.57) with tanh → 1); (7.59) for air–water. At k_m = √(ρg/σ) gravity and surface
    tension contribute equally (g/k = σk/ρ) and c_g = c.
    Parameters: sigma [N/m]; rho [kg/m³]; g [m/s²]. Returns dict(c_min [m/s], lam_m [m], k_m [rad/m]).
    Validation: V1, V2, V5, V7 — tests/test_ch07.py: test_capillary_minimum_V2_derivation,
    test_capillary_minimum_V1_numerical_minimisation, test_capillary_V5_air_water_minimum,
    test_wave_functions_V7_change_of_units. Checks: V1 equals a numerical minimisation of (7.57); σ = 0.073, ρ = 1000, g
    = 9.81 ⇒ 23.13 cm/s at 1.714 cm; V5 Wikipedia "Capillary wave" (0.23 m/s, 1.7 cm). Label: analytic, benchmark.
    """
    c_min = (4.0 * float(g) * float(sigma) / float(rho)) ** 0.25  # Eq. (7.58)
    lam_m = _TWO_PI * np.sqrt(float(sigma) / (float(rho) * float(g)))  # Eq. (7.58)
    return {"c_min": float(c_min), "lam_m": float(lam_m), "k_m": float(_TWO_PI / lam_m)}


def min_group_velocity(sigma: float = 0.0727, rho: float = 1000.0, g: float = G_BOOK) -> dict:
    """Minimum group velocity of deep-water capillary–gravity waves (the calm centre of the stone-in-a-pond train).

    Book: §7.5, text before Fig. 7.16 ("c_g,min", Exercise 7.10). Our derivation: with x = σk²/(ρg),
    c_g = (c/2)(1 + 3x)/(1 + x) and c_g² ∝ x^{−1/2}(1 + 3x)²/(1 + x); d/dx = 0 ⇒ 3x² + 6x − 1 = 0 ⇒ x = 2/√3 − 1.
    Parameters: sigma [N/m]; rho [kg/m³]; g [m/s²]. Returns dict(cg_min [m/s], k [rad/m], lam [m], c [m/s], x).
    Validation: V1, V6 — tests/test_ch07.py: test_capillary_minimum_V1_numerical_minimisation,
    test_book_V6_capillary_standing_and_group_numbers. Checks: V1 equals ``minimize_scalar`` of :func:`group_velocity`;
    V6 the book's value. Label: analytic.
    """
    x = 2.0 / np.sqrt(3.0) - 1.0
    k = np.sqrt(x * float(rho) * float(g) / float(sigma))
    cg = float(group_velocity(k, np.inf, g, sigma, rho))
    return {"cg_min": cg, "k": float(k), "lam": float(_TWO_PI / k), "c": float(phase_speed(k, np.inf, g, sigma, rho)),
            "x": float(x)}


# ======================================================================================================================
# §7.5 group velocity for any dispersion relation, beats and packets
# ======================================================================================================================
def _central4(f: Callable, x, h):
    return (-f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h) + f(x - 2 * h)) / (12.0 * h)


def _disp_fn(H=np.inf, g: float = G_BOOK, sigma: float = 0.0, rho: float = 1000.0) -> Callable:
    """The surface-wave dispersion relation (7.28)/(7.56) as a complex-safe callable ω(k) (dispersion given as data);
    ω(k) = ω(|k|) with |k| analytically continued for a complex step (:func:`_abs_k`), so dω/dk is signed for k < 0."""
    return lambda kk: omega_capillary_gravity(kk, H, sigma, rho, g)


def group_velocity_numeric(omega_fn: Callable | None = None, k=1.0, h=None, method: str = "auto", **disp):
    """c_g = dω/dk [m/s] of any dispersion relation ω(k), by the complex-step derivative Im ω(k + ih)/h (exact to round-off
    when ω is written with complex-safe numpy) or by fourth-order central differences.

    Book: §7.5, Eq. (7.67) c_g = Δω/Δk → dω/dk (tangent of ω(k), Fig. 7.14).
    Parameters: omega_fn callable ω(k) [rad/s] (None → the surface relation (7.56) built from ``**disp`` = H, g, sigma,
    rho — the dispersion as data, for explainer parity rows); k [rad/m] (array allowed); h step (default
    1e-20·max(|k|, 1) complex, 1e-3·max(|k|, 1e-12) central); method "auto" (complex step if ω returns complex output,
    else central) | "complex" | "central". Returns c_g [m/s] — the signed derivative: for the built-in surface relation
    and k < 0 it is −c_g(|k|) (a left-going wave), equal to :func:`group_velocity`. A user ``omega_fn`` used with the
    complex step must be analytic near k (write |k| as ``np.where(np.real(k) < 0, -k, k)``, not ``np.abs``).
    Validation: V1, V3 — tests/test_ch07.py: test_group_velocity_V1_complex_step_parity,
    test_group_velocity_V1_negative_k_is_signed_derivative,
    test_group_velocity_numeric_V3_central_differences_fourth_order. Checks: V1 vs analytic c_g (gravity, capillary,
    internal, Rossby-type −βk/(k² + l²)), both signs of k, deep and finite depth; V3 central differences converge at
    order 4. Label: analytic, converged.
    """
    if omega_fn is None:
        omega_fn = _disp_fn(**disp)
    k_ = _F(k)
    if method in ("auto", "complex"):
        hc = 1e-20 * np.maximum(np.abs(k_), 1.0) if h is None else _F(h)
        try:
            w = omega_fn(k_ + 1j * hc)
            if np.iscomplexobj(w):
                return _S(np.imag(w) / hc)
            if method == "complex":
                raise TypeError("omega_fn returned a real value for complex input (not complex-safe)")
        except (TypeError, ValueError):
            if method == "complex":
                raise
    if method not in ("auto", "central", "complex"):
        raise ValueError("method must be 'auto', 'complex' or 'central'")
    hr = 1e-3 * np.maximum(np.abs(k_), 1e-12) if h is None else _F(h)
    return _S(_central4(lambda kk: _F(omega_fn(kk)), k_, hr))


def group_velocity_vector(omega_fn: Callable, K, h=None, method: str = "auto") -> np.ndarray:
    """Group-velocity vector c_gi = ∂ω/∂K_i [m/s] — the gradient of ω in wavenumber space.

    Book: §7.5, unnumbered c_gi = ∂ω/∂K_i (before Fig. 7.16); §7.8, Eq. (7.143).
    Parameters: omega_fn callable ω(K) taking a 1-D array K (components) [rad/s]; K (d,) [rad/m]; h, method as in
    :func:`group_velocity_numeric`. Returns (d,) array.
    Validation: V1 — tests/test_ch07.py: test_internal_wave_velocities_V1_gradient_parity_both_signs. Checks: V1
    internal waves: equals (7.145) for k > 0 and the sign-safe form for k < 0. Label: analytic.
    """
    K_ = _F(K)
    out = np.empty_like(K_)
    for i in range(K_.size):
        def w_i(ki, i=i):
            Kc = K_.astype(np.asarray(ki).dtype, copy=True)
            Kc[i] = ki
            return omega_fn(Kc)
        out[i] = float(group_velocity_numeric(w_i, K_[i], h=h, method=method))
    return out


def beat_wave(x, t, k1: float, k2: float, a: float = 1.0, H=np.inf, g: float = G_BOOK, sigma: float = 0.0,
              rho: float = 1000.0, omega_fn: Callable | None = None, printed: bool = False) -> dict:
    """Two equal-amplitude waves of nearby wavenumber: η = 2a cos(½Δk x − ½Δω t) cos(kx − ωt) (beats).

    Book: §7.5, Eq. (7.66) and (7.67). ⚠️ The book prints the envelope argument as ½Δk x − ½Δω **x**; the next line
    and Fig. 7.13 use ½Δω **t** (sum-to-product identity) — implemented with t; ``printed=True`` returns the printed
    (wrong) form for the wrong-variant test (its envelope does not move).
    Parameters: x [m]; t [s]; k1, k2 [rad/m]; a [m]; H, g, sigma, rho the surface dispersion (7.56) (used when
    omega_fn is None); omega_fn optional callable ω(k); printed.
    Returns dict(eta, envelope (the slowly varying 2a cos(…), signed), carrier, k, omega, dk, domega, c, cg_finite
    (= Δω/Δk)). Validation: V1, V2, V6 — tests/test_ch07.py: test_beats_V2_derivation,
    test_beat_wave_V1_equals_sum_printed_fails_and_chord_to_tangent, test_book_V6_capillary_standing_and_group_numbers.
    Checks: V1 equals the plain sum a cos(k₁x − ω₁t) + a cos(k₂x − ω₂t) (the printed form must fail). Label: analytic.
    """
    wf = omega_fn if omega_fn is not None else _disp_fn(H, g, sigma, rho)
    w1, w2 = float(wf(k1)), float(wf(k2))
    dk, dw = float(k2) - float(k1), w2 - w1  # Δk = k₂ − k₁, Δω = ω₂ − ω₁
    k, w = 0.5 * (k1 + k2), 0.5 * (w1 + w2)
    x_, t_ = _F(x), _F(t)
    env_arg = 0.5 * dk * x_ - 0.5 * dw * (x_ if printed else t_)
    envl = 2.0 * float(a) * np.cos(env_arg)
    carrier = np.cos(k * x_ - w * t_)
    return {"eta": _S(envl * carrier),  # Eq. (7.66) (with ½Δω t)
            "envelope": _S(envl + 0.0 * carrier), "carrier": _S(carrier + 0.0 * envl), "k": k, "omega": w,
            "dk": dk, "domega": dw, "c": w / k, "cg_finite": dw / dk}  # Eq. (7.67): c_g = Δω/Δk


def gaussian_packet(x, t, a: float = 1.0, k0: float = 1.0, sigma_x: float = 10.0, omega_fn: Callable | None = None,
                    order: int = 2, H=np.inf, g: float = G_BOOK, sigma: float = 0.0, rho: float = 1000.0) -> dict:
    """Closed-form Gaussian wave packet η(x, 0) = a e^{−x²/2σ_x²} cos k₀x evolved with ω(k) Taylor-expanded about k₀.

    order = 2 (default): ω ≈ ω₀ + c_g κ + ½ω″κ² — the exact evolution for quadratic dispersion (a chirped, spreading
    Gaussian): η = Re{A e^{i(k₀x − ω₀t)}}, A = a σ_x/√(σ_x² + iω″t) · exp(−(x − c_g t)²/(2(σ_x² + iω″t))).
    order = 1: ω ≈ ω₀ + c_g κ gives the book's (7.68) η = a(x − c_g t) cos(k₀x − ω₀t) — the envelope rides at c_g.
    Book: §7.5, Eq. (7.68) (stated there with a pointer to Phillips 1977; written out as curation D20); Fig. 7.15.
    Our closed form (Gaussian integral), used for explainer parity without an FFT.
    Parameters: x [m]; t [s]; a [m]; k0 carrier wavenumber [rad/m] (k₀ < 0: a packet travelling to −x, the mirror
    image x → −x of the k₀ > 0 packet); sigma_x envelope half-width [m] (spectral width 1/σ_x); omega_fn ω(k)
    (default: the surface relation (7.56) from H, g, sigma, rho); order 1 or 2.
    Returns dict(eta, envelope |A|, carrier, c (= ω₀/k₀, signed), cg (signed dω/dk at k₀), omega0 (≥ 0), omega2
    (= d²ω/dk², even in k₀), x_center (= c_g t)).
    Validation: V1, V2 — tests/test_ch07.py: test_packet_envelope_V2_derivation,
    test_gaussian_packet_V1_envelope_rides_at_cg, test_gaussian_packet_V1_negative_k0_is_the_mirror_image. Checks: V1
    order 1 translates rigidly at c_g; order 2 equals :func:`linear_evolve` for quadratic dispersion; k₀ < 0 is the
    mirror image (deep and finite depth). Label: analytic.
    """
    wf = omega_fn if omega_fn is not None else _disp_fn(H, g, sigma, rho)
    k0f = float(k0)
    w0 = float(wf(k0f))
    cg = float(group_velocity_numeric(wf, k0f))
    hk = 1e-3 * abs(k0f)
    w2 = float((group_velocity_numeric(wf, k0f + hk) - group_velocity_numeric(wf, k0f - hk)) / (2 * hk)) \
        if order >= 2 else 0.0
    x_, t_ = _F(x), _F(t)
    s2 = float(sigma_x) ** 2 + 1j * w2 * t_
    A = float(a) * float(sigma_x) / np.sqrt(s2) * np.exp(-(x_ - cg * t_) ** 2 / (2.0 * s2))  # Eq. (7.68) for order 1
    carrier = np.exp(1j * (k0f * x_ - w0 * t_))
    return {"eta": _S(np.real(A * carrier)), "envelope": _S(np.abs(A)), "carrier": _S(np.real(carrier)),
            "c": w0 / k0f, "cg": cg, "omega0": w0, "omega2": w2, "x_center": _S(cg * t_)}


# ======================================================================================================================
# spectral evolution and rays
# ======================================================================================================================
def _omega_on_grid(omega_fn: Callable, kabs: np.ndarray) -> np.ndarray:
    w = np.zeros_like(kabs)
    pos = kabs > 0
    w[pos] = _F(omega_fn(kabs[pos]))
    return w


def linear_evolve(eta0, x, t, omega_fn: Callable | None = None, eta_t0=None, direction: int | None = None,
                  **disp) -> np.ndarray:
    """Evolve a linear dispersive wave field on a periodic grid: each Fourier mode moves with its own ω(k).

    Real ``eta0``: η̂(k, t) = η̂₀ cos ωt + (η̂_t0/ω) sin ωt (released from the given η and ∂η/∂t; ∂η/∂t = 0 by default,
    so a hump splits into two halves going both ways — the Cauchy–Poisson stone-in-a-pond problem); or, with
    ``direction = ±1``, one-way propagation η̂(k, t) = η̂₀ e^{∓i sgn(k) ω t} (a packet that travels only to +x or −x).
    Complex ``eta0`` (an analytic packet such as A(x)e^{ik₀x}): every mode is multiplied by e^{−iω(|k|)t} and the real
    part is returned. The k = 0 mode keeps η̂_t0·t.
    Book: §7.1 Fourier superposition (Exercise 7.3; curation N02) — our method; used for packets (7.68), the stone in a
    pond (Fig. 7.16) and crest tracking (Fig. 7.17).
    Parameters: eta0 (N,) initial elevation [m] (real, or complex analytic) on a uniform periodic grid x (N,) [m];
    t scalar or (M,) [s]; omega_fn ω(|k|) ≥ 0 [rad/s] (called with k > 0 only; None → the surface relation (7.56) from
    ``**disp`` = H, g, sigma, rho); eta_t0 (N,) [m/s] or None; direction None | +1 | −1.
    Returns η (N,) or (M, N) [m]. The domain must be long enough that the fastest component does not wrap round.
    Validation: V1, V4, V7 — tests/test_ch07.py: test_linear_evolve_V1_single_modes_and_superposition,
    test_linear_evolve_V7_nondispersive_translation, test_energy_flux_V4_packet_energy_conserved_and_carried_at_cg,
    test_gaussian_packet_V1_envelope_rides_at_cg. Checks: V1 one Fourier mode reproduces a cos(kx − ωt); V4 Σ|η̂|²-type
    energy conserved; V7 ω = ck translates any shape; V3 packet envelope speed → c_g(k₀). Label: analytic, conserved.
    """
    if omega_fn is None:
        omega_fn = _disp_fn(**disp)
    x_ = _F(x)
    N = x_.size
    dx = x_[1] - x_[0]
    kx = _TWO_PI * np.fft.fftfreq(N, d=dx)
    w = _omega_on_grid(omega_fn, np.abs(kx))
    t_ = np.atleast_1d(_F(t))
    tt = t_[:, None]
    if np.iscomplexobj(eta0):
        E0 = np.fft.fft(np.asarray(eta0, dtype=complex))
        out = np.real(np.fft.ifft(E0[None, :] * np.exp(-1j * w[None, :] * tt), axis=-1))
        return out[0] if np.ndim(t) == 0 else out
    E0 = np.fft.fft(_F(eta0))
    if direction is not None:
        if direction not in (1, -1):
            raise ValueError("direction must be None, +1 or -1")
        Et = E0[None, :] * np.exp(-1j * direction * np.sign(kx)[None, :] * w[None, :] * tt)
    else:
        Et0 = np.zeros(N, complex) if eta_t0 is None else np.fft.fft(_F(eta_t0))
        with np.errstate(divide="ignore", invalid="ignore"):
            sin_over_w = np.where(w[None, :] > 0, np.sin(w[None, :] * tt) / np.where(w > 0, w, 1.0)[None, :], tt)
        Et = E0[None, :] * np.cos(w[None, :] * tt) + Et0[None, :] * sin_over_w
    out = np.real(np.fft.ifft(Et, axis=-1))
    return out[0] if np.ndim(t) == 0 else out


def envelope(eta, axis: int = -1) -> np.ndarray:
    """Amplitude envelope |η + i H[η]| of a narrow-band signal (Hilbert transform, ``scipy.signal.hilbert``) [same units].
    Book: §7.5, the envelope a(x) of Fig. 7.15 / (7.68) — our diagnostic. Periodic signals; accurate when the band is
    narrow. Validation: V1, V4 — tests/test_ch07.py: test_envelope_V1_hilbert_modulus,
    test_energy_flux_V4_packet_energy_conserved_and_carried_at_cg. Checks: V1 envelope of a cos(kx)·e^{−x²/2s²} ≈ a
    e^{−x²/2s²}. Label: analytic."""
    from scipy.signal import hilbert

    return np.abs(hilbert(_F(eta), axis=axis))


def linear_evolve_2d(field0, x, z, t, omega_fn_K: Callable, direction_K=None, field_t0=None) -> np.ndarray:
    """Two-dimensional version of :func:`linear_evolve` on a doubly periodic (z, x) grid (``field[j, i]`` = f(z_j, x_i)).

    Each 2-D Fourier mode (k, m) moves with ω(k, m) = ``omega_fn_K(k, m)`` (e.g. internal waves N|k|/√(k² + m²), (7.138)).
    Complex ``field0`` (an analytic packet A(x, z)e^{i(k₀x + m₀z)}, the design's usage): every mode × e^{−iω(k, m)t},
    real part returned. Real ``field0``: ``direction_K = (k₀, m₀)`` gives one-way evolution (modes with K·K₀ > 0 get
    e^{−iωt}, the others the conjugate); ``None``: release from rest (cos ωt, plus field_t0 sin ωt/ω).
    Book: §7.8, Figs. 7.29, 7.32 (a packet sliding along its crests, c ⟂ c_g) — our method.
    Parameters: field0 (Nz, Nx) real or complex; x (Nx,) [m]; z (Nz,) [m]; t scalar or (M,) [s]; omega_fn_K callable
    (k, m) arrays → ω [rad/s] (called with finite arrays; K = 0 is set to ω = 0). Returns (Nz, Nx) or (M, Nz, Nx).
    Validation: V1 — tests/test_ch07.py: test_linear_evolve_2d_V1_plane_wave_and_packet_at_cg. Checks: V1 one mode
    reproduces a plane wave; V7 packet centroid moves at c_g (7.145). Label: analytic.
    """
    x_, z_ = _F(x), _F(z)
    kx = _TWO_PI * np.fft.fftfreq(x_.size, d=x_[1] - x_[0])
    kz = _TWO_PI * np.fft.fftfreq(z_.size, d=z_[1] - z_[0])
    KX, KZ = np.meshgrid(kx, kz, indexing="xy")
    Kmag = np.hypot(KX, KZ)
    w = np.zeros_like(KX)
    nz = Kmag > 0
    w[nz] = _F(omega_fn_K(KX[nz], KZ[nz]))
    t_ = np.atleast_1d(_F(t))[:, None, None]
    if np.iscomplexobj(field0):
        F0 = np.fft.fft2(np.asarray(field0, dtype=complex))
        out = np.real(np.fft.ifft2(F0[None] * np.exp(-1j * w[None] * t_), axes=(-2, -1)))
        return out[0] if np.ndim(t) == 0 else out
    F0 = np.fft.fft2(_F(field0))
    if direction_K is not None:
        s = np.sign(KX * float(direction_K[0]) + KZ * float(direction_K[1]))
        Ft = F0[None] * np.exp(-1j * s[None] * w[None] * t_)
    else:
        Ft1 = np.zeros_like(F0) if field_t0 is None else np.fft.fft2(_F(field_t0))
        with np.errstate(divide="ignore", invalid="ignore"):
            sw = np.where(w[None] > 0, np.sin(w[None] * t_) / np.where(w > 0, w, 1.0)[None], t_)
        Ft = F0[None] * np.cos(w[None] * t_) + Ft1[None] * sw
    out = np.real(np.fft.ifft2(Ft, axes=(-2, -1)))
    return out[0] if np.ndim(t) == 0 else out


def _grad_central(f: Callable, v: np.ndarray, h: np.ndarray) -> np.ndarray:
    g_ = np.empty_like(v)
    for i in range(v.size):
        e = np.zeros_like(v)
        e[i] = h[i]
        g_[i] = (-f(v + 2 * e) + 8 * f(v + e) - 8 * f(v - e) + f(v - 2 * e)) / (12.0 * h[i])
    return g_


def ray_trace(omega_fn: Callable | None = None, x0=(0.0, 0.0), k0=(0.01, 0.0), t_span=(0.0, 100.0),
              H_fn: Callable | None = None, dim: int | None = None, g: float = G_BOOK, n_out: int = 200,
              rtol: float = 1e-10, atol: float = 1e-12, t_eval=None, H_min: float = 1e-3,
              max_step: float = np.inf) -> dict:
    """Trace a ray of a slowly varying wave train: dx/dt = ∂ω/∂k (= c_g), dk/dt = −∂ω/∂x (Hamilton's equations).

    Along the ray the frequency is constant, ∂ω/∂t + c_g·∇ω = 0 (7.79), while k, c and c_g change — rays bend
    (refraction, Figs. 7.8, 7.9, 7.18). Book: §7.5, Eqs. (7.73)–(7.79); the dk/dt equation is our extension (the book
    derives only the conservation of ω; it follows from crest conservation (7.74) in the same way).
    Parameters: omega_fn ω(k_vec, x_vec) [rad/s] (default: surface gravity ω = √(g|k| tanh(|k|H(x))) with ``H_fn``);
    x0 start position (scalar or (dim,)) [m]; k0 start wavenumber (scalar or (dim,)) [rad/m]; t_span [s]; H_fn depth
    H(x_vec) [m]; dim 1 or 2 (inferred from x0); g; n_out output times evenly spaced over t_span (unless t_eval is
    given); rtol, atol (``solve_ivp`` RK45, rtol 1e-10 — our choice, tighter than the design's 1e-9: the relative ω
    drift stays ≈ 2e-9 up to the shore — 1.9e-9 measured on a 1:50 beach to H = 2 cm); H_min [m]: a terminal event
    stops the ray when H(x) < H_min (the shore); max_step [s].
    Derivatives: fourth-order central differences with steps 1e-4|k| (in k) and 1e-3/|k| (in x, a thousandth of a
    wavelength), accurate to ~1e-12 relative when the medium varies on scales ≫ λ.
    Returns dict(t (n,), x (dim, n), k (dim, n), omega (n,), omega_drift (max |ω − ω₀|/ω₀), status, message).
    Validation: V1, V2, V4 — tests/test_ch07.py: test_frequency_along_rays_V2_derivation,
    test_ray_trace_V1_homogeneous_straight_at_cg, test_ray_trace_V4_frequency_conserved_while_k_grows,
    test_ray_trace_V1_snell_closed_form_parity. Checks: V1 homogeneous: straight rays at c_g; V4 ω conserved over a
    slope (< 1e-8 asserted); V1 Snell |k| sin α constant for straight contours. Label: analytic, conserved.
    """
    x0_ = np.atleast_1d(_F(x0))
    k0_ = np.atleast_1d(_F(k0))
    d = dim if dim is not None else x0_.size
    if x0_.size != d or k0_.size != d:
        raise ValueError("x0 and k0 must have dim components")
    if omega_fn is None:
        if H_fn is None:
            raise ValueError("give omega_fn(k, x) or H_fn(x)")

        def omega_fn(kv, xv):  # noqa: F811
            kk = float(np.linalg.norm(kv))
            return float(omega_gravity(kk, float(H_fn(xv if d > 1 else xv[0])), g))  # Eq. (7.76)

    def w(kv, xv):
        return float(omega_fn(kv, xv))

    def rhs(t, y):
        xv, kv = y[:d], y[d:]
        kmag = max(float(np.linalg.norm(kv)), 1e-300)
        hk = np.full(d, 1e-4 * kmag)
        hx = np.full(d, 1e-3 / kmag)
        dwdk = _grad_central(lambda kk: w(kk, xv), kv, hk)  # ẋ = ∂ω/∂k (c_g, (7.77))
        dwdx = _grad_central(lambda xx: w(kv, xx), xv, hx)  # k̇ = −∂ω/∂x
        return np.concatenate([dwdk, -dwdx])

    events = None
    if H_fn is not None:
        def shore(t, y):
            return float(H_fn(y[:d] if d > 1 else y[0])) - H_min
        shore.terminal = True
        shore.direction = -1
        events = shore
    te = np.linspace(float(t_span[0]), float(t_span[1]), int(n_out)) if t_eval is None else _F(t_eval)
    sol = solve_ivp(rhs, t_span, np.concatenate([x0_, k0_]), method="RK45", t_eval=te, rtol=rtol, atol=atol,
                    events=events, max_step=max_step)
    X = sol.y[:d]
    Kv = sol.y[d:]
    om = np.array([w(Kv[:, i], X[:, i]) for i in range(len(sol.t))])
    w0 = w(k0_, x0_)
    return {"t": sol.t, "x": X, "k": Kv, "omega": om,  # x, k shaped (dim, n) (design contract)
            "omega_drift": float(np.max(np.abs(om - w0)) / abs(w0)) if om.size else 0.0,  # Eq. (7.79) check
            "status": sol.status, "message": sol.message}


# ======================================================================================================================
# energy
# ======================================================================================================================
def wave_energy_density(a, rho: float = 1000.0, g: float = G_BOOK, drho: float | None = None):
    """Wave energy per unit horizontal area E = ½ρga² [J/m²] (surface, (7.42)) or ½(ρ₂ − ρ₁)ga² (interface, (7.96)).

    Book: §7.2, Eq. (7.42) E = E_k + E_p = ρg⟨η²⟩ = ½ρga² (sinusoidal η); §7.7, Eq. (7.96). Parameters: a amplitude [m];
    rho [kg/m³]; g [m/s²]; drho = ρ₂ − ρ₁ [kg/m³] for an interface (then rho is ignored).
    Validation: V1 — tests/test_ch07.py: test_wave_energy_V1_quadrature_equals_closed_form,
    test_energy_flux_V1_quadrature_equals_closed_form, test_interface_energy_V1_quarter_from_direct_integration. Checks:
    V1 equals ``ch07.wave_energy`` quad; F = E c_g (7.71). Label: analytic.
    """
    r = float(rho) if drho is None else float(drho)
    return _S(0.5 * r * float(g) * _F(a) ** 2)  # Eq. (7.42) / (7.96)


def viscous_decay(a0, k, nu: float, t):
    """Viscous amplitude decay of a deep-water wave a(t) = a₀ e^{−2νk²t} [m] (Exercise 7.11; the end of the stone-in-a-pond
    train, §7.5). Book: §7.5 (pointer to Exercise 7.11; the result is the classical Lamb/Stokes decay rate, derived in
    the curation's demoted a-D68 from ε = 2νS_ijS_ij of (7.47)). Parameters: a0 [m]; k [rad/m]; nu [m²/s]; t [s].
    Validation: V2, V7 — tests/test_ch07.py: test_viscous_decay_V2_derivation, test_wave_functions_V7_change_of_units.
    Checks: V2 sympy energy-dissipation derivation. Label: analytic."""
    return _S(_F(a0) * np.exp(-2.0 * float(nu) * _F(k) ** 2 * _F(t)))


# ======================================================================================================================
# §7.7 interfaces and layers
# ======================================================================================================================
def eps2_density(rho1, rho2):
    """ε² = (ρ₂ − ρ₁)/(ρ₂ + ρ₁) [–] of (7.95) (ρ₁ upper/lighter, ρ₂ lower/heavier). Named ``eps2_density`` to avoid the
    dissipation ε of Ch. 4. Book: §7.7, Eq. (7.95). Validation: V7 — tests/test_ch07.py:
    test_interface_omega_V7_limits_and_rayleigh_taylor. Label: analytic."""
    return _S((_F(rho2) - _F(rho1)) / (_F(rho2) + _F(rho1)))


def interface_omega(k, rho1: float, rho2: float, g: float = G_BOOK):
    """Waves on the interface between two infinitely deep fluids: ω = √(gk(ρ₂ − ρ₁)/(ρ₂ + ρ₁)) = ε√(gk) [rad/s].

    Book: §7.7, Eq. (7.95) (from (7.93)–(7.94) with φ₁ = Ae^{−kz}e^{i(kx−ωt)}, φ₂ = Be^{kz}e^{i(kx−ωt)}, A = −B = iωa/k).
    Deep-water waves slowed by ε; ρ₁ → 0 recovers (7.45).
    Parameters: k [rad/m] (|k|); rho1 upper (lighter) [kg/m³]; rho2 lower (heavier) [kg/m³]; g [m/s²].
    Returns ω [rad/s]; **NaN with a RuntimeWarning when ρ₁ > ρ₂** (heavy over light: ω² < 0, Rayleigh–Taylor instability,
    Ch. 11) — never a silent real number. Assumptions: no interfacial tension, linear, irrotational in each layer.
    Validation: V1, V2, V7 — tests/test_ch07.py: test_interface_V2_derivation,
    test_interface_fields_V1_residuals_vortex_sheet_and_parity, test_interface_omega_V7_limits_and_rayleigh_taylor,
    test_interface_V1_form_published_two_fluid_dispersion. Checks: V2 sympy (7.93)–(7.94) ⇒ (7.95); V7 ρ₁ → 0 ⇒ √(gk),
    ρ₁ → ρ₂ ⇒ 0. Label: analytic, symbolic.
    """
    e2 = _F(eps2_density(rho1, rho2))
    if np.any(e2 < 0):
        warnings.warn("rho1 > rho2: heavy fluid over light — omega^2 < 0 (Rayleigh–Taylor instability, Ch. 11); "
                      "returning NaN", RuntimeWarning, stacklevel=2)
    kk = np.abs(_F(k))
    with np.errstate(invalid="ignore"):
        return _S(np.where(e2 >= 0, np.sqrt(np.maximum(e2, 0.0) * float(g) * kk), np.nan) + 0.0 * kk)  # Eq. (7.95)


def two_layer_free_surface_omega(k, H: float, rho1: float, rho2: float, g: float = G_BOOK) -> tuple:
    """The two roots of the layer-over-deep-fluid dispersion relation (7.110): barotropic ω² = gk (7.111) and baroclinic
    ω² = gk(ρ₂ − ρ₁) sinh kH/(ρ₂ cosh kH + ρ₁ sinh kH) (7.113) [rad/s].

    Book: §7.7, Eqs. (7.110), (7.111), (7.113). Upper layer ρ₁ of thickness H with a free surface over an infinitely deep
    layer ρ₂; origin at the mean free surface. (7.113) evaluated as gkΔρ tanh kH/(ρ₂ + ρ₁ tanh kH) (overflow-safe); → (7.95)
    as kH → ∞; → (7.115) as kH → 0.
    Parameters: k [rad/m] (|k|); H upper-layer thickness [m]; rho1 < rho2 [kg/m³]; g [m/s²].
    Returns the tuple (omega_bt, omega_bc) [rad/s] (baroclinic NaN + warning if ρ₁ > ρ₂).
    Validation: V2, V6, V7 — tests/test_ch07.py: test_two_layer_modes_V2_derivation,
    test_two_layer_V7_limits_and_reduced_gravity, test_book_V6_nonlinear_interface_and_internal. Checks: V2 both roots
    satisfy (7.110); V7 limits. Label: analytic, symbolic.
    """
    kk = np.abs(_F(k))
    th = np.tanh(kk * float(H))
    drho = float(rho2) - float(rho1)
    if drho < 0:
        warnings.warn("rho1 > rho2: the baroclinic mode is unstable (Ch. 11); returning NaN", RuntimeWarning,
                      stacklevel=2)
    w_bt = np.sqrt(float(g) * kk)  # Eq. (7.111): ω² = gk
    with np.errstate(invalid="ignore"):
        w2_bc = float(g) * kk * drho * th / (float(rho2) + float(rho1) * th)  # Eq. (7.113)
        w_bc = np.where(w2_bc >= 0, np.sqrt(np.maximum(w2_bc, 0.0)), np.nan) if drho >= 0 else np.full_like(kk, np.nan)
    return _S(w_bt), _S(w_bc)


def reduced_gravity_book(rho1: float, rho2: float, g: float = G_BOOK, ref: str = "lower"):
    """Reduced gravity g′ [m/s²]: ``ref="lower"`` → the book's (7.117) g′ = g(ρ₂ − ρ₁)/ρ₂; ``ref="upper"`` → ch04's
    g(ρ₂ − ρ₁)/ρ₁ (``core.similarity.reduced_gravity``); ``ref="mean"`` → g(ρ₂ − ρ₁)/ρ̄, ρ̄ = (ρ₁ + ρ₂)/2 (Boussinesq ρ₀).
    Ratios: lower/upper = ρ₁/ρ₂, lower/mean = (ρ₁ + ρ₂)/(2ρ₂); both differ from 1 by O(Δρ/ρ) (≲ 0.3 % in the ocean,
    the mean form by half as much). Book: §7.7, Eq. (7.117); ch04 §4.11 (after (4.105)).
    Parameters: rho1 upper, rho2 lower [kg/m³]; g; ref. Validation: V6, V7 — tests/test_ch07.py:
    test_two_layer_V7_limits_and_reduced_gravity, test_book_V6_nonlinear_interface_and_internal. Checks: V1 ratio
    lower/upper = ρ₁/ρ₂. Label: analytic."""
    ref_rho = {"lower": float(rho2), "upper": float(rho1), "mean": 0.5 * (float(rho1) + float(rho2))}
    if ref not in ref_rho:
        raise ValueError("ref must be 'lower' (book 7.117), 'upper' (ch04) or 'mean'")
    return float(g) * (float(rho2) - float(rho1)) / ref_rho[ref]  # Eq. (7.117) for ref="lower"


def two_layer_long_wave_speed(H: float, rho1: float, rho2: float, g: float = G_BOOK):
    """Long-wave (kH ≪ 1) speed of the baroclinic mode c = √(g′H), g′ = g(ρ₂ − ρ₁)/ρ₂ [m/s].

    Book: §7.7, Eqs. (7.115)–(7.117). Like √(gH) but reduced by √((ρ₂ − ρ₁)/ρ₂): internal waves are slow.
    Parameters: H upper-layer thickness [m]; rho1, rho2 [kg/m³]; g. Validation: V7 — tests/test_ch07.py:
    test_two_layer_V7_limits_and_reduced_gravity, test_wave_functions_V7_change_of_units. Checks: V7 = lim_{k→0} ω_bc/k
    of :func:`two_layer_free_surface_omega`. Label: analytic.
    """
    return _S(np.sqrt(reduced_gravity_book(rho1, rho2, g, "lower") * _F(H)))  # Eq. (7.116)


# ======================================================================================================================
# §7.8 internal waves (constant N)
# ======================================================================================================================
def internal_wave_omega(k, m, N: float, l=0.0):
    """Internal-gravity-wave dispersion relation ω = N√((k² + l²)/(k² + l² + m²)) = N|k_H|/K = N cos θ [rad/s].

    Book: §7.8, Eqs. (7.137), (7.138) (the printed kN/K assumes k > 0 — here the sign-safe |k|), (7.139) ω = N cos θ,
    θ = K's angle above the horizontal: frequency set by direction only, 0 ≤ ω ≤ N.
    Parameters: k, l horizontal wavenumbers [rad/m]; m vertical wavenumber [rad/m]; N buoyancy frequency [rad/s].
    Complex-safe (k² inside the square root) so complex-step c_g works for either sign of k.
    Assumptions: Boussinesq, inviscid, linear, constant N, ω ≫ Coriolis frequency.
    Validation: V1, V2, V7 — tests/test_ch07.py: test_internal_dispersion_V2_derivation,
    test_internal_wave_omega_V7_limits_and_direction_only,
    test_w_equation_residual_V1_plane_wave_right_and_wrong_frequency. Checks: V2 sympy (7.136) in (7.134) ⇒ (7.137); V7
    m = 0 ⇒ ω = N (parcel oscillation, ch01). Label: analytic, symbolic.
    """
    kh2 = k * k + l * l if (np.iscomplexobj(k) or np.iscomplexobj(m)) else _F(k) ** 2 + _F(l) ** 2
    K2 = kh2 + (m * m if np.iscomplexobj(m) else _F(m) ** 2)
    return _S(float(N) * np.sqrt(kh2 / K2))  # Eq. (7.137)


def beam_angle(omega, N: float):
    """Beam angle θ = arccos(ω/N) [rad]: the angle of the energy beams (c_g, particle motion) from the **vertical** — equal
    to the angle of K (and c) above the horizontal in (7.139).

    Book: §7.8, text before Fig. 7.33 ("four beams oriented at an angle θ with the vertical, cos θ = ω/N"). The beams
    turn toward the vertical as ω → N. Parameters: omega [rad/s] (0 ≤ ω ≤ N); N [rad/s]. Returns θ [rad] (NaN for
    ω > N: no propagating internal wave). Validation: V6, V7 — tests/test_ch07.py:
    test_internal_wave_omega_V7_limits_and_direction_only, test_book_V6_nonlinear_interface_and_internal. Checks: V1 ω/N
    = 0.71 ⇒ 44.77°. Label: analytic.
    """
    r = _F(omega) / float(N)
    with np.errstate(invalid="ignore"):
        return _S(np.where((r >= 0) & (r <= 1), np.arccos(np.clip(r, 0.0, 1.0)), np.nan))  # Eq. (7.139) inverted


def internal_wave_velocities(k, m, N: float, printed: bool = False) -> dict:
    """Phase and group velocity vectors of an internal wave in the x–z plane (l = 0) [m/s].

    c = (ω/K²)(k e_x + m e_z) (7.144); c_g = ∇_K ω = (N m/K³)(sgn(k) m e_x − |k| e_z) (sign-safe; for k > 0 this is the
    book's (7.145) (Nm/K³)(m e_x − k e_z)); c·c_g = 0 (7.146); horizontal components agree in sign, vertical ones are equal
    and opposite. ``printed=True`` evaluates the printed (7.145) literally (wrong for k < 0 — Fig. 7.29 draws k < 0).
    Book: §7.8, Eqs. (7.143)–(7.146), Figs. 7.29, 7.31.
    Parameters: k horizontal, m vertical wavenumber [rad/m]; N [rad/s].
    Returns dict(omega, c (2,), cg (2,), dot (c·c_g), theta_K (angle of K above the horizontal, from cos θ = |k|/K)).
    Validation: V1, V2 — tests/test_ch07.py: test_internal_group_velocity_V2_derivation,
    test_internal_wave_velocities_V1_gradient_parity_both_signs, test_linear_evolve_2d_V1_plane_wave_and_packet_at_cg.
    Checks: V1 equals :func:`group_velocity_vector` of :func:`internal_wave_omega` for both signs of k; c·c_g = 0 to
    round-off; c_z = −c_gz. Label: analytic.
    """
    k_, m_ = float(k), float(m)
    K = np.hypot(k_, m_)
    om = float(internal_wave_omega(k_, m_, N))
    c = om / K ** 2 * np.array([k_, m_])  # Eq. (7.144)
    if printed:
        cg = float(N) * m_ / K ** 3 * np.array([m_, -k_])  # Eq. (7.145) as printed (k > 0 only)
    else:
        cg = float(N) * m_ / K ** 3 * np.array([np.sign(k_) * m_, -abs(k_)])  # ∇_K ω, sign-safe (7.145)
    return {"omega": om, "c": c, "cg": cg, "dot": float(np.dot(c, cg)),  # Eq. (7.146)
            "theta_K": float(np.arccos(abs(k_) / K))}
