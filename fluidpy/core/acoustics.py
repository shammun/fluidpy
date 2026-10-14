"""Linear acoustics and simple waves (Kundu, Cohen & Dowling 5e, §15.2 and the wave-front picture of §15.9 — "AC" in
the analysis and curation files).

Sound is the linear limit of compressible flow: small isentropic disturbances p′, ρ′, u′ on a uniform reference state
(p₀, ρ₀, T₀, stream U) obey p′ ≅ c²ρ′ (15.5) and the convected wave equation (15.12); in one dimension every solution is
a right-going shape plus a left-going shape (d'Alembert, (15.15)–(15.17)).  The module also holds the finite-amplitude
simple wave (steepening, **ours**) and the wave fronts of a moving point source (Mach cone).

Conventions
-----------
* SI units; **subscript 0 here is the undisturbed reference state, not the stagnation state** (trap T1); primes are
  disturbances; ``c`` is the sound speed of the reference state; ``U`` > 0 is a stream in +x.
* ``f`` is the right-going shape, ``g`` the left-going shape (callables of one variable returning pressure [Pa]);
  ``None`` means "absent".
* Scalars in → floats out; arrays broadcast.

Reusable beyond Ch. 15: the same wave equation carries the arterial pressure pulse of Ch. 16 (p′ = ρ₀c u′ is the
water-hammer relation), and :func:`simple_wave_evolve` with γ = 2, c₀ = √(gH) is the shallow-water simple wave of ch07.
"""
from __future__ import annotations

from typing import Callable

import numpy as np
from scipy.integrate import quad

from ._util import as_scalar_if_0d

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d

__all__ = [
    "isentropic_pressure_perturbation", "acoustic_density", "acoustic_temperature", "condensation",
    "linearisation_error", "newton_sound_speed", "acoustic_impedance", "sound_pressure_level",
    "convected_wave_residual", "wave_equation_residual", "wave_equation_step",
    "dalembert", "plane_wave_velocity", "velocity_from_pressure", "dalembert_fg", "dalembert_initial_value",
    "dalembert_in_stream", "wave_speeds", "upstream_signal_possible", "pulse", "PULSE_KINDS", "xt_characteristics",
    "acoustic_energy", "simple_wave_speed", "simple_wave_evolve", "breaking_time",
    "source_wavefronts", "doppler_factor", "mach_cone_delay", "in_zone_of_action",
]


def _zero(s):
    return np.zeros_like(_F(s))


# ======================================================================================================================
# small amplitude: p′, ρ′, T′ and the size of what linearisation drops
# ======================================================================================================================
def isentropic_pressure_perturbation(rho_prime, rho0, p0, gamma: float = 1.4, order=1):
    """Pressure disturbance p′ that goes with a density disturbance ρ′ at constant entropy.

    Book: §15.2, the Taylor expansion p = p₀ + c²ρ′ + ½(∂²p/∂ρ²)_s ρ′² + … (unnumbered) and Eq. (15.5) p′ ≅ c²ρ′,
    with c² = γp₀/ρ₀ (15.6).  For a perfect gas (15.2) the second coefficient is (∂²p/∂ρ²)_s = γ(γ − 1)p₀/ρ₀².
    Parameters: rho_prime ρ′ [kg/m³]; rho0 [kg/m³]; p0 [Pa]; gamma [–]; order — 1 (acoustics), 2 (with the quadratic
    term) or "exact" (p₀[(1 + ρ′/ρ₀)^γ − 1]).
    Returns p′ [Pa].
    Assumptions: perfect gas, isentropic disturbance.
    Validation: V1 — tests/test_ch15.py: test_acoustic_relation_V1_linearisation_error_is_the_quadratic_term.
    Label: analytic.
    """
    r = _F(rho_prime)
    c2 = gamma * _F(p0) / _F(rho0)  # Eq. (15.6)
    if order == 1:
        return _S(c2 * r)  # Eq. (15.5)
    if order == 2:
        return _S(c2 * r + 0.5 * gamma * (gamma - 1.0) * _F(p0) / _F(rho0) ** 2 * r ** 2)
    if order == "exact":
        return _S(_F(p0) * np.expm1(gamma * np.log1p(r / _F(rho0))))  # Eq. (15.2)
    raise ValueError('order must be 1, 2 or "exact"')


def acoustic_density(p_prime, c):
    """Density disturbance of a sound wave, ρ′ = p′/c².

    Book: §15.2, Eq. (15.5)  p′ ≅ c²ρ′.
    Parameters: p_prime [Pa]; c sound speed of the reference state [m/s].   Returns ρ′ [kg/m³].
    Assumptions: small isentropic disturbance (15.7).
    Validation: V1 — tests/test_ch15.py: test_acoustic_relation_V1_linearisation_error_is_the_quadratic_term.
    Label: analytic.
    """
    return _S(_F(p_prime) / _F(c) ** 2)  # Eq. (15.5)


def acoustic_temperature(p_prime, p0, T0, gamma: float = 1.4):
    """Temperature disturbance of a sound wave, T′ = T₀ [(γ − 1)/γ] p′/p₀ (ours).

    Book: §15.2 (used in words: "if p′ > 0 then T′ > 0, so the sound speed increases within a compression") — the
    first-order form of the isentropic relation (15.2) T ∝ p^{(γ−1)/γ}.
    Parameters: p_prime [Pa]; p0 reference pressure [Pa]; T0 reference temperature [K]; gamma [–].
    Returns T′ [K].
    Assumptions: perfect gas, small isentropic disturbance.
    Validation: V1 — tests/test_ch15.py: test_acoustic_relation_V1_linearisation_error_is_the_quadratic_term.
    Label: analytic.
    """
    return _S(_F(T0) * (gamma - 1.0) / gamma * _F(p_prime) / _F(p0))


def condensation(p_prime, rho0, c):
    """Condensation ρ′/ρ₀ = p′/(ρ₀c²) — the small parameter of acoustics.

    Book: §15.2, Eq. (15.7)  p′/ρ₀c² ≪ 1 (and ρ₀c² = γp₀).
    Parameters: p_prime [Pa]; rho0 [kg/m³]; c [m/s].   Returns the condensation [–] (also equal to u′/c in a plane
    travelling wave).
    Assumptions: none (a definition).
    Validation: V1 — tests/test_ch15.py: test_acoustic_relation_V1_linearisation_error_is_the_quadratic_term.
    Label: analytic.
    """
    return _S(_F(p_prime) / (_F(rho0) * _F(c) ** 2))  # Eq. (15.7)


def linearisation_error(p_prime, rho0, c, gamma: float = 1.4):
    """Relative size of the quadratic term that (15.5) drops: ½(γ − 1) × condensation (ours).

    Book: §15.2 (the second term of the Taylor expansion against the first): ½γ(γ − 1)(p₀/ρ₀²)ρ′²/(c²ρ′).
    Parameters: p_prime [Pa]; rho0 [kg/m³]; c [m/s]; gamma [–].   Returns the ratio [–].
    Assumptions: perfect gas.
    Validation: V1 — tests/test_ch15.py: test_acoustic_relation_V1_linearisation_error_is_the_quadratic_term.
    Label: analytic.
    """
    return _S(0.5 * (gamma - 1.0) * _F(condensation(p_prime, rho0, c)))


def newton_sound_speed(T, R: float):
    """Newton's isothermal estimate of the sound speed, c = √(RT) — too low by the factor 1/√γ.

    Book: §15.2 (Newton assumed p′/p₀ = ρ′/ρ₀; Laplace's isentropic c = √(γRT) is (15.6)).
    Parameters: T [K]; R [J/(kg K)].   Returns c_Newton [m/s].
    Assumptions: isothermal compression (not what a sound wave does).
    Validation: V1 — tests/test_ch15.py: test_acoustic_relation_V1_linearisation_error_is_the_quadratic_term.
    Label: analytic.
    """
    return _S(np.sqrt(R * _F(T)))


def acoustic_impedance(rho0, c):
    """Specific acoustic impedance ρ₀c of a plane travelling wave: p′ = ±ρ₀c u′ (ours: the name).

    Book: §15.2, Eq. (15.16) (the factor 1/(ρ₀c)).
    Parameters: rho0 [kg/m³]; c [m/s].   Returns ρ₀c [Pa s/m] (about 415 for air at 15 °C).
    Assumptions: plane travelling wave in a still medium.
    Validation: V1 — tests/test_ch15.py: test_particle_velocity_V1_equals_the_time_integral_of_the_pressure_gradient.
    Label: analytic.
    """
    return _S(_F(rho0) * _F(c))


def sound_pressure_level(p_rms, p_ref: float = 2.0e-5):
    """Sound pressure level 20 log₁₀(p_rms/p_ref) in decibels (ours; not in the book).

    Book: none (used to give §15.2's "order 1 Pa or less" a familiar scale).
    Parameters: p_rms root-mean-square acoustic pressure [Pa]; p_ref reference pressure [Pa] (2e-5 Pa, the standard
    reference in air).   Returns the level [dB].
    Assumptions: none.
    Validation: V1 — tests/test_ch15.py: test_acoustic_relation_V1_linearisation_error_is_the_quadratic_term.
    Label: analytic.
    """
    return _S(20.0 * np.log10(_F(p_rms) / p_ref))


# ======================================================================================================================
# the wave equation
# ======================================================================================================================
def convected_wave_residual(p_fn: Callable, x, t, c: float, U: float = 0.0, h: float | None = None):
    """Residual of the one-dimensional convected wave equation for a given field p′(x, t).

    Book: §15.2, Eq. (15.12)  (1/c²)(∂/∂t + U_j∂/∂x_j)²p′ − ∂²p′/∂x_j∂x_j = 0, with U_j = (U, 0, 0):
    (1/c²)(p_tt + 2U p_xt + U² p_xx) − p_xx.
    Parameters: p_fn callable p′(x, t) [Pa]; x [m]; t [s]; c [m/s]; U stream speed [m/s]; h space step of the
    difference stencil [m] (default 1e-3 × max(1, |x|)); the time step is h/c.
    Returns the residual [Pa/m²] (O(h²) for a solution).
    Method: explicit second-order central differences (three-point for p_tt, p_xx; four-point cross for p_xt).
    Assumptions: uniform stream and sound speed.
    Validation: V3 — tests/test_ch15.py: test_wave_equation_V3_residual_of_a_travelling_pulse_vanishes_at_second_order.
    Label: converged.
    """
    x, t = _F(x), _F(t)
    hx = 1e-3 * np.maximum(1.0, np.abs(x)) if h is None else float(h)
    ht = hx / c
    p0_ = _F(p_fn(x, t))
    pxx = (_F(p_fn(x + hx, t)) - 2.0 * p0_ + _F(p_fn(x - hx, t))) / hx ** 2
    ptt = (_F(p_fn(x, t + ht)) - 2.0 * p0_ + _F(p_fn(x, t - ht))) / ht ** 2
    pxt = (_F(p_fn(x + hx, t + ht)) - _F(p_fn(x + hx, t - ht)) - _F(p_fn(x - hx, t + ht))
           + _F(p_fn(x - hx, t - ht))) / (4.0 * hx * ht)
    return _S((ptt + 2.0 * U * pxt + U ** 2 * pxx) / c ** 2 - pxx)  # Eq. (15.12)


def wave_equation_residual(p_fn: Callable, x, t, c: float, h: float | None = None):
    """Residual of the one-dimensional wave equation in a still fluid, (1/c²)∂²p′/∂t² − ∂²p′/∂x².

    Book: §15.2, Eq. (15.13).
    Parameters: p_fn callable p′(x, t) [Pa]; x [m]; t [s]; c [m/s]; h stencil step [m].
    Returns the residual [Pa/m²].
    Assumptions: uniform still medium.
    Validation: V3 — tests/test_ch15.py: test_wave_equation_V3_residual_of_a_travelling_pulse_vanishes_at_second_order.
    Label: converged.
    """
    return convected_wave_residual(p_fn, x, t, c, 0.0, h)  # Eq. (15.13)


def wave_equation_step(p_now, p_prev, C: float, bc: str = "periodic"):
    """One leapfrog step of the wave equation (15.13): p^{n+1}_i = 2p^n_i − p^{n−1}_i + C²(p^n_{i+1} − 2p^n_i +
    p^n_{i−1}) (**our scheme** — the book specifies none).

    Book: §15.2, Eq. (15.13).
    Parameters: p_now, p_prev — p′ at time levels n and n − 1 on a uniform grid (values at nodes) [Pa]; C = cΔt/Δx
    Courant number [–], must be ≤ 1 (raises ``ValueError`` otherwise); bc "periodic" or "fixed" (p′ = 0 at both
    ends).
    Returns p′ at level n + 1.  Second order in space and time; exact translation at C = 1.
    Assumptions: uniform still medium.
    Validation: V3 — tests/test_ch15.py: test_leapfrog_V3_second_order_and_exact_at_courant_one.
    Label: converged.
    """
    if C > 1.0 + 1e-12:
        raise ValueError(f"Courant number C = c dt/dx = {C:g} > 1: the leapfrog scheme is unstable")
    p, q = _F(p_now), _F(p_prev)
    # DEVIATION: the leapfrog scheme is our choice — the book states the wave equation but gives no numerical scheme.
    lap = np.roll(p, -1) - 2.0 * p + np.roll(p, 1)
    out = 2.0 * p - q + C ** 2 * lap
    if bc == "fixed":
        out[0] = out[-1] = 0.0
    elif bc != "periodic":
        raise ValueError('bc must be "periodic" or "fixed"')
    return out


# ======================================================================================================================
# d'Alembert
# ======================================================================================================================
def dalembert(x, t, c: float, f: Callable | None = None, g: Callable | None = None):
    """d'Alembert's solution of the wave equation: p′(x, t) = f(x − ct) + g(x + ct).

    Book: §15.2, Eq. (15.15); Fig. 15.1.
    Parameters: x [m]; t [s]; c [m/s]; f right-going shape, g left-going shape (callables of one variable [m] → [Pa];
    ``None`` = absent).
    Returns p′ [Pa].
    Assumptions: one-dimensional linear acoustics in a uniform still fluid.
    Validation: V1, V4 — tests/test_ch15.py: test_dalembert_V1_shapes_travel_unchanged,
    test_acoustic_energy_V4_is_constant_while_two_pulses_cross.
    Label: analytic, conserved.
    """
    x, t = _F(x), _F(t)
    f = _zero if f is None else f
    g = _zero if g is None else g
    return _S(_F(f(x - c * t)) + _F(g(x + c * t)))  # Eq. (15.15)


def plane_wave_velocity(x, t, c: float, rho0: float, f: Callable | None = None, g: Callable | None = None):
    """Fluid (particle) velocity of the two travelling parts: u′(x, t) = [f(x − ct) − g(x + ct)]/(ρ₀c).

    Book: §15.2, Eq. (15.16).  A compression (f > 0 or g > 0) moves the fluid in the direction the wave travels; an
    expansion moves it the other way.  u′ is not the wave speed: u′/c equals the condensation.
    Parameters: x [m]; t [s]; c [m/s]; rho0 [kg/m³]; f, g as in :func:`dalembert`.
    Returns u′ [m/s].
    Assumptions: as :func:`dalembert`; fluid at rest before the wave arrives.
    Validation: V1, V4 — tests/test_ch15.py:
    test_particle_velocity_V1_equals_the_time_integral_of_the_pressure_gradient,
    test_acoustic_energy_V4_is_constant_while_two_pulses_cross.
    Label: analytic, conserved.
    """
    x, t = _F(x), _F(t)
    f = _zero if f is None else f
    g = _zero if g is None else g
    return _S((_F(f(x - c * t)) - _F(g(x + c * t))) / (rho0 * c))  # Eq. (15.16)


def velocity_from_pressure(p_fn: Callable, x: float, t: float, rho0: float, *, t0: float = 0.0, h: float = 1e-4,
                           limit: int = 200) -> float:
    """Particle velocity from the pressure field by integrating the linear momentum law in time.

    Book: §15.2, Eq. (15.14)  u′(x, t) = −(1/ρ₀) ∫ (∂p′/∂x) dt  — an indefinite integral; the constant is fixed by
    u′ = 0 at ``t0`` (a time before the wave reaches x).
    Parameters: p_fn callable p′(x, t) [Pa]; x [m] and t [s] (scalars); rho0 [kg/m³]; t0 [s] start of the
    integration, where the fluid is at rest at x; h step of the central difference for ∂p′/∂x [m]; limit — ``quad``
    subdivisions.
    Returns u′ [m/s].
    Assumptions: linear acoustics (15.11) in a still fluid.
    Validation: V1 — tests/test_ch15.py: test_particle_velocity_V1_equals_the_time_integral_of_the_pressure_gradient.
    Label: analytic.
    """
    dpdx = lambda s: (float(p_fn(x + h, s)) - float(p_fn(x - h, s))) / (2.0 * h)  # noqa: E731
    val, _ = quad(dpdx, t0, t, limit=limit)
    return float(-val / rho0)  # Eq. (15.14)


def dalembert_fg(F: Callable, G: Callable, c: float, x_ref: float = 0.0) -> tuple[Callable, Callable]:
    """The two travelling shapes that satisfy given initial conditions p′ = F(x), ∂p′/∂t = G(x) at t = 0.

    Book: Exercise 15.1(c):  f(x) = ½[F(x) − (1/c)∫₀ˣ G dx̂],  g(x) = ½[F(x) + (1/c)∫₀ˣ G dx̂].
    Parameters: F [Pa] and G [Pa/s] callables of x [m]; c [m/s]; x_ref lower limit of the integral [m] (any value:
    the additive constant cancels in f + g).
    Returns (f, g) callables (vectorised by a loop over ``quad``).
    Assumptions: as :func:`dalembert`.
    Validation: V2, V1 — tests/test_ch15.py: test_dalembert_initial_value_V1_reproduces_both_initial_conditions.
    Label: symbolic, analytic.
    """
    def integral(s):
        return np.vectorize(lambda v: quad(lambda z: float(G(z)), x_ref, float(v))[0], otypes=[float])(_F(s))

    def f(s):
        return 0.5 * (_F(F(s)) - integral(s) / c)

    def g(s):
        return 0.5 * (_F(F(s)) + integral(s) / c)

    return f, g


def dalembert_initial_value(x, t, c: float, F: Callable, G: Callable | None = None):
    """Solution of the wave equation from initial data: p′(x, t) = ½[F(x − ct) + F(x + ct)] + (1/2c)∫_{x−ct}^{x+ct} G.

    Book: Exercise 15.1(c) (the f and g of :func:`dalembert_fg` inserted in (15.15)).
    Parameters: x [m]; t [s]; c [m/s]; F initial pressure p′(x, 0) [Pa]; G initial rate ∂p′/∂t(x, 0) [Pa/s]
    (``None`` = 0: a disturbance released from rest splits into two equal halves).
    Returns p′ [Pa].
    Assumptions: as :func:`dalembert`.
    Validation: V2, V1 — tests/test_ch15.py: test_dalembert_initial_value_V1_reproduces_both_initial_conditions.
    Label: symbolic, analytic.
    """
    x, t = np.broadcast_arrays(_F(x), _F(t))
    a, b = x - c * t, x + c * t
    out = 0.5 * (_F(F(a)) + _F(F(b)))
    if G is not None:
        out = out + np.vectorize(lambda lo, hi: quad(lambda z: float(G(z)), lo, hi)[0], otypes=[float])(a, b) / (2.0 * c)
    return _S(out)


def dalembert_in_stream(x, t, c: float, U: float, f: Callable | None = None, g: Callable | None = None):
    """Plane waves on a uniform stream: p′(x, t) = f(x − (c + U)t) + g(x + (c − U)t).

    Book: §15.2, Eq. (15.17).  For U > c both parts travel downstream: nothing propagates upstream.
    Parameters: x [m]; t [s]; c [m/s]; U stream speed [m/s] (+x); f, g as in :func:`dalembert` (shapes at t = 0).
    Returns p′ [Pa].
    Assumptions: one-dimensional linear acoustics on a uniform stream.
    Validation: V1 — tests/test_ch15.py: test_waves_in_a_stream_V1_galilean_shift_of_dalembert.
    Label: analytic.
    """
    x, t = _F(x), _F(t)
    f = _zero if f is None else f
    g = _zero if g is None else g
    return _S(_F(f(x - (c + U) * t)) + _F(g(x + (c - U) * t)))  # Eq. (15.17)


def wave_speeds(U, c):
    """Speeds of the two wave families seen from the ground: (U + c, U − c).

    Book: §15.2, Eq. (15.17).
    Parameters: U stream speed [m/s]; c sound speed [m/s].   Returns (U + c, U − c) [m/s]; the second is negative
    (upstream-running) only while U < c.
    Assumptions: uniform stream.
    Validation: V1 — tests/test_ch15.py: test_waves_in_a_stream_V1_galilean_shift_of_dalembert.
    Label: analytic.
    """
    return _S(_F(U) + _F(c)), _S(_F(U) - _F(c))


def upstream_signal_possible(U, c):
    """Can a pressure signal travel upstream?  True while |U| < c.

    Book: §15.2 (text after (15.17): "in supersonic flow only upstream disturbances may influence the flow").
    Parameters: U [m/s]; c [m/s].   Returns bool (array).
    Assumptions: uniform stream, linear waves.
    Validation: V1 — tests/test_ch15.py: test_waves_in_a_stream_V1_galilean_shift_of_dalembert.
    Label: analytic.
    """
    r = np.abs(_F(U)) < _F(c)
    return bool(r) if r.ndim == 0 else r


PULSE_KINDS = ("gaussian", "tophat", "packet", "sine")


def pulse(kind: str = "gaussian", width: float = 1.0, amplitude: float = 1.0, x0: float = 0.0) -> Callable:
    """A pulse shape s ↦ p′ for the d'Alembert functions (shared by the notebook and the sound-pulse explainer; ours).

    With ξ = (s − x0)/width and a = amplitude:
    "gaussian"  a exp(−ξ²);   "tophat"  ½a[tanh(10(ξ + ½)) − tanh(10(ξ − ½))] (a smooth box of length ``width``);
    "packet"  a exp(−ξ²) cos(6ξ) (compression and expansion in equal amounts, as in Fig. 15.1);
    "sine"  a sin(2πξ) for |ξ| ≤ ½, else 0 (one cycle).
    Book: §15.2, Fig. 15.1 (a generic waveform; the shapes are ours).
    Parameters: kind (one of :data:`PULSE_KINDS`); width [m]; amplitude [Pa]; x0 centre [m].
    Returns a callable of one variable [m] → [Pa] (vectorised).
    Assumptions: none.
    Validation: V1 — the four shapes at their defining points (a helper of the d'Alembert tests) — tests/test_ch15.py:
    test_dalembert_V1_shapes_travel_unchanged.
    Label: analytic.
    """
    if kind not in PULSE_KINDS:
        raise ValueError(f"kind must be one of {PULSE_KINDS}")

    def shape(s):
        xi = (_F(s) - x0) / width
        if kind == "gaussian":
            return amplitude * np.exp(-xi ** 2)
        if kind == "tophat":
            return 0.5 * amplitude * (np.tanh(10.0 * (xi + 0.5)) - np.tanh(10.0 * (xi - 0.5)))
        if kind == "packet":
            return amplitude * np.exp(-xi ** 2) * np.cos(6.0 * xi)
        return np.where(np.abs(xi) <= 0.5, amplitude * np.sin(2.0 * np.pi * xi), 0.0)

    return shape


def xt_characteristics(U: float, c: float, x0, t) -> dict:
    """The two families of lines in the x–t plane along which sound travels: x = x₀ + (U ± c)t.

    Book: §15.2, Eq. (15.17) (the arguments of f and g are constant on these lines).
    Parameters: U [m/s]; c [m/s]; x0 starting positions [m] (scalar or (n,)); t times [s] (scalar or (m,)).
    Returns dict: ``t`` (m,), ``x_plus`` (n, m) (the f family, speed U + c), ``x_minus`` (n, m) (the g family, speed
    U − c), ``speed_plus``, ``speed_minus`` [m/s].
    Assumptions: uniform stream.
    Validation: V1 — tests/test_ch15.py: test_waves_in_a_stream_V1_galilean_shift_of_dalembert.
    Label: analytic.
    """
    x0, t = np.atleast_1d(_F(x0)), np.atleast_1d(_F(t))
    return dict(t=t, x_plus=x0[:, None] + (U + c) * t[None, :], x_minus=x0[:, None] + (U - c) * t[None, :],
                speed_plus=U + c, speed_minus=U - c)


def acoustic_energy(p_prime, u_prime, x, rho0: float, c: float) -> float:
    """Acoustic energy per unit cross-section, ∫ ½[p′²/(ρ₀c²) + ρ₀u′²] dx (ours; not in the book).

    Book: none (a conserved quantity of (15.10)–(15.11) used as an independent check of (15.16)).
    Parameters: p_prime [Pa] and u_prime [m/s] sampled at x [m] (uniform grid, compact disturbance); rho0 [kg/m³];
    c [m/s].   Returns the energy [J/m²] (composite trapezoid).
    Assumptions: linear acoustics in a still fluid.
    Validation: V4 — tests/test_ch15.py: test_acoustic_energy_V4_is_constant_while_two_pulses_cross.
    Label: conserved.
    """
    e = 0.5 * (_F(p_prime) ** 2 / (rho0 * c ** 2) + rho0 * _F(u_prime) ** 2)
    x = _F(x)
    return float(np.sum(0.5 * (e[1:] + e[:-1]) * np.diff(x)))


# ======================================================================================================================
# finite amplitude: the simple wave (ours — the book argues in words)
# ======================================================================================================================
def simple_wave_speed(u, c0: float, gamma: float = 1.4):
    """Speed of a wavelet of a right-going finite-amplitude simple wave: dx/dt = u + c = c₀ + ½(γ + 1)u (**ours**).

    Book: §15.2 (in words: the sound speed rises inside a compression, so compressions steepen and expansions
    flatten).  In a simple wave running into still gas the Riemann invariant u − 2c/(γ − 1) = −2c₀/(γ − 1) is
    uniform, so c = c₀ + ½(γ − 1)u.  γ = 2 with c₀ = √(gH) is the shallow-water wave of ch07.
    Parameters: u particle velocity [m/s]; c0 sound speed of the undisturbed gas [m/s]; gamma [–].
    Returns the wavelet speed [m/s].
    Assumptions: perfect gas, isentropic (before any shock forms), waves of one family only.
    Validation: V4, V1 — tests/test_ch15.py: test_simple_wave_V4_riemann_invariant_and_isentrope,
    test_breaking_time_V1_first_crossing_of_neighbouring_characteristics.
    Label: analytic, conserved.
    """
    return _S(c0 + 0.5 * (gamma + 1.0) * _F(u))


def _ddx(v: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Second-order derivative on a (possibly non-uniform) grid, explicit stencils (one-sided at the ends)."""
    d = np.empty_like(v)
    hl, hr = x[1:-1] - x[:-2], x[2:] - x[1:-1]
    d[1:-1] = (-hr / (hl * (hl + hr)) * v[:-2] + (hr - hl) / (hl * hr) * v[1:-1] + hl / (hr * (hl + hr)) * v[2:])
    h0, h1 = x[1] - x[0], x[2] - x[1]
    d[0] = -(2 * h0 + h1) / (h0 * (h0 + h1)) * v[0] + (h0 + h1) / (h0 * h1) * v[1] - h0 / (h1 * (h0 + h1)) * v[2]
    h0, h1 = x[-1] - x[-2], x[-2] - x[-3]
    d[-1] = (2 * h0 + h1) / (h0 * (h0 + h1)) * v[-1] - (h0 + h1) / (h0 * h1) * v[-2] + h0 / (h1 * (h0 + h1)) * v[-3]
    return d


def breaking_time(u0, x, gamma: float = 1.4, *, max_compression_rate: float | None = None) -> float:
    """Time at which a right-going simple wave first becomes multivalued (a shock starts to form):
    t_b = 2/[(γ + 1) max(−∂u₀/∂x)] (**ours**).

    Book: §15.2 (steepening described in words; "the limiting form … is a shock wave").  Wavelets move at
    c₀ + ½(γ + 1)u₀(x₀); two neighbours meet when 1 + ½(γ + 1)(∂u₀/∂x₀)t = 0.
    Parameters: u0 initial particle velocity [m/s] sampled at x [m] (n ≥ 3); gamma [–]; max_compression_rate — give
    max(−∂u₀/∂x) [1/s] analytically instead (then u0, x are ignored).
    Returns t_b [s]; ``inf`` for a pure expansion (no part of the wave steepens).
    Method: explicit second-order differences of the samples.
    Assumptions: as :func:`simple_wave_speed`.
    Validation: V4, V1 — tests/test_ch15.py: test_simple_wave_V4_riemann_invariant_and_isentrope,
    test_breaking_time_V1_first_crossing_of_neighbouring_characteristics.
    Label: analytic, conserved.
    """
    s = float(max_compression_rate) if max_compression_rate is not None else float(np.max(-_ddx(_F(u0), _F(x))))
    return float(2.0 / ((gamma + 1.0) * s)) if s > 0.0 else float("inf")


def simple_wave_evolve(u0, x, t, c0: float, gamma: float = 1.4) -> dict:
    """Steepening of a right-going finite-amplitude wave by characteristics: each wavelet keeps its u and moves at
    c₀ + ½(γ + 1)u₀, so x(t) = x₀ + [c₀ + ½(γ + 1)u₀(x₀)]t (**ours**).

    Book: §15.2 (words only; a sine wave tends to a sawtooth).  Generalises ch07 ``simple_wave_evolve`` (γ = 2,
    c₀ = √(gH), u = 2(√(g(H + η)) − c₀)).
    Parameters: u0 initial particle velocity [m/s] at x (n,) [m]; t scalar or (m,) [s]; c0 [m/s]; gamma [–].
    Returns dict: ``x_points`` ((n,) or (m, n): where each wavelet is; also under ``x``), ``u`` (n,) [m/s]
    (unchanged along characteristics), ``c`` (n,) local sound speed c₀ + ½(γ − 1)u [m/s], ``speed`` (n,) wavelet
    speed [m/s], ``p_over_p0`` and ``rho_over_rho0`` (n,) [–] from the isentrope (15.2), ``t_break`` [s], ``valid``
    (bool or (m,): t < t_break — beyond it the curve is multivalued and a shock must be fitted).
    Assumptions: as :func:`simple_wave_speed`; gas ahead at rest.
    Validation: V4, V1 — tests/test_ch15.py: test_simple_wave_V4_riemann_invariant_and_isentrope,
    test_simple_wave_V1_gamma_two_is_the_shallow_water_wave_of_chapter_7.
    Label: analytic, conserved.
    """
    u0, x_, t_ = _F(u0), _F(x), _F(t)
    speed = c0 + 0.5 * (gamma + 1.0) * u0
    c = c0 + 0.5 * (gamma - 1.0) * u0
    tb = breaking_time(u0, x_, gamma)
    xs = x_[None, :] + speed[None, :] * np.atleast_1d(t_)[:, None]
    xp = xs[0] if t_.ndim == 0 else xs
    valid = t_ < tb
    return dict(x_points=xp, x=xp, u=u0, c=c, speed=speed, p_over_p0=(c / c0) ** (2.0 * gamma / (gamma - 1.0)),
                rho_over_rho0=(c / c0) ** (2.0 / (gamma - 1.0)), t_break=tb,
                valid=bool(valid) if valid.ndim == 0 else valid)


# ======================================================================================================================
# §15.9 — wave fronts of a moving source, Mach cone
# ======================================================================================================================
def source_wavefronts(U: float, c: float, t: float, n_fronts: int = 6, dt: float = 1.0) -> dict:
    """Wave fronts sent out at equal time intervals by a point source moving at constant speed.

    Book: §15.9, Fig. 15.17: the front emitted at time t_e is the circle |x − x_s(t_e)| = c(t − t_e).  Here the source
    moves in **+x**, x_s(t) = Ut; fronts are emitted at t_e = t − k·dt, k = 1 … n_fronts.
    Parameters: U source speed [m/s] ≥ 0; c sound speed [m/s]; t observation time [s]; n_fronts; dt emission
    interval [s].
    Returns dict: ``x_centres`` (n,) [m], ``y_centres`` (zeros), ``radii`` (n,) [m], ``t_emit`` (n,) [s],
    ``x_source`` = Ut [m], ``M`` = U/c, ``mu`` Mach angle [rad] (NaN for M < 1), ``ahead`` (n,) — how far each front
    reaches ahead of the source, radius − (x_source − centre) [m] (positive for U < c, negative for U > c: no front
    is ahead of a supersonic source).
    Assumptions: uniform still medium; weak disturbances; a plane through the source path.
    Validation: V1 — tests/test_ch15.py: test_mach_cone_V1_envelope_of_the_wave_fronts.
    Label: analytic.
    """
    k = np.arange(1, int(n_fronts) + 1, dtype=float)
    te = t - k * dt
    xc, r = U * te, c * k * dt
    M = U / c
    return dict(x_centres=xc, y_centres=np.zeros_like(xc), radii=r, t_emit=te, x_source=U * t, M=M,
                mu=float(np.arcsin(1.0 / M)) if M >= 1.0 else np.nan, ahead=r - (U * t - xc))


def doppler_factor(U, c, theta=0.0):
    """Frequency heard from a moving source over the frequency emitted, 1/(1 − M cos θ) (ours; not in the book).

    Book: none (the bunching of the fronts of Fig. 15.17b ahead of a subsonic source).
    Parameters: U source speed [m/s]; c [m/s]; theta angle between the source velocity and the direction to the
    observer at the emission time [rad].   Returns the factor [–] (inf on the Mach cone, where M cos θ = 1).
    Assumptions: uniform still medium, observer at rest.
    Validation: V1 — tests/test_ch15.py: test_mach_cone_V1_envelope_of_the_wave_fronts.
    Label: analytic.
    """
    with np.errstate(divide="ignore"):
        return _S(1.0 / (1.0 - _F(U) / _F(c) * np.cos(_F(theta))))


def mach_cone_delay(h, M, c):
    """Time between a supersonic aircraft passing overhead and its Mach cone reaching the observer:
    t = h/(U tan μ) = h√(M² − 1)/(Mc) (ours).

    Book: §15.9 ("this arrival occurs after the aircraft has passed overhead") with Eq. (15.49).
    Parameters: h height of the flight path above the observer [m]; M flight Mach number [–] > 1; c [m/s].
    Returns the delay [s] (NaN for M < 1: a subsonic aircraft is heard before it arrives).
    Assumptions: level flight at constant speed in a uniform still atmosphere.
    Validation: V1 — tests/test_ch15.py: test_mach_cone_V1_envelope_of_the_wave_fronts.
    Label: analytic.
    """
    M = _F(M)
    with np.errstate(invalid="ignore"):
        return _S(np.where(M >= 1.0, _F(h) * np.sqrt(np.where(M >= 1.0, M ** 2 - 1.0, np.nan)) / (M * _F(c)), np.nan))


def in_zone_of_action(x, y, M):
    """Is the point (x, y) inside the Mach cone of a disturbance at the origin in a stream moving in +x?

    Book: §15.9, Fig. 15.17c (zone of action inside the cone, zone of silence outside), Eq. (15.49):
    inside ⇔ x ≥ |y|√(M² − 1), i.e. |y| ≤ x tan μ.
    Parameters: x, y position relative to the source [m] (x downstream); M stream Mach number [–].
    Returns bool (array); everywhere True for M < 1 (a subsonic source is felt in the whole field).
    Assumptions: uniform stream; weak disturbances.
    Validation: V1 — tests/test_ch15.py: test_mach_cone_V1_envelope_of_the_wave_fronts.
    Label: analytic.
    """
    x, y, M = np.broadcast_arrays(_F(x), _F(y), _F(M))
    with np.errstate(invalid="ignore"):
        r = np.where(M < 1.0, True, x >= np.abs(y) * np.sqrt(np.where(M >= 1.0, M ** 2 - 1.0, 0.0)))
    return bool(r) if r.ndim == 0 else r
