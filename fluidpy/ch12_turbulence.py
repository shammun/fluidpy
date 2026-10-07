"""Chapter 12 — Turbulence (Kundu, Cohen & Dowling, *Fluid Mechanics*, 5th ed., §12.1–§12.13).

The chapter's physics as functions.  Reusable pieces live in two core modules written for this chapter:
``fluidpy.core.turbstats`` (**TS**: moments, averages, correlations, spectra, synthetic fields — §12.3–§12.6) and
``fluidpy.core.wall_turbulence`` (**WT**: wall units, law of the wall, log law, composite profiles, friction laws, surface
layer — §12.9, §12.11).  This module holds everything that is specific to the chapter:

* §12.3  Example 12.1 (window factors of a time average);
* §12.5  mean stress tensor, RANS residuals, fluxes, the displaced-parcel argument for mean(uv) < 0, the closure count,
         and the symbolic Reynolds-averaging engines (``*_sympy``);
* §12.6  isotropic two-point tensor, f and g, dissipation in its isotropic forms;
* §12.7  the two energy budgets, Kolmogorov scales, the cascade, inertial-range and model spectra;
* §12.8  plane-jet similarity (velocity, cross-flow, stress, entrainment, scalar), exponents of the free shear flows;
* §12.10 eddy viscosity, mixing length (wall profile, channel), k–ε pieces;
* §12.11 Richardson numbers, Monin–Obukhov length, the stratified surface layer, temperature spectra;
* §12.12 Taylor's theory of dispersion, Langevin particles, random walk, eddy diffusivity.

Conventions (see also :func:`conventions`)
------------------------------------------
* Over-bar = **ensemble average** unless a name says otherwise; lower-case u, v, w and primed T′ are *fluctuations*.
* ``kappa`` = von Kármán constant; ``kappa_th`` = thermal diffusivity; ``e`` = turbulent kinetic energy per unit mass
  (the "k" of k–ε); ``eps`` = its dissipation rate; ``k1``, ``K`` = wavenumbers; ``u2`` = variance of **one** velocity
  component; ``h`` = full channel height, ``delta`` = half-height or layer thickness.
* Spectra are two-sided in angular frequency or wavenumber (book normalisation) unless ``one_sided``/``two_sided`` says
  otherwise; the three-dimensional spectrum E(K) is defined on K ≥ 0 with ∫E dK = ē.
* **Stratification (§12.11).**  Computations use Kundu's ``Γ ≡ dT/dz`` (adiabatic value Γ_a = −g/C_p ≈ −9.8 K/km); the
  standard meteorological ``Γ ≡ −dT/dz`` is reported alongside wherever a verdict is given.  The "temperature" in every
  buoyancy term of the chapter is **potential** temperature; :func:`gradient_richardson_thermal` takes the in-situ gradient
  plus Γ_a.  Heat fluxes ``wT`` [K m/s] and ``H`` [W/m²] are positive upward.
* **Empirical constants are never silent defaults**: κ, B, the wake strength, jet and free-shear constants are required
  keywords or cited presets.  :data:`FREE_SHEAR_CONSTANTS` is empty on purpose.  Numbers quoted by the book stay in the
  private ``tests/book_values_ch12.json``.
* Printed slips of the book are coded in their corrected form; where a test must be able to show that the printed form
  fails, a ``printed=True`` switch reproduces it (:func:`inertial_spectrum_1d`, :func:`rans_eddy_viscosity_residual`,
  :func:`eddy_diffusivity_asymptote`).  :func:`book_slips` lists them all.
* **Re-exports.**  Every public name of ``core.turbstats`` and ``core.wall_turbulence`` is also a name of this module
  (``ch12.log_law is WT.log_law``).
* **Return shapes** follow the lesson design (``analysis/ch12_design.md`` Part C): floats, tuples and dicts of floats
  for scalar input, so that explainer parity rows can index them (``ch12.kolmogorov_scales(nu, eps)[0]``).
* Sampled quantities take a ``seed``; synthetic fields are kinematic (no cascade) — see ``core.turbstats``.
"""
from __future__ import annotations

from fractions import Fraction
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.integrate import cumulative_trapezoid, quad, solve_ivp
from scipy.special import erf, erfinv

from .core import turbstats as TS
from .core import wall_turbulence as WT
# Re-export every public name of the two core modules (the ch11 pattern): ``ch12.log_law`` *is* ``WT.log_law`` and
# ``ch12.autocorrelation`` *is* ``TS.autocorrelation``, so notebook cells and explainer parity rows need ``ch12`` only.
from .core.turbstats import *  # noqa: F401,F403
from .core.wall_turbulence import *  # noqa: F401,F403
from .core._util import as_scalar_if_0d
from .core.jets import jet_momentum_flux
from .core.stratification import adiabatic_lapse_rate, lapse_rate_stability
from .core.thermo import CP_AIR, G0

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d

#: Three-dimensional Kolmogorov constant C of E(K) = C ε^{2/3} K^{−5/3} — the usual value (the book quotes "about 1.5").
KOLMOGOROV_C: float = 1.5

#: Standard k–ε constants (Launder & Sharma 1974; https://www.cfd-online.com/Wiki/Standard_k-epsilon_model).
#: Five constants: C_mu, C_eps1, C_eps2, sigma_e (the "sigma_k" of the literature) and sigma_eps.
K_EPSILON_CONSTANTS: dict[str, float] = {"C_mu": 0.09, "C_eps1": 1.44, "C_eps2": 1.92, "sigma_e": 1.0, "sigma_eps": 1.3}

#: Constants of the free-shear similarity laws (amplitudes and half-widths of Table 12.1).  **Deliberately empty**: the
#: book's set is book-quoted (private JSON) and no public table has been confirmed yet.  Callers pass their own dict
#: ``{"C_U": …, "C_Y": …, "xi_half_U": …, "xi_half_Y": …}``; examples label theirs "illustrative".
FREE_SHEAR_CONSTANTS: dict[str, dict] = {}

FREE_SHEAR_FLOWS = ("plane_jet", "round_jet", "plane_wake", "round_wake", "plane_plume", "round_plume", "shear_layer")


def _root() -> Path:
    from .core.project import repo_root

    return repo_root()


# ======================================================================================================================
# §12.3  Example 12.1
# ======================================================================================================================
def time_average_exp_cos(t, A: float, B: float, tau: float, omega: float, window: float):
    """Sliding time average of ``u = A exp(−t/τ) + B cos ωt`` over a window Δt, in closed form.

    Book: §12.3, Example 12.1 (Eq. (12.2) applied to this signal):
    ``ū(t) = [sinh(Δt/2τ)/(Δt/2τ)] A e^{−t/τ} + [sin(ωΔt/2)/(ωΔt/2)] B cos ωt``.
    Parameters: t [s]; A, B amplitudes [unit of u]; tau decay time [s]; omega angular frequency [rad/s]; window Δt [s].
    Returns (average [unit], mean_factor [–] ≥ 1, fluctuation_factor [–] in (−0.22, 1]).  Scalar-callable.
    Reading: the decaying mean is recovered when Δt ≪ τ (first factor → 1) and the cosine is suppressed when ωΔt ≫ 1
    (second factor → 0, exactly 0 for a whole number of periods): a good window needs 1 ≪ ωΔt ≪ ωτ.
    Validation: V2 sympy integral of the signal over the window; V7 Δt → 0 and ωΔt = 2π; V1 against
    ``TS.time_average`` of a sampled signal.  Time average, not ensemble.
    Assumptions: the signal is exactly A exp(-t/tau) + B cos(omega t); centred window.
    """
    x = 0.5 * float(window) / float(tau)
    mean_factor = np.sinh(x) / x if abs(x) > 1e-6 else 1.0 + x * x / 6.0
    fluct_factor = float(np.sinc(omega * window / (2.0 * np.pi)))  # sin(ωΔt/2)/(ωΔt/2)
    tt = _F(t)
    avg = mean_factor * A * np.exp(-tt / tau) + fluct_factor * B * np.cos(omega * tt)  # Example 12.1
    return _S(avg), float(mean_factor), fluct_factor


# ======================================================================================================================
# §12.1, §12.4, §12.5  fields: divergence, frozen-field probe
# ======================================================================================================================
def _spectral_derivative(f: np.ndarray, dx: float, axis: int) -> np.ndarray:
    n = f.shape[axis]
    k = 2.0 * np.pi * np.fft.fftfreq(n, d=dx)
    if n % 2 == 0:
        k[n // 2] = 0.0  # the Nyquist mode has no well-defined odd derivative
    shape = [1] * f.ndim
    shape[axis] = n
    return np.real(np.fft.ifft(1j * k.reshape(shape) * np.fft.fft(f, axis=axis), axis=axis))


def _periodic_central(f: np.ndarray, dx: float, axis: int) -> np.ndarray:
    return (np.roll(f, -1, axis=axis) - np.roll(f, 1, axis=axis)) / (2.0 * dx)


def divergence_rms(u, v, dx: float, w=None, method: str = "spectral") -> float:
    """Root-mean-square divergence ``∂u/∂x + ∂v/∂y (+ ∂w/∂z)`` of a periodic velocity field.

    Book: §12.1 (a turbulent velocity field satisfies continuity; a random vector field does not) and §12.5,
    Eq. (12.28) ``∂u_i/∂x_i = 0`` for the fluctuation.
    Parameters: u, v (and w) periodic arrays [m/s] with x on the last axis, y on the second-to-last (z third-to-last);
    dx grid spacing [m] (same in every direction); method "spectral" (FFT derivative) or "central" (second-order).
    Returns the rms divergence [1/s]: at round-off for ``TS.synthetic_solenoidal_field``, of order u_rms/dx for noise.
    Validation: V4 solenoidal field ≈ 1e-12·(u_rms/dx); white noise O(1).
    Assumptions: periodic box, uniform grid.
    """
    comps = [_F(u), _F(v)] + ([] if w is None else [_F(w)])
    d = _spectral_derivative if method == "spectral" else _periodic_central
    nd = comps[0].ndim
    div = sum(d(c, float(dx), nd - 1 - i) for i, c in enumerate(comps))
    return float(np.sqrt(np.mean(div ** 2)))


def mean_divergence(samples_u, samples_v, dx: float, dy: float, method: str = "spectral", report: bool = False):
    """Divergence field ``∂U/∂x + ∂V/∂y`` of the ensemble-mean velocity of a two-dimensional ensemble.

    Book: §12.5, Eqs. (12.27) ``∂U_i/∂x_i = 0`` and (12.28) ``∂u_i/∂x_i = 0``: averaging commutes with differentiation
    (12.8) ``mean(∂ũ/∂x) = ∂ū/∂x``, so if every member is divergence-free, so are the mean and each fluctuation.
    Parameters
    ----------
    samples_u, samples_v : (N, ny, nx) periodic fields [m/s], members on axis 0, x on the last axis.
    dx, dy : grid spacings [m].
    method : "spectral" (exact for the fields of ``TS.synthetic_solenoidal_field``) or "central" (second-order
        differences: a spectrally solenoidal member then has a truncation-level divergence, but the commutation still
        holds to round-off).
    report : False → the field; True → a dict of summary numbers.
    Returns
    -------
    (ny, nx) array [1/s], the divergence of the mean field — or, with ``report=True``, dict(field, mean_rms,
    fluctuation_rms, member_rms [1/s], commutation_residual = max |div(mean) − mean(div)| [1/s], round-off for either
    operator because both are linear).
    Assumptions: periodic box, uniform grid; "mean" is the average over the members supplied.
    Validation: V4 round-off for a solenoidal ensemble (spectral); V1 the commutation residual is round-off for both
    operators; O(u_rms/dx) for white noise.
    """
    su, sv = _F(samples_u), _F(samples_v)
    d = _spectral_derivative if method == "spectral" else _periodic_central
    div = lambda a, b: d(a, dx, a.ndim - 1) + d(b, dy, b.ndim - 2)  # noqa: E731
    U, V = su.mean(axis=0), sv.mean(axis=0)
    field = div(U, V)                                               # Eq. (12.27): divergence of the mean field
    if not report:
        return field
    rms = lambda a: float(np.sqrt(np.mean(a ** 2)))  # noqa: E731
    div_members = div(su, sv)
    return dict(field=field, mean_rms=rms(field), fluctuation_rms=rms(div(su - U, sv - V)),  # Eq. (12.28)
                member_rms=rms(div_members),
                commutation_residual=float(np.max(np.abs(field - div_members.mean(axis=0)))))


def frozen_field_probe(field, U0: float, u_rms: float, dt: float, dx: float = 1.0, n_samples: int | None = None,
                       n_realizations: int = 64, seed: int = 0) -> dict:
    """How well Taylor's frozen-turbulence hypothesis recovers a spatial pattern, as a function of u_rms/U0.

    Book: §12.4 (after (12.23)): replacing t by x/U0 is accurate when u_rms/U0 is small.
    Model (ours, "random sweeping"): a 1-D periodic pattern F(x) passes the probe at speed U0 + u′, where u′ is a
    large-eddy velocity, constant during the record, drawn from N(0, u_rms²) for each realization.  The probe records
    F((U0 + u′) t); the frozen hypothesis reads it as the pattern at x = U0 t.
    Parameters: field (n,) periodic pattern [unit], spacing dx [m]; U0 [m/s] > 0; u_rms [m/s] ≥ 0; dt sampling
    interval [s]; n_samples (default: one pass through the box); n_realizations; seed of ``default_rng`` (default 0).
    Returns dict(record (n_samples,) = what the probe of the first realization measures, reconstruction_error = rms of
    (u_probe − u_true)/rms(F) over all realizations [–], t, x_frozen = U0 t, u_true = F(U0 t), u_probe
    (n_realizations, n_samples), error_rms (the same number as reconstruction_error), ratio = u_rms/U0).  The error
    vanishes for u_rms = 0 and grows in proportion to the ratio while it is small.
    Assumptions: one-dimensional periodic pattern, linear interpolation between its nodes; the sweeping velocity is
    constant during a record (large eddies evolve slowly); a kinematic illustration, not a simulation.
    Validation: V7 error(0) = 0 to round-off; slope 1 ± 0.2 of log error against log(u_rms/U0).
    """
    F = _F(field).ravel()
    n = F.size
    Lbox = n * float(dx)
    ns = int(Lbox / (U0 * dt)) if n_samples is None else int(n_samples)
    t = np.arange(ns) * float(dt)
    xg = np.arange(n + 1) * float(dx)
    Fp = np.append(F, F[0])
    interp = lambda x: np.interp(np.mod(x, Lbox), xg, Fp)  # noqa: E731
    rng = np.random.default_rng(seed)
    up = u_rms * rng.standard_normal(int(n_realizations))
    u_probe = interp((U0 + up)[:, None] * t[None, :])
    u_true = interp(U0 * t)
    err = float(np.sqrt(np.mean((u_probe - u_true) ** 2)) / np.sqrt(np.mean(F ** 2)))
    return dict(record=u_probe[0], reconstruction_error=err, t=t, x_frozen=U0 * t, u_true=u_true, u_probe=u_probe,
                error_rms=err, ratio=float(u_rms / U0))


# ======================================================================================================================
# §12.5  mean stress, RANS residuals, fluxes
# ======================================================================================================================
def mean_stress_tensor(P, gradU, mu: float, rho0: float, uu):
    """Mean stress tensor of a turbulent flow: ``τ̄_ij = −P δ_ij + 2μ S̄_ij − ρ0 mean(u_i u_j)``.

    Book: §12.5, Eq. (12.30), with ``S̄_ij = ½(∂U_i/∂x_j + ∂U_j/∂x_i)``.
    Parameters: P mean pressure [Pa]; gradU (d, d) with ``gradU[i, j] = ∂U_i/∂x_j`` [1/s]; mu [Pa s]; rho0 [kg/m³];
    uu (d, d) velocity covariance mean(u_i u_j) [m²/s²].  Returns (d, d) array [Pa], symmetric.
    Validation: V1 simple shear: τ̄_12 = μ dU/dy − ρ0 mean(uv); isotropic uu adds to the pressure only.
    Assumptions: Boussinesq (constant rho0, mu); gradU is the gradient of the mean velocity.
    """
    G, R = _F(gradU), _F(uu)
    S = 0.5 * (G + G.T)
    return -float(P) * np.eye(G.shape[0]) + 2.0 * mu * S - rho0 * R  # Eq. (12.30)


def _d1(f: Callable, x: np.ndarray, i: int, h: float):
    e = np.zeros_like(x)
    e[i] = h
    return (_F(f(x + e)) - _F(f(x - e))) / (2.0 * h)


def rans_momentum_residual(fields: dict, x, *, nu: float, rho0: float = 1.0, g: float = 0.0, alpha: float = 0.0,
                           T0: float = 0.0, h: float = 1e-4, ht: float = 1e-4, t: float = 0.0) -> np.ndarray:
    """Residual of the Reynolds-averaged momentum equation at a point, for fields given as callables.

    Book: §12.5, Eq. (12.30):
    ``∂U_i/∂t + U_j ∂U_i/∂x_j = −g[1 − α(T̄ − T0)]δ_i3 + (1/ρ0) ∂/∂x_j(−P δ_ij + 2μ S̄_ij − ρ0 mean(u_i u_j))``.
    With g = 0 (the default) this is the constant-density form in which P is the deviation from hydrostatic.
    Parameters: fields dict of callables of (x, t) with x a (d,) point [m]: "U" → (d,) [m/s], "P" → scalar [Pa],
    "uu" → (d, d) [m²/s²], optional "T" → scalar [K]; x (d,) point, d = 2 or 3 (gravity acts along the last axis);
    nu [m²/s]; rho0 [kg/m³]; g [m/s²]; alpha [1/K]; T0 [K]; h, ht difference steps [m], [s]; t time [s].
    Returns (d,) residual left − right [m/s²].  Second-order central differences; the stress divergence is formed as
    the difference of the stress tensor itself (conservative form).
    Validation: V1 manufactured mean flow + stress field satisfying (12.30) ⇒ residual O(h²), order 2 ± 0.15.
    Assumptions: Boussinesq fluid with constant nu; fields smooth on the scale of the difference step.
    """
    x = _F(x)
    d = x.size
    U = lambda p: _F(fields["U"](p, t))  # noqa: E731

    def stress(p):
        grad = np.stack([_d1(U, p, j, h) for j in range(d)], axis=1)  # grad[i, j] = ∂U_i/∂x_j
        return (-float(fields["P"](p, t)) * np.eye(d) + rho0 * nu * (grad + grad.T) - rho0 * _F(fields["uu"](p, t)))

    gradU = np.stack([_d1(U, x, j, h) for j in range(d)], axis=1)
    dUdt = (_F(fields["U"](x, t + ht)) - _F(fields["U"](x, t - ht))) / (2.0 * ht)
    div_tau = sum(_d1(stress, x, j, h)[:, j] for j in range(d))
    body = np.zeros(d)
    if g != 0.0:
        Tbar = float(fields["T"](x, t)) if "T" in fields else T0
        body[-1] = -g * (1.0 - alpha * (Tbar - T0))
    return dUdt + gradU @ U(x) - body - div_tau / rho0  # Eq. (12.30)


def rans_2d_residual(U: Callable, V: Callable, P: Callable, uu: Callable, uv: Callable, vv: Callable, x, y, nu: float,
                     rho: float, h: float = 1e-4, return_continuity: bool = False) -> tuple:
    """Residuals of the steady two-dimensional constant-density mean-flow equations, for fields given as callables.

    Book: §12.8, Eqs. (12.58) ``∂U/∂x + ∂V/∂y = 0``, (12.59)
    ``U U_x + V U_y = −(1/ρ)P_x + ν(U_xx + U_yy) − ∂mean(u²)/∂x − ∂mean(uv)/∂y`` and (12.60)
    ``U V_x + V V_y = −(1/ρ)P_y + ν(V_xx + V_yy) − ∂mean(uv)/∂x − ∂mean(v²)/∂y``.
    Parameters
    ----------
    U, V : callables (x, y) → mean velocities [m/s], accepting arrays.   P : callable → mean pressure [Pa].
    uu, uv, vv : callables → mean(u²), mean(uv), mean(v²) [m²/s²].
    x, y : evaluation points [m].   nu : kinematic viscosity [m²/s].   rho : density [kg/m³].
    h : difference step [m].   return_continuity : also return the residual of (12.58).
    Returns
    -------
    (x_residual, y_residual) [m/s²] = left − right of (12.59) and (12.60) — or, with ``return_continuity=True``,
    (x_residual, y_residual, continuity_residual [1/s]).  Floats for scalar points.
    Assumptions: steady mean flow, constant density, two-dimensional mean; fields smooth on the scale h.
    Numerics: second-order central differences (five-point Laplacian); error O(h²).
    Validation: V1 manufactured fields; V3 order 2 ± 0.15 in h.
    """
    x, y = _F(x), _F(y)
    dx = lambda f: (_F(f(x + h, y)) - _F(f(x - h, y))) / (2.0 * h)  # noqa: E731
    dy = lambda f: (_F(f(x, y + h)) - _F(f(x, y - h))) / (2.0 * h)  # noqa: E731
    lap = lambda f: (_F(f(x + h, y)) + _F(f(x - h, y)) + _F(f(x, y + h)) + _F(f(x, y - h)) - 4.0 * _F(f(x, y))) / h ** 2  # noqa: E731
    Uv, Vv = _F(U(x, y)), _F(V(x, y))
    cont = dx(U) + dy(V)                                                                            # Eq. (12.58)
    rx = Uv * dx(U) + Vv * dy(U) - (-dx(P) / rho + nu * lap(U) - dx(uu) - dy(uv))                   # Eq. (12.59)
    ry = Uv * dx(V) + Vv * dy(V) - (-dy(P) / rho + nu * lap(V) - dx(uv) - dy(vv))                   # Eq. (12.60)
    if return_continuity:
        return _S(rx), _S(ry), _S(cont)
    return _S(rx), _S(ry)


def rans_eddy_viscosity_residual(U: Callable, P: Callable, nu_T: Callable, e: Callable, x, *, nu: float, rho: float,
                                 h: float = 1e-4, printed: bool = False) -> np.ndarray:
    """Residual of the RANS momentum equation closed with an eddy viscosity, for steady mean fields.

    Book: §12.10, Eq. (12.97) (constant density), obtained by putting (12.94) into (12.30):
    ``∂U_i/∂t + U_j ∂U_i/∂x_j = −(1/ρ) ∂P/∂x_i + ∂/∂x_j([ν + ν_T](∂U_i/∂x_j + ∂U_j/∂x_i) − (2/3) ē δ_ij)``.
    **The page prints the pressure gradient as ∂P/∂x_j; the free index of the equation is i (slip #6) — coded with i.**
    ``printed=True`` reproduces the page literally, so that a test can show it fails: with j repeated on the right-hand
    side the summation convention turns the term into ``−(1/ρ) Σ_j ∂P/∂x_j``, the same number in every component.
    Parameters
    ----------
    U : callable x → (d,) mean velocity [m/s].   P : callable x → mean pressure [Pa].
    nu_T : callable x → eddy viscosity [m²/s].   e : callable x → turbulent kinetic energy ē [m²/s²].
    x : (d,) point [m].   nu : [m²/s].   rho : [kg/m³].   h : difference step [m].
    printed : False (corrected index, the default) or True (the page's index — a planted wrong variant).
    Returns
    -------
    (d,) residual left − right [m/s²].  With ``printed=False`` it equals :func:`rans_momentum_residual` with the stress
    of :func:`eddy_viscosity_stress` — the test of the substitution (12.94) → (12.30).
    Assumptions: constant density; steady mean fields (the term ∂U_i/∂t of (12.97) is zero for the callables accepted
    here); the turbulent-viscosity hypothesis (12.94).
    Validation: V1 manufactured field, order 2; identity with the (12.30) residual to 1e-8; ``printed=True`` leaves an
    O(1) residual on the same manufactured field whenever the pressure gradient is not the same in every direction.
    """
    x = _F(x)
    d = x.size
    Uf = lambda p: _F(U(p))  # noqa: E731

    def flux(p):
        grad = np.stack([_d1(Uf, p, j, h) for j in range(d)], axis=1)
        return (nu + float(nu_T(p))) * (grad + grad.T) - (2.0 / 3.0) * float(e(p)) * np.eye(d)

    gradU = np.stack([_d1(Uf, x, j, h) for j in range(d)], axis=1)
    gradP = np.array([float(_d1(lambda p: _F(P(p)), x, i, h)) for i in range(d)])
    div = sum(_d1(flux, x, j, h)[:, j] for j in range(d))
    # DEVIATION: ∂P/∂x_i, not the printed ∂P/∂x_j (slip #6) — j is a summed index on the right-hand side, the free
    # index of the equation is i.  printed=True keeps the page's form available so that a test can show it fails.
    if printed:
        gradP = np.full(d, float(np.sum(gradP)))  # the page's ∂P/∂x_j, j summed: one number for every component
    return gradU @ Uf(x) + gradP / rho - div  # Eq. (12.97), corrected index (unless printed=True)


def mean_heat_flux(gradT, uT, k_th: float, rho0: float, cp: float):
    """Mean heat flux of a turbulent flow: ``Q_j = −k ∂T̄/∂x_j + ρ0 C_p mean(u_j T′)``.

    Book: §12.5, Eq. (12.32) (k = ρ0 C_p κ is the thermal conductivity — ``k_th`` in fluidpy, because k is also a
    wavenumber and the "k" of k–ε in this chapter; the second term is the turbulent heat flux).
    Parameters
    ----------
    gradT : (d,) mean temperature gradient [K/m].   uT : (d,) velocity–temperature correlation mean(u_j T′) [K m/s].
    k_th : thermal conductivity [W/(m K)].   rho0 : reference density [kg/m³].   cp : specific heat [J/(kg K)].
    Returns
    -------
    Q : (d,) mean heat flux [W/m²] (a float for scalar input).
    Assumptions: Boussinesq fluid with constant k_th, ρ0 and C_p.
    Validation: V1 a hand value; units; the turbulent part alone equals :func:`turbulent_heat_flux`.
    """
    return _S(-k_th * _F(gradT) + rho0 * cp * _F(uT))  # Eq. (12.32)


def mean_scalar_flux(gradY, uY, kappa_m: float):
    """Mean flux of a passive scalar per unit density: ``−κ_m ∂Ȳ/∂x_j + mean(u_j Y′)``.

    Book: §12.5, Eq. (12.34) ``∂Ȳ/∂t + U_j ∂Ȳ/∂x_j = ∂/∂x_j(κ_m ∂Ȳ/∂x_j − mean(u_j Y′))`` — the right side is minus the
    divergence of this flux.  Parameters: gradY (d,) [1/m]; uY (d,) [m/s]; kappa_m molecular diffusivity [m²/s].
    Returns (d,) [m/s] (multiply by the mixture density for kg/(m² s)).  Assumes constant mixture density.
    Validation: V1 a hand value; zero gradient leaves the turbulent flux alone.
    """
    return _S(-kappa_m * _F(gradY) + _F(uY))  # Eq. (12.34)


def turbulent_heat_flux(w_rms, T_rms, r_wT, rho: float, cp: float):
    """Turbulent heat flux from the two rms values and their correlation coefficient: ``ρ C_p r w_rms T_rms``.

    Book: §12.5 (the term ρ0 C_p mean(u_j T′) of (12.32)) with the coefficient of (12.14):
    mean(wT′) = r_wT w_rms T_rms.  Parameters: w_rms [m/s]; T_rms [K]; r_wT in [−1, 1]; rho [kg/m³]; cp [J/(kg K)].
    Returns H [W/m²], positive upward when r_wT > 0 (daytime heating).  Scalar-callable.
    Assumptions: constant rho and C_p over the fluctuation; r_wT is the correlation coefficient (12.14).
    Validation: V1 a hand value (0.5 m/s, 0.2 K, r = 0.5 in air gives 60 W/m2); sign follows r_wT.
    """
    return _S(rho * cp * _F(r_wT) * _F(w_rms) * _F(T_rms))


def mixture_density(v, rho_s: float, rho: float):
    """Density of a binary mixture from the volume fraction v of the scalar fluid: ``ρ_m = v ρ_s + (1 − v) ρ``.

    Book: §12.5 (before (12.33)).  v in [0, 1]; rho_s, rho [kg/m³].  Returns ρ_m [kg/m³].
    Assumptions: ideal mixing: volumes add (no volume change on mixing).
    Validation: V1 limits v = 0 and v = 1; linear in v.
    """
    return _S(_F(v) * rho_s + (1.0 - _F(v)) * rho)


def mass_fraction_from_volume_fraction(v, rho_s: float, rho: float):
    """Mass fraction of the scalar fluid from its volume fraction: ``Y = v ρ_s/ρ_m``.  Book: §12.5 (before (12.33)).

    Returns Y [-] (kg of scalar fluid per kg of mixture); a float for float input.
    Assumptions: ideal mixing: volumes add.
    Validation: V1 Y = 0 and 1 at v = 0 and 1; Y = v when the densities are equal.
    """
    return _S(_F(v) * rho_s / _F(mixture_density(v, rho_s, rho)))


# ======================================================================================================================
# §12.5  why mean(uv) < 0 in a positive shear; the closure count
# ======================================================================================================================
def parcel_uv_expected(dUdy: float, l_rms: float, v_rms: float, correlation: float = 1.0) -> float:
    """Expected Reynolds shear correlation of the displaced-parcel model: ``mean(uv) = −c v_rms l_rms dU/dy``.

    Book: §12.5, Fig. 12.6 (a parcel displaced upward by ℓ keeps the mean speed of its origin, so u ≈ −ℓ dU/dy while
    v > 0) and §12.10 (the mixing-length estimate).  c is the correlation coefficient between v and ℓ.
    Parameters: dUdy [1/s]; l_rms [m]; v_rms [m/s]; correlation c in [0, 1].  Returns mean(uv) [m²/s²]: negative for
    positive shear, zero when v and ℓ are uncorrelated.  Scalar-callable (the closed form for explainer parity).
    Assumptions: first-order Taylor expansion of U over a displacement; v and the displacement jointly Gaussian with correlation c.
    Validation: V1 the mean of displaced_parcel_uv agrees within 5 standard errors; sign flips with dU/dy.
    """
    return float(-correlation * v_rms * l_rms * dUdy)


def displaced_parcel_uv(dUdy: float, l_rms: float, v_rms: float, n: int = 100000, seed: int = 0,
                        correlation: float = 1.0, u_extra_rms: float = 0.0) -> dict:
    """Monte-Carlo version of the displaced-parcel argument for the sign of the Reynolds shear stress.

    Book: §12.5, Fig. 12.6 and the momentum-flux reading that follows: each parcel arrives at the reference level with
    wall-normal velocity v after a displacement ℓ and carries the velocity deficit u = −ℓ dU/dy (first-order Taylor
    expansion of U).  Up-moving parcels (v > 0, ℓ > 0) are slow, down-moving ones fast: every exchange gives uv < 0.
    Parameters: dUdy mean shear [1/s] (either sign); l_rms displacement scale [m]; v_rms [m/s]; n parcels;
    seed of ``default_rng`` (default 0); correlation between v and ℓ (1 = a parcel keeps all of its momentum,
    0 = it forgets where it came from); u_extra_rms [m/s] optional uncorrelated part of u (fattens the scatter cloud).
    Returns dict(uv = mean(uv) [m²/s²], estimate = −mean(vℓ) dU/dy, expected (closed form,
    :func:`parcel_uv_expected`), r_uv [–], stderr = standard error of uv (std of the products/√n), nu_T = −uv/(dU/dy)
    [m²/s], samples = (u, v, ℓ) three (n,) arrays — also available as "u", "v", "l").  ``uv == estimate`` to round-off
    when u_extra_rms = 0.
    Assumptions: first-order Taylor expansion of U over a displacement; v and ℓ jointly Gaussian; a kinematic
    illustration of the sign argument, not a model of a real flow.
    Validation: V1 uv within 5 standard errors of :func:`parcel_uv_expected`; sign flips with dU/dy; zero for c = 0.
    """
    v, ell = TS.correlated_pair(int(n), float(correlation), v_rms, l_rms, seed=seed)
    u = -ell * dUdy
    if u_extra_rms:
        u = u + u_extra_rms * np.random.default_rng(seed + 1).standard_normal(int(n))
    prod = u * v
    uv = float(np.mean(prod))
    return dict(uv=uv, estimate=float(-np.mean(v * ell) * dUdy),
                expected=parcel_uv_expected(dUdy, l_rms, v_rms, correlation),
                r_uv=float(TS.correlation_coefficient(u - u.mean(), v - v.mean())) if np.any(u != 0) else float("nan"),
                stderr=float(np.std(prod, ddof=1) / np.sqrt(n)),
                nu_T=float(-uv / dUdy) if dUdy != 0 else float("nan"), samples=(u, v, ell), u=u, v=v, l=ell)


def closure_count(level: int = 3) -> list:
    """The closure problem as a count of equations against unknown moments (three space dimensions, isothermal).

    Book: §12.5 (after (12.35)): the equations for the first moments contain second moments, those for the second
    moments contain third moments, and so on.
    Returns a list of dicts, one per closure level n = 1 … ``level``, with cumulative counts:
    ``equations`` — continuity + 3 momentum (n = 1), + 6 for mean(u_i u_j) (n = 2), + 10 for mean(u_i u_j u_k) (n = 3), …;
    ``unknowns`` — U_i (3), P (1) and every velocity moment up to order n + 1 (a symmetric tensor of order m has
    (m + 1)(m + 2)/2 components: 6, 10, 15, …); ``extra`` names the further unknown correlations that appear at that
    level (pressure–velocity and dissipation tensors), counted in ``unknowns_with_extra`` for n = 2.
    Level 1: 4 equations, 10 unknowns.  Level 2: 10 equations, 20 velocity unknowns + 12 others.  Never closed.
    Validation: V1 the counts are re-derived in the test from the component formula.
    Assumptions: three space dimensions, isothermal, incompressible; symmetric moment tensors.
    """
    comps = lambda m: (m + 1) * (m + 2) // 2  # noqa: E731  symmetric tensor of order m in 3-D
    rows = []
    eq, unk = 4, 4
    for n in range(1, int(level) + 1):
        if n > 1:
            eq += comps(n)
        unk += comps(n + 1)
        extra = {1: [], 2: ["pressure–velocity-gradient correlation (6)", "dissipation tensor 2ν mean(∂u_i/∂x_k ∂u_j/∂x_k) (6)"]}.get(
            n, ["pressure and dissipation correlations of order %d (not counted)" % n])
        rows.append(dict(level=n, equations=eq, unknowns=unk, unknowns_with_extra=unk + (12 if n >= 2 else 0),
                         new_moment="mean of %d velocity fluctuations (%d components)" % (n + 1, comps(n + 1)),
                         extra=extra, closed=False))
    return rows


# ======================================================================================================================
# §12.6  homogeneous isotropic turbulence
# ======================================================================================================================
def isotropy_report(u, v, w=None, dx: float = 1.0) -> dict:
    """How isotropic is a sampled periodic velocity field?  Normal stresses, shear correlations, gradient moments.

    Book: §12.6, Eqs. (12.36) (equal normal stresses, equal same-direction gradient moments) and (12.37) (equal
    cross-direction gradient moments), with n = 2; isotropy also requires mean(uv) = 0 (Fig. 12.8).
    Parameters: u, v (and w) periodic arrays [m/s], x on the last axis; dx [m].  **Volume averages.**
    Returns dict(variances, normal_stress_ratios (each over their mean), shear_coefficients (r_12, …), longitudinal
    (mean (∂u_i/∂x_i)² for each i), transverse (mean (∂u_i/∂x_j)², i ≠ j), transverse_over_longitudinal (2 in 3-D
    isotropic incompressible turbulence, by (12.43); 3 for a 2-D solenoidal field)).  Spectral derivatives.
    Assumptions: periodic box; volume average of one realization in place of the ensemble average.
    Validation: V1 a 3-D synthetic solenoidal field gives ratios within sampling error of the isotropic values; a sheared field does not.
    """
    comps = [_F(u), _F(v)] + ([] if w is None else [_F(w)])
    nd = comps[0].ndim
    fl = [c - c.mean() for c in comps]
    var = np.array([np.mean(c ** 2) for c in fl])
    shear = {f"r_{i + 1}{j + 1}": float(np.mean(fl[i] * fl[j]) / np.sqrt(var[i] * var[j]))
             for i in range(len(fl)) for j in range(i + 1, len(fl))}
    lon, tra = [], []
    for i, c in enumerate(fl):
        for j in range(len(fl)):
            m2 = float(np.mean(_spectral_derivative(c, float(dx), nd - 1 - j) ** 2))
            (lon if i == j else tra).append(m2)
    return dict(variances=var, normal_stress_ratios=var / var.mean(), shear_coefficients=shear,
                longitudinal=np.array(lon), transverse=np.array(tra),
                transverse_over_longitudinal=float(np.mean(tra) / np.mean(lon)))


def _derivative_of(f, r, fprime=None, h: float = 1e-5):
    if fprime is not None:
        return _F(fprime(r))
    if hasattr(f, "fprime"):
        return _F(f.fprime(r))
    r = _F(r)
    return (_F(f(r + h)) - _F(f(np.abs(r - h)))) / (2.0 * h)  # f is even, so f(−x) = f(x)


def transverse_from_longitudinal(r, f, dim: int = 3, fprime=None):
    """Transverse correlation g(r) of incompressible isotropic turbulence from the longitudinal one f(r).

    Book: §12.6, the consequence of (12.40)–(12.41) (Exercise 12.18): ``g = f + (r/2) df/dr`` — **three dimensions**.
    ``dim=2`` gives the two-dimensional analogue ``g = f + r df/dr = d(r f)/dr`` (ours), which is what a 2-D solenoidal
    field obeys; do not mix the two.
    Parameters: r ≥ 0 [m]; f either a callable r → f (derivative analytic if ``fprime`` or ``f.fprime`` is given, else a
    central difference) or samples on the grid r (second-order differences, ``np.gradient``).  Returns g [–].
    Validation: V2 Gaussian f = exp(−r²/L²) ⇒ g = (1 − r²/L²) f (3-D); V1 a 3-D synthetic field.
    Assumptions: homogeneous, isotropic, incompressible turbulence; f smooth.
    """
    r_ = _F(r)
    if callable(f):
        fv, fp = _F(f(r_)), _derivative_of(f, r_, fprime)
    else:
        fv = _F(f)
        fp = np.gradient(fv, r_, edge_order=2)
    c = {3: 0.5, 2: 1.0}[int(dim)]
    return _S(fv + c * r_ * fp)  # g = f + (r/2) f'   (from Eq. (12.41))


def isotropic_correlation_tensor(rvec, f: Callable, g: Callable | None = None, u2: float = 1.0,
                                 incompressible: bool = True, fprime: Callable | None = None) -> np.ndarray:
    """Two-point velocity correlation tensor R_ij(r) of homogeneous isotropic turbulence.

    Book: §12.6, Eq. (12.40) ``R_ij = F(r) r_i r_j + G(r) δ_ij`` with ``F = u2 (f − g)/r²``, ``G = u2 g`` (general), and
    Eq. (12.41) ``R_ij = u2 {f δ_ij + (r/2)(df/dr)(δ_ij − r_i r_j/r²)}`` (incompressible).
    Parameters: rvec (3,) separation [m]; f, g callables r → coefficient; u2 one-component variance [m²/s²];
    incompressible True → (12.41) (g ignored; df/dr from ``fprime``/``f.fprime`` or a central difference),
    False → (12.40) (g required).  Returns (3, 3) array [m²/s²].
    Checks: along r the component is u2·f, across it u2·g; R_ii(0) = 3 u2 = 2ē.
    Validation: V2 divergence ∂R_ij/∂r_j = 0 (:func:`isotropic_tensor_divergence_sympy`); V1 components.
    Assumptions: homogeneous isotropic turbulence (and incompressible for (12.41)).
    """
    rv = _F(rvec)
    r = float(np.linalg.norm(rv))
    eye = np.eye(3)
    if r == 0.0:
        return u2 * float(f(0.0)) * eye
    nn = np.outer(rv, rv) / r ** 2
    if incompressible:
        fp = float(_derivative_of(f, r, fprime))
        return u2 * (float(f(r)) * eye + 0.5 * r * fp * (eye - nn))  # Eq. (12.41)
    if g is None:
        raise ValueError("the general form (12.40) needs g")
    return u2 * ((float(f(r)) - float(g(r))) * nn + float(g(r)) * eye)  # Eq. (12.40)


def isotropic_scales(r, f) -> dict:
    """Integral scales and Taylor microscales of an incompressible isotropic field from its longitudinal correlation.

    Book: §12.6, Eq. (12.39) ``Λ_f = ∫f dr``, ``Λ_g = ∫g dr``, ``λ_f² = −2/f″(0)``, ``λ_g² = −2/g″(0)`` and the results
    ``Λ_g = Λ_f/2``, ``λ_g = λ_f/√2`` stated after (12.41) (three dimensions).
    Parameters: r (n,) grid from 0 to where f has decayed [m]; f samples (n,) or a callable.
    Returns dict(Lambda_f, Lambda_g, lambda_f, lambda_g [m], Lambda_ratio, lambda_ratio, g (n,)).  Integrals over the
    whole range given ("all"); microscales by the even fit of ``TS.taylor_microscale``.
    Validation: V1 Gaussian f: Λ_g/Λ_f = ½ and λ_g/λ_f = 1/√2 to 1e-6 on a fine grid.
    Assumptions: homogeneous isotropic incompressible turbulence in three dimensions; f has decayed at the end of the grid.
    """
    r_ = _F(r)
    fv = _F(f(r_)) if callable(f) else _F(f)
    g = _F(transverse_from_longitudinal(r_, f if callable(f) else fv, dim=3))
    Lf, Lg = TS.integral_scale(r_, fv, "all"), TS.integral_scale(r_, g, "all")
    lf, lg = TS.taylor_microscale(r_, fv), TS.taylor_microscale(r_, g)
    return dict(Lambda_f=Lf, Lambda_g=Lg, lambda_f=lf, lambda_g=lg, Lambda_ratio=Lg / Lf, lambda_ratio=lg / lf, g=g)


def velocity_gradient_samples(u, v, w, dx: float) -> np.ndarray:
    """All nine velocity gradients of a periodic three-dimensional field, one sample per grid node.

    Book: §12.6, Eq. (12.42) needs ∂u_i/∂x_j of the fluctuating field.  Parameters: u, v, w periodic arrays [m/s] with
    x on the last axis, y on the second-to-last, z on the first; dx grid spacing [m].
    Returns (N, 3, 3) array with ``[n, i, j] = ∂u_i/∂x_j`` [1/s] (spectral derivatives), ready for
    :func:`dissipation_rate`.  Volume samples of one realization (a spatial average stands in for the ensemble).
    Assumptions: periodic box, uniform grid, the same spacing in the three directions.
    Validation: V1 exact derivatives of a single Fourier mode; the trace vanishes for a solenoidal field.
    """
    comps = (_F(u), _F(v), _F(w))
    out = np.empty((comps[0].size, 3, 3))
    for i, c in enumerate(comps):
        for j in range(3):
            out[:, i, j] = _spectral_derivative(c, float(dx), 2 - j).ravel()
    return out


def dissipation_rate(grad_u_samples, nu: float) -> float:
    """Mean dissipation rate of turbulent kinetic energy from samples of the fluctuating velocity gradient.

    Book: §12.6, Eq. (12.42) ``ε̄ = (ν/2) mean((∂u_i/∂x_j + ∂u_j/∂x_i)²)`` (sum over i and j) = 2ν mean(S′_ij S′_ij).
    Parameters: grad_u_samples (N, d, d) with ``[n, i, j] = ∂u_i/∂x_j`` [1/s]; nu [m²/s].  Returns ε̄ [m²/s³] ≥ 0.
    Validation: V1 pure strain and pure rotation (rotation dissipates nothing); V4 equals 2ν∫K²E dK for a synthetic
    3-D field.
    Assumptions: Newtonian fluid, incompressible fluctuations; the over-bar is estimated by the average over the samples supplied (an estimate with sampling error ~ N^{-1/2}).
    """
    G = _F(grad_u_samples)
    S2 = (G + np.swapaxes(G, -1, -2)) ** 2
    return float(0.5 * nu * np.mean(np.sum(S2, axis=(-1, -2))))  # Eq. (12.42)


def dissipation_isotropic(nu: float, u2: float | None = None, lambda_f: float | None = None,
                          lambda_g: float | None = None, dudx_sq: float | None = None):
    """Dissipation rate of isotropic turbulence from one measurable quantity.

    Book: §12.6, Eq. (12.43): ``ε̄ = 30 ν u2/λ_f² = 15 ν u2/λ_g²`` and, from its first form with the gradient moments in
    the ratio 2 : 4 : −1 (Exercise 12.19), ``ε̄ = 15 ν mean((∂u_1/∂x_1)²)``.
    Parameters: nu [m²/s]; **exactly one** of lambda_f [m] (longitudinal Taylor microscale), lambda_g [m] (transverse),
    dudx_sq = mean((∂u_1/∂x_1)²) [1/s²]; u2 one-component variance [m²/s²] (needed with a microscale).
    Returns ε̄ [m²/s³].  Scalar-callable.
    Validation: V1 the three forms agree for a Gaussian f (λ_g = λ_f/√2, mean((∂u_1/∂x_1)²) = 2u2/λ_f²); the
    λ_f-for-λ_g mutant is off by a factor 2.
    Assumptions: locally isotropic small scales (high Reynolds number).
    """
    given = [a is not None for a in (lambda_f, lambda_g, dudx_sq)]
    if sum(given) != 1:
        raise ValueError("give exactly one of lambda_f, lambda_g, dudx_sq")
    if dudx_sq is not None:
        return _S(15.0 * nu * _F(dudx_sq))                     # Eq. (12.43), one-gradient form
    if u2 is None:
        raise ValueError("u2 (one-component variance) is needed with a Taylor microscale")
    if lambda_f is not None:
        return _S(30.0 * nu * _F(u2) / _F(lambda_f) ** 2)      # Eq. (12.43)
    return _S(15.0 * nu * _F(u2) / _F(lambda_g) ** 2)          # Eq. (12.43)


def taylor_reynolds_number(u2, lam, nu):
    """Taylor-scale Reynolds number ``R_λ = λ sqrt(u2)/ν``.

    Book: §12.6, Eq. (12.44) (λ may be λ_f or λ_g — a factor √2: say which).  u2 one-component variance [m²/s²];
    lam [m]; nu [m²/s].  Returns R_λ [–].
    Assumptions: u2 is the variance of one component; say which microscale is used.
    Validation: V1 a hand value (471.4 for u2 = 1, lambda_g = 0.01/sqrt(2), nu = 1.5e-5).
    """
    return _S(_F(lam) * np.sqrt(_F(u2)) / _F(nu))  # Eq. (12.44)


# ======================================================================================================================
# §12.7  energy budgets
# ======================================================================================================================
def shear_production(uu, gradU):
    """Shear production of turbulent kinetic energy ``P = −mean(u_i u_j) ∂U_i/∂x_j`` (double contraction).

    Book: §12.7, the fifth term of (12.46) and of (12.47): the same quantity leaves the mean flow and enters the
    turbulence.  Parameters: uu (..., d, d) covariance [m²/s²]; gradU (..., d, d) with [i, j] = ∂U_i/∂x_j [1/s].
    Returns P [m²/s³] (positive = gain of turbulent energy).  For U(y): P = −mean(uv) dU/dy.
    Validation: V1 simple shear; isotropic stress with a divergence-free mean flow gives 0.
    Assumptions: none beyond the definition (a double contraction of two given tensors).
    """
    return _S(-np.einsum("...ij,...ij->...", _F(uu), _F(gradU)))  # Eqs. (12.46)–(12.47)


def reynolds_stress_production(uu, gradU) -> np.ndarray:
    """Production tensor of the Reynolds-stress budget: ``P_ij = −mean(u_i u_k) ∂U_j/∂x_k − mean(u_j u_k) ∂U_i/∂x_k``.

    Book: §12.5, Eq. (12.35) (first two terms on the right).  Returns (d, d) [m²/s³]; ½ trace = :func:`shear_production`.
    In simple shear U(y) all of it goes into the streamwise component: P_11 = −2 mean(uv) dU/dy, P_22 = P_33 = 0.
    Assumptions: none beyond the definition; gradU[i, j] = dU_i/dx_j.
    Validation: V1 half the trace equals shear_production; simple shear puts all of it in the 11 component.
    """
    R, G = _F(uu), _F(gradU)
    return -(R @ G.T) - (R @ G.T).T  # Eq. (12.35)


def buoyant_production(wT, alpha: float, g: float = G0):
    """Buoyant production (or destruction) of turbulent kinetic energy ``g α mean(w T′)``.

    Book: §12.7, the last term of (12.47).  wT [K m/s] upward positive (**potential** temperature fluctuation); alpha
    thermal expansion coefficient [1/K]; g [m/s²].  Returns [m²/s³]: positive for an upward heat flux (unstable),
    negative for a downward one (stable).  Scalar-callable.
    Assumptions: Boussinesq; wT is the flux of potential temperature.
    Validation: V1 a hand value; sign follows the heat flux.
    """
    return _S(g * alpha * _F(wT))  # Eq. (12.47)


def mean_to_turbulent_dissipation_ratio(Re, urms_over_U: float = 1.0):
    """Order of magnitude of (direct viscous dissipation of the mean flow)/(shear production): ``1/Re`` when u_rms ∼ U.

    Book: §12.7 (after (12.46)): ``2ν S̄_ij S̄_ij / [mean(u_i u_j) ∂U_i/∂x_j] ∼ ν(U/L)²/(u_rms² U/L) ∼ ν/(UL) = 1/Re``.
    Parameters: Re = UL/ν; urms_over_U (the book takes it of order one).  Returns the ratio [–].
    Assumptions: an order-of-magnitude scaling law (the constant is an argument).
    Validation: V1 arithmetic; the model channel shows the ratio falling with distance from the wall.
    """
    return _S(1.0 / (_F(Re) * urms_over_U ** 2))


def _integral(f, y):
    return float(np.trapezoid(f, y))


def mean_energy_budget(y, U, uv, nu: float, dPdx: float | None = None, rho: float = 1.0) -> dict:
    """Terms of the mean-flow kinetic-energy budget for a steady unidirectional mean flow U(y).

    Book: §12.7, Eq. (12.46).  For U = U(y) e_x with statistics independent of x it reduces to (per unit mass)
    ``0 = −U (dP/dx)/ρ  +  d/dy(ν U dU/dy)  +  d/dy(−mean(uv) U)  −  ν (dU/dy)²  +  mean(uv) dU/dy``:
    pressure work, viscous transport, turbulent transport, direct viscous dissipation (2ν S̄_ij S̄_ij), loss to turbulence.
    Parameters: y (n,) [m]; U (n,) [m/s]; uv (n,) mean(uv) [m²/s²]; nu [m²/s]; dPdx [Pa/m] (None: taken from the
    momentum balance, (1/ρ)dP/dx = [ν U′ − uv] end-to-end difference over the height); rho [kg/m³].
    Returns dict(pressure_work, viscous_transport, turbulent_transport, viscous_dissipation, loss_to_turbulence,
    residual (their sum) — arrays [m²/s³] — and ``integrals`` with the trapezoid integral of each over y [m³/s³]).
    Second-order differences (``np.gradient``, ``edge_order=2``).
    Reading: the transport terms integrate to zero between no-slip walls; what the pressure gradient puts in leaves
    through direct dissipation and — far more, by a factor ∼ Re — through the loss to turbulence.
    Validation: V4 laminar limit (uv = 0): work in = dissipation for plane Poiseuille and Couette; with the model
    channel the integral identity closes; the exchange term is minus the production of :func:`tke_budget`.
    Assumptions: steady unidirectional mean flow U(y), statistics independent of x and z; constant density (no buoyancy term).
    """
    y, U, uv = _F(y), _F(U), _F(uv)
    dU = np.gradient(U, y, edge_order=2)
    flux = nu * dU - uv
    if dPdx is None:
        dPdx_rho = (flux[-1] - flux[0]) / (y[-1] - y[0])
    else:
        dPdx_rho = float(dPdx) / rho
    terms = dict(pressure_work=-U * dPdx_rho,
                 viscous_transport=np.gradient(nu * U * dU, y, edge_order=2),
                 turbulent_transport=np.gradient(-uv * U, y, edge_order=2),
                 viscous_dissipation=-nu * dU ** 2,
                 loss_to_turbulence=uv * dU)  # Eq. (12.46)
    terms["residual"] = sum(terms.values())
    terms["integrals"] = {k: _integral(v, y) for k, v in terms.items()}
    terms["dPdx_over_rho"] = float(dPdx_rho)
    return terms


def tke_budget(y, U, uv, eps, transport=None, buoyancy=None) -> dict:
    """Terms of the turbulent-kinetic-energy budget for a steady mean shear flow U(y).

    Book: §12.7, Eq. (12.47) reduced to U(y): ``0 = transport − ε̄ − mean(uv) dU/dy + g α mean(w T′)``.
    (The label under the left side of the printed equation names the mean-flow energy; it is the turbulent ē — slip #2.)
    Parameters: y (n,) [m]; U (n,) [m/s]; uv (n,) [m²/s²]; eps (n,) dissipation rate [m²/s³] ≥ 0; transport (n,) the
    divergence term [m²/s³] (None → returned as the residual that closes the budget, i.e. **inferred, not measured**);
    buoyancy (n,) g α mean(wT′) [m²/s³] (None → 0).
    Returns dict(production, dissipation (= −eps), buoyancy, transport, residual, transport_inferred (bool), integrals).
    Validation: V1 production equals minus the "loss_to_turbulence" of :func:`mean_energy_budget` point by point.
    Assumptions: steady mean shear flow U(y), statistics independent of x and z; a transport term passed as None is inferred, not measured.
    """
    y, U, uv, eps = _F(y), _F(U), _F(uv), _F(eps)
    prod = -uv * np.gradient(U, y, edge_order=2)  # Eq. (12.47), fifth term
    buoy = np.zeros_like(prod) if buoyancy is None else _F(buoyancy)
    inferred = transport is None
    tr = -(prod - eps + buoy) if inferred else _F(transport)
    out = dict(production=prod, dissipation=-eps, buoyancy=buoy, transport=tr, residual=prod - eps + buoy + tr)
    out["integrals"] = {k: _integral(v, y) for k, v in out.items()}
    out["transport_inferred"] = inferred
    return out


def mixing_potential_energy_change(z, T_initial, alpha: float, rho0: float, g: float = G0) -> dict:
    """Change of potential energy when a column is mixed to a uniform temperature.

    Book: §12.7, Fig. 12.10 and the reading of the buoyant term of (12.47): in an unstable layer the heat flux is
    upward, the mean potential energy falls and feeds the turbulence; mixing a stable layer costs energy.
    Parameters: z (n,) heights [m]; T_initial (n,) **potential** temperature [K]; alpha [1/K]; rho0 [kg/m³]; g.
    Model: ρ = ρ0[1 − α(T − T_ref)]; the final state is the depth-mean temperature (heat conserved, Boussinesq).
    Returns dict(dPE [J/m²] = PE_final − PE_initial = −ρ0 g α ∫(T_f − T_i) z dz, T_final [K], heat_flux_sign (+1 upward,
    −1 downward, 0), stable_initially (bool: dT/dz ≥ 0 everywhere)).  dPE < 0 for an unstable column.
    Validation: V1 linear profile T = T0 + Gz over depth H: dPE = ρ0 g α G H³/12.
    Assumptions: Boussinesq, linear equation of state; complete mixing at constant heat content.
    """
    z, T = _F(z), _F(T_initial)
    H = z[-1] - z[0]
    Tf = float(np.trapezoid(T, z) / H)
    dPE = float(-rho0 * g * alpha * np.trapezoid((Tf - T) * z, z))
    return dict(dPE=dPE, T_final=Tf, heat_flux_sign=int(-np.sign(dPE)) if dPE != 0 else 0,
                stable_initially=bool(np.all(np.diff(T) >= 0)))


# ======================================================================================================================
# §12.2, §12.7  scales of the cascade
# ======================================================================================================================
def dissipation_outer_scaling(dU, L, c: float = 1.0):
    """Dissipation rate fixed by the large eddies: ``ε̄ ∼ c (ΔU)³/L``.

    Book: §12.7, Eqs. (12.48) ``Ẇ ∼ (ΔU)²[ΔU/L] = (ΔU)³/L`` (energy input to the large eddies) and (12.49) ``Ẇ = ε̄``
    (stationary turbulence).  **No viscosity appears.**
    Parameters: dU outer velocity difference [m/s]; L outer length [m]; c order-one constant (a scaling law: c = 1 by
    default, an argument so that its effect can be shown).  Returns ε̄ [m²/s³] = [W/kg].  Scalar-callable.
    Assumptions: stationary high-Reynolds-number turbulence fed by the large eddies; an order-of-magnitude scaling law (the constant is an argument).
    Validation: V2 units m2/s3; independent of viscosity by construction.
    """
    return _S(c * _F(dU) ** 3 / _F(L))  # Eqs. (12.48)–(12.49)


def kolmogorov_scales(nu, eps) -> tuple:
    """Kolmogorov length, velocity and time scales of the dissipating eddies.

    Book: §12.7, Eq. (12.50) ``η = (ν³/ε̄)^{1/4}``, ``u_K = (ν ε̄)^{1/4}``; the time ``τ_η = η/u_K = (ν/ε̄)^{1/2}`` is ours.
    Their Reynolds number ``η u_K/ν`` is exactly one.
    Parameters
    ----------
    nu : kinematic viscosity [m²/s] > 0.   eps : mean dissipation rate ε̄ [m²/s³] > 0.
    Returns
    -------
    (eta [m], u_K [m/s], tau_eta [s]) — a 3-tuple of floats for scalar input (arrays for array input), so that
    ``eta, uK, tau_eta = kolmogorov_scales(nu, eps)`` and ``kolmogorov_scales(nu, eps)[0]`` both work.
    Assumptions: the smallest eddies depend on ν and ε̄ only (Kolmogorov's first similarity hypothesis; high Reynolds
    number, locally isotropic small scales).
    Validation: V2 exponents (¾, −¼), (¼, ¼), (½, −½) from ``core.dimensional.solve_exponents``; V1 η u_K/ν = 1.
    """
    nu_, e_ = _F(nu), _F(eps)
    eta = (nu_ ** 3 / e_) ** 0.25   # Eq. (12.50)
    uK = (nu_ * e_) ** 0.25         # Eq. (12.50)
    return _S(eta), _S(uK), _S(np.sqrt(nu_ / e_))


def kolmogorov_exponents() -> dict:
    """The exponents of (12.50) by exact dimensional analysis: η = ν^a ε^b etc., as ``Fraction`` pairs (a, b).

    Book: §12.7, Eq. (12.50).  Uses ``core.dimensional.solve_exponents`` (ch01).  Returns dict(eta, u_K, tau_eta).
    Assumptions: only nu and eps matter at the smallest scales.
    Validation: V2 exact fractions (3/4, -1/4), (1/4, 1/4), (1/2, -1/2).
    """
    from .core.dimensional import solve_exponents

    out = {}
    for name, unit in (("eta", "m"), ("u_K", "m/s"), ("tau_eta", "s")):
        g = solve_exponents(name, ("nu", "eps"), {name: unit, "nu": "m**2/s", "eps": "m**2/s**3"})
        out[name] = (-g["nu"], -g["eps"])  # name · nu^p eps^q dimensionless ⇒ name = nu^{−p} eps^{−q}
    return out


def scale_separation(Re_L, c_eps: float = 1.0, u2_over_dU2: float = 1.0, microscale: str = "f") -> dict:
    """Ratios of the small scales to the outer scales as functions of the outer Reynolds number Re_L = ΔU L/ν.

    Book: §12.7, Eqs. (12.51) ``η/L ∼ Re_L^{−3/4}`` and (12.52) ``λ_T/L ∝ Re_L^{−1/2}``, obtained by putting
    ε̄ = c_eps (ΔU)³/L (12.49) into (12.50) and into (12.43).
    Parameters: Re_L [–]; c_eps the constant of (12.49); u2_over_dU2 = (one-component variance)/(ΔU)²; microscale
    "f" (factor 30 of (12.43)) or "g" (factor 15).
    Returns dict(eta_over_L = c^{−1/4} Re^{−3/4}, lambdaT_over_L = sqrt(30 or 15 · u2_over_dU2/(c Re)),
    uK_over_dU = (c/Re)^{1/4}, tau_eta_over_T = (c Re)^{−1/2} with the outer time T = L/ΔU, R_lambda = λ·sqrt(u2)/ν in
    the same units) — floats for a float Re_L.  (``lambda_over_L`` and ``tau_eta_dU_over_L`` are the same two numbers
    under their earlier names.)
    Assumptions: ε̄ = c_eps (ΔU)³/L (stationary, high Reynolds number); isotropic relation (12.43) for the microscale.
    Validation: V1 identical to kolmogorov_scales(ν, c ΔU³/L)/L for any ΔU, L, ν; slopes −¾ and −½ on log axes.
    """
    Re = _F(Re_L)
    k = {"f": 30.0, "g": 15.0}[microscale]
    lam = np.sqrt(k * u2_over_dU2 / (c_eps * Re))  # Eq. (12.52) with the factor of (12.43)
    tau = (c_eps * Re) ** -0.5
    return dict(eta_over_L=_S(c_eps ** -0.25 * Re ** -0.75),  # Eq. (12.51)
                lambdaT_over_L=_S(lam), uK_over_dU=_S((c_eps / Re) ** 0.25), tau_eta_over_T=_S(tau),
                R_lambda=_S(lam * np.sqrt(u2_over_dU2) * Re), lambda_over_L=_S(lam), tau_eta_dU_over_L=_S(tau))


def inertial_range_decades(Re_L, c_eps: float = 1.0):
    """Width of the range of scales between L and η in decades: ``log10(L/η) = ¾ log10 Re_L`` (+ ¼ log10 c_eps).

    Book: §12.7, Eq. (12.51).  Returns decades [–] (5.25 for Re_L = 10⁷).  Scalar-callable.
    Assumptions: eps = c_eps dU^3/L; high Reynolds number.
    Validation: V1 5.25 for Re_L = 1e7; slope 3/4 in log10 Re_L.
    """
    return _S(-np.log10(_F(scale_separation(Re_L, c_eps)["eta_over_L"])))


def dns_grid_points(Re_L, c_eps: float = 1.0, points_per_eta: float = 1.0):
    """Grid points needed to resolve every scale from L down to η in three dimensions: ``(L/η)³ ∝ Re_L^{9/4}``.

    Book: §12.7, Eq. (12.51) (the cost estimate is ours; it is why direct simulation is one of three responses to the
    closure problem, §12.5).  Returns the number of points [–].  Scalar-callable.
    Assumptions: uniform grid of spacing eta/points_per_eta in a cube of side L (an estimate of the cost, not a resolution rule).
    Validation: V1 3.16e13 for Re_L = 1e6; slope 9/4 in log Re_L.
    """
    return _S((points_per_eta / _F(scale_separation(Re_L, c_eps)["eta_over_L"])) ** 3)


def scale_ordering(Re_L, urms_over_dU: float = 1.0, c_eps: float = 1.0) -> dict:
    """The ordering η < λ_T < L of the length scales at a given outer Reynolds number.

    Book: §12.7 (after (12.52)): at high Reynolds number η < λ_T < Λ < L, and R_λ ∼ Re_L^{1/2}.
    Parameters: Re_L; urms_over_dU (one-component rms)/ΔU; c_eps.  Returns dict(eta_over_L, lambda_over_L (λ_g),
    R_lambda, ordered (bool: η < λ_g < L)).  The ordering fails at low Re_L (λ_g/L > 1 for Re_L < 15 u2/ΔU²).
    Assumptions: eps = c_eps dU^3/L and the isotropic relation (12.43) for lambda_g.
    Validation: V7 ordered for Re_L >= 100; fails below Re_L = 15 with u_rms = dU.
    """
    s = scale_separation(Re_L, c_eps, urms_over_dU ** 2, microscale="g")
    ordered = (_F(s["eta_over_L"]) < _F(s["lambda_over_L"])) & (_F(s["lambda_over_L"]) < 1.0)
    return dict(eta_over_L=s["eta_over_L"], lambda_over_L=s["lambda_over_L"], R_lambda=s["R_lambda"],
                ordered=bool(ordered) if np.ndim(ordered) == 0 else ordered)


def cascade_tiers(L: float, dU: float, nu: float, ratio: float = 2.0, c_eps: float = 1.0) -> dict:
    """The tiers of the energy cascade from the outer scale L down to the Kolmogorov scale η.

    Book: §12.7 (Richardson's cascade): each tier of eddies is strained by the next larger one; the transfer is inviscid
    while the eddy Reynolds number u′l′/ν ≫ 1 and ends where it is of order one.  The inertial-range estimate
    ``u′(l′) ∼ (ε̄ l′)^{1/3}`` (the same ε̄ passes through every tier) is ours.
    Parameters: L [m]; dU [m/s]; nu [m²/s]; ratio size ratio between successive tiers (> 1; a picture, not a law);
    c_eps constant of (12.49).
    Returns dict(size l′ [m], velocity u′ [m/s], turnover time l′/u′ [s], Re eddy Reynolds number u′l′/ν — four arrays,
    largest tier first, last tier ≥ η — and the floats eps [m²/s³], eta [m], n_tiers, Re_L).  Re falls as l′^{4/3} and
    reaches ≈ 1 at η.
    Assumptions: stationary cascade with the same flux ε̄ through every tier (inertial range); a discrete ladder is a
    cartoon of a continuous range of scales.
    Validation: V1 velocity³/size = ε̄ on every tier; Re of the last tier between 1 and ratio^{4/3}; V7 the number of
    tiers grows as ¾ log(Re_L)/log(ratio).
    """
    eps = float(dissipation_outer_scaling(dU, L, c_eps))
    eta = float(kolmogorov_scales(nu, eps)[0])
    n = int(np.floor(np.log(L / eta) / np.log(ratio))) + 1
    l = L / ratio ** np.arange(max(n, 1))
    u = (eps * l) ** (1.0 / 3.0)
    return dict(size=l, velocity=u, turnover=l / u, Re=u * l / nu, eps=eps, eta=eta, n_tiers=int(l.size),
                Re_L=dU * L / nu)


def richardson_diffusivity(l, eps, c: float = 1.0):
    """Richardson's four-thirds law: the effective diffusivity of a patch of size l, ``K ∼ c ε̄^{1/3} l^{4/3}``.

    Book: §12.2 (K ∝ l^{4/3}); the dimensional form with ε̄ follows from the inertial-range scales of §12.7.
    Parameters: l patch size [m]; eps [m²/s³]; c order-one constant.  Returns K [m²/s].  Scalar-callable.
    Assumptions: patch size inside the inertial range (eta << l << L); an order-of-magnitude scaling law (the constant is an argument).
    Validation: V2 units m2/s; ratio 10^(4/3) per decade of l.
    """
    return _S(c * _F(eps) ** (1.0 / 3.0) * _F(l) ** (4.0 / 3.0))


SCALE_CASES: dict[str, dict] = {
    # our own illustrative cases (L [m], dU [m/s], nu [m²/s]) — round numbers chosen by us; none is a book example
    "kitchen_mixer": dict(L=0.05, dU=1.0, nu=1.0e-6, fluid="water"),
    "wind_tunnel": dict(L=0.05, dU=1.0, nu=1.5e-5, fluid="air"),
    "atmospheric_boundary_layer": dict(L=1000.0, dU=5.0, nu=1.5e-5, fluid="air"),
    "ocean_thermocline": dict(L=10.0, dU=0.1, nu=1.0e-6, fluid="water"),
}


def scale_table(case=None, c_eps: float = 1.0, u2_over_dU2: float = 1.0):
    """All the scales of §12.7 for one flow: Re_L, ε̄, η, u_K, τ_η, λ_T, R_λ, decades of scales, DNS grid points.

    Book: §12.7, Eqs. (12.43) ``ε̄ = 15 ν u2/λ_g²``, (12.44) ``R_λ = λ sqrt(u2)/ν``, (12.49) ``ε̄ ∼ (ΔU)³/L``, (12.50)
    ``η = (ν³/ε̄)^{1/4}``, (12.51) ``η/L ∼ Re_L^{−3/4}`` (pure composition of the functions above).
    Parameters
    ----------
    case : a key of :data:`SCALE_CASES` — "kitchen_mixer", "wind_tunnel", "atmospheric_boundary_layer",
        "ocean_thermocline" (our own illustrative cases) — or a dict(L [m], dU [m/s], nu [m²/s]), or None (→ the list of
        all four cases).
    c_eps : constant of (12.49) [–].   u2_over_dU2 : (one-component variance)/(ΔU)² [–].
    Returns
    -------
    dict of floats (or a list of such dicts): name, fluid, L [m], dU [m/s], nu [m²/s], Re_L, eps [m²/s³], eta [m],
    u_K [m/s], tau_eta [s], lambda_g [m], R_lambda, decades (log10 L/η), grid_points ((L/η)³).
    Assumptions: order-of-magnitude scaling (c_eps = 1, u_rms = ΔU by default); isotropic small scales.
    Validation: V1 equals the component functions; V7 the atmospheric case has η of a fraction of a millimetre.
    """
    if case is None:
        return [scale_table(k, c_eps, u2_over_dU2) for k in SCALE_CASES]
    if isinstance(case, str) and case not in SCALE_CASES:
        raise ValueError(f"case must be one of {tuple(SCALE_CASES)} or a dict(L, dU, nu)")
    p = dict(SCALE_CASES[case]) if isinstance(case, str) else dict(case)
    L, dU, nu = float(p["L"]), float(p["dU"]), float(p["nu"])
    Re = dU * L / nu
    eps = float(dissipation_outer_scaling(dU, L, c_eps))
    eta, uK, tau_eta = kolmogorov_scales(nu, eps)
    u2 = u2_over_dU2 * dU ** 2
    lam_g = float(np.sqrt(15.0 * nu * u2 / eps))  # Eq. (12.43) solved for λ_g
    return dict(name=case if isinstance(case, str) else p.get("name", "custom"), fluid=p.get("fluid", ""), L=L, dU=dU,
                nu=nu, Re_L=Re, eps=eps, eta=eta, u_K=uK, tau_eta=tau_eta, lambda_g=lam_g,
                R_lambda=float(taylor_reynolds_number(u2, lam_g, nu)), decades=float(inertial_range_decades(Re, c_eps)),
                grid_points=float(dns_grid_points(Re, c_eps)))


# ======================================================================================================================
# §12.7  spectra
# ======================================================================================================================
def kolmogorov_constants(C: float = KOLMOGOROV_C) -> dict:
    """The three constants of the −5/3 law and how they are related.

    Book: §12.7 (after (12.55)): the three-dimensional spectrum ``E(K) = C ε̄^{2/3} K^{−5/3}`` with ē = ∫_0^∞ E dK, and
    the one-dimensional longitudinal spectrum of (12.54).  For isotropic turbulence the one-sided one-dimensional
    constant is ``C_1 = (18/55) C`` (ours: :func:`one_dimensional_from_3d` applied to a pure power law); the two-sided
    one (the book's normalisation (12.55)) is half of that.
    Returns dict(C, C1_one_sided, C1_two_sided).  For C = 1.5: 0.491 and 0.245.
    Assumptions: isotropic turbulence; C is an empirical constant (about 1.5: Sreenivasan, Phys. Fluids 7, 2778 (1995)).
    Validation: V1 one_dimensional_from_3d of a pure -5/3 law returns (18/55) C to 1e-9.
    """
    c1 = 18.0 / 55.0 * C
    return dict(C=float(C), C1_one_sided=c1, C1_two_sided=0.5 * c1)


def inertial_spectrum_1d(k1, eps, C1: float | None = None, two_sided: bool = True, printed: bool = False):
    """Kolmogorov's inertial-range form of the one-dimensional longitudinal spectrum: ``S_11 = C_1 ε̄^{2/3} k_1^{−5/3}``.

    Book: §12.7, Eq. (12.54), valid for 2π/L ≪ k_1 ≪ 2π/η.  **The page prints the exponent of k_1 as +5/3 (slip #1).**
    The corrected −5/3 is what the next sentence of the book calls the law, what Fig. 12.12 shows, and what dimensions
    require: ε̄^{2/3} k_1^{−5/3} has the m³/s² of S_11, ε̄^{2/3} k_1^{+5/3} does not.  ``printed=True`` reproduces the page
    so that a test can show it fails.
    Parameters: k1 [rad/m] > 0; eps [m²/s³]; C1 **one-sided** constant (None → (18/55)·1.5 from
    :func:`kolmogorov_constants`); two_sided True → the book's two-sided density (12.55), i.e. C1/2; False → one-sided.
    Returns S_11 [m³/s²].  Scalar-callable.
    Validation: V1 slope −5/3 on log axes; two-sided = one-sided/2; V2 units (``printed=True`` fails the check).
    Assumptions: locally isotropic turbulence at high Reynolds number; k1 inside the inertial range.
    """
    c = kolmogorov_constants()["C1_one_sided"] if C1 is None else float(C1)
    if two_sided:
        c *= 0.5                                   # Eq. (12.55)
    # DEVIATION: exponent −5/3, not the printed +5/3 (slip #1) — the printed form has the wrong units and the wrong
    # slope; printed=True keeps it available so that a test can show it fails.
    exponent = 5.0 / 3.0 if printed else -5.0 / 3.0
    return _S(c * _F(eps) ** (2.0 / 3.0) * _F(k1) ** exponent)  # Eq. (12.54), exponent corrected


def inertial_spectrum_3d(K, eps, C: float = KOLMOGOROV_C):
    """Three-dimensional inertial-range spectrum ``E(K) = C ε̄^{2/3} K^{−5/3}`` (ē = ∫_0^∞ E dK).

    Book: §12.7 (the three-dimensional form of (12.54), stated after (12.55)).  K [rad/m]; eps [m²/s³].  Returns [m³/s²].
    Assumptions: locally isotropic turbulence at high Reynolds number; K inside the inertial range.
    Validation: V2 units m3/s2; V1 slope -5/3.
    """
    return _S(C * _F(eps) ** (2.0 / 3.0) * _F(K) ** (-5.0 / 3.0))


def kolmogorov_normalize_spectrum(k1, S11, nu: float, eps: float):
    """Put a one-dimensional spectrum into Kolmogorov scaling: ``(k_1 η, S_11/(u_K² η))``.

    Book: §12.7, Eq. (12.53) ``S_11/(ν^{5/4} ε̄^{1/4}) = Φ(k_1 ν^{3/4}/ε̄^{1/4})``, i.e. ``S_11/(u_K² η) = Φ(k_1 η)`` for
    k_1 ≫ 2π/L: spectra from different flows collapse at high wavenumber (Fig. 12.12, which plots 2S_11).
    Parameters: k1 [rad/m]; S11 [m³/s²]; nu [m²/s]; eps [m²/s³].  Returns (k1·eta, S11/(u_K²·eta)) [–].
    In these variables the inertial range is Φ = C_1 (k_1 η)^{−5/3} whatever ν and ε̄ are.
    Assumptions: nu and eps of the flow are known (a change of variables).
    Validation: V1 the inertial law becomes C1 (k1 eta)^(-5/3) for any nu, eps.
    """
    eta, uK, _ = kolmogorov_scales(nu, eps)
    return _S(_F(k1) * eta), _S(_F(S11) / (uK ** 2 * eta))  # Eq. (12.53)


def fit_inertial_range(k, S, band: tuple, eps: float | None = None) -> tuple:
    """Least-squares power law ``S = a k^p`` over a band of wavenumbers (straight line on log–log axes).

    Book: §12.7, Eq. (12.54) ``S_11 = C_1 ε̄^{2/3} k_1^{−5/3}`` (exponent corrected, slip #1) and Fig. 12.12 (the −5/3
    line).
    Parameters
    ----------
    k : (n,) wavenumbers [rad/m] > 0.   S : (n,) spectrum [m³/s²].   band : (k_min, k_max) [rad/m], inclusive.
    eps : dissipation rate [m²/s³] or None.
    Returns
    -------
    (slope, constant): slope p of the fitted line [–]; constant = the prefactor a of ``S = a k^p`` [units of S·k^{−p}]
    when ``eps`` is None, or — when ``eps`` is given — the dimensionless Kolmogorov-type constant of a forced −5/3 law,
    ``geometric mean of S k^{5/3} over the band / ε̄^{2/3}``.
    Assumptions: a power law holds across the band (the caller chooses the band); points with S ≤ 0 are ignored.
    Validation: V1 slope −5/3 and constant C recovered from :func:`inertial_spectrum_3d` to 1e-12; slope within 0.02 of
    −5/3 on :func:`model_spectrum` well inside its inertial range.
    """
    k, S = _F(k), _F(S)
    m = (k >= band[0]) & (k <= band[1]) & (S > 0)
    if np.count_nonzero(m) < 2:
        raise ValueError("fewer than two points in the band")
    p, lna = np.polyfit(np.log(k[m]), np.log(S[m]), 1)
    if eps is None:
        return float(p), float(np.exp(lna))
    return float(p), float(np.exp(np.mean(np.log(S[m] * k[m] ** (5.0 / 3.0)))) / eps ** (2.0 / 3.0))  # Eq. (12.54)


def model_spectrum(K, eps: float, nu: float, L: float | None = None, kind: str = "pao", C: float = KOLMOGOROV_C,
                   c_L: float = 6.78, p0: float = 2.0):
    """A model three-dimensional energy spectrum E(K) for the whole range of scales.

    Book: §12.7 names Pao's (1965) form as the curve that fits the roll-off of Fig. 12.12; the formulas are from the
    literature, not from the book:
    ``E(K) = C ε̄^{2/3} K^{−5/3} f_L(KL) exp{−(3/2) C (Kη)^{4/3}}`` (kind "pao"), with
    ``f_L = (KL/sqrt((KL)² + c_L))^{5/3 + p0}`` when an outer scale L is given (a von Kármán-type energy-containing range,
    E ∝ K^{p0} at small K; the form used by Pope 2000, §6.5) and f_L = 1 when L is None.
    Parameters: K [rad/m] > 0; eps [m²/s³]; nu [m²/s]; L [m] or None; kind "pao" (Pao's dissipation range; the
    energy-range factor f_L is applied whenever L is given) or "pope" (the same curve under the name the notebook and
    the explainer use for "Pao's roll-off with Pope's energy-range factor"; **L is then required**); C
    three-dimensional Kolmogorov constant; c_L
    roll-off constant of the large scales (**a shape parameter, not a measured constant**: 6.78 is the value attributed to
    Pope's model spectrum and is unverified from an open source — it only sets where the spectrum bends over); p0
    low-wavenumber exponent.
    Returns E [m³/s²].  Scalar-callable.
    Assumptions: isotropic turbulence; one universal shape (no bottleneck, no intermittency correction).
    Property (exact for L = None): the dissipation integral ``2ν ∫_0^∞ K² E dK`` equals ε̄ — the exponential cut-off is
    the one that makes the model consistent with its own ε̄.
    Validation: V4 that integral to 1e-8 (L = None) and within 1 % with L when L/η ≳ 10³; V1 −5/3 slope in between.
    Label for the large-scale shape: qualitative.
    """
    if kind not in ("pao", "pope"):
        raise ValueError("kind must be 'pao' or 'pope'")
    if kind == "pope" and L is None:
        raise ValueError("kind='pope' (Pao's roll-off with the energy-range factor) needs the outer scale L")
    K_ = _F(K)
    eta = (nu ** 3 / eps) ** 0.25
    E = C * eps ** (2.0 / 3.0) * K_ ** (-5.0 / 3.0) * np.exp(-1.5 * C * (K_ * eta) ** (4.0 / 3.0))
    if L is not None:
        kl = K_ * L
        E = E * (kl / np.sqrt(kl ** 2 + c_L)) ** (5.0 / 3.0 + p0)
    return _S(E)


def one_dimensional_from_3d(k1, E_fn: Callable, K_max: float = np.inf, one_sided: bool = True):
    """One-dimensional longitudinal spectrum of isotropic turbulence from its three-dimensional spectrum.

    Book: §12.6–§12.7 use S_11(k_1) (12.45) and E(K) side by side; the link is the standard isotropic relation (ours):
    ``E_11(k_1) = ∫_{k_1}^∞ (E(K)/K)(1 − k_1²/K²) dK`` (one-sided: ∫_0^∞ E_11 dk_1 = mean(u_1²)); the book's two-sided
    S_11 is half of it.
    Parameters: k1 scalar or array [rad/m] > 0; E_fn callable K → E [m³/s²]; K_max upper limit; one_sided.
    Returns E_11 (or S_11) [m³/s²].  ``quad`` in ln K (the integrand is smooth on a log axis), rtol 1e-10.
    Validation: V1 a pure K^{−5/3} law returns (18/55) C ε̄^{2/3} k_1^{−5/3}; ∫E_11 dk_1 = (2/3)∫E dK.
    Assumptions: isotropic turbulence; E_fn decays fast enough for the integral to converge.
    """
    k1a = np.atleast_1d(_F(k1)).astype(float)
    out = np.empty(k1a.shape)
    for idx, k in np.ndenumerate(k1a):
        f = lambda s, k=k: float(E_fn(np.exp(s))) * (1.0 - k * k * np.exp(-2.0 * s))  # noqa: E731  dK/K = ds
        hi = np.log(K_max) if np.isfinite(K_max) else np.log(k) + 60.0
        val, _ = quad(f, np.log(k), hi, epsabs=0.0, epsrel=1e-10, limit=400)
        out[idx] = val
    if not one_sided:
        out = 0.5 * out
    return _S(out.reshape(np.shape(k1)))


# ======================================================================================================================
# §12.8  free turbulent shear flows — plane jet
# ======================================================================================================================
_LN2 = float(np.log(2.0))


def gaussian_profile(xi, xi_half: float):
    """Gaussian fit of a self-similar profile specified by its half-width: ``F(ξ) = exp{−ln 2 · ξ²/ξ_½²}``.

    Book: §12.8 (unnumbered, after (12.73)); ξ_½ is the one-sided half-width — the value of ξ where F = ½.
    Parameters: xi similarity variable [–]; xi_half > 0 [–] (**required**: an empirical, flow-dependent number).
    Returns F [–]; F(0) = 1, F(ξ_½) = ½, even in ξ.  Scalar-callable.
    Assumptions: a fit to measured profiles (it is not a solution of the similarity equation).
    Validation: V1 F(0) = 1, F(xi_half) = 1/2, even.
    """
    return _S(np.exp(-_LN2 * (_F(xi) / xi_half) ** 2))


def profile_integrals(xi_half_U: float, xi_half_Y: float | None = None) -> dict:
    """Closed-form integrals of the Gaussian profiles that appear in the invariants of the plane jet.

    Book: §12.8, the dimensionless integrals of (12.65) ∫F² dξ, (12.68) ∫F dξ and (12.70) ∫HF dξ (limits ±∞).
    With F = exp(−a ξ²), a = ln 2/ξ_½²: ∫F = sqrt(π/a) = ξ_½ sqrt(π/ln 2), ∫F² = sqrt(π/2a) = ξ_½ sqrt(π/(2 ln 2)),
    ∫HF = sqrt(π/(a_U + a_Y)).
    Parameters: xi_half_U half-width of the velocity profile F [–] (**required**, empirical); xi_half_Y half-width of
    the scalar profile H [–] or None.
    Returns dict(I1 = ∫F dξ, I2 = ∫F² dξ, IHF = ∫HF dξ (None without xi_half_Y)) — floats [–].
    Assumptions: Gaussian profile shapes (a fit, not a solution of the similarity equation).
    Validation: V1 against ``scipy.integrate.quad`` to 1e-12; I2 = I1/√2.
    """
    aU = _LN2 / xi_half_U ** 2
    out = dict(I1=float(np.sqrt(np.pi / aU)), I2=float(np.sqrt(np.pi / (2.0 * aU))), IHF=None)
    if xi_half_Y is not None:
        out["IHF"] = float(np.sqrt(np.pi / (aU + _LN2 / xi_half_Y ** 2)))
    return out


def _plane_jet_C5(C5, xi_half: float) -> float:
    if isinstance(C5, str):
        if C5 != "from_invariant":
            raise ValueError("C5 must be a number or 'from_invariant'")
        return 1.0 / np.sqrt(profile_integrals(xi_half)["I2"])  # J_s = ρ C5² (J_s/ρ) ∫F²dξ   (12.65)–(12.66)
    return float(C5)


def plane_jet_centerline_velocity(x, Js: float, rho: float, *, C5, xi_half: float, x0: float = 0.0):
    """Centreline mean velocity of the plane turbulent jet: ``U_CL = C_5 (J_s/ρ)^{1/2} (x − x0)^{−1/2}``.

    Book: §12.8, Eq. (12.66) at y = 0.  Parameters as :func:`plane_jet_mean_velocity`.  Returns [m/s].
    Assumptions: far field of a high-Reynolds-number plane jet in still fluid of the same density; Gaussian profile shape.
    Validation: V1 slope -1/2 in log x; equals plane_jet_mean_velocity at y = 0.
    """
    return _S(_plane_jet_C5(C5, xi_half) * np.sqrt(Js / rho) * (_F(x) - x0) ** -0.5)


def plane_jet_mean_velocity(x, y, Js: float, rho: float, *, C5, xi_half: float, x0: float = 0.0):
    """Far-field mean streamwise velocity of the plane turbulent jet: ``U = C_5 (J_s/ρ)^{1/2} x^{−1/2} F(y/x)``.

    Book: §12.8, Eqs. (12.56) (self-preserving form) and (12.66) (after the invariant (12.62) fixes the exponent
    2γ + 1 = 0 in (12.65)); δ = x by the convention C_2 − C_1 = 1, the virtual origin x0 restored as an argument.
    Parameters: x [m] > x0 downstream distance; y [m] cross-stream; Js momentum flux per unit span [N/m] = [kg/s²];
    rho ambient density [kg/m³]; C5 amplitude constant — a number (**required**, empirical) or "from_invariant", which
    sets C5 = (∫F²dξ)^{−1/2} so that the Gaussian profile carries exactly J_s; xi_half half-width of F in ξ = y/x
    (**required**); x0 virtual origin [m].
    Returns U [m/s].  Profile shape: Gaussian (:func:`gaussian_profile`).  Scalar-callable.
    Assumptions: high Reynolds number, far field, quiescent surroundings of the same density.
    Validation: V4 ρ∫U²dy independent of x and equal to J_s (with "from_invariant"); V1 U_CL ∝ x^{−1/2}, width ∝ x.
    """
    xr = _F(x) - x0
    return _S(_F(plane_jet_centerline_velocity(x, Js, rho, C5=C5, xi_half=xi_half, x0=x0))
              * _F(gaussian_profile(_F(y) / xr, xi_half)))  # Eq. (12.66)


def _gauss_int(xi, xi_half):
    """I(ξ) = ∫_0^ξ F dξ for the Gaussian profile."""
    a = _LN2 / xi_half ** 2
    return 0.5 * np.sqrt(np.pi / a) * erf(np.sqrt(a) * _F(xi))


def plane_jet_cross_velocity(x, y, Js: float, rho: float, *, C5, xi_half: float, x0: float = 0.0):
    """Mean cross-stream velocity of the plane jet from continuity: ``V = −∫_0^y ∂U/∂x dy = U_CL [ξF − ½∫_0^ξ F dξ]``.

    Book: §12.8, Eq. (12.58) integrated with V(0) = 0 (the step before (12.63)); closed form for (12.66), ours.
    Parameters as :func:`plane_jet_mean_velocity`.  Returns V [m/s]: outward near the axis, inward (entrainment) at the
    edges, V(±∞) = ∓¼ U_CL ∫F dξ.  Scalar-callable.
    Validation: V1 ∂U/∂x + ∂V/∂y = 0 to O(h²) (:func:`rans_2d_residual`); V(∞) = −½ dV̇/dx.
    Assumptions: far field of a high-Reynolds-number plane jet in still fluid of the same density; Gaussian profile shape.
    """
    xr = _F(x) - x0
    xi = _F(y) / xr
    Ucl = _F(plane_jet_centerline_velocity(x, Js, rho, C5=C5, xi_half=xi_half, x0=x0))
    return _S(Ucl * (xi * _F(gaussian_profile(xi, xi_half)) - 0.5 * _gauss_int(xi, xi_half)))


def plane_jet_stress_profile(xi, F: Callable | None = None, C3: float = 1.0, xi_half: float | None = None):
    """Profile function G(ξ) of the Reynolds shear stress of the plane jet, from the similarity equation.

    Book: §12.8, Eq. (12.63) ``{δU′_CL/U_CL} F² − {δU′_CL/U_CL + δ′} F′ ∫_0^ξ F dξ = {Ψ/U_CL²} G′`` with (12.64) and the
    plane-jet values δ = x, δU′_CL/U_CL = −½, δU′_CL/U_CL + δ′ = +½, Ψ/U_CL² = C_3.  The left side is then
    −½ d(F∫F)/dξ, so one integration with G(0) = 0 gives (ours)  **C_3 G(ξ) = −½ F(ξ) ∫_0^ξ F dξ**.
    Sign: −mean(uv) = Ψ G is negative for ξ > 0, where ∂U/∂y < 0 — the stress has the sign of the mean shear, as an eddy
    viscosity would give.  G is odd, G(0) = 0, G → 0.
    Parameters: xi [–]; F callable ξ → F (its integral by ``quad``) or None with ``xi_half`` (Gaussian, closed form);
    C3 = Ψ/U_CL² (it only rescales G: the stress C_3 G is fixed by F).
    Returns G [–].  Validation: V2 (12.63) residual 0 by sympy for a Gaussian F; V1 G odd.
    Assumptions: thin-layer similarity equation (12.63) of the plane jet with delta = x; viscous stress neglected.
    """
    xi_ = _F(xi)
    if F is None:
        if xi_half is None:
            raise ValueError("give F or xi_half")
        Fv, Iv = _F(gaussian_profile(xi_, xi_half)), _gauss_int(xi_, xi_half)
    else:
        Fv = _F(F(xi_))
        Iv = np.vectorize(lambda s: quad(lambda q: float(F(q)), 0.0, s, epsabs=1e-13, epsrel=1e-12)[0])(xi_)
    # DEVIATION (from analysis/ch12.md row 121, not from the book): the analysis note writes C_3 G = +½ F∫F; carrying
    # the signs of (12.63) through (c1 = −½, c2 = +½) gives −½ F∫F, confirmed by plane_jet_similarity_sympy and by the
    # thin-layer residual of thin_shear_layer_terms.  −mean(uv) then has the sign of ∂U/∂y, as it must.
    return _S(-0.5 * Fv * Iv / C3)  # Eq. (12.63) integrated once


def plane_jet_reynolds_stress(x, y, Js: float, rho: float, *, C5, xi_half: float, C3: float = 1.0, x0: float = 0.0):
    """Far-field Reynolds shear stress of the plane jet, ``−mean(uv) = C_3 U_CL² G(y/x)`` [m²/s²].

    Book: §12.8, Eqs. (12.57) ``−mean(uv) = Ψ(x) G(ξ)`` and (12.67) ``−mean(uv) = C_3 C_5² (J_s/ρ) x^{−1} G(y/x)``; with
    G from :func:`plane_jet_stress_profile` the product C_3 G does not depend on C_3:
    ``−mean(uv) = −½ U_CL² F ∫_0^ξ F``.
    Parameters as :func:`plane_jet_mean_velocity`, plus C3 = Ψ/U_CL² [–] (it scales G by 1/C_3 and the amplitude by
    C_3, so **the returned stress is the same for every C3** — kept as an argument to show exactly that).
    Returns −mean(uv) [m²/s²] (negative above the axis, positive below; zero on the axis; decays as 1/x).
    Scalar-callable.
    Assumptions: as :func:`plane_jet_mean_velocity` (far field, high Reynolds number, Gaussian F, thin layer).
    Validation: V1 thin-layer momentum residual (12.61) → 0; independent of C3 to round-off; odd in y.
    """
    xr = _F(x) - x0
    xi = _F(y) / xr
    Ucl = _F(plane_jet_centerline_velocity(x, Js, rho, C5=C5, xi_half=xi_half, x0=x0))
    return _S(C3 * Ucl ** 2 * _F(plane_jet_stress_profile(xi, xi_half=xi_half, C3=C3)))  # Eq. (12.67)


def plane_jet_volume_flux(x, Js: float, rho: float, *, C5, xi_half: float, x0: float = 0.0):
    """Volume flux per unit span of the plane jet: ``V̇ = C_5 (J_s/ρ)^{1/2} x^{1/2} ∫F dξ`` — it grows downstream.

    Book: §12.8, Eq. (12.68).  Returns V̇ [m²/s].  The growth ∝ x^{1/2} is the entrainment of ambient fluid.
    Assumptions: far field of a high-Reynolds-number plane jet in still fluid of the same density; Gaussian profile shape.
    Validation: V1 equals the trapezoid integral of plane_jet_mean_velocity to 1e-10; doubles from x to 4x.
    """
    return _S(_plane_jet_C5(C5, xi_half) * np.sqrt(Js / rho) * (_F(x) - x0) ** 0.5 * profile_integrals(xi_half)["I1"])  # Eq. (12.68)


def plane_jet_entrainment_velocity(x, Js: float, rho: float, *, C5, xi_half: float, x0: float = 0.0):
    """Speed at which ambient fluid flows toward the plane jet on each side: ``v_e = ½ dV̇/dx = −V(+∞)``.

    Book: §12.8, Eq. (12.68) differentiated (ours).  Returns v_e [m/s] > 0, ∝ x^{−1/2}: a fixed fraction
    (¼∫F dξ) of the centreline velocity.
    Assumptions: far field of a high-Reynolds-number plane jet in still fluid of the same density; Gaussian profile shape.
    Validation: V1 equals -V(y -> inf) of plane_jet_cross_velocity and half the x-derivative of plane_jet_volume_flux.
    """
    return _S(0.25 * _plane_jet_C5(C5, xi_half) * np.sqrt(Js / rho) * (_F(x) - x0) ** -0.5 * profile_integrals(xi_half)["I1"])


def plane_jet_mass_fraction(x, y, Ms: float, Js: float, rho: float, *, C6, xi_half_Y: float, C5=None,
                            xi_half: float | None = None, x0: float = 0.0):
    """Far-field mean mass fraction of slot fluid in the plane jet: ``Ȳ = C_6 (Ṁ_s/sqrt(ρ J_s)) x^{−1/2} H(y/x)``.

    Book: §12.8, Eqs. (12.69) ``Ȳ = Y_CL(x) H(y/x)``, (12.70) ``Ṁ_s ≅ ρ ∫ȲU dy = ρ Y_CL C_5 (J_s/ρ)^{1/2} x^{1/2} ∫HF dξ``
    and (12.71).  (The printed subscript y = 0 on the source integral of (12.70) should read x = 0 — slip #3.)
    Parameters
    ----------
    x : downstream distance [m] > x0.   y : cross-stream coordinate [m].
    Ms : slot mass flux per unit span Ṁ_s [kg/(m s)].   Js : momentum flux per unit span [N/m].   rho : ambient
        density [kg/m³].
    C6 : amplitude constant [–] — a number (**required keyword**, empirical) or "from_invariant", which sets
        ``C_6 = 1/(C_5 ∫HF dξ)`` so that (12.70) holds exactly for Gaussian F and H (then ``C5`` and ``xi_half`` of the
        velocity field must be given too; ``C5`` may itself be "from_invariant").
    xi_half_Y : half-width of H in ξ = y/x [–] (**required keyword**).
    C5, xi_half : velocity-field constants, needed only with ``C6="from_invariant"``.   x0 : virtual origin [m].
    Returns
    -------
    Ȳ [–] (kg of slot fluid per kg of mixture); a float for scalar input.
    Assumptions: far field, high Reynolds number; streamwise turbulent scalar flux mean(uY′) neglected against ȲU, as
    in (12.70); Gaussian H.
    Validation: V4 ρ∫ȲU dy = Ṁ_s at every x (with "from_invariant"); V1 Y_CL ∝ x^{−1/2}.
    """
    if isinstance(C6, str):
        if C6 != "from_invariant":
            raise ValueError("C6 must be a number or 'from_invariant'")
        if C5 is None or xi_half is None:
            raise ValueError("C6='from_invariant' needs C5 and xi_half of the velocity field")
        c6 = 1.0 / (_plane_jet_C5(C5, xi_half) * profile_integrals(xi_half, xi_half_Y)["IHF"])  # Eq. (12.70)
    else:
        c6 = float(C6)
    xr = _F(x) - x0
    return _S(c6 * Ms / np.sqrt(rho * Js) * xr ** -0.5 * _F(gaussian_profile(_F(y) / xr, xi_half_Y)))  # Eq. (12.71)


def jet_momentum_flux_per_span(y, U, rho: float) -> float:
    """Momentum flux per unit span ``J = ρ ∫U² dy`` [N/m] of a sampled profile (trapezoid) — the jet's invariant.

    Book: §12.8, Eq. (12.62) ``J_s ≅ ρ ∫U² dy = const``.  Wrapper of ``core.jets.jet_momentum_flux`` (ch09 (9.58)).
    Returns J [N/m] (float).
    Assumptions: the profile has decayed at both ends of the grid; constant density.
    Validation: V4 independent of x for plane_jet_mean_velocity; V1 equals core.jets.jet_momentum_flux.
    """
    return jet_momentum_flux(y, U, rho)  # Eq. (12.62)


def scalar_flux_per_span(y, U, Y, rho: float) -> float:
    """Flux of slot fluid per unit span ``ρ ∫ȲU dy`` [kg/(m s)] of sampled profiles (trapezoid).  Book: §12.8, Eq. (12.70).

    Returns the flux [kg/(m s)] (float).
    Assumptions: the profiles have decayed at both ends of the grid; streamwise turbulent flux neglected.
    Validation: V4 equals the slot mass flux at every x for the similarity solution.
    """
    return float(rho * np.trapezoid(_F(Y) * _F(U), _F(y)))  # Eq. (12.70)


def slot_momentum_flux(rho_s: float, U0: float, d: float) -> float:
    """Momentum flux per unit span of a uniform slot exit: ``J_s = ρ_s U0² d`` [N/m].  Book: §12.8 (before (12.72)).

    Returns J_s [N/m] (float).
    Assumptions: uniform (top-hat) exit velocity across a slot of width d.
    Validation: V1 units N/m; a hand value.
    """
    return float(rho_s * U0 ** 2 * d)


def slot_mass_flux(rho_s: float, U0: float, d: float) -> float:
    """Mass flux per unit span of a uniform slot exit: ``Ṁ_s = ρ_s U0 d`` [kg/(m s)].  Book: §12.8 (before (12.72)).

    Returns M_s [kg/(m s)] (float).
    Assumptions: uniform (top-hat) exit velocity across a slot of width d.
    Validation: V1 units kg/(m s); a hand value.
    """
    return float(rho_s * U0 * d)


def virtual_origin_fit(x, half_width) -> tuple:
    """Virtual origin and spreading rate from measured half-widths: fit ``y_½ = S (x − x0)``.

    Book: §12.8 (after (12.64)): ``δ(x) = (C_2 − C_1)(x − x_o)``, with x_o of the order of the slot width.
    Parameters: x (n,) stations [m]; half_width (n,) measured half-widths y_½ [m].
    Returns (slope S = dy_½/dx [–], x0 [m]) — the spreading rate and the virtual origin, floats.
    Assumptions: the stations lie in the self-similar far field (linear growth); ordinary least squares.
    Validation: V1 recovers (S, x0) of an exact straight line to 1e-12.
    """
    S, c = np.polyfit(_F(x), _F(half_width), 1)
    return float(S), float(-c / S)


def thin_shear_layer_terms(x: float, y, Js: float, rho: float, nu: float, *, C5, xi_half: float, h: float = 1e-5) -> dict:
    """Sizes of the terms kept and dropped in the thin-layer momentum equation of the plane jet.

    Book: §12.8, Eq. (12.61) ``U U_x + V U_y ≅ −∂mean(uv)/∂y`` obtained from (12.59) for U ≫ V, ∂/∂y ≫ ∂/∂x, no imposed
    pressure gradient and negligible viscous stress.
    Parameters: x [m] station; y (n,) [m]; Js, rho, nu; C5, xi_half as before; h relative difference step.
    Returns dict(advection_x, advection_y, stress_gradient (= −∂mean(uv)/∂y), residual = advection − stress_gradient
    (zero up to differencing error: the similarity solution satisfies (12.61) exactly), viscous_y = ν U_yy,
    viscous_x = ν U_xx, viscous_over_stress (max ratio — of order 1/Re_local), ratio_V_over_U).
    Assumptions: far field of a high-Reynolds-number plane jet in still fluid of the same density; Gaussian profile shape.
    Validation: V1 residual of (12.61) at the differencing error (~1e-8 of the stress gradient); viscous/turbulent ratio ~ 1/Re.
    """
    y = _F(y)
    kw = dict(C5=C5, xi_half=xi_half)
    U = lambda xx, yy: _F(plane_jet_mean_velocity(xx, yy, Js, rho, **kw))  # noqa: E731
    V = lambda xx, yy: _F(plane_jet_cross_velocity(xx, yy, Js, rho, **kw))  # noqa: E731
    muv = lambda xx, yy: _F(plane_jet_reynolds_stress(xx, yy, Js, rho, **kw))  # noqa: E731  (−uv)
    hx, hy = h * x, h * x
    Ux = (U(x + hx, y) - U(x - hx, y)) / (2 * hx)
    Uy = (U(x, y + hy) - U(x, y - hy)) / (2 * hy)
    Uxx = (U(x + hx, y) - 2 * U(x, y) + U(x - hx, y)) / hx ** 2
    Uyy = (U(x, y + hy) - 2 * U(x, y) + U(x, y - hy)) / hy ** 2
    stress = (muv(x, y + hy) - muv(x, y - hy)) / (2 * hy)  # ∂(−uv)/∂y
    adv_x, adv_y = U(x, y) * Ux, V(x, y) * Uy
    return dict(advection_x=adv_x, advection_y=adv_y, stress_gradient=stress, residual=adv_x + adv_y - stress,
                viscous_y=nu * Uyy, viscous_x=nu * Uxx,
                viscous_over_stress=float(np.max(np.abs(nu * Uyy)) / np.max(np.abs(stress))),
                ratio_V_over_U=float(np.max(np.abs(V(x, y))) / np.max(np.abs(U(x, y)))))


def plane_jet_eddy_viscosity_profile(xi, xi_half: float) -> dict:
    """Profile of the plane jet closed with an eddy viscosity that is uniform across the layer: ``F = sech²(a ξ)``.

    Book: §12.8 with §12.10 (12.94), (12.98) (ours): ν_T = ν̂ U_CL δ put into (12.63) gives −½ F∫F = ν̂ F′, whose
    solution is the laminar-jet shape sech² (ch09) — but with the turbulent exponents, because ν_T grows as x^{1/2}.
    Parameters: xi; xi_half half-width (**required**).  Returns dict(F, a = arccosh(√2)/ξ_½, nu_hat = ν_T/(U_CL x) =
    1/(4a²)).  Compared with the Gaussian fit it has fatter tails.
    Assumptions: eddy viscosity uniform across the layer (a closure assumption, section 12.10).
    Validation: V2 F = sech^2 satisfies -F int F / 2 = nu_hat F' exactly; F(xi_half) = 1/2.
    """
    a = float(np.arccosh(np.sqrt(2.0)) / xi_half)
    return dict(F=_S(1.0 / np.cosh(a * _F(xi)) ** 2), a=a, nu_hat=1.0 / (4.0 * a * a))


def jet_tke_budget(xi, xi_half: float, C3: float | None = None, model: str = "eddy_viscosity", e_hat: float = 0.05,
                   b: float = 1.0, C_mu: float = 0.09) -> dict:
    """A *model* of the turbulent-kinetic-energy budget across the plane jet (label: qualitative).

    Book: §12.8, Eq. (12.75) ``0 = −U ∂ē/∂x − V ∂ē/∂y − mean(uv) ∂U/∂y − ∂/∂y(mean(pv)/ρ0 + ½ mean(u_i² v)) − ε̄`` and
    Fig. 12.15.  (The printed triple-correlation term is ½ mean(ē v); with (12.47) it is ½ mean(u_i² v) — slip #10.)
    Only the production term follows from the similarity solution; ē and ε̄ are **assumed**:
    ``ē = e_hat U_CL² E(ξ)``, E = F(ξ)(1 + b ξ²/ξ_½²) (off-axis maximum, as in Fig. 12.14), and
    ``ε̄ = C_mu ē²/ν_T`` with ν_T = ν̂ U_CL x from :func:`plane_jet_eddy_viscosity_profile`; transport is the residual.
    Parameters: xi (n,); xi_half (**required**); C3 unused (kept for the signature of the plan); e_hat, b, C_mu model
    numbers (illustrative, no claim).
    Returns dict(advection, production, dissipation, transport — in units of U_CL³/x — and e (ē/U_CL²)); the four terms
    sum to zero by construction; production vanishes on the axis and peaks near the largest mean shear.
    Assumptions: a model: the shapes of e and eps are assumed, the transport term is the residual.
    Validation: qualitative (only the production term follows from the similarity solution; the four terms sum to zero by construction).
    """
    if model != "eddy_viscosity":
        raise ValueError("only model='eddy_viscosity' is implemented")
    xi_ = _F(xi)
    a = _LN2 / xi_half ** 2
    F = np.exp(-a * xi_ ** 2)
    Fp = -2.0 * a * xi_ * F
    I = _gauss_int(xi_, xi_half)
    s = b / xi_half ** 2
    E = F * (1.0 + s * xi_ ** 2)
    Ep = Fp * (1.0 + s * xi_ ** 2) + F * 2.0 * s * xi_
    nu_hat = plane_jet_eddy_viscosity_profile(0.0, xi_half)["nu_hat"]
    adv = e_hat * (F * E + 0.5 * I * Ep)          # −U ∂ē/∂x − V ∂ē/∂y  with ē ∝ x^{−1} E(ξ)
    prod = (-0.5 * F * I) * Fp                    # (−uv) ∂U/∂y
    diss = -C_mu * (e_hat * E) ** 2 / nu_hat
    return dict(xi=xi_, advection=adv, production=prod, dissipation=diss, transport=-(adv + prod + diss), e=e_hat * E,
                label="qualitative (model)")


# ======================================================================================================================
# §12.8  free shear flows in general: exponents and Table 12.1 forms
# ======================================================================================================================
def general_similarity_check(delta_fn: Callable, UCL_fn: Callable, Psi_fn: Callable, x, h: float = 1e-6) -> dict:
    """The three coefficients of the plane-jet similarity equation for trial functions δ(x), U_CL(x), Ψ(x).

    Book: §12.8, Eq. (12.63) coefficients ``c1 = δ U′_CL/U_CL``, ``c2 = δ U′_CL/U_CL + δ′``, ``c3 = Ψ/U_CL²``; simple
    similarity (12.64) needs each constant, the more general (12.74) only that they share one x-dependence.
    Parameters: three callables of x; x scalar or array [m]; h relative step of the central differences.
    Returns dict(c1, c2, c3, c1_over_c2, c1_over_c3, UCL2_delta = U_CL² δ [m³/s²] (∝ the momentum flux (12.62): it must
    not vary with x), momentum_flux_exponent = d ln(U_CL² δ)/d ln x) — floats for a float x.
    Assumptions: the trial functions are smooth; central differences with relative step h.
    Validation: V1 power family → (n, m + n, 1)·x^{m−1} to 1e-8; exponential family → c2 = 0 to 1e-8 and a non-zero
    momentum-flux exponent; V2 the same coefficients from :func:`plane_jet_similarity_sympy`.
    Reading: the power family δ ∼ x^m, U_CL ∼ x^n, Ψ ∼ x^{2n+m−1} gives (n, m + n, 1)·x^{m−1}: proportional, and the
    momentum flux is constant only for m + 2n = 0.  **The exponential family δ ∼ e^{ax}, U_CL ∼ e^{−ax} that the page
    also offers makes c2 identically zero** — the three are then not proportional — and breaks (12.62) (slip #5).
    """
    x_ = _F(x)
    hx = h * np.maximum(np.abs(x_), 1.0)
    d1 = lambda f: (_F(f(x_ + hx)) - _F(f(x_ - hx))) / (2.0 * hx)  # noqa: E731
    dl, U, Ps = _F(delta_fn(x_)), _F(UCL_fn(x_)), _F(Psi_fn(x_))
    c1 = dl * d1(UCL_fn) / U
    c2 = c1 + d1(delta_fn)
    c3 = Ps / U ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        out = dict(c1=_S(c1), c2=_S(c2), c3=_S(c3), c1_over_c2=_S(c1 / c2), c1_over_c3=_S(c1 / c3),
                   UCL2_delta=_S(U ** 2 * dl),
                   momentum_flux_exponent=_S(x_ * (2.0 * d1(UCL_fn) / U + d1(delta_fn) / dl)))
    return out


def wrong_exponent_fluxes(x, n: float, m: float = 1.0, x_ref: float = 1.0) -> tuple:
    """Momentum and volume flux of a profile family ``U_CL ∝ x^n``, ``δ ∝ x^m`` relative to their values at x_ref.

    Book: §12.8, Eqs. (12.62) ``J_s = ρ ∫U² dy = const`` with (12.65) — momentum flux ∝ U_CL² δ ∝ x^{2n+m} — and (12.68)
    — volume flux ∝ U_CL δ ∝ x^{n+m}.  Only m + 2n = 0 keeps the momentum flux flat (plane jet: m = 1, n = −½).
    Parameters: x station [m] (in units of x_ref); n decay exponent of U_CL [–]; m growth exponent of δ [–];
    x_ref reference station [m].
    Returns (momentum, volume) = ((x/x_ref)^{2n+m}, (x/x_ref)^{n+m}) [–] — a pair of floats for a float x.
    Assumptions: self-similar profiles of fixed shape (plane geometry: one cross-stream direction).
    Validation: V1 (1, x^{1/2}) for the plane jet (n, m) = (−½, 1); momentum ≠ 1 for any other n with m = 1.
    """
    r = _F(x) / x_ref
    return _S(r ** (2.0 * n + m)), _S(r ** (n + m))  # (12.62)/(12.65); (12.68)


def _solve2(a11, a12, b1, a21, a22, b2):
    det = a11 * a22 - a12 * a21
    return (b1 * a22 - a12 * b2) / det, (a11 * b2 - a21 * b1) / det


def free_shear_exponents(flow: str, return_equations: bool = False) -> dict:
    """Exact power-law exponents of the self-similar free shear flows from each flow's invariant and growth law.

    Book: §12.8, Eqs. (12.64)–(12.65) for the plane jet and Table 12.1 for the others (the two-equation derivation for
    each flow is ours, Exercises 12.25–12.28 in outline).  Width δ ∝ x^m; velocity scale (centreline velocity, or the
    centreline deficit of a wake) ∝ x^n; scalar ∝ x^s.  Two linear equations fix (m, n), solved in exact rationals:

    * plane jet:  momentum flux U²δ constant ⇒ m + 2n = 0;  linear growth δ′ = const ⇒ m = 1.      → (1, −½)
    * round jet:  U²δ² constant ⇒ 2m + 2n = 0;  m = 1.                                             → (1, −1)
    * plane wake: momentum deficit U_∞ΔU δ constant ⇒ m + n = 0;  dδ/dx ∼ ΔU/U_∞ ⇒ m − 1 = n.      → (½, −½)
    * round wake: U_∞ΔU δ² constant ⇒ 2m + n = 0;  m − 1 = n.                                      → (⅓, −⅔)
    * plane plume: m = 1; buoyancy flux U Y δ constant and d(U²δ)/dx ∼ g′δ with g′ ∝ Y
      ⇒ n + s = −1 and 2n = s + 1.                                                                 → (1, 0), s = −1
    * round plume: m = 1; U Y δ² constant, d(U²δ²)/dx ∼ g′δ² ⇒ n + s = −2 and 2n − 1 = s.          → (1, −⅓), s = −5/3
    * shear layer: velocity difference fixed ⇒ n = 0; δ′ = const ⇒ m = 1.                          → (1, 0)

    Scalar exponent of jets and wakes from the constant scalar flux (U_conv Y δ^{1 or 2}).
    Returns dict(flow, width m, velocity n, scalar s (None for the shear layer), reynolds = m + n (exponent of the local
    Reynolds number U_scale δ/ν), invariant (text)); with ``return_equations=True`` also ``equations`` — the two linear
    equations as text.  All exponents are ``fractions.Fraction``.
    Validation: V2 the table of exact values; V4 the invariants evaluated with :func:`free_shear_flow` do not vary with x.
    Assumptions: self-similar far field at high Reynolds number; each flow's invariant and growth law as listed.
    """
    Fr = Fraction
    one, zero = Fr(1), Fr(0)
    if flow == "plane_jet":
        m, n = _solve2(one, Fr(2), zero, one, zero, one)
        s, inv = -n - m, "momentum flux per span  rho*int(U^2 dy)  [(12.62)]"
        eqs = ["m + 2n = 0   (U_CL^2 * delta independent of x: momentum flux, (12.65))", "m = 1   (linear growth, delta' = constant)"]
    elif flow == "round_jet":
        m, n = _solve2(Fr(2), Fr(2), zero, one, zero, one)
        s, inv = -n - 2 * m, "momentum flux  rho*int(U^2 2*pi*r dr)"
        eqs = ["2m + 2n = 0   (U_CL^2 * delta^2 independent of x)", "m = 1   (linear growth)"]
    elif flow == "plane_wake":
        m, n = _solve2(one, one, zero, one, -one, one)
        s, inv = -m, "momentum deficit per span  rho*U_inf*int(dU dy) = drag per span"
        eqs = ["m + n = 0   (U_inf * dU_CL * delta independent of x: drag)", "m - n = 1   (d delta/dx ~ dU_CL/U_inf)"]
    elif flow == "round_wake":
        m, n = _solve2(Fr(2), one, zero, one, -one, one)
        s, inv = -2 * m, "momentum deficit  rho*U_inf*int(dU 2*pi*r dr) = drag"
        eqs = ["2m + n = 0   (U_inf * dU_CL * delta^2 independent of x: drag)", "m - n = 1   (d delta/dx ~ dU_CL/U_inf)"]
    elif flow == "plane_plume":
        m = one
        n, s = _solve2(one, one, -m, Fr(2), -one, one)
        inv = "buoyancy flux per span  g*int((rho - rho_local)/rho * U dy)"
        eqs = ["n + s = -1   (U_CL * Y_CL * delta independent of x: buoyancy flux, with m = 1)",
               "2n - s = 1   (d(U_CL^2 delta)/dx ~ Y_CL * delta: buoyancy drives the momentum flux)"]
    elif flow == "round_plume":
        m = one
        n, s = _solve2(one, one, -2 * m, Fr(2), -one, one)
        inv = "buoyancy flux  g*int((rho - rho_local)/rho * U 2*pi*r dr)"
        eqs = ["n + s = -2   (U_CL * Y_CL * delta^2 independent of x: buoyancy flux, with m = 1)",
               "2n - s = 1   (d(U_CL^2 delta^2)/dx ~ Y_CL * delta^2)"]
    elif flow == "shear_layer":
        m, n, s, inv = one, zero, None, "velocity difference U1 - U2 (imposed)"
        eqs = ["n = 0   (the velocity difference is imposed)", "m = 1   (linear growth)"]
    else:
        raise ValueError(f"flow must be one of {FREE_SHEAR_FLOWS}")
    out = dict(flow=flow, width=m, velocity=n, scalar=s, reynolds=m + n, invariant=inv)
    if return_equations:
        out["equations"] = eqs
    return out


def local_reynolds_number_exponent(flow: str) -> Fraction:
    """Exponent of x in the local Reynolds number U_scale δ/ν of a free shear flow: m + n.

    Book: §12.8 (Table 12.1 exponents).  Plane jet ½ (grows), round jet 0 (constant), plane wake 0, round wake −⅓
    (decays: the far wake eventually relaminarises), plumes 1 and ⅔, shear layer 1.
    Returns a fractions.Fraction (m + n).
    Assumptions: as free_shear_exponents.
    Validation: V2 exact fractions for the seven flows.
    """
    return free_shear_exponents(flow)["reynolds"]


def free_shear_profile(flow: str, xi, *, xi_half: float):
    """Mean-profile shape of a free shear flow in its similarity variable (Gaussian fit).

    Book: §12.8, Table 12.1 (nomenclature: F and H are approximately ``exp{−ln 2 · ξ²/ξ_½²}``).  For jets, plumes and
    wakes this is F(ξ) (the wake's is the shape of the *deficit*); for the shear layer the table's
    ``∫_{−∞}^{ξ} F/∫_{−∞}^{∞} F`` is returned, which for a Gaussian F is ``½[1 + erf(sqrt(ln 2) ξ/ξ_½)]`` and runs from
    0 to 1.
    Parameters: flow one of :data:`FREE_SHEAR_FLOWS`; xi similarity variable [–] (y/x or r/x for jets and plumes,
    y/sqrt(θx) or r/(θ²x)^{1/3} for wakes, (y − y_CL)/x for the shear layer); xi_half half-width [–] (**required
    keyword**: an empirical, flow-dependent number).
    Returns the profile value [–] (a float for a float ξ).
    Assumptions: Gaussian shape (a fit to data, not a solution).
    Validation: V1 value ½ at ξ = ξ_½ and 1 at ξ = 0 (jets, wakes, plumes); ½ at ξ = 0 and limits 0, 1 (shear layer).
    """
    if flow not in FREE_SHEAR_FLOWS:
        raise ValueError(f"flow must be one of {FREE_SHEAR_FLOWS}")
    if flow == "shear_layer":
        return _S(0.5 * (1.0 + erf(np.sqrt(_LN2) * _F(xi) / xi_half)))
    return gaussian_profile(xi, xi_half)


def _need(constants: dict | None, *keys):
    if constants is None:
        raise ValueError("free-shear constants are empirical and flow dependent: pass constants={...} "
                         "(FREE_SHEAR_CONSTANTS has no default set)")
    missing = [k for k in keys if k not in constants]
    if missing:
        raise ValueError(f"constants is missing {missing}")
    return [float(constants[k]) for k in keys]


def free_shear_centerline(flow: str, x, *, constants: dict, d: float | None = None, U0: float | None = None,
                          rho_s: float | None = None, rho: float | None = None, Y0: float = 1.0, g: float = G0,
                          theta: float | None = None, U_inf: float | None = None) -> dict:
    """Centreline (amplitude) laws of the self-similar free shear flows in nozzle or wake variables.

    Book: §12.8, Eqs. (12.72)–(12.73) for the plane jet and the forms of Table 12.1 (the numbers of the table are not
    used: ``constants`` must supply C_U and, for the scalar, C_Y):

    * plane_jet:   U_CL = C_U U0 (ρ_s/ρ)^{1/2} (x/d)^{−1/2};  Y_CL = C_Y Y0 (ρ_s/ρ)^{1/2} (x/d)^{−1/2}
    * round_jet:   U_CL = C_U U0 (ρ_s/ρ)^{1/2} (x/d)^{−1};    Y_CL = C_Y Y0 (ρ_s/ρ)^{1/2} (x/d)^{−1}
    * plane_plume: U_CL = C_U (g(ρ − ρ_s)U0 d/ρ)^{1/3};       Y_CL = C_Y Y0 (ρU0²/(g(ρ − ρ_s)d))^{1/3} (x/d)^{−1}
    * round_plume: U_CL = C_U (g(ρ − ρ_s)U0 d/ρ)^{1/3} (x/d)^{−1/3};  Y_CL = C_Y Y0 (ρU0²/(g(ρ − ρ_s)d))^{1/3} (x/d)^{−5/3}
    * plane_wake:  ΔU_CL = C_U U_∞ (x/θ)^{−1/2},  θ = drag/(ρU_∞² span)   [momentum thickness, m]
    * round_wake:  ΔU_CL = C_U U_∞ (x/θ)^{−2/3},  θ² = drag/(ρU_∞²)

    Parameters: x [m]; d slot width or nozzle diameter [m]; U0 exit speed [m/s]; rho_s exit density, rho ambient density
    [kg/m³]; Y0 exit mass fraction; theta wake momentum thickness [m]; U_inf [m/s]; constants dict with "C_U" (and
    "C_Y" if a scalar is wanted).
    Returns dict(U_CL [m/s] (the deficit ΔU_CL for wakes), Y_CL [–] or None, exponents).
    Validation: V7 each amplitude has the units of a velocity; the x-exponents equal :func:`free_shear_exponents`.
    Assumptions: self-similar far field; the amplitude constants are empirical and must be supplied.
    """
    (CU,) = _need(constants, "C_U")
    CY = float(constants["C_Y"]) if "C_Y" in constants else None
    x_ = _F(x)
    ex = free_shear_exponents(flow)
    n = float(ex["velocity"])
    Y = None
    if flow in ("plane_jet", "round_jet"):
        r = np.sqrt(rho_s / rho)
        U = CU * U0 * r * (x_ / d) ** n                               # Eq. (12.72) / Table 12.1
        if CY is not None:
            Y = CY * Y0 * r * (x_ / d) ** float(ex["scalar"])          # Eq. (12.73) / Table 12.1
    elif flow in ("plane_plume", "round_plume"):
        gp = g * (rho - rho_s)
        U = CU * (gp * U0 * d / rho) ** (1.0 / 3.0) * (x_ / d) ** n
        if CY is not None:
            Y = CY * Y0 * (rho * U0 ** 2 / (gp * d)) ** (1.0 / 3.0) * (x_ / d) ** float(ex["scalar"])
    elif flow in ("plane_wake", "round_wake"):
        U = CU * U_inf * (x_ / theta) ** n
    else:
        raise ValueError("free_shear_centerline: use free_shear_flow for the shear layer")
    return dict(U_CL=_S(U), Y_CL=None if Y is None else _S(Y), exponents=ex)


def free_shear_flow(flow: str, x, cross, *, constants: dict, d: float | None = None, U0: float | None = None,
                    rho_s: float | None = None, rho: float | None = None, Y0: float = 1.0, g: float = G0,
                    theta: float | None = None, U_inf: float | None = None, U1: float | None = None,
                    U2: float | None = None, y_CL: float = 0.0) -> dict:
    """Self-similar far-field mean velocity (and scalar) of one of the seven free turbulent shear flows.

    Book: §12.8, Table 12.1 (forms; see :func:`free_shear_centerline`) with Gaussian profiles
    (:func:`free_shear_profile`) in the similarity variable: ξ = y/x or r/x (jets, plumes), (y − y_CL)/x (shear layer),
    y/sqrt(θx) (plane wake), r/(θ²x)^{1/3} (round wake).
    Parameters: flow (one of :data:`FREE_SHEAR_FLOWS`); x [m]; cross = y or r [m]; constants (**required, no default**):
    "C_U", "xi_half_U" (+ "C_Y", "xi_half_Y" for the scalar); for the shear layer "dxi80_coeff", the coefficient of
    ``Δξ_80 = coeff·(U1 − U2)/(½(U1 + U2))`` (the span in ξ of the central 80 % of the velocity difference).
    Returns dict(U [m/s], Y (or None), width [m] (cross-stream distance where the velocity — or the wake deficit — is
    half its centreline value; half of the 10–90 % thickness for the shear layer), centreline [m/s] (U_CL; the deficit
    ΔU_CL for wakes; the mean speed ½(U1 + U2) for the shear layer), xi, Y_CL, exponents; "half_width" and "U_CL" repeat
    width and centreline).  For wakes U = U_∞ − ΔU_CL F; for the shear layer U = U2 + (U1 − U2)·profile.
    Assumptions: far field, high Reynolds number; the density difference matters only through the source terms.
    Validation: V7 exponents and half-width definition; V4 invariants independent of x; V6 book example (private
    constants).  The constants themselves carry the label "qualitative" until a public table is confirmed.
    """
    x_, c_ = _F(x), _F(cross)
    ex = free_shear_exponents(flow)
    if flow == "shear_layer":
        (coef,) = _need(constants, "dxi80_coeff")
        dxi80 = coef * (U1 - U2) / (0.5 * (U1 + U2))
        s = dxi80 / (2.0 * erfinv(0.8))                # ½[1 + erf(ξ/s)] passes 0.1 and 0.9 at ξ = ∓ s erfinv(0.8)
        xi = (c_ - y_CL) / x_
        prof = 0.5 * (1.0 + erf(xi / s))
        w = _S(0.5 * dxi80 * x_)
        return dict(U=_S(U2 + (U1 - U2) * prof), Y=None, width=w, centreline=0.5 * (U1 + U2), xi=_S(xi),
                    U_CL=0.5 * (U1 + U2), Y_CL=None, half_width=w, dxi80=float(dxi80), exponents=ex)
    CUxi, = _need(constants, "xi_half_U")
    cl = free_shear_centerline(flow, x_, constants=constants, d=d, U0=U0, rho_s=rho_s, rho=rho, Y0=Y0, g=g, theta=theta,
                               U_inf=U_inf)
    if flow == "plane_wake":
        scale = np.sqrt(theta * x_)
    elif flow == "round_wake":
        scale = (theta ** 2 * x_) ** (1.0 / 3.0)
    else:
        scale = x_
    xi = c_ / scale
    F = _F(gaussian_profile(xi, CUxi))
    if flow.endswith("wake"):
        U = U_inf - _F(cl["U_CL"]) * F
    else:
        U = _F(cl["U_CL"]) * F
    Y = None
    if cl["Y_CL"] is not None and "xi_half_Y" in constants:
        Y = _F(cl["Y_CL"]) * _F(gaussian_profile(xi, float(constants["xi_half_Y"])))
    w = _S(CUxi * scale)
    return dict(U=_S(U), Y=None if Y is None else _S(Y), width=w, centreline=cl["U_CL"], xi=_S(xi), U_CL=cl["U_CL"],
                Y_CL=cl["Y_CL"], half_width=w, exponents=ex)


def stoichiometric_mass_fraction(fuel_M: float, oxidiser_M: float, moles_O2_per_fuel: float, x_O2: float) -> dict:
    """Fuel mass fraction of a stoichiometric fuel–oxidiser mixture.

    Book: §12.8, Example 12.2 (method only; the numbers there are private).  One mole of fuel needs
    ``moles_O2_per_fuel`` moles of O2, and the oxidiser stream holds a mole (= volume) fraction ``x_O2`` of O2, so
    ``v_fuel = x_O2 v_ox/moles_O2_per_fuel`` with v_fuel + v_ox = 1; then Y = v_f M_f/(v_f M_f + v_ox M_ox).
    Parameters: fuel_M, oxidiser_M molar masses [kg/kmol]; moles_O2_per_fuel [–]; x_O2 [–].
    Returns dict(Y_fuel, v_fuel, v_oxidiser, M_mixture [kg/kmol]).  Perfect gases at equal temperature and pressure.
    Assumptions: perfect gases at equal temperature and pressure (mole fraction = volume fraction); complete single-step reaction.
    Validation: V1 the three fractions sum to one; hydrogen in air as a hand check with public molar masses.
    """
    ratio = x_O2 / moles_O2_per_fuel
    v_ox = 1.0 / (1.0 + ratio)
    v_f = 1.0 - v_ox
    Mm = v_f * fuel_M + v_ox * oxidiser_M
    return dict(Y_fuel=float(v_f * fuel_M / Mm), v_fuel=float(v_f), v_oxidiser=float(v_ox), M_mixture=float(Mm))


def round_jet_distance_for_mass_fraction(Y_target, d: float, rho_s: float, rho: float, Y0: float, C_Y: float):
    """Distance along the axis of a round jet at which the centreline mass fraction has fallen to a target value.

    Book: §12.8, Table 12.1 (round jet, scalar): ``Y_CL = C_Y Y0 (ρ_s/ρ)^{1/2} (x/d)^{−1}`` solved for x (Example 12.2).
    Parameters: Y_target [–]; d nozzle diameter [m]; rho_s, rho [kg/m³]; Y0 exit mass fraction; C_Y (**required**,
    empirical).  Returns x [m].  Scalar-callable.
    Assumptions: self-similar far field of a round jet; the result must be many diameters from the nozzle to be meaningful.
    Validation: V1 inverse of the centreline law of free_shear_centerline to round-off.
    """
    return _S(C_Y * (Y0 / _F(Y_target)) * np.sqrt(rho_s / rho) * d)


# ======================================================================================================================
# §12.10  turbulence modelling: eddy viscosity and mixing length
# ======================================================================================================================
def eddy_viscosity_stress(gradU, nu_T: float, e: float) -> np.ndarray:
    """Turbulent-viscosity hypothesis: ``mean(u_i u_j) = (2/3) ē δ_ij − ν_T (∂U_i/∂x_j + ∂U_j/∂x_i)``.

    Book: §12.10, Eq. (12.94).  Parameters: gradU (d, d) with [i, j] = ∂U_i/∂x_j [1/s]; nu_T eddy viscosity [m²/s]
    (a property of the flow, not of the fluid); e turbulent kinetic energy ē [m²/s²].
    Returns the modelled covariance mean(u_i u_j) (d, d) [m²/s²].  For d = 3 and a divergence-free mean flow its trace is
    2ē; in simple shear mean(uv) = −ν_T dU/dy.
    Assumptions: turbulent-viscosity hypothesis: the anisotropic stress is aligned with the mean strain rate (a model).
    Validation: V1 trace 2e for a divergence-free mean flow in 3-D; simple shear gives uv = -nu_T dU/dy.
    """
    G = _F(gradU)
    return (2.0 / 3.0) * e * np.eye(G.shape[0]) - nu_T * (G + G.T)  # Eq. (12.94)


def eddy_viscosity_from_data(uv, dUdy):
    """Eddy viscosity implied by measured stress and shear: ``ν_T = −mean(uv)/(dU/dy)``.

    Book: §12.10, Eq. (12.94) for simple shear (and (12.99)).  Returns ν_T [m²/s]; NaN where dU/dy = 0 (the hypothesis
    says nothing there — on a jet axis the stress vanishes with the shear).
    Assumptions: simple shear; the hypothesis (12.94) is taken as a definition of nu_T.
    Validation: V1 inverse of eddy_viscosity_stress in simple shear; NaN at zero shear.
    """
    uv_, s = _F(uv), _F(dUdy)
    with np.errstate(divide="ignore", invalid="ignore"):
        return _S(np.where(s != 0, -uv_ / np.where(s != 0, s, 1.0), np.nan))


def gradient_diffusion_flux(grad, K):
    """Gradient-diffusion hypothesis for a turbulent flux: ``mean(u_i φ′) = −K ∂φ̄/∂x_i``.

    Book: §12.10, Eqs. (12.95) (heat: K = κ_T, φ = T) and (12.96) (passive scalar: K = κ_mT, φ = Y).
    Parameters: grad mean gradient [unit of φ/m]; K eddy diffusivity [m²/s] ≥ 0.  Returns the flux [unit of φ · m/s],
    always down the mean gradient.  Scalar-callable.
    Assumptions: gradient-diffusion hypothesis: the flux runs down the mean gradient (a model; it fails for counter-gradient transport).
    Validation: V1 sign opposite to the gradient; linear in K.
    """
    return _S(-_F(K) * _F(grad))  # Eqs. (12.95)–(12.96)


def eddy_diffusivity_estimate(l_T, u_T, c: float = 1.0):
    """Eddy viscosity or diffusivity as a turbulent length times a turbulent velocity: ``ν_T, κ_T, κ_mT ∼ c l_T u_T``.

    Book: §12.10, Eq. (12.98).  l_T [m]; u_T [m/s]; c order-one constant.  Returns [m²/s].  Scalar-callable.
    Assumptions: an order-of-magnitude scaling law (the constant is an argument).
    Validation: V2 units m2/s.
    """
    return _S(c * _F(l_T) * _F(u_T))  # Eq. (12.98)


def mixing_length_stress(dUdy, l_T):
    """Mixing-length estimate of the Reynolds shear stress, ``−mean(uv) = l_T² |dU/dy| dU/dy``.

    Book: §12.10 (unnumbered, with Fig. 12.20): ``−mean(uv) = ν_T dU/dy ∼ l_T u_T dU/dy ∼ l_T² (dU/dy)²``.
    The book's (dU/dy)² is positive whatever the sign of the shear; the stress is ν_T dU/dy with ν_T > 0, so it must
    change sign with dU/dy: the square is coded as |dU/dy| dU/dy (identical for dU/dy > 0).
    Parameters: dUdy [1/s]; l_T mixing length [m].  Returns −mean(uv) [m²/s²], with the sign of dU/dy.
    Assumptions: mixing-length model: one length l_T, velocity scale l_T |dU/dy|.
    Validation: V1 odd in dU/dy; equals mixing_length_eddy_viscosity times dU/dy.
    """
    s = _F(dUdy)
    # DEVIATION: |dU/dy| dU/dy instead of the printed (dU/dy)² — the printed square loses the sign of the stress when
    # the shear is negative (upper half of a channel, upper half of a jet).
    return _S(_F(l_T) ** 2 * np.abs(s) * s)


def mixing_length_eddy_viscosity(dUdy, l_T):
    """Eddy viscosity of the mixing-length model: ``ν_T = l_T² |dU/dy|`` [m²/s] (≥ 0).

    Book: §12.10 (unnumbered relation before (12.100): u_T ∼ l_T dU/dy in (12.98)).  Near a wall l_T = κ y.
    Returns nu_T [m2/s]; a float for float input.
    Assumptions: mixing-length model.
    Validation: V1 non-negative; units m2/s.
    """
    return _S(_F(l_T) ** 2 * np.abs(_F(dUdy)))


def _ml_length(yp, kappa, damping, A_plus):
    l = kappa * yp
    if damping is None:
        return l
    if damping != "van_driest":
        raise ValueError("damping must be None or 'van_driest'")
    return l * (-np.expm1(-yp / A_plus))


def _ml_slope(yp, kappa, damping, A_plus, tau=1.0):
    """Positive root of s + l⁺² s² = τ⁺ for s = dU⁺/dy⁺, in the cancellation-free form."""
    l = _ml_length(yp, kappa, damping, A_plus)
    return 2.0 * tau / (1.0 + np.sqrt(1.0 + 4.0 * l ** 2 * tau))


def _ml_uplus_undamped(yp, kappa):
    s = 2.0 * kappa * _F(yp)
    with np.errstate(divide="ignore", invalid="ignore"):
        tail = np.where(s > 1e-4, (np.sqrt(1.0 + s * s) - 1.0) / np.where(s > 0, s, 1.0), 0.5 * s - s ** 3 / 8.0)
    return (np.arcsinh(s) - tail) / kappa


def mixing_length_intercept(kappa: float, A_plus: float | None = None) -> float:
    """Additive constant B of the logarithmic law that the wall mixing-length model implies.

    Book: §12.10, Eqs. (12.100)–(12.101): the model returns ``U⁺ = (1/κ) ln y⁺ + const``; this is the constant, i.e.
    ``lim_{y⁺→∞} [U⁺ − (1/κ) ln y⁺]``.
    Without damping (l_T = κy all the way to the wall, the book's model) it is exact: ``B = (ln 4κ − 1)/κ`` (≈ −1.2 for
    κ = 0.41) — far below measured values, because the model lets the turbulence act inside the sublayer.
    With van Driest's wall damping, l_T = κy[1 − exp(−y⁺/A⁺)] (van Driest 1956; A⁺ ≈ 26 is his value, ours to pass),
    B = B_undamped + ∫_0^∞ (s_damped − s_undamped) dy⁺ (``quad``; the integrand decays like exp(−y⁺/A⁺)).
    Parameters: kappa (**required**); A_plus (None or 0 → no damping).  Returns B [–].  Scalar-callable.
    Validation: V2 the undamped closed form by sympy; V5 damped value within 5 % of a public log-law intercept.
    Assumptions: constant-stress layer over a smooth wall; l_T = kappa y, optionally with van Driest damping.
    """
    B0 = (np.log(4.0 * kappa) - 1.0) / kappa
    if not A_plus:
        return float(B0)
    f = lambda y: float(_ml_slope(y, kappa, "van_driest", A_plus) - _ml_slope(y, kappa, None, A_plus))  # noqa: E731
    Y = 60.0 * A_plus
    val, _ = quad(f, 0.0, Y, points=[0.1 * A_plus, A_plus, 5.0 * A_plus, 20.0 * A_plus], epsabs=1e-12, epsrel=1e-11, limit=400)
    return float(B0 + val)


def mixing_length_wall_profile(yplus, kappa: float, damping: str | None = None, A_plus: float = 26.0,
                               method: str = "quad") -> dict:
    """Mean velocity near a wall from the mixing-length model with a constant total stress.

    Book: §12.10, Eq. (12.100) and its first integral ``ν dU/dy + κ² y² (dU/dy)² = τ0/ρ``; in wall units
    ``s + l⁺² s² = 1`` with s = dU⁺/dy⁺ and l⁺ = κ y⁺.  Its positive root (ours) is
    ``dU⁺/dy⁺ = 2/(1 + sqrt(1 + 4 l⁺²))`` → 1 in the sublayer (12.82) and → 1/(κ y⁺) outside it (12.101).
    Without damping the integral is closed-form: ``U⁺ = (1/κ)[asinh(2κy⁺) − (sqrt(1 + 4κ²y⁺²) − 1)/(2κy⁺)]``.
    Parameters: yplus ≥ 0 scalar or array; kappa (**required**); damping None (the book's model) or "van_driest"
    (l⁺ = κy⁺[1 − exp(−y⁺/A⁺)]); A_plus damping constant (used only with damping); method "quad" (adaptive quadrature
    of the damped correction, accurate to ~1e-10) or "trapezoid" (cumulative trapezoid on the *given* grid from the
    first point, which must be 0 — second order, for convergence studies).
    Returns dict(yplus, Uplus, slope (= dU⁺/dy⁺ = the viscous share of the stress), lT_plus, nuT_over_nu (= l⁺² s),
    uv_plus (= −mean(uv)⁺ = 1 − s), B (the intercept this model implies, :func:`mixing_length_intercept`)).
    Scalar-callable (a dict of floats for a float y⁺).
    Assumptions: constant-stress layer (τ̄ = τ0, i.e. y ≪ δ), smooth wall, mixing length l_T = κy (optionally damped).
    Validation: V1 the root satisfies the quadratic to round-off; limits; V3 trapezoid converges at order 2.
    """
    yp = _F(yplus)
    s = _ml_slope(yp, kappa, damping, A_plus)
    l = _ml_length(yp, kappa, damping, A_plus)
    if method == "trapezoid":
        ya = np.atleast_1d(yp)
        if ya[0] != 0.0:
            raise ValueError("method='trapezoid' needs a grid that starts at yplus = 0")
        U = cumulative_trapezoid(np.atleast_1d(s), ya, initial=0.0)
    elif method == "quad":
        U = _ml_uplus_undamped(yp, kappa)
        if damping is not None:
            f = lambda q: float(_ml_slope(q, kappa, damping, A_plus) - _ml_slope(q, kappa, None, A_plus))  # noqa: E731
            corr = np.vectorize(lambda Y: quad(f, 0.0, Y, epsabs=1e-12, epsrel=1e-11, limit=200,
                                               points=[p for p in (0.1 * A_plus, A_plus, 5 * A_plus, 20 * A_plus) if p < Y])[0]
                                if Y > 0 else 0.0)(yp)
            U = U + corr
    else:
        raise ValueError("method must be 'quad' or 'trapezoid'")
    return dict(yplus=_S(yp), Uplus=_S(U), slope=_S(s), lT_plus=_S(l), nuT_over_nu=_S(l ** 2 * s), uv_plus=_S(1.0 - s),
                B=mixing_length_intercept(kappa, A_plus if damping else None))


def shear_flow_eddy_viscosity_solve(y, nu_T_fn: Callable, dPdx: float, rho: float, nu: float, bc: tuple = (0.0, 0.0),
                                    tol: float = 1e-10, max_iter: int = 500, relax: float = 0.5) -> dict:
    """Solve the unidirectional mean-flow equation with an eddy viscosity between two walls.

    Book: §12.10, Eq. (12.99) ``0 = −(1/ρ) dP/dx + d/dy([ν + ν_T] dU/dy)``.
    Parameters: y (n,) nodes [m] (any spacing), first and last on the walls; nu_T_fn callable (y_face, dUdy_face) →
    ν_T [m²/s] at the cell faces (a constant, a function of y, or a mixing-length rule that uses the shear); dPdx
    [Pa/m]; rho [kg/m³]; nu [m²/s]; bc wall velocities (U(y[0]), U(y[-1])) [m/s].
    Method (ours — the book prescribes none): conservative second-order differences with the viscosity at the faces,
    tridiagonal solve, Picard iteration on ν_T with under-relaxation ``relax`` until max|ΔU| < tol·max|U|.
    Returns dict(y, U, y_face, nu_T (faces), dUdy (faces), iterations, converged, residual (max change of the last
    iteration)).  ν_T ≡ 0 returns the laminar parabola (ch08).
    Validation: V1 ν_T = 0 → Poiseuille to 1e-10 on any grid; V3 order 2 with a smooth ν_T(y).
    Assumptions: steady, fully developed unidirectional mean flow between two walls; constant density.
    """
    from scipy.linalg import solve_banded

    y = _F(y)
    n = y.size
    hf = np.diff(y)
    yf = 0.5 * (y[1:] + y[:-1])
    U = np.linspace(bc[0], bc[1], n)
    rhs_const = float(dPdx) / rho
    nuT = np.zeros(n - 1)
    it, change, converged = 0, np.inf, False
    for it in range(1, int(max_iter) + 1):
        dU = np.diff(U) / hf
        new = np.maximum(_F(nu_T_fn(yf, dU)) * np.ones(n - 1), 0.0)
        nuT = new if it == 1 else relax * new + (1.0 - relax) * nuT
        a = (nu + nuT) / hf                       # face conductance
        ab = np.zeros((3, n))
        b = np.zeros(n)
        hc = 0.5 * (hf[1:] + hf[:-1])
        ab[1, 1:-1] = -(a[1:] + a[:-1]) / hc
        ab[0, 2:] = a[1:] / hc
        ab[2, :-2] = a[:-1] / hc
        b[1:-1] = rhs_const                        # Eq. (12.99)
        ab[1, 0] = ab[1, -1] = 1.0
        b[0], b[-1] = bc
        Unew = solve_banded((1, 1), ab, b)
        change = float(np.max(np.abs(Unew - U)))
        U = Unew
        if change <= tol * max(float(np.max(np.abs(U))), 1e-300):
            converged = True
            break
    return dict(y=y, U=U, y_face=yf, nu_T=nuT, dUdy=np.diff(U) / hf, iterations=it, converged=converged, residual=change)


def _wall_grid(Re_tau: float, n: int, first: float = 0.1) -> np.ndarray:
    """n nodes from 0 to Re_tau, geometrically stretched from a first spacing of ``first`` wall units."""
    if Re_tau <= first * (n - 1):
        return np.linspace(0.0, Re_tau, n)
    from scipy.optimize import brentq

    g = lambda r: first * (r ** (n - 1) - 1.0) / (r - 1.0) - Re_tau  # noqa: E731
    r_hi = min(2.0, float(np.exp(600.0 / (n - 1))))   # keeps r**(n − 1) below the float64 range for any n
    r = brentq(g, 1.0 + 1e-12, r_hi)
    yp = first * (r ** np.arange(n) - 1.0) / (r - 1.0)
    yp[-1] = Re_tau
    return yp


def _channel_point(yp, Re_tau: float, kappa: float, A_plus: float | None, core_cap: float | None):
    """(τ⁺, l⁺, s = dU⁺/dy⁺) of the mixing-length channel at height(s) y⁺ — the one place the model is written."""
    yp = _F(yp)
    tau = 1.0 - yp / Re_tau                                   # linear total stress, (12.76)–(12.77) with (12.90)
    l = kappa * yp if core_cap is None else np.minimum(kappa * yp, core_cap * Re_tau)
    if A_plus:
        l = l * (-np.expm1(-yp / A_plus))                     # van Driest damping
    s = 2.0 * tau / (1.0 + np.sqrt(1.0 + 4.0 * l ** 2 * tau))  # positive root of s + l⁺² s² = τ⁺   (12.99)–(12.100)
    return tau, l, s


def channel_mixing_length(Re_tau: float, kappa: float, A_plus: float | None = 26.0, n: int = 400,
                          core_cap: float | None = 0.09) -> dict:
    """A whole turbulent channel from the mixing length and the exact linear total stress (the smallest closed model).

    Book: §12.10 (12.99)–(12.100) combined with §12.9 (12.76)–(12.77), ours.  In wall units over the lower half-channel
    (δ⁺ = Re_τ): ``τ⁺ = 1 − y⁺/Re_τ = s + l⁺² s²``, s = dU⁺/dy⁺, so ``s = 2τ⁺/(1 + sqrt(1 + 4 l⁺² τ⁺))``, with
    ``l⁺ = min(κ y⁺, core_cap·Re_τ)·[1 − exp(−y⁺/A⁺)]`` (van Driest damping; the cap keeps the length finite in the
    core, a standard outer-layer choice l/δ ≈ 0.09; both are options: A_plus None/0 → no damping, core_cap None → no cap).
    Parameters: Re_tau = δ u_*/ν; kappa (**required**); A_plus; n grid points (geometric stretching from Δy⁺ = 0.1);
    core_cap.
    Returns dict(yplus, y_over_delta, Uplus, dUdy_plus (viscous stress), uv_plus (= −mean(uv)⁺, Reynolds stress),
    total (τ⁺), production (= uv_plus·dUdy_plus, in u_*⁴/ν), lT_plus, nuT_over_nu, U_bulk_plus, U_cl_plus,
    Re_bulk = U_bulk·h/ν = 2 Re_τ U_bulk⁺ (h = full height), Cf = 2/U_bulk⁺², yplus_peak_production).
    U⁺ by cumulative trapezoid (second order in the grid).
    Production peaks where the viscous and Reynolds stresses are equal (s = τ⁺/2), with value τ⁺²/4 ≈ ¼, near y⁺ ≈ 12.
    A model, not a simulation: "approximate" against DNS.
    Validation: V1 total stress linear to round-off; V3 order 2; V5 U⁺ and C_f against public channel DNS in a stated band.
    Assumptions: fully developed channel, smooth walls; mixing-length closure with optional damping and core cap (a model).
    """
    yp = _wall_grid(float(Re_tau), int(n))
    tau, l, s = _channel_point(yp, Re_tau, kappa, A_plus, core_cap)
    U = cumulative_trapezoid(s, yp, initial=0.0)
    uv = tau - s
    prod = uv * s
    Ub = float(np.trapezoid(U, yp) / Re_tau)
    k = int(np.argmax(prod))
    return dict(yplus=yp, y_over_delta=yp / Re_tau, Uplus=U, dUdy_plus=s, uv_plus=uv, total=tau, production=prod,
                lT_plus=l, nuT_over_nu=l ** 2 * s, U_bulk_plus=Ub, U_cl_plus=float(U[-1]), Re_bulk=2.0 * Re_tau * Ub,
                Cf=2.0 / Ub ** 2, yplus_peak_production=float(yp[k]), Re_tau=float(Re_tau), kappa=float(kappa))


def channel_energy_budget(Re_tau: float, kappa: float, A_plus: float | None = 26.0, n: int = 400,
                          core_cap: float | None = 0.09) -> dict:
    """The two kinetic-energy budgets across a model channel, term by term, in wall units (u_*⁴/ν).

    Book: §12.7, Eqs. (12.46) and (12.47) for U(y), evaluated on :func:`channel_mixing_length`.  With the momentum
    balance 0 = 1/Re_τ + dτ⁺/dy⁺ multiplied by U⁺:
    mean flow — ``pressure_work`` U⁺/Re_τ (gain), ``transport`` d(U⁺τ⁺)/dy⁺ (moves energy toward the wall),
    ``viscous_dissipation`` −(dU⁺/dy⁺)² (loss), ``loss_to_turbulence`` −P⁺ (loss);
    turbulence — ``production`` +P⁺ = −mean(uv)⁺ dU⁺/dy⁺ (the same term, opposite sign), and
    ``dissipation_plus_transport`` = −P⁺, **the residual that closes (12.47), not a model of ε̄** (the mixing-length
    closure carries no information on how that sink splits into dissipation and transport).
    Parameters
    ----------
    Re_tau : friction Reynolds number δ u_*/ν [–] (δ = half-height).   kappa : von Kármán constant [–] (no default).
    A_plus : van Driest damping constant [–] (None or 0 → no damping).   n : grid points over the half-channel.
    core_cap : cap of the mixing length as a fraction of δ (None → no cap).
    Returns
    -------
    The :func:`channel_mixing_length` dict (arrays yplus, Uplus, uv_plus = −mean(uv)⁺, dUdy_plus, production, …) plus

    * flat term arrays, **all as sizes (≥ 0 where the name is a gain or a loss)**: ``pressure_work`` W = U⁺/Re_τ,
      ``mean_dissipation`` s² = (dU⁺/dy⁺)², ``production`` P = −mean(uv)⁺ s, ``turb_sink`` = P (the dissipation +
      transport that must remove it: a residual, **labelled model**), ``transport`` = d(U⁺τ⁺)/dy⁺ (signed), so that
      ``pressure_work + transport − mean_dissipation − production = 0`` (12.46) and ``production − turb_sink = 0``
      (12.47);
    * ``mean`` and ``turbulence`` — the same terms with their signs in the two budgets (dicts of arrays);
    * ``integrals`` over the half-channel: pressure_work (= work), mean_dissipation (= dissipation), production,
      transport; ``identity_residual`` = (work − dissipation − production)/work, which must vanish (transport
      integrates to zero): the function raises if |residual| > 1e-3; ``direct_fraction`` = dissipation/work.
    Assumptions: fully developed channel, smooth walls; mixing-length closure (a model, not a simulation); wall units.
    Validation: V4 the integral identity; V1 mean-flow terms sum to zero pointwise (exactly, by construction of the
    transport term) and agree with :func:`mean_energy_budget` to O(Δy²); V3 order 2 in n.
    """
    ch = channel_mixing_length(Re_tau, kappa, A_plus, n, core_cap)
    yp, U, s, tau, P = ch["yplus"], ch["Uplus"], ch["dUdy_plus"], ch["total"], ch["production"]
    W, tr = U / Re_tau, tau * s - U / Re_tau                               # d(Uτ)/dy = τ U' + U τ' exactly
    mean = dict(pressure_work=W, transport=tr, viscous_dissipation=-s ** 2, loss_to_turbulence=-P)  # Eq. (12.46)
    turb = dict(production=P, dissipation_plus_transport=-P)               # Eq. (12.47)
    integ = dict(work=float(np.trapezoid(W, yp)), dissipation=float(np.trapezoid(s ** 2, yp)),
                 production=float(np.trapezoid(P, yp)), transport=float(np.trapezoid(tr, yp)))
    integ.update(pressure_work=integ["work"], mean_dissipation=integ["dissipation"])
    res = (integ["work"] - integ["dissipation"] - integ["production"]) / integ["work"]
    if abs(res) > 1e-3:
        raise ArithmeticError(f"channel energy identity not closed (relative residual {res:.2e}); refine n")
    ch.update(pressure_work=W, mean_dissipation=s ** 2, turb_sink=P.copy(), transport=tr, mean=mean, turbulence=turb,
              integrals=integ, identity_residual=float(res), direct_fraction=integ["dissipation"] / integ["work"])
    return ch


def channel_energy_budget_at(yplus: float, Re_tau: float, kappa: float, A_plus: float | None,
                             core_cap: float | None = 0.09) -> dict:
    """The two kinetic-energy budgets of the model channel at **one height**, as a dict of floats (wall units).

    Book: §12.7, Eqs. (12.46) ``DĒ/Dt = ∂/∂x_j(−U_j P/ρ0 + 2ν U_i S̄_ij − mean(u_i u_j) U_i) − 2ν S̄_ij S̄_ij +
    mean(u_i u_j) ∂U_i/∂x_j − (g/ρ0) ρ̄ U_3`` and (12.47) for U(y) (production ``−mean(uv) dU/dy``), on the
    mixing-length channel of :func:`channel_mixing_length` (§12.9 (12.76)–(12.77), §12.10 (12.99)–(12.100)).
    The scalar twin of :func:`channel_energy_budget` — same model, same keys — for explainer parity and the notebook's
    worked number.  Everything is closed-form at the height except U⁺, which is the integral of the slope from the
    wall (adaptive quadrature, relative tolerance 1e-10).
    Parameters
    ----------
    yplus : height y u_*/ν [–], 0 ≤ y⁺ ≤ Re_τ.   Re_tau : δ u_*/ν [–].   kappa : von Kármán constant [–].
    A_plus : van Driest constant [–] (None or 0 → no damping).   core_cap : mixing-length cap as a fraction of δ.
    Returns
    -------
    dict of floats, in units of u_*⁴/ν (stresses in τ0, velocity in u_*): yplus, Uplus, total (τ⁺ = 1 − y⁺/Re_τ),
    lT_plus, slope (dU⁺/dy⁺ = viscous stress), uv_plus (−mean(uv)⁺ = τ⁺ − slope), pressure_work (U⁺/Re_τ),
    mean_dissipation (slope²), production ((τ⁺ − slope)·slope), turb_sink (= production; dissipation + transport as
    the residual of (12.47), labelled model), transport (d(U⁺τ⁺)/dy⁺ = τ⁺·slope − U⁺/Re_τ), viscous_transport
    (d(U⁺·slope)/dy⁺; the slope's derivative by a central difference of the closed form), reynolds_transport (the
    rest), dissipation_over_production.  ``pressure_work + transport − mean_dissipation − production = 0`` exactly.
    Assumptions: as :func:`channel_energy_budget`.
    Validation: V1 the budget closes to round-off; equals the interpolated arrays of :func:`channel_energy_budget` to
    O(Δy²) (U⁺ to ~1e-4 at n = 400, the local terms to round-off); production ≤ τ⁺²/4 with equality where the viscous
    and Reynolds stresses are equal.
    """
    y, Re = float(yplus), float(Re_tau)
    if not 0.0 <= y <= Re:
        raise ValueError("yplus must lie in [0, Re_tau]")
    tau, l, s = (float(a) for a in _channel_point(y, Re, kappa, A_plus, core_cap))
    slope = lambda q: float(_channel_point(q, Re, kappa, A_plus, core_cap)[2])  # noqa: E731
    kink = [] if core_cap is None else [core_cap * Re / kappa]
    pts = [p for p in ([0.1 * A_plus, A_plus, 5.0 * A_plus] if A_plus else []) + [1.0, 10.0, 100.0] + kink if 0 < p < y]
    U = quad(slope, 0.0, y, points=sorted(set(pts)) or None, epsabs=1e-13, epsrel=1e-10, limit=400)[0] if y > 0 else 0.0
    P = (tau - s) * s                                            # production −mean(uv)⁺ dU⁺/dy⁺   (12.46), (12.47)
    W = U / Re                                                   # pressure work  −U (dP/dx)/ρ in wall units
    tr = tau * s - W                                             # d(U⁺ τ⁺)/dy⁺
    h = 1e-4 * max(y, 1.0)
    lo, hi = max(y - h, 0.0), min(y + h, Re)
    visc_tr = s * s + U * (slope(hi) - slope(lo)) / (hi - lo)    # d(U⁺ s)/dy⁺
    return dict(yplus=y, Uplus=float(U), total=tau, lT_plus=l, slope=s, uv_plus=tau - s, pressure_work=float(W),
                mean_dissipation=s * s, production=P, turb_sink=P, transport=float(tr), viscous_transport=float(visc_tr),
                reynolds_transport=float(tr - visc_tr),
                dissipation_over_production=s * s / P if P > 0 else float("inf"))


def convective_velocity_scale(L, dT, T, g: float = G0):
    """Free-fall velocity of a buoyant fluctuation across a layer of depth L: ``w ∼ sqrt(g L ΔT/T)``.

    Book: §12.10, Eq. (12.102) ``Dw/Dt ∼ g α T′ ∼ g ΔT/T`` (α = 1/T for a perfect gas) with ``w²/L ∼ Dw/Dt``.
    Parameters: L layer depth [m]; dT temperature difference [K]; T mean absolute temperature [K].  Returns w [m/s].
    Assumptions: perfect gas (alpha = 1/T); free convection without mean shear; an order-of-magnitude scaling law (the constant is an argument).
    Validation: V2 units m/s; a hand value.
    """
    return _S(np.sqrt(g * _F(L) * _F(dT) / _F(T)))  # Eq. (12.102)


def convective_eddy_diffusivity(L, dT, T, g: float = G0):
    """Thermal eddy diffusivity of a convecting layer: ``κ_T ∼ w L`` with w from :func:`convective_velocity_scale`.

    Book: §12.10 (after (12.102), an application of (12.98)).  Returns κ_T [m²/s] — typically 10⁴–10⁵ times molecular.
    Assumptions: as convective_velocity_scale.
    Validation: V2 units m2/s; equals w L.
    """
    return _S(_F(convective_velocity_scale(L, dT, T, g)) * _F(L))


# ======================================================================================================================
# §12.10  one- and two-equation models
# ======================================================================================================================
def one_equation_closure(e, l_T, c: float, C_eps: float, sigma_e: float) -> dict:
    """Ingredients of a one-equation model: velocity scale from ē, modelled dissipation and transport coefficient.

    Book: §12.10 (unnumbered, before (12.103)): ``u_T = c sqrt(ē)``, ``ε̄ = C_ε ē^{3/2}/l_T`` and the gradient-diffusion
    model of the transport terms with coefficient ν_T/σ_e.  (The viscous transport there is printed with u_j where
    (12.47) has u_i — slip #8.)
    Parameters: e [m²/s²]; l_T prescribed length [m]; c, C_eps, sigma_e model constants (**required**).
    Returns dict(u_T [m/s], nu_T = l_T u_T [m²/s] (12.98), eps [m²/s³], transport_diffusivity = ν_T/σ_e [m²/s]).
    Assumptions: a one-equation model with a prescribed length l_T (a closure assumption).
    Validation: V1 arithmetic; units of each entry.
    """
    e_, l_ = _F(e), _F(l_T)
    uT = c * np.sqrt(e_)
    return dict(u_T=_S(uT), nu_T=_S(l_ * uT), eps=_S(C_eps * e_ ** 1.5 / l_), transport_diffusivity=_S(l_ * uT / sigma_e))


def k_epsilon_eddy_viscosity(e, eps, C_mu: float = K_EPSILON_CONSTANTS["C_mu"]):
    """Eddy viscosity of the k–ε model: ``ν_T = C_μ ē²/ε̄``.

    Book: §12.10, Eq. (12.104) (l_T = ē^{3/2}/ε̄, u_T = ē^{1/2} in (12.98)).  e [m²/s²]; eps [m²/s³].  Returns [m²/s].
    Assumptions: standard high-Reynolds-number k-epsilon model.
    Validation: V2 units m2/s; a hand value.
    """
    return _S(C_mu * _F(e) ** 2 / _F(eps))  # Eq. (12.104)


def k_epsilon_length_scale(e, eps):
    """Turbulent length scale of the k–ε model: ``l_T = ē^{3/2}/ε̄`` [m].  Book: §12.10 (before (12.104)).

    Returns l_T [m]; a float for float input.
    Assumptions: standard k-epsilon model.
    Validation: V2 units m.
    """
    return _S(_F(e) ** 1.5 / _F(eps))


def k_epsilon_rhs(e, eps, production, constants: dict = K_EPSILON_CONSTANTS):
    """Source terms of the k–ε model equations (everything except advection and the diffusion terms).

    Book: §12.10, Eqs. (12.103) ``Dē/Dt = diffusion − ε̄ − mean(u_i u_j) ∂U_i/∂x_j`` and (12.105)
    ``Dε̄/Dt = diffusion − C_ε1 (mean(u_i u_j) ∂U_i/∂x_j)(ε̄/ē) − C_ε2 ε̄²/ē``.  With P ≡ −mean(u_i u_j) ∂U_i/∂x_j:
    ``dē/dt = P − ε̄``,  ``dε̄/dt = C_ε1 P ε̄/ē − C_ε2 ε̄²/ē``.
    Parameters: e [m²/s²]; eps [m²/s³]; production P [m²/s³]; constants dict with C_eps1, C_eps2.
    Returns (de_dt [m²/s³], deps_dt [m²/s⁴]).  Note: in equilibrium P = ε̄ the ε-source is (C_ε1 − C_ε2) ε̄²/ē ≠ 0 — it is
    balanced by diffusion in the log layer (:func:`k_epsilon_loglayer_kappa`).
    Assumptions: standard k-epsilon model; advection and diffusion terms are not included.
    Validation: V1 with production = 0 the pair integrates to k_epsilon_decay; signs of each term.
    """
    e_, ep, P = _F(e), _F(eps), _F(production)
    return _S(P - ep), _S(constants["C_eps1"] * P * ep / e_ - constants["C_eps2"] * ep ** 2 / e_)  # (12.103), (12.105)


def k_epsilon_decay(e0: float, eps0: float, t, C_eps2: float = K_EPSILON_CONSTANTS["C_eps2"], method: str = "closed") -> tuple:
    """Decay of homogeneous turbulence without mean shear according to the k–ε model.

    Book: §12.10, Eqs. (12.103), (12.105) with no gradients: ``dē/dt = −ε̄``, ``dε̄/dt = −C_ε2 ε̄²/ē`` (ours: the standard
    calibration case).  Solution: ``ē = e0 (1 + t/t0)^{−n}``, ``ε̄ = eps0 (1 + t/t0)^{−(n+1)}``, ``n = 1/(C_ε2 − 1)``,
    ``t0 = n e0/eps0`` — so C_ε2 is fixed by the measured decay exponent of grid turbulence.
    Parameters: e0 [m²/s²]; eps0 [m²/s³]; t scalar or array [s] ≥ 0; C_eps2; method "closed" or "ivp"
    (``solve_ivp``, RK45, rtol 1e-10 — the independent cross-check).
    Returns (e [m²/s²], eps [m²/s³], n [–], t0 [s]) — ē(t) and ε̄(t) with the shape of t (floats for a float t), the decay
    exponent and the virtual time origin; so ``e, eps, n, t0 = k_epsilon_decay(...)``.  The eddy viscosity and length
    scale along the decay are ``k_epsilon_eddy_viscosity(e, eps)`` and ``k_epsilon_length_scale(e, eps)``.
    Assumptions: homogeneous turbulence, no mean shear, no transport (all gradients zero); high Reynolds number (the
    model has no final viscous period).
    Validation: V1 closed form vs solve_ivp to 1e-8; ε̄ = −dē/dt; V2 substitution of the power law into the pair.
    """
    tt = _F(t)
    n = 1.0 / (C_eps2 - 1.0)
    t0 = n * e0 / eps0
    if method == "closed":
        e = e0 * (1.0 + tt / t0) ** (-n)
        ep = eps0 * (1.0 + tt / t0) ** (-(n + 1.0))
    elif method == "ivp":
        ta = np.atleast_1d(tt)
        sol = solve_ivp(lambda _t, q: [-q[1], -C_eps2 * q[1] ** 2 / q[0]], (0.0, float(ta.max()) or 1.0), [e0, eps0],
                        t_eval=np.sort(ta), rtol=1e-10, atol=1e-14)
        order = np.argsort(np.argsort(ta))
        e, ep = sol.y[0][order].reshape(tt.shape), sol.y[1][order].reshape(tt.shape)
    else:
        raise ValueError("method must be 'closed' or 'ivp'")
    return _S(e), _S(ep), float(n), float(t0)  # solution of (12.103), (12.105) without gradients


def k_epsilon_loglayer_kappa(C_mu: float = K_EPSILON_CONSTANTS["C_mu"], C_eps1: float = K_EPSILON_CONSTANTS["C_eps1"],
                             C_eps2: float = K_EPSILON_CONSTANTS["C_eps2"],
                             sigma_eps: float = K_EPSILON_CONSTANTS["sigma_eps"]) -> float:
    """von Kármán constant implied by the k–ε constants in a constant-stress logarithmic layer.

    Book: §12.10, Eqs. (12.103)–(12.105) applied to the log layer of (12.88) (ours, the standard consistency relation):
    with P = ε̄ = u_*³/(κy), ē = u_*²/sqrt(C_μ) and ν_T = κ u_* y, the ε-equation requires
    ``κ² = sqrt(C_μ) (C_ε2 − C_ε1) σ_ε``.
    Returns κ [–] (≈ 0.433 for the standard constants — the model's own von Kármán constant).
    Validation: V2 sympy substitution of the log-layer forms into (12.105).
    Assumptions: constant-stress logarithmic layer with production = dissipation.
    """
    return float(np.sqrt(np.sqrt(C_mu) * (C_eps2 - C_eps1) * sigma_eps))


def k_epsilon_channel(Re_tau: float, n: int = 200, constants: dict = K_EPSILON_CONSTANTS, wall_yplus: float = 30.0,
                      *, B: float, max_iter: int = 20000, tol: float = 1e-9, relax: float = 0.5) -> dict:
    """Fully developed channel flow with the standard high-Reynolds-number k–ε model and wall functions (optional demo).

    Book: §12.10, the closed set (12.94), (12.99), (12.103)–(12.105) with the wall-function boundary condition the text
    describes (the log law (12.88) applied at the first node).  In wall units over the half-channel, y⁺ from
    ``wall_yplus`` to Re_τ: momentum is integrated once exactly, ``(1 + ν_T⁺) dU⁺/dy⁺ = 1 − y⁺/Re_τ``; ē and ε̄ obey
    ``0 = d/dy((1 + ν_T/σ) dφ/dy) + source``.  Wall node: ē⁺ = 1/sqrt(C_μ), ε̄⁺ = 1/(κ y⁺), U⁺ = ln(y⁺)/κ + B with the
    model's own κ (:func:`k_epsilon_loglayer_kappa`); centreline: zero gradient.
    Parameters: Re_tau (≥ ~300 so that a log region exists); n nodes (uniform in ln y⁺); constants; wall_yplus; B
    additive log-law constant (**required keyword**); max_iter, tol, relax of the iteration.
    Method (ours): second-order conservative differences, sinks treated implicitly (−ε̄ = −(ε̄/ē) ē, −C_ε2 ε̄²/ē =
    −(C_ε2 ε̄/ē) ε̄), tridiagonal solves, under-relaxed Picard iteration.
    Returns dict(yplus, Uplus, e_plus, eps_plus, nuT_over_nu, U_bulk_plus, Cf, kappa_model, iterations, converged,
    residual).  Label: qualitative unless the verifier's grid study and DNS comparison say otherwise.
    Assumptions: fully developed channel; standard k-epsilon model with wall functions at the first node.
    Validation: qualitative (a model with wall functions; grid study and comparison with public DNS are the verifier's).
    """
    from scipy.linalg import solve_banded

    Cmu, C1, C2 = constants["C_mu"], constants["C_eps1"], constants["C_eps2"]
    sk, se = constants["sigma_e"], constants["sigma_eps"]
    kap = k_epsilon_loglayer_kappa(Cmu, C1, C2, se)
    y = np.exp(np.linspace(np.log(wall_yplus), np.log(Re_tau), int(n)))
    hf = np.diff(y)
    hc = np.empty(y.size)
    hc[1:-1] = 0.5 * (hf[1:] + hf[:-1])
    tau = 1.0 - y / Re_tau
    k = np.full(y.size, 1.0 / np.sqrt(Cmu)) * np.maximum(tau, 0.05)
    eps = np.maximum(tau, 0.05) ** 1.5 / (kap * np.minimum(y, 0.3 * Re_tau))
    k_w, eps_w = 1.0 / np.sqrt(Cmu), 1.0 / (kap * y[0])

    def solve(phi, sigma, nuT, src, sink, wall):
        D = 1.0 + 0.5 * (nuT[1:] + nuT[:-1]) / sigma
        a = D / hf
        ab = np.zeros((3, y.size))
        b = np.zeros(y.size)
        ab[1, 1:-1] = (a[1:] + a[:-1]) / hc[1:-1] + sink[1:-1]
        ab[0, 2:] = -a[1:] / hc[1:-1]
        ab[2, :-2] = -a[:-1] / hc[1:-1]
        b[1:-1] = src[1:-1]
        ab[1, 0], b[0] = 1.0, wall
        # centreline: zero gradient (half control volume)
        ab[1, -1] = a[-1] / (0.5 * hf[-1]) + sink[-1]
        ab[2, -2] = -a[-1] / (0.5 * hf[-1])
        b[-1] = src[-1]
        return solve_banded((1, 1), ab, b)

    it, change, converged = 0, np.inf, False
    for it in range(1, int(max_iter) + 1):
        nuT = Cmu * k ** 2 / eps
        dU = tau / (1.0 + nuT)
        P = nuT * dU ** 2
        k_new = np.maximum(solve(k, sk, nuT, P, eps / k, k_w), 1e-12)
        eps_new = np.maximum(solve(eps, se, nuT, C1 * P * eps / k, C2 * eps / k, eps_w), 1e-14)
        change = float(max(np.max(np.abs(k_new - k) / np.max(k)), np.max(np.abs(eps_new - eps) / np.max(eps))))
        k = relax * k_new + (1.0 - relax) * k
        eps = relax * eps_new + (1.0 - relax) * eps
        if change < tol:
            converged = True
            break
    nuT = Cmu * k ** 2 / eps
    dU = tau / (1.0 + nuT)
    U = np.log(y[0]) / kap + B + cumulative_trapezoid(dU, y, initial=0.0)
    # bulk velocity: below the first node the log law is integrated analytically
    inner = y[0] * (np.log(y[0]) / kap + B - 1.0 / kap)
    Ub = float((inner + np.trapezoid(U, y)) / Re_tau)
    return dict(yplus=y, Uplus=U, e_plus=k, eps_plus=eps, nuT_over_nu=nuT, U_bulk_plus=Ub, Cf=2.0 / Ub ** 2,
                kappa_model=kap, iterations=it, converged=converged, residual=change, label="qualitative (model)")


# ======================================================================================================================
# §12.11  turbulence in a stratified medium
# ======================================================================================================================
def flux_richardson(wT, uw, dUdz, alpha: float, g: float = G0):
    """Flux Richardson number ``Rf = (−g α mean(wT′))/(−mean(uw) dU/dz)`` = buoyant destruction / shear production.

    Book: §12.11, Eq. (12.107).  Parameters: wT heat-flux correlation [K m/s] (upward positive; T′ = potential
    temperature fluctuation); uw Reynolds shear correlation mean(uw) [m²/s²] (negative for dU/dz > 0); dUdz [1/s];
    alpha [1/K]; g [m/s²].
    Returns Rf [–]: **negative for an upward heat flux (unstable), positive for a downward one (stable)**.
    Degenerate input: zero shear production gives ±inf (sign of the buoyancy term) or NaN (0/0) — no shear, no Rf.
    Scalar-callable.  Validation: V1 hand value; sign test with both signs of wT; the flipped-sign mutant fails.
    Assumptions: horizontally uniform mean flow U(z); Boussinesq.
    """
    num = -g * alpha * _F(wT) + 0.0   # "+ 0.0" turns −0.0 into 0.0 for a zero heat flux
    den = -_F(uw) * _F(dUdz) + 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        return _S(num / den)  # Eq. (12.107)


def turbulence_regime(Rf, Rf_cr: float = 0.25):
    """Regime of stratified shear turbulence from the flux Richardson number.

    Book: §12.11 (after (12.107)): Rf < 0 — buoyancy adds to shear production ("convective"; large −Rf = convection
    dominates); 0 ≤ Rf < Rf_cr — "shear-driven" turbulence that stratification weakens (Rf = 0 is the neutral case);
    Rf ≥ Rf_cr — turbulence cannot sustain itself ("decaying").  Rf_cr ≈ 0.25 is an **observed** value (less than 1
    because dissipation takes most of the shear production: steady (12.106) without transport gives ε̄ = P(1 − Rf)) —
    not the Ri > ¼ *theorem* of linear stability (ch11 (11.67)).
    Parameters: Rf scalar or array; Rf_cr threshold.  Returns "convective", "shear-driven", "decaying" or "undefined"
    (NaN input) — a str for scalar input, an array of str otherwise.  The boundary values belong to the upper regime
    (Rf = 0 → "shear-driven", Rf = Rf_cr → "decaying").
    Assumptions: Rf_cr is an observed value, an argument.
    Validation: V7 both sides of each boundary and the boundaries themselves; NaN input.
    """
    r = _F(Rf)
    out = np.where(np.isnan(r), "undefined", np.where(r < 0, "convective", np.where(r < Rf_cr, "shear-driven", "decaying")))
    return str(out) if out.ndim == 0 else out


def stratified_tke_budget(z, U, uw, wT, eps, alpha: float, g: float = G0, dUdz=None) -> dict:
    """Terms of the turbulent-kinetic-energy budget of a horizontally uniform stratified shear flow.

    Book: §12.11, Eq. (12.106) ``∂ē/∂t + U ∂ē/∂x = −∂/∂z(mean(pw)/ρ0 + mean(ew)) − mean(uw) ∂U/∂z + g α mean(wT′) − ε̄``.
    (The triple correlation is written mean(ew) here, ½mean(ev) in (12.75), ½mean(u_i² u_j) in (12.47) — slip #10; it is
    ½ mean(u_i² w) throughout.)
    Parameters: z (n,) heights [m]; U (n,) mean wind [m/s]; uw (n,) mean(uw) [m²/s²]; wT (n,) heat-flux correlation
    [K m/s], upward positive (potential temperature); eps (n,) dissipation rate [m²/s³]; alpha [1/K]; g [m/s²];
    dUdz (n,) or scalar mean shear [1/s] — optional: when given it is used instead of differentiating U (then z and U
    may be single heights, e.g. for a worked example at one level; a profile needs at least two heights otherwise).
    Returns dict(shear_production, buoyancy, dissipation (= −eps), residual (= minus their sum: what transport and
    unsteadiness must supply), Rf (12.107), regime) — arrays over the heights (floats and a str for one height).
    Assumptions: horizontally uniform, Boussinesq; the residual stands for transport and unsteadiness.
    Validation: V1 terms against hand values; Rf equals flux_richardson.
    """
    z_, U_, uw_, wT_, eps_ = (np.atleast_1d(_F(a)) for a in (z, U, uw, wT, eps))
    if dUdz is not None:
        dUdz = np.atleast_1d(_F(dUdz)) * np.ones_like(uw_)
    elif z_.size < 2:
        raise ValueError("one height only: pass the mean shear dUdz (a gradient cannot be formed from a single point)")
    else:
        dUdz = np.gradient(U_, z_, edge_order=2) if z_.size > 2 else np.gradient(U_, z_)
    P = -uw_ * dUdz                                         # shear production, third term of (12.106)
    Bq = g * alpha * wT_                                    # buoyant production (destruction if negative)
    Rf = flux_richardson(wT_, uw_, dUdz, alpha, g)
    one = np.ndim(uw) == 0 and np.ndim(wT) == 0 and dUdz.size == 1
    f = (lambda a: float(a[0])) if one else (lambda a: a)
    return dict(shear_production=f(P), buoyancy=f(Bq), dissipation=f(-eps_ * np.ones_like(P)), residual=f(-(P + Bq - eps_)),
                Rf=f(np.atleast_1d(Rf)), regime=turbulence_regime(f(np.atleast_1d(Rf))))  # Eq. (12.106)


def gradient_richardson_thermal(dTdz, dUdz, alpha: float, g: float = G0, Gamma_a: float = 0.0,
                                convention: str = "kundu", tol: float = 1e-12) -> dict:
    """Gradient Richardson number of a thermally stratified shear flow, from the **in-situ** temperature gradient.

    Convention Γ ≡ dT/dz (Kundu's sign; the adiabatic value is NEGATIVE, ≈ −9.8 K/km in dry air; stable when
    dT/dz > Γ_a).  Meteorology convention Γ ≡ −dT/dz: stable when Γ < Γ_d ≈ +9.8 K/km.  Both are reported.

    Book: §12.11, Eq. (12.108) ``Ri ≡ N²/(dU/dz)² = α g (dT̄/dz)/(dU/dz)²``, where — as the start of the section says —
    the adiabatic part has been subtracted, i.e. T̄ is the potential temperature.  With a thermometer's gradient that
    reads ``N² = g α (dT/dz − Γ_a)``: an atmosphere cooling upward at 6.5 K/km is *stable*, an isothermal layer strongly
    so; (12.108) fed the in-situ gradient without Γ_a would call both wrong.
    Parameters: dTdz in-situ temperature gradient [K/m], **Kundu sign** (negative when T falls with height); dUdz
    [1/s]; alpha [1/K] (1/T for a perfect gas); g; Gamma_a adiabatic gradient [K/m], Kundu sign, **negative** for air
    (``core.stratification.adiabatic_lapse_rate()``; the default 0.0 means ``dTdz`` is already a potential-temperature
    gradient, or a liquid with negligible Γ_a); convention of the returned ``text`` ("kundu" or "meteorology");
    tol neutral band of N² [1/s²].
    Returns dict(Ri [–], N2 [1/s²], dthetadz (= dT/dz − Γ_a, the potential-temperature gradient) [K/m],
    verdict_kundu (the criterion with its numbers in Kundu's convention, e.g. "stable ⇔ dT/dz > Γa: …"), verdict_met
    (the same statement in the meteorological convention Γ ≡ −dT/dz), verdict (the bare word
    "stable"/"neutral"/"unstable"), text (= the verdict string of the chosen ``convention``), text_other (the other
    one), convention; "dtheta_dz" repeats dthetadz).  Ri = +inf/−inf/NaN at dU/dz = 0 for N² >, <, = 0 (the ch11
    rule).  Array input returns arrays for Ri, N2, dthetadz and ``None`` for the text fields (one sentence per layer).
    Verdict and texts come from ``core.stratification.lapse_rate_stability`` (ch01) — never re-worded here.
    Assumptions: Boussinesq; dry air or a liquid (no moisture effects); α and Γ_a uniform over the layer.
    Validation: V1 Ri from (in-situ, Γ_a) equals Ri from dθ/dz; an isothermal layer is stable; texts identical after
    ``lapse_rate_convention``; the no-Γ_a mutant fails on the standard atmosphere.
    """
    dT, dU = _F(dTdz), _F(dUdz)
    dth = dT - Gamma_a
    N2 = g * alpha * dth                                   # Eq. (12.108) with the potential-temperature gradient
    from .ch11_instability import gradient_richardson  # reuse of ch11 (11.66), including its zero-shear rule

    N2b, dUb = np.broadcast_arrays(N2, dU)
    Ri = _F(gradient_richardson(0.0, N2=N2b, dUdz=dUb))    # Ri = N²/(dU/dz)²; ±inf or NaN where dU/dz = 0
    Ri = np.where((dUb == 0) & (np.abs(N2b) <= tol), np.nan, Ri)  # neutral band of N² at a shear-free level
    out = dict(Ri=_S(Ri), N2=_S(N2), dthetadz=_S(dth), verdict_kundu=None, verdict_met=None, verdict=None, text=None,
               text_other=None, convention=convention, dtheta_dz=_S(dth))
    if dT.ndim == 0:
        is_kundu = convention.strip().lower() == "kundu"
        k = lapse_rate_stability(float(dT), Gamma_a=float(Gamma_a), convention="kundu")
        m = lapse_rate_stability(float(dT), Gamma_a=float(Gamma_a), convention="meteorology")
        out.update(verdict_kundu=k.text, verdict_met=m.text, verdict=k.verdict, text=k.text if is_kundu else m.text,
                   text_other=m.text if is_kundu else k.text)
    return out


def turbulent_prandtl(nu_T, kappa_T):
    """Turbulent Prandtl number ``Pr_T = ν_T/κ_T`` (≈ 1 neutral; > 1 stable; < 1 unstable).  Book: §12.11, Eq. (12.109).

    Returns Pr_T [-].
    Assumptions: both eddy coefficients defined by (12.94), (12.95).
    Validation: V1 arithmetic.
    """
    return _S(_F(nu_T) / _F(kappa_T))


def flux_from_gradient_richardson(Ri, Pr_T: float = 1.0):
    """Flux Richardson number from the gradient one: ``Rf = Ri/Pr_T``.

    Book: §12.11, Eq. (12.109) ``Ri = (ν_T/κ_T) Rf`` (from (12.94), (12.95) in (12.107)).  So Rf_cr = ¼ corresponds to
    Ri = Pr_T/4, which exceeds ¼ — even 1 — in stable conditions.  Scalar-callable.
    Returns Rf [-].
    Assumptions: eddy-coefficient closures (12.94), (12.95).
    Validation: V1 Ri = Pr_T Rf round trip.
    """
    return _S(_F(Ri) / Pr_T)  # Eq. (12.109)


def monin_obukhov_length(u_star, wT, alpha: float | None = None, T: float | None = None, *, kappa: float, g: float = G0):
    """Monin–Obukhov length ``L_M = −u_*³/(κ α g mean(wT′))``.

    Book: §12.11, Eq. (12.110).  Parameters: u_star [m/s]; wT surface heat-flux correlation [K m/s], **upward
    positive**; alpha [1/K] (None → perfect gas, α = 1/T, T required [K]); kappa von Kármán constant (**required
    keyword**); g.
    Returns L_M [m]: **positive for a downward heat flux (stable), negative for an upward one (unstable)**; ``+inf`` for
    zero heat flux (neutral — the limit from either side is the same profile).  Scalar-callable.
    Meteorology writes L = −u_*³ θ̄/(κ g mean(w′θ′)): the same quantity with α = 1/θ̄.
    Validation: V1 units and sign; the no-minus mutant fails; Rf = z/L_M reproduces (12.107) with the log-law gradient.
    Assumptions: constant-flux surface layer; fluxes evaluated at the surface.
    """
    if alpha is None:
        if T is None:
            raise ValueError("give alpha, or T for a perfect gas (alpha = 1/T)")
        alpha = 1.0 / float(T)
    q = _F(wT)
    with np.errstate(divide="ignore"):
        L = np.where(q != 0, -_F(u_star) ** 3 / (kappa * alpha * g * np.where(q != 0, q, 1.0)), np.inf)  # Eq. (12.110)
    return _S(L)


def monin_obukhov_from_fluxes(tau, H, rho: float, cp: float, T: float, *, kappa: float, g: float = G0,
                              alpha: float | None = None):
    """Monin–Obukhov length from the surface stress and the surface heat flux in engineering units.

    Book: §12.11, Eq. (12.110) with ``u_* = sqrt(τ/ρ)`` (12.81) and ``mean(wT′) = H/(ρ C_p)`` (12.32).
    Parameters: tau surface stress [Pa]; H sensible heat flux [W/m²], upward positive; rho [kg/m³]; cp [J/(kg K)];
    T mean absolute temperature [K] (α = 1/T unless ``alpha`` is given); kappa (**required keyword**).
    Returns L_M [m] (sign as :func:`monin_obukhov_length`).  Scalar-callable.
    Assumptions: constant-flux surface layer; dry air as a perfect gas unless alpha is given.
    Validation: V1 equals monin_obukhov_length with u_* = sqrt(tau/rho), wT = H/(rho cp); 83.0 m for the notebook's night case.
    """
    return monin_obukhov_length(np.sqrt(_F(tau) / rho), _F(H) / (rho * cp), alpha=alpha, T=T, kappa=kappa, g=g)


def flux_richardson_surface_layer(z, L_M):
    """Flux Richardson number in the logarithmic surface layer: ``Rf = z/L_M``.

    Book: §12.11, Eq. (12.111) (from (12.107) with dU/dz = u_*/(κz) and −mean(uw) = u_*²).  So |L_M| is the height at
    which buoyancy matters as much as shear.  Returns Rf [–] (0 for L_M = ±inf).  Scalar-callable.
    Assumptions: logarithmic surface layer with constant stress and heat flux (the neutral shear u_*/(kappa z) is used).
    Validation: V1 equals flux_richardson with the log-law gradient and -uw = u_*^2.
    """
    return _S(_F(z) / _F(L_M))  # Eq. (12.111)


def gradient_richardson_surface_layer(z, L_M, Pr_T: float = 1.0, beta: float = 5.0, consistent: bool = False):
    """Gradient Richardson number in the surface layer.

    Book: §12.11, Eqs. (12.109) and (12.111): ``Ri = Pr_T Rf = Pr_T z/L_M`` (``consistent=False``, the book's two
    equations as they stand, both built on the neutral logarithmic gradient).
    ``consistent=True`` (ours) uses the shear of the log-linear profile itself, φ_m = 1 + β z/L_M: then
    ``Rf = (z/L_M)/φ_m`` and ``Ri = Pr_T (z/L_M)/φ_m`` — which saturates at Pr_T/β (= 0.2 for β = 5) however stable the
    layer; NaN where φ_m ≤ 0.
    Parameters: z [m]; L_M [m]; Pr_T turbulent Prandtl number; beta.  Returns Ri [–].  Scalar-callable.
    Assumptions: as flux_richardson_surface_layer; constant turbulent Prandtl number.
    Validation: V1 Pr_T times flux_richardson_surface_layer; the consistent form saturates at Pr_T/beta.
    """
    zeta = _F(z) / _F(L_M)
    if not consistent:
        return _S(Pr_T * zeta)
    return _S(Pr_T * zeta / _F(WT.dimensionless_shear(zeta, beta)))


def surface_layer_regime(z, L_M, crossover: float = 1.0):
    """Which part of the stratified surface layer a height lies in.

    Book: §12.11 (after (12.111) ``Rf = z/L_M``) and Fig. 12.21: L_M is the height at which buoyant destruction (or
    production) and shear production are of the same order.  For z ≪ |L_M| — stable or unstable — stratification is
    slight, the profile is logarithmic and the turbulence is mechanically forced: the "forced convection" region.  For
    z ≫ −L_M in an unstable layer buoyancy generates the turbulence: "free convection".  For z ≫ L_M in a stable layer
    stratification dominates.
    Parameters: z height [m] > 0; L_M Monin–Obukhov length [m] (sign of (12.110); ±inf = neutral); crossover = the
    |z/L_M| at which "≪" is taken to end (nominal 1, an argument: the book gives an order of magnitude, not a number).
    Returns one of four strings: "neutral" (z/L_M = 0, i.e. no heat flux), "forced convection" (0 < |z/L_M| < crossover,
    either sign), "free convection" (z/L_M ≤ −crossover), "stable" (z/L_M ≥ +crossover).  A str for scalar input, an
    array of str otherwise.
    Assumptions: constant-flux surface layer; the logarithmic-layer estimate (12.111) of Rf.
    Validation: V7 both sides of each boundary, both signs of L_M, and L_M = ±inf.
    """
    zeta = _F(z) / _F(L_M)                                     # Eq. (12.111)
    out = np.where(zeta == 0.0, "neutral",
                   np.where(zeta <= -crossover, "free convection",
                            np.where(zeta >= crossover, "stable", "forced convection")))
    return str(out) if out.ndim == 0 else out


def surface_layer_state(u_star: float, H: float, T: float, z0: float, z: float, rho: float, cp: float, kappa: float,
                        beta: float = 5.0, Pr_T: float = 1.0, convention: str = "kundu", Rf_cr: float = 0.25,
                        g: float = G0, unstable: str = "log_linear") -> dict:
    """Everything the stratified surface layer says at one height, from the surface stress and heat flux.

    Convention Γ ≡ dT/dz (Kundu's sign, Γ_a = −g/C_p NEGATIVE); the meteorological Γ ≡ −dT/dz verdict text is returned
    alongside (``text_other``).

    Book: §12.11, Eqs. (12.107)–(12.111) and the log-linear wind profile; one entry point for the notebook and the
    explainer.
    Parameters: u_star [m/s]; H surface sensible heat flux [W/m²], **upward positive**; T mean absolute temperature [K]
    (α = 1/T, perfect gas); z0 roughness length [m]; z height [m]; rho [kg/m³]; cp [J/(kg K)]; kappa von Kármán
    constant (no default); beta coefficient of the log-linear profile; Pr_T turbulent Prandtl number; convention of
    ``text``; Rf_cr; g; unstable: "log_linear" (the book's formula on both sides — its linear correction makes the shear
    vanish at z/L_M = −1/β; beyond that ``valid`` is False, φ_m and the gradients are NaN and U is the formula's value,
    which no longer increases with height) or "businger_dyer"
    (φ_m = (1 − 16ζ)^{−1/4} for ζ < 0, usable for any unstable ζ; see ``WT.dimensionless_shear``).
    Returns dict:
      wT = H/(ρ C_p) [K m/s];  L_M (12.110) [m];  zeta = z/L_M;  Rf = z/L_M (12.111);  Ri = Pr_T·Rf (12.109);
      phi_m = 1 + β ζ;  Rf_profile = ζ/φ_m and Ri_profile = Pr_T ζ/φ_m (ours: with the shear of the log-linear profile);
      U (log-linear) and U_neutral (12.93) [m/s];  regime (:func:`turbulence_regime` of Rf);  layer
      (:func:`surface_layer_regime`);  z_crit = Rf_cr·L_M [m] (height where Rf reaches Rf_cr; NaN unless stable);
      dtheta_dz = −wT Pr_T φ_m/(κ u_* z) [K/m] (flux–gradient relation (12.95) with κ_T = κ u_* z/(φ_m Pr_T));
      dTdz = dθ/dz + Γ_a [K/m] (in-situ, Kundu sign);  Gamma_a = −g/C_p;
      verdict_kundu and verdict_met — the static-stability criterion of the layer at z with its numbers, in Kundu's
      convention (Γ ≡ dT/dz) and in the meteorological one (Γ ≡ −dT/dz), both always returned, from
      ``lapse_rate_stability``;  verdict (the bare word), text (= the one of ``convention``), text_other;  valid (bool:
      φ_m > 0).  "z_Rf_cr" repeats z_crit.
      The same φ_m is used for heat and momentum apart from the factor Pr_T (a simplification, ours).
    Assumptions: horizontally uniform, steady, constant-flux layer; dry air; moderate |z/L_M|.
    Validation: V1 composition of the tested pieces; neutral limit H = 0; sign flips with H.
    """
    wT = H / (rho * cp)
    L = float(monin_obukhov_length(u_star, wT, T=T, kappa=kappa, g=g))
    zeta = 0.0 if np.isinf(L) else z / L
    phi = float(WT.dimensionless_shear(zeta, beta, unstable))
    valid = bool(phi == phi)
    Rf = zeta
    Ga = float(adiabatic_lapse_rate(cp=cp, g=g))
    dth = -wT * Pr_T * phi / (kappa * u_star * z)
    dTdz = dth + Ga
    is_kundu = convention.strip().lower() == "kundu"
    if np.isfinite(dTdz):
        k = lapse_rate_stability(dTdz, Gamma_a=Ga, convention="kundu", cp=cp, g=g)
        m = lapse_rate_stability(dTdz, Gamma_a=Ga, convention="meteorology", cp=cp, g=g)
        verdict, v_kundu, v_met = k.verdict, k.text, m.text
    else:
        verdict = v_kundu = v_met = None
    z_crit = Rf_cr * L if (L > 0 and np.isfinite(L)) else float("nan")
    return dict(wT=wT, L_M=L, zeta=zeta, Rf=Rf, Ri=Pr_T * Rf, phi_m=phi,
                Rf_profile=zeta / phi if phi == phi else float("nan"), Ri_profile=Pr_T * zeta / phi if phi == phi else float("nan"),
                U=float(WT.surface_layer_wind(z, u_star, z0, L, kappa=kappa, beta=beta, unstable=unstable)),
                U_neutral=float(WT.rough_wall_log_law(z, u_star, z0, kappa=kappa)),
                regime=turbulence_regime(Rf, Rf_cr), layer=surface_layer_regime(z, L), z_crit=z_crit,
                dtheta_dz=dth, dTdz=dTdz, Gamma_a=Ga, verdict_kundu=v_kundu, verdict_met=v_met, verdict=verdict,
                text=v_kundu if is_kundu else v_met, text_other=v_met if is_kundu else v_kundu,
                convention=convention, valid=valid, z_Rf_cr=z_crit)


def temperature_variance_budget(z, T_mean, wT, eps_T, transport=None) -> dict:
    """Terms of the budget of temperature variance ½ mean(T′²) in a horizontally uniform flow.

    Book: §12.11, Eq. (12.112): production ``−mean(wT′) dT̄/dz`` (positive when the flux runs down the mean gradient),
    turbulent and molecular transport ``−∂/∂z(½ mean(T′² w) − κ ∂(½ mean(T′²))/∂z)``, dissipation
    ``ε̄_T = κ mean((∂T′/∂x_j)²)``.  (The page prints the molecular transport without the ½ — slip #15; the derived form
    is used, see :func:`temperature_variance_sympy`.)
    Parameters: z (n,) [m]; T_mean (n,) [K] (potential temperature); wT (n,) [K m/s]; eps_T (n,) [K²/s]; transport (n,)
    [K²/s] or None (→ the residual that closes the steady budget, inferred).
    Returns dict(production, dissipation (= −eps_T), transport, residual, transport_inferred).
    Assumptions: horizontally uniform flow; steady unless the residual is read as the unsteady term.
    Validation: V1 production positive for a down-gradient flux; inferred transport closes the budget to round-off.
    """
    z_, Tm, q, eT = _F(z), _F(T_mean), _F(wT), _F(eps_T)
    # DEVIATION (documentation only): the molecular transport belonging to this budget is κ ∂(½ mean(T′²))/∂z, not the
    # printed κ ∂mean(T′²)/∂z (slip #15); it is not computed here (it is lumped into ``transport``).
    prod = -q * np.gradient(Tm, z_, edge_order=2)  # Eq. (12.112)
    inferred = transport is None
    tr = -(prod - eT) if inferred else _F(transport)
    return dict(production=prod, dissipation=-eT, transport=tr, residual=prod - eT + tr, transport_inferred=inferred)


def batchelor_scale(nu, kappa_th, eps):
    """Batchelor scale ``η_T = η (κ/ν)^{1/2}``: where temperature gradients are smeared out when ν/κ ≫ 1.

    Book: §12.11 (before (12.114)).  nu [m²/s]; kappa_th thermal (or scalar) diffusivity [m²/s]; eps [m²/s³].
    Returns η_T [m] (< η for a Prandtl number above one, e.g. water).  Scalar-callable.
    Assumptions: Prandtl (Schmidt) number above one; locally isotropic small scales.
    Validation: V1 eta_T/eta = (kappa_th/nu)^(1/2).
    """
    return _S(_F(kolmogorov_scales(nu, eps)[0]) * np.sqrt(_F(kappa_th) / _F(nu)))


def scalar_spectrum_exponents() -> dict:
    """Exponents of the inertial-range temperature spectrum S_T ∝ ε̄_T^a ε̄^b K^c by exact dimensional analysis.

    Book: §12.11, Eq. (12.113): S_T [K² m] from ε̄_T [K²/s], ε̄ [m²/s³], K [1/m] ⇒ (a, b, c) = (1, −⅓, −5/3).
    Returns dict(eps_T, eps, K) of ``Fraction``.
    Assumptions: only eps_T, eps and K matter in the inertial-convective range.
    Validation: V2 exact fractions (1, -1/3, -5/3).
    """
    from .core.dimensional import solve_exponents

    grp = solve_exponents("S_T", ("eps_T", "eps", "K"), {"S_T": "K**2*m", "eps_T": "K**2/s", "eps": "m**2/s**3", "K": "1/m"})
    return dict(eps_T=-grp["eps_T"], eps=-grp["eps"], K=-grp["K"])


def scalar_spectrum(K, eps: float, eps_T: float, nu: float, kappa_th: float, C_T: float = 1.0, cutoff: bool = False):
    """Spectrum of temperature (scalar) fluctuations: inertial −5/3 range, then Batchelor's −1 range when ν/κ ≫ 1.

    Book: §12.11, Eqs. (12.113) ``S_T ∝ ε̄_T ε̄^{−1/3} K^{−5/3}`` for 2π/L ≪ K ≪ 2π/η and (12.114) ``S_T ∝ K^{−1}`` for
    2π/η ≪ K ≪ 2π/η_T; normalisation mean(T′²) = ∫_0^∞ S_T dK.
    Parameters: K [rad/m] > 0; eps [m²/s³]; eps_T [K²/s]; nu, kappa_th [m²/s]; C_T order-one constant (a proportionality
    in the book: 1 by default, an argument); cutoff: multiply by exp(−(K η_T)²) to end the −1 range at the Batchelor
    scale (ours, for plotting).
    Returns S_T [K² m].  The two power laws are joined at K = 1/η, where they are equal for the natural amplitude of
    the second, ``C_T ε̄_T (ν/ε̄)^{1/2} K^{−1}`` (ours: the strain rate (ε̄/ν)^{1/2} is the only time scale left).  For
    ν/κ ≤ 1 only the −5/3 branch is returned (the book treats ν/κ ≫ 1 only).  Scalar-callable.
    Validation: V2 exponents (1, −⅓, −5/3) (:func:`scalar_spectrum_exponents`); continuity at the junction; slopes.
    Assumptions: locally isotropic turbulence at high Reynolds and Peclet numbers; passive scalar.
    """
    K_ = _F(K)
    eta = (nu ** 3 / eps) ** 0.25
    inertial = C_T * eps_T * eps ** (-1.0 / 3.0) * K_ ** (-5.0 / 3.0)   # Eq. (12.113)
    if nu / kappa_th > 1.0:
        visc = C_T * eps_T * np.sqrt(nu / eps) / K_                     # Eq. (12.114)
        S = np.where(K_ * eta < 1.0, inertial, visc)
    else:
        S = inertial
    if cutoff:
        S = S * np.exp(-(K_ * eta * np.sqrt(kappa_th / nu)) ** 2)
    return _S(S)


# ======================================================================================================================
# §12.12  Taylor's theory of turbulent dispersion
# ======================================================================================================================
def langevin_particles(n: int, t, u_rms: float, Lambda_t: float, seed: int = 0, dim: int = 1) -> tuple:
    """Particles leaving a point source in stationary homogeneous turbulence with an exponential velocity memory.

    Book: §12.12 (Fig. 12.24: paths X(t) from the origin, zero mean velocity).  Model (ours): each velocity component is
    an Ornstein–Uhlenbeck process (Langevin equation), whose Lagrangian autocorrelation is ``r(τ) = exp(−|τ|/Λ_t)`` — so
    Λ_t is exactly the integral time scale of (12.122).
    Parameters: n particles; t (nt,) times [s] starting at 0, increasing; u_rms [m/s] (each component); Lambda_t [s];
    seed of ``default_rng`` (default 0); dim number of independent components.
    Returns (X, u): positions (n, nt) — or (n, nt, dim) for dim > 1 — [m] and velocities of the same shape [m/s], at
    the times ``t`` supplied (particle index first), so ``X, u = langevin_particles(...)``.
    Assumptions: stationary, homogeneous turbulence with zero mean velocity; Gaussian velocities; independent
    components; a model of the Lagrangian velocity (kinematic — no flow field is solved).
    Method: the **exact joint update** of (u, X) over each step h (no time-step error): with ρ = e^{−h/Λ},
    u′ = ρu + ξ_u, X′ = X + Λ(1 − ρ)u + ξ_x, Var ξ_u = σ²(1 − ρ²), Var ξ_x = σ²Λ²(2h/Λ − 3 + 4ρ − ρ²),
    Cov(ξ_u, ξ_x) = σ²Λ(1 − ρ)².
    Validation: V1 mean(X²) within 5 s.e. of :func:`taylor_dispersion_exponential` at every output time for any step;
    V4 d mean(X²)/dt = 2 mean(Xu) (12.116).
    """
    t = _F(t)
    rng = np.random.default_rng(seed)
    shape = (int(n), int(dim))
    s2, Lam = float(u_rms) ** 2, float(Lambda_t)
    u = u_rms * rng.standard_normal(shape)
    X = np.zeros(shape)
    Xs, us = np.empty((t.size,) + shape), np.empty((t.size,) + shape)
    Xs[0], us[0] = X, u
    for k in range(t.size - 1):
        hh = t[k + 1] - t[k]
        rho = np.exp(-hh / Lam)
        var_u = s2 * (-np.expm1(-2.0 * hh / Lam))
        x = hh / Lam
        var_x = s2 * Lam ** 2 * ((2.0 * x - 3.0 + 4.0 * rho - rho * rho) if x > 1e-3 else (2.0 * x ** 3 / 3.0 - x ** 4 / 2.0))
        cov = s2 * Lam * (1.0 - rho) ** 2
        g1, g2 = rng.standard_normal(shape), rng.standard_normal(shape)
        a = np.sqrt(var_u)
        b = cov / a if a > 0 else 0.0
        c = np.sqrt(max(var_x - b * b, 0.0))
        X = X + Lam * (1.0 - rho) * u + b * g1 + c * g2
        u = rho * u + a * g1
        Xs[k + 1], us[k + 1] = X, u
    Xo, uo = np.moveaxis(Xs, 0, 1), np.moveaxis(us, 0, 1)
    if dim == 1:
        Xo, uo = Xo[..., 0], uo[..., 0]
    return Xo, uo


def dispersion_rate_from_particles(t, X, u) -> tuple:
    """Both sides of ``d mean(X²)/dt = 2 mean(X u)`` evaluated on an ensemble of particle paths.

    Book: §12.12, Eqs. (12.115) ``dX_α/dt = u_α`` and (12.116) ``d mean(X_α²)/dt = 2 mean(X_α u_α)`` (ensemble average
    over particles; no sum over the Greek index).
    Parameters: t (nt,) times [s]; X, u (n, nt) one component of position [m] and velocity [m/s], particles on axis 0.
    Returns (lhs, rhs): lhs = d mean(X²)/dt by second-order differences of the sampled mean(X²) [m²/s], rhs =
    2 mean(X u) [m²/s] — two (nt,) arrays; half of rhs is the eddy diffusivity D_T of (12.127).
    The identity is exact for every path; the two columns differ by the differencing error of the sampled mean(X²),
    which for N particles is not smooth in time (its noise, of relative size ~ N^{−1/2}, has the roughness of the
    velocity) — expect agreement at the per-cent level for N ~ 10⁴, improving with N.
    Assumptions: all particles released at X = 0 at t[0]; the average is over the particles supplied.
    Validation: V4 the two sides agree within sampling error for :func:`langevin_particles`; V1 rhs/2 within 5
    standard errors of :func:`eddy_diffusivity_exponential`.
    """
    t_, X_, u_ = _F(t), _F(X), _F(u)
    X2 = np.mean(X_ ** 2, axis=0)
    rhs = 2.0 * np.mean(X_ * u_, axis=0)                       # Eq. (12.116), right side
    return np.gradient(X2, t_, edge_order=2), rhs              # Eq. (12.116), left side (with (12.115))


def _quad_vec(f: Callable, t) -> np.ndarray:
    ta = np.atleast_1d(_F(t)).astype(float)
    out = np.array([quad(f, 0.0, ti, epsabs=1e-13, epsrel=1e-11, limit=400)[0] if ti > 0 else 0.0 for ti in ta.ravel()])
    return out.reshape(ta.shape)


def taylor_dispersion_rate(t, r_fn: Callable, u2: float):
    """Rate of growth of the mean-square displacement: ``d mean(X²)/dt = 2 u2 ∫_0^t r(τ) dτ``.

    Book: §12.12, Eq. (12.117).  (The page writes r(t′ − t) under the first integral; it equals r(t − t′) because r is
    even — the substitution τ = t − t′ gives the second form.)
    Parameters: t [s] ≥ 0; r_fn callable τ → Lagrangian autocorrelation coefficient; u2 velocity variance [m²/s²].
    Returns [m²/s].  ``quad``.  Scalar-callable.
    Assumptions: stationary homogeneous turbulence, zero mean velocity; r is the Lagrangian autocorrelation.
    Validation: V1 exponential r gives 2 u2 Lambda (1 - exp(-t/Lambda)) to 1e-10.
    """
    f = lambda tau: float(r_fn(tau))  # noqa: E731
    return _S((2.0 * u2 * _quad_vec(f, t)).reshape(np.shape(t)))  # Eq. (12.117)


def taylor_dispersion(t, r_fn: Callable, u2: float, form: str = "single"):
    """Mean-square displacement of a particle in stationary homogeneous turbulence (Taylor 1921).

    Book: §12.12, Eq. (12.119) ``mean(X²)(t) = 2 u2 t ∫_0^t (1 − τ/t) r(τ) dτ`` (``form="single"``) and Eq. (12.118)
    ``= 2 u2 ∫_0^t dt′ ∫_0^{t′} r(τ) dτ`` (``form="double"``); the two are related by an integration by parts.
    Parameters: t [s] ≥ 0 scalar or array; r_fn callable τ → r; u2 [m²/s²].
    Returns mean(X²) [m²].  Limits: u2 t² for t ≪ Λ_t (12.120), 2 u2 Λ_t t for t ≫ Λ_t (12.122).
    Validation: V1 the two forms agree to 1e-10; exponential and Gaussian closed forms; both limits.
    Assumptions: stationary homogeneous turbulence, zero mean velocity; r is the Lagrangian autocorrelation.
    """
    ta = np.atleast_1d(_F(t)).astype(float)
    out = np.zeros(ta.shape)
    for idx, ti in np.ndenumerate(ta):
        if ti <= 0:
            continue
        if form == "single":
            val = quad(lambda tau: (ti - tau) * float(r_fn(tau)), 0.0, ti, epsabs=1e-14, epsrel=1e-11, limit=400)[0]  # t(1 − τ/t) r
        elif form == "double":
            inner = lambda tp: quad(lambda tau: float(r_fn(tau)), 0.0, tp, epsabs=1e-14, epsrel=1e-11, limit=200)[0]  # noqa: E731
            val = quad(inner, 0.0, ti, epsabs=1e-13, epsrel=1e-10, limit=200)[0]
        else:
            raise ValueError("form must be 'single' or 'double'")
        out[idx] = 2.0 * u2 * val  # Eq. (12.119) / (12.118)
    return _S(out.reshape(np.shape(t)))


def taylor_dispersion_exponential(t, u2: float, Lambda_t: float):
    """Taylor's formula in closed form for ``r = exp(−τ/Λ_t)``: ``mean(X²) = 2 u2 Λ_t² [t/Λ_t − 1 + exp(−t/Λ_t)]``.

    Book: §12.12, Eq. (12.119) evaluated for an exponential correlation (ours).  Parameters: t [s] ≥ 0; u2 [m²/s²];
    Lambda_t [s].  Returns mean(X²) [m²].  Scalar-callable.
    Limits: u2 t² (12.120) and 2 u2 Λ_t t − 2 u2 Λ_t² (the diffusive law (12.122) with a constant offset).
    Evaluated with ``expm1`` and, for t/Λ_t < 1e-3, the series x²/2 − x³/6 + x⁴/24 (no cancellation).
    Assumptions: as taylor_dispersion with r = exp(-tau/Lambda_t).
    Validation: V1 equals taylor_dispersion (quad) to 1e-10; both limits; Langevin particles within 5 standard errors.
    """
    x = _F(t) / Lambda_t
    br = np.where(x < 1e-3, x ** 2 / 2.0 - x ** 3 / 6.0 + x ** 4 / 24.0, x + np.expm1(-x))
    return _S(2.0 * u2 * Lambda_t ** 2 * br)


def taylor_dispersion_gaussian(t, u2: float, t_c: float):
    """Taylor's formula in closed form for ``r = exp(−τ²/t_c²)`` (integral scale Λ_t = √π t_c/2).

    Book: §12.12, Eq. (12.119) (the Gaussian case is the one Exercise 12.38 asks about; closed form ours):
    ``mean(X²) = 2 u2 [ t (√π t_c/2) erf(t/t_c) − (t_c²/2)(1 − exp(−t²/t_c²)) ]``.
    Parameters: t [s]; u2 [m²/s²]; t_c [s].  Returns [m²].  Scalar-callable.
    Assumptions: as taylor_dispersion with r = exp(-tau^2/t_c^2).
    Validation: V1 equals taylor_dispersion (quad) to 1e-10; limits u2 t^2 and 2 u2 Lambda_t t.
    """
    tt = _F(t)
    return _S(2.0 * u2 * (tt * 0.5 * np.sqrt(np.pi) * t_c * erf(tt / t_c) + 0.5 * t_c ** 2 * np.expm1(-(tt / t_c) ** 2)))


def dispersion_local_slope(t, Lambda_t: float):
    """Local logarithmic slope ``d ln mean(X²)/d ln t`` for the exponential correlation: 2 (ballistic) → 1 (diffusive).

    Book: §12.12, Eqs. (12.120), (12.122) (the two limits); in between ``x(1 − e^{−x})/(x − 1 + e^{−x})``, x = t/Λ_t.
    Scalar-callable.  Equals 1.3 at x ≈ 3.2 and 2 − x/3 for small x.
    Returns the slope [-]; a float for float input.
    Assumptions: exponential Lagrangian correlation.
    Validation: V1 equals the numerical d ln/d ln of taylor_dispersion_exponential; limits 2 and 1.
    """
    x = _F(t) / Lambda_t
    with np.errstate(divide="ignore", invalid="ignore"):
        full = x * (-np.expm1(-x)) / np.where(x < 1e-3, 1.0, x + np.expm1(-x))
    return _S(np.where(x < 1e-3, 2.0 - x / 3.0, full))


def dispersion_regime(t, Lambda_t: float, lo: float = 0.3, hi: float = 3.0):
    """Name of the dispersion regime at time t: "ballistic" (t ≤ lo·Λ_t), "transition", "diffusive" (t ≥ hi·Λ_t).

    Book: §12.12, Eqs. (12.121) ``X_rms = u_rms t`` for t ≪ Λ_t and (12.123) ``X_rms = u_rms sqrt(2 Λ_t t)`` for t ≫ Λ_t.
    The thresholds are nominal (arguments).  str for scalar input, array of str otherwise.
    Returns 'ballistic', 'transition' or 'diffusive' (str; an array of str for array input).
    Assumptions: nominal thresholds (0.3 and 3 integral times), inclusive: t = 0.3 Lambda_t is still 'ballistic', t = 3 Lambda_t already 'diffusive'.
    Validation: V7 both sides of each threshold and the thresholds themselves.
    """
    x = _F(t) / Lambda_t
    out = np.where(x <= lo, "ballistic", np.where(x >= hi, "diffusive", "transition"))
    return str(out) if out.ndim == 0 else out


def random_walk(n_steps: int, n_walkers: int, L: float = 1.0, dim: int = 2, persistence: float = 0.0,
                seed: int = 0) -> np.ndarray:
    """Random walk of fixed step length L in a random direction: after n uncorrelated steps ``R_rms = L sqrt(n)``.

    Book: §12.12, Eqs. (12.124) ``mean(R_n²) = mean(R_{n−1}²) + L² + 2 mean(R_{n−1}·L)`` (the last term vanishes for
    uncorrelated directions) and (12.125).
    Parameters: n_steps; n_walkers; L step length [m]; dim 1, 2 or 3; persistence p in (−1, 1) — ours, the bridge to
    Taylor's theory: with probability |p| a step repeats (p > 0) or reverses (p < 0) the previous direction, otherwise
    it is drawn afresh, so successive directions have correlation p^k and mean(R_n²) → n L² (1 + p)/(1 − p);
    seed of ``default_rng`` (default 0).
    Returns the positions R, an array (n_walkers, n_steps + 1, dim) [m] with R[:, 0] = 0; the mean-square distance
    after n steps is ``(R**2).sum(axis=2).mean(axis=0)[n]``, to be compared with n L² (12.125) — or, with persistence,
    with the large-n law n L² (1 + p)/(1 − p).
    Assumptions: steps of equal length; directions uniform on the circle/sphere (±1 in one dimension).
    Validation: V1 mean(R_n²) = n L² within 5 s.e. in 1, 2 and 3 dimensions; persistence law at large n.
    """
    rng = np.random.default_rng(seed)
    nw, ns = int(n_walkers), int(n_steps)

    def directions(size):
        if dim == 1:
            return rng.choice([-1.0, 1.0], size=(size, 1))
        v = rng.standard_normal((size, dim))
        return v / np.linalg.norm(v, axis=1, keepdims=True)

    steps = np.empty((nw, ns, dim))
    prev = directions(nw)
    for k in range(ns):
        fresh = directions(nw)
        if persistence != 0.0 and k > 0:
            keep = rng.random(nw) < abs(persistence)
            fresh[keep] = np.sign(persistence) * prev[keep]
        steps[:, k] = fresh
        prev = fresh
    return np.concatenate([np.zeros((nw, 1, dim)), np.cumsum(L * steps, axis=1)], axis=1)  # R_n = R_{n−1} + L  (12.124)


def smoke_plume_width(x, U: float, w_rms: float, Lambda_t: float):
    """Root-mean-square half-width of a time-averaged smoke plume in a uniform wind: Taylor's formula with t = x/U.

    Book: §12.12, Fig. 12.27 with (12.121), (12.123): ``Z_rms ∝ x`` close to the source (x ≪ U Λ_t) and ``∝ sqrt(x)``
    far away — a parabola with a pointed vertex.  (The caption of the printed figure names the two regimes the other
    way round — slip #13; the text, the equations and the labels on the figure agree with what is coded.)
    Parameters: x distance downwind [m]; U wind speed [m/s]; w_rms cross-wind velocity rms [m/s]; Lambda_t Lagrangian
    integral time [s].  Returns Z_rms [m] (exponential correlation, :func:`taylor_dispersion_exponential`).
    Scalar-callable.  Assumes u_rms ≪ U (x plays the role of time).
    Validation: V1 local slope d ln Z/d ln x = 1 near the source and 1/2 far away (the printed caption has them swapped, slip #13).
    """
    return _S(np.sqrt(_F(taylor_dispersion_exponential(_F(x) / U, w_rms ** 2, Lambda_t))))


def plume_concentration(x, z, Q: float, U: float, sigma_z):
    """Gaussian cross-section of a plume from a line source: ``c = Q/(U sqrt(2π) σ_z) exp(−z²/(2σ_z²))``.

    Book: §12.12 (Fig. 12.27; the Gaussian shape is ours — the standard closure of the time-averaged plume).
    Parameters: x [m] (unused except for broadcasting; σ_z carries the x-dependence); z cross-wind coordinate [m];
    Q source strength per unit span [kg/(m s)]; U wind speed [m/s]; sigma_z = Z_rms(x) [m].
    Returns concentration [kg/m³]; ∫c dz = Q/U at every x.
    Assumptions: steady line source in a uniform wind; Gaussian cross-section; no ground reflection.
    Validation: V4 the cross-wind integral equals Q/U at every x.
    """
    s = _F(sigma_z) + 0.0 * _F(x)
    return _S(Q / (U * np.sqrt(2.0 * np.pi) * s) * np.exp(-0.5 * (_F(z) / s) ** 2))


def diffusivity_from_variance(t, variance):
    """A diffusivity from the growth of a variance: ``D = ½ dσ²/dt``.

    Book: §12.12, Eq. (12.126) ``ν = ½ dσ²/dt`` (a Gaussian that spreads by constant diffusivity has σ² = 2νt; note
    that ch03's vortex core radius uses 4νt).  Parameters: t (n,) [s]; variance (n,) [m²].  Returns (n,) [m²/s]
    (second-order differences).
    Assumptions: variance sampled smoothly in time (second-order differences).
    Validation: V1 sigma^2 = 2 D t returns D to round-off.
    """
    return 0.5 * np.gradient(_F(variance), _F(t), edge_order=2)  # Eq. (12.126)


def eddy_diffusivity_taylor(t, r_fn: Callable, u2: float):
    """Effective (eddy) diffusivity of turbulent dispersion: ``D_T ≡ ½ d mean(X²)/dt = u2 ∫_0^t r(τ) dτ``.

    Book: §12.12, Eq. (12.127).  It **grows with time** (12.128) before it saturates (12.129): turbulent spreading from
    a source is not diffusion with a bigger constant.  Returns D_T [m²/s].  Scalar-callable.
    Assumptions: as taylor_dispersion_rate.
    Validation: V1 half of taylor_dispersion_rate; tends to u2 Lambda_t.
    """
    return _S(0.5 * _F(taylor_dispersion_rate(t, r_fn, u2)))  # Eq. (12.127)


def eddy_diffusivity_exponential(t, u2: float, Lambda_t: float):
    """Eddy diffusivity for an exponential correlation: ``D_T = u2 Λ_t (1 − exp(−t/Λ_t))``.

    Book: §12.12, Eq. (12.127) with r = exp(−τ/Λ_t) (ours).  Limits: u2 t (12.128) and u2 Λ_t (12.129).
    Scalar-callable.
    Returns D_T [m2/s]; a float for float input.
    Assumptions: exponential Lagrangian correlation.
    Validation: V1 equals eddy_diffusivity_taylor (quad) to 1e-10; limits u2 t and u2 Lambda_t.
    """
    return _S(-u2 * Lambda_t * np.expm1(-_F(t) / Lambda_t))


def eddy_diffusivity_asymptote(t, u2: float, Lambda_t: float, which: str = "long", printed: bool = False):
    """The two limiting forms of the eddy diffusivity, each returned only where its stated condition holds.

    Book: §12.12, Eqs. (12.128) ``D_T ≅ u2 t for t ≪ Λ_t`` and (12.129) ``D_T ≅ u2 Λ_t for t ≫ Λ_t``.
    **The page prints the condition of (12.129) as t ≪ Λ_t as well (slip #12)**; the constant value is the long-time
    limit (it follows from (12.122)).  ``printed=True`` applies the constant under the printed condition so that a test
    can show it fails (at t = 0.01 Λ_t it is 100 times the exact value).
    Parameters: t [s]; u2 [m²/s²]; Lambda_t [s]; which "short" or "long".  The conditions ≪, ≫ are coded as t < Λ_t,
    t > Λ_t.  Returns D_T [m²/s] where the condition holds, NaN elsewhere.
    Assumptions: the conditions 'much less/greater than' are coded as t < Lambda_t and t > Lambda_t.
    Validation: V7 short form within 1 % of the exact D_T for t < 0.02 Lambda_t, long form for t > 5 Lambda_t; printed=True is 100 times too large at t = 0.01 Lambda_t.
    """
    tt = _F(t)
    if which == "short":
        return _S(np.where(tt < Lambda_t, u2 * tt, np.nan))            # Eq. (12.128)
    if which != "long":
        raise ValueError("which must be 'short' or 'long'")
    # DEVIATION: condition t ≫ Λ_t, not the printed t ≪ Λ_t (slip #12) — the constant value is the long-time limit
    # that follows from (12.122); printed=True applies it under the printed condition so a test can show it fails.
    cond = (tt < Lambda_t) if printed else (tt > Lambda_t)
    return _S(np.where(cond, u2 * Lambda_t, np.nan))                   # Eq. (12.129), condition corrected


# ======================================================================================================================
# symbolic engines: Reynolds averaging with an explicit averaging operator (sympy)
# ======================================================================================================================
_FLUCT_NAMES = frozenset({"u1", "u2", "u3", "p", "Tp", "Yp"})
_AVG_CLASS = None


def _avg_class():
    """The symbolic averaging wrapper ``Avg(...)`` (a sympy Function whose derivative moves inside: (12.6), (12.8))."""
    global _AVG_CLASS
    if _AVG_CLASS is None:
        import sympy as sp

        class Avg(sp.Function):
            nargs = 1

            def _eval_derivative(self, s):  # ∂/∂s mean(a) = mean(∂a/∂s): averaging commutes with differentiation
                return _average(sp.diff(self.args[0], s))

            def _latex(self, printer):
                return r"\overline{%s}" % printer._print(self.args[0])

        _AVG_CLASS = Avg
    return _AVG_CLASS


def _fluct_degree(f) -> int:
    """How many fluctuation factors a product factor carries (0 = deterministic: means, constants, averages)."""
    import sympy as sp
    from sympy.core.function import AppliedUndef

    Avg = _avg_class()
    if isinstance(f, Avg) or f.is_number or f.is_Symbol:
        return 0
    if isinstance(f, AppliedUndef):
        return 1 if f.func.__name__ in _FLUCT_NAMES else 0
    if isinstance(f, sp.Derivative):
        return _fluct_degree(f.expr)
    if isinstance(f, sp.Pow):
        d = _fluct_degree(f.base)
        if d and not (f.exp.is_Integer and f.exp > 0):
            raise ValueError("non-polynomial dependence on a fluctuation")
        return d * int(f.exp) if d else 0
    if not any(a.func.__name__ in _FLUCT_NAMES for a in f.atoms(AppliedUndef)):
        return 0
    raise ValueError(f"cannot classify factor {f!r}; expand the expression first")


def _average(expr):
    """Ensemble average of a polynomial expression in means and fluctuations.

    Rules (book §12.3 and §12.5): linear (12.4), (12.5); the mean of a mean quantity is itself; a term with exactly
    one fluctuation factor averages to zero (12.26); a product of two or more fluctuation factors is kept as
    ``Avg(product)``; derivatives pass through the average (12.6), (12.8) (handled by ``Avg._eval_derivative``).
    """
    import sympy as sp

    Avg = _avg_class()
    expr = sp.expand(sp.sympify(expr).doit())
    out = sp.Integer(0)
    for term in sp.Add.make_args(expr):
        coeff, fl, deg = sp.Integer(1), sp.Integer(1), 0
        for fac in sp.Mul.make_args(term):
            d = _fluct_degree(fac)
            if d == 0:
                coeff *= fac
            else:
                fl *= fac
                deg += d
        if deg == 0:
            out += coeff
        elif deg >= 2:
            out += coeff * Avg(fl)
    return out


def _renormalize(expr):
    import sympy as sp

    Avg = _avg_class()
    return sp.expand(expr.replace(lambda e: isinstance(e, Avg), lambda e: _average(e.args[0])))


def _use_continuity(expr, vel, x):
    """Eliminate ∂w/∂x3 (and every derivative of it) with ∂u_i/∂x_i = 0, inside and outside averages."""
    import sympy as sp

    w, z = vel[2], x[2]

    def cond(e):
        return isinstance(e, sp.Derivative) and e.expr == w and z in e.variables

    def repl(e):
        rest = list(e.variables)
        rest.remove(z)
        return -sum(sp.Derivative(vel[i], x[i], *rest) for i in range(2))

    return _renormalize(expr.replace(cond, repl))


def _sym_fields():
    import sympy as sp

    x = sp.symbols("x1 x2 x3", real=True)
    t = sp.Symbol("t", real=True)
    args = (*x, t)
    mk = lambda name: sp.Function(name)(*args)  # noqa: E731
    c = dict(zip(("nu", "rho0", "g", "alpha", "T0", "kappa_th", "kappa_m"),
                 sp.symbols("nu rho_0 g alpha T_0 kappa kappa_m", positive=True)))
    return dict(sp=sp, x=x, t=t, U=[mk(f"U{i}") for i in (1, 2, 3)], u=[mk(f"u{i}") for i in (1, 2, 3)], P=mk("P"),
                p=mk("p"), Tb=mk("Tbar"), Tp=mk("Tp"), Yb=mk("Ybar"), Yp=mk("Yp"), **c)


def _zero(expr, F) -> bool:
    """True if expr vanishes once both continuity equations (12.27), (12.28) are used."""
    e = _use_continuity(_use_continuity(F["sp"].expand(expr), F["u"], F["x"]), F["U"], F["x"])
    return F["sp"].simplify(e) == 0


def _instantaneous_momentum(F):
    """Left minus right of the Boussinesq momentum equation in flux form, ũ = U + u etc. (the form of (12.29))."""
    sp, x, t = F["sp"], F["x"], F["t"]
    ut = [F["U"][i] + F["u"][i] for i in range(3)]
    pt, Tt = F["P"] + F["p"], F["Tb"] + F["Tp"]
    eqs = []
    for i in range(3):
        e = sp.diff(ut[i], t) + sum(sp.diff(ut[j] * ut[i], x[j]) for j in range(3)) + sp.diff(pt, x[i]) / F["rho0"] \
            - F["nu"] * sum(sp.diff(ut[i], x[j], 2) for j in range(3))
        if i == 2:
            e += F["g"] * (1 - F["alpha"] * (Tt - F["T0"]))
        eqs.append(sp.expand(e))
    return eqs


def rans_sympy() -> dict:
    """Derive the Reynolds-averaged equations symbolically and compare them with the book's forms, term by term.

    Book: §12.5, Eqs. (12.27) mean continuity, (12.28) fluctuation continuity, (12.29) → (12.30) mean momentum with the
    mean stress tensor, (12.31) mean temperature, (12.34) mean passive scalar (constant mixture density).
    Method: every field is written mean + fluctuation (12.24); an explicit averaging operator (linear; kills single
    fluctuations; keeps products as ``Avg(...)``; commutes with derivatives) is applied to the Boussinesq set in flux
    form.  Three space dimensions, all components.
    Returns dict(mean_continuity, mean_momentum (list of 3: left − right of the averaged equation as derived),
    fluctuation_momentum (list of 3: total − mean), mean_temperature, mean_scalar, steps (list of (label, expression)
    for printing), checks: dict(eq_12_27, eq_12_28, eq_12_30, eq_12_31, eq_12_34) — each True when the derived equation
    minus the book's form simplifies to 0).
    Validation: V2 — this *is* the symbolic evidence; the product-rule mutant ("mean of a product = product of means")
    removes every Avg term and fails eq_12_30.
    Assumptions: Boussinesq fluid; the averaging operator obeys (12.4)-(12.9).
    """
    F = _sym_fields()
    sp, x, t, U, u = F["sp"], F["x"], F["t"], F["U"], F["u"]
    cont_total = sum(sp.diff(U[i] + u[i], x[i]) for i in range(3))
    cont_mean = _average(cont_total)                                   # Eq. (12.27)
    cont_fluct = sp.expand(cont_total.doit() - cont_mean)              # Eq. (12.28)
    E = _instantaneous_momentum(F)                                     # Eq. (12.29)
    M = [_average(e) for e in E]
    fluct = [sp.expand(E[i].doit() - M[i]) for i in range(3)]
    uu = lambda i, j: _average(u[i] * u[j])  # noqa: E731
    mu = F["rho0"] * F["nu"]
    book = []
    for i in range(3):
        stress_div = sum(sp.diff(-F["P"] * sp.KroneckerDelta(i, j) + mu * (sp.diff(U[i], x[j]) + sp.diff(U[j], x[i]))
                                 - F["rho0"] * uu(i, j), x[j]) for j in range(3))
        b = sp.diff(U[i], t) + sum(U[j] * sp.diff(U[i], x[j]) for j in range(3)) - stress_div / F["rho0"]
        if i == 2:
            b += F["g"] * (1 - F["alpha"] * (F["Tb"] - F["T0"]))
        book.append(b)                                                 # Eq. (12.30), left − right
    Tt = F["Tb"] + F["Tp"]
    heat = sp.diff(Tt, t) + sum(sp.diff((U[j] + u[j]) * Tt, x[j]) for j in range(3)) \
        - F["kappa_th"] * sum(sp.diff(Tt, x[j], 2) for j in range(3))
    heat_mean = _average(heat)
    heat_book = sp.diff(F["Tb"], t) + sum(U[j] * sp.diff(F["Tb"], x[j]) for j in range(3)) \
        + sum(sp.diff(_average(u[j] * F["Tp"]), x[j]) for j in range(3)) \
        - F["kappa_th"] * sum(sp.diff(F["Tb"], x[j], 2) for j in range(3))          # Eq. (12.31)
    Yt = F["Yb"] + F["Yp"]
    scal = sp.diff(Yt, t) + sum(sp.diff((U[j] + u[j]) * Yt, x[j]) for j in range(3)) \
        - F["kappa_m"] * sum(sp.diff(Yt, x[j], 2) for j in range(3))                 # Eq. (12.33), constant ρ_m
    scal_mean = _average(scal)
    scal_book = sp.diff(F["Yb"], t) + sum(U[j] * sp.diff(F["Yb"], x[j]) for j in range(3)) \
        - sum(sp.diff(F["kappa_m"] * sp.diff(F["Yb"], x[j]) - _average(u[j] * F["Yp"]), x[j]) for j in range(3))  # (12.34)
    checks = dict(eq_12_27=sp.simplify(cont_mean - sum(sp.diff(U[i], x[i]) for i in range(3))) == 0,
                  eq_12_28=sp.simplify(cont_fluct - sum(sp.diff(u[i], x[i]) for i in range(3))) == 0,
                  eq_12_30=all(_zero(M[i] - book[i], F) for i in range(3)),
                  eq_12_31=_zero(heat_mean - heat_book, F), eq_12_34=_zero(scal_mean - scal_book, F))
    steps = [("mean continuity (12.27)", cont_mean), ("fluctuation continuity (12.28)", cont_fluct),
             ("averaged x1-momentum, flux form (before (12.30))", M[0]), ("mean temperature (12.31)", heat_mean),
             ("mean scalar (12.34)", scal_mean)]
    return dict(mean_continuity=cont_mean, fluctuation_continuity=cont_fluct, mean_momentum=M, fluctuation_momentum=fluct,
                mean_temperature=heat_mean, mean_scalar=scal_mean, steps=steps, checks=checks, fields=F)


def reynolds_stress_budget_sympy(components=((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))) -> dict:
    """Derive the transport equation of mean(u_i u_j) from the fluctuation equation and compare with the book's.

    Book: §12.5, Eq. (12.35) (stated; derivation left to Exercise 12.16).  Method: fluctuation equation = total − mean
    (from :func:`rans_sympy`); form ``mean(u_j·eq_i + u_i·eq_j)``; compare with (12.35) written out, using
    ∂u_k/∂x_k = 0 (12.28) and ∂U_k/∂x_k = 0 (12.27).
    Parameters: components — index pairs (0-based) to derive (default: all six).
    Returns dict(budgets {(i, j): left − right as derived}, checks {(i, j): bool}, all_ok, trace_half (the ē-equation
    residual, ½ Σ of the diagonal budgets — only if all three diagonal components were requested)).
    Validation: V2 residual 0 for every component.
    Assumptions: Boussinesq fluid; both continuity equations (12.27), (12.28) are used.
    """
    R = rans_sympy()
    F = R["fields"]
    sp, x, t, U, u = F["sp"], F["x"], F["t"], F["U"], F["u"]
    fl = R["fluctuation_momentum"]
    A = lambda e: _average(e)  # noqa: E731
    budgets, checks = {}, {}
    for (i, j) in components:
        derived = _average(u[j] * fl[i] + u[i] * fl[j])
        uiuj = A(u[i] * u[j])
        book = (sp.diff(uiuj, t) + sum(U[k] * sp.diff(uiuj, x[k]) for k in range(3))
                + sum(sp.diff(A(u[i] * u[j] * u[k]), x[k]) for k in range(3))
                + sum(A(u[i] * u[k]) * sp.diff(U[j], x[k]) + A(u[j] * u[k]) * sp.diff(U[i], x[k]) for k in range(3))
                + (A(u[i] * sp.diff(F["p"], x[j])) + A(u[j] * sp.diff(F["p"], x[i]))) / F["rho0"]
                + 2 * F["nu"] * sum(A(sp.diff(u[i], x[k]) * sp.diff(u[j], x[k])) for k in range(3))
                - F["nu"] * sum(sp.diff(uiuj, x[k], 2) for k in range(3))
                - F["g"] * F["alpha"] * (A(u[j] * F["Tp"]) * (1 if i == 2 else 0) + A(u[i] * F["Tp"]) * (1 if j == 2 else 0)))  # (12.35)
        budgets[(i, j)] = derived
        checks[(i, j)] = _zero(derived - book, F)
    out = dict(budgets=budgets, checks=checks, all_ok=all(checks.values()), fields=F)
    if all((k, k) in budgets for k in range(3)):
        out["trace_half"] = sp.expand(sum(budgets[(k, k)] for k in range(3)) / 2)
    return out


def tke_budget_sympy() -> dict:
    """Derive the turbulent-kinetic-energy budget as half the trace of the Reynolds-stress budget.

    Book: §12.7, Eq. (12.47): ``Dē/Dt = ∂/∂x_j(−mean(p u_j)/ρ0 + 2ν mean(u_i S′_ij) − ½ mean(u_i² u_j)) − 2ν mean(S′_ij
    S′_ij) − mean(u_i u_j) ∂U_i/∂x_j + g α mean(u_3 T′)``, ē = ½ mean(u_i²).
    Returns dict(derived (left − right from the trace of (12.35)), book (left − right of (12.47)), check (bool),
    viscous_identity (bool: ν∇²ē − ν mean(∂u_i/∂x_k ∂u_i/∂x_k) = 2ν ∂_j mean(u_i S′_ij) − 2ν mean(S′_ij S′_ij), the step
    that needs (12.28))).
    Validation: V2 residual 0.
    Assumptions: as reynolds_stress_budget_sympy.
    """
    Rb = reynolds_stress_budget_sympy(((0, 0), (1, 1), (2, 2)))
    F = Rb["fields"]
    sp, x, t, U, u = F["sp"], F["x"], F["t"], F["U"], F["u"]
    A = lambda e: _average(e)  # noqa: E731
    Sp = [[(sp.diff(u[i], x[j]) + sp.diff(u[j], x[i])) / 2 for j in range(3)] for i in range(3)]
    e = sum(A(u[i] * u[i]) for i in range(3)) / 2
    transport = sum(sp.diff(-A(F["p"] * u[j]) / F["rho0"] + 2 * F["nu"] * sum(A(u[i] * Sp[i][j]) for i in range(3))
                            - sum(A(u[i] * u[i] * u[j]) for i in range(3)) / 2, x[j]) for j in range(3))
    diss = 2 * F["nu"] * sum(A(Sp[i][j] * Sp[i][j]) for i in range(3) for j in range(3))
    prod = -sum(A(u[i] * u[j]) * sp.diff(U[i], x[j]) for i in range(3) for j in range(3))
    buoy = F["g"] * F["alpha"] * A(u[2] * F["Tp"])
    book = sp.diff(e, t) + sum(U[j] * sp.diff(e, x[j]) for j in range(3)) - (transport - diss + prod + buoy)  # (12.47)
    visc_l = F["nu"] * sum(sp.diff(e, x[k], 2) for k in range(3)) \
        - F["nu"] * sum(A(sp.diff(u[i], x[k]) ** 2) for i in range(3) for k in range(3))
    visc_r = 2 * F["nu"] * sum(sp.diff(A(u[i] * Sp[i][j]), x[j]) for i in range(3) for j in range(3)) - diss
    return dict(derived=Rb["trace_half"], book=book, check=_zero(Rb["trace_half"] - book, F),
                viscous_identity=_zero(visc_l - visc_r, F), production=prod, dissipation=diss, fields=F)


def mean_energy_budget_sympy() -> dict:
    """Derive the kinetic-energy budget of the mean flow by multiplying the mean momentum equation by U_i.

    Book: §12.7, Eq. (12.46): ``DĒ/Dt = ∂/∂x_j(−U_j P/ρ0 + 2ν U_i S̄_ij − mean(u_i u_j) U_i) − 2ν S̄_ij S̄_ij +
    mean(u_i u_j) ∂U_i/∂x_j − (g/ρ0) ρ̄ U_3``, Ē = ½ U_i², with ρ̄ = ρ0[1 − α(T̄ − T0)] (Exercise 12.15).
    Returns dict(derived = Σ U_i × (12.30), book = left − right of (12.46), check (bool), exchange (the term
    mean(u_i u_j) ∂U_i/∂x_j that reappears with the opposite sign in (12.47))).
    Validation: V2 residual 0; the exchange term plus the production of :func:`tke_budget_sympy` is 0.
    Assumptions: Boussinesq fluid; mean continuity (12.27).
    """
    F = _sym_fields()
    sp, x, t, U, u = F["sp"], F["x"], F["t"], F["U"], F["u"]
    A = lambda e: _average(e)  # noqa: E731
    mu = F["rho0"] * F["nu"]
    S = [[(sp.diff(U[i], x[j]) + sp.diff(U[j], x[i])) / 2 for j in range(3)] for i in range(3)]
    mom = []
    for i in range(3):
        stress_div = sum(sp.diff(-F["P"] * sp.KroneckerDelta(i, j) + 2 * mu * S[i][j] - F["rho0"] * A(u[i] * u[j]), x[j])
                         for j in range(3))
        b = sp.diff(U[i], t) + sum(U[j] * sp.diff(U[i], x[j]) for j in range(3)) - stress_div / F["rho0"]
        if i == 2:
            b += F["g"] * (1 - F["alpha"] * (F["Tb"] - F["T0"]))
        mom.append(b)                                                                         # Eq. (12.30)
    derived = sp.expand(sum(U[i] * mom[i] for i in range(3)))
    Ebar = sum(U[i] ** 2 for i in range(3)) / 2
    rho_bar = F["rho0"] * (1 - F["alpha"] * (F["Tb"] - F["T0"]))
    transport = sum(sp.diff(-U[j] * F["P"] / F["rho0"] + 2 * F["nu"] * sum(U[i] * S[i][j] for i in range(3))
                            - sum(A(u[i] * u[j]) * U[i] for i in range(3)), x[j]) for j in range(3))
    diss = 2 * F["nu"] * sum(S[i][j] ** 2 for i in range(3) for j in range(3))
    exchange = sum(A(u[i] * u[j]) * sp.diff(U[i], x[j]) for i in range(3) for j in range(3))
    book = sp.diff(Ebar, t) + sum(U[j] * sp.diff(Ebar, x[j]) for j in range(3)) \
        - (transport - diss + exchange - F["g"] / F["rho0"] * rho_bar * U[2])                 # Eq. (12.46)
    return dict(derived=derived, book=book, check=_zero(derived - book, F), exchange=exchange, fields=F)


def temperature_variance_sympy() -> dict:
    """Derive the budget of temperature variance and test the printed molecular-transport term.

    Book: §12.11, Eq. (12.112).  Method: T′-equation = total − mean; multiply by T′; average.  General result (ours):
    ``D(½ mean(T′²))/Dt = −mean(u_j T′) ∂T̄/∂x_j − ∂/∂x_j(½ mean(T′² u_j)) + κ ∇²(½ mean(T′²)) − κ mean((∂T′/∂x_j)²)``,
    which for a horizontally uniform flow is (12.112) **with κ ∂(½ mean(T′²))/∂z inside the transport bracket — the page
    prints κ ∂mean(T′²)/∂z, without the ½ (slip #15)**.
    Returns dict(derived, book_corrected, check (bool: derived = corrected form), printed_check (bool: False),
    printed_minus_derived (what the printed form has too much: −(κ/2) ∇²q, written with the symbol q ≡ mean(T′²))).
    Negligible in practice (molecular transport of variance), but it is what the algebra gives.
    Validation: V2.
    Assumptions: constant thermal diffusivity; both continuity equations.
    """
    F = _sym_fields()
    sp, x, t, U, u = F["sp"], F["x"], F["t"], F["U"], F["u"]
    A = lambda e: _average(e)  # noqa: E731
    Tp, Tb, kap = F["Tp"], F["Tb"], F["kappa_th"]
    Tt = Tb + Tp
    heat = sp.diff(Tt, t) + sum(sp.diff((U[j] + u[j]) * Tt, x[j]) for j in range(3)) - kap * sum(sp.diff(Tt, x[j], 2) for j in range(3))
    fl = sp.expand(heat.doit() - _average(heat))
    derived = _average(Tp * fl)
    half = A(Tp * Tp) / 2

    def form(molecular_factor):
        return (sp.diff(half, t) + sum(U[j] * sp.diff(half, x[j]) for j in range(3))
                + sum(A(u[j] * Tp) * sp.diff(Tb, x[j]) for j in range(3))
                + sum(sp.diff(A(Tp * Tp * u[j]) / 2 - kap * molecular_factor * sp.diff(A(Tp * Tp), x[j]), x[j]) for j in range(3))
                + kap * sum(A(sp.diff(Tp, x[j]) ** 2) for j in range(3)))                     # Eq. (12.112)

    corrected, printed = form(sp.Rational(1, 2)), form(1)
    return dict(derived=derived, book_corrected=corrected, check=_zero(derived - corrected, F),
                printed_check=_zero(derived - printed, F),
                printed_minus_derived=-kap / 2 * sum(sp.Derivative(sp.Function("q")(*x, t), x[j], 2) for j in range(3)),
                fields=F)


def overlap_matching_sympy() -> dict:
    """The matching argument that gives the logarithmic law, step by step.

    Book: §12.9, Eqs. (12.85) ``dU/dy = (u_*²/ν) df/dy⁺``, (12.86) ``−dU/dy = (u_*/δ) dF/dξ``, (12.87)
    ``−ξ dF/dξ = y⁺ df/dy⁺`` (= 1/κ, a function of ξ equal to a function of y⁺), (12.88) ``f = (1/κ) ln y⁺ + B`` and
    (12.89) ``F = −(1/κ) ln ξ + A``.
    Returns dict(eq_12_85, eq_12_86 (sympy expressions of dU/dy from each law), matching (left and right of (12.87)),
    f, F_outer (the integrated laws), checks: dict(gradients_match, f_satisfies, F_satisfies, friction_law) — the last
    one is the sum of the two laws, U_∞/u_* = (1/κ) ln δ⁺ + A + B (y eliminated)).
    Validation: V2 all checks True.
    Assumptions: an overlap region exists where the inner law (12.80) and the defect law (12.84) both hold.
    """
    import sympy as sp

    y, us, nu, dl, kap, A_, B_, Uinf = sp.symbols("y u_* nu delta kappa A B U_inf", positive=True)
    f, Fo = sp.Function("f"), sp.Function("F")
    yp, xi = sp.symbols("y_plus xi", positive=True)
    dU_inner = sp.diff(us * f(y * us / nu), y)                         # Eq. (12.85)
    dU_outer = sp.diff(Uinf - us * Fo(y / dl), y)                      # Eq. (12.86): dU/dy = −(u_*/δ) F′
    f_sol = sp.log(yp) / kap + B_                                      # Eq. (12.88)
    F_sol = -sp.log(xi) / kap + A_                                     # Eq. (12.89)
    lhs = (y / us) * dU_inner.subs(f, sp.Lambda(yp, f_sol)).doit()
    rhs = (y / us) * dU_outer.subs(Fo, sp.Lambda(xi, F_sol)).doit()
    checks = dict(gradients_match=sp.simplify(lhs - rhs) == 0,
                  f_satisfies=sp.simplify(yp * sp.diff(f_sol, yp) - 1 / kap) == 0,
                  F_satisfies=sp.simplify(-xi * sp.diff(F_sol, xi) - 1 / kap) == 0)
    inner_U = us * f_sol.subs(yp, y * us / nu)
    outer_U = Uinf - us * F_sol.subs(xi, y / dl)
    fr = sp.solve(sp.Eq(inner_U, outer_U), Uinf)[0] / us
    checks["friction_law"] = sp.simplify(sp.expand_log(fr - (sp.log(dl * us / nu) / kap + A_ + B_), force=True)) == 0
    return dict(eq_12_85=dU_inner, eq_12_86=dU_outer, matching=(-xi * sp.Derivative(Fo(xi), xi), yp * sp.Derivative(f(yp), yp)),
                f=f_sol, F_outer=F_sol, friction_law=fr, checks=checks)


def plane_jet_similarity_sympy() -> dict:
    """Substitute the self-preserving forms into the thin-layer momentum equation and read off the three coefficients.

    Book: §12.8, the two unnumbered equations before (12.63), Eq. (12.63) itself and the trial families of (12.74).
    Method: with I(ξ) = ∫_0^ξ F (so F = I′), continuity gives V = −∂/∂x[U_CL δ I(y/δ)]; the residual
    ``U U_x + V U_y + ∂mean(uv)/∂y`` with −mean(uv) = Ψ G is multiplied by δ/U_CL² and compared with
    ``{δU′/U}F² − {δU′/U + δ′}F′I − {Ψ/U²}G′``.
    Returns dict(equation (the reduced left − right), coefficients (c1, c2, c3 as expressions), check (bool),
    power_family (c1, c2, c3 for δ = x^m, U_CL = x^n, Ψ = x^{2n+m−1}: (n, m + n, 1)·x^{m−1}),
    exponential_family (for δ = e^{ax}, U_CL = e^{−ax}, Ψ = e^{−ax}: c2 = 0 — slip #5), momentum_flux_exponent
    (2n + m), stress_integral_check (bool: with c1 = −½, c2 = ½ the equation integrates to C_3 G = −½ F I)).
    Validation: V2.
    Assumptions: thin-layer equation (12.61); self-preserving forms (12.56), (12.57).
    """
    import sympy as sp

    x, y, xi = sp.symbols("x y xi", positive=True)
    dl, Ucl, Psi = sp.Function("delta")(x), sp.Function("U_CL")(x), sp.Function("Psi")(x)
    I, G = sp.Function("I"), sp.Function("G")
    eta = y / dl
    U = Ucl * sp.diff(I(xi), xi).subs(xi, eta)                         # Eq. (12.56): F = I′
    V = -sp.diff(Ucl * dl * I(eta), x)                                 # Eq. (12.58) integrated, V(0) = 0
    muv = Psi * G(eta)                                                 # Eq. (12.57): −mean(uv)
    resid = U * sp.diff(U, x) + V * sp.diff(U, y) - sp.diff(muv, y)    # Eq. (12.61), left − right
    resid = (resid * dl / Ucl ** 2).doit().subs(y, xi * dl)
    resid = sp.simplify(resid)
    F_ = sp.diff(I(xi), xi)
    c1 = dl * sp.diff(Ucl, x) / Ucl
    c2 = c1 + sp.diff(dl, x)
    c3 = Psi / Ucl ** 2
    target = c1 * F_ ** 2 - c2 * sp.diff(F_, xi) * I(xi) - c3 * sp.diff(G(xi), xi)   # Eq. (12.63)
    check = sp.simplify(sp.expand(resid.doit() - target.doit())) == 0
    m, n, a = sp.symbols("m n a", real=True)
    fam = lambda d, uu_, ps: tuple(sp.simplify(c.subs({Psi: ps, Ucl: uu_, dl: d}).doit()) for c in (c1, c2, c3))  # noqa: E731
    power = fam(x ** m, x ** n, x ** (2 * n + m - 1))
    expo = fam(sp.exp(a * x), sp.exp(-a * x), sp.exp(-a * x))
    Fs = sp.Function("F")
    lhs = -sp.Rational(1, 2) * Fs(xi) ** 2 - sp.Rational(1, 2) * sp.diff(Fs(xi), xi) * sp.Integral(Fs(xi), xi)
    stress_ok = sp.simplify(lhs - sp.diff(-sp.Rational(1, 2) * Fs(xi) * sp.Integral(Fs(xi), xi), xi).doit()) == 0
    return dict(equation=resid, coefficients=dict(c1=c1, c2=c2, c3=c3), check=check, power_family=power,
                exponential_family=expo, momentum_flux_exponent=2 * n + m, stress_integral_check=stress_ok)


def gradient_moments_isotropic_sympy() -> dict:
    """Velocity-gradient moments of isotropic turbulence from the two-point tensor, and the factor 30 of (12.43).

    Book: §12.6, Eq. (12.43) and the route of Exercise 12.19:
    ``mean(∂u_i/∂x_k ∂u_j/∂x_l) = −(∂²R_ij/∂r_k∂r_l)_{r=0}`` with R_ij from (12.41) and f = 1 − r²/λ_f² + O(r⁴).
    Returns dict(m11 = mean((∂u_1/∂x_1)²), m12 = mean((∂u_1/∂x_2)²), cross = mean((∂u_1/∂x_2)(∂u_2/∂x_1)) — in units
    of u2/λ_f²: (2, 4, −1), **three different numbers** — and eps_factor = 6·(2 + 4 − 1) = 30, eps_factor_g = 15,
    ratio_lambda (λ_g²/λ_f² = ½), eps_from_m11 (= 15: ε̄ = 15 ν mean((∂u_1/∂x_1)²))).
    Validation: V2 exact rationals.
    Assumptions: homogeneous isotropic incompressible turbulence; f expanded to fourth order in r.
    """
    import sympy as sp

    r1, r2, r3, lam, u2, c4 = sp.symbols("r1 r2 r3 lambda_f u2 c4", real=True)
    rv = (r1, r2, r3)
    rr = r1 ** 2 + r2 ** 2 + r3 ** 2
    f = 1 - rr / lam ** 2 + c4 * rr ** 2                    # f(r) to O(r⁴); (r/2) f′ = −r²/λ² + 2 c4 r⁴
    half_r_fp = -rr / lam ** 2 + 2 * c4 * rr ** 2

    def R(i, j):                                            # Eq. (12.41), written without 1/r² (polynomial)
        d = 1 if i == j else 0
        return u2 * (f * d + half_r_fp * d - (-1 / lam ** 2 + 2 * c4 * rr) * rv[i] * rv[j])

    at0 = {r1: 0, r2: 0, r3: 0}
    mom = lambda i, k, j, l: sp.simplify(-sp.diff(R(i, j), rv[k], rv[l]).subs(at0) * lam ** 2 / u2)  # noqa: E731
    m11, m12, cross = mom(0, 0, 0, 0), mom(0, 1, 0, 1), mom(0, 1, 1, 0)
    g = f + half_r_fp                                       # g = f + (r/2) f′
    lam_g2 = sp.simplify(-2 / sp.diff(g.subs({r2: 0, r3: 0}), r1, 2).subs(at0))
    div = [sp.simplify(sum(sp.diff(R(i, j), rv[j]) for j in range(3))) for i in range(3)]
    return dict(m11=m11, m12=m12, cross=cross, eps_factor=6 * (m11 + m12 + cross), eps_factor_g=sp.simplify(
                6 * (m11 + m12 + cross) * lam_g2 / lam ** 2), ratio_lambda=sp.simplify(lam_g2 / lam ** 2),
                eps_from_m11=sp.simplify(6 * (m11 + m12 + cross) / m11), divergence_free=all(d == 0 for d in div))


def isotropic_tensor_divergence_sympy(f=None) -> dict:
    """Check that the incompressible isotropic tensor is divergence-free: ``∂R_ij/∂r_j = 0`` for a given f(r).

    Book: §12.6, Eq. (12.41) (obtained from (12.40) by imposing ∂u_j/∂x_j = 0, Exercise 12.18).
    Parameters: f a sympy expression in the symbol ``r`` (default exp(−r²/L²)); any smooth even f works.
    Returns dict(divergence (list of 3 simplified expressions), check (bool), g (= f + (r/2) f′), general_fails (bool:
    the general form (12.40) with an *independent* g, here g = f, is not divergence-free)).
    Validation: V2 three trial f (Gaussian, exp(−r), (1 + r²)^{−1}).
    Assumptions: homogeneous isotropic turbulence; f a smooth function of r.
    """
    import sympy as sp

    r = sp.Symbol("r", positive=True)
    Lc = sp.Symbol("L", positive=True)
    if f is None:
        f = sp.exp(-r ** 2 / Lc ** 2)
    r1, r2, r3 = sp.symbols("r1 r2 r3", positive=True)
    rv = (r1, r2, r3)
    rad = sp.sqrt(r1 ** 2 + r2 ** 2 + r3 ** 2)
    fr = f.subs(r, rad)
    fp = sp.diff(f, r).subs(r, rad)

    def R(i, j, incompressible=True):
        d = 1 if i == j else 0
        if incompressible:
            return fr * d + rad / 2 * fp * (d - rv[i] * rv[j] / rad ** 2)   # Eq. (12.41) / u2
        return fr * d                                                          # Eq. (12.40) with g = f (F = 0)

    div = [sp.simplify(sum(sp.diff(R(i, j), rv[j]) for j in range(3))) for i in range(3)]
    div_gen = [sp.simplify(sum(sp.diff(R(i, j, False), rv[j]) for j in range(3))) for i in range(3)]
    return dict(divergence=div, check=all(d == 0 for d in div), g=sp.simplify(f + r / 2 * sp.diff(f, r)),
                general_fails=any(d != 0 for d in div_gen))


SYMPY_ENGINES = ("rans", "reynolds_stress_budget", "tke_budget", "mean_energy_budget", "temperature_variance",
                 "overlap_matching", "plane_jet_similarity", "gradient_moments_isotropic", "isotropic_tensor_divergence")


def sympy_summary(name: str, cache: bool = True) -> dict:
    """Run one symbolic engine and return a small, JSON-serialisable summary (cached on disk for the notebook).

    Book: the equations each engine derives (see the engine's docstring).  Parameters: name — one of
    :data:`SYMPY_ENGINES`; cache — read/write ``outputs/ch12/cache/sympy_<name>.json`` (invalidate by deleting the
    folder after changing an engine).
    Returns dict(name, checks {label: bool}, latex {label: LaTeX string of the main expressions}, seconds).
    Assumptions: as the engine named.
    Validation: V2 (inherits the engine's checks; the cache holds only booleans and LaTeX strings).
    """
    import json
    import time

    import sympy as sp

    if name not in SYMPY_ENGINES:
        raise ValueError(f"name must be one of {SYMPY_ENGINES}")
    path = _root() / "outputs" / "ch12" / "cache" / f"sympy_{name}.json"
    if cache and path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    t0 = time.time()
    res = globals()[f"{name}_sympy"]()
    checks, latex = {}, {}

    def scan(prefix, obj):
        if isinstance(obj, (bool, np.bool_)) or obj in (sp.true, sp.false):
            checks[prefix] = bool(obj)
        elif isinstance(obj, dict):
            for k, v in obj.items():
                if k != "fields":
                    scan(f"{prefix}.{k}" if prefix else str(k), v)
        elif isinstance(obj, (list, tuple)):
            for k, v in enumerate(obj):
                scan(f"{prefix}[{k}]", v[1] if (isinstance(v, tuple) and len(v) == 2 and isinstance(v[0], str)) else v)
        elif isinstance(obj, sp.Basic):
            s = sp.latex(obj)
            latex[prefix] = s if len(s) < 4000 else s[:4000] + r"\;\dots"

    scan("", res)
    out = dict(name=name, checks=checks, latex=latex, seconds=round(time.time() - t0, 2))
    if cache:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(out, indent=1), encoding="utf-8")
    return out


# ======================================================================================================================
# chapter tables: conventions, printed slips, explainer tables
# ======================================================================================================================
def conventions():
    """The symbol and normalisation conventions of chapter 12 as a table (for the opening block of the notebook).

    Book: §12.3–§12.12 (the chapter re-uses many letters; this fixes what each one means in fluidpy).
    Returns a ``pandas.DataFrame`` with columns symbol, meaning, code_name, unit, note.
    Assumptions: none (a table).
    Validation: V1 every code_name that is a function exists in this module.
    """
    import pandas as pd

    rows = [
        ("over-bar", "ensemble average (12.1); time (12.2) or volume (12.3) average only where a name says so", "ensemble_average", "—",
         "measurements use time averages: ergodicity is assumed"),
        ("ũ = U + u", "total = mean + fluctuation (12.24); lower case u, v, w, primed T′ are fluctuations", "reynolds_decompose", "m/s", "u is NOT the velocity in this chapter"),
        ("mean(u²)", "variance of ONE velocity component in §12.6–§12.7 (ē = 3/2 of it when isotropic)", "u2", "m²/s²", ""),
        ("ē", "turbulent kinetic energy per unit mass, ½ mean(u_i u_i) — the 'k' of k–ε", "e", "m²/s²", "ch01's e was internal energy"),
        ("ε̄", "dissipation rate of ē (12.42)", "eps", "m²/s³ = W/kg", ""),
        ("κ (wall)", "von Kármán constant (12.88)", "kappa", "—", "required keyword, no default"),
        ("κ (heat)", "thermal diffusivity (12.31)", "kappa_th", "m²/s", "κ_m, κ_T, κ_mT: molecular scalar, eddy heat, eddy scalar"),
        ("k_1, K", "streamwise and three-dimensional wavenumber", "k1, K", "rad/m", "k is also conductivity in (12.32)"),
        ("S_e(ω), S_11(k_1)", "two-sided spectra, 1/2π in the forward transform, ∫_{−∞}^{∞} = variance (12.20)–(12.22), (12.55)", "two_sided=True", "unit²·s, unit²·m",
         "Fig. 12.12 plots 2S_11 (one-sided)"),
        ("E(K) (book: S(K))", "three-dimensional spectrum on K ≥ 0, ∫E dK = ē", "inertial_spectrum_3d", "m³/s²", ""),
        ("Λ_t, Λ_f, Λ_g", "integral scales (12.18), (12.39)", "integral_scale", "s or m", "ch11's Λ was a dissipation"),
        ("λ_t, λ_f, λ_g, λ_T", "Taylor microscales (12.19), (12.39); λ_g = λ_f/√2", "taylor_microscale", "s or m", "earlier chapters: wavelength"),
        ("η, u_K, η_T", "Kolmogorov length and velocity (12.50); Batchelor scale", "kolmogorov_scales", "m, m/s", "earlier chapters: η = similarity variable"),
        ("f, g", "longitudinal / transverse correlation (12.38); f also the law-of-the-wall function (12.80)", "f, g", "—", ""),
        ("F, G", "profile functions of jet velocity and stress (12.56)–(12.57); F also the defect law (12.84)", "gaussian_profile", "—", ""),
        ("δ", "jet width (δ = x by convention), boundary-layer thickness, channel half-height", "delta", "m", ""),
        ("h, d, a", "FULL channel height (12.90); pipe diameter (12.91) or slot width; pipe radius", "h, d, a", "m", "ch08 used half-heights in places"),
        ("u_*, l_ν, y⁺, U⁺", "friction velocity (12.81), viscous length, wall units (12.80)", "wall_units", "m/s, m, —", ""),
        ("Π", "Coles' wake strength", "Pi", "—", "ch01's Π were dimensionless groups"),
        ("Γ, Γ_a", "Kundu: Γ ≡ dT/dz, Γ_a = −g/C_p ≈ −9.8 K/km; meteorology: Γ ≡ −dT/dz, Γ_d ≈ +9.8 K/km", "dTdz, Gamma_a", "K/m",
         "computed in Kundu's sign, both conventions reported"),
        ("T̄, T′ (§12.11)", "POTENTIAL temperature in every buoyancy term", "wT", "K", "N² = gα(dT/dz − Γ_a) with the in-situ gradient"),
        ("mean(wT′), H", "kinematic heat flux and heat flux, positive upward; H = ρ C_p mean(wT′)", "wT, H", "K m/s, W/m²", ""),
        ("Rf, Ri", "flux (12.107) and gradient (12.108) Richardson numbers; Ri = Pr_T·Rf", "flux_richardson, gradient_richardson_thermal", "—",
         "Rf_cr ≈ ¼ is an observation; Ri > ¼ of ch11 is a theorem"),
        ("L_M", "Monin–Obukhov length (12.110): positive stable, negative unstable", "L_M", "m", ""),
        ("α", "thermal expansion coefficient; also the no-sum Greek index of §12.12", "alpha", "1/K", ""),
        ("σ² (12.126)", "variance of a spreading Gaussian, σ² = 2νt", "variance", "m²", "ch03's vortex core used 4νt"),
        ("D_T", "eddy diffusivity of dispersion (12.127), a function of time", "eddy_diffusivity_taylor", "m²/s", ""),
    ]
    return pd.DataFrame(rows, columns=["symbol", "meaning", "code_name", "unit", "note"])


def book_slips() -> list:
    """The printed slips of chapter 12 that fluidpy corrects, each with the way a test tells the two versions apart.

    Book: §12.5–§12.12 (found by reading the page images against dimensions, limits and the book's own neighbouring
    equations).  Returns a list of 15 dicts(id, where, printed, corrected, how_to_tell, coded_in, test) — ``test`` names
    the planted wrong variant a test must fail ("—" for pure wording slips).  Descriptions are ours;
    ``pd.DataFrame(book_slips())`` gives the table of the notebook's opening block.
    Assumptions: none (a table).   Validation: V1 every ``printed=True`` switch named here exists and fails its check.
    """
    rows = [
        (1, "(12.54)", "exponent of k_1 is +5/3", "−5/3", "units: ε^{2/3} k^{−5/3} has m³/s²; slope of the spectrum",
         "inertial_spectrum_1d(printed=True)"),
        (2, "(12.47)", "label under the left side names the mean-flow energy", "the turbulent ē", "it is half the trace of (12.35)", "tke_budget_sympy"),
        (3, "(12.70)", "source integral evaluated at y = 0", "at x = 0, as in (12.62)", "the integral runs over y", "plane_jet_mass_fraction"),
        (4, "after (12.67)", "undetermined constants called C_3, C_4; later C_5, C_6 said to be tabulated", "C_3 and C_5; the table's pair multiplies (12.72) and (12.73)",
         "C_5 = C_4 (ρ/J_s)^{1/2}", "free_shear_centerline"),
        (5, "(12.74)", "exponential family δ ~ e^{ax}, U_CL ~ e^{−ax}, Ψ ~ e^{−ax} offered as a solution",
         "its middle coefficient vanishes and it breaks (12.62); the power family with m + 2n = 0 is the valid one",
         "evaluate the three coefficients", "general_similarity_check, plane_jet_similarity_sympy"),
        (6, "(12.97)", "pressure gradient ∂P/∂x_j", "∂P/∂x_i (free index)", "index balance", "rans_eddy_viscosity_residual"),
        (7, "after (12.89)", "control-volume proof of (12.90)–(12.91) cited as Exercise 12.31", "Exercise 12.32", "—", "—"),
        (8, "before (12.103)", "modelled viscous transport written with u_j S′_ij", "u_i S′_ij as in (12.47)", "index balance", "one_equation_closure"),
        (9, "after (12.105)", "'five' constants but six printed entries (one repeated)", "five: C_μ, C_ε1, C_ε2, σ_e, σ_ε", "count", "K_EPSILON_CONSTANTS"),
        (10, "(12.75), (12.106)", "triple correlation as ½ mean(ē v) and mean(e w)", "½ mean(u_i² u_j) throughout", "compare with (12.47)", "jet_tke_budget, stratified_tke_budget"),
        (11, "after (12.121)", "refers to (11.119)", "(12.119)", "—", "—"),
        (12, "(12.129)", "condition t ≪ Λ_t", "t ≫ Λ_t", "at t = 0.01 Λ_t the constant is 100× the exact D_T", "eddy_diffusivity_asymptote(printed=True)"),
        (13, "Fig. 12.27 caption", "width ∝ x^{1/2} near the source, ∝ x far away", "linear near, square root far", "(12.121), (12.123); local slope of Z_rms(x)", "smoke_plume_width"),
        (14, "Exercise 12.18a", "cites (12.39) for R_ij", "(12.40)", "—", "isotropic_correlation_tensor"),
        (15, "(12.112)", "molecular transport κ ∂mean(T′²)/∂z", "κ ∂(½ mean(T′²))/∂z", "derive it: T′ × (T′-equation), averaged", "temperature_variance_sympy"),
    ]
    # "test": the call (or check) with which a test shows that the printed form fails; "—" = a wording slip, no test.
    tests = {1: "inertial_spectrum_1d(printed=True): slope +5/3 and wrong units",
             3: "scalar_flux_per_span: the flux is an integral over y at fixed x (it has no meaning 'at y = 0')",
             5: "general_similarity_check on the exponential family: c2 = 0 and U_CL^2*delta varies with x",
             6: "rans_eddy_viscosity_residual(printed=True): non-zero residual on a manufactured solution",
             12: "eddy_diffusivity_asymptote(which='long', printed=True) against eddy_diffusivity_exponential "
                 "(dispersion_regime names the regime)",
             13: "smoke_plume_width: local slope d ln Z/d ln x is 1 near the source and 1/2 far away",
             15: "temperature_variance_sympy()['printed_check'] is False"}
    keys = ("id", "where", "printed", "corrected", "how_to_tell", "coded_in")
    return [dict(zip(keys, r), test=tests.get(r[0], "—")) for r in rows]


def explainer_tables(Re_taus=(180.0, 550.0, 1000.0, 5200.0), kappa: float = 0.41, A_plus: float = 26.0,
                     n_out: int = 96, n: int = 1600, vd_decades: tuple = (-1.0, 4.5), vd_points: int = 112) -> dict:
    """Tables computed by this module for the explainers: the model channel (E5, E8) and the van Driest wall profile (E8).

    Book: §12.7 (12.46)–(12.47), §12.9 (12.76)–(12.77), §12.10 (12.99)–(12.101) — see :func:`channel_energy_budget` and
    :func:`mixing_length_wall_profile`.  κ = 0.41 is the "classical" preset of ``WT.LOG_LAW_CONSTANTS``; A⁺ = 26 is van
    Driest's damping constant.  **Ours, computed by fluidpy — no book data, no external benchmark.**
    Parameters: Re_taus friction Reynolds numbers of the channel table; kappa; A_plus; n_out points per channel profile
    (y⁺ = 0 plus a geometric grid from 0.2 to Re_τ); n grid points of the underlying solution; vd_decades (log10 of the
    first and last y⁺) and vd_points of the wall-profile table.
    Returns dict(note, kappa, A_plus,
    channel {"<Re_tau>": dict(yplus, Uplus, uv_plus, dUdy_plus, production, pressure_work, mean_dissipation, turb_sink,
    transport, U_bulk_plus, U_cl_plus, Cf, yplus_peak_production, production_peak, integrals)},
    van_driest dict(yplus, Uplus (damped, A⁺), Uplus_undamped, B (intercept with damping), B_undamped)) — lists of
    floats rounded to 6 significant digits.  The channel columns are evaluated point by point with
    :func:`channel_energy_budget_at` (U⁺ by adaptive quadrature), so they do not carry the grid error of the arrays.
    Assumptions: as the two functions named.   Validation: V1 each row reproduces ``channel_energy_budget_at`` /
    ``mixing_length_wall_profile`` to the 6 digits stored; the JSON on disk equals a fresh call.
    """
    sig = lambda a: [float(f"{v:.6g}") for v in np.atleast_1d(a)]  # noqa: E731
    out = dict(note="ours, computed by fluidpy (ch12_turbulence.explainer_tables): mixing-length channel and van Driest "
                    "wall profile; a model, not DNS",
               kappa=kappa, A_plus=A_plus, channel={})
    cols = ("Uplus", "uv_plus", "slope", "production", "pressure_work", "mean_dissipation", "turb_sink", "transport")
    for Re in Re_taus:
        b = channel_energy_budget(Re, kappa, A_plus, n=n)
        yq = np.concatenate([[0.0], np.geomspace(0.2, Re, n_out - 1)])
        yq[-1] = Re
        pts = [channel_energy_budget_at(float(q), Re, kappa, A_plus) for q in yq]
        entry = {("dUdy_plus" if c == "slope" else c): sig([p[c] for p in pts]) for c in cols}
        entry.update(yplus=sig(yq), U_bulk_plus=float(f"{b['U_bulk_plus']:.6g}"), U_cl_plus=float(f"{pts[-1]['Uplus']:.6g}"),
                     Cf=float(f"{b['Cf']:.6g}"), yplus_peak_production=float(f"{b['yplus_peak_production']:.4g}"),
                     production_peak=float(f"{float(np.max(b['production'])):.6g}"),
                     integrals={k: float(f"{v:.6g}") for k, v in b["integrals"].items()})
        out["channel"][str(int(Re))] = entry
    yv = np.logspace(vd_decades[0], vd_decades[1], int(vd_points))
    out["van_driest"] = dict(yplus=sig(yv),
                             Uplus=sig(mixing_length_wall_profile(yv, kappa, damping="van_driest", A_plus=A_plus)["Uplus"]),
                             Uplus_undamped=sig(mixing_length_wall_profile(yv, kappa)["Uplus"]),
                             B=float(f"{mixing_length_intercept(kappa, A_plus):.6g}"),
                             B_undamped=float(f"{mixing_length_intercept(kappa):.6g}"))
    return out


def write_reference_tables(dest=None) -> dict:
    """Write our computed tables to ``reference/ch12/explainer_tables.json`` (or ``dest``) and return the paths.

    Book: §12.7 (12.46)–(12.47), §12.9–§12.10 (the model channel and wall profile of :func:`explainer_tables`).
    These are fluidpy's own numbers (no book data, no external benchmark): the explainers load them and prove parity
    with :func:`channel_energy_budget_at` and :func:`mixing_length_wall_profile` at two or more points.
    Parameters: dest folder (default ``reference/ch12``).  Returns {"explainer_tables": Path}.  Called by
    ``scripts/ch12_tables.py``; re-running it regenerates the file byte for byte (deterministic, no random numbers).
    Assumptions: none beyond those of the two functions.   Validation: V1 file content equals ``explainer_tables()``.
    """
    import json

    folder = Path(dest) if dest is not None else _root() / "reference" / "ch12"
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / "explainer_tables.json"
    path.write_text(json.dumps(explainer_tables(), indent=1), encoding="utf-8")
    return {"explainer_tables": path}
