"""Statistics of fluctuating signals and fields: moments, averages, correlations, spectra, synthetic turbulence-like data.

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e, §12.3 (moments and the three averages, Eqs. 12.1–12.11), §12.4 (correlations
and spectra, Eqs. 12.12–12.23), §12.5 (Reynolds decomposition and stress, Eqs. 12.24–12.26, 12.30), §12.6 (longitudinal and
transverse correlations, Eqs. 12.38–12.39, 12.45).  First used in chapter 12; written as a core module because every later
chapter with a fluctuating quantity needs the same operations (eddy fluxes and spectra in ch13, unsteady loads in ch14).

Conventions (fixed for the whole module)
----------------------------------------
* **Ensembles carry the member index on axis 0**: ``samples[n, ...]`` is realization ``n``.  Time (or space) is the last axis.
* **Averages.**  The over-bar of the book is an ensemble average (12.1) unless the function name says ``time_`` (12.2) or
  ``volume_`` (12.3).  Which one a function uses is stated in its docstring.
* **Spectra are two-sided in angular frequency ω [rad/s] or wavenumber k [rad/m]** with the 1/2π in the forward transform,
  the book's normalisation: ``S(ω) = (1/2π)∫R(τ)e^{−iωτ}dτ`` (12.20), ``R(τ) = ∫S(ω)e^{+iωτ}dω`` (12.21), so that
  ``∫_{−∞}^{∞} S dω = variance`` (12.22).  Because S is even, functions return it on ω ≥ 0 only; the variance is then
  ``2∫_0^∞ S dω`` (:func:`spectrum_variance` does this).  ``one_sided=True`` returns ``2S`` (integrates to the variance over
  ω ≥ 0).  numpy's FFT works in cyclic frequency without the 1/2π; :func:`periodogram` converts.
* **Random numbers** always come from ``np.random.default_rng(seed)``; the seed is an argument and is named in the docstring.
* **Synthetic fields are kinematic**: random-phase fields with a prescribed spectrum are divergence-free and have the requested
  two-point statistics, but they have no energy cascade (zero velocity-gradient skewness).  Never read dynamics off them.
* **Two vs three dimensions.**  ``g = f + (r/2) f'`` (12.41) and the factor 15 of (12.43) are three-dimensional results; a 2-D
  solenoidal isotropic field obeys ``g = d(r f)/dr``.  Use ``dim=3`` whenever (12.41)–(12.43) or the cascade is the subject.
"""
from __future__ import annotations

from typing import Callable

import numpy as np
from scipy.signal import lfilter

from ._util import as_scalar_if_0d

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d

__all__ = [
    "make_ensemble", "ensemble_average", "moment", "central_moment", "statistics", "standard_error",
    "standard_error_of_mean", "time_average", "volume_average", "check_averaging_rules", "product_average_split",
    "reynolds_decompose", "correlation", "correlation_coefficient", "correlated_pair", "autocorrelation",
    "cross_correlation", "integral_scale", "correlation_time", "effective_samples", "taylor_microscale",
    "integral_scale_from_spectrum", "spectrum_from_correlation", "correlation_from_spectrum", "spectrum_variance",
    "periodogram", "taylor_frozen", "frequency_to_wavenumber_spectrum", "correlation_spectrum_pair", "smooth_signal",
    "synthetic_solenoidal_field", "white_noise_field", "shell_spectrum", "spatial_correlation",
    "longitudinal_transverse_correlation", "velocity_covariance", "reynolds_stress", "anisotropy_tensor",
    "turbulent_kinetic_energy",
]


# ======================================================================================================================
# §12.3  ensembles, moments, the three averages
# ======================================================================================================================
def make_ensemble(n_members: int, t, mean_fn, sigma: float, tau_c: float, seed: int = 0) -> np.ndarray:
    """Ensemble of ``n_members`` records u(t:n) = mean(t) + Ornstein–Uhlenbeck fluctuation (variance σ², memory τ_c).

    Book: §12.3 (realization, ensemble; Figs. 12.2–12.3).  The fluctuation is a stationary Gaussian process with
    autocorrelation ``R(τ) = σ² exp(−|τ|/τ_c)``, generated with the *exact* update
    ``u_{k+1} = ρ_k u_k + σ sqrt(1 − ρ_k²) ξ_k``, ``ρ_k = exp(−(t_{k+1} − t_k)/τ_c)`` (no time-step error).

    Parameters
    ----------
    n_members : number of realizations N.
    t : (nt,) times [s], increasing (uniform or not).
    mean_fn : callable t → mean [unit of u], or an array of shape (nt,), or a scalar (stationary mean).
    sigma : standard deviation of the fluctuation [unit of u], ≥ 0.
    tau_c : memory (integral) time of the fluctuation [s], > 0.
    seed : seed of ``np.random.default_rng`` (default 0).

    Returns
    -------
    samples : (n_members, nt) array; member index on axis 0.

    Assumptions: Gaussian, stationary fluctuation; members independent.  The OU process has a cusp in r(τ) at τ = 0, so
    its Taylor microscale (12.19) does not exist — use :func:`smooth_signal` for that.
    Validation: V1 ensemble mean → mean_fn, variance → σ², lag correlation → exp(−τ/τ_c) (5 standard errors); V7 σ = 0.
    """
    t = _F(t)
    if t.ndim != 1 or t.size < 1:
        raise ValueError("t must be a 1-D array")
    if tau_c <= 0 or sigma < 0:
        raise ValueError("need tau_c > 0 and sigma >= 0")
    rng = np.random.default_rng(seed)
    mean = _F(mean_fn(t)) if callable(mean_fn) else _F(mean_fn)
    mean = np.broadcast_to(mean, t.shape)
    xi = rng.standard_normal((int(n_members), t.size))
    u = np.empty_like(xi)
    u[:, 0] = sigma * xi[:, 0]  # stationary start
    dt = np.diff(t)
    if t.size > 1 and np.allclose(dt, dt[0], rtol=1e-10, atol=0.0):
        rho = np.exp(-dt[0] / tau_c)
        amp = sigma * np.sqrt(-np.expm1(-2.0 * dt[0] / tau_c))
        zi = (rho * u[:, 0])[:, None]
        u[:, 1:], _ = lfilter([1.0], [1.0, -rho], amp * xi[:, 1:], axis=1, zi=zi)  # exact OU update
    else:
        for k in range(t.size - 1):
            rho = np.exp(-dt[k] / tau_c)
            u[:, k + 1] = rho * u[:, k] + sigma * np.sqrt(-np.expm1(-2.0 * dt[k] / tau_c)) * xi[:, k + 1]
    return u + mean


def moment(samples, m: int, axis: int = 0):
    """m-th moment as an ensemble average, ``mean(u^m)`` over the member axis.

    Book: §12.3, Eq. (12.1) ``<u^m> = lim (1/N) Σ_n u(x, t:n)^m`` (finite N: the over-bar estimate).
    Parameters: samples [unit of u] with members on ``axis``; m integer ≥ 0.  Returns the moment [unit^m].
    Validation: V1 Gaussian moments (σ², 0, 3σ⁴) within 5 s.e.; m = 1 equals :func:`ensemble_average`.
    Assumptions: the over-bar is estimated by the average over the samples supplied (an estimate with sampling error ~ N^{-1/2}).
    """
    return _S(np.mean(_F(samples) ** m, axis=axis))  # Eq. (12.1)


def ensemble_average(samples, axis: int = 0):
    """Ensemble average (first moment) over the member axis.

    Book: §12.3, Eq. (12.10) ``ū(x, t) = (1/N) Σ_n u(x, t:n)`` (Eq. (12.1) with m = 1).
    Parameters: samples [unit of u], members on ``axis``.  Returns the mean [unit of u].
    Validation: V1 against a hand loop; mean of an OU ensemble → the prescribed mean (5 s.e.).
    Assumptions: the over-bar is estimated by the average over the samples supplied (an estimate with sampling error ~ N^{-1/2}); members are independent realizations.
    """
    return _S(np.mean(_F(samples), axis=axis))  # Eq. (12.10)


def central_moment(samples, m: int, axis: int = 0):
    """m-th central moment ``mean((u − ū)^m)`` over the member axis (ū = sample mean replaces <u>).

    Book: §12.3, Eq. (12.11).  m = 2, 3, 4 are the book's variance, skewness and kurtosis (un-normalised).
    Parameters: samples [unit], m integer.  Returns [unit^m].
    Validation: V1 identities ``mean((u−ū)²) = mean(u²) − ū²`` and the third-moment expansion to round-off.
    Assumptions: the over-bar is estimated by the average over the samples supplied (an estimate with sampling error ~ N^{-1/2}); the mean removed is the sample mean (biased by 1/N, no Bessel correction).
    """
    s = _F(samples)
    return _S(np.mean((s - np.mean(s, axis=axis, keepdims=True)) ** m, axis=axis))  # Eq. (12.11)


def statistics(samples, normalized: bool = False, axis: int = 0) -> dict:
    """Mean, variance, skewness, kurtosis, standard deviation and root mean square of an ensemble.

    Book: §12.3, Eqs. (12.10)–(12.11) and the definitions after them.  **The book's skewness and kurtosis are the raw
    third and fourth central moments** [unit³, unit⁴]; ``normalized=True`` divides them by σ³ and σ⁴ (the usual
    convention: 0 and 3 for a Gaussian).
    Returns dict(mean, variance, skewness, kurtosis, std, rms) — rms = sqrt(mean(u²)), equal to std when the mean is 0.
    Validation: V1 Gaussian samples (normalised skewness 0, kurtosis 3 within 5 s.e.); cross-check with scipy.stats.
    Assumptions: the over-bar is estimated by the average over the samples supplied (an estimate with sampling error ~ N^{-1/2}); skewness and kurtosis are the book's un-normalised central moments unless normalized=True.
    """
    s = _F(samples)
    mean = np.mean(s, axis=axis)
    var = central_moment(s, 2, axis)
    sk, ku = central_moment(s, 3, axis), central_moment(s, 4, axis)
    std = np.sqrt(var)
    if normalized:
        with np.errstate(divide="ignore", invalid="ignore"):
            sk, ku = sk / std ** 3, ku / std ** 4
    return dict(mean=_S(mean), variance=_S(var), skewness=_S(sk), kurtosis=_S(ku), std=_S(std),
                rms=_S(np.sqrt(np.mean(s ** 2, axis=axis))))


def standard_error(samples, axis: int = 0):
    """Standard error of the N-member ensemble mean, ``s/sqrt(N)`` with the unbiased sample standard deviation s.

    Book: §12.3, Fig. 12.3 (the scatter of the N-member mean shrinks with N; the N^{−1/2} law is ours, the standard result
    for independent samples).  Parameters: samples [unit], members on ``axis``.  Returns [unit].
    Validation: V1 slope −0.5 of the measured scatter against N on log axes.
    Assumptions: independent members (for a time record use effective_samples instead of N).
    """
    s = _F(samples)
    n = s.shape[axis]
    return _S(np.std(s, axis=axis, ddof=1) / np.sqrt(n))


def standard_error_of_mean(sigma, N):
    """Closed form σ/sqrt(N): scatter of a mean of N *independent* samples of standard deviation σ (scalar-callable).

    Book: §12.3 (Fig. 12.3) and §12.4 (N ≈ Δt/t_c independent samples in a record).  Returns [unit of σ].
    Assumptions: N independent samples of standard deviation sigma.
    Validation: V1 closed form; halves when N quadruples; matches the scatter of make_ensemble means within 5 %.
    """
    return _S(_F(sigma) / np.sqrt(_F(N)))


def _cumulative_integral(t: np.ndarray, y: np.ndarray, tq: np.ndarray) -> np.ndarray:
    """∫_{t0}^{tq} y dt for the piecewise-linear interpolant of (t, y) — exact (quadratic within each cell)."""
    h = np.diff(t)
    cum = np.concatenate([np.zeros(y.shape[:-1] + (1,)), np.cumsum(0.5 * (y[..., 1:] + y[..., :-1]) * h, axis=-1)], axis=-1)
    j = np.clip(np.searchsorted(t, tq, side="right") - 1, 0, t.size - 2)
    s = tq - t[j]
    slope = (y[..., j + 1] - y[..., j]) / h[j]
    return cum[..., j] + y[..., j] * s + 0.5 * slope * s ** 2


def time_average(t, u, window: float, m: int = 1):
    """Centred sliding time average of u^m over a window Δt: ``(1/Δt) ∫_{t−Δt/2}^{t+Δt/2} u^m dt``.

    Book: §12.3, Eq. (12.2).  **Time average** (not ensemble).
    Parameters: t (nt,) [s] increasing; u (..., nt) [unit], time on the last axis; window Δt [s] > 0; m power.
    Returns (..., nt) array [unit^m]; **NaN where the window leaves the record** (within Δt/2 of either end).
    Method: exact integral of the piecewise-linear interpolant of u^m (trapezoid rule when the window ends fall on
    samples); second order in the sample spacing.
    Validation: V1 against :func:`fluidpy.ch12_turbulence.time_average_exp_cos` (Example 12.1 closed form); a window of
    exactly one period removes a cosine; window → 0 returns the signal.
    Assumptions: statistically stationary record, long against the correlation time; uniform sampling; trapezoid rule inside the window; NaN where the window leaves the record.
    """
    t, y = _F(t), _F(u) ** m
    if window <= 0:
        raise ValueError("window must be > 0")
    lo, hi = t - 0.5 * window, t + 0.5 * window
    ok = (lo >= t[0] - 1e-12 * window) & (hi <= t[-1] + 1e-12 * window)
    out = np.full(y.shape, np.nan)
    if np.any(ok):
        a = _cumulative_integral(t, y, np.clip(lo[ok], t[0], t[-1]))
        b = _cumulative_integral(t, y, np.clip(hi[ok], t[0], t[-1]))
        out[..., ok] = (b - a) / window  # Eq. (12.2)
    return out


def volume_average(field, weights=None, m: int = 1):
    """Volume (spatial) average of field^m, ``(1/V) ∫ u^m dV``.

    Book: §12.3, Eq. (12.3).  **Spatial average** for a homogeneous field.
    Parameters: field (any shape) [unit]; weights = cell volumes [m³] of the same shape (None: uniform cells); m power.
    Returns a float [unit^m].  Validation: V1 uniform weights = np.mean; weighted mean of a linear field.
    Assumptions: uniform grid unless weights are given.
    """
    f = _F(field) ** m
    return float(np.mean(f)) if weights is None else float(np.sum(f * _F(weights)) / np.sum(_F(weights)))  # Eq. (12.3)


def reynolds_decompose(samples, axis: int = 0):
    """Reynolds decomposition ``ũ = U + u``: ensemble mean and fluctuations of every member.

    Book: §12.5, Eqs. (12.24)–(12.26): the mean is the expected value, the fluctuation has zero mean.
    Parameters: samples [unit], members on ``axis``.  Returns (mean, fluctuations); ``mean(fluctuations)`` is zero to
    round-off (12.26) and ``mean + fluctuations`` reproduces the samples.
    Validation: V1 both identities to 1e-12.
    Assumptions: the over-bar is estimated by the average over the samples supplied (an estimate with sampling error ~ N^{-1/2}); the fluctuations then average to zero exactly.
    """
    s = _F(samples)
    mean = np.mean(s, axis=axis)                      # Eq. (12.25)
    return mean, s - np.expand_dims(mean, axis)       # Eq. (12.24); its mean vanishes, Eq. (12.26)


def product_average_split(samples_u, samples_v, axis: int = 0):
    """The average of a product is not the product of the averages: ``mean(ũṽ) = Ū V̄ + mean(uv)``.

    Book: §12.3 (the relation after (12.9): ``mean(uv) ≠ ū v̄``) — the split itself is ours and is what makes the
    Reynolds stress appear in (12.30).
    Parameters: samples_u, samples_v [units a, b], members on ``axis``.
    Returns (mean_product, product_of_means, covariance) [unit a·b]; first = second + third to round-off.
    Validation: V1 identity to 1e-12; covariance equals ``np.cov(..., bias=True)``.
    Assumptions: the over-bar is estimated by the average over the samples supplied (an estimate with sampling error ~ N^{-1/2}); the identity holds to round-off for any sample set.
    """
    a, b = _F(samples_u), _F(samples_v)
    ma, mb = np.mean(a, axis=axis, keepdims=True), np.mean(b, axis=axis, keepdims=True)
    cov = np.mean((a - ma) * (b - mb), axis=axis)
    return _S(np.mean(a * b, axis=axis)), _S(np.squeeze(ma * mb, axis=axis)), _S(cov)


def check_averaging_rules(samples_u, samples_v, t, A: float = 2.0, m: int = 1, x=None) -> dict:
    """Residuals of the commutation rules of ensemble averaging, evaluated on two ensembles.

    Book: §12.3, Eqs. (12.4) sum, (12.5) constant factor, (12.6) time derivative, (12.7) time integral, (12.8) space
    derivative, (12.9) space integral, and ``mean(mean(u)) = mean(u)``.
    Parameters: samples_u, samples_v shaped (N, nt) or (N, nx, nt) [unit]; t (nt,) [s]; A constant; m power; x (nx,) [m]
    for 3-D input (for 2-D input the two spatial rules are evaluated along the last axis as a stand-in coordinate).
    Returns dict rule → max |left side − right side| (each rule uses the **same discrete operator** — second-order
    differences, trapezoid — on the members and on the mean, so every residual is at round-off: the rules hold because
    averaging is a linear operation), plus ``"product"`` = max |mean(uv) − ū v̄|, which is *not* small.
    Validation: V1 all six residuals < 1e-12 for random ensembles; "product" equals the covariance.
    Assumptions: the same discrete operators are applied to the members and to their mean (so linear rules hold to round-off).
    """
    u, v, t = _F(samples_u) ** m, _F(samples_v) ** m, _F(t)
    avg = lambda a: np.mean(a, axis=0)  # noqa: E731
    out = {"sum": float(np.max(np.abs(avg(u + v) - (avg(u) + avg(v))))),                      # Eq. (12.4)
           "constant": float(np.max(np.abs(avg(A * u) - A * avg(u)))),                        # Eq. (12.5)
           "d/dt": float(np.max(np.abs(avg(np.gradient(u, t, axis=-1)) - np.gradient(avg(u), t, axis=-1)))),  # Eq. (12.6)
           "integral dt": float(np.max(np.abs(avg(np.trapezoid(u, t, axis=-1)) - np.trapezoid(avg(u), t, axis=-1))))}  # (12.7)
    if u.ndim >= 3 and x is not None:
        xs, ax = _F(x), 1
    else:
        xs, ax = t, -1
    out["d/dx"] = float(np.max(np.abs(avg(np.gradient(u, xs, axis=ax)) - np.gradient(avg(u), xs, axis=ax - (ax > 0)))))  # (12.8)
    out["integral dx"] = float(np.max(np.abs(avg(np.trapezoid(u, xs, axis=ax)) - np.trapezoid(avg(u), xs, axis=ax - (ax > 0)))))  # (12.9)
    mean_field = np.broadcast_to(avg(u), u.shape)
    out["average of average"] = float(np.max(np.abs(avg(mean_field) - avg(u))))
    out["product"] = float(np.max(np.abs(avg(u * v) - avg(u) * avg(v))))
    return out


# ======================================================================================================================
# §12.4  correlations
# ======================================================================================================================
def correlation(a, b, axis: int = 0):
    """Correlation ``mean(a b)`` of two zero-mean variables over an ensemble of pairs.

    Book: §12.4, Eq. (12.12) ``R_ij = mean(u_i(x1, t1) u_j(x2, t2))``; Eq. (12.13) when b is the same variable at another
    point or time.  The inputs are used as given (no mean removed): the book's u_i are fluctuations.
    Parameters: a, b [units], pairs along ``axis``.  Returns [unit a·b].  Validation: V1 hand sum; R(a, a) = variance.
    Assumptions: inputs are fluctuations (zero mean): they are used as given, no mean is removed; the over-bar is estimated by the average over the samples supplied (an estimate with sampling error ~ N^{-1/2}).
    """
    return _S(np.mean(_F(a) * _F(b), axis=axis))  # Eq. (12.12)


def correlation_coefficient(a, b, axis: int = 0):
    """Correlation coefficient ``r = mean(ab)/sqrt(mean(a²) mean(b²))`` of two zero-mean variables, in [−1, 1].

    Book: §12.4, Eqs. (12.14) (cross) and (12.15) (auto); the bound is the Schwartz inequality (12.16) (on |mean(ab)|).
    Parameters: a, b fluctuations [units], pairs along ``axis``.  Returns r [–]; NaN if either variance is zero.
    A violation of |r| ≤ 1 beyond round-off raises (it would mean a bug); round-off overshoot is clipped.
    Validation: V1 r = ±1 for b = ±c a; invariant under scaling; V7 |r| ≤ 1 on random pairs; equals np.corrcoef for
    zero-mean data.
    Assumptions: both variances non-zero; the over-bar is estimated by the average over the samples supplied (an estimate with sampling error ~ N^{-1/2}).
    """
    a_, b_ = np.atleast_1d(_F(a)), np.atleast_1d(_F(b))
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.asarray(np.mean(a_ * b_, axis=axis)
                       / np.sqrt(np.mean(a_ ** 2, axis=axis) * np.mean(b_ ** 2, axis=axis)))  # Eq. (12.14)
    if np.any(np.abs(r[np.isfinite(r)]) > 1.0 + 1e-9):
        raise ArithmeticError("Schwartz inequality (12.16) violated beyond round-off")
    return _S(np.clip(r, -1.0, 1.0))


def correlated_pair(n: int, r: float, sigma_u: float = 1.0, sigma_v: float = 1.0, seed: int = 0):
    """n samples of a zero-mean Gaussian pair (u, v) with standard deviations σ_u, σ_v and correlation coefficient r.

    Book: §12.6, Fig. 12.8 (round cloud: mean(uv) = 0; tilted cloud: mean(uv) < 0).
    Parameters: n; r in [−1, 1]; sigma_u, sigma_v [unit]; seed of ``default_rng`` (default 0).  Returns (u, v), each (n,).
    Method: Cholesky factor of the 2 × 2 covariance: ``v = σ_v (r ξ1 + sqrt(1 − r²) ξ2)``.
    Validation: V1 sample r within 5 s.e. ((1 − r²)/sqrt(n)) of the request; r = ±1 exact.
    Assumptions: jointly Gaussian samples (Cholesky factor of the 2 x 2 covariance); |r| <= 1.
    """
    if not -1.0 <= r <= 1.0:
        raise ValueError("r must lie in [-1, 1]")
    rng = np.random.default_rng(seed)
    x1, x2 = rng.standard_normal(int(n)), rng.standard_normal(int(n))
    return sigma_u * x1, sigma_v * (r * x1 + np.sqrt(1.0 - r * r) * x2)


def autocorrelation(u, dt: float, max_lag: int | None = None, method: str = "fft", unbiased: bool = True,
                    demean: bool = True):
    """Autocorrelation function ``R(τ) = mean(u(t) u(t + τ))`` of one stationary record, for lags τ ≥ 0.

    Book: §12.4, Eq. (12.17) (time average standing in for the ensemble average; R is even, so τ ≥ 0 suffices).
    Parameters: u (n,) record [unit], uniformly sampled; dt sample spacing [s] (or [m]); max_lag number of lags
    (default n // 2); method "fft" (zero-padded FFT, O(n log n)) or "direct" (explicit sum, O(n·lags));
    unbiased: divide the lag-k sum by n − k (True) or by n (False, the biased estimate, smaller variance at large lag);
    demean: subtract the record mean first (the book's u is a fluctuation).
    Returns (lag [s], R [unit²]) with R[0] = variance of the record.
    Validation: V1 the two methods agree to 1e-10; u = cos ωt → ½cos ωτ; V3 OU record → σ²exp(−τ/τ_c) with error
    ∝ (record length)^{−1/2}.
    Assumptions: statistically stationary record, long against the correlation time; uniform sampling; time average used in place of the ensemble average (ergodicity).
    """
    x = _F(u).ravel()
    n = x.size
    if demean:
        x = x - x.mean()
    L = n // 2 if max_lag is None else int(max_lag)
    if not 0 <= L < n:
        raise ValueError("max_lag must satisfy 0 <= max_lag < len(u)")
    if method == "fft":
        nfft = 1 << int(np.ceil(np.log2(2 * n)))
        X = np.fft.rfft(x, nfft)
        acf = np.fft.irfft(X * np.conj(X), nfft)[: L + 1]
    elif method == "direct":
        acf = np.array([np.dot(x[: n - k], x[k:]) for k in range(L + 1)])
    else:
        raise ValueError("method must be 'fft' or 'direct'")
    norm = (n - np.arange(L + 1)) if unbiased else n
    return np.arange(L + 1) * float(dt), acf / norm  # Eq. (12.17)


def cross_correlation(u, v, dt: float, max_lag: int | None = None, demean: bool = True, unbiased: bool = True):
    """Cross-correlation ``R_uv(τ) = mean(u(t) v(t + τ))`` for lags of both signs.

    Book: §12.4, Eq. (12.17) with i ≠ j and Fig. 12.4: the maximum sits at the lag that aligns the two records; not even
    in τ, but ``R_uv(τ) = R_vu(−τ)``.
    Parameters: u, v (n,) records [units], same uniform sampling; dt [s]; max_lag (default n // 2).
    Returns (lag [s] from −max_lag·dt to +max_lag·dt, R [unit u·v]).  If v(t) = u(t − t0) the peak is at τ = +t0.
    Validation: V1 delayed copy peaks at the delay; R_uv(τ) = R_vu(−τ) to round-off; u = v reduces to
    :func:`autocorrelation`.
    Assumptions: statistically stationary record, long against the correlation time; uniform sampling; both records share the sampling interval.
    """
    a, b = _F(u).ravel(), _F(v).ravel()
    if a.size != b.size:
        raise ValueError("u and v must have the same length")
    n = a.size
    if demean:
        a, b = a - a.mean(), b - b.mean()
    L = n // 2 if max_lag is None else int(max_lag)
    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    c = np.fft.irfft(np.conj(np.fft.rfft(a, nfft)) * np.fft.rfft(b, nfft), nfft)  # c[k] = Σ a[j] b[j + k]
    k = np.arange(-L, L + 1)
    vals = c[k % nfft]
    norm = (n - np.abs(k)) if unbiased else n
    return k * float(dt), vals / norm


def correlation_time(lag, r):
    """Correlation time t_c: the first zero crossing of the autocorrelation coefficient (linear interpolation).

    Book: §12.4, Fig. 12.5 (t_c marked where r first reaches zero).
    Parameters: lag (n,) [s] ≥ 0; r (n,) correlation (coefficient or function).  Returns t_c [s]; ``inf`` if r never
    reaches zero on the given lags.
    Validation: V1 r = cos ωτ·envelope: first zero at π/2ω within one sample.
    Assumptions: r is sampled finely enough for linear interpolation of its first zero.
    """
    lag, r = _F(lag), _F(r)
    idx = np.where(r <= 0.0)[0]
    if idx.size == 0:
        return float("inf")
    k = int(idx[0])
    if k == 0:
        return float(lag[0])
    return float(lag[k - 1] + (lag[k] - lag[k - 1]) * r[k - 1] / (r[k - 1] - r[k]))


def integral_scale(lag, r, upto: str = "first_zero") -> float:
    """Integral scale ``Λ = ∫_0^∞ r dτ`` (time [s] or length [m]): the width of the rectangle of unit height and equal area.

    Book: §12.4, Eq. (12.18) (temporal, the "memory" of the turbulence); §12.6, Eq. (12.39) (Λ_f, Λ_g with r for τ).
    Parameters: lag (n,) ≥ 0; r (n,) correlation function or coefficient (normalised by r[0] internally);
    upto "first_zero" (integrate to the first zero crossing — robust for measured, noisy tails) or "all" (the whole
    range given — right for an exact, fully decayed correlation).
    **Warning: ``upto="first_zero"`` (the default) is not the Λ = ∫_0^∞ r dτ of (12.18) when r has a negative lobe** —
    it leaves the negative area out and so over-estimates Λ.  The transverse correlation g(r) of isotropic
    turbulence always has one (g = f + (r/2) f′ (12.41) integrates to Λ_g = Λ_f/2 only *with* its negative lobe);
    so does a correlation with an oscillating tail.  Use ``upto="all"`` on a range over which r has decayed for
    those; "first_zero" is the practical estimate for a positive, noisy r (f(r), a measured autocorrelation).
    Returns Λ in the unit of ``lag``.  Trapezoid rule.
    Validation: V1 exponential → τ_c; Gaussian exp(−τ²/t_c²) → sqrt(π) t_c/2 (1e-6 on fine grids).
    Assumptions: r decays within the range given; trapezoid rule; the result depends on where the integral is stopped (stated by 'upto').
    """
    lag, r = _F(lag), _F(r)
    r = r / r[0]
    if upto == "all":
        return float(np.trapezoid(r, lag))  # Eq. (12.18)
    if upto != "first_zero":
        raise ValueError("upto must be 'first_zero' or 'all'")
    tc = correlation_time(lag, r)
    if not np.isfinite(tc):
        return float(np.trapezoid(r, lag))
    k = int(np.searchsorted(lag, tc))
    return float(np.trapezoid(np.append(r[:k], 0.0), np.append(lag[:k], tc)))  # Eq. (12.18), truncated at t_c


def effective_samples(record_length, t_c):
    """Equivalent number of independent samples in a record of length Δt: ``N ≈ Δt/t_c``.

    Book: §12.4 (the relation after (12.18)).  Parameters: record_length Δt [s]; t_c correlation time [s] (the book's
    choice; using 2Λ_t instead is the usual statistical estimate).  Returns N [–].  Scalar-callable.
    Assumptions: samples further apart than the correlation time are treated as independent (an order-of-magnitude count).
    Validation: V1 arithmetic; the scatter of time_average over a record matches sigma/sqrt(N_eff) within a factor 2.
    """
    return _S(_F(record_length) / _F(t_c))


def taylor_microscale(lag, r, fit_points: int = 5) -> float:
    """Taylor microscale from the curvature of the correlation peak: ``λ² = −2/[d²r/dτ²]_{τ=0}``.

    Book: §12.4, Eq. (12.19) (temporal, λ_t); §12.6, Eq. (12.39) (λ_f, λ_g).  λ is where the osculating parabola
    ``1 − τ²/λ²`` crosses zero.
    Parameters: lag (n,) ≥ 0 starting at 0; r (n,) correlation (normalised by r[0] internally); fit_points number of
    leading samples (≥ 3) used in the least-squares even fit ``r ≈ 1 + a τ² + b τ⁴`` (b only when fit_points ≥ 4).
    Returns λ = sqrt(−1/a) in the unit of ``lag``; NaN if the fitted curvature is not negative.
    Assumptions: r is smooth at the origin — an exponential (OU) correlation has a cusp and **no** Taylor microscale;
    the fit then returns a value that shrinks with the sample spacing instead of converging.
    Validation: V1 Gaussian exp(−τ²/t_c²) → λ = t_c (1e-6 on fine grids).
    """
    lag, r = _F(lag), _F(r)
    n = max(3, int(fit_points))
    x, y = lag[:n], r[:n] / r[0] - 1.0
    cols = [x ** 2, x ** 4] if n >= 4 else [x ** 2]
    a = np.linalg.lstsq(np.stack(cols, axis=1), y, rcond=None)[0][0]
    return float(np.sqrt(-1.0 / a)) if a < 0 else float("nan")  # Eq. (12.19): r'' = 2a


def integral_scale_from_spectrum(S0, variance):
    """Integral scale from the zero-frequency value of the two-sided spectrum: ``Λ = π S(0)/variance``.

    Book: §12.4, the relation after (12.22): ``S_e(0) = (variance/π) Λ_t``.  S0 [unit² s] two-sided; variance [unit²].
    Returns Lambda [unit of the lag: s or m] = pi S(0)/variance; a float for float input.
    Assumptions: two-sided spectrum in the book normalisation (12.20), so S(0) = variance * Lambda/pi.
    Validation: V1 exact for the three closed-form pairs of correlation_spectrum_pair.
    """
    return _S(np.pi * _F(S0) / _F(variance))


# ======================================================================================================================
# §12.4  spectra
# ======================================================================================================================
def _cos_transform(x: np.ndarray, y: np.ndarray, w: np.ndarray, rule: str = "filon") -> np.ndarray:
    """∫ y(x) cos(w x) dx over the range of x: "filon" = exact for the piecewise-linear interpolant of y (any w);
    "trapezoid" = plain trapezoid rule on y·cos (spectrally accurate for smooth, decayed y while w·h << π)."""
    w = np.atleast_1d(w).astype(float).ravel()
    out = np.empty(w.shape)
    h = np.diff(x)
    slope = np.diff(y) / h
    if rule == "trapezoid":
        small = np.ones(w.shape, dtype=bool)
    elif rule == "filon":
        small = np.abs(w) * np.max(h) < 1e-4
    else:
        raise ValueError("rule must be 'filon' or 'trapezoid'")
    if np.any(small):
        ws = w[small][:, None]
        out[small] = np.trapezoid(y[None, :] * np.cos(ws * x[None, :]), x, axis=1)
    if np.any(~small):
        wb = w[~small][:, None]
        s, c = np.sin(wb * x[None, :]), np.cos(wb * x[None, :])
        term1 = (y[-1] * s[:, -1] - y[0] * s[:, 0]) / wb[:, 0]
        term2 = np.sum(slope[None, :] * (c[:, 1:] - c[:, :-1]), axis=1) / wb[:, 0] ** 2
        out[~small] = term1 + term2
    return out


def spectrum_from_correlation(lag, R, omega, rule: str = "filon"):
    """Two-sided spectrum as the Fourier transform of the correlation: ``S(ω) = (1/2π) ∫ R(τ) e^{−iωτ} dτ``.

    Book: §12.4, Eq. (12.20) (frequency, S_e(ω)); §12.6, Eq. (12.45) (wavenumber, S_11(k_1), with r_1 for τ and k_1 for ω).
    Because R is real and even the transform is real: ``S(ω) = (1/π) ∫_0^∞ R cos ωτ dτ``.
    Parameters: lag (n,) [s] (or [m]) — either τ ≥ 0 only (the even half is implied) or a symmetric range; R (n,)
    [unit²]; omega scalar or array [rad/s] (or [rad/m]).
    Returns S [unit² s] (or [unit² m]), **two-sided, book normalisation** (∫_{−∞}^{∞} S dω = R(0)).
    Method: ``rule="filon"`` (default) integrates the piecewise-linear interpolant of R against cos ωτ exactly (no
    aliasing at large ω; second order in the lag spacing); ``rule="trapezoid"`` is the plain trapezoid rule on R cos ωτ
    (spectrally accurate for a smooth correlation as long as ω·Δτ << π).  R must have decayed at the end of the lag range.
    Validation: V1 exponential R ↔ Lorentzian (σ²τ_c/π)/(1 + ω²τ_c²); Gaussian ↔ Gaussian; S(0) = σ²Λ/π; V2 round trip
    with :func:`correlation_from_spectrum`.
    Assumptions: R even in the lag and negligible beyond the last lag given (otherwise truncation ripples appear); trapezoid quadrature.
    """
    lag, R = _F(lag), _F(R)
    factor = 1.0 / np.pi if lag[0] >= 0.0 else 1.0 / (2.0 * np.pi)
    S = factor * _cos_transform(lag, R, _F(omega), rule)  # Eq. (12.20)
    return _S(S.reshape(np.shape(omega)))


def correlation_from_spectrum(omega, S, lag, rule: str = "filon"):
    """Correlation as the inverse transform of the two-sided spectrum: ``R(τ) = ∫ S(ω) e^{+iωτ} dω``.

    Book: §12.4, Eq. (12.21).  Parameters: omega (n,) [rad/s] — ω ≥ 0 only (even half implied) or a symmetric range;
    S (n,) two-sided spectrum [unit² s]; lag scalar or array [s].  Returns R [unit²]; R(0) is the variance (12.22).
    Method and validation as :func:`spectrum_from_correlation` (note there is no 1/2π in this member of the pair).
    Assumptions: S even in omega and negligible beyond the last frequency given; trapezoid quadrature.
    Validation: V1 round trip with spectrum_from_correlation for the exponential and Gaussian pairs (error set by the truncation).
    """
    omega, S = _F(omega), _F(S)
    factor = 2.0 if omega[0] >= 0.0 else 1.0
    R = factor * _cos_transform(omega, S, _F(lag), rule)  # Eq. (12.21)
    return _S(R.reshape(np.shape(lag)))


def spectrum_variance(omega, S, two_sided: bool = True, weights=None) -> float:
    """Variance carried by a spectrum: ``∫_{−∞}^{∞} S dω`` (two-sided) or ``∫_0^∞ S dω`` (one-sided).

    Book: §12.4, Eq. (12.22); §12.7, Eq. (12.55) ``∫_{−∞}^{∞} S_11 dk_1 = mean(u_1²) = ∫_0^∞ 2 S_11 dk_1``.
    Parameters: omega (n,) [rad/s] or k [rad/m]; S (n,) [unit² s]; two_sided: S is the two-sided density (book
    normalisation) — if the samples cover ω ≥ 0 only, the even half is added (factor 2); False: S is one-sided;
    weights: the bin weights returned by ``periodogram(..., return_weights=True)`` — then the exact discrete sum
    ``Σ w S`` is returned (Parseval to round-off) and the other arguments are ignored.
    Returns the variance [unit²] (trapezoid rule unless ``weights`` is given).
    Validation: V1 Lorentzian and Gaussian spectra integrate to σ² (1e-8 on wide fine grids); factor-2 mutant fails.
    Assumptions: the spectrum is resolved on the grid and has decayed at its end.
    """
    S = _F(S)
    if weights is not None:
        return float(np.sum(_F(weights) * S))
    omega = _F(omega)
    integral = float(np.trapezoid(S, omega))
    return 2.0 * integral if (two_sided and omega[0] >= 0.0) else integral  # Eq. (12.22) / (12.55)


def periodogram(u, d: float, two_sided: bool = True, segments: int = 1, window: str = "boxcar",
                one_sided: bool | None = None, demean: bool = True, return_weights: bool = False):
    """Spectrum of one record by the finite-window Fourier transform (periodogram), in the book's normalisation.

    Book: §12.4, the alternative to (12.20) quoted there (Exercise 12.8):
    ``S(ω) = lim (1/2πT) |∫_{−T/2}^{T/2} u(t) e^{−iωt} dt|²``; §12.6 (12.45) for a spatial record (pass d = dx).
    Parameters: u (n,) record [unit], uniform spacing d [s] (or [m]); two_sided True → the two-sided density S (book),
    False → one-sided 2S (``one_sided=True`` is a synonym and wins when given); segments: average this many equal,
    non-overlapping segments (Bartlett; reduces the scatter of each estimate from 100 % to 100 %/sqrt(segments));
    window "boxcar" or "hann"; demean: subtract the record mean first; return_weights: also return the bin weights w.
    Returns (omega [rad/s] or k [rad/m] ≥ 0, S [unit² s]) (and w): the angular frequencies ``2π j/(m d)`` of a segment
    of m samples.  ``Σ w S`` is the variance exactly: for one boxcar segment this is Parseval's theorem; with several
    segments or a taper the estimate is rescaled so that the sum equals the variance of the samples used (stated, so it
    is a normalisation there, not evidence).  ``spectrum_variance(omega, S)`` (trapezoid) agrees to O(1/m).
    Validation: V4 Parseval to 1e-12; V1 OU record → Lorentzian within 5 s.e. per bin after averaging; a sine puts its
    variance in the right bin; agreement with ``scipy.signal.welch`` after the unit conversion S(ω) = P(f)/(4π).
    Assumptions: statistically stationary record, long against the correlation time; uniform sampling; one record gives a noisy estimate (relative scatter ~ 1 per bin, reduced as 1/sqrt(segments)).
    """
    x = _F(u).ravel()
    if one_sided is not None:
        two_sided = not one_sided
    nseg = int(segments)
    m = x.size // nseg
    if m < 2:
        raise ValueError("record too short for the requested number of segments")
    x = x[: m * nseg]
    if demean:
        x = x - x.mean()
    variance = float(np.mean(x ** 2))
    if window == "boxcar":
        win = np.ones(m)
    elif window == "hann":
        win = 0.5 - 0.5 * np.cos(2.0 * np.pi * np.arange(m) / m)
    else:
        raise ValueError("window must be 'boxcar' or 'hann'")
    X = np.fft.rfft(x.reshape(nseg, m) * win, axis=1)
    T = m * float(d)
    # (1/2πT)|∫u e^{-iωt}dt|² with ∫ ≈ d·X; a taper is compensated by mean(win²)
    S = np.mean(np.abs(X) ** 2, axis=0) * float(d) ** 2 / (2.0 * np.pi * T * np.mean(win ** 2))
    omega = 2.0 * np.pi * np.fft.rfftfreq(m, d=float(d))
    dw = 2.0 * np.pi / T
    w = np.full(omega.size, 2.0 * dw)   # ±ω bins
    w[0] = dw
    if m % 2 == 0:
        w[-1] = dw                      # the Nyquist bin is its own mirror image
    if (nseg > 1 or window != "boxcar") and variance > 0.0:
        # DEVIATION: the book's formula is a limit T → ∞; for segment averages and tapers the finite-record estimate is
        # rescaled so that Σ w S equals the sample variance exactly (a normalisation, not evidence of Parseval).
        S = S * variance / float(np.sum(w * S))
    if not two_sided:
        S, w = 2.0 * S, 0.5 * w
    return (omega, S, w) if return_weights else (omega, S)


def taylor_frozen(t, u, U0: float):
    """Taylor's frozen-turbulence hypothesis: turn a time record into a space record, ``x = U0 t``.

    Book: §12.4 (after (12.23)): a probe moving at speed U0 (or a mean flow sweeping the field past a fixed probe) sees
    a pattern that does not change while it passes, accurate as u_rms/U0 → 0.
    Parameters: t (n,) [s]; u (n,) [unit]; U0 probe or sweeping speed [m/s] (> 0).  Returns (x [m], u [unit]).
    Validation: V7 a frozen sinusoid of wavenumber k sampled at U0 has frequency ω = k U0.
    Assumptions: Taylor's frozen-turbulence hypothesis: u_rms << U0, so the pattern does not evolve while it passes the probe.
    """
    return _F(t) * float(U0), _F(u)


def frequency_to_wavenumber_spectrum(omega, S, U0: float):
    """Convert a frequency spectrum to a streamwise wavenumber spectrum under the frozen hypothesis.

    Book: §12.4 (Taylor's hypothesis) with §12.6, Eq. (12.45): ``k_1 = ω/U0`` and ``S_11(k_1) = U0 S_e(ω)`` (the Jacobian
    of the change of variable, ours), so that ``∫S_11 dk_1 = ∫S_e dω`` — the variance is preserved.
    Parameters: omega [rad/s]; S [unit² s]; U0 [m/s] > 0.  Returns (k_1 [rad/m], S_11 [unit² m]).  Scalar-callable.
    Validation: V1 variance preserved to 1e-12.
    Assumptions: Taylor's frozen-turbulence hypothesis (u_rms << U0); the Jacobian U0 keeps the integral equal to the variance.
    """
    return _S(_F(omega) / float(U0)), _S(_F(S) * float(U0))


def correlation_spectrum_pair(kind: str, sigma: float, tau_c: float, omega0: float = 0.0) -> dict:
    """Closed-form correlation ↔ spectrum pairs (book normalisation) with their exact scales.

    Book: §12.4, Eqs. (12.17)–(12.22) and the relation ``S_e(0) = variance·Λ_t/π``.  The three pairs are ours:

    * ``"exponential"``:  r = exp(−|τ|/τ_c);  S = (σ²τ_c/π)/(1 + ω²τ_c²);  Λ = τ_c;  no Taylor microscale (cusp).
    * ``"gaussian"``:     r = exp(−τ²/τ_c²);  S = (σ²τ_c/2√π) exp(−ω²τ_c²/4);  Λ = √π τ_c/2;  λ = τ_c.
    * ``"damped_cosine"``: r = exp(−|τ|/τ_c) cos ω0τ;  S = (σ²τ_c/2π)[1/(1 + (ω−ω0)²τ_c²) + 1/(1 + (ω+ω0)²τ_c²)];
      Λ = τ_c/(1 + ω0²τ_c²);  first zero t_c = π/(2ω0);  no Taylor microscale (cusp).

    Parameters: kind; sigma standard deviation [unit]; tau_c [s] > 0; omega0 [rad/s] (damped cosine only).
    Returns dict(r = callable τ → coefficient, R = callable τ → σ²r, S = callable ω → two-sided spectrum [unit² s],
    Lambda_t [s], lambda_t [s] (NaN where it does not exist), S0 [unit² s], t_c [s] (inf where r has no zero),
    variance).  All callables accept scalars and arrays.
    Validation: V2 sympy transform of each r; V1 ∫S dω = σ², S(0) = σ²Λ/π.
    Assumptions: stationary signal with the stated analytic correlation; book normalisation (two-sided, 1/2pi forward).
    """
    s2 = float(sigma) ** 2
    tc = float(tau_c)
    if tc <= 0:
        raise ValueError("tau_c must be > 0")
    if kind == "exponential":
        r = lambda tau: _S(np.exp(-np.abs(_F(tau)) / tc))  # noqa: E731
        S = lambda w: _S(s2 * tc / np.pi / (1.0 + (_F(w) * tc) ** 2))  # noqa: E731
        Lam, lam, tzero = tc, float("nan"), float("inf")
    elif kind == "gaussian":
        r = lambda tau: _S(np.exp(-(_F(tau) / tc) ** 2))  # noqa: E731
        S = lambda w: _S(s2 * tc / (2.0 * np.sqrt(np.pi)) * np.exp(-(_F(w) * tc) ** 2 / 4.0))  # noqa: E731
        Lam, lam, tzero = np.sqrt(np.pi) * tc / 2.0, tc, float("inf")
    elif kind == "damped_cosine":
        w0 = float(omega0)
        r = lambda tau: _S(np.exp(-np.abs(_F(tau)) / tc) * np.cos(w0 * _F(tau)))  # noqa: E731
        S = lambda w: _S(s2 * tc / (2.0 * np.pi) * (1.0 / (1.0 + ((_F(w) - w0) * tc) ** 2)  # noqa: E731
                                                     + 1.0 / (1.0 + ((_F(w) + w0) * tc) ** 2)))
        Lam, lam = tc / (1.0 + (w0 * tc) ** 2), float("nan")
        tzero = np.pi / (2.0 * abs(w0)) if w0 != 0.0 else float("inf")
    else:
        raise ValueError("kind must be 'exponential', 'gaussian' or 'damped_cosine'")
    return dict(kind=kind, r=r, R=lambda tau: _S(s2 * _F(r(tau))), S=S, Lambda_t=float(Lam), lambda_t=float(lam),
                S0=float(s2 * Lam / np.pi), t_c=float(tzero), variance=s2)


def smooth_signal(n: int, dt: float, spectrum: Callable, seed: int = 0) -> np.ndarray:
    """Random Gaussian record with a prescribed smooth two-sided spectrum S(ω) (so its Taylor microscale exists).

    Book: §12.4, Eq. (12.19) needs a correlation that is smooth at the origin, which an Ornstein–Uhlenbeck record is not.
    Parameters: n samples; dt spacing [s]; spectrum callable ω [rad/s] → two-sided S [unit² s] (e.g. the ``S`` of
    ``correlation_spectrum_pair("gaussian", ...)``); seed of ``default_rng`` (default 0).
    Returns u (n,) [unit], zero mean, periodic with period n·dt; its expected variance is Σ S(ω_j) Δω over the resolved
    frequencies (= ∫S dω when S has decayed below the Nyquist frequency π/dt).
    Method: white noise → rFFT → multiply by sqrt(S(ω_j) Δω n) → inverse rFFT.
    Validation: V1 its autocorrelation and periodogram return the prescribed pair within 5 s.e.
    Assumptions: Gaussian random-phase signal, periodic over the record; the prescribed spectrum must decay fast enough for the Taylor microscale to exist.
    """
    rng = np.random.default_rng(seed)
    W = np.fft.rfft(rng.standard_normal(int(n)))
    omega = 2.0 * np.pi * np.fft.rfftfreq(int(n), d=float(dt))
    dw = 2.0 * np.pi / (int(n) * float(dt))
    amp = np.sqrt(np.maximum(_F(spectrum(omega)), 0.0) * dw * int(n))
    amp[0] = 0.0
    return np.fft.irfft(W * amp, int(n))


# ======================================================================================================================
# §12.1, §12.6  synthetic fields and spatial correlations
# ======================================================================================================================
def _wavenumbers(n: int, L: float, dim: int):
    k1 = 2.0 * np.pi * np.fft.fftfreq(n, d=L / n)
    ks = np.meshgrid(*([k1] * dim), indexing="ij")
    return ks, np.sqrt(sum(k ** 2 for k in ks))


def synthetic_solenoidal_field(n: int, L: float, spectrum: Callable, seed: int = 0, dim: int = 2):
    """Random divergence-free periodic velocity field with a prescribed energy spectrum E(K) — **kinematic, no cascade**.

    Book: §12.1 (a turbulent velocity field obeys ∇·u = 0, a random vector field need not), §12.6 (isotropic fields for
    Eqs. 12.36–12.43).
    Parameters: n grid points per side; L box size [m]; spectrum callable K [rad/m] → E(K) [m³/s²] with
    ``∫_0^∞ E dK = ½ mean(u_i u_i)``; seed of ``default_rng`` (default 0); dim 2 or 3.
    Returns a tuple of ``dim`` arrays: dim = 2 → (u, v) with shape (n, n) = [y, x]; dim = 3 → (u, v, w) with shape
    (n, n, n) = [z, y, x] (x on the last axis, the ``core.grids`` layout).  Node values, spacing L/n.
    Method: dim 2 — random stream function ψ̂(k) = sqrt(2 e_k)/K · noise, u = ∂ψ/∂y, v = −∂ψ/∂x; dim 3 — Gaussian
    Fourier vectors projected perpendicular to k, û = (I − k k/K²) ŵ sqrt(e_k); e_k = E(K)(2π/L)^d/(S_d K^{d−1}) is the
    energy of one mode (S_2 = 2π, S_3 = 4π), so the shell-summed energy reproduces E(K) in expectation.
    **Use dim = 3 for anything about (12.41)–(12.43) or the cascade**: in 2-D g = d(rf)/dr, not f + (r/2)f'.
    Validation: V4 spectral divergence at round-off; V1 shell spectrum → E(K) and kinetic energy → Σ e_k (5 s.e.);
    V7 statistics invariant under a 90° rotation.
    Assumptions: kinematic field: prescribed spectrum, random phases, periodic box - it satisfies continuity but no dynamics (no cascade, no skewness).
    """
    if dim not in (2, 3):
        raise ValueError("dim must be 2 or 3")
    rng = np.random.default_rng(seed)
    ks, K = _wavenumbers(int(n), float(L), dim)
    Ksafe = np.where(K > 0, K, 1.0)
    surf = 2.0 * np.pi if dim == 2 else 4.0 * np.pi
    e_k = np.where(K > 0, _F(spectrum(Ksafe)) * (2.0 * np.pi / L) ** dim / (surf * Ksafe ** (dim - 1)), 0.0)
    # remove the unpaired Nyquist modes (they cannot be both real and solenoidal)
    nyq = np.zeros(K.shape, dtype=bool)
    if n % 2 == 0:
        for k in ks:
            nyq |= np.isclose(np.abs(k), np.pi * n / L)
    e_k = np.where(nyq, 0.0, e_k)
    norm = float(n) ** (dim / 2.0)
    if dim == 2:
        ky, kx = ks
        psi_hat = np.fft.fft2(rng.standard_normal((n, n))) / norm * np.sqrt(2.0 * e_k) / Ksafe
        u = np.real(np.fft.ifft2(1j * ky * psi_hat)) * norm ** 2
        v = np.real(np.fft.ifft2(-1j * kx * psi_hat)) * norm ** 2
        return u, v
    kz, ky, kx = ks
    w_hat = [np.fft.fftn(rng.standard_normal((n, n, n))) / norm for _ in range(3)]  # order: x, y, z components
    kvec = (kx, ky, kz)
    kdotw = sum(k * w for k, w in zip(kvec, w_hat))
    comps = []
    for k, w in zip(kvec, w_hat):
        u_hat = (w - k * kdotw / Ksafe ** 2) * np.sqrt(e_k)
        comps.append(np.real(np.fft.ifftn(u_hat)) * norm ** 2)
    return tuple(comps)


def white_noise_field(n: int, seed: int = 0, dim: int = 2):
    """Independent Gaussian noise in every velocity component at every node — a random vector field that is *not* a flow.

    Book: §12.1 (turbulence is not mere randomness: it satisfies the equations of motion, here ∇·u = 0).
    Parameters: n points per side; seed of ``default_rng`` (default 0); dim 2 or 3.  Returns ``dim`` arrays of unit
    variance, layout as :func:`synthetic_solenoidal_field`.  Validation: V4 its divergence is O(1/dx), not round-off.
    Assumptions: independent Gaussian values at every node: not divergence-free, no spatial correlation (a deliberate counter-example).
    """
    rng = np.random.default_rng(seed)
    return tuple(rng.standard_normal((int(n),) * dim) for _ in range(dim))


def shell_spectrum(components, L: float, n_bins: int | None = None):
    """Energy spectrum E(K) of a periodic field by summing ½|û|² over spherical (circular) wavenumber shells.

    Book: §12.7, the normalisation ``ē = ∫_0^∞ S(K) dK`` stated after (12.55).
    Parameters: components tuple of 2 or 3 periodic arrays [m/s] (layout of :func:`synthetic_solenoidal_field`);
    L box size [m]; n_bins shells of width 2π/L (default n // 2).
    Returns (K [rad/m] shell centres, E [m³/s²]) with Σ E ΔK = ½ mean(u_i u_i) over the shells kept.
    Validation: V1 returns the spectrum given to :func:`synthetic_solenoidal_field` within sampling error.
    Assumptions: periodic box, isotropic binning in |K| (meaningful for a statistically isotropic field).
    """
    comps = [_F(c) for c in components]
    dim, n = comps[0].ndim, comps[0].shape[0]
    _, K = _wavenumbers(n, float(L), dim)
    e = sum(0.5 * np.abs(np.fft.fftn(c) / n ** dim) ** 2 for c in comps)
    dk = 2.0 * np.pi / L
    nb = n // 2 if n_bins is None else int(n_bins)
    idx = np.rint(K / dk).astype(int).ravel()
    E = np.bincount(idx, weights=e.ravel(), minlength=nb + 1)[1: nb + 1] / dk
    return dk * np.arange(1, nb + 1), E


def spatial_correlation(a, b, dx: float, axis: int = -1):
    """Two-point correlation ``R(r) = mean(a(x) b(x + r e_axis))`` of periodic fields, averaged over the whole field.

    Book: §12.4, Eq. (12.23) ``R_ij(r) = mean(u_i(x) u_j(x + r))`` (homogeneous field; **volume average**).
    Parameters: a, b periodic arrays of equal shape [units]; dx grid spacing [m] along ``axis``.
    Returns (r [m] for 0 … n/2 points, R [unit a·b]).  Periodic FFT correlation, exact for the discrete field.
    Validation: V1 against a direct shifted product; R(0) = mean(ab).
    Assumptions: homogeneous along the chosen axis and periodic; spatial average in place of the ensemble average.
    """
    A, B = _F(a), _F(b)
    n = A.shape[axis]
    c = np.fft.irfft(np.conj(np.fft.rfft(A, axis=axis)) * np.fft.rfft(B, axis=axis), n, axis=axis) / n
    c = np.moveaxis(c, axis, -1).reshape(-1, n).mean(axis=0)
    m = n // 2 + 1
    return np.arange(m) * float(dx), c[:m]  # Eq. (12.23)


def longitudinal_transverse_correlation(u, v, dx: float, w=None):
    """Longitudinal f(r) and transverse g(r) correlation coefficients of a periodic, homogeneous, isotropic field.

    Book: §12.6, Eq. (12.38): f uses the velocity component parallel to the separation r, g a component perpendicular
    to it; f(0) = g(0) = 1.  **Volume average**, then an average over the equivalent directions.
    Parameters: u, v (and w) periodic arrays, x on the last axis, y on the second-to-last (z on the third-to-last);
    dx grid spacing [m] (same in every direction).
    Returns (r [m], f, g).  2-D field: f from u along x and v along y, g from u along y and v along x.  3-D: all three
    longitudinal and all six transverse pairs.
    Remember: for an incompressible field g = f + (r/2) f' in 3-D (12.41) but g = d(r f)/dr in 2-D.
    Validation: V1 f(0) = g(0) = 1; 3-D field of 64³ obeys (12.41) within sampling error.
    Assumptions: homogeneous periodic field; separation along x; spatial average in place of the ensemble average.
    """
    comps = [_F(u), _F(v)] + ([] if w is None else [_F(w)])
    dim = comps[0].ndim
    if dim != len(comps):
        raise ValueError("pass (u, v) for a 2-D field and (u, v, w) for a 3-D field")
    axes = [dim - 1 - i for i in range(dim)]  # component i is along array axis axes[i]
    f_list, g_list = [], []
    r = None
    for i, c in enumerate(comps):
        for j, ax in enumerate(axes):
            r, R = spatial_correlation(c, c, dx, axis=ax)
            (f_list if i == j else g_list).append(R)
    u2 = np.mean([np.mean(c ** 2) for c in comps])
    return r, np.mean(f_list, axis=0) / u2, np.mean(g_list, axis=0) / u2  # Eq. (12.38)


# ======================================================================================================================
# §12.5  Reynolds stress
# ======================================================================================================================
def velocity_covariance(u_samples) -> np.ndarray:
    """Velocity correlation tensor ``mean(u_i u_j)`` of an ensemble of velocity vectors (means removed).

    Book: §12.5, the tensor in (12.30).  Parameters: u_samples (N, d) [m/s], d = 2 or 3 components per sample.
    Returns (d, d) symmetric, positive semi-definite array [m²/s²]; its trace is 2ē.
    Validation: V1 against explicit sums and ``np.cov(bias=True)``; rotation of the samples by C gives CᵀRC.
    Assumptions: the over-bar is estimated by the average over the samples supplied (an estimate with sampling error ~ N^{-1/2}); the sample mean of each component is removed.
    """
    s = _F(u_samples)
    if s.ndim != 2:
        raise ValueError("u_samples must have shape (N, components)")
    f = s - s.mean(axis=0)
    return f.T @ f / s.shape[0]


def reynolds_stress(u_samples, rho0: float = 1.0) -> np.ndarray:
    """Reynolds stress tensor ``−ρ0 mean(u_i u_j)`` from an ensemble of velocity vectors.

    Book: §12.5, Eq. (12.30) and the definition after it: symmetric, six independent components; the diagonal entries
    are normal stresses (they add to the mean pressure), the off-diagonal ones shear stresses.
    Parameters: u_samples (N, d) total or fluctuating velocities [m/s] (the ensemble mean is removed); rho0 density
    [kg/m³] (default 1 → the kinematic stress, i.e. minus the covariance, [m²/s²]).
    Returns (d, d) array [Pa] (or [m²/s²] for rho0 = 1).  **Sign:** a positive mean shear dU/dy with mean(uv) < 0 gives a
    positive shear stress component, like the viscous stress μ dU/dy.  Use :func:`velocity_covariance` for mean(u_i u_j).
    Validation: V1 symmetric; −trace/ρ0 = 2ē; sign test with a correlated pair.
    Assumptions: the over-bar is estimated by the average over the samples supplied (an estimate with sampling error ~ N^{-1/2}); constant density rho0 (Boussinesq).
    """
    return -float(rho0) * velocity_covariance(u_samples)  # Eq. (12.30): the last term of the mean stress


def turbulent_kinetic_energy(u_samples) -> float:
    """Turbulent kinetic energy per unit mass ``ē = ½ mean(u_i u_i)`` (half the trace of the covariance).

    Book: §12.6, ``R_ii(0) = mean(u_i u_i) = 2ē``.  Parameters: u_samples (N, d) [m/s].  Returns ē [m²/s²] (sum over the
    components supplied: pass all three for the full ē).  Validation: V1 isotropic Gaussian samples → (d/2)σ².
    Assumptions: the over-bar is estimated by the average over the samples supplied (an estimate with sampling error ~ N^{-1/2}); all velocity components supplied are summed (2 components give a 2-D energy).
    """
    return float(0.5 * np.trace(velocity_covariance(u_samples)))


def anisotropy_tensor(uu) -> np.ndarray:
    """Normalised anisotropy tensor ``b_ij = mean(u_i u_j)/(2ē) − δ_ij/3`` (ours; zero for isotropic turbulence).

    Book: §12.6 (isotropy: equal normal stresses, no shear stress).  Parameters: uu (3, 3) covariance [m²/s²].
    Returns (3, 3) traceless symmetric tensor [–].  Validation: V1 trace 0; isotropic input → 0.
    Assumptions: uu is a full 3 x 3 covariance with non-zero trace.
    """
    R = _F(uu)
    if R.shape != (3, 3):
        raise ValueError("anisotropy_tensor needs the full 3 x 3 covariance")
    return R / np.trace(R) - np.eye(3) / 3.0
