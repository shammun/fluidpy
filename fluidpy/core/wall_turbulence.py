"""Mean-flow summaries of wall-bounded turbulence: wall units, the law of the wall, the logarithmic law, composite profiles,
total stress, friction laws, rough walls and the stratified surface layer.

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e, §12.9 (Eqs. 12.76–12.93 and the unnumbered profile and friction
formulas of that section) and the log-linear wind profile of §12.11.  Written as a core module because every later chapter
with a wall or a ground needs the same summaries (surface layer and bulk drag in ch13, turbulent skin friction in ch14).

Conventions (fixed for the whole module)
----------------------------------------
* ``y`` (or ``z``) is the distance from the wall [m].  ``h`` is the **full** channel height (walls at y = 0 and y = h),
  ``delta`` the half-height h/2, the pipe radius or the boundary-layer thickness, ``d`` a pipe diameter, ``a`` its radius.
* Wall units: ``u_star = sqrt(tau0/rho)`` (12.81), ``l_nu = nu/u_star``, ``y+ = y/l_nu``, ``U+ = U/u_star``,
  ``Re_tau = delta+ = delta u_star/nu``.
* ``kappa`` is the **von Kármán constant** (never the thermal diffusivity, which fluidpy calls ``kappa_th``).
* **Empirical constants are never silent defaults.**  ``kappa``, the additive constant ``B``, the wake strength ``Pi`` and a
  roughness length are keyword arguments the caller must supply; :data:`LOG_LAW_CONSTANTS` offers cited public presets to
  pass explicitly (``c = LOG_LAW_CONSTANTS["classical"]; log_law(y, kappa=c["kappa"], B=c["B"])``).  The values the book
  tabulates per flow type are book-quoted numbers and live only in the private ``tests/book_values_ch12.json``.
* Stratified surface layer: ``L_M`` is the Monin–Obukhov length with the sign of (12.110) — positive stable, negative
  unstable, ``inf`` neutral.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from ._util import as_scalar_if_0d

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d

#: Cited public presets for the log-law constants, flow → dict(kappa, B, citation) (pass kappa and B explicitly; there
#: is no default set; ``pandas.DataFrame(LOG_LAW_CONSTANTS)`` prints the table).
#: "classical": the textbook smooth-wall pair, e.g. https://en.wikipedia.org/wiki/Law_of_the_wall (κ ≈ 0.41, C+ ≈ 5.0;
#: read 2026-10-06) — also ``tools/benchmarks.py::WALL_LAW``.
#: "channel_dns_lee_moser_2015": κ = 0.384 ± 0.004 from the abstract of Lee & Moser, J. Fluid Mech. 774, 395 (2015),
#: arXiv:1410.7809 (read 2026-10-06); the abstract gives no additive constant, so ``B`` is None — obtain one from a fit
#: (:func:`fit_log_law`) or from the correlation (12.92) (:func:`nagib_chauhan_B`) and say which.
LOG_LAW_CONSTANTS: dict[str, dict] = {
    "classical": {"kappa": 0.41, "B": 5.0,
                  "citation": "common textbook smooth-wall pair (kappa ~ 0.41, C+ ~ 5.0), e.g. "
                              "https://en.wikipedia.org/wiki/Law_of_the_wall (read 2026-10-06)"},
    "channel_dns_lee_moser_2015": {"kappa": 0.384, "B": None,
                                   "citation": "Lee & Moser, J. Fluid Mech. 774, 395 (2015), arXiv:1410.7809, abstract: "
                                               "kappa = 0.384 +- 0.004 at Re_tau ~ 5200 (no additive constant quoted)"},
}

__all__ = [
    "LOG_LAW_CONSTANTS", "friction_velocity", "viscous_length", "friction_reynolds_number", "wall_units",
    "from_wall_units", "law_of_the_wall_groups", "defect_law_groups", "viscous_sublayer", "log_law", "log_law_defect",
    "log_law_crossing", "fit_log_law", "log_law_indicator", "layer_name", "friction_law_from_overlap", "spalding_yplus",
    "spalding_uplus", "spalding_slope", "coles_wake", "composite_profile_plus", "composite_profile", "velocity_defect",
    "total_stress", "channel_total_stress", "stress_partition", "channel_momentum_residual",
    "channel_pressure_gradient", "pipe_pressure_gradient", "wall_stress_from_pressure_gradient",
    "boundary_layer_stress_from_profile", "zpg_boundary_layer", "skin_friction_zpg", "nagib_chauhan_kappa",
    "nagib_chauhan_B", "rough_wall_log_law", "friction_velocity_from_wind", "drag_coefficient_neutral",
    "pipe_bulk_velocity_loglaw", "pipe_friction_factor_turbulent", "dimensionless_shear", "surface_layer_wind",
]


# ======================================================================================================================
# wall units and the two scaling laws
# ======================================================================================================================
def friction_velocity(tau0, rho):
    """Friction velocity ``u_* = sqrt(tau0/rho)``.

    Book: §12.9, Eq. (12.81) ``u_*² ≡ τ0/ρ``.  Parameters: tau0 wall shear stress [Pa] ≥ 0; rho density [kg/m³].
    Returns u_* [m/s].  Validation: V1 units and a hand value; the mutant τ0/ρ (no square root) fails the units check.
    Assumptions: tau0 >= 0 is the mean wall shear stress; constant density.
    """
    return _S(np.sqrt(_F(tau0) / _F(rho)))  # Eq. (12.81)


def viscous_length(nu, u_star):
    """Viscous wall unit ``l_ν = ν/u_*`` [m].  Book: §12.9, the definition after Eq. (12.80).  nu [m²/s]; u_star [m/s].

    Returns l_nu [m]; a float for float input.
    Assumptions: u_star > 0.
    Validation: V1 units and a hand value; y+ = y/l_nu round trip with wall_units.
    """
    return _S(_F(nu) / _F(u_star))


def friction_reynolds_number(u_star, delta, nu):
    """Friction Reynolds number ``Re_τ = δ⁺ = δ u_*/ν`` — the ratio of the outer to the viscous length.

    Book: §12.9 (δ⁺ = δ/l_ν, the number that differs between the profiles of Fig. 12.18).
    Parameters: u_star [m/s]; delta half-height, radius or layer thickness [m]; nu [m²/s].  Returns Re_τ [–].
    Assumptions: delta is the outer length of the flow in hand (half-height, pipe radius or layer thickness).
    Validation: V1 hand value (0.5 m/s, 0.03 m, 1.5e-5 m2/s gives 1000).
    """
    return _S(_F(u_star) * _F(delta) / _F(nu))


def wall_units(y, U, tau0, rho, nu):
    """Convert a mean profile to wall units: ``y⁺ = y u_*/ν``, ``U⁺ = U/u_*``.

    Book: §12.9, Eqs. (12.80)–(12.81) ``U⁺ ≡ U/u_* = f(y u_*/ν) = f(y⁺)``, ``u_*² = τ0/ρ``.
    Parameters: y [m] from the wall; U [m/s]; tau0 [Pa]; rho [kg/m³]; nu [m²/s].
    Returns (yplus, Uplus, u_star [m/s], l_nu [m]).  Scalar-callable.
    Validation: V1 round trip with :func:`from_wall_units` to 1e-14.
    Assumptions: fully turbulent flow over a smooth wall; constant density (the scaling itself is a change of variables).
    """
    us = np.sqrt(_F(tau0) / _F(rho))      # Eq. (12.81)
    lnu = _F(nu) / us
    return _S(_F(y) / lnu), _S(_F(U) / us), _S(us), _S(lnu)  # Eq. (12.80)


def from_wall_units(yplus, Uplus, u_star, nu):
    """Inverse of :func:`wall_units`: ``y = y⁺ ν/u_*`` [m], ``U = U⁺ u_*`` [m/s].  Book: §12.9, Eq. (12.80).

    Returns (y [m], U [m/s]); floats for float input.
    Assumptions: u_star > 0.
    Validation: V1 round trip with wall_units to 1e-14.
    """
    return _S(_F(yplus) * _F(nu) / _F(u_star)), _S(_F(Uplus) * _F(u_star))


def law_of_the_wall_groups() -> list:
    """The two dimensionless groups of the inner layer, by the Π theorem: ``U/u_*`` and ``y u_*/ν``.

    Book: §12.9, Eq. (12.79) ``U = U(ρ, τ0, ν, y)`` ⇒ Eq. (12.80).  Five variables, three dimensions, two groups.
    Returns the list of exponent dicts from ``core.dimensional.pi_groups`` with repeating variables (rho, tau0, nu):
    ``U ρ^{1/2} τ0^{−1/2}`` and ``y ρ^{−1/2} τ0^{1/2} ν^{−1}`` (exact ``Fraction`` exponents).
    Validation: V2 exponents exact; both groups dimensionless.
    Assumptions: close to a smooth wall the outer length does not matter: U depends on rho, tau0, nu and y only.
    """
    from .dimensional import pi_groups

    variables = {"U": "m/s", "y": "m", "rho": "kg/m**3", "tau0": "Pa", "nu": "m**2/s"}
    return pi_groups(variables, solution="U", repeating=("rho", "tau0", "nu"))  # Eq. (12.79) → (12.80)


def defect_law_groups() -> list:
    """The two dimensionless groups of the outer layer: ``U/u_*`` (used as the defect (U_∞ − U)/u_*) and ``y/δ``.

    Book: §12.9, Eq. (12.83) ``U = U(ρ, τ0, δ, y)`` ⇒ Eq. (12.84).  Repeating variables (rho, tau0, delta).
    Returns the list of exponent dicts (exact ``Fraction`` exponents).  Validation: V2 exact exponents.
    Assumptions: far from the wall viscosity does not matter directly: the defect depends on rho, tau0, delta and y only.
    """
    from .dimensional import pi_groups

    variables = {"U": "m/s", "y": "m", "rho": "kg/m**3", "tau0": "Pa", "delta": "m"}
    return pi_groups(variables, solution="U", repeating=("rho", "tau0", "delta"))  # Eq. (12.83) → (12.84)


def viscous_sublayer(yplus):
    """Mean velocity in the viscous sublayer: ``U⁺ = y⁺`` (valid to about y⁺ ≈ 5).

    Book: §12.9, Eq. (12.82) ``μ dU/dy = τ0 → U = τ0 y/μ``, i.e. U⁺ = y⁺.  Returns U⁺ [–].  Scalar-callable.
    Assumptions: y+ < ~5: the Reynolds stress is negligible and the stress equals tau0 (smooth wall).
    Validation: V1 slope dU+/dy+ = 1; equals the y+ -> 0 limit of spalding_uplus.
    """
    return _S(_F(yplus) * 1.0)  # Eq. (12.82)


def log_law(yplus, *, kappa: float, B: float):
    """Logarithmic law in inner variables: ``U⁺ = (1/κ) ln y⁺ + B``.

    Book: §12.9, Eq. (12.88) (overlap layer, large y⁺ and small y/δ).
    Parameters: yplus > 0 [–]; kappa von Kármán constant [–]; B additive constant [–] — **required keywords** (the
    values depend on the flow; pass the kappa and B of a cited preset such as ``LOG_LAW_CONSTANTS["classical"]``).
    Returns U⁺ [–].  Natural logarithm.  Scalar-callable.
    Validation: V1 ``y⁺ dU⁺/dy⁺ = 1/κ``; the log10-for-ln mutant fails the slope test.
    Assumptions: fully turbulent flow over a smooth wall; constant density; overlap layer (y+ >~ 30, y/delta <~ 0.15); kappa and B are empirical.
    """
    return _S(np.log(_F(yplus)) / kappa + B)  # Eq. (12.88)


def log_law_defect(xi, *, kappa: float, A: float):
    """Logarithmic law in outer (defect) variables: ``F(ξ) = (U_∞ − U)/u_* = −(1/κ) ln ξ + A``.

    Book: §12.9, Eq. (12.89), ξ = y/δ.  Parameters: xi in (0, 1]; kappa, A [–] required keywords.  Returns F [–].
    Assumptions: overlap layer, outer variables; A depends on the flow (pressure gradient, geometry).
    Validation: V1 -xi dF/dxi = 1/kappa; log_law + log_law_defect = friction_law_from_overlap.
    """
    return _S(-np.log(_F(xi)) / kappa + A)  # Eq. (12.89)


def friction_law_from_overlap(Re_tau, *, kappa: float, A: float, B: float):
    """Friction law obtained by adding the inner and outer logarithmic laws: ``U_∞/u_* = (1/κ) ln δ⁺ + A + B``.

    Book: §12.9, Eqs. (12.88) + (12.89) (the sum eliminates y; the step is ours).
    Parameters: Re_tau = δ⁺ [–]; kappa, A, B [–] required keywords.  Returns U_∞⁺ = sqrt(2/C_f) [–].
    Validation: V1 equals log_law(y⁺) + log_law_defect(y⁺/δ⁺) for any y⁺.
    Assumptions: an overlap layer exists (large Re_tau); kappa, A, B empirical.
    """
    return _S(np.log(_F(Re_tau)) / kappa + A + B)


def log_law_crossing(*, kappa: float, B: float) -> float:
    """y⁺ where the sublayer line ``U⁺ = y⁺`` (12.82) meets the logarithmic law (12.88).

    Book: §12.9, Fig. 12.18 (the two laws drawn on one semi-log plot; the crossing lies in the buffer layer).
    Returns y⁺ [–] (≈ 10.8 for κ = 0.41, B = 5.0).  brentq on y − ln(y)/κ − B in [1, 100]; residual asserted.
    Raises ValueError when the two curves do not cross exactly once in 1 < y⁺ < 100 — a single crossing needs
    1 < B < 100 − ln(100)/κ (e.g. the undamped mixing-length intercept B = (ln 4κ − 1)/κ ≈ −1.23 gives a log line
    that lies below the sublayer line everywhere: no crossing to report).
    Assumptions: the two laws are extrapolated to their intersection (neither holds there: it lies in the buffer layer).
    Validation: V1 residual < 1e-9; 10.80 for (0.41, 5.0); moves outward when B rises.
    """
    g = lambda y: y - np.log(y) / kappa - B  # noqa: E731
    if not (kappa > 0 and g(1.0) < 0.0 < g(100.0)):
        raise ValueError(f"no single crossing of U+ = y+ with the log law in 1 < y+ < 100 for kappa = {kappa}, B = {B}: "
                         "the sublayer line lies below the log line at y+ = 1 only if B > 1, and above it at y+ = 100 "
                         "only if B < 100 - ln(100)/kappa")
    root = brentq(g, 1.0, 100.0, xtol=1e-13, rtol=1e-13)
    assert abs(g(root)) < 1e-9
    return float(root)


def fit_log_law(yplus, Uplus, window: tuple = (30.0, 0.15), Re_tau: float | None = None) -> dict:
    """Least-squares fit of ``U⁺ = (1/κ) ln y⁺ + B`` to a profile inside a window.

    Book: §12.9, Eq. (12.88) and Fig. 12.18 (the straight part of the profile on semi-log axes).
    Parameters: yplus, Uplus (n,) profile in wall units; window = (lower y⁺, upper): an upper value ≤ 1 is a fraction of
    δ⁺ = ``Re_tau`` (default: the largest y⁺ given), a value > 1 is a y⁺.
    Returns dict(kappa, B, n_points, yplus_min, yplus_max, rms_residual).
    Validation: V1 recovers (κ, B) from a synthetic log profile to 1e-10; V5 public channel DNS at Re_τ ≈ 5200.
    Assumptions: a logarithmic region exists inside the window (Re_tau >~ 10^3); ordinary least squares in ln y+.
    """
    yp, Up = _F(yplus), _F(Uplus)
    ret = float(np.max(yp)) if Re_tau is None else float(Re_tau)
    hi = window[1] * ret if window[1] <= 1.0 else float(window[1])
    m = (yp >= window[0]) & (yp <= hi)
    if np.count_nonzero(m) < 2:
        raise ValueError("fewer than two points inside the fitting window")
    slope, intercept = np.polyfit(np.log(yp[m]), Up[m], 1)
    res = Up[m] - (slope * np.log(yp[m]) + intercept)
    return dict(kappa=float(1.0 / slope), B=float(intercept), n_points=int(np.count_nonzero(m)), yplus_min=float(window[0]),
                yplus_max=float(hi), rms_residual=float(np.sqrt(np.mean(res ** 2))))


def log_law_indicator(yplus, Uplus):
    """Indicator function ``y⁺ dU⁺/dy⁺ = dU⁺/d ln y⁺``: flat and equal to 1/κ where the logarithmic law holds.

    Book: §12.9, Eq. (12.87) ``y⁺ df/dy⁺ = 1/κ`` (the matching condition of the overlap layer).
    Parameters: yplus (n,) > 0 increasing; Uplus (n,).  Returns (n,) array [–].  Second-order differences in ln y⁺
    (``np.gradient`` with ``edge_order=2``).  Validation: V1 exact 1/κ on a log profile; equals y⁺ on U⁺ = y⁺.
    Assumptions: profile sampled smoothly enough for second-order differences in ln y+.
    """
    return np.gradient(_F(Uplus), np.log(_F(yplus)), edge_order=2)  # Eq. (12.87)


def layer_name(yplus, y_over_delta, sublayer: float = 5.0, buffer: float = 30.0, wake: float = 0.15):
    """Name of the layer a point lies in: "viscous sublayer", "buffer layer", "logarithmic layer" or "wake region".

    Book: §12.9, Fig. 12.18 and the text around it: sublayer y⁺ < 5, buffer 5 ≤ y⁺ < 30, wake beyond y/δ ≈ 0.15–0.2,
    the logarithmic layer in between.  The thresholds are nominal and are arguments.
    Parameters: yplus, y_over_delta scalars or arrays [–].  Returns a str (scalar input) or an array of str.
    Rule: the two near-wall names are decided by y⁺ first (viscosity matters there whatever y/δ is); otherwise the point
    is in the wake when y/δ > ``wake`` and in the logarithmic layer if not.  At low Re_τ (δ⁺ < buffer/wake = 200) no
    logarithmic layer is left.  Boundaries belong to the outer of the two layers (y⁺ = 5 → buffer, y⁺ = 30 → log).
    Validation: V7 both sides of each boundary and the boundary itself.
    Assumptions: nominal boundaries (y+ = 5, 30; y/delta = 0.15) - conventions, not sharp transitions.
    """
    yp, xi = np.broadcast_arrays(_F(yplus), _F(y_over_delta))
    out = np.where(yp < sublayer, "viscous sublayer",
                   np.where(yp < buffer, "buffer layer", np.where(xi > wake, "wake region", "logarithmic layer")))
    return str(out) if out.ndim == 0 else out


# ======================================================================================================================
# composite profiles (Spalding inner law, Coles wake)
# ======================================================================================================================
def spalding_yplus(Uplus, *, kappa: float, B: float):
    """Spalding's single formula for the whole inner layer, in its explicit direction y⁺(U⁺).

    Book: §12.9 (unnumbered "inner profile"):
    ``y⁺ = U⁺ + e^{−κB}[exp(κU⁺) − 1 − κU⁺ − (κU⁺)²/2 − (κU⁺)³/6]``.
    Parameters: Uplus ≥ 0; kappa, B required keywords.  Returns y⁺ [–].
    Limits: y⁺ → U⁺ for small U⁺ (the bracket starts at (κU⁺)⁴/24) and U⁺ → ln(y⁺)/κ + B for large y⁺.
    The bracket is evaluated as a series for κU⁺ < 0.1 to avoid cancellation.
    Validation: V1 both limits; the bracket's series to 4th order.
    Assumptions: fully turbulent flow over a smooth wall; constant density; zero or mild pressure gradient; inner layer only (no wake).
    """
    U = _F(Uplus)
    z = kappa * U
    zs = np.where(z < 0.1, z, 0.0)
    series = zs ** 4 / 24.0 * (1.0 + zs / 5.0 * (1.0 + zs / 6.0 * (1.0 + zs / 7.0 * (1.0 + zs / 8.0))))
    with np.errstate(over="ignore"):
        direct = np.exp(z) - 1.0 - z - z ** 2 / 2.0 - z ** 3 / 6.0
    return _S(U + np.exp(-kappa * B) * np.where(z < 0.1, series, direct))  # Spalding (1961), book §12.9


def spalding_slope(Uplus, *, kappa: float, B: float):
    """Slope ``dU⁺/dy⁺`` of Spalding's profile at a given U⁺ (the viscous share of the stress in a constant-stress layer).

    Book: §12.9 (derivative of the inner profile, ours): ``dy⁺/dU⁺ = 1 + κ e^{−κB}[exp(κU⁺) − 1 − κU⁺ − (κU⁺)²/2]``.
    Returns dU⁺/dy⁺ in (0, 1].
    Assumptions: as spalding_yplus.
    Validation: V1 equals the numerical derivative of spalding_uplus to 1e-6; 1 at the wall, 1/(kappa y+) far from it.
    """
    z = kappa * _F(Uplus)
    zs = np.where(z < 0.1, z, 0.0)
    series = zs ** 3 / 6.0 * (1.0 + zs / 4.0 * (1.0 + zs / 5.0 * (1.0 + zs / 6.0 * (1.0 + zs / 7.0))))
    with np.errstate(over="ignore"):
        direct = np.exp(z) - 1.0 - z - z ** 2 / 2.0
    return _S(1.0 / (1.0 + kappa * np.exp(-kappa * B) * np.where(z < 0.1, series, direct)))


def spalding_uplus(yplus, *, kappa: float, B: float):
    """Spalding's inner profile U⁺(y⁺): sublayer, buffer layer and logarithmic layer in one curve (implicit inverse).

    Book: §12.9 (unnumbered "inner profile", see :func:`spalding_yplus`).
    Parameters: yplus ≥ 0 scalar or array; kappa, B required keywords.  Returns U⁺ [–].  Scalar-callable.
    Method: ``brentq`` on ``spalding_yplus(U) − y⁺`` in the bracket [0, max(y⁺, ln(y⁺)/κ + B) + 5] (the profile lies below
    both the sublayer line and, by a bounded amount, above/below the log line); the residual is asserted < 1e-9 y⁺.
    Validation: V1 U⁺ → y⁺ (y⁺ → 0) and → log law (relative gap < 1e-2 at y⁺ = 10³); inverse residual.
    **Approximate** against channel DNS (Lee & Moser 2015, ``reference/ch12``; (κ, B) = (0.41, 5.0), inner layer
    y/δ < 0.15): largest |ΔU⁺| = 0.82, 0.79, 0.77 wall units at Re_τ ≈ 1000, 2000, 5200 and 0.85 at 550, all in the
    buffer layer; 1.19 at Re_τ ≈ 180, where there is hardly a logarithmic layer.  A fitted curve, not a law.
    Assumptions: as spalding_yplus.
    """
    yp = np.atleast_1d(_F(yplus)).astype(float)
    out = np.empty(yp.shape)
    for idx, y in np.ndenumerate(yp):
        if y < 0:
            raise ValueError("yplus must be >= 0")
        if y == 0.0:
            out[idx] = 0.0
            continue
        hi = max(y, np.log(max(y, 1.0)) / kappa + B) + 5.0
        g = lambda U, y=y: float(spalding_yplus(U, kappa=kappa, B=B)) - y  # noqa: E731
        root = brentq(g, 0.0, hi, xtol=1e-14, rtol=1e-14)
        assert abs(g(root)) <= 1e-9 * max(y, 1.0)
        out[idx] = root
    return _S(out.reshape(np.shape(yplus)))


def coles_wake(xi, kind: str = "cubic"):
    """Coles' wake function W(ξ), rising smoothly from W(0) = 0 to W(1) = 1.

    Book: §12.9 (unnumbered): the two simple choices ``W = 3ξ² − 2ξ³`` ("cubic") and ``W = sin²(πξ/2)`` ("sin2").
    Parameters: xi = y/δ, clipped to [0, 1].  Returns W [–].  Both have zero slope at ξ = 0 and ξ = 1.
    Assumptions: an empirical shape; any smooth function with W(0) = 0, W(1) = 1 and zero end slopes serves.
    Validation: V1 end values and end slopes; the two kinds differ by < 0.02.
    """
    x = np.clip(_F(xi), 0.0, 1.0)
    if kind == "cubic":
        return _S(3.0 * x ** 2 - 2.0 * x ** 3)
    if kind == "sin2":
        return _S(np.sin(0.5 * np.pi * x) ** 2)
    raise ValueError("kind must be 'cubic' or 'sin2'")


def composite_profile_plus(yplus, Re_tau: float, *, kappa: float, B: float, Pi: float, wake: str = "cubic"):
    """Composite mean profile in wall units: Spalding's inner law plus Coles' wake, ``U⁺ = U⁺_inner + (2Π/κ) W(y⁺/δ⁺)``.

    Book: §12.9 (the unnumbered inner and outer profiles: the outer one is the log law plus the wake term; using
    Spalding's curve in place of the bare logarithm joins the two, a standard construction, ours).
    Parameters: yplus ≥ 0; Re_tau = δ⁺; kappa, B, Pi (wake strength, ≥ 0; 0 = no wake) required keywords; wake kind.
    Returns U⁺ [–]; held at its edge value for y⁺ > δ⁺.  Scalar-callable.
    Validation: V1 reduces to :func:`spalding_uplus` for Pi = 0; edge value ln(δ⁺)/κ + B + 2Π/κ at large δ⁺;
    approximate (not a benchmark) against public channel DNS.
    Assumptions: fully turbulent flow over a smooth wall; constant density; the wake strength Pi is empirical and flow dependent (0 for none).
    """
    yp = np.minimum(_F(yplus), float(Re_tau))
    U = _F(spalding_uplus(yp, kappa=kappa, B=B)) + 2.0 * Pi / kappa * _F(coles_wake(yp / float(Re_tau), wake))
    return _S(U)


def composite_profile(y, delta: float, u_star: float, nu: float, *, kappa: float, B: float, Pi: float,
                      wake: str = "cubic"):
    """Dimensional composite mean profile U(y) [m/s] of a wall layer of thickness δ (see :func:`composite_profile_plus`).

    Book: §12.9 (inner and outer profile formulas).  Parameters: y [m] from the wall; delta [m]; u_star [m/s]; nu [m²/s];
    kappa, B, Pi required keywords.  Returns U [m/s].
    Assumptions: as composite_profile_plus.
    Validation: V1 identical to composite_profile_plus after wall scaling.
    """
    lnu = nu / u_star
    return _S(u_star * _F(composite_profile_plus(_F(y) / lnu, delta / lnu, kappa=kappa, B=B, Pi=Pi, wake=wake)))


def velocity_defect(y, U, U_inf, u_star, delta):
    """Outer scaling of a profile: ``ξ = y/δ`` and the defect ``(U_∞ − U)/u_*``.

    Book: §12.9, Eq. (12.84) ``(U_∞ − U)/u_* = F(y/δ)``.  Parameters: y [m]; U, U_inf, u_star [m/s]; delta [m].
    Returns (xi, defect) [–].
    Assumptions: delta and u_star are known for the profile (the scaling is a change of variables).
    Validation: V1 round trip; a log-law profile gives log_law_defect.
    """
    return _S(_F(y) / _F(delta)), _S((_F(U_inf) - _F(U)) / _F(u_star))  # Eq. (12.84)


# ======================================================================================================================
# total stress, pressure gradient
# ======================================================================================================================
def total_stress(y, U, uv, mu: float, rho: float):
    """Total mean shear stress ``τ̄ = μ dU/dy − ρ mean(uv)`` of a unidirectional mean flow.

    Book: §12.9, the definition under Eq. (12.76).  Parameters: y (n,) [m]; U (n,) [m/s]; uv (n,) Reynolds shear
    correlation mean(uv) [m²/s²]; mu [Pa s]; rho [kg/m³].  Returns dict(total, viscous, reynolds) [Pa].
    dU/dy by second-order differences (``np.gradient``, ``edge_order=2``).
    Assumptions: unidirectional mean flow U(y); constant mu and rho.
    Validation: V1 laminar limit (uv = 0) gives mu dU/dy; parts sum to the total by construction.
    """
    visc = mu * np.gradient(_F(U), _F(y), edge_order=2)
    rey = -rho * _F(uv)
    return dict(total=visc + rey, viscous=visc, reynolds=rey)


def channel_total_stress(y, h, tau0):
    """Total stress across a fully developed channel: the straight line ``τ̄(y) = τ0 (1 − 2y/h)``.

    Book: §12.9, Eqs. (12.76)–(12.77) (∂τ̄/∂y = dP/dx is the same at every y) with (12.90) dP/dx = −2τ0/h; the linear
    form is written out by us.  **h is the full channel height** (walls at y = 0 and y = h): τ̄ = τ0 at the lower wall,
    0 on the centreline y = h/2, −τ0 at the upper wall.
    Parameters: y [m]; h [m]; tau0 wall shear stress [Pa].  Returns τ̄ [Pa].
    Validation: V1 satisfies (12.76) with (12.90) exactly; the half-height mutant τ0(1 − y/h) fails at y = h/2.
    Assumptions: fully developed channel (statistics independent of x), smooth or rough walls - the result is a force balance.
    """
    return _S(_F(tau0) * (1.0 - 2.0 * _F(y) / _F(h)))


def channel_momentum_residual(y, tau, dPdx):
    """Residual of the channel mean-momentum balance ``0 = −dP/dx + dτ̄/dy``.

    Book: §12.9, Eq. (12.76) (first equation); by (12.77) dP/dx does not depend on y.
    Parameters: y (n,) [m]; tau (n,) total stress [Pa]; dPdx [Pa/m] (scalar).  Returns (n,) residual [Pa/m].
    Assumptions: fully developed channel flow.
    Validation: V1 zero to round-off for channel_total_stress with channel_pressure_gradient.
    """
    return np.gradient(_F(tau), _F(y), edge_order=2) - float(dPdx)  # Eq. (12.76)


def channel_pressure_gradient(tau0, h):
    """Pressure gradient of fully developed channel flow: ``dP/dx = −2 τ0/h`` (h = full height).

    Book: §12.9, Eq. (12.90).  Returns dP/dx [Pa/m] (negative for τ0 > 0).
    Assumptions: fully developed flow between parallel walls a distance h apart; equal stress on both walls.
    Validation: V1 force balance on a slab: dP/dx * h = -2 tau0.
    """
    return _S(-2.0 * _F(tau0) / _F(h))  # Eq. (12.90)


def pipe_pressure_gradient(tau0, d):
    """Pressure gradient of fully developed pipe flow: ``dP/dx = −4 τ0/d`` (d = diameter).

    Book: §12.9, Eq. (12.91) (the turbulent mean obeys the same force balance as laminar flow, ch08 Eq. (8.8)).
    Returns dP/dx [Pa/m] (negative for tau0 > 0).
    Assumptions: fully developed flow in a round pipe.
    Validation: V1 force balance on a plug: dP/dx * pi d^2/4 = -tau0 * pi d; agrees with core.laminar.pipe_wall_stress.
    """
    return _S(-4.0 * _F(tau0) / _F(d))  # Eq. (12.91)


def wall_stress_from_pressure_gradient(dPdx, d):
    """Wall shear stress of a pipe from its pressure gradient: ``τ0 = −(d/4) dP/dx``.  Book: §12.9, Eq. (12.91).

    Returns tau0 [Pa].
    Assumptions: fully developed flow in a round pipe.
    Validation: V1 inverse of pipe_pressure_gradient to round-off.
    """
    return _S(-_F(dPdx) * _F(d) / 4.0)  # Eq. (12.91)


def stress_partition(yplus, Re_tau: float, *, kappa: float, B: float) -> dict:
    """Viscous and Reynolds parts of the linear total stress of a channel, in wall units (a model partition).

    Book: §12.9, Fig. 12.17a and Eqs. (12.76)–(12.77): the total τ̄⁺ = 1 − y⁺/Re_τ (Re_τ = half-height in wall units) is
    exact; the viscous part is dU⁺/dy⁺ and the Reynolds part −mean(uv)⁺ is the rest.
    Model: dU⁺/dy⁺ is taken from Spalding's inner profile (:func:`spalding_slope`), capped by the total stress — so the
    Reynolds part is ≥ 0 and the result is trustworthy in the inner layer (y/δ ≲ 0.2), indicative beyond.
    Parameters: yplus in [0, Re_tau]; kappa, B required keywords.  Returns dict(total, viscous, reynolds) [–].
    Validation: V1 the parts sum to the linear law; viscous = total = 1 at the wall.  **Approximate (model)**: the two
    shares are equal at y⁺ ≈ 9.9 for (κ, B) = (0.41, 5.0) at Re_τ = 1000 (10.3 at Re_τ = 180, 9.8 at 5200), in the
    buffer layer; the viscous share is within 0.08 τ0 of channel DNS for y/δ < 0.2 (tested at Re_τ ≈ 5200).
    Assumptions: fully developed channel; the viscous share is modelled with Spalding's profile (trustworthy for y/delta <~ 0.2).
    """
    yp = _F(yplus)
    total = 1.0 - yp / float(Re_tau)
    Up = spalding_uplus(yp, kappa=kappa, B=B)
    visc = np.minimum(_F(spalding_slope(Up, kappa=kappa, B=B)), np.maximum(total, 0.0))
    return dict(total=_S(total), viscous=_S(visc), reynolds=_S(total - visc))


def boundary_layer_stress_from_profile(x, y, U, V, dPdx, rho: float):
    """Total stress across a boundary layer by integrating the mean momentum equation from the wall.

    Book: §12.9, Eq. (12.78) ``U ∂U/∂x + V ∂U/∂y = −(1/ρ) ∂P/∂x + (1/ρ) ∂τ̄/∂y`` (Fig. 12.17b).
    Parameters: x (nx,) [m]; y (ny,) [m] starting at the wall; U, V (ny, nx) mean velocities [m/s] (``U[j, i]`` =
    U(y_j, x_i)); dPdx scalar or (nx,) [Pa/m]; rho [kg/m³].
    Returns ``tau_minus_tau0`` (ny, nx) [Pa] = τ̄(x, y) − τ0(x) = ∫_0^y [ρ(U U_x + V U_y) + dP/dx] dy (trapezoid; second-
    order differences).  Add the wall stress to get τ̄; at the edge τ̄ → 0, which returns τ0 (the momentum integral).
    Assumptions: steady two-dimensional boundary layer; mean fields on a tensor grid; normal-stress and streamwise viscous terms neglected (thin layer).
    Validation: V1 Blasius profile (ch09): the result equals mu dU/dy - tau0 to the differencing error; V3 order 2.
    """
    x, y, U, V = _F(x), _F(y), _F(U), _F(V)
    rhs = rho * (U * np.gradient(U, x, axis=1, edge_order=2) + V * np.gradient(U, y, axis=0, edge_order=2)) + _F(dPdx)
    out = np.zeros_like(U)
    out[1:] = np.cumsum(0.5 * (rhs[1:] + rhs[:-1]) * np.diff(y)[:, None], axis=0)  # Eq. (12.78) integrated in y
    return out


# ======================================================================================================================
# zero-pressure-gradient boundary layer and skin friction
# ======================================================================================================================
def zpg_boundary_layer(x, U_inf, nu, *, kappa: float) -> dict:
    """Thicknesses and skin friction of a smooth flat-plate turbulent boundary layer in zero pressure gradient.

    Book: §12.9 (unnumbered fits attributed to Monkewitz et al. 2007, for Re_x > 10⁶):
    ``θ ≈ 0.016 x Re_x^{−0.15}``; ``δ* ≈ θ exp{7.11 κ/ln Re_θ}``; ``δ_99 = 0.2 δ*[κ^{−1} ln Re_δ* + 3.30]``;
    ``C_f = τ0/(½ρU_∞²) ≅ 2.0/[κ^{−1} ln Re_δ* + 3.30]²``.
    Parameters: x distance from the (virtual) origin [m]; U_inf [m/s]; nu [m²/s]; kappa von Kármán constant (required).
    Returns dict(Re_x, theta, delta_star, delta99 [m], H = δ*/θ, Cf, U_inf_plus = sqrt(2/C_f) (the bracket),
    u_star [m/s], Re_theta, Re_delta_star, Re_tau = δ_99 u_*/ν, valid = Re_x > 1e6).
    Assumptions: smooth wall, zero pressure gradient, fully turbulent from the origin (the book warns that the virtual
    origin is not accounted for).
    Validation: V2 C_f = 2/(U_∞⁺)² identity; V7 θ < δ* < δ_99, 1.2 < H < 1.5 for 10⁶ ≤ Re_x ≤ 10⁹; V6 book example
    (private).
    """
    x, U, n = _F(x), _F(U_inf), _F(nu)
    Re_x = U * x / n
    theta = 0.016 * x * Re_x ** (-0.15)                    # momentum thickness fit
    Re_th = U * theta / n
    dstar = theta * np.exp(7.11 * kappa / np.log(Re_th))   # displacement thickness fit
    Re_ds = U * dstar / n
    bracket = np.log(Re_ds) / kappa + 3.30                 # = U_inf/u_*
    d99 = 0.2 * dstar * bracket                            # 99 % thickness fit
    Cf = 2.0 / bracket ** 2                                # skin-friction fit
    us = U / bracket
    return dict(Re_x=_S(Re_x), theta=_S(theta), delta_star=_S(dstar), delta99=_S(d99), H=_S(dstar / theta), Cf=_S(Cf),
                U_inf_plus=_S(bracket), u_star=_S(us), Re_theta=_S(Re_th), Re_delta_star=_S(Re_ds),
                Re_tau=_S(d99 * us / n), valid=bool(np.all(Re_x > 1.0e6)))


def skin_friction_zpg(Re_x, law: str = "monkewitz", kappa: float | None = None):
    """Local skin-friction coefficient ``C_f = τ0/(½ρU_∞²)`` of a smooth flat plate, turbulent, by a named correlation.

    Book: §12.9 (unnumbered): "monkewitz" — the fit of :func:`zpg_boundary_layer` (needs ``kappa``);
    "schultz_grunow" — ``C_f ≅ 0.370 (log10 Re_x)^{−2.584}``; "white" — ``C_f ≅ 0.455/[ln(0.06 Re_x)]²``.
    Also "power_fifth" — ``0.8·0.074 Re_x^{−1/5}``, the local form of the one-seventh-power drag law used in ch09
    (``core.boundary_layer.plate_drag_coefficient``; C_f = d(x C_D)/dx; coefficient unverified, see there).
    Parameters: Re_x = U_∞ x/ν [–] (> 10⁶ for the three book laws).  Returns C_f [–].  Scalar-callable.
    Validation: V7 the three book laws exceed the laminar 0.664/sqrt(Re_x) and agree with each other to the following
    spread (max/min − 1), which depends on the κ passed to "monkewitz": with κ = 0.384 within 10 % for
    10⁶ ≤ Re_x ≤ 10⁸ (8.7 % at 10⁶, 7.7 % at 10⁷, 8.8 % at 10⁸) and 12.0 % at 10⁹; with κ = 0.41, 6.7 % at 10⁶, 10.1 % at 10⁷,
    13.7 % at 10⁸ and 17.3 % at 10⁹.  They are fits to different data sets, not one law.
    Assumptions: fully turbulent flow over a smooth wall; constant density; zero pressure gradient; Re_x > 1e6; empirical correlations.
    """
    Re = _F(Re_x)
    if law == "monkewitz":
        if kappa is None:
            raise ValueError("law='monkewitz' needs kappa (the von Kármán constant has no silent default)")
        return zpg_boundary_layer(Re, 1.0, 1.0, kappa=kappa)["Cf"]
    if law == "schultz_grunow":
        return _S(0.370 * np.log10(Re) ** (-2.584))
    if law == "white":
        return _S(0.455 / np.log(0.06 * Re) ** 2)
    if law == "power_fifth":
        return _S(0.8 * 0.074 * Re ** (-0.2))
    raise ValueError("law must be 'monkewitz', 'schultz_grunow', 'white' or 'power_fifth'")


def nagib_chauhan_kappa(B):
    """von Kármán constant paired with an additive constant B by the empirical correlation ``κB = 1.6[exp(0.1663 B) − 1]``.

    Book: §12.9, Eq. (12.92) (Nagib & Chauhan 2008; stated for 0.15 < κ < 0.80 and −4 < B < 12).
    Parameters: B [–], B ≠ 0 (the limit B → 0 is κ = 1.6·0.1663).  Returns κ [–]; NaN outside −4 < B < 12.
    Assumptions: an empirical correlation across flows and pressure gradients; valid for -4 < B < 12.
    Validation: V1 the pair (kappa, B) satisfies kappa B = 1.6[exp(0.1663 B) - 1] to round-off; B -> 0 limit.
    """
    b = _F(B)
    with np.errstate(divide="ignore", invalid="ignore"):
        k = np.where(np.abs(b) < 1e-12, 1.6 * 0.1663, 1.6 * np.expm1(0.1663 * b) / np.where(b == 0, 1.0, b))  # Eq. (12.92)
    return _S(np.where((b > -4.0) & (b < 12.0), k, np.nan))


def nagib_chauhan_B(kappa: float) -> float:
    """Additive constant B paired with a von Kármán constant κ by (12.92) (inverse of :func:`nagib_chauhan_kappa`).

    Book: §12.9, Eq. (12.92).  κ(B) is monotone increasing on −4 < B < 12; brentq on that branch, residual asserted.
    Raises ValueError if κ is outside the range the correlation reaches there (≈ 0.19 … 0.85).
    Returns B [-] (float).
    Assumptions: as nagib_chauhan_kappa; kappa on the monotone branch.
    Validation: V1 round trip nagib_chauhan_kappa(nagib_chauhan_B(kappa)) = kappa to 1e-10.
    """
    g = lambda b: float(nagib_chauhan_kappa(b)) - kappa  # noqa: E731
    lo, hi = -4.0 + 1e-9, 12.0 - 1e-9
    if g(lo) * g(hi) > 0:
        raise ValueError("kappa is outside the range of the correlation (12.92) for -4 < B < 12")
    root = brentq(g, lo, hi, xtol=1e-13, rtol=1e-13)
    assert abs(g(root)) < 1e-10
    return float(root)


# ======================================================================================================================
# rough walls
# ======================================================================================================================
def rough_wall_log_law(y, u_star, y0, *, kappa: float):
    """Logarithmic profile over a hydrodynamically rough wall: ``U = (u_*/κ) ln(y/y0)``.

    Book: §12.9, Eq. (12.93); y0 is the roughness length — the height at which the logarithm extrapolates to U = 0.
    Viscosity does not appear, so the friction coefficient no longer depends on the Reynolds number.
    Parameters: y [m]; u_star [m/s]; y0 [m] > 0; kappa required keyword.  Returns U [m/s]; **NaN for y < y0** (inside
    the roughness elements the law has no meaning); U = 0 at y = y0.  Scalar-callable.
    Validation: V1 identical to :func:`log_law` with B = −(1/κ) ln y0⁺.
    Assumptions: fully rough wall (roughness elements much higher than the viscous sublayer); neutral stratification.
    """
    yy = _F(y)
    with np.errstate(divide="ignore", invalid="ignore"):
        U = _F(u_star) / kappa * np.log(yy / _F(y0))  # Eq. (12.93)
    return _S(np.where(yy >= _F(y0), U, np.nan))


def friction_velocity_from_wind(U_ref, z_ref, z0, *, kappa: float):
    """Friction velocity from one wind measurement over a rough surface (neutral): ``u_* = κ U_ref/ln(z_ref/z0)``.

    Book: §12.9, Eq. (12.93) solved for u_*.  U_ref [m/s] at height z_ref [m]; z0 roughness length [m].  Returns [m/s].
    Assumptions: neutral surface layer over a fully rough surface; z_ref inside the logarithmic layer.
    Validation: V1 inverse of rough_wall_log_law to round-off.
    """
    return _S(kappa * _F(U_ref) / np.log(_F(z_ref) / _F(z0)))  # Eq. (12.93)


def drag_coefficient_neutral(z_ref, z0, *, kappa: float):
    """Neutral bulk drag coefficient ``C_D = u_*²/U(z_ref)² = [κ/ln(z_ref/z0)]²`` (ours, from (12.93)).

    Book: §12.9, Eq. (12.93).  Returns C_D [–] (surface stress τ0 = ρ C_D U_ref²); falls as z0 falls.
    Assumptions: neutral surface layer over a fully rough surface.
    Validation: V1 equals (u_star/U)^2 from rough_wall_log_law; falls as z0 falls.
    """
    return _S((kappa / np.log(_F(z_ref) / _F(z0))) ** 2)


# ======================================================================================================================
# pipe friction
# ======================================================================================================================
def pipe_bulk_velocity_loglaw(a, u_star, nu, *, kappa: float, B: float):
    """Cross-section averaged speed of a smooth pipe when the log law is assumed to hold across the whole pipe.

    Book: §12.9, Exercise 12.34 (form only): ``U_av ≅ u_*[(1/κ) ln(a u_*/ν) + B − 3/(2κ)]`` from
    ``U_av = (2/a²) ∫_0^a U(y)(a − y) dy`` with (12.88).
    Parameters: a pipe radius [m]; u_star [m/s]; nu [m²/s]; kappa, B required keywords.  Returns U_av [m/s].
    Validation: V1 equals the numerical integral of the log law over the cross-section to 1e-10.
    Assumptions: the log law is assumed to hold over the whole cross-section (it neglects the sublayer and the wake: good to a few per cent).
    """
    return _S(_F(u_star) * (np.log(_F(a) * _F(u_star) / _F(nu)) / kappa + B - 1.5 / kappa))


def pipe_friction_factor_turbulent(Re_d, kappa: float | None = None, B: float | None = None):
    """Darcy friction factor f of fully developed turbulent flow in a smooth pipe (implicit friction law).

    Book: §12.9, Exercise 12.34 (form only).  With f = 8(u_*/U_av)² and Re_d = U_av d/ν the bulk-velocity formula of
    :func:`pipe_bulk_velocity_loglaw` becomes
    ``1/sqrt(f) = (1/sqrt(8)) [(1/κ) ln(Re_d sqrt(f)/(2 sqrt(8))) + B − 3/(2κ)]``.
    Parameters: Re_d [–] (turbulent: ≳ 4000); kappa, B: both given → the law above; **both None → Prandtl's law**
    ``1/sqrt(f) = 2.0 log10(Re_d sqrt(f)) − 0.8`` (Prandtl 1935, as quoted by McKeon, Zagarola & Smits, J. Fluid Mech.
    538, 429 (2005) — constants adjusted to data, which is why they differ slightly from the derived ones).
    Returns f [–].  brentq on x = 1/sqrt(f) in [1, 100] (f between 1 and 10⁻⁴); residual asserted < 1e-12.
    Validation: V1 residual; f decreases with Re_d and exceeds the laminar 64/Re_d for Re_d ≥ 10⁴; V5 derived
    constants from (κ, B) = (0.41, 5.0) are 1.99 and −1.02 against Prandtl's 2.0 and −0.8.
    Assumptions: smooth pipe, fully developed turbulent flow (Re_d >~ 4000).
    """
    if (kappa is None) != (B is None):
        raise ValueError("give both kappa and B, or neither (Prandtl's law)")
    Re = np.atleast_1d(_F(Re_d)).astype(float)
    out = np.empty(Re.shape)
    for idx, R in np.ndenumerate(Re):
        if kappa is None:
            g = lambda x, R=R: x - (2.0 * np.log10(R / x) - 0.8)  # noqa: E731
        else:
            g = lambda x, R=R: x - (np.log(R / (x * 2.0 * np.sqrt(8.0))) / kappa + B - 1.5 / kappa) / np.sqrt(8.0)  # noqa: E731
        x = brentq(g, 1.0, 100.0, xtol=1e-14, rtol=1e-14)
        assert abs(g(x)) < 1e-12
        out[idx] = 1.0 / x ** 2
    return _S(out.reshape(np.shape(Re_d)))


# ======================================================================================================================
# stratified surface layer
# ======================================================================================================================
def dimensionless_shear(zeta, beta: float = 5.0, unstable: str = "log_linear"):
    """Dimensionless wind shear ``φ_m = (κ z/u_*) dU/dz`` as a function of ζ = z/L_M.

    Book: §12.11, the log-linear profile ``U = (u_*/κ)[ln(z/z0) + 5 z/L_M]`` ⇒ ``φ_m = 1 + β ζ`` with the book's β = 5.
    ``unstable="businger_dyer"`` replaces the unstable side (ζ < 0) by ``φ_m = (1 − 16 ζ)^{−1/4}``, the
    "Businger–Dyer" form.  **The coefficients are literature values (unstable 16 as coded; the stable side stays the
    book's 1 + β ζ with the default β = 5 — a caller wanting the 4.7 found in part of that literature passes
    ``beta=4.7``), attribution not verified first-hand — see reference/ch12/SOURCES.md.**  Nothing here is a benchmark:
    the tests only prove that the wind profile is the integral of the coded φ_m.
    Parameters: zeta [–]; beta [–].  Returns φ_m [–]; φ_m(0) = 1; **NaN where the log-linear form would give φ_m ≤ 0**
    (ζ ≤ −1/β: the linear correction is a small-|ζ| formula and has no meaning there).
    Assumptions: constant-flux surface layer; empirical stability functions, valid for moderate |z/L_M|.
    Validation: V1 phi_m(0) = 1 on both branches; equals (kappa z/u_*) dU/dz of surface_layer_wind by differences.
    """
    z = _F(zeta)
    phi = 1.0 + beta * z
    if unstable == "businger_dyer":
        phi = np.where(z < 0, (1.0 - 16.0 * np.minimum(z, 0.0)) ** -0.25, phi)
    elif unstable != "log_linear":
        raise ValueError("unstable must be 'log_linear' or 'businger_dyer'")
    return _S(np.where(phi > 0, phi, np.nan))


def _psi_businger_dyer(zeta):
    x = (1.0 - 16.0 * np.minimum(_F(zeta), 0.0)) ** 0.25
    return 2.0 * np.log(0.5 * (1.0 + x)) + np.log(0.5 * (1.0 + x * x)) - 2.0 * np.arctan(x) + 0.5 * np.pi


def surface_layer_wind(z, u_star, z0, L_M=np.inf, *, kappa: float, beta: float = 5.0, unstable: str = "log_linear",
                       psi_at_z0: bool = False):
    """Mean wind in the stratified surface layer: the log-linear profile ``U = (u_*/κ)[ln(z/z0) + β z/L_M]``.

    Book: §12.11 (unnumbered, after (12.111)), with the book's coefficient β = 5; ``L_M = inf`` gives the neutral
    rough-wall law (12.93) ``U = (u_*/κ) ln(z/z0)``.  Stable (L_M > 0): more shear than neutral; unstable (L_M < 0): less.
    Parameters
    ----------
    z : height [m] ≥ z0.   u_star : friction velocity [m/s].   z0 : roughness length [m].
    L_M : Monin–Obukhov length [m] (sign of (12.110): positive stable, negative unstable, ±inf neutral).
    kappa : von Kármán constant [–] (required keyword).   beta : coefficient of the linear correction [–].
    unstable : "log_linear" (the book's formula on both sides of neutral) or "businger_dyer": for L_M < 0 the integral
        of φ_m = (1 − 16ζ)^{−1/4}, ``U = (u_*/κ)[ln(z/z0) − ψ_m(z/L_M)]`` with
        ``ψ_m = 2 ln((1 + x)/2) + ln((1 + x²)/2) − 2 arctan x + π/2``, ``x = (1 − 16ζ)^{1/4}`` (not from the book; the
        coefficient 16 is a literature value whose attribution is not verified first-hand — see
        reference/ch12/SOURCES.md; ψ_m itself is proved by quadrature to be ∫_0^ζ (1 − φ_m)/ζ′ dζ′ of the coded φ_m);
        stable side unchanged.  With "log_linear" the result is **NaN where ζ = z/L_M ≤ −1/β** (φ_m ≤ 0: the linear
        correction has no meaning there), exactly where :func:`dimensionless_shear` returns NaN.
    psi_at_z0 : add Paulson's lower-limit term ``+ψ_m(z0/L_M)`` (so that U(z0) = 0 exactly).  Default False — the
        common form, which neglects it because z0 ≪ |L_M| (ψ_m(z0/L_M) ≈ −4 z0/L_M, of order z0/|L_M|).
    Returns
    -------
    U [m/s]; NaN below z0, and — with ``unstable="log_linear"`` — NaN where z/L_M ≤ −1/β.  A float for scalar input.
    Assumptions: horizontally uniform, steady surface layer with constant fluxes; z ≪ boundary-layer depth; the linear
    correction holds for moderate |z/L_M| only (for z/L_M < −1/β the log-linear formula would give a wind that
    decreases with height and can turn negative: outside its range, hence NaN — use "businger_dyer" there).
    Validation: V1 neutral limit = (12.93); V7 sign test with both signs of L_M; ψ_m is the numerical integral of
    (1 − φ_m)/ζ.
    """
    zz, L = _F(z), _F(L_M)
    with np.errstate(divide="ignore", invalid="ignore"):
        zeta = np.where(np.isinf(L), 0.0, zz / np.where(np.isinf(L), 1.0, L))
        log = np.log(zz / _F(z0))
        if unstable == "businger_dyer":
            psi0 = 0.0
            if psi_at_z0:
                psi0 = _psi_businger_dyer(np.where(np.isinf(L), 0.0, _F(z0) / np.where(np.isinf(L), 1.0, L)))
            corr = np.where(zeta < 0, -(_psi_businger_dyer(zeta) - psi0), beta * zeta)
        elif unstable == "log_linear":
            # outside the range of the linear correction (φ_m = 1 + βζ ≤ 0, i.e. ζ ≤ −1/β) the formula gives a wind
            # that decreases with height and can turn negative: NaN there, as in dimensionless_shear
            corr = np.where(1.0 + beta * zeta > 0, beta * zeta, np.nan)
        else:
            raise ValueError("unstable must be 'log_linear' or 'businger_dyer'")
        U = _F(u_star) / kappa * (log + corr)  # log-linear profile, §12.11
    return _S(np.where(zz >= _F(z0), U, np.nan))

