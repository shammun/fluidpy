"""Finite differences as objects: stencils and their order, 1-D convection–diffusion schemes, von Neumann stability, CFL,
cell Péclet number, truncation errors (modified equations), the steady convection–diffusion layer and its wiggles.

Book: Kundu, Cohen & Dowling 5e, Ch. 10 §10.2 Eqs. (10.1)–(10.31), §10.4 Eqs. (10.84)–(10.94), the one-sided stencils of
§10.5 (printed p. 450 and (10.148)–(10.150)); Exercise 10.3's heated-rod series (10.199) is used (with our own numbers) in
:func:`lax_demo`.  Every equation was transcribed from the rendered page images (chapters/pages/ch10/p450–p456,
p464–p466, p477, p481, p497).

Grids.  Everything in this module is **node-based**: T_i^n ≈ T(x_i, t_n) with x_i = iΔx (Fig. 10.1).  Periodic grids do
**not** repeat the end node (x = x0 + Δx·arange(N)).  2-D arrays use the project layout ``f[j, i]`` = f(y_j, x_i).

Notation (analysis §9).  FTCS numbers α = uΔt/(2Δx), β = DΔt/Δx² (10.11); Courant number C = uΔt/Δx = 2α.  Fourier modes:
the book writes e^{iπkx_i} so that θ = kπΔx (10.26); every function here takes θ [rad per cell] directly.
Printed slip R11 (10.29)–(10.30): the upwind side and the CFL number depend on the sign of u — coded with |u| and the side
from sign(u); ``scheme="upwind_printed"`` keeps the book's stencil T_i − T_{i−1} for any sign (a test must see it fail for
u < 0).  Slip R3 (below (10.93), "forward difference"): ``steady_cd_fd(scheme="forward")`` is the true forward (downwind)
difference: its root r = 1/(1 − R_cell) is positive (monotone, anti-diffusive) for R_cell < 1 and negative — wiggles — only
for R_cell > 1.
Units: x, L, Δx [m]; t, Δt [s]; u [m/s]; D [m²/s]; T any unit (a temperature, a concentration).
Also used by Ch. 11 (discretised stability operators), Ch. 12 (numerical diffusion vs eddy viscosity), Ch. 13 (gravity-wave
CFL, upwind tracer advection), Ch. 15 (Lax–Wendroff/MacCormack).
"""
from __future__ import annotations

from functools import lru_cache
from typing import Callable, Sequence

import numpy as np

from ._util import as_scalar_if_0d, require_nonnegative, require_positive

__all__ = [
    # stencils
    "fd_weights", "fd_weights_float", "taylor_table_sympy", "fd_leading_error", "stencil_taylor_coefficients",
    "test_function", "stencil_derivative", "stencil_error", "fd_derivative", "fd_apply", "mixed_derivative_onesided",
    # model problem and schemes
    "transport_rhs", "advected_gaussian", "ftcs_coefficients", "cfl_number", "diffusion_number", "cell_peclet",
    "transport_1d_step", "solve_transport_1d", "error_norm", "observed_order", "convergence_study", "convergence_rates",
    "richardson_extrapolate", "grid_convergence_index", "richardson_three",
    "truncation_error_sympy", "truncation_terms", "ode_scheme_order", "ode_scheme_errors",
    # stability
    "propagate_error", "fourier_mode", "amplification_factor", "amplification_modulus", "amplification_curve",
    "max_amplification", "worst_theta", "is_von_neumann_stable", "ftcs_amplification_modulus2", "ftcs_stable",
    "stability_verdict", "phase_error", "numerical_diffusivity", "ftcs2d_max_amplification", "advect_periodic",
    "rod_heating_exact", "lax_demo",
    # steady convection–diffusion
    "steady_cd_exact", "cd_layer_thickness", "stretched_grid", "steady_cd_fd", "steady_cd_discrete_exact",
    "discrete_root", "wiggle_indicator", "SCHEMES_1D",
]

SCHEMES_1D = ("ftcs", "btcs", "cn", "upwind", "upwind_printed", "lax_wendroff")


# ======================================================================================================================
# Stencils from Taylor series — (10.4)–(10.8), p. 450, (10.148)–(10.150)
# ======================================================================================================================
@lru_cache(maxsize=256)
def _fd_weights_cached(offsets: tuple, m: int) -> tuple:
    import sympy as sp

    n = len(offsets)
    if m >= n:
        raise ValueError(f"need more than m = {m} points for an m-th derivative, got offsets {offsets}")
    # Taylor matching: sum_k w_k s_k^p / p! = delta_{p m} for p = 0..n-1  (Vandermonde system, exact rationals)
    A = sp.Matrix([[sp.Rational(s) ** p / sp.factorial(p) for s in offsets] for p in range(n)])
    b = sp.Matrix([1 if p == m else 0 for p in range(n)])
    w = A.LUsolve(b)
    return tuple(sp.nsimplify(x) for x in w)


def fd_weights(offsets: Sequence[int], m: int = 1, symbolic: bool = True):
    """Weights of the finite-difference stencil for the m-th derivative on the points x₀ + s·h, s ∈ ``offsets``.

    Book: §10.2 — the stencils (10.6)–(10.8) are obtained by combining the Taylor expansions (10.4)–(10.5) so that the
    low-order terms cancel; the one-sided second-order stencils of printed p. 450 and (10.148)–(10.149) are the same
    construction on (0, 1, 2) and (0, −1, −2, −3).  Here the combination is found in general by solving the
    Taylor-matching (Vandermonde) system Σ_k w_k s_k^p/p! = δ_pm, p = 0 … n−1, in exact rationals.
    Parameters: offsets (sequence of int, units of h); m (derivative order, 0 ≤ m < len(offsets)); symbolic (True →
    ``sympy.Rational`` list, False → ndarray of floats).
    Returns w with f^(m)(x₀) ≈ h^{−m} Σ_k w_k f(x₀ + s_k h) (non-dimensional weights).
    Examples: (−1, 0, 1), m = 1 → (−½, 0, ½) (centred (10.6)); m = 2 → (1, −2, 1) (10.7); (0, 1, 2), m = 1 → (−3/2, 2, −½)
    (p. 450); (0, −1, −2, −3), m = 2 → (2, −5, 4, −1) (10.149).
    Validation: V1 exact rational weights for the book stencils; V2 leading truncation term by :func:`fd_leading_error`;
    V3 observed orders on sin x.  Label: analytic, symbolic, converged.
    """
    w = _fd_weights_cached(tuple(int(s) for s in offsets), int(m))
    if symbolic:
        return list(w)
    return np.array([float(x) for x in w])


def fd_weights_float(offsets: Sequence[int], m: int = 1) -> np.ndarray:
    """Float version of :func:`fd_weights` as an ndarray.  Book: §10.2.  Label: analytic."""
    return fd_weights(offsets, m, symbolic=False)


def taylor_table_sympy(offsets: Sequence[int], order: int = 5):
    """Taylor table of a stencil: row k lists the coefficients of h^p f^(p)(x₀) in f(x₀ + s_k h), p = 0 … order.

    Book: §10.2, Eqs. (10.4)–(10.5) (the rows for s = ±1 are (10.4) and (10.5)).
    Returns sympy.Matrix (len(offsets), order + 1) with entries s_k^p/p!.  Label: symbolic.
    """
    import sympy as sp

    return sp.Matrix([[sp.Rational(s) ** p / sp.factorial(p) for p in range(order + 1)] for s in offsets])


def fd_leading_error(offsets: Sequence[int], m: int = 1, extra: int = 4):
    """Order of accuracy p and leading truncation coefficient c of the stencil of :func:`fd_weights`:
    h^{−m} Σ w_k f(x₀ + s_k h) = f^(m) + c h^p f^(m+p)(x₀) + O(h^{p+1}).

    Book: §10.2 — the O(Δx), O(Δx²) labels of (10.6)–(10.7) and of the one-sided formulas (p. 450, (10.148)–(10.149)).
    Returns (p, c) (c a sympy Rational).  Examples: forward (0, 1) → (1, ½); centred (−1, 0, 1), m = 1 → (2, 1/6);
    m = 2 → (2, 1/12).  Label: symbolic.
    """
    import sympy as sp

    w = fd_weights(offsets, m, symbolic=True)
    n = len(offsets)
    for p in range(m + 1, n + extra + m + 1):
        c = sp.nsimplify(sum(wk * sp.Rational(s) ** p / sp.factorial(p) for wk, s in zip(w, offsets)))
        if c != 0:
            return p - m, c
    return None, 0


_KIND = {
    # (kind, m) → offsets
    ("forward", 1): (0, 1), ("backward", 1): (-1, 0), ("central", 1): (-1, 0, 1), ("central2", 2): (-1, 0, 1),
    ("central", 2): (-1, 0, 1), ("onesided2", 1): (0, 1, 2), ("onesided2_backward", 1): (-2, -1, 0),
    ("forward", 2): (0, 1, 2), ("backward", 2): (-2, -1, 0), ("onesided2", 2): (0, 1, 2, 3),
    ("onesided2_backward", 2): (-3, -2, -1, 0),
}


def _kind_offsets(kind: str, m: int):
    if kind == "central2":
        m = 2
    key = (kind, int(m))
    if key not in _KIND:
        raise ValueError(f"unknown stencil {key}; choose from {sorted(_KIND)}")
    return _KIND[key], (2 if kind == "central2" else int(m))


def stencil_taylor_coefficients(kind: str = "central", m: int = 1, n_terms: int = 6) -> dict:
    """Taylor bookkeeping of a named stencil: the coefficient of h^{k−m} f^(k) in the stencil, k = 0 … n_terms − 1.

    Book: §10.2, Eqs. (10.4)–(10.7) (our D01 in table form; the explainer's term bars).  kinds: "forward", "backward",
    "central" (first derivative), "central2" (second derivative (10.7)), "onesided2" ([0, 1, 2], p. 450).
    Returns dict(offsets, weights (floats), coeffs (floats: Σ_j w_j s_j^k/k!), order (first k > m with a nonzero coefficient,
    minus m), leading (that coefficient)).  forward: coeffs (0, 1, ½, 1/6, 1/24, 1/120), order 1, leading ½; central:
    (0, 1, 0, 1/6, 0, 1/120), order 2, leading 1/6; central2: (0, 0, 1, 0, 1/12, 0), order 2, leading 1/12.
    Validation: V2 exact rationals.  Label: symbolic.
    """
    import sympy as sp

    offs, mm = _kind_offsets(kind, m)
    w = fd_weights(offs, mm, symbolic=True)
    coeffs = [sp.nsimplify(sum(wk * sp.Rational(s) ** k / sp.factorial(k) for wk, s in zip(w, offs))) for k in range(n_terms)]
    order, lead = None, 0.0
    for k in range(mm + 1, n_terms):
        if coeffs[k] != 0:
            order, lead = k - mm, float(coeffs[k])
            break
    return dict(offsets=list(offs), weights=[float(x) for x in w], coeffs=[float(c) for c in coeffs], order=order,
                leading=lead, m=mm)


def test_function(name: str, x, deriv: int = 0):
    """Smooth test functions with exact derivatives 0–5: "sin" (sin x), "exp" (eˣ), "gauss" (e^{−x²}).

    Book: §10.2 (our test fields for the stencils of (10.6)–(10.7)).  x [any unit, treated as dimensionless].
    Returns f^(deriv)(x).  Label: analytic.
    """
    x = np.asarray(x, dtype=float)
    if name == "sin":
        out = [np.sin, np.cos, lambda z: -np.sin(z), lambda z: -np.cos(z)][deriv % 4](x)
    elif name == "exp":
        out = np.exp(x)
    elif name == "gauss":
        g = np.exp(-x ** 2)
        herm = [1.0 + 0 * x, -2 * x, 4 * x ** 2 - 2, -8 * x ** 3 + 12 * x, 16 * x ** 4 - 48 * x ** 2 + 12,
                -32 * x ** 5 + 160 * x ** 3 - 120 * x]
        if deriv > 5:
            raise ValueError("deriv ≤ 5 for gauss")
        out = herm[deriv] * g
    else:
        raise ValueError("name must be sin, exp or gauss")
    return as_scalar_if_0d(out)


def stencil_derivative(func, x0: float, h: float, kind: str = "central", m: int = 1) -> float:
    """Apply a named stencil at one point: h^{−m} Σ_j w_j f(x₀ + s_j h).

    Book: §10.2, Eqs. (10.6)–(10.7), p. 450 (one-sided second order).  ``func``: a callable or a :func:`test_function` name.
    Returns the approximation (unit of f per unit x^m).  Label: analytic.
    """
    f = (lambda z: test_function(func, z)) if isinstance(func, str) else func
    offs, mm = _kind_offsets(kind, m)
    w = fd_weights_float(offs, mm)
    return float(sum(wk * f(x0 + s * h) for wk, s in zip(w, offs)) / h ** mm)  # Eq. (10.6)/(10.7)


def stencil_error(name: str, x0: float, h: float, kind: str = "central", m: int = 1) -> float:
    """Signed error (stencil − exact) of a named stencil on a :func:`test_function`.

    Book: §10.2 (the O(Δx) / O(Δx²) remainders of (10.6)–(10.7)).  E.g. sin at x₀ = 1, h = 0.1: forward −4.294e-2,
    central −9.001e-4, central2 +7.010e-4, onesided2 +1.585e-3.  Validation: V3 slopes 1 and 2.  Label: converged.
    """
    _, mm = _kind_offsets(kind, m)
    return stencil_derivative(name, x0, h, kind, m) - float(test_function(name, x0, mm))


def fd_apply(f, dx: float, offsets: Sequence[int], m: int = 1, axis: int = -1, periodic: bool = False) -> np.ndarray:
    """Apply the stencil of :func:`fd_weights` along ``axis`` of node samples with spacing ``dx``.

    Book: §10.2 (10.6)–(10.8), p. 450.  Returns an array shaped like ``f``; nodes where the stencil would leave the grid are
    NaN unless ``periodic``.  Explicit slicing — never ``np.gradient`` (first order at the edges).  Label: analytic.
    """
    f = np.moveaxis(np.asarray(f, dtype=float), axis, -1)
    w = fd_weights_float(offsets, m)
    n = f.shape[-1]
    if periodic:
        out = np.zeros_like(f)
        for wk, s in zip(w, offsets):
            out += wk * np.roll(f, -int(s), axis=-1)
    else:
        out = np.full_like(f, np.nan)
        lo, hi = -min(0, min(offsets)), n - max(0, max(offsets))
        acc = np.zeros(f.shape[:-1] + (hi - lo,))
        for wk, s in zip(w, offsets):
            acc += wk * f[..., lo + int(s): hi + int(s)]
        out[..., lo:hi] = acc
    return np.moveaxis(out / dx ** m, -1, axis)


def fd_derivative(f, dx: float, kind: str = "central", m: int = 1, axis: int = -1, periodic: bool = False) -> np.ndarray:
    """First or second derivative of node samples by a named stencil (explicit, never ``np.gradient``).

    Book: §10.2, Eq. (10.6) (forward/backward O(Δx), centred O(Δx²)), Eq. (10.7) (centred second derivative O(Δx²)),
    Eq. (10.8) (the same formulas in time), p. 450 one-sided second-order wall stencils ("onesided2" uses i, i + 1, i + 2;
    "onesided2_backward" uses i, i − 1, i − 2).
    Parameters: f (samples, any unit); dx (spacing [m] or Δt [s]); kind ("forward" | "backward" | "central" | "central2" |
    "onesided2" | "onesided2_backward"); m (1 or 2); axis; periodic.
    Returns the derivative (NaN where the stencil leaves the grid).  Validation: V1 exact on low-degree polynomials; V3
    observed orders 1 and 2.  Label: analytic, converged.
    """
    offs, mm = _kind_offsets(kind, m)
    return fd_apply(f, dx, offs, mm, axis=axis, periodic=periodic)


def mixed_derivative_onesided(f, dx: float, dy: float, side: str = "backward_x") -> np.ndarray:
    """Second-order mixed derivative ∂²f/∂x∂y, one-sided in one direction and centred in the other.

    Book: §10.5, Eq. (10.150) (block front face: backward in x over i, i−1, i−2, centred in y):
    (f_xy)_{i,j} = −1/(4ΔxΔy) [−(f_{i−2,j+1} − f_{i−2,j−1}) + 4(f_{i−1,j+1} − f_{i−1,j−1}) − 3(f_{i,j+1} − f_{i,j−1})].
    ``side``: "backward_x" (10.150), "forward_x", "backward_y", "forward_y".  Layout f[j, i]; dx, dy [m].
    Returns f_xy (NaN where the stencil leaves the grid).  Validation: V1 exact on bi-quadratics; V3 order 2.
    Label: analytic, converged.
    """
    f = np.asarray(f, dtype=float)
    if side in ("backward_x", "forward_x"):
        dyc = fd_apply(f, dy, (-1, 0, 1), 1, axis=0)
        return fd_apply(dyc, dx, (0, -1, -2) if side == "backward_x" else (0, 1, 2), 1, axis=1)  # Eq. (10.150)
    if side in ("backward_y", "forward_y"):
        dxc = fd_apply(f, dx, (-1, 0, 1), 1, axis=1)
        return fd_apply(dxc, dy, (0, -1, -2) if side == "backward_y" else (0, 1, 2), 1, axis=0)
    raise ValueError("side must be backward_x, forward_x, backward_y or forward_y")


# ======================================================================================================================
# The model problem (10.1)–(10.3) and its schemes (10.9)–(10.13), (10.29)
# ======================================================================================================================
def transport_rhs(T, u: float, D: float, dx: float, periodic: bool = True) -> np.ndarray:
    """Semi-discrete right-hand side −u T_x + D T_xx of Eq. (10.1) with the centred stencils (10.6)–(10.7).

    Book: §10.2, Eq. (10.1) ∂T/∂t + u ∂T/∂x = D ∂²T/∂x² (constant u [m/s], D [m²/s]).  Node-based grid, spacing dx [m].
    Returns dT/dt [unit of T per s]; end nodes NaN unless periodic.  Label: analytic.
    """
    T = np.asarray(T, dtype=float)
    return (-u * fd_derivative(T, dx, "central", 1, periodic=periodic)
            + D * fd_derivative(T, dx, "central2", 2, periodic=periodic))  # Eq. (10.1)


def advected_gaussian(x, t, u: float, D: float, x0: float = 0.3, s0: float = 0.05, L: float | None = 1.0,
                      n_images: int = 3):
    """Exact solution of Eq. (10.1) for a Gaussian initial condition on a periodic domain: it translates and spreads.

    Book: §10.2, Eq. (10.1) (test field of ours, the heat kernel in a moving frame; not printed):
    T(x, t) = (s₀/s) Σ_{m=−n_images}^{n_images} exp[−(x − x₀ − ut − mL)²/(2s²)],  s² = s₀² + 2Dt.
    Parameters: x [m], t [s] (≥ 0), u [m/s], D [m²/s] (≥ 0), x0 [m], s0 [m] (> 0), L [m] (period; None → infinite line),
    n_images (images on each side; 3 is round-off exact while s ≪ L).
    Returns T (dimensionless amplitude, 1 at t = 0).  Validation: V2 satisfies (10.1) (sympy); V4 ∫T dx conserved.
    Label: analytic.
    """
    require_positive("s0", s0)
    require_nonnegative("D", D)
    x = np.asarray(x, dtype=float)
    s = np.sqrt(s0 ** 2 + 2.0 * D * t)
    xc = x0 + u * t
    if L is None:
        return as_scalar_if_0d((s0 / s) * np.exp(-((x - xc) ** 2) / (2.0 * s ** 2)))
    if L is not None:
        xc = (xc % L)
    out = np.zeros_like(x)
    for mm in range(-int(n_images), int(n_images) + 1):
        out = out + np.exp(-((x - xc - mm * L) ** 2) / (2.0 * s ** 2))
    return as_scalar_if_0d((s0 / s) * out)


def ftcs_coefficients(u: float, D: float, dx: float, dt: float):
    """The two numbers of the FTCS scheme, α = uΔt/(2Δx) and β = DΔt/Δx².

    Book: §10.2, Eq. (10.11).  u [m/s], D [m²/s] (≥ 0), dx [m], dt [s].  Returns (alpha, beta) (non-dimensional);
    C = uΔt/Δx = 2α.  Example (0.1, 0.01, 0.1, 0.2) → (0.1, 0.2).  Label: analytic.
    """
    require_positive("dx", dx)
    require_positive("dt", dt)
    require_nonnegative("D", D)
    alpha = u * dt / (2.0 * dx)  # Eq. (10.11)
    beta = D * dt / dx ** 2  # Eq. (10.11)
    return as_scalar_if_0d(alpha), as_scalar_if_0d(beta)


def cfl_number(u, dt: float, dx: float):
    """Courant number C = |u|Δt/Δx — the CFL condition (10.30) requires C ≤ 1.

    Book: §10.2, Eq. (10.30) (printed uΔt/Δx ≤ 1, i.e. for u > 0; slip R11: |u| — a particle must not cross more than one
    cell per step, whichever way it moves).  u [m/s], dt [s], dx [m].  Returns C.  Label: analytic.
    """
    return as_scalar_if_0d(np.abs(np.asarray(u, dtype=float)) * dt / dx)  # Eq. (10.30)


def diffusion_number(D, dt: float, dx: float):
    """β = DΔt/Δx² (Eq. (10.11)); pure diffusion FTCS needs β ≤ ½ (Eq. (10.28)).  Label: analytic."""
    return as_scalar_if_0d(np.asarray(D, dtype=float) * dt / dx ** 2)


def cell_peclet(u, dx: float, D: float, signed: bool = False):
    """Cell Péclet (cell Reynolds) number R_cell = |u|Δx/D; wiggle-free centred differencing needs R_cell ≤ 2.

    Book: §10.2, Eq. (10.31); §10.4 below (10.91): R_cell = R/n, R = uL/D (10.87).  u [m/s], dx [m], D [m²/s] (> 0);
    ``signed`` keeps the sign of u.  Example (1, 0.01, 0.005) → 2.  Label: analytic.
    """
    require_positive("D", D)
    uu = np.asarray(u, dtype=float)
    return as_scalar_if_0d((uu if signed else np.abs(uu)) * dx / D)  # Eq. (10.31)


def _neighbours(T, periodic):
    if periodic:
        return np.roll(T, -1), np.roll(T, 1)
    Tp = np.empty_like(T)
    Tm = np.empty_like(T)
    Tp[:-1], Tp[-1] = T[1:], np.nan
    Tm[1:], Tm[0] = T[:-1], np.nan
    return Tp, Tm


def _explicit(T, Tp, Tm, alpha, beta, scheme):
    lap = Tp - 2.0 * T + Tm
    if scheme == "ftcs":
        return T - alpha * (Tp - Tm) + beta * lap  # Eq. (10.10)
    if scheme in ("upwind", "upwind_printed"):
        C = 2.0 * alpha
        if scheme == "upwind_printed" or C >= 0:
            return T - C * (T - Tm) + beta * lap  # Eq. (10.29) (+ centred diffusion): information from the left
        return T - C * (Tp - T) + beta * lap  # R11: u < 0, the upwind neighbour is i + 1
    if scheme == "lax_wendroff":
        return T - alpha * (Tp - Tm) + (2.0 * alpha ** 2 + beta) * lap  # = MacCormack for linear advection (D18)
    raise ValueError(scheme)


def _implicit_matrix(n: int, alpha: float, beta: float, theta: float, periodic: bool, right: str | None):
    """Sparse matrix of the implicit part I + θ[α(T₊ − T₋) − β(T₊ − 2T + T₋)] (BTCS θ = 1, Crank–Nicolson θ = ½)."""
    import scipy.sparse as sps

    lo = -theta * (alpha + beta)  # coefficient of T_{i−1} — Eq. (10.13) with θ = 1
    di = 1.0 + 2.0 * theta * beta
    up = theta * (alpha - beta)  # coefficient of T_{i+1}
    A = sps.lil_matrix((n, n))
    A.setdiag(np.full(n, di))
    A.setdiag(np.full(n - 1, lo), -1)
    A.setdiag(np.full(n - 1, up), 1)
    if periodic:
        A[0, n - 1] = lo
        A[n - 1, 0] = up
    else:
        A[0, :] = 0.0
        A[0, 0] = 1.0  # Dirichlet row at x = 0
        if right == "neumann":  # ghost T_N = T_{N−2} + 2Δx q: the ghost's coefficient moves onto T_{N−2}
            A[n - 1, n - 2] = lo + up
        else:
            A[n - 1, :] = 0.0
            A[n - 1, n - 1] = 1.0
    return A.tocsc()


def transport_1d_step(T, alpha: float, beta: float, scheme: str = "ftcs", periodic: bool = True, g: float | None = None,
                      q: float | None = None, dx: float | None = None, T_L: float | None = None, _lu=None):
    """Advance Eq. (10.1) by one time step with a named scheme.

    Book: §10.2 — FTCS Eq. (10.10), BTCS Eq. (10.13) (one tridiagonal solve per step), first-order upwind Eq. (10.29) (here
    with the centred diffusion term β(T₊ − 2T + T₋) added; side from sign(α)), "upwind_printed" (the printed stencil
    T_i − T_{i−1} for any sign of u — slip R11, wrong for u < 0), Crank–Nicolson "cn" (ours, θ = ½) and Lax–Wendroff
    "lax_wendroff" (ours: what MacCormack (10.101)–(10.102) reduces to for linear advection).
    Parameters
    ----------
    T : node values at t_n.   alpha, beta : Eq. (10.11) (α carries the sign of u).   scheme : see above.
    periodic : wrap-around grid (default).  Non-periodic: Dirichlet T(0) = g (10.2) (default T[0]) and, at x = L, a Neumann
        slope ∂T/∂x = q (10.2) by a ghost node T_N = T_{N−2} + 2Δx q (second order; needs ``dx``) — or, if ``q`` is None, a
        Dirichlet value ``T_L`` (default T[-1]).
    Returns the node values at t_{n+1}.  Example: FTCS on [0, 0, 1, 0, 0] (periodic) with (α, β) = (0.1, 0.2) → [0, 0.1,
    0.6, 0.3, 0]; upwind α = 0.25 → [0, 0, 0.5, 0.5, 0].
    Assumptions: constant u, D; uniform node-based grid.
    Validation: V1 u = 0 reproduces ``core.diffusion.ftcs_diffusion_1d``; V3 orders (FTCS (2, 1), BTCS (2, 1), CN (2, 2), upwind
    1); V7 u → −u mirrors the solution.  Label: analytic, converged.
    """
    T = np.asarray(T, dtype=float)
    s = scheme.lower()
    if s not in SCHEMES_1D:
        raise ValueError(f"scheme must be one of {SCHEMES_1D}")
    n = T.size
    gval = float(T[0] if g is None else g)
    neumann = (not periodic) and (q is not None)
    if neumann and dx is None:
        raise ValueError("a Neumann right end needs dx")
    TLval = float(T[-1] if T_L is None else T_L)
    if s in ("ftcs", "upwind", "upwind_printed", "lax_wendroff"):
        Tp, Tm = _neighbours(T, periodic)
        if neumann:
            Tp[-1] = T[-2] + 2.0 * dx * q  # ghost node (T_N − T_{N−2})/(2Δx) = q
        Tn = _explicit(T, Tp, Tm, alpha, beta, s)
        if not periodic:
            Tn[0] = gval
            if not neumann:
                Tn[-1] = TLval
        return Tn
    theta = 1.0 if s == "btcs" else 0.5
    if theta < 1.0:
        Tp, Tm = _neighbours(T, periodic)
        if neumann:
            Tp[-1] = T[-2] + 2.0 * dx * q
        rhs = _explicit(T, Tp, Tm, 0.5 * alpha, 0.5 * beta, "ftcs")  # explicit half of Crank–Nicolson
    else:
        rhs = T.copy()
    if not periodic:
        rhs[0] = gval
        if neumann:
            rhs[-1] -= theta * (alpha - beta) * 2.0 * dx * q  # the new ghost's 2Δx q moves to the right-hand side
        else:
            rhs[-1] = TLval
    if _lu is None:
        from scipy.sparse.linalg import splu

        _lu = splu(_implicit_matrix(n, alpha, beta, theta, periodic, "neumann" if neumann else "dirichlet"))
    return _lu.solve(rhs)  # Eq. (10.13): T^n + α(T^n_{i+1} − T^n_{i−1}) − β(T^n_{i+1} − 2T^n_i + T^n_{i−1}) = T^{n−1}


def _explicit_stable(s, alpha, beta):
    if s == "ftcs":
        return ftcs_stable(alpha, beta)
    if s == "upwind":
        return abs(2 * alpha) + 2 * beta <= 1.0 + 1e-12
    if s == "upwind_printed":
        return 0 <= 2 * alpha and 2 * alpha + 2 * beta <= 1.0 + 1e-12
    if s == "lax_wendroff":
        return 4 * alpha ** 2 + 2 * beta <= 1.0 + 1e-12 and beta >= 0
    return True


def solve_transport_1d(T0, x, u: float, D: float, dt: float, nsteps: int, scheme: str = "ftcs", g: float | None = None,
                       q: float | None = None, periodic: bool = False, save_every: int | None = None,
                       check_stability: bool = True, growth_cap: float = 1e6, T_L: float | None = None) -> dict:
    """March Eq. (10.1) with boundary conditions (10.2) and initial condition (10.3) by a named scheme.

    Book: §10.2, Eqs. (10.1)–(10.3), (10.10), (10.13), (10.29); stability (10.27) (FTCS) and (10.30) (upwind, with |u|).
    Parameters
    ----------
    T0 : initial node values (10.3).   x : node positions [m], uniform (periodic: no repeated end node).
    u [m/s], D [m²/s], dt [s], nsteps.   scheme : see :func:`transport_1d_step`.
    g : Dirichlet value at x = 0 (default T0[0]).   q : Neumann slope at x = L (None → Dirichlet end ``T_L``, default T0[-1]).
    periodic : wrap-around grid.   save_every : store every k-th state (None → the first and last only).
    check_stability : raise ValueError before marching if an explicit scheme violates its von Neumann limit (FTCS (10.27);
        upwind |C| + 2β ≤ 1; Lax–Wendroff C² + 2β ≤ 1).  False only for blow-up demonstrations.
    growth_cap : stop as soon as max|T| > growth_cap × max|T0| (no NaN plots); the step is reported.
    Returns dict(x, t (final time [s]), T (final), history (list of arrays), times (list), blew_up, step_blown (int or None),
    max_abs (per saved state), alpha, beta, courant).
    Validation: V3 observed orders; V7 blow-up just outside (10.27)/(10.30), bounded just inside.  Label: converged.
    """
    x = np.asarray(x, dtype=float)
    T = np.array(T0, dtype=float)
    dx = float(x[1] - x[0])
    alpha, beta = ftcs_coefficients(u, D, dx, dt)
    s = scheme.lower()
    if check_stability and not _explicit_stable(s, alpha, beta):
        raise ValueError(f"{s} unstable for alpha = {alpha:.4g}, beta = {beta:.4g} (C = {2 * alpha:.4g}); reduce dt or pass "
                         "check_stability=False to demonstrate the blow-up")
    gval = float(T[0] if g is None else g)
    TLval = float(T[-1] if T_L is None else T_L)
    lu = None
    if s in ("btcs", "cn"):
        from scipy.sparse.linalg import splu

        lu = splu(_implicit_matrix(T.size, alpha, beta, 1.0 if s == "btcs" else 0.5, periodic,
                                   "neumann" if (q is not None and not periodic) else "dirichlet"))
    scale = max(float(np.max(np.abs(T))), 1e-300)
    hist, times = [T.copy()], [0.0]
    blew, bstep, k = False, None, 0
    for k in range(1, nsteps + 1):
        T = transport_1d_step(T, alpha, beta, s, periodic, gval, q, dx, TLval, _lu=lu)
        if save_every and k % save_every == 0:
            hist.append(T.copy())
            times.append(k * dt)
        if not np.all(np.isfinite(T)) or np.max(np.abs(T)) > growth_cap * scale:
            blew, bstep = True, k
            break
    if times[-1] != k * dt:
        hist.append(T.copy())
        times.append(k * dt)
    return dict(x=x, t=k * dt, T=T, history=hist, times=times, blew_up=blew, step_blown=bstep,
                max_abs=[float(np.max(np.abs(h))) for h in hist], alpha=alpha, beta=beta, courant=2 * alpha)


# ======================================================================================================================
# Convergence (10.14)–(10.15), Richardson, grid convergence index
# ======================================================================================================================
def error_norm(num, exact, kind: str = "rms") -> float:
    """A measure ‖e‖ of the solution error e_i^n = T_i^n − T(x_i, t_n).

    Book: §10.2, Eq. (10.14) (error) and Eq. (10.15) ("the root mean square of the solution error on all the grid points").
    kind: "rms" (the book's), "max", "l1" (mean absolute).  Returns the same unit as the inputs.  Label: analytic.
    """
    e = np.asarray(num, dtype=float) - np.asarray(exact, dtype=float)  # Eq. (10.14)
    if kind == "rms":
        return float(np.sqrt(np.mean(e ** 2)))
    if kind == "max":
        return float(np.max(np.abs(e)))
    if kind == "l1":
        return float(np.mean(np.abs(e)))
    raise ValueError("kind must be rms, max or l1")


def observed_order(h, err) -> float:
    """Least-squares slope of log(err) vs log(h) — a rate of Eq. (10.15) ‖e‖ ≤ K Δx^a Δt^b.  Errors ≤ 1e-14 are dropped.
    Book: §10.2.  Label: analytic."""
    h = np.asarray(h, dtype=float)
    err = np.asarray(err, dtype=float)
    keep = (err > 1e-14) & np.isfinite(err)
    if keep.sum() < 2:
        raise ValueError("need at least two usable error values")
    return float(np.polyfit(np.log(h[keep]), np.log(err[keep]), 1)[0])


def convergence_rates(error_of_h: Callable[[float], float], h_list: Sequence[float]) -> dict:
    """Generic refinement study: run ``error_of_h(h)`` for each spacing.  Returns dict(h, err, order, pairwise).
    Book: §10.2, (10.15).  Label: converged."""
    hs = [float(h) for h in h_list]
    es = [float(error_of_h(h)) for h in hs]
    pw = [float(np.log(es[i] / es[i + 1]) / np.log(hs[i] / hs[i + 1])) for i in range(len(hs) - 1)]
    return dict(h=hs, err=es, order=observed_order(hs, es), pairwise=pw)


def convergence_study(scheme: str = "ftcs", n_list=(20, 40, 80, 160), u: float = 0.5, D: float = 0.01, t_end: float = 0.2,
                      dt_rule: str = "diffusive", L: float = 1.0, beta: float = 0.25, C: float = 0.5) -> dict:
    """Observed convergence rate of a scheme for (10.1) on the periodic advected Gaussian (:func:`advected_gaussian`).

    Book: §10.2, Eqs. (10.14)–(10.15) ‖e‖ ≤ K Δx^a Δt^b.  ``dt_rule``: "diffusive" (Δt = βΔx²/D — measures a + 2b, so FTCS
    shows 2) or "advective" (Δt = CΔx/|u| — measures min(a, b), so upwind shows 1).  The step count is rounded so that the
    run ends exactly at ``t_end``.
    Returns dict(dx, dt, err (rms at t_end), order (observed, in Δx), pairwise).  Asymptotic orders: FTCS diffusive 2.0,
    upwind advective 1.0, BTCS diffusive 2.0 (± 0.15) — reached only on finer grids for the last two: on the default
    n_list = (20, 40, 80, 160) the fitted orders are pre-asymptotic (upwind advective ≈ 0.76, BTCS diffusive ≈ 1.80, the
    Gaussian of width 0.05 being barely resolved at n = 20); use e.g. n_list = (160, 320, 640, 1280) to see 1 and 2.
    Label: converged.
    """
    dxs, dts, errs = [], [], []
    for n in n_list:
        dx = L / n
        x = np.arange(n) * dx
        if dt_rule == "diffusive":
            dt0 = beta * dx ** 2 / D
        else:  # advective, but never beyond the diffusion limit β ≤ beta (otherwise FTCS-type diffusion blows up on fine grids)
            dt0 = min(C * dx / max(abs(u), 1e-300), beta * dx ** 2 / D if D > 0 else np.inf)
        ns = max(1, int(np.ceil(t_end / dt0)))
        dt = t_end / ns
        r = solve_transport_1d(advected_gaussian(x, 0.0, u, D, L=L), x, u, D, dt, ns, scheme, periodic=True,
                               check_stability=False)
        errs.append(error_norm(r["T"], advected_gaussian(x, t_end, u, D, L=L)))
        dxs.append(dx)
        dts.append(dt)
    pw = [float(np.log(errs[i] / errs[i + 1]) / np.log(dxs[i] / dxs[i + 1])) for i in range(len(dxs) - 1)]
    return dict(dx=dxs, dt=dts, err=errs, order=observed_order(dxs, errs), pairwise=pw)


def richardson_extrapolate(f_coarse: float, f_fine: float, r: float = 2.0, p: float = 2.0) -> float:
    """Richardson extrapolation f₀ ≈ f_fine + (f_fine − f_coarse)/(r^p − 1) (our D23 from the error model (10.15)).
    Book: §10.5–10.6 (grid independence).  Label: analytic."""
    return float(f_fine + (f_fine - f_coarse) / (r ** p - 1.0))


def grid_convergence_index(f1: float, f2: float, f3: float, r: float = 2.0, p: float | None = None, Fs: float = 1.25) -> dict:
    """Observed order, Richardson extrapolate and Roache's grid convergence index from three grids (f1 finest).

    Book: §10.5 (Fig. 10.13 grid convergence) and §10.6 ("a mesh-refinement test is normally a must"); our D23 from the error
    model (10.15) f_h = f₀ + K h^p: p = ln[(f₃ − f₂)/(f₂ − f₁)]/ln r, f_ext = f₁ + (f₁ − f₂)/(r^p − 1),
    GCI_fine = F_s |(f₂ − f₁)/f₁|/(r^p − 1), GCI_coarse = r^p GCI_fine, asymptotic ratio GCI_coarse/(r^p GCI_fine) (≈ 1).
    ``p`` given → used instead of the observed order.  Returns dict(p, f_ext, gci_fine, gci_coarse, asymptotic_ratio,
    monotone (differences of one sign)).  Example f = 1 + h² at h = 0.025, 0.05, 0.1 → p = 2, f_ext = 1.
    Degenerate cases: a zero difference, or p = 0 (|f₃ − f₂| = |f₂ − f₁|, so r^p − 1 = 0) → f_ext and the GCIs NaN
    (p is NaN for a zero difference and the observed 0.0 for the p = 0 case), no exception.
    Validation: V1 exact on f₀ + K h^p.  Label: analytic.
    """
    d21, d32 = f2 - f1, f3 - f2
    monotone = bool(d21 * d32 > 0)
    if p is None:
        p = float(np.log(abs(d32 / d21)) / np.log(r)) if d21 != 0 and d32 != 0 else np.nan
    rp = r ** p
    if np.isfinite(p) and abs(rp - 1.0) < 1e-12:
        # p = 0 (|f₃ − f₂| = |f₂ − f₁|, e.g. equal and opposite differences): no convergence, r^p − 1 = 0 —
        # the extrapolation and the GCI are undefined, not infinite
        return dict(p=float(p), f_ext=np.nan, gci_fine=np.nan, gci_coarse=np.nan, asymptotic_ratio=np.nan,
                    monotone=monotone)
    f_ext = float(f1 + (f1 - f2) / (rp - 1.0)) if np.isfinite(p) else np.nan
    gf = float(Fs * abs(d21 / f1) / (rp - 1.0)) if np.isfinite(p) and f1 != 0 else np.nan
    gc = float(Fs * abs(d32 / f2) / (rp - 1.0)) if np.isfinite(p) and f2 != 0 else np.nan
    return dict(p=float(p), f_ext=f_ext, gci_fine=gf, gci_coarse=gc,
                asymptotic_ratio=float(gc / (rp * gf)) if gf and np.isfinite(gc) else np.nan, monotone=monotone)


def richardson_three(f1: float, f2: float, f3: float, r: float = 2.0) -> dict:
    """dict(p, f_ext) part of :func:`grid_convergence_index` (f1 finest) — scalar-callable for explainer parity rows.
    Book: §10.5 (Fig. 10.13).  Label: analytic."""
    g = grid_convergence_index(f1, f2, f3, r)
    return dict(p=g["p"], f_ext=g["f_ext"])


# ======================================================================================================================
# Consistency: truncation errors (modified equations) — (10.16)–(10.17), (10.94)
# ======================================================================================================================
def _taylor_shift(a, b, order, dsym):
    import sympy as sp

    return sum(a ** p * b ** q / (sp.factorial(p) * sp.factorial(q)) * dsym[(p, q)]
               for p in range(order + 1) for q in range(order + 1 - p))


def truncation_error_sympy(scheme: str = "ftcs", order: int = 5, substitute_pde: bool = False) -> dict:
    """Truncation error E of a scheme: what is left when the exact solution is put into the scheme, minus the PDE.

    Book: §10.2, Eqs. (10.16)–(10.17) (FTCS: [T_t + uT_x − DT_xx]_i^n + E_i^n = 0 with E = (Δt/2)T_tt + u(Δx²/6)T_xxx −
    D(Δx²/12)T_xxxx + O(Δt², Δx⁴)); §10.4, Eq. (10.94) (steady upwind: E = −(uΔx/2)T_xx, i.e. uT_x = D(1 + 0.5R_cell)T_xx).
    Method (our D03/D17): substitute the two-variable Taylor series of T about (x_i, t_n) into the stencil, divide by the
    time-step factor, subtract the PDE; keep the lowest power of Δt and the lowest power of Δx (what (10.17) writes out).
    schemes: "ftcs" (10.9), "btcs" (10.12), "upwind" (10.29 + centred diffusion, u > 0), "lax_wendroff",
    "upwind_steady" (10.93), "central_steady" (10.90).
    substitute_pde: also replace time derivatives by space derivatives using the PDE (T_t = −uT_x + DT_xx, repeatedly).
    Returns dict(E (sympy expression), E_full, order_t, order_x (powers of the leading terms), modified_equation (sympy Eq
    for the steady schemes: u T_x = (D + D_num) T_xx), symbols, leading).
    Validation: V2 reproduces (10.17) and (10.94).  Label: symbolic.
    """
    import sympy as sp

    dt, dx, u, D = sp.symbols("Delta_t Delta_x u D", positive=True)
    names, dsym = {}, {}
    for p in range(order + 1):
        for q in range(order + 1 - p):
            nm = "T" + ("_" + "x" * p + "t" * q if p + q else "")
            dsym[(p, q)] = sp.Symbol(nm)
            names[nm] = dsym[(p, q)]
    T = lambda a, b: _taylor_shift(a, b, order, dsym)  # noqa: E731
    lap = lambda: (T(dx, 0) - 2 * T(0, 0) + T(-dx, 0)) / dx ** 2  # noqa: E731
    s = scheme.lower()
    steady = s.endswith("_steady")
    if s == "ftcs":
        S = (T(0, dt) - T(0, 0)) / dt + u * (T(dx, 0) - T(-dx, 0)) / (2 * dx) - D * lap()  # Eq. (10.9)
    elif s == "btcs":
        S = (T(0, 0) - T(0, -dt)) / dt + u * (T(dx, 0) - T(-dx, 0)) / (2 * dx) - D * lap()  # Eq. (10.12)
    elif s == "upwind":
        S = (T(0, dt) - T(0, 0)) / dt + u * (T(0, 0) - T(-dx, 0)) / dx - D * lap()  # Eq. (10.29)
    elif s == "lax_wendroff":
        S = (T(0, dt) - T(0, 0)) / dt + u * (T(dx, 0) - T(-dx, 0)) / (2 * dx) - (u ** 2 * dt / 2 + D) * lap()
    elif s == "upwind_steady":
        S = u * (T(0, 0) - T(-dx, 0)) / dx - D * lap()  # Eq. (10.93) divided by Δx²/D
    elif s == "central_steady":
        S = u * (T(dx, 0) - T(-dx, 0)) / (2 * dx) - D * lap()  # Eq. (10.90)
    else:
        raise ValueError("unknown scheme")
    pde = u * dsym[(1, 0)] - D * dsym[(2, 0)] + (0 if steady else dsym[(0, 1)])
    E_full = sp.expand(S - pde)
    terms = sp.Add.make_args(E_full)
    degs = [(sp.degree(tm, dt) if tm.has(dt) else 0, sp.degree(tm, dx) if tm.has(dx) else 0) for tm in terms]
    pt = min([a for a, b in degs if a > 0 and b == 0], default=None)
    px = min([b for a, b in degs if b > 0 and a == 0], default=None)
    keep = [tm for tm, (a, b) in zip(terms, degs)
            if (a > 0 and b == 0 and a == pt) or (b > 0 and a == 0 and b == px) or (a == 0 and b == 0)]
    E = sp.Add(*keep)
    if substitute_pde and not steady:
        rep = {}
        for (p, q), sym in dsym.items():
            if q > 0:
                poly = {p: sp.Integer(1)}
                for _ in range(q):
                    new = {}
                    for k, c in poly.items():
                        new[k + 1] = new.get(k + 1, 0) - u * c
                        new[k + 2] = new.get(k + 2, 0) + D * c
                    poly = new
                rep[sym] = sum(c * dsym.get((k, 0), sp.Symbol("T_" + "x" * k)) for k, c in poly.items())
        E = sp.expand(E.subs(rep))
    mod = None
    if steady:
        Dnum = sp.simplify(-E.coeff(dsym[(2, 0)])) if E.has(dsym[(2, 0)]) else sp.Integer(0)
        mod = sp.Eq(u * dsym[(1, 0)], sp.factor(D + Dnum) * dsym[(2, 0)])
    lead = {str(k): v for k, v in sp.collect(E, list(E.free_symbols - {dt, dx, u, D}), evaluate=False).items()}
    return dict(E=E, E_full=E_full, order_t=pt, order_x=px, modified_equation=mod,
                symbols=dict(dt=dt, dx=dx, u=u, D=D, **names), leading=lead)


@lru_cache(maxsize=4)
def _gaussian_derivs():
    import sympy as sp

    x, t, u, D, x0, s0 = sp.symbols("x t u D x0 s0", real=True)
    s2 = s0 ** 2 + 2 * D * t
    T = s0 / sp.sqrt(s2) * sp.exp(-(x - x0 - u * t) ** 2 / (2 * s2))
    return {k: sp.lambdify((x, t, u, D, x0, s0), e, "numpy") for k, e in
            dict(T=T, T_tt=sp.diff(T, t, 2), T_xx=sp.diff(T, x, 2), T_xxx=sp.diff(T, x, 3), T_xxxx=sp.diff(T, x, 4)).items()}


def truncation_terms(u: float, D: float, dx: float, dt: float, x, t: float = 0.0, x0: float = 0.3, s0: float = 0.05,
                     L: float | None = 1.0, n_images: int = 3) -> dict:
    """The three leading truncation-error terms of FTCS, Eq. (10.17), evaluated on the exact periodic advected Gaussian.

    Book: §10.2, Eq. (10.17): E ≈ (Δt/2)T_tt + u(Δx²/6)T_xxx − D(Δx²/12)T_xxxx.  ``measured`` is the one-step FTCS residual of
    the exact solution, [T(x, t + Δt) − T]/Δt + u[T(x + Δx) − T(x − Δx)]/(2Δx) − D[T(x + Δx) − 2T + T(x − Δx)]/Δx² — the
    truncation error itself; it approaches ``total`` as Δx, Δt → 0.
    Parameters: u [m/s], D [m²/s], dx [m], dt [s], x [m] (float or array), t [s], x0, s0, L (period, images summed).
    Returns dict(time, conv, diff, total, measured) — floats for a float x [unit of T per s].
    Validation: V2 agrees with :func:`truncation_error_sympy`; V3 measured − total = O(Δt², Δx⁴).  Label: symbolic, converged.
    """
    f = _gaussian_derivs()
    xx = np.asarray(x, dtype=float)
    images = [0] if L is None else range(-n_images, n_images + 1)

    def val(key, xv, tv):
        return sum(f[key](xv - m * (L or 0.0), tv, u, D, x0, s0) for m in images)

    time = 0.5 * dt * val("T_tt", xx, t)
    conv = u * dx ** 2 / 6.0 * val("T_xxx", xx, t)
    diff = -D * dx ** 2 / 12.0 * val("T_xxxx", xx, t)
    T0 = val("T", xx, t)
    meas = ((val("T", xx, t + dt) - T0) / dt + u * (val("T", xx + dx, t) - val("T", xx - dx, t)) / (2 * dx)
            - D * (val("T", xx + dx, t) - 2 * T0 + val("T", xx - dx, t)) / dx ** 2)
    out = dict(time=time, conv=conv, diff=diff, total=time + conv + diff, measured=meas)
    return {k: as_scalar_if_0d(v) for k, v in out.items()}


def ode_scheme_errors(scheme: str = "forward", lam: float = -1.0, t_end: float = 1.0, n_list=(20, 40, 80, 160),
                      alpha_t: float = 1.0, beta_t: float = 0.0) -> dict:
    """Errors of a time-difference formula applied to y′ = λy, y(0) = 1 (exact e^{λt}), on several step sizes.

    Book: §10.2, Eq. (10.8) — "forward" (explicit Euler, O(Δt)), "backward" (implicit Euler, O(Δt)), "leapfrog" (centred
    O(Δt²), started from the exact y(Δt)), "trapezoidal"; §10.5, Eq. (10.163) ("alpha_beta": α(y^{n+1} − y^n)/Δt − β y′(t_n)
    = λy^{n+1}; α = 1, β = 0 backward Euler; α = 2, β = 1 trapezoidal, second order).
    Returns dict(dt, err (at t_end), order).  Label: converged.
    """
    errs, dts = [], []
    for n in n_list:
        h = t_end / n
        y = 1.0
        if scheme == "forward":
            for _ in range(n):
                y = y + h * lam * y  # Eq. (10.8) forward
        elif scheme == "backward":
            for _ in range(n):
                y = y / (1 - h * lam)  # Eq. (10.8) backward
        elif scheme == "leapfrog":
            y0, y1 = 1.0, np.exp(lam * h)
            for _ in range(n - 1):
                y0, y1 = y1, y0 + 2 * h * lam * y1  # Eq. (10.8) centred
            y = y1
        elif scheme == "trapezoidal":
            for _ in range(n):
                y = y * (1 + 0.5 * h * lam) / (1 - 0.5 * h * lam)
        elif scheme == "alpha_beta":
            dydt = lam * y
            for _ in range(n):
                y1 = (alpha_t * y / h + beta_t * dydt) / (alpha_t / h - lam)  # Eq. (10.163) with y′(t_{n+1}) = λ y^{n+1}
                dydt = alpha_t * (y1 - y) / h - beta_t * dydt
                y = y1
        else:
            raise ValueError("scheme must be forward, backward, leapfrog, trapezoidal or alpha_beta")
        errs.append(abs(y - np.exp(lam * t_end)))
        dts.append(h)
    return dict(dt=dts, err=errs, order=observed_order(dts, errs))


def ode_scheme_order(scheme: str = "forward", lam: float = -1.0, t_end: float = 1.0, n_list=(20, 40, 80, 160)) -> float:
    """Observed order of a time difference of Eq. (10.8) on y′ = λy — 1 (forward, backward), 2 (leapfrog, trapezoidal).
    Book: §10.2, Eq. (10.8); §10.5, Eq. (10.163) (trapezoidal = α = 2, β = 1).  Returns float.  Label: converged."""
    return float(ode_scheme_errors(scheme, lam, t_end, n_list)["order"])


# ======================================================================================================================
# Stability: von Neumann (10.18)–(10.28), upwind and CFL (10.29)–(10.30)
# ======================================================================================================================
def propagate_error(xi0, alpha: float, beta: float, nsteps: int, scheme: str = "ftcs", periodic: bool = True,
                    growth_cap: float = 1e6) -> dict:
    """March a disturbance ξ with the homogeneous scheme — the error equation (10.19).

    Book: §10.2, Eq. (10.18) ξ = T − T̄ (difference of two solutions of the discrete system) and Eq. (10.19)
    ξ_i^{n+1} = (α + β)ξ_{i−1}^n + (1 − 2β)ξ_i^n + (β − α)ξ_{i+1}^n (FTCS); any scheme of :func:`transport_1d_step`
    (non-periodic: ξ = 0 at both ends).  Stops when max|ξ| > growth_cap·max|ξ₀|.
    Returns dict(xi (final), max_abs (per step, ndarray), history (steps + 1, N), steps_done, blew_up).
    Validation: V1 a single Fourier mode is multiplied by G(θ) of :func:`amplification_factor` every step.  Label: analytic.
    """
    xi = np.array(xi0, dtype=float)
    scale = max(float(np.max(np.abs(xi))), 1e-300)
    hist = [xi.copy()]
    lu = None
    if scheme in ("btcs", "cn"):
        from scipy.sparse.linalg import splu

        lu = splu(_implicit_matrix(xi.size, alpha, beta, 1.0 if scheme == "btcs" else 0.5, periodic, "dirichlet"))
    blew = False
    for _ in range(nsteps):
        xi = transport_1d_step(xi, alpha, beta, scheme, periodic, 0.0, None, None, 0.0, _lu=lu)  # Eq. (10.19)
        hist.append(xi.copy())
        if np.max(np.abs(xi)) > growth_cap * scale:
            blew = True
            break
    H = np.array(hist)
    return dict(xi=xi, max_abs=np.max(np.abs(H), axis=1), history=H, steps_done=len(hist) - 1, blew_up=blew)


def fourier_mode(x, k: float, convention: str = "book"):
    """One Fourier component of the error: e^{iπkx} (book, Eq. (10.21): k counts half-waves per unit length, θ = kπΔx) or
    e^{ikx} (``convention="standard"``, θ = kΔx).  Book: §10.2, Eqs. (10.20)–(10.22).  Returns complex array.
    Label: analytic."""
    x = np.asarray(x, dtype=float)
    if convention == "book":
        return np.exp(1j * np.pi * k * x)  # Eq. (10.21)
    if convention == "standard":
        return np.exp(1j * k * x)
    raise ValueError("convention must be 'book' or 'standard'")


def amplification_factor(theta, alpha: float = 0.0, beta: float = 0.0, scheme: str = "ftcs"):
    """Von Neumann amplification factor G = g^{n+1}/g^n of a linear scheme for Eq. (10.1).

    Book: §10.2, Eq. (10.24) (FTCS): G = (α + β)e^{−iθ} + (1 − 2β) + (β − α)e^{iθ}, θ = kπΔx (10.26).
    Others (our derivations D07, D08, D18 — the book states the results only):
    "btcs" G = 1/(1 + 4β sin²(θ/2) + 2iα sin θ);  "cn" G = (1 − 2β sin²(θ/2) − iα sin θ)/(1 + 2β sin²(θ/2) + iα sin θ);
    "upwind" (10.29) with centred diffusion: G = 1 − C(1 − e^{−iθ}) − 4β sin²(θ/2) for C = 2α ≥ 0 and the mirror
    1 − C(e^{iθ} − 1) − 4β sin²(θ/2) for C < 0 ("upwind_printed" keeps the first form for any sign — slip R11);
    "lax_wendroff" / "maccormack": G = 1 − iC sin θ − C²(1 − cos θ) − 4β sin²(θ/2).
    Parameters: theta [rad] (float or array), alpha, beta (Eq. (10.11)), scheme.  Returns complex G (Python complex for a
    scalar θ).  Examples: FTCS (π/2, 0.1, 0.2) → 0.6 − 0.2i; BTCS (π, 0, 100) → 1/401.
    Validation: V1 |G|² of (10.24) equals (10.26); V3 a single-mode run grows by |G|ⁿ.  Label: analytic.
    """
    th = np.asarray(theta, dtype=float)
    s = scheme.lower()
    sh2 = np.sin(th / 2) ** 2
    if s == "ftcs":
        G = (alpha + beta) * np.exp(-1j * th) + (1 - 2 * beta) + (beta - alpha) * np.exp(1j * th)  # Eq. (10.24)
    elif s == "btcs":
        G = 1.0 / (1 + 4 * beta * sh2 + 2j * alpha * np.sin(th))
    elif s == "cn":
        G = (1 - 2 * beta * sh2 - 1j * alpha * np.sin(th)) / (1 + 2 * beta * sh2 + 1j * alpha * np.sin(th))
    elif s in ("upwind", "upwind_printed"):
        C = 2 * alpha
        if s == "upwind_printed" or C >= 0:
            G = 1 - C * (1 - np.exp(-1j * th)) - 4 * beta * sh2  # from Eq. (10.29)
        else:
            G = 1 - C * (np.exp(1j * th) - 1) - 4 * beta * sh2
    elif s in ("lax_wendroff", "maccormack"):
        C = 2 * alpha
        G = 1 - 1j * C * np.sin(th) - C ** 2 * (1 - np.cos(th)) - 4 * beta * sh2
    else:
        raise ValueError("scheme must be ftcs, btcs, cn, upwind, upwind_printed, lax_wendroff or maccormack")
    return complex(G) if np.ndim(G) == 0 else G


def amplification_modulus(theta, alpha: float = 0.0, beta: float = 0.0, scheme: str = "ftcs"):
    """|G(θ)| of :func:`amplification_factor` (the real-valued twin for parity rows).  Book: §10.2, Eq. (10.25).
    Label: analytic."""
    return as_scalar_if_0d(np.abs(amplification_factor(theta, alpha, beta, scheme)))


def amplification_curve(scheme: str = "ftcs", alpha: float = 0.0, beta: float = 0.0, n_theta: int = 361) -> dict:
    """G on θ ∈ [0, π] — the curve an explainer draws in the complex plane.  Book: §10.2, (10.24)–(10.26).
    Returns dict(theta [rad], G (complex array), modulus).  Label: analytic."""
    th = np.linspace(0.0, np.pi, int(n_theta))
    G = amplification_factor(th, alpha, beta, scheme)
    return dict(theta=th, G=G, modulus=np.abs(G))


def max_amplification(scheme: str = "ftcs", alpha: float = 0.0, beta: float = 0.0, n_theta: int = 2001) -> float:
    """max over θ ∈ [0, π] of |G(θ)| (Eq. (10.25) needs this ≤ 1); |G(−θ)| = |G(θ)| so [0, π] suffices.
    Example ("ftcs", 0, 0.51) → 1.04 (at θ = π).  Book: §10.2.  Label: analytic."""
    th = np.linspace(0.0, np.pi, int(n_theta))
    return float(np.max(np.abs(amplification_factor(th, alpha, beta, scheme))))


def worst_theta(scheme: str = "ftcs", alpha: float = 0.0, beta: float = 0.0, n_theta: int = 2001) -> float:
    """The θ ∈ [0, π] [rad] where |G| is largest (π = the two-cell zigzag, small θ = long waves).  Book: §10.2.
    Label: analytic."""
    th = np.linspace(0.0, np.pi, int(n_theta))
    return float(th[np.argmax(np.abs(amplification_factor(th, alpha, beta, scheme)))])


def is_von_neumann_stable(scheme: str = "ftcs", alpha: float = 0.0, beta: float = 0.0, n_theta: int = 2001,
                          tol: float = 1e-12) -> bool:
    """True if |G(θ)| ≤ 1 + tol for every sampled θ — Eq. (10.25).  Book: §10.2.  Label: analytic."""
    return max_amplification(scheme, alpha, beta, n_theta) <= 1.0 + tol


def ftcs_amplification_modulus2(theta, alpha: float, beta: float):
    """Closed form of |G|² for FTCS: (1 − 4β sin²(θ/2))² + (2α sin θ)².

    Book: §10.2, Eq. (10.26), θ = kπΔx (our D05: Euler's formula and 1 − cos θ = 2 sin²(θ/2)).  Example (π/2, 0.1, 0.2) →
    0.40.  Validation: V1 equals |(10.24)|² to 1e-12.  Label: analytic.
    """
    th = np.asarray(theta, dtype=float)
    return as_scalar_if_0d((1 - 4 * beta * np.sin(th / 2) ** 2) ** 2 + (2 * alpha * np.sin(th)) ** 2)  # Eq. (10.26)


def ftcs_stable(alpha: float, beta: float, tol: float = 1e-12) -> bool:
    """Noye's closed-form stability region of FTCS: 0 ≤ 4α² ≤ 2β ≤ 1.

    Book: §10.2, Eq. (10.27) (cited to Noye 1983; derived in our D06).  Pure diffusion (α = 0): β ≤ ½, Eq. (10.28); pure
    convection (β = 0, α ≠ 0): never.  Examples: (0.1, 0.2) True; (0.4, 0.2) False; (0, 0.51) False.
    Validation: V1 = brute-force max_θ|G| ≤ 1; V7 runs blow up just outside.  Label: analytic.
    """
    return bool((0.0 <= 4 * alpha ** 2 <= 2 * beta + tol) and (2 * beta <= 1.0 + tol))  # Eq. (10.27)


def stability_verdict(scheme: str = "ftcs", alpha: float = 0.0, beta: float = 0.0) -> dict:
    """Stable or not, and why — the explainers' status line.

    Book: §10.2, (10.25)–(10.30).  reason is one of the exact ASCII strings: "stable", "unstable: 4a^2 > 2b, long waves
    grow", "unstable: 2b > 1, the zigzag grows" (FTCS), "unstable: pure convection, every wave grows", "unstable: C > 1, the
    characteristic leaves the stencil", "unstable: |C| + 2b > 1, the zigzag grows" (upwind with diffusion: the limit is
    |C| + 2β ≤ 1), "unstable: C^2 + 2b > 1, the zigzag grows" (Lax–Wendroff / MacCormack with diffusion: C² + 2β ≤ 1),
    "stable for every dt (implicit)".  (C = 2a = uΔt/Δx, b = β = DΔt/Δx².)
    Returns dict(stable, reason, Gmax, theta_worst [rad]).  Label: analytic.
    """
    s = scheme.lower()
    gmax = max_amplification(s, alpha, beta)
    thw = worst_theta(s, alpha, beta)
    stable = gmax <= 1.0 + 1e-12
    if s in ("btcs", "cn"):
        reason = "stable for every dt (implicit)"
    elif stable:
        reason = "stable"
    elif s == "ftcs":
        if beta == 0 and alpha != 0:
            reason = "unstable: pure convection, every wave grows"
        elif 2 * beta > 1.0 + 1e-12:
            reason = "unstable: 2b > 1, the zigzag grows"
        else:
            reason = "unstable: 4a^2 > 2b, long waves grow"
    elif s in ("upwind", "upwind_printed", "lax_wendroff", "maccormack"):
        if abs(2 * alpha) > 1.0 + 1e-12 or (s == "upwind_printed" and alpha < 0):
            reason = "unstable: C > 1, the characteristic leaves the stencil"
        elif s in ("upwind", "upwind_printed"):
            # |G(π)| = |1 − 2|C| − 4β| ≤ 1  ⇔  |C| + 2β ≤ 1 (C = 2α)
            reason = "unstable: |C| + 2b > 1, the zigzag grows"
        else:
            # Lax–Wendroff / linear MacCormack with diffusion: |G(π)| = |1 − 2C² − 4β| ≤ 1  ⇔  C² + 2β ≤ 1
            reason = "unstable: C^2 + 2b > 1, the zigzag grows"
    else:
        reason = "unstable"
    return dict(stable=bool(stable), reason=reason, Gmax=gmax, theta_worst=thw)


def phase_error(scheme: str, C: float, theta, beta: float = 0.0):
    """Relative phase speed of a mode: (numerical phase shift per step)/(exact shift Cθ) = −arg G/(Cθ).

    Book: §10.2 (our companion to the amplitude |G| of (10.25); separates dispersion from dissipation).  C = uΔt/Δx,
    θ [rad] ∈ (0, π].  Examples: Lax–Wendroff (C = 0.5, θ = π/2) → 0.7487; upwind (C = 0.5) → 1.  Label: analytic.
    """
    th = np.asarray(theta, dtype=float)
    G = amplification_factor(th, 0.5 * C, beta, scheme)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = -np.angle(G) / (C * th)
    return as_scalar_if_0d(out)


def numerical_diffusivity(u, dx: float, D: float = 0.0, scheme: str = "upwind_steady", C: float | None = None):
    """Artificial (numerical) diffusivity added by a scheme's leading truncation term [m²/s].

    Book: §10.4, Eq. (10.94): steady first-order upwind solves uT_x = D(1 + 0.5R_cell)T_xx, i.e. D_num = 0.5R_cell D =
    |u|Δx/2 ("upwind_steady").  Unsteady upwind (ours, D08 note: modified equation with T_tt = u²T_xx): |u|Δx(1 − C)/2
    ("upwind"); FTCS: −u²Δt/2 = −|u|Δx C/2 (negative — anti-diffusion, why pure-convection FTCS fails); Lax–Wendroff,
    MacCormack and centred steady: 0 (their leading error is dispersive).
    Parameters: u [m/s], dx [m], D [m²/s] (unused, kept for the ratio D_num/D in callers), scheme (default
    "upwind_steady", the (10.94) value, which needs no C), C (Courant number uΔt/Δx — required for the unsteady schemes
    "upwind" and "ftcs"; a ValueError is raised when it is omitted, since C = 0 would silently return the steady value).
    Example (1, 0.01, C = 0.5, "upwind") → 0.0025.  Label: symbolic, converged.
    """
    uu = np.abs(np.asarray(u, dtype=float))
    if scheme in ("upwind", "ftcs") and C is None:
        raise ValueError(f'numerical_diffusivity(scheme="{scheme}") needs the Courant number C = u dt/dx; '
                         'use scheme="upwind_steady" for the steady (10.94) value |u| dx/2')
    if scheme == "upwind_steady":
        out = 0.5 * uu * dx  # Eq. (10.94): 0.5 R_cell D
    elif scheme == "upwind":
        out = 0.5 * uu * dx * (1.0 - C)
    elif scheme == "ftcs":
        out = -0.5 * uu * dx * C
    elif scheme in ("central", "central_steady", "lax_wendroff", "maccormack"):
        out = 0.0 * uu
    else:
        raise ValueError("scheme must be upwind_steady, upwind, ftcs, central, lax_wendroff or maccormack")
    return as_scalar_if_0d(out)


def ftcs2d_max_amplification(u: float, v: float, Re: float, dx: float, dt: float, n: int = 721,
                             dy: float | None = None) -> float:
    """max |G| of the explicit 2-D convection–diffusion step of the MAC predictor (forward Euler, centred convection).

    Book: §10.4 — the explicit step (10.116) with the MAC limits (10.127) ½(u² + v²)ΔtRe ≤ 1 and (10.128) 4Δt/(ReΔx²) ≤ 1
    (cited to Peyret & Taylor 1983; checked here by a von Neumann scan, ours):
    G = 1 − i(Cx sin θx + Cy sin θy) − 4dx sin²(θx/2) − 4dy sin²(θy/2), Cx = uΔt/Δx, dx = Δt/(ReΔx²).
    Parameters: u, v (non-dimensional), Re, dx, dt, n (samples per direction), dy (default dx).  Returns max |G|.
    Label: analytic.
    """
    dy = dx if dy is None else dy
    th = np.linspace(0.0, np.pi, int(n))
    TX, TY = np.meshgrid(th, th, indexing="xy")
    G = (1 - 1j * (u * dt / dx * np.sin(TX) + v * dt / dy * np.sin(TY)) - 4 * dt / (Re * dx ** 2) * np.sin(TX / 2) ** 2
         - 4 * dt / (Re * dy ** 2) * np.sin(TY / 2) ** 2)
    return float(np.max(np.abs(G)))


def _profile(name):
    if callable(name):
        return name
    if name == "square":
        return lambda s: (((s % 1.0) >= 0.2) & ((s % 1.0) < 0.4)).astype(float)
    if name in ("gauss", "gaussian"):
        return lambda s: np.exp(-((((s - 0.3) + 0.5) % 1.0 - 0.5) ** 2) / (2 * 0.05 ** 2))
    raise ValueError("profile must be square, gauss or a callable")


def advect_periodic(profile="square", C: float = 0.8, n_cells: int = 100, n_rev: float = 1.0, scheme: str = "upwind",
                    growth_cap: float = 1e6) -> dict:
    """Advect a profile round the periodic unit domain at Courant number C (pure convection T_t + uT_x = 0, u = 1).

    Book: §10.2, Eqs. (10.29)–(10.30) (upwind), (10.10) (FTCS: unstable for D = 0); §10.4, (10.101)–(10.102) MacCormack
    ("maccormack" = forward predictor/backward corrector, "maccormack_bf" the reverse; both equal Lax–Wendroff for linear
    advection).  Grid: x_i = i/N, i = 0 … N − 1.  Profiles: "square" (T = 1 for 0.2 ≤ x < 0.4) and "gauss"
    (e^{−(x−0.3)²/(2·0.05²)}, nearest periodic image).  One revolution = round(N/C) steps; C is adjusted to N·n_rev/steps
    so the pulse returns exactly.  Schemes: "upwind", "maccormack", "maccormack_bf", "lax_wendroff", "ftcs", "exact".
    Returns dict(x, T, T0, exact, amplitude_ratio (max T/max T0), rms_error, blew_up, steps, C (the adjusted value)).
    Validation: V1 C = 1 upwind and MacCormack are exact shifts; V3 orders 1 and 2.  Label: converged.
    """
    x = np.arange(n_cells) / n_cells
    f = _profile(profile)
    T0 = f(x)

    def shifted(k_steps: int) -> np.ndarray:
        # exact solution T0(x − k·C·Δx).  When k·C is a whole number of cells the answer is an integer node shift of T0
        # (np.roll, exact); otherwise evaluate the profile at the shifted argument rounded to 12 decimals so that
        # floating-point modulo (e.g. (0.2 − 1.0) % 1 = 0.19999999999999996) cannot move a square-pulse edge by a node.
        m = k_steps * Ce
        mi = int(round(m))
        if abs(m - mi) < 1e-9:
            return np.roll(T0, mi)
        return f(np.round((x - m / n_cells) % 1.0, 12))

    steps = max(1, int(round(n_rev * n_cells / abs(C))))
    Ce = np.sign(C) * n_rev * n_cells / steps
    T = T0.copy()
    scale = float(np.max(np.abs(T0)))
    blew, k = False, 0
    for k in range(1, steps + 1):
        if scheme == "maccormack":
            Ts = T - Ce * (np.roll(T, -1) - T)  # predictor, forward (10.101)
            T = 0.5 * (T + Ts - Ce * (Ts - np.roll(Ts, 1)))  # corrector, backward (10.102)
        elif scheme == "maccormack_bf":
            Ts = T - Ce * (T - np.roll(T, 1))
            T = 0.5 * (T + Ts - Ce * (np.roll(Ts, -1) - Ts))
        elif scheme == "exact":
            T = shifted(k)
        else:
            Tp, Tm = np.roll(T, -1), np.roll(T, 1)
            T = _explicit(T, Tp, Tm, 0.5 * Ce, 0.0, scheme)
        if np.max(np.abs(T)) > growth_cap * scale or not np.all(np.isfinite(T)):
            blew = True
            break
    ex = shifted(k)
    return dict(x=x, T=T, T0=T0, exact=ex, amplitude_ratio=float(np.max(T) / np.max(T0)) if not blew else float("inf"),
                rms_error=error_norm(T, ex) if not blew else float("inf"), blew_up=blew, steps=k, C=float(Ce))


def rod_heating_exact(x, t, D: float, T_wall: float, T_init: float = 0.0, L: float = 1.0, tol: float = 1e-12,
                      max_terms: int = 100000):
    """Exact temperature of a rod 0 ≤ x ≤ L at T_init whose two ends are raised to T_wall at t = 0 (Fourier sine series).

    Book: Exercise 10.3, Eq. (10.199) (the series form, here with our own wall and initial temperatures and length):
    T = T_w + (T_init − T_w) Σ_{m≥1} [4/((2m − 1)π)] sin[(2m − 1)πx/L] exp[−D(2m − 1)²π²t/L²].
    Terms are added until the term's envelope 4/((2m − 1)π)·exp[−D(2m − 1)²π²t/L²] (its bound over all x) is < tol, so
    the error of the sum is < tol·|T_w − T_init| times a factor ≤ 1 + O(e^{−8Dπ²t/L²}) — the envelope, not the actual
    term, because sin[(2m − 1)πx/L] vanishes for some m at rational x/L (x = L/5, L/3, …) and would stop the sum early;
    at t = 0 the series converges slowly (Gibbs) — the initial value is returned directly there.
    Parameters: x [m], t [s] (≥ 0), D [m²/s], T_wall, T_init (any temperature unit), L [m], tol, max_terms.
    Returns T (same unit).  Validation: V1 satisfies the PDE and BCs term by term (sympy); V3 FTCS/CN converge to it.
    Label: analytic.
    """
    x = np.asarray(x, dtype=float)
    dT = T_init - T_wall
    if t <= 0:
        out = np.where((x <= 0) | (x >= L), T_wall, T_init).astype(float)
        return as_scalar_if_0d(out)
    s = np.zeros_like(x)
    for m in range(1, max_terms + 1):
        k = (2 * m - 1) * np.pi / L
        env = 4.0 / ((2 * m - 1) * np.pi) * np.exp(-D * k ** 2 * t)  # envelope of the m-th term (|sin| ≤ 1)
        s = s + env * np.sin(k * x)  # Eq. (10.199)
        if env < tol:
            break
    return as_scalar_if_0d(T_wall + dT * s)


def lax_demo(betas=(0.45, 0.50, 0.51), n_cells: int = 20, nsteps: int = 2000, T_wall: float = 1.0, D: float = 1.0,
             L: float = 1.0, growth_cap: float = 1e6) -> dict:
    """The Lax equivalence theorem made visible: FTCS on the heated rod at several diffusion numbers β.

    Book: §10.2 — Lax equivalence (consistency + stability ⇔ convergence, for a well-posed linear problem), FTCS pure
    diffusion stable iff β ≤ ½ (10.28), the rod of Exercise 10.3 with our numbers (walls at T_wall, rod initially 0, exact
    (10.199) from :func:`rod_heating_exact`).  Δt = βΔx²/D; both ends Dirichlet.
    Returns dict β → dict(max_error_history (max |T − T_exact| per step), blew_up, step_blown, final_error (rms at the last
    step), t_final).  Expect β = 0.45, 0.50 bounded; 0.51 blown (the zigzag θ = π grows by |G| = |1 − 4β| = 1.04 per step;
    with the defaults the growth cap 1e6 is reached at step 690).
    Validation: the blow-up step and the bounded errors are asserted; no order is measured here.
    Label: qualitative (demonstration).
    """
    dx = L / n_cells
    x = np.linspace(0.0, L, n_cells + 1)
    out = {}
    for b in betas:
        dt = b * dx ** 2 / D
        T = np.zeros(n_cells + 1)
        T[0] = T[-1] = T_wall
        errs = []
        blew, sb, k = False, None, 0
        scale = abs(T_wall)
        for k in range(1, nsteps + 1):
            T = transport_1d_step(T, 0.0, b, "ftcs", periodic=False, g=T_wall, T_L=T_wall)
            e = np.abs(T - rod_heating_exact(x, k * dt, D, T_wall, 0.0, L))
            errs.append(float(np.max(e)))
            if np.max(np.abs(T)) > growth_cap * scale or not np.all(np.isfinite(T)):
                blew, sb = True, k
                break
        ex = rod_heating_exact(x, k * dt, D, T_wall, 0.0, L)
        out[float(b)] = dict(max_error_history=np.array(errs), blew_up=blew, step_blown=sb,
                             final_error=error_norm(T, ex) if not blew else float("inf"), t_final=k * dt)
    return out


# ======================================================================================================================
# Steady convection–diffusion: the layer, centred wiggles and upwind numerical diffusion — (10.84)–(10.94)
# ======================================================================================================================
def steady_cd_exact(x, R: float, L: float = 1.0, limit: str | None = None, printed: bool = False):
    """Exact solution of the steady convection–diffusion problem u T_x = D T_xx, T(0) = 0, T(L) = 1.

    Book: §10.4, Eqs. (10.84)–(10.87): T = (e^{Rx/L} − 1)/(e^R − 1), R = uL/D; Eq. (10.88) large-R form T = e^{−R(1−x/L)}
    (``limit="large_R"``).
    Overflow: the printed form overflows for R ≳ 709 (e^R > 1.8e308).  We evaluate the identical scaled form
    T = e^{−R(1−ξ)}(1 − e^{−Rξ})/(1 − e^{−R}) for R > 0 (``expm1``), expm1(Rξ)/expm1(R) for R < 0, ξ for R = 0;
    ``printed=True`` evaluates (10.86) literally (inf/NaN for large R — the overflow demonstration).
    Parameters: x [m] (0 ≤ x ≤ L), R (global Péclet number, any sign), L [m].  Returns T (0 … 1).
    Example R = 4 at x = 0.25, 0.5, 0.75: 0.03206, 0.1192, 0.3561.
    Validation: V2 satisfies (10.84)–(10.85) (sympy); V1 scaled = printed to 1e-14 for |R| ≤ 50, finite for R = 1e4.
    Label: analytic, symbolic.
    """
    xi = np.asarray(x, dtype=float) / L
    if limit == "large_R":
        return as_scalar_if_0d(np.exp(-R * (1.0 - xi)))  # Eq. (10.88)
    if printed:
        with np.errstate(over="ignore", invalid="ignore"):
            return as_scalar_if_0d((np.exp(R * xi) - 1.0) / (np.exp(R) - 1.0))  # Eq. (10.86) as printed
    if R == 0:
        return as_scalar_if_0d(xi.copy())
    if R > 0:
        out = np.exp(-R * (1.0 - xi)) * np.expm1(-R * xi) / np.expm1(-R)  # Eq. (10.86), scaled (no overflow)
    else:
        out = np.expm1(R * xi) / np.expm1(R)
    return as_scalar_if_0d(out)


def cd_layer_thickness(R: float, L: float = 1.0, level: float = float(np.exp(-1.0))) -> float:
    """Distance from the wall x = L at which the exact solution (10.86) has fallen to ``level`` (default e⁻¹ ≈ 0.368).

    Book: §10.4, Eq. (10.89) δ/L = O(1/|R|) and the levels quoted below it (computed here, not quoted): for R ≫ 1 the
    distance → L ln(1/level)/R.  Overflow-safe closed form (ours): 1 − ξ = −ln(level + (1 − level)e^{−R})/R, R > 0.
    Parameters: R (> 0), L [m], level ∈ (0, 1).  Returns δ [m].  Example R = 100 → 0.01.  Label: analytic.
    """
    require_positive("R", R)
    return float(-L * np.log(level + (1.0 - level) * np.exp(-R)) / R)


def stretched_grid(n: int, L: float = 1.0, beta_s: float = 2.0) -> np.ndarray:
    """Nodes x_j = L tanh(β_s j/n)/tanh β_s, j = 0 … n — clustered toward the wall x = L (the slope dx/dj is smallest there).

    Book: §10.4 remark after (10.92) ("nonuniform grids with local fine grid spacing inside the boundary layer"); the map is
    our choice.  β_s → 0 gives the uniform grid.
    # DEVIATION: design C.10.3 wrote x_j = L[1 − tanh(β_s(1 − j/n))/tanh β_s], which clusters the nodes at x = 0 (its slope is
    smallest at j = 0), away from the layer; we use the map that refines at x = L.
    It refines only at x = L, where the layer of (10.86)–(10.88) sits for R > 0.  For R < 0 the layer is at x = 0 and this
    grid puts its coarse cells there (worse than uniform); mirror it, x → L − x[::-1], for that case — neither
    :func:`steady_cd_fd` nor this function does so automatically.  Label: analytic.
    """
    j = np.arange(n + 1) / n
    if beta_s <= 1e-12:
        return L * j
    return L * np.tanh(beta_s * j) / np.tanh(beta_s)


def steady_cd_fd(n: int, R: float, scheme: str = "central", grid: str = "uniform", L: float = 1.0, beta_s: float = 2.0) -> dict:
    """Finite-difference solution of the steady problem (10.84)–(10.85) on n cells.

    Book: §10.4, Eqs. (10.90)–(10.91) centred: 0.5R_cell(T_{j+1} − T_{j−1}) = T_{j+1} − 2T_j + T_{j−1}, R_cell = R/n;
    Eq. (10.93) first-order upwind: R_cell(T_j − T_{j−1}) = T_{j+1} − 2T_j + T_{j−1} (slip R3: the book calls this
    "forward"; it is the backward = upwind difference for u > 0; for R < 0 the upwind side flips); "forward" (ours, R3
    illustration): R_cell(T_{j+1} − T_j) = … — the downwind difference, root r = 1/(1 − R_cell): monotone (and
    anti-diffusive, effective D(1 − 0.5R_cell)) for R_cell < 1, wiggles only for R_cell > 1.
    ``grid="stretched"``: :func:`stretched_grid` nodes (refined at x = L only — meant for R > 0) with the second-order nonuniform 3-point stencils
    (centred) or the one-sided differences (upwind / forward).  The stretched grid *shrinks* the centred wiggles (R = 80,
    n = 20: min T from −0.35 to below 1e-4 in magnitude, visibly gone) but does not remove them: the sign alternation stays
    in the coarse upstream cells, where the local R_cell ≈ 8 (:func:`wiggle_indicator` still reports it).
    Parameters: n (cells), R (uL/D), scheme, grid, L [m], beta_s.
    Returns dict(x [m], T) with T_0 = 0, T_n = 1 (tridiagonal ``scipy.linalg.solve_banded``).  Examples n = 4: R = 4 central
    → [0, 0.025, 0.1, 0.325, 1]; R = 16 central → [0, −0.05, 0.1, −0.35, 1]; R = 16 upwind → [0, 0.00641, 0.03846, 0.1987, 1].
    Validation: V1 equals :func:`steady_cd_discrete_exact` (uniform) to 1e-13; V3 orders 2 (centred) and 1 (upwind).
    Label: analytic, converged.
    """
    from scipy.linalg import solve_banded

    x = L * np.arange(n + 1) / n if grid == "uniform" else stretched_grid(n, L, beta_s)
    hm = np.diff(x)[:-1]
    hp = np.diff(x)[1:]
    a = R / L  # u/D [1/m]:  u T_x = D T_xx  ⇔  a T_x − T_xx = 0
    d2m, d2c, d2p = 2.0 / (hm * (hm + hp)), -2.0 / (hm * hp), 2.0 / (hp * (hm + hp))
    zero = 0.0 * hm
    if scheme == "central":
        den = hp * hm * (hp + hm)
        d1m, d1c, d1p = -hp ** 2 / den, (hp ** 2 - hm ** 2) / den, hm ** 2 / den  # uniform: (T_{j+1} − T_{j−1})/(2Δx), Eq. (10.90)
    elif scheme == "upwind":
        if R >= 0:
            d1m, d1c, d1p = -1.0 / hm, 1.0 / hm, zero  # Eq. (10.93): (T_j − T_{j−1})/Δx
        else:
            d1m, d1c, d1p = zero, -1.0 / hp, 1.0 / hp
    elif scheme == "forward":
        d1m, d1c, d1p = zero, -1.0 / hp, 1.0 / hp  # downwind for u > 0 (R3 illustration)
    else:
        raise ValueError("scheme must be central, upwind or forward")
    lo, di, up = a * d1m - d2m, a * d1c - d2c, a * d1p - d2p
    m = n - 1
    ab = np.zeros((3, m))
    ab[0, 1:] = up[:-1]
    ab[1, :] = di
    ab[2, :-1] = lo[1:]
    rhs = np.zeros(m)
    rhs[-1] = -up[-1]  # T_n = 1 moved to the right-hand side
    T = np.empty(n + 1)
    T[0], T[-1] = 0.0, 1.0
    T[1:-1] = solve_banded((1, 1), ab, rhs)
    return dict(x=x, T=T)


def discrete_root(R_cell: float, scheme: str = "central") -> float:
    """The nontrivial root r of the geometric trial T_j = r^j in the discrete equations (our D16).

    Book: §10.4, Eqs. (10.91), (10.93): centred r = (1 + R_cell/2)/(1 − R_cell/2) (negative ⇔ R_cell > 2: alternating signs,
    the wiggles; infinite at R_cell = 2); upwind r = 1 + R_cell (> 0 always).  Examples: (1) → 3, (4) → −3,
    ("upwind", 4) → 5.  Label: analytic.
    """
    if scheme == "central":
        if abs(R_cell - 2.0) < 1e-15:
            return float("inf")
        return float((1.0 + 0.5 * R_cell) / (1.0 - 0.5 * R_cell))
    if scheme == "upwind":
        return float(1.0 + R_cell) if R_cell >= 0 else float(1.0 / (1.0 - R_cell))
    raise ValueError("scheme must be central or upwind")


def steady_cd_discrete_exact(n: int, R: float, scheme: str = "central", L: float = 1.0) -> dict:
    """Closed-form solution of the discrete equations (10.91) / (10.93): T_j = (r^j − 1)/(r^n − 1), r = :func:`discrete_root`.

    Book: §10.4 (the book argues from a heat balance at j = n − 1; the closed form is our D16).  R_cell = 2 (centred) is
    degenerate: T_j = 0 for j < n (the limit r → ∞).  Overflow-safe for |r| > 1 (r^{j−n}(1 − r^{−j})/(1 − r^{−n})).
    Returns dict(x [m], T, r).  Validation: V1 equals the tridiagonal solve of :func:`steady_cd_fd` to 1e-13.
    Label: analytic.
    """
    Rc = R / n
    j = np.arange(n + 1, dtype=float)
    x = L * j / n
    r = discrete_root(Rc, scheme)
    if not np.isfinite(r):
        T = np.zeros(n + 1)
        T[-1] = 1.0
    elif abs(r - 1.0) < 1e-15:
        T = j / n
    elif abs(r) > 1.0:
        ri = 1.0 / r
        T = ri ** (n - j) * (1.0 - ri ** j) / (1.0 - ri ** n)
    else:
        T = (r ** j - 1.0) / (r ** n - 1.0)
    return dict(x=x, T=T, r=float(r))


def wiggle_indicator(T, tol: float = 1e-9) -> dict:
    """Detect node-to-node oscillations in a profile that should be monotone.

    Book: §10.4 (below (10.92): "an unphysical oscillatory solution with T_j < 0 … may be a virtue since it provides a
    warning").  Differences smaller than tol × (max T − min T) are ignored (round-off level).
    Returns dict(min (smallest value), n_sign_changes (of the node-to-node differences), wiggles (bool), has_negative).
    Label: analytic.
    """
    T = np.asarray(T, dtype=float)
    d = np.diff(T)
    d = d[np.abs(d) > tol * max(float(T.max() - T.min()), 1e-300)]
    sc = int(np.sum(np.sign(d[1:]) != np.sign(d[:-1]))) if d.size > 1 else 0
    return dict(min=float(T.min()), n_sign_changes=sc, wiggles=bool(sc > 0), has_negative=bool(T.min() < -tol))
