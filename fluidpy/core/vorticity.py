"""Vortex kinematics and vorticity dynamics: vortex lines and tubes, material circulation and Kelvin's theorem, the
baroclinic torque, Helmholtz's frozen-in lines, and the vorticity equation split into its terms (numeric and sympy).

Book: Ch. 5 §5.1–5.6, Eqs. (5.3)–(5.4), (5.8)–(5.13), (5.18)–(5.33), Figs. 5.1, 5.4, 5.6, 5.7, 5.9, 5.10 (rendered
pages chapters/pages/ch05/p199–p214). Every equation was transcribed from the page images.

Conventions (ch03/ch04, unchanged; design Part C): a field is a callable ``F(x, t)`` with ``x`` of shape ``(d,)`` or
``(d, N)`` (components on axis 0, points last); tensors G[i, j] = ∂u_i/∂x_j; stencils are explicit second-order central
differences (``core._stencil``) with a length step ``h`` [m] and a time step ``ht`` [s] (never the same default).
**ω is the vorticity** everywhere in this module (twice the local angular velocity); plane vorticity and circulation
are counterclockwise positive about +z. SI units throughout; no non-dimensionalisation (the chapter has none).

Closed loops are ``(d, N)`` point arrays, **counterclockwise and not closed** (the last point is not repeated),
sampled uniformly in a periodic label s_k = k/N (the material labels of the particles that form the loop, D04). Line
integrals ∮ f·dx = Σ_k f(x_k)·x′(s_k)/N use the spectral (FFT) derivative x′(s) and the periodic trapezoid rule,
which converge exponentially for smooth loops (polygons: use :func:`square_loop_points`, whose corners are smoothed in
the parametrisation).

Reuse: Ch. 6 (irrotational flow stays irrotational), Ch. 8 (vortex decay, Burgers), Ch. 9 (wall vorticity), Ch. 11
(baroclinic and shear-layer instability), Ch. 12 (stretching, enstrophy), **Ch. 13 (Bjerknes/Kelvin in a rotating
frame, baroclinic generation, planetary stretching → potential vorticity)**, Ch. 14 (trailing vortices).
"""
from __future__ import annotations

import inspect
from dataclasses import dataclass, field
from typing import Callable, NamedTuple

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

from . import _stencil as st
from ._util import as_scalar_if_0d
from .integral_theorems import Loop, Surface, _plane_basis, circulation, curl_flux, gauss_legendre_nodes
from .kinematics import acceleration, as_coord_field, streamline, velocity_gradient_at, vorticity_from_gradient

__all__ = ["planar_field_3d", "stencil_gradient_fn", "circle_loop_points", "square_loop_points", "loop_tangent",
           "loop_line_integral", "loop_circulation", "loop_vector_area", "loop_length",
           "vorticity_field", "vortex_line", "TubeStrength", "vortex_tube_strength", "Tube", "vortex_tube",
           "TubeFlux", "tube_flux_budget", "tube_flux_budget_traced",
           "material_loop", "material_circulation", "KelvinRate", "kelvin_rate_terms", "KelvinForces",
           "kelvin_force_terms", "pressure_torque_on_element", "frozen_in_check",
           "vorticity_equation_sym", "VorticityTerms", "vorticity_terms", "vorticity_terms_sym", "VorticityBudget",
           "vorticity_budget", "vorticity_budget_sym", "VORTICITY_BUDGET_PRESETS", "vorticity_budget_preset",
           "vorticity_divergence", "rotating_lamb_form_terms", "baroclinic_term", "baroclinic_term_sym",
           "stretching_tilting_split", "planetary_vorticity_terms", "absolute_circulation",
           "KELVIN_SCENARIOS", "kelvin_scenario", "kelvin_scenario_circulation", "kelvin_scenario_rate"]

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d


# ======================================================================================================================
# Small helpers: 3-D embedding, presets, periodic loops
# ======================================================================================================================
def planar_field_3d(u2: Callable) -> Callable:
    """Embed a plane field u(x, t) → (2, …) into 3-D: u₃(x, t) = (u, v, 0) evaluated at (x, y) (z-independent).

    Lets the 3-D vorticity-equation tools act on the plane fields of Ch. 3–4 (Lamb–Oseen, Taylor–Green, cellular flow).
    Label: analytic (re-packaging). Book: §5.4 (plane flow: ω = ω_z e_z, no stretching or tilting).
    """
    def u(x, t=0.0):
        x_ = _F(x)
        return st.pad3(_F(u2(x_[:2], t)))
    return u


def stencil_gradient_fn(f: Callable, h: float = 1e-4) -> Callable:
    """∇f(x, t) of a scalar field by 2nd-order central differences (a callable with the field convention).
    Label: converged (order 2). Book: helper for (5.10) and (5.28)."""
    return lambda x, t=0.0: st.grad(f, x, t, h)


def _omega_callable(omega, **p) -> Callable:
    """A vorticity callable, or a preset name of ``core.vortices`` ("gaussian_tube", "broken", "abc", "burgers") →
    the vorticity callable built with the keywords ``p``."""
    if callable(omega):
        return omega
    from . import vortices as VX
    if omega == "gaussian_tube":
        return VX.gaussian_tube_field(**p)
    if omega == "broken":
        return VX.broken_tube_field(**p)
    if omega == "abc":
        return VX.abc_flow_field(**p)  # Beltrami: ω = u
    if omega == "burgers":
        q = dict(p)
        h = q.pop("h", 1e-5)
        return vorticity_field(VX.burgers_vortex_field(**q), h)
    raise ValueError(f"unknown vorticity preset {omega!r}")


def circle_loop_points(center, radius: float, n: int = 256, normal=None) -> np.ndarray:
    """A circle sampled uniformly in θ ∈ [0, 2π) (θ_k = 2πk/n), counterclockwise about ``normal`` (default +z for a
    2-vector centre). Returns (d, n) points [m]. Label: analytic. Book: helper (the loop C of (5.8)–(5.11))."""
    c = _F(center)
    th = 2.0 * np.pi * np.arange(n) / n
    if c.size == 2 and normal is None:
        return np.stack([c[0] + radius * np.cos(th), c[1] + radius * np.sin(th)])
    nvec = np.array([0.0, 0.0, 1.0]) if normal is None else _F(normal)
    c3 = c if c.size == 3 else np.array([c[0], c[1], 0.0])
    e1, e2, _ = _plane_basis(nvec)
    return c3[:, None] + radius * (np.cos(th)[None, :] * e1[:, None] + np.sin(th)[None, :] * e2[:, None])


def square_loop_points(center, side: float, n: int = 512) -> np.ndarray:
    """A square of side ``side`` about ``center`` (plane), counterclockwise, n points (a multiple of 4). Each side is
    parametrised by τ(u) = u − sin(2πu)/(2π), u ∈ [0, 1) (dx/ds and d²x/ds² vanish at the corners), so the periodic
    loop is smooth in its label and the spectral quadrature of :func:`loop_line_integral` stays accurate.
    Returns (2, n) [m]. Label: analytic. Book: helper (the square loop of E3's baroclinic scenario)."""
    c = _F(center)
    m = n // 4
    u = np.arange(m) / m
    tau = u - np.sin(2.0 * np.pi * u) / (2.0 * np.pi)
    h = 0.5 * side
    corners = [(-h, -h), (h, -h), (h, h), (-h, h)]
    pts = []
    for k in range(4):
        a, b = _F(corners[k]), _F(corners[(k + 1) % 4])
        pts.append(a[:, None] + (b - a)[:, None] * tau[None, :])
    return np.concatenate(pts, axis=1) + c[:, None]


def loop_tangent(pts) -> np.ndarray:
    """dx/ds of a closed loop sampled uniformly in the label s ∈ [0, 1) (s_k = k/N): spectral derivative (FFT),
    shape (d, N). Book: helper for the loop integrals of (5.9)–(5.11) and (5.33). Label: analytic (exact for
    trigonometric polynomials; exponentially convergent for smooth loops)."""
    P_ = _F(pts)
    N = P_.shape[-1]
    k = 2.0 * np.pi * np.fft.fftfreq(N, d=1.0 / N)
    if N % 2 == 0:
        k[N // 2] = 0.0  # the Nyquist mode has no well-defined derivative on a real grid
    return np.real(np.fft.ifft(1j * k * np.fft.fft(P_, axis=-1), axis=-1))


def loop_line_integral(F_vals, pts) -> float:
    """∮ F·dx = ∫₀¹ F(x(s))·x′(s) ds ≈ Σ_k F_k·x′(s_k)/N (periodic trapezoid rule); F_vals (d, N) at the loop points.
    Label: analytic, converged (spectral). Book: the line integrals of (5.9)–(5.11)."""
    F_ = _F(F_vals)
    T = loop_tangent(pts)
    d = min(F_.shape[0], T.shape[0])
    return float(np.sum(F_[:d] * T[:d]) / T.shape[-1])


def loop_circulation(u, pts, t: float = 0.0) -> float:
    """Circulation Γ = ∮_C u·dx of the loop pts (d, N) at time t [m²/s], Eq. (3.18) (the Γ of (5.8)).
    Label: analytic, converged. Book: §5.2, Eq. (5.8) (definition (3.18))."""
    P_ = _F(pts)
    return loop_line_integral(_F(u(P_, t)), P_)


def loop_vector_area(pts) -> np.ndarray:
    """Vector area A_vec = ½∮ x × dx of a closed loop (3,) [m²]; a counterclockwise loop in the (x, y) plane gives
    (0, 0, A). Book: §5.6, Eq. (5.33) (∫_A Ω·n dA = Ω·A_vec for uniform Ω; our D19).
    Validation: V1 planar circle: (0, 0, πr²). Label: analytic."""
    P3 = st.pad3(_F(pts))
    T3 = st.pad3(loop_tangent(_F(pts)))
    return 0.5 * np.sum(np.cross(P3.T, T3.T), axis=0) / P3.shape[-1]  # ½ ∮ x × dx


def loop_length(pts) -> float:
    """Length ∮|dx| of a closed loop [m] (spectral tangent, trapezoid rule) — how much a material loop has been
    stretched (§5.2, Fig. 5.4). Book: helper for (5.8). Label: analytic."""
    T = loop_tangent(pts)
    return float(np.sum(np.linalg.norm(T, axis=0)) / T.shape[-1])


# ======================================================================================================================
# §5.1 vortex lines and tubes, (5.3)–(5.4)
# ======================================================================================================================
def vorticity_field(u, h: float = 1e-4) -> Callable:
    """The vorticity ω = ∇ × u of a velocity field as a callable ω(x, t) [1/s] (2nd-order stencils; the stencil
    version of ``core.kinematics.vorticity``). 3-D points give (3, …); plane (2-D) points give ω_z with shape (1, …).

    Book: §5.1 (ω = ∇ × u, twice the particle spin; (3.15)–(3.16)). Parameters: u u(x, t) [m/s]; h step [m].
    Validation: V1 u = b × x → 2b; V3 order 2. Label: analytic, converged.
    """
    def omega(x, t=0.0):
        x_ = _F(x)
        w = vorticity_from_gradient(st.grad_vector(u, x_, t, h))  # ω = ∇ × u via G = ∇u
        return w[2:3] if x_.shape[0] == 2 else w
    return omega


def vortex_line(omega, x0, t: float = 0.0, s_max: float = 1.0, n: int = 400, both: bool = True,
                rtol: float = 1e-10, atol: float = 1e-12, **p):
    """A vortex line through x0: the curve tangent to ω everywhere, Eq. (5.3) dx/ω_x = dy/ω_y = dz/ω_z.

    Book: §5.1, Eq. (5.3) — the streamline equation (3.7) with ω in place of u. Integrated as the arc-length ODE
    dx/ds = ω/|ω| (``core.kinematics.streamline`` applied to the vorticity field, DOP853); the line stops where
    |ω| → 0 (no vortex lines in irrotational fluid).

    Parameters
    ----------
    omega : ω(x, t) [1/s] callable, or a preset name ("gaussian_tube", "broken", "abc", "burgers") with keywords ``p``
    x0 : seed point [m] (3,);  t : frozen time [s];  s_max : arc length [m] each way;  n : samples;  both : both ways

    Returns
    -------
    (s, X) — arc length s [m] (0 at the point nearest x0) and points X (3, m).

    Validation: V1 ABC flow (ω = u): the vortex line equals the streamline; helical swirl u_φ = aRz: lines keep zR²
    constant; V3 tolerance study. Label: analytic, converged.
    """
    W = _omega_callable(omega, **p)
    X = streamline(W, x0, t, s_max, both, n, rtol, atol)
    ds = np.linalg.norm(np.diff(X, axis=1), axis=0)
    s = np.concatenate([[0.0], np.cumsum(ds)])
    i0 = int(np.argmin(np.linalg.norm(X - _F(x0)[:, None], axis=0)))
    return s - s[i0], X


class TubeStrength(NamedTuple):
    """Strength of a vortex tube by the two routes of §5.1: ``circulation`` ∮u·dx and ``flux`` ∫ω·n dA [m²/s]."""

    circulation: float
    flux: float
    diff: float


def vortex_tube_strength(u=None, loop: Loop | None = None, surface: Surface | None = None, omega=None,
                         t: float = 0.0, fd_step: float = 1e-4) -> TubeStrength:
    """Strength of a vortex tube: the circulation on a circuit round the tube and the vorticity flux through a
    cross-section, Γ = ∮_C u·dx = ∫_A ω·n dA (Stokes, (2.34)); dΓ = ω·n dA ↔ dQ = u·n dA.

    Book: §5.1 (the strength of a vortex tube is the circulation on a circuit lying on the tube and encircling it once;
    by Stokes' theorem it equals the vorticity integrated over the cross-section; Fig. 5.1). Orientation: the loop runs
    counterclockwise about the surface normal n (n_c into A, ch02).

    Parameters
    ----------
    u : velocity u(x, t) [m/s] (circulation route; also the flux route if ``omega`` is None)
    loop : ``core.integral_theorems.Loop`` (``planar_loop``);  surface : the matching ``Surface`` (``planar_disc``)
    omega : exact vorticity ω(x, t) [1/s] for the flux route (plane fields may return (1, …) or (3, …))
    t : time [s];  fd_step : step of the 4th-order differences when ``omega`` is None

    Returns
    -------
    TubeStrength(circulation, flux, diff) [m²/s]; a missing route is NaN.

    Validation: V1 Lamb–Oseen circle r = σ: both routes = Γ(1 − e⁻¹) = 0.632121Γ; reversing the normal flips both.
    Label: analytic.
    """
    circ = flux = np.nan
    if u is not None and loop is not None:
        circ = circulation(as_coord_field(u, t), loop)  # Γ = ∮ u·dx
    if surface is not None and (omega is not None or u is not None):
        if omega is not None:
            wf = as_coord_field(omega, t)
            if surface.dim == 2:
                curl_fn = lambda *c: _F(wf(*c))[-1]  # noqa: E731  plane: (∇×u)₃ = ω_z
            else:
                curl_fn = wf
            flux = curl_flux(as_coord_field(u, t) if u is not None else wf, surface, curl_fn=curl_fn)
        else:
            flux = curl_flux(as_coord_field(u, t), surface, fd_step=fd_step)  # ∫ ω·n dA
    return TubeStrength(float(circ), float(flux), float(circ - flux))


@dataclass(frozen=True)
class Tube:
    """A vortex tube traced from a seed circle: ``lines`` — list of (3, n) vortex lines sampled at the same arc lengths
    ``s`` (n,); ``lines[j][:, k]`` over j is the closed cross-section curve number k (counterclockwise about ω)."""

    lines: list
    seed_center: np.ndarray
    seed_normal: np.ndarray
    seed_radius: float
    s: np.ndarray = field(default_factory=lambda: np.zeros(0))

    @property
    def array(self) -> np.ndarray:
        """The lines stacked as (L, 3, n)."""
        return np.stack(self.lines)


def vortex_tube(omega, seed_center, seed_normal, seed_radius: float, n_lines: int = 24, s_max: float = 1.0,
                n: int = 200, t: float = 0.0, both: bool = True, rtol: float = 1e-10, atol: float = 1e-12,
                **p) -> Tube:
    """Trace a vortex tube: the vortex lines (5.3) through ``n_lines`` points of a seed circle (Fig. 5.1).

    Book: §5.1, Fig. 5.1 (the vortex lines through a closed curve form a tubular surface, the vortex tube). Each line
    is integrated along ω (and against it when ``both``) for arc length ``s_max`` (``core.kinematics.streamline``);
    cross-section k is the closed curve through the k-th sample of every line.

    Parameters: omega ω(x, t) [1/s] or a preset name (see :func:`vortex_line`) with keywords ``p``; seed_center [m]
    (3,); seed_normal (3,) (ω should point along it); seed_radius [m]; n_lines; s_max [m]; n samples per line; t [s];
    both. Returns :class:`Tube`.
    Raises ValueError if a line stops early (|ω| → 0 — the seed circle reaches irrotational fluid).
    Label: converged.
    """
    W = _omega_callable(omega, **p)
    c = _F(seed_center)
    e1, e2, nn = _plane_basis(seed_normal)
    th = 2.0 * np.pi * np.arange(n_lines) / n_lines
    lines = []
    for a in th:
        x0 = c + seed_radius * (np.cos(a) * e1 + np.sin(a) * e2)
        X = streamline(W, x0, t, s_max, both, n, rtol, atol)
        if X.shape[1] != n:
            raise ValueError("a vortex line stopped early (|ω| ≈ 0): choose a smaller seed radius or s_max")
        lines.append(X)
    s = np.linspace(-s_max, s_max, n) if both else np.linspace(0.0, s_max, n)
    return Tube(lines, c, nn, float(seed_radius), s)


class TubeFlux(NamedTuple):
    """Gauss budget (5.4) on a piece of vortex tube, all **outward** fluxes ∫ω·n dA [m²/s]: ``lower`` (= −Γ_lower end),
    ``side`` (the tube wall, ≈ 0), ``upper`` (= +Γ_upper end), ``total`` = lower + side + upper = ∫_V ∇·ω dV."""

    lower: float
    side: float
    upper: float
    total: float


def tube_flux_budget(field="gaussian_tube", z_a: float = 0.0, z_b: float = 1.0, R0: float = 0.1, n: int = 48,
                     R_of_z: Callable | None = None, t: float = 0.0, **p) -> TubeFlux:
    """Gauss' theorem on a piece of an **axisymmetric** vortex tube, Eq. (5.4):
    ∫_V ∇·ω dV = ∫_A ω·n dA = −Γ_lower end + Γ_upper end = 0.

    Book: §5.1, Eq. (5.4), Fig. 5.1: the curved side carries no flux (ω is tangent to it), so the strength is the same
    at both ends — tubes cannot end inside the fluid (kinematic: ∇·ω = ∇·(∇ × u) = 0). The tube wall is the surface of
    revolution R = R_t(z) through (R0, z = 0). ``lower`` = −2π∫₀^{R_t(z_a)} ω_z R dR (outward normal −e_z),
    ``upper`` = +2π∫₀^{R_t(z_b)} ω_z R dR, ``side`` = 2π∫ R_t[ω_R − R_t′ω_z] dz over the wall (n outward), all by
    Gauss–Legendre with ``n`` nodes; ``total`` = the sum.

    Parameters
    ----------
    field : "gaussian_tube" (``core.vortices.gaussian_tube_field``: R_t = R0 e^{−z/2L}), "lamb_oseen" (straight tube,
        keywords Gamma, sigma; R_t = R0), "broken" (``broken_tube_field``, R_t = R0 — **not** solenoidal) or a
        callable ω(x, t) (Cartesian, axisymmetric about the z-axis; then ``R_of_z`` is required)
    z_a, z_b : end heights [m];  R0 : wall radius at z = 0 [m];  n : Gauss–Legendre nodes per direction;  t [s]
    R_of_z : wall radius R_t(z) [m] for a callable field;  p : keywords of the preset (Gamma, a0, L, twist, sigma)

    Returns
    -------
    :class:`TubeFlux` [m²/s].

    Validation: V4 "gaussian_tube" (R0 = a0 = 0.1, Γ = 1, L = 1): lower −0.632121, upper +0.632121, side ≈ 0,
    total ≤ 1e-8; "broken": upper 0.232544 = 0.632121·e⁻¹, total −0.399577 (fails, as it must). Label: conserved.
    """
    from . import vortices as VX

    if callable(field):
        W = field
        if R_of_z is None:
            raise ValueError("a callable field needs R_of_z")
        Rt = R_of_z
    elif field == "gaussian_tube":
        L = float(p.get("L", 1.0))
        W = VX.gaussian_tube_field(**p)
        Rt = lambda z: R0 * np.exp(-_F(z) / (2.0 * L))  # noqa: E731  R/a(z) = const on the tube wall
    elif field == "broken":
        W = VX.broken_tube_field(**p)
        Rt = lambda z: R0 + 0.0 * _F(z)  # noqa: E731
    elif field == "lamb_oseen":
        Gam, sig = float(p.get("Gamma", 1.0)), float(p.get("sigma", 0.1))

        def W(x, tt=0.0):
            X = _F(x)
            r2 = X[0] ** 2 + X[1] ** 2
            return np.stack([0.0 * r2, 0.0 * r2, Gam / (np.pi * sig ** 2) * np.exp(-r2 / sig ** 2)])
        Rt = lambda z: R0 + 0.0 * _F(z)  # noqa: E731
    else:
        raise ValueError('field must be "gaussian_tube", "lamb_oseen", "broken" or a callable')

    def cap(z):
        R1 = float(Rt(z))
        r, wr = gauss_legendre_nodes(0.0, R1, n)
        pts = np.stack([r, 0.0 * r, np.full_like(r, z)])
        wz = _F(W(pts, t))[2]
        return 2.0 * np.pi * float(np.sum(wz * r * wr))  # 2π ∫ ω_z R dR

    zz, wz_ = gauss_legendre_nodes(z_a, z_b, n)
    R = _F(Rt(zz))
    dz = 1e-6 * max(1.0, abs(z_b - z_a))
    dR = (_F(Rt(zz + dz)) - _F(Rt(zz - dz))) / (2.0 * dz)  # R_t′(z)
    om = _F(W(np.stack([R, 0.0 * R, zz]), t))  # at φ = 0: ω_R = ω_x, ω_z
    side = 2.0 * np.pi * float(np.sum(R * (om[0] - dR * om[2]) * wz_))  # wall flux with n outward
    lower, upper = -cap(z_a), cap(z_b)
    return TubeFlux(lower, side, upper, lower + side + upper)  # Eq. (5.4): −Γ_lower + Γ_upper (+ side = 0)


# Dunavant degree-5 rule on a triangle (barycentric coordinates, weights summing to 1)
_TRI_B = np.array([[1 / 3, 1 / 3, 1 / 3],
                   [0.059715871789770, 0.470142064105115, 0.470142064105115],
                   [0.470142064105115, 0.059715871789770, 0.470142064105115],
                   [0.470142064105115, 0.470142064105115, 0.059715871789770],
                   [0.797426985353087, 0.101286507323456, 0.101286507323456],
                   [0.101286507323456, 0.797426985353087, 0.101286507323456],
                   [0.101286507323456, 0.101286507323456, 0.797426985353087]])
_TRI_W = np.array([0.225, *([0.132394152788506] * 3), *([0.125939180544827] * 3)])


def _refine_triangles(A, B, C, n_sub: int):
    """Split triangles (3, T) into n_sub² similar sub-triangles with the same orientation → three (3, T·n_sub²)."""
    if n_sub <= 1:
        return A, B, C
    P = lambda i, j: A + (i / n_sub) * (B - A) + (j / n_sub) * (C - A)  # noqa: E731
    As, Bs, Cs = [], [], []
    for i in range(n_sub):
        for j in range(n_sub - i):
            As.append(P(i, j)), Bs.append(P(i + 1, j)), Cs.append(P(i, j + 1))
            if i + j <= n_sub - 2:
                As.append(P(i + 1, j)), Bs.append(P(i + 1, j + 1)), Cs.append(P(i, j + 1))
    return np.concatenate(As, axis=1), np.concatenate(Bs, axis=1), np.concatenate(Cs, axis=1)


def _triangle_flux(W, A, B, C, t, n_sub=1):
    """Σ ∫ ω·n dA over triangles (A, B, C) (3, T), n along (B − A) × (C − A); returns (flux, vector-area sum)."""
    A, B, C = _refine_triangles(_F(A), _F(B), _F(C), n_sub)
    nA = 0.5 * np.cross((B - A).T, (C - A).T).T  # vector area of each triangle (3, T)
    Wq = np.zeros_like(A)
    for b, w in zip(_TRI_B, _TRI_W):
        Wq += w * _F(W(b[0] * A + b[1] * B + b[2] * C, t))
    return float(np.sum(Wq * nA)), nA.sum(axis=1)


def tube_flux_budget_traced(omega, tube: Tube, k_a: int = 0, k_b: int | None = None, t: float = 0.0,
                            n_sub: int = 4, **p) -> dict:
    """Gauss' budget (5.4) on a **traced** (not necessarily axisymmetric) vortex tube from :func:`vortex_tube`.

    Book: §5.1, Eq. (5.4), Fig. 5.1. The closed surface is a polyhedron built from the traced lines: the two ends are
    fans of triangles from the centroid of cross-sections k_a and k_b, the side is the band of triangles between
    neighbouring lines. Each triangle is split into ``n_sub``² and integrated with the 7-point degree-5 Dunavant rule.
    Because the polyhedron is closed, ``total`` → 0 up to quadrature error for any solenoidal ω; ``side`` → 0 as
    ``n_lines`` grows (chords → the tube).

    Parameters: omega ω(x, t) or a preset name (+ keywords ``p``); tube; k_a < k_b cross-section indices (default:
    first and last); t [s]; n_sub. Returns dict(flux = :class:`TubeFlux` (outward), area_lower, area_upper [m²],
    mean_omega_lower, mean_omega_upper [1/s]).

    Validation: V4 twisted narrowing Gaussian tube: total ≈ 1e-12, lower = −upper; the "broken" field fails.
    Label: conserved.
    """
    W = _omega_callable(omega, **p)
    L = tube.array
    k_b = L.shape[2] - 1 if k_b is None else int(k_b)
    nl = L.shape[0]
    jn = np.roll(np.arange(nl), -1)

    def cap(k):
        P = L[:, :, k].T  # (3, nl)
        c = P.mean(axis=1, keepdims=True)
        return _triangle_flux(W, np.repeat(c, nl, axis=1), P, P[:, jn], t, n_sub)  # n along +ω

    Fa, Aa = cap(k_a)
    Fb, Ab = cap(k_b)
    side = 0.0
    for k in range(k_a, k_b):
        P0, P1 = L[:, :, k].T, L[:, :, k + 1].T
        f1, _ = _triangle_flux(W, P0, P0[:, jn], P1[:, jn], t, n_sub)  # outward: e_θ × e_s = e_R
        f2, _ = _triangle_flux(W, P0, P1[:, jn], P1, t, n_sub)
        side += f1 + f2
    aa, ab = float(np.linalg.norm(Aa)), float(np.linalg.norm(Ab))
    return dict(flux=TubeFlux(-Fa, side, Fb, -Fa + side + Fb), area_lower=aa, area_upper=ab,
                mean_omega_lower=Fa / aa, mean_omega_upper=Fb / ab)


# ======================================================================================================================
# §5.2 Kelvin's circulation theorem, (5.8)–(5.11)
# ======================================================================================================================
def material_loop(u, pts0, t_eval, rtol: float = 1e-12, atol: float = 1e-14, method: str = "DOP853") -> np.ndarray:
    """Advect every point of a closed loop with the flow (a material contour C, Fig. 5.4): dx/dt = u(x, t).

    Book: §5.2, Fig. 5.4 and the text after (5.8) (the closed curve is made of the same fluid elements at all times).
    All N points are integrated together in **one** ``solve_ivp`` system (state flattened, ``y.reshape(d, N)``); the
    labels s_k stay fixed, so the loop remains sampled uniformly in its material label and the spectral line integrals
    stay accurate while the loop is smooth. Long runs in strong strain stretch the loop (report :func:`loop_length`).

    Parameters: u u(x, t) [m/s] (vectorised over (d, N)); pts0 (d, N) [m]; t_eval (T,) [s] (first entry = start);
    rtol, atol, method: ``solve_ivp`` settings (tight defaults: the circulation error of a stretched loop is set by the
    position error, ~4e-11 m²/s at rtol 1e-10 and ~3e-13 at 1e-12 for the cellular loop over 6 turnovers).
    Returns (T, d, N) positions [m].
    Label: converged.
    """
    P0 = _F(pts0)
    d, N = P0.shape
    te = np.atleast_1d(_F(t_eval))

    def rhs(t, y):
        return _F(u(y.reshape(d, N), t)).ravel()  # dx/dt = u(x, t) for every loop point

    if te.size == 1:
        return P0[None]
    sol = solve_ivp(rhs, (te[0], te[-1]), P0.ravel(), method=method, rtol=rtol, atol=atol, t_eval=te)
    if not sol.success:
        raise RuntimeError(f"material_loop integration failed: {sol.message}")
    return sol.y.T.reshape(te.size, d, N)


def material_circulation(u, pts0, t_eval, rtol: float = 1e-12, atol: float = 1e-14, return_loops: bool = False):
    """Circulation Γ(t) = ∮_C(t) u·dx of a material loop, the quantity of Kelvin's theorem DΓ/Dt = 0, Eq. (5.8).

    Book: §5.2, Eq. (5.8): in an inviscid, barotropic flow with conservative body forces, seen from a nonrotating
    frame, the circulation round a closed curve moving with the fluid is constant. The loop is advected by
    :func:`material_loop`; Γ at each time by the periodic trapezoid rule in the material label.

    Parameters: as :func:`material_loop`. Returns Γ (T,) [m²/s] (and the loops (T, d, N) if ``return_loops``).

    Validation: V4 steady Euler flows (cellular ψ = sin x sin y; Gaussian vortex; ABC) — Γ constant to ~1e-12 while the
    loop stretches; V1 Lamb–Oseen circle: Γ₀(1 − e^{−r²/4ν(t+t₀)}) ((5.11) at work); V3 loop resolution. Label:
    conserved, analytic, converged.
    """
    te = np.atleast_1d(_F(t_eval))
    loops = material_loop(u, pts0, te, rtol, atol)
    G = np.array([loop_circulation(u, loops[i], te[i]) for i in range(te.size)])  # Γ(t) = ∮ u·dx
    return (G, loops) if return_loops else G


class KelvinRate(NamedTuple):
    """The two integrals of (5.9) [m²/s²]: ``acceleration`` ∮(Du_i/Dt)dx_i, ``contour`` ∮u_i D(dx_i)/Dt = ∮u_i du_i;
    ``total`` = their sum = DΓ/Dt."""

    acceleration: float
    contour: float
    total: float


def kelvin_rate_terms(u, pts, t: float = 0.0, h: float = 1e-4, ht: float | None = None) -> KelvinRate:
    """Rate of change of the circulation of a material loop split as in Eq. (5.9).

    Book: §5.2, Eq. (5.9): DΓ/Dt = D/Dt ∮_C u_i dx_i = ∮_C (Du_i/Dt) dx_i + ∮_C u_i (D/Dt)(dx_i); with the element
    kinematics D(dx)/Dt = du (Fig. 5.4) the second integral is ∮ u_i du_i = ∮ d(½u_i²) = 0 on a closed loop. Here
    Du/Dt comes from ``core.kinematics.acceleration`` (stencils) and du = (dx·∇)u = G·dx (stencils), so the contour
    term is computed, not assumed.

    Parameters: u u(x, t) [m/s]; pts (d, N); t [s]; h [m], ht [s] stencil steps. Returns :class:`KelvinRate`.
    Validation: V1 contour ≈ 0 (round-off × stencil); total = dΓ/dt of :func:`material_circulation` (central
    differences). Label: analytic, converged.
    """
    P = _F(pts)
    a = _F(acceleration(u, P, t, h, ht).a)  # Du/Dt at the loop points
    acc = loop_line_integral(a, P)  # ∮ (Du_i/Dt) dx_i
    G = _F(velocity_gradient_at(u, P, t, h))  # (d, d, N)
    T = loop_tangent(P)
    du = np.einsum("ijn,jn->in", G, T)  # D(dx)/Dt = du = G·dx  (Fig. 5.4)
    U = _F(u(P, t))
    con = float(np.sum(U * du) / P.shape[-1])  # ∮ u_i du_i = ∮ d(½u²) = 0
    return KelvinRate(acc, con, acc + con)  # Eq. (5.9)


class KelvinForces(NamedTuple):
    """The three line integrals of (5.10) [m²/s²]: ``pressure`` −∮dp/ρ, ``body`` ∮g·dx (−∮dΦ when conservative),
    ``viscous`` ∮(1/ρ)(∂σ_ij/∂x_j)dx_i; the property ``total`` is their sum."""

    pressure: float
    body: float
    viscous: float

    @property
    def total(self) -> float:
        return self.pressure + self.body + self.viscous


def _field_or_const(f, P, t):
    return _F(f(P, t)) if callable(f) else float(f)


def kelvin_force_terms(pts, grad_p, rho, grad_Phi=None, visc_force=None, t: float = 0.0, g=None) -> KelvinForces:
    """The force integrals of Kelvin's proof, Eq. (5.10): ∮(Du_i/Dt)dx_i = −∮dp/ρ − ∮dΦ + ∮(1/ρ)(∂σ_ij/∂x_j)dx_i.

    Book: §5.2, Eqs. (5.10)–(5.11). The pressure integral vanishes for a **barotropic** fluid because dp/ρ(p) is the
    exact differential of the pressure function (4.67) — "ρ and p single valued" alone is not enough (a baroclinic
    field gives −∮dp/ρ = ∫(∇ρ × ∇p/ρ²)·n dA ≠ 0 by Stokes); the body-force integral vanishes when g = −∇Φ; what is
    left is (5.11).

    Parameters
    ----------
    pts : (d, N) loop [m]
    grad_p : ∇p(x, t) [Pa/m] callable (or None → 0; :func:`stencil_gradient_fn` builds one from p)
    rho : ρ(x, t) [kg/m³] callable or constant
    grad_Phi : ∇Φ(x, t) [m/s²] of a conservative body force g = −∇Φ (4.18), or None
    visc_force : net viscous force per unit **mass** (1/ρ)∂σ_ij/∂x_j (x, t) [m/s²] (e.g. ν∇²u), or None
    t : time [s];  g : any body force per unit mass g(x, t) [m/s²] (conservative or not), or None

    Returns
    -------
    :class:`KelvinForces` [m²/s²] (``.total`` their sum).

    Validation: V1 barotropic ρ(p) → pressure 0; baroclinic → −∫(∇ρ × ∇p/ρ²)·n dA (Stokes); Lamb–Oseen circle:
    viscous = ∂Γ/∂t of (5.11). Label: analytic.
    """
    P = _F(pts)
    d = P.shape[0]
    rho_ = _field_or_const(rho, P, t)
    pres = 0.0
    if grad_p is not None:
        pres = -loop_line_integral(_F(grad_p(P, t))[:d] / rho_, P)  # −∮ (1/ρ)(∂p/∂x_i) dx_i = −∮ dp/ρ
    body = 0.0
    if grad_Phi is not None:
        body += -loop_line_integral(_F(grad_Phi(P, t))[:d], P)  # −∮ dΦ
    if g is not None:
        body += loop_line_integral(_F(g(P, t))[:d], P)  # ∮ g·dx
    visc = 0.0
    if visc_force is not None:
        visc = loop_line_integral(_F(visc_force(P, t))[:d], P)  # ∮ (1/ρ)(∂σ_ij/∂x_j) dx_i   (5.11)
    return KelvinForces(pres, body, visc)  # Eq. (5.10)


def _plane_scalar(f):
    """A plane scalar field given as f(x, y) (two required arguments) or in the field convention f(x, t) → f(x, t)."""
    try:
        params = [q for q in inspect.signature(f).parameters.values()
                  if q.default is inspect.Parameter.empty and q.kind in (q.POSITIONAL_ONLY, q.POSITIONAL_OR_KEYWORD)]
    except (TypeError, ValueError):
        params = [None]
    if len(params) >= 2:
        return lambda x, t=0.0: _F(f(_F(x)[0], _F(x)[1]))
    return f


def pressure_torque_on_element(p="linear", rho="linear", center=(0.0, 0.0), radius: float = 0.01, n: int = 4000,
                               nr: int = 400, grad_rho=(10.0, 0.0), grad_p=(0.0, -9810.0), rho0: float = 1000.0,
                               p0: float = 1.0e5, t: float = 0.0, h: float | None = None) -> dict:
    """Net pressure force, centre of mass and pressure torque on a small disc of fluid (Fig. 5.6): the baroclinic
    spin-up 2·torque/I_G tends to (∇ρ × ∇p)/ρ² as the disc shrinks.

    Book: §5.2, Fig. 5.6 and restriction (3): when isopycnals are parallel to isobars (barotropic) the resultant
    pressure force passes through the centre of mass G; when they cross (baroclinic) it misses G and its torque changes
    the vorticity. Our D06: for linear fields the pressure force −∇p·A acts at the geometric centre while G sits at
    R²∇ρ/4ρ₀ from it, so torque = πR⁴(∇ρ × ∇p)/4ρ₀, I_G ≈ πρ₀R⁴/2 and dω/dt = 2 dΩ/dt = 2·torque/I_G
    → (∇ρ × ∇p)/ρ₀² — the baroclinic term (5.28) of the vorticity equation.

    Parameters
    ----------
    p, rho : "linear" (p = p₀ + ∇p·(x − c), ρ = ρ₀ + ∇ρ·(x − c)) or callables — f(x, y) or the field form f(x, t)
    center : disc centre (2,) [m];  radius : R [m]
    n : rim points (midpoint rule, spectral) for force and torque;  nr : Gauss–Legendre radial nodes (4·nr angles) for
    the mass integrals;  grad_rho [kg/m⁴], grad_p [Pa/m], rho0 [kg/m³], p0 [Pa] : the "linear" fields;  t [s]
    h : stencil step for the point value of ∇ρ × ∇p/ρ² (default R/100)

    Returns
    -------
    dict: ``force`` (2,) [N/m], ``x_G`` (2,) [m], ``offset`` x_G − centre [m], ``torque`` about G (z) [N m/m],
    ``mass`` [kg/m], ``I_G`` [kg m²/m], ``spin_up`` = 2·torque/I_G [1/s²], ``baroclinic`` (∇ρ × ∇p)_z/ρ(centre)²
    [1/s²], ``ratio`` spin_up/baroclinic (NaN when the baroclinic term is 0).

    Validation: V1 defaults (R = 0.01, ∇ρ = (10, 0), ∇p = (0, −9810)): x_G = (2.5e-7, 0), torque −7.7047e-7, ratio
    1 + 1.25e-9 = 1/(1 − |∇ρ|²R²/8ρ₀²); barotropic → torque 0; V3 ratio − 1 = O(R²). Label: analytic, converged.
    """
    c = _F(center)
    if isinstance(p, str):
        gp_, p0_ = _F(grad_p), float(p0)
        pf = lambda x, tt=0.0: p0_ + gp_[0] * (_F(x)[0] - c[0]) + gp_[1] * (_F(x)[1] - c[1])  # noqa: E731
    else:
        pf = _plane_scalar(p)
    if isinstance(rho, str):
        gr_, r0_ = _F(grad_rho), float(rho0)
        rf = lambda x, tt=0.0: r0_ + gr_[0] * (_F(x)[0] - c[0]) + gr_[1] * (_F(x)[1] - c[1])  # noqa: E731
    elif callable(rho):
        rf = _plane_scalar(rho)
    else:
        rv = float(rho)
        rf = lambda x, tt=0.0: rv + 0.0 * _F(x)[0]  # noqa: E731
    phi = 2.0 * np.pi * (np.arange(n) + 0.5) / n
    nrm = np.stack([np.cos(phi), np.sin(phi)])
    X = c[:, None] + radius * nrm
    ds = 2.0 * np.pi * radius / n
    pc = float(_F(pf(c, t)))  # a uniform pressure exerts no net force or torque: subtract it (round-off)
    dF = -(_F(pf(X, t)) - pc) * nrm * ds  # pressure force −p n ds on each rim element (n outward)
    force = dF.sum(axis=1)
    r, wr = gauss_legendre_nodes(0.0, radius, nr)
    nth = 4 * nr
    th = 2.0 * np.pi * (np.arange(nth) + 0.5) / nth
    Rr, TH = np.meshgrid(r, th, indexing="ij")
    Xi = np.stack([c[0] + Rr * np.cos(TH), c[1] + Rr * np.sin(TH)]).reshape(2, -1)
    w = (wr[:, None] * Rr * (2.0 * np.pi / nth)).ravel()
    rho_i = _F(rf(Xi, t))
    M = float(np.sum(rho_i * w))
    xG = np.sum(rho_i * w * Xi, axis=1) / M
    rel = X - xG[:, None]
    torque = float(np.sum(rel[0] * dF[1] - rel[1] * dF[0]))  # z-component of Σ (x − x_G) × dF
    IG = float(np.sum(rho_i * w * np.sum((Xi - xG[:, None]) ** 2, axis=0)))
    spin = 2.0 * torque / IG  # dω/dt = 2 dΩ/dt = 2 torque / I_G
    if isinstance(p, str) and isinstance(rho, str):  # linear fields: the exact gradients (no stencil round-off)
        gr_, gp_ = _F(grad_rho), _F(grad_p)
        bz = float((gr_[0] * gp_[1] - gr_[1] * gp_[0]) / float(rho0) ** 2)
    else:
        hh = radius / 100.0 if h is None else float(h)
        bz = float(_F(baroclinic_term(rf, pf, c, t, hh))[2])
    return dict(force=force, x_G=xG, offset=xG - c, torque=torque, mass=M, I_G=IG, spin_up=spin, baroclinic=bz,
                ratio=spin / bz if bz != 0.0 else np.nan)


def frozen_in_check(u, x0, delta0=None, t_span=(0.0, 1.0), t_eval=None, h: float = 1e-5, nu: float = 0.0,
                    rtol: float = 1e-11, atol: float = 1e-13, h_visc: float | None = None) -> dict:
    """Helmholtz's first theorem as a computation: a short material element δx that starts along ω stays along the
    vorticity, with |ω|/|δx| constant, when the flow is inviscid and barotropic.

    Book: §5.3, Helmholtz theorem (1), "vortex lines move with the fluid" (proved from Kelvin with Fig. 5.7; the book
    notes the same follows from the field equation): with ν = 0, (5.13) gives Dω/Dt = (ω·∇)u = Gω along a particle
    path — the same linear equation as D(δx)/Dt = Gδx (ch03 (3.10)); equal directions at t₀ ⇒ proportional for ever.
    One ``solve_ivp`` integrates the particle x, the element δx and ω (+ ν∇²ω of the Eulerian field by stencils when
    ν > 0), G by stencils. The Eulerian vorticity at x(t) is returned too (``omega_field``).

    Parameters
    ----------
    u : u(x, t) [m/s] (3-D);  x0 : start point [m];  delta0 : initial element (3,) [m] (default ω(x0)/|ω(x0)|·1e-3)
    t_span : (t0, t1) [s];  t_eval : output times (default 21 points);  h : stencil step [m];  nu : ν [m²/s]
    rtol, atol : ODE tolerances (rtol is raised to ≥ 1e-8 when ν > 0: the nested-stencil ν∇²ω carries ~1e-9 noise
    that would stall a tighter step control);  h_visc : step [m] of that ν∇²ω stencil (default 1e-2·max(|x0|, 1e-3))

    Returns
    -------
    dict: ``t``, ``x`` (3, T), ``delta`` (3, T), ``omega`` (3, T) (integrated), ``omega_field`` (3, T), ``angle`` (T,)
    [rad] between δx and ω, ``ratio`` (T,) |ω|/|δx| normalised to 1 at t₀, ``ratio_field`` (same with the field ω).

    Validation: V1 ABC flow and the inviscid stretched Gaussian vortex: angle < 1e-8, ratio − 1 < 1e-8 (both routes);
    viscous Burgers vortex with ν: ratio drifts (wrong-hypothesis check). Label: analytic, conserved.
    """
    t0, t1 = float(t_span[0]), float(t_span[1])
    te = np.linspace(t0, t1, 21) if t_eval is None else np.atleast_1d(_F(t_eval))
    x0 = _F(x0)
    W = vorticity_field(u, h)
    w0 = _F(W(x0, t0))
    d0 = (1e-3 * w0 / np.linalg.norm(w0)) if delta0 is None else _F(delta0)
    lap_h = 1e-2 * max(float(np.linalg.norm(x0)), 1e-3) if h_visc is None else float(h_visc)
    W_lap = vorticity_field(u, lap_h)  # same step inside and outside: noise ~ eps/h³ stays small
    if nu > 0.0:
        rtol = max(rtol, 1e-8)
        atol = max(atol, 1e-11)

    def rhs(t, y):
        x, dx, om = y[:3], y[3:6], y[6:]
        G = _F(velocity_gradient_at(u, x, t, h))
        dom = G @ om  # (ω·∇)u = G ω
        if nu > 0.0:
            dom = dom + nu * _F(st.laplacian(W_lap, x, t, lap_h, (3,)))  # + ν∇²ω of the field
        return np.concatenate([_F(u(x, t)), G @ dx, dom])  # dx/dt = u,  D(δx)/Dt = G δx,  Dω/Dt

    sol = solve_ivp(rhs, (t0, t1), np.concatenate([x0, d0, w0]), method="DOP853", rtol=rtol, atol=atol, t_eval=te)
    if not sol.success:
        raise RuntimeError(sol.message)
    X, D, Om = sol.y[:3], sol.y[3:6], sol.y[6:]
    Omf = np.stack([_F(W(X[:, i], te[i])) for i in range(te.size)], axis=1)
    nd, nw, nf = np.linalg.norm(D, axis=0), np.linalg.norm(Om, axis=0), np.linalg.norm(Omf, axis=0)
    ang = np.arctan2(np.linalg.norm(np.cross(D.T, Om.T), axis=1), np.sum(D * Om, axis=0))
    ratio = (nw / nd) / (nw[0] / nd[0])
    ratio_f = (nf / nd) / (nf[0] / nd[0])
    return dict(t=te, x=X, delta=D, omega=Om, omega_field=Omf, angle=np.abs(ang), ratio=ratio, ratio_field=ratio_f)


# ======================================================================================================================
# §5.4 and §5.6 the vorticity equation (5.12)–(5.13), (5.18)–(5.30)
# ======================================================================================================================
def _sym_curl(F, X):
    return sp.Matrix([sp.diff(F[2], X[1]) - sp.diff(F[1], X[2]), sp.diff(F[0], X[2]) - sp.diff(F[2], X[0]),
                      sp.diff(F[1], X[0]) - sp.diff(F[0], X[1])])


def _sym_grad(f, X):
    return sp.Matrix([sp.diff(f, xi) for xi in X])


def _sym_div(F, X):
    return sum(sp.diff(F[i], X[i]) for i in range(3))


def _sym_dir(a, F, X):
    """(a·∇)F for a vector F."""
    return sp.Matrix([sum(a[j] * sp.diff(F[i], X[j]) for j in range(3)) for i in range(3)])


def _sym_lap(F, X):
    return sp.Matrix([sum(sp.diff(F[i], xj, 2) for xj in X) for i in range(3)])


def vorticity_equation_sym(u_exprs, coords, t: sp.Symbol, nu=None, p_expr=None, Phi_expr=None, rho=None) -> dict:
    """The curl of the incompressible Navier–Stokes equation, term by term, Eq. (5.12) → (5.13) (sympy).

    Book: §5.4, Eq. (5.12): ∇ × {Du/Dt = −(1/ρ)∇p + g + ν∇²u}; for constant ρ and g = −∇Φ the pressure and gravity
    curls vanish (gradients); the Lamb identity gives ∇ × {∂u/∂t + (u·∇)u} = ∂ω/∂t + ∇ × (ω × u) (the book writes
    ∇(u·u) for ∇(½u·u) — harmless, its curl is 0); ∇ × ∇²u = ∇²(∇ × u); (B.3.10)
    ∇ × (ω × u) = (u·∇)ω − (ω·∇)u + ω(∇·u) − u(∇·ω), with ∇·u = ∇·ω = 0, gives (5.13).

    Parameters: u_exprs three sympy expressions (Cartesian) in ``coords`` (x, y, z) and ``t``; nu ν (default Symbol
    "nu"); p_expr, Phi_expr optional scalars; rho constant (default Symbol "rho").

    Returns
    -------
    dict of sympy Matrices (simplified): ``curl_local`` ∇×∂u/∂t = ∂ω/∂t, ``curl_advective`` ∇×(u·∇)u,
    ``curl_pressure`` ∇×(−∇p/ρ) (→ 0), ``curl_gravity`` ∇×(−∇Φ) (→ 0), ``curl_viscous`` ∇×(ν∇²u), ``lamb_curl``
    ∇×(ω×u), ``identity_B310`` (general (B.3.10) residual, → 0 for any fields), ``residual_513`` (curl of NS minus
    (5.13): → 0 when ∇·u = 0); also ``omega``, ``lhs_513`` Dω/Dt, ``rhs_513`` (ω·∇)u + ν∇²ω, ``div_u``, ``div_omega``.
    Label: symbolic.
    """
    X = list(coords)
    nu = sp.Symbol("nu", positive=True) if nu is None else nu
    rho = sp.Symbol("rho", positive=True) if rho is None else rho
    u = sp.Matrix(u_exprs)
    w = _sym_curl(u, X)
    out = dict(omega=w)
    out["curl_local"] = sp.simplify(_sym_curl(sp.diff(u, t), X))
    out["curl_advective"] = sp.simplify(_sym_curl(_sym_dir(u, u, X), X))
    out["curl_pressure"] = (sp.simplify(_sym_curl(-_sym_grad(p_expr, X) / rho, X)) if p_expr is not None
                            else sp.zeros(3, 1))
    out["curl_gravity"] = (sp.simplify(_sym_curl(-_sym_grad(Phi_expr, X), X)) if Phi_expr is not None
                           else sp.zeros(3, 1))
    out["curl_viscous"] = sp.simplify(_sym_curl(nu * _sym_lap(u, X), X))
    out["lamb_curl"] = sp.simplify(_sym_curl(w.cross(u), X))
    out["identity_B310"] = sp.simplify(out["lamb_curl"] - (_sym_dir(u, w, X) - _sym_dir(w, u, X)
                                                           + w * _sym_div(u, X) - u * _sym_div(w, X)))  # (B.3.10)
    out["lhs_513"] = sp.simplify(sp.diff(w, t) + _sym_dir(u, w, X))  # Dω/Dt
    out["rhs_513"] = sp.simplify(_sym_dir(w, u, X) + nu * _sym_lap(w, X))  # (ω·∇)u + ν∇²ω   (5.13)
    curl_ns = (out["curl_local"] + out["curl_advective"] - out["curl_pressure"] - out["curl_gravity"]
               - out["curl_viscous"])
    out["residual_513"] = sp.simplify(curl_ns - (out["lhs_513"] - out["rhs_513"]))
    out["div_u"], out["div_omega"] = sp.simplify(_sym_div(u, X)), sp.simplify(_sym_div(w, X))
    return out


class VorticityTerms(NamedTuple):
    """Terms of (5.13) at points [1/s²]: Dω/Dt = ``local`` + ``advective`` = ``stretching_tilting`` + ``diffusion``;
    ``residual`` = local + advective − stretching_tilting − diffusion (→ 0 for a Navier–Stokes solution)."""

    local: np.ndarray
    advective: np.ndarray
    stretching_tilting: np.ndarray
    diffusion: np.ndarray
    residual: np.ndarray


class VorticityBudget(NamedTuple):
    """Terms of (5.30) at points [1/s²]: ``local`` ∂ω/∂t, ``advective`` (u·∇)ω, ``stretching_tilting`` (ω·∇)u,
    ``planetary`` 2(Ω·∇)u, ``baroclinic`` ∇ρ × ∇p/ρ², ``diffusion`` ν∇²ω; ``residual`` = local + advective −
    stretching_tilting − planetary − baroclinic − diffusion."""

    local: np.ndarray
    advective: np.ndarray
    stretching_tilting: np.ndarray
    planetary: np.ndarray
    baroclinic: np.ndarray
    diffusion: np.ndarray
    residual: np.ndarray


def vorticity_budget(u, x, t: float = 0.0, nu: float = 0.0, Omega=(0.0, 0.0, 0.0), rho=None, p=None,
                     h: float = 1e-3, ht: float | None = None, omega_fn: Callable | None = None) -> VorticityBudget:
    """Every term of the vorticity equation in a rotating frame with baroclinic generation, Eq. (5.30), at points.

    Book: §5.6, Eq. (5.30): Dω/Dt = (ω + 2Ω)·∇u + (1/ρ²)∇ρ × ∇p + ν∇²ω — u, ω relative to the frame rotating at the
    constant rate Ω; the stretching/tilting acts on the **absolute** vorticity ω + 2Ω (planetary part 2(Ω·∇)u, the
    component equations Dω_z/Dt = 2Ω∂w/∂z … after (5.32)); the baroclinic term (5.28) vanishes when ∇ρ ∥ ∇p;
    molecular diffusion ν∇²ω. With Ω = 0 and no baroclinic term it is (5.13). The "Boussinesq" (5.30) keeps the full
    1/ρ in the pressure term with ∇·u = 0 and constant ν (a mixed approximation; to leading order the baroclinic term
    is ∇ρ′ × g/ρ₀, cf. (4.86)).

    Parameters
    ----------
    u : u(x, t) [m/s], 3 components (plane fields: wrap with :func:`planar_field_3d`)
    x : point(s) (3,) or (3, N) [m];  t [s];  nu : ν [m²/s];  Omega : frame rotation Ω (3,) [rad/s]
    rho, p : ρ(x, t) [kg/m³], p(x, t) [Pa] callables for the baroclinic term (both or neither)
    h, ht : stencil steps [m], [s] (nested stencils: ω from u, then ∇ω and ∇²ω — use h ≈ 1e-3 × the flow scale)
    omega_fn : exact ω(x, t) if known (removes one level of nesting)

    Returns
    -------
    :class:`VorticityBudget` (arrays (3,) or (3, N)) [1/s²].

    Validation: V1 residual → 0 for Lamb–Oseen (stretching 0, local = diffusion), Taylor–Green, Burgers (steady:
    advective = stretching + diffusion, each nonzero), Hill (ν = 0), rotating stretched column (planetary);
    V3 order 2 in h; V2 :func:`vorticity_budget_sym`. Label: analytic, converged, symbolic.
    """
    x_ = _F(x)
    W = vorticity_field(u, h) if omega_fn is None else omega_fn
    U = st.ev(u, x_, t, (3,))
    Om = st.ev(W, x_, t, (3,))
    G = st.grad_vector(u, x_, t, h)  # ∂u_i/∂x_j
    GW = st.grad_vector(W, x_, t, h)  # ∂ω_i/∂x_j
    local = st.ddt(W, x_, t, ht, (3,))  # ∂ω/∂t
    adv = np.einsum("ij...,j...->i...", GW, U)  # (u·∇)ω
    stretch = np.einsum("ij...,j...->i...", G, Om)  # (ω·∇)u
    Omg = _F(Omega).reshape((3,) + (1,) * (x_.ndim - 1))
    plan = np.einsum("ij...,j...->i...", G, 2.0 * Omg * np.ones_like(Om))  # 2(Ω·∇)u
    if rho is not None and p is not None:
        baro = _F(baroclinic_term(rho, p, x_, t, h))  # (1/ρ²) ∇ρ × ∇p   (5.28)
    else:
        baro = np.zeros_like(Om)
    diff = float(nu) * st.laplacian(W, x_, t, h, (3,))  # ν∇²ω
    res = local + adv - stretch - plan - baro - diff  # Eq. (5.30)
    return VorticityBudget(local, adv, stretch, plan, baro, diff, res)


def vorticity_terms(u, x, t: float = 0.0, nu: float = 0.0, h: float = 1e-3, ht: float | None = None,
                    omega_fn: Callable | None = None) -> VorticityTerms:
    """Terms of the vorticity equation for constant density in an inertial frame, Eq. (5.13), at points.

    Book: §5.4, Eq. (5.13): Dω/Dt = (ω·∇)u + ν∇²ω — vorticity is carried (u·∇)ω, stretched and tilted (ω·∇)u and
    diffused ν∇²ω; pressure and gravity do not appear (they act through the centre of mass). Special case of
    :func:`vorticity_budget` (Ω = 0, no baroclinic term). Parameters and validation as there. Label: analytic,
    converged.
    """
    b = vorticity_budget(u, x, t, nu, (0.0, 0.0, 0.0), None, None, h, ht, omega_fn)
    return VorticityTerms(b.local, b.advective, b.stretching_tilting, b.diffusion,
                          b.local + b.advective - b.stretching_tilting - b.diffusion)  # Eq. (5.13)


def vorticity_budget_sym(u_exprs, coords, t: sp.Symbol, nu=None, Omega=(0, 0, 0), rho_expr=None,
                         p_expr=None) -> dict:
    """Symbolic terms of (5.30): dict(``terms`` — sympy Matrices local, advective, stretching_tilting, planetary,
    baroclinic, diffusion; ``residual_530`` (simplified); ``reduces_to_513`` — with Ω → 0 and no ρ, p the residual
    equals that of (5.13) (checked symbolically)); the six terms are also top-level keys.
    Book: §5.6, Eq. (5.30) and §5.4, (5.13). Label: symbolic."""
    X = list(coords)
    nu = sp.Symbol("nu", positive=True) if nu is None else nu
    u = sp.Matrix(u_exprs)
    w = _sym_curl(u, X)
    Om = sp.Matrix(Omega)
    terms = dict(local=sp.diff(w, t), advective=_sym_dir(u, w, X), stretching_tilting=_sym_dir(w, u, X),
                 planetary=_sym_dir(2 * Om, u, X), diffusion=nu * _sym_lap(w, X),
                 baroclinic=(baroclinic_term_sym(rho_expr, p_expr, X) if rho_expr is not None and p_expr is not None
                             else sp.zeros(3, 1)))
    res = sp.simplify(terms["local"] + terms["advective"] - terms["stretching_tilting"] - terms["planetary"]
                      - terms["baroclinic"] - terms["diffusion"])
    res513 = sp.simplify(terms["local"] + terms["advective"] - terms["stretching_tilting"] - terms["diffusion"])
    # Ω → 0 and ∇ρ × ∇p → 0: the (5.30) residual must become the (5.13) residual
    red = (res + terms["planetary"] + terms["baroclinic"]).subs({o: 0 for o in Om.free_symbols})
    reduces = bool(sp.simplify(red - res513) == sp.zeros(3, 1))
    return dict(terms=terms, residual_530=res, reduces_to_513=reduces, **terms)


def vorticity_terms_sym(u_exprs, coords, t: sp.Symbol, nu=None) -> dict:
    """Symbolic terms of (5.13) (Ω = 0, no baroclinic term): dict local, advective, stretching_tilting, diffusion,
    residual. Book: §5.4, Eq. (5.13). Label: symbolic."""
    b = vorticity_budget_sym(u_exprs, coords, t, nu)
    out = {k: b[k] for k in ("local", "advective", "stretching_tilting", "diffusion")}
    out["residual"] = b["residual_530"]
    return out


VORTICITY_BUDGET_PRESETS = ("burgers", "lamb_oseen", "taylor_green", "rotating_column", "lock_exchange", "hill")


def _rotating_column_field(Omega: float, alpha: float) -> Callable:
    """u = −(α/2)(x, y, 0) + (ζ(t)/2) e_z × x + (0, 0, αz), ζ(t) = 2Ω(e^{αt} − 1): a column at rest at t = 0 in a frame
    rotating at Ω, stretched axially — ζ + 2Ω = 2Ωe^{αt} (inviscid (5.30) with uniform vorticity)."""
    def u(x, t=0.0):
        X = _F(x)
        z_t = 2.0 * Omega * np.expm1(alpha * t)
        return np.stack([-0.5 * alpha * X[0] - 0.5 * z_t * X[1], -0.5 * alpha * X[1] + 0.5 * z_t * X[0], alpha * X[2]])
    return u


def vorticity_budget_preset(name: str, x=None, y=None, z=None, t: float = 0.0, component: int = 2, **p) -> dict:
    """One component of every term of (5.30) for the scenes of E5, E7 and the C09 figure (scalar floats, for parity).

    Book: §5.4, Eq. (5.13); §5.6, Eq. (5.30). Presets (defaults; x, y, z = None → the preset's probe point):
    * "burgers" — Γ = 1e-3 m²/s, α = 1 s⁻¹, ν = 1e-6 m²/s, probe R = 1 mm: steady, advective = stretching + diffusion;
    * "lamb_oseen" — Γ = 0.01, ν = 1e-6, t₀ = 10 s, probe r = 5 mm: local = diffusion (no stretching);
    * "taylor_green" — u = U₀(cos x sin y, −sin x cos y)e^{−2νt}, U₀ = 1, ν = 0.01, probe (0.3, 0.2): local =
      diffusion;
    * "rotating_column" — Ω = 0.5 rad/s, α = 0.2 s⁻¹: a column at rest in the rotating frame at t = 0, stretched
      axially; local = planetary = 2Ωα = 0.2 s⁻²;
    * "lock_exchange" — ρ₁ = 1000, ρ₂ = 1025, δ = 0.1 m, H = 1 m, u = 0 at t = 0, probe (0, H/2): baroclinic
      = 2.4222 s⁻²; ``local`` is the tendency that (5.30) then predicts (= baroclinic), so residual = 0 — the u = 0
      snapshot cannot supply ∂ω/∂t itself;
    * "hill" — A = 1, a = 1, ν = 0, probe (0.5, 0, 0.3) (use component 1 = ω_φ at φ = 0): advective = stretching.

    Returns dict(local, advective, stretching_tilting, planetary, baroclinic, diffusion, residual) floats [1/s²].
    Label: analytic, converged.
    """
    from . import vortices as VX

    Om = (0.0, 0.0, 0.0)
    rho_f = p_f = None
    nu = float(p.get("nu", 0.0))
    if name == "burgers":
        Gam, al, nu = float(p.get("Gamma", 1e-3)), float(p.get("alpha", 1.0)), float(p.get("nu", 1e-6))
        u = VX.burgers_vortex_field(Gam, al, nu)
        probe, h = (1e-3, 0.0, 0.0), 1e-6
    elif name == "lamb_oseen":
        Gam, nu, t0 = float(p.get("Gamma", 0.01)), float(p.get("nu", 1e-6)), float(p.get("t0", 10.0))
        u = planar_field_3d(VX.lamb_oseen_field(Gam, nu, t0))
        probe, h = (5e-3, 0.0, 0.0), 5e-6
    elif name == "taylor_green":
        U0, nu = float(p.get("U0", 1.0)), float(p.get("nu", 0.01))

        def u(x, tt=0.0):
            X = _F(x)
            e = U0 * np.exp(-2.0 * nu * tt)
            return np.stack([e * np.cos(X[0]) * np.sin(X[1]), -e * np.sin(X[0]) * np.cos(X[1]), 0.0 * X[2]])
        probe, h = (0.3, 0.2, 0.0), 1e-3
    elif name == "rotating_column":
        Omz, al = float(p.get("Omega", 0.5)), float(p.get("alpha", 0.2))
        u = _rotating_column_field(Omz, al)
        Om = (0.0, 0.0, Omz)
        probe, h = (0.3, 0.2, 0.1), 1e-3
    elif name == "lock_exchange":
        r1, r2 = float(p.get("rho1", 1000.0)), float(p.get("rho2", 1025.0))
        dl, H, g = float(p.get("delta", 0.1)), float(p.get("H", 1.0)), float(p.get("g", 9.81))
        rb = 0.5 * (r1 + r2)
        rho_f = lambda x, tt=0.0: rb - 0.5 * (r2 - r1) * np.tanh(2.0 * _F(x)[0] / dl)  # noqa: E731
        p_f = lambda x, tt=0.0: rb * g * (H - _F(x)[1]) + 0.0 * _F(x)[0]  # noqa: E731  mean-density hydrostatic
        u = lambda x, tt=0.0: 0.0 * _F(x)  # noqa: E731
        probe, h = (0.0, 0.5 * H, 0.0), 1e-5
    elif name == "hill":
        A, a = float(p.get("A", 1.0)), float(p.get("a", 1.0))
        u = VX.hill_spherical_vortex_field(A, a)
        probe, h, nu = (0.5, 0.0, 0.3), 1e-3, 0.0
    else:
        raise ValueError(f"unknown preset {name!r}; choose from {VORTICITY_BUDGET_PRESETS}")
    pt = np.array([probe[0] if x is None else x, probe[1] if y is None else y, probe[2] if z is None else z], float)
    b = vorticity_budget(u, pt, t, nu, Om, rho_f, p_f, float(p.get("h", h)))
    out = {k: float(_F(getattr(b, k))[component]) for k in VorticityBudget._fields}
    if name == "lock_exchange":
        out["local"] = out["baroclinic"]  # the tendency (5.30) predicts at t = 0⁺ (the u = 0 snapshot has none)
        out["residual"] = (out["local"] + out["advective"] - out["stretching_tilting"] - out["planetary"]
                           - out["baroclinic"] - out["diffusion"])
    return out


def vorticity_divergence(u, x, t: float = 0.0, h: float = 1e-3):
    """∇·ω = ∇·(∇ × u) by nested central differences — zero to round-off, Eq. (5.18).

    Book: §5.6, Eq. (5.18): ω_{i,i} = (ε_inq u_{q,n})_{,i} = ε_inq u_{q,ni} = 0 (antisymmetric ε against the symmetric
    u_{q,ni}); the vorticity field is solenoidal even for compressible, unsteady flow. The mixed difference operators
    commute exactly, so the discrete identity also holds to round-off. Returns [1/(s m)]. Label: analytic.
    """
    return _S(st.div(vorticity_field(u, h), _F(x), t, h))


def rotating_lamb_form_terms(u, p, rho, Phi, x, t: float = 0.0, nu: float = 0.0, Omega=(0.0, 0.0, 0.0),
                             h: float = 1e-4, ht: float | None = None) -> dict:
    """The five terms of the rotating-frame momentum equation in Lamb form, Eq. (5.25), at points [m/s²].

    Book: §5.6, Eq. (5.25): ∂u_i/∂t + (½u_j² + Φ)_{,i} − ε_ijk u_j(ω_k + 2Ω_k) = −(1/ρ)p_{,i} − νε_ijk ω_{k,j}, obtained
    from (5.20) with (5.21) (Lamb identity), (5.23) (ν∇²u = −ν∇ × ω, needs ∇·u = 0), (5.24) and g = −∇Φ.

    Parameters: u (3-D) [m/s]; p [Pa]; rho [kg/m³] (callable or constant); Phi Φ(x, t) [m²/s²] (effective gravity
    potential, or None); x; t; nu; Omega (3,); h, ht.
    Returns dict ``local`` ∂u/∂t, ``grad_B`` ∇(½u² + Φ), ``lamb_abs`` −u × (ω + 2Ω), ``pressure`` −(1/ρ)∇p,
    ``viscous_curl`` −ν∇ × ω (the right-side terms carry their signs), ``residual`` = local + grad_B + lamb_abs −
    pressure − viscous_curl (equals the (5.20) residual for ∇·u = 0). Label: analytic, converged.
    """
    x_ = _F(x)
    W = vorticity_field(u, h)
    U = st.ev(u, x_, t, (3,))
    Om = st.ev(W, x_, t, (3,))
    local = st.ddt(lambda X, T: st.ev(u, X, T, (3,)), x_, t, ht, (3,))
    Phi_ = (lambda X, T: 0.0 * _F(X)[0]) if Phi is None else Phi
    B = lambda X, T: 0.5 * np.sum(st.ev(u, X, T, (3,)) ** 2, axis=0) + st.ev(Phi_, X, T)  # noqa: E731
    gB = st.grad(B, x_, t, h)
    Omg = _F(Omega).reshape((3,) + (1,) * (x_.ndim - 1)) * np.ones_like(Om)
    lamb = -np.cross(U, Om + 2.0 * Omg, axis=0)  # −ε_ijk u_j(ω_k + 2Ω_k)
    rho_ = st.ev(rho, x_, t) if callable(rho) else float(rho)
    pres = -st.grad(p, x_, t, h) / rho_  # −(1/ρ) p_{,i}
    visc = -float(nu) * _F(vorticity_from_gradient(st.grad_vector(W, x_, t, h)))  # −ν ε_ijk ω_{k,j} = −ν(∇ × ω)_i
    res = local + gB + lamb - pres - visc  # Eq. (5.25): left − right
    return dict(local=local, grad_B=gB, lamb_abs=lamb, pressure=pres, viscous_curl=visc, residual=res)


def baroclinic_term(rho, p, x, t: float = 0.0, h: float = 1e-4) -> np.ndarray:
    """Baroclinic vorticity generation (1/ρ²)∇ρ × ∇p at points (3,) or (3, N) [1/s²], Eq. (5.28).

    Book: §5.6, Eq. (5.28): −ε_nqi(p_{,i}/ρ)_{,q} = −(1/ρ)ε_nqi p_{,iq} + (1/ρ²)ε_nqi ρ_{,q} p_{,i} = 0 + (1/ρ²)[∇ρ × ∇p]_n;
    zero for a barotropic fluid (ρ = ρ(p) ⇒ ∇ρ ∥ ∇p). Plane fields (x of shape (2, …)) give (0, 0, ρ_x p_y − ρ_y p_x)/ρ².
    Order matters: ∇ρ × ∇p, not ∇p × ∇ρ (lock exchange with the heavy fluid on the left spins counterclockwise).

    Parameters: rho ρ(x, t), p p(x, t) (callables); x [m]; t [s]; h stencil step [m].
    Validation: V1 ρ = ρ(p) → 0 for arbitrary p; lock-exchange closed form; V2 :func:`baroclinic_term_sym`; V3 order 2.
    Label: analytic, converged, symbolic.
    """
    x_ = _F(x)
    gr = st.pad3(st.grad(rho, x_, t, h))
    gp = st.pad3(st.grad(p, x_, t, h))
    r = st.ev(rho, x_, t)
    return np.cross(gr, gp, axis=0) / r ** 2  # (1/ρ²) ∇ρ × ∇p   (5.28)


def baroclinic_term_sym(rho_expr, p_expr, coords) -> sp.Matrix:
    """Symbolic (1/ρ²)∇ρ × ∇p, Eq. (5.28) (3 Cartesian coordinates). Book: §5.6, Eq. (5.28). Label: symbolic."""
    X = list(coords)
    gr, gp = _sym_grad(rho_expr, X), _sym_grad(p_expr, X)
    return sp.simplify(gr.cross(gp) / rho_expr ** 2)


def stretching_tilting_split(omega, G) -> dict:
    """Split (ω·∇)u = Gω into vortex **stretching** along ω and **tilting** across it, Eqs. (5.31)–(5.32).

    Book: §5.6, Eq. (5.31): (ω·∇)u = ω ∂u/∂s (s along the vortex line, ω = |ω|); its component along e_s,
    ω ∂u_s/∂s, stretches the line ((5.32) first), the components along e_n, e_m, ω ∂u_n/∂s and ω ∂u_m/∂s, tilt it.
    Here ∂u/∂s = G e_s, so stretching = (e_s·G e_s) ω and tilting = Gω − stretching.

    Parameters: omega ω (3,) [1/s]; G velocity gradient (3, 3) [1/s] (G[i, j] = ∂u_i/∂x_j).
    Returns dict: ``e_s`` (3,), ``rate`` e_s·G e_s [1/s] (the stretching rate), ``stretching`` (3,) [1/s²] (∥ ω),
    ``tilting`` (3,) [1/s²] (⟂ ω), ``total`` Gω = (ω·∇)u, ``tilt_rate`` |tilting|/|ω| [1/s].
    Validation: V1 stretching + tilting = Gω, stretching ∥ ω, tilting ⟂ ω; plane flow (ω ∥ e_z, G planar) → both 0;
    V7 invariant under a common rotation of ω and G. Label: analytic.
    """
    w = _F(omega)
    G_ = _F(G)
    mag = float(np.linalg.norm(w))
    total = G_ @ w  # (ω·∇)u = G ω
    if mag == 0.0:
        z = np.zeros(3)
        return dict(e_s=z, rate=0.0, stretching=z, tilting=z, total=total, tilt_rate=0.0)
    es = w / mag
    rate = float(es @ G_ @ es)  # ∂u_s/∂s
    stretching = rate * w  # ω ∂u_s/∂s e_s
    tilting = total - stretching  # ω (∂u_n/∂s e_n + ∂u_m/∂s e_m)
    return dict(e_s=es, rate=rate, stretching=stretching, tilting=tilting, total=total,
                tilt_rate=float(np.linalg.norm(tilting)) / mag)


def planetary_vorticity_terms(G, Omega) -> np.ndarray:
    """Planetary stretching and tilting 2(Ω·∇)u = 2GΩ (3,) [1/s²] (§5.6 after (5.32): with Ω = Ωe_z, Dω_z/Dt =
    2Ω∂w/∂z, Dω_x/Dt = 2Ω∂u/∂z, Dω_y/Dt = 2Ω∂v/∂z — stretching of vertical *fluid* lines, Fig. 5.10).
    Parameters: G (3, 3) [1/s]; Omega (3,) [rad/s]. Book: §5.6. Label: analytic."""
    return 2.0 * _F(G) @ _F(Omega)  # 2(Ω·∇)u


def absolute_circulation(u, pts, Omega, t: float = 0.0) -> tuple[float, float]:
    """Relative and absolute circulation of a loop, Eq. (5.33): Γ_a = ∫_A(ω + 2Ω)·n dA = Γ + 2∫_A Ω·n dA.

    Book: §5.6, Eq. (5.33) (Kelvin in a rotating frame; the book points to "Exercise 5.11", the proof is Exercise
    5.10): DΓ_a/Dt = 0 for inviscid, barotropic, conservatively forced flow seen from a frame rotating at constant Ω.
    For uniform Ω, ∫_A Ω·n dA = Ω·A_vec with A_vec = ½∮x × dx (:func:`loop_vector_area`).

    Parameters: u relative velocity u(x, t) [m/s]; pts (d, N) [m]; Omega (3,) [rad/s] (a scalar = Ω e_z); t [s].
    Returns (Γ, Γ_a) [m²/s].
    Validation: V4 the rotating scenario: Γ changes, Γ_a = 2ΩA₀ constant; V1 fluid at rest in the inertial frame:
    Γ = −2ΩA, Γ_a = 0. Label: analytic, conserved.
    """
    Om = _F(Omega)
    Om3 = np.array([0.0, 0.0, float(Om)]) if Om.ndim == 0 else st.pad3(Om)
    Gm = loop_circulation(u, pts, t)
    return Gm, Gm + 2.0 * float(Om3 @ loop_vector_area(pts))  # Eq. (5.33)


# ======================================================================================================================
# E3 scenarios: the flows of Kelvin's theorem and Helmholtz's first theorem
# ======================================================================================================================
KELVIN_SCENARIOS = ("cellular", "rankine_straddle", "lamb_oseen", "baroclinic", "rotating", "gaussian", "helmholtz",
                    "helmholtz_abc")


def kelvin_scenario(name: str, **p) -> dict:
    """The flows of explainer E3 (``kelvin_material_loop``) and the C03 animation, each with a material loop, the
    fields needed by (5.10) and the circulation Kelvin's theorem predicts.

    Book: §5.2, Eqs. (5.8)–(5.11), the three sources of vorticity and the four restrictions; §5.3 Fig. 5.7
    ("helmholtz"); §5.6, Eq. (5.33) ("rotating"). Scenarios (SI units; keywords override):

    * ``"cellular"`` — steady inviscid cells ψ = Uℓ sin(x/ℓ) sin(y/ℓ) (U = 1 m/s, ℓ = 1 m), ω = 2ψ/ℓ²; circle r = 0.5
      about (π/2, 0.9), 2048 points; t_end = 6 turnover times 6·2π ≈ 37.7 s (the loop grows 3.14 → 21.0 m, spacing
      ≤ 0.021 m); Γ = 1.155130 m²/s constant to ~3e-13 (Kelvin holds). ρ = 1000, p = −ρU²(ψ̂² + |u|²/2U²) (steady Euler,
      ψ̂ = ψ/Uℓ).
    * ``"rankine_straddle"`` — Rankine vortex Γ = 2π m²/s, core σ_c = 1 m; circle r = 0.5 about (0.8, 0) straddling the
      core edge: Γ = 2 × (lens area) = 1.0982124 m²/s, constant (the kink of u at the core edge makes the loop
      quadrature only algebraically convergent: 1024 points by default). ("rankine" is an alias.)
    * ``"lamb_oseen"`` — Γ₀ = 0.01 m²/s, ν = 1e-6 m²/s, σ² = 4ν(t + t₀), t₀ = 10 s, circle r = 5 mm about the axis
      (material: u_r = 0): Γ(t) = Γ₀(1 − e^{−r²/4ν(t+t₀)}) decays — the net viscous force of (5.11); broken "viscous".
    * ``"baroclinic"`` — lock exchange at rest (ρ₁ = 1000 right, ρ₂ = 1025 left, δ = 0.1 m, H = 1 m, hydrostatic p
      with the mean density); square loop of side 0.2 m about (0, H/2): Γ(t) ≈ t·dΓ/dt(0), dΓ/dt = −∮dp/ρ;
      broken "baroclinic".
    * ``"rotating"`` — frame Ω = 0.5 rad/s; axisymmetric convergence α = 0.2 s⁻¹: u = −(α/2)x + (ζ(t)/2)e_z × x with
      ζ(t) = 2Ω(e^{αt} − 1) (at rest at t = 0, w = αz); unit circle: Γ(t) = 2ΩA₀(1 − e^{−αt}) while Γ_a = 2ΩA₀ is
      conserved ((5.33)); broken "inertial".
    * ``"gaussian"`` — steady inviscid Gaussian vortex (Γ = 2π, σ = 1), circle r = 1 about (0.8, 0): Γ constant.
    * ``"helmholtz"`` — inviscid stretched Gaussian vortex tube (3-D; Γ = 1, σ₀ = 1, α = 0.5): a circle of radius 0.3
      drawn on the tube wall R = 1: its vorticity flux is 0 and stays 0 (Helmholtz 1, Fig. 5.7).
    * ``"helmholtz_abc"`` — the ABC flow (3-D steady Euler, ω = u; ``core.vortices.abc_flow_field``) with a circle of
      radius 0.01 m about the origin in the plane spanned by ω(0)/|ω(0)| and (1, −1, 0)/√2 (normal (1, 1, −2)/√6 ⟂ ω):
      a patch on the local tube wall of a non-axisymmetric flow; t_end = 2 s; Γ = 8.016e-10 m²/s constant and the loop's
      vector area stays ⟂ ω (C05's quantitative demo).

    Returns
    -------
    dict with the contract keys ``u``, ``p``, ``rho``, ``grad_p``, ``visc_force``, ``nu``, ``Omega`` (3,), ``pts0``
    (d, N), ``t_end`` [s], ``label``, ``broken`` (None, "viscous", "baroclinic", "inertial"), plus ``Gamma_exact``
    (callable t → Γ predicted), ``Gamma_a_exact`` (rotating) and ``dim``.
    Label: analytic, conserved.
    """
    from .vortices import gaussian_vortex, lamb_oseen_field, rankine_vortex, stretched_gaussian_vortex_field

    n = int(p.get("n") or (1024 if name in ("rankine_straddle", "rankine") else 2048 if name == "cellular" else 256))
    base = dict(name=name, p=None, grad_p=None, visc_force=None, nu=0.0, Omega=np.zeros(3), rho=1000.0,
                Gamma_a_exact=None, dim=2, broken=None)
    if name == "cellular":
        U, ell, rho = float(p.get("U", 1.0)), float(p.get("ell", 1.0)), float(p.get("rho", 1000.0))
        c = p.get("center", (np.pi / 2, 0.9))
        r = float(p.get("radius", 0.5))

        def u(x, t=0.0):
            X = _F(x)
            return U * np.stack([np.sin(X[0] / ell) * np.cos(X[1] / ell), -np.cos(X[0] / ell) * np.sin(X[1] / ell)])

        def pf(x, t=0.0):
            X = _F(x)
            psi = np.sin(X[0] / ell) * np.sin(X[1] / ell)  # ψ/(Uℓ)
            return -rho * U ** 2 * (psi ** 2 + 0.5 * np.sum((u(X) / U) ** 2, axis=0))  # B = −U²ψ̂²: steady Euler
        pts0 = circle_loop_points(c, r, n)
        G0 = loop_circulation(u, pts0)
        return {**base, "u": u, "p": pf, "grad_p": stencil_gradient_fn(pf, 1e-5), "rho": rho, "pts0": pts0,
                "t_end": float(p.get("t_end", 6.0 * 2.0 * np.pi * ell / U)),
                "label": "cellular flow (inviscid): Kelvin holds, Γ constant",
                "Gamma_exact": lambda t: G0 + 0.0 * _F(t)}
    if name in ("rankine_straddle", "rankine", "gaussian"):
        Gam, sig = float(p.get("Gamma", 2 * np.pi)), float(p.get("sigma", 1.0))
        prof = gaussian_vortex if name == "gaussian" else rankine_vortex

        def u(x, t=0.0):
            X = _F(x)
            r_ = np.hypot(X[0], X[1])
            ut = _F(prof(r_, Gam, sig)[0])
            with np.errstate(divide="ignore", invalid="ignore"):
                f = np.where(r_ > 0.0, ut / np.where(r_ > 0.0, r_, 1.0), Gam / (2.0 * np.pi * sig ** 2))
            return np.stack([-f * X[1], f * X[0]])
        pts0 = circle_loop_points(p.get("center", (0.8 * sig, 0.0)),
                                  float(p.get("radius", 0.5 if name != "gaussian" else 1.0)), n)
        G0 = loop_circulation(u, pts0)
        return {**base, "u": u, "pts0": pts0, "t_end": float(p.get("t_end", 20.0)),
                "label": "vortex with the loop across its core: Kelvin holds, Γ constant",
                "Gamma_exact": lambda t: G0 + 0.0 * _F(t)}
    if name == "lamb_oseen":
        Gam, nu = float(p.get("Gamma", 0.01)), float(p.get("nu", 1e-6))
        t0, r0 = float(p.get("t0", 10.0)), float(p.get("radius", 5e-3))
        u = lamb_oseen_field(Gam, nu, t0)
        visc = lambda x, t=0.0: nu * st.laplacian(u, _F(x), t, 1e-4 * r0, (2,))  # noqa: E731  ν∇²u
        return {**base, "u": u, "nu": nu, "visc_force": visc, "pts0": circle_loop_points((0.0, 0.0), r0, n),
                "t_end": float(p.get("t_end", 60.0)), "broken": "viscous",
                "label": "Lamb–Oseen (viscous): the net viscous force on C, (5.11), drains Γ",
                "Gamma_exact": lambda t: Gam * (-np.expm1(-r0 ** 2 / (4.0 * nu * (_F(t) + t0))))}
    if name == "baroclinic":
        r1, r2 = float(p.get("rho1", 1000.0)), float(p.get("rho2", 1025.0))
        dl, H, g = float(p.get("delta", 0.1)), float(p.get("H", 1.0)), float(p.get("g", 9.81))
        side = float(p.get("side", 0.2))
        rb = 0.5 * (r1 + r2)
        rho = lambda x, t=0.0: rb - 0.5 * (r2 - r1) * np.tanh(2.0 * _F(x)[0] / dl)  # noqa: E731  ρ₂ left
        pf = lambda x, t=0.0: rb * g * (H - _F(x)[1]) + 0.0 * _F(x)[0]  # noqa: E731  hydrostatic, mean density
        u = lambda x, t=0.0: 0.0 * _F(x)  # noqa: E731  at rest at the instant the barrier is removed
        pts0 = square_loop_points((0.0, 0.5 * H), side, max(n, 1024))
        gpf = stencil_gradient_fn(pf, 1e-5)
        rate = kelvin_force_terms(pts0, gpf, rho).pressure  # −∮ dp/ρ
        return {**base, "u": u, "p": pf, "grad_p": gpf, "rho": rho, "pts0": pts0, "t_end": float(p.get("t_end", 1.0)),
                "broken": "baroclinic", "rate": rate,
                "label": "lock exchange (baroclinic): −∮dp/ρ ≠ 0 spins the interface up (initial tendency)",
                "Gamma_exact": lambda t: rate * _F(t)}
    if name == "rotating":
        Om, al, r0 = float(p.get("Omega", 0.5)), float(p.get("alpha", 0.2)), float(p.get("radius", 1.0))

        def u(x, t=0.0):
            X = _F(x)
            z_t = 2.0 * Om * np.expm1(al * t)  # ζ(t) = 2Ω(e^{αt} − 1)
            return np.stack([-0.5 * al * X[0] - 0.5 * z_t * X[1], -0.5 * al * X[1] + 0.5 * z_t * X[0]])
        A0 = np.pi * r0 ** 2
        return {**base, "u": u, "Omega": np.array([0.0, 0.0, Om]), "pts0": circle_loop_points((0.0, 0.0), r0, n),
                "t_end": float(p.get("t_end", 10.0)), "broken": "inertial",
                "label": "rotating frame: Γ changes, the absolute circulation Γ_a of (5.33) is conserved",
                "Gamma_exact": lambda t: 2.0 * Om * A0 * (-np.expm1(-al * _F(t))),
                "Gamma_a_exact": lambda t: 2.0 * Om * A0 + 0.0 * _F(t)}
    if name == "helmholtz":
        al, s0, Gam = float(p.get("alpha", 0.5)), float(p.get("sigma0", 1.0)), float(p.get("Gamma", 1.0))
        R0, eps = float(p.get("R0", 1.0)), float(p.get("radius", 0.3))
        u = stretched_gaussian_vortex_field(Gam, s0, al, 0.0)
        th = 2.0 * np.pi * np.arange(n) / n
        phi = eps / R0 * np.cos(th)
        pts0 = np.stack([R0 * np.cos(phi), R0 * np.sin(phi), eps * np.sin(th)])  # a circle drawn on R = R0
        return {**base, "u": u, "pts0": pts0, "dim": 3, "t_end": float(p.get("t_end", 2.0)),
                "label": "Helmholtz 1: a loop on a tube wall encloses no vorticity flux, and never will",
                "Gamma_exact": lambda t: 0.0 * _F(t)}
    if name == "helmholtz_abc":
        from .vortices import abc_flow_field

        u = abc_flow_field(float(p.get("A", 1.0)), float(p.get("B", 1.0)), float(p.get("C", 1.0)))
        w0 = _F(u(np.zeros(3)))
        e1 = w0 / np.linalg.norm(w0)
        e2 = np.array([1.0, -1.0, 0.0]) / np.sqrt(2.0)
        e2 = e2 - (e2 @ e1) * e1
        e2 /= np.linalg.norm(e2)
        th = 2.0 * np.pi * np.arange(n) / n
        pts0 = float(p.get("radius", 0.01)) * (np.cos(th)[None, :] * e1[:, None] + np.sin(th)[None, :] * e2[:, None])
        G0 = loop_circulation(u, pts0)
        return {**base, "u": u, "pts0": pts0, "dim": 3, "t_end": float(p.get("t_end", 2.0)),
                "label": "Helmholtz 1 (ABC flow): a patch on a tube wall keeps its tiny Γ and stays ⟂ ω",
                "Gamma_exact": lambda t: G0 + 0.0 * _F(t)}
    raise ValueError(f"unknown Kelvin scenario {name!r}; choose from {KELVIN_SCENARIOS}")


def kelvin_scenario_circulation(name: str, t: float, n: int | None = None, **p) -> float:
    """Γ(t) of an E3 scenario's material loop [m²/s] (scalar, for parity rows): closed form for "lamb_oseen" and
    "rotating", the prediction rate·t for "baroclinic", otherwise :func:`material_circulation` from t = 0.
    ``n`` loop points (None: the scenario's default, 1024 for the Rankine loop whose velocity has a kink).
    Book: §5.2, (5.8)–(5.11); §5.6, (5.33). Label: analytic, conserved."""
    s = kelvin_scenario(name, n=n, **p)
    if name in ("lamb_oseen", "rotating", "baroclinic") or float(t) == 0.0:
        return float(s["Gamma_exact"](float(t)))
    return float(material_circulation(s["u"], s["pts0"], [0.0, float(t)])[-1])


def kelvin_scenario_rate(name: str, t: float = 0.0, **p) -> dict:
    """The terms of (5.9) and (5.10) for an E3 scenario's material loop at time t [m²/s²] (scalar floats):
    ``acceleration`` ∮(Du/Dt)·dx, ``contour`` ∮u·du (≈ 0), ``pressure`` −∮dp/ρ, ``body`` (0: gravity is conservative),
    ``viscous`` ∮ν∇²u·dx, ``coriolis`` ∮(−2Ω × u)·dx (rotating frame only), ``total`` = acceleration + contour (= DΓ/Dt).
    For "baroclinic" (u ≡ 0 snapshot) the acceleration is the pressure term, the tendency (5.10) predicts at t = 0⁺.
    Book: §5.2, Eqs. (5.9)–(5.11); §5.6 (the Coriolis term that (5.33) absorbs). Label: analytic, converged."""
    s = kelvin_scenario(name, **p)
    P = s["pts0"] if float(t) == 0.0 else material_loop(s["u"], s["pts0"], [0.0, float(t)])[-1]
    scale = float(np.max(np.abs(P))) or 1.0
    kr = kelvin_rate_terms(s["u"], P, float(t), h=1e-4 * scale, ht=1e-4)
    kf = kelvin_force_terms(P, s["grad_p"], s["rho"], None, s["visc_force"], float(t))
    cor = 0.0
    Om = _F(s["Omega"])
    if np.any(Om != 0.0):
        U3 = st.pad3(_F(s["u"](P, float(t))))
        cor = loop_line_integral(-2.0 * np.cross(np.broadcast_to(Om[:, None], U3.shape), U3, axis=0)[:P.shape[0]], P)
    acc = kf.pressure if name == "baroclinic" else kr.acceleration
    return dict(acceleration=acc, contour=kr.contour, pressure=kf.pressure, body=kf.body, viscous=kf.viscous,
                coriolis=cor, total=acc + kr.contour)
