"""Chapter 6 — Ideal flow: the ideal-flow equations and where they hold, ψ and φ as Laplace problems with point
singularities, superposition (half-body, Rankine oval, cylinder with and without circulation), images and Example 6.1's
wall-pressure signal, the complex potential and Cauchy–Riemann, Blasius and Kutta–Zhukhovsky, conformal mapping and
the Zhukhovsky ellipse, finite-difference Laplace with Gauss–Seidel (Example 6.2), axisymmetric flow (Stokes stream
function, sphere, airship, the axial singularity method) and the accelerating sphere with its added mass — plus
re-exports of the reusable primitives the chapter introduced in ``fluidpy.core``.

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e (2012), Ch. 6, §§6.1–6.10, Eqs. (6.1)–(6.109), Figs. 6.1–6.30,
Examples 6.1–6.2, Exercises (test fields and optional targets only). Every equation was transcribed from the rendered
page images (chapters/pages/ch06/p225–p267).

Where the physics lives (all public names are re-exported here, so ``ch06.<name>`` reaches every callable)
--------------------------------------------------------------------------------------------------------
* ``core.potential`` (PF) — plane elements (Uniform, Source, Vortex, Doublet, Corner) and flows (Flow, FunctionFlow)
  with w, dw/dz, u, φ, ψ, C_p, p, stagnation points; half-body, cylinder (Γ_cw/Γ_ccw), images, circle theorem, Blasius
  force, Laurent coefficients; axisymmetric/3-D elements, sphere, the moving sphere (6.96)–(6.108), added mass.
* ``core.conformal`` (CM) — maps, angle preservation, grid images, Zhukhovsky map and its branch-safe inverse,
  mapped flows.
* ``core.laplace_solvers`` (LS) — five-point Laplacian, node update, Jacobi/Gauss–Seidel/SOR sweeps, solver with
  residual history, sparse direct.
* ``core.panels`` (PN) — the axial singularity method, 2-D source panels (our extension).
* this module — the chapter's worked results, checks and explainer state functions.

Conventions (analysis §9, curation §8): **Γ counterclockwise-positive in project code**; body helpers take the book's
clockwise Γ as keyword-only ``Gamma_cw=`` or the project's ``Gamma_ccw=`` (Γ_ccw = −Γ_cw), so L = ρUΓ_cw; 2-D doublet
vector from sink to source; θ from +x (downstream) in the plane formulas; §6.8 z is the horizontal symmetry axis along
the stream; §6.9 ξ = x − x_s; SI units; g = 9.81 m/s² (``G_BOOK``) in this chapter's functions.

Book slips handled (never coded as printed; printed variants kept only as labelled options for wrong-variant tests):
(6.61) 1/z² coefficient printed Ud/π − Γ²/4π² (correct −(Ud/π + Γ²/4π²)) and an extra outer square
(:func:`kutta_zhukhovsky_sym` reports both); (6.104) bracket sign (``core.potential.moving_sphere_surface_velocity``);
(6.108) stray dφ; §6.3 source velocities "from (6.8)" → (6.15); §6.7 "first-order" central differences are second-order
accurate; Example 6.2's FORTRAN loop over I that sets S(6, J) should run over J, Δψ in m²/s (not m²); (6.82) 1/r → 1/r²
(App. B form coded in :func:`spherical_continuity_residual`). Book-quoted numbers stay in the git-ignored
``tests/book_values_ch06.json``.
"""
from __future__ import annotations

from typing import Callable, Sequence

import numpy as np
import sympy as sp
from scipy import integrate
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
from scipy.special import beta as beta_fn

from .core import _stencil as st
from .core import conformal as CM  # noqa: F401
from .core import laplace_solvers as LS  # noqa: F401
from .core import panels as PN  # noqa: F401
from .core import potential as PF  # noqa: F401
from .core._util import as_scalar_if_0d
from .core.conformal import *  # noqa: F401,F403
from .core.laplace_solvers import *  # noqa: F401,F403
from .core.panels import *  # noqa: F401,F403
from .core.potential import *  # noqa: F401,F403
from .core.potential import (AxisymFlow, AxisymUniform, BlasiusForce, ComplexFlow, Doublet, Doublet3D, Flow,
                             FunctionFlow, LineSource3D, PointSource3D, Source, Uniform, Vortex, blasius_force,
                             element_from_spec, gamma_ccw_from, laplacian_residual, mirror, polygon_signed_area)
from .core.thermo import G_BOOK

G = G_BOOK  #: default g [m/s²] of this chapter's functions (the book's 9.81)
_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_C = lambda a: np.asarray(a, dtype=complex)  # noqa: E731
_S = as_scalar_if_0d
_TWO_PI = 2.0 * np.pi
_FOUR_PI = 4.0 * np.pi


# ======================================================================================================================
# §6.1 the ideal-flow equations and where they hold
# ======================================================================================================================
def flow_field_callables(flow, rho: float = 1000.0, p_inf: float = 0.0):
    """A plane :class:`~fluidpy.core.potential.ComplexFlow` as kinematics callables u(x, t) → (2, …) and
    p(x, t) (steady Bernoulli (6.18)) for the residual tools. Label: analytic."""
    def u_fn(X, t=0.0):
        X = _F(X)
        return np.stack([_F(c) for c in flow.velocity(X[0], X[1])])

    def p_fn(X, t=0.0):
        X = _F(X)
        return _F(flow.pressure(X[0], X[1], rho, p_inf))

    return u_fn, p_fn


def _residual_preset(name: str, rho: float, mu: float, **p):
    """(u_fn, p_fn) of the presets of :func:`ideal_flow_residuals`."""
    if name == "cylinder":
        fl = PF.cylinder(p.get("U", 1.0), p.get("a", 0.1))
        return flow_field_callables(fl, rho)
    if name == "poiseuille":
        G_, h0 = p.get("G", 1.0), p.get("h0", 0.01)  # G = −dp/dx [Pa/m], half-gap h0 [m]

        def u(X, t=0.0):
            X = _F(X)
            return np.stack([G_ / (2.0 * mu) * (h0 ** 2 - X[1] ** 2), 0.0 * X[1]])

        return u, (lambda X, t=0.0: -G_ * _F(X)[0])
    if name == "corner":
        A = p.get("A", 1.0)

        def u(X, t=0.0):
            X = _F(X)
            return np.stack([2.0 * A * X[0], -2.0 * A * X[1]])

        return u, (lambda X, t=0.0: -0.5 * rho * (4.0 * A ** 2) * (_F(X)[0] ** 2 + _F(X)[1] ** 2))
    raise ValueError('preset must be "cylinder", "poiseuille" or "corner"')


def ideal_flow_residuals(u_fn, p_fn=None, x=(0.2, 0.15), t: float = 0.0, rho: float = 1000.0, mu: float = 1.0e-3,
                         h: float = 1e-4, ht: float | None = None, **p) -> dict:
    """Residuals of the ideal-flow equations (6.1) ∇·u = 0, ρDu/Dt = −∇p at points x, and the two terms that make
    them hold even for μ ≠ 0: the vorticity and the net viscous force μ∇²u = −μ∇×ω (4.40).

    Parameters
    ----------
    u_fn : velocity field u(x, t) → (d,)/(d, N) [m/s] (ch03 field convention); a plane flow object (then p_fn
        defaults to steady Bernoulli (6.18) with p∞ = 0); or a preset name: "cylinder" (U = 1 m/s, a = 0.1 m, p from
        Bernoulli), "poiseuille" (plane Poiseuille u = (G/2μ)(h₀² − y²), G = −dp/dx = 1 Pa/m, gap 2h₀ = 0.02 m — the
        rotational control), "corner" (u = (2Ax, −2Ay), A = 1 1/s); preset parameters via ``**p``.
    p_fn : pressure p(x, t) [Pa].   x : (d,) or (d, N) points [m] (default (0.2, 0.15));  t [s];  ρ [kg/m³];
    μ [Pa s];  h [m], ht [s] finite-difference steps (2nd order).

    Returns
    -------
    dict(continuity ∇·u [1/s], euler = ρ(∂u/∂t + (u·∇)u) + ∇p [N/m³], viscous_force = μ∇²u [N/m³], vorticity (ω_z in
    2-D, ω in 3-D) [1/s], inertia = |ρ(u·∇)u| [N/m³] (scale for relative checks)).

    Book: §6.1, Eq. (6.1) (with (4.10), (4.38)–(4.40)). Assumptions: constant ρ, p measured from its hydrostatic value.
    Validation (planned): V1 cylinder and corner — all residuals ≈ 0; Poiseuille — viscous force μd²u/dy² = −G, and its
    Euler residual ρDu/Dt + ∇p equals +μ∇²u (the viscous term the ideal equations drop); V2 μ∇²∇φ ≡ 0 for harmonic φ. Label: analytic.
    """
    if isinstance(u_fn, str):
        u_fn, p_def = _residual_preset(u_fn, rho, mu, **p)
        p_fn = p_def if p_fn is None else p_fn
    elif hasattr(u_fn, "velocity") and not callable(u_fn):
        flow = u_fn
        u_fn, p_def = flow_field_callables(flow, rho)
        p_fn = p_def if p_fn is None else p_fn
    X = _F(x)
    d = X.shape[0]
    uv = st.ev(u_fn, X, t, (d,))
    cont = st.div(u_fn, X, t, h)
    Gm = st.grad_vector(u_fn, X, t, h)  # G[i, j] = ∂u_i/∂x_j
    adv = np.einsum("ij...,j...->i...", Gm, uv)
    dudt = st.ddt(u_fn, X, t, ht, (d,))
    gp = st.grad(p_fn, X, t, h) if p_fn is not None else np.zeros_like(uv)
    euler = rho * (dudt + adv) + gp  # Eq. (6.1)
    visc = mu * st.laplacian(u_fn, X, t, h, (d,))  # μ∇²u = −μ∇×ω, (4.40)
    if d == 2:
        vort = Gm[1, 0] - Gm[0, 1]
    else:
        vort = np.stack([Gm[2, 1] - Gm[1, 2], Gm[0, 2] - Gm[2, 0], Gm[1, 0] - Gm[0, 1]])
    return {"continuity": _S(cont), "euler": euler, "viscous_force": visc, "vorticity": _S(vort),
            "inertia": _S(np.sqrt(np.sum((rho * adv) ** 2, axis=0)))}


_REGIONS = {
    "outer": "outside thin attached boundary layers: ideal flow applies",
    "boundary_layer": "inside a boundary layer: vorticity and viscous stress are not negligible",
    "wake": "in a wake: the fluid has passed through boundary layers and carries vorticity",
    "separated": "in a separated region: a rotational bubble that does not thin as Re grows",
    "duct": "an interior (pipe/duct) flow: viscous effects fill the cross-section",
    "pipe": "an interior (pipe/duct) flow: viscous effects fill the cross-section",
    "turbulent": "a turbulent region: rotational by definition",
}


def ideal_flow_applicability(Re: float, M: float = 0.0, baroclinic: bool = False, region: str = "outer",
                             Re_min: float = 1.0e3, M_max: float = 0.3) -> dict:
    """Decision table: may the ideal-flow equations (6.1) be used here? (§6.1)

    Checks (i) Re = ρUL/μ ≥ ``Re_min`` (rule-of-thumb default 10³, a parameter) so vorticity is confined to thin
    boundary layers; (ii) M = U/c ≤ ``M_max`` (0.3, ch04's incompressibility limit) for constant density;
    (iii) no baroclinic density field — else Kelvin's theorem fails (ch05 ``kelvin_hypotheses``) and vorticity is
    created; (iv) the region is the outer flow (not a boundary layer, wake, separated bubble, duct or turbulent region).
    Returns dict(ok, reasons (list of failed checks, each a sentence), verdict (one sentence), kelvin (ch05 verdict)).
    Book: §6.1 (applicability, exclusions and inclusions). Validation (planned): V1 table cases (Re = 50, M = 0.5,
    baroclinic, wake); wrong variant: dropping the Kelvin check passes a baroclinic case. Label: analytic.
    """
    from .ch05_vorticity_dynamics import kelvin_hypotheses_text

    if region not in _REGIONS:
        raise ValueError(f"region must be one of {sorted(_REGIONS)}")
    reasons = []
    if float(Re) < Re_min:
        reasons.append(f"Re = {float(Re):.3g} < {Re_min:.3g}: boundary layers are not thin")
    if float(M) > M_max:
        reasons.append(f"M = {float(M):.3g} > {M_max:.3g}: density changes (compressible flow)")
    kel = kelvin_hypotheses_text(inviscid=True, barotropic=not baroclinic)
    if baroclinic:
        reasons.append("baroclinic density field: " + kel)
    if region != "outer":
        reasons.append(_REGIONS[region])
    ok = not reasons
    verdict = ("Ideal flow applies: " + _REGIONS["outer"] + ".") if ok else ("Ideal flow does not apply: "
                                                                             + "; ".join(r.rstrip(".") for r in reasons) + ".")
    return {"ok": ok, "reasons": reasons, "verdict": verdict, "kelvin": kel}


# ======================================================================================================================
# §6.2 ψ and φ: ω = −∇²ψ, Laplace with point singularities, orthogonality
# ======================================================================================================================
def rankine_vortex_psi(x, y, Gamma: float = 1.0, a: float = 0.1):
    """Stream function of a Rankine vortex (counterclockwise Γ, core radius a) centred at the origin:
    ψ = −Γr²/(4πa²) for r ≤ a (solid-body core), ψ = −(Γ/2π)[ln(r/a) + ½] for r > a (the ideal vortex (6.8)) —
    continuous with a continuous slope at r = a. ω = −∇²ψ = Γ/πa² inside, 0 outside (6.4).
    x, y [m]; Γ [m²/s]; a [m]. Returns ψ [m²/s]. Book: §6.2 (6.4) test field (the Rankine vortex of Ch. 3 §3.5).
    Label: analytic."""
    r2 = _F(x) ** 2 + _F(y) ** 2
    G_, a = float(Gamma), float(a)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where(r2 <= a ** 2, -G_ * r2 / (_FOUR_PI * a ** 2),
                       -G_ / _TWO_PI * (0.5 * np.log(np.maximum(r2, 1e-300) / a ** 2) + 0.5))
    return _S(out)


def vorticity_from_psi(psi_fn, x, y, h: float = 1e-3, **p):
    """Vorticity of a stream-function flow ω_z = ∂v/∂x − ∂u/∂y = −∇²ψ (6.4), by the fourth-order 9-point stencil.
    psi_fn(x, y) [m²/s] or ``"rankine"`` (:func:`rankine_vortex_psi`, keywords ``Gamma``, ``a``); x, y [m]; h [m].
    Returns ω_z [1/s]. Book: §6.2, Eq. (6.4).
    Validation (planned): V1 ψ = sin x sin y (ω = 2 sin x sin y), Rankine core Γ/πa²; V3 order 4.
    Label: analytic, converged."""
    if isinstance(psi_fn, str):
        if psi_fn != "rankine":
            raise ValueError('psi_fn must be a callable or "rankine"')
        G_, a = p.get("Gamma", 1.0), p.get("a", 0.1)
        psi_fn = lambda X, Y: rankine_vortex_psi(X, Y, G_, a)  # noqa: E731
    return _S(-_F(laplacian_residual(psi_fn, x, y, h)))  # Eq. (6.4)


def vorticity_from_psi_sym(psi_expr, x: sp.Symbol, y: sp.Symbol):
    """sympy twin: ω_z = −(ψ_xx + ψ_yy), Eq. (6.4). Label: symbolic."""
    return sp.simplify(-(sp.diff(psi_expr, x, 2) + sp.diff(psi_expr, y, 2)))


def delta_flux_check(fn, center=(0.0, 0.0), radii=(0.01, 0.1, 1.0, 10.0, 100.0), n: int = 256, kind: str = "psi",
                     h: float = 1e-5, **p) -> np.ndarray:
    """Flux of ∇f out of circles |x − x′| = r: ∮∇ψ·n ds = −Γ for a point vortex (6.6), ∮∇φ·n ds = m for a source
    (6.13) — the same for every radius (the δ at the centre is the only "source" of the Laplacian; 0 if the circle does
    not enclose it).

    ``fn``: a plane flow object (∇φ = (u, v), ∇ψ = e_z × ∇φ = (−v, u) exactly), a callable f(x, y) (4th-order
    differences, step h), or ``"vortex"`` (keyword ``Gamma``, default 2π) / ``"source"`` (keyword ``m``, default 2π) at
    the origin. ``kind`` "psi" or "phi"; periodic trapezoid with n nodes (spectral). Returns an array of fluxes [m²/s],
    one per radius. Book: §6.2, Eqs. (6.6), (6.13); Exercises 6.1, 6.3. Validation (planned): V1, V4 radius
    independence. Label: analytic.
    """
    if isinstance(fn, str):
        if fn == "vortex":
            fn = Flow([Vortex(p.get("Gamma", _TWO_PI))])
        elif fn == "source":
            fn = Flow([Source(p.get("m", _TWO_PI))])
        else:
            raise ValueError('fn must be a flow, a callable, "vortex" or "source"')
    c = PF.as_complex_point(center)
    th = np.linspace(0.0, _TWO_PI, int(n), endpoint=False)
    out = []
    for r in np.atleast_1d(_F(radii)):
        z = c + r * np.exp(1j * th)
        x, y = z.real, z.imag
        if hasattr(fn, "velocity"):
            u, v = (_F(q) for q in fn.velocity(x, y))
            gx, gy = (u, v) if kind == "phi" else (-v, u)
        else:
            cc = [1.0, -8.0, 8.0, -1.0]
            ss = [-2.0, -1.0, 1.0, 2.0]
            gx = sum(ci * _F(fn(x + si * h, y)) for ci, si in zip(cc, ss)) / (12 * h)
            gy = sum(ci * _F(fn(x, y + si * h)) for ci, si in zip(cc, ss)) / (12 * h)
        out.append(float(np.sum(gx * np.cos(th) + gy * np.sin(th)) * r * _TWO_PI / int(n)))
    return np.array(out)


def _pts_xy(pts, y=None):
    if y is not None:
        return _F(pts), _F(y)
    if isinstance(pts, tuple) and len(pts) == 2 and not np.iscomplexobj(np.asarray(pts[0])):
        return _F(pts[0]), _F(pts[1])
    a = np.asarray(pts)
    if np.iscomplexobj(a):
        return a.real, a.imag
    a = _F(a)
    return a[0], a[1]


def orthogonality_report(flow, pts, y=None, h: float = 1e-6) -> dict:
    """dict(max_dot_rel = max|∇φ·∇ψ|/|∇φ|², max_ratio_dev = max| |∇φ|/|∇ψ| − 1 |) with both gradients by central
    differences of ``flow.phi`` and ``flow.psi`` (independent of the velocity). Keep points away from branch cuts.
    Book: §6.2 (6.10). Label: analytic."""
    if isinstance(flow, (list, tuple)):
        flow = flow_from_spec(flow)
    x, y = _pts_xy(pts, y)
    px = (_F(flow.phi(x + h, y)) - _F(flow.phi(x - h, y))) / (2 * h)
    py = (_F(flow.phi(x, y + h)) - _F(flow.phi(x, y - h))) / (2 * h)
    sx = (_F(flow.psi(x + h, y)) - _F(flow.psi(x - h, y))) / (2 * h)
    sy = (_F(flow.psi(x, y + h)) - _F(flow.psi(x, y - h))) / (2 * h)
    g2 = px ** 2 + py ** 2
    return {"max_dot_rel": float(np.nanmax(np.abs(px * sx + py * sy) / g2)),
            "max_ratio_dev": float(np.nanmax(np.abs(np.sqrt(g2 / (sx ** 2 + sy ** 2)) - 1.0)))}


def orthogonality_check(flow, pts, y=None, h: float = 1e-6) -> float:
    """Equipotentials cross streamlines at right angles (§6.2 after (6.10); Exercise 6.2): max|∇φ·∇ψ|/|u|² at the
    points, both gradients by central differences of φ and ψ. ``flow`` a Flow or a spec list (:func:`flow_from_spec`);
    ``pts`` complex array, (x, y) tuple of arrays or (2, N) array (or x with ``y``). Returns a float (≈ 1e-10).
    Book: §6.2 (6.10), D04. Label: analytic."""
    return orthogonality_report(flow, pts, y, h)["max_dot_rel"]


def harmonic_polynomials(degree: int) -> list:
    """Harmonic polynomials Re zⁿ and Im zⁿ for n = 1 … degree (sympy, x and y real), each asserted harmonic: a list of
    (Re zⁿ, Im zⁿ) tuples — (φ, ψ) of w = zⁿ. Degree 2 gives (x² − y², 2xy): φ of (6.27) and ψ of (6.24) (A = 1); the
    rotated family (6.25)/(6.26) is (Re, Im) of i z². Book: §6.3 before (6.24); Exercises 6.6–6.7. Label: symbolic."""
    x, y = sp.symbols("x y", real=True)
    out = []
    for n in range(1, int(degree) + 1):
        w = sp.expand((x + sp.I * y) ** n)
        ph, ps = sp.re(w), sp.im(w)
        for f in (ph, ps):
            assert sp.simplify(sp.diff(f, x, 2) + sp.diff(f, y, 2)) == 0
        out.append((ph, ps))
    return out


# ======================================================================================================================
# §6.3 elementary flows: doublet limit, half-body, Rankine oval, superposition state
# ======================================================================================================================
def _pair_flow(m, eps):
    return Flow([Source(m, complex(-eps, 0.0)), Source(-m, complex(eps, 0.0))])


def source_sink_pair(x, y, m: float = 1.0, eps: float = 0.1):
    """Source +m at (−ε, 0) and sink −m at (+ε, 0) (6.28): φ = (m/2π)[ln√((x + ε)² + y²) − ln√((x − ε)² + y²)],
    dipole vector d = −2mε e_x (from sink to source). m [m²/s], ε [m]. Returns (phi, u, v).
    Book: §6.3, Eq. (6.28). Label: analytic."""
    fl = _pair_flow(m, eps)
    u, v = fl.velocity(x, y)
    return fl.phi(x, y), u, v


def doublet_limit_error(eps_list, d: float = 2.0, pts=((1.0, 0.0), (0.0, 1.0), (0.7, 0.7))) -> np.ndarray:
    """Relative error of the source–sink pair with m = d/2ε (2mε = |d| held fixed, dipole −d e_x) against the doublet
    (6.29) φ = |d| cos θ/2πr at the points ``pts``: max|φ_pair − φ_doublet|/max|φ_doublet| → 0 as ε² (the next term of
    the log expansion). eps_list [m]; d [m³/s]. Returns an array, one error per ε. Book: §6.3, (6.28)–(6.29).
    Validation (planned): V3 observed order 2. Label: converged."""
    P = _F(pts).reshape(-1, 2)
    x, y = P[:, 0], P[:, 1]
    ref = _F(Doublet.from_book_scalar(d).phi(x, y))
    out = []
    for e in np.atleast_1d(_F(eps_list)):
        ph = _F(_pair_flow(d / (2.0 * e), e).phi(x, y))
        out.append(float(np.max(np.abs(ph - ref)) / np.max(np.abs(ref))))
    return np.array(out)


def doublet_limit_frames(eps_list=(0.5, 0.3, 0.2, 0.1, 0.05, 0.02), d: float = 2.0, n: int = 161, xlim=(-1.5, 1.5),
                         ylim=(-1.2, 1.2)) -> list:
    """ψ grids of the source–sink pair (m = d/2ε) shrinking to the doublet, for the C04 animation (Figs. 6.5–6.6).
    Returns a list of dict(X, Y, psi, eps, err) (err = :func:`doublet_limit_error` at that ε).
    Book: §6.3 (6.28)–(6.29). Label: analytic."""
    xs = np.linspace(*xlim, int(n))
    ys = np.linspace(*ylim, int(n))
    X, Y = np.meshgrid(xs, ys, indexing="xy")
    return [{"X": X, "Y": Y, "psi": _F(_pair_flow(d / (2.0 * e), e).psi(X, Y)), "eps": float(e),
             "err": float(doublet_limit_error([e], d)[0])} for e in eps_list]


def half_body_numbers(U: float = 1.0, m: float = 2.0 * np.pi) -> dict:
    """Half-body numbers from (6.31): stagnation distance a = m/2πU (x_stag = −a), body streamline ψ = m/2,
    far-downstream half-width h_max = m/2U (mass balance 2Uh_max = m). U [m/s], m [m²/s] → lengths [m], ψ [m²/s].
    Book: §6.3, (6.31), Fig. 6.7. Label: analytic."""
    U, m = float(U), float(m)
    a = m / (_TWO_PI * U)
    return {"a": a, "x_stag": -a, "psi_body": 0.5 * m, "h_max": m / (2.0 * U)}


def half_body_shape(U: float, m: float, theta):
    """Half-body surface from ψ = m/2 in (6.31): r = m(π − θ)/(2πU sin θ), i.e. half-width h = m(π − θ)/2πU.
    θ ∈ (0, 2π) [rad] from +x (lower half θ > π). Returns (x, y) [m]. Book: §6.3. Label: analytic."""
    th = _F(theta)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = float(m) * (np.pi - th) / (_TWO_PI * float(U) * np.sin(th))
        r = np.where(np.isclose(th, np.pi), float(m) / (_TWO_PI * float(U)), r)
    return _S(r * np.cos(th)), _S(r * np.sin(th))


def half_body_surface_cp(theta):
    """Surface pressure coefficient of the half-body — our reduction of (6.32) on ψ = m/2:
    C_p = −(2k cos θ + k²), k = sin θ/(π − θ) (independent of U and m). +1 at the nose (θ → π), 0 where
    tan θ = −2(π − θ) (≈ 113°), → 0⁻ far downstream. θ [rad] from +x. Book: §6.3, Fig. 6.8 (our derivation D08).
    Label: analytic."""
    th = _F(theta)
    with np.errstate(divide="ignore", invalid="ignore"):
        k = np.sin(th) / (np.pi - th)
        cp = -(2.0 * k * np.cos(th) + k ** 2)
    return _S(np.where(np.isclose(th, np.pi), 1.0, cp))


def half_body_cp_zero_angle() -> float:
    """Angle θ [rad] (from +x at the source) where the half-body surface C_p = 0: the root of
    sin θ + 2(π − θ) cos θ = 0 in (π/2, π) (``brentq``). Book: §6.3, Fig. 6.8. Label: analytic."""
    return float(brentq(lambda t: np.sin(t) + 2.0 * (np.pi - t) * np.cos(t), 0.5 * np.pi, np.pi - 1e-3,
                        xtol=1e-15))


def half_body_net_force(U: float = 1.0, m: float = 2.0 * np.pi, x_end: float = 10.0, rho: float = 1.0,
                        n: int = 4001) -> dict:
    """Pressure force per unit depth on the half-body from its nose to the station x = x_end (Exercise 6.13 idea, our
    computation): D = −∫(p − p∞) dy, L = ∫(p − p∞) dx along the body (counterclockwise: upper surface from x_end to the
    nose, then the lower surface back), with p − p∞ = ½ρU²C_p from the surface C_p (our reduction of (6.32)). The
    gauge pressure is used because the contour is open. D → 0 as x_end → ∞ (the whole body feels no force; the tail's
    C_p → 0), L = 0 by symmetry. U [m/s], m [m²/s], x_end [m] > 0, ρ [kg/m³]; n nodes (midpoint rule on the segments).
    Returns dict(D, L) [N/m]. Book: §6.3, N29 (Exercise 6.13). Label: analytic, converged."""
    U, m = float(U), float(m)

    def xb(t):
        return float(half_body_shape(U, m, t)[0])

    t_end = brentq(lambda t: xb(t) - float(x_end), 1e-12, 0.5 * np.pi, xtol=1e-15)
    th = np.linspace(t_end, _TWO_PI - t_end, int(n))
    x, y = (_F(q) for q in half_body_shape(U, m, th))
    tm = 0.5 * (th[:-1] + th[1:])
    pg = 0.5 * rho * U ** 2 * _F(half_body_surface_cp(tm))
    D = -float(np.sum(pg * np.diff(y)))
    L = float(np.sum(pg * np.diff(x)))
    return {"D": D, "L": L}


def rankine_oval(U: float = 1.0, m: float = 2.0 * np.pi, a: float = 1.0) -> dict:
    """Rankine oval (Exercise 6.19): stream U + source m at (−a, 0) + sink −m at (+a, 0) (zero net source ⇒ closed body).
    Stagnation points x = ±L with L² = a² + ma/πU; body ψ = 0; half-width h from h = (m/πU) tan⁻¹(a/h) (``brentq``).
    Returns dict(half_length L, half_width h, psi_fn (ψ(x, y) [m²/s]), stagnation (complex (−L, +L)), flow (inside =
    the oval), x_stag, psi_body 0, closed True). U [m/s], m [m²/s], a [m]. Book: §6.3 (superposition; Exercise 6.19).
    Validation (planned): V1 ψ = 0 on the body, V2 root equation (sympy). Label: analytic."""
    U, m, a = float(U), float(m), float(a)
    base = Flow([Uniform(U), Source(m, complex(-a, 0.0)), Source(-m, complex(a, 0.0))])

    def inside(x, y):
        return _F(base.psi(_F(x), np.abs(_F(y)))) < -1e-12 * m

    L = float(np.sqrt(a ** 2 + m * a / (np.pi * U)))
    h = float(brentq(lambda hh: hh - m / (np.pi * U) * np.arctan(a / hh), 1e-12 * a, m / (2 * U) + a, xtol=1e-15))
    fl = Flow(base.elements, inside=inside, label="Rankine oval")
    return {"half_length": L, "half_width": h, "psi_fn": base.psi, "stagnation": np.array([-L, L], dtype=complex),
            "flow": fl, "x_stag": (-L, L), "psi_body": 0.0, "closed": True}


def flow_from_spec(spec) -> Flow:
    """A :class:`Flow` from a list of plain dicts (the explainers' element kit): {"kind": "uniform", "U", "V"},
    {"kind": "source", "m", "x", "y"}, {"kind": "sink", "m", "x", "y"} (a source of −m), {"kind": "vortex", "Gamma",
    "x", "y"} (counterclockwise), {"kind": "doublet", "dx", "dy", "x", "y"} (dipole vector from sink to source),
    {"kind": "corner", "A", "n"}; tuple specs of ``core.potential.element_from_spec`` and element objects also work.
    Book: §6.2–6.4 (superposition). Label: analytic."""
    if isinstance(spec, ComplexFlow):
        return spec
    els = []
    for e in spec:
        if not isinstance(e, dict):
            els.append(element_from_spec(e))
            continue
        k = e["kind"]
        z0 = complex(e.get("x", 0.0), e.get("y", 0.0))
        if k == "uniform":
            els.append(Uniform(e.get("U", 1.0), e.get("V", 0.0)))
        elif k == "source":
            els.append(Source(e.get("m", 1.0), z0))
        elif k == "sink":
            els.append(Source(-abs(e.get("m", 1.0)), z0))
        elif k == "vortex":
            els.append(Vortex(e.get("Gamma", 1.0), z0))
        elif k == "doublet":
            els.append(Doublet((e.get("dx", -1.0), e.get("dy", 0.0)), z0))
        elif k == "corner":
            els.append(PF.Corner(e.get("A", 1.0), e.get("n", 2.0)))
        else:
            raise ValueError(f"unknown element kind {k!r}")
    return Flow(els)


def superposition_state(spec, x=None, y=None, U_ref: float | None = None, box=(-4.0, 4.0, -3.0, 3.0),
                        probe=None) -> dict:
    """Everything E1's term bars and status need in one call, at the probe (x, y) [m]: the total u, v, ψ, φ, C_p; the
    per-element u_parts, v_parts; the stagnation points (list of [x, y]); ψ of the dividing streamline (through the first
    stagnation point); Σm; ``closed`` (|Σm| < 1e-12 and a body exists: a stream plus sources/sinks/doublets) and a
    status sentence.

    ``spec``: list of dicts (:func:`flow_from_spec`), tuple specs, elements, or a Flow; ``U_ref`` the C_p reference speed
    (default the free-stream speed); ``box`` the stagnation-point search box. (x may also be an (x, y) tuple.)
    Returns dict(u, v, psi, phi, cp, u_parts, v_parts, stagnation, net_source, closed, psi_dividing, status,
    contributions). Book: §6.2–6.3 (superposition, (6.16), (6.31), (6.32)). Label: analytic."""
    fl = flow_from_spec(spec)
    if probe is not None:
        x, y = probe
    elif y is None and x is not None and np.size(x) == 2:
        x, y = x
    px, py = (1.0, 1.0) if x is None else (float(x), float(y))
    contrib = fl.contributions(px, py) if isinstance(fl, Flow) else []
    u, v = fl.velocity(px, py)
    Uref = fl.U_inf if U_ref is None else float(U_ref)
    cp = 1.0 - (u ** 2 + v ** 2) / Uref ** 2 if Uref > 0 else float("nan")
    stag = fl.stagnation_points(box=box)
    psi_div = float(fl.psi(stag[0].real, stag[0].imag)) if len(stag) else float("nan")
    net = fl.net_source if isinstance(fl, Flow) else 0.0
    sing = isinstance(fl, Flow) and any(isinstance(e, (Source, Doublet)) for e in fl.elements)
    body = sing and fl.U_inf > 0
    if body and abs(net) < 1e-12:
        closed, status = True, "closed body: Σm = 0"
    elif body and net > 0:
        closed, status = False, f"open body: net source m = {net:.4g} m²/s ⇒ extends downstream"
    else:
        closed, status = False, "no body"
    return {"u": float(u), "v": float(v), "psi": float(fl.psi(px, py)), "phi": float(fl.phi(px, py)), "cp": float(cp),
            "u_parts": [float(c["u"]) for c in contrib], "v_parts": [float(c["v"]) for c in contrib],
            "stagnation": [[float(q.real), float(q.imag)] for q in stag], "net_source": net, "closed": closed,
            "psi_dividing": psi_div, "status": status, "contributions": contrib}


# ======================================================================================================================
# §6.3 the cylinder: d'Alembert, circulation, stagnation points, lift
# ======================================================================================================================
def cylinder_surface_speed(theta, U: float = 1.0, a: float = 1.0, *, Gamma_cw: float | None = None,
                           Gamma_ccw: float | None = None):
    """Surface velocity u_θ(a, θ) = −2U sin θ − Γ_cw/2πa (6.37) ((6.34) at r = a for Γ = 0) [m/s].
    θ [rad] from +x; keyword-only Γ_cw (book) or Γ_ccw (project) [m²/s]. Book: §6.3. Label: analytic."""
    Gcw = -gamma_ccw_from(Gamma_cw, Gamma_ccw)
    return _S(-2.0 * float(U) * np.sin(_F(theta)) - Gcw / (_TWO_PI * float(a)))  # Eq. (6.37)


def cylinder_surface_cp(theta, U: float = 1.0, a: float = 1.0, *, Gamma_cw: float | None = None,
                        Gamma_ccw: float | None = None):
    """Surface pressure coefficient of the cylinder C_p = 1 − (u_θ/U)², u_θ from (6.37); for Γ = 0,
    C_p = 1 − 4 sin²θ (6.35) (+1 at the stagnation points, −3 at θ = ±π/2). Book: §6.3, (6.35), (6.39).
    Validation (planned): V1 vs (6.35); V1 form cross-check Wikipedia "Potential flow around a circular cylinder".
    Label: analytic."""
    ut = _F(cylinder_surface_speed(theta, U, a, Gamma_cw=Gamma_cw, Gamma_ccw=Gamma_ccw))
    return _S(1.0 - (ut / float(U)) ** 2)  # Eq. (6.35) for Γ = 0


def cylinder_surface_pressure(theta, U: float = 1.0, a: float = 1.0, *, Gamma_cw: float | None = None,
                              Gamma_ccw: float | None = None, rho: float = 1.2, p_inf: float = 0.0):
    """Surface pressure p(a, θ) = p∞ + ½ρ[U² − (−2U sin θ − Γ_cw/2πa)²] (6.39) [Pa]. Book: §6.3. Label: analytic."""
    ut = _F(cylinder_surface_speed(theta, U, a, Gamma_cw=Gamma_cw, Gamma_ccw=Gamma_ccw))
    return _S(p_inf + 0.5 * rho * (float(U) ** 2 - ut ** 2))  # Eq. (6.39)


def cylinder_stagnation_points(U: float = 1.0, a: float = 1.0, *, Gamma_cw: float | None = None,
                               Gamma_ccw: float | None = None, include_inside: bool = False) -> np.ndarray:
    """Stagnation points of the cylinder with circulation (6.38): for |Γ| < 4πaU two surface points with
    sin θ = −Γ_cw/4πaU; at |Γ| = 4πaU one (θ = −π/2 for Γ_cw > 0); for |Γ| > 4πaU one in the flow on the −y axis
    (Γ_cw > 0) at r = [Γ + √(Γ² − (4πaU)²)]/4πU (the other root r₋ = a²/r₊ lies inside; ``include_inside``).
    Returns complex positions [m] (sorted by angle). Book: §6.3, (6.38) and the text after it, Fig. 6.12.
    Validation (planned): V1 |z| = a; Newton cross-check with ``Flow.stagnation_points``; r₊r₋ = a². Label: analytic."""
    U, a = float(U), float(a)
    Gcw = -gamma_ccw_from(Gamma_cw, Gamma_ccw)
    crit = _FOUR_PI * a * U
    sgn = -1.0 if Gcw >= 0 else 1.0  # points move down for clockwise Γ
    if abs(Gcw) < crit * (1.0 - 1e-12):
        s = -Gcw / crit  # Eq. (6.38)
        t1 = float(np.arcsin(s))
        pts = [a * np.exp(1j * t1), a * np.exp(1j * (np.pi - t1))]
    elif abs(Gcw) <= crit * (1.0 + 1e-12):
        pts = [a * np.exp(1j * sgn * 0.5 * np.pi)]
    else:
        G = abs(Gcw)
        rp = (G + np.sqrt(G ** 2 - crit ** 2)) / (_FOUR_PI * U)
        pts = [1j * sgn * rp]
        if include_inside:
            pts.append(1j * sgn * a ** 2 / rp)
    pts = np.array(pts, dtype=complex)
    return pts[np.argsort(np.angle(pts))]


def lift_per_span(rho: float = 1.0, U: float = 1.0, *, Gamma_cw: float | None = None,
                  Gamma_ccw: float | None = None) -> float:
    """Kutta–Zhukhovsky lift per unit span L = ρUΓ (6.40), (6.62) with the book's clockwise Γ (= −ρUΓ_ccw) [N/m].
    ρ [kg/m³], U [m/s]. Book: §6.3 (6.40), §6.5 (6.62). Validation (planned): V1 pressure integral and Blasius;
    V1 form cross-check Wikipedia "Kutta–Joukowski theorem". Label: analytic."""
    return float(rho) * float(U) * (-gamma_ccw_from(Gamma_cw, Gamma_ccw))  # Eq. (6.40)


def _contour_z(contour):
    return PF._contour(contour)


def _p_values(p, z):
    if not callable(p):
        return _F(p)
    import inspect
    try:
        npar = len([q for q in inspect.signature(p).parameters.values()
                    if q.default is q.empty and q.kind in (q.POSITIONAL_ONLY, q.POSITIONAL_OR_KEYWORD)])
    except (TypeError, ValueError):
        npar = 2
    return _F(p(z)) if npar == 1 else _F(p(z.real, z.imag))


def contour_force(p, contour, dcontour=None) -> BlasiusForce:
    """Pressure force per unit depth on a body from its surface pressure (6.56): D = −∮p dy, L = ∮p dx on the
    counterclockwise contour C (outward n = (e_x dy − e_y dx)/ds).

    ``p``: callable p(z) (one complex argument) or p(x, y) [Pa], or values at the vertices; ``contour``: closed polygon
    (complex or (2, N)), last ≠ first, counterclockwise (ValueError otherwise). With ``dcontour`` = dz/dτ at uniform
    τ ∈ [0, 2π) the periodic trapezoid Σ p (dz/dτ)(2π/N) is used (spectral for smooth bodies); else the trapezoid on the
    segments (second order). Returns (D, L) [N/m]. Book: §6.5, Eqs. (6.55)–(6.56).
    Validation (planned): V1 uniform p → 0; p = y → L = −A. Label: analytic."""
    z = _contour_z(contour)
    if polygon_signed_area(z) <= 0:
        raise ValueError("the contour must be counterclockwise (signed area > 0), the orientation of (6.56)")
    pv = _p_values(p, z)
    if dcontour is not None:
        dz = _C(dcontour) * (_TWO_PI / z.size)
        D = -float(np.sum(pv * dz.imag))  # Eq. (6.56): D = −∮ p dy
        L = float(np.sum(pv * dz.real))  # L = ∮ p dx
    else:
        dz = np.roll(z, -1) - z
        pm = 0.5 * (pv + np.roll(pv, -1))
        D = -float(np.sum(pm * dz.imag))
        L = float(np.sum(pm * dz.real))
    return BlasiusForce(D, L)


def complex_force_from_pressure(p, contour, dcontour=None) -> complex:
    """Complex force D − iL = −i∮_C p dz* (6.57) [N/m] (same quadrature as :func:`contour_force`). Book: §6.5.
    Label: analytic."""
    F = contour_force(p, contour, dcontour)
    return complex(F.D, -F.L)


def _circle_tangent(z):
    """If the closed polygon z samples a circle uniformly (counterclockwise), return dz/dτ = i(z − c); else None."""
    c = z.mean()
    r = np.abs(z - c)
    if r.min() <= 0 or np.ptp(r) > 1e-9 * r.mean():
        return None
    ang = np.unwrap(np.angle(z - c))
    d = np.diff(ang)
    if np.ptp(d) > 1e-9 or d.mean() <= 0:
        return None
    return 1j * (z - c)


def surface_pressure_force(flow, contour=None, rho: float = 1.2, p_inf: float = 0.0, U: float | None = None,
                           R: float | None = None, center=0j, n: int = 64) -> BlasiusForce:
    """(D, L) [N/m] on a body by integrating the steady Bernoulli pressure p = p∞ + ½ρ(U² − |u|²) (6.18) round its
    surface (6.56): ``contour`` = the body's points (complex, counterclockwise; a uniformly sampled circle is detected
    and integrated spectrally), or the circle |z − c| = R with n nodes (``R`` keyword, or a plain number as
    ``contour``). For the cylinder: D = 0 (d'Alembert), L = ρUΓ_cw (6.40). Returns (D, L) — index 0 is D.
    Book: §6.3 (Fig. 6.13, (6.39)–(6.40)), §6.5 (6.55)–(6.56). Label: analytic."""
    if contour is not None and np.ndim(contour) == 0 and not isinstance(contour, complex):
        R, contour = float(contour), None
    if contour is None:
        if R is None:
            raise ValueError("give the body contour or the circle radius R")
        th = np.linspace(0.0, _TWO_PI, int(n), endpoint=False)
        e = np.exp(1j * th)
        z = PF.as_complex_point(center) + float(R) * e
        pv = _F(flow.pressure(z.real, z.imag, rho, p_inf, U))
        return contour_force(pv, z, 1j * float(R) * e)
    z = _contour_z(contour)
    pv = _F(flow.pressure(z.real, z.imag, rho, p_inf, U))
    return contour_force(pv, z, _circle_tangent(z))


def cv_force_on_body(flow, R_outer: float, rho: float = 1.2, p_inf: float = 0.0, U: float | None = None,
                     n: int = 512, center=0j) -> BlasiusForce:
    """Force on the body from a control volume between the body and a circle of radius R (6.54) (Exercise 6.27 route):
    F_body = −∮_outer [p n + ρu(u·n)] dl, with p from Bernoulli (6.18) — independent of R.
    Returns (D, L) [N/m]. Book: §6.5, Eq. (6.54) ((4.17), ch04 ``momentum_budget``). Label: analytic."""
    th = np.linspace(0.0, _TWO_PI, int(n), endpoint=False)
    nx, ny = np.cos(th), np.sin(th)
    c = PF.as_complex_point(center)
    x, y = c.real + R_outer * nx, c.imag + R_outer * ny
    u, v = (_F(q) for q in flow.velocity(x, y))
    p = _F(flow.pressure(x, y, rho, p_inf, U))
    un = u * nx + v * ny
    dl = R_outer * _TWO_PI / int(n)
    Fx = -float(np.sum(p * nx + rho * u * un) * dl)  # Eq. (6.54) solved for the force on the body
    Fy = -float(np.sum(p * ny + rho * v * un) * dl)
    return BlasiusForce(Fx, Fy)


def circulation_family_check(U: float = 1.0, a: float = 1.0, Gammas_cw=(0.0, 1.0, 2.0, 4.0), R_far: float | None = None,
                             radii=None, n: int = 256) -> dict:
    """Non-uniqueness round a body (§6.3): every member of (6.36) satisfies the same boundary conditions — zero normal
    velocity on r = a and u → Ue_x far away — while the circulation on any loop enclosing the body is −Γ_cw (the same
    for every radius). ``Gammas_cw`` [m²/s] (the book's clockwise Γ); U [m/s]; a [m]; ``R_far`` the far-field radius
    (default 1000a); ``radii`` of the loops (default 1.5a, 3a, 10a).
    Returns dict(Gammas_cw, max_normal (array, one per Γ) [m/s], far_field (array: max|u − U| on r = R_far) [m/s],
    circulation (array: loop circulation on r = radii[0]) [m²/s], circulation_by_radius (Γ × radii), radii).
    Book: §6.3 (uniqueness, Fig. 6.12). Validation (planned): V1, V4. Label: analytic."""
    from .core.vorticity import loop_circulation

    radii = (1.5 * a, 3.0 * a, 10.0 * a) if radii is None else radii
    R_far = 1000.0 * a if R_far is None else float(R_far)
    th = np.linspace(0.0, _TWO_PI, int(n), endpoint=False)
    body = a * np.exp(1j * th)
    mn, ff, circ = [], [], []
    for G_ in Gammas_cw:
        fl = PF.cylinder(U, a, Gamma_cw=float(G_))
        mn.append(PF.normal_velocity_on(fl, body))
        ff.append(PF.far_field_check(fl, R_far, n=64))

        def ufn(P, t=0.0, fl=fl):
            return np.stack([_F(q) for q in fl.velocity(P[0], P[1])])
        circ.append([float(loop_circulation(ufn, np.stack([r * np.cos(th), r * np.sin(th)]))) for r in radii])
    circ = np.array(circ)
    return {"Gammas_cw": np.array(Gammas_cw, dtype=float), "radii": list(radii), "max_normal": np.array(mn),
            "far_field": np.array(ff), "circulation": circ[:, 0], "circulation_by_radius": circ}


def cylinder_circulation_state(U: float = 10.0, a: float = 0.1, *, Gamma_cw: float | None = None,
                               Gamma_ccw: float | None = None, rho: float = 1.2, p_inf: float = 0.0,
                               n: int = 128) -> dict:
    """E2's status, force bars and the C07 slider in one call (Γ_cw defaults to 2 m²/s if neither Γ is given):
    stagnation angles θ₁ = arcsin(−Γ_cw/4πaU), θ₂ = 180° − θ₁ (wrapped to (−180°, 180°]; NaN when the stagnation point
    is off the body), r_free (the off-body point's distance [m], NaN otherwise), the regime ("two surface points",
    "merged at the bottom" (or "... top" for Γ_cw < 0), "free stagnation point"), D and L by the surface-pressure
    quadrature (n nodes), L_KJ = ρUΓ_cw (6.40), the surface speeds at the top and bottom (6.37), min C_p.
    Returns dict(theta1_deg, theta2_deg, r_free, regime, D, L, L_KJ, speed_top, speed_bottom, cp_min, stagnation,
    stagnation_angles_deg, ratio, rhoUGamma). Book: §6.3 (6.36)–(6.40), Fig. 6.12. Label: analytic."""
    if Gamma_cw is None and Gamma_ccw is None:
        Gamma_cw = 2.0
    Gcw = -gamma_ccw_from(Gamma_cw, Gamma_ccw)
    fl = PF.cylinder(U, a, Gamma_cw=Gcw)
    stag = cylinder_stagnation_points(U, a, Gamma_cw=Gcw)
    ratio = Gcw / (_FOUR_PI * a * U)
    t1 = t2 = rfree = float("nan")
    if abs(abs(ratio) - 1.0) <= 1e-9:
        regime = "merged at the bottom" if Gcw > 0 else "merged at the top"
        t1 = t2 = -90.0 if Gcw > 0 else 90.0
    elif abs(ratio) < 1.0:
        regime = "two surface points"
        t1 = float(np.degrees(np.arcsin(-ratio)))
        t2 = 180.0 - t1
        t2 = t2 - 360.0 if t2 > 180.0 else t2
    else:
        regime = "free stagnation point"
        rfree = float(abs(stag[0]))
    F = surface_pressure_force(fl, R=a, rho=rho, p_inf=p_inf, n=n)
    th = np.linspace(0.0, _TWO_PI, 721)
    return {"theta1_deg": t1, "theta2_deg": t2, "r_free": rfree, "regime": regime, "D": F.D, "L": F.L,
            "L_KJ": lift_per_span(rho, U, Gamma_cw=Gcw),
            "speed_top": float(abs(cylinder_surface_speed(0.5 * np.pi, U, a, Gamma_cw=Gcw))),
            "speed_bottom": float(abs(cylinder_surface_speed(-0.5 * np.pi, U, a, Gamma_cw=Gcw))),
            "cp_min": float(np.min(cylinder_surface_cp(th, U, a, Gamma_cw=Gcw))),
            "stagnation": stag, "stagnation_angles_deg": np.degrees(np.angle(stag)), "ratio": float(ratio),
            "rhoUGamma": lift_per_span(rho, U, Gamma_cw=Gcw)}


def force_on_held_singularity(kind: str = "source", U: float = 1.0, strength: float = 1.0,
                              rho: float = 1.0) -> dict:
    """Force per unit depth on a source (m) or a counterclockwise vortex (Γ) held fixed in a stream U (Exercise 6.10),
    by Blasius (6.60) on a circle round it: source D = −ρmU (pushed upstream), vortex L = −ρUΓ_ccw.
    Returns dict(D, L). Book: §6.5 (6.60); Exercise 6.10 (our derivation). Label: analytic."""
    el = Source(strength) if kind == "source" else Vortex(strength) if kind == "vortex" else None
    if el is None:
        raise ValueError('kind must be "source" or "vortex"')
    F = blasius_force(Flow([Uniform(U), el]), R=1.0, rho=rho, n=64)
    return {"D": F.D, "L": F.L}


# ======================================================================================================================
# §6.3 images and Example 6.1
# ======================================================================================================================
def two_sources(m: float = 1.0, a: float = 1.0) -> Flow:
    """Source m at (a, 0) and its image m at (−a, 0) (6.41), (6.53): source near the wall x = 0, or slit flow into a
    right-angled corner. Built with :func:`~fluidpy.core.potential.mirror`. Returns a :class:`Flow`.
    Book: §6.3 (6.41), Fig. 6.16; §6.4 (6.53). Label: analytic."""
    return Flow(mirror([Source(m, complex(a, 0.0))], "x=0"), label="two sources")


def two_source_streamline(psi: float, m: float, a: float, x, both: bool = False):
    """The streamline ψ of (6.41) as y(x): from x² − y² − 2xy cot(2πψ/m) = a², the quadratic
    y² + 2x cot(c) y + (a² − x²) = 0 with c = 2πψ/m: y = −x cot c ± √(x² cot²c + x² − a²) (NaN where no real root).
    Returns y₊ (or (y₊, y₋) with ``both``) [m]. Book: §6.3, the equation after (6.41). Label: analytic."""
    x = _F(x)
    cot = 1.0 / np.tan(_TWO_PI * float(psi) / float(m))
    disc = x ** 2 * cot ** 2 + x ** 2 - float(a) ** 2
    with np.errstate(invalid="ignore"):
        rt = np.sqrt(disc)
    yp, ym = -x * cot + rt, -x * cot - rt
    return (_S(yp), _S(ym)) if both else _S(yp)


def example_6_1_times(Gamma: float = 1.0, h: float = 1.0, rho: float = 1000.0) -> dict:
    """Key instants of Example 6.1's wall-pressure signal (our reading of the book's closed form): strongest suction
    p − p∞ = −ρΓ²/4π²h² at t = 0; zero crossing at t = 4πh²/Γ (η = h); largest over-pressure ρΓ²/32π²h² at
    t = 4√3πh²/Γ (η = √3h). Γ [m²/s] (the book's: vortex of strength −Γ), h [m]. Label: analytic."""
    G, h = float(Gamma), float(h)
    return {"t_zero": _FOUR_PI * h ** 2 / G, "t_max": _FOUR_PI * np.sqrt(3.0) * h ** 2 / G,
            "p_min": -rho * G ** 2 / (_FOUR_PI * np.pi * h ** 2), "p_max": rho * G ** 2 / (32.0 * np.pi ** 2 * h ** 2)}


def example_6_1(t, Gamma: float = 1.0, h: float = 1.0, rho: float = 1000.0, p_inf: float = 0.0,
                route: str = "closed", dt: float | None = None) -> dict:
    """Example 6.1: a free vortex of strength −Γ (clockwise Γ > 0) starting at (h, 0) beside the wall x = 0 in still
    fluid: trajectory ξ(t) = (h, Γt/4πh) and the pressure at the wall origin
    (p(0, 0, t) − p∞)/ρ = (Γ²/4π²)(η² − h²)/(η² + h²)², η = Γt/4πh — suction, then over-pressure.

    ``route="closed"``: the book's closed forms (∂φ/∂t = −Γ²/4π²(h² + η²), v(0, 0, t) = Γh/π(h² + η²));
    ``route="numeric"``: the path from ``core.biot_savart.point_vortex_evolve(boundary="wall")`` (rotated so the wall is
    y′ = 0), the velocity from the vortex + image :class:`Flow`, ∂φ/∂t by a central difference in t (step dt) and the
    unsteady Bernoulli equation (4.83) with g = 0.
    Parameters: t [s] (scalar or array); Γ [m²/s]; h [m]; ρ [kg/m³]; p∞ [Pa].
    Returns dict(xi (2,)/(2, T), xi_x, xi_y [m], p_origin [Pa], v_origin [m/s], dphidt [m²/s²], unsteady = −ρ∂φ/∂t
    [Pa], speed_part = −½ρv² [Pa] (p_origin = p∞ + unsteady + speed_part), route).
    Book: §6.3, Example 6.1 (Fig. 5.14). Assumptions: self-induced velocity of the ideal vortex taken as zero.
    Validation (planned): V1 closed vs numeric (O(dt²)); sign change at 4πh²/Γ; V4 ξ_x conserved. Label: analytic.
    """
    G, h = float(Gamma), float(h)
    tt = np.atleast_1d(_F(t))
    if route == "closed":
        eta = G * tt / (_FOUR_PI * h)
        s = h ** 2 + eta ** 2
        dphidt = -G ** 2 / (_FOUR_PI * np.pi * s)
        v0 = G * h / (np.pi * s)
        p = p_inf + rho * G ** 2 / (_FOUR_PI * np.pi) * (eta ** 2 - h ** 2) / s ** 2
        xi = np.stack([np.full(tt.shape, h), eta])
    elif route == "numeric":
        from .core.biot_savart import point_vortex_evolve

        # rotate +90°: (x, y) → (x′, y′) = (−y, x); the wall x = 0 becomes y′ = 0, fluid y′ > 0; Γ_ccw = −Γ
        te = np.concatenate([[0.0], tt]) if tt[0] != 0.0 else tt
        tr = point_vortex_evolve(np.array([[0.0], [h]]), [-G], te, boundary="wall", wall_y=0.0)
        tr = tr[1:] if tt[0] != 0.0 else tr
        xi = np.stack([tr[:, 1, 0], -tr[:, 0, 0]])  # back: x = y′, y = −x′
        dt_ = (1e-4 * _FOUR_PI * h ** 2 / G) if dt is None else float(dt)

        def fl_at(tq):
            eta_q = G * tq / (_FOUR_PI * h)  # the path is ξ = (h, η): checked against xi above
            return Flow([Vortex(-G, complex(h, eta_q), cut_angle=0.0), Vortex(G, complex(-h, eta_q))])

        dphidt = np.array([(float(fl_at(tq + dt_).phi(0.0, 0.0)) - float(fl_at(tq - dt_).phi(0.0, 0.0))) / (2 * dt_)
                           for tq in tt])
        uv = np.array([fl_at(tq).velocity(0.0, 0.0) for tq in tt])
        v0 = uv[:, 1]
        p = p_inf - rho * (dphidt + 0.5 * (uv[:, 0] ** 2 + uv[:, 1] ** 2))  # (4.83) with g = 0
    else:
        raise ValueError('route must be "closed" or "numeric"')
    sq = np.ndim(t) == 0
    unst = -rho * dphidt
    spd = -0.5 * rho * np.asarray(v0) ** 2

    def one(q):
        return _S(np.asarray(q)[0]) if sq else np.asarray(q)
    return {"xi": xi[:, 0] if sq else xi, "xi_x": one(xi[0]), "xi_y": one(xi[1]), "p_origin": one(p),
            "v_origin": one(v0), "dphidt": one(dphidt), "unsteady": one(unst), "speed_part": one(spd), "route": route}


def example_6_1_wall_pressure(y_wall, t: float, Gamma: float = 1.0, h: float = 1.0, rho: float = 1000.0,
                              p_inf: float = 0.0, split: bool = False):
    """Pressure along the whole wall x = 0 in Example 6.1 (our extension of the book's origin value): the unsteady part
    −ρ∂φ/∂t = ρΓ²/4π²s² and the speed part −½ρv² = −ρΓ²h²/2π²s⁴, s² = h² + (y − η)², η = Γt/4πh
    (φ(0, y, t) = (Γ/π) tan⁻¹((y − η)/h), u = 0 and v = Γh/πs² on the wall).
    y_wall [m] (scalar or array), t [s]. Returns p [Pa] (= p∞ + both parts), or with ``split=True`` dict(unsteady,
    speed, total, p (= total), speed_part, dphidt, v, eta). Book: §6.3, Example 6.1.
    Validation (planned): V1 y = 0 equals :func:`example_6_1`. Label: analytic."""
    G_, h = float(Gamma), float(h)
    eta = G_ * float(t) / (_FOUR_PI * h)
    s2 = h ** 2 + (_F(y_wall) - eta) ** 2
    dphidt = -G_ ** 2 / (_FOUR_PI * np.pi * s2)
    v = G_ * h / (np.pi * s2)
    unst = -rho * dphidt
    spd = -0.5 * rho * v ** 2
    total = p_inf + unst + spd
    if not split:
        return _S(total)
    return {"unsteady": _S(unst), "speed": _S(spd), "total": _S(total), "p": _S(total), "speed_part": _S(spd),
            "dphidt": _S(dphidt), "v": _S(v), "eta": eta}


# ======================================================================================================================
# §6.4 complex potential: Cauchy–Riemann, corners
# ======================================================================================================================
def complex_potential_family(name: str, **p) -> tuple[Callable, Callable | None]:
    """(w, dw/dz) callables of E4's families: "corner" (A=1, n=2), "uniform" (U=1), "source" (m=2π), "vortex"
    (Gamma=2π, counterclockwise), "doublet" (d=2π, (6.49)'s scalar), "cylinder" (U=1, a=1, Gamma_cw=0) and the
    non-analytic control "conj" (f = z*, dw/dz None). Book: §6.4 (6.46)–(6.52). Label: analytic."""
    if name == "corner":
        e = PF.Corner(p.get("A", 1.0), p.get("n", 2.0))
        return e.w, e.dwdz
    if name == "uniform":
        e = Uniform(p.get("U", 1.0), p.get("V", 0.0))
        return e.w, e.dwdz
    if name == "source":
        e = Source(p.get("m", _TWO_PI))
        return e.w, e.dwdz
    if name == "vortex":
        e = Vortex(p.get("Gamma", _TWO_PI))
        return e.w, e.dwdz
    if name == "doublet":
        e = Doublet.from_book_scalar(p.get("d", _TWO_PI))
        return e.w, e.dwdz
    if name == "cylinder":
        f = PF.cylinder(p.get("U", 1.0), p.get("a", 1.0), Gamma_cw=p.get("Gamma_cw", 0.0))
        return f.w, f.dwdz
    if name == "conj":
        return (lambda z: np.conj(_C(z))), None
    raise ValueError(f"unknown family {name!r}")


def _cr_quotients(w_fn, z, h):
    def W(q):
        return complex(np.asarray(w_fn(np.array([q]))).ravel()[0])
    Dx = (W(z + h) - W(z - h)) / (2 * h)
    Dy = (W(z + 1j * h) - W(z - 1j * h)) / (2j * h)
    return Dx, Dy


def cauchy_riemann_residual(w_fn, z, h: float = 1e-6, **p) -> tuple:
    """Cauchy–Riemann residuals (6.44) from difference quotients of w along x and along iy at z:
    D_x = [w(z + h) − w(z − h)]/2h = φ_x + iψ_x, D_y = [w(z + ih) − w(z − ih)]/2ih = ψ_y − iφ_y;
    r₁ = ∂φ/∂x − ∂ψ/∂y, r₂ = ∂φ/∂y + ∂ψ/∂x (both ≈ 0 iff w is analytic, then D_x = D_y = dw/dz = u − iv (6.45)).
    ``w_fn`` a callable or a :func:`complex_potential_family` name (e.g. "conj", the non-analytic control z*, which gives
    r₁ = 2). z complex [m]; h [m]. Returns (r1, r2) [1/s]. Book: §6.4, Eqs. (6.44)–(6.45). Label: analytic."""
    if isinstance(w_fn, str):
        w_fn = complex_potential_family(w_fn, **p)[0]
    Dx, Dy = _cr_quotients(w_fn, complex(z), h)
    return float(Dx.real - Dy.real), float(Dx.imag - Dy.imag)


def cauchy_riemann_residual_sym(expr, x: sp.Symbol, y: sp.Symbol | None = None) -> tuple:
    """sympy Cauchy–Riemann residuals (φ_x − ψ_y, φ_y + ψ_x) of w = φ + iψ, Eq. (6.44), simplified ((0, 0) for analytic
    w). Call as ``(expr, x, y)`` with ``expr`` written in real symbols x, y (e.g. (x + I*y)**2), or ``(expr, z)`` with
    ``expr`` a function of the complex symbol z (substituted by x + iy). Label: symbolic."""
    if y is None:
        z = x
        x, y = sp.symbols("x y", real=True)
        expr = expr.subs(z, x + sp.I * y)
    wx = sp.expand_complex(expr)
    ph, ps = sp.re(wx), sp.im(wx)
    return (sp.simplify(sp.diff(ph, x) - sp.diff(ps, y)), sp.simplify(sp.diff(ph, y) + sp.diff(ps, x)))


def corner_info(n: float) -> dict:
    """Corner flow w = Azⁿ (6.46): dict(alpha [rad], alpha_deg, exponent = n − 1, regime ("stagnation point" for α < π,
    "uniform flow", "infinite speed at the corner" for α > π)). Book: §6.4 after (6.46). Label: analytic."""
    n = float(n)
    if n < 0.5:
        raise ValueError("n must be ≥ 1/2")
    al = np.pi / n
    regime = "stagnation point" if n > 1 else ("uniform flow" if n == 1 else "infinite speed at the corner")
    return {"alpha": float(al), "alpha_deg": float(np.degrees(al)), "exponent": n - 1.0, "regime": regime}


def corner_speed_exponent(n: float) -> float:
    """Speed exponent of the corner flow w = Azⁿ (6.46): |dw/dz| ∝ r^{n−1} = r^{(π−α)/α}, α = π/n (> 0: the corner is a
    stagnation point; < 0: infinite speed). Returns n − 1. Book: §6.4 after (6.46). Label: analytic."""
    return corner_info(n)["exponent"]


def complex_potential_probe(kind: str, x: float, y: float, A: float = 1.0, n: float = 2.0, h: float = 1e-4,
                            **p) -> dict:
    """E4's probe: the complex potential of family ``kind`` ("corner" (A, n), "source", "vortex", "doublet", "cylinder",
    or the non-analytic control "conj" = z*; parameters as in :func:`complex_potential_family`) at z = x + iy:
    φ, ψ (6.42); the exact dw/dz and (u, v) = (Re, −Im) (6.45) (NaN for "conj"); the difference quotients along x and
    along iy with step h, q_x = [w(z + h) − w(z − h)]/2h and q_y = [w(z + ih) − w(z − ih)]/2ih; the Cauchy–Riemann
    residuals cr1 = Re q_x − Re q_y (= φ_x − ψ_y), cr2 = Im q_x − Im q_y (= ψ_x + φ_y) (6.44); speed |dw/dz|.
    x, y [m]; h [m]. Returns dict(phi, psi, dwdz_re, dwdz_im, u, v, qx_re, qx_im, qy_re, qy_im, cr1, cr2, speed).
    Book: §6.4 (6.42)–(6.46), D14. Label: analytic."""
    pp = dict(p)
    if kind == "corner":
        pp.update(A=A, n=n)
    wf, dwf = complex_potential_family(kind, **pp)
    z = complex(x, y)
    w = complex(np.asarray(wf(np.array([z]))).ravel()[0])
    qx, qy = _cr_quotients(wf, z, h)
    if dwf is None:
        dw = complex(np.nan, np.nan)
    else:
        dw = complex(np.asarray(dwf(np.array([z]))).ravel()[0])
    return {"phi": w.real, "psi": w.imag, "dwdz_re": dw.real, "dwdz_im": dw.imag, "u": dw.real, "v": -dw.imag,
            "qx_re": qx.real, "qx_im": qx.imag, "qy_re": qy.real, "qy_im": qy.imag,
            "cr1": float(qx.real - qy.real), "cr2": float(qx.imag - qy.imag), "speed": float(abs(dw))}


# ======================================================================================================================
# §6.5 Blasius and Kutta–Zhukhovsky
# ======================================================================================================================
def kutta_zhukhovsky_sym() -> dict:
    """sympy: square the far-field complex velocity U + iΓ/2πz − d/2πz² (clockwise Γ, book scalar d), pick the residue
    (coefficient of 1/z) and apply Blasius (6.60): D − iL = (iρ/2)·2πi·(iUΓ/π) = −iρUΓ (6.62).

    Returns dict(square (series), residue = coeff_z1 = iUΓ/π, coeff_z2 (correct: −(Ud/π + Γ²/4π²)), coeff_z2_printed (the book's
    (6.61) prints Ud/π − Γ²/4π² — a slip; it does not affect the result), DmiL, D, L). Book: §6.5 (6.61)–(6.62).
    Label: symbolic."""
    U, Gam, d, rho = sp.symbols("U Gamma d rho", positive=True)
    z = sp.symbols("z")
    f = U + sp.I * Gam / (2 * sp.pi * z) - d / (2 * sp.pi * z ** 2)
    sq = sp.expand(f ** 2)
    res = sp.simplify(sq.coeff(z, -1))
    c2 = sp.simplify(sq.coeff(z, -2))
    DmiL = sp.simplify(sp.I * rho / 2 * 2 * sp.pi * sp.I * res)  # Eq. (6.62)
    return {"square": sq, "residue": res, "coeff_z1": res, "coeff_z2": c2,
            "coeff_z2_printed": U * d / sp.pi - Gam ** 2 / (4 * sp.pi ** 2),
            "DmiL": DmiL, "D": sp.re(DmiL), "L": sp.simplify(-sp.im(DmiL))}


def _body_flow(body: str, Gamma_cw: float, U: float, a: float, b: float | None, alpha: float):
    """(flow, inside, pressure-route contour (z, dz/dτ) or None) of E5's bodies."""
    if body == "cylinder":
        fl = PF.cylinder(U, a, Gamma_cw=Gamma_cw)
        return fl, fl.inside, None
    if body in ("ellipse", "tilted_ellipse"):
        bb = a / 1.2 if b is None else float(b)
        al = alpha if (body == "ellipse" or alpha != 0.0) else np.radians(15.0)
        fl = elliptic_cylinder_flow(U, a, bb, Gamma_cw=Gamma_cw, alpha=al)
        fl._b = bb
        return fl, fl.inside, None
    if body == "rankine_oval_vortex":
        s_, m_ = 0.5 * a, _TWO_PI * U * a
        fl = Flow([Uniform(U), Source(m_, complex(-s_, 0.0)), Source(-m_, complex(s_, 0.0)), Vortex(-Gamma_cw)])
        ro = rankine_oval(U, m_, s_)
        Lh, Hh = ro["half_length"], ro["half_width"]

        def inside(x, y):  # approximate body: the ellipse through the oval's (no-vortex) axes
            return (_F(x) / Lh) ** 2 + (_F(y) / Hh) ** 2 < 1.0
        return fl, inside, None
    raise ValueError('body must be "cylinder", "ellipse", "tilted_ellipse" or "rankine_oval_vortex"')


def laurent_contributions(body="cylinder", R: float = 0.2, *, Gamma_cw: float = 2.0, U: float = 10.0, a: float = 0.1,
                          b: float | None = None, kmax: int = 4, rho: float = 1.2, n: int = 128, alpha: float = 0.0,
                          center=0j) -> dict:
    """E5's bars: the share of each power z^k of (dw/dz)² (k = 0, −1, …, −kmax) in Blasius's (iρ/2)∮(dw/dz)² dz (6.60):
    ∮z^k dz = 2πi only for k = −1, so only that bar is nonzero, and its value is D − iL = −iρUΓ (6.61)–(6.62).

    ``body``: a name ("cylinder", "ellipse" (Zhukhovsky, b default a/1.2), "tilted_ellipse", "rankine_oval_vortex"; built
    with Γ_cw [m²/s], U [m/s], a [m], b [m], alpha [rad]) — or any flow object / dw/dz callable. R [m] the circle for the
    Laurent coefficients (FFT, n samples). Returns dict(powers [0, −1, …, −kmax], contrib_re (drag part D of each),
    contrib_im (lift part L of each) [N/m], coeffs {k: c_k of dw/dz}, square {p: coefficient of z^p in (dw/dz)²}, pairs
    ((j, k, 2πi c_j c_k × multiplicity) with j + k = −1), integral (∮(dw/dz)²dz), D, L).
    Book: §6.5 (6.61)–(6.62). Label: analytic."""
    fl = _body_flow(body, Gamma_cw, U, a, b, alpha)[0] if isinstance(body, str) else body
    kk = int(kmax)
    c = PF.laurent_coefficients(fl, R, n=n, kmin=-(kk + 2), kmax=2, center=center)
    ks = sorted(c)
    sq: dict = {}
    for j in ks:
        for k in ks:
            sq[j + k] = sq.get(j + k, 0j) + c[j] * c[k]
    pairs = []
    for j in ks:
        k = -1 - j
        if k in c and j <= k:
            mult = 1 if j == k else 2
            pairs.append((j, k, complex(2j * np.pi * c[j] * c[k] * mult)))
    powers = [0] + [-q for q in range(1, kk + 1)]
    cre, cim = [], []
    for pw in powers:
        part = 0.5j * rho * (2j * np.pi * sq.get(-1, 0j)) if pw == -1 else 0j  # (iρ/2)∮ c z^k dz
        cre.append(float(part.real))
        cim.append(float(-part.imag))
    integ = complex(2j * np.pi * sq.get(-1, 0j))
    DmiL = 0.5j * rho * integ  # Eq. (6.60)
    return {"powers": powers, "contrib_re": np.array(cre), "contrib_im": np.array(cim), "coeffs": c, "square": sq,
            "pairs": pairs, "integral": integ, "D": float(DmiL.real), "L": float(-DmiL.imag)}


def blasius_state(body: str = "cylinder", *, Gamma_cw: float = 2.0, U: float = 10.0, a: float = 0.1, R: float = 0.2,
                  rho: float = 1.2, n: int = 256, b: float | None = None, alpha: float = 0.0, offset_x: float = 0.0,
                  offset_y: float = 0.0) -> dict:
    """E5's state: Blasius (6.60) on the circle of radius R centred at (offset_x, offset_y) round the chosen body, with
    the Laurent coefficients there, the validity flag, and the independent on-body pressure route (6.56).

    Bodies: "cylinder" (radius a, (6.52)); "ellipse" (Zhukhovsky image of the circle a, map constant b (default a/1.2),
    (6.68)–(6.69)); "tilted_ellipse" (the same with the stream at angle ``alpha``, default 15° if 0); "rankine_oval_vortex"
    (stream + source 2πUa at −a/2 + sink at +a/2 + clockwise vortex Γ; its body test is approximate — the ellipse through
    the no-vortex oval's axes — and it has no pressure route (NaN)). Γ_cw [m²/s] clockwise (book); U [m/s]; a, R [m];
    ρ [kg/m³]; n quadrature nodes.
    Returns dict(D, L (Blasius on the contour) [N/m], L_KJ = ρUΓ_cw, crosses_body (a contour point inside the body:
    Blasius invalid), c0_re, cm1_re, cm1_im (Laurent coefficients c₀ and c₋₁ of dw/dz on the contour: U and iΓ_cw/2π),
    D_pressure, L_pressure (surface-pressure integral on the body, spectral, n nodes)).
    Book: §6.5 (6.56), (6.60)–(6.62). Label: analytic."""
    fl, inside, _ = _body_flow(body, Gamma_cw, U, a, b, alpha)
    c = complex(offset_x, offset_y)
    F = blasius_force(fl, R=R, rho=rho, n=n, center=c)
    th = np.linspace(0.0, _TWO_PI, int(n), endpoint=False)
    cz = c + R * np.exp(1j * th)
    crosses = bool(np.any(np.asarray(inside(cz.real, cz.imag)))) if inside is not None else False
    lc = PF.laurent_coefficients(fl, R, n=max(int(n), 8), kmin=-2, kmax=0, center=c)
    Dp = Lp = float("nan")
    if body == "cylinder":
        Fp = surface_pressure_force(fl, R=a, rho=rho, n=n)
        Dp, Lp = Fp.D, Fp.L
    elif body in ("ellipse", "tilted_ellipse"):
        bb = fl._b
        zeta = a * np.exp(1j * th)
        zb = _C(CM.joukowski(zeta, bb))
        dz = 1j * (zeta - bb ** 2 / zeta)  # dz/dθ on the body
        al = alpha if (body == "ellipse" or alpha != 0.0) else np.radians(15.0)
        _, dwz = CM.circle_flow_zeta(U, a, Gamma_cw=Gamma_cw, alpha=al)
        q = _C(dwz(zeta)) / (1.0 - bb ** 2 / zeta ** 2)  # u − iv on the body (chain rule, no inside mask)
        pv = 0.5 * rho * (U ** 2 - np.abs(q) ** 2)
        Fp = contour_force(pv, zb, dz)
        Dp, Lp = Fp.D, Fp.L
    return {"D": F.D, "L": F.L, "L_KJ": float(rho * U * Gamma_cw), "crosses_body": crosses,
            "c0_re": float(lc[0].real), "cm1_re": float(lc[-1].real), "cm1_im": float(lc[-1].imag),
            "D_pressure": Dp, "L_pressure": Lp}


# ======================================================================================================================
# §6.6 conformal mapping
# ======================================================================================================================
def cot_flow(z) -> dict:
    """The book's two-step map w = ln ζ, ζ = sin z: u − iv = (dw/dζ)(dζ/dz) = cos z/sin z = cot z (chain rule) vs the
    direct derivative of ln(sin z). Returns dict(w, dwdz_chain, dwdz_direct). Book: §6.6. Label: analytic."""
    z = _C(z)
    zeta = np.sin(z)
    chain = (1.0 / zeta) * np.cos(z)
    h = 1e-6
    direct = (np.log(np.sin(z + h)) - np.log(np.sin(z - h))) / (2 * h)
    return {"w": _S(np.log(zeta)), "dwdz_chain": _S(chain), "dwdz_direct": _S(direct)}


def joukowski_ellipse(a: float, b: float) -> dict:
    """The ellipse that the Zhukhovsky map (6.65) makes of the circle |ζ| = a > b (6.66)–(6.67): semi-axes
    A = a + b²/a (along x), B = a − b²/a, foci ±√(A² − B²) = ±2b. a, b [m]. Returns dict(A, B, foci, thickness_ratio).
    Book: §6.6, (6.67), Fig. 6.21. Label: analytic."""
    a, b = float(a), float(b)
    A, B = a + b ** 2 / a, a - b ** 2 / a
    return {"A": A, "B": B, "foci": float(np.sqrt(max(A ** 2 - B ** 2, 0.0))), "thickness_ratio": B / A}


def elliptic_cylinder_flow(U: float = 1.0, a: float = 1.2, b: float = 1.0, *, Gamma_cw: float | None = None,
                           Gamma_ccw: float | None = None, alpha: float = 0.0, branch: str = "outside") -> FunctionFlow:
    """Flow past the elliptic cylinder (6.67) with clockwise circulation Γ: the circle flow (6.68) carried to the
    z-plane by ζ(z) of (6.69) (``branch="outside"``; ``"principal"`` is the wrong variant for E6), u − iv =
    (dw/dζ)(dζ/dz). ``alpha`` tilts the stream (our extension; 0 in the book). Returns a
    :class:`~fluidpy.core.potential.FunctionFlow` (inside = the ellipse; attributes ``zeta_of``, ``ellipse`` and, for
    α = 0, ``stagnation`` = the circle's stagnation points (6.38) mapped by (6.65)).
    U [m/s], a > b [m]. Book: §6.6, (6.65)–(6.69). Validation (planned): V1 u·n = 0 on the ellipse, far field U,
    Blasius L = ρUΓ_cw. Label: analytic."""
    w_z, dw_z = CM.circle_flow_zeta(U, a, Gamma_cw=Gamma_cw, Gamma_ccw=Gamma_ccw, alpha=alpha)
    zeta_of = lambda z: _C(CM.joukowski_inverse(z, b, branch))  # noqa: E731
    dzeta = lambda z: _C(CM.joukowski_inverse_derivative(z, b, branch))  # noqa: E731
    fl = CM.mapped_flow(w_z, dw_z, zeta_of, dzeta, inside_zeta=lambda s: np.abs(s) < a * (1 - 1e-12),
                        u_inf=(U * np.cos(alpha), U * np.sin(alpha)), label="elliptic cylinder")
    fl.zeta_of = zeta_of
    fl.ellipse = joukowski_ellipse(a, b)
    if alpha == 0.0:  # stagnation points of the circle flow (6.38) carried through (6.65)
        fl.stagnation = _C(CM.joukowski(cylinder_stagnation_points(U, a, Gamma_cw=Gamma_cw, Gamma_ccw=Gamma_ccw), b))
    return fl


def joukowski_state(a: float = 1.2, b: float = 1.0, U: float = 1.0, *, Gamma_cw: float = 0.0, x: float = -3.0,
                    y: float = 0.5, branch: str = "outside") -> dict:
    """E6's probe: the inverse Zhukhovsky map (6.69) of z = x + iy on the chosen branch ("outside" — correct —
    or "principal", numpy's root, wrong for Re z < 0), the velocity there of the elliptic-cylinder flow (6.68) carried
    by that branch (u − iv = (dw/dζ)(dζ/dz), not masked, so the wrong branch shows its wrong value), the ellipse (6.67)
    and whether ζ falls inside the circle |ζ| < a (inside the body — the wrong-branch symptom).
    a > b [m]; U [m/s]; Γ_cw [m²/s] clockwise. Returns dict(zeta_re, zeta_im, zeta_abs, u, v, A, B, foci, inside).
    Book: §6.6 (6.65)–(6.69). Label: analytic."""
    z = complex(x, y)
    zeta = complex(CM.joukowski_inverse(z, b, branch))
    _, dwz = CM.circle_flow_zeta(U, a, Gamma_cw=Gamma_cw)
    q = complex(np.asarray(dwz(np.array([zeta]))).ravel()[0]) / (1.0 - b ** 2 / zeta ** 2)
    E = joukowski_ellipse(a, b)
    return {"zeta_re": zeta.real, "zeta_im": zeta.imag, "zeta_abs": abs(zeta), "u": q.real, "v": -q.imag,
            "A": E["A"], "B": E["B"], "foci": E["foci"], "inside": bool(abs(zeta) < a)}


def ellipse_surface_speed(nu, U: float = 1.0, A: float = 1.0, B: float = 0.5):
    """Exact surface speed of the ellipse x = A cos ν, y = B sin ν in a stream U along x, no circulation:
    q = U(A + B)|sin ν|/√(A² sin²ν + B² cos²ν) (from the Zhukhovsky mapped flow (6.68)–(6.69)). Label: analytic."""
    nu = _F(nu)
    return _S(float(U) * (A + B) * np.abs(np.sin(nu)) / np.sqrt(A ** 2 * np.sin(nu) ** 2 + B ** 2 * np.cos(nu) ** 2))


# ======================================================================================================================
# §6.7 numerical Laplace: the 4-point system and Example 6.2
# ======================================================================================================================
def four_point_system(boundary="xy", values=None, as_grid: bool = False, sweeps: int = 0):
    """The 16-point grid of Fig. 6.23 (i, j = 1…4): four unknowns ψ22, ψ32, ψ23, ψ33 and the average rule (6.73) as
    A ψ = b; solved exactly (``np.linalg.solve``) and optionally by ``sweeps`` Gauss–Seidel sweeps from zero.

    ``boundary="xy"`` fills the 12 boundary values with ψ^B = (i − 1)(j − 1) (the discrete-harmonic xy; exact solution
    ψ22 = 1, ψ32 = ψ23 = 2, ψ33 = 4); ``values`` = dict {(i, j): ψ} overrides (book indices). For backward
    compatibility ``boundary`` may itself be such a dict or a 4 × 4 array [j − 1, i − 1].
    Returns dict(A, b, psi (4,) in the order [ψ22, ψ32, ψ23, ψ33], grid (4 × 4 [j, i]), gs_iterates) — or, with
    ``as_grid=True``, the tuple (mask, bc) (4 × 4 [j, i]) for ``core.laplace_solvers.solve_laplace``.
    Book: §6.7, Eq. (6.73), Fig. 6.23. Label: analytic."""
    Bg = np.zeros((4, 4))
    if isinstance(boundary, str):
        if boundary != "xy":
            raise ValueError('boundary must be "xy", a dict or a 4 × 4 array')
        I, J = np.meshgrid(np.arange(1, 5), np.arange(1, 5), indexing="xy")
        Bg = ((I - 1) * (J - 1)).astype(float)
    elif isinstance(boundary, dict):
        values = {**boundary, **(values or {})}
    else:
        Bg = _F(boundary).copy()
    for (i, j), val in (values or {}).items():
        Bg[j - 1, i - 1] = float(val)
    mask = np.zeros((4, 4), bool)
    mask[1:3, 1:3] = True
    if as_grid:
        bc = Bg.copy()
        bc[mask] = 0.0
        return mask, bc
    B = lambda i, j: Bg[j - 1, i - 1]  # noqa: E731
    q = 0.25
    A = np.array([[1, -q, -q, 0], [-q, 1, 0, -q], [-q, 0, 1, -q], [0, -q, -q, 1]], dtype=float)
    b = q * np.array([B(1, 2) + B(2, 1), B(4, 2) + B(3, 1), B(1, 3) + B(2, 4), B(4, 3) + B(3, 4)])  # Eq. (6.73)
    psi = np.linalg.solve(A, b)
    grid = Bg.copy()
    grid[1, 1], grid[1, 2], grid[2, 1], grid[2, 2] = psi[0], psi[1], psi[2], psi[3]
    its, cur = [], np.zeros(4)
    for _ in range(int(sweeps)):
        p22 = q * (B(1, 2) + cur[1] + B(2, 1) + cur[2])
        p32 = q * (p22 + B(4, 2) + B(3, 1) + cur[3])
        p23 = q * (B(1, 3) + cur[3] + p22 + B(2, 4))
        p33 = q * (p23 + B(4, 3) + p32 + B(3, 4))
        cur = np.array([p22, p32, p23, p33])
        its.append(cur.copy())
    return {"A": A, "b": b, "psi": psi, "grid": grid, "gs_iterates": its}


def example_6_2_geometry(refine: int = 1, Q: float = 1.0) -> dict:
    """Grid, mask and boundary values of Example 6.2's sharp contraction (Fig. 6.24) with spacing Δ = 1/refine m:
    channel 0 ≤ x ≤ 9 m, 0 ≤ y ≤ 5 m with the corner block x ≥ 5, y ≤ 2 removed; ψ = 0 on the lower walls, ψ = Q on
    the top, ψ = Qy/5 at the inlet (x = 0) and Q(y − 2)/3 at the outlet (x = 9) (uniform inlet/outlet velocity).
    Returns dict(x, y, mask ([j, i], True = unknown), bc (NaN in the solid), dx). Book: §6.7, Example 6.2 (the book's
    grid is refine = 1: 10 × 6 nodes). Q [m²/s] is ours (default 1; the book's flow rate stays in the private test
    data, and ψ scales linearly with Q). Label: analytic."""
    r = int(refine)
    dx = 1.0 / r
    x = np.arange(9 * r + 1) * dx
    y = np.arange(5 * r + 1) * dx
    X, Y = np.meshgrid(x, y, indexing="xy")
    tol = 1e-9
    solid = (X > 5 + tol) & (Y < 2 - tol)
    bc = np.full(X.shape, np.nan)
    bc[np.isclose(Y, 0) & (X <= 5 + tol)] = 0.0  # lower wall
    bc[np.isclose(X, 5) & (Y <= 2 + tol)] = 0.0  # step face
    bc[np.isclose(Y, 2) & (X >= 5 - tol)] = 0.0  # lower wall of the outlet channel
    bc[np.isclose(Y, 5)] = Q  # upper wall
    inlet = np.isclose(X, 0)
    bc[inlet] = Q * Y[inlet] / 5.0
    outlet = np.isclose(X, 9) & (Y >= 2 - tol)
    bc[outlet] = Q * (Y[outlet] - 2.0) / 3.0
    mask = (X > tol) & (X < 9 - tol) & (Y > tol) & (Y < 5 - tol) & ~((X >= 5 - tol) & (Y <= 2 + tol))
    bc[mask] = 0.0
    return {"x": x, "y": y, "mask": mask, "bc": bc, "dx": dx, "solid": solid}


def example_6_2(Q: float = 1.0, n_iter: int | None = None, tol: float = 1e-10, method: str = "gauss_seidel",
                refine: int = 1, omega: float | None = None, order: str = "book") -> dict:
    """Example 6.2: stream function in the sharp-cornered contraction by the average rule (6.72), swept Gauss–Seidel
    with the latest values in the book's order (outer loop over i, the x index) from ψ = 0 inside.

    Q [m²/s] the flow rate per unit depth (ours, default 1 — ψ scales with Q; the book's value is private test data);
    ``n_iter`` = a fixed number of sweeps (the book's style) or None = until the residual ≤ tol; ``method`` as in
    :func:`~fluidpy.core.laplace_solvers.solve_laplace`; ``refine`` divides the 1 m spacing (×8 practical with "direct"
    or "sor").
    Returns dict(psi ([j, i], NaN in the solid), X, Y (meshgrids [m]), x, y, mask, history, iterations, grid_shape,
    probe (the (j, i) index of the fixed interior point (x, y) = (2 m, 2 m) — the same physical point at every
    refine), flux (ψ_top − ψ_bottom from u = ∂ψ/∂y on every interior vertical grid line — equals Q), dx).
    Book: §6.7, Example 6.2, Figs. 6.24–6.25. Validation (planned): V6 book grid values (private Q); V3 refinement (order
    reduced near the 270° corner); V4 flux. Label: converged, conserved.
    """
    g = example_6_2_geometry(refine, Q)
    psi, hist = LS.solve_laplace(g["mask"], np.nan_to_num(g["bc"]), method=method, tol=tol, n_iter=n_iter,
                                 omega=omega, dx=g["dx"], order=order)
    psi = np.where(g["solid"], np.nan, psi)
    flux = []
    for i in range(1, psi.shape[1] - 1):
        col = psi[:, i]
        ok = np.isfinite(col)
        u = np.diff(col[ok]) / g["dx"]
        flux.append(float(np.sum(u * g["dx"])))
    X, Y = np.meshgrid(g["x"], g["y"], indexing="xy")
    r = int(refine)
    return {"psi": psi, "X": X, "Y": Y, "x": g["x"], "y": g["y"], "mask": g["mask"], "history": hist,
            "iterations": hist["iterations"], "grid_shape": psi.shape, "probe": (2 * r, 2 * r),
            "flux": np.array(flux), "dx": g["dx"]}


def relaxation_state(problem: str = "four_point", method: str = "gauss_seidel", sweeps: int = 1,
                     omega: float | None = None, refine: int = 1, order: str = "lex") -> dict:
    """E7's state: run ``sweeps`` relaxation sweeps from ψ = 0 on one of three problems and report the residual, the
    last change, the four-point values and the sweeps needed to reach a residual of 1e-10.

    problems: "four_point" (Fig. 6.23 with ψ^B = (i − 1)(j − 1), :func:`four_point_system`), "contraction" (Example
    6.2's geometry with Q = 1 m²/s, ×refine), "xy_square" (8 × 8 nodes, ψ = xy on the frame, Δ = 1 m — discrete-harmonic,
    so the converged grid is exact). method "jacobi", "gauss_seidel" or "sor" (ω default: ω_opt of the grid's bounding
    box, ``core.laplace_solvers.optimal_sor_omega``); ``order`` "lex" (i fastest, j slowest) or "book".
    Returns dict(residual (max defect of the average rule after the sweeps, ``residual_norm``), max_change (of the last
    sweep), psi22, psi32, psi23, psi33 (the four-point values; NaN for the other problems), psi_probe (ψ at the
    problem's probe node: ψ33 for four_point, Example 6.2's probe, node (4, 4) of the square), sweeps_to_tol (sweeps
    from zero to residual ≤ 1e-10, same method)). Book: §6.7 (6.72)–(6.73), N73. Label: analytic."""
    if problem == "four_point":
        mask, bc = four_point_system("xy", as_grid=True)
        probe, dx = (2, 2), 1.0
    elif problem == "contraction":
        g = example_6_2_geometry(refine, 1.0)
        mask, bc, dx = g["mask"], np.nan_to_num(g["bc"]), g["dx"]
        probe = (2 * int(refine), 2 * int(refine))
    elif problem == "xy_square":
        I, J = np.meshgrid(np.arange(8.0), np.arange(8.0), indexing="xy")
        bc = I * J
        mask = np.zeros((8, 8), bool)
        mask[1:-1, 1:-1] = True
        bc = np.where(mask, 0.0, bc)
        probe, dx = (4, 4), 1.0
    else:
        raise ValueError('problem must be "four_point", "contraction" or "xy_square"')
    if method == "sor" and omega is None:
        omega = LS.optimal_sor_omega(mask.shape[1], mask.shape[0])
    psi = np.where(mask, 0.0, bc)
    ch = 0.0
    for _ in range(int(sweeps)):
        if method == "jacobi":
            psi, ch = LS.jacobi_sweep(psi, mask, dx=dx)
        elif method == "gauss_seidel":
            psi, ch = LS.gauss_seidel_sweep(psi, mask, order=order, dx=dx)
        elif method == "sor":
            psi, ch = LS.sor_sweep(psi, mask, omega=omega, dx=dx, order=order)
        else:
            raise ValueError('method must be "jacobi", "gauss_seidel" or "sor"')
    res = LS.residual_norm(psi, mask, dx)
    _, hist = LS.solve_laplace(mask, bc, method=method, tol=1e-10, omega=omega, dx=dx,
                               order="y_outer" if order == "lex" else order)
    four = (psi[1, 1], psi[1, 2], psi[2, 1], psi[2, 2]) if problem == "four_point" else (np.nan,) * 4
    return {"residual": res, "max_change": float(ch), "psi22": float(four[0]), "psi32": float(four[1]),
            "psi23": float(four[2]), "psi33": float(four[3]), "psi_probe": float(psi[probe]),
            "sweeps_to_tol": int(hist["sweeps"])}


# ======================================================================================================================
# §6.8 axisymmetric flow: Stokes stream function, sphere, airship, axial method
# ======================================================================================================================
def _d4(f, a, b, h, axis):
    c1 = [1.0, -8.0, 8.0, -1.0]
    s1 = [-2.0, -1.0, 1.0, 2.0]
    c2 = [-1.0, 16.0, -30.0, 16.0, -1.0]
    s2 = [-2.0, -1.0, 0.0, 1.0, 2.0]
    if axis == 0:
        d1 = sum(c * _F(f(a + s * h, b)) for c, s in zip(c1, s1)) / (12 * h)
        d2 = sum(c * _F(f(a + s * h, b)) for c, s in zip(c2, s2)) / (12 * h ** 2)
    else:
        d1 = sum(c * _F(f(a, b + s * h)) for c, s in zip(c1, s1)) / (12 * h)
        d2 = sum(c * _F(f(a, b + s * h)) for c, s in zip(c2, s2)) / (12 * h ** 2)
    return d1, d2


def stokes_operator_residual(psi_fn: Callable, R, z, h: float = 1e-3):
    """Field equation of the Stokes stream function in irrotational flow (6.77):
    ∂/∂R((1/R)∂ψ/∂R) + (1/R)∂²ψ/∂z² = −ω_φ (= 0), fourth-order differences. psi_fn(R, z) [m³/s]; R > 2h, z [m].
    Returns the left side [1/s] (= −ω_φ). Not the plane Laplacian (so complex variables do not apply): ψ = Rz, harmonic
    in the plane sense, gives −z/R² ≠ 0 (a control field); note ψ = R²z *is* a solution (axisymmetric stagnation flow).
    Book: §6.8, Eqs. (6.75)–(6.77). Label: analytic."""
    R_, z_ = _F(R), _F(z)
    pR, pRR = _d4(psi_fn, R_, z_, h, 0)
    _, pzz = _d4(psi_fn, R_, z_, h, 1)
    return _S(pRR / R_ - pR / R_ ** 2 + pzz / R_)  # Eq. (6.77)


def stokes_operator_sym(psi_expr, R: sp.Symbol, z: sp.Symbol):
    """sympy twin of (6.77): ∂/∂R((1/R)∂ψ/∂R) + (1/R)∂²ψ/∂z². Label: symbolic."""
    return sp.simplify(sp.diff(sp.diff(psi_expr, R) / R, R) + sp.diff(psi_expr, z, 2) / R)


def stokes_operator_residual_sym(psi_expr, R: sp.Symbol, z: sp.Symbol):
    """sympy twin of :func:`stokes_operator_residual`: ∂/∂R((1/R)∂ψ/∂R) + (1/R)∂²ψ/∂z² (= −ω_φ, (6.77)); zero for the
    elements (6.86)–(6.89) and for R²z (axisymmetric stagnation flow), nonzero for R z. Label: symbolic."""
    return stokes_operator_sym(psi_expr, R, z)


def axisym_laplacian_residual(phi_fn: Callable, R, z, h: float = 1e-3):
    """Axisymmetric Laplace equation for φ (6.80): (1/R)∂/∂R(R∂φ/∂R) + ∂²φ/∂z² (fourth-order differences) [1/s].
    Book: §6.8 (6.80). Label: analytic."""
    R_, z_ = _F(R), _F(z)
    pR, pRR = _d4(phi_fn, R_, z_, h, 0)
    _, pzz = _d4(phi_fn, R_, z_, h, 1)
    return _S(pRR + pR / R_ + pzz)  # Eq. (6.80)


def axisym_laplacian_sym(phi_expr, R: sp.Symbol, z: sp.Symbol):
    """sympy twin of (6.80). Label: symbolic."""
    return sp.simplify(sp.diff(R * sp.diff(phi_expr, R), R) / R + sp.diff(phi_expr, z, 2))


def spherical_continuity_residual(u_r_fn: Callable, u_t_fn: Callable, r, theta, h: float = 1e-5):
    """Axisymmetric spherical continuity in the Appendix-B form (1/r²)∂(r²u_r)/∂r + (1/r sin θ)∂(u_θ sin θ)/∂θ [1/s]
    — the book's (6.82) prints 1/r for the first factor (same zero set; a residual needs App. B's 1/r² and the 1/r of the
    second term). Book: §6.8 (6.82) (corrected factor). Label: analytic."""
    r_, t_ = _F(r), _F(theta)
    d1 = ((r_ + h) ** 2 * _F(u_r_fn(r_ + h, t_)) - (r_ - h) ** 2 * _F(u_r_fn(r_ - h, t_))) / (2 * h)
    d2 = (_F(u_t_fn(r_, t_ + h)) * np.sin(t_ + h) - _F(u_t_fn(r_, t_ - h)) * np.sin(t_ - h)) / (2 * h)
    return _S(d1 / r_ ** 2 + d2 / (r_ * np.sin(t_)))


def axisym_flux_between(flow_or_psi, P1, P2, n: int = 64, h: float = 1e-6) -> tuple:
    """Volume flow rate through the surface of revolution swept by the meridian segment P1 → P2 (P = (R, z)):
    Q = ∫2πR(u·n)ds = ∫2πR(−u_R dz + u_z dR) (Gauss–Legendre, n nodes) vs 2π[ψ(P2) − ψ(P1)] (6.78).
    ``flow_or_psi``: an :class:`~fluidpy.core.potential.AxisymFlow` or a callable ψ(R, z) (velocities by differences).
    Returns (quad, two_pi_dpsi) [m³/s]. Book: §6.8 (6.78), Fig. 6.26. Validation (planned): V1 (wrong variant without
    2π fails). Label: analytic."""
    R1, z1 = (float(v) for v in P1)
    R2, z2 = (float(v) for v in P2)
    s, ws = np.polynomial.legendre.leggauss(int(n))
    s = 0.5 * (s + 1.0)
    ws = 0.5 * ws
    R = R1 + s * (R2 - R1)
    z = z1 + s * (z2 - z1)
    if hasattr(flow_or_psi, "velocity_cyl"):
        uR, uz = (_F(q) for q in flow_or_psi.velocity_cyl(R, z))
        psi = flow_or_psi.psi
    else:
        psi = flow_or_psi
        uR = -(_F(psi(R, z + h)) - _F(psi(R, z - h))) / (2 * h) / R
        uz = (_F(psi(R + h, z)) - _F(psi(R - h, z))) / (2 * h) / R
    integrand = _TWO_PI * R * (-uR * (z2 - z1) + uz * (R2 - R1))  # Eq. (6.78) with ds = (dR, dz)
    return float(np.sum(ws * integrand)), float(_TWO_PI * (float(psi(R2, z2)) - float(psi(R1, z1))))


def sphere_surface_cp(theta):
    """Surface pressure coefficient of the sphere C_p = 1 − (9/4) sin²θ (6.91) (min −5/4 at θ = 90°, speed 1.5U);
    fore–aft symmetric ⇒ no drag. θ [rad] from +z (the downstream axis). Book: §6.8 (6.91).
    Label: analytic."""
    return _S(1.0 - 2.25 * np.sin(_F(theta)) ** 2)  # Eq. (6.91)


def cylinder_vs_sphere(r_over_a: float) -> dict:
    """3-D relief in numbers: at distance r = (r/a)·a, the largest velocity perturbation relative to U — (a/r)² for the
    cylinder ((6.34): |u − U| = Ua²/r² everywhere on the circle) and (a/r)³ for the sphere ((6.90), on the axis) — and
    the surface extremes, C_p min −3 vs −5/4 ((6.35), (6.91)), max speed 2U vs 1.5U.
    Returns dict(cyl_pert, sph_pert, cyl_cp_min, sph_cp_min, cyl_speed_max, sph_speed_max) (speeds in units of U).
    Book: §6.3 (6.34)–(6.35), §6.8 (6.90)–(6.91). Label: analytic."""
    q = 1.0 / float(r_over_a)
    return {"cyl_pert": q ** 2, "sph_pert": q ** 3, "cyl_cp_min": float(cylinder_surface_cp(0.5 * np.pi)),
            "sph_cp_min": float(sphere_surface_cp(0.5 * np.pi)), "cyl_speed_max": 2.0, "sph_speed_max": 1.5}


def line_sink_stream_function(R, z, k: float = 1.0, a: float = 1.0, method: str = "closed"):
    """Stokes stream function of the book's uniform line sink of density k [m²/s] from O (z = 0) to A (z = a):
    (6.93) ψ = (k/4π)∫₀ᵃ cos α dξ, cos α = (z − ξ)/√(R² + (z − ξ)²) (``method="quad"``), and its closed form (6.94)
    ψ = (k/4π)(r − r₁) with r = |OP|, r₁ = |AP| (``"closed"``). Returns ψ [m³/s]. Book: §6.8 (6.93)–(6.94), Fig. 6.28.
    Validation (planned): V1 closed = quad. Label: analytic."""
    R_, z_ = np.broadcast_arrays(_F(R), _F(z))
    if method == "closed":
        r = np.sqrt(R_ ** 2 + z_ ** 2)
        r1 = np.sqrt(R_ ** 2 + (z_ - float(a)) ** 2)
        return _S(float(k) / _FOUR_PI * (r - r1))  # Eq. (6.94)
    if method == "quad":
        out = np.empty(R_.shape)
        for idx in np.ndindex(R_.shape):
            Rv, zv = float(R_[idx]), float(z_[idx])
            f = lambda xi: (zv - xi) / np.hypot(Rv, zv - xi)  # noqa: E731  cos α
            out[idx] = float(k) / _FOUR_PI * quad(f, 0.0, float(a), epsabs=0.0, epsrel=1e-12, limit=200)[0]  # (6.93)
        return _S(out)
    raise ValueError('method must be "closed" or "quad"')


def _axis_root(fn, lo, hi):
    return float(brentq(fn, lo, hi, xtol=1e-14, maxiter=500))


def _body_radius(psi, z, R_hi):
    """Radius of the stream surface ψ = 0 at station z (root of ψ/R² in (R_lo, R_hi))."""
    g = lambda R: float(psi(R, z)) / R ** 2  # noqa: E731
    lo = 1e-7 * R_hi
    if g(lo) * g(R_hi) > 0:
        return 0.0
    return float(brentq(g, lo, R_hi, xtol=1e-13))


def airship(U: float = 1.0, Q: float = 1.0, a: float = 1.0, n_contour: int = 121) -> dict:
    """Source–line-sink body of revolution (6.95), Fig. 6.28: point source Q at O, uniform line sink k = Q/a from O to
    A (closure Q = ak) and a stream U along z: ψ = −(Q/4π) cos θ + (Q/4πa)(r − r₁) + ½Ur² sin²θ.

    Returns dict(flow (AxisymFlow, inside = body), psi = psi_fn (callable (R, z)), psi_fluid_frame (without the
    stream), z_nose = stagnation_front, z_tail = stagnation_rear (axis stagnation points by ``brentq``), length, R_max,
    contour (z, R) of ψ = 0, closure = net_strength = Q − ak (= 0)).
    U [m/s], Q [m³/s], a [m]. Book: §6.8 (6.93)–(6.95); Exercise 6.42 (length). Label: analytic."""
    U, Q, a = float(U), float(Q), float(a)
    els = [PointSource3D(Q, 0.0), LineSource3D(-Q / a, 0.0, a)]
    free = AxisymFlow(els)
    full = AxisymFlow([AxisymUniform(U)] + els)

    def psi(R, z):
        return full.psi(R, z)

    def uz_axis(z):
        return float(full.velocity_cyl(0.0, z)[1])

    Lsc = max(a, np.sqrt(Q / (_FOUR_PI * U)))
    z_nose = _axis_root(uz_axis, -50.0 * Lsc, -1e-9 * Lsc)
    z_tail = _axis_root(uz_axis, a * (1 + 1e-9), a + 50.0 * Lsc)
    zc = z_nose + (z_tail - z_nose) * 0.5 * (1 - np.cos(np.linspace(0, np.pi, int(n_contour))))
    Rhi = 10.0 * Lsc
    Rc = np.array([_body_radius(psi, zz, Rhi) for zz in zc[1:-1]])
    Rc = np.concatenate([[0.0], Rc, [0.0]])
    zs = np.array(zc)
    Rs = np.array(Rc)

    def inside(R, z):
        R_, z_ = np.broadcast_arrays(_F(R), _F(z))
        Rb = np.interp(z_, zs, Rs, left=0.0, right=0.0)
        return (R_ < Rb * (1 - 1e-9)) & (z_ > z_nose) & (z_ < z_tail)

    body = AxisymFlow([AxisymUniform(U)] + els, inside=inside, label="airship")
    return {"flow": body, "psi": psi, "psi_fn": psi, "psi_fluid_frame": free.psi, "z_nose": z_nose,
            "z_tail": z_tail, "stagnation_front": z_nose, "stagnation_rear": z_tail, "length": z_tail - z_nose,
            "R_max": float(np.max(Rs)), "contour": (zs, Rs), "net_strength": Q - (Q / a) * a,
            "closure": Q - (Q / a) * a}


def axisym_half_body(U: float = 1.0, Q: float = 1.0) -> dict:
    """3-D half-body (Exercise 6.36): point source Q in a stream U: nose at z = −√(Q/4πU), body ψ = Q/4π,
    far-downstream radius √(Q/πU). Returns dict(flow, z_stag, psi_body, R_far). Book: §6.8 (6.86)–(6.87).
    Label: analytic."""
    U, Q = float(U), float(Q)
    return {"flow": AxisymFlow([AxisymUniform(U), PointSource3D(Q)]), "z_stag": -np.sqrt(Q / (_FOUR_PI * U)),
            "psi_body": Q / _FOUR_PI, "R_far": float(np.sqrt(Q / (np.pi * U)))}


def perturbation_radius(kind: str = "cylinder", frac: float = 0.01) -> float:
    """Distance r/a beyond which the body's velocity perturbation is below ``frac`` of U (Exercise 6.38-type):
    cylinder |u − U| = Ua²/r² ⇒ r/a = frac^{−1/2}; sphere max|u − U| = Ua³/r³ (on the axis) ⇒ r/a = frac^{−1/3}.
    Book: §6.3 (6.34), §6.8 (6.90). Label: analytic."""
    if kind == "cylinder":
        return float(frac ** -0.5)
    if kind == "sphere":
        return float(frac ** (-1.0 / 3.0))
    raise ValueError('kind must be "cylinder" or "sphere"')


def ellipsoid_linear_source_strength(U: float = 1.0, A: float = 2.0, B: float = 1.0) -> dict:
    """Prolate spheroid (semi-axes A along the stream, B) = stream U + axial line source of linear density
    k(ξ) = Kξ between the foci ±c, c = √(A² − B²) (a classical exact result); K from u_z = 0 at the nose z = −A:
    K = −4πU/[2Ac/B² − ln((A + c)/(A − c))]. Returns dict(K [1/s], c [m]). Book: §6.8 (inverse method, our V1 target).
    Label: analytic."""
    U, A, B = float(U), float(A), float(B)
    c = float(np.sqrt(A ** 2 - B ** 2))
    K = -_FOUR_PI * U / (2 * A * c / B ** 2 - np.log((A + c) / (A - c)))
    return {"K": float(K), "c": c}


def _collocation(xi, z_lo, z_hi, how):
    N = xi.size - 1
    mid = 0.5 * (xi[:-1] + xi[1:])
    if how == "mid":
        return mid
    if how == "scaled":
        return z_lo + (mid - xi[0]) / (xi[-1] - xi[0]) * (z_hi - z_lo)
    if how == "cosine":
        return np.sort(0.5 * (z_lo + z_hi) + 0.5 * (z_hi - z_lo) * np.cos(np.pi * (np.arange(1, N + 1) - 0.5) / N))
    raise ValueError('collocation must be "mid", "scaled" or "cosine"')


def axisym_body_target(kind: str = "rankine_oval", N: int = 40, U: float = 1.0, collocation: str = "mid",
                       span: str = "body", as_dict: bool = False, fineness: float | None = None, **p):
    """Targets for the axial singularity method (E8, C14): body points and axial segment nodes.

    kinds: "ellipsoid" (A=1, B=0.2; exact smooth distribution k = Kξ between the foci ±c — the V1/V3 target),
    "airship" (Q=1, a=1; exact = point source + uniform line sink on [0, a]), "rankine_oval" (Q=1, c=1: point source at
    −c, sink at +c), "sphere" (a=1: a point doublet — ill-conditioned, fenced qualitative). N segments uniform on the
    singular interval; N body points at ``collocation`` = "mid" (default: the z of each segment midpoint, the body
    point above its segment as in Fig. 6.29), "scaled" (midpoints stretched onto the body length) or "cosine"
    (clustered at the ends — ⚠️ the influence matrix then becomes ill-conditioned fast: cond ~ 10⁸ at N = 16 for a
    fineness ratio of 2). The system is a first-kind integral equation: its condition number grows with N and with
    bluntness (the solution reports it).
    ``span="body"`` (default, the design's layout): the N segments divide the body from nose to tail and the body points
    sit above their midpoints (so ``axial_singularity_solve(z_body, R_body, U, N)`` rebuilds the same segments);
    ``span="singular"``: segments on the singular interval above (foci, [0, a], [−c, c], ±0.9a). ``fineness`` = A/B of
    the ellipsoid (A = 1). Returns (z_body, R_body), or with ``as_dict=True`` dict(z_body, R_body, xi_nodes, exact_k (at
    segment midpoints or None), kind, U, label, R_fn (exact body radius R(z)), z_ends (nose, tail)).
    Book: §6.8, Fig. 6.29. Label: analytic."""
    out = _axisym_target_dict(kind, int(N), U, collocation, fineness, **p)
    if span == "body":
        zn, zt = out["z_ends"]
        xi = np.linspace(zn, zt, int(N) + 1)
        zb = 0.5 * (xi[:-1] + xi[1:])
        out.update(z_body=zb, R_body=_F(out["R_fn"](zb)), xi_nodes=xi,
                   exact_k=None if kind != "ellipsoid" else None)
    elif span != "singular":
        raise ValueError('span must be "body" or "singular"')
    return out if as_dict else (out["z_body"], out["R_body"])


def _axisym_target_dict(kind, N, U, collocation, fineness, **p) -> dict:
    N = int(N)
    if kind == "ellipsoid":
        A = float(p.get("A", 1.0))
        B = A / float(fineness) if fineness is not None else float(p.get("B", 0.2))
        e = ellipsoid_linear_source_strength(U, A, B)
        xi = np.linspace(-e["c"], e["c"], N + 1)
        zb = _collocation(xi, -A, A, collocation)
        Rb = B * np.sqrt(1 - (zb / A) ** 2)
        exact = e["K"] * 0.5 * (xi[:-1] + xi[1:])
        return {"z_body": zb, "R_body": Rb, "xi_nodes": xi, "exact_k": exact, "kind": kind, "U": U,
                "label": f"prolate spheroid A = {A:g}, B = {B:g}", "exact_net": 0.0,
                "R_fn": lambda z, A=A, B=B: B * np.sqrt(np.clip(1 - (_F(z) / A) ** 2, 0.0, None)), "z_ends": (-A, A)}
    if kind == "airship":
        Q, a = float(p.get("Q", 1.0)), float(p.get("a", 1.0))
        s = airship(U, Q, a)
        xi = np.linspace(0.0, a, N + 1)
        zb = _collocation(xi, s["z_nose"], s["z_tail"], collocation)
        Rb = np.array([_body_radius(s["psi"], zz, 10 * max(a, s["R_max"])) for zz in zb])
        return {"z_body": zb, "R_body": Rb, "xi_nodes": xi, "exact_k": None, "kind": kind, "U": U,
                "label": "airship (6.95)", "exact_net": 0.0, "z_ends": (s["z_nose"], s["z_tail"]),
                "R_fn": lambda z, s=s, a=a: np.array([_body_radius(s["psi"], q, 10 * max(a, s["R_max"]))
                                                      for q in np.atleast_1d(_F(z))])}
    if kind == "rankine_oval":
        Q, c = float(p.get("Q", 1.0)), float(p.get("c", 1.0))
        fl = AxisymFlow([AxisymUniform(U), PointSource3D(Q, -c), PointSource3D(-Q, c)])
        uz = lambda z: float(fl.velocity_cyl(0.0, z)[1])  # noqa: E731
        zt = _axis_root(uz, c * (1 + 1e-9), c + 50.0 * max(c, np.sqrt(Q / U)))
        xi = np.linspace(-c, c, N + 1)
        zb = _collocation(xi, -zt, zt, collocation)
        Rb = np.array([_body_radius(fl.psi, zz, 10 * zt) for zz in zb])
        return {"z_body": zb, "R_body": Rb, "xi_nodes": xi, "exact_k": None, "kind": kind, "U": U,
                "label": "3-D Rankine oval", "exact_net": 0.0, "z_ends": (-zt, zt),
                "R_fn": lambda z, fl=fl, zt=zt: np.array([_body_radius(fl.psi, q, 10 * zt) for q in np.atleast_1d(_F(z))])}
    if kind == "sphere":
        a = float(p.get("a", 1.0))
        xi = np.linspace(-0.9 * a, 0.9 * a, N + 1)
        zb = _collocation(xi, -a, a, collocation)
        Rb = np.sqrt(a ** 2 - zb ** 2)
        return {"z_body": zb, "R_body": Rb, "xi_nodes": xi, "exact_k": None, "kind": kind, "U": U,
                "label": "sphere (a point doublet: fenced qualitative)", "exact_net": 0.0, "z_ends": (-a, a),
                "R_fn": lambda z, a=a: np.sqrt(np.clip(a ** 2 - _F(z) ** 2, 0.0, None))}
    raise ValueError('kind must be "ellipsoid", "airship", "rankine_oval" or "sphere"')


def axisym_body_fit(kind: str = "ellipsoid", N: int = 20, U: float = 1.0, collocation: str = "mid", **p) -> dict:
    """Run the axial singularity method on a target (:func:`axisym_body_target` + ``core.panels``): returns the target,
    the solution, when an exact distribution exists k_error = max|k_n − k_exact|/max|k_exact|, and for the ellipsoid
    body_psi_error = max|ψ|/(½UB²) on 201 points of the true body (not only the N collocation points).
    Default (ellipsoid A = 1, B = 0.2, "mid"): k_error ≈ 0.14, 0.055, 0.022, 0.010 for N = 4, 8, 16, 32 (end effects of
    piecewise-constant segments, observed order → 1) while the body error falls about as N⁻³; cond grows with N.
    Book: §6.8. Label: analytic, converged."""
    tg = axisym_body_target(kind, N, U, collocation, span="singular", as_dict=True, **p)
    sol = PN.axial_singularity_solve(tg["z_body"], tg["R_body"], U, tg["xi_nodes"])
    err = None
    if tg["exact_k"] is not None:
        err = float(np.max(np.abs(sol["k"] - tg["exact_k"])) / np.max(np.abs(tg["exact_k"])))
    body_err = None
    if kind == "ellipsoid":
        A, B = float(p.get("A", 1.0)), float(p.get("B", 0.2))
        zz = np.linspace(-0.99 * A, 0.99 * A, 201)
        RR = B * np.sqrt(1 - (zz / A) ** 2)
        body_err = float(np.max(np.abs(_F(PN.axial_singularity_psi(RR, zz, sol["k"], tg["xi_nodes"], U))))
                         / (0.5 * U * B ** 2))
    return {"target": tg, "solution": sol, "k_error": err, "body_psi_error": body_err}


def axial_state(target: str = "rankine_oval", N: int = 20, U: float = 1.0, fineness: float = 3.0) -> dict:
    """E8's state: the axial singularity method on a target body (:func:`axisym_body_target`, nose-to-tail segments,
    body points above the segment midpoints) solved by ``core.panels.axial_singularity_solve``.

    target "rankine_oval", "airship", "sphere", "ellipsoid" (``fineness`` = A/B with A = 1 m); N segments; U [m/s].
    Returns dict(cond (1-norm condition number ``np.linalg.cond(A, 1)``), net_strength Σk_nΔξ [m³/s] (→ 0: closed),
    body_error [m] = max |R_c(z_q) − R_target(z_q)| at the N − 1 stations z_q halfway between consecutive body points,
    where R_c is the root of the computed ψ(R, z_q) nearest the target radius R_t: ψ is sampled at R = R_t(0.5 + j/100),
    j = 0…100, the sign change closest to R_t is refined by ``brentq`` (no sign change: R_c = 0, i.e. the full radius
    counts as error), k_max, k_min [m²/s]).
    Book: §6.8, Fig. 6.29. Label: analytic, converged (the sphere: qualitative, ill-conditioned)."""
    tg = axisym_body_target(target, int(N), U, fineness=fineness if target == "ellipsoid" else None, as_dict=True)
    sol = PN.axial_singularity_solve(tg["z_body"], tg["R_body"], U, int(N))
    zq = 0.5 * (tg["z_body"][:-1] + tg["z_body"][1:])
    Rt = _F(tg["R_fn"](zq))

    def radius_near(zz, rt):
        Rs = rt * (0.5 + np.arange(101) / 100.0)
        ps = _F(PN.axial_singularity_psi(sol, Rs, np.full(Rs.shape, zz)))
        idx = np.nonzero(np.sign(ps[:-1]) * np.sign(ps[1:]) <= 0)[0]
        if idx.size == 0:
            return 0.0
        i0 = idx[np.argmin(np.abs(0.5 * (Rs[idx] + Rs[idx + 1]) - rt))]
        if ps[i0] == 0.0:
            return float(Rs[i0])
        return float(brentq(lambda R: float(PN.axial_singularity_psi(sol, R, zz)), Rs[i0], Rs[i0 + 1], xtol=1e-14))

    Rc = np.array([radius_near(zz, rt) for zz, rt in zip(zq, Rt)])
    return {"cond": sol["cond"], "net_strength": sol["net_strength"], "body_error": float(np.max(np.abs(Rc - Rt))),
            "k_max": float(np.max(sol["k"])), "k_min": float(np.min(sol["k"]))}


def panel_cp_error(N: int = 32, body: str | None = None, U: float = 1.0, A: float = 1.0, B: float = 0.5,
                   kind: str = "ellipse") -> float:
    """Max |C_p(panels) − C_p(exact)| at the control points of N constant-strength source panels (our extension).
    "circle": vertices on the unit circle, exact C_p = 1 − 4 sin²θ at the control-point angle — the panel result agrees
    to round-off for every N (a symmetry of the regular polygon, not a convergence study); "ellipse": vertices at
    x = A cos ν, y = B sin ν, exact speed from :func:`ellipse_surface_speed` at the control point's ν
    (ν = atan2(y/B, x/A)) — algebraic convergence in N (order 2). ``body`` is an alias of ``kind``; the default is the
    ellipse because the circle's error is round-off. Label: converged."""
    kind = body if body is not None else kind
    N = int(N)
    nu = np.linspace(0, _TWO_PI, N, endpoint=False) + np.pi / N
    if kind == "circle":
        r = PN.source_panels(np.cos(nu), np.sin(nu), U)
        return float(np.max(np.abs(r["cp"] - (1 - 4 * np.sin(r["theta"]) ** 2))))
    if kind == "ellipse":
        r = PN.source_panels(A * np.cos(nu), B * np.sin(nu), U)
        nuc = np.arctan2(r["ym"] / B, r["xm"] / A)
        cp_ex = 1.0 - (_F(ellipse_surface_speed(nuc, U, A, B)) / U) ** 2
        return float(np.max(np.abs(r["cp"] - cp_ex)))
    raise ValueError('kind must be "circle" or "ellipse"')


# ======================================================================================================================
# §6.9 the accelerating sphere and added mass
# ======================================================================================================================
def moving_sphere_dphidt_sym() -> dict:
    """sympy check of (6.101): for φ = −(a³/2|ξ|³)u_s·ξ with ξ = x − x_s(t), the chain rule ∂φ/∂t =
    (∂φ/∂x_s)·u_s + (∂φ/∂u_s)·du_s/dt equals −u·u_s − (a³/2|ξ|³) ξ·du_s/dt with u = ∇φ; also direct differentiation in
    t for a sphere with x_s = (X(t), 0, 0), u_s = (X′, 0, 0); and (6.104): the book's printed bracket minus ∇φ on the
    surface (nonzero: u_s, the sign slip) and the corrected bracket minus ∇φ (zero). Returns dict(chain, direct,
    printed_bracket, correct_bracket) (differences; chain, direct, correct_bracket simplify to 0) plus the older keys
    difference_chain, difference_direct. Book: §6.9 (6.101)–(6.104). Label: symbolic."""
    x, y, zz, xs, ys, zs, us, vs, ws, ax, ay, az = sp.symbols("x y z x_s y_s z_s u_s v_s w_s a_x a_y a_z", real=True)
    a = sp.symbols("a", positive=True)
    X = sp.Matrix([x, y, zz])
    Xs = sp.Matrix([xs, ys, zs])
    Us = sp.Matrix([us, vs, ws])
    Acc = sp.Matrix([ax, ay, az])
    xi = X - Xs
    r = sp.sqrt(xi.dot(xi))
    phi = -a ** 3 / (2 * r ** 3) * Us.dot(xi)  # Eq. (6.97)
    u = sp.Matrix([sp.diff(phi, v) for v in (x, y, zz)])
    chain = sum(sp.diff(phi, s) * c for s, c in zip((xs, ys, zs), Us)) + sum(sp.diff(phi, s) * c for s, c in
                                                                              zip((us, vs, ws), Acc))
    book = -u.dot(Us) - a ** 3 / (2 * r ** 3) * xi.dot(Acc)  # Eq. (6.101)
    t = sp.symbols("t", real=True)
    Xt = sp.Function("X")(t)
    phit = -a ** 3 / (2 * sp.sqrt((x - Xt) ** 2 + y ** 2 + zz ** 2) ** 3) * sp.diff(Xt, t) * (x - Xt)
    direct = sp.diff(phit, t)
    sub = {xs: Xt, ys: 0, zs: 0, us: sp.diff(Xt, t), vs: 0, ws: 0, ax: sp.diff(Xt, t, 2), ay: 0, az: 0}
    e1, e2, e3 = sp.symbols("e_1 e_2 e_3", real=True)
    E = sp.Matrix([e1, e2, e3])
    grad = u.subs({x: xs + a * e1, y: ys + a * e2, zz: zs + a * e3})
    printed = -a ** 3 / 2 * (-3 * a * E / a ** 5 * Us.dot(a * E) - Us / a ** 3)  # (6.104) as printed
    correct = -a ** 3 / 2 * (-3 * a * E / a ** 5 * Us.dot(a * E) + Us / a ** 3)  # (6.103) at ξ = a e_ξ
    unit = {e3: sp.sqrt(1 - e1 ** 2 - e2 ** 2)}
    return {"chain": sp.simplify(chain - book), "direct": sp.simplify(direct - book.subs(sub)),
            "printed_bracket": sp.simplify((printed - grad).subs(unit)),
            "correct_bracket": sp.simplify((correct - grad).subs(unit)),
            "difference_chain": sp.simplify(chain - book), "difference_direct": sp.simplify(direct - book.subs(sub))}


def moving_sphere_state(theta, us: float = 1.0, dus_dt: float = 0.0, accel_angle: float = 0.0, a: float = 0.1,
                        rho: float = 1000.0, p_inf: float = 0.0) -> dict:
    """E9's scalar-callable state: surface pressure of a sphere moving at u_s (along x) with acceleration of magnitude
    du_s/dt at angle ``accel_angle`` to u_s (in the x–y plane), at the surface point at angle θ from u_s in that plane
    (6.105); its split into steady and acceleration parts; the surface speed (6.104); and the two force integrals
    (6.98): steady ≈ 0 (6.106), acceleration = −M du_s/dt (6.108).
    θ, accel_angle [rad]; us [m/s]; dus_dt [m/s²]; a [m]; ρ [kg/m³]. Returns dict(steady, acceleration, total [Pa],
    surface_speed [m/s], F_steady (3,), F_accel (3,) [N], M [kg]). Book: §6.9. Label: analytic."""
    th = float(theta)
    e = np.array([np.cos(th), np.sin(th), 0.0])
    U = np.array([float(us), 0.0, 0.0])
    Acc = float(dus_dt) * np.array([np.cos(accel_angle), np.sin(accel_angle), 0.0])
    sp_ = PF.moving_sphere_surface_pressure(e, U, Acc, a, rho, p_inf, split=True)
    ua = PF.moving_sphere_surface_velocity(e, U)
    Fs = PF.sphere_force_quadrature(lambda E: PF.moving_sphere_surface_pressure(E, U, np.zeros(3), a, rho, 0.0), a)
    Fa = PF.sphere_force_quadrature(lambda E: PF.moving_sphere_surface_pressure(E, np.zeros(3), Acc, a, rho, 0.0), a)
    return {"steady": sp_["steady"], "acceleration": sp_["acceleration"], "total": sp_["total"],
            "surface_speed": float(np.linalg.norm(ua)), "F_steady": Fs, "F_accel": Fa,
            "M": PF.added_mass_sphere(a, rho)}


def added_mass_state(a: float = 0.1, rho: float = 1000.0, us: float = 1.0, dus: float = 1.0, angle_deg: float = 0.0,
                     theta_s_deg: float = 0.0) -> dict:
    """E9's probe (all scalars): a sphere of radius a moving at speed u_s with acceleration du_s/dt = ``dus`` at
    ``angle_deg`` to its velocity; the surface point at θ_s from u_s in the plane of both. Returns the surface pressure
    parts (6.105) p_steady = ½ρu_s²((9/4)cos²θ_s − 5/4), p_accel = ρ(a/2) e_ξ·du_s/dt, p_total = their sum (gauge,
    p∞ = 0) [Pa]; the pressure-force integrals (6.98) by quadrature, F_steady = component along u_s (≈ 0, (6.106)) and
    F_accel = component along du_s/dt (= −M|du_s/dt|, (6.108)) [N]; M = 2πρa³/3 and M_energy (kinetic-energy route) [kg];
    cp_steady = p_steady/(½ρu_s²) = 1 − (9/4)sin²θ_s. Returns dict(p_steady, p_accel, p_total, F_steady, F_accel, M,
    M_energy, cp_steady). Book: §6.9 (6.105)–(6.108). Label: analytic."""
    th, be = np.radians(float(theta_s_deg)), np.radians(float(angle_deg))
    st_ = moving_sphere_state(th, us=us, dus_dt=dus, accel_angle=be, a=a, rho=rho, p_inf=0.0)
    eu = np.array([1.0, 0.0, 0.0])
    ea = np.array([np.cos(be), np.sin(be), 0.0])
    q = 0.5 * rho * float(us) ** 2
    return {"p_steady": float(st_["steady"]), "p_accel": float(st_["acceleration"]), "p_total": float(st_["total"]),
            "F_steady": float(st_["F_steady"] @ eu), "F_accel": float(st_["F_accel"] @ ea), "M": st_["M"],
            "M_energy": added_mass_by_energy(a, rho), "cp_steady": float(st_["steady"]) / q if q > 0 else float("nan")}


def added_mass_by_energy(a: float = 1.0, rho: float = 1000.0, U: float = 1.0, method: str = "surface",
                         epsrel: float = 1e-11) -> float:
    """Added mass of a sphere from the kinetic energy of the fluid it moves (Exercise 6.49 route, independent of the
    pressure route): T = ½ρ∫|∇φ|² dV over r > a with φ of (6.97), M = 2T/U² (= 2πρa³/3).
    ``method="surface"``: Green's identity T = −½ρ∮φ ∂φ/∂r dA on the sphere (Gauss–Legendre × trapezoid, φ and ∇φ from
    ``core.potential``); ``"volume"``: ``dblquad`` of ½ρ|∇φ|² in s = a/r ∈ (0, 1] and θ ∈ (0, π).
    Returns M [kg]. Book: §6.9 (6.108); Exercise 6.49. Label: analytic."""
    a, U = float(a), float(U)
    us = np.array([0.0, 0.0, U])
    if method == "surface":
        mu, wmu = np.polynomial.legendre.leggauss(24)
        ph = np.linspace(0.0, _TWO_PI, 16, endpoint=False)
        MU, PH = np.meshgrid(mu, ph, indexing="ij")
        S_ = np.sqrt(1.0 - MU ** 2)
        E = np.stack([S_ * np.cos(PH), S_ * np.sin(PH), MU]).reshape(3, -1)
        W = (wmu[:, None] * np.full(PH.shape, _TWO_PI / 16)).ravel() * a ** 2
        X = a * E
        phi = _F(PF.moving_sphere_potential(X, np.zeros(3), us, a))
        dphidr = np.sum(PF.moving_sphere_velocity(X, np.zeros(3), us, a) * E, axis=0)
        T = -0.5 * float(rho) * float(np.sum(phi * dphidr * W))
        return float(2.0 * T / U ** 2)
    if method != "volume":
        raise ValueError('method must be "surface" or "volume"')

    def integrand(th, s):
        r = a / s
        x = np.array([r * np.sin(th), 0.0, r * np.cos(th)])
        u = PF.moving_sphere_velocity(x, np.zeros(3), us, a)
        return float(np.dot(u, u)) * _TWO_PI * np.sin(th) * a ** 3 / s ** 4  # r² dr = a³ ds/s⁴

    val = integrate.dblquad(integrand, 0.0, 1.0, 0.0, np.pi, epsabs=0.0, epsrel=epsrel)[0]
    T = 0.5 * float(rho) * val
    return float(2.0 * T / U ** 2)


def cylinder_added_mass(a: float = 1.0, rho: float = 1000.0, method: str = "closed") -> float:
    """Added mass per unit length of a circular cylinder M′ = ρπa² (= the displaced mass per length; Exercise 6.45):
    ``"closed"`` or ``"energy"`` (½ρ∫|∇φ|²dA of w = −Ua²/z, the moving cylinder, by ``quad``). [kg/m].
    Book: §6.9 (Exercise 6.45). Label: analytic."""
    a, rho = float(a), float(rho)
    if method == "closed":
        return float(rho * np.pi * a ** 2)
    if method == "energy":
        # |dw/dz|² = U²a⁴/r⁴ for w = −Ua²/z; T = ½ρ∫∫ U²a⁴/r⁴ r dr dθ with U = 1
        val = quad(lambda r: a ** 4 / r ** 3 * _TWO_PI, a, np.inf, epsabs=0.0, epsrel=1e-12)[0]
        return float(rho * val)  # M = 2T/U²
    raise ValueError('method must be "closed" or "energy"')


def sphere_motion(m: float, a: float, rho: float = 1000.0, F_E=None, t_eval=None, g: float | None = None,
                  mode: str | None = None, amplitude: float = 0.05, omega: float = 2.0, added_mass: bool = True,
                  u0: float = 0.0, z0: float = 0.0) -> dict:
    """Newton's law for a submerged sphere with added mass (6.109) along z (z up): (m + M) du/dt = F_E(t) + (ρV − m)g
    (buoyancy − weight, when ``g`` is given), M = 2πρa³/3 (``added_mass=False`` drops M for comparison).

    modes: None or "force" (F_E constant or callable F_E(t) [N], default 0); "ball" (g defaults to 9.81); "bubble" (m defaults to 0 —
    initial acceleration 2g, not ∞; Exercise 6.46); "oscillating" (prescribed z_s = amplitude·sin(ωt): returns the
    fluid force F_s = −M du/dt and the needed F_E = (m + M)du/dt).
    Parameters: m [kg] (≥ 0), a [m], ρ [kg/m³], t_eval [s] (default 0…1 s, 201 points), u0 [m/s], z0 [m].
    Returns dict(t, x (= z, position along the motion [m]), u, du_dt (= dudt), F_s = −M du/dt [N], M, a0 (initial
    acceleration), work (∫F_net u dt [J]), kinetic = ½(m + M)(u² − u0²) [J]). ``solve_ivp`` RK45, rtol 1e-10. Book: §6.9 (6.109).
    Validation (planned): V1 constant F_E ⇒ u = F_E t/(m + M); bubble 2g; V4 work = kinetic. Label: analytic, conserved.
    """
    a = float(a)
    V = 4.0 / 3.0 * np.pi * a ** 3
    M = PF.added_mass_sphere(a, rho) if added_mass else 0.0
    t_eval = np.linspace(0.0, 1.0, 201) if t_eval is None else _F(t_eval)
    if mode == "bubble":
        m = 0.0 if m is None else float(m)
        g = G if g is None else g
    elif mode == "ball":
        g = G if g is None else g
    m = float(m)
    if m + M <= 0:
        raise ValueError("m + M must be > 0 (a massless body without added mass has no inertia)")
    if mode == "oscillating":
        z = z0 + amplitude * np.sin(omega * t_eval)
        u = amplitude * omega * np.cos(omega * t_eval)
        du = -amplitude * omega ** 2 * np.sin(omega * t_eval)
        return {"t": t_eval, "x": z, "z": z, "u": u, "du_dt": du, "dudt": du, "M": M, "a0": float(du[0]),
                "F_s": -M * du, "F_E": (m + M) * du, "work": None, "kinetic": 0.5 * (m + M) * (u ** 2 - u[0] ** 2)}
    if mode not in (None, "force", "ball", "bubble"):
        raise ValueError('mode must be None, "force", "ball", "bubble" or "oscillating"')
    F_E = 0.0 if F_E is None else F_E
    Ffn = F_E if callable(F_E) else (lambda t, c=float(F_E): c)
    grav = (rho * V - m) * float(g) if g is not None else 0.0

    def rhs(t, s):
        F = float(Ffn(t)) + grav
        return [s[1], F / (m + M), F * s[1]]

    sol = solve_ivp(rhs, (t_eval[0], t_eval[-1]), [z0, u0, 0.0], t_eval=t_eval, rtol=1e-10, atol=1e-12)
    if not sol.success:
        raise RuntimeError(sol.message)
    z, u, W = sol.y
    du = np.array([(float(Ffn(t)) + grav) / (m + M) for t in sol.t])
    return {"t": sol.t, "x": z, "z": z, "u": u, "du_dt": du, "dudt": du, "M": M, "a0": float(du[0]), "F_s": -M * du,
            "work": W, "kinetic": 0.5 * (m + M) * (u ** 2 - u0 ** 2)}


def rayleigh_collapse_time(R0: float = 1.0, rho: float = 1000.0, dp: float = 1.0e5, method: str = "closed") -> float:
    """Collapse time of an empty spherical cavity (Rayleigh; Exercise 6.44): Ṙ² = (2Δp/3ρ)(R₀³/R³ − 1) ⇒
    t_c = R₀√(3ρ/2Δp) ∫₀¹ x^{3/2}(1 − x³)^{−1/2} dx = R₀√(3ρ/2Δp)·B(5/6, 1/2)/3 ≈ 0.91468 R₀√(ρ/Δp).
    ``method``: "closed" (Beta function) or "quad" (the integral numerically). R₀ [m], ρ [kg/m³], Δp [Pa] → t_c [s].
    Book: §6.9 (Exercise 6.44). Validation (planned): V5 Wikipedia "Rayleigh–Plesset equation" 0.91468. Label: analytic."""
    pre = float(R0) * np.sqrt(3.0 * float(rho) / (2.0 * float(dp)))
    if method == "closed":
        return float(pre * beta_fn(5.0 / 6.0, 0.5) / 3.0)
    if method == "quad":
        f = lambda x: x ** 1.5 / np.sqrt(1.0 + x + x ** 2)  # noqa: E731  (1 − x³) = (1 − x)(1 + x + x²)
        val = quad(f, 0.0, 1.0, weight="alg", wvar=(0.0, -0.5), epsabs=0.0, epsrel=1e-12)[0]
        return float(pre * val)
    raise ValueError('method must be "closed" or "quad"')
