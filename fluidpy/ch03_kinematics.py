"""Chapter 3 — Kinematics: the chapter's worked examples, the flows its figures and explainers use, and re-exports of
the reusable primitives the chapter introduced in ``fluidpy.core`` (kinematics, coords, vortices, transport).

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e (2012), Ch. 3, §§3.1–3.6, Eqs. (3.1)–(3.35), Examples 3.1–3.2.
Every equation was transcribed from the rendered page images (chapters/pages/ch03/p096–p116).

Where the physics lives
-----------------------
* ``core.kinematics`` — (3.1)–(3.6) material derivative, (3.7)–(3.8) streamlines/path lines/streak lines, (3.9)
  Galilean invariance, (3.10)–(3.21) strain, rotation, volume change, principal axes; the linear-flow helpers moved
  from ch02.
* ``core.coords`` — Fig. 3.3: polar, cylindrical, spherical coordinates and velocity components.
* ``core.vortices`` — (3.22)–(3.29): solid body, line vortex, Rankine, Gaussian; polar vorticity; circulation.
* ``core.transport`` — (3.30)–(3.35): Leibniz, control volumes, the Reynolds transport theorem.
* this module — Examples 3.1 and 3.2, Fig. 3.1's section average, Fig. 3.2's cylinder in two frames, the explainer
  fields (``unsteady_flow_preset``, ``thermal_front``, ``frame_acceleration_terms``), §3.5's parallel shear flow and
  the polar-sector circulation (Fig. 3.16), Leibniz cases.

Book typos handled (analysis §9): (3.6) prints |u| ∂/∂s without F — implemented |u| ∂F/∂s; Example 3.2's "b = 0 on
the base" — only b·n = 0 is true there (implemented with the real radial b); the "[?]" in Example 3.2's integrand is
the factor z; Fig. 3.16's C is the B of (3.25). Book-quoted numbers live only in the git-ignored
``tests/book_values_ch03.json``.
"""
from __future__ import annotations

from functools import lru_cache
from typing import Callable

import numpy as np
import sympy as sp
from scipy.integrate import dblquad, quad

from .core._util import as_scalar_if_0d
from .core.coords import (  # noqa: F401
    cartesian_components, cartesian_from_cylindrical, cartesian_from_polar_coords, cartesian_from_spherical,
    cylindrical_from_cartesian, polar_from_cartesian, spherical_from_cartesian, unit_vectors_cylindrical,
    unit_vectors_polar, unit_vectors_spherical, velocity_components,
)
from .core.index_notation import expand_indices, expand_indices_str  # noqa: F401
from .core.integral_theorems import (  # noqa: F401
    circulation, curl_flux, divergence_theorem_sphere, planar_disc, planar_loop, rectangle_loop, stokes_theorem_check,
)
from .core.kinematics import (  # noqa: F401
    VELOCITY_GRADIENT_PRESETS, AccelerationTerms, MaterialDerivativeTerms, RelativeVelocitySplit, acceleration,
    as_coord_field, deform_circle, deform_sphere, deform_square, element_rotation_rate, from_coord_field,
    galilean_transform, lagrangian_to_eulerian, lagrangian_velocity_acceleration, lagrangian_velocity_acceleration_sym,
    linear_flow_map, linear_strain_rate, material_derivative, material_derivative_sym, material_derivative_terms,
    material_line_angle, material_line_rotation_rate, material_volume_ratio, measured_strain_rates, pathline,
    perpendicular_pair_rotation_rate, principal_strain_rates, relative_velocity, relative_velocity_split,
    shear_strain_rate, strain_ellipse_axes, strain_velocity_principal, streakline, streamline, streamwise_derivative,
    velocity_gradient_at, velocity_gradient_preset, velocity_potential_2d, volumetric_strain_rate, vorticity,
    vorticity_from_gradient, vorticity_in_rotating_frame,
)
from .core.operators import curl, divergence, is_irrotational, vector_gradient  # noqa: F401
from .core.tensors import (  # noqa: F401
    antisymmetric_from_vector, antisymmetric_part, cross, levi_civita, principal_axes, rotation_matrix_2d,
    rotation_tensor, strain_rate_tensor, trace, transform_tensor, vector_from_antisymmetric,
)
from .core.transport import (  # noqa: F401
    ControlVolume, GrowingCone, GrowingCylinder, GrowingSphere, LeibnizTerms, MaterialVolumeRate, MovingBox,
    MovingEllipse2D, RTTCheck, RTTTerms, leibniz_check, leibniz_terms, material_volume_rate, reynolds_transport,
    rtt_check, rtt_ellipse_2d, surface_flux_term, swept_terms, swept_terms_sphere, swept_volume_integral,
    volume_integral, volume_integral_rate_fd, volume_rate_term,
)
from .core.vortices import (  # noqa: F401
    VORTEX_KINDS, circulation_circle, gaussian_vortex, gaussian_vortex_max_radius, line_vortex, mean_vorticity_in_disc,
    polar_vorticity_z, polar_vorticity_z_sym, rankine_vortex, solid_body_rotation, vortex_profile,
    vortex_velocity_field,
)

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d


# ======================================================================================================================
# §3.1 One-, two- and three-dimensional flows; the cylinder in two frames (Figs. 3.1, 3.2)
# ======================================================================================================================
def cross_section_average(u_fn: Callable, R: float, z: float = 0.0) -> float:
    """Area average of an axisymmetric velocity profile over a pipe cross-section: ū(z) = (1/πR²)∫₀^R u(r, z) 2πr dr.

    Book: §3.1, Fig. 3.1c–d — a flow that is really 2-D (u varies across the pipe and along it) treated as 1-D by
    averaging over the cross-section (our formula for the book's words). Poiseuille u = U(1 − r²/R²) averages to U/2.

    Parameters
    ----------
    u_fn : callable u(r, z) [m/s];  R : pipe radius [m];  z : axial station [m]

    Returns
    -------
    ū [m/s] (``quad``, its error estimate asserted < 1e-10·max(1, |ū|)).

    Validation: V1 Poiseuille → U/2, uniform → U (1e-12). Label: analytic.
    """
    val, err = quad(lambda r: u_fn(r, z) * 2.0 * np.pi * r, 0.0, float(R), epsabs=1e-13, epsrel=1e-13)
    if err > 1e-10 * max(1.0, abs(val)):
        raise RuntimeError(f"quad error estimate too large: {err:g}")
    return float(val / (np.pi * float(R) ** 2))  # ū = (1/A) ∫_A u dA


def pipe_profile(r, z=0.0, U_mean: float = 1.0, R: float = 1.0, n0: float = 20.0, L_e: float = 1.0):
    """A developing pipe profile with the same mean at every station (our model for Fig. 3.1d):
    u(r, z) = U_mean (n + 2)/n (1 − (r/R)ⁿ) with n(z) = 2 + (n0 − 2) e^{−z/L_e}.

    Near the inlet (large n) the profile is almost flat; far downstream (n → 2) it is Poiseuille, u_max = 2U_mean.
    The area average is U_mean for every n (mass conservation), so the 1-D description ū(z) = U_mean is exact while
    the 2-D profile changes. Units: r, z, R, L_e [m]; U_mean [m/s]. Not a solution of the equations of motion (Ch. 8
    gives the real entrance flow) — a teaching profile. Label: analytic.
    Book: §3.1, Fig. 3.1c–d (our profile).
    """
    r_, z_ = _F(r), _F(z)
    n = 2.0 + (float(n0) - 2.0) * np.exp(-z_ / float(L_e))
    return _S(float(U_mean) * (n + 2.0) / n * (1.0 - (np.abs(r_) / float(R)) ** n))


def _cylinder_body(xr, yr, U, a):
    """Body-frame velocity of the uniform stream U e_x past a circular cylinder of radius a (uniform stream + doublet)."""
    r2 = xr ** 2 + yr ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        u = U * (1.0 - a ** 2 * (xr ** 2 - yr ** 2) / r2 ** 2)
        v = -U * a ** 2 * 2.0 * xr * yr / r2 ** 2
    inside = r2 < a ** 2
    return np.where(inside, np.nan, u), np.where(inside, np.nan, v)


def cylinder_flow(x, y, U: float, a: float, frame: str = "body", t: float = 0.0, U_frame: float | None = None):
    """Ideal (potential) flow past a circular cylinder of radius a seen by different observers (Fig. 3.2).

    Book: §3.1 and §3.3, Fig. 3.2 — (a) the cylinder at rest in a stream U (steady), (b) the cylinder moving through
    fluid at rest (unsteady); "these flow fields only differ by a Galilean transformation", u = U + u′. The flow
    itself is Ch. 6 §6.3 (uniform stream + doublet): in the body frame
    u = U[1 − a²(x² − y²)/r⁴], v = −2U a² xy/r⁴.

    Frames: ``"body"`` — observer fixed to the cylinder (fluid at infinity moves at +U e_x, steady); ``"fluid"`` —
    observer at rest with the far fluid (the cylinder moves at −U e_x, centre at x_c = −Ut: unsteady). A number
    ``U_frame`` ∈ [0, U] overrides ``frame``: the observer moves at −U_frame e_x relative to the far fluid, so
    U_frame = 0 is the fluid frame and U_frame = U the body frame (E3's slider). Observer-frame velocity at observer
    position (x, y): u_obs = u_body(x − x_c, y) − U + U_frame with x_c = (U_frame − U) t.

    Parameters
    ----------
    x, y : observer-frame coordinates [m] (floats or arrays);  U : free-stream speed [m/s];  a : radius [m];
    frame : "body" or "fluid";  t : time [s];  U_frame : observer speed as above [m/s]

    Returns
    -------
    (u, v) [m/s]; NaN inside the cylinder.

    Validation: V1 u·n = 0 on r = a; far field → (U_frame − U + U, 0); V2 ∇·u = 0 and ∇×u = 0 (sympy); frames differ
    by exactly U. Label: analytic, symbolic.
    """
    if U_frame is None:
        if frame not in ("body", "fluid"):
            raise ValueError('frame must be "body" or "fluid" (or pass U_frame)')
        U_frame = float(U) if frame == "body" else 0.0
    xc = (float(U_frame) - float(U)) * float(t)  # cylinder centre in the observer frame
    ub, vb = _cylinder_body(_F(x) - xc, _F(y), float(U), float(a))
    return _S(ub - float(U) + float(U_frame)), _S(vb)  # u_obs = u_body − U + U_frame  (u = U + u′)


def cylinder_streamfunction(x, y, U: float, a: float):
    """Body-frame stream function ψ = U y (1 − a²/r²) of the cylinder flow (u = ∂ψ/∂y, v = −∂ψ/∂x; Ch. 6 §6.3).
    Streamlines are its contours; used to check ``core.kinematics.streamline`` (ψ constant along it). NaN inside.
    Units [m²/s]. Label: analytic.
    Book: §3.1/§3.3, Fig. 3.2 (flow from Ch. 6 §6.3).
    """
    x_, y_ = _F(x), _F(y)
    r2 = x_ ** 2 + y_ ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        psi = float(U) * y_ * (1.0 - float(a) ** 2 / r2)
    return _S(np.where(r2 < float(a) ** 2, np.nan, psi))


def cylinder_velocity_field(U: float, a: float, frame: str = "body", U_frame: float | None = None) -> Callable:
    """``cylinder_flow`` as a kinematics field u(x, t) (x (2,) or (2, N)) for ``streamline``, ``pathline``,
    ``acceleration`` … Label: analytic.
    Book: §3.1/§3.3, Fig. 3.2.
    """
    def u(x, t=0.0):
        x_ = _F(x)
        uu, vv = cylinder_flow(x_[0], x_[1], U, a, frame, t, U_frame)
        return np.stack([_F(uu), _F(vv)])
    return u


def frame_acceleration_terms(x: float, y: float, U: float, a: float, U_frame: float, t: float = 0.0,
                             h: float = 1e-5) -> dict:
    """Velocity, local and advective acceleration at the observer-frame point (x, y) near the cylinder, for an
    observer moving at −U_frame e_x relative to the still far fluid (E3's term bars; Fig. 3.2 + Eq. (3.9)).

    Book: §3.3, Eq. (3.9) and the Fig. 3.2 discussion: ∂u/∂t + (u·∇)u is the same in every Galilean frame; in the
    body frame (U_frame = U) the flow is steady and the acceleration is all advective, in the fluid frame (U_frame = 0)
    both terms are present. The observer sees u_obs(X, t) = u_body(X + (U − U_frame)t, y) + (U_frame − U) e_x.

    Parameters
    ----------
    x, y : observer-frame position [m] (at t = 0 the cylinder centre is at the origin for every observer, so (x, y)
        is the same physical point for every U_frame — the comparison of (3.9));  U, a : free-stream speed [m/s],
        radius [m];  U_frame : observer speed [m/s] (0 = fluid frame, U = body frame);  t : time [s];  h : step [m]

    Returns
    -------
    dict with (2,) arrays ``u`` [m/s], ``local``, ``advective``, ``total`` [m/s²], the same as floats (``u_x``,
    ``u_y``, ``local_x``, ``local_y``, ``advective_x`` …) and ``steady`` (True in the body frame).
    Worked number (C05): (x, y) = (0, 1.5), U = a = 1, fluid frame → local_y = −0.592593 m/s².

    Validation: V1 total independent of U_frame at t = 0 to 1e-8 (Eq. 3.9); body frame: local = 0. Label: analytic.
    """
    field = cylinder_velocity_field(U, a, U_frame=U_frame)
    P = np.array([float(x), float(y)])
    a_, loc, adv = acceleration(field, P, float(t), h=h)  # Eq. (3.9) split in this frame
    uv = _F(field(P, t))
    return {"u": uv, "local": _F(loc), "advective": _F(adv), "total": _F(a_),
            "u_x": float(uv[0]), "u_y": float(uv[1]), "local_x": float(loc[0]), "local_y": float(loc[1]),
            "advective_x": float(adv[0]), "advective_y": float(adv[1]), "total_x": float(a_[0]),
            "total_y": float(a_[1]), "steady": bool(abs(float(U_frame) - float(U)) < 1e-15)}


# ======================================================================================================================
# §3.2 A Lagrangian map with a known Eulerian twin (D01); the thermal front (C02, E2)
# ======================================================================================================================
def lagrangian_map_example(X, t, alpha: float, derivative: int = 0):
    """The uniform-stretching Lagrangian map x = X e^{αt} (label X = position at t_o = 0) and its time derivatives.

    Book: §3.2, Eqs. (3.1)–(3.2) — our worked case D01: u = dx/dt = αX e^{αt}, a = α²X e^{αt}; eliminating the label
    (X = x e^{−αt}) gives the Eulerian field u(x, t) = αx, whose material derivative αu = α²x equals a.

    Parameters
    ----------
    X : label [m];  t : time [s];  alpha : stretching rate α [1/s];  derivative : 0 (position), 1 (u), 2 (a)

    Returns
    -------
    αᵏ X e^{αt} [m/sᵏ].  Worked number (X = 2 m, α = 0.5 s⁻¹, t = 1 s): 3.297 m, 1.649 m/s, 0.824 m/s².

    Validation: V1 = ``lagrangian_velocity_acceleration`` of the map; V2 ``lagrangian_to_eulerian`` gives u = αx.
    Label: analytic, symbolic.
    """
    k = int(derivative)
    return _S(float(alpha) ** k * _F(X) * np.exp(float(alpha) * _F(t)))  # dᵏ/dtᵏ (X e^{αt})


def thermal_front(x, y, t, grad_K_per_m: float, heating_K_per_s: float = 0.0, front_speed: float = 0.0,
                  width_m: float | None = None, T0: float = 288.15):
    """Temperature field of E2 and C02's worked number: a north–south gradient (or a tanh front) that moves north at
    ``front_speed`` and warms uniformly.

    T = T0 + H t − G (y − c t) (``width_m`` None) or T = T0 + H t − G w tanh((y − c t)/w) (a front of width w whose
    largest gradient is G), with G = ``grad_K_per_m`` > 0 meaning **colder to the north** (∂T/∂y = −G), H the heating
    rate and c the northward speed of the pattern.

    Book: §3.2, Eqs. (3.4)–(3.5) — our field for the local/advective split (N10: "advection" of heat by the wind).
    Worked number: G = 1e-5 K/m (1 K per 100 km), southerly wind v = 10 m/s, pattern carried by the wind (c = 10 m/s):
    v ∂T/∂y = −1e-4 K/s, ∂T/∂t = +1e-4 K/s (0.36 K/h at a station), DT/Dt = 0 for the air parcel.

    Parameters
    ----------
    x, y [m], t [s];  grad_K_per_m G [K/m];  heating_K_per_s H [K/s];  front_speed c [m/s];  width_m w [m] or None;
    T0 [K]

    Returns
    -------
    T [K] (scalar-callable; x enters only through broadcasting).  Exact terms: :func:`thermal_front_terms`.

    Validation: V1 ``material_derivative_terms`` on this field equals the analytic split. Label: analytic.
    """
    s_ = _F(y) - float(front_speed) * _F(t) + 0.0 * _F(x)
    G, H = float(grad_K_per_m), float(heating_K_per_s)
    prof = G * s_ if width_m is None else G * float(width_m) * np.tanh(s_ / float(width_m))
    return _S(float(T0) + H * _F(t) - prof)


def thermal_front_terms(x, y, t, u: float, v: float, grad_K_per_m: float, heating_K_per_s: float = 0.0,
                        front_speed: float = 0.0, width_m: float | None = None, T0: float = 288.15) -> dict:
    """The terms of Eq. (3.5) for :func:`thermal_front` in a uniform wind (u, v) [m/s], exactly.

    Returns dict: ``T`` [K], ``dTdx`` (0), ``dTdy`` (= −G, or −G sech²(s/w)) [K/m], ``local`` ∂T/∂t at the point,
    ``advective`` u ∂T/∂x + v ∂T/∂y, ``total`` DT/Dt [K/s], and ``regime``: "warm advection" when the advective term is
    negative (wind from warm to cold: a fixed station warms), "cold advection" when positive, "no advection" when
    u ⟂ ∇T. Book: §3.2, Eq. (3.5). Label: analytic.
    """
    s_ = _F(y) - float(front_speed) * _F(t) + 0.0 * _F(x)
    G, H = float(grad_K_per_m), float(heating_K_per_s)
    dTds = -G * np.ones_like(s_) if width_m is None else -G / np.cosh(s_ / float(width_m)) ** 2
    dTdy, dTdx = dTds, np.zeros_like(s_)
    loc = H - float(front_speed) * dTds  # ∂T/∂t at fixed (x, y)
    adv = float(u) * dTdx + float(v) * dTdy  # u·∇T
    a0 = float(np.ravel(adv)[0])
    regime = "no advection" if abs(a0) < 1e-15 else ("warm advection" if a0 < 0 else "cold advection")
    return {"T": thermal_front(x, y, t, G, H, front_speed, width_m, T0), "dTdx": _S(dTdx), "dTdy": _S(dTdy),
            "local": _S(loc), "advective": _S(adv), "total": _S(loc + adv), "regime": regime}  # Eq. (3.5)


# ======================================================================================================================
# §3.3 Flow lines: Example 3.1 and the explainer's unsteady fields (E1)
# ======================================================================================================================
def streamline_slope(u: Callable, x, t: float) -> float:
    """Slope dy/dx = v/u of the streamline through x at time t — the first equality of (3.7) in a plane.

    Book: §3.3, Eq. (3.7), Example 3.1 ("dy/dx = v/u"). Returns ±inf where u = 0. Units: dimensionless.
    Worked number: u = (1, 2) → 2. Label: analytic.
    """
    U = _F(u(_F(x), t))
    with np.errstate(divide="ignore"):
        return float(np.divide(U[1], U[0]))  # dy/dx = v/u


UNSTEADY_FLOW_PRESETS = ("ex31", "ex31_current", "steady", "steady_vortex", "rotating_strain")


def unsteady_flow_preset(name: str, x, y, t, xi0: float = 1.0, omega: float = 1.0, U0: float = 0.5,
                         beta: float = 0.0, s: float = 0.5, Omega: float = 0.5):
    """Plane velocity fields for E1 (flow lines) and the C03/C04 figures. Returns (u, v) [m/s].

    * ``"ex31"`` — Example 3.1: u = ωξ_o cos ωt, v = ωξ_o sin ωt (spatially uniform, unsteady; ξ_o [m], ω [rad/s]).
    * ``"ex31_current"`` — Example 3.1 plus a steady current U0 [m/s] along x (path lines become trochoids).
    * ``"steady"`` — Example 3.1's velocity frozen at the angle β [rad]: uniform (ωξ_o cos β, ωξ_o sin β), steady —
      all three lines are the same straight line.
    * ``"steady_vortex"`` — solid-body rotation u = Ω(−y, x) (steady, curved: the three lines coincide as circles).
    * ``"rotating_strain"`` — a pure strain of rate s [1/s] whose stretching axis turns at Ω [rad/s]:
      u = s[x cos 2Ωt + y sin 2Ωt], v = s[x sin 2Ωt − y cos 2Ωt] (principal axis at angle Ωt; Ω = 0 is a steady
      hyperbolic flow).

    Book: §3.3, Example 3.1 (first preset); the others are ours. Vectorised over x, y, t; scalar-callable.
    Validation: V1 "ex31" gives the Ex. 3.1 closed forms through ``streakline``/``pathline``; "rotating_strain" is
    traceless with S's first eigenvector at angle Ωt. Label: analytic.
    """
    x_, y_, t_ = np.broadcast_arrays(_F(x), _F(y), _F(t))
    if name == "ex31":
        u, v = omega * xi0 * np.cos(omega * t_), omega * xi0 * np.sin(omega * t_)  # Example 3.1
    elif name == "ex31_current":
        u, v = U0 + omega * xi0 * np.cos(omega * t_), omega * xi0 * np.sin(omega * t_)
    elif name == "steady":
        u, v = omega * xi0 * np.cos(beta) + 0.0 * t_, omega * xi0 * np.sin(beta) + 0.0 * t_
    elif name == "steady_vortex":
        u, v = -Omega * y_, Omega * x_
    elif name == "rotating_strain":
        c2, s2 = np.cos(2 * Omega * t_), np.sin(2 * Omega * t_)
        u, v = s * (c2 * x_ + s2 * y_), s * (s2 * x_ - c2 * y_)
    else:
        raise ValueError(f"unknown preset {name!r}; choose from {UNSTEADY_FLOW_PRESETS}")
    return _S(u + 0.0 * x_), _S(v + 0.0 * x_)


def preset_field(name: str, **params) -> Callable:
    """:func:`unsteady_flow_preset` as a kinematics field u(x, t) (x (2,) or (2, N); t float or (N,) array) — the
    callable twin that ``streamline``, ``pathline`` and ``streakline`` need. Label: analytic.
    Book: §3.3, Example 3.1 and E1's fields.
    """
    def u(x, t=0.0):
        x_ = _F(x)
        uu, vv = unsteady_flow_preset(name, x_[0], x_[1], t, **params)
        return np.stack([np.broadcast_to(_F(uu), x_[0].shape), np.broadcast_to(_F(vv), x_[0].shape)])
    return u


unsteady_flow_field = preset_field  # alias (curation §8 wording)


def example_3_1(t_prime: float, xi0: float = 1.0, omega: float = 1.0, n: int = 200) -> dict:
    """Example 3.1: the streamline, path line and streak line through the origin at t = t′ for
    u = ωξ_o cos ωt, v = ωξ_o sin ωt (closed forms).

    Book: §3.3, Example 3.1 and Fig. 3.7 — streamline y = x tan(ωt′); path line
    (x + ξ_o sin ωt′)² + (y − ξ_o cos ωt′)² = ξ_o², a circle of radius ξ_o centred at (−ξ_o sin ωt′, ξ_o cos ωt′);
    streak line (x − ξ_o sin ωt′)² + (y + ξ_o cos ωt′)² = ξ_o², centred at (ξ_o sin ωt′, −ξ_o cos ωt′); all three
    tangent at the origin (common slope tan ωt′). ⚠️ t′ (drawing instant) ≠ t_o (release time) ≠ ξ_o (amplitude).

    Parameters
    ----------
    t_prime : t′ [s];  xi0 : ξ_o [m];  omega : ω [rad/s];  n : points per curve

    Returns
    -------
    dict: ``streamline`` (2, n) segment of length 4ξ_o centred at the origin; ``pathline`` (2, n) = r(t; 0, t′) for
    t ∈ [t′, t′ + 2π/ω]; ``streakline`` (2, n) = positions at t′ of particles released at the origin at
    t_o ∈ [t′ − 2π/ω, t′]; ``path_center``, ``streak_center`` (2,) (aliases ``pathline_center``,
    ``streakline_center``), ``radius`` ξ_o, ``slope`` tan ωt′
    (±inf when cos ωt′ = 0), ``angle`` ωt′, and ``field`` (a kinematics callable u(x, t)).

    Validation: V1 the closed forms satisfy (3.7)/(3.8) (sympy residual 0); numeric ``streamline``, ``pathline``,
    ``streakline`` agree (1e-9); V6 private JSON ``example_3_1``. Label: analytic, symbolic.
    """
    w, xi, tp = float(omega), float(xi0), float(t_prime)
    ang = w * tp
    s = np.linspace(-2.0 * xi, 2.0 * xi, n)
    stream = np.stack([s * np.cos(ang), s * np.sin(ang)])  # y = x tan(ωt′), parametrised along its direction
    tt = np.linspace(tp, tp + 2.0 * np.pi / w, n)
    path = np.stack([xi * (np.sin(w * tt) - np.sin(ang)), xi * (-np.cos(w * tt) + np.cos(ang))])  # x = ξ[sin ωt − sin ωt′]
    to = np.linspace(tp - 2.0 * np.pi / w, tp, n)
    streak = np.stack([xi * (np.sin(ang) - np.sin(w * to)), xi * (-np.cos(ang) + np.cos(w * to))])  # at t = t′
    with np.errstate(divide="ignore"):
        slope = float(np.divide(np.sin(ang), np.cos(ang))) if abs(np.cos(ang)) > 1e-15 else float(np.sign(np.sin(ang)) * np.inf)
    return {"streamline": stream, "pathline": path, "streakline": streak,
            "pathline_center": np.array([-xi * np.sin(ang), xi * np.cos(ang)]),
            "streakline_center": np.array([xi * np.sin(ang), -xi * np.cos(ang)]),
            "path_center": np.array([-xi * np.sin(ang), xi * np.cos(ang)]),
            "streak_center": np.array([xi * np.sin(ang), -xi * np.cos(ang)]),
            "radius": xi, "slope": slope, "angle": ang, "field": preset_field("ex31", xi0=xi, omega=w)}


# ======================================================================================================================
# §3.4 Rigid motion; potential flow
# ======================================================================================================================
def rigid_body_velocity(U, Omega, x) -> np.ndarray:
    """Rigid-body velocity u = U + Ω × x (translation U [m/s] plus rotation at angular velocity Ω [rad/s]).

    Book: §3.4 ("S_ij is zero for any rigid body motion composed of translation at a spatially uniform velocity U and
    rotation at a constant rate Ω", Exercise 3.17 — written out as D10); the Ω × x term is the velocity used to read
    the second part of (3.19). x: (3,) or (3, N) [m].

    Validation: V1 its velocity gradient has S = 0 and vorticity 2Ω for 20 random (U, Ω). Label: analytic.
    """
    x_ = _F(x)
    U_ = _F(U).reshape((3,) + (1,) * (x_.ndim - 1))
    Om = np.broadcast_to(_F(Omega).reshape((3,) + (1,) * (x_.ndim - 1)), x_.shape)
    return U_ + np.cross(Om, x_, axis=0)  # u = U + Ω × x


def potential_velocity(phi_expr, coords) -> dict:
    """Velocity u = ∇φ of a velocity potential and its vorticity (zero) and divergence (∇²φ), symbolically.

    Book: §3.4, Eq. (3.17): irrotational flow ω = 0; "u can be written as the gradient of a scalar function φ".
    ⚠️ The converse needs a simply connected region (the line vortex (3.25) is irrotational off-axis but φ = Bθ is
    multivalued) — see ``velocity_potential_2d``.

    Returns (u_exprs, curl_expr): [∂φ/∂x_i] and the vorticity of ∇φ (scalar ω₃ in 2-D, list in 3-D; identically 0).
    Validation: V2 vorticity simplifies to 0 for any φ (ch02 D26 corollary). Label: symbolic.
    """
    phi = sp.sympify(phi_expr)
    X = list(coords)
    u = [sp.diff(phi, x) for x in X]  # u_i = ∂φ/∂x_i
    if len(X) == 2:
        vort = sp.simplify(sp.diff(u[1], X[0]) - sp.diff(u[0], X[1]))
    else:
        vort = [sp.simplify(sp.diff(u[2], X[1]) - sp.diff(u[1], X[2])), sp.simplify(sp.diff(u[0], X[2]) - sp.diff(u[2], X[0])),
                sp.simplify(sp.diff(u[1], X[0]) - sp.diff(u[0], X[1]))]
    return u, vort


# ======================================================================================================================
# §3.5 Kinematics of simple plane flows
# ======================================================================================================================
def parallel_shear_kinematics(gamma: float) -> dict:
    """Everything §3.5 says about the parallel shear flow u = (u₁(x₂), 0) with local shear γ = du₁/dx₂.

    Book: §3.5, Fig. 3.14: ω₃ = −γ (clockwise spin for γ > 0); line element AB (vertical) turns at −γ, BC
    (horizontal) at 0, average −γ/2 = ω₃/2 for any perpendicular pair (D14); S = [[0, γ/2], [γ/2, 0]] (pure shear
    for the aligned element ABCD); in the principal frame at 45° S̄ = diag(γ/2, −γ/2) (the element PQRS stretches and
    compresses without shear). ⚠️ γ = 2S₁₂ here, while ch02 Example 2.4 wrote Γ ≡ S₁₂.

    Parameters
    ----------
    gamma : γ [1/s]

    Returns
    -------
    dict: ``G`` (2 × 2), ``S``, ``R`` (= G − Gᵀ), ``omega3`` (−γ), ``spin`` (−γ/2), ``rate_AB`` (−γ), ``rate_BC`` (0),
    ``lam`` (ascending eigenvalues (−|γ|/2, |γ|/2)), ``axes`` (columns), ``S_bar`` (principal rates with the
    stretching axis first: diag(|γ|/2, −|γ|/2) — the book's diag(γ/2, −γ/2) for γ > 0), ``principal_angle_deg``
    (angle of the stretching axis from x₁, in (−90°, 90°]: +45° for γ > 0, −45° for γ < 0).

    Validation: V1 as listed; pair average −γ/2 for all θ; V6 private JSON ``shear_flow_3_5``. Label: analytic.
    """
    g = float(gamma)
    G = np.array([[0.0, g], [0.0, 0.0]])  # u₁ = γ x₂ locally
    S = strain_rate_tensor(G)  # [[0, γ/2], [γ/2, 0]]
    lam, B = principal_strain_rates(G)
    k = int(np.argmax(lam))  # stretching axis
    ang = float(np.degrees(np.arctan2(B[1, k], B[0, k])))
    ang = ang - 180.0 if ang > 90.0 else (ang + 180.0 if ang <= -90.0 else ang)
    return {"G": G, "S": S, "R": rotation_tensor(G), "omega3": float(vorticity_from_gradient(G)[2]),
            "spin": float(element_rotation_rate(G)[2]), "rate_AB": float(material_line_rotation_rate(G, np.pi / 2)),
            "rate_BC": float(material_line_rotation_rate(G, 0.0)), "lam": lam, "axes": B,
            "S_bar": np.diag([lam[k], lam[1 - k]]), "principal_angle_deg": ang}


def annular_sector_circulation(u_theta, r: float, dr: float, dtheta: float, u_r: Callable | None = None,
                               theta0: float = 0.0, legs: bool = False, **kw):
    """Circulation around the small polar sector r ≤ r′ ≤ r + Δr, θ0 ≤ θ ≤ θ0 + Δθ, counter-clockwise.

    Book: §3.5, Fig. 3.16: Γ_ABCD = {∫_AB + ∫_BC + ∫_CD + ∫_DA} u·ds = −[u_θ r]_r Δθ + [u_θ r]_{r+Δr} Δθ = 0 for the
    line vortex (the radial legs contribute nothing when u_r = 0; u_θ r = B is constant). Divided by the area
    r̄ Δr Δθ it tends to ω_z of (3.23) — the construction of D17. (Fig. 3.16 labels u_θ = C/r; the text's B.)

    Parameters
    ----------
    u_theta, u_r : callables f(r) or f(r, θ) [m/s] (u_r = None → 0); ``u_theta`` may be a kind name of
    ``vortex_profile`` with keywords ``Gamma``, ``sigma`` (``**kw``);  r : inner radius [m];  dr : Δr [m];
    dtheta : Δθ [rad];  theta0 : start angle [rad];  legs : also return the four leg integrals

    Returns
    -------
    Γ [m²/s] (float), or dict(total, outer, inner, radial_start, radial_end, area) if ``legs``.

    Validation: V1 line vortex 0 (1e-12); Rankine core ω_z × area; V3 Γ/area → ω_z at order 2 in size. Label: analytic,
    converged.
    """
    from .core.vortices import _polar_fn, _profile_fn

    ut, ur = _polar_fn(_profile_fn(u_theta, **kw)), _polar_fn(u_r)
    r1, r2, t1, t2 = float(r), float(r) + float(dr), float(theta0), float(theta0) + float(dtheta)
    q = lambda f, a, b: quad(f, a, b, epsabs=1e-14, epsrel=1e-13)[0]  # noqa: E731
    radial_start = q(lambda rr: float(ur(rr, t1)), r1, r2)  # along θ = θ0 outward: ∫ u_r dr
    outer = q(lambda th: float(ut(r2, th)) * r2, t1, t2)  # along r + Δr counter-clockwise: ∫ u_θ (r+Δr) dθ
    radial_end = -q(lambda rr: float(ur(rr, t2)), r1, r2)  # along θ0 + Δθ inward
    inner = -q(lambda th: float(ut(r1, th)) * r1, t1, t2)  # along r clockwise (against u_θ)
    total = radial_start + outer + radial_end + inner
    if legs:
        return {"total": total, "outer": outer, "inner": inner, "radial_start": radial_start, "radial_end": radial_end,
                "area": 0.5 * (r2 ** 2 - r1 ** 2) * float(dtheta)}
    return float(total)


# ======================================================================================================================
# §3.6 Leibniz cases and Example 3.2
# ======================================================================================================================
_LEIBNIZ_CASES = {
    # name: (F(x, t), a(t), b(t)) as sympy expressions in x, t
    "x2t": ("x**2*t", "t", "t**2"),
    "wave": ("sin(x - t)", "t/2", "2 + t"),
    "bump": ("exp(-(x - t)**2)", "0", "1 + t/2"),
}


@lru_cache(maxsize=None)
def _leibniz_case(case: str):
    case = {"power": "x2t"}.get(case, case)
    if case not in _LEIBNIZ_CASES:
        raise ValueError(f"case must be one of {tuple(_LEIBNIZ_CASES)}")
    x, t = sp.symbols("x t", real=True)
    Fe, ae, be = (sp.sympify(s, locals={"x": x, "t": t}) for s in _LEIBNIZ_CASES[case])
    I = sp.integrate(Fe, (x, ae, be))  # closed-form integral (sympy) — an independent route to d/dt ∫F dx
    rate = sp.simplify(sp.diff(I, t))
    lam = lambda e, args=(x, t): sp.lambdify(args, e, "numpy")  # noqa: E731
    return {"F": lam(Fe), "dFdt": lam(sp.diff(Fe, t)), "a": lam(ae, (t,)), "b": lam(be, (t,)),
            "adot": lam(sp.diff(ae, t), (t,)), "bdot": lam(sp.diff(be, t), (t,)), "rate": lam(rate, (t,)),
            "exprs": (Fe, ae, be, I, rate)}


def leibniz_example(t: float, case: str = "x2t") -> dict:
    """1-D cases for Leibniz's theorem (3.30), shared by N46's figure, D21's check and E7's Leibniz mode.

    Cases (ours): ``"x2t"`` F = x²t on [t, t²] (closed-form integral (t⁷ − t⁴)/3; alias "power"); ``"wave"``
    F = sin(x − t) on [t/2, 2 + t]; ``"bump"`` a moving Gaussian F = exp(−(x − t)²) on [0, 1 + t/2].

    Returns
    -------
    dict: the (3.30) terms ``interior``, ``upper``, ``lower`` (subtracted), ``total`` (from ``leibniz_terms``);
    ``exact`` = d/dt of the sympy closed-form integral (alias ``exact_rate``) — so ``total == exact`` is a test of
    (3.30); floats ``a``, ``b``, ``adot``, ``bdot`` at t; callables ``F(x, t)``, ``dFdt(x, t)``.

    Validation: V2 total = sympy rate (1e-10). Label: symbolic, analytic.
    Book: §3.6, Eq. (3.30), Fig. 3.17.
    """
    c = _leibniz_case(case)
    tt = float(t)
    a, b, ad, bd = (float(c[k](tt)) for k in ("a", "b", "adot", "bdot"))
    terms = leibniz_terms(c["F"], c["dFdt"], a, b, ad, bd, tt)
    exact = float(c["rate"](tt))
    return {"interior": terms.interior, "upper": terms.upper, "lower": terms.lower, "total": terms.total,
            "exact": exact, "exact_rate": exact, "a": a, "b": b, "adot": ad, "bdot": bd, "F": c["F"], "dFdt": c["dFdt"]}


def example_3_2(h: float, r0: float, rdot: float, n: int = 16) -> dict:
    """Example 3.2: rate of volume increase of a fixed-height cone whose base radius grows, three ways.

    Book: §3.6, Example 3.2, Fig. 3.19: V = ⅓πhr², so dV/dt = ⅔πh r_o ṙ directly; by (3.35) with F = 1 only the
    surface term survives, dV/dt = ∫_{A*} b·n dA with b = (z/h)ṙ e_R and n = e_R cos θ − e_z sin θ on the side
    (tan θ = r_o/h), dA = z tan θ dφ dz/cos θ, giving (2πṙ tan θ/h)∫₀ʰ z² dz = ⅔πh²ṙ tan θ = ⅔πh r_o ṙ.
    ⚠️ On the base b·n = 0 (the base points move radially in its plane; the book writes b = 0); the "[?]" printed in
    the integrand is the factor z.

    Parameters
    ----------
    h : height [m];  r0 : base radius now [m];  rdot : ṙ [m/s];  n : quadrature nodes for the ``GrowingCone`` route

    Returns
    -------
    dict(direct, rtt_dblquad, rtt_cv, rtt_closed_form) [m³/s] — all equal.

    Validation: V1 the three agree (1e-10); V2 sympy of D33; V6 private JSON ``example_3_2``. Label: analytic,
    symbolic.
    """
    h_, r_, rd = float(h), float(r0), float(rdot)
    th = np.arctan2(r_, h_)
    direct = 2.0 / 3.0 * np.pi * h_ * r_ * rd  # d/dt (⅓πhr²)

    def integrand(phi, z):  # b·n dA on the side: (z/h) ṙ cos θ · z tan θ / cos θ
        return (z / h_) * rd * np.cos(th) * z * np.tan(th) / np.cos(th)
    rtt_dq, _ = dblquad(integrand, 0.0, h_, 0.0, 2.0 * np.pi, epsabs=1e-14, epsrel=1e-13)
    cone = GrowingCone(h=h_, r0=r_, rdot=rd)
    rtt_cv = surface_flux_term(lambda x, t: np.ones(x.shape[1]), cone, 0.0, n)  # (3.35) with F = 1
    return {"direct": float(direct), "rtt_dblquad": float(rtt_dq), "rtt_cv": float(rtt_cv),
            "rtt_closed_form": float(2.0 / 3.0 * np.pi * h_ ** 2 * rd * np.tan(th))}


RTT_FIELDS = ("uniform", "ramp", "warming", "carried")
CARRY_SPEED = 0.3  # speed [m/s] along x₁ at which the "carried" pattern moves


def rtt_field(name: str, dim: int = 3):
    """The scalar fields of E7 and C15, as (F, dFdt) callables.

    Book: §3.6 — our fields for Leibniz's theorem (3.30) and the Reynolds transport theorem (3.35),
    $\\frac{d}{dt}\\int_{V^*(t)} F(\\mathbf x,t)\\,dV = \\int_{V^*(t)} \\frac{\\partial F(\\mathbf x,t)}{\\partial t}\\,dV
    + \\int_{A^*(t)} F(\\mathbf x,t)\\,\\mathbf b\\cdot\\mathbf n\\,dA$.

    * ``"uniform"``  F = 1
    * ``"ramp"``     F = 1 + 0.5 x₁
    * ``"warming"``  F = 1 + 0.5 x₁ + 0.3 t
    * ``"carried"``  F = 1 + 0.5 (x₁ − 0.3 t): the ramp carried along x₁ at 0.3 m/s by a uniform flow u = (0.3, 0, 0),
      so ∂F/∂t = −0.15 and DF/Dt = ∂F/∂t + u·∇F = −0.15 + 0.3·0.5 = 0. A control volume translating with the flow
      (b = u) sees the two terms of (3.35) cancel: ∫∂F/∂t dV = −∮F b·n dA.

    Parameters
    ----------
    name : one of ``RTT_FIELDS``
    dim : 1 → F(x, t) of a position x [m] (float or array, the Leibniz mode); 2 or 3 → F(x, t) of (d, N) points [m].
          t [s]. F in arbitrary units per m³ (per m² in 2-D, per m in 1-D); dFdt in the same units per s.

    Returns
    -------
    (F, dFdt) : callables (x, t); scalar-callable (explainer parity).

    Assumptions: linear in x₁ and t, independent of x₂, x₃.
    Validation: V1 the fields are closed-form; "carried": the 2-D ellipse RTT total with b = (0.3, 0) is 0 (the two
    terms cancel), and ∂F/∂t + 0.3 ∂F/∂x₁ = 0. Label: analytic.
    """
    if name not in RTT_FIELDS:
        raise ValueError(f"name must be one of {RTT_FIELDS}")
    c1 = 0.0 if name == "uniform" else 0.5  # ∂F/∂x₁ [F/m]
    c2 = {"warming": 0.3, "carried": -0.5 * CARRY_SPEED}.get(name, 0.0)  # ∂F/∂t [F/s]

    def x1(x):
        x_ = _F(x)
        return x_ if dim == 1 else x_[0]

    def F(x, t):
        return _S(1.0 + c1 * x1(x) + c2 * float(t))  # "carried": 1 + 0.5(x₁ − 0.3t)

    def dFdt(x, t):
        return _S(c2 + 0.0 * x1(x))
    return F, dFdt


def flux_through_disc(u: Callable, center, normal, radius: float, t: float = 0.0, nr: int = 32,
                      ntheta: int = 64) -> float:
    """Volume flux ∫ u·n dA through a planar disc (midpoint rule on ``core.integral_theorems.planar_disc`` nodes).

    Book: §3.3, Fig. 3.6 (stream tube: no fluid crosses the tube wall, so the flux through two cross-sections of one
    tube is the same — N15's number). Parameters: u(x, t) (3-D) [m/s]; center (3,) [m]; normal (3,) (normalised);
    radius [m]; t [s]; nr, ntheta resolution. Returns [m³/s]. Label: analytic.
    """
    from .core.integral_theorems import planar_disc

    disc = planar_disc(_F(center), _F(normal), float(radius), nr, ntheta)
    U = _F(u(disc.points.T, t))
    return float(np.sum(np.einsum("ik,ki->k", U, disc.normals) * disc.dA))  # ∫ u·n dA


_NOT_EXPORTED = {"np", "sp", "quad", "dblquad", "lru_cache", "annotations", "Callable", "as_scalar_if_0d"}
