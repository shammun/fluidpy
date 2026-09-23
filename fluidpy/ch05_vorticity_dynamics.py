"""Chapter 5 — Vorticity dynamics: the two basic vortices with pressure and viscous stress, Kelvin's and Helmholtz's
theorems in numbers, baroclinic generation, stretching and tilting, the rotating fluid column, vortex pairs, walls,
rings and sheets — plus re-exports of the reusable primitives the chapter introduced in ``fluidpy.core``.

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e (2012), Ch. 5, §§5.1–5.8, Eqs. (5.1)–(5.33), Figs. 5.1–5.16,
Exercises 5.1–5.20 (test fields only). Every equation was transcribed from the rendered page images
(chapters/pages/ch05/p198–p222). The public names follow ``analysis/ch05_design.md`` Part C (the contract with the
notebook and the explainers).

Where the physics lives (all public names are re-exported here, so ``ch05.<name>`` reaches every callable)
--------------------------------------------------------------------------------------------------------
* ``core.vorticity`` (VD) — vortex lines and tubes (5.3)–(5.4), material circulation and Kelvin (5.8)–(5.11), the
  pressure torque (Fig. 5.6), Helmholtz's frozen-in lines, the vorticity equation (5.12)–(5.13) and (5.18)–(5.30)
  (numeric, sympy, presets), stretching/tilting (5.31)–(5.32), absolute circulation (5.33), E3's Kelvin scenarios.
* ``core.biot_savart`` (BS) — Poisson/Green (5.14) (sign corrected), Biot–Savart (5.16)–(5.17), segments and rings,
  point vortices with wall/circle/channel images (§5.7), vortex sheets (§5.8), E6/E8 presets.
* ``core.vortices`` (VX) — ch03's plane vortices + Burgers' and Hill's vortices, the Lamb–Oseen field, the stretched
  Gaussian vortex, the narrowing Gaussian vortex tube and its "broken" twin, a ring core, the ABC flow.
* ``core.integral_theorems.curl_theorem_box`` — Gauss in curl form (5.15).
* this module — the pressure fields of the rotating tank (5.6) and of the line vortex (5.7), the Rankine tornado, the
  viscous stress without net force, the rotating cylinder's torque and dissipation, Lamb–Oseen circulation (5.11),
  the lock exchange, Kelvin's four restrictions, the diffusing vortex sheet, the rotating-frame residual (5.20),
  uniform-strain stretching/tilting, the helix frame, the fluid column and the ring of air (5.33), vortex pairs and
  wall drift, ring dynamics (leap-frogging, ring toward a wall), sheet strength, convergence and roll-up, and the
  explainer scenario functions.

Conventions (analysis §9): **ω is the vorticity** — in (5.1) u_θ = ωr/2 the tank turns at ω/2 (ch03's
``solid_body_rotation(r, omega0)`` takes the rotation rate ω₀ = ω/2; :func:`solid_body_from_vorticity` converts); Γ
(``Gamma``) [m²/s] is always a circulation, γ (``gamma``) [m/s] the strength of a vortex sheet; σ is the viscous stress
(ch03's core radius is ``sigma_core``/``sigma_c`` here); counterclockwise-positive plane vorticity and circulation;
z up, g = 9.81 m/s² (the book's rounded value, ``core.thermo.G_BOOK``); SI units, nothing is non-dimensionalised;
latitudes in degrees only at interfaces named ``*_deg``.

Book slips handled in code and docstrings: (5.14) printed with −1/(4π) → +1/(4π) (``core.biot_savart``); Fig. 5.16's
caption u₁ − u₂ vs the text's u₂ − u₁ (:func:`vortex_sheet_strength` follows the text, counterclockwise positive, with a
switch); Fig. 5.2's "2ω" label (the tank turns at ω/2); "single valued" is not why ∮dp/ρ = 0 (barotropy is,
``core.vorticity.kelvin_force_terms``); "Exercise 5.8" after (5.14) → 5.9, "Exercise 5.11" after (5.33) → 5.10;
"hyperboloids of the second degree" in Fig. 5.3 are cubic surfaces (c − z)r² = const. Book-quoted numbers live only in
the git-ignored ``tests/book_values_ch05.json``.
"""
from __future__ import annotations

from typing import Callable

import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp
from scipy.linalg import expm
from scipy.optimize import brentq
from scipy.special import erf

from .core import _stencil as st
from .core import biot_savart as BS  # noqa: F401
from .core import curvilinear as CU  # noqa: F401
from .core import integral_theorems as IT  # noqa: F401
from .core import vortices as VX  # noqa: F401
from .core import vorticity as VD  # noqa: F401
from .core._util import as_scalar_if_0d
from .core.biot_savart import *  # noqa: F401,F403
from .core.constitutive import dissipation_rate
from .core.integral_theorems import CurlCheck, curl_theorem_box  # noqa: F401
from .core.kinematics import acceleration
from .core.rotating import OMEGA_EARTH, coriolis_parameter  # noqa: F401
from .core.thermo import G_BOOK
from .core.vortices import *  # noqa: F401,F403
from .core.vorticity import *  # noqa: F401,F403

G = G_BOOK  #: default g [m/s²] of this chapter's functions (the book's 9.81)
_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d


# ======================================================================================================================
# §5.1 the two basic vortices (5.1), (5.2) and their pressure (5.5)–(5.7)
# ======================================================================================================================
def solid_body_from_vorticity(r, omega):
    """Solid-body rotation produced by a uniform plane-normal vorticity ω, Eq. (5.1): u_θ = ωr/2.

    Book: §5.1, Eq. (5.1). ⚠️ ω is the **vorticity**; the fluid turns at ω/2 (ch03's (3.22) u_θ = ω₀r with ω₀ = ω/2) —
    a tank turning at 1 rad/s has ω = 2 s⁻¹. Delegates to ``core.vortices.solid_body_rotation(r, omega/2)``.
    Parameters: r [m] (≥ 0); omega ω [1/s]. Returns (u_theta [m/s], omega_z [1/s] = ω everywhere).
    Validation: V1 the numerical curl of the field is ω (not 2ω — passing ω as the rotation rate fails). Label: analytic.
    """
    u = VX.solid_body_rotation(r, 0.5 * float(omega))  # Eq. (5.1): u_θ = ωr/2  (rotation rate ω/2)
    return u, _S(np.full(np.shape(r), float(omega)))


def line_vortex_gamma(r, Gamma):
    """Ideal (irrotational) line vortex of circulation Γ, Eq. (5.2): u_θ = Γ/2πr.

    Book: §5.1, Eq. (5.2) (ch03 (3.25) with B = Γ/2π). The vorticity is a line delta: 0 for r > 0, infinite on the axis.
    Parameters: r [m]; Gamma [m²/s]. Returns (u_theta [m/s], omega_z [1/s]) — both NaN at r = 0, ω_z = 0 for r > 0.
    Validation: V1 ∮u·dx = Γ on every circle; ω_z = 0 for r > 0. Label: analytic.
    """
    u = VX.line_vortex(r, float(Gamma) / (2.0 * np.pi))  # Eq. (5.2): u_θ = Γ/2πr
    r_ = _F(r)
    return u, _S(np.where(r_ > 0.0, 0.0, np.nan))


def solid_body_pressure_gradients(r, z, omega, rho: float = 1000.0, g: float = G):
    """Radial and vertical momentum balance of solid-body rotation, Eqs. (5.5a, 5.5b): ∂p/∂r = ρu_θ²/r = ρω²r/4,
    ∂p/∂z = −ρg.

    Book: §5.1, (5.5a) −ρu_θ²/r = −∂p/∂r (the pressure gradient supplies the centripetal acceleration), (5.5b)
    0 = −∂p/∂z − ρg (hydrostatic). S = 0 in rigid rotation ⇒ τ = −pδ (no viscous stress), so Cauchy (4.24) reduces to
    Euler (4.41). Parameters: r, z [m]; omega ω [1/s]; rho [kg/m³]; g [m/s²]. Returns (dp_dr, dp_dz) [Pa/m].
    Validation: V2 sympy: the cylindrical Euler equations (``core.curvilinear.advective_acceleration``) give (5.5).
    Label: analytic, symbolic.
    """
    u, _ = solid_body_from_vorticity(r, omega)
    r_ = _F(r)
    with np.errstate(divide="ignore", invalid="ignore"):
        dpdr = np.where(r_ > 0.0, float(rho) * _F(u) ** 2 / np.where(r_ > 0.0, r_, 1.0), 0.0)  # (5.5a)
    dpdz = -float(rho) * float(g) + 0.0 * _F(z) + 0.0 * r_  # (5.5b)
    return _S(dpdr), _S(dpdz)


def solid_body_pressure(r, z, omega, rho: float = 1000.0, g: float = G, p_o: float = 0.0):
    """Pressure in a steadily rotating tank (solid-body rotation of vorticity ω), Eq. (5.6):
    p(r, z) − p_o = ⅛ρω²r² − ρgz, p_o at r = 0, z = 0.

    Book: §5.1, Eq. (5.6) (integrating (5.5a) gives ρω²r²/8 + f(z), (5.5b) gives −ρgz + g(r); consistent only for
    (5.6) — our D02). Isobars z = ω²r²/8g − (p − p_o)/ρg are paraboloids of revolution (Fig. 5.2); with the tank's
    rotation rate Ω_t = ω/2 the free surface rises Ω_t²r²/2g (Exercise 5.1's hint).
    Parameters: r, z [m]; omega ω [1/s]; rho [kg/m³]; g [m/s²]; p_o [Pa]. Returns p [Pa]. (0.1 m, ω = 10 s⁻¹ → 125 Pa.)
    Validation: V2 satisfies (5.5a, b); V1 isobar heights; V6 Exercise 5.1 via :func:`rotating_tank_free_surface`.
    Label: analytic, symbolic.
    """
    return _S(float(p_o) + float(rho) * float(omega) ** 2 * _F(r) ** 2 / 8.0 - float(rho) * float(g) * _F(z))  # (5.6)


def line_vortex_pressure(r, z, Gamma, rho: float = 1000.0, g: float = G, p_inf: float = 0.0):
    """Pressure of the ideal line vortex, Eq. (5.7): p(r, z) − p_∞ = −ρΓ²/(8π²r²) − ρgz (p_∞ far away at z = 0).

    Book: §5.1, Eq. (5.7) ((5.2) substituted in (5.5) and integrated). Isobars z = −Γ²/(8π²r²g) − (p − p_∞)/ρg form a
    funnel (Fig. 5.3; the book calls them "hyperboloids of revolution of the second degree", but (c − z)r² = const is a
    cubic surface). ½u_θ² + gz + p/ρ is the same everywhere (irrotational flow, Bernoulli (4.72)).
    Worked number: Γ = 1 m²/s, water (ρ = 1000), r = 0.1 m → deficit ρΓ²/8π²r² = 1266.5 Pa, funnel depth 0.129 m.
    Parameters: r (> 0), z [m]; Gamma [m²/s]; rho; g; p_inf [Pa]. Returns p [Pa].
    Validation: V1 (5.5a) residual with u_θ = Γ/2πr; B constant. Label: analytic.
    """
    return _S(float(p_inf) - float(rho) * float(Gamma) ** 2 / (8.0 * np.pi ** 2 * _F(r) ** 2)
              - float(rho) * float(g) * _F(z))  # Eq. (5.7)


def rankine_pressure(r, z, Gamma, a, rho: float = 1000.0, g: float = G, p_inf: float = 0.0):
    """Pressure of a Rankine vortex (core radius a, circulation Γ): the funnel (5.7) outside, a paraboloid (5.6) inside.

    Book: §5.1 (the rotating cylinder u_θ = ωa²/2r "precisely matches" a Rankine vortex with Γ = πa²ω; (3.28)) and
    Exercise 5.2 (tornado). Inside, ω_c = Γ/πa² and p − p_∞ = ρω_c²r²/8 − ρΓ²/(4π²a²) − ρgz; outside (5.7). p and ∂p/∂r
    are continuous at r = a, where p − p_∞ = −ρΓ²/(8π²a²) (half the axis deficit ρΓ²/4π²a²).
    Parameters: r, z [m]; Gamma [m²/s]; a [m]; rho; g; p_inf [Pa]. Returns p [Pa].
    Validation: V1 continuity of p and ∂p/∂r at a; equals ``core.bernoulli.rankine_vortex_pressure`` at z = 0; V6
    Exercise 5.2 (private). Label: analytic.
    """
    r_ = _F(r)
    wc = float(Gamma) / (np.pi * float(a) ** 2)
    inside = float(rho) * wc ** 2 * r_ ** 2 / 8.0 - float(rho) * float(Gamma) ** 2 / (4.0 * np.pi ** 2 * float(a) ** 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        outside = -float(rho) * float(Gamma) ** 2 / (8.0 * np.pi ** 2 * np.where(r_ > 0, r_, 1.0) ** 2)
    return _S(float(p_inf) + np.where(r_ <= float(a), inside, outside) - float(rho) * float(g) * _F(z))


def tornado_circulation(p_gauge_edge: float, a: float, rho: float = 1.184) -> float:
    """Circulation of a Rankine tornado from the gauge pressure at its core edge: p − p_∞ = −ρΓ²/(8π²a²) ⇒
    Γ = 2πa√(2|p − p_∞|/ρ) [m²/s] (Exercise 5.2(a) route: Bernoulli between infinity and r = a).
    Parameters: p_gauge_edge [Pa] (< 0); a core radius [m]; rho air density [kg/m³] (default 25 °C, 1.184).
    Book: §5.1, Eq. (5.7); Exercise 5.2. Label: analytic."""
    return float(2.0 * np.pi * a * np.sqrt(2.0 * abs(p_gauge_edge) / rho))


def rankine_radius_at_pressure(p_gauge, Gamma: float, a: float, rho: float = 1.184):
    """Radius where a Rankine vortex's gauge pressure (z = 0) equals ``p_gauge`` (< 0) [m]: outside the core
    r = (Γ/2π)√(ρ/(2|p|)), inside from the paraboloid; NaN if |p| exceeds the axis deficit ρΓ²/4π²a².
    Book: §5.1, (5.6)–(5.7); Exercise 5.2(b) (a tornado passing at speed U: time = Δr/U). Label: analytic."""
    p_ = _F(p_gauge)
    edge = rho * Gamma ** 2 / (8.0 * np.pi ** 2 * a ** 2)
    axis = 2.0 * edge
    with np.errstate(divide="ignore", invalid="ignore"):
        out = Gamma / (2.0 * np.pi) * np.sqrt(rho / (2.0 * np.abs(p_)))
        wc = Gamma / (np.pi * a ** 2)
        ins = np.sqrt(np.clip((axis - np.abs(p_)) * 8.0 / (rho * wc ** 2), 0.0, None))
    return _S(np.where(np.abs(p_) <= edge, out, np.where(np.abs(p_) <= axis, ins, np.nan)))


def isobar_height(r, param, dp_over_rho_g: float = 0.0, kind: str = "solid", g: float = G, a: float | None = None):
    """Height z(r) of the isobar (p − p_ref)/(ρg) = ``dp_over_rho_g`` [m] (Figs. 5.2, 5.3).

    ``kind="solid"`` (param = ω [1/s]): z = ω²r²/8g − Δp/ρg (paraboloid, (5.6)); ``"line"`` (param = Γ [m²/s]):
    z = −Γ²/(8π²r²g) − Δp/ρg (funnel, (5.7)); ``"rankine"`` (param = Γ, core ``a``): paraboloid inside, funnel outside.
    Book: §5.1 after (5.6) and (5.7). Parameters: r [m]; dp_over_rho_g [m]; g [m/s²]. Returns z [m].
    Validation: V1 isobar_height(0.1, 1.0, kind="line") = −0.129104 m. Label: analytic.
    """
    r_ = _F(r)
    if kind == "solid":
        return _S(float(param) ** 2 * r_ ** 2 / (8.0 * g) - dp_over_rho_g)
    if kind == "line":
        return _S(-float(param) ** 2 / (8.0 * np.pi ** 2 * r_ ** 2 * g) - dp_over_rho_g)
    if kind == "rankine":
        if a is None:
            raise ValueError("kind='rankine' needs the core radius a")
        return _S(_F(rankine_pressure(r_, 0.0, param, a, 1.0, g)) / g - dp_over_rho_g)
    raise ValueError("kind must be 'solid', 'line' or 'rankine'")


def rotating_tank_free_surface(R: float, depth: float, Omega_tank: float, g: float = G, H: float | None = None,
                               closed: bool = False) -> dict:
    """Volume-conserving free surface of water spun up in a cylindrical tank: z = z_v + Ω_t²r²/2g.

    Book: §5.1, Eq. (5.6) with ω = 2Ω_tank (Fig. 5.2): the free surface is an isobar, a paraboloid with vertex z_v
    (possibly below the bottom). The water volume πR²·depth is conserved: π∫₀^{R²} clip(z_v + q/c, 0, H_lid) dq with
    q = r², c = 2g/Ω² (the clip at 0 uncovers the bottom; the clip at a lid applies only for ``closed=True``, the
    closed tank of Exercise 5.1). z_v by ``brentq`` (bracket: surface wholly below the bottom ↔ one metre above depth).

    Parameters: R tank radius [m]; depth initial water depth [m]; Omega_tank rotation rate [rad/s]; g [m/s²];
    H tank height [m] (optional); closed — a lid at H clips the surface (needs H).

    Returns
    -------
    dict: ``z_vertex`` [m], ``z_rim`` (unclipped surface height at the wall) [m], ``rim_rise`` z_rim − depth [m],
    ``centre_drop`` depth − z_vertex [m], ``r_dry`` radius of the uncovered bottom [m], ``uncovered_area`` [m²],
    ``spills`` (H given and z_rim > H: an open tank would overflow), ``r_top`` (radius where the surface meets the lid,
    NaN otherwise), ``volume_residual`` [m³].

    Validation: V1 (0.1, 0.1, 5.0) → z_v 0.093629, z_rim 0.106371 (rise = drop = Ω²R²/4g); volume residual < 1e-12;
    V6 Exercise 5.1 uncovered area with ``closed=True`` (private JSON). Label: analytic, book-value.
    """
    c = 2.0 * g / Omega_tank ** 2
    R2 = R ** 2
    Hl = float(H) if (closed and H is not None) else np.inf

    def water(zv):  # ∫₀^{R²} clip(z_v + q/c, 0, Hl) dq  (= water volume / π)
        q1 = min(max(-c * zv, 0.0), R2)
        q2 = min(max(c * (Hl - zv), q1), R2) if np.isfinite(Hl) else R2
        part = zv * (q2 - q1) + (q2 ** 2 - q1 ** 2) / (2.0 * c)
        return part + (Hl * (R2 - q2) if np.isfinite(Hl) else 0.0)

    target = R2 * depth
    zv = brentq(lambda v: water(v) - target, -R2 / c - 1.0, depth + 1.0, xtol=1e-14, rtol=1e-15)
    z_rim = zv + R2 / c
    r_dry = float(np.sqrt(min(max(-zv * c, 0.0), R2)))
    rt2 = (Hl - zv) * c if np.isfinite(Hl) else np.inf
    return dict(z_vertex=float(zv), z_rim=float(z_rim), rim_rise=float(z_rim - depth), centre_drop=float(depth - zv),
                r_dry=r_dry, uncovered_area=float(np.pi * r_dry ** 2),
                spills=bool(H is not None and z_rim > float(H)),
                r_top=float(np.sqrt(rt2)) if 0.0 < rt2 < R2 else float("nan"),
                volume_residual=float(np.pi * (water(zv) - target)))


def bernoulli_across_vortex(kind: str, r, z=0.0, param: float = 1.0, rho: float = 1000.0, g: float = G,
                            r_ref: float = 0.0, a: float | None = None):
    """Bernoulli function difference B(r) − B(r_ref), B = ½u_θ² + gz + p/ρ, across the circular streamlines.

    Book: §5.1 after (5.6) and (5.7) (compare (4.19), (4.69)–(4.72)): in the tank −½u_θ² + gz + p/ρ = const, so
    B(r) − B(r_ref) = ω²(r² − r_ref²)/4 (rotational flow: B varies across streamlines); for the line vortex
    ½u_θ² + gz + p/ρ = const, B is the same everywhere (r_ref = 0 is read as r_ref → ∞); for a Rankine vortex (core a)
    B changes only inside the core. Uses ``core.bernoulli.bernoulli_function`` with (5.6), (5.7) and the Rankine p.
    Parameters: kind "solid" (param = ω), "line" (param = Γ), "rankine" (param = Γ, core ``a``); r [m]; z [m];
    rho; g; r_ref [m]. Returns [m²/s²] (= J/kg). (solid, r = 0.1, ω = 10 → 0.25.) Label: analytic.
    """
    from .core.bernoulli import bernoulli_function

    def B(rr):
        rr_ = _F(rr)
        if kind == "solid":
            u = _F(solid_body_from_vorticity(rr_, param)[0])
            p = _F(solid_body_pressure(rr_, z, param, rho, g))
        elif kind == "line":
            with np.errstate(divide="ignore", invalid="ignore"):
                u = float(param) / (2.0 * np.pi * rr_)
            p = _F(line_vortex_pressure(rr_, z, param, rho, g))
        elif kind == "rankine":
            if a is None:
                raise ValueError("kind='rankine' needs the core radius a")
            u = _F(VX.rankine_vortex(rr_, param, a)[0])
            p = _F(rankine_pressure(rr_, z, param, a, rho, g))
        else:
            raise ValueError("kind must be 'solid', 'line' or 'rankine'")
        return _F(bernoulli_function(u, p, z, rho, g))

    if kind == "line" and r_ref == 0.0:
        return _S(B(r))  # r_ref = 0 read as r_ref → ∞: B_∞ = p_∞/ρ = 0 (and B is uniform anyway)
    return _S(B(r) - B(r_ref))


# ----------------------------------------------------------------------------------------------------------------------
# viscous stress without net force (§5.1, Exercise 5.4) and the rotating cylinder
# ----------------------------------------------------------------------------------------------------------------------
def polar_viscous_stress(u_r_expr, u_theta_expr, r: sp.Symbol, theta: sp.Symbol, mu=None) -> dict:
    """Viscous stresses σ_rr, σ_rθ, σ_θθ of a plane polar velocity field (sympy), σ = 2μS (incompressible (4.59)).

    Book: §5.1 before (5.7): σ_rθ = μ[(1/r)∂u_r/∂θ + r∂(u_θ/r)/∂r]; for (5.2) σ_rθ = −μΓ/πr² ≠ 0 (Exercise 5.4 also
    asks for σ_rr and σ_θθ, both 0). Uses ``core.curvilinear.strain_rate(system="cylindrical")`` (App. B).
    Parameters: u_r_expr, u_theta_expr sympy expressions in ``r``, ``theta``; mu (default Symbol "mu").
    Returns dict(sigma_rr, sigma_rtheta, sigma_thetatheta) (simplified, in r, θ). Label: symbolic.
    """
    mu = sp.Symbol("mu", positive=True) if mu is None else mu
    Rs, Ph, _ = CU.coordinates("cylindrical")
    sub = {r: Rs, theta: Ph}
    back = {Rs: r, Ph: theta}
    S = CU.strain_rate([sp.sympify(u_r_expr).subs(sub), sp.sympify(u_theta_expr).subs(sub), 0], "cylindrical")
    return dict(sigma_rr=sp.simplify(2 * mu * S[0, 0]).subs(back),
                sigma_rtheta=sp.simplify(2 * mu * S[0, 1]).subs(back),
                sigma_thetatheta=sp.simplify(2 * mu * S[1, 1]).subs(back))


def line_vortex_viscous_stress(r, Gamma, mu: float = 1e-3):
    """Viscous shear stress of the ideal line vortex, σ_rθ = −μΓ/(πr²) [Pa] (§5.1 before (5.7)).
    Nonzero everywhere (elements deform, Fig. 3.16) although the net viscous force is zero (Exercise 5.4, our D03).
    Parameters: r [m]; Gamma [m²/s]; mu [Pa s]. (r = 0.1, Γ = 1, μ = 1e-3 → −0.0318310 Pa.) Book: §5.1.
    Label: analytic."""
    return _S(-float(mu) * float(Gamma) / (np.pi * _F(r) ** 2))


def polar_net_viscous_force(u_theta: Callable, r, mu: float = 1e-3, h: float | None = None) -> dict:
    """Stress and net viscous force (θ-component) of an axisymmetric swirl u_θ(r), by central differences: σ_rθ =
    μr d(u_θ/r)/dr; (1/r²)d(r²σ_rθ)/dr (metric-correct); μ(u″ + u′/r − u/r²); and the naive dσ_rθ/dr (**wrong**:
    it misses the metric terms — kept for the discrimination test).
    Parameters: u_theta callable [m/s]; r [m]; mu [Pa s]; h [m] (default 1e-4 r).
    Returns dict(sigma_rtheta [Pa], force_metric, force_laplacian, force_naive [N/m³]). Book: §5.1, Exercise 5.4
    (our D03). Label: analytic."""
    r_ = _F(r)
    hh = 1e-4 * r_ if h is None else float(h)
    f = lambda q: _F(u_theta(q))  # noqa: E731
    sig = lambda q: float(mu) * q * (f(q + hh) / (q + hh) - f(q - hh) / (q - hh)) / (2.0 * hh)  # noqa: E731
    s0 = sig(r_)
    force_metric = ((r_ + hh) ** 2 * sig(r_ + hh) - (r_ - hh) ** 2 * sig(r_ - hh)) / (2.0 * hh) / r_ ** 2
    up = (f(r_ + hh) - f(r_ - hh)) / (2.0 * hh)
    upp = (f(r_ + hh) - 2.0 * f(r_) + f(r_ - hh)) / hh ** 2
    force_lap = float(mu) * (upp + up / r_ - f(r_) / r_ ** 2)
    force_naive = (sig(r_ + hh) - sig(r_ - hh)) / (2.0 * hh)
    return dict(sigma_rtheta=_S(s0), force_metric=_S(force_metric), force_laplacian=_S(force_lap),
                force_naive=_S(force_naive))


def vortex_stress_force(kind: str, r, mu: float = 1e-3, **p) -> dict:
    """σ_rθ and the net viscous force of a plane vortex by the three routes of our D03 (Exercise 5.4, N06, N09).

    Book: §5.1 ("irrotationality does not imply the absence of viscous stresses; it implies the absence of net viscous
    forces"). Routes for the θ-component of ∂σ_ij/∂x_j [N/m³]: ``force_divergence`` (1/r²)∂(r²σ_rθ)/∂r,
    ``force_laplacian`` μ(u_θ″ + u_θ′/r − u_θ/r²) = μ(∇²u)_θ, ``force_curl`` −μ(∇ × ω)_θ = μ∂ω_z/∂r ((4.40)); σ_rθ =
    μr∂(u_θ/r)/∂r [Pa]. kinds: "solid" (``omega`` = 2 s⁻¹: σ = 0, F = 0), "line" (``Gamma`` = 1: σ = −μΓ/πr², F = 0),
    "gaussian" (``Gamma`` = 1, ``sigma_c`` = 0.1: σ ≠ 0, F ≠ 0). Central differences with h = 1e-4 r.
    Returns dict(sigma_rtheta, force_divergence, force_laplacian, force_curl). Label: analytic.
    """
    r_ = float(r)
    if kind == "solid":
        om = float(p.get("omega", 2.0))
        ut = lambda q: 0.5 * om * _F(q)  # noqa: E731
        wz = lambda q: om + 0.0 * _F(q)  # noqa: E731
    elif kind == "line":
        Gm = float(p.get("Gamma", 1.0))
        ut = lambda q: Gm / (2.0 * np.pi * _F(q))  # noqa: E731
        wz = lambda q: 0.0 * _F(q)  # noqa: E731
    elif kind == "gaussian":
        Gm, sc = float(p.get("Gamma", 1.0)), float(p.get("sigma_c", 0.1))
        ut = lambda q: _F(VX.gaussian_vortex(q, Gm, sc)[0])  # noqa: E731
        wz = lambda q: _F(VX.gaussian_vortex(q, Gm, sc)[1])  # noqa: E731
    else:
        raise ValueError("kind must be 'solid', 'line' or 'gaussian'")
    d = polar_net_viscous_force(ut, r_, mu)
    hh = 1e-4 * r_
    fcurl = float(mu) * float((wz(r_ + hh) - wz(r_ - hh)) / (2.0 * hh))  # −μ(∇×ω)_θ = μ ∂ω_z/∂r
    return dict(sigma_rtheta=float(d["sigma_rtheta"]), force_divergence=float(d["force_metric"]),
                force_laplacian=float(d["force_laplacian"]), force_curl=fcurl)


def stress_vs_net_force_table(r: float = 0.05, Gamma: float = 0.01, mu: float = 1e-3, rho: float = 1000.0,
                              t: float = 100.0) -> list[dict]:
    """The three-row table of §5.1's principle at radius r with the same Γ: solid body (σ = 0, F = 0), line vortex
    (σ ≠ 0, F = 0), Lamb–Oseen with σ_c² = 4νt (σ ≠ 0, F ≠ 0). Rows: name, u_theta, omega_z, sigma_rtheta [Pa],
    net_force [N/m³]. Book: §5.1 (N09). Label: analytic."""
    nu = mu / rho
    s2 = 4.0 * nu * t
    om = Gamma / (np.pi * r ** 2)
    rows = []
    for name, fn, wz in (("solid body", lambda q: 0.5 * om * q, om),
                         ("line vortex", lambda q: Gamma / (2.0 * np.pi * q), 0.0),
                         ("Lamb–Oseen", lambda q: Gamma / (2.0 * np.pi * q) * (-np.expm1(-q ** 2 / s2)),
                          Gamma / (np.pi * s2) * np.exp(-r ** 2 / s2))):
        d = polar_net_viscous_force(fn, r, mu)
        rows.append(dict(name=name, u_theta=float(fn(r)), omega_z=float(wz), sigma_rtheta=float(d["sigma_rtheta"]),
                         net_force=float(d["force_metric"])))
    return rows


def rotating_cylinder_flow(r, a, omega):
    """Flow outside a solid cylinder of radius a turning at ω/2 in viscous fluid, u_θ = ωa²/2r (r ≥ a); with the
    cylinder included it is a Rankine vortex with Γ = πa²ω (§5.1; derived in §8.2, (8.11)).

    Parameters: r [m]; a [m]; omega ω [1/s] (the cylinder's "vorticity", rotation rate ω/2). Returns (u_theta [m/s],
    omega_z [1/s]) via ``core.vortices.rankine_vortex(r, πa²ω, a)`` (inside: the rigid cylinder).
    Validation: V1 no slip u_θ(a) = ωa/2; Γ = πa²ω; ω_z = ω inside, 0 outside. Label: analytic.
    Book: §5.1.
    """
    return VX.rankine_vortex(r, np.pi * float(a) ** 2 * float(omega), a)  # u_θ = ωa²/2r outside, Γ = πa²ω


def torque_per_length(r, Gamma, mu: float = 1e-3):
    """Torque per unit length transmitted across the circle r by the irrotational vortex's shear stress,
    T′ = 2πr·r·σ_rθ = −2μΓ [N m/m] — the same at every radius: the cylinder's torque is carried to infinity.
    Book: §5.1 (moment of momentum, §4.9; our a-D08). Parameters: r [m]; Gamma [m²/s]; mu [Pa s]. Label: analytic."""
    return _S(2.0 * np.pi * _F(r) ** 2 * _F(line_vortex_viscous_stress(r, Gamma, mu)))  # 2πr² σ_rθ = −2μΓ


def dissipation_outside_cylinder(a: float, R_out: float, Gamma: float, mu: float = 1e-3, rho: float = 1000.0) -> dict:
    """Energy budget of the irrotational vortex outside a rotating cylinder, per unit length: the viscous dissipation in
    a ≤ r ≤ R_out equals the power put in by the cylinder minus the power passed on at R_out.

    Book: §5.1 ("the resulting viscous dissipation of kinetic energy is exactly compensated by the work done at the
    surface of the cylinder"). Dissipation from ``core.constitutive.dissipation_rate`` of the Cartesian velocity gradient
    of u_θ = Γ/2πr at (r, 0) (pure strain G = [[0, −c], [−c, 0]], c = Γ/2πr²), integrated by ``quad``:
    ∫ρε 2πr dr = (μΓ²/π)(1/a² − 1/R_out²); power = |torque| × angular speed u_θ/r = μΓ²/πr² at r = a and r = R_out.
    Returns dict(dissipated, power_in, power_out, residual) [W/m]. ((0.1, 1.0, 1.0, 1e-3) → 0.0315127, 0.0318310,
    3.18310e-4.) Label: analytic, conserved.
    """
    def eps_rho(r):
        c = Gamma / (2.0 * np.pi * r ** 2)
        Gm = np.array([[0.0, -c], [-c, 0.0]])
        return rho * float(dissipation_rate(Gm, rho, mu)) * 2.0 * np.pi * r  # ρε 2πr

    diss, _ = quad(eps_rho, a, R_out, epsabs=0.0, epsrel=1e-12, limit=200)
    p_in = mu * Gamma ** 2 / (np.pi * a ** 2)  # |T′| u_θ(a)/a = 2μΓ·Γ/2πa²
    p_out = mu * Gamma ** 2 / (np.pi * R_out ** 2)
    return dict(dissipated=diss, power_in=p_in, power_out=p_out, residual=diss - (p_in - p_out))


# ======================================================================================================================
# §5.2 Kelvin's theorem: viscous decay, lock exchange, the four restrictions
# ======================================================================================================================
def lamb_oseen_circulation(r, t, Gamma: float, nu: float, t0: float = 0.0):
    """Circulation of the circle r about a Lamb–Oseen vortex and its rate: Γ(r, t) = Γ(1 − e^{−r²/4ν(t+t₀)}),
    ∂Γ/∂t = −Γ e^{−q} q/(t + t₀), q = r²/4ν(t + t₀).

    Book: §5.2, Eq. (5.11): DΓ/Dt = ∮(1/ρ)(∂σ_ij/∂x_j)dx_i — the circle is a material loop (u_r = 0) and the net viscous
    force on it drains circulation outward. Worked number (N13): Γ₀ = 0.01 m²/s, ν = 1e-6 m²/s, r = 5 mm, t = 10 s
    → Γ = 0.00464739 m²/s, ∂Γ/∂t = −3.34538e-4 m²/s². Parameters: r [m]; t [s]; Gamma [m²/s]; nu [m²/s]; t0 [s].
    Returns (Gamma_r, dGamma_dt).
    Validation: V2 sympy ∂Γ/∂t = the (5.11) line integral; V1 equals :func:`lamb_oseen_viscous_loop_integral`.
    Label: analytic, symbolic.
    """
    tau = _F(t) + float(t0)
    q = _F(r) ** 2 / (4.0 * float(nu) * tau)
    Gr = float(Gamma) * (-np.expm1(-q))
    dG = -float(Gamma) * np.exp(-q) * q / tau  # ∂Γ/∂t = Γ e^{−q} ∂q/∂t, ∂q/∂t = −q/τ
    return _S(Gr), _S(dG)


def lamb_oseen_viscous_loop_integral(r, t, Gamma: float, nu: float, t0: float = 0.0, n: int = 256) -> float:
    """The viscous term of (5.11), ∮ν∇²u·dx, on the circle r about a Lamb–Oseen vortex at time t [m²/s²] — computed
    from the velocity field (stencil Laplacian, spectral loop quadrature via ``core.vorticity.kelvin_force_terms``),
    independently of the closed form ∂Γ/∂t of :func:`lamb_oseen_circulation`, which it must equal.
    Book: §5.2, Eq. (5.11). Label: analytic, converged."""
    u = VX.lamb_oseen_field(Gamma, nu, t0)
    pts = VD.circle_loop_points((0.0, 0.0), float(r), n)
    visc = lambda x, tt=0.0: float(nu) * st.laplacian(u, _F(x), tt, 1e-3 * float(r), (2,))  # noqa: E731  ν∇²u
    return float(VD.kelvin_force_terms(pts, None, 1.0, visc_force=visc, t=float(t)).viscous)


def lock_exchange_initial_vorticity_rate(rho1: float, rho2: float, delta: float, g: float = G) -> float:
    """Initial baroclinic spin-up at the interface of a lock exchange: Dω_z/Dt = 2(ρ₂ − ρ₁)g/((ρ₂ + ρ₁)δ) [1/s²].

    Book: §5.2, Fig. 5.5 (heavy ρ₂ on the **left**, light ρ₁ on the right) and Exercise 5.5 (answer printed there): our
    D07 evaluates the baroclinic term (1/ρ²)∇ρ × ∇p of (5.28) at the interface with ∂ρ/∂x = (ρ₁ − ρ₂)/δ, ∂p/∂y = −ρ̄g,
    ρ̄ = (ρ₁ + ρ₂)/2. Positive = counterclockwise (x right, y up): the heavy fluid slumps under the light.
    Worked number (N15): 1000 / 1025 kg/m³, δ = 0.1 m → 2.422222 s⁻². Parameters: rho1, rho2 [kg/m³]; delta [m]; g.
    Validation: V1 equals ``core.vorticity.baroclinic_term`` on :func:`lock_exchange_fields` at (0, H/2); ∇p × ∇ρ
    flips the sense (wrong variant). Label: analytic.
    """
    return float(2.0 * (rho2 - rho1) * g / ((rho2 + rho1) * delta))  # D07


def lock_exchange_fields(rho1: float = 1000.0, rho2: float = 1025.0, delta: float = 0.1, g: float = G,
                         H: float = 1.0, p_top: float = 0.0, L: float | None = None, nx: int = 128,
                         ny: int = 64) -> dict:
    """Density and pressure of the lock exchange at the instant t = 0⁺ the barrier is removed (Fig. 5.5), y up,
    interface at x = 0: ρ = ρ̄ − (Δρ/2)tanh(2x/δ) (Δρ = ρ₂ − ρ₁, ρ₂ on the left) and the hydrostatic pressure of the
    mean density, p = p_top + ρ̄g(H − y) (the fluid has not moved yet). The baroclinic term is then
    (∇ρ × ∇p)_z/ρ² = −ρ′(x)ρ̄g/ρ², peaking at 2(ρ₂ − ρ₁)g/((ρ₂ + ρ₁)δ) on the interface.

    Returns dict with the contract callables of (x, y): ``rho``, ``p``, ``grad_rho``, ``grad_p`` (→ (2, …));
    plus ``rho_fn``, ``p_fn`` in the field convention f(x, t); the grid ``x``, ``y``, ``X``, ``Y`` (x ∈ [−L/2, L/2],
    L = 2H) with ``rho_grid``, ``p_grid``, ``baroclinic_z`` [1/s²]; ``rate`` (the interface value) and
    ``circulation_rate`` = ∫ baroclinic dA over the tank [m²/s²]. Book: §5.2, Fig. 5.5; Exercise 5.5. Label: analytic.
    """
    rb = 0.5 * (rho1 + rho2)
    dr = rho2 - rho1
    rho = lambda x, y: rb - 0.5 * dr * np.tanh(2.0 * _F(x) / delta) + 0.0 * _F(y)  # noqa: E731
    p = lambda x, y: p_top + rb * g * (H - _F(y)) + 0.0 * _F(x)  # noqa: E731
    grad_rho = lambda x, y: np.stack([-dr / delta / np.cosh(2.0 * _F(x) / delta) ** 2 + 0.0 * _F(y),  # noqa: E731
                                      0.0 * _F(x) + 0.0 * _F(y)])
    grad_p = lambda x, y: np.stack([0.0 * _F(x) + 0.0 * _F(y), -rb * g + 0.0 * _F(x) + 0.0 * _F(y)])  # noqa: E731
    L = 2.0 * H if L is None else float(L)
    xg = np.linspace(-L / 2, L / 2, nx)
    yg = np.linspace(0.0, H, ny)
    X, Y = np.meshgrid(xg, yg, indexing="xy")
    gr, gp, rr = grad_rho(X, Y), grad_p(X, Y), rho(X, Y)
    bz = (gr[0] * gp[1] - gr[1] * gp[0]) / rr ** 2  # (∇ρ × ∇p)_z/ρ²   (5.28)
    circ = float(np.trapezoid(np.trapezoid(bz, xg, axis=1), yg))
    return dict(rho=rho, p=p, grad_rho=grad_rho, grad_p=grad_p,
                rho_fn=lambda x, t=0.0: rho(_F(x)[0], _F(x)[1]), p_fn=lambda x, t=0.0: p(_F(x)[0], _F(x)[1]),
                x=xg, y=yg, X=X, Y=Y, rho_grid=rr, p_grid=p(X, Y), baroclinic_z=bz,
                rate=lock_exchange_initial_vorticity_rate(rho1, rho2, delta, g), circulation_rate=circ)


def baroclinic_rate_2d(grad_rho, grad_p, rho: float) -> float:
    """(∇ρ × ∇p)_z/ρ² for plane gradients [1/s²] — the baroclinic term (5.28) from two 2-vectors (E4's readout).
    ([10, 0], [0, −9810], 1000 → −0.0981.) Book: §5.6, Eq. (5.28). Label: analytic."""
    gr, gp = _F(grad_rho), _F(grad_p)
    return float((gr[0] * gp[1] - gr[1] * gp[0]) / float(rho) ** 2)


_KELVIN_TERMS = {
    "viscous": r"$\oint_C\big(\frac1\rho\frac{\partial\sigma_{ij}}{\partial x_j}\big)dx_i$ — the net viscous force "
               "along C, (5.11)",
    "baroclinic": r"$-\oint_C\frac{dp}{\rho}$ — baroclinic: $\nabla\rho\times\nabla p\neq0$, (5.10)",
    "body": r"$\oint_C\mathbf g\cdot d\mathbf x$ — a nonconservative body force, (5.10)",
    "coriolis": r"$-\oint_C(2\boldsymbol\Omega\times\mathbf u)\cdot d\mathbf x$ — frame (Coriolis) term: only the "
                r"absolute circulation $\Gamma_a$ of (5.33) is conserved",
}


def kelvin_hypotheses(inviscid: bool = True, barotropic: bool = True, conservative: bool = True,
                      inertial: bool = True) -> dict:
    """Kelvin's theorem as a decision table: which of the four restrictions hold and which term of (5.10) survives.

    Book: §5.2, Eqs. (5.8)–(5.11) and the four restrictions: (1) no net viscous force along C, (2) conservative body
    forces, (3) barotropic ρ = ρ(p), (4) inertial frame. DΓ/Dt = 0 only when all four hold; each broken one leaves its
    term (the three ways to create or destroy vorticity, plus the frame).
    Returns dict(holds: bool, surviving_terms: list ⊂ ["viscous", "baroclinic", "body", "coriolis"], verdict: str,
    text: str (each surviving term with its TeX)). Validation: V1 all 16 combinations. Label: analytic.
    """
    flags = [("viscous", inviscid), ("baroclinic", barotropic), ("body", conservative), ("coriolis", inertial)]
    surv = [k for k, ok in flags if not ok]
    holds = not surv
    verdict = kelvin_hypotheses_text(inviscid, barotropic, conservative, inertial)
    text = (r"$\frac{D\Gamma}{Dt}=0$ (5.8)" if holds
            else "surviving term(s): " + "; ".join(_KELVIN_TERMS[k] for k in surv))
    return dict(holds=holds, surviving_terms=surv, verdict=verdict, text=text)


def kelvin_hypotheses_text(inviscid: bool = True, barotropic: bool = True, conservative: bool = True,
                           inertial: bool = True) -> str:
    """The verdict sentence of :func:`kelvin_hypotheses` (plain text, for exact-text parity rows).
    Book: §5.2. Label: analytic."""
    names = [("viscous", inviscid), ("baroclinic", barotropic), ("nonconservative body force", conservative),
             ("rotating frame", inertial)]
    broken = [n for n, ok in names if not ok]
    if not broken:
        return "Kelvin holds: DΓ/Dt = 0 for every material loop."
    return "Kelvin fails: DΓ/Dt ≠ 0 because of " + ", ".join(broken) + "."


def kelvin_scenario_gamma(name: str, t, **p):
    """Alias of ``core.vorticity.kelvin_scenario_circulation`` (Γ(t) of an E3 scenario) [m²/s]. Book: §5.2.
    Label: analytic."""
    return VD.kelvin_scenario_circulation(name, float(t), **p)


# ======================================================================================================================
# §5.4 diffusion of a vortex sheet, stretching of a tube, Burgers' balance
# ======================================================================================================================
def diffusing_vortex_sheet(y, t, gamma: float, nu: float):
    """Viscous diffusion of a plane vortex sheet ω(y, 0) = γδ(y)e_z (Exercise 5.6): (5.13) reduces to ∂ω_z/∂t = ν∂²ω_z/∂y²,
    ω_z(y, t) = γ/(2√(πνt)) exp(−y²/4νt), u(y, t) = −(γ/2) erf(y/2√(νt)).

    Book: §5.4, Exercise 5.6 (the ω_z form is printed there); u is ours (ω_z = −∂u/∂y, u odd). ∫ω dy = γ for all t and
    u(±∞) = ∓γ/2 — the sheet's velocity jump γ spreads over a thickness ~√(νt) (→ Ch. 8 Stokes' first problem).
    Parameters: y [m]; t [s] (> 0); gamma γ [m/s] (counterclockwise +); nu [m²/s]. Returns (u [m/s], omega_z [1/s]).
    ((0, 1, 1, 1e-6)[1] = 282.095; (1e-3, 1, 1, 1e-6)[0] = −0.260250.)
    Validation: V2 sympy diffusion equation; V1 ∫ω dy = γ; V3 FTCS agreement (``core.diffusion``). Label: symbolic,
    analytic, converged.
    """
    y_, t_ = _F(y), _F(t)
    s = np.sqrt(float(nu) * t_)
    w = float(gamma) / (2.0 * np.sqrt(np.pi) * s) * np.exp(-y_ ** 2 / (4.0 * s ** 2))  # Exercise 5.6
    u = -0.5 * float(gamma) * erf(y_ / (2.0 * s))
    return _S(u), _S(w)


def stretched_tube(L, L0: float, omega0: float, A0: float) -> dict:
    """An inviscid vortex tube stretched from length L₀ to L: volume AL and strength ωA are conserved (Helmholtz 2, 4 and
    continuity), so A = A₀L₀/L and ω = ω₀L/L₀ — twice as long, twice the vorticity (our D18, first part).
    Parameters: L, L0 [m]; omega0 [1/s]; A0 [m²]. Returns dict(omega [1/s], A [m²], Gamma [m²/s], radius_ratio).
    ((2, 1, 10, 1e-4) → ω 20, A 5e-5, Γ 1e-3.) Book: §5.6, Eq. (5.32). Label: analytic."""
    L_ = _F(L)
    A = A0 * L0 / L_
    return dict(omega=_S(omega0 * L_ / L0), A=_S(A), Gamma=_S(omega0 * A0 + 0.0 * L_), radius_ratio=_S(np.sqrt(L0 / L_)))


def burgers_balance(R, Gamma: float = 1e-3, alpha: float = 1.0, nu: float = 1e-6) -> dict:
    """The three terms of (5.13) for ω_z in Burgers' vortex at radius R (closed forms, our D18): with
    ω_z = (αΓ/4πν)e^{−aR²}, a = α/4ν: advective u_R∂ω_z/∂R = αaR²ω_z, stretching αω_z, diffusion
    ν∇²ω_z = αω_z(aR² − 1); steady: advective = stretching + diffusion.
    Worked number (N21): α = 1 s⁻¹, ν = 1e-6 m²/s, Γ = 1e-3 m²/s → core √(4ν/α) = 2 mm, peak ω_z = 79.5775 s⁻¹.
    Returns dict(omega_z, advective, stretching, diffusion, residual [1/s²], core_radius [m]). Book: Exercise 5.12
    (§5.4). Label: analytic."""
    R_ = _F(R)
    a = alpha / (4.0 * nu)
    w = alpha * Gamma / (4.0 * np.pi * nu) * np.exp(-a * R_ ** 2)
    adv = alpha * a * R_ ** 2 * w
    stretch = alpha * w
    diff = alpha * w * (a * R_ ** 2 - 1.0)
    return dict(omega_z=_S(w), advective=_S(adv), stretching=_S(stretch), diffusion=_S(diff),
                residual=_S(adv - stretch - diff), core_radius=float(np.sqrt(4.0 * nu / alpha)))


# ======================================================================================================================
# §5.6 rotating frame: residual (5.20), absolute vorticity, stretching/tilting, the column, the ring of air
# ======================================================================================================================
def rotating_ns_residual(u, p, rho, x, t: float = 0.0, nu: float = 1e-6, Omega=(0.0, 0.0, 0.5),
                         g_eff=(0.0, 0.0, -9.81), h: float = 1e-4, ht: float | None = None):
    """Residual of the momentum equation in a steadily rotating frame, Eq. (5.20) (Boussinesq, ∇·u = 0 (5.19)):
    ∂u/∂t + (u·∇)u + 2Ω × u + (1/ρ)∇p − g − ν∇²u [m/s²] (0 for a solution).

    Book: §5.6, Eq. (5.20): ∂u_i/∂t + u_ju_{i,j} + 2ε_ijkΩ_ju_k = −(1/ρ)p_{,i} + g_i + νu_{i,jj}, g the effective gravity
    (with the centrifugal part; ch04 D18). The Coriolis term is the kinematic +2Ω × u (``core.rotating``'s
    ``coriolis_acceleration`` returns the apparent force −2Ω × u).
    Parameters: u (3-D) [m/s]; p [Pa]; rho [kg/m³] (constant or callable); x (3,) or (3, N); t; nu; Omega (3,);
    g_eff (3,) or callable [m/s²]; h, ht. Returns (3,) or (3, N).
    Validation: V1 fluid co-rotating (u = 0) with p hydrostatic + centrifugal and g_eff: 0; inertial oscillation: 0.
    Label: analytic.
    """
    x_ = _F(x)
    a = _F(acceleration(u, x_, t, h, ht).a)  # ∂u/∂t + (u·∇)u
    U = st.ev(u, x_, t, (3,))
    Om = _F(Omega).reshape((3,) + (1,) * (x_.ndim - 1)) * np.ones_like(U)
    cor = 2.0 * np.cross(Om, U, axis=0)  # 2 ε_ijk Ω_j u_k
    rho_ = st.ev(rho, x_, t) if callable(rho) else float(rho)
    gp = st.grad(p, x_, t, h) / rho_
    gv = (st.ev(g_eff, x_, t, (3,)) if callable(g_eff)
          else _F(g_eff).reshape((3,) + (1,) * (x_.ndim - 1)) * np.ones_like(U))
    lap = float(nu) * st.laplacian(lambda X, T: st.ev(u, X, T, (3,)), x_, t, h, (3,))
    return a + cor + gp - gv - lap  # Eq. (5.20)


def absolute_vorticity(omega_rel, Omega):
    """Absolute vorticity ω + 2Ω (planetary vorticity 2Ω), §5.6 after (5.30) [1/s]; inverse of ch03's
    ``vorticity_in_rotating_frame`` (ω′ = ω − 2Ω). Scalars are z-components. Book: §5.6. Label: analytic."""
    return _S(_F(omega_rel) + 2.0 * _F(Omega))  # ω_a = ω + 2Ω


STRAIN_PRESETS = ("axial_stretch", "axial_compression", "shear_tilt", "planar", "burgers")


def strain_preset(name: str, rate: float = 1.0) -> np.ndarray:
    """Steady linear background flows u = Gx for stretching/tilting (each a steady Euler flow): "axial_stretch" (and
    "burgers") diag(−r/2, −r/2, r); "axial_compression" diag(r/2, r/2, −r); "shear_tilt" w = r x (G[2, 0] = r);
    "planar" diag(r, −r, 0) (plane strain: ω ∥ e_z is neither stretched nor tilted). Returns G (3, 3) [1/s].
    Book: §5.6, (5.32). Label: analytic."""
    r = float(rate)
    if name in ("axial_stretch", "burgers"):
        return np.diag([-0.5 * r, -0.5 * r, r])
    if name == "axial_compression":
        return np.diag([0.5 * r, 0.5 * r, -r])
    if name == "shear_tilt":
        G_ = np.zeros((3, 3))
        G_[2, 0] = r
        return G_
    if name == "planar":
        return np.diag([r, -r, 0.0])
    raise ValueError(f"unknown preset {name!r}; choose from {STRAIN_PRESETS}")


def uniform_strain_vorticity(omega0, G, t, rate: float = 1.0):
    """Weak vorticity carried by a steady linear flow u = Gx with ν = 0: Dω/Dt = (ω·∇)u = Gω ⇒ ω(t) = e^{Gt}ω₀.

    Book: §5.4, Eq. (5.13) with ν = 0 and §5.6, Eq. (5.32): stretching multiplies ω along the line, tilting turns it.
    The background flows of :func:`strain_preset` have zero (or parallel) vorticity and are steady Euler solutions; the
    carried vorticity is weak enough not to change them (Cauchy's vorticity formula: ω follows the material line
    element, e^{Gt} is the deformation gradient). Axial strain → ω₀e^{αt}; shear tilt w = sx on ω₀ = ω_x0 e_x →
    ω_z = s t ω_x0; planar strain on ω ∥ e_z → no change (2-D).
    Parameters: omega0 (3,) [1/s]; G (3, 3) [1/s] or a preset name (with ``rate``); t [s] (scalar or (T,)).
    Returns ω(t) (3,) or (T, 3) [1/s]. (([0, 0, 1], "axial_stretch", 2)[2] = 7.389056; ([1, 0, 0], "shear_tilt", 2)[2]
    = 2.) Validation: V2 d/dt(e^{Gt}ω₀) = G e^{Gt}ω₀; V1 the presets. Label: analytic.
    """
    G_ = strain_preset(G, rate) if isinstance(G, str) else _F(G)
    w0 = _F(omega0)
    if np.ndim(t) == 0:
        return expm(G_ * float(t)) @ w0  # ω(t) = e^{Gt} ω₀
    return np.stack([expm(G_ * float(tt)) @ w0 for tt in _F(t)])


def stretching_tilting_scenario(preset: str = "axial_stretch", t: float = 1.0, rate: float = 1.0,
                                angle: float = 0.0, omega0: float = 1.0) -> dict:
    """E5's state: a vortex-line element of initial magnitude ω₀ at ``angle`` [rad] from e_z (toward e_x) carried by
    the preset linear flow; returns ω(t), |ω|/ω₀, the line's angle from e_z, and the stretching/tilting split of (5.31)
    at time t (``core.vorticity.stretching_tilting_split``), plus a status string.
    Book: §5.6, Eqs. (5.31)–(5.32). Label: analytic."""
    G_ = strain_preset(preset, rate)
    w0 = omega0 * np.array([np.sin(angle), 0.0, np.cos(angle)])
    w = uniform_strain_vorticity(w0, G_, t)
    sp_ = VD.stretching_tilting_split(w, G_)
    mag = float(np.linalg.norm(w))
    tilt = float(np.arccos(np.clip(w[2] / mag, -1.0, 1.0))) if mag > 0 else 0.0
    if abs(sp_["rate"]) < 1e-14 and sp_["tilt_rate"] < 1e-14:
        status = "2-D: neither stretching nor tilting"
    elif sp_["tilt_rate"] > abs(sp_["rate"]):
        status = "tilting: ω turns"
    else:
        status = "stretching: |ω| grows" if sp_["rate"] > 0 else "compression: |ω| falls"
    return dict(omega=w, growth=mag / omega0, angle=tilt, stretching=sp_["stretching"], tilting=sp_["tilting"],
                stretch_rate=sp_["rate"], tilt_rate=sp_["tilt_rate"], G=G_, status=status)


def helix(s, a: float = 1.0, c: float = 0.3):
    """Point of an arc-length-parametrised helix of radius a and pitch 2πc: x = (a cos τ, a sin τ, cτ),
    τ = s/√(a² + c²) — a curved vortex line for the natural coordinates of Fig. 5.9. Returns (3,) (or (3, N) for an
    array s) [m]. Book: §5.6, Fig. 5.9. Label: analytic."""
    tau = _F(s) / np.sqrt(a ** 2 + c ** 2)
    return np.stack([a * np.cos(tau), a * np.sin(tau), c * tau])


def helix_frame(s, a: float = 1.0, c: float = 0.3) -> dict:
    """Natural frame (e_s, e_n, e_m) of the helix of :func:`helix` at arc length s, Fig. 5.9 (closed forms):
    e_s the unit tangent; **e_n points away from the centre of curvature** (the book's convention — minus the Frenet
    normal, i.e. radially outward for a helix); e_m = e_s × e_n; curvature κ = a/(a² + c²), torsion τ = c/(a² + c²).
    (s = 0, a = 1, c = 0.3 → κ 0.917431, τ 0.275229.) Returns dict(e_s, e_n, e_m (3,), curvature [1/m], torsion
    [1/m]). Book: §5.6, Fig. 5.9, Eq. (5.31). Label: analytic."""
    k = np.sqrt(a ** 2 + c ** 2)
    th = float(s) / k
    es = np.array([-a * np.sin(th), a * np.cos(th), c]) / k
    en = np.array([np.cos(th), np.sin(th), 0.0])  # away from the axis = away from the centre of curvature
    em = np.cross(es, en)
    return dict(e_s=es, e_n=en, e_m=em, curvature=a / k ** 2, torsion=c / k ** 2)


def helical_vortex_line(a: float, c: float, s):
    """Alias of :func:`helix` with the argument order (a, c, s) [m] → (3, N). Book: §5.6, Fig. 5.9. Label: analytic."""
    return helix(s, a, c)


def frenet_frame(curve, s) -> dict:
    """Natural (s, n, m) frame along any sampled curve (Fig. 5.9), by explicit 2nd-order differences in s (one-sided at
    the ends): e_s the unit tangent, e_n away from the centre of curvature (book convention), e_m = e_s × e_n, and the
    curvature κ [1/m]. Parameters: curve (3, N) [m]; s (N,) [m] (uniform spacing). Returns dict(e_s, e_n, e_m (3, N),
    kappa (N,)). Validation: V1 helix: κ = a/(a² + c²), e_n radially outward. Book: §5.6, Fig. 5.9. Label: analytic.
    """
    X = _F(curve)
    ds = float(_F(s)[1] - _F(s)[0])

    def d1(F):
        D = np.empty_like(F)
        D[:, 1:-1] = (F[:, 2:] - F[:, :-2]) / (2 * ds)
        D[:, 0] = (-3 * F[:, 0] + 4 * F[:, 1] - F[:, 2]) / (2 * ds)
        D[:, -1] = (3 * F[:, -1] - 4 * F[:, -2] + F[:, -3]) / (2 * ds)
        return D
    r1 = d1(X)
    r2 = d1(r1)
    T = r1 / np.linalg.norm(r1, axis=0)
    cr = np.cross(r1.T, r2.T).T
    kappa = np.linalg.norm(cr, axis=0) / np.linalg.norm(r1, axis=0) ** 3
    B = cr / np.linalg.norm(cr, axis=0)
    N = np.cross(B.T, T.T).T  # Frenet normal (toward the centre of curvature)
    en = -N
    return dict(e_s=T, e_n=en, e_m=np.cross(T.T, en.T).T, kappa=kappa)


def column_relative_vorticity(h, h0: float, zeta0: float = 0.0, f: float = 1.0e-4):
    """Relative vorticity of a fluid column whose height changes from h₀ to h (inviscid, barotropic, rotating frame):
    (ζ + f)/h = (ζ₀ + f)/h₀ ⇒ ζ = (ζ₀ + f)h/h₀ − f, f = 2Ω_z the local planetary vorticity [1/s].

    Book: §5.6, Fig. 5.10 and (5.33): our D20 applies Γ_a = (ζ + f)A = const (5.33) to a thin loop round the column
    with mass ρAh fixed. Stretching (h > h₀) spins it up cyclonically, squashing anticyclonically — the seed of
    potential vorticity (Ch. 13). ((1100, 1000, 0, 1e-4) → 1.0e-5 s⁻¹.) Parameters: h, h0 [m]; zeta0 [1/s]; f [1/s].
    Returns ζ [1/s]. Validation: V4 (ζ + f)/h invariant; small-change limit = 2Ω∂w/∂z integrated. Label: analytic,
    conserved.
    """
    return _S((float(zeta0) + float(f)) * _F(h) / float(h0) - float(f))  # (ζ + f)/h = const


_DEPTH_PRESETS = ("ridge", "trough", "slope")


def column_over_slope(x, depth="ridge", lat_deg: float = 45.0, zeta0: float = 0.0, h0: float = 1000.0,
                      bump: float = 0.2, width: float = 1.0e5, Omega: float = OMEGA_EARTH) -> dict:
    """ζ of a column carried across varying depth h(x) (E7 column mode, Fig. 5.10): (ζ + f)/h conserved,
    ζ(x) = (ζ₀ + f)h(x)/h₀ − f, f = 2Ω sin φ; h₀ is the undisturbed depth where the column has ζ₀.

    ``depth``: "ridge" h = h₀(1 − bump·e^{−x²/width²}) (column squashed → anticyclonic), "trough"
    h = h₀(1 + bump·e^{−x²/width²}) (stretched → cyclonic), "slope" h = h₀(1 + bump·tanh(x/width)), or a callable h(x).
    Parameters: x [m]; lat_deg [°] (interface in degrees); zeta0 [1/s]; h0 [m]; bump [–]; width [m]; Omega [rad/s].
    Returns dict(x, h [m], zeta [1/s], f [1/s], ratio (ζ + f)/h [1/(m s)]). Book: §5.6, Fig. 5.10. Label: analytic.
    """
    x_ = _F(x)
    if callable(depth):
        h = _F(depth(x_))
    elif depth == "ridge":
        h = h0 * (1.0 - bump * np.exp(-(x_ / width) ** 2))
    elif depth == "trough":
        h = h0 * (1.0 + bump * np.exp(-(x_ / width) ** 2))
    elif depth == "slope":
        h = h0 * (1.0 + bump * np.tanh(x_ / width))
    else:
        raise ValueError(f"depth must be one of {_DEPTH_PRESETS} or a callable")
    f = float(coriolis_parameter(np.deg2rad(lat_deg), Omega))
    z = _F(column_relative_vorticity(h, h0, zeta0, f))
    return dict(x=x_, h=h, zeta=z, f=f, ratio=(z + f) / h)


def relative_circulation_after_move(Gamma0: float, A0: float, lat0_deg: float, A1: float, lat1_deg: float,
                                    Omega: float = OMEGA_EARTH) -> float:
    """Relative circulation of a horizontal ring of air after it changes area and latitude (degrees at this interface):
    Γ₁ = Γ₀ + 2Ω(A₀ sin φ₀ − A₁ sin φ₁) [m²/s] — (5.33) with the planetary vorticity's normal part 2Ω sin φ. A ring at
    rest moved poleward acquires anticyclonic (negative in the NH) circulation. ((0, π(5e5)², 30, π(5e5)², 60) →
    −4.19261e7 m²/s.) Book: §5.6, Eq. (5.33) (Exercise 5.10). Label: analytic."""
    s0, s1 = np.sin(np.deg2rad(lat0_deg)), np.sin(np.deg2rad(lat1_deg))
    return float(Gamma0 + 2.0 * Omega * (A0 * s0 - A1 * s1))  # Γ_a conserved


def helical_swirl_field(a: float = 1.0) -> Callable:
    """The swirl u = (0, aRz, 0) in (u_R, u_φ, u_z) (Cartesian u = az(−y, x, 0)): ω = (−aR, 0, 2az), vortex lines
    zR² = const — the N02 test field (the flow of Exercise 5.3; the line equation is our (5.3) integration).
    Returns u(x, t) (3-D) [m/s]; a [1/(m s)]. Book: §5.1, Eq. (5.3). Label: analytic."""
    def u(x, t=0.0):
        X = _F(x)
        return np.stack([-a * X[2] * X[1], a * X[2] * X[0], 0.0 * X[0]])
    return u


def cellular_flow(x, t: float = 0.0, U: float = 1.0, ell: float = 1.0):
    """Steady inviscid cellular flow ψ = Uℓ sin(x/ℓ) sin(y/ℓ): u = U sin(x/ℓ)cos(y/ℓ), v = −U cos(x/ℓ)sin(y/ℓ),
    ω = 2ψ/ℓ² (a function of ψ ⇒ a steady Euler solution). x (2,) or (2, N) [m]; returns (2,) or (2, N) [m/s].
    Book: §5.2 test field (E3 "cellular"). Label: analytic."""
    X = _F(x)
    return U * np.stack([np.sin(X[0] / ell) * np.cos(X[1] / ell), -np.cos(X[0] / ell) * np.sin(X[1] / ell)])


def abc_flow(A: float = 1.0, B: float = 1.0, C: float = 1.0) -> Callable:
    """The ABC flow u = (A sin z + C cos y, B sin x + A cos z, C sin y + B cos x) (x [m], u [m/s]): a steady Euler
    (Beltrami) flow with ω = ∇ × u = u — the inviscid field of the Helmholtz demos (C05, C06, N30). Returns u(x, t).
    Alias of ``core.vortices.abc_flow_field``. Book: §5.3–5.4, §5.6 (test field). Label: analytic."""
    return VX.abc_flow_field(A, B, C)


# ======================================================================================================================
# §5.7 interactions: pairs, wall, rings
# ======================================================================================================================
def vortex_pair(Gamma1: float, Gamma2: float, h: float) -> dict:
    """Two ideal line vortices a distance h apart (Figs. 5.11, 5.12): the speeds each induces at the other,
    V₁ = Γ₁/2πh (at vortex 2), V₂ = Γ₂/2πh (at vortex 1); for Γ₁ + Γ₂ ≠ 0 the pair rotates about the centre of vorticity
    G at h₁ = Γ₂h/(Γ₁ + Γ₂) from vortex 1 with angular rate (Γ₁ + Γ₂)/2πh² (Exercise 5.18, our D21); for Γ₂ = −Γ₁ it
    translates at |Γ₁|/2πh.

    Returns dict(V1, V2 [m/s], centre_from_1 [m] (inf if Γ₁ + Γ₂ = 0), rotation_rate [rad/s], period [s],
    translation_speed [m/s]). ((1, 1, 1) → 0.159155, 0.5, 0.318310, 19.7392.) Book: §5.7.
    Validation: V1 period vs ``point_vortex_evolve``; V5 García & Haziot (2023). Label: analytic, benchmark.
    """
    V1 = Gamma1 / (2.0 * np.pi * h)
    V2 = Gamma2 / (2.0 * np.pi * h)
    tot = Gamma1 + Gamma2
    if abs(tot) > 1e-14 * max(abs(Gamma1), abs(Gamma2), 1e-300):
        rate = tot / (2.0 * np.pi * h ** 2)
        return dict(V1=V1, V2=V2, centre_from_1=Gamma2 * h / tot, rotation_rate=rate,
                    period=2.0 * np.pi / abs(rate), translation_speed=0.0)
    return dict(V1=V1, V2=V2, centre_from_1=float("inf"), rotation_rate=0.0, period=float("inf"),
                translation_speed=abs(Gamma1) / (2.0 * np.pi * h))


def vortex_near_wall_speed(Gamma: float, h: float) -> float:
    """Speed of a line vortex at distance h from a plane wall, V_A = Γ/(4πh) [m/s], parallel to the wall (induced by
    its image −Γ at distance 2h; Fig. 5.14; our D22). Book: §5.7. Label: analytic."""
    return float(Gamma / (4.0 * np.pi * h))


def ring_dynamics(rings0, t_eval, wall_z: float | None = None, core: str = "uniform", rtol: float = 1e-9,
                  atol: float = 1e-12) -> dict:
    """Coaxial thin vortex rings moving under their own and each other's induced velocity, optionally with a plane wall
    z = wall_z (image rings −Γ mirrored in the wall): leap-frogging rings (Exercise 5.15) and a ring approaching a wall
    (Fig. 5.15: it widens and slows).

    Book: §5.7, Fig. 5.15 and the leap-frogging paragraph (qualitative in the book). Model (ours): each ring i has radius
    R_i, axial position z_i, circulation Γ_i (vorticity along +e_φ, self-propelled toward +z for Γ > 0) and core radius
    a_i with a_i²R_i conserved (core volume); dR_i/dt, dz_i/dt = Kelvin self-speed (``ring_self_velocity``) +
    Σ_{j≠i} ``ring_ring_velocity`` + the image rings' contributions (each ring's own image included).

    Parameters: rings0 list of dict(R, z, Gamma, a) [m, m, m²/s, m]; t_eval (T,) [s]; wall_z [m] or None; core
    ("uniform", "hollow", "gaussian"); rtol, atol.
    Returns dict(t, R (T, n), z (T, n), a (T, n), impulse (T,) = ΣΓπR² [m⁴/s]).
    Validation: V4 two rings without a wall: ΣΓπR² conserved; leap-frog passes (z₁ − z₂ changes sign) repeatedly; ring
    + wall: R increases monotonically and the approach speed decreases (Fig. 5.15b). Label: conserved, qualitative.
    """
    R0 = np.array([float(r["R"]) for r in rings0])
    Z0 = np.array([float(r["z"]) for r in rings0])
    Gm = np.array([float(r["Gamma"]) for r in rings0])
    a0 = np.array([float(r["a"]) for r in rings0])
    n = R0.size
    vol = a0 ** 2 * R0
    te = np.atleast_1d(_F(t_eval))

    def rhs(t, y):
        R, Z = y[:n], y[n:]
        a = np.sqrt(vol / R)
        dR = np.zeros(n)
        dZ = _F(BS.ring_self_velocity(R, a, 1.0, core)) * Gm
        src = [(R[j], Z[j], Gm[j], j) for j in range(n)]
        if wall_z is not None:
            src += [(R[j], 2.0 * wall_z - Z[j], -Gm[j], -1) for j in range(n)]
        for i in range(n):
            for (Rj, Zj, Gj, j) in src:
                if j == i:
                    continue
                uR, uz = BS.ring_ring_velocity(R[i], Z[i], Rj, Zj, Gj)
                dR[i] += uR
                dZ[i] += uz
        return np.concatenate([dR, dZ])

    sol = solve_ivp(rhs, (te[0], te[-1]), np.concatenate([R0, Z0]), method="DOP853", rtol=rtol, atol=atol, t_eval=te)
    if not sol.success:
        raise RuntimeError(sol.message)
    R, Z = sol.y[:n].T, sol.y[n:].T
    return dict(t=te, R=R, z=Z, a=np.sqrt(vol / R), impulse=np.sum(Gm * np.pi * R ** 2, axis=1))


# ======================================================================================================================
# §5.8 vortex sheet
# ======================================================================================================================
def vortex_sheet_strength(u_above, u_below, convention: str = "ccw"):
    """Strength of a vortex sheet = jump in tangential velocity (circulation per unit length) [m/s].

    Book: §5.8: dΓ = u₂ds + v dn − u₁ds − v dn = (u₂ − u₁)ds, Γ ≡ dΓ/ds = u₂ − u₁ (u₁ above, u₂ below, the circuit
    counterclockwise) — ``convention="ccw"`` (default; negative for the clockwise filaments of Fig. 5.16). The figure
    caption prints dΓ/ds = u₁ − u₂, the clockwise magnitude — ``convention="caption"`` (or "cw"). ⚠️ The book reuses Γ
    for this per-length quantity; here it is ``gamma`` [m/s], never the circulation ``Gamma`` [m²/s].
    ((−1, 1) → 2.0; caption → −2.0.) Validation: V1 a counterclockwise point-vortex row gives γ > 0 with u_above = −γ/2
    (wrong-convention check). Label: analytic.
    """
    # DEVIATION: Fig. 5.16's caption prints dΓ/ds = u₁ − u₂; we follow the text's u₂ − u₁ (counterclockwise positive,
    # consistent with plane vorticity) by default and keep the caption's sign as convention="caption".
    if convention == "ccw":
        return _S(_F(u_below) - _F(u_above))  # γ = u₂ − u₁  (text)
    if convention in ("caption", "cw"):
        return _S(_F(u_above) - _F(u_below))  # caption: u₁ − u₂
    raise ValueError("convention must be 'ccw' (text) or 'caption'")


def _sheet_row(N: int, L: float):
    ds = L / N
    xs = -0.5 * L + ds * (np.arange(N) + 0.5)
    return np.stack([xs, np.zeros(N)]), ds


def discrete_sheet_u(x, y, gamma: float, N: int, L: float = 1.0):
    """Tangential velocity u at (x, y) of a finite row of N line vortices of strength γL/N at the cell midpoints of
    [−L/2, L/2] (Fig. 5.16's row of filaments) [m/s]; arrays of y are accepted.
    ((0, 0.05, 2, 100) → −0.936551.) Book: §5.8, Fig. 5.16. Label: analytic."""
    S, ds = _sheet_row(int(N), float(L))
    y_ = np.atleast_1d(_F(y))
    P = np.stack([np.full(y_.shape, float(x)), y_])
    u = BS.vortex_sheet_velocity(P, S, gamma, ds)[0]
    return float(u[0]) if np.ndim(y) == 0 else u


def discrete_sheet_convergence(gamma: float, N_list=(10, 100, 1000), L: float = 1.0, y_max: float = 0.1) -> dict:
    """How a row of N filaments approaches a continuous vortex sheet: the L1 error
    ∫_{−y_max}^{y_max}|u_N(0, y) − u_cont(0, y)|dy of the tangential velocity across the sheet's middle, with u_cont
    from ``core.biot_savart.continuous_sheet_velocity``; ∝ 1/N (the discrete profile differs from the jump over a
    band of the order of the spacing L/N).

    Book: §5.8, Fig. 5.16. Adaptive ``quad`` on each side of y = 0 (the continuous profile jumps there).
    Returns dict(N (K,), l1_error (K,) [m²/s], u_above, u_below (K,) = u_N(0, ±y_max) [m/s]).
    ((2.0) → 0.043968, 0.0044131, 0.00044128.) Label: analytic, converged.
    """
    Ns = np.atleast_1d(np.asarray(N_list, dtype=int))
    l1, ua, ub = [], [], []
    for n in Ns:
        S, ds = _sheet_row(int(n), float(L))

        def err(yy):
            u_n = float(BS.vortex_sheet_velocity(np.array([0.0, yy]), S, gamma, ds)[0])
            u_c = float(BS.continuous_sheet_velocity(0.0, yy, gamma, L)[0])
            return abs(u_n - u_c)
        pts = list(np.linspace(0.0, y_max, 9)[1:-1]) if n <= 50 else None
        e1, _ = quad(err, 0.0, y_max, limit=2000, epsabs=1e-13, epsrel=1e-10, points=pts)
        e2, _ = quad(err, -y_max, 0.0, limit=2000, epsabs=1e-13, epsrel=1e-10,
                     points=[-q for q in pts] if pts else None)
        l1.append(e1 + e2)
        ua.append(discrete_sheet_u(0.0, y_max, gamma, int(n), L))
        ub.append(discrete_sheet_u(0.0, -y_max, gamma, int(n), L))
    return dict(N=Ns, l1_error=np.array(l1), u_above=np.array(ua), u_below=np.array(ub))


def sheet_rollup(N: int = 100, gamma: float = 1.0, amplitude: float = 0.01, delta: float = 0.05, t_eval=None,
                 L: float = 1.0, rtol: float = 1e-8, atol: float = 1e-11) -> dict:
    """Roll-up of a periodic vortex sheet (period L) represented by N point vortices with Krasny's δ-smoothing — the
    Kelvin–Helmholtz roll-up of a shear layer into cat's-eye vortices (forward link to Ch. 11).

    Book: §5.8 (the sheet as a row of line vortices, Fig. 5.16) + §5.7 (vortices move with the flow). Our model: vortex
    j (strength Γ_j = γL/N, counterclockwise +) starts at x = s_j = (j + ½)L/N, y = A sin(2πs_j/L); velocities from the
    periodic row kernel u − iv = (Γ_j/2iL)cot(π(z − z_j)/L), in real form u = −(Γ/2L) sinh(2πΔy/L)/D,
    v = (Γ/2L) sin(2πΔx/L)/D, D = cosh(2πΔy/L) − cos(2πΔx/L) + δ² (δ dimensionless; δ = 0 is the singular row).
    Far from the sheet u → ∓γ/2.

    Parameters: N; gamma γ [m/s]; amplitude A [m]; delta δ [–]; t_eval (T,) [s] (default 0…4L/γ, 9 frames); L [m].
    Returns dict(t, x (T, N), y (T, N), Gamma_each, mean_y (T,) (conserved: ΣΓy), N).
    Validation: V4 mean y conserved; the roll-up shape is qualitative. Label: conserved, qualitative.
    """
    te = np.linspace(0.0, 4.0 * L / max(abs(gamma), 1e-300), 9) if t_eval is None else np.atleast_1d(_F(t_eval))
    s = (np.arange(N) + 0.5) * L / N
    x0, y0 = s.copy(), amplitude * np.sin(2.0 * np.pi * s / L)
    Gj = gamma * L / N
    k = 2.0 * np.pi / L

    def rhs(t, yv):
        X, Y = yv[:N], yv[N:]
        dx = X[:, None] - X[None, :]
        dy = Y[:, None] - Y[None, :]
        D = np.cosh(k * dy) - np.cos(k * dx) + delta ** 2
        with np.errstate(divide="ignore", invalid="ignore"):
            iD = np.where(D > 0.0, 1.0 / np.where(D > 0.0, D, 1.0), 0.0)
        u = -(Gj / (2.0 * L)) * np.sum(np.sinh(k * dy) * iD, axis=1)
        v = (Gj / (2.0 * L)) * np.sum(np.sin(k * dx) * iD, axis=1)
        return np.concatenate([u, v])

    sol = solve_ivp(rhs, (te[0], te[-1]), np.concatenate([x0, y0]), method="DOP853", rtol=rtol, atol=atol, t_eval=te)
    if not sol.success:
        raise RuntimeError(sol.message)
    X, Y = sol.y[:N].T, sol.y[N:].T
    return dict(t=te, x=X, y=Y, Gamma_each=Gj, mean_y=Y.mean(axis=1), N=N)


# ======================================================================================================================
# Explainer scenario functions (scalar-callable where the explainer mirrors them)
# ======================================================================================================================
VORTEX_PRESSURE_KINDS = ("solid", "line", "rankine", "cylinder")


def vortex_pressure_scenario(kind: str, r, z: float = 0.0, rho: float = 1000.0, mu: float = 1e-3, g: float = G,
                             **p) -> dict:
    """E2 (``vortex_pressure_funnel``): the four vortices of §5.1 in one dispatcher — velocity, pressure, Bernoulli
    function, viscous stress, net viscous force and the radial momentum balance (5.5a) at radius r.

    kinds (keywords, defaults): "solid" (``omega`` = 2 s⁻¹ — a tank at 1 rad/s), "line" (``Gamma`` = 1 m²/s),
    "rankine" (``Gamma`` = 1, core ``a`` = 0.1 m), "cylinder" (``omega`` = 2, ``a`` = 0.1: Γ = πa²ω; r < a lies inside
    the solid cylinder). Reference pressure ``p_ref`` = 0 (p_o on the axis for "solid", p_∞ otherwise).
    Tornado preset: kind "rankine", rho = 1.2, Gamma = 1e4, a = 50 → u_max 31.831 m/s, core-edge deficit 607.93 Pa.

    Returns
    -------
    dict: ``u_theta`` [m/s], ``omega_z`` [1/s], ``p`` [Pa], ``B`` (= B − B_ref) [J/kg], ``sigma_rtheta`` [Pa],
    ``net_viscous_force`` [N/m³] (θ-component, (1/r²)d(r²σ_rθ)/dr), ``centripetal`` u_θ²/r [m/s²], ``dp_dr`` [Pa/m]
    (central difference of p — (5.5a) says dp_dr/ρ = centripetal), ``radial_residual`` dp_dr/ρ − centripetal,
    ``surface_z`` (free-surface/isobar height through the reference point) [m], ``torque`` 2πr²σ_rθ [N m/m],
    ``in_fluid``, ``status``.

    Book: §5.1, Eqs. (5.1), (5.2), (5.5)–(5.7), Figs. 5.2–5.3. Validation: V1 the full balance of each kind (radial
    residual ≈ 0; σ and F per kind). Label: analytic.
    """
    pref = float(p.get("p_ref", 0.0))
    r_ = float(r)
    a = 0.0
    if kind == "solid":
        om = float(p.get("omega", 2.0))
        uf = lambda q: 0.5 * om * q  # noqa: E731
        pf = lambda q: float(solid_body_pressure(q, z, om, rho, g, pref))  # noqa: E731
        wz = om
        B = float(bernoulli_across_vortex("solid", r_, z, om, rho, g))
        zs = float(isobar_height(r_, om, 0.0, "solid", g))
        status = "rotational: B grows outward, no viscous stress (S = 0)"
        in_fluid, rigid = True, True
    elif kind in ("line", "rankine", "cylinder"):
        if kind == "line":
            Gam = float(p.get("Gamma", 1.0))
        elif kind == "rankine":
            Gam, a = float(p.get("Gamma", 1.0)), float(p.get("a", 0.1))
        else:
            a = float(p.get("a", 0.1))
            Gam = np.pi * a ** 2 * float(p.get("omega", 2.0))
        if a > 0.0:
            uf = lambda q: float(VX.rankine_vortex(q, Gam, a)[0])  # noqa: E731
            pf = lambda q: float(rankine_pressure(q, z, Gam, a, rho, g, pref))  # noqa: E731
            wz = float(VX.rankine_vortex(r_, Gam, a)[1])
            B = float(bernoulli_across_vortex("rankine", r_, z, Gam, rho, g, r_ref=np.inf, a=a))  # B − B_∞
            zs = float(isobar_height(r_, Gam, 0.0, "rankine", g, a=a))
        else:
            uf = lambda q: Gam / (2.0 * np.pi * q)  # noqa: E731
            pf = lambda q: float(line_vortex_pressure(q, z, Gam, rho, g, pref))  # noqa: E731
            wz = 0.0
            B = float(bernoulli_across_vortex("line", r_, z, Gam, rho, g))
            zs = float(isobar_height(r_, Gam, 0.0, "line", g))
        rigid = r_ < a
        in_fluid = not (kind == "cylinder" and rigid)
        if rigid:
            status = ("inside the core: rigid rotation, no viscous stress" if kind == "rankine"
                      else "inside the solid cylinder (not fluid)")
        else:
            status = "irrotational: viscous stress ≠ 0 but no net viscous force"
    else:
        raise ValueError(f"kind must be one of {VORTEX_PRESSURE_KINDS}")
    hh = 1e-5 * max(r_, 1e-9)
    ut = uf(r_)
    dpdr = (pf(r_ + hh) - pf(r_ - hh)) / (2.0 * hh)
    if rigid:
        sig = fnet = 0.0  # S = 0: no viscous stress at all
    else:
        visc = polar_net_viscous_force(lambda q: np.vectorize(uf)(q), r_, mu)
        sig, fnet = float(visc["sigma_rtheta"]), float(visc["force_metric"])
    cent = ut ** 2 / r_
    return dict(kind=kind, r=r_, u_theta=float(ut), omega_z=float(wz), p=float(pf(r_)), B=B, sigma_rtheta=sig,
                net_viscous_force=fnet, centripetal=cent, dp_dr=dpdr, radial_residual=dpdr / rho - cent,
                surface_z=zs, torque=2.0 * np.pi * r_ ** 2 * sig, in_fluid=bool(in_fluid), status=status)


def baroclinic_element_scenario(tilt: float = 0.5, grad_rho: float = 5.0, R: float = 0.05, rho0: float = 1000.0,
                                g: float = G, numeric: bool = False) -> dict:
    """E4 (``baroclinic_torque``): a fluid disc of radius R in hydrostatic pressure ∇p = (0, −ρ₀g) whose isopycnals are
    tilted by ``tilt`` [rad] from the isobars (∇ρ = |∇ρ|(sin tilt, −cos tilt): tilt 0 = stable, barotropic).

    Closed forms (our D06, linear fields): centre-of-mass offset x_G − x_c = R²∇ρ/4ρ₀, torque about G
    τ = πR⁴(∇ρ × ∇p)_z/4ρ₀, I_G ≈ πρ₀R⁴/2, spin-up 2τ/I_G = (∇ρ × ∇p)_z/ρ₀² = −|∇ρ|g sin(tilt)/ρ₀ — the baroclinic term
    (5.28). ``numeric=True`` adds the quadrature route ``core.vorticity.pressure_torque_on_element`` (ratio → 1 as
    R → 0). Returns dict(grad_rho (2,), grad_p (2,), offset (2,), torque [N m/m], I_G, spin_up, baroclinic [1/s²],
    sense, status[, numeric]). Book: §5.2, Fig. 5.6; §5.6, Eq. (5.28). Label: analytic.
    """
    gr = grad_rho * np.array([np.sin(tilt), -np.cos(tilt)])
    gp = np.array([0.0, -rho0 * g])
    cz = float(gr[0] * gp[1] - gr[1] * gp[0])  # (∇ρ × ∇p)_z
    torque = np.pi * R ** 4 * cz / (4.0 * rho0)
    IG = np.pi * rho0 * R ** 4 / 2.0
    out = dict(grad_rho=gr, grad_p=gp, offset=R ** 2 * gr / (4.0 * rho0), torque=torque, I_G=IG,
               spin_up=2.0 * torque / IG, baroclinic=cz / rho0 ** 2,
               sense="counterclockwise" if cz > 0 else ("clockwise" if cz < 0 else "none"),
               status=("barotropic: no torque" if abs(cz) <= 1e-12 * grad_rho * rho0 * g
                       else "baroclinic: the pressure force misses G"))
    if numeric:
        out["numeric"] = VD.pressure_torque_on_element("linear", "linear", (0.0, 0.0), R, grad_rho=gr, grad_p=gp,
                                                       rho0=rho0)
    return out
