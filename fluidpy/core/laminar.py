"""Exact laminar solutions of the Navier–Stokes equations: parallel flows, circular Couette flow, and the unsteady
diffusion family (Stokes' first and second problems, the thickening vortex sheet, the decaying line vortex).

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e (2012), Ch. 8, §8.2 (Eqs. (8.4)–(8.12)), §8.4 (Eqs. (8.20)–(8.31),
Examples 8.4–8.6, Exercises 8.26, 8.30), §8.5 (Eqs. (8.33)–(8.38)). Every equation was transcribed from the rendered
page images (chapters/pages/ch08/p340–p365). Reused by Ch. 9 (Rayleigh problem vs Blasius, δ₉₉), Ch. 11 (Couette,
Poiseuille and Taylor–Couette base states), Ch. 12 (laminar reference profiles, linear total stress), Ch. 13 (Ekman
spin-up, Stokes/Ekman layers), Ch. 16 (arteries, Womersley).

Conventions (analysis ch08 §9 R1–R3)
* Channel (§8.2): x along the plates, **y across, walls at y = 0 (fixed) and y = h (moving at U)**. The pressure
  gradient is passed as the book's ``dpdx`` = dp/dx [Pa/m] (favourable < 0 drives +x flow); ``G=`` is accepted as a
  keyword alias meaning **−dp/dx** (the ch04 ``plane_poiseuille`` / ``exact_solution`` convention). Giving both raises.
* Pipe: cylindrical (R, φ, z), z along the axis, radius a; ``dpdz`` = dp/dz [Pa/m] (``G=`` = −dp/dz).
* Circular Couette: inner radius R1 turning at Omega1, outer R2 at Omega2 [rad/s], counter-clockwise positive.
* §8.4: η = y/√(νt) (the book's similarity variable, **not** y/(2√(νt)), although Figs. 8.13–8.14 plot the latter).
* §8.5: ω = oscillation angular frequency [rad/s]; δ_e = √(2ν/ω) (literature penetration depth) vs the book's
  δ ~ 4√(ν/ω).
* Defaults: water μ = 1.0e-3 Pa s, ρ = 1000 kg/m³, ν = 1e-6 m²/s. All closed forms are vectorised and scalar-callable.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import quad
from scipy.special import erf, erfc, erfcinv, erfinv

from ._util import as_scalar_if_0d, require_positive

__all__ = ["channel_flow", "channel_flow_rate", "channel_shear_stress", "channel_backflow_threshold",
           "couette_poiseuille_state", "pipe_poiseuille", "pipe_shear_stress", "pipe_wall_stress", "pipe_flow_rate",
           "pipe_friction_factor", "circular_couette", "circular_couette_pressure", "circular_couette_shear_stress",
           "circular_couette_power", "circular_couette_state", "similarity_variable", "stokes_first_problem",
           "stokes_first_vorticity", "diffusion_thickness", "transition_width", "stokes_first_stopped",
           "stokes_first_state", "vortex_sheet_diffusion", "temporal_bl_wall_stress", "line_vortex_decay",
           "line_vortex_spinup", "stokes_second_problem", "stokes_layer", "stokes_layer_envelope",
           "stokes_layer_state"]

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d


def _grad(value, G, name: str = "dpdx") -> float | np.ndarray:
    """Resolve the book's pressure gradient from ``dpdx``/``dpdz`` or the ch04 alias ``G`` = −dp/dx."""
    if G is not None:
        if value is not None and np.any(_F(value) != 0.0):
            raise ValueError(f"give either {name} (= dp/d·) or G (= −dp/d·), not both")
        return -_F(G)
    if value is None:
        raise ValueError(f"{name} (dp/d·, Pa/m) or G (= −dp/d·) is required")
    return _F(value)


# ======================================================================================================================
# §8.2 plane Couette–Poiseuille flow
# ======================================================================================================================
def channel_flow(y, h, U: float = 0.0, dpdx: float | None = 0.0, mu: float = 1e-3, G: float | None = None):
    """Steady, fully developed flow between parallel plates: the Couette–Poiseuille profile.

    Book: §8.2, Eq. (8.5), from the reduced momentum equations (8.4a,b) integrated twice with u(0) = 0, u(h) = U.

    Parameters
    ----------
    y : float or array_like — distance from the fixed wall [m], 0 ≤ y ≤ h (values outside are evaluated, not masked).
    h : float — plate spacing [m], > 0.
    U : float — speed of the upper plate y = h [m/s] (lower plate at rest).
    dpdx : float — constant pressure gradient dp/dx [Pa/m] (favourable < 0).
    mu : float — dynamic viscosity μ [Pa s].
    G : float, optional — keyword alias for −dp/dx [Pa/m] (ch04 convention); excludes a nonzero ``dpdx``.

    Returns
    -------
    u : float or ndarray — streamwise velocity [m/s].

    Assumptions: incompressible, constant μ, steady, fully developed (∂u/∂x = 0 ⇒ v = 0), p = p(x) only.
    Validation (planned): V1 NS residual ≈ 0 (``core.navier_stokes.ns_incompressible_terms``), BCs; V2 sympy
    μu″ = dp/dx; parity with ch04 ``exact_solution("couette", G=−dpdx)`` and ``plane_poiseuille``.
    Label: analytic, symbolic.
    """
    y_, h = _F(y), float(h)
    require_positive("h", h)
    g = _grad(dpdx, G)
    u = float(U) / h * y_ - g / (2.0 * float(mu)) * y_ * (h - y_)  # Eq. (8.5)
    return _S(u)


def channel_flow_rate(h, U: float = 0.0, dpdx: float | None = 0.0, mu: float = 1e-3, G: float | None = None):
    """Volume flow rate per unit width and mean velocity of the Couette–Poiseuille flow.

    Book: §8.2, the displays after (8.5): Q = ∫₀ʰ u dy = (Uh/2)[1 − (h²/6μU) dp/dx] and V ≡ Q/h. The printed middle
    form of V lacks the 1/h in front of the integral (analysis §9 R9; dimensions m²/s vs m/s) — V = Q/h is coded.
    Coded in the expanded form Q = Uh/2 − h³(dp/dx)/(12μ), which stays finite at U = 0 (plane Poiseuille).

    Parameters
    ----------
    h : plate spacing [m];  U : upper-plate speed [m/s];  dpdx : dp/dx [Pa/m];  mu : μ [Pa s];
    G : optional alias −dp/dx [Pa/m].

    Returns
    -------
    (Q, V) : flow rate per unit width [m²/s] and mean velocity [m/s].

    Validation (planned): V1 ``quad`` of :func:`channel_flow` equals Q; V = Q/h; wrong variant: printed V (units).
    Label: analytic.
    """
    h = float(h)
    g = _grad(dpdx, G)
    Q = float(U) * h / 2.0 - h ** 3 * g / (12.0 * float(mu))  # Q = (Uh/2)[1 − (h²/6μU) dp/dx], expanded
    # DEVIATION: V = Q/h — the printed middle form of V omits the 1/h (analysis §9 R9)
    return _S(Q), _S(Q / h)


def channel_shear_stress(y, h, U: float = 0.0, dpdx: float | None = 0.0, mu: float = 1e-3, G: float | None = None):
    """Signed shear stress τ = μ du/dy of the Couette–Poiseuille profile (8.5).

    Book: §8.2, τ = μU/h (plane Couette) and τ = −(h/2 − y) dp/dx (plane Poiseuille, the display after (8.5));
    the general form is their sum, τ(y) = μU/h − (h/2 − y) dp/dx — linear in y (Fig. 8.4d draws |τ|).

    Parameters
    ----------
    y : [m];  h : [m];  U : [m/s];  dpdx : [Pa/m];  mu : [Pa s];  G : optional −dp/dx [Pa/m].

    Returns
    -------
    tau : float or ndarray — τ_xy [Pa], signed (positive: the fluid above drags the fluid below toward +x).

    Validation (planned): V1 μ times the analytic derivative of (8.5); |τ_w| = (h/2)|dp/dx| at both walls for U = 0.
    Label: analytic.
    """
    y_, h = _F(y), float(h)
    g = _grad(dpdx, G)
    tau = float(mu) * float(U) / h - (h / 2.0 - y_) * g  # τ = μ du/dy of Eq. (8.5)
    return _S(tau)


def channel_backflow_threshold(U, h, mu: float = 1e-3):
    """Adverse pressure gradient above which the Couette–Poiseuille flow reverses near the fixed wall.

    Book: §8.2, Fig. 8.4b (drawn, not given). Ours: du/dy(0) = U/h − (h/2μ) dp/dx < 0 ⇔ dp/dx > 2μU/h².

    Parameters
    ----------
    U : upper-plate speed [m/s];  h : spacing [m];  mu : μ [Pa s].

    Returns
    -------
    dpdx_star : float — threshold dp/dx [Pa/m] (for U > 0 backflow appears iff dp/dx > dpdx_star).

    Validation (planned): V1 bisection of min u on a y grid. Label: analytic.
    """
    return _S(2.0 * float(mu) * _F(U) / float(h) ** 2)


def couette_poiseuille_state(h, U, dpdx: float | None = 0.0, mu: float = 1e-3, G: float | None = None) -> dict:
    """Every number an explainer shows for one Couette–Poiseuille flow (curation §9, E1/IF1).

    Book: §8.2, Eq. (8.5), the flow-rate displays and Fig. 8.4.

    Parameters
    ----------
    h : [m];  U : upper-plate speed [m/s];  dpdx : dp/dx [Pa/m];  mu : μ [Pa s];  G : optional −dp/dx [Pa/m].

    Returns
    -------
    dict of floats: Q [m²/s], V [m/s], Q_couette = Uh/2, Q_poiseuille = −h³(dp/dx)/12μ [m²/s] (Q = sum);
    tau_bottom = τ(0), tau_top = τ(h) [Pa] (signed τ_xy); backflow (bool: somewhere inside 0 < y < h the fluid moves
    against the plate's direction (against +x when U = 0) — for U > 0 iff dp/dx > 2μU/h², Fig. 8.4b; mirrored for
    U < 0); threshold = 2μU/h² [Pa/m]; ratio = (dp/dx)/threshold (NaN if U = 0); zero_flow_dpdx = 6μU/h² [Pa/m] (Q = 0);
    u_max, y_umax (largest u and where, walls included) [m/s, m]; u_min, y_umin; y_reversal = h − 2μU/(h dp/dx) when
    inside (0, h), else NaN [m]; dpdx [Pa/m] (the resolved gradient).

    Validation (planned): V1 against :func:`channel_flow` sampled on a fine grid. Label: analytic.
    """
    h, U, mu = float(h), float(U), float(mu)
    g = float(_grad(dpdx, G))
    Q, V = channel_flow_rate(h, U, g, mu)
    cands = [0.0, h]
    y_rev = np.nan
    if g != 0.0:
        ys = h / 2.0 + mu * U / (h * g)  # du/dy = 0
        if 0.0 < ys < h:
            cands.append(ys)
        yr = h - 2.0 * mu * U / (h * g)  # u = 0 away from y = 0
        if 0.0 < yr < h:
            y_rev = yr
    vals = [float(channel_flow(c, h, U, g, mu)) for c in cands]
    i_max, i_min = int(np.argmax(vals)), int(np.argmin(vals))
    scale = max(abs(U), abs(g) * h ** 2 / (8.0 * mu), 1e-300)
    yy = np.linspace(0.0, h, 403)[1:-1]  # interior points
    sgn = 1.0 if U >= 0.0 else -1.0  # flow direction set by the plate (+x when U = 0)
    backflow = bool(np.min(sgn * _F(channel_flow(yy, h, U, g, mu))) < -1e-12 * scale)  # flow against the plate
    thr = float(channel_backflow_threshold(U, h, mu))
    return dict(Q=float(Q), V=float(V), Q_couette=U * h / 2.0, Q_poiseuille=-h ** 3 * g / (12.0 * mu),
                tau_bottom=float(channel_shear_stress(0.0, h, U, g, mu)),
                tau_top=float(channel_shear_stress(h, h, U, g, mu)), backflow=backflow,
                threshold=thr, ratio=(g / thr) if thr != 0.0 else np.nan, zero_flow_dpdx=6.0 * mu * U / h ** 2,
                u_max=vals[i_max], y_umax=cands[i_max], u_min=vals[i_min], y_umin=cands[i_min],
                y_reversal=float(y_rev), dpdx=g)


# ======================================================================================================================
# §8.2 round pipe (circular Poiseuille flow)
# ======================================================================================================================
def pipe_poiseuille(R, a, dpdz: float | None = None, mu: float = 1e-3, G: float | None = None):
    """Hagen–Poiseuille profile in a round tube.

    Book: §8.2, Eq. (8.6), from 0 = −dp/dz + (μ/R) d/dR(R du_z/dR) integrated twice, A ln R dropped (bounded at R = 0),
    u_z(a) = 0.

    Parameters
    ----------
    R : float or array_like — cylindrical radius [m], 0 ≤ R ≤ a.
    a : float — tube radius [m], > 0.
    dpdz : float — dp/dz [Pa/m] (< 0 drives +z flow).
    mu : float — μ [Pa s].
    G : float, optional — alias −dp/dz [Pa/m] (ch04 ``exact_solution("pipe_poiseuille")`` convention).

    Returns
    -------
    u_z : float or ndarray — axial velocity [m/s].

    Assumptions: steady, fully developed, axisymmetric, no swirl, constant μ.
    Validation (planned): V1 parity with ch04 ``exact_solution("pipe_poiseuille", G=−dpdz)`` and ch03 ``pipe_profile``
    far downstream; u_z(a) = 0; V2 sympy residual. Label: analytic, symbolic.
    """
    R_, a = _F(R), float(a)
    require_positive("a", a)
    g = _grad(dpdz, G, "dpdz")
    return _S((R_ ** 2 - a ** 2) / (4.0 * float(mu)) * g)  # Eq. (8.6)


def pipe_shear_stress(R, dpdz: float | None = None, G: float | None = None):
    """Shear stress τ = τ_zR = μ ∂u_z/∂R in Poiseuille pipe flow — linear in R.

    Book: §8.2, Eq. (8.7) (with τ_zR = μ(∂u_R/∂z + ∂u_z/∂R), u_R = 0, Appendix B).

    Parameters
    ----------
    R : cylindrical radius [m];  dpdz : dp/dz [Pa/m];  G : optional −dp/dz.

    Returns
    -------
    tau : [Pa], signed (negative for dp/dz < 0: the outer fluid holds the inner fluid back).

    Validation (planned): V1 μ × derivative of (8.6). Label: analytic.
    """
    g = _grad(dpdz, G, "dpdz")
    return _S(_F(R) / 2.0 * g)  # Eq. (8.7)


def pipe_wall_stress(a, dpdz: float | None = None, G: float | None = None):
    """Wall shear stress of pipe flow, the largest |τ|: τ₀ = (a/2) dp/dz.

    Book: §8.2, Eq. (8.8); also valid for suitable averages in turbulent pipe flow (it is the CV force balance
    πa²Δp = 2πaLτ₀ — ours).

    Parameters
    ----------
    a : tube radius [m];  dpdz : dp/dz [Pa/m];  G : optional −dp/dz.

    Returns
    -------
    tau0 : [Pa], signed like (8.8).

    Validation (planned): V1 (8.7) at R = a; V4 CV momentum balance on a pipe slug. Label: analytic, conserved.
    """
    g = _grad(dpdz, G, "dpdz")
    return _S(_F(a) / 2.0 * g)  # Eq. (8.8)


def pipe_flow_rate(a, dpdz: float | None = None, mu: float = 1e-3, G: float | None = None):
    """Hagen–Poiseuille flow rate, mean velocity and centre-line velocity.

    Book: §8.2, displays after (8.8): Q = ∫₀ᵃ u 2πR dR = −(πa⁴/8μ) dp/dz, V = Q/(πa²) = −(a²/8μ) dp/dz; u_max = 2V
    (ours).

    Parameters
    ----------
    a : radius [m];  dpdz : dp/dz [Pa/m];  mu : μ [Pa s];  G : optional −dp/dz.

    Returns
    -------
    (Q, V, u_max) : [m³/s], [m/s], [m/s].

    Validation (planned): V1 ``quad`` of 2πR u_z; u_max = 2V. Label: analytic.
    """
    a, mu = float(a), float(mu)
    g = _grad(dpdz, G, "dpdz")
    Q = -np.pi * a ** 4 / (8.0 * mu) * g  # Q = −(πa⁴/8μ) dp/dz
    V = -a ** 2 / (8.0 * mu) * g  # V = Q/(πa²)
    return _S(Q), _S(V), _S(2.0 * V)


def pipe_friction_factor(Re):
    """Darcy friction factor of laminar pipe flow, f = 64/Re with Re = Vd/ν (ours; not printed in the book).

    Book: §8.2 (follows from (8.8) and V: Δp/L = f (ρV²/2)/d, τ₀ = f ρV²/8). Parameters: Re [–] (diameter, mean
    velocity). Returns f [–]. Validation (planned): V1 from :func:`pipe_wall_stress` and :func:`pipe_flow_rate`.
    Label: analytic.
    """
    return _S(64.0 / _F(Re))


# ======================================================================================================================
# §8.2 circular Couette flow
# ======================================================================================================================
def _cc_coeffs(R1, R2, Omega1, Omega2) -> tuple[float, float]:
    R1, R2, O1, O2 = float(R1), float(R2), float(Omega1), float(Omega2)
    if R1 < 0 or not R2 > R1:
        raise ValueError("need 0 <= R1 < R2 (R2 may be np.inf)")
    if np.isinf(R2):  # exact limit R2 → ∞ of the constants after (8.9)
        return O2, (O1 - O2) * R1 ** 2
    d = R2 ** 2 - R1 ** 2
    A = (O2 * R2 ** 2 - O1 * R1 ** 2) / d  # constants after Eq. (8.9)
    B = -(O2 - O1) * R1 ** 2 * R2 ** 2 / d
    return A, B


def circular_couette(R, R1, R2, Omega1, Omega2, return_coeffs: bool = False):
    """Steady flow between concentric rotating cylinders (circular Couette flow).

    Book: §8.2, Eqs. (8.9)–(8.10); limits (8.11) (R2 → ∞, Ω₂ = 0: the irrotational vortex Γ = 2πΩ₁R₁², (5.2)) and
    (8.12) (R1, Ω₁ → 0: solid-body rotation, (5.1)). ``R2=np.inf`` and ``R1=0`` are evaluated as the exact limits
    (general Ω₂: u = Ω₂R + (Ω₁ − Ω₂)R₁²/R for R2 = inf).

    Parameters
    ----------
    R : float or array_like — radius [m] (R1 ≤ R ≤ R2).
    R1, R2 : float — inner and outer cylinder radii [m] (0 ≤ R1 < R2 ≤ inf).
    Omega1, Omega2 : float — angular velocities [rad/s] (counter-clockwise +).
    return_coeffs : bool — also return A [1/s] and B [m²/s] of (8.9).

    Returns
    -------
    u_phi : [m/s]  (or (u_phi, A, B)).

    Validation (planned): V1 u(R1) = Ω₁R₁, u(R2) = Ω₂R₂; limits vs ch05 ``rotating_cylinder_flow(r, a, omega=2Ω₁)``
    (ω there is the cylinder's vorticity 2Ω₁) and ``core.vortices.solid_body_rotation``; ω_z = 2A; V2 sympy ODE.
    Label: analytic, symbolic.
    """
    A, B = _cc_coeffs(R1, R2, Omega1, Omega2)
    R_ = _F(R)
    with np.errstate(divide="ignore", invalid="ignore"):
        u = A * R_ + (np.where(R_ != 0, B / np.where(R_ != 0, R_, 1.0), 0.0) if B != 0.0 else 0.0 * R_)  # Eq. (8.9)
    return (_S(u), A, B) if return_coeffs else _S(u)


def circular_couette_pressure(R, R1, R2, Omega1, Omega2, rho: float = 1000.0, p1: float = 0.0):
    """Pressure of circular Couette flow from the radial (centripetal) balance −u_φ²/R = −(1/ρ) dp/dR.

    Book: §8.2 (R-momentum display before (8.9): "the pressure distribution can therefore be determined once u_φ(R)
    has been found" — the integral is ours): p(R) = p₁ + ρ[A²(R² − R₁²)/2 + 2AB ln(R/R₁) − B²(1/R² − 1/R₁²)/2].

    Parameters
    ----------
    R : [m];  R1, R2 : [m];  Omega1, Omega2 : [rad/s];  rho : ρ [kg/m³];  p1 : pressure at R = R1 [Pa] (at the
    centre when R1 = 0).

    Returns
    -------
    p : [Pa].

    Validation (planned): V1 finite-difference dp/dR = ρu_φ²/R. Label: analytic.
    """
    A, B = _cc_coeffs(R1, R2, Omega1, Omega2)
    R_, R1 = _F(R), float(R1)
    p = A ** 2 / 2.0 * (R_ ** 2 - R1 ** 2)
    if B != 0.0:
        p = p + 2.0 * A * B * np.log(R_ / R1) - B ** 2 / 2.0 * (1.0 / R_ ** 2 - 1.0 / R1 ** 2)
    return _S(float(p1) + float(rho) * p)


def circular_couette_shear_stress(R, R1, R2, Omega1, Omega2, mu: float = 1e-3):
    """Shear stress σ_Rφ = μ[(1/R)∂u_R/∂φ + R ∂(u_φ/R)/∂R] = −2μB/R² of circular Couette flow.

    Book: §8.2, the σ_Rφ display after (8.11) (= −2μΩ₁R₁²/R² for the rotating cylinder in an infinite bath).

    Parameters
    ----------
    R : [m];  R1, R2 : [m];  Omega1, Omega2 : [rad/s];  mu : μ [Pa s].

    Returns
    -------
    sigma_Rphi : [Pa] (stress exerted on the fluid inside radius R by the fluid outside, φ-direction, face normal +e_R).

    Validation (planned): V1 against the printed limit; torque 2πR²σ independent of R. Label: analytic.
    """
    _, B = _cc_coeffs(R1, R2, Omega1, Omega2)
    return _S(-2.0 * float(mu) * B / _F(R) ** 2)


def circular_couette_power(R1, R2, Omega1, Omega2, mu: float = 1e-3) -> dict:
    """Torques, power input and viscous dissipation per unit length of circular Couette flow.

    Book: §8.2 (after (8.11)): "the mechanical power supplied to the fluid … equals the integrated viscous dissipation"
    (Exercise 8.12). Sign (analysis §9 R15): the book writes (2πR₁)τ_Rφ u_φ with σ_Rφ < 0; the power delivered **to the
    fluid** by the inner cylinder is −2πR₁σ_Rφ(R₁)u_φ(R₁) > 0 (the traction on the fluid at R₁ acts through a face of
    outward normal −e_R). Dissipation ε = μ[R d(u_φ/R)/dR]² integrated over R₁ ≤ R ≤ R₂ (``quad``) and in closed form
    4πμB²(1/R₁² − 1/R₂²) = 4πμB(Ω₁ − Ω₂).

    Parameters
    ----------
    R1, R2 : [m] (R2 may be inf);  Omega1, Omega2 : [rad/s];  mu : μ [Pa s].

    Returns
    -------
    dict of floats: torque_inner, torque_outer (torques per unit length exerted on the fluid by each cylinder [N m/m];
    they sum to zero in steady flow); power_inner, power_outer [W/m] (power delivered to the fluid by each cylinder;
    ``power_out`` is the same number as ``power_outer``, the design Part C key); power_in [W/m] (their sum, the net); dissipation [W/m] (closed form); dissipation_quad [W/m] (``quad``); quad_error.

    Validation (planned): V4 power_in = dissipation (quadrature), torque balance; R2 = inf, Ω₂ = 0: 4πμΩ₁²R₁².
    Label: analytic, conserved.
    """
    A, B = _cc_coeffs(R1, R2, Omega1, Omega2)
    R1, R2, O1, O2, mu = float(R1), float(R2), float(Omega1), float(Omega2), float(mu)
    # DEVIATION: power into the fluid = −2πR₁σ_Rφu_φ (> 0); the printed (2πR₁)τ_Rφu_φ is negative (analysis §9 R15)
    t_in = 4.0 * np.pi * mu * B  # = −2πR₁² σ_Rφ(R₁)
    t_out = -t_in  # = +2πR₂² σ_Rφ(R₂)
    p_in, p_out = t_in * O1, t_out * O2
    inv2 = 0.0 if np.isinf(R2) else 1.0 / R2 ** 2
    diss = 4.0 * np.pi * mu * B ** 2 * (1.0 / R1 ** 2 - inv2) if R1 > 0 else (0.0 if B == 0 else np.inf)
    if R1 > 0:
        val, err = quad(lambda r: mu * (2.0 * B / r ** 2) ** 2 * 2.0 * np.pi * r, R1, R2, epsabs=0.0, epsrel=1e-12)
    else:
        val, err = diss, 0.0
    return dict(torque_inner=t_in, torque_outer=t_out, power_inner=p_in, power_outer=p_out, power_out=p_out,
                power_in=p_in + p_out, dissipation=float(diss), dissipation_quad=float(val), quad_error=float(err),
                A=A, B=B)


def circular_couette_state(R1, R2, Omega1, Omega2, mu: float = 1e-3, rho: float = 1000.0) -> dict:
    """Everything the circular-Couette figure/explainer shows (curation §9, backup B1, IF3).

    Book: §8.2, (8.9)–(8.12); Rayleigh's inviscid criterion (stable iff (Ru_φ)² increases outward) is the Ch. 11
    teaser — ours here, not in §8.2.

    Parameters
    ----------
    R1, R2 : [m];  Omega1, Omega2 : [rad/s];  mu : [Pa s];  rho : [kg/m³].

    Returns
    -------
    dict of floats: A [1/s], B [m²/s], vorticity = 2A [1/s], torque_inner, torque_outer [N m/m], power_in,
    dissipation [W/m], dp_gap (alias dp) = p(R2) − p(R1) [Pa] (NaN if R2 = inf), rayleigh_stable (bool: A(AR² + B) ≥ 0 at both walls).

    Validation (planned): V1 against the component functions. Label: analytic.
    """
    pw = circular_couette_power(R1, R2, Omega1, Omega2, mu)
    A, B = pw["A"], pw["B"]
    R1f, R2f = float(R1), float(R2)
    dp = np.nan if np.isinf(R2f) else float(circular_couette_pressure(R2f, R1, R2, Omega1, Omega2, rho, 0.0))
    ends = [R1f] + ([] if np.isinf(R2f) else [R2f])
    phi = [A * (A * r ** 2 + B) for r in ends]  # sign of d(Ru)²/dR ∝ A(AR² + B)
    if np.isinf(R2f):
        phi.append(A * A if A != 0 else 0.0)
    return dict(A=A, B=B, vorticity=2.0 * A, torque_inner=pw["torque_inner"], torque_outer=pw["torque_outer"],
                power_in=pw["power_in"], dissipation=pw["dissipation"], dp=dp, dp_gap=dp,
                rayleigh_stable=bool(min(phi) >= -1e-15 * max(1.0, max(abs(v) for v in phi))))


# ======================================================================================================================
# §8.4 Stokes' first problem and the diffusion family
# ======================================================================================================================
def similarity_variable(y, t, nu: float = 1e-6, half: bool = False):
    """Similarity variable of Stokes' first problem, η = y/√(νt) (``half=True``: y/(2√(νt)), the axis of Figs. 8.13–8.14).

    Book: §8.4, Eq. (8.25). Parameters: y [m]; t [s] (> 0); nu ν [m²/s]. Returns η [–].
    Validation (planned): V1 definition. Label: analytic.
    """
    eta = _F(y) / np.sqrt(float(nu) * _F(t))  # Eq. (8.25)
    return _S(eta / 2.0 if half else eta)


def stokes_first_problem(y, t, U: float = 1.0, nu: float = 1e-6):
    """Stokes' first problem: a plate at y = 0 set moving at U at t = 0, u = U[1 − erf(y/(2√(νt)))] = U erfc(…).

    Book: §8.4, Eq. (8.30) (derived from (8.20)–(8.29)); also §4.6/§4.11, where ch04 used it as a test field (this is
    the promoted implementation — ``ch04.stokes_first_problem`` re-exports it).

    Parameters
    ----------
    y : float or array_like — distance from the plate [m], ≥ 0.
    t : float or array_like — time since the start [s]; for t ≤ 0 the fluid is at rest (8.21) and 0 is returned.
    U : float — plate speed [m/s].
    nu : float — ν [m²/s].

    Returns
    -------
    u : [m/s].

    Assumptions: semi-infinite fluid at rest, ∂p/∂x = 0, constant ν. erfc (not 1 − erf) avoids cancellation at large η.
    Validation (planned): V2 sympy residual of (8.20); V1 BCs, collapse in η; V3 Crank–Nicolson order 2 (I25);
    parity with ch04 ``exact_solution("stokes_first")``. Label: analytic, symbolic, converged.
    """
    y_, t_ = np.broadcast_arrays(_F(y), _F(t))
    with np.errstate(divide="ignore", invalid="ignore"):
        u = np.where(t_ > 0, float(U) * erfc(y_ / (2.0 * np.sqrt(float(nu) * np.where(t_ > 0, t_, 1.0)))), 0.0)  # Eq. (8.30)
    return _S(u)


def stokes_first_vorticity(y, t, U: float = 1.0, nu: float = 1e-6):
    """Vorticity of Stokes' first problem, ω_z = −∂u/∂y = (U/√(πνt)) exp(−y²/4νt).

    Book: §8.4 (text after (8.30): a vortex sheet created at t = 0 diffuses away; ∫₀^∞ ω dy = U — the page prints −U,
    analysis §9 R10). Parameters: y [m]; t [s] (> 0); U [m/s]; nu [m²/s]. Returns ω_z [1/s].
    Validation (planned): V1 −∂u/∂y of (8.30); ∫ω dy = U (``quad``). Label: analytic.
    """
    t_ = _F(t)
    return _S(float(U) / np.sqrt(np.pi * float(nu) * t_) * np.exp(-_F(y) ** 2 / (4.0 * float(nu) * t_)))


def diffusion_thickness(t, nu: float = 1e-6, level: float = 0.01):
    """Thickness of the diffusive layer of Stokes' first problem: where u/U falls to ``level``.

    Book: §8.4, Eq. (8.31) δ₉₉ ~ 3.64√(νt) (level 0.01); here exactly δ = 2 erfc⁻¹(level)√(νt).

    Parameters
    ----------
    t : time [s];  nu : ν [m²/s];  level : u/U at the edge (0 < level < 1).

    Returns
    -------
    delta : [m].

    Validation (planned): V1 erfc(δ/2√(νt)) = level; brentq cross-check; V6 the printed 3.64 (private JSON).
    Label: analytic, book-value.
    """
    if not 0.0 < float(level) < 1.0:
        raise ValueError("level must lie in (0, 1)")
    return _S(2.0 * erfcinv(float(level)) * np.sqrt(float(nu) * _F(t)))  # Eq. (8.31), exact coefficient


def transition_width(t, nu: float = 1e-6, level: float = 0.95):
    """Width of the thickening vortex sheet between the points u = ±level·U: 4 erf⁻¹(level)√(νt).

    Book: §8.4, Example 8.5 (u = ±0.95U; the page prints η = ±2.76 but the width 5.54√(νt) = 2 × 2.772 — rounding slip,
    analysis §9 R11). Parameters: t [s]; nu [m²/s]; level (0, 1). Returns width [m].
    Validation (planned): V1 erf(w/(4√(νt))) = level. Label: analytic.
    """
    if not 0.0 < float(level) < 1.0:
        raise ValueError("level must lie in (0, 1)")
    return _S(4.0 * erfinv(float(level)) * np.sqrt(float(nu) * _F(t)))


def stokes_first_stopped(y, t, T: float, U: float = 1.0, nu: float = 1e-6):
    """Plate started at U at t = 0 and stopped at t = T (Exercise 8.30), by superposition (ours).

    Book: §8.4 (text after (8.31): stopping the plate imposes a time scale and destroys similarity; Exercise 8.30).
    u = U erfc(y/2√(νt)) for 0 < t ≤ T and U[erfc(y/2√(νt)) − erfc(y/2√(ν(t − T)))] for t > T.

    Parameters
    ----------
    y : [m];  t : [s];  T : stopping time [s] (> 0);  U : [m/s];  nu : [m²/s].

    Returns
    -------
    u : [m/s].

    Validation (planned): V1 equals (8.30) for t ≤ T, u(0, t > T) = 0, diffusion residual 0. Label: analytic.
    """
    y_, t_ = np.broadcast_arrays(_F(y), _F(t))
    u = _F(stokes_first_problem(y_, t_, U, nu))
    late = _F(stokes_first_problem(y_, t_ - float(T), U, nu))
    return _S(np.where(t_ > float(T), u - late, u))


def stokes_first_state(t, U: float = 1.0, nu: float = 1e-6, level: float = 0.01, rho: float = 1000.0) -> dict:
    """Numbers the Stokes-first-problem explainer shows at time t (curation §9, E5).

    Book: §8.4, (8.30)–(8.31), Example 8.5 (τ_w).

    Parameters
    ----------
    t : [s];  U : [m/s];  nu : [m²/s];  level : u/U defining the edge;  rho : ρ [kg/m³] (μ = ρν for τ_w).

    Returns
    -------
    dict of floats: sqrt_nut = √(νt) [m], delta = 2 erfc⁻¹(level)√(νt) [m], eta_edge = δ/√(νt) [–], tau_w = μU/√(πνt)
    [Pa] (stress on the plate from the fluid is −τ_w), vorticity_content = ∫₀^∞ ω dy = U [m/s] (constant in time),
    omega_wall = U/√(πνt) [1/s].

    Validation (planned): V1 against the component functions. Label: analytic.
    """
    t, U, nu = float(t), float(U), float(nu)
    s = np.sqrt(nu * t)
    d = float(diffusion_thickness(t, nu, level))
    return dict(sqrt_nut=s, delta=d, eta_edge=d / s, tau_w=float(rho) * nu * U / np.sqrt(np.pi * nu * t),
                vorticity_content=U, omega_wall=U / np.sqrt(np.pi * nu * t))


def vortex_sheet_diffusion(y, t, U: float = 1.0, nu: float = 1e-6):
    """Viscous thickening of a vortex sheet that separates u = +U (y > 0) from u = −U (y < 0).

    Book: §8.4, Example 8.5: ω_z = −(U/√(πνt)) exp(−y²/4νt), u = U erf(y/2√(νt)) (from (8.32a) with n = ½, δ = √(νt)).
    Parity: ch05 ``diffusing_vortex_sheet(y, t, gamma=−2U, nu)`` (γ = u_below − u_above).

    Parameters
    ----------
    y : [m];  t : [s] (> 0);  U : [m/s];  nu : [m²/s].

    Returns
    -------
    (u, omega_z) : [m/s], [1/s].

    Validation (planned): V2 diffusion residual; V4 ∫ω dy = −2U for all t; parity ch05 with γ = −2U. Label: analytic,
    symbolic, conserved.
    """
    y_, t_ = _F(y), _F(t)
    s = np.sqrt(float(nu) * t_)
    w = -float(U) / (np.sqrt(np.pi) * s) * np.exp(-y_ ** 2 / (4.0 * s ** 2))  # Example 8.5, ω_z
    u = float(U) * erf(y_ / (2.0 * s))  # Example 8.5, u
    return _S(u), _S(w)


def temporal_bl_wall_stress(t, U: float = 1.0, nu: float = 1e-6, rho: float = 1000.0) -> dict:
    """Wall stress and skin-friction coefficient of the temporally developing boundary layer (Example 8.5).

    Book: §8.4, Example 8.5: τ_w = μ(∂u/∂y)_{y=0} = μU/√(πνt), C_f = τ_w/(½ρU²) = (2/√π)√(ν/(U²t)); with Ut → x,
    C_f = (2/√π) Re_x^{−1/2} (compare Blasius 0.664 Re_x^{−1/2}, Ch. 9).

    Parameters
    ----------
    t : [s];  U : [m/s];  nu : [m²/s];  rho : [kg/m³].

    Returns
    -------
    dict of floats: tau_w [Pa], Cf [–], Re_x (alias Rex) = U²t/ν [–], Cf_coefficient = 2/√π, Cf_of_Rex = (2/√π)Re_x^{−1/2}.

    Validation (planned): V1 finite-difference ∂u/∂y at the wall of :func:`vortex_sheet_diffusion`. Label: analytic.
    """
    t, U, nu, rho = float(t), float(U), float(nu), float(rho)
    tau = rho * nu * U / np.sqrt(np.pi * nu * t)
    Cf = tau / (0.5 * rho * U ** 2)
    Rex = U ** 2 * t / nu
    return dict(tau_w=tau, Cf=Cf, Re_x=Rex, Rex=Rex, Cf_coefficient=2.0 / np.sqrt(np.pi),
                Cf_of_Rex=2.0 / np.sqrt(np.pi) / np.sqrt(Rex))


def line_vortex_decay(r, t, Gamma: float = 1.0, nu: float = 1e-6):
    """Viscous decay of a line vortex (thin spinning cylinder stopped at t = 0): the Lamb–Oseen vortex.

    Book: §8.4, Example 8.6: u_θ = (Γ/2πr)[1 − exp(−r²/4νt)] (= Gaussian vortex (3.29) with σ² = 4νt).
    Evaluated as −expm1(−r²/4νt)/r, finite at r → 0 (limit Γr/(8πνt), 0 on the axis).

    Parameters
    ----------
    r : radius [m], ≥ 0;  t : [s] (> 0);  Gamma : Γ [m²/s];  nu : [m²/s].

    Returns
    -------
    u_theta : [m/s].

    Validation (planned): V2 residual of ∂u/∂t = ν∂/∂r[(1/r)∂(ru)/∂r]; parity ``core.vortices.gaussian_vortex(σ=2√(νt))``
    and ch04 ``exact_solution("lamb_oseen")``; V4 circulation → Γ at large r. Label: analytic, symbolic.
    """
    r_, t_ = np.broadcast_arrays(_F(r), _F(t))
    a2 = 4.0 * float(nu) * t_
    with np.errstate(divide="ignore", invalid="ignore"):
        u = np.where(r_ > 0, float(Gamma) / (2.0 * np.pi) * (-np.expm1(-r_ ** 2 / a2)) / np.where(r_ > 0, r_, 1.0),
                     0.0)  # Example 8.6
    return _S(u)


def line_vortex_spinup(r, t, Gamma: float = 1.0, nu: float = 1e-6):
    """Line vortex suddenly introduced into fluid at rest (impulsive spin-up of a thin cylinder).

    Book: §8.4, end of Example 8.6 and Exercise 8.26: u_θ = (Γ/2πr) exp(−r²/4νt) (singular at r = 0: the imposed line
    vortex; returns inf there).

    Parameters
    ----------
    r : [m];  t : [s] (> 0);  Gamma : [m²/s];  nu : [m²/s].

    Returns
    -------
    u_theta : [m/s].

    Validation (planned): V2 diffusion residual; circulation at fixed r → Γ as t → ∞, → 0 as t → 0. Label: analytic,
    symbolic.
    """
    r_, t_ = np.broadcast_arrays(_F(r), _F(t))
    with np.errstate(divide="ignore"):
        u = float(Gamma) / (2.0 * np.pi * r_) * np.exp(-r_ ** 2 / (4.0 * float(nu) * t_))  # Exercise 8.26
    return _S(u)


# ======================================================================================================================
# §8.5 Stokes' second problem
# ======================================================================================================================
def stokes_second_problem(y, t, U: float = 1.0, omega: float = 2.0 * np.pi, nu: float = 1e-6):
    """Periodic flow above a plate oscillating in its own plane, u(0, t) = U cos ωt (Stokes' second problem).

    Book: §8.5, Eq. (8.38) (from (8.35)–(8.37): k = ±(1 + i)√(ω/2ν), bounded root, A = U).

    Parameters
    ----------
    y : distance from the plate [m];  t : time [s];  U : velocity amplitude [m/s];  omega : ω [rad/s];  nu : [m²/s].

    Returns
    -------
    u : [m/s] (the periodic state after transients; no initial condition).

    Validation (planned): V2 residual of (8.20); V1 u(0, t) = U cos ωt; V3 Crank–Nicolson after 10 periods.
    Label: analytic, symbolic, converged.
    """
    k = np.sqrt(float(omega) / (2.0 * float(nu)))
    y_ = _F(y)
    return _S(float(U) * np.exp(-y_ * k) * np.cos(float(omega) * _F(t) - y_ * k))  # Eq. (8.38)


def stokes_layer(nu: float = 1e-6, omega: float = 2.0 * np.pi) -> dict:
    """Length and speed scales of the Stokes layer (8.38).

    Book: §8.5 (after (8.38)): δ ~ 4(ν/ω)^{1/2}, amplitude U e^{−4/√2} ≈ 0.06U there. Ours: δ_e = √(2ν/ω) (e-folding
    depth, the literature's penetration depth), phase (crest) speed ω/k = √(2νω), wavelength 2π√(2ν/ω).

    Parameters
    ----------
    nu : ν [m²/s];  omega : ω [rad/s].

    Returns
    -------
    dict of floats: k = √(ω/2ν) [1/m], delta_e [m], delta_book [m], ratio_book_to_e = 2√2, phase_speed [m/s],
    wavelength [m], amplitude_at_delta_book = e^{−4/√2} [–].

    Validation (planned): V1 closed forms; V3 phase speed from zero-crossing tracking. Label: analytic.
    """
    nu, omega = float(nu), float(omega)
    require_positive("omega", omega)
    k = np.sqrt(omega / (2.0 * nu))
    return dict(k=k, delta_e=1.0 / k, delta_book=4.0 * np.sqrt(nu / omega), ratio_book_to_e=2.0 * np.sqrt(2.0),
                phase_speed=omega / k, wavelength=2.0 * np.pi / k, amplitude_at_delta_book=float(np.exp(-4.0 / np.sqrt(2.0))))


def stokes_layer_envelope(y, U: float = 1.0, omega: float = 2.0 * np.pi, nu: float = 1e-6):
    """Envelope ±U exp(−y√(ω/2ν)) of the Stokes-layer profiles (8.38). Book: §8.5. y [m] → amplitude [m/s].
    Label: analytic."""
    return _S(float(U) * np.exp(-_F(y) * np.sqrt(float(omega) / (2.0 * float(nu)))))


def stokes_layer_state(nu: float, omega: float, y, U: float = 1.0) -> dict:
    """Numbers the Stokes-layer explainer shows at height y (curation §9, E7).

    Book: §8.5, (8.38).

    Parameters
    ----------
    nu : [m²/s];  omega : [rad/s];  y : height [m];  U : wall amplitude [m/s].

    Returns
    -------
    dict of floats: delta_e, delta_book [m]; amplitude = U e^{−y/δ_e} [m/s]; amp_at_book_depth = e^{−2√2} [–] (the
    relative amplitude at y = 4√(ν/ω)); phase_lag = y/δ_e [rad]; time_lag = phase_lag/ω [s]; crest_speed = √(2νω)
    [m/s]; wavelength [m]; period = 2π/ω [s]; y_over_delta_e [–].

    Validation (planned): V1 against :func:`stokes_second_problem` (maximum over a period at y). Label: analytic.
    """
    s = stokes_layer(nu, omega)
    y = float(y)
    lag = y / s["delta_e"]
    return dict(delta_e=s["delta_e"], delta_book=s["delta_book"], amplitude=float(U) * np.exp(-lag),
                amp_at_book_depth=s["amplitude_at_delta_book"], phase_lag=lag, time_lag=lag / float(omega),
                crest_speed=s["phase_speed"], wavelength=s["wavelength"], period=2.0 * np.pi / float(omega),
                y_over_delta_e=lag)
