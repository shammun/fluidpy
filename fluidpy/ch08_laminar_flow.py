"""Chapter 8 — Laminar flow: exact steady solutions (plane Couette–Poiseuille, pipe, circular Couette), elementary
lubrication theory (slider bearing, Hele-Shaw, spreading bead), similarity solutions (Stokes' first problem, vortex
sheet, line vortex, the spreading bead), Stokes' second problem, and creeping flow past a sphere (Stokes, Oseen) —
plus re-exports of the reusable primitives in ``fluidpy.core.laminar``, ``core.lubrication``, ``core.creeping`` and
``core.diffusion.crank_nicolson_1d``.

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e (2012), Ch. 8, §§8.1–8.7, Eqs. (8.1)–(8.53), Examples 8.1–8.7.
Every equation was transcribed from the rendered page images (chapters/pages/ch08/p337–p374, printed pp. 310–347).

Where the physics lives (all public names are re-exported, so ``ch08.<name>`` reaches every callable)
* ``core.laminar`` (LAM) — (8.4)–(8.12), (8.20)–(8.31), Examples 8.5–8.6, (8.33)–(8.38), the ``*_state`` dicts.
* ``core.lubrication`` (LUB) — (8.14)–(8.19), Examples 8.1–8.3, 8.7.
* ``core.creeping`` (CRP) — (8.43)–(8.53), drag laws, settling, Millikan, Oseen.
* ``core.diffusion.crank_nicolson_1d`` — our implicit scheme for (8.20) with (8.22) or (8.33).
* this module — §8.1 helpers, the sympy derivation engines (``*_sympy``), the similarity engine and collapse measure,
  the synthetic Millikan experiment.

Conventions (analysis §9): walls y = 0 (fixed) and y = h (moving); ``dpdx`` = dp/dx (book sign) with ``G=`` = −dp/dx
as a keyword alias (ch04 parity); η = y/√(νt); §8.6 θ from the downstream axis, Re = 2aU/ν. Traps for earlier
chapters' functions: ch05 ``rotating_cylinder_flow(r, a, omega)`` takes ω = 2Ω₁ (the cylinder's vorticity); ch05
``diffusing_vortex_sheet(y, t, gamma, nu)`` takes γ = u_below − u_above (Example 8.5 is γ = −2U).

Book slips (coded corrected; printed forms kept as labelled options): R6 (8.13b) ∂p/∂x → ∂p/∂y; R7 (8.17a) and
Example 8.2 lack ν; R8 Example 8.1's intermediate (1 − αx/L) and final first-power denominator; R8b (8.19) U₀
placement; R9 V without 1/h; R10 ∫ω dy = −U → +U; R11 η = ±2.76 → ±2.772; R12 minimum pressure −3μU/2a; R13 Oseen
equation −∂p/∂x_i; R14 cross-references "(9.63)" → (8.43), "(9.68)" → (8.48), "(8.33) into (8.20)" → (8.35); R15 power
sign. Book-quoted numbers live only in the git-ignored ``tests/book_values_ch08.json``.
"""
from __future__ import annotations

import functools

import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_bvp
from scipy.special import erf, erfc

from . import ch01_introduction as _ch01
from .core import creeping as CRP  # noqa: F401
from .core import curvilinear as CL
from .core import laminar as LAM  # noqa: F401
from .core import lubrication as LUB  # noqa: F401
from .core import similarity as SIM
from .core._util import as_scalar_if_0d
from .core.creeping import *  # noqa: F401,F403
from .core.diffusion import couette_startup_profile, crank_nicolson_1d, ftcs_diffusion_1d  # noqa: F401
from .core.dimensional import pi_groups
from .core.laminar import *  # noqa: F401,F403
from .core.laminar import line_vortex_decay, stokes_first_problem, stokes_first_vorticity, vortex_sheet_diffusion
from .core.lubrication import *  # noqa: F401,F403
from .core.lubrication import viscous_current_similarity
from .core.similarity import reynolds_number  # noqa: F401
from .core.thermo import G0, G_BOOK, P_ATM, R_AIR  # noqa: F401  (re-exported constants)

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d
E_CHARGE = 1.602176634e-19
"""Elementary charge [C] (exact, SI 2019; NIST CODATA 2022). Label: benchmark."""


def _z(e) -> bool:
    return sp.simplify(e) == 0


# ======================================================================================================================
# §8.1
# ======================================================================================================================
def pipe_flow_regime(U, d, nu: float = 1e-6):
    """Pipe Reynolds number Re = Ud/ν and the regime of Reynolds's dye experiment.

    Book: §8.1 (transition at a fixed Re = Ud/ν ~ 2000–3000; U mean velocity, d diameter). Labels: "laminar" for
    Re < 2000, "transitional" for 2000 ≤ Re ≤ 3000, "turbulent" above 3000 (the band is the book's; the precise onset
    depends on disturbances).

    Parameters
    ----------
    U : mean velocity [m/s];  d : diameter [m];  nu : ν [m²/s].

    Returns
    -------
    (Re, label) : Re [–] and the regime string (array of strings for array input).

    Validation (planned): V1 definition, labels switch at 2000 and 3000. Label: analytic.
    """
    Re = _F(U) * _F(d) / float(nu)
    lab = np.where(Re < 2000.0, "laminar", np.where(Re <= 3000.0, "transitional", "turbulent"))
    return _S(Re), (str(lab) if lab.ndim == 0 else lab)


def inertia_viscous_scales(U, L, nu: float = 1e-6) -> dict:
    """Order-of-magnitude inertial and viscous accelerations and their ratio, the Reynolds number.

    Book: §8.1 ("inertia ~U²/L ≫ viscous ~μU/ρL²" ⇔ Re ≫ 1; the scaling of (4.100)). Parameters: U [m/s]; L [m];
    nu [m²/s]. Returns dict(inertia = U²/L, viscous = νU/L² [m/s²], ratio = UL/ν = Re). Label: analytic.
    """
    U, L, nu = float(U), float(L), float(nu)
    return dict(inertia=U ** 2 / L, viscous=nu * U / L ** 2, ratio=U * L / nu)


def momentum_diffusivity(fluid: str = "air", T: float = 293.15, p: float = P_ATM) -> float:
    """Kinematic viscosity ν = μ/ρ (the momentum diffusivity) of air or water.

    Book: §8.1 (ν of air ≈ 15× that of water at room conditions although μ_water > μ_air; ν plays the role of κ in
    (4.89)). Air: Sutherland μ(T) (ch01, USSA-1976) and ρ = p/(R_air T); water: ch01 Vogel μ(T) and Kell ρ(T).

    Parameters
    ----------
    fluid : "air" | "water";  T : temperature [K];  p : pressure [Pa] (air only).

    Returns
    -------
    nu : [m²/s].

    Validation (planned): V5 via the ch01 correlations; ratio air/water at 293 K in [14, 17]. Label: benchmark.
    """
    if fluid == "air":
        return float(_ch01.sutherland_viscosity(T)) / (float(p) / (R_AIR * float(T)))
    if fluid == "water":
        return float(_ch01.water_viscosity(T)) / float(_ch01.water_density(T))
    raise ValueError('fluid must be "air" or "water"')


def diffusion_time(L, nu: float = 1e-6):
    """Time for momentum (vorticity) to diffuse a distance L, t ~ L²/ν. Book: §8.1, §8.4 (δ ∝ √(νt)).
    Parameters: L [m]; nu [m²/s]. Returns t [s]. Label: analytic."""
    return _S(_F(L) ** 2 / float(nu))


def wall_bc_residuals(u_wall, U_s, n, t=None):
    """Residuals of the wall conditions: no through-flow n·U_s = n·u (8.2) and no slip t·U_s = t·u (8.3).

    Book: §8.1, Eqs. (8.2)–(8.3). Parameters: u_wall fluid velocity at the wall (3,) [m/s]; U_s wall velocity (3,);
    n unit normal (3,); t optional unit tangent (3,) — if None the tangential residual is the magnitude of the whole
    tangential part of u − U_s. Returns (normal, tangential) [m/s] (both 0 when the conditions hold).
    Validation (planned): V1 zero for u = U_s; separate parts for slip/through-flow cases. Label: analytic.
    """
    d = _F(u_wall) - _F(U_s)
    nn = _F(n) / np.linalg.norm(_F(n))
    normal = float(nn @ d)  # Eq. (8.2)
    if t is not None:
        tt = _F(t) / np.linalg.norm(_F(t))
        return normal, float(tt @ d)  # Eq. (8.3)
    return normal, float(np.linalg.norm(d - normal * nn))


# ======================================================================================================================
# §8.2 derivations in sympy
# ======================================================================================================================
def _pf_out(d: dict) -> dict:
    """Add the summary keys ``residual`` (sum of |residuals|, 0) and ``checks`` to a derivation dict."""
    vals = list(d["residuals"].values()) + list(d.get("bc_residuals", []))
    d["residual"] = sp.simplify(sum(sp.Abs(sp.simplify(v_)) for v_ in vals))
    d["checks"] = dict(all_zero=all(_z(v_) for v_ in vals))
    return d


def parallel_flow_sympy(case: str = "channel") -> dict:
    """The §8.2 derivations step by step in sympy, each with a residual that simplifies to 0.

    Book: §8.2 — "channel": continuity ⇒ v = 0, (8.4a,b), the separation "function of x = function of y ⇒ constant",
    the double integral with A, B, (8.5), Q; "pipe": cylindrical momentum (Appendix B via ``core.curvilinear``),
    u_z = R²(dp/dz)/4μ + A ln R + B, A = 0 (bounded), (8.6), (8.7), Q; "circular_couette": radial and azimuthal
    momentum, (8.9), constants, (8.10), limits (8.11)–(8.12), pressure (ours) and σ_Rφ.

    Parameters
    ----------
    case : "channel" | "pipe" | "circular_couette".

    Returns
    -------
    dict of sympy objects with (design Part C 4.4) continuity, momentum (reduced equations), separation (the
    "constant" lines d²p/dx² = 0 and u‴ = 0; None for circular Couette), general (integrated profile with constants),
    constants, profile, residual (sum of |residuals|, 0), bc_residuals (list of 0); pipe adds unbounded_term (A ln R);
    circular Couette adds A, B, pressure; plus Q, V, tau, sigma, limits, a sub-dict ``residuals`` (each 0) and
    ``checks`` (booleans).

    Validation (planned): V2 every residual is 0; wrong variant: keeping A ln R leaves u unbounded. Label: symbolic.
    """
    if case == "channel":
        x, y = sp.symbols("x y", real=True)
        h, U, mu, rho = sp.symbols("h U mu rho", positive=True)
        G, A, B = sp.symbols("G A B", real=True)  # G = dp/dx here (the book's symbol dp/dx)
        u, p, f = sp.Function("u")(y), sp.Function("p")(x), sp.Function("f")
        nu = mu / rho
        mom_x = sp.Eq(0, -sp.diff(p, x) / rho + nu * sp.diff(u, y, 2))  # (8.4a)
        mom_y = sp.Eq(0, -sp.diff(p, y) / rho)  # (8.4b): identically 0 once p = p(x)
        sep = sp.diff(-sp.diff(p, x) / rho + nu * sp.diff(u, y, 2), x)  # d/dx: −p''(x)/ρ must vanish
        nonconst = sp.diff(-f(x) / rho + nu * sp.diff(u, y, 2), x)  # dp/dx = f(x): requires f'(x) = 0
        integ = sp.Eq(0, -y ** 2 / 2 * G + mu * sp.Symbol("u") + A * y + B)
        u_gen = sp.solve(integ, sp.Symbol("u"))[0]
        cons = sp.solve([u_gen.subs(y, 0), u_gen.subs(y, h) - U], [A, B], dict=True)[0]
        prof = sp.simplify(u_gen.subs(cons))
        book = U / h * y - G / (2 * mu) * y * (h - y)  # (8.5)
        Q = sp.simplify(sp.integrate(prof, (y, 0, h)))
        Qb = U * h / 2 * (1 - h ** 2 / (6 * mu * U) * G)
        tau = sp.simplify(mu * sp.diff(prof, y))
        res = dict(momentum=sp.simplify(-G + mu * sp.diff(prof, y, 2)), bc0=prof.subs(y, 0),
                   bch=sp.simplify(prof.subs(y, h) - U), vs_book=sp.simplify(prof - book), Q=sp.simplify(Q - Qb),
                   A_book=sp.simplify(cons[A] - (h / 2 * G - mu * U / h)), B_book=cons[B],
                   tau_poiseuille=sp.simplify(tau.subs(U, 0) + (h / 2 - y) * G))
        vfun = sp.Function("v")(y)
        return _pf_out(dict(case=case, symbols=(x, y, h, U, mu, rho, G),
                            continuity=sp.Eq(sp.diff(u, x) + sp.diff(vfun, y), 0),
                            continuity_text="∂u/∂x = 0 ⇒ ∂v/∂y = 0, v(0) = 0 ⇒ v ≡ 0",
                            momentum=(mom_x, mom_y), momentum_x=mom_x, momentum_y=mom_y,
                            separation=(sp.Eq(sp.diff(p, x, 2), 0), sp.Eq(sp.diff(u, y, 3), 0)),
                            separation_derivative=sep, nonconstant_requires=nonconst, general=u_gen, integrated=integ,
                            constants=cons, profile=prof, Q=Q, V=sp.simplify(Q / h), tau=tau, residuals=res,
                            bc_residuals=[res["bc0"], res["bch"]]))
    if case == "pipe":
        R, ph, z = CL.coordinates("cylindrical")
        a, mu = sp.symbols("a mu", positive=True)
        G, C1, C2 = sp.symbols("G C1 C2", real=True)  # G = dp/dz
        w = sp.Function("w")(R)
        uvec = [0, 0, w]
        div = CL.divergence(uvec, "cylindrical")
        lap = CL.vector_laplacian(uvec, "cylindrical")
        adv = CL.advective_acceleration(uvec, "cylindrical")
        zmom = sp.Eq(0, -G + mu * lap[2])
        gen = sp.dsolve(zmom, w).rhs
        consts = sorted(gen.free_symbols - {R, a, mu, G}, key=str)
        # identify the ln R coefficient and the constant
        gen2 = G * R ** 2 / (4 * mu) + C1 * sp.log(R) + C2
        chk_gen = sp.simplify(-G + mu * CL.laplacian(gen2, "cylindrical"))
        bounded = gen2.subs(C1, 0)
        c2 = sp.solve(bounded.subs(R, a), C2)[0]
        prof = sp.simplify(bounded.subs(C2, c2))
        book = (R ** 2 - a ** 2) / (4 * mu) * G  # (8.6)
        tau = sp.simplify(mu * sp.diff(prof, R))
        Q = sp.simplify(sp.integrate(prof * 2 * sp.pi * R, (R, 0, a)))
        res = dict(divergence=div, R_momentum=sp.simplify(lap[0]), phi_momentum=sp.simplify(lap[1]),
                   advective=sum(sp.Abs(sp.simplify(c)) for c in adv), general=chk_gen,
                   vs_book=sp.simplify(prof - book), tau=sp.simplify(tau - R / 2 * G),
                   Q=sp.simplify(Q + sp.pi * a ** 4 / (8 * mu) * G), wall=prof.subs(R, a))
        return _pf_out(dict(case=case, symbols=(R, a, mu, G), continuity=sp.Eq(div, 0),
                            momentum=(sp.Eq(0, -sp.Symbol("dp/dR") + mu * lap[0]),
                                      sp.Eq(0, -sp.Symbol("dp/dphi") / R + mu * lap[1]), zmom),
                            separation=(sp.Eq(sp.Symbol("d2p/dz2"), 0),
                                        sp.Eq(sp.diff(sp.diff(R * w.diff(R), R) / R, R), 0)),
                            dsolve_general=gen, dsolve_constants=consts, general=gen2,
                            unbounded_term=C1 * sp.log(R), ln_term_unbounded=sp.limit(sp.log(R), R, 0, "+"),
                            constants={C1: 0, C2: c2}, B=c2, profile=prof, tau=tau, tau_wall=sp.simplify(tau.subs(R, a)),
                            Q=Q, V=sp.simplify(Q / (sp.pi * a ** 2)), residuals=res, bc_residuals=[res["wall"]]))
    if case == "circular_couette":
        R, ph, z = CL.coordinates("cylindrical")
        R1, R2, rho, mu = sp.symbols("R_1 R_2 rho mu", positive=True)
        O1, O2, A, B, p1 = sp.symbols("Omega_1 Omega_2 A B p_1", real=True)
        v = sp.Function("v")(R)
        uvec = [0, v, 0]
        adv = CL.advective_acceleration(uvec, "cylindrical")
        lap = CL.vector_laplacian(uvec, "cylindrical")
        phimom = sp.Eq(0, mu * lap[1])
        gen = sp.dsolve(phimom, v).rhs
        g2 = A * R + B / R  # (8.9)
        chk = sp.simplify(CL.vector_laplacian([0, g2, 0], "cylindrical")[1])
        cons = sp.solve([g2.subs(R, R1) - O1 * R1, g2.subs(R, R2) - O2 * R2], [A, B], dict=True)[0]
        Ab = (O2 * R2 ** 2 - O1 * R1 ** 2) / (R2 ** 2 - R1 ** 2)
        Bb = -(O2 - O1) * R1 ** 2 * R2 ** 2 / (R2 ** 2 - R1 ** 2)
        prof = g2.subs(cons)
        book = ((O2 * R2 ** 2 - O1 * R1 ** 2) * R - (O2 - O1) * R1 ** 2 * R2 ** 2 / R) / (R2 ** 2 - R1 ** 2)  # (8.10)
        lim1 = sp.limit(prof.subs(O2, 0), R2, sp.oo)  # (8.11)
        lim2 = sp.simplify(prof.subs({O1: 0}).subs(R1, 0))  # (8.12)
        pr = p1 + rho * sp.integrate((g2 ** 2 / R).subs(R, sp.Symbol("s", positive=True)),
                                     (sp.Symbol("s", positive=True), R1, R))
        ours = p1 + rho * (A ** 2 / 2 * (R ** 2 - R1 ** 2) + 2 * A * B * sp.log(R / R1) - B ** 2 / 2 * (1 / R ** 2 - 1 / R1 ** 2))
        sig = sp.simplify(mu * R * sp.diff(g2 / R, R))
        res = dict(radial_advective=sp.simplify(adv[0] + v ** 2 / R), phi_general=chk,
                   A=sp.simplify(cons[A] - Ab), B=sp.simplify(cons[B] - Bb), vs_book=sp.simplify(prof - book),
                   limit_8_11=sp.simplify(lim1 - O1 * R1 ** 2 / R), limit_8_12=sp.simplify(lim2 - O2 * R),
                   pressure=sp.simplify(sp.expand_log(pr - ours, force=True)), sigma=sp.simplify(sig + 2 * mu * B / R ** 2),
                   vorticity=sp.simplify(CL.curl([0, g2, 0], "cylindrical")[2] - 2 * A))
        prof_s = sp.simplify(prof)
        return _pf_out(dict(case=case, symbols=(R, R1, R2, O1, O2, rho, mu),
                            continuity=sp.Eq(CL.divergence(uvec, "cylindrical"), 0),
                            momentum=(sp.Eq(-v ** 2 / R, -sp.Symbol("dpdR") / rho), phimom),
                            radial_balance=sp.Eq(-v ** 2 / R, -sp.Symbol("dpdR") / rho),
                            separation=None, dsolve_general=gen, general=g2, constants=cons, A=cons[A], B=cons[B],
                            profile=prof_s, limit_R2_inf=lim1, limit_R1_zero=lim2, pressure=ours, sigma=sig,
                            residuals=res, bc_residuals=[sp.simplify(prof_s.subs(R, R1) - O1 * R1),
                                                         sp.simplify(prof_s.subs(R, R2) - O2 * R2)]))
    raise ValueError('case must be "channel", "pipe" or "circular_couette"')


def advective_acceleration_check(case: str = "channel") -> dict:
    """(u·∇)u for the three §8.2 flows (``core.curvilinear``): 0 for the channel and the pipe, −u_φ²/R e_R for circular
    Couette flow (balanced by the pressure gradient). Book: §8.2 (end: symmetry removes the nonlinear term).
    Returns dict(advective (list of 3 sympy), expected, residual (0)). Label: symbolic."""
    if case == "channel":
        x, y, z = CL.coordinates("cartesian")
        adv = CL.advective_acceleration([sp.Function("u")(y), 0, 0], "cartesian")
        exp_ = [0, 0, 0]
    elif case == "pipe":
        R, ph, z = CL.coordinates("cylindrical")
        adv = CL.advective_acceleration([0, 0, sp.Function("w")(R)], "cylindrical")
        exp_ = [0, 0, 0]
    elif case == "circular_couette":
        R, ph, z = CL.coordinates("cylindrical")
        v = sp.Function("v")(R)
        adv = CL.advective_acceleration([0, v, 0], "cylindrical")
        exp_ = [-v ** 2 / R, 0, 0]
    else:
        raise ValueError('case must be "channel", "pipe" or "circular_couette"')
    resid = sum(sp.Abs(sp.simplify(a_ - e_)) for a_, e_ in zip(adv, exp_))
    return dict(advective=adv, expected=exp_, residual=sp.simplify(resid))


# ======================================================================================================================
# §8.3 lubrication in sympy
# ======================================================================================================================
def _coeff(term, scale):
    e = sp.simplify(term / scale)
    e = e.replace(lambda q: isinstance(q, (sp.Derivative, sp.Subs)), lambda q: sp.Integer(1))
    e = e.replace(lambda q: isinstance(q, sp.core.function.AppliedUndef), lambda q: sp.Integer(1))
    return sp.simplify(e)


def lubrication_nondim_sympy(printed_8_13b: bool = False) -> dict:
    """Cached wrapper of :func:`_lubrication_nondim_sympy` (a fresh dict each call). Book: §8.3, Eqs. (8.13)–(8.19).
    Label: symbolic."""
    return dict(_lubrication_nondim_sympy(bool(printed_8_13b)))


@functools.lru_cache(maxsize=2)
def _lubrication_nondim_sympy(printed_8_13b: bool = False) -> dict:
    """The thin-gap scaling (8.14) substituted into (6.2) and (8.13a,b) — coefficients of every term, the ε²Re_L → 0
    limit (8.17) and the profile (8.18)–(8.19), in sympy.

    Book: §8.3, Eqs. (8.13)–(8.19). The y-momentum equation uses ∂p/∂y (the page prints ∂p/∂x, analysis §9 R6;
    ``printed_8_13b=True`` keeps the printed term, which gives a pressure coefficient ε/Λ and a ∂p*/∂x* derivative —
    inconsistent with (8.16b)). The limit is returned with ν (the page's (8.17a) lacks it, R7; ``printed_8_17a`` shows
    the printed form, whose two terms have different units).

    Returns
    -------
    dict (design Part C 4.6 keys first): continuity_star ((8.15) as an Eq), x_coeffs, y_coeffs (dicts inertia,
    pressure, diff_along, diff_across), printed_13b_y_coeffs (the same with the printed ∂p/∂x: pressure ε/Λ — differs),
    limit_x ((8.17a) with ν), limit_y ((8.17b)), eq_8_18; and in more detail: symbols; continuity_coefficients (both 1
    after dividing by U/L); x_coefficients, y_coefficients (dicts of
    sympy expressions in eps, Re_L, Lambda: unsteady, advect_u, advect_v, pressure, diff_xx, diff_yy); expected_x,
    expected_y (the (8.16a,b) sets); matches (booleans); limit_8_17a (dimensional, with ν), printed_8_17a;
    profile_8_18, profile_8_19 (consistent), profile_8_19_book.

    Validation (planned): V2 coefficient sets {ε²Re_L, 1/Λ, ε², 1} and {ε⁴Re_L, 1/Λ, ε⁴, ε²}. Label: symbolic.
    """
    L, h, U, rho, mu, Pa = sp.symbols("L h U rho mu P_a", positive=True)
    eps, Re, Lam = sp.symbols("epsilon Re_L Lambda", positive=True)
    x, y, t = sp.symbols("x y t", real=True)
    us, vs, ps = sp.Function("u_s"), sp.Function("v_s"), sp.Function("p_s")
    args = (x / L, y / h, U * t / L)
    u = U * us(*args)
    v = (h / L) * U * vs(*args)
    p = Pa * ps(*args)
    back = {h: eps * L, mu: rho * U * L / Re}
    back2 = {Pa: sp.Symbol("mu_") * U * L / (Lam * (eps * L) ** 2)}

    def ex(c):
        c = sp.simplify(c.subs(Pa, mu * U * L / (Lam * h ** 2)).subs(back))
        return sp.simplify(c)

    cont = [_coeff(sp.diff(u, x), U / L), _coeff(sp.diff(v, y), U / L)]
    sx = mu * U / h ** 2 / rho  # divide x-momentum (per unit mass) by νU/h²
    xt = dict(unsteady=sp.diff(u, t), advect_u=u * sp.diff(u, x), advect_v=v * sp.diff(u, y),
              pressure=sp.diff(p, x) / rho, diff_xx=mu / rho * sp.diff(u, x, 2), diff_yy=mu / rho * sp.diff(u, y, 2))
    sy = mu * U * L / (rho * h ** 3)  # makes the ∂p/∂y coefficient 1/Λ
    py = sp.diff(p, x) if printed_8_13b else sp.diff(p, y)
    yt = dict(unsteady=sp.diff(v, t), advect_u=u * sp.diff(v, x), advect_v=v * sp.diff(v, y), pressure=py / rho,
              diff_xx=mu / rho * sp.diff(v, x, 2), diff_yy=mu / rho * sp.diff(v, y, 2))
    xc = {k: ex(_coeff(val, sx)) for k, val in xt.items()}
    yc = {k: ex(_coeff(val, sy)) for k, val in yt.items()}
    exp_x = dict(unsteady=eps ** 2 * Re, advect_u=eps ** 2 * Re, advect_v=eps ** 2 * Re, pressure=1 / Lam,
                 diff_xx=eps ** 2, diff_yy=sp.Integer(1))
    exp_y = dict(unsteady=eps ** 4 * Re, advect_u=eps ** 4 * Re, advect_v=eps ** 4 * Re, pressure=1 / Lam,
                 diff_xx=eps ** 4, diff_yy=eps ** 2)
    uf, pf = sp.Function("u")(x, y), sp.Function("p")(x)
    nu = sp.Symbol("nu", positive=True)
    lim = sp.Eq(0, -sp.diff(pf, x) / rho + nu * sp.diff(uf, y, 2))  # (8.17a) with ν restored
    printed = sp.Eq(0, -sp.diff(pf, x) / rho + sp.diff(uf, y, 2))  # (8.17a) as printed
    px, A, B, Uh, U0 = sp.symbols("p_x A B U_h U_0", real=True)
    p18 = px / mu * y ** 2 / 2 + A * y + B  # (8.18)
    cons = sp.solve([p18.subs(y, 0) - U0, p18.subs(y, h) - Uh], [A, B], dict=True)[0]
    p19 = sp.simplify(p18.subs(cons))
    p19_form = -h ** 2 / (2 * mu) * px * (y / h) * (1 - y / h) + Uh * y / h + U0 * (1 - y / h)
    p19_book = -h ** 2 / (2 * mu) * px * (y / h) * (1 - y / h) + Uh * y / h + U0
    grp = lambda c: dict(inertia=c["unsteady"], pressure=c["pressure"], diff_along=c["diff_xx"],  # noqa: E731
                         diff_across=c["diff_yy"])
    if printed_8_13b:
        printed_y = yc
    else:
        printed_y = {k: ex(_coeff(val, sy)) for k, val in
                     dict(yt, pressure=sp.diff(p, x) / rho).items()}
    usf, vsf = us(x, y, t), vs(x, y, t)
    return dict(continuity_star=sp.Eq(sp.diff(usf, x) + sp.diff(vsf, y), 0), x_coeffs=grp(xc), y_coeffs=grp(yc),
                printed_13b_y_coeffs=grp(printed_y), limit_x=lim, limit_y=sp.Eq(0, -sp.diff(sp.Function("p")(x, y), y) / rho),
                eq_8_18=p18,
                symbols=(L, h, U, rho, mu, Pa, eps, Re, Lam), continuity_coefficients=cont, x_coefficients=xc,
                y_coefficients=yc, expected_x=exp_x, expected_y=exp_y,
                matches=dict(x=all(_z(xc[k] - exp_x[k]) for k in exp_x), y=all(_z(yc[k] - exp_y[k]) for k in exp_y),
                             continuity=all(_z(c - 1) for c in cont)),
                limit_8_17a=lim, printed_8_17a=printed, profile_8_18=p18, profile_8_19=p19,
                profile_8_19_matches=_z(p19 - p19_form), profile_8_19_book=p19_book,
                book_form_top_wall=sp.simplify(p19_book.subs(y, h)))


def slider_bearing_sympy() -> dict:
    """Example 8.1 (slider bearing) step by step in sympy, with the two printed slips as checked wrong variants.

    Book: §8.3, Example 8.1: moving-CV mass balance ⇒ ∫₀ʰ(u − U)dy = C₁ = −h³p′/(12μ) − Uh/2; dp/dx = −12μC₁/h³ −
    6μU/h²; integration with h = h₀(1 + αx/L); p(0) = p(L) = p_e ⇒ C₁ = −(1 + α)Uh₀/(2 + α), C₂; exact p − p_e; the
    O(α) pressure; W = αμL²U/(2h₀²); exact W(α) (ours) and its series.

    Returns
    -------
    dict of sympy objects (design Part C 4.8 keys): C1_eq, dpdx, antiderivative_3, antiderivative_2 (∫₀ˣ of
    (1 + αx/L)^{−3}, ^{−2}), C1, C2, p_exact, ode_residual_exact (0), ode_residual_book (≠ 0: the printed first-power
    denominator), bc_residuals (0, 0), p_linear, W_linear, W_exact, W_exact_series; also symbols, C1_equation,
    p_general, p_printed, W_series, residuals (all 0: C1_book, C2_book, p_exact, p_exact_ode, p_ends, p_linear,
    W_linear, W_exact), printed_final_satisfies_ode (False), printed_intermediate_satisfies_ends (False: the
    (1 − αx/L) integrands). Cached (fresh dict per call).

    Validation (planned): V2 residuals 0; the two printed forms fail. Label: symbolic.
    """
    return dict(_slider_bearing_sympy())


@functools.lru_cache(maxsize=1)
def _slider_bearing_sympy() -> dict:
    """Uncached body of :func:`slider_bearing_sympy`. Book: §8.3, Example 8.1. Label: symbolic."""
    x, y = sp.symbols("x y", real=True)
    L, h0, mu, U = sp.symbols("L h_0 mu U", positive=True)
    al = sp.Symbol("alpha", positive=True)
    C1, C2, px, pe = sp.symbols("C_1 C_2 p_x p_e", real=True)
    hs = sp.Symbol("h", positive=True)
    u = -hs ** 2 / (2 * mu) * px * (y / hs) * (1 - y / hs) + U * y / hs
    c1eq = sp.simplify(sp.integrate(u - U, (y, 0, hs)))
    dpdx_sol = sp.solve(sp.Eq(C1, c1eq), px)[0]
    hx = h0 * (1 + al * x / L)
    dpdx = dpdx_sol.subs(hs, hx)
    pg = sp.integrate(dpdx, x) + C2
    sol = sp.solve([pg.subs(x, 0) - pe, pg.subs(x, L) - pe], [C1, C2], dict=True)[0]
    pex = sp.simplify(pg.subs(sol))
    s = x / L
    ours = pe + 6 * mu * L * U / h0 ** 2 * al * s * (1 - s) / ((2 + al) * (1 + al * s) ** 2)
    printed = pe + 6 * mu * L * U / h0 ** 2 * al * s * (1 - s) / ((2 + al) * (1 + al * s))
    plin = pe + 3 * al * mu * L * U / h0 ** 2 * s * (1 - s)
    Wlin = sp.simplify(sp.integrate(plin - pe, (x, 0, L)))
    Wex = 6 * mu * U * L ** 2 / (h0 ** 2 * al ** 2) * (sp.log(1 + al) - 2 * al / (2 + al))
    Wint = sp.integrate(ours - pe, (x, 0, L))
    c1b = -(1 + al) / (2 + al) * U * h0
    c2b = pe - 6 * mu * L * U / (h0 ** 2 * al) / (2 + al)
    pint_book = 6 * mu * L / (h0 ** 2 * al) * ((C1 / h0) / (1 - al * s) ** 2 + U / (1 - al * s)) + C2
    pint_book = pint_book.subs({C1: c1b, C2: c2b})
    ode = lambda P: sp.simplify(sp.diff(P, x) - dpdx.subs(C1, c1b))  # noqa: E731
    res = dict(C1_book=sp.simplify(sol[C1] - c1b), C2_book=sp.simplify(sol[C2] - c2b),
               p_exact=sp.simplify(pex - ours), p_exact_ode=ode(ours),
               p_ends=sp.simplify(ours.subs(x, 0) - pe) + sp.simplify(ours.subs(x, L) - pe),
               p_linear=sp.simplify(sp.series(ours, al, 0, 2).removeO() - plin),
               W_linear=sp.simplify(Wlin - al * mu * L ** 2 * U / (2 * h0 ** 2)),
               W_exact=sp.simplify(sp.expand_log(Wint - Wex, force=True)))
    xi = sp.Symbol("xi", positive=True)
    anti3 = sp.simplify(sp.integrate(1 / (1 + al * xi / L) ** 3, (xi, 0, x)))  # ∫₀ˣ dx/(1 + αx/L)³
    anti2 = sp.simplify(sp.integrate(1 / (1 + al * xi / L) ** 2, (xi, 0, x)))  # ∫₀ˣ dx/(1 + αx/L)²
    ode_book = ode(printed)
    W_ser = sp.series(Wex, al, 0, 5)
    return dict(C1_eq=sp.Eq(C1, c1eq), dpdx=dpdx, antiderivative_3=anti3, antiderivative_2=anti2,
                C1=sol[C1], C2=sol[C2], p_exact=ours, ode_residual_exact=res["p_exact_ode"],
                ode_residual_book=ode_book, bc_residuals=(sp.simplify(ours.subs(x, 0) - pe),
                                                          sp.simplify(ours.subs(x, L) - pe)),
                p_linear=plin, W_linear=Wlin, W_exact=Wex, W_exact_series=W_ser,
                symbols=(x, L, h0, mu, U, al, pe), C1_equation=sp.Eq(C1, c1eq), p_general=pg, p_printed=printed,
                W_series=W_ser, residuals=res, printed_final_satisfies_ode=_z(ode_book),
                printed_intermediate_satisfies_ends=bool(_z(pint_book.subs(x, 0) - pe) and _z(pint_book.subs(x, L) - pe)))


def reynolds_equation_sympy() -> dict:
    """The 1-D Reynolds equation of lubrication from (8.19) + continuity, step by step in sympy (design D13; the book
    leaves it to Exercises 8.19–8.20 — ours).

    Book: §8.3, (8.19), (6.2), the gap flux (Example 8.1's C₁ generalised). Steps: q = ∫₀ʰ u dy; continuity integrated
    across the gap with v(0) = 0 gives v(h) = −∫₀ʰ ∂u/∂x dy; Leibniz: ∫₀ʰ ∂u/∂x dy = ∂q/∂x − u(h)∂h/∂x with u(h) = U_h;
    kinematic condition at the (material) upper wall v(h) = ∂h/∂t + U_h∂h/∂x — the U_h∂h/∂x terms cancel, leaving
    ∂h/∂t + ∂q/∂x = 0, i.e. ∂/∂x(h³/(12μ) ∂p/∂x) = ∂h/∂t + ((U₀ + U_h)/2)∂h/∂x (U₀, U_h uniform in x).

    Returns
    -------
    dict: q (sympy), leibniz_residual (0), kinematic_cancellation (0: v(h) − U_h h_x + q_x), reynolds_residual (0: the
    p-form minus h_t + q_x), reynolds_equation (Eq), steady_q_constant (Eq ∂q/∂x = 0). Cached (fresh dict per call).

    Validation (planned): V2 every residual 0. Label: symbolic.
    """
    return dict(_reynolds_equation_sympy())


@functools.lru_cache(maxsize=1)
def _reynolds_equation_sympy() -> dict:
    x, y, t = sp.symbols("x y t", real=True)
    mu = sp.Symbol("mu", positive=True)
    Uh, U0 = sp.symbols("U_h U_0", real=True)
    h = sp.Function("h")(x, t)
    p = sp.Function("p")(x, t)
    px = sp.diff(p, x)
    u = -h ** 2 / (2 * mu) * px * (y / h) * (1 - y / h) + Uh * y / h + U0 * (1 - y / h)  # (8.19), consistent form
    q = sp.simplify(sp.integrate(u, (y, 0, h)))
    ux = sp.diff(u, x)
    int_ux = sp.integrate(ux, (y, 0, h))
    leib = sp.simplify(int_ux - (sp.diff(q, x) - u.subs(y, h) * sp.diff(h, x)))
    v_h = -int_ux  # continuity with v(0) = 0
    kin = sp.simplify(v_h - Uh * sp.diff(h, x) + sp.diff(q, x))  # h_t = v(h) − U_h h_x  ⇒  h_t + q_x = 0
    lhs = sp.diff(h ** 3 / (12 * mu) * px, x)
    rhs = sp.diff(h, t) + (U0 + Uh) / 2 * sp.diff(h, x)
    rey = sp.simplify((sp.diff(h, t) + sp.diff(q, x)) - (rhs - lhs))
    return dict(q=q, leibniz_residual=leib, kinematic_cancellation=kin, reynolds_residual=rey,
                reynolds_equation=sp.Eq(lhs, rhs), steady_q_constant=sp.Eq(sp.Derivative(sp.Symbol("q"), x), 0))


# ======================================================================================================================
# §8.4 similarity engine
# ======================================================================================================================
_SIM_KEYS = ("ode", "brackets", "delta", "n", "m", "bcs", "F", "residual")


def _sim_out(d: dict) -> dict:
    """Ensure every design Part C 4.10 key is present (None where it does not apply)."""
    for k in _SIM_KEYS:
        d.setdefault(k, None)
    return d


def similarity_reduce_sympy(case: str = "stokes1") -> dict:
    """Reduce a §8.4 PDE to its similarity ODE with the ansatz (8.32a,b) and solve the exponent bookkeeping (sympy).

    Book: §8.4 — "stokes1": u = UF(η), η = y/√(νt) in (8.20) ⇒ (8.26), BCs (8.27)–(8.28), solution (8.30);
    "stokes1_delta": Example 8.4, u = UF(y/δ(t)) ⇒ brackets [δ′/δ], [ν/δ²] proportional ⇒ δ = √(2C₁νt);
    "vortex_sheet": Example 8.5, ω = At^{−n}F(y/√(νt)), the jump constraint ⇒ n = ½, F = De^{−η²/4};
    "line_vortex": Example 8.6, u_θ = (Γ/2πr)F(r/√(νt)) ⇒ (1/η − η/2)F′ = F″, F = 1 − e^{−η²/4};
    "spreading": Example 8.7, h = At^{−n}F(x/Dt^m) ⇒ 3n + 2m = 1 and m = n (volume) ⇒ n = m = 1/5.

    Parameters
    ----------
    case : "stokes1" | "stokes1_delta" | "vortex_sheet" | "line_vortex" | "spreading".

    Returns
    -------
    dict of sympy objects, always with the design Part C 4.10 keys ode (Eq in F(η)), brackets, delta, n, m, bcs, F
    (closed form), residual (0) — None where a key does not apply — plus per case reduced, solution, pde_residual
    (the dimensional solution in the PDE, 0), exponent_equations …

    Validation (planned): V2 reduced ODEs equal (8.26) and the examples' ODEs; exponents; closed forms. Label: symbolic.
    """
    eta = sp.Symbol("eta", positive=True)
    y, r, xx = sp.symbols("y r x", positive=True)
    t, nu, U, Gam = sp.symbols("t nu U Gamma", positive=True)
    F = sp.Function("F")
    if case == "stokes1":
        u = U * F(y / sp.sqrt(nu * t))
        pde = sp.diff(u, t) - nu * sp.diff(u, y, 2)
        red = sp.simplify((pde * t / U).subs(y, eta * sp.sqrt(nu * t)).doit())
        ode = sp.Eq(-eta / 2 * F(eta).diff(eta), F(eta).diff(eta, 2))  # (8.26)
        sol = 1 - sp.erf(eta / 2)
        ures = sp.simplify(sp.diff(U * (1 - sp.erf(y / (2 * sp.sqrt(nu * t)))), t)
                           - nu * sp.diff(U * (1 - sp.erf(y / (2 * sp.sqrt(nu * t)))), y, 2))
        return _sim_out(dict(case=case, reduced=red, ode=ode, brackets=None, delta=sp.sqrt(nu * t), m=sp.Rational(1, 2), F=sol,
                    matches_8_26=_z(red + (ode.lhs - ode.rhs)) or _z(red - (ode.lhs - ode.rhs)),
                    dsolve=sp.dsolve(ode), solution=sol, A=-1 / sp.sqrt(sp.pi), B=sp.Integer(1),
                    gaussian_integral=sp.integrate(sp.exp(-sp.Symbol("xi") ** 2 / 4), (sp.Symbol("xi"), 0, sp.oo)),
                    bcs=(sol.subs(eta, 0), sp.limit(sol, eta, sp.oo)),
                    residual=sp.simplify(sp.diff(sol, eta, 2) + eta / 2 * sp.diff(sol, eta)), pde_residual=ures,
                    n=sp.Integer(0)))
    if case == "stokes1_delta":
        d = sp.Function("delta")(t)
        C1 = sp.Symbol("C_1", positive=True)
        u = U * F(y / d)
        pde = sp.diff(u, t) - nu * sp.diff(u, y, 2)
        red = sp.expand(sp.simplify((pde / U).subs(y, eta * d).doit()))
        dsol = sp.dsolve(sp.Eq(d.diff(t) / d, C1 * nu / d ** 2), d)
        dsols = dsol if isinstance(dsol, list) else [dsol]
        pos = [s_.rhs for s_ in dsols]
        Cc = sp.Symbol("C1")
        delta = None
        for e in pos:
            e0 = sp.solve(sp.Eq(e.subs(t, 0), 0), Cc)
            if e0:
                cand = sp.simplify(e.subs(Cc, e0[0]))
                if cand.subs({t: 1, nu: 1, C1: 1}).is_positive:
                    delta = cand
        return _sim_out(dict(case=case, reduced=red, brackets=(-d.diff(t) / d, nu / d ** 2), n=sp.Integer(0),
                    ode=sp.Eq(-sp.Symbol("C_1") * eta * F(eta).diff(eta), F(eta).diff(eta, 2)),
                    bracket_t=-d.diff(t) / d, bracket_nu=nu / d ** 2,
                    delta_ode=sp.Eq(d.diff(t) / d, C1 * nu / d ** 2), delta=delta,
                    delta_matches=_z(delta - sp.sqrt(2 * C1 * nu * t)) if delta is not None else False,
                    C1_recovers_8_25=sp.simplify(sp.sqrt(2 * C1 * nu * t).subs(C1, sp.Rational(1, 2)) - sp.sqrt(nu * t)),
                    residual=sp.simplify(sp.sqrt(2 * C1 * nu * t).subs(C1, sp.Rational(1, 2)) - sp.sqrt(nu * t))))
    if case == "vortex_sheet":
        n, A, D = sp.symbols("n A D", real=True)
        w = A * t ** (-n) * F(y / sp.sqrt(nu * t))
        pde = sp.diff(w, t) - nu * sp.diff(w, y, 2)
        red = sp.simplify((pde * t ** (n + 1) / A).subs(y, eta * sp.sqrt(nu * t)).doit())
        ode = sp.Eq(-n * F(eta) - eta / 2 * F(eta).diff(eta), F(eta).diff(eta, 2))
        # jump constraint: −∫ω dy = −A t^{−n} √(νt) ∫F dη = 2U must be t-independent
        tpow = sp.powsimp(t ** (-n) * sp.sqrt(t), force=True)
        n_sol = sp.solve(sp.Eq(sp.log(tpow).expand(force=True).coeff(sp.log(t)), 0), n)[0]
        Fs = D * sp.exp(-eta ** 2 / 4)
        AD = sp.solve(sp.Eq(-A * D * sp.sqrt(nu) * sp.integrate(sp.exp(-eta ** 2 / 4), (eta, -sp.oo, sp.oo)), 2 * U), A)[0] * D
        wz = sp.simplify(AD * t ** (-n_sol) * sp.exp(-y ** 2 / (4 * nu * t)))
        ub = U * sp.erf(y / (2 * sp.sqrt(nu * t)))
        return _sim_out(dict(case=case, reduced=red, ode=ode, n=n_sol, solution=Fs, F=Fs, delta=sp.sqrt(nu * t),
                    m=sp.Rational(1, 2), brackets=(-n / t, -1 / (2 * t), 1 / t), bcs=("F → 0 as η → ∞", "−∫ω dy = 2U"),
                    residual=sp.simplify((ode.lhs - ode.rhs).subs(n, n_sol).subs(F(eta), Fs).doit()),
                    AD=sp.simplify(AD), omega_z=wz, u=ub, omega_matches=_z(wz + sp.diff(ub, y)),
                    pde_residual=sp.simplify(sp.diff(wz, t) - nu * sp.diff(wz, y, 2))))
    if case == "line_vortex":
        u = Gam / (2 * sp.pi * r) * F(r / sp.sqrt(nu * t))
        pde = sp.diff(u, t) - nu * sp.diff(sp.diff(r * u, r) / r, r)
        red = sp.simplify((pde * 2 * sp.pi * r * t / Gam).subs(r, eta * sp.sqrt(nu * t)).doit())
        ode = sp.Eq((1 / eta - eta / 2) * F(eta).diff(eta), F(eta).diff(eta, 2))
        sol = 1 - sp.exp(-eta ** 2 / 4)
        ut = Gam / (2 * sp.pi * r) * (1 - sp.exp(-r ** 2 / (4 * nu * t)))
        us = Gam / (2 * sp.pi * r) * sp.exp(-r ** 2 / (4 * nu * t))
        opr = lambda q: sp.simplify(sp.diff(q, t) - nu * sp.diff(sp.diff(r * q, r) / r, r))  # noqa: E731
        return _sim_out(dict(case=case, reduced=red, ode=ode, bracket=eta ** 2 / 2, brackets=(eta ** 2 / 2,), solution=sol,
                    F=sol, delta=sp.sqrt(nu * t), n=sp.Integer(1), m=sp.Rational(1, 2),
                    residual=sp.simplify((ode.lhs - ode.rhs).subs(F(eta), sol).doit()),
                    bcs=(sol.subs(eta, 0), sp.limit(sol, eta, sp.oo)), pde_residual=opr(ut), spinup_residual=opr(us),
                    u_theta=ut, u_spinup=us))
    if case == "spreading":
        n, m = sp.symbols("n m", real=True)
        eqs = [sp.Eq(-n - 1, -4 * n - 2 * m), sp.Eq(-n + m, 0)]  # ∂h/∂t vs (ρg/3μ)∂(h³h_x)/∂x; volume
        sol = sp.solve(eqs, [n, m], dict=True)[0]
        A, D, beta = sp.symbols("A D beta", positive=True)
        h = A * t ** (-sol[n]) * F(xx / (D * t ** sol[m]))
        pde = sp.diff(h, t) - beta * sp.diff(h ** 3 * sp.diff(h, xx), xx)
        red = sp.simplify((pde * t ** (sol[n] + 1) / A).subs(xx, eta * D * t ** sol[m]).doit())
        return _sim_out(dict(case=case, exponent_equations=eqs, n=sol[n], m=sol[m], reduced=red,
                    ode=sp.Eq(red, 0), delta=D * t ** sol[m], brackets=("t^(-n-1)", "t^(-4n-2m)", "volume t^(m-n)"),
                    t_free=not red.has(t), residual=sp.simplify(sum(sp.Abs((e.lhs - e.rhs).subs(sol)) for e in eqs))))
    raise ValueError('case must be "stokes1", "stokes1_delta", "vortex_sheet", "line_vortex" or "spreading"')


def similarity_ode_solve(case: str = "stokes1", eta_max: float = 12.0, n: int = 400, tol: float = 1e-10,
                         full: bool = False):
    """Solve a §8.4 similarity ODE numerically with ``scipy.integrate.solve_bvp`` (independent of the closed form).

    Book: §8.4 — "stokes1": (8.26) F″ = −(η/2)F′, F(0) = 1, F(η_max) = 0 ((8.27)–(8.28) with ∞ truncated);
    "line_vortex": Example 8.6 F″ = (1/η − η/2)F′, F(η₀) = 0 at η₀ = 1e-4 (the axis; error O(η₀²)), F(η_max) = 1.

    Parameters
    ----------
    case : "stokes1" | "line_vortex";  eta_max : truncation of ∞;  n : initial mesh points;  tol : solve_bvp tolerance.

    full : return the diagnostic dict instead of the (eta, F) pair.

    Returns
    -------
    (eta, F) : arrays (2001 points on [η₀, η_max]);  with ``full=True`` dict(eta, F, F_exact (1 − erf(η/2) or
    1 − e^{−η²/4}), max_err, eta_max, success).

    Validation (planned): V3 max_err ≤ 1e-7, insensitive to η_max ∈ {10, 14}. Label: converged.
    """
    if case == "stokes1":
        e0, bc = 0.0, (1.0, 0.0)
        rhs = lambda e, Y: np.vstack([Y[1], -e / 2.0 * Y[1]])  # noqa: E731
        exact = lambda e: erfc(e / 2.0)  # noqa: E731
    elif case == "line_vortex":
        e0, bc = 1e-4, (0.0, 1.0)
        rhs = lambda e, Y: np.vstack([Y[1], (1.0 / e - e / 2.0) * Y[1]])  # noqa: E731
        exact = lambda e: -np.expm1(-e ** 2 / 4.0)  # noqa: E731
    else:
        raise ValueError('case must be "stokes1" or "line_vortex"')
    e = np.linspace(e0, float(eta_max), int(n))
    guess = np.vstack([bc[0] + (bc[1] - bc[0]) * (1 - np.exp(-e)), (bc[1] - bc[0]) * np.exp(-e)])
    sol = solve_bvp(rhs, lambda a, b: np.array([a[0] - bc[0], b[0] - bc[1]]), e, guess, tol=tol, max_nodes=200000)
    ee = np.linspace(e0, float(eta_max), 2001)
    Fv = sol.sol(ee)[0]
    Fe = exact(ee)
    if not full:
        return ee, Fv
    return dict(eta=ee, F=Fv, F_exact=Fe, max_err=float(np.max(np.abs(Fv - Fe))), eta_max=float(eta_max),
                success=bool(sol.success))


def similarity_collapse_error(case: str = "stokes1", n: float = 0.0, m: float = 0.5, times=(1.0, 2.0, 4.0, 8.0),
                              npts: int = 201) -> float:
    """Spread of rescaled profiles for trial similarity exponents — zero at the right exponents (curation §9, E6, C10).

    Book: §8.4, the ansatz (8.32a) γ = At^{−n}F(ξ/δ), δ ∝ t^m (and (8.32b) γ = Aξ^{−n}F(ξ/δ) for the line vortex).
    Profiles of the exact solutions (ν = 1, U = 1, β = ρg/3μ = 1, unit area) at the given times are rescaled —
    ``stokes1``: u·tⁿ vs y/tᵐ (right: n = 0, m = ½); ``vortex_sheet``: ω·tⁿ vs y/tᵐ (n = m = ½); ``line_vortex``:
    rⁿu_θ vs r/tᵐ (n = 1, m = ½); ``spreading``: h·tⁿ vs x/tᵐ (n = m = 1/5, Huppert's solution) — sampled at the same
    rescaled coordinates ξ ∈ [0, ξ_max] and compared.

    Parameters
    ----------
    case : "stokes1" | "vortex_sheet" | "line_vortex" | "spreading";  n, m : trial exponents;  times : [–];
    npts : samples in ξ.

    Returns
    -------
    spread : max over ξ of (max − min across times) divided by the largest |rescaled profile| [–].

    Validation (planned): V1 ≈ 0 (1e-12) at the right exponents, > 1e-2 away from them. Label: analytic.
    """
    ts = np.asarray(times, dtype=float)
    t0 = ts[0]
    if case == "stokes1":
        f, nat, lo = (lambda X, T: _F(stokes_first_problem(X, T, 1.0, 1.0))), 6.0 * np.sqrt(t0), 0.0
    elif case == "vortex_sheet":
        f, nat, lo = (lambda X, T: _F(vortex_sheet_diffusion(X, T, 1.0, 1.0)[1])), 6.0 * np.sqrt(t0), 0.0
    elif case == "line_vortex":
        f, nat, lo = (lambda X, T: _F(line_vortex_decay(X, T, 2.0 * np.pi, 1.0))), 6.0 * np.sqrt(t0), 0.01
    elif case == "spreading":
        f = lambda X, T: _F(viscous_current_similarity(X, T, 1.0, 3.0, 1.0, 1.0))  # noqa: E731  β = 1
        nat, lo = 1.2 * float(viscous_current_similarity(0.0, t0, 1.0, 3.0, 1.0, 1.0, True)[1]), 0.0
    else:
        raise ValueError('case must be "stokes1", "vortex_sheet", "line_vortex" or "spreading"')
    xi = np.linspace(lo * nat, nat, int(npts)) / t0 ** float(m)
    P = []
    for T in ts:
        X = xi * T ** float(m)
        prof = f(X, T)
        P.append(X ** float(n) * prof if case == "line_vortex" else T ** float(n) * prof)
    P = np.array(P)
    return float(np.max(np.max(P, axis=0) - np.min(P, axis=0)) / max(np.max(np.abs(P[0])), 1e-300))


def vorticity_content(t, U: float = 1.0, nu: float = 1e-6, return_error: bool = False):
    """∫₀^∞ ω dy of Stokes' first problem by ``quad`` — equals +U for every t > 0 (the page prints −U, analysis §9 R10).

    Book: §8.4 (text after (8.30)): no vorticity is generated after t = 0. Parameters: t [s]; U [m/s]; nu [m²/s].
    Returns the value [m/s] (with ``return_error=True``: (value, quad error estimate)). Validation (planned): V1 = U.
    Label: analytic.
    """
    s = np.sqrt(float(nu) * float(t))
    val, err = quad(lambda y: float(stokes_first_vorticity(y, t, U, nu)), 0.0, np.inf, epsabs=0.0, epsrel=1e-12,
                    limit=200) if s > 0 else (0.0, 0.0)
    return (float(val), float(err)) if return_error else float(val)


def stokes_first_pi_groups():
    """Buckingham Π groups of Stokes' first problem: u/U = f(y/√(νt), y/Ut) (Eq. (8.24)) from (u, U, y, t, ν) with U, t
    repeating (``core.dimensional.pi_groups``): three groups u/U, y/(Ut), ν/(U²t) — (8.24) combines the last two into
    y/√(νt). Book: §8.4, Eq. (8.24). Returns list of dict(name → Fraction). Label: symbolic."""
    return pi_groups({"u": "m/s", "U": "m/s", "y": "m", "t": "s", "nu": "m^2/s"}, solution="u", repeating=("U", "t"))


# ======================================================================================================================
# §8.5 Stokes' second problem in sympy
# ======================================================================================================================
def stokes_second_sympy() -> dict:
    """Stokes' second problem derived in sympy: (8.35) → (8.36) → k = ±(1 + i)√(ω/2ν) → (8.37) → bounded root → (8.38).

    Book: §8.5, Eqs. (8.33)–(8.38) (the text's "substitution of (8.33) into (8.20)" means (8.35), analysis §9 R14).

    Returns
    -------
    dict: ode (8.36), k_roots, k_book (±(1 + i)√(ω/2ν)), roots_match, bounded_root, f (8.37 with B = 0, A = U), u
    (8.38), residual (sum, 0) and residuals (ode, pde, wall, vs_book: 0).

    Validation (planned): V2 residuals 0; wrong variant: the growing root is unbounded. Label: symbolic.
    """
    y, t, w, nu, U = sp.symbols("y t omega nu U", positive=True)
    f = sp.Function("f")
    k = sp.Symbol("k")
    ode = sp.Eq(sp.I * w * f(y), nu * f(y).diff(y, 2))  # (8.36)
    roots = sp.solve(sp.Eq(sp.I * w, nu * k ** 2), k)
    kb = (1 + sp.I) * sp.sqrt(w / (2 * nu))
    match = all(any(_z(sp.expand_complex(r_ - s_ * kb)) for s_ in (1, -1)) for r_ in roots)
    fy = U * sp.exp(-kb * y)
    u = sp.simplify(sp.re(sp.expand_complex(sp.exp(sp.I * w * t) * fy)))
    book = U * sp.exp(-y * sp.sqrt(w / (2 * nu))) * sp.cos(w * t - y * sp.sqrt(w / (2 * nu)))  # (8.38)
    res = dict(ode=sp.simplify(sp.I * w * fy - nu * sp.diff(fy, y, 2)), pde=sp.simplify(sp.diff(book, t) - nu * sp.diff(book, y, 2)),
               wall=sp.simplify(book.subs(y, 0) - U * sp.cos(w * t)), vs_book=sp.simplify(sp.expand_trig(u - book)))
    return dict(ode=ode, k_roots=roots, k_book=kb, roots_match=match, bounded_root=-kb, growing_root=kb, f=fy,
                u=u, u_book=book, residuals=res, residual=sp.simplify(sum(sp.Abs(v_) for v_ in res.values())))


# ======================================================================================================================
# §8.6 scaling, Oseen
# ======================================================================================================================
def low_re_scaling_sympy() -> dict:
    """Why the dynamic pressure scale gives the wrong Re → 0 limit, and the viscous scale gives the Stokes equations.

    Book: §8.6, Eqs. (8.40)–(8.43), wrapping ``core.similarity.nondimensional_ns_coefficients`` (time scale l/U,
    steady, no gravity): dynamic scale ρU² — (advective, pressure, viscous) = (1, 1, 1/Re) (8.40); × Re: (Re, Re, 1)
    (8.41) → pressure lost as Re → 0; viscous scale μU/L — (1, 1/Re, 1/Re); × Re: (Re, 1, 1) (8.42) → (8.43).

    Returns
    -------
    dict: dynamic, dynamic_times_Re, viscous, viscous_times_Re (dicts inertia (= advective)/pressure/viscous in the
    symbol Re),
    limits (Re → 0 of the ×Re sets). Label: symbolic.
    """
    l, U, rho, mu = sp.symbols("l U rho mu", positive=True)
    Re = sp.Symbol("Re", positive=True)
    out = {}
    for key in ("dynamic", "viscous"):
        c = SIM.nondimensional_ns_coefficients(pressure_scale=key, time_scale="advective")
        cc = {k: sp.simplify(c[k].subs(mu, rho * U * l / Re)) for k in ("advective", "pressure", "viscous")}
        cc["inertia"] = cc["advective"]
        out[key] = cc
        out[key + "_times_Re"] = {k: sp.simplify(v * Re) for k, v in cc.items()}
        out[key + "_limit"] = {k: sp.limit(v * Re, Re, 0) for k, v in cc.items()}
    return out


def oseen_linearisation_sympy() -> dict:
    """Oseen's linearisation of the advective term about the uniform stream (sympy).

    Book: §8.6: u = U + u′, v = v′, w = w′ ⇒ u∂u/∂x + v∂u/∂y + w∂u/∂z = U∂u′/∂x + [u′∂u′/∂x + v′∂u′/∂y + w′∂u′/∂z];
    dropping the bracket gives ρU∂u′_i/∂x = −∂p/∂x_i + μ∇²u′_i (the page drops the minus on ∂p/∂x_i, analysis §9 R13).

    Returns dict: advection_split (Eq), dropped (the quadratic bracket), oseen_x (Eq, corrected sign),
    expansion_residual (0), linear_part, quadratic_part, oseen_equation (= oseen_x),
    oseen_equation_printed. Label: symbolic.
    """
    x, y, z = sp.symbols("x y z", real=True)
    U, rho, mu = sp.symbols("U rho mu", positive=True)
    up, vp, wp, p = (sp.Function(n_)(x, y, z) for n_ in ("u'", "v'", "w'", "p"))
    u = U + up
    adv = u * sp.diff(u, x) + vp * sp.diff(u, y) + wp * sp.diff(u, z)
    lin = U * sp.diff(up, x)
    quadp = up * sp.diff(up, x) + vp * sp.diff(up, y) + wp * sp.diff(up, z)
    lap = sp.diff(up, x, 2) + sp.diff(up, y, 2) + sp.diff(up, z, 2)
    return dict(advection_split=sp.Eq(adv, lin + quadp), dropped=quadp,
                oseen_x=sp.Eq(rho * U * sp.diff(up, x), -sp.diff(p, x) + mu * lap),
                expansion_residual=sp.simplify(adv - lin - quadp), linear_part=lin, quadratic_part=quadp,
                oseen_equation=sp.Eq(rho * U * sp.diff(up, x), -sp.diff(p, x) + mu * lap),
                oseen_equation_printed=sp.Eq(rho * U * sp.diff(up, x), sp.diff(p, x) + mu * lap))


def oseen_limit_sympy() -> dict:
    """Re → 0 limit of Oseen's stream function (8.53) equals Stokes' (8.48) (the text cites "(9.68)" for (8.48)).

    Book: §8.6, (8.53) and the remark after it. Returns dict: s (the exponent argument), series (to O(Re)),
    limit_psi, stokes_psi, difference (0), and aliases psi_oseen, psi_stokes, limit, residual; first_order
    (the O(Re) correction term). Label: symbolic.
    """
    r, th, a, U = sp.symbols("r theta a U", positive=True)
    Re = sp.Symbol("Re", positive=True)
    ps = U * a ** 2 * ((r ** 2 / (2 * a ** 2) + a / (4 * r)) * sp.sin(th) ** 2
                       - 3 / Re * (1 + sp.cos(th)) * (1 - sp.exp(-Re / 4 * r / a * (1 - sp.cos(th)))))  # (8.53)
    st = U * r ** 2 * sp.sin(th) ** 2 * (sp.Rational(1, 2) - 3 * a / (4 * r) + a ** 3 / (4 * r ** 3))  # (8.48)
    ser = sp.series(ps, Re, 0, 2).removeO()
    lim = sp.limit(ps, Re, 0)
    diff0 = sp.simplify(sp.expand_trig(lim - st).rewrite(sp.cos))
    return dict(s=Re / 4 * r / a * (1 - sp.cos(th)), series=ser, limit_psi=lim, stokes_psi=st, difference=diff0,
                psi_oseen=ps, psi_stokes=st, limit=lim, residual=diff0, first_order=sp.simplify(ser.coeff(Re, 1) * Re))


# ======================================================================================================================
# §8.6 synthetic Millikan experiment
# ======================================================================================================================
def estimate_elementary_charge(q, e_min: float = 1.2e-19, e_max: float = 2.4e-19, n_grid: int = 24001) -> dict:
    """Smallest common quantum of noisy drop charges (the book's "differences … identify the minimum difference").

    Book: §8.6 (Millikan). Method (ours): grid search over e ∈ [e_min, e_max] (the window excludes e/2 and 2e) of the
    mean squared distance of q/e to the nearest integer, then the least-squares refinement e = Σnq/Σn² with the
    integers n fixed. Parameters: q charges [C]; e_min, e_max [C]; n_grid. Returns dict(e [C], n (integers),
    misfit, e_grid_best). Label: analytic.
    """
    q = _F(q)
    eg = np.linspace(float(e_min), float(e_max), int(n_grid))
    ratio = q[None, :] / eg[:, None]
    mis = np.mean((ratio - np.round(ratio)) ** 2, axis=1)
    e0 = float(eg[int(np.argmin(mis))])
    nn = np.round(q / e0)
    e = float(np.sum(nn * q) / np.sum(nn ** 2))
    return dict(e=e, n=nn.astype(int), misfit=float(np.min(mis)), e_grid_best=e0)


def synthetic_millikan(n_drops: int = 40, seed: int = 0, noise: float = 0.01, E: float = 3.2e5, rho_p: float = 886.0,
                       rho: float = 1.2, mu: float = 1.83e-5, g: float = G0, n_max: int = 6) -> dict:
    """A synthetic Millikan oil-drop experiment: drops with random radius and charge ne, fall and rise speeds from
    Stokes' law with multiplicative noise, charges recovered with :func:`core.creeping.millikan_charge`, and e by
    :func:`estimate_elementary_charge`.

    Book: §8.6, Millikan (1911), Fig. 8.18. Synthetic data (seeded ``np.random.default_rng(seed)``); the true e is the
    exact SI value E_CHARGE. Stokes' law is used both to generate and to invert, so the slip (Cunningham) correction of
    real micron drops in air is deliberately absent.

    Parameters
    ----------
    n_drops, seed, noise (relative 1σ on each speed), E [V/m], rho_p, rho [kg/m³], mu [Pa s], g [m/s²], n_max (largest
    number of elementary charges).

    Returns
    -------
    dict: q [C], e_est [C], n (the recovered integers; = n_est), a [m], n_true, U_fall, U_rise [m/s], rel_error, e_true.

    Validation (planned): V5 e recovered within 1 % (seed 0, 40 drops, 1 % noise). Label: benchmark.
    """
    rng = np.random.default_rng(int(seed))
    a = rng.uniform(0.4e-6, 0.9e-6, int(n_drops))
    n_true = rng.integers(1, int(n_max) + 1, int(n_drops))
    Wf = 4.0 / 3.0 * np.pi * a ** 3 * g * (rho_p - rho)
    Uf = Wf / (6.0 * np.pi * mu * a)
    Ur = (n_true * E_CHARGE * E - Wf) / (6.0 * np.pi * mu * a)
    keep = Ur > 0
    a, n_true, Uf, Ur = a[keep], n_true[keep], Uf[keep], Ur[keep]
    Ufm = Uf * (1.0 + noise * rng.standard_normal(Uf.size))
    Urm = Ur * (1.0 + noise * rng.standard_normal(Ur.size))
    q = _F(CRP.millikan_charge(Ufm, Urm, rho_p, rho, mu, E, g))
    est = estimate_elementary_charge(q)
    return dict(a=a, n_true=n_true, U_fall=Ufm, U_rise=Urm, q=q, e_est=est["e"], n_est=est["n"], n=est["n"],
                rel_error=est["e"] / E_CHARGE - 1.0, e_true=E_CHARGE)
