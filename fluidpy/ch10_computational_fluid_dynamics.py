"""Chapter 10 — Computational fluid dynamics: finite differences (stencils, FTCS/BTCS, consistency, von Neumann stability,
upwind and CFL, cell Péclet number), finite elements (weak form, Galerkin, hat functions, element assembly), the two
difficulties of incompressible flow (convection-dominated wiggles; the pressure constraint) with artificial compressibility,
MacCormack, operator splitting / MAC projection on a staggered grid, the Θ-scheme and mixed (Taylor–Hood) elements; the three
worked examples (lid-driven cavity, block in a channel, cylinder in a channel) and the verification habit.

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e (2012), Ch. 10 (by H. H. Hu), §§10.1–10.6, Eqs. (10.1)–(10.199).  Every
equation was transcribed from the rendered page images (chapters/pages/ch10/p448–p499, printed pp. 421–472).

Where the physics lives (every public name is re-exported, so ``ch10.<name>`` reaches every callable):
* ``core.fd`` (FD) — (10.1)–(10.31), (10.84)–(10.94), one-sided stencils (p. 450, (10.148)–(10.150)), (10.199) (lax_demo).
* ``core.fem1d`` (FEM1) — (10.32)–(10.78).
* ``core.mac`` (MAC) — (10.111)–(10.128) (staggered grid, projection, cavity, Taylor–Green, Poiseuille).
* ``core.maccormack`` (MCK) — (10.95)–(10.110), (10.138)–(10.155) (cavity and block).
* ``core.fem2d`` (FEM2) — (10.134)–(10.137), (10.156)–(10.198) (P2–P1 elements, Newton, cylinder, inf–sup).
* this module — the sympy derivation engines, operator splitting / Θ-scheme on a linear test system, the checkerboard and
  gradient null spaces, artificial compressibility, benchmark readers (Ghia 1982, Hou et al. 1995 — public data in
  ``reference/ch10/``), Strouhal from a force history, the CFL time step, the table of printed slips.

Printed slips (analysis §9) coded corrected, printed forms kept as named options that a test must fail: R1 (10.67)–(10.68)
``FEM1.shape_slopes(printed=True)``; R2 ``FEM1.connectivity(printed=True)``; R3 ``FD.steady_cd_fd(scheme="forward")``;
R5 ``MCK.weakly_compressible_step(printed_step5=True)``; R6 ``FEM2.assemble_newton_system(printed_10_172=True)``;
R11 ``FD.transport_1d_step(scheme="upwind_printed")``; R12 ``MCK.maccormack_dt_asymptotic`` vs ``maccormack_dt``; R4, R7, R8,
R9, R10 are text slips (see :func:`book_slips`).  Book-quoted numbers live only in git-ignored ``tests/book_values_ch10.json``;
every §10.5 geometry, Mach number and safety factor used here is ours (design Part C).
"""
from __future__ import annotations

import numpy as np
import sympy as sp

from .core import fd as FD  # noqa: F401
from .core import fem1d as FEM1  # noqa: F401
from .core import fem2d as FEM2  # noqa: F401
from .core import mac as MAC  # noqa: F401
from .core import maccormack as MCK  # noqa: F401
from .core.fd import *  # noqa: F401,F403
from .core.fem1d import *  # noqa: F401,F403
from .core.fem2d import *  # noqa: F401,F403
from .core.mac import *  # noqa: F401,F403
from .core.maccormack import *  # noqa: F401,F403
from .core.diffusion import crank_nicolson_1d, ftcs_diffusion_1d, stable_time_step  # noqa: F401  (ch01/ch08 recalled)
from .core.navier_stokes import exact_solution, navier_stokes_sym  # noqa: F401  (ch04 recalled)

_Z = lambda e: sp.simplify(sp.expand(e)) == 0  # noqa: E731


# ======================================================================================================================
# §10.3 weak form — sympy (D09, D10)
# ======================================================================================================================
def weak_form_sympy() -> dict:
    """Integration by parts that turns the strong form into the weak form (our D09), checked on polynomials.

    Book: §10.3, Eqs. (10.34)–(10.36): ∫₀ᴸ T_t w + u∫₀ᴸ T_x w − D∫₀ᴸ T_xx w = 0 becomes, by parts on the last integral,
    ∫T_t w + u∫T_x w + D∫T_x w_x = D[T_x w]₀ᴸ = D q w(L) when w(0) = 0 and T_x(L) = q.
    T(x, t) = a₀(t) + a₁(t)x + a₂(t)x² + a₃(t)x³ and w = x(b₁ + b₂x) (w(0) = 0) with symbolic coefficients.
    Returns dict(steps (list of str), by_parts (identity residual, 0), weak_minus_strong (0 after inserting T_x(L) = q),
    residual (0)).  Validation: V2.  Label: symbolic.
    """
    x, t, L, u, D, q = sp.symbols("x t L u D q", positive=True)
    a = [sp.Function(f"a{k}")(t) for k in range(4)]
    b1, b2 = sp.symbols("b1 b2")
    T = sum(a[k] * x ** k for k in range(4))
    w = x * (b1 + b2 * x)
    by_parts = (sp.integrate(D * sp.diff(T, x, 2) * w, (x, 0, L))
                - (D * (sp.diff(T, x) * w).subs(x, L) - D * (sp.diff(T, x) * w).subs(x, 0)
                   - D * sp.integrate(sp.diff(T, x) * sp.diff(w, x), (x, 0, L))))  # ∫D T_xx w = D[T_x w] − D∫T_x w_x
    strong = sp.integrate((sp.diff(T, t) + u * sp.diff(T, x) - D * sp.diff(T, x, 2)) * w, (x, 0, L))  # Eq. (10.34)
    weak = (sp.integrate(sp.diff(T, t) * w + u * sp.diff(T, x) * w + D * sp.diff(T, x) * sp.diff(w, x), (x, 0, L))
            - D * sp.diff(T, x).subs(x, L) * w.subs(x, L))  # Eq. (10.35) with the boundary term
    diff = sp.simplify(weak - strong)
    steps = ["multiply (10.1) by w ∈ V and integrate over (0, L): (10.34)",
             "integrate D∫T_xx w by parts: D[T_x w]_0^L − D∫T_x w_x",
             "w(0) = 0 removes the x = 0 end; T_x(L) = q turns the x = L end into D q w(L): (10.35)–(10.36)"]
    return dict(steps=steps, by_parts=sp.simplify(by_parts), weak_minus_strong=diff, residual=sp.simplify(by_parts) + diff)


def weak_to_strong_sympy() -> dict:
    """Reverse the integration by parts (our D10): the weak residual equals ∫(T_t + uT_x − DT_xx)w + D[T_x(L) − q]w(L).

    Book: §10.3, Eqs. (10.37)–(10.38): "satisfying (10.37) for all w ∈ V requires" the PDE on (0, L) and the natural
    condition T_x(L) = q; the Dirichlet condition is not recovered (it is built into S — essential).
    Checked on polynomial T (cubic in x) and w = x(b₁ + b₂x).  Returns dict(identity (0), natural_condition (sympy Eq),
    residual (0)).  Validation: V2.  Label: symbolic.
    """
    x, t, L, u, D, q = sp.symbols("x t L u D q", positive=True)
    a = [sp.Function(f"a{k}")(t) for k in range(4)]
    b1, b2 = sp.symbols("b1 b2")
    T = sum(a[k] * x ** k for k in range(4))
    w = x * (b1 + b2 * x)
    weak = sp.integrate(sp.diff(T, t) * w + u * sp.diff(T, x) * w + D * sp.diff(T, x) * sp.diff(w, x), (x, 0, L)) \
        - D * q * w.subs(x, L)  # Eq. (10.36) residual
    strong = sp.integrate((sp.diff(T, t) + u * sp.diff(T, x) - D * sp.diff(T, x, 2)) * w, (x, 0, L)) \
        + D * (sp.diff(T, x).subs(x, L) - q) * w.subs(x, L)  # Eq. (10.37)
    ident = sp.simplify(weak - strong)
    return dict(identity=ident, natural_condition=sp.Eq(sp.Derivative(sp.Function("T")(x), x).subs(x, L), q),
                residual=ident)


# ======================================================================================================================
# §10.4 — sympy engines (R09, D18, D19, D22)
# ======================================================================================================================
def compressible_ns_sympy() -> dict:
    """The 2-D compressible continuity and momentum equations in conservation form with μ_v = 0 (10.96)–(10.98) against the
    Ch. 4 Navier–Stokes residual (4.38) with Stokes' hypothesis.

    Book: §10.4, Eqs. (10.96)–(10.99) (recap of Ch. 4 D08–D13).  (ρu)_t + (ρu²)_x + (ρvu)_y + p_x − μ∇²u − (μ/3)∂_x(∇·u)
    equals ρDu/Dt − [−p_x + ∂_i(μ(∂u/∂x_i + ∂u_i/∂x)) − ⅔μ∂_x(∇·u)] + u × (continuity residual).
    Returns dict(difference_x, difference_y (0 each), residual (0)).  Validation: V2.  Label: symbolic.
    """
    X, Y, t, mu, c = sp.symbols("x y t mu c", positive=True)
    rho = sp.Function("rho")(X, Y, t)
    u = sp.Function("u")(X, Y, t)
    v = sp.Function("v")(X, Y, t)
    p = c ** 2 * rho  # Eq. (10.99)
    div = sp.diff(u, X) + sp.diff(v, Y)
    cont = sp.diff(rho, t) + sp.diff(rho * u, X) + sp.diff(rho * v, Y)  # Eq. (10.96)
    lap = lambda f: sp.diff(f, X, 2) + sp.diff(f, Y, 2)  # noqa: E731
    mx = sp.diff(rho * u, t) + sp.diff(rho * u ** 2, X) + sp.diff(rho * v * u, Y) + sp.diff(p, X) - mu * lap(u) \
        - mu / 3 * sp.diff(div, X)  # Eq. (10.97), g = 0
    my = sp.diff(rho * v, t) + sp.diff(rho * u * v, X) + sp.diff(rho * v ** 2, Y) + sp.diff(p, Y) - mu * lap(v) \
        - mu / 3 * sp.diff(div, Y)  # Eq. (10.98)
    ns = navier_stokes_sym(rho, [u, v], p, [X, Y], t, mu, mu_v=0, form="4.38")
    dx = sp.simplify(sp.expand(mx - (ns[0] + u * cont)))
    dy = sp.simplify(sp.expand(my - (ns[1] + v * cont)))
    return dict(difference_x=dx, difference_y=dy, residual=sp.simplify(dx) + sp.simplify(dy))


def maccormack_linear_sympy() -> dict:
    """MacCormack for linear advection T_t + uT_x = 0 is Lax–Wendroff (our D18): the update, the local truncation error and
    the amplification factor.

    Book: §10.4, Eqs. (10.100)–(10.102) (the book states only "second-order accurate in both time and space").
    Predictor T*_i = T_i − C(T_{i+1} − T_i), corrector T^{n+1}_i = ½[T_i + T*_i − C(T*_i − T*_{i−1})], C = uΔt/Δx ⇒
    T^{n+1} = T − (C/2)(T₊ − T₋) + (C²/2)(T₊ − 2T + T₋); the modified-equation leading error is −(uΔx²/6)(1 − C²)T_xxx
    (vanishes at C = 1: an exact shift); |G|² = 1 − 4C²(1 − C²) sin⁴(θ/2) ⇒ stable iff C ≤ 1.
    Returns dict(update, lw_difference (0), local_error (per unit time, leading), G2, G2_factorised, G2_difference (0)).
    Validation: V2.  Label: symbolic.
    """
    C, th = sp.symbols("C theta", real=True)
    Tm2, Tm1, T0, Tp1, Tp2 = sp.symbols("T_m2 T_m1 T_0 T_p1 T_p2")
    star = lambda a, b: a - C * (b - a)  # noqa: E731  (T*_i from T_i, T_{i+1})
    s0, sm1 = star(T0, Tp1), star(Tm1, T0)
    upd = sp.expand(sp.Rational(1, 2) * (T0 + s0 - C * (s0 - sm1)))  # Eqs. (10.101)–(10.102)
    lw = sp.expand(T0 - C / 2 * (Tp1 - Tm1) + C ** 2 / 2 * (Tp1 - 2 * T0 + Tm1))
    # local error: exact solution T(x − u t) → substitute Taylor series in Δx with Δt = CΔx/u
    x, h, u = sp.symbols("x h u", positive=True)
    f = sp.Function("f")
    sh = lambda k: sum((k * h) ** m / sp.factorial(m) * sp.Derivative(f(x), (x, m)) if m else f(x) for m in range(6))  # noqa: E731
    num = upd.subs({Tm1: sh(-1), T0: sh(0), Tp1: sh(1)})
    exact = sum((-C * h) ** m / sp.factorial(m) * sp.Derivative(f(x), (x, m)) if m else f(x) for m in range(6))
    err = sp.expand(sp.simplify((num - exact).doit()))
    lead = sp.factor(sp.expand(err).coeff(h, 3))
    local_error = sp.factor(lead * h ** 3 / (C * h / u))  # per unit time: divide by Δt = C h/u
    G = sp.expand(sp.Rational(1, 2) * (1 + (1 - C * (sp.exp(sp.I * th) - 1)) - C * ((1 - C * (sp.exp(sp.I * th) - 1))
                                                                                     * (1 - sp.exp(-sp.I * th)))))
    G2 = sp.simplify(sp.expand(G * sp.conjugate(G)).rewrite(sp.cos))
    G2f = 1 - 4 * C ** 2 * (1 - C ** 2) * sp.sin(th / 2) ** 4
    return dict(update=upd, lw_difference=sp.simplify(upd - lw), local_error=local_error, G2=G2, G2_factorised=G2f,
                G2_difference=sp.simplify(sp.expand_trig(sp.expand((G2 - G2f).rewrite(sp.exp)))))


def projection_sympy() -> dict:
    """The continuous projection (our D19): taking the divergence of (10.117) with (10.118) gives ∇²p^{n+1} = ∇·u*/Δt, and the
    correction −Δt∇p is curl-free.

    Book: §10.4, Eqs. (10.117)–(10.118) (the book writes only the discrete form (10.124)).  u*, v*, p generic sympy functions
    of (x, y).  Returns dict(poisson (sympy Eq ∇²p = ∇·u*/Δt), divergence_after (0 when the Poisson equation holds),
    curl (curl of ∇p, 0)).  Validation: V2.  Label: symbolic.
    """
    X, Y, dt = sp.symbols("x y Delta_t", positive=True)
    us, vs, p = (sp.Function(n)(X, Y) for n in ("u_s", "v_s", "p"))
    u = us - dt * sp.diff(p, X)  # Eq. (10.117)
    v = vs - dt * sp.diff(p, Y)
    div = sp.expand(sp.diff(u, X) + sp.diff(v, Y))  # Eq. (10.118) demands 0
    lap = sp.diff(p, X, 2) + sp.diff(p, Y, 2)
    poisson = sp.Eq(lap, (sp.diff(us, X) + sp.diff(vs, Y)) / dt)
    div_after = sp.simplify(div.subs(sp.diff(p, X, 2), poisson.rhs - sp.diff(p, Y, 2)))
    curl = sp.simplify(sp.diff(sp.diff(p, Y), X) - sp.diff(sp.diff(p, X), Y))
    return dict(poisson=poisson, divergence_after=div_after, curl=curl)


def _strain(u, v, X, Y):
    return sp.Matrix([[sp.diff(u, X), (sp.diff(u, Y) + sp.diff(v, X)) / 2], [(sp.diff(u, Y) + sp.diff(v, X)) / 2, sp.diff(v, Y)]])


def weak_ns_identity_sympy() -> dict:
    """The viscous and pressure terms of the weak Navier–Stokes form (10.134) (our D22), checked on the unit square.

    Book: §10.4, Eqs. (10.134)–(10.136): with ∇·u = 0 and ũ = 0 on the boundary, ∫(1/Re)∇²u·ũ dΩ = −(2/Re)∫D[u]:D[ũ] dΩ
    and ∫(−∇p)·ũ dΩ = ∫p ∇·ũ dΩ (D[u] = ½[∇u + (∇u)ᵀ]).  Checked with u from a polynomial stream function (divergence-free),
    ũ = x(1 − x)y(1 − y)(c₁ + c₂x, c₃ + c₄y), p polynomial.
    Returns dict(viscous (0), pressure (0), residual (0)).  Validation: V2.  Label: symbolic.
    """
    X, Y, Re = sp.symbols("x y Re", positive=True)
    c = sp.symbols("c1:5")
    k = sp.symbols("k1:5")
    psi = k[0] * X ** 2 * Y ** 3 + k[1] * X ** 3 * Y + k[2] * X ** 2 * Y ** 2 + k[3] * X * Y ** 4
    u, v = sp.diff(psi, Y), -sp.diff(psi, X)
    bub = X * (1 - X) * Y * (1 - Y)
    ut, vt = bub * (c[0] + c[1] * X), bub * (c[2] + c[3] * Y)
    p = k[0] * X * Y + k[1] * X ** 2 + k[3] * Y ** 3
    I = lambda f: sp.integrate(sp.integrate(sp.expand(f), (X, 0, 1)), (Y, 0, 1))  # noqa: E731
    lap = lambda f: sp.diff(f, X, 2) + sp.diff(f, Y, 2)  # noqa: E731
    Du, Dt_ = _strain(u, v, X, Y), _strain(ut, vt, X, Y)
    ddot = sum(Du[i, j] * Dt_[i, j] for i in range(2) for j in range(2))
    visc = sp.simplify(I((lap(u) * ut + lap(v) * vt) / Re) + 2 / Re * I(ddot))
    pres = sp.simplify(I(-(sp.diff(p, X) * ut + sp.diff(p, Y) * vt)) - I(p * (sp.diff(ut, X) + sp.diff(vt, Y))))
    return dict(viscous=visc, pressure=pres, residual=sp.simplify(visc + pres))


def weak_ns_components_sympy() -> dict:
    """The Cartesian split of the weak momentum equation (10.156) into (10.157) and (10.158).

    Book: §10.5, Eqs. (10.156)–(10.159): 2D[u]:D[ũ] = u_xũ_x + ½(u_y + v_x)(ũ_y + ṽ_x) + v_yṽ_y (times 2), so the viscous
    integrand of (10.156) equals (1/Re)[2u_xũ_x + (u_y + v_x)ũ_y] + (1/Re)[(u_y + v_x)ṽ_x + 2v_yṽ_y] — the ũ part is (10.157),
    the ṽ part (10.158); the weak continuity (10.159) carries a minus sign so that the pressure blocks are B and Bᵀ.
    Returns dict(difference (0), book_bracket (0: (10.156)'s bracket = D:D̃), residual (0)).  Validation: V2.  Label: symbolic.
    """
    X, Y = sp.symbols("x y")
    u, v, ut, vt = (sp.Function(n)(X, Y) for n in ("u", "v", "ut", "vt"))
    Du, Dt_ = _strain(u, v, X, Y), _strain(ut, vt, X, Y)
    ddot = sum(Du[i, j] * Dt_[i, j] for i in range(2) for j in range(2))
    ux, uy, vx, vy = (sp.diff(u, X), sp.diff(u, Y), sp.diff(v, X), sp.diff(v, Y))
    utx, uty, vtx, vty = (sp.diff(ut, X), sp.diff(ut, Y), sp.diff(vt, X), sp.diff(vt, Y))
    bracket = ux * utx + sp.Rational(1, 2) * (uy + vx) * (uty + vtx) + vy * vty  # Eq. (10.156) bracket
    split = (2 * ux * utx + (uy + vx) * uty) + ((uy + vx) * vtx + 2 * vy * vty)  # Eqs. (10.157) + (10.158)
    d1 = sp.simplify(sp.expand(2 * ddot - split))
    d2 = sp.simplify(sp.expand(ddot - bracket))
    return dict(difference=d1, book_bracket=d2, residual=sp.simplify(d1 + d2))


# ======================================================================================================================
# §10.4 operator splitting and the Θ-scheme on a linear test system (N60, N61, N74)
# ======================================================================================================================
_A1 = np.array([[1.0, 1.0], [0.0, 1.0]])
_A2 = np.array([[1.0, 0.0], [-1.0, 2.0]])


def split_linear_system(A1=None, A2=None, f1=None, f2=None, phi0=None) -> dict:
    """A 2 × 2 linear test problem dφ/dt + (A₁ + A₂)φ = f₁ + f₂ with non-commuting parts, and its exact solution.

    Book: §10.4, Eqs. (10.111)–(10.112) (the operator to be split).  Default (ours): A₁ = [[1, 1], [0, 1]],
    A₂ = [[1, 0], [−1, 2]] (commutator ≠ 0), f = 0, φ₀ = (1, ½).
    Returns dict(A1, A2, f1, f2, phi0, A, commutator (A₁A₂ − A₂A₁), exact (callable t → φ(t) by ``scipy.linalg.expm``)).
    Label: analytic.
    """
    from scipy.linalg import expm

    A1 = _A1 if A1 is None else np.asarray(A1, float)
    A2 = _A2 if A2 is None else np.asarray(A2, float)
    f1 = np.zeros(2) if f1 is None else np.asarray(f1, float)
    f2 = np.zeros(2) if f2 is None else np.asarray(f2, float)
    phi0 = np.array([1.0, 0.5]) if phi0 is None else np.asarray(phi0, float)
    A = A1 + A2
    f = f1 + f2
    xs = np.linalg.solve(A, f)  # steady state (A invertible for the default pair)

    def exact(t):
        return xs + expm(-A * t) @ (phi0 - xs)

    return dict(A1=A1, A2=A2, f1=f1, f2=f2, phi0=phi0, A=A, commutator=A1 @ A2 - A2 @ A1, exact=exact)


def marchuk_yanenko(A1, A2, f1, f2, phi0, dt: float, nsteps: int) -> np.ndarray:
    """Marchuk–Yanenko fractional steps: (φ^{n+1/2} − φⁿ)/Δt + A₁φ^{n+1/2} = f₁, then (φ^{n+1} − φ^{n+1/2})/Δt + A₂φ^{n+1} = f₂.

    Book: §10.4, Eqs. (10.113)–(10.114) (both substeps implicit; first order: the split adds ½Δt²[A₁, A₂] per step).
    Returns φ after nsteps (ndarray).  Validation: V3 order 1 against :func:`split_linear_system` (expm).  Label: converged.
    """
    A1, A2 = np.asarray(A1, float), np.asarray(A2, float)
    I = np.eye(A1.shape[0])
    phi = np.asarray(phi0, float).copy()
    M1, M2 = I + dt * A1, I + dt * A2
    for _ in range(int(nsteps)):
        half = np.linalg.solve(M1, phi + dt * np.asarray(f1, float))  # Eq. (10.113)
        phi = np.linalg.solve(M2, half + dt * np.asarray(f2, float))  # Eq. (10.114)
    return phi


def theta_scheme_linear(L_stokes, L_conv, y0, dt: float, nsteps: int, theta_split: float = 1 - 1 / np.sqrt(2),
                        alpha_split: float | None = None) -> np.ndarray:
    """Glowinski's Θ-scheme on a linear model y′ + (L_s + L_c)y = 0 (L_s the "Stokes" part split with weights α, β; L_c the
    convection part, implicit only in the middle step).

    Book: §10.4, Eqs. (10.129)–(10.133): step 1 (y^{n+θ} − yⁿ)/(θΔt) + αL_s y^{n+θ} = −βL_s yⁿ − L_c yⁿ; step 2
    (y^{n+1−θ} − y^{n+θ})/((1 − 2θ)Δt) + βL_s y^{n+1−θ} + L_c y^{n+1−θ} = −αL_s y^{n+θ}; step 3 as step 1 from y^{n+1−θ};
    α + β = 1, β = θ/(1 − θ) (default; ``alpha_split`` overrides α).  Second order only for θ = 1 − 1/√2 (the book's value
    is printed to five digits; we compute it).  Returns y after nsteps.  Validation: V3 orders 2 (θ = 1 − 1/√2) and 1
    (θ = ¼); V2 :func:`theta_scheme_amplification_sympy`.  Label: converged.
    """
    Ls, Lc = np.atleast_2d(np.asarray(L_stokes, float)), np.atleast_2d(np.asarray(L_conv, float))
    th = float(theta_split)
    beta = th / (1 - th) if alpha_split is None else 1.0 - alpha_split
    alpha = 1.0 - beta
    I = np.eye(Ls.shape[0])
    y = np.atleast_1d(np.asarray(y0, float)).copy()
    S1 = I / (th * dt) + alpha * Ls
    S2 = I / ((1 - 2 * th) * dt) + beta * Ls + Lc
    for _ in range(int(nsteps)):
        y1 = np.linalg.solve(S1, y / (th * dt) - beta * Ls @ y - Lc @ y)  # Eq. (10.129)
        y2 = np.linalg.solve(S2, y1 / ((1 - 2 * th) * dt) - alpha * Ls @ y1)  # Eq. (10.131)
        y = np.linalg.solve(S1, y2 / (th * dt) - beta * Ls @ y2 - Lc @ y2)  # Eq. (10.132)
    return y


def theta_scheme_amplification_sympy() -> dict:
    """The one-step factor R of the Θ-scheme on y′ + (λ₁ + λ₂)y = 0 and its Taylor error against e^{−(λ₁+λ₂)Δt}.

    Book: §10.4, Eqs. (10.129)–(10.133) and the statement that the scheme is second order for θ = 1 − 1/√2 with α + β = 1,
    β = θ/(1 − θ) (Glowinski 1991).  The Δt² coefficient of R − e^{−λΔt} factors as (…)(2θ² − 4θ + 1)(…), which vanishes for
    every (λ₁, λ₂) only at θ = 1 − 1/√2; there the Δt³ coefficient (λ₁ = 1, λ₂ = 0) is (106 − 75√2)/6.
    Returns dict(R_minus_exp (series to Δt³), z2_coeff (factored), theta_root (1 − 1/√2), z3_at_root).
    Validation: V2.  Label: symbolic.
    """
    th, h, l1, l2 = sp.symbols("theta Delta_t lambda1 lambda2", positive=True)
    beta = th / (1 - th)
    alpha = 1 - beta
    y0 = sp.Integer(1)
    y1 = (y0 / (th * h) - beta * l1 * y0 - l2 * y0) / (1 / (th * h) + alpha * l1)  # Eq. (10.129)
    y2 = (y1 / ((1 - 2 * th) * h) - alpha * l1 * y1) / (1 / ((1 - 2 * th) * h) + beta * l1 + l2)  # Eq. (10.131)
    y3 = (y2 / (th * h) - beta * l1 * y2 - l2 * y2) / (1 / (th * h) + alpha * l1)  # Eq. (10.132)
    err = sp.series(y3 - sp.exp(-(l1 + l2) * h), h, 0, 4).removeO()
    c2 = sp.factor(sp.simplify(err.coeff(h, 2)))
    root = 1 - 1 / sp.sqrt(2)
    c3 = sp.nsimplify(sp.simplify(err.coeff(h, 3).subs(th, root).subs({l1: 1, l2: 0})))
    return dict(R_minus_exp=err, z2_coeff=c2, theta_root=root, z3_at_root=c3,
                z2_at_root=sp.simplify(c2.subs(th, root)))


def splitting_order(method: str = "marchuk_yanenko", dt_list=(0.1, 0.05, 0.025, 0.0125), t_end: float = 1.0,
                    theta_split: float | None = None) -> float:
    """Observed order in Δt of a splitting scheme on the default test system of :func:`split_linear_system`.

    Book: §10.4 — Marchuk–Yanenko (10.113)–(10.114) is first order; the Θ-scheme (10.129)–(10.133) second order for
    θ = 1 − 1/√2 (A₁ taken as the "Stokes" part, A₂ as convection).  ``method``: "marchuk_yanenko" | "theta".
    Returns the least-squares slope of log(error at t_end) vs log(Δt).  Expect 1.0 (MY), 2.0 (Θ, default θ), 1.0 (Θ, θ = ¼).
    Label: converged.
    """
    sys_ = split_linear_system()
    ex = sys_["exact"](t_end)
    errs = []
    for dt in dt_list:
        n = int(round(t_end / dt))
        if method == "marchuk_yanenko":
            y = marchuk_yanenko(sys_["A1"], sys_["A2"], sys_["f1"], sys_["f2"], sys_["phi0"], dt, n)
        elif method == "theta":
            th = 1 - 1 / np.sqrt(2) if theta_split is None else theta_split
            y = theta_scheme_linear(sys_["A1"], sys_["A2"], sys_["phi0"], dt, n, th)
        else:
            raise ValueError("method must be marchuk_yanenko or theta")
        errs.append(float(np.linalg.norm(y - ex)))
    return FD.observed_order(list(dt_list), errs)


# ======================================================================================================================
# §10.4 the checkerboard: collocated vs staggered gradients (N68, N69, D21)
# ======================================================================================================================
def checkerboard(nx: int, ny: int) -> np.ndarray:
    """The zigzag pressure p_{i,j} = (−1)^{i+j} on an ny × nx grid (layout [j, i]).  Book: §10.4 (below (10.125)).
    Label: analytic."""
    J, I = np.meshgrid(np.arange(ny), np.arange(nx), indexing="ij")
    return (-1.0) ** (I + J)


def collocated_gradient(p, dx: float, dy: float, periodic: bool = True):
    """Centred pressure gradient over two cells on a collocated grid.

    Book: §10.4, Eq. (10.125): (∂p/∂x)_{i,j} = (p_{i+1,j} − p_{i−1,j})/(2Δx), (∂p/∂y)_{i,j} = (p_{i,j+1} − p_{i,j−1})/(2Δy) —
    a checkerboard p is invisible to it ("felt like a uniform one").  periodic (default) or NaN at the edges.
    Returns (px, py).  Label: analytic.
    """
    p = np.asarray(p, dtype=float)
    return (FD.fd_derivative(p, dx, "central", 1, axis=1, periodic=periodic),
            FD.fd_derivative(p, dy, "central", 1, axis=0, periodic=periodic))  # Eq. (10.125)


def gradient_null_space(nx: int, ny: int, kind: str = "collocated") -> int:
    """Dimension of the null space of the discrete pressure-gradient operator on a periodic nx × ny grid (SVD).

    Book: §10.4, (10.125) collocated vs (10.126) staggered.  Collocated on even grids: 4 (the constant, the checkerboard and
    the two one-directional zigzags); staggered: 1 (only the constant — the physical gauge freedom).  Label: analytic.
    """
    import scipy.sparse as sps

    def d1(n, h, kind_):
        if kind_ == "collocated":
            D = sps.diags([-np.ones(n - 1), np.ones(n - 1)], [-1, 1], format="lil") / (2 * h)
            D[0, n - 1] = -1 / (2 * h)
            D[n - 1, 0] = 1 / (2 * h)
        else:  # staggered: face i − ½ gets (p_i − p_{i−1})/h
            D = sps.diags([-np.ones(n - 1), np.ones(n)], [-1, 0], format="lil") / h
            D[0, n - 1] = -1 / h
        return D.tocsr()

    Gx = sps.kron(sps.identity(ny), d1(nx, 1.0, kind))
    Gy = sps.kron(d1(ny, 1.0, kind), sps.identity(nx))
    G = sps.vstack([Gx, Gy]).toarray()
    s = np.linalg.svd(G, compute_uv=False)
    return int(np.sum(s < 1e-10 * s.max()) + max(0, nx * ny - s.size))


# ======================================================================================================================
# §10.4 artificial compressibility (N52)
# ======================================================================================================================
def artificial_compressibility_channel(ny: int = 16, Re: float = 10.0, c: float = 10.0, dtau: float | None = None,
                                       n_iter: int = 20000, tol: float = 1e-10, nx: int = 4, G: float = 1.0) -> dict:
    """Pseudo-time march of Chorin's artificial-compressibility system to steady plane Poiseuille flow.

    Book: §10.4, Eq. (10.95) ∂p/∂τ + c²∇·u = 0 (pseudotransient; meaningful only at the steady state) with the momentum
    equation (10.81) (body force G = −dp/dx for the mean pressure gradient), on a staggered channel periodic in x with walls
    at y = 0, 1 (MAC operators, quadratic wall ghosts), starting from a divergent perturbation.  Steady answer:
    u = (Re G/2) y(1 − y), v = 0, ∇·u = 0.
    Parameters: ny, Re, c (artificial sound speed), dtau (default from the acoustic and viscous limits), n_iter, tol (max
    change per pseudo-time unit), nx, G.  Returns dict(y, u, exact, max_err, div_history (max|∇·u| every 100 iterations),
    iterations, converged).  Validation: V1 Poiseuille (exact for the quadratic ghosts), V4 ∇·u → 0.  Label: converged.
    """
    g = MAC.MacGrid(nx, ny, 1.0, 1.0, (True, False), dict(top=0.0, bottom=0.0, left=0.0, right=0.0), "quadratic")
    c_ = MAC.face_coordinates(g)
    u = 0.1 * np.sin(2 * np.pi * c_["xu"]) * np.sin(np.pi * c_["yu"])  # divergent start
    v = np.zeros(g.v_shape)
    v[1:-1, :] = 0.05 * np.cos(2 * np.pi * c_["xv"][1:-1]) * np.sin(np.pi * c_["yv"][1:-1])
    p = np.zeros((ny, nx))
    h = min(g.dx, g.dy)
    umax = max(Re * G / 8.0, 1e-3)
    if dtau is None:
        dtau = 0.4 * min(h / (c * np.sqrt(2.0) + umax), Re * h ** 2 / 4.0)
    divh = []
    conv, k = False, 0
    for k in range(1, int(n_iter) + 1):
        cu, cv = MAC.convective_terms(u, v, g)
        lu, lv = MAC.laplacian_faces(u, v, g)
        gx, gy = MAC.gradient(p, g)
        un = u + dtau * (-cu - gx + lu / Re + G)
        vn = v + dtau * (-cv - gy + lv / Re)
        vn[0, :] = vn[-1, :] = 0.0
        pn = p - dtau * c ** 2 * MAC.divergence(un, vn, g)  # Eq. (10.95)
        change = max(np.max(np.abs(un - u)), np.max(np.abs(vn - v))) / dtau
        u, v, p = un, vn, pn
        if k % 100 == 0:
            divh.append(float(np.max(np.abs(MAC.divergence(u, v, g)))))
        if change < tol:
            conv = True
            break
    ex = 0.5 * Re * G * g.yc * (1 - g.yc)
    return dict(y=g.yc, u=u[:, 0], exact=ex, max_err=float(np.max(np.abs(u - ex[:, None]))), div_history=divh,
                div_final=float(np.max(np.abs(MAC.divergence(u, v, g)))), iterations=k, converged=conv, dtau=dtau)


# ======================================================================================================================
# §10.5 benchmarks (public data, reference/ch10) and force histories
# ======================================================================================================================
def _ref_path(name: str):
    from .core.refdata import load

    return load("ch10", name)


def ghia_centreline(Re: int = 100) -> dict:
    """u on the vertical centreline x = ½ of the lid-driven cavity from Ghia, Ghia & Shin (1982), Table I (public data).

    Book: §10.5 (the cavity benchmark; the book compares with Hou et al. 1995 and cites Ghia et al. 1982).  Read from
    ``reference/ch10/ghia1982_table1.csv`` (cited in SOURCES.md).  Re ∈ {100, 400, 1000, 3200, 5000, 7500, 10000}.
    Returns dict(y (ascending 0 … 1, 17 points), u).  Label: benchmark.
    """
    cols = {100: 1, 400: 2, 1000: 3, 3200: 4, 5000: 5, 7500: 6, 10000: 7}
    d = np.loadtxt(_ref_path("ghia1982_table1.csv"), delimiter=",", comments="#", skiprows=2)
    return dict(y=d[:, 0], u=d[:, cols[int(Re)]])


def _centres():
    import csv

    rows = []
    with open(_ref_path("cavity_vortex_centres.csv"), encoding="utf-8") as fh:
        for r in csv.reader(line for line in fh if not line.startswith("#")):
            rows.append(r)
    return [dict(source=r[0], method=r[1], Re=int(r[2]), x=float(r[3]), y=float(r[4])) for r in rows[1:]]


def ghia_vortex_centre(Re: int = 100) -> dict:
    """Primary-vortex centre of Ghia et al. (1982) (via Hajabdollahi & Premnath, arXiv:1202.6351).  Returns dict(x, y).
    Book: §10.5.  Label: benchmark."""
    r = [c for c in _centres() if c["source"] == "Ghia1982" and c["Re"] == int(Re)][0]
    return dict(x=r["x"], y=r["y"])


def hou_centres() -> dict:
    """Primary-vortex centres of Hou et al. (1995) (lattice Boltzmann, 256²; arXiv:comp-gas/9401003).
    Returns dict(Re → (x, y)).  Book: §10.5 (the book's comparison source).  Label: benchmark."""
    return {c["Re"]: (c["x"], c["y"]) for c in _centres() if c["source"] == "Hou1995"}


def cavity_error_vs_ghia(state: dict, Re: int = 100) -> dict:
    """Deviation of a computed centreline u(y) at x = ½ from Ghia et al. (1982) at their y points (np.interp).

    Book: §10.5 (Fig. 10.8 comparison).  ``state``: a :func:`MAC.cavity` result (staggered, with "g") or a
    :func:`MCK.cavity_maccormack` result (node arrays x, y, u).  Returns dict(max_dev, rms_dev, rel_max (max_dev ÷
    max|u_Ghia| = lid speed)).  Label: benchmark.
    """
    gh = ghia_centreline(Re)
    if "g" in state or "grid" in state:
        cl = MAC.cavity_centreline(state)
        y, u = cl["y"], cl["u"]
    else:
        y = np.asarray(state["y"])
        uu = np.asarray(state["u"])
        nx = uu.shape[1]
        u = uu[:, nx // 2] if nx % 2 else 0.5 * (uu[:, nx // 2 - 1] + uu[:, nx // 2])
    d = np.interp(gh["y"], y, u) - gh["u"]
    return dict(max_dev=float(np.max(np.abs(d))), rms_dev=float(np.sqrt(np.mean(d ** 2))),
                rel_max=float(np.max(np.abs(d)) / np.max(np.abs(gh["u"]))))


def strouhal_from_period(tau_bar: float) -> float:
    """Strouhal number St = f_s d/U = 1/τ̄ from the shedding period τ̄ in units of d/U (cyclic frequency).

    Book: §10.5 (S = nd/U with the cyclic frequency n, as in Ch. 9; Ch. 4 (4.102) uses Ω = 2πn).  Label: analytic.
    """
    return float(1.0 / tau_bar)


def dominant_frequency(t, signal, method: str = "zero_crossings") -> float:
    """Dominant cyclic frequency of a (nearly) periodic signal, e.g. the lift coefficient C_L(t).

    Book: §10.5, Fig. 10.21 (the period of the force history → St).  "zero_crossings": mean spacing of the upward zero
    crossings of the mean-removed signal (linear interpolation between samples); "fft": peak of the periodogram (Hann window)
    refined by a parabola through the three highest bins.  Returns f [cycles per unit of t].
    Validation: V1 on a synthetic sine (both methods to 0.1 %).  Label: analytic.
    """
    t = np.asarray(t, dtype=float)
    s = np.asarray(signal, dtype=float) - np.mean(signal)
    if method == "zero_crossings":
        k = np.where((s[:-1] < 0) & (s[1:] >= 0))[0]
        if k.size < 2:
            raise ValueError("fewer than two upward zero crossings")
        tc = t[k] - s[k] * (t[k + 1] - t[k]) / (s[k + 1] - s[k])
        return float((tc.size - 1) / (tc[-1] - tc[0]))
    if method == "fft":
        dt = float(np.mean(np.diff(t)))
        S = np.abs(np.fft.rfft(s * np.hanning(s.size)))
        f = np.fft.rfftfreq(s.size, dt)
        i = int(np.argmax(S[1:]) + 1)
        if 1 <= i < S.size - 1:
            a, b, c = np.log(S[i - 1] + 1e-300), np.log(S[i] + 1e-300), np.log(S[i + 1] + 1e-300)
            delta = 0.5 * (a - c) / (a - 2 * b + c)
            return float(f[i] + delta * (f[1] - f[0]))
        return float(f[i])
    raise ValueError("method must be zero_crossings or fft")


def cfl_time_step(c: float, dx: float, courant: float = 1.0) -> float:
    """Largest explicit time step Δt = C Δx/c for a signal speed c (the CFL condition (10.30) with |u| → c).

    Book: §10.2, Eq. (10.30) (a disturbance must not cross more than one cell per step) — the climate note: the external
    gravity-wave speed √(gH) (Ch. 7) sets the step of explicit ocean/atmosphere models.  c [m/s], dx [m].  Returns Δt [s].
    Examples: c = √(9.81·4000) = 198.1 m/s, Δx = 100 km → 504.8 s; Δx = 25 km → 126.2 s.  Label: analytic.
    """
    return float(courant * dx / c)  # Eq. (10.30)


def book_slips() -> list:
    """The printed slips of Chapter 10 (analysis §9 R1–R12) with the corrected form and the evaluator that shows them.

    Book: §10.3–10.5.  Returns a list of dict(id, where, printed, correct, evaluator).  Label: analytic.
    """
    return [
        dict(id="R1", where="(10.67)-(10.68)", printed="dN_A/dx = -1/(x_A - x_{A-1}), dN_{A+1}/dx = +1/(x_A - x_{A-1}) on [x_{A-1}, x_A]",
             correct="dN_{A-1}/dx = -1/h, dN_A/dx = +1/h (N_1 <-> N_{A-1}, N_2 <-> N_A)", evaluator="FEM1.shape_slopes(xa, xb, printed=True)"),
        dict(id="R2", where="text after (10.73)", printed="nonzero entries need A = e or e + 1", correct="A, B in {e - 1, e} (10.78)",
             evaluator="FEM1.connectivity(n_el, printed=True)"),
        dict(id="R3", where="text after (10.93)", printed="a forward-difference scheme", correct="T_j - T_{j-1} is a backward (upwind for u > 0) difference",
             evaluator="FD.steady_cd_fd(n, R, scheme='forward') (the true downwind difference: root r = 1/(1 - R_cell), "
                       "monotone for R_cell < 1, wiggles only for R_cell > 1) vs scheme='upwind' (never wiggles)"),
        dict(id="R4", where="no-pressure-BC paragraph, p. 445", printed="p_{0,2} will not appear in equation (10.120)",
             correct="... in equation (10.124)", evaluator="text"),
        dict(id="R5", where="cavity algorithm Step 5", printed="2 rho^{n+1} = (rho^n + rho*) - a1 + [(rho u)*_{i,j} - (rho u)*_{i-1,j}] - ...",
             correct="- a1 [(rho u)*_{i,j} - (rho u)*_{i-1,j}] (compare (10.106))", evaluator="MCK.weakly_compressible_step(..., printed_step5=True)"),
        dict(id="R6", where="(10.172)", printed="second sum: u_{A'} * integral N^u_{A',y} N^p_B", correct="v_{A'} (from dv'/dy in (10.167))",
             evaluator="FEM2.assemble_newton_system(..., printed_10_172=True)"),
        dict(id="R7", where="(10.186)", printed="third expansion v' = sum p_b psi_b", correct="p' = sum p_b psi_b", evaluator="text"),
        dict(id="R8", where="p. 458, Fig. 10.14 caption", printed="the fourth example", correct="section 10.5 has three examples", evaluator="text"),
        dict(id="R9", where="p. 470", printed="confined computed St compared with an unbounded value",
             correct="different configurations: qualitative only (unbounded literature 0.164-0.165 at Re = 100)", evaluator="text"),
        dict(id="R10", where="(10.166)", printed="beta dv*/dt(t_n)", correct="beta dv/dt(t_n) (known data at t_n)", evaluator="text"),
        dict(id="R11", where="(10.29)-(10.30)", printed="T_i - T_{i-1} and u dt/dx <= 1 with no sign of u",
             correct="upwind side from sign(u); |u| dt/dx <= 1", evaluator="FD.transport_1d_step(..., scheme='upwind_printed') with u < 0"),
        dict(id="R12", where="(10.155)", printed="dt <= (sigma/sqrt 2) Ma dx 'with large grid Reynolds numbers'",
             correct="also needs Ma << 1 (|u|/dx + |v|/dy << c sqrt(2)/dx)", evaluator="MCK.maccormack_dt_asymptotic vs MCK.maccormack_dt at Ma = 0.5"),
    ]


def derive_all() -> dict:
    """Run every sympy engine of the chapter and report the residuals that must vanish.  Returns dict name → residual (0).
    Label: symbolic."""
    out = dict(weak_form=weak_form_sympy()["residual"], weak_to_strong=weak_to_strong_sympy()["residual"],
               compressible_ns=compressible_ns_sympy()["residual"], maccormack_lw=maccormack_linear_sympy()["lw_difference"],
               projection=projection_sympy()["curl"] + projection_sympy()["divergence_after"],
               weak_ns=weak_ns_identity_sympy()["residual"], weak_ns_components=weak_ns_components_sympy()["residual"],
               theta_z2_at_root=theta_scheme_amplification_sympy()["z2_at_root"],
               galerkin_interior=FEM1.galerkin_equations_sympy(3)["check_interior"])
    return {k: sp.simplify(v) for k, v in out.items()}
