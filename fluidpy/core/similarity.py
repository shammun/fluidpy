"""Dimensionless forms of the equations and dynamic similarity: reference scales, the groups that appear in front
of each term, and the named numbers of fluid mechanics.

Book: Kundu, Cohen & Dowling 5e, Ch. 4 §4.11, Eqs. (4.99)–(4.119), Fig. 4.21, Example 4.8 (rendered pages
chapters/pages/ch04/p171–p178); Mach number also §4.2 (the M < 0.3 rule).

Scales change inside the chapter — each use states its own (``Scales`` records them; never a global):
* time: t* = Ωt (4.100, imposed frequency) or t* = Ut/l (4.109, steady boundary conditions); oscillating body U = lΩ.
* pressure: (p − p∞)/ρU² (4.100) — alternatives μU/l (viscous) and ρgl (hydrostatic) (Exercise 4.59);
  C_p = (p − p∞)/½ρU² (4.106) carries the conventional ½.
* Fr, M are **square roots** of force ratios; Fr′ = U/√(g′l), Ri = 1/Fr′²; We and Bo are ratios of **forces**
  (not per volume): We = ρU²l/σ, Bo = ρgl²/σ.
g = 9.81 m/s² by default (the book's value).

Book typos handled (analysis §9 item 11): "the [,]-brackets in (4.100)" are in (4.101); "those defined in (4.106),
(4.107) becomes" (before (4.114)) → the variables are (4.109) and the equation is (4.112).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
import sympy as sp

from ._util import as_scalar_if_0d
from .thermo import G_BOOK as G0, GAMMA_AIR, R_AIR

__all__ = ["Scales", "nondimensional_ns_coefficients", "nondimensional_ns_sym", "nondimensional_energy_coefficients",
           "nondimensional_continuity_coefficient", "strouhal_number", "reynolds_number", "froude_number",
           "reduced_gravity", "internal_froude_number", "richardson_number", "gradient_richardson_number",
           "mach_number", "compressibility_parameter", "eckert_number", "prandtl_number", "eucken_prandtl",
           "weber_number", "bond_number", "capillary_number", "rossby_number", "pressure_coefficient",
           "drag_coefficient", "lift_coefficient", "reference_area", "sphere_drag_coefficient",
           "froude_scaled_speed", "model_prototype"]

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d


# ======================================================================================================================
# named numbers
# ======================================================================================================================
def strouhal_number(Omega, l, U):
    """St = Ωl/U (unsteady/advective acceleration), Eq. (4.102). Ω [rad/s or 1/s], l [m], U [m/s]. Book: §4.11. Label: analytic."""
    return _S(_F(Omega) * _F(l) / _F(U))


def reynolds_number(U, l, nu=None, rho=None, mu=None):
    """Re = Ul/ν = ρUl/μ (inertia/viscous force), Eq. (4.103). Give ν [m²/s], or ρ [kg/m³] and μ [Pa s].
    Book: §4.11, Eq. (4.103); Fig. 4.21 uses Re = ρUd/μ. ``reynolds_number(1.0, 0.01, 1e-6)`` = 1.0e4. Scalar-callable.
    Label: analytic."""
    if nu is None:
        if rho is None or mu is None:
            raise ValueError("give nu, or rho and mu")
        nu = _F(mu) / _F(rho)
    return _S(_F(U) * _F(l) / _F(nu))


def froude_number(U, l, g: float = G0):
    """Fr = U/√(gl) (square root of inertia/gravity force), Eq. (4.104). Scalar-callable. Book: §4.11. Label: analytic."""
    return _S(_F(U) / np.sqrt(float(g) * _F(l)))


def reduced_gravity(rho1, rho2, g: float = G0):
    """Reduced gravity g′ = g(ρ₂ − ρ₁)/ρ₁ [m/s²] (after (4.105)); layer 1 lighter, above. Book: §4.11. Label: analytic."""
    return _S(float(g) * (_F(rho2) - _F(rho1)) / _F(rho1))


def internal_froude_number(U, g_prime=None, l=1.0, N=None):
    """Internal Froude number Fr′ = U/√(g′l) (two layers, (4.105)) or U/(Nl) (continuous stratification, text after
    (4.105)). Give g′ [m/s²] or N [1/s]; l [m]. Validation: V1 a linear profile with g′ = N²l gives the same Fr′.
    Book: §4.11. Label: analytic."""
    if g_prime is not None:
        return _S(_F(U) / np.sqrt(_F(g_prime) * _F(l)))
    if N is not None:
        return _S(_F(U) / (_F(N) * _F(l)))
    raise ValueError("give g_prime or N")


def richardson_number(g_prime, l, U):
    """Richardson number Ri = 1/Fr′² = g′l/U² (text after (4.105)). ``richardson_number(9.81e-3, 100.0, 0.1)`` = 98.1.
    Book: §4.11. Label: analytic."""
    return _S(_F(g_prime) * _F(l) / _F(U) ** 2)


def gradient_richardson_number(N2, dU_dz):
    """Gradient Richardson number Ri = N²(z)/(dU/dz)² (text after (4.105); N² from ch01's (1.29), Kundu's convention)
    [-]. Book: §4.11. Label: analytic."""
    return _S(_F(N2) / _F(dU_dz) ** 2)


def mach_number(U, c=None, T=None, gamma: float = GAMMA_AIR, R: float = R_AIR):
    """M = U/c (square root of inertia/compressibility force), Eq. (4.111); c = √(γRT) if T [K] is given (perfect gas,
    ch01). §4.2's rule: flows with M < 0.3 (about 100 m/s in air) are nearly incompressible.
    Validation: V5 c(288.15 K) = 340.29 m/s (USSA-1976); M(100 m/s) ≈ 0.294. Book: §4.2 and §4.11, Eq. (4.111). Label: analytic, benchmark."""
    if c is None:
        if T is None:
            raise ValueError("give c or T")
        c = np.sqrt(gamma * R * _F(T))
    return _S(_F(U) / _F(c))


def compressibility_parameter(U, c):
    """U²/c² = M², the coefficient in (4.110) — the size of isentropic departures from ∇·u = 0 (M = 0.3 → 0.09).
    Book: §4.11, Eq. (4.110). Label: analytic."""
    return _S(_F(U) ** 2 / _F(c) ** 2)


def eckert_number(U, cp, dT):
    """Ec = U²/(C_p(T_w − T_o)) (kinetic/thermal energy), Eq. (4.115). U [m/s], C_p [J/(kg K)], ΔT [K]. Book: §4.11. Label: analytic."""
    return _S(_F(U) ** 2 / (_F(cp) * _F(dT)))


def prandtl_number(nu, kappa):
    """Pr = ν/κ (momentum/thermal diffusivity) = μC_p/k, Eq. (4.116) — a fluid property. ν, κ [m²/s]. Book: §4.11. Label: analytic."""
    return _S(_F(nu) / _F(kappa))


def eucken_prandtl(gamma):
    """Eucken's kinetic-theory estimate Pr = 4γ/(9γ − 5): 2/3 for a monatomic gas (the book's 0.67 for hard spheres),
    0.74 for γ = 1.4. Book: §4.11 (the 0.67 quoted after (4.116)); formula: Wikipedia "Prandtl number" (V1 form).
    Label: analytic."""
    g = _F(gamma)
    return _S(4.0 * g / (9.0 * g - 5.0))


def weber_number(rho, U, l, sigma):
    """We = ρU²l/σ (inertia/surface-tension *force*), Eq. (4.117). Book: §4.11. Label: analytic."""
    return _S(_F(rho) * _F(U) ** 2 * _F(l) / _F(sigma))


def bond_number(rho, g, l, sigma):
    """Bo = ρgl²/σ (gravity/surface-tension force), Eq. (4.118) — = (l/ℓ_c)² with ℓ_c the capillary length.
    (For two fluids use the density difference as rho.) Book: §4.11. Label: analytic."""
    return _S(_F(rho) * _F(g) * _F(l) ** 2 / _F(sigma))


def capillary_number(mu, U, sigma):
    """Ca = μU/σ (viscous/surface-tension stress), Eq. (4.119) (= We/Re). Book: §4.11. Label: analytic."""
    return _S(_F(mu) * _F(U) / _F(sigma))


def rossby_number(U, Omega, l, factor: float = 2.0):
    """Rossby number Ro = U/(2Ωl): the advective acceleration U²/l divided by the Coriolis term 2ΩU of (4.45)
    (``factor`` multiplies Ωl; the default 2 gives that ratio). A forward pointer for E9 — not defined in Ch. 4; Ch. 13
    writes Ro = U/(fL) with f = 2Ω sin φ, which equals U/(2ΩL) at the pole. U [m/s], Ω [rad/s], l [m].
    ``rossby_number(10.0, 7.292115e-5, 1.0e6)`` = 0.06857. Book: §4.7, Eq. (4.45) and §4.11. Label: analytic."""
    return _S(_F(U) / (float(factor) * _F(Omega) * _F(l)))


def pressure_coefficient(p, p_inf, rho, U):
    """C_p = (p − p∞)/(½ρU²), Eq. (4.106) (the "Euler number" of some texts, without the ½). Book: §4.11. Label: analytic."""
    return _S((_F(p) - _F(p_inf)) / (0.5 * _F(rho) * _F(U) ** 2))


def drag_coefficient(F, rho, U, A):
    """C_D = F_D/(½ρU²A), Eq. (4.107). F [N], A reference area [m²] (:func:`reference_area`). Book: §4.11. Label: analytic."""
    return _S(_F(F) / (0.5 * _F(rho) * _F(U) ** 2 * _F(A)))


def lift_coefficient(F, rho, U, A):
    """C_L = F_L/(½ρU²A), Eq. (4.108). Book: §4.11. Label: analytic."""
    return _S(_F(F) / (0.5 * _F(rho) * _F(U) ** 2 * _F(A)))


def reference_area(shape: str, **dims) -> float:
    """Reference area A of (4.107)–(4.108) [m²]: "sphere" (d) → πd²/4; "cylinder" (d, b; axis ⟂ flow) → bd;
    "plate"/"airfoil" planform: span × chord (b, c or s, l). Book: §4.11, text after (4.108). Label: analytic."""
    if shape == "sphere":
        return float(np.pi * dims["d"] ** 2 / 4.0)
    if shape == "cylinder":
        return float(dims["b"] * dims["d"])
    if shape in ("plate", "airfoil"):
        span = dims.get("b", dims.get("s"))
        chord = dims.get("c", dims.get("l"))
        return float(span * chord)
    raise ValueError("shape must be 'sphere', 'cylinder', 'plate' or 'airfoil'")


def sphere_drag_coefficient(Re, model: str = "morrison"):
    """Drag coefficient of a smooth sphere vs Re (the curve of Fig. 4.21, drawn with a published correlation).

    Book: §4.11, Fig. 4.21 (C_D ∝ 1/Re at low Re, ≈ constant for 10³–10⁵, the dip between 10⁵ and 10⁶).
    model "morrison": C_D = 24/Re + 2.6(Re/5)/(1 + (Re/5)^1.52) + 0.411(Re/2.63e5)^−7.94/(1 + (Re/2.63e5)^−8.00)
    + 0.25(Re/1e6)/(1 + Re/1e6), valid up to Re = 1e6 (F. A. Morrison, "Data Correlation for Drag Coefficient for
    Sphere", Michigan Tech 2016; *An Introduction to Fluid Mechanics*, CUP 2013); model "stokes": 24/Re (creeping flow).
    Scalar-callable (``sphere_drag_coefficient(1e4)`` = 0.3926). Validation: V5 form reproduced; Re → 0 ratio to 24/Re → 1
    (0.1 % at Re = 0.1). Label: benchmark (form), analytic (Stokes limit).
    """
    R = _F(Re)
    if model == "stokes":
        return _S(24.0 / R)
    if model != "morrison":
        raise ValueError("model must be 'morrison' or 'stokes'")
    x = R / 2.63e5
    cd = (24.0 / R + 2.6 * (R / 5.0) / (1.0 + (R / 5.0) ** 1.52) + 0.411 * x ** -7.94 / (1.0 + x ** -8.00)
          + 0.25 * (R / 1e6) / (1.0 + R / 1e6))
    return _S(cd)


def froude_scaled_speed(U_p, l_p, l_m, g_p: float = G0, g_m: float = G0):
    """Model speed that matches the Froude number: U_m = U_p√(g_m l_m/(g_p l_p)) (Example 4.8;
    ``froude_scaled_speed(10.0, 100.0, 4.0)`` = 2.0 m/s). Book: §4.11. Label: analytic."""
    return _S(_F(U_p) * np.sqrt(float(g_m) * _F(l_m) / (float(g_p) * _F(l_p))))


def _nu(fluid, T: float) -> float:
    if isinstance(fluid, Mapping):
        return float(fluid["nu"])
    from .. import ch01_introduction as ch01  # lazy: ch01's property correlations/tables
    if fluid in ("air", "water"):
        return float(ch01.fluid_properties(fluid, T)["nu"])
    if fluid not in ch01.FLUIDS:
        raise ValueError(f"unknown fluid {fluid!r}; choose from {sorted(ch01.FLUIDS)} or pass dict(nu=...)")
    return float(ch01.FLUIDS[fluid]["nu"])


def model_prototype(l_p: float, U_p: float, scale: float, fluid_p="water", fluid_m="water", match: str = "Fr",
                    T: float = 293.15, g: float = G0) -> dict:
    """Model-test design: the model speed that matches Fr (or Re), and the groups of both flows (E9).

    Book: §4.11, Example 4.8 (Froude scaling of a ship model; Re cannot be matched in the same fluid at the same time).
    Parameters: l_p [m]; U_p [m/s]; scale λ = l_m/l_p; fluids: "air"/"water" (ν(T) from ``ch01.fluid_properties``),
    other ``ch01.FLUIDS`` names (20 °C) or dicts with nu; match "Fr" or "Re"; T [K]; g [m/s²].
    Returns dict(U_m, Re_p, Re_m, Fr_p, Fr_m, Re_ratio = Re_p/Re_m, Fr_ratio, matched (list of matched groups)).
    Validation: V1 Froude match in the same fluid: Re_p/Re_m = λ^(−3/2) (125 for λ = 1/25). Label: analytic.
    """
    nup, num = _nu(fluid_p, T), _nu(fluid_m, T)
    l_m = l_p * scale
    if match == "Fr":
        U_m = float(froude_scaled_speed(U_p, l_p, l_m, g, g))
    elif match == "Re":
        U_m = float(U_p * (num / nup) * (l_p / l_m))
    else:
        raise ValueError("match must be 'Fr' or 'Re'")
    Re_p, Re_m = float(reynolds_number(U_p, l_p, nup)), float(reynolds_number(U_m, l_m, num))
    Fr_p, Fr_m = float(froude_number(U_p, l_p, g)), float(froude_number(U_m, l_m, g))
    matched = [k for k, (a, b) in dict(Re=(Re_p, Re_m), Fr=(Fr_p, Fr_m)).items() if abs(a / b - 1.0) < 1e-9]
    return {"U_m": U_m, "Re_p": Re_p, "Re_m": Re_m, "Fr_p": Fr_p, "Fr_m": Fr_m, "Re_ratio": Re_p / Re_m,
            "Fr_ratio": Fr_p / Fr_m, "matched": matched}


# ======================================================================================================================
# reference scales
# ======================================================================================================================
@dataclass(frozen=True)
class Scales:
    """One set of reference scales for a flow, with the conventions it uses (§4.11, (4.100), (4.109), (4.113)).

    Fields (SI): l [m], U [m/s], rho [kg/m³], mu [Pa s], g [m/s²], Omega [rad/s] (imposed frequency, optional),
    p_inf [Pa], time_scale ("omega" → t_ref = 1/Ω as in (4.100); "advective" → l/U as in (4.109); None → "omega" when
    Omega is given, else "advective"), pressure_scale ("dynamic" ρU², "viscous" μU/l, "hydrostatic" ρgl), and optional
    c [m/s], T_o, T_w [K], cp [J/(kg K)], k [W/(m K)].

    Properties: ``St``, ``Re``, ``Fr``, ``M``, ``Ec``, ``Pr`` (None when the needed scale is missing), ``t_ref``,
    ``p_ref``; methods ``groups()``, ``nondimensionalise(**q)``, ``redimensionalise(**q)``; classmethod
    ``from_oscillation(l, Omega, rho, mu, g)`` (U = lΩ, text after (4.106)).
    Validation: V1 round trip nondimensionalise ↔ redimensionalise (1e-15); St = 1 for from_oscillation. Label: analytic.
    """

    l: float
    U: float
    rho: float
    mu: float
    g: float = G0
    Omega: float | None = None
    p_inf: float = 0.0
    time_scale: str | None = None
    pressure_scale: str = "dynamic"
    c: float | None = None
    T_o: float | None = None
    T_w: float | None = None
    cp: float | None = None
    k: float | None = None

    def __post_init__(self):
        ts = self.time_scale or ("omega" if self.Omega is not None else "advective")
        object.__setattr__(self, "time_scale", ts)
        if ts not in ("omega", "advective"):
            raise ValueError("time_scale must be 'omega' or 'advective'")
        if ts == "omega" and self.Omega is None:
            raise ValueError("time_scale='omega' needs Omega")
        if self.pressure_scale not in ("dynamic", "viscous", "hydrostatic"):
            raise ValueError("pressure_scale must be 'dynamic', 'viscous' or 'hydrostatic'")

    @classmethod
    def from_oscillation(cls, l: float, Omega: float, rho: float, mu: float, g: float = G0, **kw) -> "Scales":
        """Imposed length and frequency only: U = lΩ (text after (4.106)), so St = 1, Re = Ωl²/ν, Fr = Ω√(l/g)."""
        return cls(l=l, U=l * Omega, rho=rho, mu=mu, g=g, Omega=Omega, time_scale="omega", **kw)

    @property
    def t_ref(self) -> float:
        return 1.0 / self.Omega if self.time_scale == "omega" else self.l / self.U

    @property
    def p_ref(self) -> float:
        return {"dynamic": self.rho * self.U ** 2, "viscous": self.mu * self.U / self.l,
                "hydrostatic": self.rho * self.g * self.l}[self.pressure_scale]

    @property
    def St(self):
        return None if self.Omega is None else float(strouhal_number(self.Omega, self.l, self.U))

    @property
    def Re(self):
        return float(reynolds_number(self.U, self.l, rho=self.rho, mu=self.mu))

    @property
    def Fr(self):
        return float(froude_number(self.U, self.l, self.g))

    @property
    def M(self):
        return None if self.c is None else float(mach_number(self.U, c=self.c))

    @property
    def Ec(self):
        if self.cp is None or self.T_o is None or self.T_w is None:
            return None
        return float(eckert_number(self.U, self.cp, self.T_w - self.T_o))

    @property
    def Pr(self):
        if self.cp is None or self.k is None:
            return None
        return float(prandtl_number(self.mu / self.rho, self.k / (self.rho * self.cp)))

    def groups(self) -> dict:
        """The available groups among St, Re, Fr, M, Ec, Pr — (4.102)–(4.104), (4.111), (4.115), (4.116)."""
        return {k: v for k, v in dict(St=self.St, Re=self.Re, Fr=self.Fr, M=self.M, Ec=self.Ec, Pr=self.Pr).items()
                if v is not None}

    def nondimensionalise(self, **q) -> dict:
        """Starred variables of (4.100)/(4.109): x → x/l, t → t/t_ref, u → u/U, p → (p − p∞)/p_ref, g → g/g,
        rho → ρ/ρ, T → (T − T_o)/(T_w − T_o). Only the keys given are returned."""
        conv = {"x": lambda v: v / self.l, "t": lambda v: v / self.t_ref, "u": lambda v: v / self.U,
                "p": lambda v: (v - self.p_inf) / self.p_ref, "g": lambda v: v / self.g, "rho": lambda v: v / self.rho,
                "T": lambda v: (v - self.T_o) / (self.T_w - self.T_o)}
        return {k: conv[k](_F(v)) for k, v in q.items()}

    def redimensionalise(self, **q) -> dict:
        """Inverse of :meth:`nondimensionalise`."""
        conv = {"x": lambda v: v * self.l, "t": lambda v: v * self.t_ref, "u": lambda v: v * self.U,
                "p": lambda v: self.p_inf + v * self.p_ref, "g": lambda v: v * self.g, "rho": lambda v: v * self.rho,
                "T": lambda v: self.T_o + v * (self.T_w - self.T_o)}
        return {k: conv[k](_F(v)) for k, v in q.items()}


# ======================================================================================================================
# sympy: the coefficients of the scaled equations
# ======================================================================================================================
def nondimensional_ns_coefficients(pressure_scale: str = "dynamic", time_scale: str = "omega") -> dict:
    """Coefficients of each term of (4.39b) after substituting the scalings (4.100) and dividing by ρU²/l — by sympy.

    Book: §4.11, Eqs. (4.100)–(4.101): [Ωl/U]∂u*/∂t* + (u*·∇*)u* = −∇*p* + [gl/U²]g* + [μ/ρUl]∇*²u*. Pressure scale
    "dynamic" ρU² (4.100), "viscous" μU/l or "hydrostatic" ρgl (Exercise 4.59); time scale "omega" (1/Ω) or
    "advective" (l/U, (4.109)).

    Method: u = U f(x/l, t/t_ref), p = p∞ + P q(x/l, t/t_ref) (one representative component); each dimensional term is
    divided by ρU²/l and the derivative of f or q is replaced by 1, leaving its coefficient.

    Returns dict(unsteady, advective, pressure, gravity, viscous) of sympy expressions in l, U, rho, mu, g, Omega.
    Validation: V2 exactly Ωl/U, 1, 1, gl/U², μ/(ρUl) for the default scalings; viscous pressure scale → pressure
    coefficient μ/(ρUl) = 1/Re. Label: symbolic.
    """
    l, U, rho, mu, g, Om = sp.symbols("l U rho mu g Omega", positive=True)
    x, t = sp.symbols("x t")
    f, qf, gs = sp.Function("f"), sp.Function("q"), sp.Symbol("g_s")
    tref = 1 / Om if time_scale == "omega" else l / U
    P = {"dynamic": rho * U ** 2, "viscous": mu * U / l, "hydrostatic": rho * g * l}[pressure_scale]
    u = U * f(x / l, t / tref)
    p = P * qf(x / l, t / tref)
    scale = rho * U ** 2 / l
    terms = {"unsteady": rho * sp.diff(u, t), "advective": rho * u * sp.diff(u, x), "pressure": sp.diff(p, x),
             "gravity": rho * g * gs, "viscous": mu * sp.diff(u, x, 2)}
    out = {}
    for k, term in terms.items():
        e = sp.sympify(sp.simplify(term / scale))
        e = sp.sympify(e.replace(lambda a: isinstance(a, (sp.Derivative, sp.Subs)), lambda a: sp.Integer(1)))
        e = sp.sympify(e.replace(lambda a: isinstance(a, sp.core.function.AppliedUndef), lambda a: sp.Integer(1)))
        out[k] = sp.simplify(e.subs(gs, 1))
    return out


def nondimensional_ns_sym(pressure_scale: str = "dynamic", time_scale: str = "omega") -> sp.Eq:
    """The whole scaled momentum equation (4.101) as a sympy Eq with symbolic operators (du*/dt*, (u*·∇*)u*, ∇*p*, g*,
    ∇*²u*) multiplied by the coefficients of :func:`nondimensional_ns_coefficients`. Book: §4.11, Eq. (4.101).
    Label: symbolic."""
    c = nondimensional_ns_coefficients(pressure_scale, time_scale)
    dudt, adv, gradp, gstar, lap = (sp.Symbol(s) for s in (r"\partial_{t^*}u^*", r"(u^*\cdot\nabla^*)u^*",
                                                           r"\nabla^*p^*", r"g^*", r"\nabla^{*2}u^*"))
    return sp.Eq(c["unsteady"] * dudt + c["advective"] * adv, -c["pressure"] * gradp + c["gravity"] * gstar
                 + c["viscous"] * lap)


def nondimensional_energy_coefficients() -> dict:
    """Coefficients of the scaled enthalpy equation (4.114) from (4.112) with (4.109) and (4.113), by sympy.

    Book: §4.11, Eqs. (4.112)–(4.116): ρ*DT*/Dt* = [U²/C_p(T_w − T_o)]Dp*/Dt* + [U²/C_p(T_w − T_o)·μ_o/ρ_oUl]ρ*ε* +
    [k_o/C_pμ_o · μ_o/ρ_oUl]∇*·(k*∇*T*) = Ec Dp*/Dt* + (Ec/Re)ρ*ε* + (1/(Pr Re))∇*·(k*∇*T*).

    Returns dict(pressure_work, dissipation, conduction) of sympy expressions in U, l, rho_o, mu_o, k_o, C_p, T_o, T_w.
    Validation: V2 equals Ec, Ec/Re, 1/(Pr·Re). Label: symbolic.
    """
    U, l, rho, mu, k, cp, To, Tw = sp.symbols("U l rho_o mu_o k_o C_p T_o T_w", positive=True)
    dT = Tw - To
    lhs = rho * cp * dT * U / l  # ρ Dh/Dt with h ≅ C_p T, t ~ l/U
    terms = {"pressure_work": rho * U ** 2 * U / l,  # Dp/Dt with p − p∞ ~ ρ_oU²
             "dissipation": rho * (mu * U ** 2 / (rho * l ** 2)),  # ρε with ε ~ μ_oU²/(ρ_ol²) (4.113)
             "conduction": k * dT / l ** 2}  # ∂(k∂T/∂x)/∂x
    return {kk: sp.simplify(v / lhs) for kk, v in terms.items()}


def nondimensional_continuity_coefficient() -> sp.Expr:
    """Coefficient U²/c² of (4.110): ∇*·u* = −[U²/c²](1/ρ*)Dp*/Dt*, from (4.9) ∇·u = −(1/ρc²)Dp/Dt with (4.109).

    Book: §4.11, Eqs. (4.9) (the §4.11 identity), (4.109)–(4.111). Validation: V2 equals M² = U²/c². Label: symbolic.
    """
    U, l, rho, c = sp.symbols("U l rho_o c", positive=True)
    return sp.simplify(((rho * U ** 2) * (U / l) / (rho * c ** 2)) / (U / l))
