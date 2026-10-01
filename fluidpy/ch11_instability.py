"""Chapter 11 — Instability: normal modes; Kelvin–Helmholtz; thermal (Rayleigh–Bénard) convection; double diffusion; the
Taylor problem; stratified shear flows (Taylor–Goldstein, Richardson, Howard); Squire and Orr–Sommerfeld; inviscid criteria
(Rayleigh, Fjørtoft, critical layers); viscous results (neutral curves, Table 11.1, the disturbance-energy budget);
Tollmien–Schlichting waves; deterministic chaos (pendulum, Hopf, Lorenz, period doubling).

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e (2012), Ch. 11, §§11.1–11.14, Eqs. (11.1)–(11.96).  Every equation was
transcribed from the rendered page images (chapters/pages/ch11/p500–p567, printed pp. 473–540).  Contract:
``analysis/ch11_design.md`` Part C (names, signatures, return types kept).

Where the physics lives (every public name is re-exported, so ``ch11.<name>`` reaches every callable):
* ``core.stability`` (ST) — Chebyshev collocation, constrained generalised eigenproblems, (11.1) vocabulary, the
  Orr–Sommerfeld (11.79), Rayleigh (11.81) and Taylor–Goldstein (11.61) solvers, sweeps, Howard's semicircle (11.72),
  inflection points, Squire (11.78), the energy budget (11.88).
* this module — closed forms (KH (11.18), Ra (11.21), (11.44), (11.46), Ta (11.52), (11.54), Ri (11.66), Lorenz (11.91) …),
  the Bénard/Taylor eigen-solvers, the sympy derivation engines, flow wrappers with cached critical points, Table 11.1, chaos
  tools and the printed slips (:func:`book_slips`).

Conventions (design Part C, analysis §9, curation §9)
------------------------------------------------------
* **Default g = G0 = 9.80665 m/s²** (``core.thermo.G0``); pass ``g=9.81`` explicitly for parity with ``core.waves``.
* **Sign of Γ (slip S10).**  (11.21) uses Γ = −dT̄/dz (> 0 heated from below); Ch. 1's Kundu lapse rate is Γ ≡ dT/dz, the
  meteorological Γ_met = −dT/dz.  No function takes a Γ: :func:`rayleigh_number` takes ``dT`` = T_bottom − T_top
  (:func:`gamma_conventions` gives all three numbers).  §11.5's Ra ≡ gαd⁴(dT̄/dz)/(νκ) is *negative* when heated from below
  (:func:`thermal_rayleigh_signed`).
* **Stream-function signs (slip S9)**: §11.7 u = ∂ψ/∂z, w = −∂ψ/∂x; §11.8 u = ∂ψ/∂y, v = −∂ψ/∂x; §11.14 (and our Bénard
  rolls) u = −∂ψ/∂z, w = ∂ψ/∂x; each function that returns velocities states its convention.
* **Eigenvalues**: σ (e^{σt}) in §11.4, §11.6, §11.14; c (e^{ik(x−ct)}) in §11.3, §11.7–§11.11; σ = −ikc.
* **Non-dimensionalisations**: §11.4 lengths d, time d²/κ (Bénard); §11.6 gap d, x = (R − R₁)/d, σ by d²/ν (Taylor);
  §11.8–§11.11 per flow (half-width and centreline speed for Poiseuille, δ* and U∞ for Blasius, L and U₀ for tanh/sech²);
  §11.14 Lorenz's scaled time and amplitudes (D24).
* Overloaded letters: ``K`` Bénard |K|, ``k`` streamwise/axial wavenumber, ``kappa``/``kappa_s``, ``alpha`` thermal expansion,
  ``mu`` = Ω₂/Ω₁ (the book's α = μ − 1), ``beta_S`` haline, ``sigma`` growth rate (``surface_tension`` for σ_s), ``U_I`` the
  speed at an inflection point, ``r`` Lorenz's Ra/Ra_c, ``Pr`` Lorenz's "σ".

Printed slips (analysis §9): coded corrected; the printed form is a named option a test must fail — S1
``benard_determinant(printed=True)``; S2 ``benard_free_free_mode(printed=True)``; S3 ``benard_free_free_sympy()["printed_root"]``;
S4/S11 ``stratified_shear_sympy(printed=True)``; S6 ``tollmien_profile(printed=True)``; S7 ``taylor_galerkin_Ta(printed=True)``;
S12 ``taylor_perturbation_sympy()["printed_continuity_units_ok"]``; S5, S8, S9, S10 are text/convention slips.
Book-quoted numbers live only in git-ignored ``tests/book_values_ch11.json``; every run parameter here is ours.
"""
from __future__ import annotations

import hashlib
import json
import math
import warnings
from pathlib import Path
from typing import Callable, Sequence

import numpy as np
import scipy.linalg as sla
import sympy as sp
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import BarycentricInterpolator, CubicSpline
from scipy.optimize import brentq, fsolve, minimize_scalar

from .core import stability as ST  # noqa: F401
from .core._util import as_scalar_if_0d
from .core.stability import *  # noqa: F401,F403
from .core.stability import _grid_matching, _parity_mask, _second_derivative
from .core.thermo import G0

RA_FREE_FREE = 27.0 * math.pi ** 4 / 4.0  #: (27/4)π⁴ = 657.511… (§11.4, p. 491) — computed, never typed


def _F(x):
    return np.asarray(x, dtype=float)


def _FC(x):
    """As _F for real input; complex input stays complex (profiles evaluated on the complex path of rayleigh_eigs_contour)."""
    a = np.asarray(x)
    return a if np.iscomplexobj(a) else a.astype(float)


def _S(x):
    return as_scalar_if_0d(x)


# ======================================================================================================================
# caches: outputs/ch11/cache (npz, keyed by a parameter hash) and reference/ch11 (our published tables)
# ======================================================================================================================
def _root() -> Path:
    from .core.project import repo_root

    return repo_root()


def _cache_path(name: str, params: dict) -> Path:
    key = hashlib.sha1(json.dumps(params, sort_keys=True, default=str).encode()).hexdigest()[:12]
    return _root() / "outputs" / "ch11" / "cache" / f"{name}_{key}.npz"


def _ref_critical(key: str) -> dict | None:
    """Our published critical point from reference/ch11/critical_points.json (committed), or None."""
    p = _root() / "reference" / "ch11" / "critical_points.json"
    try:
        d = json.loads(p.read_text(encoding="utf-8"))
        return {k: float(v) for k, v in d[key].items()} if key in d else None
    except (OSError, ValueError, KeyError, TypeError):
        return None


def _cached(name: str, params: dict, compute: Callable[[], dict], cache: bool = True) -> dict:
    """npz cache keyed by (name, parameters).  Invalidation rule: delete ``outputs/ch11/cache`` after changing a solver."""
    path = _cache_path(name, params)
    if cache and path.exists():
        with np.load(path, allow_pickle=False) as z:
            out = {k: z[k] for k in z.files}
        return {k: (v.item() if v.ndim == 0 else v) for k, v in out.items()}
    out = compute()
    if cache:
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            np.savez(path, **{k: np.asarray(v) for k, v in out.items()})
        except OSError:
            pass
    return out


# ======================================================================================================================
# §11.1 mechanical analogies (Fig. 11.1)
# ======================================================================================================================
_WELLS = {
    "bowl": (lambda x: 0.5 * np.asarray(x) ** 2, lambda x: x),
    "cap": (lambda x: -0.5 * np.asarray(x) ** 2, lambda x: -x),
    "plane": (lambda x: 0.0 * np.asarray(x), lambda x: 0.0 * x),
    "dimple": (lambda x: 0.5 * np.asarray(x) ** 2 - 0.25 * np.asarray(x) ** 4, lambda x: x - x ** 3),
}


def potential_well_demo(shape: str, x0: float, v0: float = 0.0, damping: float = 0.3, t_end: float = 20.0,
                        n: int = 400, x_escape: float = 3.0) -> dict:
    """A damped ball in a 1-D potential, ẍ = −V′(x) − γẋ — our model of the four systems of Fig. 11.1.

    Book: §11.1, Fig. 11.1 (stable bowl, unstable cap, neutral plane, a dimple stable to small but not to large disturbances).
    Our potentials (non-dimensional): bowl ½x², cap −½x², plane 0, dimple ½x² − ¼x⁴ (rim ¼ at |x| = 1).
    Parameters: shape ∈ {"bowl", "cap", "plane", "dimple"}; x0, v0 [–]; damping γ ≥ 0 [–]; t_end; n output times; x_escape
    (terminal event at |x| = x_escape — never NaN).
    Returns dict(t, x, v, V (callable), escaped).  Example: dimple x0 = 0.3 stays (|x| < 1), x0 = 1.2 escapes.
    Validation: V4 energy ½v² + V non-increasing for γ ≥ 0 (and the γ = 0 bowl returns to its start after one period);
    no step-size / tolerance convergence study is run.  Label: conserved (an illustration of ours — the potentials and the
    damping are not the book's; solve_ivp rtol 1e-10).
    """
    if shape not in _WELLS:
        raise ValueError(f"shape must be one of {sorted(_WELLS)}")
    V, dV = _WELLS[shape]

    def rhs(t, s):
        return [s[1], -dV(s[0]) - damping * s[1]]

    ev = lambda t, s: abs(s[0]) - x_escape  # noqa: E731
    ev.terminal = True
    sol = solve_ivp(rhs, (0.0, t_end), [x0, v0], t_eval=np.linspace(0.0, t_end, n), events=ev, rtol=1e-10, atol=1e-12)
    return dict(t=sol.t, x=sol.y[0], v=sol.y[1], V=V, escaped=bool(sol.t_events[0].size), shape=shape)


# ======================================================================================================================
# §11.2 normal-mode growth for the explainers (E1)
# ======================================================================================================================
def normal_mode_growth(system: str, k: float, **params) -> complex:
    """Complex growth rate σ (e^{σt} convention) of the growing root at one wavenumber.

    Book: §11.2 (11.1) σ = −i|K|c.  system "interface" (U₁ = U₂ = 0: (11.19); Rayleigh–Taylor when ρ₁ > ρ₂; params rho1,
    rho2, g, surface_tension), "kh" ((11.18) / Ex. 11.1; params U1, U2, rho1, rho2, g, surface_tension, h): σ = −ikc₊ with c₊
    the root with c_i ≥ 0 (for a neutral pair c₊ = mean + √disc); "benard_free" (free–free Bénard, the larger root of D12 with
    K = k; params Ra, Pr, n).
    Returns σ (complex).  Example: ("kh", 1, U1=6, U2=0, rho1=1, rho2=3, g=10) → 1.3229 − 1.5i.  Label: analytic.
    """
    if system in ("interface", "kh"):
        p = dict(U1=0.0, U2=0.0, rho1=1.2, rho2=1000.0, g=G0, surface_tension=0.0, h=None)
        p.update(params)
        if system == "interface":
            p["U1"] = p["U2"] = 0.0
        cp, _ = kh_phase_speed(k, p["U1"], p["U2"], p["rho1"], p["rho2"], p["g"], p["surface_tension"], p["h"])
        return complex(sigma_from_c(k, cp))
    if system == "benard_free":
        return complex(benard_free_free_sigma(k, params.get("Ra", 1000.0), params.get("Pr", 1.0), params.get("n", 1))[0])
    raise ValueError('system must be "interface", "kh" or "benard_free"')


# ======================================================================================================================
# §11.3 Kelvin–Helmholtz
# ======================================================================================================================
def _kh_parts(k, U1, U2, rho1, rho2, g, surface_tension, h):
    k = _F(k)
    if np.any(k <= 0):
        raise ValueError("k > 0 required")
    if rho1 <= 0 or rho2 <= 0:
        raise ValueError("densities must be > 0")
    if surface_tension < 0:
        raise ValueError("surface_tension must be >= 0")
    if h is None:
        C = np.ones_like(k)
    else:
        if h <= 0:
            raise ValueError("h must be > 0 (or None for an infinitely deep lower layer)")
        C = 1.0 / np.tanh(k * h)
    den = rho1 + rho2 * C
    mean = (rho1 * U1 + rho2 * U2 * C) / den
    grav = (g / k) * (rho2 - rho1) / den
    tens = surface_tension * k / den
    shear = -rho1 * rho2 * (U1 - U2) ** 2 * C / den ** 2
    return mean, grav, tens, shear


def kh_phase_speed(k, U1: float, U2: float, rho1: float, rho2: float, g: float = G0, surface_tension: float = 0.0,
                   h: float | None = None):
    """Complex phase speeds (c₊, c₋) of interface waves on a vortex sheet between two streams (Kelvin–Helmholtz).

    Book: §11.3, Eq. (11.18)
        c = (ρ₂U₂ + ρ₁U₁)/(ρ₂ + ρ₁) ± [((ρ₂ − ρ₁)/(ρ₂ + ρ₁))(g/k) − ρ₂ρ₁(U₂ − U₁)²/(ρ₂ + ρ₁)²]^{1/2},
    upper fluid 1 (ρ₁, U₁, z > 0), lower fluid 2.  With ``surface_tension`` σ_s and/or a lower depth ``h`` the general form of
    Exercise 11.1 (C = coth kh; C = 1 for h = None):
        c = (ρ₁U₁ + ρ₂U₂C)/(ρ₁ + ρ₂C) ± [((g/k)(ρ₂ − ρ₁) + σ_s k)/(ρ₁ + ρ₂C) − ρ₁ρ₂(U₁ − U₂)²C/(ρ₁ + ρ₂C)²]^{1/2}.
    Parameters: k [1/m] > 0 (array ok); U1, U2 [m/s]; rho1, rho2 [kg/m³]; g [m/s²]; surface_tension σ_s [N/m]; h [m] or None.
    Returns (c_plus, c_minus) complex [m/s] — principal square root (``np.lib.scimath.sqrt``): c₊ has c_i ≥ 0.
    Example: k = 1, U1 = 6, U2 = 0, ρ = 1, 3, g = 10 → 1.5 ± 1.3229i; U1 = 4 → 2.4142, −0.4142 (neutral).
    Assumptions: inviscid, irrotational perturbations, linear, sheet of zero thickness.
    Validation: V1 both roots satisfy ρ₁(U₁ − c)² + ρ₂(U₂ − c)² = (g/k)(ρ₂ − ρ₁) (p. 479); V2 :func:`kh_sympy`; V7 U₁ = U₂ = 0
    → (11.19) = ``core.waves.interface_omega``/k; ρ₁ = ρ₂ → (11.20); Galilean shift; h → ∞, σ_s = 0 → (11.18).
    """
    mean, grav, tens, shear = _kh_parts(k, U1, U2, rho1, rho2, g, surface_tension, h)
    root = np.lib.scimath.sqrt(grav + tens + shear)  # Eq. (11.18) (Ex. 11.1 form): principal branch, Im ≥ 0
    return _S(mean + root), _S(mean - root)


def kh_discriminant_terms(k, U1: float, U2: float, rho1: float, rho2: float, g: float = G0, surface_tension: float = 0.0,
                          h: float | None = None) -> dict:
    """The three terms under the square root of (11.18) / Ex. 11.1 and the mean speed (the E2 term bars).

    Book: §11.3, (11.18): gravity ((ρ₂ − ρ₁)/(ρ₂ + ρ₁))(g/k), surface tension σ_s k/(ρ₁ + ρ₂), shear −ρ₁ρ₂(U₂ − U₁)²/(ρ₁ + ρ₂)²
    (coth kh factors when h is given, Ex. 11.1).  Returns dict(gravity, tension, shear, total, mean) [m²/s², m/s].
    Example: k = 1, U1 = 6, U2 = 0, ρ = 1, 3, g = 10 → 5, 0, −6.75, −1.75, 1.5.  Label: analytic.
    """
    mean, grav, tens, shear = _kh_parts(k, U1, U2, rho1, rho2, g, surface_tension, h)
    return dict(gravity=_S(grav), tension=_S(tens), shear=_S(shear), total=_S(grav + tens + shear), mean=_S(mean))


def kh_growth_rate(k, U1: float, U2: float, rho1: float, rho2: float, g: float = G0, surface_tension: float = 0.0,
                   h: float | None = None):
    """Temporal growth rate k·max(c_i, 0) [1/s] of the KH mode.  Book: §11.3, (11.18); growth e^{kc_i t} (p. 507).
    Label: analytic.
    """
    cp, _ = kh_phase_speed(k, U1, U2, rho1, rho2, g, surface_tension, h)
    return _S(_F(k) * np.maximum(np.imag(cp), 0.0))


def kh_critical_k(U1: float, U2: float, rho1: float, rho2: float, g: float = G0) -> float:
    """Wavenumber above which the stratified vortex sheet is unstable (no surface tension, deep layers).

    Book: §11.3, p. 480: unstable iff g(ρ₂² − ρ₁²) < kρ₁ρ₂(U₂ − U₁)², so k_c = g(ρ₂² − ρ₁²)/(ρ₁ρ₂(U₂ − U₁)²) [1/m].
    Returns k_c; inf for U₁ = U₂ with ρ₂ ≥ ρ₁ (no shear: no k is unstable); 0.0 when ρ₁ > ρ₂ (Rayleigh–Taylor) or when
    ρ₁ = ρ₂ with U₁ ≠ U₂ (the vortex sheet (11.20)) — every k unstable.
    Degenerate input: a uniform fluid at rest relative to itself (U₁ = U₂ and ρ₁ = ρ₂) has c = U for every k — neutral, no
    interface at all — so k_c = inf (nothing is unstable), not 0.
    Example: air over water, ΔU = 5 m/s → 327 m⁻¹ (λ_c = 1.92 cm).  Validation: V1 c_i = 0 at k_c, > 0 just above.
    Label: analytic.
    """
    num = g * (rho2 ** 2 - rho1 ** 2)
    dU2 = (U2 - U1) ** 2
    if num == 0 and dU2 == 0:
        return math.inf  # uniform fluid, no shear: the discriminant of (11.18) is 0 at every k (neutral)
    if num <= 0:
        return 0.0
    if dU2 == 0:
        return math.inf
    return num / (rho1 * rho2 * dU2)


def kh_stability_boundary(k, rho1: float, rho2: float, g: float = G0, surface_tension: float = 0.0, h: float | None = None):
    """ΔU_min(k): velocity difference at which wavenumber k becomes unstable (zero discriminant of Ex. 11.1).

    Book: Ex. 11.1(b): unstable when (tanh kh + ρ₂/ρ₁)((g/k)(ρ₂ − ρ₁)/ρ₂ + σ_s k/ρ₂) < (U₁ − U₂)².
    Returns ΔU_min [m/s] (array like k).  Label: analytic.
    """
    k = _F(k)
    T = np.ones_like(k) if h is None else np.tanh(k * h)
    val = (T + rho2 / rho1) * ((g / k) * (rho2 - rho1) / rho2 + surface_tension * k / rho2)  # Ex. 11.1(b)
    return _S(np.sqrt(np.maximum(val, 0.0)))


def kh_min_shear(rho1: float, rho2: float, g: float = G0, surface_tension: float = 0.074) -> dict:
    """Minimum velocity difference for KH instability with gravity and surface tension (deep layers).

    Book: §11.3 / Ex. 11.1 (h → ∞): unstable iff (g/k)(ρ₂ − ρ₁) + σ_s k < ρ₁ρ₂ΔU²/(ρ₁ + ρ₂); the left side is least at
    k* = √(g(ρ₂ − ρ₁)/σ_s), so ΔU_min² = 2√(g(ρ₂ − ρ₁)σ_s)(ρ₁ + ρ₂)/(ρ₁ρ₂) (D04).
    Parameters: rho1 < rho2; g; surface_tension > 0 (σ_s = 0 raises — use :func:`kh_critical_k`).
    Returns dict(dU_min [m/s], k_star [1/m], wavelength [m]).  Example: air over water → 6.70 m/s at λ = 1.73 cm.
    Label: analytic.
    """
    if surface_tension <= 0:
        raise ValueError("kh_min_shear needs surface_tension > 0; without it short waves are always unstable — use "
                         "kh_critical_k(U1, U2, rho1, rho2, g) for the cut-off wavenumber")
    drho = rho2 - rho1
    if drho <= 0:
        raise ValueError("kh_min_shear needs rho2 > rho1 (light fluid on top)")
    k_star = math.sqrt(g * drho / surface_tension)
    dU2 = 2.0 * math.sqrt(g * drho * surface_tension) * (rho1 + rho2) / (rho1 * rho2)
    return dict(dU_min=math.sqrt(dU2), k_star=k_star, wavelength=2.0 * math.pi / k_star)


def kh_unstable_band(dU: float, rho1: float, rho2: float, g: float = G0, surface_tension: float = 0.0):
    """The band of unstable wavenumbers (k₁, k₂) for a velocity difference ΔU (deep layers).

    Book: §11.3 / Ex. 11.1: negative discriminant ⇔ σ_s k² − Sk + g(ρ₂ − ρ₁) < 0, S = ρ₁ρ₂ΔU²/(ρ₁ + ρ₂); roots
    k = [S ∓ √(S² − 4σ_s g(ρ₂ − ρ₁))]/(2σ_s).  Without σ_s: (k_c, inf).
    Returns (k1, k2) [1/m]; (nan, nan) when no wavenumber is unstable.  Example: air/water, ΔU = 8 m/s → [149.2, 887.4] m⁻¹;
    ΔU = 5: empty.
    Degenerate inputs: ΔU = 0 with ρ₁ = ρ₂ (a uniform fluid, neutral at every k) → (nan, nan), with or without σ_s;
    ρ₁ > ρ₂ (top-heavy) with σ_s > 0 → (0, k₂): the lower root of the quadratic is negative and is clipped to 0
    (ΔU = 0: k₂ = √(g(ρ₁ − ρ₂)/σ_s) = 2π/:func:`rayleigh_taylor_cutoff`, the Rayleigh–Taylor cut-off).
    Label: analytic.
    """
    drho = rho2 - rho1
    S = rho1 * rho2 * dU ** 2 / (rho1 + rho2)
    if surface_tension <= 0:
        if drho == 0 and S == 0:
            return float("nan"), float("nan")  # uniform fluid, no shear: neutral at every k
        if drho <= 0:
            return 0.0, math.inf
        return (g * drho / S, math.inf) if S > 0 else (float("nan"), float("nan"))
    disc = S ** 2 - 4.0 * surface_tension * g * drho
    if disc <= 0:
        return float("nan"), float("nan")
    r = math.sqrt(disc)
    return max((S - r) / (2.0 * surface_tension), 0.0), (S + r) / (2.0 * surface_tension)


def vortex_sheet_c(U1: float, U2: float):
    """Phase speeds of a vortex sheet between equal densities: c = (U₂ + U₁)/2 ± i(U₂ − U₁)/2.

    Book: §11.3, Eq. (11.20) (ρ₁ = ρ₂, sheet strength γ = U₂ − U₁, §5.8).  Growth kc_i = k|U₂ − U₁|/2 at every k; c_r = mean
    speed (Fig. 11.3).  Returns (c_plus, c_minus) complex with Im c_plus ≥ 0.  Label: analytic.
    """
    m, a = 0.5 * (U2 + U1), 0.5 * abs(U2 - U1)
    return complex(m, a), complex(m, -a)  # Eq. (11.20)


def rayleigh_taylor_cutoff(surface_tension: float, rho_heavy: float, rho_light: float, g: float = G0) -> float:
    """Longest neutrally stable wavelength of heavy fluid over light (Rayleigh–Taylor with surface tension).

    Book: Ex. 11.2 (gravity inverted in c = ±[((g/k)(ρ₂ − ρ₁) + σk)/(ρ₁ + ρ₂ coth kh)]^{1/2}); deep layers: neutral at
    σ_s k² = g(ρ_heavy − ρ_light) ⇒ λ_c = 2π√(σ_s/((ρ_heavy − ρ_light)g)) (ours; the book prints no answer).
    Returns λ_c [m] (shorter waves are stable).  Example: water over air, σ_s = 0.074 → 1.73 cm.  Label: analytic.
    """
    if rho_heavy <= rho_light:
        raise ValueError("rho_heavy must exceed rho_light")
    return 2.0 * math.pi * math.sqrt(surface_tension / ((rho_heavy - rho_light) * g))


def kh_amplitudes(k: float, c: complex, U1: float, U2: float, zeta0: complex = 1.0):
    """Potential amplitudes from the kinematic condition: A₋ = −i(U₁ − c)ζ₀ (upper), A₊ = i(U₂ − c)ζ₀ (lower).

    Book: §11.3, Eq. (11.16) −iU₁kζ₀ − kA₋ = −ikcζ₀ = −iU₂kζ₀ + kA₊ (p. 479).  Returns (A_minus, A_plus) complex.
    Label: analytic.
    """
    return -1j * (U1 - c) * zeta0, 1j * (U2 - c) * zeta0  # Eq. (11.16)


def kh_fields(x, z, t, k: float, U1: float, U2: float, rho1: float, rho2: float, g: float = G0, zeta0: float = 0.01,
              branch: str = "+") -> dict:
    """Real perturbation fields of one KH normal mode: ζ, φ₁ (upper), φ₂ (lower) and their velocities.

    Book: §11.3, Eqs. (11.14)–(11.16): φ₁ = A₋e^{ik(x−ct) − kz}, φ₂ = A₊e^{ik(x−ct) + kz}, ζ = ζ₀e^{ik(x−ct)} (the decaying
    exponentials satisfy (11.3)–(11.5)); c from (11.18) ("+" = c₊, the growing root).  φ₁ is meaningful for z ≥ 0, φ₂ for
    z ≤ 0 (both evaluated everywhere; mask with z).  Velocities u = ∂φ/∂x, w = ∂φ/∂z (perturbation).
    Parameters: x, z [m] (broadcastable); t [s]; k [1/m]; U1, U2, rho1, rho2, g; zeta0 [m]; branch "+"/"−".
    Returns dict(zeta, phi1, phi2, u1, w1, u2, w2, c, phi, u, w (combined by the sign of z)).  Label: analytic.
    """
    cp, cm = kh_phase_speed(k, U1, U2, rho1, rho2, g)
    c = complex(cp if branch == "+" else cm)
    Am, Ap = kh_amplitudes(k, c, U1, U2, zeta0)
    x, z = np.broadcast_arrays(_F(x), _F(z))
    E = np.exp(1j * k * (x - c * t))
    p1 = Am * E * np.exp(-k * z)  # Eq. (11.15)
    p2 = Ap * E * np.exp(k * z)
    up = z >= 0
    out = dict(zeta=np.real(zeta0 * np.exp(1j * k * (_F(x) - c * t))), phi1=np.real(p1), phi2=np.real(p2),
               u1=np.real(1j * k * p1), w1=np.real(-k * p1), u2=np.real(1j * k * p2), w2=np.real(k * p2), c=c)
    out["phi"] = np.where(up, out["phi1"], out["phi2"])
    out["u"] = np.where(up, out["u1"], out["u2"])
    out["w"] = np.where(up, out["w1"], out["w2"])
    return out


def kh_residuals(k: float, U1: float, U2: float, rho1: float, rho2: float, g: float = G0, n_points: int = 20,
                 seed: int = 0) -> dict:
    """Residuals of (11.3), (11.9) and (11.13) for the computed KH modes (both roots) at random (x, z, t) (seeded).

    Book: §11.3.  φ₁, φ₂, ζ from :func:`kh_fields` (ζ₀ = 1); Laplace (11.3) by the exact second derivatives; kinematic (11.9)
    −U₁ζ_x + φ₁_z = ζ_t = −U₂ζ_x + φ₂_z and dynamic (11.13) ρ₁(φ₁_t + U₁φ₁_x + gζ) = ρ₂(φ₂_t + U₂φ₂_x + gζ) at z = 0.
    Returns dict(laplace, kinematic, dynamic) — max |residual| (complex fields), all ≤ 1e-12 relative.  Label: analytic.
    """
    rng = np.random.default_rng(seed)
    out = dict(laplace=0.0, kinematic=0.0, dynamic=0.0)
    for c in kh_phase_speed(k, U1, U2, rho1, rho2, g):
        Am, Ap = kh_amplitudes(k, c, U1, U2, 1.0)
        x, z, t = rng.uniform(0, 2 * math.pi / k, n_points), rng.uniform(-1 / k, 1 / k, n_points), rng.uniform(0, 1, n_points)
        E = np.exp(1j * k * (x - c * t))
        p1, p2 = Am * E * np.exp(-k * z), Ap * E * np.exp(k * z)
        lap = np.abs((-k ** 2) * p1 + (k ** 2) * p1) + np.abs((-k ** 2) * p2 + (k ** 2) * p2)  # φ_xx + φ_zz
        E0 = np.exp(1j * k * (x - c * t))
        zeta, zx, zt = E0, 1j * k * E0, -1j * k * c * E0
        f1z, f2z = -k * Am * E0, k * Ap * E0  # ∂φ/∂z at z = 0
        kin = np.maximum(np.abs(-U1 * zx + f1z - zt), np.abs(-U2 * zx + f2z - zt))  # Eq. (11.9)
        dyn = np.abs(rho1 * (-1j * k * c * Am * E0 + U1 * 1j * k * Am * E0 + g * zeta)
                     - rho2 * (-1j * k * c * Ap * E0 + U2 * 1j * k * Ap * E0 + g * zeta))  # Eq. (11.13)
        scale = max(1.0, abs(c) * k, g, (rho1 + rho2) * g)
        out["laplace"] = max(out["laplace"], float(np.max(lap)) / max(1.0, k ** 2 * abs(Am) + k ** 2 * abs(Ap)))
        out["kinematic"] = max(out["kinematic"], float(np.max(kin)) / max(1.0, abs(c) * k))
        out["dynamic"] = max(out["dynamic"], float(np.max(dyn)) / scale)
    return out


def kh_sympy() -> dict:
    """D02–D03: linearise the interface conditions and derive (11.18) symbolically.

    Book: §11.3, Eqs. (11.6)–(11.18).  (a) The exact kinematic condition (p. 478) −(U₁ + εφ₁_x)εζ_x + εφ₁_z = εζ_t on z = εζ:
    O(ε) at z = 0 gives (11.9).  (b) The dynamic condition with ½|∇φ̃|² and gz on z = εζ: O(ε) gives (11.13).  (c) Normal
    modes (11.14)–(11.16) → A± → the quadratic ρ₁(U₁ − c)² + ρ₂(U₂ − c)² = (g/k)(ρ₂ − ρ₁) → (11.18).  (d) The discriminant
    identity (ρ₁U₁ + ρ₂U₂)² − (ρ₁ + ρ₂)(ρ₁U₁² + ρ₂U₂²) = −ρ₁ρ₂(U₁ − U₂)², and the criterion g(ρ₂² − ρ₁²) < kρ₁ρ₂(U₂ − U₁)².
    Returns dict(kinematic_lin (0), dynamic_lin (0), quadratic, roots, residual_11_18 (0, 0), identity (0), criterion (0),
    vortex_sheet (0, 0), static (0, 0), steps).  Validation: V2.  Label: symbolic.
    """
    x, z, t, eps = sp.symbols("x z t epsilon")
    k, g, r1, r2 = sp.symbols("k g rho1 rho2", positive=True)
    U1, U2, c = sp.symbols("U1 U2 c")
    zeta, f1, f2 = sp.Function("zeta")(x, t), sp.Function("phi1")(x, z, t), sp.Function("phi2")(x, z, t)
    # (a) kinematic, exact (p. 478) with ζ → εζ, φ → εφ; Taylor expand φ_z about z = 0 to O(ε)
    kin_exact = -(U1 + eps * f1.diff(x)) * eps * zeta.diff(x) + eps * f1.diff(z) - eps * zeta.diff(t)
    kin_lin = sp.expand(kin_exact.diff(eps).subs(eps, 0)).subs(z, 0)
    kinematic_lin = sp.simplify(kin_lin - (-U1 * zeta.diff(x) + f1.diff(z).subs(z, 0) - zeta.diff(t)))  # Eq. (11.9)
    # (b) dynamic: ρ₁(φ₁_t + U₁φ₁_x + ½ε|∇φ₁|² + gζ) − ρ₂(…) on z = εζ, at O(ε)
    def bern(rho, U, f):
        return rho * (eps * f.diff(t) + eps * U * f.diff(x) + sp.Rational(1, 2) * eps ** 2 * (f.diff(x) ** 2 + f.diff(z) ** 2)
                      + g * eps * zeta)
    dyn = sp.expand((bern(r1, U1, f1) - bern(r2, U2, f2)).diff(eps).subs(eps, 0)).subs(z, 0)
    dynamic_lin = sp.simplify(dyn - (r1 * (f1.diff(t) + U1 * f1.diff(x) + g * zeta)
                                     - r2 * (f2.diff(t) + U2 * f2.diff(x) + g * zeta)).subs(z, 0))  # Eq. (11.13)
    # (c) normal modes
    z0 = sp.symbols("zeta0", positive=True)
    Am, Ap = sp.symbols("A_m A_p")
    sol = sp.solve([sp.Eq(-sp.I * U1 * k * z0 - k * Am, -sp.I * k * c * z0),
                    sp.Eq(-sp.I * U2 * k * z0 + k * Ap, -sp.I * k * c * z0)], [Am, Ap], dict=True)[0]  # Eq. (11.16)
    eq17 = r1 * (-sp.I * k * c * Am + sp.I * k * U1 * Am + g * z0) - r2 * (-sp.I * k * c * Ap + sp.I * k * U2 * Ap + g * z0)
    quad_ = sp.expand(eq17.subs(sol) / (k * z0))  # (11.17) after A± and ÷ kζ₀
    target = sp.expand(r1 * (U1 - c) ** 2 + r2 * (U2 - c) ** 2 - (g / k) * (r2 - r1))
    prop = sp.simplify(quad_ / target)
    roots = sp.solve(target, c)
    disc = (r2 - r1) / (r2 + r1) * g / k - r2 * r1 / (r2 + r1) ** 2 * (U2 - U1) ** 2
    book = [(r2 * U2 + r1 * U1) / (r2 + r1) + s * sp.sqrt(disc) for s in (1, -1)]  # Eq. (11.18)
    res = [sp.simplify(target.subs(c, b)) for b in book]
    identity = sp.expand((r1 * U1 + r2 * U2) ** 2 - (r1 + r2) * (r1 * U1 ** 2 + r2 * U2 ** 2) + r1 * r2 * (U1 - U2) ** 2)
    criterion = sp.simplify(k * (r1 + r2) ** 2 * disc - (g * (r2 ** 2 - r1 ** 2) - k * r1 * r2 * (U2 - U1) ** 2))
    vs = [sp.simplify(target.subs({r2: r1}).subs(c, (U2 + U1) / 2 + s * sp.I * (U2 - U1) / 2)) for s in (1, -1)]  # (11.20)
    st = [sp.simplify(target.subs({U1: 0, U2: 0}).subs(c, s * sp.sqrt((r2 - r1) / (r2 + r1) * g / k))) for s in (1, -1)]
    return dict(kinematic_lin=kinematic_lin, dynamic_lin=dynamic_lin, quadratic=target, proportionality=prop, roots=roots,
                residual_11_18=res, identity=identity, criterion=criterion, vortex_sheet=vs, static=st,
                steps=["linearise (11.8) and the Bernoulli matching at z = 0 → (11.9), (11.13)",
                       "normal modes (11.14) → decaying exponentials (11.15) → (11.16), (11.17)",
                       "(11.16) → A₋ = −i(U₁ − c)ζ₀, A₊ = i(U₂ − c)ζ₀; into (11.17), divide by kζ₀",
                       "ρ₁(U₁ − c)² + ρ₂(U₂ − c)² = (g/k)(ρ₂ − ρ₁); quadratic formula → (11.18)"])


def kh_depth_tension_sympy() -> dict:
    """D04 (Exercise 11.1): finite lower depth h and surface tension — the dispersion relation and its limits.

    Book: Ex. 11.1: lower mode φ₂ = B cosh k(z + h)e^{ik(x−ct)} (v₂ = 0 at z = −h), upper φ₁ = Ae^{−kz}e^{ik(x−ct)}; kinematic
    conditions as (11.16); pressure jump p₁ − p₂ = σ_s ζ_xx with the linearised Bernoulli pressures.  Result: the Ex. 11.1(a)
    formula (coth kh); h → ∞ gives (11.18) + σ_s k; ΔU_min² = 2√(gΔρσ_s)(ρ₁ + ρ₂)/(ρ₁ρ₂) at k* = √(gΔρ/σ_s).
    Returns dict(dispersion (the quadratic in c), residual (coefficient-wise difference from the quadratic whose roots are the
    Ex. 11.1(a) formula, after scaling: [0, 0, 0]), limit_residual (0), coth_limit (1), dU_min_sq, k_star).
    Validation: V2.  Label: symbolic.
    """
    k, g, r1, r2, sg, h = sp.symbols("k g rho1 rho2 sigma_s h", positive=True)
    U1, U2, c = sp.symbols("U1 U2 c")
    z0 = sp.symbols("zeta0", positive=True)
    A, B = sp.symbols("A B")
    Ch = sp.cosh(k * h) / sp.sinh(k * h)
    sol = sp.solve([sp.Eq(-sp.I * U1 * k * z0 - k * A, -sp.I * k * c * z0),
                    sp.Eq(-sp.I * U2 * k * z0 + k * B * sp.sinh(k * h), -sp.I * k * c * z0)], [A, B], dict=True)[0]
    # p_j′ = −ρ_j(φ_t + U_jφ_x + gζ) at z = 0; p₁ − p₂ = σ_s ζ_xx = −σ_s k² ζ
    p1 = -r1 * (-sp.I * k * c * A + sp.I * k * U1 * A + g * z0)
    p2 = -r2 * (-sp.I * k * c * B * sp.cosh(k * h) + sp.I * k * U2 * B * sp.cosh(k * h) + g * z0)
    dispersion = sp.simplify((p1 - p2 + sg * k ** 2 * z0).subs(sol) / (k * z0))
    Cs = sp.symbols("C", positive=True)  # stands for coth kh
    den = r1 + r2 * Cs
    mean = (r1 * U1 + r2 * U2 * Cs) / den
    disc = ((g / k) * (r2 - r1) + sg * k) / den - r1 * r2 * (U1 - U2) ** 2 * Cs / den ** 2  # Ex. 11.1(a)
    Q = sp.expand((den * (c - mean) ** 2 - den * disc).subs(Cs, Ch))  # (c − mean)² = disc, cleared of fractions
    pd_, pq_ = sp.Poly(sp.expand(dispersion), c), sp.Poly(Q, c)
    ratio = pd_.LC() / pq_.LC()
    res = [sp.simplify((a - ratio * b).rewrite(sp.exp)) for a, b in zip(pd_.all_coeffs(), pq_.all_coeffs())]
    book18 = (r2 - r1) / (r2 + r1) * g / k - r2 * r1 / (r2 + r1) ** 2 * (U2 - U1) ** 2 + sg * k / (r1 + r2)  # (11.18) + σk
    coth_limit = sp.limit(Ch, h, sp.oo)  # = 1
    limit_res = sp.simplify(disc.subs(Cs, coth_limit) - book18)
    kk = sp.symbols("kk", positive=True)
    f = (g / kk) * (r2 - r1) + sg * kk
    ks = sp.solve(sp.diff(f, kk), kk)
    ks = [s for s in ks if s.is_positive is not False][0]
    dU2 = sp.simplify(f.subs(kk, ks) * (r1 + r2) / (r1 * r2))
    return dict(dispersion=dispersion, residual=res, limit_residual=limit_res, coth_limit=coth_limit, dU_min_sq=dU2, k_star=ks)


def kh_mixing_energy(U1: float, h: float, rho: float, profile: str | Callable = "linear") -> dict:
    """Kinetic energy and momentum per unit width before and after the shear layer is smeared out (p. 482, Fig. 11.7).

    Book: §11.3, p. 482: initially U = U₁ for 0 < z < h, 0 for −h < z < 0; after mixing U(z) = U₁(½ + z/(2h)) on −h ≤ z ≤ h;
    E_initial = (ρ/2)U₁²h, E_final = (ρ/2)∫U²dz = (ρ/3)U₁²h; momentum ∫ρU dz = ρU₁h unchanged.
    Parameters: U1 [m/s]; h [m] > 0; rho [kg/m³]; profile "linear" or a callable U_final(z) (should conserve ∫U dz).
    Returns dict(E_i, E_f [kg/s² per unit width], ratio (E_f/E_i), M_i, M_f [kg/(m s)], E_f_quad (quadrature cross-check)).
    Validation: V1 ratio 2/3; V4 momentum conserved; property: any ∫U-conserving smoothing lowers ∫U².  Label: analytic.
    """
    if h <= 0:
        raise ValueError("h must be > 0")
    E_i, M_i = 0.5 * rho * U1 ** 2 * h, rho * U1 * h
    if isinstance(profile, str):
        if profile != "linear":
            raise ValueError('profile must be "linear" or a callable')
        f = lambda z: U1 * (0.5 + z / (2.0 * h))  # noqa: E731  (p. 482)
        E_f, M_f = rho * U1 ** 2 * h / 3.0, rho * U1 * h
    else:
        f = profile
        E_f = 0.5 * rho * quad(lambda z: f(z) ** 2, -h, h, epsabs=1e-13)[0]
        M_f = rho * quad(f, -h, h, epsabs=1e-13)[0]
    E_q = 0.5 * rho * quad(lambda z: f(z) ** 2, -h, h, epsabs=1e-13, epsrel=1e-13)[0]
    return dict(E_i=E_i, E_f=E_f, ratio=E_f / E_i, M_i=M_i, M_f=M_f, E_f_quad=E_q)


# ======================================================================================================================
# §11.4 Bénard — scales and base state
# ======================================================================================================================
def rayleigh_number(alpha: float, dT, d: float, kappa: float, nu: float, g: float = G0):
    """Rayleigh number Ra = gαΓd⁴/(κν) with Γ = −dT̄/dz = ΔT/d, ΔT = T_bottom − T_top (= gαΔT d³/(κν)).

    Book: §11.4, Eq. (11.21) Ra = gαΓd⁴/κν, Γ = −dT̄/dz (positive when heated from below), (11.24) Γ ≡ ΔT/d.
    **Sign convention (slip S10)**: (11.21)'s Γ is −dT̄/dz — the *opposite* of Ch. 1's Kundu lapse rate Γ ≡ dT/dz and equal to
    meteorology's Γ_met.  This function never takes a Γ: ``dT`` = T_bottom − T_top (> 0 heated from below → Ra > 0).
    Parameters: alpha [1/K]; dT [K] (array ok); d [m]; kappa, nu [m²/s]; g [m/s²] (default G0; call by keyword).
    Returns Ra [–] (< 0 when heated from above).  Example: water d = 5 mm, ΔT = 2 K, α = 2.1e-4, κ = 1.4e-7, ν = 1e-6 →
    3679 (g = 9.81), 3677 (G0).  Validation: V2 pint dimensionless; V1 sign.  Label: analytic.
    """
    if d <= 0 or kappa <= 0 or nu <= 0:
        raise ValueError("d, kappa, nu must be > 0")
    Gamma = _F(dT) / d  # Γ = −dT̄/dz = ΔT/d (11.24)
    return _S(g * alpha * Gamma * d ** 4 / (kappa * nu))  # Eq. (11.21)


def gamma_conventions(dT, d: float) -> dict:
    """One layer's temperature gradient in the three sign conventions (slip S10).

    Book: §11.4 (11.21), (11.24) Γ_T = −dT̄/dz = ΔT/d; Ch. 1 (Kundu) Γ ≡ dT/dz; meteorology Γ_met = −dT/dz.
    Parameters: dT = T_bottom − T_top [K]; d [m].  Returns dict(Gamma_11_21, dTdz_kundu, Gamma_met) [K/m].
    Example: d = 5 mm, ΔT = 2 K → 400, −400, 400.  Label: analytic.
    """
    gr = _F(dT) / d
    return dict(Gamma_11_21=_S(gr), dTdz_kundu=_S(-gr), Gamma_met=_S(gr))


def benard_scales(alpha: float, dT: float, d: float, kappa: float, nu: float, g: float = G0) -> dict:
    """Reference scales of §11.4: Ra (11.21), Pr = ν/κ, Γ = ΔT/d, w ~ κ/d, t ~ d²/κ.

    Book: §11.4, p. 485–486 (w ~ κ/d; buoyant/viscous ~ Ra), the scaling above (11.31).  Parameters as
    :func:`rayleigh_number`.  Returns dict(Ra, Pr, Gamma [K/m], w_scale [m/s], t_scale [s]).  Example: water layer →
    w ~ 2.8e-5 m/s, t ~ 179 s.  Label: analytic.
    """
    return dict(Ra=rayleigh_number(alpha, dT, d, kappa, nu, g=g), Pr=nu / kappa, Gamma=dT / d, w_scale=kappa / d,
                t_scale=d ** 2 / kappa)


def benard_base_state(z, T0: float, dT: float, d: float, rho0: float = 1000.0, alpha: float = 2.1e-4, g: float = G0,
                      P0: float = 0.0) -> dict:
    """Conduction state of the Bénard layer: T̄(z), P(z) with z centred (−d/2 ≤ z ≤ d/2).

    Book: §11.4, Eq. (11.23) 0 = −∇P/ρ₀ − g[1 − α(T̄ − T₀)]e_z, 0 = κ∂²T̄/∂z²; Eq. (11.24) T̄ = T₀ − ½ΔT − Γz, Γ = ΔT/d.
    ``T0`` = **bottom** temperature (Fig. 11.8: T₀ below, T₀ − ΔT on top; T̄ = T₀ − Γ(z + d/2)).  Integrating (11.23):
    P = P₀ − ρ₀g[z + α(½ΔT z + ½Γz²)] (ours).
    Parameters: z [m]; T0 [K]; dT = T_bottom − T_top [K]; d [m]; rho0; alpha; g; P0 = P(0) [Pa].  Returns dict(T, P).
    Validation: V1 (11.23) residual; T̄(−d/2) = T₀, T̄(d/2) = T₀ − ΔT.  Label: analytic.
    """
    z = _F(z)
    Gamma = dT / d
    T = T0 - 0.5 * dT - Gamma * z  # Eq. (11.24)
    P = P0 - rho0 * g * (z + alpha * (0.5 * dT * z + 0.5 * Gamma * z ** 2))  # (11.23) integrated
    return dict(T=_S(T), P=_S(P))


# ======================================================================================================================
# §11.4 Bénard — Chebyshev eigen-solvers
# ======================================================================================================================
def _benard_ops(K: float, N: int):
    D, z = cheb(N, (-0.5, 0.5))  # z[0] = +½ (top), z[N] = −½ (bottom)
    n = N + 1
    I = np.eye(n)
    D2 = D @ D
    return D, D2, D2 - K ** 2 * I, I, z, n


def _benard_constraints(bc, D, D2, I, n, N):
    bottom, top = bc
    for b in (bottom, top):
        if b not in ("rigid", "free"):
            raise ValueError('boundary types must be "rigid" or "free" (bc = (bottom, top))')
    Zr = np.zeros(n)
    rows = [np.hstack([I[0], Zr]), np.hstack([I[N], Zr]),  # W = 0 (11.38)/(11.43)
            np.hstack([D[0] if top == "rigid" else D2[0], Zr]),  # W′ = 0 rigid (11.38) / W″ = 0 free (11.43)
            np.hstack([D[N] if bottom == "rigid" else D2[N], Zr]),
            np.hstack([Zr, I[0]]), np.hstack([Zr, I[N]])]  # T̂ = 0
    return np.array(rows), [0, N, 1, N - 1, n, n + N]


def _benard_marginal_solve(K: float, bc, N: int):
    D, D2, L, I, z, n = _benard_ops(K, N)
    Zm = np.zeros_like(I)
    # unknowns [W; T̂]; rows [W-equation; T-equation] so that row and unknown indices coincide
    A = np.block([[L @ L, Zm], [I, L]])  # Eq. (11.39): (D² − K²)²W = Ra K² T̂ ; (D² − K²)T̂ = −W
    B = np.block([[Zm, K ** 2 * I], [Zm, Zm]])
    C, elim = _benard_constraints(bc, D, D2, I, n, N)
    w, V = constrained_eig(A, B, C, elim, return_vectors=True, sort="abs")
    return w, V, z, n


def benard_marginal_Ra(K: float, bc: Sequence[str] = ("rigid", "rigid"), mode: str = "even", N: int = 40,
                       return_mode: bool = False):
    """Marginal Rayleigh number Ra(K) of the Bénard layer (σ = 0): smallest positive real eigenvalue of (11.39).

    Book: §11.4, Eq. (11.39) (d²/dz² − K²)T̂ = −W, (d²/dz² − K²)²W = RaK²T̂ (⇔ (11.40) (d²/dz² − K²)³W = −RaK²W) with rigid
    walls (11.38)/(11.41) W = W′ = T̂ = 0 or stress-free (11.43) W = W″ = T̂ = 0 on z = ±½ (lengths by d); Pr drops out.
    Parameters: K > 0 [–]; bc (bottom, top) each "rigid"/"free"; mode "even" (W symmetric: one row of cells, Fig. 11.9), "odd"
    (two rows) or "any"; N (default 40; FAST 24); return_mode (also dict(z, W, T), W max-normalised).
    Returns Ra [–].  Examples (rigid–rigid): K = 2 → 2177.41; 3.1163 → 1707.76; 5 → 2439.32.
    Validation: V5 minimum 1707.762 at 3.117 (Chandrasekhar 1961); V1 free–free (11.44); V1 :func:`benard_marginal_Ra_det`
    (independent route, 1e-9); V3 N = 24 → 48.  Label: converged.
    """
    if K <= 0:
        raise ValueError("K must be > 0")
    if mode not in ("even", "odd", "any"):
        raise ValueError('mode must be "even", "odd" or "any"')
    if mode != "any" and bc[0] != bc[1]:
        raise ValueError("parity (even/odd) needs symmetric boundaries; use mode='any' for rigid–free")
    w, V, z, n = _benard_marginal_solve(K, tuple(bc), N)
    ok = (np.abs(w.imag) <= 1e-8 * np.abs(w)) & (w.real > 0)
    if mode != "any":
        ok &= _parity_mask(V[:n], mode)
    if not ok.any():
        raise RuntimeError("benard_marginal_Ra: no positive real eigenvalue found (increase N)")
    idx = np.nonzero(ok)[0]
    j = idx[np.argmin(w.real[idx])]
    Ra = float(w.real[j])
    if not return_mode:
        return Ra
    s = V[:n, j][np.argmax(np.abs(V[:n, j]))]
    return Ra, dict(z=z, W=np.real(V[:n, j] / s), T=np.real(V[n:, j] / s))


def benard_critical(bc: Sequence[str] = ("rigid", "rigid"), mode: str = "even", N: int = 40,
                    K_bounds: Sequence[float] | None = None, xatol: float = 1e-7) -> dict:
    """Critical point = minimum over K of the marginal Ra(K).

    Book: §11.4 — rigid–rigid Ra_cr ≈ 1708 at K_cr ≈ 3.12 (λ_cr = 2πd/K_cr ≈ 2d, Fig. 11.10); free–free (27/4)π⁴ at K² = π²/2;
    rigid–free and the odd mode (Ex. 11.7) quoted (book values private).
    Returns dict(Ra_c, K_c).  Expect rigid–rigid 1707.762, 3.1163; free–free 657.511, 2.2214; rigid–free 1100.650, 2.6823;
    odd 17610.39, 5.3647.  Validation: V5 Chandrasekhar (1961); V1 free–free.  Label: converged, benchmark.
    """
    if K_bounds is None:
        K_bounds = (3.5, 7.5) if mode == "odd" else (1.5, 4.5)
    Ks = np.linspace(K_bounds[0], K_bounds[1], 9)
    vals = [benard_marginal_Ra(K, bc, mode, N) for K in Ks]
    i = int(np.argmin(vals))
    r = minimize_scalar(lambda K: benard_marginal_Ra(K, bc, mode, N), bounds=(Ks[max(i - 1, 0)], Ks[min(i + 1, 8)]),
                        method="bounded", options=dict(xatol=xatol))
    return dict(Ra_c=float(r.fun), K_c=float(r.x))


def benard_neutral_curve(Ks, bc: Sequence[str] = ("rigid", "rigid"), mode: str = "even", N: int = 40) -> np.ndarray:
    """Marginal Ra(K) on an array of K — the neutral curve of Fig. 11.10.  Book: §11.4.  Label: converged."""
    return np.array([benard_marginal_Ra(float(K), bc, mode, N) for K in np.atleast_1d(Ks)])


def benard_neutral_table(Ks=None, N: int = 40, write: bool = True, cache: bool = True) -> dict:
    """Rigid–rigid, free–free, rigid–free and odd-mode neutral curves on one K grid (explainer E3, figure F2).

    Book: §11.4, Fig. 11.10 and p. 491.  Default K grid ``np.linspace(0.5, 10, 96)`` (Part C.5 5.4).  Writes
    ``reference/ch11/benard_neutral_curves.csv`` (ours, ≤ 6 s.f.) when ``write``.  Returns dict(K, rigid, free, rigid_free, odd).
    Label: converged.
    """
    Ks = np.linspace(0.5, 10.0, 96) if Ks is None else _F(Ks)

    def compute():
        return dict(K=Ks, rigid=benard_neutral_curve(Ks, ("rigid", "rigid"), "even", N), free=_F(benard_free_free_Ra(Ks)),
                    rigid_free=benard_neutral_curve(Ks, ("rigid", "free"), "any", N),
                    odd=benard_neutral_curve(Ks, ("rigid", "rigid"), "odd", N))

    out = _cached("benard_neutral_table", dict(K=list(np.round(Ks, 10)), N=N), compute, cache)
    if write:
        p = _root() / "reference" / "ch11" / "benard_neutral_curves.csv"
        p.parent.mkdir(parents=True, exist_ok=True)
        with p.open("w", encoding="utf-8") as f:
            f.write("# ours (fluidpy.ch11_instability.benard_neutral_table), Chebyshev N = %d; not book data\n" % N)
            f.write("K,Ra_rigid_rigid,Ra_free_free,Ra_rigid_free,Ra_odd_rigid_rigid\n")
            for row in zip(out["K"], out["rigid"], out["free"], out["rigid_free"], out["odd"]):
                f.write(",".join(f"{v:.6g}" for v in row) + "\n")
    return out


def _benard_growth_solve(K, Ra, Pr, bc, N):
    D, D2, L, I, z, n = _benard_ops(K, N)
    Zm = np.zeros_like(I)
    A = np.block([[L @ L, -Ra * K ** 2 * I], [I, L]])  # Eqs. (11.36)–(11.37)
    B = np.block([[L / Pr, Zm], [Zm, I]])
    C, elim = _benard_constraints(tuple(bc), D, D2, I, n, N)
    return constrained_eig(A, B, C, elim, sort="real")


def _benard_confirm_ladder(N: int) -> tuple:
    """Resolutions tried, in order, to confirm a leading Bénard eigenvalue that the 1.5N solve did not confirm."""
    return int(math.ceil(1.25 * N)), int(math.ceil(0.75 * N)), int(math.ceil(0.6 * N))


def benard_growth_rate(K: float, Ra: float, Pr: float, bc: Sequence[str] = ("rigid", "rigid"), N: int = 40,
                       all: bool = False, return_complex: bool = False, filter: bool = True):
    """Growth rate σ (units κ/d²) of the leading Bénard normal mode at wavenumber K (or all, sorted by Re σ).

    Book: §11.4, Eqs. (11.36)–(11.37) (σ + K² − d²/dz²)T̂ = W, (σ/Pr + K² − d²/dz²)(d²/dz² − K²)W = −RaK²T̂ with (11.38)
    (rigid) or (11.43) (free); σ real for Ra > 0 (Ex. 11.6, :func:`exchange_of_stabilities_sympy`).
    Generalised eigenproblem in [W; T̂]: L²W − RaK²T̂ = (σ/Pr)LW, W + LT̂ = σT̂ (L = d²/dz² − K²).
    Boundary unknowns are eliminated (no spurious boundary eigenvalues) and, with ``filter`` (default), only σ that reappear
    in a solve at ⌈1.5N⌉ within 1e-6·max(1, |σ|) are kept (design note: naïve row replacement gives spurious huge σ at N ≥ 32).
    Parameters: K, Ra, Pr; bc; N; all (return all, sorted); return_complex; filter.  Returns float (leading real σ) or ndarray.
    Ra < 0 (heated from above, stably stratified): exchange of stabilities no longer holds — the damped modes can be complex
    pairs (decaying internal-wave oscillations).  The default then returns **Re σ of the least-damped mode and drops Im σ**
    (``all=True`` likewise returns real parts); pass ``return_complex=True`` to see it.  Measured, rigid–rigid, K = 3:
    Ra = −2000, Pr = 1 → −28.6456 + 26.9938i (default: −28.6456); Ra = −5000, Pr = 0.7 → −23.3965 ∓ 37.8442i; but
    Ra = −2000, Pr = 7 → −44.9839 (real: the least-damped mode is still a real one).  The sign of Im σ of a pair is arbitrary.
    The leading mode is never replaced by a lower one: if the largest-Re σ of the N solve is not confirmed at ⌈1.5N⌉ it is
    re-tested at ⌈1.25N⌉, ⌈0.75N⌉, ⌈0.6N⌉ (same tolerance; the round-off of the Chebyshev D⁴ grows steeply with N, so the
    finer solve can be the less accurate one — rigid–free, K = 1.84, Ra = 1.6 Ra(K), Pr = 7.14: 7.1195004794 exact,
    7.1195004809 at N = 40, 7.1195096 at N = 60), and if it is still unconfirmed a ValueError is raised.  (Before this
    check the default returned the *second* eigenvalue there, −39.61, and at 11 of 2000 points of a (K, Ra/Ra(K), Pr) scan.)
    Raises ValueError if no eigenvalue survives the N-filter (N too small to resolve any mode, e.g. N = 8), or if the
    leading eigenvalue does not while a lower one does (measured once in 1080 wide-range points: rigid–rigid K = 5,
    Ra = 5e4, Pr = 100, where round-off is 1.3e-6 of σ = 605.29) — change N or pass ``filter=False`` (the unfiltered leading
    eigenvalue).  ``all=True`` with nothing confirmed returns an empty array.  Use N ≈ 24 – 48: accuracy is ≤ 1e-7 there and
    degrades above (rigid–free, N = 96: 5e-6 relative, shared by neighbouring N and so not detectable by the filter).
    Validation: V1 free–free = :func:`benard_free_free_sigma` (1e-8); V1 (independent route) roots of the 6 × 6 determinant
    of the exponential solutions, rigid–rigid and rigid–free (≤ 5e-8 on 61 points); all |σ_i| < 1e-10 for Ra > 0; σ₁ = 0 at
    :func:`benard_marginal_Ra` and sign σ₁ = sign(Ra − Ra(K)) on 2000 points; V3.  Label: converged.
    """
    if K <= 0 or Pr <= 0:
        raise ValueError("K > 0 and Pr > 0 required")
    w = _benard_growth_solve(K, Ra, Pr, bc, N)
    if filter and len(w):
        M = int(math.ceil(1.5 * N))
        mask = converged_mask(w, _benard_growth_solve(K, Ra, Pr, bc, M), 1e-6)
        tried = [M]
        if not mask[0]:
            # The leading candidate w[0] (largest Re σ) is not confirmed at 1.5N.  It must never be skipped silently in
            # favour of a lower mode: confirm it (and its complex partner / near neighbours above the first confirmed
            # mode) at other resolutions — the round-off of the Chebyshev D⁴ grows steeply with N, so the finer solve can
            # be the less accurate one — or raise.
            n_top = int(np.argmax(mask)) if mask.any() else len(w)
            for M2 in _benard_confirm_ladder(N):
                if M2 in tried or M2 < 6:
                    continue
                tried.append(M2)
                mask[:n_top] |= converged_mask(w[:n_top], _benard_growth_solve(K, Ra, Pr, bc, M2), 1e-6)
                if mask[0]:
                    break
        if not mask[0] and (mask.any() or not all):
            what = ("no eigenvalue" if not mask.any() else
                    f"the leading eigenvalue σ = {complex(w[0]):.9g} (a lower mode at {complex(w[int(np.argmax(mask))]):.9g} "
                    "does, and is not returned in its place)")
            raise ValueError(f"benard_growth_rate: {what} of the N = {N} solve does not reappear within 1e-6·max(1, |σ|) at "
                             f"N = {tried} (resolution too low — or too high: round-off — for the N-convergence filter) — "
                             "change N or pass filter=False")
        w = w[mask]
    if len(w) == 0 and not all:
        raise ValueError(f"benard_growth_rate: the N = {N} solve has no finite eigenvalue (filter={filter}) — increase N")
    if all:
        return w if return_complex else np.real(w)
    return complex(w[0]) if return_complex else float(np.real(w[0]))


def benard_char_roots(Ra: float, K: float):
    """Roots of (q² − K²)³ = −RaK² (characteristic equation of (11.40)): (q₀, q, q*); the six roots are ±iq₀, ±q, ±q*.

    Book: §11.4, p. 488, Eq. (11.42) q² = −K²[(Ra/K⁴)^{1/3} − 1], q² = K²[1 + ½(Ra/K⁴)^{1/3}(1 ± i√3)]; p. 489
    q₀ = K[(Ra/K⁴)^{1/3} − 1]^{1/2} (real for Ra > K⁴).  q is the principal square root (Re q > 0).
    Returns (q0 (float if real, else complex), q (complex), q_conj).  Example: K = 2, Ra = 1000 → q₀ = 3.4459,
    q = 3.8822 + 1.7705i.  Label: analytic.
    """
    r = (Ra / K ** 4) ** (1.0 / 3.0)
    q0 = K * np.lib.scimath.sqrt(r - 1.0)  # p. 489
    q = np.sqrt(complex(K ** 2 * (1.0 + 0.5 * r * (1.0 + 1j * math.sqrt(3.0)))))  # Eq. (11.42), "+" root
    q0 = float(np.real(q0)) if np.imag(q0) == 0 else complex(q0)
    return q0, complex(q), complex(np.conj(q))


def benard_determinant(Ra: float, K: float, mode: str = "even", printed: bool = False) -> complex:
    """Determinant of the 3 × 3 rigid-wall system for the even (or odd) solution of (11.40).

    Book: §11.4, p. 489: W = A cos q₀z + B cosh qz + C cosh q*z; rows W(½), W′(½), (d²/dz² − K²)²W(½):
        | cos(q₀/2)               cosh(q/2)              cosh(q*/2)            |
        | −q₀ sin(q₀/2)           q sinh(q/2)            q* sinh(q*/2)         |  = 0
        | (q₀² + K²)² cos(q₀/2)   (q² − K²)² cosh(q/2)   (q*² − K²)² cosh(q*/2) |
    (z ∈ [−½, ½]; rows = BCs at z = +½; q principal, Re q > 0).  Odd (Ex. 11.7): sin q₀z, sinh qz, sinh q*z.
    Columns 2 and 3 are complex conjugates ⇒ det is purely imaginary (use Im).
    ``printed=True`` (slip S1): the derivative line on p. 489 attaches the q* term to B, i.e. row 3 becomes
    [(q₀² + K²)² cos(q₀/2), (q² − K²)² cosh(q/2) + (q*² − K²)² cosh(q*/2), 0] (the page's matrix is right).  Measured shift of
    the marginal Ra: 0.44 % at K = 3.116, 1.6 % at K = 2, 2.4 % at K = 5.
    Parameters: Ra (> K⁴ for real q₀); K; mode; printed.  Returns complex.  Example: K = 3.1163, Ra = 1707.762 → ≈ 0.0318i.
    Label: analytic.
    """
    q0, q, qs = benard_char_roots(Ra, K)
    if mode == "even":
        f0, g0 = np.cos(q0 / 2), -q0 * np.sin(q0 / 2)
        f1, g1 = np.cosh(q / 2), q * np.sinh(q / 2)
        f2, g2 = np.cosh(qs / 2), qs * np.sinh(qs / 2)
    elif mode == "odd":
        f0, g0 = np.sin(q0 / 2), q0 * np.cos(q0 / 2)
        f1, g1 = np.sinh(q / 2), q * np.cosh(q / 2)
        f2, g2 = np.sinh(qs / 2), qs * np.cosh(qs / 2)
    else:
        raise ValueError('mode must be "even" or "odd"')
    h0 = (q0 ** 2 + K ** 2) ** 2 * f0
    h1 = (q ** 2 - K ** 2) ** 2 * f1
    h2 = (qs ** 2 - K ** 2) ** 2 * f2
    if printed:  # slip S1
        M = np.array([[f0, f1, f2], [g0, g1, g2], [h0, h1 + h2, 0.0]], dtype=complex)
    else:
        M = np.array([[f0, f1, f2], [g0, g1, g2], [h0, h1, h2]], dtype=complex)  # p. 489
    return complex(np.linalg.det(M))


def benard_marginal_Ra_det(K: float, mode: str = "even", Ra_max: float = 2e4, n_scan: int = 400,
                           printed: bool = False) -> float:
    """Marginal Ra(K) for rigid walls by Brent on Im det (:func:`benard_determinant`) — independent of Chebyshev.

    Book: §11.4, p. 489 (nonzero (A, B, C) ⇔ det = 0).  Scan Ra geometrically from max(1.02K⁴, 50) (q₀ real; the degenerate
    end Ra → K⁴ has spurious sign changes) to Ra_max (raised automatically to 20K⁴ if smaller), first sign change, brentq.
    Returns Ra.  Raises ValueError if K ≤ 0 or if Im det does not change sign below the upper end of the scan — never a
    silent NaN.  Reach of the default ``Ra_max = 2e4`` (measured, K = 0.3 … 9 in steps of 0.1): the **even** root is found
    for K ≥ 0.6 (K = 0.5: root 21 009.8; K = 0.3: 56 991.6); the **odd** root only for K ≥ 4.0 (K = 3.9: 20 061.1; K = 3:
    26 146.6; K = 2: 47 005.6; K = 1.2: 115 710.9; K = 1: 163 127.6; K = 0.5: 629 152.8).  Below those K pass ``Ra_max``
    above the root; the error message computes it for you from the Chebyshev route (:func:`benard_marginal_Ra`) and names
    a sufficient ``Ra_max`` (1.2 × that value; checked to work at each K listed).  The scan is not widened automatically:
    a call that needs a larger bracket says so.
    Validation: V1 equals :func:`benard_marginal_Ra` to 1e-8 at K = 2, 3.1163, 5.
    Label: analytic.
    """
    if not K > 0:
        raise ValueError("benard_marginal_Ra_det: K > 0 is required")
    lo = max(1.02 * K ** 4, 50.0)
    hi = max(Ra_max, 20.0 * K ** 4, 2.0 * lo)
    grid = np.geomspace(lo, hi, int(n_scan))
    f = lambda R: benard_determinant(R, K, mode, printed).imag  # noqa: E731
    v = np.array([f(R) for R in grid])
    for j in range(len(grid) - 1):
        if np.isfinite(v[j]) and np.isfinite(v[j + 1]) and v[j] * v[j + 1] < 0:
            return float(brentq(f, grid[j], grid[j + 1], xtol=1e-10, rtol=1e-14))
    reach = "K ≥ 0.6" if mode == "even" else "K ≥ 4.0"
    try:  # per-mode hint: where the independent Chebyshev route puts the root
        est = float(benard_marginal_Ra(K, mode=mode))
        hint = (f"the Chebyshev route puts the {mode} root at Ra ≈ {est:.6g}: pass Ra_max ≥ {1.2 * est:.3g}"
                if np.isfinite(est) and est > hi else
                f"the Chebyshev route gives Ra ≈ {est:.6g}, inside the scan — try a larger n_scan")
    except Exception:  # the hint must never mask the real error
        hint = "pass Ra_max above the root"
    raise ValueError(f"benard_marginal_Ra_det: no sign change of Im det for K = {K:g}, mode = {mode!r} in "
                     f"Ra ∈ [{lo:.4g}, {hi:.4g}] — raise Ra_max (the default 2e4 reaches the {mode} mode only for "
                     f"{reach}; {hint})")


def benard_free_free_Ra(K, n: int = 1):
    """Marginal Rayleigh number between stress-free surfaces: Ra = (n²π² + K²)³/K².

    Book: §11.4, Eq. (11.44) (with W = A sin nπ(z + ½); slip S2: the page prints sin nπz).  Parameters: K > 0 (array ok); n.
    Returns Ra.  Examples: Ra(1) = 1284.23, Ra(2) = 667.01, Ra(3) = 746.53, Ra(4) = 1082.06.  Label: analytic.
    """
    K = _F(K)
    return _S((n ** 2 * math.pi ** 2 + K ** 2) ** 3 / K ** 2)  # Eq. (11.44)


def benard_free_free_critical() -> dict:
    """Critical point for stress-free surfaces: dict(Ra_c = 27π⁴/4 = 657.511, K_c = π/√2 = 2.2214, K_c2 = π²/2).

    Book: §11.4, p. 490–491: dRa/dK² = 3(π² + K²)²/K² − (π² + K²)³/K⁴ = 0 (slip S3: the page prints a factor 3 on the
    second term).  Label: analytic.
    """
    return dict(Ra_c=RA_FREE_FREE, K_c=math.pi / math.sqrt(2.0), K_c2=math.pi ** 2 / 2.0)


def benard_free_free_mode(z, n: int = 1, printed: bool = False):
    """Free–free eigenfunction W(z) = sin nπ(z + ½), −½ ≤ z ≤ ½ (n = 1: cos πz, as in Lorenz's (11.90)).

    Book: §11.4, p. 490.  ``printed=True`` returns the page's sin(nπz) — non-zero at z = ±½ for odd n (slip S2).
    Label: analytic.
    """
    z = _F(z)
    return _S(np.sin(n * math.pi * z) if printed else np.sin(n * math.pi * (z + 0.5)))


def benard_free_free_sigma(K, Ra: float, Pr: float, n: int = 1):
    """The two growth rates (units κ/d²) of the free–free Bénard mode n — real (σ₊ ≥ σ₋) for Ra ≥ 0.

    Book: §11.4, (11.36)–(11.37) with W = sin nπ(z + ½) (our D12): a² = n²π² + K²,
    (σ + a²)(σ/Pr + a²)a² = RaK²  ⇔  σ² + a²(1 + Pr)σ + Pr a⁴ − Pr Ra K²/a² = 0; real for Ra > 0; σ₊ = 0 at (11.44).
    Returns (sigma_plus, sigma_minus): two floats when the discriminant a⁴(1 − Pr)² + 4 Pr Ra K²/a² is ≥ 0 (always for
    Ra ≥ 0), otherwise a **complex-conjugate pair** (σ₊ carries Im > 0): for Ra < −a⁶(1 − Pr)²/(4 Pr K²) the stably
    stratified layer answers with a damped oscillation (Re σ = −a²(1 + Pr)/2 < 0).  At K = π/√2, n = 1 that threshold is
    Ra = −845.4 (Pr = 7), −21.13 (Pr = 0.7) and 0 (Pr = 1).
    Examples: K = π/√2, Ra = 2000, Pr = 1 → 11.0155, −40.6243; K = 2.2, Ra = −2000, Pr = 7 → −58.8384 ± 51.5671i;
    K = 2.2, Ra = −50, Pr = 7 → −16.0343, −101.6425 (still real).  Label: analytic.
    """
    a2 = n ** 2 * math.pi ** 2 + K ** 2
    b = a2 * (1.0 + Pr)
    c = Pr * a2 ** 2 - Pr * Ra * K ** 2 / a2
    disc = np.lib.scimath.sqrt(b ** 2 - 4.0 * c)
    r1, r2 = (-b + disc) / 2.0, (-b - disc) / 2.0  # D12
    cv = lambda v: float(np.real(v)) if np.imag(v) == 0 else complex(v)  # noqa: E731
    return cv(r1), cv(r2)


def benard_eigenfunction(K: float, bc: Sequence[str] = ("rigid", "rigid"), mode: str = "even", x=None, z=None,
                         N: int = 40) -> dict:
    """Marginal Bénard roll: profiles W(z), T̂(z) of (11.39) and the 2-D roll fields on an (x, z) grid.

    Book: §11.4, Fig. 11.9 (even: one row of cells; odd: two rows), Fig. 11.11 (rolls).  Rolls along y:
    w = W(z)cos Kx, T′ = T̂(z)cos Kx, u = −W′(z) sin(Kx)/K, ψ = W(z) sin(Kx)/K with u = −∂ψ/∂z, w = ∂ψ/∂x (the §11.14 sign).
    Normalised so max|w| = 1 (W > 0 at its maximum).  Defaults: x over two wavelengths (81 points), z ∈ [−½, ½] (41 points).
    Returns dict(x, z (2-D grids), psi, w, u, T, W_profile, T_profile, z_nodes, Ra).  Label: converged.
    """
    md = mode if bc[0] == bc[1] else "any"
    Ra, m = benard_marginal_Ra(K, bc, md, N, return_mode=True)
    zc, Wc, Tc = m["z"], m["W"], m["T"]
    x = np.linspace(0.0, 4.0 * math.pi / K, 81) if x is None else _F(x)
    z = np.linspace(-0.5, 0.5, 41) if z is None else _F(z)
    X, Zg = np.meshgrid(x, z, indexing="xy") if (x.ndim == 1 and z.ndim == 1) else np.broadcast_arrays(x, z)
    Wi, Ti = BarycentricInterpolator(zc, Wc), BarycentricInterpolator(zc, Tc)
    dWi = BarycentricInterpolator(zc, cheb(N, (-0.5, 0.5))[0] @ Wc)
    Wz, Tz, dWz = Wi(Zg), Ti(Zg), dWi(Zg)
    return dict(x=X, z=Zg, psi=Wz * np.sin(K * X) / K, w=Wz * np.cos(K * X), u=-dWz * np.sin(K * X) / K,
                T=Tz * np.cos(K * X), W_profile=Wc, T_profile=Tc, z_nodes=zc, Ra=Ra)


def planform(x, y, K: float, kind: str = "rolls"):
    """Horizontal planform f(x, y) with ∇_H²f = −K²f: rolls cos Kx; squares cos Kx + cos Ky; hexagons Σ_{j=1}^{3}cos(k_j·x)
    with |k_j| = K at 0°, 120°, 240°.

    Book: §11.4, p. 491 (linear theory fixes |K| only; Figs. 11.11–11.12).  Validation: V1 ∇_H²f + K²f = 0.  Label: analytic.
    """
    x, y = _F(x), _F(y)
    if kind == "rolls":
        return _S(np.cos(K * x) + 0.0 * y)
    if kind == "squares":
        return _S(np.cos(K * x) + np.cos(K * y))
    if kind == "hexagons":
        out = 0.0
        for th in (0.0, 2 * math.pi / 3, 4 * math.pi / 3):
            out = out + np.cos(K * (math.cos(th) * x + math.sin(th) * y))
        return _S(out)
    raise ValueError('kind must be "rolls", "squares" or "hexagons"')


# ======================================================================================================================
# §11.4 Bénard — sympy engines
# ======================================================================================================================
def benard_perturbation_sympy() -> dict:
    """D05–D09: (11.25)–(11.29), the scaling to (11.31)–(11.32), normal modes (11.34)–(11.37), (11.39) → (11.40).

    Book: §11.4.  (a) With R_i = ∂_t u_i + ∂_i p/ρ₀ − gαT′δ_{i3} − ν∇²u_i (residual of (11.26)) and D = ∇·u:
    ∇²R_z − ∂_z(∇·R) − [∂_t∇²w − gα∇_H²T′ − ν∇⁴w] = −∂_t∂_zD + ν∇²∂_zD, i.e. zero when (11.25) holds ⇒ (11.29).  The planted
    variant with ∇² in place of ∇_H² leaves a non-zero remainder.  (b) t → (d²/κ)t, x → xd turns (11.27), (11.29) into
    (11.31)–(11.32).  (c) Normal modes and W ≡ (Γd²/κ)ŵ give (11.36)–(11.37) with Ra of (11.21).  (d) σ = 0, eliminate T̂ →
    (11.40) (constant coefficients; checked on W = e^{qz}).
    Returns dict(eq_11_27 (the scaled residual, 0), eq_11_29_residual (0), planted_residual (≠ 0), eq_11_36 (0), eq_11_37 (0),
    eq_11_40_residual (0), char_eq).  Validation: V2.  Label: symbolic.
    """
    x, y, z, t = sp.symbols("x y z t", real=True)
    rho0, g, al, nu, ka, Gam, d, Pr = sp.symbols("rho0 g alpha nu kappa Gamma d Pr", positive=True)
    u, v, w, p, T = [sp.Function(n)(x, y, z, t) for n in ("u", "v", "w", "p", "T")]
    lap = lambda f: sp.diff(f, x, 2) + sp.diff(f, y, 2) + sp.diff(f, z, 2)  # noqa: E731
    lapH = lambda f: sp.diff(f, x, 2) + sp.diff(f, y, 2)  # noqa: E731
    R = [sp.diff(q, t) + sp.diff(p, s) / rho0 - nu * lap(q) for q, s in ((u, x), (v, y), (w, z))]  # Eq. (11.26)
    R[2] = R[2] - g * al * T
    Dv = sp.diff(u, x) + sp.diff(v, y) + sp.diff(w, z)  # Eq. (11.25)
    divR = sp.diff(R[0], x) + sp.diff(R[1], y) + sp.diff(R[2], z)
    extra = -sp.diff(Dv, t, z) + nu * sp.diff(lap(Dv), z)

    def remainder(planted):
        e29 = sp.diff(lap(w), t) - g * al * (lap(T) if planted else lapH(T)) - nu * lap(lap(w))  # Eq. (11.29)
        return sp.simplify(sp.expand(lap(R[2]) - sp.diff(divR, z) - e29 - extra))

    sig, K, q, Ra, m2 = sp.symbols("sigma K q Ra m2")
    lap_s = -(K ** 2 + m2)  # ∇² → −K² + d²/dz², d²/dz² → −m²
    Th, wh, W = sp.symbols("That what W")
    eq27 = (ka / d ** 2) * sig * Th - wh * Gam - ka * lap_s * Th / d ** 2  # (11.27), dimensional, normal modes
    eq31 = sig * Th - lap_s * Th - Gam * d ** 2 / ka * wh  # Eq. (11.31)
    sc_T = sp.simplify(eq27 * d ** 2 / ka - eq31)
    eq29 = (ka / d ** 2) * sig * lap_s * wh / d ** 2 - g * al * (-K ** 2) * Th / d ** 2 - nu * lap_s ** 2 * wh / d ** 4
    eq32 = (sig / Pr - lap_s) * lap_s * wh - g * al * d ** 2 / nu * (-K ** 2) * Th  # Eq. (11.32)
    sc_w = sp.simplify((eq29 * d ** 4 / nu).subs(ka, nu / Pr) - eq32)
    nm_T = sp.simplify(eq31.subs(wh, W * ka / (Gam * d ** 2)) - ((sig - lap_s) * Th - W))  # (11.36)
    eq37 = (sig / Pr - lap_s) * lap_s * W + Ra * K ** 2 * Th  # (11.37)
    nm_W = sp.simplify((eq32 * Gam * d ** 2 / ka).subs(wh, W * ka / (Gam * d ** 2)).subs(nu, Pr * ka)
                       - eq37.subs(Ra, g * al * Gam * d ** 4 / (ka * Pr * ka)))
    Lq = q ** 2 - K ** 2
    That = -W / Lq
    sixth = sp.simplify(Lq ** 2 * W - Ra * K ** 2 * That - (Lq ** 3 * W + Ra * K ** 2 * W) / Lq)
    return dict(eq_11_27=sc_T, eq_11_29_residual=remainder(False), planted_residual=remainder(True), eq_11_32=sc_w,
                eq_11_36=nm_T, eq_11_37=nm_W, eq_11_40_residual=sixth, char_eq=sp.expand(Lq ** 3 + Ra * K ** 2))


def benard_free_free_sympy() -> dict:
    """D11: W = sin nπ(z + ½) satisfies (11.40) with (11.43) iff Ra = (n²π² + K²)³/K² (11.44); dRa/dK² = 0 ⇒ K² = π²/2,
    Ra = 27π⁴/4.  Printed slips: sin nπz fails W(±½) = 0 for n = 1 (S2); dRa/dK² with the factor 3 (S3) has no positive root.

    Book: §11.4, p. 490–491.  Returns dict(bc_W4 (W, W″, W⁗ at ±½: all 0), residual (of (11.40): 0), dRa_dK2, root (π²/2),
    Ra_c (27π⁴/4), printed_dRa_dK2, printed_root ([]), printed_bc_n1 ([1, −1])).  Validation: V2.  Label: symbolic.
    """
    z, K2 = sp.symbols("z K2", positive=True)
    n = sp.symbols("n", positive=True, integer=True)
    W = sp.sin(n * sp.pi * (z + sp.Rational(1, 2)))
    L = lambda f: sp.diff(f, z, 2) - K2 * f  # noqa: E731
    Ra_n = (n ** 2 * sp.pi ** 2 + K2) ** 3 / K2  # Eq. (11.44)
    res = sp.simplify(L(L(L(W))) + Ra_n * K2 * W)  # Eq. (11.40)
    bcs = [sp.simplify(sp.diff(W, z, m).subs(z, s)) for m in (0, 2, 4) for s in (sp.Rational(1, 2), -sp.Rational(1, 2))]
    printed_bc = [sp.sin(sp.pi * s) for s in (sp.Rational(1, 2), -sp.Rational(1, 2))]
    Ra1 = Ra_n.subs(n, 1)
    dRa = sp.simplify(sp.diff(Ra1, K2))
    book_form = 3 * (sp.pi ** 2 + K2) ** 2 / K2 - (sp.pi ** 2 + K2) ** 3 / K2 ** 2  # p. 490, corrected
    roots = [r for r in sp.solve(sp.Eq(book_form, 0), K2) if r.is_positive]
    printed_form = 3 * (sp.pi ** 2 + K2) ** 2 / K2 - 3 * (sp.pi ** 2 + K2) ** 3 / K2 ** 2  # slip S3
    printed_roots = [r for r in sp.solve(sp.Eq(printed_form, 0), K2) if r.is_positive]
    return dict(bc_W4=bcs, residual=res, dRa_dK2=dRa, dRa_minus_book=sp.simplify(dRa - book_form),
                root=roots[0] if roots else None, Ra_c=sp.simplify(Ra1.subs(K2, roots[0])) if roots else None,
                printed_dRa_dK2=printed_form, printed_root=printed_roots, printed_bc_n1=printed_bc)


def exchange_of_stabilities_sympy() -> dict:
    """D08 (Exercise 11.6): σ is real for Ra > 0 — the two energy relations and the imaginary-part identity.

    Book: Ex. 11.6: (a) ×T̂*, integrate (11.36): σI₁ + I₂ = ∫T̂*W; (b) ×W*, integrate (11.37): (σ/Pr)J₁ + J₂ = RaK²∫W*T̂, with
    I₁ = ∫|T̂|², I₂ = ∫(|T̂′|² + K²|T̂|²), J₁ = ∫(|W′|² + K²|W|²), J₂ = ∫(|W″|² + 2K²|W′|² + K⁴|W|²) (> 0, by parts with (11.38));
    (c) RaK² × conj(a) − (b): imaginary part −σ_i(RaK²I₁ + J₁/Pr) = 0 ⇒ σ_i = 0.
    Built on polynomial trial functions meeting (11.38): T̂ = (z² − ¼)(t₀ + t₁z + t₂z²), W = (z² − ¼)²(w₀ + w₁z + w₂z²),
    complex coefficients.  Returns dict(relation_T (0), relation_W (0), imag_identity (0), sigma_i_solutions ([0])).
    Validation: V2.  Label: symbolic.
    """
    z = sp.symbols("z", real=True)
    K, Pr, Ra = sp.symbols("K Pr Ra", positive=True)
    sr, si = sp.symbols("sigma_r sigma_i", real=True)
    sig = sr + sp.I * si
    tc = [sp.symbols(f"t{j}r t{j}i", real=True) for j in range(3)]
    wc = [sp.symbols(f"w{j}r w{j}i", real=True) for j in range(3)]
    T = (z ** 2 - sp.Rational(1, 4)) * sum((a + sp.I * b) * z ** j for j, (a, b) in enumerate(tc))
    W = (z ** 2 - sp.Rational(1, 4)) ** 2 * sum((a + sp.I * b) * z ** j for j, (a, b) in enumerate(wc))
    Tc, Wcj = sp.conjugate(T), sp.conjugate(W)
    integ = lambda f: sp.integrate(sp.expand(f), (z, -sp.Rational(1, 2), sp.Rational(1, 2)))  # noqa: E731
    I1 = integ(Tc * T)
    I2 = integ(sp.diff(Tc, z) * sp.diff(T, z) + K ** 2 * Tc * T)
    J1 = integ(sp.diff(Wcj, z) * sp.diff(W, z) + K ** 2 * Wcj * W)
    J2 = integ(sp.diff(Wcj, z, 2) * sp.diff(W, z, 2) + 2 * K ** 2 * sp.diff(Wcj, z) * sp.diff(W, z) + K ** 4 * Wcj * W)
    lhsT = integ(Tc * (sig * T + K ** 2 * T - sp.diff(T, z, 2)))  # ∫T̂*(σ + K² − D²)T̂
    relT = sp.simplify(sp.expand(lhsT - (sig * I1 + I2)))
    LW = sp.diff(W, z, 2) - K ** 2 * W
    lhsW = integ(Wcj * ((sig / Pr) * LW - (sp.diff(LW, z, 2) - K ** 2 * LW)))  # ∫W*(σ/Pr + K² − D²)(D² − K²)W
    relW = sp.simplify(sp.expand(lhsW + (sig / Pr) * J1 + J2))  # = −[(σ/Pr)J₁ + J₂]
    i1, i2, j1, j2 = sp.symbols("I1 I2 J1 J2", positive=True)
    eq = Ra * K ** 2 * (sp.conjugate(sig) * i1 + i2) - (sig / Pr * j1 + j2)
    im = sp.simplify(sp.im(sp.expand(eq)))
    ident = sp.simplify(im + si * (Ra * K ** 2 * i1 + j1 / Pr))
    return dict(relation_T=relT, relation_W=relW, imag_identity=ident, imag_part=im,
                sigma_i_solutions=sp.solve(sp.Eq(im, 0), si))


# ======================================================================================================================
# §11.5 double diffusion
# ======================================================================================================================
def linear_eos(T, S, rho0: float = 1027.0, alpha: float = 2e-4, beta_S: float = 7.6e-4, T0: float = 10.0, S0: float = 35.0):
    """Linear equation of state ρ = ρ₀[1 − α(T − T₀) + β(S − S₀)] (reuses Ch. 1's ``seawater_density_linear``).

    Book: §11.5, p. 492.  Parameters: T and T0 in the same scale (default T₀ = 10 °C; only differences matter); S, S0 [g/kg];
    rho0 [kg/m³]; alpha [1/K]; beta_S [1/(g/kg)].  Returns ρ [kg/m³].  Label: analytic.
    """
    from .ch01_introduction import seawater_density_linear

    return seawater_density_linear(T, S, rho0=rho0, alpha_T=alpha, beta_S=beta_S, T0=T0, S0=S0)


def thermal_rayleigh_signed(dTdz, d: float, alpha: float, nu: float, kappa: float, g: float = G0):
    """§11.5's signed thermal Rayleigh number Ra ≡ gαd⁴(dT̄/dz)/(νκ) — *negative* when heated from below (slip S10).

    Book: §11.5, p. 494.  Parameters: dTdz [K/m]; d [m]; alpha [1/K]; nu, kappa [m²/s]; g.  Returns Ra.
    Example: dT/dz = 0.01 K/m, d = 5 cm (α = 2e-4, ν = 1e-6, κ = 1.4e-7) → 875.6 (G0).  Label: analytic.
    """
    return _S(g * alpha * d ** 4 * _F(dTdz) / (nu * kappa))


def salinity_rayleigh(dSdz, d: float, beta_S: float, nu: float, kappa_s: float, g: float = G0):
    """Rs ≡ gβd⁴(dS/dz)/(νκ_s) (κ_s, not κ: Rs = (κ/κ_s)Rs′).  Book: §11.5, p. 494.  Positive salty-over-fresh.
    Example: dS/dz = 0.002 (g/kg)/m, d = 5 cm, β = 7.6e-4, κ_s = 1.5e-9 → 62 110 (G0).  Label: analytic.
    """
    return _S(g * beta_S * d ** 4 * _F(dSdz) / (nu * kappa_s))


def salinity_rayleigh_prime(dSdz, d: float, beta_S: float, nu: float, kappa: float, g: float = G0):
    """Rs′ ≡ gβd⁴(dS/dz)/(νκ) of (11.45) (κ, not κ_s).  Book: §11.5, p. 494.  Label: analytic."""
    return _S(g * beta_S * d ** 4 * _F(dSdz) / (nu * kappa))


def double_diffusive_margin(Ra, Rs):
    """Rs − Ra − 27π⁴/4 (> 0 ⇒ finger-unstable, free–free).  Book: §11.5, p. 494 (T̂ = κ_sŝ/κ turns (11.45) into (11.39)
    with Ra → Rs − Ra).  Ra in the §11.5 sign, Rs with κ_s.  Label: analytic.
    """
    return _S(_F(Rs) - _F(Ra) - RA_FREE_FREE)


def salt_finger_unstable(dTdz, dSdz, d: float, alpha: float = 2e-4, beta_S: float = 7.6e-4, nu: float = 1e-6,
                         kappa: float = 1.4e-7, kappa_s: float = 1.5e-9, g: float = G0) -> dict:
    """Finger criterion (11.46) (free–free) and the static stability of the same column.

    Book: §11.5, Eq. (11.46) (gd⁴/ν)[(β/κ_s)dS/dz − (α/κ)dT̄/dz] = 27π⁴/4 (printed "657"): unstable when the left side is
    larger; statically stable when α dT̄/dz − β dS/dz > 0 (ρ̄ decreases upward).
    Parameters: dTdz [K/m]; dSdz [(g/kg)/m]; d [m]; alpha; beta_S; nu, kappa, kappa_s [m²/s]; g.
    Returns dict(unstable, lhs, margin (lhs − 27π⁴/4), density_stable, R_rho (= αT̄_z/(βS̄_z)), Ra (§11.5), Rs).
    Example (ours): dT/dz = 0.01, dS/dz = 0.002, d = 5 cm → lhs = 61 233.19 with the default g = G0 = 9.80665 (61 254.11
    with g = 9.81), density_stable, R_ρ = 1.316; thinnest finger-unstable layer ≈ 1.61 cm.  Validation: V1 sign flips at (11.46); V7 κ_s = κ → single-component threshold.
    Label: analytic.
    """
    lhs = g * d ** 4 / nu * (beta_S / kappa_s * dSdz - alpha / kappa * dTdz)  # Eq. (11.46) left side
    R_rho = (alpha * dTdz) / (beta_S * dSdz) if dSdz != 0 else math.inf
    return dict(unstable=bool(lhs > RA_FREE_FREE), lhs=float(lhs), margin=float(lhs - RA_FREE_FREE),
                density_stable=bool(alpha * dTdz - beta_S * dSdz > 0), R_rho=R_rho,
                Ra=float(thermal_rayleigh_signed(dTdz, d, alpha, nu, kappa, g)),
                Rs=float(salinity_rayleigh(dSdz, d, beta_S, nu, kappa_s, g)))


def double_diffusive_sigma(K2: float, Ra: float, Rs: float, Pr: float, tau: float, n: int = 1) -> np.ndarray:
    """The three growth rates σ (units κ/d²) of a free–free double-diffusive mode — our cubic (the book gives none).

    Book: §11.5, generalising (11.45) to σ ≠ 0, built from the dimensional equations (D13; design verifier note (i)):
    a² = n²π² + K², τ = κ_s/κ, (σ/Pr + a²)a²(σ + a²)(σ + τa²) = K²[−Ra(σ + τa²) + τRs(σ + a²)].  σ = 0 → Rs − Ra = a⁶/K².
    Parameters: K2; Ra (§11.5 sign); Rs (κ_s); Pr; tau; n.  Returns 3 complex roots sorted by descending Re σ.
    Examples (K² = π²/2, Pr = 7, τ = 1.5e-9/1.4e-7 = 0.0107143): (Ra, Rs) = (1000, 2000) → real root +0.03301 (fingers; the
    design's "+0.0308" does not satisfy this cubic — our root does to 1e-14); (−2×10⁴, −1.9×10⁶) → 9.650 ± 70.389i (diffusive,
    oscillatory; agrees with the design).
    Validation: V1 τ = 1, Rs = 0 → :func:`benard_free_free_sigma` with Ra_§11.4 = −Ra_§11.5 (11.0155 at K = π/√2, Ra = 2000,
    Pr = 1).  Label: analytic.
    """
    a2 = n ** 2 * math.pi ** 2 + K2
    lhs = np.poly1d([1.0 / Pr, a2]) * a2 * np.poly1d([1.0, a2]) * np.poly1d([1.0, tau * a2])
    rhs = K2 * (-Ra * np.poly1d([1.0, tau * a2]) + tau * Rs * np.poly1d([1.0, a2]))
    roots = np.roots((lhs - rhs).coeffs)
    return roots[np.argsort(-roots.real)]


def salt_finger_regime(dTdz, dSdz, d: float, alpha: float = 2e-4, beta_S: float = 7.6e-4, nu: float = 1e-6,
                       kappa: float = 1.4e-7, kappa_s: float = 1.5e-9, g: float = G0, Pr: float = 7.0,
                       K2: float | None = None) -> dict:
    """Classify a column: "overturning" (top-heavy *and* unstable), "fingers", "diffusive" or "stable".

    Book: §11.5, Fig. 11.13 ((a) hot salty over cold fresh → fingers, real σ; (b) cold fresh over hot salty → growing
    oscillations, complex σ) and p. 494: with T̂ = κ_sŝ/κ the set (11.45) is the Bénard set (11.39) with Ra → Rs − Ra, so
    the free–free layer is unstable only when (gd⁴/ν)[(β/κ_s)dS/dz − (α/κ)dT̄/dz] > 27π⁴/4 (11.46, printed "657").
    Classifier (ours), s = leading root of :func:`double_diffusive_sigma` at K² (default π²/2, the critical K² of (11.46)):
    * top-heavy or neutral (α dT̄/dz − β dS/dz ≤ 0): "overturning" only if the layer is actually unstable — margin
      (11.46) > 0 or Re s > 0; otherwise "stable" (viscosity and diffusion damp the top-heavy layer, exactly as a Bénard
      layer below its critical Rayleigh number).  Example (d = 5 cm, dS/dz = 0, defaults): −0.00751 < dT̄/dz ≤ 0 K/m is
      stable (Rs − Ra = 437.8 at −0.005, σ = −0.159), −0.0076 overturns (665.5, σ = +0.156).  (11.46) only locates
      the marginal state σ = 0 (where a real root changes sign); it is not a criterion for "no growth": a top-heavy layer
      with salt stabilising (dS/dz < 0) can have margin < 0 and still a growing root, labelled "overturning".  With
      margin < 0 the cubic is positive at σ = 0, so it has 0 or 2 positive real roots, and the growing root is either
      (i) one of a complex pair, e.g. (−0.015, −0.001) → σ = 3.296 ± 8.889i — the text says oscillatory; or (ii) real,
      once the pair has merged on the real axis, e.g. (−0.03, −0.0016) → margin −47 718, roots +16.890, +7.322, −142.8 —
      the text says monotonic overturning (a real root grows although σ = 0 of (11.46) is not crossed).  "Real" means
      |Im s| ≤ 1e-9·max(1, |s|) for the leading root s.
      "stable" therefore means: (11.46) not met *and* no growing root at this K² (other K² are not scanned).
    * bottom-heavy: Re s > 0 and s real ⇒ "fingers"; Re s > 0 and s complex ⇒ "diffusive"; else "stable".
    Returns dict(regime, R_rho, margin ((11.46) lhs − 27π⁴/4), sigma_max (complex, units κ/d²), density_stable, text).
    Validation: V7 four regimes; V1 the top-heavy label flips where margin (11.46) changes sign.  Label: analytic.
    """
    K2 = math.pi ** 2 / 2.0 if K2 is None else K2
    sf = salt_finger_unstable(dTdz, dSdz, d, alpha, beta_S, nu, kappa, kappa_s, g)
    s = complex(double_diffusive_sigma(K2, sf["Ra"], sf["Rs"], Pr, kappa_s / kappa)[0])
    grows = s.real > 1e-12
    is_real = abs(s.imag) <= 1e-9 * max(1.0, abs(s))
    if not sf["density_stable"]:
        if sf["margin"] > 0 or grows:  # Eq. (11.46): unstable only for Rs − Ra > 27π⁴/4 (or a growing mode at this K²)
            how = f"σ = {s.real:.4g} (real)" if is_real else f"σ = {s.real:.4g} ± {abs(s.imag):.4g}i"
            regime = "overturning"
            if sf["margin"] > 0:
                text = (f"density increases upward and Rs − Ra = {sf['lhs']:.4g} > 27π⁴/4 = 657.5: ordinary convective "
                        f"overturning, {how}")
            elif is_real:  # (11.46) not met, yet a real root grows: the complex pair has merged on the real axis
                text = (f"density increases upward; Rs − Ra = {sf['lhs']:.4g} < 27π⁴/4 = 657.5 (the marginal state σ = 0 "
                        f"of (11.46) is not crossed) yet a real root grows: monotonic overturning, {how}")
            else:  # no steady mode grows ((11.46) not met), but an oscillatory one does
                text = (f"density increases upward; Rs − Ra = {sf['lhs']:.4g} < 27π⁴/4 = 657.5 (no steady mode grows) but "
                        f"an oscillatory mode does: overturning by growing oscillations, {how}")
        else:
            regime = "stable"
            text = (f"density does not decrease upward, but Rs − Ra = {sf['lhs']:.4g} < 27π⁴/4 = 657.5: viscosity and "
                    f"diffusion win, no overturning (Re σ = {s.real:.4g})")
    elif grows and is_real:
        regime, text = "fingers", f"statically stable but finger-unstable: σ = {s.real:.4g} (real)"
    elif grows:
        regime, text = "diffusive", f"statically stable, oscillatory growth: σ = {s.real:.4g} ± {abs(s.imag):.4g}i"
    else:
        regime, text = "stable", "no growing free–free mode at this wavenumber"
    return dict(regime=regime, R_rho=sf["R_rho"], margin=sf["margin"], sigma_max=s, density_stable=sf["density_stable"],
                text=text)


# ======================================================================================================================
# §11.6 Taylor problem
# ======================================================================================================================
def ring_interchange_energy(Gamma1: float, Gamma2: float, r1: float, r2: float) -> dict:
    """Kinetic-energy change per unit mass when rings of equal mass at r₁ < r₂ swap places conserving Γ = 2πrU_θ.

    Book: §11.6, p. 496: E = U_θ²/2 = Γ²/(8π²r²); E_final = (Γ₂²/r₁² + Γ₁²/r₂²)/(8π²), E_initial = (Γ₁²/r₁² + Γ₂²/r₂²)/(8π²),
    ΔE = (Γ₂² − Γ₁²)(1/r₁² − 1/r₂²)/(8π²) (< 0: energy released ⇒ unstable, Rayleigh).
    Parameters: Gamma1, Gamma2 [m²/s]; r1 < r2 [m].  Returns dict(E_i, E_f, dE) [m²/s²].
    Example: Γ₁ = 4, Γ₂ = 2, r₁ = 1, r₂ = 2 → 0.21531, 0.10132, −0.11399.  Validation: V1 dE = E_f − E_i.  Label: analytic.
    """
    if not 0 < r1 < r2:
        raise ValueError("need 0 < r1 < r2")
    Ei = (Gamma1 ** 2 / r1 ** 2 + Gamma2 ** 2 / r2 ** 2) / (8 * math.pi ** 2)
    Ef = (Gamma2 ** 2 / r1 ** 2 + Gamma1 ** 2 / r2 ** 2) / (8 * math.pi ** 2)
    dE = (Gamma2 ** 2 - Gamma1 ** 2) * (1 / r1 ** 2 - 1 / r2 ** 2) / (8 * math.pi ** 2)  # p. 496
    return dict(E_i=Ei, E_f=Ef, dE=dE)


def rayleigh_circulation_criterion(r, U_theta) -> dict:
    """Rayleigh's inviscid criterion for swirling flow: unstable iff dΓ²/dr < 0 somewhere (Γ = 2πrU_θ).

    Book: §11.6, p. 497 ("An inviscid Couette flow is unstable if dΓ²/dr < 0"; analogue of dρ̄/dz > 0).
    Parameters: r [m] (increasing); U_theta (array on r or callable) [m/s].
    Returns dict(unstable, where (midpoints with dΓ²/dr < 0), dGamma2_dr (finite differences on the midpoints)).
    Validation: V1 solid body stable; outer-at-rest Couette unstable; flips at Ω₂/Ω₁ = (R₁/R₂)² (Ch. 8 parity).  Label: analytic.
    """
    r = _F(r)
    U = _F(U_theta(r) if callable(U_theta) else U_theta)
    G2 = (2 * math.pi * r * U) ** 2
    dG2 = np.diff(G2) / np.diff(r)
    scale = max(np.max(np.abs(G2)), 1e-300) / max(np.ptp(r), 1e-300)
    bad = dG2 < -1e-12 * scale
    return dict(unstable=bool(bad.any()), where=0.5 * (r[1:] + r[:-1])[bad], dGamma2_dr=dG2)


def couette_rayleigh_line(R1: float, R2: float) -> float:
    """μ = Ω₂/Ω₁ on Rayleigh's neutral line, μ_R = (R₁/R₂)² (Ω₁/Ω₂ = R₂²/R₁², Fig. 11.17 dashed line).

    Book: §11.6 (Ch. 8: Rayleigh-stable iff Ω₂/Ω₁ > (R₁/R₂)²).  Example: R₂/R₁ = 1.05 → 0.9070.  Label: analytic.
    """
    return (R1 / R2) ** 2


def taylor_number(Omega1: float, Omega2: float, R1: float, R2: float, nu: float) -> dict:
    """Taylor number Ta = 4[(Ω₁R₁² − Ω₂R₂²)/(R₂² − R₁²)]Ω₁d⁴/ν², d = R₂ − R₁, with Rayleigh's flag.

    Book: §11.6, Eq. (11.52) (= −4AΩ₁d⁴/ν² with A of (11.49); ≤ 0 on/beyond the Rayleigh line — design verifier note iii).
    Parameters: Omega1, Omega2 [rad/s]; R1 < R2 [m]; nu [m²/s].
    Returns dict(Ta, rayleigh_stable (μ ≥ (R₁/R₂)² for Ω₁ > 0), mu = Ω₂/Ω₁).
    Example: R₁ = 0.1, R₂ = 0.105 m, Ω₁ = 0.5 rad/s, Ω₂ = 0, ν = 1e-6 → Ta = 6098.  Validation: V1 narrow inner limit;
    V2 dimensionless.  Label: analytic.
    """
    d = R2 - R1
    if d <= 0:
        raise ValueError("R2 > R1 required")
    Ta = 4.0 * (Omega1 * R1 ** 2 - Omega2 * R2 ** 2) / (R2 ** 2 - R1 ** 2) * Omega1 * d ** 4 / nu ** 2  # Eq. (11.52)
    mu = Omega2 / Omega1 if Omega1 != 0 else math.inf
    return dict(Ta=Ta, rayleigh_stable=bool(mu >= couette_rayleigh_line(R1, R2)), mu=mu)


def taylor_number_narrow_inner(Omega1: float, R1: float, d: float, nu: float) -> float:
    """Narrow-gap, inner-cylinder-only Ta = 2(Ω₁R₁d/ν)²(d/R₁) (§11.6, below (11.52)).  Example: 6250 for the case above.
    Label: analytic.
    """
    return 2.0 * (Omega1 * R1 * d / nu) ** 2 * (d / R1)


def _taylor_ops(k: float, N: int):
    D, x = cheb(N, (0.0, 1.0))
    n = N + 1
    I = np.eye(n)
    return D, x, D @ D - k ** 2 * I, I, n


def _taylor_constraints(D, I, n, N):
    Zr = np.zeros(n)
    rows = [np.hstack([I[0], Zr]), np.hstack([I[N], Zr]), np.hstack([D[0], Zr]), np.hstack([D[N], Zr]),
            np.hstack([Zr, I[0]]), np.hstack([Zr, I[N]])]  # Eq. (11.53): û_R = dû_R/dx = û_φ = 0 at x = 0, 1
    return np.array(rows), [0, N, 1, N - 1, n, n + N]


def taylor_marginal_Ta(k: float, mu: float, N: int = 40, return_mode: bool = False):
    """Marginal Taylor number Ta(k) of the narrow-gap problem (σ = 0): smallest positive real eigenvalue.

    Book: §11.6, Eq. (11.51) at σ = 0 — (d²/dx² − k²)²û_R = (1 + αx)û_φ, (d²/dx² − k²)û_φ = −Tak²û_R — with (11.53);
    α ≡ Ω₂/Ω₁ − 1 = μ − 1, x = (R − R₁)/d ∈ [0, 1], k ≡ kd ("d/dR" in (11.51) means d/dx).
    Parameters: k > 0; mu; N (default 40; FAST 24); return_mode (also dict(x, uR, uphi), û_R max-normalised).
    Returns Ta (NaN if no positive real eigenvalue).  Example: μ = 0, k = 3.1266 → 3389.90.
    Validation: V1 μ = 1 = rigid Bénard 1707.76 at 3.117; V3; V1 (independent route) :func:`taylor_galerkin_Ta`.
    Label: converged.
    """
    if k <= 0:
        raise ValueError("k must be > 0")
    D, x, L, I, n = _taylor_ops(k, N)
    Zm = np.zeros_like(I)
    A = np.block([[L @ L, -np.diag(1.0 + (mu - 1.0) * x)], [Zm, L]])  # Eq. (11.51), σ = 0
    B = np.block([[Zm, Zm], [-k ** 2 * I, Zm]])
    C, elim = _taylor_constraints(D, I, n, N)
    w, V = constrained_eig(A, B, C, elim, return_vectors=True, sort="abs")
    ok = (np.abs(w.imag) <= 1e-8 * np.abs(w)) & (w.real > 0)
    if not ok.any():
        return (float("nan"), None) if return_mode else float("nan")
    idx = np.nonzero(ok)[0]
    j = idx[np.argmin(w.real[idx])]
    Ta = float(w.real[j])
    if not return_mode:
        return Ta
    s = V[:n, j][np.argmax(np.abs(V[:n, j]))]
    return Ta, dict(x=x, uR=np.real(V[:n, j] / s), uphi=np.real(V[n:, j] / s))


def taylor_critical(mu: float, N: int = 40, k_bounds: Sequence[float] = (2.0, 4.5), xatol: float = 1e-7) -> dict:
    """Critical (Ta_c, k_c) of the narrow-gap problem: minimum over k of :func:`taylor_marginal_Ta`.

    Book: §11.6, (11.51)–(11.54); k_cr ≈ 3.12, λ ≈ 2d.  Returns dict(Ta_c, k_c).
    Expect μ = 1: 1707.76 at 3.1163; 0.5: 2275.09 at 3.1175; 0: 3389.90 at 3.1266; −0.5: 6413.7 at 3.1985.  Label: converged.
    """
    ks = np.linspace(k_bounds[0], k_bounds[1], 9)
    vals = np.array([taylor_marginal_Ta(k, mu, N) for k in ks])
    if not np.isfinite(vals).any():
        raise RuntimeError("taylor_critical: no stationary marginal state found (μ too negative?)")
    i = int(np.nanargmin(vals))
    r = minimize_scalar(lambda k: taylor_marginal_Ta(k, mu, N), bounds=(ks[max(i - 1, 0)], ks[min(i + 1, 8)]),
                        method="bounded", options=dict(xatol=xatol))
    return dict(Ta_c=float(r.fun), k_c=float(r.x))


def taylor_critical_approx(mu: float) -> float:
    """Ta_cr = 1708/[½(1 + Ω₂/Ω₁)] — Eq. (11.54) (approximate: +0.77 % at μ = 0, +6.5 % at μ = −0.5).  Label: analytic."""
    if mu <= -1:
        raise ValueError("(11.54) needs Ω₂/Ω₁ > −1")
    return 1708.0 / (0.5 * (1.0 + mu))  # Eq. (11.54)


def taylor_critical_table(mus=None, N: int = 40, cache: bool = True) -> dict:
    """Exact narrow-gap (Ta_c, k_c) against (11.54) on a μ grid (default −0.5 … 1.0 step 0.05; cached; ours).

    Book: §11.6.  Returns dict(mu, Ta_c, k_c, approx, rel_error (approx/exact − 1)).  Label: converged.
    """
    mus = [float(v) for v in (np.round(np.arange(-0.5, 1.0001, 0.05), 10) if mus is None else mus)]

    def compute():
        res = [taylor_critical(m, N) for m in mus]
        Ta = np.array([r["Ta_c"] for r in res])
        appr = np.array([taylor_critical_approx(m) for m in mus])
        return dict(mu=np.array(mus), Ta_c=Ta, k_c=np.array([r["k_c"] for r in res]), approx=appr, rel_error=appr / Ta - 1.0)

    return _cached("taylor_critical_table", dict(mus=mus, N=N), compute, cache)


def taylor_growth_rate(k: float, Ta: float, mu: float, N: int = 40, all: bool = False):
    """Leading growth rate σ (units ν/d², inferred) of narrow-gap Taylor modes (complex; real for μ ≥ 0).

    Book: §11.6, Eq. (11.51) (d²/dx² − k² − σ)(d²/dx² − k²)û_R = (1 + αx)û_φ, (d²/dx² − k² − σ)û_φ = −Tak²û_R, (11.53):
    L²û_R − (1 + αx)û_φ = σLû_R, Lû_φ + Tak²û_R = σû_φ.  Returns complex σ (or all, sorted by Re σ).
    Validation: σ = 0 at :func:`taylor_marginal_Ta`.  Label: converged.
    """
    D, x, L, I, n = _taylor_ops(k, N)
    Zm = np.zeros_like(I)
    A = np.block([[L @ L, -np.diag(1.0 + (mu - 1.0) * x)], [Ta * k ** 2 * I, L]])  # Eq. (11.51)
    B = np.block([[L, Zm], [Zm, I]])
    C, elim = _taylor_constraints(D, I, n, N)
    w = constrained_eig(A, B, C, elim, sort="real")
    return w if all else complex(w[0])


def taylor_eigenfunction(k: float, mu: float, x=None, z=None, N: int = 40) -> dict:
    """Marginal narrow-gap Taylor vortices: û_R(x), û_φ(x) and the meridional fields over one axial wavelength.

    Book: §11.6, (11.51), Fig. 11.16.  u_R = û_R(x)cos kz, u_φ = û_φ(x)cos kz, u_z = −û_R′(x) sin(kz)/k (continuity in the
    narrow gap), stream function ψ = û_R(x) sin(kz)/k with u_R = ∂ψ/∂z, u_z = −∂ψ/∂x (convention stated here).
    Defaults: x ∈ [0, 1] (41), z over 2π/k (41).  û_R max-normalised.
    Returns dict(x, z (2-D), u_R, u_phi, u_z, psi, Ta, uR_profile, uphi_profile, x_nodes).  Label: converged.
    """
    Ta, m = taylor_marginal_Ta(k, mu, N, return_mode=True)
    x = np.linspace(0.0, 1.0, 41) if x is None else _F(x)
    z = np.linspace(0.0, 2.0 * math.pi / k, 41) if z is None else _F(z)
    X, Z = np.meshgrid(x, z, indexing="xy")
    uR = BarycentricInterpolator(m["x"], m["uR"])(X)
    uphi = BarycentricInterpolator(m["x"], m["uphi"])(X)
    duR = BarycentricInterpolator(m["x"], cheb(N, (0.0, 1.0))[0] @ m["uR"])(X)
    return dict(x=X, z=Z, u_R=uR * np.cos(k * Z), u_phi=uphi * np.cos(k * Z), u_z=-duR * np.sin(k * Z) / k,
                psi=uR * np.sin(k * Z) / k, Ta=Ta, uR_profile=m["uR"], uphi_profile=m["uphi"], x_nodes=m["x"])


def taylor_galerkin_Ta(k: float, mu: float, n_modes: int = 4, printed: bool = False, N: int = 48) -> float:
    """Marginal Ta(k) by the Galerkin route of Exercise 11.9 (an independent discretisation of (11.51)).

    Book: Ex. 11.9, Eqs. (11.92)–(11.94): û_φ = Σ_{m=1}^{M} C_m sin(mπx) meets û_φ(0) = û_φ(1) = 0; for each m solve (11.92)
    (d²/dx² − k²)²û_R^{(m)} = (1 + αx) sin(mπx) with û_R = û_R′ = 0 (Chebyshev BVP, N); project (11.93) on sin(nπx):
    ½(n²π² + k²)C_n = Tak² Σ_m M_nm C_m, M_nm = ∫₀¹ sin(nπx)û_R^{(m)}dx.  ``printed=True`` (slip S7) uses the printed squared
    operator: ½(n²π² + k²)²C_n = −Tak²(MC)_n.
    Returns the smallest positive real Ta (NaN if none — the printed form gives none).
    Validation: V1 agrees with :func:`taylor_marginal_Ta` to 1e-4 at M = 4 (3390.27 vs 3389.90 at μ = 0).  Label: converged.
    """
    D, x, L, I, n = _taylor_ops(k, N)
    al = mu - 1.0
    w_cc = clenshaw_curtis_weights(N, (0.0, 1.0))
    L2 = (L @ L).astype(float)
    elim = [0, N, 1, N - 1]
    keep = np.setdiff1d(np.arange(n), elim)
    C = np.vstack([I[0], I[N], D[0], D[N]])
    X = -np.linalg.solve(C[:, elim], C[:, keep])
    T = np.zeros((n, len(keep)))
    T[keep, np.arange(len(keep))] = 1.0
    T[elim] = X
    Ared = L2[keep] @ T
    M = np.zeros((n_modes, n_modes))
    for m in range(1, n_modes + 1):
        uR = T @ np.linalg.solve(Ared, ((1.0 + al * x) * np.sin(m * math.pi * x))[keep])  # Eq. (11.92)
        for nn in range(1, n_modes + 1):
            M[nn - 1, m - 1] = np.sum(w_cc * np.sin(nn * math.pi * x) * uR)
    lam = np.array([nn ** 2 * math.pi ** 2 + k ** 2 for nn in range(1, n_modes + 1)])
    Dg = np.diag(0.5 * (lam ** 2 if printed else lam))
    ev = sla.eig(Dg, (-1.0 if printed else 1.0) * k ** 2 * M, right=False)  # Eq. (11.93) projected
    ev = ev[np.isfinite(ev)]
    ok = (np.abs(ev.imag) <= 1e-8 * np.abs(ev)) & (ev.real > 0)
    return float(np.min(ev.real[ok])) if ok.any() else float("nan")


def taylor_stability_boundary(R2_over_R1: float = 1.05, mus=None, N: int = 40) -> dict:
    """Narrow-gap marginal curve in the (Ω₂R₂²/ν, Ω₁R₂²/ν) plane of Fig. 11.17 for our radius ratio (default 1.05).

    Book: §11.6, Fig. 11.17 (Taylor's experiment and theory; the book's ratio is private).  From (11.52) with Ω₂ = μΩ₁:
    Ta = 4Ω₁²(R₁² − μR₂²)d⁴/((R₂² − R₁²)ν²) ⇒ Ω₁R₂²/ν = (R₂²/d²)[Ta_c(μ)(R₂² − R₁²)/(4(R₁² − μR₂²))]^{1/2}
    (scale-free: depends on R₂/R₁ only); valid for μ < (R₁/R₂)² (the Rayleigh-unstable side).
    Parameters: R2_over_R1 > 1; mus (default 25 values from −1 to 0.95(R₁/R₂)²); N.
    Returns dict(x_outer = Ω₂R₂²/ν, y_inner = Ω₁R₂²/ν, rayleigh_x, rayleigh_y (the line Ω₁/Ω₂ = R₂²/R₁²), mu, Ta_c).
    Assumptions: narrow-gap theory at a finite gap (approximate).  Label: converged.
    """
    R1, R2 = 1.0, float(R2_over_R1)
    d = R2 - R1
    muR = couette_rayleigh_line(R1, R2)
    mus = np.linspace(-1.0, 0.95 * muR, 25) if mus is None else _F(mus)
    xs, ys, mm, tc = [], [], [], []
    for m in mus:
        if m >= muR:
            continue
        try:
            Ta_c = taylor_critical(float(m), N)["Ta_c"]
        except RuntimeError:
            continue
        Y = (R2 ** 2 / d ** 2) * math.sqrt(Ta_c * (R2 ** 2 - R1 ** 2) / (4.0 * (R1 ** 2 - m * R2 ** 2)))
        xs.append(m * Y)
        ys.append(Y)
        mm.append(m)
        tc.append(Ta_c)
    xs = np.array(xs)
    rx = np.linspace(0.0, max(1.0, float(np.max(np.abs(xs))) if xs.size else 1.0), 50)
    return dict(x_outer=xs, y_inner=np.array(ys), rayleigh_x=rx, rayleigh_y=rx * R2 ** 2 / R1 ** 2, mu=np.array(mm),
                Ta_c=np.array(tc))


def taylor_perturbation_sympy() -> dict:
    """D14 (linearisation and narrow-gap steps) and slip S12.

    Book: §11.6, Eqs. (11.47)–(11.52).  Substituting ũ = U + εu, p̃ = P + εp with U_φ = AR + B/R, (1/ρ)dP/dR = U_φ²/R (11.49)
    into the axisymmetric (11.47): O(1) vanishes (base state), O(ε) is (11.50).  dU_φ/dR + U_φ/R = 2A (the factor of (11.50));
    the operator identity d/dR[(1/R)d(Rf)/dR] = f″ + f′/R − f/R²; Ta of (11.52) equals −4AΩ₁d⁴/ν².
    S12: the printed continuity ∂(Rũ_R)/∂R + ∂ũ_z/∂z = 0 mixes units (lengths scaled by λ: terms scale as λ⁰ and λ⁻¹); the
    correct (1/R)∂(Rũ_R)/∂R + ∂ũ_z/∂z = 0 is homogeneous and vanishes for a Stokes-stream-function field.
    Returns dict(lin_R, lin_phi, lin_z, base_R, base_phi (all 0), term_2A (0), operator_identity (0), narrow_gap_eqs (the two
    (11.51) expressions), Ta_coefficient (−4AΩ₁d⁴/ν²), Ta_identity (0), continuity_correct (0), continuity_printed (≠ 0),
    printed_continuity_units_ok (False), correct_continuity_units_ok (True)).  Validation: V2.  Label: symbolic.
    """
    R, z, t = sp.symbols("R z t", positive=True)
    A, B, rho, nu, eps = sp.symbols("A B rho nu epsilon")
    V = A * R + B / R  # Eq. (11.49)
    P = sp.Function("P")(R)
    uR, uphi, uz, p = [sp.Function(nm)(R, z, t) for nm in ("u_R", "u_phi", "u_z", "p")]
    tR, tphi, tz, tp = eps * uR, V + eps * uphi, eps * uz, P + eps * p
    DDt = lambda f: sp.diff(f, t) + tR * sp.diff(f, R) + tz * sp.diff(f, z)  # noqa: E731
    lap = lambda f: sp.diff(f, R, 2) + sp.diff(f, R) / R + sp.diff(f, z, 2)  # noqa: E731
    eR = DDt(tR) - tphi ** 2 / R + sp.diff(tp, R) / rho - nu * (lap(tR) - tR / R ** 2)  # Eq. (11.47)
    ephi = DDt(tphi) + tR * tphi / R - nu * (lap(tphi) - tphi / R ** 2)
    ez = DDt(tz) + sp.diff(tp, z) / rho - nu * lap(tz)
    sub = {sp.Derivative(P, R): rho * V ** 2 / R}
    e50R = sp.diff(uR, t) - 2 * V * uphi / R + sp.diff(p, R) / rho - nu * (lap(uR) - uR / R ** 2)  # Eq. (11.50)
    e50phi = sp.diff(uphi, t) + (sp.diff(V, R) + V / R) * uR - nu * (lap(uphi) - uphi / R ** 2)
    e50z = sp.diff(uz, t) + sp.diff(p, z) / rho - nu * lap(uz)
    f = sp.Function("f")(R)
    op_id = sp.simplify(sp.diff(sp.diff(R * f, R) / R, R) - (f.diff(R, 2) + f.diff(R) / R - f / R ** 2))
    Om1, Om2, R1, R2, nuv = sp.symbols("Omega1 Omega2 R1 R2 nu_v", positive=True)
    A_ = (Om2 * R2 ** 2 - Om1 * R1 ** 2) / (R2 ** 2 - R1 ** 2)  # (11.49)
    d = R2 - R1
    Ta52 = 4 * (Om1 * R1 ** 2 - Om2 * R2 ** 2) / (R2 ** 2 - R1 ** 2) * Om1 * d ** 4 / nuv ** 2  # Eq. (11.52)
    Ta_coef = -4 * A * Om1 * d ** 4 / nuv ** 2
    Psi = sp.Function("Psi")(R, z)
    ur_s, uz_s = -sp.diff(Psi, z) / R, sp.diff(Psi, R) / R
    Lu, Vu = sp.symbols("L_unit V_unit", positive=True)  # dimensions: R, z ~ L; u_R, u_z ~ V; ∂/∂R, ∂/∂z ~ 1/L
    dims_pr = [(Lu * Vu) / Lu, Vu / Lu]  # ∂(Rũ_R)/∂R, ∂ũ_z/∂z  (printed)
    dims_ok = [(Lu * Vu) / Lu / Lu, Vu / Lu]  # (1/R)∂(Rũ_R)/∂R, ∂ũ_z/∂z
    units_pr = sp.simplify(dims_pr[0] / dims_pr[1]) == 1
    units_ok = sp.simplify(dims_ok[0] / dims_ok[1]) == 1
    x_, k_, sg, al, Ta = sp.symbols("x k sigma alpha Ta")
    uRf, uphif = sp.Function("u_R")(x_), sp.Function("u_phi")(x_)
    Lop = lambda f: f.diff(x_, 2) - k_ ** 2 * f  # noqa: E731
    narrow = [sp.Eq(Lop(Lop(uRf)) - sg * Lop(uRf), (1 + al * x_) * uphif),  # Eq. (11.51)
              sp.Eq(Lop(uphif) - sg * uphif, -Ta * k_ ** 2 * uRf)]
    return dict(base_R=sp.simplify(eR.subs(eps, 0).subs(sub)), base_phi=sp.simplify(ephi.subs(eps, 0)),
                lin_R=sp.simplify(sp.diff(eR, eps).subs(eps, 0) - e50R), lin_phi=sp.simplify(sp.diff(ephi, eps).subs(eps, 0) - e50phi),
                lin_z=sp.simplify(sp.diff(ez, eps).subs(eps, 0) - e50z),
                term_2A=sp.simplify(sp.diff(V, R) + V / R - 2 * A), operator_identity=op_id, narrow_gap_eqs=narrow,
                Ta_coefficient=Ta_coef, Ta_identity=sp.simplify(Ta52 - Ta_coef.subs(A, A_)),
                continuity_correct=sp.simplify(sp.diff(R * ur_s, R) / R + sp.diff(uz_s, z)),
                continuity_printed=sp.simplify(sp.diff(R * ur_s, R) + sp.diff(uz_s, z)),
                printed_continuity_units_ok=bool(units_pr), correct_continuity_units_ok=bool(units_ok))


# ======================================================================================================================
# §11.7 stratified parallel flows
# ======================================================================================================================
def stratified_shear_sympy(printed: bool = False) -> dict:
    """D16–D19: (11.58)–(11.60) from (11.57); the Taylor–Goldstein equation (11.61); (11.61) ⇔ (11.64) under (11.63);
    (11.61) ⇔ the divergence form d/dz[(U − c)²F′] − k²(U − c)²F + N²F = 0 under (11.68).

    Book: §11.7, p. 503–506.  Slips: S4 — (11.55)'s w-equation prints −(1/ρ₀)∂p/∂x: with it p̂ cannot be eliminated
    (``printed_slip4_reaches_tg`` is False); S11 — the divergence form prints −k²(U − c)F (``printed_slip11_residual`` ≠ 0).
    ``printed=True`` puts the printed versions into ``residual``.
    Returns dict(tg_residual (0), self_adjoint_residual (0), divergence_residual (0), printed_slip11_residual (≠ 0),
    printed_slip4_reaches_tg (False), normal_modes ([0, 0, 0]), printed_slip4_residual, residual).  Validation: V2.
    Label: symbolic.
    """
    x, z, t, k, g, rho0 = sp.symbols("x z t k g rho0", positive=True)
    c = sp.symbols("c")
    U, N2 = sp.Function("U")(z), sp.Function("N2")(z)
    ph, pp, rh = sp.Function("psi")(z), sp.Function("p")(z), sp.Function("rho")(z)
    E = sp.exp(sp.I * k * (x - c * t))
    psi, p, rho = ph * E, pp * E, rh * E
    e1 = sp.diff(psi, t, z) - sp.diff(psi, x) * sp.diff(U, z) + U * sp.diff(psi, x, z) + sp.diff(p, x) / rho0  # (11.57)
    e2_ok = -sp.diff(psi, t, x) - U * sp.diff(psi, x, 2) + g * rho / rho0 + sp.diff(p, z) / rho0
    e2_pr = -sp.diff(psi, t, x) - U * sp.diff(psi, x, 2) + g * rho / rho0 + sp.diff(p, x) / rho0  # slip S4
    e3 = sp.diff(rho, t) + U * sp.diff(rho, x) + rho0 * N2 / g * sp.diff(psi, x)
    n58 = sp.simplify(sp.simplify(e1 / (sp.I * k * E)) - ((U - c) * ph.diff(z) - U.diff(z) * ph + pp / rho0))  # (11.58)
    n59 = sp.simplify(sp.simplify(e2_ok / E) - (k ** 2 * (U - c) * ph + g * rh / rho0 + pp.diff(z) / rho0))  # (11.59)
    n60 = sp.simplify(sp.simplify(e3 / (sp.I * k * E)) - ((U - c) * rh + rho0 * N2 * ph / g))  # (11.60)
    tg = (U - c) * (ph.diff(z, 2) - k ** 2 * ph) - U.diff(z, 2) * ph + N2 * ph / (U - c)  # Eq. (11.61)

    def tg_from(e2x):
        a2 = sp.simplify(e2x / E)
        rho_sol = -rho0 * N2 * ph / (g * (U - c))  # (11.60)
        p_sol = -rho0 * ((U - c) * ph.diff(z) - U.diff(z) * ph)  # (11.58)
        return sp.simplify(a2.subs(rh, rho_sol).subs(pp, p_sol).doit() + tg)

    tg_ok, tg_pr = tg_from(e2_ok), tg_from(e2_pr)
    phi = sp.Function("phi")(z)
    s = sp.sqrt(U - c)
    psi_s = s * phi  # Eq. (11.63)
    tg_phi = (U - c) * (psi_s.diff(z, 2) - k ** 2 * psi_s) - U.diff(z, 2) * psi_s + N2 * psi_s / (U - c)
    f64 = (sp.diff((U - c) * phi.diff(z), z)
           - (k ** 2 * (U - c) + U.diff(z, 2) / 2 + (U.diff(z) ** 2 / 4 - N2) / (U - c)) * phi)  # Eq. (11.64)
    sa = sp.simplify(sp.expand(tg_phi / s - f64))
    F = sp.Function("F")(z)
    psi_F = (U - c) * F  # Eq. (11.68)
    tg_F = (U - c) * (psi_F.diff(z, 2) - k ** 2 * psi_F) - U.diff(z, 2) * psi_F + N2 * psi_F / (U - c)
    div_ok = sp.diff((U - c) ** 2 * F.diff(z), z) - k ** 2 * (U - c) ** 2 * F + N2 * F
    div_pr = sp.diff((U - c) ** 2 * F.diff(z), z) - k ** 2 * (U - c) * F + N2 * F  # slip S11
    d_ok = sp.simplify(sp.expand(tg_F - div_ok))
    d_pr = sp.simplify(sp.expand(tg_F - div_pr))
    return dict(tg_residual=tg_ok, self_adjoint_residual=sa, divergence_residual=d_ok, printed_slip11_residual=d_pr,
                printed_slip4_reaches_tg=bool(tg_pr == 0), printed_slip4_residual=tg_pr, normal_modes=[n58, n59, n60],
                residual=(tg_pr, d_pr) if printed else (tg_ok, d_ok), printed=printed)


def gradient_richardson(z, U=None, N2=None, dUdz=None):
    """Gradient Richardson number Ri(z) = N²/(dU/dz)²; where dU/dz = 0 it carries the sign of N².

    Book: §11.7, Eq. (11.66); N² = −(g/ρ₀)dρ̄/dz (7.127).  Parameters: z [m]; U, N2 (callables or arrays) [m/s, 1/s²]; dUdz
    (callable/array; default central difference for callables, ``np.gradient(edge_order=2)`` for arrays).
    Returns Ri [–].  At a shear-free level (dU/dz = 0) the limit of N²/(dU/dz)² is taken: +inf for N² > 0, −inf for N² < 0
    (statically unstable — never "Ri > ¼"), nan for N² = 0 (0/0: undefined; :func:`miles_howard_stable` does not count it as
    satisfying (11.67)).
    Validation: V1 U = tanh z, N² = J sech²z → Ri = J cosh²z; V7 the three dU/dz = 0 limits.  Label: analytic.
    """
    _, N2v, Up = _richardson_inputs(z, U, N2, dUdz)
    return _S(_richardson_ratio(N2v, Up))


def _richardson_inputs(z, U, N2, dUdz):
    """(z, N², dU/dz) as float arrays (shared by the two Richardson functions; not broadcast against each other)."""
    z = _F(z)
    if N2 is None:
        raise ValueError("N2 is required")
    N2v = _F(N2(z) if callable(N2) else N2)
    if dUdz is not None:
        Up = _F(dUdz(z) if callable(dUdz) else dUdz)
    elif callable(U):
        h = 1e-6 * max(1.0, float(np.max(np.abs(z))))
        Up = (U(z + h) - U(z - h)) / (2 * h)
    elif U is not None:
        Up = np.gradient(_F(U), z, edge_order=2)
    else:
        raise ValueError("give U or dUdz")
    return z, N2v, Up


def _richardson_ratio(N2v, Up):
    with np.errstate(divide="ignore", invalid="ignore"):
        # Eq. (11.66); at dU/dz = 0: sign(N²)·inf (np.sign(0)·inf = nan — the undefined 0/0 case)
        return np.where(Up == 0, np.sign(N2v) * np.inf, N2v / Up ** 2)


def miles_howard_stable(z, U=None, N2=None, dUdz=None) -> dict:
    """Miles–Howard: Ri > ¼ everywhere ⇒ linearly stable (Ri < ¼ somewhere is only *necessary* for instability).

    Book: §11.7, Eqs. (11.65)–(11.67): c_i = 0 is forced when N² > ¼(dU/dz)² everywhere, i.e. Ri ≡ N²/(dU/dz)² > ¼.
    The verdict tests the book's inequality N² > ¼(dU/dz)² pointwise (the form before the division), so a shear-free level
    is judged by the sign of N²: N² > 0 passes (Ri = +inf), N² < 0 fails (Ri = −inf, convectively unstable), and
    N² = 0 = dU/dz fails too (0 > 0 is false; Ri is nan there).
    Parameters as :func:`gradient_richardson`.
    Returns dict(guaranteed_stable, Ri_min (smallest defined Ri; nan if Ri is undefined everywhere), z_min (its level, or
    the first level where the inequality fails if every defined Ri exceeds ¼), text).  Label: analytic.
    """
    zz, N2v, Up = (np.atleast_1d(a) for a in _richardson_inputs(z, U, N2, dUdz))
    Ri = _richardson_ratio(N2v, Up)
    if zz.shape != Ri.shape:  # constant N² and dU/dz given as numbers with an array z (or the reverse)
        zz, N2v, Up, Ri = np.broadcast_arrays(zz, N2v, Up, Ri)
    holds = N2v > 0.25 * Up ** 2  # Eq. (11.67) before dividing by (dU/dz)²: N² > ¼(dU/dz)²
    ok = bool(np.all(holds))
    defined = ~np.isnan(Ri)
    j = int(np.argmin(np.where(defined, Ri, np.inf))) if defined.any() else int(np.argmin(holds))
    if not ok and holds[j]:  # every defined Ri exceeds ¼: the failure is an undefined (N² = 0 = dU/dz) level
        j = int(np.argmin(holds))
    Ri_min, z_min = float(Ri[j]), float(zz[j])
    if ok:
        text = f"Ri > 1/4 everywhere (Ri_min = {Ri_min:.4g}): stable by Miles–Howard"
    elif math.isnan(Ri_min):
        text = (f"N² = 0 and dU/dz = 0 at z = {z_min:.4g}: Ri is undefined there and N² > (1/4)(dU/dz)² does not hold — "
                "stability not guaranteed")
    elif Ri_min < 0:
        text = (f"Ri_min = {Ri_min:.4g} < 0 at z = {z_min:.4g} (N² < 0: statically unstable there) — not Ri > 1/4, "
                "stability not guaranteed")
    elif Ri_min == 0.25:
        text = f"Ri_min = 0.25 at z = {z_min:.4g}: not > 1/4, instability allowed, not guaranteed"
    else:
        text = f"Ri_min = {Ri_min:.4g} < 1/4 at z = {z_min:.4g}: instability allowed, not guaranteed"
    return dict(guaranteed_stable=ok, Ri_min=Ri_min, z_min=z_min, text=text)


def richardson_profiles(kind: str = "tanh", J: float = 0.1, R: float = 1.0) -> dict:
    """Stratified shear-layer family (callables): U = tanh z, N² = J sech²(Rz) (J = Ri(0); R = shear/density thickness ratio).

    Book: §11.7, p. 505–506 (tanh shear layers are unstable iff Ri_min < ¼).  R = 1: Ri(z) = J cosh²z, Ri_min = J; exact neutral
    curve J = k(1 − k) (:func:`tg_tanh_neutral_J`).  Returns dict(U, Up, Upp, N2, Ri, kind, J, R).  Label: analytic.
    """
    if kind != "tanh":
        raise ValueError('only kind="tanh" is implemented')
    Up = lambda z: 1.0 / np.cosh(z) ** 2  # noqa: E731
    N2 = lambda z: J / np.cosh(np.clip(R * _F(z), -350.0, 350.0)) ** 2  # noqa: E731  (clip: cosh² overflows past |z| ≈ 355)
    return dict(U=np.tanh, Up=Up, Upp=lambda z: -2.0 * np.tanh(z) / np.cosh(np.clip(z, -350.0, 350.0)) ** 2, N2=N2,
                Ri=lambda z: N2(z) / Up(z) ** 2, kind=kind, J=J, R=R)


def tg_tanh_neutral_J(k):
    """Exact neutral stratification for U = tanh z, N² = J sech²z: J = k(1 − k), 0 < k < 1 (maximum ¼ at k = ½).

    Book: §11.7 (Miles–Howard's ¼ is attained).  Our check: ψ̂ = |tanh z|^{1−k} sech^k z with c = 0 satisfies (11.61) exactly
    at J = k(1 − k) (numerical residual ~1e-16 at k = 0.3, 0.6) — the classical Drazin/Hazel result.  Label: analytic.
    """
    k = _F(k)
    return _S(k * (1.0 - k))


def tg_tanh_neutral_mode(z, k: float):
    """ψ̂ = |tanh z|^{1−k} sech^k z — the exact neutral mode (c = 0, J = k(1 − k)) of :func:`tg_tanh_neutral_J`.
    Label: analytic.
    """
    z = _F(z)
    return _S(np.abs(np.tanh(z)) ** (1.0 - k) / np.cosh(z) ** k)


def tg_growth(k: float, J: float, R: float = 1.0, N: int = 100, map_scale: float | None = None,
              y_max: float | None = None) -> float:
    """Leading growth rate kc_i of the tanh / J sech²(Rz) layer (Taylor–Goldstein, decaying BCs; 0 if no converged mode).

    Book: §11.7 (11.61).
    Parameters: k > 0 wavenumber (scaled by 1/L, L the shear-layer half-thickness); J = Ri(0) [–]; R thickness ratio [–];
    N Chebyshev degree (default 100; the N-filter's second solve uses 150); map_scale tan-map scale s (scaled by L; default
    None → the rule :func:`decay_map_scale`(k) = min(0.5, max(0.035, 0.025/k)); a number forces that scale); y_max half-width
    of the box where ψ̂ = 0 is imposed (scaled by L).
    DEVIATION (node distribution; Part C wrote s = 3.0, loops 0–1 used a fixed 0.5): near the neutral curve the mode has a
    critical layer at z = 0 that a fixed s = 0.5 does not resolve at N = 100 for k ≥ 0.65 — the filter then dropped modes with
    kc_i up to 0.033 and reported 0 in a strip 0.04–0.06 wide under J = k(1 − k).  Clustering the nodes (s·k = 0.025, floor
    0.035) resolves it at the same N.  N-doubling with the default rule (N = 100 | 200 | 240; fixed s = 0.5, N = 100 in
    brackets): (0.9, 0.04) 0.0332469 | 0.0332455 | 0.0332455 [0]; (0.9, 0.06) 0.0201206 | 0.0201192 | 0.0201192 [0];
    (0.9, 0.08) 0.0067666 | 0.0067652 | 0.0067652 [0]; (0.85, 0.08) 0.0326887 | 0.0326878 | 0.0326878 [0]; (0.8, 0.13)
    0.0215916 | 0.0215917 | 0.0215917 [0]; (0.95, 0.01) 0.0242820 | 0.0242804 | 0.0242804 [0]; (0.95, 0.04) 0.0049144 |
    0.0049133 | 0.0049133 [0]; (0.65, 0.22) 0.0062373 | 0.0062384 | 0.0062384 [0]; (0.5, 0.24) 0.0098401 (all) [0];
    (0.4, 0.1) 0.1244766 | 0.1244767 | 0.1244767 — largest difference 1.6e-6 absolute.
    DEVIATION (truncated domain): the book's layer is unbounded; ψ̂ = 0 is imposed at ±y_max.  Default None → the box scales
    with the wavelength, y_max = :func:`decay_box`(k) = max(30, 12/k) (12 e-folds of e^{−k|z|}); an explicit y_max overrides
    it.  Checked by doubling (default: N = 100, 12 e-folds | 24 e-folds | N = 160, twice the box): k = 0.05, J = 0.04: 0.0163047
    | 0.0163048 | 0.0163046; k = 0.05, J = 0.045: 0.0082157 | 0.0082158 | 0.0082156; k = 0.1, J = 0.08: 0.0211430 | 0.0211431
    | 0.0211428; k = 0.15, J = 0.125: 0.0058248 | 0.0058252 | 0.0058244; k = 0.05, J = 0.05 (outside the tongue
    J < k(1 − k) = 0.0475): 0.  A fixed y_max = 30 gave 0.0205913 (+26 %) at (0.05, 0.04) and a spurious 0.0077967 at
    (0.05, 0.05).
    Returns kc_i [–] (scaled by U₀/L).  Scalar-callable.
    Example: J = 0, k = 0.4449 → 0.1897; J = 0.1, k = 0.4 → 0.1245; J ≥ ¼ → 0.
    Near neutral (measured, R = 1, defaults): a 0 is returned although the layer is unstable only in a sliver directly under
    J = k(1 − k), where the weak mode (kc_i ≲ 0.004) fails the N-convergence filter.  Its width in J,
    found by bisection in each column, grows smoothly with k: 0.0002 (k = 0.05), 0.0004 (0.1), 0.0007 (0.2), 0.0011 (0.3),
    0.0016 (0.4), 0.0021 (0.5), 0.0027 (0.6), 0.0032 (0.7), 0.0043 (0.8), 0.0054 (0.9), 0.0059 (0.95); the growth lost there
    is at most 0.0011 … 0.0038 (2 % of the map's maximum 0.19).  Outside that sliver the value is converged in N: on the
    whole 20 × 31 default map N = 100 and N = 200 differ by ≤ 2.1e-6 absolute (2.2e-4 relative, at (0.95, 0.04)) with the
    same zero pattern.  Label: converged (0 within 0.006 of the neutral curve means "kc_i < 0.004", not "stable").
    """
    if y_max is None:
        y_max = decay_box(k)  # box scales with 1/k: the mode decays like e^{−k|z|}
    if map_scale is None:
        map_scale = decay_map_scale(k)  # nodes follow the wave: s = 0.025/k within [0.035, 0.5] (critical-layer resolution)
    pr = richardson_profiles("tanh", J, R)
    c = taylor_goldstein_eigs(k, pr["U"], pr["Upp"], pr["N2"], domain=(-1, 1), N=N, bc="decay", map_scale=map_scale,
                              y_max=y_max)
    return float(k * c[0].imag) if len(c) else 0.0


def tg_growth_map(ks=None, Js=None, R: float = 1.0, N: int = 100, cache: bool = True, map_scale: float | None = None,
                  fast: bool = False, y_max: float | None = None) -> dict:
    """kc_i(k, J) of :func:`tg_growth` on a grid (cached; written to ``reference/ch11/tg_growth_map.csv`` by the tables
    script).  Defaults (Part C.5 5.6): k = 0.05 … 1.0 step 0.05, J = 0 … 0.3 step 0.01 (fast: step 0.1, 0.05).
    y_max : box half-width; default None → each column uses :func:`decay_box`(k) = max(30, 12/k) (DEVIATION, truncated
    domain: the box scales with 1/k; checked by doubling, see :func:`tg_growth`); a number forces one fixed box.
    map_scale : tan-map scale; default None → each column uses :func:`decay_map_scale`(k) = min(0.5, max(0.035, 0.025/k))
    (DEVIATION, node distribution: resolves the critical layer of near-neutral modes at N = 100, see :func:`tg_growth`); a
    number forces one fixed scale.  N : Chebyshev degree of every column (default 100).
    Returns dict(k, J, kci (len(J) × len(k)), neutral_J (= k(1 − k) for R = 1)).  With ``cache`` and the default grid
    (R = 1, N = 100, y_max None, map_scale None) the published 4-s.f. table reference/ch11/tg_growth_map.csv is read if
    present (≈ 2 min to recompute).
    Near neutral (measured on the default grid): every grid point with J < k(1 − k) carries growth — the zero strip under
    the exact neutral curve is 0.0075 wide in the columns k = 0.05, 0.15, …, 0.95 and 0.01 in the columns k = 0.1, 0.2, …,
    0.9 (there the neutral curve passes through a grid point, where kc_i = 0 exactly), i.e. never more than the one grid
    step in J; between grid points :func:`tg_growth` still returns 0 within 0.0002–0.006 of the curve.  N = 200 gives the
    same zero pattern.  Label: converged (0 within 0.006 of the neutral curve means "kc_i < 0.004", not "stable").
    """
    default_grid = (ks is None and Js is None and not fast and R == 1.0 and N == 100 and map_scale is None
                    and y_max is None)
    p_ref = _root() / "reference" / "ch11" / "tg_growth_map.csv"
    if cache and default_grid and p_ref.exists():
        rows = [ln for ln in p_ref.read_text(encoding="utf-8").splitlines() if ln and not ln.startswith("#")]
        kk = np.array([float(h.split("=")[1]) for h in rows[0].split(",")[1:]])
        data = np.array([[float(v) for v in ln.split(",")] for ln in rows[1:]])
        return dict(k=kk, J=data[:, 0], kci=data[:, 1:], neutral_J=kk * (1.0 - kk))
    if ks is None:
        ks = np.round(np.arange(0.05, 1.0001, 0.1 if fast else 0.05), 10)
    if Js is None:
        Js = np.round(np.arange(0.0, 0.3001, 0.05 if fast else 0.01), 10)
    ks, Js = _F(ks), _F(Js)

    def compute():
        kci = np.array([[tg_growth(k, J, R, N, map_scale, y_max) for k in ks] for J in Js])
        return dict(k=ks, J=Js, kci=kci, neutral_J=ks * (1.0 - ks))

    # the "box" and "s" entries key the cache on the truncation and node rules, so maps computed with the old fixed
    # y_max = 30 or the old fixed map scale 0.5 are not reused
    return _cached("tg_growth_map", dict(ks=list(ks), Js=list(Js), R=R, N=N,
                                         s="min(0.5,max(0.035,0.025/k))" if map_scale is None else float(map_scale),
                                         box="max(30,12/k)" if y_max is None else float(y_max)), compute, cache)


def _profile_derivs(grid, f, fp=None, fpp=None):
    y = grid.y
    F0 = _F(f(y) if callable(f) else f) * np.ones(len(y))
    F1 = (_F(fp(y)) * np.ones(len(y)) if callable(fp) else (grid.D1 @ F0 if fp is None else _F(fp)))
    F2 = (_F(fpp(y)) * np.ones(len(y)) if callable(fpp) else (grid.D2 @ F0 if fpp is None else _F(fpp)))
    return F0, F1, F2


def _weights(g, w):
    return g.w if w is None else _F(w)


def richardson_identity_check(k: float, c: complex, psi, z, U, Up, Upp, N2, w=None, grid=None) -> dict:
    """Residuals of the integral identity (11.65) and of its imaginary part for a computed Taylor–Goldstein mode.

    Book: §11.7, Eq. (11.65) ∫(N² − ¼U′²)/(U − c)|φ|²dz = ∫(U − c)(|φ′|² + k²|φ|²)dz + ∫½U″|φ|²dz with φ = ψ̂/(U − c)^{1/2}
    (11.63); imaginary part c_i∫(N² − ¼U′²)/|U − c|²|φ|²dz = −c_i∫(|φ′|² + k²|φ|²)dz.  Needs c_i ≠ 0 and ψ̂ = 0 at the ends.
    Parameters: k, c; psi (ψ̂ on the nodes z); z; U, Up, Upp, N2 (callables or arrays on z); w (quadrature weights; default the
    grid's); grid (the solver's SpectralGrid, from ``return_vectors=True`` — required for mapped grids).
    Returns dict(lhs, rhs, residual (relative), imag_lhs, imag_rhs).  Validation: V4 residual ≤ 1e-8.  Label: conserved.
    """
    g = _grid_matching(z, grid)
    wq = _weights(g, w)
    U0, U1, U2 = _profile_derivs(g, U, Up, Upp)
    N2v = _F(N2(g.y) if callable(N2) else N2) * np.ones(len(g.y))
    phi = np.asarray(psi, dtype=complex) / np.sqrt(U0 - c + 0j)
    dphi = g.D1 @ phi
    a2 = np.abs(phi) ** 2
    lhs = np.sum(wq * (N2v - 0.25 * U1 ** 2) / (U0 - c) * a2)
    rhs = np.sum(wq * (U0 - c) * (np.abs(dphi) ** 2 + k ** 2 * a2)) + np.sum(wq * 0.5 * U2 * a2)  # Eq. (11.65)
    ci = complex(c).imag
    il = ci * np.sum(wq * (N2v - 0.25 * U1 ** 2) / np.abs(U0 - c) ** 2 * a2)
    ir = -ci * np.sum(wq * (np.abs(dphi) ** 2 + k ** 2 * a2))
    return dict(lhs=complex(lhs), rhs=complex(rhs), residual=float(abs(lhs - rhs) / max(abs(lhs), abs(rhs), 1e-300)),
                imag_lhs=float(il), imag_rhs=float(ir))


def howard_identity_check(k: float, c: complex, psi, z, U, N2, w=None, grid=None, Up=None) -> dict:
    """Residuals of Howard's identities (11.69)–(11.70) for a computed Taylor–Goldstein mode; c_r as a Q-weighted mean of U.

    Book: §11.7, (11.68) F = ψ̂/(U − c), Q = |F′|² + k²|F|²; (11.69) ∫[(U − c_r)² − c_i²]Q dz = ∫N²|F|²dz; (11.70)
    c_i∫(U − c_r)Q dz = 0 ⇒ c_r = ∫UQ/∫Q (lies in the range of U, (11.71)).
    Returns dict(res_69, res_70 (relative), cr_mean (∫UQ/∫Q), int_11_72 (∫[U² − c_r² − c_i²]Q ≥ 0), in_semicircle).
    Label: conserved.
    """
    g = _grid_matching(z, grid)
    wq = _weights(g, w)
    U0, _, _ = _profile_derivs(g, U, Up, None)
    N2v = _F(N2(g.y) if callable(N2) else N2) * np.ones(len(g.y))
    F = np.asarray(psi, dtype=complex) / (U0 - c)
    Q = np.abs(g.D1 @ F) ** 2 + k ** 2 * np.abs(F) ** 2
    cr, ci = complex(c).real, complex(c).imag
    l69, r69 = np.sum(wq * ((U0 - cr) ** 2 - ci ** 2) * Q), np.sum(wq * N2v * np.abs(F) ** 2)
    l70, s70 = ci * np.sum(wq * (U0 - cr) * Q), abs(ci) * np.sum(wq * np.abs(U0 - cr) * Q)
    return dict(res_69=float(abs(l69 - r69) / max(abs(l69), abs(r69), 1e-300)), res_70=float(abs(l70) / max(s70, 1e-300)),
                cr_mean=float(np.sum(wq * U0 * Q) / np.sum(wq * Q)), int_11_72=float(np.sum(wq * (U0 ** 2 - cr ** 2 - ci ** 2) * Q)),
                in_semicircle=bool(in_howard_semicircle(c, float(np.min(U0)), float(np.max(U0)))))


def stratified_energy_budget(k: float, c: complex, psi, z, U, N2, grid=None, Up=None, w=None) -> dict:
    """Energy budget of a stratified-shear normal mode: kinetic + available potential energy (Exercise 11.10).

    Book: Ex. 11.10: ½ d/dt∫(u² + w² + g²ρ²/(ρ₀²N²))dV = −∫uw ∂U/∂z dV.  With u = ∂ψ/∂z, w = −∂ψ/∂x (§11.7), û = ψ̂′,
    ŵ = −ikψ̂, gρ̂/ρ₀ = −N²ψ̂/(U − c) (11.60): KE = ¼∫(|û|² + |ŵ|²), APE = ¼∫N²|ψ̂|²/|U − c|², P = −½∫Re(ûŵ*)U′,
    dK/dt = 2kc_i KE, dP/dt = 2kc_i APE; residual = dK/dt + dP/dt − P.  Needs N² > 0 where ψ̂ ≠ 0.
    Returns dict(dKdt, dPdt, production, residual, KE, APE).  Validation: V4 residual ≈ 0.  Label: conserved.
    """
    g = _grid_matching(z, grid)
    wq = _weights(g, w)
    U0, U1, _ = _profile_derivs(g, U, Up, None)
    N2v = _F(N2(g.y) if callable(N2) else N2) * np.ones(len(g.y))
    psi = np.asarray(psi, dtype=complex)
    uh, wh = g.D1 @ psi, -1j * k * psi
    KE = 0.25 * np.sum(wq * (np.abs(uh) ** 2 + np.abs(wh) ** 2))
    APE = 0.25 * np.sum(wq * N2v * np.abs(psi) ** 2 / np.abs(U0 - c) ** 2)
    P = np.sum(wq * (-0.5) * np.real(uh * np.conj(wh)) * U1)
    gr = 2 * k * complex(c).imag
    return dict(dKdt=float(gr * KE), dPdt=float(gr * APE), production=float(P), residual=float(gr * (KE + APE) - P),
                KE=float(KE), APE=float(APE))


# ======================================================================================================================
# §11.8 Squire and Orr–Sommerfeld
# ======================================================================================================================
def os_3d_eigs(k: float, m: float, Re: float, U: Callable, Up: Callable | None = None, N: int = 60,
               domain: Sequence[float] = (-1.0, 1.0)) -> np.ndarray:
    """Eigenvalues c of the 3-D normal-mode equations in primitive variables (û, v̂, ŵ, p̂), wall-bounded flow.

    Book: §11.8, Eq. (11.77) ik(U − c)û + v̂U′ = −ikp̂ + (1/Re)[û″ − (k² + m²)û], ik(U − c)v̂ = −p̂′ + (1/Re)[v̂″ − (k² + m²)v̂],
    ik(U − c)ŵ = −imp̂ + (1/Re)[ŵ″ − (k² + m²)ŵ], ikû + v̂′ + imŵ = 0; û = v̂ = ŵ = 0 at the walls.  Only for the numerical test
    of Squire's theorem (leading c equals :func:`orr_sommerfeld_eigs` at (k̄, Re̅), (11.78)).
    **The fifth argument is U′ (the v̂U′ term), not U″** (design Part C writes ``Upp`` here — a contract slip; passing U″ would
    be wrong physics).  Up=None → U′ from the Chebyshev matrix.
    Returns c sorted by descending c_i (finite, |c| < 10).  Example: Poiseuille k = 1, m = 0.5, Re = 10⁴ equals OS at
    k̄ = 1.1180, Re̅ = 8944.3 to 1e-10.  Label: converged.
    """
    D, y = cheb(N, domain)
    n = N + 1
    I = np.eye(n)
    Z = np.zeros((n, n))
    Lv = (D @ D - (k ** 2 + m ** 2) * I) / Re
    Uy = np.asarray(U(y), dtype=float) * np.ones(n)
    Upy = D @ Uy if Up is None else np.asarray(Up(y), dtype=float) * np.ones(n)
    ikU = 1j * k * np.diag(Uy)
    A = np.block([[ikU - Lv, np.diag(Upy), Z, 1j * k * I],
                  [Z, ikU - Lv, Z, D],
                  [Z, Z, ikU - Lv, 1j * m * I],
                  [1j * k * I, D, 1j * m * I, Z]]).astype(complex)  # Eq. (11.77)
    B = np.block([[1j * k * I, Z, Z, Z], [Z, 1j * k * I, Z, Z], [Z, Z, 1j * k * I, Z], [Z, Z, Z, Z]]).astype(complex)
    elim = [0, N, n, n + N, 2 * n, 2 * n + N]
    C = np.zeros((6, 4 * n))
    for r, j in enumerate(elim):
        C[r, j] = 1.0
    w = constrained_eig(A, B, C, elim)
    return w[np.abs(w) < 10]


def os_derivation_sympy(printed_v_sign: bool = False) -> dict:
    """D20–D21: the Orr–Sommerfeld equation (11.79) from (11.77) with m = ŵ = 0, û = φ′, v̂ = −ikφ.

    Book: §11.8, p. 510.  Solve the x-equation for p̂, differentiate, insert in the y-equation: the result is
    −1 × [(11.79) residual] (U − c)(φ″ − k²φ) − U″φ − (1/(ikRe))(φ⁗ − 2k²φ″ + k⁴φ).  The planted v̂ = +ikφ leaves a remainder.
    Parameters: printed_v_sign (put the planted variant into ``residual``).
    Returns dict(residual (0), planted_residual (≠ 0), factor (−1), steps).  Validation: V2.  Label: symbolic.
    """
    y = sp.symbols("y", real=True)
    k, Re = sp.symbols("k Re", positive=True)
    c = sp.symbols("c")
    U, phi = sp.Function("U")(y), sp.Function("phi")(y)
    OS = (U - c) * (phi.diff(y, 2) - k ** 2 * phi) - U.diff(y, 2) * phi - (phi.diff(y, 4) - 2 * k ** 2 * phi.diff(y, 2)
                                                                        + k ** 4 * phi) / (sp.I * k * Re)  # Eq. (11.79)

    def remainder(sign):
        uh, vh = phi.diff(y), sign * sp.I * k * phi
        p = (-sp.I * k * (U - c) * uh - vh * U.diff(y) + (uh.diff(y, 2) - k ** 2 * uh) / Re) / (sp.I * k)  # (11.77) x
        E = sp.expand(sp.I * k * (U - c) * vh + p.diff(y) - (vh.diff(y, 2) - k ** 2 * vh) / Re)  # (11.77) y
        fac = sp.simplify(E.coeff(phi.diff(y, 4)) / sp.expand(OS).coeff(phi.diff(y, 4)))
        return fac, sp.simplify(sp.expand(E - fac * OS))

    fac, res = remainder(-1)
    _, planted = remainder(+1)
    return dict(residual=planted if printed_v_sign else res, planted_residual=planted, factor=fac,
                steps=["û = dφ/dy, v̂ = −ikφ (u = ∂ψ/∂y, v = −∂ψ/∂x)", "x-equation of (11.77) solved for p̂",
                       "d/dy of p̂ inserted into the y-equation", "collect: −[(U − c)(φ″ − k²φ) − U″φ − (φ⁗ − 2k²φ″ + k⁴φ)/(ikRe)]"])


def energy_equation_sympy() -> dict:
    """D23 (Exercise 11.13): the pointwise identity behind (11.88), 2-D, U = U(y), from (11.96).

    Book: §11.10 (11.88); Ex. 11.13, Eq. (11.96).  Subtract the basic state from (11.96), multiply by u_i (∇·u = 0 imposed by
    u = ψ_y, v = −ψ_x): u_iR_i = ∂_t(½u_iu_i) + ∂_j[½u_iu_i(U_j + u_j) + pu_i/ρ − νu_i∂_ju_i] + u_iu_j∂_jU_i + ν(∂_ju_i)².
    Over the control volume of Fig. 11.25 the divergence integrates to zero, leaving (11.88).
    Returns dict(residual (0), terms (names)).  Validation: V2.  Label: symbolic.
    """
    x, y, t, rho, nu = sp.symbols("x y t rho nu", positive=True)
    psi, p, U = sp.Function("psi")(x, y, t), sp.Function("p")(x, y, t), sp.Function("U")(y)
    u, v = psi.diff(y), -psi.diff(x)
    lap = lambda f: f.diff(x, 2) + f.diff(y, 2)  # noqa: E731
    Rx = u.diff(t) + (U + u) * u.diff(x) + v * (U.diff(y) + u.diff(y)) + p.diff(x) / rho - nu * lap(u)  # (11.96) − basic
    Ry = v.diff(t) + (U + u) * v.diff(x) + v * v.diff(y) + p.diff(y) / rho - nu * lap(v)
    q = (u ** 2 + v ** 2) / 2
    fx = q * (U + u) + p * u / rho - nu * (u * u.diff(x) + v * v.diff(x))
    fy = q * v + p * v / rho - nu * (u * u.diff(y) + v * v.diff(y))
    rhs = q.diff(t) + fx.diff(x) + fy.diff(y) + u * v * U.diff(y) + nu * (u.diff(x) ** 2 + u.diff(y) ** 2 + v.diff(x) ** 2
                                                                          + v.diff(y) ** 2)
    return dict(residual=sp.simplify(sp.expand(u * Rx + v * Ry - rhs)),
                terms=["∂t(½u_iu_i)", "divergence of energy flux (advection, pressure work, viscous diffusion)",
                       "u_iu_j ∂U_i/∂x_j (minus the production)", "ν(∂_ju_i)² (dissipation Λ)"])


# ======================================================================================================================
# §11.9 inviscid criteria
# ======================================================================================================================
def rayleigh_criterion(y, U=None, Upp=None) -> dict:
    """Rayleigh's inflection-point criterion (necessary for inviscid instability): does U″ change sign inside (y₁, y₂)?

    Book: §11.9, Eqs. (11.83)–(11.84).  Parameters as :func:`core.stability.inflection_points`.
    Returns dict(has_inflection, y_I).  Label: analytic.
    """
    yI = inflection_points(y, U=U, Upp=Upp)
    return dict(has_inflection=bool(yI.size), y_I=yI)


def fjortoft_criterion(y, U=None, Upp=None, tol: float = 1e-12) -> dict:
    """Fjørtoft: (U − U_I)U″ < 0 somewhere, U_I = U at an inflection point (stronger necessary condition).

    Book: §11.9, Eqs. (11.85)–(11.86) and their sum ∫(U − U_I)U″/|U − c|²|φ|²dy < 0 (the book writes U₁ for U_I); equivalently
    |vorticity| has a maximum inside the flow.  Parameters: y (grid); U callable; Upp callable (default central difference).
    Returns dict(satisfied, U_I (list), y_I, min_product (min of (U − U_I)U″ on the grid)).
    Validation: V1 tanh, Bickley, sin y: yes; sinh(2y)/sinh 2: Rayleigh yes, Fjørtoft no.  Label: analytic.
    """
    y = np.sort(_F(y))
    if not callable(U):
        raise ValueError("fjortoft_criterion needs U as a callable")
    f2 = Upp if Upp is not None else (lambda yy: _second_derivative(U, yy))
    yI = inflection_points(y, Upp=f2)
    if not yI.size:
        return dict(satisfied=False, U_I=[], y_I=yI, min_product=float("nan"))
    UI = [float(U(v)) for v in yI]
    prods = [(_F(U(y)) - ui) * (_F(f2(y)) * np.ones(len(y))) for ui in UI]
    scale = max(max(float(np.max(np.abs(pq))) for pq in prods), 1e-300)
    sat = any(bool(np.any(pq < -tol * scale)) for pq in prods)
    return dict(satisfied=sat, U_I=UI, y_I=yI, min_product=float(min(np.min(pq) for pq in prods)))


def rayleigh_identity_check(k: float, c: complex, phi, y, U, Upp, w=None, grid=None) -> dict:
    """Residuals of (11.83) (real part) and (11.84) (imaginary part) for a computed Rayleigh mode (c_i ≠ 0).

    Book: §11.9, Eq. (11.83) ∫(|φ′|² + k²|φ|²)dy + ∫U″/(U − c)|φ|²dy = 0; Eq. (11.84) c_i∫U″/|U − c|²|φ|²dy = 0.
    Returns dict(res_real, res_imag (relative)).  Validation: V4 ~1e-12.  Label: conserved.
    """
    g = _grid_matching(y, grid)
    wq = _weights(g, w)
    U0 = _F(U(g.y) if callable(U) else U) * np.ones(len(g.y))
    U2 = _F(Upp(g.y) if callable(Upp) else Upp) * np.ones(len(g.y))
    phi = np.asarray(phi, dtype=complex)
    a2 = np.abs(phi) ** 2
    T1 = np.sum(wq * (np.abs(g.D1 @ phi) ** 2 + k ** 2 * a2))
    T2 = np.sum(wq * U2 / (U0 - c) * a2)
    i84 = np.sum(wq * U2 / np.abs(U0 - c) ** 2 * a2)
    s84 = np.sum(wq * np.abs(U2) / np.abs(U0 - c) ** 2 * a2)
    return dict(res_real=float(abs(np.real(T1 + T2)) / max(abs(T1), 1e-300)), res_imag=float(abs(i84) / max(s84, 1e-300)))


def critical_layer(y, U: Callable, c) -> np.ndarray:
    """Critical-layer positions y_c where U(y_c) = c_r (bracketed on y, refined by Brent).

    Book: §11.9, p. 514 (neutral modes need U = c somewhere; y_c is a critical point of (11.81)).  Returns sorted ndarray.
    Label: analytic.
    """
    y = np.sort(_F(y))
    cr = complex(c).real
    f = lambda yy: U(yy) - cr  # noqa: E731
    v = f(y)
    roots = [y[j] for j in range(len(y)) if v[j] == 0]
    roots += [brentq(f, y[j], y[j + 1], xtol=1e-14) for j in range(len(y) - 1) if v[j] * v[j + 1] < 0]
    return np.array(sorted(roots))


def cats_eye_streamfunction(x, y, y_c: float = 0.0, A: float = 0.1, phi_c: float = 1.0, k: float = 1.0,
                            Uy_c: float = 1.0, U: Callable | None = None, c: float | None = None, exact: bool = False,
                            phi: Callable | None = None):
    """Kelvin cat's-eye stream function near a critical layer, in the frame moving with c = U(y_c).

    Book: §11.9, Eq. (11.87) ψ̂ = ∫(U − c)dy + Aφ(y)e^{ikx} (real part) and its expansion (p. 514)
    ψ̂ ≅ ½(y − y_c)²[dU/dy]_{y_c} + Aφ(y_c)cos(kx) (φ(y_c) real); Fig. 11.22.
    Parameters: x, y (broadcastable); y_c; A; phi_c; k; Uy_c (U′(y_c); replaced by U′ from ``U`` when U is given);
    U (callable; needed for exact); c (default U(y_c)); exact (∫_{y_c}^{y}(U − c)dy + Aφ cos kx, φ = ``phi`` or phi_c).
    Returns ψ̂.  Label: analytic.
    """
    x, y = np.broadcast_arrays(_F(x), _F(y))
    if U is not None:
        h = 1e-6
        Uy_c = (U(y_c + h) - U(y_c - h)) / (2 * h)
    if not exact:
        return _S(0.5 * (y - y_c) ** 2 * Uy_c + A * phi_c * np.cos(k * x))  # p. 514
    if U is None:
        raise ValueError("exact=True needs U")
    cc = U(y_c) if c is None else c
    ys = np.unique(y)
    vals = np.array([quad(lambda s: U(s) - cc, y_c, yy, epsabs=1e-13)[0] for yy in ys])
    ph = phi(y) if phi is not None else phi_c
    return _S(np.interp(y, ys, vals) + A * ph * np.cos(k * x))  # Eq. (11.87)


def cats_eye_width(A: float, phi_c: float, Uy_c: float) -> float:
    """Half-width of the cat's eye, centre line to separatrix: 2√(Aφ_c/U′(y_c)) (ours, from the p. 514 expansion: the
    separatrix ψ̂ = Aφ_c meets x where cos kx = −1 at |y − y_c| = 2√(Aφ_c/U′_c)).  Example: A = 0.1, φ_c = 1, U′ = 1 → 0.632.
    Label: analytic.
    """
    return 2.0 * math.sqrt(A * phi_c / Uy_c)


def piecewise_shear_layer_c(kh, dU: float = 1.0, U_mean: float = 0.0):
    """Phase speed of the piecewise-linear (constant-vorticity) shear layer of thickness h (Exercise 11.11).

    Book: Ex. 11.11(c): c₀² = ((U₁ − U₃)/(2kh))²{(kh − 1)² − e^{−2kh}}, c₀ = c − ½(U₁ + U₃); unstable for kh < 1.2785;
    small kh: c₀ ≅ ±i((U₁ − U₃)/2)√(1 − (4/3)kh + …) (part (d)).  Parameters: kh (array ok); dU = U₁ − U₃; U_mean.
    Returns c (growing root, Im ≥ 0).  Validation: V1 neutral kh = 1.278465; max kc_i h/ΔU = 0.2012 at kh = 0.797.
    Label: analytic.
    """
    kh = _F(kh)
    c02 = (dU / (2.0 * kh)) ** 2 * ((kh - 1.0) ** 2 - np.exp(-2.0 * kh))  # Ex. 11.11(c)
    return _S(U_mean + np.lib.scimath.sqrt(c02))


def piecewise_neutral_kh() -> float:
    """Largest unstable kh of Ex. 11.11(e): root of (kh − 1)² = e^{−2kh} (= 1.278465).  Label: analytic."""
    return float(brentq(lambda a: (a - 1.0) ** 2 - math.exp(-2.0 * a), 1.0, 2.0, xtol=1e-15))


# ======================================================================================================================
# §11.9–§11.11 base profiles
# ======================================================================================================================
def _sech2(y):
    # clip: cosh² overflows past |y| ≈ 355 (boxes of decay_box(k) for k < 0.034); sech²(350) ~ 1e-304 either way
    if np.iscomplexobj(y):  # complex path (rayleigh_eigs_contour): clip the real part only
        y = np.asarray(y)
        return 1.0 / np.cosh(np.clip(y.real, -350.0, 350.0) + 1j * y.imag) ** 2
    return 1.0 / np.cosh(np.clip(y, -350.0, 350.0)) ** 2


def blasius_base(y_max: float = 20.0) -> dict:
    """Blasius base flow in δ* units for Orr–Sommerfeld: U(y) = f′(δ*_η y), U′ = f″δ*_η, U″ = −½ff″δ*_η² (exact, from (9.27)).

    Book: §11.11 (Blasius stability, Re based on δ*); Ch. 9 (9.27) f‴ + ½ff″ = 0, δ*_η = lim(η − f) = 1.7208 (computed by
    ``core.boundary_layer``).  Returns dict(U, Up, Upp (callables of y = y_dim/δ*), delta_star_eta, domain (0, 1), bc
    "semi_infinite", y_max).  Label: converged.
    """
    from .core.boundary_layer import blasius_constants, blasius_profile

    ds = blasius_constants()["delta_star"]

    def U(y):
        return blasius_profile(_F(y) * ds)[1]

    def Up(y):
        return blasius_profile(_F(y) * ds)[2] * ds

    def Upp(y):
        f, fp, fpp = blasius_profile(_F(y) * ds)
        return -0.5 * f * fpp * ds ** 2  # f‴ = −½ff″ (9.27)

    return dict(U=U, Up=Up, Upp=Upp, delta_star_eta=ds, domain=(0.0, 1.0), bc="semi_infinite", y_max=y_max)


def _falkner_skan_callables(m: float, eta_max: float | None = None, n: int = 3001):
    from .core.boundary_layer import falkner_skan

    sol = falkner_skan(m, eta_max=eta_max, n=n)
    if not sol["success"]:
        raise RuntimeError(f"falkner_skan(m={m}) did not converge")
    eta = sol["eta"]
    ds = float(eta[-1] - sol["f"][-1])  # δ*/δ = lim(η − f)
    sf, sfp, sfpp = (CubicSpline(eta, sol[key]) for key in ("f", "fp", "fpp"))
    e_end = float(eta[-1])

    def _ev(y):
        e = np.clip(_F(y) * ds, 0.0, e_end)
        return sf(e), sfp(e), sfpp(e), _F(y) * ds > e_end

    def U(y):
        f, fp, fpp, far = _ev(y)
        return np.where(far, 1.0, fp)

    def Up(y):
        f, fp, fpp, far = _ev(y)
        return np.where(far, 0.0, fpp * ds)

    def Upp(y):
        f, fp, fpp, far = _ev(y)
        return np.where(far, 0.0, (-(m + 1.0) / 2.0 * f * fpp + m * fp ** 2 - m) * ds ** 2)  # (9.36) solved for f‴

    return U, Up, Upp, ds


_PROFILE_ALIASES = {"jet": "bickley", "shear_layer": "tanh", "blasius_like": "blasius"}


def parallel_profile(name: str, **p) -> dict:
    """Base profiles of §11.9–§11.11 as callables with the domain, boundary type and scales each solver needs.

    Book: §11.9 Fig. 11.21; §11.10 (tanh, sech², Poiseuille, Couette, pressure-gradient boundary layers, Table 11.1); §11.11.
    Names (y by L, U by U₀, Re = U₀L/ν): "poiseuille" 1 − y² on ±1 (half-width, centreline speed); "couette" y on ±1;
    "tanh" / "shear_layer" tanh y (unbounded, bc "decay"); "bickley" / "jet" sech²y (unbounded; parity "even" = sinuous);
    "blasius" / "blasius_like" f′ in δ* units (semi-infinite); "falkner_skan" (m) in δ* units; "sin" sin y on |y| ≤ b
    (default π/2); "wall_vorticity_max" sinh(βy)/sinh β (β = 2; inflection with the vorticity maximum at the walls, Fig. 11.21d);
    "shear_layer_walls" tanh(y/δ)/tanh(1/δ) (δ = 0.3; Fig. 11.21f).
    Overrides: y_max, map_scale, b, beta, delta, m, parity.
    Returns dict(U, Up, Upp, domain, bc, y_max, map_scale, parity, label, length, velocity, fixed_box).  ``y_max`` of an
    unbounded or semi-infinite profile is the *smallest* box (30 tanh, 40 Bickley, 20 Blasius, 25 Falkner–Skan); the
    Orr–Sommerfeld wrappers enlarge it to :func:`decay_box`(k) for long waves (see :func:`os_box_numerics`) unless
    ``y_max`` was passed here explicitly (then ``fixed_box`` is True and the box is exactly that).  Label: analytic
    (Blasius/FS: converged).
    """
    name = _PROFILE_ALIASES.get(name, name)
    if name == "poiseuille":
        d = dict(U=lambda y: 1.0 - _FC(y) ** 2, Up=lambda y: -2.0 * _FC(y), Upp=lambda y: -2.0 + 0.0 * _FC(y), domain=(-1.0, 1.0),
                 bc="wall", label="plane Poiseuille U = 1 − y²", length="half-width", velocity="centreline speed")
    elif name == "couette":
        d = dict(U=lambda y: _FC(y), Up=lambda y: 1.0 + 0.0 * _FC(y), Upp=lambda y: 0.0 * _FC(y), domain=(-1.0, 1.0), bc="wall",
                 label="plane Couette U = y", length="half-gap", velocity="wall speed")
    elif name == "tanh":
        d = dict(U=np.tanh, Up=_sech2, Upp=lambda y: -2.0 * np.tanh(y) * _sech2(y), domain=(-1.0, 1.0), bc="decay",
                 y_max=30.0, map_scale=1.0, label="shear layer U = tanh y", length="L", velocity="U₀")
    elif name == "bickley":
        d = dict(U=_sech2, Up=lambda y: -2.0 * _sech2(y) * np.tanh(y), Upp=lambda y: 4.0 * _sech2(y) - 6.0 * _sech2(y) ** 2,
                 domain=(-1.0, 1.0), bc="decay", y_max=40.0, map_scale=1.0, label="Bickley jet U = sech² y", length="L",
                 velocity="centreline speed")
    elif name == "blasius":
        b = blasius_base(20.0 if p.get("y_max") is None else p["y_max"])
        d = dict(U=b["U"], Up=b["Up"], Upp=b["Upp"], domain=(0.0, 1.0), bc="semi_infinite", y_max=b["y_max"],
                 label="Blasius boundary layer (δ* units)", length="δ*", velocity="U∞", delta_star_eta=b["delta_star_eta"])
    elif name == "falkner_skan":
        m = float(p.get("m", 0.0))
        U, Up, Upp, ds = _falkner_skan_callables(m, p.get("eta_max"))
        d = dict(U=U, Up=Up, Upp=Upp, domain=(0.0, 1.0), bc="semi_infinite", y_max=25.0,
                 label=f"Falkner–Skan m = {m:g} (δ* units)", length="δ*", velocity="U_e", delta_star_eta=ds, m=m)
    elif name == "sin":
        b = float(p.get("b", math.pi / 2))
        d = dict(U=np.sin, Up=np.cos, Upp=lambda y: -np.sin(y), domain=(-b, b), bc="wall", label=f"U = sin y, |y| ≤ {b:.4g}",
                 length="1", velocity="1", b=b)
    elif name == "wall_vorticity_max":
        be = float(p.get("beta", 2.0))
        d = dict(U=lambda y: np.sinh(be * _FC(y)) / math.sinh(be), Up=lambda y: be * np.cosh(be * _FC(y)) / math.sinh(be),
                 Upp=lambda y: be ** 2 * np.sinh(be * _FC(y)) / math.sinh(be), domain=(-1.0, 1.0), bc="wall",
                 label=f"U = sinh({be:g}y)/sinh {be:g}", length="half-width", velocity="wall speed")
    elif name == "shear_layer_walls":
        de = float(p.get("delta", 0.3))
        t1 = math.tanh(1.0 / de)
        d = dict(U=lambda y: np.tanh(_FC(y) / de) / t1, Up=lambda y: _sech2(_FC(y) / de) / (de * t1),
                 Upp=lambda y: -2.0 * np.tanh(_FC(y) / de) * _sech2(_FC(y) / de) / (de ** 2 * t1), domain=(-1.0, 1.0),
                 bc="wall", label=f"U = tanh(y/{de:g})/tanh(1/{de:g})", length="half-width", velocity="wall speed")
    else:
        raise ValueError(f"unknown profile {name!r}")
    out = dict(y_max=None, map_scale=None, parity=p.get("parity"), name=name)
    out.update(d)
    for key in ("y_max", "map_scale"):
        if p.get(key) is not None:
            out[key] = p[key]
    out["fixed_box"] = p.get("y_max") is not None  # an explicit box is honoured exactly (no wavelength scaling)
    return out


def inviscid_profile(name: str, **p) -> dict:
    """Analytic stand-ins of Fig. 11.21 and the §11.9 examples: dict(U, Up, Upp, domain, bc, label, …).

    Book: §11.9.  Names as :func:`parallel_profile` ("blasius_like", "couette", "poiseuille", "wall_vorticity_max", "jet",
    "shear_layer", "bickley", "sin" (b), "shear_layer_walls").  Label: analytic.
    """
    return parallel_profile(name, **p)


FIG_11_21 = dict(a="blasius", b="couette", c="poiseuille", d="wall_vorticity_max", e="sin", f="shear_layer_walls")
"""Our analytic stand-ins for the six profiles of Fig. 11.21 ((e): sin y on |y| ≤ π)."""


def fig_11_21_verdicts() -> list:
    """Rayleigh and Fjørtoft verdicts for our six stand-ins of Fig. 11.21 ((a)–(c) none; (d) Rayleigh only; (e), (f) both).

    Book: §11.9, Fig. 11.21.  Returns list of dict(panel, name, label, rayleigh, fjortoft, y_I).  Label: analytic.
    """
    out = []
    for panel, name in FIG_11_21.items():
        d = parallel_profile(name, b=math.pi) if name == "sin" else parallel_profile(name)
        lo, hi = d["domain"] if d["bc"] == "wall" else (0.0, 12.0)
        y = np.linspace(lo, hi, 801)[1:-1]
        ray = rayleigh_criterion(y, Upp=d["Upp"])
        fj = fjortoft_criterion(y, d["U"], d["Upp"])["satisfied"] if ray["has_inflection"] else False
        out.append(dict(panel=panel, name=name, label=d["label"], rayleigh=ray["has_inflection"], fjortoft=fj,
                        y_I=ray["y_I"]))
    return out


_OS_BOX_RULE = "decay_box(k, y_min=profile, 12); tan s = s0*clip(y_max/100, 1, 3); semi-infinite N*sqrt(y_max/y_min)"
"""Label of the default box rule of :func:`os_box_numerics` (part of every cache key that depends on it)."""


def os_box_numerics(profile: dict, k: float, N: int, n_efold: float = 12.0) -> dict:
    """Box, map scale and degree of an Orr–Sommerfeld solve of a :func:`parallel_profile` at wavenumber k.

    Book: §11.8 (11.79)–(11.80), §11.10–§11.11 (ours: where "φ → 0 as |y| → ∞" is imposed).  DEVIATION: the book has no
    truncation.  A fixed box is too short for long waves (the mode decays like e^{−k|y|}): measured, the Bickley lower branch
    in the box 40 was wrong for every k < 0.1 (Re = 14.58, k = 0.04: c_i = −0.0184 against +0.00125), the Blasius lower
    branch in the box 20 was 8 % low at Re = 6000 (0.0725 against 0.0790), the tanh c_i at k = 0.05, Re = 1 was 0.161 against
    0.282.  Rule (walls: nothing changes):
      * unbounded ("decay"): y_max = :func:`decay_box`(k, profile box, n_efold) = max(profile box, n_efold/k); tan-map scale
        :func:`decay_box_map_scale` = s₀·min(3, max(1, y_max/100)); N unchanged;
      * "semi_infinite" (linear map on [0, y_max]): the same y_max; N·√(y_max/profile box), which keeps the number of nodes
        inside the boundary layer (Chebyshev nodes cluster quadratically at the wall);
      * ``profile["fixed_box"]`` (an explicit ``y_max`` given to :func:`parallel_profile`): exactly that box, s and N.
    Parameters
    ----------
    profile : dict from :func:`parallel_profile`.   k : wavenumber, non-dimensional (scaled by 1/L), > 0.
    N : Chebyshev degree in the profile's smallest box [–].   n_efold : e-folds of e^{−k|y|} inside the box [–]; default 12.
    Returns
    -------
    dict(y_max (non-dimensional, scaled by L; None for walls), map_scale (same; None where unused), N [–]).  Scalar-callable.
    Validation: V3 box doubling and N doubling of the neutral wavenumbers (:func:`bickley_neutral_curve`,
    :func:`blasius_neutral_curve`).  Label: converged.
    """
    bc, ym0, s0 = profile["bc"], profile.get("y_max"), profile.get("map_scale")
    if bc == "wall" or ym0 is None or profile.get("fixed_box"):
        return dict(y_max=ym0, map_scale=s0, N=int(N))
    ym = decay_box(float(k), y_min=float(ym0), n_efold=n_efold)
    if bc == "semi_infinite":
        return dict(y_max=ym, map_scale=s0, N=int(math.ceil(N * math.sqrt(ym / float(ym0)))))
    return dict(y_max=ym, map_scale=decay_box_map_scale(ym, s0=1.0 if s0 is None else float(s0)), N=int(N))


def _os_lead(profile: dict, N: int, parity=None, n_efold: float = 12.0):
    def fn(k, Re):
        nm = os_box_numerics(profile, k, N, n_efold)
        c = orr_sommerfeld_eigs(k, Re, profile["U"], profile["Upp"], domain=profile["domain"], N=nm["N"], bc=profile["bc"],
                                y_max=nm["y_max"], map_scale=nm["map_scale"], filter=False,
                                parity=parity if parity is not None else profile.get("parity"))
        c = c[np.abs(c) < 5.0]
        return c[:1] if c.size else np.array([complex(0.0, -np.inf)])

    return fn


def os_leading_mode(profile: dict, k: float, Re: float, N: int = 80, far_tol: float = 0.2, box_tol: float = 1e-5,
                    n_try: int = 200, refine: bool = True) -> dict:
    """Least-damped *discrete* Orr–Sommerfeld mode of a :func:`parallel_profile` in the wavelength-scaled box, with its
    energy budget — or a statement that none was found and where the continuous spectrum starts.

    Book: §11.8 (11.79)–(11.80), (11.88); §11.10–§11.11.  On an unbounded or semi-infinite profile the spectrum of (11.79) is
    a finite set of discrete modes plus a continuous spectrum c = U∞ − i(k² + s²)/(kRe), s real (the far-field solutions
    e^{±isy}; put U = U∞, U″ = 0 in (11.79)), whose least-damped edge is c_i = −k/Re.  A truncated box turns the continuum
    into box-dependent eigenvalues that are never a mode of the flow; in a stable flow they are often the leading
    eigenvalues of the box (Blasius, Re = 200, k = 0.05: c_i = −0.0120, −0.0026, −0.00066, −0.00032 in boxes 20, 40, 80,
    160 δ*; the edge is −0.00025).  Method (ours): eigenvalues in descending c_i; the first whose eigenfunction is localised
    (:func:`far_field_fraction` < ``far_tol``) **and** which reappears within ``box_tol`` in a second solve in a box 1.5
    times as large is the mode.  If none of the first ``n_try`` localised candidates passes, ``discrete`` is False, ``c`` is
    the analytic edge U∞ − ik/Re (c_r = NaN when the two free streams differ) and there is no budget.  A discrete mode may be
    more damped than the edge (``above_edge`` False): it is still a mode (e.g. a decaying Tollmien–Schlichting wave), but
    then not the least-damped part of the spectrum.  Walls, or an explicit box (``fixed_box``): the leading eigenvalue as is.
    Parameters
    ----------
    profile : dict from :func:`parallel_profile`.   k, Re : non-dimensional (the profile's L, U₀).   N : degree in the
    profile's smallest box (:func:`os_box_numerics`).   far_tol : localisation threshold [–]; default 0.2.
    box_tol : largest |Δc| between the two boxes [U₀]; default 1e-5.   n_try : localised candidates examined [–]; default
    200 (in effect all).  Measured: discrete modes reappear in the larger box to 1e-10 … 3e-6 (the larger value in boxes of
    240–600) and have far-field fractions ≤ 0.05 (0.02–0.05 for a mode just above the edge of the continuous spectrum);
    box eigenvalues of the continuous spectrum with a fraction below 0.2 move by ≥ 2.5e-4.  A strongly damped mode that the
    degree N does not resolve (Blasius, Re = 4464, k = 0.383: found with N = 160, c = 0.3717 − 0.0472i, moving by 2e-4
    between the boxes with N = 80) is not accepted at that N.   refine : if no mode is found, search once more with 2N
    (default True; the returned ``numerics`` then holds the degree used); what is still not found is reported as such,
    never as a wrong number.
    Returns
    -------
    dict(discrete (bool), c (complex), above_edge (bool), budget (dict of :func:`disturbance_energy_budget` or None), mode
    (:func:`os_mode`-like dict or None), far (float), numerics (dict)).  Scalar-callable.
    Validation: V3 (discrete eigenvalues unchanged by doubling the box and N; rejected ones move with the box).
    Label: converged.
    """
    k, Re = float(k), float(Re)
    nm = os_box_numerics(profile, k, N)
    kw = dict(domain=profile["domain"], bc=profile["bc"], map_scale=nm["map_scale"], parity=profile.get("parity"))
    if profile["bc"] == "wall" or profile.get("fixed_box"):
        m = os_mode(k, Re, profile["U"], profile["Upp"], N=nm["N"], y_max=nm["y_max"], **kw)
        b = disturbance_energy_budget(k, m["c"], m["phi"], m["y"], profile["Up"], Re, grid=m["grid"])
        return dict(discrete=True, c=m["c"], above_edge=True, budget=b, mode=m, far=0.0, numerics=nm)
    semi = profile["bc"] == "semi_infinite"
    centre = None if semi else 0.5 * (profile["domain"][0] + profile["domain"][1])
    res = orr_sommerfeld_eigs(k, Re, profile["U"], profile["Upp"], N=nm["N"], y_max=nm["y_max"], filter=False,
                              return_vectors=True, **kw)
    ym2 = 1.5 * nm["y_max"]
    N2 = int(math.ceil(nm["N"] * math.sqrt(1.5))) if semi else nm["N"]
    s2 = nm["map_scale"] if semi else decay_box_map_scale(ym2, s0=1.0 if profile.get("map_scale") is None
                                                          else float(profile["map_scale"]))
    w2 = orr_sommerfeld_eigs(k, Re, profile["U"], profile["Upp"], N=N2, y_max=ym2, filter=False,
                             **{**kw, "map_scale": s2})
    g, tried, far_lead = res["grid"], 0, float("nan")
    for i, c in enumerate(res["c"]):
        if not np.isfinite(c) or abs(c) >= 5.0:
            continue
        far = far_field_fraction(res["phi"][:, i], g.y, nm["y_max"], centre=centre)
        if i == 0 or not np.isfinite(far_lead):
            far_lead = far
        if far >= far_tol:
            continue
        tried += 1
        if np.min(np.abs(w2 - c)) < box_tol:
            phi = res["phi"][:, i]
            dphi = g.D1 @ phi
            scale = dphi[int(np.argmax(np.abs(dphi)))]
            phi, dphi = phi / scale, dphi / scale  # max|û| = 1, as os_mode
            m = dict(c=complex(c), y=g.y, phi=phi, u_hat=dphi, v_hat=-1j * k * phi, dphi=dphi, grid=g, k=k, Re=Re)
            b = disturbance_energy_budget(k, m["c"], phi, g.y, profile["Up"], Re, grid=g)
            return dict(discrete=True, c=m["c"], above_edge=bool(m["c"].imag > -k / Re), budget=b, mode=m, far=far,
                        numerics=nm)
        if tried >= n_try:
            break
    if refine:
        return os_leading_mode(profile, k, Re, 2 * int(N), far_tol, box_tol, n_try, refine=False)
    ymx = nm["y_max"]
    hi = float(profile["U"](ymx))
    lo = hi if semi else float(profile["U"](-ymx))
    cr = hi if abs(hi - lo) < 1e-9 else float("nan")
    return dict(discrete=False, c=complex(cr, -k / Re), above_edge=False, budget=None, mode=None, far=far_lead, numerics=nm)


# ======================================================================================================================
# §11.9 Rayleigh tables
# ======================================================================================================================
def _rayleigh_lead(pr: dict, k: float, N: int, y_max=None, parity=None):
    """Leading growing Rayleigh eigenvalue of a profile dict (complex) or None: contour solver for analytic profiles, the
    real-axis solver for the semi-infinite spline profiles (Blasius-like: no inflection point, no growing mode)."""
    if pr["bc"] == "semi_infinite":
        c = rayleigh_eigs(float(k), pr["U"], pr["Upp"], domain=pr["domain"], N=N, bc=pr["bc"], y_max=y_max,
                          map_scale=pr.get("map_scale"), unstable_only=True, tol=1e-4, ci_min=1e-4, parity=parity)
    else:
        c = rayleigh_eigs_contour(float(k), pr["U"], pr["Up"], pr["Upp"], domain=pr["domain"], N=N, bc=pr["bc"],
                                  y_max=y_max, map_scale=pr.get("map_scale"), parity=parity)
    return complex(c[0]) if len(c) else None


def sin_profile_max_growth(b: float, ks=None, N: int = 80) -> float:
    """Largest inviscid growth rate kc_i of U = sin y between walls at ±b (→ 0 as 2b → π⁺; stable for 2b < π).

    Book: §11.9, p. 513 (Tollmien's counter-example: inflected yet stable for 2b < π).  Rayleigh (11.81)–(11.82) by
    :func:`rayleigh_eigs_contour` (collocation on a path below the critical layer — converged however small c_i is; the
    earlier real-axis solve with an N-filter lost the near-neutral modes and returned 0 for b = 1.6 and a value 4 % low for
    b = 1.7).  Parameters: b [–]; ks (default 25 values in [0.02, 1.0]); N (default 80).
    Returns max kc_i over the k samples (0 if no mode with c_i > 1e-4), non-dimensional.
    Measured (default ks; N = 80 | N = 120): b = 1.5 → 0 | 0 (2b < π: stable); 1.58 → 0.0002427 | 0.0002427; 1.6 → 0.0013462 |
    0.0013462; 1.7 → 0.0118854 | 0.0118854 (k = 0.224); 2 → 0.0596099 | 0.0596099; 3 → 0.1573130 | 0.1573130 (1–4 s each).
    The exact neutral wavenumber is k_n = √(1 − (π/2b)²) (φ = cos(πy/2b), c = 0), so the unstable band 0 < k < k_n shrinks
    to nothing as 2b → π⁺ and the 25 samples cannot see it once k_n < 0.02 or c_i < 1e-4.
    Label: converged.
    """
    ks = np.linspace(0.02, 1.0, 25) if ks is None else _F(ks)
    pr = parallel_profile("sin", b=float(b))
    best = 0.0
    for k in ks:
        c = _rayleigh_lead(pr, float(k), N)
        if c is not None:
            best = max(best, float(k * c.imag))
    return best


_RAYLEIGH_TABLE_NOTE = (
    "ours (fluidpy.ch11_instability.rayleigh_spectrum_table): leading growing Rayleigh eigenvalue (largest c_i), computed "
    "on a complex collocation path below the critical layer (core.stability.rayleigh_eigs_contour) and converged up to "
    "the neutral wavenumber; c_i = 0 and c_r = null mean c_i <= 1e-4: no growing mode, or one within about 0.001 in k of "
    "the exact neutral wavenumber (measured: Bickley jet, sinuous, neutral k = 2: 0 for k > 1.9991; tanh layer, neutral "
    "k = 1: k > 0.99984; sin y on |y| <= pi, neutral k = sqrt(3)/2 = 0.866025: k > 0.86597); not book data")


def rayleigh_spectrum_table(names=None, ks=None, N: int = 120, cache: bool = True, write: bool = True) -> dict:
    """Leading unstable Rayleigh eigenvalue c(k) for the §11.9 profiles (explainer E7, figure F6; our table).

    Book: §11.9, (11.81)–(11.82).  Defaults (Part C.5 5.7): names of :func:`inviscid_profile`, k = 0.1 … 2.0 step 0.1.  c_i = 0
    and c_r = NaN where no growing mode with c_i > 1e-4 exists.  Writes ``reference/ch11/rayleigh_spectra.json`` (≤ 4 s.f.).
    Extra names (not in the default list): "jet_sinuous" (φ even) and "jet_varicose" (φ odd) — "jet" is the larger of the two,
    which is the sinuous mode at every k.
    Method: :func:`rayleigh_eigs_contour` (N, delta = 0.2, filter against ⌈1.5N⌉ with delta = 0.1) for every analytic
    profile; the Blasius-like spline profile keeps the real-axis :func:`rayleigh_eigs` (it has no inflection point).
    DEVIATION (numerical method, ours): the earlier real-axis solve with the N-filter printed 0 where a growing mode exists
    close to a neutral wavenumber — Bickley jet k = 1.8 (true c = 0.63235 + 0.02432i) and 1.9 (0.64965 + 0.01163i) — and
    was up to 2.6e-4 off in c below it (k = 1.7: 0.61500 + 0.03812i against 0.61474 + 0.03806i).
    How near the neutral point a 0 can still hide a growing mode (measured: where c_i falls to the threshold 1e-4): jet
    sinuous k > 1.99910 (neutral 2), jet varicose k > 0.99945 (neutral 1), tanh layer k > 0.99984 (neutral 1),
    sin y on |y| ≤ π k > 0.86597 (neutral √3/2 = 0.866025); there c_i < 1e-4.  Past the neutral wavenumber 0 is exact.
    DEVIATION (truncated domain): for the unbounded profiles (bc "decay": shear layer, jet) φ = 0 is imposed at ±y_max with
    y_max = max(the profile's box (30 tanh, 40 Bickley), :func:`decay_box`(k) = 12/k) — the mode decays like e^{−k|y|}, so the
    fixed boxes were 5.1e-4 (tanh) and 4.4e-4 (jet) off in c at k = 0.1 (the table's 4th significant figure); k ≥ 0.4 rows
    are unchanged by the rule.
    Returns dict(name → dict(k, c_r, c_i)).
    Validation: V3 N = 120 against N = 240 with delta = 0.1: ≤ 1e-8 in c at every entry; independent shooting
    (:func:`rayleigh_shoot`): ≤ 2e-10; V1 exact neutral wavenumbers approached linearly.  Label: converged.
    """
    names = ["blasius_like", "couette", "poiseuille", "wall_vorticity_max", "jet", "shear_layer", "sin",
             "shear_layer_walls"] if names is None else list(names)
    ks = np.round(np.arange(0.1, 2.0001, 0.1), 10) if ks is None else _F(ks)
    variants = {"jet_sinuous": ("jet", "even"), "jet_varicose": ("jet", "odd")}

    def compute():
        out = {}
        for nm in names:
            base, parity = variants.get(nm, (nm, None))
            pr = parallel_profile(base, b=math.pi) if base == "sin" else parallel_profile(base)
            cr, ci = [], []
            for k in ks:
                ym = pr.get("y_max")
                if pr["bc"] == "decay":
                    ym = max(float(ym), decay_box(float(k)))  # box scales with 1/k: the mode decays like e^{−k|y|}
                c = _rayleigh_lead(pr, float(k), N, y_max=ym, parity=parity)
                cr.append(np.nan if c is None else c.real)
                ci.append(0.0 if c is None else c.imag)
            out[f"{nm}__cr"] = np.array(cr)
            out[f"{nm}__ci"] = np.array(ci)
        out["k"] = ks
        return out

    # "box" / "method" key the cache on the truncation rule and the solver, so tables computed with the old fixed boxes or
    # the real-axis solver (false zeros below the neutral wavenumber) are not reused
    raw = _cached("rayleigh_spectrum_table", dict(names=names, ks=list(ks), N=N, box="max(profile,12/k)",
                                                  method="contour(delta=0.2,filter=1.5N@delta/2,tol=1e-6,ci_min=1e-4)"),
                  compute, cache)
    res = {nm: dict(k=_F(raw["k"]), c_r=_F(raw[f"{nm}__cr"]), c_i=_F(raw[f"{nm}__ci"])) for nm in names}
    if write:
        p = _root() / "reference" / "ch11" / "rayleigh_spectra.json"
        p.parent.mkdir(parents=True, exist_ok=True)
        js = {nm: dict(k=[float(f"{v:.4g}") for v in d["k"]],
                       c_r=[None if not np.isfinite(v) else float(f"{v:.4g}") for v in d["c_r"]],
                       c_i=[float(f"{v:.4g}") for v in d["c_i"]]) for nm, d in res.items()}
        p.write_text(json.dumps(dict(note=_RAYLEIGH_TABLE_NOTE, N=N, spectra=js), indent=1), encoding="utf-8")
    return res


# ======================================================================================================================
# §11.10 viscous results — wrappers with caches
# ======================================================================================================================
def poiseuille_spectrum(k: float, Re: float, N: int = 100) -> np.ndarray:
    """The plane-Poiseuille Orr–Sommerfeld spectrum (the Y-shaped branches), all finite c with |c| < 5, sorted by c_i.

    Book: §11.8, (11.79)–(11.80) with U = 1 − y².  Unfiltered (the reader sees the whole spectrum).  Label: converged (TS mode).
    """
    c = orr_sommerfeld_eigs(k, Re, lambda y: 1.0 - y ** 2, lambda y: -2.0 + 0.0 * y, N=N, filter=False)
    return c[np.abs(c) < 5.0]


def poiseuille_critical(N: int = 100, cache: bool = True) -> dict:
    """Critical point of plane Poiseuille flow from the Orr–Sommerfeld equation (cached).

    Book: §11.10 (Table 11.1, book value private; L = half-width, U₀ = centreline speed).  Returns dict(Re_c, k_c, c_r).
    Expect 5772.22, 1.02056, 0.26400.  With ``cache`` and the default N the published value in
    reference/ch11/critical_points.json is used (≈ 12 s to recompute; ``cache=False`` always recomputes).
    Validation: V5 Orszag (1971); V3 N = 60/80/100.  Label: converged, benchmark.
    """
    if cache and N == 100 and (ref := _ref_critical("poiseuille")):
        return ref

    def compute():
        r = critical_point(_os_lead(parallel_profile("poiseuille"), N), (0.95, 1.10), (5000.0, 7000.0))
        return dict(Re_c=r["Re_c"], k_c=r["k_c"], c_r=r["c_r"])

    d = _cached("poiseuille_critical", dict(N=N), compute, cache)
    return {k: float(v) for k, v in d.items()}


def poiseuille_neutral_curve(Re_values=None, N: int = 80, cache: bool = True, k_bounds: Sequence[float] = (0.15, 1.35),
                             n_k: int = 41, fast: bool = False) -> dict:
    """Neutral curve (lower/upper branch k(Re)) of plane Poiseuille flow — the "thumb" (cached).

    Book: §11.10.  Re_values default: 24 log-spaced from 5800 to 10⁶ (fast: 8 to 10⁵).  Returns dict(Re, k_lower, k_upper,
    c_lower, c_upper).  Measured at Re = 10⁴ (N = 80, defaults): unstable k from 0.7972 to 1.0947.  Label: converged.
    """
    if Re_values is None:
        Re_values = np.geomspace(5800.0, 1e5 if fast else 1e6, 8 if fast else 24)
    Re_values = [float(r) for r in np.atleast_1d(Re_values)]

    def compute():
        nc = neutral_curve(_os_lead(parallel_profile("poiseuille"), N), Re_values, k_bounds, n_k=n_k)
        return {key: nc[key] for key in ("Re", "k_lower", "k_upper", "c_lower", "c_upper")}

    return _cached("poiseuille_neutral", dict(Re=Re_values, N=N, kb=list(k_bounds), n_k=n_k), compute, cache)


def couette_max_growth(k, Re, N: int = 80) -> float:
    """Largest c_i of plane Couette flow (U = y on ±1) over all given (k, Re) pairs — negative: linearly stable.

    Book: §11.10 (Table 11.1 "always stable").  Expect k = 1, Re = 10³: −0.1192; Re = 10⁴: −0.0521.  Label: converged.
    """
    best = -np.inf
    for kk in np.atleast_1d(k):
        for R in np.atleast_1d(Re):
            c = orr_sommerfeld_eigs(float(kk), float(R), lambda y: y, lambda y: 0.0 * y, N=N, filter=False)
            c = c[np.abs(c) < 5]
            if c.size:
                best = max(best, float(c[0].imag))
    return best


def blasius_critical(N: int = 100, y_max: float | None = None, cache: bool = True) -> dict:
    """Critical point of the (parallel) Blasius boundary layer (cached).

    Book: §11.11 (Table 11.1, Re based on δ*; Fig. 11.26), Orr–Sommerfeld (11.79)–(11.80).  Base flow :func:`blasius_base`;
    linear map on [0, y_max]; minimum over k of the neutral Re (:func:`critical_point`).
    Parameters
    ----------
    N : Chebyshev degree in the profile's smallest box (20 δ*) [–]; default 100.
    y_max : None (default) → the wavelength-scaled box of :func:`os_box_numerics`, max(20, 12/k) δ* with degree
        N·√(y_max/20) (at k_c: 39.5 δ*, 12 e-folds of e^{−ky}, degree 141); a number [δ*] → exactly that box with degree N.
    cache : with the defaults, return the published value of reference/ch11/critical_points.json (≈ 70 s to recompute);
        otherwise the npz cache of outputs/ch11/cache.
    Returns
    -------
    dict(Re_c [–, U∞δ*/ν], k_c (αδ*) [–], omega_c (= k_c c_r, units U∞/δ*), c_r [units U∞]).
    Expect 519.060, 0.30377, 0.12049, 0.39664 (Thomas via Gallagher, Griffiths & Stephen 2016: 519.2, 0.303, 0.120; their own
    n = 1 row 519.12; Jordinson 1970: 520 — ours is 2.7e-4 below Thomas's value; not tuned).
    DEVIATION (truncation): the book has no box.  The former default, the fixed box 20 δ* (only k_c·y_max = 6 e-folds),
    gave 519.0765 — 3.2e-5 (relative) above the box-converged value; it is still available as ``y_max=20.0``.
    Validation: V3, Re_c — rule (default), N = 60 / 80 / 100 / 140: 518.9724 / 519.06016 / 519.06012 / 519.06012; fixed boxes
    (N, y_max) = (100, 20) → 519.0765, (100, 40) → 519.0558 (box right, N too small: the linear map needs N ∝ √y_max),
    (142, 40) → 519.06012, (200, 40) → 519.06012, (200, 80) → 519.06012: converged to 519.0601 (2e-7 relative);
    k_c = 0.30377 (the minimum is flat: 0.303771 – 0.303772), c_r = 0.39664, ω_c = 0.12049.  V5 as above.
    Label: converged, benchmark.
    """
    if cache and N == 100 and y_max is None and (ref := _ref_critical("blasius")):
        return ref

    def compute():
        r = critical_point(_os_lead(parallel_profile("blasius", y_max=y_max), N), (0.27, 0.34), (350.0, 1200.0))
        return dict(Re_c=r["Re_c"], k_c=r["k_c"], c_r=r["c_r"], omega_c=r["k_c"] * r["c_r"])

    box = _OS_BOX_RULE if y_max is None else float(y_max)
    d = _cached("blasius_critical", dict(N=N, y_max=box), compute, cache)
    return {k: float(v) for k, v in d.items()}


def blasius_neutral_curve(Re_values=None, in_frequency: bool = False, cache: bool = True, N: int = 80,
                          y_max: float | None = None, k_bounds: Sequence[float] = (0.04, 0.42), n_k: int = 31,
                          fast: bool = False) -> dict:
    """Neutral curve of the Blasius boundary layer in (Re_δ*, kδ*) or (Re_δ*, F = ων/U∞² = kc_r/Re_δ*) (cached).

    Book: §11.11, Fig. 11.26.  Returns dict(Re, k_lower, k_upper, c_lower, c_upper, F_lower, F_upper; with in_frequency the
    keys ``lower``/``upper`` hold F, otherwise k).
    y_max : None (default) → the wavelength-scaled box of :func:`os_box_numerics` (max(20, 12/k) δ*, N·√(y_max/20)); a
    number → exactly that box with degree N (the former default 20.0 is still available this way).  DEVIATION (truncation):
    the fixed box 20 δ* holds only 1.4 e-folds of e^{−ky} on the lower branch at Re = 6000 and put it at k = 0.0725; boxes
    40 / 80 / 160 give 0.07879 / 0.07902 / 0.07902 (Re = 3026: 0.0995 → 0.1024; Re = 1051: 0.1660 → 0.1663; Re = 530:
    unchanged to 1e-4; the upper branch moves by < 5e-5 everywhere).
    Validation: V3 box and N doubling (numbers above).  Label: converged.
    """
    if Re_values is None:
        Re_values = np.geomspace(530.0, 3000.0 if fast else 6000.0, 8 if fast else 22)
    Re_values = [float(r) for r in np.atleast_1d(Re_values)]

    def compute():
        nc = neutral_curve(_os_lead(parallel_profile("blasius", y_max=y_max), N), Re_values, k_bounds, n_k=n_k)
        out = {key: nc[key] for key in ("Re", "k_lower", "k_upper", "c_lower", "c_upper")}
        out["F_lower"] = out["k_lower"] * np.real(out["c_lower"]) / out["Re"]
        out["F_upper"] = out["k_upper"] * np.real(out["c_upper"]) / out["Re"]
        return out

    box = _OS_BOX_RULE if y_max is None else float(y_max)
    d = _cached("blasius_neutral", dict(Re=Re_values, N=N, y_max=box, kb=list(k_bounds), n_k=n_k), compute, cache)
    d["lower"], d["upper"] = (d["F_lower"], d["F_upper"]) if in_frequency else (d["k_lower"], d["k_upper"])
    return d


def falkner_skan_neutral_curve(m: float, Re_values=None, cache: bool = True, N: int = 80, y_max: float | None = None,
                               k_bounds: Sequence[float] = (0.03, 0.8), n_k: int = 35, fast: bool = False) -> dict:
    """Neutral curve of a Falkner–Skan boundary layer (favourable m > 0: the loop closes; adverse m < 0: flat upper branch).

    Book: §11.10, Fig. 11.24 (kδ* vs Re_δ*); base flows from ``core.boundary_layer.falkner_skan`` (Ch. 9 (9.36)) in δ* units.
    y_max : None (default) → the wavelength-scaled box of :func:`os_box_numerics` (max(25, 12/k) δ*, N·√(y_max/25)); a number
    → exactly that box with degree N (the former default 25.0).  DEVIATION (truncation): as :func:`blasius_neutral_curve`.
    Measured, lower branch, box 25 → scaled box (N = 80; N = 120 in brackets): m = 0.1: Re = 5000: 0.1245 → 0.1249 (0.1248),
    2e4: 0.0682 → 0.0718 (0.0719), 1e5: 0.0305 → 0.0424 (0.0425); m = −0.05: Re = 1000: 0.1164 → 0.1168 (0.1168), 5000:
    0.0536 → 0.0592 (0.0592), 2e4: none found (the box-25 flow looked unstable down to k = 0.03) → 0.0369 (0.0369).
    Not a box effect and still open: at Re ≥ 2e4 the *upper* branch depends on N (m = 0.1, Re = 1e5: 0.1266 with N = 80,
    0.1248 with N = 120; m = −0.05, Re = 1e5: 0.4801, 0.4924) — the critical layer is under-resolved there, so those upper
    branch values are good to about 3 %; pass a larger N for Re ≥ 2e4.
    Returns dict(Re, k_lower, k_upper, c_lower, c_upper, m).  Label: converged (lower branch; upper branch for Re ≤ 5000),
    qualitative (upper branch at Re ≥ 2e4 with the default N).
    """
    if Re_values is None:
        Re_values = np.geomspace(100.0, 2e4 if fast else 1e5, 8 if fast else 20)
    Re_values = [float(r) for r in np.atleast_1d(Re_values)]

    def compute():
        nc = neutral_curve(_os_lead(parallel_profile("falkner_skan", m=m, y_max=y_max), N), Re_values, k_bounds, n_k=n_k)
        out = {key: nc[key] for key in ("Re", "k_lower", "k_upper", "c_lower", "c_upper")}
        out["m"] = m
        return out

    box = _OS_BOX_RULE if y_max is None else float(y_max)
    return _cached("falkner_skan_neutral", dict(m=m, Re=Re_values, N=N, y_max=box, kb=list(k_bounds), n_k=n_k), compute,
                   cache)


def tanh_max_growth(N: int = 120, y_max: float = 30.0) -> dict:
    """Most-amplified inviscid mode of the tanh shear layer (Rayleigh, c_r = 0).

    Book: §11.10 (Fig. 11.23, inviscid limit); §11.7, p. 506.  Returns dict(k, kci, c).  Expect k = 0.4449, kc_i = 0.18970
    (Michalke 1964: 0.4446, 0.1897).  Label: converged, benchmark (approximate digits).
    """
    pr = parallel_profile("tanh", y_max=y_max)
    r = max_growth(lambda k: rayleigh_eigs(k, pr["U"], pr["Upp"], N=N, bc="decay", y_max=y_max,
                                           map_scale=pr["map_scale"], unstable_only=True), (0.2, 0.8))
    return dict(k=r["k"], kci=r["growth"], c=r["c"])


def tanh_shear_layer_neutral_curve(Re_values=None, N: int = 100, cache: bool = True, y_max: float = 40.0,
                                   fast: bool = False) -> dict:
    """Upper neutral wavenumber k_u(Re) of the tanh shear layer (unstable for 0 < k < k_u; k_u → 1 as Re → ∞; Re_c = 0).

    Book: §11.10, Fig. 11.23 (kL vs Re = U₀L/ν).  OS on (−y_max, y_max), tan map; Brent in k on the leading c_i.  Tested:
    y_max = 40 and 60 agree to 1e-8 for Re ≥ 2 (k_u(2) = 0.2732, k_u(10) = 0.6544, k_u(200) = 0.9702).
    Returns dict(Re, k_upper).  Label: converged.
    """
    if Re_values is None:
        Re_values = np.geomspace(2.0, 400.0, 8 if fast else 18)
    Re_values = [float(r) for r in np.atleast_1d(Re_values)]
    fn = _os_lead(parallel_profile("tanh", y_max=y_max), N)

    def compute():
        ku = []
        for R in Re_values:
            ks = np.linspace(0.02, 1.05, 30)
            v = np.array([fn(k, R)[0].imag for k in ks])
            j = [i for i in range(len(ks) - 1) if v[i] > 0 >= v[i + 1]]
            ku.append(brentq(lambda k: fn(k, R)[0].imag, ks[j[-1]], ks[j[-1] + 1], xtol=1e-9) if j else np.nan)  # noqa: B023
        return dict(Re=np.array(Re_values), k_upper=np.array(ku))

    return _cached("tanh_neutral", dict(Re=Re_values, N=N, y_max=y_max), compute, cache)


def bickley_critical(parity: str = "sinuous", cache: bool = True, N: int = 80, y_max: float = 60.0) -> dict:
    """Critical point of the Bickley jet U = sech²y (Orr–Sommerfeld; sinuous = φ even, varicose = φ odd).

    Book: §11.10 (Table 11.1 jet).  Returns dict(Re_c, k_c, c_r).  Ours: sinuous 4.017 at k = 0.1728, c_r = 0.0433
    (y_max = 60 and 80 agree to 1e-6).  Validation: V5 Tatsumi & Kakutani (1958) ≈ 4.0 at k ≈ 0.2 (approximate); V1 inviscid
    neutral modes k = 2 (sinuous), k = 1 (varicose), c = 2/3.  With ``cache`` and the defaults (sinuous) the published value
    in reference/ch11/critical_points.json is used.  Label: converged, benchmark (approximate).
    """
    if cache and parity == "sinuous" and N == 80 and y_max == 60.0 and (ref := _ref_critical("bickley_sinuous")):
        return ref
    par = {"sinuous": "even", "varicose": "odd"}[parity]

    def compute():
        r = critical_point(_os_lead(parallel_profile("bickley", y_max=y_max), N, parity=par), (0.08, 0.32), (2.0, 15.0), n_k=7)
        return dict(Re_c=r["Re_c"], k_c=r["k_c"], c_r=r["c_r"])

    d = _cached("bickley_critical", dict(parity=parity, N=N, y_max=y_max), compute, cache)
    return {k: float(v) for k, v in d.items()}


def _os_track(profile: dict, N: int, Re: float, k0: float, c0: complex, k1: float, n: int = 24,
              n_efold: float = 12.0) -> complex:
    """Follow one Orr–Sommerfeld eigenvalue from (k0, c0) to k1 at fixed Re by continuation: n equal steps in k, at each the
    eigenvalue nearest to the linear extrapolation of the last two.  Returns c at k1 (ours; used to tell mode families apart)."""
    c_prev, c = None, complex(c0)
    for k in np.linspace(float(k0), float(k1), int(n) + 1)[1:]:
        nm = os_box_numerics(profile, float(k), N, n_efold)
        w = orr_sommerfeld_eigs(float(k), Re, profile["U"], profile["Upp"], domain=profile["domain"], N=nm["N"],
                                bc=profile["bc"], y_max=nm["y_max"], map_scale=nm["map_scale"], filter=False,
                                parity=profile.get("parity"))
        guess = c if c_prev is None else 2.0 * c - c_prev
        c_prev, c = c, complex(w[np.argmin(np.abs(w - guess))])
    return c


def bickley_neutral_curve(Re_values=None, N: int = 80, cache: bool = True, k_bounds: Sequence[float] = (0.02, 2.0),
                          n_k_long: int = 30, n_k: int = 26, y_max: float | None = None, fast: bool = False,
                          n_efold: float = 12.0) -> dict:
    """Neutral wavenumbers of the sinuous Bickley jet U = sech²y at each Re: the main unstable band and the long-wave band.

    Book: §11.10 (Table 11.1 jet; a neutral curve c_i(k, Re) = 0 as in Figs. 11.23–11.24), Orr–Sommerfeld (11.79)–(11.80)
    with φ even.  Method (ours): at each Re the sign of the leading c_i on a k grid (``n_k_long`` log-spaced points from
    k_bounds[0] to 0.3, then ``n_k`` equally spaced to k_bounds[1]) and Brent on every sign change; the roots are sorted
    into bands by the direction of the crossing.  The stability boundary needs no mode identification (the flow is unstable
    where *any* mode grows; box eigenvalues of the continuous spectrum have c_i < −k/Re < 0 and cannot fake a crossing);
    the family of each neutral mode is found afterwards by continuation in k (:func:`_os_track`).
    DEVIATION (truncation): box max(40, 12/k) by default (:func:`os_box_numerics`); ``y_max`` = a number forces that box.
    What the converged computation shows (and the fixed boxes 40 / 60 hid — there c_i at k < 0.1 was wrong, e.g.
    Re = 14.58, k = 0.04: −0.0184 (box 40), −0.0026 (60), +0.00125 (150 … 1200)):
      * up to Re ≈ 17.5 one unstable band; its lower edge falls steeply (k = 0.1257 at Re = 4.1, 0.0655 at 4.72, 0.0427 at
        5.44, 0.0296 at 6.26, 0.0211 at 7.2) and leaves the resolved range k ≥ 0.02 just above Re = 7.2: for larger Re
        ``k_lower`` is NaN — unstable down to the smallest k computed, lower edge unknown;
      * between Re = 17.4 and 17.6, at k ≈ 0.068, a stable gap opens inside the band and splits it: the main band above
        (lower edge 0.0737 at Re = 17.6, 0.0758 at 19.3, then ≈ 1.6/Re: 0.0620 at 25.6, 0.0234 at 68.7; c_r from −0.03 to
        −0.22; below k = 0.02 beyond Re ≈ 80) and a long-wave band below it (upper edge ``k_long_upper`` = 0.0624 at
        Re = 17.6, 0.0447 at 19.3, 0.0222 at 25.6 — roughly (7/Re)³; c_r from −0.014 to −0.001; below k = 0.02 beyond
        Re ≈ 26, where nothing is reported).  Up to Re = 18.6 both edges of the gap are neutral points of one mode (followed
        continuously in k); from Re = 19.0 on they belong to two different modes (``gap_same_mode``) — the "kink" of the old
        table was this gap, drawn with wrong numbers on its left.
      * the critical point (:func:`bickley_critical`, 4.017 at k = 0.173, 10 e-folds in its box 60) is not affected: at
        Re = 4.0 every k from 0.02 to 0.23 decays.
    Parameters
    ----------
    Re_values : Reynolds numbers U₀L/ν [–] (default 40 log-spaced 4.1 … 1000; fast: 8).   N : Chebyshev degree [–].
    cache : npz cache in outputs/ch11/cache.   k_bounds : scanned wavenumbers [1/L]; below 0.02 the box 12/k > 600 and
    the growth c_i ~ 1e-5 of the long-wave mode is no longer resolved (noise 2e-6) — nothing is reported there.
    n_k_long, n_k : grid sizes [–] (fast: 14, 12).   y_max : None or an explicit box [L].   fast : coarser k grid.
    Returns
    -------
    dict(Re, k_upper, c_upper, k_lower, c_lower (main band; k_lower NaN where the band reaches below k_bounds[0]),
    k_long_upper, c_long_upper, k_long_lower, c_long_lower (the long-wave band below the gap; NaN where there is no gap or
    the edge is below k_bounds[0]), gap_same_mode (1.0: the two edges of the gap are the same mode; 0.0: different modes;
    NaN: no gap), unstable_at_k_min (1.0 where c_i > 0 at k_bounds[0]), k_min).  k non-dimensional (scaled by 1/L), c by U₀.
    Assumptions: parallel base flow, 2-D sinuous disturbances.
    Validation: V3 the 40-point table recomputed with (N, n_efold) = (80, 12), (160, 12), (80, 24), (160, 24): every neutral
    k agrees to ≤ 5e-6 for k ≥ 0.03 and to ≤ 8e-6 at k ≈ 0.021–0.022 (4e-4 relative: 0.022197, 0.022189, 0.022194,
    0.022196 at Re = 25.6), so the 4th significant figure of the two entries nearest k = 0.02 is uncertain by one unit.
    ``n_efold`` : e-folds of e^{−k|y|} in the box (default 12).  Label: converged.
    """
    if Re_values is None:
        Re_values = np.geomspace(4.1, 1000.0, 8 if fast else 40)
    Re_values = [float(r) for r in np.atleast_1d(Re_values)]
    k_lo, k_hi = float(k_bounds[0]), float(k_bounds[1])
    nl, nu = (14, 12) if fast else (int(n_k_long), int(n_k))
    k_mid = min(0.3, k_hi)
    ks = np.unique(np.concatenate([np.geomspace(k_lo, k_mid, nl), np.linspace(k_mid, k_hi, nu)])) if k_lo < k_mid \
        else np.linspace(k_lo, k_hi, nu)
    pr = parallel_profile("bickley", parity="even", y_max=y_max)
    fn = _os_lead(pr, N, n_efold=n_efold)
    nan = complex(np.nan, np.nan)

    def compute():
        keys = ("k_upper", "k_lower", "k_long_upper", "k_long_lower")
        out = {key: [] for key in keys}
        out.update({"c" + key[1:]: [] for key in keys})
        out.update(gap_same_mode=[], unstable_at_k_min=[])
        for Re in Re_values:
            v = np.array([fn(k, Re)[0].imag for k in ks])
            rising, falling = [], []
            for j in range(len(ks) - 1):
                if np.isfinite(v[j]) and np.isfinite(v[j + 1]) and v[j] * v[j + 1] < 0:
                    r = brentq(lambda k: fn(k, Re)[0].imag, ks[j], ks[j + 1], xtol=1e-9)  # noqa: B023
                    (rising if v[j] < 0 else falling).append(r)
            ku = falling[-1] if falling else np.nan  # c_i turns negative for good: upper branch
            below = [r for r in rising if not np.isfinite(ku) or r < ku]
            kl = below[-1] if below else np.nan  # lower edge of the band that ends at k_upper
            gap = [r for r in falling if np.isfinite(kl) and r < kl]
            kg = gap[-1] if gap else np.nan  # upper edge of the long-wave band (lower edge of the stable gap)
            outer = [r for r in rising if np.isfinite(kg) and r < kg]
            ko = outer[-1] if outer else np.nan
            vals = dict(k_upper=ku, k_lower=kl, k_long_upper=kg, k_long_lower=ko)
            for key in keys:
                out[key].append(vals[key])
                out["c" + key[1:]].append(complex(fn(vals[key], Re)[0]) if np.isfinite(vals[key]) else nan)
            same = np.nan
            if np.isfinite(kg):
                c_end = _os_track(pr, N, Re, kg, out["c_long_upper"][-1], kl, n_efold=n_efold)
                same = float(abs(c_end - out["c_lower"][-1]) < 2e-3)
            out["gap_same_mode"].append(same)
            out["unstable_at_k_min"].append(float(v[0] > 0))
        res = {key: np.array(val) for key, val in out.items()}
        res.update(Re=np.array(Re_values), k_min=k_lo)
        return res

    box = f"{_OS_BOX_RULE}; n_efold = {float(n_efold):g}" if y_max is None else float(y_max)
    return _cached("bickley_neutral", dict(Re=Re_values, N=N, box=box, k=[float(k) for k in ks]), compute, cache)


def table_11_1(cache: bool = True) -> list:
    """Our recomputation of Table 11.1 (critical Reynolds numbers) with the published benchmarks (the book's rounded column
    stays private).

    Book: §11.10, Table 11.1.  Rows: jet (:func:`bickley_critical`), shear layer (Re_c = 0, :func:`tanh_shear_layer_neutral_curve`),
    Blasius (:func:`blasius_critical`), plane Poiseuille (:func:`poiseuille_critical`), pipe (stated: linearly stable, not
    computed — cylindrical OS not coded), plane Couette (:func:`couette_max_growth` < 0).
    Returns list of dict(flow, U, Re_c_ours, k_c_ours, benchmark, source, length_scale, remark).  Label: converged, benchmark.
    """
    P = poiseuille_critical(cache=cache)
    Bl = blasius_critical(cache=cache)
    Bj = bickley_critical(cache=cache)
    cc = couette_max_growth([0.5, 1.0, 2.0], [1e3, 1e4, 1e5])
    return [
        dict(flow="jet (Bickley)", U="sech²(y/L)", Re_c_ours=Bj["Re_c"], k_c_ours=Bj["k_c"], benchmark=4.0,
             source="Tatsumi & Kakutani (1958)", length_scale="L (jet width scale)", remark="sinuous mode; inflectional"),
        dict(flow="shear layer", U="tanh(y/L)", Re_c_ours=0.0, k_c_ours=0.0, benchmark=0.0, source="Betchov & Szewczyk (1963)",
             length_scale="L", remark="unstable at every Re for small enough k"),
        dict(flow="Blasius", U="f′(η)", Re_c_ours=Bl["Re_c"], k_c_ours=Bl["k_c"], benchmark=519.2,
             source="Thomas, via Gallagher, Griffiths & Stephen (2016); Jordinson (1970): 520", length_scale="δ*",
             remark="parallel-flow approximation"),
        dict(flow="plane Poiseuille", U="1 − (y/L)²", Re_c_ours=P["Re_c"], k_c_ours=P["k_c"], benchmark=5772.22,
             source="Orszag (1971)", length_scale="L = half-width", remark="TS waves (viscous instability)"),
        dict(flow="pipe", U="1 − (r/R)²", Re_c_ours=math.inf, k_c_ours=float("nan"), benchmark=math.inf,
             source="not computed here", length_scale="R", remark="linearly stable; transition is finite-amplitude"),
        dict(flow="plane Couette", U="y/L", Re_c_ours=math.inf, k_c_ours=float("nan"), benchmark=math.inf,
             source="Romanov (1973)", length_scale="L = half-gap", remark=f"max c_i on our (k, Re) grid = {cc:.3g} < 0"),
    ]


def ts_wave_fields(x, y, t: float, k: float, c: complex, phi, phi_y, amp: float = 0.05, U=None) -> dict:
    """Real fields of an Orr–Sommerfeld (Tollmien–Schlichting) mode on the (y, x) grid (len(y) × len(x)).

    Book: §11.8, p. 510: [u, v, ψ] = [û, v̂, φ]e^{ik(x−ct)}, û = dφ/dy, v̂ = −ikφ; u = ∂ψ/∂y, v = −∂ψ/∂x (the §11.8 sign).
    Parameters: x (1-D); y (1-D nodes of phi); t; k; c; phi (complex on y); phi_y (dφ/dy on y, e.g. ``os_mode(...)["u_hat"]``);
    amp; U (callable or array: base flow, gives u_total = U + u).
    Returns dict(psi, u, v, u_total (None without U), uv_mean (x-average of uv on y)).  Label: analytic.
    """
    x, y = _F(x), _F(y)
    phi, phi_y = np.asarray(phi, dtype=complex), np.asarray(phi_y, dtype=complex)
    E = amp * np.exp(1j * k * (x[None, :] - c * t))
    u = np.real(phi_y[:, None] * E)
    out = dict(psi=np.real(phi[:, None] * E), u=u, v=np.real(-1j * k * phi[:, None] * E), u_total=None,
               uv_mean=0.5 * np.real(phi_y * np.conj(-1j * k * phi)) * abs(amp) ** 2 * math.exp(2 * k * complex(c).imag * t))
    if U is not None:
        Uy = _F(U(y) if callable(U) else U)
        out["u_total"] = Uy[:, None] + u
    return out


def ts_mode(flow: str = "poiseuille", k: float = 1.0, Re: float = 1e4, N: int = 100) -> dict:
    """Leading Orr–Sommerfeld mode of a named flow (:func:`parallel_profile`) with its energy budget (A4, E8).

    Book: §11.8–§11.11, (11.79)–(11.80), (11.88).  Returns :func:`core.stability.os_mode`'s dict plus ``budget`` and
    ``profile``.  Example: Poiseuille Re = 10⁴, k = 1 → c = 0.23752649 + 0.00373967i, P/Λ = 1.616.
    Unbounded and semi-infinite flows ("tanh", "bickley" (all parities, as before), "blasius"): solved in the
    wavelength-scaled box of :func:`os_box_numerics` (N is the degree in the profile's smallest box) — the same numerics as
    the published grids os_grid_*.csv; the returned mode is the leading eigenvalue of that box.  Where the flow is stable
    this can be an eigenvalue of the discretised continuous spectrum, not a mode (``far`` = :func:`far_field_fraction` of φ
    is then ~1, and :func:`os_leading_mode` reports ``discrete=False``).  Channel flows are unchanged.  Label: converged.
    """
    pr = parallel_profile(flow)
    nm = os_box_numerics(pr, k, N)
    m = os_mode(k, Re, pr["U"], pr["Upp"], domain=pr["domain"], N=nm["N"], bc=pr["bc"], y_max=nm["y_max"],
                map_scale=nm["map_scale"], parity=pr.get("parity"))
    m["budget"] = disturbance_energy_budget(k, m["c"], m["phi"], m["y"], pr["Up"], Re, grid=m["grid"])
    m["profile"] = pr
    m["far"] = 0.0 if pr["bc"] == "wall" else far_field_fraction(
        m["phi"], m["y"], nm["y_max"], centre=None if pr["bc"] == "semi_infinite" else 0.5 * sum(pr["domain"]))
    return m


_OS_TABLE_FLOWS = ("poiseuille", "blasius", "tanh", "bickley")
_OS_HDR = "# ours (fluidpy.ch11_instability.neutral_curve_tables), Orr–Sommerfeld, Chebyshev; not book data\n"
_OS_GRID_NOTE = ("# leading mode at each (Re, k): c_r, c_i, production P, dissipation Lambda, energy E (mode normalised "
                 "max|u_hat| = 1)")
_OS_GRID_NOTE_OPEN = ("; box {box}; the mode is the least-damped eigenvalue whose eigenfunction is localised (|phi| beyond half "
                      "the box < 20 % of its maximum) and which is unchanged (1e-5) in a box 1.5 times as large (searched again with twice the degree if none is found) - box-dependent "
                      "eigenvalues of the continuous spectrum are never listed; where c_i < -k/Re the mode is more damped than "
                      "the edge of the continuous spectrum; rows with P = nan: no such mode found - c_i is then that edge, -k/Re "
                      "(analytic, not a mode), and c_r its free-stream speed ({cr})")
_OS_NEUTRAL_NOTE = dict(
    blasius="; box max(20, 12/k) delta* (the fixed box 20 put the lower branch up to 8 % too low at Re = 6000)",
    bickley="; sinuous mode, box max(40, 12/k), k scanned from {k_min:g} to {k_max:g}; k_lower, k_upper = edges of the "
            "unstable band that contains the fastest-growing wave; k_lower = nan with k_upper finite: that band reaches "
            "below k = {k_min:g} (its lower edge is not resolved - the flow is UNSTABLE there, not stable){gap}",
)


def _write_bickley_neutral(out: Path, nc: dict, label: str) -> dict:
    """Write os_neutral_bickley.csv (same five columns as the other flows) and os_neutral_bickley_longwave.csv."""
    Re, kg = np.asarray(nc["Re"]), np.asarray(nc["k_long_upper"])
    has = np.isfinite(kg)
    gap = ""
    if has.any():
        i0 = int(np.argmax(has))
        left = f"between Re = {Re[i0 - 1]:.4g} and {Re[i0]:.4g}" if i0 > 0 else f"below Re = {Re[i0]:.4g}"
        gap = (f"; a stable gap opens inside the band {left}: from there on k_lower is the upper edge of the gap (it "
               "jumps - do not join it across) and a second, long-wave unstable band lies below the gap "
               "(os_neutral_bickley_longwave.csv)")
    note = _OS_NEUTRAL_NOTE["bickley"].format(k_min=float(nc["k_min"]), k_max=float(nc["k_max"]), gap=gap)
    p = out / "os_neutral_bickley.csv"
    with p.open("w", encoding="utf-8") as f:
        f.write(_OS_HDR + f"# flow: {label}{note}\nRe,k_lower,k_upper,cr_lower,cr_upper\n")
        for R, kl, ku, cl, cu in zip(Re, nc["k_lower"], nc["k_upper"], nc["c_lower"], nc["c_upper"]):
            f.write(f"{R:.4g},{kl:.4g},{ku:.4g},{np.real(cl):.4g},{np.real(cu):.4g}\n")
    p2 = out / "os_neutral_bickley_longwave.csv"
    with p2.open("w", encoding="utf-8") as f:
        f.write(_OS_HDR + f"# flow: {label}; sinuous mode, box max(40, 12/k): the long-wave unstable band below the stable "
                f"gap, k_long_lower < k < k_long_upper (the gap is k_long_upper < k < k_lower of os_neutral_bickley.csv); "
                f"nan = no gap at this Re, or that edge lies below the smallest k computed ({float(nc['k_min']):g}) and is "
                "unknown; same_mode = 1: both edges of the gap are neutral points of one mode, 0: of two different modes\n"
                "Re,k_long_lower,k_long_upper,cr_long_upper,same_mode\n")
        for R, ko, kgi, cg, sm in zip(Re, nc["k_long_lower"], kg, nc["c_long_upper"], nc["gap_same_mode"]):
            f.write(f"{R:.4g},{ko:.4g},{kgi:.4g},{np.real(cg):.4g},{sm:.4g}\n")
    return {"os_neutral_bickley": str(p), "os_neutral_bickley_longwave": str(p2)}


def neutral_curve_tables(out_dir: str | Path | None = None, fast: bool = False, n_Re: int = 40, cache: bool = True,
                         flows: Sequence[str] | None = None, write_modes: bool = True) -> dict:
    """Write our Orr–Sommerfeld tables for explainer E8 / figure F7 to ``reference/ch11/`` (≤ 4 s.f., labelled "ours").

    Book: §11.8–§11.11.  Part C.5 5.8: neutral curves on ``n_Re`` log-spaced Re for Poiseuille, Blasius, tanh, Bickley
    (sinuous); leading c and the budget (P, Λ, E) on a 24 Re × 25 k grid per flow (fast: 8 × 8); mode φ, û, v̂ on 41 y points
    for the presets (Poiseuille Re = 10⁴, k = 1; Poiseuille at its critical point; Blasius Re = 1000, k = 0.25; tanh Re = 50,
    k = 0.45; Bickley Re = 20, k = 0.3).
    flows : subset of ("poiseuille", "blasius", "tanh", "bickley") to (re)write; None = all.   write_modes : write
    os_modes.json (default True).
    Truncation (DEVIATION, ours): Blasius, tanh and Bickley are solved in the wavelength-scaled box of
    :func:`os_box_numerics`; the fixed boxes used before (20, 30, 40) were too short for long waves — the Bickley lower
    branch, the Blasius lower branch above Re ≈ 1500 and the small-k columns of the three grids were wrong.  The Bickley
    neutral table comes from :func:`bickley_neutral_curve` (two unstable bands; the second one in
    os_neutral_bickley_longwave.csv).  In the grids of these three flows a point where no discrete mode lies above the
    continuous spectrum (:func:`os_leading_mode`) carries the analytic edge c_i = −k/Re and NaN for P, Λ, E — never a
    box-dependent eigenvalue.  Whatever is not resolved is NaN and the file header says what a NaN means.
    Returns dict(name → path).  Label: converged.
    """
    out = Path(out_dir) if out_dir is not None else _root() / "reference" / "ch11"
    out.mkdir(parents=True, exist_ok=True)
    paths: dict = {}
    nr = 8 if fast else n_Re
    spec_all = dict(
        poiseuille=dict(Re=np.geomspace(5800.0, 1e6, nr), kb=(0.15, 1.35), grid_Re=np.geomspace(2000.0, 1e6, 8 if fast else 24),
                        grid_k=np.linspace(0.2, 1.4, 8 if fast else 25), N=80),
        blasius=dict(Re=np.geomspace(530.0, 6000.0, nr), kb=(0.04, 0.42), grid_Re=np.geomspace(200.0, 6000.0, 8 if fast else 24),
                     grid_k=np.linspace(0.05, 0.45, 8 if fast else 25), N=80),
        tanh=dict(Re=np.geomspace(2.0, 400.0, nr), kb=(0.02, 1.05), grid_Re=np.geomspace(1.0, 400.0, 8 if fast else 24),
                  grid_k=np.linspace(0.05, 1.2, 8 if fast else 25), N=100),
        bickley=dict(Re=np.geomspace(4.1, 1000.0, nr), kb=(0.02, 2.0), grid_Re=np.geomspace(2.0, 1000.0, 8 if fast else 24),
                     grid_k=np.linspace(0.05, 1.8, 8 if fast else 25), N=80),
    )
    names = list(_OS_TABLE_FLOWS) if flows is None else [_PROFILE_ALIASES.get(n, n) for n in flows]
    for name in names:
        if name not in spec_all:
            raise ValueError(f"neutral_curve_tables: unknown flow {name!r} (choose from {_OS_TABLE_FLOWS})")
    for name in names:
        spec = spec_all[name]
        pr = parallel_profile(name, parity="even") if name == "bickley" else parallel_profile(name)
        wall = pr["bc"] == "wall"
        label = f"{pr['label']}; L = {pr['length']}, U0 = {pr['velocity']}"
        fn = _os_lead(pr, spec["N"])
        Rs = [float(r) for r in spec["Re"]]
        if name == "bickley":
            nc = dict(bickley_neutral_curve(Rs, N=spec["N"], cache=cache, k_bounds=spec["kb"], fast=fast))
            nc["k_max"] = spec["kb"][1]
            paths.update(_write_bickley_neutral(out, nc, label))
        else:
            def nc_compute(fn=fn, Rs=Rs, spec=spec):
                nc = neutral_curve(fn, Rs, spec["kb"], n_k=31)
                return {key: nc[key] for key in ("Re", "k_lower", "k_upper", "c_lower", "c_upper")}

            key = dict(Re=Rs, kb=list(spec["kb"]), N=spec["N"])
            if not wall:
                key["box"] = _OS_BOX_RULE
            nc = _cached(f"os_neutral_{name}", key, nc_compute, cache)
            p = out / f"os_neutral_{name}.csv"
            with p.open("w", encoding="utf-8") as f:
                f.write(_OS_HDR + f"# flow: {label}{_OS_NEUTRAL_NOTE.get(name, '')}\nRe,k_lower,k_upper,cr_lower,cr_upper\n")
                for R, kl, ku, cl, cu in zip(nc["Re"], nc["k_lower"], nc["k_upper"], nc["c_lower"], nc["c_upper"]):
                    f.write(f"{R:.4g},{kl:.4g},{ku:.4g},{np.real(cl):.4g},{np.real(cu):.4g}\n")
            paths[f"os_neutral_{name}"] = str(p)

        def grid_compute(pr=pr, spec=spec, wall=wall):
            GR, GK = spec["grid_Re"], spec["grid_k"]
            cr = np.full((len(GR), len(GK)), np.nan)
            ci, P, D, E = cr.copy(), cr.copy(), cr.copy(), cr.copy()
            for i, R in enumerate(GR):
                for j, k in enumerate(GK):
                    if wall:
                        try:
                            m = os_mode(float(k), float(R), pr["U"], pr["Upp"], domain=pr["domain"], N=spec["N"], bc=pr["bc"],
                                        y_max=pr.get("y_max"), map_scale=pr.get("map_scale"), parity=pr.get("parity"))
                        except RuntimeError:
                            continue
                        c = m["c"]
                        b = disturbance_energy_budget(float(k), c, m["phi"], m["y"], pr["Up"], float(R), grid=m["grid"])
                    else:
                        try:
                            lm = os_leading_mode(pr, float(k), float(R), spec["N"])
                        except RuntimeError:
                            continue
                        c, b = lm["c"], lm["budget"]
                    cr[i, j], ci[i, j] = c.real, c.imag
                    if b is not None:
                        P[i, j], D[i, j], E[i, j] = b["production"], b["dissipation"], b["E"]
            return dict(Re=GR, k=GK, cr=cr, ci=ci, P=P, D=D, E=E)

        key = dict(Re=list(spec["grid_Re"]), k=list(spec["grid_k"]), N=spec["N"])
        note = _OS_GRID_NOTE
        if not wall:
            key.update(box=_OS_BOX_RULE, select="localised (far < 0.2) and box-independent (1e-5), all candidates, one retry at 2N")
            ym0 = pr["y_max"]
            note += _OS_GRID_NOTE_OPEN.format(
                box=f"max({ym0:g}, 12/k)", cr="nan: the two streams differ, -1 and +1" if name == "tanh" else
                ("1" if pr["bc"] == "semi_infinite" else "0"))
        gd = _cached(f"os_grid_{name}", key, grid_compute, cache)
        p = out / f"os_grid_{name}.csv"
        with p.open("w", encoding="utf-8") as f:
            f.write(_OS_HDR + note + "\nRe,k,c_r,c_i,P,Lambda,E\n")
            for i, R in enumerate(gd["Re"]):
                for j, k in enumerate(gd["k"]):
                    f.write(f"{R:.4g},{k:.4g},{gd['cr'][i, j]:.4g},{gd['ci'][i, j]:.4g},{gd['P'][i, j]:.4g},"
                            f"{gd['D'][i, j]:.4g},{gd['E'][i, j]:.4g}\n")
        paths[f"os_grid_{name}"] = str(p)
    if not write_modes:
        return paths
    Pc = poiseuille_critical(cache=cache)
    presets = [("poiseuille", 1.0, 1e4), ("poiseuille", Pc["k_c"], Pc["Re_c"]), ("blasius", 0.25, 1000.0), ("tanh", 0.45, 50.0),
               ("bickley", 0.3, 20.0)]
    modes = []
    for name, k, R in presets:
        pr = parallel_profile(name, parity="even") if name == "bickley" else parallel_profile(name)
        nm = os_box_numerics(pr, k, 100 if name != "bickley" else 80)
        m = os_mode(k, R, pr["U"], pr["Upp"], domain=pr["domain"], N=nm["N"], bc=pr["bc"],
                    y_max=nm["y_max"], map_scale=nm["map_scale"], parity=pr.get("parity"))
        b = disturbance_energy_budget(k, m["c"], m["phi"], m["y"], pr["Up"], R, grid=m["grid"])
        lo, hi = (pr["domain"] if pr["bc"] == "wall" else ((0.0, 8.0) if pr["bc"] == "semi_infinite" else (-6.0, 6.0)))
        ys = np.linspace(lo, hi, 41)
        f_re = BarycentricInterpolator(m["y"], m["phi"].real) if pr["bc"] == "wall" else None
        interp = (lambda arr, f_re=f_re: np.interp(ys, m["y"][::-1], arr[::-1])) if f_re is None else (
            lambda arr: BarycentricInterpolator(m["y"], arr)(ys))
        modes.append(dict(flow=name, k=float(f"{k:.4g}"), Re=float(f"{R:.6g}"), c=[float(f"{m['c'].real:.6g}"),
                                                                                   float(f"{m['c'].imag:.4g}")],
                          P=float(f"{b['production']:.4g}"), Lambda=float(f"{b['dissipation']:.4g}"), E=float(f"{b['E']:.4g}"),
                          y=[float(f"{v:.4g}") for v in ys],
                          phi_r=[float(f"{v:.4g}") for v in interp(m["phi"].real)],
                          phi_i=[float(f"{v:.4g}") for v in interp(m["phi"].imag)],
                          u_r=[float(f"{v:.4g}") for v in interp(m["u_hat"].real)],
                          u_i=[float(f"{v:.4g}") for v in interp(m["u_hat"].imag)],
                          v_r=[float(f"{v:.4g}") for v in interp(m["v_hat"].real)],
                          v_i=[float(f"{v:.4g}") for v in interp(m["v_hat"].imag)]))
    p = out / "os_modes.json"
    p.write_text(json.dumps(dict(note="ours (fluidpy.ch11_instability.neutral_curve_tables): Orr–Sommerfeld modes, "
                                 "max|u_hat| = 1; not book data", modes=modes), indent=1), encoding="utf-8")
    paths["os_modes"] = str(p)
    return paths


# ======================================================================================================================
# §11.11 Tollmien's approximate Blasius profile
# ======================================================================================================================
def tollmien_coefficients(eta1: float = 0.2) -> dict:
    """Coefficients of U/U∞ = aη (η ≤ η₁), 1 − b(1 − η)² (η₁ ≤ η ≤ 1), 1 (η ≥ 1) from continuity of value and slope at η₁:
    a = 2/(1 + η₁), b = 1/(1 − η₁²).  Default η₁ = 0.2 (ours; the book's breakpoint is private) → a = 1.6667, b = 1.0417.
    Label: analytic.
    """
    return dict(a=2.0 / (1.0 + eta1), b=1.0 / (1.0 - eta1 ** 2), eta1=float(eta1))


def tollmien_profile(eta, eta1: float = 0.2, printed: bool = False, a: float | None = None, b: float | None = None):
    """Tollmien's piecewise approximation of the Blasius profile (first TS calculations; zero curvature at the wall).

    Book: §11.11, p. 520: U/U∞ = aη (0 ≤ η ≤ η₁), middle branch (η₁ ≤ η ≤ 1), 1 (η ≥ 1), η = y/δ.  Slip S6: the page prints the
    middle branch as 1 − b[1 − η²] (discontinuous at η₁); the corrected 1 − b(1 − η)² matches the inner branch in value and
    slope and reaches 1 with zero slope at η = 1.  ``printed=True`` uses the printed middle branch.
    DEVIATION (design C.2 7.19, CHG): a, b from continuity with our η₁ (default 0.2); the book's coefficients stay private (pass
    them via eta1, a, b in the private V6 test).  Returns U/U∞.  Validation: V1 continuity of value and slope (corrected);
    printed fails continuity.  Label: analytic.
    """
    co = tollmien_coefficients(eta1)
    a = co["a"] if a is None else a
    b = co["b"] if b is None else b
    e = _F(eta)
    mid = 1.0 - b * (1.0 - e ** 2) if printed else 1.0 - b * (1.0 - e) ** 2  # p. 520 (printed) / corrected (S6)
    return _S(np.where(e <= eta1, a * e, np.where(e <= 1.0, mid, 1.0)))


# ======================================================================================================================
# §11.14 deterministic chaos
# ======================================================================================================================
def pendulum_rhs(t, s, g_over_l: float = 1.0, damping: float = 0.0):
    """Pendulum as a first-order system Ẋ = Y, Ẏ = −(g/l) sin X (− γY with damping, ours).

    Book: §11.14, Eq. (11.89) (from Ẍ + (g/l) sin X = 0, p. 525).  Returns [Ẋ, Ẏ].  Label: analytic.
    """
    X, Y = s[0], s[1]
    return [Y, -g_over_l * np.sin(X) - damping * Y]  # Eq. (11.89)


def pendulum_energy(s, g_over_l: float = 1.0):
    """E = ½Y² − (g/l)cos X — conserved by (11.89) without damping (to 1e-9 with rtol 1e-10).  Label: analytic."""
    return _S(0.5 * _F(s[1]) ** 2 - g_over_l * np.cos(_F(s[0])))


def phase_portrait(rhs: Callable, starts, t_end: float = 20.0, n: int = 400, **kw) -> list:
    """Trajectories of a 2-D autonomous system from several starts (a phase portrait, Fig. 11.28).

    Book: §11.14.  Parameters: rhs(t, s, **kw); starts (iterable of (X, Y)); t_end; n; kw (passed to rhs, e.g. g_over_l, mu).
    Returns list of dict(t, X, Y).  Label: converged (DOP853, rtol 1e-10).
    """
    out = []
    for s0 in starts:
        sol = solve_ivp(lambda t, s: rhs(t, s, **kw), (0.0, t_end), list(s0), t_eval=np.linspace(0.0, t_end, n),
                        rtol=1e-10, atol=1e-12, method="DOP853")
        out.append(dict(t=sol.t, X=sol.y[0], Y=sol.y[1]))
    return out


def hopf_normal_form(t, s, mu: float, omega: float = 1.0):
    """Supercritical Hopf normal form ż = (μ + iω)z − |z|²z in real form (ours: draws Fig. 11.28's attractors).

    Book: §11.14, Fig. 11.28 (fixed point → repeller with a stable limit cycle growing with R − R_cr).  Returns [ẋ, ẏ].
    Label: analytic.
    """
    x, y = s[0], s[1]
    r2 = x * x + y * y
    return [mu * x - omega * y - r2 * x, omega * x + mu * y - r2 * y]


def limit_cycle_amplitude(mu):
    """Radius √μ of the Hopf limit cycle (0 for μ ≤ 0) — Fig. 11.28c branches ±√μ.  Label: analytic."""
    return _S(np.sqrt(np.maximum(_F(mu), 0.0)))


def lorenz_rhs(t, s, Pr: float = 10.0, r: float = 28.0, b: float = 8.0 / 3.0):
    """Lorenz's model: Ẋ = Pr(Y − X), Ẏ = −XZ + rX − Y, Ż = XY − bZ.

    Book: §11.14, Eq. (11.91), r = Ra/Ra_cr, b = 4π²/(π² + k²) (Lorenz 1963: Pr = 10, b = 8/3, r = 28).  Example:
    rhs(0, (1, 1, 1)) = (0, 26, −1.667).  Label: analytic.
    """
    X, Y, Z = s[0], s[1], s[2]
    return [Pr * (Y - X), -X * Z + r * X - Y, X * Y - b * Z]  # Eq. (11.91)


def lorenz_integrate(s0=(1.0, 1.0, 1.0), t_end: float = 40.0, Pr: float = 10.0, r: float = 28.0, b: float = 8.0 / 3.0,
                     n: int = 4001, rtol: float = 1e-10, atol: float = 1e-12, method: str = "DOP853",
                     fast: bool = False) -> dict:
    """Integrate (11.91) with an adaptive Runge–Kutta (DOP853, rtol 1e-10).

    Book: §11.14, Figs. 11.29–11.30.  Pointwise values beyond the predictability time (~20–30 time units) are not
    reproducible: assert invariants and statistics only.  fast → t_end ≤ 25.  Returns dict(t, X, Y, Z).  Label: converged.
    """
    t_end = min(t_end, 25.0) if fast else t_end
    sol = solve_ivp(lorenz_rhs, (0.0, t_end), list(s0), args=(Pr, r, b), t_eval=np.linspace(0.0, t_end, n), method=method,
                    rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(f"lorenz_integrate failed: {sol.message}")
    return dict(t=sol.t, X=sol.y[0], Y=sol.y[1], Z=sol.y[2])


def lorenz_rk4(s0, dt: float, n_steps: int, Pr: float = 10.0, r: float = 28.0, b: float = 8.0 / 3.0) -> np.ndarray:
    """Fixed-step classical RK4 for (11.91) — the twin of the explainer's JS integrator (parity rows use dt = 0.005).

    Book: §11.14.  Returns ndarray (n_steps + 1, 3).  Example: one step dt = 0.01 from (1, 1, 1) → (1.012567, 1.259918,
    0.984891).  Label: analytic (the scheme), converged (vs DOP853 for t ≲ 5).
    """
    out = np.empty((int(n_steps) + 1, 3))
    s = np.asarray(s0, dtype=float)
    out[0] = s
    f = lambda v: np.array(lorenz_rhs(0.0, v, Pr, r, b))  # noqa: E731
    for i in range(int(n_steps)):
        k1 = f(s)
        k2 = f(s + 0.5 * dt * k1)
        k3 = f(s + 0.5 * dt * k2)
        k4 = f(s + dt * k3)
        s = s + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        out[i + 1] = s
    return out


def lorenz_b(k):
    """b = 4π²/(π² + k²) (= 8/3 at k² = π²/2, the free–free critical wavenumber of (11.44)).  Book: §11.14.
    Label: analytic.
    """
    return _S(4.0 * math.pi ** 2 / (math.pi ** 2 + _F(k) ** 2))


def lorenz_r(Ra, k):
    """r = Ra/Ra_cr(k) = Ra k²/(π² + k²)³ (Ra_cr(k) from (11.44), free–free).  Book: §11.14 (r = Ra/Ra_cr).  Label: analytic."""
    k = _F(k)
    return _S(_F(Ra) * k ** 2 / (math.pi ** 2 + k ** 2) ** 3)


def lorenz_fixed_points(r: float, b: float = 8.0 / 3.0) -> list:
    """Steady states of (11.91): (0, 0, 0) and, for r > 1, (±√(b(r − 1)), ±√(b(r − 1)), r − 1).

    Book: §11.14, p. 528.  Returns list of 3-tuples (origin, C₊, C₋).  Example: r = 28 → (±8.4853, ±8.4853, 27).
    Label: analytic.
    """
    pts = [(0.0, 0.0, 0.0)]
    if r > 1:
        a = math.sqrt(b * (r - 1.0))
        pts += [(a, a, r - 1.0), (-a, -a, r - 1.0)]
    return pts


def lorenz_jacobian(s, Pr: float = 10.0, r: float = 28.0, b: float = 8.0 / 3.0) -> np.ndarray:
    """Jacobian of (11.91): [[−Pr, Pr, 0], [r − Z, −1, −X], [Y, X, −b]] (D25, ours).  Label: analytic."""
    X, Y, Z = s
    return np.array([[-Pr, Pr, 0.0], [r - Z, -1.0, -X], [Y, X, -b]])


def lorenz_eigs(r: float, Pr: float = 10.0, b: float = 8.0 / 3.0, which: str = "C") -> np.ndarray:
    """Eigenvalues of the Jacobian at a fixed point, sorted by descending real part ("C" = C₊ (same as C₋), "O" = origin).

    Book: §11.14 (D25, ours).  Example r = 28: C → 0.0940 ± 10.1945i, −13.8546; origin → 11.8277, −2.6667, −22.8277.
    Label: analytic.
    """
    pts = lorenz_fixed_points(r, b)
    if which == "C":
        if len(pts) < 3:
            raise ValueError("C± exist only for r > 1")
        p = pts[1]
    elif which in ("O", "origin"):
        p = pts[0]
    else:
        raise ValueError('which must be "C" or "O"')
    ev = np.linalg.eigvals(lorenz_jacobian(p, Pr, r, b))
    return ev[np.argsort(-ev.real, kind="stable")]


def lorenz_hopf_r(Pr: float = 10.0, b: float = 8.0 / 3.0) -> float:
    """r_H = Pr(Pr + b + 3)/(Pr − b − 1): above it C± are unstable (subcritical Hopf; D25, ours).

    Book: §11.14 says only "if r is large" (p. 528).  Pr = 10, b = 8/3 → 24.7368.  Validation: V1 eigenvalues of
    :func:`lorenz_jacobian` at C± are ±9.6245i at r_H; V5 Wikipedia "Lorenz system" (24.74).  Label: analytic, benchmark.
    """
    if Pr <= b + 1:
        raise ValueError("no Hopf bifurcation of the convection states for Pr <= b + 1")
    return Pr * (Pr + b + 3.0) / (Pr - b - 1.0)


def lorenz_divergence(Pr: float = 10.0, b: float = 8.0 / 3.0) -> float:
    """∇·ṡ = −(Pr + 1 + b): every phase-space volume contracts (dissipative).  Label: analytic."""
    return -(Pr + 1.0 + b)


def lorenz_separation(s0=(1.0, 1.0, 1.0), delta0: float = 1e-8, t_end: float = 40.0, Pr: float = 10.0, r: float = 28.0,
                      b: float = 8.0 / 3.0, n: int = 4001, window: Sequence[float] | None = None) -> dict:
    """Two Lorenz runs started δ₀ apart (in X): separation |δ|(t) and the slope of ln|δ| over the growth window.

    Book: §11.14, p. 528 (sensitivity to initial conditions).  Integrated together (DOP853, rtol 1e-12).  Default window:
    from the first time |δ| > 100δ₀ to the first time |δ| > 0.1 (before saturation).
    Returns dict(t, sep, slope (1/time; ≈ the largest Lyapunov exponent ~0.9, qualitative), window, a, b).
    Example: δ₀ = 1e-8 from (1, 1, 1): |δ| reaches 1 at t ≈ 29.1.  Label: converged.
    """
    def rhs6(t, y):
        return lorenz_rhs(t, y[:3], Pr, r, b) + lorenz_rhs(t, y[3:], Pr, r, b)

    y0 = list(s0) + [s0[0] + delta0, s0[1], s0[2]]
    sol = solve_ivp(rhs6, (0.0, t_end), y0, t_eval=np.linspace(0.0, t_end, n), method="DOP853", rtol=1e-12, atol=1e-14)
    d = np.linalg.norm(sol.y[:3] - sol.y[3:], axis=0)
    if window is None:
        i0 = np.nonzero(d > 100 * delta0)[0]
        i1 = np.nonzero(d > 0.1)[0]
        window = (float(sol.t[i0[0]]) if i0.size else 0.0, float(sol.t[i1[0]]) if i1.size else float(sol.t[-1]))
    msk = (sol.t >= window[0]) & (sol.t <= window[1]) & (d > 0)
    slope = float(np.polyfit(sol.t[msk], np.log(d[msk]), 1)[0]) if msk.sum() > 2 else float("nan")
    return dict(t=sol.t, sep=d, slope=slope, window=tuple(window), a=sol.y[:3], b=sol.y[3:])


def lorenz_predictability_time(delta0: float = 1e-8, threshold: float = 1.0, s0=(1.0, 1.0, 1.0), **kw) -> float:
    """First time two runs started δ₀ apart differ by more than ``threshold`` (NaN if never).  Example: ≈ 29.1 (δ₀ = 1e-8, 1).
    Label: converged.
    """
    r = lorenz_separation(s0=s0, delta0=delta0, **kw)
    j = np.nonzero(r["sep"] > threshold)[0]
    return float(r["t"][j[0]]) if j.size else float("nan")


def lorenz_largest_lyapunov(t_end: float = 200.0, renorm_dt: float = 1.0, cache: bool = True, delta0: float = 1e-8,
                            s0=(1.0, 1.0, 1.0), transient: float = 20.0, Pr: float = 10.0, r: float = 28.0,
                            b: float = 8.0 / 3.0) -> float:
    """Largest Lyapunov exponent by two-trajectory renormalisation (Benettin et al. 1980; ours, optional, cached).

    Book: §11.14 (not quantified in the book).  ≈ 0.9 at Lorenz's parameters (literature ≈ 0.906, quoted for orientation
    only — not a cited benchmark row).  A finite-time statistical estimate: it depends on t_end, the start and the
    renormalisation interval, and no convergence study in t_end is run.
    Validation: only the band 0.7 < λ < 1.2 (positive ⇒ sensitive dependence) is tested.  Label: qualitative.
    """
    def compute():
        s = solve_ivp(lorenz_rhs, (0.0, transient), list(s0), args=(Pr, r, b), method="DOP853", rtol=1e-11,
                      atol=1e-13).y[:, -1]
        p = s + np.array([delta0, 0.0, 0.0])
        acc, steps = 0.0, int(t_end / renorm_dt)
        for _ in range(steps):
            sa = solve_ivp(lorenz_rhs, (0.0, renorm_dt), s, args=(Pr, r, b), method="DOP853", rtol=1e-11, atol=1e-13).y[:, -1]
            sb = solve_ivp(lorenz_rhs, (0.0, renorm_dt), p, args=(Pr, r, b), method="DOP853", rtol=1e-11, atol=1e-13).y[:, -1]
            dd = np.linalg.norm(sb - sa)
            acc += math.log(dd / delta0)
            s, p = sa, sa + (sb - sa) * (delta0 / dd)
        return dict(lam=acc / (steps * renorm_dt))

    return float(_cached("lorenz_lyapunov", dict(t_end=t_end, dt=renorm_dt, d0=delta0, s0=list(s0), tr=transient, Pr=Pr, r=r,
                                                 b=b), compute, cache)["lam"])


def lorenz_fields(x, z, X: float, Y: float, Z: float, k: float = math.pi / math.sqrt(2.0)) -> dict:
    """Roll stream function, temperature anomaly and velocities of Lorenz's truncation (11.90) (proportionality constants 1).

    Book: §11.14, Eq. (11.90) ψ ∝ X cos(πz) sin(kx), T′ ∝ Y cos(πz) cos(kx) + Z sin(2πz), −½ ≤ z ≤ ½, with the §11.14
    convention u = −∂ψ/∂z, w = ∂ψ/∂x (slip S9).  Returns dict(psi, T, u, w).  Label: analytic.
    """
    x, z = np.broadcast_arrays(_F(x), _F(z))
    psi = X * np.cos(math.pi * z) * np.sin(k * x)  # Eq. (11.90)
    T = Y * np.cos(math.pi * z) * np.cos(k * x) + Z * np.sin(2 * math.pi * z)
    u = X * math.pi * np.sin(math.pi * z) * np.sin(k * x)  # −∂ψ/∂z
    w = X * k * np.cos(math.pi * z) * np.cos(k * x)  # ∂ψ/∂x
    return dict(psi=psi, T=T, u=u, w=w)


def lorenz_r_sweep(r_values=None, t_end: float = 30.0, n: int = 1500, cache: bool = True, s0=(1.0, 1.0, 1.0),
                   Pr: float = 10.0, b: float = 8.0 / 3.0) -> dict:
    """Precomputed X(t), Z(t) for a list of r (conduction → steady convection → chaos) for the slider figure F8 (cached).

    Book: §11.14.  Returns dict(r, t, X (n_r × n), Z).  Label: converged.
    """
    r_values = [0.5, 5.0, 10.0, 15.0, 20.0, 24.0, 25.0, 28.0, 30.0] if r_values is None else [float(v) for v in r_values]

    def compute():
        res = [lorenz_integrate(s0, t_end, Pr, rv, b, n=n) for rv in r_values]
        return dict(r=np.array(r_values), t=res[0]["t"], X=np.array([d["X"] for d in res]), Z=np.array([d["Z"] for d in res]))

    return _cached("lorenz_r_sweep", dict(r=r_values, t_end=t_end, n=n, s0=list(s0), Pr=Pr, b=b), compute, cache)


def lorenz_sympy() -> dict:
    """D24: Lorenz's system (11.91) by Galerkin projection of the 2-D Boussinesq equations onto the three modes of (11.90).

    Book: §11.14, (11.90)–(11.91) ("on substitution … Lorenz finally obtained" — the steps are ours).  Non-dimensional
    Boussinesq (lengths d, time d²/κ, T′ by ΔT, stress-free walls z = ±½), u = −ψ_z, w = ψ_x (§11.14 sign):
    (1/Pr)(∂_t∇²ψ + u·∇∇²ψ) = Ra ∂T′/∂x + ∇⁴ψ,  ∂_tT′ + u·∇T′ = ψ_x + ∇²T′.
    Truncation ψ = A cos πz sin kx, T′ = B cos πz cos kx + C sin 2πz; projection (orthogonality on [−½, ½] × one wavelength):
    Ȧ = Pr(kRaB/a² − a²A), Ḃ = kA − a²B − πkAC, Ċ = (πk/2)AB − 4π²C, a² = π² + k².  Scaling τ = a²t, A = √2a²X/(πk),
    B = √2Y/(πr), C = Z/(πr), r = k²Ra/a⁶ gives (11.91) with b = 4π²/a².
    Returns dict(dA, dB, dC, scaled_residuals ([0, 0, 0]), b_at_kc (8/3), r_formula).  Validation: V2.  Label: symbolic.
    """
    x, z, t = sp.symbols("x z t", real=True)
    k, Pr, Ra = sp.symbols("k Pr Ra", positive=True)
    A, B, C = [sp.Function(nm)(t) for nm in ("A", "B", "C")]
    psi = A * sp.cos(sp.pi * z) * sp.sin(k * x)  # Eq. (11.90)
    T = B * sp.cos(sp.pi * z) * sp.cos(k * x) + C * sp.sin(2 * sp.pi * z)
    lap = lambda f: f.diff(x, 2) + f.diff(z, 2)  # noqa: E731
    u, w = -psi.diff(z), psi.diff(x)
    adv = lambda f: u * f.diff(x) + w * f.diff(z)  # noqa: E731
    vort = (lap(psi).diff(t) + adv(lap(psi))) / Pr - Ra * T.diff(x) - lap(lap(psi))
    temp = T.diff(t) + adv(T) - psi.diff(x) - lap(T)
    Lx = 2 * sp.pi / k
    proj = lambda f, mode: sp.integrate(sp.integrate(sp.expand(f * mode), (x, 0, Lx)), (z, -sp.Rational(1, 2), sp.Rational(1, 2)))  # noqa: E731
    mA, mB, mC = sp.cos(sp.pi * z) * sp.sin(k * x), sp.cos(sp.pi * z) * sp.cos(k * x), sp.sin(2 * sp.pi * z)
    dA = sp.solve(proj(vort, mA), A.diff(t))[0]
    dB = sp.solve(proj(temp, mB), B.diff(t))[0]
    dC = sp.solve(proj(temp, mC), C.diff(t))[0]
    a2 = sp.pi ** 2 + k ** 2
    X, Y, Z, r = sp.symbols("X Y Z r", positive=True)
    al, be, ga = sp.sqrt(2) * a2 / (sp.pi * k), sp.sqrt(2) / (sp.pi * r), 1 / (sp.pi * r)
    sub = {A: al * X, B: be * Y, C: ga * Z, Ra: r * a2 ** 3 / k ** 2}
    bb = 4 * sp.pi ** 2 / a2
    resX = sp.simplify(dA.subs(sub) / (a2 * al) - Pr * (Y - X))  # dX/dτ = (dA/dt)/(a²α)
    resY = sp.simplify(dB.subs(sub) / (a2 * be) - (-X * Z + r * X - Y))
    resZ = sp.simplify(dC.subs(sub) / (a2 * ga) - (X * Y - bb * Z))
    return dict(dA=sp.simplify(dA), dB=sp.simplify(dB), dC=sp.simplify(dC), scaled_residuals=[resX, resY, resZ],
                b_at_kc=sp.simplify(bb.subs(k, sp.pi / sp.sqrt(2))), r_formula=sp.Eq(r, k ** 2 * Ra / a2 ** 3))


def logistic_map(A: float, x0: float, n: int) -> np.ndarray:
    """Iterates x_{n+1} = Ax_n(1 − x_n), returning x_0 … x_n.

    Book: Ex. 11.14 (logistic map as a transition toy; background state x* = 1 − 1/A, stable for 1 < A < 3 — D34).
    Label: analytic.
    """
    x = np.empty(int(n) + 1)
    x[0] = x0
    for i in range(int(n)):
        x[i + 1] = A * x[i] * (1.0 - x[i])
    return x


def logistic_fixed_point(A: float) -> dict:
    """x* = 1 − 1/A and its multiplier f′(x*) = 2 − A (stable iff 1 < A < 3).  Book: Ex. 11.14(a).  Label: analytic."""
    lam = 2.0 - A
    return dict(x_star=1.0 - 1.0 / A, multiplier=lam, stable=bool(abs(lam) < 1.0))


def cobweb(A: float, x0: float, n: int) -> dict:
    """Cobweb polyline of the logistic map: (x₀, 0) → (x₀, x₁) → (x₁, x₁) → (x₁, x₂) → …  Book: §11.14 (B1).
    Returns dict(xs, ys).  Label: analytic.
    """
    seq = logistic_map(A, x0, n)
    xs, ys = [seq[0]], [0.0]
    for i in range(int(n)):
        xs += [seq[i], seq[i + 1]]
        ys += [seq[i + 1], seq[i + 1]]
    return dict(xs=np.array(xs), ys=np.array(ys))


def bifurcation_diagram(A_values, n_transient: int = 500, n_keep: int = 100, x0: float = 0.5) -> dict:
    """Points (A, x) of the logistic map's attractor after a transient — the bifurcation tree of Fig. 11.31.

    Book: §11.14, Fig. 11.31.  Returns dict(A, x) (flat arrays).  Label: converged.
    """
    A = _F(A_values)
    x = np.full_like(A, x0)
    for _ in range(int(n_transient)):
        x = A * x * (1.0 - x)
    xs = np.empty((int(n_keep), A.size))
    for i in range(int(n_keep)):
        x = A * x * (1.0 - x)
        xs[i] = x
    return dict(A=np.tile(A, int(n_keep)), x=xs.ravel())


def _f_iter(A, x, p):
    for _ in range(p):
        x = A * x * (1.0 - x)
    return x


def superstable_points(n_max: int = 8) -> np.ndarray:
    """Superstable parameters S_n (x = ½ on the stable 2ⁿ-cycle): S₀ = 2, S₁ = 1 + √5, then f_A^{2ⁿ}(½) = ½ by Brent near the
    geometric prediction.  Expect 2, 3.236068, 3.498562, 3.554641, 3.566667, 3.569244, 3.569795.  Label: converged.
    """
    s = [2.0, 1.0 + math.sqrt(5.0)]
    for n in range(2, int(n_max) + 1):
        gap = (s[-1] - s[-2]) / 4.669
        pred = s[-1] + gap
        g = lambda A, n=n: _f_iter(A, 0.5, 2 ** n) - 0.5  # noqa: E731
        grid = np.linspace(pred - 0.35 * gap, pred + 0.35 * gap, 41)
        v = [g(A) for A in grid]
        roots = [brentq(g, grid[j], grid[j + 1], xtol=1e-15) for j in range(40) if v[j] * v[j + 1] < 0]
        if not roots:
            raise RuntimeError(f"superstable_points: no root near {pred}")
        s.append(min(roots, key=lambda q: abs(q - pred)))
    return np.array(s[: int(n_max) + 1])


def period_doubling_points(n_max: int = 6) -> dict:
    """Period-doubling parameters A_n (period 2^{n−1} → 2ⁿ) and superstable points S_n of the logistic map.

    Book: §11.14, Fig. 11.31 (R_n), Feigenbaum ratio (p. 530).  Method (ours): A_n is where the multiplier of the 2^{n−1}-cycle
    through (near) ½ equals −1, bracketed by S_{n−1} < A_n < S_n.  Returns dict(A_n (3, 1 + √6 = 3.449490, 3.544090,
    3.564407, …), S_n (S₀ … S_{n_max})).  Validation: V1 A₁, A₂; V5 Wikipedia "Feigenbaum constants".  Label: converged, benchmark.
    """
    s = superstable_points(n_max + 1)
    out = [3.0]
    for n in range(2, int(n_max) + 1):
        p = 2 ** (n - 1)

        def mult(A, p=p):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)  # fsolve's "xtol too small" at the converged root
                xx = fsolve(lambda q: _f_iter(A, q, p) - q, 0.5, xtol=1e-14)[0]
            lam = 1.0
            for _ in range(p):
                lam *= A * (1.0 - 2.0 * xx)
                xx = A * xx * (1.0 - xx)
            return lam + 1.0

        out.append(brentq(mult, s[n - 1] + 1e-12, s[n] - 1e-12, xtol=1e-14))
    return dict(A_n=np.array(out), S_n=s[: int(n_max) + 1])


def feigenbaum_estimate(n_max: int = 6) -> list:
    """Feigenbaum δ estimates δ_n = (S_n − S_{n−1})/(S_{n+1} − S_n) from superstable points S₀ … S_{n_max}.

    Book: §11.14, p. 530: (R_n − R_{n−1})/(R_{n+1} − R_n) → 4.6692.  Expect 4.709, 4.681, 4.663, 4.668, 4.669 (n_max = 6);
    within 1e-3 of δ = 4.669201609… by n_max = 8 (Feigenbaum 1978; Wikipedia "Feigenbaum constants").  Returns list.
    Label: converged, benchmark.
    """
    s = superstable_points(n_max)
    return [float(v) for v in (s[1:-1] - s[:-2]) / (s[2:] - s[1:-1])]


# ======================================================================================================================
# printed slips, sympy roll-call
# ======================================================================================================================
def book_slips() -> list:
    """Printed slips of Chapter 11 (S1–S12; "slip #1 … #12" in reader-facing text): key, where, printed, correct, evaluator.

    Book: analysis §9 (each read on the page image); design convention 5.  Returns list of dicts.  Label: analytic (each coded
    slip has a planted-variant test).
    """
    rows = [
        ("S1", "§11.4, p. 489, derivative line of the even solution",
         "(d²/dz² − K²)²W = A(q₀² + K²)²cos q₀z + B(q² − K²)²cosh qz + B(q*² − K²)²cosh q*z",
         "… + C(q*² − K²)²cosh q*z (the 3 × 3 matrix below it is right)", "benard_determinant(Ra, K, printed=True)"),
        ("S2", "§11.4, p. 490", "W = A sin(nπz) on −½ ≤ z ≤ ½", "W = A sin nπ(z + ½) (n = 1: cos πz); (11.44) unaffected",
         "benard_free_free_mode(z, 1, printed=True)"),
        ("S3", "§11.4, p. 490, dRa/dK²", "3(π² + K²)²/K² − 3(π² + K²)³/K⁴", "3(π² + K²)²/K² − (π² + K²)³/K⁴ ⇒ K² = π²/2",
         "benard_free_free_sympy()['printed_root'] == []"),
        ("S4", "§11.7, (11.55) w-equation", "−(1/ρ₀)∂p/∂x", "−(1/ρ₀)∂p/∂z",
         "stratified_shear_sympy()['printed_slip4_reaches_tg'] is False"),
        ("S5", "§11.3, p. 480", "see (7.96)", "Ch. 7 (7.95), the interface dispersion relation", "text box"),
        ("S6", "§11.11, p. 520, Tollmien profile", "middle branch 1 − b[1 − (y/δ)²]",
         "1 − b(1 − y/δ)² (continuous value and slope at the breakpoint)", "tollmien_profile(eta, printed=True)"),
        ("S7", "Exercise 11.9, (11.93)", "(d²/dR² − k²)²û_φ = −Tak²û_R", "(d²/dR² − k²)û_φ = −Tak²û_R ((11.51) at σ = 0)",
         "taylor_galerkin_Ta(k, mu, printed=True)"),
        ("S8", "§11.7, p. 507, Howard's chain", "unweighted ∫[U_min − U][U_max − U]dz ≤ 0 'recast' with Q; (U_min − U) < 0",
         "the weight Q ≥ 0 inside from the start; (U_min − U) ≤ 0", "text box"),
        ("S9", "§11.7, §11.8, §11.14", "three stream-function sign conventions",
         "not an error: each section states its convention", "docstrings"),
        ("S10", "§11.4 (11.21); §11.5", "Γ = −dT̄/dz (opposite to Ch. 1's Γ ≡ dT/dz); §11.5 redefines Ra with +dT̄/dz",
         "convention change: rayleigh_number takes dT = T_bottom − T_top", "gamma_conventions, thermal_rayleigh_signed"),
        ("S11", "§11.7, p. 506, divergence form", "d/dz[(U − c)²F′] − k²(U − c)F + N²F = 0", "… − k²(U − c)²F + N²F = 0",
         "stratified_shear_sympy()['printed_slip11_residual'] != 0"),
        ("S12", "§11.6, (11.47) and (11.50) continuity", "∂(Rũ_R)/∂R + ∂ũ_z/∂z = 0", "(1/R)∂(Rũ_R)/∂R + ∂ũ_z/∂z = 0",
         "taylor_perturbation_sympy()['printed_continuity_units_ok'] is False"),
    ]
    return [dict(key=k, where=w, printed=p, correct=c, evaluator=e) for k, w, p, c, e in rows]


def _zero(v) -> bool:
    if isinstance(v, (list, tuple)):
        return all(_zero(e) for e in v)
    try:
        return sp.simplify(v) == 0
    except Exception:  # noqa: BLE001
        return False


def derive_all() -> dict:
    """Run every sympy engine of the chapter and summarise its residuals (True = all expected zeros are zero).

    Book: the derivations D02–D24 of the design (Part F).  Returns dict(engine → dict(ok, detail)).  Label: symbolic.
    """
    out = {}
    r = kh_sympy()
    out["kh_sympy"] = dict(ok=_zero([r["kinematic_lin"], r["dynamic_lin"], r["residual_11_18"], r["identity"], r["criterion"]]),
                           detail="(11.9), (11.13), (11.18), identity, criterion")
    r = kh_depth_tension_sympy()
    out["kh_depth_tension_sympy"] = dict(ok=_zero([r["residual"], r["limit_residual"]]), detail="Ex. 11.1 form, h → ∞")
    r = benard_perturbation_sympy()
    out["benard_perturbation_sympy"] = dict(
        ok=_zero([r["eq_11_27"], r["eq_11_29_residual"], r["eq_11_32"], r["eq_11_36"], r["eq_11_37"], r["eq_11_40_residual"]])
        and not _zero(r["planted_residual"]), detail="(11.29), (11.31)–(11.37), (11.40); planted ∇² fails")
    r = exchange_of_stabilities_sympy()
    out["exchange_of_stabilities_sympy"] = dict(ok=_zero([r["relation_T"], r["relation_W"], r["imag_identity"]]),
                                                detail="Ex. 11.6")
    r = benard_free_free_sympy()
    out["benard_free_free_sympy"] = dict(ok=_zero([r["residual"], r["bc_W4"], r["dRa_minus_book"]]) and r["printed_root"] == [],
                                         detail="(11.44), K² = π²/2, slip S3")
    r = taylor_perturbation_sympy()
    out["taylor_perturbation_sympy"] = dict(
        ok=_zero([r["lin_R"], r["lin_phi"], r["lin_z"], r["base_R"], r["base_phi"], r["term_2A"], r["operator_identity"],
                  r["Ta_identity"], r["continuity_correct"]]) and not r["printed_continuity_units_ok"],
        detail="(11.50), 2A, Ta = −4AΩ₁d⁴/ν², slip S12")
    r = stratified_shear_sympy()
    out["stratified_shear_sympy"] = dict(
        ok=_zero([r["tg_residual"], r["self_adjoint_residual"], r["divergence_residual"], r["normal_modes"]])
        and not _zero(r["printed_slip11_residual"]) and not r["printed_slip4_reaches_tg"], detail="(11.58)–(11.64), (11.68)")
    r = os_derivation_sympy()
    out["os_derivation_sympy"] = dict(ok=_zero(r["residual"]) and not _zero(r["planted_residual"]), detail="(11.79)")
    r = energy_equation_sympy()
    out["energy_equation_sympy"] = dict(ok=_zero(r["residual"]), detail="(11.88) pointwise identity")
    r = lorenz_sympy()
    out["lorenz_sympy"] = dict(ok=_zero(r["scaled_residuals"]) and sp.simplify(r["b_at_kc"] - sp.Rational(8, 3)) == 0,
                               detail="(11.91) from (11.90)")
    return out


# ======================================================================================================================
# our published tables (reference/ch11) — written by reference/ch11/make_refs.py and scripts/ch11_tables.py
# ======================================================================================================================
def _sig(v, n: int = 4):
    if isinstance(v, (list, tuple, np.ndarray)):
        return [_sig(e, n) for e in v]
    v = float(v)
    return None if not np.isfinite(v) else float(f"{v:.{n}g}")


_REF_HDR = "# ours (fluidpy.ch11_instability), computed; not book data\n"
_TG_CSV_NOTE = ("# U = tanh z, N^2 = J sech^2 z; k c_i of the leading converged unstable Taylor-Goldstein mode (N = 100, box "
                "max(30, 12/k), tan-map scale min(0.5, max(0.035, 0.025/k))); every entry with J < k(1 - k) (the exact neutral "
                "curve) is positive, so the zero strip under the curve is 0.0075 (k = 0.05, 0.15, ..., 0.95) or 0.01 "
                "(k = 0.1, 0.2, ..., 0.9, where the curve passes through a grid point) wide - one J step at most; off the "
                "grid a weak mode (k c_i < 0.004) within 0.0002 (k = 0.05) to 0.006 (k = 0.95) of the curve is still "
                "reported as 0")


def _write_tg_growth_csv(path: Path, gm: dict) -> Path:
    """Write a :func:`tg_growth_map` result as the published CSV (4 s.f.; header states what a 0 means)."""
    with Path(path).open("w", encoding="utf-8") as f:
        f.write(_REF_HDR + _TG_CSV_NOTE + "\nJ," + ",".join(f"k={k:.4g}" for k in gm["k"]) + "\n")
        for J, row in zip(gm["J"], gm["kci"]):
            f.write(f"{J:.4g}," + ",".join(f"{v:.4g}" for v in row) + "\n")
    return Path(path)


def _tg_map_json(gm: dict) -> dict:
    """The ``tg_map`` entry of explainer_tables.json (4 s.f.) for a :func:`tg_growth_map` result."""
    return dict(k=_sig(gm["k"]), J=_sig(gm["J"]), kci=[_sig(r) for r in gm["kci"]], neutral_J=_sig(gm["neutral_J"]))


def write_reference_tables(out_dir: str | Path | None = None, fast: bool = False, cache: bool = True,
                           include_os: bool = True) -> dict:
    """Write every table the explainers and the notebook read to ``reference/ch11/`` (ours; ≤ 4 s.f. in the JSON).

    Book: §11.4–§11.14.  Files: benard_neutral_curves.csv (:func:`benard_neutral_table`), taylor_critical.csv,
    tg_growth_map.csv, rayleigh_spectra.json, critical_points.json, explainer_tables.json (Part C.5 5.4–5.7, 5.9: Bénard table,
    Taylor table + eigenfunction ψ on 41 × 41 for μ = 0, 0.5, 1, TG map, Rayleigh spectra, Lorenz r-sweep (X, Z on 300
    points)), and with ``include_os`` the Orr–Sommerfeld tables of :func:`neutral_curve_tables`.
    Runtime: ~10 min uncached (TG map dominates); every piece is cached in ``outputs/ch11/cache``.
    Returns dict(name → path).  Label: converged.
    """
    out = Path(out_dir) if out_dir is not None else _root() / "reference" / "ch11"
    out.mkdir(parents=True, exist_ok=True)
    paths: dict = {}
    hdr = "# ours (fluidpy.ch11_instability), computed; not book data\n"
    bt = benard_neutral_table(np.linspace(0.5, 10.0, 24 if fast else 96), write=True, cache=cache)
    paths["benard_neutral_curves"] = str(_root() / "reference" / "ch11" / "benard_neutral_curves.csv")
    mus = np.round(np.arange(-0.5, 1.0001, 0.25 if fast else 0.05), 10)
    tt = taylor_critical_table(mus, cache=cache)
    p = out / "taylor_critical.csv"
    with p.open("w", encoding="utf-8") as f:
        f.write(hdr + "# narrow gap (11.51)-(11.53), Chebyshev N = 40; Ta_11_54 = 1708/(0.5(1 + mu))\n"
                      "mu,Ta_c,k_c,Ta_11_54,rel_error\n")
        for row in zip(tt["mu"], tt["Ta_c"], tt["k_c"], tt["approx"], tt["rel_error"]):
            f.write(",".join(f"{v:.6g}" for v in row) + "\n")
    paths["taylor_critical"] = str(p)
    gm = tg_growth_map(fast=fast, cache=cache)
    p = _write_tg_growth_csv(out / "tg_growth_map.csv", gm)
    paths["tg_growth_map"] = str(p)
    rs = rayleigh_spectrum_table(ks=np.round(np.arange(0.1, 2.0001, 0.3 if fast else 0.1), 10), cache=cache, write=True)
    paths["rayleigh_spectra"] = str(_root() / "reference" / "ch11" / "rayleigh_spectra.json")
    tef = {}
    for mu in (0.0, 0.5, 1.0):
        kc = taylor_critical(mu)["k_c"]
        e = taylor_eigenfunction(kc, mu)
        tef[f"{mu:g}"] = dict(k_c=_sig(kc, 6), Ta=_sig(e["Ta"], 6), x=_sig(e["x"][0]), z=_sig(e["z"][:, 0]),
                              psi=[_sig(r) for r in e["psi"]])
    lz = lorenz_r_sweep(cache=cache)
    sel = np.linspace(0, len(lz["t"]) - 1, 300).astype(int)
    tables = dict(
        note="ours (fluidpy.ch11_instability.write_reference_tables), computed; values to 4 significant figures; not book data",
        benard=dict(K=_sig(bt["K"]), rigid=_sig(bt["rigid"]), free=_sig(bt["free"]), rigid_free=_sig(bt["rigid_free"]),
                    odd=_sig(bt["odd"])),
        taylor=dict(mu=_sig(tt["mu"]), Ta_c=_sig(tt["Ta_c"]), k_c=_sig(tt["k_c"]), approx=_sig(tt["approx"]),
                    eigenfunctions=tef),
        tg_map=_tg_map_json(gm),
        rayleigh={nm: dict(k=_sig(d["k"]), c_r=_sig(d["c_r"]), c_i=_sig(d["c_i"])) for nm, d in rs.items()},
        lorenz_sweep=dict(r=_sig(lz["r"]), t=_sig(lz["t"][sel]), X=[_sig(r[sel]) for r in lz["X"]],
                          Z=[_sig(r[sel]) for r in lz["Z"]]),
    )
    p = out / "explainer_tables.json"
    p.write_text(json.dumps(tables, separators=(",", ":")), encoding="utf-8")
    paths["explainer_tables"] = str(p)
    crit = dict(note="ours (fluidpy.ch11_instability), computed; published benchmarks are cited in SOURCES.md / benchmarks.json",
                benard={nm: benard_critical(bc, md) for nm, bc, md in
                        (("rigid_rigid", ("rigid", "rigid"), "even"), ("free_free", ("free", "free"), "even"),
                         ("rigid_free", ("rigid", "free"), "any"), ("rigid_rigid_odd", ("rigid", "rigid"), "odd"))},
                taylor_mu0=taylor_critical(0.0), poiseuille=poiseuille_critical(cache=cache),
                blasius=blasius_critical(cache=cache), bickley_sinuous=bickley_critical(cache=cache),
                tanh_inviscid={k: v for k, v in tanh_max_growth().items() if k != "c"},
                lorenz_r_H=lorenz_hopf_r(), piecewise_layer_neutral_kh=piecewise_neutral_kh(),
                feigenbaum_ratios=feigenbaum_estimate(8))
    p = out / "critical_points.json"
    p.write_text(json.dumps(crit, indent=1), encoding="utf-8")
    paths["critical_points"] = str(p)
    if include_os:
        paths.update(neutral_curve_tables(out, fast=fast, cache=cache))
    return paths
