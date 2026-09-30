"""The explicit MacCormack predictor–corrector scheme and its use on the weakly compressible (isothermal, low-Mach)
Navier–Stokes equations: the lid-driven cavity and the square block in a channel of §10.5.

Book: Kundu, Cohen & Dowling 5e, Ch. 10 §10.4 Eqs. (10.95)–(10.110) (artificial compressibility, compressible NS with
μ_v = 0, p = c²ρ, MacCormack (10.100)–(10.102), its NS form (10.103)–(10.109), arrangements FF/BB …, the semi-empirical
time-step limit (10.110)); §10.5 Eqs. (10.138)–(10.155) and the six-substep cavity algorithm with coefficients a₁–a₁₁
(printed pp. 451–452).  Pages chapters/pages/ch10/p466–p469, p476–p482.

Grid: **node-based** (collocated) — ρ, ρu, ρv all live at the nodes x_i = iΔx, y_j = jΔy, array layout f[j, i]; walls are
node lines (the cavity walls are i = 0, n and j = 0, n).  Non-dimensional variables as in the book's cavity: lengths L (the
cavity side or the block side), velocity U, time L/U, density ρ₀, pressure ρ₀U²; the equation of state (10.99) becomes
p = ρ/Ma², Ma = U/c; Re = ρ₀UL/μ.  The dimensional form (10.103)–(10.109) is the same code with c₁ … c₅
(``ns_coefficients``); :func:`cavity_coefficients` returns a₁ … a₁₁ = the same numbers with c = 1/Ma and μ = 1/Re.
A state is a tuple (rho, rhou, rhov) or a dict with those keys; with ``perturbation`` (the book's practice, p. 455) ``rho``
holds ρ′ = ρ − 1.
Printed slips: R5 (cavity Step 5, a stray "+" after −a₁: ``weakly_compressible_step(printed_step5=True)`` reproduces it and
breaks the conservation form of the continuity update); R12 ((10.155) also needs Ma ≪ 1).  The block-surface densities
(10.147)–(10.154), the outflow treatment and the corner averaging are heuristic closures ("it was found that the
conditions derived from the momentum equations give better results").
Conservation: on a periodic box the interior update conserves mass to round-off; the cavity is **not** mass-conserving —
the wall-density updates (10.139)–(10.146) are not in conservation form (relative drift ≈ −0.6 %, see
:func:`cavity_maccormack`).
Low Mach: on a fixed grid the weakly compressible error *grows* as Ma falls (Δt ∝ Ma·Δx, the acoustic truncation error
∝ cΔx² accumulates); the O(Ma²) compressibility error is real but hidden behind it until the grid is fine.
Reused by Ch. 15 (compressible flow: MacCormack is the classic shock-capturing predictor–corrector).
"""
from __future__ import annotations

import hashlib

import numpy as np

from ._util import as_scalar_if_0d

__all__ = [
    "pressure_isothermal", "ns_fluxes", "maccormack_step", "maccormack_advection_1d", "lax_wendroff_advection_1d",
    "ns_coefficients", "cavity_coefficients", "coefficients_from_ns", "ns_predictor", "ns_corrector",
    "maccormack_dt", "maccormack_dt_additive", "maccormack_dt_asymptotic", "cavity_density_bc", "make_cavity_bc", "weakly_compressible_step",
    "cavity_init", "cavity_maccormack", "block_wall_density", "block_geometry", "block_init", "block_channel",
    "body_forces", "wc_taylor_green", "ARRANGEMENTS",
]

ARRANGEMENTS = ("FF/BB", "BB/FF", "FB/BF", "BF/FB")


def _tuple(state):
    if isinstance(state, dict):
        return state["rho"], state["rhou"], state["rhov"]
    return tuple(state)


def _like(template, tup):
    if isinstance(template, dict):
        out = dict(template)
        out.update(rho=tup[0], rhou=tup[1], rhov=tup[2])
        return out
    return tuple(tup)


def pressure_isothermal(rho, c: float):
    """Isothermal equation of state p = c²ρ.

    Book: §10.4, Eq. (10.99) (non-dimensional form p = ρ/Ma², §10.5).  rho [kg/m³] (or ρ/ρ₀), c [m/s] (or 1/Ma).
    Returns p [Pa] (or p/(ρ₀U²)).  Label: analytic.
    """
    return as_scalar_if_0d(c ** 2 * np.asarray(rho, dtype=float))  # Eq. (10.99)


def ns_fluxes(rho, rhou, rhov, c: float):
    """Inviscid flux vectors of the 2-D isothermal compressible equations in conservation form U_t + E_x + F_y = 0.

    Book: §10.4, Eqs. (10.96)–(10.100) with p = c²ρ (10.99): U = (ρ, ρu, ρv), E = (ρu, ρu² + c²ρ, ρuv),
    F = (ρv, ρuv, ρv² + c²ρ).  (The viscous terms of (10.97)–(10.98) are added separately, centred.)
    Returns (E, F), each of shape (3, …).  SI or non-dimensional (c = 1/Ma).  Label: analytic.
    """
    rho = np.asarray(rho, dtype=float)
    u, v = rhou / rho, rhov / rho
    E = np.stack([rhou, rhou * u + c ** 2 * rho, rhou * v])  # Eq. (10.100) with (10.96)–(10.99)
    F = np.stack([rhov, rhov * u, rhov * v + c ** 2 * rho])
    return E, F


def _split_arrangement(arrangement: str):
    pred, corr = arrangement.replace(" ", "").split("/") if "/" in arrangement else (arrangement[0], arrangement[1])
    px, cx = pred[0], corr[0]
    py, cy = (pred[1], corr[1]) if len(pred) > 1 else ("F", "B")
    return px, py, cx, cy


def maccormack_step(U, flux_E, flux_F=None, dt: float = 1e-3, dx: float = 1.0, dy: float = 1.0,
                    arrangement: str = "FF/BB"):
    """One MacCormack predictor–corrector step for U_t + E(U)_x + F(U)_y = 0 on a periodic grid.

    Book: §10.4, Eqs. (10.100)–(10.102): predictor U* = Uⁿ − (Δt/Δx)(Eⁿ_{i+1} − Eⁿ_i) − (Δt/Δy)(Fⁿ_{j+1} − Fⁿ_j);
    corrector U^{n+1} = ½[Uⁿ + U* − (Δt/Δx)(E*_i − E*_{i−1}) − (Δt/Δy)(F*_j − F*_{j−1})] (FF/BB); the arrangements BB/FF,
    FB/BF, BF/FB swap the one-sided directions (letters = x then y; predictor/corrector).  1-D: ``flux_F=None`` and
    arrangement "F/B" or "FB".
    Parameters: U (nvar, ny, nx) or (nvar, nx); flux_E, flux_F callables U → flux of the same shape; dt, dx, dy.
    Returns U^{n+1}.  Second order in space and time.
    Validation: V1 for linear advection one step equals Lax–Wendroff (identity); V3 order 2.  Label: analytic, converged.
    """
    U = np.asarray(U, dtype=float)
    px, py, cx, cy = _split_arrangement(arrangement)
    ax_x, ax_y = U.ndim - 1, U.ndim - 2

    def diff(f, kind, axis):
        return (np.roll(f, -1, axis) - f) if kind == "F" else (f - np.roll(f, 1, axis))

    Us = U - dt / dx * diff(flux_E(U), px, ax_x)  # Eq. (10.101)
    if flux_F is not None:
        Us = Us - dt / dy * diff(flux_F(U), py, ax_y)
    Un = U + Us - dt / dx * diff(flux_E(Us), cx, ax_x)
    if flux_F is not None:
        Un = Un - dt / dy * diff(flux_F(Us), cy, ax_y)
    return 0.5 * Un  # Eq. (10.102)


def maccormack_advection_1d(T0, C: float, nsteps: int = 1, arrangement: str = "FB") -> np.ndarray:
    """MacCormack for linear advection T_t + uT_x = 0 on a periodic grid, Courant number C = uΔt/Δx.

    Book: §10.4, (10.101)–(10.102) with E = uT.  "FB": forward predictor, backward corrector; "BF" the reverse.
    For linear advection each step equals Lax–Wendroff (our D18), stable iff |C| ≤ 1.  Example: [0, 0, 1, 0, 0], C = 0.5, one
    step → [0, −0.125, 0.75, 0.375, 0].  Returns T after nsteps.  Label: analytic, converged.
    """
    T = np.array(T0, dtype=float)
    fb = arrangement.replace("/", "").upper().startswith("F")
    for _ in range(int(nsteps)):
        if fb:
            Ts = T - C * (np.roll(T, -1) - T)  # Eq. (10.101), forward
            T = 0.5 * (T + Ts - C * (Ts - np.roll(Ts, 1)))  # Eq. (10.102), backward
        else:
            Ts = T - C * (T - np.roll(T, 1))
            T = 0.5 * (T + Ts - C * (np.roll(Ts, -1) - Ts))
    return T


def lax_wendroff_advection_1d(T0, C: float, nsteps: int = 1) -> np.ndarray:
    """Lax–Wendroff T^{n+1} = T − (C/2)(T₊ − T₋) + (C²/2)(T₊ − 2T + T₋) — the one-line form MacCormack reduces to (D18).
    Label: analytic."""
    T = np.array(T0, dtype=float)
    for _ in range(int(nsteps)):
        Tp, Tm = np.roll(T, -1), np.roll(T, 1)
        T = T - 0.5 * C * (Tp - Tm) + 0.5 * C ** 2 * (Tp - 2 * T + Tm)
    return T


# ======================================================================================================================
# Coefficients
# ======================================================================================================================
def ns_coefficients(dt: float, dx: float, dy: float, mu: float, c: float | None = None) -> dict:
    """c₁ = Δt/Δx, c₂ = Δt/Δy, c₃ = μΔt/Δx², c₄ = μΔt/Δy², c₅ = μΔt/(12ΔxΔy).

    Book: §10.4, Eq. (10.109) (c₅: the cross-derivative stencil carries 1/(4ΔxΔy) and the μ/3 of (10.97)).
    Units: dt [s], dx, dy [m], mu [Pa·s] (dimensional (10.103)–(10.108) with ρ in kg/m³).  ``c`` (sound speed [m/s]) is
    stored with the coefficients when given, for :func:`ns_predictor`.  Label: analytic.
    """
    out = dict(c1=dt / dx, c2=dt / dy, c3=mu * dt / dx ** 2, c4=mu * dt / dy ** 2, c5=mu * dt / (12.0 * dx * dy))  # Eq. (10.109)
    if c is not None:
        out["c"] = float(c)
    return out


def cavity_coefficients(dt: float, dx: float, dy: float, Ma: float, Re: float) -> dict:
    """a₁ … a₁₁ of the non-dimensional six-substep MacCormack algorithm.

    Book: §10.5, printed p. 452: a₁ = Δt/Δx, a₂ = Δt/Δy, a₃ = Δt/(ΔxMa²), a₄ = Δt/(ΔyMa²), a₅ = 4Δt/(3ReΔx²),
    a₆ = Δt/(ReΔy²), a₇ = Δt/(ReΔx²), a₈ = 4Δt/(3ReΔy²), a₉ = Δt/(12ReΔxΔy), a₁₀ = 2(a₅ + a₆), a₁₁ = 2(a₇ + a₈).
    Equal to :func:`coefficients_from_ns` of (10.109) with c = 1/Ma, μ = 1/Re.  The dict also carries dt, dx, dy, Ma, Re so
    it can drive :func:`weakly_compressible_step`.  Label: analytic.
    """
    a = dict(a1=dt / dx, a2=dt / dy, a3=dt / (dx * Ma ** 2), a4=dt / (dy * Ma ** 2), a5=4 * dt / (3 * Re * dx ** 2),
             a6=dt / (Re * dy ** 2), a7=dt / (Re * dx ** 2), a8=4 * dt / (3 * Re * dy ** 2), a9=dt / (12 * Re * dx * dy))
    a["a10"] = 2 * (a["a5"] + a["a6"])
    a["a11"] = 2 * (a["a7"] + a["a8"])
    a.update(dt=dt, dx=dx, dy=dy, Ma=Ma, Re=Re)
    return a


def coefficients_from_ns(c_coeffs: dict, c: float | None = None) -> dict:
    """Map c₁ … c₅ of (10.109) and the sound speed c to the a₁ … a₁₁ form (a₃ = c₁c², a₅ = 4c₃/3, a₉ = c₅ …).

    Book: §10.4 (10.103)–(10.109) ↔ §10.5 Step 2 (the pressure flux c²ρ split out of (ρu² + c²ρ)).  Label: analytic.
    """
    k = c_coeffs
    c = k.get("c") if c is None else c
    if c is None:
        raise ValueError("the sound speed c is needed (pass c= or ns_coefficients(..., c=...))")
    a = dict(a1=k["c1"], a2=k["c2"], a3=k["c1"] * c ** 2, a4=k["c2"] * c ** 2, a5=4 * k["c3"] / 3, a6=k["c4"],
             a7=k["c3"], a8=4 * k["c4"] / 3, a9=k["c5"])
    a["a10"] = 2 * (a["a5"] + a["a6"])
    a["a11"] = 2 * (a["a7"] + a["a8"])
    return a


# ======================================================================================================================
# The NS predictor and corrector (10.103)–(10.108) — interior nodes
# ======================================================================================================================
def _D(f, kind: str, axis: int):
    """One-sided difference on the interior [1:-1, 1:-1]: F → f_{+1} − f_0, B → f_0 − f_{−1}."""
    if axis == 1:
        return (f[1:-1, 2:] - f[1:-1, 1:-1]) if kind == "F" else (f[1:-1, 1:-1] - f[1:-1, :-2])
    return (f[2:, 1:-1] - f[1:-1, 1:-1]) if kind == "F" else (f[1:-1, 1:-1] - f[:-2, 1:-1])


def _stage(rho, m, n, a: dict, kx: str, ky: str, one: float = 1.0):
    """Increments (Δρ, Δm, Δn) of one MacCormack stage on the interior (a-coefficient form of (10.103)–(10.108)).
    ``rho`` is the density perturbation when one = 1 (ρ = 1 + ρ′) or the density itself when one = 0."""
    R = rho + one
    u = m / R
    v = n / R
    c = (slice(1, -1), slice(1, -1))
    dr = -a["a1"] * _D(m, kx, 1) - a["a2"] * _D(n, ky, 0)  # continuity, Eq. (10.103)/(10.106)
    cross_v = v[2:, 2:] + v[:-2, :-2] - v[:-2, 2:] - v[2:, :-2]  # v_{i+1,j+1} + v_{i−1,j−1} − v_{i+1,j−1} − v_{i−1,j+1}
    cross_u = u[2:, 2:] + u[:-2, :-2] - u[:-2, 2:] - u[2:, :-2]
    dm = (-a["a3"] * _D(rho, kx, 1) - a["a1"] * _D(m * u, kx, 1) - a["a2"] * _D(m * v, ky, 0)
          - a["a10"] * u[c] + a["a5"] * (u[1:-1, 2:] + u[1:-1, :-2]) + a["a6"] * (u[2:, 1:-1] + u[:-2, 1:-1])
          + a["a9"] * cross_v)  # Eq. (10.104)/(10.107) (Step 2 / Step 5 form)
    dn = (-a["a4"] * _D(rho, ky, 0) - a["a1"] * _D(m * v, kx, 1) - a["a2"] * _D(n * v, ky, 0)
          - a["a11"] * v[c] + a["a7"] * (v[1:-1, 2:] + v[1:-1, :-2]) + a["a8"] * (v[2:, 1:-1] + v[:-2, 1:-1])
          + a["a9"] * cross_u)  # Eq. (10.105)/(10.108)
    return dr, dm, dn


def _apply_stage(src, a, kx, ky, periodic, one, base=None, printed=False):
    """One MacCormack stage.  base None → predictor U* = U + ΔU(U); base given → corrector U = ½(Uⁿ + U* + ΔU(U*)).
    ``printed`` (corrector only): continuity as printed in cavity Step 5 (slip R5)."""
    rho, m, n = (np.asarray(f) for f in src)
    rr, mm, nn = ([np.pad(f, 1, mode="wrap") for f in (rho, m, n)] if periodic else (rho, m, n))
    dr, dm, dn = _stage(rr, mm, nn, a, kx, ky, one)
    if printed and base is not None:
        # Step 5 as printed (slip R5): 2ρ^{n+1} = (ρⁿ + ρ*) − a₁ + [(ρu)*_{i,j} − (ρu)*_{i−1,j}] − a₂[…]
        dr = -a["a1"] + _D(mm, kx, 1) - a["a2"] * _D(nn, ky, 0)
    if periodic:
        out = [rho + dr, m + dm, n + dn]
        if base is not None:
            out = [0.5 * (b + o) for b, o in zip(base, out)]  # Eqs. (10.106)–(10.108)
        return tuple(out)
    new = [rho.copy(), m.copy(), n.copy()]
    if base is None:
        new[0][1:-1, 1:-1] += dr  # Eq. (10.103)
        new[1][1:-1, 1:-1] += dm  # Eq. (10.104)
        new[2][1:-1, 1:-1] += dn  # Eq. (10.105)
        return tuple(new)
    b0, b1, b2 = base
    new[0][1:-1, 1:-1] = 0.5 * (b0[1:-1, 1:-1] + rho[1:-1, 1:-1] + dr)  # Eq. (10.106)
    new[1][1:-1, 1:-1] = 0.5 * (b1[1:-1, 1:-1] + m[1:-1, 1:-1] + dm)  # Eq. (10.107)
    new[2][1:-1, 1:-1] = 0.5 * (b2[1:-1, 1:-1] + n[1:-1, 1:-1] + dn)  # Eq. (10.108)
    return tuple(new)


def ns_predictor(state, coeffs: dict, arrangement: str = "FF/BB", periodic: bool = True, c: float | None = None):
    """Predictor (10.103)–(10.105) for (ρ, ρu, ρv) with the pressure replaced by c²ρ (dimensional form).

    Book: §10.4, Eqs. (10.103)–(10.105), coefficients (10.109).  ``state``: tuple or dict (rho, rhou, rhov) [j, i] (SI);
    ``coeffs``: :func:`ns_coefficients` (with ``c``) or already the a-form; convective/pressure differences one-sided (the
    predictor letters of ``arrangement``), viscous ones centred.  ``periodic``: all nodes updated (wrap-around); otherwise
    only interior nodes (boundaries unchanged).  Returns the starred state (same type as ``state``).  Label: analytic.
    """
    a = coeffs if "a1" in coeffs else coefficients_from_ns(coeffs, c)
    px, py, _, _ = _split_arrangement(arrangement)
    return _like(state, _apply_stage(_tuple(state), a, px, py, periodic, 0.0))


def ns_corrector(state, star, coeffs: dict, arrangement: str = "FF/BB", periodic: bool = True, c: float | None = None):
    """Corrector (10.106)–(10.108): 2Uⁿ⁺¹ = Uⁿ + U* − (opposite one-sided differences of the starred fluxes) + viscous(U*).

    Book: §10.4, Eqs. (10.106)–(10.108).  Arguments as :func:`ns_predictor`.  Returns the new state.  Label: analytic.
    """
    a = coeffs if "a1" in coeffs else coefficients_from_ns(coeffs, c)
    _, _, cx, cy = _split_arrangement(arrangement)
    return _like(state, _apply_stage(_tuple(star), a, cx, cy, periodic, 0.0, base=_tuple(state)))


# ======================================================================================================================
# Time step
# ======================================================================================================================
def maccormack_dt(u, v, c: float, dx: float, dy: float, rho: float = 1.0, mu: float = 1e-2, sigma: float = 0.8) -> float:
    """Semi-empirical stability limit of explicit MacCormack for the Navier–Stokes equations.

    Book: §10.4, Eq. (10.110) (Tannehill, Anderson & Pletcher 1997, as cited):
    Δt ≤ σ/(1 + 2/Re_Δ) · [|u|/Δx + |v|/Δy + c√(1/Δx² + 1/Δy²)]⁻¹,  Re_Δ = min(ρ|u|Δx/μ, ρ|v|Δy/μ).
    Parameters: u, v — characteristic speeds [m/s] (scalars; arrays are reduced to max|u|, max|v|); c [m/s]; dx, dy [m];
    rho [kg/m³]; mu [Pa·s]; sigma (safety factor; the formula's source recommends a value a little below 1 — our default 0.8).  Returns Δt [s].  Non-dimensional use: c = 1/Ma, μ = 1/Re, ρ = 1.
    Example (1, 0, 12.5, 1/64, 1/64, 1, 0.01) → 2.935e-4 (σ = 0.8).
    # DEVIATION: a zero speed would give Re_Δ = 0 and Δt = 0 (and with arrays the book's pointwise minimum vanishes at every
    stagnation point or wall); zero speeds are left out of the minimum, and Re_Δ = ∞ if both are zero.
    Validation: V1 → (10.155) for Δx = Δy, Re_Δ → ∞, |u| + |v| ≪ c√2.  Label: analytic (semi-empirical formula).
    """
    uu = float(np.max(np.abs(u)))
    vv = float(np.max(np.abs(v)))
    re_list = [r for r in (rho * uu * dx / mu, rho * vv * dy / mu) if r > 0]
    re_d = min(re_list) if re_list else np.inf
    br = uu / dx + vv / dy + c * np.sqrt(1.0 / dx ** 2 + 1.0 / dy ** 2)
    return float(sigma / (1.0 + 2.0 / re_d) / br)  # Eq. (10.110)


def maccormack_dt_additive(u, v, c: float, dx: float, dy: float, rho: float = 1.0, mu: float = 1e-2,
                           sigma: float = 0.8) -> float:
    """Stable MacCormack step with the convective, acoustic and viscous rates added (our default for the cavity).

    Book: §10.4, Eq. (10.110) rearranged.  (10.110) is Δt ≤ σ/[(1 + 2/Re_Δ)·B], B = |u|/Δx + |v|/Δy + c√(1/Δx² + 1/Δy²);
    expanding, (1 + 2/Re_Δ)·B = B + (2μ/ρ)/(|u|Δx)·B ⊇ B + 2ν/Δx², i.e. the formula adds a viscous rate but also
    multiplies the acoustic rate by 2/Re_Δ.  Here only the viscous rate of the centred 2-D Laplacian is added:
    Δt = σ / [B + 2ν(1/Δx² + 1/Δy²)],  ν = μ/ρ.
    The viscous term alone is the exact limit of the predictor–corrector (Heun) step for pure diffusion,
    Δt·4ν(1/Δx² + 1/Δy²) ≤ 2.
    # DEVIATION: not printed in the book — an additive combination of the limits in (10.110); it is stable in our cavity
    runs up to σ = 1 at Re = 1 … 400 on 16² … 64² (Ma = 0.08), where the inviscid bracket alone blows up at low Re.
    Parameters: u, v — characteristic speeds [m/s] (arrays reduced to max |·|); c [m/s]; dx, dy [m]; rho [kg/m³];
    mu [Pa·s]; sigma (safety factor).  Returns Δt [s].  Non-dimensional use: c = 1/Ma, μ = 1/Re, ρ = 1.
    Validation: V1 → the inviscid bracket as μ → 0; → the diffusive Heun limit as c, u, v → 0.  Label: analytic
    (semi-empirical; stability demonstrated by runs).
    """
    uu = float(np.max(np.abs(u)))
    vv = float(np.max(np.abs(v)))
    br = uu / dx + vv / dy + c * np.sqrt(1.0 / dx ** 2 + 1.0 / dy ** 2)
    return float(sigma / (br + 2.0 * (mu / rho) * (1.0 / dx ** 2 + 1.0 / dy ** 2)))  # Eq. (10.110), additive form


def maccormack_dt_asymptotic(Ma: float, dx: float, sigma: float = 0.8) -> float:
    """Asymptotic MacCormack step Δt ≤ (σ/√2) Ma Δx (non-dimensional, Δx = Δy).

    Book: §10.5, Eq. (10.155) — (10.110) with a large grid Reynolds number; slip R12: it also needs |u| + |v| ≪ c√2
    (Ma ≪ 1), which the book does not state.  Example (0.08, 1/64) → 7.071e-4 (σ = 0.8).  Label: analytic.
    """
    return float(sigma / np.sqrt(2.0) * Ma * dx)  # Eq. (10.155)


# ======================================================================================================================
# Lid-driven cavity (six-substep algorithm, (10.138)–(10.146))
# ======================================================================================================================
def _onesided(f, side: str, dt_over_2h: float):
    """dt/(2h)·[−f₂ + 4f₁ − 3f₀] along the wall normal (side = left/right/bottom/top) for all wall nodes."""
    if side == "left":
        return dt_over_2h * (-f[:, 2] + 4 * f[:, 1] - 3 * f[:, 0])
    if side == "right":
        return dt_over_2h * (-f[:, -3] + 4 * f[:, -2] - 3 * f[:, -1])
    if side == "bottom":
        return dt_over_2h * (-f[2, :] + 4 * f[1, :] - 3 * f[0, :])
    if side == "top":
        return dt_over_2h * (-f[-3, :] + 4 * f[-2, :] - 3 * f[-1, :])
    raise ValueError(side)


def cavity_density_bc(state, side: str, stage: str, dt: float, dx: float, dy: float, U: float = 1.0, star=None):
    """Wall density of the cavity from the continuity equation with one-sided second-order differences.

    Book: §10.5, Eq. (10.138) (continuity on a wall, where v = 0 — or u = 0 — along it) and Eqs. (10.139)–(10.146):
    left wall (i = 0, incl. corners) ρ* = ρⁿ − (Δt/2Δx)[−(ρu)ⁿ_{i+2} + 4(ρu)ⁿ_{i+1} − 3(ρu)ⁿ_i]; right wall (i = n_x − 1, incl.
    corners) with + and i − 1, i − 2; bottom (j = 0, interior i) with ρv; lid (j = n_y − 1, interior i):
    ρ* = ρⁿ − (ΔtU/2Δx)(ρⁿ_{i+1} − ρⁿ_{i−1}) + (Δt/2Δy)[−(ρv)ⁿ_{j−2} + 4(ρv)ⁿ_{j−1} − 3(ρv)ⁿ_j].  Corrector (stage
    "corrector", needs ``star``): 2ρⁿ⁺¹ = ρⁿ + ρ* ∓ the same bracket of the starred values.
    Parameters: state (level n; tuple or dict (rho, rhou, rhov)), side ("left" | "right" | "bottom" | "top"), stage, dt, dx, dy,
    U (lid speed), star (the starred state, corrector only).  Returns the new wall values (1-D array).  Label: analytic.
    """
    rho_n, _, _ = _tuple(state)
    src = _tuple(state) if stage == "predictor" else _tuple(star)
    rs, m, n = src
    if side == "left":
        inc = _onesided(m, "left", dt / (2 * dx))  # Eq. (10.139)/(10.140)
        base, s_val = rho_n[:, 0], rs[:, 0]
    elif side == "right":
        inc = -_onesided(m, "right", dt / (2 * dx))  # Eq. (10.141)/(10.142)
        base, s_val = rho_n[:, -1], rs[:, -1]
    elif side == "bottom":
        inc = _onesided(n, "bottom", dt / (2 * dy))[1:-1]  # Eq. (10.143)/(10.144)
        base, s_val = rho_n[0, 1:-1], rs[0, 1:-1]
    elif side == "top":
        inc = dt * U / (2 * dx) * (rs[-1, 2:] - rs[-1, :-2]) - _onesided(n, "top", dt / (2 * dy))[1:-1]  # Eq. (10.145)/(10.146)
        base, s_val = rho_n[-1, 1:-1], rs[-1, 1:-1]
    else:
        raise ValueError("side must be left, right, bottom or top")
    if stage == "predictor":
        return base - inc
    return 0.5 * (base + s_val - inc)


def make_cavity_bc(wall_density: str = "continuity", corner: str = "extrapolate", corner_u: float = 1.0):
    """Boundary function of the MacCormack cavity (Steps 3 and 6): wall densities, no-slip walls, lid u = U.

    Book: §10.5.  ``wall_density``:
    "continuity" (the book's) — Eqs. (10.139)–(10.146), a predictor–corrector for ρ on each wall from the continuity equation
    with one-sided second-order differences;
    "momentum" (ours: the book's block closure (10.147)–(10.154) transferred to the cavity walls) — the wall density from the
    wall-normal momentum equation (left wall ↔ (10.152) "back", right ↔ (10.151) "front", bottom ↔ (10.153) "top",
    lid ↔ (10.154) "bottom"; on the lid v = 0 and u = U along the wall, so the same terms vanish).
    Why the option: the continuity closure fixes the lid pressure p = ρ/Ma² through Uρ_x = −(ρv)_y, so its truncation error is
    amplified by 1/Ma² — our 32²–64² runs get *worse* as Ma decreases (the book's much finer grid hides it); the momentum
    closure ties the wall density to the interior pressure gradient instead.
    ``corner``: the two top corners (lid meets side wall) — "extrapolate" (default, linear along the side wall from the two
    nodes below) or "book" ((10.139)/(10.141) applied there as printed, "including two corner points").
    # DEVIATION: with u = U at a top corner, (10.139) becomes a downwind one-sided advection of density along the lid,
    ρ_t = −U(−ρ₂ + 4ρ₁ − 3ρ₀)/(2Δx), whose self-coefficient +3U/(2Δx) makes the corner grow exponentially (with u = 0 at the
    corner it is a sink −3ρU/(2Δx) that drives ρ negative in a few steps): both blow up in our runs, so the default
    extrapolates the corner density (the book averages at the block corners for a similar reason).
    Corners with ``wall_density="momentum"`` (ours; the wall-normal momentum closure has no normal at a corner): the two
    bottom corners take the average of their two neighbouring wall nodes, ρ(0,0) = ½[ρ(1,0) + ρ(0,1)] (as the book does at
    the block corners), and the two top corners are always extrapolated linearly down the side wall, ρ_top = 2ρ_{−2} −
    ρ_{−3} — ``corner`` is ignored in this mode.
    ``corner_u``: velocity of the top corner nodes as a fraction of U (1 = the corners move with the lid; the velocity then
    changes "linearly between two grid points" down the side wall, p. 449; it enters only the cross-derivative stencils).
    Returns bcf(state_n, stage_state, star, stage, dt, dx, dy, U, one, M=…, Re=…).  Label: analytic (closures heuristic).
    """
    def bcf(state_n, stage_state, star, stage, dt, dx, dy, U, one, M=None, Re=None):
        rho, m, n = (np.array(f) for f in stage_state)
        m[:, 0] = m[:, -1] = 0.0
        m[0, :] = 0.0
        n[:, 0] = n[:, -1] = 0.0
        n[0, :] = n[-1, :] = 0.0
        if wall_density == "continuity":
            for side in ("left", "right"):
                val = cavity_density_bc(state_n, side, stage, dt, dx, dy, U, star)
                rho[:, 0 if side == "left" else -1] = val
            rho[0, 1:-1] = cavity_density_bc(state_n, "bottom", stage, dt, dx, dy, U, star)
            rho[-1, 1:-1] = cavity_density_bc(state_n, "top", stage, dt, dx, dy, U, star)
        elif wall_density == "momentum":
            R = rho + one
            m[-1, :] = R[-1, :] * U
            uu, vv = m / R, n / R
            jj = np.arange(1, rho.shape[0] - 1)
            ii = np.arange(1, rho.shape[1] - 1)
            nx1, ny1 = rho.shape[1] - 1, rho.shape[0] - 1
            rho[jj, 0] = _wall_rho(rho, uu, vv, "back", np.zeros_like(jj), jj, dx, dy, M, Re)
            rho[jj, nx1] = _wall_rho(rho, uu, vv, "front", np.full_like(jj, nx1), jj, dx, dy, M, Re)
            rho[0, ii] = _wall_rho(rho, uu, vv, "top", ii, np.zeros_like(ii), dx, dy, M, Re)
            rho[ny1, ii] = _wall_rho(rho, uu, vv, "bottom", ii, np.full_like(ii, ny1), dx, dy, M, Re)
            rho[0, 0] = 0.5 * (rho[0, 1] + rho[1, 0])  # bottom corners: average of the two neighbouring wall nodes
            rho[0, -1] = 0.5 * (rho[0, -2] + rho[1, -1])
        else:
            raise ValueError("wall_density must be continuity or momentum")
        if corner == "extrapolate" or wall_density == "momentum":
            rho[-1, 0] = 2.0 * rho[-2, 0] - rho[-3, 0]
            rho[-1, -1] = 2.0 * rho[-2, -1] - rho[-3, -1]
        m[-1, :] = (rho[-1, :] + one) * U  # lid row: u = U (Steps 3 and 6)
        m[-1, 0] *= corner_u
        m[-1, -1] *= corner_u
        return rho, m, n

    bcf.needs_MRe = wall_density == "momentum"
    return bcf


_cavity_bcf = make_cavity_bc()


def _wc_step_bounded(state, dt, dx, dy, M, Re, bcf, arrangement, printed, one, U, mask=None):
    """Bounded-domain MacCormack step (the six substeps) with a boundary function."""
    a = cavity_coefficients(dt, dx, dy, M, Re)
    px, py, cx, cy = _split_arrangement(arrangement)
    kw = dict(M=M, Re=Re) if getattr(bcf, "needs_MRe", False) else {}
    star = _apply_stage(state, a, px, py, False, one)  # Step 2 (Step 1 inside: u = ρu/ρ)
    if mask is not None:
        star = tuple(np.where(mask, s0, s) for s0, s in zip(state, star))
    star = bcf(state, star, None, "predictor", dt, dx, dy, U, one, **kw)  # Step 3
    new = _apply_stage(star, a, cx, cy, False, one, base=state, printed=printed)  # Steps 4–5
    if mask is not None:
        new = tuple(np.where(mask, s0, s) for s0, s in zip(state, new))
    return bcf(state, new, star, "corrector", dt, dx, dy, U, one, **kw)  # Step 6


def weakly_compressible_step(state, coeffs: dict, bc_fn="cavity", arrangement: str = "FF/BB", dtype=np.float64,
                             perturbation: bool = True, printed_step5: bool = False, U: float = 1.0):
    """One step of the non-dimensional six-substep MacCormack algorithm.

    Book: §10.5, printed pp. 451–452: Step 1 velocities from momenta; Step 2 interior predictor with a₁ … a₁₁; Step 3
    boundary conditions at t_{n+1} for the starred variables; Step 4 starred velocities; Step 5 interior corrector; Step 6
    boundary conditions.  (Equivalent to (10.103)–(10.108) with c = 1/Ma, μ = 1/Re.)
    Parameters
    ----------
    state : tuple or dict (rho, rhou, rhov) node arrays [j, i]; ``rho`` is ρ′ = ρ − 1 when ``perturbation``.
    coeffs : :func:`cavity_coefficients` (carries dt, dx, dy, Ma, Re).
    bc_fn : "cavity" (walls (10.139)–(10.146), lid u = U), "periodic" (all nodes updated — tests), or a callable made by
        :func:`make_cavity_bc`.
    arrangement : "FF/BB" | "BB/FF" | "FB/BF" | "BF/FB".   dtype : float precision (the p. 455 double-precision lesson).
    printed_step5 : Step 5 continuity with the printed stray "+" (slip R5; breaks the conservation form).   U : lid speed.
    Returns the new state (same type).  Cost ~0.3 ms per step on 65² nodes.  Label: converged.
    """
    tup = tuple(np.asarray(f, dtype=dtype) for f in _tuple(state))
    one = 1.0 if perturbation else 0.0
    dt, dx, dy, Ma, Re = (coeffs[k] for k in ("dt", "dx", "dy", "Ma", "Re"))
    px, py, cx, cy = _split_arrangement(arrangement)
    if bc_fn == "periodic":
        star = _apply_stage(tup, coeffs, px, py, True, one)
        new = _apply_stage(star, coeffs, cx, cy, True, one, base=tup, printed=printed_step5)
        return _like(state, new)
    bcf = _cavity_bcf if bc_fn in (None, "cavity") else bc_fn
    return _like(state, _wc_step_bounded(tup, dt, dx, dy, Ma, Re, bcf, arrangement, printed_step5, one, U))


def cavity_init(n: int, dtype=np.float64, perturbation: bool = True, U: float = 1.0):
    """Fluid at rest (ρ = 1) in the (n + 1) × (n + 1)-node cavity with the lid row moving.  Returns (rho, rhou, rhov)."""
    rho = np.zeros((n + 1, n + 1), dtype=dtype) + (0.0 if perturbation else 1.0)
    m = np.zeros_like(rho)
    m[-1, :] = U
    return rho, m, np.zeros_like(rho)


def _node_streamfunction(u, dy):
    psi = np.zeros_like(u)
    psi[1:, :] = np.cumsum(0.5 * (u[1:, :] + u[:-1, :]), axis=0) * dy  # trapezoid in y from the bottom wall
    return psi


def cavity_maccormack(Re: float = 100.0, Ma: float = 0.08, n: int = 32, t_end: float = 20.0, arrangement: str = "FF/BB",
                      cycle: bool = False, tol_steady: float = 1e-7, cache: bool = True, dt: float | None = None,
                      sigma: float = 0.8, dt_rule: str = "additive", dtype=np.float64, perturbation: bool = True,
                      U: float = 1.0, printed_step5: bool = False, save_every: int | None = None, check_every: int = 50,
                      wall_density: str = "continuity", corner: str = "extrapolate", fast: bool = False,
                      check_stability: bool = True) -> dict:
    """Lid-driven cavity by the explicit MacCormack scheme on the weakly compressible NS equations, marched to steady state.

    Book: §10.5 (Fig. 10.6; scales L, U, L/U, ρ₀, ρ₀U²; p = ρ/Ma²; Re = ρ₀UL/μ; wall density from continuity (10.138)–(10.146);
    the six-substep algorithm).  The book's run uses a much finer grid (values private); here n cells per side (n + 1 nodes).
    Parameters
    ----------
    Re, Ma : Reynolds and Mach numbers (our default Ma = 0.08).   n : cells per side (``fast`` halves it).   t_end : max time.
    arrangement, cycle : one-sided arrangement; ``cycle`` alternates FB/BF and BF/FB (the book's block practice).
    tol_steady : stop when max|Δ(ρu)|, max|Δ(ρv)| per unit time < tol (checked every ``check_every`` steps); 0 → run to t_end.
    cache : npz in outputs/ch10 keyed by the parameters.   dt : default by ``dt_rule`` with ``sigma`` (|u| = |v| = U):
        "additive" (default, ours — see :func:`maccormack_dt_additive`) = σ/[|u|/Δx + |v|/Δy + c√(1/Δx² + 1/Δy²)
        + (2/Re)(1/Δx² + 1/Δy²)]; "tap" = the full (10.110) (the book's; very conservative when Re_Δ is small, e.g.
        Δt ×1/33 at Re = 2 on 32²); "asymptotic" = (10.155) (large Re_Δ only); "inviscid" = (10.110) without the viscous
        factor (unstable when the viscous limit binds: Re = 5 on 64², Re = 2 on 32² blow up).
    check_stability : raise ValueError when Δt (given or by the rule) exceeds the additive limit with σ = 1 — the largest
        step our runs found stable (Re = 1 … 400, 16² … 64²; the inviscid bound at σ ≈ 1 already blows up at Re = 20, 64²).
    dtype, perturbation : precision and ρ′ = ρ − 1 storage (p. 455).   printed_step5 : slip R5 variant.
    wall_density, corner : see :func:`make_cavity_bc` ("continuity" = the book's).
    Returns dict(x, y (node coordinates), u, v, rho (full density), p (= ρ/Ma², mean removed), psi (node streamfunction),
    t, steps, converged, residual, dt, mass_drift (relative), history (list of (t, max|Δ|/Δt)), snapshots, Re, Ma, n).
    Cost (additive Δt, to t = 30): 32² ≈ 10 s, 64² ≈ 50 s, 128² ≈ 8 min — 64² and up are cached / script-only.
    Validation: V5 Ghia et al. (1982) centreline — max deviation 4.1 % of U at 32², 1.38 % at 64², 0.44 % at 128²
    (Re = 100, Ma = 0.08, t = 30), i.e. inside a 1 % benchmark tolerance only at 128²; Hou et al. (1995) centre; V3 grid
    refinement (order ≈ 1.6).  Not mass-conserving: the wall-density closures (10.139)–(10.146) are not in conservation
    form, so ``mass_drift`` is −0.68 %, −0.56 %, −0.49 % at 32², 64², 128² (Ma = 0.08) — it is reported, not zero.
    Low Mach: on a fixed grid the error does not fall as Ma falls (Δt ∝ Ma·Δx and the acoustic truncation error
    accumulate; the O(Ma²) compressibility error is hidden behind it), which is why the book's cavity uses a much finer
    grid.  Label: converged; benchmark only at 128² (1.38 % of U at 64²).
    """
    if fast:
        n = max(8, n // 2)
    dx = dy = 1.0 / n
    c = 1.0 / Ma
    if dt is None:
        if dt_rule == "tap":
            dt = maccormack_dt(U, U, c, dx, dy, 1.0, 1.0 / Re, sigma)  # Eq. (10.110)
        elif dt_rule == "asymptotic":
            dt = maccormack_dt_asymptotic(Ma, dx, sigma)  # Eq. (10.155)
        elif dt_rule == "inviscid":
            dt = sigma / (2 * U / dx + c * np.sqrt(2.0) / dx)
        elif dt_rule == "additive":
            # DEVIATION: default step = additive rates instead of the full (10.110) — (10.110) multiplies the whole
            # bracket, acoustic term included, by (1 + 2/Re_Δ), which the book itself calls "quite conservative" at small
            # mesh Reynolds numbers; adding only the viscous rate keeps the convective, acoustic and viscous limits.
            dt = maccormack_dt_additive(U, U, c, dx, dy, 1.0, 1.0 / Re, sigma)
        else:
            raise ValueError("dt_rule must be additive, tap, asymptotic or inviscid")
    if check_stability:
        lim = maccormack_dt_additive(U, U, c, dx, dy, 1.0, 1.0 / Re, 1.0)
        if dt > lim * (1 + 1e-12):
            raise ValueError(f"dt = {dt:.4g} exceeds the MacCormack stability limit {lim:.4g} (additive form of (10.110) "
                             f"with sigma = 1; Re = {Re:g}, Ma = {Ma:g}, n = {n}); use dt_rule='additive' or a smaller dt, or "
                             f"pass check_stability=False to watch it blow up")
    path = None
    if cache and not save_every:
        from .project import repo_root

        tag = (f"Re{Re}_Ma{Ma}_n{n}_t{t_end}_{arrangement}_{cycle}_tol{tol_steady}_dt{dt}_{np.dtype(dtype).name}_"
               f"{perturbation}_{U}_{printed_step5}_{wall_density}_{corner}")
        h = hashlib.md5(tag.encode()).hexdigest()[:8]
        path = repo_root() / "outputs" / "ch10" / f"mck_cavity_Re{int(Re)}_n{n}_{h}.npz"
        if path.is_file():
            d = np.load(path, allow_pickle=True)
            out = {k: d[k] for k in d.files}
            for k in ("t", "dt", "residual", "mass_drift", "Re", "Ma"):
                out[k] = float(out[k])
            out["steps"], out["n"] = int(out["steps"]), int(out["n"])
            out["converged"] = bool(out["converged"])
            out["history"] = [tuple(r) for r in out["history"]]
            out["snapshots"] = []
            out["cached"] = str(path)
            return out
    one = 1.0 if perturbation else 0.0
    bcf = make_cavity_bc(wall_density, corner)
    state = cavity_init(n, dtype, perturbation, U)
    nsteps = int(np.ceil(t_end / dt))
    hist, snaps = [], []
    mass0 = float(np.sum(state[0].astype(np.float64) + one))
    conv, res, k = False, np.inf, 0
    arrs = ("FB/BF", "BF/FB") if cycle else (arrangement,)
    prev = state
    for k in range(1, nsteps + 1):
        state = _wc_step_bounded(state, dt, dx, dy, Ma, Re, bcf, arrs[(k - 1) % len(arrs)], printed_step5, one, U)
        if k % check_every == 0 or k == nsteps:
            res = max(float(np.max(np.abs(state[1] - prev[1]))), float(np.max(np.abs(state[2] - prev[2])))) / (check_every * dt)
            hist.append((k * dt, res))
            prev = state
            if not np.isfinite(res):
                break
            if tol_steady and res < tol_steady:
                conv = True
                break
        if save_every and k % save_every == 0:
            snaps.append((k * dt, [s.copy() for s in state]))
    rho, m, nn = (f.astype(np.float64) for f in state)
    R = rho + one
    u, v = m / R, nn / R
    p = R / Ma ** 2
    x = np.linspace(0, 1, n + 1)
    out = dict(x=x, y=x.copy(), rho=R, u=u, v=v, p=p - p.mean(), psi=_node_streamfunction(u, dy), t=k * dt, steps=k,
               converged=conv, residual=float(res), dt=dt, mass_drift=float(np.sum(R) - mass0) / mass0, history=hist,
               snapshots=snaps, Re=Re, Ma=Ma, n=n)
    if path is not None:
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(path, **{k2: v2 for k2, v2 in out.items() if k2 not in ("snapshots", "history")},
                 history=np.array(hist, dtype=float).reshape(-1, 2))
    return out


def wc_taylor_green(n: int = 32, Ma: float = 0.05, Re: float = 100.0, t_end: float = 1.0, sigma: float = 0.8,
                    arrangement: str = "FF/BB") -> dict:
    """Weakly compressible MacCormack on the periodic box [0, 2π]² started from the Taylor–Green vortex.

    Book: §10.4 — the claim that the low-Mach isothermal equations "approximate the incompressible limit" (below (10.99)),
    tested against the exact incompressible decay e^{−2t/Re}.  Initial density from the incompressible pressure,
    ρ = 1 + Ma²p.  Returns dict(err_u (max), decay_measured (rms amplitude ratio), decay_exact, mass_drift, steps, dt).
    Validation: V3 (err falls ~ 7× per halving of Δx at fixed Ma); V4 mass conserved to round-off (periodic box, no walls).
    At fixed Δx the error *grows* as Ma falls (32²: 0.0062, 0.0099, 0.0197 at Ma = 0.1, 0.05, 0.025 — Δt ∝ Ma·Δx and the
    acoustic truncation error accumulate), so the O(Ma²) approach to the incompressible limit shows only on fine grids.
    Label: converged, conserved.
    """
    L = 2 * np.pi
    x = np.arange(n) * L / n
    X, Y = np.meshgrid(x, x, indexing="xy")
    dx = L / n
    u0, v0 = np.sin(X) * np.cos(Y), -np.cos(X) * np.sin(Y)
    rho = Ma ** 2 * 0.25 * (np.cos(2 * X) + np.cos(2 * Y))  # perturbation ρ′ = Ma² p
    state = (rho, (1 + rho) * u0, (1 + rho) * v0)
    dt = sigma / (2.0 / dx + np.sqrt(2.0) / (Ma * dx))
    nst = int(np.ceil(t_end / dt))
    dt = t_end / nst
    mass0 = float(np.sum(1 + rho))
    a = cavity_coefficients(dt, dx, dx, Ma, Re)
    for _ in range(nst):
        state = weakly_compressible_step(state, a, "periodic", arrangement)
    R = 1 + state[0]
    u = state[1] / R
    ke = float(np.sum(u ** 2 + (state[2] / R) ** 2))
    return dict(err_u=float(np.max(np.abs(u - u0 * np.exp(-2 * t_end / Re)))),
                decay_measured=float(np.sqrt(ke / float(np.sum(u0 ** 2 + v0 ** 2)))), decay_exact=float(np.exp(-2 * t_end / Re)),
                mass_drift=float(np.sum(R) - mass0) / mass0, steps=nst, dt=dt)


# ======================================================================================================================
# Square block in a channel (10.147)–(10.155)
# ======================================================================================================================
def _wall_rho(rho, u, v, side, i, j, dx, dy, M, Re):
    k = M ** 2 / Re
    if side == "front":  # fluid at i−1, i−2, i−3
        return ((4 * rho[j, i - 1] - rho[j, i - 2]) / 3.0
                + 8.0 / (9.0 * dx) * k * (-5 * u[j, i - 1] + 4 * u[j, i - 2] - u[j, i - 3])
                - 1.0 / (18.0 * dy) * k * (-(v[j + 1, i - 2] - v[j - 1, i - 2]) + 4 * (v[j + 1, i - 1] - v[j - 1, i - 1])
                                           - 3 * (v[j + 1, i] - v[j - 1, i])))  # Eq. (10.151)
    if side == "back":
        return ((4 * rho[j, i + 1] - rho[j, i + 2]) / 3.0
                - 8.0 / (9.0 * dx) * k * (-5 * u[j, i + 1] + 4 * u[j, i + 2] - u[j, i + 3])
                - 1.0 / (18.0 * dy) * k * (-(v[j + 1, i + 2] - v[j - 1, i + 2]) + 4 * (v[j + 1, i + 1] - v[j - 1, i + 1])
                                           - 3 * (v[j + 1, i] - v[j - 1, i])))  # Eq. (10.152)
    if side == "top":
        return ((4 * rho[j + 1, i] - rho[j + 2, i]) / 3.0
                - 8.0 / (9.0 * dy) * k * (-5 * v[j + 1, i] + 4 * v[j + 2, i] - v[j + 3, i])
                - 1.0 / (18.0 * dx) * k * (-(u[j + 2, i + 1] - u[j + 2, i - 1]) + 4 * (u[j + 1, i + 1] - u[j + 1, i - 1])
                                           - 3 * (u[j, i + 1] - u[j, i - 1])))  # Eq. (10.153)
    if side == "bottom":
        return ((4 * rho[j - 1, i] - rho[j - 2, i]) / 3.0
                + 8.0 / (9.0 * dy) * k * (-5 * v[j - 1, i] + 4 * v[j - 2, i] - v[j - 3, i])
                - 1.0 / (18.0 * dx) * k * (-(u[j - 2, i + 1] - u[j - 2, i - 1]) + 4 * (u[j - 1, i + 1] - u[j - 1, i - 1])
                                           - 3 * (u[j, i + 1] - u[j, i - 1])))  # Eq. (10.154)
    raise ValueError("side must be front, back, top or bottom")


class _Block:
    """Index bookkeeping of the block geometry (node-based grid, block side 1)."""

    def __init__(self, dx, H, ahead, behind, D=1.0):
        self.dx = dx
        self.H, self.ahead, self.behind = H, ahead, behind
        self.nx = int(round((ahead + D + behind) / dx)) + 1
        self.ny = int(round(H / dx)) + 1
        self.i0 = int(round(ahead / dx))
        self.i1 = int(round((ahead + D) / dx))
        self.j0 = int(round(0.5 * (H - D) / dx))
        self.j1 = int(round(0.5 * (H + D) / dx))
        solid = np.zeros((self.ny, self.nx), dtype=bool)
        solid[self.j0:self.j1 + 1, self.i0:self.i1 + 1] = True
        self.solid = solid
        self.x = np.arange(self.nx) * dx
        self.y = np.arange(self.ny) * dx
        jj = np.arange(self.j0 + 1, self.j1)
        ii = np.arange(self.i0 + 1, self.i1)
        self.front = (np.full(jj.size, self.i0), jj)
        self.back = (np.full(jj.size, self.i1), jj)
        self.bottom = (ii, np.full(ii.size, self.j0))
        self.top = (ii, np.full(ii.size, self.j1))

    @classmethod
    def from_mask(cls, mask, dx):
        js, is_ = np.where(mask)
        g = cls.__new__(cls)
        g.dx = dx
        g.ny, g.nx = mask.shape
        g.i0, g.i1, g.j0, g.j1 = int(is_.min()), int(is_.max()), int(js.min()), int(js.max())
        g.solid = mask
        return g


def block_geometry(dx: float = 0.125, H: float = 4.0, ahead: float = 8.0, behind: float = 20.0) -> dict:
    """Node counts and block indices of the channel of :func:`block_channel` (our geometry, block side 1).
    Returns dict(nx, ny, i0, i1, j0, j1, nodes_across_block).  Label: analytic."""
    g = _Block(dx, H, ahead, behind)
    return dict(nx=g.nx, ny=g.ny, i0=g.i0, i1=g.i1, j0=g.j0, j1=g.j1, nodes_across_block=g.i1 - g.i0 + 1)


def block_wall_density(state, side: str, dx: float, dy: float, Ma: float, Re: float, i=None, j=None, geo=None,
                       perturbation: bool = True) -> dict:
    """Density on a face of the block from the wall-normal momentum equation (u = v = 0 on the face).

    Book: §10.5, Eq. (10.147) (front face: ∂ρ/∂x = (Ma²/Re)(4/3 u_xx + 1/3 v_xy)), one-sided stencils (10.148)–(10.150), and
    the resulting Eqs. (10.151) front, (10.152) back, (10.153) top, (10.154) bottom (heuristic closures, the book's choice).
    Parameters: state (tuple/dict (rho, rhou, rhov); rho = ρ′ if ``perturbation`` — the relation has weights summing to 1,
    so the constant cancels), side ("front" | "back" | "top" | "bottom"), dx, dy, Ma, Re; the face nodes by index arrays
    (i, j) or a geometry ``geo`` (as stored by :func:`block_channel`; the face interior is used).
    Returns dict(i, j, rho (face values, same convention as the input)).  Validation: V2 sympy solve of the discretised
    (10.147) (analyst: residual 0); V1 exact for manufactured quadratic fields.  Label: analytic (heuristic closure).
    """
    rho, m, n = (np.asarray(f, dtype=float) for f in _tuple(state))
    R = rho + (1.0 if perturbation else 0.0)
    u, v = m / R, n / R
    if i is None or j is None:
        if geo is None:
            raise ValueError("give the face nodes (i, j) or a block geometry")
        i, j = getattr(geo, side)
    i, j = np.asarray(i), np.asarray(j)
    return dict(i=i, j=j, rho=_wall_rho(rho, u, v, side, i, j, dx, dy, Ma, Re))


def block_init(dx: float = 0.125, H: float = 4.0, ahead: float = 8.0, behind: float = 20.0, U: float = 1.0):
    """Initial state (fluid moving at U everywhere except on/in the block) and the geometry.  Returns (state, geo)."""
    geo = _Block(dx, H, ahead, behind)
    rho = np.zeros((geo.ny, geo.nx))
    m = np.full_like(rho, U)
    m[geo.solid] = 0.0
    return (rho, m, np.zeros_like(rho)), geo


def _block_bcf_factory(geo: _Block, M: float, Re: float, U: float):
    def bcf(state_n, stage_state, star, stage, dt, dx, dy, Uw, one):
        rho, m, n = (np.array(f) for f in stage_state)
        rn = state_n[0]
        rs, ms, ns = state_n if stage == "predictor" else star
        rstar = stage_state[0] if stage == "predictor" else star[0]

        def upd(base_vals, star_vals, incr):
            return base_vals - incr if stage == "predictor" else 0.5 * (base_vals + star_vals - incr)

        # inflow x = 0 is set below from the interior (zero normal gradient) — see the DEVIATION note in block_channel
        # sliding plates j = 0, ny − 1 (interior i): ∂(ρu)/∂x = U ∂ρ/∂x (centred) and one-sided ∂(ρv)/∂y, as in (10.145)
        incb = dt * U / (2 * dx) * (rs[0, 2:] - rs[0, :-2]) + dt / (2 * dy) * (-ns[2, 1:-1] + 4 * ns[1, 1:-1] - 3 * ns[0, 1:-1])
        rho[0, 1:-1] = upd(rn[0, 1:-1], rstar[0, 1:-1], incb)
        inct = dt * U / (2 * dx) * (rs[-1, 2:] - rs[-1, :-2]) - dt / (2 * dy) * (-ns[-3, 1:-1] + 4 * ns[-2, 1:-1] - 3 * ns[-1, 1:-1])
        rho[-1, 1:-1] = upd(rn[-1, 1:-1], rstar[-1, 1:-1], inct)
        # outflow x = L (interior j): ∂(ρu)/∂x = ∂(ρv)/∂x = 0; continuity then leaves −∂(ρv)/∂y (centred)
        inco = dt / (2 * dy) * (ns[2:, -1] - ns[:-2, -1])
        rho[1:-1, -1] = upd(rn[1:-1, -1], rstar[1:-1, -1], inco)
        for jj in (0, -1):  # outflow corners on the plates: backward one-sided ∂ρ/∂x
            drc = dt * U / (2 * dx) * (3 * rs[jj, -1] - 4 * rs[jj, -2] + rs[jj, -3])
            dvy = (dt / (2 * dy) * (-ns[2, -1] + 4 * ns[1, -1] - 3 * ns[0, -1]) if jj == 0
                   else -dt / (2 * dy) * (-ns[-3, -1] + 4 * ns[-2, -1] - 3 * ns[-1, -1]))
            rho[jj, -1] = upd(rn[jj, -1], rstar[jj, -1], drc + dvy)
        # inflow x = 0 (all j, incl. the plate corners): ∂ρ/∂x = 0 by one-sided second-order extrapolation
        rho[:, 0] = (4.0 * rho[:, 1] - rho[:, 2]) / 3.0
        R = rho + one
        m[:, 0] = R[:, 0] * U  # inflow u = U, v = 0
        n[:, 0] = 0.0
        m[0, :] = R[0, :] * U  # sliding plates u = U, v = 0
        m[-1, :] = R[-1, :] * U
        n[0, :] = n[-1, :] = 0.0
        m[1:-1, -1] = (4 * m[1:-1, -2] - m[1:-1, -3]) / 3.0  # outflow ∂(ρu)/∂x = 0 (one-sided second order)
        n[1:-1, -1] = (4 * n[1:-1, -2] - n[1:-1, -3]) / 3.0
        m[geo.solid] = 0.0
        n[geo.solid] = 0.0
        uu, vv = m / R, n / R
        for side in ("front", "back", "top", "bottom"):
            ii, jj = getattr(geo, side)
            rho[jj, ii] = _wall_rho(rho, uu, vv, side, ii, jj, dx, dy, M, Re)  # Eqs. (10.151)–(10.154)
        corners = [((geo.i0, geo.j0), ("front", "bottom")), ((geo.i0, geo.j1), ("front", "top")),
                   ((geo.i1, geo.j0), ("back", "bottom")), ((geo.i1, geo.j1), ("back", "top"))]
        for (ic, jc), (s1, s2) in corners:
            r1 = _wall_rho(rho, uu, vv, s1, np.array([ic]), np.array([jc]), dx, dy, M, Re)[0]
            r2 = _wall_rho(rho, uu, vv, s2, np.array([ic]), np.array([jc]), dx, dy, M, Re)[0]
            rho[jc, ic] = 0.5 * (r1 + r2)  # corners: the average of the two sides (p. 455)
        inner = geo.solid.copy()
        inner[geo.j0, :] = inner[geo.j1, :] = False
        inner[:, geo.i0] = inner[:, geo.i1] = False
        rho[inner] = 0.0 if one else 1.0
        return rho, m, n

    return bcf


def body_forces(state, mask, Re: float, Ma: float, dx: float | None = None, perturbation: bool = True) -> dict:
    """Drag and lift coefficients of the block from the pressure and viscous tractions on its faces.

    Book: §10.5 (C_D (4.107) and C_L (4.108) histories, Figs. 10.10, 10.13): C_D = F_x/(½ρ₀U²D), C_L = F_y/(½ρ₀U²D) per unit
    span.  Traction σ·n with σ = −pI + (1/Re)(∇u + ∇uᵀ − ⅔(∇·u)I), p = ρ/Ma², n the outward normal of the block; on a face
    where u = v = 0 only wall-normal derivatives survive (front: F_x = p − (4/3)u_x/Re, F_y = −v_x/Re; top: F_x = u_y/Re,
    F_y = −p + (4/3)v_y/Re; …); one-sided second-order derivatives into the fluid; trapezoid rule along each face.
    Parameters: state (tuple/dict, or the dict returned by :func:`block_channel`), mask (bool solid mask or the geometry),
    Re, Ma, dx (needed with a bare mask), perturbation.  Returns dict(CD, CL).
    Validation: V1 tractions exact on manufactured fields; the block C_D grid sequence is monotone but not asymptotic (no
    order asserted).  Label: qualitative (demonstration).
    """
    if isinstance(state, dict) and "geo" in state and mask is None:
        mask = state["geo"]
    rho, m, n = (np.asarray(f, dtype=float) for f in _tuple(state))
    one = 1.0 if perturbation else 0.0
    geo = mask if isinstance(mask, _Block) else _Block.from_mask(np.asarray(mask, bool), dx if dx is not None else state["dx"])
    R = rho + one
    u, v = m / R, n / R
    p = rho / Ma ** 2  # the constant part of ρ/Ma² cancels on the closed surface
    h = geo.dx
    jj = np.arange(geo.j0, geo.j1 + 1)
    ii = np.arange(geo.i0, geo.i1 + 1)
    wj = np.full(jj.size, h)
    wj[[0, -1]] = 0.5 * h
    wi = np.full(ii.size, h)
    wi[[0, -1]] = 0.5 * h
    Fx = Fy = 0.0
    i = geo.i0  # front face, outward normal −x, fluid at i − 1, i − 2
    ux = -(-u[jj, i - 2] + 4 * u[jj, i - 1] - 3 * u[jj, i]) / (2 * h)
    vx = -(-v[jj, i - 2] + 4 * v[jj, i - 1] - 3 * v[jj, i]) / (2 * h)
    Fx += np.sum(wj * (p[jj, i] - 4.0 / 3.0 * ux / Re))
    Fy += np.sum(wj * (-vx / Re))
    i = geo.i1  # back face, normal +x
    ux = (-u[jj, i + 2] + 4 * u[jj, i + 1] - 3 * u[jj, i]) / (2 * h)
    vx = (-v[jj, i + 2] + 4 * v[jj, i + 1] - 3 * v[jj, i]) / (2 * h)
    Fx += np.sum(wj * (-p[jj, i] + 4.0 / 3.0 * ux / Re))
    Fy += np.sum(wj * (vx / Re))
    j = geo.j1  # top face, normal +y
    uy = (-u[j + 2, ii] + 4 * u[j + 1, ii] - 3 * u[j, ii]) / (2 * h)
    vy = (-v[j + 2, ii] + 4 * v[j + 1, ii] - 3 * v[j, ii]) / (2 * h)
    Fx += np.sum(wi * (uy / Re))
    Fy += np.sum(wi * (-p[j, ii] + 4.0 / 3.0 * vy / Re))
    j = geo.j0  # bottom face, normal −y
    uy = -(-u[j - 2, ii] + 4 * u[j - 1, ii] - 3 * u[j, ii]) / (2 * h)
    vy = -(-v[j - 2, ii] + 4 * v[j - 1, ii] - 3 * v[j, ii]) / (2 * h)
    Fx += np.sum(wi * (-uy / Re))
    Fy += np.sum(wi * (p[j, ii] - 4.0 / 3.0 * vy / Re))
    return dict(CD=float(2.0 * Fx), CL=float(2.0 * Fy))


def block_channel(Re: float = 20.0, Ma: float = 0.06, dx: float = 0.125, t_end: float = 30.0, H: float = 4.0,
                  ahead: float = 8.0, behind: float = 20.0, cycle: bool = True, cache: bool = True, sigma: float = 0.8,
                  record_every: float = 0.05, U: float = 1.0, kick: float = 0.0, save_field: bool = True,
                  verbose: bool = False) -> dict:
    """Flow past a square block between two sliding plates by explicit MacCormack (weakly compressible, our geometry).

    Book: §10.5, Fig. 10.9 (a block of side D between plates sliding at the inflow speed U — the block moving through fluid at
    rest, seen in the block frame; outflow ∂(ρu)/∂x = ∂(ρv)/∂x = 0; continuity-based densities on the outer boundary;
    (10.147)–(10.154) on the block (corners averaged); FB/BF–BF/FB cycling and ρ′ = ρ − 1 (p. 455); Δt from (10.155)).
    Our geometry (the book's is private): channel height H = 4, 8 sides ahead, 20 behind, block centred.
    # DEVIATION: at the inflow the book updates ρ from continuity like the other outer boundaries; with u = U there that
    update is ρ_t = −U(−ρ₂ + 4ρ₁ − 3ρ₀)/(2Δx), a downwind one-sided advection whose self-coefficient +3U/(2Δx) grows
    exponentially (the plate corner blows up by t ≈ 1.3 in our runs).  We set the inflow density by zero-gradient
    extrapolation ρ₀ = (4ρ₁ − ρ₂)/3 (the usual subsonic-inflow treatment: the density/pressure is taken from the interior).
    Parameters: Re, Ma, dx (= dy), t_end, H, ahead, behind, cycle, cache (npz in outputs/ch10), sigma, record_every (time
    between force samples), U, kick (amplitude of a short transverse push behind the block at 2 < t < 3 to shorten the
    symmetric transient; 0 = none), save_field, verbose.
    Returns dict(t, CD, CL (histories), state (final (rho′, rhou, rhov)), x, y, u, v, rho, solid, dt, steps, geometry numbers).
    Cost: Δx = 1/8 ≈ 10 s to t = 30; 1/16 ≈ 1–2 min; finer grids script-only (cached).
    Grid study (Re = 20, mean C_D over 20 < t < 30 at Δx = 1/8, 1/16, 1/32: 4.27, 4.42, 4.53) is **non-asymptotic**:
    observed order p ≈ 0.4, GCI_fine ≈ 9–10 % (depending on the averaging window), so the Richardson value is not trustworthy and the result is not grid
    independent (the lift does go to zero).  Label: qualitative (non-asymptotic grid study; heuristic boundary closures).
    """
    state, geo = block_init(dx, H, ahead, behind, U)
    dt = maccormack_dt_asymptotic(Ma, dx, sigma)
    path = None
    if cache:
        from .project import repo_root

        tag = f"Re{Re}_Ma{Ma}_dx{dx}_t{t_end}_H{H}_{ahead}_{behind}_{cycle}_{sigma}_{record_every}_{U}_{kick}"
        h = hashlib.md5(tag.encode()).hexdigest()[:8]
        path = repo_root() / "outputs" / "ch10" / f"block_Re{int(Re)}_dx{dx:g}_{h}.npz"
        if path.is_file():
            d = np.load(path, allow_pickle=True)
            out = {k: d[k] for k in d.files}
            out["state"] = (out.pop("rho_p"), out.pop("rhou"), out.pop("rhov"))
            for k in ("dt", "Re", "Ma", "dx"):
                out[k] = float(out[k])
            out["steps"] = int(out["steps"])
            out["geo"] = geo
            out["cached"] = str(path)
            return out
    nst = int(np.ceil(t_end / dt))
    bcf = _block_bcf_factory(geo, Ma, Re, U)
    arrs = ("FB/BF", "BF/FB") if cycle else ("FF/BB",)
    rec = max(1, int(round(record_every / dt)))
    ts, cd, cl = [], [], []
    k = 0
    for k in range(1, nst + 1):
        state = _wc_step_bounded(state, dt, dx, dx, Ma, Re, bcf, arrs[(k - 1) % len(arrs)], False, 1.0, U, mask=geo.solid)
        if kick and 2.0 < k * dt < 3.0:
            jj = slice(geo.j1 - 1, geo.j1 + 2)
            ii = slice(geo.i1 + 2, geo.i1 + 5)
            state[2][jj, ii] += kick * dt
        if k % rec == 0:
            f = body_forces(state, geo, Re, Ma)
            ts.append(k * dt)
            cd.append(f["CD"])
            cl.append(f["CL"])
            if verbose and len(ts) % 200 == 0:
                print(f"  t = {k * dt:.2f}  CD = {f['CD']:.4f}  CL = {f['CL']:+.4f}", flush=True)
            if not np.isfinite(f["CD"]):
                break
    R = state[0] + 1.0
    out = dict(t=np.array(ts), CD=np.array(cd), CL=np.array(cl), state=state, dt=dt, steps=k, Re=Re, Ma=Ma, dx=dx,
               nx=geo.nx, ny=geo.ny, block=np.array([ahead, ahead + 1.0, 0.5 * (H - 1.0), 0.5 * (H + 1.0)]),
               x=geo.x, y=geo.y, u=state[1] / R, v=state[2] / R, rho=R, solid=geo.solid, geo=geo)
    if path is not None:
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(path, t=out["t"], CD=out["CD"], CL=out["CL"], rho_p=state[0], rhou=state[1], rhov=state[2], dt=dt,
                 steps=k, Re=Re, Ma=Ma, dx=dx, nx=geo.nx, ny=geo.ny, block=out["block"], x=geo.x, y=geo.y, u=out["u"],
                 v=out["v"], rho=R, solid=geo.solid)
    return out
