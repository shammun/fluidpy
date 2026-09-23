"""Constitutive law of a Newtonian fluid and the viscous dissipation rate.

Book: Kundu, Cohen & Dowling 5e, Ch. 4 §4.5, Eqs. (4.25)–(4.37), and §4.8, Eqs. (4.58)–(4.59) (read from the
rendered pages chapters/pages/ch04/p138–p141, p151).

Conventions (``knowledge/notation.md``)
---------------------------------------
* Velocity gradient ``G[i, j] = ∂u_i/∂x_j`` [1/s]; strain rate S = ½(G + Gᵀ) (3.12); S_mm = tr S = ∇·u.
* Stress τ_ij [Pa]: first index = normal of the face, second = direction of the force (§2.4); **tension positive**,
  so a fluid at rest has τ = −pδ (4.26). The viscous part is σ_ij (4.27) ("deviatoric" in the book — traceless only
  when μ_v = 0 or ∇·u = 0, since tr σ = 3μ_v∇·u).
* Viscosities: shear μ [Pa s], second coefficient λ [Pa s], bulk μ_v = λ + ⅔μ [Pa s] (text after (4.35)); Stokes'
  assumption (4.36) is μ_v = 0.
* Every function accepts trailing point axes: G may be (d, d) or (d, d, N); p, μ may be arrays over the points.
"""
from __future__ import annotations

import numpy as np

from ._util import as_scalar_if_0d
from .tensors import strain_rate_tensor

__all__ = ["static_stress", "total_stress", "isotropic_fourth_order", "linear_stress", "viscous_stress",
           "newtonian_stress", "mean_pressure", "thermodynamic_pressure_from_stress", "pressure_difference",
           "bulk_viscosity", "lam_from_bulk", "stokes_assumption_holds", "dissipation_rate", "stress_on_plane",
           "STRESS_LAB_PRESETS", "stress_lab_gradient", "shear_stress_parallel_flow", "deviatoric_part"]

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731


def _eye_like(G: np.ndarray) -> np.ndarray:
    """δ_ij shaped (d, d, 1, …) to broadcast against G of shape (d, d, …)."""
    d = G.shape[0]
    return np.eye(d).reshape((d, d) + (1,) * (G.ndim - 2))


def _trace(G: np.ndarray) -> np.ndarray:
    return np.trace(G, axis1=0, axis2=1)


# ======================================================================================================================
# §4.5 static stress, decomposition, isotropic linear law
# ======================================================================================================================
def static_stress(p, dim: int = 3) -> np.ndarray:
    """Stress in a fluid at rest: τ_ij = −p δ_ij, Eq. (4.26).

    Book: §4.5, Eq. (4.26). The only isotropic 2nd-order tensor is δ_ij (§2.7); the minus sign because normal stresses
    are tension-positive (Fig. 2.4) while pressure compresses.

    Parameters
    ----------
    p : thermodynamic pressure [Pa] (float or array over points);  dim : 2 or 3

    Returns
    -------
    tau : (dim, dim) or (dim, dim, N) [Pa]

    Validation: V1 traction on any plane = −p n (isotropic); invariant under rotation (``tensors.is_isotropic``).
    Label: analytic.
    """
    p_ = _F(p)
    return -p_ * np.eye(dim).reshape((dim, dim) + (1,) * p_.ndim)  # Eq. (4.26): τ_ij = −p δ_ij


def total_stress(p, sigma) -> np.ndarray:
    """Total stress = static part + viscous part: τ_ij = −p δ_ij + σ_ij, Eq. (4.27).

    Book: §4.5, Eq. (4.27). Parameters: p [Pa]; sigma (d, d[, N]) viscous stress [Pa]. Returns τ [Pa].
    Validation: V1 σ = 0 reduces to (4.26). Label: analytic.
    """
    s = _F(sigma)
    return -_F(p) * _eye_like(s) + s  # Eq. (4.27)


def isotropic_fourth_order(lam, mu, gam, dim: int = 3) -> np.ndarray:
    """The general isotropic fourth-order tensor K_ijmn = λ δ_ij δ_mn + μ δ_im δ_jn + γ δ_in δ_jm, Eq. (4.29).

    Book: §4.5, Eq. (4.29) (cited from Aris 1962; we check isotropy numerically instead of proving it).

    Parameters
    ----------
    lam, mu, gam : the scalars λ, μ, γ [Pa s] (γ here is the third Lamé-type coefficient, not a ratio of specific
        heats — ``knowledge/notation.md``);  dim : 2 or 3

    Returns
    -------
    K : (dim,)*4 array [Pa s]

    Notes
    -----
    Contracting with a *symmetric* S gives σ_ij = λS_mmδ_ij + (μ + γ)S_ij (4.28), already symmetric for any γ: only the
    sum μ + γ is physical; the book's "γ = μ" (4.30) names it 2μ (it is forced if K itself must be symmetric in i, j).
    Validation: V7 invariant under 50 random rotations (``tensors.is_isotropic`` on order 4, 1e-12); V2 sympy K:S for
    symmetric S depends on μ + γ only. Label: analytic, symbolic.
    """
    d = np.eye(dim)
    return (float(lam) * np.einsum("ij,mn->ijmn", d, d) + float(mu) * np.einsum("im,jn->ijmn", d, d)
            + float(gam) * np.einsum("in,jm->ijmn", d, d))  # Eq. (4.29)


def linear_stress(K, S) -> np.ndarray:
    """Most general linear viscous stress: σ_ij = K_ijmn S_mn, Eq. (4.28).

    Book: §4.5, Eq. (4.28). Parameters: K (d,d,d,d) [Pa s]; S (d,d[,N]) strain rate [1/s]. Returns σ [Pa].
    Validation: V1 with K from :func:`isotropic_fourth_order` (λ, μ, μ) equals :func:`viscous_stress` with
    μ_v = λ + ⅔μ. Label: analytic.
    """
    return np.einsum("ijmn,mn...->ij...", _F(K), _F(S))  # Eq. (4.28): σ_ij = K_ijmn S_mn


def viscous_stress(G, mu, mu_v=0.0) -> np.ndarray:
    """Newtonian viscous stress σ_ij = μ(∂u_i/∂x_j + ∂u_j/∂x_i) + (μ_v − ⅔μ)(∂u_m/∂x_m) δ_ij, Eq. (4.59).

    Book: §4.5, Eq. (4.37) (viscous part) and §4.8, Eq. (4.59) (written out). Equivalent to 2μS_ij + λS_mmδ_ij of
    (4.30)–(4.31) with λ = μ_v − ⅔μ.

    Parameters
    ----------
    G : velocity gradient ∂u_i/∂x_j [1/s], (d, d) or (d, d, N);  mu : shear viscosity μ [Pa s];
    mu_v : bulk viscosity μ_v [Pa s] (0 = Stokes' assumption (4.36))

    Returns
    -------
    sigma : (d, d[, N]) [Pa], symmetric.

    Assumptions: isotropic Newtonian fluid (linear in S, depends on S only — rigid rotation gives 0).
    Validation: V1 parallel shear u = (γy, 0, 0) gives σ₁₂ = μγ = ``ch01.newton_shear_stress`` (link to (1.3)); rigid
    rotation → 0; V2 sympy (4.31) ≡ (4.37) ≡ (4.59). Label: analytic, symbolic.
    """
    G_ = _F(G)
    S = strain_rate_tensor(G_) if G_.ndim == 2 else 0.5 * (G_ + np.swapaxes(G_, 0, 1))
    div = _trace(S)
    third = 1.0 / 3.0  # the book's 3-D deviator; a 2 × 2 G is a plane flow of a 3-D fluid (u₃ = 0, S₃₃ = 0)
    mu_, mv = _F(mu), _F(mu_v)
    return 2.0 * mu_ * (S - third * div * _eye_like(S)) + mv * div * _eye_like(S)  # Eq. (4.37)/(4.59)


def newtonian_stress(G, p=0.0, mu=1.0e-3, lam=None, mu_v=None, incompressible: bool = False,
                     tol: float = 1e-9) -> np.ndarray:
    """Full stress tensor of a Newtonian fluid, Eqs. (4.31), (4.35), (4.37).

    Book: §4.5 — (4.31) τ_ij = −pδ_ij + 2μS_ij + λS_mmδ_ij; (4.37) τ_ij = −pδ_ij + 2μ(S_ij − ⅓S_mmδ_ij) + μ_vS_mmδ_ij
    with μ_v = λ + ⅔μ; incompressible (4.35) τ_ij = −pδ_ij + 2μS_ij.

    Parameters
    ----------
    G : velocity gradient [1/s], (d, d) or (d, d, N) (nested lists accepted)
    p : thermodynamic (or, if incompressible, mechanical) pressure [Pa]
    mu : shear viscosity [Pa s]
    lam : second viscosity coefficient λ [Pa s] — give λ **or** μ_v (at most one; both raises ``ValueError``);
        neither → Stokes' assumption (4.36) μ_v = 0, i.e. λ = −⅔μ
    mu_v : bulk viscosity [Pa s]
    incompressible : use (4.35) and raise ``ValueError`` if |tr G| > tol·‖G‖
    tol : relative tolerance of that check [-]

    Returns
    -------
    tau : (d, d[, N]) [Pa]

    Validation: V1 parallel shear τ₁₂ = μγ (``newtonian_stress([[0, 10, 0], [0, 0, 0], [0, 0, 0]], mu=1e-3)[0, 1]``
    = 0.010 Pa); at rest → (4.26); V2 sympy (4.31) ≡ (4.37); traces (4.32)–(4.34). Label: analytic, symbolic.
    """
    G_ = _F(G)
    if lam is not None and mu_v is not None:
        raise ValueError("give lam or mu_v, not both (mu_v = lam + 2 mu/3)")
    if incompressible:
        tr = np.abs(_trace(G_))
        if np.any(tr > tol * max(float(np.linalg.norm(G_)), 1e-300)):
            raise ValueError(f"incompressible=True but tr G = ∇·u = {np.max(tr):.3e} ≠ 0")
        S = 0.5 * (G_ + np.swapaxes(G_, 0, 1))
        return -_F(p) * _eye_like(S) + 2.0 * _F(mu) * S  # Eq. (4.35)
    if lam is not None:
        S = 0.5 * (G_ + np.swapaxes(G_, 0, 1))
        return -_F(p) * _eye_like(S) + 2.0 * _F(mu) * S + _F(lam) * _trace(S) * _eye_like(S)  # Eq. (4.31)
    mv = 0.0 if mu_v is None else mu_v
    return total_stress(p, viscous_stress(G_, mu, mv))  # Eq. (4.37)


def _trace3(T: np.ndarray, tau33):
    """τ_ii summed over all three directions: a 3 × 3 tensor as is; a 2 × 2 (plane-flow) tensor needs τ₃₃ explicitly."""
    if T.shape[0] == 3 and T.shape[1] == 3:
        if tau33 is not None:
            raise ValueError("tau33 is only for a 2 × 2 (plane-flow) tensor; a 3 × 3 tensor already contains τ₃₃")
        return _trace(T)
    if T.shape[0] == 2 and T.shape[1] == 2:
        if tau33 is None:
            raise ValueError("(4.32)–(4.34) need the full trace τ₁₁ + τ₂₂ + τ₃₃: for a 2 × 2 (plane-flow) tensor pass "
                             "tau33 (e.g. −p + (μ_v − ⅔μ)∇·u for a Newtonian fluid with u₃ = 0) or use the 3 × 3 tensor")
        return _trace(T) + _F(tau33)
    raise ValueError(f"tau must be 3 × 3 (or 2 × 2 with tau33); got shape {T.shape[:2]}")


def mean_pressure(tau, tau33=None):
    """Mean (mechanical) pressure p̄ = −⅓ τ_ii (sum over i = 1, 2, 3), Eq. (4.33).

    Book: §4.5, Eq. (4.33). Parameters: tau (3, 3[, N]) [Pa]; for a 2 × 2 plane-flow tensor pass ``tau33`` [Pa] (the
    out-of-plane normal stress, which is not zero: τ₃₃ = −p + (μ_v − ⅔μ)∇·u when u₃ = 0) — otherwise ``ValueError``
    (a 2 × 2 trace alone would give a wrong p̄). Returns p̄ [Pa].
    Validation: V1 τ = −pδ gives p; V2 (4.34) p − p̄ = μ_v∇·u (G = diag(1, 0), p = 100, μ = 1, μ_v = 0.5: 0.5 Pa with
    the 3 × 3 tensor or with tau33). Label: analytic.
    """
    T = _F(tau)
    return as_scalar_if_0d(-_trace3(T, tau33) / 3.0)  # Eq. (4.33)


def thermodynamic_pressure_from_stress(tau, div_u, mu, lam, tau33=None):
    """Thermodynamic pressure recovered from the trace of the stress, p = −⅓τ_ii + (⅔μ + λ)∇·u, Eq. (4.32).

    Book: §4.5, Eq. (4.32) (the trace of (4.31) with δ_ii = 3). Parameters: tau [Pa] (3 × 3, or 2 × 2 with ``tau33`` as
    in :func:`mean_pressure`); div_u [1/s]; mu, lam [Pa s]. Returns p [Pa].
    Validation: V1 round trip with :func:`newtonian_stress`. Label: analytic.
    """
    return as_scalar_if_0d(mean_pressure(tau, tau33) + (2.0 / 3.0 * _F(mu) + _F(lam)) * _F(div_u))  # Eq. (4.32)


def pressure_difference(div_u, mu=None, lam=None, mu_v=None):
    """Thermodynamic minus mean pressure p − p̄ = (⅔μ + λ)∇·u = μ_v ∇·u, Eq. (4.34).

    Book: §4.5, Eq. (4.34). Give μ_v, or μ and λ. Returns [Pa]; zero for ∇·u = 0 or μ_v = 0 (Stokes (4.36)).
    Validation: V1 equals p − ``mean_pressure(newtonian_stress(...))``. Label: analytic.
    """
    if mu_v is None:
        if mu is None or lam is None:
            raise ValueError("give mu_v, or both mu and lam")
        mu_v = _F(lam) + 2.0 / 3.0 * _F(mu)
    return as_scalar_if_0d(_F(mu_v) * _F(div_u))  # Eq. (4.34)


def bulk_viscosity(lam, mu):
    """Bulk viscosity μ_v = λ + ⅔μ [Pa s] (text after (4.35)); Stokes' assumption (4.36) sets it to 0.

    Book: §4.5, text after Eq. (4.35), Eq. (4.36). Label: analytic.
    """
    return as_scalar_if_0d(_F(lam) + 2.0 / 3.0 * _F(mu))


def lam_from_bulk(mu_v, mu):
    """Second viscosity coefficient λ = μ_v − ⅔μ [Pa s] (inverse of :func:`bulk_viscosity`).

    Book: §4.5, text after Eq. (4.35). Label: analytic.
    """
    return as_scalar_if_0d(_F(mu_v) - 2.0 / 3.0 * _F(mu))


def stokes_assumption_holds(lam, mu, tol: float = 1e-12) -> bool:
    """True if λ + ⅔μ = 0 (Stokes' assumption, Eq. (4.36)) to within tol·|μ|.

    Book: §4.5, Eq. (4.36). Label: analytic.
    """
    return bool(abs(float(lam) + 2.0 / 3.0 * float(mu)) <= tol * max(abs(float(mu)), 1e-300))


def shear_stress_parallel_flow(mu, dudy):
    """Shear stress τ₁₂ = μ du/dy of a parallel flow u(y), the special case of (4.37) that recovers Newton's (1.3).

    Book: §4.5, text after Eq. (4.37) ("τ₁₂ = μ(∂u₁/∂x₂ + ∂u₂/∂x₁)"). Units: μ [Pa s], du/dy [1/s] → [Pa].
    Validation: V1 equals ``newtonian_stress`` off-diagonal for G = [[0, γ, 0], …]. Label: analytic.
    """
    return as_scalar_if_0d(_F(mu) * _F(dudy))


# ======================================================================================================================
# §4.8 viscous dissipation
# ======================================================================================================================
def dissipation_rate(G, rho, mu, mu_v=0.0, form: str = "sum_of_squares"):
    """Kinetic-energy dissipation rate per unit mass ε = (1/ρ)σ_ij S_ij, Eq. (4.58).

    Book: §4.8, Eq. (4.58): ε ≡ (1/ρ)σ_ijS_ij = 2ν(S_ij − ⅓S_mmδ_ij)² + (μ_v/ρ)S_mm², with σ from (4.59).

    Parameters
    ----------
    G : velocity gradient [1/s], (d, d) or (d, d, N);  rho : density [kg/m³];  mu, mu_v : viscosities [Pa s]
    form : "sum_of_squares" (right side of (4.58)), "contraction" ((1/ρ)σ_ijS_ij) or "both" (returns the pair)

    Returns
    -------
    eps : [W/kg] = [m²/s³] (≥ 0 whenever μ, μ_v ≥ 0 — a sum of squares).

    Notes
    -----
    Only S enters: σ_ij ∂u_j/∂x_i = σ_ij(S_ij + R_ji/2…) = σ_ijS_ij because σ is symmetric (text after (4.57)).
    Validation: V1 both routes agree on 1000 random G (1e-12); parallel shear ε = νγ²; rigid rotation 0; V7 invariant
    under rotations; μ < 0 can make ε < 0 (the second-law discrimination of (4.63)). Label: analytic.
    """
    G_ = _F(G)
    S = 0.5 * (G_ + np.swapaxes(G_, 0, 1))
    rho_ = _F(rho)
    third = 1.0 / 3.0  # 3-D deviator also for a plane (2 × 2) G: S₃₃ = 0 but the ⅓S_mmδ₃₃ part still counts
    div = _trace(S)
    dev = S - third * div * _eye_like(S)
    dev33 = (-third * div) if G_.shape[0] == 2 else 0.0  # the missing (3, 3) entry of the deviator in plane flow
    sq = (2.0 * _F(mu) * (np.sum(dev * dev, axis=(0, 1)) + dev33 ** 2)
          + _F(mu_v) * div ** 2) / rho_  # Eq. (4.58), right side
    if form == "sum_of_squares":
        return as_scalar_if_0d(sq)
    sig = viscous_stress(G_, mu, mu_v)
    con = np.sum(sig * S, axis=(0, 1)) / rho_  # Eq. (4.58): ε = (1/ρ) σ_ij S_ij
    if form == "contraction":
        return as_scalar_if_0d(con)
    if form == "both":
        return as_scalar_if_0d(sq), as_scalar_if_0d(con)
    raise ValueError("form must be 'sum_of_squares', 'contraction' or 'both'")


# ======================================================================================================================
# E3 helpers: the traction on a rotatable plane, preset velocity gradients
# ======================================================================================================================
STRESS_LAB_PRESETS = ("shear", "extension", "rotation", "expansion")


def stress_lab_gradient(name: str, rate: float = 1.0, dim: int = 3) -> np.ndarray:
    """Velocity gradient G [1/s] of E3's presets (our examples, not the book's).

    * "shear": u = (γy, 0, 0) → G₁₂ = rate;  * "extension": planar pure strain diag(rate, −rate, 0) (∇·u = 0);
    * "rotation": rigid rotation about z, G₁₂ = −rate, G₂₁ = rate (S = 0, ω₃ = 2·rate);
    * "expansion": isotropic expansion rate·δ (∇·u = dim·rate — the only preset where μ_v matters).

    Book: §4.5 (the chain G → S → σ of (4.28)–(4.37)); presets are ours. Label: analytic.
    """
    G = np.zeros((dim, dim))
    r = float(rate)
    if name == "shear":
        G[0, 1] = r
    elif name == "extension":
        G[0, 0], G[1, 1] = r, -r
    elif name == "rotation":
        G[0, 1], G[1, 0] = -r, r
    elif name == "expansion":
        G = r * np.eye(dim)
    else:
        raise ValueError(f"unknown preset {name!r}; choose from {STRESS_LAB_PRESETS}")
    return G


def stress_on_plane(G, p=0.0, mu=1.0e-3, mu_v=0.0, theta=0.0):
    """Normal and (signed) shear stress on a plane with normal n = (cos θ, sin θ, 0) in a Newtonian fluid:
    f_j = n_iτ_ij with τ from (4.37); σ_n = f·n (tension positive), τ_s = f·t with t = (−sin θ, cos θ, 0).

    Book: §4.5, Eq. (4.37) (τ from the strain rate) and §2.6, Eq. (2.15) f_j = n_iτ_ij (Cauchy) — E3's rotatable plane.

    Parameters
    ----------
    G : velocity gradient [1/s], (d, d) (nested lists accepted) or a preset name of :data:`STRESS_LAB_PRESETS`
        (then with rate 1 s⁻¹)
    p : pressure [Pa];  mu, mu_v : viscosities [Pa s]
    theta : angle of the plane's normal from the x axis in the x–y plane [rad]

    Returns
    -------
    (sigma_n, tau_s) [Pa], plain floats; τ_s is signed (positive along t, counter-clockwise of n).

    Validation: V1 shear G₁₂ = 10 s⁻¹, μ = 1e-3, θ = π/4 → (0.010, 0.0); θ = 0 → (−p, μγ); rotation preset: τ_s = 0,
    σ_n = −p for every θ; expansion: σ_n = −p + 3μ_v·rate for every θ (isotropic). Label: analytic.
    """
    if isinstance(G, str):
        G = stress_lab_gradient(G)
    G_ = _F(G)
    d = G_.shape[0]
    th = float(theta)
    n_ = np.array([np.cos(th), np.sin(th), 0.0][:d])
    t_ = np.array([-np.sin(th), np.cos(th), 0.0][:d])
    tau = newtonian_stress(G_, p, mu, mu_v=mu_v)
    f = np.einsum("ij,i->j", tau, n_)  # Eq. (2.15): f_j = n_i τ_ij
    return float(f @ n_), float(f @ t_)


def deviatoric_part(A) -> np.ndarray:
    """Deviator A − ⅓(tr A)δ of a 3 × 3 tensor (trailing point axes allowed); traceless by construction.

    Book: §4.5, Eq. (4.37) (S_ij − ⅓S_mmδ_ij) and §4.8, Eq. (4.58) (its square gives ε). Label: analytic.
    """
    A_ = _F(A)
    return A_ - _trace(A_) / 3.0 * _eye_like(A_)
