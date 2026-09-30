"""Bluff-body flows (Ch. 9, §9.7–§9.9): cylinder and sphere regimes, the Kármán vortex street, Strouhal shedding, the separated-flow
pressure-drag model, the drag crisis and the dynamics of sports balls.

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e, §9.7 (streamlining, form drag), §9.8 (cylinder, Figs. 9.16–9.21, Kármán street,
Strouhal number (4.102), drag crisis), §9.9 (sphere, cricket/tennis/baseball).  Most of this section of the book is descriptive; the
quantitative parts coded here are (i) the point-vortex street (b/a = arccosh(√2)/π = 0.2805 and its linear stability), (ii) the Strouhal
relation, (iii) our separated-pressure model (ideal C_p = 1 − 4 sin²φ to the separation angle, constant base pressure behind), (iv) sports-ball
kinematics y = ½at².  Regime thresholds, separation angles and base pressures are ROUNDED experimental values (arguments with defaults, flagged
'qualitative', labelled illustrative); they are not derived and not tested for accuracy.

Convention (analysis §9): ``phi`` / ``phi_deg`` is measured from the FORWARD stagnation point (ch06's θ is from +x downstream; ch08's sphere θ from
the downstream axis).  Re = U∞d/ν on the DIAMETER.  Pressure coefficient C_p = (p − p∞)/(½ρU∞²).
"""
from __future__ import annotations

import warnings

import numpy as np
from scipy.linalg import eig, expm

from ._util import as_scalar_if_0d

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d

# ----------------------------------------------------------------------------------------------------------------------
# regimes (rounded; qualitative)
# ----------------------------------------------------------------------------------------------------------------------
CYLINDER_THRESHOLDS = dict(creeping=1.0, attached_eddies=4.0, street_onset=40.0, irregular=200.0, turbulent_wake=3000.0,
                           critical=3e5, supercritical=6e5)
"""Rounded regime boundaries on the cylinder Reynolds number (Fig. 9.16–9.21 discussion); experiments give ranges, not sharp values.  Label: illustrative."""

SPHERE_THRESHOLDS = dict(creeping=1.0, steady_wake=130.0, critical=3e5, supercritical=8e5)
"""Rounded regime boundaries on the sphere Reynolds number (§9.9); illustrative."""


def cylinder_flow_regime(Re, thresholds: dict | None = None) -> dict:
    """Flow regime of a smooth circular cylinder at Re = U∞d/ν (ROUNDED thresholds — qualitative, illustrative).

    Book: §9.8 (Figs. 9.16–9.21): the book's rounded values, experimental.  Labels (ASCII, exact): Re < 1 "creeping flow: symmetric, no wake";
    1 ≤ Re < 4 "creeping to attached eddies"; 4 ≤ Re < 40 "two steady attached eddies"; 40 ≤ Re < 200 "laminar Karman street"; 200 ≤ Re < 3000
    "irregular vortices, St stays near 0.2"; 3000 ≤ Re < 3e5 "subcritical: laminar separation, wide wake"; 3e5 ≤ Re < 6e5 "critical: drag crisis";
    Re ≥ 6e5 "supercritical: turbulent separation, narrow wake".
    Parameters: Re [–]; thresholds [dict] overrides entries of :data:`CYLINDER_THRESHOLDS` (keys creeping, attached_eddies, street_onset, irregular,
    turbulent_wake, critical, supercritical).
    Returns dict(label, thresholds, separation_deg (82.0 subcritical, 125.0 supercritical, else None; from the forward stagnation point),
    St (0.2 for 40 ≤ Re < 3000, else None)).
    Validation: V7 table lookup at each boundary (a value on a threshold belongs to the upper regime).  Label: qualitative.
    """
    th = dict(CYLINDER_THRESHOLDS)
    if thresholds:
        th.update(thresholds)
    R = float(Re)
    sep, st = None, None
    if R < th["creeping"]:
        lab = "creeping flow: symmetric, no wake"
    elif R < th["attached_eddies"]:
        lab = "creeping to attached eddies"
    elif R < th["street_onset"]:
        lab = "two steady attached eddies"
    elif R < th["irregular"]:
        lab, st = "laminar Karman street", 0.2
    elif R < th["turbulent_wake"]:
        lab, st = "irregular vortices, St stays near 0.2", 0.2
    elif R < th["critical"]:
        lab, sep = "subcritical: laminar separation, wide wake", 82.0
    elif R < th["supercritical"]:
        lab = "critical: drag crisis"
    else:
        lab, sep = "supercritical: turbulent separation, narrow wake", 125.0
    return dict(label=lab, thresholds=th, separation_deg=sep, St=st)


def sphere_flow_regime(Re, thresholds: dict | None = None, sep_sub_deg: float = 80.0, sep_super_deg: float = 120.0) -> dict:
    """Flow regime of a smooth sphere at Re = U∞d/ν (ROUNDED, qualitative, illustrative).

    Book: §9.9 (Figs. 9.22–9.23): the book transplants the cylinder range for the attached eddy (sphere separation actually starts near Re ≈ 20 — slip R13);
    loops are shed above Re ≈ 130 but with NO regular street; drag crisis near 3×10⁵.  Labels (ASCII, exact): Re < 1 "creeping flow: symmetric, no wake";
    1 ≤ Re < 130 "steady wake with an attached doughnut eddy"; 130 ≤ Re < 3e5 "unsteady wake, loops shed (no regular street)"; 3e5 ≤ Re < 8e5
    "critical: drag crisis"; Re ≥ 8e5 "supercritical: turbulent separation, narrow wake".
    Parameters: Re [–]; thresholds overrides :data:`SPHERE_THRESHOLDS`; sep_sub_deg, sep_super_deg [deg] separation angles from the forward stagnation
    point before / after the crisis (rounded experimental, illustrative).
    Returns dict(label, thresholds, separation_deg (sep_sub_deg / sep_super_deg or None), St (None: no regular shedding)).  Label: qualitative.
    """
    th = dict(SPHERE_THRESHOLDS)
    if thresholds:
        th.update(thresholds)
    R = float(Re)
    sep = None
    if R < th["creeping"]:
        lab = "creeping flow: symmetric, no wake"
    elif R < th["steady_wake"]:
        lab = "steady wake with an attached doughnut eddy"
    elif R < th["critical"]:
        lab, sep = "unsteady wake, loops shed (no regular street)", float(sep_sub_deg)
    elif R < th["supercritical"]:
        lab = "critical: drag crisis"
    else:
        lab, sep = "supercritical: turbulent separation, narrow wake", float(sep_super_deg)
    return dict(label=lab, thresholds=th, separation_deg=sep, St=None)


# ----------------------------------------------------------------------------------------------------------------------
# Kármán street
# ----------------------------------------------------------------------------------------------------------------------


def karman_street_ratio() -> float:
    """Stable spacing ratio b/a = arccosh(√2)/π = 0.2805 of Kármán's staggered double row of point vortices (b = distance between rows, a = spacing in a row).

    Book: §9.8 — 'a staggered row is stable only if the lateral/longitudinal spacing is 0.28' (Kármán 1912); closed form cosh(πb/a) = √2.
    Returns b/a [–].  Validation: V1/V2 cosh(π·0.2805) = √2; literature 0.281 (Horváth 2020 JGR).  Label: analytic."""
    return float(np.arccosh(np.sqrt(2.0)) / np.pi)


def karman_street_velocity(a: float, b: float, Gamma: float) -> float:
    """Self-induced speed of a staggered street, U = (Γ/2a)tanh(πb/a)  [m/s] (Lamb §156).  a, b [m]; Gamma [m²/s] (magnitude of each circulation).
    Derived from the complex velocity of a row, u − iv = (Γ/2ia)cot(π(z − z₀)/a) (ch05 point vortices).  Label: analytic."""
    return _S(Gamma / (2.0 * a) * np.tanh(np.pi * _F(b) / a))


def _lattice_T(dz, k: float, a: float, N: int):
    """T(dz; k) = Σ_m e^{ikma}/(dz − ma)² over the integers m (dz complex, not on the real axis).

    Exact for k = 0: (π/a)²/sin²(πdz/a); exact for k = ±π/a: (π/a)² cos(πdz/a)/sin²(πdz/a) (derivatives of π/(a sin) and π cot/a); otherwise the sum is
    truncated at |m| ≤ N (the neglected oscillating tail is O(1/(N²·ka)); it is NOT corrected, so near k → 0 the error grows like 2/(N a²))."""
    th = float(k) * a
    if abs(th) < 1e-12:
        return (np.pi / a) ** 2 / np.sin(np.pi * dz / a) ** 2
    if abs(abs(th) - np.pi) < 1e-12:
        return (np.pi / a) ** 2 * np.cos(np.pi * dz / a) / np.sin(np.pi * dz / a) ** 2
    m = np.arange(-int(N), int(N) + 1)
    return np.sum(np.exp(1j * th * m) / (dz - m * a) ** 2)


def _same_row_sum(k: float, a: float):
    """Σ_{m≠0} e^{ikma}/(ma)² = (2/a²)Σ_{m≥1}cos(mθ)/m² = (2/a²)(π²/6 − πθ/2 + θ²/4), θ = ka mod 2π ∈ [0, 2π]  (exact, all k)."""
    th = np.mod(abs(float(k)) * a, 2.0 * np.pi)
    return 2.0 / a ** 2 * (np.pi ** 2 / 6.0 - np.pi * th / 2.0 + th ** 2 / 4.0)


def _street_bloch_matrix(b_over_a: float, offset: float, k: float, Gamma: float, a: float, N: int) -> np.ndarray:
    """4×4 complex matrix of the linearised double row for a perturbation with Bloch wavenumber k.

    Rows: upper (y = +b/2, circulation −Γ, vortices at na + ib/2) and lower (y = −b/2, +Γ, at (n + offset)a − ib/2).  Perturbations ζ_{r,n} = P_r e^{ikna},
    ζ̄_{r,n} = Q_r e^{ikna} (independent amplitudes of ζ and its conjugate), state (P₁, Q₁, P₂, Q₂).  From dż̄ = (1/2πi)ΣΓ_j/(z_k − z_j):
    Ṗ_r = (Γ_r/2πi)[E(0) − E(k)]Q_r + (Γ_s/2πi)[C(0)Q_r − C(k)Q_s],  Q̇_r = −(Γ_r/2πi)[E(0) − E(k)]P_r − (Γ_s/2πi)[T(0)P_r − T(k)P_s],
    with E the same-row sums, T(k) = Σ_m e^{ikma}/(dz_rs − ma)², C(k) = conj(T(−k)), dz_rs = z_r0 − z_s0."""
    b = b_over_a * a
    z0 = [0.5j * b, offset * a - 0.5j * b]
    G = [-Gamma, Gamma]
    E0, Ek = _same_row_sum(0.0, a), _same_row_sum(k, a)
    M = np.zeros((4, 4), dtype=complex)
    c = 1.0 / (2j * np.pi)
    for r in (0, 1):
        s = 1 - r
        dz = z0[r] - z0[s]
        Tk, T0 = _lattice_T(dz, k, a, N), _lattice_T(dz, 0.0, a, N)
        Ck, C0 = np.conj(_lattice_T(dz, -k, a, N)), np.conj(T0)
        same = G[r] * c * (E0 - Ek)
        M[2 * r, 2 * r + 1] = same + G[s] * c * C0
        M[2 * r, 2 * s + 1] = -G[s] * c * Ck
        M[2 * r + 1, 2 * r] = -(same + G[s] * c * T0)
        M[2 * r + 1, 2 * s] = G[s] * c * Tk
    return M


def _sort_ev(ev: np.ndarray) -> np.ndarray:
    return ev[np.lexsort((ev.imag, ev.real))]


def karman_street_spectrum(b_over_a: float, offset: float = 0.5, k: float | None = None, Gamma: float = 1.0, a: float = 1.0,
                           n_terms: int = 2000) -> np.ndarray:
    """The four eigenvalues σ [1/s] of the linearised infinite double row of point vortices for one perturbation wavenumber k (perturbation ∝ e^{σt + ikx}).

    Book: §9.8 (Kármán's point-vortex model; his result b/a = 0.2805).  Method (ours, following Lamb Art. 156): displace vortex n of each row by
    P_r e^{ikna}; linearise the induced velocity (Γ/2πi)/(z_k − z_j); the sums over the infinite row are lattice sums — exact for k = 0 and k = π/a and for the
    same-row part at any k, truncated at ±``n_terms`` otherwise.
    Parameters: b_over_a [–] row separation / streamwise spacing; offset [–] streamwise shift of the lower row in units of a (0.5 = Kármán's staggered street,
    0 = the symmetric, always-unstable row); k [1/m] wavenumber (default π/a, the alternate-vortex mode; k and k + 2π/a are the same mode); Gamma [m²/s];
    a [m]; n_terms lattice truncation for k ∉ {0, π/a}.
    Returns complex ndarray(4) [1/s] sorted by (real part, imaginary part).  For k = π/a, offset 0.5: σ² = ±… with growth rate
    (πΓ/2a²)|½ − sech²(πb/a)| (see :func:`karman_street_growth_closed`); at b/a = arccosh(√2)/π the spectrum is ±0.7854i·Γ/a² (neutral).
    Validation: V1 closed form at k = π/a (1e-12); V3 lattice truncation; independent cross-check against the finite periodic-cell eigenproblem
    :func:`karman_street_spectrum_periodic`.  Label: analytic, converged.
    """
    kk = np.pi / a if k is None else float(k)
    return _sort_ev(np.linalg.eigvals(_street_bloch_matrix(b_over_a, offset, kk, Gamma, a, n_terms)))


def karman_street_growth(b_over_a: float, offset: float = 0.5, Gamma: float = 1.0, a: float = 1.0, k: float | None = None, n_terms: int = 2000) -> float:
    """Largest growth rate max Re σ [1/s] of the double row of point vortices (≥ 0 up to noise of order 1e-8·Γ/a² at exactly neutral spacing).

    Book: §9.8.  ``k`` [1/m]: one wavenumber, or ``None`` = the maximum of Re σ over 61 wavenumbers in (0, π/a] (the most unstable mode of the street).
    Other parameters as :func:`karman_street_spectrum`.  Expect (Γ = a = 1, offset 0.5): 0.6400 (b/a = 0.1), 0.2982 (0.2), 0.1099 (0.25), ≈ 0 (0.28055),
    0.0663 (0.3), 0.2208 (0.35), 0.5359 (0.5), 0.7737 (1.0); non-staggered rows (offset 0): π/4 for every b/a.
    Validation: V1 closed form :func:`karman_street_growth_closed`; V2 :func:`ch09_boundary_layers.karman_street_sympy`.  Label: analytic, converged."""
    if k is not None:
        return float(np.max(karman_street_spectrum(b_over_a, offset, k, Gamma, a, n_terms).real))
    ks = np.linspace(np.pi / a / 61.0, np.pi / a, 61)
    return float(max(np.max(karman_street_spectrum(b_over_a, offset, kk, Gamma, a, n_terms).real) for kk in ks))


def karman_street_growth_closed(b_over_a: float, Gamma: float = 1.0, a: float = 1.0):
    """Closed-form growth rate of the alternate-vortex mode (k = π/a) of the staggered street: (πΓ/2a²)|½ − sech²(πb/a)|  [1/s].

    Book: §9.8 (Kármán's stability analysis).  Zero at cosh²(πb/a) = 2, i.e. b/a = arccosh(√2)/π = 0.28055 (:func:`karman_street_ratio`).
    b_over_a [–]; Gamma [m²/s]; a [m].  Derivation: the 4×4 matrix of :func:`karman_street_spectrum` at k = π/a factorises as
    ((μ − γ)² + σ²)((μ + γ)² + σ²) — see :func:`ch09_boundary_layers.karman_street_sympy`.  Label: analytic."""
    return _S(np.pi * Gamma / (2.0 * a ** 2) * np.abs(0.5 - 1.0 / np.cosh(np.pi * _F(b_over_a)) ** 2))


def _street_matrix(b_over_a: float, N: int, Gamma: float, a: float, stagger: float):
    """Real 4N×4N linearisation of a double row of point vortices (period cell L = Na, N vortices per row per cell)."""
    b = b_over_a * a
    L = N * a
    n = np.arange(N)
    z = np.concatenate([n * a + 0.5j * b, (n + stagger) * a - 0.5j * b])
    G = np.concatenate([-Gamma * np.ones(N), Gamma * np.ones(N)])
    dz = z[:, None] - z[None, :]
    np.fill_diagonal(dz, 1.0)
    csc2 = 1.0 / np.sin(np.pi * dz / L) ** 2
    c = G / (2j * L)
    A = -(c[None, :]) * (np.pi / L) * csc2  # A_kj = −c_j (π/L) csc²(π(z_k − z_j)/L)
    np.fill_diagonal(A, 0.0)
    Ac = np.conj(A)
    Q = -Ac
    Q[np.diag_indices(2 * N)] = Ac.sum(axis=1)  # dζ_k/dt = Σ conj(A_kj)(ζ̄_k − ζ̄_j) = (Qζ̄)_k
    Qr, Qi = Q.real, Q.imag
    return np.block([[Qr, Qi], [Qi, -Qr]]), z, G


def karman_street_spectrum_periodic(b_over_a: float, n_pairs: int = 16, Gamma: float = 1.0, a: float = 1.0, stagger: float = 0.5,
                                    fast: bool = False) -> np.ndarray:
    """All eigenvalues σ [1/s] of the linearised motion of a double row of point vortices, perturbations periodic over ``n_pairs`` spacings (finite-cell method).

    Book: §9.8; method (ours, following Lamb Art. 156): the velocity of a periodic row is (Γ/2iL)cot(π(z−z₀)/L), linearised about the street, eigenvalues by
    ``scipy.linalg.eig`` of the real 4N×4N matrix.  An independent route to :func:`karman_street_spectrum` (it contains the Bloch wavenumbers k = 2πj/(Na)).
    Rigid translation gives a defective zero eigenvalue (noise ≈ 1e-8·Γ/a²).  ``stagger`` 0.5 = Kármán's staggered street, 0 = the symmetric row.  Label: analytic."""
    M, _, _ = _street_matrix(b_over_a, 8 if fast else int(n_pairs), Gamma, a, stagger)  # fast: 8 pairs per cell
    return eig(M, right=False)


def karman_street_positions(t, b_over_a: float, eps: float = 1e-2, mode: str = "unstable", n_cells: int = 8, Gamma: float = 1.0, a: float = 1.0,
                            offset: float = 0.5) -> dict:
    """Positions of the vortices of a slightly perturbed street after time t, from the LINEAR evolution of the k = π/a (alternate-vortex) eigenmode.

    Book: §9.8 (street picture, Fig. 9.17); ours.  The 4×4 system of :func:`karman_street_spectrum` at k = π/a is advanced with the matrix exponential
    from the real initial perturbation ζ_n = eps·a·(−1)ⁿ·P_r (P_r from the eigenvector of the "unstable" (largest Re σ) or "stable" (smallest) mode,
    made real by adding its image under the reality map (P, Q) → (Q̄, P̄)).  For b/a = 0.28055 (neutral) both choices oscillate.
    Parameters: t [s] scalar; b_over_a [–]; eps [–] amplitude in units of a; mode "unstable" | "stable"; n_cells number of vortices drawn per row; Gamma [m²/s]; a [m];
    offset [–] (0.5 staggered).
    Returns dict(zA, zB [complex arrays of length n_cells, m] positions of the upper and lower rows, growth [1/s] = Re σ of the chosen mode).
    Validity: linear regime only (|ζ| ≪ b, a).  Label: qualitative (linear)."""
    if mode not in ("unstable", "stable"):
        raise ValueError('mode must be "unstable" or "stable"')
    M = _street_bloch_matrix(b_over_a, offset, np.pi / a, Gamma, a, 2000)
    w, V = np.linalg.eig(M)
    i = int(np.argmax(w.real)) if mode == "unstable" else int(np.argmin(w.real))
    v = V[:, i]
    swap = lambda x: np.array([np.conj(x[1]), np.conj(x[0]), np.conj(x[3]), np.conj(x[2])])  # noqa: E731  reality map (P, Q) → (Q̄, P̄)
    x0 = v + swap(v)
    if np.linalg.norm(x0) < 1e-8 * np.linalg.norm(v):
        x0 = 1j * v + swap(1j * v)
    x0 = x0 / np.max(np.abs(x0[[0, 2]]))
    xt = expm(M * float(t)) @ (eps * a * x0)
    n = np.arange(int(n_cells))
    z_up = n * a + 0.5j * b_over_a * a
    z_lo = (n + offset) * a - 0.5j * b_over_a * a
    sign = (-1.0) ** n
    return dict(zA=z_up + xt[0] * sign, zB=z_lo + xt[2] * sign, growth=float(w[i].real))


# ----------------------------------------------------------------------------------------------------------------------
# Strouhal
# ----------------------------------------------------------------------------------------------------------------------


def shedding_frequency(U, d, St: float = 0.2) -> dict:
    """Vortex-shedding frequency f = St·U/d  [Hz] and its angular counterpart ω = 2πf  [rad/s].

    Book: §9.8, Eq. (4.102) St = Ωd/U∞ with St ≈ 0.2 over a wide Re range.  The quoted 0.2 belongs to the CYCLIC frequency f (St = fd/U): the book's Ω in
    (4.102) is the frequency in Hz.  TRAP: if Ω were the ANGULAR frequency 2πf, the same street would have Ωd/U = 2π·0.2 ≈ 1.26.
    Parameters: U [m/s]; d [m] cylinder diameter; St [–] in the cyclic form.
    Returns dict(f [Hz], omega_rad [rad/s]) — e.g. U = 10 m/s, d = 2 mm: f = 1000 Hz, ω = 6283 rad/s.  Label: analytic."""
    f = St * _F(U) / _F(d)
    return dict(f=_S(f), omega_rad=_S(2.0 * np.pi * f))


def shedding_angular_frequency(U, d, St: float = 0.2):
    """Ω = 2πf = 2π·St·U/d  [rad/s] with the cyclic-frequency St = fd/U (the Ω-vs-f trap of :func:`shedding_frequency`)."""
    return shedding_frequency(U, d, St)["omega_rad"]


def strouhal_of_re(Re):
    """Strouhal number of a circular cylinder vs Re (Roshko's plateau).  ≈ 0.21 for 300 ≤ Re < 3×10⁵ (the subcritical range; the book quotes ≈ 0.2);
    below Re = 300 (where St still varies with Re) and above the drag crisis no fit with a primary-source citation is coded: NaN and a warning.

    Book: §9.8, Eq. (4.102) recalled.  The plateau value is from Roshko (1954/1961) via secondary sources (qualitative).  Re [–] (array or float).  Label: qualitative."""
    R = _F(Re)
    ok = (R >= 300.0) & (R < 3e5)
    if not np.all(ok):
        warnings.warn("strouhal_of_re: no cited fit outside 300 <= Re < 3e5; returning NaN there", stacklevel=2)
    return _S(np.where(ok, 0.21, np.nan))


# ----------------------------------------------------------------------------------------------------------------------
# separated-pressure model (ours)
# ----------------------------------------------------------------------------------------------------------------------


def cp_ideal_cylinder(phi_deg):
    """Ideal-flow surface pressure coefficient of a circular cylinder, C_p = 1 − 4 sin²φ  (φ from the forward stagnation point).  Book: ch06 (6.x), quoted in §9.8."""
    return _S(1.0 - 4.0 * np.sin(np.deg2rad(_F(phi_deg))) ** 2)


def separated_cp(phi_deg, phi_sep_deg: float = 82.0, cp_base: float | None = None):
    """OUR separated-flow model of the cylinder surface pressure: ideal C_p = 1 − 4 sin²φ up to the separation angle, a constant base pressure behind.

    Book: §9.7–9.8 (Fig. 9.20 shows measured C_p with a plateau in the wake; the book gives no formula).
    Parameters: phi_deg [deg] from the forward stagnation point; phi_sep_deg [deg] (≈82 laminar, ≈125 turbulent layer — rounded experimental values);
    cp_base [–] wake pressure (``None``: the ideal value at separation, a crude choice; measured subcritical ≈ −1.2).
    Returns C_p(φ) [–].  Label: qualitative (a model, not a solution of the flow)."""
    phi = _F(phi_deg)
    cpb = 1.0 - 4.0 * np.sin(np.deg2rad(phi_sep_deg)) ** 2 if cp_base is None else cp_base
    return _S(np.where(phi <= phi_sep_deg, 1.0 - 4.0 * np.sin(np.deg2rad(phi)) ** 2, cpb))


def pressure_drag_from_cp(phi_deg, cp=None, n: int = 200, breakpoints=()) -> float:
    """Pressure-drag coefficient C_D,p = ½∮C_p cos φ dφ = ∫₀^π C_p cos φ dφ of a circular cylinder (φ from the forward stagnation point; friction excluded).

    Derivation: F_x per unit span = ∮ p cos φ a dφ over 0…2π, C_D = F/(½ρU²·2a) ⇒ C_D = ½∮C_p cos φ dφ (symmetric top/bottom ⇒ ∫₀^π).
    Forms: (a) ``phi_deg`` [deg] and ``cp`` [–] SAMPLES: if the samples span 0…360° the periodic trapezoid rule of ½∮ is used (endpoint 360° dropped
    if repeated); if they span 0…180° the trapezoid of ∫₀^π (second order — a jump in C_p at φ_s costs O(Δφ)); (b) ``phi_deg`` a CALLABLE cp(φ [rad]) with
    ``cp`` omitted: Gauss–Legendre with ``n`` points per sub-interval split at ``breakpoints`` [rad] (put one at a separation angle).
    Label: analytic (quadrature)."""
    if callable(phi_deg) and cp is None:
        cp_func = phi_deg
        xg, wg = np.polynomial.legendre.leggauss(int(n))
        edges = [0.0, *sorted(breakpoints), np.pi]
        tot = 0.0
        for lo, hi in zip(edges[:-1], edges[1:]):
            ph = 0.5 * (hi - lo) * xg + 0.5 * (hi + lo)
            tot += 0.5 * (hi - lo) * np.sum(wg * _F(cp_func(ph)) * np.cos(ph))
        return float(tot)
    ph, c = np.deg2rad(_F(phi_deg)), _F(cp)
    order = np.argsort(ph)
    ph, c = ph[order], c[order]
    if ph[-1] - ph[0] > np.pi + 1e-9:  # full circle: periodic trapezoid of ½∮
        if abs(ph[-1] - ph[0] - 2.0 * np.pi) < 1e-9:
            ph, c = ph[:-1], c[:-1]
        ph_c, c_c = np.append(ph, ph[0] + 2.0 * np.pi), np.append(c, c[0])
        return float(0.5 * np.trapezoid(c_c * np.cos(ph_c), ph_c))
    return float(np.trapezoid(c * np.cos(ph), ph))


def separated_pressure_drag(phi_sep_deg: float = 82.0, cp_base: float | None = -1.2, n: int = 64) -> float:
    """Pressure drag coefficient of the separated-flow model :func:`separated_cp`, by Gauss–Legendre; closed form
    C_D,p = sin φₛ (1 − (4/3) sin²φₛ − C_b).

    Book: §9.7 (form drag from a low, nearly uniform wake pressure).  φₛ → 180° with C_b → 1 gives the ideal flow, C_D → 0
    (d'Alembert).  Our model: QUALITATIVE.  Parameters: phi_sep_deg [deg]; cp_base [–] wake (base) pressure coefficient, default −1.2 (ILLUSTRATIVE: a typical
    measured subcritical magnitude, not derived; ``None`` → the crude rule C_b = ideal value at φₛ, which gives an unphysical 2.59 at 82°); n GL points.
    Returns C_D,pressure [–] (no skin friction) — expect (82°, −1.2) → 0.8838, (125°, −0.6) → 0.5778, (90°, ideal C_b = −3) → 2.667.
    NOTE: at a FIXED base pressure the drag is not monotone in φₛ (C_b = −1.2: D(82°) = 0.884, D(90°) = 0.867, D(125°) = 1.07); a later separation
    lowers the drag only through the accompanying rise of the base pressure (e.g. (125°, −0.6) → 0.578).
    Label: qualitative; analytic (closed-form parity, d'Alembert limit)."""
    ps = np.deg2rad(phi_sep_deg)
    cpb = 1.0 - 4.0 * np.sin(ps) ** 2 if cp_base is None else float(cp_base)
    return pressure_drag_from_cp(lambda ph: np.where(ph <= ps, 1.0 - 4.0 * np.sin(ph) ** 2, cpb), n=n, breakpoints=(ps,))


def _separated_drag_closed(phi_sep_deg: float, cp_base: float | None = None) -> float:
    ps = np.deg2rad(phi_sep_deg)
    s = np.sin(ps)
    cpb = 1.0 - 4.0 * s ** 2 if cp_base is None else cp_base
    return float(s * (1.0 - 4.0 / 3.0 * s ** 2 - cpb))  # C_D,p = sin φs (1 − (4/3) sin²φs − C_b)


# illustrative model parameters for the two sides of the drag crisis (ours; typical measured base pressures)
_SUBCRITICAL = dict(phi_sep_deg=82.0, cp_base=-1.2)
_SUPERCRITICAL = dict(phi_sep_deg=125.0, cp_base=-0.6)


def drag_crisis_state(Re, rough: bool = False, Re_cr: float = 3e5) -> dict:
    """Illustrative separation angle, base pressure and model pressure drag of a cylinder across the drag crisis.

    Book: §9.8 (drag crisis at Re_cr ≈ 3×10⁵: the layer turns turbulent before separating, separation moves from ≈ 82° to ≈ 125° and C_D falls by a factor ≈ 3–4;
    surface roughness triggers the transition earlier).  Model (OURS, qualitative, illustrative): φ_s = 82° below Re_cr, 125° from 2·Re_cr up, blended across the critical
    band by a smoothstep in log₂(Re/Re_cr); base pressures C_b = −1.2 (subcritical) and −0.6 (supercritical) — typical magnitudes chosen to show the crisis, not derived;
    C_D,p from :func:`separated_pressure_drag`'s closed form.  ``rough``: the transition is tripped at Re_cr/3 (no blending band: 125° for Re ≥ Re_cr/3).
    Parameters: Re [–]; rough [bool]; Re_cr [–] smooth-cylinder critical Reynolds number.
    NOTE (sphere): the sphere's critical Reynolds number is used as an argument (default 3e5, the design contract value); the book's text says ≈ 5×10⁵ (G9) — pass ``Re_cr=5e5`` to follow it.
    Returns dict(label [ASCII, :func:`cylinder_flow_regime` with the same critical band], phi_sep_deg [deg], St (None here: see :func:`cylinder_state`), cb [–],
    cd_model [–, pressure part], blend [0 = subcritical … 1 = supercritical], qualitative = True).  Label: qualitative."""
    R = float(Re)
    if rough:
        rc = Re_cr / 3.0
        w = 1.0 if R >= rc else 0.0
        reg = cylinder_flow_regime(R, thresholds=dict(critical=rc, supercritical=rc))
    else:
        t = min(max(np.log2(max(R, 1e-300) / Re_cr), 0.0), 1.0)
        w = float(t * t * (3.0 - 2.0 * t))
        reg = cylinder_flow_regime(R, thresholds=dict(critical=Re_cr, supercritical=2.0 * Re_cr))
    phi = _SUBCRITICAL["phi_sep_deg"] + w * (_SUPERCRITICAL["phi_sep_deg"] - _SUBCRITICAL["phi_sep_deg"])
    cb = _SUBCRITICAL["cp_base"] + w * (_SUPERCRITICAL["cp_base"] - _SUBCRITICAL["cp_base"])
    return dict(label=reg["label"], phi_sep_deg=float(phi), St=reg["St"], cb=float(cb), cp_base=float(cb), cd_model=_separated_drag_closed(phi, cb),
                blend=w, qualitative=True)


def cylinder_state(Re, rough: bool = False, Re_cr: float = 3e5) -> dict:
    """Illustrative state of a circular cylinder at Re: regime label, Strouhal number, separation angle, base pressure and model pressure drag.

    Book: §9.8 (regimes, St ≈ 0.2, drag crisis, roughness).  Combines :func:`cylinder_flow_regime` (label, St) with :func:`drag_crisis_state`.  The separated
    model applies to high-Re separated flow only, so ``phi_sep_deg``, ``cb`` and ``cd_model`` are ``None`` for Re < 3000 (``turbulent_wake`` threshold).
    Parameters: Re [–]; rough [bool] roughness lowers the critical Re to Re_cr/3; Re_cr [–].
    Returns dict(label, phi_sep_deg [deg], St [–] (0.2 for 40 ≤ Re < 3000, else None), cb [–], cd_model [–], qualitative = True).  QUALITATIVE (illustrative numbers).
    Validation: V7 regime boundaries; model values from the closed form.  Label: qualitative."""
    d = drag_crisis_state(Re, rough=rough, Re_cr=Re_cr)
    R = float(Re)
    st = cylinder_flow_regime(R)["St"]
    if R < CYLINDER_THRESHOLDS["turbulent_wake"]:
        return dict(label=d["label"], phi_sep_deg=None, St=st, cb=None, cp_base=None, cd_model=None, blend=d["blend"], qualitative=True)
    return dict(d, St=st)


def drag_crisis_pair(Re_cr: float = 3e5) -> dict:
    """The model before/after the drag crisis: dict(subcritical, supercritical, ratio) each with phi_sep_deg, cp_base, cd_model (ratio = supercritical/subcritical).
    Book: §9.8 (separation 82° → 125°, C_D falls sharply, "a factor ≈ 3–4" in the measured drag).  Model values are OURS (qualitative, illustrative): the
    pressure-only model gives a WEAKER drop (ratio ≈ 0.65), because it has no skin friction or wake-width physics."""
    sub = dict(_SUBCRITICAL, cd_model=_separated_drag_closed(**_SUBCRITICAL))
    sup = dict(_SUPERCRITICAL, cd_model=_separated_drag_closed(**_SUPERCRITICAL))
    return dict(subcritical=sub, supercritical=sup, ratio=sup["cd_model"] / sub["cd_model"], Re_cr=Re_cr)


_CD_ANCHORS = np.array([[1.0, 8.0 * np.pi / 2.002], [10.0, 2.8], [40.0, 1.6], [1e3, 1.0], [1e4, 1.1], [2e5, 1.2], [3e5, 1.15], [4e5, 0.6], [6e5, 0.35],
                        [2e6, 0.5], [1e7, 0.7]])


def cylinder_cd_schematic(Re):
    """SCHEMATIC C_D(Re) of a smooth cylinder (log–log interpolation through rounded anchor points; below Re = 1 Lamb's asymptote 8π/[Re(2.002 − ln Re)]; the Re = 1 anchor is Lamb's value 8π/2.002 = 12.55, so the curve is continuous).

    Book: Fig. 9.21 (the book's data are private and no dataset was fetched — analysis §8).  QUALITATIVE: shows ≈60 at Re = 0.1, ≈1 in the middle,
    a dip near the critical Reynolds number and a slow recovery; do not use for numbers.  Re [–]; returns C_D [–].  Label: qualitative."""
    R = _F(Re)
    lamb = 8.0 * np.pi / (R * (2.002 - np.log(np.maximum(R, 1e-12))))
    lx, ly = np.log(_CD_ANCHORS[:, 0]), np.log(_CD_ANCHORS[:, 1])
    mid = np.exp(np.interp(np.log(np.maximum(R, 1e-12)), lx, ly))
    return _S(np.where(R < 1.0, lamb, mid))


# ----------------------------------------------------------------------------------------------------------------------
# sports balls
# ----------------------------------------------------------------------------------------------------------------------


def ball_swing_deflection(F_over_W: float, distance: float, U: float, g: float = 9.81):
    """Sideways deflection y = ½ a t² of a swinging ball under a constant lateral force F = (F/W)·mg  [m]; t = distance/U, a = (F/W)g.

    Book: §9.9 (cricket-ball swing: seam trips one side; constant lateral force ⇒ parabolic path).  Parameters: F_over_W [–] side force / weight;
    distance [m] travelled; U [m/s]; g [m/s²].  Our one-line kinematics; the parameters are the caller's.  Label: analytic."""
    return _S(0.5 * F_over_W * g * (_F(distance) / U) ** 2)


def magnus_sign(Re_slow, Re_fast, Re_cr: float = 3e5) -> str:
    """Sign of the Magnus effect of a spinning ball from which side is past the drag crisis — a truth table (qualitative).

    Convention (book §9.9): the Reynolds number of each side is built on the fluid velocity RELATIVE TO THE SURFACE.  ``Re_fast`` belongs to the side with the
    LARGER relative velocity — the side whose surface moves AGAINST the oncoming flow (clockwise spin: the lower side in the book's figure) — and ``Re_slow`` to
    the side whose surface moves with the flow (Re_slow ≤ Re_fast).  Book: if only the larger-relative-speed side is past the crisis (Re_slow < Re_cr ≤ Re_fast) its
    separation is delayed, giving lower pressure on that side and a force opposite to the ideal-flow Magnus force: NEGATIVE Magnus effect.  If both sides are past
    the crisis (Re_cr ≤ Re_slow) the separation point moves upstream with Re on the faster side, its pressure is higher, and the ordinary POSITIVE effect results.
    (The book's sentence prints 'Re < Re_cr' twice — slip R10; the second is Re > Re_cr.)  Both sides below Re_cr is NOT discussed by the book (both layers
    laminar) and the answer is not asserted here.
    Returns "+", "−" (Unicode minus) or "none" (undetermined: Re_slow == Re_fast, i.e. no spin; the convention Re_slow ≤ Re_fast violated; or both sides
    subcritical, which the book does not treat).  Label: qualitative."""
    if Re_slow >= Re_fast:
        return "none"
    if Re_slow < Re_cr <= Re_fast:
        return "−"
    if Re_slow >= Re_cr:
        return "+"
    return "none"  # both subcritical: not supported by the book's argument


__all__ = [n for n in dir() if not n.startswith("_") and n not in ("annotations", "warnings", "eig", "expm")]
