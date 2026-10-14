"""Chapter 16 — Introduction to Biofluid Mechanics (Kundu, Cohen & Dowling 5e, §16.1–§16.4, Eqs. (16.1)–(16.222)).

The whole book applied to blood and sap: Poiseuille flow and its resistance, a cell-free plasma layer, distensible
tubes, the Windkessel, pulse waves (Moens–Korteweg), Womersley's oscillatory pipe flow, viscous waves in an elastic
tube, Murray's law, wave reflection at junctions, Dean flow in a curved tube, collapsible tubes and flow limitation,
Casson (yield-stress) pipe flow, and water transport in plants.

**This module re-exports** every public name of the two new core modules

* ``fluidpy.core.womersley`` — the Womersley number (16.12) and rigid-tube oscillatory flow (16.69)–(16.100), complex
  Bessel functions (series and overflow-free ratios);
* ``fluidpy.core.rheology`` — Newtonian, Casson (16.5)–(16.7), (16.207)–(16.220), power-law and Herschel–Bulkley laws
  and their pipe fluxes;

and wraps the earlier chapters' tools it re-uses: ch08 ``core.laminar`` (Poiseuille (8.6)–(8.8), the Stokes layer),
ch15 ``core.acoustics`` (d'Alembert, p′ = ρc u′), ch10 ``core.maccormack`` (the nonlinear 1-D tube model), ch01
``capillary_rise``.  Everything an explainer parity row needs is reachable as ``ch16.<name>`` and is scalar-callable.

Conventions (analysis ch16 §9 — traps T1–T27; ``conventions_table()``)
* x along the vessel, r radial (cylindrical (r, θ, x)); Δp = p₁ − p₂ > 0 drives flow in +x; **G = −dp/dx**.
* Stresses in §§16.3-S/FL/CA are positive magnitudes τ = −μ du/dr (ch08 code uses signed τ = μ du/dR) — trap T5.
* SI inside every function; mmHg, CGS (poise, dyn/cm²), PRU and bar only through the converters (trap T3).
* **Two time conventions** (trap T1): rigid-tube Womersley uses e^{+iωt} (16.73); the viscous elastic-tube waves use
  e^{i(kx − ωt)} (16.116).  Every docstring says which.
* Overloaded symbols get distinct code names (trap T2): Womersley α = ``alpha``, Shapiro's area ratio A/A₀ =
  ``area_ratio``, toroidal angle = ``alpha_t``; compliance = ``K`` (Windkessel functions only), metabolic cost per
  volume = ``b_metabolic``, Dean's K = ``K_dean``; peripheral resistance ``R`` (Windkessel), reflection coefficient
  ``refl``, radius of curvature ``R_curv``; Poisson ratio = ``poisson`` (ν is kinematic viscosity ``nu``).
* Hematocrit as a fraction (trap T13).  Inputs used in examples are ours (:func:`illustrative_inputs`), never the
  book's.

**Printed slips coded corrected** (analysis §9; each changes the physics or a number) — ``printed_*`` options exist for
the tests that must fail: (16.30) ``two_layer_flow_rates(printed_16_30=True)``; (16.44)–(16.45) the a₀⁴ (pint check);
(16.54) ``tube1d_run(printed_16_54=True)``; (16.76) right side −G/μ (``womersley_ode_residual``); (16.92) wall shear
*stress* with F₂ (``womersley_wall_shear``); (16.96) ``womersley_F2_asymptote(printed_16_96=True)``; (16.140)/(16.143)
``murray_radius(printed_16_140=True)``; (16.201) ``collapsible_differentials_sympy(printed_16_201=True)``; (16.174)
advection-term sign (S34, found in verification: eliminating p₁ from (16.167)–(16.168) gives +(1/r)(…)∇²ψ on the
left) ``dean_series_sympy(printed_16_174=True)`` — only the K⁴ flux coefficient changes (``DEAN_C4_DERIVED``); the
minimum-pressure-drop example (S26), the cP unit chain (S27) and the xylem viscosity (S28) are number slips handled in
the converters and the notebook.  ``book_slips()`` lists all 34 rows.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from fractions import Fraction
from functools import lru_cache
from typing import Callable

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

from .core import acoustics as AC
from .core import laminar as LAM
from .core import womersley as WOM
from .core.rheology import *  # noqa: F401,F403
from .core.womersley import *  # noqa: F401,F403
from .core._util import as_scalar_if_0d, require_positive
from .core.acoustics import (acoustic_impedance, dalembert, pulse, velocity_from_pressure,  # noqa: F401
                             wave_equation_step, wave_speeds)
from .core.laminar import (pipe_flow_rate, pipe_poiseuille, pipe_wall_stress, stokes_layer,  # noqa: F401
                           stokes_second_problem)
from .core.maccormack import maccormack_step  # noqa: F401
from .core.similarity import reynolds_number, strouhal_number  # noqa: F401
from .core.thermo import G0, R_U
from .core.units import Q_, check_dimensions  # noqa: F401
from .ch01_introduction import capillary_rise, newton_shear_stress  # noqa: F401

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d

# ======================================================================================================================
# Units and constants (I1; trap T3, slip S27)
# ======================================================================================================================
MMHG: float = 133.322387415          #: 1 mmHg in Pa (conventional, exact by definition: 13.5951 g/cm³ × 9.80665 m/s² × 1 mm)
PRU: float = MMHG / 1.0e-6           #: peripheral resistance unit, 1 mmHg per cm³/s, in Pa s/m³
BAR: float = 1.0e5                   #: 1 bar in Pa (exact)
R_MOLAR: float = R_U / 1000.0        #: universal gas constant per mole [J/(mol K)] (CODATA 2019 exact, via core.thermo)
R_BOOK_CM3_BAR: float = 83.141       #: the book's R in cm³ bar/(mol K) for (16.221) (= R_MOLAR × 10, consistency check)

#: Blood and plasma properties used as defaults (cited public values, not the book's):
#: ρ = 1060 kg/m³ and whole-blood viscosity 3–4 mPa s at 37 °C (Wikipedia "Hemorheology", fetched 2026-10-14; we take
#: the high-shear value 3.5 mPa s); plasma viscosity ≈ 1.2 mPa s (same source: 1.10–1.30 mPa s); body temperature
#: 310.15 K.
BLOOD: dict = dict(rho=1060.0, mu=3.5e-3, mu_p=1.2e-3, T=310.15)


def _check_hematocrit(H):
    """Validate a hematocrit given as a fraction (trap T13): raises if H > 1 (a percent passed by mistake) or H ≤ 0.

    Book: §16.2 (hematocrit defined as a volume fraction; the prose quotes percents).  Returns H as float/array.
    """
    H = _F(H)
    if np.any(H > 1.0):
        raise ValueError("hematocrit must be a FRACTION (e.g. 0.45), not a percent (trap T13)")
    if np.any(H <= 0.0):
        raise ValueError("hematocrit must be > 0")
    return _S(H)


def mmHg_to_Pa(x):
    """mmHg → Pa (1 mmHg = 133.322387415 Pa).  Book: §16.2 (the conversion quoted on p. 789). Label: analytic."""
    return _S(_F(x) * MMHG)


def Pa_to_mmHg(x):
    """Pa → mmHg.  Book: §16.2. Label: analytic."""
    return _S(_F(x) / MMHG)


def cP_to_Pa_s(x):
    """centipoise → Pa s: 1 cP = 10⁻² P = 10⁻² dyn s/cm² = 10⁻³ Pa s.

    Book: §16.2 units passage, printed "1 cP = 0.01 P = 1 dyn s/cm² = 0.1 N s/m²".
    DEVIATION: slip S27 — the last two members are 100× too large; 1 cP = 0.01 dyn s/cm² = 1e-3 Pa s (pint-verified).
    Label: analytic.
    """
    return _S(_F(x) * 1.0e-3)


def poise_to_Pa_s(x):
    """poise → Pa s (1 P = 1 dyn s/cm² = 0.1 Pa s).  Book: §16.2. Label: analytic."""
    return _S(_F(x) * 0.1)


def dyn_cm2_to_Pa(x):
    """dyn/cm² → Pa (× 0.1).  Book: §16.2, after (16.5). Label: analytic."""
    return _S(_F(x) * 0.1)


def bar_to_Pa(x):
    """bar → Pa (× 1e5).  Book: §16.4, after (16.221). Label: analytic."""
    return _S(_F(x) * BAR)


def mL_to_m3(x):
    """millilitre → m³ (× 1e-6).  Book: §16.2 (volumes quoted in mL). Label: analytic."""
    return _S(_F(x) * 1.0e-6)


# ======================================================================================================================
# §16.2 circulation bookkeeping (I1, I2)
# ======================================================================================================================
def pulse_pressure(pS, pD):
    """Pulse pressure p_S − p_D (systolic minus diastolic).

    Book: §16.2 "The Heart as a Pump" (definition) and after (16.11).
    Parameters
    ----------
    pS, pD : systolic and diastolic pressures [Pa] (any one unit for both).
    Returns
    -------
    pulse pressure [same unit].
    Assumptions: none.  Validation: V1. Label: analytic.
    """
    return _S(_F(pS) - _F(pD))


def stroke_volume(EDV, ESV):
    """Stroke volume SV = EDV − ESV.

    Book: §16.2, Eq. (16.1).
    Parameters
    ----------
    EDV, ESV : end-diastolic and end-systolic ventricular volumes [m³] (or any one unit).
    Returns
    -------
    SV [same unit].
    Assumptions: none.  Validation: V1; read off :func:`pv_loop_schematic`. Label: analytic.
    """
    return _S(_F(EDV) - _F(ESV))  # Eq. (16.1)


def ejection_fraction(SV, EDV, *, percent: bool = False):
    """Ejection fraction EF = SV/EDV (the book's (SV/EDV) × 100 % as a fraction; ``percent=True`` multiplies by 100).

    Book: §16.2, Eq. (16.2).
    Parameters
    ----------
    SV, EDV : [m³] (same unit);  percent : return a percent instead of a fraction.
    Returns
    -------
    EF [–] (or [%]).
    Assumptions: none.  Validation: V1. Label: analytic.
    """
    ef = _F(SV) / _F(EDV)  # Eq. (16.2)
    return _S(100.0 * ef if percent else ef)


def cardiac_output(SV, HR):
    """Cardiac output CO = SV × HR, returned in m³/s (the book works in mL/min).

    Book: §16.2, Eq. (16.3).  (The book writes one value in "mL/m" — slip S32, mL/min.)
    Parameters
    ----------
    SV : stroke volume [m³];  HR : heart rate [beats/min].
    Returns
    -------
    CO [m³/s] (multiply by 6e4 for L/min — ×60 s/min × 1000 L/m³ — or by 6e7 for mL/min).
    Assumptions: steady beat-to-beat.  Validation: V1; V4 equals the time mean of the ejection flow. Label: analytic.
    """
    return _S(_F(SV) * _F(HR) / 60.0)  # Eq. (16.3), per second


def cardiac_index(CO, BSA):
    """Cardiac index CI = CO/BSA.  Book: §16.2, Eq. (16.4).

    Parameters: CO [m³/s];  BSA body surface area [m²].  Returns CI [m³/(s m²)].  Label: analytic.
    """
    return _S(_F(CO) / _F(BSA))  # Eq. (16.4)


def peripheral_resistance(pA_mean, Q, pV_mean=0.0):
    """Total peripheral resistance R = Δp̄/Q with Δp̄ = p̄_A − p̄_V.

    Book: §16.2, Eqs. (16.9)–(16.10); with p̄_V = 0, p̄_A = QR.  Divide by :data:`PRU` for mmHg per cm³/s.
    Parameters
    ----------
    pA_mean : time-mean arterial pressure [Pa];  Q : time-mean flow (cardiac output) [m³/s];
    pV_mean : time-mean venous (right-atrial) pressure [Pa].
    Returns
    -------
    R [Pa s/m³].
    Assumptions: time averages over a cardiac cycle.  Validation: V1. Label: analytic.
    """
    return _S((_F(pA_mean) - _F(pV_mean)) / _F(Q))  # Eqs. (16.9), (16.10)


def mean_arterial_pressure(pS, pD):
    """The one-third rule p̄_A ≈ ⅓p_S + ⅔p_D = p_D + ⅓(p_S − p_D) — an *estimate* (trap T15).

    Book: §16.2, Eq. (16.11).
    Parameters
    ----------
    pS, pD : systolic, diastolic pressures [Pa].
    Returns
    -------
    p̄_A estimate [Pa].
    Assumptions: a typical waveform.  Validation: V1; compared with :func:`time_mean` of a Windkessel waveform.
    Label: analytic.
    """
    return _S(_F(pD) + (_F(pS) - _F(pD)) / 3.0)  # Eq. (16.11)


def time_mean(p, t):
    """Time mean of a sampled signal over its span, (1/(t_N − t_0))∫p dt (trapezoid; ours).

    Book: §16.2 (the time averages of (16.9)–(16.11)).
    Parameters
    ----------
    p : samples (any unit);  t : times [s] (increasing, one period for a periodic mean).
    Returns
    -------
    mean [unit of p].
    Assumptions: samples resolve the signal.  Validation: V1 (a constant; a sine over a period). Label: analytic.
    """
    p, t = _F(p), _F(t)
    return float(np.trapezoid(p, t) / (t[-1] - t[0]))


def pv_loop_schematic(EDV, ESV, *, p_min: float = 0.6e3, p_ED: float = 1.2e3, p_open: float = 10.5e3,
                      p_peak: float = 15.5e3, p_ES: float = 13.0e3, n: int = 200) -> dict:
    """A synthetic left-ventricular pressure–volume loop A → B → C → D → A (ours; illustrative pressures).

    Book: §16.2 "Work Done by the Left Ventricle", Fig. 16.6 (described, not copied): A→B filling (mitral open,
    passive rise to the end-diastolic pressure), B→C isovolumic contraction at EDV, C→D ejection (aortic open, volume
    falls to ESV), D→A isovolumic relaxation at ESV.  Traversed **counter-clockwise** in the (V, p) plane (trap T11).
    Leg shapes (ours, design C.1): filling p = p_min + (p_ED − p_min)s², s = (V − ESV)/(EDV − ESV); ejection the
    parabola through (EDV, p_open), (½(EDV + ESV), p_peak), (ESV, p_ES); the isovolumic legs straight.
    Parameters
    ----------
    EDV, ESV : volumes [m³];  p_min : pressure at A, start of filling [Pa];  p_ED : end-diastolic pressure (B) [Pa];
    p_open : pressure at aortic-valve opening (C) [Pa];  p_peak : pressure at mid-ejection [Pa];  p_ES : end-systolic
    pressure (D) [Pa];  n : points per curved leg.
    Returns
    -------
    dict: V [m³], p [Pa] (closed polygon, last point = first), corners {A, B, C, D: (V, p)}.
    Assumptions: quasi-static chamber, uniform pressure.
    Validation: V1 :func:`pv_loop_work` equals :func:`pv_loop_schematic_work` (exact for these legs at large n).
    Label: analytic.
    """
    n = int(n)
    s = np.linspace(0.0, 1.0, n, endpoint=False)
    dV = EDV - ESV
    V1 = ESV + dV * s
    p1 = p_min + (p_ED - p_min) * s ** 2                                     # A → B filling
    m = 4
    s2 = np.linspace(0.0, 1.0, m, endpoint=False)
    V2 = np.full_like(s2, EDV)
    p2 = p_ED + (p_open - p_ED) * s2                                          # B → C isovolumic contraction
    q = np.linspace(0.0, 1.0, n, endpoint=False)                              # 0 at C, 1 at D
    V3 = EDV - dV * q
    p3 = p_open * (1 - q) * (1 - 2 * q) + 4 * p_peak * q * (1 - q) + p_ES * q * (2 * q - 1)  # C → D (Lagrange)
    V4 = np.full_like(s2, ESV)
    p4 = p_ES + (p_min - p_ES) * s2                                           # D → A isovolumic relaxation
    V = np.concatenate([V1, V2, V3, V4, [V1[0]]])
    p = np.concatenate([p1, p2, p3, p4, [p1[0]]])
    corners = dict(A=(float(ESV), float(p_min)), B=(float(EDV), float(p_ED)), C=(float(EDV), float(p_open)),
                   D=(float(ESV), float(p_ES)))
    return dict(V=V, p=p, corners=corners)


def pv_loop_schematic_work(EDV, ESV, *, p_min: float = 0.6e3, p_ED: float = 1.2e3, p_open: float = 10.5e3,
                           p_peak: float = 15.5e3, p_ES: float = 13.0e3) -> float:
    """Exact work per beat of :func:`pv_loop_schematic`: W = ΔV[⅙(p_ES + 4p_peak + p_open) − p_min − ⅓(p_ED − p_min)].

    Book: §16.2, Fig. 16.6 (W = −∮p dV, ours): Simpson's rule is exact for the parabolic ejection leg, and the filling
    leg integrates to ΔV[p_min + ⅓(p_ED − p_min)].
    Parameters as :func:`pv_loop_schematic`.  Returns W [J].  Label: analytic.
    """
    dV = float(EDV) - float(ESV)
    return dV * ((p_ES + 4.0 * p_peak + p_open) / 6.0 - p_min - (p_ED - p_min) / 3.0)


def pv_loop_work(p, V) -> dict:
    """Work done on the blood per beat, W = −∮p dV (trapezoid along the polygon) = the signed (shoelace) area of the
    loop in the (V, p) plane.

    Book: §16.2, Fig. 16.6 ("the area within the pressure–volume diagram represents the net work done by the LV on
    blood"); the sign convention −∮p dV is ours (derivation D01): for the counter-clockwise loop ∮p dV < 0 and W > 0.
    Parameters
    ----------
    p : pressures around the loop [Pa];  V : volumes [m³] (same length; closed or open polygon).
    Returns
    -------
    dict: work [J] (= −∮p dV), signed_area [J] (shoelace in (V, p), > 0 counter-clockwise), orientation
    ("counterclockwise" / "clockwise").
    Assumptions: quasi-static, the polygon resolves the loop.
    Validation: V1 a rectangle gives Δp ΔV exactly; reversing the order flips only the sign. Label: analytic.
    """
    p, V = _F(p), _F(V)
    if p[0] == p[-1] and V[0] == V[-1]:
        p, V = p[:-1], V[:-1]
    pn, Vn = np.roll(p, -1), np.roll(V, -1)
    oint = float(np.sum(0.5 * (p + pn) * (Vn - V)))  # ∮ p dV, trapezoid on each edge (exact for a polygon)
    signed_area = float(0.5 * np.sum(V * pn - Vn * p))  # shoelace
    return dict(work=-oint, signed_area=signed_area,
                orientation="counterclockwise" if signed_area > 0 else "clockwise")


def cardiac_power(W, HR):
    """Mean mechanical power of the ventricle, W × HR/60 [W].  Book: §16.2 (work per beat × beats per second; ours).

    Parameters: W work per beat [J];  HR [beats/min].  Returns power [W].  Label: analytic.
    """
    return _S(_F(W) * _F(HR) / 60.0)


# ======================================================================================================================
# §16.2 blood and vessels (I3)
# ======================================================================================================================
def equal_volume_sphere_area(V):
    """Surface area of the sphere with volume V: (36πV²)^{1/3} (ours; compare a red cell's area with it).

    Book: §16.2 (RBC dimensions described; the comparison is ours).
    Parameters: V [m³].  Returns area [m²].  Validation: V1 (r = 1 sphere). Label: analytic.
    """
    return _S(np.cbrt(36.0 * np.pi * _F(V) ** 2))


def mean_speed(Q, A_total):
    """Mean speed ū = Q/A_total across a vessel class (continuity; ours, Fig. 16.10's reading).

    Book: §16.2, Fig. 16.10.  Parameters: Q [m³/s];  A_total total cross-section of the class [m²].  Returns [m/s].
    Label: analytic.
    """
    return _S(_F(Q) / _F(A_total))


#: Illustrative vessel classes (ours, order of magnitude only — Wikipedia "Blood vessel" and "Womersley number",
#: fetched 2026-10-14; deliberately not the book's table): radius a [m], wall thickness h [m], mean speed u_mean [m/s].
VESSELS: dict = {
    "aorta": dict(a=11.0e-3, h=1.6e-3, u_mean=0.25),
    "large artery": dict(a=1.5e-3, h=0.45e-3, u_mean=0.12),
    "arteriole": dict(a=20.0e-6, h=12.0e-6, u_mean=5.0e-3),
    "capillary": dict(a=4.0e-6, h=0.5e-6, u_mean=5.0e-4),
    "venule": dict(a=12.0e-6, h=1.0e-6, u_mean=2.0e-3),
    "vein": dict(a=2.0e-3, h=0.2e-3, u_mean=0.05),
    "vena cava": dict(a=14.0e-3, h=1.1e-3, u_mean=0.12),
}
VESSELS_ILLUSTRATIVE = VESSELS  # alias


def vessel_scales_table(vessels: dict | None = None, *, HR: float = 66.0, rho: float = 1060.0, mu: float = 3.5e-3):
    """Reynolds number (diameter), Womersley number and wall-thickness ratio per vessel class (ours).

    Book: §16.2 "Nature of Blood Vessels" (classes described) with α = a√(ω/ν) (16.12) and Re = ρū·2a/μ (§16.3).  The
    inputs are :data:`VESSELS` (ours, cited orders of magnitude), never the book's values.
    Parameters
    ----------
    vessels : {name: dict(a [m], h [m], u_mean [m/s])} (default :data:`VESSELS`);  HR [beats/min];  rho [kg/m³];
    mu [Pa s].
    Returns
    -------
    pandas.DataFrame: vessel, a, h, u_mean, Re, alpha, h_over_a, regime.
    Assumptions: Newtonian blood (high-shear viscosity).  Validation: V1 per row. Label: analytic.
    """
    import pandas as pd

    vessels = VESSELS if vessels is None else vessels
    omega = 2.0 * np.pi * HR / 60.0
    rows = []
    for name, v in vessels.items():
        a, h, u = v["a"], v["h"], v["u_mean"]
        al = float(WOM.womersley_number(a, omega, mu / rho))
        rows.append(dict(vessel=name, a=a, h=h, u_mean=u, Re=rho * u * 2.0 * a / mu, alpha=al, h_over_a=h / a,
                         regime=WOM.womersley_regime(al)))
    return pd.DataFrame(rows)


# ======================================================================================================================
# §16.3 steady Poiseuille flow and resistance (I6, I7) — thin wrappers of ch08 core.laminar
# ======================================================================================================================
def poiseuille_velocity(r, a, dp, L, mu):
    """Hagen–Poiseuille profile u(r) = (a²Δp/4μL)(1 − r²/a²).

    Book: §16.3, Eqs. (16.13)–(16.15) (recap of ch08 (8.6)); calls ``core.laminar.pipe_poiseuille`` with G = Δp/L.
    Parameters
    ----------
    r [m];  a : radius [m];  dp : Δp = p₁ − p₂ [Pa];  L : length [m];  mu [Pa s].
    Returns
    -------
    u [m/s].
    Assumptions: steady, fully developed, rigid, Newtonian.
    Validation: V1 identical to ch08 (8.6). Label: analytic.
    """
    return LAM.pipe_poiseuille(r, a, mu=mu, G=_F(dp) / _F(L))  # Eqs. (16.13)–(16.15)


def poiseuille_flow_rate(a, dp, L, mu):
    """Poiseuille formula Q = πa⁴Δp/(8μL).

    Book: §16.3, Eq. (16.17) (ch08 ``pipe_flow_rate``); mean speed (16.18) V = Q/(πa²) = u_max/2, u_max (16.16).
    Parameters
    ----------
    a [m];  dp [Pa];  L [m];  mu [Pa s].
    Argument order (positional): a, dp, L, mu — :func:`vascular_resistance` takes (mu, L, a); call by keyword when in
    doubt.
    Returns
    -------
    Q [m³/s].
    Assumptions: as :func:`poiseuille_velocity`.
    Validation: V1 quadrature of :func:`poiseuille_velocity`; equals ch08. Label: analytic.
    """
    Q, _, _ = LAM.pipe_flow_rate(a, mu=mu, G=_F(dp) / _F(L))
    return Q  # Eq. (16.17)


def poiseuille_wall_stress(a, dp, L):
    """Wall shear stress magnitude τ = −μ du/dr|_a = aΔp/(2L) > 0.

    Book: §16.3, Eq. (16.19); the last member is printed −(a/2)Δp/L, it is +(a/2)Δp/L (slip S1, loose; trap T5:
    a positive magnitude here, signed in ch08's ``pipe_wall_stress``).
    Parameters
    ----------
    a [m];  dp [Pa];  L [m].
    Returns
    -------
    tau [Pa] (> 0 for Δp > 0).
    Assumptions: as :func:`poiseuille_velocity`.  Validation: V4 force balance; V1 = −ch08 value. Label: analytic.
    """
    return _S(-_F(LAM.pipe_wall_stress(a, G=_F(dp) / _F(L))))  # Eq. (16.19), corrected sign


def vascular_resistance(mu, L, a):
    """Poiseuille resistance R = Δp/Q = 8μL/(πa⁴) — halve the radius, multiply R by 16.

    Book: §16.3, Eq. (16.20) (the text says "fourth power of vessel diameter"; R ∝ a⁻⁴ ∝ d⁻⁴ either way).
    Parameters
    ----------
    mu [Pa s];  L [m];  a : radius [m].
    Argument order (positional): mu, L, a — the reverse of :func:`poiseuille_flow_rate`'s (a, dp, L, mu); swapping
    a and μ is silent and wrong by orders of magnitude, so call by keyword when in doubt.
    Returns
    -------
    R [Pa s/m³] (÷ :data:`PRU` for mmHg s/cm³).
    Assumptions: as :func:`poiseuille_velocity`.
    Validation: V1 Δp/Q of :func:`poiseuille_flow_rate`; R a⁴ constant. Label: analytic.
    """
    return _S(8.0 * _F(mu) * _F(L) / (np.pi * _F(a) ** 4))  # Eq. (16.20)


def network_resistance(R_list, *, arrangement: str):
    """Resistance of vessels in series (ΣR_i) or in parallel (1/Σ(1/R_i)) (ours; Ohm-like Δp = QR).

    Book: §16.3 after (16.20) (arterioles as the main resistance; the network rule is ours).
    Parameters: R_list sequence [Pa s/m³];  arrangement "series" | "parallel" (required keyword).
    Returns R_total [Pa s/m³].
    Validation: V1 N identical parallel tubes give R/N. Label: analytic.
    """
    R = _F(R_list)
    if arrangement == "series":
        return float(np.sum(R))
    if arrangement == "parallel":
        return float(1.0 / np.sum(1.0 / R))
    raise ValueError('arrangement must be "series" or "parallel"')


def series_pressure_drop(R_list, Q, p_in):
    """Pressures at the nodes of resistances in series carrying the same Q: p_{k+1} = p_k − QR_k (ours, Fig. 16.2).

    Book: §16.2, Fig. 16.2 (mean pressure falling along the circulation; the series model is ours).
    Parameters: R_list [Pa s/m³];  Q [m³/s];  p_in [Pa].  Returns array of the pressure after each element
    (len(R_list)) [Pa].
    Label: analytic.
    """
    return float(p_in) - float(Q) * np.cumsum(_F(R_list))


def entrance_length(Re, d, *, regime: str, allow_low_re: bool = False):
    """Entrance (inlet) length ℓ from the book's correlations: ℓ/d = 0.06 Re (laminar, Re > 50), ℓ/d = 0.693 Re^{1/4}
    (turbulent).

    Book: §16.3 "Effect of Developing Flow", Eq. (16.37).  Other sources give 0.0575 Re and 1.359 Re^{1/4} (Çengel, via
    Wikipedia "Entrance length (fluid dynamics)") — a different turbulent correlation, factor ≈ 2 (slip S31, loose).
    Parameters
    ----------
    Re : Reynolds number on the diameter [–];  d : diameter [m];  regime : "laminar" | "turbulent";
    allow_low_re : allow the laminar formula at Re ≤ 50 (outside its stated range; warns).
    Returns
    -------
    l [m].
    Assumptions: steady flow, uniform inlet.  Validation: V1; V5 within 5 % of 0.0575 Re (laminar). Label: analytic.
    """
    Re, d = _F(Re), _F(d)
    if regime == "laminar":
        if np.any(Re <= 50.0):
            if not allow_low_re:
                raise ValueError("(16.37) laminar form is stated for Re > 50 (pass allow_low_re=True to extrapolate)")
            warnings.warn("(16.37) used outside Re > 50", stacklevel=2)
        return _S(0.06 * Re * d)  # Eq. (16.37), laminar
    if regime == "turbulent":
        return _S(0.693 * Re ** 0.25 * d)  # Eq. (16.37), turbulent
    raise ValueError('regime must be "laminar" or "turbulent"')


# ======================================================================================================================
# §16.2 + §16.3 Fåhræus–Lindqvist: hematocrit bookkeeping and the two-layer model (I5)
# ======================================================================================================================
def cell_conservation(Qc, Qp, H_c, a, delta) -> dict:
    """Red-cell and volume bookkeeping of a core (cells) + cell-free plasma layer.

    Book: §16.2, Eq. (16.8): Q_F H_F = Q_c H_c, Q_c + Q_p = Q_F, H_T a² = H_c(a − δ)²; and the ratio printed after
    (16.34), H_c/H_F = 1 + Q_p/Q_c (derivation a-d03).
    Parameters
    ----------
    Qc, Qp : core and plasma-layer flow rates [m³/s];  H_c : core hematocrit [–];  a : radius [m];
    delta : cell-free layer thickness [m].
    Returns
    -------
    dict: Q_F, H_F (feed = discharge), H_T (tube), Hc_over_HF (= 1 + Q_p/Q_c), HT_over_Hc (= (1 − δ/a)²),
    HT_over_HD.
    Assumptions: steady, no cells in the layer.
    Validation: V1 H_T < H_D whenever the core moves faster than the mean. Label: analytic.
    """
    Qc, Qp, H_c, a, delta = (float(v) for v in (Qc, Qp, H_c, a, delta))
    _check_hematocrit(H_c)
    QF = Qc + Qp  # Eq. (16.8)
    HF = Qc * H_c / QF
    HT = H_c * (a - delta) ** 2 / a ** 2
    return dict(Q_F=QF, H_F=HF, H_T=HT, Hc_over_HF=1.0 + Qp / Qc, HT_over_Hc=(a - delta) ** 2 / a ** 2,
                HT_over_HD=HT / HF)


def two_layer_velocity(r, a, delta, dp, L, mu_p, mu_c):
    """Axial velocity of the two-region (core + cell-free plasma layer) Poiseuille flow.

    Book: §16.3 "Hagen-Poiseuille Flow and the Fahraeus-Lindqvist Effect", Eqs. (16.27) (plasma, a − δ ≤ r ≤ a) and
    (16.28) (core, 0 ≤ r ≤ a − δ), from (16.21)–(16.26).  (The signs printed in (16.21) and (16.24) are inconsistent —
    slip S2, loose; the profiles (16.27)–(16.28) are right.)
    Parameters
    ----------
    r [m];  a : radius [m];  delta : layer thickness [m] (0 ≤ δ < a);  dp : Δp [Pa];  L [m];
    mu_p : plasma viscosity [Pa s];  mu_c : core viscosity [Pa s].
    Returns
    -------
    u [m/s].
    Assumptions: steady, Newtonian layers, sharp interface, rigid tube.
    Validation: V2 sympy: the stress laws, regularity, stress and velocity continuity, no slip (16.21)–(16.26).
    Label: symbolic.
    """
    r = _F(r)
    a, delta, dp, L, mu_p, mu_c = (float(v) for v in (a, delta, dp, L, mu_p, mu_c))
    pre = a ** 2 / (4.0 * mu_p) * dp / L
    s, sc = r / a, (a - delta) / a
    up = pre * (1.0 - s ** 2)  # Eq. (16.27)
    uc = pre * (1.0 - sc ** 2 - mu_p / mu_c * s ** 2 + mu_p / mu_c * sc ** 2)  # Eq. (16.28)
    return _S(np.where(r <= a - delta, uc, up))


def two_layer_flow_rates(a, delta, dp, L, mu_p, mu_c, *, printed_16_30: bool = False):
    """Plasma-layer, core and total flow rates of the two-layer model.

    Book: §16.3, Eqs. (16.29) Q_p = (πΔp/(8μ_pL))[a² − (a − δ)²]², (16.30) Q_c, (16.31)
    Q = (πa⁴Δp/(8μ_pL))[1 − (1 − δ/a)⁴(1 − μ_p/μ_c)].
    DEVIATION: slip S3 — (16.30) is printed (πa²Δp/(4μ_pL))[a² − (1 − μ_p/(2μ_c))(a − δ)⁴/a²]; integrating (16.28)
    gives Q_c = (π(a − δ)²Δp/(4μ_pL))[a² − (1 − μ_p/(2μ_c))(a − δ)²], the only form that adds with (16.29) to (16.31).
    ``printed_16_30=True`` returns the printed Q_c (and the inconsistent sum) for the test that must fail.
    Parameters
    ----------
    a [m];  delta [m];  dp [Pa];  L [m];  mu_p, mu_c [Pa s];  printed_16_30 : use the printed (16.30).
    Returns
    -------
    dict: Qp, Qc, Q (= Qp + Qc) [m³/s].
    Assumptions: as :func:`two_layer_velocity`.
    Validation: V1 quadrature of :func:`two_layer_velocity`; Q equals (16.31). Label: analytic, symbolic.
    """
    a, delta, dp, L, mu_p, mu_c = (float(v) for v in (a, delta, dp, L, mu_p, mu_c))
    b = a - delta
    Qp = np.pi * dp / (8.0 * mu_p * L) * (a ** 2 - b ** 2) ** 2  # Eq. (16.29)
    if printed_16_30:
        Qc = np.pi * a ** 2 * dp / (4.0 * mu_p * L) * (a ** 2 - (1.0 - mu_p / (2.0 * mu_c)) * b ** 4 / a ** 2)  # printed
    else:
        Qc = np.pi * b ** 2 * dp / (4.0 * mu_p * L) * (a ** 2 - (1.0 - mu_p / (2.0 * mu_c)) * b ** 2)  # (16.30) corrected
    return dict(Qp=float(Qp), Qc=float(Qc), Q=float(Qp + Qc))  # Eq. (16.31) = Qp + Qc


def two_layer_total_flow_rate(a, delta, dp, L, mu_p, mu_c):
    """Total flow of the two-layer model in the book's closed form (16.31).  Label: analytic.

    Book: §16.3, Eq. (16.31).  Parameters as :func:`two_layer_flow_rates`.  Returns Q [m³/s].
    """
    return _S(np.pi * _F(a) ** 4 * _F(dp) / (8.0 * _F(mu_p) * _F(L))
              * (1.0 - (1.0 - _F(delta) / _F(a)) ** 4 * (1.0 - _F(mu_p) / _F(mu_c))))  # Eq. (16.31)


def apparent_viscosity_two_layer(a, delta, mu_p, mu_c):
    """Apparent viscosity of the two-layer flow by Poiseuille analogy:
    μ_app = μ_p[1 − (1 − δ/a)⁴(1 − μ_p/μ_c)]⁻¹.

    Book: §16.3, Eqs. (16.32)–(16.33).  μ_app < μ_c iff μ_c > μ_p: a thinner tube gives the slow plasma layer a bigger
    share (the Fåhræus–Lindqvist effect, derivation D07).
    Parameters
    ----------
    a [m];  delta [m];  mu_p, mu_c [Pa s].
    Returns
    -------
    mu_app [Pa s].
    Assumptions: as :func:`two_layer_velocity`.
    Validation: V1 = πa⁴Δp/(8QL) with :func:`two_layer_flow_rates`; V7 δ → 0 ⇒ μ_c, μ_c = μ_p ⇒ μ_p. Label: analytic.
    """
    return _S(_F(mu_p) / (1.0 - (1.0 - _F(delta) / _F(a)) ** 4 * (1.0 - _F(mu_p) / _F(mu_c))))  # Eq. (16.33)


def apparent_viscosity_thin_layer(a, delta, mu_p, mu_c):
    """Thin-layer limit μ_app ≈ μ_c[1 + 4(δ/a)(μ_c/μ_p − 1)]⁻¹ (from (1 − δ/a)⁴ ≈ 1 − 4δ/a).

    Book: §16.3, Eq. (16.34).
    Parameters: a [m];  delta [m];  mu_p, mu_c [Pa s].  Returns mu_app [Pa s].
    Validation: V7 relative error against :func:`apparent_viscosity_two_layer` ∝ (δ/a)². Label: analytic.
    """
    return _S(_F(mu_c) / (1.0 + 4.0 * _F(delta) / _F(a) * (_F(mu_c) / _F(mu_p) - 1.0)))  # Eq. (16.34)


def charm_kurland_core_viscosity(H_c, T, mu_p):
    """Core viscosity from the core hematocrit (Charm & Kurland): μ_c = μ_p/(1 − α_cH_c),
    α_c = 0.070 exp[2.49H_c + (1107/T)exp(−1.69H_c)].

    Book: §16.3, Eqs. (16.35)–(16.36) (stated valid to H_c = 0.60; T in kelvin).
    Parameters
    ----------
    H_c : core hematocrit [–] (fraction, 0 ≤ H_c ≤ 0.6 — raises above, the stated range);  T : temperature [K];
    mu_p : plasma viscosity [Pa s].
    Returns
    -------
    mu_c [Pa s].
    Assumptions: the correlation's range.  Validation: V1 hand value; H_c → 0 ⇒ μ_p. Label: analytic.
    """
    H_c = _F(H_c)
    if np.any(H_c > 1.0):
        raise ValueError("H_c is a fraction (trap T13)")
    if np.any(H_c > 0.6):
        raise ValueError("(16.36) is stated valid to H_c = 0.60")
    return _S(_charm_kurland(H_c, T, mu_p))


def _charm_kurland(H_c, T, mu_p):
    """(16.35)–(16.36) without the range check (used inside the root find of :func:`fahraeus_lindqvist_model`)."""
    H_c, T = _F(H_c), _F(T)
    alpha_c = 0.070 * np.exp(2.49 * H_c + 1107.0 / T * np.exp(-1.69 * H_c))  # Eq. (16.36)
    den = 1.0 - alpha_c * H_c
    if np.any(den <= 0.0):
        raise ValueError("alpha_c H_c >= 1: the correlation gives no finite core viscosity")
    return _F(mu_p) / den  # Eq. (16.35)


def fahraeus_lindqvist_model(a, delta, H_F, mu_p, T) -> dict:
    """Close the two-layer model (ours): find the core hematocrit H_c with H_c/H_F = 1 + Q_p/Q_c, μ_c = μ_c(H_c) from
    Charm–Kurland, then μ_app(a) from (16.33).

    Book: §16.3, text after (16.34) ("a modeling procedure can be developed") with (16.8), (16.29)–(16.31), (16.35),
    (16.36).  The book leaves δ and the procedure open; this closure and the choice of δ are ours (trap T19) — the
    measured comparison is :func:`pries_relative_viscosity`.
    Method: brentq on f(H) = H − H_F[1 + Q_p/Q_c(μ_c(H))] for H ∈ (H_F, H_max), H_max where α_cH = 1 (or 0.999); the
    residual is asserted.
    Parameters
    ----------
    a : radius [m];  delta : cell-free layer [m];  H_F : feed (= discharge) hematocrit [–];  mu_p [Pa s];  T [K].
    Returns
    -------
    dict: H_c, mu_c, mu_app, mu_rel (= μ_app/μ_p), H_T, Qp_over_Qc, residual, in_range (H_c ≤ 0.6, the stated
    validity of (16.36)).
    Assumptions: as :func:`two_layer_velocity`; the Charm–Kurland range.
    Validation: V1 residual < 1e-12; V7 μ_app falls as a decreases at fixed δ (the sign of the FL effect).
    Label: analytic (model, ours).
    """
    a, delta, H_F, mu_p, T = (float(v) for v in (a, delta, H_F, mu_p, T))
    _check_hematocrit(H_F)

    def ratio(H):
        mu_c = float(_charm_kurland(H, T, mu_p))
        q = two_layer_flow_rates(a, delta, 1.0, 1.0, mu_p, mu_c)
        return q["Qp"] / q["Qc"], mu_c

    def f(H):
        return H - H_F * (1.0 + ratio(H)[0])

    # upper end of the physical bracket: α_c H = 1
    def g(H):
        return 0.070 * np.exp(2.49 * H + 1107.0 / T * np.exp(-1.69 * H)) * H - 1.0

    H_hi = brentq(g, 1e-6, 0.999) * (1.0 - 1e-9) if g(0.999) > 0 else 0.999
    if f(H_hi) < 0:
        raise ValueError("no core hematocrit satisfies the cell balance: delta too large for this H_F")
    Hc = brentq(f, H_F, H_hi, xtol=1e-15, rtol=1e-15, maxiter=200)
    res = f(Hc)
    assert abs(res) < 1e-12, res
    r_pc, mu_c = ratio(Hc)
    mu_app = float(apparent_viscosity_two_layer(a, delta, mu_p, mu_c))
    return dict(H_c=Hc, mu_c=mu_c, mu_app=mu_app, mu_rel=mu_app / mu_p, H_T=Hc * (a - delta) ** 2 / a ** 2,
                Qp_over_Qc=r_pc, residual=res, in_range=bool(Hc <= 0.6))


def pries_relative_viscosity(D_um, H_D):
    """In-vitro relative apparent viscosity of blood in glass tubes (Pries et al. 1992), μ_rel(D, H_D) — a cited
    empirical law used as the measured counterpart of the Fåhræus–Lindqvist effect.

    μ_rel = 1 + (μ₀.₄₅ − 1)[(1 − H_D)^C − 1]/[(1 − 0.45)^C − 1],
    μ₀.₄₅ = 220 exp(−1.3D) + 3.2 − 2.44 exp(−0.06D^{0.645}),
    C = (0.8 + exp(−0.075D))(−1 + 1/(1 + 10⁻¹¹D¹²)) + 1/(1 + 10⁻¹¹D¹²), D in µm.
    Book: §16.2 (the FL effect described, Fig. 16.9); the law is not in the book.  Source: T. W. Secomb & A. R. Pries,
    C. R. Physique 14, 470–478 (2013), Eqs. (3)–(5) (read first-hand 2026-10-14; doi:10.1016/j.crhy.2013.04.002),
    reanalysing Pries et al., Am. J. Physiol. 263, H1770 (1992).  High-shear data only.
    Parameters
    ----------
    D_um : tube diameter [µm];  H_D : discharge hematocrit [–].
    Returns
    -------
    mu_rel = μ_app/μ_plasma [–].
    Assumptions: in vitro, high shear rates (ū/D ≥ 50 s⁻¹).
    Validation: V5 (the cited formula; μ_rel(·, 0.45) = μ₀.₄₅). Label: benchmark.
    """
    D, H = _F(D_um), _F(H_D)
    mu45 = 220.0 * np.exp(-1.3 * D) + 3.2 - 2.44 * np.exp(-0.06 * D ** 0.645)
    w = 1.0 / (1.0 + 1e-11 * D ** 12)
    C = (0.8 + np.exp(-0.075 * D)) * (-1.0 + w) + w
    return _S(1.0 + (mu45 - 1.0) * ((1.0 - H) ** C - 1.0) / ((1.0 - 0.45) ** C - 1.0))


# ======================================================================================================================
# §16.3 elastic tube: hoop stress, tube law, distensible Poiseuille flow, compliance (I8)
# ======================================================================================================================
def _beta_guard(beta, name="a0(p - pe)/(E h)"):
    """Tube-law singularity guard (trap T16): raise for scalars, NaN + warning for arrays."""
    beta = _F(beta)
    bad = beta >= 1.0
    if np.any(bad):
        if beta.ndim == 0:
            raise ValueError(f"{name} = {float(beta):.4g} >= 1: the Hookean tube law (16.41) has no finite radius "
                             "(radius → ∞ as this ratio → 1; trap T16)")
        warnings.warn(f"{name} >= 1 at some points: returning NaN there (trap T16)", stacklevel=3)
    return beta, bad


def hoop_stress(p, pe, a, h):
    """Circumferential (hoop) stress of a thin-walled tube, σ_θθ = (p − p_e)a/h (Laplace's law).

    Book: §16.3 "Effect of Tube Wall Elasticity on Poiseuille Flow", Eqs. (16.38)–(16.39) (half-tube force balance,
    ∫₀^π sin θ dθ = 2).
    Parameters
    ----------
    p : internal pressure [Pa];  pe : external pressure [Pa];  a : (deformed) radius [m];  h : wall thickness [m].
    Returns
    -------
    sigma [Pa].
    Assumptions: thin wall (h/a ≪ 1), radial stress neglected.  Validation: V1 the force balance by quadrature.
    Label: analytic.
    """
    return _S((_F(p) - _F(pe)) * _F(a) / _F(h))  # Eq. (16.39)


def tube_radius_hookean(p, pe, a0, E, h):
    """Pressure–radius relation (tube law) of a thin Hookean tube: a = a₀[1 − a₀(p − p_e)/(Eh)]⁻¹.

    Book: §16.3, Eqs. (16.40)–(16.41) (hoop strain (a − a₀)/a₀ = σ_θθ/E with σ_θθ from (16.39) at the deformed
    radius).  The radius → ∞ as a₀(p − p_e)/(Eh) → 1 (trap T16): raises for a scalar there, NaN for array entries.
    Parameters
    ----------
    p [Pa];  pe [Pa];  a0 : unstressed radius [m];  E : Young's modulus [Pa];  h : wall thickness [m].
    Argument order (positional): p, pe, a0, E, h — note the sibling :func:`tube_law_hookean` takes (A, A0, E, h, a0);
    call by keyword when mixing them.
    Returns
    -------
    a [m].
    Assumptions: thin wall, linear elastic (Hookean) wall, quasi-static.
    Regime: meant for distension, p ≥ p_e (a ≥ a₀).  For p < p_e the formula still returns a circular radius a < a₀,
    but a real thin tube under negative transmural pressure buckles and collapses non-circularly long before that
    (the collapsible-tube regime of §16.3, Shapiro's law (16.188) — :func:`collapsible_tube_law`); small negative
    values are a linear-elastic extrapolation only.
    Validation: V1 inverse of :func:`tube_law_hookean`; V7 E → ∞ ⇒ a₀. Label: analytic.
    """
    beta = _F(a0) * (_F(p) - _F(pe)) / (_F(E) * _F(h))
    beta, bad = _beta_guard(beta)
    a = _F(a0) / (1.0 - beta)  # Eq. (16.41)
    return _S(np.where(bad, np.nan, a))


def tube_law_hookean(A, A0, E, h, a0):
    """Transmural pressure from the area of a Hookean tube: p − p_e = (Eh/a₀)(1 − a₀/a) = (Eh/a₀)[1 − (A₀/A)^{1/2}].

    Book: §16.3 "Pulse Wave Propagation …", Eq. (16.55) (from (16.41)).
    Parameters
    ----------
    A : area [m²];  A0 : unstressed area πa₀² [m²];  E [Pa];  h [m];  a0 [m].
    Argument order (positional): A, A0, E, h, a0 (a₀ last — unlike :func:`tube_radius_hookean`'s p, pe, a0, E, h).
    Returns
    -------
    p − p_e [Pa].
    Assumptions: as :func:`tube_radius_hookean`.  Validation: V1 round trip with :func:`tube_radius_hookean`.
    Label: analytic.
    """
    return _S(_F(E) * _F(h) / _F(a0) * (1.0 - np.sqrt(_F(A0) / _F(A))))  # Eq. (16.55)


def tube_law_linear(A, A0, E, h, a0):
    """Linearised tube law p − p_e = (Eh/(2a₀A₀))(A − A₀), ∂p/∂A = Eh/(2a₀A₀).

    Book: §16.3, Eq. (16.58) (Taylor expansion of (16.55) about A₀).
    Parameters: A [m²];  A0 [m²];  E [Pa];  h [m];  a0 [m].  Returns p − p_e [Pa].
    Argument order (positional): A, A0, E, h, a0 (as :func:`tube_law_hookean`).
    Validation: V2 sympy series of (16.55). Label: symbolic.
    """
    return _S(_F(E) * _F(h) / (2.0 * _F(a0) * _F(A0)) * (_F(A) - _F(A0)))  # Eq. (16.58)


def elastic_tube_flow_rate(p1, p2, pe, a0, E, h, mu, L):
    """Steady flow through a distensible Hookean tube (local Poiseuille + tube law):
    Q = (πa₀³Eh/(24μL)){[1 − a₀(p₁ − p_e)/(Eh)]⁻³ − [1 − a₀(p₂ − p_e)/(Eh)]⁻³}.

    Book: §16.3, Eqs. (16.42)–(16.46).  DEVIATION: slip S4 — (16.44) and (16.45) print 8μ/(πa₀) where dimensions and
    (16.46) require 8μ/(πa₀⁴); (16.46) is right and is what is coded.  Evaluated in the stable form
    expm1(−3 log1p(−β₁)) − expm1(−3 log1p(−β₂)) so the rigid limit E → ∞ (πa₀⁴Δp/(8μL)) keeps its digits.
    Parameters
    ----------
    p1, p2 : inlet, outlet pressures [Pa];  pe : external pressure [Pa];  a0 [m];  E [Pa];  h [m];  mu [Pa s];  L [m].
    Returns
    -------
    Q [m³/s].
    Assumptions: steady, laminar, slowly varying radius (lubrication), thin Hookean wall.
    Validation: V2 sympy: d/dp of the bracket gives (16.44) with a₀⁴; V1 independent-route cross-check: solve_ivp of
    (16.43) along x (no refinement order asserted); V7 E → ∞ ⇒ (16.17). Label: symbolic, analytic.
    Argument order (positional): p1, p2, pe, a0, E, h, mu, L — call by keyword when in doubt.
    """
    a0, E, h, mu, L = (float(v) for v in (a0, E, h, mu, L))
    b1 = a0 * (_F(p1) - float(pe)) / (E * h)
    b2 = a0 * (_F(p2) - float(pe)) / (E * h)
    _beta_guard(np.maximum(b1, b2))
    diff = np.expm1(-3.0 * np.log1p(-b1)) - np.expm1(-3.0 * np.log1p(-b2))
    return _S(np.pi * a0 ** 3 * E * h / (24.0 * mu * L) * diff)  # Eq. (16.46)


def elastic_tube_pressure(x, Q, p1, pe, a0, E, h, mu):
    """Pressure along a distensible Hookean tube carrying Q: integrate the corrected (16.44) from 0 to x,
    [1 − β(x)]⁻³ = [1 − β₁]⁻³ − 24μQx/(πa₀³Eh), β = a₀(p − p_e)/(Eh).

    Book: §16.3, Eqs. (16.43)–(16.45) (with the a₀⁴ of slip S4 corrected).
    Parameters
    ----------
    x : position(s) from the inlet [m];  Q [m³/s];  p1 : inlet pressure [Pa];  pe [Pa];  a0 [m];  E [Pa];  h [m];
    mu [Pa s].
    Returns
    -------
    p(x) [Pa] (NaN where the tube would have collapsed to zero radius — the bracket ≤ 0).
    Assumptions: as :func:`elastic_tube_flow_rate`.
    Validation: V1 independent-route cross-check: solve_ivp of dp/dx = −8μQ/(πa(p)⁴) (16.43) (no refinement order
    asserted). Label: analytic.
    Argument order (positional): x, Q, p1, pe, a0, E, h, mu.
    """
    a0, E, h, mu, Q = (float(v) for v in (a0, E, h, mu, Q))
    b1 = a0 * (float(p1) - float(pe)) / (E * h)
    _beta_guard(b1)
    rhs = (1.0 - b1) ** -3 - 24.0 * mu * Q * _F(x) / (np.pi * a0 ** 3 * E * h)
    with np.errstate(invalid="ignore"):
        beta = 1.0 - np.where(rhs > 0.0, rhs, np.nan) ** (-1.0 / 3.0)
    return _S(float(pe) + E * h / a0 * beta)


def elastic_tube_resistance(p1, p2, pe, a0, E, h, mu, L):
    """Effective resistance (p₁ − p₂)/Q of a distensible tube — falls as the mean pressure rises (ours; the mechanism of
    Fig. 16.23).

    Book: §16.3, from (16.46); Fig. 16.23 (pulmonary resistance against pressure) described.
    Parameters as :func:`elastic_tube_flow_rate`.  Returns R_eff [Pa s/m³].
    Validation: V7 → 8μL/(πa₀⁴) as E → ∞. Label: analytic.
    """
    return _S((_F(p1) - _F(p2)) / _F(elastic_tube_flow_rate(p1, p2, pe, a0, E, h, mu, L)))


def compliance_from_tube(p, pe, a0, E, h, L):
    """Compliance K = dV/dp of a Hookean tube segment: K = 2πLa₀³/(Eh) [1 − a₀(p − p_e)/(Eh)]⁻³ (ours).

    Book: §16.3, compliance defined after (16.48) ("the change in volume for a change in pressure"); the closed form
    from V = πa²L and (16.41) is ours.
    Parameters
    ----------
    p [Pa];  pe [Pa];  a0 [m];  E [Pa];  h [m];  L : segment length [m].
    Returns
    -------
    K [m³/Pa].
    Assumptions: as :func:`tube_radius_hookean`.
    Validation: V1 central difference of πa(p)²L. Label: analytic.
    """
    beta = _F(a0) * (_F(p) - _F(pe)) / (_F(E) * _F(h))
    beta, bad = _beta_guard(beta)
    K = 2.0 * np.pi * _F(L) * _F(a0) ** 3 / (_F(E) * _F(h)) / (1.0 - beta) ** 3
    return _S(np.where(bad, np.nan, K))


# ======================================================================================================================
# §16.3 the Windkessel (I9)
# ======================================================================================================================
class HalfSineEjection:
    """Ventricular ejection Q(t) = Q_pk sin(πτ/T_s) for τ = t mod T < T_s, else 0; Q_pk = πSV/(2T_s) so each beat
    ejects SV (ours, the standard idealisation).  Callable; ``breakpoints(t0, t1)`` lists the kinks for quadrature.

    Book: §16.3 "Elasticity of the Aorta and the Windkessel Theory" (Q(t) the ejection rate in (16.47)–(16.52)); the
    waveform is ours.
    Parameters: SV [m³];  T : period [s];  Ts : systole duration [s] (0 < Ts ≤ T).
    Validation: V1 ∫ over one beat = SV. Label: analytic.
    """

    def __init__(self, SV: float, T: float, Ts: float):
        if not 0.0 < Ts <= T:
            raise ValueError("need 0 < Ts <= T")
        self.SV, self.T, self.Ts = float(SV), float(T), float(Ts)
        self.Q_peak = np.pi * self.SV / (2.0 * self.Ts)
        self.period = self.T

    def __call__(self, t):
        tau = np.mod(_F(t), self.T)
        return _S(np.where(tau < self.Ts, self.Q_peak * np.sin(np.pi * tau / self.Ts), 0.0))

    def breakpoints(self, t0: float, t1: float) -> list:
        k0, k1 = int(np.floor(t0 / self.T)) - 1, int(np.ceil(t1 / self.T)) + 1
        pts = []
        for k in range(k0, k1 + 1):
            pts += [k * self.T, k * self.T + self.Ts]
        return [p for p in pts if t0 < p < t1]


def ejection_half_sine(t, SV, T_period, Ts):
    """Half-sine ventricular ejection Q(t) (ours): Q_max sin(π(t mod T)/T_s) in systole, 0 in diastole,
    Q_max = πSV/(2T_s).

    Book: §16.3 (the Q(t) of (16.47)–(16.52)); waveform ours.
    Parameters: t [s];  SV [m³];  T_period [s];  Ts systole duration [s].  Returns Q [m³/s].
    Validation: V1 ∫ over one beat = SV. Label: analytic.
    """
    return HalfSineEjection(SV, T_period, Ts)(t)


def windkessel_rhs(t, p, Q_fn: Callable, R, K):
    """Right side of the two-element Windkessel, dp/dt = [Q(t) − p/R]/K.

    Book: §16.3, Eq. (16.49) Q(t) = K(dp/dt) + p(t)/R, from the volume balance (16.47) and the compliance chain rule
    (16.48) dv/dt = (dv/dp)(dp/dt), K = dv/dp (derivation D09).  Pressures are transmural with p̄_V = 0 (trap T12).
    Parameters
    ----------
    t [s];  p : aortic pressure [Pa];  Q_fn : inflow Q(t) [m³/s];  R : peripheral resistance [Pa s/m³];
    K : compliance [m³/Pa].
    Returns
    -------
    dp/dt [Pa/s].
    Assumptions: incompressible blood, instantaneous transmission, constant R and K.
    Validation: V1 constant inflow ⇒ exact exponential. Label: analytic.
    """
    return (float(Q_fn(t)) - p / R) / K  # Eq. (16.49) solved for dp/dt


def windkessel_ode(t, Q_fn: Callable, R, K, p0, *, rtol: float = 1e-12, atol: float | None = None,
                   max_step: float | None = None, method: str = "DOP853"):
    """Aortic pressure p(t) of the two-element Windkessel by numerical integration of (16.49) (solve_ivp, DOP853).

    Book: §16.3, Eq. (16.49).
    Parameters
    ----------
    t : output times [s] (increasing; the integration starts at t[0] with p = p0);  Q_fn : Q(t) [m³/s];  R [Pa s/m³];
    K [m³/Pa];  p0 : pressure at t[0] [Pa];  rtol, atol : solver tolerances (atol default 1e-12·max(|p0|, QR scale));
    max_step : [s] (default (t[-1] − t[0])/2000, or T/200 if Q_fn has a period);  method : solve_ivp method.
    Returns
    -------
    p [Pa] at t.
    Assumptions: as :func:`windkessel_rhs`.
    Validation: V1 independent-route cross-check: agrees with the closed form :func:`windkessel_pressure` (16.52) to
    1e-8 relative (no refinement order asserted). Label: analytic.
    """
    t = _F(t)
    if max_step is None:
        per = getattr(Q_fn, "period", None)
        max_step = per / 200.0 if per else (t[-1] - t[0]) / 2000.0
    scale = max(abs(float(p0)), abs(float(np.mean(np.atleast_1d(Q_fn(t))))) * R, 1.0)
    atol = 1e-12 * scale if atol is None else atol
    sol = solve_ivp(lambda tt, y: [windkessel_rhs(tt, y[0], Q_fn, R, K)], (t[0], t[-1]), [float(p0)], t_eval=t,
                    rtol=rtol, atol=atol, max_step=max_step, method=method)
    if not sol.success:
        raise RuntimeError(sol.message)
    return _S(sol.y[0])


def windkessel_pressure(t, Q_fn: Callable, R, K, p0, *, epsabs: float = 0.0, epsrel: float = 1e-12):
    """Aortic pressure of the Windkessel from the closed-form solution
    p(t) = (1/K)e^{−t/RK}∫₀ᵗQ(τ)e^{τ/RK}dτ + p₀e^{−t/RK}, by quadrature.

    Book: §16.3, Eqs. (16.50)–(16.52) (integrating factor; p₀ = end-diastolic aortic pressure; derivation D10).
    Evaluated by the equivalent stable recursion p_{k+1} = p_k e^{−Δ/RK} + (1/K)∫_{t_k}^{t_{k+1}} Q(τ)e^{−(t_{k+1}−τ)/RK}dτ
    (no e^{+t/RK} overflow for long runs), each step by ``quad`` split at ``Q_fn.breakpoints`` if present.
    Parameters
    ----------
    t : times [s], increasing, t[0] = 0 for the book's form;  Q_fn [m³/s];  R [Pa s/m³];  K [m³/Pa];  p0 : p(t[0]) [Pa];
    epsabs, epsrel : quad tolerances.
    Returns
    -------
    p [Pa] at t.
    Assumptions: as :func:`windkessel_rhs`.
    Validation: V1 constant Q ⇒ QR(1 − e^{−t/RK}) + p₀e^{−t/RK}; Q = 0 ⇒ p₀e^{−t/RK}; V1 = :func:`windkessel_ode`.
    Label: analytic.
    """
    t = _F(t)
    tau = float(R) * float(K)
    out = np.empty_like(t)
    out[0] = float(p0)
    bp = getattr(Q_fn, "breakpoints", None)
    for k in range(len(t) - 1):
        a, b = float(t[k]), float(t[k + 1])
        pts = bp(a, b) if bp else None
        integ, _ = quad(lambda s: float(Q_fn(s)) * np.exp(-(b - s) / tau), a, b, points=pts or None, limit=200,
                        epsabs=epsabs, epsrel=epsrel)
        out[k + 1] = out[k] * np.exp(-(b - a) / tau) + integ / float(K)  # Eq. (16.52), stepwise
    return _S(out)


def windkessel_periodic_p0_general(Q_fn: Callable, R, K, T):
    """End-diastolic pressure of the periodic Windkessel state for any periodic inflow, p₀ = I_T/(1 − e^{−T/RK}),
    I_T = response from rest over one period (ours: the fixed point p(T) = p(0) of a linear ODE).

    Book: §16.3, from (16.52) with p(T) = p₀.
    Parameters: Q_fn [m³/s] periodic with period T;  R [Pa s/m³];  K [m³/Pa];  T [s].  Returns p₀ [Pa].
    Validation: V1 = :func:`windkessel_periodic_p0` for half-sine ejection; V4 ∫Q dt = ∫p/R dt. Label: analytic.
    """
    I_T = float(np.atleast_1d(windkessel_pressure([0.0, float(T)], Q_fn, R, K, 0.0))[-1])
    return I_T / (1.0 - np.exp(-float(T) / (float(R) * float(K))))


def _wk_constants(SV, HR, systole_share, R, K):
    T = 60.0 / float(HR)
    Ts = float(systole_share) * T
    tau_c = float(R) * float(K)
    ws = np.pi / Ts
    Qmax = np.pi * float(SV) / (2.0 * Ts)
    B = float(R) * Qmax * ws * tau_c / (1.0 + (ws * tau_c) ** 2)  # design C.12 sign: B > 0
    E1, E2 = np.exp(-Ts / tau_c), np.exp(-(T - Ts) / tau_c)
    p0 = B * (1.0 + E1) * E2 / (1.0 - np.exp(-T / tau_c))
    return T, Ts, tau_c, ws, Qmax, B, p0


def windkessel_periodic_p0(*, SV, HR, systole_share, R, K) -> float:
    """End-diastolic pressure of the periodic Windkessel beat for half-sine ejection (closed form, ours):
    p₀ = B(1 + e^{−T_s/τ})e^{−(T−T_s)/τ}/(1 − e^{−T/τ}), B = RQ_maxω_sτ/(1 + (ω_sτ)²), τ = RK, ω_s = π/T_s.

    Book: §16.3, (16.49), (16.52) with p(T) = p₀ (derivation D10 step 8; design C.12).
    Parameters: SV [m³];  HR [beats/min];  systole_share T_s/T [–];  R [Pa s/m³];  K [m³/Pa].  Returns p₀ [Pa].
    Validation: V1 = :func:`windkessel_periodic_p0_general` (quadrature). Label: analytic.
    """
    return float(_wk_constants(SV, HR, systole_share, R, K)[6])


def windkessel_periodic_state(t, *, SV, HR, systole_share, R, K):
    """Closed-form periodic Windkessel pressure for half-sine ejection (ours; what the E2 explainer mirrors).

    Book: §16.3, (16.49) with Q = :class:`HalfSineEjection`.  With τ_c = RK, ω_s = π/T_s, s = t mod T and
    B = RQ_maxω_sτ_c/(1 + (ω_sτ_c)²) (design C.12): systole (s < T_s)
    p = RQ_max(sin ω_ss − ω_sτ_c cos ω_ss)/(1 + (ω_sτ_c)²) + (p₀ + B)e^{−s/τ_c}; diastole
    p = [B + (p₀ + B)e^{−T_s/τ_c}]e^{−(s−T_s)/τ_c}; periodicity p(T) = p₀ gives :func:`windkessel_periodic_p0`.
    Parameters
    ----------
    t [s];  SV [m³];  HR [beats/min];  systole_share : T_s/T [–];  R [Pa s/m³];  K [m³/Pa].
    Returns
    -------
    p [Pa].
    Assumptions: as :func:`windkessel_rhs`; the periodic state (all transients gone).
    Validation: V1 equals :func:`windkessel_ode` after many beats; V4 mean = R·SV·HR/60. Label: analytic.
    """
    T, Ts, tau_c, ws, Qmax, B, p0 = _wk_constants(SV, HR, systole_share, R, K)
    tt = _F(t)
    tau = np.atleast_1d(np.mod(tt, T)).astype(float)
    sys_ = tau < Ts
    p = np.empty_like(tau)
    # each branch is evaluated only where it applies, so its exponent is never positive (no overflow for tiny R*K)
    ts = tau[sys_]
    p[sys_] = float(R) * Qmax * (np.sin(ws * ts) - ws * tau_c * np.cos(ws * ts)) / (1.0 + (ws * tau_c) ** 2) \
        + (p0 + B) * np.exp(-ts / tau_c)  # Eq. (16.49) solved for systole
    td = tau[~sys_]
    p[~sys_] = (B + (p0 + B) * np.exp(-Ts / tau_c)) * np.exp(-(td - Ts) / tau_c)  # Eq. (16.49) solved for diastole
    return _S(p.reshape(np.shape(tt)))


def windkessel_summary(*, SV, HR, systole_share, R, K, n: int = 2001) -> dict:
    """Numbers of the periodic Windkessel beat (ours): p₀ (end-diastolic), p_S (max), p_D (min), the true mean
    (= R·Q̄, exactly), the ⅓-rule estimate (16.11) and its error, τ = RK, Q̄, time of the peak.

    Book: §16.3 (16.49) with §16.2 (16.9)–(16.11).
    Parameters as :func:`windkessel_periodic_state`; n samples per beat (the maximum is refined by a bounded scalar
    minimisation).
    Returns dict: p0, pS, pD, p_mean, map_third, map_error (= map_third − p_mean), pulse_pressure, tau, Q_mean, t_peak,
    T [Pa, s, m³/s].
    Validation: V4 sampled mean = R·SV·HR/60 to 1e-9. Label: analytic.
    """
    from scipy.optimize import minimize_scalar

    T, Ts, tau_c, ws, Qmax, B, p0 = _wk_constants(SV, HR, systole_share, R, K)
    t = np.linspace(0.0, T, int(n))
    kw = dict(SV=SV, HR=HR, systole_share=systole_share, R=R, K=K)
    p = _F(windkessel_periodic_state(t, **kw))
    i = int(np.argmax(p))
    lo, hi = t[max(i - 1, 0)], t[min(i + 1, len(t) - 1)]
    res = minimize_scalar(lambda s: -float(windkessel_periodic_state(s, **kw)), bounds=(lo, hi), method="bounded",
                          options=dict(xatol=1e-12))
    pS, t_peak = -float(res.fun), float(res.x)
    pD = float(p.min())
    Q_mean = float(SV) * float(HR) / 60.0
    p_mean = float(R) * Q_mean  # exact: ∫Q dt = ∫p/R dt over a periodic beat
    third = float(mean_arterial_pressure(pS, pD))
    return dict(p0=float(p0), pS=pS, pD=pD, p_mean=p_mean, p_mean_sampled=time_mean(p, t), map_third=third,
                map_error=third - p_mean, pulse_pressure=pS - pD, tau=tau_c, Q_mean=Q_mean, t_peak=t_peak, T=T)


windkessel_periodic_summary = windkessel_summary  # alias (curation name)


def windkessel_terms(t, *, SV, HR, systole_share, R, K) -> dict:
    """Where the inflow goes at an instant of the periodic beat: storage K dp/dt = Q − p/R and run-off p/R (E2 term
    bars; ours).

    Book: §16.3, Eq. (16.49) read term by term ((16.47): inflow = storage + outflow).
    Parameters: t [s];  SV [m³];  HR [beats/min];  systole_share [–];  R [Pa s/m³];  K [m³/Pa].
    Returns dict Q, storage, runoff [m³/s], p [Pa], dpdt [Pa/s].
    Validation: V1 storage = K dp/dt of :func:`windkessel_periodic_state` (central difference). Label: analytic.
    """
    T = 60.0 / float(HR)
    Q = _F(ejection_half_sine(t, SV, T, float(systole_share) * T))
    p = _F(windkessel_periodic_state(t, SV=SV, HR=HR, systole_share=systole_share, R=R, K=K))
    run = p / float(R)
    return dict(Q=_S(Q), storage=_S(Q - run), runoff=_S(run), p=_S(p), dpdt=_S((Q - run) / float(K)))


def windkessel_terms_general(t, p, Q_fn: Callable, R, K) -> dict:
    """:func:`windkessel_terms` for any inflow Q_fn and a pressure p on its solution (ours).  Returns dict Q, storage,
    runoff [m³/s], dpdt [Pa/s].  Book: §16.3, (16.49).  Label: analytic."""
    Q = _F(Q_fn(t))
    run = _F(p) / float(R)
    return dict(Q=_S(Q), storage=_S(Q - run), runoff=_S(run), dpdt=_S((Q - run) / float(K)))


# ======================================================================================================================
# §16.3 inviscid pulse waves (I10)
# ======================================================================================================================
def moens_korteweg_speed(E, h, rho, a0, *, midwall: bool = False):
    """Moens–Korteweg pulse-wave speed c = √(Eh/(2ρa₀)) (thin wall), or the mid-wall form √(Eh/(2ρ(a₀ + h/2))).

    Book: §16.3, Eqs. (16.61)–(16.62) (from the linearised tube equations (16.56)–(16.58); derivation D12) and (16.63)
    (Fung's mid-wall strain).  Eh/(2ρa) = Eh/(ρd) (trap T17).
    Parameters
    ----------
    E : Young's modulus [Pa];  h : wall thickness [m];  rho : blood density [kg/m³];  a0 : radius [m];
    midwall : use (16.63).
    Returns
    -------
    c [m/s].
    Assumptions: thin Hookean wall, inviscid, long waves, small amplitude.
    Validation: V1 = :func:`bramwell_hill_speed` with ∂p/∂A of (16.55) at A₀; V5 form of Wikipedia
    "Moens–Korteweg equation". Label: analytic, benchmark.
    """
    a_eff = _F(a0) + (0.5 * _F(h) if midwall else 0.0)
    return _S(np.sqrt(_F(E) * _F(h) / (2.0 * _F(rho) * a_eff)))  # Eq. (16.62) / (16.63)


def bramwell_hill_speed(A, dpdA, rho):
    """Long-wave speed from any tube law: c² = (A/ρ) dp/dA (Bramwell–Hill form).

    Book: §16.3, Eq. (16.62) (second member) and Eq. (16.185) for a general tube law P̂(A).
    Parameters: A [m²];  dpdA : d(p − p_e)/dA [Pa/m²];  rho [kg/m³].  Returns c [m/s].
    Validation: V1. Label: analytic.
    """
    return _S(np.sqrt(_F(A) * _F(dpdA) / _F(rho)))  # Eq. (16.62), (16.185)


def youngs_modulus_from_pwv(c, h, rho, a0, *, midwall: bool = False):
    """Invert Moens–Korteweg: E = 2ρa₀c²/h (or 2ρ(a₀ + h/2)c²/h) — the non-invasive stiffness estimate.

    Book: §16.3, after (16.65) (Mazumdar's transit-time idea) with (16.62)/(16.63).
    Parameters: c [m/s];  h [m];  rho [kg/m³];  a0 [m];  midwall.  Returns E [Pa].
    Validation: V1 round trip. Label: analytic.
    """
    a_eff = _F(a0) + (0.5 * _F(h) if midwall else 0.0)
    return _S(2.0 * _F(rho) * a_eff * _F(c) ** 2 / _F(h))


def linearisation_ratios(p_amp, p_mean, u, c, *, Eh_over_a0: float | None = None) -> dict:
    """The smallness ratios that justify the linearised tube equations: pressure amplitude/mean, u/c and (if the wall
    is given) the hoop strain p_amp/(Eh/a₀).

    Book: §16.3, the conditions stated before (16.56) (values there are the book's; inputs here are the reader's).
    Parameters: p_amp [Pa];  p_mean [Pa];  u [m/s];  c [m/s];  Eh_over_a0 [Pa] (optional).
    Returns dict p_ratio, u_over_c, strain (None without the wall) [–].
    Label: analytic.
    """
    strain = None if Eh_over_a0 is None else float(p_amp) / float(Eh_over_a0)
    return dict(p_ratio=float(p_amp) / float(p_mean), u_over_c=float(u) / float(c), strain=strain)


def pulse_wave(x, t, c, f1: Callable | None = None, f2: Callable | None = None):
    """d'Alembert solution in the book's *time-shape* form p = f₁(t − x/c) + f₂(t + x/c).

    Book: §16.3, Eq. (16.65) (f₂ = 0 for a wave travelling in +x only); evaluated through ch15 ``AC.dalembert`` with
    the space shapes f(s) = f₁(−s/c), g(s) = f₂(s/c).
    Parameters
    ----------
    x [m];  t [s];  c [m/s];  f1, f2 : callables of time [s] → pressure [Pa] (None = absent).
    Returns
    -------
    p [Pa].
    Assumptions: linear waves in a uniform thin elastic tube.
    Validation: V1 = (16.66) with f(x − ct) = f₁(t − x/c). Label: analytic.
    """
    f = None if f1 is None else (lambda s: f1(-_F(s) / c))
    g = None if f2 is None else (lambda s: f2(_F(s) / c))
    return AC.dalembert(x, t, c, f, g)  # Eq. (16.65)


def velocity_from_pulse(p, rho, c, *, direction: int = +1):
    """Fluid velocity carried by a travelling pressure wave: u = ±p/(ρc) (p̂₁ = ρcû₁ for +x, p̂₂ = −ρcû₂ for −x).

    Book: §16.3, Eq. (16.68) (from (16.57) applied to f(x ∓ ct); the water-hammer relation, ch15 (15.16)).
    Parameters: p : pressure of the wave [Pa];  rho [kg/m³];  c [m/s];  direction : +1 (travelling in +x) or −1.
    Returns u [m/s].
    Validation: V1 ρ∂u/∂t = −∂p/∂x (16.57) by finite differences on :func:`pulse_wave`. Label: analytic.
    """
    if direction not in (1, -1):
        raise ValueError("direction must be +1 or -1")
    return _S(direction * _F(p) / (_F(rho) * _F(c)))  # Eq. (16.68)


def velocity_from_pulse_shapes(x, t, c, rho, f1: Callable | None = None, f2: Callable | None = None):
    """Fluid velocity of the two-way d'Alembert pulse in time-shape form: u = [f₁(t − x/c) − f₂(t + x/c)]/(ρc).

    Book: §16.3, Eqs. (16.65)–(16.68); via ch15 ``AC.plane_wave_velocity``.  Parameters as :func:`pulse_wave` plus
    rho [kg/m³].  Returns u [m/s].  Label: analytic.
    """
    f = None if f1 is None else (lambda s: f1(-_F(s) / c))
    g = None if f2 is None else (lambda s: f2(_F(s) / c))
    return AC.plane_wave_velocity(x, t, c, rho, f, g)  # Eq. (16.68)


def wall_bulge(p, *, a0, E, h, pe: float = 0.0):
    """Radial wall displacement η = a − a₀ for a transmural pressure, from the Hookean tube law (16.41) (E3; ours).

    Book: §16.3, Eq. (16.41).  Small-amplitude value ≈ a₀²(p − p_e)/(Eh).
    Parameters: p [Pa];  a0 [m];  E [Pa];  h [m];  pe [Pa].  Returns eta [m].  Label: analytic.
    """
    return _S(_F(tube_radius_hookean(p, pe, a0, E, h)) - float(a0))


def pulse_state(x, t, *, c, amplitude, width, rho, a0, E, h, ic: str = "right", x0: float = 0.0,
                shape: str = "gaussian") -> dict:
    """A pressure pulse in an artery and what it does: p, u and the wall radius a from the tube law (E3 parity; ours).

    Book: §16.3, (16.65)–(16.68) and (16.41).  ``ic="right"``: one right-going pulse p = amplitude·g((x − x₀ − ct)/width),
    u = p/(ρc); ``ic="rest"``: the same bump released from rest splits into two halves travelling ±x,
    p = ½[g(x − ct) + g(x + ct)], u = ½[g(x − ct) − g(x + ct)]/(ρc).  g = ``AC.pulse(shape)`` (Gaussian by default).
    Parameters
    ----------
    x [m];  t [s];  c [m/s];  amplitude [Pa];  width [m];  rho [kg/m³];  a0 [m];  E [Pa];  h [m];  ic : "right" | "rest";
    x0 : initial centre [m];  shape : "gaussian" | "tophat" | "packet" | "sine".
    Returns
    -------
    dict p [Pa] (transmural), u [m/s], a [m].
    Regime / reading: the pulse pressure is treated as the *whole* transmural pressure (p_e = 0 and no mean pressure
    in the tube law), so ``a0`` is the radius at the mean (diastolic) pressure about which the pulse is a perturbation,
    not the unstressed radius; a(x, t) shows the wall's response to the pulse alone.  Linear acoustics: amplitude ≪ ρc².
    Validation: V1 u = p/(ρc) for the right-going pulse; a from :func:`tube_radius_hookean`. Label: analytic.
    """
    if ic == "right":
        f = AC.pulse(shape, width, amplitude, x0)
        g = None
    elif ic == "rest":
        f = AC.pulse(shape, width, 0.5 * amplitude, x0)
        g = f
    else:
        raise ValueError('ic must be "right" or "rest"')
    p = AC.dalembert(x, t, c, f, g)  # Eq. (16.66)
    u = AC.plane_wave_velocity(x, t, c, rho, f, g)  # Eq. (16.67)–(16.68)
    a = tube_radius_hookean(p, 0.0, a0, E, h)  # Eq. (16.41)
    return dict(p=p, u=u, a=a)


# ======================================================================================================================
# Nonlinear 1-D tube model (ours, I11): ∂A/∂t + ∂Q/∂x = 0, ∂Q/∂t + ∂(Q²/A + Π(A)/ρ)/∂x = 0
# ======================================================================================================================
def tube1d_pressure_flux(A, A0, E, h, a0, *, printed_16_54: bool = False):
    """Pressure part Π(A) of the momentum flux of the 1-D tube equations with the Hookean law (16.55) (ours):
    corrected A∂p/∂x = ∂Π/∂x with Π = (Eh/a₀)√(A₀A); printed (16.54) Π = (p − p_e)A = (Eh/a₀)(A − √(A₀A)).

    Book: §16.3, Eqs. (16.53)–(16.55).  DEVIATION: slip S5 — the printed right side −∂((p − p_e)A)/∂x of (16.54) is
    physically −A∂p/∂x (the wall's axial reaction cancels (p − p_e)∂A/∂x); both give c² = Eh/(2ρa₀) once linearised.
    Parameters: A [m²];  A0 [m²];  E [Pa];  h [m];  a0 [m];  printed_16_54.  Returns Π [N/m = Pa m²].
    Argument order (positional): A, A0, E, h, a0 (as :func:`tube_law_hookean`).
    Label: analytic.
    """
    A = _F(A)
    k = float(E) * float(h) / float(a0)
    if printed_16_54:
        return _S(k * (A - np.sqrt(float(A0) * A)))
    return _S(k * np.sqrt(float(A0) * A))


def tube1d_run(x, A_init, u_init, A0, E, h, a0, rho, t_end, *, inflow: Callable | None = None, CFL: float = 0.8,
               nu_art: float = 0.02, printed_16_54: bool = False, bc: str | None = None, eps_av: float = 0.0,
               snapshots=(), fast: bool = False) -> dict:
    """March the nonlinear one-dimensional elastic-tube equations with MacCormack's scheme (**ours**):
    U = (A, Q = Au), ∂U/∂t + ∂E(U)/∂x = 0, E = (Q, Q²/A + Π(A)/ρ), Π from :func:`tube1d_pressure_flux`.

    Book: §16.3, Eqs. (16.53) (continuity), (16.54) (momentum, corrected — slip S5) and (16.55) (tube law); the
    numerical method is ch10's MacCormack step (``core.maccormack.maccormack_step``, (10.100)–(10.102)), alternating
    F/B and B/F.  DEVIATION: the scheme, the conservative form and the artificial damping are ours — the book stops at
    the linear wave equation.  Damping: after each step U ← U − ν_art Δ⁴U (fourth difference: third-order small for a
    smooth pulse, so the scheme stays second order, and it conserves Σ A Δx on a periodic grid); optionally the
    pressure-switched second-difference smoothing of ch15 ``euler1d_run`` (``eps_av`` > 0) for steepened fronts.
    Inflow (optional): ``inflow(t)`` is the inlet transmural pressure [Pa]; the inlet area follows from (16.55) and the
    inlet velocity from the outgoing Riemann invariant u − 4(c₀ − c(A)) (c = c₀(A₀/A)^{1/4} for the Hookean law); the
    outlet is transmissive.
    Parameters
    ----------
    x : uniform nodes [m];  A_init [m²], u_init [m/s] : initial fields;  A0 [m²];  E [Pa];  h [m];  a0 [m];
    rho [kg/m³];  t_end [s];  inflow : callable t → inlet p − p_e [Pa] or None;  CFL : Courant number on |u| + c (≤ 1);
    nu_art : fourth-difference damping coefficient (≤ 1/16);  printed_16_54 : use the printed momentum flux;
    bc : "periodic" | "transmissive" (default periodic, or inflow/transmissive when ``inflow`` is given);
    eps_av : second-difference smoothing coefficient (0 = off);  snapshots : times [s] at which to store the fields;
    fast : unused (the caller picks the grid; common notebook signature).
    Returns
    -------
    dict: x, t (final time), A, u, p (transmural, (16.55)) at t, n_steps, mass, mass0 (Σ A Δx), momentum, momentum0,
    snapshots = dict(t, A, u, p) with arrays of shape (n_snap, n).
    Cost: ≈ 0.02 s for 400 nodes over one domain crossing.
    Assumptions: inviscid, 1-D, Hookean thin wall, p_e constant.
    Validation: V3 a small pulse travels at c(A₀) (16.62) within 0.5 %; V4 Σ A Δx conserved on a periodic domain;
    V3 order ≈ 2 for a smooth pulse. Label: converged, conserved.
    """
    from .core.maccormack import maccormack_step

    x = _F(x)
    dx = float(x[1] - x[0])
    if not np.allclose(np.diff(x), dx, rtol=1e-9, atol=0.0):
        raise ValueError("tube1d_run needs a uniform grid")
    if CFL > 1.0:
        raise ValueError("CFL must be <= 1")
    if nu_art < 0.0 or nu_art > 1.0 / 16.0:
        raise ValueError("nu_art must be in [0, 1/16] (stability of the fourth-difference damping)")
    if bc is None:
        bc = "transmissive" if inflow is not None else "periodic"
    if bc not in ("periodic", "transmissive"):
        raise ValueError('bc must be "periodic" or "transmissive"')
    A0, E, h, a0, rho = (float(v) for v in (A0, E, h, a0, rho))
    k = E * h / a0
    c0 = np.sqrt(k / (2.0 * rho))
    A = _F(A_init) + 0.0 * x
    u = _F(u_init) + 0.0 * x
    U = np.array([A, A * u])
    snap_t = sorted(float(s) for s in snapshots)
    snaps = dict(t=[], A=[], u=[], p=[])

    def flux(W):
        return np.array([W[1], W[1] ** 2 / W[0] + _F(tube1d_pressure_flux(W[0], A0, E, h, a0,
                                                                         printed_16_54=printed_16_54)) / rho])

    def speed(W):
        Aa = W[0]
        if printed_16_54:
            c2 = k / rho * (1.0 - 0.5 * np.sqrt(A0 / Aa))
        else:
            c2 = k / (2.0 * rho) * np.sqrt(A0 / Aa)
        return np.abs(W[1] / Aa) + np.sqrt(np.maximum(c2, 0.0))

    def apply_bc(W, tt):
        if bc == "transmissive":
            W[:, -1] = W[:, -2]
            if inflow is None:
                W[:, 0] = W[:, 1]
            else:  # prescribed inlet pressure; velocity from the outgoing invariant u − 4(c0 − c)
                beta = a0 * float(inflow(tt)) / (E * h)
                if beta >= 1.0:
                    raise ValueError("inflow pressure beyond the tube-law singularity (trap T16)")
                A_in = A0 / (1.0 - beta) ** 2  # (16.41) with A = πa²
                A1, u1 = W[0, 1], W[1, 1] / W[0, 1]
                c1, cin = c0 * (A0 / A1) ** 0.25, c0 * (A0 / A_in) ** 0.25
                u_in = (u1 - 4.0 * (c0 - c1)) + 4.0 * (c0 - cin)
                W[0, 0], W[1, 0] = A_in, A_in * u_in
        return W

    def store(W, tt):
        snaps["t"].append(tt)
        snaps["A"].append(W[0].copy())
        snaps["u"].append(W[1] / W[0])
        snaps["p"].append(_F(tube_law_hookean(W[0], A0, E, h, a0)))

    tot0 = U.sum(axis=1) * dx
    t, n = 0.0, 0
    U = apply_bc(U, 0.0) if bc == "transmissive" else U
    while snap_t and snap_t[0] <= 0.0:
        store(U, snap_t.pop(0))
    t_end = float(t_end)
    while t < t_end * (1.0 - 1e-14):
        if np.any(U[0] <= 0.0):
            raise FloatingPointError("tube1d_run: non-positive area (reduce CFL or raise the damping)")
        target = min(t_end, snap_t[0]) if snap_t else t_end
        dt = min(CFL * dx / float(np.max(speed(U))), target - t)
        U = maccormack_step(U, flux, None, dt, dx, arrangement="F/B" if n % 2 == 0 else "B/F")
        if nu_art > 0.0:  # DEVIATION: fourth-difference damping (ours)
            d4 = (np.roll(U, -2, axis=1) - 4.0 * np.roll(U, -1, axis=1) + 6.0 * U - 4.0 * np.roll(U, 1, axis=1)
                  + np.roll(U, 2, axis=1))
            if bc == "transmissive":
                d4[:, :2] = 0.0
                d4[:, -2:] = 0.0
            U = U - nu_art * d4
        if eps_av > 0.0:  # DEVIATION: pressure-switched smoothing (ours), as ch15 euler1d_run
            Aa = U[0]
            s = np.abs(np.roll(Aa, -1) - 2.0 * Aa + np.roll(Aa, 1)) / (np.roll(Aa, -1) + 2.0 * Aa + np.roll(Aa, 1))
            nu_ = np.minimum(0.45, eps_av * np.maximum(s, np.roll(s, -1)))
            d = nu_[None, :] * (np.roll(U, -1, axis=1) - U)
            if bc == "transmissive":
                d[:, -1] = 0.0
            U = U + d - np.roll(d, 1, axis=1)
        t += dt
        n += 1
        if bc == "transmissive":
            U = apply_bc(U, t)
        while snap_t and snap_t[0] <= t * (1.0 + 1e-12):
            store(U, snap_t.pop(0))
    tot = U.sum(axis=1) * dx
    A, Q = U
    snapshots_out = {key: np.array(val) for key, val in snaps.items()}
    return dict(x=x, A=A, u=Q / A, p=_F(tube_law_hookean(A, A0, E, h, a0)), t=t, n_steps=n, mass=tot[0],
                mass0=tot0[0], momentum=tot[1], momentum0=tot0[1], snapshots=snapshots_out)


# ======================================================================================================================
# §16.3 viscous waves in an elastic tube (I13) — time convention e^{i(kx − ωt)} (trap T1)
# ======================================================================================================================
def elastic_tube_characteristic(Xk, alpha, poisson=0.5):
    """Residual of the characteristic equation (16.136) of long viscous waves in a thin elastic tube:
    X²[2β̂J₀(β̂)/J₁(β̂) − 4] + X[4ν̂ − 1 − 2β̂J₀(β̂)/J₁(β̂)] + (1 − ν̂²), X = (k²/ω²)(Eh/(2ρa₀)), β̂ = i^{1/2}α.

    Book: §16.3 "Wave Propagation in a Viscous Liquid Contained in an Elastic Cylindrical Tube", Eqs. (16.134)–(16.136)
    (verified symbolically: the 2 × 2 determinant of (16.134)–(16.135) is (16.136), derivation D21).  Convention
    e^{i(kx − ωt)} (16.116): β = √(iω/ν), β̂ = βa₀ = e^{iπ/4}α (trap T1 — not the i^{3/2} of the rigid-tube flow).
    Parameters
    ----------
    Xk : X = k²c₀²/ω² (complex), c₀ the Moens–Korteweg speed;  alpha : Womersley number [–];
    poisson : Poisson ratio ν̂ (may be complex, (16.137)).
    Returns
    -------
    residual (complex) [–].
    Assumptions: a/λ ≪ 1, h/a ≪ 1, wall inertia neglected, linear, thin Hookean membrane.
    Validation: V2 sympy determinant; V3 :func:`elastic_tube_determinant_numeric` vanishes at its roots. Label: symbolic.
    """
    bh = WOM.I12 * float(alpha)
    Lr = 2.0 * bh * complex(WOM.bessel_ratio_j0_j1(bh))
    X = np.asarray(Xk, dtype=complex)
    nu_ = complex(poisson)
    return WOM._S_c(X ** 2 * (Lr - 4.0) + X * (4.0 * nu_ - 1.0 - Lr) + (1.0 - nu_ ** 2))  # Eq. (16.136)


def _char_roots(alpha, poisson):
    """Both roots of (16.136) for scalar or array α: array of shape (2,) + shape(α)."""
    bh = WOM.I12 * np.asarray(alpha, dtype=float)
    Lr = 2.0 * bh * np.asarray(WOM.bessel_ratio_j0_j1(bh), dtype=complex)
    nu_ = complex(poisson)
    a2, a1, a0_ = Lr - 4.0, 4.0 * nu_ - 1.0 - Lr, 1.0 - nu_ ** 2
    disc = np.sqrt(a1 * a1 - 4.0 * a2 * a0_ + 0j)
    return np.array([(-a1 + disc) / (2.0 * a2), (-a1 - disc) / (2.0 * a2)])


def _track_pressure_root(alphas_desc, poisson):
    """Follow the pressure root down a decreasing α path (vectorised roots, scalar tracking)."""
    R = _char_roots(alphas_desc, poisson)  # (2, n)
    out = np.empty(R.shape[1], dtype=complex)
    prev = R[np.argmin(np.abs(R[:, 0] - 1.0)), 0]
    for i in range(R.shape[1]):
        r0, r1 = R[0, i], R[1, i]
        prev = r0 if abs(r0 - prev) <= abs(r1 - prev) else r1
        out[i] = prev
    return out, R


def elastic_tube_roots(alpha, poisson=0.5, *, alpha_start: float = 200.0, n_steps: int = 400) -> dict:
    """Both roots X of (16.136), labelled by continuation from large α (trap T7): "pressure" → 1 as α → ∞ (the
    Moens–Korteweg wave, c → c₀) and "wall" → 0 (a fast, heavily damped wall-dominated mode).

    Book: §16.3, Eq. (16.136) and the text after it.  Root selection is ours: at α ≥ ``alpha_start`` the pressure root
    is the one nearer 1; below, follow it down a log-spaced path in α choosing the nearest root at each step.
    Parameters
    ----------
    alpha : Womersley number [–] (> 0);  poisson : ν̂ (may be complex);  alpha_start : continuation start;
    n_steps : continuation steps.
    Returns
    -------
    dict(pressure=X_p, other=X_w) complex [–] (``wall`` is an alias of ``other``).
    Assumptions: as :func:`elastic_tube_characteristic`.
    Validation: V7 X_p → 1 with 1 − c/c₀ ∝ 1/α; residuals of both roots < 1e-12. Label: analytic.
    """
    alpha = float(alpha)
    require_positive("alpha", alpha)
    if alpha >= alpha_start:
        r = _char_roots(alpha, poisson)
        xp = r[np.argmin(np.abs(r - 1.0))]
    else:
        path = np.geomspace(alpha_start, alpha, int(n_steps))
        trk, R = _track_pressure_root(path, poisson)
        xp, r = trk[-1], R[:, -1]
    xw = r[np.argmax(np.abs(r - xp))]
    return dict(pressure=complex(xp), other=complex(xw), wall=complex(xw))


def elastic_tube_wave_curve(alphas, c0: float = 1.0, poisson=0.5, *, branch: str = "pressure",
                            E_ratio: complex = 1.0 + 0.0j, alpha_start: float = 200.0, n_per_decade: int = 200) -> dict:
    """:func:`elastic_tube_wave` for a whole grid of α in one continuation (fast; for figures and sliders — ours).

    Book: §16.3, (16.136)–(16.137).  Parameters: alphas [–] (any order, > 0);  c0 [m/s];  poisson;  branch
    ("pressure" | "other");  E_ratio;  alpha_start, n_per_decade : continuation path.
    Returns dict of arrays aligned with ``alphas``: X, c_over_c0, k2_over_k1, damping_per_wavelength, c [m/s].
    Validation: V1 equals :func:`elastic_tube_wave` point by point. Label: analytic.
    """
    al = np.asarray(alphas, dtype=float)
    lo = float(al.min())
    n = max(int(n_per_decade * np.log10(alpha_start / min(lo, alpha_start))) + 2, 2)
    path = np.unique(np.concatenate([np.geomspace(alpha_start, min(lo, alpha_start), n), al[al < alpha_start]]))[::-1]
    trk, R = _track_pressure_root(path, poisson)
    Xp = np.empty(al.shape, dtype=complex)
    Xo = np.empty(al.shape, dtype=complex)
    for idx, a in np.ndenumerate(al):
        if a >= alpha_start:
            r = _char_roots(a, poisson)
            xp = r[np.argmin(np.abs(r - 1.0))]
        else:
            j = int(np.argmin(np.abs(path - a)))
            xp, r = trk[j], R[:, j]
        Xp[idx] = xp
        Xo[idx] = r[np.argmax(np.abs(r - xp))]
    X = (Xp if branch == "pressure" else Xo) / complex(E_ratio)
    s = np.sqrt(X + 0j)
    s = np.where(s.real < 0, -s, s)
    k2k1 = s.imag / s.real
    return dict(X=X, c_over_c0=1.0 / s.real, k2_over_k1=k2k1, damping_per_wavelength=np.exp(-2.0 * np.pi * k2k1),
                c=float(c0) / s.real)


def elastic_tube_wave(alpha, c0: float = 1.0, poisson=0.5, *, omega: float = 1.0, branch: str = "pressure",
                      E_ratio: complex = 1.0 + 0.0j) -> dict:
    """Wave speed and damping of a viscous wave in an elastic (or viscoelastic) tube from (16.136).

    Book: §16.3, Eq. (16.136) ("the wave speed ω/k₁ and the damping factor … from the real and imaginary parts of k/ω")
    and Eq. (16.137) for a viscoelastic wall: E* = E − iωE′, ν̂* = ν̂ − iων̂′ (Kelvin–Voigt in the e^{−iωt} convention).
    With E* the quadratic is solved for X* = k²E*h/(2ρa₀ω²) and X = X*/E_ratio, E_ratio = E*/E, so c₀ stays the
    elastic Moens–Korteweg speed.  k/ω = √X/c₀ with Re√X > 0 (a wave travelling in +x); k = k₁ + ik₂, decaying as
    e^{−k₂x}.  (The book prints |k| = √(k₁ + k₂) — slip S13; it is √(k₁² + k₂²).)
    Parameters
    ----------
    alpha : Womersley number [–];  c0 : Moens–Korteweg speed √(Eh/(2ρa₀)) [m/s];  poisson : ν̂ (or ν̂*, complex);
    omega : ω [rad/s] (for k₁, k₂; k₀ = ω/c₀);  branch : "pressure" | "other" ("wall" accepted);
    E_ratio : E*/E = 1 − iωE′/E (complex).
    Defaults: c0 = 1 and omega = 1 are placeholders that make the dimensionless outputs (X, c_over_c0, k2_over_k1,
    damping_per_wavelength) correct for any tube; c, k1, k2 are dimensional only when *both* c0 [m/s] and omega
    [rad/s] are set to the physical values (with the defaults they are in units of c₀ and ω/c₀).
    Returns
    -------
    dict: X (complex), c [m/s], c_over_c0 [–], k1, k2 [1/m], k2_over_k1 [–], damping_per_wavelength = e^{−2πk₂/k₁}
    [–], branch.
    Assumptions: as :func:`elastic_tube_characteristic`.
    Validation: V7 α → ∞: c → c₀ (1 − c/c₀ ∝ α⁻¹), k₂/k₁ ∝ α⁻¹; E_ratio = 1 recovers the elastic case exactly.
    Label: analytic.
    """
    Er = complex(E_ratio)
    roots = elastic_tube_roots(alpha, poisson)
    if branch not in roots:
        raise ValueError('branch must be "pressure" or "other"')
    Xstar = roots[branch]
    X = Xstar / Er
    s = np.sqrt(X + 0j)  # k/k₀, k₀ = ω/c₀
    if s.real < 0:
        s = -s
    c = float(c0) / s.real
    k2k1 = s.imag / s.real
    return dict(X=complex(X), c=float(c), c_over_c0=float(1.0 / s.real), k1=float(omega) * s.real / float(c0),
                k2=float(omega) * s.imag / float(c0), k2_over_k1=float(k2k1),
                damping_per_wavelength=float(np.exp(-2.0 * np.pi * k2k1)), branch=branch)


def _elastic_tube_modes(r, A, eta, k, omega, nu, E, h, a0, rho, poisson):
    """Closed-form long-wave mode shapes (û, v̂) of the viscous elastic tube, convention e^{i(kx − ωt)}:
    û(r) = AJ₀(βr) + (k/ω)(Eh/(ρa₀²))η̂ − (ν̂/(βa₀))AJ₁(βa₀),
    v̂(r) = −(ikA/β)J₁(βr) − (ik²/ω)(Eh η̂/(ρa₀²))(r/2) + (ikν̂/(βa₀))A(r/2)J₁(βa₀),  β = √(iω/ν).

    Book: §16.3, Eqs. (16.131), (16.133) (from (16.130) and (16.132)).
    Parameters: r [m];  A [m/s] (arbitrary amplitude);  eta η̂ [m];  k (complex) [1/m];  omega [rad/s];  nu [m²/s];
    E [Pa];  h [m];  a0 [m];  rho [kg/m³];  poisson ν̂.  Returns (u_hat, v_hat) complex [m/s].
    Validation: V1 residuals of (16.130) and (16.122) by finite differences. Label: analytic.
    """
    from scipy.special import jv

    r = np.asarray(r, dtype=float)
    beta = np.sqrt(1j * float(omega) / float(nu))
    J1a = jv(1, beta * a0)
    stiff = float(E) * float(h) / (float(rho) * float(a0) ** 2)
    u_hat = A * jv(0, beta * r) + k / omega * stiff * eta - complex(poisson) / (beta * a0) * A * J1a  # (16.131)
    v_hat = (-1j * k * A / beta * jv(1, beta * r) - 1j * k ** 2 / omega * stiff * eta * r / 2.0
             + 1j * k * complex(poisson) / (beta * a0) * A * r / 2.0 * J1a)  # (16.133)
    return u_hat, v_hat


def elastic_tube_determinant_numeric(Xk, alpha, poisson=0.5, *, rtol: float = 1e-11, s0: float = 1e-4) -> complex:
    """Determinant of the two wall conditions for the long-wave modes, computed **without Bessel closed forms** (ours):
    the regular solution φ of φ″ + φ′/s + β̂²φ = 0 and ∫₀ˢ sφ ds are integrated with complex ``solve_ivp`` (DOP853) on
    s = r/a₀ ∈ [s₀, 1] from the series at s₀; then û = Aφ + C, C from (16.130) (constant particular solution), v̂ from
    (16.122)/(16.132), ξ̂ from (16.129), and the two conditions (16.124), (16.125) give a homogeneous 2 × 2 system in
    (A, η̂).  Its determinant (rows normalised) vanishes exactly where (16.136) does.

    Book: §16.3, Eqs. (16.122)–(16.125), (16.128)–(16.133) (the route the book condenses "after considerable algebra").
    A dimensional set is chosen internally (ω = 2π rad/s, a₀ = 1 cm, ρ = 1000 kg/m³, c₀ = 5 m/s, ν from α); the root
    X is independent of it.
    Parameters
    ----------
    Xk : X = k²c₀²/ω² (complex);  alpha : Womersley number;  poisson : ν̂;  rtol : ODE tolerance;  s0 : series start.
    Returns
    -------
    det (complex, dimensionless after row normalisation).
    Cost: ≈ 30 ms per call.
    Assumptions: as :func:`elastic_tube_characteristic`.
    Validation: V1 independent-route cross-check: vanishes (to ≈ 1e-9 of its typical size) at the roots of
    :func:`elastic_tube_characteristic`; V3 order 6 in the series start s₀ (the truncated series at s₀ has error
    O(s₀⁶); observed order ≈ 6.0).
    Label: analytic, converged.
    """
    alpha = float(alpha)
    nu_ = complex(poisson)
    omega, a0, rho, c0 = 2.0 * np.pi, 0.01, 1000.0, 5.0
    nu = omega * a0 ** 2 / alpha ** 2
    mu = rho * nu
    Eh = 2.0 * rho * a0 * c0 ** 2
    k = np.sqrt(complex(Xk)) * omega / c0
    bh2 = 1j * alpha ** 2  # β̂² = iω a₀²/ν (e^{−iωt})

    def rhs(s, y):
        phi, dphi, I = y[0], y[1], y[2]
        return [dphi, -dphi / s - bh2 * phi, s * phi]

    y0 = [1.0 - bh2 * s0 ** 2 / 4.0 + bh2 ** 2 * s0 ** 4 / 64.0, -bh2 * s0 / 2.0 + bh2 ** 2 * s0 ** 3 / 16.0,
          s0 ** 2 / 2.0 - bh2 * s0 ** 4 / 16.0]
    sol = solve_ivp(rhs, (s0, 1.0), np.asarray(y0, dtype=complex), method="DOP853", rtol=rtol, atol=1e-14)
    phi1, dphi1, I1 = sol.y[:, -1]
    beta2 = 1j * omega / nu
    # û(r) = A φ(r/a₀) + C,  C = RHS/β²,  RHS = (ik/μ)(Eh/a₀²)η̂ + (ν̂/a₀)û′(a₀),  û′(a₀) = A φ′(1)/a₀
    dU_A = dphi1 / a0                                   # û′(a₀) per unit A
    C_eta = (1j * k / mu) * (Eh / a0 ** 2) / beta2      # C per unit η̂
    C_A = (nu_ / a0) * dU_A / beta2                     # C per unit A
    # (16.124): û(a₀) + iω ξ̂ = 0,  ξ̂ = (iν̂/(k a₀)) η̂ − ((1 − ν̂²)/(Eh k²)) μ û′(a₀)          (16.129)
    xi_eta = 1j * nu_ / (k * a0)
    xi_A = -(1.0 - nu_ ** 2) / (Eh * k ** 2) * mu * dU_A
    row1 = np.array([phi1 + C_A + 1j * omega * xi_A, C_eta + 1j * omega * xi_eta])
    # (16.125): v̂(a₀) + iω η̂ = 0,  v̂(a₀) = −(ik/a₀)∫₀^{a₀} r û dr = −(ik/a₀)[A a₀² I(1) + C a₀²/2]      (16.132)
    v_A = -(1j * k / a0) * (a0 ** 2 * I1 + C_A * a0 ** 2 / 2.0)
    v_eta = -(1j * k / a0) * (C_eta * a0 ** 2 / 2.0)
    row2 = np.array([v_A, v_eta + 1j * omega])
    row1 = row1 / np.linalg.norm(row1)
    row2 = row2 / np.linalg.norm(row2)
    return complex(row1[0] * row2[1] - row1[1] * row2[0])


# ======================================================================================================================
# §16.3 Murray's law and the optimal bifurcation (I14)
# ======================================================================================================================
def murray_cost(a, Q, L, mu, b_metabolic):
    """Power to run a vessel: pumping 8μLQ²/(πa⁴) plus upkeep b πa²L.

    Book: §16.3 "Blood Vessel Bifurcation …", Eq. (16.138) (Ê₁; the book's K is our ``b_metabolic`` [W/m³], trap T2).
    Parameters
    ----------
    a : radius [m];  Q : flow [m³/s];  L : length [m];  mu [Pa s];  b_metabolic : upkeep per unit volume [W/m³].
    Returns
    -------
    E1 [W].
    Assumptions: Poiseuille flow, cost ∝ volume.  Validation: V1 :func:`murray_terms` sum. Label: analytic.
    """
    return _S(8.0 * _F(mu) * _F(L) / (np.pi * _F(a) ** 4) * _F(Q) ** 2
              + _F(b_metabolic) * np.pi * _F(a) ** 2 * _F(L))  # Eq. (16.138)


def murray_terms(a, Q, L, mu, b_metabolic) -> dict:
    """The two terms of (16.138) and their total (E6 bars; ours).  Book: §16.3, Eq. (16.138).

    Returns dict pumping, upkeep, total [W].  At the optimum pumping = upkeep/2 (derivation D22).  Label: analytic.
    """
    a = _F(a)
    pump = 8.0 * _F(mu) * _F(L) / (np.pi * a ** 4) * _F(Q) ** 2
    up = _F(b_metabolic) * np.pi * a ** 2 * _F(L)
    return dict(pumping=_S(pump), upkeep=_S(up), total=_S(pump + up))


def murray_radius(Q, mu, b_metabolic, *, printed_16_140: bool = False):
    """Optimal radius a = [16μ/(π²b)]^{1/6} Q^{1/3} (Murray's law Q ∝ a³).

    Book: §16.3, Eqs. (16.139)–(16.140), printed a = [(16μ/π²)K]^{1/6}Q^{1/3}.
    DEVIATION: slip S15 — solving (16.139) gives a⁶ = 16μQ²/(π²K): K belongs in the **denominator** (also required
    by dimensions, K in W/m³).  ``printed_16_140=True`` returns the printed form for the test that must fail.
    Parameters
    ----------
    Q [m³/s];  mu [Pa s];  b_metabolic [W/m³];  printed_16_140.
    Returns
    -------
    a [m].
    Assumptions: as :func:`murray_cost`.
    Validation: V2 sympy ∂Ê₁/∂a = 0; V1 grid search of :func:`murray_cost`. Label: symbolic, analytic.
    """
    mu, b = _F(mu), _F(b_metabolic)
    if printed_16_140:
        return _S((16.0 * mu / np.pi ** 2 * b) ** (1.0 / 6.0) * np.cbrt(_F(Q)))  # printed (16.140)
    return _S((16.0 * mu / (np.pi ** 2 * b)) ** (1.0 / 6.0) * np.cbrt(_F(Q)))  # (16.140) corrected


def murray_min_cost(Q, L, mu, b_metabolic):
    """Minimum cost Ê₁,min = (3π/2) b L a² at the optimal radius.

    Book: §16.3, Eq. (16.141).  Parameters as :func:`murray_cost` without a.  Returns [W].
    Validation: V1 = :func:`murray_cost` at :func:`murray_radius`. Label: analytic.
    """
    a = _F(murray_radius(Q, mu, b_metabolic))
    return _S(1.5 * np.pi * _F(b_metabolic) * _F(L) * a ** 2)  # Eq. (16.141)


def murray_daughters(a0, flow_split):
    """Daughter radii at a Murray-optimal fork: a₁ = a₀s^{1/3}, a₂ = a₀(1 − s)^{1/3}, s = Q₁/Q₀ — so a₀³ = a₁³ + a₂³.

    Book: §16.3, Eqs. (16.143)–(16.144).
    Parameters: a0 parent radius [m];  flow_split s = Q₁/Q₀ ∈ [0, 1].  Returns (a1, a2) [m].
    Validation: V1 cube law. Label: analytic.
    """
    s = _F(flow_split)
    return _S(_F(a0) * np.cbrt(s)), _S(_F(a0) * np.cbrt(1.0 - s))  # Eq. (16.144)


def murray_angles(a0, a1, a2):
    """Optimal branching angles from the radii (law of cosines on the weighted triangle a₀², a₁², a₂²):
    cos θ = (a₀⁴ + a₁⁴ − a₂⁴)/(2a₀²a₁²), cos φ = (a₀⁴ − a₁⁴ + a₂⁴)/(2a₀²a₂²).

    Book: §16.3, Eq. (16.146) (also cos(θ + φ) = (a₀⁴ − a₁⁴ − a₂⁴)/(2a₁²a₂²)).  If a cosine falls outside [−1, 1] no
    interior optimum exists: NaN and ``valid=False`` are returned.
    Parameters: a0, a1, a2 [m].
    Returns dict theta, phi, total (θ + φ) [rad], valid (bool; scalar inputs) — arrays for array inputs.
    Validation: V3 :func:`bifurcation_point` reproduces them; symmetric fork θ = arccos 2^{−1/3}; the third relation
    of (16.146) holds for cube-law radii. Label: analytic.
    """
    a0, a1, a2 = (_F(v) for v in (a0, a1, a2))
    ct = (a0 ** 4 + a1 ** 4 - a2 ** 4) / (2.0 * a0 ** 2 * a1 ** 2)  # Eq. (16.146)
    cp = (a0 ** 4 - a1 ** 4 + a2 ** 4) / (2.0 * a0 ** 2 * a2 ** 2)
    bad = (np.abs(ct) > 1.0) | (np.abs(cp) > 1.0)
    with np.errstate(invalid="ignore"):
        th = np.where(bad, np.nan, np.arccos(np.clip(ct, -1.0, 1.0)))
        ph = np.where(bad, np.nan, np.arccos(np.clip(cp, -1.0, 1.0)))
    valid = ~bad
    return dict(theta=_S(th), phi=_S(ph), total=_S(th + ph), valid=bool(valid) if valid.ndim == 0 else valid)


def murray_theta(a0, a1):
    """Optimal angle of branch 1 with a₂ eliminated by the cube law:
    cos θ = [a₀⁴ + a₁⁴ − (a₀³ − a₁³)^{4/3}]/(2a₀²a₁²).

    Book: §16.3, Eq. (16.147).  Parameters a0, a1 [m] (a1 < a0).  Returns theta [rad].
    Validation: V1 = :func:`murray_angles` with a₂ = (a₀³ − a₁³)^{1/3}. Label: analytic.
    """
    a0, a1 = _F(a0), _F(a1)
    ct = (a0 ** 4 + a1 ** 4 - (a0 ** 3 - a1 ** 3) ** (4.0 / 3.0)) / (2.0 * a0 ** 2 * a1 ** 2)  # Eq. (16.147)
    return _S(np.arccos(np.clip(ct, -1.0, 1.0)))


def murray_symmetric_angle() -> float:
    """Branch angle of a symmetric Murray fork, θ = φ = arccos 2^{−1/3} ≈ 37.47° (ours, from (16.146) with
    a₁ = a₂ = 2^{−1/3}a₀).  Book: §16.3, Eqs. (16.144), (16.146).  Returns [rad].
    Validation: V1 = :func:`murray_angles` of the symmetric cube-law fork. Label: analytic."""
    return float(np.arccos(2.0 ** (-1.0 / 3.0)))


def bifurcation_cost(B, A, C, D, a0, a1, a2):
    """Weighted length Σ a_i²L_i of a planar fork with fixed ends A, C, D and branch point B — proportional to
    Ê₂,min = (3π/2)K(L₀a₀² + L₁a₁² + L₂a₂²).

    Book: §16.3, Eq. (16.142).  Parameters: B, A, C, D points (x, y) [m];  a0, a1, a2 radii [m].  Returns [m³].
    Label: analytic.
    """
    B, A, C, D = (np.asarray(v, dtype=float) for v in (B, A, C, D))
    return float(a0 ** 2 * np.linalg.norm(B - A) + a1 ** 2 * np.linalg.norm(C - B) + a2 ** 2 * np.linalg.norm(D - B))


def bifurcation_point(A, C, D, a0, a1, a2) -> dict:
    """Optimal branch point B minimising :func:`bifurcation_cost` (weighted Fermat point; ours: Nelder–Mead then
    Powell polish), with the angles it produces.

    Book: §16.3, Eqs. (16.145)–(16.146) (the variational argument "making displacements of B, one at a time"); the
    numerical minimisation is our independent route.
    Parameters: A, C, D points [m];  a0, a1, a2 [m].
    Returns dict: B (2,), cost, theta, phi [rad] (θ between the continuation of AB and BC, φ for BD).
    Validation: V1 independent-route cross-check: angles equal the closed form :func:`murray_angles` to 1e-6 rad (no
    refinement order asserted). Label: analytic.
    """
    from scipy.optimize import minimize

    A, C, D = (np.asarray(v, dtype=float) for v in (A, C, D))
    f = lambda B: bifurcation_cost(B, A, C, D, a0, a1, a2)  # noqa: E731
    x0 = (A + C + D) / 3.0
    scale = max(np.linalg.norm(C - A), np.linalg.norm(D - A))
    r1 = minimize(f, x0, method="Nelder-Mead", options=dict(xatol=1e-12 * scale, fatol=1e-16, maxiter=20000))
    r2 = minimize(f, r1.x, method="Nelder-Mead", options=dict(xatol=1e-14 * scale, fatol=1e-18, maxiter=20000,
                                                              initial_simplex=None))
    B = r2.x
    e0 = (B - A) / np.linalg.norm(B - A)
    e1 = (C - B) / np.linalg.norm(C - B)
    e2 = (D - B) / np.linalg.norm(D - B)
    th = float(np.arccos(np.clip(e0 @ e1, -1, 1)))
    ph = float(np.arccos(np.clip(e0 @ e2, -1, 1)))
    return dict(B=B, cost=float(r2.fun), theta=th, phi=ph)


# ======================================================================================================================
# §16.3 wave reflection at junctions (I15)
# ======================================================================================================================
def tube_admittance(A, rho, c):
    """Characteristic admittance Y = A/(ρc) (flow per unit pressure of a forward wave).

    Book: §16.3 "Reflection of Waves at Arterial Junctions", Eqs. (16.150)–(16.151) (from p̂ = ρcû (16.68)).
    Parameters: A [m²];  rho [kg/m³];  c [m/s].  Returns Y [m³/(Pa s)].  Label: analytic.
    """
    return _S(_F(A) / (_F(rho) * _F(c)))  # Eq. (16.150)


def tube_impedance(A, rho, c):
    """Characteristic impedance Z = ρc/A = 1/Y (ch15's acoustic ρc divided by the area).

    Book: §16.3, Eq. (16.151).  Parameters: A [m²];  rho;  c.  Returns Z [Pa s/m³].  Label: analytic.
    """
    return _S(_F(AC.acoustic_impedance(rho, c)) / _F(A))  # Eq. (16.151)


def junction_coefficients(Y_parent, Y_daughters):
    """Pressure reflection and transmission coefficients of a junction:
    R = (Y₁ − ΣY_j)/(Y₁ + ΣY_j), T = 2Y₁/(Y₁ + ΣY_j) = 1 + R.

    Book: §16.3, Eqs. (16.148), (16.152)–(16.154) (derivation D25).  The flux balance is written with magnitudes; the
    book's (16.149) with (16.158) counts the sign of Q_R twice (slip S16, loose).  ΣY_j = Y₁ matched (R = 0), < Y₁
    closed-end-like (R > 0), > Y₁ open-end-like (R < 0) (trap T18).
    Parameters: Y_parent [m³/(Pa s)];  Y_daughters sequence (or scalar total).  Returns (R, T) [–].
    Validation: V1 np.linalg.solve of (16.148), (16.152); V4 :func:`junction_energy`. Label: analytic.
    """
    Y1 = float(Y_parent)
    Ys = float(np.sum(_F(Y_daughters)))
    R = (Y1 - Ys) / (Y1 + Ys)  # Eq. (16.153)
    T = 2.0 * Y1 / (Y1 + Ys)
    return R, T


def junction_energy(Y_parent, Y_daughters):
    """Energy shares at a junction: reflected R², transmitted (ΣY_j/Y₁)T², and their sum (= 1 identically; ours).

    Book: §16.3, Eqs. (16.161)–(16.163) (energy flux pQ = p²/Z); the conservation identity is ours (derivation D26).
    Parameters as :func:`junction_coefficients`.  Returns dict reflected (R²), transmitted (ΣY_jT²/Y₁), total [–].
    Validation: V4 total = 1 to 1e-14 for random admittances. Label: conserved.
    """
    R, T = junction_coefficients(Y_parent, Y_daughters)
    trans = float(np.sum(_F(Y_daughters))) / float(Y_parent) * T ** 2  # Eq. (16.163)
    return dict(reflected=R ** 2, transmitted=trans, total=R ** 2 + trans)  # Eq. (16.162)


def junction_from_geometry(a_parent, a_daughters, c_parent, c_daughters, rho) -> dict:
    """Admittances, R, T and energy shares of a fork from radii and wave speeds (E7 presets; ours).

    Book: §16.3, (16.150)–(16.163).  Parameters: a_parent [m];  a_daughters sequence [m];  c_parent [m/s];
    c_daughters sequence or scalar [m/s];  rho [kg/m³].
    Returns dict Y1, Y (array), ratio (ΣY/Y₁), R, T, reflected (R²), transmitted (ΣY_jT²/Y₁).
    Label: analytic.
    """
    ad = np.atleast_1d(_F(a_daughters))
    cd = np.broadcast_to(_F(c_daughters), ad.shape)
    Y1 = float(tube_admittance(np.pi * float(a_parent) ** 2, rho, c_parent))
    Y = _F(tube_admittance(np.pi * ad ** 2, rho, cd))
    R, T = junction_coefficients(Y1, Y)
    en = junction_energy(Y1, Y)
    return dict(Y1=Y1, Y=Y, ratio=float(Y.sum() / Y1), R=R, T=T, reflected=en["reflected"],
                transmitted=en["transmitted"])


def junction_waveforms(x, t, f: Callable, c1, R, *, Y1: float = 1.0, P_I: float = 1.0) -> dict:
    """Pressure and flow in the parent tube (x ≤ 0) when a wave f meets a junction at x = 0:
    p = P_I[f(t − x/c₁) + R f(t + x/c₁)],  Q = Y₁P_I[f(t − x/c₁) − R f(t + x/c₁)].

    Book: §16.3, Eqs. (16.155)–(16.160) (the reflected shape is written f because R does not depend on frequency here —
    slip S17, loose).  p and Q no longer share a shape once R ≠ 0.
    Parameters: x [m] (≤ 0);  t [s];  f : time shape (max 1);  c1 [m/s];  R [–];  Y1 [m³/(Pa s)];  P_I [Pa].
    Returns dict p [Pa], Q [m³/s].
    Label: analytic.
    """
    x, t = _F(x), _F(t)
    inc = _F(f(t - x / float(c1)))
    ref = _F(f(t + x / float(c1)))
    return dict(p=_S(P_I * (inc + R * ref)), Q=_S(Y1 * P_I * (inc - R * ref)))  # Eqs. (16.159)–(16.160)


def network_pulse_run(L_parent, L_daughters, c_parent, c_daughters, Y_parent, Y_daughters, *, width: float,
                      nx: int = 400, t_end: float | None = None, amplitude: float = 1.0, t_center: float | None = None,
                      x_probe: float | None = None, n_frames: int = 120, fast: bool = False) -> dict:
    """A pressure pulse running down a parent tube into a junction with any number of daughters, computed on a grid
    (ours): the linear 1-D tube equations in each segment, solved exactly by transporting the characteristic variables
    w± = p ± ZQ (Z = 1/Y) one cell per step (Courant number 1 in every segment, Δx_j = c_jΔt), with p continuous and
    Q conserved at the node.

    Book: §16.3, (16.148)–(16.153) (junction conditions) and (16.65)–(16.68) (the d'Alembert pieces); the grid method is
    ours (DEVIATION: characteristic transport at C = 1 instead of the leapfrog of the analysis plan — both are exact at
    C = 1; this one handles the junction node directly).  Junction: p_J = (Y₁w₁⁺ + ΣY_jw_j⁻)/(Y₁ + ΣY_j), outgoing
    w₁⁻ = 2p_J − w₁⁺, w_j⁺ = 2p_J − w_j⁻.  Ends: the incident pulse p_I(t) = amplitude·exp(−((t − t_c)/width)²) enters at
    x = −L_parent (w₁⁺ = 2p_I); all far ends are non-reflecting.
    Parameters
    ----------
    L_parent [m];  L_daughters sequence [m];  c_parent [m/s];  c_daughters sequence or scalar [m/s];
    Y_parent [m³/(Pa s)];  Y_daughters sequence;  width : time width of the Gaussian pulse [s];  nx : cells in the
    parent;  t_end [s] (default: until the reflection has passed the probe);  amplitude [Pa];  t_center [s] (default
    3·width);  x_probe [m] (default −L_parent/2);  n_frames : stored frames;  fast : nx = 160.
    Returns
    -------
    dict: x = dict(parent, daughters list) [m]; t (n_frames,) [s]; p_parent (n_frames, nx + 1) [Pa]; p_daughters list
    of (n_frames, n_j + 1) [Pa]; probe_t, probe_p, probe_Q (every step) [s, Pa, m³/s]; R_measured, T_measured (peak
    ratios: reflected at the probe, transmitted in daughter 1, over incident), R_theory, T_theory (16.153), dt.
    Cost: ≈ 30 ms at nx = 400.
    Assumptions: linear, inviscid, long wave, no storage at the junction.
    Validation: V1 independent-route cross-check: R_measured, T_measured equal (16.153) to ≤ 1e-3 (no refinement
    order asserted). Label: analytic.
    """
    if fast:
        nx = 160
    Ld = np.atleast_1d(_F(L_daughters))
    nd = Ld.size
    cd = np.broadcast_to(_F(c_daughters), Ld.shape).astype(float)
    Yd = np.broadcast_to(_F(Y_daughters), Ld.shape).astype(float)
    L1, c1, Y1 = float(L_parent), float(c_parent), float(Y_parent)
    width = float(width)
    t_center = 3.0 * width if t_center is None else float(t_center)
    x_probe = -0.5 * L1 if x_probe is None else float(x_probe)
    if t_end is None:
        t_end = t_center + (L1 + abs(x_probe)) / c1 + 4.0 * width
    dt = L1 / (c1 * int(nx))
    nseg = [int(nx)] + [max(2, int(round(Ld[j] / (cd[j] * dt)))) for j in range(nd)]
    Ys = np.concatenate([[Y1], Yd])
    Zs = 1.0 / Ys
    wp = [np.zeros(n + 1) for n in nseg]
    wm = [np.zeros(n + 1) for n in nseg]
    xp = -L1 + np.arange(nseg[0] + 1) * c1 * dt
    xd = [np.arange(nseg[j + 1] + 1) * cd[j] * dt for j in range(nd)]
    i_probe = int(np.argmin(np.abs(xp - x_probe)))
    pinc = lambda tt: amplitude * np.exp(-((tt - t_center) / width) ** 2)  # noqa: E731
    nsteps = int(np.ceil(t_end / dt))
    t_arr = np.arange(nsteps + 1) * dt
    p_pr = np.zeros(nsteps + 1)
    Q_pr = np.zeros(nsteps + 1)
    p_d1 = np.zeros(nsteps + 1)
    frame_steps = np.unique(np.round(np.linspace(0, nsteps, int(n_frames))).astype(int))
    fr_t, fr_p, fr_d = [], [], [[] for _ in range(nd)]
    wp[0][0] = 2.0 * pinc(0.0)
    for k in range(0, nsteps + 1):
        if k > 0:
            tk = k * dt
            for j in range(nd + 1):  # transport one cell (C = 1)
                wp[j][1:] = wp[j][:-1].copy()
                wm[j][:-1] = wm[j][1:].copy()
            wp[0][0] = 2.0 * pinc(tk)  # incident wave enters the parent
            for j in range(1, nd + 1):
                wm[j][-1] = 0.0  # non-reflecting far ends of the daughters
            pJ = (Ys[0] * wp[0][-1] + sum(Ys[j] * wm[j][0] for j in range(1, nd + 1))) / Ys.sum()  # (16.148), (16.152)
            wm[0][-1] = 2.0 * pJ - wp[0][-1]
            for j in range(1, nd + 1):
                wp[j][0] = 2.0 * pJ - wm[j][0]
        p_pr[k] = 0.5 * (wp[0][i_probe] + wm[0][i_probe])
        Q_pr[k] = 0.5 * (wp[0][i_probe] - wm[0][i_probe]) / Zs[0]
        p_d1[k] = 0.5 * (wp[1][nseg[1] // 2] + wm[1][nseg[1] // 2])
        if k in frame_steps:
            fr_t.append(k * dt)
            fr_p.append(0.5 * (wp[0] + wm[0]))
            for j in range(nd):
                fr_d[j].append(0.5 * (wp[j + 1] + wm[j + 1]))
    # reflected wave at the probe = record − the known incident wave p_I(t − (x_probe + L)/c₁)
    x_pr = float(xp[i_probe])
    t_ret = t_arr - (x_pr + L1) / c1
    p_incident = np.where(t_ret >= -0.5 * dt, pinc(t_ret), 0.0)  # nothing is injected before t = 0
    p_reflected = p_pr - p_incident
    P_inc = amplitude
    P_ref = p_reflected[np.argmax(np.abs(p_reflected))]
    P_tr = p_d1[np.argmax(np.abs(p_d1))]
    R_th, T_th = junction_coefficients(Y1, Yd)
    return dict(x=dict(parent=xp, daughters=xd), t=np.array(fr_t), p_parent=np.array(fr_p),
                p_daughters=[np.array(v) for v in fr_d], probe_t=t_arr, probe_p=p_pr, probe_Q=Q_pr,
                R_measured=float(P_ref / P_inc), T_measured=float(P_tr / P_inc), R_theory=R_th, T_theory=T_th,
                dt=dt, x_probe=float(xp[i_probe]))


# ======================================================================================================================
# §16.3 flow in a curved tube: Dean (I16)
# ======================================================================================================================
DEAN_BOOK = dict(c2=0.0306, c4=0.0120, lam4=0.0110)  #: the printed coefficients of (16.178), (16.182)


def curvature_ratio(a, R_curv):
    """Curvature ratio δ = a/R (tube radius over radius of curvature of the centre line).

    Book: §16.3 "Flow in a Rigid-Walled Curved Tube", Eq. (16.164).  Parameters: a [m];  R_curv [m].  Returns δ [–].
    Label: analytic.
    """
    return _S(_F(a) / _F(R_curv))  # Eq. (16.164)


def dean_number(delta, U_mean, a, nu, *, convention: str):
    """Dean number in one of the chapter's conventions (trap T8):
    κ = √δ·2aW̄₀/ν (16.171), K = 2κ² (16.178), D = 4√2κ (16.179)–(16.180); "De" = Re_d√δ = κ.

    Book: §16.3, Eqs. (16.171), (16.178), (16.179)–(16.180).  (16.171) writes κ = √δ 2Re, consistent only with
    Re = aW̄₀/ν (slip S18, loose); (16.180) prints Ĝ = −8μW̄/a² with the sign of a gradient (slip S19, loose).
    Parameters
    ----------
    delta : a/R [–];  U_mean : mean axial speed W̄₀ [m/s];  a : radius [m];  nu [m²/s];
    convention : "kappa" | "De" | "K" | "D" (required keyword).
    Returns
    -------
    Dean number [–].
    Assumptions: δ ≪ 1.  Validation: V1 conversions round-trip (:func:`dean_convert`). Label: analytic.
    """
    kappa = np.sqrt(_F(delta)) * 2.0 * _F(a) * _F(U_mean) / _F(nu)  # Eq. (16.171)
    return _S(dean_convert(kappa, from_="kappa", to=convention))


def dean_convert(value, *, from_: str, to: str):
    """Convert between the Dean-number conventions κ (= De), K = 2κ², D = 4√2 κ (so K = D²/16).

    Book: §16.3, (16.171), (16.178) (K = (2a/R)(W_max a/ν)² = 2κ², W_max = 2W̄), (16.179)–(16.181).
    Parameters: value;  from_, to : "kappa" | "De" | "K" | "D".  Returns the converted value.
    Validation: V1 round trips; K = D²/16. Label: analytic.
    """
    v = _F(value)
    names = ("kappa", "De", "K", "D")
    if from_ not in names or to not in names:
        raise ValueError(f"conventions are {names}")
    kappa = {"kappa": v, "De": v, "K": np.sqrt(v / 2.0), "D": v / (4.0 * np.sqrt(2.0))}[from_]
    return _S({"kappa": kappa, "De": kappa, "K": 2.0 * kappa ** 2, "D": 4.0 * np.sqrt(2.0) * kappa}[to])


def dean_flux_ratio(K=None, *, D=None, kappa=None, exact: bool = False, exact_c4: bool = False):
    """Flux in a gently curved pipe relative to a straight one at the same pressure gradient:
    Q/(πa²W̄) = 1 − 0.0306(K/576)² + 0.0120(K/576)⁴ (= 1 − 0.0306(D/96)⁴ + 0.0120(D/96)⁸).

    Book: §16.3, Eqs. (16.178), (16.181) (Dean 1928).  ``exact=True`` replaces the printed 0.0306 by our exact
    1541/50400 = 0.0305754 (:func:`dean_series_sympy`, order 2) and keeps the printed 0.0120 (design C.12);
    ``exact_c4=True`` also replaces 0.0120 by our order-4 value 6471982981/542442700800 = 0.0119312 — 0.6 % below
    the printed coefficient (equations (16.173)–(16.176) solved exactly, with the advection sign of (16.174)
    corrected — slip S34; an independent spectral Navier–Stokes solve gives 0.0119311).
    # DEVIATION: c₄ from the consistent sign of (16.174), not the printed one (which would give 0.0121130).
    Parameters
    ----------
    K : Dean's K = 2κ² [–] (or give D or kappa);  exact, exact_c4 : use the derived coefficients.
    Returns
    -------
    ratio [–].
    Assumptions: δ ≪ 1, small K (series), fully developed, steady.
    Validation: V2 our sympy re-derivation reproduces 0.0306 (1541/50400); V3 c₄ against a spectral Navier–Stokes
    solve. Label: symbolic.
    """
    if K is None:
        K = dean_convert(D, from_="D", to="K") if D is not None else dean_convert(kappa, from_="kappa", to="K")
    q = _F(K) / 576.0
    c2 = 1541.0 / 50400.0 if exact else DEAN_BOOK["c2"]
    c4 = float(DEAN_C4_DERIVED) if exact_c4 else DEAN_BOOK["c4"]
    return _S(1.0 - c2 * q ** 2 + c4 * q ** 4)  # Eq. (16.178)


#: Our exact second flux coefficient of (16.178) from :func:`dean_series_sympy` at order 4 (≈ 0.0119312), with the
#: sign of the advection term of (16.174) corrected (slip S34; printed sign gives 821328579941/67805337600000).
DEAN_C4_DERIVED = Fraction(6471982981, 542442700800)


def curved_resistance_ratio(K):
    """Resistance of the curved pipe over the straight one, λ = (Q_c/Q_s)⁻¹ = 1 + 0.0306(K/576)² − 0.0110(K/576)⁴.

    Book: §16.3, Eq. (16.182) (the text after it swaps "straight" and "curved" — slip S20, loose).  The K⁴ coefficient
    follows from inverting (16.178): 0.0306² − 0.0120 = −0.01106.
    Parameters: K [–].  Returns λ [–].
    Validation: V1 = 1/(16.178) to O(K⁴). Label: analytic.
    """
    q = _F(K) / 576.0
    return _S(1.0 + DEAN_BOOK["c2"] * q ** 2 - DEAN_BOOK["lam4"] * q ** 4)  # Eq. (16.182)


def dean_secondary_streamfunction(r, alpha_t):
    """Leading secondary-flow stream function ψ₀ = r(1 − r²)²(4 − r²) sin α/144 (ours, from (16.174) at O(1)).

    Book: §16.3, Eqs. (16.172)–(16.177) (the book states the series; ψ₀ is our solution of ∇⁴ψ₀ = 8r(1 − r²) sin α with
    ψ₀ = ∂ψ₀/∂r = 0 at r = 1, where w₀ = 2(1 − r²)).  Physical ψ ≈ κψ₀; u = (1/r)∂ψ/∂α, v = −∂ψ/∂r (16.172).
    Coordinates: r = r′/a, α_t the angle from the outer bend (α_t = 0 outer — the centrifugal term −w²cos α of (16.167)).
    Parameters: r [–] (0…1);  alpha_t [rad].  Returns ψ₀ [–].
    Validation: V2 :func:`dean_series_sympy` (exact rationals). Label: symbolic.
    """
    r, al = _F(r), _F(alpha_t)
    return _S(r * (1.0 - r ** 2) ** 2 * (4.0 - r ** 2) * np.sin(al) / 144.0)


def dean_axial_correction(r, alpha_t):
    """First curvature correction of the axial velocity, w₁ = +r(1 − r²)(19 − 21r² + 9r⁴ − r⁶)cos α/5760 (ours, from
    (16.173) at O(κ²): ∇²w₁ = (1/2r)(∂ψ₀/∂α ∂w₀/∂r), w₁(1) = 0).

    Book: §16.3, (16.173), (16.177).  The **plus** sign (positive at α_t = 0, the outer bend) is what the exact series
    gives — the fast core is carried outward; the analysis note had a minus (design C.8 "header correction").
    Parameters: r [–];  alpha_t [rad].  Returns w₁ [–] (w = w₀ + κ²w₁ + …).
    Validation: V2 :func:`dean_series_sympy` (exact rationals). Label: symbolic.
    """
    r, al = _F(r), _F(alpha_t)
    return _S(_W1_COEF * r * (1.0 - r ** 2) * (19.0 - 21.0 * r ** 2 + 9.0 * r ** 4 - r ** 6) * np.cos(al))


#: sign and size of w₁ — checked against the exact series in :func:`dean_series_sympy` (see the test)
_W1_COEF = 1.0 / 5760.0


def dean_axial_velocity(r, alpha_t, K):
    """Axial velocity to O(κ²) in a gently curved pipe, w = w₀ + κ²w₁ with w₀ = 2(1 − r²), κ² = K/2 (ours; E8, F6).

    Book: §16.3, (16.177).  Scaled by the straight-pipe mean speed at the same gradient.
    Parameters: r [–];  alpha_t [rad];  K [–].  Returns w [–].  Label: symbolic.
    """
    return _S(2.0 * (1.0 - _F(r) ** 2) + _F(K) / 2.0 * _F(dean_axial_correction(r, alpha_t)))


def dean_secondary_velocity(r, alpha_t, K):
    """Leading secondary velocity (u, v) = κ((1/r)∂ψ₀/∂α, −∂ψ₀/∂r), κ = √(K/2) (radial, angular; ours).

    Book: §16.3, Eq. (16.172) with ψ ≈ κψ₀ (16.177).  u > 0 on α_t = 0: the core moves toward the outer bend.
    Parameters: r [–] (> 0 for u);  alpha_t [rad];  K [–].  Returns (u, v) [–] (scaled like w).  Label: symbolic.
    """
    r, al = _F(r), _F(alpha_t)
    kappa = np.sqrt(_F(K) / 2.0)
    f = (1.0 - r ** 2) ** 2 * (4.0 - r ** 2) / 144.0          # ψ₀ = r f(r) sin α
    dpsi_dr = (f + r * (-4.0 * r * (1.0 - r ** 2) * (4.0 - r ** 2) - 2.0 * r * (1.0 - r ** 2) ** 2) / 144.0) * np.sin(al)
    u = kappa * f * np.cos(al)                                 # (1/r)∂ψ₀/∂α = f cos α
    v = -kappa * dpsi_dr
    return _S(u), _S(v)


def dean_series_reliability(K, *, threshold: float = 0.1) -> dict:
    """Our own reliability criterion for the small-K series (16.178): the ratio of the K⁴ term to the K² term,
    |0.0120(K/576)⁴|/|0.0306(K/576)²|, against ``threshold`` (never the book's quoted range).

    Book: §16.3, (16.178).  Parameters: K [–];  threshold [–].  Returns dict(ratio, reliable).  Label: analytic.
    """
    q = float(K) / 576.0
    ratio = DEAN_BOOK["c4"] * q ** 2 / DEAN_BOOK["c2"] if q > 0 else 0.0
    return dict(ratio=ratio, reliable=bool(ratio < threshold))


class _FP:
    """Fourier–polynomial field Σ_m [A_m(r) cos mα + B_m(r) sin mα] with exact rational (Laurent) polynomials in r."""

    def __init__(self, t=None):
        self.t = {} if t is None else {k: dict(v) for k, v in t.items()}

    def add(self, m, kind, p, c):
        if c == 0:
            return
        if m < 0:
            m = -m
            if kind == "s":
                c = -c
        if m == 0 and kind == "s":
            return
        d = self.t.setdefault((m, kind), {})
        d[p] = d.get(p, Fraction(0)) + c
        if d[p] == 0:
            del d[p]
        if not d:
            del self.t[(m, kind)]

    def items(self):
        for (m, kd), d in self.t.items():
            for p, c in d.items():
                yield m, kd, p, c

    def __add__(self, o):
        out = _FP(self.t)
        for m, kd, p, c in o.items():
            out.add(m, kd, p, c)
        return out

    def __sub__(self, o):
        return self + o.scale(-1)

    def scale(self, s):
        out = _FP()
        for m, kd, p, c in self.items():
            out.add(m, kd, p, c * Fraction(s))
        return out

    def __mul__(self, o):
        out = _FP()
        for m1, k1, p1, c1 in self.items():
            for m2, k2, p2, c2 in o.items():
                c, p = c1 * c2 / 2, p1 + p2
                if k1 == "c" and k2 == "c":
                    out.add(m1 - m2, "c", p, c)
                    out.add(m1 + m2, "c", p, c)
                elif k1 == "s" and k2 == "s":
                    out.add(m1 - m2, "c", p, c)
                    out.add(m1 + m2, "c", p, -c)
                elif k1 == "s" and k2 == "c":
                    out.add(m1 + m2, "s", p, c)
                    out.add(m1 - m2, "s", p, c)
                else:  # cos m1 · sin m2
                    out.add(m1 + m2, "s", p, c)
                    out.add(m1 - m2, "s", p, -c)
        return out

    def dr(self):
        out = _FP()
        for m, kd, p, c in self.items():
            out.add(m, kd, p - 1, c * p)
        return out

    def da(self):
        out = _FP()
        for m, kd, p, c in self.items():
            if kd == "c":
                out.add(m, "s", p, -c * m)
            else:
                out.add(m, "c", p, c * m)
        return out

    def rpow(self, k):
        out = _FP()
        for m, kd, p, c in self.items():
            out.add(m, kd, p + k, c)
        return out

    def lap(self):
        out = _FP()
        for m, kd, p, c in self.items():
            out.add(m, kd, p - 2, c * (p * p - m * m))
        return out

    def solve_poisson(self):
        """f with ∇²f = self, f(1) = 0, regular at r = 0."""
        out = _FP()
        for (m, kd), d in self.t.items():
            val = Fraction(0)
            for p, c in d.items():
                q = p + 2
                if q * q == m * m or q < m:
                    raise ValueError("non-regular particular solution")
                cc = c / (q * q - m * m)
                out.add(m, kd, q, cc)
                val += cc
            out.add(m, kd, m, -val)
        return out

    def solve_biharmonic(self):
        """f with ∇⁴f = self, f(1) = f′(1) = 0, regular at r = 0."""
        out = _FP()
        for (m, kd), d in self.t.items():
            g = {}
            for p, c in d.items():
                q = p + 2
                if q * q == m * m or q < m:
                    raise ValueError("non-regular particular solution")
                g[q] = c / (q * q - m * m)
            fp = {}
            for p, c in g.items():
                q = p + 2
                if q * q == m * m or q < m:
                    raise ValueError("non-regular particular solution")
                fp[q] = fp.get(q, Fraction(0)) + c / (q * q - m * m)
            f1 = sum(fp.values(), Fraction(0))
            df1 = sum((c * p for p, c in fp.items()), Fraction(0))
            c2 = (m * f1 - df1) / 2
            c1 = -f1 - c2
            for p, c in fp.items():
                out.add(m, kd, p, c)
            out.add(m, kd, m, c1)
            out.add(m, kd, m + 2, c2)
        return out

    def disk_mean(self):
        """(1/π)∫₀¹∫₀^{2π} f r dr dα = 2∫₀¹ A₀(r) r dr."""
        d = self.t.get((0, "c"), {})
        return sum((2 * c / (p + 2) for p, c in d.items()), Fraction(0))

    def to_sympy(self):
        import sympy as sp

        r, al = sp.symbols("r alpha", real=True)
        e = 0
        for m, kd, p, c in self.items():
            e += sp.Rational(c.numerator, c.denominator) * r ** p * (sp.cos(m * al) if kd == "c" else sp.sin(m * al))
        return sp.expand(e)


@lru_cache(maxsize=4)
def dean_series_sympy(order: int = 2, printed_16_174: bool = False) -> dict:
    """Dean's small-κ series for fully developed flow in a gently curved pipe, solved exactly order by order (ours).

    Book: §16.3, Eqs. (16.172)–(16.177): with p₀ = −Gs and the scaling (κ/2)G = 8 (so w₀ = 2(1 − r²), mean 1),
    w = Σκ^{2n}w_n, ψ = κΣκ^{2n}ψ_n,
      ∇²w_n = −8δ_{n0} + (1/2r)Σ_{i+j=n−1}(∂ψ_i/∂α ∂w_j/∂r − ∂ψ_i/∂r ∂w_j/∂α)                         (16.173)
      2∇⁴ψ_n = −(1/r)Σ_{i+j=n−1}(∂ψ_i/∂r ∂/∂α − ∂ψ_i/∂α ∂/∂r)∇²ψ_j − 2Σ_{i+j=n} w_i(sin α ∂w_j/∂r + (cos α/r)∂w_j/∂α)
                                                                                                      (16.174)
    with ψ = ∂ψ/∂r = w = 0 at r = 1 (16.176).
    DEVIATION (slip S34): the book prints (16.174) as (2/κ)∇⁴ψ − (1/r)(…)∇²ψ = −2w(…); eliminating p₁ from
    (16.167)–(16.168) with (16.172) gives (2/κ)∇⁴ψ + (1/r)(…)∇²ψ = −2w(…), which is what is coded (the − sign above
    after moving the term to the right).  ``printed_16_174=True`` reproduces the printed sign.  Only the order-κ⁸
    (K⁴) flux coefficient is affected: ψ₀, w₁ and c₂ = 1541/50400 are identical in both versions, while
    c₄ = 6471982981/542442700800 ≈ 0.0119312 (consistent sign; an independent spectral Navier–Stokes solve gives
    0.0119311) versus 821328579941/67805337600000 ≈ 0.0121130 with the printed sign.  The book's 0.0120 lies between.  Each order is a polynomial Poisson or biharmonic problem per Fourier
    mode, solved by undetermined coefficients in exact rational arithmetic (``fractions``).  The flux ratio at fixed
    gradient is mean(w) = 1 + κ⁴m₂ + κ⁸m₄; with K = 2κ² the coefficients of (K/576)² and (K/576)⁴ are
    c₂ = −576²m₂/4 and c₄ = 576⁴m₄/16.
    Parameters
    ----------
    order : highest n of w_n (2 gives c₂; 4 gives c₄ as well).
    printed_16_174 : if True, use the printed (inconsistent) sign of the advection term of (16.174) (slip S34).
    Returns
    -------
    dict: w (list of sympy expressions w₀…w_order), psi (ψ₀…ψ_{order−1}), w0, psi0, w1, means (Fractions),
    c2 (Fraction, = 1541/50400), flux_coefficient (sympy Rational(1541, 50400)), c4 and flux_coefficient_4
    (order ≥ 4: 6471982981/542442700800 ≈ 0.0119312, 0.6 % below the printed 0.0120; with printed_16_174=True
    821328579941/67805337600000 ≈ 0.0121130), else None.
    Cost: < 1 s for order 2, a few seconds for order 4 (cached).
    Assumptions: δ ≪ 1, steady, fully developed, small κ.
    Validation: V2 each w_n, ψ_n satisfies its equation and boundary conditions (sympy); c₂ reproduces 0.0306;
    V3 c₄ matches an independent spectral Navier–Stokes solve (0.0119311).
    Label: symbolic.
    """
    order = int(order)
    printed_16_174 = bool(printed_16_174)
    one = _FP()
    one.add(0, "c", 0, Fraction(1))
    sin_a, cos_a = _FP(), _FP()
    sin_a.add(1, "s", 0, Fraction(1))
    cos_a.add(1, "c", 0, Fraction(1))
    w: list[_FP] = []
    psi: list[_FP] = []
    lap_psi: list[_FP] = []
    for n in range(order + 1):
        src = _FP()
        if n == 0:
            src = one.scale(-8)
        for i in range(n):
            j = n - 1 - i
            jac = psi[i].da() * w[j].dr() - psi[i].dr() * w[j].da()
            src = src + jac.rpow(-1).scale(Fraction(1, 2))
        w.append(src.solve_poisson())
        if n == order:
            break
        s2 = _FP()
        for i in range(n):
            j = n - 1 - i
            term = psi[i].dr() * lap_psi[j].da() - psi[i].da() * lap_psi[j].dr()
            # DEVIATION: sign of the psi-advection term of (16.174) — eliminating p1 from (16.167)-(16.168) with
            # (16.172) puts it on the left with a + sign, i.e. on the right with a − sign (slip S34); the printed
            # (16.174) has − on the left (+ on the right), reproduced with printed_16_174=True.
            s2 = (s2 + term.rpow(-1)) if printed_16_174 else (s2 - term.rpow(-1))  # Eq. (16.174), advection term
        for i in range(n + 1):
            j = n - i
            term = w[i] * (sin_a * w[j].dr() + cos_a * w[j].da().rpow(-1))
            s2 = s2 - term.scale(2)
        p_n = s2.scale(Fraction(1, 2)).solve_biharmonic()
        psi.append(p_n)
        lap_psi.append(p_n.lap())
    import sympy as sp

    means = [wn.disk_mean() for wn in w]
    c2 = -means[2] * 576 ** 2 / 4 if order >= 2 else None
    c4 = means[4] * Fraction(576) ** 4 / 16 if order >= 4 else None
    W = [x.to_sympy() for x in w]
    PS = [x.to_sympy() for x in psi]
    return dict(w=W, psi=PS, w0=W[0], psi0=PS[0] if PS else None, w1=W[1] if len(W) > 1 else None, means=means,
                c2=c2, c4=c4, flux_coefficient=sp.Rational(c2.numerator, c2.denominator) if c2 is not None else None,
                flux_coefficient_4=sp.Rational(c4.numerator, c4.denominator) if c4 is not None else None,
                _fields=(w, psi))


# ======================================================================================================================
# §16.3 collapsible tubes: Shapiro's 1-D steady theory (I17)
# ======================================================================================================================
def bending_stiffness_parameter(E, h, a0, poisson, *, coeff: float):
    """Shapiro's stiffness parameter K_p ∝ EI, I = (h/a₀)³/(1 − ν̂²): K_p = coeff·E(h/a₀)³/(1 − ν̂²).

    Book: §16.3, Eqs. (16.186)–(16.187).  The proportionality constant is not given in the book — a required keyword.
    Parameters: E [Pa];  h [m];  a0 [m];  poisson [–];  coeff [–] (required).  Returns K_p [Pa].  Label: analytic.
    """
    return _S(float(coeff) * _F(E) * (_F(h) / _F(a0)) ** 3 / (1.0 - _F(poisson) ** 2))  # Eq. (16.187)


ALPHA_MAX_COLLAPSE = 1.0  # upper end of the fitted (collapse) branch of (16.188), Fig. 16.19
_ALPHA_ROUNDOFF = 1e-12   # relative round-off allowance above α = 1 (α_in = alpha_from_tube_law(0) is exactly 1)


def _collapse_guard(al, name: str = "area_ratio", stacklevel: int = 3):
    """Range guard of Shapiro's tube law (16.188): valid for α ≤ 1 (collapse branch, P ≤ 0).

    For α > 1 (distension) Fig. 16.19 rises almost vertically and the fit fails: the shapiro form caps P at 1 and
    c² = nK_pα⁻ⁿ/ρ falls as the tube inflates (non-physical).  Scalars: RuntimeWarning, the extrapolated value is
    still returned.  Arrays: RuntimeWarning, NaN at the offending entries (the trap-T16 pattern).
    Returns (al as float array, mask of entries to blank — all False for a scalar).
    """
    al = _F(al)
    bad = al > ALPHA_MAX_COLLAPSE * (1.0 + _ALPHA_ROUNDOFF)
    if np.any(bad):
        if al.ndim == 0:
            warnings.warn(f"{name} = {float(al):.4g} > 1: Shapiro's tube law (16.188) is fitted to the collapse branch "
                          "(α ≤ 1, Fig. 16.19) and is extrapolated here", RuntimeWarning, stacklevel=stacklevel)
            bad = np.zeros_like(bad)
        else:
            warnings.warn(f"{name} > 1 at {int(np.count_nonzero(bad))} point(s): outside the fitted collapse branch of "
                          "(16.188) (α ≤ 1, Fig. 16.19); returning NaN there", RuntimeWarning, stacklevel=stacklevel)
    return al, bad


def _tube_law_raw(al, n: float, form: str):
    """Unguarded (16.188) for internal root finding (brackets may probe α > 1); P [–]."""
    if form == "shapiro":
        return 1.0 - al ** (-float(n))  # Eq. (16.188)
    if form == "pure":
        return -al ** (-float(n))
    raise ValueError('form must be "shapiro" or "pure"')


def collapsible_tube_law(area_ratio, *, n: float = 1.5, form: str = "shapiro"):
    """Shapiro's collapsible tube law P = (p − p_e)/K_p as a function of α = A/A₀:
    "shapiro" −P ≈ α⁻ⁿ − 1 (P = 1 − α⁻ⁿ), "pure" −P = α⁻ⁿ (the two fits of Fig. 16.19).

    Book: §16.3, Eqs. (16.186), (16.188) (n = 3/2); Fig. 16.19.  P < 0 means collapse (trap T9).
    **Valid for α ≤ 1 (collapse branch, Fig. 16.19).**  For α > 1 the measured curve rises almost vertically and the
    fit fails (the shapiro form caps P at 1): a scalar α > 1 gives a RuntimeWarning (value extrapolated), array
    entries α > 1 give NaN plus a RuntimeWarning.
    Parameters: area_ratio α [–] (> 0, ≤ 1);  n [–];  form.  Returns P [–].
    Assumptions: tube law depends on α only; longitudinal tension neglected.  Label: analytic.
    """
    al = _F(area_ratio)
    require_positive("area_ratio", al)
    al, bad = _collapse_guard(al)
    P = _tube_law_raw(al, n, form)  # Eq. (16.188)
    return _S(np.where(bad, np.nan, P))


def collapsible_tube_law_derivatives(area_ratio, n: float = 1.5, form: str = "shapiro"):
    """(P, dP/dα, d²P/dα²) of :func:`collapsible_tube_law` (both forms share the derivatives nα^{−n−1}, −n(n+1)α^{−n−2}).

    Book: §16.3, (16.188), used in (16.189), (16.206).  Valid for α ≤ 1 (same guard as :func:`collapsible_tube_law`:
    scalar α > 1 warns, array entries α > 1 are NaN in all three).  Returns a tuple of floats/arrays.  Label: analytic.
    """
    al = _F(area_ratio)
    require_positive("area_ratio", al)
    al, bad = _collapse_guard(al)
    n = float(n)
    P = _tube_law_raw(al, n, form)
    dP = n * al ** (-n - 1.0)
    d2P = -n * (n + 1.0) * al ** (-n - 2.0)
    return tuple(_S(np.where(bad, np.nan, v)) for v in (P, dP, d2P))


def alpha_from_tube_law(P, n: float = 1.5, form: str = "shapiro"):
    """Invert the tube law: α = (1 − P)^{−1/n} ("shapiro", P < 1) or (−P)^{−1/n} ("pure", P < 0).

    Book: §16.3, Eq. (16.188) inverted (ours).  Parameters: P [–];  n;  form.  Returns α [–].
    The fit is valid for α ≤ 1, i.e. P ≤ 0 (shapiro form): P > 0 (distension) warns for a scalar (value returned) and
    gives NaN plus a RuntimeWarning for array entries.
    Validation: V1 round trip with :func:`collapsible_tube_law`. Label: analytic."""
    P = _F(P)
    if form == "shapiro":
        if np.any(P >= 1.0):
            raise ValueError("the shapiro form has P < 1")
        al = (1.0 - P) ** (-1.0 / float(n))
    elif form == "pure":
        if np.any(P >= 0.0):
            raise ValueError("the pure form has P < 0")
        al = (-P) ** (-1.0 / float(n))
    else:
        raise ValueError('form must be "shapiro" or "pure"')
    al, bad = _collapse_guard(al, name="alpha_from_tube_law (P > 0: distended)")
    return _S(np.where(bad, np.nan, al))


def collapsible_wave_speed(area_ratio, Kp, rho, n: float = 1.5):
    """Long-wave speed in a collapsible tube, c² = (A/ρ)∂(p − p_e)/∂A = nK_pα⁻ⁿ/ρ.

    Book: §16.3, Eq. (16.189) (from (16.185) with (16.188); derivation a-d37).  Valid for α ≤ 1 (the range of
    (16.188); for α > 1 this c falls as the tube inflates, which is non-physical): scalar α > 1 warns, array entries
    α > 1 are NaN.
    Parameters: area_ratio α [–];  Kp [Pa];  rho [kg/m³];  n [–].  Returns c [m/s].
    Validation: V1 = :func:`bramwell_hill_speed` with a finite-difference dp/dA. Label: analytic.
    """
    al, bad = _collapse_guard(area_ratio)
    c = np.sqrt(float(n) * _F(Kp) * al ** (-float(n)) / _F(rho))  # Eq. (16.189)
    return _S(np.where(bad, np.nan, c))


def speed_index(u, c):
    """Shapiro's speed index S = u/c (the collapsible tube's Mach / Froude number).

    Book: §16.3, Eq. (16.195).  S < 1 subcritical, S = 1 choking (flow limitation), S > 1 supercritical.
    Parameters: u [m/s];  c [m/s].  Returns S [–].  Label: analytic.
    """
    return _S(_F(u) / _F(c))  # Eq. (16.195)


def shapiro_M(area_ratio, n: float = 1.5, *, tube_law: Callable | None = None, h: float = 1e-5):
    """Shapiro's coefficient M = 3 + α(∂²P/∂α²)/(∂P/∂α) (not a Mach number; trap T2).

    Book: §16.3, Eq. (16.206).  For the power law P = 1 − α⁻ⁿ: M = 3 − (n + 1) = 2 − n (= ½ for n = 3/2, trap T9).
    A general ``tube_law`` callable P(α) is differentiated by central differences (step h·α).
    Parameters: area_ratio α [–];  n [–];  tube_law callable (optional);  h relative step.  Returns M [–].
    Validation: V1 2 − n for the power law. Label: analytic.
    """
    al = _F(area_ratio)
    if tube_law is None:
        return _S(np.full_like(al, 2.0 - float(n)) + 0.0 * al)  # Eq. (16.206) for P = 1 − α⁻ⁿ
    d = h * al
    P1 = (_F(tube_law(al + d)) - _F(tube_law(al - d))) / (2.0 * d)
    P2 = (_F(tube_law(al + d)) - 2.0 * _F(tube_law(al)) + _F(tube_law(al - d))) / d ** 2
    return _S(3.0 + al * P2 / P1)  # Eq. (16.206)


def friction_R_per_length(area_ratio, A0, mu, *, f_L_prime: float = 8.0):
    """Laminar viscous resistance per unit length R (so that τ_ws/A = RQ) from (16.193):
    R = 4μf′_L/(α²A₀d₀²), d₀ = 2√(A₀/π); f′_L = 8 gives Poiseuille's 8μ/(πa₀⁴) at α = 1.

    Book: §16.3, Eq. (16.193) (τ_ws dx/A = (μu/d₀)(1/α)(4f′_L dx/d₀)) with u = Q/(αA₀).  DEVIATION: the book defines
    f′_L = (A/A_e)f_L with A_e undefined (slip S21, unresolved without Shapiro 1977a); we take f′_L constant (keyword).
    Parameters: area_ratio α [–];  A0 [m²];  mu [Pa s];  f_L_prime [–].  Returns R [Pa s/m⁴].
    Validation: V1 Poiseuille at α = 1. Label: analytic.
    """
    d0 = 2.0 * np.sqrt(_F(A0) / np.pi)
    return _S(4.0 * _F(mu) * f_L_prime / (_F(area_ratio) ** 2 * _F(A0) * d0 ** 2))


def wall_friction_term(u, area_ratio, d0, rho, mu, *, regime: str, f: float | None = None,
                       d_e: float | None = None, f_L_prime: float = 8.0):
    """Wall-friction pressure gradient τ_ws/A [Pa/m]: turbulent ½ρu²·4f_T/d_e (16.192) or laminar
    (μu/d₀)(1/α)(4f′_L/d₀) (16.193).

    Book: §16.3, Eqs. (16.192)–(16.193).  Parameters: u [m/s];  area_ratio [–];  d0 [m];  rho;  mu;  regime;
    f : Fanning f_T (turbulent, required);  d_e : hydraulic diameter 4A/s [m] (turbulent, default d₀√α);  f_L_prime.
    Returns [Pa/m].  Label: analytic.
    """
    u, al = _F(u), _F(area_ratio)
    if regime == "turbulent":
        if f is None:
            raise ValueError("turbulent friction needs the Fanning coefficient f")
        de = float(d0) * np.sqrt(al) if d_e is None else float(d_e)
        return _S(0.5 * float(rho) * u ** 2 * 4.0 * float(f) / de)  # Eq. (16.192)
    if regime == "laminar":
        return _S(float(mu) * u / float(d0) / al * 4.0 * f_L_prime / float(d0))  # Eq. (16.193)
    raise ValueError('regime must be "laminar" or "turbulent"')


@dataclass
class CollapsibleTube:
    """Parameters of a steady collapsible-tube flow for Shapiro's equations (ours; every field a constant or a callable
    of x [m]).

    A0 : neutral area [m²];  Kp : stiffness parameter [Pa];  pe : external pressure [Pa];  dzdx : slope of the
    elevation [–] (−1 for a vertical tube draining downward, the giraffe jugular);  rho [kg/m³];  mu [Pa s];  n : tube
    law exponent;  form : "shapiro" | "pure";  friction : include RQ;  f_L_prime : see :func:`friction_R_per_length`;
    g [m/s²];  h_fd : step for derivatives of callables [m];  extra : other keys (e.g. the length L).
    Book: §16.3, the independent quantities of Shapiro's theory listed after (16.195) (dA₀, dp_e, g dz, f dx, dK_p).
    Label: V7 (structure).
    """

    A0: float | Callable = 1.1310e-4
    Kp: float | Callable = 200.0
    pe: float | Callable = 0.0
    dzdx: float | Callable = 0.0
    rho: float = 1000.0
    mu: float = 1.0e-3
    n: float = 1.5
    form: str = "shapiro"
    friction: bool = True
    f_L_prime: float = 8.0
    g: float = G0
    h_fd: float = 1e-6
    extra: dict = field(default_factory=dict)

    def val(self, name, x):
        f = getattr(self, name)
        return float(f(x)) if callable(f) else float(f)

    def der(self, name, x):
        f = getattr(self, name)
        if not callable(f):
            return 0.0
        h = self.h_fd
        return (float(f(x + h)) - float(f(x - h))) / (2.0 * h)

    @property
    def L(self):
        return self.extra.get("L")


_TUBE_KEYS = ("A0", "Kp", "pe", "dzdx", "rho", "mu", "n", "form", "friction", "f_L_prime", "g", "h_fd")


def _as_tube(prm) -> CollapsibleTube:
    """Accept the design's ``prm`` dict (keys L, A0, Kp, n, rho, mu; optional pe, dzdx, form, friction, f_L_prime, g —
    constants or callables of x) or a :class:`CollapsibleTube`."""
    if isinstance(prm, CollapsibleTube):
        return prm
    if not isinstance(prm, dict):
        raise TypeError("prm must be a dict or a CollapsibleTube")
    kw = {k: prm[k] for k in _TUBE_KEYS if k in prm}
    extra = {k: v for k, v in prm.items() if k not in _TUBE_KEYS}
    return CollapsibleTube(**kw, extra=extra)


def collapsible_terms(area_ratio, Q, prm, *, x: float = 0.0, dA0dx: float | None = None, dpedx: float | None = None,
                      dzdx: float | None = None, dKpdx: float | None = None) -> dict:
    """The drivers of Shapiro's area equation (16.204) at one station, each already divided as in the equation
    (E9 term bars; ours):  (1 − S²)(1/α)dα/dx = S²(1/A₀)dA₀/dx − (1/ρc²)[dp_e/dx + ρg dz/dx + RQ + P dK_p/dx].

    Book: §16.3, Eq. (16.204) (derivation D30: substitute (16.200) into (16.194), collect dα).  The driver sum is
    formed **before** dividing by 1 − S² (ch15 lesson, trap T10).  Friction: laminar R = 8πμ/(αA₀)² per unit length
    (:func:`friction_R_per_length` with f′_L = 8), so that τ_ws dx/A = RQ dx.
    Valid for α ≤ 1 (the collapse branch of (16.188), Fig. 16.19): α > 1 gives a RuntimeWarning (values extrapolated).
    Parameters
    ----------
    area_ratio α [–];  Q [m³/s];  prm : dict (L, A0, Kp, n, rho, mu; optional pe, dzdx, form, friction, f_L_prime,
    g — constants or callables of x) or :class:`CollapsibleTube`;  x : station [m] (where callables are evaluated);
    dA0dx [m], dpedx [Pa/m], dzdx [–], dKpdx [Pa/m] : local gradients (override those implied by ``prm``).
    Returns
    -------
    dict: S2, c [m/s], u [m/s], P, area, external, gravity, friction, stiffness (the five terms of the right side
    [1/m]), driver (their sum), one_minus_S2, dalpha_dx [1/m] (NaN at S = 1).
    Label: analytic.
    """
    prm = _as_tube(prm)
    x = float(x)
    al = float(area_ratio)
    if al > ALPHA_MAX_COLLAPSE * (1.0 + _ALPHA_ROUNDOFF):  # (16.188) range guard (scalar fast path)
        _collapse_guard(al, stacklevel=4)
    A0, Kp = prm.val("A0", x), prm.val("Kp", x)
    n = float(prm.n)
    an = al ** (-n)
    P = 1.0 - an if prm.form == "shapiro" else -an  # (16.188), scalar fast path
    dP = n * an / al
    rc2 = al * Kp * dP  # ρc² = αK_p ∂P/∂α (16.198)
    c = (rc2 / prm.rho) ** 0.5
    u = Q / (al * A0)
    S2 = u * u / (c * c)
    # friction_R_per_length inlined: R = 4μf′_L/(α²A₀d₀²), d₀² = 4A₀/π
    R = 4.0 * prm.mu * prm.f_L_prime / (al * al * A0 * (4.0 * A0 / np.pi)) if prm.friction else 0.0
    gA0 = prm.der("A0", x) if dA0dx is None else float(dA0dx)
    gpe = prm.der("pe", x) if dpedx is None else float(dpedx)
    gz = prm.val("dzdx", x) if dzdx is None else float(dzdx)
    gKp = prm.der("Kp", x) if dKpdx is None else float(dKpdx)
    t_area = S2 / A0 * gA0
    t_ext = -gpe / rc2
    t_grav = -prm.rho * prm.g * gz / rc2
    t_fric = -R * Q / rc2
    t_stiff = -P * gKp / rc2
    drv = t_area + t_ext + t_grav + t_fric + t_stiff  # Eq. (16.204), right side
    om = 1.0 - S2
    return dict(S2=S2, c=c, u=u, P=P, area=t_area, external=t_ext, gravity=t_grav, friction=t_fric,
                stiffness=t_stiff, driver=drv, one_minus_S2=om, dalpha_dx=al * drv / om if om != 0.0 else np.nan)


def collapsible_rhs(x, y, prm, Q: float):
    """dα/dx of Shapiro's equation (16.204) for P = P(α): α·driver/(1 − S²).  y = [α].  Book: §16.3, (16.204).

    Parameters: x [m];  y = [α];  prm (dict or :class:`CollapsibleTube`);  Q [m³/s].  Returns [dα/dx] [1/m].
    Label: analytic.
    """
    return [collapsible_terms(y[0], Q, prm, x=x)["dalpha_dx"]]


def collapsible_speed_index_rhs(x, area_ratio, Q, prm) -> float:
    """(1/S²)dS²/dx from Shapiro's second equation (16.205) (for the consistency check against S² computed from α;
    ours):  (1 − S²)(1/S²)dS²/dx = (1/A₀)(dA₀/dx)[−2 + (2 − M)S²] + (M/ρc²)[dp_e/dx + ρg dz/dx + RQ]
    + (1/ρc²)(dK_p/dx)[MP − (1 − S²)α dP/dα].

    Book: §16.3, Eqs. (16.205)–(16.206).  Returns (1/S²)dS²/dx [1/m] (NaN at S = 1).  Label: analytic.
    """
    prm = _as_tube(prm)
    tr = collapsible_terms(area_ratio, Q, prm, x=x)
    al = float(area_ratio)
    A0, Kp = prm.val("A0", x), prm.val("Kp", x)
    P, dP, d2P = collapsible_tube_law_derivatives(al, prm.n, prm.form)
    M = 3.0 + al * float(d2P) / float(dP)  # Eq. (16.206)
    rc2 = al * Kp * float(dP)
    S2 = tr["S2"]
    R = float(friction_R_per_length(al, A0, prm.mu, f_L_prime=prm.f_L_prime)) if prm.friction else 0.0
    bracket = prm.der("pe", x) + prm.rho * prm.g * prm.val("dzdx", x) + R * Q
    rhs = (prm.der("A0", x) / A0 * (-2.0 + (2.0 - M) * S2) + M / rc2 * bracket
           + prm.der("Kp", x) / rc2 * (M * float(P) - (1.0 - S2) * al * float(dP)))  # Eq. (16.205)
    return rhs / (1.0 - S2) if S2 != 1.0 else np.nan


def collapsible_march(Q, alpha_in, x_span, prm, *, rtol: float = 1e-9, atol: float = 1e-12,
                      max_step: float | None = None, n_out: int = 400, S_margin: float = 1e-9) -> dict:
    """March α(x) along a collapsible tube with Shapiro's equation (16.204) (solve_ivp RK45) and stop at the critical
    point: a terminal event where S² crosses 1 ∓ ``S_margin`` (choking / flow limitation).

    Book: §16.3, Eqs. (16.203)–(16.206); Fig. 16.20 (smooth transitions need a vanishing right side at S = 1).
    Parameters
    ----------
    Q : flow [m³/s];  alpha_in : α at x_span[0] [–];  x_span : (x0, x1) [m];  prm : dict or :class:`CollapsibleTube`;
    rtol, atol, max_step : solver controls;  n_out : output points;  S_margin : event offset of S² from 1.
    Returns
    -------
    dict: x, alpha, S, p (= p_e + K_pP(α)), u, c, choked (bool), x_choke (or None).
    Assumptions: steady, 1-D, incompressible, P = P(α) (16.203).
    Validation: V1 frictionless case equals :func:`collapsible_frictionless_alpha` (Bernoulli); V1 S² from α obeys
    (16.205); signs of dα/dx for S < 1 and S > 1 (friction) (independent-route cross-checks; no refinement order
    asserted). Label: analytic.
    Regime: the tube law (16.188) is valid for α ≤ 1 (collapse branch, Fig. 16.19); states α > 1 warn.
    """
    prm = _as_tube(prm)
    x0, x1 = (float(v) for v in x_span)
    S_in2 = collapsible_terms(alpha_in, Q, prm, x=x0)["S2"]
    target = 1.0 - S_margin if S_in2 < 1.0 else 1.0 + S_margin

    def ev(x, y):
        return collapsible_terms(y[0], Q, prm, x=x)["S2"] - target

    ev.terminal = True
    max_step = (x1 - x0) / 50.0 if max_step is None else max_step
    sol = solve_ivp(lambda x, y: collapsible_rhs(x, y, prm, Q), (x0, x1), [float(alpha_in)], events=ev,
                    dense_output=True, rtol=rtol, atol=atol, max_step=max_step)
    x_end = float(sol.t[-1])
    choked = bool(sol.t_events[0].size > 0)
    xs = np.linspace(x0, x_end, int(n_out))
    al = sol.sol(xs)[0]
    S = np.empty_like(xs)
    p = np.empty_like(xs)
    u = np.empty_like(xs)
    c = np.empty_like(xs)
    for i, (xx, aa) in enumerate(zip(xs, al)):
        tr = collapsible_terms(aa, Q, prm, x=xx)
        S[i], u[i], c[i] = np.sqrt(tr["S2"]), tr["u"], tr["c"]
        p[i] = prm.val("pe", xx) + prm.val("Kp", xx) * tr["P"]  # Eq. (16.203)
    return dict(x=xs, alpha=al, S=S, p=p, u=u, c=c, choked=choked, x_choke=x_end if choked else None)


def choking_flux(area_ratio, A0, Kp, rho, n: float = 1.5):
    """Flow at which the speed index reaches 1 at a section of area ratio α: Q* = αA₀c(α) (S = 1).

    Book: §16.3, (16.189), (16.195) (flow limitation at S = 1, the Laval-throat analogy).
    Parameters: area_ratio [–];  A0 [m²];  Kp [Pa];  rho [kg/m³];  n.  Returns Q* [m³/s].  Label: analytic.
    """
    al = _F(area_ratio)
    return _S(al * _F(A0) * _F(collapsible_wave_speed(al, Kp, rho, n)))


def collapsible_frictionless_alpha(Q, H, A0, Kp, rho, n: float = 1.5, *, form: str = "shapiro",
                                   branch: str = "sub"):
    """Area ratio of steady frictionless flow from Bernoulli: K_pP(α) + ½ρQ²/(α²A₀²) = H (uniform p_e, K_p, level;
    ours — the closed-form reference for the throat transition of Fig. 16.20(a)).

    Book: §16.3, the frictionless limit of (16.204) (dα/dx from A₀(x) only).  The two branches meet at S = 1 (where
    d/dα of the left side vanishes); "sub" is α above that point (S < 1), "super" below (S > 1).
    Parameters: Q [m³/s];  H : K_pP + ½ρu² [Pa] (a constant);  A0 [m²] (local);  Kp [Pa];  rho;  n;  form;  branch.
    Returns α [–] (NaN if H is below the minimum: no steady solution at this A₀).
    Regime: the tube law (16.188) is valid for α ≤ 1 (collapse branch, Fig. 16.19).  The sub branch stays at α ≤ 1
    everywhere iff H ≤ 0 (for n = 3/2 and S = 1 at a throat of area ratio α_t: iff α_t ≤ 0.25^{2/3} = 0.397); a root
    α > 1 is returned with a RuntimeWarning (extrapolated fit).  The root search itself probes α > 1 silently.
    Label: analytic.
    """
    A0, Kp, rho, Q, H = (float(v) for v in (A0, Kp, rho, Q, H))
    n = float(n)
    if not n < 2.0:
        raise ValueError("collapsible_frictionless_alpha needs n < 2 (S² ∝ α^(n−2) must fall with α)")
    # Critical section S = 1: S² = ρQ²/(α³A₀²K_p dP/dα) with dP/dα = nα^{−n−1} for both forms of (16.188), so
    # S² = ρQ²α^{n−2}/(A₀²K_p n) and the root is closed-form — no bracket that wide upstream sections (α* ≪ 1e-6,
    # α* = α_t(A₀,t/A₀)⁴ at n = 3/2) could fall outside of.
    a_star = float((rho * Q ** 2 / (A0 ** 2 * Kp * n)) ** (1.0 / (2.0 - n)))  # Eq. (16.189) with S = 1

    def f(al):
        return Kp * _tube_law_raw(al, n, form) + 0.5 * rho * Q ** 2 / (al * A0) ** 2 - H

    def done(al):
        al = float(al)
        if al > ALPHA_MAX_COLLAPSE * (1.0 + _ALPHA_ROUNDOFF):
            _collapse_guard(al, name="collapsible_frictionless_alpha root", stacklevel=4)
        return al

    f_star = f(a_star)
    if f_star > 1e-12 * max(abs(H), Kp):
        return float("nan")
    if f_star >= 0.0:  # at the critical section itself (round-off): both branches meet at α*
        return done(a_star)
    if branch == "sub":
        hi = a_star * 2.0
        while f(hi) < 0:
            hi *= 2.0
            if hi > 1e8:
                return float("nan")
        return done(brentq(f, a_star, hi, xtol=1e-15))
    lo = a_star / 2.0
    while f(lo) < 0:  # the kinetic term ∝ α⁻² outgrows K_pα⁻ⁿ (n < 2) as α → 0, so this terminates
        lo /= 2.0
        if lo < 1e-300:
            return float("nan")
    return done(brentq(f, lo, a_star, xtol=1e-15 * a_star))


def _is_uniform(prm: CollapsibleTube) -> bool:
    return not any(callable(getattr(prm, k)) for k in ("A0", "Kp", "pe", "dzdx"))


def _alpha_star(Q, prm: CollapsibleTube) -> float:
    """Area ratio at which S = 1 for flux Q in a uniform tube (S² = ρQ²/(α³A₀²K_pP′(α)))."""
    A0, Kp, n = prm.val("A0", 0.0), prm.val("Kp", 0.0), float(prm.n)
    # S² ∝ α^{n−2} (P′ = nα^{−n−1}): closed form for the power law
    return float((prm.rho * Q ** 2 / (A0 ** 2 * Kp * n)) ** (1.0 / (2.0 - n)))


def _length_to(alpha_end, alpha_in, Q, prm: CollapsibleTube) -> float:
    """Distance from α_in to α_end in a uniform tube: x = ∫ (1 − S²)/(α·driver) dα (ours; regular at S = 1)."""
    def dx_dlna(s):  # in s = ln α the integrand (1 − S²)/driver is smooth over many decades of α
        tr = collapsible_terms(np.exp(s), Q, prm)
        return tr["one_minus_S2"] / tr["driver"]

    val, _ = quad(dx_dlna, np.log(alpha_in), np.log(alpha_end), epsabs=1e-13, epsrel=1e-11, limit=200)  # x [m]
    return float(val)


def _outlet_state(Q, p1, pe, prm: CollapsibleTube, L):
    """Outlet state of a uniform tube with constant p_e (α_in from p₁); None if it chokes (S = 1) before L.

    Uniform tube: every term of (16.204) depends on α only, so x(α) = ∫(1 − S²)/(α·driver) dα is a regular quadrature
    (the integrand vanishes at S = 1) — robust where an x-march meets the square-root singularity at the choke.
    Non-uniform tubes fall back to :func:`collapsible_march`.
    """
    al_in = float(alpha_from_tube_law((p1 - pe) / prm.val("Kp", 0.0), prm.n, prm.form))
    if collapsible_terms(al_in, Q, prm)["S2"] >= 1.0:
        return None
    if not _is_uniform(prm):
        m = collapsible_march(Q, al_in, (0.0, L), prm, n_out=3)
        return None if m["choked"] else m
    a_st = _alpha_star(Q, prm)
    if collapsible_terms(al_in, Q, prm)["driver"] >= 0.0:  # area does not fall: no choking downstream
        L_max = np.inf
    else:
        L_max = _length_to(a_st, al_in, Q, prm)
    if L_max < L:
        return None
    a_out = brentq(lambda a: _length_to(a, al_in, Q, prm) - L, a_st, al_in, xtol=1e-15, rtol=1e-13)
    tr = collapsible_terms(a_out, Q, prm)
    p_out = prm.val("pe", L) + prm.val("Kp", L) * tr["P"]
    return dict(x=np.array([0.0, L]), alpha=np.array([al_in, a_out]), S=np.array([0.0, np.sqrt(tr["S2"])]),
                p=np.array([p1, p_out]), L_max=L_max)


def _tube_length(prm: CollapsibleTube, L):
    L = prm.L if L is None else L
    if L is None:
        raise ValueError("the tube length L is required (prm['L'] or L=)")
    return float(L)


def starling_state(p1, p2, pe, prm, *, L: float | None = None, _lim: dict | None = None) -> dict:
    """Steady flow through a collapsible tube in a chamber (Starling resistor), uniform tube, laminar friction: find Q
    with outlet pressure p₂; if p₂ is below the value at which S reaches 1 at the outlet, Q stays at that limiting
    value (flow limitation, Case 1; ours).

    Book: §16.3 "Starling resistor" Case 1 (flow limitation) with Shapiro's (16.204); the computation is ours: the
    inlet area follows from p₁ − p_e = K_pP(α_in) (16.203).  For a uniform tube every term of (16.204) depends on α
    only, so the code integrates x(α) = ∫(1 − S²)/(α·driver) dα (in ln α; regular at S = 1, where the integrand
    vanishes): p₂(Q) from the α reached at x = L, and Q_lim from L_needed(Q) = L with the outlet exactly at S = 1.
    DEVIATION (ours, numerical): an x-march with a terminal event creeps along the square-root singularity at S = 1
    and overestimates Q_lim (by 2e-3 for the design's vein); it is kept only for non-uniform tubes (bisection).
    Parameters: p1 inlet pressure [Pa];  p2 outlet pressure [Pa];  pe chamber pressure [Pa];  prm dict (L, A0, Kp, n,
    rho, mu) or :class:`CollapsibleTube`;  L : length [m] (overrides prm["L"]).
    Returns dict: Q [m³/s], alpha_out [–], S_out [–], limited (bool), Q_lim [m³/s], p2_lim [Pa].
    Assumptions: as :func:`collapsible_march`; downstream of a choked outlet the flow adjusts outside the model.
    Cost: ≈ 0.3 s (one bisection for the limit, one root find).
    Regime: the tube law (16.188) is valid for α ≤ 1, so p₁ ≤ p_e (inlet not distended); p₁ > p_e warns.
    Validation: V1 the plateau Q equals the flux at which S_out → 1 (closed-form cross-check; no refinement order
    asserted). Label: analytic.
    """
    prm = _as_tube(prm)
    L = _tube_length(prm, L)
    p1, p2, pe = float(p1), float(p2), float(pe)
    lim = _flow_limit(p1, pe, prm, L) if _lim is None else _lim
    Q_lim, p2_lim = lim["Q_lim"], lim["p2_lim"]
    base = dict(Q_lim=Q_lim, p2_lim=p2_lim)
    if p2 >= p1:
        al_in = float(alpha_from_tube_law((p1 - pe) / prm.val("Kp", 0.0), prm.n, prm.form))
        return dict(Q=0.0, alpha_out=al_in, S_out=0.0, limited=False, **base)
    if p2 <= p2_lim:
        return dict(Q=Q_lim, alpha_out=lim["alpha_out"], S_out=lim["S_out"], limited=True, **base)

    def g(Q):
        m = _outlet_state(Q, p1, pe, prm, L)
        return (m["p"][-1] if m is not None else -np.inf) - p2

    Q = brentq(g, Q_lim * 1e-9, Q_lim * (1.0 - 1e-7), xtol=Q_lim * 1e-11, rtol=1e-12)
    m = _outlet_state(Q, p1, pe, prm, L)
    return dict(Q=Q, alpha_out=float(m["alpha"][-1]), S_out=float(m["S"][-1]), limited=False, **base)


def _flow_limit(p1, pe, prm, L, n_bisect: int = 50) -> dict:
    al_in = float(alpha_from_tube_law((p1 - pe) / prm.val("Kp", 0.0), prm.n, prm.form))
    Q_hi = float(choking_flux(al_in, prm.val("A0", 0.0), prm.val("Kp", 0.0), prm.rho, prm.n))
    if _is_uniform(prm):  # root of L_needed(Q) = L, the length over which friction brings S to 1 (α-integration)
        def g(Q):
            return _length_to(_alpha_star(Q, prm), al_in, Q, prm) - L

        Q_lim = brentq(g, Q_hi * 1e-6, Q_hi * (1.0 - 1e-12), xtol=Q_hi * 1e-14, rtol=1e-13)
        a_st = _alpha_star(Q_lim, prm)
        P_st = float(collapsible_tube_law(a_st, n=prm.n, form=prm.form))
        return dict(Q_lim=Q_lim, p2_lim=prm.val("pe", L) + prm.val("Kp", L) * P_st, S_out=1.0, alpha_out=a_st)
    lo, hi = 0.0, Q_hi
    for _ in range(int(n_bisect)):
        mid = 0.5 * (lo + hi)
        if _outlet_state(mid, p1, pe, prm, L) is None:
            hi = mid
        else:
            lo = mid
    m = _outlet_state(lo, p1, pe, prm, L)
    return dict(Q_lim=lo, p2_lim=float(m["p"][-1]), S_out=float(m["S"][-1]), alpha_out=float(m["alpha"][-1]))


def flow_limitation_curve(p2_list, *, p1, pe, prm, L: float | None = None) -> dict:
    """Flow against outlet pressure for a Starling resistor (Case 1): Q rises as p₂ falls, then plateaus at Q_choke
    when S = 1 at the outlet (ours).

    Book: §16.3 "Starling resistor", Case 1, with (16.204).
    Parameters: p2_list [Pa];  p1, pe [Pa];  prm dict (L, A0, Kp, n, rho, mu) or :class:`CollapsibleTube`;  L [m]
    (overrides prm["L"]).
    Returns dict: p2, Q, S_out, limited (arrays), Q_choke [m³/s], p2_choke [Pa] (Q_lim, p2_lim aliases).
    Regime: as :func:`starling_state` (p₁ ≤ p_e, α ≤ 1).
    Cost ≈ 0.1 s per point.  Validation: V1 Q_choke and p₂ at the choke equal the closed-form α-integral (no
    refinement order asserted).  Label: analytic.
    """
    prm = _as_tube(prm)
    L = _tube_length(prm, L)
    p2v = _F(p2_list)
    lim = _flow_limit(float(p1), float(pe), prm, L)
    out = [starling_state(p1, v, pe, prm, L=L, _lim=lim) for v in np.atleast_1d(p2v)]
    return dict(p2=p2v, Q=np.array([o["Q"] for o in out]), S_out=np.array([o["S_out"] for o in out]),
                limited=np.array([o["limited"] for o in out]), Q_choke=lim["Q_lim"], p2_choke=lim["p2_lim"],
                Q_lim=lim["Q_lim"], p2_lim=lim["p2_lim"])


def collapsible_influence_coefficients(S2, area_ratio, n: float = 1.5, form: str = "shapiro") -> dict:
    """Influence coefficients of Shapiro's equations (ours): dα/α and dS²/S² per unit of each independent change,
    the collapsible-tube twin of ch15's ``influence_coefficients`` (every entry has the factor 1/(1 − S²)).

    Book: §16.3, (16.204)–(16.206) rearranged (driver D = (dp_e + ρg dz + RQ dx)/(ρc²)).
    Parameters: S2 [–];  area_ratio α [–];  n;  form.
    Returns dict {driver: (coef in dα/α, coef in dS²/S²)} for drivers "dA0/A0", "D", "dKp/Kp".
    Label: analytic.
    """
    S2 = float(S2)
    al = float(area_ratio)
    P, dP, d2P = collapsible_tube_law_derivatives(al, n, form)
    M = 3.0 + al * float(d2P) / float(dP)
    k = 1.0 / (1.0 - S2)
    q = float(P) / (al * float(dP))  # P K_p/(ρc²)
    return {"dA0/A0": (S2 * k, (-2.0 + (2.0 - M) * S2) * k), "D": (-k, M * k),
            "dKp/Kp": (-q * k, (M * q - (1.0 - S2)) * k)}


collapsible_influence_matrix = collapsible_influence_coefficients  # analysis-plan name (I17)


def collapsible_differentials_sympy(printed_16_201: bool = False) -> dict:
    """Symbolic check of Shapiro's differential relations (ours):
    (1) the logarithmic derivative of c² = (α/ρ)K_p ∂P/∂α (16.198) against (16.201) — corrected last term
        (αK_p/ρc²)∂/∂x(∂P/∂α)dx, or the printed ∂/∂x(∂P/∂x) with ``printed_16_201=True`` (slip S23);
    (2) (16.204) and (16.205) re-derived from (16.194), (16.196), (16.200) (P = P(α)) and (16.202): residuals.

    Book: §16.3, Eqs. (16.194)–(16.206).  DEVIATION: slip S23 — the printed inner derivative ∂P/∂x should be ∂P/∂α.
    Returns dict: dc_difference (sympy expression; 0 for the corrected form, non-zero when P depends on x for the
    printed one), eq204_residual, eq205_residual (both 0), M_power_law (= 2 − n).
    Label: symbolic.
    """
    import sympy as sp

    al, x, rho, n = sp.symbols("alpha x rho n", positive=True)
    dal, dx = sp.symbols("dalpha dx")
    Pf = sp.Function("P")(al, x)
    Kp = sp.Function("K_p")(x)
    c2 = al / rho * Kp * sp.diff(Pf, al)  # (16.198)
    two_dc_c = sp.diff(c2, al) / c2 * dal + sp.diff(c2, x) / c2 * dx  # d(c²)/c² = 2dc/c
    inner = sp.diff(Pf, x) if printed_16_201 else sp.diff(Pf, al)
    book = ((1 + al * sp.diff(Pf, al, 2) / sp.diff(Pf, al)) * dal / al + sp.diff(Kp, x) / Kp * dx
            + al * Kp / (rho * c2) * sp.diff(inner, x) * dx)  # (16.201)
    dc_diff = sp.simplify(two_dc_c - book)
    # (2) P = P(α) only
    P1 = sp.Function("P")(al)
    A0, dA0, dpe, g, dz, RQdx, u, Kp0, dKp = sp.symbols("A0 dA0 dpe g dz RQdx u K_p dKp", real=True)
    c2a = al / rho * Kp0 * sp.diff(P1, al)
    S2 = u ** 2 / c2a
    dA_over_A = dal / al + dA0 / A0  # (16.196)
    dp = dpe + rho * c2a * dal / al + P1 * dKp  # (16.200) with ∂P/∂x = 0
    mom = dp + rho * g * dz + RQdx - rho * u ** 2 * dA_over_A  # (16.194)
    dal_sol = sp.solve(sp.Eq(mom, 0), dal)[0]
    eq204 = sp.simplify((1 - S2) * dal_sol / al - (S2 * dA0 / A0 - (dpe + rho * g * dz + RQdx + P1 * dKp) / (rho * c2a)))
    M = 3 + al * sp.diff(P1, al, 2) / sp.diff(P1, al)
    two_dc = (1 + al * sp.diff(P1, al, 2) / sp.diff(P1, al)) * dal_sol / al + dKp / Kp0
    dS2_S2 = -2 * dal_sol / al - 2 * dA0 / A0 - two_dc  # (16.202)
    rhs205 = (dA0 / A0 * (-2 + (2 - M) * S2) + M / (rho * c2a) * (dpe + rho * g * dz + RQdx)
              + 1 / (rho * c2a) * dKp * (M * P1 - (1 - S2) * al * sp.diff(P1, al)))
    eq205 = sp.simplify((1 - S2) * dS2_S2 - rhs205)
    Ppow = 1 - al ** (-n)
    Mpow = sp.simplify(3 + al * sp.diff(Ppow, al, 2) / sp.diff(Ppow, al))
    return dict(dc_difference=dc_diff, eq204_residual=eq204, eq205_residual=eq205, M_power_law=Mpow)


# ======================================================================================================================
# §16.4 plants (I19)
# ======================================================================================================================
def osmotic_pressure(c, T):
    """van 't Hoff osmotic pressure Π = RTc of a dilute solution (ours; the −RTc term of (16.221)).

    Book: §16.4, the term −RTc of Eq. (16.221).  Parameters: c solute concentration [mol/m³];  T [K].  Returns Π [Pa].
    Assumptions: dilute ideal solution.  Validation: V1 R from CODATA (``core.thermo``). Label: analytic.
    """
    return _S(R_MOLAR * _F(T) * _F(c))


def water_potential(p, c, T, z, *, rho: float = 1000.0, g: float = G0):
    """Water potential ψ = p − RTc + ρgz in SI (Pa): pressure, osmotic and gravitational terms.

    Book: §16.4, Eq. (16.221) (stated in bar with R = 83.141 cm³ bar/(mol K), c in mol/cm³, ρgz in CGS — the ρgz term
    must be divided by 10⁶ to be in bar, slip S30, loose).  SI here: no conversion needed.  ψ = 0 for pure water at
    atmospheric (gauge 0) pressure at z = 0; ψ < 0 under tension (trap T22).
    Parameters
    ----------
    p : hydrostatic (gauge) pressure [Pa];  c : solute concentration [mol/m³];  T [K];  z : height [m];
    rho [kg/m³];  g [m/s²].
    Returns
    -------
    psi [Pa].
    Assumptions: dilute solution.  Validation: V1 term by term; V1 book-unit version after conversion. Label: analytic.
    """
    return _S(_F(p) - R_MOLAR * _F(T) * _F(c) + _F(rho) * _F(g) * _F(z))  # Eq. (16.221)


def water_potential_terms(p, c, T, z, *, rho: float = 1000.0, g: float = G0) -> dict:
    """The three terms of the water potential (16.221) separately: pressure p, osmotic −RTc, gravity ρgz, and total.

    Book: §16.4, Eq. (16.221).  Parameters as :func:`water_potential`.  Returns dict pressure, osmotic, gravity,
    total [Pa].  Validation: V1 total = :func:`water_potential`. Label: analytic.
    """
    pr = _F(p) + 0.0 * _F(z)
    osm = -R_MOLAR * _F(T) * _F(c) + 0.0 * _F(z)
    gr = _F(rho) * _F(g) * _F(z)
    return dict(pressure=_S(pr), osmotic=_S(osm), gravity=_S(gr), total=_S(pr + osm + gr))


def water_potential_book_units(p_bar, c_mol_cm3, T, z_cm, rho_g_cm3: float = 1.0, g_cm_s2: float = 980.0):
    """(16.221) in the book's units, returning bar: ψ = p − RTc + ρgz/10⁶ with R = 83.141 cm³ bar/(mol K).

    Book: §16.4, Eq. (16.221) and the unit list after it.  DEVIATION: slip S30 — the printed formula adds ρgz in
    dyn/cm² to bar; we divide it by 10⁶ (1 bar = 10⁶ dyn/cm²).
    Parameters: p_bar [bar];  c_mol_cm3 [mol/cm³];  T [K];  z_cm [cm];  rho_g_cm3 [g/cm³];  g_cm_s2 [cm/s²].
    Returns psi [bar].  Validation: V1 = :func:`water_potential`/1e5 (with g = 9.80; the book's R = 83.141 is 0.004 %
    below CODATA's 83.1446 cm³ bar/(mol K)). Label: analytic.
    """
    return _S(_F(p_bar) - R_BOOK_CM3_BAR * _F(T) * _F(c_mol_cm3) + _F(rho_g_cm3) * g_cm_s2 * _F(z_cm) / 1.0e6)


def xylem_reynolds(u, a, nu):
    """Reynolds number of sap flow in a xylem conduit, Re = u·2a/ν (diameter).

    Book: §16.4 "Xylem" (Re = ud/ν; the book's quoted viscosity is 10× water's for its quoted Re — slip S28).
    Parameters: u [m/s];  a : conduit radius [m];  nu [m²/s].  Returns Re [–].  Label: analytic.
    """
    return _S(_F(u) * 2.0 * _F(a) / _F(nu))


def xylem_column(z, *, psi_soil, u, a, mu, rho: float = 1000.0, g: float = G0, c: float = 0.0,
                 T: float = 293.15) -> dict:
    """Pressure and water-potential terms up a xylem conduit with mean sap speed u (ours; E10, F8):
    p(z) = p_root − ρgz − (8μu/a²)z, with the root xylem at the soil potential, p_root = ψ_soil + RTc.

    Book: §16.4 (16.221) with Poiseuille's (16.17) in a rigid conduit (the xylem is treated as a rigid tube in quasi-
    steady Stokes flow there): Δp/L = 8μQ/(πa⁴) = 8μu/a².
    Parameters
    ----------
    z : heights [m];  psi_soil [Pa];  u : mean sap speed [m/s];  a : conduit radius [m];  mu [Pa s];  rho [kg/m³];
    g [m/s²];  c : sap solute concentration [mol/m³];  T [K].
    Returns
    -------
    dict: z, p, psi, gravity (ρgz), osmotic (−RTc), viscous (8μuz/a²) [Pa], grad_viscous (8μu/a²) [Pa/m].
    Validation: V1 ψ falls along the column by exactly the viscous drop. Label: analytic.
    """
    z = _F(z)
    osm = -float(osmotic_pressure(c, T))
    p_root = float(psi_soil) - osm
    grad = 8.0 * float(mu) * float(u) / float(a) ** 2
    visc = grad * z
    hyd = rho * g * z
    p = p_root - hyd - visc
    return dict(z=z, p=p, psi=_F(water_potential(p, c, T, z, rho=rho, g=g)), gravity=hyd, osmotic=osm + 0.0 * z,
                viscous=visc, grad_viscous=grad)


def cavitation_margin(p, *, p_cav):
    """How far the xylem pressure is above an assumed cavitation threshold, p − p_cav (> 0: the column holds; ours).

    Book: §16.4 (cavitation named).  ``p_cav`` is a required keyword (illustrative; not a book value).
    Parameters: p [Pa];  p_cav [Pa] (negative).  Returns margin [Pa].  Label: analytic.
    """
    return _S(_F(p) - float(p_cav))


def sieve_tube_pressure_drop(Q, mu, a, r_pore, L, ell, N) -> dict:
    """Pressure drop over one phloem sieve tube and one sieve plate (Rand et al. 1980):
    Δp = (8μQ/πa⁴)[L + (ℓ/N)(a/r)⁴] + 2Δp′, Δp′ = (8μQ/πa³)(a_e/r){0.57N[(a_e/r)³ − 1] − 1.5(1 − r/a_e)},
    a_e = a/√N.

    Book: §16.4 "Phloem", Eq. (16.222) (creeping conical flow through N pores; the book notes it is untested for N ≠ 1).
    Parameters
    ----------
    Q : flow [m³/s];  mu [Pa s];  a : sieve-tube radius [m];  r_pore : pore radius [m];  L : tube length [m];
    ell : plate thickness [m];  N : number of pores [–].
    Returns
    -------
    dict: total, tube (8μQL/πa⁴), plate_viscous (8μQℓ/(Nπr⁴)), plate_entry (2Δp′) [Pa], plate_share [–].
    Assumptions: creeping flow; the N pores fit in the plate, N r² ≤ a² (r ≤ a_e) — otherwise the formula's entry term
    turns negative, so the function raises (ours).
    Validation: V1 a_e = r ⇒ Δp′ = 0; ℓ = 0, N = 1, r = a ⇒ Poiseuille. Label: analytic.
    """
    Q, mu, a, r, L, ell, N = (float(v) for v in (Q, mu, a, r_pore, L, ell, N))
    ae = a / np.sqrt(N)
    if r > ae * (1.0 + 1e-12):
        raise ValueError(f"pore radius r = {r:g} m exceeds a/sqrt(N) = {ae:g} m: N pores of that size do not fit")
    tube = 8.0 * mu * Q / (np.pi * a ** 4) * L
    plate = 8.0 * mu * Q / (np.pi * a ** 4) * ell / N * (a / r) ** 4
    dpp = 8.0 * mu * Q / (np.pi * a ** 3) * (ae / r) * (0.57 * N * ((ae / r) ** 3 - 1.0) - 1.5 * (1.0 - r / ae))
    total = tube + plate + 2.0 * dpp  # Eq. (16.222)
    return dict(total=total, tube=tube, plate_viscous=plate, plate_entry=2.0 * dpp,
                plate_share=(plate + 2.0 * dpp) / total)


# ======================================================================================================================
# Tables (I20)
# ======================================================================================================================
def illustrative_inputs() -> dict:
    """Inputs used by the scripts, notebook and explainers — **ours**, chosen to be physiologically plausible and
    visibly not the book's numbers (rule 9).  Sources for the property defaults: :data:`BLOOD`.

    Book: §16.2–§16.4 (the quantities the chapter uses; values ours).  Returns dict of floats in SI.
    Label: V7 (structure).
    """
    return dict(
        rho=BLOOD["rho"], mu=BLOOD["mu"], mu_plasma=BLOOD["mu_p"], T_body=BLOOD["T"],
        HR=66.0, SV=78e-6, EDV=125e-6, ESV=47e-6, BSA=1.9, systole_share=0.35,
        R_periph=1.4e8, K_compliance=1.1e-8,
        p_ED=1.2e3, p_open=10.5e3, p_peak=15.5e3, p_ES=13.0e3, p_min=0.6e3,
        a_aorta=11e-3, h_aorta=1.6e-3, E_aorta=0.45e6, a_artery=3.5e-3, h_artery=0.6e-3, E_artery=0.8e6,
        L_vessel=0.1, a_arteriole=20e-6,
        tau_y=5.0e-3, Kc=0.055, H=0.45, C_F=0.30, delta_cfl=2.5e-6,
        G_womersley=600.0, Q_murray=5e-6, b_metabolic=2000.0,
        dean_delta=0.1, a_collapse=6e-3, Kp_collapse=200.0, L_collapse=0.3, rho_water=1000.0, mu_water=1.0e-3,
        tree_height=30.0, a_xylem=20e-6, u_xylem=1.0e-3, psi_soil=-0.3e6, T_plant=293.15, p_cav=-3.0e6,
        a_sieve=10e-6, r_pore=1.5e-6, L_sieve=300e-6, ell_plate=1.0e-6, N_pores=20.0, mu_phloem=1.5e-3,
        u_phloem=1.0e-4,
    )


def book_slips():
    """The 34 suspected slips of the chapter and what the code does (no book number is reproduced).

    Book: §16.2–§16.4 (analysis ch16 §9).  ``kind``: "slip" = an error on the page; "loose" = an inconsistency or
    imprecision.  ``changes_physics`` marks those coded corrected (with a ``printed_*`` option where an equation is
    involved).  Returns pandas.DataFrame: id, kind, where, printed, correct, code, changes_physics.
    Label: V7 (structure).
    """
    import pandas as pd

    rows = [
        ("S1", "loose", "(16.19) last member", "−(a/2)Δp/L", "+(a/2)Δp/L (τ > 0 by the book's own choice)", "poiseuille_wall_stress", False),
        ("S2", "loose", "(16.21), (16.24)", "−μ du/dr = −(r/2)Δp/L", "−μ du/dr = +(r/2)Δp/L; (16.27)–(16.28) are right", "two_layer_velocity", False),
        ("S3", "slip", "(16.30)", "πa²Δp/(4μ_pL)[a² − (1 − μ_p/2μ_c)(a − δ)⁴/a²]", "π(a − δ)²Δp/(4μ_pL)[a² − (1 − μ_p/2μ_c)(a − δ)²]", "two_layer_flow_rates(printed_16_30=True)", True),
        ("S4", "slip", "(16.44), (16.45)", "8μ/(πa₀)", "8μ/(πa₀⁴) (dimensions; (16.46) needs it)", "elastic_tube_flow_rate", True),
        ("S5", "slip", "(16.54)", "−∂((p − p_e)A)/∂x", "−A ∂p/∂x (same once linearised)", "tube1d_run(printed_16_54=True)", True),
        ("S6", "slip", "(16.76)", "right side +G/μ", "−G/μ (from (16.73)–(16.74); (16.77) consistent with it)", "womersley_ode_residual", False),
        ("S7", "loose", "(16.88), (16.89)", "√(ω/2ν); Σ to ∞", "√(nω/2ν) for harmonic n; Σ to N", "womersley_large_alpha(omega=nω)", False),
        ("S8", "slip", "(16.92)", "'shear rate' τ = ∂u/∂r with F(α_n)", "wall shear stress −μ∂u/∂r with F₂(α_n)", "womersley_wall_shear", True),
        ("S9", "loose", "(16.94)", "e^{inωt}", "e^{iωt} (single mode)", "womersley_flow_rate", False),
        ("S10", "slip", "(16.96)", "(2/(i^{1/2}α))(1 + 1/(2α))", "(2/(i^{1/2}α))(1 + i^{3/2}/(2α))", "womersley_F2_asymptote(printed_16_96=True)", True),
        ("S11", "loose", "worked example after (16.100)", "quoted X at the example α", "(16.91) gives X ≈ 3 % lower (Y within rounding); the amplitude ratio holds", "womersley_XY", False),
        ("S12", "loose", "(16.108)", "Eh/(1 − ν̂)", "Eh/(1 − ν̂²)", "—", False),
        ("S13", "loose", "(16.116)", "|k| = √(k₁ + k₂)", "√(k₁² + k₂²)", "elastic_tube_wave", False),
        ("S14", "loose", "text after (16.127)", "σ_rx ≈ μ ∂v/∂x", "σ_rx ≈ μ ∂u/∂r (what (16.128)–(16.129) use)", "elastic_tube_determinant_numeric", False),
        ("S15", "slip", "(16.140), (16.143)", "[(16μ/π²)K]^{1/6}", "[16μ/(π²K)]^{1/6}", "murray_radius(printed_16_140=True)", True),
        ("S16", "loose", "(16.149) with (16.158)", "Q_I − Q_R with Q_R = −Y₁P_R g", "the sign of Q_R counted twice; (16.152) is right", "junction_coefficients", False),
        ("S17", "loose", "(16.159), (16.160)", "reflected shape written f", "g (= f only because R is frequency-independent here)", "junction_waveforms", False),
        ("S18", "loose", "(16.171)", "κ = √δ 2Re", "consistent only with Re = aW̄₀/ν", "dean_number", False),
        ("S19", "loose", "(16.180)", "Ĝ = −8μW̄/a²", "+8μW̄/a² (so D > 0)", "dean_convert", False),
        ("S20", "loose", "after (16.182)", "'straight and curved, respectively'", "curved and straight", "curved_resistance_ratio", False),
        ("S21", "loose", "(16.193)", "f′_L = (A/A_e)f_L", "A_e undefined (unresolved)", "friction_R_per_length(f_L_prime=)", False),
        ("S22", "loose", "before (16.196)", "dA_e = (dA₀/dx)dx", "dA₀", "—", False),
        ("S23", "slip", "(16.201)", "(αK_p/ρc²)∂/∂x(∂P/∂x)dx", "(αK_p/ρc²)∂/∂x(∂P/∂α)dx", "collapsible_differentials_sympy(printed_16_201=True)", True),
        ("S24", "loose", "after (16.205)", "'the stiffness should decrease with x'", "with P < 0 the stiffness must increase (trap T9)", "collapsible_terms", False),
        ("S25", "loose", "(16.209) first member", "(r/2)dp/dx", "−(r/2)dp/dx", "casson_wall_stress", False),
        ("S26", "slip", "text after (16.209)", "the minimum pressure drop quoted in dyn/cm²", "1000× the printed figure; the printed mmHg figure is consistent", "casson_min_pressure_drop", True),
        ("S27", "slip", "units passage, §16.2", "1 cP = 0.01 P = 1 dyn s/cm² = 0.1 N s/m²", "1 cP = 0.01 P = 0.01 dyn s/cm² = 1e-3 Pa s", "cP_to_Pa_s", True),
        ("S28", "slip", "xylem paragraph", "the quoted ν with the quoted Re", "the quoted Re needs water's ν, 10× smaller", "xylem_reynolds", True),
        ("S29", "loose", "(16.5) vs the text", "the quoted τ_y range at H = 0.45", "the correlation gives 0.062–0.113 dyn/cm² for C_F = 0.21–0.46", "yield_stress_correlation", False),
        ("S30", "loose", "(16.221) units", "ρgz added in CGS to ψ in bar", "ρgz/10⁶", "water_potential_book_units", False),
        ("S31", "loose", "(16.37) turbulent", "the turbulent constant of (16.37)", "other sources use 1.359 Re^(1/4) (Çengel) — a different correlation", "entrance_length", False),
        ("S32", "loose", "§16.2 text", "a cardiac output written in mL/m", "mL/min", "cardiac_output", False),
        ("S33", "loose", "Exercise 16.4", "Shapiro (1997a)", "1977a", "—", False),
        ("S34", "slip", "(16.174)", "(2/κ)∇⁴ψ − (1/r)(∂ψ/∂r ∂/∂α − ∂ψ/∂α ∂/∂r)∇²ψ = −2w(…)",
         "+ (1/r)(…)∇²ψ (eliminating p₁ from (16.167)–(16.168) with (16.172)); changes only the K⁴ flux coefficient",
         "dean_series_sympy(printed_16_174=True), DEAN_C4_DERIVED", True),
    ]
    return pd.DataFrame(rows, columns=["id", "kind", "where", "printed", "correct", "code", "changes_physics"])


def traps():
    """The chapter's traps: places where the page is right and the reader goes wrong (analysis §9, T1–T27).

    Book: §16.2–§16.4.  Returns pandas.DataFrame: id, trap, guard.  Label: V7 (structure).
    """
    import pandas as pd

    rows = [
        ("T1", "two time conventions: e^{+iωt} (rigid Womersley) vs e^{i(kx − ωt)} (elastic tube)", "docstrings; I32 vs I12"),
        ("T2", "overloaded symbols α, K, R, δ, ψ, c, T, M, S, P, X, κ, λ, E, H, h, η, ν̂", "distinct code names"),
        ("T3", "units: SI, CGS, mmHg, bar, PRU", "converters MMHG, PRU, cP_to_Pa_s, bar_to_Pa"),
        ("T4", "G = −dp/dx (amplitude)", "G keyword as ch08"),
        ("T5", "τ = −μ du/dr > 0 here vs signed τ in ch08", "poiseuille_wall_stress returns the magnitude"),
        ("T6", "Re, St on the diameter; α² = Re·St/4; other α² = 2π Re St", "womersley_from_re_st"),
        ("T7", "two roots of (16.136)", "elastic_tube_roots by continuation"),
        ("T8", "Dean conventions κ, K = 2κ², D = 4√2κ, De", "dean_convert"),
        ("T9", "P < 0 in collapse; sign of P dK_p/dx; M = 2 − n", "collapsible_terms, shapiro_M"),
        ("T10", "S plays M's role: form the driver, divide by 1 − S² last", "collapsible_terms; terminal event"),
        ("T11", "PV loop counter-clockwise: ∮p dV < 0", "pv_loop_work"),
        ("T12", "Windkessel pressures transmural, p̄_V = 0", "windkessel_rhs"),
        ("T13", "hematocrit fraction vs percent", "yield_stress_correlation raises for H > 1"),
        ("T14", "Casson k_c (dimensionless) vs K_c = k_c√μ_p", "casson_constants"),
        ("T15", "the ⅓ rule is an estimate", "windkessel_periodic_summary rule_error"),
        ("T16", "tube-law singularity a₀(p − p_e)/Eh → 1", "tube_radius_hookean raises"),
        ("T17", "Moens–Korteweg radius vs diameter forms", "moens_korteweg_speed"),
        ("T18", "'closed end' means R = +1", "junction_coefficients"),
        ("T19", "the FL closure is ours; δ is free", "fahraeus_lindqvist_model"),
        ("T20", "Fig. 16.14 axis wall → centre, mirrored", "full-diameter figures"),
        ("T21", "Casson plug a_c ≥ a: no flow", "casson_pipe_flow_rate returns 0"),
        ("T22", "water potential < 0 under tension", "water_potential"),
        ("T23", "rigid-tube Womersley is not a wave", "docstrings"),
        ("T24", "'ours' demos labelled with cited inputs", "illustrative_inputs"),
        ("T25", "large-α Bessel overflow", "jve-scaled ratios"),
        ("T26", "equations stated without derivation", "curation §4b/§4c"),
        ("T27", "apparent viscosity: secant vs tangent", "apparent_viscosity, differential_viscosity"),
    ]
    return pd.DataFrame(rows, columns=["id", "trap", "guard"])


def conventions_table():
    """Symbols that change meaning in this chapter and the code names that keep them apart (trap T2).

    Book: §16.2–§16.4.  Returns pandas.DataFrame: symbol, meaning, code_name, where.  Label: V7 (structure).
    """
    import pandas as pd

    rows = [
        ("α", "Womersley number a√(ω/ν)", "alpha", "(16.12), (16.77)–(16.100)"),
        ("α", "area ratio A/A₀ (Shapiro)", "area_ratio", "(16.186)–(16.206)"),
        ("α", "toroidal angle (0 = outer bend)", "alpha_t", "(16.166)–(16.176)"),
        ("K", "compliance dv/dp", "K (Windkessel functions)", "(16.48)–(16.52)"),
        ("K", "metabolic cost per volume", "b_metabolic", "(16.138)–(16.143)"),
        ("K", "Dean's number 2κ²", "K (Dean functions)", "(16.178)"),
        ("R", "peripheral resistance", "R (Windkessel functions)", "(16.9), (16.49)"),
        ("R", "reflection coefficient", "refl / R in junction_coefficients", "(16.153)"),
        ("R", "radius of curvature", "R_curv", "(16.164)"),
        ("R", "viscous resistance per unit length", "friction_R_per_length", "(16.204)"),
        ("δ", "cell-free layer thickness", "delta", "(16.21)–(16.34)"),
        ("δ", "curvature ratio a/R", "delta (Dean functions)", "(16.164)"),
        ("ν, ν̂", "kinematic viscosity; Poisson ratio", "nu; poisson", "(16.12); (16.104)"),
        ("M", "Shapiro's coefficient (not Mach)", "shapiro_M", "(16.206)"),
        ("S", "speed index u/c", "speed_index", "(16.195)"),
        ("X", "Re(1 − F₂); also k²c₀²/ω²", "womersley_XY; Xk", "(16.97); (16.136)"),
        ("H", "hematocrit (fraction)", "H, H_c, H_F", "(16.5), (16.8)"),
        ("e^{iωt} / e^{−iωt}", "time convention", "I32 (rigid) / I12 (elastic)", "(16.73) / (16.116)"),
    ]
    return pd.DataFrame(rows, columns=["symbol", "meaning", "code_name", "where"])


def units_table():
    """Unit conversions used in the chapter (all exact except where noted; trap T3, slip S27).

    Book: §16.2 (mmHg, cP, PRU), §16.4 (bar, R).  Returns pandas.DataFrame.  Label: V7 (structure).
    """
    import pandas as pd

    rows = [
        ("1 mmHg", MMHG, "Pa", "conventional definition"),
        ("1 PRU (mmHg s/cm³)", PRU, "Pa s/m³", "MMHG/1e-6"),
        ("1 cP", 1e-3, "Pa s", "slip S27 corrected"),
        ("1 P = 1 dyn s/cm²", 0.1, "Pa s", ""),
        ("1 dyn/cm²", 0.1, "Pa", ""),
        ("1 bar = 10⁶ dyn/cm²", BAR, "Pa", ""),
        ("R", R_MOLAR, "J/(mol K)", "= 83.1446 cm³ bar/(mol K); the book's 83.141"),
    ]
    return pd.DataFrame(rows, columns=["unit", "value_SI", "SI_unit", "note"])
