"""Verification suite for Chapter 16 — Introduction to Biofluid Mechanics (Kundu, Cohen & Dowling 5e, §§16.1–16.4,
Eqs. (16.1)–(16.222)).

Evidence levels (``verify-implementation`` skill): V1 analytic · V2 symbolic (sympy / dimensions) · V3 convergence or an
independent numerical route · V4 conservation / invariant · V5 published benchmark (``reference/ch16/benchmarks.json``,
see SOURCES.md) · V6 book value (private, ``tests/book_values_ch16.json``, skipped when absent) · V7 limits / symmetry /
guards.  Every test name is ``test_<concept>_<level>_<what>``; the comment on the ``def`` line repeats the level and the
curation ID.

A items (CORE, >= 2 independent levels): C01 PV loop work · C02 Casson law (16.6) · C03 Womersley number (16.12) ·
C04 resistance (16.20) · C05 two-layer apparent viscosity (16.33) · C06 tube law (16.41) · C07 Windkessel (16.49) ·
C08 Moens–Korteweg (16.62) · C09 Womersley profile (16.82) · C10 amplitude and phase (16.99)–(16.100) · C11 viscous
elastic-tube waves (16.136) · C12 Murray (16.140) · C13 junction reflection (16.153) · C14 Dean (16.177)–(16.178) ·
C15 Shapiro's equations (16.204)–(16.206) · C16 Casson pipe flow (16.215)–(16.220) · C17 water potential (16.221).

Independent routes written here (not in ``fluidpy``): a Chebyshev–Fourier Newton solver of the Dean equations in the
Navier–Stokes form of the Imperial College M4A33 notes (reference/ch16/SOURCES.md); a 4 x 4 Bessel determinant of the
un-eliminated elastic-tube system (16.122)–(16.125), (16.128)–(16.129); the closed-form Starling-resistor choke
(Bernoulli with friction integrated in the area ratio); sympy re-derivations of the ★★/★★★ rows of curation §4b.

Printed slips kept as options that must FAIL an independent check: (16.30) ``printed_16_30``, (16.54) ``printed_16_54``,
(16.96) ``printed_16_96``, (16.140) ``printed_16_140``, (16.201) ``printed_16_201``, (16.174) ``printed_16_174``
(slip S34, found here in loop 1).

All inputs are ours (``ch16.illustrative_inputs()`` or visibly odd numbers); nothing here is quoted from the book.

Run: ``.venv/Scripts/python.exe -m pytest tests/test_ch16.py -q -p no:cacheprovider``  (``-m "not slow"`` skips the
script runs).
"""
from __future__ import annotations

import json
import math
import re
import subprocess
import sys
import warnings
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import numpy as np
import pytest
import sympy as sp
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq, fsolve, minimize_scalar
from scipy.special import jv

from fluidpy import ch16_biofluid_mechanics as ch16
from fluidpy.core import laminar as LAM
from fluidpy.core.units import Q_
from tools.convergence import observed_order, pairwise_orders

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference" / "ch16"
BOOK = Path(__file__).resolve().parent / "book_values_ch16.json"
book_only = pytest.mark.skipif(not BOOK.exists(), reason="book values are private; see CLAUDE.md rule 9")
slow = pytest.mark.slow
PI = math.pi
I = ch16.illustrative_inputs()
I32 = complex(np.exp(0.75j * np.pi))
I12 = complex(np.exp(0.25j * np.pi))


def book():
    return json.loads(BOOK.read_text(encoding="utf-8"))


def bench():
    return json.loads((REF / "benchmarks.json").read_text(encoding="utf-8"))


def rel(a, b):
    return abs(complex(a) - complex(b)) / max(abs(complex(b)), 1e-300)


def close(a, b, rtol=1e-10, atol=0.0):
    """``numpy.allclose`` with NO absolute tolerance unless one is written."""
    return bool(np.allclose(a, b, rtol=rtol, atol=atol))


# ======================================================================================================================
# Independent helper 1: Dean equations, Navier–Stokes form (Imperial College M4A33 notes, eq. (9)), spectral Newton
#   K J(psi, v) = 1 + lap v ;  K J(psi, Om) = lap Om - 2 K v v_z ;  Om = -lap psi ;  v = psi = psi_r = 0 at r = 1,
#   J(f, g) = f_z g_x - f_x g_z = (1/r)(f_th g_r - f_r g_th).   Straight pipe v = (1 - r^2)/4; book's K = K_lec^2 / 8.
#   ``sJ = -1`` flips the vorticity-advection sign (the form printed in (16.174), see test_dean_vorticity_V2_...).
# ======================================================================================================================
def _cheb(N):
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2, np.ones(N - 1), 2]) * (-1) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    D = np.outer(c, 1 / c) / (X - X.T + np.eye(N + 1))
    return D - np.diag(D.sum(axis=1)), x


def _fourier(M):
    h = 2 * np.pi / M
    k = np.arange(1, M)
    col = np.zeros(M)
    col[1:] = 0.5 * (-1) ** k / np.tan(k * h / 2)
    col2 = np.zeros(M)
    col2[0] = -np.pi ** 2 / (3 * h ** 2) - 1 / 6
    col2[1:] = -0.5 * (-1) ** k / np.sin(h * k / 2) ** 2
    D1 = np.array([np.roll(col, i) for i in range(M)]).T     # D1[j, i] = col[(j - i) mod M]
    D2 = np.array([np.roll(col2, i) for i in range(M)]).T
    return D1, D2, h * np.arange(M)


def _dean_ops(N=21, M=16):
    D, x = _cheb(N)
    D2 = D @ D
    n2 = (N + 1) // 2
    r = x[:n2]
    Dt, Dtt, th = _fourier(M)
    S = np.roll(np.eye(M), M // 2, axis=1)                   # theta -> theta + pi
    Dr = np.kron(D[:n2, :n2], np.eye(M)) + np.kron(D[:n2, n2:][:, ::-1], S)
    Drr = np.kron(D2[:n2, :n2], np.eye(M)) + np.kron(D2[:n2, n2:][:, ::-1], S)
    Dth = np.kron(np.eye(n2), Dt)
    rinv = np.kron(1 / r, np.ones(M))
    Lap = Drr + rinv[:, None] * Dr + (rinv ** 2)[:, None] * np.kron(np.eye(n2), Dtt)
    RR, TT = np.meshgrid(r, th, indexing="ij")
    return dict(r=r, Dr=Dr, Dth=Dth, Lap=Lap, rinv=rinv, RR=RR.ravel(), TT=TT.ravel(), n2=n2, M=M)


def _dean_solve(K, op, U0=None, sJ=1.0, tol=1e-13):
    n = op["n2"] * op["M"]
    Dr, Dth, Lap, ri = op["Dr"], op["Dth"], op["Lap"], op["rinv"]
    s, c = np.sin(op["TT"]), np.cos(op["TT"])
    bnd = np.isclose(op["RR"], 1.0)
    U = np.concatenate([(1 - op["RR"] ** 2) / 4, np.zeros(2 * n)]) if U0 is None else U0.copy()
    Id = np.eye(n)
    for _ in range(40):
        v, psi, Om = U[:n], U[n:2 * n], U[2 * n:]
        vr, vt, pr, pt, Or, Ot = Dr @ v, Dth @ v, Dr @ psi, Dth @ psi, Dr @ Om, Dth @ Om
        vz = s * vr + c * ri * vt
        F1 = K * ri * (pt * vr - pr * vt) - 1 - Lap @ v
        F2 = Om + Lap @ psi
        F3 = sJ * K * ri * (pt * Or - pr * Ot) - Lap @ Om + 2 * K * v * vz
        J11 = K * ri[:, None] * (pt[:, None] * Dr - pr[:, None] * Dth) - Lap
        J12 = K * ri[:, None] * (vr[:, None] * Dth - vt[:, None] * Dr)
        J13 = np.zeros((n, n))
        J21, J22, J23 = np.zeros((n, n)), Lap.copy(), Id.copy()
        J31 = 2 * K * (np.diag(vz) + v[:, None] * (s[:, None] * Dr + (c * ri)[:, None] * Dth))
        J32 = sJ * K * ri[:, None] * (Or[:, None] * Dth - Ot[:, None] * Dr)
        J33 = sJ * K * ri[:, None] * (pt[:, None] * Dr - pr[:, None] * Dth) - Lap
        F1[bnd], J11[bnd], J12[bnd], J13[bnd] = v[bnd], Id[bnd], 0, 0
        F2[bnd], J21[bnd], J22[bnd], J23[bnd] = psi[bnd], 0, Id[bnd], 0
        F3[bnd], J31[bnd], J32[bnd], J33[bnd] = pr[bnd], 0, Dr[bnd], 0
        dU = np.linalg.solve(np.block([[J11, J12, J13], [J21, J22, J23], [J31, J32, J33]]),
                             -np.concatenate([F1, F2, F3]))
        U = U + dU
        if np.max(np.abs(dU)) < tol:
            return U
    raise RuntimeError("Dean Newton did not converge")


def _dean_flux_ratio_numeric(U, op):
    """(mean of v over the disk)/(1/8): theta mean exact; radial integral by interpolation in s = 2r^2 - 1."""
    from numpy.polynomial import chebyshev as C
    n2, M = op["n2"], op["M"]
    vbar = U[:n2 * M].reshape(n2, M).mean(axis=1)
    ci = C.chebint(C.chebfit(2 * op["r"] ** 2 - 1, vbar, n2 - 1))
    return 8.0 * 0.5 * (C.chebval(1.0, ci) - C.chebval(-1.0, ci))


@lru_cache(maxsize=4)
def _dean_flux_samples(sJ=1.0, N=21, M=16):
    """Flux ratios of the spectral solutions at 25 Dean numbers, x = K/576 in [0.05, 0.35] (K = K_lec^2/8)."""
    op = _dean_ops(N, M)
    xs = np.linspace(0.05, 0.35, 25)
    F, U = [], None
    for x in xs:
        U = _dean_solve(math.sqrt(8 * 576 * x), op, U0=U, sJ=sJ)
        F.append(_dean_flux_ratio_numeric(U, op))
    return xs, np.array(F)


def _dean_fit(xs, F, nterms=5):
    """Least-squares fit of F - 1 = -c2 x^2 + c4 x^4 + ... with ``nterms`` even powers; returns (c2, c4)."""
    A = np.vstack([xs ** (2 * k) for k in range(1, nterms + 1)]).T
    coef, *_ = np.linalg.lstsq(A, F - 1, rcond=None)
    return -coef[0], coef[1]


@lru_cache(maxsize=4)
def _dean_coefficients_numeric(sJ=1.0, N=21, M=16):
    """Fit Q/(pi a^2 W) = 1 - c2 x^2 + c4 x^4 - ... (x = K/576, K = K_lec^2/8) to the spectral solutions."""
    return _dean_fit(*_dean_flux_samples(sJ, N, M), nterms=5)


# ======================================================================================================================
# Independent helper 2: the un-eliminated long-wave system of the viscous elastic tube, unknowns (A, eta, xi, p),
# u(r) = A J0(beta r) + k p/(omega rho)  (particular solution of (16.123)); nondimensional a0 = rho = omega = 1, c0 = 1.
# ======================================================================================================================
def _elastic_det4(X, alpha, nu_hat=0.5):
    a0 = rho = omega = 1.0
    nu = omega * a0 ** 2 / alpha ** 2
    mu, Eh = rho * nu, 2.0
    k = np.sqrt(complex(X))
    beta = np.sqrt(1j * omega / nu)
    Cp = k / (omega * rho)
    J0, J1 = jv(0, beta * a0), jv(1, beta * a0)
    du_A = -beta * J1
    M = np.zeros((4, 4), complex)
    M[0] = [1j * nu_hat * mu / (a0 * k) * du_A, -Eh / a0 ** 2, 0, 1]                        # (16.128)
    M[1] = [(1 - nu_hat ** 2) * mu / (Eh * k ** 2) * du_A, -1j * nu_hat / (k * a0), 1, 0]   # (16.129)
    M[2] = [J0, 0, 1j * omega, Cp]                                                          # (16.124)
    M[3] = [-(1j * k / a0) * a0 * J1 / beta, 1j * omega, 0, -(1j * k / a0) * Cp * a0 ** 2 / 2]  # (16.125), (16.122)
    M = M / np.abs(M).max(axis=1, keepdims=True)
    return np.linalg.det(M)


def _root4(X0, alpha, nu_hat=0.5):
    f = lambda v: [_elastic_det4(complex(*v), alpha, nu_hat).real, _elastic_det4(complex(*v), alpha, nu_hat).imag]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        s = fsolve(f, [X0.real * 1.01, X0.imag * 0.99 + 1e-3], xtol=1e-14)
    return complex(*s)


# ======================================================================================================================
# Independent helper 3: Starling resistor, uniform tube, laminar friction, power-law tube law, open inlet (alpha = 1).
# Bernoulli with friction d(p + rho u^2/2) = -R' Q dx, R' = 8 pi mu/(alpha A0)^2, p - pe = Kp (1 - alpha^-n):
#   x(alpha) = G(alpha_in) - G(alpha),  G(a) = n Kp A0^2 a^(2-n)/((2-n) 8 pi mu Q) - rho Q ln(a)/(8 pi mu).
# ======================================================================================================================
def _starling_closed(L, A0, Kp, rho, mu, n=1.5):
    G = lambda a, Q: n * Kp * A0 ** 2 * a ** (2 - n) / ((2 - n) * 8 * PI * mu * Q) - rho * Q * math.log(a) / (8 * PI * mu)
    a_star = lambda Q: (rho * Q ** 2 / (A0 ** 2 * Kp * n)) ** (1 / (2 - n))
    Qhi = A0 * math.sqrt(n * Kp / rho)
    Qc = brentq(lambda Q: G(1.0, Q) - G(a_star(Q), Q) - L, Qhi * 1e-6, Qhi * (1 - 1e-12), xtol=1e-22, rtol=1e-15)
    return dict(Q=Qc, alpha=a_star(Qc), p2=Kp * (1 - a_star(Qc) ** (-n)), G=G, a_star=a_star)


VEIN = dict(L=0.1, A0=PI * 0.004 ** 2, Kp=200.0, n=1.5, rho=1060.0, mu=3.5e-3)     # the design's vein (ours)


# ======================================================================================================================
# C01 — the PV loop and the work per beat, W = −∮p dV (ours); N04–N14 bookkeeping
# ======================================================================================================================
def _schematic_kwargs():
    return dict(p_min=I["p_min"], p_ED=I["p_ED"], p_open=I["p_open"], p_peak=I["p_peak"], p_ES=I["p_ES"])


def test_pv_loop_V1_rectangle_and_orientation():  # V1 — C01, D01 (trap T11)
    for dp, dV in ((7.3e3, 41e-6), (13.1e3, 66e-6)):
        V = np.array([50e-6, 50e-6 + dV, 50e-6 + dV, 50e-6])
        p = np.array([1e3, 1e3, 1e3 + dp, 1e3 + dp])           # counterclockwise in (V, p)
        w = ch16.pv_loop_work(p, V)
        assert rel(w["work"], dp * dV) < 1e-13 and rel(w["signed_area"], dp * dV) < 1e-13
        assert w["orientation"] == "counterclockwise"
        wr = ch16.pv_loop_work(p[::-1], V[::-1])                 # traversed the other way: only the sign changes
        assert rel(wr["work"], -dp * dV) < 1e-13 and wr["orientation"] == "clockwise"
        wc = ch16.pv_loop_work(np.append(p, p[0]), np.append(V, V[0]))   # closed polygon gives the same
        assert rel(wc["work"], dp * dV) < 1e-13


def test_pv_loop_V2_exact_work_from_the_leg_shapes():  # V2 — C01, D01: sympy ∮ over the four legs of the schematic
    EDV, ESV = I["EDV"], I["ESV"]
    kw = _schematic_kwargs()
    V, s = sp.symbols("V s", real=True)
    dV = sp.nsimplify(EDV - ESV)
    pmin, pED, popen, ppk, pES = (sp.nsimplify(kw[k]) for k in ("p_min", "p_ED", "p_open", "p_peak", "p_ES"))
    fill = pmin + (pED - pmin) * s ** 2                                    # A → B, V = ESV + s ΔV
    q = s
    eject = popen * (1 - q) * (1 - 2 * q) + 4 * ppk * q * (1 - q) + pES * q * (2 * q - 1)   # C → D, V = EDV − q ΔV
    oint = sp.integrate(fill * dV, (s, 0, 1)) + sp.integrate(eject * (-dV), (q, 0, 1))    # isovolumic legs: dV = 0
    W_exact = float(-oint)
    assert rel(ch16.pv_loop_schematic_work(EDV, ESV, **kw), W_exact) < 1e-13
    lp = ch16.pv_loop_schematic(EDV, ESV, **kw, n=4000)
    assert rel(ch16.pv_loop_work(lp["p"], lp["V"])["work"], W_exact) < 1e-7
    assert W_exact > 0 and ch16.pv_loop_work(lp["p"], lp["V"])["orientation"] == "counterclockwise"
    c = lp["corners"]
    assert c["A"] == (ESV, kw["p_min"]) and c["B"] == (EDV, kw["p_ED"]) and c["C"] == (EDV, kw["p_open"])


def test_pv_loop_V3_polygon_converges_second_order():  # V3 — C01: trapezoid along the polygon, error ∝ n⁻²
    EDV, ESV = I["EDV"], I["ESV"]
    kw = _schematic_kwargs()
    W = ch16.pv_loop_schematic_work(EDV, ESV, **kw)
    ns = np.array([25, 50, 100, 200, 400])
    err = [abs(ch16.pv_loop_work(**{k: v for k, v in ch16.pv_loop_schematic(EDV, ESV, **kw, n=int(n)).items()
                                     if k in ("p", "V")})["work"] - W) for n in ns]
    assert abs(observed_order(1.0 / ns, err) - 2.0) < 0.15


def test_circulation_bookkeeping_V1_hand_values():  # V1 — N06, N07, N09, N12–N14, N08, N10, N11
    SV = ch16.stroke_volume(131e-6, 53e-6)
    assert rel(SV, 78e-6) < 1e-12
    assert rel(ch16.ejection_fraction(SV, 131e-6), 78 / 131) < 1e-12
    assert rel(ch16.ejection_fraction(SV, 131e-6, percent=True), 7800 / 131) < 1e-12
    CO = ch16.cardiac_output(SV, 63.0)
    assert rel(CO, 78e-6 * 63 / 60) < 1e-14                               # m³/s
    assert rel(ch16.cardiac_index(CO, 1.83), CO / 1.83) < 1e-14
    assert rel(ch16.cardiac_power(1.17, 63.0), 1.17 * 63 / 60) < 1e-14
    assert ch16.pulse_pressure(16.1e3, 10.3e3) == pytest.approx(5.8e3, rel=1e-13)
    assert rel(ch16.mean_arterial_pressure(16.1e3, 10.3e3), (16.1e3 + 2 * 10.3e3) / 3) < 1e-14   # (16.11)
    R = ch16.peripheral_resistance(12.4e3, CO, pV_mean=0.6e3)
    assert rel(R, (12.4e3 - 0.6e3) / CO) < 1e-14                           # (16.9)–(16.10)
    assert ch16.MMHG == 133.322387415 and rel(ch16.PRU, 133.322387415e6) < 1e-15
    assert rel(ch16.mmHg_to_Pa(97.0), 97 * 133.322387415) < 1e-15 and rel(ch16.Pa_to_mmHg(ch16.mmHg_to_Pa(97.0)), 97) < 1e-15
    assert ch16.cP_to_Pa_s(3.3) == pytest.approx(3.3e-3, rel=1e-15)        # slip S27 corrected
    assert ch16.poise_to_Pa_s(0.033) == pytest.approx(3.3e-3, rel=1e-15)
    assert ch16.dyn_cm2_to_Pa(17.0) == pytest.approx(1.7, rel=1e-15) and ch16.bar_to_Pa(0.37) == pytest.approx(3.7e4, rel=1e-15)
    assert ch16.mL_to_m3(78.0) == pytest.approx(78e-6, rel=1e-15)
    p = ch16.series_pressure_drop([2e7, 9e7, 3e7], 9e-5, 13.2e3)
    assert close(p, [13.2e3 - 1.8e3, 13.2e3 - 9.9e3, 13.2e3 - 12.6e3], rtol=1e-13)
    assert rel(ch16.mean_speed(9e-5, 0.3), 3e-4) < 1e-14
    assert rel(ch16.equal_volume_sphere_area(4 / 3 * PI * 2.7e-6 ** 3), 4 * PI * 2.7e-6 ** 2) < 1e-13
    with pytest.raises(ValueError):
        ch16._check_hematocrit(45.0)                                        # trap T13: a percent passed
    with pytest.raises(ValueError):
        ch16._check_hematocrit(0.0)
    assert ch16._check_hematocrit(0.41) == 0.41


def test_circulation_V2_units_are_consistent():  # V2 — N07, N08, N12 (pint)
    CO = Q_(78, "mL") * Q_(63, "1/min")
    assert rel(CO.to("m**3/s").magnitude, ch16.cardiac_output(78e-6, 63.0)) < 1e-12
    assert rel(Q_(1, "mmHg").to("Pa").magnitude, ch16.MMHG) < 1e-9
    assert rel((Q_(1, "mmHg") / Q_(1, "cm**3/s")).to("Pa*s/m**3").magnitude, ch16.PRU) < 1e-9
    assert rel(Q_(1, "cP").to("Pa*s").magnitude, ch16.cP_to_Pa_s(1.0)) < 1e-12
    assert rel(Q_(1, "dyn/cm**2").to("Pa").magnitude, ch16.dyn_cm2_to_Pa(1.0)) < 1e-12


def test_time_mean_V1_and_one_third_rule_estimate():  # V1 — N14 (trap T15)
    t = np.linspace(0, 0.83, 4001)
    assert rel(ch16.time_mean(3.7 + 0 * t, t), 3.7) < 1e-14
    assert abs(ch16.time_mean(np.sin(2 * PI * t / 0.83), t)) < 1e-12
    # a waveform that sits at pS for a third of the period and pD for two thirds has exactly the ⅓-rule mean
    pS, pD = 16.1e3, 10.3e3
    tt = np.linspace(0, 3.0, 300001)
    p = np.where(tt < 1.0, pS, pD)
    assert rel(ch16.time_mean(p, tt), ch16.mean_arterial_pressure(pS, pD)) < 1e-5


# ======================================================================================================================
# C02 — the Casson law (16.6)/(16.207); N19–N22 (yield stress, apparent and differential viscosity, Whitmore)
# ======================================================================================================================
def test_casson_V1_stress_and_shear_rate_are_inverses():  # V1 — C02
    ty, Kc = I["tau_y"], I["Kc"]
    g = np.geomspace(1e-3, 1e4, 300)
    tau = ch16.casson_stress(g, ty, Kc)
    assert close((np.sqrt(tau) - np.sqrt(ty)) / np.sqrt(g), Kc, rtol=1e-12)     # √τ − √τ_y = K_c√γ̇ (16.207)
    assert close(ch16.casson_shear_rate(tau, ty, Kc), g, rtol=1e-10)
    assert np.all(ch16.casson_shear_rate(np.linspace(0, ty * 0.999, 7), ty, Kc) == 0.0)   # below yield: no shear
    assert close(ch16.casson_apparent_viscosity(g, ty, Kc), tau / g, rtol=1e-12)          # secant (N21)
    assert close(ch16.apparent_viscosity(tau, g), tau / g, rtol=1e-15)
    h = g * 1e-6                                                                          # tangent (trap T27)
    fd = (ch16.casson_stress(g + h, ty, Kc) - ch16.casson_stress(g - h, ty, Kc)) / (2 * h)
    assert close(ch16.differential_viscosity(g, ty, Kc), fd, rtol=1e-7)
    assert np.all(ch16.differential_viscosity(g, ty, Kc) < ch16.casson_apparent_viscosity(g, ty, Kc))
    assert close(ch16.newtonian_stress(g, 3.3e-3), 3.3e-3 * g, rtol=1e-15)


def test_casson_V2_derivation():  # V2 — D02: square the law, secant slope, the two limits, K_c = k_c√μ_p
    tau, ty, Kc, g, kc, mup = sp.symbols("tau tau_y K_c gamma k_c mu_p", positive=True)
    tau_g = (sp.sqrt(ty) + Kc * sp.sqrt(g)) ** 2                                  # step 1: square (16.207)
    mu_app = sp.simplify(tau_g / g)                                                # step 2: secant
    assert sp.simplify(mu_app - (sp.sqrt(ty / g) + Kc) ** 2) == 0                  # step 3: D02's form
    assert sp.limit(mu_app, g, sp.oo) == Kc ** 2                                   # step 4: high-shear limit
    assert sp.limit(mu_app, g, 0, "+") == sp.oo                                    # step 4: diverges at low shear
    # step 5: (16.6) √(τ/μ_p) = k_c√γ̇ + √(τ_y/μ_p) is (16.207) with K_c = k_c√μ_p
    lhs = sp.sqrt(tau_g.subs(Kc, kc * sp.sqrt(mup)) / mup)
    assert sp.simplify(lhs - (kc * sp.sqrt(g) + sp.sqrt(ty / mup))) == 0
    tyn, Kcn = ch16.casson_constants(1.41, 1.9, 1.3e-3)                          # τ_y = (√(τ_y/μ_p))²μ_p, K_c = k_c√μ_p
    assert rel(tyn, 1.9 ** 2 * 1.3e-3) < 1e-15 and rel(Kcn, 1.41 * math.sqrt(1.3e-3)) < 1e-15
    assert rel(ch16.casson_stress_plasma_form(37.0, 1.41, 1.9, 1.3e-3), float(tau_g.subs({ty: tyn, Kc: Kcn, g: 37.0}))) < 1e-13
    # dimensions: K_c in (Pa s)^½, τ in Pa
    Kq = Q_(0.055, "Pa**0.5*s**0.5")
    assert ((Q_(5e-3, "Pa") ** 0.5 + Kq * Q_(100, "1/s") ** 0.5) ** 2).check("[pressure]")


def test_casson_V7_limits():  # V7 — C02: τ_y → 0 Newtonian with μ = K_c²; μ_app − K_c² ∝ γ̇^(−1/2)
    Kc = 0.061
    g = np.geomspace(0.1, 1e3, 50)
    assert close(ch16.casson_stress(g, 0.0, Kc), Kc ** 2 * g, rtol=1e-14)
    gg = np.geomspace(1e6, 1e10, 9)
    ex = ch16.casson_apparent_viscosity(gg, 4e-3, Kc) - Kc ** 2
    slope = np.polyfit(np.log(gg), np.log(ex), 1)[0]
    assert abs(slope + 0.5) < 1e-3


def test_whitmore_and_yield_stress_V1():  # V1 — N22 (16.7), N19 (16.5), C02 constants (trap T14)
    mup = 1.27e-3
    w = ch16.whitmore_casson(np.array([0.5, 50.0, 800.0]), mup)
    assert rel(w["tau_y"], 4.0 * mup) < 1e-14 and rel(w["Kc"], 1.53 * math.sqrt(mup)) < 1e-14
    assert close(np.sqrt(w["tau"] / mup), 1.53 * np.sqrt([0.5, 50.0, 800.0]) + 2.0, rtol=1e-13)   # (16.7)
    assert rel(ch16.casson_apparent_viscosity(1e12, w["tau_y"], w["Kc"]) / mup, 1.53 ** 2) < 1e-5  # μ∞/μ_p = k_c²
    ty, Kc = ch16.casson_constants(1.53, 2.0, mup)
    assert rel(ty, w["tau_y"]) < 1e-14 and rel(Kc, w["Kc"]) < 1e-14
    assert close(ch16.casson_stress_plasma_form(np.array([3.0, 30.0]), 1.53, 2.0, mup),
                 ch16.casson_stress(np.array([3.0, 30.0]), ty, Kc), rtol=1e-13)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tyc = ch16.yield_stress_correlation(0.41, 0.37)
    assert rel(tyc, ((0.41 - 0.1) * (0.37 + 0.5)) ** 2 * 0.1) < 1e-14       # (16.5) in dyn/cm² → Pa
    with pytest.raises(ValueError):
        ch16.yield_stress_correlation(41.0, 0.37)                           # percent passed (trap T13)


# ======================================================================================================================
# C03 — the Womersley number (16.12); D03
# ======================================================================================================================
def test_womersley_number_V1_three_routes():  # V1 — C03, D03
    a, HR, rho, mu = 3.1e-3, 71.0, 1060.0, 3.5e-3
    om, nu = 2 * PI * HR / 60, mu / rho
    al = ch16.womersley_number(a, om, nu)
    assert rel(al, a * math.sqrt(om / nu)) < 1e-15
    assert rel(ch16.stokes_number(a, om, nu), al ** 2) < 1e-14
    u = 0.23
    Re, St = 2 * a * u / nu, 2 * a * om / u                                # the chapter's diameter forms
    assert rel(ch16.womersley_from_re_st(Re, St), al) < 1e-14              # α² = Re·St/4
    de = LAM.stokes_layer(nu, om)["delta_e"]
    assert rel(al, math.sqrt(2) * a / de) < 1e-14                          # α = √2 a/δ_e
    assert (ch16.womersley_regime(0.3), ch16.womersley_regime(4.0), ch16.womersley_regime(17.0)) == (
        "quasi-steady", "transitional", "inertial")


def test_womersley_number_V2_derivation_and_dimensions():  # V2 — D03
    a, om, nu, u = sp.symbols("a omega nu u", positive=True)
    Re, St = 2 * a * u / nu, 2 * a * om / u
    assert sp.simplify(Re * St / 4 - a ** 2 * om / nu) == 0
    de = sp.sqrt(2 * nu / om)
    assert sp.simplify((sp.sqrt(2) * a / de) ** 2 - a ** 2 * om / nu) == 0
    q = Q_(3.1e-3, "m") * (Q_(7.4, "1/s") / Q_(3.3e-6, "m**2/s")) ** 0.5
    assert q.dimensionless


def test_vessel_scales_table_V1_rows_recomputed():  # V1 — N24, C03 figure
    df = ch16.vessel_scales_table(HR=66.0)
    om = 2 * PI * 66 / 60
    for _, row in df.iterrows():
        v = ch16.VESSELS[row["vessel"]]
        assert rel(row["alpha"], v["a"] * math.sqrt(om * 1060 / 3.5e-3)) < 1e-13
        assert rel(row["Re"], 1060 * v["u_mean"] * 2 * v["a"] / 3.5e-3) < 1e-13
        assert rel(row["h_over_a"], v["h"] / v["a"]) < 1e-14
    a = df.set_index("vessel")["alpha"]
    assert a["capillary"] < 0.1 < 10 < a["aorta"]                          # V7: the two regimes the chapter names
    assert set(ch16.BLOOD) >= {"rho", "mu", "mu_p", "T"} and ch16.BLOOD["mu"] > ch16.BLOOD["mu_p"]


# ======================================================================================================================
# C04 — vascular resistance (16.20) and the Poiseuille recap R04–R07; N27 entrance length
# ======================================================================================================================
def test_resistance_V1_poiseuille_quadrature_and_ch08():  # V1 — C04, R04–R07
    a, dp, L, mu = 1.37e-3, 211.0, 0.071, 3.3e-3
    r = np.linspace(0, a, 101)
    u = ch16.poiseuille_velocity(r, a, dp, L, mu)
    assert close(u, a ** 2 * dp / (4 * mu * L) * (1 - (r / a) ** 2), rtol=1e-13, atol=1e-18)       # (16.15)
    Qq, _ = quad(lambda s: 2 * PI * s * float(ch16.poiseuille_velocity(s, a, dp, L, mu)), 0, a, epsabs=0, epsrel=1e-13)
    Q = ch16.poiseuille_flow_rate(a, dp, L, mu)
    assert rel(Q, Qq) < 1e-12 and rel(Q, PI * a ** 4 * dp / (8 * mu * L)) < 1e-14                     # (16.17)
    assert rel(ch16.vascular_resistance(mu, L, a), dp / Q) < 1e-13                                    # (16.20)
    assert rel(Q / (PI * a ** 2), 0.5 * u[0]) < 1e-13                                                 # (16.18)
    tau = ch16.poiseuille_wall_stress(a, dp, L)
    assert tau > 0 and rel(tau, a * dp / (2 * L)) < 1e-14                                             # (16.19), S1
    h = a * 1e-6                                                                                      # −μ du/dr at the wall
    assert rel(-mu * (float(ch16.poiseuille_velocity(a, a, dp, L, mu)) - float(ch16.poiseuille_velocity(a - h, a, dp, L, mu))) / h, tau) < 1e-5


def test_resistance_V4_force_balance_and_scaling():  # V4 — C04: τ_w·2πaL = πa²Δp; R a⁴ constant; networks
    a, dp, L = 2.3e-3, 151.0, 0.093
    assert rel(ch16.poiseuille_wall_stress(a, dp, L) * 2 * PI * a * L, PI * a ** 2 * dp) < 1e-14
    aa = np.geomspace(5e-6, 1e-2, 40)
    Ra = ch16.vascular_resistance(3.5e-3, 0.05, aa) * aa ** 4
    assert np.ptp(Ra) / Ra.mean() < 1e-14
    assert rel(ch16.vascular_resistance(3.5e-3, 0.05, 1e-3) / ch16.vascular_resistance(3.5e-3, 0.05, 2e-3), 16) < 1e-13
    R1 = ch16.vascular_resistance(3.5e-3, 0.002, 21e-6)
    assert rel(ch16.network_resistance([R1] * 37, arrangement="parallel"), R1 / 37) < 1e-13
    assert rel(ch16.network_resistance([R1, 2 * R1, 4 * R1], arrangement="series"), 7 * R1) < 1e-14
    with pytest.raises(ValueError):
        ch16.network_resistance([R1], arrangement="mixed")


def test_resistance_V2_derivation():  # V2 — D04: disc integral of the profile, R = Δp/Q, units
    r, a, mu, L, dp = sp.symbols("r a mu L Delta_p", positive=True)
    u = a ** 2 * dp / (4 * mu * L) * (1 - r ** 2 / a ** 2)
    assert sp.simplify(mu * (sp.diff(u, r, 2) + sp.diff(u, r) / r) + dp / L) == 0                     # momentum
    Q = sp.integrate(2 * sp.pi * r * u, (r, 0, a))
    assert sp.simplify(Q - sp.pi * a ** 4 * dp / (8 * mu * L)) == 0
    assert sp.simplify(dp / Q - 8 * mu * L / (sp.pi * a ** 4)) == 0
    Rq = 8 * Q_(3.5e-3, "Pa*s") * Q_(0.1, "m") / (PI * Q_(2e-3, "m") ** 4)
    assert Rq.check("[pressure]*[time]/[length]**3")


def test_entrance_length_V1_V5_correlations():  # V1 + V5 — N27 (16.37), slip S31
    d = 4.1e-3
    assert rel(ch16.entrance_length(733.0, d, regime="laminar"), 0.06 * 733 * d) < 1e-14
    assert rel(ch16.entrance_length(5200.0, d, regime="turbulent"), 0.693 * 5200 ** 0.25 * d) < 1e-14
    with pytest.raises(ValueError):
        ch16.entrance_length(30.0, d, regime="laminar")
    with pytest.warns(UserWarning):
        ch16.entrance_length(30.0, d, regime="laminar", allow_low_re=True)
    bm = bench()["entrance_length_wikipedia"]
    lam_ref = bm["laminar_coefficient"] * 733 * d
    assert abs(ch16.entrance_length(733.0, d, regime="laminar") / lam_ref - 1) < 0.05        # 0.06 vs 0.0575: 4.3 %
    turb_ratio = bm["turbulent_coefficient"] * 5200 ** 0.25 * d / ch16.entrance_length(5200.0, d, regime="turbulent")
    assert 1.9 < turb_ratio < 2.0                                       # S31: a different correlation, factor ≈ 1.96


# ======================================================================================================================
# C05 — the two-layer (Fåhræus–Lindqvist) model (16.21)–(16.36); N29–N36; benchmark Pries/Secomb–Pries
# ======================================================================================================================
FL = dict(a=19e-6, delta=2.3e-6, dp=37.0, L=1.3e-3, mu_p=1.27e-3, mu_c=4.1e-3)


def test_two_layer_V2_derivation_profiles_fluxes_viscosity():  # V2 — D05, D06, D07 (slips S2, S3)
    r, a, d, G, mp, mc = sp.symbols("r a delta G mu_p mu_c", positive=True)     # G = Δp/L
    C1, C2 = sp.symbols("C1 C2")
    b = a - d
    # D05: −μ du/dr = (r/2)Δp/L in each layer (sign corrected, S2); regular axis; continuity at b; no slip at a
    uc = -G * r ** 2 / (4 * mc) + C1
    up = -G * r ** 2 / (4 * mp) + C2
    sol = sp.solve([up.subs(r, a), (uc - up).subs(r, b)], [C1, C2])
    uc, up = sp.simplify(uc.subs(sol)), sp.simplify(up.subs(sol))
    assert sp.simplify(-mc * sp.diff(uc, r) - (-mp * sp.diff(up, r))) == 0                 # (16.23) stress continuous
    page_27 = a ** 2 / (4 * mp) * G * (1 - (r / a) ** 2)                                    # (16.27) as printed
    page_28 = a ** 2 / (4 * mp) * G * (1 - (b / a) ** 2 - mp / mc * (r / a) ** 2 + mp / mc * (b / a) ** 2)   # (16.28)
    assert sp.simplify(up - page_27) == 0 and sp.simplify(uc - page_28) == 0
    # D06: the layer fluxes
    Qp = sp.integrate(2 * sp.pi * r * up, (r, b, a))
    Qc = sp.integrate(2 * sp.pi * r * uc, (r, 0, b))
    assert sp.simplify(Qp - sp.pi * G / (8 * mp) * (a ** 2 - b ** 2) ** 2) == 0                              # (16.29)
    corr_30 = sp.pi * b ** 2 * G / (4 * mp) * (a ** 2 - (1 - mp / (2 * mc)) * b ** 2)
    printed_30 = sp.pi * a ** 2 * G / (4 * mp) * (a ** 2 - (1 - mp / (2 * mc)) * b ** 4 / a ** 2)
    assert sp.simplify(Qc - corr_30) == 0                                                                    # (16.30) corrected
    # S3's size: printed − correct = +πΔp a²δ(2a − δ)/(4μ_pL) (> 0; analysis §9 S3 writes δ(δ − 2a), the opposite sign)
    assert sp.simplify(printed_30 - corr_30 - sp.pi * G * a ** 2 * d * (2 * a - d) / (4 * mp)) == 0
    total_31 = sp.pi * a ** 4 * G / (8 * mp) * (1 - (1 - d / a) ** 4 * (1 - mp / mc))
    assert sp.simplify(Qp + Qc - total_31) == 0                                                              # (16.31)
    assert sp.simplify(Qp + printed_30 - total_31) != 0                                                       # printed fails
    # D07: μ_app by Poiseuille analogy (16.32) → (16.33) → thin-layer (16.34)
    mu_app = sp.simplify(sp.pi * a ** 4 * G / (8 * total_31))
    assert sp.simplify(mu_app - mp / (1 - (1 - d / a) ** 4 * (1 - mp / mc))) == 0
    e = sp.symbols("e", positive=True)                                                                      # e = δ/a
    thin = mc / (1 + 4 * e * (mc / mp - 1))
    full = mp / (1 - (1 - e) ** 4 * (1 - mp / mc))
    assert sp.simplify(sp.series(full - thin, e, 0, 2).removeO()) == 0                                       # same to O(e)
    # and the coded functions are these expressions
    num = {a: FL["a"], d: FL["delta"], G: FL["dp"] / FL["L"], mp: FL["mu_p"], mc: FL["mu_c"]}
    q = ch16.two_layer_flow_rates(FL["a"], FL["delta"], FL["dp"], FL["L"], FL["mu_p"], FL["mu_c"])
    assert rel(q["Qp"], float(Qp.subs(num))) < 1e-12 and rel(q["Qc"], float(corr_30.subs(num))) < 1e-12
    assert rel(ch16.apparent_viscosity_two_layer(FL["a"], FL["delta"], FL["mu_p"], FL["mu_c"]), float(mu_app.subs(num))) < 1e-12
    for frac in (0.0, 0.3, 0.7, 0.95, 1.0):
        rr = frac * FL["a"]
        ref = (uc if rr <= FL["a"] - FL["delta"] else up).subs(num).subs(r, rr)
        assert abs(float(ch16.two_layer_velocity(rr, *FL.values())) - float(ref)) <= 1e-13 * float(uc.subs(num).subs(r, 0))


def test_two_layer_V1_quadrature_of_the_profile_and_slip_16_30():  # V1 — C05 (slip S3 printed_16_30 must fail)
    a, d, dp, L, mp, mc = FL.values()
    rr = np.linspace(0, a, 2001)
    u = ch16.two_layer_velocity(rr, a, d, dp, L, mp, mc)
    assert abs(u[-1]) < 1e-18 and np.all(np.diff(u) <= 1e-18)                         # no slip; monotone
    b = a - d
    Qc_q, _ = quad(lambda s: 2 * PI * s * float(ch16.two_layer_velocity(s, a, d, dp, L, mp, mc)), 0, b, epsabs=0, epsrel=1e-13)
    Qp_q, _ = quad(lambda s: 2 * PI * s * float(ch16.two_layer_velocity(s, a, d, dp, L, mp, mc)), b, a, epsabs=0, epsrel=1e-13)
    q = ch16.two_layer_flow_rates(a, d, dp, L, mp, mc)
    assert rel(q["Qc"], Qc_q) < 1e-11 and rel(q["Qp"], Qp_q) < 1e-11
    assert rel(q["Q"], ch16.two_layer_total_flow_rate(a, d, dp, L, mp, mc)) < 1e-13
    assert rel(ch16.apparent_viscosity_two_layer(a, d, mp, mc), PI * a ** 4 * dp / (8 * q["Q"] * L)) < 1e-12   # (16.32)
    qp = ch16.two_layer_flow_rates(a, d, dp, L, mp, mc, printed_16_30=True)
    assert rel(qp["Qc"], Qc_q) > 0.1                                                    # printed (16.30) fails the integral
    assert rel(qp["Q"], ch16.two_layer_total_flow_rate(a, d, dp, L, mp, mc)) > 0.05


def test_two_layer_V7_limits_and_thin_layer_order():  # V7 — C05, N33, N35
    a, mp, mc = 23e-6, 1.2e-3, 3.9e-3
    assert rel(ch16.apparent_viscosity_two_layer(a, 0.0, mp, mc), mc) < 1e-14          # no layer: the core fluid
    assert rel(ch16.apparent_viscosity_two_layer(a, 3e-6, mp, mp), mp) < 1e-14          # one fluid
    assert ch16.apparent_viscosity_two_layer(a, 3e-6, mp, mc) < mc                      # the FL sign: thinner than the core
    e = np.array([0.0025, 0.005, 0.01, 0.02])                                           # asymptotic range
    err = np.abs(ch16.apparent_viscosity_thin_layer(a, e * a, mp, mc) / ch16.apparent_viscosity_two_layer(a, e * a, mp, mc) - 1)
    assert abs(observed_order(e, err) - 2.0) < 0.15                                     # (16.34) is first order exact
    D = np.array([300, 150, 80, 40, 20]) * 1e-6
    mu = ch16.apparent_viscosity_two_layer(D / 2, 2.5e-6, mp, mc)
    assert np.all(np.diff(mu) < 0)                                                      # smaller tube → smaller μ_app


def test_fahraeus_lindqvist_V1_closure_cell_balance_and_charm_kurland():  # V1 — N30, N36 (16.8), (16.35)–(16.36)
    T = 310.15
    Hc, mp = 0.47, 1.2e-3
    al = 0.070 * math.exp(2.49 * Hc + 1107 / T * math.exp(-1.69 * Hc))
    assert rel(ch16.charm_kurland_core_viscosity(Hc, T, mp), mp / (1 - al * Hc)) < 1e-14
    assert rel(ch16.charm_kurland_core_viscosity(0.0, T, mp), mp) < 1e-15
    with pytest.raises(ValueError):
        ch16.charm_kurland_core_viscosity(0.65, T, mp)
    m = ch16.fahraeus_lindqvist_model(17e-6, 2.1e-6, 0.41, mp, T)
    q = ch16.two_layer_flow_rates(17e-6, 2.1e-6, 1.0, 1.0, mp, m["mu_c"])
    assert abs(m["H_c"] - 0.41 * (1 + q["Qp"] / q["Qc"])) < 1e-12                      # H_c/H_F = 1 + Q_p/Q_c
    assert rel(m["mu_c"], ch16.charm_kurland_core_viscosity(m["H_c"], T, mp)) < 1e-13
    cc = ch16.cell_conservation(q["Qc"], q["Qp"], m["H_c"], 17e-6, 2.1e-6)
    assert rel(cc["H_F"], 0.41) < 1e-12 and rel(cc["Q_F"], q["Qc"] + q["Qp"]) < 1e-15              # (16.8)
    assert rel(cc["HT_over_Hc"], (1 - 2.1 / 17) ** 2) < 1e-13 and rel(cc["H_T"], m["H_T"]) < 1e-12
    assert cc["HT_over_HD"] < 1.0                                                                    # Fåhræus effect


def test_pries_V5_secomb_pries_form_and_minimum():  # V5 — N29 (Secomb & Pries 2013, Eqs. (3)–(5))
    bm = bench()["pries_secomb_pries_2013"]["coefficients"]
    D = np.geomspace(3.0, 500.0, 200)
    H = 0.37
    mu45 = bm["a1"] * np.exp(-bm["b1"] * D) + bm["a2"] - bm["a3"] * np.exp(-bm["b3"] * D ** bm["e3"])
    w = 1 / (1 + bm["k"] * D ** bm["m"])
    C = (bm["c1"] + np.exp(-bm["c2"] * D)) * (-1 + w) + w
    ref = 1 + (mu45 - 1) * ((1 - H) ** C - 1) / ((1 - bm["H_ref"]) ** C - 1)
    assert close(ch16.pries_relative_viscosity(D, H), ref, rtol=1e-13)
    assert close(ch16.pries_relative_viscosity(D, 0.45), mu45, rtol=1e-13)                          # μ_rel(·, 0.45) = μ₀.₄₅
    assert close(ch16.pries_relative_viscosity(D, 1e-12), 1.0, rtol=0, atol=1e-9)                  # no cells: plasma
    res = minimize_scalar(lambda x: float(ch16.pries_relative_viscosity(x, 0.45)), bounds=(3.0, 30.0), method="bounded")
    assert abs(res.x - bench()["pries_secomb_pries_2013"]["statement_min_diameter_um"]) < 1.0       # "down to about 7 µm"


def test_fahraeus_lindqvist_V7_model_trend_against_pries():  # V7 (qualitative, Open item) — trap T19
    Ds = np.array([300, 150, 80, 40, 20.0])
    mod = [ch16.fahraeus_lindqvist_model(D * 1e-6 / 2, 2.5e-6, 0.45, 1.2e-3, 310.15)["mu_rel"] for D in Ds]
    pr = ch16.pries_relative_viscosity(Ds, 0.45)
    assert np.all(np.diff(mod) < 0) and np.all(np.diff(pr) < 0)                        # both fall from 300 to 20 µm


# ======================================================================================================================
# C06 — hoop stress (16.39) and the Hookean tube law (16.41); N37–N45 (distensible tube, compliance)
# ======================================================================================================================
TUBE = dict(a0=3.7e-3, E=0.61e6, h=0.53e-3)


def test_tube_law_V2_derivation():  # V2 — D08: half-tube balance, Hooke, solve for a; (16.55), (16.58)
    p, pe, a, a0, h, E, th, sig = sp.symbols("p p_e a a_0 h E theta sigma", positive=True)
    bal = sp.Eq(2 * sig * h, sp.integrate((p - pe) * a * sp.sin(th), (th, 0, sp.pi)))        # (16.38)
    s39 = sp.solve(bal, sig)[0]
    assert sp.simplify(s39 - (p - pe) * a / h) == 0                                           # (16.39)
    a41 = sp.solve(sp.Eq(s39 / E, (a - a0) / a0), a)[0]                                      # (16.40)
    assert sp.simplify(a41 - a0 / (1 - a0 * (p - pe) / (E * h))) == 0                         # (16.41)
    A, A0, q = sp.symbols("A A_0 q", positive=True)                                            # q = 1 − a₀(p − p_e)/(Eh) > 0
    tl = E * h / a0 * (1 - sp.sqrt(A0 / A))                                                   # (16.55)
    assert sp.simplify(tl.subs(A, sp.pi * (a0 / q) ** 2).subs(A0, sp.pi * a0 ** 2) - E * h / a0 * (1 - q)) == 0
    assert sp.simplify(a41.subs(p, pe + (1 - q) * E * h / a0) - a0 / q) == 0
    lin = sp.series(tl, A, A0, 2).removeO()
    assert sp.simplify(lin - E * h / (2 * a0 * A0) * (A - A0)) == 0                           # (16.58)
    vals = {E: 0.61e6, h: 0.53e-3, a0: 3.7e-3, A0: PI * 3.7e-3 ** 2, A: 1.013 * PI * 3.7e-3 ** 2}
    assert rel(ch16.tube_law_linear(vals[A], vals[A0], vals[E], vals[h], vals[a0]), float(lin.subs(vals))) < 1e-13
    assert rel(ch16.tube_law_hookean(vals[A], vals[A0], vals[E], vals[h], vals[a0]), float(tl.subs(vals))) < 1e-13


def test_tube_law_V1_force_balance_round_trips_and_guard():  # V1 + V7 — C06 (trap T16)
    a0, E, h = TUBE.values()
    dpv = np.array([-2e3, 0.0, 5e3, 13e3, 21e3])
    a = ch16.tube_radius_hookean(dpv + 1.1e3, 1.1e3, a0, E, h)
    sig = ch16.hoop_stress(dpv + 1.1e3, 1.1e3, a, h)
    assert close((a - a0) / a0, sig / E, rtol=1e-12, atol=1e-16)                      # Hooke with the deformed radius
    for dp_, aa in zip(dpv, a):
        f, _ = quad(lambda t: dp_ * aa * math.sin(t), 0, PI)
        assert abs(f - 2 * float(ch16.hoop_stress(dp_, 0.0, aa, h)) * h) <= 1e-12 * max(abs(f), 1)   # (16.38)
    assert close(ch16.tube_law_hookean(PI * a ** 2, PI * a0 ** 2, E, h, a0), dpv, rtol=1e-9, atol=1e-7)
    assert close(ch16.wall_bulge(dpv, a0=a0, E=E, h=h), a - a0, rtol=1e-13, atol=1e-18)
    beta = a0 * 1.0 / (E * h)                                                            # small-amplitude value a₀²Δp/(Eh)(1 + β + …)
    assert abs(ch16.wall_bulge(1.0, a0=a0, E=E, h=h) / (a0 ** 2 / (E * h)) - 1 - beta) < 2 * beta ** 2
    assert rel(ch16.tube_radius_hookean(5e3, 0.0, a0, 1e15, h), a0) < 1e-10             # E → ∞: rigid
    with pytest.raises(ValueError):
        ch16.tube_radius_hookean(E * h / a0 * 1.01, 0.0, a0, E, h)
    with pytest.warns(UserWarning):
        out = ch16.tube_radius_hookean(np.array([1e3, E * h / a0 * 1.5]), 0.0, a0, E, h)
    assert np.isnan(out[1]) and np.isfinite(out[0])
    A = PI * a0 ** 2 * np.array([0.98, 1.0, 1.02])
    assert close(ch16.tube_law_linear(A, PI * a0 ** 2, E, h, a0), E * h / (2 * a0 * PI * a0 ** 2) * (A - PI * a0 ** 2), rtol=1e-13, atol=1e-12)


def test_distensible_flux_V2_derivation_and_slip_S4_units():  # V2 — N39, slip S4
    p, pe, a0, E, h, mu, Q, x = sp.symbols("p p_e a_0 E h mu Q x", positive=True)
    a = a0 / (1 - a0 * (p - pe) / (E * h))
    F = E * h / (3 * a0) * (1 - a0 * (p - pe) / (E * h)) ** -3                        # antiderivative used by (16.45)
    assert sp.simplify(sp.diff(F, p) - (1 - a0 * (p - pe) / (E * h)) ** -4) == 0
    # (16.43) dp/dx = −8μQ/(πa⁴) ⇔ [1 − β]⁻⁴dp = −8μQ dx/(πa₀⁴): the a₀⁴ of the corrected (16.44)
    assert sp.simplify((1 - a0 * (p - pe) / (E * h)) ** -4 * (-8 * mu * Q / (sp.pi * a ** 4)) + 8 * mu * Q / (sp.pi * a0 ** 4)) == 0
    # units: the corrected right side 8μQ dx/(πa₀⁴) is a pressure; the printed 8μQ dx/(πa₀) is not
    corr = 8 * Q_(3.5e-3, "Pa*s") * Q_(1e-6, "m**3/s") * Q_(0.1, "m") / (PI * Q_(3e-3, "m") ** 4)
    printed = 8 * Q_(3.5e-3, "Pa*s") * Q_(1e-6, "m**3/s") * Q_(0.1, "m") / (PI * Q_(3e-3, "m"))
    assert corr.check("[pressure]") and not printed.check("[pressure]")
    vals = {pe: 0.4e3, a0: 3.7e-3, E: 0.61e6, h: 0.53e-3, mu: 3.5e-3}
    Fp = lambda pp: float(F.subs(vals).subs(p, pp))
    Qref = PI * vals[a0] ** 4 / (8 * vals[mu] * 0.17) * (Fp(9.3e3) - Fp(2.1e3))            # ∫ of the corrected (16.44)
    assert rel(ch16.elastic_tube_flow_rate(9.3e3, 2.1e3, 0.4e3, 3.7e-3, 0.61e6, 0.53e-3, 3.5e-3, 0.17), Qref) < 1e-12


def test_distensible_flux_V1_ode_route_and_rigid_limit():  # V1 (independent-route cross-check) + V7 — N39–N41 (16.42)–(16.46)
    a0, E, h = TUBE.values()
    mu, L, pe, p1, p2 = 3.5e-3, 0.17, 0.4e3, 9.3e3, 2.1e3
    Q = ch16.elastic_tube_flow_rate(p1, p2, pe, a0, E, h, mu, L)

    def rhs(xx, y):
        aa = float(ch16.tube_radius_hookean(y[0], pe, a0, E, h))
        return [-8 * mu * Q / (PI * aa ** 4)]                                        # (16.43)

    xs = np.linspace(0, L, 6)
    sol = solve_ivp(rhs, (0, L), [p1], t_eval=xs, rtol=1e-12, atol=1e-9, method="DOP853")
    assert abs(sol.y[0, -1] - p2) < 1e-6 * p1                                         # p(L) = p₂ for that Q
    assert close(ch16.elastic_tube_pressure(xs, Q, p1, pe, a0, E, h, mu), sol.y[0], rtol=1e-9)
    assert rel(ch16.elastic_tube_resistance(p1, p2, pe, a0, E, h, mu, L), (p1 - p2) / Q) < 1e-14
    assert ch16.elastic_tube_resistance(p1, p2, pe, a0, E, h, mu, L) < ch16.vascular_resistance(mu, L, a0)   # distended
    Rs = [ch16.elastic_tube_resistance(pm + 500, pm - 500, 0.0, a0, E, h, mu, L) for pm in (1e3, 4e3, 8e3)]
    assert Rs[0] > Rs[1] > Rs[2]                                                       # Fig. 16.23 mechanism (N44)
    Qr = ch16.elastic_tube_flow_rate(p1, p2, pe, a0, 1e14, h, mu, L)
    assert rel(Qr, PI * a0 ** 4 * (p1 - p2) / (8 * mu * L)) < 1e-8                    # E → ∞: (16.17)


def test_distensible_flux_V3_quadrature_of_16_44_orders_two_and_four():  # V3 — N39–N41: refine the (16.44) integral
    """Q = (πa₀⁴/(8μL))∫_{p₂}^{p₁}(1 − a₀(p − p_e)/(Eh))⁻⁴ dp (the corrected (16.44)) by composite trapezoid and Simpson
    rules on the module's own a(p) = :func:`tube_radius_hookean` — an independent route to the closed-form
    :func:`elastic_tube_flow_rate`; the errors fall at the design orders 2 and 4 (strongly distended tube, β up to 0.34)."""
    a0, E, h = TUBE.values()
    mu, L, pe, p1, p2 = 3.5e-3, 0.17, 0.4e3, 30e3, 2.1e3
    Q = float(ch16.elastic_tube_flow_rate(p1, p2, pe, a0, E, h, mu, L))
    g = lambda p: PI * np.asarray(ch16.tube_radius_hookean(p, pe, a0, E, h), float) ** 4 / (8 * mu * L)
    hs, e_trap, e_simp = [], [], []
    for n in (8, 16, 32, 64):
        p = np.linspace(p2, p1, n + 1)
        y, d = g(p), p[1] - p[0]
        e_trap.append(abs(d * (y.sum() - 0.5 * (y[0] + y[-1])) - Q) / Q)
        e_simp.append(abs(d / 3 * (y[0] + y[-1] + 4 * y[1:-1:2].sum() + 2 * y[2:-1:2].sum()) - Q) / Q)
        hs.append(d)
    assert abs(observed_order(hs, e_trap) - 2.0) < 0.15 and abs(observed_order(hs, e_simp) - 4.0) < 0.15
    assert e_simp[-1] < 2e-8 and np.all(np.diff(e_trap) < 0) and np.all(np.diff(e_simp) < 0)


def test_compliance_V1_central_difference_of_volume():  # V1 — N45
    a0, E, h = TUBE.values()
    L = 0.21
    for p in (0.0, 6e3, 14e3):
        dp = 1e-2
        V = lambda q: PI * float(ch16.tube_radius_hookean(q, 0.0, a0, E, h)) ** 2 * L
        assert rel(ch16.compliance_from_tube(p, 0.0, a0, E, h, L), (V(p + dp) - V(p - dp)) / (2 * dp)) < 1e-7


# ======================================================================================================================
# C07 — the two-element Windkessel (16.49)–(16.52); D09, D10
# ======================================================================================================================
WK = dict(SV=I["SV"], HR=I["HR"], systole_share=I["systole_share"], R=I["R_periph"], K=I["K_compliance"])


def test_windkessel_V1_constant_inflow_and_free_decay():  # V1 — C07 (16.52)
    R, K = 1.7e8, 0.9e-8
    tau = R * K
    t = np.linspace(0, 4.1, 60)
    Q0, p0 = 8.3e-5, 9.1e3
    exact = Q0 * R * (1 - np.exp(-t / tau)) + p0 * np.exp(-t / tau)
    assert close(ch16.windkessel_pressure(t, lambda s: Q0, R, K, p0), exact, rtol=1e-11)
    assert close(ch16.windkessel_ode(t, lambda s: Q0, R, K, p0), exact, rtol=1e-9)
    assert close(ch16.windkessel_pressure(t, lambda s: 0.0, R, K, p0), p0 * np.exp(-t / tau), rtol=1e-13)
    assert rel(ch16.windkessel_rhs(0.3, 7e3, lambda s: Q0, R, K), (Q0 - 7e3 / R) / K) < 1e-15


def test_windkessel_V2_derivation_and_periodic_closed_form():  # V2 — D09, D10
    t, R, K, Qm, ws, p0, s = sp.symbols("t R K Q_m omega_s p_0 s", positive=True)
    p = sp.Function("p")
    tau = R * K
    # D10: the integrating-factor solution (16.52) satisfies (16.49) for any Q (checked for Q = Q_m sin ω_s t)
    Qf = Qm * sp.sin(ws * t)
    sol52 = sp.exp(-t / tau) / K * sp.integrate((Qm * sp.sin(ws * s)) * sp.exp(s / tau), (s, 0, t)) + p0 * sp.exp(-t / tau)
    assert sp.simplify(K * sp.diff(sol52, t) + sol52 / R - Qf) == 0
    assert sp.simplify(sol52.subs(t, 0) - p0) == 0
    # the systolic closed form of design C.12 is the same function
    B = R * Qm * ws * tau / (1 + (ws * tau) ** 2)
    sys_form = R * Qm * (sp.sin(ws * t) - ws * tau * sp.cos(ws * t)) / (1 + (ws * tau) ** 2) + (p0 + B) * sp.exp(-t / tau)
    assert sp.simplify(sp.expand(sys_form - sol52)) == 0
    # numerically: the coded periodic state satisfies (16.49), is continuous at T_s and periodic
    T = 60 / WK["HR"]
    Ts = WK["systole_share"] * T
    tt = np.concatenate([np.linspace(0.01, Ts - 0.01, 20), np.linspace(Ts + 0.01, T - 0.01, 20)])
    hh = 1e-6
    pf = lambda q: np.asarray(ch16.windkessel_periodic_state(q, **WK), float)
    dpdt = (pf(tt + hh) - pf(tt - hh)) / (2 * hh)
    Qv = np.asarray(ch16.ejection_half_sine(tt, WK["SV"], T, Ts), float)
    assert close(WK["K"] * dpdt + pf(tt) / WK["R"], Qv, rtol=0, atol=1e-8 * np.max(Qv))
    assert rel(pf(Ts - 1e-12), pf(Ts + 1e-12)) < 1e-9
    assert rel(pf(T), pf(0.0)) < 1e-12 and rel(pf(0.0), ch16.windkessel_periodic_p0(**WK)) < 1e-14
    assert rel(pf(T - 1e-12), pf(0.0)) < 1e-9                       # periodic: the diastole ends where systole starts


def test_windkessel_V1_periodic_state_reached_by_marching():  # V1 (independent-route cross-check) — C07: ODE and quadrature from rest reach the closed form
    T = 60 / WK["HR"]
    Ts = WK["systole_share"] * T
    Qf = ch16.HalfSineEjection(WK["SV"], T, Ts)
    nb = 25
    t = np.linspace(0, nb * T, nb * 40 + 1)
    p_ode = ch16.windkessel_ode(t, Qf, WK["R"], WK["K"], 0.0)
    p_q = ch16.windkessel_pressure(t, Qf, WK["R"], WK["K"], 0.0)
    assert close(p_ode, p_q, rtol=1e-8, atol=1e-6)
    last = t >= (nb - 1) * T
    assert close(p_q[last], ch16.windkessel_periodic_state(t[last], **WK), rtol=1e-6)
    assert rel(ch16.windkessel_periodic_p0_general(Qf, WK["R"], WK["K"], T), ch16.windkessel_periodic_p0(**WK)) < 1e-10


def test_windkessel_V3_crank_nicolson_march_order_two():  # V3 — C07 (16.49): refine a time march against (16.52)
    """Trapezoidal (Crank–Nicolson) march of K dp/dt + p/R = Q(t) from rest over three beats, with Δt = T/m chosen so
    that T_s falls on a node (the half-sine ejection has a kink there); max error against the closed form
    :func:`windkessel_pressure` falls at order 2."""
    T = 60 / WK["HR"]
    Ts = WK["systole_share"] * T
    Qf = ch16.HalfSineEjection(WK["SV"], T, Ts)
    R, K = WK["R"], WK["K"]
    hs, errs = [], []
    for m in (40, 80, 160, 320):
        assert abs(Ts / (T / m) - round(Ts / (T / m))) < 1e-9                       # T_s on the grid
        dt = T / m
        t = np.arange(3 * m + 1) * dt
        q = np.asarray(Qf(t), float)
        p = np.zeros_like(t)
        A = K / dt
        for j in range(t.size - 1):
            p[j + 1] = ((A - 0.5 / R) * p[j] + 0.5 * (q[j] + q[j + 1])) / (A + 0.5 / R)   # (16.49), trapezoid
        exact = np.asarray(ch16.windkessel_pressure(t, Qf, R, K, 0.0), float)
        hs.append(dt)
        errs.append(np.max(np.abs(p - exact)) / np.max(exact))
    assert abs(observed_order(hs, errs) - 2.0) < 0.15 and errs[-1] < 1e-4
    assert np.all(np.abs(np.asarray(pairwise_orders(hs, errs)) - 2.0) < 0.05)


def test_windkessel_V4_volume_balance_over_a_beat():  # V4 — C07: ∫Q dt = ∫p/R dt; mean p = R·CO; ∫Q = SV
    T = 60 / WK["HR"]
    Ts = WK["systole_share"] * T
    Qint, _ = quad(lambda s: float(ch16.ejection_half_sine(s, WK["SV"], T, Ts)), 0, T, points=[Ts], epsabs=0, epsrel=1e-13)
    assert rel(Qint, WK["SV"]) < 1e-12
    pint, _ = quad(lambda s: float(ch16.windkessel_periodic_state(s, **WK)), 0, T, points=[Ts], epsabs=0, epsrel=1e-13)
    assert rel(pint / WK["R"], Qint) < 1e-10
    sm = ch16.windkessel_summary(**WK)
    assert ch16.windkessel_periodic_summary is ch16.windkessel_summary
    assert rel(sm["p_mean"], WK["R"] * ch16.cardiac_output(WK["SV"], WK["HR"])) < 1e-14
    assert rel(sm["p_mean_sampled"], sm["p_mean"]) < 1e-5
    assert sm["pD"] <= sm["p0"] * (1 + 1e-12) and sm["pS"] > sm["p_mean"] > sm["pD"]
    assert rel(sm["map_third"], sm["pD"] + (sm["pS"] - sm["pD"]) / 3) < 1e-14 and sm["map_error"] != 0.0   # T15
    tt = np.linspace(0.02, T - 0.02, 13)
    tr = ch16.windkessel_terms(tt, **WK)
    assert close(tr["storage"] + tr["runoff"], tr["Q"], rtol=1e-13, atol=1e-18)
    hh = 1e-6
    dp = (np.asarray(ch16.windkessel_periodic_state(tt + hh, **WK)) - np.asarray(ch16.windkessel_periodic_state(tt - hh, **WK))) / (2 * hh)
    assert close(tr["storage"], WK["K"] * dp, rtol=0, atol=1e-8 * WK["SV"] / T)
    trg = ch16.windkessel_terms_general(tt, tr["p"], ch16.HalfSineEjection(WK["SV"], T, Ts), WK["R"], WK["K"])
    assert close(trg["storage"], tr["storage"], rtol=1e-13, atol=1e-18)


def test_windkessel_V5_V7_diastolic_decay_and_limits():  # V5 (Wikipedia form) + V7 — C07
    T = 60 / WK["HR"]
    Ts = WK["systole_share"] * T
    td = np.linspace(Ts, T, 30)
    p = np.asarray(ch16.windkessel_periodic_state(td, **WK))
    assert close(p, p[0] * np.exp(-(td - Ts) / (WK["R"] * WK["K"])), rtol=1e-12)     # P(t_d)e^{−(t−t_d)/RC}
    assert "P(t_d) exp(-(t - t_d)/(R C))" in bench()["windkessel_wikipedia"]["what"]
    # K → 0: p follows QR (no storage); stiffer aorta → larger pulse pressure
    pp = [ch16.windkessel_summary(**{**WK, "K": k})["pulse_pressure"] for k in (2e-8, 1e-8, 0.5e-8)]
    assert pp[0] < pp[1] < pp[2]
    tiny = {**WK, "K": 1e-14}
    tt = np.linspace(0.05, Ts - 0.05, 7)
    assert close(ch16.windkessel_periodic_state(tt, **tiny), WK["R"] * np.asarray(ch16.ejection_half_sine(tt, WK["SV"], T, Ts)), rtol=1e-4)


def test_windkessel_V7_periodic_state_warning_free_for_tiny_RK():  # V7 — C07 guard (loop-1 Should-fix 2), as -W error
    """With R·K ≪ T the unused branch's exponential used to overflow (RuntimeWarning); every warning is an error here."""
    T = 60 / WK["HR"]
    Ts = WK["systole_share"] * T
    t = np.linspace(0.0, 2 * T, 401)                                                  # both phases, two beats
    with warnings.catch_warnings():
        warnings.simplefilter("error")                                                # = python -W error
        for k in (1e-14, 1e-16, 1e-20):
            p = np.asarray(ch16.windkessel_periodic_state(t, **{**WK, "K": k}), float)
            assert np.all(np.isfinite(p)) and np.all(p >= -1e-9 * WK["R"] * WK["SV"] / T)
            ps = ch16.windkessel_periodic_state(0.3 * Ts, **{**WK, "K": k})            # scalar call stays scalar
            assert np.ndim(ps) == 0 and math.isfinite(float(ps))
            sm = ch16.windkessel_summary(**{**WK, "K": k})
            assert all(math.isfinite(float(v)) for v in (sm["pS"], sm["pD"], sm["p_mean"]))
        tau = np.mod(t, T)
        dia = tau > Ts + 1e-3
        p = np.asarray(ch16.windkessel_periodic_state(t, **{**WK, "K": 1e-14}), float)
        assert np.all(np.abs(p[dia]) < 1e-200)                                        # no storage: p = QR = 0 in diastole
    # unchanged on the ordinary branch (shape kept, matches the regression values of the V2 test)
    t2 = np.linspace(0, T, 12).reshape(3, 4)
    assert np.asarray(ch16.windkessel_periodic_state(t2, **WK)).shape == (3, 4)


# ======================================================================================================================
# C08 — the wave equation and the Moens–Korteweg speed (16.56)–(16.68)
# ======================================================================================================================
def test_moens_korteweg_V1_bramwell_hill_and_inverse():  # V1 + V5 (form) — C08, N53
    a0, E, h, rho = 11e-3, 0.45e6, 1.6e-3, 1060.0
    c = ch16.moens_korteweg_speed(E, h, rho, a0)
    A0 = PI * a0 ** 2
    dA = A0 * 1e-6
    dpdA = (ch16.tube_law_hookean(A0 + dA, A0, E, h, a0) - ch16.tube_law_hookean(A0 - dA, A0, E, h, a0)) / (2 * dA)
    assert rel(c, ch16.bramwell_hill_speed(A0, dpdA, rho)) < 1e-9                     # c² = (A/ρ)dp/dA at A₀
    assert rel(c, math.sqrt(E * h / (rho * 2 * a0))) < 1e-15                            # T17: Eh/(2ρa) = Eh/(ρd)
    assert "2 r rho" in bench()["moens_korteweg_wikipedia"]["what"]                     # V5 form: √(E h/(2 r ρ)), r = radius
    cm = ch16.moens_korteweg_speed(E, h, rho, a0, midwall=True)
    assert cm < c and rel(cm, math.sqrt(E * h / (2 * rho * (a0 + h / 2)))) < 1e-15      # (16.63)
    assert rel(ch16.youngs_modulus_from_pwv(c, h, rho, a0), E) < 1e-14
    assert rel(ch16.youngs_modulus_from_pwv(cm, h, rho, a0, midwall=True), E) < 1e-14
    hs = np.array([1e-3, 1e-4, 1e-5]) * a0
    gap = [1 - ch16.moens_korteweg_speed(E, x, rho, a0, midwall=True) / ch16.moens_korteweg_speed(E, x, rho, a0) for x in hs]
    assert gap[0] > gap[1] > gap[2] > 0                                                  # h/a → 0: forms agree


def test_moens_korteweg_V2_derivation():  # V2 — D11, D12: linearise (16.53)–(16.55), cross-differentiate
    x, t, rho, E, h, a0, eps = sp.symbols("x t rho E h a_0 epsilon", positive=True)
    A0 = sp.pi * a0 ** 2
    Ah = sp.Function("Ah")(x, t)
    uh = sp.Function("uh")(x, t)
    A = A0 + eps * Ah
    u = eps * uh
    p = E * h / a0 * (1 - sp.sqrt(A0 / A))                                              # (16.55), p_e = 0
    cont = sp.diff(A, t) + sp.diff(A * u, x)                                            # (16.53)
    mom = rho * A * (sp.diff(u, t) + u * sp.diff(u, x)) + A * sp.diff(p, x)              # (16.54) corrected
    c1 = sp.simplify(sp.diff(cont, eps).subs(eps, 0))
    m1 = sp.simplify(sp.diff(mom, eps).subs(eps, 0))
    assert sp.simplify(c1 - (sp.diff(Ah, t) + A0 * sp.diff(uh, x))) == 0                # (16.56)
    dpdA = sp.simplify(sp.diff(p, eps).subs(eps, 0) / Ah)
    assert sp.simplify(dpdA - E * h / (2 * a0 * A0)) == 0                               # (16.58)
    assert sp.simplify(m1 / A0 - (rho * sp.diff(uh, t) + dpdA * sp.diff(Ah, x))) == 0     # (16.57)
    # the printed (16.54) −∂((p − p_e)A)/∂x linearises to the same equation (S5: different only when nonlinear)
    mom_pr = rho * A * (sp.diff(u, t) + u * sp.diff(u, x)) + sp.diff(p * A, x)
    assert sp.simplify(sp.diff(mom_pr, eps).subs(eps, 0) - m1) == 0
    assert sp.simplify(sp.diff(mom_pr, eps, 2).subs(eps, 0) - sp.diff(mom, eps, 2).subs(eps, 0)) != 0
    # cross-differentiation: f(x − ct) solves both iff c² = (A₀/ρ)dp/dA = Eh/(2ρa₀)  (16.59)–(16.62)
    c = sp.symbols("c", positive=True)
    f = sp.Function("f")
    Aw = f(x - c * t)
    uw = c * Aw / A0                                                                    # from (16.56)
    res = rho * sp.diff(uw, t) + dpdA * sp.diff(Aw, x)
    csol = sp.solve(sp.simplify(res / sp.diff(f(x - c * t), x)), c)
    assert any(sp.simplify(cs ** 2 - E * h / (2 * rho * a0)) == 0 for cs in csol)


def test_pulse_wave_V1_dalembert_and_water_hammer():  # V1 — R09, R10 (16.65)–(16.68)
    c, rho = 5.3, 1060.0
    f1 = lambda s: 900.0 * np.exp(-((np.asarray(s) - 0.2) / 0.05) ** 2)
    f2 = lambda s: 300.0 * np.exp(-((np.asarray(s) - 0.3) / 0.04) ** 2)
    x = np.linspace(-1, 2, 61)
    t0 = 0.31
    p = ch16.pulse_wave(x, t0, c, f1, f2)
    assert close(p, f1(t0 - x / c) + f2(t0 + x / c), rtol=1e-13, atol=1e-10)
    hx, ht = 1e-4, 1e-5
    pxx = (ch16.pulse_wave(x + hx, t0, c, f1, f2) - 2 * p + ch16.pulse_wave(x - hx, t0, c, f1, f2)) / hx ** 2
    ptt = (ch16.pulse_wave(x, t0 + ht, c, f1, f2) - 2 * p + ch16.pulse_wave(x, t0 - ht, c, f1, f2)) / ht ** 2
    assert np.max(np.abs(ptt - c ** 2 * pxx)) < 1e-4 * np.max(np.abs(ptt))                # (16.61)
    u = ch16.velocity_from_pulse_shapes(x, t0, c, rho, f1, f2)
    assert close(u, (f1(t0 - x / c) - f2(t0 + x / c)) / (rho * c), rtol=1e-12, atol=1e-14)
    ut = (ch16.velocity_from_pulse_shapes(x, t0 + ht, c, rho, f1, f2) - ch16.velocity_from_pulse_shapes(x, t0 - ht, c, rho, f1, f2)) / (2 * ht)
    px = (ch16.pulse_wave(x + hx, t0, c, f1, f2) - ch16.pulse_wave(x - hx, t0, c, f1, f2)) / (2 * hx)
    assert np.max(np.abs(rho * ut + px)) < 1e-5 * np.max(np.abs(px))                     # (16.57)
    assert rel(ch16.velocity_from_pulse(1300.0, rho, c), 1300 / (rho * c)) < 1e-15
    assert rel(ch16.velocity_from_pulse(1300.0, rho, c, direction=-1), -1300 / (rho * c)) < 1e-15
    lr = ch16.linearisation_ratios(2.6e3, 13e3, 0.9, 5.3, Eh_over_a0=65e3)
    assert rel(lr["p_ratio"], 0.2) < 1e-14 and rel(lr["u_over_c"], 0.9 / 5.3) < 1e-14 and rel(lr["strain"], 0.04) < 1e-14


def test_pulse_state_V1_right_going_and_released_from_rest():  # V1 — E3 parity function
    a0, E, h, rho = 4e-3, 0.7e6, 0.6e-3, 1060.0
    c = ch16.moens_korteweg_speed(E, h, rho, a0)
    x = np.linspace(-0.5, 1.5, 81)
    st = ch16.pulse_state(x, 0.07, c=c, amplitude=1.9e3, width=0.06, rho=rho, a0=a0, E=E, h=h)
    g = 1.9e3 * np.exp(-((x - c * 0.07) / 0.06) ** 2)
    assert close(st["p"], g, rtol=1e-12, atol=1e-9) and close(st["u"], g / (rho * c), rtol=1e-12, atol=1e-12)
    assert close(st["a"], ch16.tube_radius_hookean(g, 0.0, a0, E, h), rtol=1e-14)
    sr = ch16.pulse_state(x, 0.07, c=c, amplitude=1.9e3, width=0.06, rho=rho, a0=a0, E=E, h=h, ic="rest")
    gl = 1.9e3 * np.exp(-((x + c * 0.07) / 0.06) ** 2)
    assert close(sr["p"], 0.5 * (g + gl), rtol=1e-12, atol=1e-9) and close(sr["u"], 0.5 * (g - gl) / (rho * c), rtol=1e-12, atol=1e-12)


def _tube1d_case(n, eps, printed=False, nu=0.02):
    a0, h, E, rho = 0.01, 1e-3, 0.5e6, 1060.0
    A0 = PI * a0 ** 2
    c0 = math.sqrt(E * h / (2 * rho * a0))
    x = np.linspace(0, 2.0, n, endpoint=False)
    g = np.exp(-((x - 0.6) / 0.08) ** 2)
    A = A0 * (1 + eps * g)
    u = 4 * (c0 - c0 * (A0 / A) ** 0.25)                      # simple right-going wave of the corrected equations
    out = ch16.tube1d_run(x, A, u, A0, E, h, a0, rho, 0.5 / c0, printed_16_54=printed, nu_art=nu)
    return x, out, A0, c0


def test_tube1d_V3_small_pulse_speed_order_and_mass():  # V3 + V4 — C08 (16.53)–(16.55), N50
    hs, es = [], []
    for n in (200, 400, 800):
        x, o, A0, c0 = _tube1d_case(n, 1e-6)
        ex = 1e-6 * np.exp(-((((x - 0.6 - c0 * o["t"]) + 1.0) % 2.0 - 1.0) / 0.08) ** 2)
        es.append(np.max(np.abs(o["A"] / A0 - 1 - ex)) / 1e-6)
        hs.append(x[1] - x[0])
        assert abs(o["mass"] / o["mass0"] - 1) < 1e-13                                   # Σ A Δx conserved
    assert abs(observed_order(hs, es) - 2.0) < 0.15
    x, o, A0, c0 = _tube1d_case(800, 1e-3)
    d = o["A"] - A0
    shift = np.sum(x * d) / np.sum(d) - 0.6
    assert abs(shift / o["t"] / c0 - 1) < 5e-3                                           # travels at c(A₀) (16.62)


def test_tube1d_V7_printed_16_54_breaks_the_simple_wave():  # V7 (structure: a simple wave leaves rest behind) — slip S5: printed_16_54 must fail
    for printed, bound in ((False, 1e-5), (True, None)):
        x, o, A0, c0 = _tube1d_case(800, 0.2, printed=printed)
        left = (x > 0.1) & (x < 0.45)                       # behind the pulse: a simple wave leaves it at rest
        back = np.max(np.abs(o["A"][left] / A0 - 1))
        if printed:
            assert back > 1e-3                              # the printed momentum flux sheds a reflected wave
        else:
            assert back < bound


# ======================================================================================================================
# C09 — Womersley's profile (16.72)–(16.88); D13–D15; N55–N59
# ======================================================================================================================
def _womersley_U_mpmath(r, a, omega, rho, mu, G):
    """(16.80) U(r) = −(iG/ωρ)(1 − J₀(i^{3/2}α r/a)/J₀(i^{3/2}α)) in 30-digit arithmetic (independent of fluidpy)."""
    import mpmath as mp
    mp.mp.dps = 30
    alpha = a * mp.sqrt(omega * rho / mu)
    z = mp.exp(3j * mp.pi / 4) * alpha
    return complex(-1j * G / (omega * rho) * (1 - mp.besselj(0, z * r / a) / mp.besselj(0, z)))


def test_womersley_profile_V1_against_30_digit_bessel():  # V1 — C09 (16.80)–(16.82)
    a, rho, G = 2.9e-3, 1060.0, 713.0
    for alpha in (0.7, 3.0, 8.0, 20.0, 45.0):
        mu = 3.5e-3
        omega = alpha ** 2 * mu / (rho * a ** 2)
        r = np.linspace(0, a, 9)
        for t in (0.0, 0.37 * 2 * PI / omega):
            ref = np.array([(_womersley_U_mpmath(x, a, omega, rho, mu, G) * np.exp(1j * omega * t)).real for x in r])
            got = ch16.womersley_velocity(r, t, a, omega, rho, mu, G)
            assert np.max(np.abs(got - ref)) < 1e-11 * np.max(np.abs(ref)) + 1e-300
        nd = ch16.womersley_velocity_nd(r / a, 1.1, alpha)
        assert close(nd, ch16.womersley_velocity(r, 1.1 / omega, a, omega, rho, mu, G) / (G * a ** 2 / (4 * mu)), rtol=1e-11, atol=1e-13)
    assert close(ch16.womersley_F1(np.array([2.9e-3]), 2.9e-3, 8.0), 1.0, rtol=1e-14)            # F₁(a) = 1
    for al_ in (1e-4, 0.5, 6.0):
        ss = np.array([0.0, 0.4, 0.9])
        ref = np.array([_womersley_U_mpmath(x, 1.0, al_ ** 2, 1.0, 1.0, 1.0) for x in ss])   # U/(Ga²/μ) with a = μ = ρ = G = 1
        assert np.max(np.abs(ch16.womersley_profile_complex(ss, 1.0, al_) - ref)) < 1e-12 * np.max(np.abs(ref))
    assert close(ch16.womersley_F1(np.linspace(0, 1, 5), 1.0, 0.0), 1.0, rtol=0)                 # α = 0


def test_womersley_profile_V2_derivation():  # V2 — D13 (★★★): every line of (16.74)–(16.82), slip S6
    r, t, a, om, rho, mu, G = sp.symbols("r t a omega rho mu G", positive=True)
    U = sp.Function("U")
    # (16.72) with −∂p/∂x = G e^{iωt} (16.73) and u = U(r)e^{iωt} (16.75): ρ∂u/∂t = G e^{iωt} + μ(u_rr + u_r/r)
    u = U(r) * sp.exp(sp.I * om * t)
    pde = rho * sp.diff(u, t) - G * sp.exp(sp.I * om * t) - mu * (sp.diff(u, r, 2) + sp.diff(u, r) / r)
    ode = sp.expand(sp.simplify(pde / (-mu * sp.exp(sp.I * om * t))))
    lhs76 = sp.diff(U(r), r, 2) + sp.diff(U(r), r) / r - sp.I * om * rho / mu * U(r)
    assert sp.simplify(ode - (lhs76 + G / mu)) == 0                    # LHS of (16.76) = −G/μ (printed +G/μ: S6)
    # (16.77): the particular solution G/(ωρi) satisfies the corrected equation and not the printed one
    Up = G / (om * rho * sp.I)
    assert sp.simplify(lhs76.subs(U(r), Up).doit() + G / mu) == 0
    assert sp.simplify(lhs76.subs(U(r), Up).doit() - G / mu) != 0
    # J₀(λr) with λ = i√(iωρ/μ) = i^{3/2}√(ωρ/μ) solves the homogeneous equation (λ² = −iωρ/μ)
    lam = sp.I * sp.sqrt(sp.I * om * rho / mu)
    assert sp.simplify(lam ** 2 + sp.I * om * rho / mu) == 0
    hom = sp.diff(sp.besselj(0, lam * r), r, 2) + sp.diff(sp.besselj(0, lam * r), r) / r - sp.I * om * rho / mu * sp.besselj(0, lam * r)
    assert abs(complex(hom.subs({r: 0.37, om: 2.3, rho: 1.1, mu: 0.07}).evalf(30))) < 1e-12        # double-precision subs
    # (16.78)–(16.80): C₁ from U(a) = 0 and the profile; (16.82): −iG/(ωρ) = Ga²/(iμα²) with α² = a²ωρ/μ
    al = a * sp.sqrt(om * rho / mu)
    C1 = sp.I * G / (om * rho) / sp.besselj(0, lam * a)                                 # (16.79)
    U80 = C1 * sp.besselj(0, lam * r) + Up
    assert sp.simplify(U80.subs(r, a)) == 0
    assert sp.simplify(U80 - (-sp.I * G / (om * rho)) * (1 - sp.besselj(0, lam * r) / sp.besselj(0, lam * a))) == 0
    assert sp.simplify(-sp.I * G / (om * rho) - G * a ** 2 / (sp.I * mu * al ** 2)) == 0
    assert abs(complex((sp.I * sp.sqrt(sp.I)).evalf(30)) - I32) < 1e-15                     # the argument i^{3/2}α
    assert sp.simplify(lam * a - sp.I * sp.sqrt(sp.I) * al) == 0
    # the coded profile satisfies the corrected (16.76): residual = −G/μ (finite differences)
    a_, om_, rho_, mu_, G_ = 3.3e-3, 7.1, 1060.0, 3.5e-3, 517.0
    res = ch16.womersley_ode_residual(np.linspace(0.1, 0.9, 9) * a_, a_, om_, rho_, mu_, G_)
    assert np.max(np.abs(res + G_ / mu_)) < 1e-5 * G_ / mu_


def test_womersley_limits_V2_derivation():  # V2 — D14 (small α, (16.87)) and D15 (large α, (16.88), slip S7)
    s, x, al, y, G, om, rho, mu, n = sp.symbols("s x alpha y G omega rho mu n", positive=True)
    z = sp.exp(3 * sp.I * sp.pi / 4) * al
    J0 = lambda w: sum((-1) ** m * (w / 2) ** (2 * m) / sp.factorial(m) ** 2 for m in range(4))   # P374, 4 terms
    ratio = sp.series(J0(z * s) / J0(z), al, 0, 5).removeO()                                    # D14 steps 1–4
    one_minus_F1 = sp.simplify(sp.expand(1 - ratio))
    lead = sp.simplify(sp.expand(one_minus_F1).coeff(al, 2))
    assert sp.simplify(lead - sp.I * (1 - s ** 2) / 4) == 0                                     # 1 − F₁ ≈ iα²(1 − s²)/4
    shape = sp.simplify(sp.expand(one_minus_F1 / (sp.I * al ** 2)).subs(al, 0))
    assert sp.simplify(shape - (1 - s ** 2) / 4) == 0                                           # (16.87): Poiseuille
    # D15: near the wall (y = a − r ≪ a) the 1/r term drops: U″ − (iωρ/μ)U = −G/μ, U(0) = 0, U bounded as y grows
    lam = sp.sqrt(n * om * rho / (2 * mu)) * (1 + sp.I)                                         # harmonic n (S7)
    assert sp.simplify(lam ** 2 - sp.I * n * om * rho / mu) == 0                                # the root with Re > 0 decays
    U = G / (sp.I * n * om * rho) * (1 - sp.exp(-lam * y))
    assert sp.simplify(sp.diff(U, y, 2) - sp.I * n * om * rho / mu * U + G / mu) == 0 and U.subs(y, 0) == 0
    assert sp.simplify(G / (sp.I * n * om * rho) - G * x ** 2 / (sp.I * mu * (x * sp.sqrt(n * om * rho / mu)) ** 2)) == 0   # Ga²/(iμα_n²)


def test_elastic_tube_V2_derivation_long_wave_modes():  # V2 — D18: (16.118)–(16.125) from normal modes (16.116)
    x, r, t, k, om, rho, mu = sp.symbols("x r t k omega rho mu")
    uh, vh, ph = sp.Function("uh")(r), sp.Function("vh")(r), sp.symbols("ph")
    E = sp.exp(sp.I * (k * x - om * t))                                                       # T1: e^{i(kx − ωt)}
    u, v, p = uh * E, vh * E, ph * E
    cont = sp.diff(u, x) + sp.diff(r * v, r) / r                                                # (16.118)
    mom = rho * sp.diff(u, t) + sp.diff(p, x) - mu * (sp.diff(u, r, 2) + sp.diff(u, r) / r)      # (16.119)
    assert sp.simplify(cont * r / E - (sp.diff(r * vh, r) + sp.I * k * r * uh)) == 0           # (16.122)
    m = sp.expand(sp.simplify(mom / (-mu * E)))
    assert sp.simplify(m - (sp.diff(uh, r, 2) + sp.diff(uh, r) / r + sp.I * om * rho / mu * uh - sp.I * k * ph / mu)) == 0   # (16.123)
    xi, eta = sp.symbols("xi eta")
    assert sp.simplify(sp.diff(xi * E, t) / E + sp.I * om * xi) == 0                            # (16.124): u = ∂ξ/∂t at the wall
    assert sp.simplify(sp.diff(eta * E, t) / E + sp.I * om * eta) == 0                          # (16.125): v = ∂η/∂t


def test_womersley_profile_V3_crank_nicolson_order_two():  # V3 — C09 (independent PDE route, (16.72))
    for alpha in (3.0, 8.0):
        a, om, rho = 1.0, 1.0, 1.0
        mu = om * a * a / alpha ** 2
        hs, es = [], []
        for nr, ns in ((20, 40), (40, 80), (80, 160), (160, 320)):
            s = ch16.womersley_pde_solve(a, om, rho, mu, 1.0, nr=nr, steps_per_period=ns)
            R, T = np.meshgrid(s["r"], s["t"])
            ex = ch16.womersley_velocity(R, T, a, om, rho, mu, 1.0)
            es.append(np.max(np.abs(s["u"] - ex)) / np.max(np.abs(ex)))
            hs.append(s["dr"])
        assert abs(observed_order(hs, es) - 2.0) < 0.15 and es[-1] < 1e-4
    s = ch16.womersley_pde_solve(1.0, 1.0, 1.0, 1 / 9, 1.0, nr=40, steps_per_period=80, periodic=False, periods=60)
    s2 = ch16.womersley_pde_solve(1.0, 1.0, 1.0, 1 / 9, 1.0, nr=40, steps_per_period=80)
    assert np.max(np.abs(s["u"] - s2["u"])) < 1e-6 * np.max(np.abs(s2["u"]))            # marching reaches the periodic map


def test_womersley_profile_V7_small_and_large_alpha():  # V7 — D14, D15, (16.87)–(16.88), N58 (slip S7)
    a, rho, mu, G = 1.0, 1.0, 1.0, 1.0
    r = np.linspace(0, a, 41)
    errs = []
    als = np.array([0.05, 0.1, 0.2, 0.4])
    for al in als:
        om = al ** 2 * mu / (rho * a ** 2)
        t = np.linspace(0, 2 * PI / om, 13)
        R, T = np.meshgrid(r, t)
        ex = ch16.womersley_velocity(R, T, a, om, rho, mu, G)
        errs.append(np.max(np.abs(ch16.womersley_small_alpha(R, T, a, om, rho, mu, G) - ex)) / (G * a ** 2 / (4 * mu)))
    assert abs(observed_order(als, errs) - 2.0) < 0.15                                  # quasi-steady error ∝ α²
    errl = []
    alsl = np.array([20.0, 40.0, 80.0, 160.0])
    for al in alsl:
        om = al ** 2 * mu / (rho * a ** 2)
        t = np.linspace(0, 2 * PI / om, 25)
        rr = np.concatenate([np.linspace(0, a * (1 - 12 / al), 21), a - a / al * np.linspace(0, 12, 121)[::-1]])   # resolve the layer
        R, T = np.meshgrid(rr, t)
        ex = ch16.womersley_velocity(R, T, a, om, rho, mu, G)
        errl.append(np.max(np.abs(ch16.womersley_large_alpha(R, T, a, om, rho, mu, G) - ex)) / (G / (rho * om)))
        core = ch16.womersley_velocity(0.0, t, a, om, rho, mu, G)                    # plug: ρ∂u/∂t = G cos ωt
        assert np.max(np.abs(core - G * np.sin(om * t) / (rho * om))) < 3.0 / al * G / (rho * om)
    assert abs(observed_order(1 / alsl, errl) - 1.0) < 0.15                              # thin-layer error ∝ 1/α
    # near the wall the large-α profile is ch08's Stokes second problem in the oscillating frame
    al = 60.0
    om = al ** 2 * mu / (rho * a ** 2)
    y = np.linspace(0, 6 * math.sqrt(2 * mu / (rho * om)), 25)
    t = 0.3
    plug = G * np.sin(om * t) / (rho * om)
    # seen from the plug, the wall moves with −plug = (G/ρω)cos(ωt + π/2): ch08's Stokes second problem, shifted by T/4
    st = np.asarray(LAM.stokes_second_problem(y, t + PI / (2 * om), U=G / (rho * om), omega=om, nu=mu / rho))
    rel_frame = ch16.womersley_large_alpha(a - y, t, a, om, rho, mu, G) - plug
    assert np.max(np.abs(rel_frame - st)) < 1e-12 * G / (rho * om)


def test_womersley_bessel_series_V1():  # V1 — P374 from-scratch J₀, J₁ (E4/E5 parity, design C.12)
    for al in (0.5, 3.0, 8.0, 12.0, 20.0, 25.0, 30.0):
        z = I32 * al
        assert rel(ch16.bessel_j0_series(z, 100), jv(0, z)) < 1e-12
        assert rel(ch16.bessel_j1_series(z, 100), jv(1, z)) < 1e-12
    zz = np.array([0.3 + 0.1j, 2.0 - 1.0j, 5.5j])
    assert close(ch16.bessel_ratio_j0_j1(zz), jv(0, zz) / jv(1, zz), rtol=1e-13)


def test_womersley_harmonics_V1_superposition():  # V1 — N57 (16.83)–(16.86)
    a, om, rho, mu = 2.1e-3, 7.3, 1060.0, 3.5e-3
    Gn = [310.0, 520.0 - 140.0j, 90.0 + 60.0j]
    r = np.linspace(0, a, 11)
    t = np.linspace(0, 2 * PI / om, 7)
    R, T = np.meshgrid(r, t)
    got = ch16.womersley_harmonics(R, T, a, om, rho, mu, Gn)
    mean = ch16.poiseuille_velocity(R, a, Gn[0], 1.0, mu)
    h1 = ch16.womersley_velocity(R, T, a, om, rho, mu, Gn[1])
    h2 = ch16.womersley_velocity(R, T, a, 2 * om, rho, mu, Gn[2])                          # α₂ = α√2 via 2ω
    assert close(got, mean + h1 + h2, rtol=1e-12, atol=1e-15)


def test_womersley_flow_rate_V7_scalar_gradient_has_one_meaning():  # V7 — C10, N57 (review Should-fix 1)
    """A bare scalar G means the n = 1 amplitude in womersley_flow_rate; the harmonics form refuses it (it could be
    read as the mean G₀ — 245× apart at these inputs); [G₀] is a steady mean, [0, G] the single harmonic."""
    a, om, rho, mu = 2.1e-3, 7.3, 1060.0, 3.5e-3
    t = np.linspace(0, 2 * PI / om, 9)
    for scalar in (100.0, np.float64(100.0), np.array(100.0), 100, 100.0 + 20.0j):
        with pytest.raises(ValueError, match="scalar"):
            ch16.womersley_flow_rate_harmonics(t, a, om, rho, mu, scalar)
    one = ch16.womersley_flow_rate(t, a, om, rho, mu, 100.0)
    assert close(one, ch16.womersley_flow_rate_harmonics(t, a, om, rho, mu, [0.0, 100.0]), rtol=1e-14, atol=1e-20)
    assert close(ch16.womersley_flow_rate_harmonics(t, a, om, rho, mu, [100.0]), PI * 100.0 * a ** 4 / (8 * mu),
                 rtol=1e-14)                                                                       # steady mean, (16.17)
    assert close(ch16.womersley_flow_rate(t, a, om, rho, mu, [100.0]), PI * 100.0 * a ** 4 / (8 * mu), rtol=1e-14)
    Qmean = PI * 100.0 * a ** 4 / (8 * mu)
    assert abs(np.max(np.abs(one)) - Qmean) > 0.3 * Qmean and abs(np.mean(one[:-1])) < 1e-12 * Qmean   # zero mean


def test_womersley_V5_wikipedia_forms():  # V5 — C09, C10 (Wikipedia "Pulsatile flow", Womersley 1955; P′ = −G)
    a, om, rho, mu = 2.7e-3, 7.9, 1060.0, 3.5e-3
    G = 450.0
    alpha = a * math.sqrt(om * rho / mu)
    Lam = alpha * I32
    Pp = -G
    r = np.linspace(0, a, 9)
    t = np.linspace(0, 2 * PI / om, 11)
    R, T = np.meshgrid(r, t)
    u_w = np.real(1j * Pp / (rho * om) * (1 - jv(0, Lam * R / a) / jv(0, Lam)) * np.exp(1j * om * T))
    assert close(ch16.womersley_velocity(R, T, a, om, rho, mu, G), u_w, rtol=1e-11, atol=1e-14)
    Q_w = np.real(PI * a ** 2 * 1j * Pp / (rho * om) * (1 - 2 / Lam * jv(1, Lam) / jv(0, Lam)) * np.exp(1j * om * t))
    assert close(ch16.womersley_flow_rate(t, a, om, rho, mu, G), Q_w, rtol=1e-11, atol=1e-18)
    tw_w = np.real(Pp * a / Lam * jv(1, Lam) / jv(0, Lam) * np.exp(1j * om * t))           # their τ_w = μ∂u/∂r|_R
    assert close(ch16.womersley_wall_shear(t, a, om, rho, mu, G), -tw_w, rtol=1e-11, atol=1e-12)   # ours is −μ∂u/∂r (S8)


# ======================================================================================================================
# C10 — flow rate, amplitude and phase lag (16.89)–(16.100); D16, D17; N60–N64
# ======================================================================================================================
def test_womersley_flow_V1_quadrature_wall_shear_amplitude_and_phase():  # V1 — C10, N62 (slip S8)
    a, rho, mu, G = 3.1e-3, 1060.0, 3.5e-3, 640.0
    for alpha in (1.5, 4.0, 8.0, 15.0):
        om = alpha ** 2 * mu / (rho * a ** 2)
        for t in (0.0, 0.41 / om, 2.2 / om):
            Qq, _ = quad(lambda s: 2 * PI * s * float(ch16.womersley_velocity(s, t, a, om, rho, mu, G)), 0, a,
                         epsabs=0, epsrel=1e-12, limit=200)
            assert rel(ch16.womersley_flow_rate(t, a, om, rho, mu, G), Qq) < 1e-10
            h = a * 1e-5
            dudr = (float(ch16.womersley_velocity(a, t, a, om, rho, mu, G)) - float(ch16.womersley_velocity(a - h, t, a, om, rho, mu, G))) / h
            d2 = (float(ch16.womersley_velocity(a, t, a, om, rho, mu, G)) - 2 * float(ch16.womersley_velocity(a - h, t, a, om, rho, mu, G))
                  + float(ch16.womersley_velocity(a - 2 * h, t, a, om, rho, mu, G))) / h ** 2
            tau_fd = -mu * (dudr + 0.5 * h * d2)                                         # second-order one-sided
            assert abs(ch16.womersley_wall_shear(t, a, om, rho, mu, G) - tau_fd) < 1e-6 * G * a
        tt = np.linspace(0, 2 * PI / om, 20001)
        Qt = ch16.womersley_flow_rate(tt, a, om, rho, mu, G)
        Q0 = PI * G * a ** 4 / (8 * mu)
        assert rel(np.max(Qt) / Q0, ch16.womersley_flow_ratio(alpha)) < 1e-7              # (16.99)/(16.95)
        X, Y = ch16.womersley_XY(alpha)
        assert rel(complex(X, Y), 1 - complex(ch16.womersley_F2(alpha))) < 1e-14           # (16.97)
        assert rel(ch16.womersley_flow_ratio(alpha), 8 * math.hypot(X, Y) / alpha ** 2) < 1e-12
        phi = ch16.womersley_phase_lag(alpha)
        assert rel(math.tan(phi), X / Y) < 1e-12                                            # (16.100)
        t_pk = tt[np.argmax(Qt)]
        assert abs(((om * t_pk - phi) + PI) % (2 * PI) - PI) < 2 * PI / 20000 * 1.01         # Q peaks φ/ω after G cos ωt
    G_n = [210.0, 640.0, -130.0 + 70.0j]
    om = 7.0
    tt = np.linspace(0, 1.0, 9)
    Qh = ch16.womersley_flow_rate_harmonics(tt, a, om, rho, mu, G_n)
    Qs = (PI * G_n[0] * a ** 4 / (8 * mu) + ch16.womersley_flow_rate(tt, a, om, rho, mu, G_n[1])
          + ch16.womersley_flow_rate(tt, a, 2 * om, rho, mu, G_n[2]))
    assert close(Qh, Qs, rtol=1e-12)
    taw = ch16.womersley_wall_shear(tt, a, om, rho, mu, G_n)
    assert close(taw, G_n[0] * a / 2 + ch16.womersley_wall_shear(tt, a, om, rho, mu, G_n[1])
                 + ch16.womersley_wall_shear(tt, a, 2 * om, rho, mu, G_n[2]), rtol=1e-12)


def test_womersley_flow_V2_derivation():  # V2 — D16, D17: ∫sJ₀ = zJ₁, (16.90)–(16.91), (16.98)–(16.100)
    z, s, X, Y, w, t = sp.symbols("z s X Y omega t", real=True)
    assert sp.simplify(sp.diff(z * sp.besselj(1, z), z) - z * sp.besselj(0, z)) == 0          # P376
    import mpmath as mp
    mp.mp.dps = 25
    for al in (2.0, 6.0, 11.0):                     # ∫₀¹ 2(1 − F₁(s))s ds = 1 − F₂  (D16, a = 1)
        zz = mp.exp(3j * mp.pi / 4) * al
        lhs = mp.quad(lambda x: 2 * x * (1 - mp.besselj(0, zz * x) / mp.besselj(0, zz)), [0, 1])
        assert abs(complex(lhs) - (1 - complex(ch16.womersley_F2(al)))) < 1e-14
        assert abs(complex(2 * mp.besselj(1, zz) / (zz * mp.besselj(0, zz))) - complex(ch16.womersley_F2(al))) < 1e-14
    # D17: Q ∝ (Z/i)e^{iωt}; its real and imaginary parts are (16.98); |Q| and tan φ (16.99)–(16.100)
    Qc = sp.expand((X + sp.I * Y) / sp.I * (sp.cos(w * t) + sp.I * sp.sin(w * t)))
    assert sp.simplify(sp.re(Qc) - (Y * sp.cos(w * t) + X * sp.sin(w * t))) == 0
    assert sp.simplify(sp.im(Qc) - (-(X * sp.cos(w * t) - Y * sp.sin(w * t)))) == 0
    phi = sp.atan2(X, Y)
    Xp, Yp = sp.symbols("X_p Y_p", positive=True)
    amp_phase = sp.sqrt(Xp ** 2 + Yp ** 2) * sp.cos(w * t - sp.atan(Xp / Yp))
    assert sp.simplify(sp.expand_trig(amp_phase) - (Yp * sp.cos(w * t) + Xp * sp.sin(w * t))) == 0


def test_womersley_flow_V7_limits_series_and_slip_16_96():  # V7 — N63, N64 (16.93), (16.96), slip S10
    als = np.array([0.02, 0.04, 0.08, 0.16])
    e2 = [abs(ch16.womersley_F2_series(x, 2) - ch16.womersley_F2(x)) for x in als]
    e4 = [abs(ch16.womersley_F2_series(x, 4) - ch16.womersley_F2(x)) for x in als]
    assert abs(observed_order(als, e2) - 4.0) < 0.15 and abs(observed_order(als, e4) - 6.0) < 0.15
    assert rel(ch16.womersley_F2_series(0.3, 2), 1 - 0.3 ** 2 * 1j / 8) < 1e-15                # (16.93)
    big = np.array([40.0, 80.0, 160.0, 320.0])
    ec = [abs(ch16.womersley_F2_asymptote(x) - ch16.womersley_F2(x)) for x in big]
    ep = [abs(ch16.womersley_F2_asymptote(x, printed_16_96=True) - ch16.womersley_F2(x)) for x in big]
    el = [abs(2 / (I12 * x) - ch16.womersley_F2(x)) for x in big]                              # leading term alone
    assert abs(observed_order(1 / big, ec) - 3.0) < 0.15
    assert abs(observed_order(1 / big, ep) - 2.0) < 0.15 and abs(observed_order(1 / big, el) - 2.0) < 0.15
    assert all(p_ > l_ for p_, l_ in zip(ep, el))                    # the printed correction is worse than none
    a = np.geomspace(0.01, 200, 300)
    ratio, phi = ch16.womersley_flow_ratio(a), ch16.womersley_phase_lag(a)
    assert np.all(np.diff(ratio) < 0) and np.all(np.diff(phi) > 0)
    assert abs(ratio[0] - 1) < 1e-8 and phi[0] < 1e-4 and abs(phi[-1] - PI / 2) < 0.01
    assert abs(ratio[-1] * a[-1] ** 2 / 8 - 1) < 0.01                                         # inertia-limited 8/α²
    assert ch16.womersley_flow_ratio(0.0) == 1.0 and ch16.womersley_phase_lag(0.0) == 0.0
    assert rel(ch16.womersley_F2(1e-4), 1 - 1e-8j / 8) < 1e-15 and rel(ch16.womersley_F2(0.0), 1.0) < 1e-15


# ======================================================================================================================
# C11 — viscous waves in an elastic tube (16.116)–(16.137); D18–D21
# ======================================================================================================================
def test_elastic_tube_V2_derivation_wall_equations_and_modes():  # V2 — D19, D20 (★★★): (16.128)–(16.133)
    p, eta, xi, Eh, nh, k, a0, mu, du = sp.symbols("p eta xi Eh nu_h k a_0 mu du")
    # (16.126)–(16.127) with wall inertia and 2μ dv/dr dropped (h/a ≪ 1, a²k² ≪ 1): solve for p̂ and ξ̂
    e126 = p - Eh / (1 - nh ** 2) * (eta / a0 ** 2 + sp.I * nh * k / a0 * xi)
    e127 = -mu * du + Eh / (1 - nh ** 2) * (-k ** 2 * xi + sp.I * nh * k / a0 * eta)
    sol = sp.solve([e126, e127], [p, xi], dict=True)[0]
    assert sp.simplify(sol[p] - (Eh / a0 ** 2 * eta - sp.I * nh / (a0 * k) * mu * du)) == 0        # (16.128)
    assert sp.simplify(sol[xi] - (sp.I * nh / (k * a0) * eta - (1 - nh ** 2) / (Eh * k ** 2) * mu * du)) == 0   # (16.129)
    # (16.130)–(16.133): û = AJ₀(βr) + C solves û″ + û′/r + β²û = β²C; v̂ from (16.122)
    r, A, om, nu, rho, beta = sp.symbols("r A omega nu rho beta", positive=True)
    uh = A * sp.besselj(0, beta * r) + sp.Symbol("C")
    assert sp.simplify(sp.diff(uh, r, 2) + sp.diff(uh, r) / r + beta ** 2 * uh - beta ** 2 * sp.Symbol("C")) == 0
    vh = -sp.I * k * A / beta * sp.besselj(1, beta * r) - sp.I * k * sp.Symbol("C") * r / 2
    assert sp.simplify(sp.diff(r * vh, r) + sp.I * k * r * uh) == 0                               # (16.122)
    # the constant: C = (ik/μ)p̂/β² with β² = iω/ν, p̂ from (16.128), û′(a₀) = −AβJ₁(βa₀) → (16.131)
    C = sp.I * k / (rho * nu) * (Eh / a0 ** 2 * eta - sp.I * nh / (a0 * k) * rho * nu * (-A * beta * sp.besselj(1, beta * a0))) / (sp.I * om / nu)
    page = k / om * Eh / (rho * a0 ** 2) * eta - nh / (beta * a0) * A * sp.besselj(1, beta * a0)
    assert sp.simplify((C - page).subs(om, sp.I * 0 + beta ** 2 * nu / sp.I)) == 0


def test_elastic_tube_V2_derivation_rows_and_characteristic():  # V2 — D21 (★★★): (16.134), (16.135) → (16.136)
    eta, A, k, om, rho, a0, Eh, nh, b, J0, J1, X = sp.symbols("eta A k omega rho a_0 Eh nu_h beta J0 J1 X")
    nu = sp.I * om / b ** 2                        # β² = iω/ν
    mu = rho * nu
    bh = b * a0
    du = -A * b * J1                               # û′(a₀)
    C = k / om * Eh / (rho * a0 ** 2) * eta - nh / bh * A * J1           # (16.131) constant
    u_a = A * J0 + C
    v_a = -sp.I * k * A / b * J1 - sp.I * k ** 2 / om * Eh * eta / (rho * a0 ** 2) * a0 / 2 + sp.I * k * nh / bh * A * a0 / 2 * J1   # (16.133)
    xi = sp.I * nh / (k * a0) * eta - (1 - nh ** 2) / (Eh * k ** 2) * mu * du                       # (16.129)
    bc124 = sp.expand(u_a + sp.I * om * xi)
    bc125 = sp.expand(v_a + sp.I * om * eta)
    row134 = eta * (om / k * nh / a0 - k * Eh / (om * rho * a0 ** 2)) - A * (J0 + J1 * (sp.I * b * om * mu * (1 - nh ** 2) / (Eh * k ** 2) - nh / bh))
    row135 = eta * (1 - k ** 2 * Eh / (om ** 2 * 2 * rho * a0)) - A * J1 * (k / (om * b) - k * nh / (2 * om * b))
    assert sp.simplify(bc124 + row134) == 0                                              # (16.134) = −(16.124)
    assert sp.simplify(bc125 * sp.I / om + row135) == 0                                  # (16.135) = −(i/ω)(16.125)
    M = sp.Matrix([[sp.diff(row134, eta), sp.diff(row134, A)], [sp.diff(row135, eta), sp.diff(row135, A)]])
    det = sp.together(sp.expand(M.det()))
    Lr = 2 * bh * J0 / J1
    char = X ** 2 * (Lr - 4) + X * (4 * nh - 1 - Lr) + (1 - nh ** 2)                      # (16.136)
    detX = sp.simplify(det.subs(Eh, 2 * rho * a0 * om ** 2 * X / k ** 2))
    ratio = sp.simplify(detX * X / char)
    assert sp.simplify(sp.diff(ratio, X)) == 0 and sp.simplify(ratio) != 0              # det = (const/X)·(16.136)
    assert not sp.simplify(ratio).has(k)
    for al_, Xv, nv in ((8.0, 1.07 + 0.11j, 0.5), (3.0, 0.6 - 0.2j, 0.3)):
        bhat = complex(np.exp(0.25j * PI)) * al_                                           # β̂ = βa₀ = e^{iπ/4}α (T1)
        val = complex(char.subs({J0: jv(0, bhat), J1: jv(1, bhat), b: bhat, a0: 1, nh: nv, X: Xv}))
        assert rel(ch16.elastic_tube_characteristic(Xv, al_, nv), val) < 1e-12


def test_elastic_tube_V1_independent_4x4_determinant():  # V1 (independent-route cross-check) — C11: un-eliminated system vs (16.136) roots
    for al in (2.0, 5.0, 8.0, 20.0):
        roots = ch16.elastic_tube_roots(al)
        for key in ("pressure", "other"):
            X0 = roots[key]
            assert abs(ch16.elastic_tube_characteristic(X0, al)) < 1e-12
            assert rel(_root4(X0, al), X0) < 1e-10
        assert abs(ch16.elastic_tube_determinant_numeric(roots["pressure"], al)) < 1e-8
        assert abs(ch16.elastic_tube_determinant_numeric(1.3 * roots["pressure"], al)) > 1e3 * abs(
            ch16.elastic_tube_determinant_numeric(roots["pressure"], al))
    assert rel(_root4(ch16.elastic_tube_roots(8.0, poisson=0.3)["pressure"], 8.0, 0.3),
               ch16.elastic_tube_roots(8.0, poisson=0.3)["pressure"]) < 1e-10


def test_elastic_tube_V3_ode_determinant_series_start_order_six():  # V3 — C11 (elastic_tube_determinant_numeric)
    """The ODE determinant starts φ from its series truncated after s⁴ at s₀, so at an exact root of (16.136) its
    residual must fall as s₀⁶ (the first dropped term β̂⁶s⁶/2304) once the ODE tolerance is far below it."""
    for al in (2.0, 8.0):
        X0 = ch16.elastic_tube_roots(al)["pressure"]
        s0s = np.array([0.2, 0.1, 0.05, 0.025])
        res = [abs(ch16.elastic_tube_determinant_numeric(X0, al, s0=s, rtol=1e-13)) for s in s0s]
        assert abs(observed_order(s0s, res) - 6.0) < 0.15
        assert res[-1] < 1e-12
        far = abs(ch16.elastic_tube_determinant_numeric(1.3 * X0, al, s0=0.025, rtol=1e-13))
        assert far > 1e8 * res[-1]                                                  # small only at the root


def test_elastic_tube_V7_inviscid_limit_and_damping():  # V7 — C11 (trap T7), N68, N72
    al = np.array([25.0, 50.0, 100.0, 200.0])
    w = [ch16.elastic_tube_wave(x, c0=5.0, omega=7.0) for x in al]
    one_minus = [1 - d["c_over_c0"] for d in w]
    k21 = [d["k2_over_k1"] for d in w]
    assert abs(observed_order(1 / al, one_minus) - 1.0) < 0.1 and abs(observed_order(1 / al, k21) - 1.0) < 0.1
    assert abs(ch16.elastic_tube_roots(400.0)["pressure"] - 1) < 0.02
    d8 = ch16.elastic_tube_wave(8.0, c0=5.0, omega=7.0)
    assert rel(d8["c"], 5.0 * d8["c_over_c0"]) < 1e-15 and rel(d8["k1"], 7.0 / d8["c"]) < 1e-14
    assert rel(d8["damping_per_wavelength"], math.exp(-2 * PI * d8["k2"] / d8["k1"])) < 1e-14
    assert d8["c_over_c0"] < 1 and d8["k2"] > 0                                          # slower and damped
    assert rel(complex(d8["X"]) ** 0.5, complex(d8["k1"], d8["k2"]) * 5.0 / 7.0) < 1e-14  # k/k₀ = √X
    ev = ch16.elastic_tube_wave(8.0, E_ratio=1.0 + 0j)
    assert ev["X"] == ch16.elastic_tube_wave(8.0)["X"]
    ve = ch16.elastic_tube_wave(8.0, E_ratio=1 - 0.2j)                                    # viscoelastic wall (16.137)
    assert ve["k2_over_k1"] > ev["k2_over_k1"]
    other = ch16.elastic_tube_wave(8.0, branch="other")
    assert other["c_over_c0"] > 2.0                                                      # T7: the other root is the fast wall mode
    curve = ch16.elastic_tube_wave_curve(np.array([3.0, 8.0, 20.0]))
    assert close(curve["c_over_c0"], [ch16.elastic_tube_wave(x)["c_over_c0"] for x in (3.0, 8.0, 20.0)], rtol=1e-12)


def test_elastic_tube_V1_mode_shapes_satisfy_their_equations():  # V1 — D20 (16.130)–(16.133)
    a0, om, nu, E, h, rho, nh = 4e-3, 7.0, 3.3e-6, 0.5e6, 0.6e-3, 1060.0, 0.5
    al = a0 * math.sqrt(om / nu)
    c0 = math.sqrt(E * h / (2 * rho * a0))
    k = om / c0 * np.sqrt(ch16.elastic_tube_roots(al)["pressure"])
    Aamp, eta = 0.01, 2e-6
    r = np.linspace(0.2, 0.8, 7) * a0
    hh = a0 * 1e-4
    u = lambda x: ch16._elastic_tube_modes(x, Aamp, eta, k, om, nu, E, h, a0, rho, nh)[0]
    v = lambda x: ch16._elastic_tube_modes(x, Aamp, eta, k, om, nu, E, h, a0, rho, nh)[1]
    upp = (u(r + hh) - 2 * u(r) + u(r - hh)) / hh ** 2
    up = (u(r + hh) - u(r - hh)) / (2 * hh)
    dua = (u(a0) - u(a0 - hh)) / hh                                                     # û′(a₀), first order
    dua = (3 * u(a0) - 4 * u(a0 - hh) + u(a0 - 2 * hh)) / (2 * hh)
    rhs = 1j * k / (rho * nu) * E * h / a0 ** 2 * eta + nh / a0 * dua                    # (16.130)
    lhs = upp + up / r + 1j * om / nu * u(r)
    assert np.max(np.abs(lhs - rhs)) < 1e-5 * np.max(np.abs(rhs))
    drv = ((r + hh) * v(r + hh) - (r - hh) * v(r - hh)) / (2 * hh)
    assert np.max(np.abs(drv + 1j * k * r * u(r))) < 1e-6 * np.max(np.abs(k * r * u(r)))   # (16.122)


# ======================================================================================================================
# C12 — Murray's law (16.138)–(16.147); D22–D24
# ======================================================================================================================
def test_murray_V2_derivation_and_slip_16_140():  # V2 — D22, D23 (slip S15)
    a, Q, L, mu, K = sp.symbols("a Q L mu K", positive=True)
    E1 = 8 * mu * L / (sp.pi * a ** 4) * Q ** 2 + K * sp.pi * a ** 2 * L               # (16.138)
    dE = sp.diff(E1, a)
    assert sp.simplify(dE - (-32 * mu * L / sp.pi * Q ** 2 * a ** -5 + 2 * K * sp.pi * L * a)) == 0   # (16.139)
    a_opt = (16 * mu / (sp.pi ** 2 * K)) ** sp.Rational(1, 6) * Q ** sp.Rational(1, 3)              # corrected (16.140)
    assert sp.simplify(dE.subs(a, a_opt)) == 0 and sp.simplify(sp.diff(E1, a, 2).subs(a, a_opt)) != 0
    assert sp.simplify(sp.diff(E1, a, 2).subs(a, a_opt)).is_positive
    a_pr = (16 * mu / sp.pi ** 2 * K) ** sp.Rational(1, 6) * Q ** sp.Rational(1, 3)                 # printed (16.140)
    assert sp.simplify(dE.subs(a, a_pr)) != 0
    assert sp.simplify(E1.subs(a, a_opt) - sp.Rational(3, 2) * sp.pi * K * L * a_opt ** 2) == 0     # (16.141)
    pump = (8 * mu * L / (sp.pi * a ** 4) * Q ** 2).subs(a, a_opt)
    assert sp.simplify(pump - (K * sp.pi * a ** 2 * L).subs(a, a_opt) / 2) == 0                      # pumping = ½ upkeep
    Q1, Q2 = sp.symbols("Q1 Q2", positive=True)                                                     # D23 cube law
    c = (16 * mu / (sp.pi ** 2 * K)) ** sp.Rational(1, 6)
    assert sp.simplify((c * (Q1 + Q2) ** sp.Rational(1, 3)) ** 3 - (c * Q1 ** sp.Rational(1, 3)) ** 3 - (c * Q2 ** sp.Rational(1, 3)) ** 3) == 0
    # units: the corrected radius is a length, the printed one is not (K in W/m³)
    kq, muq, Qq = Q_(2000, "W/m**3"), Q_(3.5e-3, "Pa*s"), Q_(5e-6, "m**3/s")
    assert (16 * muq * Qq ** 2 / (PI ** 2 * kq)).check("[length]**6")                    # a⁶ of the corrected form
    assert not (16 * muq / PI ** 2 * kq * Qq ** 2).check("[length]**6")                  # the printed a⁶ is not a length⁶
    b = 2000.0
    assert rel(ch16.murray_radius(5e-6, 3.5e-3, b), (16 * 3.5e-3 / (PI ** 2 * b)) ** (1 / 6) * (5e-6) ** (1 / 3)) < 1e-14
    assert ch16.murray_radius(5e-6, 3.5e-3, 2 * b) < ch16.murray_radius(5e-6, 3.5e-3, b)            # dearer upkeep → thinner
    assert ch16.murray_radius(5e-6, 3.5e-3, 2 * b, printed_16_140=True) > ch16.murray_radius(5e-6, 3.5e-3, b, printed_16_140=True)


def test_murray_V1_grid_search_minimum_and_terms():  # V1 — C12
    Q, L, mu, b = 7.3e-6, 0.031, 3.5e-3, 1700.0
    aa = np.linspace(1e-3, 1e-2, 200001)
    cost = ch16.murray_cost(aa, Q, L, mu, b)
    a_num = aa[np.argmin(cost)]
    a_opt = ch16.murray_radius(Q, mu, b)
    assert abs(a_num - a_opt) < 2 * (aa[1] - aa[0])
    assert rel(ch16.murray_min_cost(Q, L, mu, b), ch16.murray_cost(a_opt, Q, L, mu, b)) < 1e-13
    t = ch16.murray_terms(a_opt, Q, L, mu, b)
    assert rel(t["pumping"], t["upkeep"] / 2) < 1e-12 and rel(t["total"], ch16.murray_cost(a_opt, Q, L, mu, b)) < 1e-14
    a_pr = ch16.murray_radius(Q, mu, b, printed_16_140=True)
    assert ch16.murray_cost(a_pr, Q, L, mu, b) > 1.5 * ch16.murray_min_cost(Q, L, mu, b)          # printed is far off
    bm = bench()["murray_wikipedia"]
    assert "r^3 = r1^3 + r2^3" in bm["what"]
    a1, a2 = ch16.murray_daughters(a_opt, 0.37)
    assert rel(a1 ** 3 + a2 ** 3, a_opt ** 3) < 1e-14
    assert rel(a1, ch16.murray_radius(0.37 * Q, mu, b)) < 1e-13                           # each daughter is optimal


def test_murray_V2_branch_angles_derivation():  # V2 — D24 (★★★): weighted Fermat point → (16.146)–(16.147)
    th, ph, a0, a1, a2 = sp.symbols("theta phi a0 a1 a2", positive=True)
    # stationarity of Σa_i²L_i under any displacement of B: a₀²e₀ = a₁²e₁ + a₂²e₂ (unit vectors along the branches)
    fx = a0 ** 2 - a1 ** 2 * sp.cos(th) - a2 ** 2 * sp.cos(ph)
    fy = a1 ** 2 * sp.sin(th) - a2 ** 2 * sp.sin(ph)
    # steps 7–9 with cos and sin as polynomial symbols (c² + s² = 1): no trigonometric solve, deterministic and fast
    ct, st, cp, sp_ = sp.symbols("c_t s_t c_p s_p", real=True)
    eq = sp.expand((a0 ** 2 - a1 ** 2 * ct) ** 2 + (a1 ** 2 * st) ** 2 - a2 ** 4).subs(st ** 2, 1 - ct ** 2)   # |a₂²e₂|² = a₂⁴
    cth = sp.solve(eq, ct)[0]
    assert sp.simplify(cth - (a0 ** 4 + a1 ** 4 - a2 ** 4) / (2 * a0 ** 2 * a1 ** 2)) == 0                  # (16.146) cos θ
    eqp = sp.expand((a0 ** 2 - a2 ** 2 * cp) ** 2 + (a2 ** 2 * sp_) ** 2 - a1 ** 4).subs(sp_ ** 2, 1 - cp ** 2)
    cph = sp.solve(eqp, cp)[0]
    assert sp.simplify(cph - (a0 ** 4 - a1 ** 4 + a2 ** 4) / (2 * a0 ** 2 * a2 ** 2)) == 0                  # (16.146) cos φ
    # the third relation: |a₁²e₁ + a₂²e₂|² = a₀⁴ with e₁·e₂ = cos θ cos φ − sin θ sin φ = cos(θ + φ) = C
    C = sp.symbols("C", real=True)
    eq3 = sp.expand(a1 ** 4 + a2 ** 4 + 2 * a1 ** 2 * a2 ** 2 * C - a0 ** 4)
    assert sp.simplify(sp.solve(eq3, C)[0] - (a0 ** 4 - a1 ** 4 - a2 ** 4) / (2 * a1 ** 2 * a2 ** 2)) == 0
    # (16.147): a₂³ = a₀³ − a₁³
    c147 = cth.subs(a2, (a0 ** 3 - a1 ** 3) ** sp.Rational(1, 3))
    assert sp.simplify(c147 - (a0 ** 4 + a1 ** 4 - (a0 ** 3 - a1 ** 3) ** sp.Rational(4, 3)) / (2 * a0 ** 2 * a1 ** 2)) == 0
    # numerical: the coded angles satisfy both force components
    for s in (0.5, 0.3, 0.12):
        A0r = 1.7e-3
        b1, b2 = ch16.murray_daughters(A0r, s)
        ang = ch16.murray_angles(A0r, b1, b2)
        vals = {a0: A0r, a1: b1, a2: b2, th: ang["theta"], ph: ang["phi"]}
        assert abs(float(fx.subs(vals))) < 1e-15 and abs(float(fy.subs(vals))) < 1e-15
        assert rel(ch16.murray_theta(A0r, b1), ang["theta"]) < 1e-12
        c = (A0r ** 4 - b1 ** 4 - b2 ** 4) / (2 * b1 ** 2 * b2 ** 2)
        assert rel(math.cos(ang["total"]), c) < 1e-12


def test_murray_V1_numerical_branch_point_reproduces_angles():  # V1 (independent-route cross-check) — D24 check, C12 (bifurcation_point)
    A, C, D = np.array([0.0, 0.0]), np.array([0.09, 0.037]), np.array([0.083, -0.051])
    for s in (0.5, 0.31, 0.17):
        a0 = 1.3e-3
        a1, a2 = ch16.murray_daughters(a0, s)
        bp = ch16.bifurcation_point(A, C, D, a0, a1, a2)
        ang = ch16.murray_angles(a0, a1, a2)
        assert abs(bp["theta"] - ang["theta"]) < 1e-6 and abs(bp["phi"] - ang["phi"]) < 1e-6
        assert rel(bp["cost"], ch16.bifurcation_cost(bp["B"], A, C, D, a0, a1, a2)) < 1e-14
        for dB in (np.array([1e-4, 0]), np.array([0, 1e-4]), np.array([-7e-5, 7e-5])):
            assert ch16.bifurcation_cost(bp["B"] + dB, A, C, D, a0, a1, a2) > bp["cost"]


def test_murray_V7_symmetric_and_degenerate_forks():  # V7 — C12
    th = ch16.murray_symmetric_angle()
    assert rel(th, math.acos(2 ** (-1 / 3))) < 1e-15 and abs(math.degrees(th) - 37.467) < 1e-3
    a1, a2 = ch16.murray_daughters(1.0, 0.5)
    ang = ch16.murray_angles(1.0, a1, a2)
    assert rel(ang["theta"], th) < 1e-14 and rel(ang["phi"], th) < 1e-14 and ang["valid"] is True
    assert ch16.murray_theta(1.0, 1e-4 ** (1 / 3) * 1.0) > 1.4                           # a tiny branch leaves near 90°
    assert ch16.murray_theta(1.0, (1 - 1e-6) ** (1 / 3)) < 0.01                           # the main branch goes straight
    bad = ch16.murray_angles(1.0, 0.2, 0.2)                                              # not a cube-law fork: no optimum
    assert bad["valid"] is False and math.isnan(bad["theta"])
    assert rel(ch16.murray_cost(2e-3, 3e-6, 0.05, 3.5e-3, 1900.0),
               8 * 3.5e-3 / PI * 0.05 * 3e-6 ** 2 / 2e-3 ** 4 + PI * 1900.0 * 0.05 * 2e-3 ** 2) < 1e-14   # V5 Wikipedia form


# ======================================================================================================================
# C13 — reflection and transmission at a junction (16.148)–(16.163); D25, D26
# ======================================================================================================================
def test_junction_V1_linear_solve_and_limits():  # V1 — C13 (trap T18)
    rng = np.random.default_rng(16)
    for _ in range(20):
        Y1 = rng.uniform(0.1, 3)
        Yd = rng.uniform(0.05, 2, size=rng.integers(1, 4))
        # unknowns P_R, P_T (per unit P_I): continuity (16.148) P_I + P_R = P_T; flux (16.152) Y₁(P_I − P_R) = ΣY_j P_T
        M = np.array([[1.0, -1.0], [-Y1, -Yd.sum()]])
        PR, PT = np.linalg.solve(M, [-1.0, -Y1])
        R, T = ch16.junction_coefficients(Y1, Yd)
        assert abs(R - PR) < 1e-14 and abs(T - PT) < 1e-14 and abs(T - 1 - R) < 1e-14
    assert ch16.junction_coefficients(2.0, [1.0, 1.0]) == (0.0, 1.0)                     # matched
    R, _ = ch16.junction_coefficients(1.0, [1e-12])
    assert abs(R - 1) < 1e-11                                                           # closed end: R = +1
    R, _ = ch16.junction_coefficients(1.0, [1e12])
    assert abs(R + 1) < 1e-11                                                           # open end: R = −1
    A, rho, c = 3.1e-4, 1060.0, 5.7
    assert rel(ch16.tube_admittance(A, rho, c) * ch16.tube_impedance(A, rho, c), 1.0) < 1e-15
    assert (Q_(A, "m**2") / (Q_(rho, "kg/m**3") * Q_(c, "m/s"))).check("[length]**3/([pressure]*[time])")


def test_junction_V4_energy_is_conserved():  # V4 — D26 (16.161)–(16.163): R² + (ΣY/Y₁)T² = 1
    rng = np.random.default_rng(1601)
    for _ in range(200):
        Y1 = rng.uniform(0.01, 10)
        Yd = rng.uniform(0.01, 10, size=rng.integers(1, 5))
        e = ch16.junction_energy(Y1, Yd)
        assert abs(e["total"] - 1) < 1e-14
        R, T = ch16.junction_coefficients(Y1, Yd)
        assert rel(e["reflected"], R ** 2) < 1e-15


def test_junction_V2_derivation():  # V2 — D25, D26
    Y1, Ys, PI_, PR, PT = sp.symbols("Y1 Ys P_I P_R P_T", positive=True)
    sol = sp.solve([PI_ + PR - PT, Y1 * (PI_ - PR) - Ys * PT], [PR, PT], dict=True)[0]
    R = sp.simplify(sol[PR] / PI_)
    T = sp.simplify(sol[PT] / PI_)
    assert sp.simplify(R - (Y1 - Ys) / (Y1 + Ys)) == 0 and sp.simplify(T - 2 * Y1 / (Y1 + Ys)) == 0   # (16.153)
    Z1, Zs = 1 / Y1, 1 / Ys
    assert sp.simplify(R - (1 / Z1 - 1 / Zs) / (1 / Z1 + 1 / Zs)) == 0                                  # (16.154)
    assert sp.simplify(R ** 2 + Ys / Y1 * T ** 2 - 1) == 0                                              # energy identity


def test_junction_V1_grid_pulse_matches_theory():  # V1 (independent-route cross-check) — C13 (network_pulse_run, exact C = 1 transport)
    for Yd in ([0.3, 0.45], [1.1, 0.9, 0.4]):
        Y1 = 1.0
        out = ch16.network_pulse_run(0.8, [0.6] * len(Yd), 5.0, [6.5] + [5.5] * (len(Yd) - 1), Y1, Yd, width=0.012, nx=400)
        R, T = ch16.junction_coefficients(Y1, Yd)
        assert abs(out["R_measured"] - R) < 1e-3 and abs(out["T_measured"] - T) < 1e-3
        assert out["R_theory"] == R and out["T_theory"] == T
        out2 = ch16.network_pulse_run(0.8, [0.6] * len(Yd), 5.0, [6.5] + [5.5] * (len(Yd) - 1), Y1, Yd, width=0.012,
                                      nx=200)
        assert abs(out2["R_measured"] - out["R_measured"]) < 1e-12                      # C = 1 transport is exact:
        assert abs(out2["T_measured"] - out["T_measured"]) < 1e-12                      # nothing to refine
    f = lambda s: np.exp(-((np.asarray(s) - 0.1) / 0.02) ** 2)
    x, t = np.linspace(-1, 0, 21), 0.17
    w = ch16.junction_waveforms(x, t, f, 5.0, -0.3, Y1=2.0, P_I=900.0)
    assert close(w["p"], 900 * (f(t - x / 5) - 0.3 * f(t + x / 5)), rtol=1e-14, atol=1e-12)
    assert close(w["Q"], 2 * 900 * (f(t - x / 5) + 0.3 * f(t + x / 5)), rtol=1e-14, atol=1e-12)
    g = ch16.junction_from_geometry(4e-3, list(ch16.murray_daughters(4e-3, 0.5)), 5.0, 5.0, 1060.0)
    assert rel(g["ratio"], 2 ** (1 / 3)) < 1e-13 and g["R"] < 0                          # Murray fork with equal c: open-end-like
    assert abs(g["reflected"] + g["transmitted"] - 1) < 1e-14


# ======================================================================================================================
# C14 — Dean flow in a gently curved pipe (16.164)–(16.182); D27; N86–N93
# ======================================================================================================================
def test_dean_V2_derivation_leading_orders_and_flux_coefficient():  # V2 — D27 (★★★), steps 1–14 re-run by hand
    r, al = sp.symbols("r alpha", positive=True)
    lap = lambda F: sp.diff(F, r, 2) + sp.diff(F, r) / r + sp.diff(F, al, 2) / r ** 2      # (16.175)
    w0 = 2 * (1 - r ** 2)
    assert sp.simplify(lap(w0) + 8) == 0 and w0.subs(r, 1) == 0                              # steps 1–2
    assert sp.integrate(w0 * 2 * r, (r, 0, 1)) == 1                                          # mean speed 1
    psi0 = r / 144 * (1 - r ** 2) ** 2 * (4 - r ** 2) * sp.sin(al)                           # step 8
    forcing = -2 * w0 * (sp.sin(al) * sp.diff(w0, r) + sp.cos(al) / r * sp.diff(w0, al)) / 2  # steps 3–4 (16.174)/2
    assert sp.simplify(lap(lap(psi0)) - forcing) == 0 and sp.simplify(forcing - 8 * r * (1 - r ** 2) * sp.sin(al)) == 0
    assert psi0.subs(r, 1) == 0 and sp.diff(psi0, r).subs(r, 1) == 0                         # (16.176)
    w1 = r * (1 - r ** 2) * (19 - 21 * r ** 2 + 9 * r ** 4 - r ** 6) / 5760 * sp.cos(al)      # step 11 (plus sign)
    rhs1 = (sp.diff(psi0, al) * sp.diff(w0, r) - sp.diff(psi0, r) * sp.diff(w0, al)) / (2 * r)   # steps 9–10 (16.173)
    assert sp.simplify(lap(w1) - rhs1) == 0 and w1.subs(r, 1) == 0
    # steps 12–14: the α-mean of w₂ from ψ₀ and w₁ only; ∇²w̄₂ = mean of (1/2r)(ψ₀_α w₁_r − ψ₀_r w₁_α)
    S = sp.simplify(sp.integrate((sp.diff(psi0, al) * sp.diff(w1, r) - sp.diff(psi0, r) * sp.diff(w1, al)) / (2 * r),
                                 (al, 0, 2 * sp.pi)) / (2 * sp.pi))
    c1 = sp.symbols("c1")
    rwp = sp.integrate(r * S, (r, 0, r))                                                     # r w̄₂′ (regular on the axis)
    w2bar = sp.integrate(rwp / r, r) + c1
    w2bar = w2bar.subs(c1, sp.solve(w2bar.subs(r, 1), c1)[0])
    m2 = sp.integrate(w2bar * 2 * r, (r, 0, 1))
    f, g = r * (1 - r ** 2) ** 2 * (4 - r ** 2) / 144 / r * r, r * (1 - r ** 2) * (19 - 21 * r ** 2 + 9 * r ** 4 - r ** 6) / 5760
    assert sp.simplify(sp.integrate(r * (f / r) * r * g, (r, 0, 1)) - sp.Rational(1541, 1045094400)) == 0   # step 13 ∫r f g
    c2 = sp.nsimplify(-m2 * 576 ** 2 / 4)                                                    # step 14, K = 2κ²
    assert c2 == sp.Rational(1541, 50400)
    ds = ch16.dean_series_sympy(2)
    assert ds["flux_coefficient"] == sp.Rational(1541, 50400) and ds["c2"] == Fraction(1541, 50400)
    assert rel(ch16.dean_flux_ratio(300.0, exact=True), 1 - float(c2) * (300 / 576) ** 2 + 0.0120 * (300 / 576) ** 4) < 1e-15
    assert sp.simplify(ds["psi0"].subs(sp.Symbol("r", positive=True), r) - psi0) == 0 or abs(
        float(ch16.dean_secondary_streamfunction(0.43, 1.1)) - float(psi0.subs({r: 0.43, al: 1.1}))) < 1e-15
    for rr, aa in ((0.2, 0.3), (0.5, 0.0), (0.83, 2.6)):
        assert abs(float(ch16.dean_secondary_streamfunction(rr, aa)) - float(psi0.subs({r: rr, al: aa}))) < 1e-15
        assert abs(float(ch16.dean_axial_correction(rr, aa)) - float(w1.subs({r: rr, al: aa}))) < 1e-15


def test_dean_vorticity_V2_printed_16_174_nonlinear_sign_is_inconsistent():  # V2 — N89 (new slip, reported)
    """Eliminate p₁ from the book's own (16.167)–(16.168) with (16.172): the vorticity-advection term comes out with the
    sign opposite to the one printed in (16.174) (the centrifugal right side and the viscous term agree)."""
    r, a, k = sp.symbols("r alpha kappa", positive=True)
    psi, w, p1 = sp.Function("psi")(r, a), sp.Function("w")(r, a), sp.Function("p1")(r, a)
    u, v = sp.diff(psi, a) / r, -sp.diff(psi, r)                                              # (16.172)
    zeta = sp.diff(v, r) + v / r - sp.diff(u, a) / r
    Er = (u * sp.diff(u, r) + v / r * sp.diff(u, a) - v ** 2 / r - w ** 2 * sp.cos(a)) - (-sp.diff(p1, r) - 2 / k / r * sp.diff(zeta, a))   # (16.167)
    Ea = (u * sp.diff(v, r) + v / r * sp.diff(v, a) + u * v / r + w ** 2 * sp.sin(a)) - (-sp.diff(p1, a) / r + 2 / k * sp.diff(zeta, r))   # (16.168)
    curl = sp.simplify((sp.diff(r * Ea, r) - sp.diff(Er, a)) / r)
    assert not curl.has(p1)
    lap = lambda F: sp.diff(F, r, 2) + sp.diff(F, r) / r + sp.diff(F, a, 2) / r ** 2
    adv = lambda g: (sp.diff(psi, r) * sp.diff(g, a) - sp.diff(psi, a) * sp.diff(g, r)) / r
    rhs = -2 * w * (sp.sin(a) * sp.diff(w, r) + sp.cos(a) / r * sp.diff(w, a))
    printed = 2 / k * lap(lap(psi)) - adv(lap(psi)) - rhs                                    # (16.174) as printed
    consistent = 2 / k * lap(lap(psi)) + adv(lap(psi)) - rhs                                 # nonlinear sign flipped
    assert sp.simplify(sp.expand(curl - consistent)) == 0
    assert sp.simplify(sp.expand(curl - printed)) != 0 and sp.simplify(sp.expand(curl + printed)) != 0


def test_dean_V3_numeric_navier_stokes_leading_orders():  # V3 — C14: ψ₀, w₁ and c₂ from the independent solver
    op = _dean_ops(21, 16)
    n = op["n2"] * op["M"]
    r = op["r"]
    i = int(np.argmin(np.abs(r - 0.5)))
    d_w1, d_psi = [], []
    Ks = np.array([0.25, 0.5, 1.0])
    for Kl in Ks:                                    # κ = K_lec/4; v = w/8; u_book = 8 u_lec
        U = _dean_solve(Kl, op)
        v = U[:n].reshape(op["n2"], op["M"])
        w1n = (v[:, 0] - v[:, op["M"] // 2]) / 2 / (Kl ** 2 / 128)
        ur = (op["Dth"] @ U[n:2 * n]).reshape(op["n2"], op["M"])[:, 0] / r          # radial velocity on the outer ray
        f = ch16.dean_secondary_streamfunction(r, PI / 2) / r
        d_w1.append(abs(w1n[i] / ch16.dean_axial_correction(r[i], 0.0) - 1))
        d_psi.append(abs(ur[i] / (Kl / 32) / f[i] - 1))
    assert d_w1[0] < 1e-9 and d_psi[0] < 1e-5                                     # sign and size of w₁ (outer bend) and ψ₀
    assert abs(observed_order(Ks, d_psi) - 2.0) < 0.15                            # the next order is κ²
    c2, c4 = _dean_coefficients_numeric(sJ=1.0)
    assert rel(c2, 1541 / 50400) < 2e-6                                          # (16.178): 0.0306 = 1541/50400
    assert rel(ch16.dean_flux_ratio(576.0, exact=True), 1 - 1541 / 50400 + float(ch16.DEAN_BOOK["c4"])) < 1e-15
    _, c4_printed = _dean_coefficients_numeric(sJ=-1.0)                         # the equations as printed
    assert rel(c4_printed, float(Fraction(821328579941, 67805337600000))) < 2e-5  # = the module's exact series


def test_dean_V3_spectral_resolution_study_and_fit_truncation():  # V3 — C14: the independent solver is converged
    """Resolution study of the Chebyshev–Fourier Dean solver at the largest sampled Dean number (x = K/576 = 0.35):
    the flux ratio converges geometrically in the radial resolution N and in M; then the c₂, c₄ fit converges as more
    terms of the series are fitted, towards the exact 1541/50400 and 6471982981/542442700800."""
    K_lec = math.sqrt(8 * 576 * 0.35)
    F = {}
    for N in (9, 13, 17, 21, 25):
        op = _dean_ops(N, 16)
        F[N] = _dean_flux_ratio_numeric(_dean_solve(K_lec, op), op)
    d = [abs(F[N] - F[25]) for N in (9, 13, 17, 21)]
    assert all(d[i + 1] < 1e-1 * d[i] for i in range(3)) and d[-1] < 1e-10       # spectral (geometric) in N
    for M in (8, 12):
        op = _dean_ops(21, M)
        assert abs(_dean_flux_ratio_numeric(_dean_solve(K_lec, op), op) - F[21]) < 1e-9   # converged in M
    xs, Fs = _dean_flux_samples(1.0, 21, 16)
    c2x, c4x = 1541 / 50400, 6471982981 / 542442700800
    e2, e4 = [], []
    for nt in (3, 4, 5, 6):                                                    # (7 terms over-fits: floor ~1e-6)
        c2, c4 = _dean_fit(xs, Fs, nt)
        e2.append(abs(c2 / c2x - 1))
        e4.append(abs(c4 / c4x - 1))
    assert all(e4[i + 1] < 0.3 * e4[i] for i in range(3)) and e4[-1] < 1e-6        # fit truncation → exact c₄
    assert e2[-1] < 1e-8


def test_dean_flux_V3_order4_coefficient_is_the_navier_stokes_value():  # V3 + V2 — C14 (16.178) K⁴ coefficient
    """The K⁴ flux coefficient of (16.178) from the Navier–Stokes-consistent Dean equations is 0.0119312 (spectral
    solver here; exact 6471982981/542442700800 from the module's own series with the (16.174) advection sign corrected);
    the module's DEAN_C4_DERIVED = 0.0121130 expands the misprinted (16.174)."""
    _, c4 = _dean_coefficients_numeric(sJ=1.0)
    exact_ns = Fraction(6471982981, 542442700800)
    assert rel(c4, float(exact_ns)) < 2e-5                                        # the two independent routes agree
    assert rel(float(ch16.DEAN_C4_DERIVED), c4) < 1e-4, (
        f"DEAN_C4_DERIVED = {float(ch16.DEAN_C4_DERIVED):.7f}, Navier–Stokes value {c4:.7f} "
        f"({100 * (float(ch16.DEAN_C4_DERIVED) / c4 - 1):+.2f} %); printed (16.174) has the advection sign flipped")


def test_dean_flux_V2_V3_printed_16_174_option_reproduces_the_wrong_variant():  # V2 + V3 — C14, slip S34 (loop 2)
    """The corrected series and DEAN_C4_DERIVED are the exact Navier–Stokes value; ``printed_16_174=True`` reproduces
    the misprinted variant 821328579941/67805337600000 = 0.0121130, which must FAIL the spectral NS value."""
    ok = ch16.dean_series_sympy(4)
    pr = ch16.dean_series_sympy(4, printed_16_174=True)
    assert ok["c4"] == Fraction(6471982981, 542442700800) == ch16.DEAN_C4_DERIVED      # exact rational
    assert pr["c4"] == Fraction(821328579941, 67805337600000)
    assert abs(float(pr["c4"]) - 0.0121130) < 5e-8 and abs(float(ok["c4"]) - 0.0119312) < 5e-8
    # the leading orders are untouched by the sign (it first enters ψ₁): ψ₀, w₁, c₂ identical
    assert ok["c2"] == pr["c2"] == Fraction(1541, 50400)
    assert sp.simplify(ok["psi0"] - pr["psi0"]) == 0 and sp.simplify(ok["w1"] - pr["w1"]) == 0
    assert sp.simplify(ok["psi"][1] - pr["psi"][1]) != 0                                 # ψ₁ is where they part
    # V3: the independent spectral Navier–Stokes solve separates them
    _, c4_ns = _dean_coefficients_numeric(sJ=1.0)
    assert rel(float(ok["c4"]), c4_ns) < 2e-5
    assert rel(float(pr["c4"]), c4_ns) > 1e-2                                           # wrong variant fails (1.5 %)
    # dean_flux_ratio(exact_c4=True) carries the corrected coefficient
    K = 400.0
    q = K / 576
    assert rel(ch16.dean_flux_ratio(K, exact=True, exact_c4=True), 1 - 1541 / 50400 * q ** 2 + float(ok["c4"]) * q ** 4) < 1e-15


def test_dean_V1_conventions_and_resistance():  # V1 — N86–N93 (trap T8, slips S18–S20)
    for kap in (0.3, 7.0, 31.0):
        K = ch16.dean_convert(kap, from_="kappa", to="K")
        D = ch16.dean_convert(kap, from_="kappa", to="D")
        assert rel(K, 2 * kap ** 2) < 1e-15 and rel(D, 4 * math.sqrt(2) * kap) < 1e-15 and rel(K, D ** 2 / 16) < 1e-14
        assert rel(ch16.dean_convert(K, from_="K", to="kappa"), kap) < 1e-14
        assert rel(ch16.dean_convert(D, from_="D", to="De"), kap) < 1e-14
    assert rel(ch16.dean_number(0.09, 0.3, 4e-3, 3.3e-6, convention="kappa"), math.sqrt(0.09) * 2 * 4e-3 * 0.3 / 3.3e-6) < 1e-14
    assert rel(ch16.curvature_ratio(4e-3, 0.05), 0.08) < 1e-15
    for K in (50.0, 200.0, 400.0):
        x = K / 576
        F = ch16.dean_flux_ratio(K)
        assert rel(F, 1 - 0.0306 * x ** 2 + 0.0120 * x ** 4) < 1e-15
        assert rel(ch16.dean_flux_ratio(D=ch16.dean_convert(K, from_="K", to="D")), F) < 1e-13
        lam = ch16.curved_resistance_ratio(K)
        assert abs(1 / F - lam) < (abs(0.0306 ** 2 - 0.0120 + 0.0110) + 0.002) * x ** 4    # (16.182) = 1/(16.178) to O(K⁴)
    rel_ = ch16.dean_series_reliability(300.0, threshold=0.25)
    assert rel(rel_["ratio"], 0.0120 * (300 / 576) ** 2 / 0.0306) < 1e-14 and rel_["reliable"]
    assert not ch16.dean_series_reliability(700.0, threshold=0.25)["reliable"]
    rr, aa, K = 0.47, 0.6, 180.0
    assert rel(ch16.dean_axial_velocity(rr, aa, K), 2 * (1 - rr ** 2) + K / 2 * ch16.dean_axial_correction(rr, aa)) < 1e-15
    u, v = ch16.dean_secondary_velocity(rr, aa, K)
    h = 1e-6
    kap = math.sqrt(K / 2)
    dpa = (ch16.dean_secondary_streamfunction(rr, aa + h) - ch16.dean_secondary_streamfunction(rr, aa - h)) / (2 * h)
    dpr = (ch16.dean_secondary_streamfunction(rr + h, aa) - ch16.dean_secondary_streamfunction(rr - h, aa)) / (2 * h)
    assert rel(u, kap * dpa / rr) < 1e-8 and rel(v, -kap * dpr) < 1e-8                   # (16.172)
    assert ch16.dean_secondary_velocity(0.3, 0.0, K)[0] > 0                                # core moves to the outer bend


# ======================================================================================================================
# C15 — collapsible tubes: Shapiro's equations (16.185)–(16.206); D28–D30; Starling resistor
# ======================================================================================================================
def test_collapsible_V2_derivation_shapiro_equations():  # V2 — D29, D30 (★★★): (16.204), (16.205) from the differentials
    al, rho, Kp, A0, u = sp.symbols("alpha rho K_p A_0 u", positive=True)
    dal, dA0, dpe, gdz, RQ, dKp = sp.symbols("dalpha dA0 dpe gdz RQdx dKp", real=True)
    P = sp.Function("P")(al)
    Pa, Paa = sp.diff(P, al), sp.diff(P, al, 2)
    c2 = al / rho * Kp * Pa                                                                    # (16.198)
    S2 = u ** 2 / c2
    dA_A = dal / al + dA0 / A0                                                                 # (16.196)
    dp = dpe + rho * c2 * dal / al + P * dKp                                                   # (16.200), P = P(α)
    # (16.194) with u du = −u² dA/A (16.190):  dp + ρg dz + R Q dx − ρu² dA/A = 0
    da = sp.solve(sp.Eq(dp + rho * gdz + RQ - rho * u ** 2 * dA_A, 0), dal)[0]
    rhs204 = S2 * dA0 / A0 - (dpe + rho * gdz + RQ + P * dKp) / (rho * c2)
    assert sp.simplify((1 - S2) * da / al - rhs204) == 0                                     # (16.204)
    two_dc_c = (1 + al * Paa / Pa) * da / al + dKp / Kp                                       # (16.201), P = P(α)
    dS2_S2 = -2 * da / al - 2 * dA0 / A0 - two_dc_c                                           # (16.202)
    M = 3 + al * Paa / Pa                                                                      # (16.206)
    rhs205 = (dA0 / A0 * (-2 + (2 - M) * S2) + M / (rho * c2) * (dpe + rho * gdz + RQ)
              + dKp / (rho * c2) * (M * P - (1 - S2) * al * Pa))
    assert sp.simplify((1 - S2) * dS2_S2 - rhs205) == 0                                       # (16.205)
    # at S = 1 the two right sides are proportional: smooth passage needs both to vanish together
    assert sp.simplify((rhs205 + M * rhs204).subs(u, sp.sqrt(c2))) == 0
    # (16.201) from c² = (α/ρ)K_p ∂P/∂α with P = P(α, x): the inner derivative is ∂P/∂α (slip S23)
    x = sp.symbols("x")
    Px = sp.Function("P")(al, x)
    Kx = sp.Function("K")(x)
    cc = al / rho * Kx * sp.diff(Px, al)
    dln = sp.diff(cc, al) / cc * dal + sp.diff(cc, x) / cc * sp.Symbol("dx")
    corr = (1 + al * sp.diff(Px, al, 2) / sp.diff(Px, al)) * dal / al + sp.diff(Kx, x) / Kx * sp.Symbol("dx") \
        + al * Kx / (rho * cc) * sp.diff(sp.diff(Px, al), x) * sp.Symbol("dx")
    printed = corr - al * Kx / (rho * cc) * sp.diff(sp.diff(Px, al), x) * sp.Symbol("dx") \
        + al * Kx / (rho * cc) * sp.diff(sp.diff(Px, x), x) * sp.Symbol("dx")
    assert sp.simplify(dln - corr) == 0 and sp.simplify(dln - printed) != 0
    # the module's own engine agrees, and its printed option fails
    d = ch16.collapsible_differentials_sympy()
    assert d["eq204_residual"] == 0 and d["eq205_residual"] == 0 and sp.simplify(d["dc_difference"]) == 0
    assert sp.simplify(d["M_power_law"] - (2 - sp.Symbol("n", positive=True))) == 0
    assert sp.simplify(ch16.collapsible_differentials_sympy(printed_16_201=True)["dc_difference"]) != 0


def test_collapsible_V1_tube_law_wave_speed_and_M():  # V1 — D28 (16.185)–(16.189), (16.206) (trap T9)
    Kp, rho, A0 = 230.0, 1050.0, 1.3e-4
    al = np.array([0.15, 0.4, 0.8, 0.97, 1.0])                 # the fitted collapse branch of (16.188): α ≤ 1 only
    for form in ("shapiro", "pure"):
        P = ch16.collapsible_tube_law(al, n=1.5, form=form)
        assert close(ch16.alpha_from_tube_law(P, 1.5, form), al, rtol=1e-13)                    # round trip
        h = al * 1e-6
        # second-order one-sided (backward) stencil for ∂P/∂α everywhere: at α = 1 a central stencil would step
        # outside the fitted range (α > 1 is NaN by design)
        dP = (3 * ch16.collapsible_tube_law(al, form=form) - 4 * ch16.collapsible_tube_law(al - h, form=form)
              + ch16.collapsible_tube_law(al - 2 * h, form=form)) / (2 * h)
        c = ch16.collapsible_wave_speed(al, Kp, rho)
        assert close(c ** 2, al / rho * Kp * dP, rtol=1e-8)                                        # (16.198)
        assert close(c, ch16.bramwell_hill_speed(al * A0, Kp * dP / A0, rho), rtol=1e-8)            # (16.185)
        inner = al[al < 1]                                       # shapiro_M's central differences need α(1 + h) ≤ 1
        assert close(ch16.shapiro_M(inner, 1.5, tube_law=lambda q: ch16.collapsible_tube_law(q, form=form)), 0.5, rtol=1e-5)
    assert np.all(ch16.collapsible_tube_law(al[al < 1]) < 0)                                      # P < 0 when collapsed
    assert close(ch16.collapsible_tube_law(al), 1 - al ** -1.5, rtol=1e-15)                        # (16.188)
    assert close(ch16.shapiro_M(al, 1.5), 0.5, rtol=1e-15) and close(ch16.shapiro_M(al, 2.5), -0.5, rtol=1e-15)
    assert rel(ch16.choking_flux(0.4, A0, Kp, rho), 0.4 * A0 * ch16.collapsible_wave_speed(0.4, Kp, rho)) < 1e-15
    assert rel(ch16.speed_index(0.3, 0.75), 0.4) < 1e-15
    assert rel(ch16.bending_stiffness_parameter(2e5, 3e-4, 4e-3, 0.45, coeff=0.7), 0.7 * 2e5 * (3e-4 / 4e-3) ** 3 / (1 - 0.45 ** 2)) < 1e-14
    d0, mu = 2 * math.sqrt(A0 / PI), 3.5e-3
    for a_ in (0.3, 1.0):
        uu = 0.07
        assert rel(ch16.wall_friction_term(uu, a_, d0, rho, mu, regime="laminar"),
                   ch16.friction_R_per_length(a_, A0, mu) * uu * a_ * A0) < 1e-13                # τ_w s/A = R′Q
    assert rel(ch16.friction_R_per_length(1.0, PI * 2e-3 ** 2, mu), 8 * mu / (PI * 2e-3 ** 4)) < 1e-14   # Poiseuille
    assert rel(ch16.wall_friction_term(0.5, 1.0, d0, rho, mu, regime="turbulent", f=0.006), 0.5 * rho * 0.25 * 4 * 0.006 / d0) < 1e-14


def test_collapsible_V7_tube_law_range_guard_alpha_le_1():  # V7 — D28, (16.188) valid for α ≤ 1 (review Must-fix 1)
    """(16.188) is Shapiro's fit to the collapse branch (Fig. 16.19).  Scalars α > 1 warn and return the extrapolated
    value; array entries α > 1 are NaN with one RuntimeWarning; α ≤ 1 (and α = 1 + round-off) is silent."""
    Kp, rho, A0 = 230.0, 1050.0, 1.3e-4
    # scalar α > 1: RuntimeWarning, extrapolated value still returned
    with pytest.warns(RuntimeWarning, match="16.188"):
        assert rel(ch16.collapsible_tube_law(1.7), 1 - 1.7 ** -1.5) < 1e-15
    with pytest.warns(RuntimeWarning, match="16.188"):
        assert rel(ch16.collapsible_tube_law(1.7, form="pure"), -1.7 ** -1.5) < 1e-15
    with pytest.warns(RuntimeWarning, match="16.188"):
        P, dP, d2P = ch16.collapsible_tube_law_derivatives(1.7)
    assert rel(dP, 1.5 * 1.7 ** -2.5) < 1e-15 and rel(d2P, -3.75 * 1.7 ** -3.5) < 1e-15
    with pytest.warns(RuntimeWarning, match="16.188"):
        assert rel(ch16.collapsible_wave_speed(1.7, Kp, rho), math.sqrt(1.5 * Kp * 1.7 ** -1.5 / rho)) < 1e-15
    with pytest.warns(RuntimeWarning, match="16.188"):                     # P > 0 ⇔ α > 1 (distended)
        assert rel(ch16.alpha_from_tube_law(0.3), 0.7 ** (-1 / 1.5)) < 1e-15
    with pytest.warns(RuntimeWarning, match="16.188"):
        ch16.collapsible_terms(1.3, 1e-5, dict(A0=A0, Kp=Kp, rho=rho, mu=1e-3))
    # arrays: NaN exactly at the α > 1 entries, the others untouched
    al = np.array([0.5, 1.0, 1.7, 0.9, 2.5])
    bad = al > 1
    with pytest.warns(RuntimeWarning, match="NaN"):
        P = ch16.collapsible_tube_law(al)
    assert np.array_equal(np.isnan(P), bad) and close(P[~bad], 1 - al[~bad] ** -1.5, rtol=1e-15)
    with pytest.warns(RuntimeWarning, match="NaN"):
        trio = ch16.collapsible_tube_law_derivatives(al)
    assert all(np.array_equal(np.isnan(v), bad) for v in trio)
    with pytest.warns(RuntimeWarning, match="NaN"):
        c = ch16.collapsible_wave_speed(al, Kp, rho)
    assert np.array_equal(np.isnan(c), bad) and close(c[~bad], np.sqrt(1.5 * Kp * al[~bad] ** -1.5 / rho), rtol=1e-15)
    with pytest.warns(RuntimeWarning, match="NaN"):
        a_back = ch16.alpha_from_tube_law(np.array([-1.0, 0.0, 0.3]))
    assert np.isnan(a_back[2]) and close(a_back[:2], [2 ** (-1 / 1.5), 1.0], rtol=1e-15)
    # silent on the fitted branch, at α = 1 exactly and within the stated round-off allowance (1e-12 relative)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        ch16.collapsible_tube_law(np.linspace(0.05, 1.0, 50))
        ch16.collapsible_tube_law(1.0 + 1e-13)
        ch16.collapsible_wave_speed(1.0, Kp, rho)
        ch16.collapsible_tube_law_derivatives(np.array([0.3, 1.0]))
        assert ch16.alpha_from_tube_law(0.0) == 1.0
        ch16.collapsible_terms(1.0, 1e-5, dict(A0=A0, Kp=Kp, rho=rho, mu=1e-3))
    assert ch16.ALPHA_MAX_COLLAPSE == 1.0


def test_collapsible_V2_throat_condition_alpha_t_le_0397():  # V2 + V1 — C15, Fig. 16.20(a): H ≤ 0 ⇔ upstream α ≤ 1
    """At a throat with S = 1 (u = c(α_t)), Bernoulli's constant is H = K_pP(α_t) + ½ρc_t² = K_p[1 − (1 − n/2)α_t⁻ⁿ];
    far upstream (u → 0) K_pP(α_far) = H, so α_far ≤ 1 ⇔ H ≤ 0 ⇔ α_t ≤ (1 − n/2)^{1/n} (= 0.25^{2/3} = 0.3969 at n = 3/2).
    The subcritical branch rises monotonically towards α_far, so H ≤ 0 keeps every upstream state on the fit."""
    al, n, K, r = sp.symbols("alpha_t n K_p rho", positive=True)
    P = 1 - al ** -n                                                                          # (16.188)
    c2 = al / r * K * sp.diff(P, al)                                                          # (16.198)
    H = K * P + r * c2 / 2
    assert sp.simplify(H - K * (1 - (1 - n / 2) * al ** -n)) == 0
    a_max = (1 - n / 2) ** (1 / n)
    assert sp.simplify(sp.powdenest(H.subs(al, a_max), force=True)) == 0                  # (valid for 0 < n < 2)
    assert sp.simplify(H.subs({n: sp.Rational(3, 2)}).subs(al, sp.Rational(1, 4) ** sp.Rational(2, 3))) == 0
    assert sp.simplify(sp.diff(H, al) - K * n * (1 - n / 2) * al ** (-n - 1)) == 0          # > 0 for n < 2: monotone
    a_max15 = 0.25 ** (2 / 3)
    assert abs(float(a_max.subs(n, sp.Rational(3, 2))) - a_max15) < 1e-15 and abs(a_max15 - 0.39685) < 1e-5
    # numerically, with the module (the demo's inputs): α_t = 0.35 → H < 0, all upstream α ≤ 1, silent
    I_ = ch16.illustrative_inputs()
    rho, Kp = I_["rho_water"], I_["Kp_collapse"]
    A0 = PI * I_["a_collapse"] ** 2
    Hs = lambda a_t: Kp * float(ch16.collapsible_tube_law(a_t)) + 0.5 * rho * (
        float(ch16.choking_flux(a_t, A0, Kp, rho)) / (a_t * A0)) ** 2
    for a_t in (0.2, 0.3, 0.35, 0.39, 0.40, 0.45, 0.6):
        assert rel(Hs(a_t), Kp * (1 - 0.25 * a_t ** -1.5)) < 1e-13
        assert (Hs(a_t) <= 0) == (a_t <= a_max15)
    Qs = float(ch16.choking_flux(0.35, A0, Kp, rho))
    H = Hs(0.35)
    assert abs(H + 41.47) < 0.01 and abs(Qs * 1e6 - 47.646) < 1e-3                         # the demo's numbers
    A0f = lambda x: A0 * (1 + (x / 0.05) ** 2)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        sub = np.array([ch16.collapsible_frictionless_alpha(Qs, H, A0f(x), Kp, rho, n=1.5, branch="sub")
                        for x in np.linspace(-0.2, -1e-3, 60)])                  # the demo shows x ≥ −12 cm
        a_far = float(ch16.alpha_from_tube_law(H / Kp))
    assert np.all(np.diff(sub) < 0) and np.all(sub < a_far) and a_far < 1                   # rises towards α_far ≤ 1
    a12 = ch16.collapsible_frictionless_alpha(Qs, H, A0f(-0.12), Kp, rho, n=1.5, branch="sub")
    assert abs(a12 - 0.8758) < 5e-5                                                          # the demo's α at x = −12 cm
    # α_t above the bound: H > 0 and the far-upstream root leaves the fit → the module warns
    Qb = float(ch16.choking_flux(0.45, A0, Kp, rho))
    with pytest.warns(RuntimeWarning, match="16.188"):
        a_bad = ch16.collapsible_frictionless_alpha(Qb, Hs(0.45), 15 * A0, Kp, rho, n=1.5, branch="sub")
    assert a_bad > 1


def test_collapsible_V1_frictionless_alpha_on_wide_upstream_sections():  # V1 — C15 (loop 3 finding: a* bracket)
    """Far upstream of a throat (A₀ ≫ A₀,throat) the subcritical root tends to α_far = (1 − H/K_p)^{−1/n} (u → 0).  The
    critical area ratio is α* = (ρQ²/(A₀²K_pn))^{1/(2−n)} (S = 1), = α_t(A₀,t/A₀)⁴ at n = 3/2; the module brackets α*
    in a fixed [1e-6, 1e6], so α* < 1e-6 (here A₀/A₀,t ≳ 24) must still give the root, not a ValueError."""
    I_ = ch16.illustrative_inputs()
    rho, Kp, n = I_["rho_water"], I_["Kp_collapse"], 1.5
    A0t = PI * I_["a_collapse"] ** 2
    Qs = float(ch16.choking_flux(0.35, A0t, Kp, rho))
    H = Kp * float(ch16.collapsible_tube_law(0.35)) + 0.5 * rho * (Qs / (0.35 * A0t)) ** 2
    a_far = (1 - H / Kp) ** (-1 / n)
    got, failures = {}, []
    for ratio in (5.0, 20.0, 25.0, 30.0, 100.0):
        a_star = (rho * Qs ** 2 / ((ratio * A0t) ** 2 * Kp * n)) ** (1 / (2 - n))
        assert rel(a_star, 0.35 * ratio ** -4) < 1e-12
        try:
            a = ch16.collapsible_frictionless_alpha(Qs, H, ratio * A0t, Kp, rho, n=n, branch="sub")
        except ValueError as e:
            failures.append(f"A0/A0t = {ratio:g} (alpha* = {a_star:.2e}): ValueError {e}")
            continue
        res = Kp * (1 - a ** -n) + 0.5 * rho * (Qs / (a * ratio * A0t)) ** 2 - H          # Bernoulli residual
        assert abs(res) < 1e-9 * Kp and a < a_far
        got[ratio] = a
    assert not failures, "collapsible_frictionless_alpha raises on valid wide sections: " + "; ".join(failures)
    assert got[5.0] < got[20.0] < got[25.0] < got[30.0] < got[100.0] < a_far and a_far - got[100.0] < 1e-4


def test_collapsible_V1_starling_choke_closed_form():  # V1 — C15 (implementer point 1: Q_choke, p₂ at the choke)
    ref = _starling_closed(**VEIN)
    prm = dict(VEIN)
    fl = ch16.flow_limitation_curve(np.array([-2500.0, -1200.0, -400.0, -100.0]), p1=0.0, pe=0.0, prm=prm)
    assert rel(fl["Q_choke"], ref["Q"]) < 1e-11 and rel(fl["p2_choke"], ref["p2"]) < 1e-11
    assert abs(ref["Q"] - 1.99899e-5) < 1e-10 and abs(ref["p2"] + 946.13) < 0.01         # the design's numbers
    assert np.all(fl["limited"] == [True, True, False, False])
    assert fl["Q"][0] == fl["Q"][1] == fl["Q_choke"] and np.all(np.diff(fl["Q"][1:]) < 0)    # plateau, then falls with p₂
    G = ref["G"]
    for frac in (0.3, 0.7, 0.95):                                                       # subcritical outlet states
        Q = frac * ref["Q"]
        a_out = brentq(lambda a: G(1.0, Q) - G(a, Q) - VEIN["L"], ref["a_star"](Q), 1.0, xtol=1e-15)
        p2 = VEIN["Kp"] * (1 - a_out ** -1.5)
        st = ch16.starling_state(0.0, p2, 0.0, prm)
        assert rel(st["Q"], Q) < 1e-9 and not st["limited"] and rel(st["alpha_out"], a_out) < 1e-9
    assert ch16.starling_state(0.0, 10.0, 0.0, prm)["Q"] == 0.0


def test_collapsible_V1_march_against_closed_form_and_S2_equation():  # V1 — C15 (collapsible_march, (16.205) along the march; no order asserted)
    ref = _starling_closed(**VEIN)
    Q = 0.8 * ref["Q"]
    tube = ch16.CollapsibleTube(A0=VEIN["A0"], Kp=VEIN["Kp"], rho=VEIN["rho"], mu=VEIN["mu"])
    m = ch16.collapsible_march(Q, 1.0, (0.0, VEIN["L"]), tube, n_out=41)
    G = ref["G"]
    x_cf = np.array([G(1.0, Q) - G(a, Q) for a in m["alpha"]])
    assert close(x_cf, m["x"], rtol=0, atol=1e-8)                                        # x(α) of Bernoulli with friction
    assert not m["choked"] and np.all(np.diff(m["alpha"]) < 0) and np.all(np.diff(m["S"]) > 0)   # S < 1: area falls, S rises
    assert close(m["p"], VEIN["Kp"] * (1 - m["alpha"] ** -1.5), rtol=1e-12, atol=1e-9)
    lnS2 = np.log(m["S"][1:-1] ** 2)
    dlnS2 = np.gradient(np.log(m["S"] ** 2), m["x"])
    pred = np.array([ch16.collapsible_speed_index_rhs(xx, aa, Q, tube) for xx, aa in zip(m["x"], m["alpha"])])
    assert np.max(np.abs(dlnS2[2:-2] - pred[2:-2])) < 2e-3 * np.max(np.abs(pred))       # (16.205) along the march
    mc = ch16.collapsible_march(1.02 * ref["Q"], 1.0, (0.0, VEIN["L"]), tube)
    assert mc["choked"] and mc["x_choke"] < VEIN["L"] and abs(mc["S"][-1] - 1) < 1e-3     # chokes inside the tube
    tr = ch16.collapsible_terms(0.6, Q, tube)
    assert tr["friction"] < 0 and tr["dalpha_dx"] < 0
    sup = ch16.collapsible_terms(0.1, Q, tube)                                            # supercritical: the opposite
    assert sup["S2"] > 1 and sup["dalpha_dx"] > 0
    a_s = (VEIN["rho"] * Q ** 2 / (VEIN["A0"] ** 2 * VEIN["Kp"] * 1.5)) ** 2              # S = 1 exactly
    t1 = ch16.collapsible_terms(a_s, Q, tube)
    assert abs(t1["one_minus_S2"]) < 1e-12 and np.isfinite(t1["driver"])                  # T10: driver formed first
    t2 = ch16.collapsible_terms(0.6, Q, tube, dKpdx=50.0)                                 # P < 0 and stiffening
    assert t2["P"] < 0 and t2["stiffness"] > 0                                            # S24 / T9
    assert rel(t1["driver"], t1["area"] + t1["external"] + t1["gravity"] + t1["friction"] + t1["stiffness"]) < 1e-14
    assert rel(ch16.collapsible_rhs(0.0, [0.6], tube, Q)[0], tr["dalpha_dx"]) < 1e-15


@pytest.mark.filterwarnings("error::RuntimeWarning")
def test_collapsible_V1_frictionless_throat_chokes_at_the_throat():  # V1 — Fig. 16.20(a) (D30 at S = 1)
    Kp, rho, A0t = 200.0, 1000.0, 1.0e-4
    A0f = lambda x: A0t * (1 + 0.6 * (x - 0.05) ** 2 / 0.05 ** 2)                       # throat at x = 0.05
    a_t = 0.35                                    # ≤ 0.25^{2/3} = 0.397: H ≤ 0, every upstream state on α ≤ 1
    Q = float(ch16.choking_flux(a_t, A0t, Kp, rho))                                      # S = 1 at the throat
    H = Kp * float(ch16.collapsible_tube_law(a_t)) + 0.5 * rho * Q ** 2 / (a_t * A0t) ** 2
    assert H < 0
    for xx in (0.0, 0.02, 0.04):
        a_sub = ch16.collapsible_frictionless_alpha(Q, H, A0f(xx), Kp, rho, branch="sub")
        S2 = ch16.collapsible_terms(a_sub, Q, dict(A0=A0f(xx), Kp=Kp, rho=rho, mu=1e-3, friction=False))["S2"]
        assert S2 < 1                                                                      # subcritical upstream
    assert rel(ch16.collapsible_frictionless_alpha(Q, H, A0t, Kp, rho, branch="sub"), a_t) < 1e-6
    tube = ch16.CollapsibleTube(A0=A0f, Kp=Kp, rho=rho, mu=1e-3, friction=False)
    a_in = ch16.collapsible_frictionless_alpha(Q, H, A0f(0.0), Kp, rho, branch="sub")
    m = ch16.collapsible_march(Q, a_in, (0.0, 0.049), tube, n_out=11)
    ref = [ch16.collapsible_frictionless_alpha(Q, H, A0f(xx), Kp, rho, branch="sub") for xx in m["x"]]
    assert close(m["alpha"], ref, rtol=1e-6)                                              # (16.204) ≡ Bernoulli here
    assert abs(ch16.collapsible_terms(a_t, Q, tube, x=0.05)["driver"]) < 1e-9             # dA₀/dx = 0 where S = 1


def test_collapsible_V2_influence_coefficients():  # V2 — C15 (the ch15-style table)
    S2, al, n = 0.37, 0.62, 1.5
    ic = ch16.collapsible_influence_coefficients(S2, al, n)
    P, dP, d2P = (float(v) for v in ch16.collapsible_tube_law_derivatives(al, n))
    M = 3 + al * d2P / dP
    k = 1 / (1 - S2)
    assert rel(ic["dA0/A0"][0], S2 * k) < 1e-14 and rel(ic["dA0/A0"][1], (-2 + (2 - M) * S2) * k) < 1e-14
    assert rel(ic["D"][0], -k) < 1e-14 and rel(ic["D"][1], M * k) < 1e-14
    q = P / (al * dP)
    assert rel(ic["dKp/Kp"][0], -q * k) < 1e-14 and rel(ic["dKp/Kp"][1], (M * q - (1 - S2)) * k) < 1e-14
    assert ch16.collapsible_influence_matrix is ch16.collapsible_influence_coefficients


# ======================================================================================================================
# C16 — Casson pipe flow (16.207)–(16.220); D31–D33
# ======================================================================================================================
CAS = dict(a=0.31e-3, dp=19.0, L=0.021, tau_y=4.7e-3, Kc=0.057)


def test_casson_pipe_V2_derivation():  # V2 — D31, D32, D33 (★★★)
    s, a, c, G, K = sp.symbols("s a c G K", positive=True)            # r = s², a = A², a_c = C²; G = Δp/L
    A, C = sp.symbols("A C", positive=True)
    r = s ** 2
    tau, tauy = G * r / 2, G * C ** 2 / 2                              # D31: (16.209) τ = Δp r/(2L), plug edge (16.210)
    dudr = -(sp.sqrt(tau) - sp.sqrt(tauy)) ** 2 / K ** 2              # (16.211)
    page_212 = -(1 / (2 * K ** 2)) * G * (r + C ** 2 - 2 * sp.sqrt(r * C ** 2))
    assert sp.simplify(dudr - page_212) == 0                                                     # (16.212)
    u215 = 1 / (4 * K ** 2) * G * ((A ** 4 - r ** 2) - sp.Rational(8, 3) * C * (A ** 3 - s ** 3) + 2 * C ** 2 * (A ** 2 - r))   # (16.215)
    assert sp.simplify(sp.diff(u215, s) / sp.diff(r, s) - dudr) == 0 and sp.simplify(u215.subs(s, A)) == 0
    tw = G * A ** 2 / 2
    u216 = A ** 2 * tw / (2 * K ** 2) * ((1 - (r / A ** 2) ** 2) - sp.Rational(8, 3) * sp.sqrt(tauy / tw) * (1 - (s / A) ** 3)
                                         + 2 * tauy / tw * (1 - r / A ** 2))                     # (16.216)
    assert sp.simplify(u216 - u215) == 0
    uc = 1 / (4 * K ** 2) * G * (A - C) ** 3 * (A + C / 3)                                        # (16.218)
    assert sp.simplify(u215.subs(s, C) - uc) == 0 and sp.simplify(dudr.subs(s, C)) == 0           # plug edge, γ̇ = 0
    Q = sp.pi * C ** 4 * uc + sp.integrate(2 * sp.pi * r * u215 * sp.diff(r, s), (s, C, A))       # (16.219)
    xi = C ** 2 / A ** 2
    Q220 = sp.pi / 8 / K ** 2 * G * A ** 8 * (1 - sp.Rational(16, 7) * sp.sqrt(xi) + sp.Rational(4, 3) * xi - xi ** 4 / 21)
    assert sp.simplify(sp.expand(Q - Q220)) == 0                                                  # (16.220)
    assert sp.simplify(Q220.subs(C, 0) - sp.pi * A ** 8 * G / (8 * K ** 2)) == 0                  # ξ → 0: Poiseuille, μ = K_c²
    assert sp.simplify(Q220.subs(C, A)) == 0                                                      # ξ → 1: no flow
    a_, dp_, L_, ty_, Kc_ = CAS.values()
    ac_ = 2 * L_ * ty_ / dp_
    num = {A: math.sqrt(a_), C: math.sqrt(ac_), G: dp_ / L_, K: Kc_}
    for rr in (0.2 * a_, 0.6 * a_, 0.97 * a_):
        assert rel(ch16.casson_pipe_velocity(rr, *CAS.values()), float(u215.subs(num).subs(s, math.sqrt(rr)))) < 1e-12
    assert rel(ch16.casson_core_velocity(*CAS.values()), float(uc.subs(num))) < 1e-12
    assert rel(ch16.casson_pipe_flow_rate(*CAS.values()), float(Q220.subs(num))) < 1e-12


def test_casson_pipe_V1_quadrature_and_plug():  # V1 — C16
    a, dp, L, ty, Kc = CAS.values()
    ac = ch16.casson_plug_radius(dp, L, ty)
    assert rel(ac, 2 * L * ty / dp) < 1e-15 and 0 < ac < a
    assert rel(ac / a, ty / ch16.casson_wall_stress(dp, L, a)) < 1e-14                            # (16.210)
    uc = ch16.casson_core_velocity(a, dp, L, ty, Kc)
    r = np.linspace(0, a, 401)
    u = ch16.casson_pipe_velocity(r, a, dp, L, ty, Kc)
    assert close(u[r <= ac], uc, rtol=1e-13) and abs(u[-1]) < 1e-15 * uc
    Qq, _ = quad(lambda x: 2 * PI * x * float(ch16.casson_pipe_velocity(x, a, dp, L, ty, Kc)), 0, a, points=[ac],
                 epsabs=0, epsrel=1e-13)
    Q = ch16.casson_pipe_flow_rate(a, dp, L, ty, Kc)
    assert rel(Q, Qq) < 1e-11
    xi = ac / a
    QP = PI * a ** 4 * dp / (8 * Kc ** 2 * L)
    assert rel(Q / QP, ch16.casson_flux_ratio(xi)) < 1e-13 and rel(ch16.casson_poiseuille_flux_ratio(xi), ch16.casson_flux_ratio(xi)) < 1e-15
    rr = np.linspace(ac * 1.05, a * 0.999, 50)                  # the law holds in the annulus (FD step ≪ r − a_c)
    h = a * 1e-7
    gam = -(ch16.casson_pipe_velocity(rr + h, a, dp, L, ty, Kc) - ch16.casson_pipe_velocity(rr - h, a, dp, L, ty, Kc)) / (2 * h)
    assert close(ch16.casson_stress(gam, ty, Kc), dp * rr / (2 * L), rtol=1e-6)


def test_casson_pipe_V4_V7_force_balance_limits_and_slip_S26():  # V4 + V7 — C16, N105, N106 (trap T21)
    a, dp, L, ty, Kc = CAS.values()
    assert rel(ch16.casson_wall_stress(dp, L, a) * 2 * PI * a * L, PI * a ** 2 * dp) < 1e-14       # force balance
    dmin = ch16.casson_min_pressure_drop(L, a, ty)
    assert rel(dmin, 2 * L * ty / a) < 1e-15 and rel(ch16.casson_plug_radius(dmin, L, ty), a) < 1e-14
    assert ch16.casson_pipe_flow_rate(a, 0.999 * dmin, L, ty, Kc) == 0.0                          # a_c ≥ a: no flow
    assert np.all(ch16.casson_pipe_velocity(np.linspace(0, a, 5), a, 0.9 * dmin, L, ty, Kc) == 0)
    assert rel(ch16.casson_pipe_flow_rate(a, dp, L, 1e-30, Kc), PI * a ** 4 * dp / (8 * Kc ** 2 * L)) < 1e-10
    xs = np.array([0.9, 0.95, 0.975, 0.9875])
    q = ch16.casson_flux_ratio(xs)
    assert np.all(np.diff(q) < 0) and abs(observed_order(1 - xs, q) - 3.0) < 0.15                  # Q ∝ (1 − ξ)³ (u_c ∝ (√a − √a_c)³)
    # S26: the minimum pressure drop for L/a = 200 and τ_y = 0.05 dyn/cm² (our numbers) is 20 dyn/cm², not 0.05
    assert rel(ch16.casson_min_pressure_drop(200.0, 1.0, ch16.dyn_cm2_to_Pa(0.05)), ch16.dyn_cm2_to_Pa(20.0)) < 1e-14


def test_power_law_and_herschel_bulkley_V1():  # V1 — S01 pointer (optional functions)
    a, dp, L = 0.4e-3, 31.0, 0.017
    m = 3.7e-3
    assert rel(ch16.power_law_pipe_flow_rate(a, dp, L, m, 1.0), PI * a ** 4 * dp / (8 * m * L)) < 1e-13
    for n in (0.6, 1.4):
        Qq, _ = quad(lambda x: 2 * PI * x * float(ch16.power_law_pipe_velocity(x, a, dp, L, m, n)), 0, a, epsabs=0, epsrel=1e-12)
        assert rel(ch16.power_law_pipe_flow_rate(a, dp, L, m, n), Qq) < 1e-10
        assert rel(ch16.herschel_bulkley_pipe_flow_rate(a, dp, L, m, n, 1e-30), ch16.power_law_pipe_flow_rate(a, dp, L, m, n)) < 1e-9
    tw = a * dp / (2 * L)
    phi = 0.23
    bingham = PI * a ** 4 * dp / (8 * m * L) * (1 - 4 / 3 * phi + phi ** 4 / 3)                   # Buckingham–Reiner
    assert rel(ch16.herschel_bulkley_pipe_flow_rate(a, dp, L, m, 1.0, phi * tw), bingham) < 1e-12
    assert ch16.herschel_bulkley_pipe_flow_rate(a, dp, L, m, 1.0, 1.01 * tw) == 0.0


# ======================================================================================================================
# C17 — water potential (16.221); xylem and phloem (16.222); D34
# ======================================================================================================================
def test_water_potential_V1_terms_and_column():  # V1 — C17, D34
    p, c, T, z = -0.83e6, 41.0, 296.0, 17.0
    t = ch16.water_potential_terms(p, c, T, z)
    R = 8.31446261815324
    assert rel(t["osmotic"], -R * T * c) < 1e-15 and rel(t["gravity"], 1000 * 9.80665 * z) < 1e-15
    assert rel(t["total"], p - R * T * c + 1000 * 9.80665 * z) < 1e-15
    assert rel(ch16.water_potential(p, c, T, z), t["total"]) < 1e-15
    assert rel(ch16.osmotic_pressure(c, T), R * T * c) < 1e-15
    zz = np.linspace(0, 31, 32)
    col = ch16.xylem_column(zz, psi_soil=-0.27e6, u=1.3e-3, a=21e-6, mu=1.0e-3, c=11.0, T=T)
    grad = 8 * 1e-3 * 1.3e-3 / 21e-6 ** 2
    assert rel(col["grad_viscous"], grad) < 1e-14
    assert close(col["psi"], -0.27e6 - grad * zz, rtol=1e-12)                                     # ψ falls by the viscous drop
    assert close(col["p"], -0.27e6 + R * T * 11.0 - 1000 * 9.80665 * zz - grad * zz, rtol=1e-12)
    assert col["psi"][-1] < 0 and col["p"][-1] < 0                                                 # T22: tension
    still = ch16.water_potential(-1000 * 9.80665 * zz, 0.0, T, zz)                                  # pure water at rest
    assert np.max(np.abs(still)) < 1e-9                                                            # ψ uniform: no flow
    assert rel(ch16.cavitation_margin(-1.1e6, p_cav=-2.6e6), 1.5e6) < 1e-15
    assert rel(ch16.xylem_reynolds(1.3e-3, 21e-6, 1.0e-6), 1.3e-3 * 42e-6 / 1e-6) < 1e-15


def test_water_potential_V2_units_and_book_units_form():  # V2 — D34 (slip S30), N110
    for q in (Q_(-0.8, "MPa"), Q_(8.314, "J/(mol*K)") * Q_(300, "K") * Q_(30, "mol/m**3"),
              Q_(1000, "kg/m**3") * Q_(9.8, "m/s**2") * Q_(20, "m")):
        assert q.check("[pressure]")
    # (16.221) in bar with R = 83.141 cm³ bar/(mol K) and ρgz/10⁶ equals the SI form ÷ 1e5 (g = 980 cm/s²)
    p, c, T, z = -8.3, 4.1e-5, 296.0, 1700.0                                   # bar, mol/cm³, K, cm
    b = ch16.water_potential_book_units(p, c, T, z)
    si = ch16.water_potential(p * 1e5, c * 1e6, T, z / 100, g=9.80) / 1e5
    assert abs(b - si) < 2e-4 * abs(si)                                         # R 83.141 vs CODATA 83.1446 (4e-5 rel)
    assert rel(ch16.R_BOOK_CM3_BAR, ch16.R_MOLAR * 10) < 5e-5
    without = p - ch16.R_BOOK_CM3_BAR * T * c + 1.0 * 980 * z                 # the printed form, ρgz not divided
    assert abs(without - b) > 1e5 * abs(b - si)


def test_sieve_tube_V1_limits():  # V1 — N113 (16.222)
    Q, mu, a, L = 3e-14, 1.5e-3, 10e-6, 300e-6
    d = ch16.sieve_tube_pressure_drop(Q, mu, a, a, L, 0.0, 1)                  # ℓ = 0, N = 1, r = a: Poiseuille only
    assert rel(d["total"], 8 * mu * Q * L / (PI * a ** 4)) < 1e-14 and d["plate_entry"] == 0.0
    N, r = 9, a / 3                                                            # a_e = a/√N = r: Δp′ = 0
    d = ch16.sieve_tube_pressure_drop(Q, mu, a, r, L, 1e-6, N)
    assert abs(d["plate_entry"]) < 1e-15 * d["total"]
    assert rel(d["plate_viscous"], 8 * mu * (Q / N) * 1e-6 / (PI * r ** 4)) < 1e-13             # N pores in parallel
    assert rel(d["total"], d["tube"] + d["plate_viscous"] + d["plate_entry"]) < 1e-15
    d = ch16.sieve_tube_pressure_drop(Q, mu, a, 1.5e-6, L, 1e-6, 20)
    ae = a / math.sqrt(20)
    dpp = 8 * mu * Q / (PI * a ** 3) * (ae / 1.5e-6) * (0.57 * 20 * ((ae / 1.5e-6) ** 3 - 1) - 1.5 * (1 - 1.5e-6 / ae))
    assert rel(d["plate_entry"], 2 * dpp) < 1e-13 and 0 < d["plate_share"] < 1
    with pytest.raises(ValueError):
        ch16.sieve_tube_pressure_drop(Q, mu, a, 4e-6, L, 1e-6, 20)


# ======================================================================================================================
# Tables, slips, contract coverage, book values, scripts, references
# ======================================================================================================================
def test_tables_V7_structure_and_slips_with_printed_options():  # V7 — trap T2, T24, slips table
    sl = ch16.book_slips()
    assert len(sl) == 34 and list(sl["id"]) == [f"S{i}" for i in range(1, 35)]
    assert set(sl["kind"]) == {"slip", "loose"} and (sl["kind"] == "slip").sum() == 12
    s34 = sl[sl["id"] == "S34"].iloc[0]                                       # found in verification (loop 1)
    assert s34["kind"] == "slip" and bool(s34["changes_physics"]) and "(16.174)" in s34["where"]
    assert "printed_16_174" in s34["code"] and "DEAN_C4_DERIVED" in s34["code"]
    for opt in ("printed_16_30", "printed_16_54", "printed_16_96", "printed_16_140", "printed_16_201", "printed_16_174"):
        assert any(opt in c for c in sl["code"]), opt
        assert opt in (ch16.__doc__ or ""), opt                               # listed in the module docstring too
    assert "34 rows" in ch16.__doc__ and "34 suspected slips" in ch16.book_slips.__doc__
    assert len(ch16.traps()) == 27 and len(ch16.conventions_table()) >= 15 and len(ch16.units_table()) == 7
    inp = ch16.illustrative_inputs()
    assert all(isinstance(v, float) for v in inp.values())
    if BOOK.exists():
        forbidden = book()["_forbidden_public"]
        src = (Path(__file__)).read_text(encoding="utf-8")
        assert not [f for f in forbidden if f in src]                         # this public file quotes no book number


def test_contract_V7_every_design_part_c_function_is_exercised_here():  # V7 — design Part C (C.1–C.10)
    lines = (ROOT / "analysis" / "ch16_design.md").read_text(encoding="utf-8").splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith("### C.1 "))
    stop = next(i for i, l in enumerate(lines) if l.startswith("### C.11 "))
    names = []
    for l in lines[start:stop]:
        if l.startswith("| `"):
            names += re.findall(r"`([A-Za-z_][A-Za-z0-9_]*)`", l.split("|")[1])
    assert len(names) >= 160
    src = Path(__file__).read_text(encoding="utf-8")
    missing = [n for n in names if not hasattr(ch16, n)]
    untested = [n for n in names if not re.search(r"ch16\." + re.escape(n) + r"\b", src)]
    assert not missing and not untested, (missing, untested)


def test_docstrings_V7_no_stale_validation_text():  # V7 — ch01 lesson: stale phrases
    import inspect
    from fluidpy.core import rheology, womersley
    for mod in (ch16, womersley, rheology):
        text = inspect.getsource(mod)
        for phrase in ("not yet pinned", "benchmark pending", "TODO", "FIXME", "to be verified"):
            assert phrase not in text, (mod.__name__, phrase)


@book_only
def test_book_V6_worked_numbers():  # V6 — private book values (report states relative errors only)
    b = book()
    c = b["sec_16_2_circulation"]["worked_CO"]
    assert rel(ch16.cardiac_output(ch16.mL_to_m3(c["SV_mL"]), c["HR_per_min"]) * 60e6, c["CO_mL_per_min_printed"]) < 1e-12
    w = b["sec_16_2_vessels"]["worked_PRU"]
    R = ch16.peripheral_resistance(ch16.mmHg_to_Pa(w["pA_mean_mmHg"]), w["Q_cm3_s"] * 1e-6) / ch16.PRU
    assert round(R, 1) == w["R_PRU"] and rel(R, b["sec_16_2_vessels"]["derived_by_analyst"]["R_PRU_exact"]) < 1e-4
    wa = b["sec_16_3_womersley_and_steady"]["worked_alpha"]
    al = ch16.womersley_number(wa["a_cm"] / 100, 2 * PI * wa["HR_per_min"] / 60, ch16.poise_to_Pa_s(wa["mu_poise"]) / (wa["rho_g_cm3"] * 1000))
    assert round(al) == wa["alpha"] and abs(al / wa["alpha"] - 1) < 0.006
    X, Y = ch16.womersley_XY(b["sec_16_3_womersley_rigid"]["example_alpha"])
    dv = b["sec_16_3_womersley_rigid"]["derived_by_analyst"]
    assert abs(X - dv["X"]) < 5e-5 and abs(Y - dv["Y"]) < 5e-5
    assert rel(Y, b["sec_16_3_womersley_rigid"]["Y"]) < 0.006                       # within the printed rounding
    assert 0.025 < 1 - X / b["sec_16_3_womersley_rigid"]["X"] < 0.04                 # slip S11: ≈ 3 % below the print
    assert abs(ch16.womersley_flow_ratio(b["sec_16_3_womersley_rigid"]["example_alpha"]) - 0.1) < 0.006
    cs = b["sec_16_3_casson"]["worked_min_dp"]
    dmin = ch16.casson_min_pressure_drop(cs["L_over_a"], 1.0, ch16.dyn_cm2_to_Pa(cs["tau_y_dyn_cm2"]))
    assert rel(dmin, b["sec_16_3_casson"]["derived_by_analyst"]["correct_Pa"]) < 1e-12
    assert abs(ch16.Pa_to_mmHg(dmin) - cs["printed_mmHg"]) < 0.005                    # the printed mmHg is right
    assert rel(dmin / ch16.dyn_cm2_to_Pa(cs["printed_dyn_cm2"]), 1000) < 1e-12        # S26: the dyn/cm² figure is 1000× small
    xy = b["sec_16_4_plants"]["xylem"]
    Re_pr = ch16.xylem_reynolds(xy["speed_cm_s"] / 100, xy["radius_mm"] / 1000, xy["nu_cm2_s_printed"] * 1e-4)
    Re_w = ch16.xylem_reynolds(xy["speed_cm_s"] / 100, xy["radius_mm"] / 1000, 1e-6)
    assert rel(Re_w, xy["Re_printed"]) < 1e-12 and rel(Re_pr, xy["Re_printed"] / 10) < 1e-12   # S28
    assert rel(ch16.R_MOLAR * 10, b["sec_16_4_plants"]["R_cm3_bar_per_mol_K"]) < 5e-5
    assert rel(ch16.MMHG, b["sec_16_2_circulation"]["mmHg_to_Pa_printed"]) < 0.003
    ws, bl = b["sec_16_3_womersley_and_steady"], b["sec_16_2_blood"]
    Re = ch16.reynolds_number(ws["aorta_mean_speed_cm_s"] / 100, ws["aorta_root_diameter_mm"] / 1000,
                              rho=bl["density_kg_m3"], mu=ch16.cP_to_Pa_s(bl["asymptotic_viscosity_cP"]))
    assert rel(Re, ws["aorta_Re"]) < 0.005
    dean = b["sec_16_3_curved_tube"]
    assert rel(1541 / 50400, dean["coefficients_178"][0]) < 0.001
    assert rel(ch16.dean_convert(dean["valid_D_max"], from_="D", to="K"), dean["derived_by_analyst"]["K_from_D_at_100"]) < 1e-12
    assert abs(dean["coefficients_182"][1] - (dean["coefficients_178"][0] ** 2 - dean["coefficients_178"][1])) < 1e-4   # (16.182) = 1/(16.178)
    rng = b["sec_16_2_blood"]["yield_stress_range_dyn_cm2_H45"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        lo = ch16.yield_stress_correlation(0.45, 0.21) / 0.1
        hi = ch16.yield_stress_correlation(0.45, 0.46) / 0.1
    assert lo > rng[0] and hi > rng[1]                                               # S29: the correlation's range is higher


@slow
def test_scripts_V7_every_chapter_script_runs_headless(tmp_path):  # V7 — the demos exit 0; make_refs idempotent
    ran = 0
    for s in sorted((ROOT / "scripts").glob("ch16_*.py")):
        if s.name == "ch16_common.py":
            continue
        r = subprocess.run([sys.executable, str(s), "--no-show", "--fast", "--out", str(tmp_path)], cwd=str(ROOT),
                           capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=900)
        assert r.returncode == 0, (s.name, r.stderr[-2000:])
        ran += 1
    assert ran == 13 and len(list(tmp_path.glob("*.png"))) >= 13
    before = (REF / "benchmarks.json").read_bytes()
    r = subprocess.run([sys.executable, str(REF / "make_refs.py")], cwd=str(ROOT), capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=120)
    assert r.returncode == 0 and (REF / "benchmarks.json").read_bytes() == before


def test_reference_V5_sources_file_cites_every_entry():  # V5 — data-and-benchmarks policy
    src = (REF / "SOURCES.md").read_text(encoding="utf-8")
    bm = bench()
    keys = [k for k in bm if not k.startswith("_")]
    assert len(keys) == 7
    for key in keys:
        assert f"`{key}`" in src and bm[key]["url"] in src and bm[key]["verified"] in src
