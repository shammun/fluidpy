"""Verification suite for Chapter 1 — Introduction (Kundu, Cohen & Dowling 5e, §§1.1–1.11).

Evidence levels (``verify-implementation`` skill): V1 analytic · V2 symbolic / dimensional · V3 convergence ·
V4 conservation / invariant · V5 published benchmark (``reference/ch01/``, see SOURCES.md) · V6 book value (private,
``tests/book_values_ch01.json``, skipped when absent) · V7 limits / symmetry / monotonicity. Every test name is
``test_<concept>_<level>_<what>``; the comment on the ``def`` line repeats the level.

A items (CORE, ≥ 2 independent levels): C06 C12 C20 C25 C35 C36 C40 C45 C50 C51 C54 C55 C64 C67 C69.
Derivations re-derived with sympy: D14, D18, D19 (★★★), D28 (★★★), D34, D36 (all ★★/★★★) plus D10, D20.

Run: ``.venv/Scripts/python.exe -m pytest tests/test_ch01.py -q -p no:cacheprovider``.
"""
from __future__ import annotations

import csv
import itertools
import json
import math
import re
import warnings
from fractions import Fraction
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import numpy as np  # noqa: E402
import pytest  # noqa: E402
import sympy as sp  # noqa: E402

from fluidpy import ch01_introduction as ch01  # noqa: E402
from fluidpy.core.units import Q_, ureg  # noqa: E402
from tools.convergence import observed_order, pairwise_orders  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference" / "ch01"
BOOK = Path(__file__).resolve().parent / "book_values_ch01.json"
needs_ref = pytest.mark.skipif(not (REF / "benchmarks.json").exists(), reason="run reference/ch01/make_refs.py")
book_only = pytest.mark.skipif(not BOOK.exists(), reason="book values are private; see CLAUDE.md rule 9")


def rel(a, b):
    """Maximum relative error |a − b| / |b| (b must be nonzero)."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return float(np.max(np.abs(a - b) / np.abs(b)))


def book():
    return json.loads(BOOK.read_text(encoding="utf-8"))


def ref_json(name):
    return json.loads((REF / name).read_text(encoding="utf-8"))


def ref_csv(name):
    with open(REF / name, newline="", encoding="utf-8") as f:
        return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]


def geopotential_m(Z_km, r0_km):
    """USSA-1976 geopotential altitude H = r0 Z/(r0 + Z) in metres (NASA-TM-X-74335, r0 from Table 2 + errata)."""
    return 1e3 * r0_km * Z_km / (r0_km + Z_km)


G = ch01.G0
M_AIR = ch01.molecular_mass(ch01.M_W_AIR)
N_AIR = ch01.number_density(1.225, ch01.M_W_AIR)


# ====================================================================================================================
# C06 · Continuum hypothesis (§1.4) — A item; with C04, C05 (D34), C07, N03, C24, D35
# ====================================================================================================================
def test_continuum_V1_relative_noise_follows_poisson_law():  # V1 (Poisson counting, D35)
    L = np.logspace(-8, -6, 7)  # 25 … 2.5e7 molecules per box
    mean, std = ch01.sample_density(L, N_AIR, M_AIR, n_samples=4000, seed=1)
    measured = std / mean
    expected = ch01.density_noise_expected(L, N_AIR)
    se = 1.0 / np.sqrt(2 * (4000 - 1))  # relative standard error of a sample std
    assert np.all(np.abs(measured / expected - 1) < 5 * se + 0.01), measured / expected
    slope = np.polyfit(np.log(L), np.log(measured), 1)[0]
    assert abs(slope + 1.5) < 0.05, slope  # noise ∝ L^(−3/2)


def test_continuum_V2_noise_law_symbolic_D35():  # V2 (D35 re-derived)
    n, L, m = sp.symbols("n L m", positive=True)
    Nbar = n * L**3  # step: expected count
    sigma_N = sp.sqrt(Nbar)  # Poisson: Var(N) = N̄
    rel_noise = sp.simplify((m * sigma_N / L**3) / (m * Nbar / L**3))  # steps 3-4
    assert sp.simplify(rel_noise - (n * L**3) ** sp.Rational(-1, 2)) == 0
    assert sp.simplify(sp.diff(sp.log(rel_noise), L) * L) == sp.Rational(-3, 2)  # log–log slope −3/2 (step 5)
    # worked number: 10 µm cube of sea-level air
    Nb = N_AIR * (1e-5) ** 3
    assert 2.4e10 < Nb < 2.6e10 and abs(ch01.density_noise_expected(1e-5, N_AIR) - Nb**-0.5) < 1e-18


def test_continuum_V1_from_scratch_uniform_positions_agree_with_sampler():  # V1 (independent mechanism: binomial)
    rng = np.random.default_rng(3)
    n_tot, s, reps = 20000, 0.1, 600  # 20000 molecules in a unit cube, sub-cube of side 0.1 → mean 20 molecules
    counts = np.array([np.sum(np.all(rng.random((n_tot, 3)) < s, axis=1)) for _ in range(reps)])
    rel_scratch = counts.std(ddof=1) / counts.mean()
    mean, std = ch01.sample_density(s, float(n_tot), 1.0, n_samples=reps, seed=4)  # n = n_tot per unit volume
    se = 1.0 / np.sqrt(2 * (reps - 1))
    assert abs(rel_scratch / (std / mean) - 1) < 5 * se * 1.5, (rel_scratch, std / mean)
    assert abs(counts.mean() / (n_tot * s**3) - 1) < 0.02


def test_continuum_V1_box_average_matches_quadrature_and_limits():  # V1 + V7 (plateau → drift)
    rho0, eps, Lf = 1.2, 0.1, 1.0
    x0 = 0.25
    for L in (1e-4, 0.1, 0.37, 1.0, 2.5):
        x = np.linspace(x0 - L / 2, x0 + L / 2, 20001)
        quad = np.trapezoid(rho0 * (1 + eps * np.sin(2 * np.pi * x / Lf)), x) / L
        assert rel(ch01.box_average_density(L, rho0, eps, Lf), quad) < 1e-8
    assert rel(ch01.box_average_density(1e-9, rho0, eps, Lf), rho0 * (1 + eps)) < 1e-12  # L → 0: point value
    assert rel(ch01.box_average_density(Lf, rho0, eps, Lf), rho0) < 1e-12  # whole wavelength: mean
    assert ch01.density_expected(0.5, N_AIR, M_AIR) == pytest.approx(N_AIR * M_AIR, rel=1e-14)
    assert ch01.density_expected(0.5, N_AIR, M_AIR, gradient=0.2) == pytest.approx(N_AIR * M_AIR * 1.05, rel=1e-14)


def test_continuum_V1_sample_mean_tracks_box_average():  # V1 (statistical: sampler mean = deterministic box average)
    L = np.array([0.05, 0.3, 0.8])
    mean, std = ch01.sample_density(L, 1e6, 1.0, n_samples=400, variation=0.2, L_flow=1.0, seed=7)
    target = ch01.box_average_density(L, 1e6, 0.2, 1.0)
    assert np.all(np.abs(mean - target) < 5 * std / np.sqrt(400) + 1e-9)
    m2, s2 = ch01.sample_density(1e-3, 1e15, 1.0, n_samples=50, gradient=0.3, L_flow=1.0, seed=2)
    assert abs(m2 / (1e15 * (1 + 0.3 * 1e-3 / 2)) - 1) < 1e-3


def test_kinetic_pressure_V2_derivation_D34():  # V2 (D34 re-derived: momentum flux of a Maxwellian gas)
    u, m, n, kB, T = sp.symbols("u m n k_B T", positive=True)
    sigma2 = kB * T / m
    f = sp.exp(-u**2 / (2 * sigma2)) / sp.sqrt(2 * sp.pi * sigma2)  # Gaussian u_x distribution
    p = 2 * m * n * sp.integrate(u**2 * f, (u, 0, sp.oo))  # steps 1-5: 2 m u_x per hit, flux n u_x f, u_x > 0 only
    assert sp.simplify(p - n * kB * T) == 0  # step 8
    ux2 = sp.integrate(u**2 * f, (u, -sp.oo, sp.oo))
    assert sp.simplify(3 * m * ux2 - 3 * kB * T) == 0  # ⟨|u|²⟩ = 3⟨u_x²⟩ = 3k_BT/m (step 7)


def test_kinetic_pressure_V1_sampled_molecules_give_nkT():  # V1 (statistical check of D34: momentum flux = n k_B T)
    T = 288.15
    u = ch01.maxwellian_velocities(1_000_000, T, M_AIR, seed=0)
    p_ref = N_AIR * ch01.K_B * T
    assert rel(ch01.molecular_pressure(N_AIR, M_AIR, u), p_ref) < 0.01  # rel SE 0.08 %
    assert rel(ch01.wall_impact_pressure(u, M_AIR, N_AIR), p_ref) < 0.012  # rel SE 0.22 %
    assert rel(np.var(u, axis=0), ch01.K_B * T / M_AIR) < 0.01
    assert rel(ch01.molecular_gas_pressure(N_AIR, 1.0, T), p_ref) < 1e-14
    c_bar = np.mean(np.linalg.norm(u, axis=1))
    assert rel(ch01.mean_molecular_speed(T, M_AIR), c_bar) < 0.005


def test_molecular_spacing_V1_liquid_vs_gas():  # V1 (C04)
    s_w = ch01.mean_molecular_spacing(1000.0, ch01.MOLAR_MASS["H2O"])
    s_a = ch01.mean_molecular_spacing(1.225, ch01.M_W_AIR)
    assert 0.30e-9 < s_w < 0.32e-9 and 3.3e-9 < s_a < 3.5e-9 and 10 < s_a / s_w < 12
    assert rel(ch01.number_density(1.225, ch01.M_W_AIR) ** (-1 / 3), s_a) < 1e-14


def test_knudsen_V5_jennings_mean_free_path():  # V5 (Jennings 1988 via Tsalikis et al. 2024) + V2 dims
    ref = ref_json("benchmarks.json")["jennings"]
    lam = ch01.mean_free_path_air(300.0)
    assert rel(lam * 1e9, ref["mean_free_path_air_300K_1atm_nm"]) < 0.01, lam
    q = np.sqrt(np.pi / 8) * Q_(1.85e-5, "Pa*s") / (ref["phi"] * (Q_(1.18, "kg/m**3") * Q_(101325, "Pa")) ** 0.5)
    assert q.check("[length]")
    assert rel(ch01.mean_free_path_jennings(1.85e-5, 1.18, 101325.0), q.to("m").magnitude) < 1e-12
    assert ch01.knudsen_number(67e-9, 1.0) == pytest.approx(6.7e-8) and ch01.knudsen_number(1e-7, 1e-7) == 1.0
    with pytest.raises(ValueError):
        ch01.knudsen_number(1e-7, 0.0)
    # altitude: l ∝ 1/p at fixed T (V7)
    assert rel(ch01.mean_free_path_air(250.0, 5e3) / ch01.mean_free_path_air(250.0, 5e4), 10.0) < 1e-12


@book_only
def test_knudsen_V6_book_mean_free_path_order_of_magnitude():  # V6 (book quotes an approximate value)
    b = book()["continuum"]["mean_free_path_air_room_nm"]["value"]
    lam = ch01.mean_free_path_air(293.15) * 1e9
    assert abs(math.log10(lam / b)) < 0.5  # same order of magnitude (reported as relative difference in the report)


def test_fluid_particle_V1_collision_time():  # V1 (C23/C24 numbers)
    tc = ch01.collision_time(ch01.mean_free_path_air(300.0), 300.0, M_AIR)
    assert 1e-10 < tc < 2e-10
    assert rel(tc, ch01.mean_free_path_air(300.0) / ch01.mean_molecular_speed(300.0, M_AIR)) < 1e-14


# ====================================================================================================================
# C12 · Newton's law of viscosity (§1.5) — A item; with C02, C08, C10, C11, N05, C13, C14
# ====================================================================================================================
def test_newton_viscosity_V1_quadratic_profile_exact_stress():  # V1 (τ = −2μU0y/b² exactly for a parabola)
    b, U0, mu = 0.01, 0.3, 1.0e-3
    y = np.linspace(-b, b, 41)
    u = U0 * (1 - (y / b) ** 2)
    tau = ch01.shear_stress_profile(u, y, mu)
    assert np.max(np.abs(tau - (-2 * mu * U0 * y / b**2))) < 1e-10 * mu * U0 / b
    assert ch01.newton_shear_stress(1e-3, 1000.0) == pytest.approx(1.0)  # worked number: water Couette 1 mm, 1 m/s
    with pytest.raises(ValueError):
        ch01.newton_shear_stress(-1.0, 1.0)


def test_newton_viscosity_V3_stress_stencils_second_order():  # V3
    hs, errs = [], []
    for N in (21, 41, 81, 161):
        y = np.linspace(0.0, 1.0, N)
        tau = ch01.shear_stress_profile(np.sin(3 * y), y, 2.0)
        errs.append(np.max(np.abs(tau - 6 * np.cos(3 * y))))  # includes the one-sided end stencils
        hs.append(y[1])
    p = observed_order(hs, errs)
    assert abs(p - 2.0) < 0.15, (p, pairwise_orders(hs, errs))
    ynu = np.array([0.0, 0.1, 0.3, 0.6, 1.0, 1.5])  # non-uniform grid: still exact for a quadratic
    assert np.max(np.abs(ch01.derivative_2nd_order(ynu**2, ynu) - 2 * ynu)) < 1e-12


def test_newton_viscosity_V2_units_of_transport_laws():  # V2 (dimensional homogeneity of (1.1)–(1.4))
    tau = Q_(1e-3, "Pa*s") * Q_(1000.0, "1/s")
    assert tau.check("[pressure]") and rel(ch01.newton_shear_stress(1e-3, 1000.0), tau.to("Pa").magnitude) < 1e-14
    nu = Q_(1.8e-5, "Pa*s") / Q_(1.2, "kg/m**3")
    assert nu.check("[length]**2/[time]") and rel(ch01.kinematic_viscosity(1.8e-5, 1.2), nu.to("m**2/s").magnitude) < 1e-14
    J = -Q_(1.2, "kg/m**3") * Q_(2.6e-5, "m**2/s") * Q_(3.0, "1/m")
    assert J.check("[mass]/[length]**2/[time]")
    assert rel(ch01.fick_mass_flux(1.2, 2.6e-5, 3.0), J.to("kg/(m**2*s)").magnitude) < 1e-14
    q = -Q_(0.6, "W/(m*K)") * Q_(50.0, "K/m")
    assert q.check("[power]/[length]**2") and rel(ch01.fourier_heat_flux(0.6, 50.0), q.to("W/m**2").magnitude) < 1e-14
    kap = Q_(0.6, "W/(m*K)") / (Q_(998, "kg/m**3") * Q_(4182, "J/(kg*K)"))
    assert kap.check("[length]**2/[time]") and rel(ch01.thermal_diffusivity(0.6, 998, 4182), kap.to("m**2/s").magnitude) < 1e-14


def test_transport_laws_V1_fluxes_point_down_gradient():  # V1 (C10, C11; second law iii)
    g = np.array([0.3, -1.0, 2.0])
    assert np.dot(ch01.fick_mass_flux(1.2, 2e-5, g), g) < 0 and np.dot(ch01.fourier_heat_flux(0.03, g), g) < 0
    with pytest.raises(ValueError):
        ch01.fick_mass_flux(1.2, -1e-5, g)
    with pytest.raises(ValueError):
        ch01.fourier_heat_flux(-0.1, g)


def test_couette_startup_V3_ftcs_converges_to_series_second_order():  # V3 (FTCS vs separation-of-variables series)
    h, U, nu, t_end, r = 1.0e-3, 1.0, 1.0e-6, 0.05, 0.4  # t in units: t ν/h² = 0.05
    hs, errs = [], []
    for N in (11, 21, 41, 81):
        y = np.linspace(0.0, h, N)
        dy = y[1] - y[0]
        nsteps = int(round(t_end / (r * dy**2 / nu)))
        u0 = np.zeros(N)
        u0[-1] = U
        F = ch01.ftcs_diffusion_1d(u0, nu, dy, t_end / nsteps, nsteps, save_every=nsteps)
        exact = ch01.couette_startup_profile(y, t_end, U, h, nu, nterms=None)
        errs.append(np.max(np.abs(F[-1] - exact)))
        hs.append(dy)
    p = observed_order(hs, errs)
    assert abs(p - 2.0) < 0.15, (p, errs, pairwise_orders(hs, errs))


def test_couette_startup_V1_steady_state_uniform_stress_and_limits():  # V1 + V7
    h, U, nu, mu = 1e-3, 1.0, 1e-6, 1e-3
    y = np.linspace(0, h, 51)
    assert np.max(np.abs(ch01.couette_startup_profile(y, 50.0, U, h, nu) - U * y / h)) < 1e-12  # t ≫ h²/ν
    u_init = ch01.couette_startup_profile(y, 0.0, U, h, nu)
    assert u_init[-1] == U and np.all(u_init[:-1] == 0)
    a = ch01.couette_startup_profile(y, 0.2, U, h, nu, nterms=200)
    b = ch01.couette_startup_profile(y, 0.2, U, h, nu, nterms=None)
    assert np.max(np.abs(a - b)) < 1e-12
    dy = y[1] - y[0]
    F = ch01.ftcs_diffusion_1d(np.where(np.isclose(y, h), U, 0.0), nu, dy, ch01.stable_time_step(nu, dy), 12000,
                               save_every=100)  # 12000 steps = 2.2 h²/ν: first mode decayed by e^(−π²·2.2)
    tb, tt = ch01.wall_shear_history(F, dy, mu)
    assert rel(tb[-1], mu * U / h) < 1e-6 and rel(tt[-1], mu * U / h) < 1e-6  # steady Couette: τ = μU/h at both walls
    assert np.all(np.diff(tb) >= -1e-12)  # bottom stress rises monotonically as momentum arrives
    assert ch01.diffusion_time(h, nu) == pytest.approx(1.0)  # worked number h²/ν = 1 s for water


def test_couette_startup_V1_initial_state_any_gap_width():  # V1 (t = 0 exact initial state; relative plate tolerance)
    for h in (1e-12, 1e-6, 1.0, 1e3):
        y = np.linspace(0.0, h, 11)
        u = ch01.couette_startup_profile(y, 0.0, 2.5, h, 1e-6)
        assert u[-1] == 2.5 and np.all(u[:-1] == 0.0), (h, u)  # only the moving plate has speed U
        assert ch01.couette_startup_profile(h, 0.0, 2.5, h, 1e-6) == 2.5
        assert ch01.couette_startup_profile(0.9 * h, 0.0, 2.5, h, 1e-6) == 0.0


def test_ftcs_V7_stability_limit_enforced():  # V7 (stability limit r = 1/2)
    N, D, dy = 41, 1.0, 1.0 / 40
    rng = np.random.default_rng(0)
    f0 = rng.standard_normal(N)
    f0[0] = f0[-1] = 0.0
    with pytest.raises(ValueError):
        ch01.ftcs_diffusion_1d(f0, D, dy, 0.51 * dy**2 / D, 10)
    stable = ch01.ftcs_diffusion_1d(f0, D, dy, 0.49 * dy**2 / D, 3000, save_every=3000)[-1]
    unstable = ch01.ftcs_diffusion_1d(f0, D, dy, 0.51 * dy**2 / D, 3000, save_every=3000, check_stability=False)[-1]
    assert np.max(np.abs(stable)) < np.max(np.abs(f0)) and np.max(np.abs(unstable)) > 1e3 * np.max(np.abs(f0))
    assert ch01.stable_time_step(D, dy, 0.9) * D / dy**2 == pytest.approx(0.45)
    assert ch01.ftcs_stable_time_step is ch01.stable_time_step


def test_ftcs_V4_zero_flux_conserves_integral_and_gaussian_order():  # V4 + V3 (heat-kernel exact solution)
    y = np.linspace(-3, 3, 121)
    dy = y[1] - y[0]
    f0 = np.exp(-(y - 0.5) ** 2 / 0.1)
    F = ch01.ftcs_diffusion_1d(f0, 0.5, dy, 0.45 * dy**2 / 0.5, 2000, bc=("neumann", "neumann"), save_every=500)
    integ = [np.trapezoid(f, dx=dy) for f in F]
    assert np.max(np.abs(np.array(integ) - integ[0])) < 1e-12 * integ[0]
    Fp = ch01.ftcs_diffusion_1d(f0, 0.5, dy, 0.45 * dy**2 / 0.5, 500, bc=("periodic", "periodic"), save_every=500)
    assert abs(np.sum(Fp[-1]) - np.sum(f0)) < 1e-10
    hs, errs = [], []
    for N in (61, 121, 241, 481):
        yy = np.linspace(-3, 3, N)
        d = yy[1] - yy[0]
        ns = int(round(0.02 / (0.25 * d**2)))
        f = ch01.ftcs_diffusion_1d(ch01.gaussian_spreading(yy, 0.0, 1.0, t0=0.01), 1.0, d, 0.02 / ns, ns,
                                   bc=("neumann", "neumann"), save_every=ns)[-1]
        errs.append(np.max(np.abs(f - ch01.gaussian_spreading(yy, 0.02, 1.0, t0=0.01))))
        hs.append(d)
    assert abs(observed_order(hs, errs) - 2.0) < 0.15, (errs, pairwise_orders(hs, errs))
    with pytest.raises(ValueError):
        ch01.ftcs_diffusion_1d(f0, 1.0, dy, 1e-6, 1, bc=("periodic", "dirichlet"))


def test_ftcs_V1_boundary_names_validated():  # V1 (unknown bc raises; names case-insensitive; periodic node convention)
    y = np.linspace(0.0, 1.0, 21)
    f0 = np.sin(np.pi * y) + 0.3
    dy, dt = y[1], 0.4 * y[1] ** 2
    for bad in (("dirichlet", "robin"), ("neuman", "neumann"), ("", "dirichlet"), "neumann", ("dirichlet",)):
        with pytest.raises(ValueError):
            ch01.ftcs_diffusion_1d(f0, 1.0, dy, dt, 5, bc=bad)
    ref = ch01.ftcs_diffusion_1d(f0, 1.0, dy, dt, 50, bc=("neumann", "dirichlet"))
    assert np.array_equal(ch01.ftcs_diffusion_1d(f0, 1.0, dy, dt, 50, bc=(" Neumann", "DIRICHLET")), ref)
    assert np.all(ref[:, -1] == f0[-1])  # Dirichlet end held
    # periodic: node N−1 neighbours node 0; a sine sampled without the repeated endpoint decays as exp(−D k'² t) with
    # the discrete FTCS factor (1 − 4r sin²(k dy/2)) per step, k = 2π/(N dy)
    N = 32
    yp = dy * np.arange(N)
    k = 2 * np.pi / (N * dy)
    r = 0.4
    Fp = ch01.ftcs_diffusion_1d(np.sin(k * yp), 1.0, dy, r * dy**2, 40, bc=("periodic", "periodic"), save_every=40)
    amp = (1 - 4 * r * np.sin(k * dy / 2) ** 2) ** 40
    assert np.max(np.abs(Fp[-1] - amp * np.sin(k * yp))) < 1e-13


def test_viscosity_temperature_V5_sutherland_matches_ussa_table():  # V5 (USSA-1976 Table 2 via PDAS) + V7
    for row in ref_csv("ussa1976_table2.csv"):
        T_row = [r for r in ref_csv("ussa1976_table1.csv") if r["Z_km"] == row["Z_km"]]
        if not T_row:
            continue
        T = T_row[0]["T_K"]
        assert rel(ch01.sutherland_viscosity(T), row["mu_Pa_s"]) < 1e-4, row
        assert rel(ch01.kinematic_viscosity(row["mu_Pa_s"], T_row[0]["rho_kg_m3"]), row["nu_m2_s"]) < 2e-4, row
    c = ref_json("ussa1976_constants.json")
    assert ch01.sutherland_viscosity(300.0) == ch01.sutherland_viscosity(300.0, c["sutherland_beta_kg_per_m_s_sqrtK"],
                                                                         c["sutherland_S_K"])
    T = np.linspace(200, 400, 21)  # gas: μ rises with T; liquid: μ falls
    assert np.all(np.diff(ch01.sutherland_viscosity(T)) > 0)
    assert np.all(np.diff(ch01.water_viscosity(np.linspace(275, 370, 20))) < 0)
    assert ch01.viscosity_power_law(300.0, 1.8e-5, 300.0) == 1.8e-5
    assert rel(ch01.viscosity_power_law(400.0, 1.0, 100.0), 2.0) < 1e-14  # n = 1/2


def test_water_viscosity_V5_iapws_2008_check_point():  # V5 (IAPWS R12-08 Table 4, 298.15 K)
    chk = ref_json("benchmarks.json")["iapws_viscosity_check"][0]
    assert rel(ch01.water_viscosity(chk["T_K"]), chk["mu_microPa_s"] * 1e-6) < 0.01  # fit claims ±1 %
    with pytest.raises(ValueError):
        ch01.water_viscosity(100.0)


def test_fluid_properties_V5_air_matches_ussa_sea_level():  # V5 + V1 (C14: air ν ≈ 15× water ν)
    t1 = ref_csv("ussa1976_table1.csv")[0]
    t2 = ref_csv("ussa1976_table2.csv")[0]
    air = ch01.fluid_properties("air", T=288.15, p=101325.0)
    assert rel(air["rho"], t1["rho_kg_m3"]) < 1e-4 and rel(air["mu"], t2["mu_Pa_s"]) < 1e-4
    assert rel(air["nu"], t2["nu_m2_s"]) < 2e-4
    w = ch01.fluid_properties("water")
    assert rel(w["nu"], w["mu"] / w["rho"]) < 1e-14 and "sigma" in w
    with pytest.raises(ValueError):
        ch01.fluid_properties("mercury")
    for name, f in ch01.FLUIDS.items():
        assert rel(f["nu"], f["mu"] / f["rho"]) < 1e-14 and f["kappa_m"] > 0 and f["k"] > 0 and f["cp"] > 0, name
    assert 14 < ch01.FLUIDS["air"]["nu"] / ch01.FLUIDS["water"]["nu"] < 16
    assert ch01.FLUIDS["honey"]["mu"] / ch01.FLUIDS["air"]["mu"] > 1e5


def test_fluid_vs_solid_V1_deformation_histories():  # V1 + V7 (C02, N01)
    t = np.linspace(-1, 10, 111)
    tau, G, mu = 2.0, 50.0, 4.0
    sol = ch01.shear_deformation_history(t, tau, "solid", G=G, t_off=5.0)
    flu = ch01.shear_deformation_history(t, tau, "fluid", mu=mu, t_off=5.0)
    assert np.allclose(sol[(t >= 0) & (t < 5)], tau / G) and np.all(sol[t >= 5] == 0) and np.all(sol[t < 0] == 0)
    assert np.allclose(flu[(t >= 0) & (t <= 5)], tau * t[(t >= 0) & (t <= 5)] / mu)
    assert np.allclose(flu[t >= 5], tau * 5.0 / mu)  # the fluid keeps its strain (no memory)
    assert np.all(ch01.shear_deformation_history(t, 1.0, "bingham", mu=mu, tau_y=1.5) == 0)  # below yield
    mx = ch01.shear_deformation_history(np.array([4.9, 5.1]), tau, "maxwell", G=G, mu=mu, t_off=5.0)
    assert mx[0] - mx[1] == pytest.approx(tau / G, rel=1e-6, abs=0.05)  # elastic part recovers
    assert ch01.shear_deformation_history(3.0, tau, "fluid", mu=1e12) < 1e-10  # μ → ∞: no flow
    for bad in (dict(kind="solid"), dict(kind="fluid"), dict(kind="glass", mu=1.0)):
        with pytest.raises(ValueError):
            ch01.shear_deformation_history(1.0, 1.0, **bad)


def test_normal_shear_stress_V1_traction_split():  # V1 (C03)
    rng = np.random.default_rng(5)
    t, n = rng.standard_normal((50, 3)), rng.standard_normal((50, 3))
    s_n, tau_vec, tau_mag = ch01.traction_components(t, n)
    nh = n / np.linalg.norm(n, axis=1, keepdims=True)
    assert np.max(np.abs(np.sum(tau_vec * nh, axis=1))) < 1e-12
    assert np.max(np.abs(s_n[:, None] * nh + tau_vec - t)) < 1e-12
    s1, _, m1 = ch01.traction_components([0, 0, -5.0], [0, 0, 2.0])  # pure pressure on a face
    assert s1 == -5.0 and m1 == 0.0


def test_mass_fraction_V1_dry_air_composition():  # V1 (C09) + V5 (USSA Table 3 mean molecular weight)
    X = {k: v[0] for k, v in ch01.DRY_AIR.items()}
    M = {k: v[1] for k, v in ch01.DRY_AIR.items()}
    Y = ch01.mass_fractions(X, M)
    assert abs(sum(Y.values()) - 1) < 1e-14
    Mmix = 1.0 / sum(Y[k] / M[k] for k in Y)
    assert abs(Mmix - 28.96) < 0.01
    arr = ch01.mass_fractions(np.array([1.0, 1.0]), np.array([2.0, 6.0]))
    assert np.allclose(arr, [0.25, 0.75])
    comp = ref_json("ussa1976_constants.json")["composition"]
    M0 = sum(Mi * Fi for Mi, Fi in comp.values()) / sum(Fi for _, Fi in comp.values())
    assert rel(ch01.M_W_AIR, M0) < 5e-6, M0  # USSA sea-level mean molecular weight 28.9644 kg/kmol


# ====================================================================================================================
# §1.6 Surface tension (B items under C20)
# ====================================================================================================================
def test_surface_tension_V5_iapws_table():  # V5 (IAPWS R1-76(2014) Table 1) + V7
    for row in ref_csv("iapws_sigma.csv"):
        s = ch01.surface_tension_water(row["t_C"] + 273.15) * 1e3
        assert abs(s - row["sigma_calc_mN_m"]) <= 0.005 + 1e-9, row  # calculated column, 2-decimal rounding
        assert abs(s - row["sigma_exp_mN_m"]) <= row["uncertainty_mN_m"], row  # experimental ± uncertainty
    T = np.linspace(273.16, 647.0, 200)
    assert np.all(np.diff(ch01.surface_tension_water(T)) < 0) and ch01.surface_tension_water(647.096) == 0.0
    with pytest.raises(ValueError):
        ch01.surface_tension_water(700.0)


def test_laplace_jump_V1_special_cases_and_units():  # V1 + V2 (C16, N07)
    s, R = 0.0728, 1e-3
    assert rel(ch01.laplace_pressure_jump(s, R), 2 * s / R) < 1e-14
    assert ch01.laplace_pressure_jump(s, np.inf, np.inf) == 0.0
    assert abs(ch01.laplace_pressure_jump(s, R, -R)) < 1e-12
    assert rel(ch01.laplace_pressure_jump(s, R, np.inf), s / R) < 1e-14
    dp = Q_(s, "N/m") * (1 / Q_(R, "m") + 1 / Q_(R, "m"))
    assert dp.check("[pressure]") and 140 < dp.to("Pa").magnitude < 150  # worked number ≈ 146 Pa


def test_capillary_rise_V1_force_balance_and_laplace_consistency():  # V1 (Ex. 1.1) + V7
    s, rho, R, alpha = 0.0728, 998.0, 1e-3, np.radians(70.0)
    h = ch01.capillary_rise(s, alpha, rho, R)
    assert abs(s * 2 * np.pi * R * np.sin(alpha) - rho * G * h * np.pi * R**2) < 1e-14
    p_E = ch01.P_ATM - ch01.laplace_pressure_jump(s, R / np.sin(alpha))  # meniscus radius R/sin α
    assert rel(ch01.hydrostatic_pressure_uniform(h, ch01.P_ATM, rho), p_E) < 1e-12
    assert ch01.capillary_rise(s, 0.0, rho, R) == 0.0
    assert rel(ch01.capillary_rise(s, alpha, rho, R / 3), 3 * h) < 1e-14
    assert rel(ch01.capillary_rise_deg(s, 70.0, rho, R), h) < 1e-14
    assert ch01.alpha_from_contact_angle(0.0) == pytest.approx(np.pi / 2)
    h_full = ch01.capillary_rise(s, ch01.alpha_from_contact_angle(0.0), rho, R)
    assert 0.0145 < h_full < 0.0152  # worked number ≈ 15 mm
    assert ch01.contact_angle_to_alpha is ch01.alpha_from_contact_angle


# ====================================================================================================================
# C20 · Hydrostatic law (§1.7) — A item; with C17, C18, C19, C21, N11, C22, D05, D37
# ====================================================================================================================
def test_hydrostatics_V2_uniform_solution_residual_and_units():  # V2
    z, p0, rho, g = sp.symbols("z p0 rho g", real=True)
    p = p0 - rho * g * z  # (1.9)
    assert sp.simplify(sp.diff(p, z) + rho * g) == 0 and p.subs(z, 0) == p0  # (1.8) and p(0) = p0
    f = sp.lambdify((z, p0, rho, g), p, "numpy")
    zz = np.linspace(-50, 10, 61)
    assert np.max(np.abs(ch01.hydrostatic_pressure_uniform(zz, 101325.0, 1025.0) - f(zz, 101325.0, 1025.0, G))) < 1e-9
    q = Q_(101325, "Pa") - Q_(1000, "kg/m**3") * Q_(G, "m/s**2") * Q_(-10, "m")
    assert q.check("[pressure]")
    assert rel(ch01.hydrostatic_pressure_uniform(-10.0, 101325.0, 1000.0), q.to("Pa").magnitude) < 1e-14
    assert 0.97e5 < ch01.hydrostatic_pressure_uniform(-10.0, 0.0, 1000.0) < 0.99e5  # 10 m of water ≈ 1 atm


def test_hydrostatics_V1_integrator_matches_closed_forms():  # V1 (uniform, isothermal, linear lapse)
    z = np.linspace(0.0, 11000.0, 23)
    p_u = ch01.integrate_hydrostatic(-z / 1000, lambda zz, p: 1000.0, 101325.0, z0=0.0)
    assert rel(p_u, ch01.hydrostatic_pressure_uniform(-z / 1000, 101325.0, 1000.0)) < 1e-10
    p_iso = ch01.integrate_hydrostatic(z, lambda zz, p: p / (ch01.R_AIR * 250.0), 101325.0)
    assert rel(p_iso, ch01.isothermal_pressure(z, 101325.0, 250.0)) < 1e-8
    p, rho, T = ch01.atmosphere_from_temperature(z, lambda zz: 288.15 - 6.5e-3 * np.asarray(zz), 101325.0)
    assert rel(p, ch01.linear_lapse_pressure(z, 101325.0, 288.15, -6.5e-3)) < 1e-8
    assert rel(rho, p / (ch01.R_AIR * T)) < 1e-14
    zmid = np.linspace(-2000, 2000, 9)  # integrate both ways from a middle level
    pm = ch01.integrate_hydrostatic(zmid, lambda zz, pp: pp / (ch01.R_AIR * 270.0), 9e4, z0=0.0)
    assert rel(pm, ch01.isothermal_pressure(zmid, 9e4, 270.0)) < 1e-8
    assert rel(ch01.linear_lapse_pressure(z, 1e5, 280.0, -1e-9), ch01.isothermal_pressure(z, 1e5, 280.0)) < 1e-5  # Γ → 0


def test_hydrostatics_V3_from_scratch_euler_first_order():  # V3 (from-scratch Euler march, design §7)
    H, T = 3000.0, 260.0
    hs, errs = [], []
    for n in (50, 100, 200, 400):
        z = np.linspace(0.0, H, n + 1)
        dz = z[1] - z[0]
        p = np.empty_like(z)
        p[0] = 1e5
        for k in range(n):
            p[k + 1] = p[k] - p[k] / (ch01.R_AIR * T) * G * dz  # p[k+1] = p[k] − ρ(z_k) g Δz
        errs.append(abs(p[-1] - ch01.integrate_hydrostatic(z, lambda zz, pp: pp / (ch01.R_AIR * T), 1e5)[-1]))
        hs.append(dz)
    assert abs(observed_order(hs, errs) - 1.0) < 0.15


@needs_ref
def test_hydrostatics_V5_ussa1976_table_profiles():  # V5 (PDAS USSA-1976 Table 1, geometric → geopotential)
    c = ref_json("ussa1976_constants.json")
    assert np.allclose(ch01.USSA_BASES / 1e3, c["layers_H_km"][:-1]) and np.allclose(ch01.USSA_LAPSE * 1e3,
                                                                                      c["layers_L_K_per_km"])
    assert ch01.USSA_Z_TOP / 1e3 == pytest.approx(c["layers_H_km"][-1])
    worst = {"T": 0.0, "p": 0.0, "rho": 0.0}
    for row in ref_csv("ussa1976_table1.csv"):
        H = geopotential_m(row["Z_km"], c["r0_km"])
        T, p, rho = ch01.standard_atmosphere(H)
        assert abs(T - row["T_K"]) < 1e-3, row
        assert rel(p, row["p_Pa"]) < 1e-3 and rel(rho, row["rho_kg_m3"]) < 1e-3, (row, p, rho)
        worst = {"T": max(worst["T"], abs(T - row["T_K"])), "p": max(worst["p"], rel(p, row["p_Pa"])),
                 "rho": max(worst["rho"], rel(rho, row["rho_kg_m3"]))}
    print("USSA worst errors", worst)
    for Zk in (5.0, 10.0):  # the integrator itself (C49) through the 0–11 km layer
        row = [r for r in ref_csv("ussa1976_table1.csv") if r["Z_km"] == Zk][0]
        H = geopotential_m(Zk, c["r0_km"])
        p, _, _ = ch01.atmosphere_from_temperature(np.array([0.0, H]), lambda zz: 288.15 - 6.5e-3 * np.asarray(zz))
        assert rel(p[-1], row["p_Pa"]) < 1e-3
    with pytest.raises(ValueError):
        ch01.standard_atmosphere(9e4)


@needs_ref
def test_standard_atmosphere_V5_geometric_height_option():  # V5 (PDAS USSA-1976 Table 1 is tabulated at geometric Z)
    c = ref_json("ussa1976_constants.json")
    assert ch01.USSA_R0 == pytest.approx(c["r0_km"] * 1e3, rel=1e-15)  # r0 = 6356.766 km (Table 2 + NASA errata)
    rows = ref_csv("ussa1976_table1.csv")
    Z = np.array([r["Z_km"] * 1e3 for r in rows])
    T, p, rho = ch01.standard_atmosphere(Z, geometric=True)
    assert np.max(np.abs(T - [r["T_K"] for r in rows])) < 1e-3
    assert rel(p, [r["p_Pa"] for r in rows]) < 1e-3 and rel(rho, [r["rho_kg_m3"] for r in rows]) < 1e-3
    T10 = ch01.standard_atmosphere(10000.0, geometric=True)[0]
    assert abs(T10 - 223.252) < 1e-3  # Table 1 at Z = 10 km geometric
    assert abs(ch01.standard_atmosphere(10000.0)[0] - 223.25) > 0.1  # 10 km geopotential is a different level (223.15 K)
    Tg, pg, rg = ch01.standard_atmosphere(geopotential_m(Z / 1e3, c["r0_km"]))  # same as converting by hand
    assert rel(T, Tg) < 1e-14 and rel(p, pg) < 1e-12 and rel(rho, rg) < 1e-12
    assert np.isfinite(ch01.standard_atmosphere(86000.0, geometric=True)[0])  # top of the model maps inside the range
    for bad in (-1.0, 86001.0):
        with pytest.raises(ValueError):
            ch01.standard_atmosphere(bad, geometric=True)


def test_hydrostatics_V1_layered_tank_and_gauge():  # V1 (C17, C21)
    z = np.linspace(-6, 1, 71)
    one = ch01.layered_pressure(z, [10.0], [1000.0])
    assert rel(one[z <= 0], ch01.hydrostatic_pressure_uniform(z[z <= 0], ch01.P_ATM, 1000.0)) < 1e-14
    assert np.all(one[z > 0] == ch01.P_ATM)
    two = ch01.layered_pressure(z, [2.0, 10.0], [800.0, 1000.0])
    zi = np.array([-2.0 - 1e-9, -2.0 + 1e-9])
    assert abs(np.diff(ch01.layered_pressure(zi, [2.0, 10.0], [800.0, 1000.0]))[0]) < 1e-4  # continuous at interface
    slope = np.gradient(two, z)
    assert abs(slope[10] + 1000 * G) < 1e-6 and abs(slope[55] + 800 * G) < 1e-6
    assert ch01.gauge_pressure(ch01.P_ATM) == 0 and ch01.gauge_pressure(0.0) == -ch01.P_ATM
    assert ch01.absolute_pressure(ch01.gauge_pressure(2.5e5)) == 2.5e5
    with pytest.raises(ValueError):
        ch01.layered_pressure(z, [1.0], [1.0, 2.0])


def test_buoyancy_V1_face_integral_equals_archimedes_D37():  # V1 (D37) + V4 (Pascal: horizontal forces cancel)
    rho = 1025.0
    pfn = lambda x, y, zz: ch01.hydrostatic_pressure_uniform(zz, ch01.P_ATM, rho)  # noqa: E731
    box = (0.1, 0.3, -0.2, 0.3, -5.0, -4.6)
    V = 0.2 * 0.5 * 0.4
    F = ch01.net_pressure_force_on_box(pfn, box)
    assert rel(F[2], ch01.buoyancy_force(rho, V)) < 1e-10
    assert np.max(np.abs(F[:2])) < 1e-9 * ch01.P_ATM * 0.2
    assert rel(ch01.buoyancy_force(1000.0, 1e-3), 9.80665) < 1e-14  # 1 L in water


def test_buoyancy_V2_derivation_D37():  # V2 (D37 steps 4–6)
    p0, rho, g, z1, z2, A = sp.symbols("p0 rho g z1 z2 A", real=True)
    p = lambda z: p0 - rho * g * z  # noqa: E731
    F_net = (p(z1) - p(z2)) * A
    assert sp.simplify(F_net - rho * g * A * (z2 - z1)) == 0


def test_pressure_isotropy_V1_wedge_closed_form_vanishes_linearly():  # V1 (C18: closed-form wedge balance, ∝ dz; not a discretisation)
    dz = np.array([1e-1, 5e-2, 2.5e-2, 1.25e-2])
    for th in (0.2, 0.7, 1.3):
        d = ch01.wedge_pressure_difference(1000.0, dz, th)
        assert np.all(d["p3_minus_p1"] == 0)
        assert abs(observed_order(dz, d["p2_minus_p1"]) - 1.0) < 0.05
        assert rel(d["p2_minus_p1"], 0.5 * 1000.0 * G * dz) < 1e-14
        f = ch01.wedge_face_forces(1e5, 1000.0, dz, th)
        assert np.max(np.abs(f["res_x"])) < 1e-8 and np.max(np.abs(f["res_z"])) < 1e-8
        assert abs(observed_order(dz, f["weight_to_face_force"]) - 1.0) < 0.05


# ====================================================================================================================
# C25 · First law (§1.8) — A item; with C23, C26, C27, C28
# ====================================================================================================================
def test_first_law_V1_isothermal_and_isochoric_paths():  # V1
    R, cv, T = ch01.R_AIR, ch01.CV_AIR, 300.0
    d = ch01.process_path("isothermal", (0.8, T), (1.6, T), n=4001)
    out = ch01.process_heat_work(d["v"], d["T"])
    assert rel(out["w"][-1], -R * T * np.log(2)) < 1e-6 and rel(out["q"][-1], R * T * np.log(2)) < 1e-6
    assert abs(out["de"][-1]) < 1e-12 and -59.8e3 < out["w"][-1] < -59.6e3  # worked number −59.7 kJ/kg
    d2 = ch01.process_path("isochoric", (0.8, 300.0), (0.8, 450.0), n=11)
    o2 = ch01.process_heat_work(d2["v"], d2["T"])
    assert np.all(o2["w"] == 0) and rel(o2["q"][-1], cv * 150.0) < 1e-14
    assert rel(o2["p"], R * d2["T"] / 0.8) < 1e-14
    with pytest.raises(ValueError):
        ch01.process_path("isothermal", (0.8, 300.0), (1.6, 310.0))
    with pytest.raises(ValueError):
        ch01.process_heat_work([1.0], [300.0])


def test_first_law_V3_trapezoid_work_second_order():  # V3
    hs, errs = [], []
    for n in (11, 21, 41, 81):
        d = ch01.process_path("isentropic", (0.8, 300.0), (1.6, 300.0 * 0.5**0.4), n=n)
        w = ch01.process_heat_work(d["v"], d["T"])["w"][-1]
        errs.append(abs(w - ch01.path_heat_work_totals("isentropic", 0.8, 300.0, 1.6, 300.0 * 0.5**0.4)["w"]))
        hs.append(0.8 / (n - 1))
    assert abs(observed_order(hs, errs) - 2.0) < 0.15, pairwise_orders(hs, errs)


def test_first_law_V4_state_vs_path_functions():  # V4 (Δe, Δs path-independent; q, w not)
    v1, T1, v2, T2 = 0.8, 300.0, 1.6, 420.0
    kinds = ["isochoric-isobaric", "isobaric-isochoric", "isentropic-isochoric", "isothermal-isochoric"]
    tots = {k: ch01.path_heat_work_totals(k, v1, T1, v2, T2) for k in kinds}
    de = [t["de"] for t in tots.values()]
    ds = [t["ds"] for t in tots.values()]
    assert np.ptp(de) < 1e-9 and np.ptp(ds) < 1e-12
    assert rel(ds[0], ch01.perfect_gas_entropy_change(T1, v1, T2, v2)) < 1e-12
    assert np.ptp([t["q"] for t in tots.values()]) > 1e4 and np.ptp([t["w"] for t in tots.values()]) > 1e4
    for k in kinds:  # closed forms = trapezoid on fine sampled paths (and their q, w, Δs)
        d = ch01.process_path(k, (v1, T1), (v2, T2), n=4001)
        o = ch01.process_heat_work(d["v"], d["T"])
        assert d["corner"] == 4000 and rel(o["q"][-1], tots[k]["q"]) < 1e-5 and rel(o["w"][-1], tots[k]["w"]) < 1e-5
        s = ch01.entropy_change_reversible(o["q"], d["T"])
        assert rel(s, tots[k]["ds"]) < 1e-5
    assert ch01.path_heat_work_totals("isobaric", 0.8, 300.0, 1.2, 450.0)["w"] == pytest.approx(-ch01.R_AIR * 150.0)
    assert ch01.path_heat_work_totals("isochoric", 0.8, 300.0, 0.8, 450.0)["w"] == 0.0
    iso = ch01.path_heat_work_totals("isothermal", 0.8, 300.0, 1.6, 300.0)
    assert iso["q"] == pytest.approx(-iso["w"]) and iso["de"] == 0.0
    with pytest.raises(ValueError):
        ch01.path_heat_work_totals("zigzag", v1, T1, v2, T2)


def test_equation_of_state_V1_two_properties_fix_state():  # V1 (C28) + V7
    p, rho, T = ch01.perfect_gas_state(p=9e4, T=270.0)
    assert rel(ch01.perfect_gas_state(p=p, rho=rho)[2], 270.0) < 1e-14
    assert rel(ch01.perfect_gas_state(rho=rho, T=T)[0], 9e4) < 1e-14
    for bad in ({}, {"p": 1.0}, {"p": 1.0, "rho": 1.0, "T": 1.0}):
        with pytest.raises(ValueError):
            ch01.perfect_gas_state(**bad)
    assert ch01.specific_volume(4.0) == 0.25
    e = ch01.perfect_gas_internal_energy(300.0)
    h = ch01.perfect_gas_enthalpy(300.0)
    assert rel(h - e, ch01.R_AIR * 300.0) < 1e-12 and rel(ch01.enthalpy(e, 1e5, ch01.R_AIR * 300.0 / 1e5), h) < 1e-12


# ====================================================================================================================
# C35 · Gibbs relations (§1.8) — A item; with C29–C34, N15
# ====================================================================================================================
def test_gibbs_V2_derivation_D10_and_residuals():  # V2 (D10 steps + perfect-gas residuals of both forms)
    de, dq, dv, dp, p, v, T, ds = sp.symbols("de dq dv dp p v T ds")
    dq_ = de + p * dv  # step 1
    step2 = sp.Eq(T * ds, dq_)  # step 2
    dh = de + p * dv + v * dp  # steps 3–4 (product rule)
    assert sp.simplify(step2.rhs - (dh - v * dp)) == 0  # steps 5–6: T ds = dh − v dp
    Ts, vs, ps, cv, R = sp.symbols("T v p c_v R", positive=True)
    s_Tv = cv * sp.log(Ts) + R * sp.log(vs)
    e = cv * Ts
    p_Tv = R * Ts / vs
    assert sp.simplify(Ts * sp.diff(s_Tv, Ts) - sp.diff(e, Ts)) == 0  # dT coefficient of T ds − de − p dv
    assert sp.simplify(Ts * sp.diff(s_Tv, vs) - p_Tv) == 0  # dv coefficient
    cp = cv + R
    s_Tp = cp * sp.log(Ts) - R * sp.log(ps)
    assert sp.simplify(Ts * sp.diff(s_Tp, Ts) - cp) == 0 and sp.simplify(Ts * sp.diff(s_Tp, ps) + R * Ts / ps) == 0


def test_gibbs_V1_both_forms_agree_numerically():  # V1 (C35 worked number)
    ds_v = ch01.perfect_gas_entropy_change(300.0, 0.861, 600.0, 1.722)  # constant p ⇒ v doubles with T
    ds_p = ch01.perfect_gas_entropy_change_p(300.0, 1e5, 600.0, 1e5)
    assert rel(ds_v, ds_p) < 1e-12 and rel(ds_p, ch01.CP_AIR * np.log(2)) < 1e-14 and 696 < ds_p < 697
    T1, v1 = 300.0, 0.9
    v2 = 0.5
    T2 = T1 * (v1 / v2) ** 0.4  # isentrope from (1.25)/(1.26)
    assert abs(ch01.perfect_gas_entropy_change(T1, v1, T2, v2)) < 1e-12


def test_gibbs_V4_closed_reversible_cycle_returns_zero_entropy():  # V4
    legs = [ch01.process_path("isothermal", (0.8, 300.0), (1.6, 300.0), n=2001),
            ch01.process_path("isochoric", (1.6, 300.0), (1.6, 450.0), n=2001),
            ch01.process_path("isobaric", (1.6, 450.0), (0.8, 225.0), n=2001),
            ch01.process_path("isochoric", (0.8, 225.0), (0.8, 300.0), n=2001)]
    v, T = ch01.join_paths(*[(d["v"], d["T"]) for d in legs])
    o = ch01.process_heat_work(v, T)
    s = ch01.entropy_change_reversible(o["q"], T, cumulative=True)
    assert abs(s[-1]) < 1e-5 and abs(o["de"][-1]) < 1e-9 and abs(o["q"][-1] + o["w"][-1]) < 1e-9
    assert abs(o["w"][-1]) > 1e3  # net work around the cycle ≠ 0 (path function)
    with pytest.raises(ValueError):
        ch01.join_paths((v[:5], T[:5]), (v[10:], T[10:]))


def test_clausius_duhem_V1_irreversible_processes():  # V1 + V4 (C34)
    st = ch01.irreversible_process("stirring", 300.0, 0.8, T2=330.0)
    fe = ch01.irreversible_process("free_expansion", 300.0, 0.8, v2=1.6)
    assert st["q"] == 0 and rel(st["ds"], ch01.CV_AIR * np.log(1.1)) < 1e-14 and st["ds"] > st["int_dq_over_T"]
    assert rel(fe["ds"], ch01.R_AIR * np.log(2)) < 1e-14 and fe["ds"] > fe["int_dq_over_T"] == 0.0
    assert ch01.free_expansion(0.8, 1.6, 300.0) == fe
    sp_ = ch01.stirred_isochoric_process(0.8, 300.0, 330.0)
    assert rel(sp_["s"][-1], st["ds"]) < 1e-12 and np.all(sp_["q"] == 0) and ch01.stirred_isochoric_path is ch01.stirred_isochoric_process
    for bad in (("stirring", dict(T2=290.0)), ("free_expansion", dict(v2=0.1)), ("teleport", {})):
        with pytest.raises(ValueError):
            ch01.irreversible_process(bad[0], 300.0, 0.8, **bad[1])


def test_specific_heats_V1_partial_derivative_definitions():  # V1 (C30, C31) + V3 (central difference order 2)
    assert rel(ch01.specific_heat_cp(lambda T, p: ch01.perfect_gas_enthalpy(T), 300.0, 1e5), ch01.CP_AIR) < 1e-6
    a, b = ch01.VDW_CO2["a"], ch01.VDW_CO2["b"]
    Rc = ch01.gas_constant(ch01.MOLAR_MASS["CO2"])
    cv = 650.0
    assert rel(ch01.specific_heat_cv(lambda T, v: ch01.van_der_waals_internal_energy(T, v, a, cv), 300.0, 0.02), cv) < 1e-6
    steps = [1e-1, 5e-2, 2.5e-2, 1.25e-2]
    errs = [abs(ch01.partial_derivative(lambda x, y: np.sin(x) * y, "x", {"x": 0.7, "y": 2.0}, step=s) - 2 * np.cos(0.7))
            for s in steps]
    assert abs(observed_order(steps, errs) - 2.0) < 0.15
    # van der Waals contrast (C41): (∂e/∂v)_T = T(∂p/∂T)_v − p from Gibbs; zero for a = 0
    T, v = 300.0, 0.02
    lhs = ch01.partial_derivative(lambda T, v: ch01.van_der_waals_internal_energy(T, v, a, cv), "v", {"T": T, "v": v})
    rhs = T * ch01.partial_derivative(lambda T, v: ch01.van_der_waals_pressure(T, v, a, b, Rc), "T", {"T": T, "v": v}) \
        - ch01.van_der_waals_pressure(T, v, a, b, Rc)
    assert rel(lhs, rhs) < 1e-5 and rel(lhs, a / v**2) < 1e-6
    assert rel(ch01.van_der_waals_pressure(T, v, 0.0, 0.0, Rc), ch01.perfect_gas_pressure(1 / v, T, Rc)) < 1e-14
    am, bm = ch01.van_der_waals_constants_per_mass(0.3640, 4.267e-5, ch01.MOLAR_MASS["CO2"])
    assert rel(am, a) < 1e-12 and rel(bm, b) < 1e-12


# ====================================================================================================================
# C36 · Speed of sound (§1.8) — A item; with C37
# ====================================================================================================================
def test_sound_speed_V1_polytropic_eos_exact():  # V1
    K, gam = 101325.0 / 1.225**1.4, 1.4
    for rho in (0.3, 1.225, 5.0):
        c = ch01.sound_speed_from_eos(lambda r, s: K * r**gam, rho)
        assert rel(c, np.sqrt(gam * K * rho**gam / rho)) < 1e-8
    assert rel(ch01.sound_speed_from_eos(lambda r, s: K * r**gam, 1.225), 340.29) < 1e-4  # worked number ≈ 340 m/s
    with pytest.raises(ValueError):
        ch01.sound_speed_from_eos(lambda r, s: -r, 1.0)


def test_sound_speed_V5_ussa_table_and_perfect_gas():  # V5 (PDAS USSA-1976 c column) + V1 cross-check
    for row in ref_csv("ussa1976_table1.csv"):
        assert rel(ch01.perfect_gas_sound_speed(row["T_K"]), row["c_m_s"]) < 2e-4, row
    T = 288.15
    p = 101325.0
    rho = ch01.perfect_gas_density(p, T)
    c_eos = ch01.sound_speed_from_eos(lambda r, s: ch01.isentropic_pressure(r, p, rho), rho)
    assert rel(c_eos, ch01.perfect_gas_sound_speed(T)) < 1e-8
    assert rel(ch01.perfect_gas_sound_speed(T) / np.sqrt(ch01.R_AIR * T), np.sqrt(1.4)) < 1e-14  # not Newton's √(RT)


def test_sound_speed_V7_incompressible_limit_tait():  # V7 + V1
    rho0 = 1000.0
    cs = [ch01.sound_speed_from_eos(lambda r, s, K0=K0: ch01.tait_pressure(r, rho0, K0), rho0) for K0 in (2.2e7, 2.2e9, 2.2e11)]
    assert rel(cs[1], np.sqrt(2.2e9 / rho0)) < 1e-6 and 1470 < cs[1] < 1500
    assert cs[0] < cs[1] < cs[2] and cs[2] / cs[1] == pytest.approx(10.0, rel=1e-6)  # c ∝ √K0 → ∞
    assert ch01.tait_pressure(rho0) == ch01.P_ATM
    with pytest.raises(ValueError):
        ch01.tait_pressure(1000.0, K0=0.0)


def test_thermal_expansion_V1_perfect_gas_linear_eos_and_water():  # V1 (C37, C48)
    for T in (200.0, 300.0, 400.0):
        a = ch01.thermal_expansion_coefficient(lambda T, p: ch01.perfect_gas_density(p, T), T, 1e5)
        assert rel(a, ch01.perfect_gas_expansion_coefficient(T)) < 1e-6 and rel(a, 1 / T) < 1e-6
    a0, T0 = 2e-4, 283.15
    lin = lambda T, p: 1027.0 * (1 - a0 * (T - T0))  # noqa: E731
    assert rel(ch01.thermal_expansion_coefficient(lin, 293.15, 1e5), a0 / (1 - a0 * 10)) < 1e-6
    wd = lambda T, p: ch01.water_density(T)  # noqa: E731
    assert ch01.thermal_expansion_coefficient(wd, 275.15, 1e5) < 0 < ch01.thermal_expansion_coefficient(wd, 281.15, 1e5)
    T = np.linspace(273.15, 283.15, 2001)
    assert abs(T[np.argmax(ch01.water_density(T))] - 277.13) < 0.05  # maximum near 3.98 °C
    assert 998.1 < ch01.water_density(293.15) < 998.3


# ====================================================================================================================
# C40 · Perfect-gas law (§1.9) — A item; with R03, C38, C39, C41, C48, C62, C63, D11
# ====================================================================================================================
@needs_ref
def test_perfect_gas_V5_codata_constants_and_ussa_sea_level():  # V5 (CODATA; USSA-1976)
    c = ref_json("constants.json")
    assert ch01.K_B == c["k_B_J_per_K"] and ch01.N_A == c["N_A_per_mol"]
    assert rel(ch01.R_U / 1e3, c["R_J_per_mol_K"]) < 1e-10
    u = ref_json("ussa1976_constants.json")
    assert ch01.G0 == u["g0_m_per_s2"] and ch01.P_ATM == u["P0_Pa"] and ch01.GAMMA_AIR == u["gamma"]
    assert rel(ch01.R_AIR, u["R_star_J_per_kmol_K"] / ch01.M_W_AIR) < 5e-5  # CODATA R_u vs USSA R* (1.7e-5)
    t1 = ref_csv("ussa1976_table1.csv")
    assert rel(ch01.perfect_gas_density(101325.0, 288.15), t1[0]["rho_kg_m3"]) < 1e-4
    for row in t1:  # p/(ρT) of every tabulated state = R of air
        assert rel(row["p_Pa"] / (row["rho_kg_m3"] * row["T_K"]), ch01.R_AIR) < 2e-4, row
    for row in ref_csv("ussa1976_table2.csv"):
        t = [r for r in t1 if r["Z_km"] == row["Z_km"]][0]
        assert rel(ch01.number_density(t["rho_kg_m3"], ch01.M_W_AIR), row["n_per_m3"]) < 2e-4, row
        assert rel(ch01.mean_molecular_speed(t["T_K"], M_AIR), row["V_m_s"]) < 1e-4, row
        assert rel(ch01.scale_height(t["T_K"], g=t["g_m_s2"]), row["Hp_m"]) < 2e-4, row


def test_perfect_gas_V1_constants_chain_D11():  # V1 (D11: p = n k_B T/V → ρ R T)
    assert rel(ch01.R_U, ch01.K_B * ch01.N_A_KMOL) < 1e-15
    assert rel(ch01.K_B / ch01.molecule_mass(ch01.M_W_AIR), ch01.R_AIR) < 1e-12  # step 4: k_B/m = k_B A_o/M_w
    assert rel(ch01.gas_constant(ch01.M_W_AIR), ch01.R_AIR) < 1e-15
    n_over_V = 2.547e25
    T = 288.15
    assert rel(ch01.molecular_gas_pressure(n_over_V, 1.0, T),
               ch01.perfect_gas_pressure(n_over_V * M_AIR, T)) < 1e-12  # steps 1–6
    assert rel(ch01.perfect_gas_density(101325.0, 288.15), 101325.0 / (287.058 * 288.15)) < 1e-5  # worked number
    assert 286.9 < ch01.R_U / 28.96 < 287.2  # R_air = 8314/28.96 ≈ 287
    assert rel(ch01.number_density(1.225, ch01.M_W_AIR), 101325.0 / (ch01.K_B * 288.15)) < 1e-4


def test_perfect_gas_V2_dimensions_and_kmol_trap():  # V2
    p = Q_(1.225, "kg/m**3") * Q_(ch01.R_AIR, "J/(kg*K)") * Q_(288.15, "K")
    assert p.check("[pressure]") and rel(ch01.perfect_gas_pressure(1.225, 288.15), p.to("Pa").magnitude) < 1e-14
    with pytest.warns(UserWarning):
        vec = ch01.dimension_vector("J/(kmol*K)")
    assert list(vec) == [1, 2, -2, -1]
    Ru_mol = Q_(8.314462618, "J/(mol*K)").to("J/(kmol*K)").magnitude
    assert rel(Ru_mol, ch01.R_U) < 1e-10  # the factor 1000 between mol and kmol


def test_isothermal_atmosphere_V1_scale_height():  # V1 (C62, C63) + V2; the book's own example is in the private book-value test
    z, p0, g, R, T = sp.symbols("z p0 g R T", positive=True)
    p = p0 * sp.exp(-g * z / (R * T))
    assert sp.simplify(sp.diff(p, z) + p * g / (R * T)) == 0
    T0 = 288.15  # USSA sea-level temperature (not the book's example temperature)
    H = ch01.scale_height(T0)
    assert rel(H, ch01.R_AIR * T0 / G) < 1e-14  # H = R T/g
    assert rel(ch01.isothermal_pressure(H, 1e5, T0), 1e5 / np.e) < 1e-14
    assert rel(ch01.isothermal_density(3000.0, 1.2, T0), 1.2 * ch01.isothermal_pressure(3000.0, 1.0, T0)) < 1e-14
    assert rel(ch01.scale_height(2 * T0), 2 * H) < 1e-14


# ====================================================================================================================
# C45 · Isentropic perfect gas (§1.9) — A item; with C42–C44, N16, C46, C47, D14
# ====================================================================================================================
def test_isentropic_law_V2_derivation_D14():  # V2 (D14 re-derived step by step)
    v, R, cv = sp.symbols("v R c_v", positive=True)
    C1 = sp.Symbol("C1", positive=True)
    cp = cv + R  # (1.23)
    gam = cp / cv  # (1.24)
    Tf = sp.Function("T")
    sol = sp.dsolve(sp.Eq(cv * Tf(v).diff(v), -R * Tf(v) / v), Tf(v))  # step 2: C_v dT = −p dv, p = RT/v
    T = C1 * v ** (-R / cv)
    assert sp.simplify(sol.rhs / v ** (-R / cv)).free_symbols <= {sp.Symbol("C1")}  # dsolve gives this family
    assert sp.simplify(cv * sp.diff(T, v) + R * T / v) == 0
    p = R * T / v
    assert sp.simplify(cp * sp.diff(T, v) - v * sp.diff(p, v)) == 0  # step 3: C_p dT = v dp holds on the same path
    assert sp.simplify(sp.diff(p, v) / p + gam / v) == 0  # step 6: dp/p = −γ dv/v
    rho = sp.Symbol("rho", positive=True)
    assert sp.simplify(sp.diff(1 / rho, rho) / (1 / rho) + 1 / rho) == 0  # step 7: dv/v = −dρ/ρ
    assert sp.simplify(sp.diff(p * v**gam, v)) == 0  # step 10: p/ρ^γ = p v^γ = const


def test_isentropic_law_V3_from_scratch_integration_converges():  # V3 (design §7 from-scratch dp/dρ = γp/ρ)
    rho0, p0, rho1, gam = 1.225, 101325.0, 2.45, 1.4
    exact = ch01.isentropic_pressure(rho1, p0, rho0, gam)
    hs, errs = [], []
    for n in (100, 200, 400, 800):
        rho = np.linspace(rho0, rho1, n + 1)
        p = p0
        for k in range(n):
            p += gam * p / rho[k] * (rho[k + 1] - rho[k])  # explicit Euler
        errs.append(abs(p - exact))
        hs.append(rho[1] - rho[0])
    assert abs(observed_order(hs, errs) - 1.0) < 0.15
    assert rel(exact, p0 * 2**1.4) < 1e-14


def test_isentropic_ratios_V1_consistency_and_limits():  # V1 + V7 (C46, C47)
    r = np.array([0.1, 0.5, 2.0, 8.0])
    Tr, rr = ch01.isentropic_ratios(r, 1.4)
    assert np.max(np.abs(Tr / rr**0.4 - 1)) < 1e-14 and np.max(np.abs(rr * Tr / r - 1)) < 1e-14  # p = ρRT
    T1, _ = ch01.isentropic_ratios(r, 1.0 + 1e-9)
    assert np.max(np.abs(T1 - 1)) < 1e-8  # γ → 1: isothermal
    Tp, rp = ch01.isentropic_ratios(2.0)
    assert 351.0 < 288.15 * Tp < 351.6 and abs(rp - 1.641) < 1e-3  # bicycle pump
    assert rel(ch01.isentropic_pressure(1.225 * rp, 1e5, 1.225), 2e5) < 1e-12
    with pytest.raises(ValueError):
        ch01.isentropic_ratios(0.0)


def test_specific_heat_relations_V1_round_trips():  # V1 (C42, C43, N16)
    assert rel(ch01.CP_AIR - ch01.CV_AIR, ch01.R_AIR) < 1e-13 and rel(ch01.CP_AIR / ch01.CV_AIR, 1.4) < 1e-14
    assert rel(ch01.gamma_from_cp(ch01.cp_from_gamma(5 / 3)), 5 / 3) < 1e-13
    assert rel(ch01.cv_from_cp(ch01.CP_AIR), ch01.CV_AIR) < 1e-13
    assert rel(ch01.cp_from_gamma(1.4) - ch01.cv_from_gamma(1.4), ch01.R_AIR) < 1e-13
    assert abs(ch01.gamma_from_cp(1005.0) - 1.40) < 0.002
    assert ch01.GAMMA_BY_ATOMICITY["diatomic"] == 1.4 and ch01.GAMMA_IDEAL is ch01.GAMMA_BY_ATOMICITY


@book_only
def test_perfect_gas_V6_book_constants():  # V6 (book rounds constants; ≤ 0.5 %)
    b = book()["perfect_gas"]
    ours = {"k_B_J_per_K": ch01.K_B, "A0_per_kmol": ch01.N_A_KMOL, "R_u_J_per_kmol_K": ch01.R_U,
            "M_w_dry_air_kg_per_kmol": ch01.M_W_AIR, "R_air_J_per_kg_K": ch01.R_AIR, "gamma_air": ch01.GAMMA_AIR,
            "cp_air_J_per_kg_K": ch01.CP_AIR}
    for k, v in ours.items():
        assert rel(v, b[k]["value"]) < 0.005, (k, v, b[k]["value"])
    u = book()["units"]
    assert u["celsius_offset_K"]["value"] == ch01.KELVIN_OFFSET
    assert rel(ch01.P_ATM / 1e3, u["p_atm_kPa"]["value"]) < 0.005 and rel(ch01.P_ATM / 1e5, u["p_atm_bar"]["value"]) < 0.005


# ====================================================================================================================
# C50 · Displaced-parcel equation (§1.10) — A item; D18
# ====================================================================================================================
def test_parcel_equation_V2_derivation_D18():  # V2 (D18 re-derived: steps 3–12)
    zeta, g, rho0, a, b, t = sp.symbols("zeta g rho0 a b t", real=True)
    rho_p = rho0 + a * zeta  # step 4 (a = dρ_a/dz)
    rho_e = rho0 + b * zeta  # step 5 (b = dρ/dz)
    acc = -g * (rho_p - rho_e) / rho_p  # step 3
    assert sp.simplify(rho_p - rho_e - (a - b) * zeta) == 0  # step 6
    lin = sp.series(acc, zeta, 0, 2).removeO()  # step 7: drop O(ζ²)
    N2 = -g / rho0 * (b - a)  # step 10, Eq. (1.29)
    assert sp.simplify(lin + N2 * zeta) == 0  # step 8–9: ζ'' + N²ζ = 0
    s = sp.Symbol("s", positive=True)
    z0 = sp.Symbol("zeta0", real=True)
    for sol, n2 in ((z0 * sp.cos(s * t), s**2), (z0 * sp.cosh(s * t), -s**2), (z0 + 0 * t, 0)):  # step 12
        assert sp.simplify(sp.diff(sol, t, 2) + n2 * sol) == 0 and sol.subs(t, 0) == z0
        assert sp.diff(sol, t).subs(t, 0) == 0
    f = sp.lambdify((g, rho0, b, a), N2)
    for args in ((9.81, 1025.0, -0.01, 0.0), (9.80665, 1.2, -1e-4, -1.1e-4)):
        assert rel(ch01.brunt_vaisala_sq(args[1], args[2], args[3], g=args[0]), f(*args)) < 1e-14


def test_parcel_equation_V1_analytic_solution_satisfies_ode():  # V1 (numerical 2nd derivative of the code)
    t = np.linspace(0, 3000, 30001)
    dt = t[1] - t[0]
    for N2, w0 in ((1e-4, 0.0), (-2e-5, 0.01), (0.0, 0.02), (3e-4, -0.05)):
        z = ch01.parcel_displacement(t, 5.0, N2, w0=w0)
        acc = (z[2:] - 2 * z[1:-1] + z[:-2]) / dt**2
        scale = max(abs(N2) * np.max(np.abs(z)), 1e-12)
        assert np.max(np.abs(acc + N2 * z[1:-1])) < 1e-5 * scale + 1e-9
        assert z[0] == pytest.approx(5.0) and (z[1] - z[0]) / dt == pytest.approx(w0, abs=1e-3)


def test_parcel_equation_V7_nonlinear_parcel_converges_to_linear():  # V7 (small-amplitude limit: error ∝ ζ0², ocean)
    t = np.linspace(0, 1500, 31)
    rho0, b, a = 1025.0, -0.01, -0.004
    N2 = ch01.brunt_vaisala_sq(rho0, b, a)
    z0s = [8.0, 4.0, 2.0, 1.0]
    errs = [np.max(np.abs(ch01.parcel_ode_from_gradients(t, z0, rho0, b, a) - ch01.parcel_displacement(t, z0, N2)))
            for z0 in z0s]
    p = observed_order(z0s, errs)
    assert abs(p - 2.0) < 0.15, (p, errs)
    ts, zs = ch01.parcel_ode((0.0, 1500.0), 4.0, lambda z: rho0 + b * z, lambda zeta: rho0 + a * zeta, t_eval=t)
    assert np.max(np.abs(zs - ch01.parcel_ode_from_gradients(t, 4.0, rho0, b, a))) < 1e-8
    assert ch01.parcel_ode_from_gradients(0.0, 3.0, rho0, b, a) == 3.0


def test_parcel_equation_V7_atmospheric_parcel_converges_to_linear():  # V7 (small-amplitude limit, dry parcel)
    t = np.linspace(0, 1500, 31)
    T0, G_ = 288.15, -6.5e-3
    N2 = ch01.brunt_vaisala_sq_from_lapse(T0, G_)
    z0s = [200.0, 100.0, 50.0, 25.0]
    errs = [np.max(np.abs(ch01.parcel_ode_atmosphere(t, z0, T0, G_) - ch01.parcel_displacement(t, z0, N2))) for z0 in z0s]
    assert abs(observed_order(z0s, errs) - 2.0) < 0.15, errs
    zn = ch01.parcel_ode_atmosphere(t, 100.0, T0, ch01.adiabatic_lapse_rate())  # neutral: stays put
    assert np.max(np.abs(zn - 100.0)) < 1e-8
    eps = 1e-3
    da = (ch01.parcel_acceleration_atmosphere(eps, T0, G_) - ch01.parcel_acceleration_atmosphere(-eps, T0, G_)) / (2 * eps)
    assert rel(-da, N2) < 1e-6 and ch01.parcel_acceleration_atmosphere(0.0, T0, G_) == 0.0


def _parcel_potential(zeta, g, rho0, drho_dz, drho_a_dz, denominator_gradient):
    """Exact potential V(ζ) with dV/dζ = −g (b − a) ζ/(ρ0 + c ζ) for linear profiles (b = dρ/dz, a = dρ_a/dz).

    c = a divides by the PARCEL density ρ_p = ρ0 + aζ (Newton II for the parcel's own mass); c = b would divide by the
    environment density. ∫ ζ/(ρ0 + cζ) dζ = ζ/c − (ρ0/c²) ln(1 + cζ/ρ0).
    """
    c = denominator_gradient
    return -g * (drho_dz - drho_a_dz) * (zeta / c - rho0 / c**2 * np.log1p(c * zeta / rho0))


def test_parcel_equation_V4_energy_first_integral_divides_by_parcel_density():  # V4 (first integral) + V1 (turning point, period)
    # Book §1.10 before linearising: (ρ_p V) ζ'' = −ρ_p V g + ρ_e V g ⇒ ζ'' = −g (ρ_p − ρ_e)/ρ_p. With linear profiles
    # ρ_p = ρ0 + aζ, ρ_e = ρ0 + bζ this is ζ'' = g (b − a) ζ/(ρ0 + aζ) = −dV/dζ, so E = ½ζ'² + V(ζ) is conserved exactly.
    # A strongly non-linear case (|a ζ0|/ρ0 = 0.2, |b ζ0|/ρ0 = 0.4) separates ÷ρ_p from ÷ρ_e by ~20 m in the turning
    # point, while both share the same N² and small-amplitude limit (so the V7 test cannot tell them apart).
    from scipy.integrate import quad
    from scipy.optimize import brentq

    g, rho0, b, a, z0 = G, 1.0, -0.01, 0.005, 40.0
    env = lambda z: rho0 + b * z  # noqa: E731
    par = lambda zeta: rho0 + a * zeta  # noqa: E731
    Vp = lambda z: _parcel_potential(z, g, rho0, b, a, a)  # noqa: E731  (÷ parcel density: the correct law)
    Ve = lambda z: _parcel_potential(z, g, rho0, b, a, b)  # noqa: E731  (÷ environment density: the wrong law)
    z1 = brentq(lambda z: Vp(z) - Vp(z0), -99.0, -1e-9)  # lower turning point, V(z1) = V(z0)
    z1_env = brentq(lambda z: Ve(z) - Ve(z0), -99.0, -1e-9)
    assert abs(z1 - z1_env) > 10.0  # the test discriminates the two laws by metres, not by round-off
    c, d = 0.5 * (z0 + z1), 0.5 * (z0 - z1)  # period by quadrature, ζ = c + d sin θ removes the endpoint singularities
    P, P_err = quad(lambda th: 2 * d * np.cos(th) / np.sqrt(2 * (Vp(z0) - Vp(c + d * np.sin(th)))), -np.pi / 2,
                    np.pi / 2, limit=200)
    assert P_err < 1e-10 * P
    t = np.linspace(0.0, 1.25 * P, 25001)
    dt = t[1] - t[0]
    ts, zs = ch01.parcel_ode((0.0, t[-1]), z0, env, par, t_eval=t, zeta_max=1e6)
    assert ts.size == t.size
    w = (zs[2:] - zs[:-2]) / (2 * dt)  # central difference, O(dt²)
    E = 0.5 * w**2 + Vp(zs[1:-1])
    assert np.ptp(E) < 1e-6 * abs(Vp(z0)), np.ptp(E) / abs(Vp(z0))  # V4: first integral conserved
    E_env = 0.5 * w**2 + Ve(zs[1:-1])
    assert np.ptp(E_env) > 1e-2 * abs(Ve(z0))  # the ÷ρ_e invariant is visibly not conserved by this trajectory

    def extremum(i):  # parabola through three samples → vertex (time, value)
        y0, y1, y2 = zs[i - 1], zs[i], zs[i + 1]
        s = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
        return t[i] + s * dt, y1 - 0.25 * (y0 - y2) * s

    i_min = int(np.argmin(zs))
    assert abs(extremum(i_min)[1] - z1) < 1e-6 * abs(z1), (extremum(i_min)[1], z1, z1_env)  # V1 turning point
    late = t > 0.75 * P
    i_max = int(np.argmax(np.where(late, zs, -np.inf)))
    t_ret, z_ret = extremum(i_max)
    assert abs(t_ret - P) < 1e-6 * P and abs(z_ret - z0) < 1e-6 * z0, (t_ret, P, z_ret)  # V1 period, return to ζ0
    assert abs(extremum(i_min)[0] - P / 2) < 1e-6 * P  # asymmetric potential, but half-period to the lower turn
    # the gradient wrapper runs the same law
    assert np.max(np.abs(ch01.parcel_ode_from_gradients(t[::500], z0, rho0, b, a) - zs[::500])) < 1e-6 * z0


def test_parcel_ode_V7_runaway_capped_with_nan():  # V7 (N² < 0: NaN after the cap, never overflow)
    t = np.linspace(0, 5000, 501)
    z = ch01.parcel_ode_from_gradients(t, 1.0, 1000.0, 0.01, 0.0)  # denser above: unstable
    fin = np.isfinite(z)
    assert np.any(~fin) and np.all(np.isnan(z[~fin])) and not np.any(np.isinf(z))
    assert np.all(fin[: np.argmin(fin)]) and not np.any(fin[np.argmin(fin):])  # finite then NaN
    assert np.nanmax(np.abs(z)) <= 1000.0 * 1.0 + 1e-6
    z2 = ch01.parcel_ode_from_gradients(t, 1.0, 1000.0, 0.01, 0.0, zeta_max=50.0)
    assert np.nanmax(np.abs(z2)) <= 50.0 + 1e-9 and np.isnan(z2[-1])
    za = ch01.parcel_ode_atmosphere(t, 1.0, 288.15, -0.012)
    assert np.isnan(za[-1]) and not np.any(np.isinf(za)) and np.nanmax(np.abs(za)) <= 1000.0 + 1e-6
    zs = ch01.parcel_ode_atmosphere(t, 1.0, 288.15, -0.012, zeta_max=1e7)  # cap off: the too-cold guard stops it
    assert np.isnan(zs[-1]) and not np.any(np.isinf(zs))
    early = t < 400
    lin = ch01.parcel_displacement(t[early], 1.0, ch01.brunt_vaisala_sq(1000.0, 0.01, 0.0))
    assert rel(z[early], lin) < 1e-3  # before the cap the run follows cosh


# ====================================================================================================================
# C51 · Brunt–Väisälä frequency (§1.10) — A item; with C52, C60, N21, C61, D36
# ====================================================================================================================
def test_brunt_vaisala_V1_thermocline_worked_number_and_timescales():  # V1
    N2 = ch01.brunt_vaisala_sq(1025.0, -0.01, 0.0, g=9.81)
    assert rel(N2, 9.81 * 0.01 / 1025.0) < 1e-14 and 9.5e-5 < N2 < 9.6e-5
    kind, per = ch01.stability_timescale(N2)
    assert kind == "period" and 640 < per < 645
    t = np.linspace(0, per, 200001)
    z = ch01.parcel_displacement(t, 1.0, N2)
    first_zero = t[np.argmax(z <= 0)]
    assert abs(first_zero - per / 4) < 2 * (t[1] - t[0])
    assert ch01.stability_timescale(-1e-4) == ("efold", pytest.approx(100.0))
    assert ch01.stability_timescale(0.0) == ("none", float("inf"))
    assert ch01.brunt_vaisala_sq(1000.0, 0.0, 0.0) == 0.0  # uniform incompressible: neutral
    assert ch01.brunt_vaisala_sq(1.2, -1e-4, -1e-4) == 0.0  # isentropic: neutral
    with pytest.raises(ValueError):
        ch01.brunt_vaisala_sq(0.0, 1.0, 1.0)


def test_brunt_vaisala_V2_isothermal_atmosphere_symbolic():  # V2 (N² = g²/(C_p T))
    z, p0, g, R, T, gam = sp.symbols("z p0 g R T gamma", positive=True)
    p = p0 * sp.exp(-g * z / (R * T))
    rho = p / (R * T)
    drho_a = -rho * g / (gam * R * T)  # −ρg/c², c² = γRT
    N2 = -g / rho * (sp.diff(rho, z) - drho_a)
    cp = gam * R / (gam - 1)
    assert sp.simplify(N2 - g**2 / (cp * T)) == 0
    val = ch01.brunt_vaisala_sq_from_lapse(250.0, 0.0)
    assert rel(val, G**2 / (ch01.CP_AIR * 250.0)) < 1e-14 and 3.8e-4 < val < 3.9e-4


def test_brunt_vaisala_V1_density_form_equals_lapse_and_theta_forms():  # V1 (consistency: three independent routes to N²)
    z = np.linspace(0.0, 6000.0, 6001)
    Tfun = lambda zz: 288.0 - 6e-3 * np.asarray(zz) + 3.0 * np.sin(np.asarray(zz) / 700.0)  # noqa: E731
    dTdz = -6e-3 + 3.0 / 700.0 * np.cos(z / 700.0)
    p, rho, T = ch01.atmosphere_from_temperature(z, Tfun, 101325.0)
    drho = np.gradient(rho, z, edge_order=2)
    c = ch01.perfect_gas_sound_speed(T)
    N2_rho = ch01.brunt_vaisala_sq(rho, drho, ch01.isentropic_density_gradient(rho, c))  # (1.29)
    N2_lapse = ch01.brunt_vaisala_sq_from_lapse(T, dTdz)
    theta = ch01.potential_temperature(T, p)
    N2_theta = ch01.brunt_vaisala_sq_from_theta(theta, np.gradient(theta, z, edge_order=2))
    sl = slice(5, -5)
    scale = np.max(np.abs(N2_lapse))
    assert np.max(np.abs(N2_rho[sl] - N2_lapse[sl])) < 1e-4 * scale
    assert np.max(np.abs(N2_theta[sl] - N2_lapse[sl])) < 1e-4 * scale
    assert np.array_equal(ch01.classify_stability(N2_lapse), ch01.classify_stability(N2_rho))


def test_stability_classification_V7_labels_and_tolerance():  # V7 (C52)
    assert list(ch01.classify_stability(np.array([1e-4, 0.0, -1e-4, 5e-13, -5e-13]))) == \
        ["stable", "neutral", "unstable", "neutral", "neutral"]
    assert ch01.classify_stability(2e-12) == "stable" and ch01.classify_stability(-1e-3, tol=1e-2) == "neutral"


def test_ocean_criterion_V1_sign_and_limits():  # V1 + V7 (N21, C61)
    rho, c = 1025.0, 1500.0
    assert rel(ch01.isentropic_density_gradient(rho, c), -rho * G / c**2) < 1e-14
    assert 4.4e-3 < rho * G / c**2 < 4.6e-3  # worked number
    rng = np.random.default_rng(11)
    for drho in rng.uniform(-0.02, 0.02, 200):
        N2 = ch01.brunt_vaisala_sq(rho, drho, ch01.isentropic_density_gradient(rho, c))
        ind = ch01.ocean_potential_density_gradient(drho, rho, c)
        assert np.sign(ind) == -np.sign(N2)
    assert rel(ch01.ocean_potential_density_gradient(-0.01, rho, 1e12), -0.01) < 1e-12  # c → ∞
    # perfect gas: −ρg/c² equals the isentropic parcel gradient from (1.26) and hydrostatics
    T, p = 260.0, 7e4
    rho_g = ch01.perfect_gas_density(p, T)
    drho_a = rho_g / (ch01.GAMMA_AIR * p) * (-rho_g * G)  # dρ_a/dz = (ρ/γp) dp/dz
    assert rel(ch01.isentropic_density_gradient(rho_g, ch01.perfect_gas_sound_speed(T)), drho_a) < 1e-12


def test_seawater_eos_V1_linear_fit_consistent_with_eos80():  # V1 (C60; self-consistency only, see SOURCES.md)
    rho0, aT, bS = ch01.seawater_linear_coefficients()
    assert abs(rho0 - 1027.0) < 0.1 and abs(aT - 1.67e-4) < 2e-6 and abs(bS - 7.6e-4) < 2e-6
    assert ch01.seawater_density_linear(283.15, 35.0) == 1027.0
    for T in np.linspace(278.15, 288.15, 5):
        for S in np.linspace(33, 37, 5):
            assert abs(ch01.seawater_density_linear(T, S) - ch01.seawater_density_eos80(T, S)) < 0.5
    assert ch01.seawater_density_eos80(283.15, 35.0) > ch01.seawater_density_eos80(283.15, 34.0)


@needs_ref
def test_seawater_eos_V5_unesco_eos80_check_values():  # V5 (Fofonoff & Millard 1983, Unesco Tech. Pap. Mar. Sci. 44, p. 19)
    chk = ref_json("benchmarks.json")["unesco_eos80_check"]
    rows = [r for r in chk["rows"] if r["p_dbar"] == 0.0]  # the one-atmosphere equation (p = 0 dbar gauge)
    assert len(rows) == 4
    for r in rows:
        T90_K = r["t68_C"] / 1.00024 + ch01.KELVIN_OFFSET  # the function takes ITS-90 kelvin and converts to t68
        rho = ch01.seawater_density_eos80(T90_K, r["S"])
        assert abs(rho - r["rho_kg_m3"]) <= 5e-6, (r, rho)  # half a unit of the 5th printed decimal
        assert rel(1.0 / rho, r["V_1e-3_m3_kg"] * 1e-3) < 1e-8, (r, 1.0 / rho)  # specific-volume column [m^3/kg]
        wrong_scale = ch01.seawater_density_eos80(r["t68_C"] + ch01.KELVIN_OFFSET, r["S"])  # t68 fed in as if ITS-90
        assert abs(wrong_scale - r["rho_kg_m3"]) > 5e-6  # the table resolves the temperature-scale conversion


# ====================================================================================================================
# C54 · Adiabatic lapse rate (§1.10) and the two lapse-rate conventions — A item; with C49, C53, N17, D19
# ====================================================================================================================
def test_adiabatic_lapse_rate_V2_derivation_D19():  # V2 (D19 ★★★: any fluid, then a perfect gas)
    T, p, g, R, cp = sp.symbols("T p g R c_p", positive=True)
    Gf = sp.Function("G")(T, p)  # any smooth Gibbs free energy per unit mass
    v = sp.diff(Gf, p)  # step 6
    s = -sp.diff(Gf, T)  # step 6
    h = Gf + T * s  # g = h − Ts
    Cp = sp.diff(h, T)  # (1.14)
    alpha = sp.diff(v, T) / v  # step 8: (∂v/∂T)_p = vα
    assert sp.simplify(sp.diff(v, T) + sp.diff(s, p)) == 0  # step 7: Maxwell relation
    assert sp.simplify(sp.diff(h, p) - (T * sp.diff(s, p) + v)) == 0  # step 4
    assert sp.simplify(sp.diff(h, p) - v * (1 - alpha * T)) == 0  # step 10
    dTdp_s = -sp.diff(s, p) / sp.diff(s, T)  # isentrope slope
    assert sp.simplify(dTdp_s - v * alpha * T / Cp) == 0  # step 12: C_p dT = vαT dp
    Gamma_a = dTdp_s * (-g / v)  # steps 13–14: dp/dz = −ρg = −g/v
    assert sp.simplify(Gamma_a + g * alpha * T / Cp) == 0  # step 15: (1.30)
    Gpg = cp * (T - T * sp.log(T)) + R * T * sp.log(p)  # perfect gas, constant C_p
    subs = lambda e: sp.simplify(e.subs(Gf, Gpg).doit())  # noqa: E731
    assert sp.simplify(subs(alpha) - 1 / T) == 0 and sp.simplify(subs(Cp) - cp) == 0
    assert sp.simplify(subs(Gamma_a) + g / cp) == 0  # step 16: −g/C_p
    alpha_s, T_s, cp_s = 1.5e-4, 283.15, 4190.0  # water number of the check
    assert rel(ch01.adiabatic_lapse_rate(T=T_s, cp=cp_s, alpha=alpha_s), -G * alpha_s * T_s / cp_s) < 1e-14


def test_adiabatic_lapse_rate_V1_perfect_gas_sign_and_general_form():  # V1
    Ga = ch01.adiabatic_lapse_rate()
    assert Ga < 0 and rel(Ga, -G / ch01.CP_AIR) < 1e-14 and abs(Ga * 1e3 + 9.7607) < 1e-4
    T = np.array([220.0, 288.15, 310.0])
    assert rel(ch01.adiabatic_lapse_rate(T=T, alpha=1 / T), -G / ch01.CP_AIR) < 1e-14
    assert ch01.adiabatic_lapse_rate(T=300.0, alpha=0.0) == 0.0  # α = 0: no adiabatic temperature change
    assert -0.11e-3 < ch01.adiabatic_lapse_rate(T=283.15, cp=4190.0, alpha=1.5e-4) < -0.09e-3  # water
    with pytest.raises(ValueError):
        ch01.adiabatic_lapse_rate(alpha=1e-3)
    with pytest.warns(UserWarning, match="alpha"):  # T without alpha: perfect gas assumed, T ignored — say so
        Gw = ch01.adiabatic_lapse_rate(T=283.15)
    assert Gw == Ga
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # no warning on the normal calls
        ch01.adiabatic_lapse_rate()
        ch01.adiabatic_lapse_rate(T=283.15, cp=4190.0, alpha=1.5e-4)
    q = Q_(G, "m/s**2") * Q_(1 / 288.15, "1/K") * Q_(288.15, "K") / Q_(ch01.CP_AIR, "J/(kg*K)")
    assert q.check("[temperature]/[length]")


def test_adiabatic_lapse_rate_V1_from_scratch_lifted_parcel_slope_at_release():  # V1 (design §7 from-scratch)
    T0, Gam = 288.15, -6.5e-3
    z = np.linspace(0.0, 2000.0, 2001)
    pe = ch01.linear_lapse_pressure(z, 101325.0, T0, Gam)  # hydrostatic environment
    Te = T0 + Gam * z
    i0 = 1000  # release height 1000 m
    Tp = Te[i0] * (pe / pe[i0]) ** ((ch01.GAMMA_AIR - 1) / ch01.GAMMA_AIR)  # parcel on its isentrope, (1.26)
    slope = np.gradient(Tp, z, edge_order=2)
    assert rel(slope[i0], ch01.adiabatic_lapse_rate()) < 1e-6  # −g/C_p exactly at the release height
    local = ch01.adiabatic_lapse_rate() * Tp / Te  # elsewhere −(g/C_p) T_p/T_e
    assert np.max(np.abs(slope[5:-5] / local[5:-5] - 1)) < 1e-6
    assert abs(slope[100] / ch01.adiabatic_lapse_rate() - 1) > 5e-3  # …which is not −g/C_p away from release
    met = -slope[i0]
    assert rel(met, ch01.lapse_rate_convention(ch01.adiabatic_lapse_rate(), "meteorology")) < 1e-6


@needs_ref
def test_adiabatic_lapse_rate_V5_ams_glossary_value():  # V5 (AMS: g/c_pd ≈ 9.8 °C/km, rate of decrease)
    ams = ref_json("benchmarks.json")["ams_lapse_rate"]
    met_per_km = ch01.lapse_rate_convention(ch01.adiabatic_lapse_rate(), "meteorology") * 1e3
    assert met_per_km > 0 and round(met_per_km, 1) == ams["dry_adiabatic_lapse_rate_C_per_km"]
    assert rel(met_per_km, ams["dry_adiabatic_lapse_rate_C_per_km"]) < 0.01


@book_only
def test_adiabatic_lapse_rate_V6_book_rounded_value():  # V6 (book prints a 1-significant-figure value)
    s = book()["stratification"]
    b = s["adiabatic_lapse_rate_C_per_km"]["value"]
    ours = ch01.adiabatic_lapse_rate() * 1e3
    assert np.sign(ours) == np.sign(b) and round(ours) == b  # same sign convention (Γ ≡ dT/dz), same 1-s.f. value
    lab = s["lab_scale_temperature_change"]
    dT = ch01.adiabatic_lapse_rate() * lab["dz_m"]
    assert float(f"{dT:.0e}") == lab["dT_C"]  # lab scale: T and θ differ negligibly
    assert ch01.P_REF == s["standard_pressure_p0_kPa_approx"]["value"] * 1e3


def test_lapse_rate_convention_V1_sign_flip_involution_and_errors():  # V1 (C53, N17)
    Ga = ch01.adiabatic_lapse_rate()
    assert ch01.lapse_rate_convention(Ga, "meteorology") == pytest.approx(G / ch01.CP_AIR, rel=1e-14)
    assert ch01.lapse_rate_convention(Ga, "kundu") == Ga and ch01.lapse_rate_convention(Ga) == Ga
    x = np.linspace(-15e-3, 10e-3, 26)
    assert np.array_equal(ch01.lapse_rate_convention(ch01.lapse_rate_convention(x, "meteorology"), "meteorology"), x)
    assert np.array_equal(ch01.lapse_rate_convention(x, "  Meteorology "), -x)
    for bad in ("foo", "", "met", "Kundu's", None):
        with pytest.raises(ValueError):
            ch01.lapse_rate_convention(1e-3, bad)
    with pytest.raises(ValueError):
        ch01.lapse_rate_stability(-6.5e-3, convention="meteo")
    assert ch01.LAPSE_RATE_CONVENTIONS == ("kundu", "meteorology")


def test_lapse_rate_stability_V1_exact_texts_standard_troposphere():  # V1 (exact strings, both conventions)
    k = ch01.lapse_rate_stability(-6.5e-3, convention="kundu")
    m = ch01.lapse_rate_stability(-6.5e-3, convention="meteorology")
    assert (k.verdict, k.text, k.code) == ("stable", "stable ⇔ dT/dz > Γa: −6.5 > −9.8 K/km", 1)
    assert (m.verdict, m.text, m.code) == ("stable", "stable ⇔ Γ < Γa: 6.5 < 9.8 K/km", 1)
    assert k[:2] == ("stable", "stable ⇔ dT/dz > Γa: −6.5 > −9.8 K/km") and isinstance(k, ch01.LapseStability)
    assert k.Gamma == -6.5e-3 and m.Gamma == 6.5e-3 and k.Gamma_a < 0 < m.Gamma_a and k.margin == m.margin
    assert ch01.lapse_rate_stability(-12e-3).text == "stable ⇔ dT/dz > Γa: −12.0 < −9.8 K/km"
    assert ch01.lapse_rate_stability(-12e-3, convention="meteorology").text == "stable ⇔ Γ < Γa: 12.0 > 9.8 K/km"
    assert ch01.lapse_rate_stability(-6.5e-3, prefix=False).text == "−6.5 > −9.8 K/km"
    assert ch01.lapse_rate_stability(-6.5e-3, prefix=False, ascii_minus=True).text == "-6.5 > -9.8 K/km"
    assert ch01.lapse_rate_stability(-6.5e-3, prefix=False, per_km=False).text == "−0.00650 > −0.00976 K/m"
    neutral = ch01.lapse_rate_stability(ch01.adiabatic_lapse_rate())
    assert neutral.code == 0 and neutral.verdict == "neutral" and " = " in neutral.text
    assert ch01.lapse_rate_stability(-1e-6, prefix=False).text == "0.0 > −9.8 K/km"  # never prints "−0.0"
    assert ch01.lapse_rate_stability(0.0, prefix=False).text == "0.0 > −9.8 K/km"
    with pytest.raises(ValueError):
        ch01.lapse_rate_stability(np.array([-6.5e-3, -7e-3]))


def test_lapse_rate_stability_V1_auto_decimals_near_adiabat():  # V1 (auto-decimals rule)
    assert ch01.lapse_rate_stability(-9.8e-3).text == "stable ⇔ dT/dz > Γa: −9.80 < −9.76 K/km"
    assert ch01.lapse_rate_stability(-9.8e-3, convention="meteorology").text == "stable ⇔ Γ < Γa: 9.80 > 9.76 K/km"
    assert ch01.lapse_rate_stability(-9.76e-3, prefix=False).text == "−9.760 > −9.761 K/km"
    Ga = ch01.adiabatic_lapse_rate()
    for conv in ("kundu", "meteorology"):
        for dm in (1e-8, -1e-8, 3e-7, -3e-7, 4e-5, -4e-5, 2e-3, -2e-3):  # margins ≥ 10 tol
            r = ch01.lapse_rate_stability(Ga + dm, convention=conv, prefix=False, ascii_minus=True)
            a_txt, op, b_txt = r.text.split()[:3]
            assert a_txt != b_txt, r.text  # two different numbers never print alike
            a, b = float(a_txt), float(b_txt)
            assert {">": a > b, "<": a < b}[op], r.text  # the printed relation is true for the printed numbers
            assert len(a_txt.split(".")[1]) <= 6


def test_lapse_rate_stability_V7_same_code_both_conventions_and_sign_of_N2():  # V7 (invariance under the sign convention; sweep −15…+10 K/km)
    Ga = ch01.adiabatic_lapse_rate()
    sweep = np.concatenate([np.linspace(-15e-3, 10e-3, 2501), [Ga, Ga + 1e-7, Ga - 1e-7]])
    to_code = {"stable": 1, "neutral": 0, "unstable": -1}
    for x in sweep:
        k = ch01.lapse_rate_stability(float(x), convention="kundu")
        m = ch01.lapse_rate_stability(float(x), convention="meteorology")
        assert k.code == m.code and k.verdict == m.verdict and k.margin == m.margin
        assert k.Gamma == pytest.approx(-m.Gamma) and k.Gamma_a == pytest.approx(-m.Gamma_a)
        N2 = ch01.brunt_vaisala_sq_from_lapse(260.0, float(x))
        if abs(k.margin) > 1e-9 or k.margin == 0.0:
            assert to_code[ch01.classify_stability(N2, tol=0.0 if k.margin == 0 else 1e-12)] == k.code, x
        assert np.sign(ch01.potential_temperature_gradient(260.0, float(x), theta=280.0)) == k.code or k.code == 0
    custom = ch01.lapse_rate_stability(-3e-3, Gamma_a=-2e-3)  # a moist-like reference gradient
    assert custom.code == -1 and custom.Gamma_a == -2e-3


def test_lapse_rate_V3_profile_derivative_second_order():  # V3 + V1 (lapse_rate, parcel_temperature)
    hs, errs = [], []
    for N in (21, 41, 81, 161):
        z = np.sort(np.linspace(0, 3000, N) + np.r_[0, np.sin(np.arange(1, N - 1)) * 3000 / (N * 8), 0])
        T = 288.0 + 5 * np.cos(z / 600)
        errs.append(np.max(np.abs(ch01.lapse_rate(T, z) + 5 / 600 * np.sin(z / 600))))
        hs.append(3000 / (N - 1))
    assert abs(observed_order(hs, errs) - 2.0) < 0.25, pairwise_orders(hs, errs)
    z = np.linspace(0, 1000, 11)
    assert np.allclose(ch01.lapse_rate(290 - 6.5e-3 * z, z), -6.5e-3, rtol=1e-12)
    Tp = ch01.parcel_temperature(290.0, z, z0=200.0)
    assert np.allclose(np.diff(Tp) / np.diff(z), ch01.adiabatic_lapse_rate(), rtol=1e-12) and Tp[2] == 290.0


# ====================================================================================================================
# C55 · Potential temperature (§1.10) — A item; with C56–C59, N18, D20, D36
# ====================================================================================================================
def test_potential_temperature_V1_definition_and_worked_number():  # V1 (D20)
    assert ch01.potential_temperature(273.0, ch01.P_REF) == 273.0
    th = ch01.potential_temperature(250.0, 5e4)
    assert rel(th, 250.0 * 2 ** (2 / 7)) < 1e-12 and 304.7 < th < 304.9
    T = np.array([220.0, 260.0, 300.0])
    p = np.array([2.5e4, 6e4, 1.02e5])
    assert rel(ch01.temperature_from_potential(ch01.potential_temperature(T, p), p), T) < 1e-14
    with pytest.raises(ValueError):
        ch01.potential_temperature(250.0, 0.0)


def test_potential_temperature_V2_derivation_D20_and_eq_1_32():  # V2 (D20 steps 4–5; (1.32) residual)
    gam, R, cv, T, p, po = sp.symbols("gamma R c_v T p p_o", positive=True)
    cp = cv + R
    assert sp.simplify(((cp / cv) - 1) / (cp / cv) - R / cp) == 0  # (γ−1)/γ = R/C_p
    z, g = sp.symbols("z g", positive=True)
    Tz = sp.Function("T")(z)
    pz = sp.Function("p")(z)
    kap = R / cp
    theta = Tz * (po / pz) ** kap
    lhs = Tz / theta * sp.diff(theta, z)
    rhs = sp.diff(Tz, z) + g / cp
    hydro = {sp.Derivative(pz, z): -pz * g / (R * Tz)}  # (1.8) with ρ = p/RT
    assert sp.simplify(lhs.subs(hydro) - rhs) == 0  # (1.32)


def test_potential_temperature_V4_constant_along_dry_adiabat():  # V4 (invariant along an independently built column)
    z = np.linspace(0, 10000, 41)
    p, rho, T = ch01.atmosphere_from_temperature(z, lambda zz: 300.0 - G / ch01.CP_AIR * np.asarray(zz))
    th = ch01.potential_temperature(T, p)
    assert np.ptp(th) / th[0] < 1e-8
    N2 = ch01.brunt_vaisala_sq_from_theta(th, np.gradient(th, z))
    assert np.max(np.abs(N2)) < 1e-10


def test_potential_temperature_V2_derivation_D36():  # V2 (D36 ★★: N² = (g/θ)dθ/dz = (g/T)(Γ − Γ_a))
    z, g, R, gam, po = sp.symbols("z g R gamma p_o", positive=True)
    Tz = sp.Function("T")(z)
    pz = sp.Function("p")(z)
    rho = pz / (R * Tz)
    dlnrho = sp.diff(rho, z) / rho
    assert sp.simplify(dlnrho - (sp.diff(pz, z) / pz - sp.diff(Tz, z) / Tz)) == 0  # step 2
    dlnrho_a = sp.diff(pz, z) / (gam * pz)  # step 3
    theta = Tz * (po / pz) ** ((gam - 1) / gam)
    dlntheta = sp.diff(theta, z) / theta
    assert sp.simplify(dlntheta - (sp.diff(Tz, z) / Tz - (gam - 1) / gam * sp.diff(pz, z) / pz)) == 0  # step 5
    N2_rho = -g * (dlnrho - dlnrho_a)  # step 1 (ρ_a = ρ at rest height)
    assert sp.simplify(N2_rho - g * dlntheta) == 0  # steps 6–7
    cp = gam * R / (gam - 1)
    hydro = {sp.Derivative(pz, z): -rho * g}
    assert sp.simplify((g * dlntheta).subs(hydro) - g / Tz * (sp.diff(Tz, z) + g / cp)) == 0  # step 8


def test_potential_temperature_gradient_V1_matches_finite_differences():  # V1 (C56; identity θρ_θ = p_ref/R)
    z = np.linspace(0, 11000, 11001)
    T, p, rho = ch01.standard_atmosphere(z)
    th = ch01.potential_temperature(T, p)
    dth = ch01.potential_temperature_gradient(T, np.full_like(z, -6.5e-3), p=p)
    assert np.max(np.abs(dth[2:-2] / np.gradient(th, z)[2:-2] - 1)) < 1e-6
    assert rel(ch01.potential_temperature_gradient(T[5], -6.5e-3, theta=th[5]), dth[5]) < 1e-14
    with pytest.raises(ValueError):
        ch01.potential_temperature_gradient(280.0, -6.5e-3)
    rth = ch01.potential_density(rho, p)
    assert np.max(np.abs(th * rth / (ch01.P_REF / ch01.R_AIR) - 1)) < 1e-12
    lhs = -np.gradient(rth, z)[2:-2] / rth[2:-2]  # (1.34)
    rhs = np.gradient(th, z)[2:-2] / th[2:-2]
    assert np.max(np.abs(lhs - rhs)) < 1e-6 * np.max(np.abs(rhs))
    assert ch01.potential_density(1.2, ch01.P_REF) == 1.2


def test_synthetic_column_V1_layers_consistent():  # V1 (N18; default mixed layer = dry adiabat, Fig. 1.9 neutral)
    z = np.linspace(0, 2000, 201)
    T = ch01.synthetic_boundary_layer_profile(z)
    col = ch01.synthetic_boundary_layer_column(z)
    assert rel(col["T"], T) < 1e-14
    Ga = ch01.adiabatic_lapse_rate()
    assert ch01.MIXED_LAYER_DT_DZ == Ga and rel(Ga, -G / ch01.CP_AIR) < 1e-14  # default slope is −g/C_p, not −9.8e-3
    mids = {"mixed": 400.0, "inversion": 900.0, "upper": 1500.0}
    slopes = {"mixed": ch01.MIXED_LAYER_DT_DZ, "inversion": 0.01, "upper": -0.0045}
    for k, zz in mids.items():
        i = int(np.argmin(np.abs(z - zz)))
        assert col["dT_dz"][i] == slopes[k]
        assert abs(col["N2"][i] - ch01.brunt_vaisala_sq_from_lapse(col["T"][i], slopes[k])) \
            <= 1e-10 * abs(ch01.brunt_vaisala_sq_from_lapse(col["T"][i], slopes[k])) + 1e-15
    assert abs(ch01.synthetic_boundary_layer_profile(800.0 - 1e-9) - ch01.synthetic_boundary_layer_profile(800.0 + 1e-9)) < 1e-8
    below = z < 790
    assert np.all(col["stability"][below] == "neutral"), set(col["stability"][below])  # mixed layer neutral …
    assert np.max(np.abs(col["N2"][below])) < 1e-15  # … because N² = 0 (scale of N² here ~1e-4 s⁻²)
    assert np.ptp(col["theta"][below]) < 1e-9  # θ constant along the dry adiabat
    assert np.all(col["stability"][(z > 810) & (z < 990)] == "stable") and np.all(col["stability"][z > 1010] == "stable")
    steep = ch01.synthetic_boundary_layer_column(z, mixed_dT_dz=-9.8e-3)  # the old default: steeper than Γa → unstable
    assert np.all(steep["stability"][below] == "unstable")
    assert np.max(np.abs(col["theta"] * col["rho_theta"] / (ch01.P_REF / ch01.R_AIR) - 1)) < 1e-12
    p_int = ch01.integrate_hydrostatic(z, lambda zz, pp: pp / (ch01.R_AIR * float(ch01.synthetic_boundary_layer_profile(zz))),
                                       ch01.P_ATM, rtol=1e-11)
    assert rel(col["p"], p_int) < 1e-7  # closed-form layers = direct integration of (1.8)+(1.22)


@book_only
def test_scale_height_V6_book_value():  # V6 (≤ 0.5 %; the book's worked scale-height example, private values)
    s = book()["stratification"]
    H = ch01.scale_height(s["scale_height_T_K"])
    assert rel(H / 1e3, s["scale_height_km"]["value"]) < 0.005
    assert rel(ch01.isothermal_pressure(H, 1e5, s["scale_height_T_K"]), 1e5 / np.e) < 1e-14
    band = s["isothermal_band_T_K"]
    z = np.linspace(0, band["height_km"] * 1e3, 701)
    T = ch01.standard_atmosphere(z)[0]
    dev = np.max(np.abs(T / band["value"] - 1)) * 100
    assert round(dev) <= band["within_percent"], dev  # USSA stays within the band to the book's rounding


# ====================================================================================================================
# C64 · Dimensional homogeneity (§1.11) — A item; with R01, R02, R03
# ====================================================================================================================
def test_dimensional_homogeneity_V2_pint_checks_of_chapter_laws():  # V2
    checks = {
        "(1.9)": (Q_(1e5, "Pa") - Q_(1000, "kg/m**3") * Q_(9.8, "m/s**2") * Q_(3, "m"), "[pressure]"),
        "(1.22)": (Q_(1.2, "kg/m**3") * Q_(287, "J/(kg*K)") * Q_(290, "K"), "[pressure]"),
        "(1.29)": (Q_(9.8, "m/s**2") / Q_(1025, "kg/m**3") * Q_(0.01, "kg/m**4"), "1/[time]**2"),
        "(1.30)": (Q_(9.8, "m/s**2") * Q_(3e-3, "1/K") * Q_(290, "K") / Q_(1005, "J/(kg*K)"), "[temperature]/[length]"),
        "(1.5)": (Q_(0.07, "N/m") / Q_(1, "mm"), "[pressure]"),
        "Ex1.1": (Q_(0.07, "N/m") / (Q_(1000, "kg/m**3") * Q_(9.8, "m/s**2") * Q_(1, "mm")), "[length]"),
        "H=RT/g": (Q_(287, "J/(kg*K)") * Q_(250, "K") / Q_(9.8, "m/s**2"), "[length]"),
        "Ex1.4": (Q_(1.2, "kg/m**3") * Q_(100, "m") ** 5 / Q_(0.02, "s") ** 2, "[energy]"),
        "c": ((Q_(1.4, "") * Q_(1e5, "Pa") / Q_(1.2, "kg/m**3")) ** 0.5, "[velocity]"),
    }
    for name, (q, dim) in checks.items():
        assert ch01.check_dimensions(q, dim), name
    out = ch01.dimensional_check(lambda rho, U: 0.5 * rho * U**2, "pressure", rho=Q_(998, "kg/m**3"), U=Q_(2, "m/s"))
    assert out.to("Pa").magnitude == pytest.approx(1996.0)
    with pytest.raises(AssertionError):  # a wrong law ρ g z² fails the test
        ch01.dimensional_check(lambda rho, g, z: rho * g * z**2, "pressure", rho=Q_(1000, "kg/m**3"),
                               g=Q_(9.8, "m/s**2"), z=Q_(2, "m"))


def test_dimensional_homogeneity_V1_dimension_vectors():  # V1 (§1.11 Step 2)
    cases = {"Pa": [1, -1, -2, 0], "J/(kg*K)": [0, 2, -2, -1], "m/s": [0, 1, -1, 0], "Pa*s": [1, -1, -1, 0],
             "W/m**2": [1, 0, -3, 0], "rad": [0, 0, 0, 0], "dimensionless": [0, 0, 0, 0]}
    for u, vec in cases.items():
        assert list(ch01.dimension_vector(u)) == vec, u
        d = ureg.parse_expression(u).dimensionality  # independent route: pint dimensionality dict
        assert [d.get(k, 0) for k in ("[mass]", "[length]", "[time]", "[temperature]")] == vec
    assert list(ch01.dimension_vector({"M": 1, "L": -3})) == [1, -3, 0, 0]
    assert list(ch01.dimension_vector([0, 2, -1])) == [0, 2, -1, 0] and list(ch01.dimension_vector(3.0)) == [0, 0, 0, 0]
    assert list(ch01.dimension_vector("m/s", basis=("L", "T"))) == [1, -1]
    for bad in (dict(q="A"), dict(q="Pa", basis=("L", "T")), dict(q={"X": 1}), dict(q=[1, 2])):
        with pytest.raises(ValueError):
            ch01.dimension_vector(**bad)
    with pytest.raises(ValueError):
        ch01.dimension_vector("mol", drop_substance=False)


def test_dimensional_homogeneity_V7_values_change_groups_do_not():  # V7 (C64: invariance under a change of units)
    v = {"dp": 100.0, "dx": 2.0, "d": 0.05, "eps": 1e-4, "U": 0.1, "rho": 1000.0, "mu": 1e-3}
    cgs = ch01.rescale_units(v, ch01.PIPE, "cgs")
    assert cgs["dp"] == pytest.approx(Q_(100, "Pa").to("dyn/cm**2").magnitude)  # independent: pint conversion
    assert cgs["mu"] == pytest.approx(Q_(1e-3, "Pa*s").to("poise").magnitude)
    imp = ch01.rescale_units(v, ch01.PIPE, "imperial")
    assert imp["dp"] == pytest.approx(Q_(100, "Pa").to("lb/(ft*s**2)").magnitude, rel=1e-12)
    assert imp["U"] == pytest.approx(Q_(0.1, "m/s").to("ft/s").magnitude, rel=1e-12)
    for g in ch01.pi_groups(ch01.PIPE, "dp", ("U", "d", "rho")):
        for sys_ in ("cgs", "imperial"):
            assert rel(ch01.group_value(g, ch01.rescale_units(v, ch01.PIPE, sys_)), ch01.group_value(g, v)) < 1e-12
    Pi1 = ch01.group_value({"dp": 1, "U": -2, "rho": -1}, v)
    assert Pi1 == pytest.approx(10.0)  # worked number: Π₁ = 10 in SI and cgs
    th = ch01.rescale_units({"T0": 300.0}, {"T0": "K"}, "imperial")["T0"]
    assert th == pytest.approx(Q_(300, "K").to("degR").magnitude)


def test_units_V1_celsius_kelvin_round_trip_and_pint():  # V1 + V7 (R03)
    T = np.array([-40.0, 0.0, 20.0, 100.0])
    assert np.allclose(ch01.kelvin_to_celsius(ch01.celsius_to_kelvin(T)), T, atol=1e-12)
    assert ch01.celsius_to_kelvin(0.0) == 273.15 and isinstance(ch01.celsius_to_kelvin(0.0), float)
    assert np.allclose(ch01.celsius_to_kelvin(T), Q_(T, "degC").to("K").magnitude, atol=1e-12)


# ====================================================================================================================
# C67 · Dimensional matrix and rank (§1.11) — A item; with C65, C66, C68
# ====================================================================================================================
def test_dimensional_matrix_V1_pipe_matrix_and_minors():  # V1 (Eq. 1.39, Step 3)
    A, names, rows = ch01.dimensional_matrix(ch01.PIPE)
    expected = np.array([[1, 0, 0, 0, 0, 1, 1], [-1, 1, 1, 1, 1, -3, -1], [-2, 0, 0, 0, -1, 0, -1]])
    assert np.array_equal(A, expected) and names == ["dp", "dx", "d", "eps", "U", "rho", "mu"] and rows == ["M", "L", "T"]
    assert ch01.minor_determinant(A, (0, 1, 2), (0, 1, 2)) == 0
    assert ch01.minor_determinant(A, (0, 1, 2), (4, 5, 6)) == -1
    assert ch01.minor_determinant(A, (0, 1, 2), (2, 4, 5)) == -1  # columns (d, U, ρ): the design's worked number
    assert ch01.rank_by_minors(A) == (3, (0, 1, 2), (0, 1, 4))
    assert ch01.rank_by_minors(np.zeros((2, 3))) == (0, (), ())
    with pytest.raises(ValueError):
        ch01.minor_determinant(A, (0, 1), (0, 1, 2))
    A2, _, rows2 = ch01.dimensional_matrix(ch01.PIPE, drop_zero_rows=False)
    assert rows2 == ["M", "L", "T", "Θ"] and not np.any(A2[3])
    # from-scratch cofactor expansion over every 3×3 minor (design §7)
    def det3(M):
        return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
                + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
    for cols in itertools.combinations(range(7), 3):
        assert det3(A[:, cols].tolist()) == ch01.minor_determinant(A, (0, 1, 2), cols)


def test_dimensional_matrix_V2_rank_agrees_with_numpy_and_sympy():  # V2 (independent linear-algebra routes)
    rng = np.random.default_rng(0)
    for k in range(200):
        m, n = (3, int(rng.integers(2, 8))) if k % 2 else (4, int(rng.integers(3, 7)))
        M = rng.integers(-3, 4, size=(m, n))
        if k % 5 == 0 and m > 2:
            M[2] = M[0] + 2 * M[1]  # force a dependent row
        r, rows, cols = ch01.rank_by_minors(M)
        assert r == np.linalg.matrix_rank(M) == sp.Matrix(M.tolist()).rank()
        if r:
            assert ch01.minor_determinant(M, rows, cols) != 0
            assert round(np.linalg.det(M[np.ix_(rows, cols)].astype(float))) == ch01.minor_determinant(M, rows, cols)
    Fr = ch01.minor_determinant([[Fraction(1, 2), 1], [1, 4]], (0, 1), (0, 1))
    assert Fr == 1


@book_only
def test_dimensional_matrix_V6_book_matrices_ranks_and_exponents():  # V6
    da = book()["dimensional_analysis"]
    A, names, _ = ch01.dimensional_matrix(ch01.PIPE)
    pm = da["pipe_matrix_1_39"]
    assert np.array_equal(A, np.array([pm["M"], pm["L"], pm["T"]]))
    assert ch01.minor_determinant(A, (0, 1, 2), (0, 1, 2)) == da["pipe_det_first_three_columns"]
    assert ch01.minor_determinant(A, (0, 1, 2), (4, 5, 6)) == da["pipe_det_last_three_columns"]
    assert ch01.rank_by_minors(A)[0] == da["pipe_rank"] and len(names) == da["pipe_n"]
    groups = ch01.pi_groups(ch01.PIPE, "dp", ("U", "d", "rho"))
    assert len(groups) == da["pipe_n_groups"]
    e = da["pipe_exponents_Pi1"]
    assert (groups[0]["U"], groups[0]["d"], groups[0]["rho"]) == (e["a_U"], e["b_d"], e["c_rho"])
    rx = da["rank_example_matrix"]
    assert ch01.rank_by_minors(rx["rows"])[0] == rx["rank"]
    presets = {"example_1_2": (ch01.SCALE_HEIGHT, ["M", "L", "T", "theta"]), "example_1_3": (ch01.PYTHAGORAS, ["L"]),
               "example_1_4": (ch01.BLAST, ["M", "L", "T"]), "example_1_5": (ch01.RAYLEIGH, ["M", "L", "T"])}
    for key, (preset, rows) in presets.items():
        ex = da[key]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            Am, _, _ = ch01.dimensional_matrix(preset)
            assert np.array_equal(Am, np.array([ex[r] for r in rows])), key
            assert ch01.rank_by_minors(Am)[0] == ex["rank"] and Am.shape[1] == ex["n"]
            assert len(ch01.pi_groups(preset, **{k: v for k, v in ch01.PRESET_INFO[key_map(key)].items()
                                                  if k in ("solution", "repeating")})) == ex["n_groups"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        g = ch01.pi_groups(ch01.SCALE_HEIGHT, "H", ("T0", "Mw", "g", "Ru"))[0]
    ex = da["example_1_2"]["exponents"]
    assert (g["T0"], g["Mw"], g["g"], g["Ru"]) == (ex["a_To"], ex["b_Mw"], ex["c_g"], ex["d_Ru"])


def key_map(example_key):
    return {"example_1_2": "scale_height", "example_1_3": "pythagoras", "example_1_4": "blast",
            "example_1_5": "rayleigh"}[example_key]


# ====================================================================================================================
# C69 · Buckingham's Π theorem (§1.11) — A item; with N22–N25, C70–C76, D28
# ====================================================================================================================
def test_buckingham_V2_derivation_D28_nullspace_rank_nullity():  # V2 (D28 ★★★)
    A = sp.Matrix([[1, 0, 0, 0, 0, 1, 1], [-1, 1, 1, 1, 1, -3, -1], [-2, 0, 0, 0, -1, 0, -1]])
    ns = A.nullspace()  # step 4
    assert A.rank() == 3 and len(ns) == A.shape[1] - A.rank() == 4  # step 5: rank–nullity
    assert all((A * k).is_zero_matrix for k in ns)  # step 3
    k1 = sp.Matrix([1, 0, 0, 0, -2, -1, 0])
    assert (A * k1).is_zero_matrix and sp.Matrix.hstack(*ns, k1).rank() == len(ns)  # step 7: in the span
    Anum, names, _ = ch01.dimensional_matrix(ch01.PIPE)
    groups = ch01.pi_groups(ch01.PIPE, "dp", ("U", "d", "rho"))
    E = sp.Matrix([[sp.Rational(Fraction(g.get(nm, 0)).numerator, Fraction(g.get(nm, 0)).denominator) for nm in names]
                   for g in groups])
    assert (sp.Matrix(Anum.tolist()) * E.T).is_zero_matrix  # every code group lies in null(A)
    assert sp.Matrix.hstack(*ns, *[E.row(i).T for i in range(4)]).rank() == 4  # and spans the same space
    # steps 9–11: choose base-unit scales that make the repeating variables 1; the rest become their Π values
    rng = np.random.default_rng(8)
    vals = dict(zip(names, rng.uniform(0.1, 10.0, 7)))
    Arep = Anum[:, [names.index(k) for k in ("U", "d", "rho")]].astype(float)
    ln_lam = np.linalg.solve(Arep.T, -np.log([vals["U"], vals["d"], vals["rho"]]))
    new = ch01.rescale_units(vals, ch01.PIPE, dict(zip(("M", "L", "T"), np.exp(ln_lam))) | {"Θ": 1.0})
    for k in ("U", "d", "rho"):
        assert abs(new[k] - 1.0) < 1e-12
    for g in groups:
        head = next(iter(g))
        assert rel(new[head], ch01.group_value(g, vals)) < 1e-12  # step 11: q_j′ = Π_j


def test_buckingham_V1_exponent_solve_matches_linear_algebra():  # V1 (from-scratch np.linalg.solve, design §7)
    A, names, _ = ch01.dimensional_matrix(ch01.PIPE)
    R = A[:, [names.index(k) for k in ("U", "d", "rho")]].astype(float)
    sol = np.linalg.solve(R, -A[:, names.index("dp")].astype(float))
    g = ch01.solve_exponents("dp", ("U", "d", "rho"), ch01.PIPE)
    assert np.allclose(sol, [float(g["U"]), float(g["d"]), float(g["rho"])]) and (g["U"], g["d"], g["rho"]) == (-2, 0, -1)
    assert isinstance(g["U"], Fraction)
    with pytest.raises(ValueError):
        ch01.solve_exponents("dp", ("d", "eps", "dx"), ch01.PIPE)  # singular repeating set (pure lengths)
    with pytest.raises(ValueError):
        ch01.solve_exponents("dp", ("U", "d"), ch01.PIPE)  # cannot cancel M
    assert ch01.group_latex(ch01.pi_groups(ch01.PIPE, "dp", ("U", "d", "rho"))[0]) == r"\frac{\Delta p}{U^{2} \rho}"
    expr = ch01.group_expression({"dp": 1, "U": -2, "rho": -1})
    dp, U, rho = sp.symbols("dp U rho", positive=True)
    assert sp.simplify(expr - dp / (U**2 * rho)) == 0
    assert ch01.group_latex({"lam": Fraction(1, 2)}) == r"\lambda^{1/2}"


def test_buckingham_V1_groups_count_dimensionless_independent_all_presets():  # V1 (property tests: n − r, zero dimension, independence)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        for key, info in ch01.PRESET_INFO.items():
            var = info["variables"]
            A, names, _ = ch01.dimensional_matrix(var)
            r = np.linalg.matrix_rank(A.astype(float))
            for rep in (info["repeating"], None):  # the book's set and the automatic one
                gs = ch01.pi_groups(var, info["solution"], rep)
                assert len(gs) == len(names) - r, key
                assert info["solution"] in gs[0] and gs[0][info["solution"]] == 1
                for g in gs:
                    assert not np.any(ch01.group_dimension(g, var)), (key, g)
                    qty = Q_(1.0, "")
                    for nm, e in g.items():  # independent dimension check with pint (kmol kept as a unit)
                        qty = qty * Q_(1.0, var[nm]) ** float(e)
                    dims = dict(qty.dimensionality)
                    dims.pop("[substance]", None)
                    assert all(v == 0 for v in dims.values()), (key, g, dims)
                assert ch01.groups_independent(gs, names)
                assert np.linalg.matrix_rank(ch01.exponent_matrix(gs, names)) == len(gs)
    g = ch01.pi_groups(ch01.PIPE, "dp", ("U", "d", "rho"))
    combo = {k: g[0].get(k, 0) - 2 * g[3].get(k, 0) for k in set(g[0]) | set(g[3])}  # Π₁/Π₄²
    assert not ch01.groups_independent([g[0], g[3], combo]) and ch01.groups_independent(g)
    assert {"theta0": 1} in ch01.pi_groups(ch01.PENDULUM, "tau", ("Lp", "m", "g"))  # dimensionless variable = own group
    with pytest.raises(ValueError):
        ch01.pi_groups(ch01.PIPE, "dp", ("U", "d"))
    with pytest.raises(ValueError):
        ch01.pi_groups(ch01.PIPE, "dp", ("dp", "d", "rho"))
    with pytest.raises(ValueError):
        ch01.pi_groups(ch01.PIPE, "zz")


def test_buckingham_V1_poiseuille_data_collapse():  # V1 (C72: Π₁ = 32 Π₂ Π₄ for laminar data)
    rng = np.random.default_rng(2)
    g = ch01.pi_groups(ch01.PIPE, "dp", ("U", "d", "rho"))
    for _ in range(50):
        v = {"dx": rng.uniform(0.5, 5), "d": rng.uniform(0.005, 0.05), "eps": 1e-5, "U": rng.uniform(0.01, 0.2),
             "rho": rng.uniform(800, 1200), "mu": rng.uniform(1e-3, 5e-2)}
        v["dp"] = ch01.poiseuille_pressure_drop(v["mu"], v["U"], v["dx"], v["d"])
        assert rel(ch01.group_value(g[0], v), 32 * ch01.group_value(g[1], v) * ch01.group_value(g[3], v)) < 1e-12


def test_examples_V1_blast_rayleigh_pythagoras():  # V1 + V7 (C73–C76)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        g = ch01.pi_groups(ch01.SCALE_HEIGHT, "H", ("T0", "Mw", "g", "Ru"))
    assert len(g) == 1
    vals = {"H": ch01.scale_height(288.15), "T0": 288.15, "Mw": ch01.M_W_AIR, "g": G, "Ru": ch01.R_U}
    assert rel(ch01.group_value(g[0], vals), 1.0) < 1e-12  # the constant of Ex. 1.2 is 1
    E, t, rho = 8.4e13, 0.025, 1.2
    D = ch01.blast_radius(E, t, rho)
    assert rel(ch01.blast_energy(D, t, rho), E) < 1e-12
    ts = np.logspace(-3, -1, 9)
    assert abs(np.polyfit(np.log(ts), np.log(ch01.blast_radius(E, ts, rho, K=ch01.TAYLOR_K_GAMMA14)), 1)[0] - 0.4) < 1e-12
    bg = ch01.pi_groups(ch01.BLAST, "E", ("D", "rho", "t"))
    assert bg == [{"E": 1, "D": -5, "rho": -1, "t": 2}]
    ratio = ch01.rayleigh_scattering_ratio(1e-24, 10.0, 450e-9) / ch01.rayleigh_scattering_ratio(1e-24, 10.0, 700e-9)
    assert rel(ratio, (700 / 450) ** 4) < 1e-12 and 5.8 < ratio < 5.9
    assert rel(ch01.rayleigh_scattering_ratio(2e-24, 10.0, 5e-7), 4 * ch01.rayleigh_scattering_ratio(1e-24, 10.0, 5e-7)) < 1e-12
    assert rel(ch01.rayleigh_scattering_ratio(1e-24, 20.0, 5e-7), ch01.rayleigh_scattering_ratio(1e-24, 10.0, 5e-7) / 4) < 1e-12
    rg = ch01.pi_groups(ch01.RAYLEIGH, "S", ("I", "lam"))
    assert len(rg) == 4 and not np.any(ch01.group_dimension({"S": 1, "I": -1, "d": 2, "lam": 4, "V": -2}, ch01.RAYLEIGH))
    beta = np.random.default_rng(1).uniform(0.1, 1.4, 20)
    C = 3.0
    A_, B_ = C * np.cos(beta), C * np.sin(beta)
    phi = ch01.pythagoras_phi(beta)
    assert np.allclose(A_**2 * phi + B_**2 * phi, C**2 * phi) and np.allclose(phi * C**2, 0.5 * A_ * B_)
    r, gcol, bcol = ch01.wavelength_to_rgb(450.0)
    assert bcol > r and ch01.wavelength_to_rgb(650.0)[0] == 1.0 and ch01.wavelength_to_rgb(900.0) == (0.0, 0.0, 0.0)
    assert [c.shape for c in ch01.wavelength_to_rgb(np.array([400.0, 500.0]))] == [(2,), (2,), (2,)]


@needs_ref
def test_blast_V5_taylor_constant():  # V5 (Taylor 1950 via Díaz 2020)
    tb = ref_json("benchmarks.json")["taylor_blast"]
    assert ch01.TAYLOR_K_GAMMA14 == tb["S_gamma_minus5"]
    assert rel(tb["S_gamma"] ** -5, ch01.TAYLOR_K_GAMMA14) < 0.003  # S = 1.032 is rounded to 4 s.f.
    assert rel(ch01.blast_energy(100.0, 0.02, 1.2, K=ch01.TAYLOR_K_GAMMA14),
               0.856 * ch01.blast_energy(100.0, 0.02, 1.2)) < 1e-14


def test_blast_V1_hemisphere_equals_free_sphere_of_twice_the_energy():  # V1 (image argument: ground burst E ≡ sphere 2E)
    K, rho = ch01.TAYLOR_K_GAMMA14, 1.2
    D = np.array([50.0, 100.0, 140.0])
    t = np.array([0.006, 0.016, 0.025])
    E_h = ch01.blast_energy(D, t, rho, K=K, geometry="hemisphere")
    assert rel(E_h, K * rho * D**5 / (2 * t**2)) < 1e-14  # E = K ρ D⁵/(2t²)
    assert rel(2 * E_h, ch01.blast_energy(D, t, rho, K=K)) < 1e-14  # the sphere formula returns the 2E equivalent
    E = 8.4e13
    ts = np.logspace(-3, -1, 7)
    assert rel(ch01.blast_radius(E, ts, rho, K=K, geometry="hemisphere"), ch01.blast_radius(2 * E, ts, rho, K=K)) < 1e-14
    assert rel(ch01.blast_energy(ch01.blast_radius(E, ts, rho, K=K, geometry="Hemisphere "), ts, rho, K=K,
                                 geometry="hemisphere"), E) < 1e-12  # round trip; name case/space-insensitive
    assert rel(ch01.blast_radius(E, 0.01, rho, geometry="hemisphere") / ch01.blast_radius(E, 0.01, rho), 2 ** 0.2) < 1e-14
    assert ch01.blast_energy(100.0, 0.02, rho, geometry="sphere") == ch01.blast_energy(100.0, 0.02, rho)  # default
    for fn, args in ((ch01.blast_energy, (100.0, 0.02, rho)), (ch01.blast_radius, (E, 0.02, rho))):
        with pytest.raises(ValueError):
            fn(*args, geometry="cylinder")
    with pytest.raises(TypeError):
        ch01.blast_energy(100.0, 0.02, rho, 1.0, "hemisphere")  # keyword-only


# ====================================================================================================================
# Machinery and drawing helpers called by the notebook (design Part C.0, C.7) — smoke tests
# ====================================================================================================================
def test_machinery_smoke_style_anim_interact_embed(tmp_path, capsys):  # smoke
    import matplotlib.pyplot as plt

    from fluidpy.core import anim, embed, interact, style

    fast = style.setup_notebook(fast=True)
    assert fast is True and style.setup_notebook(fast=False) is False and "accent" in style.COLORS
    fig, ax = plt.subplots()
    (ln,) = ax.plot([0, 1], [0, 1])
    path = style.savefig(fig, "ch01", "smoke_test", root=tmp_path)
    assert path.exists() and path.stat().st_size > 0
    a = anim.animate(lambda i: ln.set_ydata([0, i]), frames=3, fig=fig)
    html = anim.animation_html(a, player="frames")
    assert "<" in html and "script" in html.lower()
    fig2, ax2 = plt.subplots()
    (ln2,) = ax2.plot([0, 1], [0, 1])
    anim.show_animation(anim.animate(lambda i: ln2.set_ydata([0, i]), frames=2, fig=fig2), player="frames")
    sf = interact.slider_figure(lambda H: {"p": ([0, 1], [1.0, np.exp(-1 / H)])}, "H", [1.0, 2.0, 4.0], unit="km")
    assert len(sf.data) == 3 and len(sf.layout.sliders[0].steps) == 3
    af = interact.animate_figure(lambda t: {"u": ([0, 1], [0.0, t])}, [0.0, 0.5, 1.0])
    assert len(af.frames) == 3
    w = interact.live(lambda N2=1e-4: N2, N2=(-1e-4, 1e-4, 1e-5))
    assert w is not None
    embed.show_viz("ch01", "continuum_averaging_volume")
    out = capsys.readouterr().out
    assert "fluidpy-viz" in out or "HTML" in out
    assert "ch01/continuum_averaging_volume" in embed.viz_html("ch01", "continuum_averaging_volume", mode="srcdoc")
    plt.close("all")


def test_drawing_helpers_smoke():  # smoke (scripts/ch01_drawings.py, scripts/ch01_book_map.py)
    import matplotlib.pyplot as plt

    from scripts.ch01_book_map import draw_book_map
    from scripts.ch01_drawings import draw_capillary, draw_cube_forces, draw_parcel_column, draw_triangle_split, draw_wedge

    for fn, kw in ((draw_cube_forces, {}), (draw_wedge, dict(dz=1.0, theta=0.61)), (draw_capillary, dict(R=1.0, alpha=1.22, h=3.0)),
                   (draw_parcel_column, dict(z0=500.0, zeta=300.0)), (draw_triangle_split, dict(beta=0.52))):
        fig, ax = plt.subplots()
        assert fn(ax, **kw) is ax and len(ax.get_children()) > 5
        plt.close(fig)
    ax = draw_book_map(highlight=("ch13",))
    assert ax is not None and len(ax.get_children()) > 10
    plt.close("all")
