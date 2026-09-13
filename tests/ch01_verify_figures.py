"""Chapter 1 verification figures and report metrics (not collected by pytest).

Reproduces the chapter's key plots with our own code into ``outputs/ch01/verify/`` (git-ignored) and prints the numbers
quoted in ``reports/ch01_verification.md`` (convergence orders, benchmark errors, invariant residuals).

Run: ``.venv/Scripts/python.exe tests/ch01_verify_figures.py``
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fluidpy import ch01_introduction as ch01  # noqa: E402
from fluidpy.core.style import COLORS, use_style  # noqa: E402
from tools.convergence import observed_order, pairwise_orders  # noqa: E402

OUT = ROOT / "outputs" / "ch01" / "verify"
REF = ROOT / "reference" / "ch01"
G = ch01.G0


def rows(name):
    with open(REF / name, newline="", encoding="utf-8") as f:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"figure: outputs/ch01/verify/{name}.png")


def fig_1_9():
    z = np.linspace(0, 2000, 401)
    col = ch01.synthetic_boundary_layer_column(z)
    fig, (a, b) = plt.subplots(1, 2, figsize=(9, 4.5), sharey=True)
    colors = {"stable": COLORS["teal"], "neutral": COLORS["amber"], "unstable": COLORS["rose"]}
    for lab, c in colors.items():
        m = col["stability"] == lab
        a.scatter(col["T"][m] - 273.15, z[m], s=4, color=c, label=lab)
        b.scatter(col["theta"][m] - 273.15, z[m], s=4, color=c)
    for T0 in np.arange(5, 30, 5):  # neutral reference lines: dry adiabats
        a.plot(T0 + ch01.adiabatic_lapse_rate() * z, z, color=COLORS["muted"], lw=0.7)
        b.axvline(T0, color=COLORS["muted"], lw=0.7)
    a.set_xlabel("T [°C]"); a.set_ylabel("z [m]"); b.set_xlabel("θ [°C]")
    a.set_xlim(0, 25); b.set_xlim(10, 30); a.legend(markerscale=4)
    a.set_title("synthetic T(z) + dry adiabats"); b.set_title("θ(z)")
    save(fig, "fig1_9_T_theta_synthetic")
    i = int(np.argmin(np.abs(z - 400)))
    print(f"synthetic mixed layer: dT/dz = {col['dT_dz'][i]*1e3:.2f} K/km, Γa = {ch01.adiabatic_lapse_rate()*1e3:.4f} K/km,"
          f" N² = {col['N2'][i]:.3e} 1/s², label = {col['stability'][i]}")


def fig_parcel():
    t = np.linspace(0, 3000, 601)
    fig, ax = plt.subplots(figsize=(7, 4))
    for (drho, lab, c) in ((-0.01, "stable", COLORS["teal"]), (0.0, "neutral", COLORS["amber"]), (0.01, "unstable", COLORS["rose"])):
        N2 = ch01.brunt_vaisala_sq(1025.0, drho, 0.0)
        ax.plot(t, ch01.parcel_displacement(t, 5.0, N2), color=c, lw=3, alpha=0.35)
        ax.plot(t, ch01.parcel_ode_from_gradients(t, 5.0, 1025.0, drho, 0.0), color=c, ls="--", label=f"{lab} (nonlinear)")
    ax.set_ylim(-20, 60); ax.set_xlabel("t [s]"); ax.set_ylabel("ζ [m]"); ax.legend()
    ax.set_title("parcel: linear (thick) vs nonlinear (dashed); unstable run ends at the cap")
    save(fig, "parcel_regimes")


def fig_couette():
    h, U, nu = 1e-3, 1.0, 1e-6
    y = np.linspace(0, h, 41)
    dy = y[1] - y[0]
    dt = ch01.stable_time_step(nu, dy, 0.8)
    times = [0.01, 0.05, 0.15, 0.5]
    nsteps = int(round(times[-1] / dt))
    F = ch01.ftcs_diffusion_1d(np.where(np.isclose(y, h), U, 0.0), nu, dy, dt, nsteps)
    fig, (a, b) = plt.subplots(1, 2, figsize=(10, 4))
    for tt in times:
        k = int(round(tt / dt))
        a.plot(F[k], y * 1e3, "o", ms=3, color=COLORS["accent"])
        a.plot(ch01.couette_startup_profile(y, k * dt, U, h, nu, nterms=None), y * 1e3, color=COLORS["ink"], lw=1)
    a.set_xlabel("u [m/s]"); a.set_ylabel("y [mm]"); a.set_title("FTCS (dots) vs series (lines), t = 0.01…0.5 s")
    hs, errs = [], []
    for N in (11, 21, 41, 81, 161):
        yy = np.linspace(0, h, N); d = yy[1] - yy[0]
        ns = int(round(0.05 / (0.4 * d**2 / nu)))
        u0 = np.zeros(N); u0[-1] = U
        FF = ch01.ftcs_diffusion_1d(u0, nu, d, 0.05 / ns, ns, save_every=ns)
        errs.append(np.max(np.abs(FF[-1] - ch01.couette_startup_profile(yy, 0.05, U, h, nu, nterms=None)))); hs.append(d)
    b.loglog(hs, errs, "o-", color=COLORS["accent"]); b.loglog(hs, errs[0] * (np.array(hs) / hs[0]) ** 2, "--", color=COLORS["muted"], label="slope 2")
    b.set_xlabel("Δy [m]"); b.set_ylabel("max error [m/s]"); b.legend()
    save(fig, "couette_ftcs_vs_series")
    print(f"FTCS vs Couette (r = 0.4, t = 0.05 h²/ν): errors {['%.3e' % e for e in errs]}, order {observed_order(hs, errs):.3f},"
          f" pairwise {[round(p, 3) for p in pairwise_orders(hs, errs)]}")


def fig_continuum():
    n, m = ch01.number_density(1.225, ch01.M_W_AIR), ch01.molecular_mass(ch01.M_W_AIR)
    L = np.logspace(-8.5, -5, 15)
    mean, std = ch01.sample_density(L, n, m, n_samples=2000, seed=1)
    fig, (a, b) = plt.subplots(1, 2, figsize=(10, 4))
    a.loglog(L, std / mean, "o", color=COLORS["accent"], label="sample_density")
    a.loglog(L, ch01.density_noise_expected(L, n), color=COLORS["ink"], label="(nL³)^−1/2")
    a.set_xlabel("box side L [m]"); a.set_ylabel("relative noise"); a.legend()
    LL = np.logspace(-4, 1, 200)
    b.semilogx(LL, ch01.box_average_density(LL, 1.225, 0.1, 1.0), color=COLORS["teal"])
    b.axhline(1.225 * 1.1, ls=":", color=COLORS["muted"]); b.axhline(1.225, ls=":", color=COLORS["muted"])
    b.set_xlabel("box side L [m]"); b.set_ylabel("noiseless box density [kg/m³]")
    save(fig, "continuum_noise_and_drift")
    slope = np.polyfit(np.log(L), np.log(std / mean), 1)[0]
    print(f"continuum noise slope = {slope:.4f} (design −1.5)")


def fig_lapse():
    x = np.linspace(-15e-3, 10e-3, 501)
    N2 = ch01.brunt_vaisala_sq_from_lapse(260.0, x)
    codes = np.array([ch01.lapse_rate_stability(float(v)).code for v in x])
    codes_m = np.array([ch01.lapse_rate_stability(float(v), convention="meteorology").code for v in x])
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(x * 1e3, N2 * 1e4, color=COLORS["ink"])
    ax.scatter(x * 1e3, codes * 3, s=3, color=COLORS["accent"], label="lapse_rate_stability code ×3 (Kundu)")
    ax.axvline(ch01.adiabatic_lapse_rate() * 1e3, color=COLORS["rose"], ls="--", label="Γa (Kundu)")
    ax.axhline(0, color=COLORS["muted"], lw=0.7)
    sec = ax.secondary_xaxis("top", functions=(lambda v: -v, lambda v: -v)); sec.set_xlabel("meteorology Γ = −dT/dz [K/km]")
    ax.set_xlabel("dT/dz [K/km]"); ax.set_ylabel("N² [1e−4 s⁻²]"); ax.legend(loc="upper left")
    save(fig, "lapse_rate_sweep")
    print(f"lapse sweep: codes identical in both conventions: {np.array_equal(codes, codes_m)};"
          f" sign(N²) == code: {np.array_equal(np.sign(np.round(N2, 15)), codes)}")


def fig_ussa():
    c = json.loads((REF / "ussa1976_constants.json").read_text(encoding="utf-8"))
    t1 = rows("ussa1976_table1.csv")
    Z = np.array([r["Z_km"] for r in t1])
    H = 1e3 * c["r0_km"] * Z / (c["r0_km"] + Z)
    T, p, rho = ch01.standard_atmosphere(H)
    zz = np.linspace(0, 50000, 500)
    Tc, pc, _ = ch01.standard_atmosphere(zz)
    fig, (a, b) = plt.subplots(1, 2, figsize=(9, 4.5), sharey=True)
    a.plot(Tc, zz / 1e3, color=COLORS["ink"]); a.plot([r["T_K"] for r in t1], H / 1e3, "o", color=COLORS["orange"])
    b.semilogx(pc, zz / 1e3, color=COLORS["ink"]); b.semilogx([r["p_Pa"] for r in t1], H / 1e3, "o", color=COLORS["orange"])
    a.set_xlabel("T [K]"); b.set_xlabel("p [Pa]"); a.set_ylabel("geopotential H [km]")
    a.set_title("standard_atmosphere (line) vs PDAS USSA-1976 (dots)")
    save(fig, "ussa1976_vs_table")
    eT = np.max(np.abs(T - np.array([r["T_K"] for r in t1])))
    ep = np.max(np.abs(p / np.array([r["p_Pa"] for r in t1]) - 1))
    er = np.max(np.abs(rho / np.array([r["rho_kg_m3"] for r in t1]) - 1))
    ec = np.max(np.abs(ch01.perfect_gas_sound_speed(np.array([r["T_K"] for r in t1])) / np.array([r["c_m_s"] for r in t1]) - 1))
    print(f"USSA Table 1 (0–50 km): max |ΔT| = {eT:.2e} K, max rel p = {ep:.2e}, rho = {er:.2e}, c = {ec:.2e}")
    t2 = rows("ussa1976_table2.csv")
    Td = {r["Z_km"]: r for r in t1}
    emu = max(abs(ch01.sutherland_viscosity(Td[r["Z_km"]]["T_K"]) / r["mu_Pa_s"] - 1) for r in t2)
    en = max(abs(ch01.number_density(Td[r["Z_km"]]["rho_kg_m3"], ch01.M_W_AIR) / r["n_per_m3"] - 1) for r in t2)
    eV = max(abs(ch01.mean_molecular_speed(Td[r["Z_km"]]["T_K"], ch01.molecular_mass(ch01.M_W_AIR)) / r["V_m_s"] - 1) for r in t2)
    eH = max(abs(ch01.scale_height(Td[r["Z_km"]]["T_K"], g=Td[r["Z_km"]]["g_m_s2"]) / r["Hp_m"] - 1) for r in t2)
    enu = max(abs(ch01.kinematic_viscosity(r["mu_Pa_s"], Td[r["Z_km"]]["rho_kg_m3"]) / r["nu_m2_s"] - 1) for r in t2)
    print(f"USSA Table 2: rel μ {emu:.2e}, ν {enu:.2e}, n {en:.2e}, mean speed {eV:.2e}, H_p {eH:.2e}")
    print(f"perfect_gas_density(101325, 288.15) rel err vs 1.2250: {ch01.perfect_gas_density(101325.0, 288.15)/1.2250-1:.2e}")
    comp = c["composition"]
    M0 = sum(a * b for a, b in comp.values()) / sum(b for _, b in comp.values())
    print(f"USSA Table 3 mean molecular weight {M0:.5f} vs M_W_AIR {ch01.M_W_AIR} (rel {ch01.M_W_AIR/M0-1:.2e})")


def fig_pi_collapse():
    rng = np.random.default_rng(2)
    g = ch01.pi_groups(ch01.PIPE, "dp", ("U", "d", "rho"))
    x, y = [], []
    for _ in range(80):
        v = {"dx": rng.uniform(0.5, 5), "d": rng.uniform(0.005, 0.05), "eps": 1e-5, "U": rng.uniform(0.01, 0.2),
             "rho": rng.uniform(800, 1200), "mu": rng.uniform(1e-3, 5e-2)}
        v["dp"] = ch01.poiseuille_pressure_drop(v["mu"], v["U"], v["dx"], v["d"])
        x.append(ch01.group_value(g[1], v) * ch01.group_value(g[3], v)); y.append(ch01.group_value(g[0], v))
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.loglog(x, y, "o", color=COLORS["accent"]); xx = np.array([min(x), max(x)])
    ax.loglog(xx, 32 * xx, color=COLORS["ink"], label="Π₁ = 32 Π₂Π₄")
    ax.set_xlabel("Π₂ Π₄ = (Δx/d)(μ/ρUd)"); ax.set_ylabel("Π₁ = Δp/ρU²"); ax.legend()
    save(fig, "pi_collapse_pipe")


def fig_sigma():
    t = rows("iapws_sigma.csv")
    T = np.linspace(273.16, 647.0, 300)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(T - 273.15, ch01.surface_tension_water(T) * 1e3, color=COLORS["blue"], label="surface_tension_water")
    ax.errorbar([r["t_C"] for r in t], [r["sigma_exp_mN_m"] for r in t], yerr=[r["uncertainty_mN_m"] for r in t], fmt="o",
                color=COLORS["orange"], label="IAPWS Table 1 (exp.)")
    ax.set_xlabel("t [°C]"); ax.set_ylabel("σ [mN/m]"); ax.legend()
    save(fig, "iapws_surface_tension")
    e = max(abs(ch01.surface_tension_water(r["t_C"] + 273.15) * 1e3 - r["sigma_calc_mN_m"]) for r in t)
    print(f"IAPWS σ: max |ours − calculated column| = {e:.4f} mN/m")


def metrics():
    b = json.loads((REF / "benchmarks.json").read_text(encoding="utf-8"))
    lam = ch01.mean_free_path_air(300.0) * 1e9
    print(f"Jennings mean free path at 300 K: {lam:.3f} nm vs 67.3 nm (rel {lam/67.3-1:.2e})")
    mu = ch01.water_viscosity(298.15) * 1e6
    print(f"water viscosity 298.15 K: {mu:.3f} µPa s vs IAPWS 889.7351 (rel {mu/889.7351-1:.2e})")
    met = -ch01.adiabatic_lapse_rate() * 1e3
    print(f"dry adiabatic lapse rate (meteorology) {met:.4f} K/km vs AMS 9.8 (rel {met/9.8-1:.2e})")
    print(f"Taylor: 1.032^-5 = {1.032**-5:.4f} vs 0.856 (rel {1.032**-5/0.856-1:.2e})")
    # convergence studies
    def study(name, hs, errs, design):
        print(f"order[{name}] = {observed_order(hs, errs):.3f} (design {design}); pairwise {[round(p, 3) for p in pairwise_orders(hs, errs)]}")
    steps = [1e-1, 5e-2, 2.5e-2, 1.25e-2]
    study("partial_derivative", steps, [abs(ch01.partial_derivative(lambda x, y: np.sin(x) * y, "x", {"x": 0.7, "y": 2.0}, step=s) - 2 * np.cos(0.7)) for s in steps], 2)
    hs, errs = [], []
    for n in (11, 21, 41, 81):
        d = ch01.process_path("isentropic", (0.8, 300.0), (1.6, 300.0 * 0.5**0.4), n=n)
        errs.append(abs(ch01.process_heat_work(d["v"], d["T"])["w"][-1] - ch01.path_heat_work_totals("isentropic", 0.8, 300.0, 1.6, 300.0 * 0.5**0.4)["w"])); hs.append(0.8 / (n - 1))
    study("process_heat_work trapezoid", hs, errs, 2)
    hs, errs = [], []
    for N in (21, 41, 81, 161):
        yy = np.linspace(0, 1, N)
        errs.append(np.max(np.abs(ch01.shear_stress_profile(np.sin(3 * yy), yy, 2.0) - 6 * np.cos(3 * yy)))); hs.append(yy[1])
    study("shear_stress_profile stencils", hs, errs, 2)
    hs, errs = [], []
    for N in (61, 121, 241, 481):
        yy = np.linspace(-3, 3, N); d = yy[1] - yy[0]; ns = int(round(0.02 / (0.25 * d**2)))
        f = ch01.ftcs_diffusion_1d(ch01.gaussian_spreading(yy, 0.0, 1.0, t0=0.01), 1.0, d, 0.02 / ns, ns, bc=("neumann", "neumann"), save_every=ns)[-1]
        errs.append(np.max(np.abs(f - ch01.gaussian_spreading(yy, 0.02, 1.0, t0=0.01)))); hs.append(d)
    study("FTCS Gaussian (zero-flux)", hs, errs, 2)
    t = np.linspace(0, 1500, 31)
    z0s = [8.0, 4.0, 2.0, 1.0]
    N2 = ch01.brunt_vaisala_sq(1025.0, -0.01, -0.004)
    study("parcel_ode_from_gradients → linear (vs ζ0)", z0s, [np.max(np.abs(ch01.parcel_ode_from_gradients(t, z, 1025.0, -0.01, -0.004) - ch01.parcel_displacement(t, z, N2))) for z in z0s], 2)
    z0a = [200.0, 100.0, 50.0, 25.0]
    N2a = ch01.brunt_vaisala_sq_from_lapse(288.15, -6.5e-3)
    study("parcel_ode_atmosphere → linear (vs ζ0)", z0a, [np.max(np.abs(ch01.parcel_ode_atmosphere(t, z, 288.15, -6.5e-3) - ch01.parcel_displacement(t, z, N2a))) for z in z0a], 2)
    hs, errs = [], []
    for n in (50, 100, 200, 400):
        z = np.linspace(0.0, 3000.0, n + 1); dz = z[1]; p = 1e5
        for k in range(n):
            p = p - p / (ch01.R_AIR * 260.0) * G * dz
        errs.append(abs(p - ch01.isothermal_pressure(3000.0, 1e5, 260.0))); hs.append(dz)
    study("from-scratch Euler hydrostatic", hs, errs, 1)
    hs, errs = [], []
    for n in (100, 200, 400, 800):
        rho = np.linspace(1.225, 2.45, n + 1); p = 101325.0
        for k in range(n):
            p += 1.4 * p / rho[k] * (rho[k + 1] - rho[k])
        errs.append(abs(p - ch01.isentropic_pressure(2.45, 101325.0, 1.225))); hs.append(rho[1] - rho[0])
    study("from-scratch Euler isentrope", hs, errs, 1)
    dz = np.array([1e-1, 5e-2, 2.5e-2, 1.25e-2])
    study("wedge p2 − p1", dz, ch01.wedge_pressure_difference(1000.0, dz, 0.7)["p2_minus_p1"], 1)
    # invariants
    y = np.linspace(-3, 3, 121); dy = y[1] - y[0]
    F = ch01.ftcs_diffusion_1d(np.exp(-(y - 0.5) ** 2 / 0.1), 0.5, dy, 0.45 * dy**2 / 0.5, 2000, bc=("neumann", "neumann"), save_every=500)
    I = [np.trapezoid(f, dx=dy) for f in F]
    print(f"FTCS zero-flux integral drift: {np.max(np.abs(np.array(I) - I[0])) / I[0]:.2e}")
    zz = np.linspace(0, 10000, 41)
    pz, rz, Tz = ch01.atmosphere_from_temperature(zz, lambda q: 300.0 - G / ch01.CP_AIR * np.asarray(q))
    th = ch01.potential_temperature(Tz, pz)
    print(f"θ along dry adiabat: ptp/θ = {np.ptp(th)/th[0]:.2e}")
    T, p, rho = ch01.standard_atmosphere(np.linspace(0, 11000, 1101))
    print(f"θ ρθ − p_ref/R (USSA 0–11 km): {np.max(np.abs(ch01.potential_temperature(T, p)*ch01.potential_density(rho, p)/(ch01.P_REF/ch01.R_AIR)-1)):.2e}")
    legs = [ch01.process_path("isothermal", (0.8, 300.0), (1.6, 300.0), n=2001), ch01.process_path("isochoric", (1.6, 300.0), (1.6, 450.0), n=2001),
            ch01.process_path("isobaric", (1.6, 450.0), (0.8, 225.0), n=2001), ch01.process_path("isochoric", (0.8, 225.0), (0.8, 300.0), n=2001)]
    v, Tp = ch01.join_paths(*[(d["v"], d["T"]) for d in legs])
    o = ch01.process_heat_work(v, Tp)
    print(f"closed reversible cycle: Δs = {ch01.entropy_change_reversible(o['q'], Tp):.2e} J/(kg K), Δe = {o['de'][-1]:.2e}, net w = {o['w'][-1]:.1f} J/kg")
    box = ch01.net_pressure_force_on_box(lambda x, yy, zq: ch01.hydrostatic_pressure_uniform(zq, ch01.P_ATM, 1025.0), (0.1, 0.3, -0.2, 0.3, -5.0, -4.6))
    print(f"box force: Fx, Fy = {box[0]:.2e}, {box[1]:.2e} N; Fz/ρgV − 1 = {box[2]/ch01.buoyancy_force(1025.0, 0.04)-1:.2e}")
    zg = np.linspace(0.0, 6000.0, 6001)
    Tf = lambda q: 288.0 - 6e-3 * np.asarray(q) + 3.0 * np.sin(np.asarray(q) / 700.0)  # noqa: E731
    pg, rg, Tg = ch01.atmosphere_from_temperature(zg, Tf)
    N2l = ch01.brunt_vaisala_sq_from_lapse(Tg, -6e-3 + 3.0 / 700.0 * np.cos(zg / 700.0))
    N2r = ch01.brunt_vaisala_sq(rg, np.gradient(rg, zg, edge_order=2), ch01.isentropic_density_gradient(rg, ch01.perfect_gas_sound_speed(Tg)))
    thg = ch01.potential_temperature(Tg, pg)
    N2t = ch01.brunt_vaisala_sq_from_theta(thg, np.gradient(thg, zg, edge_order=2))
    s = np.max(np.abs(N2l))
    print(f"N² three routes: max|ρ-form − lapse| / max|N²| = {np.max(np.abs(N2r - N2l)[5:-5])/s:.2e}, θ-form = {np.max(np.abs(N2t - N2l)[5:-5])/s:.2e}")


if __name__ == "__main__":
    use_style()
    fig_1_9(); fig_parcel(); fig_couette(); fig_continuum(); fig_lapse(); fig_ussa(); fig_pi_collapse(); fig_sigma()
    metrics()
