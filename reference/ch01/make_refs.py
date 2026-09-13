"""Write the cited public benchmark data used by ``tests/test_ch01.py`` (V5 evidence) into ``reference/ch01/``.

Every number below was read from the cited public source on 2026-09-12/13 (see ``SOURCES.md`` in this folder for the
exact document, page/table and how it was obtained). Nothing here comes from the textbook; book-quoted numbers live in
the git-ignored ``tests/book_values_ch01.json``.

Run: ``.venv/Scripts/python.exe reference/ch01/make_refs.py``  (idempotent; overwrites the files).
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

# --- CODATA 2018/2022 exact defining constants (NIST CUU pages "Value?k", "Value?na", "Value?r") ----------------------
CODATA = {
    "k_B_J_per_K": 1.380649e-23,          # exact
    "N_A_per_mol": 6.02214076e23,         # exact
    "R_J_per_mol_K": 8.314462618,         # exact (printed 8.314 462 618...)
    "_source": "NIST CODATA, https://physics.nist.gov/cgi-bin/cuu/Value?k , ?na , ?r (fetched 2026-09-12)",
}

# --- U.S. Standard Atmosphere 1976, NASA-TM-X-74335 (NTRS 19770009539), Tables 2-4 (read from the scanned PDF) --------
USSA_CONSTANTS = {
    "g0_m_per_s2": 9.80665,               # Table 2, Category II
    "P0_Pa": 101325.0,                    # Table 2 (1.013250e5)
    "T0_K": 288.15,                       # Table 2
    "R_star_J_per_kmol_K": 8314.32,       # p. 3, paragraph on R* (Table 2 prints 8.31432e-3: NASA errata sheet)
    "r0_km": 6356.766,                    # Table 2 prints 6.356766e6 km; NASA errata sheet: too large by 1000
    "sutherland_beta_kg_per_m_s_sqrtK": 1.458e-6,   # Table 2
    "sutherland_S_K": 110.4,              # Table 2 prints 110 K; NASA errata sheet (p. 242 of the PDF): 110.4 K
    "gamma": 1.40,                        # Table 2
    # Table 3: molecular weight M_i [kg/kmol] and fractional volume F_i of sea-level dry air (first 10 species)
    "composition": {
        "N2": [28.0134, 0.78084], "O2": [31.9988, 0.209476], "Ar": [39.948, 0.00934], "CO2": [44.00995, 0.000314],
        "Ne": [20.183, 0.00001818], "He": [4.0026, 0.00000524], "Kr": [83.80, 0.00000114],
        "Xe": [131.30, 0.000000087], "CH4": [16.04303, 0.000002], "H2": [2.01594, 0.0000005],
    },
    # Table 4: geopotential base heights H_b [km'] and molecular-scale temperature gradients L_M,b [K/km']
    "layers_H_km": [0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.8520],
    "layers_L_K_per_km": [-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0],
    "_source": "NASA-TM-X-74335 (1976), https://ntrs.nasa.gov/api/citations/19770009539/downloads/19770009539.pdf, "
               "doc pp. 2-3 (PDF pp. 18-19) and errata sheet (PDF p. 242); fetched 2026-09-12",
}

# PDAS "Tables of the U.S. Standard Atmosphere, 1976", Table 1 (SI) and Table 2 (SI), geometric altitude Z [km],
# columns copied from the HTML page (https://www.pdas.com/bigtables.html, fetched 2026-09-12).
# Z_km, H_km(rounded), T_K, p_Pa, rho_kg_m3, c_m_s, g_m_s2
USSA_TABLE1 = [
    (0, 0.0, 288.150, 1.0132e5, 1.2250e0, 340.29, 9.8066),
    (5, 5.0, 255.676, 5.4048e4, 7.3643e-1, 320.55, 9.7912),
    (10, 10.0, 223.252, 2.6500e4, 4.1351e-1, 299.53, 9.7759),
    (15, 15.0, 216.650, 1.2112e4, 1.9476e-1, 295.07, 9.7605),
    (20, 19.9, 216.650, 5.5293e3, 8.8910e-2, 295.07, 9.7452),
    (25, 24.9, 221.552, 2.5492e3, 4.0084e-2, 298.39, 9.7300),
    (30, 29.9, 226.509, 1.1970e3, 1.8410e-2, 301.71, 9.7147),
    (40, 39.7, 250.350, 2.8714e2, 3.9957e-3, 317.19, 9.6844),
    (45, 44.7, 264.164, 1.4910e2, 1.9663e-3, 325.82, 9.6693),
    (50, 49.6, 270.650, 7.9779e1, 1.0269e-3, 329.80, 9.6542),
]
# Z_km, H_km, mu_Pa_s, nu_m2_s, Hp_m, n_per_m3, V_m_s(mean particle speed), mean_free_path_m, M_kg_kmol
USSA_TABLE2 = [
    (0, 0.0, 1.7894e-5, 1.4607e-5, 8435, 2.5470e25, 458.94, 6.6332e-8, 28.964),
    (5, 5.0, 1.6282e-5, 2.2110e-5, 7496, 1.5312e25, 432.31, 1.1034e-7, 28.964),
    (10, 10.0, 1.4577e-5, 3.5250e-5, 6555, 8.5976e24, 403.97, 1.9651e-7, 28.964),
    (20, 19.9, 1.4216e-5, 1.5989e-4, 6382, 1.8486e24, 397.95, 9.1393e-7, 28.964),
    (30, 29.9, 1.4753e-5, 8.0132e-4, 6693, 3.8278e23, 406.91, 4.4137e-6, 28.964),
    (50, 49.6, 1.7037e-5, 1.6590e-2, 8047, 2.1350e22, 444.79, 7.9130e-5, 28.964),
]

# IAPWS R1-76(2014), Table 1: t [°C], experimental sigma [mN/m], uncertainty [mN/m], calculated sigma [mN/m]
IAPWS_SIGMA = [
    (0.01, 75.64, 0.38, 75.65), (5, 74.94, 0.37, 74.94), (10, 74.23, 0.37, 74.22), (15, 73.49, 0.37, 73.49),
    (20, 72.74, 0.36, 72.74), (25, 71.98, 0.36, 71.97), (30, 71.19, 0.36, 71.19), (40, 69.59, 0.35, 69.60),
    (50, 67.93, 0.34, 67.94), (60, 66.24, 0.33, 66.24), (80, 62.68, 0.31, 62.67), (100, 58.92, 0.29, 58.91),
    (160, 46.58, 0.23, 46.59), (200, 37.68, 0.22, 37.67),
]
IAPWS_SIGMA_EQ = {"B_mN_per_m": 235.8, "b": -0.625, "mu": 1.256, "Tc_K": 647.096,
                  "form": "sigma = B tau^mu (1 + b tau), tau = 1 - T/Tc",
                  "_note": "the PDF text layer drops the minus sign of b; the tabulated calculated column fixes b = -0.625"}

# IAPWS R12-08 (2008 viscosity formulation), Table 4 computer-program verification points
IAPWS_VISCOSITY_CHECK = [
    {"T_K": 298.15, "rho_kg_m3": 998.0, "mu_microPa_s": 889.735100},
    {"T_K": 373.15, "rho_kg_m3": 1000.0, "mu_microPa_s": 307.883622},
]

# Jennings (1988) via Tsalikis, Mavrantzas & Pratsinis (2024), open access (CC BY 4.0)
JENNINGS = {
    "mean_free_path_air_300K_1atm_nm": 67.3,
    "phi": 0.4987445,
    "formula": "lambda = sqrt(pi/8) * mu / (phi * sqrt(rho * P))  (their Eq. 14)",
    "MD_value_nm": 38.5, "MD_uncertainty_nm": 1.0,
    "_source": "Tsalikis et al., Aerosol Sci. Technol. 58(8) (2024), doi:10.1080/02786826.2024.2333859; "
               "ETH research collection doi:10.3929/ethz-b-000669211",
}

# Taylor (1950) via Díaz (2020), arXiv:2009.05674
TAYLOR = {"gamma": 1.4, "S_gamma_minus5": 0.856, "S_gamma": 1.032,
          "_source": "J. S. Díaz, 'Explosion analysis from images: Trinity and Beirut', arXiv:2009.05674, p. 3"}

# AMS Glossary of Meteorology, 'adiabatic lapse rate' (definition 1, dry-adiabatic lapse rate)
AMS = {"dry_adiabatic_lapse_rate_C_per_km": 9.8, "approximate": True, "convention": "rate of DECREASE, -dT/dz",
       "definition": "g/c_pd", "_source": "https://glossary.ametsoc.org/wiki/Dry-adiabatic_lapse_rate (fetched 2026-09-13)"}


# UNESCO EOS-80 check values: Fofonoff & Millard (1983), Unesco Tech. Pap. Mar. Sci. 44, p. 19 (PDF p. 23), quoting
# Unesco Report No. 38 p. 191. S = practical salinity (PSS-78), t = IPTS-68 temperature [°C], p = gauge pressure
# [dbar], rho [kg/m^3], V = specific volume [1e-3 m^3/kg]. The source warns that the last decimal may differ by round-off.
UNESCO_EOS80_CHECK = {
    "rows": [
        {"S": 0.0, "t68_C": 5.0, "p_dbar": 0.0, "rho_kg_m3": 999.96675, "V_1e-3_m3_kg": 1.000033251},
        {"S": 0.0, "t68_C": 5.0, "p_dbar": 10000.0, "rho_kg_m3": 1044.12802, "V_1e-3_m3_kg": 0.957736964},
        {"S": 0.0, "t68_C": 25.0, "p_dbar": 0.0, "rho_kg_m3": 997.04796, "V_1e-3_m3_kg": 1.00296078},
        {"S": 0.0, "t68_C": 25.0, "p_dbar": 10000.0, "rho_kg_m3": 1037.90204, "V_1e-3_m3_kg": 0.963482064},
        {"S": 35.0, "t68_C": 5.0, "p_dbar": 0.0, "rho_kg_m3": 1027.67547, "V_1e-3_m3_kg": 0.973069835},
        {"S": 35.0, "t68_C": 5.0, "p_dbar": 10000.0, "rho_kg_m3": 1069.48914, "V_1e-3_m3_kg": 0.935025857},
        {"S": 35.0, "t68_C": 25.0, "p_dbar": 0.0, "rho_kg_m3": 1023.34306, "V_1e-3_m3_kg": 0.977189409},
        {"S": 35.0, "t68_C": 25.0, "p_dbar": 10000.0, "rho_kg_m3": 1062.53817, "V_1e-3_m3_kg": 0.941142660},
    ],
    "t90_from_t68": "t90 = t68 / 1.00024",
    "_source": "N. P. Fofonoff and R. C. Millard Jr., 'Algorithms for computation of fundamental properties of "
               "seawater', Unesco Technical Papers in Marine Science 44 (1983), p. 19, check-value table for the "
               "EOS-80 (from Unesco Report 38, p. 191); scanned PDF via WHOI/MBL darchive, page image read 2026-09-13",
}


def main() -> int:
    HERE.mkdir(parents=True, exist_ok=True)
    (HERE / "constants.json").write_text(json.dumps(CODATA, indent=2), encoding="utf-8")
    (HERE / "ussa1976_constants.json").write_text(json.dumps(USSA_CONSTANTS, indent=2), encoding="utf-8")
    with open(HERE / "ussa1976_table1.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["Z_km", "H_km", "T_K", "p_Pa", "rho_kg_m3", "c_m_s", "g_m_s2"])
        w.writerows(USSA_TABLE1)
    with open(HERE / "ussa1976_table2.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["Z_km", "H_km", "mu_Pa_s", "nu_m2_s", "Hp_m", "n_per_m3", "V_m_s", "mfp_m", "M_kg_kmol"])
        w.writerows(USSA_TABLE2)
    with open(HERE / "iapws_sigma.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["t_C", "sigma_exp_mN_m", "uncertainty_mN_m", "sigma_calc_mN_m"])
        w.writerows(IAPWS_SIGMA)
    other = {"iapws_sigma_equation": IAPWS_SIGMA_EQ, "iapws_viscosity_check": IAPWS_VISCOSITY_CHECK,
             "jennings": JENNINGS, "taylor_blast": TAYLOR, "ams_lapse_rate": AMS,
             "unesco_eos80_check": UNESCO_EOS80_CHECK}
    (HERE / "benchmarks.json").write_text(json.dumps(other, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"wrote reference data to {HERE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
