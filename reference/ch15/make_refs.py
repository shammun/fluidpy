"""Write ``reference/ch15/benchmarks.json`` — the cited public values and forms the chapter-15 tests compare with.

Only statements read first-hand in a public source on the date given are stored (see SOURCES.md for the full citation
and what was read).  Nothing here comes from the textbook.  The primary source is NACA Report 1135 (Ames Research
Staff, 1953), a U.S. government report in the public domain: its table rows were read from the rendered page images of
the NASA copy and typed by hand; its equations are stored as *forms* (a formula as the report writes it, with the
report's equation number) and the tests evaluate our functions against the form at our own inputs.  Numbers are stored
exactly as the source prints them; ``decimals`` records how many decimals each column prints.

Run: ``.venv/Scripts/python.exe reference/ch15/make_refs.py``  (idempotent: rewrites the same bytes).
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
NACA = "Ames Research Staff, 'Equations, Tables, and Charts for Compressible Flow', NACA Report 1135 (1953)"
NACA_URL = "https://www.nasa.gov/wp-content/uploads/2023/03/equations-tables-charts-compressibleflow-report-1135.pdf"
TODAY = "2026-10-13"

BENCHMARKS = {
    "_note": "Cited public values and forms for compressible flow; our own numbers are computed, never typed. See SOURCES.md.",
    "naca1135_table1_subsonic": {
        "what": "Table I, subsonic flow, gamma = 7/5 (report p. 633): rows M = 0.50, 0.55, 0.60",
        "columns": ["M", "p/pt", "rho/rhot", "T/Tt", "beta", "q/pt", "A/A*", "V/a*"],
        "decimals": [2, 4, 4, 4, 4, 4, 4, 5],
        "rows": [
            [0.50, 0.8430, 0.8852, 0.9524, 0.8660, 0.1475, 1.3398, 0.53452],
            [0.55, 0.8142, 0.8634, 0.9430, 0.8352, 0.1724, 1.2550, 0.58506],
            [0.60, 0.7840, 0.8405, 0.9328, 0.8000, 0.1976, 1.1882, 0.63481],
        ],
        "source": NACA + ", Table I", "url": NACA_URL,
        "how_obtained": "page rendered at 4x and read as an image; typed by hand; every entry then reproduced by our functions",
        "verified": TODAY,
    },
    "naca1135_table2_supersonic": {
        "what": "Table II, supersonic flow, gamma = 7/5 (report p. 634): rows M = 1.60, 1.65, 1.70, 2.00, 2.05, 2.10; "
                "beta = sqrt(M^2 - 1), q = dynamic pressure, nu = Prandtl-Meyer angle, mu = Mach angle (degrees), "
                "subscript 2 = behind a normal shock, t = total (stagnation)",
        "columns": ["M", "p/pt", "rho/rhot", "T/Tt", "beta", "q/pt", "A/A*", "V/a*", "nu_deg", "mu_deg", "M2", "p2/p1",
                    "rho2/rho1", "T2/T1", "pt2/pt1", "p1/pt2"],
        "decimals": [2, 4, 4, 4, 3, 4, 3, 5, 3, 2, 4, 3, 3, 3, 4, 4],
        "rows": [
            [1.60, 0.2353, 0.3557, 0.6614, 1.249, 0.4216, 1.250, 1.42539, 14.861, 38.68, 0.6684, 2.820, 2.032, 1.388, 0.8952, 0.2628],
            [1.65, 0.2184, 0.3373, 0.6475, 1.312, 0.4162, 1.292, 1.45439, 16.338, 37.31, 0.6540, 3.010, 2.115, 1.423, 0.8760, 0.2493],
            [1.70, 0.2026, 0.3197, 0.6337, 1.375, 0.4098, 1.338, 1.48247, 17.810, 36.03, 0.6405, 3.205, 2.198, 1.458, 0.8557, 0.2368],
            [2.00, 0.1278, 0.2300, 0.5556, 1.732, 0.3579, 1.688, 1.63299, 26.380, 30.00, 0.5774, 4.500, 2.667, 1.688, 0.7209, 0.1773],
            [2.05, 0.1182, 0.2176, 0.5433, 1.790, 0.3478, 1.760, 1.65530, 27.748, 29.20, 0.5691, 4.736, 2.740, 1.729, 0.6975, 0.1695],
            [2.10, 0.1094, 0.2058, 0.5313, 1.847, 0.3376, 1.837, 1.67687, 29.097, 28.44, 0.5613, 4.978, 2.812, 1.770, 0.6742, 0.1622],
        ],
        "source": NACA + ", Table II", "url": NACA_URL,
        "how_obtained": "page rendered at 3.2x and read as an image; typed by hand; every entry then reproduced by our functions",
        "verified": TODAY,
    },
    "naca1135_forms": {
        "what": "equations of the report, by its own equation numbers; xi = p2/p1, theta = shock-wave angle, delta = "
                "deflection (our sigma, delta), tilde u and tilde v = velocity components normal and parallel to the shock",
        "eq_105": "M1^2 = [(gamma+1) xi + (gamma-1)] / (2 gamma)",
        "eq_106": "rho2/rho1 = u1/u2 = [(gamma+1) xi + (gamma-1)] / [(gamma-1) xi + (gamma+1)]",
        "eq_111": "pt2/pt1 = exp(-Delta s/R) = xi^(-1/(gamma-1)) {[(gamma+1) xi + (gamma-1)] / [(gamma-1) xi + (gamma+1)]}^(gamma/(gamma-1))",
        "weak_shock_in_xi": "pt2/pt1 = 1 - (gamma+1)/(12 gamma^2) (xi-1)^3 + (gamma+1)/(8 gamma^2) (xi-1)^4 + ...",
        "eq_124": "u1~ u2~ = a*^2 - (gamma-1)/(gamma+1) v~^2   (Prandtl's relation for an oblique shock)",
        "eq_127": "(p2 - p1)/(rho2 - rho1) = gamma (p2 + p1)/(rho2 + rho1)",
        "eq_132": "M2^2 = {(gamma+1)^2 M1^4 sin^2(theta) - 4 (M1^2 sin^2(theta) - 1)(gamma M1^2 sin^2(theta) + 1)} / "
                  "{[2 gamma M1^2 sin^2(theta) - (gamma-1)] [(gamma-1) M1^2 sin^2(theta) + 2]}",
        "eq_137": "V2^2/V1^2 = 1 - 4 (M1^2 sin^2(theta) - 1)(gamma M1^2 sin^2(theta) + 1) / [(gamma+1)^2 M1^4 sin^2(theta)]",
        "eq_138": "cot(delta) = tan(theta) [ (gamma+1) M1^2 / (2 (M1^2 sin^2(theta) - 1)) - 1 ]",
        "eq_139b": "tan(delta) = [M1^2 sin(2 theta) - 2 cot(theta)] / [2 + M1^2 (gamma + cos(2 theta))]",
        "eq_142": "pt2/pt1 = {(gamma+1) M1^2 sin^2(theta) / [(gamma-1) M1^2 sin^2(theta) + 2]}^(gamma/(gamma-1)) "
                  "{(gamma+1) / [2 gamma M1^2 sin^2(theta) - (gamma-1)]}^(1/(gamma-1))",
        "eq_145": "(p2 - p1)/q1 = 4 (M1^2 sin^2(theta) - 1) / [(gamma+1) M1^2]",
        "eq_146_gamma_7_5": {"cubic": [35, 216], "quartic": [245, 864],
                             "form": "pt2/pt1 = 1 - (35/216) m^3 + (245/864) m^4 + ..., m = M1^2 sin^2(theta) - 1"},
        "eq_147": "Delta s/R = [1/(gamma-1)] Delta s/c_v = 2 gamma/(3 (gamma+1)^2) m^3 - 2 gamma^2/(gamma+1)^3 m^4 + ..., m = M1^2 sin^2(theta) - 1",
        "eq_167": "sin^2(theta*) = {(gamma+1) M1^2 - (3-gamma) + sqrt((gamma+1)[(gamma+1) M1^4 - 2(3-gamma) M1^2 + (gamma+9)])} / (4 gamma M1^2)   (M2 = 1)",
        "eq_168": "sin^2(theta_dmax) = {(gamma+1) M1^2 - 4 + sqrt((gamma+1)[(gamma+1) M1^4 + 8(gamma-1) M1^2 + 16])} / (4 gamma M1^2)   (maximum deflection)",
        "eq_169": {"coefficient_gamma_7_5": 0.2946, "form": "delta* ~ (M1^2 - 1)^(3/2) / (sqrt(2) (gamma+1) M1^2)   (M1 close to 1, radians)"},
        "eq_170": {"coefficient_gamma_7_5": 0.3208, "form": "delta_max ~ 4 (M1^2 - 1)^(3/2) / (3 sqrt(3) (gamma+1) M1^2)   (M1 close to 1, radians)"},
        "wedge_limit_deg_gamma_7_5": 45.6,
        "wedge_limit_form": "the bow wave is attached for a wedge of semivertex angle less than asin(1/gamma) = 45.6 deg for gamma = 7/5",
        "eq_171c": "nu = sqrt((gamma+1)/(gamma-1)) atan(sqrt((gamma-1)/(gamma+1) (M^2-1))) - atan(sqrt(M^2-1))",
        "expansion_series": "p2/p1 = 1 - gamma M1^2/sqrt(M1^2-1) (Dnu) + gamma M1^2 [(gamma+1) M1^4 - 4 (M1^2-1)] / [4 (M1^2-1)^2] (Dnu)^2 - ...   (Dnu in radians)",
        "source": NACA + ", report pp. 620-626", "url": NACA_URL,
        "how_obtained": "pages rendered and read as images (equations transcribed by hand)",
        "verified": TODAY,
    },
    "normal_shock_tables": {
        "what": "normal-shock table for gamma = 1.4: M1, M2, p2/p1, rho2/rho1, T2/T1, p02/p01",
        "columns": ["M1", "M2", "p2/p1", "rho2/rho1", "T2/T1", "p02/p01"],
        "decimals": [4, 4, 4, 4, 4, 4],
        "rows": [
            [1.5, 0.7011, 2.4583, 1.8621, 1.3202, 0.9298],
            [2.0, 0.5774, 4.5000, 2.6667, 1.6875, 0.7209],
            [2.5, 0.5130, 7.1250, 3.3333, 2.1375, 0.4990],
            [3.0, 0.4752, 10.3333, 3.8571, 2.6790, 0.3283],
        ],
        "source": "Wikipedia, 'Normal shock tables' (secondary)", "url": "https://en.wikipedia.org/wiki/Normal_shock_tables",
        "how_obtained": "page fetched and the four rows extracted; the M1 = 2 row agrees with NACA Report 1135, Table II",
        "verified": TODAY,
    },
    "choked_flow": {
        "what": "critical pressure ratio p*/p0 = (2/(gamma+1))^(gamma/(gamma-1)) ('p* = 0.528 p0' for gamma = 1.4) and the "
                "page's table of minimum upstream/downstream pressure ratios for choked flow (the reciprocal); mass flow "
                "mdot = C_d A sqrt(gamma rho0 P0 (2/(gamma+1))^((gamma+1)/(gamma-1)))",
        "critical_ratio_air": 0.528,
        "min_upstream_over_downstream": [["dry air", 1.400, 1.893], ["helium", 1.660, 2.049], ["hydrogen", 1.410, 1.899],
                                         ["methane", 1.307, 1.837], ["propane", 1.131, 1.729], ["butane", 1.096, 1.708],
                                         ["chlorine", 1.355, 1.866], ["carbon monoxide", 1.404, 1.895]],
        "source": "Wikipedia, 'Choked flow' (secondary)", "url": "https://en.wikipedia.org/wiki/Choked_flow",
        "how_obtained": "page fetched and read",
        "verified": TODAY,
    },
    "fanno_flow": {
        "what": "adiabatic constant-area flow with friction, ratios to the sonic state of the same process",
        "differential": "dM^2/M^2 = gamma M^2/(1 - M^2) (1 + (gamma-1)/2 M^2) 4 f dx/D_h",
        "length": "4 f L*/D_h = (1 - M^2)/(gamma M^2) + (gamma+1)/(2 gamma) ln[ M^2 / ( (2/(gamma+1)) (1 + (gamma-1)/2 M^2) ) ]",
        "ratios": "with B = (2/(gamma+1))(1 + (gamma-1)/2 M^2): p/p* = 1/(M sqrt(B)), rho/rho* = sqrt(B)/M, T/T* = 1/B, V/V* = M/sqrt(B), "
                  "p0/p0* = B^((gamma+1)/(2(gamma-1)))/M",
        "entropy": "Delta s/c_p = ln[ M^((gamma-1)/gamma) B^(-(gamma+1)/(2 gamma)) ]",
        "source": "Wikipedia, 'Fanno flow' (secondary)", "url": "https://en.wikipedia.org/wiki/Fanno_flow",
        "how_obtained": "page fetched and the formulas read",
        "verified": TODAY,
    },
    "rayleigh_flow": {
        "what": "frictionless constant-area flow with heat exchange, ratios to the sonic state of the same process",
        "differential": "dM^2/M^2 = (1 + gamma M^2)/(1 - M^2) (1 + (gamma-1)/2 M^2) dT0/T0",
        "ratios": "p/p* = (gamma+1)/(1 + gamma M^2), rho/rho* = (1 + gamma M^2)/((gamma+1) M^2), T/T* = (gamma+1)^2 M^2/(1 + gamma M^2)^2, "
                  "v/v* = (gamma+1) M^2/(1 + gamma M^2), T0/T0* = 2 (gamma+1) M^2 (1 + (gamma-1)/2 M^2)/(1 + gamma M^2)^2, "
                  "p0/p0* = (gamma+1)/(1 + gamma M^2) [ (2/(gamma+1)) (1 + (gamma-1)/2 M^2) ]^(gamma/(gamma-1))",
        "mach_of_maximum_static_enthalpy": 0.845,
        "statement": "for subsonic flow the maximum of H occurs at M = 0.845; between there and M = 1 heat is added but T/T* decreases",
        "source": "Wikipedia, 'Rayleigh flow' (secondary)", "url": "https://en.wikipedia.org/wiki/Rayleigh_flow",
        "how_obtained": "page fetched and read (the page states 0.845 without the formula 1/sqrt(gamma); gamma = 1.4 is implied by its charts)",
        "verified": TODAY,
    },
}


def main() -> int:
    out = HERE / "benchmarks.json"
    out.write_text(json.dumps(BENCHMARKS, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {out} ({len([k for k in BENCHMARKS if not k.startswith('_')])} entries)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
