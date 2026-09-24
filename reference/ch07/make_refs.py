"""Write the cited public benchmark values and forms used by ``tests/test_ch07.py`` into ``reference/ch07/``.

Chapter 7 (gravity waves) needs no datasets: almost every check is Tier 1 (exact linear-wave solutions, sympy
re-derivations, convergence orders of our own discretisations). The public sources below pin the published *numbers* the
chapter's code reproduces from its own implementation — the maximum error and optimal exponent of the Fenton & McKee
(1990) explicit wavelength formula, the maximum error of Guo's (2002) formula, the air–water capillary minimum, the
steepness of the limiting Stokes wave, the depth decay of deep-water Stokes drift (true V5 benchmarks) — and the
standard *forms* of the capillary–gravity dispersion relation, Stokes drift, the Bélanger jump relations, the KdV
equation, the cnoidal and solitary waves and the Stokes amplitude dispersion (V1 "form" cross-checks: a published closed
form reproduced identically is not an independent number). Every entry was read from the cited page on 2026-09-24 (see
``SOURCES.md``). Nothing here comes from the textbook; book-quoted numbers and printed closed forms live in the
git-ignored ``tests/book_values_ch07.json``.

Run: ``.venv/Scripts/python.exe reference/ch07/make_refs.py``  (idempotent; overwrites ``benchmarks.json``).
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

BENCHMARKS = {
    "_meta": {
        "chapter": "ch07",
        "written_by": "reference/ch07/make_refs.py",
        "verified": "2026-09-24",
        "note": "Public sources only (pages fetched on the date above, see SOURCES.md). No book values.",
    },
    # J. D. Fenton & W. D. McKee, "On calculating the lengths of water waves", Coastal Engineering 14, 499-513 (1990),
    # Eq. (21): L = (gT^2/2π){tanh[(2π√(d/g)/T)^{3/2}]}^{2/3}; abstract: "exact in the limits of short and long waves, and
    # in the intermediate range has an accuracy always better than 1.7%"; p. 507: "The minimum value of the maximum
    # error over all wave lengths in approximating Eq. 14 was for ν = 1.49. Rounding to ν = 3/2 ... this approximation
    # has a maximum error of only 1.7%". Family: kd = x^2 [coth(x^ν)]^{1/ν}, x = ω√(d/g).
    "fenton_mckee_1990": {
        "formula": "k*d = x**2 * (coth(x**nu))**(1/nu), x = omega*sqrt(d/g), nu = 3/2 (Eq. 21: L = L0*tanh((2*pi*sqrt(d/g)/T)**1.5)**(2/3))",
        "max_error_wavelength_percent": 1.7,
        "max_error_statement": "accuracy always better than 1.7% (abstract); maximum error of only 1.7% (p. 507)",
        "optimal_nu": 1.49,
        "limits": "exact for long (shallow) and short (deep) waves",
        "source": "https://johndfenton.com/Papers/Fenton90c+McKee-On-calculating-the-lengths-of-water-waves.pdf "
                  "(PDF fetched 2026-09-24, text extracted with pymupdf) — Coastal Engineering 14, 499-513 (1990), "
                  "doi:10.1016/0378-3839(90)90032-R",
        "label": "V5",
    },
    # J. D. Fenton, "A note on two approximations to the linear dispersion relation for surface gravity water waves"
    # (October 26, 2006): Eq. (2) "accurate to within 1.5% over all wavelengths"; Eq. (3) (Guo with exponent 5/2) "has a
    # maximum error of about 0.7%, about half the error of equation (2)"; both exact in both limits.
    "fenton_2006_note": {
        "fm_statement_percent": 1.5,
        "guo_5_2_statement_percent": 0.7,
        "guo_over_fm": "about half",
        "guo_5_2_formula": "k*d = x**2 * (1 - exp(-x**2.5))**(-2/5)",
        "source": "https://johndfenton.com/Papers/Dispersion-Relation.pdf (PDF fetched 2026-09-24, text extracted)",
        "label": "V5 (loose statements; see SOURCES.md for the measured values)",
    },
    # J. Guo, "Simple and explicit solution of wave dispersion equation", Coastal Engineering 45, 71-74 (2002),
    # doi:10.1016/S0378-3839(02)00039-X: shape parameter β = 2.4908 and "a maximum relative error of 0.75%" — read from
    # search-engine summaries of the abstract (ScienceDirect, ResearchGate, academia.edu returned 403 on fetch).
    "guo_2002": {
        "formula": "k*h = x**2 * (1 - exp(-x**beta))**(-1/beta), x = omega*sqrt(h/g)",
        "beta": 2.4908,
        "max_rel_error_percent": 0.75,
        "source": "doi:10.1016/S0378-3839(02)00039-X (abstract via search results 2026-09-24; publisher page 403)",
        "label": "V5 (search-confirmed)",
    },
    # Wikipedia, "Capillary wave": ω² = |k|(((ρ − ρ')/(ρ + ρ'))g + σk²/(ρ + ρ')) (two deep fluids); for air–water the
    # minimum phase speed c_m = 0.23 m/s at λ_m = 1.7 cm, λ_m = 2π√(σ/((ρ − ρ')g)).
    "capillary_wave": {
        "dispersion": "omega**2 = abs(k)*(((rho - rho_p)/(rho + rho_p))*g + sigma/(rho + rho_p)*k**2)",
        "lambda_m": "2*pi*sqrt(sigma/((rho - rho_p)*g))",
        "c_min_m_s": 0.23,
        "lambda_m_cm": 1.7,
        "digits": "as printed (2 significant figures)",
        "source": "https://en.wikipedia.org/wiki/Capillary_wave (fetched 2026-09-24)",
        "label": "V5 (numbers) + V1 form",
    },
    # IAPWS R1-76(2014): σ(20 °C) = 72.74 mN/m (already reproduced by ch01.surface_tension_water, V5 in tests/test_ch01.py).
    "iapws_surface_tension": {
        "sigma_20C_mN_m": 72.74,
        "source": "https://iapws.org/public/documents/CH-L9/Surf-H2O-2014.pdf (reused from ch01; read 2026-09-24 by the "
                  "concept analyst)",
        "label": "V5 (input property)",
    },
    # HandWiki, "Physics:Stokes wave": "The maximum wave steepness, for periodic and propagating deep-water waves, is
    # H/λ = 0.1410633 ± 4·10⁻⁷ ... a sharp wave crest – with an angle of 120°"; the page cites Schwartz & Fenton for
    # the number; third-order deep water c = ω/k = (1 + ½(ka)²)√(g/k) + O((ka)⁴).
    "stokes_limit": {
        "H_over_lambda": 0.1410633,
        "uncertainty": 4e-7,
        "crest_angle_deg": 120.0,
        "attribution_on_page": "Schwartz & Fenton (cited by HandWiki); crest angle 2π/3 also in Dyachenko, Lushnikov & "
                               "Korotkevich, arXiv:1507.02784 (Stud. Appl. Math. 137, 419-472, 2016)",
        "phase_speed_third_order": "c = (1 + (k*a)**2/2)*sqrt(g/k)",
        "source": "https://handwiki.org/wiki/Physics:Stokes_wave (fetched 2026-09-24); https://arxiv.org/abs/1507.02784",
        "label": "V5 (constant) + V1 form (speed)",
    },
    # Wikipedia, "Stokes drift": deep water ū_S ≈ ωka²e^{2kz} = (4π²a²/λT)e^{4πz/λ}; "at a depth of a quarter wavelength,
    # z = −λ/4, it is about 4% of its value at the mean free surface".
    "stokes_drift": {
        "deep_form": "omega*k*a**2*exp(2*k*z)",
        "quarter_wavelength_fraction_percent": 4,
        "source": "https://en.wikipedia.org/wiki/Stokes_drift (fetched 2026-09-24)",
        "label": "V5 (the 4 % number) + V1 form",
    },
    # Wikipedia, "Hydraulic jumps in rectangular channels": y₂/y₁ = ½(√(1 + 8Fr₁²) − 1); ΔE = (y₂ − y₁)³/(4y₁y₂).
    "belanger": {
        "depth_ratio": "(sqrt(1 + 8*Fr1**2) - 1)/2",
        "head_loss": "(y2 - y1)**3/(4*y1*y2)",
        "source": "https://en.wikipedia.org/wiki/Hydraulic_jumps_in_rectangular_channels (fetched 2026-09-24)",
        "label": "V1 form",
    },
    # Wikipedia, "Cnoidal wave": ∂tη + √(gh)∂xη + (3/2)√(g/h)η∂xη + (1/6)h²√(gh)∂x³η = 0; η = η₂ + H cn²((x − ct)/Δ | m),
    # η₂ = (H/m)(1 − m − E/K), Δ = h√(4mh/(3H)), c = √(gh)[1 + (H/(mh))(1 − m/2 − (3/2)E/K)], λ = 2ΔK(m);
    # solitary limit η = H sech²((x − ct)/Δ), c = √(gh)(1 + H/(2h)), Δ = h√(4h/(3H)); Ursell U = Hλ²/h³.
    "cnoidal": {
        "kdv": "eta_t + sqrt(g*h)*eta_x + 3/2*sqrt(g/h)*eta*eta_x + 1/6*h**2*sqrt(g*h)*eta_xxx = 0",
        "trough": "H/m*(1 - m - E/K)",
        "Delta": "h*sqrt(4*m*h/(3*H))",
        "c": "sqrt(g*h)*(1 + H/(m*h)*(1 - m/2 - 3/2*E/K))",
        "wavelength": "2*Delta*K",
        "solitary": "H*sech((x - c*t)/Delta)**2, c = sqrt(g*h)*(1 + H/(2*h)), Delta = h*sqrt(4*h/(3*H))",
        "ursell": "H*lambda**2/h**3",
        "note": "Wikipedia's H (height) is the book's a for the solitary wave; K, E complete elliptic integrals of "
                "parameter m",
        "source": "https://en.wikipedia.org/wiki/Cnoidal_wave (fetched 2026-09-24)",
        "label": "V1 form",
    },
    # Mowbray & Rarity, J. Fluid Mech. 28, 1-16 (1967): four beams (St Andrew's cross) from an oscillating body in a
    # linearly stratified fluid; angle to the vertical arccos(ω/N) — citation (pedigree), not a number.
    "st_andrews_cross": {
        "citation": "D. E. Mowbray & B. S. H. Rarity, J. Fluid Mech. 28, 1-16 (1967)",
        "url": "https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/theoretical-and-"
               "experimental-investigation-of-the-phase-configuration-of-internal-waves-of-small-amplitude-in-a-"
               "density-stratified-liquid/BD4159D2420DA60E9640520CD06B68AB",
        "label": "pedigree (qualitative)",
    },
}


def main() -> int:
    out = HERE / "benchmarks.json"
    out.write_text(json.dumps(BENCHMARKS, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
