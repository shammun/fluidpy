"""Write the cited public benchmark values and forms used by ``tests/test_ch08.py`` into ``reference/ch08/``.

Chapter 8 (laminar flow) needs no datasets: almost every check is Tier 1 (exact viscous solutions, sympy re-derivations,
convergence orders of our own discretisations). The public sources below pin the published *numbers* that the chapter's
code reproduces from its own implementation — the optimum taper and maximum load of the inclined slider bearing and the
taper of maximum peak pressure (San Andrés), the front constant η_N of the planar viscous gravity current (Huppert 1982
as summarised by Ball & Huppert), the exact elementary charge (NIST CODATA 2022) — and the standard *forms* of the Oseen
and Proudman–Pearson drag laws, Stokes' law and its pressure field, Stokes' second problem, the Hagen–Poiseuille law and
circular Couette flow (V1 "form" cross-checks: a published closed form reproduced identically is not an independent
number). Every entry was read from the cited page on 2026-09-28 (see ``SOURCES.md``). Nothing here comes from the
textbook; book-quoted numbers and printed closed forms live in the git-ignored ``tests/book_values_ch08.json``.

Run: ``.venv/Scripts/python.exe reference/ch08/make_refs.py``  (idempotent; overwrites ``benchmarks.json``).
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

BENCHMARKS = {
    "_meta": {
        "chapter": "ch08",
        "written_by": "reference/ch08/make_refs.py",
        "verified": "2026-09-28",
        "note": "Public sources only (pages fetched on the date above, see SOURCES.md). No book values.",
    },
    # L. San Andrés, Modern Lubrication, Notes 2 Appendix "One-dimensional fluid film bearings" (Texas A&M, © 2009, rev.
    # Aug 2012), pp. 5-6: load w = (6 μ U L² B / h2²) W(α), Eq. (8b) W(α) = [ln α + 2(1 − α)/(1 + α)]/(1 − α)², α = h1/h2
    # (inlet/exit film ratio); "αopt := 2.1889", "W(αopt) = 0.0267"; peak pressure P_max(α) = (α − 1)/(4α(1 + α))
    # (normalised by 6μUL/h2²) with "αPmax = 2.414", "Pmax(αPmax) = 0.043". PDF text extracted with PyMuPDF.
    "san_andres_slider": {
        "W_form": "W(K) = (log(K) + 2*(1 - K)/(1 + K))/(1 - K)**2, K = h_inlet/h_exit, w = 6*mu*U*L**2*B/h_exit**2 * W",
        "K_opt": 2.1889,
        "W_opt": 0.0267,
        "Pmax_form": "Pmax(K) = (K - 1)/(4*K*(1 + K)), normalised by 6*mu*U*L/h_exit**2",
        "K_Pmax": 2.414,
        "Pmax_max": 0.043,
        "mapping": "K = 1 + alpha (book gap h = h0(1 + alpha x/L), pad moving +x: inlet h0(1+alpha) at x = L, exit h0)",
        "source": "https://rotorlab.tamu.edu/me626/Notes_pdf/Notes02_App_1D_bearings.pdf (PDF fetched 2026-09-28, "
                  "text extracted with PyMuPDF, pp. 5-6 read)",
        "digits": "as printed (K_opt 5 significant figures; W, Pmax 3 figures; K_Pmax 4 figures)",
        "label": "V5",
    },
    # T. V. Ball & H. E. Huppert, "Similarity solutions and viscous gravity current adjustment times" (preprint, JFM
    # style), Appendix (a) "Two-dimensional viscous gravity currents" (Huppert 1982): h_t − β(h³h_x)_x = 0 (A1),
    # ∫_0^{x_N} h dx = A (A2), h = (3/10)^{1/3} η_N^{2/3} (A²/β)^{1/5} t^{−1/5}(1 − y²)^{1/3} (A4),
    # x_N = η_N(βA³)^{1/5}t^{1/5} (A5), η_N = [(1/5)(3/10)^{1/3} π^{1/2} Γ(1/3)/Γ(5/6)]^{−3/5} = 1.411… (A6a,b).
    "huppert_planar_current": {
        "pde": "h_t - beta*(h**3*h_x)_x = 0, integral_0^{x_N} h dx = A (half-area)",
        "h_form": "(3/10)**(1/3) * eta_N**(2/3) * (A**2/beta)**(1/5) * t**(-1/5) * (1 - (x/x_N)**2)**(1/3)",
        "xN_form": "eta_N * (beta*A**3)**(1/5) * t**(1/5)",
        "eta_N_form": "((1/5)*(3/10)**(1/3)*sqrt(pi)*gamma(1/3)/gamma(5/6))**(-3/5)",
        "eta_N": 1.411,
        "digits": "printed as 1.411... (4 significant figures)",
        "source": "https://warwick.ac.uk/fac/sci/maths/people/staff/tball/publications/similarity_preprint.pdf "
                  "(PDF fetched 2026-09-28, text extracted with PyMuPDF, Appendix a, Eqs. A1-A7); citation of record "
                  "H. E. Huppert, J. Fluid Mech. 121, 43-58 (1982)",
        "label": "V5",
    },
    # Wikipedia, "Oseen equations": F = 6πμau(1 + (3/8)Re), Re = ρua/μ (radius); C_d = (12/Re)(1 + (3/8)Re);
    # Proudman & Pearson (1957): F = 6πμaU(1 + (3/8)Re + (9/40)Re² ln Re + O(Re²)).
    "oseen_drag": {
        "force_form": "F = 6*pi*mu*a*U*(1 + 3/8*Re_a), Re_a = rho*U*a/mu (radius)",
        "cd_form": "C_d = 12/Re_a*(1 + 3/8*Re_a)",
        "proudman_pearson_form": "F = 6*pi*mu*a*U*(1 + 3/8*Re_a + 9/40*Re_a**2*log(Re_a) + O(Re_a**2))",
        "oseen_coeff_radius": 0.375,
        "pp_coeff_radius": 0.225,
        "source": "https://en.wikipedia.org/wiki/Oseen_equations (fetched 2026-09-28)",
        "label": "V1 form (coefficients 3/8, 9/40 are exact rationals of the theory)",
    },
    # Wikipedia, "Stokes' law": F_d = −6πμRv; settling v = (2/9)(ρ_p − ρ_f)gR²/μ; p(r, θ) = −(3μRu/2) cos θ/r².
    "stokes_law": {
        "drag_form": "F = 6*pi*mu*R*v",
        "settling_form": "v = 2/9*(rho_p - rho_f)/mu*g*R**2",
        "pressure_form": "p - p_inf = -3*mu*R*u/2*cos(theta)/r**2 (theta from the direction of the stream)",
        "source": "https://en.wikipedia.org/wiki/Stokes%27_law (fetched 2026-09-28)",
        "label": "V1 form",
    },
    # Wikipedia, "Stokes problem": u = U e^{−√(ω/2ν) y} cos(ωt − √(ω/2ν) y); penetration depth δ = √(2ν/ω).
    "stokes_second_problem": {
        "u_form": "U*exp(-sqrt(omega/(2*nu))*y)*cos(omega*t - sqrt(omega/(2*nu))*y)",
        "penetration_depth": "sqrt(2*nu/omega)",
        "source": "https://en.wikipedia.org/wiki/Stokes_problem (fetched 2026-09-28)",
        "label": "V1 form",
    },
    # Wikipedia, "Hagen–Poiseuille equation": Δp = 8μLQ/(πR⁴); u = (G/4μ)(R² − r²); Darcy f = 64/Re, Re = ρvd/μ.
    "hagen_poiseuille": {
        "dp_form": "Delta_p = 8*mu*L*Q/(pi*R**4)",
        "u_form": "G/(4*mu)*(R**2 - r**2), G = -dp/dz",
        "darcy_f": "64/Re, Re = rho*V*d/mu (mean velocity, diameter)",
        "source": "https://en.wikipedia.org/wiki/Hagen%E2%80%93Poiseuille_equation (fetched 2026-09-28)",
        "label": "V1 form",
    },
    # Wikipedia, "Taylor–Couette flow": v_θ = Ar + B/r, A = Ω1(μ − η²)/(1 − η²), B = Ω1R1²(1 − μ)/(1 − η²),
    # μ = Ω2/Ω1, η = R1/R2; Rayleigh: stable iff (r v_θ)² increases outward, i.e. μ > η² for co-rotation.
    "taylor_couette": {
        "A_form": "Omega1*(mu_r - eta**2)/(1 - eta**2)",
        "B_form": "Omega1*R1**2*(1 - mu_r)/(1 - eta**2)",
        "rayleigh": "stable iff d(r*v)**2/dr >= 0 everywhere; for the Couette profile: mu_r > eta**2",
        "source": "https://en.wikipedia.org/wiki/Taylor%E2%80%93Couette_flow (fetched 2026-09-28)",
        "label": "V1 form",
    },
    # NIST CODATA 2022: e = 1.602 176 634 × 10⁻¹⁹ C (exact).
    "elementary_charge": {
        "e_C": 1.602176634e-19,
        "status": "exact (SI 2019), CODATA 2022",
        "source": "https://physics.nist.gov/cgi-bin/cuu/Value?e (fetched 2026-09-28)",
        "label": "V5",
    },
    # Reused (verified in reference/ch04/ on 2026-09-23): Morrison (2016) sphere-drag correlation, 24/Re for Re < 2.
    "morrison_sphere_drag": {
        "reuse": "reference/ch04/benchmarks.json → morrison_sphere_drag (core.similarity.sphere_drag_coefficient)",
        "claim_tested": "the correlation lies between Stokes 24/Re and Oseen (24/Re)(1 + 3Re/16) for 0.1 <= Re <= 5",
        "label": "V5 (correlation, 3 % class)",
    },
    # Reused (verified in reference/ch01/ on 2026-09-12): USSA-1976 Sutherland constants β = 1.458e-6, S = 110.4 K.
    "ussa1976_sutherland": {
        "reuse": "reference/ch01/ussa1976_constants.json (sutherland_beta_kg_per_m_s_sqrtK, sutherland_S_K)",
        "label": "V5 (input property)",
    },
}


def main() -> None:
    (HERE / "benchmarks.json").write_text(json.dumps(BENCHMARKS, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("wrote", HERE / "benchmarks.json")


if __name__ == "__main__":
    main()
