"""Write the cited public benchmark values used by ``tests/test_ch09.py`` into ``reference/ch09/benchmarks.json``.

Chapter 9 (boundary layers) is almost entirely Tier 1 (exact similarity solutions, sympy re-derivations, convergence
orders of our own marching scheme).  The public numbers below pin the *values* that our code reproduces from its own
implementation: the Blasius wall-shear constant, the Falkner-Skan wall-shear table (including the zero-shear separation
member and the second, reversed-flow branch), the Hiemenz stagnation-point constant, the Bickley plane-jet constants,
Thwaites' linear fit and the flat-plate friction laws.  Every entry was read from the cited page or PDF on 2026-09-30
(see ``SOURCES.md``).  Nothing here comes from the textbook; book-quoted numbers live in the git-ignored
``tests/book_values_ch09.json``.

Normalisation warning (checked on the source page): Belden et al. write the Falkner-Skan equation as
f''' + f f'' + beta(1 - f'^2) = 0 with eta = y sqrt((m+1)/2 * U/(nu x)), m = beta/(2 - beta); the book (and fluidpy) write
f''' + ((n+1)/2) f f'' - n f'^2 + n = 0 with eta = y sqrt(U_e/(nu x)), n = m.  The two f''(0) differ by
f''_book(0) = sqrt((m+1)/2) * kappa (the tests convert with ``belden_to_book``).

Run: ``.venv/Scripts/python.exe reference/ch09/make_refs.py``  (idempotent; overwrites ``benchmarks.json``).
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

BENCHMARKS = {
    "_meta": {
        "chapter": "ch09",
        "written_by": "reference/ch09/make_refs.py",
        "verified": "2026-09-30",
        "note": "Public sources only (pages fetched on the date above, see SOURCES.md). No book values.",
    },
    # E. R. Belden, Z. A. Dickman, S. J. Weinstein, A. D. Archibee, E. Burroughs, N. S. Barlow, "Asymptotic Approximant for the
    # Falkner-Skan Boundary-Layer equation", arXiv:1907.09912 (17 Jul 2019), Fig. 2 caption: (beta, kappa = f''(0), B) triples from a
    # shooting method with RK4 step 1e-4 and domain-length convergence L = 10..40; equation (1) f''' + f f'' + beta(1 - f'^2) = 0.
    "belden_falkner_skan": {
        "equation": "f''' + f f'' + beta (1 - f'^2) = 0, f(0) = f'(0) = 0, f'(inf) = 1; eta = y sqrt((m + 1)/2 * U/(nu x)), m = beta/(2 - beta)",
        "points": [
            {"tag": "I", "beta": 0.5, "kappa": 0.927680039836653},
            {"tag": "II (Blasius)", "beta": 0.0, "kappa": 0.469599988361},
            {"tag": "III", "beta": -0.12, "kappa": 0.28176052424},
            {"tag": "IV (separation)", "beta": -0.198837735, "kappa": 0.0},
            {"tag": "V (reversed branch)", "beta": -0.12, "kappa": -0.1429351943576},
            {"tag": "VI (reversed branch)", "beta": -0.02, "kappa": -0.065168585542904},
        ],
        "digits": "as printed in the caption (kappa 12-15 significant figures; beta_sep 9 figures)",
        "source": "https://arxiv.org/abs/1907.09912 (PDF text extracted with PyMuPDF, Figure 2 caption and Eq. (1)-(2))",
        "label": "V5",
    },
    # Wikipedia "Blasius boundary layer": small-eta expansion f(eta) = alpha eta^2/2 + O(eta^5) with alpha = 0.332057336215196 (for
    # f''' + f f''/2 = 0).  The same page also prints 0.332043934904293 where it converts the boundary condition at infinity; the two
    # numbers are presented as distinct on the page and differ by 4e-5 relative; only the first is the converged constant.
    "blasius_wikipedia": {
        "fpp0": 0.332057336215196,
        "fpp0_other_number_on_the_same_page": 0.332043934904293,
        "delta_star_coeff_rounded": 1.72,
        "theta_coeff_rounded": 0.665,
        "note": "0.665 is a rounding of 0.6641 (2 f''(0)); the exact value is 0.664115",
        "source": "https://en.wikipedia.org/wiki/Blasius_boundary_layer (read 2026-09-30)",
        "label": "V5 (rounded thicknesses: V1-level cross-check only)",
    },
    # P. Weidman & M. R. Turner, "Stagnation-point flows with stretching surfaces: A unified formulation and new results" (preprint,
    # U. Colorado / U. Surrey), Section 5.1: at sigma = 0 "pure Hiemenz flow" F0''' + F0 F0'' - F0'^2 + 1 = 0 (this is the book's
    # Falkner-Skan equation at n = 1), leading-order wall shear F''(0) = 1.232588 + 0.18533 sigma.
    "hiemenz_weidman_turner": {
        "equation": "F''' + F F'' - F'^2 + 1 = 0, F(0) = F'(0) = 0, F'(inf) = 1  (= book (9.36) at n = 1)",
        "Fpp0": 1.232588,
        "digits": "7 significant figures",
        "source": "https://personalpages.surrey.ac.uk/m.turner/Crane_Revision_FINAL.pdf (PDF text extracted with PyMuPDF, Eq. (5.3) and (5.5))",
        "label": "V5",
    },
    # Wikipedia "Stagnation point flow": asymptotic behaviour section prints the Hiemenz displacement thickness delta* = 0.6479 delta.
    "hiemenz_delta_star_wikipedia": {
        "delta_star_over_delta": 0.6479,
        "source": "https://en.wikipedia.org/wiki/Stagnation_point_flow (read 2026-09-30; only this number, no f''(0), is on the page)",
        "label": "V5 (4 significant figures)",
    },
    # Wikipedia "Bickley jet" (Bickley, Phil. Mag. 23 (1937) 727-731): u = 0.4543 (M^2/(nu rho^2 x))^(1/3) sech^2(xi),
    # xi = 0.2752 (M/(nu^2 rho))^(1/3) y / x^(2/3), Q = 2 rho int_0^inf u dy = 3.3019 (M nu rho^2 x)^(1/3),
    # v = 0.5503 (M nu/(rho x^2))^(1/3) (2 xi sech^2 xi - tanh xi).   M = momentum flux per unit span (the book's J).
    "bickley_jet_wikipedia": {
        "u0_coeff": 0.4543,
        "xi_coeff": 0.2752,
        "Q_coeff": 3.3019,
        "v_coeff": 0.5503,
        "source": "https://en.wikipedia.org/wiki/Bickley_jet (read 2026-09-30)",
        "label": "V5 (4-5 significant figures)",
    },
    # Agrawal, Bose, Griffin, Moin, "An extension of Thwaites' method for turbulent boundary layers" (J. Fluid Mech., under
    # consideration; arXiv:2310.16337), Sec. 2: L(m) ~= 0.45 + 6 m (2.4), where m = (theta^2/U_e) d2U/dn2|wall = -(theta^2/nu) dU_e/ds
    # = -lambda (Thwaites' m, by (2.2) nu d2U/dn2|wall = dP_e/ds = -U_e dU_e/ds) and L = 2 Re_theta dtheta/ds = 2l - 2(2+H)lambda (2.3);
    # theta^2 = 0.45 nu/U^6 int U^5 dr (2.5); laminar separation "when m ~= 0.09".
    "thwaites_agrawal": {
        "L_fit": "L = 0.45 + 6 m with m = -lambda (paper: m = (theta^2/Ue) d2U/dn2 at the wall = -(theta^2/nu) dUe/ds)",
        "theta_sq_form": "theta^2 = 0.45 nu / Ue^6 * int_0^s Ue^5 dr",
        "separation_abs_m": 0.09,
        "sign_note": "m = -lambda, so 0.45 + 6 m = 0.45 - 6 lambda, the book's linear fit",
        "source": "https://arxiv.org/abs/2310.16337 (PDF text extracted with PyMuPDF, Eqs. 2.4-2.5)",
        "label": "V5 (form)",
    },
    # Wikipedia "Skin friction drag" (read 2026-09-30): laminar c_f = 0.664 Re_x^(-1/2) (Blasius); turbulent c_f = 0.0576 Re_x^(-1/5)
    # "Prandtl's one-seventh-power law" (local coefficient).  Averaging c_f over the plate multiplies by 5/4: C_f,L = 0.0720 Re_L^(-1/5).
    "skin_friction_wikipedia": {
        "laminar_cf_coeff": 0.664,
        "turbulent_local_cf_coeff": 0.0576,
        "turbulent_mean_coeff_derived": 0.0720,
        "note": "fluidpy.plate_drag_coefficient uses 0.074 (a commonly quoted Prandtl mean-value coefficient); 2.8 % above 5/4 x 0.0576",
        "source": "https://en.wikipedia.org/wiki/Skin_friction_drag (read 2026-09-30)",
        "label": "V5 secondary (form and coefficient to 3 %)",
    },
    # Karman (1911/12) stable staggered street: b/a = arccosh(sqrt 2)/pi ~ 0.281 (secondary sources: Horvath et al., J. Geophys. Res.
    # Atmos. 125 (2020), doi:10.1029/2019JD032121; "A Stroll down Karman Street", Resonance 10(8), ias.ac.in; arXiv:1807.00203).
    "karman_ratio": {
        "b_over_a": 0.281,
        "closed_form": "arccosh(sqrt 2)/pi",
        "source": "https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JD032121 (search snippet quoted: r_K = cosh^-1(sqrt2)/pi ~ 0.281)",
        "label": "V5 (3 significant figures; the exact value is V1)",
    },
}


def main() -> None:
    out = HERE / "benchmarks.json"
    out.write_text(json.dumps(BENCHMARKS, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("wrote", out)


if __name__ == "__main__":
    main()
