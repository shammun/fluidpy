"""Published reference values for the V5 evidence level, each with its source.

Every entry carries a citation. **Rule: nothing goes in here from memory.** Before a value is used in a test, the
verifier confirms it against the cited source (or the source's own table in ``reference/chNN/``) and records the check
in ``reference/chNN/SOURCES.md`` with the date. Values marked ``VERIFY`` below are the usual textbook figures supplied
as a starting point — confirm and then delete the marker.

Import in tests::

    from tools.benchmarks import BLASIUS, WALL_LAW, STABILITY, cite

    assert abs(fpp0 - BLASIUS["fpp0"]) < 1e-5, cite("BLASIUS")
"""
from __future__ import annotations

SOURCES: dict[str, str] = {
    "BLASIUS": "Howarth, Proc. R. Soc. Lond. A 164, 547 (1938); tabulated in Schlichting & Gersten, "
               "Boundary-Layer Theory, 9th ed., Ch. 6.",
    "FALKNER_SKAN": "Falkner & Skan (1931); separation value tabulated in Schlichting & Gersten, Ch. 7.",
    "WALL_LAW": "Kim, Moin & Moser, J. Fluid Mech. 177, 133 (1987); Moser, Kim & Mansour, Phys. Fluids 11, 943 (1999).",
    "STABILITY": "Orszag, J. Fluid Mech. 50, 689 (1971) (plane Poiseuille); Chandrasekhar, Hydrodynamic and "
                 "Hydromagnetic Stability (1961) (Rayleigh-Benard).",
    "GHIA": "Ghia, Ghia & Shin, J. Comput. Phys. 48, 387 (1982), doi:10.1016/0021-9991(82)90058-4 — tables I-II.",
    "SOD": "Sod, J. Comput. Phys. 27, 1 (1978) — or an exact Riemann solver written for this project.",
    "PIPE": "Moody, Trans. ASME 66, 671 (1944); Colebrook, J. Inst. Civ. Eng. 11, 133 (1939).",
    "DRAG": "Schlichting & Gersten, Boundary-Layer Theory; Roshko, NACA TN 3169 (1954) (Strouhal number).",
}


def cite(key: str) -> str:
    """The citation string for an assertion message, so a failing test says where the number came from."""
    return f"{key}: {SOURCES[key]}"


# --- flat-plate laminar boundary layer -------------------------------------------------------------------------------
BLASIUS = {
    "fpp0": 0.332057,          # f''(0)                                        VERIFY against Howarth's table
    "delta99_coeff": 5.0,      # delta_99 = 5.0 * x / sqrt(Re_x)               VERIFY (4.91 with a different definition)
    "delta_star_coeff": 1.7208,  # displacement thickness coefficient          VERIFY
    "theta_coeff": 0.664,      # momentum thickness coefficient                VERIFY
    "cf_coeff": 0.664,         # c_f = 0.664 / sqrt(Re_x)                      VERIFY
    "cd_coeff": 1.328,         # C_D = 1.328 / sqrt(Re_L) (one side)           VERIFY
}

FALKNER_SKAN = {
    "beta_separation": -0.198838,   # wedge parameter at which f''(0) = 0      VERIFY
}

# --- wall-bounded turbulence -----------------------------------------------------------------------------------------
WALL_LAW = {
    "kappa": 0.41,             # von Karman constant (0.38-0.41 in the literature - state which you use)
    "B": 5.0,                  # log-law intercept (5.0-5.2)
    "y_plus_viscous": 5.0,     # edge of the viscous sublayer
    "y_plus_log_start": 30.0,  # start of the log layer
}

# --- linear stability ------------------------------------------------------------------------------------------------
STABILITY = {
    "poiseuille_Re_c": 5772.22,     # plane Poiseuille, based on centreline velocity and half-gap   VERIFY (Orszag)
    "poiseuille_alpha_c": 1.02056,  # critical streamwise wavenumber                                VERIFY
    "rayleigh_benard_rigid": 1707.762,   # critical Rayleigh number, rigid-rigid                    VERIFY
    "rayleigh_benard_free": 657.511,     # critical Rayleigh number, free-free                      VERIFY
    "pipe_transition_Re": 2300.0,        # conventional transition value (not a stability limit)
}

# --- drag ------------------------------------------------------------------------------------------------------------
DRAG = {
    "stokes_sphere_coeff": 24.0,   # C_D = 24/Re for Re << 1
    "cylinder_strouhal": 0.20,     # St over roughly 300 < Re < 2e5                                 VERIFY (Roshko)
}

# --- pipe flow -------------------------------------------------------------------------------------------------------
PIPE = {
    "laminar_friction_coeff": 64.0,   # f = 64/Re (Darcy friction factor)
}

TABLES = {
    # name -> (chapter, file under reference/<chNN>/, what it contains).  make_refs.py creates the file; SOURCES.md
    # records where it came from. Keep the raw numbers out of this module.
    "ghia_re100": ("ch??", "ghia_re100.csv", "u on the vertical centreline, v on the horizontal centreline, Re=100"),
    "sod_exact": ("ch??", "sod_exact.npz", "exact Riemann solution of Sod's problem at t=0.2"),
}
