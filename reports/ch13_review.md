# Chapter 13 review — Geophysical Fluid Dynamics (2026-10-07)

Independent read-only review of `fluidpy/core/gfd.py`, `core/vertical_modes.py`, `core/shallow_water.py`,
`fluidpy/ch13_geophysical_fluid_dynamics.py`, the `core/stability.py` diff, `scripts/ch13_*.py` and `tests/test_ch13.py`
against the rendered book pages. Suite re-run by the reviewer: 89 passed. Reviewer's own mutants: 34 of 36 killed
(the two survivors are test gaps; the code is right).

## Must fix
1. **`barotropic_run` (`core/shallow_water.py`) — invariants drift when `zeta0` is not band-limited.** The 2/3 mask is
   applied to the advective product only, never to the state. Inviscid, 64², t = 4: energy +18 %, enstrophy +130 % for
   white-noise input as given; −4e-5 / −2e-4 once the input is truncated first. Shipped caches and tests start from a
   pre-truncated field, so they are unaffected, but the docstring label `conserved` is false in general.
   Fix: project the initial state onto the retained modes when `dealias=True`; document; add a test with a
   non-band-limited start.
2. **`wkb_vertical_structure` docstring** says the error falls as 1/(Hm)²; measured and tested order is 1. Doc only.

## Should fix (orchestrator rulings in brackets)
- Westward flow over a step, §13.13: the reviewer and the verifier independently find the printed description
  inconsistent with potential-vorticity conservation. [Promote to slip #14, keeping the "true for a finite ridge" caveat.]
- Slip #9: four of the five wrong cross-references are confirmed; one of them is in fact correct. [Remove that one.]
- Slip #13: the row misdescribes the attribution. [Reword; keep it unasserted.]
- Slips #10 and #11 are order-of-magnitude statements that are loose or inconsistent, not false. [Keep the ids so
  numbering stays stable; mark both `kind = "loose"` and teach them as inconsistencies, not errors.]
- New trap: the text names an acceleration term of the vertical momentum equation as the vertical Coriolis force; the
  sign of the force is the opposite. [Add as trap T18.]
- Test gaps: `poincare_orbit` and `inertial_oscillation` never check dy/dt = v; `kelvin_state(wall="north")` is never
  exercised.
- Docstrings: `shallow_water_omega` sign statement holds only for k > 0; `rayleigh_kuo_eigs` compares two box values
  (the radiating far field is the issue; shooting gives 0.0575); Chebyshev barotropic digits; `barotropic_run` drift
  order is 5; `absolute_vorticity_gradient` and `ekman_residual` claim a convergence order no test asserts; one script
  comment about the wind changing sign.
- Argument order: `eady_max_growth_rate` / `eady_time_scale` take (f, N) while siblings take (N, f);
  `poincare_omega(K, f, c)` vs `shallow_water_omega(k, l, c, f0, beta)`. [Signatures stay; notebook and explainers call
  them by keyword.]
- NaN band: keep planetary-scale external-mode examples at 35°.
- Public-repo hygiene (the automatic check passes): a storyboard example, a few test inputs and one notebook-builder
  line sit too close to the book's own "typical" values or wording. [Replace with our illustrative inputs; book-only
  tests read the private file.]

## Verified against the rendered pages
- Slips #1–#8 false as printed; slip #5 re-read; `ekman_surface_printed` is a trap (true for f > 0).
- §13.4–13.5 thin-shell equations, geostrophy, thermal wind; Ekman layers (surface, bottom, finite depth, pumping,
  Sverdrup) for both signs of f; vertical modes, including weight-1 orthogonality for both lids; the shallow-water
  set, the cubic and its limits, C-grid stencils, 1-D potential-vorticity conservation, adjustment energy; internal
  waves and lee waves; Rossby waves and the Rayleigh–Kuo operator (the `stability.py` change is correct); the Eady
  matrix entry by entry, tilt and heat-flux direction in both hemispheres; Fjørtoft's argument; step flow.
- Not re-read on the page: two late equations of §13.16 (standard form assumed); the online benchmark source.

## Status
Must-fix items sent to the implementer; test additions to the verifier. See the verification report for the re-run.
