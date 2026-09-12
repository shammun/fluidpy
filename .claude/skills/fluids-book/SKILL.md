---
name: fluids-book
description: Conventions and context for the fluidpy project (Kundu, Cohen & Dowling, Fluid Mechanics 5e, turned into tested Python, teaching notebooks and interactive explainers) - how the book's chapters build on each other, file naming, docstring contract, units and notation, validation vocabulary, curation tiers and what "done" means for a chapter. Load whenever working on any chapter, module, test, notebook, explainer or knowledge file in this repo.
---

# fluids-book — project conventions

## The book in one paragraph (our summary)
Kundu, Cohen & Dowling (5th ed., Academic Press 2012; 16 chapters, printed pages 1–852 = PDF pages 28–879). Ch. 1 sets
the vocabulary (continuum, viscosity, statics, thermodynamics, stratification, dimensional analysis); Ch. 2 the maths
(Cartesian tensors, Gauss and Stokes); Ch. 3 kinematics (Lagrangian/Eulerian views, strain and rotation, Reynolds
transport); Ch. 4 the conservation laws and Navier–Stokes (with rotating frames, Bernoulli, Boussinesq, similarity) — the
trunk everything else grows from. Ch. 5 vorticity dynamics; Ch. 6 ideal (potential) flow; Ch. 7 gravity waves; Ch. 8
exact laminar solutions; Ch. 9 boundary layers; Ch. 10 CFD; Ch. 11 instability; Ch. 12 turbulence; Ch. 13 geophysical
fluid dynamics (rotation + stratification: geostrophy, Ekman layers, shallow water, Kelvin and Rossby waves, baroclinic
instability — the chapter Shammunul's climate-dynamics PhD leans on most, drawing on 4, 5, 7, 8, 11, 12); Ch. 14
aerodynamics; Ch. 15 compressible flow; Ch. 16 biofluid mechanics. Appendices A–D (constants and properties, maths
tools, founders, visual resources) are reference only.

## What is different about this project
- The book has **no code**: every function is new; the chapter text and its equation *images* are the only specification.
- Correctness comes from mathematics and physics (`verify-implementation`), not from matching an author's output.
- We teach selectively: the `concept-curator` tiers every item CORE / SUPPORT / NOTE / SKIP; CORE ideas get the full
  treatment (and possibly an explainer), SKIP items are listed with the reason it is safe to skip them.

## File naming
- Chapter ids: `ch01` … `ch16`; slug from `book.yaml` (`ch07` → `gravity_waves`).
- Chapter module: `fluidpy/ch07_gravity_waves.py`; tests `tests/test_ch07.py`; notebook builder
  `notebooks/build_ch07.py` → `notebooks/ch07_gravity_waves.ipynb` (+ `_colab.ipynb`, `.html`); explainers
  `viz/ch07/<slug>.html`; reports `reports/ch07_{verification,review,viz}.md`; knowledge `knowledge/ch07.md`.
- Reusable primitives: `fluidpy/core/<topic>.py` (grids, operators, ode, potential, waves, similarity, thermo, …).
  Machinery (not physics): `project`, `embed`, `anim`, `interact`, `style`, `units`, `refdata`.
- Scripts mirror examples and figures: `scripts/ch07_particle_orbits.py`. Figures → `outputs/ch07/` (git-ignored).

## Notation and units (the most common source of bugs)
- SI everywhere in dimensional code; every docstring lists each symbol with its unit.
- Non-dimensional groups keep the book's names (`Re`, `Fr`, `Ro`, `Ri`, `Ra`, `Pr`, `Ma`, `We`, `St`) and the docstring
  says which reference scales the book uses *in that chapter*.
- Angles in radians inside functions; degrees only at the interface with `_deg` in the name.
- Coordinates: the book usually takes z (or y) positive upward with the free surface/wall at a stated level — record it
  per chapter in `knowledge/notation.md` and never assume it carries over.
- Pressure absolute unless the name says `_gauge`; kinematic viscosity `nu` vs dynamic `mu` never interchanged.

## Docstring contract (every implemented function)
```python
def phase_speed(k: float | np.ndarray, H: float, g: float = 9.81) -> np.ndarray:
    """Phase speed of linear surface gravity waves on water of depth H.

    Book: §7.2, Eq. (7.36)–(7.37)  (derived from the linearised free-surface conditions (7.12), (7.18), (7.21)).
    Parameters
    ----------
    k : wavenumber [1/m], > 0.   H : undisturbed depth [m], > 0.   g : gravitational acceleration [m/s^2].
    Returns
    -------
    c : phase speed [m/s].
    Assumptions: inviscid, irrotational, constant density, small amplitude (ka << 1), no surface tension.
    Validation: V2 sympy (satisfies the dispersion relation); V7 deep (kH→∞: c→sqrt(g/k)) and shallow (c→sqrt(gH)) limits.
    Label: symbolic.
    """
```
(The equation numbers above are illustrative — always take them from the rendered page.)

## Validation vocabulary (exact words in reports and `knowledge/concept_map.md`)
`analytic` (V1) · `symbolic` (V2) · `converged` (V3) · `conserved` (V4) · `benchmark` (V5, cited) · `book-value` (V6,
private) · `qualitative` (must appear in Open items) · `unverified` (blocks PASS unless justified).

## What "done" means for a chapter
1. `analysis/chNN.md` (inventory + dependency graph), `analysis/chNN_curation.md` (spine, tiers, ≤5 explainers),
   `analysis/chNN_design.md` (storyboards).
2. `fluidpy/chNN_<slug>.py` implements every CORE/SUPPORT item and every function in design Part C; scripts run.
3. `tests/test_chNN.py` passes; CORE ≥ 2 evidence levels, SUPPORT ≥ 1; `reports/chNN_verification.md` PASS;
   `reports/chNN_review.md` has no open Must-fix.
4. `viz/chNN/*.html` (≤ 5) pass lint + shot at every size with parity rows; `reports/chNN_viz.md` PASS.
5. The notebook executes headlessly with 0 errors; it follows the teaching style and embeds every explainer once.
6. `knowledge/` updated; `progress.json` all `pass`/`done`; published: page, Colab twin, gallery and index live,
   `tools/check_public.py` clean.

## Public-repo rule
The repo and site are public. Never commit or publish the book's text, page images, figure crops or scans, transcribed
tables, exercise text, or book-quoted numbers (those live in git-ignored `tests/book_values_chNN.json`). Equations may
be shown with their numbers; explanations are ours; figures are generated by our code.
