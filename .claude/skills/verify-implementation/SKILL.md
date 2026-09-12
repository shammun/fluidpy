---
name: verify-implementation
description: How to prove that Python written from a textbook's equations is actually correct when there is no reference implementation — the V1-V7 evidence ladder (analytic, symbolic, convergence, conservation, benchmark, book value, limits), tolerances, test skeletons, the benchmark catalogue and the verification report format. Load when testing, verifying or reviewing any chapter.
---

# verify-implementation — evidence that the physics is right

The book ships no code, so there is nothing to diff against. Correctness is established by making the mathematics test
itself. Evidence is proportional to the curation tier (`analysis/chNN_curation.md` §2): every computable **CORE** item
(every new idea of the chapter) needs **at least two independent evidence levels**, at least one of them V1, V2, V3 or
V5; a conceptual CORE item with no computable output gets a demonstration test (the notebook's code for it runs and
asserts its point); coded **NOTE** items at least one level; every function the notebook or an explainer calls (design
Part C) at least a smoke test with a physical sanity assertion. RECAP items were tested in their own chapter; SKIP
items are not coded.
"Independent" means a bug that fools one level must not fool the other: an analytic check and a symbolic check of the
*same* wrong equation are not independent; an analytic check plus a conservation check are.

---

## V1 — Analytic truth (the strongest, use it wherever a closed form exists)
Compare the code to an exact solution of the *same* problem on a field of points, not at one point.

Classic exact solutions to test against (state which one, and its assumptions, in the test name):
Couette and plane/pipe Poiseuille flow · Stokes' first (Rayleigh) and second problems (erf / exponentially damped
oscillation) · starting/decaying flows by separation of variables · Lamb–Oseen and Rankine vortices · potential-flow
superpositions (source+sink+uniform = Rankine oval, doublet = cylinder, with and without circulation; the Kutta–
Joukowski lift) · hydrostatics and manometry · 1-D isentropic relations, normal-shock jump conditions, Prandtl–Meyer ·
the Sod shock-tube exact solution · Stokes flow past a sphere (drag 6πμaU) · Womersley pulsatile flow · linear
gravity/acoustic wave dispersion relations · Ekman and Stokes layers.

```python
def test_stokes_first_problem_matches_erf():
    y = np.linspace(0, 0.05, 200); t, nu, U = 2.0, 1e-6, 1.0
    u = ch05.stokes_first_problem(y, t, nu, U)
    exact = U * erfc(y / (2 * np.sqrt(nu * t)))          # Book Eq. (5.31)
    assert np.max(np.abs(u - exact)) < 1e-10
```

## V2 — Symbolic re-derivation (sympy) and dimensional homogeneity
Two distinct checks:
1. **Residual check** — substitute the coded solution into the governing equation and simplify to zero:
```python
import sympy as sp
y, h, mu, dpdx = sp.symbols("y h mu dpdx", positive=True)
u = dpdx / (2*mu) * (y**2 - h**2)                        # Book Eq. (5.14)
assert sp.simplify(mu * sp.diff(u, y, 2) - dpdx) == 0    # governing equation
assert sp.simplify(u.subs(y,  h)) == 0                   # no-slip, Eq. (5.15)
```
2. **Dimensional homogeneity** — every term of every implemented equation must have the same dimensions. Use
   `fluidpy/core/units.py` (a `pint` registry) and assert the function's output unit:
```python
from fluidpy.core.units import Q_, dimensional_check
def test_dynamic_pressure_units():
    dimensional_check(lambda rho, U: 0.5*rho*U**2, "pressure", rho=Q_(998, "kg/m**3"), U=Q_(2, "m/s"))
```
A dimensional-homogeneity test catches ν↔μ, missing ρ and missing ½ faster than any other test.

## V3 — Convergence / order of accuracy (for every discretisation)
Refine and measure; do not eyeball. `tools/convergence.py` provides `observed_order(h, err)` (least-squares slope of
log err vs log h) and `richardson(f_coarse, f_fine, p)`.
* With an **exact solution**: errors at h, h/2, h/4 → observed order within ±0.15 of the design order.
* Without one: **method of manufactured solutions** — pick a smooth analytic field, substitute it into the operator to
  get a source term, solve with that source, and measure the error against the manufactured field.
* Also assert the *stability* limit the book states (CFL, diffusion number): the scheme must blow up just above it and
  be stable just below.
```python
def test_laplacian_second_order():
    hs, errs = [], []
    for n in (32, 64, 128, 256):
        x = np.linspace(0, 1, n+1); u = np.sin(2*np.pi*x)
        num = core.operators.laplacian_1d(u, x[1]-x[0])
        errs.append(np.max(np.abs(num[1:-1] + (2*np.pi)**2*u[1:-1]))); hs.append(x[1]-x[0])
    assert abs(observed_order(hs, errs) - 2.0) < 0.15
```

## V4 — Conservation laws and invariants
Assert the budget, not the vibe: total mass, momentum and energy; `div(u)` to machine precision for an incompressible
solver; circulation (Kelvin's theorem) for an inviscid barotropic flow; vorticity in 2-D inviscid flow; entropy
non-decreasing across a shock; Bernoulli constant along a streamline in steady inviscid flow; the Reynolds transport
identity between the integral and differential forms the book derives.
Report the residual as a number in the verification report; a solver whose mass drifts by 1e-3 per step is broken even
if the picture looks right.

## V5 — Published benchmarks (with a citation, never from memory)
Store the numbers in `reference/chNN/` via `make_refs.py`, and record each in `reference/chNN/SOURCES.md`
(value, source, DOI/URL, date verified). `tools/benchmarks.py` holds the ones already digitised.
Candidates — **verify each against the cited source before using it**:

| Benchmark | Value to reproduce | Usual source |
|---|---|---|
| Blasius flat plate | f''(0) = 0.332057; δ99 = 5.0 x/√Re_x; c_f = 0.664/√Re_x | Howarth (1938); Schlichting |
| Falkner–Skan | separation at β = −0.198838 | Schlichting, *Boundary-Layer Theory* |
| Lid-driven cavity | centreline u and v profiles at Re = 100/400/1000 | Ghia, Ghia & Shin, JCP 48 (1982) |
| Sod shock tube | densities/pressures at t = 0.2 | Sod, JCP 27 (1978) — or your own exact Riemann solver |
| Pipe friction | f = 64/Re (laminar); Colebrook–White implicit form; Moody chart | Moody (1944); Colebrook (1939) |
| Rayleigh–Bénard | Ra_c = 1707.762 (rigid–rigid), 657.51 (free–free) | Chandrasekhar (1961) |
| Plane Poiseuille stability | Re_c = 5772.22 at α = 1.02056 | Orszag, JFM 50 (1971) |
| Taylor–Couette | critical Taylor number, narrow gap | Taylor (1923); Chandrasekhar |
| Sphere/cylinder drag | C_D = 24/Re (Stokes); standard drag curve; Strouhal ≈ 0.2 in the vortex-shedding range | Schlichting; Roshko (1954) |
| Turbulent channel | law of the wall u+ = y+ and u+ = (1/κ)ln y+ + B, κ ≈ 0.41, B ≈ 5.0–5.2; DNS profiles at Re_τ = 180 | Kim, Moin & Moser, JFM 177 (1987) |
| Fluid properties | ρ, μ, c_p of water/air vs T | NIST / CoolProp / the ICAO standard atmosphere |

## V6 — Numbers and figures the book prints
Grep `chapters/chNN.txt` for `= `, `%`, `Re =`, "approximately", worked-example answers and table captions. Reproduce
each and tabulate `book value | our value | relative error | comment`. Because these are book-derived, keep the raw
table **private**: put the values in `tests/book_values_chNN.json` (git-ignored) and guard the tests with

```python
BOOK = Path("tests/book_values_ch05.json")
book_only = pytest.mark.skipif(not BOOK.exists(), reason="book values are private; see CLAUDE.md rule 9")
```
The published report may say "reproduces the book's Example 5.3 to 0.2 %" — not the book's table itself.

## V7 — Limits, symmetry and invariance (cheap, catches sign errors)
Asymptotic limits (Re → 0 and → ∞, μ → ∞, Ma → 0, γ → 1, gap → 0) · symmetry (mirror the geometry, the solution must
mirror) · Galilean invariance (add a uniform velocity; the pressure field must not change) · rotational invariance of
tensor quantities · monotonicity (drag must increase with speed) · positivity (density, absolute temperature,
turbulent kinetic energy) · reciprocity in Stokes flow.

---

## Default tolerances
| Quantity | Pass |
|---|---|
| closed-form comparison, analytic function | ≤ 1e-10 relative |
| ODE/BVP solution vs exact (rtol=1e-8 solver) | ≤ 1e-6 relative |
| observed order of accuracy | within ±0.15 of design order (±0.25 if the error range spans < 2 decades) |
| conservation residual, incompressible divergence | ≤ 1e-12 absolute (machine) or the stated truncation error |
| published benchmark table | ≤ 1 % of the reported value unless the source states its own uncertainty |
| digitised-from-a-plot benchmark | ≤ 3 %, and say "digitised" in the report |
| book-quoted value | ≤ 0.5 %, or explain the difference (rounding in the book, different property table) |
Tighten these where the algorithm is deterministic. **Never loosen one to make a test pass.**

## Test file skeleton (`tests/test_chNN.py`)
```python
import numpy as np, pytest, sympy as sp
from pathlib import Path
from tools.convergence import observed_order
from tools.benchmarks import BLASIUS
from fluidpy import ch09_boundary_layers as ch7   # (illustrative)

REF = Path("reference/ch07")
needs_ref = pytest.mark.skipif(not REF.exists(), reason="run reference/ch07/make_refs.py")

def test_blasius_wall_shear_benchmark():            # V5
    eta, f = ch7.blasius_similarity()
    assert abs(f[2, 0] - BLASIUS["fpp0"]) < 1e-5    # Howarth (1938), see reference/ch07/SOURCES.md

def test_blasius_satisfies_ode():                   # V2 (independent of the benchmark)
    eta, f = ch7.blasius_similarity()
    res = np.gradient(f[2], eta) + 0.5 * f[0] * f[2]
    assert np.max(np.abs(res[10:-10])) < 1e-6

def test_blasius_far_field_limit():                 # V7
    eta, f = ch7.blasius_similarity(eta_max=12.0)
    assert abs(f[1, -1] - 1.0) < 1e-8
```
Run with `.venv/Scripts/python.exe -m pytest tests/test_chNN.py -q -p no:cacheprovider`.

## Report format (`reports/chNN_verification.md`)
```
# Chapter N verification — <title>                     date, commit hash
## Environment: python x.y, numpy, scipy, sympy versions
## Validation table
| Concept / Eq. | Tier | fluidpy target | Evidence (V-levels) | Numbers (error, order, residual) | Label | Notes |
## Functions used by the notebook and explainers (design Part C) — test name each
## Convergence studies (scheme | grids | observed order | design order)
## Conservation / invariant residuals
## Benchmarks used (value | our value | source + DOI | date verified)
## Numbers from the text (book vs ours)  [private values redacted to percentages]
## Figures reproduced with our code (figure number | outputs/chNN/verify/<png> (local) | one-sentence visual verdict)
## Deviations & justifications (every `# DEVIATION` in the code appears here)
## Open items (anything qualitative/unverified, missing benchmark, needs user)
## Verdict: PASS / FAIL
```

## Failure loop
If a test fails: (1) check the *test* first — wrong reference length, degrees vs radians, a benchmark quoted at a
different Reynolds number, comparing a non-dimensional result to a dimensional one; (2) then the code; (3) rerun.
Maximum 3 rounds, then write the discrepancy in Open items with your best hypothesis and stop.

Hypotheses that explain most failures in this domain, in the order they actually occur:
sign convention (pressure gradient, stress tensor, y-axis direction) · a factor of 2 or ½ (dynamic pressure, half-gap
vs full gap, half-angle) · ν vs μ · the wrong reference length or velocity in a Reynolds number · boundary condition
applied at the node instead of the face (half-cell offset) · degrees vs radians · gauge vs absolute pressure ·
per-unit-depth vs total force · a non-dimensionalisation that changed between book sections · an equation the book
states for one regime being used in another.

## Public-repo rule
Benchmark tables are committed only with a citable public source. Book tables, page crops and figure comparisons stay
local (`tests/book_values_*.json`, `reports/**/figures/` — both git-ignored). `tools/check_public.py` enforces it.
