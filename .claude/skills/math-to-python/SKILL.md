---
name: math-to-python
description: How to turn textbook mathematics (definitions, theorems, derivations, numbered equations) into faithful, testable Python — transcription discipline, symbolic-first workflow, choosing the numerical method, grid and boundary conventions, non-dimensionalisation, units, and the numerical pitfalls that bite in fluid dynamics. Load before writing or reviewing any implemented code.
---

# math-to-python — coding a book that has no code

## 0. The order of work (do not skip a step)
1. **Transcribe** the equation from the page image into LaTeX in `analysis/chNN.md`. Never from the extracted text
   alone: in this book's extraction `¼` is `=`, `ð…Þ` are parentheses, minus signs vanish, ω prints as `u`, fractions
   are split across lines. Render the defining page: `.venv/Scripts/python.exe tools/render_pages.py chNN --eq N.M`.
2. **Symbolically** state it in sympy and check it against the book's stated special cases (the book almost always
   gives one: "for a flat plate this reduces to …"). If the special case does not come out, your transcription is wrong.
3. **Decide the numerical method** (table below) and what the book leaves unspecified.
4. **Implement** with the book's own symbols and a `# Eq. (n.m)` comment on the line that is the equation.
5. **Validate** (see `verify-implementation`) before you move on. An unvalidated function is a liability for every later
   chapter that reuses it.

## 1. Transcription discipline
- One function per equation, named after what it computes, not after the equation number.
- Keep the book's symbol names in the signature (`nu`, `rho`, `dpdx`, `U_inf`, `Re_x`); add units in the docstring.
- If the book writes a result in a compact form and its expanded form elsewhere, implement the form it *derives* and
  assert the compact form in a test — this catches transcription errors in both places.
- If you suspect a typo in the book (dimensions do not balance, the limit is wrong), do **not** silently fix it: code
  what you believe is correct, and record in the analysis and the report "book prints X; dimensional analysis requires Y;
  implemented Y" with the reasoning.

## 2. Choosing the numerical method
| The book gives you | Implement with | Watch for |
|---|---|---|
| a closed-form expression | direct numpy evaluation, vectorised | cancellation; evaluate the numerically stable form |
| an implicit algebraic equation (Colebrook, normal-shock M₂, Kutta condition) | `scipy.optimize.brentq` with a bracket derived from physics; `newton` only with an analytic derivative | brackets that fail at extreme parameters; always assert the residual |
| an integral (streamfunction, drag, energy) | `scipy.integrate.quad` (adaptive, returns its error estimate — assert it), `simpson` on given samples | improper/singular integrands: substitute or split at the singularity |
| an ODE initial-value problem | `scipy.integrate.solve_ivp(..., rtol=1e-8, atol=1e-10, dense_output=True)`; default `RK45`, `Radau`/`BDF` when stiff | stiffness (chemistry, thin layers); events for shock/separation detection |
| an ODE boundary-value problem (Blasius, Falkner–Skan, Orr–Sommerfeld) | shooting (`brentq` on the missing initial condition) **and** `scipy.integrate.solve_bvp` — implement one, cross-check with the other | infinite domains: truncate at η_max and show the answer is insensitive to it |
| an eigenvalue problem (stability, normal modes) | `scipy.linalg.eig`/`eigh` on the discretised operator; Chebyshev collocation for Orr–Sommerfeld | spurious modes: refine and keep only the modes that converge |
| a PDE with a scheme the book specifies | write that scheme explicitly; do not substitute a library solver | boundary treatment; stability limit |
| a PDE with no scheme specified | pick the simplest that resolves it (2nd-order central in space, RK4 or Crank–Nicolson in time), document the choice, and show grid-independence | say in the docstring that this is *your* choice, not the book's |
| a linear system | `scipy.sparse.linalg.spsolve`/`cg` with a sparse assembly; dense only for n ≲ 2000 | assemble with `lil`/`coo`, solve in `csr`; check the condition number when a result looks odd |
| a series solution | sum with an explicit truncation test (add terms until the increment is < tol), never a fixed term count | slow convergence near the domain edges |
| a statistical/turbulence quantity | numpy/pandas with explicit averaging windows | ensemble vs time vs spatial averaging — say which |

## 3. Grids, fields and boundaries
- Standard layout: `x` varies along the **last** axis, `y` along the second-to-last (`u[j, i]` = `u(y_j, x_i)`), and
  `np.meshgrid(x, y, indexing="xy")`. State it once in `fluidpy/core/grids.py` and never deviate.
- Say in every docstring whether a value lives at a **node** or a **cell face/centre**. Half-cell offsets are the single
  most common discretisation bug; a staggered (MAC) arrangement is worth it for incompressible solvers.
- Boundary conditions are part of the operator, not an afterthought: implement Dirichlet/Neumann/periodic explicitly in
  `core/operators.py` with a `bc=` argument, and unit-test each with a function whose derivative is known.
- `np.gradient` is second-order in the interior and **first-order at the edges** — do not use it inside a scheme whose
  order you are about to measure. Write the difference stencils explicitly.
- Time stepping: expose `dt` and assert the stability criterion (`dt <= CFL*dx/|u|max`, `dt <= dx^2/(2*nu)`), raising a
  clear error rather than producing NaNs.

## 4. Non-dimensionalisation and units
- Implement the non-dimensional equations (that is what the book analyses) and provide a thin dimensional wrapper
  built from the reference scales; never mix.
- Put the reference scales in one place per chapter (a small dataclass `Scales(L, U, rho, mu)` with a
  `nondimensionalise()`/`redimensionalise()` pair) so a change of reference length cannot silently leak.
- `fluidpy/core/units.py` wraps a `pint` registry. Use it in tests (dimensional homogeneity) and at script/notebook
  boundaries; keep the hot numerical paths as plain floats.
- Temperatures in kelvin inside functions. Pressures absolute unless the name says `_gauge`.

## 5. Numerical pitfalls that actually occur in fluid dynamics
- **Cancellation**: `sqrt(1+x) - 1` for small x, `(p - p_inf)/p_inf` for weak flows, `1 - M^2` near M = 1 → use the
  algebraically rearranged form (`x/(sqrt(1+x)+1)`), `np.expm1`, `np.log1p`.
- **Singularities**: potential-flow sources at their own location, the leading edge of a boundary layer, the pipe
  centreline in cylindrical coordinates (`r = 0`), the stagnation point in a log-law. Handle them analytically and test
  the limit explicitly.
- **Stiffness**: near-wall grids, chemical source terms, very low Mach numbers → implicit solver, and say so.
- **Transcendental branches**: `arctan2`, not `arctan`, for angles; `np.unwrap` for phase; watch the branch cut in
  conformal maps (Joukowski).
- **Interpolation**: books plot; you compute. Use `scipy.interpolate` with the same order as the scheme, and never
  interpolate across a shock or a discontinuity.
- **Random anything**: seed it (`np.random.default_rng(0)`) and put the seed in the docstring, or the notebook will not
  reproduce.
- **Floating-point comparison**: no `==` on floats; compare with a tolerance tied to the quantity's scale.
- **Performance**: vectorise first; a 200×200 grid stepped 5000 times in pure Python loops takes minutes and will blow
  the notebook budget. If a loop is unavoidable, isolate it behind a function so it can later be `numba.njit`-ed, and
  cache expensive results to `outputs/`/`reference/` with a documented invalidation rule.

## 6. Figures that match the book
- Reproduce the book's axes, ranges, non-dimensional groups and curve labels — that is what makes the comparison
  meaningful — but draw them with our code and our style: `fluidpy.core.style.use_style()` (palette shared with the
  explainers, explicit axis labels **with units**, legend entries naming the parameter value). Never copy a book figure.
- Never set colours per point for a continuous field without a colourbar.
- Save with `fluidpy.core.style.savefig(fig, "chNN", "fig_<n>_<slug>")` (→ `outputs/chNN/`, git-ignored), and return the
  figure so the notebook can show it. Animations and interactive figures: skill `python-viz`.
- Functions that an explainer mirrors must accept plain floats (scalar-callable) so `selftest()` parity rows can call
  them as `chNN.fn(1.0, 10.0)`.

## 7. Lessons learned (append one line per new pitfall; the knowledge-keeper promotes chapter lessons here)
- (none yet)
