# ch08 — Laminar Flow: derivation review

Independent, fresh-context, read-only review of `fluidpy/core/{laminar,lubrication,creeping,diffusion}.py`,
`fluidpy/ch08_laminar_flow.py`, `scripts/ch08_*.py` and `tests/test_ch08.py` against the rendered page images
(p338–p374). At review time: `pytest tests/test_ch08.py -m "not slow"` → 96 passed (57.5 s); all 12 scripts exit 0
with `--no-show`.

## Must fix

| # | Where | Finding | Fix | Status |
|---|---|---|---|---|
| M1 | `fluidpy/core/laminar.py` `couette_poiseuille_state` | Velocity-extremum location has the wrong sign. From (8.5) $u = \frac{U}{h}y - \frac{1}{2\mu}\frac{dp}{dx}\,y(h-y)$, $\frac{du}{dy} = \frac{U}{h} - \frac{1}{2\mu}\frac{dp}{dx}(h-2y) = 0$ gives $y^* = \frac{h}{2} - \frac{\mu U}{h\,dp/dx}$; the code used `+`. Wrong `u_max`, `y_umax`, `u_min`, `y_umin` whenever $U \ne 0$ (e.g. $h$ = 1 cm, $U$ = 1 cm/s, $dp/dx$ = −1 Pa/m: returned 0.016 m/s at 4 mm, true 0.018 m/s at 6 mm). The test only covered $U = 0$, where $y^* = h/2$ either way. | `ys = h/2 - mu*U/(h*g)`; add a dense-grid test for $U \ne 0$, both signs of $dp/dx$ | **fixed** (`laminar.py:199`; agrees with a 2 000 001-point grid to 9–10 digits); test with the verifier |

Should-fix 2, 3 (83 docstrings), 4, 5, 6, 7a fixed in `fluidpy/` (ch08 tests 96/96, full suite 977 passed with
`-m "not slow"`); 1 and 7b (tests) with the verifier.

## Should fix

1. Oseen wake direction is not pinned by a public test (a θ ↔ π − θ swap in both functions would pass). The coded
   ψ (8.53) is right: it satisfies $E^4\psi = \frac{U}{\nu}\,\partial_x(E^2\psi)$ to ~1e-17 with the + sign. Add that
   V2 test plus "fluid-frame |ψ| larger at θ < π/2".
2. Validation labels that overstate: `thin_film_spread` "converged" (no order asserted) · `sphere_drag_quadrature`
   "converged" (Gauss–Legendre exact from n = 1 → analytic) · `similarity_ode_solve` "converged" (no order study) ·
   Oseen / Proudman–Pearson "benchmark" (a published formula, not data → analytic form cross-check).
3. Verification O1: 83 docstrings still say "Validation (planned)" — replace with the test names.
4. `couette_poiseuille_state` with $U = 0$, $dp/dx > 0$ returns `backflow=True` for a pure Poiseuille flow running
   in −x — a flow direction, not backflow (teaching risk in E1).
5. `inertia_viscous_ratio` docstring: say which Re (book: radius; module: diameter) and add the verification O3
   finding. **Resolved differently:** the θ-part of O3 is right (no $O(U^2a/r^2)$ term in the θ-component of
   $\mathbf{u}\cdot\nabla\mathbf{u}$; radial term $\frac{3U^2a}{4r^2}(2 - 3\sin^2\theta)$), but its "prefactor ≈ 1/8,
   crossover 10–20 $a/Re_a$" is wrong. With the viscous term $|\nu\nabla^2\mathbf{u}| = |\nabla p|/\rho =
   \frac{3\nu U a}{2r^3}(4\cos^2\theta + \sin^2\theta)^{1/2}$ the axis ratio is
   $\frac{1.5U^2a/r^2}{3\nu U a/r^3} = \tfrac12 Re_a \frac{r}{a}$, so the crossover is $r/a \approx 2/Re_a = 4/Re$ (checked by
   hand by the orchestrator; the function gives 0.498 at θ = 0, 0.499 at π/2, r/a = 400). The book's "$r/a \sim 1/Re$"
   is right as an order of magnitude.
6. `settling_state["valid"]` uses Re < 0.1 — our choice; the book says Stokes/Oseen are "fairly accurate for Re < 5".
7. Small: `stokes_first_vorticity(y, 0)` returns NaN (no t ≤ 0 guard) · test comment at `tests/test_ch08.py` ~1383
   says "stokes_first(2U)" while the (correct) assertion uses U.

## Printed-slip corrections (analysis §9) — all judged right

- R6 (8.13b) prints $\partial p/\partial x$ in the y-equation; (8.16b) has $\partial p^*/\partial y^*$.
- R7 (8.17a) and Example 8.2 lack ν; units do not balance.
- R8 Example 8.1 re-derived by hand: $p - p_e = \frac{6\mu L U}{h_0^2}\,\frac{\alpha\,s(1-s)}{(2+\alpha)(1+\alpha s)^2}$,
  $s = x/L$ — the square is correct; $C_1$, $C_2$, the O(α) pressure and W are right; the page's intermediate
  integrals use $(1 - \alpha x/L)$.
- R8b printed (8.19) gives $u(h) = U_h + U_0$. · R9 $V = \int u\,dy$ lacks $1/h$. · R10 $-\int \partial u/\partial y\,dy = +U$.
- R11 $2\,\mathrm{erf}^{-1}(0.95) = 2.772$. · R12 rear minimum $-3\mu U/2a$. · R13 Oseen needs $-\partial p/\partial x_i$.
- R14 "(9.63)", "(9.68)" are slips for (8.43), (8.48).
- R15 power into the fluid $= -2\pi R_1 \sigma_{R\varphi} u_\varphi$ (the fluid's face at $R_1$ has normal $-\mathbf{e}_R$).
- Exact slider-load series re-derived: $\alpha/2 - 3\alpha^2/4 + 33\alpha^3/40 - 13\alpha^4/16 + 171\alpha^5/224$ — matches.

## Verified (checked hardest)

- §8.2: (8.5) with the book's A, B; Q; τ; (8.6)–(8.8); pipe Q and V; f = 64/Re from (8.8). Circular Couette (8.9),
  (8.10), the $R_2 = \infty$ and $R_1 = 0$ branches, the pressure integral, $\sigma_{R\varphi} = -2\mu B/R^2$, torques
  $\pm 4\pi\mu B$, dissipation $= 4\pi\mu B(\Omega_1 - \Omega_2)$ = power in; Rayleigh flag switches at $\Omega_2/\Omega_1 = \eta^2$.
- §8.3: scaling coefficients vs (8.16a,b) incl. $\Lambda = \mu U L/(P_a h^2)$ (book's oil: $\varepsilon^2 Re_L = 10^{-3}$);
  Hele-Shaw; thin-film flux $q = -\frac{\rho g}{3\mu} h^3 h_x$; Newton Jacobian of `thin_film_spread` entry by entry;
  Huppert $x_N$, $h_c$; pad-frame Reynolds solver; inlet backflow switches at α = 1.
- §8.4: (8.30); $\omega = -\partial u/\partial y$; $\int\omega\,dy = +U$; δ = 3.643$\sqrt{\nu t}$; Example 8.5 (Galilean map
  uses U, not 2U); Example 8.6 as Lamb–Oseen ($\sigma^2 = 4\nu t$) and spin-up; Example 8.7 $n = m = 1/5$.
- §8.5: (8.38) with the bounded root; $e^{-4/\sqrt2} = 0.0591$ at $4\sqrt{\nu/\omega}$; CN θ-scheme with BCs at
  $t_{n+1}$ and two backward-Euler start-up steps.
- §8.6: $E^2$; (8.48)–(8.50) with θ from the downstream axis; surface tractions ($3\mu U/2a$ uniform x-traction);
  running drag $\pi\mu a U(1-\cos^3\theta)$ and $\pi\mu a U(2 - 3\cos\theta + \cos^3\theta)$; (8.52) with $Re = 2aU/\nu$;
  Oseen ψ and velocities; $C_D = \frac{24}{Re}\left(1 + \frac{3}{16}Re\right)$; Proudman–Pearson (radius form); settling and
  Millikan force balance.
- Reuse: `ch04.stokes_first_problem` re-exports the core function (same signature and values; returns 0 instead of
  NaN for t ≤ 0 — ch04 calls it with t > 0 only).

## Not caught by the tests as they stood

The M1 extremum for $U \ne 0$; an Oseen wake on the wrong side; explainer parity rows built on
`couette_poiseuille_state` keys never checked against a grid.

## Found later (explainer build)

- `slider_bearing_state["inlet_backflow"]` tested only x = L; recirculation next to the pad occurs where
  $h > 1.5\,h_m$ ($h_m$ = gap at the pressure peak), i.e. wide/narrow gap ratio > 2: α > 1 (wide end x = L) **or**
  α < −½ (wide end x = 0). **Fixed** in `fluidpy/core/lubrication.py` (flag now checks the wide end; new keys
  `backflow_x`, `backflow_any`); brute-force 201 × 399 pad-frame scan agrees in every case; ch08 tests 99 passed.
  Follow-up: `viz/ch08/slider_bearing.html` Explain §7 still says the flag tests x = L only (for the viz review).
