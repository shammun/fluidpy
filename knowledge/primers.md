# Primer register — concepts and tools already explained in a notebook

Appended by the knowledge-keeper after every chapter from the notebook's `metadata.fluidpy.primers` and the design's
prerequisite ledger. Later chapters do not repeat a primer: they write a one-sentence reminder ("primed in Ch. 1,
P44") and point here. IDs are the notebook's `P` numbers (not in numeric order inside ch01).

| Term (maths / physics / Python) | Explained in (chapter · notebook section · CORE block) | One-line gist (our words) |
|---|---|---|
| **Python and tooling** | | |
| matplotlib figures (P01) | ch01 · §1.1 · front | `fig, ax = plt.subplots()`, `ax.plot`, labels with units; every figure gets What you see / How to read it / What would change if |
| pint quantities (P02) | ch01 · §1.2 · R01 | `Q_(value, "unit")`, `.to()`, `.dimensionality`; °C is an offset unit, formulas use kelvin |
| numpy arrays (P03) | ch01 · §1.2 · R03 | arithmetic acts on every element at once (vectorised), no loops |
| f-strings (P04) | ch01 · §1.2 | `f"{x:.3g} Pa"` inserts a value with 3 significant figures and a unit |
| np.linspace and np.logspace (P06) | ch01 · §1.3 | evenly spaced numbers on a linear axis, or with a constant factor between neighbours on a log axis |
| np.random.default_rng (P10) | ch01 · §1.4 · C06 | a seeded generator makes "random" results reproducible; `normal`, `poisson`, `random` |
| tuple unpacking (P14) | ch01 · §1.4 · C06 | `a, b = f()` takes a function's several return values apart in one line |
| assert np.allclose (P15) | ch01 · §1.4 · C06 | stops with an error if two results differ beyond `rtol`; our proof that a from-scratch version matches the library |
| animate and show_animation (P16) | ch01 · §1.4 · C06 | `update(i)` moves already-drawn artists per frame; `player="frames"` steps, `player="video"` plays smoothly |
| slider_figure (P17) | ch01 · §1.4 · C06 | computes every slider position in Python up front, so the plotly figure keeps working on the published page |
| show_viz (P18) | ch01 · §1.4 · C06 | embeds a chapter explainer (local file in Jupyter, Pages copy in Colab, full-window block on the page) |
| Python dictionaries (P23) | ch01 · §1.5 · C12 | a dict maps names to values; fluidpy returns several named results this way |
| np.gradient (P22) | ch01 · §1.5 · C12 | central differences inside, one-sided at the ends; use `edge_order=2` for second order at the edges |
| functions as arguments and lambda (P29) | ch01 · §1.7 · C20 | functions are values; `lambda z, p: 1000.0` is a one-line density function passed to `integrate_hydrostatic` |
| scipy.integrate.solve_ivp (P31) | ch01 · §1.7 · C20 | adaptive Runge–Kutta that picks its own steps to meet a tolerance |
| sympy (P40) | ch01 · §1.8 · C35 | algebra with symbols: `symbols`, `Function`, `diff`, `simplify`; an expression that simplifies to 0 is an identity |
| plotly 3-D surface (P41) | ch01 · §1.9 · C40 | `go.Surface` draws a rotatable surface, `go.Scatter3d` adds points; stays interactive on the page |
| np.where and np.select (P46) | ch01 · §1.10 · C51 | element-wise choice between values; several regimes at once (how "stable/neutral/unstable" labels are assigned) |
| live widgets (P47) | ch01 · §1.10 · C51 | ipywidgets sliders that re-run a function while a kernel runs; frozen on the page, so pair with a static figure |
| itertools.combinations (P55) | ch01 · §1.11 · C67 | every way to choose k of n items, unordered, no repeats (35 ways to pick 3 of 7 columns) |
| np.linalg.det and np.linalg.matrix_rank (P56) | ch01 · §1.11 · C67 | floating-point determinant (round before comparing to −1) and rank with a tiny tolerance |
| sympy Matrix and nullspace (P61) | ch01 · §1.11 · C69 | exact integer/rational matrices; `.rank()` and `.nullspace()` without round-off |
| fractions.Fraction (P60) | ch01 · §1.11 · C69 | exact rationals so an exponent ½ never becomes 0.49999 |
| np.linalg.solve (P57) | ch01 · §1.11 · C69 | solves B·x = b for a square nonsingular B |
| **Physics vocabulary** | | |
| stress (P05) | ch01 · §1.3 | force per area, split into normal (push/pull) and shear (slide) parts; shear strain γ is the lean angle, dγ/dt its rate |
| mole, kilomole, molecular weight and Avogadro's number (P07) | ch01 · §1.3 | a kmol is 6.022×10²⁶ molecules; M_w = mass of a kmol; one molecule weighs M_w/A_o; n = ρA_o/M_w |
| temperature as molecular kinetic energy (P08) | ch01 · §1.4 · C06 | ½m⟨u²⟩ = (3/2)k_BT: hotter gas means faster molecules |
| Newton's second law and momentum (P09) | ch01 · §1.4 · C06 | net force = rate of change of momentum; draw a free-body diagram, add forces with signs |
| Gaussian velocity components and mean square speed (P11) | ch01 · §1.4 · C06 | each component is Gaussian with variance k_BT/m, so ⟨u²⟩ = 3⟨u_x²⟩ |
| boundary conditions (P20) | ch01 · §1.5 · C12 | what the solution does at the edges; no-slip at a wall; steady vs transient |
| weight and gravitational acceleration (P24) | ch01 · §1.7 · C20 | W = mg downward, g = 9.80665 m/s²; fluid weight ρVg |
| net force from pressure (P28) | ch01 · §1.7 · C20 | pressure pushes along the inward normal; the net force is −∮p n̂ dA, for a box a sum over six faces |
| internal and kinetic energy per unit mass (P32) | ch01 · §1.8 · C25 | e is the molecules' random-motion energy [J/kg]; bulk motion adds u²/2; thermodynamics here uses e |
| perfect-gas law as a working model (P33) | ch01 · §1.8 · C25 | until §1.9 derives it: pv = RT with R = 287 and e = C_vT with C_v = 718 J/(kg K) |
| **Mathematics** | | |
| Poisson counting statistics (P12) | ch01 · §1.4 · C06 | rare independent hits: variance = mean, relative scatter 1/√N̄ |
| power laws and log–log plots (P13) | ch01 · §1.4 · C06 | y = ax^k is a straight line of slope k on log–log axes; `np.polyfit` on the logs measures k |
| exponent rules (P43) | ch01 · §1.4 · C06 | a^m a^n = a^(m+n), (a^m)^n = a^(mn), (ab)^m = a^m b^m; e.g. (L³)^(−1/2) = L^(−3/2) |
| ordinary derivative as a slope (P19) | ch01 · §1.5 · C12 | du/dy = change of u per metre of y, unit (m/s)/m = 1/s |
| partial derivative (P25) | ch01 · §1.5 · C12 | slope in one variable with the others held fixed; ∂²f/∂y² = curvature; ∇ points uphill, fluxes run "down the gradient" |
| finite differences (P21) | ch01 · §1.5 · C12 | central difference is second order; FTCS steps ∂f/∂t = D∂²f/∂y² and is stable only for r = DΔt/Δy² ≤ ½ |
| first-order Taylor expansion (P26) | ch01 · §1.7 · C20 | f(z + dz) ≈ f(z) + f′dz; the neglected terms shrink like dz² and vanish after dividing by dz |
| definite integral (P27) | ch01 · §1.7 · C20 | adds pieces f dz between limits; ∫c dz = c(b − a); `np.trapezoid` approximates from samples |
| explicit stepping (P30) | ch01 · §1.7 · C20 | Euler y_{k+1} = y_k + fΔx (error ∝ Δx); Euler–Cromer updates velocity first and keeps oscillator energy bounded |
| differentials (P34) | ch01 · §1.8 · C25 | d = change of a state quantity (route-independent); δ = small amount transferred along a route (route-dependent) |
| line integral along a path (P35) | ch01 · §1.8 · C25 | ∫p dv adds p × dv along a route; the area under the route on a p–v diagram; different routes, different areas |
| natural logarithm and exponential (P36) | ch01 · §1.8 · C25 | ln undoes exp; ∫dv/v = ln(v₂/v₁); shrinking by e⁻¹ ≈ 0.368 is one e-folding |
| trapezoid rule (P37) | ch01 · §1.8 · C25 | integral from samples by trapezoids, error ∝ Δx²; `np.trapezoid`, `cumulative_trapezoid` |
| partial derivative with a variable held fixed (P39) | ch01 · §1.8 · C35 | (∂h/∂T)_p: the subscript names what is held fixed, and changing it changes the answer |
| product rule for differentials (P38) | ch01 · §1.8 · C35 | d(pv) = p dv + v dp; dropping one term is the classic slip |
| separation of variables (P42) | ch01 · §1.9 · C45 | put each variable on its own side and integrate: dp/p = γ dρ/ρ → ln p = γ ln ρ + C |
| square root of a negative number (P45) | ch01 · §1.10 · C50 | i² = −1; Euler's formula e^(iθ) = cos θ + i sin θ; cos(iσt) = cosh σt, so an imaginary frequency means growth |
| linear second-order ODE (P44) | ch01 · §1.10 · C50 | try e^(λt): λ² + N² = 0; N² > 0 cos/sin, N² < 0 cosh/sinh growth, N² = 0 double root ζ = A + Bt |
| inequalities under a sign change (P48) | ch01 · §1.10 · C54 | multiplying by −1 reverses an inequality — exactly the move between the two lapse-rate conventions |
| chain rule (P49) | ch01 · §1.10 · C54 | df = (∂f/∂x)_y dx + (∂f/∂y)_x dy; divide by dz for the rate along a path |
| Gibbs free energy (P50) | ch01 · §1.10 · C54 | g ≡ h − Ts, a device: dg = v dp − s dT involves only dp and dT |
| exact differentials and Maxwell relations (P51) | ch01 · §1.10 · C54 | equal mixed second derivatives of g give (∂v/∂T)_p = −(∂s/∂p)_T: an entropy derivative from measurable v(T) |
| logarithmic differentiation (P52) | ch01 · §1.10 · C55 | d ln f/dz = (1/f) df/dz; logs turn products of powers into weighted sums of relative rates |
| matrices, determinants and minors (P53) | ch01 · §1.11 · C67 | determinant is zero exactly when columns are dependent; a minor is the determinant of a square piece; 3×3 by cofactors |
| linear independence and rank (P54) | ch01 · §1.11 · C67 | rank = number of independent columns = size of the largest nonzero minor |
| null space and rank–nullity theorem (P58) | ch01 · §1.11 · C69 | all k with A·k = 0; rank + null-space dimension = number of columns |
| scaling the base units (P59) | ch01 · §1.11 · C69 | changing m, kg, s by factors multiplies a quantity's number by λ_M^a λ_L^b λ_T^c; the quantity is unchanged |

Glosses (one sentence where used, no demo) in ch01: elastic shear modulus G (§1.3), vapour pressure and cavitation
(§1.3), κ = k/ρC_p (§1.5), dot product and unit normal (§1.3, §1.5), radius of curvature and principal radii (§1.6),
limits and orders of smallness (§1.7), similar triangles (Ex. 1.3), light wavelength and intensity (Ex. 1.5), similarity
solution (Ex. 1.4), boundary layer (§1.4, §1.5), Prandtl number (§1.5), thermocline and inversion (§1.10), potential
density forward gloss in C60 (§1.10), implicit-function rule in the D19 sympy cell (§1.10).
